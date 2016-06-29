# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.

from __future__ import print_function, absolute_import

import os

import yaml
import numpy as np

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.image.utils as afwImageUtils
from lsst.afw.table import SourceCatalog, SchemaMapper, Field
from lsst.afw.table import MultiMatch, SimpleRecord, GroupView
from lsst.afw.fits.fitsLib import FitsError
import lsst.daf.persistence as dafPersist
import lsst.pipe.base as pipeBase
from lsst.utils import getPackageDir

# from .base import ValidateErrorNoStars
# from .calcSrd import calcAM1, calcAM2, calcAM3, calcPA1, calcPA2
# from . import calcSrd
# from .check import checkAstrometry, checkPhotometry, positionRms
# from .plot import plotAstrometry, plotPhotometry, plotPA1, plotAMx
# from .print import printPA1, printPA2, printAMx
from .srdSpec import srdSpec, loadSrdRequirements
from .util import getCcdKeyName, repoNameToPrefix, loadParameters
from .io import JobSerializer, persist_job
from .matchreduce import (MatchedMultiVisitDataset, AnalyticPhotometryModel,
                          AnalyticAstrometryModel, positionRms)
from .calcsrd import (PA1Measurement, PA2Measurement, PF1Measurement,
                      AMxMeasurement)
from .base import Metric


def loadAndMatchData(repo, dataIds,
                     matchRadius=afwGeom.Angle(1, afwGeom.arcseconds),
                     verbose=False):
    """Load data from specific visit.  Match with reference.

    Parameters
    ----------
    repo : string
        The repository.  This is generally the directory on disk
        that contains the repository and mapper.
    dataIds : list of dict
        List of `butler` data IDs of Image catalogs to compare to reference.
        The `calexp` cpixel image is needed for the photometric calibration.
    matchRadius :  afwGeom.Angle().
        Radius for matching.
    verbose : bool, optional
        Output additional information on the analysis steps.

    Returns
    -------
    afw.table.GroupView
        An object of matched catalog.
    """

    # Following
    # https://github.com/lsst/afw/blob/tickets/DM-3896/examples/repeatability.ipynb
    butler = dafPersist.Butler(repo)
    dataset = 'src'

    # 2016-02-08 MWV:
    # I feel like I could be doing something more efficient with
    # something along the lines of the following:
    #    dataRefs = [dafPersist.ButlerDataRef(butler, vId) for vId in dataIds]

    ccdKeyName = getCcdKeyName(dataIds[0])

    schema = butler.get(dataset + "_schema", immediate=True).schema
    mapper = SchemaMapper(schema)
    mapper.addMinimalSchema(schema)
    mapper.addOutputField(Field[float]('base_PsfFlux_snr', "PSF flux SNR"))
    mapper.addOutputField(Field[float]('base_PsfFlux_mag', "PSF magnitude"))
    mapper.addOutputField(Field[float]('base_PsfFlux_magerr', "PSF magnitude uncertainty"))
    newSchema = mapper.getOutputSchema()

    # Create an object that can match multiple catalogs with the same schema
    mmatch = MultiMatch(newSchema,
                        dataIdFormat={'visit': int, ccdKeyName: int},
                        radius=matchRadius,
                        RecordClass=SimpleRecord)

    # create the new extented source catalog
    srcVis = SourceCatalog(newSchema)

    for vId in dataIds:
        try:
            calexpMetadata = butler.get("calexp_md", vId, immediate=True)
        except FitsError as fe:
            print(fe)
            print("Could not open calibrated image file for ", vId)
            print("Skipping %s " % repr(vId))
            continue
        except TypeError as te:
            # DECam images that haven't been properly reformatted
            # can trigger a TypeError because of a residual FITS header
            # LTV2 which is a float instead of the expected integer.
            # This generates an error of the form:
            #
            # lsst::pex::exceptions::TypeError: 'LTV2 has mismatched type'
            #
            # See, e.g., DM-2957 for details.
            print(te)
            print("Calibration image header information malformed.")
            print("Skipping %s " % repr(vId))
            continue

        calib = afwImage.Calib(calexpMetadata)

        oldSrc = butler.get('src', vId, immediate=True)
        print(len(oldSrc), "sources in ccd %s  visit %s" % (vId[ccdKeyName], vId["visit"]))

        # create temporary catalog
        tmpCat = SourceCatalog(SourceCatalog(newSchema).table)
        tmpCat.extend(oldSrc, mapper=mapper)
        tmpCat['base_PsfFlux_snr'][:] = tmpCat['base_PsfFlux_flux'] / tmpCat['base_PsfFlux_fluxSigma']
        with afwImageUtils.CalibNoThrow():
            (tmpCat['base_PsfFlux_mag'][:], tmpCat['base_PsfFlux_magerr'][:]) = \
                calib.getMagnitude(tmpCat['base_PsfFlux_flux'],
                                   tmpCat['base_PsfFlux_fluxSigma'])

        srcVis.extend(tmpCat, False)
        mmatch.add(catalog=tmpCat, dataId=vId)

    # Complete the match, returning a catalog that includes
    # all matched sources with object IDs that can be used to group them.
    matchCat = mmatch.finish()

    # Create a mapping object that allows the matches to be manipulated
    # as a mapping of object ID to catalog of sources.
    allMatches = GroupView.build(matchCat)

    return allMatches


def analyzeData(allMatches, safeSnr=50.0, verbose=False):
    """Calculate summary statistics for each star.

    Parameters
    ----------
    allMatches : afw.table.GroupView
        GroupView object with matches.
    safeSnr : float, optional
        Minimum median SNR for a match to be considered "safe".
    verbose : bool, optional
        Output additional information on the analysis steps.

    Returns
    -------
    pipeBase.Struct containing:
    - mag: mean PSF magnitude for good matches
    - magerr: median of PSF magnitude for good matches
    - magrms: standard deviation of PSF magnitude for good matches
    - snr: median PSF flux SNR for good matches
    - dist: RMS RA/Dec separation, in milliarcsecond
    - goodMatches: all good matches, as an afw.table.GroupView;
        good matches contain only objects whose detections all have
          * a PSF Flux measurement with S/N > 1
          * a finite (non-nan) PSF magnitude
            - This separate check is largely to reject failed zeropoints.
          * and do not have flags set for bad, cosmic ray, edge or saturated
    - safeMatches: safe matches, as an afw.table.GroupView;
        safe matches are good matches that are sufficiently bright and sufficiently compact
    """

    # Filter down to matches with at least 2 sources and good flags
    flagKeys = [allMatches.schema.find("base_PixelFlags_flag_%s" % flag).key
                for flag in ("saturated", "cr", "bad", "edge")]
    nMatchesRequired = 2

    psfSnrKey = allMatches.schema.find("base_PsfFlux_snr").key
    psfMagKey = allMatches.schema.find("base_PsfFlux_mag").key
    psfMagErrKey = allMatches.schema.find("base_PsfFlux_magerr").key
    extendedKey = allMatches.schema.find("base_ClassificationExtendedness_value").key

    def goodFilter(cat, goodSnr=3):
        if len(cat) < nMatchesRequired:
            return False
        for flagKey in flagKeys:
            if cat.get(flagKey).any():
                return False
        if not np.isfinite(cat.get(psfMagKey)).all():
            return False
        psfSnr = np.median(cat.get(psfSnrKey))
        # Note that this also implicitly checks for psfSnr being non-nan.
        return psfSnr >= goodSnr

    goodMatches = allMatches.where(goodFilter)

    # Filter further to a limited range in S/N and extendedness
    # to select bright stars.
    safeMaxExtended = 1.0

    def safeFilter(cat):
        psfSnr = np.median(cat.get(psfSnrKey))
        extended = np.max(cat.get(extendedKey))
        return psfSnr >= safeSnr and extended < safeMaxExtended

    safeMatches = goodMatches.where(safeFilter)

    # Pass field=psfMagKey so np.mean just gets that as its input
    goodPsfSnr = goodMatches.aggregate(np.median, field=psfSnrKey)  # SNR
    goodPsfMag = goodMatches.aggregate(np.mean, field=psfMagKey)  # mag
    goodPsfMagRms = goodMatches.aggregate(np.std, field=psfMagKey)  # mag
    goodPsfMagErr = goodMatches.aggregate(np.median, field=psfMagErrKey)
    # positionRms knows how to query a group so we give it the whole thing
    #   by going with the default `field=None`.
    dist = goodMatches.aggregate(positionRms)

    return pipeBase.Struct(
        mag = goodPsfMag,
        magerr = goodPsfMagErr,
        magrms = goodPsfMagRms,
        snr = goodPsfSnr,
        dist = dist,
        goodMatches = goodMatches,
        safeMatches = safeMatches,
    )


def didThisRepoPass(repo, dataIds, configFile, **kwargs):
    """Convenience function for calling didIPass using the standard conventions for output filenames.

    Parameters
    ----------
    repo : str
        Path name of repository
    dataIds : list
        Data Ids that were analyzed
    configFile : str
        Configuration file with requirements specified as a dict.  E.g.,

        requirements: {'PA1': 25, 'PA2': 35}

    Returns
    -------
    bool
        Did all of the measured and required metrics pass.

    Raises
    ------
    AttributeError
        If the configuration file does not contain a `requirements` dict.

    See Also
    --------
    didIPass : The key function that does the work.
    """
    outputPrefix = repoNameToPrefix(repo)
    filters = set(d['filter'] for d in dataIds)
    try:
        requirements = loadParameters(configFile).requirements
    except AttributeError as ae:
        print("Configuration file %s does not contain a `requirements` dict." % configFile)
        raise(ae)

    return didIPass(outputPrefix, filters, requirements, **kwargs)


def didThisRepoPassSrd(repo, dataIds, level='design', **kwargs):
    """Convenience function for calling didIPass using the LSST SRD requirements.

    Parameters
    ----------
    repo : str
        Path name of repository
    dataIds : list
        Data Ids that were analyzed

    Returns
    -------
    bool
        Did all of the measured and required metrics pass.

    See Also
    --------
    didIPass : The key function that does the work.
    """
    outputPrefix = repoNameToPrefix(repo)
    filters = set(d['filter'] for d in dataIds)

    requirements = loadSrdRequirements(srdSpec, level=level)

    return didIPass(outputPrefix, filters, requirements, **kwargs)


def didIPass(*args, **kwargs):
    """Did this set pass.

    Returns
    -------
    bool
        Did all of the measured and requiremd metrics pass.
    """
    passedScores = scoreMetrics(*args, **kwargs)

    didAllPass = True
    for (metric, filter), passed in passedScores.iteritems():
        if not passed:
            print("Failed metric, filter: %s, %s" % (metric, filter))
            didAllPass = False

    return didAllPass


def scoreMetrics(outputPrefix, filters, requirements, verbose=False):
    """Score Key Performance metrics.  Returns dict((metric, filter), Pass/Fail)

    Parameters
    ----------
    outputPrefix : str
        The starting name for the output JSON files with the results
    filters : list, str, or None
        The filters in the analysis.  Output JSON files will be searched as
            "%s%s" % (outputPrefix, filters[i])
        If `None`, then JSON files will be searched for as just
            "%s" % outputPrefix.
    requirements : dict
        The requirements on each of the Key Performance Metrics
        Skips measurements for any metric without an entry in `requirements`.

    Returns
    -------
    dict of (str, str) -> bool
        A dictionary of results.  (metricName, filter) : True/False


    We provide the ability to check against configured standards
    instead of just the srdSpec because
    1. Different data sets may not lend themselves to satisfying the SRD.
    2. The pipeline continues to improve.
       Specifying a set of standards and updating that allows for a natural tightening of requirements.

    Note that there is no filter-dependence for the requirements.
    """
    if isinstance(filters, str):
        filters = list(filters)

    fileSnippet = dict(
        zip(
            ("PA1", "PF1", "PA2", "AM1", "AF1", "AM2", "AF2", "AM3", "AF3"),
            ("PA1", "PA2", "PA2", "AM1", "AM1", "AM2", "AM2", "AM3", "AM3")
        )
    )
    lookupKeyName = dict(
        zip(
            ("PA1", "PF1", "PA2", "AM1", "AF1", "AM2", "AF2", "AM3", "AF3"),
            ("PA1", "PF1", "PA2", "AMx", "AFx", "AMx", "AFx", "AMx", "AFx")
        )
    )
    metricsToConsider = ("PA1", "PF1", "PA2",
                         "AM1", "AF1", "AM2", "AF2", "AM3", "AF3")

    if verbose:
        print("{:16s}   {:13s} {:20s}".format("Measured", "Required", "Passes"))

    passed = {}
    for f in filters:
        if f:
            thisPrefix = "%s%s_" % (outputPrefix, f)
        else:
            thisPrefix = outputPrefix
        # get json files
        # Multiple metrics are sometimes stored in a file.
        # The names in those files may be generic ("AMx" instead of "AM1")
        # so we have three different, almost identical tuples here.
        for metricName in metricsToConsider:
            jsonFile = "%s%s.%s" % (thisPrefix, fileSnippet[metricName], 'json')

            metricNameKey = lookupKeyName[metricName]

            metricUnitsKey = metricNameKey.lower()+'Units'

            try:
                metricResults = loadKpmFromJson(jsonFile).getDict()
            except IOError:
                print("No results available for %s" % metricName)
                continue

            if metricName not in requirements:
                if verbose:
                    print("No requirement specified for %s.  Skipping." % metricName)
                continue

            # Check values against configured standards
            passed[(metricName, f)] = metricResults[metricNameKey] <= requirements[metricName]

            if verbose:
                kpmInfoToPrint = {
                    "name": metricName,
                    "value": metricResults[metricNameKey],
                    "units": metricResults[metricUnitsKey],
                    "spec": requirements[metricName],
                    "result": passed[(metricName, f)],
                }
                kpmInfoFormat = "{name:4s}: {value:5.2f} {units:4s} < {spec:5.2f} {units:4s} == {result}"
                print(kpmInfoFormat.format(**kpmInfoToPrint))

    return passed


####
def run(repo, dataIds, outputPrefix=None, level="design", verbose=False, **kwargs):
    """Main executable.

    Runs multiple filters, if necessary, through repeated calls to `runOneFilter`.
    Assesses results against SRD specs at specified `level`.

    Inputs
    ------
    repo : string
        The repository.  This is generally the directory on disk
        that contains the repository and mapper.
    dataIds : list of dict
        List of `butler` data IDs of Image catalogs to compare to reference.
        The `calexp` cpixel image is needed for the photometric calibration.
    outputPrefix : str, optional
        Specify the beginning filename for output files.
        The name of each filter will be appended to outputPrefix.
    level : str
        The level of the specification to check: "design", "minimum", "stretch"
    verbose : bool
        Provide detailed output.

    Outputs
    -------
    Names of plot files or JSON file are generated based on repository name,
    unless overriden by specifying `ouputPrefix`.
    E.g., Analyzing a repository "CFHT/output"
        will result in filenames that start with "CFHT_output_".
    The filter name is added to this prefix.  If the filter name has spaces,
        there will be annoyance and sadness as those spaces will appear in the filenames.
    """

    allFilters = set([d['filter'] for d in dataIds])

    if outputPrefix is None:
        outputPrefix = repoNameToPrefix(repo)

    for filt in allFilters:
        # Do this here so that each outputPrefix will have a different name for each filter.
        thisOutputPrefix = "%s_%s_" % (outputPrefix.rstrip('_'), filt)
        theseVisitDataIds = [v for v in dataIds if v['filter'] == filt]
        runOneFilter(repo, theseVisitDataIds, outputPrefix=thisOutputPrefix, verbose=verbose, filterName=filt,
                     **kwargs)

    # if verbose:
    #     print("==============================")
    #     print("Comparison against *LSST SRD*.")

    # SRDrequirements = {}
    # for k, v in srdSpec.getDict().iteritems():
    #     if isinstance(v, dict):
    #         SRDrequirements[k] = v[level]
    #     else:
    #         SRDrequirements[k] = v

    # scoreMetrics(outputPrefix, allFilters, SRDrequirements, verbose=verbose)


def runOneFilter(repo, visitDataIds, brightSnr=100,
                 medianAstromscatterRef=25, medianPhotoscatterRef=25, matchRef=500,
                 makePrint=True, makePlot=True, makeJson=True,
                 filterName=None, outputPrefix=None,
                 verbose=False,
                 **kwargs):
    """Main executable for the case where there is just one filter.

    Plot files and JSON files are generated in the local directory
        prefixed with the repository name (where '_' replace path separators),
    unless overriden by specifying `outputPrefix`.
    E.g., Analyzing a repository "CFHT/output"
        will result in filenames that start with "CFHT_output_".

    Parameters
    ----------
    repo : string
        The repository.  This is generally the directory on disk
        that contains the repository and mapper.
    dataIds : list of dict
        List of `butler` data IDs of Image catalogs to compare to reference.
        The `calexp` cpixel image is needed for the photometric calibration.
    brightSnr : float, optional
        Minimum SNR for a star to be considered bright
    medianAstromscatterRef : float, optional
        Expected astrometric RMS [mas] across visits.
    medianPhotoscatterRef : float, optional
        Expected photometric RMS [mmag] across visits.
    matchRef : int, optional
        Expectation of the number of stars that should be matched across visits.
    makePrint : bool, optional
        Print calculated quantities (to stdout).
    makePlot : bool, optional
        Create plots for metrics.  Saved to current working directory.
    makeJson : bool, optional
        Create JSON output file for metrics.  Saved to current working directory.
    outputPrefix : str, optional
        Specify the beginning filename for output files.
    filterName : str, optional
        Name of the filter (bandpass).
    verbose : bool, optional
        Output additional information on the analysis steps.

    """
    # Cache the YAML definitions of metrics (optional)
    yamlPath = os.path.join(getPackageDir('validate_drp'),
                            'metrics.yaml')
    with open(yamlPath) as f:
        yamlDoc = yaml.load(f)

    if outputPrefix is None:
        outputPrefix = repoNameToPrefix(repo)

    matchedDataset = MatchedMultiVisitDataset(repo, visitDataIds,
                                              verbose=verbose)
    photomModel = AnalyticPhotometryModel(matchedDataset)
    astromModel = AnalyticAstrometryModel(matchedDataset)

    measurements = []

    PA1 = PA1Measurement(matchedDataset, bandpass=filterName, verbose=verbose,
                         linkedBlobs=[photomModel, astromModel])
    measurements.append(PA1)
    print('filterName', filterName, PA1.bandpass)
    print('design', PA1.checkSpec('design'), filterName)

    pa2Metric = Metric.fromYaml('PA2', yamlDoc=yamlDoc)
    for specName in pa2Metric.getSpecNames(bandpass=filterName):
        PA2 = PA2Measurement(matchedDataset, bandpass=filterName,
                             specName=specName, verbose=verbose)
        measurements.append(PA2)
        print('PA2', specName, PA2.checkSpec(specName))

    pf1Metric = Metric.fromYaml('PF1', yamlDoc=yamlDoc)
    for specName in pf1Metric.getSpecNames(bandpass=filterName):
        PF1 = PF1Measurement(matchedDataset, bandpass=filterName,
                             specName=specName, verbose=verbose)
        measurements.append(PF1)
        print('PF1', specName, PF1.checkSpec(specName))

    AM1 = AMxMeasurement(1, matchedDataset,
                         bandpass=filterName, verbose=verbose,
                         linkedBlobs=[photomModel, astromModel])
    measurements.append(AM1)

    AM2 = AMxMeasurement(2, matchedDataset,
                         bandpass=filterName, verbose=verbose,
                         linkedBlobs=[photomModel, astromModel])
    measurements.append(AM2)

    AM3 = AMxMeasurement(3, matchedDataset,
                         bandpass=filterName, verbose=verbose,
                         linkedBlobs=[photomModel, astromModel])
    measurements.append(AM3)

    job_serializer = JobSerializer(
        measurements=measurements,
        blobs=[matchedDataset, photomModel, astromModel])
    persist_job(job_serializer, outputPrefix.rstrip('_') + '.json')

    # if makePlot:
    #     plotAstrometry(dist, magavg, struct.snr,
    #                    fit_params=astromStruct.params,
    #                    brightSnr=brightSnr, outputPrefix=outputPrefix)
    #     plotPhotometry(magavg, struct.snr, mmagerr, mmagrms,
    #                    fit_params=photStruct.params,
    #                    brightSnr=brightSnr, filterName=filterName, outputPrefix=outputPrefix)

    # magKey = allMatches.schema.find("base_PsfFlux_mag").key

    # AM1, AM2, AM3 = [calcOrNone(func, safeMatches, ValidateErrorNoStars, verbose=verbose)
    #                  for func in (calcAM1, calcAM2, calcAM3)]
    # PA1, PA2 = [func(safeMatches, magKey, verbose=verbose) for func in (calcPA1, calcPA2)]

    # blob = MultiVisitStarBlobSerializer.init_from_structs(
    #     filterName, struct, astromStruct, photStruct)
    # json.dumps(blob.json)

    # measurement_serializers = []

    # # Serialize AMx, AFx, ADx
    # for AMx in (AM1, AM2, AM3):
    #     if AMx is None:
    #         continue
    #     x = AMx.x

    #     AMx_serializer = MeasurementSerializer(
    #         metric=calcSrd.AMxSerializer(x=x),
    #         value=DatumSerializer(
    #             value=AMx.AMx,
    #             units='milliarcsecond',
    #             label='AM{0:d}'.format(x),
    #             description='Median RMS of the astrometric distance '
    #                         'distribution for stellar pairs with separation '
    #                         'of D arcmin (repeatability)'),
    #         parameters=calcSrd.AMxParamSerializer(
    #             D=DatumSerializer(
    #                 value=AMx.D,
    #                 units='arcmin',
    #                 label='D',
    #                 description='Fiducial distance between two objects to '
    #                             'consider'),
    #             annulus=DatumSerializer(
    #                 value=AMx.annulus,
    #                 units='arcmin',
    #                 label='annulus',
    #                 description='Annulus for selecting pairs of stars'),
    #             mag_range=DatumSerializer(
    #                 value=AMx.magRange,
    #                 units='mag',
    #                 label='mag range',
    #                 description='(bright, faint) magnitude selection range')),
    #         blob_id=blob.id)
    #     measurement_serializers.append(AMx_serializer)

    #     # So... only one spec level is computed???
    #     AFx_serializer = MeasurementSerializer(
    #         metric=calcSrd.AFxSerializer(x=x, level=AMx.level),
    #         value=DatumSerializer(
    #             value=AMx.AFx,
    #             units='',
    #             label='AF{0:d}'.format(x),
    #             description='Fraction of pairs that deviate more than AD{0:d} '
    #                         'from median AM{0:d} ({1})'.format(x, AMx.level)),
    #         parameters=calcSrd.AFxParamSerializer(
    #             ADx=DatumSerializer(
    #                 value=AMx.ADx_spec,
    #                 units='milliarcsecond',
    #                 label='AD{0:d}'.format(x),
    #                 description='Deviation from median RMS AM{0:d} '
    #                             'containing AF{0:d} of sample'.format(x)),
    #             AMx=DatumSerializer(
    #                 value=AMx.AMx,
    #                 units='milliarcsecond',
    #                 label='AM{0:d}'.format(x),
    #                 description='Median RMS of the astrometric distance '
    #                             'distribution for stellar pairs with '
    #                             'separation of D arcmin (repeatability)'),
    #             D=DatumSerializer(
    #                 value=AMx.D,
    #                 units='arcmin',
    #                 label='D',
    #                 description='Fiducial distance between two objects to '
    #                             'consider'),
    #             annulus=DatumSerializer(
    #                 value=AMx.annulus,
    #                 units='arcmin',
    #                 label='annulus',
    #                 description='Annulus for selecting pairs of stars'),
    #             mag_range=DatumSerializer(
    #                 value=AMx.magRange,
    #                 units='mag',
    #                 label='mag range',
    #                 description='(bright, faint) magnitude selection range')),
    #         blob_id=blob.id)
    #     measurement_serializers.append(AFx_serializer)

    # # Serialize PA1
    # PA1_serializer = MeasurementSerializer(
    #     metric=calcSrd.PA1Serializer(),
    #     value=DatumSerializer(
    #         value=PA1.rms,
    #         units='millimag',
    #         label='PA1',
    #         description='Median RMS of visit-to-visit relative photometry. '
    #                     'LPM-17.'),
    #     parameters=calcSrd.PA1ParamSerializer(
    #         num_random_shuffles=50),
    #     blob_id=blob.id)
    # json.dumps(PA1_serializer.json)
    # measurement_serializers.append(PA1_serializer)
    # # FIXME need to include the rest of PA1's measurement struct in a blob

    # # Serialize PA2 with each level of PF1
    # for level in srdSpec.levels:
    #     PA2_serializer = MeasurementSerializer(
    #         metric=calcSrd.PA2Serializer(spec_level=level),
    #         value=DatumSerializer(
    #             value=PA2.PA2_measured[level],
    #             units='millimag',
    #             label='PA2',
    #             description='Mags from mean relative photometric RMS that '
    #                         'encompasses PF1 of measurements.'),
    #         parameters=calcSrd.PA2ParamSerializer(
    #             num_random_shuffles=50,  # FIXME
    #             PF1=DatumSerializer(
    #                 value=PA2.PF1_spec[level],  # input for PA2
    #                 units='',
    #                 label='PF1',
    #                 description='Fraction of measurements between PA1 and '
    #                             'PF2, {0} spec'.format(level))),
    #         blob_id=blob.id)
    #     json.dumps(PA2_serializer.json)
    #     measurement_serializers.append(PA2_serializer)

    # # Serialize PF1 with each level of PA2
    # for level in srdSpec.levels:
    #     PF1_serializer = MeasurementSerializer(
    #         metric=calcSrd.PF1Serializer(spec_level=level),
    #         value=DatumSerializer(
    #             value=PA2.PF1_measured[level],
    #             units='',
    #             label='PF1',
    #             description='Fraction of measurements between PA1 and PF2, '
    #                         '{0} spec'.format(level)),
    #         parameters=calcSrd.PF1ParamSerializer(
    #             PA2=DatumSerializer(
    #                 value=PA2.PA2_spec[level],
    #                 units='millimag',
    #                 label='PA2',
    #                 description='Mags from mean relative photometric RMS that '
    #                             'encompasses PF1 of measurements at '
    #                             '{0} spec'.format(level))),
    #         blob_id=blob.id)
    #     json.dumps(PF1_serializer.json)
    #     measurement_serializers.append(PF1_serializer)

    # # Wrap measurements in a Job
    # job_serializer = JobSerializer(
    #     measurements=measurement_serializers,
    #     blobs=[blob])
    # persist_job(job_serializer, outputPrefix.rstrip('_') + '.json')

    # if makePrint:
    #     print("=============================================")
    #     print("Detailed comparison against SRD requirements.")
    #     print("The LSST SRD is at:  http://ls.st/LPM-17")
    #     printPA1(PA1)
    #     printPA2(PA2)
    #     for metric in (AM1, AM2, AM3):
    #         if metric:
    #             print("--")
    #             printAMx(metric)

    # if makePlot:
    #     plotPA1(PA1, outputPrefix=outputPrefix)
    #     for metric in (AM1, AM2, AM3):
    #         if metric:
    #             plotAMx(metric, outputPrefix=outputPrefix)

    # if makeJson:
    #     for struct in (astromStruct, photStruct):
    #         outfile = outputPrefix + "%s.json" % struct.name
    #         saveKpmToJson(struct, outfile)

    #     for metric in (AM1, AM2, AM3, PA1, PA2):
    #         if metric:
    #             outfile = outputPrefix + "%s.json" % metric.name
    #             saveKpmToJson(metric, outfile)
