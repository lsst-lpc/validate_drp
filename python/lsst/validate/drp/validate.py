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

import json
import numpy as np

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.image.utils as afwImageUtils
from lsst.afw.table import SourceCatalog, SchemaMapper, Field
from lsst.afw.table import MultiMatch, SimpleRecord, GroupView
from lsst.afw.fits.fitsLib import FitsError
import lsst.daf.persistence as dafPersist
import lsst.pipe.base as pipeBase
import pickle
import shelve

def sourceFlux(sourceFluxField = 'base_PsfFlux'): # 'base_CircularApertureFlux_6_0'):#'base_PsfFlux'):#
    return sourceFluxField


from .base import ValidateErrorNoStars
from .calcSrd import calcAM1, calcAM2, calcAM3, calcPA1, calcPA2
from . import calcSrd
from .check import checkAstrometry, checkPhotometry, positionRms
from .plot import plotAstrometry, plotPhotometry, plotPA1, plotAMx,  plotVisitVsTime, plotAstromPhotRMSvsTimeCcd
from .print import printPA1, printPA2, printAMx
from .srdSpec import srdSpec, loadSrdRequirements
from .util import getCcdKeyName, repoNameToPrefix, calcOrNone, loadParameters
from .io import (saveKpmToJson, loadKpmFromJson, MultiVisitStarBlobSerializer,
                 MeasurementSerializer, DatumSerializer, JobSerializer,
                 persist_job)

import lsst.pex.config as pexConfig

sourceFluxField = sourceFlux()


sourceFluxFields = ['base_PsfFlux', 'base_CircularApertureFlux_6_0']
commentsFluxFields = ['PSF','CircularAperture 6_0']



def loadAndMatchData(repo, dataIds,
                     matchRadius=afwGeom.Angle(1, afwGeom.arcseconds),
                     verbose=False, MetaData=False):
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
    # mapper.addOutputField(Field[float]('base_PsfFlux_snr', "PSF flux SNR"))
    # mapper.addOutputField(Field[float]('base_PsfFlux_mag', "PSF magnitude"))
    # mapper.addOutputField(Field[float]('base_PsfFlux_magerr', "PSF magnitude uncertainty"))
    for i in range(len(sourceFluxFields)): ### sourceflux
        mapper.addOutputField(Field[float](sourceFluxFields[i]+'_snr', commentsFluxFields[i]+" flux SNR"))
        mapper.addOutputField(Field[float](sourceFluxFields[i]+'_mag', commentsFluxFields[i]+" magnitude"))
        mapper.addOutputField(Field[float](sourceFluxFields[i]+'_magerr', commentsFluxFields[i]+" magnitude uncertainty"))
    mapper.addOutputField(Field[float]('MJD-OBS', "Date observation (mjd)"))
    mapper.addOutputField(Field[float]('FLUXMAG0', "zero point flux"))
    mapper.addOutputField(Field[float]('FLUXMAG0ERR', "zero point flux error" ))
    mapper.addOutputField(Field[float]('PSF-FWHM', "Full width at half maximum (FWHM) of the psf (arcseconds)"))
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
            # print('dataIds vId',vId)
            calexp  =  butler.get('calexp',  vId, immediate=True)
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

        fluxmag0 = calexpMetadata.get('FLUXMAG0')
        fluxmag0_err = calexpMetadata.get('FLUXMAG0ERR')
      
        tmpCat['MJD-OBS'][:] = calexpMetadata.get('MJD-OBS')
        tmpCat['FLUXMAG0'][:] = fluxmag0
        tmpCat['FLUXMAG0ERR'][:] = fluxmag0_err

        sigmaToFwhm  =  2.0*np.sqrt(2.0*np.log(2.0))
 
        wcs = calexp.getWcs()
        psf = calexp.getPsf() # make a setup meas_extensions_psfex to import this function
        shape = psf.computeShape()
        psf_fwhm = shape.getDeterminantRadius()* wcs.pixelScale().asArcseconds()  *  sigmaToFwhm
        # print('PSF-FWHM :',   psf_fwhm)
        tmpCat['PSF-FWHM'][:] = psf_fwhm

        # tmpCat['base_PsfFlux_snr'][:] = tmpCat['base_PsfFlux_flux'] / tmpCat['base_PsfFlux_fluxSigma']
        # with afwImageUtils.CalibNoThrow():
        #     (tmpCat['base_PsfFlux_mag'][:], tmpCat['base_PsfFlux_magerr'][:]) = \
        #         calib.getMagnitude(tmpCat['base_PsfFlux_flux'],
        #                            tmpCat['base_PsfFlux_fluxSigma'])
        for sourceFluxField in sourceFluxFields:
            tmpCat[sourceFluxField+'_snr'][:] = tmpCat[sourceFluxField+'_flux'] / tmpCat[sourceFluxField+'_fluxSigma']
            with afwImageUtils.CalibNoThrow():
                (tmpCat[sourceFluxField+'_mag'][:], tmpCat[sourceFluxField+'_magerr'][:]) = \
                    calib.getMagnitude(tmpCat[sourceFluxField+'_flux'],
                                       tmpCat[sourceFluxField+'_fluxSigma'])
        srcVis.extend(tmpCat, False)
        mmatch.add(catalog=tmpCat, dataId=vId)

    # Complete the match, returning a catalog that includes
    # all matched sources with object IDs that can be used to group them.
    matchCat = mmatch.finish()

    # matchCat_schema =  matchCat.getSchema()
    # print('schema_matchCat.getOrderedNames()',matchCat_schema.getOrderedNames())

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

 #   flagKeys = [allMatches.schema.find("base_PixelFlags_flag_%s" % flag).key 
 #               for flag in ("saturated", "cr", "bad", "edge")] #anciens flags
 #   nMatchesRequired = 2

    flagKeys = [allMatches.schema.find("base_%s" % flag).key
                for flag in ("PixelFlags_flag_saturated", "PixelFlags_flag_cr", "PixelFlags_flag_bad", "PixelFlags_flag_edge", "SdssCentroid_flag", "SdssShape_flag")] # considering flags from jointcal
    nMatchesRequired = 2 


    psfSnrKey = allMatches.schema.find(sourceFluxField+"_snr").key  # il faudrait renommer les "psf*" en autre chose, vu que cela depend maintenant de sourceFluxField
    psfMagKey = allMatches.schema.find(sourceFluxField+"_mag").key
    psfMagErrKey = allMatches.schema.find(sourceFluxField+"_magerr").key
    # ###
    centroid = "base_SdssCentroid"
    shape = "base_SdssShape"

    fluxKey = allMatches.schema.find(sourceFluxField+"_flux").key
    fluxErrKey = allMatches.schema.find(sourceFluxField+"_fluxSigma").key
    parentKey = allMatches.schema.find("parent").key
    # ##
 #   schemaAllMatches=allMatches.schema#
 #   print('=============================== schemaAllMatches',schemaAllMatches.getOrderedNames())

    # psfSnrKey = allMatches.schema.find("base_PsfFlux_snr").key
    # psfMagKey = allMatches.schema.find("base_PsfFlux_mag").key
    # psfMagErrKey = allMatches.schema.find("base_PsfFlux_magerr").key
    extendedKey = allMatches.schema.find("base_ClassificationExtendedness_value").key


    
    def goodFilter(cat, goodSnr=3):#, maxMag=22.5):
        if len(cat) < nMatchesRequired:
            return False
        for flagKey in flagKeys: # ne pas considerer les sources avec des bad flags
            if cat.get(flagKey).any():
                return False
        if not np.isfinite(cat.get(psfMagKey)).all():
            return False
        # ### Selections in JointCal
       # print('(dir(cat)',dir(cat))
        for src in cat:
            # print('(dir(src)',dir(src))
            # Reject negative flux
            flux = src.get(fluxKey)
            
            if flux < 0:
                print('flux negatif', flux)
                return False
            # Reject objects with too large magnitude
           # maxMag=22.5
           # mag = src.get(psfMagKey)
           # magErr = src.get(psfMagErrKey)
           # fluxErr = src.get(fluxErrKey)
            #if mag > maxMag or magErr > 0.1 or flux/fluxErr < 10: # 
          #  if magErr > 0.1 or flux/fluxErr < 10:
                #print('mag or magERR rejected')
            #    return False

            # cut snr 100 direct (pour le test dans /test_validation_phVisits_cutssnr100direct/ ) :
           # SNR = src.get(psfSnrKey)
           # if SNR<100:
           #     return False

            # Reject blends (?)
            if src.get(parentKey) != 0:
                print("parentKey Rejection")
                return False
         #   footprint = src.getFootprint() # ???
         #   if footprint is not None and len(footprint.getPeaks()) > 1:
         #       return False
          # Check consistency of variances and second moments
            vx = np.square(src.get(centroid + "_xSigma"))
            vy = np.square(src.get(centroid + "_ySigma"))
            mxx = src.get(shape + "_xx")
            myy = src.get(shape + "_yy")
            mxy = src.get(shape + "_xy")
            vxy = mxy*(vx+vy)/(mxx+myy)

            if vxy*vxy > vx*vy or np.isnan(vx) or np.isnan(vy):
                print('shapevxy*vxy > vx*vy or np.isnan(vx) or np.isnan(vy) rejection')
                return False

        # ## fin nouveaux ajouts
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
    goodPsfSnr = goodMatches.aggregate(np.median, field=psfSnrKey)
    goodPsfMag = goodMatches.aggregate(np.mean, field=psfMagKey)
    goodPsfMagRms = goodMatches.aggregate(np.std, field=psfMagKey)
    goodPsfMagErr = goodMatches.aggregate(np.median, field=psfMagErrKey)
    # positionRms knows how to query a group so we give it the whole thing
    #   by going with the default `field=None`.
    dist = goodMatches.aggregate(positionRms)

    # ajout d'un ndarray contenant toutes les sources
    print('======================================================================')

    data_titles=['visit', 'ccd', 'coord_ra', 'coord_dec', 'MJD-OBS', sourceFluxField+'_snr', sourceFluxField+'_flux', sourceFluxField+'_fluxSigma', sourceFluxField+'_mag', sourceFluxField+'_magerr', 'Nb_group','id', 'PSF-FWHM', 'FLUXMAG0', 'FLUXMAG0ERR',  'MeanGrpRa', 'MeanGrpDec', 'MedGrpSnr', 'base_ClassificationExtendedness_value']
    
    Dtype=[]
    nbsources=100

    for title in data_titles:
        if title =='visit' or title=='ccd' or title=='Nb_group' or title=='id':
            Dtype.append((title,np.int))
        else:
            Dtype.append((title,np.float))

    sources_ndarray = np.zeros(nbsources, dtype = Dtype)
 
    Nb_group=0
    compteur_sources=0
    
    for group in goodMatches.groups:
        snr = group.get(sourceFluxField+"_snr")
        RA = group.get('coord_ra')
        Dec = group.get('coord_dec')
        psf_fwhm = group.get('PSF-FWHM')
        MeanPsf_fwhm = np.mean(psf_fwhm)
        MeanRA = np.mean(RA)
        MeanDec = np.mean(Dec)
        FluxMag0 = group.get('FLUXMAG0')
        FluxMag0Err = group.get('FLUXMAG0ERR')
        visit= group.get('visit')
        ccd = group.get('ccd')
        MedSnr = np.median(snr)

        for src in group:# range(len(RA)):
           # print( i)
           # print(' RA[i]', RA[i])
            sources_ndarray['id'][compteur_sources] = src.get('id')
            sources_ndarray['visit'][compteur_sources] = src.get('visit')
            sources_ndarray['coord_ra'][compteur_sources] = src.get('coord_ra')
            sources_ndarray['coord_dec'][compteur_sources] = src.get('coord_dec')
            sources_ndarray['ccd'][compteur_sources] = src.get('ccd')
            sources_ndarray['MJD-OBS'][compteur_sources] = src.get('MJD-OBS')
            sources_ndarray['PSF-FWHM'][compteur_sources] = src.get('PSF-FWHM')
            sources_ndarray[sourceFluxField+'_snr'][compteur_sources] = src.get(sourceFluxField+'_snr')
            sources_ndarray[sourceFluxField+'_flux'][compteur_sources] = src.get(sourceFluxField+'_flux')
            sources_ndarray[sourceFluxField+'_fluxSigma'][compteur_sources] = src.get(sourceFluxField+'_fluxSigma')
            sources_ndarray[sourceFluxField+'_mag'][compteur_sources] = src.get(sourceFluxField+'_mag')
            sources_ndarray[sourceFluxField+'_magerr'][compteur_sources] = src.get(sourceFluxField+'_magerr')
            sources_ndarray['base_ClassificationExtendedness_value'][compteur_sources] = src.get('base_ClassificationExtendedness_value')
            sources_ndarray['Nb_group'][compteur_sources] =  Nb_group
            sources_ndarray['FLUXMAG0'][compteur_sources] = src.get('FLUXMAG0')
            sources_ndarray['FLUXMAG0ERR'][compteur_sources] = src.get('FLUXMAG0ERR')
            
            sources_ndarray['MeanGrpRa'][compteur_sources] =  MeanRA 
            sources_ndarray['MeanGrpDec'][compteur_sources] =  MeanDec
            sources_ndarray['MedGrpSnr'][compteur_sources] =  MedSnr 
            # sources_ndarray[][compteur_sources] = src.get()

            compteur_sources+=1
            if compteur_sources >= len(sources_ndarray):
                sources_ndarray = np.resize(sources_ndarray, compteur_sources+100)

        Nb_group+=1

    sources_ndarray = np.resize(sources_ndarray, compteur_sources)
    #print('Dtype', Dtype)
    #print('set(visits)',set(sources_ndarray['visit']))
    print('Dtype sources_ndarray.dtype', sources_ndarray.dtype)
    #print('shape ndarray', np.shape(sources_ndarray))

    print(' ndarray', sources_ndarray)


    # extraction du ndarray dans un fichier pkl
    filename='sources_ndarray_grp_'+str(len(set(sources_ndarray['visit'])))+'visits.pkl'
    output_file=open(filename, 'wb')
    pickle.dump(sources_ndarray, output_file)
    output_file.close()
    print('sources_ndarray sauved \n')
    print( '======================================================================')

    return pipeBase.Struct(
       # name="structgoodMatches", #ajout pour le json
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

    if verbose:
        print("==============================")
        print("Comparison against *LSST SRD*.")

    SRDrequirements = {}
    for k, v in srdSpec.getDict().iteritems():
        if isinstance(v, dict):
            SRDrequirements[k] = v[level]
        else:
            SRDrequirements[k] = v

    scoreMetrics(outputPrefix, allFilters, SRDrequirements, verbose=verbose)


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

    if outputPrefix is None:
        outputPrefix = repoNameToPrefix(repo)

    filterName = set([dId['filter'] for dId in visitDataIds]).pop()
    allMatches = loadAndMatchData(repo, visitDataIds, verbose=verbose)
   # allMatches, metaData = loadAndMatchData(repo, visitDataIds, verbose=verbose, MetaData=True)
    struct = analyzeData(allMatches, brightSnr,  verbose=verbose)

    magavg = struct.mag
    magerr = struct.magerr
    magrms = struct.magrms
    dist = struct.dist
    match = len(struct.goodMatches)
    safeMatches = struct.safeMatches

    mmagerr = 1000*magerr
    mmagrms = 1000*magrms

    astromStruct = \
        checkAstrometry(struct.snr, dist, match,
                        brightSnr=brightSnr,
                        medianRef=medianAstromscatterRef, matchRef=matchRef)
    photStruct = \
        checkPhotometry(struct.snr, magavg, mmagerr, mmagrms, dist, match,
                        brightSnr=brightSnr,
                        medianRef=medianPhotoscatterRef, matchRef=matchRef)
    if makePlot:
        plotAstrometry(dist, magavg, struct.snr,
                       fit_params=astromStruct.params,
                       brightSnr=brightSnr, outputPrefix=outputPrefix)
        
        plotVisitVsTime(struct.goodMatches,
                        outputPrefix=outputPrefix)
 
        plotAstromPhotRMSvsTimeCcd(dist, magavg, struct.snr, struct.goodMatches, mmagrms,
                                   fit_params=astromStruct.params, srcFluxField=sourceFluxField,
                                   brightSnr=brightSnr, outputPrefix=outputPrefix)

        plotAstromPhotRMSvsTimeCcd(dist, magavg, struct.snr, struct.goodMatches, mmagrms,
                                   fit_params=astromStruct.params, srcFluxField=sourceFluxField,
                                   brightSnr=brightSnr, outputPrefix=outputPrefix, zoom=True, dico=True)
        print('')
        print('astromStruct.params,',astromStruct.params)
        print('')
        print('photStruct.params',photStruct.params)
        print('')
        plotPhotometry(magavg, struct.snr, mmagerr, mmagrms,
                       fit_params=photStruct.params,
                       brightSnr=brightSnr, filterName=filterName, outputPrefix=outputPrefix)
   # if makeJson: #test
     #   outfile=open('structgoodMatches.pkl', 'w')

   #     structure= pipeBase.Struct(name="structgoodMatches",
                     #              model_name="photErrModel",
                      #             doc=photErrModel.__doc__,
                      #             params=fit_params,
                      #             photRmsScatter=photScatter)

        #outfile = outputPrefix + "%sAAAAA.json" % struct.name
        #saveKpmToJson(struct, outfile) # marche pas "is not JSON serializable"
     #   save_obj(struct, 'GOODMATCHES')  # marche pas, "TypeError: can't pickle SwigPyObject objects"

       # cPickle.dump(struct.goodMatches, outfile)
       ## pickle#.dump(struct, outfile)
       # outfile.close()
        #outfile = outputPrefix + "structgoodMatches.shelve"
        #saveKpmToJson(struct.goodMatches, outfile)
      #  shelve_goodMatches = shelve.open('structgoodMatches.shelve')
       # shelve_goodMatches['good'] = struct.goodMatches
      #  shelve_goodMatches.close()
        
    # magKey = allMatches.schema.find("base_PsfFlux_mag").key
    magKey = allMatches.schema.find(sourceFluxField+"_mag").key

    AM1, AM2, AM3 = [calcOrNone(func, safeMatches, ValidateErrorNoStars, verbose=verbose)
                     for func in (calcAM1, calcAM2, calcAM3)]
    PA1, PA2 = [func(safeMatches, magKey, verbose=verbose) for func in (calcPA1, calcPA2)]

    blob = MultiVisitStarBlobSerializer.init_from_structs(
        filterName, struct, astromStruct, photStruct)
    json.dumps(blob.json)

    measurement_serializers = []

    # Serialize AMx, AFx, ADx
    for AMx in (AM1, AM2, AM3):
        if AMx is None:
            continue
        x = AMx.x

        AMx_serializer = MeasurementSerializer(
            metric=calcSrd.AMxSerializer(x=x),
            value=DatumSerializer(
                value=AMx.AMx,
                units='milliarcsecond',
                label='AM{0:d}'.format(x),
                description='Median RMS of the astrometric distance '
                            'distribution for stellar pairs with separation '
                            'of D arcmin (repeatability)'),
            parameters=calcSrd.AMxParamSerializer(
                D=DatumSerializer(
                    value=AMx.D,
                    units='arcmin',
                    label='D',
                    description='Fiducial distance between two objects to '
                                'consider'),
                annulus=DatumSerializer(
                    value=AMx.annulus,
                    units='arcmin',
                    label='annulus',
                    description='Annulus for selecting pairs of stars'),
                mag_range=DatumSerializer(
                    value=AMx.magRange,
                    units='mag',
                    label='mag range',
                    description='(bright, faint) magnitude selection range')),
            blob_id=blob.id)
        measurement_serializers.append(AMx_serializer)

        # So... only one spec level is computed???
        AFx_serializer = MeasurementSerializer(
            metric=calcSrd.AFxSerializer(x=x, level=AMx.level),
            value=DatumSerializer(
                value=AMx.AFx,
                units='',
                label='AF{0:d}'.format(x),
                description='Fraction of pairs that deviate more than AD{0:d} '
                            'from median AM{0:d} ({1})'.format(x, AMx.level)),
            parameters=calcSrd.AFxParamSerializer(
                ADx=DatumSerializer(
                    value=AMx.ADx_spec,
                    units='milliarcsecond',
                    label='AD{0:d}'.format(x),
                    description='Deviation from median RMS AM{0:d} '
                                'containing AF{0:d} of sample'.format(x)),
                AMx=DatumSerializer(
                    value=AMx.AMx,
                    units='milliarcsecond',
                    label='AM{0:d}'.format(x),
                    description='Median RMS of the astrometric distance '
                                'distribution for stellar pairs with '
                                'separation of D arcmin (repeatability)'),
                D=DatumSerializer(
                    value=AMx.D,
                    units='arcmin',
                    label='D',
                    description='Fiducial distance between two objects to '
                                'consider'),
                annulus=DatumSerializer(
                    value=AMx.annulus,
                    units='arcmin',
                    label='annulus',
                    description='Annulus for selecting pairs of stars'),
                mag_range=DatumSerializer(
                    value=AMx.magRange,
                    units='mag',
                    label='mag range',
                    description='(bright, faint) magnitude selection range')),
            blob_id=blob.id)
        measurement_serializers.append(AFx_serializer)

    # Serialize PA1
    PA1_serializer = MeasurementSerializer(
        metric=calcSrd.PA1Serializer(),
        value=DatumSerializer(
            value=PA1.rms,
            units='millimag',
            label='PA1',
            description='Median RMS of visit-to-visit relative photometry. '
                        'LPM-17.'),
        parameters=calcSrd.PA1ParamSerializer(
            num_random_shuffles=50),
        blob_id=blob.id)
    json.dumps(PA1_serializer.json)
    measurement_serializers.append(PA1_serializer)
    # FIXME need to include the rest of PA1's measurement struct in a blob

    # Serialize PA2 with each level of PF1
    for level in srdSpec.levels:
        PA2_serializer = MeasurementSerializer(
            metric=calcSrd.PA2Serializer(spec_level=level),
            value=DatumSerializer(
                value=PA2.PA2_measured[level],
                units='millimag',
                label='PA2',
                description='Mags from mean relative photometric RMS that '
                            'encompasses PF1 of measurements.'),
            parameters=calcSrd.PA2ParamSerializer(
                num_random_shuffles=50,  # FIXME
                PF1=DatumSerializer(
                    value=PA2.PF1_spec[level],  # input for PA2
                    units='',
                    label='PF1',
                    description='Fraction of measurements between PA1 and '
                                'PF2, {0} spec'.format(level))),
            blob_id=blob.id)
        json.dumps(PA2_serializer.json)
        measurement_serializers.append(PA2_serializer)

    # Serialize PF1 with each level of PA2
    for level in srdSpec.levels:
        PF1_serializer = MeasurementSerializer(
            metric=calcSrd.PF1Serializer(spec_level=level),
            value=DatumSerializer(
                value=PA2.PF1_measured[level],
                units='',
                label='PF1',
                description='Fraction of measurements between PA1 and PF2, '
                            '{0} spec'.format(level)),
            parameters=calcSrd.PF1ParamSerializer(
                PA2=DatumSerializer(
                    value=PA2.PA2_spec[level],
                    units='millimag',
                    label='PA2',
                    description='Mags from mean relative photometric RMS that '
                                'encompasses PF1 of measurements at '
                                '{0} spec'.format(level))),
            blob_id=blob.id)
        json.dumps(PF1_serializer.json)
        measurement_serializers.append(PF1_serializer)

    # Wrap measurements in a Job
    job_serializer = JobSerializer(
        measurements=measurement_serializers,
        blobs=[blob])
    persist_job(job_serializer, outputPrefix.rstrip('_') + '.json')

    if makePrint:
        print("=============================================")
        print("Detailed comparison against SRD requirements.")
        print("The LSST SRD is at:  http://ls.st/LPM-17")
        printPA1(PA1)
        printPA2(PA2)
        for metric in (AM1, AM2, AM3):
            if metric:
                print("--")
                printAMx(metric)

    if makePlot:
        plotPA1(PA1, outputPrefix=outputPrefix)
        for metric in (AM1, AM2, AM3):
            if metric:
                plotAMx(metric, outputPrefix=outputPrefix)

    if makeJson:
        for struct in (astromStruct, photStruct):
            outfile = outputPrefix + "%s.json" % struct.name
            saveKpmToJson(struct, outfile)

        for metric in (AM1, AM2, AM3, PA1, PA2):
            if metric:
                outfile = outputPrefix + "%s.json" % metric.name
                saveKpmToJson(metric, outfile)













# Selections in JointCal (Pierre Astier):

#class StarSelectorConfig(pexConfig.Config):
#
#    badFlags = pexConfig.ListField(
#        doc = "List of flags which cause a source to be rejected as bad",
#        dtype = str,
#        default = ["base_PixelFlags_flag_saturated",
#                   "base_PixelFlags_flag_cr",
#                   "base_PixelFlags_flag_interpolated",
#                   "base_SdssCentroid_flag",
#                   "base_SdssShape_flag"],
#    )
#    sourceFluxField = pexConfig.Field(
#        doc = "Type of source flux",
#        dtype = str,
#        default = "slot_CalibFlux"
#    )
#    maxMag = pexConfig.Field(
#        doc = "Maximum magnitude for sources to be included in the fit",
#        dtype = float,
#        default = 22.5,
#    )
#    coaddName = pexConfig.Field(
#        doc = "Type of coadd",
#        dtype = str,
#        default = "deep"
#    )
#    centroid = pexConfig.Field(
#        doc = "Centroid type for position estimation",
#        dtype = str,
#        default = "base_SdssCentroid",
#    )
#    shape = pexConfig.Field(
#        doc = "Shape for error estimation",
#        dtype = str,
#        default = "base_SdssShape",
#    )
#
#
#class StarSelector(object):
#
#    ConfigClass = StarSelectorConfig
#
#    def __init__(self, config):
#        """Construct a star selector
#        @param[in] config: An instance of StarSelectorConfig
#        """
#        self.config = config
#
#    def select(self, srcCat, calib):
#        """Return a catalog containing only reasonnable stars / galaxies."""
#
#        schema = srcCat.getSchema()
#        newCat = afwTable.SourceCatalog(schema)
#        fluxKey = schema[self.config.sourceFluxField+"_flux"].asKey()
#        fluxErrKey = schema[self.config.sourceFluxField+"_fluxSigma"].asKey()
#        parentKey = schema["parent"].asKey()
#        flagKeys = []
#        for f in self.config.badFlags:
#            key = schema[f].asKey()
#            flagKeys.append(key)
#        fluxFlagKey = schema[self.config.sourceFluxField+"_flag"].asKey()
#        flagKeys.append(fluxFlagKey)
#
#        for src in srcCat:
#            # Do not consider sources with bad flags
#            for f in flagKeys:
#                rej = 0
#                if src.get(f):
#                    rej = 1
#                    break
#            if rej == 1:
#                continue
#            # Reject negative flux
#            flux = src.get(fluxKey)
#            if flux < 0:
#                continue
#            # Reject objects with too large magnitude
#            fluxErr = src.get(fluxErrKey)
#            mag, magErr = calib.getMagnitude(flux, fluxErr)
#            if mag > self.config.maxMag or magErr > 0.1 or flux/fluxErr < 10:
#                continue
#            # Reject blends
#            if src.get(parentKey) != 0:
#                continue
#            footprint = src.getFootprint()
#            if footprint is not None and len(footprint.getPeaks()) > 1:
#                continue
#
#            # Check consistency of variances and second moments
#            vx = np.square(src.get(self.config.centroid + "_xSigma"))
#            vy = np.square(src.get(self.config.centroid + "_ySigma"))
#            mxx = src.get(self.config.shape + "_xx")
#            myy = src.get(self.config.shape + "_yy")
#            mxy = src.get(self.config.shape + "_xy")
#            vxy = mxy*(vx+vy)/(mxx+myy)
#
#            if vxy*vxy > vx*vy or np.isnan(vx) or np.isnan(vy):
#                continue
#
#            newCat.append(src)
#
#        return newCat
