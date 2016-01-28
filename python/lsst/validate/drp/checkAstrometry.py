#!/usr/bin/env python

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

from __future__ import print_function, division

import os.path
import sys

import matplotlib.pylab as plt
import numpy as np

import lsst.daf.persistence as dafPersist
import lsst.pipe.base as pipeBase
from lsst.afw.table import matchRaDec, SimpleCatalog, SourceCatalog, SchemaMapper, Field
from lsst.afw.table import MultiMatch, SimpleRecord, GroupView

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.coord as afwCoord


def getCcdKeyName(dataid):
    """Return the key in a dataId that's referring to the CCD or moral equivalent.

    @param dataid  A dictionary that will be searched for a key that matches
        an entry in the hardcoded list of possible names for the CCD field.

    Different camera mappings use different keys
      to indicate the different amps/ccds in the same exposure.
    We here look through the reference dataId to determine which one this.
    """
    possibleCcdFieldNames = ['ccd', 'ccdnum', 'camcol']

    for name in possibleCcdFieldNames:
        if name in dataid:
            return name
    else:
        return None


def isExtended(source, extendedKey, extendedThreshold=1.0):
    """Is the source extended attribute above the threshold.

    Higher values of extendedness indicate a resolved object
    that is larger than a point source.
    """
    return source.get(extendedKey) >= extendedThreshold

def averageRaDec(cat):
    """Calculate the RMS for RA, Dec for a set of observations an object."""
    ra = np.mean(cat.get('coord_ra'))
    dec = np.mean(cat.get('coord_dec'))
    return ra, dec

# Some thoughts from Paul Price on how to do the coordinate differences correctly:
#    mean = sum([afwGeom.Extent3D(coord.toVector()) 
#                for coord in coordList, afwGeom.Point3D(0, 0, 0)])
#    mean /= len(coordList)
#    mean = afwCoord.IcrsCoord(mean)

def positionDiff(cat):
    """Calculate the diff RA, Dec from the mean for a set of observations an object for each observation.

    @param[in]  cat -- Collection with a .get method 
         for 'coord_ra', 'coord_dec' that returns radians.

    @param[out]  pos_median -- median diff of positions in milliarcsec.  Float.

    This is WRONG!
    Doesn't do wrap-around
    """
    ra_avg, dec_avg = averageRaDec(cat)
    ra, dec = cat.get('coord_ra'), cat.get('coord_dec')
    # Approximating that the cos(dec) term is roughly the same
    #   for all observations of this object.
    ra_diff  = (ra - ra_avg) * np.cos(dec)**2
    dec_diff = dec - dec_avg
    pos_diff = np.sqrt(ra_diff**2 + dec_diff**2)  # radians
    pos_diff = [afwGeom.radToMas(p) for p in pos_diff]  # milliarcsec

    return pos_diff

def positionRms(cat):
    """Calculate the RMS for RA, Dec for a set of observations an object.

    @param[in]  cat -- Collection with a .get method 
         for 'coord_ra', 'coord_dec' that returns radians.

    @param[out]  pos_rms -- RMS of positions in milliarcsecond.  Float.

    This is WRONG!
    Doesn't do wrap-around
    """
    ra_avg, dec_avg = averageRaDec(cat)
    ra, dec = cat.get('coord_ra'), cat.get('coord_dec')
    # Approximating that the cos(dec) term is roughly the same
    #   for all observations of this object.
    ra_var = np.var(ra) * np.cos(dec_avg)**2
    dec_var = np.var(dec)
    pos_rms = np.sqrt(ra_var + dec_var)  # radians
    pos_rms = afwGeom.radToMas(pos_rms)  # milliarcsec

    return pos_rms


def loadAndMatchData(repo, visitDataIds,
                     matchRadius=afwGeom.Angle(1, afwGeom.arcseconds)):
    """Load data from specific visit.  Match with reference.

    @param repo  The repository.  This is generally the directory on disk
                    that contains the repository and mapper.
    @param visitDataIds   Butler Data ID of Image catalogs to compare to reference.
           The actual pixel image is also needed for now for the photometric calibration.
           List.

    @param matchRadius    Radius for matching.  afwGeom.Angle().

    Return a pipeBase.Struct with mag, dist, and number of matches.
    """

    # Following 
    # https://github.com/lsst/afw/blob/tickets/DM-3896/examples/repeatability.ipynb
    butler = dafPersist.Butler(repo)
    dataset = 'src'

    dataRefs = [dafPersist.ButlerDataRef(butler, vId) for vId in visitDataIds]

    # retrieve the schema of the source catalog and extend it in order to add a field to record the ccd number
    ccdKeyName = getCcdKeyName(visitDataIds[0])

    schema = butler.get(dataset + "_schema", immediate=True).schema
    mapper = SchemaMapper(schema)
    mapper.addMinimalSchema(schema)
#    mapper.addOutputField(Field[int](ccdKeyName, "CCD number"))
    mapper.addOutputField(Field[float]('base_PsfFlux_mag', "PSF magnitude"))
    mapper.addOutputField(Field[float]('base_PsfFlux_magerr', "PSF magnitude uncertainty"))
    newSchema = mapper.getOutputSchema()

    # Create an object that can match multiple catalogs with the same schema
    mmatch = MultiMatch(newSchema, 
                        dataIdFormat = {'visit': int, ccdKeyName: int},
                        radius=matchRadius,
                        RecordClass=SimpleRecord)

    if False:
        # Create visit catalogs by merging those from constituent CCDs. We also convert from BaseCatalog to SimpleCatalog, but in the future we'll probably want to have the source transformation tasks generate SimpleCatalogs (or SourceCatalogs) directly.
        byVisit = {}
        for dataRef in dataRefs:
            catalog = byVisit.setdefault(dataRef.dataId["visit"], SimpleRecord.Catalog(schema))
            catalog.extend(dataRef.get(dataset, immediate=True), deep=True)

        # Loop over visits, adding them to the match.
        for visit, catalog in byVisit.iteritems():
            mmatch.add(catalog, dict(visit=visit))

    # create the new extented source catalog
    srcVis = SourceCatalog(newSchema)

    for vId in visitDataIds:
        calib = afwImage.Calib(butler.get("calexp_md", vId, immediate=True))
        calib.setThrowOnNegativeFlux(False)
        oldSrc = butler.get('src', vId, immediate=True)
        print(len(oldSrc), "sources in ccd: ", vId[ccdKeyName])

        # create temporary catalog
        tmpCat = SourceCatalog(SourceCatalog(newSchema).table)
        tmpCat.extend(oldSrc, mapper=mapper)
        (tmpCat['base_PsfFlux_mag'][:],
         tmpCat['base_PsfFlux_magerr'][:]) = \
             calib.getMagnitude(tmpCat['base_PsfFlux_flux'],
                                tmpCat['base_PsfFlux_fluxSigma'])
        srcVis.extend(tmpCat, False)
        mmatch.add(catalog=tmpCat, dataId=vId)

    # Complete the match, returning a catalog that includes all matched sources with object IDs that can be used to group them.
    matchCat = mmatch.finish()

    # Create a mapping object that allows the matches to be manipulated as a mapping of object ID to catalog of sources.
    allMatches = GroupView.build(matchCat)
    
    # Filter down to matches with at least 2 sources and good flags
    flagKeys = [allMatches.schema.find("base_PixelFlags_flag_%s" % flag).key
                for flag in ("saturated", "cr", "bad", "edge")]
    nMatchesRequired = 2


    def goodFilter(cat):
        if len(cat) < nMatchesRequired:
            return False
        for flagKey in flagKeys:
            if cat.get(flagKey).any():
                return False
        return True

    goodMatches = allMatches.where(goodFilter)

    psfMagKey = allMatches.schema.find("base_PsfFlux_mag").key
#    apMagKey = allMatches.schema.find("base_CircularApertureFlux_12_0_mag").key
    extendedKey = newSchema["base_ClassificationExtendedness_value"].asKey()

    # Filter further to a limited range in magnitude and extendedness to select bright stars (plotted below). The upturn at the bright end is due to the brighter-fatter effect on the sensors, but shouldn't matter for this test because it will affect all visits in approximately the same way (since they all have the same exposure time).
    safeMinMag = 16.5
#    safeMaxMag = 18.0
    safeMaxMag = 21.0
    safeMaxExtended = 1.0
    def safeFilter(cat):
        psfMag = np.mean(cat.get(psfMagKey))
#        modelMag = np.mean(cat.get(modelMagKey))
#        extend = np.abs(psfMag - modelMag)
        extended = np.max(cat.get(extendedKey))
        return psfMag >= safeMinMag and extended < safeMaxExtended

    safeMatches = goodMatches.where(safeFilter)

#    SUMMARY_STATISTICS = False
    SUMMARY_STATISTICS = True
    if SUMMARY_STATISTICS:
        # Pass field=psfMagKey so np.mean just gets that as its input
        goodPsfMag = goodMatches.aggregate(np.mean, field=psfMagKey)
        goodPsfMagRms = goodMatches.aggregate(np.std, field=psfMagKey)  # mmag
        # positionRms knows how to query a group so we give it the whole thing
        #   by going with the default `field=None`.
        dist = goodMatches.aggregate(positionRms)
    else:
        goodPsfMag = goodMatches.apply(lambda x: x, field=psfMagKey)
        dist = goodMatches.apply(positionDiff)

    return pipeBase.Struct(
        mag = goodPsfMag,
        mmagrms = 1000*goodPsfMagRms,  # mag -> mmag
        dist = dist,
        match = len(dist)
    )


def plotAstrometry(mag, mmagrms, dist, match, good_mag_limit=19.5):
    """Plot angular distance between matched sources from different exposures."""

    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['mathtext.default'] = 'regular'

    fig, ax = plt.subplots(ncols=2, nrows=3, figsize=(18, 22))
    ax[0][0].hist(dist, bins=80, histtype='stepfilled')
    ax[0][1].scatter(mag, dist, s=10, color='b')
    ax[0][0].set_xlim([0., 900.])
    ax[0][0].set_xlabel("Distance in mas", fontsize=20)
    ax[0][0].tick_params(labelsize=20)
    ax[0][0].set_title("Median : %.1f mas" % (np.median(dist)), fontsize=20, x=0.6, y=0.88)
    ax[0][1].set_xlabel("Magnitude", fontsize=20)
    ax[0][1].set_ylabel("Distance in mas", fontsize=20)
    ax[0][1].set_ylim([0., 900.])
    ax[0][1].tick_params(labelsize=20)
    ax[0][1].set_title("Number of matches : %d" % match, fontsize=20)

    ax[1][0].hist(dist, bins=150, histtype='stepfilled')
    ax[1][0].set_xlim([0., 400.])
    ax[1][1].scatter(mag, dist, s=10, color='b')
    ax[1][1].set_ylim([0., 400.])
    ax[1][0].set_xlabel("Distance in mas", fontsize=20)
    ax[1][1].set_xlabel("Magnitude", fontsize=20)
    ax[1][1].set_ylabel("Distance in mas", fontsize=20)
    ax[1][0].tick_params(labelsize=20)
    ax[1][1].tick_params(labelsize=20)

    idxs = np.where(np.asarray(mag) < good_mag_limit)

    ax[2][0].hist(np.asarray(dist)[idxs], bins=100, histtype='stepfilled')
    ax[2][0].set_xlabel("Distance in mas - mag < %.1f" % good_mag_limit, fontsize=20)
    ax[2][0].set_xlim([0, 200])
    ax[2][0].set_title("Median (mag < %.1f) : %.1f mas" %
                       (good_mag_limit, np.median(np.asarray(dist)[idxs])),
                       fontsize=20, x=0.6, y=0.88)
    ax[2][1].scatter(np.asarray(mag)[idxs], np.asarray(dist)[idxs], s=10, color='b')
    ax[2][1].set_xlabel("Magnitude", fontsize=20)
    ax[2][1].set_ylabel("Distance in mas - mag < %.1f" % good_mag_limit, fontsize=20)
    ax[2][1].set_ylim([0., 200.])
    ax[2][0].tick_params(labelsize=20)
    ax[2][1].tick_params(labelsize=20)

    plt.suptitle("Astrometry Check", fontsize=30)
    plotPath = "check_astrometry.png"
    plt.savefig(plotPath, format="png")


def checkAstrometry(mag, mmagrms, dist, match,
                    good_mag_limit=19.5,
                    medianRef=100, matchRef=500):
    """Print out the astrometric scatter for all stars, and for good stars.
    @param mag    Magnitude.  List or numpy.array.
    @param dist   Separation from reference.  List of numpy.array
    @param match  Number of stars matched.  Integer.

    @param medianRef  Median reference astrometric scatter in arcseconds.
    @param matchRef   Should match at least matchRef stars.

    Return the astrometric scatter (RMS, arcsec) for all good stars.

    Notes:
       The scatter and match defaults are appropriate to SDSS are stored here.
    For SDSS, stars with mag < 19.5 should be completely well measured.
    """

    assert len(mag) == len(dist)
    print("Median value of the astrometric scatter - all magnitudes:",
          np.median(dist), "mas")
    print("Median value of the photometric scatter - all magnitudes:",
          np.median(mmagrms), "mmag")

    idxs = np.where(np.asarray(mag) < good_mag_limit)
    astromScatter = np.median(np.asarray(dist)[idxs])
    photomScatter = np.median(np.asarray(mmagrms)[idxs])
    print("Astrometric scatter (median) - mag < %.1f : %.1f %s" %
          (good_mag_limit, astromScatter, "mas"))
    print("Photometric scatter (median) - mag < %.1f : %.1f %s" %
          (good_mag_limit, photomScatter, "mmag"))

    if astromScatter > medianRef:
        print("Median astrometric scatter %.1f %s is larger than reference : %.1f %s " %
              (astromScatter, "mas", medianRef, "mas"))
    if match < matchRef:
        print("Number of matched sources %d is too small (shoud be > %d)" % (match, matchRef))

    return astromScatter


def run(repo, visitDataIds, good_mag_limit, medianRef, matchRef):
    """Main executable.
    """

    struct = loadAndMatchData(repo, visitDataIds)
    mag = struct.mag
    mmagrms = struct.mmagrms
    dist = struct.dist
    match = struct.match
    checkAstrometry(mag, mmagrms, dist, match,
                    good_mag_limit=good_mag_limit,
                    medianRef=medianRef, matchRef=matchRef)
    plotAstrometry(mag, mmagrms, dist, match, good_mag_limit=good_mag_limit)


def defaultData(repo):
    """Example of loading dataIds for use by checkAstrometry.

    This example is based on the CFHT data in validation_data_cfht
    and is provided here for reference.
    For general usage, write your own equivalent to defaultData 
    and then pass to the `checkAstrometry.run` method.
    See the same `__main___` below for an example.
    """
    # List of visits to be considered
    visits = [176837, 176846, 176850]

    # List of CCD to be considered (source catalogs will be concateneted)
    ccd = [10]  # , 12, 14, 18]
    filter = 'z'

    # Reference values for the median astrometric scatter and the number of matches
    good_mag_limit = 21
    medianRef = 25
    matchRef = 5000

    visitDataIds = [{'visit': v, 'filter': filter, 'ccdnum': c} for v in visits
                    for c in ccd]

    return visitDataIds, good_mag_limit, medianRef, matchRef


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("""Usage: check_astrometry repo
where repo is the path to a repository containing the output of processCcd
""")
        sys.exit(1)

    repo = sys.argv[1]
    if not os.path.isdir(repo):
        print("Could not find repo %r" % (repo,))
        sys.exit(1)

    args = defaultData(repo)
    run(repo, *args)
