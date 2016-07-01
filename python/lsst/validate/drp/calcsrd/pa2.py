# LSST Data Management System
# Copyright 2016 AURA/LSST.
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

import numpy as np

from ..base import MeasurementBase, Metric
from ..util import getRandomDiffRmsInMas


class PA2Measurement(MeasurementBase):
    """Measurement of PA2: millimag from median RMS (see PA1) of which
    PF1 of the samples can be found.

    Parameters
    ----------
    matchedDataset : lsst.validate.drp.matchreduce.MatchedMultiVisitDataset
    bandpass : str
        Bandpass (filter name) used in this measurement (e.g., `'r'`).
    specName : str
        Name of a specification level to measure against (e.g., design,
        minimum, stretch).
    verbose : bool, optional
        Output additional information on the analysis steps.
    linkedBlobs : list, optional
        A `list` of additional blobs (subclasses of BlobSerializerBase) that
        can provide additional context to the measurement, though aren't
        direct dependencies of the computation (e.g., `matchedDataset).

    Notes
    -----
    The LSST Science Requirements Document (LPM-17) is commonly referred
    to as the SRD.  The SRD puts a limit that no more than PF1 % of difference
    will vary by more than PA2 millimag.  The design, minimum, and stretch
    goals are PF1 = (10, 20, 5) % at PA2 = (15, 15, 10) millimag following
    LPM-17 as of 2011-07-06, available at http://ls.st/LPM-17.
    """

    metric = None
    value = None
    units = 'millimag'
    label = 'PA2'
    schema = 'pa2-1.0.0'

    # FIXME somehow this isn't depdendent on a median RMS (e.g., PA1 measure)
    def __init__(self, matchedDataset, bandpass, specName, verbose=False,
                 linkedBlobs=None, metricYamlDoc=None, metricYamlPath=None):
        MeasurementBase.__init__(self)
        self.bandpass = bandpass
        self.specName = specName  # spec-dependent measure because of PF1 dep
        self.metric = Metric.fromYaml(self.label,
                                      yamlDoc=metricYamlDoc,
                                      yamlPath=metricYamlPath)

        # TODO register input parameters for serialization
        # note that matchedDataset is treated as a blob, separately

        self.matchedDataset = matchedDataset

        # Add external blob so that links will be persisted with
        # the measurement
        if linkedBlobs is not None:
            for blob in linkedBlobs:
                self.linkBlob(blob)
        self.linkBlob(self.matchedDataset)

        matches = matchedDataset.safeMatches
        magKey = matchedDataset.magKey
        magDiffs = matches.aggregate(getRandomDiffRmsInMas, field=magKey)
        # FIXME Get magdiffs from a blob supplied by PA1

        pf1Val = self.metric.getSpec(specName, bandpass=self.bandpass).\
            PF1.getSpec(specName, bandpass=self.bandpass).value
        pf1Percentile = 100. - pf1Val
        self.value = np.percentile(np.abs(magDiffs), pf1Percentile)
