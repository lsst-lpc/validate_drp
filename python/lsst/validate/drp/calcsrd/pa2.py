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

import lsst.pipe.base as pipeBase

from ..base import MeasurementBase, Metric
from .utils import (
    getRandomDiffRmsInMas, computeWidths, radiansToMilliarcsec,
    arcminToRadians, sphDist, matchVisitComputeDistance)


class PA2Measurement(MeasurementBase):
    """Docstring for PA2Measurement. """

    metric = None
    value = None
    units = 'millimag'
    label = 'PA2'
    schema = 'pa1-1.0.0'

    def __init__(self, matchedDataset, bandpass,
                 numRandomShuffles=50, verbose=False,
                 linkedBlobs=None, metricYamlDoc=None, metricYamlPath=None):
        MeasurementBase.__init__(self)
        self.bandpass = bandpass
        self.metric = Metric.fromYaml(self.label,
                                      yamlDoc=metricYamlDoc,
                                      yamlPath=metricYamlPath)

        # register input parameters for serialization
        # note that matchedDataset is treated as a blob, separately
        self.registerParameter('num_random_shuffles', numRandomShuffles)

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
        PF1_percentiles = 100. - np.asarray([srdSpec.PF1[l]
                                             for l in srdSpec.levels])
        PA2_measured = dict(zip(srdSpec.levels,
                                np.percentile(np.abs(magDiffs),
                                              PF1_percentiles)))

        PF1_measured = {l: 100*np.mean(np.asarray(magDiffs) > srdSpec.PA2[l])
                        for l in srdSpec.levels}

        return pipeBase.Struct(name='PA2', pa2Units='mmag', pf1Units='%',
                            PA2=PA2_measured['design'], PF1=PF1_measured['design'],
                            PA2_measured=PA2_measured,
                            PF1_measured=PF1_measured,
                            PF1_spec=srdSpec.PF1, PA2_spec=PA2_spec)
