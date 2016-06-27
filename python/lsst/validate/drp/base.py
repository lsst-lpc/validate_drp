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

import abc
import numpy as np
import astropy.units


class ValidateError(Exception):
    """Base classes for exceptions in validate_drp."""
    pass


class ValidateErrorNoStars(ValidateError):
    """To be raised by tests that find no stars satisfying a set of criteria.

    Some example cases that might return such an error:
    1. There are no stars between 19-21 arcmin apart.
    2. There are no stars in the given magnitude range.
    """
    pass


class ValidateErrorSpecification(ValidateError):
    """Indicates an error with accessing or using requirement specifications."""
    pass


class ValidateErrorUnknownSpecificationLevel(ValidateErrorSpecification):
    """Indicates the requested level of requirements is unknown."""
    pass


class JsonSerializationMixin(object):
    """Mixin that provides serialization support"""

    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def json(self):
        """a `dict` that can be serialized as semantic SQuaSH json."""
        pass

    @staticmethod
    def jsonify_dict(d):
        """Recursively render JSON on all values in a dict."""
        json_dict = {}
        for k, v in d.iteritems():
            json_dict[k] = JsonSerializationMixin._jsonify_value(v)
        return json_dict

    @staticmethod
    def jsonify_list(lst):
        json_array = []
        for v in lst:
            json_array.append(JsonSerializationMixin._jsonify_value(v))
        return json_array

    @staticmethod
    def _jsonify_value(v):
        if isinstance(v, JsonSerializationMixin):
            return v.json
        elif isinstance(v, dict):
            return JsonSerializationMixin.jsonify_dict(v)
        elif isinstance(v, (list, tuple)):
            return JsonSerializationMixin.jsonify_list(v)
        else:
            return v


class Datum(JsonSerializationMixin):
    """A value annotated with units, a plot label and description.

    Parameters
    ----------
    value : `str`, `int`, `float` or 1-D iterable.
        Value of the datum.
    units : `str`
        Astropy-compatible units string. See
        http://docs.astropy.org/en/stable/units/
    label : `str`, optional
        Label suitable for plot axes (without units).
    description : `str`, optional
        Extended description.
    """
    def __init__(self, value, units, label=None, description=None):
        self._doc = {}
        self.value = value
        self.units = units
        self.label = label
        self.description = description

    @property
    def json(self):
        """Datum as a `dict` compatible with overall Job JSON schema."""
        # Copy the dict so that the serializer is immutable
        v = self.value
        if isinstance(v, np.ndarray):
            v = v.tolist()
        d = {
            'value': v,
            'units': self.units,
            'label': self.label,
            'description': self.description
        }
        return d

    @property
    def value(self):
        """Value of the datum (`str`, `int`, `float` or 1-D iterable.)."""
        return self._doc['value']

    @value.setter
    def value(self, v):
        self._doc['value'] = v

    @property
    def units(self):
        """Astropy-compatible unit string."""
        return self._doc['units']

    @units.setter
    def units(self, value):
        # verify that Astropy can parse the unit string
        if value is not None and value != 'millimag':
            print(value)
            astropy.units.Unit(value, parse_strict='raise')
        self._doc['units'] = value

    @property
    def label(self):
        """Label for plotting (without units)."""
        return self._doc['label']

    @label.setter
    def label(self, value):
        assert isinstance(value, basestring) or None
        self._doc['label'] = value

    @property
    def description(self):
        """Extended description of Datum."""
        return self._doc['description']

    @description.setter
    def description(self, value):
        assert isinstance(value, basestring) or None
        self._doc['description'] = value


class Metric(JsonSerializationMixin):
    """Docstring for Metric. """
    def __init__(self, name, description, operatorStr, specs=None,
                 referenceDoc=None, referenceUrl=None, referencePage=None):
        self.name = name
        self.description = description
        self.referenceDoc = referenceDoc
        self.referenceUrl = referenceUrl
        self.referencePage = referencePage

        self._operatorStr = operatorStr
        self.operator = Metric.convertOperatorStr(operatorStr)

        if specs is None:
            self.specs = []
        else:
            self.specs = specs

    def fromYamlDoc(self, yamlDoc, metricName, resolveDependencies=True):
        metricDoc = dict(yamlDoc[metricName])

        m = Metric(
            metricName,
            description=metricDoc.pop('description'),
            operatorStr=metricDoc.pop('operatorStr'),
            referenceDoc=metricDoc['reference'].pop('doc'),
            referenceUrl=metricDoc['reference'].pop('url'),
            referencePage=metricDoc['reference'].pop('page'))

        for specDoc in metricDoc['specs']:
            if 'dependencies' in specDoc and resolveDependencies:
                deps = {}
                for depItem in specDoc['dependencies']:
                    if isinstance(depItem, basestring):
                        # This is a metric
                        name = depItem
                        d = Metric.fromYamlDoc(yamlDoc, depItem,
                                               resolveDependencies=False)
                    elif isinstance(depItem, dict):
                        # Likely a Datum
                        # in yaml, wrapper object is dict with single key-val
                        name = depItem.keys()[0]
                        depItem = dict(depItem[name])
                        v = depItem.pop('value')
                        units = depItem.pop('units')
                        d = Datum(v, units, **depItem)
                    else:
                        raise RuntimeError(
                            'Cannot process dependency %r' % depItem)
                    deps[d.name] = d
            spec = Specification(name=depItem.pop('level'),
                                 value=depItem.pop('value'),
                                 units=depItem.pop('units'),
                                 filters=depItem.pop('filters'),
                                 dependencies=deps)
            m.specs.append(spec)

    def getSpec(self, name, bandpass=None):
        """Get a specification by name, and other qualitifications."""
        candidates = []

        # First collect candidate specifications by name
        for spec in self.specs:
            if spec.name == name:
                candidates.append(spec)
        if len(candidates) == 1:
            return candidates[0]

        # Filter down by bandpass
        if bandpass is not None:
            for spec in candidates:
                if hasattr(spec, 'bandpasses') and bandpass not in spec:
                    candidates.remove(spec)
        if len(candidates) == 1:
            return candidates[0]

        raise RuntimeError(
            'No spec found for name={0} bandpass={1}'.format(name, bandpass))

    def checkSpec(self, value, specName, bandpass=None):
        """Compare a measurement `value` against a named specification level
        (:class:`SpecLevel`).

        Returns
        -------
        passed : bool
            `True` if a value meets the specification, `False` otherwise.
        """
        spec = self.getSpec(specName, bandpass=None)

        # NOTE: assumes units are the same
        return self.operator(value, spec.value)

    @property
    def json(self):
        """Render metric as a JSON object (`dict`)."""
        return JsonSerializationMixin.jsonify_dict({
            'name': self.name,
            'reference': self.referenceDoc,
            'description': self.description,
            'specifications': self.specs})


class Specification(JsonSerializationMixin):
    """A specification level or threshold associated with a Metric."""
    def __init__(self, name, value, units, bandpasses=None, dependencies=None):
        self.name = name
        self.value = value
        self.units = units
        self.bandpasses = bandpasses
        self.dependencies = dependencies

    def __getattr__(self, key):
        """Access dependencies with keys as attributes."""
        return self.dependencies[key]

    @property
    def json(self):
        return JsonSerializationMixin.jsonify_dict({
            'name': self.level,
            'value': Datum(self.value, self.units),
            'filters': self.bandpasses,
            'dependencies': self.dependencies})


class MetricBase(JsonSerializationMixin):
    """Baseclass for Metric definition classes."""
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def name(self):
        """Name of the metric."""
        pass

    @abc.abstractproperty
    def reference(self):
        """Identifier for document that defines the metric."""
        pass

    @abc.abstractproperty
    def description(self):
        """Short description of the metric."""
        pass

    @abc.abstractproperty
    def operator(self):
        """Function from the `operator` module (e.g.
        :py:class:`operator.gt`:."""
        pass

    @abc.abstractproperty
    def specs(self):
        """`dict` of :class:`SpecLevel` objects, keyed by role ."""
        pass

    def checkSpec(self, value, specName):
        """Compare a measurement `value` against a named specification level
        (:class:`SpecLevel`).

        Returns
        -------
        passed : bool
            `True` if a value meets the specification, `False` otherwise.
        """
        specLevel = self.specs[specName]
        return self.operator(value, specLevel)

    @property
    def json(self):
        """Render metric as a JSON object (`dict`)."""
        # FIXME add spec_level for dependent metrics
        # For example, PF1 is dependent on the specification level
        # to determine PA2.
        return {
            'name': self.name,
            'reference': self.reference,
            'description': self.description}


class SpecLevel(object):
    """Specification Level definition classes.

    Attributes
    ----------
    value : float or int
        A scalar value that defines a specification threshold level for a
        Metric.
    units : str
        Astropy-compatible unit string annotating a value.

    Parameters
    ----------
    value : float or int
        A scalar value that defines a specification threshold level for a
        Metric.
    units : str
        Astropy-compatible unit string annotating a value.
    """
    def __init__(self, value, units):
        self.value = value
        self.units = units


class MeasurementBase(JsonSerializationMixin):
    """Baseclass for Measurement classes.

    Attributes
    ----------
    metric
    units
    label
    json
    schema
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self):
        self._linkedBlobs = []
        self._params = {}

    def linkBlob(self, blob):
        """Add a blob so that it will be linked to this measurement in
        the serialized JSON. This does not cause the blob itself to be
        persisted in a Job document, however.
        """
        self._linkedBlobs.append(blob)

    def registerParameter(self, paramKey, paramValue):
        """Register a measurement input parameter so that it can be persisted
        with the measurement.
        """
        self._params[paramKey] = paramValue

    @abc.abstractproperty
    def metric(self):
        """An instance derived from
        :class:`~lsst.validate.drp.base.MetricBase`.
        """
        pass

    @abc.abstractproperty
    def value(self):
        """Metric measurement value."""
        pass

    @abc.abstractproperty
    def units(self):
        """Astropy-compatible units string. (`str`)."""
        pass

    @abc.abstractproperty
    def label(self):
        """Lable (`str`) suitable for plot axes; without units."""
        pass

    @abc.abstractproperty
    def schema(self):
        """A `str` identifying the schema of this measurement and parameters.
        """
        pass

    @property
    def json(self):
        """a `dict` that can be serialized as semantic SQuaSH json."""
        object_doc = {'metric': self.metric,
                      'value': self.value,
                      'parameters': self._params,
                      'blobs': self._linkedBlobs,
                      'schema': self.schema}
        json_doc = JsonSerializationMixin.jsonify_dict(object_doc)
        return json_doc


class BlobBase(JsonSerializationMixin):
    """Baseclass for Blob classes. Blobs are flexible containers of data
    that are serialized to JSON.
    """

    def __init__(self):
        self._doc = {}

    @abc.abstractproperty
    def schema(self):
        """A `str` identify the schema of this blob."""
        pass

    @property
    def json(self):
        json_doc = JsonSerializationMixin.jsonify_dict(self._doc)
        return json_doc
