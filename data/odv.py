# -*- coding: utf-8 -*-
"""ODV data manager module.

This module manages Ocean Data View's spreadsheed formated text files.

Disclaimer
----------
This software may be used, copied, or redistributed as long as it is
not sold and this copyright notice is reproduced on each copy made.
This routine is provided as is without any express or implied
warranties whatsoever.

Author
------
Sebastian Krieger (sebastian.krieger@usp.br)

Revision
--------
1 (2016-04-21 12:23 -0300)

"""
from __future__ import division, unicode_literals

__version__ = '$Revision: 1 $'
# $Source$

import codecs
import re

from collections import OrderedDict
from datetime import datetime
from lxml import etree
from matplotlib.dates import date2num, num2date, strpdate2num
from numpy import arange, array, genfromtxt, ma
from StringIO import StringIO

import atlantis.data


class Sequence(atlantis.data.Sequence):
    """
    Manages ODV-like data files.

    Parameters
    ----------
    url : string, optional
        Gives the full path for data file. If file exists, tries to
        read its content and rebuild dataset. Otherwise initializes a
        blank database.

    Examples
    --------
    >> dat = odv.Sequence('odv_file.tsv')

    """
    # Valid attribute list.
    _attr_list = ['__name__', 'url', 'time', 'fields', 'history']

    def __init__(self, url=None, **kwargs):
        # Resets fields dictionary.
        self.fields = OrderedDict()
        if url is None:
            # Creates default data fields:
            #   cruise, station, type, yyyy-mm-dd, hh:mm, dec_year, longitude,
            #   latitude, bot. depth, depth
            fields = ['Cruise', 'Station', 'Type', 'YYYY-MM-DD', 'hh:mm',
                'Decimal year', 'Julian day', 'Longitude', 'Latitude',
                'Bot. depth', 'Depth']
            units = [None, None, None, None, None, None, None, 'degrees_east',
                'degrees_north', 'm', 'm']
            fmt = ['{}', '{}', '{}', '{}', '{}', '{:.6f}', '{:.4f}', '{:.6f}', '{:.6f}', '{:.1f}', '{:.2f}']
            #
            for c, u, f in zip(fields, units, fmt):
                var = atlantis.data.Variable(
                    standard_name=c,
                    units=u,
                    string_format=f
                )
                self.__setitem__(c, var)
            #
            self.__setattr__('time', None)
        else:
            self.read(url)
        #
        self.__setattr__('url', url)
        #
        return None


    def __len__(self):
        try:
            return len(self.time)
        except:
            return 0


    def __getattr__(self, name):
        # Checks whether attribute name is valid:
        if name not in self._attr_list:
            raise ValueError('Invalid attribute `{}`.'.format(name))
        # Returns attribute value or None if value was not assigned.
        try:
            return self.__dict__[name]
        except:
            return None


    def __setattr__(self, name, value):
        # Assigns value to class' dict of attributes.
        if name == 'time':
            value = self._check_attribute_time(value)
        self.__dict__[name] = value


    def __getitem__(self, key):
        if isinstance(key, list):
            return OrderedDict([(item, self[item]) for item in key])
        elif isinstance(key, basestring):
            return self.fields[key]
        else:
            raise ValueError('Unable to get field `{}`.'.format(key))


    def __setitem__(self, key, value):
        # Makes sure that value is a Variable object
        if not isinstance(value, atlantis.data.Variable):
            raise ValueError('`{}` is not a Variable!'.format(key))
        if (~isinstance(value.data, ma.MaskedArray)) & \
                    (value.data is not None):
            value.data = ma.masked_invalid(value.data)
        if (~isinstance(value.standard_deviation, ma.MaskedArray)) & \
                    (value.standard_deviation is not None):
            value.standard_deviation = \
                ma.masked_invalid(value.standard_deviation)
        # Assigns value to class' list of items (fields)
        self.fields[key] = value


    def __delitem__(self, key):
        del self.fields[key]


    def __iter__(self):
        return iter(self.fields)


    def __str__(self):
        return self._dump().encode('utf-8')


    def _check_attribute_time(self, value, update=True):
        """
        Checks new time values and update helper variables
        if appropriate.

        """
        if value is None:
            return None
        # Checks instance type of first parameter, and assumes the rest is of
        # same type.
        if isinstance(value[0], datetime):
            T = value
            value = date2num(T)
        elif isinstance(value[0], float):
            T = num2date(value)
        # If `update` is True, then changes the values of helper fields.
        if update:
            D = ['{:04d}-{:02d}-{:02d}'.format(t.year, t.month, t.day)
                for t in T]
            H = ['{:02d}:{:02d}'.format(t.hour, t.minute) for t in T]
            DY, JD = self._convert_datenum(value)
            if 'YYYY-MM-DD' in self.fields:
                self['YYYY-MM-DD'].data = D
            if 'hh:mm' in self.fields:
                self['hh:mm'].data = H
            if 'Decimal year' in self.fields:
                self['Decimal year'].data = DY
            if 'Julian day' in self.fields:
                self['Julian day'].data = JD
        #
        return value


    def _dump_header(self, s, dtype='info', result='list'):
        """Generates ODV file header."""
        dump = []
        for key, value in s.items():
            if isinstance(value, basestring):
                dump.append(u'//<{0}>{1}</{0}>'.format(key, value))
            elif isinstance(value, atlantis.data.Variable):
                if value.string_format is None:
                    value.string_format = '{:.2f}'
                s = '//<{}>'.format('DataVariable')
                if value.units is not None:
                    s += 'label=\"{} [{}]\"'.format(key, value.units)
                elif value.canonical_units is not None:
                    s += 'label=\"{} [{}]\"'.format(key, value.canonical_units)
                else:
                    s += 'label=\"{}\"'.format(key)
                if value.standard_name:
                    s += ' standard_name=\"{}\"'.format(value.standard_name)
                if value.long_name:
                    s += ' long_name=\"{}\"'.format(value.long_name)
                if value.string_format:
                    s += ' string_format=\"{}\"'.format(value.string_format)
                if value.symbol:
                    s += ' symbol=\"{}\"'.format(value.symbol)
                #if value.description:
                #    s += 'description=\"{}\"'.format(value.description)
                s += '</{}>'.format('DataVariable')
                dump.append(s)
        #
        if result == 'list':
            return dump
        elif result == 'string':
            return '\n'.join(dump)
        else:
            raise ValueError('Invalid type `{}`'.format(result))


    def _dump(self, fields=None):
        """Returns human-readable data dump."""
        # Checks input parameters.
        if fields is None:
            fields = self.fields
        # Header information
        file_header = OrderedDict([
            ('Encoding', 'UTF-8'),
            ('Creator', 'Atlantis ODV data manager {}'.format(__version__)),
            ('CreateTime', '{} UTC'.format(datetime.now().isoformat())),
            ('MissingValueIndicators', None),
            ('Software', None),
            ('Source', None),
            ('SourceLastModified', None),
            ('Version', 'ODV Spreasheet V4.6')
        ])
        data_header = OrderedDict([
            ('DataField', 'Ocean'),
            ('DataType', None),
            ('History', self.history)
        ])
        # Initializes dump
        dump = []
        # Appends file header to dump.
        dump += self._dump_header(file_header)
        dump += self._dump_header(data_header)
        if fields is None:
            dump += self._dump_header(self.fields)
        else:
            dump += self._dump_header(
                OrderedDict([(field, self[field]) for field in fields])
            )
        # Creates data header and reads string format for each field.
        fields_labels = []
        string_formats = {}
        for f in fields:
            units = self[f].units
            if (units is not None) & (units != ''):
                fields_labels.append('{} [{}]'.format(f, units))
            else:
                fields_labels.append(f)
            #
            s = self[f].string_format
            if (s is None) | (s == ''):
                string_formats[f] = '{:.2f}'
            else:
                string_formats[f] = s
        #
        dump.append('\t'.join(fields_labels))
        #
        # Walks through each data entry
        n = len(self)
        for i in xrange(n):
            line = []
            for f in fields:
                s = string_formats[f]
                try:
                    line.append(s.format(self[f].data[i]))
                except:
                    line.append('')
            # Appends line to data dump.
            dump.append('\t'.join(line))
        #
        return '\n'.join(dump)


    def _convert_datenum(self, t):
        """Returns decimal year and Julian day."""
        T = num2date(t)
        Y = array([_t.year for _t in T], dtype=float)
        T0 = array(date2num([datetime(year=_t.year, month=1, day=1, hour=0,
            minute=0) for _t in T]), dtype=float)
        T1 = array(date2num([datetime(year=_t.year+1, month=1, day=1, hour=0,
            minute=0) for _t in T]), dtype=float)
        julian = (t - T0) + 1
        decimal_year = Y + (t - T0) / (T1 - T0)
        #
        return decimal_year, julian


    def _read_text_file(self, url):
        """Returns all content of a given text file."""
        try:
            with codecs.open(url, 'rb', encoding='utf-8') as f:
                return f.read()
        except:
            with open(url, 'rb') as f:
                return f.read()


    def _get_header(self, s, comments='//', delimiter='\t'):
        """Analyzes the header of ODV file."""
        # Makes sure that file content is given as a list.
        if isinstance(s, basestring):
            s = s.splitlines()
        # Lenght of the comment character(s)
        comments_len = len(comments)
        # Initializes header dictionary
        header = []
        # Walks through each line to extract header data.
        for line_n, line in enumerate(s):
            # If first character(s) is a comment, get the header data it
            # contains.
            if line[:comments_len] == comments:
                header.append(line[comments_len:])
            # Otherwise exits loop
            else:
                break

        # Initializes fields and values lists.
        fields = OrderedDict()

        # The first data line contains the data columns and here we extract
        # field names and units.
        columns = s[line_n].split(delimiter)
        units = [''] * len(columns)
        for i, column in enumerate(columns):
            columns[i], units[i] = self._field_unit_from_column(column)
            fields[columns[i]] = atlantis.data.Variable(units=units[i])

        # Walks through each header line and retrieves variable parameters. The
        # header parameters are treated as XML data.
        xml = '<root>\n' + '\n'.join(header) + '</root>\n'
        root = etree.fromstring(xml)
        # Regular expression pattern for key="value".
        pattern = '([^=]+)=\"([^\"]+)\" ?'
        # Alias dictionary from ODV DataVariable parameters and data.Variable
        # attributes: 'long_name', 'standard_name', 'units',
        # 'canonical_units', 'description', 'grib', 'amip', 'id',
        # 'missing_value',  'scale_factor', 'add_offset', 'grids', 'data',
        # 'flag_values',  'flag_meanings', 'string_format', 'symbol')
        attrib_keys = dict(standard_name='standard_name',
            long_name='long_name', units='units',
            string_format='string_format', description='description',
            symbol='symbol')
        # Walks through each childeren in root.
        for child in root.iterchildren():
            if child.tag == 'DataVariable':
                attributes = dict()
                for key, value in re.findall(pattern, child.text):
                    if key == 'label':
                        label, units = self._field_unit_from_column(column)
                        attributes['units'] = units
                    else:
                        attributes[attrib_keys[key]] = value
                # Appends attributes to variable
                if label in fields.keys():
                    for key, value in attributes.items():
                        setattr(fields[label], key, value)
                else:
                    fields[label] = atlantis.data.Variable(**attributes)

        # Treats each header entry as XML.
        return header, fields, line_n


    def _field_unit_from_column(self, column):
        """Extracts field name and unit from column header."""
        pattern = '(\S+) \[(.+)\]'
        try:
            return re.match(pattern, column).groups()
        except:
            return column, ''


    def read(self, url=None, fields=None, delimiter='\t'):
        """
        Reads data from ODV formatted file.

        Parameters
        ----------
        url : string, optional
            Full path and file name to save the data. If omitted,
            assumes path indicated at sequence initialization.
        fields : sequence, optional
            Sets the fields to be saved. Default is to save all fields
            in dataset.

        Returns
        -------
        dat : array like
            Structured array of file contents.

        """
        if fields is not None:
            raise ValueError('Not implemented yet.')
        # Reads the content of file.
        if url is None:
            url = self.url
        f = self._read_text_file(url)
        # Splits all content lines and analyzes the header to extract
        # additional information
        header, fields, skip_header = self._get_header(f,
            comments='//', delimiter=delimiter)
        keys = fields.keys()
        # Sets data converters according to field names.
        converters = dict()
        for i, key in enumerate(keys):
            if key == 'YYYY-MM-DD':
                converters[i] = strpdate2num('%Y-%m-%d')
            elif key == 'hh:mm':
                converters[i] = strpdate2num('%H:%M')
        # Converts data content in structured array using numpy.genfromtxt.
        dat_keys = [b'{:d}'.format(a) for a in range(len(keys))]
        dat = genfromtxt(url, delimiter=delimiter, skip_header=skip_header+1,
            dtype=None, names=dat_keys, converters=converters)
        # Sets data in field container.
        for dat_key, key in zip(dat_keys, keys):
            fields[key].data = dat[dat_key]
        # Updates class containers
        self.fields = fields
        # Update date and time.
        T0 = 693596.0 #strpdate2num('%H:%M')('00:00')
        self.time = fields['YYYY-MM-DD'].data + fields['hh:mm'].data - T0
        # Returns data structured array.
        return dat


    def write(self, url=None, fields=None):
        """
        Saves data to ODV formatted file.

        Parameters
        ----------
        url : string, optional
            Full path and file name to save the data. If omitted,
            assumes path indicated at sequence initialization.
        fields : sequence, optional
            Sets the fields to be saved. Default is to save all fields
            in dataset.

        Returns
        -------
        Nothing

        """
        if url is None:
            url = self.url
        else:
            self.url = url
        #
        with file(url, 'w') as f:
            f.write(self._dump(fields=fields).encode('utf-8'))

        return None


    def append_fields(self, names, data, units=None, string_format=None,
        duplicate='error'):
        """
        Appends fields to dataset.

        Parameters
        ----------
        names : string, Sequence
            String or sequence of strings correspondig to the names of
            the new fields.
        data : array or sequence of arrays
            Array or sequence of arrays storing the fields to add to
            the base. Data can also be stored as a data.Sequence
            object.
        units : string or sequence, optional
            String or sequence of strings assigning data units to the
            new fields.
        string_format : string or sequence, optional
            String or sequence of strings assigning the format for data
            output.
        duplicate : string, optional
            Sets behaviour when appending duplicated fields. Valid
            options are `error` (default), `skip`, `rename`,
            `overwrite`.

        Returns
        -------
        Nothing

        """
        # Checks if input parameters are. Assumes that all parameters have
        # same structure.
        if (not isinstance(names, list)) & (not isinstance(names, tuple)):
            names = [names]
        if (not isinstance(data, list)) & (not isinstance(data, tuple)):
            data = [data]
        if (not isinstance(units, list)) & (not isinstance(units, tuple)):
            units = [units]
        if ((not isinstance(string_format, list)) &
                (not isinstance(string_format, tuple))):
            string_format = [string_format]
        #
        nvars = len(names)
        if len(data) != nvars:
            raise ValueError('Data size does not match name size.')
        if (len(units) == 1) & (nvars > 1):
            units = units * nvars
        elif len(units) != nvars:
            raise ValueError('Units size does not match name size.')
        if (len(string_format) == 1) & (nvars > 1):
            string_format = string_format * nvars
        elif len(string_format) != nvars:
            raise ValueError('String format size does not match name size.')

        # Walks through each input parameter and appends data.
        for n, d, u, sf in zip(names, data, units, string_format):
            if n in self.fields:
                if duplicate == 'error':
                    raise ValueError('Field `{}` already exists in ' 'database'.format(n))
                elif duplicate == 'skip':
                    continue
                elif duplicate == 'rename':
                    raise ValueError('TODO.')
                elif duplicate == 'overwrite':
                    print ('*** Warning *** : duplicated entry `{}` ' 'overwritten.'.format(n))
                else:
                    raise ValueError('Invalid option `{}` for duplicated ' 'fields.'.format(duplicate))
            if isinstance(d, atlantis.data.Variable):
                # Makes sure that data is one-dimensional
                d.data = d.data.flatten()
                var = d
                d = var.data
            else:
                var = atlantis.data.Variable(
                    standard_name=n,
                    units=u,
                    string_format=sf,
                    data=d.flatten()
                )
            # Checks if data has the same size as time array
            if d.shape != self.time.shape:
                raise ValueError('New data has not same size as currently '
                    'stored data.')
            # Finally, appends the data.
            self.__setitem__(n, var)
        #
        return None
