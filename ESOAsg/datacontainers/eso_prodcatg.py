r"""Class that defines the properties of an ESO PRODCATG

The current implementation follows the ESO Science Data Products Standard Version 6.
More information in the `Document ESO-044286 <https://www.eso.org/sci/observing/phase3/p3sdpstd.pdf>`_.
"""

import collections
import pkg_resources
from astropy.io import ascii

from ESOAsg import msgs
from ESOAsg.ancillary import cleaning_lists

# Ordered dictionary that defines all properties of different PRODCATG
PRODCATG = collections.OrderedDict()
PRODCATG['SCIENCE.IMAGE'] = {'data_format': 'Single image stored in the primary HDU.'}
PRODCATG['SCIENCE.MEFIMAGE'] = {'data_format': 'Multiple images stored in multi-extension ' +
                                               'FITS format (MEF).'}
PRODCATG['SCIENCE.IMAGE.FLUXMAP'] = {'data_format': 'Flux map, single image stored in the primary HDU'}
PRODCATG['SCIENCE.SPECTRUM'] = {'data_format': 'Single target one-dimensional spectrum ' +
                                               'in binary table format (single row array).'}
PRODCATG['SCIENCE.SRCTBL'] = {'data_format': 'Source list stored as FITS binary table resulting ' +
                                             'from the detection of sources on an image (both single ' +
                                             'image or MEF image). Merged  multi-band source lists ' +
                                             'do not qualify as source lists but should be submitted ' +
                                             'as catalogue data instead.'}
PRODCATG['SCIENCE.CUBE.IFS'] = {'data_format': '3D cube for IFS with two spatial axes and one ' +
                                               'spectral axis.'}
PRODCATG['SCIENCE.VISIBILITY'] = {'data_format': 'Visibility and closure phase for optical interferometry.'}
PRODCATG['SCIENCE.CATALOG'] = {'data_format': 'Scientific catalogue in single FITS binary table.'}
PRODCATG['SCIENCE.MCATALOG'] = {'data_format': 'Metadata definition file of a catalogue submitted ' +
                                               'in a tile-by-tile fashion, i.e. partitioned in' +
                                               'multiple FITS binary tables.'}
PRODCATG['SCIENCE.CATALOGTILE'] = {'data_format': 'Data file for a catalogue submitted in a ' +
                                                  'tile-by-tile fashion, i.e. partitioned in ' +
                                                  'multiple FITS binary tables.'}

# Ordered Dictionary that defines the conditions
HEADER_KEYWORDS_TABLE_LEGEND = collections.OrderedDict()
HEADER_KEYWORDS_TABLE_LEGEND['M'] = 'Keyword mandatory in the primary header.'
HEADER_KEYWORDS_TABLE_LEGEND['Mext'] = 'Keyword mandatory in the extension header.'
HEADER_KEYWORDS_TABLE_LEGEND['Meso'] = 'Keyword mandatory for ESO files.'
HEADER_KEYWORDS_TABLE_LEGEND['Mapp'] = 'Keyword mandatory in the primary header, ' + \
                                       'where applicable, ' + \
                                       'refer to the notes for more information.'
HEADER_KEYWORDS_TABLE_LEGEND['Mextapp'] = 'Keyword mandatory in the extension header, ' + \
                                          'where applicable, ' + \
                                          'refer to the notes for more information.'
HEADER_KEYWORDS_TABLE_LEGEND['R'] = 'Keyword recommended; it means applicable ' + \
                                    'unless you have a good reason to not use it.'
HEADER_KEYWORDS_TABLE_LEGEND['Rext'] = 'Keyword recommended in the extension header; ' + \
                                       'it means applicable ' + \
                                       'unless you have a good reason to not use it.'
HEADER_KEYWORDS_TABLE_LEGEND['NotAlw'] = 'Keyword not allowed: the Phase 3 system ' + \
                                         'will reject the file.'
HEADER_KEYWORDS_TABLE_LEGEND['O'] = 'It makes sense to provide the value, but it ' + \
                                    'is optional.'
HEADER_KEYWORDS_TABLE_LEGEND['Oext'] = 'It makes sense to provide the value in the extension ' + \
                                       'header, but it is optional.'
HEADER_KEYWORDS_TABLE_LEGEND['N/A'] = 'Keyword not applicable, i.e. it is not ' + \
                                      'relevant in the context of the PRODCATG (or ' + \
                                      'it is computed by the P3 system) but it is ' + \
                                      'tolerated, that is it is not checked by our system.'
HEADER_KEYWORDS_TABLE_LEGEND['Reserved'] = 'Reserved keyword. If present, modified ' + \
                                           'by the Phase 3 system during the archiving process.'
HEADER_KEYWORDS_TABLE_LEGEND['Contact'] = 'Please contact usd-help@eso.org subject: `phase 3` ' + \
                                          'if you intend to submit polarisation data.'
HEADER_KEYWORDS_TABLE_LEGEND['None'] = 'Not in use for the selected PRODCATG'


def _read_header_keywords_notes_table() -> dict:
    r"""Read the notes of the table describing the `PRODCATG` keywords in the ESO Science Data Products Standard

    Returns:
        dict: Dict containing the notes

    """
    header_keywords_notes = pkg_resources.resource_filename('ESOAsg',
                                                            'datacontainers/data/p3sdpstd_header_keywords_notes.csv')
    header_keywords_table_notes = ascii.read(header_keywords_notes, format='csv', delimiter=';')
    notes_dict = {}
    for idx in range(0, len(header_keywords_table_notes)):
        notes_dict[str(header_keywords_table_notes['Number'][idx])] = header_keywords_table_notes['Note'][idx]
    return notes_dict


def _read_header_keywords_table(prodcatg_type) -> dict:
    r"""Read the table describing the `PRODCATG` header keywords in the ESO Science Data Products Standard

    For a given `PRODCATG` it reads the table containing the information regarding the header keywords and set it in
    a dictionary. This contains a dictionary for each `keyword` with keys:
        * `iterable`: if it is True it means that the keyword is actually a list of entries. For instance: `CRVALi` to
            indicate `CRVAL1`, `CRVAL2`, etc.
        * `condition`: describe if the header keyword is mandatory (`M`), mandatory for ESO files (`Meso`), etc. The
            full list is in the dictionary: `HEADER_KEYWORDS_TABLE_LEGEND`
        * `header`: indicates the header where the keyword apply
        * `note`: notes relative to a specific keyword

    Returns:
        dict: dictionary containing `iterable`, `condition`, `header`, and `note` for each `keyword`

    """
    header_keywords = pkg_resources.resource_filename('ESOAsg',
                                                      'datacontainers/data/p3sdpstd_header_keywords.csv')
    header_keywords_table = ascii.read(header_keywords, format='csv', delimiter=';', comment='#')
    header_keywords_notes = _read_header_keywords_notes_table()

    keyword_list = []
    iterable_list = []
    condition_list = []
    header_list = []
    note_list = []

    for col_name in header_keywords_table.colnames:
        if col_name.startswith(prodcatg_type):
            for idx_keyword in range(0, len(header_keywords_table['KEYWORD', col_name])):
                # setting keyword list
                _keyword_split = header_keywords_table['KEYWORD'][idx_keyword].split(' ')
                keyword_length = len(_keyword_split)
                if keyword_length == 1:
                    keyword_list.append(_keyword_split[0])
                elif keyword_length == 2:
                    keyword_list.append(_keyword_split[0])
                    keyword_list.append(_keyword_split[1])
                else:
                    msgs.error('Unexpected entry in the keyword list while reading: \n {}'.format(header_keywords))
                # setting the condition list
                _condition_split = header_keywords_table[col_name][idx_keyword].split(' ')
                for idx_condition in range(0, keyword_length):
                    _keyword_notes_list = []
                    _keyword_conditions_list = []
                    _keyword_header_list = []
                    for _condition_single in _condition_split:
                        if _condition_single in list(header_keywords_notes.keys()):
                            _keyword_notes_list.append(_condition_single)
                        elif _condition_single in list(HEADER_KEYWORDS_TABLE_LEGEND.keys()):
                            _keyword_conditions_list.append(_condition_single)
                            if 'ext' in _condition_single:
                                _keyword_header_list.append('EXTENSION')
                            else:
                                _keyword_header_list.append('PRIMARY')
                        else:
                            msgs.error('{} is an unexpected entry in the condition'.format(_condition_single) +
                                       'list in: \n {}'.format(header_keywords))
                    note_list.append(_keyword_notes_list)
                    condition_list.append(_keyword_conditions_list)
                    header_list.append(_keyword_header_list)
            # find iterable headers
            for _key in keyword_list:
                if _key.endswith('i') | _key.endswith('ia') | _key.endswith('i_j'):
                    iterable_list.append(True)
                else:
                    iterable_list.append(False)
    keywords_dict = {}
    for keyword, iterable, condition, header, note in zip(keyword_list, iterable_list, condition_list, header_list,
                                                          note_list):
        keywords_dict[keyword] = {'iterable': iterable,
                                  'condition': condition,
                                  'header': header,
                                  'note': note}
    return keywords_dict


def _get_header_table_legend(short_name) -> str:
    r"""Read the table

    Args:
        short_name (str):

    Returns:
        str: string containing the explicit value of `short_name`

    """
    if short_name in list(HEADER_KEYWORDS_TABLE_LEGEND.keys()):
        return HEADER_KEYWORDS_TABLE_LEGEND[short_name]
    else:
        raise ValueError('{} is not present in the legend table.'.format(short_name)
                         + 'Possible values are: \n {}'.format(HEADER_KEYWORDS_TABLE_LEGEND.keys()))


class ProdCatg:
    r"""Class that defines the data product categories of the ESO Archive
    
    Args:
        prodcatg_type (str): Name of the data product category

    Raises:
        TypeError: if `prodcatg` is not a `str`
        ValueError: if `prodcatg` is not a valid category accordingly to the ESO standard

    """

    def __init__(self, prodcatg_type=None):
        self.prodcatg = prodcatg_type

    def __str__(self):
        return 'ESO PRODCATG: {}'.format(self.prodcatg)

    @property
    def prodcatg(self):
        r"""Name of the `prodcatg`.

        This set also the `descriptor` dictionary containing the relevant properties of the `prodcatg`.

        """
        return self.__prodcatg

    @prodcatg.setter
    def prodcatg(self, prodcatg_type):
        if prodcatg_type is None:
            self.__prodcatg = None
            self.__descriptor = None
            self.__header_keywords = None
        elif not isinstance(prodcatg_type, str):
            raise TypeError('The type of the input prodcatg={} is {}. '.format(prodcatg_type, type(prodcatg_type)) +
                            'Only \'str\' are accepted')
        elif prodcatg_type.upper() in PRODCATG.keys():
            self.__prodcatg = prodcatg_type.upper()
            self.__descriptor = PRODCATG[prodcatg_type.upper()]
            self.__header_keywords = _read_header_keywords_table(prodcatg_type.upper())
        else:
            raise ValueError('The input prodcatg={} is not valid. '.format(prodcatg_type.upper()) +
                             'Accepted values are: {}'.format(list(PRODCATG.keys())))

    @property
    def descriptor(self):
        r"""Dictionary that dictates the property of a `prodcatg`

        """
        return self.__descriptor

    @property
    def header_keywords(self):
        r"""Dictionary that dictates the property of the header keywords of a `prodcatg`

        """
        return self.__header_keywords

    def _get_descriptor_value(self, descriptor_key):
        """Get values corresponding to a key in the descriptor dict.

        If the `descriptor_key` is not present in the descriptor dict a warning is printed and the code returns
        a `None`

        Args:
            descriptor_key (str):  key of the descriptor dictionary

        Returns:
            any: value corresponding to the descriptor_key in the descriptor dictionary

        """
        if descriptor_key in self.descriptor.keys():
            return self.descriptor[descriptor_key]
        else:
            msgs.warning('Key {} not present in the descriptor of {}'.format(descriptor_key, self.prodcatg))
            return None

    def data_format(self):
        r"""Returns a string describing the data format associated to the `prodcatg`.

        Returns:
            str: string containing the general description of the data format of the given `prodcatg`.

        """
        descriptor_key = 'data_format'
        return self._get_descriptor_value(descriptor_key)

    def get_header_keyword_dictionary(self, header_keyword):
        r"""

        Args:
            header_keyword:

        Returns:

        """
        if header_keyword in self.header_keywords.keys():
            return self.header_keywords[header_keyword]
        else:
            msgs.warning('Keyword {} is not characterized for {}. Possible values are: \n'.format(header_keyword,
                                                                                                  self.prodcatg) +
                         '{}'.format(self.header_keywords.keys()))
            return None

    def _show_header_keyword_info(self, header_keyword):
        r"""Given the name of a header keyword, the method prints its status for the assigned `prodcatg`

        Args:
            header_keyword (str):

        """
        keyword_dictionary = self.get_header_keyword_dictionary(header_keyword)
        msgs.info('The keyword {} for a {} is:'.format(header_keyword, self.prodcatg))
        for condition_value in keyword_dictionary['condition']:
            msgs.pre_indent('{} - {} '.format(condition_value,
                                              _get_header_table_legend(condition_value)))

        return None

    def show_header_keywords_info(self, header_keywords):
        r"""Given the name of a header keywords, the method prints its status for the assigned `prodcatg`

        Args:
            header_keywords (list):

        """
        header_keywords_list = cleaning_lists.make_list_of_strings(header_keywords)
        for header_keyword in header_keywords_list:
            keyword_dictionary = self.get_header_keyword_dictionary(header_keyword)
            msgs.info('The keyword {} for a {} is:'.format(header_keyword, self.prodcatg))
            for which_header in keyword_dictionary['header']:
                msgs.pre_indent('Applicable to the {} header'.format(which_header))
            for condition_value in keyword_dictionary['condition']:
                msgs.pre_indent('{} - {} '.format(condition_value,
                                                  _get_header_table_legend(condition_value)))

        return None
