"""
Module to compare lists
"""

import numpy as np
from astropy.io import ascii

from ESOAsg import msgs
from ESOAsg.core import fitsfiles


def convert_to_numpy_array(input_list):
    r""" Convert an input to a numpy array
    
    Args:
        input_list:
            to be converted into a numpy array. At the moment the supported type for the input are: `list`,

    Returns:
        output_np_array (`np.array`):
            output array with the same dimension of `input_list`
    """
    if type(input_list) is list:
        output_np_array = np.copy(np.array(input_list))
        msgs.warning('Converting list to np.array with dimension {}.'.format(np.shape(output_np_array)))
    else:
        output_np_array = np.copy(input_list)
    return output_np_array


class Lists:
    r"""
    A class used to check lists.

    This allows to perform some checks on an input list (could also be given with a text files or from a fits file
    header).
    
    Attributes:
        cards (`np.array`):
            The list you want to test.
        values (`np.array`)
            If not `None` these are the values to be associated to `cards`.

    Methods:
        get_cards:
            returns the cards of the list


    """

    def __init__(self, cards=None, values=None, from_fits=None, which_hdu=0, from_txt=None,
                 data_start=None, data_end=None):
        r"""

        Given the variety of possible format of the input, particular care should be used with the `from_txt` option.
        The options `data_start` and `data_end` can help with this, but if the format is not a simple plain text
        with space/tab separated column, reading the text file with another module and then store the data with
        cards and values could be a more safe option.

        Args:
            cards (`np.array` or `list`):
                array that you want to use as cards
            values (`np.array` or `list`):
                values associated to each element of `cards`. This implies that `cards` and `values` needs to have the
                same length. But multiple `values` could be associated to `cards`.
            from_fits (`str`):
                fits file name form which read the header.
            which_hdu (`numpy.int`):
                select from which HDU you are getting the header. See `fitsfiles` in `ESOAsg.core` for further
                details.
            from_txt (`str`):
                Ascii file from which cards and values will be read. The assumption here is that `cards` are stored
                in the first column, and `values` in the second one. The options `data_start` and `data_end` will
                select only the (`data_end`-`data_start`) lines following the line `data_start`. This is the same
                option present in `ascii.read` from astropy.
            data_start (`np.int`):
                First line to be read from the file: `from_txt`.
            data_end (`np.int`):
                Last line to be read from the file: `from_txt`.

        Returns:
            a list object containing cards and corresponding values
        """

        # Loading from input
        if cards is not None:
            # Loading cards
            msgs.work('Loading cards.')
            _cards = convert_to_numpy_array(cards)
            self.cards = _cards
            # Loading values
            if values is not None:
                _values = convert_to_numpy_array(values)
                if np.shape(_values) == np.shape(_cards):
                    msgs.work('Loading values.')
                    self.values = _values
                else:
                    if np.ndim(_values) == 1:
                        if np.shape(_values[:]) == np.shape(_cards):
                            self.values = _values
                        else:
                            msgs.error('Values and Cards should have the same length.')
                    elif np.ndim(_values) == 2:
                        for index in np.arange(np.shape(_values)[0], dtype=np.int_):
                            if np.shape(_values[index, :]) != np.shape(_cards):
                                msgs.error('Values and Cards should have the same length.')
                        self.values = _values
                    else:
                        msgs.error('Values and Cards should have the same length.')
            else:
                msgs.work('Values is empty.')
                self.values = np.array([])

        # Loading from fits file
        elif from_fits is not None:
            msgs.work('Loading header from fits file: {}'.format(from_fits))
            _hdu = fitsfiles.header_from_file(from_fits, which_hdu=which_hdu)
            _cards, _values = np.zeros_like(0), np.zeros_like(0)
            for index in list(_hdu.keys()):
                if 'COMMENT' not in index:
                    _cards = np.append(_cards, index)
                    _values = np.append(_values, _hdu[index])
            self.cards = np.array(_cards)
            self.values = np.array(_values)

        # Loading from text file
        elif from_txt is not None:
            msgs.work('Loading list from text file: {}'.format(from_txt))
            _full_table = ascii.read(from_txt, guess=True, data_start=data_start, data_end=data_end)
            self.cards = np.array(_full_table[_full_table.colnames[0]].data)
            if np.size(_full_table.colnames) < 2:
                self.values = np.array([])
            elif np.size(_full_table.colnames) >= 2:
                _values = np.array(_full_table[_full_table.colnames[1]].data)
                if np.size(_full_table.colnames) > 2:
                    for index in _full_table.colnames[2:]:
                        _values = np.vstack((_values, _full_table[index].data))
                self.values = _values

        # Creating empty object
        else:
            msgs.work('Creating empty lists.Lists object.')
            self.cards = np.array([])
            self.values = np.array([])

        return

    def get_cards(self, check_cards=None):
        r"""Returns the `cards` present in a list object. If `cards` is not `None`, it checks that such cards are
        contained in the list and return them. If not present an error is raised.

        Args:
            check_cards (`np.array`):
                List of cards that needs to be checked.

        Returns:
            selected_cards (`np.array`):
                Cards present in the list.
            missing_cards (`np.array`):
                Subset of `check_cards` not present in the list. It is empty if all cards are present.
        """
        if check_cards is not None:
            _check_cards = np.isin(self.cards, convert_to_numpy_array(check_cards))
            selected_cards = np.copy(self.cards[_check_cards])
            if np.size(selected_cards) < 1:
                msgs.error('None of the cards in input are present in the list.')
            else:
                msgs.info('There are {} occurrences of the {} cards in input.'.format(np.size(selected_cards),
                                                                                      np.size(check_cards)))
                if np.all(np.isin(check_cards, selected_cards)):
                    msgs.info('All cards in input are present in the list.')
                    missing_cards = np.array([])
                else:
                    msgs.info('Not all the cards in input are present in the list.')
                    missing_cards = check_cards[np.logical_not(np.isin(check_cards, selected_cards))]
                    msgs.info('The missing cards are:')
                    for missing_card in missing_cards:
                        msgs.info(' - {}'.format(missing_card))
        else:
            selected_cards = np.copy(self.cards)
            missing_cards = np.array([])
        return selected_cards, missing_cards

    def get_values(self, cards=None):
        r"""

        Args:
            cards:

        Returns:

        """
        return self.get_cards(self, cards=cards), np.array(self.values)
