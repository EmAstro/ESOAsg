"""
Module to compare lists
"""


import numpy as np
import os
from astropy.io import ascii

from ESOAsg import msgs
from ESOAsg.core import fitsfiles


def _convert_to_numpy_array(input_list):
    r""" Convert an input to a numpy array
    
    Args:
        input_list:
            to be converted into a numpy array. At the moment the supported type for the input are: `list`,

    Returns:
        output_np_array (`np.array`):
            output array with the same dimension of `input_list`
    """
    if isinstance(input_list, list):
        output_np_array = np.copy(np.array(input_list))
        msgs.work('Converting list to np.array with dimension {}.'.format(np.shape(output_np_array)))
    else:
        output_np_array = np.copy(input_list)
    return output_np_array


def _find_max_str_length(vector):
    r"""Given a `np.array` finds the length of the longest element in it.

    Args:
        vector (`np.array`):
            input `np.array`
    Returns:
        max_length (`np.int`):
            maximum length of the string in an array
    """
    if np.size(vector) == 0:
        max_length = np.int(0)
    else:
        max_length = len(max(vector, key=len))
    return max_length


def _get_y_dimension(array):
    r"""Returns the y-dimension of a 2D-array. If is 1D it returns 1.

    Args:
        array (`np.array`):
            input array

    Returns:
        y_dimension (`np.int`):
            size of the y dimension of the array

    """
    if np.ndim(array) == 1:
        y_dimension = 1
    else:
        y_dimension, _ = np.shape(array)
    return y_dimension


def _print_1cards_2values(cards, values1, values2, on_terminal=True, on_file=None):
    r""" Print on terminal (if `on_terminal` is `True`) and/or on a text file (if `on_file` is not `None`) cards
    and values for

    Args:
        cards (`np.array`):
            cards to be printed
        values1 (`np.array`):
            first values to be printed
        values2 (`np.array`):
            second values to be printed
        on_terminal (`bool`):
            if `True`, the `cards` and the two `values` are printed on the terminal.
        on_file (`str`, `None`):
            if not `None`, `cards` and `values` will be stored in this text file

    Returns:
        Nothing is returned, but if `on_file` is set, `cards` and `values` will be saved on a text file, and if
        `on_terminal` is `True` it will be also displayed on the terminal.
    """
    # the maximum number of chars required to display any item in list
    max_length = np.max((_find_max_str_length(cards), _find_max_str_length(values1),
                         _find_max_str_length(values2), 10))
    if np.ndim(cards) != 1:
        msgs.error('Cards should be a 1D array.')
    else:
        spaces_values1 = _get_y_dimension(values1)
        spaces_values2 = _get_y_dimension(values2)
        header_line = 'CARDS' + ' '*(max_length-4) + 'VALUES_1st' + ' '*((spaces_values1-1)*(max_length+1)) + \
                      ' '*(max_length-9) + 'VALUES_2nd' + ' '*((spaces_values2-1)*(max_length+1)) + \
                      ' '*(max_length-9)+'\n' + '-----' + ' '*(max_length-4)+'----------' + \
                      ' '*((spaces_values1-1)*(max_length+1)) + ' '*(max_length-9) + '----------'
        complete_matrix = np.vstack((cards, values1, values2))
        if on_file is None:
            filename = '_temp.txt'
        else:
            filename = on_file
        np.savetxt(filename, complete_matrix.T, fmt='%-'+str(max_length)+'s', header=header_line)
        if on_terminal:
            print(" ")
            with open(filename) as text_file:
                line_number = np.int(0)
                for line in text_file:
                    if line_number < 2:
                        print(line.strip())
                    else:
                        print('  '+line.strip())
                    line_number += 1
            print(" ")
        os.remove(filename) if os.path.exists(filename) and on_file is None else None


def _print_1cards_1values(cards, values1, on_terminal=True, on_file=None):
    r""" Print on terminal (if `on_terminal` is `True`) and/or on a text file (if `on_file` is not `None`) cards
    and values for

    Args:
        cards (`np.array`):
            cards to be printed
        values1 (`np.array`):
            first values to be printed
        on_terminal (`bool`):
            if `True`, the `cards` and the `values1` are printed on the terminal.
        on_file (`str`, `None`):
            if not `None`, `cards` and `values1` will be stored in this text file

    Returns:
        Nothing is returned, but if `on_file` is set, `cards` and `values1` will be saved on a text file, and if
        `on_terminal` is `True` it will be also displayed on the terminal.
    """
    # the maximum number of chars required to display any item in list
    max_length = np.max((_find_max_str_length(cards), _find_max_str_length(values1), 10))
    if np.ndim(cards) != 1:
        msgs.error('Cards should be a 1D array.')
    else:
        spaces_values1 = _get_y_dimension(values1)
        header_line = 'CARDS' + ' '*(max_length-4) + 'VALUES' + ' '*((spaces_values1-1)*(max_length+1)) + \
                      ' '*(max_length-10)+'\n' + '-----' + ' '*(max_length-4)+'------' + \
                      ' '*((spaces_values1-1)*max_length)
        complete_matrix = np.vstack((cards, values1))
        if on_file is None:
            filename = '_temp.txt'
        else:
            filename = on_file
        np.savetxt(filename, complete_matrix.T, fmt='%-'+str(max_length)+'s', header=header_line)
        if on_terminal:
            print(" ")
            with open(filename) as text_file:
                line_number = np.int(0)
                for line in text_file:
                    if line_number < 2:
                        print(line.strip())
                    else:
                        print('  '+line.strip())
                    line_number += 1
            print(" ")
        os.remove(filename) if os.path.exists(filename) and on_file is None else None


class Lists:
    r"""A class used to check lists.

    This allows to perform some checks on an input list (could also be given with a text files or from a fits file
    header). A Lists object is constituted by an array of cards and an array of corresponding values. Multiple values
    can be assigned to the same cards. But the length of values should be the same of the length of cards.
    For instance:
    Cards:      Values1:    Values2:
    -----       -------     -------
    JE          ICH             IO
    TU          DU              TU
    ELLES       SIE         LEI
    NOUS        WIR             NOI
    VOUS        IHR             VOI
    ELLES       SIE             ESSE
    -----       -------     -------

    Attributes:
        cards (`np.array`):
            The list you want to test.
        values (`np.array`)
            If not `None` these are the values to be associated to `cards`.

    Methods:
        get_cards:
            returns the cards from a list
        get_values:
            returns the values from a list
        compare_with:
            compare two lists
        difference_with:
            find differences between two lists

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
            from_txt (`str`, `None`):
                Ascii file from which cards and values will be read. The assumption here is that `cards` are stored
                in the first column, and `values` in the second one. The options `data_start` and `data_end` will
                select only the (`data_end`-`data_start`) lines following the line `data_start`. This is the same
                option present in `ascii.read` from astropy.
            data_start (`np.int`, `None`):
                First line to be read from the file: `from_txt`.
            data_end (`np.int`, `None`):
                Last line to be read from the file: `from_txt`.

        Returns:
            a list object containing cards and corresponding values
        """

        # Loading from input
        if cards is not None:
            # Loading cards
            msgs.work('Loading cards.')
            _cards = _convert_to_numpy_array(cards)
            self.cards = _cards
            # Loading values
            if values is not None:
                _values = _convert_to_numpy_array(values)
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
            _hdu = fitsfiles.header_from_fits_file(from_fits, which_hdu=which_hdu)
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
        r"""Returns the `cards` present in a list object. If `check_cards` is not `None`, it checks that such cards are
        contained in the list and return them. If not present an empty `np.vector` is returned.

        Args:
            check_cards (`np.array`, `None`):
                List of cards that needs to be checked.

        Returns:
            selected_cards (`np.array`):
                Cards present in the list.
            missing_cards (`np.array`):
                Subset of `check_cards` not present in the list. It is empty if all cards are present.
        """
        if check_cards is not None:
            _check_cards = np.isin(self.cards, _convert_to_numpy_array(check_cards))
            selected_cards = np.copy(self.cards[_check_cards])
            if np.size(selected_cards) < 1:
                msgs.warning('None of the cards in input are present in the list.')
                selected_cards = np.array([])
                missing_cards = np.array([])
            else:
                msgs.work('There are {} occurrences of the {} cards in input.'.format(np.size(selected_cards),
                                                                                      np.size(check_cards)))
                if np.all(np.isin(check_cards, selected_cards)):
                    msgs.work('All cards in input are present in the list.')
                    missing_cards = np.array([])
                else:
                    msgs.work('Not all the cards in input are present in the list.')
                    missing_cards = check_cards[np.logical_not(np.isin(check_cards, selected_cards))]
                    msgs.work('The missing cards are:')
                    for missing_card in missing_cards:
                        msgs.work(' - {}'.format(missing_card))
        else:
            selected_cards = np.copy(self.cards)
            missing_cards = np.array([])
        return selected_cards, missing_cards

    def get_values(self, check_cards=None):
        r"""Returns the values associated to the cards in a list. If `check_cards` is defined, only the subset of values
        selected will be returned.

        Args:
            check_cards (`np.array`):
                List of cards that needs to be returned.

        Returns:
            selected_cards (`np.array`):
                Cards present in the list.
            selected_values (`np.array`):
                Values associated to `check_cards`
            missing_cards (`np.array`):
                Subset of `check_cards` not present in the list. It is empty if all cards are present.
        """
        selected_cards, missing_cards = self.get_cards(check_cards=check_cards)
        if check_cards is not None:
            _check_cards = np.isin(self.cards, _convert_to_numpy_array(check_cards))
            if np.size(self.values) == 0:
                msgs.warning('No values present in the list.')
                selected_values = np.array([])
            else:
                if np.ndim(self.values) < 2:
                    selected_values = np.copy(self.values[_check_cards])
                else:
                    selected_values = np.copy(self.values[0, :][_check_cards])
                    for index in np.arange(1, np.size(self.values[:, 0])):
                        selected_values = np.vstack((selected_values, np.copy(self.values[index, :][_check_cards])))
        else:
            selected_values = np.array([])
        return selected_cards, selected_values, missing_cards

    def compare_with(self, second_list, check_cards=None, on_terminal=True, on_file=None):
        r"""Compare two list objects. The code will go through all cards and values present in the
        input list (or the `check_cards` subset if set not to None) and compare with cards and values in the
        `second_list`.

        Args:
            second_list (`Lists`):
                List that you want to compare with the input list.
            check_cards (`np.array`):
                List of cards that needs to be returned.
            on_terminal (`bool`):
                if `True`, the `cards` and respective `values` are printed on the terminal.
            on_file (`str`, `None`):
                if not `None`, `cards` and respective `values` will be stored in this text file

        Returns:
            first_cards (`np.array`):
                `cards` present in both the 1st and the 2nd lists.
            first_values (`np.array`):
                `values` from the 1st list associated to the `cards`
            second_values (`np.array`):
                `values` from the 2nd list associated to the `cards`
        """

        if not isinstance(second_list, Lists):
            msgs.error('The second list is not a Lists object.')
        if not isinstance(on_terminal, bool):
            msgs.error('The on_terminal option should be a bool.')

        if check_cards is None:
            check_cards, _ = self.get_cards(check_cards)

        # loading the first list:
        _first_cards, _first_values, _first_missing = self.get_values(check_cards=check_cards)
        msgs.work('Checking for cards in the second_list')
        # comparing with values from the second list:
        _second_cards, _second_values, _second_missing = second_list.get_values(check_cards=check_cards)

        # define where cards are present in both lists
        _overlap = np.isin(_first_cards, _second_cards)
        first_cards, first_values, first_missing = self.get_values(check_cards=_first_cards[_overlap])
        second_cards, second_values, second_missing = second_list.get_values(check_cards=_first_cards[_overlap])

        # Printing results
        if on_terminal:
            msgs.info('Cards and values from both the first and the second list:')
        _print_1cards_2values(first_cards, first_values, second_values, on_terminal=on_terminal, on_file=on_file)

        if on_terminal:
            # Case of all cards present in the first and in the second lists
            if np.size(_first_missing) == 0 and np.size(_second_missing) == 0:
                msgs.info('All cards are present in both lists')
            # Case of some cards missing in the first list but present in the second:
            elif np.size(_first_missing) != 0 and np.size(_second_missing) == 0:
                _only_second = np.isin(check_cards, _first_missing)
                _only_second_cards, _only_second_values, _only_second_missing = second_list.get_values(
                    check_cards=check_cards[_only_second])
                msgs.warning('{} Cards are missing in the first list,'.format(np.size(_first_missing)))
                msgs.warning('but present in the second:')
                _print_1cards_1values(_only_second_cards, _only_second_values, on_terminal=on_terminal, on_file=None)
            # Case of some cards missing in the second list but present in the first:
            elif np.size(_first_missing) == 0 and np.size(_second_missing) != 0:
                _only_first = np.isin(check_cards, _second_missing)
                _only_first_cards, _only_first_values, _only_first_missing = self.get_values(
                    check_cards=check_cards[_only_first])
                msgs.warning('{} Cards are missing in the second list,'.format(np.size(_second_missing)))
                msgs.warning('but present in the first:')
                _print_1cards_1values(_only_first_cards, _only_first_values, on_terminal=on_terminal, on_file=None)
            # Case of some cards missing in the first and in the second list:
            else:
                _first_missing = ~np.isin(check_cards, _first_cards)
                _second_missing = ~np.isin(check_cards, _second_cards)
                _both_missing = _first_missing * _second_missing
                _first_only_missing = _first_missing * ~_both_missing
                _second_only_missing = _second_missing * ~_both_missing
                if np.size(check_cards[_both_missing]) != 0:
                    msgs.warning('{} Cards are missing from both lists:'.format(np.size(check_cards[_both_missing])))
                    for both_missing in check_cards[_both_missing]:
                        print(' ----------> '+both_missing)
                if np.size(check_cards[_second_only_missing]) != 0:
                    _only_second_cards, _only_second_values, _only_second_missing = self.get_values(
                        check_cards=check_cards[_second_only_missing])
                    msgs.warning('{} Cards are missing in the first list,'.format(
                        np.size(check_cards[_second_only_missing])))
                    msgs.warning('but present in the second:')
                    _print_1cards_1values(_only_second_cards, _only_second_values, on_terminal=on_terminal,
                                          on_file=None)
                if np.size(check_cards[_first_only_missing]) != 0:
                    _only_first_cards, _only_first_values, _only_first_missing = second_list.get_values(
                        check_cards=check_cards[_first_only_missing])
                    msgs.warning('{} Cards are missing in the second list,'.format(
                        np.size(check_cards[_first_only_missing])))
                    msgs.warning('but present in the first:')
                    _print_1cards_1values(_only_first_cards, _only_first_values, on_terminal=on_terminal,
                                          on_file=None)

        return first_cards, first_values, second_values

    def difference_with(self, second_list, check_cards=None, on_terminal=True, on_file=None):
        r"""Finds differences between two list objects. The code will go through all cards and values present in the
        input list (or the `check_cards` subset if set not to None) and compare with cards and values in the
        `second_list`.

        Args:
            second_list (`Lists`):
                List that you want to compare with the input list.
            check_cards (`np.array`):
                List of cards that needs to be returned.
            on_terminal (`bool`):
                if `True`, the `cards` and the respective `values` are printed on the terminal.
            on_file (`str`, `None`):
                if not `None`, `cards` and respective `values` will be stored in this text file

        Returns:
            diff_first_cards (`np.array`):
                `cards` present in both the 1st and the 2nd lists that that have different values
            diff_first_values (`np.array`):
                `values` from the 1st list associated to the `cards`
            diff_second_cards (`np.array`):
                `values` from the 2nd list associated to the `cards`
        """
        first_cards, first_values, second_values = self.compare_with(second_list, check_cards=check_cards,
                                                                     on_terminal=False, on_file=None)
        if np.ndim(first_values) == 1 and np.ndim(second_values) == 1:
            _differences_1st_2nd = np.isin(first_values, second_values)
            diff_first_cards, diff_first_values, diff_first_missing = self.get_values(
                check_cards=first_cards[~_differences_1st_2nd])
            diff_second_cards, diff_second_values, diff_second_missing = second_list.get_values(
                check_cards=first_cards[~_differences_1st_2nd])
            _print_1cards_2values(diff_first_cards, diff_first_values, diff_second_values,
                                  on_terminal=on_terminal, on_file=on_file)
            return diff_first_cards, diff_first_values, diff_second_values
        else:
            print('WORK IN PROGRESS')

