r""" Module that performs some simple operations on headers

.. include:: ./include/links.rst

"""

from ESOAsg import msgs
from ESOAsg.ancillary import cleaning_lists
from ESOAsg.ancillary import checks


def delete_history_keywords(input_header, in_place=True, verbose=True):
    r"""Remove `HISTORY` keywords from an header

    Given an header as input, all the occurrences of the `HISTORY` keywords are removed. The default is to modify the
    header in place. If `in_place` is set to `False` the code will create copy of the header and return it.

    Args:
        input_header (`header`_): header from which the `HISTORY` keywords will be removed
        in_place (bool): if set to `True` the header will be modified in place, while if set to `False` a copy will be
            created
        verbose (bool): print on screen additional information on the process

    Returns:
        `header`_: header with history keywords removed
    """
    if not checks.header_is_valid(input_header):
        msgs.warning('The input header is not valid')
        return input_header

    # Primary Header
    if 'HISTORY' in input_header.keys():
        if verbose:
            history_cards = [history_card for history_card in input_header
                             if history_card.startswith('HISTORY')]
            history_values = [input_header[history_card] for history_card in input_header
                              if history_card.startswith('HISTORY')]
            for history_card, history_value in zip(history_cards, history_values):
                msgs.work('Cleaning cards: {} = {}'.format(history_card, history_value))
        if in_place:
            del input_header['HISTORY']
            return input_header
        else:
            output_header = input_header.copy()
            del output_header['HISTORY']
            return output_header
    else:
        if verbose:
            msgs.info('No HISTORY cards in the header')
        return input_header


def clean_non_ascii_history_keywords(header):
    r"""

    Args:
        header:

    Returns:

    """

    # Data Header
    if 'HISTORY' in header.keys():
        history_values_hdr1 = header['HISTORY'][:]
        for history_number in range(0, len(history_values_hdr1)):
            clean_history = cleaning_lists.remove_non_ascii(history_values_hdr1[history_number])
            if len(clean_history) > 0:
                header['HISTORY'][history_number] = str(clean_history)
            else:
                header['HISTORY'][history_number] = str(' ')
    return header
