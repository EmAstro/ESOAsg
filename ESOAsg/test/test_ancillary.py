r""" Module to run tests on the ancillary module
"""

from astropy.io.fits import Header
from ESOAsg.ancillary import cleaning_headers


def test_delete_history():
    r"""Test for cleaning_headers.delete_history_keywords

    The test creates an header, fill it with HISTORY keywords and test if delete_history_keywords is able to remove them
    """
    input_header = Header({'SIMPLE': True})
    input_header.add_history('Testing the capabilities of removing HISTORY from header')
    input_header.add_history(' - HISTORY line 1')
    input_header.add_history(' - HISTORY line 2')
    input_header.add_history(' - HISTORY line 3')
    output_header = cleaning_headers.delete_history_keywords(input_header, verbose=False, in_place=False)
    assert len(output_header) == 1
