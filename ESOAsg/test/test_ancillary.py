r""" Module to run tests on the ancillary module
"""

from astropy.io.fits import Header
from ESOAsg.ancillary import cleaning_headers


def test_delete_history():
    r"""Test for cleaning_headers.delete_history_keywords

    The test creates an header, fill it with HISTORY keywords and test if delete_history_keywords is
    able to remove them. The option `in_place` is set to False, so the code will create a copy of the input header.
    """
    input_header = Header({'SIMPLE': True})
    input_header.add_history('Testing the capabilities of removing HISTORY from header')
    output_header = cleaning_headers.history_keywords_delete(input_header, verbose=False, in_place=False)
    assert len(output_header) == 1


def test_delete_history_in_place():
    r"""Test for cleaning_headers.delete_history_keywords

    The test creates an header, fill it with HISTORY keywords and test if delete_history_keywords is
    able to remove them. The option `in_place` is set to True, so the code will change the header in place.
    """
    input_header = Header({'SIMPLE': True})
    input_header.add_history('Testing the capabilities of removing HISTORY from header')
    _ = cleaning_headers.history_keywords_delete(input_header, verbose=False, in_place=True)
    assert len(input_header) == 1
