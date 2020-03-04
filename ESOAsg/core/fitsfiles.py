"""
Module to hack a file.fits
This is heavily using:
`astropy.io.fits <https://docs.astropy.org/en/stable/io/fits/api/files.html#astropy.io.fits.open>`_
"""

from astropy.io import fits
import numpy as np
import os
import re

from ESOAsg import msgs


def get_hdul(fits_name, mode='readonly', checksum=True):
    r"""
    Wrapper for astropy `fits.open`. It checks if the file exists and in case returns its HDUList.

    Args:
        fits_name (`str`):
            fits file name
        mode (`str`):
            Open mode for the file. Possibilities are: `readonly’, `update`, `append`, `denywrite`, or `ostream`.
        checksum (`bool`):
            If True, verifies that both `DATASUM` and `CHECKSUM` card values (when present in the HDU header)
            match the header and data of all HDU’s in the file. Updates to a file that already has a checksum
            will preserve and update the existing checksums unless this argument is given a value of `remove`,
            in which case the `CHECKSUM` and `DATASUM` values are not checked, and are removed when saving
            changes to the file.

    Returns:
        hdul (`hdul`):
            list-like collection of HDU objects.
    """
    if fits_name is None:
        msgs.error('No file selected.')
    elif not os.path.exists(fits_name):
        msgs.error('File does not exist.')
    else:
        hdul = fits.open(fits_name, mode=mode, checksum=checksum)
        msgs.info('The fits file {} contains {} HDUs'.format(fits_name, len(hdul)))
        msgs.info('Summary:')
        hdul.info()
        return hdul


def header_from_file(fits_name, which_hdu=0, mode='readonly', checksum=True):
    r"""
    Load an header with the information from a fits file

    Args:
        fits_name (`str`):
            fits file name
        which_hdu (`numpy.int`):
            select from which HDU you are getting the header
        mode (`str`):
            Open mode for the file. Possibilities are: `readonly’, `update`, `append`, `denywrite`, or `ostream`.
        checksum (`bool`):
            If True, verifies that both `DATASUM` and `CHECKSUM` card values (when present in the HDU header)
            match the header and data of all HDU’s in the file. Updates to a file that already has a checksum
            will preserve and update the existing checksums unless this argument is given a value of `remove`,
            in which case the `CHECKSUM` and `DATASUM` values are not checked, and are removed when saving
            changes to the file.

    Returns:
         header (`hdu.header`):
             the header corresponding to `which_hdu` from `fits_name`
    """
    hdul = get_hdul(fits_name, mode=mode, checksum=checksum)
    return hdul[which_hdu].header


def header_from_txt_file(txt_file):
    r"""Load an header from a text file into an `astropy` `hdu.header` object.

    The text file in input should contain the information in the standard header format:
        SIMPLE  =                T / Standard FITS format
        BITPIX  =               16 / # of bits storing pix values
        NAXIS   =                2 / # of axes in frame
        NAXIS1  =             2148 / # pixels/axis
        NAXIS2  =             1365 / # pixels/axis
        ORIGIN  = 'ESO'            / European Southern Observatory.
        DATE    = '2003-07-25T05:41:32.569' / UT date when this file was written
        CRVAL1  =              1.0 / value of ref pixel
        CRPIX1  =              1.0 / Ref. pixel of center of rotation
        CDELT1  =              2.0 / Binning factor
        etc..
    Cards will be read only if there is a value associated (i.e. if they are followed by a = sign). Note that
    COMMENTS will not be propagated.
    In case the file does not exist, an empty header will be returned and a warning statement will be
    raised.

    Args:
        txt_file (`str`):
            txt file name

    Returns:
         header_from_txt (`hdu.header`):
             an header object
    """

    # Checks for txt_name
    assert isinstance(txt_file, str), '`txt_name` needs to be a str'

    # creating hdu object
    hdu = fits.PrimaryHDU()
    header_from_txt = hdu.header
    if not os.path.isfile(txt_file):
        msgs.warning('File {} not exists. Returning empty header'.format(txt_file))
    else:
        with open(txt_file, 'r') as txt_header:
            for line in txt_header:
                card, value, comment = from_line_to_header_card(line)
                if card is None:
                    msgs.warning('The following line will not be added to the header\n {}'.format(line))
                else:
                    add_header_card(header_from_txt, card, value, comment=comment)
    for card in header_from_txt.cards:
        print(card)
    return header_from_txt


def from_line_to_header_card(line):  # written by Ema 04.03.2020
    r"""Given a line of text from an header, it returns card, value (and comment, if present).

    This is a tool to read in headers that have been saved in ascii format. Typically a line will be in the form:
    `DATE    = '2003-07-25T05:41:32.569' / UT date when this file was written`
    the code will try to divide the line into:
     - card = DATE
     - value = '2003-07-25T05:41:32.569'
     - comment = UT date when this file was written
     care is taken for cases in which values and/or comment contains characters like `=` or `/`
    In the possible case that `line` could no be processed, None, None, None will be returned and a warning statement
    will be raised.

    Args:
        line (`str`):
            input string to be divided into card, value, comment

    Returns:
        card, value, comment (`str`, `int`, `float`):
            if there are no comment, a `None` will be returned.
    """

    # checks for line
    assert isinstance(line, str), '`line` needs to be a str'

    if '=' not in line:
        msgs.warning('The following line could not be interpreted as header:\n {}'.format(line))
        card, value, comment = None, None, None
    else:
        # taking as card, everything that appears before the first occurrence of `=`
        card = line[0:re.search(r'=', line).start()].strip()
        line_leftover = line[re.search(r'=', line).start()+1:].strip()
        if line_leftover.count('/') == 0:
            # no comment present
            value = line_leftover
            comment = ' '
            value = check_value(value)
        elif line_leftover.count('/') == 1:
            # comment present after /
            value = line_leftover[0:re.search(r'/', line_leftover).start()].strip()
            comment = line_leftover[re.search(r'/', line_leftover).start()+1:].strip()
        else:
            test_strings = ['pixels/axis', 'PIXELS/AXIS', 'arcsec/mm', 'ARCSEC/MM',
                            'nm/mm', 'NM/MM', 'grooves/nm', 'GROOVES/NM']
            if any(test_string in line_leftover for test_string in test_strings):
                # Checking for typical values in the comments
                value = line_leftover[0:re.search(r'/', line_leftover).start()].strip()
                comment = line_leftover[re.search(r'/', line_leftover).start() + 1:].strip()
            elif line_leftover.startswith("'") and line_leftover.count("'") == 2:
                # Checking for strings in the values
                value = line_leftover[1:re.search(r"'", line_leftover[1:]).start()].strip()
                line_leftover = line_leftover[re.search(r"'", line_leftover[1:]).start()+1:].strip()
                comment = line_leftover[re.search(r'/', line_leftover).start()+1:].strip()
            else:
                # Troubles
                msgs.warning('The following line could not be interpreted as header:\n {}'.format(line))
                card, value, comment = None, None, None
        if len(comment) == 0:
            comment = None
    return card, check_value(value), comment


def check_value(value):  # written by Ema 04.03.2020
    r"""Guess for the best type of header values

    Args:
        value (`str`):
            input string value

    Returns:
        value (`str`,`int`,`float`,`bool`)
            output value with (hopefully) the correct type

    """

    if value is not None:
        if value == 'T':
            value = np.bool(True)
        elif value == 'F':
            value = np.bool(False)
        elif value.startswith("'") and value.endswith("'"):
            value = str(value[1:-1])
        elif any(character.isalpha() for character in value):
            if value.count('E') == 1 or value.count('e') == 1:
                value_test = value
                for exponent in ['e+', 'E+', 'e-', 'E-', 'e', 'E', '.']:
                    value_test = value_test.replace(exponent, '')
                    if all(character.isdigit() for character in value_test):
                        value = np.float(value)
            else:
                value = str(value)
        elif all(character.isdigit() for character in value):
            value = np.int(value)
        elif value.startswith('+') or value.startswith('-'):
            if all(character.isdigit() for character in value[1:]):
                value = np.int(value)
        else:
            value = np.float(value)

    return value


def new_fits_like(source_fits, which_hdul, output_fits, overwrite=True):
    r"""
    Create a fits file called `fits_new` that has the standard HDUL[0] and the HDUL[`which_hdul`] from
    `source_fits` appended one after the other.

    Args:
        source_fits (`str`):
            input fits file name
        which_hdul (`numpy.array`):
            select from which HDUL will be copied in the `source_fits`. It needs to be an array of integer. So
            for a single HDUL use [0], while for multiple: [0,1,2]..
        output_fits (`str`):
            output fits file name
        overwrite (`bool`):
            if `True` overwrite the `output_fits` file

    Returns:
        The code creates a new fits file with the same HDUL[0] of the input file.
    """
    source_hdul = get_hdul(source_fits, mode='readonly', checksum=True)
    output_hdu = fits.PrimaryHDU()
    output_hdul = fits.HDUList([output_hdu])
    for output_which_hdul in which_hdul:
        output_hdul.append(source_hdul[output_which_hdul])
    output_hdul.writeto(output_fits, overwrite=overwrite, checksum=True)
    source_hdul.close()


def transfer_header_cards(source_header, output_header, source_cards, output_cards=None,
                          with_comment=True, delete_card=True):
    r"""
    Copy cards, values (and optionally comments, if `with_comment`=`True`) from the header `source_header` to the
    header `output_header`. `source_cards` is a list containing all the cards you would like to transfer.
    If `output_cards` is defined, the cards from `source_cards[i]` will be saved in the `output_header` has
    `output_cards[i]`. If `delete_card`=`True` the card will be removed from the `source_header`.

    ..note ::
        Both `source_header` and `output_header` are modified in place. I.e. there is no backup option for the
        original values.

    Args:
        source_header (`hdu.header'):
            Header from which the cards will be taken.
        output_header (`hdu.header'):
            Header that will be modified with cards from `source_header`.
        source_cards (`list`):
            List of cards you want to transfer from `source_header` to `output_header`.
        output_cards (`list`):
            If not `None` the cards in `output_header` will be saved with the new names listed here.
        with_comment (`bool`):
            if true, also the associated comment will be copied
        delete_card (`bool`):
            if true, the card will be removed from the `source_header`

    Returns:
        `source_header` and `output_header` with update values.
    """
    if output_cards is None:
        output_cards = source_cards
    if len(output_cards) != len(source_cards):
        msgs.error("Incompatible length between output and source cards lists.")

    for source_card, output_card in zip(source_cards, output_cards):
        msgs.work("Transferring header card {} to {}.".format(source_card, output_card))
        if with_comment:
            add_header_card(output_header, output_card, source_header[source_card],
                            comment=source_header.comments[source_card])
        else:
            add_header_card(output_header, output_card, source_header[source_card],
                            comment=None)
        if delete_card:
            del source_header[source_card]


def add_header_card(header, card, value, comment=None):
    r"""
    Add `card` to `header` with `value`. If `comment` is not `None`, a comment will be also included.

    ..note ::
        `header` is modified in place. I.e. there is no backup option for the original value.

    Args:
        header (`hdu.header'):
            Header to which add the card.
        card (`str`):
            Card to add to the header.
        value (`str`):
            Value to associate to the card.
        comment (`str`, `None`):
            If not `None`, comment to be added to the card.

    Returns:
        `header` with update values.
    """
    if comment is None:
        header[card] = value
    else:
        header[card] = value, comment
