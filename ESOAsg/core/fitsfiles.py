"""
Module to hack a file.fits
This is heavily using:
`astropy.io.fits <https://docs.astropy.org/en/stable/io/fits/api/files.html#astropy.io.fits.open>`_
"""

from astropy.io import fits
import numpy as np
import os
import re
import ast
# import time

from ESOAsg import msgs
from ESOAsg.ancillary import checks


def get_hdul(fits_name, mode='readonly', checksum=True):  # Written by Ema 05.03.2020
    r"""Wrapper for astropy `fits.open`. It checks if the file exists and in case returns its HDUList.

    Args:
        fits_name (`str`):
            fits file name
        mode (`str`):
            Open mode for the file. Possibilities are: `readonly’, `update`, `append`, `denywrite`, or `ostream`
        checksum (`bool`):
            If True, verifies that both `DATASUM` and `CHECKSUM` card values (when present in the HDU header)
            match the header and data of all HDU’s in the file. Updates to a file that already has a checksum
            will preserve and update the existing checksums unless this argument is given a value of `remove`,
            in which case the `CHECKSUM` and `DATASUM` values are not checked, and are removed when saving
            changes to the file

    Returns:
        hdul (`hdul`):
            list-like collection of HDU objects

    """
    if not checks.fits_file_is_valid(fits_name):
        msgs.error('Fits file not valid')
    else:
        hdul = fits.open(fits_name, mode=mode, checksum=checksum)
        msgs.info('The fits file {} contains {} HDUs'.format(fits_name, len(hdul)))
        msgs.info('Summary:')
        hdul.info()
        return hdul


def header_from_fits_file(fits_name, which_hdu=0):  # Written by Ema 05.03.2020
    r"""Load an header with the information from a fits file

    Args:
        fits_name (`str`):
            fits file name
        which_hdu (`numpy.int`):
            select from which HDU you are getting the header. Default = 1

    Returns:
         header (`hdu.header`):
             the header corresponding to `which_hdu` from `fits_name`

    """

    assert isinstance(which_hdu, (int, np.int_)), 'which_hdu must be an int'
    if not checks.fits_file_is_valid(fits_name):
        msgs.error('Fits file not valid')
    else:
        header = fits.getheader(fits_name, which_hdu)
        return header


def header_from_txt_file(txt_file):  # written by Ema 05.03.2020
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
    Cards will be read only if there is a value associated (i.e. if they are followed by a = sign). In case the file
    does not exist, an empty header will be returned and a warning statement will be raised.

    Args:
        txt_file (`str`):
            txt file name

    Returns:
         header_from_txt (`hdu.header`):
             an header object

    """

    # Checks for txt_name
    assert isinstance(txt_file, (str, np.str)), '`txt_name` needs to be a str'

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
                    if 'END' not in line and 'END\n' not in line:
                        msgs.warning('The following line will not be added to the header\n {}'.format(line))
                else:
                    add_header_card(header_from_txt, card, value, comment=comment)

    return header_from_txt


def from_line_to_header_card(line):  # written by Ema 05.03.2020
    r"""Given a line of text from an header, it returns card, value (and comment, if present).

    This is a tool to read in headers that have been saved in ascii format. Typically a line will be in the form:
    `DATE    = '2003-07-25T05:41:32.569' / UT date when this file was written`
    the code will try to divide the line into:
     - card = DATE
     - value = '2003-07-25T05:41:32.569'
     - comment = UT date when this file was written
     Care is taken for cases in which values and/or comment contains characters like `=` or `/`
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
    assert isinstance(line, (str, np.str)), '`line` needs to be a str'

    if '=' not in line:
        if 'END' in line and 'END\n' not in line:
            msgs.warning('The following line could not be interpreted as header:\n {}'.format(line))
        card, value, comment = None, None, None
    else:
        # taking as card, everything that appears before the first occurrence of `=`
        card, line_leftover = re.split("=", line, maxsplit=1)
        card = card.strip()
        # now check how many occurrences of ` \ ` are in line_leftover and split values from comments
        if line_leftover.count(' / ') == 0:
            # good splitting and no comment present
            value, comment = line_leftover.strip(), None
        elif line_leftover.count(' / ') == 1:
            # good splitting
            value, comment = re.split(" / ", line_leftover, maxsplit=1)
            value, comment = value.strip(), comment.strip()
            if len(comment) == 0:
                comment = None
        else:
            # Troubles but fingers crossed
            value, comment = re.split(" / ", line_leftover, maxsplit=1)
            value, comment = value.strip(), comment.strip()
            msgs.warning('The following line should be double checked:\n {}'.format(line))

    return card, check_value(value), comment


def check_value(value):  # written by Ema 05.03.2020
    r"""Guess for the best type of header values.

    This is based on `ast.literal_eval`
    Args:
        value (`str`):
            input string value

    Returns:
        value (`str`,`int`,`float`,`bool`)
            output value with (hopefully) the correct type

    """
    special_char = re.compile('[@_!#$%^&*()<>?/\\|}{~:]')
    if value is not None:
        if value == 'T':
            value = np.bool(True)
        elif value == 'F':
            value = np.bool(False)
        elif special_char.search(value) is not None:
            value = str(value)
        else:
            value = ast.literal_eval(value)
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
    r"""Transfer header cards from one header to another

    Cards, values (and optionally comments, if `with_comment`=`True`) from the header `source_header`  will be
    transfer to the header `output_header`.
    `source_cards` is a list containing all the cards that needs to be transfer. If `output_cards` is defined, the cards
    from `source_cards[i]` will be saved in the `output_header` has `output_cards[i]`. If `delete_card`=`True` the
    card will be removed from the `source_header`.

    ..note ::
        Both `source_header` and `output_header` are modified in place. I.e. there is no backup option for the
        original values.
        If a card is not present in `source_header` a warning is raised and will not be transferred to `output_header`

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
        if source_card not in source_header:
            msgs.warning('{} not present in `source_header`. The card will not be transferred'.format(source_card))
            continue
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
