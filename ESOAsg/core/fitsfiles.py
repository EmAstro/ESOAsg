"""
Module to hack a file.fits
"""

from astropy.io import fits
from os import path

from ESOAsg import msgs


def get_hdul(fits_name):
    r"""
    Wrapper for astropy `fits.open`

    Args:
        fits_name (`str`):
            fits file name
    Returns:
        hdul object
    """
    if fits_name is None:
        msgs.error('No file selected.')
    elif not path.exists(fits_name):
        msgs.warning('File does not exist.')
    else:
        hdul = fits.open(fits_name)
        msgs.info('The fits file {} contains {} HDUs'.format(fits_name, len(hdul)))
        return hdul


def header_from_file(fits_name, which_hdu=0):
    r"""
    Load an header with the information from a fits file

    Args:
        fits_name (`str`):
            fits file name
        which_hdu (`numpy.int`):
            select from which HDU you are getting the header

    Returns:
         header (`fits.header`):
             the header corresponding to `which_hdu` from `fits_name`
    """
    hdul = get_hdul(fits_name)
    return hdul[which_hdu].header


def new_fits_like(fits_original, which_hdul, fits_new, overwrite=True):
    r"""
    Create a fits file called `fits_new` that has the standard HDUL[0] and the HDUL[`which_hdul`] from
    `fits_original` appended one after the other.

    Args:
        fits_original (`str`):
            input fits file name
        which_hdul (`numpy.array`):
            select from which HDUL will be copied in the `fits_new`
        fits_new (`str`):
            output fits file name
        overwrite (`bool`):
            if `True` overwrite the `fits_new` file

    Returns:
        The code creates a new fits file with the same HDUL[0] of the input file.
    """
    hdul_original = get_hdul(fits_original)
    hdu_new = fits.PrimaryHDU()
    hdul_new = fits.HDUList([hdu_new])
    for which_hdul_new in which_hdul:
        hdul_new.append(hdul_original[which_hdul_new])
    hdul_new.writeto(fits_new, overwrite=overwrite, checksum=True)
    hdul_original.close()


def transfer_header_cards(source_header, output_header, source_cards, output_cards=None,
                          with_comment=True, delete_card=True):
    r"""

    Args:
        source_header
        output_header
        source_cards
        output_cards
        with_comment (`bool`):
            if true, also the associated comment will be copied
        delete_card (`bool`):
            if true, the card will be removed from the `source_header`

    Returns:

    """
    if output_cards is None:
        output_cards = source_cards

    for source_card, output_card in zip(source_cards, output_cards):
        msgs.info("Transferring header card {} to {}.".format(source_card, output_card))
        if with_comment:
            add_header_card(output_header, output_card, source_header[source_card], comment=source_header.comments[source_card])
        else:
            add_header_card(output_header, output_card, source_header[source_card], comment=None)
        if delete_card:
            del source_header[source_card]
    return output_header


def add_header_card(header, card, value, comment=None):
    r"""

    Args:
        header
        card
        value
        comment

    Returns:

    """
    if comment is None:
        header[card] = value
    else:
        header[card] = value, comment

