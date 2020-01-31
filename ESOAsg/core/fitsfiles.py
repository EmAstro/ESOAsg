"""
Module to hack a file.fits
This is heavily using:
`astropy.io.fits <https://docs.astropy.org/en/stable/io/fits/api/files.html#astropy.io.fits.open>`_
"""

from astropy.io import fits
from os import path

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
    elif not path.exists(fits_name):
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
