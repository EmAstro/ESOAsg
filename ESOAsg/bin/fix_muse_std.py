#!/usr/bin/env python3

import argparse

import numpy as np

from astropy.io import fits
from astropy.coordinates import name_resolve
from astropy.coordinates import FK5
# from astropy import units as u

from ESOAsg import msgs
from ESOAsg.core import fitsfiles
# from ESOAsg.core import download_archive

# from IPython import embed


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=r"""
        This macro make the spectra of standard stars collected with MUSE compatible 
        with the Phase3 standard. For more details on the project, please see:
        `https://www.eso.org/sci/publications/messenger/archive/no.178-dec19/messenger-no178-17-18.pdf`_      

        This uses ESOAsg version {:s}
        """.format(msgs._version),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument('-i', '--input_fits', nargs='+', type=str, default=None,
                        help='Original fits file')
    parser.add_argument('-o', '--output_fits', nargs='+', type=str, default=None,
                        help='Output fits file with corrected header and structure')
    parser.add_argument('-v', '--version', action='version', version=msgs._version)
    return parser.parse_args()


EXAMPLES = r"""
        Example:
        fix_muse_std.py --input_fits 
        """

if __name__ == '__main__':
    args = parse_arguments()

    # getting fits names
    input_fits = np.str(args.input_fits[0])
    if args.output_fits is None:
        output_fits = input_fits.replace('.fits', '_fixed.fits')
    else:
        output_fits = np.str(args.output_fits[0])

    msgs.start()

    # Copy relevant information from input into output file
    fitsfiles.new_fits_like(input_fits, [0], output_fits, overwrite=True)

    # This file will be modified in place
    hdul = fits.open(output_fits, 'update', checksum=True)
    hdr0 = hdul[0].header
    hdr1 = hdul[1].header

    msgs.work('Refactoring data structure.')
    FLUX = hdul[1].data
    WAVE = hdr1['CRVAL1'] + (hdr1['CDELT1']*np.arange(0., len(FLUX)))

    # double precision:
    col_format = 'D'
    col_dtype = np.float64

    col1 = fits.Column(name='WAVE', format=col_format, unit='angstrom',
                       array=np.array(WAVE, dtype=col_dtype))
    col2 = fits.Column(name='FLUX', format=col_format, unit='erg.s**(-1).cm**(-2).angstrom**(-1)',
                       array=np.array(FLUX, dtype=col_dtype))
    col3 = fits.Column(name='ERR', format=col_format, unit='erg.s**(-1).cm**(-2).angstrom**(-1)',
                       array=np.array(FLUX-FLUX, dtype=col_dtype))

    hdul[1] = fits.BinTableHDU.from_columns([col1, col2, col3])

    msgs.work('Refactoring Header.')

    # Now working on HDU1
    change1_cards, change1_values, change1_comments = [], [], []
    # Object names
    object_name: str = input_fits.replace('_av.fits', '').split('/')[-1]
    # few adjustments
    if object_name.startswith('HD'):
        object_name = object_name.replace('HD', 'HD ')
    elif object_name.startswith('IRAS'):
        object_name = object_name.replace('IRAS', 'IRAS ')
        object_name = object_name.replace('p', '+')
    elif object_name.startswith('B86'):
        object_name = object_name.replace('B86', '[B86] ')
    else:
        pass
    change1_cards.append('OBJECT')
    change1_values.append(object_name)
    change1_comments.append('Target designation')
    msgs.warning('The object name is taken from the file_name and is: {}'.format(object_name))

    # RA Dec
    radec = name_resolve.get_icrs_coordinates(object_name).transform_to(FK5(equinox='J2000'))
    msgs.warning('The coordinates of the object are taken from an online archive:')
    msgs.warning('RA={}, Dec={}'.format(radec.ra, radec.dec))
    change1_cards.append('RA')
    change1_values.append(radec.ra.value)
    change1_comments.append('[deg] Spectroscopic target position (J2000.0)')
    change1_cards.append('DEC')
    change1_values.append(radec.dec.value)
    change1_comments.append('[deg] Spectroscopic target position (J2000.0)')

    # EXTNAME
    change1_cards.append('EXTNAME')
    change1_values.append('SPECTRUM')
    change1_comments.append('FITS Extension name ')

    """
    # NAXIS1
    change1_cards.append('NAXIS1')
    change1_values.append(2)
    change1_comments.append('length of dimension 1')
    
    # NAXIS2
    change1_cards.append('NAXIS2')
    change1_values.append(1)
    change1_comments.append('length of dimension 2')
    """

    # TUCD1
    msgs.warning("AIR WAVE?!?!?!?!?!")
    change1_cards.append('TUCD1')
    change1_values.append('em.wl;obs.atmos')
    change1_comments.append('UCD of field 1')
    # TUTYP1
    change1_cards.append('TUTYP1')
    change1_values.append('spec:Data.SpectralAxis.Value')
    change1_comments.append(' ')
    # TDMIN1
    change1_cards.append('TDMIN1')
    change1_values.append(np.nanmin(WAVE))
    change1_comments.append('Start in spectral coordinates')
    # TDMAX1
    change1_cards.append('TDMAX1')
    change1_values.append(np.nanmax(WAVE))
    change1_comments.append('Stop in spectral coordinates')

    # TUCD2
    change1_cards.append('TUCD2')
    change1_values.append('phot.flux.density;em.wl;meta.main')
    change1_comments.append('UCD of field 2')
    # TUTYP2
    change1_cards.append('TUTYP2')
    change1_values.append('spec:Data.FluxAxis.Value')
    change1_comments.append(' ')

    # TUCD3
    change1_cards.append('TUCD3')
    change1_values.append('stat.error;phot.flux.density;em.wl;meta.main')
    change1_comments.append('UCD of field 3')
    # TUTYP3
    change1_cards.append('TUTYP3')
    change1_values.append('spec:Data.FluxAxis.Accuracy.StatError')
    change1_comments.append(' ')

    change1_cards.append('APERTURE')
    if object_name == '[B86] 133':
        change1_values.append(4./3600.)
    if object_name == 'HD 193256':
        change1_values.append(4.6/3600.)
    else:
        change1_values.append(6./3600.)
    change1_comments.append('[deg] Aperture width')

    for cards, values, comments in zip(change1_cards, change1_values, change1_comments):
        msgs.work('HDUL[1], Updating {}'.format(cards))
        fitsfiles.add_header_card(hdul[1].header, cards, values, comment=comments)

    # Now working on HDU0
    change0_cards, change0_values, change0_comments = [], [], []

    change0_cards.append('ORIGIN')
    change0_values.append('ESO-PARANAL')
    change0_comments.append('European Southern Observatory')

    change0_cards.append('TELESCOP')
    change0_values.append('ESO-VLT-U4')
    change0_comments.append('ESO telescope designation')

    change0_cards.append('INSTRUME')
    change0_values.append('MUSE')
    change0_comments.append('Instrument name')

    change0_cards.append('OBJECT')
    change0_values.append(object_name)
    change0_comments.append('Target designation')

    change0_cards.append('RA')
    change0_values.append(radec.ra.value)
    change0_comments.append('[deg] Spectroscopic target position (J2000.0)')

    change0_cards.append('DEC')
    change0_values.append(radec.dec.value)
    change0_comments.append('[deg] Spectroscopic target position (J2000.0)')

    change0_cards.append('RADECSYS')
    change0_values.append('FK5')
    change0_comments.append('Coordinate reference frame')

    change0_cards.append('EQUINOX')
    change0_values.append(2000.0)
    change0_comments.append('[y] Standard FK5')

    change0_cards.append('PRODCATG')
    change0_values.append('SCIENCE.SPECTRUM')
    change0_comments.append('Data product category')

    change0_cards.append('WAVELMIN')
    change0_values.append(np.nanmin(WAVE))
    change0_comments.append('[angstrom] Minimum wavelength')

    change0_cards.append('WAVELMAX')
    change0_values.append(np.nanmax(WAVE))
    change0_comments.append('[angstrom] Maximum wavelength')

    change0_cards.append('SPEC_RES')
    change0_values.append(2750.)
    change0_comments.append('Spectral resolution')

    change0_cards.append('SPEC_BIN')
    change0_values.append(hdr1['CDELT1'])
    change0_comments.append('[angstrom] Wavelength bin sizes')

    change0_cards.append('REFERENC')
    change0_values.append('doi:10.1051/0004-6361/201936178')
    change0_comments.append(' ')

    change0_cards.append('OBSTECH')
    change0_values.append('IFU')
    change0_comments.append('Technique of observation')

    change0_cards.append('EXT_OBJ')
    change0_values.append(False)
    change0_comments.append(' ')

    for cards, values, comments in zip(change0_cards, change0_values, change0_comments):
        msgs.work('HDUL[0], Updating {}'.format(cards))
        fitsfiles.add_header_card(hdul[0].header, cards, values, comment=comments)

    msgs.work('Updating checksum and datasum')
    hdul[0].add_datasum()
    hdul[1].add_datasum()
    hdul[0].add_checksum(override_datasum=True)
    hdul[1].add_checksum(override_datasum=True)
    msgs.work('Flushing changes.')
    hdul.flush()
    hdul.close()

    msgs.newline()
    msgs.info('File {} produced.'.format(output_fits))
    msgs.end()
