#!/usr/bin/env python3

import argparse

import numpy as np
import os

from astropy.io import fits

# from astropy import coordinates
# from astropy import units as u

from ESOAsg import __version__
from ESOAsg import msgs
from ESOAsg.core import fitsfiles
from ESOAsg.ancillary import checks

# from IPython import embed


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=r"""
        Manipulates SINFONI cubes to become compliant with the Phase 3 standard
        
        This uses ESOAsg version {:s}
        """.format(__version__),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument('input_fits', nargs='+', type=str,
                        help='Sinfoni fits file to be modified. This may contain wildcards.')
    parser.add_argument('-out', '--output', nargs='+', type=str, default=None,
                        help='Name of the output modified file')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    return parser.parse_args()

EXAMPLES = r"""
        Example:
        fix_sinfoni.py input_fits.fits 
        """

if __name__ == '__main__':
    args = parse_arguments()

    # File_names
    for file_name in args.input_fits:
        if not checks.fits_file_is_valid(file_name):
            msgs.error('File {} is not a valid fits file'.format(file_name))
    input_fits_files = args.input_fits

    # output
    if args.output is not None:
        overwrite = False
        if len(input_fits_files) == 1:
            output_fits_files = args.output
        else:
            if np.str(args.output[0]).endswith('.fits'):
                output_ending = np.str(args.output[0])
            else:
                output_ending = np.str(args.output[0])+'.fits'
            output_fits_files = [output_temp.replace('.fits', output_ending) for output_temp in input_fits_files]
    else:
        output_fits_files = input_fits_files
        overwrite = True

    msgs.start()

    for fits_in, fits_out in zip(input_fits_files, output_fits_files):
        if os.path.exists(fits_out):
            os.rename(fits_out, fits_out.replace('.fits', '_old.fits'))
            msgs.warning('{} already exists. Backup created.'.format(fits_out))

        # 1. create a copy of the file where there is a primary HDU and data are in the 'DATA" HDU
        fitsfiles.new_fits_like(fits_in, [0], fits_out, overwrite=overwrite)
        hdul = fitsfiles.get_hdul(fits_out, 'update', checksum=True)
        hdr0 = hdul[0].header
        hdr1 = hdul[1].header

        # 2. update cards for primary header:
        not_to_be_transfer = [hdr1_card for hdr1_card in hdr1 if
                              hdr1_card.startswith('COMMENT') or
                              hdr1_card.startswith('NAXIS') or
                              hdr1_card.startswith('CRPIX') or
                              hdr1_card.startswith('CRVAL') or
                              hdr1_card.startswith('CTYPE') or
                              hdr1_card.startswith('CD1_') or
                              hdr1_card.startswith('CD2_') or
                              hdr1_card.startswith('CD3_') or
                              hdr1_card.startswith('CUNIT') or
                              hdr1_card.startswith('CSYER') or
                              hdr1_card.startswith('HDUCLAS')]
        cards_to_be_transfer = [hdr1_card for hdr1_card in hdr1 if hdr1_card not in not_to_be_transfer]
        fitsfiles.transfer_header_cards(hdr1, hdr0, cards_to_be_transfer, with_comment=True, delete_card=True)

        EXPTIME_sec = np.float32(hdr0['EXPTIME'])
        EXPTIME_day = EXPTIME_sec / (60. * 60. * 24.)
        MJDEND = np.float32(hdr0['MJD-OBS']) + EXPTIME_day
        fitsfiles.add_header_card(hdr0, 'MJD-END', MJDEND, 'End of observation')

        # 3. update cards for `DATA` header
        fitsfiles.add_header_card(hdr1, 'EXTNAME', 'DATA', 'This extension contains data value')
        fitsfiles.add_header_card(hdr1, 'XTENSION', 'IMAGE', 'IMAGE extension')

        # 4. create white light image

        # 5. update checksum and datasum
        msgs.work('Updating checksum and datasum')
        hdul[0].add_datasum()
        hdul[1].add_datasum()
        hdul[0].add_checksum(override_datasum=True)
        hdul[1].add_checksum(override_datasum=True)

        """
        msgs.work('Updating MJD-END')


        msgs.work('Updating TEXPTIME')
        fitsfiles.add_header_card(hdr0, 'TEXPTIME', EXPTIME_sec)

        msgs.work('Updating WAVELMIN')
        WAVELMIN = np.float32(hdr1['CDELT3']) * (1. - np.float32(hdr1['CRPIX3'])) + np.float32(hdr1['CRVAL3'])
        fitsfiles.add_header_card(hdr0, 'WAVELMIN', WAVELMIN)

        msgs.work('Updating WAVELMAX')
        zMax, _, _ = np.shape(hdul[1].data)
        WAVELMAX = np.float32(hdr1['CDELT3']) * (np.float32(zMax) - np.float32(hdr1['CRPIX3'])) + np.float32(
            hdr1['CRVAL3'])
        fitsfiles.add_header_card(hdr0, 'WAVELMAX', WAVELMAX)

        msgs.work('Updating REFERENC')
        fitsfiles.add_header_card(hdr0, 'REFERENC', ' ')

        # Updating values with the same CARD in the header
        CARDS_LIST = ['DATE', 'ORIGIN', 'TELESCOP', 'INSTRUME', 'RA', 'DEC',
                      'RADECSYS', 'EQUINOX', 'DATE-OBS', 'MJD-OBS',
                      'SPECSYS', 'PRODCATG']
        fitsfiles.transfer_header_cards(hdr1, hdr0, CARDS_LIST, delete_card=True)

        # Updating values without removing the original ones
        CARDS_LIST = ['OBJECT', 'EXPTIME']
        fitsfiles.transfer_header_cards(hdr1, hdr0, CARDS_LIST, delete_card=False)

        # Updating values with different CARD in the header
        CARDS_INPUT = ['ESO OBS ID', 'ESO OBS PROG ID', 'ESO PRO ANCESTOR', 'ESO PRO DATANCOM',
                       'ESO PRO TECH', 'ESO PRO REC1 PIPE ID']
        CARDS_OUTPUT = ['OBID1', 'PROG_ID', 'PROV1', 'NCOMBINE',
                        'OBSTECH', 'PROCSOFT']
        fitsfiles.transfer_header_cards(hdr1, hdr0, CARDS_INPUT, output_cards=CARDS_OUTPUT, delete_card=True)

        # CARDS are already present in hdr0
        # HDUCLASS', 'HDUCLAS1', 'HDUCLAS2', 'HDUDOC', 'HDUVERS'
        """
        hdul.flush()
        hdul.close()
        msgs.info('File {} produced.'.format(fits_out))

    msgs.end()
