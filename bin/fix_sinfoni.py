#!/usr/bin/env python3

import argparse

import numpy as np
import os

from astropy.io import fits
from astropy import units as u

# from astropy import coordinates

from ESOAsg import __version__
from ESOAsg import msgs
from ESOAsg import images
# from ESOAsg.core import astro
from ESOAsg.core import fitsfiles
from ESOAsg.core import download_archive
from ESOAsg.ancillary import checks


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=r"""
        Manipulates SINFONI cubes to become compliant with the Phase 3 standard
        
        This code takes care of a series of simple rearrangements that transform SINFONI cube coming from the pipeline
        into (almost) Phase 3 compliant cubes (almost) ready to be ingested on the ESO Archive. The code trys its best
        to get all the keywords properly, however, the user should always double check the outcome to be sure that 
        everything is on place.
        
        To summarize the steps:
            1. Creates a PrimaryHDU and put the data in the `DATA` extension (similar to MUSE datacubes).
            2. Places properly header cards into the PrimaryHDU and in `DATA`
            3. Creates the white-light image 
            4. Performs few consistency checks for header values:
                * EXPTIME > 0 : If not, an attempt is made to recover the error from the PROV-j entries
        
        The following header keywords can be added by user. We refer to the `ESO Science Data Products Standard` for a
        exhaustive explanation of meanings and rules associated to them.
            * PrimaryHDU:
                * `FLUXCAL`  : either `ABSOLUTE` (`default`) or `UNCALIBRATED`
                * `REFERENC` : DOI to the related scientific publication
                * `ABMAGLIM` : 5-sigma detection limit derived from the white-light image
        
        Note that, if `output` is not defined, the fits file will be overwritten. 
        
        This uses ESOAsg version {__version__:s}
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument('input_fits', nargs='+', type=str,
                        help='Sinfoni fits file to be modified. This may contain wildcards.')
    parser.add_argument('-out', '--output', nargs='+', type=str, default=None,
                        help='Name of the output modified file')
    parser.add_argument('-cal', '--fluxcal', type=str, default='ABSOLUTE',
                        help='Value for the header keyword: `FLUXCAL`')
    parser.add_argument('-r', '--referenc', nargs='+', type=str, default=None,
                        help='DOI of the related scientific publication')
    parser.add_argument('-mlim', '--abmaglim', nargs='+', type=float, default=None,
                        help='5-sigma detection limit in the AB system')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    return parser.parse_args()


EXAMPLES = r"""
        Example:
        fix_sinfoni.py input_fits.fits 
        fix_sinfoni.py input_fits.fits -out _cleaned -cal ABSOLUTE
        """

if __name__ == '__main__':
    args = parse_arguments()

    # File_names
    input_fits_files = []
    for file_name in args.input_fits:
        if not checks.fits_file_is_valid(file_name):
            msgs.warning('File {} is not a valid fits file. Skipping the procedure'.format(file_name))
        elif len(fitsfiles.get_hdul(file_name)) > 1:
            msgs.warning('{} cube already processed? Skipping the procedure'.format(file_name))
        else:
            input_fits_files.append(file_name)
    if len(input_fits_files) == 0:
        msgs.error('No valid file to process')

    # output
    if args.output is not None:
        overwrite = False
        if len(input_fits_files) == 1:
            output_fits_files = [output_fits_temp if output_fits_temp.endswith('.fits') else np.str(
                output_fits_temp)+'.fits' for output_fits_temp in [args.output[0]]]
            output_fits_images = [output_fits_image.replace('.fits', '_whitelight.fits') for output_fits_image
                                  in output_fits_files]
        else:
            if np.str(args.output[0]).endswith('.fits'):
                output_ending = np.str(args.output[0])
            else:
                output_ending = np.str(args.output[0])+'.fits'
            output_fits_files = [output_temp.replace('.fits', output_ending) for output_temp in input_fits_files]
            output_fits_images = [output_fits_image.replace('.fits', '_whitelight'+output_ending) for output_fits_image
                                  in input_fits_files]

    else:
        output_fits_files = input_fits_files.copy()
        output_fits_images = [output_fits_image.replace('.fits', '_whitelight.fits') for output_fits_image in
                              input_fits_files]
        overwrite = True

    # reference
    if args.referenc is not None:
        reference = str(args.referenc[0])
    else:
        reference = str(' ')

    # fluxcal
    if args.fluxcal == 'ABSOLUTE':
        fluxcal = 'ABSOLUTE'
    elif args.fluxcal == 'UNCALIBRATED':
        fluxcal = 'UNCALIBRATED'
    else:
        msgs.error('Possible values for fluxcal are: `ABSOLUTE` or `UNCALIBRATED`')

    # abmaglim
    if args.abmaglim is not None:
        abmaglim = args.abmaglim
        assert isinstance(abmaglim, (int, np.float_)), 'ABMAGLIM must be a float'
        if abmaglim < 0:
            msgs.error('ABMAGLIM must be positive')
    else:
        abmaglim = np.float_(-1.)

    msgs.start()

    for fits_in, fits_out, image_out in zip(input_fits_files, output_fits_files, output_fits_images):
        if os.path.exists(fits_out):
            os.rename(fits_out, fits_out.replace('.fits', '_old.fits'))
            msgs.warning('{} already exists. Backup created.'.format(fits_out))
        if os.path.exists(image_out):
            os.rename(image_out, image_out.replace('.fits', '_old.fits'))
            msgs.warning('{} already exists. Backup created.'.format(image_out))

        # 1. create a copy of the file where there is a primary HDU and data are in the 'DATA" HDU
        fitsfiles.new_fits_like(fits_in, [0], fits_out, overwrite=overwrite)
        hdul = fitsfiles.get_hdul(fits_out, 'update', checksum=True)
        hdr0 = hdul[0].header
        hdr1 = hdul[1].header

        # 1.1 Check for HISTORY

        if 'HISTORY' in hdr0.keys():
            history_cards_hdr0 = [history_card_hdr0 for history_card_hdr0 in hdr0
                                  if history_card_hdr0.startswith('HISTORY')]
            history_values_hdr0 = [hdr0[history_card_hdr0] for history_card_hdr0 in hdr0
                                   if history_card_hdr0.startswith('HISTORY')]
            del hdr0['HISTORY'][:]

        if 'HISTORY' in hdr1.keys():
            history_values_hdr1 = hdr1['HISTORY'][:]
            for history_number in range(0, len(history_values_hdr1)):
                clean_history = checks.remove_non_ascii(history_values_hdr1[history_number])
                if len(clean_history) > 0:
                    hdr1['HISTORY'][history_number] = str(clean_history)
                else:
                    hdr1['HISTORY'][history_number] = str(' ')

        # 2. update cards for headers:

        # Updating values with different CARD in the header
        CARDS_INPUT = ['ESO OBS ID', 'ESO OBS PROG ID', 'ESO PRO ANCESTOR', 'ESO PRO DATANCOM',
                       'ESO PRO TECH', 'ESO PRO REC1 PIPE ID']
        CARDS_OUTPUT = ['OBID1', 'PROG_ID', 'PROV1', 'NCOMBINE',
                        'OBSTECH', 'PROCSOFT']
        fitsfiles.transfer_header_cards(hdr1, hdr0, CARDS_INPUT, output_cards=CARDS_OUTPUT, delete_card=True)

        # Transfer cards from HDU1 to the PrimaryHDU
        not_to_be_transfer = [hdr1_card for hdr1_card in hdr1 if
                              hdr1_card.startswith('COMMENT') or
                              hdr1_card.startswith('NAXIS') or
                              hdr1_card.startswith('CRPIX') or
                              hdr1_card.startswith('CRVAL') or
                              hdr1_card.startswith('CDELT') or
                              hdr1_card.startswith('CTYPE') or
                              hdr1_card.startswith('CD1_') or
                              hdr1_card.startswith('CD2_') or
                              hdr1_card.startswith('CD3_') or
                              hdr1_card.startswith('CUNIT') or
                              hdr1_card.startswith('CSYER') or
                              hdr1_card.startswith('HDUCLAS') or
                              hdr1_card.startswith('XTENSION') or
                              hdr1_card.startswith('PCOUNT') or
                              hdr1_card.startswith('GCOUNT') or
                              hdr1_card.startswith('HDUDOC') or
                              hdr1_card.startswith('HDUVER') or
                              hdr1_card.startswith('HISTORY')]
        cards_to_be_transfer = [hdr1_card for hdr1_card in hdr1 if hdr1_card not in not_to_be_transfer]
        fitsfiles.transfer_header_cards(hdr1, hdr0, cards_to_be_transfer, with_comment=True, delete_card=True)

        # Updating the FITS file definition comment line
        hdr0['COMMENT'] = "  FITS (Flexible Image Transport System) format is defined in 'Astronomy" \
                          + "  and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H"

        # Including specific cards
        msgs.work('Updating SPEC_RES')
        if hdr0['HIERARCH ESO INS GRAT1 NAME'].strip() == 'J':
            hdr0['SPEC_RES'] = 2000.
        elif hdr0['HIERARCH ESO INS GRAT1 NAME'].strip() == 'H':
            hdr0['SPEC_RES'] = 3000.
        elif hdr0['HIERARCH ESO INS GRAT1 NAME'].strip() == 'K':
            hdr0['SPEC_RES'] = 4000.
        elif hdr0['HIERARCH ESO INS GRAT1 NAME'].strip() == 'H + K' or hdr0['HIERARCH ESO INS GRAT1 NAME'].strip() == \
                'H+K':
            hdr0['SPEC_RES'] = 1500.
        else:
            msgs.error('GRAT1 name: {} not recognized'.format(hdr0['HIERARCH ESO INS GRAT1 NAME']))

        msgs.work('Updating PRODCATG')
        hdr0['PRODCATG'] = str('SCIENCE.CUBE.IFS')

        if abmaglim > 0.:
            msgs.work('Updating ABMAGLIM')
            hdr0['ABMAGLIM'] = abmaglim

        msgs.work('Updating REFERENC')
        hdr0['REFERENC'] = reference

        msgs.work('Updating PROVENANCE')
        prov_to_be_transfer = [hdr0_card for hdr0_card in hdr0 if hdr0_card.startswith('ESO PRO REC')
                               and 'RAW' in hdr0_card and 'NAME' in hdr0_card and hdr0[hdr0_card].startswith('SINFO.')
                               and hdr0[hdr0_card.replace('NAME', 'CATG')].startswith('OBJECT')]

        for prov, prov_number in zip(prov_to_be_transfer, range(1, len(prov_to_be_transfer), 1)):
            hdr0['PROV'+str(prov_number)] = hdr0[prov]
            msgs.work('    - {}'.format(hdr0[prov]))
            msgs.work('      {}'.format(hdr0[prov.replace('NAME', 'CATG')]))

        msgs.work('Updating MJD-END')
        EXPTIME_sec = np.float32(hdr0['EXPTIME'])
        EXPTIME_day = EXPTIME_sec / (60. * 60. * 24.)
        MJDEND = np.float32(hdr0['MJD-OBS']) + EXPTIME_day
        fitsfiles.add_header_card(hdr0, 'MJD-END', MJDEND, 'End of observation')

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

        msgs.work('Updating FLUXCAL')
        fitsfiles.add_header_card(hdr0, 'FLUXCAL', fluxcal)

        # update cards for `DATA` header
        fitsfiles.add_header_card(hdr1, 'EXTNAME', 'DATA', 'This extension contains data value')
        fitsfiles.add_header_card(hdr1, 'XTENSION', 'IMAGE', 'IMAGE extension')

        #  Dealing with douplicate FITS comments
        all_comments = [comment for comment in hdr1['COMMENT']]
        del hdr1['COMMENT']
        for comments_to_be_kept in all_comments:
            if '(Flexible Image Transport System)' not in comments_to_be_kept \
                    and '2001A&A...376..359H' not in comments_to_be_kept:
                hdr1['COMMENT'] = str('  ' + comments_to_be_kept)

        # 3. create white light image

        # Updating ASSON value in the header
        msgs.work('Updating ASSON')
        hdr0['ASSON1'] = str(image_out)

        # Actually creating the white-light image
        image_hdu = fits.PrimaryHDU()
        image_hdul = fits.HDUList([image_hdu])
        if hdr1['CUNIT3'].startswith('um'):
            to_ang = 1000.
        else:
            msgs.error('Spectral unit: {} not recognized'.format(hdr1['CUNIT3']))
        delta_wave_bin = hdr1['CDELT3']
        image_hdul.append(fits.ImageHDU(to_ang*delta_wave_bin*np.nansum(hdul[1].data, axis=0, dtype=np.float_)))
        image_hdr0 = image_hdul[0].header
        image_hdr1 = image_hdul[1].header
        card_for_image0 = ['WAVELMIN', 'WAVELMAX', 'OBJECT', 'TELESCOP', 'INSTRUME', 'RADECSYS', 'RA', 'DEC',
                           'EQUINOX']
        fitsfiles.transfer_header_cards(hdr0, image_hdr0, card_for_image0, with_comment=True, delete_card=False)
        image_hdr0['PRODCATG'] = str('ANCILLARY.IMAGE.WHITELIGHT')

        card_for_image1 = ['CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2', 'CDELT1', 'CDELT2', 'CUNIT1', 'CUNIT2',
                            'NAXIS1', 'NAXIS2', 'EXTNAME', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'CTYPE1', 'CTYPE2']
        fitsfiles.transfer_header_cards(hdr1, image_hdr1, card_for_image1, with_comment=True, delete_card=False)

        # 5. few checks

        # EXPTIME
        if hdr0['EXPTIME'] <= 0:
            get_prov = [hdr0[prov] for prov in hdr0 if prov.startswith('PROV')]
            if len(get_prov) > 0:
                msgs.warning('Inconsiste value for EXPTIME: {}. Tying to fix it with provenances'.format(hdr0[
                                                                                                            'EXPTIME']))
                msgs.work('Downloading {} headers'.format(len(get_prov)))
                EXPTIME_sec = 0
                for prov_number in hdr0['PROV*']:
                    file_id = hdr0[prov_number].replace('.fits', '')
                    download_archive.get_header_from_archive(file_id, text_file=file_id+'.hdr')
                    hdr_prov = fitsfiles.header_from_txt_file(file_id+'.hdr')
                    msgs.work('{} has EXPTIME of: {}'.format(file_id, hdr_prov['EXPTIME']))
                    EXPTIME_sec = EXPTIME_sec+hdr_prov['EXPTIME']
                hdr0['EXPTIME'] = np.float(EXPTIME_sec)
                msgs.warning('Updating value for EXPTIME to: {}'.format(EXPTIME_sec))
                fitsfiles.add_header_card(hdr0, 'TEXPTIME', EXPTIME_sec)
                EXPTIME_day = np.float(EXPTIME_sec) / (60. * 60. * 24.)
                MJDEND = np.float32(hdr0['MJD-OBS']) + EXPTIME_day
                fitsfiles.add_header_card(hdr0, 'MJD-END', MJDEND, 'End of observation')
                msgs.warning('Updating value for MJD-END to: {}'.format(MJDEND))
            else:
                msgs.error('Inconsiste value for EXPTIME: {}'.format(hdr0['EXPTIME']))

        # ABMAGLIM
        if hdr0['FLUXCAL'] is 'ABSOLUTE':
            if 'ABMAGLIM' not in hdr0.keys():
                msgs.warning('Missing ABMAGLIM keyword')
                msgs.warning('Trying to estimate it from whitelight image')
                whitelight = images.Images(data=image_hdul[1].data)
                whitelight.calc_background(method='median', nsigma=3., find_sources=True, src_nsigma=3.,
                                           src_npixels=5, src_dilate_size=2)
                whitelight_mean, whitelight_median, whitelight_std = whitelight.get_clean_stats(nsigma=3.)
                # ToDO
                # Improve on this
                n_pixels_psf = np.pi*(3)**2.
                five_sigma_nu = 5.*3.34e4*np.power((WAVELMAX+WAVELMIN)*to_ang/2.,
                                                   2.)*whitelight_std*n_pixels_psf/(delta_wave_bin*to_ang)
                hdr0['ABMAGLIM'] = -2.5*np.log10(five_sigma_nu/3631.)
                msgs.warning('ABMAGLIM={}. This is most probably not correct at the moment'.format(hdr0['ABMAGLIM'] ))

        # 6. update checksum and datasum
        msgs.work('Updating checksum and datasum')
        hdul[0].add_datasum()
        hdul[1].add_datasum()
        hdul[0].add_checksum(override_datasum=True)
        hdul[1].add_checksum(override_datasum=True)
        image_hdul[0].add_datasum()
        image_hdul[1].add_datasum()
        image_hdul[0].add_checksum(override_datasum=True)
        image_hdul[1].add_checksum(override_datasum=True)

        """
        # CARDS are already present in hdr0
        # HDUCLASS', 'HDUCLAS1', 'HDUCLAS2', 'HDUDOC', 'HDUVERS'
        """

        image_hdul.writeto(image_out)
        hdul.flush()
        hdul.close()
        msgs.info('File {} produced.'.format(fits_out))
        msgs.info('Image {} produced.'.format(image_out))

    msgs.end()
