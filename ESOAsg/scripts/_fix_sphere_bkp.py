r"""
fix_sphere
==========
Script to make data produced by the `SPHERE Data Center <https://sphere.osug.fr>`_ (almost) phase 3 compliant

.. topic:: Inputs:

    - **input_fits** - Input fits file

"""


import argparse

from ESOAsg import __version__

SUPPORTED_INSTRUMENT = ['IFS', 'IRDIS']

EXAMPLES = str(r"""EXAMPLES:""" + """\n""" + """\n""" +
               r"""Given a IFS file in input from the `SPHERE Data Center` transform it in a """ +
               r"""fits file with the same name and suffix `_phase3` that is """ +
               r"""(almost) compliant with the """ +
               r"""format of a PRODCATG = 'SCIENCE.CUBE.IFS' and create the """ +
               r"""corresponding whitelight image:""" + """\n""" +
               r""">>> fix_sphere  BD11_IFS_CUBE_cADI.fits""" +
               r"""--suffix _phase3 """ +
               r"""--whitelight """ + """\n""" +
               r""" """)


def parser(options=None):
    parser = argparse.ArgumentParser(
        description=r"""Manipulate SPHERE data produced by the `SPHERE Data Center` (https://sphere.osug.fr) """ +
                    r"""and make them compliant with the  """ +
                    r"""`ESO Phase 3 Standard` (https://www.eso.org/sci/observing/phase3/p3sdpstd.pdf). """ +
                    """\n""" + """\n""" +
                    r"""This code takes care of a series of rearrangements that transform SPHERE """ +
                    r"""data into (almost) Phase 3 compliant data (almost) ready to be ingested on the ESO """ +
                    r"""Archive. The code tries its best to get all the keywords properly. """ +
                    r"""However, it is user's responsibility to check the outcome """ +
                    r"""and to make sure that everything is on place. """ +
                    """\n""" + """\n""" +
                    r"""To summarize the steps for IFS data: """ +
                    """\n""" +
                    r"""- Creates a PrimaryHDU and put the data in the `DATA` extension; """ +
                    """\n""" +
                    r"""- Places properly header cards into the PrimaryHDU and in `DATA` """ +
                    """\n""" +
                    r"""- Fix astrometry based on object name """ +
                    """\n""" +
                    r"""- Try to reconstruct some missing header keywords """ +
                    """\n""" +
                    r"""- Make files fits compliant """ +
                    """\n""" +
                    r"""- Creates the white-light images """ +
                    """\n""" + """\n""" +
                    r"""To summarize the steps for IRDIS data: """ +
                    """\n""" +
                    r"""- Creates a PrimaryHDU and put the data in the `DATA` extension; """ +
                    """\n""" +
                    r"""- Places properly header cards into the PrimaryHDU and in `DATA` """ +
                    """\n""" +
                    r"""- Fix astrometry based on object name """ +
                    """\n""" +
                    r"""- Try to reconstruct some missing header keywords """ +
                    """\n""" +
                    r"""- Make files fits compliant """ +
                    """\n""" + """\n""" +
                    r"""The following header keywords can be added by user. """ +
                    r"""We refer to the `ESO Science Data Products Standard` for an exhaustive explanation """ +
                    r"""of meanings and rules associated to them. """ +
                    """\n""" +
                    r"""PrimaryHDU: """ +
                    """\n""" +
                    r"""* `FLUXCAL`  : either `ABSOLUTE` (`default`) or `UNCALIBRATED` """ +
                    """\n""" +
                    r"""* `REFERENC` : DOI to the related scientific publication """ +
                    """\n""" +
                    r"""* `ABMAGLIM` : 5-sigma detection limit derived from the white-light image """ +
                    """\n""" + """\n""" +
                    r"""Note that, if `suffix` is not defined, the fits file will be overwritten.  """ +
                    """\n""" + """\n""" +
                    r"""This uses ESOAsg version {:s}""".format(__version__),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument("input_fits", nargs="+", type=str,
                        help=r"Input SPHERE fits file")
    parser.add_argument("-s", "--suffix", nargs="+", type=str, default=None,
                        help=r"Suffix to be added to file names the will be used as output." +
                             r"If it is not set, the input files will be overwritten")
    parser.add_argument("-wl", "--whitelight", action="store_true", default=False,
                        help=r"Create the white light image for IFS cubes")
    parser.add_argument("-v", "--version", action="version", version=__version__)
    return parser.parse_args()


def main(args):
    import numpy as np
    import os
    import shutil
    from astropy.coordinates import SkyCoord
    from astropy.coordinates import name_resolve
    from astropy import units as u
    from astropy.io import fits
    from ESOAsg.ancillary import cleaning_lists
    from ESOAsg.ancillary import cleaning_headers
    from ESOAsg.core import fitsfiles
    from ESOAsg import msgs

    # Cleaning input lists
    input_fits_files = cleaning_lists.make_list_of_fits_files(args.input_fits)
    # Make whitelight images
    if args.whitelight:
        make_whitelight_image = True
    else:
        make_whitelight_image = False
    # Creating output list
    if args.suffix is None:
        overwrite = False
        msgs.warning('The file will overwrite the input files')
    else:
        overwrite = True
    suffix_string = cleaning_lists.make_string(args.suffix)
    output_fits_files = cleaning_lists.make_list_of_fits_files_and_append_suffix(input_fits_files,
                                                                                 suffix=suffix_string)
    if make_whitelight_image:
        output_whitelight_files = cleaning_lists.make_list_of_fits_files_and_append_suffix(input_fits_files,
                                                                                           suffix=suffix_string + '_WL')
    else:
        output_whitelight_files = [None] * len(input_fits_files)
    '''
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
    '''

    msgs.start()

    for fits_in, fits_out, image_out in zip(input_fits_files, output_fits_files, output_whitelight_files):
        if os.path.exists(fits_out):
            shutil.copy(fits_out, fits_out.replace('.fit', '_old.fit'))
            msgs.warning('{} already exists. Backup created.'.format(fits_out))
        if image_out is not None:
            if os.path.exists(image_out):
                shutil.copy(image_out, image_out.replace('.fit', '_old.fit'))
                msgs.warning('{} already exists. Backup created.'.format(image_out))

        full_hdul = fitsfiles.get_hdul(fits_in)
        try:
            instrument = full_hdul[0].header['HIERARCH ESO SEQ ARM']
        except KeyError:
            msgs.error('Failed to read the keyword HIERARCH ESO SEQ ARM from the primary header')
        finally:
            if instrument in SUPPORTED_INSTRUMENT:
                msgs.info('The input file is from SPHERE/{}'.format(instrument))
            else:
                msgs.warning('Instrument SPHERE/{} not supported'.format(instrument))

        # ToDo
        # These needs to be transformed in objects
        if instrument.startswith('IFS'):
            msgs.work('Fixing header for SPHERE/{} file {}'.format(instrument, fits_in))

            # Create a copy of the file where there is a primary HDU and data are in the 'DATA" HDU
            msgs.work('Reshaping cube into PrimaryHEADER and Data Header')
            fitsfiles.new_fits_like(fits_in, [0], fits_out, overwrite=overwrite, fix_header=True)
            hdul = fitsfiles.get_hdul(fits_out, 'update', checksum=True)
            hdr0 = hdul[0].header
            hdr1 = hdul[1].header

            # Clean HISTORY keywords
            _ = cleaning_headers.history_keywords_delete(hdr0, verbose=True, in_place=True)
            _ = cleaning_headers.history_keywords_delete(hdr1, verbose=True, in_place=True)

            # Update cards for headers:
            # Updating values with different CARD in the header
            cards_input = ['CRPIX4', 'CRVAL4', 'CTYPE4', 'CUNIT4', 'CD4_4', 'CD1_4', 'CD2_4', 'CD4_1', 'CD4_2']
            cards_output = ['CRPIX3', 'CRVAL3', 'CTYPE3', 'CUNIT3', 'CD3_3', 'CD1_3', 'CD2_3', 'CD3_1', 'CD3_2']
            fitsfiles.transfer_header_cards(hdr1, hdr1, cards_input, output_cards=cards_output, delete_card=True)
            # Remove not used values
            cards_to_be_removed_hdr1 = ['CD4_3', 'CD3_4']
            for card_to_be_removed_hdr1 in cards_to_be_removed_hdr1:
                hdr1.remove(card_to_be_removed_hdr1, ignore_missing=True)

            # Transfer cards from HDU1 to the PrimaryHDU
            not_to_be_transfer = [hdr1_card for hdr1_card in hdr1 if
                                  hdr1_card.startswith('COMMENT') or
                                  hdr1_card.startswith('EXTNAME') or
                                  hdr1_card.startswith('BITPIX') or
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
                                  hdr1_card.startswith('HDUVER')]
            cards_to_be_transfer = []
            for hdr1_card in hdr1:
                if hdr1_card not in not_to_be_transfer:
                    cards_to_be_transfer.append(hdr1_card)
            fitsfiles.transfer_header_cards(hdr1, hdr0, cards_to_be_transfer, with_comment=True, delete_card=True)

            # Try to guess coordinates
            if 'CRVAL1' not in hdr1.keys():
                msgs.warning('CRVAL position keywords not preset')
                if 'OBJECT' in hdr0.keys():
                    try:
                        object_coordinate = SkyCoord.from_name(str(hdr0['OBJECT']).strip())
                        ra_obj, dec_obj = object_coordinate.ra.degree, object_coordinate.dec.degree
                        if 'RA' in hdr0.keys() and 'DEC' in hdr0.keys():
                            pointing_coordinate = SkyCoord(float(hdr0['RA']), float(hdr0['DEC']), unit='deg')
                            msgs.work('Testing from separation from pointing position')
                            separation = object_coordinate.separation(pointing_coordinate).arcsec
                            if separation < 120.:
                                msgs.info('Object - Pointing separation is {}'.format(separation))
                                msgs.info('Updating CRVAL1 = {}'.format(ra_obj))
                                msgs.info('Updating CRVAL2 = {}'.format(dec_obj))
                                hdr1['CRVAL1'] = ra_obj
                                hdr1['CRVAL2'] = dec_obj
                                msgs.work('Updating CUNIT')
                                hdr1['CUNIT1'] = 'deg'
                                hdr1['CUNIT2'] = 'deg'
                                msgs.work('Updating CTYPE')
                                hdr1['CTYPE1'] = 'RA---TAN'
                                hdr1['CTYPE2'] = 'DEC--TAN'
                                msgs.work('Updating CRPIX')
                                hdr1['CRPIX1'] = float(hdul[1].data.shape[2])/2.
                                hdr1['CRPIX2'] = float(hdul[1].data.shape[1])/2.
                                msgs.info('Updating CD1 and CD2')
                                hdr1['CD1_1'] = 2.06E-06
                                hdr1['CD2_2'] = 2.06E-06
                                hdr1['CD1_2'] = 0.
                                hdr1['CD2_1'] = 0.
                                msgs.work('Updating RA, DEC')
                                hdr0['RA'] = ra_obj
                                hdr0.comments['RA'] = object_coordinate.ra.to_string(u.hour)
                                hdr0['DEC'] = dec_obj
                                hdr0.comments['DEC'] = object_coordinate.dec.to_string(u.degree, alwayssign=True)
                            else:
                                msgs.warning('Object - Pointing separation is {}'.format(separation))
                                msgs.warning('This is suspicious, CRVAL not updated')
                    except name_resolve.NameResolveError:
                        msgs.warning('Object {} not recognized'.format(str(hdr0['OBJECT']).strip()))
                        msgs.warning('CRVAL not updated')

            # Updating file prodcatg
            msgs.work('Updating PRODCATG to SCIENCE.CUBE.IFS')
            hdr0['PRODCATG'] = str('SCIENCE.CUBE.IFS')
            # Some more updates
            msgs.work('Setting NAXIS = 0 in primary header')
            hdr0['NAXIS'] = 0
            if 'OBSTECH' not in hdr0.keys():
                msgs.warning('OBSTECH missing')
                if 'ESO PRO TECH' in hdr0.keys():
                    msgs.info('Deriving OBSTECH from HIERARCH ESO PRO TECH')
                    msgs.work('Updating OBSTECH to {}'.format(str(hdr0['HIERARCH ESO PRO TECH'])))
                    hdr0['OBSTECH'] = str(hdr0['HIERARCH ESO PRO TECH'])
            if 'EXPTIME' not in hdr0.keys():
                msgs.warning('EXPTIME missing')
                if 'ESO DET SEQ1 REALDIT' in hdr0.keys() and 'ESO DET NDIT' in hdr0.keys():
                    msgs.info('Deriving EXPTIME and TEXPTIME as REALDIT * DIT')
                    hdr0['EXPTIME'] = hdr0['HIERARCH ESO DET SEQ1 REALDIT'] * hdr0['HIERARCH ESO DET NDIT']
                    hdr0['TEXPTIME'] = hdr0['HIERARCH ESO DET SEQ1 REALDIT'] * hdr0['HIERARCH ESO DET NDIT']
                    msgs.work('Updating EXPTIME to {}'.format(str(hdr0['EXPTIME'])))
                    msgs.work('Updating TEXPTIME to {}'.format(str(hdr0['TEXPTIME'])))
            if 'WAVELMIN' not in hdr0.keys():
                msgs.warning('WAVELMIN missing')
                z_pixel = np.arange(int(hdul[1].data.shape[0]))
                z_wave = float(hdr1['CRVAL3']) + (z_pixel * float(hdr1['CD3_3']))
                if str(hdr1['CUNIT3']).strip().upper() == 'MICRONS':
                    msgs.info('Deriving WAVELMIN and WAVELMAX from CRVAL1')
                    z_wave = z_wave * 1000. # convert to nanometers
                    hdr0['WAVELMIN'] = np.nanmin(z_wave)
                    hdr0['WAVELMAX'] = np.nanmax(z_wave)
                    msgs.work('Updating WAVELMIN to {}'.format(str(hdr0['WAVELMIN'])))
                    msgs.work('Updating WAVELMAX to {}'.format(str(hdr0['WAVELMAX'])))
                else:
                    msgs.warning('Unknown units {}. WAVELMIN and WAVELMAX not calculated'.format(str(hdr1['CUNIT3'])))
            if 'SPEC_RES' not in hdr0.keys():
                msgs.warning('SPEC_RES missing')
                if 'WAVELMAX' in hdr0.keys():
                    msgs.info('Deriving SPEC_RES from WAVELMAX')
                    if (float(hdr0['WAVELMAX']) > 1300.) and (float(hdr0['WAVELMAX']) < 1400.):
                        hdr0['SPEC_RES'] = 50.
                        msgs.work('Updating SPEC_RES to {}'.format(str(hdr0['SPEC_RES'])))
                    elif (float(hdr0['WAVELMAX']) > 1600.) and (float(hdr0['WAVELMAX']) < 1700.):
                        hdr0['SPEC_RES'] = 30.
                        msgs.work('Updating SPEC_RES to {}'.format(str(hdr0['SPEC_RES'])))
                    else:
                        msgs.warning('WAVELMAX = {} is not in the expected ' /
                                     + 'range of possible values'.format(str(hdr0['WAVELMAX'])))
            if 'PROGID' not in hdr0.keys():
                msgs.warning('PROG_ID missing')
                if 'ESO OBS PROG ID' in hdr0.keys():
                    msgs.info('Deriving PROG_ID from HIERARCH ESO OBS PROG ID')
                    msgs.work('Updating PROG_ID to {}'.format(str(hdr0['HIERARCH ESO OBS PROG ID'])))
                    hdr0['PROG_ID'] = str(hdr0['HIERARCH ESO OBS PROG ID'])
            if 'MJD-END' not in hdr0.keys():
                msgs.warning('MJD-END missing')
                if 'TEXPTIME' in hdr0.keys():
                    msgs.info('Deriving MJD-END from MJD-OBS and TEXPTIME')
                    texptime_sec = float(hdr0['TEXPTIME'])
                    texptime_day = texptime_sec / (60. * 60. * 24.)
                    mjdend = float(hdr0['MJD-OBS']) + texptime_day
                    fitsfiles.add_header_card(hdr0, 'MJD-END', mjdend, 'End of observation')
                    msgs.work('MJD-OBS = {} and TEXPTIME = {} days'.format(str(hdr0['MJD-OBS']), str(texptime_day)))
                    msgs.work('Updating MJD-END to {}'.format(str(hdr0['MJD-END'])))

            # Remove not used values
            cards_to_be_removed_hdr0 = ['ERRDATA', 'QUALDATA', 'SCIDATA']
            for card_to_be_removed_hdr0 in cards_to_be_removed_hdr0:
                hdr0.remove(card_to_be_removed_hdr0, ignore_missing=True)
            cards_to_be_removed_hdr1 = ['HDUCLASS3']
            for card_to_be_removed_hdr1 in cards_to_be_removed_hdr1:
                hdr1.remove(card_to_be_removed_hdr1, ignore_missing=True)

            # Updating the FITS file definition comment line
            hdr0.add_comment("  FITS (Flexible Image Transport System) format is defined in 'Astronomy" + "  and "
                             "Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H", after='EXTEND')
            if 'COMMENT' in hdr1.keys():
                comment_values_hdr1 = hdr1['COMMENT'][:]
                for index, comment_value_hdr1 in enumerate(comment_values_hdr1):
                    msgs.work('Removing COMMENT card : {}'.format(comment_value_hdr1))
                hdr1.remove('COMMENT', ignore_missing=True, remove_all=True)

            # Creating white light image keyword:
            if make_whitelight_image:
                fitsfiles.add_header_card(hdr0, 'ASSON1', image_out.split('/')[-1],
                                          'ANCILLARY.IMAGE.WHITELIGHT filename')
                msgs.work('Updating ASSON1 to {}'.format(hdr0['ASSON1']))

            # Actually creating the white-light image
            if make_whitelight_image:
                msgs.info('Making white light image')
                image_hdu = fits.PrimaryHDU()
                image_hdul = fits.HDUList([image_hdu])
                if str(hdr1['CUNIT3']).strip().upper() == 'MICRONS':
                    to_ang = 10000.
                else:
                    msgs.error('Spectral unit: {} not recognized'.format(hdr1['CUNIT3']))
                delta_wave_bin = hdr1['CD3_3']
                image_hdul.append(fits.ImageHDU(to_ang * delta_wave_bin * np.nansum(hdul[1].data,
                                                                                    axis=0, dtype=np.float_)))
                image_hdr0 = image_hdul[0].header
                image_hdr1 = image_hdul[1].header
                card_for_image0 = ['WAVELMIN', 'WAVELMAX', 'OBJECT', 'TELESCOP', 'INSTRUME', 'RADECSYS', 'RA', 'DEC',
                                   'EQUINOX']
                fitsfiles.transfer_header_cards(hdr0, image_hdr0, card_for_image0, with_comment=True, delete_card=False)
                image_hdr0['PRODCATG'] = str('ANCILLARY.IMAGE.WHITELIGHT')

                card_for_image1 = ['CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2', 'CUNIT1', 'CUNIT2',
                                   'NAXIS1', 'NAXIS2', 'EXTNAME', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'CTYPE1', 'CTYPE2']
                fitsfiles.transfer_header_cards(hdr1, image_hdr1, card_for_image1, with_comment=True, delete_card=False)

            # Update checksum and datasum
            msgs.work('Updating checksum and datasum')
            hdul[0].add_checksum(override_datasum=False)
            hdul[1].add_datasum()
            hdul[1].add_checksum(override_datasum=True)
            hdul.flush(output_verify='fix')
            hdul.close()
            msgs.info('File {} produced.'.format(fits_out))
            if make_whitelight_image:
                image_hdul[0].add_datasum()
                image_hdul[1].add_datasum()
                image_hdul[0].add_checksum(override_datasum=True)
                image_hdul[1].add_checksum(override_datasum=True)
                image_hdul.writeto(image_out, overwrite=True, output_verify='fix')
            msgs.info('Image {} produced.'.format(image_out))

        elif instrument.startswith('IRDIS'):
            msgs.work('Fixing header for SPHERE/{} file {}'.format(instrument, fits_in))
            hdr = fitsfiles.header_from_fits_file(fits_in)
            if 'ESO DPR TECH' in hdr.keys():
                if str(hdr['ESO DPR TECH']).strip() == 'IMAGE,DUAL,CORONOGRAPHY':
                    msgs.work('Working with {} as observing technique'.format(str(hdr['ESO DPR TECH']).strip()))
                elif 'DUAL' in str(hdr['ESO DPR TECH']).strip():
                    msgs.error('{} needs to be tested'.format(str(hdr['ESO DPR TECH']).strip()))
                else:
                    msgs.error('Only DUAL imaging currently implemented')
            else:
                msgs.error('Cannot recognize the observing technique')

            # defining the two fits_out files:
            fits_out_index = [0, 1]
            fits_out_files = []
            for index in fits_out_index:
                fits_out_file = fits_out.replace('.fit', '_' + str(index) + '.fit')
                if os.path.exists(fits_out_file):
                    shutil.copy(fits_out_file, fits_out_file.replace('.fit', '_old.fit'))
                    msgs.warning('{} already exists. Backup created.'.format(fits_out))
                fitsfiles.new_fits_like(fits_in, [0], fits_out_file, overwrite=overwrite,
                                        fix_header=True, empty_primary_hdu=False)
                fits_out_files.append(fits_out_file)

            for index, fits_out_file in zip(fits_out_index, fits_out_files):
                hdul = fitsfiles.get_hdul(fits_out_file, 'update', checksum=True)
                hdr0 = hdul[0].header
                hdul[0].data = hdul[0].data[index, :, :]

                # Clean HISTORY keywords
                _ = cleaning_headers.history_keywords_delete(hdr0, verbose=True, in_place=True)
                # Try to guess coordinates
                if 'CRVAL1' not in hdr0.keys():
                    msgs.warning('CRVAL position keywords not preset')
                    if 'OBJECT' in hdr0.keys():
                        try:
                            object_coordinate = SkyCoord.from_name(str(hdr0['OBJECT']).strip())
                            ra_obj, dec_obj = object_coordinate.ra.degree, object_coordinate.dec.degree
                            if 'RA' in hdr0.keys() and 'DEC' in hdr0.keys():
                                pointing_coordinate = SkyCoord(float(hdr0['RA']), float(hdr0['DEC']), unit='deg')
                                msgs.work('Testing from separation from pointing position')
                                separation = object_coordinate.separation(pointing_coordinate).arcsec
                                if separation < 120.:
                                    msgs.info('Object - Pointing separation is {}'.format(separation))
                                    msgs.info('Updating CRVAL1 = {}'.format(ra_obj))
                                    msgs.info('Updating CRVAL2 = {}'.format(dec_obj))
                                    hdr0['CRVAL1'] = ra_obj
                                    hdr0['CRVAL2'] = dec_obj
                                    msgs.work('Updating CUNIT')
                                    hdr0['CUNIT1'] = 'deg'
                                    hdr0['CUNIT2'] = 'deg'
                                    msgs.work('Updating CTYPE')
                                    hdr0['CTYPE1'] = 'RA---TAN'
                                    hdr0['CTYPE2'] = 'DEC--TAN'
                                    msgs.work('Updating CRPIX')
                                    hdr0['CRPIX1'] = float(hdul[0].data.shape[1]) / 2.
                                    hdr0['CRPIX2'] = float(hdul[0].data.shape[0]) / 2.
                                    msgs.info('Updating CD1 and CD2')
                                    hdr0['CD1_1'] = hdr0['PIXSCAL'] * 2.778E-4 / 1000.
                                    hdr0['CD2_2'] = hdr0['PIXSCAL'] * 2.778E-4 / 1000.
                                    hdr0['CD1_2'] = 0.
                                    hdr0['CD2_1'] = 0.
                                    msgs.work('Updating RA, DEC')
                                    hdr0['RA'] = ra_obj
                                    hdr0.comments['RA'] = object_coordinate.ra.to_string(u.hour)
                                    hdr0['DEC'] = dec_obj
                                    hdr0.comments['DEC'] = object_coordinate.dec.to_string(u.degree, alwayssign=True)
                                else:
                                    msgs.warning('Object - Pointing separation is {}'.format(separation))
                                    msgs.warning('This is suspicious, CRVAL not updated')
                        except name_resolve.NameResolveError:
                            msgs.warning('Object {} not recognized'.format(str(hdr0['OBJECT']).strip()))
                            msgs.warning('CRVAL not updated')

                # Updating file prodcatg
                msgs.work('Updating PRODCATG to SCIENCE.IMAGE')
                hdr0['PRODCATG'] = str('SCIENCE.IMAGE')

                if 'PROGID' not in hdr0.keys():
                    msgs.warning('PROG_ID missing')
                    if 'ESO OBS PROG ID' in hdr0.keys():
                        msgs.info('Deriving PROG_ID from HIERARCH ESO OBS PROG ID')
                        msgs.work('Updating PROG_ID to {}'.format(str(hdr0['HIERARCH ESO OBS PROG ID'])))
                        hdr0['PROG_ID'] = str(hdr0['HIERARCH ESO OBS PROG ID'])

                # Update checksum and datasum
                msgs.work('Updating checksum and datasum')
                hdul[0].add_checksum(override_datasum=False)
                hdul.flush(output_verify='fix')
                hdul.close()
        else:
            msgs.warning('The Instrument {} is not supported \nThe file {} will not be processed'.format(instrument,
                                                                                                         fits_in))

    msgs.end()

    '''

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
            hdr0['PROV' + str(prov_number)] = hdr0[prov]
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
        image_hdul.append(fits.ImageHDU(to_ang * delta_wave_bin * np.nansum(hdul[1].data, axis=0, dtype=np.float_)))
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
                    download_archive.get_header_from_archive(file_id, text_file=file_id + '.hdr')
                    hdr_prov = fitsfiles.header_from_txt_file(file_id + '.hdr')
                    msgs.work('{} has EXPTIME of: {}'.format(file_id, hdr_prov['EXPTIME']))
                    EXPTIME_sec = EXPTIME_sec + hdr_prov['EXPTIME']
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
                n_pixels_psf = np.pi * (3) ** 2.
                five_sigma_nu = 5. * 3.34e4 * np.power((WAVELMAX + WAVELMIN) * to_ang / 2.,
                                                       2.) * whitelight_std * n_pixels_psf / (delta_wave_bin * to_ang)
                hdr0['ABMAGLIM'] = -2.5 * np.log10(five_sigma_nu / 3631.)
                msgs.warning('ABMAGLIM={}. This is most probably not correct at the moment'.format(hdr0['ABMAGLIM']))

 
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
    '''
