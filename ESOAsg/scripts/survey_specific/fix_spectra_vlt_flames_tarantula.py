r"""
fix_spectra_vlt_flames_tarantula
================================
Script to make data produced by the `VLT/FLAMES TARANTULA group <https://www.roe.ac.uk/~cje/tarantula/>`_
(almost) phase 3 compliant

.. topic:: Inputs:

    - **input_fits** - Input fits file

"""

import argparse

from ESOAsg import __version__

EXAMPLES = str(r"""EXAMPLES:""" + """\n""" + """\n""" +
               r"""TO BE DONE """ + """\n""" +
               r""" """)


def parser(options=None):
    parser = argparse.ArgumentParser(
        description=r"""Manipulate VLT/FLAMES data produced by the """ +
                    r"""VLT/FLAMES TARANTULA group (https://www.roe.ac.uk/~cje/tarantula/) """ +
                    r"""and make them compliant with the  """ +
                    r"""`ESO Phase 3 Standard` (https://www.eso.org/sci/observing/phase3/p3sdpstd.pdf). """ +
                    """\n""" + """\n""" +
                    r"""Note that, if `suffix` is not defined, the fits file will be overwritten.  """ +
                    """\n""" + """\n""" +
                    r"""This uses ESOAsg version {:s}""".format(__version__),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument("input_fits", nargs="+", type=str,
                        help=r"Input spectra")
    parser.add_argument("-s", "--suffix", nargs="+", type=str, default=None,
                        help=r"Suffix to be added to file names the will be used as output." +
                             r"If it is not set, the input files will be overwritten")

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
    from ESOAsg.core import fitsfiles
    from ESOAsg import msgs

    # Cleaning input lists
    input_fits_files = cleaning_lists.make_list_of_fits_files(args.input_fits)
    if args.suffix is None:
        overwrite = False
        msgs.warning('The file will overwrite the input files')
    else:
        overwrite = True
    suffix_string = cleaning_lists.make_string(args.suffix)
    output_fits_files = cleaning_lists.make_list_of_fits_files_and_append_suffix(input_fits_files,
                                                                                 suffix=suffix_string)

    for fits_in, fits_out in zip(input_fits_files, output_fits_files):
        if os.path.exists(fits_out):
            shutil.copy(fits_out, fits_out.replace('.fit', '_old.fit'))
            msgs.warning('{} already exists. Backup created.'.format(fits_out))

        msgs.work('Fixing file: {}'.format(fits_in))

        # Create a copy of the file where there is a primary HDU and data are in the 'DATA" HDU
        msgs.work('Reshaping spectra into PrimaryHEADER and Data Header')

        hdul_in = fitsfiles.get_hdul(fits_in, 'readonly', checksum=True)
        hdu0_out = fits.PrimaryHDU(header=hdul_in[0].header.copy())
        hdul_out = fits.HDUList([hdu0_out])
        FLUX = hdul_in['SKY SUBTRACTED'].data
        ERR = hdul_in['ORIGINAL ERROR'].data
        CRVAL1 = hdul_in['SKY SUBTRACTED'].header['CRVAL1']
        CRPIX1 = hdul_in['SKY SUBTRACTED'].header['CRPIX1']
        CDELT1 = hdul_in['SKY SUBTRACTED'].header['CRPIX1']
        WAVE = (CDELT1*np.arange(len(FLUX))) + CRVAL1

        # double precision:
        col_dtype = np.float64
        col_format = str(len(np.array(WAVE, dtype=col_dtype))) + 'D'

        col1 = fits.Column(name='WAVE', format=col_format, unit='angstrom',
                           array=[np.array(WAVE, dtype=col_dtype)])
        col2 = fits.Column(name='FLUX', format=col_format, unit='adu',
                           array=[np.array(FLUX, dtype=col_dtype)])
        col3 = fits.Column(name='ERR', format=col_format, unit='adu',
                           array=[np.array(ERR, dtype=col_dtype)])

        hdu1_out = fits.BinTableHDU.from_columns([col1, col2, col3], nrows=1)
        hdu1_out.header['EXTNAME'] = 'SPECTRUM'
        hdu1_out.header['TUCD1'] = 'em.wl'
        hdu1_out.header['TUCD2'] = 'phot.flux.density;em.wl;stat.uncalib'
        hdu1_out.header['TUCD3'] = 'stat.error;phot.flux.density;em.wl;stat.uncalib'

        hdu1_out.header['EXTNAME'] = 'SPECTRUM'

        hdul_out.append(hdu1_out)

        msgs.work('Updating checksum and datasum')
        hdul_out[0].add_datasum()
        hdul_out[1].add_datasum()
        hdul_out[0].add_checksum(override_datasum=True)
        hdul_out[1].add_checksum(override_datasum=True)
        msgs.work('Saving changes.')
        hdul_out.writeto(fits_out, overwrite=True)
        hdul_out.close()
        hdul_in.close()

        msgs.newline()
        msgs.info('File {} produced.'.format(fits_out))
        msgs.end()


def not_to_be_used():
    '''
    msgs.work('Refactoring data structure.')
    msgs.work('Taking data from the file {}'.format(str(input_fits.replace('.fits', '.dat'))))
    spectra = ascii.read(input_fits.replace('.fits', '.dat'))
    WAVE = spectra['col1']
    FLUX = spectra['col2']
    ERR = spectra['col3']



    msgs.work('Refactoring Header.')

    # Now working on HDU1
    change1_cards, change1_values, change1_comments = [], [], []
    # Object names
    object_name: str = input_fits.replace('_av.fits', '').split('/')[-1]
    if object_name == 'HD193281A':
        begin_name = 'HD193281'
    elif object_name == 'HD193281B':
        begin_name = 'HD193281'
    elif object_name == 'HD193256':
        begin_name = 'HD193281'
    else:
        begin_name = object_name
    msgs.work('Searching for missing KEYWORDS')
    fullHeader = ascii.read('PROV.txt')
    this_object = np.zeros(len(fullHeader['FILE']), dtype=np.bool)
    for file_name, i in zip(fullHeader['FILE'], range(0, len(fullHeader['FILE']))):
        this_object[i] = file_name.startswith(begin_name)
    print(fullHeader['FILE'][this_object].data)


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

    """
    result_from_query = download_archive.query_from_radec(radec, radius=600., instrument='MUSE')
    download_archive.download(result_from_query['dp_id'])
    """

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
    change1_cards.append('TUCD1')
    change1_values.append('em.wl;obs.atmos')
    change1_comments.append('Air wavelength')
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
    change1_values.append('phot.flux.density;em.wl; ;meta.main')
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

    #SPEC_VAL
    change1_cards.append('SPEC_VAL')
    change1_values.append(0.1*(np.nanmax(WAVE)+np.nanmin(WAVE))/2.)
    change1_comments.append('[nm]')

    #SPEC_BW
    change1_cards.append('SPEC_BW')
    change1_values.append(0.1*(np.nanmax(WAVE)-np.nanmin(WAVE)))
    change1_comments.append('[nm]')

    #VOCLASS
    change1_cards.append('VOCLASS')
    change1_values.append('SPECTRUM v1.0')
    change1_comments.append('VO Data Model')

    #VOPUB
    change1_cards.append('VOPUB')
    change1_values.append('ESO/SAF')
    change1_comments.append(' ')

    #TITLE
    change1_cards.append('TITLE')
    change1_values.append(object_name + str(' MUSE spectrum'))
    change1_comments.append(' ')

    #NELEM
    change1_cards.append('NELEM')
    change1_values.append(np.int(len(np.array(WAVE, dtype=col_dtype))))
    change1_comments.append(' ')

    change1_cards.append('TELAPSE')
    change1_values.append(86400.*np.float_(np.nanmax(fullHeader['MJD-END'][this_object].data)-
                                           np.nanmin(fullHeader['MJD-OBS'][this_object].data)))
    change1_comments.append('[s]')

    change1_cards.append('TMID')
    change1_values.append(0.5*(np.nanmax(fullHeader['MJD-END'][this_object].data)+np.nanmin(fullHeader['MJD-OBS'][
                                                                                           this_object].data)))
    change1_comments.append(' ')


    # Diameter of the aperture used to extract the spectra in deg.
    change1_cards.append('APERTURE')
    if object_name == '[B86] 133':
        change1_values.append(2.*4./3600.)
    if object_name == 'HD 193256':
        change1_values.append(2.*4.6/3600.)
    if object_name == 'HD 193281B':
        change1_values.append(2.*1.2/3600.)
    else:
        change1_values.append(2.*6./3600.)
    change1_comments.append('[deg] Aperture width')

    for cards, values, comments in zip(change1_cards, change1_values, change1_comments):
        msgs.work('HDUL[1], Updating {}'.format(cards))
        fitsfiles.add_header_card(hdul[1].header, cards, values, comment=comments)

    # Now working on HDU0
    change0_cards, change0_values, change0_comments = [], [], []

    change0_cards.append('MJD-OBS')
    change0_values.append(np.nanmin(fullHeader['MJD-OBS'][this_object].data))
    change0_comments.append(' ')

    change0_cards.append('MJD-END')
    change0_values.append(np.nanmax(fullHeader['MJD-END'][this_object].data))
    change0_comments.append(' ')

    change0_cards.append('EXPTIME')
    change0_values.append(np.nansum(fullHeader['TEXPTIME'][this_object].data))
    change0_comments.append(' ')

    change0_cards.append('TEXPTIME')
    change0_values.append(np.nansum(fullHeader['TEXPTIME'][this_object].data))
    change0_comments.append(' ')

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
    change0_values.append(0.1*np.nanmin(WAVE))
    change0_comments.append('[nm] Minimum wavelength')

    change0_cards.append('WAVELMAX')
    change0_values.append(0.1*np.nanmax(WAVE))
    change0_comments.append('[nm] Maximum wavelength')

    change0_cards.append('SPEC_BIN')
    change0_values.append(0.1*hdr1['CDELT1'])
    change0_comments.append('[nm] Wavelength bin sizes')

    change0_cards.append('REFERENC')
    change0_values.append('doi:10.1051/0004-6361/201936178')
    change0_comments.append(' ')

    change0_cards.append('OBSTECH')
    change0_values.append('IFU')
    change0_comments.append('Technique of observation')

    change0_cards.append('PROCSOFT')
    change0_values.append('muse/2.6')
    change0_comments.append(' ')

    change0_cards.append('SPECSYS')
    change0_values.append('BARYCENT')
    change0_comments.append(' ')

    change0_cards.append('PROG_ID')
    change0_values.append('099.D-0623(A)')
    change0_comments.append(' ')

    change0_cards.append('FLUXCAL')
    change0_values.append('ABSOLUTE')
    change0_comments.append(' ')

    change0_cards.append('EXT_OBJ')
    change0_values.append(False)
    change0_comments.append(' ')

    change0_cards.append('M_EPOCH')
    change0_values.append(False)
    change0_comments.append(' ')

    change0_cards.append('CONTNORM')
    change0_values.append(False)
    change0_comments.append(' ')

    change0_cards.append('TOT_FLUX')
    change0_values.append(False)
    change0_comments.append(' ')

    change0_cards.append('FLUXERR')
    change0_values.append(-2)
    change0_comments.append(' ')

    change0_cards.append('SNR')
    msgs.info('Calculated SNR={}'.format(np.nanmedian(FLUX[FLUX > 0.]/ERR[FLUX > 0.])))
    change0_values.append(np.nanmedian(FLUX[FLUX > 0.]/ERR[FLUX > 0.]))
    change0_comments.append(' ')

    change0_cards.append('NCOMBINE')
    change0_values.append(len(fullHeader['MJD-END'][this_object].data))
    change0_comments.append(' ')

    change0_cards.append('DATE')
    change0_values.append(fullHeader['DATE'][this_object].data[0])
    change0_comments.append(' ')

    #SPEC_RES
    change0_cards.append('SPEC_RES')
    change0_values.append(np.int(3000))
    change0_comments.append('Spectral resolution')

    # add ancillary file
    ancillary_file = input_fits.replace('_av.fits', '.jpg')

    change0_cards.append('ASSON1')
    change0_values.append(ancillary_file)
    change0_comments.append(' ')

    change0_cards.append('ASSOM1')
    change0_values.append(hashlib.md5(open(ancillary_file, 'rb').read()).hexdigest())
    change0_comments.append(' ')

    change0_cards.append('ASSOC1')
    change0_values.append('ANCILLARY.PREVIEW')
    change0_comments.append(' ')

    # PROVENANCEs
    INDEX = 1
    for filename in fullHeader['FILE'][this_object].data:
        PROV = filename.split("/")[1]
        change0_cards.append('PROV'+str(INDEX))
        change0_values.append(PROV)
        change0_comments.append(' ')
        INDEX = INDEX+1

    INDEX = 1
    for OBID in np.unique(np.array(fullHeader['OBID1'][this_object].data)):
        change0_cards.append('OBID'+str(INDEX))
        change0_values.append(OBID)
        change0_comments.append(' ')
        INDEX = INDEX+1

    for cards, values, comments in zip(change0_cards, change0_values, change0_comments):
        msgs.work('HDUL[0], Updating {}'.format(cards))
        fitsfiles.add_header_card(hdul[0].header, cards, values, comment=comments)


    '''