r"""Class to work on images

This includes some photutils code: `https://photutils.readthedocs.io/en/stable/index.html`
"""

from astropy.stats import sigma_clipped_stats
from photutils import make_source_mask
import numpy as np

from ESOAsg import msgs
from ESOAsg.core import fitsfiles
from ESOAsg.ancillary import checks


class Images:
    r"""A class used to define and make simple operations on images

    This allows to perform some basic checks on a 2D image.

    Attributes:

    Methods:

    """

    def __init__(self, data=None, errors=None, background=None, quality=None, bad_pixels_mask=None, sources_mask=None):
        r"""Instantiate the class image

        """
        self.data = data
        self.errors = errors
        self.background = background
        self.quality = quality
        self.bad_pixels_mask = bad_pixels_mask
        self.sources_mask = sources_mask

        # check that all attributes are images
        if not self._check_image():
            msgs.errors('Not a valid Image')

    def from_fits(self, fits_name=None, data_hdu=None, errors_hdu=None, background_hdu=None, quality_hdu=None,
                  bad_pixels_mask_hdu=None, sources_mask_hdu=None):
        r"""Load an Image from a fits file

        Args:
            fits_name (`str`):
                fits file name where the image is stored
            data_hdu (`int'):
                hdu where the data are
            errors_hdu (`int'):
                hdu where the error are
            background_hdu (`int'):
                hdu where the background is
            quality_hdu (`int'):
                hdu where the quality is
            bad_pixels_mask_hdu (`int'):
                hdu where the bad_pixel_mask is
            sources_mask_hdu (`int'):
                hdu where the sources_mask is

        """

        if not checks.fits_file_is_valid(fits_name):
            msgs.error('{} not a valid fits file'.format(fits_name))
        hdul = fitsfiles.get_hdul(fits_name)

        for attribute, value in zip(self.__dict__.keys(), self.__dict__.values()):
            which_hdu = vars()[attribute + '_hdu']
            if which_hdu is not None:
                assert isinstance(which_hdu, (int, np.int_)), r'{} must be an integer'.format(
                    attribute + '_hdu')
                setattr(self, attribute, hdul[which_hdu].data)

        if not self._check_image():
            msgs.errors('Not a valid Image')

    def _check_image(self):
        r"""Check that all entries are valid
        """
        good_image = True
        for attribute, value in zip(self.__dict__.keys(), self.__dict__.values()):
            if value is not None:
                if not checks.image2d_is_valid(value):
                    msgs.warning('No valid {} present'.format(attribute))
                    good_image = False
        return good_image

    def calc_background(self, method='median', find_sources=True, sigma=3.):
        r"""

        Args:
            sigma:
            find_sources
            method:

        """
        if self.data is None:
            msgs.warning('There are no data to calculate the background')
            return False
        if self.background is not None:
            msgs.warning('The background will be overwritten')
        if find_sources:
            if self.sources_mask is not None:
                msgs.warning('The sources_mask will be overwritten')
                self.find_sources()
            # ToDo
            msgs.error('Not implemented yet')
            msgs.work('Not impleSearching for bright sources with a {}-sigma clipping algorithm'.format(sigma))
        msgs.work('Finding background with a {}-sigma clipping algorithm'.format(sigma))
        mean, median, std = sigma_clipped_stats(self.data, mask=self.get_full_mask(), sigma=sigma)
        msgs.info('Mean={}, Median={}, Std={}'.format(mean, median, std))
        return mean, median, std

    def get_full_mask(self, mask_bad_pixels=True, mask_quality=True, mask_sources=True, mask_nans=True):
        r"""

        Returns:

        """
        if self.data is None:
            msgs.warning('data is empty')
            return None
        else:
            full_mask = np.zeros_like(self.data, dtype=np.bool)
            if mask_bad_pixels:
                if self.bad_pixels_mask is not None:
                    full_mask[self.bad_pixels_mask] = True
            if mask_quality:
                if self.quality is not None:
                    full_mask[self.quality] = True
            if mask_sources:
                if self.sources_mask is not None:
                    full_mask[self.sources_mask] = True
            if mask_nans:
                full_mask[np.isnan(self.data)] = True
        return full_mask

    def find_sources(self, sigma=3., remove_background=True, mask_radius=2.):
        r"""Filling source mask with a sigma detection algorithm

        """
        if self.data is None:
            msgs.warning('data is empty, mask_sources set to None')
            return None

