r"""Class to work on images

This includes some photutils code: `https://photutils.readthedocs.io/en/stable/index.html`
"""

from astropy.stats import sigma_clipped_stats
from astropy.stats import SigmaClip
from photutils import make_source_mask
from photutils import Background2D, MedianBackground
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

    def calc_background(self, method='median', nsigma=3., find_sources=True, src_nsigma=3.,
                        src_npixels=5, src_dilate_size=2):
        r"""Calculate background given data.

        Args:
            method (`str`):
                method used to calculate the background. `median`, `average` , or `bg2d`
            nsigma (`float`):
                rejection sigma for the background estimation
            find_sources ('bool'):
                if `True` the sources_mask will be updated
            src_nsigma (`int`):
                sources detection sigma limit
            src_npixels (`int`):
                number of pixels to consider a source detected
            src_dilate_size (`int`):
                pixels buffer around sources to be included in the sources_mask

        """
        if self.data is None:
            msgs.warning('There are no data to calculate the background')
            return False
        if self.background is not None:
            msgs.warning('The background will be overwritten')
        if find_sources:
            self.find_sources(nsigma=src_nsigma, npixels=src_npixels, dilate_size=src_dilate_size)

        if method == 'bg2d':
            sigma_clip = SigmaClip(sigma=nsigma)
            bkg_estimator = MedianBackground()
            self.background = Background2D(self.data, (30, 30), filter_size=(3, 3), sigma_clip=sigma_clip,
                                           bkg_estimator=bkg_estimator, mask=self.get_full_mask)
        elif method == 'median' or method == 'mean':
            bg_mean, bg_median, bg_std = self.get_clean_stats(nsigma=nsigma)
            if method == 'mean' :
                self.background = np.full_like(self.data, fill_value=bg_mean, dtype=np.float_)
            else:
                self.background = np.full_like(self.data, fill_value=bg_median, dtype=np.float_)
        else:
            msgs.error('Not a valid method for background subtraction')

        return True

    def find_sources(self, nsigma=3., npixels=5, dilate_size=2):
        r"""Filling source mask with photutils.make_source_mask

        """
        if self.data is None:
            msgs.warning('data is empty, mask_sources set to None')
            return False
        if self.sources_mask is not None:
            msgs.warning('sources_mask will be overwritten')

        if self.background is not None:
            data_clean = np.copy(self.data - self.background)
        else:
            data_clean = np.copy(self.data)

        msgs.work('Searching for sources {}-sigma above the background'.format(nsigma))
        self.sources_mask = make_source_mask(data_clean, nsigma=nsigma, npixels=npixels, dilate_size=dilate_size)

        return True

    def get_full_mask(self, include_bad_pixels=True, include_quality=True, include_sources=True, include_nans=True):
        r""" Sum over all the masks of an Image

        Args:
            include_bad_pixels (`bool`):
                include mask_bad_pixels in the final mask
            include_quality (`bool`):
                include mask_quality in the final mask
            include_sources (`bool`):
                include sources_mask in the final mask
            include_nans (`bool`):
                mask `NaNs` in the data and add them to the final mask


        Returns:
            full_mask (`np.array`, `bool`):
                the mask resulting from the combination of all the `Trues` switches.
        """
        if self.data is None:
            msgs.warning('data is empty')
            return None
        else:
            full_mask = np.zeros_like(self.data, dtype=np.bool)
            if include_bad_pixels:
                if self.bad_pixels_mask is not None:
                    full_mask[self.bad_pixels_mask] = True
            if include_quality:
                if self.quality is not None:
                    full_mask[self.quality] = True
            if include_sources:
                if self.sources_mask is not None:
                    full_mask[self.sources_mask] = True
            if include_nans:
                full_mask[np.isnan(self.data)] = True
        return full_mask

    def get_clean_stats(self, nsigma=3.):
        """Calculate stats of data removing all the masked pixels
        """
        if self.data is None:
            msgs.error('data is empty')
        mean, median, std = sigma_clipped_stats(self.data, mask=self.get_full_mask(), sigma=nsigma)
        return mean, median, std
