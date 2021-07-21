from astropy.coordinates import SkyCoord
from astropy.coordinates import name_resolve
from ESOAsg import msgs


def coordinates_from_object_name(object_name, strip_object_name=True):
    r"""Returns RA and Dec of an object given its name

    This is a wrapper around the the `astropy` class `SkyCoord`. The frame is fixed to ICRS and the
    coordinates are in degrees.

    Args:
        object_name (str): Name of the object that will be used to query
        strip_object_name (bool): if `True` the object name will be stripped, i.e. any trailing space
            will be removed

    Returns:
        tuple:
            tuple containing:
            * object_name(`str`_): stripped name of the object
            * object_ra (float): ra of the object in degrees
            * object_dec (float): dec of the object in degrees
            * object_coordinate (`SkyCoord`_): SkyCoord of the object

    """
    if strip_object_name:
        object_name = object_name.strip()
    object_coordinate = SkyCoord.from_name(object_name, frame='icrs')
    object_ra, object_dec = object_coordinate.ra.degree, object_coordinate.dec.degree
    return object_name, object_ra, object_dec, object_coordinate

