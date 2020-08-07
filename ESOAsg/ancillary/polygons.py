"""
This is a direct copy of `spherical geometry <https://github.com/spacetelescope/spherical_geometry>`_
It is copied here to avoid to deal with the C routines required for the full
package. This will change in the future.
* qd-library 2.3.7
The `spherical_geometry.polygon` module defines the `SphericalPolygon` class for
managing polygons on the unit sphere.
"""
# ToDo
# just simply import teh spherical_geometry in a give package
# This is used in:
# - ancillary.astro

import numpy as np
from numpy.core.umath_tests import inner1d

__all__ = ['SingleSphericalPolygon', 'SphericalPolygon']


def contours_to_polygons(contours, max_vertices=30):
    r"""Converts a contour into a polygon

    The resulting `polygons` is a list with N elements (with N matching the number of `contours`). Each elements
    contains a string defining the location in the sky of the polygon with RA, Dec, separated by commas and with the
    first RA, Dec pair that matches the last one (to close the polygon)

    Args:
        contours (list):
        max_vertices (int):

    Returns:
        list: list of strings defining each polygon

    """
    polygons = []
    for iii, contour in enumerate(contours):
        if len(contour) > max_vertices:
            contour_clean = contour[0: len(contour): int(len(contour) / max_vertices + 1)]
        else:
            contour_clean = contour
        # Construct the polygon as a string in the right format
        polygon = ' '.join(['%.4f, %.4f,' % (ra, dec) for ra, dec in contour_clean])[
                  :-1]  # Remove the last character, which is an unwanted extra comma
        polygons.append(polygon)
    return polygons


def two_d(vec):
    """
    Reshape a one dimensional vector so it has a second dimension
    """
    shape = list(vec.shape)
    shape.append(1)
    shape = tuple(shape)
    return np.reshape(vec, shape)


def midpoint(a, b):
    """
    Returns the midpoint on the great circle arc between *A* and *B*.

    Parameters
    ----------
    a, b : (*x*, *y*, *z*) triples or Nx3 arrays of triples
        The endpoints of the great circle arc.  It is assumed that
        these points are already normalized.

    Returns
    -------
    midpoint : (*x*, *y*, *z*) triple or Nx3 arrays of triples
        The midpoint between *A* and *B*, normalized on the unit
        sphere.
    """
    pp = (a + b) / 2.0
    # Now normalize...
    ll = np.sqrt(np.sum(pp * pp, axis=-1))
    ll = two_d(ll)
    return pp / ll


def triple_product(a, b, c):
    return inner1d(c, np.cross(a, b))


def normalize_vector(xyz, output=None):
    r"""
    Normalizes a vector so it falls on the unit sphere.

    Parameters
    ----------
    xyz : Nx3 array of vectors
        The input vectors

    output : Nx3 array of vectors, optional
        The array to store the results in.  If `None`, a new array
        will be created and returned.

    Returns
    -------
    output : Nx3 array of vectors
    """
    xyz = np.asanyarray(xyz, dtype=np.float64)

    if output is None:
        output = np.empty(xyz.shape, dtype=np.float64)

    ll = np.sqrt(np.sum(xyz * xyz, axis=-1))

    output = xyz / two_d(ll)

    return output


def lonlat_to_vector(lon, lat, degrees=True):
    r"""Converts a location on the unit sphere from longitude and latitude to an x, y, z vector

    This is taken from `spherical geometry <https://github.com/spacetelescope/spherical_geometry>`_

    Args:
        lon (any): longitude as scalar or array
        lat (any): latitude as scalar or array
        degrees (bool, optional): if `True` `lon` and `lat` are in decimal degrees, otherwise in radians

    Returns:
        tuple: x, y, z in format of scalars or 1-D arrays of the same length

    .. notes::

       Where longitude is *l* and latitude is *b*:

    .. math::
        x = \cos l \cos b

        y = \sin l \cos b

        z = \sin b
    """
    lon = np.asanyarray(lon)
    lat = np.asanyarray(lat)

    if degrees:
        lon_rad = np.deg2rad(lon)
        lat_rad = np.deg2rad(lat)
    else:
        lon_rad = lon
        lat_rad = lat

    cos_lat = np.cos(lat_rad)

    return (
        np.cos(lon_rad) * cos_lat,
        np.sin(lon_rad) * cos_lat,
        np.sin(lat_rad))


def vector_to_lonlat(x, y, z, degrees=True):
    r"""
    Converts a vector to longitude and latitude.

    Parameters
    ----------
    x, y, z : scalars or 1-D arrays
        The input vectors

    degrees : bool, optional
        If `True` (default) the result is returned in decimal degrees,
        otherwise radians.

    Returns
    -------
    lon, lat : tuple of scalars or arrays of the same length

    Notes
    -----
    Where longitude is *l* and latitude is *b*:

    .. math::
        l = \arctan2(y, x)

        b = \arctan2(z, \sqrt{x^2 + y^2})
    """
    x = np.asanyarray(x, dtype=np.float64)
    y = np.asanyarray(y, dtype=np.float64)
    z = np.asanyarray(z, dtype=np.float64)

    lon = np.arctan2(y, x)
    lon = np.remainder(lon, 2.0 * np.pi)

    lat = np.arctan2(z, np.sqrt(x ** 2 + y ** 2))
    result = (lon, lat)

    if degrees:
        return np.rad2deg(result[0]), np.rad2deg(result[1])
    else:
        return result


class SingleSphericalPolygon(object):
    r"""
    Polygons are represented by both a set of points (in Cartesian
    (*x*, *y*, *z*) normalized on the unit sphere), and an inside
    point.  The inside point is necessary, because both the inside and
    outside of the polygon are finite areas on the great sphere, and
    therefore we need a way of specifying which is which.
    """

    def __init__(self, points, inside=None):
        r"""
        Parameters
        ----------
        points : An Nx3 array of (*x*, *y*, *z*) triples in vector space
            These points define the boundary of the polygon.
            It may contain zero points, in which it defines the null
            polygon.  It may not contain one, two or three points.
            Four points are needed to define a triangle, since the
            polygon must be closed.
        inside : An (*x*, *y*, *z*) triple, optional
            This point must be inside the polygon.  If not provided, an
            interior point will be calculated
        """
        if len(points) == 0:
            # handle special case of initializing with an empty list of
            # vertices (ticket #1079).
            self._inside = np.zeros(3)
            self._points = np.asanyarray(points)
            return

        if not np.array_equal(points[0], points[-1]):
            points = list(points[:])
            points.append(points[0])

        if len(points) < 3:
            raise ValueError("Polygon made of too few points")

        self._points = points = np.asanyarray(points)
        new_inside = self._find_new_inside()

        if inside is None:
            self._inside = np.asanyarray(new_inside)

            if not self.is_clockwise():
                self._points = points[::-1]
        else:
            self._inside = np.asanyarray(inside)

            if self.contains_point(new_inside) != self.is_clockwise():
                self._points = points[::-1]

    def __len__(self):
        return 1

    def __repr__(self):
        return '%s(%r, %r)' % (self.__class__.__name__,
                               self.points, self.inside)

    def __iter__(self):
        """
        Iterate over all base polygons that make up this multi-polygon
        set.
        """
        yield self

    @property
    def points(self):
        """
        The points defining the polygon.  It is an Nx3 array of
        (*x*, *y*, *z*) vectors.  The polygon will be explicitly
        closed, i.e., the first and last points are the same.
        """
        return self._points

    @property
    def inside(self):
        """
        Get the inside point of the polygon.
        """
        return self._inside

    def is_clockwise(self):
        """
        Return True if the points in this polygon are in clockwise order.
        The normal vector to the two arcs containing a vertes points outward
        from the sphere if the angle is clockwise and inward if the angle is
        counter-clockwise. The sign of the inner product of the normal vector
        with the vertex tells you this. The polygon is ordered clockwise if
        the vertices are predominantly clockwise and counter-clockwise if
        the reverse.
        """

        points = np.vstack((self._points, self._points[1]))
        a = points[:-2]
        b = points[1:-1]
        c = points[2:]
        orient = triple_product(a-b, c-b, b)
        return np.sum(orient) > 0.0

    def to_lonlat(self):
        """
        Convert `SingleSphericalPolygon` footprint to longitude and latitutde.
        Returns
        -------
        lon, lat : list of float
            List of *lon* and *lat* in degrees corresponding
            to `points`.
        """
        if len(self.points) == 0:
            return np.array([])
        return vector_to_lonlat(self.points[:, 0], self.points[:, 1],
                                self.points[:, 2], degrees=True)

    @classmethod
    def from_lonlat(cls, lon, lat, center=None, degrees=True):
        r"""
        Create a new `SingleSphericalPolygon` from a list of
        (*longitude*, *latitude*) points.
        Parameters
        ----------
        lon, lat : 1-D arrays of the same length
            The vertices of the polygon in longitude and
            latitude.
        center : (*lon*, *lat*) pair, optional
            A point inside of the polygon to define its inside.
        degrees : bool, optional
            If `True`, (default) *lon* and *lat* are in decimal degrees,
            otherwise in radians.
        Returns
        -------
        polygon : `SingleSphericalPolygon` object
        """
        # Convert to Cartesian
        x, y, z = lonlat_to_vector(lon, lat, degrees=degrees)

        points = np.dstack((x, y, z))[0]

        if center is not None:
            center = lonlat_to_vector(*center, degrees=degrees)

        return cls(points, center)

    def _find_new_inside(self):
        """
        Finds an acceptable inside point inside of *points* that is
        also inside of *polygons*.
        """
        npoints = len(self._points)
        if npoints > 4:
            points = np.vstack((self._points, self._points[1]))
            a = points[:-2]
            b = points[1:-1]
            c = points[2:]
            orient = triple_product(a-b, c-b, b)
            if np.sum(orient) < 0.0:
                orient = -1.0 * orient
            midpoint_ac = midpoint(a, c)
            candidate = max(zip(orient, midpoint_ac), key=lambda x: x[0])
            inside = candidate[1]
        else:
            # Fall back on computing the mean point
            inside = self._points.mean(axis=0)
            normalize_vector(inside, output=inside)
        return inside


class SphericalPolygon(SingleSphericalPolygon):
    r"""
    Polygons are represented by both a set of points (in Cartesian
    (*x*, *y*, *z*) normalized on the unit sphere), and an inside
    point.  The inside point is necessary, because both the inside and
    outside of the polygon are finite areas on the great sphere, and
    therefore we need a way of specifying which is which.
    This class contains a list of disjoint closed polygons.
    """

    def __init__(self, init, inside=None):
        r"""
        Parameters
        ----------
        init : object
            May be either:
               - A list of disjoint `SphericalPolygon` objects.
               - An Nx3 array of (*x*, *y*, *z*) triples in Cartesian
                 space.  These points define the boundary of the
                 polygon.
                 It may contain zero points, in which it defines the
                 null polygon.  It may not contain one or two points.
        inside : An (*x*, *y*, *z*) triple, optional
            If *init* is an array of points, this point must be inside
            the polygon.  If it is not provided, one will be created.
        """
        for polygon in init:
            if not isinstance(polygon, (SphericalPolygon, SingleSphericalPolygon)):
                break
        else:
            self._polygons = tuple(init)
            return

        self._polygons = (SingleSphericalPolygon(init, inside),)

    def __iter__(self):
        """
        Iterate over all base polygons that make up this multi-polygon
        set.
        """
        for polygon in self._polygons:
            for subpolygon in polygon:
                yield subpolygon

    @property
    def points(self):
        """
        The points defining the polygons.  It is an iterator over
        disjoint closed polygons, where each element is an Nx3 array
        of (*x*, *y*, *z*) vectors.  Each polygon is explicitly
        closed, i.e., the first and last points are the same.
        """
        for polygon in self:
            yield polygon.points

    @property
    def inside(self):
        """
        Iterate over the inside point of each of the polygons.
        """
        for polygon in self:
            yield polygon.inside

    @property
    def polygons(self):
        """
        Get a sequence of all of the subpolygons.  Each subpolygon may
        itself have subpolygons.  To get a flattened sequence of all
        base polygons, use `iter_polygons_flat`.
        """
        return self._polygons

    def to_lonlat(self):
        """
        Convert the `SphericalPolygon` footprint to longitude and latitude
        coordinates.
        Returns
        -------
        polyons : iterator
            Each element in the iterator is a tuple of the form (*lon*,
            *lat*), where each is an array of points.
        """
        for polygon in self:
            yield polygon.to_lonlat()

    @classmethod
    def from_lonlat(cls, lon, lat, center=None, degrees=True):
        # TODO Move into SingleSphericalPolygon
        r"""
        Create a new `SphericalPolygon` from a list of
        (*longitude*, *latitude*) points.
        Parameters
        ----------
        lon, lat : 1-D arrays of the same length
            The vertices of the polygon in longitude and
            latitude.
        center : (*lon*, *lat*) pair, optional
            A point inside of the polygon to define its inside.
        degrees : bool, optional
            If `True`, (default) *lon* and *lat* are in decimal degrees,
            otherwise in radians.
        Returns
        -------
        polygon : `SphericalPolygon` object
        """
        polygon = SingleSphericalPolygon.from_lonlat(lon, lat,
                                                     center, degrees)
        return cls((polygon,))



