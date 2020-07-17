from ESOAsg import msgs


def from_polygons(polygons=None, open_link=False, show_link=False):
    if polygons is not None:
        for iii, polygon in enumerate(polygons):
            url = 'http://archive.eso.org/scienceportal/home?' + 'poly=' + polygon + '&sort=-obs_date'
            if show_link:
                msgs.info('ASP link to region N.{} is:\n {}\n'.format(np.str(iii + 1), url))
            if open_link:
                webbrowser.open(url)


def from_radec(positions, radius=None, open_link=False, show_link=False):
    r"""Query the ESO ASP service given a position in RA and Dec.

     The `positions` value (or list) needs to be given as an `astropy.coordinates.SkyCoord` object.

    Args:
        positions (`astropy.coordinates.SkyCoord`):
            Coordinates (or list of coordinates) of the sky you want to query in the format of an
            `astropy.coordinates.SkyCoord` object. For further detail see here:
            `astropy coordinates <https://docs.astropy.org/en/stable/coordinates/>`_
        radius (`float`):
            Search radius you want to query in arcseconds. Note that in case `None` is given, the query will be
            performed with the `CONTAINS(POINT('',RA,Dec), s_region)` clause instead of the
            `CONTAINS(s_region,CIRCLE('',RA,Dec,radius/3600.))` one. See here for further examples:
            `tap obs examples <http://archive.eso.org/tap_obs/examples>`_
        open_link (`bool`):
            open a link to the ASP page
        show_link (`bool`):
            show the link on the terminal

    """
    # Check inputs:
    # Working on positions
    if isinstance(positions, list):
        positions_list = positions
    else:
        positions_list = [positions]
    for position in positions_list:
        assert isinstance(position, coordinates.SkyCoord), r'Input positions not a SkyCoord object'
    # Working on radius
    if radius is not None:
        if isinstance(radius, int):
            radius = float(radius)
        else:
            assert isinstance(radius, float), r'Input radius is not a number'

    for iii, position in enumerate(positions_list):
        position.transform_to(ICRS)
        radec_string = np.str(np.float_(position.ra.degree)) + ',' + np.str(np.float_(position.dec.degree))
        radius_string = np.str(radius / 3600.)
        url = 'http://archive.eso.org/scienceportal/home?' + 'pos=' + radec_string + '&r=' + radius_string + \
              '&sort=-obs_date'
        if show_link:
            msgs.info('ASP link to region N.{} is:\n {}\n'.format(np.str(iii + 1), url))
        if open_link:
            webbrowser.open(url)



































