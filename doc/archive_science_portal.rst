.. _archive-science-portal:


======================
Archive Science Portal
======================

The `Archive Science Portal <http://archive.eso.org/scienceportal/home>`_  provides the primary entry point 
for interactive searching and browsing of the ESO Science Archive using an intuitive interactive user interface 
that supports iterative queries.

The archive contains both pipeline-processed data streams of many La Silla-Paranal instruments and a large 
collection of coherent data sets provided by users of ESO telescopes' time including the results of Public Surveys. 

The module `archive_science_portal` provides some simple functions to create an `url` corresponding to your search.

Please, when you use data from the archive, follow the :doc:`datapolicy`.

Overview
========

The module `archive_science_portal` uses `core.asp_queries` to create the `url` that links
to your query at ASP.
The basic idea is to start from the `base url <http://archive.eso.org/scienceportal/home?>`_ and keep adding conditions until your search criteria is satisfied.
For instance, to search for all the `FORS2` and `XSHOOTER` spectra present in the archive:
::

    from ESOAsg.core import asp_queries
    # define the lists of instruments and data_types to be queried
    instruments = ['FORS2', 'XSHOOOTER']
    data_types = ['spectrum']
    # create the url
    url = '{0}{1}{2}'.format(asp_queries.base_url(),
                             asp_queries.condition_instruments(instruments),
                             asp_queries.condition_data_types(data_types, connector='&'))
    # run the query
    asp_queries.run_query(url)

This will open `this link <http://archive.eso.org/scienceportal/home?ins_id=FORS2,XSHOOOTER&dp_type=SPECTRUM>`_ on your browser.

The module `archive_science_portal` provides a set of functions to avoid the actual creation of the url by the user and to just obtain the link(s) to the page she/he is interested in.


Some examples
=============

query to a specific position in the sky
---------------------------------------

In early 2019 the red supergiant star `Betelgeuse <http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=Betelgeuse>`_ began to dim (see, `ESO science release #2003 <https://www.eso.org/public/news/eso2003/>`_). 
To open a page containing all the reduced spectra and data cubes associated to the star, you can run:
::

    from astropy.coordinates import SkyCoord
    from ESOAsg import archive_science_portal
    # resolve the star location from the name
    star_position = SkyCoord.from_name('Betelgeuse')
    # run the query with a search radius of 5 arcseconds anse selection spectrum and cube
    # as data types
    archive_science_portal.query_from_radec(star_position, radius=5., open_link=True, 
                                            data_types=['spectrum', 'cube'])

This will open `this page <http://archive.eso.org/scienceportal/home?pos=88.79293899,7.407064&r=0.001388888888888889&dp_type=SPECTRUM,CUBE&sort=-obs_date>`_ on your broswer.
From there you can select the data to download (currently, Aug. 2020, 57 spectra from `HARPS`, `UVES`, and `XSHOOTER` and 10 datacubes from `ALMA` could be retrieved).

query to a polygon in the sky
-----------------------------

The `version 2.0 of the archive science portal <https://archive.eso.org/cms/eso-archive-news/eso-archive-science-portal-2-0-released.html>`_  offers the possibility of searching data withing irregular polygons.
This is a relevant tool to detect electromagnetic counterpart of gravitational wave detections. 

This possibility is provided by the `archive_science_portal.query_from_polygons` function. 
For instance, following the procedure described in this `jupyter notebook <https://github.com/EmAstro/ESOAsg/blob/master/doc/notebooks/HOWTO_getDataFromGWContours.ipynb>`_, you can investigate which data are present within the 50% confidence contours of the location of the `GW superevent s191205ah <https://gracedb.ligo.org/superevents/S191205ah/>`_.
The code will open the following 3 pages: 
`Polygon-1 <http://archive.eso.org/scienceportal/home?poly=348.0469,%2020.7425,%20351.5625,%2020.1056,%20355.0781,%2018.8395,%20358.5938,%2018.2100,%200.7031,%2015.0948,%202.1094,%2012.0247,%202.8125,%208.9893,%203.5156,%205.9792,%203.5156,%202.9855,%203.5156,%200.0000,%203.5156,%20-2.9855,%203.5156,%20-5.9792,%203.5156,%20-8.9893,%202.8125,%20-12.0247,%200.0000,%20-12.0247,%20356.4844,%20-14.4776,%20352.9687,%20-13.2481,%20349.4531,%20-11.4152,%20345.9375,%20-10.2000,%20343.1250,%20-7.7827,%20341.7187,%20-4.7802,%20341.0156,%20-1.7908,%20340.3125,%201.1938,%20340.3125,%204.1815,%20341.0156,%207.1808,%20341.7187,%2010.2000,%20343.1250,%2013.2481,%20344.1797,%2016.0242,%20346.2891,%2019.1551&amp;sort=-obs_date>`_,
`Polygon-2 <http://archive.eso.org/scienceportal/home?poly=167.8711,%20-2.5374,%20166.9922,%20-2.0894,%20166.1133,%20-1.3430,%20165.9375,%20-0.5968,%20165.7617,%200.1492,%20165.5859,%200.8953,%20165.4102,%201.6415,%20165.2344,%202.3880,%20165.7617,%203.1349,%20166.2891,%203.8824,%20166.4648,%204.6305,%20166.6406,%205.3794,%20167.5195,%205.8292,%20168.3984,%205.6792,%20169.2773,%204.9300,%20169.4531,%204.1815,%20169.6289,%203.4338,%20169.8047,%202.6867,%20169.6289,%201.9401,%20169.4531,%201.1938,%20169.2773,%200.4476,%20169.1016,%20-0.2984,%20168.9258,%20-1.0445,%20168.7500,%20-1.7908,%20167.8711,%20-2.5374&amp;sort=-obs_date>`_, and
`Polygon-3 <http://archive.eso.org/scienceportal/home?poly=325.1953,%20-5.6792,%20325.3711,%20-4.9300,%20326.2500,%20-4.1815,%20327.1289,%20-3.7328,%20328.0078,%20-3.8824,%20328.8867,%20-4.6305,%20329.0625,%20-5.3794,%20329.2383,%20-6.1292,%20329.4141,%20-6.8801,%20329.2383,%20-7.6322,%20329.0625,%20-8.3856,%20328.8867,%20-9.1404,%20328.7109,%20-9.8969,%20327.8320,%20-10.6551,%20327.6562,%20-11.4152,%20326.7773,%20-11.5675,%20325.8984,%20-11.7198,%20325.3711,%20-10.9589,%20324.8437,%20-10.2000,%20325.0195,%20-9.4428,%20325.1953,%20-8.6873,%20325.0195,%20-7.9334,%20324.8437,%20-7.1808,%20325.0195,%20-6.4294,%20325.1953,%20-5.6792&amp;sort=-obs_date>`_.

The `ESOAsg` package also provides the `get_data_from_gw_event` script to fully exploit this functionality (see :ref:`archive-access-scripts` for futher details).
