.. _archive-observations:

====================
Archive Observations
====================

The `ESO archive <http://archive.eso.org/cms.html>`_ currently (June 2020) contains more than `1.7 million spectra <https://archive.eso.org/scienceportal/home?data_release_date=*:2020-06-24&dp_type=SPECTRUM&sort=-obs_date&s=P%2fDSS2%2fcolor&f=177.115919&fc=-1,-1&cs=J2000&av=true&ac=false&c=9,10,11,12,13,14,15,16,17,18,19,20&ta=RES&dts=true&sdtm=%7b%22SPECTRUM%22%3atrue%7d&at=119.452774,-60.30286&sr=i>`_, more than `650,000 images <https://archive.eso.org/scienceportal/home?data_release_date=*:2020-06-24&dp_type=IMAGE&sort=-obs_date&s=P%2fDSS2%2fcolor&f=177.115919&fc=-1,-1&cs=J2000&av=true&ac=false&c=9,10,11,12,13,14,15,16,17,18,19,20&ta=RES&dts=true&sdtm=%7b%22IMAGE%22%3atrue%7d&at=160.465004,19.501825&sr=i>`_, and more than `240,000 cubes <https://archive.eso.org/scienceportal/home?data_release_date=*:2020-06-24&dp_type=CUBE&sort=-obs_date&s=P%2fDSS2%2fcolor&f=177.115919&fc=-1,-1&cs=J2000&av=true&ac=false&c=9,10,11,12,13,14,15,16,17,18,19,20&ta=RES&dts=true&sdtm=%7b%22CUBE%22%3atrue%7d&at=239.591811,-14.166308&sr=i>`_.

There are three main ways to access the vaste amount of information present in the `ESO archive <http://archive.eso.org/cms.html>`_:

* the `Raw Data query form <http://archive.eso.org/eso/eso_archive_main.html>`_
* the `Science Portal <http://archive.eso.org/scienceportal/home>`_ to browse and access the processed data (see also ref::`archive-science-portal`)
* the `Programmatic and Tools access <http://archive.eso.org/programmatic/>`_ which permits direct database access to both raw and processed data, and to the ambient condition measurements

In addition, the `archive_observations` provides simple wrappers to efficiently embed the access to the ESO archive into `python` routines.

Overview
========

Similarly to the module `archive_catalogues` (see :ref:`archive-catalogues`), the base for `archive_observations` is the `query_observations.ESOObservations` class (a child of the `query.Query` class) that has the following attributes:

* `query` -- A string that contains the query to be perfomed via the TAP Service `tap_obs`
* `result_from_query` -- A table containing the result of the query
* `maxrec` -- An integer that the define the maximum number of records that will be returned for a query
* `type_of_query` -- A string that defines if the query will be run `synchronously` or `asynchronously`

After defining a `query` you can fill the `result_from_query` attribute using the method `run_query()`, for instance:
::

    from ESOAsg.queries import query_observations
    # define a query to obtain all ALMA observations in the ESO Archive
    query = "SELECT target_name, dp_id, s_ra, s_dec, instrument_name FROM ivoa.ObsCore WHERE instrument_name = 'ALMA'"
    # instantiate the class adn limit to the first 100 results
    alma_query = query_observations.ESOObservations(query=query, maxrec=100)
    # run the query
    alma_query.run_query()
    # print the result on terminal
    alma_query.result_from_query.pprint()

The module `archive_observations` provides a set of functions to help a user to search and download for data without explicitly creating a TAP query.


Some examples
=============

query to a specific position in the sky
---------------------------------------

The access to data in the ESO archive that intersects a cone with a certain radius (for instance around the star HD 057060: RA=07:18:40.38, Dec.=-24:33:31.3, radius=5arcsec) it is possible either with the `archive science portal <https://archive.eso.org/scienceportal/home?pos=109.66825,-24.5587&r=0.00138888888>`_ (see also :ref:`archive-science-portal`), or by running the `TAP query to ivoa.ObsCore <http://archive.eso.org/programmatic/#TAP?f=text&m=200&q=SELECT%0A%20%20%20%20%20%20%20%20target_name%2C%20dp_id%2C%20s_ra%2C%20s_dec%2C%20t_exptime%2C%20em_min%2C%20em_max%2C%0A%20%20%20%20%20%20%20%20dataproduct_type%2C%20instrument_name%2C%20obstech%2C%20abmaglim%2C%0A%20%20%20%20%20%20%20%20proposal_id%2C%20obs_collection%0AFROM%0A%20%20%20%20%20%20%20%20ivoa.ObsCore%0AWHERE%0A%20%20%20%20%20%20%20%20INTERSECTS(CIRCLE('ICRS'%2C109.66824871%2C-24.55869773%2C5.%2F3600.)%2C%20s_region)%20%3D%201&>`_:
::

    SELECT
        target_name, dp_id, s_ra, s_dec, t_exptime, em_min, em_max,
        dataproduct_type, instrument_name, obstech, abmaglim,
        proposal_id, obs_collection
    FROM
        ivoa.ObsCore
    WHERE
        INTERSECTS(CIRCLE('ICRS',109.66824871,-24.55869773,0.0014), s_region) = 1

The same result, could be obtained using the function `query_from_radec()`:
::

    from ESOAsg import archive_observations
    from astropy.coordinates import SkyCoord
    star_position = SkyCoord.from_name('HD 057060')
    eso_data = archive_observations.query_from_radec(star_position, radius=5.)

`eso_data` is an `astropy.table` containing the information on the available data products.
It is possible to download all this data feeding the `dp_id` column to the `download` function:
::

    archive_observations.download(eso_data['dp_id']

A step-by-step explanation for the procedure to select and download data is provided in `this jupyter notebook <http://localhost:8889/notebooks/HOWTO_getDataFromRaDec.ipynb>`_.

The `ESOAsg` package also provides the `get_data_from_radec` script to fully exploit this functionality (see :ref:`archive-access-scripts` for futher details).
