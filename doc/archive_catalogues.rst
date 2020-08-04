==================
Archive Catalogues
==================

The ESO catalogue facility provides access to the collection of data that were produced by PIs of ESO programmes and then integrated into the ESO science archive through the Phase 3 process.
The full list of available catalogue could be found `here <https://www.eso.org/qi/>`_.

To access the data you can use the `programmatic access <http://archive.eso.org/programmatic/#TAP>`_ via the `tap_cat` TAP Service. 
The module `archive_catalogues` provides some simple `python` wrapper around this.

More examples on how to use the TAP services for the ESO archive are provided `in these notebooks <http://archive.eso.org/programmatic/HOWTO/>`_.

Please, when you use data from the archive, follow the :doc:`datapolicy`.

Overview
========

The module `archive_catalogue` is based on the `query.ESOCatalogues` class (a child of the `query.Query` class) that has the following attributes:

* `query` -- A string that contains the query to be perfomed via the TAP Service `tap_cat`
* `result_from_query` -- A table containing the result of the query
* `maxrec` -- An integer that the define the maximum number of records that will be returned for a query
* `type_of_query` -- A string that defines if the query will be run `synchronously` or `asynchronously`

After defining a `query` the `result_from_query` attribute is automatically filled by the method `run_query()`, for instance:
::

    from ESOAsg.queries import query_catalogues
    # define a query to obtain all `table_name` in the ESO Archive
    query = 'SELECT schema_name, table_name from TAP_SCHEMA.tables'
    # instantiate the class
    catalogue_list = query_catalogues.ESOCatalogues(query=query)
    # run the query
    catalogue_list.run_query()
    # print the result on terminal
    catalogue_list.result_from_query.pprint()

The module `archive_catalogue` provides a set of additional quality-of-life functions to circumvent the actual creation of the queries.
These allow one to directly obtain the information needed in a `table` format.

Some examples
=============

.. note::
   The way the queries are created allows one to set as input either `collections` or `tables`.
   However, we strongly discurage to use both at the same times.
   Given that the connector between the two conditions is an `AND` this may give rise to an un-expected behaviour

Which catalogues are available
------------------------------

In general, it is possible to explore which catalogue are availabe in the ESO archive either with the `query interface webpage <https://www.eso.org/qi/>`_, or with the `archive science portal <https://archive.eso.org/scienceportal/home?data_release_date=*:2020-07-23&dp_type=CATALOG&sort=-obs_date&s=P%2fDSS2%2fcolor&f=134.496111&fc=-1,-1&cs=J2000&av=true&ac=false&c=9,10,11,12,13,14,15,16,17,18,19,20&ta=RES&dts=true&sdtm=%7b%22CATALOG%22%3atrue%7d&at=0,0&sr=i>`_, or by running the `TAP query to obtain all versions of all catalogues <http://archive.eso.org/programmatic/#TAP?e=1&f=text&m=200&q=SELECT%20cat_id%2C%20collection%2C%20table_name%2C%20title%2C%20number_rows%2C%20number_columns%2C%20version%2C%20acknowledgment%20FROM%20TAP_SCHEMA.tables%20WHERE%20schema_name%20%3D%20'safcat'%0A&>`_:
::

    SELECT
        cat_id, collection, table_name, title, number_rows, number_columns, 
        version, acknowledgment
    FROM
        TAP_SCHEMA.tables 
    WHERE 
        schema_name = 'safcat'

or this one to obtain `only the latest version of the catalogues <http://archive.eso.org/programmatic/#TAP?e=1&f=text&m=200&q=SELECT%20t1.cat_id%2C%20t1.collection%2C%20t1.table_name%2C%20t1.title%2C%20t1.number_rows%2C%20t1.number_columns%2C%20t1.version%2C%20t1.acknowledgment%20FROM%20tables%20t1%20LEFT%20OUTER%20JOIN%20tables%20t2%20ON%20(t1.title%20%3D%20t2.title%20AND%20t1.version%20%3C%20t2.version)%20WHERE%20t2.title%20IS%20null%20AND%20t1.cat_id%20IS%20NOT%20null%20AND%20t1.schema_name%20%3D%20'safcat'%0A&>`_:
::

    SELECT
        t1.cat_id, t1.collection, t1.table_name, t1.title,
        t1.number_rows, t1.number_columns,
        t1.version, t1.acknowledgment
    FROM
        tables t1
    LEFT OUTER JOIN 
        tables t2 ON (t1.title = t2.title AND t1.version < t2.version)
    WHERE
       t2.title IS null AND t1.cat_id IS NOT null AND t1.schema_name = 'safcat'

Alternatively, it is possible to obtain an `astropy.table` containing information on all catalogues and all their versions using (note that the query is more complicated of the ones above because more information are collected):
::

    from ESOAsg import archive_catalogues
    archive_catalogues.all_catalogues(all_versions=True)

This returns an `astropy.table` containing:

+------------------+-----------------------------------------------------------------------------------------+
| Column name      | Description                                                                             |
+==================+=========================================================================================+
| collection       | Name of the Phase 3 collection the catalog belongs to                                   |
+------------------+-----------------------------------------------------------------------------------------+
| title            | Title of the catalog                                                                    |
+------------------+-----------------------------------------------------------------------------------------+
| version          | Version of the catalog                                                                  |
+------------------+-----------------------------------------------------------------------------------------+
| table_name       | The fully qualified table name                                                          |
+------------------+-----------------------------------------------------------------------------------------+
| filter           | Name(s) of the filter bandpasses the original data were aquired with                    |
+------------------+-----------------------------------------------------------------------------------------+
| instrument       | Name(s) of the instrument(s) the original data were acquired with                       |
+------------------+-----------------------------------------------------------------------------------------+
| telescope        | Name(s) of the telescope(s) the original data were acquired with                        |
+------------------+-----------------------------------------------------------------------------------------+
| publication_date | The data the catalog was published                                                      |
+------------------+-----------------------------------------------------------------------------------------+
| description      | Describes tables in the tableset                                                        |
+------------------+-----------------------------------------------------------------------------------------+
| number_rows      | Number of rows present in this version of the catalog                                   |
+------------------+-----------------------------------------------------------------------------------------+
| number_columns   | Number of columns present in this version of the catalog                                |
+------------------+-----------------------------------------------------------------------------------------+
| rel_descr_url    | Location of the release description document (typically a pdf)                          |
+------------------+-----------------------------------------------------------------------------------------+
| acknowledgment   | It provides the sentence to be used in your publication when making use of this catalog |
+------------------+-----------------------------------------------------------------------------------------+
| cat_id           | Internal catalog identifier,                                                            |
+------------------+-----------------------------------------------------------------------------------------+
| mjd_obs          | The observational data this catalog is based were taken between mjd_obs and mjd_end     |
+------------------+-----------------------------------------------------------------------------------------+
| mjd_end          | The observational data this catalog is based were taken between mjd_obs and mjd_end     |
+------------------+-----------------------------------------------------------------------------------------+
| skysqdeg         | Area of the sky (in square degrees) covered by this catalog                             |
+------------------+-----------------------------------------------------------------------------------------+
| bibliography     | Bibliographic reference in the form of either a BIBCODE or a DOI                        |
+------------------+-----------------------------------------------------------------------------------------+
| document_id      | Internal identifier of the release description document [#foot_cat]_                    |
+------------------+-----------------------------------------------------------------------------------------+
| from_column      | Column in the from_table                                                                |
+------------------+-----------------------------------------------------------------------------------------+
| target_table     | The table with the primary key                                                          |
+------------------+-----------------------------------------------------------------------------------------+
| target_column    | Column in the target_table                                                              |
+------------------+-----------------------------------------------------------------------------------------+
| last_version     | True if this is the latest version of the catalog                                       |
+------------------+-----------------------------------------------------------------------------------------+
| RA_id            | Identifier for RA in the catalog                                                        |
+------------------+-----------------------------------------------------------------------------------------+
| Dec_id           | Identifier for Dec in the catalog                                                       |
+------------------+-----------------------------------------------------------------------------------------+

.. note::
   At first sight it may seem that not all catalogs have the `RA_id` and `Dec_id`.
   This is because the catalogue is spreaded into more than one table.
   To identify the same source among the differnt tables of a catalogue the `target_table` and `target_column` should be used.

Which columns are in a catalog
------------------------------

It is possible to get information on all columns present in a catalogue by running the following `TAP query <http://archive.eso.org/programmatic/#TAP?e=1&f=text&m=200&q=SELECT%20table_name%2C%20column_name%2C%20ucd%2C%20datatype%2C%20description%2C%20unit%0AFROM%20TAP_SCHEMA.columns%0AWHERE%20table_name%20%3D%20'viking_er5_zyjj_1j_2hks_catMetaData_fits_V4'%0A&>`_ for the `VIKING DR4 <https://www.eso.org/rm/api/v1/public/releaseDescriptions/135>`_ catalogue:
::

    SELECT 
        table_name, column_name, ucd, datatype, description, unit
    FROM 
        TAP_SCHEMA.columns
    WHERE 
        table_name = 'viking_er5_zyjj_1j_2hks_catMetaData_fits_V4'

A similar result can be obtained running:
::

    archive_catalogues.columns_info(tables='viking_er5_zyjj_1j_2hks_catMetaData_fits_V4')

where the result is stored in an `astropy.table`.


.. rubric:: Footnotes

.. [#foot_cat] The web user interface for this catalog is reachable via the URL computed appending the `cat_id` to the string: https://www.eso.org/qi/catalogQuery/index/
