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



Some examples
=============

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

Alternatively, it is possible to obtain an `astropy.table` containing information on all catalogues (and all their versions) using:
::

    from ESOAsg import archive_catalogues
    archive_catalogues.all_catalogues(all_versions=True)

