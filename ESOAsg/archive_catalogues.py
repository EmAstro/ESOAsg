import numpy as np
from astropy.table import MaskedColumn

from ESOAsg import msgs
from ESOAsg.ancillary import checks
from ESOAsg.core import tap_queries
from ESOAsg.queries import query_catalogues


def all_catalogues_info(all_versions=False, verbose=False):
    r"""Load an `astropy.table` with information on all catalogues present in the ESO archive

    The output table will contain the following columns:
    `collection`, `title`, `version`, `table_name`, `filter`, `instrument`, `telescope`, `publication_date`,
    `description`, `number_rows`, `number_columns`, `rel_descr_url`, `acknowledgment`, `cat_id`, `mjd_obs`,
    `mjd_end`, `skysqdeg`, `bibliography`, `document_id`, `from_column`, `target_table`, `target_column`,
    `last_version`

    For further information check the `ESO catalogue facility <https://www.eso.org/qi/>`_

    .. note::
        This is analogue to run:

        >>> catalogues_info(collections=None, tables=None)

        with the difference that, given that no constraints are set, this returns the master table with all catalogues
        present in the ESO archive

    Args:
        verbose (bool): if set to `True` additional info will be displayed
        all_versions (bool): if set to `True` also obsolete versions of the catalogues are listed

    Returns:
        astropy.table: table containing the information on all catalogues present in the ESO archive

    """
    return catalogues_info(all_versions=all_versions, collections=None, tables=None, verbose=verbose)


def catalogues_info(all_versions=False, collections=None, tables=None, verbose=False):
    r"""Load an `astropy.table` with information on the selected catalogues

    Specific catalogues can be selected by selecting a list of collections or a list of table_names. If both
    `collections` and `tables` are set to `None` information on all ESO catalogues will be returned. For further
    information check the `ESO catalogue facility <https://www.eso.org/qi/>`_

    The output table will contain the following columns:
    `collection`, `title`, `version`, `table_name`, `filter`, `instrument`, `telescope`, `publication_date`,
    `description`, `number_rows`, `number_columns`, `rel_descr_url`, `acknowledgment`, `cat_id`, `mjd_obs`,
    `mjd_end`, `skysqdeg`, `bibliography`, `document_id`, `from_column`, `target_table`, `target_column`,
    `last_version`

    .. note::
       The way the query is created is to set as input or `collections` or `tables`. Particular attention should be
       given if both `collections` and `tables` are not `None`. Given that the connector between the two conditions is
       an `AND` this may give rise to an un-expected behaviour

    Args:
        all_versions (bool): if set to `True` also obsolete versions of the catalogues are listed
        collections (any, optional): list of `str` containing the names of the collections (or single `str`) for which
            the query will be limited
        tables (any, optional): list of `str` containing the table_name of the tables (or single `str`) for which the
            query will be limited
        verbose (bool): if set to `True` additional info will be displayed

    Returns:
        astropy.table: table containing the information on the selected catalogues

    """
    # Test on collections
    clean_collections = _is_collection_list_at_eso(collections)
    # Test for tables
    clean_tables = _is_table_list_at_eso(tables)
    # instantiate ESOCatalogues
    query_for_catalogues = query_catalogues.ESOCatalogues(query=tap_queries.create_query_all_catalogues(
        all_versions=all_versions, collections=clean_collections, tables=clean_tables))
    # rise warning
    if (collections is not None) and (tables is not None):
        msgs.warning('Setting conditions for both `collections` and `tables`. Please check that this is the wanted '
                     'behaviour:')
    # Print query
    if verbose:
        query_for_catalogues.print_query()
    # Obtaining query results
    query_for_catalogues.run_query(to_string=True)
    # Sorting
    query_for_catalogues.result_from_query.sort(['collection', 'table_name', 'version'])
    # Checking for obsolete and creating last_version column
    # this is redundant in case `all_version` = Tru
    query_for_catalogues.set_last_version(update=True)
    catalogues_table = query_for_catalogues.get_result_from_query()
    # Get info on source ID RA and DEC for all collections
    id_ra_dec_table = _get_id_ra_dec_from_columns(collections=clean_collections)
    ra_id, dec_id = [], []
    for t_name in catalogues_table.iterrows('table_name'):
        ra_id_table = id_ra_dec_table[(id_ra_dec_table['table_name'] == t_name) &
                                      (id_ra_dec_table['ucd'] == 'pos.eq.ra;meta.main')]['column_name'].data.data
        dec_id_table = id_ra_dec_table[(id_ra_dec_table['table_name'] == t_name) &
                                       (id_ra_dec_table['ucd'] == 'pos.eq.dec;meta.main')]['column_name'].data.data
        if len(ra_id_table) == 1 and len(dec_id_table) == 1:
            ra_id.append(ra_id_table[0]), dec_id.append(dec_id_table[0])
        elif len(ra_id_table) > 1 and len(dec_id_table) > 1:
            msgs.warning('Impossible to identify RA and Dec columns in table: {}'.format(t_name))
        else:
            ra_id.append(None), dec_id.append(None)
    catalogues_table.add_column(MaskedColumn(data=np.asarray(ra_id), name='RA_id', dtype=str,
                                             description='Identifier for RA in the catalog'))
    catalogues_table.add_column(MaskedColumn(data=np.asarray(dec_id), name='Dec_id', dtype=str,
                                             description='Identifier for Dec in the catalog'))
    return catalogues_table


def columns_info(collections=None, tables=None, verbose=False):
    r"""Load a query that get names (and corresponding ucd) of the columns present in a collection

    If `collections` and `tables` are `None` the query for the column of all collections and tables in the ESO
    archive is returned.

    .. note::
       The way the query is created is to set as input or `collections` or `tables`. Particular attention should be
       given if both `collections` and `tables` are not `None`. Given that the connector between the two conditions is
       an `AND` this may give rise to an un-expected behaviour

    Args:
        collections (any, optional): list of `str` containing the names of the collections (or single `str`) from which
            the columns will be extracted
        tables (any, optional): list of `str` containing the names of the tables (or single `str`) from which the
            columns will be extracted
        verbose (bool): if set to `True` additional info will be displayed

    Returns:
        astropy.table: table of all columns present in a table/collection. Information are stored in `table_name`,
            `column_name`, `ucd`, `datatype`, `description`, and `unit`

    """
    # test on collections
    clean_collections = _is_collection_list_at_eso(collections)
    # test on tables
    clean_tables = _is_table_list_at_eso(tables)
    # instantiate ESOCatalogues
    query_all_columns_info = query_catalogues.ESOCatalogues(query=tap_queries.create_query_all_columns(
        collections=clean_collections, tables=clean_tables))
    # rise warning
    if (collections is not None) and (tables is not None):
        msgs.warning('Setting conditions for both `collections` and `tables`. Please check that this is the wanted '
                     'behaviour:')
    # Print query
    if verbose or ((collections is not None) and (tables is not None)):
        query_all_columns_info.print_query()
    # Obtaining query results
    query_all_columns_info.run_query(to_string=True)
    all_columns_table = query_all_columns_info.get_result_from_query()
    return all_columns_table


def get_catalogues(collections=None, tables=None, columns=None, type_of_query='sync', all_versions=False,
                   maxrec=None, verbose=False):
    r"""Query the ESO tap_cat service for specific catalogues

    There are two ways to select the catalogues you are interested in. Either you select directly the table_name (or the
    list of table_names) that you want to query, or you select a collection (or a list of collections). If you select
    this latter option, what happens in the background is that the code is going to search for the table(s)
    corresponding to the given collection and query them

    If you are asking for more than one table, the result will be listed in a list of `astropy.tables` one per
    table

    Args:
        collections (any): list of `str` containing the names (or a single `str`) of the collections for
            which the query will be limited
        tables (any): list of `str` containing the table_name of the tables for which the query will be limited
        columns (any): list of the `column_name` that you want to download. The full list of the columns in a
            table can be found by running `columns_info()`
        all_versions (bool): if set to `True` also obsolete versions of the catalogues are searched in case
            `collections` is given
        type_of_query (str): type of query to be run
        maxrec (int, optional): define the maximum number of entries that a single query can return. If it is `None` the
            value is set by the limit of the service.
        verbose (bool): if set to `True` additional info will be displayed

    Returns:
        any: `astropy.table` or `list` of `astropy.tables` containing the queried catalogues

    """
    # ToDo EMA: conditions to select properties from catalogues should be added in this. Both with ANDs and with ORs
    # test on collections
    clean_collections = _is_collection_list_at_eso(collections)
    # test on tables
    clean_tables = _is_table_list_at_eso(tables)
    # the query works on
    if clean_tables is None:
        clean_tables = []
    if clean_collections is not None:
        for clean_collection in clean_collections:
            clean_tables.append(_get_tables_from_collection(clean_collection, all_versions=all_versions))
    # This removes possible duplicates and removes None
    clean_tables = list(filter(None, list(set(clean_tables))))
    # if max_rec is set to None, the entire length of the catalogue is considered:
    if maxrec is None:
        maxrec_list = []
        for table in clean_tables:
            maxrec_list.append(_get_catalogue_length_from_table(table, all_versions=all_versions))
    else:
        maxrec_list = [maxrec] * len(clean_tables)

    list_of_catalogues = []
    for table_name, maxrec_for_table in zip(clean_tables, maxrec_list):
        # test for columns
        columns_in_table = _is_column_list_in_catalogues(columns, tables=table_name)
        # instantiate ESOCatalogues
        query_table = query_catalogues.ESOCatalogues(query=tap_queries.create_query_table(table_name,
                                                                                          columns=columns_in_table),
                                                     type_of_query=type_of_query, maxrec=maxrec_for_table)
        if verbose:
            query_table.print_query()
        query_table.run_query(to_string=True)
        catalogue = query_table.get_result_from_query()
        list_of_catalogues.append(catalogue)
        msgs.info('The query to {} returned {} entries (with a limit set to maxrec={})'.format(table_name,
                                                                                               len(catalogue),
                                                                                               maxrec_for_table))
    if len(list_of_catalogues) == 0:
        return None
    elif len(list_of_catalogues) == 1:
        return list_of_catalogues[0]
    else:
        return list_of_catalogues


def _get_id_ra_dec_from_columns(collections=None):
    r"""Returns the column names corresponding to source ID, RA, and DEC from a list of collections

    This is base on the tokens:
    * `meta.id;meta.main` -> source ID
    * `pos.eq.ra;meta.main` -> RA
    * `pos.eq.dec;meta.main` -> Dec

    Args:
        collections (any): list of `str` (or a single `str`) containing the names of the collections from which the
            name of the source ID, RA, and Dec will be extracted

    Returns:
        astropy.table: table containing the column names corresponding to source ID, RA, and Dec

    """
    all_columns_table = columns_info(collections)
    filter_for_tokens = ((all_columns_table['ucd'].data == 'meta.id;meta.main') | (
                          all_columns_table['ucd'].data == 'pos.eq.ra;meta.main') | (
                          all_columns_table['ucd'].data == 'pos.eq.dec;meta.main'))
    return all_columns_table[filter_for_tokens]


def _get_tables_from_collection(collection, all_versions=False):
    r"""Returns the table_name corresponding to a given collection

    Args:
        collection (str): name of the collection for which the tables will be extracted
        all_versions (bool): if set to `True` also obsolete versions of the catalogues are listed

    Returns:
        list: list containing the `table_name` corresponding to the selected `collection`

    """
    if not _is_collection_at_eso(collection):
        return None
    table_all_catalogues = all_catalogues_info(verbose=False, all_versions=all_versions)
    table_selected_catalogues = table_all_catalogues[(table_all_catalogues['collection'].data == collection)]
    list_selected_tables = table_selected_catalogues['table_name'].data.data.to_list()
    return list_selected_tables


def _get_catalogue_length_from_table(table_name, all_versions=False):
    r"""Returns the length of a catalogue given a `table_name`

    Args:
        table_name (str): name of the collection for which the tables will be extracted
        all_versions (bool): if set to `True` also obsolete versions of the catalogues are listed

    Returns:
        int: `number_rows` corresponding to the selected `table_name`

    """
    if not _is_table_at_eso(table_name):
        return None
    table_all_catalogues = all_catalogues_info(verbose=False, all_versions=all_versions)
    table_selected_catalogues = table_all_catalogues[(table_all_catalogues['table_name'].data == table_name)]
    selected_number_columns = int(table_selected_catalogues['number_rows'].data.data)
    return selected_number_columns


def _is_collection_list_at_eso(collections):
    r"""Check if a given list of collections is present in the ESO archive

    Args:
        collections (any): str or list of collections to be tested

    Returns:
        list: same of `collections` but collection not present at ESO are removed

    """
    assert collections is None or isinstance(collections, (str, list)), r'`collections` must be `None`, ' \
                                                                        r'or a `str` or a `list`'
    collections_list = checks.from_string_to_list(collections)
    if collections_list is not None:
        clean_collections = []
        for collection in collections_list:
            if _is_collection_at_eso(collection):
                clean_collections.append(collection)
    else:
        clean_collections = None
    return clean_collections


def _is_collection_at_eso(collection):
    r"""Check if a given collection is present in the ESO archive

    Args:
        collection (str): collection to be tested

    Returns:
        bool: `True` if the table is present in the ESO archive, `False` and warning raised otherwise

    """
    is_at_eso = True
    table_all_catalogues = all_catalogues_info(verbose=False, all_versions=False)
    all_collections_list = np.unique(table_all_catalogues['collection'].data.data).tolist()
    if collection not in all_collections_list:
        msgs.warning('Collection: {} not recognized. Possible values are:\n{}'.format(collection, all_collections_list))
        is_at_eso = False
    return is_at_eso


def _is_table_list_at_eso(tables):
    r"""Check if a given list of table_names is present in the ESO archive

    Args:
        tables (any): `list` of table_names (or a single `str`) to be tested

    Returns:
        list: same of `tables` but tables not present at ESO are removed

    """
    assert tables is None or isinstance(tables, (str, list)), r'`tables` must be `None` or a `str` or a `list`'
    tables_list = checks.from_string_to_list(tables)
    if tables is not None:
        clean_tables = []
        for table in tables_list:
            if _is_table_at_eso(table):
                clean_tables.append(table)
    else:
        clean_tables = None
    return clean_tables


def _is_table_at_eso(table_name):
    r"""Check if a given table is present at ESO

    Args:
        table_name (str): table to be tested.

    Returns:
        bool: `True` if the table is present in the ESO archive. `False` and warning raised otherwise.

    """
    is_at_eso = True
    table_all_catalogues = all_catalogues_info(verbose=False, all_versions=True)
    all_table_list = table_all_catalogues['table_name'].data.data.tolist()
    last_version_list = table_all_catalogues['last_version'].data.data.tolist()
    if table_name not in all_table_list:
        msgs.warning('Table: {} not recognized. Possible values are:\n{}'.format(table_name, table_all_catalogues))
        is_at_eso = False
    else:
        if not last_version_list[all_table_list.index(table_name)]:
            msgs.warning('{} is not the most recent version of the queried catalogue'.format(table_name))
    return is_at_eso


def _is_column_list_in_catalogues(columns, collections=None, tables=None):
    r"""Check if a given list of columns is present in the ESO archive

    It is possible to test for given collection (or table) by setting the appropriate values as input

    Args:
        columns (any): list of string containing the column_name (or the single `str`) to be tested
        collections (any): list of `str` containing the names of the collections (or a single `str`) from which the
            columns will be extracted
        tables (any): list of `str`  (or a single `str`) containing the names of the tables from which the columns
            will be extracted

    Returns:
        list: same of `columns` but columns not present at in the collections/tables are removed

    """
    assert columns is None or isinstance(columns, (str, list)), r'`columns` must be `None` or a `str` or a `list`'
    columns_list = checks.from_string_to_list(columns)
    if columns is not None:
        # test if it is a valid column
        clean_columns = []
        for column in columns_list:
            if _is_column_in_catalogues(column, collections=collections, tables=tables):
                clean_columns.append(column)
    else:
        clean_columns = None
    return clean_columns


def _is_column_in_catalogues(column_name, collections=None, tables=None):
    r"""Check if a given column is present in the ESO archive

    Args:
        column_name (str): column to be tested
        collections (any): list of `str` containing the names of the collections (or a single `str`) from which the
            columns will be extracted
        tables (any): list of `str`  (or a single `str`) containing the names of the tables from which the columns
            will be extracted

    Returns:
        bool: `True` if the column is present in the selected collections/tables. `False` and warning raised
            otherwise

    """
    is_at_eso = True
    table_all_columns = columns_info(collections=collections, tables=tables, verbose=False)
    all_column_list = table_all_columns['column_name'].data.data.tolist()
    if column_name not in all_column_list:
        msgs.warning('Column: {} not recognized. Possible values are:\n{}'.format(column_name, all_column_list))
        is_at_eso = False
    return is_at_eso


'''
def query_from_radec(table_names, position, radius=None, maxrec=default.get_value('maxrec')):
    r"""Query the ESO tap_cat service for a specific a specific location.

    The `position` needs to be given as an `astropy.coordinates.SkyCoord` object.

    Args:
        table_names (`list`):
            list of the tables to be queried. To check the full list of catalogues run `all_catalogues()`
        position (`astropy.coordinates.SkyCoord`):
            Coordinates of the sky you want to query in the format of an `astropy.coordinates.SkyCoord` object. Note
            that at the moment it works for one target at the time. For further detail see here:
            `astropy coordinates <https://docs.astropy.org/en/stable/coordinates/>`_
        radius (`float`):
            Search radius you want to query in arcseconds. Note that in case `None` is given, the query will be
            performed with the `CONTAINS(POINT('',RA,Dec), s_region)` clause instead of the
            `CONTAINS(s_region,CIRCLE('',RA,Dec,radius/3600.))` one. See here for further examples:
            `tap obs examples <http://archive.eso.org/tap_obs/examples>`_
        maxrec (`numpy.int`):
            Define the maximum number of file that a single query can return from the ESO archive. You probably never
            need this. By default is set by the `default.txt` file.

    Returns:
        result_from_query (`astropy.Table`):
            Result from the query to the TAP service
    """
'''