import numpy as np
from astropy.table import MaskedColumn

from ESOAsg import msgs
from ESOAsg import default
from ESOAsg.ancillary import checks
from ESOAsg.core import tap_queries
from ESOAsg.queries import query_catalogues


def all_catalogues(verbose=False, all_versions=False):
    r"""Load a list with all ESO catalogues

    For further information check `https://www.eso.org/qi/`

    Args:
        verbose (`bool`):
            if set to `True` additional info will be displayed
        all_versions (`bool`):
            if set to `True` also obsolete versions of the catalogues are listed. If `False` only catalogues with the
            `last_version` = `True` are loaded.

    Returns:
        all_catalogues_table (`astropy.table`):
            `astropy.table` containing: `cat_id`, `collection`, `table_name`, `title`, `number_rows`, `version`,
            `acknowledgment` of all catalogues currently present at ESO. In addition the column `last_version` is added.
            This is an attempt to remove obsolete catalogues based on the version number and the title of the catalogue.
    """
    query_all_catalogues = query_catalogues.ESOCatalogues(query=tap_queries.create_query_all_catalogues(
        all_versions=all_versions))
    # Obtaining query results
    query_all_catalogues.run_query()
    # Sorting
    query_all_catalogues.result_from_query.sort(['collection', 'table_name', 'version'])
    # Checking for obsolete
    unique_titles = np.unique(query_all_catalogues.result_from_query['title'].data).tolist()
    last_version = np.zeros_like(query_all_catalogues.result_from_query['version'].data, dtype=bool)
    for unique_title in unique_titles:
        most_recent_version = np.nanmax(
            query_all_catalogues.result_from_query['version'].data[(query_all_catalogues.result_from_query[
                                                                                  'title'].data == unique_title)])
        last_version[(query_all_catalogues.result_from_query['title'].data == unique_title) &
                     (query_all_catalogues.result_from_query['version'].data == most_recent_version)] = True
    query_all_catalogues.result_from_query.add_column(MaskedColumn(data=last_version, name='last_version', dtype=bool,
                                                 description='True if this is the latest version of the catalog'))
    query_all_catalogues.result_from_query['table_name'].data.data[:] = checks.from_bytes_to_string(query_all_catalogues.result_from_query[
                                                                                      'table_name'].data.data)
    return query_all_catalogues.result_from_query



