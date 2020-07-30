r"""Child class to run queries to the ESO catalogue
"""

import numpy as np
from astropy.table import MaskedColumn

from ESOAsg import default
from ESOAsg import msgs
from ESOAsg.core import tap_queries
from ESOAsg.queries import query


class ESOCatalogues(query.Query):
    r"""This class is designed to query the ESO archive for raw, reduced, and ambient data.

    """
    def __init__(self, query=None, result_from_query=None, maxrec=default.get_value('maxrec')):
        # assign tap_service
        super().__init__(tap_service=tap_queries.define_tap_service('eso_tap_cat'), query=query,
                         result_from_query=result_from_query, maxrec=maxrec)

    def set_last_version(self, update=True):
        r"""Set the `last_version` column to the `result_from_query` attribute

        `last_version` is a column of `bool` where `False` means that there is a more update version of a catalogue

        This works only if `result_from_query` contains the columns: `version` and `title`. In case the `last_version`
        column is already present, a warning is raised.

        Args:
            update (bool): in case the `last_version` column is already present the code will update the value only
                if `update` is set to `True`

        Returns:
            None

        """
        # Require that title and version are present in result_from_query
        for check_column in ['title', 'version']:
            if check_column not in self.which_columns():
                msgs.warning('{} column not present, `last_version` will not be created'.format(check_column))

        # Check that last_version is not present in result_from_query
        if 'last_version' in self.which_columns():
            if update:
                msgs.warning('`last_version` column already present')
            else:
                msgs.warning('`last_version` column already present and it will not be updated')
                return

        # Get last versions
        unique_titles = np.unique(self.result_from_query['title'].data).tolist()
        last_version = np.zeros_like(self.result_from_query['version'].data, dtype=bool)
        for unique_title in unique_titles:
            most_recent_version = np.nanmax(self.result_from_query['version'].data[
                                                (self.result_from_query['title'].data == unique_title)])
            last_version[(self.result_from_query['title'].data == unique_title) &
                         (self.result_from_query['version'].data == most_recent_version)] = True
        self.result_from_query.add_column(MaskedColumn(data=last_version, name='last_version', dtype=bool,
                                                       description='True if this is the latest version of the catalog'))
        return

