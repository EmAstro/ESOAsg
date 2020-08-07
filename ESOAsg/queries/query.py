from ESOAsg import msgs
from ESOAsg.ancillary import cleaning_lists
from ESOAsg.core import tap_queries

__all__ = ['Query']


class Query:
    r"""Base class that dictate the general behaviour of a query

    Attributes:
        tap_service (pyvo.dal.tap.TAPService): TAP service that will be used for the query
        query (str): String containing the query
        type_of_query (str): type of query to be run
        maxrec (int, optional): define the maximum number of entries that a single query can return
        result_from_query (astropy.table.Table): result from the query to the TAP service
    """

    def __init__(self, tap_service=None, query=None, type_of_query='sync', result_from_query=None, maxrec=None):
        self.tap_service = tap_service
        self.query = query
        self.result_from_query = result_from_query
        self.maxrec = maxrec
        if type_of_query not in tap_queries.TAP_QUERY_TYPES:
            msgs.warning('{} not a valid entry for the type of TAP query. Possibilities are: {}'.format(
                type_of_query, tap_queries.TAP_QUERY_TYPES))
            msgs.warning('The `type_of_query` attribute will be set to `sync`')
            self.type_of_query = 'sync'
        else:
            self.type_of_query = type_of_query

    def clean_query(self):
        r"""Set the `query` attribute to None
        """
        self.query = None

    def clean_result_from_query(self):
        r"""Set the `result_from_query` attribute to None
        """
        self.result_from_query = None

    def get_result_from_query(self):
        r"""Returns a copy of `result_from_query`
        """
        if self.result_from_query is None:
            copy_of_result_from_query = None
        else:
            copy_of_result_from_query = self.result_from_query.copy()
        return copy_of_result_from_query

    def run_query(self, to_string=True):
        r"""Run the query and store results in the `result_from_query` attribute

        Args:
            to_string (`bool`, optional): if set to True, if a column is in `bytes` format it transform it to `str`

        """
        self.result_from_query = tap_queries.run_query(self.tap_service, self.query, self.type_of_query,
                                                       maxrec=self.maxrec)
        if to_string:
            for column_id in self.which_columns():
                self.result_from_query[column_id].data.data[:] = cleaning_lists.from_bytes_to_string(
                    self.result_from_query[column_id].data.data)
        return

    def print_query(self):
        r"""Print the query on the terminal
        """
        tap_queries.print_query(self.query)
        return

    def which_service(self):
        r"""Return the `tap_service` that is used together with a description of it

        Returns:
            None
        """
        return tap_queries.which_service(self.tap_service)

    def which_columns(self):
        r"""Return a list with all the names of the columns in the attribute `result_from_query`

        Returns:
            list_of_columns (`list`): List of all the names of the columns
        """
        if self.result_from_query is None:
            msgs.warning('`result_from_query` is empty. You may want to run the query first')
            list_of_columns = []
        else:
            list_of_columns = self.result_from_query.colnames
        return list_of_columns
