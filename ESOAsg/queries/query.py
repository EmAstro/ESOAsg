from ESOAsg import default
from ESOAsg import msgs
from ESOAsg.ancillary import checks
from ESOAsg.core import tap_queries


class Query:
    r"""Base class that dictate the general behaviour of a query

    Attributes:
        tap_service (`str`)
        query (`str`)
        result_from_query (`str`)
        maxrec (`int`)
    """

    def __init__(self, tap_service=None, query=None, result_from_query=None, maxrec=default.get_value('maxrec')):
        self.tap_service = tap_service
        self.query = query
        self.result_from_query = result_from_query
        self.maxrec = maxrec

    def clean_query(self):
        r"""Set the `query` attribute to None
        """
        self.query = None

    def clean_result_from_query(self):
        r"""Set the `result_from_query` attribute to None
        """
        self.result_from_query = None

    def run_query(self, to_string=True):
        r"""Run the query and store results in the `result_from_query` attribute

        Args:
            to_string (`bool`, optional): if set to True, if a column is in `bytes` format it transform it to `str`

        """
        self.result_from_query = tap_queries.run_query(self.tap_service, self.query)
        if to_string:
            for column_id in self.which_columns():
                self.result_from_query[column_id].data.data[:] = checks.from_bytes_to_string(self.result_from_query[
                    column_id].data.data)
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
            list_of_columns = self.which_columns.colnames
        return list_of_columns