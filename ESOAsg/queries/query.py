from ESOAsg import default
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

    def run_query(self):
        r"""Run the query and store results in result_from_query
        """
        self.result_from_query = tap_queries.run_query(self.tap_service, self.query)

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
        tap_queries.which_service(self.tap_service)
        return
