from ESOAsg import msgs
from ESOAsg import default
from ESOAsg.core import tap_queries
from ESOAsg.queries import query


class ESOObservations(query.Query):
    r"""This class is designed to query the ESO archive for raw, reduced, and ambient data.

    """

    def __init__(slef, query=None, result_from_query=None, maxrec=default.get_value('maxrec')):
        # assign tap_service
        super().__init__(tap_service=tap_queries.define_tap_service('eso_tap_obs'), query=query,
                         result_from_query=result_from_query, maxrec=maxrec)
    '''
    def __init__(self, tap_service, query):
        self.tap_service = tap_queries.define_tap_service('eso_tap_obs')
        self.query = tap_queries.query_obscore_base()

    @property
    def query(self):
        r"""Return the query
        """
        return self._query

    @property
    def tap_service(self):
        r"""Return the `TAP services<http://archive.eso.org/programmatic/#TAP>`_ used
        """
        return self._tap_service

    def print_query(self):
        r"""Print the query on the terminal
        """
        tap_queries(self.query)
    '''
