from ESOAsg import msgs
from ESOAsg import default
from ESOAsg.core import tap_queries
from ESOAsg.queries import query


class ESOObservations(query.Query):
    r"""This class is designed to query the ESO archive for raw, reduced, and ambient data.

    """
    tap_service = tap_queries.define_tap_service('eso_tap_obs')

    def __init__(self):

        super(ESOObservations, self).__init__()

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
