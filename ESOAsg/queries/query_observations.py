from ESOAsg import msgs
from ESOAsg import default
from ESOAsg.core import tap_queries
from ESOAsg.queries import query


class ESOObservations(query.Query):
    r"""This class is designed to query the ESO archive for raw, reduced, and ambient data.

    """

    def __init__(self,  query=None, result_from_query=None, type_of_query='sync', maxrec=default.get_value('maxrec')):
        # assign tap_service
        super().__init__(tap_service=tap_queries.define_tap_service('eso_tap_obs'), query=query,
                         result_from_query=result_from_query,  type_of_query=type_of_query, maxrec=maxrec)