r""" Module to run tests on queries
"""

from ESOAsg.queries import query_observations


def test_query_observations():
    r"""Test for from ESOAsg.queries import query_observations
    """
    query = "SELECT instrument_name FROM ivoa.ObsCore WHERE instrument_name = 'ALMA'"
    alma_query = query_observations.ESOObservations(query=query, maxrec=1)
    alma_query.run_query()
    assert alma_query.result_from_query['instrument_name'] == 'ALMA'
