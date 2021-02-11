r""" Module to run tests on queries
"""

from ESOAsg.queries import query_observations
from ESOAsg import archive_catalogues

# from ESOAsg.queries import query_catalogues
# from ESOAsg import archive_observations


def test_query_observations():
    r"""Test for from ESOAsg.queries import query_observations
    """
    query = "SELECT instrument_name FROM ivoa.ObsCore WHERE instrument_name = 'OMEGACAM'"
    lpo_query = query_observations.ESOObservations(query=query, maxrec=1)
    lpo_query.run_query()
    assert lpo_query.result_from_query['instrument_name'] == 'OMEGACAM', r'Unable to find LPO data'


def test_all_catalogues():
    r"""Test for archive_catalogues
    
    Running all_catalogues_info and check that a catalogue is returned
    """
    all_catalogues = archive_catalogues.all_catalogues_info(all_versions=True, verbose=False)
    assert 'GAIAESO' in all_catalogues['collection'].data.data.tolist(), r'GAIAESO not in all_catalogues'
