from ESOAsg.filters import filter_curves


def test_all_filters():
    r"""Test for data access to all filters

    Reading the file containing all filter information and test if `SPHERE` is in the instrument_name
    """
    all_filters_table = filter_curves._read_filters_file()
    assert 'SPHERE' in all_filters_table['instrument_name'].data.tolist(), r'SPHERE not in all_filters'
