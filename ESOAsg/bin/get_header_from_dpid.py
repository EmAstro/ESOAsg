#!/usr/bin/env python3

import argparse
# from astropy import coordinates
# from astropy import units as u
import numpy as np
import os.path
# import glob

# from ESOAsg.core import download_archive
from ESOAsg import msgs
from ESOAsg import __version__
from ESOAsg.core import fitsfiles
from ESOAsg.core import archive_observations
from ESOAsg.ancillary import cleaning_lists

# from IPython import embed


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=r"""
        This macro reads the gets the header of a data product (DP) stored in the ESO archive.
        In case the --cards argument is not included, it will download the full header of the data
        product on your disk.
        
        This uses ESOAsg version {:s}
        """.format(__version__),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES)

    parser.add_argument('input_dpid', nargs='+', type=str,
                        help='data product ID')
    parser.add_argument('-car', '--cards', nargs='+', type=str, default=None,
                        help='Header cards you want to check')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    return parser.parse_args()


EXAMPLES = r"""
        Examples:
        Get the the value of "HIERARCH ESO OBS ID" from the dpid "OMEGA.2016-08-28T07:14:36.769"
        get_header_from_dpid.py 'OMEGA.2016-08-28T07:14:36.769' -car 'HIERARCH ESO OBS ID'
        """

if __name__ == '__main__':
    # input arguments and tests
    args = parse_arguments()

    # File_names
    input_dpids = cleaning_lists.make_list_of_strings(args.input_dpid)

    # cards
    if args.cards is not None:
        input_cards = cleaning_lists.make_list_of_strings(args.cards)

    # Starting to read the header(s)
    msgs.start()

    if args.cards is not None:
        dpid_names = []
        card_names = []
        card_values = []

    for input_dpid in input_dpids:
        archive_observations.get_header_from_archive(input_dpid, text_file=input_dpid+'.hdr')
        if args.cards is not None:
            hdr_from_dpid = fitsfiles.header_from_txt_file(input_dpid+'.hdr')
            for input_card in input_cards:
                dpid_names.append(input_dpid)
                card_names.append(input_card)
                card_values.append(hdr_from_dpid[input_card])

    if args.cards is not None:
        msgs.info('Summary:')
        msgs.info('DPID   CARD   VALUE')
        for dpid_name, card_name, card_value in zip(dpid_names, card_names, card_values):
            print('{}   {}   {}'.format(dpid_name, card_name, card_value))

    msgs.end()
