# from ESOAsg.core import eso_headers

# from collections import OrderedDict


class EsoHeaderParam:
    """Define a class used to hold ESO header parameters
    """

    def __init__(self, cards=None, descriptions=None, prodcatg=None, mandatory=None):
        self.cards = cards
        self.descriptions = descriptions
        self.prodcatg = prodcatg
        self.mandatory = mandatory