# from ESOAsg import msgs

import collections

from ESOAsg.datacontainers import eso_prodcatg

__all__ = ['DataContainer']


class DataContainer:
    r"""Base class that dictate the general behaviour of a data container

    Attributes:
        prodcatg (:obj:`eso_prodcatg.ProdCatg`)

    """

    def __init__(self, prodcatg_type=None):
        r"""Instantiate the class DataContainer

        """
        self.prodcatg = eso_prodcatg.ProdCatg(prodcatg_type=prodcatg_type)




