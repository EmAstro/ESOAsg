# from ESOAsg.core import eso_headers
from astropy.io.fits import header
from ESOAsg.datacontainers import eso_prodcatg
from ESOAsg import msgs


class ESOHeader(header.Header):

    def __init__(self, cards=[], copy=False, is_primary=False, prodcatg_type=None):
        header.Header.__init__(self, cards=cards, copy=copy)
        self.is_primary = is_primary
        self.prodcatg = prodcatg_type

    @property
    def prodcatg(self):
        r"""Name of the `prodcatg`.

        This set also the `descriptor` dictionary containing the relevant properties of the `prodcatg`.

        """
        return self.__prodcatg

    @prodcatg.setter
    def prodcatg(self, prodcatg_type):
        self.__prodcatg = eso_prodcatg.ProdCatg(prodcatg_type=prodcatg_type)
        if self.is_primary is True:
            _prodcatg_value = self.get('PRODCATG', default=None)
            if _prodcatg_value is None:
                msgs.info('Added PRODCATG = {} to the header'.format(prodcatg_type))
                self.set('PRODCATG', prodcatg_type)
            elif _prodcatg_value != prodcatg_type:
                msgs.warning('Updating value fo PRODCATG from {} to {}'.format(_prodcatg_value, prodcatg_type))
                self.set('PRODCATG', prodcatg_type)
            elif _prodcatg_value == prodcatg_type:
                msgs.info('PRODCATG = {}'.format(_prodcatg_value))
            else:
                msgs.error('Cannot set the value of PRODCATG')

    def data_format(self):
        r"""Returns a string describing the data format associated to the `prodcatg`.

        Returns:
            str: string containing the general description of the data format of the given `prodcatg`.

        """
        return self.prodcatg.data_format()
