r""" Load default values
"""

import pkg_resources
from ESOAsg import msgs


class Default:
    r"""Set default values for the ESOAsg module
    """

    def __init__(self):
        self.default_dict = {}
        self._load_from_file()

    def _load_from_file(self):
        r"""Load default dict with values in default.txt

        """
        default_file = pkg_resources.resource_filename('ESOAsg', 'default.txt')
        default_list = [line.strip().replace('==', '::')
                        for line in open(default_file) if not line.strip().startswith('#') and line.strip() != '']

        for default in default_list:
            default_quantity, default_value = default.split('::')
            self.default_dict[default_quantity] = default_value

    def get_value(self, card_name):  # written by Ema 05.03.2020
        r""" Gets the default value from the parameter: card_name
        Args:
            card_name (`str`):
        Returns:
            the default value of `card_name`. If not existing an error is raised
        """
        assert isinstance(card_name, str), '{} not in string format'.format(card_name)
        if card_name in self.default_dict:
            return self.default_dict[card_name]
        else:
            msgs.error('{} Key not present in the default dictionary'.format(card_name))