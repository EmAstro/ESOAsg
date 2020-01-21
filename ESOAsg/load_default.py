"""
Load default values
"""

import pkg_resources

class Default:
    """Set default values ESOAsg
    """

    def __init__(self):
        self.default_dict = {}

    def _load_from_file(self):
        """Load default dict with values in default.txt
        """

        default_file = pkg_resources.resource_filename('ESOAsg', 'default.txt')
        default_list = [line.strip().replace('=', ':') for line in open(default_file)
                   if not line.strip().startswith('#') and line.strip() != '']

        for default in default_list:
            default_quantity, default_value = default.split(':')
            self.default_dict[default_quantity] = default_value

    def get_value(self, card_name):
        """Gets the default value from the parameter: card_name
        """
        return(self.default_dict[card_name])
