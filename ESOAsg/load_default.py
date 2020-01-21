"""
Load default values
"""

import pkg_resources

from ESOAsg import msgs

class Default:
    """Set default values ESOAsg
    """

    def __init__(self):
        self.default = {}

    def _load_from_file(self):
        """Load default dict with values in default.txt
        """

        default_file = pkg_resources.resource_filename('ESOAsg', 'default.txt')
        default_list = [line.strip().replace('=', ':') for line in open(default_file)
                   if not line.strip().startswith('#') and line.strip() != '']

        msgs.info("Loading default variables")
        for default in default_list:
            default_quantity, default_value = default.split(':')
            self.default[default_quantity] = default_value

