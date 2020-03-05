r""" ESOAsg package initialization.

The current main purpose of this is to provide package-level globals.
"""

# Set version
__version__ = '0.00'

from ESOAsg import check_requirements  # THIS IMPORT DOES THE CHECKING.  KEEP IT
from ESOAsg import msgs
from ESOAsg import default

# Import and instantiate the logger
msgs = msgs.Messages()

# Define default values
default = default.Default()
