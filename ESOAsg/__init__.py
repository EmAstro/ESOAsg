# Set version
__version__ = '0.1.0'

from ESOAsg import check_requirements  # THIS IMPORT DOES THE CHECKING.  KEEP IT
from ESOAsg import msgs
from ESOAsg import default

# Import and instantiate the logger
msgs = msgs.Messages()

# Define default values
default = default.Default()
