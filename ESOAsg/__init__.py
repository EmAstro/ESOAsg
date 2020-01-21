"""
ESOAsg package initialization.
The current main purpose of this is to provide package-level globals.
"""

# Imports for signal and log handling
import os
import sys
import signal
import warnings


# Set version
__version__ = '0.00'


# Import and instantiate the logger
from pypeit import pypmsgs
msgs = pypmsgs.Messages()

from ESOAsg import check_requirements  # THIS IMPORT DOES THE CHECKING.  KEEP IT

