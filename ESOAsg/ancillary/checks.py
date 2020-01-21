"""
Module that performs basic checks on your machine
"""

import sys
import shutil

# ESOAsg imports
from ESOAsg import msgs
from ESOAsg import default

def check_disk_space():
    """Given a limit in GB
    """
    total, used, free = shutil.disk_usage("./")
    total = total / (1024**3)
    used  = used / (1024**3)
    free  = free / (1024**3)
    msgs.info("Your disk shows:  Total: {} GB - Used: {} GB - Free: {} GB".format(total, used, free))
