"""
Module that performs basic checks on your machine
"""

import sys
import shutil
import numpy as np

# ESOAsg imports
from ESOAsg import msgs
from ESOAsg import default

def check_disk_space(min_disk_space=np.float32(default.get_value('min_disk_space'))):
    """Given a limit in GB in the variable min_disk_space the macro returns
    true if there is enough space and rises an error otherwise.

    Parameters
    ----------
    check_disk_space : np.float32
        Size of free space on disk required

    Returns
    -------
    enough_space : np.bool
        True if there is enough space on disk
    """
    total, used, free = shutil.disk_usage("./")
    total = total / (1024**3)
    used  = used / (1024**3)
    free  = free / (1024**3)
    msgs.info('Your disk shows: Total: {0:.2f} GB, Used: {0:.2f} GB, Free: {0:.2f} GB'.format(total, used, free))
    if (free>min_disk_space):
        enough_space = np.bool(1)
    else:
        enough_space = np.bool(0)
        msgs.error('Not enough space on disk')
    return enough_space
