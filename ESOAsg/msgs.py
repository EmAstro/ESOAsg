"""
Module for terminal log
This is inspired by PypeIt: 
https://github.com/pypeit/PypeIt
"""

import sys

from ESOAsg import __version__ 


class Messages:
    r"""Create coloured text for messages printed to screen.
    For further details on colours see the following example:
    http://ascii-table.com/ansi-escape-sequences.php

    Args:
        colors (`bool`):
            If true, the screen output will have colors, otherwise
            normal screen output will be displayed
    """

    def __init__(self, colors=True):

        # Initialize variables
        self._version = __version__

        # Use colors?
        self._start = None
        self._end = None
        self._black_CL = None
        self._yellow_CL = None
        self._blue_CL = None
        self._green_CL = None
        self._red_CL = None
        self._white_RD = None
        self._white_GR = None
        self._white_BK = None
        self._white_BL = None
        self._black_YL = None
        self._yellow_BK = None

        self.disable_colors()
        if colors:
            self.enable_colors()

    @staticmethod
    def _print(premsg, msg):
        r"""Print to standard error
        """
        _msg = premsg+msg
        print(_msg, file=sys.stderr)

    def error(self, msg):
        r"""Print an error message
        """
        premsg = '\n'+self._start + self._white_RD + '[ERROR]   ::' + self._end + ' '
        self._print(premsg, msg)
        self.newline()
        sys.exit(1)

    def info(self, msg):
        r"""Print an information message
        """
        premsg = self._start + self._green_CL + '[INFO]    ::' + self._end + ' '
        self._print(premsg, msg)

    def start(self):
        r"""Print a staring message
        """
        premsg = self._start + self._green_CL + '[START]   ::' + self._end + ' '
        self._print(premsg, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

    def end(self):
        r"""Print a ending message
        """
        premsg = self._start + self._red_CL + '[END]     ::' + self._end + ' '
        self._print(premsg, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

    def warning(self, msg):
        r"""Print a warning message
        """
        premsg = self._start + self._red_CL + '[WARNING] ::' + self._end + ' '
        self._print(premsg, msg)

    def bug(self, msg):
        r"""Print a bug message
        """
        premsg = self._start + self._white_BK + '[BUG]     ::' + self._end + ' '
        self._print(premsg, msg)

    def work(self, msg):
        r"""Print a work in progress message
        """
        premsg = self._start + self._black_CL + '[WORKING] ::' + self._end + ' '
        self._print(premsg, msg)

    def pre_indent(self, msg):
        r"""Print an indent
        """
        premsg = '             '
        self._print(premsg, msg)

    @staticmethod
    def newline():
        r"""Return a text string containing a newline to be used with messages
        """
        return '\n             '

    @staticmethod
    def indent():
        r"""Return a text string containing an indent to be used with messages
        """
        return '             '

    # Set the colors
    def enable_colors(self):
        r"""Enable colored output text
        """

        # Start and end coloured text
        self._start = '\x1B['
        self._end = '\x1B[' + '0m'

        # Clear Backgrounds
        self._black_CL = '1;30m'
        self._yellow_CL = '1;33m'
        self._blue_CL = '1;34m'
        self._green_CL = '1;32m'
        self._red_CL = '1;31m'

        # Coloured Backgrounds
        self._white_RD = '1;37;41m'
        self._white_GR = '1;37;42m'
        self._white_BK = '1;37;40m'
        self._white_BL = '1;37;44m'
        self._black_YL = '1;37;43m'
        self._yellow_BK = '1;33;40m'

    def disable_colors(self):
        r"""Disable colored output text
        """

        # Start and end coloured text
        self._start = ''
        self._end = ''

        # Clear Backgrounds
        self._black_CL = ''
        self._yellow_CL = ''
        self._blue_CL = ''
        self._green_CL = ''
        self._red_CL = ''

        # Coloured Backgrounds
        self._white_RD = ''
        self._white_GR = ''
        self._white_BK = ''
        self._white_BL = ''
        self._black_YL = ''
        self._yellow_BK = ''
