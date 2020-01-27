"""
Module to compare lists
"""


class Lists:
    r"""
    A class used to check lists.

    This allows to perform some checks on an input list (could also be given with a text files or from a fits file
    header).
    
    Attributes:
        cards (`list`):
            The list you want to test.
        values (`list`)
            If not `None` these are the values to be associated to `cards`.

    Methods:


    """

    def __init__(self, cards=None, values=None, from_header=None,
                 from_txt=None):
        r"""

        Args:
            cards
            values
            from_header
            from_txt

        Returns:

        """

        if cards is not None:
            self.cards = list(cards)
        elif from_txt is not None:
            msgs.working('Loading list from text file {}'.format(from_txt))
            from_txt_file=open(from_txt,"r")
            lines=from_txt_file.readlines()
            full=[]
            for line in lines:
                full.append(line.split(' ')[1])
            from_txt_file.close()
        elif from_header is not None:
            msgs.working('Loading header from file {}'.format(from_header))
        else:
            msgs.working('Creating empty list object.')
            self.cards = []
        return
