# from ESOAsg import msgs

__all__ = ['DataContainer']

PRODCATG = ['SCIENCE.IMAGE',
            'SCIENCE.MEFIMAGE',
            'SCIENCE.IMAGE.FLUXMAP',
            'SCIENCE.SPECTRUM',
            'SCIENCE.SRCTBL',
            'SCIENCE.CUBE.IFS',
            'SCIENCE.VISIBILITY',
            'SCIENCE.CATALOG',
            'SCIENCE.MCATALOG',
            'SCIENCE.CATALOGTILE'
            ]


class DataContainer:
    r"""Base class that dictate the general behaviour of a data container

    Attributes:

    """

    def __init__(self, primary_header=None, header=None, data=None, error=None, background=None, quality=None,
                 mask=None, others=None):
        r"""Instantiate the class LightCurves

        """
        self.primary_header = primary_header
        self.header = header
        self.data = data
        self.error = error
        self.background = background
        self.quality = quality
        self.mask = mask
        self.others = others


