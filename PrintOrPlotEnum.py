from enum import Enum

class PRINTORPLOT(Enum):
    """
    Enum for specifying type of output.
    """
    PRINT = 1
    PLOT = 2
    PLOTSAVE = 3
    PLOTGRAPHONLY = 4
