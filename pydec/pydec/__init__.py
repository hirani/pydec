"""PyDEC: Software and Algorithms for Discrete Exterior Calculus
"""

from .version import version as __version__

from .dec import *
from .fem import *
from .math import *
from .io import *
from .mesh import *
from .util import *
from .vis import *

__all__ = list(filter(lambda s:not s.startswith('_'),dir()))
#__all__ += ['test', '__version__']
__all__ += ['__version__']

# TODO: replace testing framework with pytest (?) since Tester of nose/numpy is deprecated
#from pydec.testing import Tester
#test = Tester().test
