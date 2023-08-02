"PyDEC meshes"

from .info import __doc__

from .simplex import *
from .regular_cube import *
from .generation import *
from .subdivision import *
from .ncube import *

__all__ = list(filter(lambda s:not s.startswith('_'),dir()))

# TODO: replace testing framework with pytest (?) since Tester is deprecated
#from pydec.testing import Tester
#test = Tester().test

