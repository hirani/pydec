"General utility functions"

from info import __doc__

from util import *

__all__ = filter(lambda s:not s.startswith('_'),dir())

