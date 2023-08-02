"Finite Element matrix creation"

from .info import __doc__

from .innerproduct import *

__all__ = list(filter(lambda s:not s.startswith('_'),dir()))

