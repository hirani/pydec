"PyDEC mesh and array IO"

from .info import __doc__

from .meshio import *
from .arrayio import *

__all__ = list(filter(lambda s:not s.startswith('_'),dir()))

