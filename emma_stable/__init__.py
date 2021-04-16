"""Extended Module Materials Assembly (EMMA)"""

from distutils.version import LooseVersion

import ase
import numpy

from emma_stable.mc_emma import *
from emma_stable.all import *
from emma_stable.restart_gulp import *
from emma_stable.gulp import *

__all__ = ['mc_emma','all','restart_gulp','gulp','make_random_structure']

__version__ = '4.00'

if LooseVersion(np.__version__) < '1.9':
    # Make isinstance(x, numbers.Integral) work also for np.intxx:
    import numbers
    numbers.Integral.register(np.integer)
    del numbers
