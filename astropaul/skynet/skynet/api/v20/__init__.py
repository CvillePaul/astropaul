"""
SkyNet API: Python interface (client-side) modules
"""

import os.path
import re

# Automatically detect the current API version
__api_version__ = re.search(
    r'skynet{0}api{0}v(\d\d+){0}'.format('\\' + os.path.sep),
    __file__).groups()[0]
