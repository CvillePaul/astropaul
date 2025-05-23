"""
SkyNet API: client-side part and modules common to both client and server sides

API users requiring a specific API version should explicitly import from the
corresponding subpackage, e.g.

from skynet.api.v20.errors import ApiError

Users may also choose to import without an explicit API version specification,
like

from skynet.api.errors import ApiError

In the latter case, modules within skynet.api itself (that use absolute vs
relative imports, like "from .errors import ApiError") get the same API version
as they belong to, as well as modules within web2py.applications.api*. For such
modules, importing without explicit version specification is the preferred way
since it makes maintenance easier (no need to change imports when developing a
new API version).

Modules outside the SkyNet API (like user scripts) may do that as well, then
they always get the most recent API version. However, this is only recommended
for features that are not likely to change their interface on API upgrades.
Beware that compatibility can be broken at any point on major API upgrades, so
it is up to you to maintain compatibility.
"""

from __future__ import absolute_import, division, print_function

import sys
import os.path
import re
from glob import glob

try:
    # noinspection PyCompatibility
    from importlib.abc import MetaPathFinder
except ImportError:
    MetaPathFinder = object


auto_api_version_regex = re.compile(r'skynet\.api\.(?!v\d+)')
api_package_regex = re.compile(r'v(\d\d+)')
web2py_api_regex = re.compile(
    r'applications{0}api(\d\d+){0}'.format('\\' + os.path.sep))
skynet_api_regex = re.compile(
    r'skynet{0}api{0}v(\d\d+){0}'.format('\\' + os.path.sep))


# Current API version - guess from subpackage names
__current_api_version__ = sorted(
    [api_package_regex.match(_p).groups()[0]
     for _p in [os.path.basename(__p)
                for __p in glob(os.path.join(os.path.dirname(__file__), 'v*'))]
     if api_package_regex.match(_p)],
    key=lambda p: (int(p[0]), int(p[1:])))[-1]


class SkyNetAPIModuleFinder(MetaPathFinder):
    """
    SkyNet API import hook class

    Uses SkyNetAPIModuleLoader to load modules for
    "from skynet.api.errors"-style imports. An instance of this class should be
    inserted in sys.meta_path.
    """
    api_ver = None

    def find_module(self, fullname, _=None):
        """
        Find the given module; for modules within skynet.api (excluding those
        with an explicit API version spec, like skynet.api.v20.errors), use
        SkyNetAPIModuleLoader; for all other modules, fall back to the normal
        Python module loader

        :param fullname: fully-qualified module name, e.g. skynet.api.errors
        :param _: module search path (unused)

        :return: SkyNetAPIModuleLoader instance for API modules without an
            explicit version spec; None for all other modules
        """
        if auto_api_version_regex.match(fullname):
            # By default, assume that the caller needs the most recent API
            self.api_ver = __current_api_version__

            # Guess the caller's API version by traversing the stack until me
            # encounter a module with either "applications/api##" or
            # "skynet/api/v##" in its path. If we reach the bottom without
            # reaching such module, then we're called not from within an
            # API-specific part of the library, but from a user script. We are
            # here only when the user doesn't import the specific API (e.g.
            # "from skynet.api.v20.import errors"). Then it is reasonable to
            # assume that the user wants the most recent API, so we may leave
            # api_ver as is
            # noinspection PyProtectedMember
            frame = sys._getframe(1)
            try:
                while frame:
                    p = frame.f_code.co_filename
                    try:
                        # In web2py/applications/api*?
                        self.api_ver = web2py_api_regex.search(p).groups()[0]
                    except AttributeError:
                        try:
                            # In skynet.api.v*?
                            self.api_ver = \
                                skynet_api_regex.search(p).groups()[0]
                        except AttributeError:
                            pass
                        else:
                            # Called from skynet.api.v*
                            break
                    else:
                        # Called from web2py/applications/api*
                        break

                    frame = frame.f_back
            finally:
                del frame

            return self

    def load_module(self, fullname):
        """
        Load the given SkyNet API module

        :param fullname: fully-qualified module name, e.g. skynet.api.errors;
            the appropriate API version is appended automatically

        :return: module instance
        """
        # Adjust the module name to include the appropriate API version
        mod = __import__(
            auto_api_version_regex.sub(
                'skynet.api.v' + self.api_ver + '.', fullname),
            fromlist=['__dict__'])

        # Also add module to cache with its generic name (w/o API version
        # number) - required by Python 3
        sys.modules[fullname] = mod

        return mod


sys.meta_path.append(SkyNetAPIModuleFinder())
