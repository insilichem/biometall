"""
biometall
Software to predict feasible metal binding areas into proteins
"""
from .biometall import run

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
