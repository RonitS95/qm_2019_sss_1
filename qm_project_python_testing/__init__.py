"""
qm_2019_sss_1_cookie
qm_2019_sss_1
"""

# Add imports here
from .qm_project_python_testing import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
