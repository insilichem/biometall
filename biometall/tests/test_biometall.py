"""
Unit and regression test for the biometall package.
"""

# Import package, test suite, and other packages as needed
import biometall
import pytest
import sys

def test_biometall_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "biometall" in sys.modules
