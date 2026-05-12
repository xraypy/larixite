"""Shared pytest configuration for the larixite test suite."""

import warnings


def pytest_configure(config):
    """Suppress pymatgen rounding warnings that occur inside Flask handlers
    (where pytest.warns can't capture them). The warnings are explicitly
    tested in test_webapp.py::TestPymatgenRoundingWarning."""
    warnings.filterwarnings("ignore",
                            message=".*fractional coordinates rounded.*",
                            category=UserWarning,
                            module="larixite.cif_cluster")
    warnings.filterwarnings("ignore",
                            message=".*fractional coordinates rounded.*",
                            category=UserWarning,
                            module="larixite.struct")
