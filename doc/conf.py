"""Sphinx configuration file for an LSST stack package.

This configuration only affects single-package Sphinx documentation builds.
"""

from documenteer.sphinxconfig.stackconf import build_package_configs
import lsst.obs.hsc.sims


_g = globals()
_g.update(build_package_configs(
    project_name='obs_hsc_sims',
    version=lsst.obs.hsc.sims.version.__version__))
