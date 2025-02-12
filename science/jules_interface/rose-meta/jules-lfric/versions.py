import re
import sys

from metomi.rose.upgrade import MacroUpgrade


class UpgradeError(Exception):
    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__


"""
Copy this template and complete to add your macro
class vnXX_txxx(MacroUpgrade):
    # Upgrade macro for <TICKET> by <Author>
    BEFORE_TAG = "vnX.X"
    AFTER_TAG = "vnX.X_txxx"
    def upgrade(self, config, meta_config=None):
        # Add settings
        return config, self.reports
"""


class vn20_t249(MacroUpgrade):
    """Upgrade macro for ticket #249 by Denis Sergeev."""

    BEFORE_TAG = "vn2.0"
    AFTER_TAG = "vn2.0_t249"

    def upgrade(self, config, meta_config=None):
        # Commands From: science/jules_interface/rose-meta/jules-lfric
        nml = "namelist:specified_surface"
        self.add_setting(config, [nml, "surf_temp_forcing"], "'none'")
        self.add_setting(config, [nml, "internal_flux_method"], "'uniform'")
        self.add_setting(config, [nml, "internal_flux_value"], "0.0")

        return config, self.reports
