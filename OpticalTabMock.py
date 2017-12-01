#!/usr/bin/env python2
# -*- coding:utf-8 -*-
import mock
from Strata import Strata
pyqt5 = False
if pyqt5:
    from PyQt5.QtCore import *
    from PyQt5.QtGui import *
    from PyQt5.QtWidgets import *
    # QString = unicode
else:
    from PyQt4.QtCore import *
    from PyQt4.QtGui import *

class OpticalTab(QWidget):
    dirty = pyqtSignal()

    def __init__(self, parent=None):
        super(OpticalTab, self).__init__(parent)
        self.strata = Strata()
        self.editOpticalParametersBox = mock.Mock()
        self.opticalCanvas = mock.Mock()
        self.optimization1DCanvas = mock.Mock()

    def update_stratum_inputBoxes(self):
        pass

    def stratumTable_refresh(self):
        pass

    def update_opticalCanvas(self):
        pass

# vim: ts=4 sw=4 sts=4 expandtab
