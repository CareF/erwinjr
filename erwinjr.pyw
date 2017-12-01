#!/usr/bin/env python2
# -*- coding:utf-8 -*-

# ===========================================================================
# ErwinJr is a simulation program for quantum semiconductor lasers.
# Copyright (C) 2012 Kale J. Franz, PhD
# Copyright (C) 2017 Ming Lyu (CareF)
#
# A portion of this code is Copyright (c) 2011, California Institute of
# Technology ("Caltech"). U.S. Government sponsorship acknowledged.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ===========================================================================

# TODO:
# replace np.hstack (done, but self.strata part to test)
# find replacement for psyco
# try to seperate this file to smaller ones (done for quantumtab
# check unnecessary function call
# Ctrl+z support
# add status bar

from __future__ import division

__pyqt5__ = False
__USE_MATPLOTLIB__ = True

import os
import sys
import traceback
from functools import partial
import time

from QCLayers import QCLayers, cst
from Strata import Strata
import SaveLoad

if __pyqt5__:
    from PyQt5.QtCore import *
    from PyQt5.QtGui import *
    from PyQt5.QtWidgets import *
    # QString is automatically unicode
    # qsettings.value(...).toxxxx -> qsettings.value(...)
else:
    from PyQt4.QtCore import *
    from PyQt4.QtGui import *

if __USE_MATPLOTLIB__:
    from QuantumTabMatplotlib import QuantumTab
    from OpticalTab import OpticalTab
else:
    from QuantumTab import QuantumTab
    from OpticalTab import OpticalTab

# ===========================================================================
# Version
# ===========================================================================
ejVersion = 171109
majorVersion = '3.4.0'

# ===========================================================================
# Debug options
# ===========================================================================
DEBUG = 1


class MainWindow(QMainWindow):
    def __init__(self, fileName=None, parent=None):
        super(MainWindow, self).__init__(parent)

        self.filename = fileName

        self.create_main_frame()
        self.create_Quantum_menu()

        qsettings = QSettings(parent=self)
        self.recentFiles = qsettings.value("RecentFiles").toStringList()
        self.restoreGeometry(
                qsettings.value("MainWindow/Geometry").toByteArray())
        self.restoreState(qsettings.value("MainWindow/State").toByteArray())
        self.updateFileMenu()

        self.dirty = False

        if self.filename:
            self.fileOpen(self.filename)
        else:
            QTimer.singleShot(0, self.loadInitialFile)

        self.dirty = False
        self.update_windowTitle()

# ===========================================================================
# Create Main Frame
# ===========================================================================

    def create_main_frame(self):

        self.mainTabWidget = QTabWidget()

        # ###############################
        #
        # Thermal Tab
        #
        # ###############################

#        vBox1 = QVBoxLayout()

#        thermalTable = QTableWidget()
#        thermalTable.setSelectionBehavior(QTableWidget.SelectRows)
#        thermalTable.setSelectionMode(QTableWidget.SingleSelection)
#        thermalTable.setMaximumWidth(380)
#        thermalTable.setMinimumWidth(380)
#        thermalTable.itemChanged.connect(self.stratumTable_itemChanged)
#        thermalTable.itemSelectionChanged.connect(self.stratumTable_itemSelectionChanged)
#        vBox1.addWidget(thermalTable)

#        thermalWidget = QWidget()
#        thermalWidget.setLayout(vBox1)
#        thermalWidget.setAutoFillBackground(True)
#        thermalWidget.setBackgroundRole(QPalette.Window)

        # ##########################
        #
        # Quantum Tab
        #
        # ##########################
        self.quantumWidget = QuantumTab(self)
        self.quantumWidget.dirty.connect(self.winUpdate)
        self.quantumWidget.toOpticalParamsButton.clicked.connect(
                self.transfer_optical_parameters)
        self.mainTabWidget.addTab(self.quantumWidget, 'Quantum')

        # ##########################
        #
        # Optical Tab
        #
        # ##########################
        self.opticalWidget = OpticalTab(self)
        self.mainTabWidget.addTab(self.opticalWidget, 'Optical')
        #  self.mainTabWidget.addTab(thermalWidget, 'Thermal')
        self.mainTabWidget.currentChanged.connect(self.change_main_tab)
        #  self.opticalWidget.dirty.connect(self.winUpdate)
        self.connect(self.opticalWidget,
                     SIGNAL('dirty'),
                     self.winUpdate)

        self.setCentralWidget(self.mainTabWidget)

    def winUpdate(self):
        """ SLOT connected to self.quantumWidget.dirty and
        self.opticalWidget.dirty """
        self.dirty = True
        self.update_windowTitle()

    def transfer_optical_parameters(self):
        """ SLOT connected to quantumWidget.toOpticalParamsButton.clicked()
        Transfer informatin in quantum tab to optical tab"""
        self.quantumWidget.transfer_params(self.opticalWidget.strata)

        # GUI settings
        self.quantumWidget.toOpticalParamsButton.setEnabled(False)
        self.opticalWidget.editOpticalParametersBox.setChecked(False)

        # update all the input boxes
        self.opticalWidget.update_stratum_inputBoxes()

        # update the stratumTable
        self.opticalWidget.stratumTable_refresh()

# ===========================================================================
# General Menu Functions
# ===========================================================================
    def change_main_tab(self, tabIdx):
        self.menuBar().clear()
        if tabIdx == 0:
            self.create_Quantum_menu()
        elif tabIdx == 1:
            self.create_Optical_menu()
        elif tabIdx == 2:
            pass
        else:
            assert 1==2

    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def create_action(self, text, slot=None, shortcut=None, icon=None,
                      tip=None, checkable=False, ischecked=False):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon("images/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            action.triggered.connect(slot)
        if checkable:
            action.setCheckable(True)
        if ischecked:
            action.setChecked(True)
        return action

    def create_Quantum_menu(self):
        # file menu
        self.file_menu = self.menuBar().addMenu("&File")
        newFileAction = self.create_action(
            "&New...", slot=self.fileNew, shortcut=QKeySequence.New,
            icon="filenew", tip="New ErwinJr file")
        openFileAction = self.create_action(
            "&Open", shortcut="Ctrl+O", slot=self.fileOpen,
            tip="Open ErwinJr file", icon="fileopen")
        saveFileAction = self.create_action(
            "&Save", shortcut="Ctrl+S", slot=self.fileSave,
            tip="Save ErwinJr file", icon="filesave")
        saveAsFileAction = self.create_action(
            "S&ave As", shortcut="Ctrl+W", slot=self.fileSaveAs,
            tip="Save ErwinJr file as", icon="filesaveas")
        exportQuantumCanvasAction = self.create_action(
                "Export Band Diagram Image",
                slot=self.exportBandDiagram,
                tip="Export Band Diagram Image")
        exportBandCSVAction = self.create_action(
            "Export Band Diagram Data", slot=self.export_band_diagram_data,
            tip="Export Band Diagram Data")
        quit_action = self.create_action(
            "&Quit", slot=self.close, shortcut="Ctrl+Q",
            tip="Close the application", icon="filequit")
        self.fileMenuActions = (
            newFileAction, openFileAction, saveFileAction, saveAsFileAction,
            None, exportBandCSVAction, exportQuantumCanvasAction, None,
            quit_action)
        self.file_menu.aboutToShow.connect(self.updateFileMenu)

        # edit menu
        self.edit_menu = self.menuBar().addMenu("&Edit")
        temperatureAction = self.create_action(
            "&Temperature", slot=self.set_temperature, tip="Set temperature")
        bumpLayerAction = self.create_action(
            "&Bump First Layer", slot=self.quantumWidget.bump_first_layer,
            tip="Move zeroth layer to first layer")
        copyStructureAction = self.create_action(
            "&Copy Structure", slot=self.quantumWidget.copy_structure,
            tip="Copy Layer Structure to Clipboard")
        self.add_actions(self.edit_menu, (temperatureAction,
                                          bumpLayerAction, None,
                                          copyStructureAction))

        # view menu
        self.view_menu = self.menuBar().addMenu("&View")
        VXBandAction = self.create_action(
                "X Valley Conduction Band",
                checkable=True, ischecked=self.quantumWidget.plotVX,
                slot=self.quantumWidget.view_VXBand)
        VLBandAction = self.create_action(
                "L Valley Conduction Band",
                checkable=True, ischecked=self.quantumWidget.plotVL,
                slot=self.quantumWidget.view_VLBand)
        LHBandAction = self.create_action(
                "Light Hole Valence Band",
                checkable=True, ischecked=self.quantumWidget.plotLH,
                slot=self.quantumWidget.view_LHBand)
        SOBandAction = self.create_action(
                "Split Off Valence Band",
                checkable=True, ischecked=self.quantumWidget.plotSO,
                slot=self.quantumWidget.view_SOBand)
        self.add_actions(self.view_menu, (VXBandAction,
                                          VLBandAction,
                                          LHBandAction,
                                          SOBandAction))

        # help menu
        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action("&About", shortcut='F1',
                                          slot=self.on_about)
        licenses_action = self.create_action("&License",
                                             slot=self.on_licenses)
        tutorialAction = self.create_action("&Tutorial",
                                            slot=self.on_tutorial)
        self.add_actions(self.help_menu, (tutorialAction,
                                          about_action,
                                          licenses_action))

    def create_Optical_menu(self):
        # file menu
        self.file_menu = self.menuBar().addMenu("&File")
        newFileAction = self.create_action(
            "&New...", self.fileNew, QKeySequence.New,
            "filenew", "New ErwinJr file")
        openFileAction = self.create_action(
            "&Open", shortcut="Ctrl+O", slot=self.fileOpen,
            tip="Open ErwinJr file", icon="fileopen")
        saveFileAction = self.create_action(
            "&Save", shortcut="Ctrl+S", slot=self.fileSave,
            tip="Save ErwinJr file", icon="filesave")
        saveAsFileAction = self.create_action(
            "S&ave As", shortcut="Ctrl+W", slot=self.fileSaveAs,
            tip="Save ErwinJr file as", icon="filesaveas")
        quit_action = self.create_action(
            "&Quit", slot=self.close, shortcut="Ctrl+Q",
            tip="Close the application", icon="filequit")
        self.fileMenuActions = (newFileAction, openFileAction,
                                saveFileAction, saveAsFileAction,
                                None, quit_action)
        self.file_menu.aboutToShow.connect(self.updateFileMenu)
        # self.add_actions(self.file_menu, (newFileAction, openFileAction, saveFileAction, saveAsFileAction, None, quit_action))

        # help menu
        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action("&About", shortcut='F1',
                                          slot=self.on_about)
        licenses_action = self.create_action("&License",
                                             slot=self.on_licenses)
        tutorialAction = self.create_action("&Tutorial",
                                            slot=self.on_tutorial)
        self.add_actions(self.help_menu, (tutorialAction,
                                          about_action, licenses_action))


# ===========================================================================
# File Menu Items
# ===========================================================================
    def updateFileMenu(self):
        """SLOT connected to self.file_menu.aboutToShow()
        Update for recent files"""
        self.file_menu.clear()
        self.add_actions(self.file_menu, self.fileMenuActions[:-1])
        current = (QString(self.filename)
                   if self.filename is not None else None)
        recentFiles = []
        for fname in self.recentFiles:
            if fname != current and QFile.exists(fname):
                recentFiles.append(fname)
        if recentFiles:
            self.file_menu.addSeparator()
            for i, fname in enumerate(recentFiles):
                action = QAction(
                        "&{0}  {1}".format(i + 1, QFileInfo(
                            fname).fileName()), self)
                action.setData(QVariant(fname))
                action.triggered.connect(partial(self.fileOpen, fname))
                self.file_menu.addAction(action)
        self.file_menu.addSeparator()
        self.file_menu.addAction(self.fileMenuActions[-1])

    def addRecentFile(self, fname):
        if fname is None:
            return
        if __pyqt5__:
            if fname not in self.recentFiles:
                self.recentFiles.insert(0, fname)
                while self.recentFiles.count() > 9:
                    self.recentFiles.pop()
        else:
            if not self.recentFiles.contains(fname):
                self.recentFiles.prepend(QString(fname))
                while self.recentFiles.count() > 9:
                    self.recentFiles.takeLast()

    def loadInitialFile(self):
        qsettings = QSettings()
        fname = unicode(qsettings.value("LastFile").toString())
        if fname and QFile.exists(fname):
            if fname.split('.')[-1] == 'qcl':
                self.qclLoad(fname)

            self.quantumWidget.reloaded()
            self.opticalWidget.strata.populate_rIndexes()
            self.opticalWidget.update_stratum_inputBoxes()
            self.opticalWidget.stratumTable_refresh()
            self.opticalWidget.opticalCanvas.clear()
            self.opticalWidget.update_opticalCanvas()

        self.filename = fname
        self.addRecentFile(fname)
        self.dirty = False
        self.update_windowTitle()

    def fileNew(self):
        if not self.okToContinue():
            return False

        self.filename = None
        #  self.quantumWidget.quantumCanvas.clear()
        self.opticalWidget.opticalCanvas.clear()
        self.opticalWidget.optimization1DCanvas.clear()

        self.quantumWidget.qclayers = QCLayers()
        self.quantumWidget.reloaded()

        self.opticalWidget.strata = Strata()

        self.opticalWidget.update_stratum_inputBoxes()

        self.opticalWidget.update_opticalCanvas()
        self.opticalWidget.stratumTable_refresh()

        self.dirty = False
        self.update_windowTitle()

        return True

    def okToContinue(self):
        if self.dirty:
            reply = QMessageBox.question(
                self, "ErwinJr " + str(majorVersion) + " - Unsaved Changes",
                "Save unsaved changes?",
                QMessageBox.Yes|QMessageBox.No|QMessageBox.Cancel)
            if reply == QMessageBox.Cancel:
                return False
            elif reply == QMessageBox.Yes:
                self.fileSave()
        return True

    def update_windowTitle(self):
        if self.filename is not None:
            self.setWindowTitle("ErwinJr " + str(majorVersion) + " - %s[*]" %
                                os.path.basename(str(self.filename)))
        else:
            self.setWindowTitle("ErwinJr " + str(majorVersion) + "[*]")
        self.setWindowModified(self.dirty)

    def fileOpen(self, fname=None):
        # clear all old data, also calls self.okToContinue()
        if not self.fileNew():
            return False
        if not fname:
            dir = os.path.dirname(str(self.filename)) if self.filename else "."
            fname =unicode(QFileDialog.getOpenFileName(
                self, "ErwinJr - Choose file", dir,
                "ErwinJr files (*.qcl)\nAll files (*.*)"))
        # open file and determine if it is from the Matlab version of ErwinJr
        filehandle = open(fname, 'rU')
        firstLine = filehandle.readline()
        filehandle.close()
        if fname:
            if firstLine.split(':')[0] == 'Description':
                QMessageBox.warning(
                    self, 'ErwinJr Error',
                    'Older .qcl format is no longer supported for Ver>3.0.')
                #  self.qclPtonLoad(fname)
            elif firstLine == 'ErwinJr Data File\n':
                self.qclLoad(fname)
            else:
                QMessageBox.warning(self, 'ErwinJr Error',
                                    'Could not recognize input file.')
                return
            self.quantumWidget.reloaded()

            # if firstLine == 'ErwinJr Data File\n':
            self.opticalWidget.strata.populate_rIndexes()
            self.opticalWidget.update_stratum_inputBoxes()
            self.opticalWidget.stratumTable_refresh()
            self.opticalWidget.opticalCanvas.clear()
            self.opticalWidget.update_opticalCanvas()

        self.filename = fname
        self.addRecentFile(fname)
        self.dirty = False
        self.update_windowTitle()

        return True

    def qclLoad(self, fname):
        #  print "Loading "+fname
        try:
            with open(fname, 'rU') as f:
                SaveLoad.qclLoad(f, self.quantumWidget.qclayers,
                                 self.opticalWidget.strata)
        except Exception as err:
            QMessageBox.warning(self, "ErwinJr - Warning",
                                "Could not load *.qcl file.\n"+
                                traceback.format_exc())

    def fileSave(self):
        if self.filename is None:
            return self.fileSaveAs()
        else:
            # os.path.extsep
            if self.filename.split('.')[-1] == 'qcl':
                if self.qclSave(self.filename):
                    self.dirty = False
                    self.update_windowTitle()
                    return True
                else:
                    return False
            else:
                raise IOError('The *.' + self.filename.split('.')[-1] +
                              ' extension is not supported.')
                return False

    def fileSaveAs(self):
        fname = self.filename if self.filename is not None else "."
        typeString = "ErwinJr 2.x file (*.qcl)\nAll files (*.*)"
        fname = unicode(QFileDialog.getSaveFileName(
            self, "ErwinJr - Save File", QString(fname), typeString))
        if fname:
            if "." not in fname:
                fname += ".qcl"
            self.addRecentFile(fname)
            self.filename = fname
            return self.fileSave()
        return False

    def qclSave(self, fname):
        try:
            with open(fname, 'w') as f:
                f.write("ErwinJr Data File\n")
                f.write("Version:" + str(ejVersion) + '\n')
                SaveLoad.qclSave(f, self.quantumWidget.qclayers,
                                 self.opticalWidget.strata)
        except Exception as err:
            QMessageBox.warning(self, "ErwinJr - Warning",
                                "Could not save *.qcl file.\n"+
                                traceback.format_exc())
        return True

    def closeEvent(self, event):
        if self.okToContinue():
            qsettings = QSettings()
            filename = (QVariant(QString(self.filename)) if self.filename
                        else QVariant())
            qsettings.setValue("LastFile", filename)
            recentFiles = (QVariant(self.recentFiles) if self.recentFiles
                           else QVariant())
            qsettings.setValue("RecentFiles", recentFiles)
            qsettings.setValue(
                    "MainWindow/Geometry", QVariant(self.saveGeometry()))
            qsettings.setValue(
                    "MainWindow/State", QVariant(self.saveState()))
        else:
            event.ignore()


# ===========================================================================
# Export Functions
# ===========================================================================
    def exportBandDiagram(self):
        if __USE_MATPLOTLIB__:
            self.quantumWidget.export_quantumCanvas(
                    self.filename.split('.')[0])
        else:
            fname = unicode(QFileDialog.getSaveFileName(
                self, "ErwinJr - Export Band Structure Image",
                self.filename.split('.')[0],
                "Portable Network Graphics file (*.png)"))
            if not fname:
                return

            # set background color to white and save presets
            bgRole = self.mainTabWidget.backgroundRole()
            self.mainTabWidget.setBackgroundRole(QPalette.Base)
            self.quantumWidget.export_quantumCanvas(fname)
            self.mainTabWidget.setBackgroundRole(bgRole)

    def export_band_diagram_data(self):
        fname = unicode(QFileDialog.getSaveFileName(
            self, "ErwinJr - Export Band Structure Data",
            self.filename.split('.')[0],
            "Comma-Separated Value file (*.csv)"))
        if fname != '':
            # if user doesn't click cancel
            self.quantumWidget.export_band_data(fname)


# ===========================================================================
# Edit Menu Items
# ===========================================================================
    def set_temperature(self):
        # TODO: make it explict
        nowTemp = cst.Temperature
        newTemp, buttonResponse = QInputDialog.getDouble(
            self, 'ErwinJr Input Dialog', 'Set Temperature',
            value=nowTemp, min=0)
        if buttonResponse:
            cst.set_temperature(newTemp)
            self.quantumWidget.set_temperature(newTemp)


# ===========================================================================
# Help Menu Items
# ===========================================================================
    def on_about(self):
        msg = """ ErwinJr 3.x Authors and Contributors

         * Ming Lyu
            minglyu@princeton.edu

        ErwinJr 2.x Authors and Contributors

         * Kale J. Franz, PhD (Jet Propulsion Laboratory)
            kfranz@alumni.princeton.edu
            www.kalefranz.com

With contributions from:
         * Yamac Dikmelik (Johns Hopkins University)
         * Yu Song (Princeton University)
        """
        QMessageBox.about(self, "ErwinJr " + str(ejVersion), msg.strip())

    def on_licenses(self):
        copyright1 = """
#=======================================
# ErwinJr is a simulation program for quantum semiconductor lasers.
# Copyright (C) 2012 Kale J. Franz, PhD
# Copyright (C) 2017 Ming Lyu
#
# A portion of this code is Copyright (c) 2011, California Institute of
# Technology ("Caltech"). U.S. Government sponsorship acknowledged.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#=======================================
"""
        QMessageBox.about(self, "ErwinJr " + str(ejVersion),
                          copyright1.strip())

    def on_tutorial(self):
        if os.name == "nt":
            os.startfile("tutorial.pdf")
        elif os.name == "posix":
            os.system("/usr/bin/xdg-open tutorial.pdf")


def main():
    app = QApplication(sys.argv)
    app.setOrganizationName("JPL")
    app.setOrganizationDomain("erwinjr.org")
    app.setApplicationName("ErwinJr")
    qsettingsSystem = QSettings(QSettings.SystemScope, "JPL", "ErwinJr")
    installDirectory = str(qsettingsSystem.value(
        'installDirectory').toString())
    if installDirectory:
        os.chdir(installDirectory)

    # Create and display the splash screen
    splash_pix = QPixmap('images/erwinjr_splash.png')
    splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
    splash.setMask(splash_pix.mask())
    splash.show()
    app.processEvents()

    time.sleep(1)

    app.setWindowIcon(QIcon('images/EJpng48x48.png'))

    # this block handles a filename passed in by command line
    try:
        fileName = sys.argv[1]
        name, ext = os.path.splitext(fileName)
        assert ext == ".qcl"
        assert os.path.exists(fileName)
        fileName = os.path.abspath(fileName)
    except (IndexError, AssertionError):
        fileName = None

    form = MainWindow(fileName)
    form.show()
    splash.finish(form)

    qsettings = QSettings()
    if not qsettings.value('firstRun').toInt()[1]:
        if not installDirectory:
            qsettingsSystem.setValue("installDirectory", QVariant(os.getcwd()))
        firstRunBox = QMessageBox(
            QMessageBox.Question, 'EwrinJr '+str(majorVersion),
            ("Welcome to ErwinJr!\n"
             "Since this is your first time running the program, "
             "would you like to open an example file or a blank file?"),
            parent=form)
        firstRunBox.addButton("Blank File", QMessageBox.NoRole)
        firstRunBox.addButton("Example File", QMessageBox.YesRole)
        ansr =firstRunBox.exec_()
        if ansr:
            form.fileOpen('examples/NPhoton PQLiu.qcl')
        else:
            form.fileNew()
        qsettings.setValue("firstRun", 1)

    app.exec_()


if __name__ == "__main__":
    main()

# vim: ts=4 sw=4 sts=4 expandtab
