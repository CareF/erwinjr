#!/usr/bin/env python2
# -*- coding:utf-8 -*-

#===============================================================================
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
#===============================================================================

# TODO: 
# replace np.hstack (done, but self.strata part to test)
# find replacement for psyco
# try to seperate this file to smaller ones (done for quantumtab
# check unnecessary function call
# Ctrl+z support
# add status bar

from __future__ import division
import os, sys
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import PyQt4.Qwt5 as Qwt
import numpy as np
from numpy import pi, sqrt
from functools import partial
import time

import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D

import settings
import SupportClasses
from QCLayers import QCLayers, cst
from Strata import Strata
from QCLayers import h, c0, e0
import SaveLoad

from QuantumTab import QuantumTab

#============================================================================
# Version
#============================================================================
ejVersion = 171101
majorVersion = '3.2.0'

#============================================================================
# Debug options
#============================================================================
DEBUG = 1

class MainWindow(QMainWindow):
    def __init__(self, fileName=None, parent=None):
        super(MainWindow, self).__init__(parent)

        self.newLineChar = '\n'
        if sys.platform == "darwin":
            self.newLineChar = '\n'

        self.strata = Strata()

        self.stratumMaterialsList = ['Active Core', 
                                     'InP',
                                     'GaAs',
                                     'InGaAs', 
                                     'InAlAs', 
                                     'Au', 
                                     'SiNx', 
                                     'SiO2', 
                                     'Air']
        self.waveguideFacetsList = ['as-cleaved + as-cleaved',
                                    'as-cleaved + perfect HR',
                                    'as-cleaved + perfect AR',
                                    'perfect AR + perfect HR',
                                    'custom coating + as-cleaved',
                                    'custom coating + perfect HR',
                                    'custom coating + perfect AR']

        self.filename = fileName


        self.create_main_frame()
        self.create_Quantum_menu()

        self.stratumTable_refresh()

        self.update_stratum_inputBoxes()

        qsettings = QSettings()
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



#===============================================================================
# Create Main Frame        
#===============================================================================

    def create_main_frame(self):

        self.mainTabWidget = QTabWidget()

        # ##########################
        #
        # Optical Tab
        #
        # ##########################

        #vBox1
        vBox1Grid = QGridLayout()

        self.editOpticalParametersBox = QCheckBox('Edit Parameters')
        self.editOpticalParametersBox.setChecked(False)
        self.connect(self.editOpticalParametersBox, 
                SIGNAL('toggled(bool)'), 
                self.edit_optical_parameters)
        vBox1Grid.addWidget(self.editOpticalParametersBox, 
                0,0, 1,2, Qt.AlignCenter)

        wlLabel = QLabel('<b><center>Wavelength</center</b>')
        vBox1Grid.addWidget(wlLabel, 1,0, 1,2)
        self.wavelengthBox = QDoubleSpinBox()
        self.wavelengthBox.setValue(self.strata.wavelength)
        self.wavelengthBox.setSuffix(u' \u03BCm')
        self.wavelengthBox.setDecimals(3)
        self.wavelengthBox.setSingleStep(0.1)
        self.wavelengthBox.setRange(0.0,30.0)
        self.wavelengthBox.setReadOnly(True)
        self.wavelengthBox.setStyleSheet('color:gray')
        self.connect(self.wavelengthBox, SIGNAL("valueChanged(double)"), 
                self.input_wavelength)
        vBox1Grid.addWidget(self.wavelengthBox, 2,0, 1,2)

        vBox1Grid.addItem(QSpacerItem(20,20), 3,0, 1,2)

        vBox1Grid.addWidget(QLabel(
            '<b><center>Operating<br>Field</center></b>'), 4,0, 1,1)
        self.operatingFieldBox = QDoubleSpinBox()
        self.operatingFieldBox.setDecimals(1)
        self.operatingFieldBox.setRange(0.,300.)
        self.operatingFieldBox.setSingleStep(1)
        self.operatingFieldBox.setSuffix(u' kV/cm')
        self.operatingFieldBox.setReadOnly(True)
        self.operatingFieldBox.setStyleSheet('color:gray')
        self.connect(self.operatingFieldBox, SIGNAL("valueChanged(double)"), 
                self.input_operatingField)
        vBox1Grid.addWidget(self.operatingFieldBox, 5,0, 1,1)

        vBox1Grid.addWidget(QLabel(
            '<b><center>Active Core<br>Period Length</center></b>'), 4,1, 1,1)
        self.ACPeriodLengthBox = QDoubleSpinBox()
        self.ACPeriodLengthBox.setDecimals(1)
        self.ACPeriodLengthBox.setRange(0.,10000.)
        self.ACPeriodLengthBox.setSingleStep(1)
        self.ACPeriodLengthBox.setSuffix(u' \u212B')
        self.connect(self.ACPeriodLengthBox, SIGNAL("valueChanged(double)"), 
                self.input_ACPeriodLength)
        self.ACPeriodLengthBox.setReadOnly(True)
        self.ACPeriodLengthBox.setStyleSheet('color:gray')
        vBox1Grid.addWidget(self.ACPeriodLengthBox, 5,1, 1,1)

        vBox1Grid.addWidget(QLabel(
            '<b><center>Active Core<br>Periods</center></b>'), 6,0, 1,1)
        self.ACPeriodsBox = QSpinBox()
        self.ACPeriodsBox.setRange(1,99)
        self.connect(self.ACPeriodsBox, SIGNAL("valueChanged(int)"), 
                self.input_ACPeriods)
        vBox1Grid.addWidget(self.ACPeriodsBox, 7,0, 1,1)

        vBox1Grid.addWidget(QLabel(
            '<b><center>Operating<br>Voltage</center></b>'), 6,1, 1,1)
        self.OperatingVoltageBox = QLineEdit()
        self.OperatingVoltageBox.setMaximumWidth(150)
        self.OperatingVoltageBox.setReadOnly(True)
        self.OperatingVoltageBox.setStyleSheet('color:gray')
        vBox1Grid.addWidget(self.OperatingVoltageBox, 7,1, 1,1)

        vBox1Grid.addItem(QSpacerItem(20,20), 8,0, 1,2)

        vBox1Grid.addWidget(QLabel(
            u'<b><center><i>\u03B1<sub>core</sub></i></center></b>'), 9,0, 1,1)
        self.aCoreBox = QLineEdit()
        self.aCoreBox.setMaximumWidth(150)
        self.aCoreBox.setReadOnly(True)
        self.aCoreBox.setStyleSheet('color:gray')
        self.connect(self.aCoreBox, SIGNAL("editingFinished()"), 
                self.input_aCore)
        vBox1Grid.addWidget(self.aCoreBox, 10,0, 1,1)

        vBox1Grid.addWidget(QLabel(
            u'<center><b>\u00F1<i><sub>core</sub></i></b></center>'), 9,1, 1,1)
        self.nCoreBox = QLineEdit()
        self.nCoreBox.setMaximumWidth(150)
        self.nCoreBox.setReadOnly(True)
        self.nCoreBox.setStyleSheet('color:gray')
        vBox1Grid.addWidget(self.nCoreBox, 10,1, 1,1)

        vBox1Grid.addItem(QSpacerItem(20,20), 11,0, 1,2)

        vBox1Grid.addWidget(QLabel(
            u'<center><b>Transition Broadening</b></center>'), 12,0, 1,2)
        self.transitionBroadeningBox = QDoubleSpinBox()
        self.transitionBroadeningBox.setMaximumWidth(85)
        self.transitionBroadeningBox.setDecimals(1)
        self.transitionBroadeningBox.setRange(0.,1000.)
        self.transitionBroadeningBox.setSingleStep(1)
        self.transitionBroadeningBox.setSuffix(u' meV')
        self.transitionBroadeningBox.setReadOnly(True)
        self.transitionBroadeningBox.setStyleSheet('color:gray')
        self.connect(self.transitionBroadeningBox, 
                SIGNAL("valueChanged(double)"), 
                self.input_transitionBroadening)
        vBox1Grid.addWidget(self.transitionBroadeningBox, 13,0, 1,2, Qt.AlignCenter)

        vBox1subGrid1 = QGridLayout()

        vBox1subGrid1.addWidget(QLabel(
            u'<b><center><i>\u03C4<sub>upper</sub></i></center></b>'), 
            0,0, 1,1)
        self.tauUpperBox = QDoubleSpinBox()
        self.tauUpperBox.setDecimals(3)
        self.tauUpperBox.setRange(0.,99.)
        self.tauUpperBox.setSingleStep(1)
        self.tauUpperBox.setSuffix(u' ps')
        self.tauUpperBox.setReadOnly(True)
        self.tauUpperBox.setStyleSheet('color:gray')
        #self.tauUpperBox.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed))
        self.connect(self.tauUpperBox, SIGNAL("valueChanged(double)"), 
                self.input_tauUpper)
        vBox1subGrid1.addWidget(self.tauUpperBox, 1,0, 1,1)

        vBox1subGrid1.addWidget(QLabel(
            u'<b><center><i>\u03C4<sub>lower</sub></i></b></center></b>'), 
            0,1, 1,1)
        self.tauLowerBox = QDoubleSpinBox()
        self.tauLowerBox.setDecimals(3)
        self.tauLowerBox.setRange(0.,99.)
        self.tauLowerBox.setSingleStep(1)
        self.tauLowerBox.setSuffix(u' ps')
        self.tauLowerBox.setReadOnly(True)
        self.tauLowerBox.setStyleSheet('color:gray')
        #self.tauLowerBox.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed))
        self.connect(self.tauLowerBox, SIGNAL("valueChanged(double)"), 
                self.input_tauLower)
        vBox1subGrid1.addWidget(self.tauLowerBox, 1,1, 1,1)

        vBox1subGrid1.addWidget(QLabel(u'<b><center><i>\u03C4<sub>upper,lower</sub></i></b></center></b>'), 0,2, 1,1)
        self.tauUpperLowerBox = QDoubleSpinBox()
        self.tauUpperLowerBox.setDecimals(3)
        self.tauUpperLowerBox.setRange(0.,99.)
        self.tauUpperLowerBox.setSingleStep(1)
        self.tauUpperLowerBox.setSuffix(u' ps')
        self.tauUpperLowerBox.setReadOnly(True)
        self.tauUpperLowerBox.setStyleSheet('color:gray')
        #self.tauUpperLowerBox.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed))
        self.connect(self.tauUpperLowerBox, SIGNAL("valueChanged(double)"), 
                self.input_tauUpperLower)
        vBox1subGrid1.addWidget(self.tauUpperLowerBox, 1,2, 1,1)

        vBox1Grid.addLayout(vBox1subGrid1, 14,0, 1,2)

        vBox1Grid.addWidget(QLabel(
            u'<b><center>optical dipole</center></b>'), 15,0, 1,1)
        self.opticalDipoleBox = QDoubleSpinBox()
        self.opticalDipoleBox.setDecimals(1)
        self.opticalDipoleBox.setRange(0.,10000.)
        self.opticalDipoleBox.setSingleStep(1)
        self.opticalDipoleBox.setSuffix(u' \u212B')
        self.opticalDipoleBox.setReadOnly(True)
        self.opticalDipoleBox.setStyleSheet('color:gray')
        self.connect(self.opticalDipoleBox, 
                SIGNAL("valueChanged(double)"), 
                self.input_opticalDipole)
        vBox1Grid.addWidget(self.opticalDipoleBox, 16,0, 1,1)

        vBox1Grid.addWidget(QLabel(
            u'<center><b>Figure of Merit</b></center>'), 15,1, 1,1)
        self.FoMBox = QLineEdit()
        self.FoMBox.setMaximumWidth(150)
        self.FoMBox.setReadOnly(True)
        self.FoMBox.setStyleSheet('color:gray')
        vBox1Grid.addWidget(self.FoMBox, 16,1, 1,1)

        vBox1Grid.addItem(QSpacerItem(20,20), 17,0, 1,2)

        vBox1Grid.addWidget(QLabel(
            u'<center><b>Waveguide Facets</b></center>'), 18,0, 1,2)
        self.waveguideFacetsBox = QComboBox();
        self.waveguideFacetsBox.addItems(self.waveguideFacetsList)
        self.connect(self.waveguideFacetsBox, 
                SIGNAL("currentIndexChanged(const QString &)"), 
                self.input_waveguideFacets)
        vBox1Grid.addWidget(self.waveguideFacetsBox, 19,0, 1,2)

        vBox1Grid.addWidget(QLabel(
            u'<b><center>Waveguide<br>Length</center></b>'), 20,0, 1,1)
        self.waveguideLengthBox = QDoubleSpinBox()
        self.waveguideLengthBox.setDecimals(1)
        self.waveguideLengthBox.setRange(0.,20.)
        self.waveguideLengthBox.setSingleStep(1)
        self.waveguideLengthBox.setSuffix(u' mm')
        self.connect(self.waveguideLengthBox, 
                SIGNAL("valueChanged(double)"), 
                self.input_waveguideLength)
        vBox1Grid.addWidget(self.waveguideLengthBox, 21,0, 1,1)


        self.customFacetBoxLabel = QLabel(
                u'<b><center>Custom<br>Reflectivity</center></b>')
        self.customFacetBoxLabel.setStyleSheet('color:gray')
        vBox1Grid.addWidget(self.customFacetBoxLabel, 20,1, 1,1)
        self.customFacetBox = QDoubleSpinBox()
        self.customFacetBox.setDecimals(1)
        self.customFacetBox.setRange(0.,100.)
        self.customFacetBox.setSingleStep(1)
        self.customFacetBox.setSuffix(u'%')
        self.customFacetBox.setEnabled(False)
        self.connect(self.customFacetBox, SIGNAL("valueChanged(double)"), 
                self.input_customFacet)
        vBox1Grid.addWidget(self.customFacetBox, 21,1, 1,1)


        vBox1GridWidget = QWidget()
        vBox1GridWidget.setLayout(vBox1Grid)
        vBox1GridWidget.setContentsMargins(0,0,0,0)
        vBox1GridWidget.setMaximumWidth(235)
        vBox1 = QVBoxLayout()
        vBox1.addWidget(vBox1GridWidget)
        vBox1.setSpacing(0)
        vBox1.addStretch()


        #set up stratumTable
        self.stratumTable = QTableWidget()
        self.stratumTable.setSelectionBehavior(QTableWidget.SelectRows)
        self.stratumTable.setSelectionMode(QTableWidget.SingleSelection)
        self.stratumTable.setMinimumWidth(450)
        self.stratumTable.setMaximumWidth(450)
        if sys.platform == 'darwin':
            self.stratumTable.setMinimumWidth(550)
            self.stratumTable.setMaximumWidth(550)
        self.stratumTable.setMinimumHeight(450)
        self.stratumTable.setMaximumHeight(650)
        self.connect(self.stratumTable, 
                SIGNAL("itemChanged(QTableWidgetItem*)"),
                self.stratumTable_itemChanged)
        self.connect(self.stratumTable, SIGNAL("itemSelectionChanged()"),
                self.stratumTable_itemSelectionChanged)


        insertStratumAboveButton = QPushButton("Insert Stratum Above")
        self.connect(insertStratumAboveButton, SIGNAL("clicked()"), 
                self.insert_stratumAbove)
        insertStratumBelowButton = QPushButton("Insert Stratum Below")
        self.connect(insertStratumBelowButton, SIGNAL("clicked()"), 
                self.insert_stratumBelow)
        deleteStratumButton = QPushButton("Delete Stratum")
        self.connect(deleteStratumButton, SIGNAL("clicked()"), self.delete_stratum)

        #vBox2
        vBox2Grid = QGridLayout()
        vBox2Grid.addWidget(insertStratumAboveButton, 1,0, 1,1)
        vBox2Grid.addWidget(insertStratumBelowButton, 1,1, 1,1)
        vBox2Grid.addWidget(deleteStratumButton, 1,2, 1,1)
        vBox2Grid.addWidget(self.stratumTable, 2,0, 1,3)

        #Optimization
        optiFrameLayout = QGridLayout()

        self.opti1DChoiceBox = QComboBox()
        self.opti1DChoiceBox.addItems(['Thickness','Doping'])
        self.connect(self.opti1DChoiceBox, 
                SIGNAL("currentIndexChanged(const QString &)"), 
                self.input_opti1DChoice)
        optiFrameLayout.addWidget(self.opti1DChoiceBox, 1,0, 1,1)

        opti1DLayerBoxLabel = QLabel('<b><center>1<sup>st</sup> Dimension<br>'
                'Strata Number(s)</center></b>')
        opti1DLayerBoxLabel.setToolTip('Ex: 6,8')
        optiFrameLayout.addWidget(opti1DLayerBoxLabel, 0,1, 1,1)
        self.opti1DLayerBox = QLineEdit()
        self.connect(self.opti1DLayerBox, SIGNAL('editingFinished()'), 
                self.input_opti1DLayer)
        self.opti1DLayerBox.setToolTip('Ex: 6,8')
        optiFrameLayout.addWidget(self.opti1DLayerBox, 1,1, 1,1)

        opti1DRangeBoxLabel = QLabel('<b><center>Optimization<br>Range</center></b>')
        opti1DRangeBoxLabel.setToolTip('Ex: 1:0.1:3')
        optiFrameLayout.addWidget(opti1DRangeBoxLabel, 0,2, 1,1)
        self.opti1DRangeBox = QLineEdit()
        self.opti1DRangeBox.setToolTip('Ex: 1:0.1:3')
        self.connect(self.opti1DRangeBox, SIGNAL('editingFinished()'), 
                self.input_opti1DRange)
        optiFrameLayout.addWidget(self.opti1DRangeBox, 1,2, 1,1)

        self.opti1DRunButton = QPushButton('Optimize 1D')
        self.connect(self.opti1DRunButton, SIGNAL('clicked(bool)'), self.run_opti1D)
        optiFrameLayout.addWidget(self.opti1DRunButton, 1,3, 1,1)

        self.opti2DChoiceBox = QComboBox()
        self.opti2DChoiceBox.addItems(['Thickness','Doping'])
        self.connect(self.opti2DChoiceBox, 
                SIGNAL("currentIndexChanged(const QString &)"), 
                self.input_opti2DChoice)
        optiFrameLayout.addWidget(self.opti2DChoiceBox, 3,0, 1,1)

        opti2DLayerBoxLabel = QLabel('<b><center>2<sup>nd</sup> Dimension<br>'
                'Strata Number(s)</center></b>')
        opti2DLayerBoxLabel.setToolTip('Ex: 2 5 7')
        optiFrameLayout.addWidget(opti2DLayerBoxLabel, 2,1, 1,1)
        self.opti2DLayerBox = QLineEdit()
        self.opti2DLayerBox.setToolTip('Ex: 2 5 7')
        self.connect(self.opti2DLayerBox, SIGNAL('editingFinished()'), 
                self.input_opti2DLayer)
        optiFrameLayout.addWidget(self.opti2DLayerBox, 3,1, 1,1)

        opti2DRangeBoxLabel = QLabel('<b><center>Optimization<br>Range</center></b>')
        opti2DRangeBoxLabel.setToolTip('Ex: 1:5')
        optiFrameLayout.addWidget(opti2DRangeBoxLabel, 2,2, 1,1)
        self.opti2DRangeBox = QLineEdit()
        self.opti2DRangeBox.setToolTip('Ex: 1:5')
        self.connect(self.opti2DRangeBox, SIGNAL('editingFinished()'), 
                self.input_opti2DRange)
        optiFrameLayout.addWidget(self.opti2DRangeBox, 3,2, 1,1)

        self.opti2DRunButton = QPushButton('Optimize 2D')
        self.connect(self.opti2DRunButton, SIGNAL('clicked(bool)'), self.run_opti2D)
        optiFrameLayout.addWidget(self.opti2DRunButton, 3,3, 1,1)

        self.optiFrame = QGroupBox('Optimization')
        self.optiFrame.setMaximumWidth(450)
        self.optiFrame.setMinimumWidth(450)
        self.optiFrame.setLayout(optiFrameLayout)
        vBox2Grid.addWidget(self.optiFrame, 3,0, 1,3)

        vBox2 = QVBoxLayout()
        vBox2.addLayout(vBox2Grid)
        vBox2.addStretch()



        #vBox3
        self.plotModeButton = QPushButton("Plot Mode")
        self.connect(self.plotModeButton, SIGNAL("clicked()"), self.solve_mode)
        vBox3 = QVBoxLayout()
        vBox3.addWidget(self.plotModeButton)

        vBox3.addWidget(QLabel(u'<center><b><i>\u03B2<sub>eff</sub></i></b></center>'))
        self.betaEffBox = QLineEdit()
        self.betaEffBox.setMaximumWidth(150)
        self.betaEffBox.setEnabled(False)
        vBox3.addWidget(self.betaEffBox)

        self.modeCalculationsBox = QTextEdit('')
        self.modeCalculationsBox.setReadOnly(True)
        self.modeCalculationsBox.setSizePolicy(QSizePolicy(
            QSizePolicy.Fixed,QSizePolicy.Fixed))
        self.modeCalculationsBox.setMaximumHeight(175)
        self.modeCalculationsBox.setMaximumWidth(150)
        vBox3.addWidget(self.modeCalculationsBox)

        vBox3.addStretch()


        #vBox4

        #set up opticalCanvas for stratum / mode plot
        self.opticalCanvas = Qwt.QwtPlot(self)
        self.opticalCanvas.setCanvasBackground(Qt.white)
        self.opticalCanvas.canvas().setCursor(Qt.ArrowCursor)

        #optical optimization canvas
        self.optimization1DCanvas = Qwt.QwtPlot(self)
        self.optimization1DCanvas.setCanvasBackground(Qt.white)
        self.optimization1DCanvas.canvas().setCursor(Qt.ArrowCursor)
        self.optimization1DCanvas.setVisible(False)

        #2D optical optimization canvas
        #optimization2DFig = Figure((5.0, 4.0), dpi=dpi)
        self.optimization2DFig = Figure()
        self.optimization2DCanvas = FigureCanvas(self.optimization2DFig)
        #self.optimization2DAxes = self.optimization2DFig.add_subplot(111, projection='3d')
        margins = [0.05,0.05,0.95,0.95]
        self.optimization2DAxes = self.optimization2DFig.add_axes(
                margins, projection='3d')
        self.optimization2DAxes.autoscale(enable=True, axis='both', tight=True)
        #get the background color of the central widget
        #bgColor = self.mainTabWidget.palette().brush(QPalette.Window).color().name()
        bgColorRed = self.mainTabWidget.palette().brush(QPalette.Window).color().red()
        bgColorBlue = self.mainTabWidget.palette().brush(QPalette.Window).color().blue()
        bgColorGreen = self.mainTabWidget.palette().brush(QPalette.Window).color().green()
        self.bgColor = (bgColorRed/255.0, bgColorGreen/255.0, bgColorBlue/255.0, 1)
        self.optimization2DAxes.patch.set_color(self.bgColor)
        self.optimization2DFig.patch.set_color(self.bgColor)
        self.optimization2DCanvas.setVisible(False)



        vBox4 = QVBoxLayout()
        vBox4.addWidget(self.opticalCanvas)
        vBox4.addWidget(self.optimization1DCanvas)
        vBox4.addWidget(self.optimization2DCanvas)

        opticalLayout = QHBoxLayout()
        opticalLayout.addLayout(vBox1)
        opticalLayout.addLayout(vBox2)
        opticalLayout.addLayout(vBox3)
        opticalLayout.addLayout(vBox4)  



        self.opticalWidget = QWidget()
        self.opticalWidget.setLayout(opticalLayout)
        self.opticalWidget.setAutoFillBackground(True)
        self.opticalWidget.setBackgroundRole(QPalette.Window)        

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
#        self.connect(thermalTable,SIGNAL("itemChanged(QTableWidgetItem*)"),self.stratumTable_itemChanged)
#        self.connect(thermalTable,SIGNAL("itemSelectionChanged()"),self.stratumTable_itemSelectionChanged)
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
        self.connect(self.quantumWidget, 
                SIGNAL('dirty'), 
                self.quantumUpdate)
        self.connect(self.quantumWidget.toOpticalParamsButton, 
                SIGNAL("clicked()"), 
                self.transfer_optical_parameters)
        self.mainTabWidget.addTab(self.quantumWidget, 'Quantum')

        self.mainTabWidget.addTab(self.opticalWidget, 'Optical')
#        self.mainTabWidget.addTab(thermalWidget, 'Thermal')
        self.connect(self.mainTabWidget, SIGNAL('currentChanged(int)'), 
                self.change_main_tab)

        self.setCentralWidget(self.mainTabWidget)

        self.quantumWidget.layerTable.selectRow(0)
        self.quantumWidget.layerTable.setFocus()




#===============================================================================
# Optical Tab Input Controls
#===============================================================================

    def update_stratum_inputBoxes(self):
        self.wavelengthBox.setValue(self.strata.wavelength)
        self.operatingFieldBox.setValue(self.strata.operatingField)
        self.ACPeriodLengthBox.setValue(self.strata.Lp)
        self.ACPeriodsBox.setValue(self.strata.Np)
        self.OperatingVoltageBox.setText('{0:.1f} V'.format(
                    self.strata.Np*self.strata.operatingField/self.strata.Lp))
        self.aCoreBox.setText(
                '{0:.3f} cm^-1'.format(self.strata.aCore))
        self.nCoreBox.setText(
                '{0.real:2.3f}+{0.imag:1.3e}j'.format(self.strata.nCore))
        self.transitionBroadeningBox.setValue(
                self.strata.transitionBroadening * 1000) #display in meV
        self.tauUpperBox.setValue(self.strata.tauUpper)
        self.tauLowerBox.setValue(self.strata.tauLower)
        self.tauUpperLowerBox.setValue(self.strata.tauUpperLower)
        self.opticalDipoleBox.setValue(self.strata.opticalDipole)
        self.FoMBox.setText(u'{0:4.0f} ps \u212B^2'.format(self.strata.FoM))
        self.waveguideFacetsBox.setCurrentIndex(
                self.waveguideFacetsList.index(self.strata.waveguideFacets))
        self.waveguideLengthBox.setValue(self.strata.waveguideLength)
        self.customFacetBox.setValue(self.strata.customFacet * 100) #display in percent

        self.strata.updateFacets()

    def update_modeCalculations_box(self):

        self.strata.calculate_performance_parameters()
        reportString = (u"\u0393: <b>%3.1f%%</b><br>"
                u"<i>\u03B1<sub>wg</sub></i> : %3.1f cm<sup>-1</sup><br>"
                u"<i>\u03B1<sub>mirror</sub></i> : %3.1f cm<sup>-1</sup><br>"
                u"gain: %3.3f cm/A<br>"
                u"<i>J<sub>th0</sub></i> : <b>%3.3f kA/cm<sup>2</sup></b><br>"
                u"<i>I<sub>th0</sub></i> : %5.1f mA<br>"
                u"<i>V<sub>op</sub></i> : %3.1f V<br>"
                u"<i>\u03B7<sub>voltage</sub></i> : %3.1f%%<br>"
                u"<i>\u03B7<sub>extraction</sub></i> : %3.1f%%<br>"
                u"<i>\u03B7<sub>inversion</sub></i> : %3.1f%%<br>"
                u"<i>\u03B7<sub>modal</sub></i> : %3.1f%%<br>")%(
                        self.strata.confinementFactor*100, 
                        self.strata.waveguideLoss, 
                        self.strata.mirrorLoss,
                        self.strata.gain, 
                        self.strata.Jth0, 
                        self.strata.Ith0*1000, 
                        self.strata.operatingVoltage, 
                        self.strata.voltageEfficiency*100,
                        self.strata.extractionEfficiency*100,
                        self.strata.inversionEfficiency*100, 
                        self.strata.modalEfficiency*100)

        #  reportString = u""

        #  #confinement factor
        #  reportString += u"\u0393: <b>{0:3.1f}%</b><br>".format(self.strata.confinementFactor*100)
        #  #waveguide loss
        #  reportString += u"<i>\u03B1<sub>wg</sub></i> : {0:3.1f} cm<sup>-1</sup><br>".format(self.strata.waveguideLoss)
        #  #mirror loss
        #  reportString += u"<i>\u03B1<sub>mirror</sub></i> : {0:3.1f} cm<sup>-1</sup><br>".format(self.strata.mirrorLoss)
        #  #gain
        #  reportString += u"gain: {0:3.3f} cm/A<br>".format(self.strata.gain)
        #  #Jth0
        #  reportString += u"<i>J<sub>th0</sub></i> : <b>{0:3.3f} kA/cm<sup>2</sup></b><br>".format(self.strata.Jth0)
        #  #Ith0
        #  reportString += u"<i>I<sub>th0</sub></i> : {0:5.1f} mA<br>".format(self.strata.Ith0*1000)

        #  #Voltage
        #  reportString += u"<i>V<sub>op</sub></i> : {0:3.1f} V<br>".format(self.strata.operatingVoltage)
        #  #Voltage Efficiency
        #  reportString += u"<i>\u03B7<sub>voltage</sub></i> : {0:3.1f}%<br>".format(self.strata.voltageEfficiency*100)
        #  #Extraction Efficiency
        #  reportString += u"<i>\u03B7<sub>extraction</sub></i> : {0:3.1f}%<br>".format(self.strata.extractionEfficiency*100)
        #  #Inversion Efficiency
        #  reportString += u"<i>\u03B7<sub>inversion</sub></i> : {0:3.1f}%<br>".format(self.strata.inversionEfficiency*100)
        #  #Modal Efficiency
        #  reportString += u"<i>\u03B7<sub>modal</sub></i> : {0:3.1f}%<br>".format(self.strata.modalEfficiency*100)

        self.modeCalculationsBox.setText(reportString)

    def edit_optical_parameters(self, toggleState):
        if toggleState == True:
            self.wavelengthBox.setReadOnly(False)
            self.wavelengthBox.setStyleSheet('color:black')
            self.operatingFieldBox.setReadOnly(False)
            self.operatingFieldBox.setStyleSheet('color:black')
            self.ACPeriodLengthBox.setReadOnly(False)
            self.ACPeriodLengthBox.setStyleSheet('color:black')
            self.aCoreBox.setReadOnly(False)
            self.aCoreBox.setStyleSheet('color:black')
            self.tauUpperBox.setReadOnly(False)
            self.tauUpperBox.setStyleSheet('color:black')
            self.tauUpperLowerBox.setReadOnly(False)
            self.tauUpperLowerBox.setStyleSheet('color:black')
            self.tauLowerBox.setReadOnly(False)
            self.tauLowerBox.setStyleSheet('color:black')
            self.opticalDipoleBox.setReadOnly(False)
            self.opticalDipoleBox.setStyleSheet('color:black')
            self.transitionBroadeningBox.setReadOnly(False)
            self.transitionBroadeningBox.setStyleSheet('color:black')
        else:
            self.wavelengthBox.setReadOnly(True)
            self.wavelengthBox.setStyleSheet('color:gray')
            self.operatingFieldBox.setReadOnly(True)
            self.operatingFieldBox.setStyleSheet('color:gray')
            self.ACPeriodLengthBox.setReadOnly(True)
            self.ACPeriodLengthBox.setStyleSheet('color:gray')
            self.aCoreBox.setReadOnly(True)
            self.aCoreBox.setStyleSheet('color:gray')
            self.tauUpperBox.setReadOnly(True)
            self.tauUpperBox.setStyleSheet('color:gray')
            self.tauLowerBox.setReadOnly(True)
            self.tauLowerBox.setStyleSheet('color:gray')
            self.tauUpperLowerBox.setReadOnly(True)
            self.tauUpperLowerBox.setStyleSheet('color:gray')
            self.opticalDipoleBox.setReadOnly(True)
            self.opticalDipoleBox.setStyleSheet('color:gray')
            self.transitionBroadeningBox.setReadOnly(True)
            self.transitionBroadeningBox.setStyleSheet('color:gray')

    def input_wavelength(self, value):
        self.strata.wavelength = value

        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()

    def input_operatingField(self, value):
        self.strata.operatingField = value

        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()

    def input_ACPeriodLength(self, value):
        self.strata.Lp = value

        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()

    def input_ACPeriods(self, value):
        self.strata.Np = value

        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()

    def input_aCore(self):
        initialText = unicode(self.aCoreBox.text())
        txt = initialText.split()[0]
        try:
            value = float(txt)
            self.strata.aCore = value
            kCore = 1/(4*pi) * self.strata.aCore * self.strata.wavelength*1e-4 
            # See Def of acore
            # 1e-4: aCore in cm-1, wl in um
            self.strata.nCore = self.quantumWidget.get_nCore(
                    self.strata.wavelength) + 1j*kCore


            self.dirty = True
            self.update_windowTitle()       
            self.update_stratum_inputBoxes()
            self.stratumTable_refresh()
        except ValueError:
            self.aCore.setText(initialText)            

    def input_transitionBroadening(self, value):
        self.strata.transitionBroadening = value / 1000

        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()

    def input_tauUpper(self, value):
        self.strata.tauUpper = value
        self.strata.FoM = self.strata.opticalDipole**2 * self.strata.tauUpper * (1- self.strata.tauLower/self.strata.tauUpperLower)

        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()

    def input_tauLower(self, value):
        self.strata.tauLower = value
        self.strata.FoM = self.strata.opticalDipole**2 * self.strata.tauUpper * (1- self.strata.tauLower/self.strata.tauUpperLower)

        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()

    def input_tauUpperLower(self, value):
        self.strata.tauUpperLower = value
        self.strata.FoM = self.strata.opticalDipole**2 * self.strata.tauUpper * (1- self.strata.tauLower/self.strata.tauUpperLower)

        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()

    def input_opticalDipole(self, value):
        self.strata.opticalDipole = value
        self.strata.FoM = self.strata.opticalDipole**2 * self.strata.tauUpper * (1- self.strata.tauLower/self.strata.tauUpperLower)

        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()

    def input_waveguideFacets(self, selection):
        self.strata.waveguideFacets = selection
        if selection == 'as-cleaved + as-cleaved':
            self.customFacetBoxLabel.setStyleSheet('color:gray')
            self.customFacetBox.setEnabled(False)
        elif selection == 'as-cleaved + perfect HR':
            self.customFacetBoxLabel.setStyleSheet('color:gray')
            self.customFacetBox.setEnabled(False)
        elif selection == 'as-cleaved + perfect AR':
            self.customFacetBoxLabel.setStyleSheet('color:gray')
            self.customFacetBox.setEnabled(False)
        elif selection == 'perfect AR + perfect HR':
            self.customFacetBoxLabel.setStyleSheet('color:gray')
            self.customFacetBox.setEnabled(False)
        elif selection == 'custom coating + as-cleaved':
            self.customFacetBoxLabel.setStyleSheet('color:black')
            self.customFacetBox.setEnabled(True)
        elif selection == 'custom coating + perfect HR':
            self.customFacetBoxLabel.setStyleSheet('color:black')
            self.customFacetBox.setEnabled(True)
        elif selection == 'custom coating + perfect AR':
            self.customFacetBoxLabel.setStyleSheet('color:black')
            self.customFacetBox.setEnabled(True)

        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()

    def input_waveguideLength(self, value):
        self.strata.waveguideLength = value

        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()

    def input_customFacet(self, value):
        self.strata.customFacet = value / 100.0

        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()

    def input_opti1DChoice(self, selectionString):
        pass

    def input_opti1DLayer(self):
        pass

    def input_opti1DRange(self):
        pass

    def input_opti2DChoice(self, selectionString):
        pass

    def input_opti2DLayer(self):
        pass

    def input_opti2DRange(self):
        pass




#===============================================================================
# Optical Tab Strata Table Control
#===============================================================================

    def stratumTable_refresh(self):
        #calculate index for each layer
        self.strata.populate_rIndexes()

        #set properties for Active Core Layer
        for q in xrange(self.strata.stratumDopings.size):
            if self.strata.stratumMaterials[q] == 'Active Core':
                self.strata.stratumThicknesses[q] = \
                        self.strata.Np * self.strata.Lp * 1e-4
                self.strata.stratumDopings[q] = self.strata.nD
                self.strata.stratumRIndexes[q] = self.strata.nCore

        #update table
        self.stratumTable.clear()
        self.stratumTable.setColumnCount(6)
        self.stratumTable.setRowCount(self.strata.stratumDopings.size) #need to change
        self.stratumTable.setHorizontalHeaderLabels(
                ['Material', 'Mole Frac', 'Thickness', 'Doping', 
                    'Refractive Index', 'Loss'])
        vertLabels = []
        for n in xrange(self.strata.stratumDopings.size):
            vertLabels.append(str(n+1))
        self.stratumTable.setVerticalHeaderLabels(vertLabels)

        for q in xrange(self.strata.stratumDopings.size):
            #Stratum Material Setup
            materialWidget = QComboBox()
            materialWidget.addItems(self.stratumMaterialsList)
            materialWidget.setCurrentIndex(
                    self.stratumMaterialsList.index(self.strata.stratumMaterials[q]))
            self.connect(materialWidget, SIGNAL("currentIndexChanged(int)"), 
                    partial(self.stratumTable_materialChanged, q))
            self.stratumTable.setCellWidget(q, 0, materialWidget)

            #Stratum Composition Setup
            composition = QTableWidgetItem()
            if self.strata.stratumMaterials[q] not in self.strata.needsCompositionList:
                composition.setFlags(Qt.NoItemFlags)
            else:
                composition.setData(0,'{0:3.3f}'.format(
                    self.strata.stratumCompositions[q]))
                composition.setTextAlignment(Qt.AlignCenter)
            self.stratumTable.setItem(q, 1, composition)

            #Stratum Thickness Setup
            thickness = QTableWidgetItem(unicode(self.strata.stratumThicknesses[q]))
            thickness.setTextAlignment(Qt.AlignCenter)
            self.stratumTable.setItem(q, 2, thickness)
            if self.strata.stratumMaterials[q] == 'Active Core':
                thickness.setFlags(Qt.NoItemFlags)

            #Stratum Doping Setup
            doping = QTableWidgetItem()
            if self.strata.stratumMaterials[q] in self.strata.notDopableList:
                doping.setFlags(Qt.NoItemFlags)
            else:
                doping.setData(0,'{0:3.2f}'.format(self.strata.stratumDopings[q]))
                doping.setTextAlignment(Qt.AlignCenter)
                if self.strata.stratumMaterials[q] == 'Active Core':
                    doping.setFlags(Qt.NoItemFlags)
            self.stratumTable.setItem(q, 3, doping)

            #Stratum RIndex Setup
            rIndex = QTableWidgetItem(
                    '{0.real:2.3f}+{0.imag:1.3e}j'.format(
                        self.strata.stratumRIndexes[q]))
            rIndex.setTextAlignment(Qt.AlignCenter)
            rIndex.setFlags(Qt.NoItemFlags)
            self.stratumTable.setItem(q, 4, rIndex)

            #Stratum Loss Setup
            loss = self.strata.stratumRIndexes[q].imag*4*pi/self.strata.wavelength/1e-4
            alpha = QTableWidgetItem('{0:3.2f}'.format(loss))
            alpha.setTextAlignment(Qt.AlignCenter)
            alpha.setFlags(Qt.NoItemFlags)
            self.stratumTable.setItem(q, 5, alpha)

        self.stratumTable.resizeColumnsToContents()

        self.update_opticalCanvas()

    def stratumTable_itemChanged(self, item):
        column = self.stratumTable.currentColumn()
        row = self.stratumTable.currentRow()
        if column == -1: #column == -1 on GUI initialization
            return
        elif column == 0:
            return
        elif column == 1:
            xFrac = float(item.text())
            if xFrac < 0 or xFrac > 1:
                QMessageBox.warning(self,
                        'ErwinJr Error',
                        'Mole Fraction must be between 0 and 1')
            else:
                self.strata.stratumCompositions[row] = xFrac
        elif column == 2:
            self.strata.stratumThicknesses[row] = float(item.text())
        elif column == 3:
            self.strata.stratumDopings[row] = float(item.text())

        self.stratumTable_refresh()
        self.stratumTable.selectRow(row)

        self.dirty = True
        self.update_windowTitle()   

    def stratumTable_itemSelectionChanged(self):
        self.strata.stratumSelected = self.stratumTable.currentRow()
        if self.strata.stratumSelected >= 0 and \
                self.strata.stratumSelected < self.quantumWidget.layerNum():
            self.strata.populate_x()
            self.update_opticalCanvas()

    def stratumTable_materialChanged(self, row, selection):
        self.strata.stratumMaterials[row] = self.stratumMaterialsList[selection]

        self.stratumTable_refresh()
        self.stratumTable.selectRow(row)

        self.dirty = True
        self.update_windowTitle()  

    def insert_stratumAbove(self):
        row = self.stratumTable.currentRow()
        if row == -1:
            return

        #  if row == 0:
            #  self.strata.stratumMaterials.insert(row, self.strata.stratumMaterials[row])
            #  self.strata.stratumCompositions = np.hstack([self.strata.stratumCompositions[row], self.strata.stratumCompositions[row:,]])
            #  self.strata.stratumThicknesses = np.hstack([self.strata.stratumThicknesses[row], self.strata.stratumThicknesses[row:,]])
            #  self.strata.stratumDopings = np.hstack([self.strata.stratumDopings[row], self.strata.stratumDopings[row:,]])
            #  self.strata.stratumRIndexes = np.hstack([self.strata.stratumRIndexes[row], self.strata.stratumRIndexes[row:,]])

        #  else:
            #  self.strata.stratumMaterials.insert(row, self.strata.stratumMaterials[row])
            #  self.strata.stratumCompositions = np.hstack([self.strata.stratumCompositions[0:row], self.strata.stratumCompositions[row], self.strata.stratumCompositions[row:,]])
            #  self.strata.stratumThicknesses = np.hstack([self.strata.stratumThicknesses[0:row], self.strata.stratumThicknesses[row], self.strata.stratumThicknesses[row:,]])
            #  self.strata.stratumDopings = np.hstack([self.strata.stratumDopings[0:row], self.strata.stratumDopings[row], self.strata.stratumDopings[row:,]])
            #  self.strata.stratumRIndexes = np.hstack([self.strata.stratumRIndexes[0:row], self.strata.stratumRIndexes[row], self.strata.stratumRIndexes[row:,]])

        self.strata.stratumMaterials.insert(row, self.strata.stratumMaterials[row])

        self.strata.stratumCompositions = np.insert(
                self.strata.stratumCompositions, row, 
                self.strata.stratumCompositions[row])
        self.strata.stratumThicknesses = np.insert(
                self.strata.stratumThicknesses, row, 
                self.strata.stratumThicknesses[row])
        self.strata.stratumDopings = np.insert(
                self.strata.stratumDopings, row, 
                self.strata.stratumDopings[row])
        self.strata.stratumRIndexes = np.insert(
                self.strata.stratumRIndexes, row, 
                self.strata.stratumRIndexes[row])

        self.stratumTable_refresh()
        self.stratumTable.selectRow(row)
        self.stratumTable.setFocus()

        self.dirty = True
        self.update_windowTitle()

    def insert_stratumBelow(self):
        row = self.stratumTable.currentRow()
        if row == -1:
            return

        #  if row == self.strata.stratumDopings.size-1:
            #  self.strata.stratumMaterials.insert(row, self.strata.stratumMaterials[row])
            #  self.strata.stratumCompositions = np.hstack([self.strata.stratumCompositions[:], self.strata.stratumCompositions[row]])
            #  self.strata.stratumThicknesses = np.hstack([self.strata.stratumThicknesses[:], self.strata.stratumThicknesses[row]])
            #  self.strata.stratumDopings = np.hstack([self.strata.stratumDopings[:], self.strata.stratumDopings[row]])
            #  self.strata.stratumRIndexes = np.hstack([self.strata.stratumRIndexes[:], self.strata.stratumRIndexes[row]])

        #  else:
            #  self.strata.stratumMaterials.insert(row, self.strata.stratumMaterials[row])
            #  self.strata.stratumCompositions = np.hstack([self.strata.stratumCompositions[0:row+1], self.strata.stratumCompositions[row], self.strata.stratumCompositions[row+1:,]])
            #  self.strata.stratumThicknesses = np.hstack([self.strata.stratumThicknesses[0:row+1], self.strata.stratumThicknesses[row], self.strata.stratumThicknesses[row+1:,]])
            #  self.strata.stratumDopings = np.hstack([self.strata.stratumDopings[0:row+1], self.strata.stratumDopings[row], self.strata.stratumDopings[row+1:,]])
            #  self.strata.stratumRIndexes = np.hstack([self.strata.stratumRIndexes[0:row+1], self.strata.stratumRIndexes[row], self.strata.stratumRIndexes[row+1:,]])

        self.strata.stratumMaterials.insert(row, self.strata.stratumMaterials[row])
        self.strata.stratumCompositions = np.insert(
                self.strata.stratumCompositions, row+1, 
                self.strata.stratumCompositions[row])
        self.strata.stratumThicknesses = np.insert(
                self.strata.stratumThicknesses, row+1, 
                self.strata.stratumThicknesses[row])
        self.strata.stratumDopings = np.insert(
                self.strata.stratumDopings, row+1, 
                self.strata.stratumDopings[row])
        self.strata.stratumRIndexes = np.insert(
                self.strata.stratumRIndexes, row+1, 
                self.strata.stratumRIndexes[row])

        self.stratumTable_refresh()
        self.stratumTable.selectRow(row+1)
        self.stratumTable.setFocus()

        self.dirty = True
        self.update_windowTitle()

    def delete_stratum(self):
        #don't delete last stratum
        if self.strata.stratumDopings.size == 1:
            return

        row = self.stratumTable.currentRow()
        if row == -1:
            return

        self.strata.stratumMaterials.pop(row)
        #  self.strata.stratumCompositions = np.hstack([self.strata.stratumCompositions[0:row], self.strata.stratumCompositions[row+1:,]])
        #  self.strata.stratumThicknesses = np.hstack([self.strata.stratumThicknesses[0:row], self.strata.stratumThicknesses[row+1:,]])
        #  self.strata.stratumDopings = np.hstack([self.strata.stratumDopings[0:row], self.strata.stratumDopings[row+1:,]])
        #  self.strata.stratumRIndexes = np.hstack([self.strata.stratumRIndexes[0:row], self.strata.stratumRIndexes[row+1:,]])
        self.strata.stratumCompositions = np.delete(
                self.strata.stratumCompositions, row)
        self.strata.stratumThicknesses = np.delete(
                self.strata.stratumThicknesses, row)
        self.strata.stratumDopings = np.delete(
                self.strata.stratumDopings, row)
        self.strata.stratumRIndexes = np.delete(
                self.strata.stratumRIndexes, row)

        #if current row was last row (now deleted)
        if row+1 > self.strata.stratumThicknesses.size:
            self.strata.stratumSelected -= 1
            row -= 1

        self.stratumTable.selectRow(row)
        self.stratumTable_refresh()
        self.stratumTable.selectRow(row)
        self.stratumTable.setFocus()

        self.dirty = True
        self.update_windowTitle()        




#===============================================================================
# Optical Tab Plotting and Plot Control
#===============================================================================

    def update_opticalCanvas(self):
        self.strata.populate_x()

        self.opticalCanvas.clear()

        self.curvenR = Qwt.QwtPlotCurve()
        self.curvenR.setData(self.strata.xPoints,self.strata.xn.real)
        self.curvenR.setPen(QPen(Qt.black, 1.5))
        if settings.antialiased:
            self.curvenR.setRenderHint(Qwt.QwtPlotItem.RenderAntialiased)
        self.curvenR.attach(self.opticalCanvas)
        self.opticalCanvas.setAxisTitle(Qwt.QwtPlot.yLeft, 'Refractive Index')

        if self.strata.stratumSelected >= 0 and \
                self.strata.stratumSelected < self.strata.stratumThicknesses.size:
            mask = ~np.isnan(self.strata.xStratumSelected)
            self.stratumSelection = SupportClasses.MaskedCurve(
                    self.strata.xPoints, self.strata.xStratumSelected,
                    mask)
            self.stratumSelection.setPen(QPen(Qt.blue, 2))
            if settings.antialiased:
                self.stratumSelection.setRenderHint(Qwt.QwtPlotItem.RenderAntialiased)
            self.stratumSelection.attach(self.opticalCanvas)

        #plot Intensity
        if hasattr(self.strata,'xI'):
            self.curvexI = Qwt.QwtPlotCurve()
            self.curvexI.setData(self.strata.xPoints, 
                    self.strata.xI*self.strata.stratumRIndexes.real.max())
            self.curvexI.setPen(QPen(Qt.red, 1.5))
            if settings.antialiased:
                self.curvexI.setRenderHint(Qwt.QwtPlotItem.RenderAntialiased)
            self.curvexI.attach(self.opticalCanvas)
            self.opticalCanvas.setAxisTitle(Qwt.QwtPlot.yLeft, 
                    'Refractive Index, Mode Intensity')

        self.opticalCanvas.setAxisTitle(Qwt.QwtPlot.xBottom, u'Position (\u03BCm)')
        self.opticalCanvas.replot()

    def solve_mode(self):
        betaInit = self.betaEffBox.text()
        if betaInit == '':
            betaInit = None
        else:
            betaInit = complex(str(betaInit))
        self.strata.beta = self.strata.beta_find(betaInit)
        self.betaEffBox.setText(
                '{0.real:2.3f}+{0.imag:1.3e}j'.format(self.strata.beta))
        self.strata.mode_plot()
        self.update_modeCalculations_box()
        self.update_opticalCanvas()

    def run_opti1D(self):
        #get initial parameters
        try:
            optiType1D  = self.opti1DChoiceBox.currentText()
            strata1D    = np.array(SupportClasses.matlab_range(
                self.opti1DLayerBox.text()), dtype=int)
            strata1D   -= 1 #indexing starts at 0
            optiRange1D = np.array(SupportClasses.matlab_range(
                self.opti1DRangeBox.text()))
        except ValueError:
            QMessageBox.warning(self,"ErwinJr Error", "Invalid entry.")
            return

        #set up GUI
        self.optiFrame.setEnabled(False)
        self.plotModeButton.setEnabled(False)

        Jth0Array = np.zeros(optiRange1D.size)*np.NaN
        ylabel = '<i>J<sub>th0</sub></i>'

        stratumThicknessesInitial = self.strata.stratumThicknesses.copy()
        stratumDopingsInitial = self.strata.stratumDopings.copy()

        for q, rangeValue in enumerate(optiRange1D):
            if optiType1D == 'Thickness':
                self.strata.stratumThicknesses[strata1D] = rangeValue
                xlabel = u'Thickness (\u03BCm)'
            elif optiType1D == 'Doping':
                self.strata.stratumDopings[strata1D] = rangeValue
                xlabel = 'Doping (x10<sup>17</sup> cm<sup>-3</sup>)'
            elif optiType1D == 'Active Core Periods':
                pass
            elif optiTyp1D == 'deltaE':
                pass
            elif optiTyp1D == 'E3c':
                pass
            elif optiTyp1D == 'Custom Facet':
                pass
            elif optiTyp1D == 'Waveguide Length':
                pass
            elif optiTyp1D == 'Ridge Width':
                pass
            elif optiType1D == 'Tsink':
                pass
            self.stratumTable_refresh()
            self.update_stratum_inputBoxes()
            self.solve_mode()
            Jth0Array[q] = self.strata.Jth0
            self.plot_on_optimization1DCanvas(optiRange1D, xlabel, Jth0Array, ylabel)

        #reset initial values
        self.strata.stratumThicknesses = stratumThicknessesInitial
        self.strata.stratumDopings = stratumDopingsInitial
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()
        self.solve_mode()

        #reset GUI
        self.optiFrame.setEnabled(True)
        self.plotModeButton.setEnabled(True)

    def plot_on_optimization1DCanvas(self, xVals, xlabel, yVals, ylabel):
        self.optimization2DCanvas.setVisible(False)
        self.optimization1DCanvas.setVisible(True)

        self.optimization1DCanvas.clear()

        mask = ~np.isnan(yVals)
        optiCurve = SupportClasses.MaskedCurve(xVals, yVals, mask)
        optiCurve.setPen(QPen(Qt.blue, 1.5))
        if settings.antialiased:
            optiCurve.setRenderHint(Qwt.QwtPlotItem.RenderAntialiased)
        optiCurve.attach(self.optimization1DCanvas)

        self.optimization1DCanvas.setAxisTitle(Qwt.QwtPlot.xBottom, xlabel)
        self.optimization1DCanvas.setAxisTitle(Qwt.QwtPlot.yLeft, ylabel)
        self.optimization1DCanvas.replot()

    def run_opti2D(self):
        #get initial parameters
        try:
            optiType1D  = self.opti1DChoiceBox.currentText()
            strata1D    = np.array(SupportClasses.matlab_range(
                self.opti1DLayerBox.text()), dtype=int)
            strata1D   -= 1 #indexing starts at 0
            optiRange1D = np.array(SupportClasses.matlab_range(
                self.opti1DRangeBox.text()))
            optiType2D  = self.opti2DChoiceBox.currentText()
            strata2D    = np.array(SupportClasses.matlab_range(
                self.opti2DLayerBox.text()), dtype=int)
            strata2D   -= 1 #indexing starts at 0
            optiRange2D = np.array(SupportClasses.matlab_range(
                self.opti2DRangeBox.text()))
        except ValueError:
            QMessageBox.warning(self,"ErwinJr Error", "Invalid entry.")
            return

        Jth0Array = np.NaN * np.zeros((optiRange1D.size, optiRange2D.size))
        zlabel = '$J_{th0}$'

        stratumThicknessesInitial = self.strata.stratumThicknesses.copy()
        stratumDopingsInitial = self.strata.stratumDopings.copy()

        for qq, rangeValue2D in enumerate(optiRange2D):
            if optiType2D == 'Thickness':
                self.strata.stratumThicknesses[strata2D] = rangeValue2D
                xlabel = u'Thickness ($\mu m$)'
            elif optiType2D == 'Doping':
                self.strata.stratumDopings[strata2D] = rangeValue2D
                xlabel = 'Doping ($x10^{17} cm^{-3}$)'
            for q, rangeValue1D in enumerate(optiRange1D):
                if optiType1D == 'Thickness':
                    self.strata.stratumThicknesses[strata1D] = rangeValue1D
                    ylabel = u'Thickness ($\mu m$)'
                elif optiType1D == 'Doping':
                    self.strata.stratumDopings[strata1D] = rangeValue1D
                    ylabel = 'Doping ($x10^{17} cm^{-3}$)'

                self.stratumTable_refresh()
                self.update_stratum_inputBoxes()
                self.solve_mode()
                Jth0Array[q,qq] = self.strata.Jth0
                self.plot_on_optimization2DCanvas(optiRange1D, xlabel, 
                        optiRange2D, ylabel, Jth0Array, zlabel)
                QCoreApplication.processEvents()

        #reset initial values
        self.strata.stratumThicknesses = stratumThicknessesInitial
        self.strata.stratumDopings = stratumDopingsInitial
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()
        self.solve_mode()

        #reset GUI
        self.optiFrame.setEnabled(True)
        self.plotModeButton.setEnabled(True)

    def plot_on_optimization2DCanvas(self, 
            xVals, xlabel, yVals, ylabel, zVals, zlabel):
        self.optimization1DCanvas.setVisible(False)
        self.optimization2DCanvas.setVisible(True)

        X,Y = meshgrid(yVals, xVals)
        Z = zVals

        self.optimization2DAxes.cla()
        self.optimization2DAxes.patch.set_color(self.bgColor)
        self.optimization2DAxes.mouse_init()

        normd = matplotlib.colors.Normalize(
                np.nanmin(np.nanmin(Z)), np.nanmax(np.nanmax(Z)))
        self.optimization2DAxes.plot_surface(X, Y, Z, 
                cstride=1, rstride=1, norm=normd, cmap=matplotlib.cm.Blues_r, 
                linewidth=0, antialiased=False, shade=False)
        self.optimization2DAxes.set_zlim(0.95*np.nanmin(np.nanmin(Z)), 
                1.05*np.nanmax(np.nanmax(Z)))
        self.optimization2DAxes.set_xlabel(xlabel)
        self.optimization2DAxes.set_ylabel(ylabel)
        self.optimization2DAxes.set_zlabel(zlabel)
        self.optimization2DCanvas.draw()


    def quantumUpdate(self):
        """ SLOT connected to self.quantumWidget.dirty """
        self.dirty = True
        self.update_windowTitle()


    def transfer_optical_parameters(self):
        """ SLOT connected to SIGNAL self.quantumWidget.toOpticalParamsButton.clicked()
        Transfer informatin in quantum tab to optical tab"""
        self.quantumWidget.transfer_params(self.strata)

        #GUI settings
        self.quantumWidget.toOpticalParamsButton.setEnabled(False)
        self.editOpticalParametersBox.setChecked(False)

        #update all the input boxes
        self.update_stratum_inputBoxes()

        #update the stratumTable
        self.stratumTable_refresh()



#===============================================================================
# General Menu Functions
#===============================================================================

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

    def create_action(self, text, slot=None, shortcut=None, 
                        icon=None, tip=None, checkable=False, 
                        signal="triggered()", ischecked=False):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon("images/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        if ischecked:
            action.setChecked(True)
        return action

    def create_Quantum_menu(self):
        #file menu
        self.file_menu = self.menuBar().addMenu("&File")
        newFileAction      = self.create_action("&New...", 
                slot=self.fileNew, shortcut=QKeySequence.New, 
                icon="filenew", tip="New ErwinJr file")
        openFileAction     = self.create_action("&Open", 
                shortcut="Ctrl+O", slot=self.fileOpen, 
                tip="Open ErwinJr file", icon="fileopen")
        saveFileAction     = self.create_action("&Save", 
                shortcut="Ctrl+S", slot=self.fileSave, 
                tip="Save ErwinJr file", icon="filesave")
        saveAsFileAction   = self.create_action("S&ave As", 
                shortcut="Ctrl+W", slot=self.fileSaveAs, 
                tip="Save ErwinJr file as", icon="filesaveas")
        exportQuantumCanvasAction = self.create_action("Export Band Diagram Image", 
                slot=self.exportBandDiagram, 
                tip="Export Band Diagram Image")
        exportBandCSVAction = self.create_action("Export Band Diagram Data", 
                slot=self.export_band_diagram_data, 
                tip="Export Band Diagram Data")
        quit_action = self.create_action("&Quit", 
                slot=self.close, shortcut="Ctrl+Q", 
                tip="Close the application", icon="filequit")
        self.fileMenuActions = (newFileAction, openFileAction, 
                saveFileAction, saveAsFileAction, None, 
                exportBandCSVAction, exportQuantumCanvasAction, None, 
                quit_action)
        self.connect(self.file_menu, SIGNAL("aboutToShow()"), 
                self.updateFileMenu)
        #  self.add_actions(self.file_menu, 
                #  (newFileAction, openFileAction, 
                    #  saveFileAction, saveAsFileAction, None, quit_action))

        #edit menu
        self.edit_menu = self.menuBar().addMenu("&Edit")
        temperatureAction = self.create_action("&Temperature", 
                slot=self.set_temperature, tip="Set temperature")
        bumpLayerAction = self.create_action("&Bump First Layer", 
                slot=self.quantumWidget.bump_first_layer, 
                tip="Move zeroth layer to first layer")
        copyStructureAction = self.create_action("&Copy Structure", 
                slot=self.quantumWidget.copy_structure, 
                tip="Copy Layer Structure to Clipboard")
        self.add_actions(self.edit_menu, 
                (temperatureAction, bumpLayerAction, None, copyStructureAction))

        #view menu
        self.view_menu = self.menuBar().addMenu("&View")
        VXBandAction = self.create_action("X Valley Conduction Band", 
                checkable=True, ischecked=self.quantumWidget.plotVX,
                slot=self.quantumWidget.view_VXBand)
        VLBandAction = self.create_action("L Valley Conduction Band", 
                checkable=True, ischecked=self.quantumWidget.plotVL,
                slot=self.quantumWidget.view_VLBand)
        LHBandAction = self.create_action("Light Hole Valence Band", 
                checkable=True, ischecked=self.quantumWidget.plotLH,
                slot=self.quantumWidget.view_LHBand)
        SOBandAction = self.create_action("Split Off Valence Band", 
                checkable=True, ischecked=self.quantumWidget.plotSO,
                slot=self.quantumWidget.view_SOBand)
        self.add_actions(self.view_menu, 
                (VXBandAction,VLBandAction,LHBandAction,SOBandAction))        

        #help menu
        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action("&About",shortcut='F1', 
                slot=self.on_about)
        licenses_action = self.create_action("&License", 
                slot=self.on_licenses)
        tutorialAction = self.create_action("&Tutorial", 
                slot=self.on_tutorial)
        self.add_actions(self.help_menu, 
                (tutorialAction,about_action,licenses_action))

    def create_Optical_menu(self):
        #file menu
        self.file_menu = self.menuBar().addMenu("&File")
        newFileAction      = self.create_action("&New...", 
                self.fileNew, QKeySequence.New, 
                "filenew", "New ErwinJr file")
        openFileAction     = self.create_action("&Open", 
                shortcut="Ctrl+O", slot=self.fileOpen, 
                tip="Open ErwinJr file", icon="fileopen")
        saveFileAction     = self.create_action("&Save", 
                shortcut="Ctrl+S", slot=self.fileSave, 
                tip="Save ErwinJr file", icon="filesave")
        saveAsFileAction   = self.create_action("S&ave As", 
                shortcut="Ctrl+W", slot=self.fileSaveAs, 
                tip="Save ErwinJr file as", icon="filesaveas")
        quit_action = self.create_action("&Quit", 
                slot=self.close, shortcut="Ctrl+Q", 
                tip="Close the application", icon="filequit")
        self.fileMenuActions = (newFileAction, openFileAction, 
                saveFileAction, saveAsFileAction, None, quit_action)
        self.connect(self.file_menu, SIGNAL("aboutToShow()"), 
                self.updateFileMenu)
        #self.add_actions(self.file_menu, (newFileAction, openFileAction, saveFileAction, saveAsFileAction, None, quit_action))

        #help menu
        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action("&About",shortcut='F1', slot=self.on_about)
        licenses_action = self.create_action("&License", slot=self.on_licenses)
        tutorialAction = self.create_action("&Tutorial", slot=self.on_tutorial)
        self.add_actions(self.help_menu, (tutorialAction,about_action,licenses_action))




#===============================================================================
# File Menu Items
#===============================================================================

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
                self.connect(action, SIGNAL("triggered()"),
                             partial(self.fileOpen,fname))
                self.file_menu.addAction(action)
        self.file_menu.addSeparator()
        self.file_menu.addAction(self.fileMenuActions[-1])

    def addRecentFile(self, fname):
        if fname is None:
            return
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
            self.strata.populate_rIndexes()
            self.update_stratum_inputBoxes()
            self.stratumTable_refresh()
            self.opticalCanvas.clear()
            self.update_opticalCanvas()

        self.filename = fname
        self.addRecentFile(fname)
        self.dirty = False
        self.update_windowTitle()       

    def fileNew(self):
        if not self.okToContinue():
            return False

        self.filename = None
        self.quantumWidget.plotDirty = False
        self.quantumWidget.quantumCanvas.clear()
        self.opticalCanvas.clear()
        self.optimization1DCanvas.clear()

        self.quantumWidget.qclayers = QCLayers()
        self.quantumWidget.reloaded()

        self.strata = Strata()

        self.update_stratum_inputBoxes()

        self.update_opticalCanvas()
        self.stratumTable_refresh()

        self.dirty = False
        self.update_windowTitle()

        return True

    def okToContinue(self):
        if self.dirty:
            reply = QMessageBox.question(self, 
                    "ErwinJr " + str(majorVersion) + " - Unsaved Changes", 
                    "Save unsaved changes?", 
                    QMessageBox.Yes|QMessageBox.No|QMessageBox.Cancel)
            if reply == QMessageBox.Cancel:
                return False
            elif reply == QMessageBox.Yes:
                self.fileSave()
        return True

    def update_windowTitle(self):
        if self.filename is not None:
            self.setWindowTitle("ErwinJr " + str(majorVersion) + 
                    " - %s[*]" % os.path.basename(str(self.filename)))
        else:
            self.setWindowTitle("ErwinJr " + str(majorVersion) + "[*]")
        self.setWindowModified(self.dirty)

    def fileOpen(self, fname = None):
        #if not self.okToContinue():
        #    return False

        #clear all old data, also calls self.okToContinue()
        if not self.fileNew(): 
            return False
        if fname is None:
            dir = os.path.dirname(str(self.filename)) if self.filename else "."
            fname =unicode(QFileDialog.getOpenFileName(self,
                "ErwinJr - Choose file", dir, 
                "ErwinJr files (*.qcl)\nAll files (*.*)"))
        #open file and determine if it is from the Matlab version of ErwinJr
        filehandle = open(fname, 'r')
        firstLine = filehandle.readline()
        filehandle.close()
        if fname:
            if firstLine.split(':')[0] == 'Description':
                QMessageBox.warning(self, 'ErwinJr Error', 
                        'Older .qcl format is no longer supported '
                        'for Ver>3.0.')
                #  self.qclPtonLoad(fname)
            elif firstLine == 'ErwinJr Data File' + self.newLineChar:
                self.qclLoad(fname)
            else:
                QMessageBox.warning(self,'ErwinJr Error',
                        'Could not recognize input file.')
                return
            self.quantumWidget.reloaded()

            #if firstLine == 'ErwinJr Data File\n':
            self.strata.populate_rIndexes()
            self.update_stratum_inputBoxes()
            self.stratumTable_refresh()
            self.opticalCanvas.clear()
            self.update_opticalCanvas()

        self.filename = fname
        self.addRecentFile(fname)
        self.dirty = False
        self.update_windowTitle()

        return True


    def qclLoad(self, fname):
        #  print "Loading "+fname
        try:
            with open(fname, 'r') as f:
                SaveLoad.qclLoad(f, self.quantumWidget.qclayers, self.strata)
        except Exception as err:
            QMessageBox.warning(self,"ErwinJr - Warning",
                             "Could not load *.qcl file.\n"+str(err))

    def fileSave(self):
        if self.filename is None:
            return self.fileSaveAs()
        else:
            #os.path.extsep
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
        fname = unicode(QFileDialog.getSaveFileName(self,
            "ErwinJr - Save File", QString(fname), typeString))
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
                SaveLoad.qclSave(f, self.quantumWidget.qclayers, self.strata)
        except Exception as err:
            QMessageBox.warning(self,"ErwinJr - Warning",
                             "Could not save *.qcl file.\n"+str(err))
        return True

    def closeEvent(self, event):
        if self.okToContinue():
            qsettings = QSettings()
            filename = QVariant(QString(self.filename)) if self.filename \
                    else QVariant()
            qsettings.setValue("LastFile", filename)
            recentFiles = QVariant(self.recentFiles) if self.recentFiles \
                    else QVariant()
            qsettings.setValue("RecentFiles", recentFiles)
            qsettings.setValue("MainWindow/Geometry", QVariant(self.saveGeometry()))
            qsettings.setValue("MainWindow/State", QVariant(self.saveState()))
        else:
            event.ignore()




#===============================================================================
# Export Functions
#===============================================================================

    def exportBandDiagram(self):
        fname = unicode(QFileDialog.getSaveFileName(self,
            "ErwinJr - Export Band Structure Image",
            self.filename.split('.')[0], 
            "Portable Network Graphics file (*.png)"))
        if not fname:
            return

        #set background color to white and save presets
        bgRole = self.mainTabWidget.backgroundRole()
        self.mainTabWidget.setBackgroundRole(QPalette.Base)
        self.quantumWidget.export_quantumCanvas(fname)
        self.mainTabWidget.setBackgroundRole(bgRole)

    def export_band_diagram_data(self):
        fname = unicode(QFileDialog.getSaveFileName(self,
            "ErwinJr - Export Band Structure Data",
            self.filename.split('.')[0], 
            "Comma-Separated Value file (*.csv)"))
        if fname != '': #if user doesn't click cancel
            self.quantumWidget.export_band_data(fname)


#===============================================================================
# Edit Menu Items
#===============================================================================
    def set_temperature(self):
        # TODO: make it explict 
        nowTemp = cst.Temperature
        newTemp, buttonResponse = QInputDialog.getDouble(self, 
                'ErwinJr Input Dialog', 'Set Temperature', 
                value=nowTemp, min=0)
        if buttonResponse:
            cst.set_temperature(newTemp)
            self.quantumWidget.set_temperature(newTemp)




#===============================================================================
# Help Menu Items
#===============================================================================

    def on_about(self):
        msg = """ ErwinJr 2.x Authors and Contributors

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
        QMessageBox.about(self, "ErwinJr " + str(ejVersion), copyright1.strip())

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
    qsettingsSystem = QSettings(QSettings.SystemScope,"JPL","ErwinJr")
    installDirectory = str(qsettingsSystem.value('installDirectory').toString())
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

    #this block handles a filename passed in by command line
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

    # Import Psyco if available
    #  try:
        #  import psyco
        #  psyco.full()
    #  except ImportError:
        #  noPsycoBox = QMessageBox(QMessageBox.Question, 'EwrinJr '+str(majorVersion), "Psyco could not be loaded.\nExecution will be slowed.")
        #  noPsycoBox.exec_()

    qsettings = QSettings()
    if not qsettings.value('firstRun').toInt()[1]:
        if not installDirectory:
            qsettingsSystem.setValue("installDirectory", QVariant(os.getcwd()))
        firstRunBox = QMessageBox(QMessageBox.Question, 
                'EwrinJr '+str(majorVersion), 
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
