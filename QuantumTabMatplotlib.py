#!/usr/bin/env python2
# -*- coding:utf-8 -*-

#============================================================================
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
#============================================================================

# TODO: 
# In plot controls, add "show one period" 
# save and load pickle for qclayers
# replace qwt by matplotlib

from PyQt4.QtCore import *
from PyQt4.QtGui import *
#  from matplotlib.backend_bases import key_press_handler
#  from matplotlib.backends.backend_qt4agg \
        #  import NavigationToolbar2QT as NavigationToolbar
import sys
import numpy as np
from numpy import pi, sqrt
from functools import partial, wraps

import settings
from QCLayers import QCLayers, cst
from QCLayers import h, c0, e0
from EJcanvas import EJcanvas, EJplotControl

#============================================================================
# Debug options
#============================================================================
DEBUG = 2
if DEBUG >=3: 
    import pickle
    import inspect

class QuantumTab(QWidget):
    """The Quantum Tab of ErwinJr. This is designed to be a GUI wrapper of
    the class QCLayers
    Member viable (! label public interface): 
        !qclayers: A class to describe the physics of QC layers
        numMaterials: copy from qclayers, the number of difference kinds of 
                    materials supported
        materialList: for Material column in the layerTable
        colors: colors for different wavefunctions

        customized SIGNAL: 'dirty'

        --- state selection ---
        stateHolder: the list containing the index of seleted states. the
                    index is defined in qclayers
        pairSelected: Boolean flag for if a pair of states is selected
        --- plot flags ---
        !plotDirty: the quantumCanvas is to be updated (To check)
        solveType: 'basis' or 'whole', decide different kinds of solving
        !plotVX, !plotVL, !plotLH, !plotSO: show X point of conduction band (L
                    point, LH band, SO band)

        --- GUI widget ---
        1st column (settingBox):
            inputSubstrateBox
            inputEFieldBox
            inputHorzResBox
            inputVertResBox
            inputRepeatsBox
            inputARInjectorCheck
            inputInjectorARCheck
            LpFirstSpinbox          LpLastSpinbox
            targetWL_box
            goalFuncBox
            GlobalOptButton

        2nd column (layerBox):
            insertLayerAboveButton  deleteLayerButton
            OptimizeFoMButton       OptimizeDipoleButton
            layerTable

        3rd column (solveBox): 
            solveBasisButton
            solveWholeButton
            DescriptionBox
            mtrl_well               mtrl_barr
            MoleFracWellBox[n]      MoleFracBarrBox[n]     offsetLabel[n]
            strainDescription
            LOPhononDescription
            wellSelectButton
            zoominButton            zoomOutButton
            panButton               clearWFsButton
            pairSelectButton
            FoMButton               !toOpticalParamsButton
            pairSelectString

        4th column (figureBox):
            !quantumCanvas**

        --- *quantumCanvas ---
        zoomer
        picker
        picker2
        panner

        --- Optimization ---
        OptGoalsName: Name displayed for global optimization target functions
        OptGoalsFunc: the target functions for optimization
        OptGoal: the target function to optimize

    Member method (! label public interface): 
        !reloaded
        !get_nCore
        !layerNum
        !transfer_params
        !bump_first_layer
        !copy_structure
        !view_VXBand(VLBand, LHBand, SOBand)
        !export_quantumCanvas
        !export_band_data
        !set_temperature
        (private is omitted)
        """
    def __init__(self, parent=None):
        super(QuantumTab, self).__init__(parent)
        self.qclayers = QCLayers()
        self.numMaterials = self.qclayers.numMaterials #=8, to improve (TODO)

        # colors for different wavefunctions
        self.colors = ((0.584, 0.450, 0.701), (0.431, 0.486, 0.745),
                (0.576, 0.694, 0.517), (0.682, 0.780, 0.321),
                (0.501, 0.501, 0.509), (0.854, 0.741, 0.247), 
                (0.874, 0.607, 0.290), (0.823, 0.341, 0.278),
                (0.725, 0.321, 0.623), (0.411, 0.741, 0.270), 
                (0.078, 0.078, 0.078), (0.431, 0.803, 0.870), 
                (0.223, 0.321, 0.643))

        # Global Optimization target functions (goal)
        self.OptGoalsName = ('FoM', 'Dipole')
        self.OptGoalsFunc = (self.qclayers.figure_of_merit, 
                self.qclayers.dipole)

        self.updating = False
        self.plotDirty = False
        self.solveType = None
        self.plotVX = False
        self.plotVL = False
        self.plotLH = False
        self.plotSO = False

        self.stateHolder = []
        self.pairSelected = False

        # Platform dependent settings, eg. layerout size settings
        if sys.platform == 'win32': 
            layerTableSize = 340
            DescriptionBoxWidth = 190
            LpStringBoxWidth=135
        elif sys.platform == 'darwin':
            layerTableSize = 405
            DescriptionBoxWidth = 285
            LpStringBoxWidth=130
        elif sys.platform == 'linux2':
            layerTableSize = 365
            DescriptionBoxWidth = 240
            LpStringBoxWidth=150
        else:
            QMessageBox.warning(self, 'ErwinJr - Warning', 
                    'Platform %s not tested.'%sys.platform)
            layerTableSize = 340
            DescriptionBoxWidth = 190
            LpStringBoxWidth=135
        pairSelectStringWidth = DescriptionBoxWidth


        self.updating = True
        ########################################################
        # settingBox, containing all setting parameter
        ########################################################
        settingBox = QVBoxLayout()

        settingBox.addWidget(QLabel(
            "<center><b>Substrate</b></center>"))
        self.inputSubstrateBox = QComboBox()
        self.inputSubstrateBox.addItems(self.qclayers.substratesList)
        self.connect(self.inputSubstrateBox,
                SIGNAL("currentIndexChanged(const QString)"), 
                self.input_substrate)
        settingBox.addWidget(self.inputSubstrateBox)

        settingBox.addWidget(QLabel(
            '<center><b><i>E<sub>field</sub></i></b></center>'))
        self.inputEFieldBox = QDoubleSpinBox()
        self.inputEFieldBox.setDecimals(1)
        self.inputEFieldBox.setSuffix(' kV/cm')
        self.inputEFieldBox.setRange(0.0,250.0)
        self.connect(self.inputEFieldBox, 
                SIGNAL("valueChanged(double)"), 
                self.input_EField)
        settingBox.addWidget(self.inputEFieldBox)

        settingBox.addWidget(QLabel(
            '<center><b>Horizontal<br>Resolution</b></center>'))
        self.inputHorzResBox = QComboBox();
        # TODO: add unit here (angstrom)
        self.inputHorzResBox.addItems(['1.0','0.5','0.25','0.2','0.1'])
        # remove 0.25? because to 0.1
        self.connect(self.inputHorzResBox, 
                SIGNAL("currentIndexChanged(int)"), 
                self.input_horzRes)
        settingBox.addWidget(self.inputHorzResBox)

        settingBox.addWidget(QLabel(
            '<center><b>Vertical<br>Resolution</b></center>'))
        self.inputVertResBox = QDoubleSpinBox()
        self.inputVertResBox.setDecimals(2)
        self.inputVertResBox.setValue(0.5)
        self.inputVertResBox.setRange(0.0,10.0)
        self.inputVertResBox.setSingleStep(0.1)
        self.inputVertResBox.setSuffix(' meV')
        self.connect(self.inputVertResBox,
                SIGNAL("valueChanged(double)"), 
                self.input_vertRes)
        # TODO: check the SLOT if it can make use of the input
        settingBox.addWidget(self.inputVertResBox)

        settingBox.addWidget(QLabel(
            '<center><b>Structure Repeats</b></center>'))
        self.inputRepeatsBox = QSpinBox()
        self.inputRepeatsBox.setValue(1)
        self.inputRepeatsBox.setRange(1,5)
        self.connect(self.inputRepeatsBox, 
                SIGNAL("valueChanged(int)"), 
                self.input_repeats) 
        # TODO: check the SLOT if it can make use of the input
        settingBox.addWidget(self.inputRepeatsBox)

        # Basis solver devider setting
        basis_groupBox = QGroupBox("Basis Divisions")
        self.inputARInjectorCheck = QCheckBox("AR->Injector")
        self.inputInjectorARCheck = QCheckBox("Injector->AR")
        self.inputARInjectorCheck.setChecked(True)
        self.inputInjectorARCheck.setChecked(True)
        basisLayout = QVBoxLayout()
        basisLayout.addWidget(self.inputARInjectorCheck)
        basisLayout.addWidget(self.inputInjectorARCheck)
        basis_groupBox.setLayout(basisLayout)
        self.connect(self.inputARInjectorCheck,
                SIGNAL("stateChanged(int)"), 
                self.input_basis)
        self.connect(self.inputInjectorARCheck,
                SIGNAL("stateChanged(int)"), 
                self.input_basis)
        settingBox.addWidget(basis_groupBox)

        # Period information groupbox
        LpLayout_groupBox = QGroupBox("Period Info")
        self.LpFirstSpinbox = QSpinBox()
        self.LpFirstSpinbox.setValue(1)
        self.LpFirstSpinbox.setRange(1,1)
        self.connect(self.LpFirstSpinbox, 
                SIGNAL("valueChanged(int)"), 
                self.update_inputBoxes) 
        # TODO: change this slot to prevent unnecessary data updates
        self.LpLastSpinbox  = QSpinBox()
        self.LpLastSpinbox.setValue(1)
        self.LpLastSpinbox.setRange(1,1)
        self.connect(self.LpLastSpinbox,
                SIGNAL("valueChanged(int)"), 
                self.update_inputBoxes)
        # TODO: change this slot to prevent unnecessary data updates
        self.LpStringBox = QTextEdit('')
        self.LpStringBox.setReadOnly(True)
        self.LpStringBox.setSizePolicy(
                QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed))
        self.LpStringBox.setMaximumHeight(95)
        self.LpStringBox.setMaximumWidth(LpStringBoxWidth)
        LpLayout = QGridLayout()
        LpLayout.addWidget(QLabel('<b>first</b>'), 0,0)
        LpLayout.addWidget(QLabel('<b>last</b>'), 0,1)
        LpLayout.addWidget(self.LpFirstSpinbox, 1,0)
        LpLayout.addWidget(self.LpLastSpinbox, 1,1)
        LpLayout.addWidget(self.LpStringBox, 2,0, 1,2)
        LpLayout_groupBox.setLayout(LpLayout)
        settingBox.addWidget(LpLayout_groupBox)

        # Global Optimization groupbox
        GlobalOptLayout_groupBox = QGroupBox("Global Optimization")
        self.targetWL_box = QLineEdit('')
        self.targetWL_box.setValidator(QDoubleValidator(0,100,1))
        self.targetWL_box.setSizePolicy(QSizePolicy(
            QSizePolicy.Minimum,QSizePolicy.Ignored))
        #  self.targetWL_box.setMaximumWidth(50)
        self.connect(self.targetWL_box,
                SIGNAL("editingFinished()"), 
                self.set_targetWL)
        self.goalFuncBox = QComboBox()
        self.goalFuncBox.addItems(self.OptGoalsName)
        self.OptGoal = self.OptGoalsFunc[self.goalFuncBox.currentIndex()]
        self.connect(self.goalFuncBox, 
                SIGNAL("currentIndexChanged(int)"), 
                self.set_goal)
        self.GlobalOptButton = QPushButton("Optimize")
        self.connect(self.GlobalOptButton,
                SIGNAL("clicked()"),
                self.GlobalOptimization)
        GlobalOptLayout = QGridLayout()
        GlobalOptLayout.addWidget(QLabel(u"<b>\u03BB</b>:"), 0, 0)
        GlobalOptLayout.addWidget(self.targetWL_box, 0, 1)
        GlobalOptLayout.addWidget(QLabel('um'), 0, 2)
        GlobalOptLayout.addWidget(QLabel("<b>Target function</b>"), 1,0,1,3)
        GlobalOptLayout.addWidget(self.goalFuncBox, 2, 0, 1, 3)
        GlobalOptLayout.addWidget(self.GlobalOptButton, 3, 0, 1, 3)
        GlobalOptLayout_groupBox.setLayout(GlobalOptLayout)
        settingBox.addWidget(GlobalOptLayout_groupBox)

        settingBox.addStretch()
        #vBox1: settingBox end


        ########################################################
        # layerBox, containing layer table and related buttons
        ########################################################
        layerBox = QGridLayout()
        self.insertLayerAboveButton = QPushButton("Insert Layer Above")
        self.connect(self.insertLayerAboveButton, 
                SIGNAL("clicked()"),
                self.insert_layerAbove)
        layerBox.addWidget(self.insertLayerAboveButton, 0,0)
        self.deleteLayerButton = QPushButton("Delete Layer")
        self.connect(self.deleteLayerButton, 
                SIGNAL("clicked()"),
                self.delete_layer)
        layerBox.addWidget(self.deleteLayerButton, 0,1)        

        self.OptimizeFoMButton = QPushButton("Optimize Width (FoM)")
        self.OptimizeFoMButton.setEnabled(False)
        self.connect(self.OptimizeFoMButton, 
                SIGNAL("clicked()"), 
                partial(self.OptimizeLayer, 
                    goal = self.qclayers.figure_of_merit))
        self.OptimizeDipoleButton = QPushButton("Optimize Width (Dipole)")
        self.connect(self.OptimizeDipoleButton, 
                SIGNAL("clicked()"),
                partial(self.OptimizeLayer,
                    goal = self.qclayers.dipole))
        self.OptimizeDipoleButton.setEnabled(False)
        layerBox.addWidget(self.OptimizeFoMButton, 1,0)
        layerBox.addWidget(self.OptimizeDipoleButton, 1,1)

        #set up layerTable
        self.layerTable = QTableWidget()
        self.layerTable.setSelectionBehavior(QTableWidget.SelectRows)
        self.layerTable.setSelectionMode(QTableWidget.SingleSelection)
        self.layerTable.setMaximumWidth(layerTableSize)
        self.layerTable.setMinimumWidth(layerTableSize)
        self.connect(self.layerTable,
                SIGNAL("itemChanged(QTableWidgetItem*)"),
                self.layerTable_itemChanged)
        self.connect(self.layerTable,
                SIGNAL("itemSelectionChanged()"),
                self.layerTable_itemSelectionChanged)
        layerBox.addWidget(self.layerTable, 2,0, 1,2)
        #vBox2: layerBox end


        ########################################################
        # figureBox, containing the band and wavefunc figure
        ########################################################
        #set up quantumCanvas for band structure plot
        self.quantumCanvas = EJcanvas(xlabel=u'Position (Å)', 
                ylabel='Energy (eV)', parent=self)
        #  self.quantumCanvas.mpl_connect('key_press_event',
                #  partial(key_press_handler, canvas=self.quantumCanvas,
                    #  toolbar=self.mpl_toolbar))
        #  self.mpl_toolbar = NavigationToolbar(self.quantumCanvas, self)
        self.plotControl = EJplotControl(self.quantumCanvas, self)
        #  self.quantumCanvas.test()
        #  self.plotControl.save_figure()
        figureBox = QVBoxLayout()
        #  figureBox.addWidget(self.self.mpl_toolbar)
        figureBox.addWidget(self.quantumCanvas)
        #vBox4: figureBox end


        ########################################################
        # solveBox, containing all eigensolver and calculation control
        ########################################################
        solveBox = QVBoxLayout()
        self.solveBasisButton = QPushButton("Solve Basis")
        self.connect(self.solveBasisButton,
                SIGNAL("clicked()"),
                self.solve_basis)
        solveBox.addWidget(self.solveBasisButton)
        self.solveWholeButton = QPushButton("Solve Whole")
        self.connect(self.solveWholeButton,
                SIGNAL("clicked()"),
                self.solve_whole)
        solveBox.addWidget(self.solveWholeButton)

        #set up description box
        self.DescriptionBox = QTextEdit('')
        self.DescriptionBox.setReadOnly(False)
        self.DescriptionBox.setSizePolicy(QSizePolicy(
            QSizePolicy.Fixed,QSizePolicy.Fixed))
        self.DescriptionBox.setMaximumHeight(40)
        self.DescriptionBox.setMaximumWidth(DescriptionBoxWidth)
        self.connect(self.DescriptionBox, 
                SIGNAL("textChanged()"), 
                self.input_description)
        DescLayout = QVBoxLayout()
        DescLayout.addWidget(self.DescriptionBox)
        DescLayout_groupBox = QGroupBox("Description")
        DescLayout_groupBox.setLayout(DescLayout)
        solveBox.addWidget(DescLayout_groupBox)

        #set up material composition inputs
        #  self.mtrl_well = QLabel(
                #  '<center><b>In<sub>x</sub>Ga<sub>1-x</sub>As</b></center>')
        #  self.mtrl_barr = QLabel(
                #  '<center><b>Al<sub>1-x</sub>In<sub>x</sub>As</b></center>')
        self.mtrl_well = QLabel()
        self.mtrl_barr = QLabel()
        self.input_substrate('InP')
        # TODO: set these according to substrate dict
        self.MoleFracWellBox = []
        self.MoleFracBarrBox = []
        self.offsetLabel = []
        for n in range(self.numMaterials//2):
            self.MoleFracWellBox.append(QDoubleSpinBox())
            self.MoleFracWellBox[n].setDecimals(3)
            self.MoleFracWellBox[n].setValue(0.53)
            self.MoleFracWellBox[n].setRange(0.0, 1.0)
            self.MoleFracWellBox[n].setSingleStep(0.001)
            self.connect(self.MoleFracWellBox[n],
                    SIGNAL("editingFinished()"), 
                    partial(self.input_moleFrac, 2*n))
            self.MoleFracBarrBox.append(QDoubleSpinBox())
            self.MoleFracBarrBox[n].setDecimals(3)
            self.MoleFracBarrBox[n].setValue(0.52)
            self.MoleFracBarrBox[n].setRange(0.0, 1.0)
            self.MoleFracBarrBox[n].setSingleStep(0.001)
            self.connect(self.MoleFracBarrBox[n],
                    SIGNAL("editingFinished()"), 
                    partial(self.input_moleFrac, 2*n+1))
            self.offsetLabel.append(QLabel(''))
        self.strainDescription = QLabel('')
        self.LOPhononDescription = QLabel('')
        mtrl_grid = QGridLayout()
        mtrl_grid.addWidget(QLabel(
            '<center><b>Mole Fractions</b></center>'), 0,0, 1,4)
        mtrl_grid.addWidget(self.mtrl_well, 1,1)
        mtrl_grid.addWidget(self.mtrl_barr, 1,2)
        mtrl_grid.addWidget(QLabel('<center><b>Offset</b></center>'), 1,3)
        for n in range(self.numMaterials//2):
            mtrl_grid.addWidget(QLabel(
                '<center><b>#%d</b></center>'%(n+1)), 2+n, 0)
            mtrl_grid.addWidget(self.MoleFracWellBox[n], 2+n, 1)
            mtrl_grid.addWidget(self.MoleFracBarrBox[n], 2+n, 2)
            mtrl_grid.addWidget(self.offsetLabel[n], 2+n, 3)
        mtrl_grid.addWidget(QLabel('<center>(well)</center>'), 6,1)
        mtrl_grid.addWidget(QLabel('<center>(barrier)</center>'), 6,2)
        mtrl_grid.addWidget(self.strainDescription, 7,0, 1,4)
        mtrl_grid.addWidget(self.LOPhononDescription, 8,0, 1,4)
        mtrl_groupBox = QGroupBox()
        mtrl_groupBox.setLayout(mtrl_grid)
        solveBox.addWidget(mtrl_groupBox)

        #set up plot control inputs
        self.zoominButton = QPushButton("Zoom In")
        self.plotControl.set_action('zoom', self.zoominButton)
        self.zoomOutButton = QPushButton("Zoom Out")
        self.plotControl.set_action('home', self.zoomOutButton)
        self.panButton = QPushButton("Pan") # to move
        self.plotControl.set_action('pan', self.panButton)
        self.wellSelectButton = QPushButton("Layer Select")
        self.plotControl.set_custom('wellselect', 
                self.wellSelectButton, 
                self.well_select)
        self.connect(self.wellSelectButton, SIGNAL("clicked()"), 
                partial(self.plotControl.custom, 'wellselect'))
        self.clearWFsButton = QPushButton("Clear")
        self.connect(self.clearWFsButton,
                SIGNAL("clicked()"),
                self.clear_WFs)
        plotControl_grid = QGridLayout()
        plotControl_grid.addWidget(self.wellSelectButton, 0,0, 1,2)
        plotControl_grid.addWidget(self.zoominButton, 1,0, 1,1)
        plotControl_grid.addWidget(self.zoomOutButton, 1,1, 1,1)
        plotControl_grid.addWidget(self.panButton, 2,0, 1,1)
        plotControl_grid.addWidget(self.clearWFsButton, 2,1, 1,1)
        plotControl_groupBox = QGroupBox("Plot Controls")
        plotControl_groupBox.setLayout(plotControl_grid)
        solveBox.addWidget(plotControl_groupBox)

        #set up Calculate controls
        self.pairSelectButton = QPushButton("Pair Select")
        self.plotControl.set_custom('pairselect', 
                self.pairSelectButton,
                self.state_pick)
        self.connect(self.pairSelectButton, SIGNAL("clicked()"), 
                partial(self.plotControl.custom, 'pairselect'))
        self.FoMButton = QPushButton("Figure of Merit")
        self.FoMButton.setEnabled(False)
        self.connect(self.FoMButton, 
                SIGNAL("clicked()"), 
                self.updateFoM)
        self.toOpticalParamsButton = QPushButton("-> Optical Params")
        self.toOpticalParamsButton.setEnabled(False)
        # SLOT connected to SIGNAL of tranferOpticalParametersButton should
        # be defined out side this tab, in the main window
        self.pairSelectString = QTextEdit('')
        self.pairSelectString.setReadOnly(True)
        self.pairSelectString.setMaximumWidth(pairSelectStringWidth)
        self.pairSelectString.setSizePolicy(QSizePolicy(
            QSizePolicy.Fixed,QSizePolicy.Fixed))
        calculateControl_grid = QGridLayout()
        calculateControl_grid.addWidget(self.pairSelectButton, 0,0, 1,2)
        calculateControl_grid.addWidget(self.FoMButton, 1,0, 1,1)
        calculateControl_grid.addWidget(self.toOpticalParamsButton, 1,1, 1,1)
        calculateControl_grid.addWidget(self.pairSelectString, 2,0, 1,2)
        calculateControl_groupBox = QGroupBox("Calculate")
        calculateControl_groupBox.setLayout(calculateControl_grid)
        solveBox.addWidget(calculateControl_groupBox)

        solveBox.addStretch()
        #vBox3: solveBox end


        quantumLayout = QHBoxLayout()
        quantumLayout.addLayout(settingBox)
        quantumLayout.addLayout(layerBox)
        quantumLayout.addLayout(solveBox)
        quantumLayout.addLayout(figureBox)

        self.setLayout(quantumLayout)
        self.setAutoFillBackground(True)
        self.setBackgroundRole(QPalette.Window)

        self.layerTable_refresh()
        self.layerTable.selectRow(1)
        self.updating = False

        self.update_inputBoxes()
        self.layerTable.setFocus()
    # __init__ end

    def settingslot(fn):
        @wraps(fn)
        def wrapper(self, *args, **kwargs):
            if self.updating:
                return
            else:
                #  print args
                #  print kwargs
                return fn(self, *args, **kwargs)
        return wrapper

    def reloaded(self):
        self.update_Lp_limits()
        self.update_inputBoxes()
        self.layerTable_refresh()
        self.updating = True
        self.layerTable.selectRow(1)
        self.layerTable.setFocus()
        self.updating = False
        self.qclayers.populate_x()
        self.plotDirty = True
        self.update_quantumCanvas()

#==========================================================================
# Input Controls
#==========================================================================
    def update_inputBoxes(self):
        """ Update all input boxes except for the layerTable. 
        SLOT connected to LpFirstSpinbox/LpLastSpinBox.valueChanged(int)
        TODO: add GlobalOpt boxes
        """
        self.updating = True
        try:
            self.inputSubstrateBox.setCurrentIndex(
                    self.qclayers.substratesList.index(
                        self.qclayers.substrate))
        except Exception as err:
            QMessageBox.warning(self,"ErwinJr - Warning",
                             "Substrate data wrong.\n"+str(err))

        self.qclayers.update_alloys()
        self.qclayers.update_strain()
        self.qclayers.populate_x()
        for n in range(self.numMaterials//2):
            self.MoleFracWellBox[n].setValue(self.qclayers.moleFrac[2*n])
            self.MoleFracBarrBox[n].setValue(self.qclayers.moleFrac[2*n+1])
            self.offsetLabel[n].setText("%6.0f meV" %
                    ((self.qclayers.EcG[2*n+1]-self.qclayers.EcG[2*n])*1000))

        self.DescriptionBox.setText(self.qclayers.description)
        strainString = ("<center>Net Strain: <b>%6.3f%%</b></center>" %
                self.qclayers.netStrain)
        self.strainDescription.setText(strainString)
        hwLOString = ("<center>LO phonon: <b>%4.1f ~ %4.1f meV</b></center>" %
                (min(self.qclayers.hwLO)*1000, max(self.qclayers.hwLO)*1000))
        self.LOPhononDescription.setText(hwLOString)

        self.inputVertResBox.setValue(self.qclayers.vertRes)
        self.inputEFieldBox.setValue(self.qclayers.EField)
        self.inputHorzResBox.setCurrentIndex(self.inputHorzResBox.findText( 
            QString(unicode(self.qclayers.xres))))
        self.inputRepeatsBox.setValue(self.qclayers.repeats)
        self.update_Lp_box()

        self.emit(SIGNAL('dirty'))
        self.updating = False

    def input_substrate(self, substrateType):
        """ SLOT connected to inputSubstrateBox.currentIndexChanged(QString)
        update substrate chosen """
        # TODO put the following infor mation to a dict
        if substrateType == 'InP':
            self.qclayers.substrate = 'InP'
            self.materialList = ['InGaAs/AlInAs #1', 
                    'InGaAs/AlInAs #2', 
                    'InGaAs/AlInAs #3', 
                    'InGaAs/AlInAs #4']
            self.mtrl_well.setText( '<center><b>\
                    In<sub>x</sub>Ga<sub>1-x</sub>As\
                    </b></center>')
            self.mtrl_barr.setText('<center><b>\
                    Al<sub>1-x</sub>In<sub>x</sub>As\
                    </b></center')

        elif substrateType == 'GaAs':
            self.qclayers.substrate = 'GaAs'
            self.materialList = ['AlGaAs/AlGaAs #1', 
                    'AlGaAs/AlGaAs #2', 
                    'AlGaAs/AlGaAs #3', 
                    'AlGaAs/AlGaAs #4']
            self.mtrl_well.setText('<center><b>\
                    Al<sub>x</sub>Ga<sub>1-x</sub>As\
                    </b></center')
            self.mtrl_barr.setText('<center><b>\
                    Al<sub>x</sub>Ga<sub>1-x</sub>As\
                    </b></center')

        elif substrateType == 'GaSb':
            self.qclayers.substrate = 'GaSb'
            self.materialList = ['InAsSb/AlGaSb #1', 
                    'InAsSb/AlGaSb #2', 
                    'InAsSb/AlGaSb #3', 
                    'InAsSb/AlGaSb #4']
            self.mtrl_well.setText('<center><b>\
                    InAs<sub>y</sub>Sb<sub>1-y</sub>\
                    </b></center')
            self.mtrl_barr.setText('<center><b>\
                    Al<sub>x</sub>Ga<sub>1-x</sub>Sb\
                    </b></center')

        elif substrateType == 'GaN':
            #  self.input_substrate(self.qclayers.substrate)
            QMessageBox.information(self, 'ErwinJr Error', 
                    'III-Nitride substrates have not yet been implemented.')
            self.inputSubstrateBox.setCurrentIndex(
                    self.inputSubstrateBox.findText(self.qclayers.substrate))
            return

        else:
            raise TypeError('substrate selection not allowed')
            return

        if not self.updating:
            self.update_Lp_limits()
            self.update_inputBoxes()
            self.layerTable_refresh()
            self.qclayers.populate_x()
            self.update_quantumCanvas()

    @settingslot
    def input_EField(self, *args):
        """ SLOT connected to inputEFieldBox.valueChanged(double)
        update external E field in unit kV/cm """
        self.qclayers.EField = float(self.inputEFieldBox.value())
        self.qclayers.populate_x()
        self.update_quantumCanvas()
        self.emit(SIGNAL('dirty'))

    @settingslot
    def input_horzRes(self, *args):
        """ SLOT connected to inputHorzResBox.currentIndexChanged(int)
        update position resolution (xres), in angstrom """
        horzRes = unicode(self.inputHorzResBox.currentText())
        horzRes = float(horzRes)
        self.qclayers.set_xres(horzRes)
        self.qclayers.populate_x()
        self.update_quantumCanvas()
        self.emit(SIGNAL('dirty'))

    @settingslot
    def input_vertRes(self, *args):
        """ SLOT connected to inputVertResBox.valueChanged
        Update initial energy resolution for eigensolver. Set this too small
        may result in loss of some eigenvalue. """
        self.qclayers.vertRes = float(self.inputVertResBox.value())
        self.emit(SIGNAL('dirty'))

    @settingslot
    def input_repeats(self, *args):
        """ SLOT connected to SINGAL self.inputRepeatsBox.valueChanged(int)
        update number of repeats for the whole structure."""
        self.qclayers.repeats = int(self.inputRepeatsBox.value())
        self.qclayers.populate_x()
        self.update_quantumCanvas()
        self.emit(SIGNAL('dirty'))

    def input_basis(self):
        """ SLOT connected to self.inputARInjectorCheck.stateChanged(int) and
        self.inputInjectorARCheck.stateChanged(int)
        update dividers info
        """
        self.qclayers.basisARInjector = self.inputARInjectorCheck.isChecked()
        self.qclayers.basisInjectorAR = self.inputInjectorARCheck.isChecked()
        self.emit(SIGNAL('dirty'))

    def update_Lp_limits(self):
        """
        Update Lp select range in the Period Info box (GUI)
        """
        self.LpFirstSpinbox.setRange(1,self.qclayers.layerWidth.size-1)
        self.LpFirstSpinbox.setValue(1)
        self.LpLastSpinbox.setRange(1,self.qclayers.layerWidth.size-1)
        self.LpLastSpinbox.setValue(self.qclayers.layerWidth.size-1)

    def update_Lp_box(self):
        """
        Update Lp box in the Period Info box (GUI): 
            Lp:total length
            well: persentage of well material 
            nD: average doping (cm-3)
            ns: 2D carrier density in 1E11 cm-2
        """
        LpFirst = self.LpFirstSpinbox.value()
        LpLast = self.LpLastSpinbox.value()+1 
            #+1 because range is not inclusive of last value
        # total length of the layers (1 period)
        Lp = sum(self.qclayers.layerWidth[LpFirst:LpLast]) *self.qclayers.xres
        Lp_string  = u"Lp: %g \u212B<br>" % Lp
        # total length of well (1 period)
        Lw = sum((1-self.qclayers.layerBarriers[LpFirst:LpLast])
                *self.qclayers.layerWidth[LpFirst:LpLast]) *self.qclayers.xres
        if Lp == 0: 
            Lp_string += u"wells: NA%%<br>" 
            # average doping of the layers
            Lp_string += (u"n<sub>D</sub>: NA\u00D710<sup>17</sup>"
                    u"cm<sup>-3</sup><br>") 
        else: 
            Lp_string += u"wells: %6.1f%%<br>" % (100.0*Lw/Lp)
            # average doping of the layers
            nD = self.qclayers.xres * sum(
                    self.qclayers.layerDopings[LpFirst:LpLast]
                    *self.qclayers.layerWidth[LpFirst:LpLast])/Lp
            Lp_string += (u"n<sub>D</sub>: %6.3f\u00D710<sup>17</sup>"
                    u"cm<sup>-3</sup><br>") % nD
        # 2D carrier density in 1E11cm-2
        ns = self.qclayers.xres * sum(
                self.qclayers.layerDopings[LpFirst:LpLast]
                *self.qclayers.layerWidth[LpFirst:LpLast])*1e-2
        Lp_string += (u"n<sub>s</sub>: %6.3f\u00D710<sup>11</sup>"
            u"cm<sup>-2</sup") % ns
        self.LpStringBox.setText(Lp_string)


    @settingslot
    def input_description(self):
        """ SLOT connected to self.DescriptionBox.textChanged()
        Change description string for the design"""
        self.qclayers.description = self.DescriptionBox.toPlainText()
        self.emit(SIGNAL('dirty'))

    @settingslot
    def input_moleFrac(self, boxID):
        """ SLOT connected to self.MoleFracWellBox[n].editingFinished() 
            and MoleFracBarrBox[n].editingFinished()
        Update moleFrac for material."""
        self.qclayers.moleFrac[boxID] = float(
                self.MoleFracWellBox[boxID//2].value() if boxID % 2 == 0
                else self.MoleFracBarrBox[(boxID-1)//2].value())
        self.emit(SIGNAL('dirty'))

        self.update_inputBoxes()
        self.qclayers.update_alloys()
        self.qclayers.update_strain()
        self.qclayers.populate_x()
        self.update_quantumCanvas()

    @settingslot
    def set_temperature(self, T): 
        self.qclayers.Temperature = T
        self.qclayers.populate_x()
        self.qclayers.populate_x_band()
        self.update_quantumCanvas()


#============================================================================
# Layer Table Control
#============================================================================
    def layerNum(self):
        return self.qclayers.layerWidth.size

    def layerTable_refresh(self):
        """Refresh layer table, called every time after data update"""
        # Block itemChanged SIGNAL while refreshing
        #  self.clear_WFs()
        self.layerTable.blockSignals(True) 
        self.layerTable.clear()
        self.layerTable.setColumnCount(6)
        self.layerTable.setRowCount(self.qclayers.layerWidth.size+1)
        self.layerTable.setHorizontalHeaderLabels(['Width', 'ML', 'Brr', 
            'AR', 'Doping', 'Material'])
        #  vertLabels = []
        #  for n in xrange(self.qclayers.layerWidth.size+1):
            #  vertLabels.append(str(n))
        vertLabels = [str(n) for n in
                range(self.qclayers.layerWidth.size+1)]
        self.layerTable.setVerticalHeaderLabels(vertLabels)        

        #color for barrier layers
        gray = QColor(230,230,240)  # for Barrier layers
        gray2 = QColor(230,230,230) # for unchangable background

        for q, layerWidth in enumerate(self.qclayers.layerWidth):
            #Width Setup
            width = QTableWidgetItem("%5.1f" %
                    (layerWidth*self.qclayers.xres))
            width.setTextAlignment(Qt.AlignCenter)
            if bool(self.qclayers.layerBarriers[q]):
                width.setBackgroundColor(gray)
            self.layerTable.setItem(q, 0, width)
            if q == 0:
                width.setFlags(Qt.NoItemFlags)
                width.setBackgroundColor(gray2)

            #ML Setup
            numML = self.qclayers.xres*layerWidth/self.qclayers.MLThickness[q]
            item = QTableWidgetItem("%5.1f" % numML)
            item.setTextAlignment(Qt.AlignCenter)
            if bool(self.qclayers.layerBarriers[q]):
                item.setBackgroundColor(gray)
            self.layerTable.setItem(q, 1, item)
            if q == 0:
                item.setFlags(Qt.NoItemFlags)
                item.setBackgroundColor(gray2)

            #Barrier Layer Setup
            item = QTableWidgetItem()
            #  item.setCheckState(int(self.qclayers.layerBarriers[q])*2)
            item.setCheckState(Qt.Checked if
                    self.qclayers.layerBarriers[q]==1 else Qt.Unchecked)
            if bool(self.qclayers.layerBarriers[q]):
                item.setBackgroundColor(gray)
            self.layerTable.setItem(q, 2, item)
            if q == 0:
                item.setFlags(Qt.NoItemFlags)
                item.setBackgroundColor(gray2)

            #Active Region Layer Setup
            item = QTableWidgetItem()
            #  item.setCheckState(int(self.qclayers.layerARs[q])*2)
            item.setCheckState(Qt.Checked if
                    self.qclayers.layerARs[q]==1 else Qt.Unchecked)
            if bool(self.qclayers.layerBarriers[q]):
                item.setBackgroundColor(gray)
            self.layerTable.setItem(q, 3, item)
            if q == 0:
                item.setFlags(Qt.NoItemFlags)
                item.setBackgroundColor(gray2)

            #Layer Doping Setup
            doping = QTableWidgetItem(unicode(self.qclayers.layerDopings[q]))
            doping.setTextAlignment(Qt.AlignCenter)
            if bool(self.qclayers.layerBarriers[q]):
                doping.setBackgroundColor(gray)
            self.layerTable.setItem(q, 4, doping)
            if q == 0:
                doping.setFlags(Qt.NoItemFlags)
                doping.setBackgroundColor(gray2)

            #Material Setup
            if q == 0:
                item = QTableWidgetItem(unicode(self.materialList[
                    int(self.qclayers.layerMaterials[q])-1]))
                #TODO: reformat layerMaterials to int begin at 0
                item.setBackgroundColor(gray2)
                item.setFlags(Qt.NoItemFlags)
                self.layerTable.setItem(q, 5, item)
            else:
                materialWidget = QComboBox()
                materialWidget.addItems(self.materialList)
                materialWidget.setCurrentIndex(
                        self.qclayers.layerMaterials[q]-1)
                self.connect(materialWidget, 
                        SIGNAL("currentIndexChanged(int)"), 
                        partial(self.layerTable_materialChanged, q))
                self.layerTable.setCellWidget(q, 5, materialWidget)

        self.layerTable.resizeColumnsToContents()

        self.layerTable.blockSignals(False)


    def insert_layerAbove(self):
        """ SLOT connected to self.insertLayerAboveButton.clicked()"""
        row = self.layerTable.currentRow()
        if row == -1:
            return

        self.qclayers.layerWidth = np.insert(
                self.qclayers.layerWidth, row, 0)
        self.qclayers.layerBarriers = np.insert(
                self.qclayers.layerBarriers, row,
                0 if self.qclayers.layerBarriers[row] == 1 else 1)
        self.qclayers.layerARs = np.insert(
                self.qclayers.layerARs, row, 
                self.qclayers.layerARs[row])
        self.qclayers.layerMaterials = np.insert(
                self.qclayers.layerMaterials, row,
                self.qclayers.layerMaterials[row])
        self.qclayers.layerDopings = np.insert(
                self.qclayers.layerDopings, row,
                self.qclayers.layerDopings[row])
        self.qclayers.layerDividers = np.insert(
                self.qclayers.layerDividers, row,
                self.qclayers.layerDividers[row])

        self.update_Lp_limits()
        self.update_Lp_box()

        self.layerTable_refresh()
        self.qclayers.populate_x()
        self.layerTable.selectRow(row)
        self.layerTable.setFocus()

        self.emit(SIGNAL('dirty'))


    def delete_layer(self):
        """ SLOT connected to self.deleteLayerButton.clicked()"""
        #don't delete last layer
        if self.qclayers.layerWidth.size == 1:
            return
        row = self.layerTable.currentRow()
        if row == -1 or row >= self.qclayers.layerWidth.size:
            return

        self.qclayers.layerWidth = np.delete(
                self.qclayers.layerWidth, row)
        self.qclayers.layerBarriers = np.delete(
                self.qclayers.layerBarriers, row)
        self.qclayers.layerARs = np.delete(
                self.qclayers.layerARs, row)
        self.qclayers.layerMaterials = np.delete(
                self.qclayers.layerMaterials, row)
        self.qclayers.layerDopings = np.delete(
                self.qclayers.layerDopings, row)
        self.qclayers.layerDividers = np.delete(
                self.qclayers.layerDividers, row)

        if row == self.qclayers.layerWidth.size: #if row == last_row
            #make first item the same as last item
            self.qclayers.layerWidth[0] = self.qclayers.layerWidth[-1]
            self.qclayers.layerBarriers[0] = self.qclayers.layerBarriers[-1]
            self.qclayers.layerARs[0] = self.qclayers.layerARs[-1]
            self.qclayers.layerMaterials[0] = self.qclayers.layerMaterials[-1]
            self.qclayers.layerDopings[0] = self.qclayers.layerDopings[-1]
            self.qclayers.layerDividers[0] = self.qclayers.layerDividers[-1]

        self.update_Lp_limits()
        self.update_Lp_box()

        self.qclayers.update_strain()
        self.layerTable_refresh()
        self.qclayers.populate_x()
        self.layerTable.selectRow(row)
        self.layerTable.setFocus()

        self.emit(SIGNAL('dirty'))


    def layerTable_itemChanged(self, item):
        """SLOT connected to layerTable.itemChanged(QTableWidgetItem*)
        Update layer profile after user input"""
        #TODO: redo illegal input
        #  column = self.layerTable.currentColumn()
        #  row = self.layerTable.currentRow()
        #  print "---debug, itemChanged--- (%d, %d)"%(column, row)
        #  print "--debug, itemChanged (%d, %d)"%(item.column(), item.row())
        #  print item.text()
        #  if column == -1: #column == -1 on GUI initialization
            #  return
        column = item.column()
        row = item.row()
        if column == 0: #column == 0 for item change in Widths column
            new_width = float(item.text())
            new_width_int = int(np.round(new_width/self.qclayers.xres))
            #  if np.mod(new_width, self.qclayers.xres) != 0 \
                    #  and self.qclayers.xres != 0.1:
            if np.abs(new_width_int * self.qclayers.xres-new_width) > 1E-9:
                # TODO: bug to fix, np.mod is not good for xres < 0.5
                # potential solution is to change internal length to int
                # times xres
                QMessageBox.warning(self,"ErwinJr - Warning", 
                        ("You entered a width that is not compatible with "
                        "the minimum horizontal resolution. "
                        "%f %% %f = %f"%(new_width, self.qclayers.xres,
                            np.mod(new_width, self.qclayers.xres))))
                return
            if row == self.qclayers.layerWidth.size: #add row at end of list
                self.qclayers.layerWidth = np.append(
                        self.qclayers.layerWidth, new_width_int)
                self.qclayers.layerBarriers = np.append(
                        self.qclayers.layerBarriers, 
                        0 if self.qclayers.layerBarriers[-1] == 1 else 1)
                self.qclayers.layerARs = np.append(
                        self.qclayers.layerARs, 
                        self.qclayers.layerARs[-1])
                self.qclayers.layerMaterials = np.append(
                        self.qclayers.layerMaterials, 
                        self.qclayers.layerMaterials[-1])
                self.qclayers.layerDopings = np.append(
                        self.qclayers.layerDopings, 
                        self.qclayers.layerDopings[-1])
                self.qclayers.layerDividers = np.append(
                        self.qclayers.layerDividers, 
                        self.qclayers.layerDividers[-1])
                row += 1 #used so that last (blank) row is again selected

                #make first item the same as last item
                self.qclayers.layerWidth[0] = self.qclayers.layerWidth[-1]
                self.qclayers.layerBarriers[0] = self.qclayers.layerBarriers[-1]
                self.qclayers.layerARs[0] = self.qclayers.layerARs[-1]
                self.qclayers.layerMaterials[0] = self.qclayers.layerMaterials[-1]
                self.qclayers.layerDopings[0] = self.qclayers.layerDopings[-1]
                self.qclayers.layerDividers[0] = self.qclayers.layerDividers[-1]
                self.update_Lp_limits()

            elif row == self.qclayers.layerWidth.size-1:
                self.qclayers.layerWidth[row] = new_width_int
                #make first item the same as last item
                self.qclayers.layerWidth[0] = self.qclayers.layerWidth[-1]
                #  self.qclayers.layerBarriers[0] = self.qclayers.layerBarriers[-1]
                #  self.qclayers.layerARs[0] = self.qclayers.layerARs[-1]
                #  self.qclayers.layerMaterials[0] = self.qclayers.layerMaterials[-1]
                #  self.qclayers.layerDopings[0] = self.qclayers.layerDopings[-1]
                #  self.qclayers.layerDividers[0] = self.qclayers.layerDividers[-1]  

            else: #change Width of selected row in-place
                self.qclayers.layerWidth[row] = new_width_int

        elif column == 1: #column == 1 for ML
            if self.qclayers.xres != 0.1:
                QMessageBox.warning(self,"ErwinJr - Warning", 
                        (u"Horizontal Resolution of 0.1 \u212B required" 
                        u"when setting monolayer thicknesses."))
                return
            if row == self.qclayers.layerWidth.size: #add row at end of list
                pass
            elif row == self.qclayers.layerWidth.size-1:
                self.qclayers.layerWidth[row] = int(np.round( 
                    self.qclayers.MLThickness[row] * float(item.text())
                    / self.qclayers.xres))

                #make first item the same as last item
                self.qclayers.layerWidth[0] = self.qclayers.layerWidth[-1]
                #  self.qclayers.layerBarriers[0] = self.qclayers.layerBarriers[-1]
                #  self.qclayers.layerARs[0] = self.qclayers.layerARs[-1]
                #  self.qclayers.layerMaterials[0] = self.qclayers.layerMaterials[-1]
                #  self.qclayers.layerDopings[0] = self.qclayers.layerDopings[-1]
                #  self.qclayers.layerDividers[0] = self.qclayers.layerDividers[-1]

                self.update_Lp_limits()

            else: #change Width of selected row in-place
                self.qclayers.layerWidth[row] = int(np.round(
                        self.qclayers.MLThickness[row] * float(item.text()) 
                        / self.qclayers.xres ))
        elif column == 2: #column == 2 for item change in Barrier column
            if row == self.qclayers.layerWidth.size: 
                #don't do anything if row is last row
                return
            #  self.qclayers.layerBarriers[row] = int(item.checkState())//2
            self.qclayers.layerBarriers[row] = (item.checkState() == Qt.Checked)
            if row == self.qclayers.layerWidth.size-1:
                self.qclayers.layerBarriers[0] = self.qclayers.layerBarriers[-1]

        elif column == 3: #column == 3 for item change in AR column
            if row == self.qclayers.layerWidth.size: 
                #don't do anything if row is last row
                return
            #  self.qclayers.layerARs[row] = int(item.checkState())//2
            self.qclayers.layerARs[row] = (item.checkState() == Qt.Checked)
            if row == self.qclayers.layerWidth.size-1:
                self.qclayers.layerARs[0] = self.qclayers.layerARs[-1]

        elif column == 4: #column == 4 for item change in Doping column
            if row == self.qclayers.layerWidth.size: 
                #don't do anything if row is last row
                return
            self.qclayers.layerDopings[row] = float(item.text())
            if row == self.qclayers.layerWidth.size-1:
                self.qclayers.layerDopings[0] = self.qclayers.layerDopings[-1]

        elif column == 5: #column == 5 for item change in Materials column
            # See layerTable_materialChanged for more information
            #self.qclayers.layerWidth[row] = int(item.text[row])
            if row == self.qclayers.layerWidth.size-1:
                self.qclayers.layerMaterials[0] = self.qclayers.layerMaterials[-1]
        else:
            pass
        self.layerTable_refresh()
        #self.qclayers.populate_x()
        #self.layerTable.selectRow(row)
        self.layerTable.setCurrentCell(row,column)
        self.layerTable.setFocus()

        self.update_Lp_box()

        self.emit(SIGNAL('dirty'))


    def layerTable_itemSelectionChanged(self):
        """SLOT connected to layerTable.itemSelectionChanged()"""
        #This is the primary call to update_quantumCanvas
        self.qclayers.layerSelected = self.layerTable.currentRow()
        if not self.updating and self.qclayers.layerSelected >= 0 and \
                self.qclayers.layerSelected < self.qclayers.layerWidth.size:
            self.qclayers.populate_x()
            self.update_quantumCanvas()


    def layerTable_materialChanged(self, row, selection):
        """SLOT as partial(self.layerTable_materialChanged, q)) connected to 
        materialWidget.currentIndexChanged(int) """
        self.qclayers.layerMaterials[row] = selection+1
        #self.layerTable_refresh()
        self.qclayers.populate_x()
        self.layerTable.selectRow(row)

        self.emit(SIGNAL('dirty'))


    def bump_first_layer(self):
        """Move zeroth layer to first layer"""
        self.qclayers.layerWidth = np.insert(self.qclayers.layerWidth, 
                0, self.qclayers.layerWidth[-1])
        self.qclayers.layerBarriers = np.insert(self.qclayers.layerBarriers, 
                0, self.qclayers.layerBarriers[-1])
        self.qclayers.layerARs = np.insert(self.qclayers.layerARs, 
                0, self.qclayers.layerARs[-1])
        self.qclayers.layerMaterials = np.insert(self.qclayers.layerMaterials, 
                0, self.qclayers.layerMaterials[-1])
        self.qclayers.layerDopings = np.insert(self.qclayers.layerDopings, 
                0, self.qclayers.layerDopings[-1])
        self.qclayers.layerDividers = np.insert(self.qclayers.layerDividers, 
                0, self.qclayers.layerDividers[-1])

        self.update_inputBoxes()
        self.layerTable_refresh()
        self.layerTable.setCurrentCell(1,0)
        self.layerTable.setFocus()
        self.update_quantumCanvas()
        self.emit(SIGNAL('dirty'))


    def copy_structure(self):
        clipboard = QApplication.clipboard()
        string = ''
        for layer in self.qclayers.layerWidth[1:]:
            string += '%g\n' % (layer*self.qclayers.xres)
        clipboard.setText(string)


    def OptimizeLayer(self, goal):
        """ SLOT connected to SINGAL 
                self.OptimizeFomButton/OptimizeDipoleButton.clicked()
        Optimize the thickness of selected layer to maximaze the goal
        function using Newton's method.. 
        Support only for whole solve """
        row = self.layerTable.currentRow()
        if row == -1 or row >= self.qclayers.layerWidth.size:
            QMessageBox.warning(self, "ErwinJr Error", 
                "Invalid layer selection.")
            return

        self.Calculating(True)

        try:
            xres = self.qclayers.xres
            step = 1 # * xres
            upper = self.stateHolder[1]
            lower = self.stateHolder[0]        
            old_width = -1
            origin_width = new_width = self.qclayers.layerWidth[row]
            if DEBUG >= 1:
                print "--debug-- width optimization"

            goals = np.empty(3)
            goal_old = goals[1] = np.abs(goal(upper,lower))
            width_tried = [origin_width]
            goal_tried = [goal_old]
            #  while abs(old_width - new_width) >= 0.7*xres:
            while old_width != new_width :
                # Solve for values of goal near old_width
                # improve: only solve for eigen states near selection
                goal_old = goals[1]
                self.qclayers.layerWidth[row] = new_width - step
                self.qclayers.populate_x()
                self.qclayers.populate_x_band()
                self.qclayers.solve_psi()
                goals[0] = np.abs(goal(upper,lower))

                self.qclayers.layerWidth[row] = new_width + step
                self.qclayers.populate_x()
                self.qclayers.populate_x_band()
                self.qclayers.solve_psi()
                goals[2] = np.abs(goal(upper,lower))
                diff = (goals[2] - goals[0])/2
                diff2 = goals[0] + goals[2] - 2*goals[1]

                step_cutoff = 0.5E3 
                old_width = new_width
                # set a cutoff s.t. Newton's method won't go too far
                if -diff2 < 1/step_cutoff:
                    # When Newton's method is not a good one
                    new_width += int(np.round(step *
                        step_cutoff*diff/goals[1]))
                else:
                    new_width += -int(np.round(step * diff/diff2))
                if new_width <= 0:
                    new_width = 1
                if DEBUG >= 1:
                    print "Layer # %d width = %.1f; goal = %f"%(
                            row, old_width*xres, goals[1])
                    print "\tdiff = %f; diff2 = %f, new_width= %.1f"%(
                            diff, diff2, new_width*xres)

                self.qclayers.layerWidth[row] = new_width
                self.qclayers.populate_x()
                self.qclayers.populate_x_band()
                self.qclayers.solve_psi()
                goal_new = np.abs(goal(upper,lower))
                E_i = self.qclayers.EigenE[upper]
                E_j = self.qclayers.EigenE[lower]
                wavelength = h*c0/(e0*np.abs(E_i-E_j))*1e6
                if DEBUG >= 1:
                    print "\tgoal_new = %f, wl = %.1f um"%(goal_new, wavelength)

                while goal_new < goal_old*0.95: 
                    #  So a step will not go too far
                    #  new_width = (old_width + new_width)/2
                    new_width = int( (old_width + new_width)/(2))
                    if DEBUG >= 1:
                        print "\tGoing too far, back a little bit: "
                        print "\tnew_width=%.1f"%(new_width*xres)
                    self.qclayers.layerWidth[row] = new_width
                    self.qclayers.populate_x()
                    self.qclayers.populate_x_band()
                    self.qclayers.solve_psi()
                    goal_new = np.abs(goal(upper,lower))
                    E_i = self.qclayers.EigenE[upper]
                    E_j = self.qclayers.EigenE[lower]
                    wavelength = h*c0/(e0*np.abs(E_i-E_j))*1e6
                    if DEBUG >= 1:
                        print "\tgoal_new = %f, wl = %.1f um"%(goal_new, wavelength)

                goal_old = goals[1]
                goals[1] = goal_new
                if new_width in width_tried:
                    print "new_width has been tried"
                    break
                width_tried.append(new_width)
                goal_tried.append(goal_new)

            self.qclayers.layerWidth[row] = new_width
        finally:
            self.Calculating(False)
            if self.qclayers.layerWidth[row] != origin_width: 
                self.clear_WFs()
                self.layerTable_refresh()
                self.qclayers.populate_x()
                self.qclayers.populate_x_band()
                self.emit(SIGNAL('dirty'))
                self.plotDirty = True
                self.update_quantumCanvas()
        if DEBUG >= 1:
            print "done"


#==========================================================================
# Quantum Tab Plotting and Plot Control
#==========================================================================
    def update_quantumCanvas(self):
        #  print "update "+inspect.stack()[1][3]
        #  print inspect.stack()[2][3]
        # TODO: check plotDirty usage
        if self.plotDirty: 
            self.quantumCanvas.clear()
        self.quantumCanvas.axes.plot(self.qclayers.xPoints, 
                self.qclayers.xVc, 'k', linewidth=1)

        # plot Conduction Band L-Valley/X-Valley, Light Hole Valence Bnad and
        # Spin-Orbit coupling Valence Band
        for bandFlag, xv, conf in ((self.plotVL, self.qclayers.xVL, 'g--'), 
                (self.plotVX, self.qclayers.xVX, 'm-.'), 
                (self.plotLH, self.qclayers.xVLH, 'k'),
                (self.plotSO, self.qclayers.xVSO, 'r--')):
            if bandFlag:
                self.quantumCanvas.axes.plot(self.qclayers.xPoints, xv,
                        conf, linewidth=1)

        # highlight selected layer & make AR layers bold
        self.quantumCanvas.axes.plot(self.qclayers.xPoints,
                self.qclayers.xARs, 'k', linewidth=1.5)
        if self.qclayers.layerSelected >= 0 and \
                self.qclayers.layerSelected < self.qclayers.layerWidth.size: 
            self.quantumCanvas.axes.plot(self.qclayers.xPoints,
                    self.qclayers.xLayerSelected, 'b', 
                    linewidth = 1.5 if self.qclayers.layerARs[
                        self.qclayers.layerSelected] == 1 else 1)

        if self.plotDirty and hasattr(self.qclayers, 'EigenE'):
            self.curveWF = []
            for n in range(self.qclayers.EigenE.size):
                if n in self.stateHolder:
                    curve, = self.quantumCanvas.axes.plot( 
                            self.qclayers.xPointsPost, 
                            self.qclayers.xyPsiPsi[:,n]+self.qclayers.EigenE[n], 
                            'k', linewidth=2)
                else:
                    curve, = self.quantumCanvas.axes.plot( 
                            self.qclayers.xPointsPost, 
                            self.qclayers.xyPsiPsi[:,n]+self.qclayers.EigenE[n], 
                            color = self.colors[np.mod(n, len(self.colors))])
                self.curveWF.append(curve)

        self.quantumCanvas.draw()
        self.plotDirty = False


    def well_select(self, event):
        """ callback registered in plotControl when it's in wellselect mode. 
        It's mpl_connect to button_release_event of quantumCanvas """
        if event.button == 1: # left button clicked
            x = event.xdata
            xLayerNum = np.argmin((self.qclayers.xPoints-x)**2)
            layerNum = self.qclayers.xLayerNums[xLayerNum]
            self.layerTable.selectRow(layerNum)
            self.layerTable.setFocus()
        #  if on:
            #  self.picker2.setEnabled(True)
            #  #  self.quantumCanvas.canvas().setCursor(Qt.PointingHandCursor)
        #  else:
            #  self.picker2.setEnabled(False)
            #  #  self.quantumCanvas.canvas().setCursor(Qt.ArrowCursor)


    def clear_WFs(self):
        self.plotDirty = False
        if hasattr(self.qclayers, 'EigenE'):
            delattr(self.qclayers,'EigenE')
        self.update_quantumCanvas()


#============================================================================
# Export Functions
#============================================================================
    def export_quantumCanvas(self, filename=None):
        self.plotControl.save_figure(
                "ErwinJr - Export Band Structure Image", 
                filename, u'png')

    def export_band_data(self, fname):
        np.savetxt(fname.split('.')[0] + '_CB' + '.csv', 
                np.column_stack([self.qclayers.xPoints,self.qclayers.xVc]), 
                delimiter=',')

        try: self.qclayers.xyPsiPsi
        except AttributeError: pass #band structure hasn't been solved yet
        else:
            xyPsiPsiEig = np.zeros(self.qclayers.xyPsiPsi.shape)
            for q in xrange(self.qclayers.EigenE.size):
                xyPsiPsiEig[:,q] = self.qclayers.xyPsiPsi[:,q] + \
                        self.qclayers.EigenE[q]
            np.savetxt(fname.split('.')[0] + '_States' + '.csv', 
                    np.column_stack([self.qclayers.xPointsPost, xyPsiPsiEig]), 
                    delimiter=',')


#============================================================================
# View Band Items
#============================================================================
    def view_Band(self, band):
        # TODO: combine following slot to this one
        self.bandFlags[band] = not self.bandFlags[band]
        self.plotDirty = True
        self.update_quantumCanvas()

    def view_VXBand(self):
        if self.plotVX:
            self.plotVX = False
        else:
            self.plotVX = True
        self.plotDirty = True
        self.update_quantumCanvas()

    def view_VLBand(self):
        if self.plotVL:
            self.plotVL = False
        else:
            self.plotVL = True
        self.plotDirty = True
        self.update_quantumCanvas()

    def view_LHBand(self):
        if self.plotLH:
            self.plotLH = False
        else:
            self.plotLH = True
        self.plotDirty = True
        self.update_quantumCanvas()

    def view_SOBand(self):
        if self.plotSO:
            self.plotSO = False
        else:
            self.plotSO = True
        self.plotDirty = True
        self.update_quantumCanvas()


#==========================================================================
# Calculations
#==========================================================================
    def Calculating(self, is_doing):
        """UI repaint for doing calculating """
        for button in (self.solveWholeButton, self.solveBasisButton,
                self.pairSelectButton):
            button.setEnabled(not is_doing)
            button.repaint()
        if self.pairSelected: 
            self.FoMButton.setEnabled(not is_doing)
            self.FoMButton.repaint()
            if self.solveType == 'whole':
                for button in (self.OptimizeFoMButton, self.OptimizeDipoleButton):
                    button.setEnabled(not is_doing)
                    button.repaint()


    def solve_whole(self):  #solves whole structure
        """SLOT connected to solveWholeButton.clicked()
        Whole solver """
        if hasattr(self.qclayers, "EigenE"):
            self.clear_WFs()
        self.pairSelected = False
        self.Calculating(True)

        self.qclayers.populate_x_band()
        try:
            self.stateHolder = []
            self.qclayers.solve_psi()
            self.plotDirty = True
            self.solveType = 'whole'
            self.update_quantumCanvas()
            if DEBUG >= 4: 
                with open('qclayer.pkl','wb') as f:
                    pickle.dump(self.qclayers, f, pickle.HIGHEST_PROTOCOL)
                print self.qclayers.EigenE
        except (IndexError,TypeError) as err:
            QMessageBox.warning(self, 'ErwinJr - Error', str(err))

        self.Calculating(False)


    def solve_basis(self):  #solves structure with basis
        """SLOT connected to solveBasisButton.clicked()
        Basis solver """
        self.Calculating(True)

        try:
            self.stateHolder = []
            self.dCL = self.qclayers.basisSolve()
            self.qclayers.convert_dCL_to_data(self.dCL)
            self.solveType = 'basis'        
            self.plotDirty = True
            self.update_quantumCanvas()
        except (ValueError,IndexError) as err:
            QMessageBox.warning(self,"ErwinJr - Error", str(err))

        self.Calculating(False)


    def state_pick(self, event):
        """ callback registered in plotControl when it's in pairselect mode. 
        It's mpl_connect to button_release_event of quantumCanvas """
        if not hasattr(self.qclayers, 'EigenE'):
            # Not yet solved 
            return

        if event.button == 1: # left button clicked: select a state
            if len(self.stateHolder) >= 2:
                # start new pair selection
                self.stateHolder = []
                self.pairSelected = False
                for button in (self.FoMButton, self.OptimizeFoMButton,
                        self.OptimizeDipoleButton):
                    # TODO: move all similiar codes to a button update
                    button.setEnabled(False)
                    button.repaint()

            # TODO: replace follwoing by matplotlib lines picker
            x = event.xdata
            y = event.ydata
            xData = np.tile(self.qclayers.xPointsPost,
                    (self.qclayers.xyPsiPsi.shape[1],1)).T
            yData = self.qclayers.xyPsiPsi + self.qclayers.EigenE

            width, height = self.quantumCanvas.axes.bbox.size
            xmin, xmax = self.quantumCanvas.axes.get_xlim()
            ymin, ymax = self.quantumCanvas.axes.get_ylim()
            xScale = (xmax - xmin)/width
            yScale = (ymax - ymin)/height

            r = np.nanmin( sqrt(
                ((xData-x)/xScale)**2+((yData-y)/yScale)**2), axis=0)
            ss = np.nanargmin(r)
            self.stateHolder.append(ss)
            self.plotDirty = True
            #  self.curveWF[ss].set_color('black')
            #  self.curveWF[ss].set_linewidth(2)
        elif event.button == 3: # right button clicked: remove last selection
            if len(self.stateHolder) == 2:
                self.pairSelected = False
                for button in (self.FoMButton, self.OptimizeFoMButton,
                        self.OptimizeDipoleButton):
                    button.setEnabled(False)
                    button.repaint()
            self.stateHolder.pop()
            self.plotDirty = True
        else:
            return
        self.update_quantumCanvas()

        if len(self.stateHolder) == 1:
            self.pairString  = (u"selected: %d, ..<br>"%self.stateHolder[0])
        elif len(self.stateHolder) == 2:
            self.pairSelected = True
            #TODO: put these enablement to a functions
            self.FoMButton.setEnabled(True)
            if self.solveType == 'whole':
                self.OptimizeFoMButton.setEnabled(True)
                self.OptimizeDipoleButton.setEnabled(True)
            E_i = self.qclayers.EigenE[self.stateHolder[0]]
            E_j = self.qclayers.EigenE[self.stateHolder[1]]
            if E_i > E_j:
                upper = self.stateHolder[0]
                lower = self.stateHolder[1]
            else:
                upper = self.stateHolder[1]
                lower = self.stateHolder[0]

            self.eDiff = 1000*(E_i-E_j)
            self.wavelength = h*c0/(e0*np.abs(E_i-E_j))*1e6

            if self.solveType is 'basis':
                couplingEnergy = self.qclayers.coupling_energy(
                        self.dCL, upper, lower)
                self.transitionBroadening = self.qclayers.broadening_energy(
                        upper, lower)
                self.qclayers.populate_x_band()
                self.opticalDipole = self.qclayers.dipole(upper, lower)            
                self.tauUpperLower = 1/self.qclayers.lo_transition_rate(
                        upper, lower)
                self.pairString  = (u"selected: %d, %d<br>" 
                        u"energy diff: <b>%6.1f meV</b> (%6.1f um)<br>" 
                        u"coupling: %6.1f meV<br>broadening: %6.1f meV<br>" 
                        u"dipole: <b>%6.1f \u212B</b>" 
                        u"<br>LO scattering: <b>%6.2g ps</b><br>") % ( 
                                self.stateHolder[0], 
                                self.stateHolder[1], 
                                self.eDiff, self.wavelength,
                                couplingEnergy, 
                                self.transitionBroadening, 
                                self.opticalDipole, 
                                self.tauUpperLower)

            elif self.solveType is 'whole':
                self.qclayers.populate_x_band()
                self.opticalDipole = self.qclayers.dipole(upper, lower)
                self.tauUpperLower = 1/self.qclayers.lo_transition_rate(
                        upper, lower)
                self.transitionBroadening = 0.1 * self.eDiff #TODO
                self.pairString = (u"selected: %d, %d<br>"
                        u"energy diff: <b>%6.1f meV</b> (%6.1f um)<br>"
                        u"dipole: %6.1f \u212B<br>" 
                        u"LO scattering: %6.2g ps<br>") % (
                                self.stateHolder[0], 
                                self.stateHolder[1], 
                                self.eDiff, self.wavelength, 
                                self.opticalDipole, 
                                self.tauUpperLower)
            else:
                self.FoMButton.setEnabled(False)

        self.pairSelectString.clear()
        self.pairSelectString.setText(self.pairString)        


    def updateFoM(self):
        """ SLOT connected to FoMButton.clicked()
        Calculate Figure of merit.  """
        if len(self.stateHolder) < 2:
            return
        self.Calculating(True)
        self.FoMButton.setEnabled(False)
        self.FoMButton.repaint()

        upper = self.stateHolder[1]
        lower = self.stateHolder[0]        
        if upper < lower:
            upper, lower = lower, upper

        self.tauLower = self.qclayers.lo_life_time(lower)
        self.tauUpper = self.qclayers.lo_life_time(upper)
        self.FoM = self.opticalDipole**2 * self.tauUpper \
                * (1- self.tauLower/self.tauUpperLower)
        # tauUpperLower is the inverse of transition rate (lifetime)
        self.alphaISB = self.qclayers.alphaISB(upper, lower)

        self.FoMString  = (u"<i>\u03C4<sub>upper</sub></i> : %6.2f ps<br>"
                u"<i>\u03C4<sub>lower</sub></i> : %6.2f ps"
                u"<br>FoM: <b>%6.0f ps \u212B<sup>2</sup></b>"
                u"<br><i>\u03B1<sub>ISB</sub></i> : %.2f cm<sup>2</sup>") % (
                        self.tauUpper, self.tauLower, self.FoM, self.alphaISB)
        self.pairSelectString.setText(self.pairString + self.FoMString)

        self.Calculating(False)
        self.FoMButton.setEnabled(True)
        self.toOpticalParamsButton.setEnabled(True)


    def get_nCore(self, wl):
        """Get overall active core complex reflaction index (imag part being
        decay), by average over width. Used for optical mode calculation. 
        wl for wavelength"""
        n = np.empty(self.numMaterials)
        for q in range(self.numMaterials):
            n[q] = self.qclayers.moleFrac[q] *\
                        cst[self.qclayers.Mat1[q]].rIndx(wl)\
                    + (1-self.qclayers.moleFrac[q]) *\
                        cst[self.qclayers.Mat2[q]].rIndx(wl)
        # Average n?
        nCore = np.sum(self.qclayers.MaterialWidth*n)/np.sum(self.qclayers.MaterialWidth) 
        return nCore

    def transfer_params(self, strata):
        """ transfer parameters to strata. """
        #set wavelength
        strata.wavelength = 1.24/np.abs(self.eDiff)*1000
        #set operating field
        strata.operatingField = self.qclayers.EField
        #set Lp
        LpFirst = self.LpFirstSpinbox.value()
        LpLast = self.LpLastSpinbox.value()+1 
            #+1 because range is not inclusive of last value
        strata.Lp = self.qclayers.xres * np.sum(
                self.qclayers.layerWidth[LpFirst:LpLast])
        #set nD doping sheet density
        strata.nD = np.sum(self.qclayers.layerDopings[LpFirst:LpLast] *
                self.qclayers.layerWidth[LpFirst:LpLast]) / \
                np.sum(self.qclayers.layerWidth[LpFirst:LpLast])
        #set aCore
        strata.aCore = self.alphaISB
        #set nCore
        kCore = 1/(4*pi) * strata.aCore * strata.wavelength*1e-4 
        # See Def of acore
        # 1e-4: aCore in cm-1, wl in um
        strata.nCore = self.get_nCore(strata.wavelength) + 1j*kCore

        #set tauUpper
        strata.tauUpper = self.tauUpper
        #set tauLower
        strata.tauLower = self.tauLower
        #set tauUpperLower
        strata.tauUpperLower = self.tauUpperLower
        #set optical dipole
        strata.opticalDipole = self.opticalDipole
        #set figure of merit
        strata.FoM = self.FoM
        #2gamma transition broadening
        strata.transitionBroadening = self.transitionBroadening / 1000 
            # store in eV


#==========================================================================
# Global Optimization
#===========================================================================
    def set_targetWL(self):
        """ SLOT connected to self.targetWL_box.editingFinished()
        To set target wavelength."""
        try:
            wl = float(self.targetWL_box.text())
        except ValueError:
            QMessageBox.warning(self, 'ErwinJr Error', 
                'Invalid input:%s'%(self.targetWL_box.text()))
            self.targetWL_box.setText('')
        self.targetWL = wl
        self.targetWL_box.setText('%.1f'%self.targetWL)


    def set_goal(self, goal): 
        """ SLOT connected to self.goalFuncBox.currentIndexChanged(int)
        To set goal function for global optimization."""
        self.OptGoal = self.OptGoalsFunc[goal]


    def GlobalOptimization(self):
        """ SLOT connected to self.GlobalOptButton.cliecked()
        To conduct global Optimization (TBD)."""
        if not hasattr(self, 'targetWL'):
            QMessageBox.warning(self, "ErwinJr Error", 
                "Target wavelength is not set.")
            return
        if DEBUG >= 1:
            print "Global Optimization for %s"%self.OptGoal.__name__
        Jaco = 0


# vim: ts=4 sw=4 sts=4 expandtab
