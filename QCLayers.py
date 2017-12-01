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
# material related codes should be moved to MaterialConstants
# try use dict type for substrate restriction on material
# Add solve well (helping tight binding model design
# Performance improve: lo_transition_rate and solve_psi:
#               parallel support: Done
#               reduce reductant calculation
#               Graphic card support?
# Try matrix eigen solver in solve_psi - It's O(n^3), compared to what we
# have for O(n^2) (or O(mn) where m for E resolution and n for x)
# replace CLIB by Cython

from __future__ import division

__LOG__ = False
__DEBUG__ = 1

__USE_CLIB__ = True
__MORE_INTERPOLATION__ = True  # One more time interpolation for eigen solver
__MULTI_PROCESSING__ = True

import copy
import sys
import numpy as np
from numpy import sqrt, exp, sin, cos, log, pi, conj, real, imag
from scipy import interpolate

from settings import (wf_scale, wf_min_height, pretty_plot_factor,
                      plot_decimate_factor, phonon_integral_factor)
import MaterialConstantsDict
cst = MaterialConstantsDict.MaterialConstantsDict()

if __LOG__:
    import pickle
    logcount = 0
if __DEBUG__ >= 1:
    from time import time

# TODO: replace CLIB by Cython
if __USE_CLIB__:
    from ctypes import *
    if __MULTI_PROCESSING__:
        if sys.platform in ('linux2', 'darwin', 'cygwin'):
            cQ = CDLL('./cQCLayersMP.so')
        elif sys.platform == 'win32':
            cQ = CDLL('cQCLayersMP.dll')
    else:
        if sys.platform in ('linux2', 'darwin', 'cygwin'):
            cQ = CDLL('./cQCLayers.so')
        elif sys.platform == 'win32':
            cQ = CDLL('cQCLayers.dll')

# ===========================================================================
# Global Variables
# ===========================================================================
from scipy.constants import (e as e0, epsilon_0 as eps0,
                             electron_mass as m0, c as c0)
from scipy.constants import h, hbar

ANG = 1e-10    # angstrom to meter
KVpCM = 1e5    # KV/cm to V/m
meV = 1e-3     # meV to eV

INV_INF = 1e-20  # for infinit small decay rate (ns-1)
PAD_HEAD = 100  # width padded in the head of the given region for basis solver
PAD_TAIL = 30

# ===========================================================================
# Reference
# [0]Kale Franz's thesis
# [1]Handbook of Optics, Vol.2, ISBN: 0070479747
# [2]Van de Walle C G. Band lineups and deformation potentials in the
#    model-solid theory[J]. Physical review B, 1989, 39(3): 1871.
# [3]Peter Qiang Liu's thesis
# ===========================================================================

# for In0.53Ga0.47As, EcG = 0.22004154
#    use this as a zero point baseline
bandBaseln = 0.22004154


class QCLayers(object):
    """Class for QCLayers
    Member variables:
        parameters for each layer, np.array type, with len = No. of layers:
            layerWidth -int in pixel (xres per pixel), width of each layer
            layerBarriers -boolean(TBD), if the layer is barrier or not
            layerARs -boolean(TBD), if the layer is active region or not
                      only affect basis solver (negelet some coupling
            layerMaterials -int(TBD), label of material, depending on
                        substrate, the material is defined in erwinjr.pyw
            layerDopings -Doping per volumn in unit 1e17 cm-3
            layerDividers -??? seems not used
        xres: position resolution, in angstrom
        vertRes: vertical/energy resolution, in meV
        EField: external (static) electrical field, in kV/cm = 1e5 V/m
        layerSelected: (for GUI) a label which layer is selected in GUI,
                        with default -1 indicating non selected.
        repeats: (int) is the number of repeat times for the given structure
        solver: ??? seems not used
        Temperature: Temperature of the device, affecting material property
                      seems not used
        TempFoM: ??? seems not used (Figure of Merit?
        diffLength: ??? seems not used
        substrate: The substrate material for the device, which determined
                      the well and barrier material
                      substrate   | well            | barrier
                      InP         | In_xGa_{1-x}As  | Al_{1-x}In_xAs
                      GaAs        | Al_xGa_{1-x}As  | Al_xGa_{1-x}As
                      GaSb        | InAs_ySb_{1-y}  | Al_xGa_{1-x}Sb
                      GaN  (TBD)
        basisARInjector & basisInjectorAR: where should the layer separate
                    for basis solver
        moleFrac: mole fraction for each possible layer material, in format
                    [well, barrier]*4
    """
    def __init__(self):
        self.layerWidth = np.array([1, 1])      # pix
        self.layerBarriers = np.array([0, 0])   # boolean
        self.layerARs = np.array([0, 0])        # boolean
        self.layerMaterials = np.array([1, 1])  # label
        self.layerDopings = np.array([0., 0.])  # 1e17 cm-3
        self.layerDividers = np.array([0, 0])   # Basis solver divider

        self.xres = 0.5          # angstrom
        self.EField = 0          # kV/cm = 1e5 V/m
        self.layerSelected = -1  # int
        self.vertRes = 0.5       # meV
        self.repeats = 2         # repeats n times for the given structure

        self.description = ""
        self.solver = "SolverH"  # ?
        self.Temperature = 300
        self.TempFoM = 300       # ?
        self.diffLength = 0      # ?
        self.basisARInjector = True
        self.basisInjectorAR = True
        #  self.designByAngs = True
        #  self.designByML = False
        self.substratesList = ('InP', 'GaAs', 'GaSb', 'GaN')
        self.substrate = 'InP'

        self.moleFrac = [0.53, 0.52, 0.53, 0.52, 0.53, 0.52, 0.53, 0.52]

        self.update_alloys()
        self.update_strain()
        self.populate_x()

    def set_xres(self, res):
        for n in range(self.layerWidth.size):
            self.layerWidth[n] = int(np.round(
                self.layerWidth[n] * self.xres / res))
        self.xres = res

    def populate_x(self):
        """Extend layer information to position functions
        Layer data: layerWidth
                with len = # of layers and each value repr. a layer
        Position data (OUTPUT/update member variables):
                    xPoints
                       position grid
                    xBarriers: from layerBarriers, is barrier layer
                        should be boolean (TBD)
                    xARs: from layerARs and xVc, xVc if is active region
                        (layerARs == 1) otherwise np.NaN
                    xMaterials: from layerMaterials label/index of material
                        should be int starting from 0 (TBD)
                    xDopings: from layerDopings, doping per volumn
                    xLayerNums
                       at xPoints[q] it's xLayerNums[q]-th layer
                    xLayerSelected: from layerSelected and xVc, xVc if it's
                        the layer indicated by layerSelected, otherwise NaN
        """
        #  print "-----debug----- QCLayers populate_x called"
        #  print self.layerBarriers
        self.xPoints = self.xres * np.arange(0, self.layerWidth.sum())

        layerNumCumSum = np.concatenate(([0], self.layerWidth.cumsum()))
        self.xBarriers = np.zeros(self.xPoints.shape)
        self.xARs = np.zeros(self.xPoints.shape)
        self.xMaterials = np.zeros(self.xPoints.shape)
        self.xDopings = np.zeros(self.xPoints.shape)
        self.xLayerNums = np.zeros(self.xPoints.shape)

        # extend layer data for all xpoints
        for q in xrange(0, self.layerWidth.size):
            self.xBarriers[layerNumCumSum[q]:
                           layerNumCumSum[q+1]] = self.layerBarriers[q]
            if self.layerARs[q] == 1:
                self.xARs[layerNumCumSum[q]-1:
                          layerNumCumSum[q+1]+1] = 1
            self.xMaterials[layerNumCumSum[q]:
                            layerNumCumSum[q+1]] = self.layerMaterials[q]
            self.xDopings[layerNumCumSum[q]:
                          layerNumCumSum[q+1]] = self.layerDopings[q]
            self.xLayerNums[layerNumCumSum[q]:
                            layerNumCumSum[q+1]] = q

        # duplicate layer based on user input repeats
        #  repeats = self.repeats
        if self.repeats >= 2:
            self.xPoints = np.arange(
                0, self.xres*(self.layerWidth.sum() +
                              self.layerWidth[1:].sum()*(self.repeats-1)),
                self.xres)
            self.xBarriers = np.concatenate((self.xBarriers, np.tile(
                self.xBarriers[layerNumCumSum[1]:], self.repeats-1)))
            self.xARs = np.concatenate((self.xARs, np.tile(
                self.xARs[layerNumCumSum[1]:], self.repeats-1)))
            self.xMaterials = np.concatenate((self.xMaterials, np.tile(
                self.xMaterials[layerNumCumSum[1]:], self.repeats-1)))
            self.xDopings = np.concatenate((self.xDopings, np.tile(
                self.xDopings[layerNumCumSum[1]:], self.repeats-1)))
            self.xLayerNums = np.concatenate((self.xLayerNums, np.tile(
                self.xLayerNums[layerNumCumSum[1]:], self.repeats-1)))

        # this hack is needed because sometimes
        # self.xPoints is one element too big
        self.xPoints = self.xPoints[0:self.xBarriers.size]

        self.update_strain()
        # Following are equiv. elec potential for different bands
        # external field is included
        # xVX, xVL, xVLH and xVSO are used for checking if there's indrect
        # bandgap, s.t. we can prevent its effect
        self.xVc = np.zeros(self.xPoints.size)
        self.xVX = np.zeros(self.xPoints.size)
        self.xVL = np.zeros(self.xPoints.size)
        self.xVLH = np.zeros(self.xPoints.size)
        self.xVSO = np.zeros(self.xPoints.size)
        for MLabel in range(1, 5):
            indx = np.nonzero(self.xMaterials == MLabel)[0]
            if indx.size != 0:
                material = np.where(self.xBarriers[indx] == 1,
                                    MLabel*2-1, (MLabel-1)*2)
                self.xVc[indx] = (self.EcG[material] - self.xPoints[indx] *
                                  ANG * self.EField * KVpCM)
                self.xVX[indx] = (self.EcX[material] - self.xPoints[indx] *
                                  ANG * self.EField * KVpCM)
                self.xVL[indx] = (self.EcL[material] - self.xPoints[indx] *
                                  ANG * self.EField * KVpCM)
                self.xVLH[indx] = (self.EvLH[material] - self.xPoints[indx] *
                                   ANG * self.EField * KVpCM)
                self.xVSO[indx] = (self.EvSO[material] - self.xPoints[indx] *
                                   ANG * self.EField * KVpCM)

        # make array to show selected layer in mainCanvas
        try:
            self.xLayerSelected = np.zeros(self.xPoints.shape)*np.NaN
            layerSelected = self.layerSelected
            if layerSelected != -1:
                if layerSelected == 0:
                    # row for first layer is selected
                    for repeat in xrange(1, self.repeats+1):
                        base = layerNumCumSum[-1] * (repeat - 1)
                        self.xLayerSelected[
                            base+layerNumCumSum[layerSelected]:
                            base+layerNumCumSum[layerSelected+1]+1
                        ] = self.xVc[
                            base+layerNumCumSum[layerSelected]:
                            base+layerNumCumSum[layerSelected+1]+1
                        ]
                elif self.layerSelected == self.layerWidth.size:
                    # last (blank) layer row is selected
                    pass
                else:
                    for repeat in xrange(1, self.repeats+1):
                        base = np.sum(self.layerWidth[1:])*(repeat-1)
                        self.xLayerSelected[
                            base+layerNumCumSum[layerSelected]-1:
                            base+layerNumCumSum[layerSelected+1]+1
                        ] = self.xVc[
                            base+layerNumCumSum[layerSelected]-1:
                            base+layerNumCumSum[layerSelected+1]+1
                        ]
        except IndexError:
            # index error happens in SolveBasis when the selected layer is
            # greater than the number of layers in the solve swath
            # however, xLayerSelected is not used for the SolveBasis function
            #  print ("Index Error for layer selection at function"
                #  "qclayer.populate_x")
            pass

        self.xARs[np.nonzero(self.xARs == 0)[0]] = np.NaN
        self.xARs *= self.xVc

    def populate_x_band(self):
        """Extend layer information to position functions for band parameter
        OUTPUT/update member variables):
            xEg, xMc, xESO, xEp, xF,
            whose value is determeined by the layer material
        TODO: tell why we need populate_x and populate_x_band
        """
        #  print "------debug------- QCLayers populate_x_band called"
        # Following parameters can be looked up in cQCLayers.c
        self.xEg = np.zeros(self.xPoints.size)
        self.xMc = np.zeros(self.xPoints.size)  # Seems not to be used
        self.xESO = np.zeros(self.xPoints.size)
        self.xEp = np.zeros(self.xPoints.size)
        self.xF = np.zeros(self.xPoints.size)
        for MLabel in range(1, 5):
            indx = np.nonzero(self.xMaterials == MLabel)[0]
            if indx.size != 0:
                material = np.where(self.xBarriers[indx] == 1,
                                    MLabel*2-1, (MLabel-1)*2)
                self.xEg[indx] = self.EgLH[material]
                self.xMc[indx] = self.me[material]
                self.xESO[indx] = self.ESO[material]
                self.xEp[indx] = self.Ep[material]
                self.xF[indx] = self.F[material]

    def update_alloys(self):  # c is a Material_Constant class instance
        """ update material parameter for the alloy used.
        (Always followed by update_strain)
        OUTPUT/update member variable:
            all parameters listed in variables: np.array with len=numMaterials
                        labeled by sequence [well, barrier]*4
            self.numMaterials: Number of differenc types of material
                               supported
            self.epsrho: ???
        """
        variables = [
            'EgG', 'EgL', 'EgX', 'VBO', 'DSO',  # unit eV
            'me0',                # seems not used
            'acG', 'acL', 'acX',  # Pikus-Bir interaction parameter
            'Ep', 'F',  # effective mass parameter, unit eV (Ep) and 1 (F)
            'XiX',      # strain correction to band at X point, unit eV
            'b', 'av', 'alG',  # strain correction to bands at Gamma, unit eV
            'beG', 'alL',      # Varsh correction
            # 'beL', 'alX', 'beX',  # seems not used
            'epss', 'epsInf',  # static and high-freq permitivity
            'hwLO',            # LO phonon energy, unit eV
            'alc',             # lattice const, unit angstrom
            'c11', 'c12'      # elestic stiffness constants
        ]
        #  print "----debug--- substrate is "+self.substrate
        # substrate restriction on layer material,
        # see doc string of QCLayers class
        # Material are labeled by sequence [well, barrier]*4
        if self.substrate == 'InP':
            self.numMaterials = 8
            self.Mat1 = ['InAs']*8
            self.Mat2 = ['GaAs', 'AlAs']*4
            MatCross = ['InGaAs', 'AlInAs']*4
        elif self.substrate == 'GaAs':
            self.numMaterials = 8
            self.Mat1 = ['AlAs']*8
            self.Mat2 = ['GaAs']*8
            MatCross = ['AlGaAs']*8  # Note EgG_AlGaAs's moleFrac deps
        elif self.substrate == 'GaSb':
            self.numMaterials = 8
            self.Mat1 = ['InAs', 'AlSb']*4
            self.Mat2 = ['InSb', 'GaSb']*4
            MatCross = ['InAsSb', 'AlGaSb']*4
            # Note EgG's bowing moleFrac deps
        else:
            raise TypeError('substrate selection not allowed')

        for item in variables:
            setattr(self, item, np.empty(self.numMaterials))
            para = getattr(self, item)
            for n in range(self.numMaterials):
                para[n] = self.moleFrac[n]*getattr(cst[self.Mat1[n]], item) \
                    + (1-self.moleFrac[n])*getattr(cst[self.Mat2[n]], item)
                if MatCross[n] in cst and hasattr(cst[MatCross[n]], item):
                    # bowing parameter
                    para[n] -= self.moleFrac[n]*(1-self.moleFrac[n]) \
                                * getattr(cst[MatCross[n]], item)

        # See MaterialConstantsDict.py...
        # TODO: move to MaterialConstantsDict.py
        if self.substrate == 'GaAs':
            for n in range(self.numMaterials):
                EgG_AlGaAs = cst['AlGaAs'].EgG + 1.310*self.moleFrac[n]
                self.EgG[n] = self.moleFrac[n]*cst['AlAs'].EgG \
                    + (1-self.moleFrac[n])*cst['GaAs'].EgG \
                    - self.moleFrac[n]*(1-self.moleFrac[n])*EgG_AlGaAs
        elif self.substrate == 'GaSb':
            for n in range(1, self.numMaterials, 2):
                EgG_AlGaSb = cst['AlGaSb'].EgG + 1.22*self.moleFrac[n]
                self.EgG[n] = self.moleFrac[n]*cst['AlSb'].EgG \
                    + (1-self.moleFrac[n])*cst['GaSb'].EgG \
                    - self.moleFrac[n]*(1-self.moleFrac[n])*EgG_AlGaSb
        #  print "----debug----"
        #  for item in variables:
            #  ll = copy.copy(getattr(self,item))
            #  #  print item, getattr(self, item)
            #  setattr(self, item+"_new", ll)

        # set this once the others are set ???
        self.epsrho = 1 / (1/self.epsInf - 1/self.epss)

    def update_strain(self):  # c is a Material_Constant class instance
        """Update strain and strain related parameters inside each layers
        (Always called after update_alloys)
        OUTPUT/update member variables:
            (all below are np.array with len=numMaterials)
            (Material are labeled by sequence [well, barrier]*4)
            self.a_parallel: lattice const. within/parallel to the layer plane
            self.eps_parallel: strain tensor within/parallel to the layer plane
            self.a_perp: lattice const. perpendicular to the layer plane
            self.eps_perp: strain tensor perpendicular to the layer plane
            self.MaterialWidth: total width of a each material
            self.netStrain: spacial average of eps_perp in unit of percentage
            self.avghwLO: spacial average of hwLO in unit of eV
            self.MLThickness: monolayer thickness? shown in GUI as
                    layerWidth/MLThickness. This is an average number of
                    layers and actually the edge is rough
            self.Pec, self.Pe, self.Qe, self.Varsh: correction terms on bands,
                                See Kales's thesis, sec2
            self.ESO: spin-orbit splitting, including strain correction
            self.EgLH, self.EgSO: band bottom/top at Gamma Epoints respect to
                                conduction band
            self.me: effective mass ignoring energy dependence, in unit m0
                        only used as a backup and for broadening_energy
                        may be deleted for further update
            self.EcG, self.EcL, self.EcX: conduction band bottom at
                                Gamma, L and X points, respect to a give
                                baseline
            self.EvLH, self.EvSO: valence band (LH/SO) top at Gamma point
            (EcL, EcX, EvLH, EvSO are only used for plotting?)
        """
        if self.substrate in cst.substrateSet:
            self.a_parallel = cst[self.substrate].alc
            # parallel littice constant depends on substrate
        else:
            raise TypeError('substrate selection not allowed')

        # [2]Walle eqn 1b
        self.eps_parallel = self.a_parallel / self.alc - 1
        # [2]Walle eqn 2a and 4a
        self.a_perp = self.alc * (1 - 2 * self.c12 / self.c11 *
                                  self.eps_parallel)
        # [2]Walle eqn 2b
        self.eps_perp = self.a_perp/self.alc - 1
        #             = -2*self.c12/self.c11*self.eps_parallel

        # total width of different material
        self.MaterialWidth = np.zeros(self.numMaterials)
        for i in range(4):
            # Note that material are labeled by sequence [well, barrier]*4
            # [1:] because first layer doesn't count
            #  indx = np.nonzero(self.layerMaterials[1:] == i+1)[0]
            #  # print i+1, self.layerMaterials
            # (BUG FIXED: self.layerMaterials[1:] results in material index
            # mismatch by 1)
            # self.layerWidth includes an extra layer to promise first=last
            indx = (self.layerMaterials == i+1)
            indx[0] = False  # s.t. 1st layer doesn't count
            self.MaterialWidth[2*i+1] = self.xres*np.sum(
                self.layerWidth[indx] * self.layerBarriers[indx])
            self.MaterialWidth[2*i] = self.xres*np.sum(
                self.layerWidth[indx]) - self.MaterialWidth[2*i+1]
        #  print "------debug-----", np.sum(self.MaterialWidth)
        self.netStrain = 100 * np.sum(
            self.MaterialWidth*self.eps_perp
        ) / np.sum(self.MaterialWidth)  # in percentage
        self.avghwLO = np.sum(self.MaterialWidth*self.hwLO
                              ) / np.sum(self.MaterialWidth)

        self.MLThickness = np.zeros(self.layerMaterials.size)
        for n, (MLabel, BLabel) in enumerate(zip((1, 1, 2, 2, 3, 3, 4, 4),
                                                 (0, 1) * 4)):
            # MLThickness is monolayer thickness of the material
            self.MLThickness[
                (self.layerMaterials == MLabel) &
                (self.layerBarriers == BLabel)
            ] = self.a_perp[n] / 2.0

        # Pikus-Bir interaction correction to bands offset,
        # According to Kale's, Eq.(2.14),
        # Pec for \delta E_{c} and Pe for \delta E_{v}
        self.Pec = (2*self.eps_parallel+self.eps_perp) * (self.acG)
        self.Pe = (2*self.eps_parallel+self.eps_perp) * (self.av)
        # Kale's Thesis, Eq.(2.16)
        self.Qe = - self.b * (self.c11+2*self.c12) /\
            self.c11 * self.eps_parallel

        # corrections to the method used to calculate band edges,
        # thanks to Yu Song
        # conduction band edge at different point, Eq.(2.7)
        # Varsh correction is added later
        self.EcG = self.VBO + self.EgG + self.Pec - bandBaseln
        # band edge at L and X?
        # only used in diagram..
        self.EcL = self.VBO + self.EgL + (
            2*self.eps_parallel+self.eps_perp) * (self.acL+self.av) \
            - bandBaseln
        self.EcX = self.VBO + self.EgX + (
            2*self.eps_parallel+self.eps_perp)*(self.acX+self.av) \
            + 2/3 * self.XiX * (self.eps_perp-self.eps_parallel) \
            - bandBaseln

        self.ESO = sqrt(9*self.Qe**2+2*self.Qe*self.DSO+self.DSO**2)
        self.EgLH = self.EgG + self.Pec + self.Pe \
            - 1/2*(self.Qe - self.DSO + self.ESO)
        self.EgSO = self.EgG + self.Pec + self.Pe \
            - 1/2*(self.Qe - self.DSO - self.ESO)

        # Varsh correction comes here
        # temperature correction to conduction band edge, Eq.(2.10) in Kale's
        self.Varsh = - self.alG*cst.Temperature**2/(cst.Temperature+self.beG)
        # the Varsh correction should be part conduction band, part valence
        # 1st MAJOR assumption:
        #   Varshney contribution to band edge is in proportion to percent
        #   of band offset
        # 2nd major assumption:
        #   Temperature affects sattelite valleys in the same way it does
        #   the Gamma valley
        # ??But why is it working on offset between wells and barriers, not
        # within the same material? This is not reasonable! (TODO)
        barrs = np.array([1, 3, 5, 7])
        wells = np.array([0, 2, 4, 6])
        CBOffset = self.EcG[barrs] - self.EcG[wells]
        VBOffset = (self.EcG[barrs] - self.EgLH[barrs]) \
            - (self.EcG[wells] - self.EgLH[wells])
        percentCB = CBOffset / (CBOffset + VBOffset)
        percentCB = np.column_stack([percentCB, percentCB]).flatten()
        # applies percent CV to both well and barrier slots
        self.EcG += percentCB * self.Varsh
        self.EcL += percentCB * self.Varsh
        self.EcX += percentCB * self.Varsh
        self.EvLH = self.EcG - self.EgLH - ((1-percentCB) * self.Varsh)
        self.EvSO = self.EcG - self.EgSO - ((1-percentCB) * self.Varsh)

        # Eq.(2.20) in Kale's, with Eq=0. Note that E(C-SO) = EgSO = ESO+EgLH
        self.me = 1 / ((1+2*self.F) + self.Ep/self.EgLH * (
            self.EgLH+2/3*self.ESO)/(self.EgLH + self.ESO))

    def eff_mass(self, E):
        """Calculate effective mass according to energy E,
        according to Eq.(2.20) in Kale's thesis
        """
        #  xMcE = self.xMc * (1 - (self.xVc - E) / self.xEg)
        xMcE = 1 / (1+2*self.xF + self.xEp/3 * (
            2 / ((E-self.xVc)+self.xEg) + 1 / (
                 (E-self.xVc)+self.xEg+self.xESO)))
        return xMcE

    def solve_psi(self):
        """ solve eigen modes
        OUTPUT: (doesn't return, but update member variables
            self.EigenE is the eignenergy of the layer structure
            self.xPointsPost[x] is a shorter version of self.xPoints
            self.xyPsi[x, n] is the wave function at position
                    self.xPointsPost[x] corresiponding to the
                    eigenenergy EigenE[n], and without solutions near zero
            self.xyPsiPsi[x, n] is the scaled norm of xyPsi
             -- and the above two also cut long zero heads and tials --
             -- for better plot --
            self.xyPsiPsi2[x, n] is a more precise version corresponding to
                    position self.xPoints[x]
        TODO: try matrix eigen solver?
        """
        Epoints = np.arange(min(self.xVc),
                            max(self.xVc-115*self.EField*1e-5),  # ?115e-5?
                            self.vertRes/1000)
        xMcE = np.zeros(self.xPoints.shape)
        xPsi = np.zeros(self.xPoints.shape)
        psiEnd = np.zeros(Epoints.size)

        # TODO: add adaptive spacing for Eq
        if __USE_CLIB__:
            # Call C function to get boundary dependence of energy EPoints[n],
            # the return value is psiEnd[n]
            # for n with psiEnd[n]=0, EPoints[n] is eigenenergy
            cQ.psiFnEnd(Epoints.ctypes.data_as(c_void_p),
                        int(Epoints.size), int(xPsi.size),
                        c_double(self.xres), c_double(self.EField),
                        self.xVc.ctypes.data_as(c_void_p),
                        self.xEg.ctypes.data_as(c_void_p),
                        self.xF.ctypes.data_as(c_void_p),
                        self.xEp.ctypes.data_as(c_void_p),
                        self.xESO.ctypes.data_as(c_void_p),
                        self.xMc.ctypes.data_as(c_void_p),
                        xMcE.ctypes.data_as(c_void_p),
                        xPsi.ctypes.data_as(c_void_p),
                        psiEnd.ctypes.data_as(c_void_p))
        else:
            for p, Eq in enumerate(Epoints):
                xMcE = m0 * eff_mass(Eq)
                xMcE[0:-1] = 0.5 * (xMcE[0:-1]+xMcE[1:])
                xPsi[0] = 0
                xPsi[1] = 1
                for q in xrange(1, xPsi.size-1):
                    xPsi[q+1] = xMcE[q] * ((2 * (self.xres*1e-10 / hbar)**2 *
                                           (self.xVc[q] - Eq)*e0 +
                                           1 / xMcE[q] + 1 / xMcE[q-1]) *
                                           xPsi[q] - xPsi[q-1] / xMcE[q-1])
                psiEnd[p] = xPsi[-1]
        if __LOG__:
            global logcount
            with file("EpointsLog%d.pkl" % logcount, 'w') as logfile:
                pickle.dump((Epoints, psiEnd), logfile)
            logcount += 1
            print "log saved for Epoints and psiEnd (%d)" % logcount
        # TODO: maybe improved
        tck = interpolate.splrep(Epoints, psiEnd)
        self.EigenE = interpolate.sproot(tck, mest=len(Epoints))

        if __MORE_INTERPOLATION__:
            # Near the above approximation result,
            # try to get a more precise result
            xnear = np.empty(3*len(self.EigenE))
            fxnear = np.empty(3*len(self.EigenE))
            for q in xrange(self.EigenE.size):
                # 100000 is an estimate for the precision of above
                # approximation
                # TODO: change the three calls of psiFn to loop
                approxwidth = self.vertRes/100000
                xnear[3*q:3*q+3] = self.EigenE[q] + approxwidth * \
                    np.arange(-1, 2)

            for n, Eq in enumerate(xnear):
                if __USE_CLIB__:
                    cQ.psiFn(c_double(Eq), int(1), int(xPsi.size),
                             c_double(self.xres),
                             self.xVc.ctypes.data_as(c_void_p),
                             self.xEg.ctypes.data_as(c_void_p),
                             self.xF.ctypes.data_as(c_void_p),
                             self.xEp.ctypes.data_as(c_void_p),
                             self.xESO.ctypes.data_as(c_void_p),
                             self.xMc.ctypes.data_as(c_void_p),
                             xMcE.ctypes.data_as(c_void_p),
                             xPsi.ctypes.data_as(c_void_p))
                else:
                    xMcE = self.eff_mass(Eq)
                    xMcE[1:] = (xMcE[:-1] + xMcE[1:])/2
                    xPsi[0] = 0
                    for x in range(1, xPsi.size):
                        # not tested.. C files are better
                        xPsi[x] = (xPsi[x] * (
                            2*(self.xres*ANG/hbar)**2 * (self.xVc[x] - Eq) *
                            e0 + 1/xMcE[x] + 1/xMcE[n-1]) -
                                   xPsi[x-1]/xMcE[x-1]) * xMcE[x]
                fxnear[n] = xPsi[-1]
            idxs = 3*np.arange(len(self.EigenE))+1
            if __USE_CLIB__:
                cQ.inv_quadratic_interp(xnear.ctypes.data_as(c_void_p),
                                        fxnear.ctypes.data_as(c_void_p),
                                        idxs.ctypes.data_as(POINTER(c_int)),
                                        int(idxs.size),
                                        self.EigenE.ctypes.data_as(c_void_p))
            else:
                for q, idx in enumerate(idxs):
                    # do quadratic interpolation
                    x0 = xnear[idx-1]
                    fx0 = fxnear[idx-1]
                    x1 = xnear[idx]
                    fx1 = fxnear[idx]
                    x2 = xnear[idx+1]
                    fx2 = fxnear[idx+1]
                    d1 = (fx1-fx0)/(x1-x0)
                    d2 = (fx2-fx1)/(x2-x1)
                    # inverse quadratic interpolation
                    x3 = x0*fx1*fx2/(fx0-fx1)/(fx0-fx2) +\
                        x1*fx0*fx2/(fx1-fx0)/(fx1-fx2) +\
                        x2*fx0*fx1/(fx2-fx0)/(fx2-fx1)
                    self.EigenE[q] = x3

        # make array for Psi and fill it in
        if __USE_CLIB__:
            # with eigenenregy EigenE, here call C function to get wave
            # function
            self.xyPsi = np.zeros(self.xPoints.size*self.EigenE.size)
            cQ.psiFill(int(xPsi.size), c_double(self.xres),
                       int(self.EigenE.size),
                       self.EigenE.ctypes.data_as(c_void_p),
                       self.xVc.ctypes.data_as(c_void_p),
                       self.xEg.ctypes.data_as(c_void_p),
                       self.xF.ctypes.data_as(c_void_p),
                       self.xEp.ctypes.data_as(c_void_p),
                       self.xESO.ctypes.data_as(c_void_p),
                       self.xMc.ctypes.data_as(c_void_p),
                       xMcE.ctypes.data_as(c_void_p),
                       self.xyPsi.ctypes.data_as(c_void_p))
            # self.xyPsi.resize(a.xPoints.size, a.EigenE.size)
            self.xyPsi = self.xyPsi.reshape(self.xPoints.size,
                                            self.EigenE.size, order='F')
        else:
            self.xyPsi = np.zeros((self.xPoints.size, self.EigenE.size))
            for p, Eq in enumerate(self.EigenE):
                xMcE = m0 * eff_mass(self.EigenE)
                xMcE[0:-1] = 0.5 * (xMcE[0:-1]+xMcE[1:])
                xPsi[1] = 1
                for q in xrange(2, xPsi.size-1):
                    xPsi[q+1] = ((2 * (self.xres*1e-10 / hbar)**2 *
                                  (self.xVc[q] - Eq)*e0 + 1 / xMcE[q] +
                                  1 / xMcE[q-1]) *
                                 xPsi[q] - xPsi[q-1] / xMcE[q-1]) * xMcE[q]
                psiInt = np.sum(xPsi**2 * (
                    1+(Eq-self.xVc)/(Eq-self.xVc+self.xEg)))
                A = 1 / sqrt(self.xres * 1e-10 * psiInt)
                self.xyPsi[:, p] = A * xPsi

        # remove states that come from oscillating end points
        # TODO: change to remove non-bounded states, with user options
        #       wf with non-negeligiable amplitudes higher than barrier
        #       should be removed
        if True:
            # looks like we should change -1 to -2 (following)???
            #  psiEnd = self.xyPsi[-1,:]
            #  idxs = abs(psiEnd)<10
            #  idxs = np.nonzero(abs(psiEnd)<10)[0]
            psiEnd = self.xyPsi[-2, :]
            idxs = np.abs(psiEnd) < 200/self.xres
            # 200 depends on how precise we want about eigenenergy solver
            # (TODO: more analysis and test about this value
            self.EigenE = self.EigenE[idxs]
            self.xyPsi = self.xyPsi[:, idxs]

        # 4.5e-10 scales size of wavefunctions, arbitrary for nice plots
        self.xyPsiPsi = self.xyPsi*self.xyPsi*wf_scale

        # remove states that are smaller than minimum height (remove zero
        # solutions?)-test case not showing any effect
        # addresses states high above band edge
        # 0.014 is arbitrary; changes if 4.5e-10 changes
        #  idxs = np.nonzero(self.xyPsiPsi.max(0) >
        #                   wf_scale * wf_min_height)[0]
        idxs = self.xyPsiPsi.max(0) > wf_scale*wf_min_height
        self.EigenE = self.EigenE[idxs]
        self.xyPsi = self.xyPsi[:, idxs]
        self.xyPsiPsi = self.xyPsiPsi[:, idxs]
        self.xyPsiPsi2 = copy.deepcopy(self.xyPsiPsi)

        # implement pretty plot:
        # remove long zero head and tail of the wave functions
        # test case shows on "solve whole"
        for q in xrange(self.EigenE.size):
            # 0.0005 is arbitrary
            prettyIdxs = np.nonzero(
                self.xyPsiPsi[:, q] >
                wf_scale * pretty_plot_factor)[0]
            self.xyPsiPsi[0:prettyIdxs[0], q] = np.NaN
            self.xyPsiPsi[prettyIdxs[-1]:, q] = np.NaN

        # decimate plot points: seems for better time and memory performance?
        idxs = np.arange(0, self.xPoints.size, int(
            plot_decimate_factor/self.xres), dtype=int)
        self.xyPsiPsiDec = np.zeros((idxs.size, self.EigenE.size))
        for q in xrange(self.EigenE.size):
            self.xyPsiPsiDec[:, q] = self.xyPsiPsi[idxs, q]
        self.xyPsiPsi = self.xyPsiPsiDec
        self.xPointsPost = self.xPoints[idxs]

    def basisSolve(self):
        """ solve basis for the QC device, with each basis being eigen mode of
        a seperate part of the layer structure
        OUTPUT:
            dCL: a list, each element is a QCLayers class, with layer structure
                  limited within a seperate sigle active/injection area, and
                  layer structure in dCL also includes pedding at head/tail
                  with same material as the first/last layer and barrier type
        """
        # self.basisInjectorAR is 0-to-1
        # self.basisARInjector is 1-to-0

        # find all 0-to-1 & 1-to-0 transition points
        # (1 for active region, 0 for injection retion, and active regions are
        # required to start and end by a well layer)
        # where for all n,
        # self.layerARs[zeroTOone[n]] = self.layerARs[oneTOzero[n]] = 0
        # but self.layerARs[zeroTOone[n]+1] = self.layerARs[oneTOzero[n]-1] = 1
        # TODO: try always at left
        zeroTOone = []
        oneTOzero = []
        for q in xrange(0, self.layerARs.size-1):
            if self.layerARs[q] == 0 and self.layerARs[q+1] == 1:
                zeroTOone.append(q)
            if self.layerARs[q] == 1 and self.layerARs[q+1] == 0:
                oneTOzero.append(q+1)

        dividers = [0, self.layerARs.size-1]
        if self.basisInjectorAR:
            dividers += zeroTOone
        if self.basisARInjector:
            dividers += oneTOzero
        # This converts the list into a set, thereby removing duplicates,
        # and then back into a list.
        dividers = list(set(dividers))
        dividers.sort()

        # this is dataClassesList.
        # it holds all of the Data classes for each individual solve section
        dCL = []
        # for first period only
        # this handles all of the solving
        for n in range(len(dividers)-1):
            dCL.append(copy.deepcopy(self))
            dCL[n].repeats = 1

            # substitute proper layer characteristics into dCL[n], hear/tail
            #  padding
            layer = range(dividers[n], dividers[n+1]+1)
            dCL[n].layerWidth = np.concatenate(
                    ([int(PAD_HEAD/self.xres)],
                        self.layerWidth[layer],
                        [int(PAD_TAIL/self.xres)]))
            dCL[n].layerBarriers = np.concatenate(
                    ([1], self.layerBarriers[layer], [1]))
            dCL[n].layerARs = np.concatenate(
                    ([0], self.layerARs[layer], [0]))
            dCL[n].layerMaterials = np.concatenate(
                    ([self.layerMaterials[layer][0]],
                        self.layerMaterials[layer],
                        [self.layerMaterials[layer][-1]]))
            dCL[n].layerDopings = np.concatenate(
                    ([0], self.layerDopings[layer], [0]))
            dCL[n].layerDividers = np.concatenate(
                    ([0], self.layerDividers[layer], [0]))

            # update and solve
            dCL[n].update_alloys()
            dCL[n].update_strain()
            dCL[n].populate_x()
            dCL[n].populate_x_band()
            dCL[n].solve_psi()

            # caculate offsets
            dCL[n].widthOffset = self.xres * np.sum(
                    self.layerWidth[range(0, dividers[n])])
            dCL[n].fieldOffset = -(dCL[n].widthOffset-PAD_HEAD) *\
                ANG * dCL[n].EField * KVpCM

        # create dCL's and offsets for repeat periods
        period = len(dCL)
        counter = period
        if self.repeats > 1:
            for q in xrange(1, self.repeats):
                for p in xrange(0, period):
                    dCL.append(copy.deepcopy(dCL[p]))
                    dCL[counter].widthOffset = self.xres * np.sum(
                            self.layerWidth[1:])*q + dCL[p].widthOffset
                    dCL[counter].fieldOffset = -(
                        dCL[counter].widthOffset-100) * \
                        ANG * dCL[counter].EField * KVpCM
                    counter += 1
        return dCL

    def convert_dCL_to_data(self, dCL):
        """ post processng of dCL list
        INPUT:
            self: orginal QCLayers class
            dCL: result of basisSolve(self)
        OUPUT:
            get wave functions (dCL[n].xyPsi) and eigenenrgies (dCL[n].EigenE)
            in dCL and update them in self; format them in length compatibale
            for self and update self.xyPsiPsi
            self.moduleID: moduleID[n] is the label of the position area for
                    mode self.eigenE[n] and self.xyPsi[n]
        """
        # count number of wavefunctions
        numWFs = sum([dC.EigenE.size for dC in dCL])

        self.xPointsPost = np.arange(-PAD_HEAD, self.xPoints[-1] + PAD_TAIL +
                                     self.xres, self.xres)
        self.xyPsi = np.zeros((self.xPointsPost.size, numWFs))
        self.xyPsiPsi = np.NaN*np.zeros(self.xyPsi.shape)
        self.EigenE = np.zeros(numWFs)
        self.moduleID = np.zeros(numWFs, dtype=np.int8)
        counter = 0
        for n, dC in enumerate(dCL):
            for q in xrange(dC.EigenE.size):
                self.EigenE[counter] = dC.EigenE[q] + dC.fieldOffset
                self.moduleID[counter] = n
                wf = np.zeros(self.xPointsPost.size)
                begin = int(dC.widthOffset/self.xres)
                end = begin + dC.xyPsi[:, q].size
                wf[begin:end] = dC.xyPsi[:, q]
                self.xyPsi[:, counter] = wf
                self.xyPsiPsi[:, counter] = wf**2 * wf_scale
                counter += 1

        # cut head and tial to promise the figure is in the right place?
        head = int(PAD_HEAD/self.xres)
        tail = -int(PAD_TAIL/self.xres)
        self.xPointsPost = self.xPointsPost[head:tail]
        self.xyPsi = self.xyPsi[head:tail]
        self.xyPsiPsi = self.xyPsiPsi[head:tail]

        # implement pretty plot
        # remove long zero head and tail of the wave functions
        for q in xrange(self.EigenE.size):
            prettyIdxs = np.nonzero(
                self.xyPsiPsi[:, q] >
                wf_scale * pretty_plot_factor)[0]
            if prettyIdxs.size != 0:
                self.xyPsiPsi[0:prettyIdxs[0], q] = np.NaN
                self.xyPsiPsi[prettyIdxs[-1]:, q] = np.NaN
            else:
                self.xyPsiPsi[:, q] = np.NaN

        # sort by ascending energy
        sortID = np.argsort(self.EigenE)
        self.EigenE = self.EigenE[sortID]
        self.xyPsi = self.xyPsi[:, sortID]
        self.xyPsiPsi = self.xyPsiPsi[:, sortID]
        self.moduleID = self.moduleID[sortID]

        #  #decimate plot points
        # idxs = np.arange(0,self.xPoints.size, int(
        #     plot_decimate_factor/self.xres), dtype=int)
        #  self.xyPsiPsiDec = np.zeros([idxs.size, self.EigenE.size])
        #  for q in xrange(self.EigenE.size):
        #      self.xyPsiPsiDec[:,q] = self.xyPsiPsi[idxs,q]
        #  self.xyPsiPsi = self.xyPsiPsiDec
        #  self.xPointsPost = self.xPoints[idxs]

    def lo_transition_rate(self, upper, lower):
        """ LO phonon scattering induced decay life time calculator
        INPUT:
            upper: the higher energy state index
            lower: the lower energy state index
        OUTPUT
            T1 decay life time between upper and lower states induced by LO
            phonon scattering
        """
        if upper < lower:
            upper, lower = lower, upper

        psi_i = self.xyPsi[:, upper]
        psi_j = self.xyPsi[:, lower]
        E_i = self.EigenE[upper]
        E_j = self.EigenE[lower]
        #  print "---debug---"
        #  print "E_i = %f; E_j=%f"%(E_i, E_j)

        if E_i-E_j-self.hwLO[0] < 0:
            # energy difference is smaller than a LO phonon
            # LO phonon scatering doesn't happen
            return INV_INF

        # zero head and tail cut off
        idxs_i = np.nonzero(
            psi_i > wf_scale * phonon_integral_factor)[0]
        idxs_j = np.nonzero(
            psi_j > wf_scale * phonon_integral_factor)[0]
        idx_first = (idxs_i[0], idxs_j[0])
        idx_last = (idxs_i[-1], idxs_j[-1])
        if max(idx_first) > min(idx_last):
            # wavefunction not overlap
            return INV_INF

        idx_first = min(idx_first)
        idx_last = max(idx_last)
        psi_i = psi_i[idx_first:idx_last]
        psi_j = psi_j[idx_first:idx_last]
        xPoints = self.xPoints[idx_first:idx_last]

        xMcE_j = self.eff_mass(E_j)
        # weight non-parabolic effective mass by probability density
        McE_j = m0*np.sum(xMcE_j[idx_first:idx_last] * psi_j**2) /\
            np.sum(psi_j**2)
        xMcE_i = self.eff_mass(E_i)
        # weight non-parabolic effective mass by probability density
        McE_i = m0*np.sum(xMcE_i[idx_first:idx_last] * psi_i**2) /\
            np.sum(psi_i**2)
        #  print McE_i, McE_j

        # Kale's thesis Eq.(2.68)
        kl = sqrt(2*McE_j/hbar**2 * (E_i-E_j-self.hwLO[0])*e0)
        if __USE_CLIB__:
            inv_tau_int = cQ.inv_tau_int
            inv_tau_int.restype = c_double
            Iij = inv_tau_int(xPoints.size, c_double(self.xres),
                              c_double(kl), xPoints.ctypes.data_as(c_void_p),
                              psi_i.ctypes.data_as(c_void_p),
                              psi_j.ctypes.data_as(c_void_p))
        else:
            dIij = np.empty(xPoints.size)
            for n in xrange(xPoints.size):
                x1 = xPoints[n]*ANG
                x2 = xPoints*ANG
                # first integral for eq.(2.69)
                dIij[n] = np.sum(psi_i*psi_j * exp(-kl*abs(x1-x2)) *
                                 psi_i[n]*psi_j[n] * (self.xres*ANG)**2)
            Iij = np.sum(dIij)
            #  psi_corr = self.xres*ANG * psi_i * psi_j
            #  x = np.outer(xPoints * ANG, np.ones(xPoints.size))
            #  dIij = np.outer(psi_corr, psi_corr) * exp(-kl * np.abs(x-x.T))
            #  Iij = np.sum(dIij)
        # looks similiar with eq.(2.69) but not exact in detail
        inverse_tau = sqrt(McE_j*McE_i) * e0**2 * self.hwLO[0]*e0/hbar * \
            Iij / (4 * hbar**2 * self.epsrho[0]*eps0 * kl)
        #  print "rate = %f"%(inverse_tau/1e12)
        return inverse_tau/1e12  # to ps

    def lo_life_time(self, state):
        """ return the life time due to LO phonon scattering of the
        given state(label)
        TODO: ?what if state is a lower state and there's no coupled lower
        states?"""
        rate = [self.lo_transition_rate(state, q)
                for q in range(state)]
        #  print "---debug---", rate
        rate = sum(rate)
        return 1/rate

    def dipole(self, upper, lower):
        """ Return optical dipole between self's upper level state
        and lower level state, in unit angstrom
        z = i\hbar/(2\Delta E) \<\psi_i|(m*^{-1} P_z + P_z m*^{-1})|\psi_j\>
        """
        # TODO: improve performance
        if upper < lower:
            upper, lower = lower, upper
        psi_i = self.xyPsi[:, upper]
        psi_j = self.xyPsi[:, lower]
        E_i = self.EigenE[upper]
        E_j = self.EigenE[lower]

        #  self.populate_x_band()
        # This energy dependence can be as large as -70%/+250%...
        xMcE_i = self.eff_mass(E_i)
        #  print max(xMcE_i/self.xMc), min(xMcE_i/self.xMc)
        xMcE_j = self.eff_mass(E_j)
        #  print xMcE_j/self.xMc
        xMcE_j_avg = 0.5 * (xMcE_j[0:-1]+xMcE_j[1:])
        psi_i_avg = 0.5 * (psi_i[0:-1]+psi_i[1:])
        # Kale's (2.43) and (2.47), however for varying eff mass model, this
        # should start with (2.36)
        z = np.sum(psi_i_avg * np.diff(psi_j/xMcE_i) +
                   1/xMcE_j_avg * (psi_i_avg * np.diff(psi_j)))
        z *= hbar**2/(2*(E_i-E_j)*e0*m0) / ANG  # e0 transform eV to J
        return z

    def coupling_energy(self, dCL, upper, lower):
        """Calculate the coupling energy between upper level and lower level
        with levels(basis) defined in dCL
        coupling energy = <upper|H|lower> with H = H0 + V1 + V2,
        H0 + V1 |upper> = E(upper) |upper>; H0 + V2 |lower> = E(lower) |lower>;
        so <upper|H|lower> = <upper|V2|lower> + E(upper) <upper|lower>
                           = <upper|V1|lower> + E(lower) <upper|lower>
        while H0 includes potential without wells, V1 and V2 are wells for
        module/dCL corresponds to upper and lower respectively
        The result is only used in the calculate box..
        Old version only includes <upper|V1|lower>
        """
        # here, psi_i is the left-most wavefunction,
        # not the wf with the highest energy
        # but does it matter?..
        module_i = self.moduleID[upper]
        module_j = self.moduleID[lower]
        if module_i > module_j:
            module_i, module_j = module_j, module_i
            psi_i = self.xyPsi[:, lower]
            psi_j = self.xyPsi[:, upper]
            Ej = self.EigenE[upper]
        else:
            psi_i = self.xyPsi[:, upper]
            psi_j = self.xyPsi[:, lower]
            Ej = self.EigenE[lower]

        # old version of coupling calculation
        #  DeltaV = np.ones(self.xPointsPost.size)
        #  first = int(dCL[module_i].widthOffset/self.xres)
        # last = first + dCL[module_i].xBarriers[
        #     int(PAD_HEAD/self.xres):].size
        #  print "---debug--- coupling_energy"
        #  print first,last
        # DeltaV[first:last] = dCL[module_i].xBarriers[
        #     int(PAD_HEAD/self.xres):]
        #  DeltaV = 1 - DeltaV #=is well
        #  couplingEnergy = np.sum(psi_i * DeltaV * psi_j) * self.xres * ANG \
        #        * abs(self.EcG[1] - self.EcG[0]) /meV #* (1-self.xBarriers)
        # DeltaV * (self.EcG[1] (barrier) - dta.Ecg[0] (well)) = Vi(wells)

        # Ming's version for calculating coupling, 08.23.2017
        if module_j - module_i != 1:
            return 0
        DeltaV = np.ones(self.xPointsPost.size)
        first = int(dCL[module_i].widthOffset/self.xres)
        last = first + dCL[module_i].xBarriers[
                int(PAD_HEAD/self.xres):int(PAD_TAIL/self.xres)].size
        #  print "---debug--- coupling_energy"
        #  print first,last
        DeltaV[first:last] = dCL[module_i].xBarriers[
                int(PAD_HEAD/self.xres):int(PAD_TAIL/self.xres)]
        DeltaV = 1 - DeltaV  # =is well
        jMat = int(self.xMaterials[last+1])
        DeltaV *= (self.EcG[2*jMat-1] - self.EcG[2*(jMat-1)])/meV  # unit meV
        couplingEnergy = (np.sum(psi_i * (DeltaV + Ej) * psi_j)) * \
            self.xres * ANG
        return couplingEnergy  # unit meV

    def broadening_energy(self, upper, lower):
        """interface roughness induced broadening: Khurgin, yentings thesis"""
        if upper < lower:
            upper, lower = lower, upper
        psi_i = self.xyPsi[:, upper]
        psi_j = self.xyPsi[:, lower]

        transitions = np.bitwise_xor(self.xBarriers[0:-1].astype(bool),
                                     self.xBarriers[1:].astype(bool))
        transitions = np.append(transitions, False)  # last element is not
        psi2int2 = np.sum((psi_i[transitions]**2-psi_j[transitions]**2)**2)
        DeltaLambda = 0.76 * 1e-9 * 1e-9  # 0.79nm^2
        # effective mass (self.me) update?
        twogamma = pi*self.me[0]*m0*e0**2/hbar**2 * DeltaLambda**2 *\
            (self.EcG[1] - self.EcG[0])**2 * psi2int2
        twogamma /= meV*e0  # convert to meV
        return twogamma

    def alphaISB(self, stateR, lower):
        """intersubband transition.. etc."""
        statesQ = []
        dipoles = []
        gammas = []
        energies = []
        for q in xrange(stateR+1, self.EigenE.size):
            dp = self.dipole(q, stateR)
            if abs(dp) > 1e-6:
                statesQ.append(q)
                dipoles.append(dp)
        for q in statesQ:
            gamma = self.broadening_energy(q, stateR)/2
            gammas.append(gamma)
            energies.append(self.EigenE[q]-self.EigenE[stateR])

        dipoles = np.array(dipoles)*1e-10   # in m
        gammas = np.array(gammas)/1000      # from meV to eV
        energies = abs(np.array(energies))  # in eV

        neff = 3
        Lp = self.xres * np.sum(self.layerWidth[1:]) * 1e-10  # in m
        Nq = np.sum(self.layerDopings[1:]*self.layerWidth[1:]) /\
            np.sum(self.layerWidth[1:])
        Nq *= 100**3  # convert from cm^-3 to m^-3
        Ns = self.xres * np.sum(self.layerDopings[1:] *
                                self.layerWidth[1:]) * 1e11  # in cm^-2
        Ns *= 100**2  # from cm^-2 to m^-2
        hw = self.EigenE[stateR] - self.EigenE[lower]

        # hw = np.arange(0.15, 0.5, 0.01)
        # for enerG in hw:
        #     alphaISB = np.sum(energies * dipoles**2 * gammas / (
        #         #(energies - enerG)**2 + gammas**2))

        alphaISB = np.sum(energies*e0/h/c0 * dipoles**2 * gammas /
                          ((energies - hw)**2 + gammas**2))
        alphaISB *= 4*pi*e0**2 / (eps0*neff) * pi/(2*Lp) * Ns
        alphaISB /= e0*100

        if False:  # plot loss diagram
            hw = np.arange(0.15, 0.5, 0.001)
            alphaISBw = np.zeros(hw.size)
            for q, enerG in enumerate(hw):
                alphaISBw[q] = np.sum(energies*e0/h/c0 * dipoles**2 * gammas /
                                      ((energies - enerG)**2 + gammas**2))
            alphaISBw *= 4*pi*e0**2 / (eps0*neff) * pi/(2*Lp) * Ns
            alphaISBw /= e0*100
            plt.figure()
            plt.plot(hw, alphaISBw)
            plt.show()

        return alphaISB

    def figure_of_merit(self, upper, lower):
        if upper < lower:
            upper, lower = lower, upper
        tauLower = self.lo_life_time(lower)
        tauUpper = self.lo_life_time(upper)
        tauUpperLower = 1/self.lo_transition_rate(upper, lower)
        opticalDipole = self.dipole(upper, lower)
        return opticalDipole**2 * tauUpper * (1 - tauLower/tauUpperLower)


if __name__ == "__main__":
    if __USE_CLIB__:
        print "1/alpha ~ %d" % cQ.inv_alpha()
        print "It's one of the greatest damn mysteries of physics"
    else:
        print "QCLayers.py called"

# vim: ts=4 sw=4 sts=4 expandtab
