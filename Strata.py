#!/usr/bin/env python
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

from __future__ import division

__USE_CLIB__ = True

import numpy as np
from numpy import sqrt, exp, sin, cos, log, pi, conj, real, imag
from scipy import interpolate
import copy
import sys

from QCLayers import cst

# TODO: replace CLIB by Cython
from ctypes import *
try:
    if sys.platform in ('linux2', 'darwin', 'cygwin'):
        cS = CDLL('./cStrata.so')
    elif sys.platform == 'win32':
        cS = CDLL('cStrata.dll')
except:
    print "unable to load cStrata"

# =============================================================================
# Global Variables
# =============================================================================
from scipy.constants import h
from scipy.constants import (e as e0, epsilon_0 as eps, c as c0)

# ===========================================================================
# Reference
# [0]Kale Franz's thesis
# [1]Handbook of Optics, Vol.2, ISBN: 0070479747
# ===========================================================================


def zero_find(xVals, yVals):
    """To find zero points for function y(x) using iterpolation
    """
    # TODO: may be improved for near degenerate states
    # For eigen energy solver, psiEnd's dependence on energy is significant
    # near eigenenergy
    tck = interpolate.splrep(xVals.real, yVals.real)
    #  print "------debug------ Here zero_find is called"
    return interpolate.sproot(tck, mest=len(xVals))


class Strata(object):
    """Strata property for optical mode solver
    """
    def __init__(self):
        self.stratumMaterials = ['InP']
        self.stratumCompositions = np.array([0.])
        self.stratumThicknesses = np.array([0.])
        self.stratumDopings = np.array([0.])

        self.wavelength = 4.7  # unit micron
        self.operatingField = 0
        self.Lp = 1
        self.Np = 1

        # aCore is the alpha defined in [1], Sec 36.3, (10), measured in cm-1
        # representing decay rate in the material.. wl independent?
        self.aCore = 0  # = 4\pi k /\lambda
        self.nCore = 4  # index of the active core?
        self.nD = 0     # ?not used

        self.tauUpper = 0.0
        self.tauLower = 0.0
        self.tauUpperLower = 1.0e-3
        self.opticalDipole = 0.0
        self.FoM = 0.0
        self.transitionBroadening = 1.0e-5  # eV
        self.waveguideFacets = 'as-cleaved + as-cleaved'
        # can be combination of "as-cleaved", "perfect HR", "perfect AR",
        # "custom coating"
        self.customFacet = 0.0
        self.waveguideLength = 3.0  # unit?

        self.frontFacet = 0
        self.backFacet = 0

        self.beta = 3+0j

        self.xres = 0.01  # um -> angstrom?
        self.stratumSelected = 0

        self.notDopableList = ['Air', 'Au', 'SiO2', 'SiNx']
        self.needsCompositionList = ['InGaAs', 'InAlAs']

        self.populate_rIndexes()

    def populate_rIndexes(self):
        """ Matrial reflection index for GaAs, InAs, AlAs and InP """
        wl = self.wavelength  # unit um, see [1] Table22

        self.stratumRIndexes = np.zeros(self.stratumDopings.size,
                                        dtype=np.complex128)
        for q, material in enumerate(self.stratumMaterials):
            # Calculate reflection index (complex for decay) for each stratum?
            # TODO: combine codes for different materials
            if material == 'Active Core':
                self.stratumRIndexes[q] = self.nCore
            elif material == 'InP':
                # 0.95 < wl < 10
                n_InP = cst['InP'].rIndx(wl)
                nue = 1
                me0 = cst['InP'].me0
                a = 8.97E-5*wl**2/me0*self.stratumDopings[q]
                eps = n_InP**2 - a / (1+1j*5.305e-3*wl**2*nue)
                n_InPd = sqrt(0.5 * (abs(eps) + eps.real))
                k_InPd = sqrt(0.5 * (abs(eps) - eps.real))
                self.stratumRIndexes[q] = n_InPd + 1j*k_InPd
            elif material == 'GaAs':
                # 1.4 < wl < 11
                n_GaAs = cst['GaAs'].rIndx(wl)
                nue = 1
                me0 = cst['GaAs'].me0
                a = 8.97E-5*wl**2/me0*self.stratumDopings[q]
                eps = n_GaAs**2 - a / (1+1j*5.305e-3*wl**2*nue)
                n_GaAsd = sqrt(0.5 * (abs(eps) + eps.real))
                k_GaAsd = sqrt(0.5 * (abs(eps) - eps.real))
                self.stratumRIndexes[q] = n_GaAsd + 1j*k_GaAsd
            elif material == 'InGaAs':
                # 3.7 < wl < 31.3
                n_InAs = cst['InAs'].rIndx(wl)
                # 1.4 < wl < 11
                n_GaAs = cst['GaAs'].rIndx(wl)
                xFrac = self.stratumCompositions[q]
                # bowing parameters not found: negeleted
                n_InGaAs = xFrac*n_InAs + (1-xFrac)*n_GaAs
                nue = 1
                me0 = xFrac*cst['InAs'].me0 + (1-xFrac)*cst['GaAs'].me0
                a = 8.97E-5*wl**2/me0*self.stratumDopings[q]
                eps = n_InGaAs**2 - a / (1+1j*5.305e-3*wl**2*nue)
                n_InGaAs = sqrt(0.5 * (abs(eps) + eps.real))
                k_InGaAs = sqrt(0.5 * (abs(eps) - eps.real))
                self.stratumRIndexes[q] = n_InGaAs + 1j*k_InGaAs
            elif material == 'InAlAs':
                xFrac = self.stratumCompositions[q]
                # 3.7 < wl < 31.3
                n_InAs = cst['InAs'].rIndx(wl)
                # 0.56 < wl < 2.2
                n_AlAs = cst['AlAs'].rIndx(wl)
                # bowing parameters not found: negeleted
                n_AlInAs = (1-xFrac)*n_AlAs + xFrac*n_InAs
                nue = 1
                me0 = (1-xFrac)*cst['AlAs'].me0 + xFrac*cst['InAs'].me0
                a = 8.97E-5*wl**2/me0*self.stratumDopings[q]
                eps = n_AlInAs**2 - a / (1+1j*5.305e-3*wl**2*nue)
                n_AlInAs = sqrt(0.5 * (abs(eps) + eps.real))
                k_AlInAs = sqrt(0.5 * (abs(eps) - eps.real))
                self.stratumRIndexes[q] = n_AlInAs + 1j*k_AlInAs
            elif material == 'Au':
                C1 = -0.1933
                C2 = 0.3321
                C3 = 0.0938
                D1 = -0.382
                D2 = 6.8522
                D3 = -0.1289
                n_Au = C1+wl*C2+wl*C3**2
                k_Au = D1+wl*D2+wl*D3**2
                self.stratumRIndexes[q] = n_Au+k_Au*1j
            elif material == 'SiNx':
                # from Jean Nguyen's Thesis
                C1 = 2.0019336
                C2 = 0.15265213
                C3 = 4.0495557
                D0 = -0.00282
                D1 = 0.003029
                D2 = -0.0006982
                D3 = -0.0002839
                D4 = 0.0001816
                D5 = -3.948e-005
                D6 = 4.276e-006
                D7 = -2.314e-007
                D8 = 4.982e-009
                n_SiNx = C1 + C2/wl**2 + C3/wl**4
                k_SiNx = D0 + D1*wl + D2*wl**2 + D3*wl**3 + D4*wl**4 \
                    + D5*wl**5 + D6*wl**6 + D7*wl**7 + D8*wl**8
                k_SiNx *= 100
                self.stratumRIndexes[q] = n_SiNx+k_SiNx*1j
            elif material == 'SiO2':
                # from Jean Nguyen's Thesis
                C1 = 1.41870
                C2 = 0.12886725
                C3 = 2.7573641e-5
                n_SiO2 = C1 + C2/wl**2 + C3/wl**4
                # this is a 4 peak Lorentzian fit to her data
                y0=-797.4627
                xc1=2.83043; w1=6.083822; A1=10881.9438
                xc2=8.95338; w2=1.38389113; A2=9167.662815
                xc3=12.3845492; w3=3.9792077; A3=12642.72911
                xc4=15.6387213; w4=0.6057751177; A4=3292.325272
                alpha = y0 + 2*A1/pi*w1/(4*(wl-xc1)**2+w1**2) \
                    + 2*A2/pi*w2/(4*(wl-xc2)**2+w2**2) \
                    + 2*A3/pi*w3/(4*(wl-xc3)**2+w3**2) \
                    + 2*A4/pi*w4/(4*(wl-xc4)**2+w4**2)
                k_SiO2 = alpha * wl*1e-4 / (4*pi)
                self.stratumRIndexes[q] = n_SiO2 + k_SiO2*1j
            elif material == 'Air':
                self.stratumRIndexes[q] = 1

    def populate_x(self):
        """Extend layer information to position functions?
        Layer data: stratumThicknesses
                    stratumThickNum
                    stratumMaterials
                    stratumRIndexes
                with len = # of layers and each value repr. a layer
        Position data: xPoints - position
                       xn
                       xAC
                       xStratumSelected"""
        # use rounding to work with selected resolution
        self.stratumThickNum = np.round(self.stratumThicknesses /
                                        self.xres).astype(np.int64)
        self.stratumThicknesses = self.stratumThickNum * self.xres

        # convert to int to prevent machine rounding errors
        self.xPoints = self.xres * np.arange(0, self.stratumThickNum.sum(), 1)

        stratumThickNumCumSum = np.concatenate(
            ([0], self.stratumThickNum.cumsum()))
        self.xn = np.zeros(self.xPoints.size, dtype=np.complex128)
        # binary designation for Active Core
        self.xAC = np.zeros(self.xPoints.size)

        # extend layer data for all xpoints
        for q in xrange(0, self.stratumThickNum.size):
            self.xn[stratumThickNumCumSum[q]:
                    stratumThickNumCumSum[q+1]] = self.stratumRIndexes[q]
            if self.stratumMaterials[q] == 'Active Core':
                self.xAC[stratumThickNumCumSum[q]:
                         stratumThickNumCumSum[q+1]] = 1

        # make array to show selected stratum in mainCanvas
        # ? what is this used for?
        self.xStratumSelected = np.zeros(self.xPoints.shape) * np.NaN
        if self.stratumSelected != -1:  # if not no row selected
            self.xStratumSelected[
                stratumThickNumCumSum[self.stratumSelected]:
                stratumThickNumCumSum[self.stratumSelected+1]] \
                   = self.xn.real[
                       stratumThickNumCumSum[self.stratumSelected]:
                       stratumThickNumCumSum[self.stratumSelected+1]]

    def chi_find(self, beta):
        # ?... beta is a float number
        # This function is not actually called
        print "----debug---- chi_find("+beta
        z0 = 0.003768
        k = 2*pi/self.wavelength

        alpha = sqrt(self.stratumRIndexes**2-beta**2)
        if alpha[0].imag < 0:
            alpha[0] = conj(alpha[0])
            # alpha[0] = -alpha[0]
        if alpha[-1].imag < 0:
            alpha[-1] = conj(alpha[-1])
            # alpha[-1] = -alpha[-1]
        gamma = z0*alpha/self.stratumRIndexes**2
        phi = k*self.stratumThicknesses*alpha
        # zeta  = k*self.stratumThicknesses/z0

        Mj = []
        M = np.array([[1+0j, 0], [0, 1+0j]])
        for q in xrange(alpha.size):
            Mj.append(np.array([[cos(phi[q]), -1j/gamma[q]*sin(phi[q])],
                                [-1j*gamma[q]*sin(phi[q]), cos(phi[q])]]))

        Mj.reverse()
        for mj in Mj:
            M = np.dot(mj, M)

        gammas = gamma[0]
        gammac = gamma[-1]

        chi = gammac*M[0, 0] + gammac*gammas*M[0, 1] + M[1, 0] + gammas*M[1, 1]
        return chi

    def beta_find(self, betaInit=None):
        # ? Seems to relate to EMF mode
        if True:  # betaInit == None:
            betaMax = max(self.stratumRIndexes.real)
            betaMin = min(self.stratumRIndexes.real)

            betas = np.arange(betaMin.real+0.01, betaMax.real, 0.01)

            if __USE_CLIB__:  # do chi_find in c
                chiImag = np.zeros(len(betas), dtype=float)
                betasReal = betas.real
                betasImag = betas.imag
                stratumRIndexesReal = self.stratumRIndexes.real.copy()
                stratumRIndexesImag = self.stratumRIndexes.imag.copy()
                cS.chiImag_array(
                    c_double(self.wavelength),
                    self.stratumThicknesses.ctypes.data_as(c_void_p),
                    stratumRIndexesReal.ctypes.data_as(c_void_p),
                    stratumRIndexesImag.ctypes.data_as(c_void_p),
                    int(self.stratumRIndexes.size),
                    betasReal.ctypes.data_as(c_void_p),
                    betasImag.ctypes.data_as(c_void_p), int(betasReal.size),
                    chiImag.ctypes.data_as(c_void_p))
                beta0s = zero_find(betas.real, chiImag)
            else:
                chi = np.zeros(betas.size, dtype=np.complex128)
                for p, beta in enumerate(betas):
                    chi[p] = self.chi_find(beta)
                beta0s = zero_find(betas.real, chi.imag)
            beta = max(beta0s)+1j*min(self.stratumRIndexes.imag)
        else:
            beta = betaInit

        beta_find_precision = 1e-5
        if True:  # setting to True makes the function stall in Mac OS X
            betaIn = beta
            stratumRIndexesReal = self.stratumRIndexes.real.copy()
            stratumRIndexesImag = self.stratumRIndexes.imag.copy()
            betaOut = np.array([0.0, 0.0])
            beta = cS.beta_find(
                c_double(self.wavelength),
                self.stratumThicknesses.ctypes.data_as(c_void_p),
                stratumRIndexesReal.ctypes.data_as(c_void_p),
                stratumRIndexesImag.ctypes.data_as(c_void_p),
                int(self.stratumRIndexes.size), c_double(betaIn.real),
                c_double(betaIn.imag), c_double(beta_find_precision),
                betaOut.ctypes.data_as(c_void_p))
            beta = betaOut[0] + 1j*betaOut[1]
        else:
            rInc = 0.0001
            iInc = 1j * 1e-6
            abschiNew = 1
            while True:
                betas = [beta, beta+rInc, beta-rInc, beta+iInc, beta-iInc,
                         beta+rInc+iInc, beta-rInc-iInc, beta+rInc-iInc,
                         beta-rInc+iInc]
                if True:
                    chi = np.zeros(len(betas), dtype=np.complex128)
                    for p, betaIn in enumerate(betas):
                        chi[p] = self.chi_find(betaIn)
                else:  # do chi_find in c
                    chi = np.zeros(len(betas), dtype=float)
                    abschi_find = cS.abschi_find
                    abschi_find.restype = c_double
                    for p, betaIn in enumerate(betas):
                        stratumRIndexesReal = self.stratumRIndexes.real.copy()
                        stratumRIndexesImag = self.stratumRIndexes.imag.copy()
                        chi[p] = abschi_find(
                            c_double(self.wavelength),
                            self.stratumThicknesses.ctypes.data_as(c_void_p),
                            stratumRIndexesReal.ctypes.data_as(c_void_p),
                            stratumRIndexesImag.ctypes.data_as(c_void_p),
                            int(self.stratumRIndexes.size),
                            c_double(betaIn.real), c_double(betaIn.imag))
                abschiOld = abschiNew
                abschiNew = min(abs(chi))
                idx = argmin(abs(chi))
                beta = betas[idx]
                if abs(abschiOld - abschiNew)/abschiOld < beta_find_precision:
                    break
        return beta

    def mode_plot(self):
        #  print "-----debug----- mode_plot is called"
        n = copy.copy(self.stratumRIndexes)[::-1]
        thicknesses = copy.copy(self.stratumThicknesses)[::-1]
        ThickNum = copy.copy(self.stratumThickNum)[::-1]
        #  xres = self.xres

        z0 = 0.003768
        # z0 = 376.8
        k = 2*pi/self.wavelength
        M = np.array([[1+0j, 0], [0, 1+0j]])

        alpha = sqrt(n**2-self.beta**2)
        if alpha[0].imag < 0:
            alpha[0] = conj(alpha[0])
        if alpha[-1].imag < 0:
            alpha[-1] = conj(alpha[-1])
        gamma = z0*alpha/n**2
        phi = k*thicknesses*alpha
        # zeta  = k*thicknesses/z0

        ncs = stratumThickNumCumSum = np.concatenate(([0], ThickNum.cumsum()))
        xI = np.zeros(self.xPoints.size, dtype=np.complex128)

        for q in xrange(ThickNum.size-1, -1, -1):
            xvec = copy.copy(self.xPoints[ncs[q]:ncs[q+1]])[::-1]
            if len(xvec) == 0:  # make sure xvec isn't empty
                continue
            xvec -= min(xvec)
            field = np.dot(M, np.array([1, gamma[-1]]))
            U = field[0]
            V = field[1]
            if q == 0 or q == self.stratumThicknesses.size-1:
                xI[ncs[q]:ncs[q+1]] = real(U)*exp(1j*k*alpha[q]*xvec) / n[q]**2
            else:
                xI[ncs[q]:ncs[q+1]] = U*cos(-k*alpha[q]*xvec) \
                        + 1j/gamma[q] * V*sin(-k*alpha[q]*xvec)
                xI[ncs[q]:ncs[q+1]] /= n[q]**2
                Mj = np.array([[cos(phi[q]), -1j/gamma[q]*sin(phi[q])],
                               [-1j*gamma[q]*sin(phi[q]), cos(phi[q])]])
                M = np.dot(Mj, M)

        xI = abs(xI)**2
        xI = xI / max(xI)
        self.xI = xI[::-1]

        # calculate confinement factor
        self.confinementFactor = np.sum(self.xI * self.xAC) / np.sum(self.xI)

    def calculate_performance_parameters(self):
        # waveguide loss
        self.waveguideLoss = 4 * pi * self.beta.imag \
                / (self.wavelength * 1e-6) * 1e-2

        # mirror loss
        self.mirrorLoss = -1 / (2 * self.waveguideLength * 0.1) \
            * log(self.frontFacet * self.backFacet)

        # transition cross-section
        Eph = h * c0 / (self.wavelength * 1e-6)
        neff = self.beta.real
        z = self.opticalDipole * 1e-10
        deltaE = 0.1*Eph
        sigma0 = 4*pi*e0**2 / (h*c0*eps0*neff) * Eph/deltaE * z**2

        # gain
        tauEff = self.tauUpper * (1 - self.tauLower /
                                  self.tauUpperLower) * 1e-12
        Lp = self.Lp * 1e-10
        self.gain = sigma0 * tauEff / (e0 * Lp)  # m/A
        self.gain *= 100  # cm/A

        # threshold current density
        self.Jth0 = (self.waveguideLoss + self.mirrorLoss) \
            / (self.gain * self.confinementFactor)  # A/cm^2
        self.Jth0 *= 1e-3  # kA/cm^2

        # threshold current
        self.Ith0 = self.Jth0*1e3 * (self.Np * self.Lp*1e-8) \
            * self.waveguideLength*1e-1

        # operating voltage
        self.operatingVoltage = self.operatingField*1e3 * self.Lp*1e-8 \
            * self.Np

        # voltage efficiency
        self.voltageEfficiency = 1.24/self.wavelength * self.Np \
            / self.operatingVoltage

        # extraction efficiency
        self.extractionEfficiency = self.mirrorLoss \
            / (self.mirrorLoss + self.waveguideLoss)

        # population inverstion efficiency
        tauEff = self.tauUpper * (1 - self.tauLower / self.tauUpperLower)
        self.inversionEfficiency = tauEff / (tauEff + self.tauLower)

        # modal efficiency
        xI = self.xI/max(self.xI)
        U = xI[np.nonzero(self.xAC)[0]]

        # interoplate over U_AC for each Np at the point xbar for Ubar
        # this is Faist's version
        numACs = self.stratumMaterials.count('Active Core')
        try:
            xVals = np.arange(self.xres,
                              numACs*self.Np*self.Lp*1e-4,
                              self.xres)
            assert xVals.size == U.size
        except AssertionError:
            try:
                xVals = np.arange(0, numACs*self.Np*self.Lp*1e-4, self.xres)
                assert xVals.size == U.size
            except AssertionError:
                xVals = np.arange(self.xres,
                                  numACs*self.Np*self.Lp*1e-4-self.xres,
                                  self.xres)
                assert xVals.size == U.size
        tck = interpolate.splrep(xVals, U, s=0)
        minx = 0.5*self.Lp*1e-4
        maxx = numACs*self.Np*self.Lp*1e-4-0.5*self.Lp*1e-4
        xbar = np.linspace(minx, maxx, numACs*self.Np)
        Ubar = interpolate.splev(xbar, tck, der=0)
        self.modalEfficiency = np.sum(Ubar)**2 / (
            numACs * self.Np * np.sum(Ubar**2))

        # Kale's version
        modalEfficiency = np.sum(Ubar) / self.Np
        # since Ubar taken from normalized xI
        # I guess we'll go with Faist's version since he probably knows
        # better than I do. 

    def updateFacets(self):
        if self.waveguideFacets == 'as-cleaved + as-cleaved':
            self.frontFacet = reflectivity(self.beta)
            self.backFacet = reflectivity(self.beta)
        elif self.waveguideFacets == 'as-cleaved + perfect HR':
            self.frontFacet = ThePhysics.reflectivity(self.beta)
            self.backFacet = 1
        elif self.waveguideFacets == 'as-cleaved + perfect AR':
            self.frontFacet = 1e-9
            self.backFacet = ThePhysics.reflectivity(self.beta)
        elif self.waveguideFacets == 'perfect AR + perfect HR':
            self.frontFacet = 1e-9
            self.backFacet = 1
        elif self.waveguideFacets == 'custom coating + as-cleaved':
            self.frontFacet = self.customFacet
            self.backFacet = ThePhysics.reflectivity(self.beta)
        elif self.waveguideFacets == 'custom coating + perfect HR':
            self.frontFacet = self.customFacet
            self.backFacet = 1
        elif self.waveguideFacets == 'custom coating + perfect AR':
            self.frontFacet = 1e-9
            self.backFacet = self.customFacet


def reflectivity(beta):
    # should be member method for strata
    beta = beta.real
    return ((beta - 1) / (beta + 1))**2


if __name__ == "__main__":
    print ('Answer to the Ultimate Question of Life,'
           'The Universe, and Everything is'), cS.answer()

# vim: ts=4 sw=4 sts=4 expandtab
