#===============================================================================
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
#===============================================================================

# TODO: 
# material related codes should be moved to MaterialConstants
# try separate this file to smaller ones

from __future__ import division
#  from numpy import *
import numpy as np
from numpy import sqrt, exp, sin, cos, log, pi, conj, real, imag
from scipy import interpolate
import copy
from multiprocessing import Process, Queue

#import pylab as plt

import settings

#  import MaterialConstants
#  create global variable c for material constants
#  c = MaterialConstants.MaterialConstants()

import MaterialConstantsDict
cst = MaterialConstantsDict.MaterialConstantsDict()

from ctypes import *
try:
    cFunctions=CDLL('./cFunctions.so')
except WindowsError:
    cFunctions=CDLL('cFunctions.dll')

USE_CLIB = True

#===============================================================================
# Global Variables
#===============================================================================
e0 = 1.60217653e-19  #electron charge
eps0 = 8.854187e-12
m0 = 9.10938188e-31   #free electron mass (kg)
h = 6.6260693e-34
hbar = 6.6260693e-34/(2*pi) #Planck's constant (J s)
kb = 1.386505e-23 / e0 #eV/K
c0 = 299792458
ANG = 1e-10 # angstrom to meter
KVpCM = 1e5 # KV/cm to V/m

#===============================================================================
# Reference
# [1]Handbook of Optics, Vol.2, ISBN: 0070479747
#===============================================================================

def zero_find(xVals, yVals):
    """To find zero points for function y(x) using iterpolation
    """
    tck = interpolate.splrep(xVals.real,yVals.real)
    #  print "------debug------ Here zero_find is called"
    return interpolate.sproot(tck, mest=len(xVals))
    #  xnew = np.linspace(min(xVals.real),max(xVals.real),5e5)
    #  ynew = interpolate.splev(xnew,tck,der=0)

    #  #This routine looks for all of the zero crossings, and then picks each one out
    #  gtz = ynew > 0
    #  ltz = ynew < 0
    #  overlap1 = np.bitwise_and(gtz[0:-1],ltz[1:])
    #  overlap2 = np.bitwise_and(gtz[1:],ltz[0:-1])
    #  overlap  = np.bitwise_or(overlap1, overlap2)
    #  idxs = np.nonzero(overlap == True)[0]
    #  zeroCrossings = np.zeros(idxs.size)

    #  if False:
        #  cFunctions.inv_quadratic_interp(xnew.ctypes.data_as(c_void_p), 
                #  ynew.ctypes.data_as(c_void_p),
                #  idxs.ctypes.data_as(c_void_p), int(zeroCrossings.size),
                #  zeroCrossings.ctypes.data_as(c_void_p))
    #  else:
        #  for q, idx in enumerate(idxs): # do quadratic interpolation
            #  x0=xnew[idx-1]; fx0=ynew[idx-1]
            #  x1=xnew[idx];   fx1=ynew[idx]
            #  x2=xnew[idx+1]; fx2=ynew[idx+1]
            #  d1=(fx1-fx0)/(x1-x0); d2=(fx2-fx1)/(x2-x1)
            #  #inverse quadratic interpolation
            #  x3 = x0*fx1*fx2/(fx0-fx1)/(fx0-fx2) \
                    #  + x1*fx0*fx2/(fx1-fx0)/(fx1-fx2) \
                    #  + x2*fx0*fx1/(fx2-fx0)/(fx2-fx1)
            #  zeroCrossings[q] = x3

    #  print zeroCrossings
    #  return zeroCrossings
    

class Strata(object):
    """Strata property for optical mode solver
    """
    def __init__(self):
        self.stratumMaterials = ['InP']
        self.stratumCompositions = np.array([0.])
        self.stratumThicknesses = np.array([0.])
        self.stratumDopings = np.array([0.])
        
        self.wavelength = 4.7 # unit? micron?
        self.operatingField = 0
        self.Lp = 1
        self.Np = 1

        # aCore is the alpha defined in [1], Sec 36.3, (10), measured in cm-1
        # representing decay rate in the material.. wl independent?
        self.aCore = 0 # = 4\pi k /\lambda
        self.nCore = 4 # index of the active core?
        self.nD = 0    # ?not used
        
        self.tauUpper = 0.0
        self.tauLower = 0.0
        self.tauUpperLower = 1.0e-3
        self.opticalDipole = 0.0
        self.FoM = 0.0 # ?
        self.transitionBroadening = 1.0e-5 # eV
        self.waveguideFacets = 'as-cleaved + as-cleaved'  
        # can be combination of "as-cleaved", "perfect HR", "perfect AR",
        # "custom coating"
        self.customFacet = 0.0  
        self.waveguideLength = 3.0 #unit?
        
        self.frontFacet = 0
        self.backFacet = 0
        
        self.beta = 3+0j
        
        self.xres = 0.01 #um -> angstrom?
        self.stratumSelected = 0
        
        self.notDopableList = ['Air', 'Au', 'SiO2', 'SiNx']
        self.needsCompositionList = ['InGaAs','InAlAs']

        self.populate_rIndexes()
        
    def populate_rIndexes(self):
        """ Matrial reflection index for GaAs, InAs, AlAs and InP """
        # Repeated codes, should be moved outside 1
        # n = sqrt(c1 + c2 * wl**2/(wl**2-c3**2) + c4 * wl**2/(wl**2-c5**2) )
        wl = self.wavelength
        n_GaAs = sqrt( cst['GaAs'].C1 
                + cst['GaAs'].C2*wl**2/(wl**2-cst['GaAs'].C3**2) 
                + cst['GaAs'].C4*wl**2/(wl**2-cst['GaAs'].C5**2) )
        n_InAs = sqrt( cst['InAs'].C1 
                + cst['InAs'].C2*wl**2/(wl**2-cst['InAs'].C3**2) 
                + cst['InAs'].C4*wl**2/(wl**2-cst['InAs'].C5**2) )
        n_AlAs = sqrt( cst['AlAs'].C1 
                + cst['AlAs'].C2*wl**2/(wl**2-cst['AlAs'].C3**2) 
                + cst['AlAs'].C4*wl**2/(wl**2-cst['AlAs'].C5**2) )
        n_InP  = sqrt( cst['InP'].C1 
                + cst['InP'].C2*wl**2/(wl**2-cst['InP'].C3**2) 
                + cst['InP'].C4*wl**2/(wl**2-cst['InP'].C5**2) )
        
        self.stratumRIndexes = np.zeros(self.stratumDopings.size, dtype=complex)
        for q, material in enumerate(self.stratumMaterials):
            # Calculate reflection index (complex for decay) for each stratum?
            #TODO: combine codes for different materials
            if material == 'Active Core':
                self.stratumRIndexes[q] = self.nCore
            elif material == 'InP':
                nue = 1
                me0 = cst['InP'].me0
                a = 8.97E-5*wl**2/me0*self.stratumDopings[q]
                eps = n_InP**2 - a / (1+1j*5.305e-3*wl**2*nue)
                n_InPd = sqrt( 0.5 *(abs(eps) + eps.real))
                k_InPd = sqrt( 0.5 *(abs(eps) - eps.real))
                self.stratumRIndexes[q] = n_InPd + 1j*k_InPd
            elif material == 'GaAs':
                nue = 1
                me0 = cst['GaAs'].me0
                a = 8.97E-5*wl**2/me0*self.stratumDopings[q]
                eps = n_GaAs**2 - a / (1+1j*5.305e-3*wl**2*nue)
                n_GaAsd = sqrt( 0.5 *(abs(eps) + eps.real))
                k_GaAsd = sqrt( 0.5 *(abs(eps) - eps.real))
                self.stratumRIndexes[q] = n_GaAsd + 1j*k_GaAsd
            elif material == 'InGaAs':
                xFrac = self.stratumCompositions[q]
                # TODO: bowing parameters?
                n_InGaAs = xFrac*n_InAs + (1-xFrac)*n_GaAs
                nue = 1
                me0 = xFrac*cst['InAs'].me0 + (1-xFrac)*cst['GaAs'].me0
                a = 8.97E-5*wl**2/me0*self.stratumDopings[q]
                eps = n_InGaAs**2 - a / (1+1j*5.305e-3*wl**2*nue)
                n_InGaAs = sqrt( 0.5 *(abs(eps) + eps.real))
                k_InGaAs = sqrt( 0.5 *(abs(eps) - eps.real))
                self.stratumRIndexes[q] = n_InGaAs + 1j*k_InGaAs
            elif material == 'InAlAs':
                xFrac = self.stratumCompositions[q]
                # TODO: bowing parameters?
                n_AlInAs = (1-xFrac)*n_AlAs + xFrac*n_InAs
                nue = 1
                me0 = (1-xFrac)*cst['AlAs'].me0 + xFrac*cst['InAs'].me0
                a = 8.97E-5*wl**2/me0*self.stratumDopings[q]
                eps = n_AlInAs**2 - a / (1+1j*5.305e-3*wl**2*nue)
                n_AlInAs = sqrt( 0.5 *(abs(eps) + eps.real))
                k_AlInAs = sqrt( 0.5 *(abs(eps) - eps.real))
                self.stratumRIndexes[q] = n_AlInAs + 1j*k_AlInAs
            elif material == 'Au':
                C1=-0.1933; C2=0.3321; C3=0.0938;
                D1=-0.382; D2=6.8522; D3=-0.1289;
                n_Au = C1+wl*C2+wl*C3**2
                k_Au = D1+wl*D2+wl*D3**2
                self.stratumRIndexes[q] = n_Au+k_Au*1j
            elif material == 'SiNx':
                #from Jean Nguyen's Thesis
                C1 = 2.0019336; C2 = 0.15265213; C3 = 4.0495557
                D0=-0.00282; D1=0.003029; D2=-0.0006982; D3=-0.0002839;
                D4=0.0001816; D5=-3.948e-005; D6=4.276e-006; D7=-2.314e-007; D8=4.982e-009
                n_SiNx = C1 + C2/wl**2 + C3/wl**4
                k_SiNx = D0 + D1*wl + D2*wl**2 + D3*wl**3 + D4*wl**4 \
                        + D5*wl**5 + D6*wl**6 + D7*wl**7 + D8*wl**8
                k_SiNx *= 100
                self.stratumRIndexes[q] = n_SiNx+k_SiNx*1j
            elif material == 'SiO2':
                #from Jean Nguyen's Thesis
                C1 = 1.41870; C2 = 0.12886725; C3 = 2.7573641e-5
                n_SiO2 = C1 + C2/wl**2 + C3/wl**4
                
                #this is a 4 peak Lorentzian fit to her data
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
                
    def get_nCore(self, data):
        """Get equiv. overall core relection index and decay? for InAs/GaAs layers"""

        # Matrial reflection index for GaAs, InAs, AlAs
        # Repeated codes, should be moved outside 1
        wl = self.wavelength
        n_GaAs = sqrt(cst['GaAs'].C1 + cst['GaAs'].C2*wl**2/(wl**2-cst['GaAs'].C3**2) 
                + cst['GaAs'].C4*wl**2/(wl**2-cst['GaAs'].C5**2) )
        n_InAs = sqrt(cst['InAs'].C1 + cst['InAs'].C2*wl**2/(wl**2-cst['InAs'].C3**2) 
                + cst['InAs'].C4*wl**2/(wl**2-cst['InAs'].C5**2) )
        n_AlAs = sqrt(cst['AlAs'].C1 + cst['AlAs'].C2*wl**2/(wl**2-cst['AlAs'].C3**2) 
                + cst['AlAs'].C4*wl**2/(wl**2-cst['AlAs'].C5**2) )
        
        n=np.zeros(8)
        # TODO: bowing parameters? 
        for i in range(0, 8, 2):
            n[i] = data.moleFrac[i]*n_InAs + (1-data.moleFrac[i])*n_GaAs
        for i in range(1, 8, 2):
            n[i] = data.moleFrac[i]*n_InAs + (1-data.moleFrac[i])*n_AlAs
        nCore = sum(data.h*n)/sum(data.h) # Average n?
        
        kCore = 1/(4*pi) * self.aCore * wl*1e-4 
        # See Def of acore
        # 1e-4: aCore in cm-1, wl in um
        
        return nCore+kCore*1j
        
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
        #use rounding to work with selected resolution
        self.stratumThickNum = np.round(self.stratumThicknesses
                /self.xres).astype(np.int64)
        self.stratumThicknesses = self.stratumThickNum * self.xres
        
        #convert to int to prevent machine rounding errors
        self.xPoints = self.xres * np.arange(0, self.stratumThickNum.sum(), 1)
        
        stratumThickNumCumSum = np.concatenate( ([0], 
            self.stratumThickNum.cumsum()) )
        self.xn = np.zeros(self.xPoints.size, dtype=complex)
        self.xAC = np.zeros(self.xPoints.size) #binary designation for Active Core

        #extend layer data for all xpoints
        for q in xrange(0,self.stratumThickNum.size):
            self.xn[stratumThickNumCumSum[q] :
                    stratumThickNumCumSum[q+1]] = self.stratumRIndexes[q]
            if self.stratumMaterials[q] == 'Active Core':
                self.xAC[stratumThickNumCumSum[q] : 
                        stratumThickNumCumSum[q+1]] = 1
                
        # make array to show selected stratum in mainCanvas
        # ? what is this used for? 
        self.xStratumSelected = np.zeros(self.xPoints.shape) * np.NaN
        if self.stratumSelected != -1: #if not no row selected
            self.xStratumSelected[stratumThickNumCumSum[self.stratumSelected]
                        : stratumThickNumCumSum[self.stratumSelected+1] ] \
                   = self.xn.real[stratumThickNumCumSum[self.stratumSelected]
                        : stratumThickNumCumSum[self.stratumSelected+1] ]

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
        phi   = k*self.stratumThicknesses*alpha
        #zeta  = k*self.stratumThicknesses/z0
        
        Mj = []
        M = np.array([[1+0j,0],[0,1+0j]])
        for q in xrange(alpha.size):
            Mj.append(np.array([[cos(phi[q]), -1j/gamma[q]*sin(phi[q])],[-1j*gamma[q]*sin(phi[q]), cos(phi[q])]]))
        
        Mj.reverse()
        for mj in Mj:
            M = np.dot(mj,M)
            
        gammas = gamma[0]
        gammac = gamma[-1]
        
        chi = gammac*M[0,0] + gammac*gammas*M[0,1] + M[1,0] + gammas*M[1,1]
        return chi
        
    def beta_find(self, betaInit = None):
        #? Seems to relate to EMF mode
        if True: #betaInit == None:
            betaMax  = max(self.stratumRIndexes.real)
            betaMin  = min(self.stratumRIndexes.real)

            betas = np.arange(betaMin.real+0.01,betaMax.real,0.01)

            if True: #do chi_find in c
                chiImag = np.zeros(len(betas),dtype=float)
                betasReal = betas.real
                betasImag = betas.imag
                stratumRIndexesReal = self.stratumRIndexes.real.copy()
                stratumRIndexesImag = self.stratumRIndexes.imag.copy()
                cFunctions.chiImag_array(c_double(self.wavelength), 
                        self.stratumThicknesses.ctypes.data_as(c_void_p), 
                        stratumRIndexesReal.ctypes.data_as(c_void_p), 
                        stratumRIndexesImag.ctypes.data_as(c_void_p), 
                        int(self.stratumRIndexes.size), 
                        betasReal.ctypes.data_as(c_void_p), 
                        betasImag.ctypes.data_as(c_void_p),
                        int(betasReal.size), chiImag.ctypes.data_as(c_void_p))
                beta0s = zero_find(betas.real, chiImag)
            else:
                chi=np.zeros(betas.size, dtype=complex)
                for p, beta in enumerate(betas):
                    chi[p] = self.chi_find(beta)
                beta0s = zero_find(betas.real, chi.imag)
            beta = max(beta0s)+1j*min(self.stratumRIndexes.imag)
        else:
            beta = betaInit
        
        beta_find_precision = 1e-5
        if True: #setting to True makes the function stall in Mac OS X
            betaIn = beta
            stratumRIndexesReal = self.stratumRIndexes.real.copy()
            stratumRIndexesImag = self.stratumRIndexes.imag.copy()
            betaOut = np.array([0.0, 0.0])
            beta = cFunctions.beta_find(c_double(self.wavelength), 
                    self.stratumThicknesses.ctypes.data_as(c_void_p), 
                    stratumRIndexesReal.ctypes.data_as(c_void_p), 
                    stratumRIndexesImag.ctypes.data_as(c_void_p), 
                    int(self.stratumRIndexes.size), c_double(betaIn.real), 
                    c_double(betaIn.imag), c_double(beta_find_precision), 
                    betaOut.ctypes.data_as(c_void_p))
            beta = betaOut[0] + 1j*betaOut[1]
        else:
            rInc = 0.0001; iInc = 1j*1e-6
            abschiNew=1
            while True:
                betas = [beta, beta+rInc, beta-rInc, beta+iInc, beta-iInc, 
                         beta+rInc+iInc, beta-rInc-iInc, beta+rInc-iInc,
                         beta-rInc+iInc]
                if True:
                    chi = np.zeros(len(betas),dtype=complex)
                    for p, betaIn in enumerate(betas):
                        chi[p] = self.chi_find(betaIn)
                else: #do chi_find in c
                    chi = np.zeros(len(betas),dtype=float)
                    abschi_find = cFunctions.abschi_find
                    abschi_find.restype = c_double
                    for p, betaIn in enumerate(betas):
                        stratumRIndexesReal = self.stratumRIndexes.real.copy()
                        stratumRIndexesImag = self.stratumRIndexes.imag.copy()
                        chi[p] = abschi_find(c_double(self.wavelength), 
                                self.stratumThicknesses.ctypes.data_as(c_void_p), 
                                stratumRIndexesReal.ctypes.data_as(c_void_p), 
                                stratumRIndexesImag.ctypes.data_as(c_void_p), 
                                int(self.stratumRIndexes.size), c_double(betaIn.real), 
                                c_double(betaIn.imag))
                abschiOld = abschiNew
                abschiNew = min(abs(chi))
                idx=argmin(abs(chi))
                beta=betas[idx]
                if abs(abschiOld -abschiNew)/abschiOld < beta_find_precision:
                    break
        return beta
        
    def mode_plot(self):
        #  print "-----debug----- mode_plot is called"
        n=copy.copy(self.stratumRIndexes)[::-1]
        thicknesses = copy.copy(self.stratumThicknesses)[::-1]
        ThickNum = copy.copy(self.stratumThickNum)[::-1]
        #  xres = self.xres        
        
        z0 = 0.003768
        #z0 = 376.8
        k = 2*pi/self.wavelength
        M = np.array([[1+0j,0],[0,1+0j]])
        
        alpha = sqrt(n**2-self.beta**2)
        if alpha[0].imag < 0:
            alpha[0] = conj(alpha[0])
        if alpha[-1].imag < 0:
            alpha[-1] = conj(alpha[-1])
        gamma = z0*alpha/n**2
        phi   = k*thicknesses*alpha
        #zeta  = k*thicknesses/z0
        
        ncs = stratumThickNumCumSum = np.concatenate( ([0], ThickNum.cumsum()) )
        xI = np.zeros(self.xPoints.size, dtype=complex)

        for q in xrange(ThickNum.size-1, -1, -1):
            xvec  = copy.copy(self.xPoints[ncs[q]:ncs[q+1]])[::-1]
            if len(xvec) == 0: #make sure xvec isn't empty
                continue
            xvec -= min(xvec)
            field = np.dot(M,np.array([1,gamma[-1]]))
            U = field[0]
            V = field[1]
            if q == 0 or q == self.stratumThicknesses.size-1:
                xI[ncs[q]:ncs[q+1]] = real(U)*exp(1j*k*alpha[q]*xvec) / n[q]**2
            else:
                xI[ncs[q]:ncs[q+1]]  = U*cos(-k*alpha[q]*xvec) \
                        + 1j/gamma[q] * V*sin(-k*alpha[q]*xvec)
                xI[ncs[q]:ncs[q+1]] /= n[q]**2
                Mj = np.array([[cos(phi[q]), -1j/gamma[q]*sin(phi[q])],
                    [-1j*gamma[q]*sin(phi[q]), cos(phi[q])]])
                M = np.dot(Mj,M)
        
        xI = abs(xI)**2 
        xI = xI / max(xI)
        self.xI = xI[::-1]
        
        #calculate confinement factor
        self.confinementFactor = sum(self.xI * self.xAC) / sum(self.xI)
        
    def calculate_performance_parameters(self):
        #waveguide loss
        self.waveguideLoss = 4 * pi * self.beta.imag \
                / (self.wavelength * 1e-6) * 1e-2
        
        #mirror loss
        self.mirrorLoss = -1 / (2 * self.waveguideLength * 0.1) \
                * log(self.frontFacet * self.backFacet)
        
        #transition cross-section
        Eph = h * c0 / (self.wavelength * 1e-6)
        neff = self.beta.real
        z = self.opticalDipole * 1e-10
        deltaE = 0.1*Eph
        sigma0 = 4*pi*e0**2 / (h*c0*eps0*neff) * Eph/deltaE * z**2
        
        #gain
        tauEff = self.tauUpper * (1 - self.tauLower 
                / self.tauUpperLower) * 1e-12
        Lp = self.Lp * 1e-10
        self.gain = sigma0 * tauEff / (e0 * Lp) #m/A
        self.gain *= 100 #cm/A
        
        #threshold current density
        self.Jth0 = (self.waveguideLoss + self.mirrorLoss) \
                / (self.gain * self.confinementFactor) #A/cm^2
        self.Jth0 *= 1e-3 #kA/cm^2
        
        #threshold current
        self.Ith0 = self.Jth0*1e3 * (self.Np * self.Lp*1e-8) \
                * self.waveguideLength*1e-1
        
        #operating voltage
        self.operatingVoltage = self.operatingField*1e3 * self.Lp*1e-8 * self.Np
        
        #voltage efficiency
        self.voltageEfficiency = 1.24/self.wavelength * self.Np \
                / self.operatingVoltage
        
        #extraction efficiency
        self.extractionEfficiency = self.mirrorLoss \
                / (self.mirrorLoss + self.waveguideLoss)
        
        #population inverstion efficiency
        tauEff = self.tauUpper * (1 - self.tauLower / self.tauUpperLower)
        self.inversionEfficiency = tauEff / (tauEff + self.tauLower)
        
        #modal efficiency
        xI = self.xI/max(self.xI)
        U = xI[np.nonzero(self.xAC)[0]]
        
        #interoplate over U_AC for each Np at the point xbar for Ubar
        #this is Faist's version
        numACs = self.stratumMaterials.count('Active Core')
        try:
            xVals = np.arange(self.xres,numACs*self.Np*self.Lp*1e-4,self.xres)
            assert xVals.size == U.size
        except AssertionError:
            try:
                xVals = np.arange(0,numACs*self.Np*self.Lp*1e-4,self.xres)
                assert xVals.size == U.size
            except AssertionError:
                xVals = np.arange(self.xres,
                        numACs*self.Np*self.Lp*1e-4-self.xres,self.xres)
                assert xVals.size == U.size
        tck = interpolate.splrep(xVals,U,s=0)
        minx = 0.5*self.Lp*1e-4
        maxx = numACs*self.Np*self.Lp*1e-4-0.5*self.Lp*1e-4
        xbar = np.linspace(minx, maxx, numACs*self.Np)
        Ubar = interpolate.splev(xbar,tck,der=0)
        self.modalEfficiency = sum(Ubar)**2 / (numACs * self.Np * sum(Ubar**2))
        
        #Kale's version
        modalEfficiency = sum(Ubar) / self.Np #since Ubar taken from normalized xI
        #I guess we'll go with Faist's version since he probably knows better than I do.
       

class QCLayers(object):
    """Class for QCLayers
    Member variables: 
        parameters for each layer, np.array type, with len = No. of layers: 
            layerWidths -float in angstrom, width of each layer
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
        TempFoM: ??? seems not used
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
        moleFrac: mole fraction for each possible layer material, in format [
    """
    def __init__(self):
        self.layerWidths = np.array([1.,1.]) # angstrom
        self.layerBarriers = np.array([0,0]) # boolean
        self.layerARs = np.array([0,0])      # boolean
        self.layerMaterials = np.array([1,1]) #label
        self.layerDopings = np.array([0.,0.]) #1e17 cm-3
        self.layerDividers = np.array([0,0]) #?
        
        self.xres = 0.5 # angstrom
        self.EField = 0 # kV/cm = 1e5 V/m
        self.layerSelected = -1 # int
        self.vertRes = 0.5 # meV
        self.repeats = 2 # repeats n times for the given structure
        
        self.description = ""
        self.solver = "SolverH" #?
        self.Temperature = 300
        self.TempFoM = 300 # ?
        self.diffLength = 0 # ?
        self.basisARInjector = True
        self.basisInjectorAR = True
        #  self.designByAngs = True
        #  self.designByML = False
        self.substrate = 'InP'
        
        #  self.moleFrac1 = 0.53
        #  self.moleFrac2 = 0.52
        #  self.moleFrac3 = 0.53
        #  self.moleFrac4 = 0.52
        #  self.moleFrac5 = 0.53
        #  self.moleFrac6 = 0.52
        #  self.moleFrac7 = 0.53
        #  self.moleFrac8 = 0.52
        self.moleFrac = [0.53, 0.52, 0.53, 0.52, 0.53, 0.52, 0.53, 0.52]
        
        self.update_alloys()
        self.update_strain()
        self.populate_x()

    def populate_x(self):
        """Extend layer information to position functions
        Layer data: layerWidth
                    layerNum
                with len = # of layers and each value repr. a layer
        Position data (OUTPUT/update member variables): 
                    xPoints 
                       position grid
                    xBarriers: from layerBarriers, is barrier layer
                        should be boolean (TBD)
                    xARs: from layerARs, is active region
                        should be boolean (TBD)
                    xMaterials: from layerMaterials label/index of material
                        should be int starting from 0 (TBD)
                    xDopings: from layerDopings, doping per volumn
                    xLayerNums 
                       at xPoints[q] it's xLayerNums[q]-th layer"""
        #  print "-----debug----- QCLayers populate_x called"
        #  print self.layerBarriers
        #use rounding to work with selected resolution
        self.layerNum = np.round(self.layerWidths /
                self.xres).astype(np.int64)
        self.layerWidths = self.layerNum * self.xres
        
        #convert to int to prevent machine rounding errors
        self.xPoints = self.xres * np.arange(0, self.layerNum.sum())
        
        #  layerWidthsCumSum = np.concatenate([[0.],self.layerWidths.cumsum()])
        layerNumCumSum = np.concatenate( ([0], self.layerNum.cumsum()) )
        self.xBarriers = np.zeros(self.xPoints.shape)
        self.xARs = np.zeros(self.xPoints.shape)
        self.xMaterials = np.zeros(self.xPoints.shape)
        self.xDopings = np.zeros(self.xPoints.shape)
        self.xLayerNums = np.zeros(self.xPoints.shape)

        #extend layer data for all xpoints
        for q in xrange(0,self.layerWidths.size):
            self.xBarriers[ layerNumCumSum[q] : 
                    layerNumCumSum[q+1] ] = self.layerBarriers[q]
            if self.layerARs[q] == 1:
                self.xARs[ layerNumCumSum[q]-1 : layerNumCumSum[q+1]+1] = 1
            self.xMaterials[layerNumCumSum[q] : 
                    layerNumCumSum[q+1] ] = self.layerMaterials[q]
            self.xDopings[layerNumCumSum[q] : 
                    layerNumCumSum[q+1] ] = self.layerDopings[q]
            self.xLayerNums[layerNumCumSum[q] : layerNumCumSum[q+1] ] = q

            
        #plt.plot(self.xPoints, self.xBarriers,'o')

        #duplicate layer based on user input repeats
        #  repeats = self.repeats
        if self.repeats >= 2:
            self.xPoints = np.arange(0, self.layerWidths.sum()
                + self.layerWidths[1:].sum()*(self.repeats-1), self.xres)
            self.xBarriers = np.hstack(   [self.xBarriers, np.tile(
                self.xBarriers[layerNumCumSum[1]:], self.repeats-1)  ])
            self.xARs = np.hstack(        [self.xARs, np.tile(
                self.xARs[layerNumCumSum[1]:], self.repeats-1)       ])
            self.xMaterials = np.hstack(  [self.xMaterials, np.tile(
                self.xMaterials[layerNumCumSum[1]:], self.repeats-1) ])
            self.xDopings = np.hstack(    [self.xDopings, np.tile(
                self.xDopings[layerNumCumSum[1]:], self.repeats-1)   ])
            self.xLayerNums = np.hstack(  [self.xLayerNums, np.tile(
                self.xLayerNums[layerNumCumSum[1]:], self.repeats-1) ])
        
        
        #this hack is needed because sometimes self.xPoints is one element too big
        self.xPoints = self.xPoints[0 : self.xBarriers.size]
        
        self.update_strain()
        # Following are equiv. elec potential for different bands
        # external field is included
        self.xVc  = np.zeros(self.xPoints.size)
        self.xVX  = np.zeros(self.xPoints.size)
        self.xVL  = np.zeros(self.xPoints.size)
        self.xVLH = np.zeros(self.xPoints.size)
        self.xVSO = np.zeros(self.xPoints.size)
        for MLabel in range(1,5):
            indx = np.nonzero(self.xMaterials == MLabel)[0]
            if indx.size != 0:
                material = np.where(self.xBarriers[indx] == 1, 
                        MLabel*2-1, (MLabel-1)*2)
                self.xVc[indx]  = self.EcG[material] \
                        - self.xPoints[indx] * ANG * self.EField * KVpCM
                self.xVX[indx]  = self.EcX[material] \
                        - self.xPoints[indx] * ANG * self.EField * KVpCM
                self.xVL[indx]  = self.EcL[material] \
                        - self.xPoints[indx] * ANG * self.EField * KVpCM
                self.xVLH[indx] = self.EvLH[material] \
                        - self.xPoints[indx] * ANG * self.EField * KVpCM
                self.xVSO[indx] = self.EvSO[material] \
                        - self.xPoints[indx] * ANG * self.EField * KVpCM
            
        #  self.xVc = self.xBarriers * 0.520 - self.xPoints * self.EField * 1e-5
        
        # make array to show selected layer in mainCanvas
        try:
            self.xLayerSelected = np.zeros(self.xPoints.shape)*np.NaN
            layerSelected = self.layerSelected
            if layerSelected != -1:
                if layerSelected == 0: #row for first layer is selected
                    for repeat in xrange(1,self.repeats+1):
                        base = layerNumCumSum[-1] * (repeat - 1)
                        self.xLayerSelected[base+layerNumCumSum[layerSelected]
                                :base+layerNumCumSum[layerSelected+1]+1] \
                            = self.xVc[     base+layerNumCumSum[layerSelected]
                                :base+layerNumCumSum[layerSelected+1]+1]
                elif self.layerSelected == self.layerWidths.size: 
                    #last (blank) layer row is selected
                    pass
                else: 
                    for repeat in xrange(1,self.repeats+1):
                        base = sum(self.layerNum[1:])*(repeat-1)
                        self.xLayerSelected[base+layerNumCumSum[layerSelected]-1
                                :base+layerNumCumSum[layerSelected+1]+1]\
                            = self.xVc[base+layerNumCumSum[layerSelected]-1
                                :base+layerNumCumSum[layerSelected+1]+1]
        except IndexError:
            #index error happens in SolveBasis when the selected layer is greater than the number of layers in the solve swath
            # however, xLayerSelected is not used for the SolveBasis function
            print "Index Error for layer selection at function \
            qclayer.populate_x"

        self.xARs[np.nonzero(self.xARs==0)[0]] = np.NaN
        self.xARs *= self.xVc

    def populate_x_band(self):
        """Extend layer information to position functions for band parameter
        OUTPUT/update member variables):
            xEg, xMc, xESO, xEp, xF, 
            whose value is determeined by the layer material
        """
        #  print "------debug------- QCLayers populate_x_band called"
        # Following parameters can be looked up in cFunctions.c
        self.xEg = np.zeros(self.xPoints.size)  
        self.xMc = np.zeros(self.xPoints.size)  # Seems not to be used
        self.xESO = np.zeros(self.xPoints.size) 
        self.xEp = np.zeros(self.xPoints.size)
        self.xF = np.zeros(self.xPoints.size)
        for MLabel in range(1,5):
            indx = np.nonzero(self.xMaterials == MLabel)[0]
            if indx.size !=0 :
                material = np.where(self.xBarriers[indx] == 1, 
                        MLabel*2-1, (MLabel-1)*2)
                self.xEg[indx] = self.EgLH[material]
                self.xMc[indx] = self.me[material]
                self.xESO[indx]= self.ESO[material]
                self.xEp[indx] = self.Ep[material]
                self.xF[indx]  = self.F[material]
        #  self.xMc = self.xMc #??????

    def update_alloys(self):  # c is a Material_Constant class instance
        """
        update material parameter for the alloy used.
        """
        variables = ['EgG', 'EgL', 'EgX', 'VBO', 'DSO', 'me0', 'acG', 'acL', 'acX', 
                 'Ep', 'F', 'XiX', 'b', 'av', 'alG', 'beG', 'alL', 'beL', 'alX', 
                 'beX', 'epss', 'epsInf', 'hwLO', 'alc', 'c11', 'c12']
        #  print "----debug--- substrate is "+self.substrate
        if self.substrate == 'InP':
            self.numMaterials = 8
            Mat1 = ['InAs']*8
            Mat2 = ['GaAs', 'AlAs']*4
            MatCross = ['InGaAs', 'AlInAs']*4
        elif self.substrate == 'GaAs': 
            self.numMaterials = 8
            Mat1 = ['AlAs']*8
            Mat2 = ['GaAs']*8
            MatCross = ['AlGaAs']*8 # Note EgG_AlGaAs's moleFrad deps
        elif self.substrate == 'GaSb': 
            self.numMaterials = 8
            Mat1 = ['InAs', 'AlSb']*4
            Mat2 = ['InSb', 'GaSb']*4
            MatCross = ['InAsSb', 'AlGaSb']*4 # Note EgG's bowing moleFrad deps
        else: 
            raise TypeError('substrate selection not allowed')

        self.h = np.zeros(self.numMaterials)
        for item in variables:
            setattr(self, item, np.empty(self.numMaterials))
            para = getattr(self, item)
            for n in range(self.numMaterials):
                para[n] = self.moleFrac[n]*getattr(cst[Mat1[n]], item) \
                    + (1-self.moleFrac[n])*getattr(cst[Mat2[n]], item)
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

        #set this once the others are set        
        self.epsrho = 1 / (1/self.epsInf - 1/self.epss)     

    def update_strain(self):  # c is a Material_Constant class instance

        if self.substrate == 'InP':
            self.a_parallel = cst['InP'].alc
        elif self.substrate == 'GaSb':
            self.a_parallel = cst['GaSb'].alc
        elif self.substrate == 'GaAs':
            self.a_parallel = cst['GaAs'].alc
        else:
            raise TypeError('substrate selection not allowed')
            
        for i in range(4): 
            indx = np.nonzero(self.layerMaterials[1:] == i+1)[0]
            # [1:] because first layer doesn't count
            self.h[2*i+1] = sum(self.layerWidths[indx]
                    * self.layerBarriers[indx])
            self.h[2*i] = sum(self.layerWidths[indx]) - self.h[2*i+1]
        
        self.eps_parallel = self.a_parallel / self.alc - 1; #Walle eqn 1b
        
        self.a_perp   = self.alc * (1 - 2* self.c12 / self.c11 * self.eps_parallel); #Walle eqn 2a
        self.eps_perp = self.a_perp/self.alc - 1; #Walle eqn 2b
        #             = -2*self.c12/self.c11*self.eps_parallel
        
        #  print "------debug-----", sum(self.h)
        self.mismatch = 100 * sum(self.h*self.eps_perp) / sum(self.h);
        
        self.MLThickness = np.zeros(self.layerMaterials.size)
        for n, (MLabel, BLabel) in enumerate( zip((1,1,2,2,3,3,4,4),
            (0,1)*4)): 
            #? what's MLThickness???
            self.MLThickness[np.nonzero( 
                  (self.layerMaterials == MLabel) 
                & (self.layerBarriers == BLabel)
                )[0] ] = self.a_perp[n] / 2.0
    
        # for In0.53Ga0.47As, EcG = 0.22004154
        #    use this as a zero point baseline
        baseln = 0.22004154
        
        # Pikus-Bir interaction correction to bands offset, 
        # According to Kale's, Eq.(2.14), Pe for \delta E_{v}
        self.Pec = (2*self.eps_parallel+self.eps_perp) * (self.acG)
        #  self.Pe = 2*self.av * (self.c11-self.c12) / self.c11 * self.eps_parallel
        self.Pe  = (2*self.eps_parallel+self.eps_perp) * (self.av)
        # Kale's Thesis, Eq.(2.16)
        self.Qe = - self.b * (self.c11+2*self.c12) / self.c11 * self.eps_parallel
        # temperature correction to conduction band edge, Eq.(2.10) in Kale's 
        self.Varsh = - self.alG*cst.Temperature**2/(cst.Temperature+self.beG)
        
        # calculations for effective mass
        #  in part following Sugawara, PRB 48, 8102 (1993)
        
        # Eq.(2.15) in Kale's
        self.ESO  = sqrt(9*self.Qe**2+2*self.Qe*self.DSO+self.DSO**2);
        self.EgLH = self.EgG + self.Pec + self.Pe + self.Varsh \
                - 1/2*(self.Qe - self.DSO + self.ESO);
        self.EgSO = self.EgG + self.Pec + self.Pe + self.Varsh \
                - 1/2*(self.Qe - self.DSO - self.ESO);
        
        #  self.alpha = (self.Qe + self.DSO + 
                #  sqrt(9*self.Qe.^2+2*self.Qe*self.DSO+self.DSO.^2)) \
                #  / (2*sqrt(9*self.Qe.^2+2*self.Qe*self.DSO+self.DSO.^2));
        #  self.beta = (-self.Qe - self.DSO + 
                #  sqrt(9*self.Qe.^2+2*self.Qe*self.DSO+self.DSO.^2)) \
                #  / (2*sqrt(9*self.Qe.^2+2*self.Qe*self.DSO+self.DSO.^2));
        #  self.alpha = self.alpha / sqrt(self.alpha.^2 + self.beta.^2);
        #  self.beta = self.beta / sqrt(self.alpha.^2 + self.beta.^2);
        self.me = 1 / ((1+2*self.F) + self.Ep/self.EgLH*(self.EgLH+2/3*self.ESO)/(self.EgLH + self.ESO));
        #  self.me2 = 1 / ((1+2*self.F) + self.Ep/3 *( 
            #  (sqrt(2)*self.alpha-self.beta).^2 / self.EgLH 
            #  + (sqrt(2)*self.beta+self.alpha).^2 / (self.EgLH + self.ESO))); #using Igor's Ep instead of the calculated P2
        # self.me3 = 1 / ((1+2*self.F) + self.Ep/3 *( 2 / self.EgLH + 1 / (self.EgLH + self.ESO)));
        # self.me4 = 1 / ((1+2*self.F) + self.Ep/3 *( 2 / self.EgLH + 1 / (self.EgSO)));
        
        #corrections to the method used to calculate band edges, thanks to Yu Song
        # band edge at different point, Eq.(2.7)
        self.EcG = self.VBO + self.EgG + self.Pec - baseln; # Varsh?
        # band edge at L and X?
        self.EcL = self.VBO + self.EgL \
                + (2*self.eps_parallel+self.eps_perp) * (self.acL+self.av) - baseln;
        self.EcX = self.VBO + self.EgX \
                + (2*self.eps_parallel+self.eps_perp)*(self.acX+self.av) \
                + 2/3 * self.XiX * (self.eps_perp-self.eps_parallel) - baseln;
        
        # Same as previous lines except for Varsh.. but why?
        self.EgLH = self.EgG + self.Pec + self.Pe \
                - 1/2*(self.Qe - self.DSO + self.ESO);
        self.EgSO = self.EgG + self.Pec + self.Pe \
                - 1/2*(self.Qe - self.DSO - self.ESO);
        
        #1st MAJOR assumption: 
        #   Varshney contribution to band edge is in proportion to percent 
        #   of band offset
        #2nd major assumption: 
        #   Temperature affects sattelite valleys in the same way it does 
        #   the Gamma valley
        barrs = np.array([1,3,5,7])
        wells = np.array([0,2,4,6])
        CBOffset = self.EcG[barrs] - self.EcG[wells]
        VBOffset = (self.EcG[barrs] - self.EgLH[barrs]) - (self.EcG[wells] - self.EgLH[wells])
        percentCB = CBOffset / (CBOffset + VBOffset)
        percentCB = np.column_stack([percentCB,percentCB]).flatten() 
        #applies percent CV to both well and barrier slots
        
        self.EcG += percentCB * self.Varsh
        self.EcL += percentCB * self.Varsh
        self.EcX += percentCB * self.Varsh
        self.EvLH = self.EcG - self.EgLH - ((1-percentCB) * self.Varsh)
        self.EvSO = self.EcG - self.EgSO - ((1-percentCB) * self.Varsh)
        
        
        #other code that I've junked for some reason or other
        # self.EgAV = self.EgG - self.alG*c.Temperature.^2/(c.Temperature+self.beG) + (2*self.eps_parallel+self.eps_perp) * (self.acG+self.av);
        # self.EgEG = self.EgAV + 1/3*self.ESO;
        # self.P2 = 1/2 * (1/self.me0 - (1+2*self.F)/1) * self.EgG * (self.EgG+self.DSO) / (self.EgG+2/3*self.DSO); #strain independent
        # self.P2 = 1/2 * (1/self.me0 - (1+2*self.F)/1) * self.EgLH * (self.EgLH+self.ESO) / (self.EgLH+2/3*self.ESO);
        # Pe = 2*self.acG * (self.c11-self.c12) / self.c11 * self.eps_parallel;
        # Elh = -Pe + 1/2 * (self.Qe - self.DSO + sqrt(9*self.Qe.^2+2*self.Qe*self.DSO+self.DSO.^2))
        # Elh = -Pe + 1/2 * (self.Qe - self.DSO + sqrt(9*self.Qe.^2+2*self.Qe*self.DSO+self.DSO.^2));
        # Eso = -Pe + 1/2 * (self.Qe - self.DSO - sqrt(9*self.Qe.^2+2*self.Qe*self.DSO+self.DSO.^2));
        # self.me = 1 / ((1+self.F) + 2*self.P2/3 *( (sqrt(2)*self.alpha-self.beta).^2 / self.EgLH + (sqrt(2)*self.beta+self.alpha).^2 / self.EgSO));
        # self.me3 = 1 / ((2+2*self.F) + self.Ep*(self.EgLH + 2/3*self.ESO) / self.EgLH / (self.EgLH + self.ESO) ); #using Igor's effective mass formula
        # self.me4 = self.EgAV / self.Ep;

    def solve_psi(self):
        """ solve eigen mode
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
        """
        Epoints = np.arange(min(self.xVc),
                max(self.xVc-115*self.EField*1e-5), #?115e-5?
                self.vertRes/1000)
        xMcE = np.zeros(self.xPoints.shape)
        xPsi = np.zeros(self.xPoints.shape)
        psiEnd = np.zeros(Epoints.size)
        
        #TODO: add adaptive spacing for Eq
        #TODO: convert nested for loop to C
        if USE_CLIB:
            # Call C function to get boundary dependence of energy EPoints[n], 
            # the return value is psiEnd[n]
            # for n with psiEnd[n]=0, EPoints[n] is eigenenergy
            cFunctions.psiFnEnd(Epoints.ctypes.data_as(c_void_p), 
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
            #  psiFnEnd(double *eEq, int eEqSize, int xPsiSize, double xres,
                    #  double *xVc, double *xEg, double *xF, double *xEp, double
                    #  *xESO, double *xMcE, double *xPsi, double *xPsiEnd) 
            # my_sum.sum(a.ctypes.data_as(c_void_p), int(10))
        else:
            for p, Eq in enumerate(Epoints):
                if True:
                    xMcE = m0 / (1+2*self.xF + self.xEp/3 * 
                            (2 / ((Eq-self.xVc)+self.xEg) + 
                                1 / ((Eq-self.xVc)+self.xEg+self.xESO) ))
                else:
                    xMcE = self.xMc * (1 - (self.xVc - Eq) / self.xEg)
                xMcE[0:-1] = 0.5 * (xMcE[0:-1]+xMcE[1:])
                xPsi[0] = 0
                xPsi[1] = 1
                for q in xrange(1,xPsi.size-1):
                    xPsi[q+1] = xMcE[q] *((2 * (self.xres*1e-10 / hbar)**2 
                        * (self.xVc[q] - Eq)*e0 + 1 / xMcE[q] + 1 / xMcE[q-1]) 
                        * xPsi[q] - xPsi[q-1] / xMcE[q-1]) 
                psiEnd[p] = xPsi[-1]
                
        #interpolate between solved-for E points        
        tck = interpolate.splrep(Epoints,psiEnd,s=0)
        #adds 100 points per solved-for E point
        xnew = np.linspace(Epoints[0],Epoints[-1],Epoints.size*1e2) 
        ynew = interpolate.splev(xnew,tck,der=0)
        
#        #plot interpolated points over solved-for E points
#        plt.figure()
#        plt.plot(Epoints,psiEnd,'x',xnew,ynew,'o-')
#        plt.show()
       
        #find Eigen Energies
        #This routine looks for all of the zero crossings, and then picks each one out
        #TODO: try to replace this part by zero_find function
        gtz = ynew > 0
        ltz = ynew < 0
        overlap1 = np.bitwise_and(gtz[0:-1],ltz[1:])
        overlap2 = np.bitwise_and(gtz[1:],ltz[0:-1])
        overlap  = np.bitwise_or(overlap1, overlap2)
        idxs = np.nonzero(overlap == True)[0]
        #need this to maintain compatibility with 32-bit and 64-bit systems
        idxs = idxs.astype(float) 
        self.EigenE = np.zeros(idxs.size)

        if USE_CLIB:
            # use inverse quadratic to get an approximation of zeros
            cFunctions.inv_quadratic_interp(xnew.ctypes.data_as(c_void_p), 
                    ynew.ctypes.data_as(c_void_p), 
                    idxs.ctypes.data_as(c_void_p), 
                    int(idxs.size), self.EigenE.ctypes.data_as(c_void_p))
        else:
            for q, idx in enumerate(idxs): # do quadratic interpolation
                x0=xnew[idx-1]; fx0=ynew[idx-1]
                x1=xnew[idx];   fx1=ynew[idx]
                x2=xnew[idx+1]; fx2=ynew[idx+1]
                d1=(fx1-fx0)/(x1-x0); d2=(fx2-fx1)/(x2-x1)
                #inverse quadratic interpolation
                x3 = x0*fx1*fx2/(fx0-fx1)/(fx0-fx2) \
                        + x1*fx0*fx2/(fx1-fx0)/(fx1-fx2) \
                        + x2*fx0*fx1/(fx2-fx0)/(fx2-fx1)
                self.EigenE[q] = x3
#                if abs(d1) > 1e15 and abs(d2) > 1e15:
#                    self.EigenE[q] = 0

        if True:
            # Near the above approximation result, 
            # try to get a more precise result
            for q in xrange(self.EigenE.size):
                # 100000 is an estimate for the precision of above
                # approximation
                # TODO: change the three calls of psiFn to loop
                approxwidth = self.vertRes/100000
                x0=self.EigenE[q]-approxwidth 
                x1=self.EigenE[q]
                x2=self.EigenE[q]+approxwidth
                
                cFunctions.psiFn(c_double(x0), int(1), int(xPsi.size), 
                        c_double(self.xres), 
                        self.xVc.ctypes.data_as(c_void_p),
                        self.xEg.ctypes.data_as(c_void_p), 
                        self.xF.ctypes.data_as(c_void_p), 
                        self.xEp.ctypes.data_as(c_void_p), 
                        self.xESO.ctypes.data_as(c_void_p), 
                        self.xMc.ctypes.data_as(c_void_p), 
                        xMcE.ctypes.data_as(c_void_p), 
                        xPsi.ctypes.data_as(c_void_p))
                fx0 = xPsi[-1]
                
                cFunctions.psiFn(c_double(x1), int(1), int(xPsi.size), 
                        c_double(self.xres), 
                        self.xVc.ctypes.data_as(c_void_p), 
                        self.xEg.ctypes.data_as(c_void_p), 
                        self.xF.ctypes.data_as(c_void_p), 
                        self.xEp.ctypes.data_as(c_void_p), 
                        self.xESO.ctypes.data_as(c_void_p), 
                        self.xMc.ctypes.data_as(c_void_p), 
                        xMcE.ctypes.data_as(c_void_p), 
                        xPsi.ctypes.data_as(c_void_p))
                fx1 = xPsi[-1]
                
                cFunctions.psiFn(c_double(x2), int(1), int(xPsi.size), 
                        c_double(self.xres), 
                        self.xVc.ctypes.data_as(c_void_p), 
                        self.xEg.ctypes.data_as(c_void_p), 
                        self.xF.ctypes.data_as(c_void_p), 
                        self.xEp.ctypes.data_as(c_void_p), 
                        self.xESO.ctypes.data_as(c_void_p), 
                        self.xMc.ctypes.data_as(c_void_p), 
                        xMcE.ctypes.data_as(c_void_p), 
                        xPsi.ctypes.data_as(c_void_p))
                fx2 = xPsi[-1]
                
                #  psiFn(double Eq, int startpoint, int xPsiSize, double xres, 
                        #  double *xVc, double *xEg, double *xF, double *xEp, 
                        #  double *xESO, double *xMc, double *xMcE, double *xPsi)
                
               
                d1=(fx1-fx0)/(x1-x0); d2=(fx2-fx1)/(x2-x1);
                #inverse quadratic interpolation
                x3 = x0*fx1*fx2/(fx0-fx1)/(fx0-fx2) \
                        + x1*fx0*fx2/(fx1-fx0)/(fx1-fx2) \
                        + x2*fx0*fx1/(fx2-fx0)/(fx2-fx1)
                self.EigenE[q] = x3

        #make array for Psi and fill it in
        if USE_CLIB:
            # with eigenenregy EigenE, here call C function to get wave
            # function
            self.xyPsi = np.zeros(self.xPoints.size*self.EigenE.size)
            cFunctions.psiFill(int(xPsi.size), c_double(self.xres),
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
            #self.xyPsi.resize(a.xPoints.size, a.EigenE.size)
            self.xyPsi = self.xyPsi.reshape(self.xPoints.size, 
                    self.EigenE.size, order='F')
        else:
            self.xyPsi = np.zeros((self.xPoints.size,self.EigenE.size))
            for p, Eq in enumerate(self.EigenE):
                if True:
                    xMcE = m0 / (1+2*self.xF + self.xEp/3 * (
                        2 / ((Eq-self.xVc)+self.xEg) + 1 / (
                             (Eq-self.xVc)+self.xEg+self.xESO) ))
                else:
                    xMcE = self.xMc * (1 - (self.xVc - Eq) / self.xEg)
                xMcE[0:-1] = 0.5 * (xMcE[0:-1]+xMcE[1:])
                xPsi[1] = 1
                for q in xrange(2,xPsi.size-1):
                    xPsi[q+1] = ((2 * (self.xres*1e-10 / hbar)**2 
                        * (self.xVc[q] - Eq)*e0 + 1 / xMcE[q] + 
                            1 / xMcE[q-1]) 
                        * xPsi[q] - xPsi[q-1] / xMcE[q-1]) * xMcE[q]
                psiInt = sum(xPsi**2 * (1+(Eq-self.xVc)/(Eq-self.xVc+self.xEg)))
                A = 1 / sqrt( self.xres * 1e-10 * psiInt)
                self.xyPsi[:,p] = A * xPsi
        
        #remove states that come from oscillating end points
        # looks like we should change -1 to -2 (following)
        psiEnd = self.xyPsi[-1,:]
        idxs = abs(psiEnd)<10
        #  idxs = np.nonzero(abs(psiEnd)<10)[0]
        #  psiEnd = self.xyPsi[-2,:]
        #  idxs = np.abs(psiEnd)<200/self.xres 
        # 200 depends on how precise we want about eigenenergy solver 
        # (TODO: more analysis and test about this value
        self.EigenE = self.EigenE[idxs]
        self.xyPsi = self.xyPsi[:,idxs]

        #4.5e-10 scales size of wavefunctions, arbitrary for nice plots
        self.xyPsiPsi = self.xyPsi*self.xyPsi*settings.wf_scale 

        #remove states that are smaller than minimum height (remove zero
        # solutions?)-test case not showing any effect
        # addresses states high above band edge
        #0.014 is arbitrary; changes if 4.5e-10 changes
        #  idxs = np.nonzero(self.xyPsiPsi.max(0) > 
                #  settings.wf_scale * settings.wf_min_height)[0] 
        #  idxs = self.xyPsiPsi.max(0) > settings.wf_scale*settings.wf_min_height
        #  self.EigenE = self.EigenE[idxs]
        #  self.xyPsi = self.xyPsi[:,idxs]
        #  self.xyPsiPsi = self.xyPsiPsi[:,idxs]
        
        self.xyPsiPsi2 = copy.deepcopy(self.xyPsiPsi)

        # implement pretty plot: 
        # remove long zero head and tail of the wave functions
        # test case shows on "solve whole"
        for q in xrange(self.EigenE.size):
            #0.0005 is arbitrary
            prettyIdxs = np.nonzero(self.xyPsiPsi[:,q] > 
                    settings.wf_scale * settings.pretty_plot_factor)[0] 
            #  prettyIdxs = self.xyPsiPsi[:, q] > \
                    #  settings.wf_scale * settings.pretty_plot_factor
            self.xyPsiPsi[0:prettyIdxs[0],q] = np.NaN
            self.xyPsiPsi[prettyIdxs[-1]:,q] = np.NaN
            #  print q,prettyIdxs
        #  print self.xyPsiPsi

        #decimate plot points: seems for better time and memory performance?
        idxs = np.arange(0, self.xPoints.size, 
                int(settings.plot_decimate_factor/self.xres), dtype=int)
        self.xyPsiPsiDec = np.zeros((idxs.size, self.EigenE.size))
        for q in xrange(self.EigenE.size):
            self.xyPsiPsiDec[:,q] = self.xyPsiPsi[idxs,q]
        self.xyPsiPsi = self.xyPsiPsiDec
        self.xPointsPost = self.xPoints[idxs]


def basisSolve(data):
    """ solve basis for the QC device, with each basis being eigen mode of 
    a seperate part of the layer structure
    INPUT: 
        data: a QCLayers class
    OUTPUT: 
        dCL: a list, each element is a QCLayers class, with layer structure 
              limited within a seperate sigle active/injection area
    """
    #TODO: move to QCLayers.. should be a member method
    #self.data.basisInjectorAR is 0-to-1
    #self.data.basisARInjector is 1-to-0
        
    #find all 0-to-1 & 1-to-0 transitions
    # where for all n, 
    # data.layerARs[zeroTOone[n]] = data.layerARs[oneTOzero[n]] = 0
    # but data.layerARs[zeroTOone[n]+1] = data.layerARs[oneTOzero[n]-1] = 1
    zeroTOone = []
    oneTOzero = []
    layerAR = np.insert(data.layerARs, 0, data.layerARs[-1])
    for q in xrange(0,layerAR.size-1):
        if layerAR[q] == 0 and layerAR[q+1] == 1:
            zeroTOone.append(q-1)
        if layerAR[q] == 1 and layerAR[q+1] == 0:
            oneTOzero.append(q)

    dividers = [0, data.layerARs.size-1]
    if data.basisInjectorAR:
        dividers += zeroTOone
    if data.basisARInjector:
        dividers += oneTOzero
    # This converts the list into a set, thereby removing duplicates,
    # and then back into a list.
    dividers = list(set(dividers)) 
    dividers.sort()
    
    #the first region is always the "wrap around" region
    #  layers = [range(dividers[q], dividers[q+1]+1) for q in
            #  range(len(dividers)-1)]
    
    #this is dataClassesList. 
    # it holds all of the Data classes for each individual solve section
    dCL = [] 

    #for first period only
    #this handles all of the solving        
    for n in range(len(dividers)-1):
        dCL.append(copy.deepcopy(data))
        dCL[n].repeats = 1
        
        #substitute proper layer characteristics into dCL[n]
        layer = range(dividers[n], dividers[n+1]+1)
        dCL[n].layerWidths = np.hstack([100, data.layerWidths[layer], 30])
        dCL[n].layerBarriers = np.hstack([1, data.layerBarriers[layer], 1])
        dCL[n].layerARs = np.hstack([0, data.layerARs[layer], 0])
        dCL[n].layerMaterials = np.hstack([data.layerMaterials[layer][0], 
            data.layerMaterials[layer], data.layerMaterials[layer][-1]])
        dCL[n].layerDopings = np.hstack([0, data.layerDopings[layer], 0])
        dCL[n].layerDividers = np.hstack([0, data.layerDividers[layer], 0])
       
        #update and solve
        dCL[n].update_alloys()
        dCL[n].update_strain()
        dCL[n].populate_x()
        dCL[n].populate_x_band()
        dCL[n].solve_psi()
        
        #caculate offsets
        dCL[n].widthOffset = sum(data.layerWidths[range(0,dividers[n])]) #- 100/data.xres
        dCL[n].fieldOffset = -(dCL[n].widthOffset-100) * ANG \
                * dCL[n].EField * KVpCM            
    
    solvePeriods = len(dCL)
    
    #create dCL's and offsets for repeat periods
    counter = solvePeriods
    if data.repeats > 1:
        for q in xrange(1,data.repeats):
            for p in xrange(0,solvePeriods):
                dCL.append(copy.deepcopy(dCL[p]))
                dCL[counter].widthOffset = sum(data.layerWidths[1:])*q \
                        + dCL[p].widthOffset #- 100/data.xres
                dCL[counter].fieldOffset = -(dCL[counter].widthOffset-100)*ANG \
                        * dCL[counter].EField * KVpCM
                counter += 1
    return dCL

def convert_dCL_to_data(data, dCL):
    """ post processng of dCL list
    INPUT: 
        data: orginal QCLayers class
        dCL: result of basisSolve(data)
    OUPUT:
        get wave function (dCL[n].xyPsi) and eigenenrgy (dCL[n].EigenE) 
        in dCL and update them in data, and format it in length for data
        and update data.xyPsiPsi
    """
    #count number of wavefunctions
    numWFs = sum([dC.EigenE.size for dC in dCL])

    data.xPointsPost = np.arange(-100, 
            data.xPoints[-1] + 30 + data.xres, 
            data.xres)

    data.xyPsi = np.zeros((data.xPointsPost.size, numWFs))
    data.xyPsiPsi = np.NaN*np.zeros(data.xyPsi.shape)
    data.EigenE = np.zeros(numWFs)
    data.moduleID = np.zeros(numWFs)
    counter = 0
    for n, dC in enumerate(dCL):
        for q in xrange(dC.EigenE.size):
            #  wf = np.NaN*np.zeros(data.xPointsPost.size)
            #  begin = int(dC.widthOffset/data.xres)
            #  end = begin + dC.xyPsiPsi2[:,q].size 
            #  wf[begin:end] = dC.xyPsiPsi2[:,q]
            #  data.xyPsiPsi[:,counter] = wf
            data.EigenE[counter] = dC.EigenE[q] + dC.fieldOffset
            
            data.moduleID[counter] = n
            wf = np.zeros(data.xPointsPost.size)
            begin = int(dC.widthOffset/data.xres)
            end = begin + dC.xyPsi[:,q].size
            wf[begin:end] = dC.xyPsi[:,q]
            data.xyPsi[:,counter] = wf
            data.xyPsiPsi[:,counter] = wf**2 * settings.wf_scale
            
            counter += 1
    # cut head and tial to promise the figure is in the right place?
    head = int(100/data.xres)
    tail = -int(30/data.xres)
    data.xPointsPost = data.xPointsPost[head:tail]
    data.xyPsi = data.xyPsi[head:tail]
    data.xyPsiPsi = data.xyPsiPsi[head:tail]
    
    #implement pretty plot
    # remove long zero head and tail of the wave functions
    for q in xrange(data.EigenE.size):
        prettyIdxs = np.nonzero(data.xyPsiPsi[:,q] > 
                settings.wf_scale * settings.pretty_plot_factor)[0] 
        #0.0005 is arbitrary
        data.xyPsiPsi[0:prettyIdxs[0],q] = np.NaN
        data.xyPsiPsi[prettyIdxs[-1]:,q] = np.NaN
        
    #sort by ascending energy
    sortID = np.argsort(data.EigenE)
    data.EigenE = data.EigenE[sortID]
    data.xyPsi = data.xyPsi[:,sortID]
    data.xyPsiPsi = data.xyPsiPsi[:,sortID]
    data.moduleID = data.moduleID[sortID]

#    #decimate plot points
#    idxs = np.arange(0,data.xPoints.size, int(settings.plot_decimate_factor/data.xres), dtype=int)
#    data.xyPsiPsiDec = np.zeros([idxs.size, data.EigenE.size])
#    for q in xrange(data.EigenE.size):
#        data.xyPsiPsiDec[:,q] = data.xyPsiPsi[idxs,q]
#    data.xyPsiPsi = data.xyPsiPsiDec
#    data.xPointsPost = data.xPoints[idxs]

def lo_phonon_time(data, upper, lower):
    if upper < lower:
        temp = upper
        upper = lower
        lower = temp
        
    psi_i = data.xyPsi[:,upper]
    psi_j = data.xyPsi[:,lower]
    E_i = data.EigenE[upper]
    E_j = data.EigenE[lower]

    if E_i-E_j-data.hwLO[0] < 0:
        return 1e20
        
    idxs_i = np.nonzero(psi_i >
            settings.wf_scale * settings.phonon_integral_factor)[0]
    idxs_j = np.nonzero(psi_j >
            settings.wf_scale * settings.phonon_integral_factor)[0]
    idx_first = [idxs_i[0], idxs_j[0]]
    idx_last  = [idxs_i[-1], idxs_j[-1]]
    if max(idx_first) > min(idx_last):
        return 1e20
    else:
        idx_first = min([idxs_i[0], idxs_j[0]])
        idx_last  = max([idxs_i[-1], idxs_j[-1]])
        psi_i = psi_i[idx_first:idx_last]
        psi_j = psi_j[idx_first:idx_last]
        xPoints = data.xPoints[idx_first:idx_last]

        xMcE_j = data.xMc * (1 - (data.xVc - E_j) / data.xEg)        
        #weight non-parabolic effective mass by probability density
        McE_j = m0*sum(xMcE_j[idx_first:idx_last] * psi_j**2) / sum(psi_j**2) 
        
        kl = sqrt(2*McE_j/hbar**2 * (E_i-E_j-data.hwLO[0])*e0)
        dIij = np.zeros(xPoints.size)
        for n in xrange(xPoints.size):
            x1 = xPoints[n]*1e-10
            x2 = xPoints*1e-10
            dIij[n] = sum(psi_i*psi_j * exp(-kl*abs(x1-x2)) 
                    * psi_i[n]*psi_j[n] * (data.xres*1e-10)**2)
        Iij = sum(dIij)
        inverse_tau = McE_j * e0**2 * data.hwLO[0]*e0/hbar * Iij \
                / (4 * hbar**2 * data.epsrho[0]*eps0 * kl)
        tau = 1e12/inverse_tau
        return tau

def dipole(data, upper, lower):
    if upper < lower:
        temp = upper
        upper = lower
        lower = temp
    psi_i = data.xyPsi[:,upper]
    psi_j = data.xyPsi[:,lower]
    E_i = data.EigenE[upper]
    E_j = data.EigenE[lower]
    
    data.populate_x_band()
    xMcE_i = data.xMc * (1 - (data.xVc - E_i) / data.xEg)
    xMcE_j = data.xMc * (1 - (data.xVc - E_j) / data.xEg)
    xMcE_j_avg = 0.5 * (xMcE_j[0:-1]+xMcE_j[1:])
    psi_i_avg = 0.5 * (psi_i[0:-1]+psi_i[1:])
    z = sum(psi_i_avg * np.diff(psi_j/xMcE_i) 
            + 1/xMcE_j_avg * (psi_i_avg * np.diff(psi_j)))
    z *= hbar**2/(2*(E_i-E_j)*e0) * 1e10/m0
    return z

def coupling_energy(data, dCL, upper, lower):
    #here, psi_i is the left-most wavefunction, not the wf with the highest energy
    module_i = data.moduleID[upper]
    module_j = data.moduleID[lower]
    psi_i = data.xyPsi[:,upper]
    psi_j = data.xyPsi[:,lower]
    if module_i > module_j:
        module_i = data.moduleID[lower]
        module_j = data.moduleID[upper]
        psi_i = data.xyPsi[:,lower]
        psi_j = data.xyPsi[:,upper]

    DeltaV = np.ones(data.xPointsPost.size)
    first = int(dCL[int(module_i)].widthOffset/data.xres)
    last = first + dCL[int(module_i)].xBarriers[int(100/data.xres):].size
    #  print "---debug--- coupling_energy"
    #  print first,last
    DeltaV[first:last] = dCL[int(module_i)].xBarriers[int(100/data.xres):]
    DeltaV = 1 - DeltaV
    couplingEnergy = sum(psi_i * DeltaV * psi_j) * data.xres*1e-10 \
            * abs(data.EcG[1] - data.EcG[0]) * 1000 #* (1-data.xBarriers)
    return couplingEnergy

def broadening_energy(data, upper, lower):
    if upper < lower:
        temp = upper
        upper = lower
        lower = temp
    psi_i = data.xyPsi[:,upper]
    psi_j = data.xyPsi[:,lower]
    
    transitions = np.bitwise_xor(data.xBarriers[0:-1].astype(bool),
            data.xBarriers[1:].astype(bool))
    transitionIdxs = np.nonzero(transitions == True)[0]# + int(100/data.xres)
    psi2int2 = sum((psi_i[transitionIdxs]**2-psi_j[transitionIdxs]**2)**2)
    DeltaLambda = 0.76 * 1e-9 * 1e-9
    twogamma = pi*data.me[0]*m0*e0**2/hbar**2 * DeltaLambda**2 \
            * (data.EcG[1] - data.EcG[0])**2 * psi2int2
    twogamma *=  1e3/e0 #convert to meV
    return twogamma

def alphaISB(data, stateR, lower):
    statesQ = []
    dipoles = []
    gammas = []
    energies = []
    for q in xrange(stateR+1, data.EigenE.size):
        dp = dipole(data, q, stateR)
        if abs(dp) > 1e-6:
            statesQ.append(q)
            dipoles.append(dp)
    for q in statesQ:
        gamma = broadening_energy(data, q, stateR)/2
        gammas.append(gamma)
        energies.append(data.EigenE[q]-data.EigenE[stateR])
        
    dipoles = np.array(dipoles)*1e-10 #in m
    gammas = np.array(gammas)/1000 #from meV to eV
    energies = abs(np.array(energies)) #in eV
    
    neff = 3
    Lp = sum(data.layerWidths[1:]) * 1e-10 #in m
    Nq = sum(data.layerDopings[1:]*data.layerWidths[1:]) / sum(data.layerWidths[1:])
    Nq *= 100**3 #convert from cm^-3 to m^-3
    Ns = sum(data.layerDopings[1:]*1e17 * data.layerWidths[1:]*1e-8) #in cm^-2
    Ns *= 100**2 #from cm^-2 to m^-2
    hw = data.EigenE[stateR] - data.EigenE[lower]
    
#    hw = np.arange(0.15,0.5,0.01)
#    for enerG in hw:
#        alphaISB = sum(energies * dipoles**2 * gammas / ((energies - enerG)**2 + gammas**2))

    h = 4.1356675e-15 #eV s
    #eps0 = 5.527e7 #1/eV-m
    eps0 = 8.854e-12 #C/V-m
    #print energies/h/c0 * dipoles**2
    #print gammas / ((energies - hw)**2 + gammas**2)
    
    alphaISB = sum(energies/h/c0 * dipoles**2 * gammas 
            / ((energies - hw)**2 + gammas**2))
    alphaISB *= 4*pi*e0**2 / (eps0*neff) * pi/(2*Lp) * Ns
    alphaISB /= e0*100
    
    if False: #plot loss diagram
        hw = np.arange(0.15,0.5,0.001)
        alphaISBw = np.zeros(hw.size)
        for q, enerG in enumerate(hw):
            alphaISBw[q] = sum(energies/h/c0 * dipoles**2 * gammas 
                    / ((energies - enerG)**2 + gammas**2))
        alphaISBw *= 4*pi*e0**2 / (eps0*neff) * pi/(2*Lp) * Ns
        alphaISBw /= e0*100
        plt.figure()
        plt.plot(hw,alphaISBw)
        plt.show()
        
    return alphaISB

def reflectivity(beta):
    beta = beta.real
    return ((beta - 1) / (beta + 1))**2



if __name__  == "__main__":
    print 'hi', cFunctions.returnme()

    
