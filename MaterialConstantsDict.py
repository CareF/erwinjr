#!/usr/bin/env python
# -*- coding:utf-8 -*-

#===============================================================================
# ErwinJr is a simulation program for quantum semiconductor lasers.
# Copyright (C) 2012 Kale J. Franz, PhD
# Copyright (C) 2017 Ming Lyu, converting it to new format
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

#===============================================================================
# Reference
# [1]Handbook of Optics, Vol.2, ISBN: 0070479747
# [2]Herve P J L, Vandamme L K J. Empirical temperature dependence of the 
#    refractive index of semiconductors[J]. Journal of Applied Physics, 
#    1995, 77(10): 5476-5477.
#===============================================================================
class Material(object):
    def __init__(self, Name, Temperature=300):
        self.Name = Name
        self.set_temperature(Temperature)
        
    def set_temperature(self, Temperature):
        self.Temperature = Temperature

class MaterialConstantsDict(dict):
    def __init__(self, Temperature=300):
        self.set_constants()
        self.set_temperature(Temperature)
        
    def set_temperature(self, Temperature):
        self.Temperature = Temperature
        self["GaAs"].alc = 5.65325 + 3.88e-5*(Temperature-300) #Angs at 80 K
        self['InAs'].alc = 6.0583 + 2.74e-5*(Temperature-300)
        self['AlAs'].alc = 5.6611 + 2.90e-5*(Temperature-300)
        self['AlSb'].alc = 6.1355 + 2.60e-5*(Temperature-300)
        self['GaSb'].alc = 6.0959 + 4.72e-5*(Temperature-300)
        self['InSb'].alc = 6.4794 + 3.48e-5*(Temperature-300)
        self['InP'].alc = 5.8697+2.79e-5*(Temperature-300)
        
    def set_constants(self):
        
        # GaAs constants
        # from Vurgaftman
        self['GaAs'] = Material("GaAs")
        self['GaAs'].c11 = 1221 #GPa
        self['GaAs'].c12 = 566  #GPa
        self['GaAs'].EgG = 1.519 #eV at 0 K 
        self['GaAs'].EgL = 1.815 #eV
        self['GaAs'].EgX = 1.981 #eV
        self['GaAs'].VBO = -0.80 #eV
        self['GaAs'].DSO = 0.341 #eV Delta_SO from top of valence band
        self['GaAs'].acG = -7.17 #eV
        self['GaAs'].acL = -4.91 #eV, NextNano DB
        self['GaAs'].acX = -0.16 #eV, NaxtNano DB
        self['GaAs'].av  = -1.16 #eV, valence band deformation potential
        self['GaAs'].b   = -2.0  #eV, Pikus Burr Uniaxial Strain deformation potential
        self['GaAs'].XiG = 0 #eV, uniaxial deformation potential, NextNano DB
        self['GaAs'].XiL = 6.5 #eV
        self['GaAs'].XiX = 14.26 #eV
        self['GaAs'].me0 = 0.067 #1/m0: effective mass
        self['GaAs'].Ep  = 28.8 #eV
        self['GaAs'].F   = -1.94
        self['GaAs'].alG = 0.5405e-3 #eV/K, Varshni alpha(Gamma)
        self['GaAs'].beG = 204 #K, Varshni beta(Gamma)
        self['GaAs'].alX = 0.460e-3
        self['GaAs'].beX = 204
        self['GaAs'].alL = 0.605e-3
        self['GaAs'].beL = 204
        self['GaAs'].epss = 12.9
        self['GaAs'].epsInf = 10.86
        self['GaAs'].hwLO = 35.3*1e-3
        # Handbook of Optics, 2nd edition, Vol. 2. McGraw-Hill 1994
        # Table 22 Room-temperature Dispersion Formulas for Crystals
        # Sellmeier dispersion formula: 
        # n^2 = permitivity = 2 Lorenzians
        #     = c1 + c2 * wl**2/(wl**2-c3**2) + c4 * wl**2/(wl**2-c5**2) 
        # Temperature deps ~5E-5 [2] which can be ignored
        # TODO: material dependent formula should be here
        self['GaAs'].C1 = 3.5 
        self['GaAs'].C2 = 7.4969
        self['GaAs'].C3 = 0.4082
        self['GaAs'].C4 = 1.9347
        self['GaAs'].C5 = 37.17
        
        # InAs constants
        # from Vurgaftman
        self['InAs'] = Material("InAs")
        self['InAs'].c11 = 832.9
        self['InAs'].c12 = 452.6
        self['InAs'].EgG = 0.417
        self['InAs'].EgL = 1.133
        self['InAs'].EgX = 1.433
        self['InAs'].VBO = -0.59
        self['InAs'].DSO = 0.39
        self['InAs'].acG = -5.08
        self['InAs'].acL = -3.89 #eV, NextNano DB
        self['InAs'].acX = -0.08 #eV, NextNano DB
        self['InAs'].av  = -1.00
        self['InAs'].b   = -1.8
        self['InAs'].XiG = 0
        self['InAs'].XiL = 11.35
        self['InAs'].XiX = 3.7
        self['InAs'].me0 = 0.026
        self['InAs'].Ep  = 21.5
        self['InAs'].F   = -2.9
        self['InAs'].alG = 0.276e-3 
        self['InAs'].beG = 93 
        self['InAs'].alX = 0.276e-3
        self['InAs'].beX = 93
        self['InAs'].alL = 0.276e-3
        self['InAs'].beL = 93
        self['InAs'].epss = 14.3
        self['InAs'].epsInf = 11.6
        self['InAs'].hwLO = 29.93*1e-3
        self['InAs'].C1 = 11.1 #Handbook of Optics, 2nd edition, Vol. 2. McGraw-Hill 1994
        self['InAs'].C2 = 0.71
        self['InAs'].C3 = 2.551
        self['InAs'].C4 = 2.75
        self['InAs'].C5 = 45.66
        
        # AlAs constants
        # from Vurgaftman
        self['AlAs'] = Material('AlAs')
        self['AlAs'].c11 = 1250
        self['AlAs'].c12 = 534
        self['AlAs'].EgG = 3.099
        self['AlAs'].EgL = 2.46
        self['AlAs'].EgX = 2.24
        self['AlAs'].VBO = -1.33
        self['AlAs'].DSO = 0.28
        self['AlAs'].acG = -5.64
        self['AlAs'].acL = -3.07 #NextNano DB
        self['AlAs'].acX = 2.54  #NextNano DB
        self['AlAs'].av  = -2.47
        self['AlAs'].b   = -2.3
        self['AlAs'].XiG = 0
        self['AlAs'].XiL = 11.35
        self['AlAs'].XiX = 6.11
        self['AlAs'].me0 = 0.15
        self['AlAs'].Ep  = 21.1
        self['AlAs'].F   = -0.48
        self['AlAs'].alG = 0.855e-3
        self['AlAs'].beG = 530
        self['AlAs'].alX = 0.70e-3
        self['AlAs'].beX = 530
        self['AlAs'].alL = 0.605e-3
        self['AlAs'].beL = 204
        self['AlAs'].epss = 10.06
        self['AlAs'].epsInf = 8.16
        self['AlAs'].hwLO = 49.8*1e-3
        self['AlAs'].C1 = 2.0792 #Handbook of Optics, 2nd edition, Vol. 2. McGraw-Hill 1994
        self['AlAs'].C2 = 6.0840
        self['AlAs'].C3 = 0.2822
        self['AlAs'].C4 = 1.900
        self['AlAs'].C5 = 27.62
        
        # AlSb constants
        # from Vurgaftman
        self['AlSb'] = Material('AlSb')
        self['AlSb'].c11 = 876.9
        self['AlSb'].c12 = 434.1
        self['AlSb'].EgG = 2.386
        self['AlSb'].EgL = 2.329
        self['AlSb'].EgX = 1.696
        self['AlSb'].VBO = -0.41
        self['AlSb'].DSO = 0.676
        self['AlSb'].acG = -4.5
        self['AlSb'].acL = 0 #NextNano DB
        self['AlSb'].acX = 2.54  #NextNano DB
        self['AlSb'].av  = -1.4
        self['AlSb'].b   = -1.35
        self['AlSb'].XiG = 0
        self['AlSb'].XiL = 11.35
        self['AlSb'].XiX = 6.11
        self['AlSb'].me0 = 0.14
        self['AlSb'].Ep  = 18.7
        self['AlSb'].F   = -0.56
        self['AlSb'].alG = 0.42e-3
        self['AlSb'].beG = 140
        self['AlSb'].alX = 0.39e-3
        self['AlSb'].beX = 140
        self['AlSb'].alL = 0.58e-3
        self['AlSb'].beL = 140
        self['AlSb'].epss = 12.04 #ISBN 0849389127
        self['AlSb'].epsInf = 10.24 #ISBN 0849389127
        self['AlSb'].hwLO = 42.7 #http://prb.aps.org/pdf/PRB/v43/i9/p7231_1
        
        # GaSb constants
        # from Vurgaftman
        self['GaSb'] = Material('GaSb')
        self['GaSb'].c11 = 884.2
        self['GaSb'].c12 = 402.6
        self['GaSb'].EgG = 0.812
        self['GaSb'].EgL = 0.875
        self['GaSb'].EgX = 1.141
        self['GaSb'].VBO = -0.03
        self['GaSb'].DSO = 0.76
        self['GaSb'].acG = -7.5
        self['GaSb'].acL = 0 #unknown
        self['GaSb'].acX = 0 #unknown
        self['GaSb'].av  = -0.8
        self['GaSb'].b   = -2.0
        self['GaSb'].XiG = 0 #unknown
        self['GaSb'].XiL = 0 #unknown
        self['GaSb'].XiX = 0 #unknown
        self['GaSb'].me0 = 0.039
        self['GaSb'].Ep  = 27.0
        self['GaSb'].F   = -1.63
        self['GaSb'].alG = 0.417e-3
        self['GaSb'].beG = 140
        self['GaSb'].alX = 0.475e-3
        self['GaSb'].beX = 94
        self['GaSb'].alL = 0.597e-3
        self['GaSb'].beL = 140
        self['GaSb'].epss = 0 #unknown
        self['GaSb'].epsInf = 0 #unknown
        self['GaSb'].hwLO = 0 #unknown
     
        # InSb constants
        # from Vurgaftman
        self['InSb'] = Material('InSb')
        self['InSb'].c11 = 684.7
        self['InSb'].c12 = 373.5
        self['InSb'].EgG = 0.235
        self['InSb'].EgL = 0.93
        self['InSb'].EgX = 0.63
        self['InSb'].VBO = 0
        self['InSb'].DSO = 0.81
        self['InSb'].acG = -6.94
        self['InSb'].acL = 0 #unknown
        self['InSb'].acX = 0 #unknown
        self['InSb'].av  = -0.36
        self['InSb'].b   = -2.0
        self['InSb'].XiG = 0 #unknown
        self['InSb'].XiL = 0 #unknown
        self['InSb'].XiX = 0 #unknown
        self['InSb'].me0 = 0.0135
        self['InSb'].Ep  = 23.3
        self['InSb'].F   = -0.23
        self['InSb'].alG = 0.32e-3
        self['InSb'].beG = 170
        self['InSb'].alX = 0e-3 #unknown
        self['InSb'].beX = 0 #unknown
        self['InSb'].alL = 0e-3 #unknown
        self['InSb'].beL = 0 #unknown
        self['InSb'].epss = 0 #unknown
        self['InSb'].epsInf = 0 #unknown
        self['InSb'].hwLO = 0 #unknown
     
        # InP constants
        self['InP'] = Material('InP')
        self['InP'].me0 = 0.0795
        self['InP'].C1 = 7.255 #Handbook of Optics, 2nd edition, Vol. 2. McGraw-Hill 1994
        self['InP'].C2 = 2.316
        self['InP'].C3 = 0.6263
        self['InP'].C4 = 2.765
        self['InP'].C5 = 32.935
        
        # TODO: consider transform bowing parameter to function
        # InGaAs constants
        self['InGaAs'] = Material('InGaAs')
        self['InGaAs'].EgG = 0.477
        self['InGaAs'].EgL = 0.33
        self['InGaAs'].EgX = 1.4
        self['InGaAs'].VBO = -0.38
        self['InGaAs'].DSO = 0.15
        self['InGaAs'].acG = 2.61
        self['InGaAs'].acL = 2.61 #NextNano DB
        self['InGaAs'].acX = 2.61 #NextNano DB
        self['InGaAs'].me0 = 0.0091
        self['InGaAs'].Ep  = -1.48
        self['InGaAs'].F   = 1.77
        
        # AlInAs constants
        self['AlInAs'] = Material('AlInAs')
        self['AlInAs'].EgG = 0.70
        self['AlInAs'].EgL = 0
        self['AlInAs'].EgX = 0
        self['AlInAs'].VBO = -0.64
        self['AlInAs'].DSO = 0.15
        self['AlInAs'].acG = -1.4
        self['AlInAs'].acL = -1.4 #NextNano DB
        self['AlInAs'].acX = -1.4 #NextNano DB
        self['AlInAs'].me0 = 0.049
        self['AlInAs'].Ep  = -4.81
        self['AlInAs'].F   = -4.44
        
        # AlGaAs constants
        self['AlGaAs'] = Material('AlGaAs')
        self['AlGaAs'].EgG = -0.127 #+ 1.310*x(Al)
        # To describe the band gap bowing at the Gamma point in AlxGa1-xAs,
        # a linear or even a quadratic interpolation is not sufficient. 
        # Here, the bowing parameter is not constant, 
        # it depends on the alloy composition x.
        self['AlGaAs'].EgL = 0.055
        self['AlGaAs'].EgX = 0
        self['AlGaAs'].VBO = 0
        self['AlGaAs'].DSO = 0
        self['AlGaAs'].acG = 0
        self['AlGaAs'].acL = 0
        self['AlGaAs'].acX = 0
        self['AlGaAs'].me0 = 0
        self['AlGaAs'].Ep  = 0
        self['AlGaAs'].F   = 0
        
        # AlAsSb constants
        self['AlAsSb'] = Material('AlGaAs')
        self['AlAsSb'].EgG = 0.8
        self['AlAsSb'].EgL = 0.28
        self['AlAsSb'].EgX = 0.28
        self['AlAsSb'].DSO = 0.15
        self['AlAsSb'].VBO = -1.71
        
        # AlGaSb constants
        self['AlGaSb'] = Material('AlGaAs')
        self['AlGaSb'].EgG = -0.044 #+ 1.22*x(Al), same as AlGaAs
        self['AlGaSb'].EgL = 0
        self['AlGaSb'].EgX = 0
        self['AlGaSb'].VBO = 0
        self['AlGaSb'].DSO = 0.2
        self['AlGaSb'].acG = 0
        self['AlGaSb'].acL = 0 #NextNano DB
        self['AlGaSb'].acX = 0 #NextNano DB
        self['AlGaSb'].me0 = 0
        self['AlGaSb'].Ep  = 0
        self['AlGaSb'].F   = 0
        
        # InAsSb constants
        self['InAsSb'] = Material('InAsSb')
        self['InAsSb'].EgG = 0.67
        self['InAsSb'].EgL = 0.6
        self['InAsSb'].EgX = 0.6
        self['InAsSb'].VBO = 0
        self['InAsSb'].DSO = 1.2
        self['InAsSb'].acG = 0
        self['InAsSb'].acL = 0 #NextNano DB
        self['InAsSb'].acX = 0 #NextNano DB
        self['InAsSb'].me0 = 0.035
        self['InAsSb'].Ep  = 0
        self['InAsSb'].F   = 0

# vim: ts=4 sw=4 sts=4 expandtab
