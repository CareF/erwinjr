#!/usr/bin/env python2
# -*- coding:utf-8 -*-

#===============================================================================
# ErwinJr is a simulation program for quantum semiconductor lasers.
# Copyright (C) 2012 Kale J. Franz, PhD
#               2017 Ming Lyu
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

from __future__ import division
import numpy as np
import sys
from QCLayers import QCLayers
from Strata import Strata
newLineChar = '\n'
if sys.platform == 'darwin':
    newLineChar = '\n'

numMaterials = 8
def qclLoad(filehandle, qclayers=None, strata=None):
    """Load file with filename 'fname' into qclayers and strata"""
    if not isinstance(qclayers, QCLayers) and not isinstance(strata, Strata):
        raise TypeError("qclLoad: Nothing to load.." 
                "Both QCLayers and Strata are not valid type")
    valDict = {}
    filehandle.readline() #throw away 'ErwinJr Data File'
    while True:
        line = filehandle.readline()
        if line == '# QC layers #'+newLineChar:
            break
        line = line.split(':')
        valDict[line[0]] = line[1].strip()
    
    if isinstance(qclayers, QCLayers):
        qclayers.description = valDict['Description']
        qclayers.substrate   = valDict['Substrate']
        qclayers.EField      = float(valDict['Efield'])
        qclayers.xres        = float(valDict['xres'])
        qclayers.vertRes     = float(valDict['Eres'])
        for n in range(numMaterials):
            qclayers.moleFrac[n] = float(valDict['moleFrac%d'%(n+1)])
        qclayers.solver      = valDict['Solver']
        qclayers.Temperature = float(valDict['Temp'])
        qclayers.TempFoM     = float(valDict['TempFoM'])
        qclayers.repeats     = int(valDict['PlotPeriods'])
        qclayers.diffLength  = float(valDict['DiffLeng'])

        # QC layers
        lines = []
        while True:
            line = filehandle.readline()
            if line == '# Optical strata #'+newLineChar:
                break
            lines.append(line)
        rows = len(lines)
        qclayers.layerWidth = np.empty(rows, np.int_)
        for item in ('layerBarriers', 'layerARs', 
                'layerDopings', 'layerMaterials', 'layerDividers'):
            setattr(qclayers, item, np.zeros(rows))
        for q, line in enumerate(lines):
            line = line.split('\t')
            qclayers.layerWidth[q]     = int(np.round(float(line[1])
                                                  /qclayers.xres))
            qclayers.layerBarriers[q]  = float(line[2])
            qclayers.layerARs[q]       = float(line[3])
            qclayers.layerMaterials[q] = float(line[4])
            qclayers.layerDopings[q]   = float(line[5])
            qclayers.layerDividers[q]  = float(line[6])

    if isinstance(strata, Strata):
        strata.wavelength           = float(valDict['Wavelength'])
        strata.operatingField       = float(valDict['StratumField'])
        strata.Lp                   = float(valDict['Lp'])
        strata.Np                   = float(valDict['Np'])
        strata.aCore                = float(valDict['alphaCore'])
        strata.nCore                = complex(valDict['nCore'])
        strata.nD                   = float(valDict['nD'])
        strata.transitionBroadening = float(valDict['transitionBroadening'])        
        strata.tauUpper             = float(valDict['tauUpper'])
        strata.tauLower             = float(valDict['tauLower'])
        strata.tauUpperLower        = float(valDict['tauUpperLower'])
        strata.opticalDipole        = float(valDict['opticalDipole'])
        strata.FoM                  = float(valDict['FoM'])
        strata.waveguideFacets      = valDict['waveguideFacets']
        strata.waveguideLength      = float(valDict['waveguideLength'])
        strata.customFacet          = float(valDict['customFacet'])

        # Optical strata
        lines = filehandle.readlines()
        rows = len(lines)
        for item in ('stratumCompositions', 'stratumThicknesses',
                'stratumDopings'):
            setattr(strata, item, np.zeros(rows))
        strata.stratumMaterials = []
        for q, line in enumerate(lines):
            line = line.split('\t')
            strata.stratumMaterials.append(str(line[1]))
            strata.stratumCompositions[q] = float(line[2])
            strata.stratumThicknesses[q] = float(line[3])
            strata.stratumDopings[q]     = float(line[4])

    return True

def qclSave(filehandle, qclayers, strata):
    """Save file with filename 'fname' from qclayers and strata"""
    if not isinstance(qclayers, QCLayers) and not isinstance(strata, Strata):
        raise TypeError("qclSave: Nothing to save.." 
                "Both QCLayers and Strata are not valid type")

    if isinstance(qclayers, QCLayers):
        filehandle.write("Description:" + qclayers.description + '\n')
        filehandle.write("Substrate:" + qclayers.substrate + '\n')
        filehandle.write("Efield:" + str(qclayers.EField) + '\n')
        filehandle.write("xres:" + str(qclayers.xres) + '\n')
        filehandle.write("Eres:" + str(qclayers.vertRes) + '\n')
        for n in range(numMaterials ):
            filehandle.write("moleFrac%d:"%(n+1) +
                    str(qclayers.moleFrac[n]) + '\n')
        filehandle.write("Solver:" + qclayers.solver + '\n')
        filehandle.write("Temp:" + str(qclayers.Temperature) + '\n')
        filehandle.write("TempFoM:" + str(qclayers.TempFoM) + '\n')
        filehandle.write("PlotPeriods:" + str(qclayers.repeats) + '\n')
        filehandle.write("DiffLeng:" + str(qclayers.diffLength) + '\n')

    if isinstance(strata, Strata):
        filehandle.write("Wavelength:" + str(strata.wavelength) + '\n')
        filehandle.write("StratumField:" + str(strata.operatingField) + '\n')
        filehandle.write("Lp:" + str(strata.Lp) + '\n')
        filehandle.write("Np:" + str(strata.Np) + '\n')
        filehandle.write("alphaCore:" + str(strata.aCore) + '\n')
        filehandle.write("nCore:" + str(strata.nCore) + '\n')
        filehandle.write("nD:" + str(strata.nD) + '\n')
        filehandle.write("transitionBroadening:" + 
                str(strata.transitionBroadening) + '\n')
        filehandle.write("tauUpper:" + str(strata.tauUpper) + '\n')
        filehandle.write("tauLower:" + str(strata.tauLower) + '\n')
        filehandle.write("tauUpperLower:" + str(strata.tauUpperLower) + '\n')
        filehandle.write("opticalDipole:" + str(strata.opticalDipole) + '\n')
        filehandle.write("FoM:" + str(strata.FoM) + '\n')
        filehandle.write("waveguideFacets:" + strata.waveguideFacets + '\n')
        filehandle.write("waveguideLength:" + str(strata.waveguideLength) + '\n')
        filehandle.write("customFacet:" + str(strata.customFacet) + '\n')

    filehandle.write("# QC layers #\n")
    if isinstance(qclayers, QCLayers):
        for row in xrange(qclayers.layerWidth.size):
            string = "%d\t%f\t%d\t%d\t%d\t%f\t%d\n" % (row+1, 
                    qclayers.xres * qclayers.layerWidth[row], 
                    qclayers.layerBarriers[row], 
                    qclayers.layerARs[row], 
                    qclayers.layerMaterials[row], 
                    qclayers.layerDopings[row], 
                    qclayers.layerDividers[row])
            filehandle.write(string)

    filehandle.write("# Optical strata #\n")
    if isinstance(strata, Strata):
        for row in xrange(strata.stratumDopings.size):
            string = "%d\t%s\t%f\t%f\t%f\n" % (row+1, 
                    strata.stratumMaterials[row], 
                    strata.stratumCompositions[row], 
                    strata.stratumThicknesses[row], 
                    strata.stratumDopings[row])
            filehandle.write(string)

    return True

# vim: ts=4 sw=4 sts=4 expandtab
