#!/usr/bin/env python2
# -*- coding:utf-8 -*-
from __future__ import division
import numpy as np
from numpy import pi, sqrt
import ThePhysics
from ThePhysics import h, c0, e0
import sys
import cProfile

newLineChar = '\n'
if sys.platform == 'darwin':
    newLineChar = '\n'

def qclLoad(fname):
    qclayers = ThePhysics.QCLayers()
    #  print "Loading "+fname
    valDict = {}
    filehandle = open(fname, 'r')
    filehandle.readline() #throw away 'ErwinJr Data File'
    while True:
        line = filehandle.readline()
        if line == '# QC layers #'+newLineChar:
            break
        line = line.split(':')
        valDict[line[0]] = line[1].strip()
    
    qclayers.description = valDict['Description']
    qclayers.substrate   = valDict['Substrate']
    qclayers.EField      = float(valDict['Efield'])
    qclayers.xres        = float(valDict['xres'])
    qclayers.vertRes     = float(valDict['Eres'])
    qclayers.moleFrac[0]   = float(valDict['moleFrac1'])
    qclayers.moleFrac[1]   = float(valDict['moleFrac2'])
    qclayers.moleFrac[2]   = float(valDict['moleFrac3'])
    qclayers.moleFrac[3]   = float(valDict['moleFrac4'])
    qclayers.moleFrac[4]   = float(valDict['moleFrac5'])
    qclayers.moleFrac[5]   = float(valDict['moleFrac6'])
    qclayers.moleFrac[6]   = float(valDict['moleFrac7'])
    qclayers.moleFrac[7]   = float(valDict['moleFrac8'])
    qclayers.solver      = valDict['Solver']
    qclayers.Temperature = float(valDict['Temp'])
    qclayers.TempFoM     = float(valDict['TempFoM'])
    qclayers.repeats     = int(valDict['PlotPeriods'])
    qclayers.diffLength  = float(valDict['DiffLeng'])
    
    #  self.strata.wavelength           = float(valDict['Wavelength'])
    #  self.strata.operatingField       = float(valDict['StratumField'])
    #  self.strata.Lp                   = float(valDict['Lp'])
    #  self.strata.Np                   = float(valDict['Np'])
    #  self.strata.aCore                = float(valDict['alphaCore'])
    #  self.strata.nCore                = complex(valDict['nCore'])
    #  self.strata.nD                   = float(valDict['nD'])
    #  self.strata.transitionBroadening = float(valDict['transitionBroadening'])        
    #  self.strata.tauUpper             = float(valDict['tauUpper'])
    #  self.strata.tauLower             = float(valDict['tauLower'])
    #  self.strata.tauUpperLower        = float(valDict['tauUpperLower'])
    #  self.strata.opticalDipole        = float(valDict['opticalDipole'])
    #  self.strata.FoM                  = float(valDict['FoM'])
    #  self.strata.waveguideFacets      = valDict['waveguideFacets']
    #  self.strata.waveguideLength      = float(valDict['waveguideLength'])
    #  self.strata.customFacet          = float(valDict['customFacet'])
    
    lines = []
    while True:
        line = filehandle.readline()
        if line == '# Optical strata #'+newLineChar:
            break
        lines.append(line)
    rows = len(lines)
    variables = ['layerWidths', 'layerBarriers', 'layerARs', 'layerDopings', 
                 'layerMaterials', 'layerDividers']
    for item in variables:
        setattr(qclayers, item, np.zeros(rows))
    for q, line in enumerate(lines):
        line = line.split('\t')
        qclayers.layerWidths[q]    = float(line[1])
        qclayers.layerBarriers[q]  = float(line[2])
        qclayers.layerARs[q]       = float(line[3])
        qclayers.layerMaterials[q] = float(line[4])
        qclayers.layerDopings[q]   = float(line[5])
        qclayers.layerDividers[q]  = float(line[6])
    
    #  lines = filehandle.readlines()
    #  rows = len(lines)
    #  variables = ['stratumCompositions', 'stratumThicknesses', 'stratumDopings']
    #  for item in variables:
        #  setattr(self.strata, item, np.zeros(rows))
    #  self.strata.stratumMaterials = []
    #  for q, line in enumerate(lines):
        #  line = line.split('\t')
        #  self.strata.stratumMaterials.append(str(line[1]))
        #  self.strata.stratumCompositions[q] = float(line[2])
        #  self.strata.stratumThicknesses[q] = float(line[3])
        #  self.strata.stratumDopings[q]     = float(line[4])
    
    filehandle.close()
    qclayers.update_alloys()
    qclayers.update_strain()
    qclayers.populate_x()
    qclayers.populate_x_band()
    return qclayers

def check_class(a, b):
    for key, item in a.__dict__.items():
        if isinstance(item, np.ndarray):
            equal = (item == getattr(b, key)).all()
        else:
            equal = (item == getattr(b,key))
        print key, equal

if __name__  == "__main__":
    if len(sys.argv) != 3:
        print "Usage: python2 %s input_filename output_name"%sys.argv[0]
    f_input = sys.argv[1]
    f_output = sys.argv[2]
    qclayers = qclLoad(f_input)
    qclayers.solve_psi()
    upper = 19
    lower = 15
    FoM = qclayers.figure_of_merit(upper, lower)
    print qclayers.EigenE
    print FoM
    
# vim: ts=4 sw=4 sts=4 expandtab
