#!/usr/bin/env python2
# -*- coding:utf-8 -*-

#============================================================================
# ErwinJr is a simulation program for quantum semiconductor lasers.
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

# This is a command line performance test, and can be extended to a command
# line tool for ErwinJr

from __future__ import division
import numpy as np
from QCLayers import QCLayers
import sys
import cProfile
import SaveLoad

newLineChar = '\n'
if sys.platform == 'darwin':
    newLineChar = '\n'

def qclLoad(fname):
    qclayers = QCLayers()
    #  print "Loading "+fname
    filehandle = open(fname, 'rU')
    SaveLoad.qclLoad(filehandle, qclayers, None)
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

def main(qclayers):
    qclayers.solve_psi()
    upper = 19
    lower = 15
    FoM = qclayers.figure_of_merit(upper, lower)
    print FoM

if __name__  == "__main__":
    if not len(sys.argv) in (2,3):
        print "Usage: python2 %s input_filename [output_name]"%sys.argv[0]
    import cProfile, pstats
    f_input = sys.argv[1]
    qclayers = qclLoad(f_input)
    if len(sys.argv) == 3: 
        f_output = sys.argv[2]
    else: 
        import time
        f_output = 'PerformanceTest_'+time.strftime("%m-%d-%y-%H:%M:%S",
                time.localtime())
    cProfile.run('main(qclayers)', filename = f_output)
    p = pstats.Stats(f_output)
    p.strip_dirs().sort_stats('tottime').print_stats(5)
    if len(sys.argv) == 2: 
        import subprocess
        subprocess.Popen(["rm", f_output])
    
# vim: ts=4 sw=4 sts=4 expandtab
