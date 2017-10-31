#!/usr/bin/env python2
# -*- coding:utf-8 -*-
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
    filehandle = open(fname, 'r')
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
    if len(sys.argv) == 3: 
        f_output = sys.argv[2]
    else: 
        import time
        f_output = 'PerformanceTest_'+time.strftime("%m-%d-%y-%H:%M:%S",
                time.localtime())
    qclayers = qclLoad(f_input)
    cProfile.run('main(qclayers)', filename = f_output)
    p = pstats.Stats(f_output)
    p.strip_dirs().sort_stats('tottime').print_stats(5)
    
# vim: ts=4 sw=4 sts=4 expandtab
