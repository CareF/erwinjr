#!/usr/bin/env python
# -*- coding:utf-8 -*-
import numpy as np
import sys, os, re
from QCLayers import QCLayers
import SaveLoad
import json
qcMaterial = {"InP":  ["InGaAs", "AlInAs"], 
              "GaAs": ["AlGaAs", "AlGaAs"], 
              "GaSb": ["InAsSb", "AlGaSb"]
              }
class QCEncoder(json.JSONEncoder):
    def default(self, o):
        if not isinstance(o, QCLayers):
            raise TypeError("Input type is not QCLayers")
        usedMaterial = int(np.max(o.layerMaterials))
        materialList = (2 * (o.layerMaterials - 1) + o.layerBarriers)
        materialList = materialList.astype(int)
        layerWidth = o.layerWidth * o.xres
        output = {"FileType":    "ErwinJr2 Data File",
                  "Version":     "181106",
                  "Descrption":  o.description,
                  "Substrate":   o.substrate, 
                  "EField":      o.EField,
                  "x resolution":        o.xres,
                  "E resolution":        o.vertRes, 
                  "Solver":      o.solver, 
                  "Temperature": o.Temperature, 
                  "Repeats":     o.repeats,
                  "Materials": {"Compostion":    (qcMaterial[o.substrate] * 
                                                  usedMaterial), 
                                "Mole Fraction": o.moleFrac[:2*usedMaterial]
                                },
                  "QC Layers": {"Material":      materialList[1:].tolist(),
                                "Width":         layerWidth[1:].tolist(), 
                                "Doping":        o.layerDopings[1:].tolist(), 
                                "Active Region": o.layerARs[1:].tolist()
                                }}
        return output 

JSONTemplate = """{
    "FileType": "ErwinJr2 Data File", 
    "Version": "181107", 
    "Descrption": %s, 
    "Substrate": %s, 
    "EField": %s, 
    "x resolution": %s, 
    "E resolution": %s, 
    "Solver": %s, 
    "Temperature": %s, 
    "Repeats": %s, 
    "Materials": {
        "Compostion": %s, 
        "Mole Fraction": %s
    }, 
    "QC Layers": {
        "Material": %s, 
        "Width": %s, 
        "Doping": %s, 
        "Active Region": %s
    }
}"""

def qclSaveJSON(fhandle, qclayers):
    "Save file with filename as json"
    if not isinstance(qclayers, QCLayers):
        raise TypeError("qclSave: Nothing to save.."
                        "QCLayers not valid type")
    #  json.dump(qclayers, fhandle, cls=QCEncoder, indent=4)
    o=qclayers
    usedMaterial = int(np.max(o.layerMaterials))
    materialList = (2 * (o.layerMaterials - 1) + o.layerBarriers)
    materialList = materialList.astype(int)
    layerWidth = o.layerWidth * o.xres
    parameters = [json.encoder.encode_basestring(o.description)] + [
        repr(s).replace("'","\"") for s in (
            o.substrate, o.EField, o.xres, o.vertRes, o.solver, 
            o.Temperature, o.repeats, 
            qcMaterial[o.substrate] * usedMaterial, 
            o.moleFrac[:2*usedMaterial], materialList[1:].tolist(),
            layerWidth[1:].tolist(), o.layerDopings[1:].tolist(), 
            o.layerARs[1:].tolist())]
    fhandle.write(JSONTemplate%tuple(parameters))

if __name__ == "__main__":
    if len(sys.argv) == 3: 
        qclayers = QCLayers()
        with open(sys.argv[1], 'rU') as fin, open(sys.argv[2], 'w') as fout:
            SaveLoad.qclLoad(fin, qclayers)
            qclSaveJSON(fout, qclayers)
    else: 
        print "Call this converter by python %s <orignal> <new>"%__file__

# vim: ts=4 sw=4 sts=4 expandtab
