#setup script for py2exe
#run at command line:  python setup.py py2exe

from distutils.core import setup
import py2exe
import matplotlib
import os

filesList = matplotlib.get_py2exe_datafiles()
filesList.append(('', ['EJico.ico']))
filesList.append(('', ['tutorial.pdf']))
filesList.append(('', ['cQCLayers.so']))
filesList.append(('', ['cStrata.so']))
filesList.append(('', ['license.txt']))
filesList.append(('', ['EJico.ico']))
filesList.append(('src', ['cQCLayers.c']))
filesList.append(('src', ['cStrata.c']))
filesList.append(('src', ['complex.h']))
filesList.append(('src', ['erwinjr.pyw']))
filesList.append(('src', ['MaterialConstantsDict.py']))
filesList.append(('src', ['settings.py']))
filesList.append(('src', ['readme.txt']))
filesList.append(('src', ['create_exe.py']))
filesList.append(('src', ['SupportClasses.py']))
filesList.append(('src', ['QCLayers.py']))
filesList.append(('src', ['Strata.py']))
filesList.append(('src', ['QuantumTab.py']))
filesList.append(('src', ['setup_script.iss']))
dirPaths = ['images', 'examples']
for dirPath in dirPaths:
    for files in os.listdir(dirPath):
        f1 = dirPath + '/' + files
        if os.path.isfile(f1): # skip directories
            f2 = dirPath, [f1]
            filesList.append(f2)

setup(windows=[{"script": "erwinjr.pyw",
                "icon_resources": [(1,"EJico.ico")]}],
      data_files=filesList,
      options={"py2exe": {"skip_archive": True, "includes": ["sip","PyQt4.QtSvg","matplotlib.backends.backend_qt4agg"]}}
      )
#http://www.riverbankcomputing.com/pipermail/pyqt/2010-February/025806.html
