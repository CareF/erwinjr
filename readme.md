ABOUT
=====

ErwinJr is an open source design and simulation program for quantum semiconductor devices including quantum cascade lasers. It is useful for modeling intersubband semiconductor devices.  The code base (Version before 3.0) was written at the NASA Jet Propulsion Laboratory. And now it's updating and maintained by minglyu@princeton. 

ErwinJr is a multi-platform application that runs on most desktop operating systems including Windows, Mac OS X, and Linux.  This is because it is written in Python and uses the Qt graphical user interface framework.

### WINDOWS INSTALLATION ###

(Update: GUI under Win is not tested for Ver>=3.0 because qwt support for new Windows installation is out of date.)

1) Do a FULL install of pythonxy from http://www.pythonxy.com/  (no longer recommended)

1*) (Recommended) if you are not going to use optical cavity design tab, try `noOptic` branch (which gets rid of qwt dependence), and it should be able to work under any Python distribution with pyqt4, numpy, scipy and matplotlib.. (And mock package is also needed for blocking optical tab). Anaconda is recommended. 

~~2) Under Start|Run type cmd~~

3) type the following at the command prompt (See following section about compiling using Visual Studio)

		cd \erwinjr
		make

4) open (double click) on the file erwinjr.pyw, or in command line, `python erwinjr.pyw`

### Multiprocessing Support ###

The multiprocessing feature (cQCLayersMP) requires openmp support, tested under Linux with openmp (ver>=5.0) and Windows with Visual Studio 2017. 

To compile for cQCLayersMP under *nix, using following command: 

	make cQCLayersMP.so

See Makefile for building detail. 

### Visual Studio support ###

Visual Studio support is tested under Windows 10 and Visual Studio 2017 community version. 

Choosing from `Solution Configurations` to compile for non-multiprocessing (`Release`) or multiprocessing (`ReleaseMP`), and press `Ctrl+Shift+B` to build. 

### Python Dependences ###

PyQt4, PyQwt5, numpy, scipy, matplotlib

~~psyco is also useful~~ (psyco is out of date)

Built and tested on Python 2.7.

### For Arch Linux User ###

*UPDATE*: Bug fixed. Compiling locally is no longer required. 

`community/python2-pyqwt` package on the official source has some compiling bug, which will result in `Segmentation fault (core dumped)` when import PyQt4.Qwt5. See [bug report](https://bugs.archlinux.org/task/53918?project=5&cat%5B0%5D=33&string=python2-pyqwt) 

This can be fixed by building the package locally, using following commands: 

	asp export pacman community/python2-pyqwt
	cd python2-pyqwt
	makepkg
	sudo pacman -U pyton2-pyqwt-5.2.0-2-x86_64.pkg.tar.xz

