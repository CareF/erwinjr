ABOUT
=====

ErwinJr is an open source design and simulation program for quantum semiconductor devices including quantum cascade lasers. It is useful for modeling intersubband semiconductor devices.  The code base was written at the NASA Jet Propulsion Laboratory.

ErwinJr is a multi-platform application that runs on most desktop operating systems including Windows, Mac OS X, and Linux.  This is because it is written in Python and uses the Qt graphical user interface framework.


== WINDOWS INSTALLATION ==

1) Do a FULL install of pythonxy from http://www.pythonxy.com/

2) Under Start|Run type cmd

3) type the following at the command prompt
   cd \erwinjr
   gcc -c -fPIC cFunctions.c
   gcc -shared -fPIC -o cFunctions.so cFunctions.o
   
4) open (double click) on the file erwinjr.pyw


== Python Dependences ==

PyQt4, PyQwt5, numpy, scipy, matplotlib
psyco is also useful

Built and tested on Python 2.6.

== For Arch Linux User ==

community/python2-pyqwt package on the official source has some compiling bug, which will result in `Segmentation fault (core dumped)' when import PyQt4.Qwt5. See [bug report](https://bugs.archlinux.org/task/53918?project=5&cat%5B0%5D=33&string=python2-pyqwt) 

This can be fixed by building the package locally, using following commands: 
	asp export pacman community/python2-pyqwt
	cd python2-pyqwt
	makepkg
	sudo pacman -U pyton2-pyqwt-5.2.0-2-x86_64.pkg.tar.xz
