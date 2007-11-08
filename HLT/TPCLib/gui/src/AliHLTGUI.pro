TEMPLATE	= app
LANGUAGE	= C++

CONFIG	+= qt warn_on release thread debug

HEADERS	+= AliHLTGUI.h

SOURCES	+= main.cxx \
	   AliHLTGUI.cxx

FORMS	=  AliHLTGUIMainForm.ui

IMAGES	=  AliceLogo.gif \
           KIPLogo.gif

include (flags)

#/************************************************************************
#**
#** This file is property of and copyright by the Computer Science/Computer 
#** Engineering Group, Kirchhoff Institute for Physics, Ruprecht-Karls-
#** University, Heidelberg, Germany, 2005
#** This file has been written by Jochen Thaeder, 
#** thaeder@kip.uni-heidelberg.de
#**
#**
#** See the file license.txt for details regarding usage, modification,
#** distribution and warranty.
#** Important: This file is provided without any warranty, including
#** fitness for any particular purpose.
#**
#**
#** Newer versions of this file's package will be made available from 
#** http://www.kip.uni-heidelberg.de/ti/HLT/
#** or the corresponding page of the Heidelberg Alice HLT group.
#**
#*************************************************************************/

TARGET   = AliHLTGUI









includeFile = $(QTROOTSYSDIR)/include/rootcint.pri
exists ($$includeFile) {
  include ($$includeFile)
}
!exists ($$includeFile) {
  includeFile = $$ROOTINCDIR/rootcint.pri
  exists ($$includeFile) {
    include ($$includeFile)
  }
  !exists ($$includeFile) {
    message (" ")
    message ("WARNING:  The $$inlcudeFile was not found !!!")
    message ("Please update your Qt layer version from http://root.bnl.gov ")
    message (" ")
  }
}

# Loading HLT libraries 
LIBS += $$ALIHLT_LDFLAGS $$ALIHLT_LIBS

# Loading ALICE TPC libraries 
LIBS += $$ALIROOT_LDFLAGS $$ALIROOT_LIBS

# Loading ROOT libraries
LIBS += -L$$ROOTLIBDIR

unix {
#
#  generate the link to the proper version of the ROOT resource file
#
  rootrc.target   = .rootrc
  ROOTRESOURCEFILE=rootrcqtgui
  !exists ($$ROOTLIBDIR/libQtGui.$$QMAKE_EXTENSION_SHLIB) {
      message ("No ROOT Qt Extension was found. Use Qt-layer instead")
      ROOTRESOURCEFILE = rootrcqt
  }
  rootrc.commands = @rm -rf .rootrc; ln -s $$ROOTRESOURCEFILE $$rootrc.target

  QMAKE_EXTRA_UNIX_TARGETS += rootrc
  PRE_TARGETDEPS  += $$rootrc.target
  QMAKE_CLEAN     += $$rootrc.target
}

unix {
  ROOTCINT_DIR = .rootcint
  UI_DIR = .ui
  MOC_DIR = .moc
  OBJECTS_DIR = .obj
}
