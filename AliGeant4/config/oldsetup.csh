# $Id$
# ----------------------------------------------------------------
# This script sets the default environment variables for
# Alice Geant4 based prototype
# Options: -g     verbose mode
#          local  recomended local configuration
#   
# csh version by I.Gonzalez  18.2.2000



#
# ==================================================================
# Alice configuration options: Please define your selection with the
# variables below.
# ==================================================================
# 

#
# ====== AG4_VERSION
# Geant4 version
# If set: the provided Geant4 version and not the default one is set
#setenv AG4_VERSION 2.0_test

#
# ====== AG4_VISUALIZE
# Set/Unset to get/avoid Geant4 visualisation.
setenv AG4_VISUALIZE 1
#unsetenv AG4_VISUALIZE

#
# ====== AG4_OPACS
# Set/Unset to get OPACS as a Geant4 visualisation driver.
#setenv AG4_OPACS 1
unsetenv AG4_OPACS

#
# ====== AG4_MAKESHLIB
# If set: shared libraris are created
setenv AG4_MAKESHLIB 1
#unsetenv AG4_MAKESHLIB


#
# Resolve input parameters
#
set VERBOSE = "NO"
set LOCAL = "NO"
set SILENT = "NO"
foreach param ( $* )
  switch ( $param )
    case -g:
      set VERBOSE="YES"; shift; breaksw;
    case local:
      set LOCAL="YES"; shift; breaksw;
    case silent:
      set SILENT="YES"; shift; breaksw;
  endsw
end

#
# ==================================================================
# Key path variables
# ==================================================================
#

if ( "$LOCAL" == "NO" ) then
  # 
  # AFS
  #

  # ====== ALICE_BASE
  # ALICE base directory
  set ALICE_BASE = /afs/cern.ch/alice/offline

  # ====== G4_BASE
  # Geant4 base directory
  set G4_BASE = ${ALICE_BASE}/geant4

  # ====== G4MC_BASE
  # Geant4_mc base directory
  set G4MC_BASE = ${ALICE_BASE}

  # ====== LHCXX_BASE
  # LHC++ base directory
  set LHCXX_BASE = /afs/cern.ch/sw/lhcxx/specific/@sys

  # ====== ROOT_BASE
  # Root base directory  
  set ROOT_BASE = /afs/cern.ch/alice/library/root
  # to be removed after aCC support will be defined
  # in the offline in a standard way    
  if ( "`uname`" == "HP-UX" ) then 
   set ROOT_BASE = /afs/cern.ch/alice/library/.hp_ux102/root.2.23.aCC
  endif 

  # ====== IRST_BASE
  # IRST code check tool base directory
  setenv IRST_BASE /afs/cern.ch/alice/offline/geant4/tools/IRST

else
  # 
  # recommended local installation
  #

  # ====== ALICE_BASE
  # ALICE base directory
  set ALICE_BASE = $ALICE_INSTALL/$ALICE_LEVEL

  # ====== G4_BASE
  # Geant4 base directory
  set G4_BASE = $ALICE_BASE

  # ====== G4MC_BASE
  # Geant4_mc base directory
  set G4MC_BASE = $ALICE_BASE

  # ====== LHCXX_BASE
  # LHC++ base directory
  set LHCXX_BASE = $ALICE_BASE

  # ====== ROOT_BASE
  # Root base directory  
  set ROOT_BASE = $ALICE_BASE

  # ====== IRST_BASE
  # IRST code check tool base directory
  #setenv IRST_BASE $HOME/dev/tools/IRST

endif


#....................................... SKIP ................................
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# You should not need to change 
# the lines below this one
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

# ====== AG4_INSTALL
# Alice Geant4 install directory
#
setenv AG4_INSTALL $ALICE_ROOT/AliGeant4

# Set variables depending on other variables
# if opacs is selected then select also visualize
if ( "$?AG4_OPACS" == 1 ) then
  setenv AG4_VISUALIZE 1
endif

# Start the output
if ( "$VERBOSE" == "YES" ) then
  echo " "
  echo "    =========================================="
  echo "      ALICE Geant4 environment configuration"
  echo "    =========================================="

  #
  # Some output on the selections... 
  #

  if ("$?AG4_VISUALIZE" == 1) then
    echo "Geant4 visualisation is     selected."
  else
    echo "Geant4 visualisation is NOT selected."
  endif
  if ("$?AG4_OPACS" == 1) then
    echo "OPACS driver         is     selected."
  else
    echo "OPACS driver         is NOT selected."
  endif
  if ("$?AG4_NOPHYSICS" == 1) then
    echo "Only geantino or charged geantino can be shooted."
  else
    echo "Full physics has been selected."
  endif
  if ("$?AG4_ODBMS" == 1) then
    echo "The environment for using Objectivity will be set."
  else
    echo "No Geant4 Persistency."
  endif

endif

#
# ROOT Alice definitions & options
# ==================================
#
if ( "$?ROOTSYS" == 0 ) then
  setenv ROOTSYS ${ROOT_BASE}
endif
if ( "$VERBOSE" == "YES" ) then
  echo " "
  echo "ROOT"
  echo "===="
  echo ROOTSYS set to $ROOTSYS
endif

#
# Geant4  
# ==================================
#
if ( "$VERBOSE" == "YES" ) then
  echo " "
  echo "Geant4 env. variables..."
  echo "============================"
endif

if ("$?AG4_VERSION" == 0) then
  setenv G4INSTALL ${G4_BASE}/geant4
else
  setenv G4INSTALL ${G4_BASE}/g4dev/geant4.${AG4_VERSION}
endif

if ("$?AG4_MAKESHLIB" == 0) then
  unsetenv G4LIB_BUILD_SHARED
else  
  setenv G4LIB_BUILD_SHARED 1
endif  

# path to data files needed by hadronic processes
setenv G4LEVELGAMMADATA ${G4INSTALL}/data/PhotonEvaporation

# This is architecture dependent...
set SYSTEM = `uname`
if ( $SYSTEM == "HP-UX" ) then
  setenv G4SYSTEM "HP-aCC"
  #setenv G4USE_OSPACE 1        # compiling with Object Space STL
endif 
if ( $SYSTEM == "Linux" ) then
  # check compiler version
  set COMPILER = "g++"
  set TMPFILE = `mktemp -u /tmp/g4compiler.XXXXXX`
  ${COMPILER} -v >& ${TMPFILE}
  if ( "`cat ${TMPFILE} | grep 2.96 `" != "" ) then
    echo "WARNING: Found version of 'g++' (2.96.XXX) is known to be buggy!"
  endif   
  if ( "`cat ${TMPFILE} | grep egcs `" != "" ) then
    setenv G4SYSTEM "Linux-egcs"
  else   
    setenv G4SYSTEM "Linux-g++"
  endif  
  if ( "`cat ${TMPFILE} | grep 3.1 `" != "" ) then
    setenv GCC3 1
  endif
  rm ${TMPFILE}
endif
if ( $SYSTEM == "OSF1" ) then
  setenv G4SYSTEM "DEC-cxx"
  #setenv G4NO_STD_NAMESPACE 1  # compiling witn non ISO/ANSI setup
endif 
if ( $SYSTEM == "SunOS" ) then
  setenv G4SYSTEM "SUN-CC"
  setenv G4USE_OSPACE 1         # compiling with Object Space STL
endif 
if ( "$VERBOSE" == "YES" ) then
  echo "Architecture is $SYSTEM"
  echo "Geant4 is istalled in $G4INSTALL"
  echo "Geant4 architecture type is $G4SYSTEM"
  if ( "$?G4USE_OSPACE" == 1 ) then
    echo "ObjectSpace will be used for STL"
  endif
endif

#
# Geant4_mc  
# ==================================
#
if ( "$VERBOSE" == "YES" ) then
  echo " "
  echo "Geant4_mc env. variables..."
  echo "============================"
endif
setenv MCINSTALL ${G4MC_BASE}/geant4_mc
if ( "$VERBOSE" == "YES" ) then
  echo "geant4_mc base directory: $MCINSTALL"
endif

#
# CLHEP
# ==================================
#

if ( -d $LHCXX_BASE/CLHEP/1.8.0.0 ) then
  setenv CLHEP_BASE_DIR $LHCXX_BASE/CLHEP/1.8.0.0
else
  echo "WARNING: CLHEP has not been found in the default path."
  if ( "$VERBOSE" == "YES" ) then
    echo "         Please set the variable CLHEP_BASE_DIR to its base path"
    echo "         Example: setenv  CLHEP_BASE_DIR /afs/cern.ch/sw/lhcxx/specific/@sys/CLHEP/pro"
  endif
endif
if ( "$VERBOSE" == "YES" ) then
  echo "CLHEP base directory: $CLHEP_BASE_DIR"
endif


#
# Visualization
# ==================================
#

if ( "$?AG4_VISUALIZE" == 1 ) then
  if ( "$VERBOSE" == "YES" ) then
    echo "G4 Visualization env. variables..."
  endif

  #
  # tcsh UI
  #
  if ( "$VERBOSE" == "YES" ) then
    echo "* tcsh UI..."
  endif
  setenv G4UI_USE_TCSH         1

  #
  # Xm UI
  #
  if ( "$VERBOSE" == "YES" ) then
    echo "* X11 with Motif..."
  endif
  setenv G4UI_BUILD_XM_DRIVER  1
  setenv G4UI_BUILD_XM_SESSION 1
  setenv G4UI_USE_XM           1

  #
  # Fukui Renderer
  #
  if ( "$VERBOSE" == "YES" ) then
    echo "* Fukui Renderer (DAWN)..."
  endif
  setenv G4VIS_BUILD_DAWN_DRIVER     1
  setenv G4VIS_USE_DAWN              1
  #setenv G4DAWNFILE_VIEWER   david
  setenv DAWN_HOME ${G4_BASE}/../tools/bin
  if ( "`echo ${PATH} | grep ${DAWN_HOME} `" == "" ) then
    setenv PATH "${PATH}:${DAWN_HOME}"
    rehash
  endif
  setenv G4DAWN_MULTI_WINDOW 1
  if ( $SYSTEM == "Linux" ) then
    setenv G4DAWN_NAMED_PIPE 1
  endif

  if ( "$VERBOSE" == "YES" ) then
    if ("$?G4VIS_USE_DAWN" == 1) then
      echo "  Dawn driver activated"
    endif
    if ("$?G4DAWNFILE_VIEWER" == 1) then
      echo "  Dawn file viewer set to ${G4DAWNFILE_VIEWER}"
    endif
    if ("$?DAWN_HOME" == 1) then
      echo "  Dawn home path set to ${DAWN_HOME}"
    endif
    if ("$?G4DAWN_MULTI_WINDOW" == 1) then
      echo "  Dawn multi window selected"
    endif
    if ("$?G4DAWN_NAMED_PIPE" == 1) then
      echo "  Dawn named pipe selected"
    endif
  endif

  # David flags
  # Set colors for overlappings
  setenv DAVID_RAINBOW_FLAG 1
  #setenv DAVID_HIGHLIGHT_COLOR_R  r
  #setenv DAVID_HIGHLIGHT_COLOR_G  g
  #setenv DAVID_HIGHLIGHT_COLOR_B  b

  # If set volumes names are shown
  setenv DAVID_DISPLAY_PVNAME   1
  # If set supresses the call to dawn
  #setenv DAVID_NO_VIEW  1
  setenv DAVID_EPSILON_3D  0.001

  if ( "$VERBOSE" == "YES" ) then
    if ("$?DAVID_RAINBOW_FLAG" == 1) then
      echo "  DAVID Rainbow mode is ON"
    endif
    if ("$?DAVID_HIGHLIGHT_COLOR_R" == 1) then
      echo "  DAVID Highlight color (Red) set to ${DAVID_HIGHLIGHT_COLOR_R}"
    endif
    if ("$?DAVID_HIGHLIGHT_COLOR_G" == 1) then
      echo "  DAVID Highlight color (Green) set to ${DAVID_HIGHLIGHT_COLOR_G}"
    endif
    if ("$?DAVID_HIGHLIGHT_COLOR_B" == 1) then
      echo "  DAVID Highlight color (Blue) set to ${DAVID_HIGHLIGHT_COLOR_B}"
    endif
    if ("$?DAVID_DISPLAY_PVNAME" == 1) then
      echo "  DAVID will display intersected volumes name"
    endif
    if ("$?DAVID_DISPLAY_PVNAME" == 1) then
      echo "  Dawn will not be called from DAVID"
    endif
    if ("$?DAVID_EPSILON_3D" == 1) then
      echo "  DAVID tolerance set to ${DAVID_EPSILON_3D}"
    endif
  endif

  #
  # OpenGL
  #
  if ( "$VERBOSE" == "YES" ) then
    echo "* OpenGL..."
  endif
  setenv G4VIS_BUILD_OPENGLX_DRIVER  1
  setenv G4VIS_BUILD_OPENGLXM_DRIVER 1
  setenv G4VIS_USE_OPENGLX           1
  setenv G4VIS_USE_OPENGLXM          1
  #setenv OGLHOME /usr/local
  #setenv OGLLIBS "-L$OGLHOME/lib -lMesaGLU -lMesaGL"
  setenv OGLHOME /usr
  setenv OGLLIBS "-L$OGLHOME/lib -lGLU -lGL"
  if ( $SYSTEM == "HP-UX" ) then
    setenv OGLLIBS "-L/usr/lib ${OGLLIBS}"
  endif
  if ( $SYSTEM == "OSF1" ) then
    # temporarily excluded
    # due to problems with Root
    unsetenv G4VIS_BUILD_OPENGLX_DRIVER
    unsetenv G4VIS_BUILD_OPENGLXM_DRIVER
    unsetenv G4VIS_USE_OPENGLX
    unsetenv G4VIS_USE_OPENGLXM
    unsetenv OGLHOME
    unsetenv OGLLIBS
  endif
  if ( "$VERBOSE" == "YES" ) then
    if ("$?G4VIS_USE_OPENGLX" == 1) then
      echo "  OpenGL and  X11 driver activated"
    endif
    if ("$?G4VIS_USE_OPENGLXM" == 1) then
      echo "  OpenGL with Motif extension driver activated"
    endif
    if ("$?OGLHOME" == 1) then
      echo "  OpenGL path set to ${OGLHOME}"
    endif
    if ("$?OGLLIBS" == 1) then
      echo "  OpenGL libraries set to ${OGLLIBS}"
    endif
  endif

  #
  # OpenInventor
  #
  if ( "$VERBOSE" == "YES" ) then
    echo "* OpenInventor..."
  endif
  #setenv G4VIS_USE_OPENINVENTOR 1
  #setenv OIHOME whatever
  #setenv HEPVISDIR something
  if ( "$VERBOSE" == "YES" ) then
    if ("$?G4VIS_USE_OPENINVENTOR" == 1) then
      echo "  OpenInventor driver activated"
      echo "  OpenInventor path is ${OIHOME}"
      echo "  HepVis path is ${HEPVISDIR}"
    endif
  endif

  #
  # VRML1
  #
  if ( "$VERBOSE" == "YES" ) then
    echo "* VRML..."
  endif
  setenv G4VIS_BUILD_VRML_DRIVER        1
  setenv G4VIS_USE_VRML                 1
  #Set preferred vrml viewer to be invoked in this mode
  setenv G4VRMLFILE_VIEWER vrweb
  #Set host name for VRML1 visualization.... the g4vrmlview machine!
  setenv G4VRML_HOST_NAME nohost
  #Command to run java listener for VRML driver
  #alias javavrml "java -classpath $G4_BASE/tools/bin/java g4vrmlview vrweb"

  if ( "$VERBOSE" == "YES" ) then
    if ("$?G4VIS_USE_VRML" == 1) then
      echo "  VRML driver activated"
      echo "  Host Name for remote visualization is ${G4VRML_HOST_NAME}"
    endif
  endif

  #
  # GAG
  #
  if ( "$VERBOSE" == "YES" ) then
    echo "* Geant Adaptative GUI (GAG)..."
  endif
  #setenv G4UI_USE_GAG                   1
  #setenv MOMOPATH     ${G4_BASE}/tools/GAG/tcltk
  #if ( "`echo ${PATH} | grep ${MOMOPATH} `" == "" ) then
  #  setenv PATH "${PATH}:${MOMOPATH}"
  #  rehash
  #endif
  #set    NCLASSPATH = ".:${G4_BASE}/tools/swing-1.0.3/swingall.jar:${G4_BASE}/tools/GAG/java/GAG.jar"
  #if ("$?CLASSPATH" == 0) then
  #  setenv CLASSPATH $NCLASSPATH
  #else
  #  if ( "`echo ${CLASSPATH} | grep ${NCLASSPATH} `" == "" ) then
  #    setenv CLASSPATH "${CLASSPATH}:${NCLASSPATH}"
  #  endif
  #endif

  if ( "$VERBOSE" == "YES" ) then
    if ("$?G4UI_USE_GAG" == 1) then
      echo "  GAG UI activated"
      echo "  Momo path set to $MOMOPATH"
      echo "    NOTE: Run "\'tmomo\' "to use Momo (TK/Tcl version)"
      echo "  The path to the java GAG code was updated"
      echo "    NOTE: Run "\'java gag\'" to use GAG (java version)"
    endif
  endif

else
  if ( "$VERBOSE" == "YES" ) then
    echo Unsetting G4 Visualization env. variables...
  endif

  #XM
  unsetenv G4UI_BUILD_XM_DRIVER
  unsetenv G4UI_BUILD_XM_SESSION
  unsetenv G4UI_USE_XM

  # Dawn
  unsetenv G4VIS_BUILD_DAWN_DRIVER
  unsetenv G4VIS_USE_DAWN
  unsetenv G4DAWNFILE_VIEWER
  unsetenv DAWN_HOME
  unsetenv G4DAWN_MULTI_WINDOW
  if ( $SYSTEM == "Linux" ) then
    unsetenv G4DAWN_NAMED_PIPE
  endif

  # David
  unsetenv DAVID_RAINBOW_FLAG
  unsetenv DAVID_HIGHLIGHT_COLOR_R
  unsetenv DAVID_HIGHLIGHT_COLOR_G
  unsetenv DAVID_HIGHLIGHT_COLOR_B
  unsetenv DAVID_DISPLAY_PVNAME
  unsetenv DAVID_NO_VIEW
  unsetenv DAVID_EPSILON_3D
  
  # OpenGL
  unsetenv G4VIS_BUILD_OPENGLX_DRIVER
  unsetenv G4VIS_BUILD_OPENGLXM_DRIVER
  unsetenv G4VIS_USE_OPENGLX
  unsetenv G4VIS_USE_OPENGLXM
  unsetenv OGLHOME
  unsetenv OGLLIBS

  # OpenInventor
  unsetenv G4VIS_USE_OPENINVENTOR

  # VRML1
  unsetenv G4VIS_BUILD_VRML_DRIVER
  unsetenv G4VIS_USE_VRML
  unsetenv G4VRMLFILE_VIEWER
  unsetenv G4VRML_HOST_NAME

  # GAG
  unsetenv G4UI_USE_GAG
  unsetenv MOMOPATH

endif

#
# OPACS
#
if ( "$?AG4_OPACS" == 1 ) then
  if ( "$VERBOSE" == "YES" ) then
    echo "* OPACS..."
  endif
  setenv AG4_VISUALIZE 1

  #
  # OpenGL: needed by OPACS
  #
  setenv G4VIS_BUILD_OPENGLX_DRIVER 1
  setenv G4VIS_USE_OPENGLX          1
  setenv OGLHOME /usr/local
  setenv OGLLIBS "-L$OGLHOME/lib -lMesaGLU -lMesaGL"
  if ( $SYSTEM == "HP-UX" ) then
    setenv OGLLIBS "-L/usr/lib ${OGLLIBS}"
  endif

  #
  # OPACS
  #
  setenv G4VIS_BUILD_OPACS_DRIVER 1
  setenv G4UI_BUILD_WO_DRIVER     1
  setenv G4UI_BUILD_WO_SESSION    1
  setenv G4VIS_USE_OPACS          1
  setenv G4UI_USE_WO              1
  setenv OPACS_HOME $G4_BASE/tools/OPACS
  #setenv OPACS_HOME /afs/cern.ch/rd44/dev/OPACS
  if ( $SYSTEM == "Linux" ) then
    setenv G4_OPACS_WIDGET_SET lesstif
  else
    setenv G4_OPACS_WIDGET_SET Xm
  endif
  source $OPACS_HOME/OPACS/v3/setup.csh
  setenv WOENVIRONMENT $AG4_INSTALL/bin/Alice.odb
  setenv OPATH "$OPATH $G4INSTALL/environments/OPACS/usr"
  if ( "$VERBOSE" == "YES" ) then
    if ("$?G4VIS_USE_OPACS" == 1) then
      echo "  OPACS driver activated"
      echo "  OPACS path set to $OPACS_HOME"
    endif
  endif
else   
  if ( "$VERBOSE" == "YES" ) then
    echo "* Unsetting OPACS driver env. variables ..."
  endif
  unsetenv G4VIS_BUILD_OPACS_DRIVER
  unsetenv G4UI_BUILD_WO_DRIVER
  unsetenv G4UI_BUILD_WO_SESSION
  unsetenv G4VIS_USE_OPACS
  unsetenv G4UI_USE_WO
  unsetenv OPACS_HOME
  unsetenv G4_OPACS_WIDGET_SET
  unsetenv G4OROOT
  unsetenv WOENVIRONMENT
  unsetenv OPATH
endif

#
# path to AliGeant4 config scripts
#  
if ( "`echo ${PATH} | grep ${AG4_INSTALL}/config `" == "" ) then
  if ( "$VERBOSE" == "YES" ) then
    echo Adding ${AG4_INSTALL}/config to the path...
  endif
  setenv PATH "${PATH}:${AG4_INSTALL}/config"
endif

#
# path to shared libraries
# 
if ( $SYSTEM == "HP-UX" ) then
  set SHLIBVAR = $SHLIB_PATH
  set SHLIBVARNAME = SHLIB_PATH
endif 
if ( $SYSTEM == "Linux" ) then
  set SHLIBVAR = $LD_LIBRARY_PATH
  set SHLIBVARNAME = LD_LIBRARY_PATH
endif
if ( $SYSTEM == "OSF1" ) then
  set SHLIBVAR = $LD_LIBRARY_PATH
  set SHLIBVARNAME = LD_LIBRARY_PATH
endif
if ( $SYSTEM == "SunOS" ) then
  set SHLIBVAR = $LD_LIBRARY_PATH
  set SHLIBVARNAME = LD_LIBRARY_PATH
endif

if ( "`echo ${SHLIBVAR} | grep ${G4INSTALL}/lib/${G4SYSTEM} `" == "" ) then
  if ( "$VERBOSE" == "YES" ) then
    echo Adding ${G4INSTALL}/lib/${G4SYSTEM} to the shared libraries path...
  endif
  set SHLIBVAR="${G4INSTALL}/lib/${G4SYSTEM}:${SHLIBVAR}"
endif
if ( "`echo ${SHLIBVAR} | grep ${MCINSTALL}/lib/${G4SYSTEM} `" == "" ) then
  if ( "$VERBOSE" == "YES" ) then
    echo Adding ${MCINSTALL}/lib/${G4SYSTEM} to the shared libraries path...
  endif
  set SHLIBVAR="${MCINSTALL}/lib/${G4SYSTEM}:${SHLIBVAR}"
endif
if ( "`echo ${SHLIBVAR} | grep ${CLHEP_BASE_DIR}/lib `" == "" ) then
  if ( "$VERBOSE" == "YES" ) then
    echo Adding ${CLHEP_BASE_DIR}/lib to the shared libraries path...
  endif
  set SHLIBVAR="${CLHEP_BASE_DIR}/lib:${SHLIBVAR}"
endif

setenv $SHLIBVARNAME $SHLIBVAR


#
# Remove unneeded variables. If this is not done the vars. remain in the env.
#

unset ALICE_BASE
unset G4_BASE
unset LHCXX_BASE
unset SYSTEM
unset NCLASSPATH
unset SHLIBVAR
unset SHLIBVARNAME
unset LOCAL
unset VERBOSE

if ( "$SILENT" == "NO" ) then
  echo "Default ALICE environment for GEANT4 has been set."
endif
