# $Id$
# ----------------------------------------------------------------
# This script sets the default environment variables for
# Alice Geant4 based prototype
# Options: -g     verbose mode
#          local  recomended local configuration
#          silent no output
#   
# by I. Hrivnacova, 18.8.1998
#
# sh version modified by I. Gonzalez 2.3.2000

#
# ==================================================================
# Alice configuration options: 
# Please define your selection with the variables below.
# ==================================================================
# 

#
# ====== AG4_VERSION
# Geant4 version
# If set: the provided Geant4 version and not the default one is set
#export AG4_VERSION=2.0_opt_global

#
# ====== AG4_VISUALIZE
# Set/Unset to get/avoid Geant4 visualisation.
#export AG4_VISUALIZE=1
unset  AG4_VISUALIZE

#
# ====== AG4_OPACS
# Set/Unset to get OPACS as a Geant4 visualisation driver.
#export AG4_OPACS=1
unset AG4_OPACS

#
# ====== AG4_MAKESHLIB
# If set: shared libraris are created
export AG4_MAKESHLIB=1
#unset AG4_MAKESHLIB


#
# Resolve input parameters
#
VERBOSE="NO"
LOCAL="NO"
SILENT="NO"
for param in $*
do
  case $param in
    -g) VERBOSE="YES"; shift 1;;
    local) LOCAL="YES"; shift 1;;
    silent) SILENT="YES"; shift 1;;
  esac
done


#
# ==================================================================
# Key path variables
# ==================================================================
#

if [ "$LOCAL" = "NO" ]; then
  # 
  # AFS
  #

  # ====== ALICE_BASE
  # ALICE base directory
  ALICE_BASE=/afs/cern.ch/alice/offline

  # ====== G4_BASE
  # Geant4 base directory
  G4_BASE=${ALICE_BASE}/geant4

  # ====== G4MC_BASE
  # Geant4_mc base directory
  G4MC_BASE=${ALICE_BASE}

  # ====== LHCXX_BASE
  # LHC++ base directory
  LHCXX_BASE=/afs/cern.ch/sw/lhcxx/specific/@sys

  # ====== ROOT_BASE
  # Root base directory  
  ROOT_BASE=/afs/cern.ch/alice/library/root
  # to be removed after aCC support will be defined
  # in the offline in a standard way  
  if [ `uname` = "HP-UX" ]; then 
   ROOT_BASE=/afs/cern.ch/alice/library/.hp_ux102/root.2.23.aCC
  fi 

  # ====== IRST_BASE
  # IRST code check tool base directory
  export IRST_BASE=/afs/cern.ch/alice/offline/geant4/tools/IRST

else
  # 
  # recommended local installation
  #

  # ====== ALICE_BASE
  # ALICE base directory
  ALICE_BASE=$ALICE_INSTALL/$ALICE_LEVEL

  # ====== G4_BASE
  # Geant4 base directory
  G4_BASE=$ALICE_BASE

  # ====== G4MC_BASE
  # Geant4_mc base directory
  G4MC_BASE=$ALICE_BASE

  # ====== LHCXX_BASE
  # LHC++ base directory
  LHCXX_BASE=$ALICE_BASE

  # ====== ROOT_BASE
  # Root base directory  
  ROOT_BASE=$ALICE_BASE

  # ====== IRST_BASE
  # IRST code check tool base directory
  #export IRST_BASE=$HOME/dev/tools/IRST

fi


#....................................... SKIP ................................
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# You should not need to change 
# the lines below this one
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

# ====== AG4_INSTALL
# Alice Geant4 install directory
#
export AG4_INSTALL=$ALICE_ROOT/AliGeant4

# Set variables depending on other variables
# if opacs is selected then select also visualize
if [ $AG4_OPACS ]; then
  export AG4_VISUALIZE=1
fi

# Start the output
if [ "$VERBOSE" = "YES" ]; then
  echo " "
  echo "    =========================================="
  echo "      ALICE Geant4 environment configuration"
  echo "    =========================================="

  #
  # Some output on the selections... 
  #

  if [ $AG4_VISUALIZE ]; then
    echo "Geant4 visualisation is     selected."
  else
    echo "Geant4 visualisation is NOT selected."
  fi
  if [ $AG4_OPACS ]; then
    echo "OPACS driver         is     selected."
  else
    echo "OPACS driver         is NOT selected."
  fi
  if [ $AG4_NOPHYSICS ]; then
    echo "Only geantino or charged geantino can be shooted."
  else
    echo "Full physics has been selected."
  fi
  if [ $AG4_ODBMS ]; then
    echo "The environment for using Objectivity will be set."
  else
    echo "No Geant4 Persistency."
  fi

fi

#
# ROOT Alice definitions & options
# ==================================
#
if [ "$ROOTSYS" = "" ]; then
  export ROOTSYS=$ROOT_BASE
fi
if [ "$VERBOSE" = "YES" ]; then
  echo " "
  echo "ROOT"
  echo "===="
  echo "ROOTSYS set to $ROOTSYS"
fi

#
# Geant4
# ==================================
#
if [ "$VERBOSE" = "YES" ]; then
  echo " "
  echo "Geant4 env. variables..."
  echo "============================"
fi
if [ "$AG4_VERSION" = "" ]; then
  export G4INSTALL=$G4_BASE/geant4
else
  export G4INSTALL=$G4_BASE/g4dev/geant4.$AG4_VERSION
fi
export MCINSTALL=$ALICE_ROOT/geant4_mc

if [ "$AG4_MAKESHLIB" = "" ]; then
  unset G4LIB_BUILD_SHARED
else
  export G4LIB_BUILD_SHARED=1
fi  

# path to data files needed by hadronic processes
export G4LEVELGAMMADATA=$G4INSTALL/data/PhotonEvaporation

# This is architecture dependent...
SYSTEM=`uname`
if [ "$SYSTEM" = "HP-UX" ]; then
  export G4SYSTEM="HP-aCC"
  #export G4USE_OSPACE=1        # compiling with Object Space STL
fi  
if [ "$SYSTEM" = "Linux" ]; then
  # check compiler version
  COMPILER="g++"
  TMPFILE=`mktemp -u /tmp/g4compiler.XXXXXX`
  ${COMPILER} -v 2> ${TMPFILE}
  if [ "`cat ${TMPFILE} | grep 2.96 `" != "" ]; then
    echo "WARNING: Found version of 'g++' (2.96.XXX) is known to be buggy!"
  fi   
  if [ "`cat ${TMPFILE} | grep egcs `" != "" ]; then
    export G4SYSTEM="Linux-egcs"
  else   
    export G4SYSTEM="Linux-g++"
  fi 
  if [ "`cat ${TMPFILE} | grep 3.1 `" != "" ]; then
    export GCC3=1
  fi  
  rm ${TMPFILE}
fi
if [ "$SYSTEM" = "OSF1" ]; then
  export G4SYSTEM="DEC-cxx"
  #export G4NO_STD_NAMESPACE=1  # compiling with non ISO/ANSI setup
fi
if [ "$SYSTEM" = "SunOS" ]; then
  export G4SYSTEM="SUN-CC"
  export G4USE_OSPACE=1         # compiling with Object Space STL
fi
if [ "$VERBOSE" = "YES" ]; then
  echo "Architecture is $SYSTEM"
  echo "Geant4 is istalled in $G4INSTALL"
  echo "Geant4 architecture type is $G4SYSTEM"
  if [ $G4USE_OSPACE ]; then
    echo "ObjectSpace will be used for STL"
  fi
fi

#
# Geant4_mc  
# ==================================
#
if [ "$VERBOSE" = "YES" ]; then
  echo " "
  echo "Geant4_mc env. variables..."
  echo "============================"
fi
export MCINSTALL=$G4MC_BASE/geant4_mc
if [ "$VERBOSE" = "YES" ]; then
  echo "geant4_mc base directory: $MCINSTALL"
fi

#
# CLHEP
# ==================================
#

if [ -d $LHCXX_BASE/CLHEP/1.8.0.0 ]; then
  export CLHEP_BASE_DIR=$LHCXX_BASE/CLHEP/1.8.0.0
else
  echo "WARNING: CLHEP has not been found in the default path."
  if [ "$VERBOSE" = "YES" ]; then
    echo "         Please set the variable CLHEP_BASE_DIR to its base path"
    echo "         Example: export CLHEP_BASE_DIR=/afs/cern.ch/sw/lhcxx/specific/@sys/CLHEP/pro"
  fi
fi
if [ "$VERBOSE" = "YES" ]; then
  echo "CLHEP base directory: $CLHEP_BASE_DIR"
fi


#
# Visualization
# ==================================
#

if [ $AG4_VISUALIZE ]; then
  if [ "$VERBOSE" = "YES" ]; then
    echo "G4 Visualization env. variables..."
  fi

  #
  # tcsh UI
  #
  if [ "$VERBOSE" = "YES" ]; then
    echo "* tcsh UI..."
  fi
  export G4UI_USE_TCSH=1

  #
  # Xm UI
  #
  if [ "$VERBOSE" = "YES" ]; then
    echo "* X11 with Motif..."
  fi
  export G4UI_BUILD_XM_DRIVER=1
  export G4UI_BUILD_XM_SESSION=1
  export G4UI_USE_XM=1

  #
  # Fukui Renderer
  #
  if [ "$VERBOSE" = "YES" ]; then
    echo "* Fukui Renderer (DAWN)..."
  fi
  export G4VIS_BUILD_DAWN_DRIVER=1
  export G4VIS_USE_DAWN=1
  #export G4DAWNFILE_VIEWER=david
  export DAWN_HOME=${G4_BASE}/tools/bin
  if [ "`echo ${PATH} | grep ${DAWN_HOME} `" = "" ]; then
    export PATH=$PATH:$DAWN_HOME
  fi
  export G4DAWN_MULTI_WINDOW=1
  if [ `uname` = "Linux" ]; then
    export G4DAWN_NAMED_PIPE=1
  fi

  if [ "$VERBOSE" = "YES" ]; then
    if [ $G4VIS_USE_DAWN ]; then
      echo "  Dawn driver activated"
    fi
    if [ $G4DAWNFILE_VIEWER ]; then
      echo "  Dawn file viewer set to ${G4DAWNFILE_VIEWER}"
    fi
    if [ $DAWN_HOME ]; then
      echo "  Dawn home path set to ${DAWN_HOME}"
    fi
    if [ $G4DAWN_MULTI_WINDOW ]; then
      echo "  Dawn multi window selected"
    fi
    if [ $G4DAWN_NAMED_PIPE ]; then
      echo "  Dawn named pipe selected"
    fi
  fi


  # David flags
  # Set colors for overlappings
  export DAVID_RAINBOW_FLAG=1
  #export DAVID_HIGHLIGHT_COLOR_R=r
  #export DAVID_HIGHLIGHT_COLOR_G=g
  #export DAVID_HIGHLIGHT_COLOR_B=b

  # If set volumes names are shown
  export DAVID_DISPLAY_PVNAME=1
  # If set supresses the call to dawn
  #export DAVID_NO_VIEW=1
  export DAVID_EPSILON_3D=0.001

  if [ "$1" = "-g"  ]; then
    if [ $DAVID_RAINBOW_FLAG ]; then
      echo "  DAVID Rainbow mode is ON"
    fi
    if [ $DAVID_HIGHLIGHT_COLOR_R ]; then
      echo "  DAVID Highlight color (Red) set to ${DAVID_HIGHLIGHT_COLOR_R}"
    fi
    if [ $DAVID_HIGHLIGHT_COLOR_G ]; then
      echo "  DAVID Highlight color (Green) set to ${DAVID_HIGHLIGHT_COLOR_G}"
    fi
    if [ $DAVID_HIGHLIGHT_COLOR_B ]; then
      echo "  DAVID Highlight color (Blue) set to ${DAVID_HIGHLIGHT_COLOR_B}"
    fi
    if [ $DAVID_DISPLAY_PVNAME ]; then
      echo "  DAVID will display intersected volumes name"
    fi
    if [ $DAVID_DISPLAY_PVNAME ]; then
      echo "  Dawn will not be called from DAVID"
    fi
    if [ $DAVID_EPSILON_3D ]; then
      echo "  DAVID tolerance set to ${DAVID_EPSILON_3D}"
    fi
  fi

  #
  # OpenGL
  #
  if [ "$VERBOSE" = "YES" ]; then
    echo "* OpenGL..."
  fi
  export G4VIS_BUILD_OPENGLX_DRIVER=1
  export G4VIS_BUILD_OPENGLXM_DRIVER=1
  export G4VIS_USE_OPENGLX=1
  export G4VIS_USE_OPENGLXM=1
  #export OGLHOME=/usr/local
  export OGLHOME=/export/alice/tools
  #export OGLLIBS="-L$OGLHOME/lib -lMesaGLU -lMesaGL"
  export OGLLIBS="-L$OGLHOME/lib -lGLU -lGL"
  if [ "$SYSTEM" = "HP-UX" ]; then
    export OGLLIBS="-L/usr/lib $OGLLIBS"
  fi
  if [ "$SYSTEM" = "OSF1" ]; then
    # temporarily excluded
    # due to problems with Root
    unset G4VIS_BUILD_OPENGLX_DRIVER
    unset G4VIS_BUILD_OPENGLXM_DRIVER
    unset G4VIS_USE_OPENGLX
    unset G4VIS_USE_OPENGLXM
    unset OGLHOME
    unset OGLLIBS
  fi
  if [ "$VERBOSE" = "YES" ]; then
    if [ $G4VIS_USE_OPENGLX ]; then
      echo "  OpenGL and  X11 driver activated"
    fi
    if [ $G4VIS_USE_OPENGLXM ]; then
      echo "  OpenGL with Motif extension driver activated"
    fi
    if [ $OGLHOME ]; then
      echo "  OpenGL path set to ${OGLHOME}"
    fi
    if [ $OGLLIBS ]; then
      echo "  OpenGL libraries set to ${OGLLIBS}"
    fi
  fi

  #
  # OpenInventor
  #
  if [ "$VERBOSE" = "YES" ]; then
    echo "* OpenInventor..."
  fi
  #export G4VIS_USE_OPENINVENTOR=1
  #export OIHOME=whatever
  #export HEPVISDIR=something
  if [ "$VERBOSE" = "YES" ]; then
    if [ $G4VIS_USE_OPENINVENTOR ]; then
      echo "  OpenInventor driver activated"
      echo "  OpenInventor path is ${OIHOME}"
      echo "  HepVis path is ${HEPVISDIR}"
    fi
  fi

  #
  # VRML1
  #
  if [ "$VERBOSE" = "YES" ]; then
    echo "* VRML..."
  fi
  export G4VIS_BUILD_VRML_DRIVER=1    
  export G4VIS_USE_VRML=1
  #Set preferred vrml viewer to be invoked in this mode
  export G4VRMLFILE_VIEWER=vrweb
  #Set host name for VRML1 visualization.... the g4vrmlview machine!
  export G4VRML_HOST_NAME=nohost
  #Command to run java listener for VRML driver
  #alias javavrml "java -classpath $G4_BASE/tools/bin/java g4vrmlview vrweb"

  if [ "$VERBOSE" = "YES" ]; then
    if [ $G4VIS_USE_VRML ]; then
      echo "  VRML driver activated"
      echo "  Host Name for remote visualization is ${G4VRML_HOST_NAME}"
    fi
  fi

  #
  # GAG
  #
  if [ "$VERBOSE" = "YES" ]; then
    echo "* Geant Adaptative GUI (GAG)..."
  fi
  export G4UI_USE_GAG=1
  export MOMOPATH=${G4_BASE}/tools/GAG/tcltk
  if [ "`echo ${PATH} | grep ${MOMOPATH} `" = "" ]; then
    export PATH=$PATH:$MOMOPATH
  fi
  NCLASSPATH=".:${G4_BASE}/tools/swing-1.0.3/swingall.jar:${G4_BASE}/tools/GAG/java/GAG.jar"
  if [ "$CLASSPATH" = "" ]; then
    export CLASSPATH=$NCLASSPATH
  else
    if [ "`echo ${CLASSPATH} | grep ${NCLASSPATH} `" = "" ]; then
      export CLASSPATH="${CLASSPATH}:${NCLASSPATH}"
    fi
  fi

  if [ "$VERBOSE" = "YES" ]; then
    if [ $G4UI_USE_GAG ]; then
      echo "  GAG UI activated"
      echo "  Momo path set to $MOMOPATH"
      echo "    NOTE: Run "\'tmomo\' "to use Momo (TK/Tcl version)"
      echo "  The path to the java GAG code was updated"
      echo "    NOTE: Run "\'java gag\'" to use GAG (java version)"
    fi
  fi


else
  if [ "$VERBOSE" = "YES" ]; then
    echo Unsetting G4 Visualization env. variables...
  fi

  #XM
  unset G4UI_BUILD_XM_DRIVER
  unset G4UI_BUILD_XM_SESSION
  unset G4UI_USE_XM

  # Dawn
  unset G4VIS_BUILD_DAWN_DRIVER
  unset G4VIS_USE_DAWN
  unset G4DAWNFILE_VIEWER
  unset DAWN_HOME
  unset G4DAWN_MULTI_WINDOW
  if [ "$SYSTEM" = "Linux" ]; then
    unset G4DAWN_NAMED_PIPE
  fi

  # David
  unset DAVID_RAINBOW_FLAG
  unset DAVID_HIGHLIGHT_COLOR_R
  unset DAVID_HIGHLIGHT_COLOR_G
  unset DAVID_HIGHLIGHT_COLOR_B
  unset DAVID_DISPLAY_PVNAME
  unset DAVID_NO_VIEW
  unset DAVID_EPSILON_3D
  
  # OpenGL
  unset G4VIS_BUILD_OPENGLX_DRIVER
  unset G4VIS_BUILD_OPENGLXM_DRIVER
  unset G4VIS_USE_OPENGLX
  unset G4VIS_USE_OPENGLXM
  unset OGLHOME
  unset OGLLIBS

  # OpenInventor
  #unset G4VIS_USE_OPENINVENTOR

  # VRML1
  unset G4VIS_BUILD_VRML_DRIVER
  unset G4VIS_USE_VRML
  unset G4VRMLFILE_VIEWER
  unset G4VRML_HOST_NAME

  # GAG
  unset G4UI_USE_GAG
  unset MOMOPATH

fi

#
# OPACS
#
if [ $AG4_OPACS ]; then
  if [ "$VERBOSE" = "YES" ]; then
    echo "* OPACS..."
  fi
  export AG4_VISUALIZE=1

  #
  # OpenGL: needed by OPACS
  #
  export G4VIS_BUILD_OPENGLX_DRIVER=1
  export G4VIS_USE_OPENGLX=1
  export OGLHOME=/usr/local
  export OGLLIBS="-L$OGLHOME/lib -lMesaGLU -lMesaGL"
  if [ "$SYSTEM" = "HP-UX" ]; then
    export OGLLIBS="-L/usr/lib $OGLLIBS"
  fi
    
  #
  # OPACS
  #
  export G4VIS_BUILD_OPACS_DRIVER=1
  export G4UI_BUILD_WO_DRIVER=1
  export G4UI_BUILD_WO_SESSION=1
  export G4VIS_USE_OPACS=1
  export G4UI_USE_WO=1
  export OPACS_HOME=$G4_BASE/tools/OPACS
  #export OPACS_HOME=/afs/cern.ch/rd44/dev/OPACS
  if [ "$SYSTEM" = "Linux" ]; then
    export G4_OPACS_WIDGET_SET=lesstif
  else
    export G4_OPACS_WIDGET_SET=Xm
  fi
  . $OPACS_HOME/OPACS/v3/setup.sh
  export WOENVIRONMENT=$AG4_INSTALL/bin/Alice.odb
  export OPATH="$OPATH $G4INSTALL/environments/OPACS/usr"
  if [ "$VERBOSE" = "YES" ]; then
    if [ $G4VIS_USE_OPACS ]; then
      echo "  OPACS driver activated"
      echo "  OPACS path set to $OPACS_HOME"
    fi
  fi
else   
  if [ "$VERBOSE" = "YES" ]; then
    echo "* Unsetting OPACS driver env. variables ..."
  fi
  unset G4VIS_BUILD_OPACS_DRIVER
  unset G4UI_BUILD_WO_DRIVER
  unset G4UI_BUILD_WO_SESSION
  unset G4VIS_USE_OPACS
  unset G4UI_USE_WO
  unset OPACS_HOME
  unset G4_OPACS_WIDGET_SET
  unset G4OROOT
  unset WOENVIRONMENT
  unset OPATH
fi

#
# path to AliGeant4 config scripts
#  
if [ "`echo ${PATH} | grep ${AG4_INSTALL}/config `" = "" ]; then
  if [ "$VERBOSE" = "YES" ]; then
    echo Adding ${AG4_INSTALL}/config to the path...
  fi
  export PATH=${PATH}:${AG4_INSTALL}/config
fi

#
# path to shared libraries
# 
if [ "$SYSTEM" = "HP-UX" ]; then
  SHLIBVAR=$SHLIB_PATH
  SHLIBVARNAME=SHLIB_PATH
fi  
if [ "$SYSTEM" = "Linux" ]; then
  SHLIBVAR=$LD_LIBRARY_PATH
  SHLIBVARNAME=LD_LIBRARY_PATH
fi
if [ "$SYSTEM" = "OSF1" ]; then
  SHLIBVAR=$LD_LIBRARY_PATH
  SHLIBVARNAME=LD_LIBRARY_PATH
fi
if [ "$SYSTEM" = "SunOS" ]; then
  SHLIBVAR=$LD_LIBRARY_PATH
  SHLIBVARNAME=LD_LIBRARY_PATH
fi

if [ "`echo ${SHLIBVAR} | grep ${G4INSTALL}/lib/${G4SYSTEM} `" = "" ]; then
  if [ "$VERBOSE" = "YES" ]; then
    echo Adding ${G4INSTALL}/lib/${G4SYSTEM} to the shared libraries path...
  fi
  SHLIBVAR="${G4INSTALL}/lib/${G4SYSTEM}:${SHLIBVAR}"
fi
if [ "`echo ${SHLIBVAR} | grep ${MCINSTALL}/lib/${G4SYSTEM} `" = "" ]; then
  if [ "$VERBOSE" = "YES" ]; then
    echo Adding ${MCINSTALL}/lib/${G4SYSTEM} to the shared libraries path...
  fi
  SHLIBVAR="${MCINSTALL}/lib/${G4SYSTEM}:${SHLIBVAR}"
fi
if [ "`echo ${SHLIBVAR} | grep ${CLHEP_BASE_DIR}/lib `" = "" ]; then
  if [ "$VERBOSE" = "YES" ]; then
    echo Adding ${CLHEP_BASE_DIR}/lib to the shared libraries path...
  fi
  SHLIBVAR="${CLHEP_BASE_DIR}/lib:${SHLIBVAR}"
fi

export $SHLIBVARNAME=$SHLIBVAR


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

if [ "$SILENT" = "NO" ]; then
  echo "Default ALICE environment for GEANT4 has been set."
fi
