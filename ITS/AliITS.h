#ifndef ITS_H
#define ITS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//           Manager and hits classes for set: ITS                    //
////////////////////////////////////////////////////////////////////////

#include "TObjArray.h"
#include "AliDetector.h"
#include "AliITSgeom.h"
#include "AliITSdigit.h"
#include "AliITSmodule.h"

class AliITS : public AliDetector {
////////////////////////////////////////////////////////////////////////
//
//      An overview of the basic philosophy of the ITS code development
// and analysis is show in the figure below.
//Begin_Html
/*
<img src="picts/ITS/ITS_Analysis_schema.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>Roberto Barbera is in charge of the ITS Offline code (1999).
<a href="mailto:roberto.barbera@ct.infn.it">Roberto Barbera</a>.
</font>
<pre>
*/
//End_Html
//
// Version: 0
// Written by Rene Brun, Federico Carminati, and Roberto Barbera
//
// Version: 1
// Modified and documented by Bjorn S. Nilsen
// July 11 1999
//
// AliITS is the general base class for the ITS. Also see AliDetector for
// futher information.
//
// Data members:
//
// AliITSgeom *fITSgeom
//     All of the geometry information about the active volumes that
// make up the ITS are described in the AliITSgeom class. This includes
// the transformation functions between the local and global coordinate
// systems and the like. See the full description found in the AliITSgeom
// class. Here in the AliITS class is kept the pointer to the geometry
// used in the simulations or that thought to be the correct one for the
// data. Until a more general class is define and a more appropriate
// place is found to keep this geometry information, it is kept here in
// AliITS.
//
// TObjArray *fITSpoints
//     This is a pointer to the points, to be used by the tracking algorithms
// for example, found in the detectors of the ITS. To allow for the most
// general points structure it is defined to be a pointer to a TObjArray where
// each array element would be one point found in the ITS detectors. An
// Addpoints function is defined below. By default an array of 16 TObjects are
// defined during the initialization of AliITS. This is automatically expanded
// when necessary by the Addpoints function.
//
// Bool_t fEuclidOut
//     This is a flag used to indicate that an Euclid compatible CAD
// file will be created upon the creation of the ITS Monte Carlo
// geometry definition, in the function CreatGeometry. If fEuclidOut is
// true, then a file called ITSgeometry.euc will be created.
//
// Int_t fIdN
//     This variable contains the number of layers defined in the ITS
// geometry. It is primarily used as a size indicator for fIdSens and
// fIdName described below. In general the number of layers, ladders, or
// detectors should be gotten from the AliITSgeom functions. Upon
// creating the AliITS object it is set to zero.
//
// Int_t *fIdSens
//     This is a pointer to an array containing the Monte Carlo volume
// numbers for the different layers of the ITS. These numbers are needed
// by the StepManager function to determine what layer a hit was on. It
// is sized and initialized in the Init function and the AliITSv? Init
// function, called after a call to CreateGeometry. Upon creating the
// AliITS object it points to zero. This variable is made a pointer
// in order to keep the maximum flexibility at this level of the code.
//
// char **fIdName
//     This is a pointer to an array of characters containing the names of
// the different ITS layers as defined in the Monte Carlo geometry data-
// base. It is sized and filled in the AliITSv? Init function, called
// after a call to CreatGeometry. Upon creating the AliITS object it
// points to zero. This variable is make a pointer in order to keep the
// maximum flexibility at this level of the code.
//
// Member Functions:
//
// AliITS()
//     The default constructor of the AliITS class. In addition to
// creating the AliITS class it zeros the variables fIshunt (a member
// of AliDetector class), fEuclidOut, and fIdN, and zeros the pointers
// fITSpoints, fIdSens, and fIdName. The AliDetector default constructor
// is also called.
//
// AliITS(const char *name, const char *title)
//     The constructor of the AliITS class. In addition to creating the
// AliITS class, it allocates memory for the TClonesArrays fHits and
// fDigits, and for the TObjArray fITSpoints. It also zeros the variables
// fIshunt (a member of AliDetector class), fEuclidOut, and fIdN, and zeros
// the pointers fIdSens and fIdName. To help in displaying hits via the ROOT
// macro display.C AliITS also sets the marker color to red. The variables
// passes with this constructor, const char *name and *title, are used by
// the constructor of AliDetector class. See AliDetector class for a
// description of these parameters and its constructor functions.
//
// ~AliITS()
//     The default destructor of the AliITS class. In addition to deleting
// the AliITS class it deletes the memory pointed to by the fHits, fDigits,
// fIdSens, fIdName, and fITSpoints.
//
// AddHit(Int_t track, Int_t *vol, Float_t *hits)
//     The function to add information to the AliITShit class. See the
// AliITShit class for a full description. This function allocates the
// necessary new space for the hit information and passes the variable
// track, and the pointers *vol and *hits to the AliITShit constructor
// function.
//
// AddDigit(Int_t *track, Int_t *digits)
//     The function to add information to the AliITSdigits class. See the
// AliITSdigits class for a full description. This function allocates the
// necessary new space for the digits information and passes the pointers
// *track and *digits to the AliITSdigits constructor function.
//
// BuildGeometry()
//     This function builds a simple ITS geometry used by the ROOT macro
// display.C. In general the geometry as coded is wrong.
//
// CreateGeometry()
//     This function builds the detailed geometry used by the Geant
// Monte Carlo. As defined here it is a dummy routine to be replaced
// by the version coded up in AliITSv? where the specific geometry to
// be used by the simulation is defined. See the definition of AliITSv5
// or the other routines for a complete definition.
//
// CreateMaterials()
//     This function defines the default materials used in the Geant
// Monte Carlo simulations. In general it is automatically replaced by
// the CreatMaterials routine defined in AliITSv?. Should the function
// CreateMaterials not exist for the geometry version you are using this
// one is used. See the definition found in AliITSv5 or the other routine
// for a complete definition.
//
// IsVersion()
//     Returns the version number of the AliITS class. At present it is
// version 1.
//
// DistancetoPrimitive(Int_t x, Int_t y)
//     A dummy routine used by the ROOT macro display.C to allow for the
// use of the mouse (pointing device) in the macro. In general this should
// never be called. If it is it returns the number 9999 for any value of
// x and y.
//
// Init()
//     This routine initializes the AliITS class. It is intended to be called
// from the Init function in AliITSv?. Besides displaying a banner
// indicating that it has been called it initializes the array fIdSens.
// Therefore it should be called after a call to CreateGeometry.
//
// MakeBranch(Option_t *Opt=" ")
//     Creates the TTree branch where the class AliITS is kept.
//
// SetEUCLID(bool_t euclid=1)
//     Sets the flag fEuclidOut to true (default) of false (euclid=0).
// By setting or clearing the fEuclidOut flag you can controls whether
// or not a euclid formatted output file of the ITS geometry is written.
// If fEuclidOut is set true then a file called ITSgeometry.euc will be
// written after the ITS geometry is defined in the Monte Carlo. If
// fEuclidOut is set false then no file is created.
//
// StepManager()
//     Dummy routine which is replaced by the routine StepManager() defined
// in AliITSv?. If no such routine exist then this routine returns zero.
// See AliITSv? for a detailed description of the step manager routines.
//
// GetITSgeom()
//     Returns the value of the pointer fITSgeom. This is used to get
// access to the ITS geometry stored in the file. See AliITSgeom for a
// full description of the geometry package.
//
// GetITSpoints()
//     Returns the value of the pointer fITSpoints. This is used to get
// access to the ITS cluster objects, if filled, stored in the file. See
// AliITSCluster for a full description of the cluster data.
////////////////////////////////////////////////////////////////////////
 protected:
    AliITSgeom  *fITSgeom;    // Pointer to ITS geometry
    TObjArray   *fITSmodules; // Pointer to ITS modules
    // Defined here since it doesn't have a place in AliDetector like fDigit
    TObjArray   *fITSpoints;  // Pointer to ITS points

    Bool_t fEuclidOut; // Flag to write out geometry in euclid format
    Int_t  fIdN; // the number of layers
    Int_t  *fIdSens; char **fIdName;   //layer identifier
    // Geometry and Stepmanager version numbers used.
    Int_t fMajorVersion,fMinorVersion;

  public:
                          AliITS();
                          AliITS(const char *name, const char *title);
           virtual        ~AliITS();

           virtual void   AddHit(Int_t, Int_t*, Float_t*);
           virtual void   AddDigit(Int_t*, Int_t*);
           virtual Int_t  AddDigit(AliITSdigit *d);
//         virtual void   AddPoint(); // yet to be defined

           virtual void   BuildGeometry();
           virtual void   CreateGeometry() {};
           virtual void   CreateMaterials();

    inline virtual TObjArray* GetModules() {return fITSmodules;}
    inline virtual TObjArray* GetPoints(){return fITSpoints;}

    inline void GetGeometryVersion(Int_t &a,Int_t &b)
      {a = fMajorVersion;b=fMinorVersion;return;}
    inline virtual Int_t  IsVersion() {return 1;}
                   Int_t  DistancetoPrimitive(Int_t px, Int_t py);
           virtual void   Init();
           virtual void   MakeBranch(Option_t *opt=" ");
    inline virtual void   SetEUCLID(Bool_t euclid=1){fEuclidOut = euclid;}
           virtual void   StepManager()=0;
    //
    // ITS geometry functions
    inline virtual AliITSgeom *GetITSgeom(){return fITSgeom;}
    inline virtual TObjArray  *GetITSpoints(){return fITSpoints;}

    ClassDef(AliITS,1)
};
#endif
