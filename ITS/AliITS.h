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

           virtual void   BuildGeometry() {};
           virtual void   CreateGeometry() {};
           virtual void   CreateMaterials() {};

           virtual TObjArray* GetModules() const {return fITSmodules;}
           virtual TObjArray* GetPoints() const {return fITSpoints;}

           void GetGeometryVersion(Int_t &a,Int_t &b) const 
      {a = fMajorVersion;b=fMinorVersion;return;}
           virtual Int_t  IsVersion() const {return 1;}
                   Int_t  DistancetoPrimitive(Int_t px, Int_t py);
           virtual void   Init();
           virtual void   MakeBranch(Option_t *opt=" ");
           virtual void   SetEUCLID(Bool_t euclid=1) {fEuclidOut = euclid;}
           virtual void   StepManager()=0;
    //
    // ITS geometry functions
           virtual AliITSgeom *GetITSgeom() const {return fITSgeom;}
           virtual TObjArray  *GetITSpoints() const {return fITSpoints;}

    ClassDef(AliITS,1)
};
#endif
