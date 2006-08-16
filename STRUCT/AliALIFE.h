#ifndef ALIALIFE_H
#define ALIALIFE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TObject.h"
#include "TString.h"

class AliALIFE : public  TObject {
 public:
    AliALIFE(const char *name1, const char *name2);
    AliALIFE();    
    AliALIFE(const AliALIFE &rhs);

    virtual ~AliALIFE() {}
    void Cylinder(Float_t rmin, Float_t rmax,
		  Float_t zmin, Float_t zmax,
		  Float_t pos[3],
		  const char* Material, 
		  const char* Field="MF", 
		  const char* Cuts="$UNSHIELDED");
    void OnionCylinder(Float_t* r, Int_t nr, Float_t zmin, Float_t zmax,
		       Float_t pos[3],
		       const char** Materials,
		       const char** Fields=0,
		       const char** Cuts=0);
    
    void Cone(Float_t rmin1, Float_t rmin2,
	      Float_t rmax1, Float_t rmax2,
	      Float_t zmin, Float_t zmax,
	      Float_t pos[3],
	      const char* Material, 
	      const char* Field="MF", 
	      const char* Cuts="$UNSHIELDED");
    
    void OnionCone(Float_t* r1, Float_t* r2, Int_t nr,
		   Float_t zmin, Float_t zmax,
		   Float_t pos[3],
		   const char** Materials,
		   const char** Fields=0,
		   const char** Cuts=0);

    void PolyCone(Float_t* rmin, Float_t* rmax, Float_t* z, Int_t nz,
		  Float_t pos[3], 
		  const char* Material, 
		  const char* Field="MF", 
		  const char* Cuts="$UNSHIELDED");

    void OnionPolyCone(Float_t** r , Float_t* z, Int_t nr, Int_t nz,
		       Float_t pos[3], 
		       const char** Materials,
		       const char** Fields=0,
		       const char** Cuts=0);
    
    void Comment(const char* Comment);

    void Finish(Bool_t remove = 0);

    void SetDefaultVolume(TString vol1, TString vol2) 
	{fDefaultVolume1=vol1; fDefaultVolume2=vol2;}
    
    void SetDefaultVolume(TString vol) 
	{fDefaultVolume1=vol;}

 protected:
    Int_t        fNBodies;          // current number of bodies
    Int_t        fNVolumes;         // current number of volumes
    TString       fBodyFile;        // File for Fluka bodies
    TString       fVolumeFile;      // File for Fluka volumes
    FILE         *fFile1;           // ! output file for fluka geometry in ALIFE format  
    FILE         *fFile2;           // ! scratch file
    TString      fDefaultVolume1;   // default external volume 1 
    TString      fDefaultVolume2;   // default external volume 2   
 private:
    void BodyHeader();
    void VolumeHeader();
    void Copy(TObject &alife) const;

    AliALIFE & operator = (const AliALIFE &rhs) {
	rhs.Copy(*this);
	return (*this);
    }
   ClassDef(AliALIFE,1)
};
#endif















