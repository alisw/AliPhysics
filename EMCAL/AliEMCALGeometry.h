#ifndef ALIEMCALGEOMETRY_H
#define ALIEMCALGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Geometry class  for EMCAL : singleton
// EMCAL consists of a layers of scintillator, and lead.
//                  
//*-- Author: Sahal Yacoob (LBL / UCT)
//*--   and : Yves Schutz (Subatech)
        
#include <assert.h> 

// --- ROOT system ---
#include "TString.h"
#include "TObjArray.h"
#include "TVector3.h"

//class TObjArray ;  
//class TVector3; 
//class TMatrix ; 

// --- AliRoot header files ---

#include "AliGeometry.h"


class AliEMCALGeometry : public AliGeometry {

public: 

  AliEMCALGeometry() {
    // default ctor 
    // must be kept public for root persistency purposes, but should never be called by the outside world
  } ;  

  AliEMCALGeometry(const AliEMCALGeometry & geom) {
    // cpy ctor requested by Coding Convention but not yet needed
    assert(0==1) ;
  } 
  
  virtual ~AliEMCALGeometry(void) ; 
  static AliEMCALGeometry * GetInstance(const Text_t* name, const Text_t* title="") ; 
  static AliEMCALGeometry * GetInstance() ; 

  AliEMCALGeometry & operator = (const AliEMCALGeometry  & rvalue) const {
    // assignement operator requested by coding convention but not needed
    assert(0==1) ;
    return *(GetInstance()) ; 
  }
  virtual void GetGlobal(const AliRecPoint *, TVector3 &, TMatrix &) const {}
  virtual void GetGlobal(const AliRecPoint *, TVector3 &) const {}
  // General

  Bool_t  IsInitialized(void) const { return fgInit ; }  
                                                                       
  // Return EMCA geometrical parameters

  // geometry
  const Float_t GetAirGap() const { return fAirGap ; }
  const Float_t GetArm1PhiMin() const { return fArm1PhiMin ; }
  const Float_t GetArm1PhiMax() const { return fArm1PhiMax ; }
  const Float_t GetIPDistance()   const { return  fIPDistance  ; } 
  const Float_t GetEnvelop(Int_t index) const { return fEnvelop[index] ; }  
  const Float_t GetShellThickness() const { return fShellThickness ; }
  const Float_t GetZLength() const { return fZLength ; } 
  const Float_t GetGap2Active() const {return  fGap2Active ; }
  const Int_t   GetNLayers() const {return fNLayers ;}
  const Int_t   GetNZ() const {return fNZ ;}
  const Int_t   GetNPhi() const {return fNPhi ;}
  
protected:

  AliEMCALGeometry(const Text_t* name, const Text_t* title="") : AliGeometry(name, title) { 
    // ctor only for internal usage (singleton)
    Init() ; 
  }
  void Init(void) ;            // initializes the parameters of EMCAL 

private:

  static AliEMCALGeometry * fgGeom ; // pointer to the unique instance of the singleton 
  static Bool_t fgInit ;             // Tells if geometry has been succesfully set up 

  // geometry
  Float_t fAirGap ;                  // Distance between envelop and active material 
  Float_t fArm1PhiMin ;              // Minimum angular position of EMCAL in Phi (degrees)
  Float_t fArm1PhiMax ;              // Maximum angular position of EMCAL in Phi (degrees)

// It is assumed that Arm1 and Arm2 have the same following parameters
  Float_t fEnvelop[3] ;              // the GEANT TUB for the detector 
  Float_t fIPDistance ;              // Distance of the inner surface to the interaction point
  Float_t fShellThickness ;          // Total thickness in (x,y) direction
  Float_t fZLength ;                 // Total length in z direction
  Float_t fGap2Active ;              // Gap between the envelop and the active material
  Int_t fNLayers ;                  // Number of layers of material in the R direction
  Int_t fNZ ;                      // Number of Towers in the Z direction
  Int_t fNPhi ;                    //Number of Towers in the Phi Direction
 
  ClassDef(AliEMCALGeometry,2)       // EMCAL geometry class 

} ;

#endif // AliEMCALGEOMETRY_H
