#ifndef ALIEMCALGEOMETRY_H
#define ALIEMCALGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Geometry class  for EMCAL : singleton
// EMCAL consists of a shell of Pb
//                  
//*-- Author: Yves Schutz (SUBATECH)

#include <assert.h> 

// --- ROOT system ---

class TObjArray ;  
class TVector3; 
class TMatrix ; 

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
  const Float_t GetArm2PhiMin() const { return fArm2PhiMin ; }
  const Float_t GetArm2PhiMax() const { return fArm2PhiMax ; }
  const Float_t GetIPDistance()   const { return  fIPDistance  ; } 
  const Float_t GetEnvelop(Int_t index) const { return fEnvelop[index] ; }  
  const Float_t GetShellThickness() const { return fShellThickness ; }
  const Float_t GetZLength() const { return fZLength ; } 

  // material 
  const Float_t GetAmat()   const { return  fAmat ; }  
  const Float_t GetZmat()   const { return  fZmat ; }   
  const Float_t GetDmat()   const { return  fDmat ; }  
  const Float_t GetRmat()   const { return  fRmat ; }  
  const Float_t GetEmat()   const { return  fEmat ; }  
  const Float_t GetLmat()   const { return  fEmat * fRmat ; }  

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
  Float_t fArm1PhiMin ;              // Minimum phi angle covered by Arm 1 
  Float_t fArm1PhiMax ;              // Maximum phi angle covered by Arm 1       
  Float_t fArm2PhiMin ;              // Minimum phi angle covered by Arm 2        
  Float_t fArm2PhiMax ;              // Maximum phi angle covered by Arm 2
  // It is assumed that Arm1 and Arm2 have the same following parameters
  Float_t fEnvelop[3] ;              // the GEANT TUB that contains the 2 arms
  Float_t fIPDistance ;              // Distance of the inner surface to the interaction point
  Float_t fShellThickness ;          // Total thickness in (x,y) direction
  Float_t fZLength ;                 // Total length in z direction

  //material
  Float_t fAmat ;  // average atomic weight of the active material
  Float_t fZmat ;  // average atomic number of the active material 
  Float_t fDmat ;  // average density of the active material
  Float_t fRmat ;  // average radiation length of the active material
  Float_t fEmat ;  // thickness of the active material in radiation length units
 
  ClassDef(AliEMCALGeometry,1)       // EMCAL geometry class 

} ;

#endif // AliEMCALGEOMETRY_H
