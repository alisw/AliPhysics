#ifndef ALIPHOSVFAST_H
#define ALIPHOSVFAST_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  Manager class  for PHOS                   //
//  Version SUBATECH                          //
//  Author  Y. Schutz SUBATECH                //
//       This is the class to be used for     //  
//       fast simulations                     //
////////////////////////////////////////////////

/* $Id$ */

// --- ROOT system ---
#include "TClonesArray.h"

// --- AliRoot header files ---
#include "AliPHOS.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSPID.h"


class AliPHOSvFast : public AliPHOS {

public:

  AliPHOSvFast(void) ;
  AliPHOSvFast(const char *name, const char *title="") ;
  virtual ~AliPHOSvFast(void) ;

  void           AddRecParticle(Int_t primary) ;                    // adds primary particle to the RecParticles list
  virtual void   BuildGeometry(void) ;                              // creates the geometry for the ROOT display
  virtual void   CreateGeometry(void) ;                             // creates the geometry for GEANT
  Float_t        GetBigBox(Int_t index) ;                             
  virtual AliPHOSGeometry * GetGeometry() { return fGeom ; }  
  virtual void   Init(void) ;                                       // does nothing
  Int_t   IsVersion(void) const { return -1 ; }
  void           MakeBranch(Option_t* opt) ;
  RecParticlesList * RecParticles() { return fRecParticles ; }      // gets TClonesArray of reconstructed particles
  void           SetBigBox(Int_t index, Float_t value) ;                             
  virtual void   StepManager(void) ;                                // does the tracking through PHOS and a preliminary digitalization
  
private:
  
  Float_t fBigBoxX ;                    // main box containing all PHOS (EMC+PPSD)
  Float_t fBigBoxY ;                    // main box containing all PHOS (EMC+PPSD)
  Float_t fBigBoxZ ;                    // main box containing all PHOS (EMC+PPSD)
  AliPHOSGeometry * fGeom ;             // geometry definition
  Int_t fNRecParticles ;                // number of detected particles
  RecParticlesList * fRecParticles ;    // list of detected particles 

  ClassDef(AliPHOSvFast,1)  // PHOS main class , version for fast simulation

};

#endif // AliPHOSVFAST_H
