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
#include "TRandom.h"

// --- AliRoot header files ---
#include "AliPHOS.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSPID.h"
#include "AliPHOSFastRecParticle.h"

class AliPHOSvFast : public AliPHOS {

public:

  AliPHOSvFast(void) ;
  AliPHOSvFast(const char *name, const char *title="") ;
  virtual ~AliPHOSvFast(void) ;

  void           AddRecParticle(const AliPHOSFastRecParticle & rp) ; // adds primary particle to the RecParticles list
  virtual void   BuildGeometry(void) ;                               // creates the geometry for the ROOT display
  virtual void   CreateGeometry(void) ;                              // creates the geometry for GEANT
  Float_t        GetBigBox(Int_t index) ;                             
  virtual AliPHOSGeometry * GetGeometry() { return fGeom ; }  
  virtual void   Init(void) ;                                        // does nothing
  Int_t   IsVersion(void) const { return -1 ; }
  void    MakeBranch(Option_t* opt) ;
  Double_t MakeEnergy(const Double_t energy) ;                       // makes the detected energy    
  TVector3 MakePosition(const Double_t energy, const TVector3 pos, const Double_t th, const Double_t ph) ; 
                                                                     // makes the detected position
  void MakeRecParticle(const Int_t modid, const TVector3 pos, AliPHOSFastRecParticle & rp) ;  // makes a reconstructes particle from primary
  Int_t   MakeType(AliPHOSFastRecParticle & rp) ;                    // gets the detected type of particle
  FastRecParticlesList * FastRecParticles() { return fFastRecParticles ; } // gets TClonesArray of reconstructed particles
  virtual void ResetPoints() ; 
  void         ResetFastRecParticles() ; 
  void         SetBigBox(Int_t index, Float_t value) ;                             
  Double_t     SigmaE(Double_t energy) ;    // calulates the energy resolution at a given Energy                           
  Double_t     SigmaP(Double_t energy, Int_t inc) ; // calulates the position resolution at a given Energy at a given incidence                           
  virtual void StepManager(void) ;          // does the tracking through PHOS and a preliminary digitalization
  
private:
  
  Float_t fBigBoxX ;                         // main box containing all PHOS (EMC+PPSD)
  Float_t fBigBoxY ;                         // main box containing all PHOS (EMC+PPSD)
  Float_t fBigBoxZ ;                         // main box containing all PHOS (EMC+PPSD)
  FastRecParticlesList * fFastRecParticles ; // list of particles modified by the response function 
  AliPHOSGeometry * fGeom ;                  // geometry definition
  Int_t fNRecParticles ;                     // number of detected particles
  TRandom fRan ;                             // random number generator
  Double_t fResPara1 ;                       // parameter for the energy resolution dependence  
  Double_t fResPara2 ;                       // parameter for the energy resolution dependence  
  Double_t fResPara3 ;                       // parameter for the energy resolution dependence 
  Double_t fPosParaA0 ;                      // parameter for the position resolution
  Double_t fPosParaA1 ;                      // parameter for the position resolution 
  Double_t fPosParaB0 ;                      // parameter for the position resolution 
  Double_t fPosParaB1 ;                      // parameter for the position resolution 
  Double_t fPosParaB2 ;                      // parameter for the position resolution

  ClassDef(AliPHOSvFast,1)  // PHOS main class , version for fast simulation

};

#endif // AliPHOSVFAST_H
