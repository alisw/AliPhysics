#ifndef ALIEMCAL_H
#define ALIEMCAL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
//  Base Class for EMCAL     
//  holds all geant information of
//  materials, etc.
//                  
//*-- Author: Yves Schutz (SUBATECH) 

// --- ROOT system ---

class TString ;
class TFolder ;
class TRandom ; 
class TGraph;
class TF1;

// --- AliRoot header files ---
class AliRawReader;
#include "AliDetector.h"
#include "AliEMCALGeometry.h" 
#include "AliEMCALRawUtils.h"
#include "AliReconstructor.h"
class AliEMCALTriggerData;

class AliEMCAL : public AliDetector {

 public:
  
  AliEMCAL(); 
  AliEMCAL(const char* name, const char* title="", const Bool_t checkGeoAndRun = kTRUE);

  virtual ~AliEMCAL() ; 
  virtual void   AddHit(Int_t, Int_t*, Float_t *) {
    Fatal("AddHit(Int_t, Int_t*, Float_t *", "not to be used: use AddHit( Int_t shunt, Int_t primary, Int_t track,Int_t id, Float_t *hits )") ;  
  }
  virtual AliDigitizer* CreateDigitizer(AliDigitizationInput* digInput) const;
  virtual void  CreateMaterials() ;   
  virtual void  Init() ;   
  virtual void  Digits2Raw();
  
  virtual void  FinishRun() {}                  
  virtual AliEMCALGeometry * GetGeometry() const ;
   // {return AliEMCALGeometry::GetInstance(GetTitle(),"") ;  }   
  virtual void    Hits2SDigits();
  virtual Int_t   IsVersion(void) const = 0 ;   
  
   //  
  virtual AliLoader* MakeLoader(const char* topfoldername);
  virtual const TString Version() const {return TString(" ") ; }   

  virtual void  SetCheckRunNumberAndGeoVersion(Bool_t check) { fCheckRunNumberAndGeoVersion = check ; }

  Bool_t Raw2SDigits(AliRawReader* rawReader);
  
protected:
  void InitConstants();  //initializes some params

  Int_t    fBirkC0; // constants for Birk's Law implementation
  Double_t fBirkC1; // constants for Birk's Law implementation
  Double_t fBirkC2; // constants for Birk's Law implementation

  AliEMCALGeometry* fGeometry;              //!
  Bool_t   fCheckRunNumberAndGeoVersion;    // Check if run number corresponds to the requested geometry and V1 is used
  
  //For embedding
  static AliEMCALRawUtils    * fgRawUtils;        // raw utilities class, for embedding 
  TClonesArray               * fTriggerData;      // Trigger parameters data container

private:
  AliEMCAL(const AliEMCAL& emcal);
  AliEMCAL & operator = (const AliEMCAL & /*rvalue*/);

  ClassDef(AliEMCAL,13) // Electromagnetic calorimeter (base class)
    
} ;

#endif // ALIEMCAL_H
