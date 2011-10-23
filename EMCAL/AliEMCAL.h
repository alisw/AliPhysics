#ifndef ALIEMCAL_H
#define ALIEMCAL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */
/* History of cvs commits:
 *
 * $Log$
 * Revision 1.43  2007/03/10 22:19:01  pavlinov
 * move one varibels from AliEMCALv2 to AliEMCAL
 *
 * Revision 1.42  2007/02/24 20:42:35  pavlinov
 * fixed error of Geant3 parameters initialisation
 *
 * Revision 1.41  2007/02/05 10:43:25  hristov
 * Changes for correct initialization of Geant4 (Mihaela)
 *
 * Revision 1.40  2006/12/05 17:19:26  gustavo
 * Updated AliEMCAL::Digits2Raw, reads first provisional RCU mapping files to make Raw data with new AliCaloAltroMapping and AliCaloRawStream
 *
 *
 */
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
#include "AliEMCALTrigger.h" 
class AliEMCALRawUtils;
#include "AliReconstructor.h"
class AliEMCALTriggerData;

class AliEMCAL : public AliDetector {

 public:
  
  AliEMCAL(); 
  AliEMCAL(const char* name, const char* title="");

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
  
  virtual AliTriggerDetector* CreateTriggerDetector() const 
    { return new AliEMCALTrigger(); }

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
  AliEMCALTriggerData        * fTriggerData;      // Trigger parameters data container

private:
  AliEMCAL(const AliEMCAL& emcal);
  AliEMCAL & operator = (const AliEMCAL & /*rvalue*/);

  ClassDef(AliEMCAL,12) // Electromagnetic calorimeter (base class)
    
} ;

#endif // ALIEMCAL_H
