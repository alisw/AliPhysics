#ifndef ALIEMCAL_H
#define ALIEMCAL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

///////////////////////////////////////////////////////////////////////////////
///
/// \class AliEMCAL
/// \ingroup EMCALbase
/// \brief Base Class for EMCAL description
/// 
/// This class contains material definitions    
/// for the EMCAL - It does not place the detector in ALICE (who does? AliEMCALv0?)
///
/// \author Yves Schutz (SUBATECH) 
/// \author Sahal Yacoob (LBNL/UCT)
/// \author Alexei Pavlinov (WSU) 
//////////////////////////////////////////////////////////////////////////////

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

class AliEMCAL : public AliDetector 
{

 public:
  
  AliEMCAL(); 
  AliEMCAL(const char* name, const char* title="", const Bool_t checkGeoAndRun = kTRUE);

  virtual ~AliEMCAL() ; 
  
  /// See in AliEMCALv2 
  virtual void  AddHit(Int_t, Int_t*, Float_t *) {
    Fatal("AddHit(Int_t, Int_t*, Float_t *", "not to be used: use AddHit( Int_t shunt, Int_t primary, Int_t track,Int_t id, Float_t *hits )") ;  
  }
  
  virtual AliDigitizer* CreateDigitizer(AliDigitizationInput* digInput) const;
  
  virtual void  CreateMaterials() ;   
  
  virtual void  Init() ;   
  
  virtual void  Digits2Raw();
  
  virtual void  FinishRun() {}                  
  
  virtual AliEMCALGeometry * GetGeometry() const ;
   // {return AliEMCALGeometry::GetInstance(GetTitle(),"") ;  }   
  
  virtual void  Hits2SDigits();
  
  virtual       Int_t   IsVersion(void) const = 0 ;   
  virtual const TString   Version()     const { return TString(" ") ; }   

   //  
  virtual AliLoader* MakeLoader(const char* topfoldername);

  virtual void  SetCheckRunNumberAndGeoVersion(Bool_t check) { fCheckRunNumberAndGeoVersion = check ; }

  Bool_t        Raw2SDigits(AliRawReader* rawReader);
  
protected:
  
  void          InitConstants();  

  Int_t                     fBirkC0;      ///<  Constant 0 for Birk's Law implementation
  Double_t                  fBirkC1;      ///<  Constant 1 for Birk's Law implementation
  Double_t                  fBirkC2;      ///<  Constant 2 for Birk's Law implementation

  Bool_t                    fCheckRunNumberAndGeoVersion; ///< Check if run number corresponds to the requested geometry and V1 is used
  
  AliEMCALGeometry        * fGeometry;    //!<! EMCal geometry access
  
  // For embedding
  static AliEMCALRawUtils * fgRawUtils;   ///<  Raw utilities class, for embedding 
  TClonesArray            * fTriggerData; ///<  Trigger parameters data container

private:
  
  AliEMCAL              (const AliEMCAL & emcal);
  AliEMCAL & operator = (const AliEMCAL & /*rvalue*/);

  /// \cond CLASSIMP
  ClassDef(AliEMCAL,13) ;
  /// \endcond

} ;

#endif // ALIEMCAL_H
