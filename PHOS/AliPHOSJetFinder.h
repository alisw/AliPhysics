#ifndef ALIPHOSJETFINDER_H
#define ALIPHOSJETFINDER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
//                  
//*-- Author: D.Peressounko


// --- ROOT system ---
#include "TNamed.h"
class TClonesArray ;
class TObjArray ;
class AliPHOSDigit ;

// --- Standard library ---

// --- AliRoot header files ---

class AliPHOSJetFinder : public TNamed {

public:
  AliPHOSJetFinder() ;          // ctor

  virtual ~AliPHOSJetFinder() ; // dtor

  void FindJetsFromParticles(const TClonesArray * plist,TObjArray * jetslist) ; //Do the job
  void FindJetsFromDigits(const TClonesArray * digits,TObjArray * jetslist) ; //Do the job

  void Print(Option_t * option = "") ;

  void SetEtSeed(Double_t etseed){fEtSeed = etseed ;} ;
  void SetEtMin(Double_t etmin){fEtMin = etmin ;} ;
  void SetConRad(Double_t r){fConeRad = r ;} ;
  void SetMaxConeMove(Double_t move){fMaxConeMove=move ; } ;
  void SetMinConeMove(Double_t move){fMinConeMove=move ; } ;
  void SetStatusCode(Int_t stc = 1){fStatusCode=stc ;} ;
  
private:
  Double_t Calibrate(const AliPHOSDigit * digit) ;
  void    CalculateEEtaPhi(const AliPHOSDigit * d,Double_t &e, Double_t &Eta, Double_t &phi);

private:
  Int_t     fNJets ;
  Int_t     fStatusCode ; //Status code of particles to handle
  Int_t     fMode  ;   //Mode for background calculation

  Double_t  fConeRad ;   //Maximal radius of jet
  Double_t  fMaxConeMove ;
  Double_t  fMinConeMove ;
  Double_t  fEtSeed ;
  Double_t  fEtMin ;   
  Double_t  fPrecBg ;
  Double_t  fSimGain ;
  Double_t  fSimPedestal ;


  TClonesArray * fParticles ;
  TObjArray *    fJets ;

  ClassDef(AliPHOSJetFinder,1)  //Class to find Jets

};

#endif // AliPHOSJETFINDER_H
