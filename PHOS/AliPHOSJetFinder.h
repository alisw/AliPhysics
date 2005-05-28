#ifndef ALIPHOSJETFINDER_H
#define ALIPHOSJETFINDER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 */

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
  AliPHOSJetFinder(const AliPHOSJetFinder & jet) : TNamed(jet) {
    // copy ctor: no implementation yet
    Fatal("cpy ctor", "not implemented") ;
  }
  virtual ~AliPHOSJetFinder() ; // dtor

  void FindJetsFromParticles(const TClonesArray * plist,TObjArray * jetslist) ; //Do the job
  void FindJetsFromDigits(const TClonesArray * digits,TObjArray * jetslist) ; //Do the job

  void Print(const Option_t * = "") const ;

  void SetEtSeed(Double_t etseed){fEtSeed = etseed ;} ;
  void SetEtMin(Double_t etmin){fEtMin = etmin ;} ;
  void SetConRad(Double_t r){fConeRad = r ;} ;
  void SetMaxConeMove(Double_t move){fMaxConeMove=move ; } ;
  void SetMinConeMove(Double_t move){fMinConeMove=move ; } ;
  void SetStatusCode(Int_t stc = 1){fStatusCode=stc ;} ;
  AliPHOSJetFinder & operator = (const AliPHOSJetFinder & /*rvalue*/)  {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented") ; return *this ; 
  }
  
private:
  Double_t Calibrate(const AliPHOSDigit * digit) ;
  void    CalculateEEtaPhi(const AliPHOSDigit * d,Double_t &e, Double_t &Eta, Double_t &phi);

private:
  Int_t     fNJets ; //Number of jets
  Int_t     fStatusCode ; //Status code of particles to handle
  Int_t     fMode  ;   //Mode for background calculation

  Double_t  fConeRad ;   //Maximal radius of jet
  Double_t  fMaxConeMove ; //Maximal cone movement
  Double_t  fMinConeMove ; //Minimum cone movement
  Double_t  fEtSeed ;      //Transversal energy seed
  Double_t  fEtMin ;       //Minimal transversal energy
  Double_t  fPrecBg ;      //Precision due to background?  
  Double_t  fSimGain ;     //Simulated digit gain
  Double_t  fSimPedestal ; //Simulated digit pedestal


  TClonesArray * fParticles ; //Particles array
  TObjArray *    fJets ;      //Jets array

  ClassDef(AliPHOSJetFinder,1)  //Class to find Jets

};

#endif // AliPHOSJETFINDER_H
