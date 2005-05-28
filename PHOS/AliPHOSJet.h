#ifndef ALIPHOSJET_H
#define ALIPHOSJET_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Base Class for Jets in ALICE
//                  
//*-- Author: D.Peressounko


// --- ROOT system ---
#include "TObject.h"
class TParticle ;
#include "TArrayI.h"

// --- Standard library ---

// --- AliRoot header files ---

class AliPHOSJet : public TObject {

public:
  AliPHOSJet() ;          // ctor
  AliPHOSJet(const AliPHOSJet & jet) : TObject(jet){
    // copy ctor: no implementation yet
    Fatal("cpy ctor", "not implemented") ;
  }
  virtual ~AliPHOSJet() ; 

  void AddDigit(Double_t e, Double_t eta, Double_t phi, Int_t index) ;
  void AddParticle(const TParticle * p, Int_t index) ;
  //adds particle p to jet. index: index of p in list of all particles in event

  Double_t DistanceToJet(const TParticle *p)const ;
  //calculates distance to Jet in accordance with some scheam: 
  //geometrical, inv mass, etc

  void CalculateAll(void) ;
  //calculate final Energy, Eta & phi from intermediate ones.

  const Int_t * Indexs(Int_t & nIndexs)const{nIndexs = fNpart; return fList->GetArray() ;}

  Bool_t IsInCone(const TParticle * p)const ;
  Bool_t IsInCone(Double_t eta, Double_t phi)const ;
  Bool_t AcceptConeDeviation(const TParticle *p)const ;
  Bool_t AcceptConeDeviation(Double_t e, Double_t eta, Double_t phi)const ;

  void SetConeRadius(Double_t r){fConeRad = r ;} ;
  void SetMaxConeMove(Double_t max = 0.15){fMaxConeMove = max ;} ;
  void SetMinConeMove(Double_t min = 0.05){fMinConeMove = min ;} ;

  Double_t Energy(void)const{if(fEnergy) return fEnergy ;
                        else return fSumEnergy ;}
  Double_t Eta(void)const{if(fEta) return fEta; 
                          else return fSumEta/fSumEnergy ;}
  Double_t Phi(void)const{if(fPhi) return fPhi; 
                          else return fSumPhi/fSumEnergy ;}
  Int_t GetNJetParticles(void)const{return fNpart;}

  void Print(const Option_t * = "") const ;
  AliPHOSJet & operator = (const AliPHOSJet & /*rvalue*/)  {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented") ; return *this ; 
  }

private:
  Int_t      fNpart ;     //Number of particles in jet
  TArrayI *  fList ;      //Indexes of particles in list

  Double_t   fConeRad ; // Cone radius 
  Double_t   fMaxConeMove ;// Maximum Cone movement
  Double_t   fMinConeMove ;// Minimum Cone movement 
  Double_t   fSumEnergy ; //! Intermediate energy
  Double_t   fSumEta ;    //! Intermediate eta
  Double_t   fSumPhi ;    //! Intermediate phi
  Double_t   fEnergy ;    //Energy of the jet
  Double_t   fEta ;       //Eta directtion of the jet
  Double_t   fPhi ;       //Phi direction of the jet
  Double_t   fLEnergy ;    //Energy of Leading particle of jet
  Double_t   fLEta ;       //Eta directtion of Leading particle of the jet
  Double_t   fLPhi ;       //Phi direction of leading particles of the jet

  ClassDef(AliPHOSJet,1)  // description 

};

//////////////////////////////////////////////////////

#endif // ALIPHOSJET_H





