#ifndef ALICFPI0_H
#define ALICFPI0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: $ */
 
//_________________________________________________________________________
// Class to collect neutral pion kinematics and information
//
//-- Author: Paul Batzing 
//Class contains minimal information for collecting trigger pi0s.
//

#include "AliVParticle.h"
#include "TLorentzVector.h"


class AliCFPI0 :public AliVParticle{
 
  public:
    
    AliCFPI0() ;
    AliCFPI0(TLorentzVector V);
    ~AliCFPI0(){}
    
    
    Double_t Px() const {return fPx;}
    Double_t Py() const {return fPy;}
    Double_t Pz() const {return fPz;}
    Double_t Pt() const {return fPt;}
    Double_t P()  const {return fP;}

    Double_t E() const {return fE;}
    Double_t Phi() const {return fPhi;}
    Double_t Theta() const {return fTheta;}
    Double_t Eta() const {return fEta;}
    
    Bool_t IsPHOS() const {return fIsPhos ;}
    Bool_t IsEMCAL() const {return fIsEmcal ;}
    
    void SetPx(Double_t lPx) {fPx = lPx;}
    void SetPy(Double_t lPy) {fPy = lPy;}
    void SetPz(Double_t lPz) {fPz = lPz;}
    void SetPt(Double_t lPt) {fPt = lPt;}
    void SetE(Double_t lE) {fE = lE;}
    void SetPhi(Double_t lPhi) {fPhi = lPhi;}
    void SetTheta(Double_t lTheta) {fTheta = lTheta;}
    void SetEta(Double_t lEta ) {fEta = lEta;}
    
    void SetIsPhos(Bool_t lIsPhos) {fIsPhos = lIsPhos;}
    void SetIsEmcal(Bool_t lIsEmcal) {fIsEmcal = lIsEmcal;}
    
    
    
  //do I need these? All virtual functions must be implementet?
  Bool_t   PxPyPz(Double_t p[3]) const {return p[0]==fPx&&p[1]==fPy&&p[2]==fPz; }
  Double_t Xv() const {return 0; }
  Double_t Yv() const {return 0; }
  Double_t Zv() const {return 0;}
  Bool_t   XvYvZv(Double_t x[3]) const {return x[0]==0&&x[1]==0&&x[2]==0;}  
  Double_t OneOverPt()  const {return 0;}
  Double_t M()          const {return 0;}
  Double_t Y()          const {return 0;}
  Short_t Charge()      const {return 0;}
  Int_t   GetLabel()    const {return 0;}
  // PID
  Int_t   PdgCode()     const {return 0;}
  const Double_t *PID() const {return 0;} // return PID object (to be defined, still)
    
    
  private:
  Double_t fPx;
  Double_t fPy;
  Double_t fPz;
  Double_t fPt;
  Double_t fP;
  Double_t fPhi;
  Double_t fTheta;
  Double_t fEta;
  Double_t fE;

  
  Bool_t fIsPhos;
  Bool_t fIsEmcal;

   ClassDef(AliCFPI0, 1);

  
};

#endif
