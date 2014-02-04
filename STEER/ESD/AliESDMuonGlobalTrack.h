#ifndef AliESDMuonGlobalTrack_H
#define AliESDMuonGlobalTrack_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//====================================================================================================================================================
//
//      ESD description of an ALICE muon forward track, combining the information of the Muon Spectrometer and the Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "TMath.h"
#include "TMatrixD.h"
#include "TDatabasePDG.h"
#include "TArrayI.h"
#include "TLorentzVector.h"

#include "AliVParticle.h"

class AliESDEvent;
class TClonesArray;

//====================================================================================================================================================

class AliESDMuonGlobalTrack : public AliVParticle {

public:

  AliESDMuonGlobalTrack();
  AliESDMuonGlobalTrack(Double_t px, Double_t py, Double_t pz);
  virtual ~AliESDMuonGlobalTrack() {;}
  AliESDMuonGlobalTrack(const AliESDMuonGlobalTrack& esdTrack);
  AliESDMuonGlobalTrack& operator=(const AliESDMuonGlobalTrack& esdTrack);
  virtual void Copy(TObject &obj) const;

  void  SetCharge(Int_t charge) { fCharge = charge; } 
  Short_t GetCharge() const { return fCharge; }

  /* Double_t GetOffset(Double_t x, Double_t y, Double_t z); */
  /* Double_t GetOffsetX(Double_t x, Double_t z); */
  /* Double_t GetOffsetY(Double_t y, Double_t z); */

  // Set and Get methods for kinematics at primary vertex
  void SetPxPyPz(Double_t px, Double_t py, Double_t pz);

  // Get and Set methods for global tracking info
  Double_t GetChi2(void) const {return fChi2;}        // chi2/ndf
  void     SetChi2(Double_t Chi2) {fChi2 = Chi2;}     // chi2/ndf
  
  // Get and Set methods for trigger matching
  void  SetMatchTrigger(Int_t matchTrigger) { fMatchTrigger = matchTrigger; }
  Int_t GetMatchTrigger() { return fMatchTrigger; }

  // Kinematics
  Double_t Pt()       const { return fPt;  }
  Double_t Eta()      const { return fEta; }
  Double_t Rapidity() const { return fRapidity; }
  Double_t Px()       const { return fPx; }
  Double_t Py()       const { return fPy; }
  Double_t Pz()       const { return fPz; }
  Double_t P()        const { return fP;  }

  Bool_t   PxPyPz(Double_t p[3]) const { p[0] = Px(); p[1] = Py(); p[2] = Pz(); return kTRUE; }

  // Additional methods to comply with AliVParticle
  Double_t Xv() const {return -999.;} // put reasonable values here
  Double_t Yv() const {return -999.;} //
  Double_t Zv() const {return -999.;} //
  Bool_t   XvYvZv(Double_t x[3]) const { x[0] = Xv(); x[1] = Yv(); x[2] = Zv(); return kTRUE; }  
  Double_t OneOverPt() const { return (Pt() != 0.) ? 1./Pt() : FLT_MAX; }
  Double_t Phi() const { return TMath::Pi()+TMath::ATan2(-Py(), -Px()); }
  Double_t Theta() const { return TMath::ATan2(Pt(), Pz()); }
  Double_t E() const { return TMath::Sqrt(M()*M() + P()*P()); }
  Double_t M() const { return TDatabasePDG::Instance()->GetParticle("mu-")->Mass(); }
  Double_t Y() const { return Rapidity(); }
  Short_t  Charge() const { return fCharge; }

  // Dummy
  const Double_t *PID() const { return (Double_t*)0x0; }
  Int_t PdgCode() const { return 0; }
  
  // Set the corresponding MC track number
  void  SetLabel(Int_t label) { fLabel = label; }
  // Return the corresponding MC track number
  Int_t GetLabel() const { return fLabel; }

  AliESDEvent* GetESDEvent() const { return fESDEvent; }
  void         SetESDEvent(AliESDEvent* evt) { fESDEvent = evt; }  
  
protected:

  Short_t fCharge, fMatchTrigger;

  // kinematics at vertex
  Double_t fPx, fPy, fPz, fPt, fP, fEta, fRapidity;

  // global tracking info
  Double_t fChi2;                 //  chi2 in the MUON+MFT track fit
  Double_t fChi2MatchTrigger;     //  chi2 of trigger/track matching

  Int_t fLabel;                   //  point to the corresponding MC track

  AliESDEvent *fESDEvent;         //! Pointer back to event to which the track belongs
  
  ClassDef(AliESDMuonGlobalTrack,1) // MUON+MFT ESD track class 

};

//====================================================================================================================================================

#endif 
