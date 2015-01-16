#ifndef ALIMUONINFOSTORERD_H
#define ALIMUONINFOSTORERD_H

/* Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//***********************************************************
// Class AliMuonInfoStoreRD
// class used to extract and store reco info of muon track
// Author: X-M. Zhang, zhang@clermont.in2p3.fr
//                     zhangxm@iopp.ccnu.edu.cn
//***********************************************************

#include <TObject.h>
#include <TVector3.h>
#include <TString.h>
#include <TMath.h>

class AliAODTrack;
class AliESDMuonTrack;

class AliMuonInfoStoreRD : public TObject {
 public:

  AliMuonInfoStoreRD();
  AliMuonInfoStoreRD(AliAODTrack     *trk, UInt_t selMask);
  AliMuonInfoStoreRD(AliESDMuonTrack *trk, UInt_t selMask);
  AliMuonInfoStoreRD(const AliMuonInfoStoreRD &src);
  AliMuonInfoStoreRD& operator=(const AliMuonInfoStoreRD &src);
  virtual ~AliMuonInfoStoreRD();

  TVector3 MomentumAtVtx() const { return fMomentumAtVtx; }
  TVector3 MomentumAtDCA() const { return fMomentumAtDCA; }
  TVector3 MomentumUncor() const { return fMomentumUncor; }

  void XYZAtDCA(Double_t dca[3]) const { for (Int_t i=3; i--;) dca[i]=fDCA[i]; }
  Double_t DCA() const  { return TMath::Sqrt(fDCA[0]*fDCA[0]+fDCA[1]*fDCA[1]); }

  Short_t  Charge()       const { return fCharge;           }
  Int_t    MatchTrigger() const { return fMatchTrigger;     }
  Double_t Chi2Tracker()  const { return fChi2FitMomentum;  }
  Double_t Chi2Trigger()  const { return fChi2MatchTrigger; }
  Double_t RabsEnd()      const { return fRabsEnd;          }
  UInt_t   SelMask()      const { return fSelMask;          }
  Bool_t   IsSelected(const UInt_t filter) { return ((fSelMask & filter) == filter); }
  static const char* StdBranchName()       { return fgkStdBranchName.Data();       }

 private:

  void FillMuonInfo(AliAODTrack *trk);
  void FillMuonInfo(AliESDMuonTrack *trk);

  void SetMomentumAtVtx(Double_t p[3])    { fMomentumAtVtx.SetXYZ(p[0],p[1],p[2]); }
  void SetMomentumAtDCA(Double_t p[3])    { fMomentumAtDCA.SetXYZ(p[0],p[1],p[2]); }
  void SetMomentumUncor(Double_t p[3])    { fMomentumUncor.SetXYZ(p[0],p[1],p[2]); }
  void SetDCA(Double_t dca[3])            { for (Int_t i=3; i--;) fDCA[i]=dca[i];  }
  void SetCharge(Short_t charge)          { fCharge           = charge;  }
  void SetMatchTrigger(Int_t trigger)     { fMatchTrigger     = trigger; }
  void SetChi2FitMomentum(Double_t chi2)  { fChi2FitMomentum  = chi2;    }
  void SetChi2MatchTrigger(Double_t chi2) { fChi2MatchTrigger = chi2;    }
  void SetRabsEnd(Double_t rAbsEnd)       { fRabsEnd          = rAbsEnd; }

  static const TString fgkStdBranchName; // Standard branch name

  TVector3 fMomentumAtVtx; // momentum corrected w vtx
  TVector3 fMomentumAtDCA; // momentum at DCA in vtx plane
  TVector3 fMomentumUncor; // momentum at first station

  Double_t fDCA[3];           // distance of closet approach
  Short_t  fCharge;           // track charge
  Int_t    fMatchTrigger;     // type of match trigger
  Double_t fChi2FitMomentum;  // chi2/NDF of momentum fit
  Double_t fChi2MatchTrigger; // chi2 of trigger matching
  Double_t fRabsEnd;          // position at the end of front absorber
  UInt_t   fSelMask;          // mask of single muon selection

  ClassDef(AliMuonInfoStoreRD, 8);
};

#endif
