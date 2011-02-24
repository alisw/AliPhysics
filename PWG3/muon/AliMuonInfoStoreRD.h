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

class AliAODTrack;
class AliESDMuonTrack;

class AliMuonInfoStoreRD : public TObject {
 public:

  AliMuonInfoStoreRD();
  AliMuonInfoStoreRD(AliAODTrack *trk);
  AliMuonInfoStoreRD(AliESDMuonTrack *trk);
  AliMuonInfoStoreRD(const AliMuonInfoStoreRD &src);
  AliMuonInfoStoreRD& operator=(const AliMuonInfoStoreRD &src);
  virtual ~AliMuonInfoStoreRD();

  TVector3 Momentum()      const { return fMomentum; }

  void XYZAtDCA(Double_t dca[3]) const { for (Int_t i=3; i--;) dca[i]=fDCA[i]; }
  Double_t DCA() const  { return TMath::Sqrt(fDCA[0]*fDCA[0]+fDCA[1]*fDCA[1]); }

  Short_t  Charge()       const { return fCharge;           }
  Int_t    MatchTrigger() const { return fMatchTrigger;     }
  Double_t Chi2Tracker()  const { return fChi2FitMomentum;  }
  Double_t Chi2Trigger()  const { return fChi2MatchTrigger; }
  Double_t RabsEnd()      const { return fRabsEnd;          }

  Bool_t IsSelected();

  static const char* StdBranchName()              { return fgkStdBranchName.Data(); }
  static void SelectionCust(Double_t cuts[16])    { for (Int_t i=16; i--;) cuts[i]=fgCuts[i]; }
  static void SetSelectionCuts(Double_t cuts[16]) { for (Int_t i=16; i--;) fgCuts[i]=cuts[i]; }

 private:

  void FillMuonInfo(AliAODTrack *trk);
  void FillMuonInfo(AliESDMuonTrack *trk);

  void SetMomentum(Double_t p[3])      { fMomentum.SetXYZ(p[0],p[1],p[2]);      }

  void SetDCA(Double_t dca[3]) { for (Int_t i=3; i--;) fDCA[i]=dca[i];    }
  void SetCharge(Short_t charge)           { fCharge           = charge;  }
  void SetMatchTrigger(Int_t trigger)      { fMatchTrigger     = trigger; }
  void SetChi2FitMomentum(Double_t chi2)   { fChi2FitMomentum  = chi2;    }
  void SetChi2MatchTrigger(Double_t chi2)  { fChi2MatchTrigger = chi2;    }
  void SetRabsEnd(Double_t rAbsEnd)        { fRabsEnd          = rAbsEnd; }

  static const TString fgkStdBranchName;  // Standard branch name
  static Double_t fgCuts[16];  // 0, min of 3-momentum
                               // 1, max of 3-momentum
                               // 2, pt_Min
                               // 3, pt_Max
                               // 4, eta_Min
                               // 5, eta_Max
                               // 6, dca_Min
                               // 7, dca_Max
                               // 8, about trigger matching
                               // 9, about trigger matching
                               //10, rAbs_Min
                               //11, rAbs_Max
                               //12, chi2Tracker Min
                               //13, chi2Tracker Max
                               //14, chi2Trigger Min
                               //15, chi2Trigger Max

  TVector3 fMomentum;       // momentum corrected w vtx

  Double_t fDCA[3];            // distance of closet approach
  Short_t  fCharge;            // track charge
  Int_t    fMatchTrigger;      // type of match trigger
  Double_t fChi2FitMomentum;   // chi2/NDF of momentum fit
  Double_t fChi2MatchTrigger;  // chi2 of trigger matching
  Double_t fRabsEnd;  // position at the end of front absorber

  ClassDef(AliMuonInfoStoreRD, 5);
};

#endif
