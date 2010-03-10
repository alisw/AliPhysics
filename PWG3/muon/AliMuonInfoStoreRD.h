#ifndef ALIMUONINFOSTORERD_H
#define ALIMUONINFOSTORERD_H

/* Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

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
  TVector3 MomentumAtDCA() const { return fMomentumAtDCA; }

  void XYZAtDCA(Double_t dca[3]) const { for (Int_t i=3; i--;) dca[i]=fDCA[i]; }
  Double_t DCA() const  { return TMath::Sqrt(fDCA[0]*fDCA[0]+fDCA[1]*fDCA[1]); }

  Short_t  Charge()           const { return fCharge;           }
  Int_t    MatchTrigger()     const { return fMatchTrigger;     }
  Int_t    NClusters()        const { return fNClusters;        }
  UInt_t   MUONClusterMap()   const { return fMUONClusterMap;   }
  Double_t Chi2FitMomentum()  const { return fChi2FitMomentum;  }
  Double_t Chi2MatchTrigger() const { return fChi2MatchTrigger; }

  Bool_t MuonSelection();

  static const char* StdBranchName()                  { return fgkStdBranchName.Data(); }
  static const void SelectionCust(Double_t cuts[10])  { for (Int_t i=10; i--;) cuts[i]=fgCuts[i]; }
  static void SetSelectionCuts(Double_t cuts[10]) { for (Int_t i=10; i--;) fgCuts[i]=cuts[i]; }

 private:

  void FillMuonInfo(AliAODTrack *trk);
  void FillMuonInfo(AliESDMuonTrack *trk);

  void SetMomentum(Double_t p[3])      { fMomentum.SetXYZ(p[0],p[1],p[2]);      }
  void SetMomentumAtDCA(Double_t p[3]) { fMomentumAtDCA.SetXYZ(p[0],p[1],p[2]); }

  void SetDCA(Double_t dca[3]) { for (Int_t i=3; i--;) fDCA[i]=dca[i]; }
  void SetCharge(Short_t charge)           { fCharge        = charge;  }
  void SetNClusters(Int_t ncls)            { fNClusters     = ncls;    }
  void SetMatchTrigger(Int_t trigger)      { fMatchTrigger  = trigger; }
  void SetMUONClusterMap(UInt_t clMap)     { fMUONClusterMap  = clMap; }
  void SetChi2FitMomentum(Double_t chi2)   { fChi2FitMomentum  = chi2; }
  void SetChi2MatchTrigger(Double_t chi2)  { fChi2MatchTrigger = chi2; }

  static const TString fgkStdBranchName;  // Standard branch name
  static Double_t fgCuts[10];  // 0, min of 3-momentum
                               // 1, max of 3-momentum
                               // 2, pt_Min
                               // 3, pt_Max
                               // 4, eta_Min
                               // 5, eta_Max
                               // 6, dca_Min
                               // 7, dca_Max
                               // 8, about trigger matching
                               // 9, about trigger matching

  TVector3 fMomentum;       // momentum corrected w vtx
  TVector3 fMomentumAtDCA;  // momentum at DCA in vtx plane

  Double_t fDCA[3];            // distance of closet approach
  Short_t  fCharge;            // track charge
  Int_t    fMatchTrigger;      // type of match trigger
  Int_t    fNClusters;         // number of clusters in the track
  UInt_t   fMUONClusterMap;    // map of MUON clusters
  Double_t fChi2FitMomentum;   // chi2/NDF of momentum fit
  Double_t fChi2MatchTrigger;  // chi2 of trigger matching

  ClassDef(AliMuonInfoStoreRD, 4);
};

#endif
