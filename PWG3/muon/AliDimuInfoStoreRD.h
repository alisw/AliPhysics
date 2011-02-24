#ifndef ALIDIMUINFOSTORERD_H
#define ALIDIMUINFOSTORERD_H

/* Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//***********************************************************
// Class AliDimuInfoStoreRD
// class used to extract and store reco info of dimu candidate
// Author: X-M. Zhang, zhang@clermont.in2p3.fr
//                     zhangxm@iopp.ccnu.edu.cn
//***********************************************************

#include <TObject.h>
#include <TRef.h>
#include <TVector3.h>
#include <TString.h>

#include "AliMuonInfoStoreRD.h"

class AliDimuInfoStoreRD : public TObject {
 public:

  AliDimuInfoStoreRD();
  AliDimuInfoStoreRD(AliMuonInfoStoreRD* const trk0, AliMuonInfoStoreRD* const trk1);
  AliDimuInfoStoreRD(const AliDimuInfoStoreRD &src);
  AliDimuInfoStoreRD& operator=(const AliDimuInfoStoreRD &src);
  virtual ~AliDimuInfoStoreRD();

  AliMuonInfoStoreRD* Muon(Int_t i) const { return (i<2 ? (AliMuonInfoStoreRD*)(fMuonRef[i].GetObject()) : 0x0); }

  TVector3 Momentum() const { return fMomentum; }
  Short_t  Charge()   const { return fCharge;   }
  Double_t InvM()     const { return fInvM;     }

  Bool_t IsSelected();

  static const char* StdBranchName() { return fgkStdBranchName.Data(); }
  static void SetSelectionCuts(Double_t cuts[16]) { for (Int_t i=16; i--;) fgCutd[i]=cuts[i]; }

 protected:

  void FillDimuInfo();
  TRef fMuonRef[2];  // ref to the two corresponding muon tracks

 private:

  static const TString fgkStdBranchName;  // Standard branch name
  static Double_t fgCutd[16];             // single muon cuts for dimuon selection

  TVector3 fMomentum;  // 3-momentum of dimuon
  Short_t  fCharge;    // charge of dimuon
  Double_t fInvM;      // invariance mass of dimuon

  ClassDef(AliDimuInfoStoreRD, 5);
};

#endif
