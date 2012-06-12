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
  AliDimuInfoStoreRD(AliMuonInfoStoreRD* const trk0, AliMuonInfoStoreRD* const trk1, UInt_t selMask);
  AliDimuInfoStoreRD(const AliDimuInfoStoreRD &src);
  AliDimuInfoStoreRD& operator=(const AliDimuInfoStoreRD &src);
  virtual ~AliDimuInfoStoreRD();

  AliMuonInfoStoreRD* Muon(Int_t i) const { return (i<2 ? (AliMuonInfoStoreRD*)(fMuonRef[i].GetObject()) : 0x0); }

  TVector3 Momentum() const { return fMomentum; }
  Short_t  Charge()   const { return fCharge;   }
  Double_t InvM()     const { return fInvM;     }
  UInt_t   SelMask()  const { return fSelMask;  }
  Bool_t IsSelected(const UInt_t filter) { return ((fSelMask & filter) == filter); }

  static const char* StdBranchName() { return fgkStdBranchName.Data(); }

 protected:

  void FillDimuInfo();
  UInt_t fSelMask;   // dimuon selection mask
  TRef fMuonRef[2];  // ref to the two corresponding muon tracks

 private:

  static const TString fgkStdBranchName;  // Standard branch name

  TVector3 fMomentum;  // 3-momentum of dimuon
  Short_t  fCharge;    // charge of dimuon
  Double_t fInvM;      // invariance mass of dimuon

  ClassDef(AliDimuInfoStoreRD, 7);
};

#endif
