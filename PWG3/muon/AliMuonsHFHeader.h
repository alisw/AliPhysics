#ifndef ALIMUONSHFHEADER_H
#define ALIMUONSHFHEADER_H

/* Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//***********************************************************
// Class AliMuonsHFHeader
// class used to extract and store info at event level
// Author: X-M. Zhang, zhang@clermont.in2p3.fr
//                     zhangxm@iopp.ccnu.edu.cn
//***********************************************************

#include <TNamed.h>
#include <TString.h>

class TList;
class AliAODEvent;
class AliESDEvent;

class AliMuonInfoStoreRD;
class AliDimuInfoStoreRD;
class AliMuonInfoStoreMC;
class AliDimuInfoStoreMC;

class AliMuonsHFHeader : public TNamed {
 public :

  AliMuonsHFHeader();
  AliMuonsHFHeader(const AliMuonsHFHeader &src);
  AliMuonsHFHeader& operator=(const AliMuonsHFHeader &src);
  ~AliMuonsHFHeader();

  void GetXYZ(Double_t *vtx)  const { for (Int_t i=3; i--;) vtx[i]=fVtx[i]; }
  Double_t Vx()               const { return fVtx[0]; }
  Double_t Vy()               const { return fVtx[1]; }
  Double_t Vz()               const { return fVtx[2]; }
  Double_t Vt()               const { return TMath::Sqrt(fVtx[0]*fVtx[0] + fVtx[1]*fVtx[1]); }
  Int_t VtxContrsN()          const { return fVtxContrsN; }
  TString FiredTriggerClass() const { return fFiredTriggerClass; }
  Double_t Centrality()       const { return fCentrality; }

  void SetEvent(AliVVertex *vertex);
  void SetFiredTriggerClass(TString trigger) { fFiredTriggerClass=trigger; }
  Bool_t EventSelection();

  void CreateHistograms(TList *list);
  void FillHistosEvnH(TList *list);
  void FillHistosMuon(TList *list, AliMuonInfoStoreRD* infoStore, Int_t src=-1);
  void FillHistosDimu(TList *list, AliDimuInfoStoreRD* infoStore, Int_t src=-1);

  static const char* StdBranchName()             { return fgkStdBranchName.Data();          }
  static void SetAnaMode(Int_t anaMode=0)        { fgAnaMode=anaMode;                       }
  static void SetIsMC(Int_t isMC=kFALSE)         { fgIsMC   =isMC;                          }
  static void SetSelectionCuts(Double_t cuts[3]) { for (Int_t i=3; i--;) fgCuts[i]=cuts[i]; }

 private :

  void CreateHistosEvnH(TList *list);
  void CreateHistosMuon(TList *list, TString sName="");
  void CreateHistosDimu(TList *list, TString sName="");

  static const TString fgkStdBranchName;  // Standard branch name
  static Int_t  fgAnaMode;                // analysis mode
  static Bool_t fgIsMC;                   // flag to use MC
  static Double_t fgCuts[3];  // 0, low limit of num. of vtx contributors
                              // 1, up limit of vz
                              // 2, up limit of vt

  Double_t fVtx[3];   // position of vtx
  Int_t fVtxContrsN;  // num. of contributors of vtx rec

  TString fFiredTriggerClass; // trigger class

  Double_t fCentrality;  // event centrality class

  ClassDef(AliMuonsHFHeader, 3)
};

#endif
