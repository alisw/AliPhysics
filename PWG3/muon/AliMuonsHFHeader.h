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

  ULong64_t TriggerMask()     const { return fTriggerMask;        }
  TString FiredTrigger()      const { return fFiredTrigger;       }
  Int_t  NFiredTrigger()      const { return fNFiredTrigger;      }
  Int_t  EventType()          const { return fEventType;          }
  Bool_t IsPhysicsTriggered() const { return fIsPhysicsTriggered; }
  Bool_t IsPhysicsAccepted()  const { return fIsPhysicsAccepted;  }
  Bool_t IsTriggerFired(TString trigger);

  void GetXYZ(Double_t *vtx) const { for (Int_t i=3; i--;) vtx[i]=fVtx[i]; }
  Double_t Vx() const { return fVtx[0]; }
  Double_t Vy() const { return fVtx[1]; }
  Double_t Vz() const { return fVtx[2]; }
  Double_t Vt() const { return TMath::Sqrt(fVtx[0]*fVtx[0] + fVtx[1]*fVtx[1]); }
  Bool_t  IsUnrecoVertex() const { return fUnrecoVertex;  }
  Int_t NVtxContributors() const { return fNContributors; }

  void GetXYZSPD(Double_t *vtx) const { for (Int_t i=3; i--;) vtx[i]=fVtxSPD[i]; }
  Double_t VxSPD() const { return fVtxSPD[0]; }
  Double_t VySPD() const { return fVtxSPD[1]; }
  Double_t VzSPD() const { return fVtxSPD[2]; }
  Double_t VtSPD() const { return TMath::Sqrt(fVtxSPD[0]*fVtxSPD[0] + fVtxSPD[1]*fVtxSPD[1]); }
  Bool_t  IsUnrecoVtxSPD()    const { return fUnrecoVtxSPD;     }
  Int_t NVtxContributorsSPD() const { return fNContributorsSPD; }
  Int_t NTrackletsSPD()       const { return fNTrackletsSPD;    }

  Double_t Centrality() const { return fCentrality; }

  void SetEvent(AliAODEvent *event);
  void SetEvent(AliESDEvent *event);

  void EventSelection(TString triggerName);
  void CreateHistograms(TList *listEvent=0, TList *listMuon=0, TList *listDimu=0);
  void CreateHistosEventH(TList *list);
  void CreateHistosMuonRD(TList *list);
  void CreateHistosDimuRD(TList *list);
  void CreateHistosMuonMC(TList *list);
  void CreateHistosDimuMC(TList *list);

  void FillHistosEventH(TList *list);
  void FillHistosMuonRD(TList *list, AliMuonInfoStoreRD* const muonStoreRD);
  void FillHistosDimuRD(TList *list, AliDimuInfoStoreRD* const dimuStoreRD);
  void FillHistosMuonMC(TList *list, AliMuonInfoStoreMC* const muonStoreMC);
  void FillHistosDimuMC(TList *list, AliDimuInfoStoreMC* const dimuStoreMC);

  static const char* StdBranchName()    { return fgkStdBranchName.Data(); }
  static const Bool_t IsEventSelected() { return fgIsEventSelected; }
  static void SetAnaMode(Int_t anaMode=0)        { fgAnaMode=anaMode; }
  static void SetIsMC(Int_t isMC=kFALSE)         { fgIsMC   =isMC;    }
  static void SetSelectionCuts(Double_t cuts[3]) { for (Int_t i=3; i--;) fgCuts[i]=cuts[i]; }

 private :

  void SetFiredTrigger(TString str);
  void PhysicsTriggerAna(const AliESDEvent *esd);
  void EventSelection();

  static const TString fgkStdBranchName;  // Standard branch name
  static Bool_t fgIsEventSelected;        // flag for event selection
  static Int_t  fgAnaMode;                // analysis mode
  static Bool_t fgIsMC;                   // flag to use MC
  static Double_t fgCuts[3];  // 0, low limit of num. of vtx contributors
                              // 1, up limit of vz
                              // 2, up limit of vt

  ULong64_t fTriggerMask;      // trigger mask
  TString fFiredTrigger;       // fired of trigger class of event
  Int_t fNFiredTrigger;        // num. of fired trigger class
  Bool_t fIsPhysicsTriggered;  // flag of final physics trigger from AliPhysicsSelection
  Bool_t fIsPhysicsAccepted;   // flag of physiscs selection w/ BKG Id
  Int_t fEventType;            // event type

  Double32_t fVtx[3];      // position of vtx
  Bool_t fUnrecoVertex;    // flag for unreco vtx
  Int_t fNContributors;    // num. of contributors of vtx rec

  Double32_t fVtxSPD[3];   // position of vtx
  Bool_t fUnrecoVtxSPD;    // flag for unreco vtx
  Int_t fNContributorsSPD; // num. of contributors of vtx rec
  Int_t fNTrackletsSPD;  // num. of SPD tracklets

  Double_t fCentrality;    // event centrality class

  ClassDef(AliMuonsHFHeader, 2)
};

#endif
