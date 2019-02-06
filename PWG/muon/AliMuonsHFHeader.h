#ifndef ALIMUONSHFHEADER_H
#define ALIMUONSHFHEADER_H

/* Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//***********************************************************
// Class AliMuonsHFHeader
// class used to extract and store info at event level
// Author: X-M. Zhang, zhang@clermont.in2p3.fr
//                     zhangxm@iopp.ccnu.edu.cn
//***********************************************************

#include <TNamed.h>
#include <TString.h>

class TList;
class AliMCEvent;
class AliVEventHandler;

class AliMuonInfoStoreRD;
class AliDimuInfoStoreRD;
class AliMuonInfoStoreMC;
class AliDimuInfoStoreMC;

class AliMuonsHFHeader : public TNamed {
 public:

  AliMuonsHFHeader();
  AliMuonsHFHeader(const AliMuonsHFHeader &src);
  AliMuonsHFHeader& operator=(const AliMuonsHFHeader &src);
  ~AliMuonsHFHeader();

  void       GetVMC(Double_t *vtx) const { for (Int_t i=3; i--;) vtx[i]=fVMC[i]; }
  void       GetXYZ(Double_t *vtx) const { for (Int_t i=3; i--;) vtx[i]=fVtx[i]; }
  Double_t   Vx()                  const { return fVtx[0]; }
  Double_t   Vy()                  const { return fVtx[1]; }
  Double_t   Vz()                  const { return fVtx[2]; }
  Double_t   Vt()                  const { return TMath::Sqrt(fVtx[0]*fVtx[0] + fVtx[1]*fVtx[1]); }
  Int_t      VtxContrsN()          const { return fVtxContrsN; }
  TString    FiredTriggerClass()   const { return fFiredTriggerClass; }
  UInt_t     SelectionMask()       const { return fSelMask; }
  Bool_t     IsMB()                const { return fIsMB; }
  Bool_t     IsMU()                const { return fIsMU; }
  Int_t      CentQA()              const { return fCentQA;}
  Double32_t EventPlane()          const { return fEventPlane; }
  Bool_t     IsSelected();

  void SetEventInfo(AliVEventHandler* const handler);

  void CreateHistograms(TList *list);
  void FillHistosEvnH(TList *list);
  void FillHistosMuon(TList *list, AliMuonInfoStoreRD* infoStore, Int_t src=0);
  void FillHistosDimu(TList *list, AliDimuInfoStoreRD* infoStore, Int_t src=0);

  static const char* StdBranchName()             { return fgkStdBranchName.Data();          }
  static void SetSelectionCuts(Double_t cuts[5]) { for (Int_t i=5; i--;) fgCuts[i]=cuts[i]; }

  void    SetAnaMode(Int_t anaMode)       { fAnaMode = anaMode;       }
  void    SetIsMC(Int_t isMC)             { fIsMC    = isMC;          }
  Int_t   GetAnaMode()              const { return fAnaMode;          }
  Int_t   GetNumOfTracklets()       const { return fNumOfTrklets;     }
  UInt_t  GetTrgInpts()             const { return fTrgInpts;         }
  Float_t GetImpParam()             const { return fImpParam;         }

  Bool_t  IsMC()                    const { return fIsMC;             }
  Bool_t  IsEvtInChunk()            const { return fIsEvtInChunk;     }
  Bool_t  IsVtxSeled2013pA()        const { return fIsVtxSeled2013pA; }

  enum { kV0M=0, kV0A, kV0C, kCL1, kZNA, kZNC };
  Float_t Centrality(Int_t centEstor=0);

  enum {
    kPUc1z1 = BIT(0),  // (minContributor=3, minZdis=0.5)
    kPUc1z2 = BIT(1),  // (3, 0.6)
    kPUc1z3 = BIT(2),  // (3, 0.8)
    kPUc1z4 = BIT(3),  // (3, 0.9)
    kPUc2z1 = BIT(4),  // (minContributor=4, minZdis=0.5)
    kPUc2z2 = BIT(5),  // (4, 0.6)
    kPUc2z3 = BIT(6),  // (4, 0.8)
    kPUc2z4 = BIT(7),  // (4, 0.9)
    kPUc3z1 = BIT(8),  // (minContributor=5, minZdis=0.5)
    kPUc3z2 = BIT(9),  // (5, 0.6)
    kPUc3z3 = BIT(10), // (5, 0.8)
    kPUc3z4 = BIT(11), // (5, 0.9)
    kPUc4z1 = BIT(12), // (minContributor=6, minZdis=0.5)
    kPUc4z2 = BIT(13), // (6, 0.6)
    kPUc4z3 = BIT(14), // (6, 0.8)
    kPUc4z4 = BIT(15)  // (6, 0.9)
  };
  UInt_t GetPUMask() const { return fPUMask; }

 private :

  void CreateHistosEvnH(TList *list, TString sName="");
  void CreateHistosMuon(TList *list, TString sName="");
  void CreateHistosDimu(TList *list, TString sName="");

  UInt_t CollectPUMask(AliVEvent *event); // collect the mask for the plie-up identification 

  static const TString fgkStdBranchName; // Standard branch name
  static Double_t fgCuts[5]; // 0, low limit of num. of vtx contributors
                             // 1, up limit of vz
                             // 2, up limit of vt
                             // 3, centrality max
                             // 4, centrality min

  Int_t      fAnaMode;           // analysis mode
  Bool_t     fIsMC;              // flag to use MC
  UInt_t     fSelMask;           // mask of physics selection
  Bool_t     fIsMB;              // is min. bias triggered event (for real data)
  Bool_t     fIsMU;              // is MUON triggered event (for real data)
  Double_t   fVtx[3];            // position of vtx
  Double_t   fVMC[3];            // position of vtx in MC
  Int_t      fVtxContrsN;        // num. of contributors of vtx rec
  TString    fFiredTriggerClass; // trigger class
  Int_t      fCentQA;            // quality of centrality determination
  Double32_t fEventPlane;        // event plane angle
  Bool_t     fIsEvtInChunk;      // flag denotes whether the event is in the chunk
  Bool_t     fIsVtxSeled2013pA;  // vertex selection based on 2013pA criteria
  Int_t      fNumOfTrklets;      // num. of tracklets
  UInt_t     fTrgInpts;          // trigger inputs (ID)
  Float_t    fImpParam;          // impact parameter
  
  Float_t fCentralityV0M; //centrality determinated via employed estimators
  Float_t fCentralityV0A;
  Float_t fCentralityV0C;
  Float_t fCentralityCL1;
  Float_t fCentralityZNA;
  Float_t fCentralityZNC;
  
  UInt_t fPUMask; // mask for tagging pile-up event

  ClassDef(AliMuonsHFHeader, 8)
};

#endif
