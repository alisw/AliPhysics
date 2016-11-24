#ifndef ALIPHOSPPBPI0HEADER_H
#define ALIPHOSPPBPI0HEADER_H

/* Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliPHOSpPbPi0Header
// class used to extract ,store info and fill histograms at event level
// Author: H-S. Zhu, hongsheng.zhu@cern.ch
//                   hszhu@iopp.ccnu.edu.cn
//*************************************************************************

#include <TNamed.h>
#include <TString.h>

class TList;
class TParticle;
class TClonesArray;

class AliVEventHandler;
class AliStack;
class AliVEvent;
class AliVCaloCells;
class AliPHOSGeoUtils;

class AliCaloClusterInfo;

class AliPHOSpPbPi0Header : public TNamed {
 public :

  AliPHOSpPbPi0Header();
  AliPHOSpPbPi0Header(const AliPHOSpPbPi0Header &src);
  AliPHOSpPbPi0Header& operator=(const AliPHOSpPbPi0Header &src);
  ~AliPHOSpPbPi0Header();

  void     GetXYZ(Double_t *vtx)   const { for (Int_t i=3; i--;)   vtx[i]=fVtx[i]; }
  Double_t Vx()                    const { return fVtx[0];                         }
  Double_t Vy()                    const { return fVtx[1];                         }
  Double_t Vz()                    const { return fVtx[2];                         }
  TString  FiredTriggerClass()     const { return fFiredTriggerClass;              }
  UInt_t   SelectionMask()         const { return fSelMask;                        }
  Bool_t   IsVertexOK()            const { return fIsVertexOK;                     }
  Bool_t   IsPileup()              const { return fIsPileup;                       }
  Float_t  Centrality()            const { return fCentrality;                     }

  Bool_t   IsSelected();

  void SetEventInfo(AliVEventHandler* const handler);

  void CreateHistograms(TList *listQA, TList *listRD, TList *listMC);
  void FillHistosEvent(TList *listQA);
  void FillHistosCaloCellsQA(TList *listQA, AliVCaloCells* const cells, AliPHOSGeoUtils* const phosGeo);
  void FillHistosCaloCluster(TList *listQA, TClonesArray* const caloClArr, Int_t cent);
  void FillHistosPi0(TList *listRD, TClonesArray* const caloClArr, Int_t cent);
  void FillHistosMixPi0(TList *listRD, TClonesArray* const caloClArr, TList* const eventlist, Int_t cent);
  void FillHistosMC(TList *listMC, AliStack* const stack, TClonesArray* const caloClArr, AliPHOSGeoUtils* const phosGeo, Int_t cent);

  static void SetIsMC(Bool_t isMC=kFALSE)         { fgIsMC                          = isMC;    }
  static void SetUseFiducialCut(Bool_t fc=kFALSE) { fgUseFiducialCut                = fc;      }
  static void SetNCent(Int_t ncent=10)            { fgNCent                         = ncent;   }
  static void SetSelectionCuts(Double_t cuts[3])  { for (Int_t i=3; i--;) fgCuts[i] = cuts[i]; }

 private :

  void CreateHistosEvent(TList *listQA);
  void CreateHistosCaloCellsQA(TList *listQA);
  void CreateHistosCaloCluster(TList *listQA);
  void CreateHistosPi0(TList *listRD);
  void CreateHistosMixPi0(TList *listRD);
  void CreateHistosMC(TList *listMC);

  Bool_t    CheckEventVertex(AliVEvent* const event);
  TString   ClassifyMCPi0(Int_t index, AliStack* const stack);
  Int_t     HitPHOSModule(TParticle* const pMC, AliPHOSGeoUtils* const phosGeo);

  static Bool_t   fgIsMC;           // flag to use MC
  static Bool_t   fgIspARun;        // flag to use pA vertex cut
  static Bool_t   fgUseFiducialCut; // flag to use fiducial cut
  static Int_t    fgNCent;          // # of centrality bins
  static Double_t fgCuts[3];        // 0, up limit of vz
                                    // 1, centrality max
                                    // 2, centrality min

  enum { kAll,      kCpv,        kDisp,       kBoth,          kCpv2,      kDisp2,      kBoth2,     kPIDs                   };   // PID
  enum { kPtClu,    kEtaClu,     kPhiClu,     kM02Clu,        kM20Clu,    kTOFClu,     kNCellsClu, kNClustersClu, kVarsClu };   // clusters
  enum { kPtPi0,    kEtaPi0,     kPhiPi0,     kAsyPi0,        kAnglePi0,  kInvMassPi0, kVarsPi0                            };   // pi0
  enum { kPtMixPi0, kEtaMixPi0,  kPhiMixPi0,  kInvMassMixPi0, kVarsMixPi0                                                  };   // Mixed pi0
  enum { kPtMC,     kRapidityMC, kRadiusMC,   kPhiMC,         kInvMassMC, kVarsMC                                          };   // MC

  Double_t fVtx[3];                 // position of vtx
  TString  fFiredTriggerClass;      // trigger class
  UInt_t   fSelMask;                // mask of physics selection
  Bool_t   fIsVertexOK;             // is vertex OK
  Bool_t   fIsPileup;               // is Pileup from SPD
  Float_t  fCentrality;             // event certrality

  ClassDef(AliPHOSpPbPi0Header, 3)
};

#endif
