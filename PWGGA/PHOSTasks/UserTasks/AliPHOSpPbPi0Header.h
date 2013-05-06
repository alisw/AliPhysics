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

class AliInputEventHander;
class AliMCEvent;
class AliStack;
class AliAODMCParticle;
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
  Int_t    VtxContrsN()            const { return fVtxContrsN;                     }
  Bool_t   IspAVertexOK()          const { return fIspAVertexOK;                   }
  Bool_t   IsPileupSPD()           const { return fIsPileupSPD;                    }
  Float_t  Centrality()            const { return fCentrality;                     }
  Double_t MagneticField()         const { return fMagneticField;                  }
  Bool_t   IsSelected();

  void SetEventInfo(AliInputEventHandler* const handler);

  void CreateHistograms(TList *listEvent, TList *listCaloCl, TList *listPi0, TList *listMC);
  void FillHistosEvent(TList *list);
  void FillHistosCaloCellsQA(TList *list, AliVCaloCells* const cells, AliPHOSGeoUtils* const phosGeo);
  void FillHistosCaloCluster(TList *list, TClonesArray* const caloClArr, Int_t cent);
  void FillHistosPi0(TList *list, TClonesArray* const caloClArr, Int_t cent);
  void FillHistosMixPi0(TList *list, TClonesArray* const caloClArr, TList* const eventlist, Int_t cent);
  void FillHistosMC(TList *list, AliMCEvent* const mcEvent, AliPHOSGeoUtils* const phosGeo, Int_t cent);
  void FillHistosMC(TList *list, AliStack*   const stack,   AliPHOSGeoUtils* const phosGeo, Int_t cent);

  static void SetIsMC(Bool_t isMC=kFALSE) { fgIsMC  = isMC;  }
  static void SetNCent(Int_t ncent=10)    { fgNCent = ncent; }

 private :

  void CreateHistosEvent(TList *list);
  void CreateHistosCaloCellsQA(TList *list);
  void CreateHistosCaloCluster(TList *list);
  void CreateHistosPi0(TList *list);
  void CreateHistosMixPi0(TList *list);
  void CreateHistosMC(TList *list);

  Bool_t IspAVertexOK(AliAODEvent* const event);
  Bool_t IspAVertexOK(AliESDEvent* const event);

  Int_t  HitPHOSModule(AliAODMCParticle* const pMC, AliPHOSGeoUtils* const phosGeo);
  Int_t  HitPHOSModule(TParticle* const pMC, AliPHOSGeoUtils* const phosGeo);
  Double_t PrimaryParticleWeight(Int_t pdg, Double_t pt);

  static Bool_t fgIsMC;            // flag to use MC
  static Int_t  fgNCent;           // # of centrality bins

  enum { kAll,      kCpv,       kDisp,       kFiducial,      kBoth,      kBothPlusFiducial, kPIDs                               };   // PID
  enum { kPtClu,    kEtaClu,    kPhiClu,     kM02Clu,        kM20Clu,    kTOFClu,           kNCellsClu, kNClustersClu, kVarsClu };   // clusters
  enum { kPtPi0,    kEtaPi0,    kPhiPi0,     kAsyPi0,        kAnglePi0,  kInvMassPi0,       kVarsPi0                            };   // pi0
  enum { kPtMixPi0, kEtaMixPi0, kPhiMixPi0,  kInvMassMixPi0, kVarsMixPi0                                                        };   // Mixed pi0
  enum { kVertexMC, kPtMC,      kRapidityMC, kPhiMC,         kWeightMC,  kVarsMC                                                };   // MC

  Double_t fVtx[3];                       // position of vtx
  TString  fFiredTriggerClass;            // trigger class
  UInt_t   fSelMask;                      // mask of physics selection
  Int_t    fVtxContrsN;                   // num. of contributors of vtx rec
  Bool_t   fIspAVertexOK;                 // is pA collision vertex OK
  Bool_t   fIsPileupSPD;                  // is Pileup from SPD
  Float_t  fCentrality;                   // event certrality
  Double_t fMagneticField;                // magnetic field

  ClassDef(AliPHOSpPbPi0Header, 1)
};

#endif
