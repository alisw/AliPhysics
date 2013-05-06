#ifndef ALIANALYSISTASKSEPHOSPPBPI0_H
#define ALIANALYSISTAKSSEPHOSPPBPI0_H

/* Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskSEPHOSpPbPi0
// AliAnalysisTaskSE for the gamma and pi0 from pPb collision analysis
// Author: H-S. Zhu, hongsheng.zhu@cern.ch
//                   hszhu@iopp.ccnu.edu.cn
//*************************************************************************

#include "AliAnalysisTaskSE.h"

class TH2I;
class TList;
class TString;
class TArray;
class TClonesArray;
class AliESDEvent;
class AliPHOSGeoUtils;
class AliPHOSpPbPi0Header;

class AliAnalysisTaskSEPHOSpPbPi0 : public AliAnalysisTaskSE {
 public:

  AliAnalysisTaskSEPHOSpPbPi0();
  AliAnalysisTaskSEPHOSpPbPi0(const char *name);
  virtual ~AliAnalysisTaskSEPHOSpPbPi0();
  
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  void SetPHOSBadMap(Int_t mod,TH2I *hMap)
  {
    if(fPHOSBadMap[mod]) delete fPHOSBadMap[mod];
    fPHOSBadMap[mod]=new TH2I(*hMap);
    printf("Set %s \n",fPHOSBadMap[mod]->GetName());
  }

  void SetUseMC(Bool_t isMC=kFALSE)                           { fIsMC          = isMC;                         }
  void SetXBins(const TArrayF& tCent, const TArrayI& tBuffer) { fCentralityBin = tCent; fBufferSize = tBuffer; }
  void SetLogWeight(Float_t logWeight)                  const { AliCaloClusterInfo::SetLogWeight(logWeight);   }

  static void SetMinNCells(Int_t ncells=2)                    { fgMinNCells                         = ncells;  }
  static void SetMinClusterEnergy(Double_t energy=0.3)        { fgMinClusterEnergy                  = energy;  }
  static void SetMinM02(Double_t m02=0.2)                     { fgMinM02                            = m02;     }
  static void SetMinDistToBad(Double_t dist=2.5)              { fgMinDistToBad                      = dist;    }

 private:

  AliAnalysisTaskSEPHOSpPbPi0(const AliAnalysisTaskSEPHOSpPbPi0&);            // not implemented
  AliAnalysisTaskSEPHOSpPbPi0& operator=(const AliAnalysisTaskSEPHOSpPbPi0&); // not implemented

  void PHOSInitialize(AliESDEvent* const esd);
  void FillCaloClusterInfo(AliESDEvent* const esd);

  Bool_t IsGoodCaloCluster(Int_t iMod, Int_t cellX, Int_t cellZ);

  static Int_t    fgMinNCells;
  static Double_t fgMinClusterEnergy;
  static Double_t fgMinM02;
  static Double_t fgMinDistToBad;

  Bool_t               fIsMC;               // flag of whether the input is MC
  TArrayF              fCentralityBin;      // Centrality bin
  TArrayI              fBufferSize;         // Buffer size for event mixing

  Int_t                fRunNumber;          // Run Number
  TH2I                *fPHOSBadMap[5];      // Container of PHOS bad channels map
  AliPHOSGeoUtils     *fPHOSGeo;            // PHOS geometry

  TList               *fEventList[10][10];  // Event list for mixing
  TList               *fListEvent;          // output list of Event
  TList               *fListCaloCl;         // output list of Calo Cluster histograms
  TList               *fListPi0;            // output list of Pi0 histograms
  TList               *fListMC;             // output list of MC histograms
  AliPHOSpPbPi0Header *fHeader;             // output for info at event level
  TClonesArray        *fCaloClArr;          // Container of Calo clusters Info
   
  ClassDef(AliAnalysisTaskSEPHOSpPbPi0, 1);
};

#endif
