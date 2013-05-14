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
class AliAODEvent;
class AliESDEvent;
class AliPHOSGeoUtils;
class AliPHOSCalibData;
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
  void SetEventCuts(Double_t cuts[4])                   const { AliPHOSpPbPi0Header::SetSelectionCuts(cuts);   }

  static void SetRemovePileup(Bool_t rm=kFALSE)               { fgRemovePileup                   = rm;         }
  static void SetUseFiducialCut(Bool_t fc=kFALSE)             { fgUseFiducialCut                 = fc;         }
  static void SetDecaliWidth(Double_t width)                  { fgDecaliWidth                    = width;      }
  static void SetCaloClCuts(Double_t cuts[5])                 { for (Int_t i=5; i--; ) fgCuts[i] = cuts[i];    }          

 private:

  AliAnalysisTaskSEPHOSpPbPi0(const AliAnalysisTaskSEPHOSpPbPi0&);            // not implemented
  AliAnalysisTaskSEPHOSpPbPi0& operator=(const AliAnalysisTaskSEPHOSpPbPi0&); // not implemented

  void PHOSInitialize(AliESDEvent* const esd);
  void FillCaloClusterInfo(AliAODEvent* const aod, AliESDEvent* const esd);

  Bool_t IsGoodCaloCluster(Int_t iMod, Int_t cellX, Int_t cellZ);

  static const Float_t kLogWeight = 4.5;

  static Bool_t        fgRemovePileup;      // flag of remove pileup events
  static Bool_t        fgUseFiducialCut;    // flag of use fiducial cut
  static Double_t      fgDecaliWidth;       // decalibration width
  static Double_t      fgCuts[5];           // 0, min of cluster Energy
                                            // 1, min of NCells
                                            // 2, min of M02
                                            // 3, min of distance to bad channels
                                            // 4, max of the cluster TOF 

  Bool_t               fIsMC;               // flag of whether the input is MC
  TArrayF              fCentralityBin;      // Centrality bin
  TArrayI              fBufferSize;         // Buffer size for event mixing

  Int_t                fRunNumber;          // Run Number
  TH2I                *fPHOSBadMap[5];      // Container of PHOS bad channels map
  AliPHOSGeoUtils     *fPHOSGeo;            // PHOS geometry
  AliPHOSCalibData    *fPHOSCalibData;      // PHOS calibration

  TList               *fEventList[10][10];  // Event list for mixing
  TList               *fListQA;             // output list of QA histograms
  TList               *fListRD;             // output list of Real Data histograms
  TList               *fListMC;             // output list of MC histograms
  AliPHOSpPbPi0Header *fHeader;             // info at event level
  TClonesArray        *fCaloClArr;          // Container of Calo clusters Info
   
  ClassDef(AliAnalysisTaskSEPHOSpPbPi0, 1);
};

#endif
