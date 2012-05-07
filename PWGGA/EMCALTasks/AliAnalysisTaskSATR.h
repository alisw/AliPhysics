/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#ifndef ALIANALYSISTASKSATR_H
#define ALIANALYSISTASKSATR_H

class TH1F;
class TH2F;
class TH1I;
class TH2I;
class TList;
class AliCaloCalibPedestal;
class AliAnalysisTaskEMCALClusterizeFast;

const Float_t   L0Calib     = .065;
const Int_t     L0Ampbins   = 100;
const Float_t   L0Amplow    = 0;
const Float_t   L0Ampup     = 100;
const Int_t     L0Timebins  = 20;
const Float_t   L0Timelow   = 0;
const Float_t   L0Timeup    = 20;
const Int_t     Ebins       = 100;
const Float_t   Elow        = 0;
const Float_t   Eup         = 25;
const Int_t     Ptbins      = 200;
const Float_t   Ptlow       = 0;
const Float_t   Ptup        = 100;
const Int_t     TOFbins     = 100;
const Float_t   TOFlow      = 0;
const Float_t   TOFup       = 1e-6;
const Float_t   L1Amplow    = 0;
const Float_t   L1Ampup     = 400;
const Int_t     L1Ampbins   = 400;
const Int_t     Indexesbins = 1440;
const Int_t     Indexeslow  = 0;
const Int_t     Indexesup   = 2880;
const Int_t     nPhibins    = 60;
const Int_t     nPhilow     = 0;
const Int_t     nPhiup      = 60;
const Int_t     nEtabins    = 48;
const Int_t     nEtalow     = 0;
const Int_t     nEtaup      = 48;
const Int_t     RowTrgbins  = 60;
const Int_t     RowTrglow   = 0;
const Int_t     RowTrgup    = 60;
const Int_t     ColTrgbins  = 48;
const Int_t     ColTrglow   = 0;
const Int_t     ColTrgup    = 48;

#include <AliAnalysisTaskSE.h>

class AliAnalysisTaskSATR : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskSATR();
  AliAnalysisTaskSATR(const char *name);
  virtual ~AliAnalysisTaskSATR();
  
  virtual void     UserCreateOutputObjects();
  virtual void     UserExec(Option_t *option);
  virtual void     Terminate(Option_t *);
  virtual void     Init();
  
  Bool_t                  GetCheckDeadClusters()                                      const { return fCheckDeadClusters;                        }
  AliCaloCalibPedestal*   GetPedestal()                                               const { return fPedestal;                                 }
  Bool_t                  GetTimeCutOn()                                              const { return fTimeCutOn;                                }
  
  void                    SetCheckDeadClusters(Bool_t c)                                    { fCheckDeadClusters = c;                           }
  void                    SetPedestal(AliCaloCalibPedestal *pds)                            { fPedestal = pds;                                  }
  void                    SetCaloClustersName(TString name)                                 { fCaloClustersName = name;                         }
  void                    SetTriggerClustersName(TString name)                              { fTriggerClustersName = name;                      }
  void                    SetTimeCutOn(Bool_t yes)                                          { fTimeCutOn = yes;                                 }
  void                    SetCutL0Amp(Float_t min = -1, Float_t max = 1000)                 { fMinCutL0Amp = min; fMaxCutL0Amp = max;           }
  void                    SetCutClusEnergy(Float_t min = -1, Float_t max = 100)             { fMinCutClusEnergy = min; fMaxCutClusEnergy = max; }
  void                    SetMinDistanceFromBadTower(Float_t d)                             { fMinDistanceFromBadTower = d;                     }
  void                    SetClusterizer(AliAnalysisTaskEMCALClusterizeFast *c)             { fClusterizer = c;                                 }
  void                    SetTriggerClusterizer(AliAnalysisTaskEMCALClusterizeFast *c)      { fTriggerClusterizer = c;                          }
  void                    SetLoadPed(Bool_t yes)                                            { fLoadPed = yes;                                   }
  void                    SetOCDBPath(const TString &path)                                  { fOCDBpath = path;                                 }
  void                    SetL0TimeCut(Int_t min, Int_t max)                                { fMinL0Time = min; fMaxL0Time = max;               }
  void                    SetClusTimeCut(Float_t min, Float_t max)                          { fMinClusTime = min; fMaxClusTime = max;           }
  
 protected:
  AliVCluster*                           GetClusterFromId(TClonesArray *caloClusters, Int_t id);
  
  TList                                 *fOutput;                         // Output list
  TH1F                                  *fHistEclus;                      // Energy spectrum of clusters
  TH1F                                  *fHistEmaxClus;
  TH2I                                  *fHistEtavsPhiMaxClus;
  TH2F                                  *fHistEtavsEmaxClus;
  TH2F                                  *fHistPhivsEmaxClus;
  TH2F                                  *fHistTOFvsEclus;
  TH2F                                  *fHistTOFvsEclusC;
  TH2F                                  *fHistNcellsvsEclus; 
  
  TH1F                                  *fHistAmpTClus;
  TH1F                                  *fHistAmpMaxTClus;
  TH2I                                  *fHistEtavsPhiMaxTClus;
  
  TH2F                                  *fHistEmaxClusvsAmpMaxTClus;
  TH2F                                  *fHistEmaxClusvsAmpMatchedTClus;
  TH1F                                  *fHistEmaxClusNotMatchingTClus;
  TH2I                                  *fHistEtavsPhiMaxClusNotMatchingTClus;
  TH2F                                  *fHistEmatchedClusvsAmpMaxTClus;
  TH1F                                  *fHistAmpMaxTClusNotMatchingClus;
  TH2I                                  *fHistEtavsPhiMaxTClusNotMatchingClus;
  TH2I                                  *fHistIdxMaxClusvsIdxMaxTClus;
  TH2I                                  *fHistPhiMaxClusvsPhiMaxTClus;
  TH2I                                  *fHistEtaMaxClusvsEtaMaxTClus;
  TH2F                                  *fHistTOFmaxClusvsTimeMaxTClus;
  TH2F                                  *fHistEmatchedClusvsAmpMatchedTClus;
  TH1F                                  *fHistEmatchedClus;
  TH1F                                  *fHistEmaxMatchedClus;
	
  TH1F                                  *fHistAmpL1TimeSum;
  TH1F                                  *fHistAmpMaxL1TimeSum;
  TH2F                                  *fHistAmpMaxL1TimeSumVScent;
  
  TH2F                                  *fHistAmpFastORvsAmpL1TimeSum;
  
  TH1F                                  *fHistAmpFastOR;
  TH1F                                  *fHistAmpMaxFastOR;
  TH1F                                  *fHistTimeFastOR;
  TH2I                                  *fHistEtavsPhiFastOR;
  TH2I                                  *fHistEtavsPhiMaxFastOR;
  TH1F                                  *fHistTimeDispFastOR;
  TH2F                                  *fHistTimevsL0TimeFastOR;
  TH1I                                  *fHistNtimesFastOR;

  TH1F                                  *fHistEcells;
  TH1F                                  *fHistEmaxCell;
  TH2F                                  *fHistTOFvsEcells;
  TH2F                                  *fHistTOFvsEcellsC;
  
  TH2F                                  *fHistEmaxCellvsAmpFastOR;

  TString                                fCaloClustersName;
  TString                                fTriggerClustersName;
  Float_t                                fMinCutL0Amp;
  Float_t                                fMaxCutL0Amp;
  Float_t                                fMinCutClusEnergy;
  Float_t                                fMaxCutClusEnergy;
  Bool_t                                 fTimeCutOn;
  Int_t                                  fMinL0Time;
  Int_t                                  fMaxL0Time;
  Float_t                                fMinClusTime;
  Float_t                                fMaxClusTime;
  Bool_t                                 fCheckDeadClusters;
  AliCaloCalibPedestal                  *fPedestal;
  Bool_t                                 fLoadPed;
  TString                                fOCDBpath;                               // path with OCDB location
  Float_t                                fMinDistanceFromBadTower;
  AliAnalysisTaskEMCALClusterizeFast    *fClusterizer;
  AliAnalysisTaskEMCALClusterizeFast    *fTriggerClusterizer;
  Int_t                                  fRun;
  
private:
  
  AliAnalysisTaskSATR (const AliAnalysisTaskSATR&);           // not implemented
  AliAnalysisTaskSATR operator=(const AliAnalysisTaskSATR&);  // not implemented
  
  ClassDef(AliAnalysisTaskSATR, 1);
};

#endif

