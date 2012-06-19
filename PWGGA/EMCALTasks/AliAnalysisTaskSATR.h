#ifndef ALIANALYSISTASKSATR_H
#define ALIANALYSISTASKSATR_H

// $Id$

class TH1F;
class TH2F;
class TH1I;
class TH2I;
class TList;
class AliCaloCalibPedestal;
class AliAnalysisTaskEMCALClusterizeFast;
class AliEMCALGeometry;

#include "AliAnalysisTaskSE.h"

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

  Float_t                                fL0Calib                                    ; // L0 amplitude calibration  
  TString                                fCaloClustersName                           ; // Calo cluster collection name
  TString                                fTriggerClustersName                        ; // Trigger cluster collection name
  Float_t                                fMinCutL0Amp                                ; // Min L0 amplitude
  Float_t                                fMaxCutL0Amp                                ; // Max L0 amplitude
  Float_t                                fMinCutClusEnergy                           ; // Min cluster energy
  Float_t                                fMaxCutClusEnergy                           ; // Max cluster energy
  Bool_t                                 fTimeCutOn                                  ; // True = time cut on
  Int_t                                  fMinL0Time                                  ; // Min L0 time
  Int_t                                  fMaxL0Time                                  ; // Max L0 time
  Float_t                                fMinClusTime                                ; // Min clus time
  Float_t                                fMaxClusTime                                ; // Max clus time
  Bool_t                                 fCheckDeadClusters                          ; // True = check for dead clusters
  AliCaloCalibPedestal                  *fPedestal                                   ; // Calo calib pedestal object
  Bool_t                                 fLoadPed                                    ; // True = load pedesta
  TString                                fOCDBpath                                   ; // Path with OCDB location
  Float_t                                fMinDistanceFromBadTower                    ; // Min distance from bad tower
  AliAnalysisTaskEMCALClusterizeFast    *fClusterizer                                ; // Clusterizer
  AliAnalysisTaskEMCALClusterizeFast    *fTriggerClusterizer                         ; // Trigger clusterizer

  AliEMCALGeometry                      *fGeom                                       ; //!Pointer to emcal geometry object
  Int_t                                  fRun                                        ; //!Current run
  TList                                 *fOutput                                     ; //!Output list
  TH1F                                  *fHistEclus                                  ; //!Energy spectrum of clusters
  TH1F                                  *fHistEmaxClus                               ; //!Energy of max cluster per event
  TH2I                                  *fHistEtavsPhiMaxClus                        ; //!Position (eta-phi) of max cluster per event
  TH2F                                  *fHistEtavsEmaxClus                          ; //!Eta vs. energy of max cluster per event
  TH2F                                  *fHistPhivsEmaxClus                          ; //!Phi vs. energy of max cluster per event
  TH2F                                  *fHistTOFvsEclus                             ; //!TOF vs. energy of clusters
  TH2F                                  *fHistTOFvsEclusC                            ; //!Output histogram
  TH2F                                  *fHistNcellsvsEclus                          ; //!Output histogram
  
  TH1F                                  *fHistAmpTClus                               ; //!Output histogram
  TH1F                                  *fHistAmpMaxTClus                            ; //!Output histogram
  TH2I                                  *fHistEtavsPhiMaxTClus                       ; //!Output histogram
  
  TH2F                                  *fHistEmaxClusvsAmpMaxTClus                  ; //!Output histogram
  TH2F                                  *fHistEmaxClusvsAmpMatchedTClus              ; //!Output histogram
  TH1F                                  *fHistEmaxClusNotMatchingTClus               ; //!Output histogram
  TH2I                                  *fHistEtavsPhiMaxClusNotMatchingTClus        ; //!Output histogram
  TH2F                                  *fHistEmatchedClusvsAmpMaxTClus              ; //!Output histogram
  TH1F                                  *fHistAmpMaxTClusNotMatchingClus             ; //!Output histogram
  TH2I                                  *fHistEtavsPhiMaxTClusNotMatchingClus        ; //!Output histogram
  TH2I                                  *fHistIdxMaxClusvsIdxMaxTClus                ; //!Output histogram
  TH2I                                  *fHistPhiMaxClusvsPhiMaxTClus                ; //!Output histogram
  TH2I                                  *fHistEtaMaxClusvsEtaMaxTClus                ; //!Output histogram
  TH2F                                  *fHistTOFmaxClusvsTimeMaxTClus               ; //!Output histogram
  TH2F                                  *fHistEmatchedClusvsAmpMatchedTClus          ; //!Output histogram
  TH1F                                  *fHistEmatchedClus                           ; //!Output histogram
  TH1F                                  *fHistEmaxMatchedClus                        ; //!Output histogram
	
  TH1F                                  *fHistAmpL1TimeSum                           ; //!Output histogram
  TH1F                                  *fHistAmpMaxL1TimeSum                        ; //!Output histogram
  TH2F                                  *fHistAmpMaxL1TimeSumVScent                  ; //!Output histogram
  
  TH2F                                  *fHistAmpFastORvsAmpL1TimeSum                ; //!Output histogram
  
  TH1F                                  *fHistAmpFastOR                              ; //!Output histogram
  TH1F                                  *fHistAmpMaxFastOR                           ; //!Output histogram
  TH1F                                  *fHistTimeFastOR                             ; //!Output histogram
  TH2I                                  *fHistEtavsPhiFastOR                         ; //!Output histogram
  TH2I                                  *fHistEtavsPhiMaxFastOR                      ; //!Output histogram
  TH1F                                  *fHistTimeDispFastOR                         ; //!Output histogram
  TH2F                                  *fHistTimevsL0TimeFastOR                     ; //!Output histogram
  TH1I                                  *fHistNtimesFastOR                           ; //!Output histogram

  TH1F                                  *fHistEcells                                 ; //!Output histogram
  TH1F                                  *fHistEmaxCell                               ; //!Output histogram
  TH2F                                  *fHistTOFvsEcells                            ; //!Output histogram
  TH2F                                  *fHistTOFvsEcellsC                           ; //!Output histogram
  
  TH2F                                  *fHistEmaxCellvsAmpFastOR                    ; //!Output histogram
  
private:
  
  AliAnalysisTaskSATR (const AliAnalysisTaskSATR&);           // not implemented
  AliAnalysisTaskSATR operator=(const AliAnalysisTaskSATR&);  // not implemented
  
  ClassDef(AliAnalysisTaskSATR, 2);
};

#endif

