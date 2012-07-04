#ifndef ALIANALYSISTASKEMCAL_H
#define ALIANALYSISTASKEMCAL_H

// $Id: AliAnalysisTaskEmcal.h 56756 2012-05-30 05:03:02Z loizides $

class TClonesArray;
class TString;
class TList;
class AliEmcalParticle;
class AliMCParticle;
class AliVCluster;
class AliVTrack;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEmcal : public AliAnalysisTaskSE {
 public:
  
  enum EmcalAnaType {
    kTPC       = 0,     // TPC only analysis
    kEMCAL     = 1,     // EMCal + TPC analysis
    kTPCSmall  = 2,     // TPC only in EMCal acceptance
    kEMCALOnly = 3,     // EMCal only analysis
  };

  enum BeamType {
    kNA       = -1,
    kpp       = 0,
    kAA       = 1,
    kpA       = 2
  };

  AliAnalysisTaskEmcal();
  AliAnalysisTaskEmcal(const char *name, Bool_t histo=kFALSE); 
  virtual ~AliAnalysisTaskEmcal();

  void                        UserExec(Option_t *option);

  void                        SetAnaType(EmcalAnaType type)                         { fAnaType        = type        ; }
  void                        SetCentRange(Double_t min, Double_t max)              { fMinCent = min; fMaxCent = max; }
  void                        SetClusName(const char *n)                            { fCaloName       = n           ; }
  void                        SetClusPtCut(Double_t cut)                            { fClusPtCut      = cut         ; }
  void                        SetClusTimeCut(Double_t min, Double_t max)            { fClusTimeCutLow = min; fClusTimeCutUp = max;      }
  void                        SetHistoBins(Int_t nbins, Double_t min, Double_t max) { fNbins = nbins; fMinBinPt = min; fMaxBinPt = max; }
  void                        SetOffTrigger(UInt_t t)                               { fOffTrigger    = t;            }
  void                        SetPtCut(Double_t cut)                                { SetClusPtCut(cut); SetTrackPtCut(cut);            }
  void                        SetTrackPtCut(Double_t cut)                           { fTrackPtCut     = cut        ; }
  void                        SetTracksName(const char *n)                          { fTracksName     = n          ; }
  void                        SetVzRange(Double_t min, Double_t max)                { fMinVz = min; fMaxVz   = max;  }

 protected:
  Bool_t                      AcceptCluster(AliVCluster        *clus,  Bool_t acceptMC = kFALSE) const;
  Bool_t                      AcceptEmcalPart(AliEmcalParticle *part,  Bool_t acceptMC = kFALSE) const;
  Bool_t                      AcceptTrack(AliVTrack            *track, Bool_t acceptMC = kFALSE) const;
  virtual void                ExecOnce();
  virtual Bool_t              FillHistograms()                                     { return fCreateHisto; }
  BeamType                    GetBeamType();
  TClonesArray               *GetArrayFromEvent(const char *name, const char *clname=0);
  virtual Bool_t              IsEventSelected() const;
  virtual Bool_t              RetrieveEventObjects();
  virtual Bool_t              Run()                                                { return kTRUE                 ; }
  void                        SetInitialized(Bool_t ini = kTRUE)                   { fInitialized    = ini        ; }

  EmcalAnaType                fAnaType;                    // analysis type
  Bool_t                      fInitialized;                // whether or not the task has been already initialized
  Bool_t                      fCreateHisto;                // whether or not create histograms
  TString                     fTracksName;                 // name of track collection
  TString                     fCaloName;                   // name of calo cluster collection
  Double_t                    fMinCent;                    // min centrality for event selection
  Double_t                    fMaxCent;                    // max centrality for event selection
  Double_t                    fMinVz;                      // min vertex for event selection
  Double_t                    fMaxVz;                      // max vertex for event selection
  UInt_t                      fOffTrigger;                 // offline trigger for event selection
  Int_t                       fNbins;                      // no. of pt bins
  Double_t                    fMinBinPt;                   // min pt in histograms
  Double_t                    fMaxBinPt;                   // max pt in histograms
  Double_t                    fClusPtCut;                  // cut on cluster pt
  Double_t                    fTrackPtCut;                 // cut on track pt
  Double_t                    fClusTimeCutLow;             // low time cut for clusters
  Double_t                    fClusTimeCutUp;              // up time cut for clusters
  TClonesArray               *fTracks;                     //!tracks
  TClonesArray               *fCaloClusters;               //!clusters
  Double_t                    fCent;                       //!event centrality
  Int_t                       fCentBin;                    //!event centrality bin
  Double_t                    fEPV0;                       //!event plane V0
  Double_t                    fEPV0A;                      //!event plane V0A
  Double_t                    fEPV0C;                      //!event plane V0C
  Double_t                    fVertex[3];                  //!event vertex
  BeamType                    fBeamType;                   //!event beam type
  TList                      *fOutput;                     //!output list

 private:
  AliAnalysisTaskEmcal(const AliAnalysisTaskEmcal&);            // not implemented
  AliAnalysisTaskEmcal &operator=(const AliAnalysisTaskEmcal&); // not implemented

  ClassDef(AliAnalysisTaskEmcal, 4) // EMCAL base analysis task
};
#endif
