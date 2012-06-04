#ifndef ALIANALYSISTASKEMCAL_H
#define ALIANALYSISTASKEMCAL_H

// $Id: AliAnalysisTaskEmcal.h 56756 2012-05-30 05:03:02Z loizides $

class TClonesArray;
class TString;
class TList;
class AliMCParticle;
class AliVTrack;
class AliVCluster;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEmcal : public AliAnalysisTaskSE {
 public:
  
  enum EmcalAnaType {
    kTPC       = 0,     // TPC only analysis
    kEMCAL     = 1,     // EMCal + TPC analysis
  };

  enum BeamType {
    kNA       = -1,
    kpp       = 0,
    kAA       = 1,
    kpA       = 2
  };

  AliAnalysisTaskEmcal();
  AliAnalysisTaskEmcal(const char *name);
  AliAnalysisTaskEmcal(const char *name, Bool_t histo); 
  virtual ~AliAnalysisTaskEmcal();

  virtual void                UserCreateOutputObjects();
  virtual void                UserExec(Option_t *option);
  virtual void                Terminate(Option_t *option);

  void                        SetTracksName(const char *n)                         { fTracksName     = n          ; }
  void                        SetClusName(const char *n)                           { fCaloName       = n          ; }
  void                        SetAnaType(EmcalAnaType type)                        { fAnaType        = type       ; }
  void                        SetPtCut(Float_t cut)                                { fPtCut          = cut        ; }
  void                        SetHistoBins(Int_t nbins, Float_t min, Float_t max)  { fNbins = nbins; fMinBinPt = min; fMaxBinPt = max; }

 protected:

  Bool_t                      AcceptTrack(AliVTrack* track, Bool_t acceptMC = kFALSE)              const;
  Bool_t                      AcceptCluster(AliVCluster* clus, Bool_t acceptMC = kFALSE)           const;
  BeamType                    GetBeamType()                                                             ;

  virtual Bool_t              RetrieveEventObjects();
  virtual Bool_t              Run()                        { return kTRUE ; }
  virtual Bool_t              FillHistograms()             { return kFALSE; }

  EmcalAnaType                fAnaType;                    // analysis type
  Bool_t                      fInitialized;                // whether or not the task has been already initialized
  Bool_t                      fCreateHisto;                // whether or not create histograms

  TString                     fTracksName;                 // name of track collection
  TString                     fCaloName;                   // name of calo cluster collection
  Int_t                       fNbins;                      // no. of pt bins
  Float_t                     fMinBinPt;                   // min pt in histograms
  Float_t                     fMaxBinPt;                   // max pt in histograms
  Float_t                     fPtCut;                      // cut on particle pt

  TClonesArray               *fTracks;                     //!tracks
  TClonesArray               *fCaloClusters;               //!clusters
  Float_t                     fCent;                       //!event centrality
  Int_t                       fCentBin;                    //!event centrality bin
  Double_t                    fVertex[3];                  //!event vertex
  BeamType                    fBeamType;                   //!event beam type

  TList                      *fOutput;                     //!output list

 private:
  AliAnalysisTaskEmcal(const AliAnalysisTaskEmcal&);            // not implemented
  AliAnalysisTaskEmcal &operator=(const AliAnalysisTaskEmcal&); // not implemented

  ClassDef(AliAnalysisTaskEmcal, 3) // EMCAL base analysis task
};
#endif
