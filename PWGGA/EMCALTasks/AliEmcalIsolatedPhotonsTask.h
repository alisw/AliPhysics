#ifndef ALIEMCALISOLATEDPHOTONSTASK_H
#define ALIEMCALISOLATEDPHOTONSTASK_H

// $Id: AliEmcalIsolatedPhotonsTask.h $

class TClonesArray;
class TString;
class AliESDtrackCuts;
class AliVTrack;
class AliVCluster;
class TList;
class TH1F;
class TH2F;
class AliEmcalJet;

#include "AliAnalysisTaskSE.h"

class AliEmcalIsolatedPhotonsTask : public AliAnalysisTaskSE {
 public:
  AliEmcalIsolatedPhotonsTask();
  AliEmcalIsolatedPhotonsTask(const char *name);
  virtual ~AliEmcalIsolatedPhotonsTask();

  void                        UserCreateOutputObjects();
  void                        UserExec(Option_t *option);
  void                        Terminate(Option_t *option);

  void                        SetClusName(const char *n)                    { fCaloName      = n          ; }
  void                        SetJetsName(const char *n)                    { fJetsName      = n          ; }
  void                        SetTracksName(const char *n)                  { fTracksName    = n          ; }
  virtual void                SetTrackCuts(AliESDtrackCuts *cuts)           { fESDTrackCuts = cuts        ; }
  virtual AliESDtrackCuts    *GetTrackCuts()                          const { return fESDTrackCuts        ; }
  virtual void                SetAODFilterBit(const Int_t b)                { fFilterBit = b              ; }
  virtual Int_t               GetAODFilterBit()                       const { return fFilterBit           ; }
  virtual void                SetSkimmedESD(const Bool_t s)                 { fSkimmedESD = s             ; }
  virtual Bool_t              GetSkimmedESD()                         const { return fSkimmedESD          ; }


 protected:

  AliVTrack                  *GetTrack(const Int_t i)          const;
  Int_t                       GetNumberOfTracks()              const;
  AliVCluster                *GetCaloCluster(const Int_t i)    const;
  Int_t                       GetNumberOfCaloClusters()        const;
  AliEmcalJet                *GetJet(const Int_t i)            const;
  Int_t                       GetNumberOfJets()                const;
  void                        FillHistograms()                      ;
  void                        RetrieveEventObjects()                ;
  Bool_t                      AcceptTrack(AliVTrack *track)         ;

  TList                      *fOutput;                 // Output list

  TString                     fTracksName;             // name of track collection
  TString                     fCaloName;               // name of calo cluster collection
  TString                     fJetsName;               // name of jet collection
  Bool_t                      fSkimmedESD;             // flag if skimmed ESD
  AliESDtrackCuts            *fESDTrackCuts;           // Track cuts
  Int_t                       fFilterBit;              // AOD filter bit
  TClonesArray               *fTracks;                 //!Tracks
  TClonesArray               *fCaloClusters;           //!Clusters
  TClonesArray               *fJets;                   //!Jets
  TH1F                       *fHistTracksPt;           // Pt spectrum of tracks
  TH1F                       *fHistClustersEnergy;     // Energy spectrum of clusters
  TH2F                       *fHistEPcorrelation;      // Energy-momentum correlation
  TH1F                       *fHistJetsEnergy;         // Energy spectrum of jets
  TH1F                       *fHistJetsNE;             // Jet neutral energy spectrum
  TH1F                       *fHistJetsNEF;            // Jet neutral energy fraction
  TH1F                       *fHistJetsZ;              // Constituent Pt over Jet E ratio
  TH2F                       *fHistTrPhiEta;           // Phi-Eta distribution of tracks
  TH2F                       *fHistClusPhiEta;         // Phi-Eta distribution of clusters
  TH2F                       *fHistJetPhiEta;          // Phi-Eta distribution of jets

  Int_t                       Ptbins;                  // No. of pt bins
  Float_t                     Ptlow;                   // Min pt
  Float_t                     Ptup;                    // Max pt
  Int_t                       Ebins;                   // No. of e bins
  Float_t                     Elow;                    // Min e
  Float_t                     Eup;                     // Max e

 private:
  AliEmcalIsolatedPhotonsTask(const AliEmcalIsolatedPhotonsTask&);            // not implemented
  AliEmcalIsolatedPhotonsTask &operator=(const AliEmcalIsolatedPhotonsTask&); // not implemented

  ClassDef(AliEmcalIsolatedPhotonsTask, 1) // Isolated photons task
};
#endif
