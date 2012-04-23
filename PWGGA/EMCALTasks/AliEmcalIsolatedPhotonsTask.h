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
  void                        SetTrgClusName(const char *n)                 { fTrgClusName   = n          ; }
  virtual void                SetTrackCuts(AliESDtrackCuts *cuts)           { fESDTrackCuts = cuts        ; }
  virtual AliESDtrackCuts    *GetTrackCuts()                          const { return fESDTrackCuts        ; }
  virtual void                SetAODFilterBit(const Int_t b)                { fFilterBit = b              ; }
  virtual Int_t               GetAODFilterBit()                       const { return fFilterBit           ; }

 protected:

  AliVTrack                  *GetTrack(const Int_t i)          const;
  Int_t                       GetNumberOfTracks()              const;
  AliVCluster                *GetCaloCluster(const Int_t i)    const;
  Int_t                       GetNumberOfCaloClusters()        const;
  AliEmcalJet                *GetJet(const Int_t i)            const;
  Int_t                       GetNumberOfJets()                const;
  AliVCluster                *GetTrgCluster(const Int_t i)     const;
  Int_t                       GetNumberOfTrgClusters()         const;
  void                        FillHistograms()                      ;
  void                        RetrieveEventObjects()                ;
  Bool_t                      AcceptTrack(AliVTrack *track)         ;

  TList                      *fOutput;                 // Output list

  TString                     fTracksName;             // name of track collection
  TString                     fCaloName;               // name of calo cluster collection
  TString                     fJetsName;               // name of jet collection
  TString                     fTrgClusName;            // name of trg clus name
  AliESDtrackCuts            *fESDTrackCuts;           // Track cuts
  Int_t                       fFilterBit;              // AOD filter bit
  TClonesArray               *fTracks;                 //!Tracks
  TClonesArray               *fCaloClusters;           //!Clusters
  TClonesArray               *fJets;                   //!Jets
  TClonesArray               *fTrgClusters;            //!Trg Clusters
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
  TH1F                       *fHistMaxTrgCluster;      // Energy distribution of max trigger clusters
  TH1F                       *fHistTrackPhi[3];        // Phi distribution of hybrid tracks
  TH1F                       *fHistTrackEta[3];        // Eta distribution of hybrid tracks

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
