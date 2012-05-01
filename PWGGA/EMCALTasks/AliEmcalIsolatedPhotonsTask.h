#ifndef ALIEMCALISOLATEDPHOTONSTASK_H
#define ALIEMCALISOLATEDPHOTONSTASK_H

// $Id: AliEmcalIsolatedPhotonsTask.h $

class TClonesArray;
class TString;
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
  Bool_t                      AcceptTrack(AliVTrack* /*track*/)     ;

  TList                      *fOutput;                 // Output list

  TString                     fTracksName;             // name of track collection
  TString                     fCaloName;               // name of calo cluster collection
  TString                     fJetsName;               // name of jet collection
  TString                     fTrgClusName;            // name of trg clus name
  TClonesArray               *fTracks;                 //!Tracks
  TClonesArray               *fCaloClusters;           //!Clusters
  TClonesArray               *fJets;                   //!Jets
  TClonesArray               *fTrgClusters;            //!Trg Clusters
  AliCentrality              *fCent;                   // Event centrality
  TH1F                       *fHistCentrality;         // Event centrality distribution
  TH1F                       *fHistJetsE[4];           // Jet energy spectrum
  TH1F                       *fHistJetsNE[4];          // Jet neutral energy spectrum
  TH1F                       *fHistJetsNEF[4];         // Jet neutral energy fraction
  TH1F                       *fHistJetsZ[4];           // Constituent Pt over Jet E ratio
  TH1F                       *fHistLeadingJetE[4];     // Leading jet energy spectrum
  TH1F                       *fHistTracksPtLJ[4];      // Pt spectrum of tracks
  TH1F                       *fHistClusELJ[4];         // Energy spectrum of clusters
  TH1F                       *fHistTracksPtBkg[4];     // Pt spectrum of tracks
  TH1F                       *fHistClusEBkg[4];        // Energy spectrum of clusters
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
