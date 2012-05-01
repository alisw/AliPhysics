#ifndef ALIEMCALQATASK_H
#define ALIEMCALQATASK_H

// $Id$

class TClonesArray;
class TString;
class AliVTrack;
class AliVCluster;
class TList;
class TH1F;
class TH2F;
class AliEmcalJet;

#include "AliAnalysisTaskSE.h"

class AliEmcalQATask : public AliAnalysisTaskSE {
 public:
  AliEmcalQATask();
  AliEmcalQATask(const char *name);
  virtual ~AliEmcalQATask();

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
  TH2F                       *fHistTracksCent;         // Number of tracks vs. centrality
  TH2F                       *fHistClusCent;           // Number of clusters vs. centrality
  TH1F                       *fHistTracksPt;           // Pt spectrum of tracks
  TH1F                       *fHistClustersEnergy;     // Energy spectrum of clusters
  TH2F                       *fHistEPcorrelation;      // Energy-momentum correlation
  TH1F                       *fHistJetsEnergy;         // Energy spectrum of jets
  TH2F                       *fHistTrPhiEta;           // Phi-Eta distribution of tracks
  TH2F                       *fHistClusPhiEta;         // Phi-Eta distribution of clusters
  TH2F                       *fHistJetPhiEta;          // Phi-Eta distribution of jets
  TH1F                       *fHistMaxTrgCluster;      // Energy distribution of max trigger clusters
  TH1F                       *fHistTrackPhi[5];        // Phi distribution of hybrid tracks
  TH1F                       *fHistTrackEta[5];        // Eta distribution of hybrid tracks

  Int_t                       Ptbins;                  // No. of pt bins
  Float_t                     Ptlow;                   // Min pt
  Float_t                     Ptup;                    // Max pt
  Int_t                       Ebins;                   // No. of e bins
  Float_t                     Elow;                    // Min e
  Float_t                     Eup;                     // Max e

 private:
  AliEmcalQATask(const AliEmcalQATask&);            // not implemented
  AliEmcalQATask &operator=(const AliEmcalQATask&); // not implemented

  ClassDef(AliEmcalQATask, 1) // Quality task for Emcal analysis
};
#endif
