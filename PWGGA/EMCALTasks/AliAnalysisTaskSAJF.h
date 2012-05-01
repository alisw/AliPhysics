#ifndef ALIANALYSISTASKSAJF_H
#define ALIANALYSISTASKSAJF_H

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

class AliAnalysisTaskSAJF : public AliAnalysisTaskSE {
 public:
  
  enum SAJFAnaType {
    kFullAcceptance  = 0,     // Full acceptance
    kEMCAL           = 1,     // EMCal acceptance only
    kEMCALFiducial   = 2      // EMCal fiduacial region only
  };

  AliAnalysisTaskSAJF();
  AliAnalysisTaskSAJF(const char *name);
  virtual ~AliAnalysisTaskSAJF();

  void                        UserCreateOutputObjects();
  void                        UserExec(Option_t *option);
  void                        Terminate(Option_t *option);

  void                        SetClusName(const char *n)                    { fCaloName      = n          ; }
  void                        SetJetsName(const char *n)                    { fJetsName      = n          ; }
  void                        SetTracksName(const char *n)                  { fTracksName    = n          ; }
  void                        SetTrgClusName(const char *n)                 { fTrgClusName   = n          ; }
  void                        SetAnaType(SAJFAnaType type)                  { fAnaType       = type       ; }

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
  Bool_t                      AcceptJet(AliEmcalJet* jet)           ;

  SAJFAnaType                 fAnaType;                // analysis type
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
  TH2F                       *fHistJetPhiEta;          // Phi-Eta distribution of jets
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
  AliAnalysisTaskSAJF(const AliAnalysisTaskSAJF&);            // not implemented
  AliAnalysisTaskSAJF &operator=(const AliAnalysisTaskSAJF&); // not implemented

  ClassDef(AliAnalysisTaskSAJF, 1) // Isolated photons task
};
#endif
