#ifndef ALIANALYSISTASKSAJF_H
#define ALIANALYSISTASKSAJF_H

// $Id$

class TClonesArray;
class TString;
class TH1F;
class TH2F;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskSAJF : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskSAJF();
  AliAnalysisTaskSAJF(const char *name);
  virtual ~AliAnalysisTaskSAJF();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);
  void                        Init();

  void                        SetJetMinRC2LJ(Float_t d)                            { fMinRC2LJ                = d          ; } 
  void                        SetEmbJetsName(const char *n)                        { fEmbJetsName             = n          ; }
  void                        SetRandTracksName(const char *n)                     { fRandTracksName          = n          ; }
  void                        SetRandClusName(const char *n)                       { fRandCaloName            = n          ; }
  void                        SetRhoName(const char *n)                            { fRhoName                 = n          ; }
  void                        SetSkipEventsWithNoJets(Bool_t s)                    { fSkipEventsWithNoJets    = s          ; } 

 protected:

  Bool_t                      RetrieveEventObjects()                                                                        ;
  Bool_t                      FillHistograms()                                                                              ;
  void                        GetLeadingJets(Int_t &maxJetIndex, Int_t &max2JetIndex)                                       ;
  void                        DoJetLoop()                                                                                   ;
  void                        DoEmbJetLoop(AliEmcalJet* &embJet, TObject* &maxPart)                                         ;
  void                        DoTrackLoop(Int_t maxJetIndex)                                                                ;
  void                        DoClusterLoop(Int_t maxJetIndex)                                                              ;
  void                        GetRigidCone(Float_t &pt, Float_t &eta, Float_t &phi, Bool_t acceptMC = kFALSE, 
					   AliEmcalJet *jet = 0, TClonesArray* tracks = 0, TClonesArray* clusters = 0) const;

  Float_t                     fMinRC2LJ;                   // Minimum distance random cone to leading jet
  TString                     fEmbJetsName;                // Name of embedded jets collection
  TString                     fRandTracksName;             // Name of randomized track collection
  TString                     fRandCaloName;               // Name of randomized calo cluster collection
  TString                     fRhoName;                    // Name of rho object
  Bool_t                      fSkipEventsWithNoJets;       // Whether or not skip events with no jets

  TClonesArray               *fEmbJets;                    //!Embedded Jets
  TClonesArray               *fRandTracks;                 //!Randomized tracks
  TClonesArray               *fRandCaloClusters;           //!Randomized clusters
  Float_t                     fRho;                        //!Event rho

  // General histograms
  TH1F                       *fHistCentrality;             //!Event centrality distribution
  TH1F                       *fHistRejectedEvents;         //!Rejected events
  TH2F                       *fHistJetPhiEta[4];           //!Phi-Eta distribution of jets
  TH1F                       *fHistJetsPt[4];              //!Inclusive jet pt spectrum
  TH2F                       *fHistJetsPtArea[4];          //!Jet pt vs. area
  TH1F                       *fHistLeadingJetPt[4];        //!Leading jet pt spectrum
  TH1F                       *fHist2LeadingJetPt[4];       //!Second leading jet pt spectrum
  TH2F                       *fHistJetsNEFvsPt[4];         //!Jet neutral energy fraction vs. jet pt
  TH2F                       *fHistJetsZvsPt[4];           //!Constituent Pt over Jet Pt ratio vs. jet pt

  // Rho
  TH1F                       *fHistRho[4];                 //!Rho distribution
  TH2F                       *fHistRhoVSleadJetPt;         //!Area(leadjetarea) * rho vs. leading jet pt
  TH1F                       *fHistCorrJetsPt[4];          //!Corrected inclusive jet pt spectrum
  TH1F                       *fHistCorrLeadingJetPt[4];    //!Corrected leading jet pt spectrum

  // Random cones
  TH2F                       *fHistRCPhiEta;               //!Phi-Eta distribution of random cones
  TH1F                       *fHistRCPt[4];                //!Random cone pt
  TH1F                       *fHistRCPtExLJ[4];            //!Random cone pt, imposing min distance from leading jet
  TH1F                       *fHistRCPtRand[4];            //!Random cone pt, randomized particles
  TH2F                       *fHistRCPtExLJVSDPhiLJ;       //!Random cone pt, imposing min distance from leading jet, vs. deltaPhi leading jet
  TH2F                       *fHistRhoVSRCPt;              //!Rho vs. Pt(RCExLJ) / Area(RCExLJ)
  TH1F                       *fHistDeltaPtRC[4];           //!deltaPt = Pt(RC) - A * rho
  TH1F                       *fHistDeltaPtRCExLJ[4];       //!deltaPt = Pt(RC) - A * rho, imposing min distance from leading jet
  TH1F                       *fHistDeltaPtRCRand[4];       //!deltaPt = Pt(RC) - A * rho, randomzied particles

  // Jet embedding
  TH1F                       *fHistEmbJets[4];             //!Pt distribution of embedded jets
  TH1F                       *fHistEmbPart[4];             //!Pt distribution of embedded particle
  TH2F                       *fHistEmbJetPhiEta;           //!Phi-Eta distribution of embedded jets
  TH2F                       *fHistEmbPartPhiEta;          //!Phi-Eta distribution of embedded particles
  TH2F                       *fHistRhoVSEmbBkg;            //!Area(embjet) * rho vs. Pt(embjet) - Pt(embtrack)
  TH1F                       *fHistDeltaPtEmb[4];          //!deltaPt = Pt(embjet) - Area(embjet) * rho - Pt(embtrack)

 private:
  AliAnalysisTaskSAJF(const AliAnalysisTaskSAJF&);            // not implemented
  AliAnalysisTaskSAJF &operator=(const AliAnalysisTaskSAJF&); // not implemented

  ClassDef(AliAnalysisTaskSAJF, 6) // jet analysis task
};
#endif
