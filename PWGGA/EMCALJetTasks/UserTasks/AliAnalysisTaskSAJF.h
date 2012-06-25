#ifndef ALIANALYSISTASKSAJF_H
#define ALIANALYSISTASKSAJF_H

// $Id$

class TClonesArray;
class TString;
class TH1F;
class TH2F;
class AliRhoParameter;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskSAJF : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskSAJF();
  AliAnalysisTaskSAJF(const char *name);
  virtual ~AliAnalysisTaskSAJF();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  void                        SetJetMinRC2LJ(Float_t d)                            { fMinRC2LJ                = d          ; } 
  void                        SetEmbJetsName(const char *n)                        { fEmbJetsName             = n          ; }
  void                        SetRandTracksName(const char *n)                     { fRandTracksName          = n          ; }
  void                        SetRandClusName(const char *n)                       { fRandCaloName            = n          ; }
  void                        SetEmbTracksName(const char *n)                      { fEmbTracksName           = n          ; }
  void                        SetEmbClusName(const char *n)                        { fEmbCaloName             = n          ; }
  void                        SetRhoName(const char *n)                            { fRhoName                 = n          ; }        

 protected:
  void                        ExecOnce()                                                                                    ;
  Bool_t                      RetrieveEventObjects()                                                                        ;
  Bool_t                      FillHistograms()                                                                              ;
  void                        GetLeadingJets(Int_t &maxJetIndex, Int_t &max2JetIndex)                                       ;
  void                        DoJetLoop()                                                                                   ;
  void                        DoEmbJetLoop(AliEmcalJet* &embJet, TObject* &embPart)                                         ;
  void                        DoTrackLoop()                                                                                 ;
  void                        DoClusterLoop()                                                                               ;
  void                        GetRigidCone(Float_t &pt, Float_t &ptrigid, Float_t &eta, Float_t &phi, 
					   AliEmcalJet *jet = 0, TClonesArray* tracks = 0, TClonesArray* clusters = 0) const;

  Float_t                     fMinRC2LJ;                   // Minimum distance random cone to leading jet
  TString                     fEmbJetsName;                // Name of embedded jet collection
  TString                     fEmbTracksName;              // Name of embedded track collection
  TString                     fEmbCaloName;                // Name of embedded calo cluster collection
  TString                     fRandTracksName;             // Name of randomized track collection
  TString                     fRandCaloName;               // Name of randomized calo cluster collection
  TString                     fRhoName;                    // Name of rho object

  TClonesArray               *fEmbJets;                    //!Embedded Jets
  TClonesArray               *fEmbTracks;                  //!Embedded tracks
  TClonesArray               *fEmbCaloClusters;            //!Embedded clusters  
  TClonesArray               *fRandTracks;                 //!Randomized tracks
  TClonesArray               *fRandCaloClusters;           //!Randomized clusters
  AliRhoParameter            *fRho;                        //!Event rho
  Double_t                    fRhoVal;                     //!Event rho value
  Int_t                       fEmbeddedClusterId;          //!Embedded cluster id
  Int_t                       fEmbeddedTrackId;            //!Embedded track id

  // General histograms
  TH1F                       *fHistCentrality;             //!Event centrality distribution
  TH1F                       *fHistEvents[4];              //!Events accepted/rejected
  TH1F                       *fHistTracksPt[4];            //!Inclusive track pt spectrum
  TH1F                       *fHistClustersPt[4];          //!Inclusive clusters pt spectrum
  TH2F                       *fHistJetPhiEta[4];           //!Phi-Eta distribution of jets
  TH1F                       *fHistJetsPt[4];              //!Inclusive jet pt spectrum
  TH2F                       *fHistJetsPtArea[4];          //!Jet pt vs. area
  TH1F                       *fHistLeadingJetPt[4];        //!Leading jet pt spectrum
  TH1F                       *fHist2LeadingJetPt[4];       //!Second leading jet pt spectrum
  TH2F                       *fHistJetsNEFvsPt[4];         //!Jet neutral energy fraction vs. jet pt
  TH2F                       *fHistJetsZvsPt[4];           //!Constituent Pt over Jet Pt ratio vs. jet pt
  TH2F                       *fHistMaxTrackPtvsJetPt[4];   //!Max constituent track pt vs. jet pt
  TH2F                       *fHistMaxClusPtvsJetPt[4];    //!Max constituent cluster pt vs. jet pt
  TH2F                       *fHistMaxPartPtvsJetPt[4];    //!Max constituent particle (track or cluster) pt vs. jet pt
  TH2F                       *fHistMaxTrackPtvsJetCorrPt[4];   //!Max constituent track pt vs. jet pt
  TH2F                       *fHistMaxClusPtvsJetCorrPt[4];    //!Max constituent cluster pt vs. jet pt
  TH2F                       *fHistMaxPartPtvsJetCorrPt[4];    //!Max constituent particle (track or cluster) pt vs. jet pt
  TH1F                       *fHistDeltaVectorPt;              //!Delta Pt between vector and scalar sum

  // Rho
  TH1F                       *fHistRho[4];                    //!Rho distribution
  TH2F                       *fHistRhoVSleadJetPt;            //!Area(leadjetarea) * rho vs. leading jet pt
  TH1F                       *fHistCorrJetsPt[4];             //!Corrected inclusive jet pt spectrum
  TH2F                       *fHistCorrJetsPtArea[4];         //!Jet pt vs. area
  TH1F                       *fHistCorrLeadingJetPt[4];       //!Corrected leading jet pt spectrum

  // Random cones
  TH2F                       *fHistRCPhiEta;               //!Phi-Eta distribution of random cones
  TH1F                       *fHistRCPtRigid[4];           //!Random cone pt, rigid
  TH1F                       *fHistRCPt[4];                //!Random cone pt
  TH1F                       *fHistRCPtExLJ[4];            //!Random cone pt, imposing min distance from leading jet
  TH1F                       *fHistRCPtRand[4];            //!Random cone pt, randomized particles
  TH2F                       *fHistRCPtExLJVSDPhiLJ;       //!Random cone pt, imposing min distance from leading jet, vs. deltaPhi leading jet
  TH2F                       *fHistRhoVSRCPt;              //!Rho vs. Pt(RCExLJ) / Area(RCExLJ)
  TH1F                       *fHistDeltaPtRCRigid[4];      //!deltaPt = Pt(RC) - A * rho, rigid
  TH1F                       *fHistDeltaPtRC[4];           //!deltaPt = Pt(RC) - A * rho
  TH1F                       *fHistDeltaPtRCExLJ[4];       //!deltaPt = Pt(RC) - A * rho, imposing min distance from leading jet
  TH1F                       *fHistDeltaPtRCRand[4];       //!deltaPt = Pt(RC) - A * rho, randomzied particles

  // Jet embedding
  TH2F                       *fHistEmbNotFoundPhiEta;      //!Phi-Eta of "not found" embedded particles
  TH1F                       *fHistEmbJetsPt[4];           //!Pt distribution of embedded jets
  TH1F                       *fHistEmbJetsCorrPt[4];       //!Pt distribution of embedded jets
  TH1F                       *fHistEmbPart[4];             //!Pt distribution of embedded particle
  TH2F                       *fHistEmbJetPhiEta;           //!Phi-Eta distribution of embedded jets
  TH2F                       *fHistEmbPartPhiEta;          //!Phi-Eta distribution of embedded particles
  TH2F                       *fHistRhoVSEmbBkg;            //!Area(embjet) * rho vs. Pt(embjet) - Pt(embtrack)
  TH1F                       *fHistDeltaPtEmb[4];          //!deltaPt = Pt(embjet) - Area(embjet) * rho - Pt(embtrack)

 private:
  AliAnalysisTaskSAJF(const AliAnalysisTaskSAJF&);            // not implemented
  AliAnalysisTaskSAJF &operator=(const AliAnalysisTaskSAJF&); // not implemented

  ClassDef(AliAnalysisTaskSAJF, 9) // jet analysis task
};
#endif
