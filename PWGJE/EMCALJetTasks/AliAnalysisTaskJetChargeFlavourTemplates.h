#ifndef ALIANALYSISTASKJETCHARGETEMPLATES_H
#define ALIANALYSISTASKJETCHARGETEMPLATES_H

class TTree;
class AliAnalysisManager;

#include "AliAnalysisTaskEmcalJet.h"


class AliAnalysisTaskJetChargeFlavourTemplates : public AliAnalysisTaskEmcalJet {
 public:


  AliAnalysisTaskJetChargeFlavourTemplates();
  AliAnalysisTaskJetChargeFlavourTemplates(const char *name);
  virtual ~AliAnalysisTaskJetChargeFlavourTemplates();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetJetContainer(Int_t c)            { fContainer = c; }
  void SetJetPtThreshold(Float_t f)        { fPtThreshold = f; }
  void SetCentralitySelectionOn(Bool_t t)  { fCentSelectOn = t; }
  void SetMinCentrality(Float_t t)         { fCentMin = t; }
  void SetMaxCentrality(Float_t t)         { fCentMax = t; }
  void SetJetRadius(Double_t t)            { fJetRadius = t; }

 protected:
  // Obligatory Functions
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();
  // Essential variables
  // Jet container to be analysed: 0 for raw, 1 for subtracted
  Int_t                               fContainer;
  // Jet lower Pt threshold
  Float_t                             fPtThreshold;
  // Switch on/off centrality selection
  Bool_t                              fCentSelectOn;
  // Min centrality value
  Float_t                             fCentMin;
  // Max centrality value
  Float_t                             fCentMax;
  // Jet radius

  Double_t                            fJetRadius;
  // Splits Jets into too high and Lower Pt Bins
  // User histograms and output tree
  // Histograms first


  TH1F                                *fhJetPt;
  TH1F                                *fhJetPhi;
  TH1F                                *fhJetEta;

  TH1F                                *JC;

  TH1F                                *JCUp;
  TH1F                                *JCDown;
  TH1F                                *JCGluon;
  TH1F                                *JCOther;
  TH1F                                *JCUnmatched;

  TH1F                                *JCLow;

  TH1F                                *JCUpLow;
  TH1F                                *JCDownLow;
  TH1F                                *JCGluonLow;
  TH1F                                *JCOtherLow;
  TH1F                                *JCUnmatchedLow;

  TH1F                                *JCMid;

  TH1F                                *JCUpMid;
  TH1F                                *JCDownMid;
  TH1F                                *JCGluonMid;
  TH1F                                *JCOtherMid;
  TH1F                                *JCUnmatchedMid;

  TH1F                                *JCHigh;

  TH1F                                *JCUpHigh;
  TH1F                                *JCDownHigh;
  TH1F                                *JCGluonHigh;
  TH1F                                *JCOtherHigh;
  TH1F                                *JCUnmatchedHigh;
  // Here is the TTree
  TTree                               *fTreeJets;
  // These are the branch variables; there are nBranches of them
  static const Int_t nBranchesJetChargeFlavourTemplates = 27;
  Double_t                            fTreeBranch[nBranchesJetChargeFlavourTemplates];
  TChain                              *pChain;

 private:
  AliAnalysisTaskJetChargeFlavourTemplates(const AliAnalysisTaskJetChargeFlavourTemplates&);            // not implemented
  AliAnalysisTaskJetChargeFlavourTemplates &operator=(const AliAnalysisTaskJetChargeFlavourTemplates&); // not implemented

  ClassDef(AliAnalysisTaskJetChargeFlavourTemplates, 2)
};
#endif
