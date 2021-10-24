#ifndef AliAnalysisTaskJetCharge_H
#define AliAnalysisTaskJetCharge_H

class TTree;
class AliAnalysisManager;

#include "AliAnalysisTaskEmcalJet.h"


class AliAnalysisTaskJetCharge : public AliAnalysisTaskEmcalJet {
 public:



  AliAnalysisTaskJetCharge();
  AliAnalysisTaskJetCharge(const char *name);
  virtual ~AliAnalysisTaskJetCharge();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetJetContainer(Int_t c)            { fContainer = c; }
  void SetJetPtMinThreshold(Float_t f)        { fPtThreshold = f; }
  void SetCentralitySelectionOn(Bool_t t)  { fCentSelectOn = t; }
  void SetMinCentrality(Float_t t)         { fCentMin = t; }
  void SetMaxCentrality(Float_t t)         { fCentMax = t; }
  void SetJetRadius(Double_t t)            { fJetRadius = t; }
  void SetK(Double_t t)                    { JetChargeK = t; }
  void SetJetMidPt(Double_t t)             { JetMidPt = t; }
  void SetJetHighPt(Double_t t)            { JetHighPt = t; }

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
  // Value of K the scaling factor on jet charge default k = 0.5
  Double_t                            JetChargeK = 0.5;                     //
  // Set the Mid Pt jet threshold Defult = 40
  Double_t                            JetMidPt = 40;                        //
  // Set the Highest jet Pt Band threshold Default = 40
  Double_t                            JetHighPt = 80;                       //
  // Splits Jets into too high and Lower Pt Bins
  // User histograms and output tree
  // Histograms first

  TH1F                                *fhJetPt;                             //!
  TH1F                                *fhJetPhi;                            //!
  TH1F                                *fhJetEta;                            //!

  TH1F                                *fhJetCharge;                         //!
  TH1F                                *fhJetChargeLow;                      //!
  TH1F                                *fhJetChargeMid;                      //!
  TH1F                                *fhJetChargeHigh;                     //!

  // Here is the TTree
  TTree                               *fTreeJets;
  // These are the branch variables; there are nBranchesJetCharge of them
  static const Int_t nBranchesJetCharge = 7;
  Double_t                            fTreeBranch[nBranchesJetCharge];
  TChain                              *pChain;

 private:
  AliAnalysisTaskJetCharge(const AliAnalysisTaskJetCharge&);            // not implemented
  AliAnalysisTaskJetCharge &operator=(const AliAnalysisTaskJetCharge&); // not implemented

  ClassDef(AliAnalysisTaskJetCharge, 3)
};
#endif
