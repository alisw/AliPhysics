#ifndef ALIANALYSISTASKJETCharge_H
#define ALIANALYSISTASKJETCharge_H

class TH1;
class TH2;
class TH3;
class TH3F;
class TTree;
class THnSparse;
class TClonesArray;
class TArrayI;
class AliAnalysisManager;
class AliJetContainer;
class AliEmcalJetFinder;
class AliFJWrapper;

#include "AliAnalysisTaskEmcalJet.h"
#include "AliFJWrapper.h"
#include "AliClusterContainer.h"
#include <TChain.h>
//#include <fastjet/PseudoJet.hh>
//#include <fastjet/SharedPtr.hh>



class AliAnalysisTaskJetCharge : public AliAnalysisTaskEmcalJet {
 public:
   AliEventCuts fEventCuts; // Event Cuts




  enum JetShapeSub {
    kNoSub = 0,
    kConstSub = 1
  };

  AliAnalysisTaskJetCharge();
  AliAnalysisTaskJetCharge(const char *name);
  virtual ~AliAnalysisTaskJetCharge();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetJetContainer(Int_t c)            { fContainer = c; }
  void SetJetPtThreshold(Float_t f)        { fPtThreshold = f; }
  void SetCentralitySelectionOn(Bool_t t)  { fCentSelectOn = t; }
  void SetMinCentrality(Float_t t)         { fCentMin = t; }
  void SetMaxCentrality(Float_t t)         { fCentMax = t; }
  void SetJetShapeSub(JetShapeSub t)       { fJetShapeSub = t; }
  void SetJetRadius(Double_t t)            { fJetRadius = t; }

 protected:
  // Obligatory Functions
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();
  // Essential variables
  // Jet container to be analysed: 0 for raw, 1 for subtracted
  Int_t                               fContainer;
  // Jet subtraction method
  JetShapeSub                         fJetShapeSub;
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

  TH1F                                *fhJetCharge;
  TH1F                                *fhJetChargeLow;
  TH1F                                *fhJetChargeMid;
  TH1F                                *fhJetChargeHigh;

  // Here is the TTree
  TTree                               *fTreeJets;
  const Int_t nBranches = 9;
  // These are the branch variables; there are nBranches of them
  Double_t                            fTreeBranch[nBranches];
  TChain                              *pChain;

 private:
  AliAnalysisTaskJetCharge(const AliAnalysisTaskJetCharge&);            // not implemented
  AliAnalysisTaskJetCharge &operator=(const AliAnalysisTaskJetCharge&); // not implemented

  ClassDef(AliAnalysisTaskJetCharge, 1)
};
#endif
