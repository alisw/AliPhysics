#ifndef ALIANALYSISTASKJETCHARGE_H
#define ALIANALYSISTASKJETCHARGE_H


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
  void SetJetPtMinThreshold(Float_t f)     { fPtThreshold = f; }
  void SetCentralitySelectionOn(Bool_t t)  { fCentSelectOn = t; }
  void SetMinCentrality(Float_t t)         { fCentMin = t; }
  void SetMaxCentrality(Float_t t)         { fCentMax = t; }
  void SetJetRadius(Float_t t)            { fJetRadius = t; }
  void SetK(Float_t t)                    { JetChargeK = t; }

 protected:
  // Obligatory Functions
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();
  // Essential variables
  // Jet container to be analysed: 0 for raw, 1 for subtracted
  Int_t                               fContainer;                           //
  // Jet lower Pt threshold
  Float_t                             fPtThreshold;                         //
  // Switch on/off centrality selection
  Bool_t                              fCentSelectOn;                        //
  // Min centrality value
  Float_t                             fCentMin;                             //
  // Max centrality value
  Float_t                             fCentMax;                             //
  // Jet radius
  Float_t                             fJetRadius;
  // Value of K the scaling factor on jet charge default k = 0.5
  Float_t                             JetChargeK = 0.5;                     //

  // Here is the TTree
  TTree                               *fTreeJets;                           //!
  // These are the branch variables; there are nBranches of them

  Float_t Pt;                     //
  Float_t Phi;                    //
  Float_t Eta;                    //
  Float_t JetCharge;              //
  Float_t JetChargeMid;              //
  Float_t JetChargeHigh;              //
  Float_t LeadingTrackPt;         //


  TChain                              *pChain;                                                     //



  ClassDef(AliAnalysisTaskJetCharge, 6)
};
#endif
