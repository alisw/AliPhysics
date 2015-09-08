#ifndef ALIANALYSISTASKEMCALJETTAGGER_H
#define ALIANALYSISTASKEMCALJETTAGGER_H

class TH1;
class TH2;
class TH3;
class AliJetContainer;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskEmcalJetTagger : public AliAnalysisTaskEmcalJet {
 public:
  enum JetTaggingMethod {
    kGeo      = 0,
    kFraction = 1
  };

  enum JetTaggingType {
    kTag      = 0,
    kClosest  = 1
  };

  AliAnalysisTaskEmcalJetTagger();
  AliAnalysisTaskEmcalJetTagger(const char *name);
  virtual ~AliAnalysisTaskEmcalJetTagger();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetJetContainerBase(Int_t c)                             { fContainerBase = c;}
  void SetJetContainerTag(Int_t c)                              { fContainerTag  = c;}

  void SetJetTaggingType(JetTaggingType t)                      { fJetTaggingType = t;}
  void SetJetTaggingMethod(JetTaggingMethod m)                  { fJetTaggingMethod = m;}

  void SetMinFractionShared(Double_t f)                         { fMinFractionShared = f; }

  void SetUseSumw2(Bool_t b)                                    { fUseSumw2 = b;}
  
  void SetTypeAcceptance(Int_t type)                            { fTypeAcc = type; /*see Init()*/}
  void SetMaxDistance(Double_t dist)                            { fMaxDist = dist; }
  
 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();
  void                                Init();
  Double_t GetDeltaPhi(const AliEmcalJet* jet1, const AliEmcalJet* jet2);
  Double_t GetDeltaPhi(Double_t phi1,Double_t phi2);

  void     MatchJetsGeo(Int_t c1 = -1, Int_t c2 = -1, Int_t iDebug = 0, Float_t maxDist = 0.3, Int_t type = 2, Bool_t bReset = kTRUE);
  void     ResetTagging(const Int_t c);
  
 private:
  JetTaggingType                      fJetTaggingType;             // jet matching type
  JetTaggingMethod                    fJetTaggingMethod;           // jet matching method
  Int_t                               fContainerBase;              // jets to be tagged
  Int_t                               fContainerTag;               // jets used for tagging
  Double_t                            fMinFractionShared;          // only fill histos for jets if shared fraction larger than X
  Bool_t                              fUseSumw2;                   // activate sumw2 for output histograms
  Bool_t                              fMatchingDone;               // flag to indicate if matching is done or not
  Int_t                               fTypeAcc;                    // acceptance cut for the jet containers, see method MatchJetsGeo in .cxx for possibilities
  Double_t                            fMaxDist;                    // distance allowed for two jets to match
  Bool_t                              fInit;                       // true when the containers are initialized
  TH3F            **fh3PtJet1VsDeltaEtaDeltaPhi;  //!pt jet 1 vs deta vs dphi
  TH2F            **fh2PtJet1VsDeltaR;            //!pt jet 1 vs dR
  TH2F            **fh2PtJet2VsFraction;          //!pt jet 1 vs shared fraction
  
  TH2F            **fh2PtJet1VsLeadPtAllSel;      //!all jets after std selection
  TH2F            **fh2PtJet1VsLeadPtTagged;      //!tagged jets
  TH2F            **fh2PtJet1VsPtJet2;            //!pT of base jet vs tagged jet
  TH2F            **fh2PtJet2VsRelPt;             //!pT of tagged jet vs pt base jet / pt tagged jet
  
  TH3F             *fh3PtJetDEtaDPhiConst;        //!pt jet vs delta eta vs delta phi of constituents
  TH3F             *fh3PtJetAreaDRConst;          //!pt jet vs Area vs delta R of constituents

  AliAnalysisTaskEmcalJetTagger(const AliAnalysisTaskEmcalJetTagger&);            // not implemented
  AliAnalysisTaskEmcalJetTagger &operator=(const AliAnalysisTaskEmcalJetTagger&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetTagger, 7)
};
#endif

