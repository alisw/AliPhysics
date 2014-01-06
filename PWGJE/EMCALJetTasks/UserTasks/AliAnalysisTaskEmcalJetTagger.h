#ifndef ALIANALYSISTASKEMCALJETTAGGER_H
#define ALIANALYSISTASKEMCALJETTAGGER_H

class TH1;
class TH2;
class TH3;
class TH3F;
class THnSparse;
class TClonesArray;
class TArrayI;
class AliAnalysisManager;
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

 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();


  Double_t GetDeltaPhi(const AliEmcalJet* jet1, const AliEmcalJet* jet2);
  Double_t GetDeltaPhi(Double_t phi1,Double_t phi2);
  Double_t GetDeltaR(const AliEmcalJet* jet1, const AliEmcalJet* jet2) const;
  Double_t GetFractionSharedPt(const AliEmcalJet *jet1, const AliEmcalJet *jet2) const;

  void     MatchJetsGeo(Int_t c1 = -1, Int_t c2 = -1, Int_t iDebug = 0, Float_t maxDist = 0.3, Int_t type = 2);
  void     ResetTagging(const Int_t c);
  
  JetTaggingType                      fJetTaggingType;             // jet matching type
  JetTaggingMethod                    fJetTaggingMethod;           // jet matching method

  Int_t                               fContainerBase;              // jets to be tagged
  Int_t                               fContainerTag;               // jets used for tagging
  
 private:
  Bool_t            fMatchingDone;                // flag to indicate if matching is done or not
  TH3F            **fh3PtJet1VsDeltaEtaDeltaPhi;  //!pt jet 1 vs deta vs dphi
  TH2F            **fh2PtJet1VsDeltaR;            //!pt jet 1 vs dR
  
  TH2F            **fh2PtJet1VsLeadPtAllSel;      //!all jets after std selection
  TH2F            **fh2PtJet1VsLeadPtTagged;      //!tagged jets
  TH2F            **fh2PtJet1VsPtJet2;            //!pT of base jet vs tagged jet

  TH3F             *fh3PtJetDEtaDPhiConst;        //!pt jet vs delta eta vs delta phi of constituents
  TH2F             *fh2PtJetDRConst;              //!pt jet vs delta R of constituents
  TH3F             *fh3PtJetAreaDRConst;          //!pt jet vs Area vs delta R of constituents

  AliAnalysisTaskEmcalJetTagger(const AliAnalysisTaskEmcalJetTagger&);            // not implemented
  AliAnalysisTaskEmcalJetTagger &operator=(const AliAnalysisTaskEmcalJetTagger&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetTagger, 2)
};
#endif

