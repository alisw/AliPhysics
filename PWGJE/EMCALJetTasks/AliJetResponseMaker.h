#ifndef ALIJETRESPONSEMAKER_H
#define ALIJETRESPONSEMAKER_H

// $Id$

class AliGenPythiaEventHeader;
class TClonesArray;
class TH1;
class TH2;
class TProfile;
class AliNamedArrayI;

#include "AliEmcalJet.h"
#include "AliAnalysisTaskEmcalJet.h"

class AliJetResponseMaker : public AliAnalysisTaskEmcalJet {
 public:
  AliJetResponseMaker();
  AliJetResponseMaker(const char *name);
  virtual ~AliJetResponseMaker();

  enum MatchingType{
    kNoMatching = 0,
    kGeometrical = 1,
    kMCLabel = 2,
    kSameCollections = 3
  };

  void                        UserCreateOutputObjects();
  Bool_t                      UserNotify();

  void                        SetJets2Name(const char *n)                                     { fJets2Name         = n         ; }
  void                        SetTracks2Name(const char *n)                                   { fTracks2Name       = n         ; }
  void                        SetClus2Name(const char *n)                                     { fCalo2Name         = n         ; }
  void                        SetJet2EtaLimits(Float_t min=-999, Float_t max=-999)            { fJet2MinEta = min, fJet2MaxEta = max ; }
  void                        SetJet2PhiLimits(Float_t min=-999, Float_t max=-999)            { fJet2MinPhi = min, fJet2MaxPhi = max ; }
  void                        SetRho2Name(const char *n)                                      { fRho2Name          = n         ; }
  void                        SetPtBiasJet2Clus(Float_t b)                                    { fPtBiasJet2Clus    = b         ; }
  void                        SetPtBiasJet2Track(Float_t b)                                   { fPtBiasJet2Track   = b         ; }
  void                        SetMatching(MatchingType t, Double_t p1=1, Double_t p2=1)       { fMatching = t; fMatchingPar1 = p1; fMatchingPar2 = p2; }
  void                        SetPtHardBin(Int_t b)                                           { fSelectPtHardBin   = b         ; }
  void                        SetAreMCCollections(Bool_t f1, Bool_t f2)                       { fAreCollections1MC = f1; fAreCollections2MC = f2; }
  void                        SetIsEmbedded(Bool_t i)                                         { fIsEmbedded        = i         ; }
  void                        SetIsPythia(Bool_t i)                                           { fIsPythia          = i         ; }

 protected:
  Bool_t                      PythiaInfoFromFile(const char* currFile, Float_t &fXsec, Float_t &fTrials, Int_t &pthard);
  Bool_t                      IsEventSelected();
  Bool_t                      AcceptJet(AliEmcalJet* jet) const;
  Bool_t                      AcceptBiasJet2(AliEmcalJet *jet) const;
  void                        ExecOnce();
  void                        DoJetLoop(Bool_t order);
  Bool_t                      FillHistograms();
  Bool_t                      RetrieveEventObjects();
  Bool_t                      Run();
  Bool_t                      DoJetMatching();
  void                        SetMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, MatchingType matching);
  void                        GetGeometricalMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d) const;
  void                        GetMCLabelMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d1, Double_t &d2) const;
  void                        GetSameCollectionsMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d1, Double_t &d2) const;

  TString                     fTracks2Name;                   // name of second track collection
  TString                     fCalo2Name;                     // name of second cluster collection
  TString                     fJets2Name;                     // name of second jet collection
  TString                     fRho2Name;                      // name of second jet collection
  Float_t                     fPtBiasJet2Track;               // select jets 2 with a minimum pt track
  Float_t                     fPtBiasJet2Clus;                // select jets 2 with a minimum pt cluster
  Bool_t                      fAreCollections1MC;             // collections 1 MC
  Bool_t                      fAreCollections2MC;             // collections 1 MC
  MatchingType                fMatching;                      // matching type
  Double_t                    fMatchingPar1;                  // matching parameter for jet1-jet2 matching
  Double_t                    fMatchingPar2;                  // matching parameter for jet2-jet1 matching
  Float_t                     fJet2MinEta;                    // minimum eta jet 2 acceptance
  Float_t                     fJet2MaxEta;                    // maximum eta jet 2 acceptance
  Float_t                     fJet2MinPhi;                    // minimum phi jet 2 acceptance
  Float_t                     fJet2MaxPhi;                    // maximum phi jet 2 acceptance  
  Int_t                       fSelectPtHardBin;               // select one pt hard bin for analysis
  Bool_t                      fIsEmbedded;                    // trigger, embedded signal
  Bool_t                      fIsPythia;                      // trigger, if it is a PYTHIA production

  AliGenPythiaEventHeader    *fPythiaHeader;                  //!event Pythia header
  Int_t                       fPtHardBin;                     //!event pt hard bin
  Int_t                       fNTrials;                       //!event trials
  TClonesArray               *fTracks2;                       //!Tracks 2
  TClonesArray               *fCaloClusters2;                 //!Clusters 2
  TClonesArray               *fJets2;                         //!Jets 2
  AliRhoParameter            *fRho2;                          //!Event rho 2
  Double_t                    fRho2Val;                       //!Event rho 2 value 
  AliNamedArrayI             *fTracks2Map;                    //!MC particle map
  // General histograms
  TH1                        *fHistTrialsAfterSel;            //!total number of trials per pt hard bin after selection
  TH1                        *fHistEventsAfterSel;            //!total number of events per pt hard bin after selection
  TH1                        *fHistTrials;                    //!trials from pyxsec.root
  TProfile                   *fHistXsection;                  //!x section from pyxsec.root
  TH1                        *fHistEvents;                    //!total number of events per pt hard bin
  // Jets 1
  TH2                        *fHistJets1PhiEta;               //!phi-eta distribution of jets 1
  TH2                        *fHistJets1PtArea;               //!inclusive jet pt vs area histogram 1
  TH2                        *fHistJets1CorrPtArea;           //!inclusive jet pt vs. area histogram 1
  TH2                        *fHistLeadingJets1PtArea;        //!leading jet pt vs area histogram 1
  TH2                        *fHistLeadingJets1CorrPtArea;    //!leading jet pt vs. area histogram 1
  TH2                        *fHistJets1NEFvsPt;              //!Jet neutral energy fraction vs. jet pt 1
  TH2                        *fHistJets1CEFvsCEFPt;           //!Jet charged energy fraction vs. charged jet pt 1
  TH2                        *fHistJets1ZvsPt;                //!Constituent Pt over Jet Pt ratio vs. jet pt 1
  // Jets 2
  TH2                        *fHistJets2PhiEta;               //!phi-eta distribution of jets 2
  TH2                        *fHistJets2PtArea;               //!inclusive jet pt vs. area histogram 2
  TH2                        *fHistJets2CorrPtArea;           //!inclusive jet pt vs. area histogram 2
  TH2                        *fHistLeadingJets2PtArea;        //!leading jet pt vs. area histogram 2
  TH2                        *fHistLeadingJets2CorrPtArea;    //!leading jet pt vs. area histogram 2
  TH2                        *fHistJets2PhiEtaAcceptance;     //!phi-eta distribution of jets 2 using jet 1 cuts (acceptance, leading hadron bias, ...)
  TH2                        *fHistJets2PtAreaAcceptance;     //!inclusive jet pt vs. area histogram 2 using jet 1 cuts (acceptance, leading hadron bias, ...)
  TH2                        *fHistJets2CorrPtAreaAcceptance; //!inclusive jet pt vs. area histogram 2 using jet 1 cuts (acceptance, leading hadron bias, ...)
  TH2                        *fHistLeadingJets2PtAreaAcceptance;     //!leading jet pt vs. area histogram 2 using jet 1 cuts (acceptance, leading hadron bias, ...)
  TH2                        *fHistLeadingJets2CorrPtAreaAcceptance; //!leading jet pt vs. area histogram 2 using jet 1 cuts (acceptance, leading hadron bias, ...)
  TH2                        *fHistJets2NEFvsPt;              //!Jet neutral energy fraction vs. jet pt 2
  TH2                        *fHistJets2CEFvsCEFPt;           //!Jet charged energy fraction vs. charged jet pt 2
  TH2                        *fHistJets2ZvsPt;                //!Constituent Pt over Jet Pt ratio vs. jet pt 2
  // Jet1-Jet2 matching
  TH2                        *fHistCommonEnergy1vsJet1Pt;              //!common energy 1 (%) vs jet 1 pt
  TH2                        *fHistCommonEnergy2vsJet2Pt;              //!common energy 2 (%) vs jet 2 pt
  TH2                        *fHistDistancevsJet1Pt;                   //!distance vs jet 1 pt
  TH2                        *fHistDistancevsJet2Pt;                   //!distance vs jet 2 pt
  TH2                        *fHistDistancevsCommonEnergy1;            //!distance vs common energy 1 (%)
  TH2                        *fHistDistancevsCommonEnergy2;            //!distance vs common energy 2 (%)
  TH2                        *fHistJet2PtOverJet1PtvsJet2Pt;           //!jet 2 pt over jet 1 pt vs jet 2 pt
  TH2                        *fHistJet1PtOverJet2PtvsJet1Pt;           //!jet 1 pt over jet 2 pt vs jet 1 pt
  TH2                        *fHistDeltaEtaPhi;                        //!delta eta-phi between matched jets
  TH2                        *fHistDeltaPtvsJet1Pt;                    //!delta pt between matched jets vs jet 1 pt
  TH2                        *fHistDeltaPtvsJet2Pt;                    //!delta pt between matched jets vs jet 2 pt
  TH2                        *fHistDeltaPtvsMatchingLevel;             //!delta pt between matched jets vs matching level
  TH2                        *fHistDeltaCorrPtvsJet1Pt;                //!delta pt corr between matched jets vs jet 1 pt
  TH2                        *fHistDeltaCorrPtvsJet2Pt;                //!delta pt corr between matched jets vs jet 2 pt
  TH2                        *fHistDeltaCorrPtvsMatchingLevel;         //!delta pt corr between matched jets vs matching level
  TH2                        *fHistNonMatchedJets1PtArea;              //!non-matched jet 1 pt distribution
  TH2                        *fHistNonMatchedJets2PtArea;              //!non-matched jet 2 pt distribution
  TH2                        *fHistNonMatchedJets1CorrPtArea;          //!non-matched jet pt distribution
  TH2                        *fHistNonMatchedJets2CorrPtArea;          //!non-matched jet pt distribution
  TH2                        *fHistJet1PtvsJet2Pt;                     //!correlation jet 1 pt vs jet 2 pt
  TH2                        *fHistJet1CorrPtvsJet2CorrPt;             //!correlation jet 1 corr pt vs jet 2 corr pt
  TH2                        *fHistJet1MCPtvsJet2Pt;                   //!correlation jet 1 MC pt vs jet 2 pt
  TH2                        *fHistMissedJets2PtArea;                  //!jets 2 not found in jet 1 collection


 private:
  AliJetResponseMaker(const AliJetResponseMaker&);            // not implemented
  AliJetResponseMaker &operator=(const AliJetResponseMaker&); // not implemented

  ClassDef(AliJetResponseMaker, 13) // Jet response matrix producing task
};
#endif
