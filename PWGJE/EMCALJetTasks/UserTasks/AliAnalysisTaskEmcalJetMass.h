#ifndef ALIANALYSISTASKEMCALJETMASS_H
#define ALIANALYSISTASKEMCALJETMASS_H

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

class AliAnalysisTaskEmcalJetMass : public AliAnalysisTaskEmcalJet {
 public:
  enum JetMassType {
    kRaw   = 0,  //mass form anti-kt 4-vector
    kDeriv = 1   //area based subtracted jet mass
  };

  AliAnalysisTaskEmcalJetMass();
  AliAnalysisTaskEmcalJetMass(const char *name);
  virtual ~AliAnalysisTaskEmcalJetMass();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetJetContainerBase(Int_t c)                                  { fContainerBase     = c   ; }
  void SetJetContainerUnsub(Int_t c)                                 { fContainerUnsub    = c   ; }
  void SetMinFractionShared(Double_t f, Bool_t useUnsubJet = kFALSE) { fMinFractionShared = f   ; SetUseUnsubJet(useUnsubJet); }
  void SetUseUnsubJet(Bool_t b)                                      { fUseUnsubJet       = b   ; }
  void SetJetMassType(JetMassType t)                                 { fJetMassType       = t   ; }
  void SetUseSumw2(Bool_t b)                                         { fUseSumw2          = b   ; }

 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();

  Double_t                            GetJetMass(AliEmcalJet *jet);
 
 private:
  Int_t                               fContainerBase;              // jets to be analyzed
  Int_t                               fContainerUnsub;             // unsubtracted jets
  Double_t                            fMinFractionShared;          // only fill histos for jets if shared fraction larger than X
  Bool_t                              fUseUnsubJet;                // calc fraction of unsubtracted jet (relevant for constituent subtraction)
  JetMassType                         fJetMassType;                // jet mass type to be used
  Bool_t                              fUseSumw2;                   // activate sumw2 for output histograms

  TH3F            **fh3PtJet1VsMassVsLeadPtAllSel;         //!all jets after std selection pt vs mass vs leading track pt
  TH3F            **fh3PtJet1VsMassVsLeadPtTagged;         //!tagged jets pt vs mass vs leading track pt
  TH3F            **fh3PtJet1VsMassVsLeadPtTaggedMatch;    //!tagged jets pt vs mass vs leading track pt matched to MC
  TProfile        **fpPtVsMassJet1All;                     //!pT vs avg mass of all jets
  TProfile        **fpPtVsMassJet1Tagged;                  //!pT vs avg mass of tagged jets
  TProfile        **fpPtVsMassJet1TaggedMatch;             //!pT vs avg mass of tagged jets matched to MC
  TH2F            **fh2MassVsAreaJet1All;                  //!mass vs area of all jets
  TH2F            **fh2MassVsAreaJet1Tagged;               //!mass vs area of tagged jets
  TH2F            **fh2MassVsAreaJet1TaggedMatch;          //!mass vs area of tagged jets matched to MC
  TH2F            **fh2MassVsNConstJet1All;                //!mass vs number of constituents of all jets
  TH2F            **fh2MassVsNConstJet1Tagged;             //!mass vs number of constituents of tagged jets
  TH2F            **fh2MassVsNConstJet1TaggedMatch;        //!mass vs number of constituents of tagged jets matched to MC

  TH3F            *fh3PtJet1VsMassVsCentAllSel;            //!all jets after std selection: pt vs mass vs centrality
  TH3F            *fh3PtJet1VsMassVsCentTagged;            //!tagged jets: pt vs mass vs centrality
  TH3F            *fh3PtJet1VsMassVsCentTaggedMatch;       //!tagged jets: pt vs mass vs centrality matched to MC

  TH3F            **fh3PtJet1VsRatVsLeadPtAllSel;          //!all jets after std selection pt vs mass/pt vs leading track pt
  TH3F            **fh3PtJet1VsRatVsLeadPtTagged;          //!tagged jets pt vs mass/pt vs leading track pt
  TH3F            **fh3PtJet1VsRatVsLeadPtTaggedMatch;     //!tagged jets pt vs mas/pts vs leading track pt matched to MC
  TProfile        **fpPtVsRatJet1All;                      //!pT vs avg mass/pt of all jets
  TProfile        **fpPtVsRatJet1Tagged;                   //!pT vs avg mass/pt of tagged jets
  TProfile        **fpPtVsRatJet1TaggedMatch;              //!pT vs avg mass/pt of tagged jets matched to MC
  TH2F            **fh2RatVsAreaJet1All;                   //!mass/pt vs area of all jets
  TH2F            **fh2RatVsAreaJet1Tagged;                //!mass/pt vs area of tagged jets
  TH2F            **fh2RatVsAreaJet1TaggedMatch;           //!mass/pt vs area of tagged jets matched to MC
  TH2F            **fh2RatVsNConstJet1All;                 //!mass/pt vs number of constituents of all jets
  TH2F            **fh2RatVsNConstJet1Tagged;              //!mass/pt vs number of constituents of tagged jets
  TH2F            **fh2RatVsNConstJet1TaggedMatch;         //!mass/pt vs number of constituents of tagged jets matched to MC

  TH3F            **fh3JetPtVsMassVsEPRelAllSel;           //!jet pt vs mass vs (phi-ep) for all jets
  TH3F            **fh3JetPtVsMassVsEPRelTagged;           //!jet pt vs mass vs (phi-ep) for tagged jets
  TH3F            **fh3JetPtVsMassVsEPRelTaggedMatch;      //!jet pt vs mass vs (phi-ep) for tagged matched jets
  THnSparse       *fhnJetPtMassAndBkg;                     //! correlation between leading jet kinematics and event background
  
  AliAnalysisTaskEmcalJetMass(const AliAnalysisTaskEmcalJetMass&);            // not implemented
  AliAnalysisTaskEmcalJetMass &operator=(const AliAnalysisTaskEmcalJetMass&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetMass, 11)
};
#endif

