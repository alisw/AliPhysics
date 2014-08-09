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
  void SetJetContainerBase(Int_t c)                             { fContainerBase     = c   ; }
  void SetMinFractionShared(Double_t f)                         { fMinFractionShared = f   ; }
  void SetJetMassType(JetMassType t)                            { fJetMassType       = t   ; }

 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();

  Double_t                            GetJetMass(AliEmcalJet *jet);
 
  Int_t                               fContainerBase;              // jets to be analyzed
  Double_t                            fMinFractionShared;          // only fill histos for jets if shared fraction larger than X
  JetMassType                         fJetMassType;                // jet mass type to be used

  TH3F            **fh3PtJet1VsMassVsLeadPtAllSel;         //!all jets after std selection pt vs mass vs leading track pt
  TH3F            **fh3PtJet1VsMassVsLeadPtTagged;         //!tagged jets pt vs mass vs leading track pt
  TH3F            **fh3PtJet1VsMassVsLeadPtTaggedMatch;    //!tagged jets pt vs mass vs leading track pt matched to MC
  TProfile        **fpPtVsMassJet1All;               //!pT vs avg mass of all jets
  TProfile        **fpPtVsMassJet1Tagged;            //!pT vs avg mass of tagged jets
  TProfile        **fpPtVsMassJet1TaggedMatch;       //!pT vs avg mass of tagged jets matched to MC
  TH2F            **fh2MassVsAreaJet1All;            //!mass vs area of all jets
  TH2F            **fh2MassVsAreaJet1Tagged;         //!mass vs area of tagged jets
  TH2F            **fh2MassVsAreaJet1TaggedMatch;    //!mass vs area of tagged jets matched to MC
  TH2F            **fh2MassVsNConstJet1All;          //!mass vs number of constituents of all jets
  TH2F            **fh2MassVsNConstJet1Tagged;       //!mass vs number of constituents of tagged jets
  TH2F            **fh2MassVsNConstJet1TaggedMatch;  //!mass vs number of constituents of tagged jets matched to MC
  TH2F            **fh2EtMassOverEtRSq;           //!Et vs (M/Et*R)^2

 private:
  AliAnalysisTaskEmcalJetMass(const AliAnalysisTaskEmcalJetMass&);            // not implemented
  AliAnalysisTaskEmcalJetMass &operator=(const AliAnalysisTaskEmcalJetMass&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetMass, 6)
};
#endif

