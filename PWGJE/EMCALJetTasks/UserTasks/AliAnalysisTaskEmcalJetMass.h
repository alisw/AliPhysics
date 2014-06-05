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

  AliAnalysisTaskEmcalJetMass();
  AliAnalysisTaskEmcalJetMass(const char *name);
  virtual ~AliAnalysisTaskEmcalJetMass();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetJetContainerBase(Int_t c)                             { fContainerBase     = c   ; }
  void SetMinFractionShared(Double_t f)                         { fMinFractionShared = f   ; }
  void SetJetMassAverage(Double_t m)                            { fJetMassAvg        = m   ; }

 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();

  Double_t                            GetJetMass(AliEmcalJet *jet);
  TLorentzVector                      GetSubtractedVector(AliEmcalJet *jet);
  TLorentzVector                      GetBkgVector(AliEmcalJet *jet, AliJetContainer *cont);
 
  Int_t                               fContainerBase;              // jets to be analyzed
  Double_t                            fMinFractionShared;          // only fill histos for jets if shared fraction larger than X
  Double_t                            fJetMassAvg;                 // average jet mass

  TH2F            **fh2PtJet1VsLeadPtAllSel;         //!all jets after std selection vs leading track pt
  TH2F            **fh2PtJet1VsLeadPtTagged;         //!tagged jets vs leading track pt
  TH2F            **fh2PtJet1VsLeadPtTaggedMatch;    //!tagged jets vs leading track pt matched to MC
  TH2F            **fh2PtVsMassJet1All;              //!pT vs mass of all jets
  TH2F            **fh2PtVsMassJet1Tagged;           //!pT vs mass of tagged jets
  TH2F            **fh2PtVsMassJet1TaggedMatch;      //!pT vs mass of tagged jets matched to MC
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

  ClassDef(AliAnalysisTaskEmcalJetMass, 4)
};
#endif

