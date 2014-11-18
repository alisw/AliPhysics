#ifndef ALIANALYSISTASKEMCALHJETMASS_H
#define ALIANALYSISTASKEMCALHJETMASS_H

class TH1;
class TH2;
class TH3;
class TH3F;
class THnSparse;
class TClonesArray;
class TArrayI;
class AliAnalysisManager;
class AliJetContainer;
class AliEmcalJet;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskEmcalHJetMass : public AliAnalysisTaskEmcalJet {
 public:
  enum JetMassType {
    kRaw   = 0,  //mass form anti-kt 4-vector
    kDeriv = 1   //area based subtracted jet mass
  };

  AliAnalysisTaskEmcalHJetMass();
  AliAnalysisTaskEmcalHJetMass(const char *name);
  virtual ~AliAnalysisTaskEmcalHJetMass();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetJetContainerBase(Int_t c)                             { fContainerBase     = c   ; }
 void SetJetContainerUnsub(Int_t c)                                 { fContainerUnsub    = c   ; }
  void SetMinFractionShared(Double_t f, Bool_t useUnsubJet = kFALSE) { fMinFractionShared = f   ; fUseUnsubJet = useUnsubJet; }
  void SetJetMassType(JetMassType t)                            { fJetMassType       = t   ; }
  void SetMaxDeltaPhi(Double_t dphi)                            { fDPhiHJetMax       = dphi; }

 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHJetHistograms(Double_t pt, const AliEmcalJet *jet);

  Double_t                            GetJetMass(const AliEmcalJet *jet) const;
  Double_t                            GetDeltaPhi(const AliVParticle *vp, const AliEmcalJet* jet) const;
  Double_t                            GetDeltaPhi(Double_t phi1,Double_t phi2) const; 

  Int_t                               fContainerBase;              // jets to be analyzed
  Int_t                               fContainerUnsub;             // unsubtracted jets
  Double_t                            fMinFractionShared;          // only fill histos for jets if shared fraction larger than X
  Bool_t                              fUseUnsubJet;                // calc fraction of unsubtracted jet (relevant for constituent subtraction
  JetMassType                         fJetMassType;                // jet mass type to be used
  Double_t                            fDPhiHJetMax;                // maximum delta phi between hadron and jet

  TH1F            **fh1PtHadron;                        //!pt of hadrons
  TH3F            **fh3PtJet1VsMassVsHPtAllSel;         //!all jets after std selection pt vs mass vs track pt
  TH3F            **fh3PtJet1VsMassVsHPtTagged;         //!tagged jets pt vs mass vs track pt
  TH3F            **fh3PtJet1VsMassVsHPtTaggedMatch;    //!tagged jets pt vs mass vs track pt matched to MC

  TH3F            **fh3PtJet1VsRatVsHPtAllSel;          //!all jets after std selection pt vs mass/pt vs track pt
  TH3F            **fh3PtJet1VsRatVsHPtTagged;          //!tagged jets pt vs mass/pt vs track pt
  TH3F            **fh3PtJet1VsRatVsHPtTaggedMatch;     //!tagged jets pt vs mas/pts vs track pt matched to MC

 private:
  AliAnalysisTaskEmcalHJetMass(const AliAnalysisTaskEmcalHJetMass&);            // not implemented
  AliAnalysisTaskEmcalHJetMass &operator=(const AliAnalysisTaskEmcalHJetMass&); // not implemented

  ClassDef(AliAnalysisTaskEmcalHJetMass, 1)
};
#endif

