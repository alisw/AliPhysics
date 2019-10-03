#ifndef ALIANALYSISTASKDELTAPTJEMB_H
#define ALIANALYSISTASKDELTAPTJEMB_H

// $Id$

class TClonesArray;
class TString;
class TH1;
class TH2;
class TH3;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskDeltaPtJEmb : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskDeltaPtJEmb();
  AliAnalysisTaskDeltaPtJEmb(const char *name);
  virtual ~AliAnalysisTaskDeltaPtJEmb() {;}

  void                        UserCreateOutputObjects();

  void SetMinFractionShared(Double_t f)                         { fMinFractionShared = f   ; }


 protected:
  void                        ExecOnce()                                                                                    ;
  Bool_t                      FillHistograms()                                                                              ;

  AliJetContainer            *fJetsCont;                   //!Jets
  AliParticleContainer       *fTracksCont;                 //!Tracks
  AliClusterContainer        *fCaloClustersCont;           //!Clusters  


  Double_t                    fMinFractionShared;          // only fill histos for jets if shared fraction larger than X

  // Jet embedding
  TH3                       **fHistEmbJetsPtArea;          //!Pt vs. area of embedded jets
  TH3                       **fHistEmbJetsCorrPtArea;      //!Pt-rho*A vs. area of embedded jets
  TH2                       **fHistEmbPartPtvsJetPt;       //!MC jet pt total jet pt
  TH2                       **fHistEmbPartPtvsJetCorrPt;   //!MC jet pt total jet pt - rho*A
  TH2                       **fHistJetPtvsJetCorrPt;       //!Pt vs jet pt - rho*A
  TH1                       **fHistDistLeadPart2JetAxis;   //!Distance between leading particle and jet axis
  TH2                       **fHistEmbBkgArea;             //!Pt(embjet) - Pt(embtrack) vs. area of embedded jets
  TH2                       **fHistRhoVSEmbBkg;            //!Area(embjet) * rho vs. Pt(embjet) - Pt(embtrack)
  TH2                       **fHistDeltaPtEmbArea;         //!deltaPt = Pt(embjet) - Area(embjet) * rho - Pt(embtrack) vs. Area(embjet)
  TH2                       **fHistDeltaPtEmbvsEP;         //!deltaPt = Pt(embjet) - Area(embjet) * rho - Pt(embtrack) vs. event plane
  TH2                       **fHistDeltaPtEmbPtProbe;      //!deltaPt = Pt(embjet) - Area(embjet) * rho - Pt(embtrack) vs. Pt(probe)
  TH2                        *fHistEmbJetsPhiEta;          //!Phi-Eta distribution of embedded jets
  TH2                        *fHistLeadPartPhiEta;         //!Phi-Eta distribution of the leading particle of embedded jets

 private:
  AliAnalysisTaskDeltaPtJEmb(const AliAnalysisTaskDeltaPtJEmb&);            // not implemented
  AliAnalysisTaskDeltaPtJEmb &operator=(const AliAnalysisTaskDeltaPtJEmb&); // not implemented

  ClassDef(AliAnalysisTaskDeltaPtJEmb, 1) // deltaPt jet embedding analysis task
};
#endif
