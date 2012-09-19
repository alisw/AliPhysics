#ifndef ALIANALYSISTASKSAJF_H
#define ALIANALYSISTASKSAJF_H

// $Id$

class TClonesArray;
class TString;
class TH1F;
class TH2F;
class TH3F;
class AliRhoParameter;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskSAJF : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskSAJF();
  AliAnalysisTaskSAJF(const char *name);
  virtual ~AliAnalysisTaskSAJF();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  void                        SetLeadingHadronType(Int_t t)     { fLeadingHadronType = t; }

 protected:
  Bool_t                      FillHistograms()                                              ;
  Int_t                       DoJetLoop()                                                   ;
  Float_t*                    GenerateFixedBinArray(Int_t n, Float_t min, Float_t max) const;

  Int_t                       fLeadingHadronType;          // 0 = charged, 1 = neutral, 2 = both

  // General histograms
  TH1F                       *fHistEvents[4];              //!Events accepted/rejected
  TH1F                       *fHistLeadingJetPt[4];        //!Leading jet pt spectrum
  TH1F                       *fHist2LeadingJetPt[4];       //!Second leading jet pt spectrum
  TH1F                       *fHistLeadingJetCorrPt[4];    //!Corrected leading jet pt spectrum
  TH2F                       *fHistRhoVSleadJetPt[4];      //!Area(leadjet) * rho vs. leading jet pt
  TH2F                       *fNjetsVsCent;                //!No. of jets vs. centrality

  // Inclusive jets histograms
  TH3F                       *fHistJetPhiEta[4];           //!Phi-Eta distribution of jets
  TH3F                       *fHistJetsPtArea[4];          //!Jet pt vs. area
  TH3F                       *fHistJetsCorrPtArea[4];      //!Jet corr pt vs. area
  TH3F                       *fHistJetsNEFvsPt[4];         //!Jet neutral energy fraction vs. jet pt
  TH3F                       *fHistJetsZvsPt[4];           //!Constituent Pt over Jet Pt ratio vs. jet pt
  TH2F                       *fHistConstituents[4];        //!x axis = constituents pt; y axis = no. of constituents
  TH2F                       *fHistTracksJetPt[4];         //!Track pt vs. jet pt
  TH2F                       *fHistClustersJetPt[4];       //!Cluster pt vs. jet pt
  TH3F                       *fHistJetNconstVsPt[4];       //!Jet no. of constituents vs. pt

 private:
  AliAnalysisTaskSAJF(const AliAnalysisTaskSAJF&);            // not implemented
  AliAnalysisTaskSAJF &operator=(const AliAnalysisTaskSAJF&); // not implemented

  ClassDef(AliAnalysisTaskSAJF, 13) // jet analysis task
};
#endif
