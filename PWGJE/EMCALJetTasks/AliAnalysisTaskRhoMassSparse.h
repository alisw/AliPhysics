#ifndef ALIANALYSISTASKRHOMASSSPARSE_H
#define ALIANALYSISTASKRHOMASSSPARSE_H

// $Id$

#include "AliAnalysisTaskRhoMassBase.h"

class AliAnalysisTaskRhoMassSparse : public AliAnalysisTaskRhoMassBase {

 public:
  AliAnalysisTaskRhoMassSparse();
  AliAnalysisTaskRhoMassSparse(const char *name, Bool_t histo=kFALSE);
  virtual ~AliAnalysisTaskRhoMassSparse() {}

  enum JetRhoMassType {
    kMd     = 0,            //rho_m from arXiv:1211.2811
    kMdP    = 1,            //rho_m using P instead of pT
    kMd4    = 2             //rho_m using addition of 4-vectors
  };

  void             UserCreateOutputObjects();

  void             SetExcludeLeadJets(UInt_t n)     { fNExclLeadJets  = n   ; }
  void             SetRhoMassType(JetRhoMassType t) { fJetRhoMassType = t   ; }
  void             SetPionMassForClusters(Bool_t b) { fPionMassClusters = b ; }
  void             SetRhoCMS(Bool_t cms)           { fRhoCMS = cms ; }
  Bool_t           IsJetOverlapping(AliEmcalJet* jet1, AliEmcalJet* jet2);
  Bool_t           IsJetSignal(AliEmcalJet* jet1);


 protected:
  Bool_t           Run();

  Double_t         GetSumMConstituents(AliEmcalJet *jet);
  Double_t         GetSumPtConstituents(AliEmcalJet *jet);
  Double_t         GetMd(AliEmcalJet *jet);

  UInt_t           fNExclLeadJets;                 // number of leading jets to be excluded from the median calculation
  Bool_t           fRhoCMS;                        // flag to run CMS method
  JetRhoMassType   fJetRhoMassType;                // method for rho_m calculation
  Bool_t           fPionMassClusters;              // assume pion mass for clusters

  TH2F            *fHistMdAreavsCent;              //! Md/Area vs cent for all kt clusters
  TH2F            *fHistOccCorrvsCent;             //!occupancy correction vs. centrality

  AliAnalysisTaskRhoMassSparse(const AliAnalysisTaskRhoMassSparse&);             // not implemented
  AliAnalysisTaskRhoMassSparse& operator=(const AliAnalysisTaskRhoMassSparse&);  // not implemented
  
  ClassDef(AliAnalysisTaskRhoMassSparse, 1); // Rho_m task
};
#endif
