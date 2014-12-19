#ifndef ALIANALYSISTASKRHOMASS_H
#define ALIANALYSISTASKRHOMASS_H

// $Id$

#include "AliAnalysisTaskRhoMassBase.h"

class AliAnalysisTaskRhoMass : public AliAnalysisTaskRhoMassBase {

 public:
  AliAnalysisTaskRhoMass();
  AliAnalysisTaskRhoMass(const char *name, Bool_t histo=kFALSE);
  virtual ~AliAnalysisTaskRhoMass() {}

  enum JetRhoMassType {
    kMd     = 0,            //rho_m from arXiv:1211.2811
    kMdP    = 1,            //rho_m using P instead of pT
    kMd4    = 2             //rho_m using addition of 4-vectors
  };

  void             UserCreateOutputObjects();

  void             SetExcludeLeadJets(UInt_t n)     { fNExclLeadJets  = n   ; }
  void             SetRhoMassType(JetRhoMassType t) { fJetRhoMassType = t   ; }
  void             SetPionMassForClusters(Bool_t b) { fPionMassClusters = b ; }

 protected:
  Bool_t           Run();

  Double_t         GetSumMConstituents(AliEmcalJet *jet);
  Double_t         GetSumPtConstituents(AliEmcalJet *jet);
  Double_t         GetMd(AliEmcalJet *jet);

  UInt_t           fNExclLeadJets;                 // number of leading jets to be excluded from the median calculation
  JetRhoMassType   fJetRhoMassType;                // method for rho_m calculation
  Bool_t           fPionMassClusters;              // assume pion mass for clusters

  TH2F            *fHistMdAreavsCent;              //! Md/Area vs cent for all kt clusters

  AliAnalysisTaskRhoMass(const AliAnalysisTaskRhoMass&);             // not implemented
  AliAnalysisTaskRhoMass& operator=(const AliAnalysisTaskRhoMass&);  // not implemented
  
  ClassDef(AliAnalysisTaskRhoMass, 2); // Rho_m task
};
#endif
