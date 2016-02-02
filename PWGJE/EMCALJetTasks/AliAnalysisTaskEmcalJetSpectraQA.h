#ifndef ALIANALYSISTASKEMCALJETSPECTRAQA_H
#define ALIANALYSISTASKEMCALJETSPECTRAQA_H

class TH2;
class TH3;
class THnSparse;

#include <TH3F.h>

#include "AliTLorentzVector.h"
#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskEmcalJetSpectraQA : public AliAnalysisTaskEmcalJet {
 public:

  struct AliEmcalJetInfo : public AliTLorentzVector {
    AliEmcalJetInfo();
    AliEmcalJetInfo(const AliEmcalJet& jet);

    Double_t fArea;
    Double_t fMCPt;
    Double_t fNConstituents;
    Double_t fNEF;
    Double_t fCent;
    Double_t fEP;
    Double_t fCorrPt;
    Double_t fZ;
    Double_t fLeadingPt;
  };

  AliAnalysisTaskEmcalJetSpectraQA();
  AliAnalysisTaskEmcalJetSpectraQA(const char *name);
  virtual ~AliAnalysisTaskEmcalJetSpectraQA() {;}

  void                        UserCreateOutputObjects();

  void                        SetHistoType(Int_t t)               { fHistoType             = t; }
  void                        SetDefaultClusterEnergy(Int_t d)    { fDefaultClusterEnergy  = d; }
  void                        SetJetEPaxis(Bool_t b)              { fJetEPaxis             = b; }

 protected:
  void                        AllocateTHX();
  void                        AllocateTHnSparse();

  Bool_t                      FillHistograms();
  void                        FillJetHisto(const AliEmcalJetInfo& jetInfo);

  Int_t                       fHistoType;                   // histogram type (0=TH2, 1=THnSparse)
  Int_t                       fDefaultClusterEnergy;        // default cluster energy
  Bool_t                      fJetEPaxis;                   // whether a EP-jet axis should be included in the THnSparse

  TH2                       **fHistRejectionReason;         //!Rejection reason vs. jet pt
  TH2                       **fHistTracksJetPt;             //!Track pt vs. jet pt
  TH2                       **fHistClustersJetPt;           //!Cluster pt vs. jet pt
  TH2                       **fHistTracksPtDist;            //!Track pt vs. distance form jet axis
  TH2                       **fHistClustersPtDist;          //!Cluster pt vs. distance form jet axis
  TH3                       **fHistTracksZJetPtJetConst;    //!Track z vs. jet pt vs. no. of jet const
  TH3                       **fHistClustersZJetPtJetConst;  //!Cluster z vs. jet pt vs. no. of jet const
  TH2                        *fHistRhoVsCent;               //!Background vs. centrality

  // Inclusive jets histograms
  THnSparse                  *fHistJetObservables;          //!Jet observables

  // TH2/TH3 versions
  TH3                       **fHistJetPtEtaPhi;             //!Jet Pt vs. Eta vs. Phi
  TH2                       **fHistJetPtArea;               //!Jet Pt vs. Area
  TH2                       **fHistJetPtEP;                 //!Jet Pt vs. event plane
  TH2                       **fHistJetPtNEF;                //!Jet Pt vs. neutral energy fraction
  TH2                       **fHistJetPtZ;                  //!Jet Pt vs. z
  TH2                       **fHistJetPtLeadingPartPt;      //!Jet Pt vs. leading particle pt
  TH3                       **fHistJetCorrPtEtaPhi;         //!Jet corrPt vs. Eta vs. Phi
  TH2                       **fHistJetCorrPtArea;           //!Jet corrPt vs. Area
  TH2                       **fHistJetCorrPtEP;             //!Jet corrPt vs. event plane
  TH2                       **fHistJetCorrPtNEF;            //!Jet corrPt vs. neutral energy fraction
  TH2                       **fHistJetCorrPtZ;              //!Jet corrPt vs. z
  TH2                       **fHistJetCorrPtLeadingPartPt;  //!Jet corrPt vs. leading particle pt
  TH2                       **fHistJetPtCorrPt;             //!Jet Pt vs. corrPt
  TH2                       **fHistJetPtMCPt;               //!Jet Pt vs. MCPt
  TH2                       **fHistJetMCPtCorrPt;           //!Jet MCPt vs. corrPt

 private:
  AliAnalysisTaskEmcalJetSpectraQA(const AliAnalysisTaskEmcalJetSpectraQA&);            // not implemented
  AliAnalysisTaskEmcalJetSpectraQA &operator=(const AliAnalysisTaskEmcalJetSpectraQA&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetSpectraQA, 1) // jet spectra QA task
};
#endif
