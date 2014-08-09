#ifndef ALIANALYSISTASKSAJF_H
#define ALIANALYSISTASKSAJF_H

class TH2;
class THnSparse;

#include <TH3F.h>

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskSAJF : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskSAJF();
  AliAnalysisTaskSAJF(const char *name);
  virtual ~AliAnalysisTaskSAJF() {;}

  void                        UserCreateOutputObjects();

  void                        SetHistoType(Int_t t) { fHistoType = t; }

 protected:
  void                        AllocateTHX();
  void                        AllocateTHnSparse();

  Bool_t                      FillHistograms();
  void                        FillJetHisto(Double_t cent, Double_t ep, Double_t eta, Double_t phi, Double_t pt, Double_t MCpt, Double_t corrpt, Double_t area, 
					   Double_t NEF, Double_t z, Int_t n, Double_t leadingpt);

  Int_t                       fHistoType;                      // histogram type (0=TH2, 1=THnSparse)

  TH2                       **fHistTracksJetPt;             //!Track pt vs. jet pt
  TH2                       **fHistClustersJetPt;           //!Cluster pt vs. jet pt
  TH2                       **fHistTracksPtDist;            //!Track pt vs. distance form jet axis
  TH2                       **fHistClustersPtDist;          //!Cluster pt vs. distance form jet axis

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
  AliAnalysisTaskSAJF(const AliAnalysisTaskSAJF&);            // not implemented
  AliAnalysisTaskSAJF &operator=(const AliAnalysisTaskSAJF&); // not implemented

  ClassDef(AliAnalysisTaskSAJF, 17) // jet analysis task
};
#endif
