#ifndef ALIANALYSISTASKSAJF_H
#define ALIANALYSISTASKSAJF_H

// $Id$

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

  // Inclusive jets histograms
  THnSparse                  *fHistJetObservables;             //!Jet-wise observables

  // TH2/TH3 versions
  TH3                        *fHistJetPtEtaPhi[4];             //!Jet Pt vs. Eta vs. Phi
  TH2                        *fHistJetPtArea[4];               //!Jet Pt vs. Area
  TH2                        *fHistJetPtEP[4];                 //!Jet Pt vs. event plane
  TH2                        *fHistJetPtNEF[4];                //!Jet Pt vs. neutral energy fraction
  TH2                        *fHistJetPtZ[4];                  //!Jet Pt vs. z
  TH2                        *fHistJetPtLeadingPartPt[4];      //!Jet Pt vs. leading particle pt
  TH3                        *fHistJetCorrPtEtaPhi[4];         //!Jet corrPt vs. Eta vs. Phi
  TH2                        *fHistJetCorrPtArea[4];           //!Jet corrPt vs. Area
  TH2                        *fHistJetCorrPtEP[4];             //!Jet corrPt vs. event plane
  TH2                        *fHistJetCorrPtNEF[4];            //!Jet corrPt vs. neutral energy fraction
  TH2                        *fHistJetCorrPtZ[4];              //!Jet corrPt vs. z
  TH2                        *fHistJetCorrPtLeadingPartPt[4];  //!Jet corrPt vs. leading particle pt
  TH2                        *fHistJetPtCorrPt[4];             //!Jet Pt vs. corrPt
  TH2                        *fHistJetPtMCPt[4];               //!Jet Pt vs. MCPt
  TH2                        *fHistJetMCPtCorrPt[4];           //!Jet MCPt vs. corrPt

  TH2                        *fHistTracksJetPt[4];             //!Track pt vs. jet pt
  TH2                        *fHistClustersJetPt[4];           //!Cluster pt vs. jet pt
  TH2                        *fHistTracksPtDist[4];            //!Track pt vs. distance form jet axis
  TH2                        *fHistClustersPtDist[4];          //!Cluster pt vs. distance form jet axis

 private:
  AliAnalysisTaskSAJF(const AliAnalysisTaskSAJF&);            // not implemented
  AliAnalysisTaskSAJF &operator=(const AliAnalysisTaskSAJF&); // not implemented

  ClassDef(AliAnalysisTaskSAJF, 17) // jet analysis task
};
#endif
