#ifndef ALIANALYSISTASKSAJF_H
#define ALIANALYSISTASKSAJF_H

// $Id$

#include <TH3F.h>

class TH2;
class THnSparse;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskSAJF : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskSAJF();
  AliAnalysisTaskSAJF(const char *name);
  virtual ~AliAnalysisTaskSAJF() {;}

  void                        UserCreateOutputObjects();

 protected:
  Bool_t                      FillHistograms()                                              ;
  void                        FillJetHisto(Double_t cent, Double_t ep, Double_t eta, Double_t phi, Double_t pt, Double_t MCpt, Double_t corrpt, Double_t area, 
					   Double_t NEF, Double_t z, Int_t n, Double_t leadingpt);

  // Inclusive jets histograms
  THnSparse                  *fHistJetObservables;         //!Jet-wise observables
  TH2                        *fHistTracksJetPt[4];         //!Track pt vs. jet pt
  TH2                        *fHistClustersJetPt[4];       //!Cluster pt vs. jet pt
  TH2                        *fHistTracksPtDist[4];        //!Track pt vs. distance form jet axis
  TH2                        *fHistClustersPtDist[4];      //!Cluster pt vs. distance form jet axis

 private:
  AliAnalysisTaskSAJF(const AliAnalysisTaskSAJF&);            // not implemented
  AliAnalysisTaskSAJF &operator=(const AliAnalysisTaskSAJF&); // not implemented

  ClassDef(AliAnalysisTaskSAJF, 15) // jet analysis task
};
#endif
