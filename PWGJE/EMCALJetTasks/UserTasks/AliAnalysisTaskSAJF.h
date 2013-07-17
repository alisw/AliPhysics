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
  void                        DoJetLoop(const Int_t nacc, const Int_t* sortedJets)          ;

  // General histograms
  THnSparse                  *fHistEventObservables;       //!Event-wise observables (centrality, event plane, rho, number of jets, leading jet pt, leading jet corr pt, 
                                                           // leading jet area)

  // Inclusive jets histograms
  THnSparse                  *fHistJetObservables;         //!Jet-wise observables (centrality, event plane, phi, eta, pt, MCPt, area, corr pt, NEF, Z, no. of constituents
                                                           // leading hadron pt)
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
