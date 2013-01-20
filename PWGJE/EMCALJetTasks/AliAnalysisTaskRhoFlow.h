#ifndef ALIANALYSISTASKRHOFLOW_H
#define ALIANALYSISTASKRHOFLOW_H

// $Id$

#include "AliAnalysisTaskRho.h"

class AliAnalysisTaskRhoFlow : public AliAnalysisTaskRho {

 public:
  AliAnalysisTaskRhoFlow();
  AliAnalysisTaskRhoFlow(const char *name);
  virtual ~AliAnalysisTaskRhoFlow() {}

  void             UserCreateOutputObjects();

 protected:
  Bool_t           Run();
  Bool_t           FillHistograms();

  AliAnalysisTaskRhoFlow(const AliAnalysisTaskRhoFlow&);             // not implemented
  AliAnalysisTaskRhoFlow& operator=(const AliAnalysisTaskRhoFlow&);  // not implemented

  Double_t               fRhoNearSide;                   //!Rho in the near side
  Double_t               fRhoAwaySide;                   //!Rho in the away side
  Double_t               fRhoPerpSide1;                  //!Rho in the perpendicular side 1
  Double_t               fRhoPerpSide2;                  //!Rho in the perpendicular side 2
  TH2F                  *fHistRhoNearVsCent;             //!Near side rho vs. centrality
  TH2F                  *fHistDeltaRhoNearVsCent;        //!Rho - rho_near vs. centrality
  TH2F                  *fHistRhoAwayVsCent;             //!Away side rho vs. centrality
  TH2F                  *fHistDeltaRhoAwayVsCent;        //!Rho - rho_away vs. centrality
  TH2F                  *fHistRhoPerp1VsCent;            //!Perpendicualr side 1 rho vs. centrality
  TH2F                  *fHistDeltaRhoPerp1VsCent;       //!Rho - rho_perp1 vs. centrality
  TH2F                  *fHistRhoPerp2VsCent;            //!Perpendicualr side 2 rho vs. centrality
  TH2F                  *fHistDeltaRhoPerp2VsCent;       //!Rho - rho_perp2 vs. centrality
  
  ClassDef(AliAnalysisTaskRhoFlow, 1); // Rho task for flow bias study
};
#endif
