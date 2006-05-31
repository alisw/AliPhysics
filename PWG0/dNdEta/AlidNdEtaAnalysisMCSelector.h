/* $Id$ */

#ifndef ALIDNDETAANALYSISSELECTORMC_H
#define ALIDNDETAANALYSISSELECTORMC_H

#include "AlidNdEtaAnalysisSelector.h"

class TH3F;
class TH1F;

class AlidNdEtaAnalysisMCSelector : public AlidNdEtaAnalysisSelector {
  public:
    AlidNdEtaAnalysisMCSelector();
    virtual ~AlidNdEtaAnalysisMCSelector();

    virtual void    Init(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    Terminate();

 protected:

 private:
    TH3F* fVertex;  //! vertex of counted particles
    TH1F* fPartEta; //! counted particles as function of eta
    Int_t fEvents;  //! number of processed events

    ClassDef(AlidNdEtaAnalysisMCSelector, 0);
};

#endif
