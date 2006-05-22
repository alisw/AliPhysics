#ifndef ALIDNDETAANALYSISSELECTORMC_H
#define ALIDNDETAANALYSISSELECTORMC_H

#include "AlidNdEtaAnalysisSelector.h"

class AlidNdEtaAnalysisMCSelector : public AlidNdEtaAnalysisSelector {
  public:
    AlidNdEtaAnalysisMCSelector(TTree *tree=0);
    virtual ~AlidNdEtaAnalysisMCSelector();

    virtual void    Init(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);

 protected:

 private:

  ClassDef(AlidNdEtaAnalysisMCSelector, 0);
};

#endif
