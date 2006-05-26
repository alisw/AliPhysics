/* $Id$ */

#ifndef ALIDNDETAANALYSISSELECTOR_H
#define ALIDNDETAANALYSISSELECTOR_H

#include "AliSelector.h"

class dNdEtaAnalysis;

class AlidNdEtaAnalysisSelector : public AliSelector {
  public:
    AlidNdEtaAnalysisSelector();
    virtual ~AlidNdEtaAnalysisSelector();

    virtual void    SlaveBegin(TTree *tree);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:
    virtual void WriteObjects();

    dNdEtaAnalysis* fdNdEtaAnalysis;      // contains the intermediate histograms (on each slave)
    dNdEtaAnalysis* fdNdEtaAnalysisFinal; // contains the final histograms

 private:
    ClassDef(AlidNdEtaAnalysisSelector, 0);
};

#endif
