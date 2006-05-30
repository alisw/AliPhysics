/* $Id$ */

#ifndef ALIDNDETAANALYSISSELECTOR_H
#define ALIDNDETAANALYSISSELECTOR_H

#include "AliSelector.h"

class dNdEtaAnalysis;
class dNdEtaCorrection;

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
    dNdEtaCorrection* fdNdEtaCorrection;  // correction map

 private:
    ClassDef(AlidNdEtaAnalysisSelector, 0);
};

#endif
