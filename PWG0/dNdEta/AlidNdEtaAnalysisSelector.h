#ifndef ALIDNDETAANALYSISSELECTOR_H
#define ALIDNDETAANALYSISSELECTOR_H

#include "AliSelector.h"

class AliESDtrackCuts;
class dNdEtaCorrection;
class dNdEtaAnalysis;

class AlidNdEtaAnalysisSelector : public AliSelector {
  public:
    AlidNdEtaAnalysisSelector(TTree *tree=0);
    virtual ~AlidNdEtaAnalysisSelector();

    virtual void    SlaveBegin(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:
  AliESDtrackCuts*  fEsdTrackCuts;          // Object containing the parameters of the esd track cuts

  dNdEtaAnalysis* fdNdEtaAnalysis;      // contains the intermediate histograms (on each slave)

  dNdEtaCorrection* fdNdEtaCorrection; // correction map
  dNdEtaAnalysis* fdNdEtaAnalysisFinal; // contains the final histograms

 private:

  ClassDef(AlidNdEtaAnalysisSelector, 0);
};

#endif
