/* $Id$ */

#ifndef ALIDNDETAANALYSISESDSELECTOR_H
#define ALIDNDETAANALYSISESDSELECTOR_H

#include "AliSelector.h"

class AliESDtrackCuts;
class dNdEtaAnalysis;

class AlidNdEtaAnalysisESDSelector : public AliSelector {
  public:
    AlidNdEtaAnalysisESDSelector();
    virtual ~AlidNdEtaAnalysisESDSelector();

    virtual void    SlaveBegin(TTree *tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:
    dNdEtaAnalysis* fdNdEtaAnalysis;      // contains the target histograms
    AliESDtrackCuts*  fEsdTrackCuts;     // Object containing the parameters of the esd track cuts

 private:

  ClassDef(AlidNdEtaAnalysisESDSelector, 0);
};

#endif
