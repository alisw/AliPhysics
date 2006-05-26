/* $Id$ */

#ifndef ALIDNDETAANALYSISESDSELECTOR_H
#define ALIDNDETAANALYSISESDSELECTOR_H

#include "AlidNdEtaAnalysisSelector.h"

class AliESDtrackCuts;
class dNdEtaCorrection;

class AlidNdEtaAnalysisESDSelector : public AlidNdEtaAnalysisSelector {
  public:
    AlidNdEtaAnalysisESDSelector();
    virtual ~AlidNdEtaAnalysisESDSelector();

    virtual void    SlaveBegin(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);

 protected:
    virtual void WriteObjects();

    AliESDtrackCuts*  fEsdTrackCuts;     // Object containing the parameters of the esd track cuts
    dNdEtaCorrection* fdNdEtaCorrection; // correction map

 private:

  ClassDef(AlidNdEtaAnalysisESDSelector, 0);
};

#endif
