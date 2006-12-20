/* $Id$ */

#ifndef ALIDNDETAANALYSISESDSELECTOR_H
#define ALIDNDETAANALYSISESDSELECTOR_H

#include "AliSelectorRL.h"

class AliESDtrackCuts;
class dNdEtaAnalysis;
class AlidNdEtaCorrection;
class TH1F;

class AlidNdEtaAnalysisESDSelector : public AliSelectorRL {
  public:
    AlidNdEtaAnalysisESDSelector();
    virtual ~AlidNdEtaAnalysisESDSelector();

    virtual void    Begin(TTree* tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:
    void ReadUserObjects(TTree* tree);

    dNdEtaAnalysis* fdNdEtaAnalysis;        // contains the uncorrected histograms
    TH1F*           fMult;                  // raw multiplicity histogram (control histogram)

    AliESDtrackCuts*  fEsdTrackCuts;     // Object containing the parameters of the esd track cuts

 private:
    AlidNdEtaAnalysisESDSelector(const AlidNdEtaAnalysisESDSelector&);
    AlidNdEtaAnalysisESDSelector& operator=(const AlidNdEtaAnalysisESDSelector&);

  ClassDef(AlidNdEtaAnalysisESDSelector, 0);
};

#endif
