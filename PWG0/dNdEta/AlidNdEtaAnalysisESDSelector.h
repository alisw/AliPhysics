/* $Id$ */

#ifndef ALIDNDETAANALYSISESDSELECTOR_H
#define ALIDNDETAANALYSISESDSELECTOR_H

#include "AliSelector.h"

class AliESDtrackCuts;
class dNdEtaAnalysis;
class AlidNdEtaCorrection;

class AlidNdEtaAnalysisESDSelector : public AliSelector {
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

    dNdEtaAnalysis* fdNdEtaAnalysisMBVtx;   // contains the histograms for the triggered events with vertex
    dNdEtaAnalysis* fdNdEtaAnalysisMB;      // contains the histograms corrected with vtx recon eff
    dNdEtaAnalysis* fdNdEtaAnalysis;        // contains the histograms corrected with vtx recon eff and trigger bias eff

    AliESDtrackCuts*  fEsdTrackCuts;     // Object containing the parameters of the esd track cuts

    AlidNdEtaCorrection* fdNdEtaCorrection; // correction maps

 private:

  ClassDef(AlidNdEtaAnalysisESDSelector, 0);
};

#endif
