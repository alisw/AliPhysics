/* $Id$ */

#ifndef ALIDNDETAANALYSISSELECTORMC_H
#define ALIDNDETAANALYSISSELECTORMC_H

#include "AliSelectorRL.h"

class TH3F;
class TH1F;
class dNdEtaAnalysis;
class AliESDtrackCuts;

class AlidNdEtaAnalysisMCSelector : public AliSelectorRL {
  public:
    AlidNdEtaAnalysisMCSelector();
    virtual ~AlidNdEtaAnalysisMCSelector();

    virtual void    Begin(TTree *tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual void    SlaveTerminate();
    virtual void    Init(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    Terminate();

  protected:
    void ReadUserObjects(TTree* tree);

  private:
    AliESDtrackCuts* fEsdTrackCuts;     // Object containing the parameters of the esd track cuts, needed for fdNdEtaAnalysisTracks

    dNdEtaAnalysis* fdNdEtaAnalysis;        // contains the dndeta from the full sample
    dNdEtaAnalysis* fdNdEtaAnalysisTr;      // contains the dndeta from the triggered events
    dNdEtaAnalysis* fdNdEtaAnalysisTrVtx;   // contains the dndeta from the triggered events with vertex
    dNdEtaAnalysis* fdNdEtaAnalysisTracks;  // contains the dndeta from the triggered events with vertex counted from the mc particles associated to the tracks (comparing this to the raw values from the esd shows the effect of the detector resolution)

    // the following are control histograms to check the dNdEtaAnalysis class
    TH3F* fVertex;     // vertex of counted particles
    TH1F* fPartEta[3]; // counted particles as function of eta (full vertex range, below 0 range, above 0 range)
    TH1F* fPartPt;     // counted particles as function of pt
    TH1F* fEvents;     // events counted as function of vtx

    AlidNdEtaAnalysisMCSelector(const AlidNdEtaAnalysisMCSelector&);
    AlidNdEtaAnalysisMCSelector& operator=(const AlidNdEtaAnalysisMCSelector&);

    ClassDef(AlidNdEtaAnalysisMCSelector, 0);
};

#endif
