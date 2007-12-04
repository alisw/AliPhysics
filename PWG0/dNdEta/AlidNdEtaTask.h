/* $Id$ */

#ifndef AlidNdEtaTask_H
#define AlidNdEtaTask_H

#include "AliAnalysisTask.h"
#include <TString.h>

class AliESDtrackCuts;
class dNdEtaAnalysis;
class TH1F;
class TH3F;
class AliESDEvent;

class AlidNdEtaTask : public AliAnalysisTask {
  public:
    enum AnalysisMethod { kSPD = 0, kTPC };

    AlidNdEtaTask(const char* opt = "");
    virtual ~AlidNdEtaTask();

    virtual void   ConnectInputData(Option_t *);
    virtual void   CreateOutputObjects();
    virtual void   Exec(Option_t*);
    virtual void   Terminate(Option_t*);

    void SetTrackCuts(AliESDtrackCuts* cuts) { fEsdTrackCuts = cuts; }
    void SetAnalysisMode(AnalysisMethod mode) { fAnalysisMode = mode; }
    void SetReadMC(Bool_t flag = kTRUE) { fReadMC = flag; }

 protected:
    AliESDEvent *fESD;    //! ESD object
    TList* fOutput;                  //! list send on output slot 0

    TString fOption;      // option string
    AnalysisMethod fAnalysisMode; // detector that is used for analysis

    Bool_t  fReadMC;       // if true reads MC data (to build correlation maps)

    AliESDtrackCuts* fEsdTrackCuts;         // Object containing the parameters of the esd track cuts

    // Gathered from ESD
    dNdEtaAnalysis* fdNdEtaAnalysisESD;     //! contains the dndeta from the ESD
    // control hists
    TH1F* fMult;                            //! raw multiplicity histogram (control histogram)
    TH1F* fMultVtx;                            //! raw multiplicity histogram of evts with vtx (control histogram)
    TH1F* fPartEta[3];            //! counted particles as function of eta (full vertex range, below 0 range, above 0 range)
    TH1F* fEvents;                //! events counted as function of vtx

    // Gathered from MC (when fReadMC is set)
    dNdEtaAnalysis* fdNdEtaAnalysis;        //! contains the dndeta from the full sample
    dNdEtaAnalysis* fdNdEtaAnalysisTr;      //! contains the dndeta from the triggered events
    dNdEtaAnalysis* fdNdEtaAnalysisTrVtx;   //! contains the dndeta from the triggered events with vertex
    dNdEtaAnalysis* fdNdEtaAnalysisTracks;  //! contains the dndeta from the triggered events with vertex counted from the mc particles associated to the tracks (comparing this to the raw values from the esd shows the effect of the detector resolution)
    // the following are control histograms to check the dNdEtaAnalysis class
    TH3F* fVertex;                //! vertex of counted particles
    TH1F* fPartPt;                //! counted particles as function of pt

 private:
    AlidNdEtaTask(const AlidNdEtaTask&);
    AlidNdEtaTask& operator=(const AlidNdEtaTask&);

  ClassDef(AlidNdEtaTask, 1);
};

#endif
