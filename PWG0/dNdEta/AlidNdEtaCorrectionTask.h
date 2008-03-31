/* $Id$ */

#ifndef AlidNdEtaCorrectionTask_H
#define AlidNdEtaCorrectionTask_H

#include "AliAnalysisTask.h"
#include <TString.h>
#include "AliPWG0Helper.h"

class AliESDtrackCuts;
class dNdEtaAnalysis;
class AlidNdEtaCorrection;
class TH1F;
class AliESDEvent;
class TParticlePDG;
class TH2F;
class TProfile;

class AlidNdEtaCorrectionTask : public AliAnalysisTask {
  public:
    AlidNdEtaCorrectionTask(const char* opt = "");
    virtual ~AlidNdEtaCorrectionTask();

    virtual void   ConnectInputData(Option_t *);
    virtual void   CreateOutputObjects();
    virtual void   Exec(Option_t*);
    virtual void   Terminate(Option_t *);

    void SetTrackCuts(AliESDtrackCuts* cuts) { fEsdTrackCuts = cuts; }
    void SetAnalysisMode(AliPWG0Helper::AnalysisMode mode) { fAnalysisMode = mode; }
    void SetOnlyPrimaries(Bool_t flag = kTRUE) { fOnlyPrimaries = flag; }

 protected:
    Bool_t SignOK(TParticlePDG* particle);

    AliESDEvent *fESD;               //! ESD object
    TList* fOutput;                  //! list send on output slot 0

    TString fOption;                 // option string
    AliPWG0Helper::AnalysisMode fAnalysisMode;    // detector that is used for analysis
    Int_t fSignMode;                 // if 0 process all particles, if +-1 process only particles with that sign
    Bool_t fOnlyPrimaries;           // only process primaries (syst. studies)

    AliESDtrackCuts*  fEsdTrackCuts;             // Object containing the parameters of the esd track cuts

    AlidNdEtaCorrection* fdNdEtaCorrection;      //! contains the intermediate histograms (on each slave)

    dNdEtaAnalysis* fdNdEtaAnalysisMC;           //! analysis from MC (only triggered, vertex events)
    dNdEtaAnalysis* fdNdEtaAnalysisESD;          //! analysis from ESD (not yet corrected!)

    // control histograms
    TH1F* fPIDParticles;                         //! pid of primary particles
    TH1F* fPIDTracks;                            //! pid of reconstructed tracks
 
    TH2F* fVertexCorrelation;                    //! ESD z-vtx vs MC z-vtx
    TProfile* fVertexProfile;                    //! Profile of MC z-vtx - ESD z-vtx vs. MC z-vtx
    TH1F* fVertexShiftNorm;                      //! (MC z-vtx - ESD z-vtx) / (sigma_ESD-z-vtx) histogrammed

    TH2F* fEtaCorrelation;                       //! ESD eta vs MC eta
    TProfile* fEtaProfile;                       //! Profile of MC eta - ESD eta vs. MC eta

    // histograms for systematic studies (must be enabled with option)

    TH1F* fSigmaVertexTracks;                    //! (accepted tracks) vs (n of sigma to vertex cut)
    TH1F* fSigmaVertexPrim;                      //! (accepted primaries) vs (n of sigma to vertex cut)
                                                 // enable with option: sigma-vertex

    TH1F* fMultAll; //! primary particles  in |eta| < 1 and pT > 0.2 in all events
    TH1F* fMultTr; //! primary particles  in |eta| < 1 and pT > 0.2 in triggered events
    TH1F* fMultVtx; //! primary particles  in |eta| < 1 and pT > 0.2 in triggered events with vertex

    TH1F* fDeltaPhi[8]; //! delta phi of primaries, secondaries, other (= unclear cases)

    AlidNdEtaCorrection* fdNdEtaCorrectionProcessType[3]; //! correction for specific process type (ND, SD, DD)
                                                          // enable with option: process-types

 private:
    AlidNdEtaCorrectionTask(const AlidNdEtaCorrectionTask&);
    AlidNdEtaCorrectionTask& operator=(const AlidNdEtaCorrectionTask&);

  ClassDef(AlidNdEtaCorrectionTask, 1);
};

#endif
