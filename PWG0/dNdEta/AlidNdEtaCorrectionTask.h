/* $Id$ */

#ifndef AlidNdEtaCorrectionTask_H
#define AlidNdEtaCorrectionTask_H

#include "AliAnalysisTask.h"
#include <TString.h>
#include "AliPWG0Helper.h"

class AliESDtrackCuts;
class dNdEtaAnalysis;
class AlidNdEtaCorrection;
class TH1;
class TH1F;
class AliESDEvent;
class TParticlePDG;
class TH2F;
class TProfile;

class AlidNdEtaCorrectionTask : public AliAnalysisTask {
  public:
    AlidNdEtaCorrectionTask();
    AlidNdEtaCorrectionTask(const char* opt);
    virtual ~AlidNdEtaCorrectionTask();

    virtual void   ConnectInputData(Option_t *);
    virtual void   CreateOutputObjects();
    virtual void   Exec(Option_t*);
    virtual void   Terminate(Option_t *);

    void SetTrackCuts(AliESDtrackCuts* cuts) { fEsdTrackCuts = cuts; }
    void SetAnalysisMode(AliPWG0Helper::AnalysisMode mode) { fAnalysisMode = mode; }
    void SetOnlyPrimaries(Bool_t flag = kTRUE) { fOnlyPrimaries = flag; }
    void SetTrigger(AliPWG0Helper::Trigger trigger) { fTrigger = trigger; }
    void SetFillPhi(Bool_t flag = kTRUE) { fFillPhi = flag; }
    void SetDeltaPhiCut(Float_t cut) { fDeltaPhiCut = cut; }

    void SetOption(const char* opt) { fOption = opt; }

 protected:
    Bool_t SignOK(TParticlePDG* particle);

    AliESDEvent *fESD;               //! ESD object
    TList* fOutput;                  //! list send on output slot 0

    TString fOption;                 // option string
    AliPWG0Helper::AnalysisMode fAnalysisMode;    // detector that is used for analysis
    AliPWG0Helper::Trigger fTrigger; // trigger used in the analysis
    Bool_t fFillPhi;                           // if true phi is filled as 3rd coordinate in all maps
    Float_t fDeltaPhiCut;                      // cut in delta phi (only SPD)

    Int_t fSignMode;                 // if 0 process all particles, if +-1 process only particles with that sign
    Bool_t fOnlyPrimaries;           // only process primaries (syst. studies)
    Int_t fStatError;                // statistical error evaluation: if set to 1 we only count unique primaries (binomial errors are valid), for 2 all the rest

    AliESDtrackCuts*  fEsdTrackCuts;             // Object containing the parameters of the esd track cuts

    AlidNdEtaCorrection* fdNdEtaCorrection;      //! contains the intermediate histograms (on each slave)

    dNdEtaAnalysis* fdNdEtaAnalysisMC;           //! analysis from MC (only triggered, vertex events)
    dNdEtaAnalysis* fdNdEtaAnalysisESD;          //! analysis from ESD (not yet corrected!)

    // control histograms
    TH1F* fPIDParticles;                         //! pid of primary particles
    TH1F* fPIDTracks;                            //! pid of reconstructed tracks

    TH2F* fVertexCorrelation;                    //! ESD z-vtx vs MC z-vtx
    TH2F* fVertexCorrelationShift;               //! (MC z-vtx - ESD z-vtx) vs MC z-vtx
    TProfile* fVertexProfile;                    //! Profile of MC z-vtx - ESD z-vtx vs. MC z-vtx
    TH1F* fVertexShift;                          //! (MC z-vtx - ESD z-vtx) in +- 10 cm
    TH1F* fVertexShiftNorm;                      //! (MC z-vtx - ESD z-vtx) / (sigma_ESD-z-vtx) histogrammed

    TH2F* fEtaCorrelation;                       //! ESD eta vs MC eta
    TH2F* fEtaCorrelationShift;                  //! (MC eta - ESD eta) vs MC eta
    TProfile* fEtaProfile;                       //! Profile of MC eta - ESD eta vs. MC eta
    TH1F* fEtaResolution;                        //! MC eta - ESD eta in |eta| < 1
    TH2F* fDeltaPhiCorrelation;                  //! delta phi ESD vs. MC

    TH2F* fpTResolution;                         //! (MC pT - ESD pT) / MC pT vs. MC pT in |eta| < 0.9

    AliESDtrackCuts*  fEsdTrackCutsPrim;         //! control histograms for primaries
    AliESDtrackCuts*  fEsdTrackCutsSec;          //! control histograms for secondaries

    // histograms for systematic studies (must be enabled with option)

    TH1* fTemp1;                                 //! temp histogram for quick study of variables
    TH1* fTemp2;                                 //! temp histogram for quick study of variables

    TH1F* fMultAll; //! primary particles  in |eta| < 1 and pT > 0.2 in all events
    TH1F* fMultTr; //! primary particles  in |eta| < 1 and pT > 0.2 in triggered events
    TH1F* fMultVtx; //! primary particles  in |eta| < 1 and pT > 0.2 in triggered events with vertex

    TH2* fDeltaPhi[8]; //! delta phi of primaries, secondaries, other (= unclear cases)

    TH2F* fEventStats;  //! some stats on number of events, see CreateOutputObjects for a detailed definition

    AlidNdEtaCorrection* fdNdEtaCorrectionSpecial[4];   //! correction maps used for systematic studies, may contain:
                                                        // for specific process type (ND, SD, DD), enable with option: process-types
                                                        // for particle species (pi, K, p, rest), enable with: particle-species

 private:
    AlidNdEtaCorrectionTask(const AlidNdEtaCorrectionTask&);
    AlidNdEtaCorrectionTask& operator=(const AlidNdEtaCorrectionTask&);

  ClassDef(AlidNdEtaCorrectionTask, 1);
};

#endif
