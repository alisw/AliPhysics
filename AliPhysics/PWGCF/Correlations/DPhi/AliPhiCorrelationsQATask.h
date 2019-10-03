/* $Id: AliPhiCorrelationsQATask.h 44748 2010-10-29 18:27:12Z jgrosseo $ */

#ifndef AliPhiCorrelationsQATask_H
#define AliPhiCorrelationsQATask_H

#include "AliAnalysisTaskSE.h"

class AliESDtrackCuts;
class TH1F;
class TH2F;
#include <THn.h>
class TH3F;
class AliESDEvent;
class TH1D;

class AliPhiCorrelationsQATask : public AliAnalysisTaskSE {
  public:
    AliPhiCorrelationsQATask(const char* opt = "");
    virtual ~AliPhiCorrelationsQATask();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t*);
    virtual void   Terminate(Option_t*);
    
    void SetUseUncheckedCentrality() { fUseUncheckedCentrality = kTRUE; }
    void SetTrackCuts(AliESDtrackCuts* cuts, AliESDtrackCuts* cuts2 ) { fEsdTrackCuts = cuts; fEsdTrackCuts2 = cuts2; }
    void SetTPCOnly(Bool_t flag) { fTPCOnly = flag; }
    
    void SetOption(const char* opt) { fOption = opt; }

 protected:
    TList* fOutput;                         //! list send on output slot 0

    TString fOption;                        // option string
    
    Bool_t fTPCOnly;                        // tpc only track cuts
    
    AliESDtrackCuts* fEsdTrackCuts;         // Object containing the parameters of the esd track cuts
    AliESDtrackCuts* fEsdTrackCuts2;        // Object containing the parameters of the esd track cuts
    AliESDtrackCuts* fCheckITS;             // Object containing the parameters of the esd track cuts
    AliESDtrackCuts* fGlobalTracks;         // Object containing the parameters of the esd track cuts
    
    TH2F* fCentralityCorrelation;           // correlation of SPD and V0 centrality estimators
    THnF* fDCAPrimaries;                    // DCA distribution of primaries
    THnF* fDCASecondaries;                  // DCA distribution of secondaries
    
    Bool_t fUseUncheckedCentrality;         // for MC!

 private:
    AliPhiCorrelationsQATask(const AliPhiCorrelationsQATask&);
    AliPhiCorrelationsQATask& operator=(const AliPhiCorrelationsQATask&);

  ClassDef(AliPhiCorrelationsQATask, 1);
};

#endif
