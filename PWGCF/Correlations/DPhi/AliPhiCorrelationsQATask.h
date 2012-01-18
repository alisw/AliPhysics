/* $Id: AliPhiCorrelationsQATask.h 44748 2010-10-29 18:27:12Z jgrosseo $ */

#ifndef AliPhiCorrelationsQATask_H
#define AliPhiCorrelationsQATask_H

#include "AliAnalysisTaskSE.h"

class AliESDtrackCuts;
class TH1F;
class TH2F;
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
    
    void SetOption(const char* opt) { fOption = opt; }

 protected:
    TList* fOutput;                         //! list send on output slot 0

    TString fOption;                        // option string
    
    AliESDtrackCuts* fEsdTrackCuts;         // Object containing the parameters of the esd track cuts
    AliESDtrackCuts* fEsdTrackCuts2;        // Object containing the parameters of the esd track cuts
    AliESDtrackCuts* fCheckITS;             // Object containing the parameters of the esd track cuts
    AliESDtrackCuts* fGlobalTracks;         // Object containing the parameters of the esd track cuts
    
    TH2F* fCentralityCorrelation;           // correlation of SPD and V0 centrality estimators
    TH2F* fDCAPrimaries;                    // DCA distribution of primaries
    TH2F* fDCASecondaries;                  // DCA distribution of secondaries
    
    Bool_t fUseUncheckedCentrality;         // for MC!

 private:
    AliPhiCorrelationsQATask(const AliPhiCorrelationsQATask&);
    AliPhiCorrelationsQATask& operator=(const AliPhiCorrelationsQATask&);

  ClassDef(AliPhiCorrelationsQATask, 1);
};

#endif
