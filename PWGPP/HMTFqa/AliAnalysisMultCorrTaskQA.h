#ifndef AliAnalysisMultCorrTaskQA_H
#define AliAnalysisMultCorrTaskQA_H

class TList;
class TH1F;
class TH2F;
class TH3F;
class TProfile;
class TVector3;

class AliESDpid;
class AliESDtrackCuts;
class AliAnalysisUtils;
class AliESDEvent;
class AliPhysicsSelection;

class AliESDtrack;
class AliPPVsMultUtils;
class AliAnalysisFilter;
//#include "TString.h"
//#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisMultCorrTaskQA : public AliAnalysisTaskSE {
public:
    AliAnalysisMultCorrTaskQA();
    AliAnalysisMultCorrTaskQA(const char *name);
    virtual ~AliAnalysisMultCorrTaskQA();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    //virtual void   SetCuts(AliESDtrackCuts *); //cuts=0){fCuts = cuts;}
    virtual void   SetCuts(AliESDtrackCuts* cuts=0){fCuts = cuts;}
    virtual void   SetTrackFilterSKE(AliAnalysisFilter* trackF) {fTrackFilterSKE = trackF;}
    
    void   LoopESD(AliESDEvent *);
private:
    // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
    // your data member object is created on the worker nodes and streaming is not needed.
    // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
    TList  *fListHist;  //! List of histograms
    AliESDtrackCuts* fCuts; //!
    AliAnalysisFilter* fTrackFilterSKE; //!    //Added by me

//===========================================================================================
//   Histograms
//===========================================================================================
    
    TH1D *fHistEventCounter; //!
    TH1D *fHistRefMult08; //!
    TH1D *fHistRefMult05; //!
    TH1D *fHistV0M; //!
    TH1D *fHistV0A; //!
    TH1D *fHistV0C; //!
    TH1D *fHistV0Aamp; //!
    TH1D *fHistV0Camp; //!
    TH1D *fHistV0Mamp; //!
    TH1D *fdNdeta; //!
    TH1D *fPMult08; //!
    TH1D *fdNdeta05; //!
    TH1D *fPMult05; //!
    TH2D *fcorrRef05Ref08; //!
    TH2D *fcorrV0ARef08; //!
    TH2D *fcorrV0CRef08; //!
    TH2D *fcorrV0MRef08; //!
    TH2D *fcorrV0AampRef08; //!
    TH2D *fcorrV0CampRef08; //!
    TH2D *fcorrV0MampRef08; //!
    TProfile *fcorrRef05Ref08pfx; //!
    TProfile *fcorrV0ARef08pfx; //!
    TProfile *fcorrV0CRef08pfx; //!
    TProfile *fcorrV0MRef08pfx; //!
    TProfile *fcorrV0AampRef08pfx; //!
    TProfile *fcorrV0CampRef08pfx; //!
    TProfile *fcorrV0MampRef08pfx; //!
    TH2D *fModulesV0; //!
    
    AliPPVsMultUtils *fPPVsMultUtils;

    AliAnalysisMultCorrTaskQA(const AliAnalysisMultCorrTaskQA&);            // not implemented
    AliAnalysisMultCorrTaskQA& operator=(const AliAnalysisMultCorrTaskQA&); // not implemented

    ClassDef(AliAnalysisMultCorrTaskQA, 11); //11
};

#endif
