#ifndef AliAnalysisTaskMultspec_FDcorr_H
#define AliAnalysisTaskMultspec_FDcorr_H

//class AliESDpid;
//class AliESDEvent;

#include "TString.h"
#include "AliPIDResponse.h"
#include "AliAnalysisTaskSE.h"
#include "THistManager.h"
#include "TH1D.h"
#include "TH2D.h"
#include "AliMCEvent.h"

struct LamFiller_FD{

    Double32_t Pt_lam;
    Double32_t Pt_xi;
    Double32_t DcaV0Daught;     //[0,2.54,8]
    Double32_t DcaPosToPV;      //[0,2.54,8]
    Double32_t DcaNegToPV;      //[0,2.54,8]
    Double32_t V0Rad;           //[0,25.4,8]
    Double32_t V0CosPA;         //[0.95,1,16]
    Double32_t LeastCRawsOvF;   //[0,2.54,8]
    Double32_t DistOverTotP;    //[0,254,8]
    Double32_t NSigPos;         //[-10,10,8]
    Double32_t NSigNeg;         //[-10,10,8]
    int    LeastCRaws;          //[0,254,8]
    bool   TOFmatch;
    bool   ITSmatch;
    bool   IsPrimary;
    bool   IsParticle;

};


class AliAnalysisTaskMultspec_FDcorr : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskMultspec_FDcorr();
    AliAnalysisTaskMultspec_FDcorr(const char *name, TString lExtraOptions = "");
    virtual ~AliAnalysisTaskMultspec_FDcorr();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);
    
  private:
    
    //outputs
    THistManager *fHistos_misc;                               //!<! Output histos
    TTree        *fTree;                                      //!<! Output Tree
    
    //fillers
    LamFiller_FD *ffillV0 = nullptr;                               //!<! Transient V0filler
    
    //objects retreived from input handler
    AliPIDResponse *fPIDResponse;                              //!
    
    //AliEventCuts object
    AliEventCuts fEventCuts;                                   //

    //functions to allow flushing part of code out of UserExec
    void DataPosting();
    
    AliAnalysisTaskMultspec_FDcorr(const AliAnalysisTaskMultspec_FDcorr&);            // not implemented
    AliAnalysisTaskMultspec_FDcorr& operator=(const AliAnalysisTaskMultspec_FDcorr&); // not implemented

    ClassDef(AliAnalysisTaskMultspec_FDcorr, 1);


};

#endif
