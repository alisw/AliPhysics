#ifndef AliAnalysisBGMonitorQAHLT_H
#define AliAnalysisBGMonitorQAHLT_H
#include "AliAnalysisTaskSE.h"
class AliESDEvent;
class AliESDfriend;
class AliAnalysisCuts;
class TH1D;
class TH1F;
class TH2F;
class TH2D;
class AliAnalysisBGMonitorQAHLT : public AliAnalysisTaskSE {
    public:
    AliAnalysisBGMonitorQAHLT();
    AliAnalysisBGMonitorQAHLT(const char* name);
    virtual                            ~AliAnalysisBGMonitorQAHLT();
    virtual void                    UserCreateOutputObjects();
    virtual void                    Exec(Option_t* option);
    virtual void                    Terminate(Option_t* option);
    virtual void                    DrawHist(Int_t* ftrigger, Int_t fSpdT, Int_t fSpdC1, Int_t fSpdC2, Int_t* BBFlagC, Int_t* BBFlagA);
    // virtual void   Terminate(Option_t *);
    private:
    AliESDEvent* fESD;        //! ESD event
    AliESDfriend* fESDfriend; //! ESDfriend
    TTree *fTreeTrack2;        //! tree
    TList *fList;             //! list
    Int_t fUseTree;
    Int_t runNumber;
    Int_t ftrigger[kMaxUShort];
    Int_t fSpdC1;
    Int_t fSpdC2;
    Int_t fSpdT;
    Int_t ntracks;
    Int_t nV0A;
    Int_t nV0C;
    Int_t nV0ABG;
    Int_t nV0CBG;
    Int_t bgID;
    Int_t BGFlagA[kMaxUShort];
    Int_t BGFlagC[kMaxUShort];
    Int_t BBFlagA[kMaxUShort];
    Int_t BBFlagC[kMaxUShort];
    Int_t ADBGFlagA[kMaxUShort];
    Int_t ADBGFlagC[kMaxUShort];
    Int_t ADBBFlagA[kMaxUShort];
    Int_t ADBBFlagC[kMaxUShort];
    UShort_t ntr;
    UShort_t nbunch;
    AliAnalysisBGMonitorQAHLT(const AliAnalysisBGMonitorQAHLT&); // not implemented
    AliAnalysisBGMonitorQAHLT &operator=(const AliAnalysisBGMonitorQAHLT&); // not implemented
    ClassDef(AliAnalysisBGMonitorQAHLT,2);// example of analysis
};
#endif