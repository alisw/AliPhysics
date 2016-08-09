#ifndef AliAnalysisBGMonitorQAHLT_H
#define AliAnalysisBGMonitorQAHLT_H
#include "AliAnalysisTask.h"
class AliVEvent;
class AliVfriendEvent;
class AliVMultiplicity;
class AliVVZERO;
class AliVVZEROfriend;
class TList;

class AliAnalysisBGMonitorQAHLT : public AliAnalysisTask {
    public:
    AliAnalysisBGMonitorQAHLT();
    AliAnalysisBGMonitorQAHLT(const char* name);
    virtual                            ~AliAnalysisBGMonitorQAHLT();
    virtual void                    CreateOutputObjects();
    virtual void                    Exec(Option_t* option);
    virtual void                    Terminate(Option_t* option);
    virtual Bool_t                  ResetOutputData();
    
    static Int_t ResetHistograms(TList* list);
    static void CreateHistograms(TList*& list, Option_t* options=NULL);
    static void FillHistograms(TList* list,
                               std::string* firedTriggerClasses,
                               AliVMultiplicity* mult,
                               AliVVZERO* vzero,
                               AliVVZEROfriend* vzeroFriend );

    private:
    AliVEvent* fESD;        //! ESD event
    AliVfriendEvent* fESDfriend; //! ESDfriend
    TList *fList;             //! list
    Int_t fUseTree;
    Int_t triggerType;
    Int_t fSpdClusters;
    Int_t fSpdTracklets;
    Int_t ntracks;
    Int_t nV0A;
    Int_t nV0C;
    Int_t nV0ABG;
    Int_t nV0CBG;
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
