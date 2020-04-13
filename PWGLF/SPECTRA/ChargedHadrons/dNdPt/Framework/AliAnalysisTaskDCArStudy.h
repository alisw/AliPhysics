#ifndef AliAnalysisTaskDCArStudy_H
#define AliAnalysisTaskDCArStudy_H

#include "AliAnalysisTaskMKBase.h"

class AliESDtrackCuts;
class AliVEvent;
class AliESDEvent;
class AliAODEvent;
class AliMCEvent;
class AliStack;
class AliHeader;
class AliGenEventHeader;
class AliESDtrack;
class AliMCParticle;

class AliAnalysisTaskDCArStudy : public AliAnalysisTaskMKBase
{
    public:
                                AliAnalysisTaskDCArStudy();
                                AliAnalysisTaskDCArStudy(const char *name);
        virtual                 ~AliAnalysisTaskDCArStudy();

        virtual void            AddOutput();                     //called at the beginning
        virtual Bool_t          IsEventSelected();               //called for each event
        virtual void            AnaEvent();                      //called once for every selected event        
        virtual void            AnaTrackDATA(Int_t flag = 0);    //called once for every track in DATA event
        virtual void            AnaTrackMC(Int_t flag = 0);      //called once for every track in MC event
        
        static AliAnalysisTaskDCArStudy* AddTaskDCArStudy(const char* name = "TaskDCArStudy", const char* outfile = 0);

    protected:     
        THnSparseD*             fHistDCA;    //-> dca hist  
        THnSparseD*             fHistDCATPC; //-> dca hist with tpc tracks              
        
    private:
        AliAnalysisTaskDCArStudy(const AliAnalysisTaskDCArStudy&); // not implemented
        AliAnalysisTaskDCArStudy& operator=(const AliAnalysisTaskDCArStudy&); // not implemented

        ClassDef(AliAnalysisTaskDCArStudy, 1);
};

#endif