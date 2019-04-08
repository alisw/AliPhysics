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

        virtual void            AddOutput(); //called at the beginning
        virtual void            AnaTrack();  //called for every track
        virtual void            AnaEvent();  //called for every track        
        
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