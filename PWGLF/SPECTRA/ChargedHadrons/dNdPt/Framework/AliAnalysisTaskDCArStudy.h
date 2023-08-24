#ifndef AliAnalysisTaskDCArStudy_H
#define AliAnalysisTaskDCArStudy_H

#include "AliAnalysisTaskMKBase.h"
#include "AliAnalysisHelpersHist.h"
#include "THn.h"

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
class AliMCSpectraWeights;

class AliAnalysisTaskDCArStudy : public AliAnalysisTaskMKBase
{
    public:
                                AliAnalysisTaskDCArStudy();
                                AliAnalysisTaskDCArStudy(const char *name);
        virtual                 ~AliAnalysisTaskDCArStudy();

        virtual void            AddOutput();                     //called at the beginning
        virtual Bool_t          IsEventSelected();               //called for each event
        virtual void            AnaEventDATA();                  //called once for every selected DATA event
        virtual void            AnaEventMC();                    //called once for every selected mc event
        virtual void            AnaTrackDATA(Int_t flag = 0);    //called once for every track in DATA track
        virtual void            AnaTrackMC(Int_t flag = 0);      //called once for every track in MC track
        
        static AliAnalysisTaskDCArStudy* AddTaskDCArStudy(const char* name = "TaskDCArStudy", const char* outfile = 0);

    protected:
        Hist<THnF>        fHistDCA; //!<! hist for DCA distros
        Hist<THnF>        fHistDCAPCC; //!<! hist for DCA distros with PCC scaling
        Hist<THnF>        fHistDCAPCCSysUp; //!<! hist for DCA distros with PCC scaling                                    //!
        Hist<THnF>        fHistDCAPCCSysDown; //!<! hist for DCA distros with PCC scaling                                    //!
        Hist<THnF>        fHistSecWeights;//!
        Hist<THnF>        fEventHist;//!
        AliMCSpectraWeights*    fMCSpectraWeights; //!<! object to determine efficiency scaling
        
    private:
        AliAnalysisTaskDCArStudy(const AliAnalysisTaskDCArStudy&); // not implemented
        AliAnalysisTaskDCArStudy& operator=(const AliAnalysisTaskDCArStudy&); // not implemented

        ClassDef(AliAnalysisTaskDCArStudy, 1);
};

#endif
