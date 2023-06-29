/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliSPDppAnalysisTaskSim_H
#define AliSPDppAnalysisTaskSim_H

#include "AliAnalysisTaskSE.h"

#include "AliAnalysisUtils.h"
#include "AliGenEventHeader.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "TParticle.h"
#include "AliEventCuts.h"

class AliSPDppAnalysisTaskSim : public AliAnalysisTaskSE  
{
    public:
                                AliSPDppAnalysisTaskSim();
                                AliSPDppAnalysisTaskSim(const char *name);
        virtual                 ~AliSPDppAnalysisTaskSim();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
        virtual void            ProcessMCParticles();
        virtual void            ProcessData();
    
        void                    UsekINT1(bool use_int1) { fUseINT1 = use_int1; }
    

    private:
        
        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list
        TH1F*                   MultDist05;     //! multiplicity distribution histogram for eta<.5>
        TH1F*                   RecMultDist05;     //! multiplicity distribution histogram for eta<.5>
        TH2F*                   responseMatrix;  //!
        TH1F*                   MultDist05inelgr0;     //! multiplicity distribution histogram for eta<.5>
        TH1F*                   RecMultDist05inelgr0;     //! multiplicity distribution histogram for eta<.5>
        TH2F*                   responseMatrixinelgr0;  //!
    
    
        AliVMultiplicity*       fMultiplicity=nullptr;
    
        AliMCEvent*             fMCEvent;       //! corresponding MC event
    
        AliEventCuts            fEventCuts;
        AliAODVZERO             *fAODV0;
    
        bool fUseINT1;
    
        TList                   *fQAList;       //!

        
        Int_t eventcount1 = 0;
    
        AliSPDppAnalysisTaskSim(const AliSPDppAnalysisTaskSim&); // not implemented
        AliSPDppAnalysisTaskSim& operator=(const AliSPDppAnalysisTaskSim&); // not implemented

        ClassDef(AliSPDppAnalysisTaskSim, 1);
};

#endif
