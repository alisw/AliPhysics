/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskJetFlow_H
#define AliAnalysisTaskJetFlow_H

#include <AliAnalysisTaskSE.h>
#include <AliFlowTrackCuts.h>
#include <AliFlowEvent.h>
#include <TString.h>

class AliAnalysisTaskJetFlow : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskJetFlow();
                                AliAnalysisTaskJetFlow(const char *name);
        virtual                 ~AliAnalysisTaskJetFlow();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
        // setters
        void                    SetCutsRP(AliFlowTrackCuts* c)          {fCutsRP = c;}
        void                    SetCutsPOI(AliFlowTrackCuts* c)         {fCutsPOI = c;}
        void                    SetCutsNull(AliFlowTrackCuts* c)        {fCutsNull = c;}
        void                    SetJetCollectionName(TString jets)      {fJetsName = jets;}
        void                    SetDebugMode(Int_t d)                   {fDebug = d;}

    private:

        // technical stuff
        Int_t                   fDebug;         // debug level (0 none, 1 fcn calls, 2 verbose)
        TString                 fJetsName;      // name of jet list
        TList*                  fOutputList;    //! output list
        // cut objects
        AliFlowTrackCuts*       fCutsRP;        // rp cuts
        AliFlowTrackCuts*       fCutsPOI;       // poi cuts
        AliFlowTrackCuts*       fCutsNull;      // empty cuts
        // input, output
        AliFlowEvent*           fFlowEvent;     //! flow event
        // histograms
        /* TH1F*                   fHistPt;        //! dummy histogram */

        AliAnalysisTaskJetFlow(const AliAnalysisTaskJetFlow&); // not implemented
        AliAnalysisTaskJetFlow& operator=(const AliAnalysisTaskJetFlow&); // not implemented

        ClassDef(AliAnalysisTaskJetFlow, 1);
};

#endif
