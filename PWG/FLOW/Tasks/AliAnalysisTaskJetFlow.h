/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskJetFlow_H
#define AliAnalysisTaskJetFlow_H

#include <AliAnalysisTaskSE.h>
#include <AliFlowTrackCuts.h>
#include <AliFlowEvent.h>
#include <TString.h>
#include <AliVEvent.h>

class AliAnalysisTaskJetFlow : public AliAnalysisTaskSE
{
    public:
        // enumerators
        enum dataType           {kESD, kAOD, kESDMC, kAODMC };  // data type
        // constructors, destructor
                                AliAnalysisTaskJetFlow();
                                AliAnalysisTaskJetFlow(const char *name);
        virtual                 ~AliAnalysisTaskJetFlow();
        // virtual methods
        virtual void            LocalInit();
        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
        // setters
        void                    SetCutsRP(AliFlowTrackCuts* c)          {fCutsRP        = c;}
        void                    SetCutsPOI(AliFlowTrackCuts* c)         {fCutsPOI       = c;}
        void                    SetCutsNull(AliFlowTrackCuts* c)        {fCutsNull      = c;}
        void                    SetJetCollectionName(TString jets)      {fJetsName      = jets;}
        void                    SetDebugMode(Int_t d)                   {fDebug         = d;}
        void                    SetPtBump(Float_t b)                    {fPtBump        = b;}
        void                    SetMinMaxCentrality(Float_t min, Float_t max)       {fCentralityMin = min; fCentralityMax = max; }
        // analysis details
        Bool_t                  PassesCuts(AliVEvent* event); 

    private:

        // analysis flags and task setup specifics
        Int_t                   fDebug;                 // debug level (0 none, 1 fcn calls, 2 verbose)
        TString                 fJetsName;              // name of jet list
        TList*                  fOutputList;            //! output list
        dataType                fDataType;              //! data type
        // members
        Float_t                 fPtBump;                // track pt += ptbump
        Float_t                 fCentralityMin;         // minimium centrality
        Float_t                 fCentralityMax;         // maximum centrality
        // cut objects
        AliFlowTrackCuts*       fCutsRP;                // rp cuts
        AliFlowTrackCuts*       fCutsPOI;               // poi cuts
        AliFlowTrackCuts*       fCutsNull;              // empty cuts
        // containers, setup
        AliFlowEvent*           fFlowEvent;             //! container for flow analysis
        // histograms
        TH1F*                   fHistAnalysisSummary;   //! analysis summary

        AliAnalysisTaskJetFlow(const AliAnalysisTaskJetFlow&);                  // not implemented
        AliAnalysisTaskJetFlow& operator=(const AliAnalysisTaskJetFlow&);       // not implemented

        ClassDef(AliAnalysisTaskJetFlow, 2);
};

#endif
