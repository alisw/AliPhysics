/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskJetFlow_H
#define AliAnalysisTaskJetFlow_H

// includes
#include <AliAnalysisTaskSE.h>
// forward declarations
class TString;
class TList;
class TArrayD;
class AliFlowTrackCuts;
class AliFlowEventCuts;
class AliFlowEvent;
class TH1;

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
        void                    SetCutsRP(AliFlowTrackCuts* tpc, AliFlowTrackCuts* vzero)  {fCutsRP_TPC = tpc; fCutsRP_VZERO = vzero; }

        void                    SetCutsPOI(AliFlowTrackCuts* c)         {fCutsPOI       = c;}
        void                    SetCutsNull(AliFlowTrackCuts* c)        {fCutsNull      = c;}
        void                    SetJetCollectionName(TString jets)      {fJetsName      = jets;}
        void                    SetDebugMode(Int_t d)                   {fDebug         = d;}
        void                    SetPtBump(Float_t b)                    {fPtBump        = b;}
        void                    SetCCMaxPt(Float_t m)                   {fCCMaxPt       = m;}
        void                    SetCCBinsInPt(Int_t b)                  {fCCBinsInPt    = b;}
        void                    SetMinMaxCentrality(Float_t min, Float_t max)   {fCentralityMin = min; fCentralityMax = max; }
        void                    SetMinMaxPOIPt(Float_t min, Float_t max)        {fPOIPtMin = min; fPOIPtMax = max; }        
        void                    SetDoVParticleAnalysis(Bool_t d)        {fVParticleAnalysis = d; }
        void                    SetDoTestFlowAnalysis(Bool_t t, TArrayD* pt = 0x0)
                {fDoTestFlowAnalysis = t; 
                 fPtBins = pt;  }
        // analysis details
        Bool_t                  PassesCuts(AliVEvent* event); 
        void                    DoTestFlowAnalysis();

    private:

        // analysis flags and task setup specifics
        Int_t                   fDebug;                 // debug level (0 none, 1 fcn calls, 2 verbose)
        TString                 fJetsName;              // name of jet list
        TList*                  fOutputList;            //! output list
        dataType                fDataType;              //! data type
        Bool_t                  fVParticleAnalysis;     // do the analysis on vparticles instead of jets
        Bool_t                  fDoTestFlowAnalysis;    // do a quick and dirty crude flow estimate
        Bool_t                  fInitialized;           //! check if the analysis is initialized
        // members
        Float_t                 fPtBump;                // track pt += ptbump
        Float_t                 fCCMaxPt;               // max pt for flow analysis (common constants)
        Float_t                 fCCBinsInPt;            // bins in pt for flow analysis (common constants)
        Float_t                 fCentralityMin;         // minimium centrality
        Float_t                 fCentralityMax;         // maximum centrality
        Float_t                 fPOIPtMin;              // minimum pt for poi's
        Float_t                 fPOIPtMax;              // maximum pt for poi's
        TArrayD*                fPtBins;                // pt bins for flow analysis
        // cut objects
        AliFlowTrackCuts*       fCutsRP_TPC;            // rp cuts for tpc
        AliFlowTrackCuts*       fCutsRP_VZERO;          // rp cuts for fzero
        AliFlowTrackCuts*       fCutsPOI;               // poi cuts
        AliFlowTrackCuts*       fCutsNull;              // empty cuts
        AliFlowEventCuts*       fCutsEvent;             // event cuts
        // containers, setup
        AliFlowEvent*           fFlowEvent_TPC;         //! container for flow analysis
        AliFlowEvent*           fFlowEvent_VZERO;       //! container for flow analysis
        // histograms
        TH1F*                   fHistAnalysisSummary;   //! analysis summary
        TH1F*                   fCentralitySelection;   //! centrality selection
        TProfile*               fv2VZEROA;              //! v2 from VZEROA
        TProfile*               fv2VZEROC;              //! v2 from VZEROC

        AliAnalysisTaskJetFlow(const AliAnalysisTaskJetFlow&);                  // not implemented
        AliAnalysisTaskJetFlow& operator=(const AliAnalysisTaskJetFlow&);       // not implemented

        ClassDef(AliAnalysisTaskJetFlow, 3);
};

#endif
