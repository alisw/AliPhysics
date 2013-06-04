/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskJetFlow_H
#define AliAnalysisTaskJetFlow_H

// root includes
#include <TMath.h>
//aliroot includes
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
                                AliAnalysisTaskJetFlow(
                                        const char *name,
                                        AliFlowTrackCuts* rpCutsTPC,
                                        AliFlowTrackCuts* rpCutsVZERO,
                                        TString jetName,
                                        TString picoName);
        virtual                 ~AliAnalysisTaskJetFlow();
        // virtual methods
        virtual void            LocalInit();
        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
        // setters
        void                    SetExplicitOutlierCut(Int_t c)          {fExplicitOutlierCut = c;}
        void                    SetDebugMode(Int_t d)                   {fDebug         = d;}
        void                    SetPtBump(Float_t b)                    {fPtBump        = b;}
        void                    SetCCMaxPt(Float_t m)                   {fCCMaxPt       = m;}
        void                    SetCCBinsInPt(Int_t b)                  {fCCBinsInPt    = b;}
        void                    SetMinMaxCentrality(Float_t min, Float_t max)   {fCentralityMin = min; fCentralityMax = max; }
        void                    SetMinMaxPOIPt(Float_t min, Float_t max)        {fPOIPtMin = min; fPOIPtMax = max; }        
        void                    SetDoVParticleAnalysis(Bool_t d)        {fVParticleAnalysis = d; }
        void                    SetMinimizeDiffBins(Bool_t b)           {fMinimizeDiffBins = b; }
        void                    SetDoTestFlowAnalysis(Bool_t t, TArrayD* pt = 0x0)
                {fDoTestFlowAnalysis = t; 
                 fPtBins = pt;  }
        void                    SetDoMultWeight(Bool_t m)               {fDoMultWeight = m; }
        // analysis details
        Bool_t                  PassesCuts(AliVEvent* event);
        Bool_t                  PassesCuts(Int_t year); 
        void                    DoTestFlowAnalysis();
        /* inline */    Double_t PhaseShift(Double_t x) const {  
            while (x>=TMath::TwoPi())x-=TMath::TwoPi();
            while (x<0.)x+=TMath::TwoPi();
            return x; }
        /* inline */    Double_t PhaseShift(Double_t x, Double_t n) const {
            x = PhaseShift(x);
            if(TMath::Nint(n)==2) while (x>TMath::Pi()) x-=TMath::Pi();
            if(TMath::Nint(n)==3) {
                if(x>2.*TMath::TwoPi()/n) x = TMath::TwoPi() - x;
                if(x>TMath::TwoPi()/n) x = TMath::TwoPi()-(x+TMath::TwoPi()/n);
            }
            return x; }
    private:

        // analysis flags and task setup specifics
        Int_t                   fDebug;                 // debug level (0 none, 1 fcn calls, 2 verbose)
        Int_t                   fExplicitOutlierCut;    // cut on multiplicity ourliers explicitely (slow)
        TString                 fJetsName;              // name of jet list
        TString                 fTracksName;            // name of track list
        TList*                  fOutputList;            //! output list
        dataType                fDataType;              //! data type
        Bool_t                  fVParticleAnalysis;     // do the analysis on vparticles instead of jets
        Bool_t                  fMinimizeDiffBins;      // minimize variables (for low statistics)
        Bool_t                  fDoTestFlowAnalysis;    // do a quick and dirty crude flow estimate
        Bool_t                  fDoMultWeight;          // weight events with multiplicity
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
        AliFlowTrackCuts*       fCutsNull;              // empty cuts
        AliFlowEventCuts*       fCutsEvent;             // event cuts
        // containers, setup
        AliFlowEvent*           fFlowEvent_TPC;         //! container for flow analysis
        AliFlowEvent*           fFlowEvent_VZERO;       //! container for flow analysis
        // histograms
        TH1F*                   fHistAnalysisSummary;   //! analysis summary
        TH1F*                   fCentralitySelection;   //! centrality selection
        TH1F*                   fVZEROAEP;              //! VZEROA EP
        TH1F*                   fVZEROCEP;              //! VZEROC EP
        TProfile*               fv2VZEROA;              //! v2 from VZEROA
        TProfile*               fv2VZEROC;              //! v2 from VZEROC
        TProfile*               fTempA;                 //! internal bookkeeping
        TProfile*               fTempC;                 //! internal bookkeeping

        AliAnalysisTaskJetFlow(const AliAnalysisTaskJetFlow&);                  // not implemented
        AliAnalysisTaskJetFlow& operator=(const AliAnalysisTaskJetFlow&);       // not implemented

        ClassDef(AliAnalysisTaskJetFlow, 4);
};

#endif
