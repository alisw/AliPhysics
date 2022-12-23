/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskDensity_H
#define AliAnalysisTaskDensity_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliAODMCParticle.h"


// root basics
#include "TString.h"
#include "TRandom.h"
#include "TList.h"
#include "TH1D.h"

// ali containers
#include "AliPtSubEventContainer.h"

class TList;
class AliAODEvent;
class AliAODTrack;
class AliMCEvent;

class AliPtSubEventContainer;


class AliAnalysisTaskDensity : public AliAnalysisTaskSE  
{
    public:
        AliAnalysisTaskDensity();
        AliAnalysisTaskDensity(const char *name);
        virtual ~AliAnalysisTaskDensity();

        virtual void UserCreateOutputObjects();
        virtual void UserExec(Option_t* option);
        virtual void Terminate(Option_t* option);

        void ProcessMCParticles();

        void SetDefaultSettings();
        void SetMode(TString mode){fMode = mode;};
        void SetCorrelationOrder(Int_t m){ fMPar = m;};
        void SetMultSelectionMethod(TString method){ fMultSelMethod = method;};

        void SetMaxPileup(Int_t npilup){ fMaxPileup = npilup;};
        void SetMaxPrimaryVz(Double_t vz){ fMaxPrimaryVertexZ = vz; };
        void SetTightEventCuts(Double_t zvertex, Int_t pileup);

        void SetFilterBit(UInt_t bit){ fFilterBit = bit; };
        void SetTPCMinCls(UShort_t tpcmin){fTPCMinCls = tpcmin;};
        void SetPtRange(Double_t low, Double_t high){ fPtLow = low; fPtHigh = high;};
        void SetMaxDCAz(Double_t dcaz){fMaxDCAz = dcaz;};
        void SetTightTrackCuts(UInt_t filterbit, UShort_t ncls, Double_t dcaz);

        bool AcceptAODTrack(AliAODTrack* track);
        bool AcceptMCTrack(AliAODMCParticle* track);
        bool AcceptAODEvent(AliAODEvent* event);

        void FillTrackSelection(AliAODTrack* track);

        void ProcessEventCorrelation(Int_t id);
        void FillWPCounter(vector<vector<double>> &inarr, double w, double p);
        void ClearWPCounter();

    protected:
        // QA
        AliEventCuts fEventCuts;

    private:
        // containers
        TList*                  fPtSampleList;          //! output list for BS samples        
        TList*                  fQAList;                //! QA AliEvents cuts 
        TList*                  fQATrackList;           //! QA for track cuts
        // classes
        AliAODEvent*            fAOD;                    //! input event
        AliMCEvent*             fMCEvent;                //! mc event
        AliPtSubEventContainer** PtSubContainer;        
        // event settings
        Double_t    fCentrality;
        Double_t    fMaxPrimaryVertexZ;
        Int_t       fMaxPileup;
        Double_t    fPtLow;
        Double_t    fPtHigh;
        // track settings
        UInt_t      fFilterBit;
        UShort_t    fTPCMinCls;
        Int_t       fMPar;
        Int_t       fNOfTPCClusters;
        Double_t    fMaxDCAz;

        TString     fMultSelMethod;
        TString     fMode;

        // histos
        TH1D*       fhDCAzDistribution;
        TH1D*       fhDCAxyDistribution;
        TH1D*       fhEtaDistribution;
        TH1D*       fhPtDistribution;

        // Other
        std::vector<std::vector<std::vector<Double_t>>> wpSubEvent;
        std::vector<std::vector<Double_t>> wp;
        TRandom* rndGenerator;


        AliAnalysisTaskDensity(const AliAnalysisTaskDensity&); // not implemented
        AliAnalysisTaskDensity& operator=(const AliAnalysisTaskDensity&); // not implemented

        ClassDef(AliAnalysisTaskDensity, 1);
};

#endif
