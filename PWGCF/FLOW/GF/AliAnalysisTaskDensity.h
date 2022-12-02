/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskDensity_H
#define AliAnalysisTaskDensity_H

#include "AliAnalysisTaskSE.h"

// root basics
#include "TString.h"
#include "TRandom.h"
#include "TList.h"

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

        void SetDefaultCut();
        void SetMode(TString mode){fMode = mode;};
        void SetTPCMinCls(UShort_t tpcmin){fTPCMinCls = tpcmin;};
        void SetMaxPileup(Int_t npilup){ fMaxPileup = npilup;};
        void SetCorrelationOrder(Int_t m){ fMPar = m;};
        void SetPtRange(Double_t low, Double_t high){ fPtLow = low; fPtHigh = high;};
        void SetMultSelectionMethod(TString method){ fMultSelMethod = method;};
        void SetFilterBit(UInt_t bit){ fFilterBit = bit; };
        void SetPrimaryVertexZ(Double_t vz){ fLimitVertexZ = vz; };
 
        bool AcceptAODTrack(AliAODTrack* track);
             
        void FillWPCounter(vector<vector<double>> &inarr, double w, double p);
        void ClearWPCounter();

    private:
        // containers
        TList*                  fPtSampleList;          //! output list for samples

        // classes
        AliAODEvent*            fAOD;                    //! input event
        AliMCEvent*             fMCEvent;                //! mc event
        AliPtSubEventContainer** PtSubContainer;        

        // events and track variables
        Double_t    fCentrality;
        Double_t    fLimitVertexZ;
        Double_t    fPtLow;
        Double_t    fPtHigh;

        UInt_t      fFilterBit;
        UShort_t    fTPCMinCls;
        Int_t       fMPar;
        Int_t       fNOfTPCClusters;
        Int_t       fMaxPileup;

        TString     fMultSelMethod;
        TString     fMode;


        // Other
        std::vector<std::vector<std::vector<Double_t>>> wpSubEvent;
        std::vector<std::vector<Double_t>> wp;
        TRandom* rndGenerator;


        AliAnalysisTaskDensity(const AliAnalysisTaskDensity&); // not implemented
        AliAnalysisTaskDensity& operator=(const AliAnalysisTaskDensity&); // not implemented

        ClassDef(AliAnalysisTaskDensity, 1);
};

#endif
