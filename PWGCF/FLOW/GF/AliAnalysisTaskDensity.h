/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskDensity_H
#define AliAnalysisTaskDensity_H

// ali basic
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliAODMCParticle.h"
#include "AliGFWCuts.h"

#include <iostream>
#include <utility>
#include <vector>

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
class AliGFWCuts;
class TProfile;
class TH2D;

class AliPtSubEventContainer;

class AliAnalysisTaskDensity : public AliAnalysisTaskSE  
{
    public:
        AliAnalysisTaskDensity();
        AliAnalysisTaskDensity(const char *name, Bool_t bUseEff);
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

        void SetFilterBit(UInt_t bit){ fFilterBit = bit; };
        void SetTPCMinCls(UShort_t tpcmin){fTPCMinCls = tpcmin;};
        void SetPtRange(Double_t low, Double_t high){ fPtLow = low; fPtHigh = high;};
        void SetMaxDCAz(Double_t dcaz){fMaxDCAz = dcaz;};
        void SetEtaGap(Double_t gap) { fEtaGap = gap; };
        void SetSystFlag(Int_t iflag){ if(!fGFWSelection) fGFWSelection = new AliGFWCuts(); fGFWSelection->SetupCuts(iflag);};
        void SetDCAxyFunctionalForm(TString iform) { fDCAxyFunctionalForm = iform; }; //Call after SystFlag
        void SetAbsEta(Double_t ieta) {fAbsEta = ieta;};

        bool CheckTrigger(Double_t lCent);
        bool AcceptAODTrack(AliAODTrack* track);
        bool AcceptMCTrack(AliAODMCParticle* track);
        bool AcceptAODEvent(AliAODEvent* event);

        void ProcessTrack(Double_t lweight, Double_t lpt, Double_t leta);
        void ProcessEventCorrelation(Int_t id);
        void FillWPCounter(std::vector<std::vector<Double_t>> &inarr, Double_t w, Double_t p);
        void FillEventQA(Double_t lCent, Int_t lNch);
        void ClearWPCounter();

        Double_t WeightPt(Double_t pt);

    protected:
    // QA
        AliEventCuts fEventCuts;

    private:
    // input/output containers
        TList*                  fPtSampleList;          //! output list for BS samples        
        TList*                  fInputListEfficiency;   
        TList*                  fQAEventList;           //! QA AliEvents cuts + Costum QA 
        TList*                  fQATrackList;           //! QA for track cuts
    // ali classes
        AliGFWCuts*                 fGFWSelection;
        AliAODEvent*                fAOD;                  //! input event
        AliMCEvent*                 fMCEvent;              //! mc event
        AliPtSubEventContainer**    PtSubContainer;         
    // task setting
        Int_t       fSystFlag;

    // event settings
        Bool_t      fUseEffeciency;
        Double_t    fCentrality;
        Double_t    fMaxPrimaryVertexZ;
        Int_t       fMaxPileup;
        Double_t    fPtLow;
        Double_t    fPtHigh;
        Double_t    fEtaGap;
        Double_t    fAbsEta;
    // track settings
        UInt_t      fFilterBit;
        UShort_t    fTPCMinCls;
        Int_t       fMPar;
        Int_t       fNOfTPCClusters;
        Double_t    fMaxDCAz;
        TString     fMultSelMethod;
        TString     fMode;
        TString     fDCAxyFunctionalForm;
        Double_t    fChi2PerClusterTPC;
    // histos
        TH1D*       fhDCAzDistribution;     //!
        TH1D*       fhDCAxyDistribution;    //!
        TH1D*       fhEtaDistribution;      //!
        TH2D*       fhPtDistribution;       //!
        TH2D*       fhCentSelected;         //!
        TH1D*       fhPrimaryVzSelected;    //!
        TH1D*       fhCrossedRowsTPC;       //!
        TH1D*       fhChiPerTPCCls;         //!
        TH1D*       fhChiPerITSCls;         //!
        TH1D*       fhNofPileupSelected;    //!
        TH2D*       fhSubEtaDistribution;   //!
        TH2D*       fhMultSelected;         //!

        TH1D*       fhEffeciencyPt;         //!

    // Other
        std::vector<std::vector<std::vector<Double_t>>> wpThreeSubEvent;
        std::vector<std::vector<std::vector<Double_t>>> wpTwoSubEvent;
        std::vector<std::vector<Double_t>> wp;
        TRandom* rndGenerator;
        
        ClassDef(AliAnalysisTaskDensity, 2);
};

#endif
