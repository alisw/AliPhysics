#ifndef ALIJCDIJETTASK_H
#define ALIJCDIJETTASK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Analysis task for providing various dijet informations 
// author: O. Saarimaki, D.J. Kim(dong.jo.kim@cern.ch)
// ALICE Group University of Jyvaskyla 
// Finland 
//
// Fill the analysis containers for ESD or AOD
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <TDirectory.h>
#include <TComplex.h>
#include <AliLog.h>
#include <AliAnalysisTaskSE.h>
#include <AliJCatalystTask.h>
#include "AliJCDijetHistos.h"
#include "AliJCDijetAna.h"
#include "AliMCEvent.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliEventCuts.h"
#include "AliYAMLConfiguration.h"



using namespace std;
//==============================================================
class TClonesArray;
class AliJCDijetHistos;
class AliJCDijetAna;

class AliJCDijetTask : public AliAnalysisTaskSE {

    public:
        AliJCDijetTask();
        AliJCDijetTask(const char *name,  TString inputformat);
        AliJCDijetTask(const AliJCDijetTask& ap);   
        AliJCDijetTask& operator = (const AliJCDijetTask& ap);
        virtual ~AliJCDijetTask();

        static AliAnalysisTask *AddTaskJCDijetTask(TString taskName,
                                    Bool_t isMC,
                                    TString sJCatalyst        = "JCatalystTask",
                                    TString sJCatalystDetMC   = "JCatalystDetMCTask",
                                    UInt_t flags              = 0,
                                    TString centBins          = "0.0 5.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0",
                                    TString sDijetMBins       = "0, 20, 40, 45, 55, 65, 75, 85, 100, 120, 150, 250, 400, 500, 100000",
                                    double jetCone            = 0.4,
                                    double ktjetCone          = 0.4,
                                    int ktScheme              = 1,
                                    int antiktScheme          = 1,
                                    Bool_t usePionMass        = false,
                                    Bool_t useDeltaPhiBGSubtr = true,
                                    double particleEtaCut     = 0.8,
                                    double particlePtCut      = 0.15,
                                    double leadingJetCut      = 20.0,
                                    double subleadingJetCut   = 20.0,
                                    double minJetPt           = 10.0,
                                    double constituentCut     = 5.0,
                                    double deltaPhiCut        = 2.0,
                                    double matchingR          = 0.2,
                                    double trackingIneff      = 0.0,
                                    TString sAnchorPeriodForTracking = "",
                                    AliJCDijetAna::jetClasses lUnfJetClassTrue = AliJCDijetAna::iAcc,
                                    AliJCDijetAna::jetClasses lUnfJetClassDet = AliJCDijetAna::iAcc,
                                    Bool_t useCoveredAreaRho  = false);


        // methods to fill from AliAnalysisTaskSE
        virtual void UserCreateOutputObjects(); 
        virtual void Init();   
        virtual void LocalInit() { Init(); }
        virtual void UserExec(Option_t *option);
        virtual void Terminate(Option_t* );
        AliJCatalystTask *GetJCatalystTask() {return fJCatalystTask;}
        void    SetCentralityBins( vector<double> centralityBins ) {fcentralityBins=centralityBins; }
        void    SetDijetMBins( TString dijetMBins ) {fsDijetMBins=dijetMBins; }
        void    SetJetConeSize(double jetCone, double ktjetCone) {fjetCone=jetCone; fktJetCone=ktjetCone; }
        void    SetBGSubtrSettings(int ktScheme, int antiktScheme, Bool_t usePionMass, Bool_t useDeltaPhiBGSubtr, Bool_t luseCrho) {fktScheme=ktScheme; fantiktScheme=antiktScheme; fusePionMass=usePionMass; fuseDeltaPhiBGSubtr=useDeltaPhiBGSubtr; bUseCrho=luseCrho;}
        void    SetUnfoldingJetSets(AliJCDijetAna::jetClasses lJetClassTrue, AliJCDijetAna::jetClasses lJetClassDet) { iUnfJetClassTrue = lJetClassTrue; iUnfJetClassDet = lJetClassDet;}
        Bool_t  IsMC()const{ return fIsMC; }
        void    SetIsMC(Bool_t b) { fIsMC=b; }
        void    SetCuts(double particleEta,
                        double particlePt,
                        double leadingJet,
                        double subleadingJet,
                        double constituent,
                        double deltaPhi,
                        double matchingR,
                        double trackingIneff,
                        double minJetPt) {
                        fparticleEtaCut=particleEta;
                        fparticlePtCut=particlePt;
                        fleadingJetCut=leadingJet;
                        fsubleadingJetCut=subleadingJet;
                        fconstituentCut=constituent;
                        fdeltaPhiCut=deltaPhi;
                        fmatchingR = matchingR;
                        ftrackingIneff = trackingIneff;
                        fMinJetPt = minJetPt;
        }
        void AddFlags(UInt_t nflags){flags |= nflags;}
        enum{
            DIJET_VERTEX13PA      = 0x1,
            DIJET_PILEUPSPD       = 0x2,
            DIJET_UTILSPILEUPSPD  = 0x4,
            DIJET_ALIEVENTCUT     = 0x8,
            DIJET_CATALYST        = 0x100
        };

        // Methods specific for this class
        void SetJCatalystTaskName(TString name){ fJCatalystTaskName=name; } // Setter for filter task name
        void SetJCatalystTaskNameDetMC(TString name){ fJCatalystDetMCTaskName=name; } // Setter for filter task name
        void AddArtificialTrackingEfficiencyConfig(TString sAnchorPeriod);
        void SetArtificialTrackingEfficiencyFromYAML();

        AliEventCuts fEventCuts; // Event cut object

    private:

        AliJCatalystTask *fJCatalystTask;  //
        AliJCatalystTask *fJCatalystDetMCTask;  //
        TString           fJCatalystTaskName; // Name for JCatalyst task
        TString           fJCatalystDetMCTaskName; // Name for JCatalyst task
        vector<double> fcentralityBins;
        TString fsDijetMBins;
        double fjetCone;
        double fktJetCone;
        int  fktScheme;
        int  fantiktScheme;
        Bool_t fusePionMass;
        Bool_t fuseDeltaPhiBGSubtr;
        Bool_t fIsMC;       // MC data or real data
        double fparticleEtaCut;
        double fparticlePtCut;
        double fleadingJetCut;
        double fsubleadingJetCut;
        double fMinJetPt;
        double fconstituentCut;
        double fdeltaPhiCut;
        double fmatchingR;
        double ftrackingIneff;
        AliJCDijetHistos *fhistos;
        AliJCDijetHistos *fhistosDetMC;
        AliJCDijetAna *fana;
        AliJCDijetAna *fanaMC;
        int fCBin;
        int fCBinDetMC;
        TDirectory     *fOutput; // Output directory
        UInt_t flags; //
        AliAnalysisUtils *fUtils; //!
        double fptHardBin;
        double fPythiaSigma;
        double fPythiaTrial;
        int fDetMCFlag;
        AliJCDijetAna::jetClasses iUnfJetClassTrue;
        AliJCDijetAna::jetClasses iUnfJetClassDet;
        bool bUseCrho;
        bool bGoodEvent;
        bool bGoodMCEvent;
        TH1D *fTrackEfficiencyHistogram; ///< This will define the tracking efficiency uncertainty in pt function. Implemented similarly as AliEmcalJetTask. At the moment only implemented for MB runs of pp and p--Pb.
        PWG::Tools::AliYAMLConfiguration fYAMLConfig; ///< yaml configuration
        TString fsAnchorPeriod;

        ClassDef(AliJCDijetTask, 3); 
};
#endif // ALIJCDIJETTASK_H
