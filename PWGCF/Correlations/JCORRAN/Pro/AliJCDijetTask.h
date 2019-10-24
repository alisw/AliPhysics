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

        // methods to fill from AliAnalysisTaskSE
        virtual void UserCreateOutputObjects(); 
        virtual void Init();   
        virtual void LocalInit() { Init(); }
        virtual void UserExec(Option_t *option);
        virtual void Terminate(Option_t* );
        AliJCatalystTask *GetJCatalystTask() {return fJCatalystTask;}
        void    SetCentralityBins( vector<double> centralityBins ) {fcentralityBins=centralityBins; }
        void    SetJetConeSize(double jetCone, double ktjetCone) {fjetCone=jetCone; fktJetCone=ktjetCone; }
        void    SetBGSubtrSettings(int ktScheme, int antiktScheme, Bool_t usePionMass, Bool_t useDeltaPhiBGSubtr) {fktScheme=ktScheme; fantiktScheme=antiktScheme; fusePionMass=usePionMass; fuseDeltaPhiBGSubtr=useDeltaPhiBGSubtr; }
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
            DIJET_UTILSPILEUPSPD  = 0x4
        };

        // Methods specific for this class
        void SetJCatalystTaskName(TString name){ fJCatalystTaskName=name; } // Setter for filter task name
        void SetJCatalystTaskNameDetMC(TString name){ fJCatalystDetMCTaskName=name; } // Setter for filter task name

    private:

        AliJCatalystTask *fJCatalystTask;  //
        AliJCatalystTask *fJCatalystDetMCTask;  //
        TString           fJCatalystTaskName; // Name for JCatalyst task
        TString           fJCatalystDetMCTaskName; // Name for JCatalyst task
        vector<double> fcentralityBins;
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

        ClassDef(AliJCDijetTask, 1); 
};
#endif // ALIJCDIJETTASK_H
