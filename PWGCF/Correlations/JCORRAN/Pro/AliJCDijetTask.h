#ifndef ALIJCDIJETTASK_H
#define ALIJCDIJETTASK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Analysis task for providing various flow informations 
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


using namespace std;
//==============================================================
class TClonesArray;
class AliJCDijetHistos;

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
  void    SetJetConeSize(double jetCone) {fjetCone=jetCone; }
  Bool_t  IsMC()const{ return fIsMC; }
  void    SetIsMC(Bool_t b) { fIsMC=b; }
  void    SetCuts(double particleEta, double particlePt, double leadingJet, double subleadingJet, double constituent, double deltaPhi) {fparticleEtaCut=particleEta; fparticlePtCut=particlePt; fleadingJetCut=leadingJet; fsubleadingJetCut=subleadingJet; fconstituentCut=constituent; fdeltaPhiCut=deltaPhi;}
  //double  GetLeadingJetCut() {return fleadingJetCut;}
  //double  GetSubleadingJetCut() {return fsubleadingJetCut;}
  //double  GetConstituentCut() {return fconstituentCut;}
  void CalculateJetsDijets(TClonesArray *inList,
                           int lDebug,
                           int lCBin,
                           double lParticleEtaCut,
                           double lParticlePtCut,
                           double lJetCone,
                           double lConstituentCut,
                           double lLeadingJetCut,
                           double lSubleadingJetCut,
                           double lDeltaPhiCut);



  // Methods specific for this class
  void SetJCatalystTaskName(TString name){ fJCatalystTaskName=name; } // Setter for filter task name

 private:

  AliJCatalystTask *fJCatalystTask;  //
  TString           fJCatalystTaskName; // Name for JCatalyst task
  vector<double> fcentralityBins;
  double fjetCone;
  Bool_t      fIsMC;       // MC data or real data
  double fparticleEtaCut;
  double fparticlePtCut;
  double fleadingJetCut;
  double fsubleadingJetCut;
  double fconstituentCut;
  double fdeltaPhiCut;
  AliJCDijetHistos *fhistos;
  int fCBin;
  TDirectory     *fOutput; // Output directory

  ClassDef(AliJCDijetTask, 1); 
};
#endif // ALIJCDIJETTASK_H
