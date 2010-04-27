#ifndef ALIANALYSISTASKDIJETS_H
#define ALIANALYSISTASKDIJETS_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>

class AliAnalysisTaskDiJets : public AliAnalysisTaskSE
{
 public:
    AliAnalysisTaskDiJets();
    AliAnalysisTaskDiJets(const char* name);
    virtual ~AliAnalysisTaskDiJets() {;}

    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);

    void SetFillAOD(Bool_t fill) { fFillAOD=fill; }

 private:
  AliAnalysisTaskDiJets(const AliAnalysisTaskDiJets &det);
  AliAnalysisTaskDiJets &operator=(const AliAnalysisTaskDiJets &det);

 private:
  TClonesArray* fDiJets;    // Array of dijets
  TClonesArray* fDiJetsIn;  // Array of dijets

  Bool_t        fFillAOD;   // option to fill AOD branch

  TList        *fHistList;  // Output list

  TH1F         *fH1DeltaPt;  // Pt difference
  TH1F         *fH1DeltaPhi; // delta phi plot in (0,pi)
  TH1F         *fH1PhiImbal; // phi imbalance (-pi,pi)
  TH1F         *fH1Asym;     // asymmetry of the dijet
  TH2F         *fH2Pt2vsPt1; // scatter plot with the two jets' Pt
  TH2F         *fH2DifvsSum; // Pt difference vs Pt sum

  ClassDef(AliAnalysisTaskDiJets, 2); // Analysis task for standard dijet analysis
};

#endif
