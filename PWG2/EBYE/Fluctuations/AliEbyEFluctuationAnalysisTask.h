#ifndef AliEbyEFluctuationAnalysisTask_cxx
#define AliEbyEFluctuationAnalysisTask_cxx

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-------------------------------------------------------------------------
//                 AliEbyEFluctuationAnalysisTask
//   This is the class to deal with the EbyE Charge Fluctuation  analysis
//               Origin: Satyajit Jena, sjena@cern.ch
//-------------------------------------------------------------------------

class TH1F;
class TH2F;
class TList;
class TTree;
class AliESDVertex;
class AliESDtrackCuts;
class AliESDEvent;

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif


class AliEbyEFluctuationAnalysisTask : public AliAnalysisTaskSE {
 public:
    AliEbyEFluctuationAnalysisTask();
    AliEbyEFluctuationAnalysisTask(const char *name);
    virtual ~AliEbyEFluctuationAnalysisTask();
    
    virtual void     UserCreateOutputObjects();
    virtual void     UserExec(Option_t *option);
    virtual void     Terminate(Option_t *);
    
    void SetAnalysisCutObject(AliESDtrackCuts *const trackCuts) {
      fTrackCut = trackCuts;}
    void SetCentralityType(const char *centtype, 
			   const char *centmode, 
			   Int_t cbin) { //To be removed
      fCentType = centtype; 
      fCentMode = centmode; 
      fCentBin = cbin;
    } 

    void SetAnalysisType(const char *analysistype) { fAnalysisType = analysistype;}
    void SetAnalysisMode(const char *analysisMode) { fAnalysisMode = analysisMode;}
    void SetCentralityEstimator(const char *centestimator) { fCentEstimator = centestimator;}
    void SetCentType(const char *centtype) { fCentType = centtype;}

    void SetVertex(Double_t vx, Double_t vy, Double_t vz) {fVx = vx; fVy = vy; fVz = vz; }

    Int_t GetCentrality(AliESDEvent *esd) const;
    Int_t FindCentralityESD(Double_t mult) const;
    Int_t FindCentralityMC(Double_t mult) const;
 
    const AliESDVertex *GetVertex(AliESDEvent* esd);

 private:
    AliESDtrackCuts *fTrackCut;//cut object
    
    TList   *fListQA;//! Output QA list
    TList   *fListResults;//!Output results list
    TString  fAnalysisType;//"ESD","MC","AOD"
    TString  fAnalysisMode;//"TPC","Global","Hybrid"
    TString  fCentType;//to be removed
    Int_t    fCentBin;
    TString  fCentMode;
    TString  fCentEstimator;//Centrality estimator "V0M","TKL","TRK","FMD","CL0","CL1"
    Double_t fVx;//Vx
    Double_t fVy;//Vy
    Double_t fVz;//Vz
     
    TH1F  *fHistEventStats;//event stats
    TH2F  *fhVertex; // 1-vx,2-vy,3-vz,4-vxacc,5-vyAcc,6-vzAcc
    TH2F  *fhVzeroMult;//V0 amplitude
    TH2F  *fCentVsMult; //centraltiy vs multiplicity
    TH2F  *fhVnT;//
    TH2F  *fHCentrality[20];//nPos vs nNeg vs centrality

    //_________________________________________________
    AliEbyEFluctuationAnalysisTask(const AliEbyEFluctuationAnalysisTask&); // not implemented
    AliEbyEFluctuationAnalysisTask& operator=(const AliEbyEFluctuationAnalysisTask&); // not implemented
    
    ClassDef(AliEbyEFluctuationAnalysisTask, 1); 
};

#endif

