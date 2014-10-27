#ifndef ALIANALYSISTASKCPQA_H
#define ALIANALYSISTASKCPQA_H

class TList;
class TH1F;
class TH2F;
class TH3F;
class AliTriggerAnalysis;


#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCPQA : public AliAnalysisTaskSE {
 public:

  AliAnalysisTaskCPQA(const char *name="<default name>");
  virtual ~AliAnalysisTaskCPQA() ;// { /*if (fOutputList) delete fOutputList;*/}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t* option);
  virtual void   Terminate(Option_t *);

  void LoopESD();
  void LoopESDMC();

  
  void UseMC(Bool_t useMC=kTRUE) { fUseMC = useMC;}
  
 private:
  Bool_t       fUseMC; //MC flag
  AliESDEvent *fESD; // esd event
  TList	      *fOutputList; //! output list

  TH1F *fhEvent;//!


//  Double_t fEtaMaxM;
//  Double_t fEtaMaxD;
//  Double_t fVtxZmax;

  TH2F *fhV0A[4];//!
  TH2F *fhV0C[4];//!
  TH2F *fhV0online[4];//!
  TH2F *fhV0offline[4];//!
  TH1F *fhSPDFiredChip[4];//!
  TH1F *fhSPDFastOrChip[4];//!
  TH1F *fhReferenceMultiplicity[4];//!
  TH3F *fhVtxTrack[4];//!

  AliTriggerAnalysis * fTriggerAnalysis; // trigger analysis object, to get the offline triggers

  TH1F* Hist1D(const char* name, Int_t nBins, Double_t xMin, Double_t xMax,  const char* xLabel="", Int_t color=1, Int_t ls=1, const char* yLabel="");
  TH2F *Hist2D(const char* name, Int_t nBinsx, Double_t xMin, Double_t xMax, Int_t nBinsy, Double_t yMin, Double_t yMax, const char* xLabel="", const char* yLabel="", Int_t color=1);
  TH3F *Hist3D(const char* name, Int_t nBinsx, Double_t xMin, Double_t xMax, Int_t nBinsy, Double_t yMin, Double_t yMax,  Int_t nBinsz, Double_t zMin, Double_t zMax, const char* xLabel="", const char* yLabel="", const char *zLabel="");

  // public:
    AliAnalysisTaskCPQA(const AliAnalysisTaskCPQA&); // not implemented
    AliAnalysisTaskCPQA& operator=(const AliAnalysisTaskCPQA&); // not implemented
  
  ClassDef(AliAnalysisTaskCPQA, 1);
};

#endif
