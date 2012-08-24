#ifndef AliEbyEHigherMomentsTask_cxx
#define AliEbyEHigherMomentsTask_cxx

//=========================================================================//
//                                                                         //
//           Analysis Task for Net-Charge Higher Moment Analysis           //
//              Author: Satyajit Jena || Nirbhay K. Behera                 //
//                      sjena@cern.ch || nbehera@cern.ch                   //
//                               V0.0 23/08/2012                           //
//                                                                         //
//=========================================================================//


class TH1D;
class TH2F;
class TString;
class AliAODEvent;
class TList;


#include "AliAnalysisTaskSE.h"

class AliEbyEHigherMomentsTask: public AliAnalysisTaskSE {
 public:
  AliEbyEHigherMomentsTask( const char *name = "HigherMomentAnalysis");
  virtual ~AliEbyEHigherMomentsTask();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetVertexDiamond(Double_t vx, Double_t vy, Double_t vz) {fVxMax = vx;fVyMax = vy; fVzMax = vz;}
  void SetCentralityEstimator(const char* centralityEstimator) { fCentralityEstimator = centralityEstimator;}
  void SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
  void SetDCA(Double_t xy, Double_t z) { fDCAxy = xy; fDCAz = z; }
  void SetPtRange(Double_t ptl, Double_t pth){fPtLowerLimit = ptl; fPtHigherLimit = pth;}
  void SetEta(Double_t eta){fEtaLowerLimit=-eta; fEtaHigherLimit=eta;}
  void SetTPCNclus(Int_t nclus) { fTPCNClus = nclus;}
  void SetAODtrackCutBit(Int_t bit){
    nAODtrackCutBit = bit;
  }
  
  void SetKinematicsCutsAOD(Double_t ptl, Double_t pth, Double_t eta){
    fPtLowerLimit = ptl;
    fPtHigherLimit = pth;
    fEtaLowerLimit = -eta;
    fEtaHigherLimit = eta;
    
  }
 
 
 private:
  TList *fListOfHistos;
  
  TString          fAnalysisType;          // "MC", "ESD", "AOD"
  TString          fCentralityEstimator;   // "V0M","TRK","TKL","ZDC","FMD"

  Double_t fVxMax;               //vxmax
  Double_t fVyMax;//vymax
  Double_t fVzMax;//vzmax

  Double_t fDCAxy;
  Double_t fDCAz;
  Double_t fPtLowerLimit;
  Double_t fPtHigherLimit;
  Double_t fEtaLowerLimit;
  Double_t fEtaHigherLimit;
  Int_t fTPCNClus;
  Int_t nAODtrackCutBit;//track cut bit from track selection (only used for AODs)
  TH1D *fEventCounter;

  TH2F *fhNplusNminus[91]; //Data 
  
  TH1D *fHistQA[11];
  
  AliEbyEHigherMomentsTask(const AliEbyEHigherMomentsTask&);
  AliEbyEHigherMomentsTask& operator = (const AliEbyEHigherMomentsTask&);//Not implimented..
  ClassDef(AliEbyEHigherMomentsTask, 1);

};

#endif

 
