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
class TH2D;
class TH3D;
class THnSparse;
class TString;
class AliAODEvent;
class AliAODTrack;
class TList;

#include "TParticle.h"
#include "AliAnalysisTaskSE.h"

class AliEbyEHigherMomentsTask: public AliAnalysisTaskSE {
 public:
  AliEbyEHigherMomentsTask( const char *name = "HigherMomentAnalysis");
  virtual ~AliEbyEHigherMomentsTask();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   doAODEvent();
  virtual void   doMCAODEvent();
  
  virtual void   Terminate(Option_t *);
  
  void SetVertexDiamond(Double_t vx, Double_t vy, Double_t vz) {fVxMax = vx;fVyMax = vy; fVzMax = vz;}
  void SetCentralityEstimator(const char* centralityEstimator) { fCentralityEstimator = centralityEstimator;}
  void SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
  void SetAODtrackCutBit(Int_t bit){ fAODtrackCutBit = bit;}
  void SetKinematicsCutsAOD(Double_t ptl, Double_t pth, Double_t eta){
    fPtLowerLimit = ptl;
    fPtHigherLimit = pth;
    fEtaLowerLimit = -1.*eta;
    fEtaHigherLimit = eta;
    
  }
  void SetNumberOfPtBins(Int_t nPtBins){ fNptBins = nPtBins;}
  
  
  
 private:
  
  Bool_t ProperVertex(AliAODEvent *fAOD) const;
  Bool_t AcceptTrack(AliAODTrack* track) const;
  Int_t  GetPtBin(Double_t pt);


  TList *fListOfHistos;
  TClonesArray          *fArrayMC;
  TString          fAnalysisType;          // "MC", "ESD", "AOD"
  TString          fCentralityEstimator;   // "V0M","TRK","TKL","ZDC","FMD"
  
  Int_t fCentrality;
  Double_t fVxMax;               //vxmax
  Double_t fVyMax;//vymax
  Double_t fVzMax;//vzmax
  Double_t fPtLowerLimit;
  Double_t fPtHigherLimit;
  Int_t fNptBins;
  Int_t fBin;
  Double_t fEtaLowerLimit;
  Double_t fEtaHigherLimit;
  Int_t fAODtrackCutBit;//track cut bit from track selection (only used for AODs)
  TH1D *fEventCounter;
  
  TH1D *fHistQA[4];
  TH2D *fHistVxVy;
  
  
  THnSparse *fTHnCentNplusNminusCh;
  THnSparse *fTHnCentNplusNminusChTruth;
  THnSparse *fPtBinNplusNminusCh;
  THnSparse *fPtBinNplusNminusChTruth;
  
  
  
  AliEbyEHigherMomentsTask(const AliEbyEHigherMomentsTask&);
  AliEbyEHigherMomentsTask& operator = (const AliEbyEHigherMomentsTask&);//Not implimented..
  ClassDef(AliEbyEHigherMomentsTask, 1);

};

#endif

 
