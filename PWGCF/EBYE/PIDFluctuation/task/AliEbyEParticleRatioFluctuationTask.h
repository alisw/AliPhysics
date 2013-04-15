#ifndef AliEbyEParticleRatioFluctuationTask_cxx
#define AliEbyEParticleRatioFluctuationTask_cxx

//=========================================================================//
//                                                                         //
//             AliEbyE Analysis for Particle Ratio Fluctuation             //
//              Author:   Deepika Rathee  || Satyajit Jena                 //
//                        drathee@cern.ch || sjena@cern.ch                 //
//                                                                         //
//=========================================================================//

class TH1D;
class TH2F;
class TH3F;
class TString;
class AliAODEvent;
class AliAODTrack;
class AliAODMCParticle;
class TList;
class AliESDtrackCuts;
class AliHelperPID;

#include "AliAnalysisTaskSE.h"
#include "AliPID.h"
#include "THnSparse.h"

class AliEbyEParticleRatioFluctuationTask: public AliAnalysisTaskSE {
 public:
  AliEbyEParticleRatioFluctuationTask( const char *name = "HigherMomentAnalysis");
  virtual ~AliEbyEParticleRatioFluctuationTask();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  static const Int_t kNCentralityBins = 20; //! N centrality bins
  static const Int_t kNSparseData = 14;     //! N Sparse bins
   
  void SetVertexDiamond(Double_t vx, Double_t vy, Double_t vz) {fVxMax = vx;fVyMax = vy; fVzMax = vz;}
  void SetKinematicsCuts(Double_t ptl, Double_t pth, Double_t eta) {fPtLowerLimit = ptl; fPtHigherLimit = pth; fEtaLowerLimit = -eta; fEtaHigherLimit = eta; }
  void SetAODtrackCutBit(Int_t bit) {fAODtrackCutBit = bit; }
  void SetDCA(Double_t xy, Double_t z) { fDCAxy = xy; fDCAz = z; }
  void SetTPCNclus(Int_t nclus) { fTPCNClus = nclus;}
  void SetCentralityEstimator(const char* cent) { fCentralityEstimator = cent;}
  void SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
  void SetAnalysisData(const char* analysisData) {fAnalysisData = analysisData;}
  void RunQA() {isQA = kTRUE;}
  void Debug() {fDebug = kTRUE;}
  void SetHelperPID(AliHelperPID* pid){fHelperPID = pid;}

 private:

  Bool_t       AcceptEvent(AliAODEvent *event) const; //! accept eventc
  Bool_t       AcceptTrack(AliAODTrack *track) const; //! accept track
  Bool_t       AcceptMCTrack(AliAODMCParticle *track) const; //! accept track
  TList        *fThnList;
  TString      fAnalysisType;          //! "PbPb", "pp", "pA"
  TString      fAnalysisData;          //! "AOD", "AODMC"
  TString      fCentralityEstimator;   //! "V0M","TRK","TKL","ZDC","FMD"  
  Double_t     fVxMax;                 //!  vxmax
  Double_t     fVyMax;                 //!  vymax
  Double_t     fVzMax;                 //!  vzmax
  Double_t     fDCAxy;
  Double_t     fDCAz;
  Double_t     fPtLowerLimit;
  Double_t     fPtHigherLimit;
  Double_t     fEtaLowerLimit;
  Double_t     fEtaHigherLimit;
  Int_t        fTPCNClus;
  Int_t        fAODtrackCutBit;
  Bool_t       isQA;
  Bool_t       fDebug;
  AliHelperPID *fHelperPID;
  TH1D         *fEventCounter;
  TH1D         *fHistQA[14];
  THnSparseI   *fHistoCorrelation; 
  
  //________________________________
  AliEbyEParticleRatioFluctuationTask(const AliEbyEParticleRatioFluctuationTask&);
  AliEbyEParticleRatioFluctuationTask& operator = (const AliEbyEParticleRatioFluctuationTask&);
  ClassDef(AliEbyEParticleRatioFluctuationTask, 1);

};

#endif

 
