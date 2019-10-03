#ifndef AliEbyEHigherMomentsTaskPID_cxx
#define AliEbyEHigherMomentsTaskPID_cxx

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
class AliPIDResponse;
class AliHelperPID;
class TList;

#include "TParticle.h"
#include "AliAnalysisTaskSE.h"

class AliEbyEHigherMomentsTaskPID: public AliAnalysisTaskSE {
 public:
  AliEbyEHigherMomentsTaskPID( const char *name = "HigherMomentAnalysis");
  virtual ~AliEbyEHigherMomentsTaskPID();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   doAODEvent();
  virtual void   doMCAODEvent();
  
  virtual void   Terminate(Option_t *);
  
  void SetVertexDiamond(Double_t vx, Double_t vy, Double_t vz) {fVxMax = vx;fVyMax = vy; fVzMax = vz;}
  void SetCentralityEstimator(const char* centralityEstimator) { fCentralityEstimator = centralityEstimator;}
  void SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
  void SetRapidityCut(Double_t rapidity){ fRapidityCut = rapidity;}
  void SetNSigmaCut(Double_t nsigma){ fNSigmaCut = nsigma;}
  void SetParticleSpecies(AliPID::EParticleType pid) {fParticleSpecies = pid;}
  void SetAODtrackCutBit(Int_t bit){ fAODtrackCutBit = bit;}
  void SetHelperPID(AliHelperPID* pid){fHelperPID = pid;}
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
  AliPIDResponse	*fPIDResponse;
  AliPID::EParticleType fParticleSpecies;
 
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
  Double_t fRapidityCut;
  Double_t fNSigmaCut;
  Int_t fAODtrackCutBit;//track cut bit from track selection (only used for AODs)
  AliHelperPID *fHelperPID;
  TH1D *fEventCounter;
  
  TH1D *fHistQA[4];
  TH2D *fHistVxVy;
  
  THnSparse *fTHnCentNplusNminusPid;
  THnSparse *fTHnCentNplusNminusPidTruth;
  THnSparse *fPtBinNplusNminusPid;
  THnSparse *fPtBinNplusNminusPidTruth;
  
  AliEbyEHigherMomentsTaskPID(const AliEbyEHigherMomentsTaskPID&);
  AliEbyEHigherMomentsTaskPID& operator = (const AliEbyEHigherMomentsTaskPID&);//Not implimented..
  ClassDef(AliEbyEHigherMomentsTaskPID, 1);

};

#endif

 
