//--- Task for investigation of the DoubleHyperHydrogen4 ---
//---     Author: Janik Ditzel; janik.ditzel@cern.ch     ---

#ifndef ALIANALYSISTASKS3PARTICLEYIELDS_H
#define ALIANALYSISTASKS3PARTICLEYIELDS_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"

class AliAnalysisTaskS3ParticleYields : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskS3ParticleYields();
  AliAnalysisTaskS3ParticleYields(const char *name);
  virtual ~AliAnalysisTaskS3ParticleYields();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(const Option_t*);
  void SelectPIDcheckOnly(Bool_t pidch = kFALSE) {fPIDCheckOnly = pidch;};
  void SetTriggerMask(UInt_t triggerMask = AliVEvent::kINT7) {fTriggerMask = triggerMask;};
  enum PdgCodeType_t {
    kPDGPionPlus,
    kPDGPionMinus,
    kPDGProton,
    kPDGAntiProton,
    kPDGKaonMinus,
    kPDGKaonPlus,
    kPDGDeuteron,
    kPDGAntiDeuteron,
    kPDGTriton,
    kPDGAntiTriton,
    kPDGHelium3,
    kPDGAntiHelium3,
    kPDGHelium4,
    kPDGAntiHelium4,
    kPDGLambda,
    kPDGAntiLambda,
    kPDGDeltaMinus,
    kPDGOmegaMinus,
    kPDGOmegaPlus,
    kPDGXiMinus,
    kPDGXiPlus,
    kPDGOmegaOmega,
    kPDGAntiOmegaOmega,
    kPDGLambdaNeutronNeutron,
    kPDGAntiLambdaNeutronNeutron,
    kPDGXi0Proton,
    kPDGAntiXi0Proton,
    kPDGLambdaN,
    kPDGAntiLambdaN,
    kPDGOmegaProton,
    kPDGAntiOmegaProton,
    kPDGLambdaLambda,
    kPDGAntiLambdaLambda,
    kPDGHyperHelium4,
    kPDGAntiHyperHelium4,
    kPDGHyperHelium5,
    kPDGAntiHyperHelium5,
    kPDGDoubleHyperHydrogen4,
    kPDGAntiDoubleHyperHydrogen4,
    kPDGHyperHydrogen3,
    kPDGAntiHyperHydrogen3,
    kPDGHyperHydrogen4,
    kPDGAntiHyperHydrogen4,
 
    kPdgCodeMax  // Number of enum entries                                                                                                                                          
  };
  static const Int_t fgkPdgCode[];
  
 private:
  AliESDInputHandler    *fInputHandler;
  AliESDpid             *fPID;
  AliESDEvent           *fESDevent;
  AliStack              *fStack;
  AliESDv0              *fV0;
  AliMCEvent            *mcEvent;
  AliESDtrackCuts       *trackCutsV0;
  TH2F                  *fHistdEdx; 
  TH2F                  *fHistdEdxV0;  
  TH1F                  *fHistNumEvents; 
  TH1F			*fHistTrigger;	
  TH1F			*fHistV0;
  TH1F                  *fHistEvents;
  TTree                 *fTree;   
  TTree                 *hTree; 
  TTree                 *fTreeGen; 
  TList                 *fHistogramList; 
  TVector3              fPrimaryVertex;
  Double_t              fMagneticField;
  Int_t                 fNV0Cand;
  Bool_t                fPIDCheckOnly; 
  Bool_t                fMCtrue;
  AliEventCuts          fEventCuts; 
  UInt_t		fTriggerMask;	
  Int_t		        fPeriod;   
  Bool_t                fBetheSplines;   
  Double_t              fBetheParamsHe[6];   
  Double_t              fBetheParamsT[6]; 
  Int_t                 fonTheFly;
  Int_t                 frunnumber;
Int_t ftrig, fz, fmc;
  Float_t fmLambda, fpLambda, fptLambda, fctLambda, fdcaLambda, fcosLambda, fyLambda;
  Float_t fpy, fhe3y, fpLy, fpiP, fhe3P, fhe3Pt, fpP, fpPt, fpLP, fpiy, fpchi2, fhe3chi2;
  Float_t fpDcaSec, fpiDcaSec, fpiDca, fpDca, fpLDca, fpLDcaSec, fpDcaz, fhe3Dcaz, fhe3Dca;
  Float_t fpiNcls, fhe3Ncls, fpNcls, fpLNcls, fpiNclsITS, fhe3NclsITS, fpNclsITS, fpLNclsITS;
  Float_t fpiDedxSigma, fhe3DedxSigma, fpDedxSigma, fpLDedxSigma, fpiDedx, fhe3Dedx, fpDedx, fpLDedx;
  Float_t farmalpha, farmpt;
  Float_t fthetaP, fthetaN, fEtaHe3, fEtaP, fEtaPL, fEtaPi, fPhiHe3, fPhiP, fPhiPL, fPhiPi;
  Float_t fGeoLengthHe3, fGeoLengthP, fGeoLengthPi, fGeoLengthPL, fTOFSignalHe3, fTOFSignalP, fTOFSignalPi, fTOFSignalPL;  
  Float_t fpHe3Gen, fyHe3Gen, fpPGen, fyPGen, fpLambdaGen, fyLambdaGen, fmLambdaGen;
Int_t fMCtrueHe3, fisPrimaryHe3, fisWeakHe3, fisMaterialHe3, fisfromHypertriton, fisPrimaryP, fisWeakP, fisMaterialP, fMCtrueP, fMCtrueL, fisWeakL, fisMaterialL, fisPrimaryL;
  Int_t fisPrimaryGenHe3, fisSecondaryGenHe3, fisPrimaryGenP, fisMaterialGenP, fisSecondaryGenP, fisMaterialGenHe3, fisWeakGenL, fisMaterialGenL, fisPrimaryGenL;
  Int_t fHe3Charge, fPCharge, fLambdaCharge;

  TVector3              fVertexPosition; 
  Int_t              fNumberV0s;
  Int_t                 fCentrality;
  Int_t              fTrigger; 
  TString               fTriggerClasses;

  Int_t fMultV0M, fMultOfV0M, fMultSPDTracklet, fMultSPDCluster, fMultRef05, fMultRef08, tSPDCluster, tSPDTracklets, tSPDFiredChips0, tSPDFiredChips1, tV0Multiplicity;

  Bool_t TriggerSelection();
  Double_t Bethe(const AliESDtrack& track, Double_t mass, Int_t charge, Double_t* params);
  Double_t GeoLength(const AliESDtrack& track);
  Double_t TOFSignal(const AliESDtrack& track);
  void dEdxCheck();
  void LambdaV0s();
  void ProtonTracks();
  void MCGenerated();
  void SetMultiplicity();
  void SetBetheBlochParams(Int_t runnumber);
  void ResetVals(TString mode);
  Int_t GetLabel(Int_t labelFirstMother, Int_t particlePdgCode);
  Float_t GetInvPtDevFromBC(Int_t b, Int_t c);
  AliAnalysisTaskS3ParticleYields(const AliAnalysisTaskS3ParticleYields&);
  AliAnalysisTaskS3ParticleYields &operator=(const AliAnalysisTaskS3ParticleYields&);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskS3ParticleYields, 1);
  /// \endcond
};



#endif
