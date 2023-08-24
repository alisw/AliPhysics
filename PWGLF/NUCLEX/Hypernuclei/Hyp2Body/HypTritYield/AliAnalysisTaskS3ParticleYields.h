//---     Author: Janik Ditzel; janik.ditzel@cern.ch     ---

#ifndef ALIANALYSISTASKS3PARTICLEYIELDS_H
#define ALIANALYSISTASKS3PARTICLEYIELDS_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliTRDonlineTrackMatching.h"
#include "AliCDBManager.h"
#include "AliGeomManager.h"

class AliAnalysisTaskS3ParticleYields : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskS3ParticleYields();
  AliAnalysisTaskS3ParticleYields(const char *name);
  virtual ~AliAnalysisTaskS3ParticleYields();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(const Option_t*);
  void SelectPIDcheckOnly(Bool_t pidch = kFALSE) {fPIDCheckOnly = pidch;};
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
  TTree                 *fTree;   
  TTree                 *hTree; 
  TTree                 *fTreeGen; 
  TList                 *fHistogramList; 
  TVector3              fPrimaryVertex;  
  Bool_t                fPIDCheckOnly; 
  Bool_t                fMCtrue;
  AliEventCuts          fEventCuts; 	   
  Double_t              fBetheParamsHe[6];   
  Double_t              fBetheParamsT[6]; 
  Int_t                 MB, HMV0, HMSPD, HNU, HQU;
  Int_t                 fYear;  // Jahr für  OCDB
  //saved in Tree
  Int_t                 fNumberV0s; 
  Int_t                 fNTracks; 
  Int_t                 fonTheFly;
  Int_t                 frunnumber;
  Int_t                 fCharge;
  Int_t			fTrigMB, fTrigHMV0, fTrigHMSPD, fTrigHNU, fTrigHQU;
  Int_t                 fMultV0M, fMultOfV0M, fMultSPDTracklet, fMultSPDCluster, fMultRef05, fMultRef08, fSPDCluster, fSPDTracklets, fSPDFiredChips0, fSPDFiredChips1, fV0Multiplicity;
  Int_t                 fisPrimaryP, fisWeakP, fisMaterialP, fMCtrueP, fMCtrueL, fisWeakL, fisMaterialL, fisPrimaryL;
  Int_t                 fKinkP, fTPCrefitP, fITSrefitP, fKinkPi, fTPCrefitPi, fITSrefitPi, fpiNcls, fpiNclsITS, fpNcls, fpNclsITS;
  Float_t               fITSchi2P, fITSchi2Pi;
  Float_t               fLambdaM, fLambdaE, fLambdaP, fLambdaPx, fLambdaPy, fLambdaPz, fLambdaPt, fLambdaCt, fLambdaDca, fLambdaCos, fLambdaY;
  Float_t               fpiP, fpiPt, fpiPx, fpiPy, fpiPz, fpiE, fpiy, fpichi2, fpiDcaSec, fpiDca, fpiDcaz, fpiDedxSigma, fpiDedx, fEtaPi, fPhiPi, fGeoLengthPi, fTOFSignalPi, fTOFnSigmaPi, fpiSigmaYX, fpiSigmaXYZ, fpiSigmaZ;
  Float_t               fpy, fpP, fpPx, fpPy, fpPz, fpPt, fpE, fpDcaSec, fpDca, fpDcaz, fpchi2, fpDedxSigma, fpDedx, fEtaP, fGeoLengthP, fPhiP, fTOFSignalP, fTOFnSigmaP, fpSigmaYX, fpSigmaXYZ, fpSigmaZ;
  Float_t               farmalpha, farmpt, fthetaP, fthetaN; 
  Int_t                 fTRDvalidP; // has valid TRD track
  Int_t                 fTRDtrigHNUP; // HNU fired by track
  Int_t                 fTRDtrigHQUP; // HQU fired by track
  Int_t                 fTRDPidP;
  Int_t                 fTRDnTrackletsP;
  Int_t                 fTRDPtP;
  Int_t                 fTRDLayerMaskP;
  Float_t               fTRDSagittaP;
  Int_t                 fTRDvalidPi; // has valid TRD track
  Int_t                 fTRDtrigHNUPi; // HNU fired by track
  Int_t                 fTRDtrigHQUPi; // HQU fired by track
  Int_t                 fTRDPidPi;
  Int_t                 fTRDnTrackletsPi;
  Int_t                 fTRDPtPi;
  Int_t                 fTRDLayerMaskPi;
  Float_t               fTRDSagittaPi;
  //Functions
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
  Double_t TRDtrack(AliESDtrack* esdTrack, Int_t particleID);
  //
  AliAnalysisTaskS3ParticleYields(const AliAnalysisTaskS3ParticleYields&);
  AliAnalysisTaskS3ParticleYields &operator=(const AliAnalysisTaskS3ParticleYields&);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskS3ParticleYields, 1);
  /// \endcond
};



#endif
