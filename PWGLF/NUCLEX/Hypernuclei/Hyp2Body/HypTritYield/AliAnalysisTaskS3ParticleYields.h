//--- Task for investigation of the DoubleHyperHydrogen4 ---
//---     Author: Janik Ditzel; janik.ditzel@cern.ch     ---

#ifndef ALIANALYSISTASKS3PARTICLEYIELDS_H
#define ALIANALYSISTASKS3PARTICLEYIELDS_H

class TH1F;
class TH2F;
class AliESDEvent;
class AliESDpid;
class AliESDtrackCuts;
class AliESDv0;
class AliESDVertex;
class AliESDInputHandler;
class AliESDtrack;

#include "AliAnalysisTaskSE.h"
#include "AliStack.h"
#include "AliEventCuts.h"
#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <TClonesArray.h>
#include "AliPID.h"
#include "AliVertexerTracks.h"
#include "THnSparse.h"


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
  AliESDInputHandler    *fInputHandler;        //!<! Input handler
  AliESDpid             *fPID;                 //!<! ESD pid
  AliESDEvent           *fESDevent;            //!<! ESD event
  AliStack              *fStack;               //!<! MC stack
  AliESDv0              *fV0;                  //!<! ESD v0 - He4 + pi
  TH2F                  *fHistdEdx;            //<   Histogram of Tpc dEdx for pid qa
  THnSparseF		*fHistData;
  THnSparseF		*fHistMC;
  TH2F                  *fHistdEdxV0;          //<   Histogram of Tpc dEdx for pid qa
  TH1F                  *fHistNumEvents;       //<   Histogram of number of events
  TH1F			            *fHistTrigger;	 	//<   Histogram of trigger for all events 
  TH1F			            *fHistV0;	 	//<   Histogram of trigger for all V0s 
  TH1F                  *fHistEvents;
  TTree                 *fTree;                //<   Tree containing reduced events
  TTree                 *gTree;                //<   Tree containing reduced events
  TTree                 *hTree;                //<   Tree containing reduced events
  TTree                 *fTreeGen;                //<   Tree containing reduced events
  TList                 *fHistogramList;       //<   List of histograms
  TVector3              fPrimaryVertex;       //!<! Vector of primary vertex of collision
  Double_t              fMagneticField;       //!<! Magnetic field
  Int_t                 fNV0Cand;             //!<! Number of V0 candidates in a event
  Bool_t                fPIDCheckOnly;        //< Flag to reduce the task only to PID check for Hypertriton daughters
  Bool_t                fMCtrue;              //< Flag for activating MC analysis (set automatically)
  AliEventCuts          fEventCuts;           //< 2015 event cuts as advised by PDG (AliEventCuts)
  UInt_t		            fTriggerMask;		//< Triggermask for event cuts
  Int_t		              fPeriod;              //< Data period for centrality selector
  Bool_t                fBetheSplines;        //< Switch between built in BetheSplines and personal Fit
  Double_t              fBetheParamsHe[6];    //< Bethe Aleph He3 Parameter + TPC sigma: [0][i] he3 [2][i] t
  Double_t              fBetheParamsT[6];     //< Bethe Aleph He3 Parameter + TPC sigma: [0][i] he3 [2][i] t
  Int_t                 fonTheFly;
  Int_t                 frunnumber;
  Float_t fmLambda, fpLambda, fptLambda, fctLambda, fdcaLambda, fcosLambda, fyLambda;
  Float_t fpy, fhe3y, fpLy, fpiP, fhe3P, fhe3Pt, fpP, fpPt, fpLP, fpiy, fpchi2, fhe3chi2;
  Float_t fpDcaSec, fpiDcaSec, fpiDca, fpDca, fpLDca, fpLDcaSec, fpDcaz, fhe3Dcaz, fhe3Dca;
  Float_t fpiNcls, fhe3Ncls, fpNcls, fpLNcls, fpiNclsITS, fhe3NclsITS, fpNclsITS, fpLNclsITS;
  Float_t fpiDedxSigma, fhe3DedxSigma, fpDedxSigma, fpLDedxSigma, fpiDedx, fhe3Dedx, fpDedx, fpLDedx;
  Float_t farmalpha, farmpt;
  Int_t ftrig, fz, fmc;
  Float_t fthetaP, fthetaN, fEtaHe3, fEtaP, fEtaPL, fEtaPi, fPhiHe3, fPhiP, fPhiPL, fPhiPi;
  Float_t fGeoLengthHe3, fGeoLengthP, fGeoLengthPi, fGeoLengthPL, fTOFSignalHe3, fTOFSignalP, fTOFSignalPi, fTOFSignalPL;
  Int_t fMCtrueHe3, fisPrimaryHe3, fisWeakHe3, fisMaterialHe3, fisfromHypertriton, fisPrimaryP, fisWeakP, fisMaterialP, fMCtrueP, fMCtrueL;
  Float_t fpHe3Gen, fyHe3Gen, fpPGen, fyPGen, fpLambdaGen, fyLambdaGen;
  Int_t fisPrimaryGenHe3, fisSecondaryGenHe3, fisPrimaryGenP, fisMaterialGenP, fisSecondaryGenP, fisMaterialGenHe3;
  Int_t fHe3Charge, fPCharge, fLambdaCharge;

  TVector3              fVertexPosition; //< position of primary vertex
  Int_t              fNumberV0s;      //< number of v0s in event
  Int_t                 fCentrality;     //< centrality of event
  Int_t              fTrigger;        //< array of Triggers
  TString               fTriggerClasses; //< fired trigger classes

  Int_t fMultV0M, fMultOfV0M, fMultSPDTracklet, fMultSPDCluster, fMultRef05, fMultRef08, tSPDCluster, tSPDTracklets, tSPDFiredChips0, tSPDFiredChips1, tV0Multiplicity;

  Bool_t TriggerSelection();
  Double_t Bethe(const AliESDtrack& track, Double_t mass, Int_t charge, Double_t* params);
  Double_t GeoLength(const AliESDtrack& track);
  Double_t TOFSignal(const AliESDtrack& track);
  void dEdxCheck();
  void V0Analysis(AliESDtrackCuts trackCutsV0, AliMCEvent* mcEvent);
  void He3PYields(AliESDtrackCuts trackCutsV0, AliMCEvent* mcEvent);
  void MCGenerated(AliMCEvent* mcEvent);
  void SetMultiplicity();
  void SetBetheBlochParams(Int_t runnumber);
  Float_t GetInvPtDevFromBC(Int_t b, Int_t c);
  AliAnalysisTaskS3ParticleYields(const AliAnalysisTaskS3ParticleYields&);
  AliAnalysisTaskS3ParticleYields &operator=(const AliAnalysisTaskS3ParticleYields&);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskS3ParticleYields, 1);
  /// \endcond
};



#endif
