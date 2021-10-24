//--- Task for investigation of the DoubleHyperHydrogen4 ---
//---     Author: Janik Ditzel; janik.ditzel@cern.ch     ---

#ifndef ALIANALYSISTASKDOUBLEHYPNUCTREE_H
#define ALIANALYSISTASKDOUBLEHYPNUCTREE_H

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
#include "AliMCParticle.h"
#include "AliESDEvent.h"


class AliAnalysisTaskDoubleHypNucTree : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskDoubleHypNucTree();
  AliAnalysisTaskDoubleHypNucTree(const char *name);
  virtual ~AliAnalysisTaskDoubleHypNucTree();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(const Option_t*);
  void SelectPIDcheckOnly(Bool_t pidch = kFALSE) {fPIDCheckOnly = pidch;};
  void SetTriggerMask(UInt_t triggerMask = AliVEvent::kAny) {fTriggerMask = triggerMask;};

enum PdgCodeType_t {
    kPDGPionPlus,
    kPDGPionMinus,
    kPDGProton,
    kPDGAntiProton,    
    kPDGDeuteron,
    kPDGAntiDeuteron,
    kPDGTriton,
    kPDGAntiTriton,
    kPDGHelium3,
    kPDGAntiHelium3,
    kPDGHelium4,
    kPDGAntiHelium4,
    kPDGHyperHydrogen3,
    kPDGAntiHyperHydrogen3,
    kPDGHyperHydrogen4,
    kPDGAntiHyperHydrogen4,
    kPDGHyperHelium4,
    kPDGAntiHyperHelium4,
    kPDGHyperHelium5,
    kPDGAntiHyperHelium5,
    kPDGDoubleHyperHydrogen4,
    kPDGAntiDoubleHyperHydrogen4,
    kPdgCodeMax  // Number of enum entries                                                                                                                                          
  };

  static const Int_t fgkPdgCode[];

 private:
  //-- General variables --
  AliESDInputHandler    *fInputHandler;        
  AliESDpid             *fPID;                 
  AliESDEvent           *fESDevent;     
  AliMCEvent            *mcEvent;
  AliStack              *fStack;     
  AliEventCuts          fEventCuts; 
  AliESDtrackCuts       *trackCutsNuclei;
  AliESDtrackCuts       *trackCutsStrong;
  AliESDtrackCuts       *trackCutsSoft;
  TTree                 *fTree;
  TTree                 *gTree;
  TTree                 *fTreeGen;
  TTree                 *gTreeGen;
  TList                 *fHistogramList; 
  TH2F                  *fHistdEdx;            
  TH1F                  *fHistNumEvents;       
  TH1F			*fHistTrigger;	
  TH1F 	                *hInvMass4LHe;
  TH1F 	                *hInvMass4LLH;
  Double_t              fBetheParamsHe[6];  
  Double_t              fBetheParamsT[6];
  Bool_t                fPIDCheckOnly;  
  Int_t                 fMCtrue; 
  UInt_t		fTriggerMask; 
  Int_t                 fvariante; 
  //-- Trigger --
  Int_t                 MB, HMV0, HMSPD, HNU, HQU;
  //-- Tracks -- 
  AliESDtrack           *track1;
  AliESDtrack           *track2;
  AliESDtrack           *track3;
  AliESDtrack           *track4;
  Double_t              ptot1, sign1, ptot2, sign2, ptot3, sign3, ptot4, sign4;
  Bool_t                He3Pos1, He3Neg1, pPos2, pNeg2, piPos3, piNeg3, piPos4, piNeg4;
  //-- Full track arrays --
  std::vector < int >   He3PosArray;
  std::vector < int >   PPosArray;
  std::vector < int >   PiPosArray;
  std::vector < int >   PiPosSecArray;
  std::vector < int >   He3NegArray;
  std::vector < int >   PNegArray;
  std::vector < int >   PiNegArray;
  std::vector < int >   PiNegSecArray;
  //-- Reduced track arrays --
  std::vector < int >   He3PosArrayAcc;
  std::vector < int >   PPosArrayAcc;
  std::vector < int >   PiPosArrayAcc;
  std::vector < int >   He3NegArrayAcc;
  std::vector < int >   PNegArrayAcc;
  std::vector < int >   PiNegArrayAcc;
  //-- Array index counter --
  Int_t                 He3PosCounter;
  Int_t                 He3NegCounter;
  Int_t                 PPosCounter;
  Int_t                 PNegCounter;
  Int_t                 PiPosCounter;
  Int_t                 PiNegCounter;
  Int_t                 PiPosSecCounter;
  Int_t                 PiNegSecCounter;
  Int_t                 He3PosCounterAcc;
  Int_t                 He3NegCounterAcc;
  Int_t                 PPosCounterAcc;
  Int_t                 PNegCounterAcc;
  Int_t                 PiPosCounterAcc;
  Int_t                 PiNegCounterAcc;
  //-- Vertex Reco --
  AliESDVertex          *primVertex;
  const AliESDVertex    *vertex;
  AliESDVertex          *secVertex;
  AliESDVertex          *tertVertex;
  AliVertexerTracks     *secvertexer;
  AliVertexerTracks     *tertvertexer;   
  AliExternalTrackParam *exTrack;
  TObjArray             *trkArray;
  TObjArray             *trkArray1;
  Double_t              PrimVertex[3];
  Double_t              SecVertex[3];
  Double_t              TertVertex[3];
  Double_t              TertVertex2[3];
  Double_t              cov[21];
  Double_t              cov0[21];
  Double_t              cov1[21];
  Double_t              cov2[21];
  UShort_t              id[2];
  Double_t              pxpypz[3];
  Double_t              xyz[3];
  Double_t              dd[3];  
  Short_t               sign; 
  //-- MC --
  Int_t                 noCombinatoricsMC;
  Int_t                 lessCombinatoricsMC;
  Int_t                 stackN;
  Int_t                 label1;
  Int_t                 label2;
  Int_t                 label3;
  Int_t                 label4;   
  Int_t                 labelMother1;
  Int_t                 labelMother2;
  Int_t                 labelMother3;
  Int_t                 labelMother4;  
  Int_t                 labelGrandMother1;
  Int_t                 labelGrandMother2;
  Int_t                 labelGrandMother3;
  Int_t                 labelGrandMother4;
  Long_t                PDGCodeMother;
  AliMCParticle         *FirstDaughter;
  AliMCParticle         *SecondDaughter;
  AliMCParticle         *ThirdDaughter;
  AliMCParticle         *FourthDaughter;
  AliMCParticle         *ParticleMother;
  AliMCParticle         *ParticleMother1;
  AliMCParticle         *ParticleMother2;
  AliMCParticle         *ParticleMother3;
  AliMCParticle         *ParticleMother4;
  AliMCParticle         *ParticleGrandMother1;
  AliMCParticle         *ParticleGrandMother2;
  AliMCParticle         *ParticleGrandMother3;
  AliMCParticle         *ParticleGrandMother4;
  //-- Particle --
  TLorentzVector        *sublorentzsum;
  TLorentzVector        *sublorentzsum2;
  TLorentzVector        *particle1;
  TLorentzVector        *particle2;
  TLorentzVector        *particle3;
  TLorentzVector        *particle4;
  TLorentzVector        *lorentzsum;
  TLorentzVector        *lorentzsum2;
  TVector3              *h;
  Double_t              xthiss; 
  Double_t              xpp; 
  //-- saved in Tree --
  Int_t		        fPeriod; 
  Int_t                 fMagneticField;
  Int_t                 fCentrality;
  Int_t                 frunnumber;
  Int_t                 mctruth;
  Int_t                 feventclass;
  Int_t                 fTrigMB, fTrigHMV0, fTrigHMSPD, fTrigHNU, fTrigHQU;
  Int_t                 fPDGMother, fChargeMother;
  Int_t                 fNclsDaughter, fNclsITSDaughter, fNclsDaughter1, fNclsITSDaughter1, fNclsDaughter2, fNclsITSDaughter2, fNclsDaughter3, fNclsITSDaughter3;
  Int_t                 fPropDCADaughter, fPropDCADaughter1, fPropDCADaughter2, fPropDCADaughter3, fPropDCADaughter4;
  Float_t               fEDaughter, fpDaughter, fptDaughter,  fpxDaughter,  fpyDaughter,  fpzDaughter, fyDaughter, fdEdxDaughter, fdEdxSigmaDaughter, fDcaDaughter, fDcaDaughtero, fDcazDaughter, fDcaSecDaughter,  fChi2Daughter,  fEtaDaughter, fPhiDaughter, fGeoLengthDaughter, fTOFSignalDaughter;
  Float_t               fEDaughter1, fpDaughter1, fptDaughter1, fpxDaughter1, fpyDaughter1, fpzDaughter1, fyDaughter1, fdEdxDaughter1, fdEdxSigmaDaughter1, fDcaDaughter1, fDcaDaughter1o, fDcaSecDaughter1, fDcazDaughter1, fChi2Daughter1, fEtaDaughter1, fPhiDaughter1, fGeoLengthDaughter1, fTOFSignalDaughter1;
  Float_t               fEDaughter2, fpDaughter2, fptDaughter2, fpxDaughter2, fpyDaughter2, fpzDaughter2, fyDaughter2, fdEdxDaughter2, fdEdxSigmaDaughter2, fDcaDaughter2, fDcaDaughter2o, fDcaSecDaughter2, fDcazDaughter2, fChi2Daughter2, fEtaDaughter2, fPhiDaughter2, fGeoLengthDaughter2, fTOFSignalDaughter2;
  Float_t               fEDaughter3, fpDaughter3, fptDaughter3, fpxDaughter3, fpyDaughter3, fpzDaughter3, fyDaughter3, fdEdxDaughter3, fdEdxSigmaDaughter3, fDcaDaughter3, fDcaDaughter3o, fDcaSecDaughter3, fDcazDaughter3, fChi2Daughter3, fEtaDaughter3, fPhiDaughter3, fGeoLengthDaughter3, fTOFSignalDaughter3;
  Float_t               fEDaughter4, fpDaughter4, fptDaughter4, fpxDaughter4, fpyDaughter4, fpzDaughter4, fyDaughter4, fDcaDaughter4, fDcazDaughter4, fDcaSecDaughter4;
  Float_t               fSigmaYXDaughter, fSigmaXYZDaughter, fSigmaZDaughter, fSigmaYXDaughter1, fSigmaXYZDaughter1, fSigmaZDaughter1, fSigmaYXDaughter2, fSigmaXYZDaughter2, fSigmaZDaughter2, fSigmaYXDaughter3, fSigmaXYZDaughter3, fSigmaZDaughter3, fSigmaYXDaughter4, fSigmaXYZDaughter4, fSigmaZDaughter4, fPtUncertDaughter, fPtUncertDaughter1, fPtUncertDaughter2, fPtUncertDaughter3, fPtUncertDaughter4;
  Int_t                 fTPCRefitDaughter, fITSRefitDaughter, fTPCRefitDaughter1, fITSRefitDaughter1, fTPCRefitDaughter2, fITSRefitDaughter2, fTPCRefitDaughter3, fITSRefitDaughter3;
  Float_t               fImParDaughter, fImPar2Daughter, fImParDaughter1, fImPar2Daughter1, fImParDaughter2, fImPar2Daughter2, fImParDaughter3, fImPar2Daughter3, fImParDaughter4, fImPar2Daughter4, fImParzDaughter, fImParz2Daughter, fImParzDaughter1, fImParz2Daughter1, fImParzDaughter2, fImParz2Daughter2, fImParzDaughter3, fImParz2Daughter3, fImParzDaughter4, fImParz2Daughter4;
  Float_t               fdEdxSigmaPion, fdEdxSigmaDeuteron, fdEdxSigmaTriton, fdEdxSigmaAlpha;
  Float_t               fDCA2B, fDCA2Bo, fDCA3B1, fDCA3B2, fDCA3B3, fPA, fSubPA, fSubPA2, fSubPA3, fSubPA4, fDecAngle, farmalpha, farmpt;
  Float_t               fPrimVertexX, fPrimVertexY, fPrimVertexZ,fSecVertexX, fSecVertexY, fSecVertexZ, fTertVertexX, fTertVertexY, fTertVertexZ, fTertVertex2X, fTertVertex2Y, fTertVertex2Z;
  Float_t               fmMother, fmMother2, fEMother, fpxMother, fpyMother, fpzMother, fptMother, fpMother, fyMother, fctMother;
  Float_t               fmSubMother, fESubMother, fpxSubMother, fpySubMother, fpzSubMother, fptSubMother, fpSubMother, fySubMother, fctSubMother;
  Float_t               fthetaP, fthetaN;
  //-- Functions --
  //-- only PID check --
  void dEdxCheck();
  //-- init track arrays --
  void InitArrays();
  //-- main function for track finding --
  void FindTracks();
  //-- pos 4LHe reco --
  void CalcPosDaughterHypNuc();
  //-- neg 4LHe reco --
  void CalcNegDaughterHypNuc();
  //-- pos 4LLH reco --
  void CalcPosMotherHypNuc();
  //-- neg 4LLH reco
  void CalcNegMotherHypNuc();
  //-- function for the cascade vertex reco --
  void CreateSecVertex();
  //-- filling track information into tree --
  void SetDaughterInformation(TString Mother);
  //-- set Bethe_Bloch Params by runnumber --
  void SetBetheBlochParams(Int_t runNumber);
  //-- Main MC function --
  void MCGenerated();                                      
  //-- gen 4LLH --
  void MCFourBodyDecay (Int_t stackN, AliMCParticle *tparticleMother, Long_t PDGMother, Long_t PDGFirstDaughter, Long_t PDGSecondDaughter, Long_t PDGThirdDaughter, Long_t PDGFourthDaughter, Double_t massFirstDaughter, Double_t massSecondDaughter, Double_t massThirdDaughter, Double_t massFourthDaughter);
  //-- gen 4LHe --
  void MCThreeBodyDecay (Int_t stackN, AliMCParticle *tparticleMother, Long_t PDGMother, Long_t PDGFirstDaughter, Long_t PDGSecondDaughter, Long_t PDGThirdDaughter, Double_t massFirstDaughter, Double_t massSecondDaughter, Double_t massThirdDaughter); 
  //-- stack label calculation --
  Int_t GetLabel(Int_t labelFirstMother, Int_t particlePdgCode, Int_t motherparticlePdgCode);
  Int_t GetLabel(Int_t labelFirstMother, Int_t particlePdgCode);
  //-- trigger selection and simulation --
  Bool_t TriggerSelection();
  //-- see AliExternalTrackParams --
  Float_t GetInvPtDevFromBC(Int_t b, Int_t c);
  //-- clac of dEdx --
  Double_t Bethe(const AliESDtrack& track, Double_t mass, Int_t charge, Double_t* params);
  //
  Double_t GeoLength(const AliESDtrack& track);
  //
  Double_t GetTOFSignal(const AliESDtrack& track);
  //
  void ResetVals(TString mode);
  //
  AliAnalysisTaskDoubleHypNucTree(const AliAnalysisTaskDoubleHypNucTree&);
  AliAnalysisTaskDoubleHypNucTree &operator=(const AliAnalysisTaskDoubleHypNucTree&);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskDoubleHypNucTree, 3);
  /// \endcond
};
#endif
