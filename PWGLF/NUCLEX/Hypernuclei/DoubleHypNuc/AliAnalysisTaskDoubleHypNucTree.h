//--- Task for investigation of the DoubleHyperHydrogen4 ---
//---     Author: Janik Ditzel; janik.ditzel@cern.ch     ---

#ifndef ALIANALYSISTASKDOUBLEHYPNUCTREE_H
#define ALIANALYSISTASKDOUBLEHYPNUCTREE_H


#define HomogeneousField ///homogenous field in z direction
#include <algorithm>
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliInputEventHandler.h"
#include "AliAODTrack.h"
#include "AliPIDResponse.h"
#include "AliAODVertex.h"
#include "AliAODEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliEventCuts.h"
#include "AliKalmanTrack.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliMCEvent.h"
#include "AliMCVertex.h"
#include "AliPID.h"
#include "AliStack.h"
#include "AliVertexerTracks.h"
#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include <fstream>
#include "KFParticle.h"
#include "KFVertex.h"
#include "KFPVertex.h"
#include "KFPTrack.h"
#include <vector>
#include "TCanvas.h"
#include "TChain.h"
#include <TClonesArray.h>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TObject.h"
#include "TPDGCode.h"
#include "TLorentzVector.h"
#include "TRandom2.h"
#include "TTree.h"
#include "TVector3.h"

class TH1F;
class TH2F;
class AliESDEvent;
class AliESDpid;
class AliESDVertex;
class AliESDTrack;
class AliAODVertex;
class AliAODTrack;
class AliAODEvent;
class AliVEvent;
class KFParticle;
class KFVertex;
class AliKalmanTrack;
class AliITStrackerMI;

using namespace std;
// _________________________________________________ //
class KFParticleHypNuc : public KFParticle {
public:
  Bool_t CheckDaughter(KFParticle daughter) {
    Float_t m[8], mV[36], D[3][3];
    if (KFParticleBase::GetMeasurement(daughter, m, mV, D)) return kTRUE;
    return kFALSE;
  }
};
// _________________________________________________ //
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
        kPDGHyperHydrogen4Star,
        kPDGAntiHyperHydrogen4Star,
        kPDGHyperHelium4Star,
        kPDGAntiHyperHelium4Star,
        kPDGHyperHelium5,
        kPDGAntiHyperHelium5,
        kPDGDoubleHyperHydrogen4,
        kPDGAntiDoubleHyperHydrogen4,
        kPdgCodeMax  // Number of enum entries
    };

  static const Int_t fgkPdgCode[];

private:
  //-- General variables --
  AliInputEventHandler    *fInputHandler;
  AliPIDResponse        *fPID;
  AliESDEvent           *fESDevent;
  AliAODEvent           *fAODevent;
  AliVHeader            *eventHeader;
  AliVEvent             *fevent;
  AliMCEvent            *mcEvent;
  AliStack              *fStack;     
  AliEventCuts          fEventCuts; 
  TTree                 *fTree;
  TTree                 *fTreeKF;
  TTree                 *gTree;
  TTree                 *gTreeKF;
  TTree                 *hTree;
  TTree                 *hTreeKF;
  TTree                 *iTree;
  TTree                 *iTreeKF;
  TTree                 *fTreeGen;
  TTree                 *gTreeGen;
  TList                 *fHistogramList; 
  TH2F                  *fHistdEdx;            
  TH1F                  *fHistNumEvents;       
  TH1F			*fHistTrigger;
  TH1F                  *fHistTrigger1;
  TH1F                  *fHistTrigger2;
  TH1F                  *fHistCentrality1;
  TH1F                  *fHistCentrality2;
  TH1F                  *fHistNtracks1;
  TH1F                  *fHistNtracks2;
  Double_t              fBetheParamsHe[6];  
  Double_t              fBetheParamsT[6];
  Bool_t                fPIDCheckOnly;
  Int_t                 fMCtrue; 
  UInt_t		        fTriggerMask; 
  Int_t                 kVariante; 
  Int_t                 kStandardReco;
  Int_t                 kKFReco;
  Int_t                 kSkip3LHAnalysis;
  Int_t                 kSkip4LHAnalysis;
  Int_t                 kAOD, kESD;
  //-- Trigger --
  Int_t                 MB, HMV0, HMSPD, HNU, HQU, Central, SemiCentral;
  //-- Tracks -- 
  AliAODTrack           *track1;
  AliAODTrack           *track2;
  AliAODTrack           *track3;
  AliAODTrack           *track4;
  //-- Full track arrays --
  std::vector < int >   He4PosArray;
  std::vector < int >   He3PosArray;
  std::vector < int >   PPosArray;
  std::vector < int >   PiPosArray;
  std::vector < int >   PiPosSecArray;
  std::vector < int >   He4NegArray;
  std::vector < int >   He3NegArray;
  std::vector < int >   PNegArray;
  std::vector < int >   PiNegArray;
  std::vector < int >   PiNegSecArray;
  //-- Array index counter --
  Int_t                 He4PosCounter;
  Int_t                 He3PosCounter;
  Int_t                 He4NegCounter;
  Int_t                 He3NegCounter;
  Int_t                 PPosCounter;
  Int_t                 PNegCounter;
  Int_t                 PiPosCounter;
  Int_t                 PiNegCounter;
  Int_t                 PiPosSecCounter;
  Int_t                 PiNegSecCounter;
  //-- Vertex Reco --
  AliESDVertex          *primVertex;
  KFVertex              primKFVertex;
  const AliAODVertex    *vertex;
  AliESDVertex          *secVertex;
  AliESDVertex          *tertVertex;
  AliESDVertex          *KFtertVertex;
  AliESDVertex          *KFsecVertex;
  AliVertexerTracks     *secvertexer;
  AliVertexerTracks     *tertvertexer;   
  AliExternalTrackParam *exTrack;
  AliExternalTrackParam *exTrack1;
  AliExternalTrackParam *exTrack2;
  AliExternalTrackParam *exTrack3;
  AliExternalTrackParam *exTrack4;
  TObjArray             *trkArray;
  TObjArray             *trkArray1;
  Double_t              PrimVertex[3];
  Double_t              SecVertex[3];
  Double_t              TertVertex[3];
  Double_t              PrimVertexKF[3];
  Double_t              SecVertexKF[3];
  Double_t              SecVertexErrKF[3];
  Double_t              TertVertexKF[3];
  Double_t              TertVertexErrKF[3];
  Double_t              cov[21];
  Double_t              cov0[21];
  Double_t              cov1[21];
  Double_t              cov2[21];
  Double_t              pxpypz[3];
  Double_t              xyz[3];
  Double_t              dd[3];
  Double_t              kMagF;
  Short_t               sign;
  KFParticle            KFtrack1;
  KFParticle            KFtrack2;
  KFParticle            KFtrack3;
  KFParticle            KFtrack4;
  //-- Cuts --
  Double_t              kDCATracksCut;
  Double_t              kTrackPtUncertCut;
  Double_t              kPointingAngleCut;
  Double_t              kPIDSelecCut;
  Double_t              kMin4LLHMass;
  Double_t              kMax4LLHMass;
  Double_t              kMin4LHeMass;
  Double_t              kMax4LHeMass;
  Double_t              kMin4LHMass;
  Double_t              kMax4LHMass;
  Double_t              kMin3LHMass;
  Double_t              kMax3LHMass;
  //-- MC --
  Bool_t                konlyBG;  
  Bool_t                konlySig;
  Bool_t                kMCPIDCheck;
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
  AliMCParticle         *CheckParticle;
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
  Double_t               xthiss; 
  Double_t               xpp; 
  //-- saved in Tree --
  Int_t                  fEventID;
  Int_t                  fParticleID;
  Int_t		         fPeriod; 
  Int_t                  fMagneticField;
  Int_t                  fCentrality;
  Int_t                  frunnumber;
  Int_t                  fmctruth;
  Int_t                  feventclass;
  Int_t                  fPrimary4LHe;
  Int_t                  fisExcited;
  Int_t                  fDecayChannel;
  Int_t                  fRecoMethod;
  Int_t                  fisWOTrackVertex;
  const Char_t*          fTriggerString;
  Int_t                  fisOnlineV0_13;
  Int_t                  fisOnlineV0_14;
  Int_t                  fisOnlineV0_23;
  Int_t                  fisOnlineV0_24;
  Double_t               fptDaughterUnProp, fpxDaughterUnProp, fpyDaughterUnProp, fpzDaughterUnProp;
  Double_t               fptDaughter1UnProp, fpxDaughter1UnProp, fpyDaughter1UnProp, fpzDaughter1UnProp;
  Double_t               fptDaughter2UnProp, fpxDaughter2UnProp, fpyDaughter2UnProp, fpzDaughter2UnProp;
  Double_t               fptDaughter3UnProp, fpxDaughter3UnProp, fpyDaughter3UnProp, fpzDaughter3UnProp;
  Int_t                  fTrigMB, fTrigHMV0, fTrigHMSPD, fTrigHNU, fTrigHQU, fTrigkCentral, fTrigkSemiCentral;
  Int_t                  fPDGMother, fChargeMother;
  Int_t                  fPrimVertNDF, fSecVertNDF, fTertVertNDF;  
  Int_t                  fNclsDaughter, fNclsITSDaughter, fNclsDaughter1, fNclsITSDaughter1, fNclsDaughter2, fNclsITSDaughter2, fNclsDaughter3, fNclsITSDaughter3;
  Int_t                  fPropDCADaughter, fPropDCADaughter1, fPropDCADaughter2, fPropDCADaughter3, fPropDCADaughter4;
  Int_t                  fTPCRefitDaughter, fITSRefitDaughter, fTPCRefitDaughter1, fITSRefitDaughter1, fTPCRefitDaughter2, fITSRefitDaughter2, fTPCRefitDaughter3, fITSRefitDaughter3;
  Int_t                  fITSLayer1Daughter, fITSLayer2Daughter, fITSLayer3Daughter, fITSLayer4Daughter, fITSLayer5Daughter, fITSLayer6Daughter, fITSLayer1Daughter1, fITSLayer2Daughter1, fITSLayer3Daughter1, fITSLayer4Daughter1, fITSLayer5Daughter1, fITSLayer6Daughter1, fITSLayer1Daughter2, fITSLayer2Daughter2, fITSLayer3Daughter2, fITSLayer4Daughter2, fITSLayer5Daughter2, fITSLayer6Daughter2, fITSLayer1Daughter3, fITSLayer2Daughter3, fITSLayer3Daughter3, fITSLayer4Daughter3, fITSLayer5Daughter3, fITSLayer6Daughter3;
  Int_t                  fTrackPIDDaughter, fTrackPIDDaughter1, fTrackPIDDaughter2, fTrackPIDDaughter3;
  Double_t               fEDaughter, fpDaughter, fptDaughter,  fpxDaughter,  fpyDaughter,  fpzDaughter, fyDaughter, fdEdxDaughter, fdEdxSigmaDaughter, fDcaDaughter, fDcaDaughtero, fDcazDaughter, fDcaSecDaughter,  fChi2Daughter,  fEtaDaughter, fPhiDaughter, fGeoLengthDaughter, fTOFSignalDaughter;
  Double_t               fEDaughter1, fpDaughter1, fptDaughter1, fpxDaughter1, fpyDaughter1, fpzDaughter1, fyDaughter1, fdEdxDaughter1, fdEdxSigmaDaughter1, fDcaDaughter1, fDcaDaughter1o, fDcaSecDaughter1, fDcazDaughter1, fChi2Daughter1, fEtaDaughter1, fPhiDaughter1, fGeoLengthDaughter1, fTOFSignalDaughter1;
  Double_t               fEDaughter2, fpDaughter2, fptDaughter2, fpxDaughter2, fpyDaughter2, fpzDaughter2, fyDaughter2, fdEdxDaughter2, fdEdxSigmaDaughter2, fDcaDaughter2, fDcaDaughter2o, fDcaSecDaughter2, fDcazDaughter2, fChi2Daughter2, fEtaDaughter2, fPhiDaughter2, fGeoLengthDaughter2, fTOFSignalDaughter2;
  Double_t               fEDaughter3, fpDaughter3, fptDaughter3, fpxDaughter3, fpyDaughter3, fpzDaughter3, fyDaughter3, fdEdxDaughter3, fdEdxSigmaDaughter3, fDcaDaughter3, fDcaDaughter3o, fDcaSecDaughter3, fDcazDaughter3, fChi2Daughter3, fEtaDaughter3, fPhiDaughter3, fGeoLengthDaughter3, fTOFSignalDaughter3;
  Double_t               fEDaughter4, fpDaughter4, fptDaughter4, fpxDaughter4, fpyDaughter4, fpzDaughter4, fyDaughter4, fDcaDaughter4, fDcazDaughter4, fDcaSecDaughter4;
  Double_t               fSigmaYXDaughter, fSigmaXYZDaughter, fSigmaZDaughter, fSigmaYXDaughter1, fSigmaXYZDaughter1, fSigmaZDaughter1, fSigmaYXDaughter2, fSigmaXYZDaughter2, fSigmaZDaughter2, fSigmaYXDaughter3, fSigmaXYZDaughter3, fSigmaZDaughter3, fSigmaYXDaughter4, fSigmaXYZDaughter4, fSigmaZDaughter4, fPtUncertDaughter, fPtUncertDaughter1, fPtUncertDaughter2, fPtUncertDaughter3, fPtUncertDaughter4;  
  Double_t               fImParDaughter,  fImParDaughter1,  fImParDaughter2,  fImParDaughter3,  fImParDaughter4,  fImParzDaughter,  fImParzDaughter1,  fImParzDaughter2,  fImParzDaughter3,  fImParzDaughter4;
  Double_t               fdEdxSigmaPion, fdEdxSigmaDeuteron, fdEdxSigmaTriton, fdEdxSigmaAlpha;
  //
  Double_t               fDCA2B, fDCA2Bo, fDCA3B1, fDCA3B2, fDCA3B3, fPA, fSubPA, fSubPA2, fDecAngle, farmalpha, farmpt;
  Double_t               fPrimVertexX, fPrimVertexY, fPrimVertexZ,fSecVertexX, fSecVertexY, fSecVertexZ, fTertVertexX, fTertVertexY, fTertVertexZ, fV0VertexX_13, fV0VertexY_13, fV0VertexZ_13, fV0VertexX_14, fV0VertexY_14, fV0VertexZ_14, fV0VertexX_23, fV0VertexY_23, fV0VertexZ_23, fV0VertexX_24, fV0VertexY_24, fV0VertexZ_24;
  Double_t               fPrimVertChi2, fSecVertChi2, fTertVertChi2;  
  Double_t               fmMother, fmMother2, fEMother, fpxMother, fpyMother, fpzMother, fptMother, fpMother, fyMother, fctMother;
  Double_t               fmSubMother, fESubMother, fpxSubMother, fpySubMother, fpzSubMother, fptSubMother, fpSubMother, fySubMother, fctSubMother;
  //-- KF --
  Int_t                  fPrimVertNDFKF, fSecVertNDFKF, fTertVertNDFKF;
  Double_t               fEDaughter1KF, fLabelDaughter1KF, fptDaughter1KF,  fpxDaughter1KF,  fpyDaughter1KF,  fpzDaughter1KF, fyDaughter1KF;
  Double_t               fEDaughter2KF, fLabelDaughter2KF, fptDaughter2KF,  fpxDaughter2KF,  fpyDaughter2KF,  fpzDaughter2KF, fyDaughter2KF;
  Double_t               fEDaughter3KF, fLabelDaughter3KF, fptDaughter3KF,  fpxDaughter3KF,  fpyDaughter3KF,  fpzDaughter3KF, fyDaughter3KF;
  Double_t               fEDaughterKF, fLabelDaughterKF, fptDaughterKF,  fpxDaughterKF,  fpyDaughterKF,  fpzDaughterKF, fyDaughterKF;
  Double_t               fDCA2BXYKF, fDCA3B1XYKF, fDCA3B2XYKF, fDCA3B3XYKF, fPAKF, fSubPAKF, fSubPA2KF;
  Double_t               fDCA2BZKF, fDCA3B1ZKF, fDCA3B2ZKF, fDCA3B3ZKF;
  Double_t               fPrimVertexXKF, fPrimVertexYKF, fPrimVertexZKF, fSecVertexXKF, fSecVertexYKF, fSecVertexZKF, fTertVertexXKF, fTertVertexYKF, fTertVertexZKF;
  Double_t               fPrimVertChi2KF, fSecVertChi2KF, fTertVertChi2KF;  
  Double_t               fmMotherKF, fEMotherKF, fpxMotherKF, fpyMotherKF, fpzMotherKF, fptMotherKF, fpMotherKF, fyMotherKF, fctMotherKF;
  Double_t               fmSubMotherKF, fESubMotherKF, fpxSubMotherKF, fpySubMotherKF, fpzSubMotherKF, fptSubMotherKF, fpSubMotherKF, fySubMotherKF, fctSubMotherKF;
  Double_t               fDcaDaughterXYKF, fDcaDaughter1XYKF, fDcaDaughter2XYKF, fDcaDaughter3XYKF, fDcaDaughter4XYKF, fDcaSecDaughterXYKF, fDcaSecDaughter1XYKF, fDcaSecDaughter2XYKF, fDcaSecDaughter3XYKF, fDcaSecDaughter4XYKF;
  Double_t               fDcaDaughterZKF, fDcaDaughter1ZKF, fDcaDaughter2ZKF, fDcaDaughter3ZKF, fDcaDaughter4ZKF, fDcaSecDaughterZKF, fDcaSecDaughter1ZKF, fDcaSecDaughter2ZKF, fDcaSecDaughter3ZKF, fDcaSecDaughter4ZKF;
  Double_t               fPrimVertexXErrKF, fPrimVertexYErrKF, fPrimVertexZErrKF, fSecVertexXErrKF, fSecVertexYErrKF, fSecVertexZErrKF, fTertVertexXErrKF, fTertVertexYErrKF, fTertVertexZErrKF;
  Double_t               fmMotherErrKF, fEMotherErrKF, fpxMotherErrKF, fpyMotherErrKF, fpzMotherErrKF, fptMotherErrKF, fpMotherErrKF, fctMotherErrKF;
  Double_t               fmSubMotherErrKF, fESubMotherErrKF, fpxSubMotherErrKF, fpySubMotherErrKF, fpzSubMotherErrKF, fptSubMotherErrKF, fpSubMotherErrKF, fctSubMotherErrKF;
  Double_t               fmSubMotherCheck, fmMotherCheck;
  Double_t               fCovMatrixTrack[21], fCovMatrixTrack1[21], fCovMatrixTrack2[21], fCovMatrixTrack3[21], fCovMatrixTrack4[21];
  Double_t               fTrackPar[7], fTrackPar1[7], fTrackPar2[7], fTrackPar3[7], fTrackPar4[7];
  //
  Double_t               fthetaP, fthetaN;  
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
  void CalcPosMotherHypNuc(Int_t Track1Entry, Int_t Track2Entry, Int_t Track3Entry, Int_t kDecayChannel);
  //-- neg 4LLH reco
  void CalcNegMotherHypNuc(Int_t Track1Entry, Int_t Track2Entry, Int_t Track3Entry, Int_t kDecayChannel);
  //-- function for the vertex reco --
  Int_t StandardReconstruction(TString Mother, Int_t kDecayChannel, Int_t ksign);
  //-- function for the  vertex reco with KFParticle--
  Int_t KFReconstruction(TString Mother, Int_t kDecayChannel, Int_t ksign);
  //-- function for the CORRECT initialization of a KFParticle--
  KFParticle CreateKFParticle(AliExternalTrackParam& track, Double_t Mass, Int_t Charge); 
  //-- function for the CORRECT initialization of a KFVertex--
  KFVertex CreateKFVertex(const AliVVertex& vertex);
  //-- check whether two used tracks are from a V0 --
  void GetV0Status();
  //-- filling track information into tree --
  void SetDaughterInformation(TString Mother, Int_t kDecayChannel);
  void SetDaughterInformationKF(KFParticle daughter1, KFParticle daughter2, KFParticle daughter3, Int_t kDecayChannel);
  void SetDaughterInformationKF(KFParticle daughter1, KFParticle daughter2, KFParticle daughter3, KFParticle daughter4, Int_t kDecayChannel);
  //-- set Bethe_Bloch Params by runnumber --
  void SetBetheBlochParams(Int_t runNumber);
  //-- Main MC function --
  void MCGenerated();                                      
  //-- gen 4LLH --
  void MCFourBodyDecay (Int_t stackN, AliMCParticle *tparticleMother, Long_t PDGMother, Int_t kDecayChannel, Long_t PDGFirstDaughter, Long_t PDGSecondDaughter, Long_t PDGThirdDaughter, Long_t PDGFourthDaughter, Double_t massFirstDaughter, Double_t massSecondDaughter, Double_t massThirdDaughter, Double_t massFourthDaughter);
  //-- gen 4LHe --
  void MCThreeBodyDecay (Int_t stackN, AliMCParticle *tparticleMother, Long_t PDGMother, Long_t PDGFirstDaughter, Long_t PDGSecondDaughter, Long_t PDGThirdDaughter, Double_t massFirstDaughter, Double_t massSecondDaughter, Double_t massThirdDaughter);
  //--gen 3LH
  void MCTwoBodyDecay(Int_t stackN, AliMCParticle* ParticleMother, Long_t PDGMother, Long_t PDGFirstDaughter, Long_t PDGSecondDaughter, Double_t massFirstDaughter, Double_t massSecondDaughter);
  //-- stack label calculation --
  Int_t GetLabel(Int_t labelFirstMother, Int_t particlePdgCode, Int_t motherparticlePdgCode);
  Int_t GetLabel(Int_t labelFirstMother, Int_t particlePdgCode);
  //Convert AOD to ESD Vertex
  AliESDVertex* AODToESDVertex(const AliAODVertex &aodVert);
  //-- trigger selection and simulation -
  Bool_t TriggerSelection();
  //-- see AliExternalTrackParams --
  Float_t GetInvPtDevFromBC(Int_t b, Int_t c);
  //-- clac of dEdx --
  Double_t Bethe(const AliAODTrack& track, Double_t mass, Int_t charge, Double_t* params);
  //
  Int_t CustomTrackCut(const AliAODTrack& track, Int_t particle);  
  //
  Double_t GeoLength(const AliAODTrack& track);
  //
  Double_t GetLengthInActiveZone(const AliExternalTrackParam  *paramT, Double_t deltaY, Double_t deltaZ, Double_t bz, Double_t exbPhi);
  //
  Double_t GetTOFSignal(const AliAODTrack& track);  
  //
  void ResetVals(TString mode);
  //
  AliAnalysisTaskDoubleHypNucTree(const AliAnalysisTaskDoubleHypNucTree&);
  AliAnalysisTaskDoubleHypNucTree &operator=(const AliAnalysisTaskDoubleHypNucTree&);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskDoubleHypNucTree, 1);
  /// \endcond
};
#endif
