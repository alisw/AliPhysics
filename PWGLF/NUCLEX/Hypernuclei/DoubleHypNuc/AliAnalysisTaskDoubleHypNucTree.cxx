//--- Task for investigation of the {}^{4}_{#Lambda#Lambda}H ---
//---     Author: Janik Ditzel; janik.ditzel@cern.ch     ---

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom2.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDv0.h"
#include "AliESDVertex.h"
#include "AliVertexerTracks.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCVertex.h"
#include "AliStack.h"
#include "TPDGCode.h"
#include "AliEventCuts.h"
#include "AliAnalysisTaskDoubleHypNucTree.h"
#include "TLorentzVector.h"
#include <TClonesArray.h>
#include "TObject.h"

using namespace std;

ClassImp(AliAnalysisTaskDoubleHypNucTree)
// _________________________________________________ //
// Default Constructor
AliAnalysisTaskDoubleHypNucTree::AliAnalysisTaskDoubleHypNucTree()
:AliAnalysisTaskSE("AliAnalysisTaskDoubleHypNucTree"),
  fPIDCheckOnly(kFALSE),
  fInputHandler(0),
  fPID(0),
  fESDevent(0),
  mcEvent(0),
  fStack(),
  fEventCuts(),  
  fTriggerMask(),
  MB(0),
  HMV0(0),
  HMSPD(0),
  HNU(0),
  HQU(0),
  trackCutsNuclei(0),
  trackCutsStrong(0),
  trackCutsSoft(0),
  fTree(0),
  gTree(0),
  fTreeGen(0),
  gTreeGen(0), 
  fvariante(1),
  He3PosArray(),
  PPosArray(),
  PiPosArray(),
  PiPosSecArray(),
  He3NegArray(),
  PNegArray(),
  PiNegArray(),
  PiNegSecArray(),
  He3PosCounter(),
  He3NegCounter(),
  PPosCounter(),
  PNegCounter(),
  PiPosCounter(),
  PiNegCounter(),
  PiPosSecCounter(),
  PiNegSecCounter(),
  He3PosArrayAcc(),
  PPosArrayAcc(),
  PiPosArrayAcc(),
  He3NegArrayAcc(),
  PNegArrayAcc(),
  PiNegArrayAcc(),
  He3PosCounterAcc(),
  He3NegCounterAcc(),
  PPosCounterAcc(),
  PNegCounterAcc(),
  PiPosCounterAcc(),
  PiNegCounterAcc(),
  fBetheParamsHe(),
  fBetheParamsT(),  
  fMCtrue(0),
  vertex(),
  primVertex(),
  secvertexer(),
  secVertex(),
  tertvertexer(),
  tertVertex(),
  PrimVertex(),
  SecVertex(),
  TertVertex(),
  TertVertex2(),
  track1(),
  track2(),
  track3(),
  track4(),
  exTrack(),
  trkArray(),
  trkArray1(),
  ptot1(0),
  sign1(0),
  ptot2(0),
  sign2(0),
  ptot3(0),
  sign3(0),
  ptot4(0),
  sign4(0),
  He3Pos1(0),
  He3Neg1(0),
  pPos2(0),
  pNeg2(0),
  piPos3(0),
  piNeg3(0),
  piPos4(0),
  piNeg4(0),    
  cov(),
  cov0(),
  cov1(),
  cov2(),
  id(),
  pxpypz(),
  xyz(),
  sign(0),
  dd(),
  xthiss(),
  xpp(),
  h(),
  fthetaP(-99), 
  fthetaN(-99),
//
  noCombinatoricsMC(0),
  lessCombinatoricsMC(0),
  stackN(0),
  PDGCodeMother(),
  ParticleMother(),
  FirstDaughter(),
  SecondDaughter(),
  ThirdDaughter(),
  FourthDaughter(),
  label1(),
  labelMother1(),
  labelGrandMother1(),
  ParticleMother1(),
  ParticleGrandMother1(),
  label2(),
  labelMother2(),
  labelGrandMother2(),
  ParticleMother2(),
  ParticleGrandMother2(),
  label3(),
  labelMother3(),
  labelGrandMother3(),
  ParticleMother3(),
  ParticleGrandMother3(),
  label4(),
  labelMother4(),
  labelGrandMother4(),
  ParticleMother4(),
  ParticleGrandMother4(),
//
  lorentzsum(NULL),
  lorentzsum2(NULL),
  sublorentzsum(NULL),
  sublorentzsum2(NULL),
  particle1(NULL),
  particle2(NULL),
  particle3(NULL),
  particle4(NULL),
//Histos
  fHistogramList(NULL),
  fHistdEdx(0),
  fHistNumEvents(0),
  fHistTrigger(0),
  hInvMass4LHe(0),
  hInvMass4LLH(0),
//saved in tree
  fCentrality(-99),
  frunnumber(-99),
  fPeriod(0),
  fMagneticField(0),
  feventclass(0),
  mctruth(0),
  fTrigMB(-99),			
  fTrigHMV0(-99),
  fTrigHMSPD(-99),
  fTrigHNU(0),
  fTrigHQU(0),
  fDCA2B(-99),
  fDCA2Bo(-99),
  fDCA3B1(-99),
  fDCA3B2(-99),
  fDCA3B3(-99),
  fDecAngle(-99),
  farmalpha(-99),
  farmpt(-99),
  fPA(-99),
  fSubPA(-99),
  fSubPA2(-99),
  fSubPA3(-99),
  fSubPA4(-99),
  fPrimVertexX(-99),
  fPrimVertexY(-99),
  fPrimVertexZ(-99),
  fSecVertexX(-99),
  fSecVertexY(-99),
  fSecVertexZ(-99),
  fTertVertexX(-99),
  fTertVertexY(-99),
  fTertVertexZ(-99),
  fTertVertex2X(-99),
  fTertVertex2Y(-99),
  fTertVertex2Z(-99),
//
  fPDGMother(-99),
  fChargeMother(-99),
  fmMother(-99),
  fmMother2(-99),
  fEMother(-99),
  fpxMother(-99),
  fpyMother(-99),
  fpzMother(-99),
  fptMother(-99),
  fpMother(-99),
  fyMother(-99),
  fctMother(-99),
//
  fmSubMother(-99),
  fESubMother(-99),
  fpxSubMother(-99),
  fpySubMother(-99),
  fpzSubMother(-99),
  fptSubMother(-99),
  fpSubMother(-99),
  fySubMother(-99),
  fctSubMother(-99),
//
  fEDaughter(-99),
  fpDaughter(-99),
  fptDaughter(-99),
  fpxDaughter(-99),
  fpyDaughter(-99),
  fpzDaughter(-99),
  fyDaughter(-99),
  fdEdxDaughter(-99),
  fdEdxSigmaDaughter(-99),
  fDcaDaughter(-99),
  fDcaDaughtero(-99),
  fDcazDaughter(-99),
  fDcaSecDaughter(-99),
  fImParDaughter(-99),      
  fImPar2Daughter(-99),
  fImParzDaughter(-99),
  fImParz2Daughter(-99),
  fNclsDaughter(-99),
  fChi2Daughter(-99),
  fNclsITSDaughter(-99),
  fEtaDaughter(-99),
  fPhiDaughter(-99),
  fGeoLengthDaughter(-99),
  fTOFSignalDaughter(-99),
  fSigmaYXDaughter(-99),
  fSigmaXYZDaughter(-99),
  fSigmaZDaughter(-99),
  fPtUncertDaughter(-99),
  fTPCRefitDaughter(-99),
  fITSRefitDaughter(-99),
  fPropDCADaughter(-1),
//
  fEDaughter1(-99),
  fpDaughter1(-99),
  fptDaughter1(-99),
  fpxDaughter1(-99),
  fpyDaughter1(-99),
  fpzDaughter1(-99),
  fyDaughter1(-99),
  fdEdxDaughter1(-99),
  fdEdxSigmaDaughter1(-99),
  fDcaDaughter1(-99),
  fDcaDaughter1o(-99),
  fDcazDaughter1(-99),
  fDcaSecDaughter1(-99),
  fImParDaughter1(-99),
  fImPar2Daughter1(-99),
  fImParzDaughter1(-99),
  fImParz2Daughter1(-99),
  fNclsDaughter1(-99),
  fChi2Daughter1(-99),
  fNclsITSDaughter1(-99),
  fEtaDaughter1(-99),
  fPhiDaughter1(-99),
  fGeoLengthDaughter1(-99),
  fTOFSignalDaughter1(-99),
  fSigmaYXDaughter1(-99),
  fSigmaXYZDaughter1(-99),
  fSigmaZDaughter1(-99),
  fPtUncertDaughter1(-99),
  fTPCRefitDaughter1(-99),
  fITSRefitDaughter1(-99),
  fPropDCADaughter1(-1),
//
  fEDaughter2(-99),
  fpDaughter2(-99),
  fptDaughter2(-99),
  fpxDaughter2(-99),
  fpyDaughter2(-99),
  fpzDaughter2(-99),
  fyDaughter2(-99),
  fdEdxDaughter2(-99),
  fdEdxSigmaDaughter2(-99),
  fDcaDaughter2(-99),
  fDcaDaughter2o(-99),
  fDcazDaughter2(-99),
  fDcaSecDaughter2(-99),
  fImParDaughter2(-99),
  fImPar2Daughter2(-99),
  fImParzDaughter2(-99),
  fImParz2Daughter2(-99),
  fNclsDaughter2(-99),
  fChi2Daughter2(-99),
  fNclsITSDaughter2(-99),
  fEtaDaughter2(-99),
  fPhiDaughter2(-99),
  fGeoLengthDaughter2(-99),
  fTOFSignalDaughter2(-99),
  fSigmaYXDaughter2(-99),
  fSigmaXYZDaughter2(-99),
  fSigmaZDaughter2(-99),
  fPtUncertDaughter2(-99),
  fTPCRefitDaughter2(-99),
  fITSRefitDaughter2(-99),
  fPropDCADaughter2(-1),
//
  fEDaughter3(-99),
  fpDaughter3(-99),
  fptDaughter3(-99),
  fpxDaughter3(-99),
  fpyDaughter3(-99),
  fpzDaughter3(-99),
  fyDaughter3(-99),
  fdEdxDaughter3(-99),
  fdEdxSigmaDaughter3(-99),
  fDcaDaughter3(-99),
  fDcaDaughter3o(-99),
  fDcazDaughter3(-99),
  fDcaSecDaughter3(-99),
  fImParDaughter3(-99),
  fImPar2Daughter3(-99),
  fImParzDaughter3(-99),
  fImParz2Daughter3(-99),
  fNclsDaughter3(-99),
  fChi2Daughter3(-99),
  fNclsITSDaughter3(-99),
  fEtaDaughter3(-99),
  fPhiDaughter3(-99),
  fGeoLengthDaughter3(-99),
  fTOFSignalDaughter3(-99),
  fSigmaYXDaughter3(-99),
  fSigmaXYZDaughter3(-99),
  fSigmaZDaughter3(-99),
  fPtUncertDaughter3(-99),
  fTPCRefitDaughter3(-99),
  fITSRefitDaughter3(-99),
  fPropDCADaughter3(-1),
//
  fEDaughter4(-99),
  fpDaughter4(-99),
  fptDaughter4(-99),
  fpxDaughter4(-99),
  fpyDaughter4(-99),
  fpzDaughter4(-99),
  fyDaughter4(-99),
  fDcaDaughter4(-99),
  fDcazDaughter4(-99),
  fDcaSecDaughter4(-99),
  fImParDaughter4(-99),
  fImPar2Daughter4(-99),
  fImParzDaughter4(-99),
  fImParz2Daughter4(-99),
  fSigmaYXDaughter4(-99),
  fSigmaXYZDaughter4(-99),
  fSigmaZDaughter4(-99),
  fPtUncertDaughter4(-99),
  fPropDCADaughter4(-1),
  fdEdxSigmaPion(-99),
  fdEdxSigmaDeuteron(-99),
  fdEdxSigmaTriton(-99),
  fdEdxSigmaAlpha(-99)
{

}
// _________________________________________________ //
// Constructor
AliAnalysisTaskDoubleHypNucTree::AliAnalysisTaskDoubleHypNucTree(const char *name)
  :AliAnalysisTaskSE(name),
   fPIDCheckOnly(kFALSE),
   fInputHandler(0),
   fPID(0),
   fESDevent(0),
   mcEvent(0),
   fStack(),
   fEventCuts(),  
   fTriggerMask(),
   MB(0),
   HMV0(0),
   HMSPD(0),
   HNU(0),
   HQU(0),
   trackCutsNuclei(0),
   trackCutsStrong(0),
   trackCutsSoft(0),
   fTree(0),
   gTree(0),
   fTreeGen(0),
   gTreeGen(0), 
   fvariante(1),
   He3PosArray(),
   PPosArray(),
   PiPosArray(),
   PiPosSecArray(),
   He3NegArray(),
   PNegArray(),
   PiNegArray(),
   PiNegSecArray(),
   He3PosCounter(),
   He3NegCounter(),
   PPosCounter(),
   PNegCounter(),
   PiPosCounter(),
   PiNegCounter(),
   PiPosSecCounter(),
   PiNegSecCounter(),
   He3PosArrayAcc(),
   PPosArrayAcc(),
   PiPosArrayAcc(),
   He3NegArrayAcc(),
   PNegArrayAcc(),
   PiNegArrayAcc(),
   He3PosCounterAcc(),
   He3NegCounterAcc(),
   PPosCounterAcc(),
   PNegCounterAcc(),
   PiPosCounterAcc(),
   PiNegCounterAcc(),
   fBetheParamsHe(),
   fBetheParamsT(),  
   fMCtrue(0),
   vertex(),
   primVertex(),
   secvertexer(),
   secVertex(),
   tertvertexer(),
   tertVertex(),
   PrimVertex(),
   SecVertex(),
   TertVertex(),
   TertVertex2(),
   track1(),
   track2(),
   track3(),
   track4(),
   exTrack(),
   trkArray(),
   trkArray1(),
   ptot1(0),
   sign1(0),
   ptot2(0),
   sign2(0),
   ptot3(0),
   sign3(0),
   ptot4(0),
   sign4(0),
   He3Pos1(0),
   He3Neg1(0),
   pPos2(0),
   pNeg2(0),
   piPos3(0),
   piNeg3(0),
   piPos4(0),
   piNeg4(0),    
   cov(),
   cov0(),
   cov1(),
   cov2(),
   id(),
   pxpypz(),
   xyz(),
   sign(0),
   dd(),
   xthiss(),
   xpp(),
   h(),
   fthetaP(-99), 
   fthetaN(-99),
   //
   noCombinatoricsMC(0),
   lessCombinatoricsMC(0),
   stackN(0),
   PDGCodeMother(),
   ParticleMother(),
   FirstDaughter(),
   SecondDaughter(),
   ThirdDaughter(),
   FourthDaughter(),
   label1(),
   labelMother1(),
   labelGrandMother1(),
   ParticleMother1(),
   ParticleGrandMother1(),
   label2(),
   labelMother2(),
   labelGrandMother2(),
   ParticleMother2(),
   ParticleGrandMother2(),
   label3(),
   labelMother3(),
   labelGrandMother3(),
   ParticleMother3(),
   ParticleGrandMother3(),
   label4(),
   labelMother4(),
   labelGrandMother4(),
   ParticleMother4(),
   ParticleGrandMother4(),
   //
   lorentzsum(NULL),
   lorentzsum2(NULL),
   sublorentzsum(NULL),
   sublorentzsum2(NULL),
   particle1(NULL),
   particle2(NULL),
   particle3(NULL),
   particle4(NULL),
   //Histos
   fHistogramList(NULL),
   fHistdEdx(0),
   fHistNumEvents(0),
   fHistTrigger(0),
   hInvMass4LHe(0),
   hInvMass4LLH(0),
   //saved in tree
   fCentrality(-99),
   frunnumber(-99),
   fPeriod(0),
   fMagneticField(0),
   feventclass(0),
   mctruth(0),
   fTrigMB(-99),			
   fTrigHMV0(-99),
   fTrigHMSPD(-99),
   fTrigHNU(0),
   fTrigHQU(0),
   fDCA2B(-99),
   fDCA2Bo(-99),
   fDCA3B1(-99),
   fDCA3B2(-99),
   fDCA3B3(-99),
   fDecAngle(-99),
   farmalpha(-99),
   farmpt(-99),
   fPA(-99),
   fSubPA(-99),
   fSubPA2(-99),
   fSubPA3(-99),
   fSubPA4(-99),
   fPrimVertexX(-99),
   fPrimVertexY(-99),
   fPrimVertexZ(-99),
   fSecVertexX(-99),
   fSecVertexY(-99),
   fSecVertexZ(-99),
   fTertVertexX(-99),
   fTertVertexY(-99),
   fTertVertexZ(-99),
   fTertVertex2X(-99),
   fTertVertex2Y(-99),
   fTertVertex2Z(-99),
   //
   fPDGMother(-99),
   fChargeMother(-99),
   fmMother(-99),
   fmMother2(-99),
   fEMother(-99),
   fpxMother(-99),
   fpyMother(-99),
   fpzMother(-99),
   fptMother(-99),
   fpMother(-99),
   fyMother(-99),
   fctMother(-99),
   //
   fmSubMother(-99),
   fESubMother(-99),
   fpxSubMother(-99),
   fpySubMother(-99),
   fpzSubMother(-99),
   fptSubMother(-99),
   fpSubMother(-99),
   fySubMother(-99),
   fctSubMother(-99),
   //
   fEDaughter(-99),
   fpDaughter(-99),
   fptDaughter(-99),
   fpxDaughter(-99),
   fpyDaughter(-99),
   fpzDaughter(-99),
   fyDaughter(-99),
   fdEdxDaughter(-99),
   fdEdxSigmaDaughter(-99),
   fDcaDaughter(-99),
   fDcaDaughtero(-99),
   fDcazDaughter(-99),
   fDcaSecDaughter(-99),
   fImParDaughter(-99),
   fImPar2Daughter(-99),
   fImParzDaughter(-99),
   fImParz2Daughter(-99),
   fNclsDaughter(-99),
   fChi2Daughter(-99),
   fNclsITSDaughter(-99),
   fEtaDaughter(-99),
   fPhiDaughter(-99),
   fGeoLengthDaughter(-99),
   fTOFSignalDaughter(-99),
   fSigmaYXDaughter(-99),
   fSigmaXYZDaughter(-99),
   fSigmaZDaughter(-99),
   fPtUncertDaughter(-99),
   fTPCRefitDaughter(-99),
   fITSRefitDaughter(-99),
   fPropDCADaughter(-1),
   //
   fEDaughter1(-99),
   fpDaughter1(-99),
   fptDaughter1(-99),
   fpxDaughter1(-99),
   fpyDaughter1(-99),
   fpzDaughter1(-99),
   fyDaughter1(-99),
   fdEdxDaughter1(-99),
   fdEdxSigmaDaughter1(-99),
   fDcaDaughter1(-99),
   fDcaDaughter1o(-99),
   fDcazDaughter1(-99),
   fDcaSecDaughter1(-99),
   fImParDaughter1(-99),
   fImPar2Daughter1(-99),
   fImParzDaughter1(-99),
   fImParz2Daughter1(-99),
   fNclsDaughter1(-99),
   fChi2Daughter1(-99),
   fNclsITSDaughter1(-99),
   fEtaDaughter1(-99),
   fPhiDaughter1(-99),
   fGeoLengthDaughter1(-99),
   fTOFSignalDaughter1(-99),
   fSigmaYXDaughter1(-99),
   fSigmaXYZDaughter1(-99),
   fSigmaZDaughter1(-99),
   fPtUncertDaughter1(-99),
   fTPCRefitDaughter1(-99),
   fITSRefitDaughter1(-99),
   fPropDCADaughter1(-1),
   //
   fEDaughter2(-99),
   fpDaughter2(-99),
   fptDaughter2(-99),
   fpxDaughter2(-99),
   fpyDaughter2(-99),
   fpzDaughter2(-99),
   fyDaughter2(-99),
   fdEdxDaughter2(-99),
   fdEdxSigmaDaughter2(-99),
   fDcaDaughter2(-99),
   fDcaDaughter2o(-99),
   fDcazDaughter2(-99),
   fDcaSecDaughter2(-99),
   fImParDaughter2(-99),
   fImPar2Daughter2(-99),
   fImParzDaughter2(-99),
   fImParz2Daughter2(-99),
   fNclsDaughter2(-99),
   fChi2Daughter2(-99),
   fNclsITSDaughter2(-99),
   fEtaDaughter2(-99),
   fPhiDaughter2(-99),
   fGeoLengthDaughter2(-99),
   fTOFSignalDaughter2(-99),
   fSigmaYXDaughter2(-99),
   fSigmaXYZDaughter2(-99),
   fSigmaZDaughter2(-99),
   fPtUncertDaughter2(-99),
   fTPCRefitDaughter2(-99),
   fITSRefitDaughter2(-99),
   fPropDCADaughter2(-1),
   //
   fEDaughter3(-99),
   fpDaughter3(-99),
   fptDaughter3(-99),
   fpxDaughter3(-99),
   fpyDaughter3(-99),
   fpzDaughter3(-99),
   fyDaughter3(-99),
   fdEdxDaughter3(-99),
   fdEdxSigmaDaughter3(-99),
   fDcaDaughter3(-99),
   fDcaDaughter3o(-99),
   fDcazDaughter3(-99),
   fDcaSecDaughter3(-99),
   fImParDaughter3(-99),
   fImPar2Daughter3(-99),
   fImParzDaughter3(-99),
   fImParz2Daughter3(-99),
   fNclsDaughter3(-99),
   fChi2Daughter3(-99),
   fNclsITSDaughter3(-99),
   fEtaDaughter3(-99),
   fPhiDaughter3(-99),
   fGeoLengthDaughter3(-99),
   fTOFSignalDaughter3(-99),
   fSigmaYXDaughter3(-99),
   fSigmaXYZDaughter3(-99),
   fSigmaZDaughter3(-99),
   fPtUncertDaughter3(-99),
   fTPCRefitDaughter3(-99),
   fITSRefitDaughter3(-99),
   fPropDCADaughter3(-1),
   //
   fEDaughter4(-99),
   fpDaughter4(-99),
   fptDaughter4(-99),
   fpxDaughter4(-99),
   fpyDaughter4(-99),
   fpzDaughter4(-99),
   fyDaughter4(-99),
   fDcaDaughter4(-99),
   fDcazDaughter4(-99),
   fDcaSecDaughter4(-99),
   fImParDaughter4(-99),
   fImPar2Daughter4(-99),
   fImParzDaughter4(-99),
   fImParz2Daughter4(-99),
   fSigmaYXDaughter4(-99),
   fSigmaXYZDaughter4(-99),
   fSigmaZDaughter4(-99),
   fPtUncertDaughter4(-99),
   fPropDCADaughter4(-1),
   fdEdxSigmaPion(-99),
   fdEdxSigmaDeuteron(-99),
   fdEdxSigmaTriton(-99),
   fdEdxSigmaAlpha(-99)
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
  DefineOutput(4, TTree::Class());
  DefineOutput(5, TTree::Class());

  // _________________________________________________ //
  // __ track cut settings __ //
  // __ 3He tracks __ //
  trackCutsNuclei = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");
  trackCutsNuclei->SetEtaRange(-0.9,0.9);
  trackCutsNuclei->SetAcceptKinkDaughters(kFALSE);
  trackCutsNuclei->SetRequireTPCRefit(kTRUE);
  trackCutsNuclei->SetMaxChi2PerClusterTPC(5);
  trackCutsNuclei->SetMinNClustersTPC(80);
  //trackCutsNuclei->SetMaxRel1PtUncertainty(0.1);
  //trackCutsNuclei->SetPtRange(0.0, 10.0);

  // __ pion tracks __ //
  trackCutsStrong = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");
  trackCutsStrong->SetEtaRange(-0.9,0.9);
  trackCutsStrong->SetAcceptKinkDaughters(kFALSE);
  trackCutsStrong->SetRequireTPCRefit(kTRUE);
  trackCutsStrong->SetMaxChi2PerClusterTPC(5);
  trackCutsStrong->SetMinNClustersTPC(60);
  //trackCutsStrong->SetMaxRel1PtUncertainty(0.2);
  trackCutsStrong->SetPtRange(0.0, 1.0);

  // __ proton tracks __ //
  trackCutsSoft = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");
  trackCutsSoft->SetEtaRange(-0.9,0.9);
  trackCutsSoft->SetAcceptKinkDaughters(kFALSE);
  trackCutsSoft->SetRequireTPCRefit(kTRUE);
  trackCutsSoft->SetMaxChi2PerClusterTPC(5);
  trackCutsSoft->SetMinNClustersTPC(60);
  //trackCutsSoft->SetMaxRel1PtUncertainty(0.1);
  //trackCutsSoft->SetPtRange(0.0, 5.0);
}
// _________________________________________________ //
// __ Destructor __ //
AliAnalysisTaskDoubleHypNucTree::~AliAnalysisTaskDoubleHypNucTree() {

  ResetVals("Event");

}
// _________________________________________________ //
const Int_t AliAnalysisTaskDoubleHypNucTree::fgkPdgCode[] = {
  211,                //PionPlus
  -211,               //PionMinus
  2212,               //Proton
  -2212,              //Anti-Proton
  1000010020,         //Deuteron
  -1000010020,        //Anti-Deuteron
  1000010030,         //Triton
  -1000010030,        //Anti-Triton
  1000020030,         //Helium 3
  -1000020030,        //Anti-Helium 3
  1000020040,         //Helium 4
  -1000020040,        //Anti-Helium 4
  1010010030,	      //HyperHydrogen 3
  -1010010030,	      //Anti-HyperHydrogen3
  1010010040,	      //HyperHydrogen 4
  -1010010040,	      //Anti-HyperHydrogen 4
  1010020040, 	      //HyperHelium 4
  -1010020040, 	      //Anti-HyperHelium 4
  1010020050,	      //HyperHelium 5
  -1010020050,	      //Anti-HyperHelium 5
  1020010040,	      //DoubleHyperHydrogen 4
  -1020010040	      //Anti-DoubleHyperHydrogen 4 
};
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::UserCreateOutputObjects() {
  fInputHandler = dynamic_cast<AliESDInputHandler*>
    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!fInputHandler) {
    AliError("Could not get ESD InputHandler.\n");
    return;
  }
  fPID = fInputHandler->GetESDpid();
  if (!fPID) {
    AliError("Could not get PID response.\n");
    return;
  }
  // _________________________________________________ //
  // __ Histograms __ //
  fHistdEdx = new TH2F("fHistdEdX","dE/dx;#frac{#it{p}}{z} (GeV/#it{c});TPC Signal (a.u.)",1000,-5.0,5.0,1000,0.0,2000);
	
  fHistNumEvents = new TH1F("fHistNumEvents","Number of Events",2,0,2);
  fHistNumEvents->GetXaxis()->SetBinLabel(1,"before PhysSel");
  fHistNumEvents->GetXaxis()->SetBinLabel(2,"after PhysSel");

  fHistTrigger = new TH1F("fHistTrigger","Trigger",7,0,7);
  fHistTrigger->GetXaxis()->SetBinLabel(1,"other");
  fHistTrigger->GetXaxis()->SetBinLabel(2,"kINT7");
  fHistTrigger->GetXaxis()->SetBinLabel(3,"kHighMultV0");
  fHistTrigger->GetXaxis()->SetBinLabel(4,"kHighMultSPD");
  fHistTrigger->GetXaxis()->SetBinLabel(5,"HNU");
  fHistTrigger->GetXaxis()->SetBinLabel(6,"HQU");
  fHistTrigger->GetXaxis()->SetBinLabel(7,"HJT");

  hInvMass4LHe = new TH1F("hInvMass4LHe", "", 120, 3.87, 3.99);
  hInvMass4LHe->GetXaxis()->SetTitle("Inv. Mass (GeV/#it{c}^{2})");
  hInvMass4LHe->GetYaxis()->SetTitle("Count / (2 MeV/#it{c}^{2})");

  hInvMass4LLH = new TH1F("hInvMass4LLH", "", 120, 4.04, 4.16);
  hInvMass4LLH->GetXaxis()->SetTitle("Inv. Mass (GeV/#it{c}^{2})");
  hInvMass4LLH->GetYaxis()->SetTitle("Count / (2 MeV/#it{c}^{2})");

  fHistogramList = new TList();
  fHistogramList->SetOwner(kTRUE);
  fHistogramList->SetName(GetName());
  fHistogramList->Add(fHistdEdx);
  fHistogramList->Add(fHistNumEvents);
  fHistogramList->Add(fHistTrigger);
  fHistogramList->Add(hInvMass4LHe);
  fHistogramList->Add(hInvMass4LLH);
  fEventCuts.AddQAplotsToList(fHistogramList);
  // _________________________________________________ //
  // __ associated Tree __ //
  fTree = new TTree("fTree", "fTree");    
  fTree->Branch("fPeriod",             &fPeriod,             "fPeriod/I");
  fTree->Branch("frunnumber",          &frunnumber,          "frunnumber/I");
  fTree->Branch("feventclass",         &feventclass,         "feventclass/I");
  fTree->Branch("fCentrality",         &fCentrality,         "fCentrality/I");
  fTree->Branch("fMagneticField",      &fMagneticField,      "fMagneticField/I");  	
  fTree->Branch("fPrimVertexX",        &fPrimVertexX,        "fPrimVertexX/F");
  fTree->Branch("fPrimVertexY",        &fPrimVertexY,        "fPrimVertexY/F");
  fTree->Branch("fPrimVertexZ",        &fPrimVertexZ,        "fPrimVertexZ/F");
  fTree->Branch("fSecVertexX",         &fSecVertexX,         "fSecVertexX/F");
  fTree->Branch("fSecVertexY",         &fSecVertexY,         "fSecVertexY/F");
  fTree->Branch("fSecVertexZ",         &fSecVertexZ,         "fSecVertexZ/F");
  fTree->Branch("fTertVertexX",        &fTertVertexX,        "fTertVertexX/F");
  fTree->Branch("fTertVertexY",        &fTertVertexY,        "fTertVertexY/F");
  fTree->Branch("fTertVertexZ",        &fTertVertexZ,        "fTertVertexZ/F");
  fTree->Branch("fTertVertex2X",       &fTertVertex2X,       "fTertVertex2X/F");
  fTree->Branch("fTertVertex2Y",       &fTertVertex2Y,       "fTertVertex2Y/F");
  fTree->Branch("fTertVertex2Z",       &fTertVertex2Z,       "fTertVertex2Z/F");
  fTree->Branch("fTrigMB",             &fTrigMB,             "fTrigMB/I");
  fTree->Branch("fTrigHMV0",           &fTrigHMV0,           "fTrigHMV0/I");
  fTree->Branch("fTrigHMSPD",          &fTrigHMSPD,          "fTrigHMSPD/I");
  fTree->Branch("fTrigHNU",            &fTrigHNU,            "fTrigHNU/I");
  fTree->Branch("fTrigHQU",            &fTrigHQU,            "fTrigHQU/I");
  fTree->Branch("fDCA2B",              &fDCA2B,              "fDCA2B/F");
  fTree->Branch("fDCA2Bo",             &fDCA2Bo,             "fDCA2Bo/F");
  fTree->Branch("fDCA3B1",             &fDCA3B1,             "fDCA3B1/F");
  fTree->Branch("fDCA3B2",             &fDCA3B2,             "fDCA3B2/F");
  fTree->Branch("fDCA3B3",             &fDCA3B3,             "fDCA3B3/F");
  fTree->Branch("fPA",                 &fPA,                 "fPA/F");
  fTree->Branch("fSubPA",              &fSubPA,              "fSubPA/F");
  fTree->Branch("fSubPA2",             &fSubPA2,             "fSubPA2/F");
  fTree->Branch("fSubPA3",             &fSubPA3,             "fSubPA3/F");
  fTree->Branch("fSubPA4",             &fSubPA4,             "fSubPA4/F");
  fTree->Branch("fDecAngle",           &fDecAngle,           "fDecAngle/F");
  fTree->Branch("farmalpha",           &farmalpha,           "farmalpha/F");
  fTree->Branch("farmpt",              &farmpt,              "farmpt/F");
  fTree->Branch("fPDGMother",          &fPDGMother,          "fPDGMother/I");
  fTree->Branch("fChargeMother",       &fChargeMother,       "fChargeMother/I");
  fTree->Branch("mctruth",             &mctruth,             "mctruth/I"); 
  fTree->Branch("fmMother",            &fmMother,            "fmMother/F");
  fTree->Branch("fmMother2",           &fmMother2,           "fmMother2/F");
  fTree->Branch("fEMother",            &fEMother,            "fEMother/F");
  fTree->Branch("fpxMother",           &fpxMother,           "fpxMother/F");
  fTree->Branch("fpyMother",           &fpyMother,           "fpyMother/F");
  fTree->Branch("fpzMother",           &fpzMother,           "fpzMother/F");
  fTree->Branch("fptMother",           &fptMother,           "fptMother/F");
  fTree->Branch("fpMother",            &fpMother,            "fpMother/F");
  fTree->Branch("fyMother",            &fyMother,            "fyMother/F");
  fTree->Branch("fctMother",           &fctMother,           "fctMother/F");
  fTree->Branch("fmSubMother",         &fmSubMother,         "fmSubMother/F");
  fTree->Branch("fESubMother",         &fESubMother,         "fESubMother/F");
  fTree->Branch("fpxSubMother",        &fpxSubMother,        "fpxSubMother/F");
  fTree->Branch("fpySubMother",        &fpySubMother,        "fpySubMother/F");
  fTree->Branch("fpzSubMother",        &fpzSubMother,        "fpzSubMother/F");
  fTree->Branch("fptSubMother",        &fptSubMother,        "fptSubMother/F");
  fTree->Branch("fpSubMother",         &fpSubMother,         "fpSubMother/F");
  fTree->Branch("fySubMother",         &fySubMother,         "fySubMother/F");
  fTree->Branch("fctSubMother",        &fctSubMother,        "fctSubMother/F");
  fTree->Branch("fEDaughter",          &fEDaughter,          "fEDaughter/F");
  fTree->Branch("fpDaughter",          &fpDaughter,          "fpDaughter/F");
  fTree->Branch("fptDaughter",         &fptDaughter,         "fptDaughter/F");
  fTree->Branch("fpxDaughter",         &fpxDaughter,         "fpxDaughter/F");
  fTree->Branch("fpyDaughter",         &fpyDaughter,         "fpyDaughter/F");
  fTree->Branch("fpzDaughter",         &fpzDaughter,         "fpzDaughter/F");
  fTree->Branch("fyDaughter",          &fyDaughter,          "fyDaughter/F");
  fTree->Branch("fdEdxDaughter",       &fdEdxDaughter,       "fdEdxDaughter/F");
  fTree->Branch("fdEdxSigmaDaughter",  &fdEdxSigmaDaughter,  "fdEdxSigmaDaughter/F");
  fTree->Branch("fDcaDaughter",        &fDcaDaughter,        "fDcaDaughter/F");
  fTree->Branch("fDcaDaughtero",       &fDcaDaughtero,       "fDcaDaughtero/F");
  fTree->Branch("fDcazDaughter",       &fDcazDaughter,       "fDcazDaughter/F");
  fTree->Branch("fDcaSecDaughter",     &fDcaSecDaughter,     "fDcaSecDaughter/F");
  fTree->Branch("fImParDaughter",      &fImParDaughter,      "fImParDaughter/F");
  fTree->Branch("fImPar2Daughter",     &fImPar2Daughter,     "fImPar2Daughter/F");
  fTree->Branch("fImParzDaughter",     &fImParzDaughter,     "fImParzDaughter/F");
  fTree->Branch("fImParz2Daughter",    &fImParz2Daughter,    "fImParz2Daughter/F");
  fTree->Branch("fNclsDaughter",       &fNclsDaughter,       "fNclsDaughter/I");
  fTree->Branch("fChi2Daughter",       &fChi2Daughter,       "fChi2Daughter/F");
  fTree->Branch("fNclsITSDaughter",    &fNclsITSDaughter,    "fNclsITSDaughter/I");
  fTree->Branch("fEtaDaughter",        &fEtaDaughter,        "fEtaDaughter/F");
  fTree->Branch("fPhiDaughter",        &fPhiDaughter,        "fPhiDaughter/F");
  fTree->Branch("fGeoLengthDaughter",  &fGeoLengthDaughter,  "fGeoLengthDaughter/F");
  fTree->Branch("fTOFSignalDaughter",  &fTOFSignalDaughter,  "fTOFSignalDaughter/F");
  fTree->Branch("fSigmaYXDaughter",    &fSigmaYXDaughter,    "fSigmaYXDaughter/F");
  fTree->Branch("fSigmaXYZDaughter",   &fSigmaXYZDaughter,   "fSigmaXYZDaughter/F");
  fTree->Branch("fSigmaZDaughter",     &fSigmaZDaughter,     "fSigmaZDaughter/F");
  fTree->Branch("fPtUncertDaughter",   &fPtUncertDaughter,   "fPtUncertDaughter/F");
  fTree->Branch("fTPCRefitDaughter",   &fTPCRefitDaughter,   "fTPCRefitDaughter/I");
  fTree->Branch("fITSRefitDaughter",   &fITSRefitDaughter,   "fITSRefitDaughter/I");
  fTree->Branch("fPropDCADaughter",    &fPropDCADaughter,    "fPropDCADaughter/I");
  fTree->Branch("fEDaughter1",         &fEDaughter1,         "fEDaughter1/F");
  fTree->Branch("fpDaughter1",         &fpDaughter1,         "fpDaughter1/F");
  fTree->Branch("fptDaughter1",        &fptDaughter1,        "fptDaughter1/F");
  fTree->Branch("fpxDaughter1",        &fpxDaughter1,        "fpxDaughter1/F");
  fTree->Branch("fpyDaughter1",        &fpyDaughter1,        "fpyDaughter1/F");
  fTree->Branch("fpzDaughter1",        &fpzDaughter1,        "fpzDaughter1/F");
  fTree->Branch("fyDaughter1",         &fyDaughter1,         "fyDaughter1/F");
  fTree->Branch("fdEdxDaughter1",      &fdEdxDaughter1,      "fdEdxDaughter1/F");
  fTree->Branch("fdEdxSigmaDaughter1", &fdEdxSigmaDaughter1, "fdEdxSigmaDaughter1/F");
  fTree->Branch("fDcaDaughter1",       &fDcaDaughter1,       "fDcaDaughter1/F");
  fTree->Branch("fDcaDaughter1o",      &fDcaDaughter1o,      "fDcaDaughter1o/F");
  fTree->Branch("fDcazDaughter1",      &fDcazDaughter1,      "fDcazDaughter1/F");
  fTree->Branch("fDcaSecDaughter1",    &fDcaSecDaughter1,    "fDcaSecDaughter1/F");
  fTree->Branch("fImParDaughter1",     &fImParDaughter1,     "fImParDaughter1/F");
  fTree->Branch("fImPar2Daughter1",    &fImPar2Daughter1,    "fImPar2Daughter1/F");
  fTree->Branch("fImParzDaughter1",    &fImParzDaughter1,    "fImParzDaughter1/F");
  fTree->Branch("fImParz2Daughter1",   &fImParz2Daughter1,   "fImParz2Daughter1/F");
  fTree->Branch("fNclsDaughter1",      &fNclsDaughter1,      "fNclsDaughter1/I");
  fTree->Branch("fChi2Daughter1",      &fChi2Daughter1,      "fChi2Daughter1/F");
  fTree->Branch("fNclsITSDaughter1",   &fNclsITSDaughter1,   "fNclsITSDaughter1/I");
  fTree->Branch("fEtaDaughter1",       &fEtaDaughter1,       "fEtaDaughter1/F");
  fTree->Branch("fPhiDaughter1",       &fPhiDaughter1,       "fPhiDaughter1/F");
  fTree->Branch("fGeoLengthDaughter1", &fGeoLengthDaughter1, "fGeoLengthDaughter1/F");
  fTree->Branch("fTOFSignalDaughter1", &fTOFSignalDaughter1, "fTOFSignalDaughter1/F");
  fTree->Branch("fSigmaYXDaughter1",   &fSigmaYXDaughter1,   "fSigmaYXDaughter1/F");
  fTree->Branch("fSigmaXYZDaughter1",  &fSigmaXYZDaughter1,  "fSigmaXYZDaughter1/F");
  fTree->Branch("fSigmaZDaughter1",    &fSigmaZDaughter1,    "fSigmaZDaughter1/F");
  fTree->Branch("fPtUncertDaughter1",  &fPtUncertDaughter1,  "fPtUncertDaughter1/F");
  fTree->Branch("fTPCRefitDaughter1",  &fTPCRefitDaughter1,  "fTPCRefitDaughter1/I");
  fTree->Branch("fITSRefitDaughter1",  &fITSRefitDaughter1,  "fITSRefitDaughter1/I");
  fTree->Branch("fPropDCADaughter1",   &fPropDCADaughter1,   "fPropDCADaughter1/I");
  fTree->Branch("fEDaughter2",         &fEDaughter2,         "fEDaughter2/F");
  fTree->Branch("fpDaughter2",         &fpDaughter2,         "fpDaughter2/F");
  fTree->Branch("fptDaughter2",        &fptDaughter2,        "fptDaughter2/F");
  fTree->Branch("fpxDaughter2",        &fpxDaughter2,        "fpxDaughter2/F");
  fTree->Branch("fpyDaughter2",        &fpyDaughter2,        "fpyDaughter2/F");
  fTree->Branch("fpzDaughter2",        &fpzDaughter2,        "fpzDaughter2/F");
  fTree->Branch("fyDaughter2",         &fyDaughter2,         "fyDaughter2/F");
  fTree->Branch("fdEdxDaughter2",      &fdEdxDaughter2,      "fdEdxDaughter2/F");
  fTree->Branch("fdEdxSigmaDaughter2", &fdEdxSigmaDaughter2, "fdEdxSigmaDaughter2/F");
  fTree->Branch("fDcaDaughter2",       &fDcaDaughter2,       "fDcaDaughter2/F");
  fTree->Branch("fDcaDaughter2o",      &fDcaDaughter2o,      "fDcaDaughter2o/F");
  fTree->Branch("fDcazDaughter2",      &fDcazDaughter2,      "fDcazDaughter2/F");
  fTree->Branch("fDcaSecDaughter2",    &fDcaSecDaughter2,    "fDcaSecDaughter2/F");
  fTree->Branch("fImParDaughter2",     &fImParDaughter2,     "fImParDaughter2/F");
  fTree->Branch("fImPar2Daughter2",    &fImPar2Daughter2,    "fImPar2Daughter2/F");
  fTree->Branch("fImParzDaughter2",    &fImParzDaughter2,    "fImParzDaughter2/F");
  fTree->Branch("fImParz2Daughter2",   &fImParz2Daughter2,   "fImParz2Daughter2/F");
  fTree->Branch("fNclsDaughter2",      &fNclsDaughter2,      "fNclsDaughter2/I");
  fTree->Branch("fChi2Daughter2",      &fChi2Daughter2,      "fChi2Daughter2/F");
  fTree->Branch("fNclsITSDaughter2",   &fNclsITSDaughter2,   "fNclsITSDaughter2/I");
  fTree->Branch("fEtaDaughter2",       &fEtaDaughter2,       "fEtaDaughter2/F");
  fTree->Branch("fPhiDaughter2",       &fPhiDaughter2,       "fPhiDaughter2/F");
  fTree->Branch("fGeoLengthDaughter2", &fGeoLengthDaughter2, "fGeoLengthDaughter2/F");
  fTree->Branch("fTOFSignalDaughter2", &fTOFSignalDaughter2, "fTOFSignalDaughter2/F");
  fTree->Branch("fSigmaYXDaughter2",   &fSigmaYXDaughter2,   "fSigmaYXDaughter2/F");
  fTree->Branch("fSigmaXYZDaughter2",  &fSigmaXYZDaughter2,  "fSigmaXYZDaughter2/F");
  fTree->Branch("fSigmaZDaughter2",    &fSigmaZDaughter2,    "fSigmaZDaughter2/F");
  fTree->Branch("fPtUncertDaughter2",  &fPtUncertDaughter2,  "fPtUncertDaughter2/F");
  fTree->Branch("fTPCRefitDaughter2",  &fTPCRefitDaughter2,  "fTPCRefitDaughter2/I");
  fTree->Branch("fITSRefitDaughter2",  &fITSRefitDaughter2,  "fITSRefitDaughter2/I");
  fTree->Branch("fPropDCADaughter2",   &fPropDCADaughter2,   "fPropDCADaughter2/I");
  fTree->Branch("fEDaughter3",         &fEDaughter3,         "fEDaughter3/F");
  fTree->Branch("fpDaughter3",         &fpDaughter3,         "fpDaughter3/F");
  fTree->Branch("fptDaughter3",        &fptDaughter3,        "fptDaughter3/F");
  fTree->Branch("fpxDaughter3",        &fpxDaughter3,        "fpxDaughter3/F");
  fTree->Branch("fpyDaughter3",        &fpyDaughter3,        "fpyDaughter3/F");
  fTree->Branch("fpzDaughter3",        &fpzDaughter3,        "fpzDaughter3/F");
  fTree->Branch("fyDaughter3",         &fyDaughter3,         "fyDaughter3/F");
  fTree->Branch("fdEdxDaughter3",      &fdEdxDaughter3,      "fdEdxDaughter3/F");
  fTree->Branch("fdEdxSigmaDaughter3", &fdEdxSigmaDaughter3, "fdEdxSigmaDaughter3/F");
  fTree->Branch("fDcaDaughter3",       &fDcaDaughter3,       "fDcaDaughter3/F");
  fTree->Branch("fDcaDaughter3o",      &fDcaDaughter3o,      "fDcaDaughter3o/F");
  fTree->Branch("fDcazDaughter3",      &fDcazDaughter3,      "fDcazDaughter3/F");
  fTree->Branch("fDcaSecDaughter3",    &fDcaSecDaughter3,    "fDcaSecDaughter3/F");
  fTree->Branch("fImParDaughter3",     &fImParDaughter3,     "fImParDaughter3/F");
  fTree->Branch("fImPar2Daughter3",    &fImPar2Daughter3,    "fImPar2Daughter3/F");
  fTree->Branch("fImParzDaughter3",    &fImParzDaughter3,    "fImParzDaughter3/F");
  fTree->Branch("fImParz2Daughter3",   &fImParz2Daughter3,   "fImParz2Daughter3/F");
  fTree->Branch("fNclsDaughter3",      &fNclsDaughter3,      "fNclsDaughter3/I");
  fTree->Branch("fChi2Daughter3",      &fChi2Daughter3,      "fChi2Daughter3/F");
  fTree->Branch("fNclsITSDaughter3",   &fNclsITSDaughter3,   "fNclsITSDaughter3/I");
  fTree->Branch("fEtaDaughter3",       &fEtaDaughter3,       "fEtaDaughter3/F");
  fTree->Branch("fPhiDaughter3",       &fPhiDaughter3,       "fPhiDaughter3/F");
  fTree->Branch("fGeoLengthDaughter3", &fGeoLengthDaughter3, "fGeoLengthDaughter3/F");
  fTree->Branch("fTOFSignalDaughter3", &fTOFSignalDaughter3, "fTOFSignalDaughter3/F");
  fTree->Branch("fSigmaYXDaughter3",   &fSigmaYXDaughter3,   "fSigmaYXDaughter3/F");
  fTree->Branch("fSigmaXYZDaughter3",  &fSigmaXYZDaughter3,  "fSigmaXYZDaughter3/F");
  fTree->Branch("fSigmaZDaughter3",    &fSigmaZDaughter3,    "fSigmaZDaughter3/F");
  fTree->Branch("fPtUncertDaughter3",  &fPtUncertDaughter3,  "fPtUncertDaughter3/F");
  fTree->Branch("fTPCRefitDaughter3",  &fTPCRefitDaughter3,  "fTPCRefitDaughter3/I");
  fTree->Branch("fITSRefitDaughter3",  &fITSRefitDaughter3,  "fITSRefitDaughter3/I");
  fTree->Branch("fPropDCADaughter3",   &fPropDCADaughter3,   "fPropDCADaughter3/I");
  fTree->Branch("fEDaughter4",         &fEDaughter4,         "fEDaughter4/F");
  fTree->Branch("fpDaughter4",         &fpDaughter4,         "fpDaughter4/F");
  fTree->Branch("fptDaughter4",        &fptDaughter4,        "fptDaughter4/F");
  fTree->Branch("fpxDaughter4",        &fpxDaughter4,        "fpxDaughter4/F");
  fTree->Branch("fpyDaughter4",        &fpyDaughter4,        "fpyDaughter4/F");
  fTree->Branch("fpzDaughter4",        &fpzDaughter4,        "fpzDaughter4/F");
  fTree->Branch("fyDaughter4",         &fyDaughter4,         "fyDaughter4/F");
  fTree->Branch("fDcaDaughter4",       &fDcaDaughter4,       "fDcaDaughter4/F");
  fTree->Branch("fDcazDaughter4",      &fDcazDaughter4,      "fDcazDaughter4/F");
  fTree->Branch("fDcaSecDaughter4",    &fDcaSecDaughter4,    "fDcaSecDaughter4/F");
  fTree->Branch("fImParDaughter4",     &fImParDaughter4,     "fImParDaughter4/F");
  fTree->Branch("fImPar2Daughter4",    &fImPar2Daughter4,    "fImPar2Daughter4/F");
  fTree->Branch("fImParzDaughter4",    &fImParzDaughter4,    "fImParzDaughter4/F");
  fTree->Branch("fImParz2Daughter4",   &fImParz2Daughter4,   "fImParz2Daughter4/F");
  fTree->Branch("fSigmaYXDaughter4",   &fSigmaYXDaughter4,   "fSigmaYXDaughter4/F");
  fTree->Branch("fSigmaXYZDaughter4",  &fSigmaXYZDaughter4,  "fSigmaXYZDaughter4/F");
  fTree->Branch("fSigmaZDaughter4",    &fSigmaZDaughter4,    "fSigmaZDaughter4/F");
  fTree->Branch("fPtUncertDaughter4",  &fPtUncertDaughter4,  "fPtUncertDaughter4/F");
  fTree->Branch("fPropDCADaughter4",   &fPropDCADaughter4,   "fPropDCADaughter4/I");
  fTree->Branch("fdEdxSigmaPion",      &fdEdxSigmaPion,      "fdEdxSigmaPion/F");
  fTree->Branch("fdEdxSigmaDeuteron",  &fdEdxSigmaDeuteron,  "fdEdxSigmaDeuteron/F");
  fTree->Branch("fdEdxSigmaTriton",    &fdEdxSigmaTriton,    "fdEdxSigmaTriton/F");
  fTree->Branch("fdEdxSigmaAlpha",     &fdEdxSigmaAlpha,     "fdEdxSigmaAlpha/F");
  // _________________________________________________ //
  // __ associated 4LHe tree __ //
  gTree = new TTree("gTree", "gTree");    
  gTree->Branch("fPeriod",             &fPeriod,             "fPeriod/I");
  gTree->Branch("frunnumber",          &frunnumber,          "frunnumber/I");
  gTree->Branch("feventclass",         &feventclass,         "feventclass/I");
  gTree->Branch("fCentrality",         &fCentrality,         "fCentrality/I");
  gTree->Branch("fMagneticField",      &fMagneticField,      "fMagneticField/I");  	
  gTree->Branch("fPrimVertexX",        &fPrimVertexX,        "fPrimVertexX/F");
  gTree->Branch("fPrimVertexY",        &fPrimVertexY,        "fPrimVertexY/F");
  gTree->Branch("fPrimVertexZ",        &fPrimVertexZ,        "fPrimVertexZ/F");
  gTree->Branch("fTertVertexX",        &fTertVertexX,        "fTertVertexX/F");
  gTree->Branch("fTertVertexY",        &fTertVertexY,        "fTertVertexY/F");
  gTree->Branch("fTertVertexZ",        &fTertVertexZ,        "fTertVertexZ/F");
  gTree->Branch("fTrigMB",             &fTrigMB,             "fTrigMB/I");
  gTree->Branch("fTrigHMV0",           &fTrigHMV0,           "fTrigHMV0/I");
  gTree->Branch("fTrigHMSPD",          &fTrigHMSPD,          "fTrigHMSPD/I");
  gTree->Branch("fTrigHNU",            &fTrigHNU,            "fTrigHNU/I");
  gTree->Branch("fTrigHQU",            &fTrigHQU,            "fTrigHQU/I");
  gTree->Branch("fDCA3B1",             &fDCA3B1,             "fDCA3B1/F");
  gTree->Branch("fDCA3B2",             &fDCA3B2,             "fDCA3B2/F");
  gTree->Branch("fDCA3B3",             &fDCA3B3,             "fDCA3B3/F");
  gTree->Branch("fSubPA",              &fSubPA,              "fSubPA/F");
  gTree->Branch("fPDGMother",          &fPDGMother,          "fPDGMother/I");
  gTree->Branch("fChargeMother",       &fChargeMother,       "fChargeMother/I");
  gTree->Branch("mctruth",             &mctruth,             "mctruth/I"); 
  gTree->Branch("fmSubMother",         &fmSubMother,         "fmSubMother/F");
  gTree->Branch("fESubMother",         &fESubMother,         "fESubMother/F");
  gTree->Branch("fpxSubMother",        &fpxSubMother,        "fpxSubMother/F");
  gTree->Branch("fpySubMother",        &fpySubMother,        "fpySubMother/F");
  gTree->Branch("fpzSubMother",        &fpzSubMother,        "fpzSubMother/F");
  gTree->Branch("fptSubMother",        &fptSubMother,        "fptSubMother/F");
  gTree->Branch("fpSubMother",         &fpSubMother,         "fpSubMother/F");
  gTree->Branch("fySubMother",         &fySubMother,         "fySubMother/F");
  gTree->Branch("fctSubMother",        &fctSubMother,        "fctSubMother/F");
  gTree->Branch("fEDaughter",          &fEDaughter,          "fEDaughter/F");
  gTree->Branch("fpDaughter",          &fpDaughter,          "fpDaughter/F");
  gTree->Branch("fptDaughter",         &fptDaughter,         "fptDaughter/F");
  gTree->Branch("fpxDaughter",         &fpxDaughter,         "fpxDaughter/F");
  gTree->Branch("fpyDaughter",         &fpyDaughter,         "fpyDaughter/F");
  gTree->Branch("fpzDaughter",         &fpzDaughter,         "fpzDaughter/F");
  gTree->Branch("fyDaughter",          &fyDaughter,          "fyDaughter/F");
  gTree->Branch("fdEdxDaughter",       &fdEdxDaughter,       "fdEdxDaughter/F");
  gTree->Branch("fdEdxSigmaDaughter",  &fdEdxSigmaDaughter,  "fdEdxSigmaDaughter/F");
  gTree->Branch("fDcaDaughter",        &fDcaDaughter,        "fDcaDaughter/F");
  gTree->Branch("fDcaDaughtero",       &fDcaDaughtero,       "fDcaDaughtero/F");
  gTree->Branch("fDcazDaughter",       &fDcazDaughter,       "fDcazDaughter/F");
  gTree->Branch("fDcaSecDaughter",     &fDcaSecDaughter,     "fDcaSecDaughter/F");
  gTree->Branch("fImParDaughter",      &fImParDaughter,      "fImParDaughter/F");
  gTree->Branch("fImParzDaughter",     &fImParzDaughter,     "fImParzDaughter/F");
  gTree->Branch("fNclsDaughter",       &fNclsDaughter,       "fNclsDaughter/I");
  gTree->Branch("fChi2Daughter",       &fChi2Daughter,       "fChi2Daughter/F");
  gTree->Branch("fNclsITSDaughter",    &fNclsITSDaughter,    "fNclsITSDaughter/I");
  gTree->Branch("fEtaDaughter",        &fEtaDaughter,        "fEtaDaughter/F");
  gTree->Branch("fPhiDaughter",        &fPhiDaughter,        "fPhiDaughter/F");
  gTree->Branch("fGeoLengthDaughter",  &fGeoLengthDaughter,  "fGeoLengthDaughter/F");
  gTree->Branch("fTOFSignalDaughter",  &fTOFSignalDaughter,  "fTOFSignalDaughter/F");
  gTree->Branch("fSigmaYXDaughter",    &fSigmaYXDaughter,    "fSigmaYXDaughter/F");
  gTree->Branch("fSigmaXYZDaughter",   &fSigmaXYZDaughter,   "fSigmaXYZDaughter/F");
  gTree->Branch("fSigmaZDaughter",     &fSigmaZDaughter,     "fSigmaZDaughter/F");
  gTree->Branch("fPtUncertDaughter",   &fPtUncertDaughter,   "fPtUncertDaughter/F");
  gTree->Branch("fTPCRefitDaughter",   &fTPCRefitDaughter,   "fTPCRefitDaughter/I");
  gTree->Branch("fITSRefitDaughter",   &fITSRefitDaughter,   "fITSRefitDaughter/I");
  gTree->Branch("fPropDCADaughter",    &fPropDCADaughter,    "fPropDCADaughter/I");
  gTree->Branch("fEDaughter1",         &fEDaughter1,         "fEDaughter1/F");
  gTree->Branch("fpDaughter1",         &fpDaughter1,         "fpDaughter1/F");
  gTree->Branch("fptDaughter1",        &fptDaughter1,        "fptDaughter1/F");
  gTree->Branch("fpxDaughter1",        &fpxDaughter1,        "fpxDaughter1/F");
  gTree->Branch("fpyDaughter1",        &fpyDaughter1,        "fpyDaughter1/F");
  gTree->Branch("fpzDaughter1",        &fpzDaughter1,        "fpzDaughter1/F");
  gTree->Branch("fyDaughter1",         &fyDaughter1,         "fyDaughter1/F");
  gTree->Branch("fdEdxDaughter1",      &fdEdxDaughter1,      "fdEdxDaughter1/F");
  gTree->Branch("fdEdxSigmaDaughter1", &fdEdxSigmaDaughter1, "fdEdxSigmaDaughter1/F");
  gTree->Branch("fDcaDaughter1",       &fDcaDaughter1,       "fDcaDaughter1/F");
  gTree->Branch("fDcaDaughter1o",      &fDcaDaughter1o,      "fDcaDaughter1o/F");
  gTree->Branch("fDcazDaughter1",      &fDcazDaughter1,      "fDcazDaughter1/F");
  gTree->Branch("fDcaSecDaughter1",    &fDcaSecDaughter1,    "fDcaSecDaughter1/F");
  gTree->Branch("fImParDaughter1",     &fImParDaughter1,     "fImParDaughter1/F");
  gTree->Branch("fImParzDaughter1",    &fImParzDaughter1,    "fImParzDaughter1/F");
  gTree->Branch("fNclsDaughter1",      &fNclsDaughter1,      "fNclsDaughter1/I");
  gTree->Branch("fChi2Daughter1",      &fChi2Daughter1,      "fChi2Daughter1/F");
  gTree->Branch("fNclsITSDaughter1",   &fNclsITSDaughter1,   "fNclsITSDaughter1/I");
  gTree->Branch("fEtaDaughter1",       &fEtaDaughter1,       "fEtaDaughter1/F");
  gTree->Branch("fPhiDaughter1",       &fPhiDaughter1,       "fPhiDaughter1/F");
  gTree->Branch("fGeoLengthDaughter1", &fGeoLengthDaughter1, "fGeoLengthDaughter1/F");
  gTree->Branch("fTOFSignalDaughter1", &fTOFSignalDaughter1, "fTOFSignalDaughter1/F");
  gTree->Branch("fSigmaYXDaughter1",   &fSigmaYXDaughter1,   "fSigmaYXDaughter1/F");
  gTree->Branch("fSigmaXYZDaughter1",  &fSigmaXYZDaughter1,  "fSigmaXYZDaughter1/F");
  gTree->Branch("fSigmaZDaughter1",    &fSigmaZDaughter1,    "fSigmaZDaughter1/F");
  gTree->Branch("fPtUncertDaughter1",  &fPtUncertDaughter1,  "fPtUncertDaughter1/F");
  gTree->Branch("fTPCRefitDaughter1",  &fTPCRefitDaughter1,  "fTPCRefitDaughter1/I");
  gTree->Branch("fITSRefitDaughter1",  &fITSRefitDaughter1,  "fITSRefitDaughter1/I");
  gTree->Branch("fPropDCADaughter1",   &fPropDCADaughter1,   "fPropDCADaughter1/I");
  gTree->Branch("fEDaughter2",         &fEDaughter2,         "fEDaughter2/F");
  gTree->Branch("fpDaughter2",         &fpDaughter2,         "fpDaughter2/F");
  gTree->Branch("fptDaughter2",        &fptDaughter2,        "fptDaughter2/F");
  gTree->Branch("fpxDaughter2",        &fpxDaughter2,        "fpxDaughter2/F");
  gTree->Branch("fpyDaughter2",        &fpyDaughter2,        "fpyDaughter2/F");
  gTree->Branch("fpzDaughter2",        &fpzDaughter2,        "fpzDaughter2/F");
  gTree->Branch("fyDaughter2",         &fyDaughter2,         "fyDaughter2/F");
  gTree->Branch("fdEdxDaughter2",      &fdEdxDaughter2,      "fdEdxDaughter2/F");
  gTree->Branch("fdEdxSigmaDaughter2", &fdEdxSigmaDaughter2, "fdEdxSigmaDaughter2/F");
  gTree->Branch("fDcaDaughter2",       &fDcaDaughter2,       "fDcaDaughter2/F");
  gTree->Branch("fDcaDaughter2o",      &fDcaDaughter2o,      "fDcaDaughter2o/F");
  gTree->Branch("fDcazDaughter2",      &fDcazDaughter2,      "fDcazDaughter2/F");
  gTree->Branch("fDcaSecDaughter2",    &fDcaSecDaughter2,    "fDcaSecDaughter2/F");
  gTree->Branch("fImParDaughter2",     &fImParDaughter2,     "fImParDaughter2/F");
  gTree->Branch("fImParzDaughter2",    &fImParzDaughter2,    "fImParzDaughter2/F");
  gTree->Branch("fNclsDaughter2",      &fNclsDaughter2,      "fNclsDaughter2/I");
  gTree->Branch("fChi2Daughter2",      &fChi2Daughter2,      "fChi2Daughter2/F");
  gTree->Branch("fNclsITSDaughter2",   &fNclsITSDaughter2,   "fNclsITSDaughter2/I");
  gTree->Branch("fEtaDaughter2",       &fEtaDaughter2,       "fEtaDaughter2/F");
  gTree->Branch("fPhiDaughter2",       &fPhiDaughter2,       "fPhiDaughter2/F");
  gTree->Branch("fGeoLengthDaughter2", &fGeoLengthDaughter2, "fGeoLengthDaughter2/F");
  gTree->Branch("fTOFSignalDaughter2", &fTOFSignalDaughter2, "fTOFSignalDaughter2/F");
  gTree->Branch("fSigmaYXDaughter2",   &fSigmaYXDaughter2,   "fSigmaYXDaughter2/F");
  gTree->Branch("fSigmaXYZDaughter2",  &fSigmaXYZDaughter2,  "fSigmaXYZDaughter2/F");
  gTree->Branch("fSigmaZDaughter2",    &fSigmaZDaughter2,    "fSigmaZDaughter2/F");
  gTree->Branch("fPtUncertDaughter2",  &fPtUncertDaughter2,  "fPtUncertDaughter2/F");
  gTree->Branch("fTPCRefitDaughter2",  &fTPCRefitDaughter2,  "fTPCRefitDaughter2/I");
  gTree->Branch("fITSRefitDaughter2",  &fITSRefitDaughter2,  "fITSRefitDaughter2/I");
  gTree->Branch("fPropDCADaughter2",   &fPropDCADaughter2,   "fPropDCADaughter2/I");
  gTree->Branch("fdEdxSigmaPion",      &fdEdxSigmaPion,      "fdEdxSigmaPion/F");
  gTree->Branch("fdEdxSigmaDeuteron",  &fdEdxSigmaDeuteron,  "fdEdxSigmaDeuteron/F");
  gTree->Branch("fdEdxSigmaTriton",    &fdEdxSigmaTriton,    "fdEdxSigmaTriton/F");
  gTree->Branch("fdEdxSigmaAlpha",     &fdEdxSigmaAlpha,     "fdEdxSigmaAlpha/F");
  // _________________________________________________ //
  // __ generated tree __ //
  fTreeGen = new TTree("fTreeGen", "fTreeGen");  
  fTreeGen->Branch("fPeriod",         &fPeriod,              "fPeriod/I");
  fTreeGen->Branch("frunnumber",      &frunnumber,           "frunnumber/I");
  fTreeGen->Branch("fCentrality",     &fCentrality,          "fCentrality/I");
  fTreeGen->Branch("fMagneticField",  &fMagneticField,       "fMagneticField/I");
  fTreeGen->Branch("fTrigMB",         &fTrigMB,              "fTrigMB/I");
  fTreeGen->Branch("fTrigHMV0",       &fTrigHMV0,            "fTrigHMV0/I");
  fTreeGen->Branch("fTrigHMSPD",      &fTrigHMSPD,           "fTrigHMSPD/I");
  fTreeGen->Branch("fTrigHNU",        &fTrigHNU,             "fTrigHNU/I");
  fTreeGen->Branch("fTrigHQU",        &fTrigHQU,             "fTrigHQU/I");
  fTreeGen->Branch("fPDGMother",      &fPDGMother,           "fPDGMother/I");
  fTreeGen->Branch("fChargeMother",   &fChargeMother,        "fChargeMother/I");
  fTreeGen->Branch("fmMother",        &fmMother,             "fmMother/F");
  fTreeGen->Branch("fptMother",       &fptMother,            "fptMother/F");
  fTreeGen->Branch("fyMother",        &fyMother,             "fyMother/F");
  fTreeGen->Branch("fctMother",       &fctMother,            "fctMother/F");
  fTreeGen->Branch("fmSubMother",     &fmSubMother,          "fmSubMother/F");
  fTreeGen->Branch("fptSubMother",    &fptSubMother,         "fptSubMother/F");
  fTreeGen->Branch("fySubMother",     &fySubMother,          "fySubMother/F");
  fTreeGen->Branch("fctSubMother",    &fctSubMother,         "fctSubMother/F");
  fTreeGen->Branch("fPDGMother",      &fPDGMother,           "fPDGMother/I");
  fTreeGen->Branch("fChargeMother",   &fChargeMother,        "fChargeMother/I");
  fTreeGen->Branch("fptDaughter",     &fptDaughter,          "fptDaughter/F");
  fTreeGen->Branch("fyDaughter",      &fyDaughter,           "fyDaughter/F");
  fTreeGen->Branch("fptDaughter1",    &fptDaughter1,         "fptDaughter1/F");
  fTreeGen->Branch("fyDaughter1",     &fyDaughter1,          "fyDaughter1/F");
  fTreeGen->Branch("fptDaughter2",    &fptDaughter2,         "fptDaughter2/F");
  fTreeGen->Branch("fyDaughter2",     &fyDaughter2,          "fyDaughter2/F");
  fTreeGen->Branch("fptDaughter3",    &fptDaughter3,         "fptDaughter3/F");
  fTreeGen->Branch("fyDaughter3",     &fyDaughter3,          "fyDaughter3/F");

  // __ generated tree __ //
  gTreeGen = new TTree("gTreeGen", "gTreeGen");  
  gTreeGen->Branch("fPeriod",         &fPeriod,              "fPeriod/I");
  gTreeGen->Branch("frunnumber",      &frunnumber,           "frunnumber/I");
  gTreeGen->Branch("fCentrality",     &fCentrality,          "fCentrality/I");
  gTreeGen->Branch("fMagneticField",  &fMagneticField,       "fMagneticField/I");
  gTreeGen->Branch("fTrigMB",         &fTrigMB,              "fTrigMB/I");
  gTreeGen->Branch("fTrigHMV0",       &fTrigHMV0,            "fTrigHMV0/I");
  gTreeGen->Branch("fTrigHMSPD",      &fTrigHMSPD,           "fTrigHMSPD/I");
  gTreeGen->Branch("fTrigHNU",        &fTrigHNU,             "fTrigHNU/I");
  gTreeGen->Branch("fTrigHQU",        &fTrigHQU,             "fTrigHQU/I");
  gTreeGen->Branch("fPDGMother",      &fPDGMother,           "fPDGMother/I");
  gTreeGen->Branch("fChargeMother",   &fChargeMother,        "fChargeMother/I");
  gTreeGen->Branch("fmSubMother",     &fmSubMother,          "fmSubMother/F");
  gTreeGen->Branch("fptSubMother",    &fptSubMother,         "fptSubMother/F");
  gTreeGen->Branch("fySubMother",     &fySubMother,          "fySubMother/F");
  gTreeGen->Branch("fctSubMother",    &fctSubMother,         "fctSubMother/F");
  gTreeGen->Branch("fPDGMother",      &fPDGMother,           "fPDGMother/I");
  gTreeGen->Branch("fChargeMother",   &fChargeMother,        "fChargeMother/I");
  gTreeGen->Branch("fptDaughter",     &fptDaughter,          "fptDaughter/F");
  gTreeGen->Branch("fyDaughter",      &fyDaughter,           "fyDaughter/F");
  gTreeGen->Branch("fptDaughter1",    &fptDaughter1,         "fptDaughter1/F");
  gTreeGen->Branch("fyDaughter1",     &fyDaughter1,          "fyDaughter1/F");
  gTreeGen->Branch("fptDaughter2",    &fptDaughter2,         "fptDaughter2/F");
  gTreeGen->Branch("fyDaughter2",     &fyDaughter2,          "fyDaughter2/F");

  PostData(1, fHistogramList);
  PostData(2, fTree);
  PostData(3, gTree);
  PostData(4, fTreeGen);
  PostData(5, gTreeGen);


  //*******--------******** Info ********--------*******//
  //  EParticleType :                                   //
  //  kElectron = 0, kMuon = 1, kPion = 2, kKaon = 3,   //  
  //  kProton = 4, kDeuteron = 5, kTriton = 6, kHe3 = 7 //
  //  kAlpha = 8, kPhoton = 9, kPi0 = 10, kNeutron = 11 //
  //  kKaon0 = 12, kEleCon = 13, kUnknown = 14          //
  //                                                    //
  //  Masses: 4Li: 3.74958 GeV/c^2                      //
  //          4He: 3.727379 GeV/c^2                     //
  //          3He: 2.80923 GeV/c^2                      //
  //            t: 2.80925 GeV/c^2                      //
  //            d: 1.875613 GeV/c^2                     //
  //            p: 0.93827 GeV/c^2                      //
  //           pi: 0.13957 GeV/c^2                      //
  //          3LH: 2.99131 GeV/c^2                      //
  //          4LH: 3.931   GeV/c^2                      //
  //         4LHe: 3.929   GeV/c^2                      //
  //         5LHe: 4.841   GeV/c^2                      //
  //         4LLH: 4.106   GeV/c^2                      //
  //*******-------********--------********-------*******//

}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::UserExec(Option_t *) {
  //__ MC __ //
  fMCtrue = kTRUE;
  AliMCEventHandler* mcEventHandler = dynamic_cast<AliMCEventHandler*>
    (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!mcEventHandler) {
    fMCtrue = kFALSE;
  }
  mcEvent = 0x0;
  if (mcEventHandler) mcEvent = mcEventHandler->MCEvent();
  if (!mcEvent) {
    if (fMCtrue) return;
  }
  if (fMCtrue) {
    fStack = mcEvent->Stack();
    if (!fStack) return;
  }
  // __ Data __ //
  fESDevent = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESDevent || !fPID) {
    return;
  } 
  
  fHistNumEvents->Fill(0);

  // __ runnumber & period __ //
  fPeriod    = fESDevent->GetPeriodNumber();
  frunnumber = fESDevent->GetRunNumber();

  // __ EventCuts & TimeRangeCut __ //
  if(!fMCtrue && (frunnumber == 297219 || frunnumber == 297194 || frunnumber == 297029
		  || frunnumber == 296890 || frunnumber == 296849 || frunnumber == 296750
		  || frunnumber == 296749 || frunnumber == 297481)) fEventCuts.UseTimeRangeCut();
  
  // __ classify events __ //
  if(!fMCtrue && fEventCuts.AcceptEvent(fESDevent)) feventclass = 1;
  else if(fMCtrue) feventclass = 1;
  else feventclass = 0;

  // __ Number of acc events __ //
  fHistNumEvents->Fill(1);

  // __ centrality, 0 = V0M __ //
  Float_t centrality = -1;
  AliMultSelection *fMultSelection = (AliMultSelection*)fESDevent->FindListObject("MultSelection");
  centrality = fMultSelection->GetMultiplicityPercentile("V0M");
  fCentrality = centrality;

  // __ Trigger Selection __ //
  fEventCuts.OverrideAutomaticTriggerSelection(fTriggerMask);
  AliAnalysisTaskDoubleHypNucTree::TriggerSelection();

  //__ SetBetheBlochParams __ //
  AliAnalysisTaskDoubleHypNucTree::SetBetheBlochParams(frunnumber);   

  // __ MagneticField __ //
  fMagneticField  = fESDevent->GetMagneticField();

  // __ MC Settings __ //
  if(fMCtrue){
    // __ each track is checked for 4LLH as mother __ //
    noCombinatoricsMC = 1;
    // __ only pion tracks are checked for 4LLH as mother __ //
    lessCombinatoricsMC = 0; 
  }

  // __ setting up vertex reconstruction __ //
  // 1 = reconstruction with coord of prim vtx __ //
  // 2 = reconstruction with coord of tert vtx __ //
  fvariante = 1;
  // __ copy prim vtx __ //
  vertex = fESDevent->GetPrimaryVertexSPD();
  // __ coord of prim vtx __ //
  PrimVertex[0] = vertex->GetX();
  PrimVertex[1] = vertex->GetY();
  PrimVertex[2] = vertex->GetZ();

  // __ only PID __ //
  if(fPIDCheckOnly){
    AliAnalysisTaskDoubleHypNucTree::dEdxCheck();
  }// __ Analysis __ //
  else{
    if(fMCtrue) AliAnalysisTaskDoubleHypNucTree::MCGenerated();
    AliAnalysisTaskDoubleHypNucTree::FindTracks();    
  }

  PostData(1, fHistogramList);
  PostData(2, fTree);
  PostData(3, gTree);
  PostData(4, fTreeGen);
  PostData(5, gTreeGen);
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::dEdxCheck(){
  AliESDtrackCuts* trackCutsPid = new AliESDtrackCuts("trackCutsPid", "trackCutsPid");
  trackCutsPid = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  trackCutsPid->SetEtaRange(-0.9,0.9);
  for (Int_t itrack = 0; itrack < fESDevent->GetNumberOfTracks(); itrack++) {
    AliESDtrack* track = fESDevent->GetTrack(itrack);
    if (!trackCutsPid->AcceptTrack(track)) continue;
    Double_t momentum = track->GetInnerParam()->GetP();
    fHistdEdx->Fill(momentum * track->GetSign(), track->GetTPCsignal());
  }
  delete trackCutsPid;
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::InitArrays(){

  const Int_t nTracksEvent = fESDevent->GetNumberOfTracks();
  // __ Track Array for identified particles __ //
  // __ track counter __ //
  He3PosCounter   = 0;
  He3NegCounter   = 0;
  PPosCounter     = 0;
  PNegCounter     = 0;
  PiPosCounter    = 0;
  PiNegCounter    = 0;
  PiPosSecCounter = 0;
  PiNegSecCounter = 0;
  // __ resize vector arrays to number of tracks __ //
  He3PosArray.resize(nTracksEvent);
  He3NegArray.resize(nTracksEvent);
  PPosArray.resize(nTracksEvent);
  PNegArray.resize(nTracksEvent);
  PiPosArray.resize(nTracksEvent);
  PiNegArray.resize(nTracksEvent);
  PiPosSecArray.resize(nTracksEvent);
  PiNegSecArray.resize(nTracksEvent);
  // __ set all elements to zero __ //
  for(int i = 0; i<nTracksEvent;i++){
    He3PosArray[i]   = 0;
    He3NegArray[i]   = 0;
    PPosArray[i]     = 0;
    PNegArray[i]     = 0;
    PiPosArray[i]    = 0;
    PiNegArray[i]    = 0;
    PiPosSecArray[i] = 0;
    PiNegSecArray[i] = 0;
  }

  //__ Reduced Track Array for 4LLH __ //   
  // __ track counter __ //
  // __ starts at 1 due to array checks in the loop __ //
  He3PosCounterAcc = 1;
  He3NegCounterAcc = 1;
  PPosCounterAcc   = 1;
  PNegCounterAcc   = 1;
  PiPosCounterAcc  = 1;
  PiNegCounterAcc  = 1;
  // __ resize vector arrays to number of tracks __ //
  He3PosArrayAcc.resize(nTracksEvent+1);
  He3NegArrayAcc.resize(nTracksEvent+1);
  PPosArrayAcc.resize(nTracksEvent+1);
  PNegArrayAcc.resize(nTracksEvent+1);
  PiPosArrayAcc.resize(nTracksEvent+1);
  PiNegArrayAcc.resize(nTracksEvent+1);
  // __ set all elements to zero __ //
  for(int i = 0; i<nTracksEvent+1;i++){
    He3PosArrayAcc[i] = 0;
    He3NegArrayAcc[i] = 0;
    PPosArrayAcc[i]   = 0;
    PNegArrayAcc[i]   = 0;
    PiPosArrayAcc[i]  = 0;
    PiNegArrayAcc[i]  = 0;
  }

}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::FindTracks(){

  const Int_t nTracksEvent = fESDevent->GetNumberOfTracks();

  // __ init track arrays for identified particles __ //
  AliAnalysisTaskDoubleHypNucTree::InitArrays();

  // __ track loop for PID __ //
  for (Int_t qTracks = 0; qTracks < nTracksEvent; qTracks++) {

    AliESDtrack *trackq = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(qTracks));

    if (!trackq->GetInnerParam()) continue;

    // __ fill energy loss plot __ //
    fHistdEdx->Fill(trackq->GetInnerParam()->GetP()*trackq->GetSign(), trackq->GetTPCsignal());

    // __ positive tracks __ //
    if(trackq->GetSign() > 0){      
      if(TMath::Abs(AliAnalysisTaskDoubleHypNucTree::Bethe(*trackq, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) <= 4) {
	if (!trackCutsNuclei->AcceptTrack(trackq)) continue; 
	He3PosArray[He3PosCounter] = qTracks;	
	He3PosCounter++;	
      }
      if(TMath::Abs(AliAnalysisTaskDoubleHypNucTree::Bethe(*trackq, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT)) <= 4) {
	if (!trackCutsSoft->AcceptTrack(trackq)) continue; 
	PPosArray[PPosCounter] = qTracks;
	PPosCounter++;	
      }
      if(TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kPion)) <= 3) {
	if (!trackCutsStrong->AcceptTrack(trackq)) continue; 
	PiPosArray[PiPosCounter] = qTracks;
	PiPosCounter++;
	PiPosSecArray[PiPosSecCounter] = qTracks;
	PiPosSecCounter++;	
      }
    }
    // __ negative tracks __ //
    if(trackq->GetSign() < 0){      
      if(TMath::Abs(AliAnalysisTaskDoubleHypNucTree::Bethe(*trackq, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) <= 4) {
	if (!trackCutsNuclei->AcceptTrack(trackq)) continue; 
	He3NegArray[He3NegCounter] = qTracks;
	He3NegCounter++;
      }
      if(TMath::Abs(AliAnalysisTaskDoubleHypNucTree::Bethe(*trackq, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT)) <= 4) {
	if (!trackCutsSoft->AcceptTrack(trackq)) continue; 
	PNegArray[PNegCounter] = qTracks;
	PNegCounter++;	
      }
      if(TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kPion)) <= 3) {
	if (!trackCutsStrong->AcceptTrack(trackq)) continue; 
	PiNegArray[PiNegCounter] = qTracks;
	PiNegCounter++;
	PiNegSecArray[PiNegSecCounter] = qTracks;
	PiNegSecCounter++;	
      }
    }
  }
  
  //__ positive 4LHe __ //
  AliAnalysisTaskDoubleHypNucTree::CalcPosDaughterHypNuc();
  AliAnalysisTaskDoubleHypNucTree::ResetVals("");
  // __ positive 4LLH __ //
  AliAnalysisTaskDoubleHypNucTree::CalcPosMotherHypNuc();
  AliAnalysisTaskDoubleHypNucTree::ResetVals("");
  //__ negative 4LHe __ //
  AliAnalysisTaskDoubleHypNucTree::CalcNegDaughterHypNuc();
  AliAnalysisTaskDoubleHypNucTree::ResetVals("");
  // __ negative 4LLH __ //
  AliAnalysisTaskDoubleHypNucTree::CalcNegMotherHypNuc();
  AliAnalysisTaskDoubleHypNucTree::ResetVals("");
  // __ reset __ //
  ResetVals("Event");
  
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::CalcPosDaughterHypNuc(){

  // ____________________________ first track _______________________________ //
  for (const Int_t &He3PosTracks : He3PosArray) {

    if(He3PosTracks == 0) continue;
    track1 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(He3PosTracks));

    ptot1 = track1->GetInnerParam()->GetP();
    sign1 = track1->GetSign();    

    // __ MC part __ //
    if(fMCtrue && noCombinatoricsMC){
	
      label1 = track1->GetLabel();
		  
      labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));
      labelGrandMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother1));

      AliMCParticle *part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(label1))->Particle());
      if(part->PdgCode() != fgkPdgCode[kPDGHelium3])continue;
      if(part)delete part;
	      
      ParticleGrandMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());

      if(ParticleGrandMother1->PdgCode() != fgkPdgCode[kPDGHyperHelium4]) continue;
      if(ParticleGrandMother1) delete ParticleGrandMother1;
    }
    // ____________________________ second track _______________________________ //
    for (const Int_t &PPosTracks : PPosArray) {

      if(PPosTracks == 0) continue;

      if(He3PosTracks == PPosTracks) continue;

      track2 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(PPosTracks));

      ptot2 = track2->GetInnerParam()->GetP();
      sign2 = track2->GetSign();          

      // __ MC part __ //
      if(fMCtrue && noCombinatoricsMC){
	
	label2 = track2->GetLabel();
		  
	labelMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label2));
	labelGrandMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother2));

	AliMCParticle *part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(label2))->Particle());
	if(part->PdgCode() != fgkPdgCode[kPDGProton]) continue;
	if(part) delete part;
	      
	ParticleGrandMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother2))->Particle());

	if(ParticleGrandMother2->PdgCode() != fgkPdgCode[kPDGHyperHelium4]) continue;
	if(ParticleGrandMother2) delete ParticleGrandMother2;
      }

      // __ dca rejection bet tracks __ //
      if(TMath::Abs(track2->GetDCA(track1, fMagneticField, xthiss, xpp)) > 1.2)  continue;

      // ____________________________ third track _______________________________ //
      for (const Int_t &PiNegTracks : PiNegArray) {

	if(PiNegTracks == 0) continue;

	// __ reject using same track twice __ //
	if (PPosTracks == PiNegTracks || He3PosTracks == PiNegTracks) continue;

	track3 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(PiNegTracks));        	
	
	ptot3 = track3->GetInnerParam()->GetP();
	sign3 = track3->GetSign();

	// __ MC part __ //
	if(fMCtrue && (lessCombinatoricsMC || noCombinatoricsMC)){
	
	  label3 = track3->GetLabel();
		  
	  labelMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label3));
	  labelGrandMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother3));

	  AliMCParticle *part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(label3))->Particle());
	  if(part->PdgCode() != fgkPdgCode[kPDGPionMinus]) continue;
	  if(part) delete part;
	      
	  ParticleGrandMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());

	  if(ParticleGrandMother3->PdgCode() != fgkPdgCode[kPDGHyperHelium4]) continue;
	  if(ParticleGrandMother3) delete ParticleGrandMother3;
	}

	// __ dca rejection betw tracks __ //
	if(TMath::Abs(track2->GetDCA(track1, fMagneticField, xthiss, xpp)) > 1.2
	   || TMath::Abs(track3->GetDCA(track1, fMagneticField, xthiss, xpp)) > 1.2
	   || TMath::Abs(track3->GetDCA(track2, fMagneticField, xthiss, xpp)) > 1.2) continue;

	
	//______ daughter hypernucleus mass _______ //
	sublorentzsum = new TLorentzVector(0.,0.,0.,0.);
	particle1     = new TLorentzVector(0.,0.,0.,0.);
	particle2     = new TLorentzVector(0.,0.,0.,0.);
	particle3     = new TLorentzVector(0.,0.,0.,0.);
	h             = new TVector3(0., 0., 0.);
	particle1->SetXYZM(2.*track1->Px(), 2.*track1->Py(), 2.*track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
	particle2->SetXYZM(   track2->Px(),    track2->Py(),    track2->Pz(), AliPID::ParticleMass(AliPID::kProton));
	particle3->SetXYZM(   track3->Px(),    track3->Py(),    track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
	sublorentzsum->SetXYZM(0.,0.,0.,0.);
	*sublorentzsum = *particle1 + *particle2 + *particle3;

	if(sublorentzsum->M() > 3.99 || sublorentzsum->M() < 3.87) {//3.90, 3.96
	  if(sublorentzsum) delete sublorentzsum;
	  if(particle1)     delete particle1;
	  if(particle2)     delete particle2;
	  if(particle3)     delete particle3;
	  if(h)             delete h;
	  AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	  continue;
	}

	// __ MC part __ //
	if(fMCtrue){
	  label1 = 0;
	  label2 = 0;
	  label3 = 0;
		  
	  label1 = track1->GetLabel();
	  label2 = track2->GetLabel();
	  label3 = track3->GetLabel();	
		  
	  labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));
	  labelGrandMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother1));
	      
	  labelMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label2));
	  labelGrandMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother2));
	      
	  labelMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label3));
	  labelGrandMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother3));

	  ParticleGrandMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	  ParticleGrandMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother2))->Particle());
	  ParticleGrandMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());
	}
	
	// _________________________________ //
	Double_t impar[2];
	// __ init AliVertexerTracks __ //
	tertvertexer = new AliVertexerTracks(fESDevent->GetMagneticField());
	// __ get prim vtx __ //
	primVertex = new AliESDVertex(*vertex);
	// __ tert vtx __ //
	trkArray = new TObjArray(3);
	trkArray->AddAt(track1,0);
	trkArray->AddAt(track2,1);
	trkArray->AddAt(track3,2);
	// __ start at prim vtx __ //
	tertvertexer->SetVtxStart(primVertex);
	tertVertex = (AliESDVertex*)tertvertexer->VertexForSelectedESDTracks(trkArray);
	if(trkArray) delete trkArray;
	TertVertex[0] = tertVertex->GetX();
	TertVertex[1] = tertVertex->GetY();
	TertVertex[2] = tertVertex->GetZ();
	// _________________________________ //
	// __ propagate tracks and get impact parameters __ //
	if(track1->PropagateToDCA(tertVertex, fMagneticField, 10, impar)) fPropDCADaughter = 1;
	fImParDaughter  = impar[0];
	fImParzDaughter = impar[1];
	if(track2->PropagateToDCA(tertVertex, fMagneticField, 10, impar)) fPropDCADaughter1 = 1;
	fImParDaughter1  = impar[0];
	fImParzDaughter1 = impar[1];
	if(track3->PropagateToDCA(tertVertex, fMagneticField, 10, impar)) fPropDCADaughter2 = 1;
	fImParDaughter2  = impar[0];
	fImParzDaughter2 = impar[1];
	if(tertVertex) delete tertVertex;    
	// _________________________________ //	
	if (tertvertexer) delete tertvertexer;
	if (primVertex)   delete primVertex;
	// _________________________________________________ //
	//********** 4LHe **********
	// _________________________________________________ //
	// __ daughter hypernucleus cos(PA) to prim vtx __ //
	dd[0] = PrimVertex[0] - TertVertex[0];
	dd[1] = PrimVertex[1] - TertVertex[1];
	dd[2] = PrimVertex[2] - TertVertex[2];
	h->SetXYZ(-dd[0],-dd[1],-dd[2]);

	fSubPA = TMath::Cos(sublorentzsum->Angle(*h));
	if(fSubPA < 0.99) {
	  if(sublorentzsum) delete sublorentzsum;
	  if(particle1)     delete particle1;
	  if(particle2)     delete particle2;
	  if(particle3)     delete particle3;
	  if(h)             delete h;
	  AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	  continue;
	}
	// _________________________________________________ //
	// __ trigger __ //
	fTrigMB    = MB;		
	fTrigHMV0  = HMV0;
	fTrigHMSPD = HMSPD;
	fTrigHNU   = HNU;
	fTrigHQU   = HQU;
	// _________________________________________________ //	   
	// __ coord of prim vtx __ //
	fPrimVertexX = PrimVertex[0];
	fPrimVertexY = PrimVertex[1];
	fPrimVertexZ = PrimVertex[2];
	// _________________________________________________ //
	// __ coord of tert vtx __ //
	// __ reco starts at prim vtx __ //
	fTertVertexX = TertVertex[0];
	fTertVertexY = TertVertex[1];
	fTertVertexZ = TertVertex[2];
	// _________________________________________________ //
	// __ dca between tracks __ //
	fDCA3B1 = TMath::Abs(track2->GetDCA(track1,  fMagneticField, xthiss, xpp));
	fDCA3B2 = TMath::Abs(track3->GetDCA(track1,  fMagneticField, xthiss, xpp));
	fDCA3B3 = TMath::Abs(track3->GetDCA(track2,  fMagneticField, xthiss, xpp));	      
	// __ cuts from 04/2020 __ //
	fDcaDaughtero  = track1->GetD(fPrimVertexX, fPrimVertexY, fMagneticField);
	fDcaDaughter1o = track2->GetD(fPrimVertexX, fPrimVertexY, fMagneticField);
	fDcaDaughter2o = track3->GetD(fPrimVertexX, fPrimVertexY, fMagneticField);
	// _________________________________________________ //
	// __ track information __ //
	AliAnalysisTaskDoubleHypNucTree::SetDaughterInformation("4LHe");		       
	// _________________________________________________ //	  	
	// __ daughter hypernucleus information __ //
	fPDGMother       = fgkPdgCode[kPDGHyperHelium4];
	fChargeMother    = 2;
	fmSubMother      = sublorentzsum->M();
	fESubMother      = sublorentzsum->E();
	fpxSubMother     = sublorentzsum->Px();
	fpySubMother     = sublorentzsum->Py();
	fpzSubMother     = sublorentzsum->Pz();
	fptSubMother     = sublorentzsum->Pt();
	fpSubMother      = sublorentzsum->P();
	fySubMother      = sublorentzsum->Rapidity();
	hInvMass4LHe->Fill(fmSubMother);	  	
	// _________________________________________________ //
	// __ MC part __ //
	if(fMCtrue){
	  if(ParticleGrandMother1->PdgCode() == fgkPdgCode[kPDGHyperHelium4]
	     && ParticleGrandMother2->PdgCode() == fgkPdgCode[kPDGHyperHelium4]
	     && ParticleGrandMother3->PdgCode() == fgkPdgCode[kPDGHyperHelium4]
	     && labelGrandMother1 == labelGrandMother2
	     && labelGrandMother2 == labelGrandMother3){
	    if(TMath::Abs(label1) == TMath::Abs(GetLabel(labelGrandMother1, fgkPdgCode[kPDGHelium3]))
	       && TMath::Abs(label2) == TMath::Abs(GetLabel(labelGrandMother1, fgkPdgCode[kPDGProton]))
	       && TMath::Abs(label3) == TMath::Abs(GetLabel(labelGrandMother1, fgkPdgCode[kPDGPionMinus]))){
	      mctruth = 1;		      	     
	    }	  
	  }
	}
	// __ fill tree __ //
	gTree->Fill();
	// __ fill reduced arrays __ //
	if(He3PosArrayAcc[He3PosCounterAcc-1] != He3PosTracks){
	  He3PosArrayAcc[He3PosCounterAcc] = He3PosTracks;
	  He3PosCounterAcc++;
	}
	if(PPosArrayAcc[PPosCounterAcc-1] != PPosTracks){
	  PPosArrayAcc[PPosCounterAcc] = PPosTracks;
	  PPosCounterAcc++;
	}
	if(PiNegArrayAcc[PiNegCounterAcc-1] != PiNegTracks){
	  PiNegArrayAcc[PiNegCounterAcc] = PiNegTracks;
	  PiNegCounterAcc++;
	}
	// __ reset __ //
	AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	if(sublorentzsum)  delete sublorentzsum;
	if(particle1)      delete particle1;
	if(particle2)      delete particle2;
	if(particle3)      delete particle3;
	if(h)              delete h;	
      }
    }
  } 
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::CalcNegDaughterHypNuc(){

  // ____________________________ first track _______________________________ //
  for (const Int_t &He3NegTracks : He3NegArray) {

    if(He3NegTracks == 0) continue;
    
    track1 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(He3NegTracks));

    ptot1 = track1->GetInnerParam()->GetP();
    sign1 = track1->GetSign();

    // __ MC part __ //
    if(fMCtrue && noCombinatoricsMC){
	
      label1 = track1->GetLabel();
		  
      labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));
      labelGrandMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother1));

      AliMCParticle *part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(label1))->Particle());
      if(part->PdgCode() != fgkPdgCode[kPDGAntiHelium3])continue;
      if(part)delete part;
	      
      ParticleGrandMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());

      if(ParticleGrandMother1->PdgCode() != fgkPdgCode[kPDGAntiHyperHelium4]) continue;
      if(ParticleGrandMother1) delete ParticleGrandMother1;
    }
    // ____________________________ second track _______________________________ //
    for (const Int_t &PNegTracks : PNegArray) {

      if(PNegTracks == 0) continue;

      if(He3NegTracks == PNegTracks) continue;

      track2 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(PNegTracks));

      ptot2 = track2->GetInnerParam()->GetP();
      sign2 = track2->GetSign();          

      // __ MC part __ //
      if(fMCtrue && noCombinatoricsMC){
	
	label2 = track2->GetLabel();
		  
	labelMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label2));
	labelGrandMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother2));

	AliMCParticle *part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(label2))->Particle());
	if(part->PdgCode() != fgkPdgCode[kPDGAntiProton]) continue;
	if(part) delete part;
	      
	ParticleGrandMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother2))->Particle());

	if(ParticleGrandMother2->PdgCode() != fgkPdgCode[kPDGAntiHyperHelium4]) continue;
	if(ParticleGrandMother2) delete ParticleGrandMother2;
      }

      // __ dca rejection bet tracks __ //
      if(TMath::Abs(track2->GetDCA(track1, fMagneticField, xthiss, xpp)) > 1.2)  continue;

      // ____________________________ third track _______________________________ //
      for (const Int_t &PiPosTracks : PiPosArray) {

	if(PiPosTracks == 0) continue;

	// __ reject using same track twice __ //
	if (PNegTracks == PiPosTracks || He3NegTracks == PiPosTracks) continue;

	track3 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(PiPosTracks));        	
	
	ptot3 = track3->GetInnerParam()->GetP();
	sign3 = track3->GetSign();

	// __ MC part __ //
	if(fMCtrue && (lessCombinatoricsMC || noCombinatoricsMC)){
	
	  label3 = track3->GetLabel();
		  
	  labelMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label3));
	  labelGrandMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother3));

	  AliMCParticle *part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(label3))->Particle());
	  if(part->PdgCode() != fgkPdgCode[kPDGPionPlus]) continue;
	  if(part) delete part;
	      
	  ParticleGrandMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());

	  if(ParticleGrandMother3->PdgCode() != fgkPdgCode[kPDGAntiHyperHelium4]) continue;
	  if(ParticleGrandMother3) delete ParticleGrandMother3;
	}
	
	// __ dca rejection betw tracks __ //
	if(TMath::Abs(track2->GetDCA(track1, fMagneticField, xthiss, xpp)) > 1.2
	   || TMath::Abs(track3->GetDCA(track1, fMagneticField, xthiss, xpp)) > 1.2
	   || TMath::Abs(track3->GetDCA(track2, fMagneticField, xthiss, xpp)) > 1.2) continue;

	
	//______ daughter hypernucleus mass _______ //
	sublorentzsum = new TLorentzVector(0.,0.,0.,0.);
	particle1     = new TLorentzVector(0.,0.,0.,0.);
	particle2     = new TLorentzVector(0.,0.,0.,0.);
	particle3     = new TLorentzVector(0.,0.,0.,0.);
	h             = new TVector3(0., 0., 0.);
	particle1->SetXYZM(2.*track1->Px(), 2.*track1->Py(), 2.*track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
	particle2->SetXYZM(   track2->Px(),    track2->Py(),    track2->Pz(), AliPID::ParticleMass(AliPID::kProton));
	particle3->SetXYZM(   track3->Px(),    track3->Py(),    track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
	sublorentzsum->SetXYZM(0.,0.,0.,0.);
	*sublorentzsum = *particle1 + *particle2 + *particle3;

	if(sublorentzsum->M() > 3.99 || sublorentzsum->M() < 3.87) {
	  if(sublorentzsum) delete sublorentzsum;
	  if(particle1)     delete particle1;
	  if(particle2)     delete particle2;
	  if(particle3)     delete particle3;
	  if(h)             delete h;
	  AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	  continue;
	}

	// __ MC part __ //
	if(fMCtrue){
	  label1 = 0;
	  label2 = 0;
	  label3 = 0;
		  
	  label1 = track1->GetLabel();
	  label2 = track2->GetLabel();
	  label3 = track3->GetLabel();	
		  
	  labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));
	  labelGrandMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother1));
	      
	  labelMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label2));
	  labelGrandMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother2));
	      
	  labelMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label3));
	  labelGrandMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother3));

	  ParticleGrandMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	  ParticleGrandMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother2))->Particle());
	  ParticleGrandMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());
	}
	
	// _________________________________ //
	Double_t impar[2];
	// __ init AliVertexerTracks __ //
	tertvertexer = new AliVertexerTracks(fESDevent->GetMagneticField());
	// __ get prim vtx __ //
	primVertex = new AliESDVertex(*vertex);
	// __ tert vtx __ //
	trkArray = new TObjArray(3);
	trkArray->AddAt(track1,0);
	trkArray->AddAt(track2,1);
	trkArray->AddAt(track3,2);
	// __ start at prim vtx __ //
	tertvertexer->SetVtxStart(primVertex);
	tertVertex = (AliESDVertex*)tertvertexer->VertexForSelectedESDTracks(trkArray);
	if(trkArray) delete trkArray;
	TertVertex[0] = tertVertex->GetX();
	TertVertex[1] = tertVertex->GetY();
	TertVertex[2] = tertVertex->GetZ();
	// _________________________________ //
	// __ propagate tracks and get impact parameters __ //
	if(track1->PropagateToDCA(tertVertex, fMagneticField, 10, impar)) fPropDCADaughter = 1;
	fImParDaughter  = impar[0];
	fImParzDaughter = impar[1];
	if(track2->PropagateToDCA(tertVertex, fMagneticField, 10, impar)) fPropDCADaughter1 = 1;
	fImParDaughter1  = impar[0];
	fImParzDaughter1 = impar[1];
	if(track3->PropagateToDCA(tertVertex, fMagneticField, 10, impar)) fPropDCADaughter2 = 1;
	fImParDaughter2  = impar[0];
	fImParzDaughter2 = impar[1];
	if(tertVertex) delete tertVertex;    
	// _________________________________ //
	if (tertvertexer) delete tertvertexer;
	if (primVertex)   delete primVertex;
	// _________________________________________________ //
	//********** 4LHe **********
	// _________________________________________________ //
	// __ daughter hypernucleus cos(PA) to prim vtx __ //
	dd[0] = PrimVertex[0] - TertVertex[0];
	dd[1] = PrimVertex[1] - TertVertex[1];
	dd[2] = PrimVertex[2] - TertVertex[2];
	h->SetXYZ(-dd[0],-dd[1],-dd[2]);

	fSubPA = TMath::Cos(sublorentzsum->Angle(*h));
	if(fSubPA < 0.99) {
	  if(sublorentzsum) delete sublorentzsum;
	  if(particle1)     delete particle1;
	  if(particle2)     delete particle2;
	  if(particle3)     delete particle3;
	  if(h)             delete h;
	  AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	  continue;
	}
	// _________________________________________________ //
	// __ trigger __ //
	fTrigMB    = MB;		
	fTrigHMV0  = HMV0;
	fTrigHMSPD = HMSPD;
	fTrigHNU   = HNU;
	fTrigHQU   = HQU;
	// _________________________________________________ //	   
	// __ coord of prim vtx __ //
	fPrimVertexX = PrimVertex[0];
	fPrimVertexY = PrimVertex[1];
	fPrimVertexZ = PrimVertex[2];
	// _________________________________________________ //
	// __ coord of tert vtx __ //
	// __ reco starts at prim vtx __ //
	fTertVertexX = TertVertex[0];
	fTertVertexY = TertVertex[1];
	fTertVertexZ = TertVertex[2];
	// _________________________________________________ //
	// __ dca between tracks __ //
	fDCA3B1 = TMath::Abs(track2->GetDCA(track1,  fMagneticField, xthiss, xpp));
	fDCA3B2 = TMath::Abs(track3->GetDCA(track1,  fMagneticField, xthiss, xpp));
	fDCA3B3 = TMath::Abs(track3->GetDCA(track2,  fMagneticField, xthiss, xpp));	      
	// __ cuts from 04/2020 __ //
	fDcaDaughtero  = track1->GetD(fPrimVertexX, fPrimVertexY, fMagneticField);
	fDcaDaughter1o = track2->GetD(fPrimVertexX, fPrimVertexY, fMagneticField);
	fDcaDaughter2o = track3->GetD(fPrimVertexX, fPrimVertexY, fMagneticField);
	// _________________________________________________ //
	// __ track information __ //
	AliAnalysisTaskDoubleHypNucTree::SetDaughterInformation("4LHe");		       
	// _________________________________________________ //	  	
	// __ daughter hypernucleus information __ //
	fPDGMother       = fgkPdgCode[kPDGAntiHyperHelium4];
	fChargeMother    = -2;
	fmSubMother      = sublorentzsum->M();
	fESubMother      = sublorentzsum->E();
	fpxSubMother     = sublorentzsum->Px();
	fpySubMother     = sublorentzsum->Py();
	fpzSubMother     = sublorentzsum->Pz();
	fptSubMother     = sublorentzsum->Pt();
	fpSubMother      = sublorentzsum->P();
	fySubMother      = sublorentzsum->Rapidity();
	hInvMass4LHe->Fill(fmSubMother);	  
	// _________________________________________________ //	
	// __ MC part __ //
	if(fMCtrue){
	  if(ParticleGrandMother1->PdgCode() == fgkPdgCode[kPDGAntiHyperHelium4]
	     && ParticleGrandMother2->PdgCode() == fgkPdgCode[kPDGAntiHyperHelium4]
	     && ParticleGrandMother3->PdgCode() == fgkPdgCode[kPDGAntiHyperHelium4]
	     && labelGrandMother1 == labelGrandMother2
	     && labelGrandMother2 == labelGrandMother3){
	    if(TMath::Abs(label1) == TMath::Abs(GetLabel(labelGrandMother1, fgkPdgCode[kPDGAntiHelium3]))
	       && TMath::Abs(label2) == TMath::Abs(GetLabel(labelGrandMother1, fgkPdgCode[kPDGAntiProton]))
	       && TMath::Abs(label3) == TMath::Abs(GetLabel(labelGrandMother1, fgkPdgCode[kPDGPionPlus]))){
	      mctruth = 1;		      	     
	    }	  
	  }
	}
	// __ fill tree __ //
	gTree->Fill();
	// __ fill reduced arrays __ //
	if(He3NegArrayAcc[He3NegCounterAcc-1] != He3NegTracks){
	  He3NegArrayAcc[He3NegCounterAcc] = He3NegTracks;
	  He3NegCounterAcc++;
	}
	if(PNegArrayAcc[PNegCounterAcc-1] != PNegTracks){
	  PNegArrayAcc[PNegCounterAcc] = PNegTracks;
	  PNegCounterAcc++;
	}
	if(PiPosArrayAcc[PiPosCounterAcc-1] != PiPosTracks){
	  PiPosArrayAcc[PiPosCounterAcc] = PiPosTracks;
	  PiPosCounterAcc++;
	}
	// __ reset __ //
	AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	if(sublorentzsum)  delete sublorentzsum;
	if(particle1)      delete particle1;
	if(particle2)      delete particle2;
	if(particle3)      delete particle3;
	if(h)              delete h;	
      }
    }
  }  
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::CalcPosMotherHypNuc(){
  // ____________________________ first track _______________________________ //
  for (const Int_t &He3PosTracks : He3PosArrayAcc) {

    if(He3PosTracks == 0) continue;

    track1 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(He3PosTracks));
    
    ptot1 = track1->GetInnerParam()->GetP();
    sign1 = track1->GetSign();

    // __ MC part __ //
    if(fMCtrue && noCombinatoricsMC){
	
      label1 = track1->GetLabel();
		  
      labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));
      labelGrandMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother1));

      AliMCParticle *part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(label1))->Particle());
      if(part->PdgCode() != fgkPdgCode[kPDGHelium3])continue;
      if(part)delete part;
	      
      ParticleGrandMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother1))->Particle());

      if(ParticleGrandMother1->PdgCode() != fgkPdgCode[kPDGDoubleHyperHydrogen4]) continue;
      if(ParticleGrandMother1) delete ParticleGrandMother1;
    }

    // ____________________________ second track _______________________________ //
    for (const Int_t &PPosTracks : PPosArrayAcc) {

      if(PPosTracks == 0) continue;
      // __ reject using same track twice __ //
      if(He3PosTracks == PPosTracks) continue;

      track2 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(PPosTracks));

      ptot2 = track2->GetInnerParam()->GetP();
      sign2 = track2->GetSign();       

      // __ MC part __ //
      if(fMCtrue && noCombinatoricsMC){
	
	label2 = track2->GetLabel();
		  
	labelMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label2));
	labelGrandMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother2));

	AliMCParticle *part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(label2))->Particle());
	if(part->PdgCode() != fgkPdgCode[kPDGProton]) continue;
	if(part) delete part;
	      
	ParticleGrandMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother2))->Particle());

	if(ParticleGrandMother2->PdgCode() != fgkPdgCode[kPDGDoubleHyperHydrogen4]) continue;
	if(ParticleGrandMother2) delete ParticleGrandMother2;
      }

      // __ dca rejection bet tracks __ //
      if(TMath::Abs(track2->GetDCA(track1, fMagneticField, xthiss, xpp)) > 1.2)  continue;

      // ____________________________ third track _______________________________ //
      for (const Int_t &PiNegTracks : PiNegArrayAcc) {

	if(PiNegTracks == 0) continue;

	// __ reject using same track twice __ //
	if (PPosTracks == PiNegTracks || He3PosTracks == PiNegTracks) continue;

	track3 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(PiNegTracks));        	
	
	ptot3 = track3->GetInnerParam()->GetP();
	sign3 = track3->GetSign();

	// __ MC part __ //
	if(fMCtrue && (lessCombinatoricsMC || noCombinatoricsMC)){
	
	  label3 = track3->GetLabel();
		  
	  labelMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label3));
	  labelGrandMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother3));

	  AliMCParticle *part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(label3))->Particle());
	  if(part->PdgCode() != fgkPdgCode[kPDGPionMinus]) continue;
	  if(part) delete part;
	      
	  ParticleGrandMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother3))->Particle());

	  if(ParticleGrandMother3->PdgCode() != fgkPdgCode[kPDGDoubleHyperHydrogen4]) continue;
	  if(ParticleGrandMother3) delete ParticleGrandMother3;
	}

	// __ dca rejection betw tracks __ //
	if(TMath::Abs(track2->GetDCA(track1, fMagneticField, xthiss, xpp)) > 1.2
	   || TMath::Abs(track3->GetDCA(track1, fMagneticField, xthiss, xpp)) > 1.2
	   || TMath::Abs(track3->GetDCA(track2, fMagneticField, xthiss, xpp)) > 1.2) continue;

	// ____________________________ fourth track _______________________________ //
	for (const Int_t &PiNegSecTracks : PiNegSecArray) {

	  if(PiNegSecTracks == 0) continue;

	  // __ reject using same track twice __ //
	  if (PPosTracks == PiNegSecTracks || He3PosTracks == PiNegSecTracks || PiNegTracks == PiNegSecTracks) continue;

	  track4 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(PiNegSecTracks));	  	  

	  ptot4 = track4->GetInnerParam()->GetP();
	  sign4 = track4->GetSign();


	  // __ MC part __ //
	  if(fMCtrue && (lessCombinatoricsMC || noCombinatoricsMC)){
	
	    label4 = track4->GetLabel();
		  
	    labelMother4 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label4));
	      
	    ParticleMother4 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother4))->Particle());

	    AliMCParticle *part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(label4))->Particle());
	    if(part->PdgCode() != fgkPdgCode[kPDGPionMinus]) continue;
	    if(part)delete part;

	    if(ParticleMother4->PdgCode() != fgkPdgCode[kPDGDoubleHyperHydrogen4]) continue;

	    if(ParticleMother4) delete ParticleMother4;	    
	  }

	  // _________________________________________________ //
	  lorentzsum     = new TLorentzVector(0.,0.,0.,0.);
	  particle1      = new TLorentzVector(0.,0.,0.,0.);
	  particle2      = new TLorentzVector(0.,0.,0.,0.);
	  particle3      = new TLorentzVector(0.,0.,0.,0.);
	  particle4      = new TLorentzVector(0.,0.,0.,0.);		
	  // _________________________________________________ //
	  particle1->SetXYZM(2.*track1->Px(), 2.*track1->Py(), 2.*track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
	  particle2->SetXYZM(   track2->Px(),    track2->Py(),    track2->Pz(), AliPID::ParticleMass(AliPID::kProton));
	  particle3->SetXYZM(   track3->Px(),    track3->Py(),    track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
	  particle4->SetXYZM(   track4->Px(),    track4->Py(),    track4->Pz(), AliPID::ParticleMass(AliPID::kPion));
	  lorentzsum->SetXYZM(0.,0.,0.,0.);
	  *lorentzsum     = *particle1 + *particle2 + *particle3 + *particle4;
	  // _________________________________________________ //
	  // __ select mother hypernucleus by inv mass __ //
	  if(lorentzsum->M() > 4.16 || lorentzsum->M() < 4.04) {
	    if(lorentzsum)    delete lorentzsum;
	    if(particle1)     delete particle1;
	    if(particle2)     delete particle2;
	    if(particle3)     delete particle3;
	    if(particle4)     delete particle4;
	    AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	    continue;
	  }
	  else{
	    if(lorentzsum)    delete lorentzsum;
	    if(particle1)     delete particle1;
	    if(particle2)     delete particle2;
	    if(particle3)     delete particle3;
	    if(particle4)     delete particle4;
	  }

	  // __ MC part __ //
	  if(fMCtrue){
	    label1 = 0;
	    label2 = 0;
	    label3 = 0;
	    label4 = 0;
		  
	    label1 = track1->GetLabel();
	    label2 = track2->GetLabel();
	    label3 = track3->GetLabel();	
	    label4 = track4->GetLabel();
		  
	    labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));
	    labelGrandMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother1));
	      
	    labelMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label2));
	    labelGrandMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother2));
	      
	    labelMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label3));
	    labelGrandMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother3));
	      
	    labelMother4 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label4));
	    labelGrandMother4 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother4));
	      
	    ParticleMother4 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother4))->Particle());

	    ParticleGrandMother4 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother4))->Particle());
	    ParticleGrandMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother3))->Particle());
	    ParticleGrandMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother2))->Particle());
	    ParticleGrandMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother1))->Particle());
	  }
	  // _________________________________________________ //
	  //********** 4LLH **********
	  // _________________________________________________ //
	  // __ trigger __ //
	  fTrigMB    = MB;		
	  fTrigHMV0  = HMV0;
	  fTrigHMSPD = HMSPD;
	  fTrigHNU   = HNU;
	  fTrigHQU   = HQU;
	  // _________________________________________________ //
	  sublorentzsum  = new TLorentzVector(0.,0.,0.,0.);
	  sublorentzsum2 = new TLorentzVector(0.,0.,0.,0.);
	  lorentzsum     = new TLorentzVector(0.,0.,0.,0.);
	  lorentzsum2    = new TLorentzVector(0.,0.,0.,0.);
	  particle1      = new TLorentzVector(0.,0.,0.,0.);
	  particle2      = new TLorentzVector(0.,0.,0.,0.);
	  particle3      = new TLorentzVector(0.,0.,0.,0.);
	  particle4      = new TLorentzVector(0.,0.,0.,0.);		
	  h              = new TVector3(0., 0., 0.);
	  // _________________________________________________ //
	  particle1->SetXYZM(2.*track1->Px(), 2.*track1->Py(), 2.*track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
	  particle2->SetXYZM(   track2->Px(),    track2->Py(),    track2->Pz(), AliPID::ParticleMass(AliPID::kProton));
	  particle3->SetXYZM(   track3->Px(),    track3->Py(),    track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
	  particle4->SetXYZM(   track4->Px(),    track4->Py(),    track4->Pz(), AliPID::ParticleMass(AliPID::kPion));
	  sublorentzsum->SetXYZM(0.,0.,0.,0.);
	  sublorentzsum2->SetXYZM(0.,0.,0.,0.);
	  lorentzsum->SetXYZM(0.,0.,0.,0.);
	  lorentzsum2->SetXYZM(0.,0.,0.,0.);
	  *lorentzsum     = *particle1 + *particle2 + *particle3 + *particle4;
	  *sublorentzsum  = *particle1 + *particle2 + *particle3;
	  sublorentzsum2->SetXYZM(sublorentzsum->X(), sublorentzsum->Y(), sublorentzsum->Z(), 3.929);	    
	  *lorentzsum2    = *sublorentzsum2 + *particle4;
	  // _________________________________________________ //	   
	  // __ coord of prim vtx __ //
	  fPrimVertexX = PrimVertex[0];
	  fPrimVertexY = PrimVertex[1];
	  fPrimVertexZ = PrimVertex[2];
	  // _________________________________________________ //
	  // __ create sec / tert vtx __ //
	  AliAnalysisTaskDoubleHypNucTree::CreateSecVertex();
	  // _________________________________________________ //
	  // __ coord of sec vtx __ //
	  fSecVertexX = SecVertex[0];
	  fSecVertexY = SecVertex[1];
	  fSecVertexZ = SecVertex[2];
	  // _________________________________________________ //
	  // __ coord of tert vtx __ //
	  // __ reco starts at prim vtx __ //
	  fTertVertexX = TertVertex[0];
	  fTertVertexY = TertVertex[1];
	  fTertVertexZ = TertVertex[2];
	  // __ reco starts at sec vtx __ //
	  fTertVertex2X = TertVertex2[0];
	  fTertVertex2Y = TertVertex2[1];
	  fTertVertex2Z = TertVertex2[2];
	  // _________________________________________________ //
	  // __ dca between tracks __ //
	  fDCA2B  = TMath::Abs(track4->GetDCA(exTrack, fMagneticField, xthiss, xpp));
	  fDCA3B1 = TMath::Abs(track2->GetDCA(track1,  fMagneticField, xthiss, xpp));
	  fDCA3B2 = TMath::Abs(track3->GetDCA(track1,  fMagneticField, xthiss, xpp));
	  fDCA3B3 = TMath::Abs(track3->GetDCA(track2,  fMagneticField, xthiss, xpp));	      
	  // __ cuts from 04/2020 __ //
	  fDCA2Bo        = TMath::Abs(track4->GetD(sublorentzsum->X(), sublorentzsum->Y(), fMagneticField));
	  fDcaDaughtero  = track1->GetD(fPrimVertexX, fPrimVertexY, fMagneticField);
	  fDcaDaughter1o = track2->GetD(fPrimVertexX, fPrimVertexY, fMagneticField);
	  fDcaDaughter2o = track3->GetD(fPrimVertexX, fPrimVertexY, fMagneticField);
	  fDcaDaughter3o = track4->GetD(fPrimVertexX, fPrimVertexY, fMagneticField);
	  // _________________________________________________ //
	  // __ track information __ //
	  AliAnalysisTaskDoubleHypNucTree::SetDaughterInformation("4LLH");		
	  // _________________________________________________ //
	  // __ information of hypernuclei track __ //
	  Float_t xv[2],  yv[3];
	  exTrack->GetImpactParameters(xv,yv);
	  fpDaughter4        = exTrack->GetP();
	  fptDaughter4       = sublorentzsum->Pt();
	  fpxDaughter4       = sublorentzsum->Px();
	  fpyDaughter4       = sublorentzsum->Py();
	  fpzDaughter4       = sublorentzsum->Pz();
	  fEDaughter4        = sublorentzsum->E();
	  fyDaughter4        = sublorentzsum->Rapidity();	    
	  fDcaDaughter4      = xv[0];
	  fDcazDaughter4     = xv[1];
	  fSigmaYXDaughter4  = yv[0];
	  fSigmaXYZDaughter4 = yv[1];
	  fSigmaZDaughter4   = yv[2];
	  fPtUncertDaughter4 = TMath::Sqrt(exTrack->GetSigma1Pt2())*fptDaughter4;
	  fDcaSecDaughter4   = TMath::Abs(exTrack->GetD(SecVertex[0], SecVertex[1], fMagneticField));
	  // _________________________________________________ //
	  // __ mother hypernucleus information __ //
	  fPDGMother       = fgkPdgCode[kPDGDoubleHyperHydrogen4];
	  fChargeMother    = 1;
	  fmMother         = lorentzsum->M();
	  fmMother2        = lorentzsum2->M();
	  fEMother         = lorentzsum->E();
	  fpxMother        = lorentzsum->Px();
	  fpyMother        = lorentzsum->Py();
	  fpzMother        = lorentzsum->Pz();
	  fptMother        = lorentzsum->Pt();
	  fpMother         = lorentzsum->P();
	  fyMother         = lorentzsum->Rapidity();
	  hInvMass4LLH->Fill(fmMother);
	  // _________________________________________________ //
	  // __ daughter hypernucleus information __ //
	  fmSubMother      = sublorentzsum->M();
	  fESubMother      = sublorentzsum->E();
	  fpxSubMother     = sublorentzsum->Px();
	  fpySubMother     = sublorentzsum->Py();
	  fpzSubMother     = sublorentzsum->Pz();
	  fptSubMother     = sublorentzsum->Pt();
	  fpSubMother      = sublorentzsum->P();
	  fySubMother      = sublorentzsum->Rapidity();
	  // _________________________________________________ //
	  // __ mother hypernucleus ctau and cos(PA) __ //
	  dd[0] = PrimVertex[0] - SecVertex[0];
	  dd[1] = PrimVertex[1] - SecVertex[1];
	  dd[2] = PrimVertex[2] - SecVertex[2];
	  h->SetXYZ(-dd[0],-dd[1],-dd[2]);	 
	      
	  fctMother        = (fmMother*TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2)))/lorentzsum->P();
	  fPA              = TMath::Cos(lorentzsum->Angle(*h));	      	     
	  // _________________________________________________ //
	  // __ daughter hypernucleus ctau and cos(PA) to sec vtx __ //
	  dd[0] = SecVertex[0] - TertVertex[0];
	  dd[1] = SecVertex[1] - TertVertex[1];
	  dd[2] = SecVertex[2] - TertVertex[2];
	  h->SetXYZ(-dd[0],-dd[1],-dd[2]);
	      
	  fctSubMother     = (fmSubMother*TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2)))/sublorentzsum->P();
	  fSubPA           = TMath::Cos(sublorentzsum->Angle(*h));	  
	  // _________________________________________________ //
	  // __ daughter hypernucleus cos(PA) to prim vtx __ //
	  dd[0] = PrimVertex[0] - TertVertex[0];
	  dd[1] = PrimVertex[1] - TertVertex[1];
	  dd[2] = PrimVertex[2] - TertVertex[2];
	  h->SetXYZ(-dd[0],-dd[1],-dd[2]);

	  fSubPA2          = TMath::Cos(sublorentzsum->Angle(*h));
	  fDecAngle        = sublorentzsum->Angle(particle4->Vect());
	  if(fSubPA2 < 0.99) {
	    if(sublorentzsum)  delete sublorentzsum;
	    if(sublorentzsum2) delete sublorentzsum2;
	    if(lorentzsum)     delete lorentzsum;
	    if(lorentzsum2)    delete lorentzsum2;
	    if(particle1)      delete particle1;
	    if(particle2)      delete particle2;
	    if(particle3)      delete particle3;
	    if(particle4)      delete particle4;
	    if(exTrack)        delete exTrack;
	    if(h)              delete h;
	    AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	    continue;
	  }
	  // _________________________________________________ //
	  // __ daughter hypernucleus cos(PA) to sec vtx (2) __ //
	  dd[0] = SecVertex[0] - TertVertex2[0];
	  dd[1] = SecVertex[1] - TertVertex2[1];
	  dd[2] = SecVertex[2] - TertVertex2[2];
	  h->SetXYZ(-dd[0], -dd[1], -dd[2]);

	  fSubPA3 = TMath::Cos(sublorentzsum->Angle(*h));
	  // _________________________________________________ //
	  // __ daughter hypernucleus cos(PA) to prim vtx (2) __ //
	  dd[0] = PrimVertex[0] - TertVertex2[0];
	  dd[1] = PrimVertex[1] - TertVertex2[1];
	  dd[2] = PrimVertex[2] - TertVertex2[2];
	  h->SetXYZ(-dd[0], -dd[1], -dd[2]);

	  fSubPA4 = TMath::Cos(sublorentzsum->Angle(*h));
	  // _________________________________________________ //
	  // __ armenteros podolanski __ //
	  TVector3 vecN(0.,0.,0.);
	  TVector3 vecP(0.,0.,0.);
	  TVector3 vecM(0.,0.,0.); 

	  vecP      = sublorentzsum->Vect();
	  vecN      = particle4->Vect();
	  vecM      = lorentzsum->Vect();

	  fthetaP   = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
	  fthetaN   = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));
	  farmalpha = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
	  farmpt    = vecP.Mag()*sin(fthetaP);
	  // _________________________________________________ //
	  // __ MC part __ //
	  if(fMCtrue){
	    if(ParticleGrandMother1->PdgCode() == fgkPdgCode[kPDGDoubleHyperHydrogen4]
	       && ParticleGrandMother2->PdgCode() == fgkPdgCode[kPDGDoubleHyperHydrogen4]
	       && ParticleGrandMother3->PdgCode() == fgkPdgCode[kPDGDoubleHyperHydrogen4]
	       && ParticleMother4->PdgCode() == fgkPdgCode[kPDGDoubleHyperHydrogen4]
	       && labelGrandMother1 == labelGrandMother2
	       && labelGrandMother2 == labelGrandMother3
	       && labelGrandMother3 == labelMother4){
	      if(TMath::Abs(label1) == TMath::Abs(GetLabel(labelGrandMother1, fgkPdgCode[kPDGHelium3], fgkPdgCode[kPDGHyperHelium4]))
		 && TMath::Abs(label2) == TMath::Abs(GetLabel(labelGrandMother1, fgkPdgCode[kPDGProton], fgkPdgCode[kPDGHyperHelium4]))
		 && TMath::Abs(label3) == TMath::Abs(GetLabel(labelGrandMother1, fgkPdgCode[kPDGPionMinus], fgkPdgCode[kPDGHyperHelium4]))
		 && TMath::Abs(label4) == TMath::Abs(GetLabel(labelGrandMother1, fgkPdgCode[kPDGPionMinus]))){
		mctruth = 1;		      	     
	      }	  
	    }
	  }
	  // __ fill tree __ //
	  fTree->Fill();
	  // __ reset __ //
	  AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	  if(sublorentzsum)  delete sublorentzsum;
	  if(sublorentzsum2) delete sublorentzsum2;
	  if(lorentzsum)     delete lorentzsum;
	  if(lorentzsum2)    delete lorentzsum2;
	  if(particle1)      delete particle1;
	  if(particle2)      delete particle2;
	  if(particle3)      delete particle3;
	  if(particle4)      delete particle4;
	  if(exTrack)        delete exTrack;
	  if(h)              delete h;
	  // _________________________________________________ //
	}//track4
      }//track3
    }//track2
  }//track1
}
void AliAnalysisTaskDoubleHypNucTree::CalcNegMotherHypNuc(){
  // ____________________________ first track _______________________________ //
  for (const Int_t &He3NegTracks : He3NegArrayAcc) {

    if(He3NegTracks == 0) continue;

    track1 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(He3NegTracks));

    ptot1 = track1->GetInnerParam()->GetP();
    sign1 = track1->GetSign();

    // __ MC part __ //
    if(fMCtrue && noCombinatoricsMC){
	
      label1 = track1->GetLabel();
		  
      labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));
      labelGrandMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother1));

      AliMCParticle *part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(label1))->Particle());
      if(part->PdgCode() != fgkPdgCode[kPDGAntiHelium3])continue;
      if(part)delete part;
	      
      ParticleGrandMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother1))->Particle());

      if(ParticleGrandMother1->PdgCode() != fgkPdgCode[kPDGAntiDoubleHyperHydrogen4]) continue;
      if(ParticleGrandMother1) delete ParticleGrandMother1;
    }


    // ____________________________ second track _______________________________ //
    for (const Int_t &PNegTracks : PNegArrayAcc) {

      if(PNegTracks == 0) continue;
      // __ reject using same track twice __ //
      if(He3NegTracks == PNegTracks) continue;

      track2 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(PNegTracks));

      ptot2 = track2->GetInnerParam()->GetP();
      sign2 = track2->GetSign();

      // __ MC part __ //
      if(fMCtrue && noCombinatoricsMC){
	
	label2 = track2->GetLabel();
		  
	labelMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label2));
	labelGrandMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother2));

	AliMCParticle *part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(label2))->Particle());
	if(part->PdgCode() != fgkPdgCode[kPDGAntiProton]) continue;
	if(part) delete part;
	      
	ParticleGrandMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother2))->Particle());

	if(ParticleGrandMother2->PdgCode() != fgkPdgCode[kPDGAntiDoubleHyperHydrogen4]) continue;
	if(ParticleGrandMother2) delete ParticleGrandMother2;
      }
      // __ dca rejection betw tracks __ //
      if(TMath::Abs(track2->GetDCA(track1, fMagneticField, xthiss, xpp)) > 1.2) continue;

      // ____________________________ third track _______________________________ //
      for (const Int_t &PiPosTracks : PiPosArrayAcc) {

	if(PiPosTracks == 0) continue;

	// __ reject using same track twice __ //
	if (PNegTracks == PiPosTracks || He3NegTracks == PiPosTracks) continue;

	track3 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(PiPosTracks));        	
	
	ptot3 = track3->GetInnerParam()->GetP();
	sign3 = track3->GetSign();

	// __ MC part __ //
	if(fMCtrue && (lessCombinatoricsMC || noCombinatoricsMC)){
	
	  label3 = track3->GetLabel();
		  
	  labelMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label3));
	  labelGrandMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother3));

	  AliMCParticle *part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(label3))->Particle());
	  if(part->PdgCode() != fgkPdgCode[kPDGPionPlus]) continue;
	  if(part) delete part;
	      
	  ParticleGrandMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother3))->Particle());

	  if(ParticleGrandMother3->PdgCode() != fgkPdgCode[kPDGAntiDoubleHyperHydrogen4]) continue;
	  if(ParticleGrandMother3) delete ParticleGrandMother3;
	}

	// __ dca rejection betw tracks __ //
	if(TMath::Abs(track2->GetDCA(track1, fMagneticField, xthiss, xpp)) > 1.2
	   || TMath::Abs(track3->GetDCA(track1, fMagneticField, xthiss, xpp)) > 1.2
	   || TMath::Abs(track3->GetDCA(track2, fMagneticField, xthiss, xpp)) > 1.2) continue;

	// ____________________________ fourth track _______________________________ //
	for (const Int_t &PiPosSecTracks : PiPosSecArray) {

	  if(PiPosSecTracks == 0) continue;

	  // __ reject using same track twice __ //
	  if (PNegTracks == PiPosSecTracks || He3NegTracks == PiPosSecTracks || PiPosTracks == PiPosSecTracks) continue;

	  track4 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(PiPosSecTracks));	  	  

	  ptot4 = track4->GetInnerParam()->GetP();
	  sign4 = track4->GetSign();

	  // __ MC part __ //
	  if(fMCtrue && (lessCombinatoricsMC || noCombinatoricsMC)){
	
	    label4 = track4->GetLabel();
		  
	    labelMother4 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label4));
	      
	    ParticleMother4 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother4))->Particle());

	    AliMCParticle *part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(label4))->Particle());
	    if(part->PdgCode() != fgkPdgCode[kPDGPionPlus]) continue;
	    if(part)delete part;

	    if(ParticleMother4->PdgCode() != fgkPdgCode[kPDGAntiDoubleHyperHydrogen4]) continue;

	    if(ParticleMother4) delete ParticleMother4;	    
	  }

	  // _________________________________________________ //
	  lorentzsum     = new TLorentzVector(0.,0.,0.,0.);
	  particle1      = new TLorentzVector(0.,0.,0.,0.);
	  particle2      = new TLorentzVector(0.,0.,0.,0.);
	  particle3      = new TLorentzVector(0.,0.,0.,0.);
	  particle4      = new TLorentzVector(0.,0.,0.,0.);		
	  // _________________________________________________ //
	  particle1->SetXYZM(2.*track1->Px(), 2.*track1->Py(), 2.*track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
	  particle2->SetXYZM(   track2->Px(),    track2->Py(),    track2->Pz(), AliPID::ParticleMass(AliPID::kProton));
	  particle3->SetXYZM(   track3->Px(),    track3->Py(),    track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
	  particle4->SetXYZM(   track4->Px(),    track4->Py(),    track4->Pz(), AliPID::ParticleMass(AliPID::kPion));
	  lorentzsum->SetXYZM(0.,0.,0.,0.);
	  *lorentzsum     = *particle1 + *particle2 + *particle3 + *particle4;
	  // _________________________________________________ //
	  // __ select mother hypernucleus by inv mass __ //
	  if(lorentzsum->M() > 4.16 || lorentzsum->M() < 4.04) {
	    if(lorentzsum)    delete lorentzsum;
	    if(particle1)     delete particle1;
	    if(particle2)     delete particle2;
	    if(particle3)     delete particle3;
	    if(particle4)     delete particle4;
	    AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	    continue;
	  }
	  else{
	    if(lorentzsum)    delete lorentzsum;
	    if(particle1)     delete particle1;
	    if(particle2)     delete particle2;
	    if(particle3)     delete particle3;
	    if(particle4)     delete particle4;
	  }

	  // __ MC part __ //
	  if(fMCtrue){
	    label1 = 0;
	    label2 = 0;
	    label3 = 0;
	    label4 = 0;
		  
	    label1 = track1->GetLabel();
	    label2 = track2->GetLabel();
	    label3 = track3->GetLabel();	
	    label4 = track4->GetLabel();
		  
	    labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));
	    labelGrandMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother1));
	      
	    labelMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label2));
	    labelGrandMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother2));
	      
	    labelMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label3));
	    labelGrandMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother3));
	      
	    labelMother4 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label4));
	    labelGrandMother4 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother4));
	      
	    ParticleMother4 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother4))->Particle());

	    ParticleGrandMother4 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother4))->Particle());
	    ParticleGrandMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother3))->Particle());
	    ParticleGrandMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother2))->Particle());
	    ParticleGrandMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother1))->Particle());
	  }
	  // _________________________________________________ //
	  //********** #bar{4LLH} **********
	  // _________________________________________________ //
	  // __ trigger __ //
	  fTrigMB    = MB;		
	  fTrigHMV0  = HMV0;
	  fTrigHMSPD = HMSPD;
	  fTrigHNU   = HNU;
	  fTrigHQU   = HQU;
	  // _________________________________________________ //
	  sublorentzsum  = new TLorentzVector(0.,0.,0.,0.);
	  sublorentzsum2 = new TLorentzVector(0.,0.,0.,0.);
	  lorentzsum     = new TLorentzVector(0.,0.,0.,0.);
	  lorentzsum2    = new TLorentzVector(0.,0.,0.,0.);
	  particle1      = new TLorentzVector(0.,0.,0.,0.);
	  particle2      = new TLorentzVector(0.,0.,0.,0.);
	  particle3      = new TLorentzVector(0.,0.,0.,0.);
	  particle4      = new TLorentzVector(0.,0.,0.,0.);		
	  h              = new TVector3(0., 0., 0.);
	  // _________________________________________________ //
	  particle1->SetXYZM(2.*track1->Px(), 2.*track1->Py(), 2.*track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
	  particle2->SetXYZM(   track2->Px(),    track2->Py(),    track2->Pz(), AliPID::ParticleMass(AliPID::kProton));
	  particle3->SetXYZM(   track3->Px(),    track3->Py(),    track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
	  particle4->SetXYZM(   track4->Px(),    track4->Py(),    track4->Pz(), AliPID::ParticleMass(AliPID::kPion));
	  sublorentzsum->SetXYZM(0.,0.,0.,0.);
	  sublorentzsum2->SetXYZM(0.,0.,0.,0.);
	  lorentzsum->SetXYZM(0.,0.,0.,0.);
	  lorentzsum2->SetXYZM(0.,0.,0.,0.);
	  *lorentzsum     = *particle1 + *particle2 + *particle3 + *particle4;
	  *sublorentzsum  = *particle1 + *particle2 + *particle3;
	  sublorentzsum2->SetXYZM(sublorentzsum->X(), sublorentzsum->Y(), sublorentzsum->Z(), 3.929);	    
	  *lorentzsum2    = *sublorentzsum2 + *particle4;	    
	  // _________________________________________________ //
	  // __ coord of prim vtx __ //
	  fPrimVertexX = PrimVertex[0];
	  fPrimVertexY = PrimVertex[1];
	  fPrimVertexZ = PrimVertex[2];
	  // _________________________________________________ //
	  // __ create sec / tert vtx __ //
	  AliAnalysisTaskDoubleHypNucTree::CreateSecVertex();
	  // _________________________________________________ //
	  // __ coord of sec vtx __ //
	  fSecVertexX = SecVertex[0];
	  fSecVertexY = SecVertex[1];
	  fSecVertexZ = SecVertex[2];
	  // _________________________________________________ //
	  // __ coord of tert vtx __ //
	  // __ reco starts at prim vtx __ //
	  fTertVertexX = TertVertex[0];
	  fTertVertexY = TertVertex[1];
	  fTertVertexZ = TertVertex[2];
	  // __ reco starts at sec vtx __ //
	  fTertVertex2X = TertVertex2[0];
	  fTertVertex2Y = TertVertex2[1];
	  fTertVertex2Z = TertVertex2[2];
	  // _________________________________________________ //
	  // __ dca between tracks __ //
	  fDCA2B  = TMath::Abs(track4->GetDCA(exTrack, fMagneticField, xthiss, xpp));
	  fDCA3B1 = TMath::Abs(track2->GetDCA(track1,  fMagneticField, xthiss, xpp));
	  fDCA3B2 = TMath::Abs(track3->GetDCA(track1,  fMagneticField, xthiss, xpp));
	  fDCA3B3 = TMath::Abs(track3->GetDCA(track2,  fMagneticField, xthiss, xpp));	      
	  // __ cuts from 04/2020 __ //
	  fDCA2Bo        = TMath::Abs(track4->GetD(sublorentzsum->X(), sublorentzsum->Y(), fMagneticField));
	  fDcaDaughtero  = track1->GetD(fPrimVertexX, fPrimVertexY, fMagneticField);
	  fDcaDaughter1o = track2->GetD(fPrimVertexX, fPrimVertexY, fMagneticField);
	  fDcaDaughter2o = track3->GetD(fPrimVertexX, fPrimVertexY, fMagneticField);
	  fDcaDaughter3o = track4->GetD(fPrimVertexX, fPrimVertexY, fMagneticField);
	  // _________________________________________________ //
	  // __ track information __ //
	  AliAnalysisTaskDoubleHypNucTree::SetDaughterInformation("4LLH");		
	  // _________________________________________________ //
	  // __ information of hypernuclei track __ //
	  Float_t xv[2],  yv[3];
	  exTrack->GetImpactParameters(xv,yv);
	  fpDaughter4        = exTrack->GetP();
	  fptDaughter4       = sublorentzsum->Pt();
	  fpxDaughter4       = sublorentzsum->Px();
	  fpyDaughter4       = sublorentzsum->Py();
	  fpzDaughter4       = sublorentzsum->Pz();
	  fEDaughter4        = sublorentzsum->E();
	  fyDaughter4        = sublorentzsum->Rapidity();	    
	  fDcaDaughter4      = xv[0];
	  fDcazDaughter4     = xv[1];
	  fSigmaYXDaughter4  = yv[0];
	  fSigmaXYZDaughter4 = yv[1];
	  fSigmaZDaughter4   = yv[2];
	  fPtUncertDaughter4 = TMath::Sqrt(exTrack->GetSigma1Pt2())*fptDaughter4;
	  fDcaSecDaughter4   = TMath::Abs(exTrack->GetD(SecVertex[0], SecVertex[1], fMagneticField));
	  // _________________________________________________ //
	  // __ mother hypernucleus information __ //
	  fPDGMother       = fgkPdgCode[kPDGAntiDoubleHyperHydrogen4];
	  fChargeMother    = -1;
	  fmMother         = lorentzsum->M();
	  fmMother2        = lorentzsum2->M();
	  fEMother         = lorentzsum->E();
	  fpxMother        = lorentzsum->Px();
	  fpyMother        = lorentzsum->Py();
	  fpzMother        = lorentzsum->Pz();
	  fptMother        = lorentzsum->Pt();
	  fpMother         = lorentzsum->P();
	  fyMother         = lorentzsum->Rapidity();
	  hInvMass4LLH->Fill(fmMother);
	  // _________________________________________________ //
	  // __ daughter hypernucleus information __ //
	  fmSubMother      = sublorentzsum->M();
	  fESubMother      = sublorentzsum->E();
	  fpxSubMother     = sublorentzsum->Px();
	  fpySubMother     = sublorentzsum->Py();
	  fpzSubMother     = sublorentzsum->Pz();
	  fptSubMother     = sublorentzsum->Pt();
	  fpSubMother      = sublorentzsum->P();
	  fySubMother      = sublorentzsum->Rapidity();
	  // _________________________________________________ //
	  // __ mother hypernucleus ctau and cos(PA) __ //
	  dd[0] = PrimVertex[0] - SecVertex[0];
	  dd[1] = PrimVertex[1] - SecVertex[1];
	  dd[2] = PrimVertex[2] - SecVertex[2];
	  h->SetXYZ(-dd[0],-dd[1],-dd[2]);	 
	      
	  fctMother        = (fmMother*TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2)))/lorentzsum->P();
	  fPA              = TMath::Cos(lorentzsum->Angle(*h));	      	     
	  // _________________________________________________ //
	  // __ daughter hypernucleus ctau and cos(PA) to sec vtx __ //
	  dd[0] = SecVertex[0] - TertVertex[0];
	  dd[1] = SecVertex[1] - TertVertex[1];
	  dd[2] = SecVertex[2] - TertVertex[2];
	  h->SetXYZ(-dd[0],-dd[1],-dd[2]);
	      
	  fctSubMother     = (fmSubMother*TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2)))/sublorentzsum->P();
	  fSubPA           = TMath::Cos(sublorentzsum->Angle(*h));	  
	  // _________________________________________________ //
	  // __ daughter hypernucleus cos(PA) to prim vtx __ //
	  dd[0] = PrimVertex[0] - TertVertex[0];
	  dd[1] = PrimVertex[1] - TertVertex[1];
	  dd[2] = PrimVertex[2] - TertVertex[2];
	  h->SetXYZ(-dd[0],-dd[1],-dd[2]);

	  fSubPA2          = TMath::Cos(sublorentzsum->Angle(*h));
	  fDecAngle        = sublorentzsum->Angle(particle4->Vect());
	  if(fSubPA2 < 0.99) {
	    if(sublorentzsum)  delete sublorentzsum;
	    if(sublorentzsum2) delete sublorentzsum2;
	    if(lorentzsum)     delete lorentzsum;
	    if(lorentzsum2)    delete lorentzsum2;
	    if(particle1)      delete particle1;
	    if(particle2)      delete particle2;
	    if(particle3)      delete particle3;
	    if(particle4)      delete particle4;
	    if(exTrack)        delete exTrack;
	    if(h)              delete h;
	    AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	    continue;
	  }
	  // _________________________________________________ //
	  // __ daughter hypernucleus cos(PA) to sec vtx (2) __ //
	  dd[0] = SecVertex[0] - TertVertex2[0];
	  dd[1] = SecVertex[1] - TertVertex2[1];
	  dd[2] = SecVertex[2] - TertVertex2[2];
	  h->SetXYZ(-dd[0], -dd[1], -dd[2]);

	  fSubPA3 = TMath::Cos(sublorentzsum->Angle(*h));
	  // _________________________________________________ //
	  // __ daughter hypernucleus cos(PA) to prim vtx (2) __ //
	  dd[0] = PrimVertex[0] - TertVertex2[0];
	  dd[1] = PrimVertex[1] - TertVertex2[1];
	  dd[2] = PrimVertex[2] - TertVertex2[2];
	  h->SetXYZ(-dd[0], -dd[1], -dd[2]);

	  fSubPA4 = TMath::Cos(sublorentzsum->Angle(*h));
	  // _________________________________________________ //
	  // __ armenteros podolanski __ //
	  TVector3 vecN(0.,0.,0.);//
	  TVector3 vecP(0.,0.,0.);
	  TVector3 vecM(0.,0.,0.); 

	  vecN      = sublorentzsum->Vect();
	  vecP      = particle4->Vect();
	  vecM      = lorentzsum->Vect();

	  fthetaP   = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
	  fthetaN   = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));
	  farmalpha = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
	  farmpt    = vecP.Mag()*sin(fthetaP);
	  // _________________________________________________ //
	  // __ MC part __ //
	  if(fMCtrue){
	    if(ParticleGrandMother1->PdgCode() == fgkPdgCode[kPDGAntiDoubleHyperHydrogen4]
	       && ParticleGrandMother2->PdgCode() == fgkPdgCode[kPDGAntiDoubleHyperHydrogen4]
	       && ParticleGrandMother3->PdgCode() == fgkPdgCode[kPDGAntiDoubleHyperHydrogen4]
	       && ParticleMother4->PdgCode() == fgkPdgCode[kPDGAntiDoubleHyperHydrogen4]
	       && labelGrandMother1 == labelGrandMother2
	       && labelGrandMother2 == labelGrandMother3
	       && labelGrandMother3 == labelMother4){
	      if(TMath::Abs(label1) == TMath::Abs(GetLabel(labelGrandMother1, fgkPdgCode[kPDGAntiHelium3], fgkPdgCode[kPDGAntiHyperHelium4]))
		 && TMath::Abs(label2) == TMath::Abs(GetLabel(labelGrandMother1, fgkPdgCode[kPDGAntiProton], fgkPdgCode[kPDGAntiHyperHelium4]))
		 && TMath::Abs(label3) == TMath::Abs(GetLabel(labelGrandMother1, fgkPdgCode[kPDGPionPlus], fgkPdgCode[kPDGAntiHyperHelium4]))
		 && TMath::Abs(label4) == TMath::Abs(GetLabel(labelGrandMother1, fgkPdgCode[kPDGPionPlus]))){
		mctruth = 1;		      
	      }
	    }
	  }
	  // __ fill tree __ //
	  fTree->Fill();
	  // __ reset __ //
	  AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	  if(sublorentzsum)  delete sublorentzsum;
	  if(sublorentzsum2) delete sublorentzsum2;
	  if(lorentzsum)     delete lorentzsum;
	  if(lorentzsum2)    delete lorentzsum2;
	  if(particle1)      delete particle1;
	  if(particle2)      delete particle2;
	  if(particle3)      delete particle3;
	  if(particle4)      delete particle4;
	  if(exTrack)        delete exTrack;
	  if(h)              delete h;	  
	  // _________________________________________________ //
	}//track4
      }//track3
    }//track2
  }//track1  
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::CreateSecVertex(){
  // ******************************************************************************
  // * Two different methods:                                                     *
  // * fvariante = 2 first creates the tertiary vertex                            *
  // * which coordinates are taken to reconstruct the daughter hypernucleus track *
  // * MC checks showed, that the Pointing resolution is worse in this case       *
  // * fvariante = 1 uses the coordinates of the primary vertex to reconstruct    *
  // * the daughter hypernucleus track                                            *
  // ******************************************************************************
  Double_t impar[2];
  // __ daughter hypernucleus track __ //
  exTrack = new AliExternalTrackParam();  
  // __ init AliVertexerTracks __ //
  secvertexer = new AliVertexerTracks(fESDevent->GetMagneticField());
  tertvertexer = new AliVertexerTracks(fESDevent->GetMagneticField());
  //secvertexer->SetITSMode();
  //tertvertexer->SetITSMode();
  // __ get prim vtx __ //
  primVertex = new AliESDVertex(*vertex);
  // _________________________________ //
  if(fvariante == 2) {
    // __ tert vtx __ //
    trkArray = new TObjArray(3);
    trkArray->AddAt(track1,0);
    trkArray->AddAt(track2,1);
    trkArray->AddAt(track3,2);
    // __start at prim vtx __ //
    tertvertexer->SetVtxStart(primVertex);
    tertVertex = (AliESDVertex*)tertvertexer->VertexForSelectedESDTracks(trkArray);
    if(trkArray) delete trkArray;
    TertVertex[0] = tertVertex->GetX();
    TertVertex[1] = tertVertex->GetY();
    TertVertex[2] = tertVertex->GetZ();
    // _________________________________ //
    // __ propagate tracks and get impact parameters __ //
    if(track1->PropagateToDCA(tertVertex, fMagneticField, 10, impar)) fPropDCADaughter = 1;
    fImParDaughter = impar[0];
    fImParzDaughter = impar[1];
    if(track2->PropagateToDCA(tertVertex, fMagneticField, 10, impar)) fPropDCADaughter1 = 1;
    fImParDaughter1 = impar[0];
    fImParzDaughter1 = impar[1];
    if(track3->PropagateToDCA(tertVertex, fMagneticField, 10, impar)) fPropDCADaughter2 = 1;
    fImParDaughter2 = impar[0];
    fImParzDaughter2 = impar[1];
    if(tertVertex) delete tertVertex;		
    // _________________________________ //
    // __ Setting up daughter hypernucleus track here __ //
    // __ cov matrix calculated from daughter tracks __ //
    track1->CheckCovariance(); track2->CheckCovariance(); track3->CheckCovariance();
    track1->GetCovarianceXYZPxPyPz(cov0); track2->GetCovarianceXYZPxPyPz(cov1); track3->GetCovarianceXYZPxPyPz(cov2);
    for(int i = 0; i<21; i++) cov[i] = TMath::Sqrt(cov0[i]*cov0[i] + cov1[i]*cov1[i] + cov2[i]*cov2[i]);
    // __ daughter hypernucleus sign = helium track sign __ //
    sign = sign1;
    // __ use coord from tert vtx for daughter hypernucleus track __ //
    xyz[0]=TertVertex[0]; xyz[1]=TertVertex[1]; xyz[2]=TertVertex[2];
    // __ daughter hypernucleus momentum __ //
    pxpypz[0] = sublorentzsum->Px(); pxpypz[1] = sublorentzsum->Py(); pxpypz[2] = sublorentzsum->Pz();
    exTrack->Set(xyz,pxpypz,cov,sign);
    // __ check cov __ //
    exTrack->CheckCovariance();
    track4->CheckCovariance();
    // _________________________________ //
    // __ create sec vtx __ //
    trkArray1 = new TObjArray(2);
    trkArray1->AddAt(exTrack,0);
    trkArray1->AddAt(track4,1);
    id[0]=0;id[1]=1;
    // __ start at prim vtx __ //
    secvertexer->SetVtxStart(primVertex);
    secVertex = (AliESDVertex*)secvertexer->VertexForSelectedTracks(trkArray1, id, kTRUE, kTRUE, kFALSE);
    if(trkArray1) delete trkArray1;
    SecVertex[0] = secVertex->GetX();
    SecVertex[1] = secVertex->GetY();
    SecVertex[2] = secVertex->GetZ();	      	    		  
    // _________________________________ //
    // __ propagate tracks and get impact parameters __ //
    if(track4->PropagateToDCA(secVertex, fMagneticField, 10, impar)) fPropDCADaughter3 = 1;
    fImParDaughter3 = impar[0];
    fImParzDaughter3 = impar[1];
    if(exTrack->PropagateToDCA(secVertex, fMagneticField, 10, impar)) fPropDCADaughter4 = 1;
    fImParDaughter4 = impar[0];
    fImParzDaughter4 = impar[1];
    if(secVertex) delete secVertex;
  }
  // _________________________________ //
  else {
    // __ Setting up daughter hypernucleus track here __ //
    // __ cov matrix calculated from daughter tracks __ //
    track1->CheckCovariance(); track2->CheckCovariance(); track3->CheckCovariance();
    track1->GetCovarianceXYZPxPyPz(cov0); track2->GetCovarianceXYZPxPyPz(cov1); track3->GetCovarianceXYZPxPyPz(cov2);
    for(int i = 0; i<21; i++) cov[i] = TMath::Sqrt(cov0[i]*cov0[i] + cov1[i]*cov1[i] + cov2[i]*cov2[i]);
    // __ daughter hypernucleus sign = helium track sign __ //
    sign = sign1;
    // __ use coord from prim vtx for daughter hypernucleus track __ //
    xyz[0]=PrimVertex[0]; xyz[1]=PrimVertex[1]; xyz[2]=PrimVertex[2];
    // __ daughter hypernucleus momentum __ //
    pxpypz[0] = sublorentzsum->Px(); pxpypz[1] = sublorentzsum->Py(); pxpypz[2] = sublorentzsum->Pz();	      
    exTrack->Set(xyz,pxpypz,cov,sign);
    // __ check cov __ //
    exTrack->CheckCovariance();
    track4->CheckCovariance();
    // _________________________________ //
    // __ create sec vtx __ //
    trkArray1 = new TObjArray(2);
    trkArray1->AddAt(exTrack,0);
    trkArray1->AddAt(track4,1);
    id[0]=0;id[1]=1;
    // __ start at prim vtx __ //
    secvertexer->SetVtxStart(primVertex);
    secVertex = (AliESDVertex*)secvertexer->VertexForSelectedTracks(trkArray1, id,  kTRUE, kTRUE, kFALSE);
    if(trkArray1) delete trkArray1;
    SecVertex[0] = secVertex->GetX();
    SecVertex[1] = secVertex->GetY();
    SecVertex[2] = secVertex->GetZ();
    // _________________________________ //
    // __ propagate tracks and get impact parameters __ //
    if(track4->PropagateToDCA(secVertex, fMagneticField, 10, impar)) fPropDCADaughter3 = 1;
    fImParDaughter3  = impar[0];
    fImParzDaughter3 = impar[1];
    if(exTrack->PropagateToDCA(secVertex, fMagneticField, 10, impar)) fPropDCADaughter4 = 1;
    fImParDaughter4  = impar[0];
    fImParzDaughter4 = impar[1];
    // _________________________________ //
    // __ tert vtx __ //
    trkArray = new TObjArray(3);
    trkArray->AddAt(track1,0);
    trkArray->AddAt(track2,1);
    trkArray->AddAt(track3,2);
    // __ start at prim vtx __ //
    tertvertexer->SetVtxStart(primVertex);
    tertVertex = (AliESDVertex*)tertvertexer->VertexForSelectedESDTracks(trkArray);
    if(trkArray) delete trkArray;
    TertVertex[0] = tertVertex->GetX();
    TertVertex[1] = tertVertex->GetY();
    TertVertex[2] = tertVertex->GetZ();
    // _________________________________ //
    // __ propagate tracks and get impact parameters __ //
    if(track1->PropagateToDCA(tertVertex, fMagneticField, 10, impar)) fPropDCADaughter = 1;
    fImParDaughter  = impar[0];
    fImParzDaughter = impar[1];
    if(track2->PropagateToDCA(tertVertex, fMagneticField, 10, impar)) fPropDCADaughter1 = 1;
    fImParDaughter1  = impar[0];
    fImParzDaughter1 = impar[1];
    if(track3->PropagateToDCA(tertVertex, fMagneticField, 10, impar)) fPropDCADaughter2 = 1;
    fImParDaughter2  = impar[0];
    fImParzDaughter2 = impar[1];
    if(tertVertex) delete tertVertex;    
    // _________________________________ //
    //__ do it again for vertex starting at sec vertex __ //
    tertvertexer->SetVtxStart(secVertex);
    if (secVertex) delete secVertex;
    // __ tert vtx __ //
    trkArray = new TObjArray(3);
    trkArray->AddAt(track1, 0);
    trkArray->AddAt(track2, 1);
    trkArray->AddAt(track3, 2);
    // __ start at sec vtx __ //
    tertVertex = (AliESDVertex*)tertvertexer->VertexForSelectedESDTracks(trkArray);
    if (trkArray) delete trkArray;
    TertVertex2[0] = tertVertex->GetX();
    TertVertex2[1] = tertVertex->GetY();
    TertVertex2[2] = tertVertex->GetZ();
    // _________________________________ //
    // __ propagate tracks and get impact parameters __ //
    track1->PropagateToDCA(tertVertex, fMagneticField, 10, impar);
    fImPar2Daughter  = impar[0];
    fImParz2Daughter = impar[1];
    track2->PropagateToDCA(tertVertex, fMagneticField, 10, impar);
    fImPar2Daughter1  = impar[0];
    fImParz2Daughter1 = impar[1];
    track3->PropagateToDCA(tertVertex, fMagneticField, 10, impar);
    fImPar2Daughter2  = impar[0];
    fImParz2Daughter2 = impar[1];
    if (tertVertex) delete tertVertex;
  }
  if(secvertexer)  delete secvertexer;
  if(tertvertexer) delete tertvertexer;
  if(primVertex)   delete primVertex;
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::SetDaughterInformation(TString Mother) {

  Float_t xv[2],  yv[3];
  if(Mother == "4LHe" || Mother == "4LLH"){
    // _________________________________ //
    // __ He3 track __ //
    fpDaughter          = track1->GetInnerParam()->GetP();
    fpxDaughter         = particle1->Px();
    fpyDaughter         = particle1->Py();
    fpzDaughter         = particle1->Pz();
    fptDaughter         = particle1->Pt();
    fEDaughter          = particle1->E();
    fyDaughter          = particle1->Rapidity();
    fdEdxDaughter       = track1->GetTPCsignal();
    fdEdxSigmaDaughter  = AliAnalysisTaskDoubleHypNucTree::Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
    track1->GetImpactParameters(xv,yv);
    fDcaDaughter        = xv[0];
    fDcazDaughter       = xv[1];
    fSigmaYXDaughter    = yv[0];
    fSigmaXYZDaughter   = yv[1];
    fSigmaZDaughter     = yv[2];
    fDcaSecDaughter     = TMath::Abs(track1->GetD(TertVertex[0], TertVertex[1], fMagneticField));
    fNclsDaughter       = track1->GetTPCNcls();
    fNclsITSDaughter    = track1->GetNumberOfITSClusters();
    fChi2Daughter       = track1->GetTPCchi2() / (Float_t) track1->GetTPCclusters(0);
    fEtaDaughter        = track1->Eta();
    fPhiDaughter        = track1->Phi();
    fGeoLengthDaughter  = AliAnalysisTaskDoubleHypNucTree::GeoLength(*track1);
    fTOFSignalDaughter  = AliAnalysisTaskDoubleHypNucTree::GetTOFSignal(*track1);
    fPtUncertDaughter   = TMath::Sqrt(track1->GetSigma1Pt2())*fptDaughter;
    fTPCRefitDaughter   = (track1->GetStatus() & AliESDtrack::kTPCrefit) != 0;
    fITSRefitDaughter   = (track1->GetStatus() & AliESDtrack::kITSrefit) != 0;  
    fdEdxSigmaTriton    = AliAnalysisTaskDoubleHypNucTree::Bethe(*track1, AliPID::ParticleMass(AliPID::kTriton), 1, fBetheParamsT);
    fdEdxSigmaAlpha     = AliAnalysisTaskDoubleHypNucTree::Bethe(*track1, AliPID::ParticleMass(AliPID::kAlpha),  2, fBetheParamsHe);
    // _________________________________ //
    // __ p track __ //
    fpDaughter1         = track2->GetInnerParam()->GetP();
    fpxDaughter1        = particle2->Px();
    fpyDaughter1        = particle2->Py();
    fpzDaughter1        = particle2->Pz();
    fptDaughter1        = particle2->Pt();
    fEDaughter1         = particle2->E();
    fyDaughter1         = particle2->Rapidity();
    fdEdxDaughter1      = track2->GetTPCsignal();
    fdEdxSigmaDaughter1 =  AliAnalysisTaskDoubleHypNucTree::Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
    track2->GetImpactParameters(xv,yv);
    fDcaDaughter1       = xv[0];
    fDcazDaughter1      = xv[1];
    fSigmaYXDaughter1   = yv[0];
    fSigmaXYZDaughter1  = yv[1];
    fSigmaZDaughter1    = yv[2];
    fDcaSecDaughter1    = TMath::Abs(track2->GetD(TertVertex[0], TertVertex[1], fMagneticField));
    fNclsDaughter1      = track2->GetTPCNcls();
    fNclsITSDaughter1   = track2->GetNumberOfITSClusters();
    fChi2Daughter1      = track2->GetTPCchi2() / (Float_t) track2->GetTPCclusters(0);
    fEtaDaughter1       = track2->Eta();
    fPhiDaughter1       = track2->Phi();
    fGeoLengthDaughter1 = AliAnalysisTaskDoubleHypNucTree::GeoLength(*track2);
    fTOFSignalDaughter1 = AliAnalysisTaskDoubleHypNucTree::GetTOFSignal(*track2);
    fPtUncertDaughter1  = TMath::Sqrt(track2->GetSigma1Pt2())*fptDaughter1;
    fTPCRefitDaughter1  = (track2->GetStatus() & AliESDtrack::kTPCrefit) != 0;
    fITSRefitDaughter1  = (track2->GetStatus() & AliESDtrack::kITSrefit) != 0;
    fdEdxSigmaPion      =  fPID->NumberOfSigmasTPC(track2, AliPID::kPion);
    fdEdxSigmaDeuteron  =  AliAnalysisTaskDoubleHypNucTree::Bethe(*track2, AliPID::ParticleMass(AliPID::kDeuteron), 1, fBetheParamsT);
    // _________________________________ //
    // __ pion track __ //
    fpDaughter2         = track3->GetInnerParam()->GetP();
    fpxDaughter2        = particle3->Px();
    fpyDaughter2        = particle3->Py();
    fpzDaughter2        = particle3->Pz();
    fptDaughter2        = particle3->Pt();
    fEDaughter2         = particle3->E();
    fyDaughter2         = particle3->Rapidity();
    fdEdxDaughter2      = track3->GetTPCsignal();
    fdEdxSigmaDaughter2 = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
    track3->GetImpactParameters(xv,yv);
    fDcaDaughter2       = xv[0];
    fDcazDaughter2      = xv[1];
    fSigmaYXDaughter2   = yv[0];
    fSigmaXYZDaughter2  = yv[1];
    fSigmaZDaughter2    = yv[2];
    fDcaSecDaughter2    = TMath::Abs(track3->GetD(TertVertex[0], TertVertex[1], fMagneticField));
    fNclsDaughter2      = track3->GetTPCNcls();
    fNclsITSDaughter2   = track3->GetNumberOfITSClusters();
    fChi2Daughter2      = track3->GetTPCchi2() / (Float_t) track3->GetTPCclusters(0);
    fEtaDaughter2       = track3->Eta();
    fPhiDaughter2       = track3->Phi();
    fGeoLengthDaughter2 = AliAnalysisTaskDoubleHypNucTree::GeoLength(*track3);
    fTOFSignalDaughter2 = AliAnalysisTaskDoubleHypNucTree::GetTOFSignal(*track3);
    fPtUncertDaughter2  = TMath::Sqrt(track3->GetSigma1Pt2())*fptDaughter2;
    fTPCRefitDaughter2  = (track3->GetStatus() & AliESDtrack::kTPCrefit) != 0;
    fITSRefitDaughter2  = (track3->GetStatus() & AliESDtrack::kITSrefit) != 0;
  }
  if(Mother == "4LLH"){
    // _________________________________ //
    // __ pion track __ //
    fpDaughter3         = track4->GetInnerParam()->GetP();
    fpxDaughter3        = particle4->Px();
    fpyDaughter3        = particle4->Py();
    fpzDaughter3        = particle4->Pz();
    fptDaughter3        = particle4->Pt();
    fEDaughter3         = particle4->E();
    fyDaughter3         = particle4->Rapidity();
    fdEdxDaughter3      = track4->GetTPCsignal();
    fdEdxSigmaDaughter3 = fPID->NumberOfSigmasTPC(track4, AliPID::kPion);
    track4->GetImpactParameters(xv,yv);
    fDcaDaughter3       = xv[0];
    fDcazDaughter3      = xv[1];
    fSigmaYXDaughter3   = yv[0];
    fSigmaXYZDaughter3  = yv[1];
    fSigmaZDaughter3    = yv[2];
    fDcaSecDaughter3    = TMath::Abs(track4->GetD(SecVertex[0], SecVertex[1], fMagneticField));
    fNclsDaughter3      = track4->GetTPCNcls();
    fNclsITSDaughter3   = track4->GetNumberOfITSClusters();
    fChi2Daughter3      = track4->GetTPCchi2() / (Float_t) track4->GetTPCclusters(0);
    fEtaDaughter3       = track4->Eta();
    fPhiDaughter3       = track4->Phi();
    fGeoLengthDaughter3 = AliAnalysisTaskDoubleHypNucTree::GeoLength(*track4);
    fTOFSignalDaughter3 = AliAnalysisTaskDoubleHypNucTree::GetTOFSignal(*track4);
    fPtUncertDaughter3  = TMath::Sqrt(track4->GetSigma1Pt2())*fptDaughter3;
    fTPCRefitDaughter3  = (track4->GetStatus() & AliESDtrack::kTPCrefit) != 0;
    fITSRefitDaughter3  = (track4->GetStatus() & AliESDtrack::kITSrefit) != 0;
  }
  return;
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::Terminate(const Option_t*) {
  if (!GetOutputData(0)) return;
}
/// Set trigger information in reduced event
/// \return returns kTRUE is successful.
// _________________________________________________ //
Bool_t AliAnalysisTaskDoubleHypNucTree::TriggerSelection() {
  //******************************
  //*   get trigger information  *
  //******************************

  MB    = 0;
  HMV0  = 0;
  HMSPD = 0;
  HNU   = 0;
  HQU   = 0;

  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7)        MB    = kTRUE;
  if (fInputHandler->IsEventSelected() & AliVEvent::kHighMultV0)  HMV0  = kTRUE;
  if (fInputHandler->IsEventSelected() & AliVEvent::kHighMultSPD) HMSPD = kTRUE;
	
  Int_t nTrdTracks = fESDevent->GetNumberOfTrdTracks();
  if (!fMCtrue){
    // __ Data: get TRD trigger information from trigger classes __ //
    TString classes = fESDevent->GetFiredTriggerClasses();   
    if (classes.Contains("HNU")) HNU = 1;
    if (classes.Contains("HQU")) HQU = 1; 
		
  } else {
    // __ MC: simulate TRD trigger __ //
    if (nTrdTracks > 0) {
      for (Int_t iTrack = 0; iTrack < nTrdTracks; ++iTrack) {
	AliESDTrdTrack* trdTrack = fESDevent->GetTrdTrack(iTrack);
	if (!trdTrack) continue;					
	// __ simulate HNU __ //
	if((trdTrack->GetPID() >= 255 && trdTrack->GetNTracklets() == 4) || 
	   (trdTrack->GetPID() >= 235 && trdTrack->GetNTracklets() > 4)) {	
	  HNU = 1;
	}
	// __ simulate HQU __ //
	if (TMath::Abs(trdTrack->GetPt()) >= 256 &&
	    trdTrack->GetPID() >= 130 && trdTrack->GetNTracklets() >= 5 && (trdTrack->GetLayerMask() & 1) ){	
	  Float_t sag = AliAnalysisTaskDoubleHypNucTree::GetInvPtDevFromBC(trdTrack->GetB(), trdTrack->GetC());
	  if (sag < 0.2 && sag > -0.2) {
	    HQU = 1;
	  }
	}
      }     
    }
  } 
  // fill histogram
  fHistTrigger->Fill(0);
  if (MB)    fHistTrigger->Fill(1);
  if (HMV0)  fHistTrigger->Fill(2);
  if (HMSPD) fHistTrigger->Fill(3);
  if (HNU)   fHistTrigger->Fill(4);
  if (HQU)   fHistTrigger->Fill(5);
  Bool_t isTriggered = kFALSE;
  if(MB || HMV0 || HMSPD || HNU || HQU) isTriggered = kTRUE;
  return isTriggered;
}
// _________________________________________________ //
// __ calculate stack label of specific daughter __ //
Int_t AliAnalysisTaskDoubleHypNucTree::GetLabel(Int_t labelGrandMother, Int_t particlePdgCode){

  Int_t labelFirstDaughter = mcEvent->GetLabelOfParticleFirstDaughter(TMath::Abs(labelGrandMother));
  Int_t labelLastDaughter  = mcEvent->GetLabelOfParticleLastDaughter(TMath::Abs(labelGrandMother));
  Int_t diff               = TMath::Abs(labelLastDaughter - labelFirstDaughter) + 1;  
  Int_t returnval          = -99;
  
  for(Int_t i = 0; i<diff; i++){
      
    Int_t labelDaughter     = labelFirstDaughter + i;
    AliMCParticle *Daughter = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelDaughter))->Particle());
      
    if(Daughter->PdgCode() == particlePdgCode){
      returnval = labelDaughter;	 
    }
    if(Daughter) delete Daughter;
  }
    
  return returnval;
}
// _________________________________________________ //
// __ calculate stack label of specific granddaughter __ //
Int_t AliAnalysisTaskDoubleHypNucTree::GetLabel(Int_t labelGrandMother, Int_t particlePdgCode, Int_t motherparticlePdgCode){
  Int_t labelFirstDaughter = mcEvent->GetLabelOfParticleFirstDaughter(TMath::Abs(labelGrandMother));
  Int_t labelLastDaughter  = mcEvent->GetLabelOfParticleLastDaughter(TMath::Abs(labelGrandMother));
  Int_t diff               = TMath::Abs(labelLastDaughter - labelFirstDaughter) + 1;
  Int_t returnval          = -99;
  
  for(Int_t i = 0; i<diff;i++){
    
    Int_t labelDaughter     = labelFirstDaughter + i;
    AliMCParticle *Daughter = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelDaughter))->Particle());
    
    if(Daughter->PdgCode() == motherparticlePdgCode){
      
      Int_t labelFirstEnkel =  mcEvent->GetLabelOfParticleFirstDaughter(TMath::Abs(labelDaughter));
      Int_t labelLastEnkel  = mcEvent->GetLabelOfParticleLastDaughter(TMath::Abs(labelDaughter));
      Int_t diffEnkel        = TMath::Abs(labelLastEnkel - labelFirstEnkel) + 1;
      
      for(Int_t j = 0; j<diffEnkel; j++){
	
	Int_t labelEnkel     = labelFirstEnkel + j;
	AliMCParticle *Enkel =  new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelEnkel))->Particle());
	
	if(Enkel->PdgCode() == particlePdgCode){
	  returnval = labelEnkel;			     
	}
	if(Enkel) delete Enkel;
      }
    }		
    if(Daughter) delete Daughter;
  }
  return returnval;
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::MCGenerated() {

  // Monte Carlo for genenerated particles                                                                                                                                          
  stackN = 0;
  for(stackN = 0; stackN < mcEvent->GetNumberOfTracks(); stackN++) //loop over stack                                                                                                      
    {		 	

      ParticleMother = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(stackN))->Particle());
      PDGCodeMother  = ParticleMother->PdgCode();        		
	
      //DoubleHyperHydrogen4
      if(PDGCodeMother == fgkPdgCode[kPDGDoubleHyperHydrogen4]) //check mother PDG                                                                                                        
	{			 
	  AliAnalysisTaskDoubleHypNucTree::MCFourBodyDecay(stackN, ParticleMother, PDGCodeMother, fgkPdgCode[kPDGHelium3],
							   fgkPdgCode[kPDGProton], fgkPdgCode[kPDGPionMinus], fgkPdgCode[kPDGPionMinus],
							   AliPID::ParticleMass(AliPID::kHe3), AliPID::ParticleMass(AliPID::kProton), AliPID::ParticleMass(AliPID::kPion), AliPID::ParticleMass(AliPID::kPion));
	}
      //AntiDoubleHyperHydrogen4
      if(PDGCodeMother == fgkPdgCode[kPDGAntiDoubleHyperHydrogen4]) //check mother PDG
	{		 
		   	  
	  AliAnalysisTaskDoubleHypNucTree::MCFourBodyDecay(stackN, ParticleMother, PDGCodeMother, fgkPdgCode[kPDGAntiHelium3],
							   fgkPdgCode[kPDGAntiProton], fgkPdgCode[kPDGPionPlus], fgkPdgCode[kPDGPionPlus],
							   AliPID::ParticleMass(AliPID::kHe3), AliPID::ParticleMass(AliPID::kProton), AliPID::ParticleMass(AliPID::kPion), AliPID::ParticleMass(AliPID::kPion));
	}
      if(PDGCodeMother == fgkPdgCode[kPDGHyperHelium4]) //check mother PDG                                                                                                        
	{			 
	  AliAnalysisTaskDoubleHypNucTree::MCThreeBodyDecay(stackN, ParticleMother, PDGCodeMother, fgkPdgCode[kPDGHelium3],
							    fgkPdgCode[kPDGProton], fgkPdgCode[kPDGPionMinus],
							    AliPID::ParticleMass(AliPID::kHe3), AliPID::ParticleMass(AliPID::kProton), AliPID::ParticleMass(AliPID::kPion));
	}
      //AntiDoubleHyperHydrogen4
      if(PDGCodeMother == fgkPdgCode[kPDGAntiHyperHelium4]) //check mother PDG
	{		 
		   	  
	  AliAnalysisTaskDoubleHypNucTree::MCThreeBodyDecay(stackN, ParticleMother, PDGCodeMother, fgkPdgCode[kPDGAntiHelium3],
							    fgkPdgCode[kPDGAntiProton], fgkPdgCode[kPDGPionPlus],
							    AliPID::ParticleMass(AliPID::kHe3), AliPID::ParticleMass(AliPID::kProton), AliPID::ParticleMass(AliPID::kPion));
	}
      if(ParticleMother) delete ParticleMother;
    }//end loop over stack                                                                                                                                                                              
}
//_________________________________________________ 
void AliAnalysisTaskDoubleHypNucTree::MCFourBodyDecay(Int_t stackN, AliMCParticle *ParticleMother, Long_t PDGMother, Long_t PDGFirstDaughter, Long_t PDGSecondDaughter, Long_t PDGThirdDaughter, Long_t PDGFourthDaughter, Double_t massFirstDaughter, Double_t massSecondDaughter, Double_t massThirdDaughter, Double_t massFourthDaughter){

  fTrigMB    = MB;		
  fTrigHMV0  = HMV0;
  fTrigHMSPD = HMSPD;
  fTrigHNU   = HNU;
  fTrigHQU   = HQU;

  sign = PDGMother/TMath::Abs(PDGMother);

  FirstDaughter  = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(AliAnalysisTaskDoubleHypNucTree::GetLabel(stackN, PDGFirstDaughter,  sign*fgkPdgCode[kPDGHyperHelium4])))->Particle());
  SecondDaughter = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(AliAnalysisTaskDoubleHypNucTree::GetLabel(stackN, PDGSecondDaughter, sign*fgkPdgCode[kPDGHyperHelium4])))->Particle());
  ThirdDaughter  = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(AliAnalysisTaskDoubleHypNucTree::GetLabel(stackN, PDGThirdDaughter,  sign*fgkPdgCode[kPDGHyperHelium4])))->Particle());
  FourthDaughter = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(AliAnalysisTaskDoubleHypNucTree::GetLabel(stackN, PDGFourthDaughter)))->Particle());

  sublorentzsum = new TLorentzVector(0.,0.,0.,0.);
  lorentzsum    = new TLorentzVector(0.,0.,0.,0.);
  particle1     = new TLorentzVector(0.,0.,0.,0.);
  particle2     = new TLorentzVector(0.,0.,0.,0.);
  particle3     = new TLorentzVector(0.,0.,0.,0.);
  particle4     = new TLorentzVector(0.,0.,0.,0.);		
	      
  particle1->SetXYZM(FirstDaughter->Px(),   FirstDaughter->Py(),  FirstDaughter->Pz(),  massFirstDaughter);
  particle2->SetXYZM(SecondDaughter->Px(),  SecondDaughter->Py(), SecondDaughter->Pz(), massSecondDaughter);
  particle3->SetXYZM(ThirdDaughter->Px(),   ThirdDaughter->Py(),  ThirdDaughter->Pz(),  massThirdDaughter);
  particle4->SetXYZM(FourthDaughter->Px(),  FourthDaughter->Py(), FourthDaughter->Pz(), massFourthDaughter);
  sublorentzsum->SetXYZM(0.,0.,0.,0.);
  lorentzsum->SetXYZM(0.,0.,0.,0.);
  *lorentzsum    = *particle1 + *particle2 + *particle3 + *particle4;
  *sublorentzsum = *particle1 + *particle2 + *particle3;
	  
  if(FirstDaughter->PdgCode() == PDGFirstDaughter){
    if(SecondDaughter->PdgCode() == PDGSecondDaughter){
      if(ThirdDaughter->PdgCode() == PDGThirdDaughter){
	if(FourthDaughter->PdgCode() == PDGFourthDaughter){

	  dd[0] = ParticleMother->Xv() - FourthDaughter->Xv();
	  dd[1] = ParticleMother->Yv() - FourthDaughter->Yv();
	  dd[2] = ParticleMother->Zv() - FourthDaughter->Zv();
	  
	  fPDGMother    = (Int_t)PDGMother;
	  fChargeMother = sign;
	  fmMother      = lorentzsum->M();
	  fptMother     = lorentzsum->Pt();
	  fyMother      = lorentzsum->Rapidity();
	  fctMother     = (lorentzsum->M()*TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2)))/lorentzsum->P();

	  dd[0] = FourthDaughter->Xv() - FirstDaughter->Xv();
	  dd[1] = FourthDaughter->Yv() - FirstDaughter->Yv();
	  dd[2] = FourthDaughter->Zv() - FirstDaughter->Zv();
	  
	  fmSubMother   = sublorentzsum->M();
	  fptSubMother  = sublorentzsum->Pt();
	  fySubMother   = sublorentzsum->Rapidity();
	  fctSubMother  = (sublorentzsum->M()*TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2)))/sublorentzsum->P();
	      
	  fptDaughter   = particle1->Pt();
	  fyDaughter    = particle1->Rapidity();
	  fptDaughter1  = particle2->Pt();
	  fyDaughter1   = particle2->Rapidity();
	  fptDaughter2  = particle3->Pt();
	  fyDaughter2   = particle3->Rapidity();
	  fptDaughter3  = particle4->Pt();
	  fyDaughter3   = particle4->Rapidity();
	  fTreeGen->Fill();
	  AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	}
      }
    }
  }
  if(FirstDaughter)  delete FirstDaughter;
  if(SecondDaughter) delete SecondDaughter;
  if(ThirdDaughter)  delete ThirdDaughter;
  if(FourthDaughter) delete FourthDaughter;
  if(sublorentzsum) delete sublorentzsum;
  if(lorentzsum)    delete lorentzsum;
  if(particle1)     delete particle1;
  if(particle2)     delete particle2;
  if(particle3)     delete particle3;
  if(particle4)     delete particle4;
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::MCThreeBodyDecay(Int_t stackN, AliMCParticle *ParticleMother, Long_t PDGMother, Long_t PDGFirstDaughter, Long_t PDGSecondDaughter, Long_t PDGThirdDaughter, Double_t massFirstDaughter, Double_t massSecondDaughter, Double_t massThirdDaughter){

  fTrigMB    = MB;		
  fTrigHMV0  = HMV0;
  fTrigHMSPD = HMSPD;
  fTrigHNU   = HNU;
  fTrigHQU   = HQU;

  sign = PDGMother/TMath::Abs(PDGMother);

  FirstDaughter  = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(AliAnalysisTaskDoubleHypNucTree::GetLabel(stackN, PDGFirstDaughter)))->Particle());
  SecondDaughter = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(AliAnalysisTaskDoubleHypNucTree::GetLabel(stackN, PDGSecondDaughter)))->Particle());
  ThirdDaughter  = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(AliAnalysisTaskDoubleHypNucTree::GetLabel(stackN, PDGThirdDaughter)))->Particle());

  sublorentzsum = new TLorentzVector(0.,0.,0.,0.);
  particle1     = new TLorentzVector(0.,0.,0.,0.);
  particle2     = new TLorentzVector(0.,0.,0.,0.);
  particle3     = new TLorentzVector(0.,0.,0.,0.);		
	      
  particle1->SetXYZM(FirstDaughter->Px(),   FirstDaughter->Py(),  FirstDaughter->Pz(),  massFirstDaughter);
  particle2->SetXYZM(SecondDaughter->Px(),  SecondDaughter->Py(), SecondDaughter->Pz(), massSecondDaughter);
  particle3->SetXYZM(ThirdDaughter->Px(),   ThirdDaughter->Py(),  ThirdDaughter->Pz(),  massThirdDaughter);
  sublorentzsum->SetXYZM(0.,0.,0.,0.);
  *sublorentzsum = *particle1 + *particle2 + *particle3;
	  
  if(FirstDaughter->PdgCode() == PDGFirstDaughter){
    if(SecondDaughter->PdgCode() == PDGSecondDaughter){
      if(ThirdDaughter->PdgCode() == PDGThirdDaughter){

	dd[0] = ParticleMother->Xv() - FirstDaughter->Xv();
	dd[1] = ParticleMother->Yv() - FirstDaughter->Yv();
	dd[2] = ParticleMother->Zv() - FirstDaughter->Zv();
	  
	fPDGMother    = (Int_t)PDGMother;
	fChargeMother = sign;

	fmSubMother   = sublorentzsum->M();
	fptSubMother  = sublorentzsum->Pt();
	fySubMother   = sublorentzsum->Rapidity();
	fctSubMother  = (sublorentzsum->M()*TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2)))/sublorentzsum->P();
	      
	fptDaughter   = particle1->Pt();
	fyDaughter    = particle1->Rapidity();
	fptDaughter1  = particle2->Pt();
	fyDaughter1   = particle2->Rapidity();
	fptDaughter2  = particle3->Pt();
	fyDaughter2   = particle3->Rapidity();
	gTreeGen->Fill();
	AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      }
    }
  }
  if(FirstDaughter)  delete FirstDaughter;
  if(SecondDaughter) delete SecondDaughter;
  if(ThirdDaughter)  delete ThirdDaughter;
  if(sublorentzsum) delete sublorentzsum;
  if(particle1)     delete particle1;
  if(particle2)     delete particle2;
  if(particle3)     delete particle3;
}
// _________________________________________________ //
Double_t AliAnalysisTaskDoubleHypNucTree::Bethe(const AliESDtrack& track, Double_t mass, Int_t charge, Double_t* params){
  Double_t expected = charge*charge*AliExternalTrackParam::BetheBlochAleph(charge*track.GetInnerParam()->GetP()/mass,params[0],params[1],params[2],params[3],params[4]);
  Double_t sigma    = expected*params[5];
  if (TMath::IsNaN(expected)) return -999;
  return (track.GetTPCsignal() - expected) / sigma;
}
// _________________________________________________ //
Double_t AliAnalysisTaskDoubleHypNucTree::GetTOFSignal(const AliESDtrack& track){
  Float_t mass = 0;
  Float_t time = -1;
  Float_t beta = 0;
  Float_t gamma = 0;
  Float_t length = 0;
  Float_t time0 = 0;
  length = track.GetIntegratedLength();
  time0  = fPID->GetTOFResponse().GetStartTime(track.P());//fPID->GetTOFResponse().GetTimeZero();
  time   = track.GetTOFsignal() - time0;
  if (time > 0) {
    beta  = length / (2.99792457999999984e-02 * time);
    gamma = 1/TMath::Sqrt(1 - beta*beta);
    mass  = (track.GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
  }
  return mass;
}
// _________________________________________________ //
Double_t AliAnalysisTaskDoubleHypNucTree::GeoLength(const AliESDtrack& track) {
  Double_t deadZoneWidth = 3.0;
  Double_t lengthInActiveZone = track.GetLengthInActiveZone(1, deadZoneWidth, 220, track.GetESDEvent()->GetMagneticField(),0,0);
  return lengthInActiveZone;
}
// _________________________________________________ //
Float_t AliAnalysisTaskDoubleHypNucTree::GetInvPtDevFromBC(Int_t b, Int_t c) {
  //returns d(1/Pt) in c/GeV
  //in case of no gtu simulation -> return maximum 0.5
  if(b==0 && c==0) return 0.5;
  Int_t tmp = (((b & 0xfff) << 12) ^ 0x800000) - 0x800000;
  tmp += (c & 0xfff);
  Float_t invPtDev = tmp * 0.000001;
  return invPtDev;
}
// _________________________________________________ //
// __ reset variables __ //
void AliAnalysisTaskDoubleHypNucTree::ResetVals(TString mode){
  if(mode == "Event"){
    fPeriod             = -1;
    fMCtrue             = 0;
    frunnumber          = -1;
    fMagneticField      = -1;
    fCentrality         = -1;
    feventclass         = -1;
    MB                  = 0;
    HMV0                = 0;
    HMSPD               = 0;
    HNU                 = 0;
    HQU                 = 0;
  }
  if(mode == ""){
    fPrimVertexX        = -99;
    fPrimVertexY        = -99;
    fPrimVertexZ        = -99;
    fSecVertexX         = -99;
    fSecVertexY         = -99;
    fSecVertexZ         = -99;
    fTertVertexX        = -99;
    fTertVertexY        = -99;
    fTertVertexZ        = -99;
    fTertVertex2X       = -99;
    fTertVertex2Y       = -99;
    fTertVertex2Z       = -99;
    fTrigMB             = -99;
    fTrigHMV0           = -99;
    fTrigHMSPD          = -99;
    fTrigHNU            = -99;
    fTrigHQU            = -99;
    fDCA2B              = -99;
    fDCA2Bo             = -99;
    fDCA3B1             = -99;
    fDCA3B2             = -99;
    fDCA3B3             = -99;
    fPA                 = -99;
    fSubPA              = -99;
    fSubPA2             = -99;
    fSubPA3             = -99;
    fSubPA4             = -99;
    fDecAngle           = -99;
    farmalpha           = -99;
    farmpt              = -99;
    fPDGMother          = -99;
    fChargeMother       = -99;
    mctruth             = 0;
    fmMother            = -99;
    fmMother2           = -99;
    fEMother            = -99;
    fpxMother           = -99;
    fpyMother           = -99;
    fpzMother           = -99;
    fptMother           = -99;
    fpMother            = -99;
    fyMother            = -99;
    fctMother           = -99;
    fmSubMother         = -99;
    fESubMother         = -99;
    fpxSubMother        = -99;
    fpySubMother        = -99;
    fpzSubMother        = -99;
    fptSubMother        = -99;
    fpSubMother         = -99;
    fySubMother         = -99;
    fctSubMother        = -99;
    fEDaughter          = -99;
    fpDaughter          = -99;
    fptDaughter         = -99;
    fpxDaughter         = -99;
    fpyDaughter         = -99;
    fpzDaughter         = -99;
    fyDaughter          = -99;
    fdEdxDaughter       = -99;
    fdEdxSigmaDaughter  = -99;
    fDcaDaughter        = -99;
    fDcaDaughtero       = -99;
    fDcazDaughter       = -99;
    fDcaSecDaughter     = -99;
    fImParDaughter      = -99;
    fImParzDaughter     = -99;
    fImPar2Daughter     = -99;
    fImParz2Daughter    = -99;
    fNclsDaughter       = -99;
    fChi2Daughter       = -99;
    fNclsITSDaughter    = -99;
    fEtaDaughter        = -99;
    fPhiDaughter        = -99;
    fGeoLengthDaughter  = -99;
    fTOFSignalDaughter  = -99;
    fSigmaYXDaughter    = -99;
    fSigmaXYZDaughter   = -99;
    fSigmaZDaughter     = -99;
    fPtUncertDaughter   = -1;
    fTPCRefitDaughter   = -1;
    fITSRefitDaughter   = -1;
    fPropDCADaughter    = -1;
    fEDaughter1         = -99;
    fpDaughter1         = -99;
    fptDaughter1        = -99;
    fpxDaughter1        = -99;
    fpyDaughter1        = -99;
    fpzDaughter1        = -99;
    fyDaughter1         = -99;
    fdEdxDaughter1      = -99;
    fdEdxSigmaDaughter1 = -99;
    fDcaDaughter1       = -99;
    fDcaDaughter1o      = -99;
    fDcazDaughter1      = -99;
    fDcaSecDaughter1    = -99;
    fImParDaughter1     = -99;
    fImParzDaughter1    = -99;
    fImPar2Daughter1    = -99;
    fImParz2Daughter1   = -99;
    fNclsDaughter1      = -99;
    fChi2Daughter1      = -99;
    fNclsITSDaughter1   = -99;
    fEtaDaughter1       = -99;
    fPhiDaughter1       = -99;
    fGeoLengthDaughter1 = -99;
    fTOFSignalDaughter1 = -99;
    fSigmaYXDaughter1   = -99;
    fSigmaXYZDaughter1  = -99;
    fSigmaZDaughter1    = -99;
    fPtUncertDaughter1  = -1;
    fTPCRefitDaughter1  = -1;
    fITSRefitDaughter1  = -1;
    fPropDCADaughter1   = -1;
    fEDaughter2         = -99;
    fpDaughter2         = -99;
    fptDaughter2        = -99;
    fpxDaughter2        = -99;
    fpyDaughter2        = -99;
    fpzDaughter2        = -99;
    fyDaughter2         = -99;
    fdEdxDaughter2      = -99;
    fdEdxSigmaDaughter2 = -99;
    fDcaDaughter2       = -99;
    fDcaDaughter2o      = -99;
    fDcazDaughter2      = -99;
    fDcaSecDaughter2    = -99;
    fImParDaughter2     = -99;
    fImParzDaughter2    = -99;
    fImPar2Daughter2    = -99;
    fImParz2Daughter2   = -99;
    fNclsDaughter2      = -99;
    fChi2Daughter2      = -99;
    fNclsITSDaughter2   = -99;
    fEtaDaughter2       = -99;
    fPhiDaughter2       = -99;
    fGeoLengthDaughter2 = -99;
    fTOFSignalDaughter2 = -99;
    fSigmaYXDaughter2   = -99;
    fSigmaXYZDaughter2  = -99;
    fSigmaZDaughter2    = -99;
    fPtUncertDaughter2  = -1;
    fTPCRefitDaughter2  = -1;
    fITSRefitDaughter2  = -1;
    fPropDCADaughter2   = -1;
    fEDaughter3         = -99;
    fpDaughter3         = -99;
    fptDaughter3        = -99;
    fpxDaughter3        = -99;
    fpyDaughter3        = -99;
    fpzDaughter3        = -99;
    fyDaughter3         = -99;
    fdEdxDaughter3      = -99;
    fdEdxSigmaDaughter3 = -99;
    fDcaDaughter3       = -99;
    fDcaDaughter3o      = -99;
    fDcazDaughter3      = -99;
    fDcaSecDaughter3    = -99;
    fImParDaughter3     = -99;
    fImParzDaughter3    = -99;
    fImPar2Daughter3    = -99;
    fImParz2Daughter3   = -99;
    fNclsDaughter3      = -99;
    fChi2Daughter3      = -99;
    fNclsITSDaughter3   = -99;
    fEtaDaughter3       = -99;
    fPhiDaughter3       = -99;
    fGeoLengthDaughter3 = -99;
    fTOFSignalDaughter3 = -99;
    fSigmaYXDaughter3   = -99;
    fSigmaXYZDaughter3  = -99;
    fSigmaZDaughter3    = -99;
    fPtUncertDaughter3  = -1;
    fTPCRefitDaughter3  = -1;
    fITSRefitDaughter3  = -1;
    fPropDCADaughter3   = -1;
    fEDaughter4         = -99;
    fpDaughter4         = -99;
    fptDaughter4        = -99;
    fpxDaughter4        = -99;
    fpyDaughter4        = -99;
    fpzDaughter4        = -99;
    fyDaughter4         = -99;
    fDcaDaughter4       = -99;
    fDcazDaughter4      = -99;
    fDcaSecDaughter4    = -99;
    fImParDaughter4     = -99;
    fImParzDaughter4    = -99;
    fImPar2Daughter4    = -99;
    fImParz2Daughter4   = -99;
    fSigmaYXDaughter4   = -99;
    fSigmaXYZDaughter4  = -99;
    fSigmaZDaughter4    = -99;
    fPtUncertDaughter4  = -1;
    fPropDCADaughter4   = -1;
    fthetaP             = -99; 
    fthetaN             = -99;
    fdEdxSigmaPion      = -99;
    fdEdxSigmaDeuteron  = -99;
    fdEdxSigmaTriton    = -99;
    fdEdxSigmaAlpha     = -99;

  }
  return;
}
// _________________________________________________ //
// __ set Bethe-Bloch parameter __ //
void AliAnalysisTaskDoubleHypNucTree::SetBetheBlochParams(Int_t runNumber) {  
  if(runNumber == 170593){                                                                      
    fBetheParamsHe[0] = 3.18506;
    fBetheParamsHe[1] = 16.6883;
    fBetheParamsHe[2] = -0.200774;
    fBetheParamsHe[3] = 1.8954;
    fBetheParamsHe[4] = 1.33783;
    fBetheParamsHe[5] = 0.06;

    fBetheParamsT[0] = 3.67773;
    fBetheParamsT[1] = 14.5622;
    fBetheParamsT[2] = 5.72337;
    fBetheParamsT[3] = 1.91099;
    fBetheParamsT[4] = 2.16825;
    fBetheParamsT[5] = 0.06;
  }
  else if(runNumber < 246994) return;
  // __ 2015 PbPb __ //
  if(runNumber >= 244917 && runNumber <= 246994){
    fBetheParamsHe[0] = 2.45605;
    fBetheParamsHe[1] = 19.8067;
    fBetheParamsHe[2] = -0.77472;
    fBetheParamsHe[3] = 1.96279;
    fBetheParamsHe[4] = 0.172695;
    fBetheParamsHe[5] = 0.06;
      
    fBetheParamsT[0] = 2.32603;
    fBetheParamsT[1] = 19.2492;
    fBetheParamsT[2] = 30.7943;
    fBetheParamsT[3] = 2.1697;
    fBetheParamsT[4] = -8.11114;
    fBetheParamsT[5] = 0.06;
  }
  // __ 2018 PbPb pass 1__ //
  /* else if(runNumber >= 295581 && runNumber <= 297624){
     fBetheParamsT[0] = 0.669634;
     fBetheParamsT[1] = 53.1497;
     fBetheParamsT[2] =-1.32853e-08;
     fBetheParamsT[3] = 2.5775;
     fBetheParamsT[4] = 17.7607;
     fBetheParamsT[5] = 0.06;
      
     fBetheParamsHe[0] = 1.50582;
     fBetheParamsHe[1] = 33.7232;
     fBetheParamsHe[2] = -0.0923749;
     fBetheParamsHe[3] = 2.00901;
     fBetheParamsHe[4] = 2.28772;
     fBetheParamsHe[5] = 0.06;
     }*/
  // __ 2018 PbPb pass 3__ //
  else if(runNumber >= 295581 && runNumber <= 297624){
    fBetheParamsT[0] = 0.648689;
    fBetheParamsT[1] = 56.6706;
    fBetheParamsT[2] =-1.63243e-10;
    fBetheParamsT[3] = 2.46921;
    fBetheParamsT[4] = 16.8531;
    fBetheParamsT[5] = 0.06;
      
    fBetheParamsHe[0] = 1.70184;
    fBetheParamsHe[1] = 28.4426;
    fBetheParamsHe[2] = 3.21871e-12;
    fBetheParamsHe[3] = 2.06952;
    fBetheParamsHe[4] = 2.77971;
    fBetheParamsHe[5] = 0.06;
  }
  //__ 2016 pp __
  else if(runNumber >= 252235 && runNumber <= 264347 ) { 
    if(!fMCtrue) { // Data
      // LHC16 + LHC18
      // He3
      fBetheParamsHe[0] = 1.81085;
      fBetheParamsHe[1] = 29.4656;
      fBetheParamsHe[2] = 0.0458225;
      fBetheParamsHe[3] = 2.08689;
      fBetheParamsHe[4] = 2.28772;
      fBetheParamsHe[5] = 0.06;
      // Triton
      fBetheParamsT[0] = 1.58385;
      fBetheParamsT[1] = 25.8334;
      fBetheParamsT[2] = 0.00908038;
      fBetheParamsT[3] = 2.24769;
      fBetheParamsT[4] = 2.87755;
      fBetheParamsT[5] = 0.06;
    } else { // MC
      if (runNumber >= 262424 || runNumber <= 256418 ) {
	//LHC18a2b (->LHC16)
	// He3
	fBetheParamsHe[0] = 3.05245;
	fBetheParamsHe[1] = 15.7252;
	fBetheParamsHe[2] = -0.00453331;
	fBetheParamsHe[3] = 2.17241;
	fBetheParamsHe[4] = 2.88422;
	fBetheParamsHe[5] = 0.0834274;
	// Triton
	fBetheParamsT[0] = 2.74259;
	fBetheParamsT[1] = 18.3295;
	fBetheParamsT[2] = 5.91594;
	fBetheParamsT[3] = 1.93471;
	fBetheParamsT[4] = 0.292147;
	fBetheParamsT[5] = 0.0728241;
      }
      if (runNumber >= 256941 && runNumber <= 258537 ) {
	// LHC18a2b2 (LHC16k)
	// He3
	fBetheParamsHe[0] = 2.80527;
	fBetheParamsHe[1] = 14.2379;
	fBetheParamsHe[2] = 0.0232811;
	fBetheParamsHe[3] = 2.11464;
	fBetheParamsHe[4] = 1.615;
	fBetheParamsHe[5] = 0.0815227;
	// Triton
	fBetheParamsT[0] = 1.31603;
	fBetheParamsT[1] = 36.1798;
	fBetheParamsT[2] = 493.036;
	fBetheParamsT[3] = 2.10841;
	fBetheParamsT[4] = 7.43391;
	fBetheParamsT[5] = 0.0769041;
      }
      if (runNumber >= 258962 && runNumber <= 259888 ) {
	//LHC18a2b3 (->LHC16l)
	// He3
	fBetheParamsHe[0] = 2.80121;
	fBetheParamsHe[1] = 14.2397;
	fBetheParamsHe[2] = 0.0100894;
	fBetheParamsHe[3] = 2.10396;
	fBetheParamsHe[4] = 1.41608;
	fBetheParamsHe[5] = 0.0817429;
	// Triton
	fBetheParamsT[0] = 4.80597;
	fBetheParamsT[1] = 13.8813;
	fBetheParamsT[2] = 189.651;
	fBetheParamsT[3] = 2.05969;
	fBetheParamsT[4] = 4.38013;
	fBetheParamsT[5] = 0.077593;
      }
    }
  }
  // __ 2017 pp __ //
  if (runNumber >= 270581 && runNumber <= 282704) { 
    if(!fMCtrue) {
      //LHC17
      // He3
      fBetheParamsHe[0] = 3.20025;
      fBetheParamsHe[1] = 16.4971;
      fBetheParamsHe[2] = -0.0116571;
      fBetheParamsHe[3] = 2.3152;
      fBetheParamsHe[4] = 3.11135;
      fBetheParamsHe[5] = 0.06;
      // Triton
      fBetheParamsT[0] = 1.69461;
      fBetheParamsT[1] = 27.6917;
      fBetheParamsT[2] = 0.372214;
      fBetheParamsT[3] = 2.05305;
      fBetheParamsT[4] = -1.25037;
      fBetheParamsT[5] = 0.06;
    } else {
      // LHC18a2a (->LHC17)
      // He3
      fBetheParamsHe[0] = 3.12796;
      fBetheParamsHe[1] = 16.1359;
      fBetheParamsHe[2] = -0.00682978;
      fBetheParamsHe[3] = 2.26624;
      fBetheParamsHe[4] = 2.58652;
      fBetheParamsHe[5] = 0.0847009;
      // Triton
      fBetheParamsT[0] = 2.8303;
      fBetheParamsT[1] = 15.4337;
      fBetheParamsT[2] = 3.18352;
      fBetheParamsT[3] = 2.20975;
      fBetheParamsT[4] = 0.218244;
      fBetheParamsT[5] = 0.0780191;
    }
  }
  // __ 2018 pp __ //
  if (runNumber >= 285009 && runNumber <= 294925) { 
    if(!fMCtrue) {
      // LHC16 + LHC18
      // He3
      fBetheParamsHe[0] = 1.81085;
      fBetheParamsHe[1] = 29.4656;
      fBetheParamsHe[2] = 0.0458225;
      fBetheParamsHe[3] = 2.08689;
      fBetheParamsHe[4] = 2.28772;
      fBetheParamsHe[5] = 0.06;
      // Triton
      fBetheParamsT[0] = 1.58385;
      fBetheParamsT[1] = 25.8334;
      fBetheParamsT[2] = 0.00908038;
      fBetheParamsT[3] = 2.24769;
      fBetheParamsT[4] = 2.87755;
      fBetheParamsT[5] = 0.06;
    } else {
      //LHC18a2d (->LHC18)
      // He3
      fBetheParamsHe[0] = 3.07104;
      fBetheParamsHe[1] = 15.8085;
      fBetheParamsHe[2] = 0.0150992;
      fBetheParamsHe[3] = 2.13909;
      fBetheParamsHe[4] = 2.59495;
      fBetheParamsHe[5] = 0.0865179;
      // Triton
      fBetheParamsT[0] = 2.54486;
      fBetheParamsT[1] = 17.1203;
      fBetheParamsT[2] = -0.0452007;
      fBetheParamsT[3] = 2.00988;
      fBetheParamsT[4] = 0.849292;
      fBetheParamsT[5] = 0.0768715;
    }
  }
  return;
}
