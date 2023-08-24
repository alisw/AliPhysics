//--- Task for investigation of the {}^{4}_{#Lambda#Lambda}H ---
//---     Author: Janik Ditzel; janik.ditzel@cern.ch     ---


#include "AliAnalysisTaskDoubleHypNucTreeLS.h"

ClassImp(AliAnalysisTaskDoubleHypNucTreeLS)
// _________________________________________________ //
// Default Constructor
AliAnalysisTaskDoubleHypNucTreeLS::AliAnalysisTaskDoubleHypNucTreeLS()
:AliAnalysisTaskSE("AliAnalysisTaskDoubleHypNucTreeLS"),
  fPIDCheckOnly(kFALSE), fInputHandler(0), fPID(0), fESDevent(0), mcEvent(0), fStack(), fEventCuts(), fAODevent(), fevent(), eventHeader(),
//
  fTriggerMask(), MB(0), HMV0(0), HMSPD(0), HNU(0), HQU(0), Central(0), SemiCentral(0),
//
  fTree(0), fTreeKF(0), gTree(0), gTreeKF(0), hTree(0), hTreeKF(0), iTree(0), iTreeKF(0), fTreeGen(0), gTreeGen(0),
//
  kVariante(1), kStandardReco(0), kKFReco(0), kSkip3LHAnalysis(0), kSkip4LHAnalysis(0), kAOD(0), kESD(0),
  //
  He4PosArray(), He4NegArray(), He3PosArray(), PPosArray(), PiPosArray(), PiPosSecArray(), He3NegArray(), PNegArray(), PiNegArray(), PiNegSecArray(), He3PosCounter(), He3NegCounter(), PPosCounter(), PNegCounter(), PiPosCounter(), PiNegCounter(), PiPosSecCounter(), PiNegSecCounter(), He4PosCounter(), He4NegCounter(),
//
  fBetheParamsHe(), fBetheParamsT(),
//
  fMCtrue(0), konlyBG(0), konlySig(0), kMCPIDCheck(0), fPrimary4LHe(0), stackN(0), fisExcited(0),
//
  vertex(), primVertex(), secvertexer(), secVertex(), tertvertexer(), tertVertex(), PrimVertex(), SecVertex(), TertVertex(),
//
  track1(), track2(), track3(), track4(), exTrack(), exTrack1(), exTrack2(), exTrack3(), exTrack4(), trkArray(), trkArray1(),
//
  kDCATracksCut(0), kTrackPtUncertCut(0), kPointingAngleCut(0), kPIDSelecCut(0), kMin4LLHMass(0), kMax4LLHMass(0), kMin4LHeMass(0), kMax4LHeMass(0), kMin4LHMass(0), kMax4LHMass(0), kMin3LHMass(0), kMax3LHMass(0),
//
  kMagF(0),
//
  cov(), cov0(), cov1(), cov2(), pxpypz(), xyz(), sign(0), dd(), xthiss(), xpp(), h(), fthetaP(-99), fthetaN(-99),
//  
  PDGCodeMother(), ParticleMother(NULL), FirstDaughter(NULL), SecondDaughter(NULL), ThirdDaughter(NULL), FourthDaughter(NULL), label1(), labelMother1(), labelGrandMother1(), ParticleMother1(NULL), ParticleGrandMother1(NULL), label2(), labelMother2(), labelGrandMother2(), ParticleMother2(NULL), ParticleGrandMother2(NULL), label3(), labelMother3(), labelGrandMother3(), ParticleMother3(NULL), ParticleGrandMother3(NULL), label4(), labelMother4(), labelGrandMother4(), ParticleMother4(NULL), ParticleGrandMother4(NULL), CheckParticle(NULL),
//
  lorentzsum(NULL), lorentzsum2(NULL), sublorentzsum(NULL), sublorentzsum2(NULL), particle1(NULL), particle2(NULL), particle3(NULL), particle4(NULL),
//Histos
  fHistogramList(NULL), fHistdEdx(0), fHistNumEvents(0), fHistTrigger(0), fHistCentrality1(), fHistTrigger1(), fHistCentrality2(), fHistTrigger2(), fHistNtracks1(), fHistNtracks2(),
//saved in tree
  fEventID(-1), fParticleID(-1), fCentrality(-99), frunnumber(-99), fPeriod(0), fMagneticField(0), feventclass(0), fmctruth(0), fDecayChannel(0), fRecoMethod(0), fisWOTrackVertex(), fTriggerString(),
  fTrigMB(-99), fTrigHMV0(-99), fTrigHMSPD(-99), fTrigHNU(0), fTrigHQU(0), fTrigkCentral(0), fTrigkSemiCentral(0),
//
  fisOnlineV0_13(0), fisOnlineV0_14(0), fisOnlineV0_23(0), fisOnlineV0_24(0),
//
  fDCA2B(-99), fDCA2Bo(-99), fDCA3B1(-99), fDCA3B2(-99), fDCA3B3(-99), fDecAngle(-99), farmalpha(-99), farmpt(-99), fPA(-99), fSubPA(-99), fSubPA2(-99),
//
  fPrimVertexX(-99), fPrimVertexY(-99), fPrimVertexZ(-99), fSecVertexX(-99), fSecVertexY(-99), fSecVertexZ(-99), fTertVertexX(-99), fTertVertexY(-99), fTertVertexZ(-99), fTertVertChi2(-99), fTertVertNDF(-99), fSecVertChi2(-99), fSecVertNDF(-99), fPrimVertChi2(-99), fPrimVertNDF(-99), fV0VertexX_13(-99), fV0VertexY_13(-99), fV0VertexZ_13(-99), fV0VertexX_14(-99), fV0VertexY_14(-99), fV0VertexZ_14(-99), fV0VertexX_23(-99), fV0VertexY_23(-99), fV0VertexZ_23(-99), fV0VertexX_24(-99), fV0VertexY_24(-99), fV0VertexZ_24(-99),
//
  fPDGMother(-99), fChargeMother(-99), fmMother(-99), fmMother2(-99), fEMother(-99), fpxMother(-99), fpyMother(-99), fpzMother(-99), fptMother(-99), fpMother(-99), fyMother(-99), fctMother(-99), fmMotherCheck(-99), fmSubMotherCheck(-99),
//
  fmSubMother(-99), fESubMother(-99), fpxSubMother(-99), fpySubMother(-99), fpzSubMother(-99), fptSubMother(-99), fpSubMother(-99), fySubMother(-99), fctSubMother(-99),
//
  fEDaughter(-99), fpDaughter(-99), fptDaughter(-99), fpxDaughter(-99), fpyDaughter(-99), fpzDaughter(-99), fyDaughter(-99), fdEdxDaughter(-99), fdEdxSigmaDaughter(-99), fDcaDaughter(-99), fDcaDaughtero(-99), fDcazDaughter(-99), fDcaSecDaughter(-99), fImParDaughter(-99), fImParzDaughter(-99), fNclsDaughter(-99), fChi2Daughter(-99), fNclsITSDaughter(-99), fEtaDaughter(-99), fPhiDaughter(-99), fGeoLengthDaughter(-99), fTOFSignalDaughter(-99), fSigmaYXDaughter(-99), fSigmaXYZDaughter(-99), fSigmaZDaughter(-99), fPtUncertDaughter(-99), fTPCRefitDaughter(-99), fITSRefitDaughter(-99), fPropDCADaughter(0),
//
  fEDaughter1(-99), fpDaughter1(-99), fptDaughter1(-99), fpxDaughter1(-99), fpyDaughter1(-99), fpzDaughter1(-99), fyDaughter1(-99), fdEdxDaughter1(-99), fdEdxSigmaDaughter1(-99), fDcaDaughter1(-99), fDcaDaughter1o(-99), fDcazDaughter1(-99), fDcaSecDaughter1(-99), fImParDaughter1(-99), fImParzDaughter1(-99), fNclsDaughter1(-99), fChi2Daughter1(-99), fNclsITSDaughter1(-99), fEtaDaughter1(-99), fPhiDaughter1(-99), fGeoLengthDaughter1(-99), fTOFSignalDaughter1(-99), fSigmaYXDaughter1(-99), fSigmaXYZDaughter1(-99), fSigmaZDaughter1(-99), fPtUncertDaughter1(-99), fTPCRefitDaughter1(-99), fITSRefitDaughter1(-99), fPropDCADaughter1(0),
//
  fEDaughter2(-99), fpDaughter2(-99), fptDaughter2(-99), fpxDaughter2(-99), fpyDaughter2(-99), fpzDaughter2(-99), fyDaughter2(-99), fdEdxDaughter2(-99), fdEdxSigmaDaughter2(-99), fDcaDaughter2(-99), fDcaDaughter2o(-99), fDcazDaughter2(-99), fDcaSecDaughter2(-99), fImParDaughter2(-99), fImParzDaughter2(-99), fNclsDaughter2(-99), fChi2Daughter2(-99), fNclsITSDaughter2(-99), fEtaDaughter2(-99), fPhiDaughter2(-99), fGeoLengthDaughter2(-99), fTOFSignalDaughter2(-99), fSigmaYXDaughter2(-99), fSigmaXYZDaughter2(-99), fSigmaZDaughter2(-99), fPtUncertDaughter2(-99), fTPCRefitDaughter2(-99), fITSRefitDaughter2(-99), fPropDCADaughter2(0),
//
  fEDaughter3(-99), fpDaughter3(-99), fptDaughter3(-99), fpxDaughter3(-99), fpyDaughter3(-99), fpzDaughter3(-99), fyDaughter3(-99), fdEdxDaughter3(-99), fdEdxSigmaDaughter3(-99), fDcaDaughter3(-99), fDcaDaughter3o(-99), fDcazDaughter3(-99), fDcaSecDaughter3(-99), fImParDaughter3(-99), fImParzDaughter3(-99), fNclsDaughter3(-99), fChi2Daughter3(-99), fNclsITSDaughter3(-99), fEtaDaughter3(-99), fPhiDaughter3(-99), fGeoLengthDaughter3(-99), fTOFSignalDaughter3(-99), fSigmaYXDaughter3(-99), fSigmaXYZDaughter3(-99), fSigmaZDaughter3(-99), fPtUncertDaughter3(-99), fTPCRefitDaughter3(-99), fITSRefitDaughter3(-99), fPropDCADaughter3(0),
//
  fEDaughter4(-99), fpDaughter4(-99), fptDaughter4(-99), fpxDaughter4(-99), fpyDaughter4(-99), fpzDaughter4(-99), fyDaughter4(-99), fDcaDaughter4(-99), fDcazDaughter4(-99), fDcaSecDaughter4(-99), fImParDaughter4(-99), fImParzDaughter4(-99), fSigmaYXDaughter4(-99), fSigmaXYZDaughter4(-99), fSigmaZDaughter4(-99), fPtUncertDaughter4(-99), fPropDCADaughter4(0),
//
  fCovMatrixTrack(), fCovMatrixTrack1(), fCovMatrixTrack2(), fCovMatrixTrack3(), fCovMatrixTrack4(), fTrackPar(), fTrackPar1(), fTrackPar2(), fTrackPar3(), fTrackPar4(),
//
  fTrackPIDDaughter(-99), fTrackPIDDaughter1(-99), fTrackPIDDaughter2(-99), fTrackPIDDaughter3(-99),
//
  fITSLayer1Daughter(0), fITSLayer2Daughter(0), fITSLayer3Daughter(0), fITSLayer4Daughter(0), fITSLayer5Daughter(0), fITSLayer6Daughter(0), fITSLayer1Daughter1(0), fITSLayer2Daughter1(0), fITSLayer3Daughter1(0), fITSLayer4Daughter1(0), fITSLayer5Daughter1(0), fITSLayer6Daughter1(0), fITSLayer1Daughter2(0), fITSLayer2Daughter2(0), fITSLayer3Daughter2(0), fITSLayer4Daughter2(0), fITSLayer5Daughter2(0), fITSLayer6Daughter2(0), fITSLayer1Daughter3(0), fITSLayer2Daughter3(0), fITSLayer3Daughter3(0), fITSLayer4Daughter3(0), fITSLayer5Daughter3(0), fITSLayer6Daughter3(0),
//
  fdEdxSigmaPion(-99), fdEdxSigmaDeuteron(-99), fdEdxSigmaTriton(-99), fdEdxSigmaAlpha(-99),
//
  PrimVertexKF(), SecVertexKF(), TertVertexKF(), KFtertVertex(), KFsecVertex(), SecVertexErrKF(), TertVertexErrKF(), fPrimVertexXKF(-99), fPrimVertexYKF(-99), fPrimVertexZKF(-99), fSecVertexXKF(-99), fSecVertexYKF(-99), fSecVertexZKF(-99), fTertVertexXKF(-99), fTertVertexYKF(-99), fTertVertexZKF(-99), fPrimVertexXErrKF(-99), fPrimVertexYErrKF(-99), fPrimVertexZErrKF(-99), fSecVertexXErrKF(-99), fSecVertexYErrKF(-99), fSecVertexZErrKF(-99), fTertVertexXErrKF(-99), fTertVertexYErrKF(-99), fTertVertexZErrKF(-99), fTertVertChi2KF(-99), fTertVertNDFKF(-99), fSecVertChi2KF(-99), fSecVertNDFKF(-99), fPrimVertChi2KF(-99), fPrimVertNDFKF(-99),
//
  KFtrack1(), KFtrack2(), KFtrack3(), KFtrack4(), primKFVertex(),
//
  fPAKF(-99), fSubPAKF(-99), fSubPA2KF(-99),
//
  fmMotherKF(-99), fmMotherErrKF(-99), fEMotherKF(-99), fEMotherErrKF(-99), fpxMotherKF(-99), fpxMotherErrKF(-99), fpyMotherKF(-99), fpyMotherErrKF(-99), fpzMotherKF(-99), fpzMotherErrKF(-99), fptMotherKF(-99), fptMotherErrKF(-99), fpMotherKF(-99), fpMotherErrKF(-99), fyMotherKF(-99), fctMotherKF(-99), fctMotherErrKF(-99),
//
  fmSubMotherKF(-99), fmSubMotherErrKF(-99), fESubMotherKF(-99), fESubMotherErrKF(-99), fpxSubMotherKF(-99), fpxSubMotherErrKF(-99), fpySubMotherKF(-99), fpySubMotherErrKF(-99), fpzSubMotherKF(-99), fpzSubMotherErrKF(-99), fptSubMotherKF(-99), fptSubMotherErrKF(-99), fpSubMotherKF(-99), fpSubMotherErrKF(-99), fySubMotherKF(-99), fctSubMotherKF(-99), fctSubMotherErrKF(-99),
//  
  fEDaughterKF(-99), fLabelDaughterKF(-99), fptDaughterKF(-99), fpxDaughterKF(-99), fpyDaughterKF(-99), fpzDaughterKF(-99), fyDaughterKF(-99),
  fEDaughter1KF(-99), fLabelDaughter1KF(-99), fptDaughter1KF(-99), fpxDaughter1KF(-99), fpyDaughter1KF(-99), fpzDaughter1KF(-99), fyDaughter1KF(-99),
  fEDaughter2KF(-99), fLabelDaughter2KF(-99), fptDaughter2KF(-99), fpxDaughter2KF(-99), fpyDaughter2KF(-99), fpzDaughter2KF(-99), fyDaughter2KF(-99),
  fEDaughter3KF(-99), fLabelDaughter3KF(-99), fptDaughter3KF(-99), fpxDaughter3KF(-99), fpyDaughter3KF(-99), fpzDaughter3KF(-99), fyDaughter3KF(-99),
//  
  fptDaughterUnProp(-99), fpxDaughterUnProp(-99), fpyDaughterUnProp(-99), fpzDaughterUnProp(-99), 
  fptDaughter1UnProp(-99), fpxDaughter1UnProp(-99), fpyDaughter1UnProp(-99), fpzDaughter1UnProp(-99), 
  fptDaughter2UnProp(-99), fpxDaughter2UnProp(-99), fpyDaughter2UnProp(-99), fpzDaughter2UnProp(-99), 
  fptDaughter3UnProp(-99), fpxDaughter3UnProp(-99), fpyDaughter3UnProp(-99), fpzDaughter3UnProp(-99), 
//
  fDcaDaughterXYKF(-99), fDcaDaughter1XYKF(-99), fDcaDaughter2XYKF(-99), fDcaDaughter3XYKF(-99), fDcaDaughter4XYKF(-99), fDcaSecDaughterXYKF(-99), fDcaSecDaughter1XYKF(-99), fDcaSecDaughter2XYKF(-99), fDcaSecDaughter3XYKF(-99), fDcaSecDaughter4XYKF(-99), fDcaDaughterZKF(-99), fDcaDaughter1ZKF(-99), fDcaDaughter2ZKF(-99), fDcaDaughter3ZKF(-99), fDcaDaughter4ZKF(-99), fDcaSecDaughterZKF(-99), fDcaSecDaughter1ZKF(-99), fDcaSecDaughter2ZKF(-99), fDcaSecDaughter3ZKF(-99), fDcaSecDaughter4ZKF(-99),
//
  fDCA2BXYKF(-99), fDCA3B1XYKF(-99), fDCA3B2XYKF(-99), fDCA3B3XYKF(-99), fDCA2BZKF(-99), fDCA3B1ZKF(-99), fDCA3B2ZKF(-99), fDCA3B3ZKF(-99)
{

}
// _________________________________________________ //
// Constructor
AliAnalysisTaskDoubleHypNucTreeLS::AliAnalysisTaskDoubleHypNucTreeLS(const char* name)
  :AliAnalysisTaskSE(name),
   fPIDCheckOnly(kFALSE), fInputHandler(0), fPID(0), fESDevent(0), mcEvent(0), fStack(), fEventCuts(), fAODevent(), fevent(), eventHeader(),
   //
   fTriggerMask(), MB(0), HMV0(0), HMSPD(0), HNU(0), HQU(0), Central(0), SemiCentral(0),
   //
   fTree(0), fTreeKF(0), gTree(0), gTreeKF(0), hTree(0), hTreeKF(0), iTree(0), iTreeKF(0), fTreeGen(0), gTreeGen(0),
   //
   kVariante(1), kStandardReco(0), kKFReco(0), kSkip3LHAnalysis(0), kSkip4LHAnalysis(0), kAOD(0), kESD(0),
   //
   He4PosArray(), He4NegArray(), He3PosArray(), PPosArray(), PiPosArray(), PiPosSecArray(), He3NegArray(), PNegArray(), PiNegArray(), PiNegSecArray(), He3PosCounter(), He3NegCounter(), PPosCounter(), PNegCounter(), PiPosCounter(), PiNegCounter(), PiPosSecCounter(), PiNegSecCounter(), He4PosCounter(), He4NegCounter(),
   //
   fBetheParamsHe(), fBetheParamsT(),
   //
   fMCtrue(0), konlyBG(0), konlySig(0), kMCPIDCheck(0), fPrimary4LHe(0), stackN(0), fisExcited(0),
   //
   vertex(), primVertex(), secvertexer(), secVertex(), tertvertexer(), tertVertex(), PrimVertex(), SecVertex(), TertVertex(),
   //
   track1(), track2(), track3(), track4(), exTrack(), exTrack1(), exTrack2(), exTrack3(), exTrack4(), trkArray(), trkArray1(),
   //
   kDCATracksCut(0), kTrackPtUncertCut(0), kPointingAngleCut(0), kPIDSelecCut(0), kMin4LLHMass(0), kMax4LLHMass(0), kMin4LHeMass(0), kMax4LHeMass(0), kMin4LHMass(0), kMax4LHMass(0), kMin3LHMass(0), kMax3LHMass(0),
   //
   kMagF(0),
   //
   cov(), cov0(), cov1(), cov2(), pxpypz(), xyz(), sign(0), dd(), xthiss(), xpp(), h(), fthetaP(-99), fthetaN(-99),
   //  
   PDGCodeMother(), ParticleMother(NULL), FirstDaughter(NULL), SecondDaughter(NULL), ThirdDaughter(NULL), FourthDaughter(NULL), label1(), labelMother1(), labelGrandMother1(), ParticleMother1(NULL), ParticleGrandMother1(NULL), label2(), labelMother2(), labelGrandMother2(), ParticleMother2(NULL), ParticleGrandMother2(NULL), label3(), labelMother3(), labelGrandMother3(), ParticleMother3(NULL), ParticleGrandMother3(NULL), label4(), labelMother4(), labelGrandMother4(), ParticleMother4(NULL), ParticleGrandMother4(NULL), CheckParticle(NULL),
   //
   lorentzsum(NULL), lorentzsum2(NULL), sublorentzsum(NULL), sublorentzsum2(NULL), particle1(NULL), particle2(NULL), particle3(NULL), particle4(NULL),
   //Histos
   fHistogramList(NULL), fHistdEdx(0), fHistNumEvents(0), fHistTrigger(0), fHistCentrality1(), fHistTrigger1(), fHistCentrality2(), fHistTrigger2(), fHistNtracks1(), fHistNtracks2(),
   //saved in tree
   fEventID(-1), fParticleID(-1), fCentrality(-99), frunnumber(-99), fPeriod(0), fMagneticField(0), feventclass(0), fmctruth(0), fDecayChannel(0), fRecoMethod(0), fisWOTrackVertex(), fTriggerString(),
   fTrigMB(-99), fTrigHMV0(-99), fTrigHMSPD(-99), fTrigHNU(0), fTrigHQU(0), fTrigkCentral(0), fTrigkSemiCentral(0),
   //
   fisOnlineV0_13(0), fisOnlineV0_14(0), fisOnlineV0_23(0), fisOnlineV0_24(0),
   //
   fDCA2B(-99), fDCA2Bo(-99), fDCA3B1(-99), fDCA3B2(-99), fDCA3B3(-99), fDecAngle(-99), farmalpha(-99), farmpt(-99), fPA(-99), fSubPA(-99), fSubPA2(-99),
   //
   fPrimVertexX(-99), fPrimVertexY(-99), fPrimVertexZ(-99), fSecVertexX(-99), fSecVertexY(-99), fSecVertexZ(-99), fTertVertexX(-99), fTertVertexY(-99), fTertVertexZ(-99), fTertVertChi2(-99), fTertVertNDF(-99), fSecVertChi2(-99), fSecVertNDF(-99), fPrimVertChi2(-99), fPrimVertNDF(-99), fV0VertexX_13(-99), fV0VertexY_13(-99), fV0VertexZ_13(-99), fV0VertexX_14(-99), fV0VertexY_14(-99), fV0VertexZ_14(-99), fV0VertexX_23(-99), fV0VertexY_23(-99), fV0VertexZ_23(-99), fV0VertexX_24(-99), fV0VertexY_24(-99), fV0VertexZ_24(-99),
   //
   fPDGMother(-99), fChargeMother(-99), fmMother(-99), fmMother2(-99), fEMother(-99), fpxMother(-99), fpyMother(-99), fpzMother(-99), fptMother(-99), fpMother(-99), fyMother(-99), fctMother(-99), fmMotherCheck(-99), fmSubMotherCheck(-99),
   //
   fmSubMother(-99), fESubMother(-99), fpxSubMother(-99), fpySubMother(-99), fpzSubMother(-99), fptSubMother(-99), fpSubMother(-99), fySubMother(-99), fctSubMother(-99),
   //
   fEDaughter(-99), fpDaughter(-99), fptDaughter(-99), fpxDaughter(-99), fpyDaughter(-99), fpzDaughter(-99), fyDaughter(-99), fdEdxDaughter(-99), fdEdxSigmaDaughter(-99), fDcaDaughter(-99), fDcaDaughtero(-99), fDcazDaughter(-99), fDcaSecDaughter(-99), fImParDaughter(-99), fImParzDaughter(-99), fNclsDaughter(-99), fChi2Daughter(-99), fNclsITSDaughter(-99), fEtaDaughter(-99), fPhiDaughter(-99), fGeoLengthDaughter(-99), fTOFSignalDaughter(-99), fSigmaYXDaughter(-99), fSigmaXYZDaughter(-99), fSigmaZDaughter(-99), fPtUncertDaughter(-99), fTPCRefitDaughter(-99), fITSRefitDaughter(-99), fPropDCADaughter(0),
   //
   fEDaughter1(-99), fpDaughter1(-99), fptDaughter1(-99), fpxDaughter1(-99), fpyDaughter1(-99), fpzDaughter1(-99), fyDaughter1(-99), fdEdxDaughter1(-99), fdEdxSigmaDaughter1(-99), fDcaDaughter1(-99), fDcaDaughter1o(-99), fDcazDaughter1(-99), fDcaSecDaughter1(-99), fImParDaughter1(-99), fImParzDaughter1(-99), fNclsDaughter1(-99), fChi2Daughter1(-99), fNclsITSDaughter1(-99), fEtaDaughter1(-99), fPhiDaughter1(-99), fGeoLengthDaughter1(-99), fTOFSignalDaughter1(-99), fSigmaYXDaughter1(-99), fSigmaXYZDaughter1(-99), fSigmaZDaughter1(-99), fPtUncertDaughter1(-99), fTPCRefitDaughter1(-99), fITSRefitDaughter1(-99), fPropDCADaughter1(0),
   //
   fEDaughter2(-99), fpDaughter2(-99), fptDaughter2(-99), fpxDaughter2(-99), fpyDaughter2(-99), fpzDaughter2(-99), fyDaughter2(-99), fdEdxDaughter2(-99), fdEdxSigmaDaughter2(-99), fDcaDaughter2(-99), fDcaDaughter2o(-99), fDcazDaughter2(-99), fDcaSecDaughter2(-99), fImParDaughter2(-99), fImParzDaughter2(-99), fNclsDaughter2(-99), fChi2Daughter2(-99), fNclsITSDaughter2(-99), fEtaDaughter2(-99), fPhiDaughter2(-99), fGeoLengthDaughter2(-99), fTOFSignalDaughter2(-99), fSigmaYXDaughter2(-99), fSigmaXYZDaughter2(-99), fSigmaZDaughter2(-99), fPtUncertDaughter2(-99), fTPCRefitDaughter2(-99), fITSRefitDaughter2(-99), fPropDCADaughter2(0),
   //
   fEDaughter3(-99), fpDaughter3(-99), fptDaughter3(-99), fpxDaughter3(-99), fpyDaughter3(-99), fpzDaughter3(-99), fyDaughter3(-99), fdEdxDaughter3(-99), fdEdxSigmaDaughter3(-99), fDcaDaughter3(-99), fDcaDaughter3o(-99), fDcazDaughter3(-99), fDcaSecDaughter3(-99), fImParDaughter3(-99), fImParzDaughter3(-99), fNclsDaughter3(-99), fChi2Daughter3(-99), fNclsITSDaughter3(-99), fEtaDaughter3(-99), fPhiDaughter3(-99), fGeoLengthDaughter3(-99), fTOFSignalDaughter3(-99), fSigmaYXDaughter3(-99), fSigmaXYZDaughter3(-99), fSigmaZDaughter3(-99), fPtUncertDaughter3(-99), fTPCRefitDaughter3(-99), fITSRefitDaughter3(-99), fPropDCADaughter3(0),
   //
   fEDaughter4(-99), fpDaughter4(-99), fptDaughter4(-99), fpxDaughter4(-99), fpyDaughter4(-99), fpzDaughter4(-99), fyDaughter4(-99), fDcaDaughter4(-99), fDcazDaughter4(-99), fDcaSecDaughter4(-99), fImParDaughter4(-99), fImParzDaughter4(-99), fSigmaYXDaughter4(-99), fSigmaXYZDaughter4(-99), fSigmaZDaughter4(-99), fPtUncertDaughter4(-99), fPropDCADaughter4(0),
   //
   fCovMatrixTrack(), fCovMatrixTrack1(), fCovMatrixTrack2(), fCovMatrixTrack3(), fCovMatrixTrack4(), fTrackPar(), fTrackPar1(), fTrackPar2(), fTrackPar3(), fTrackPar4(),
   //
   fTrackPIDDaughter(-99), fTrackPIDDaughter1(-99), fTrackPIDDaughter2(-99), fTrackPIDDaughter3(-99),
   //
   fITSLayer1Daughter(0), fITSLayer2Daughter(0), fITSLayer3Daughter(0), fITSLayer4Daughter(0), fITSLayer5Daughter(0), fITSLayer6Daughter(0), fITSLayer1Daughter1(0), fITSLayer2Daughter1(0), fITSLayer3Daughter1(0), fITSLayer4Daughter1(0), fITSLayer5Daughter1(0), fITSLayer6Daughter1(0), fITSLayer1Daughter2(0), fITSLayer2Daughter2(0), fITSLayer3Daughter2(0), fITSLayer4Daughter2(0), fITSLayer5Daughter2(0), fITSLayer6Daughter2(0), fITSLayer1Daughter3(0), fITSLayer2Daughter3(0), fITSLayer3Daughter3(0), fITSLayer4Daughter3(0), fITSLayer5Daughter3(0), fITSLayer6Daughter3(0),
   //
   fdEdxSigmaPion(-99), fdEdxSigmaDeuteron(-99), fdEdxSigmaTriton(-99), fdEdxSigmaAlpha(-99),
   //
   PrimVertexKF(), SecVertexKF(), TertVertexKF(), KFtertVertex(), KFsecVertex(), SecVertexErrKF(), TertVertexErrKF(), fPrimVertexXKF(-99), fPrimVertexYKF(-99), fPrimVertexZKF(-99), fSecVertexXKF(-99), fSecVertexYKF(-99), fSecVertexZKF(-99), fTertVertexXKF(-99), fTertVertexYKF(-99), fTertVertexZKF(-99), fPrimVertexXErrKF(-99), fPrimVertexYErrKF(-99), fPrimVertexZErrKF(-99), fSecVertexXErrKF(-99), fSecVertexYErrKF(-99), fSecVertexZErrKF(-99), fTertVertexXErrKF(-99), fTertVertexYErrKF(-99), fTertVertexZErrKF(-99), fTertVertChi2KF(-99), fTertVertNDFKF(-99), fSecVertChi2KF(-99), fSecVertNDFKF(-99), fPrimVertChi2KF(-99), fPrimVertNDFKF(-99),
   //
   KFtrack1(), KFtrack2(), KFtrack3(), KFtrack4(), primKFVertex(),
   //
   fPAKF(-99), fSubPAKF(-99), fSubPA2KF(-99),
   //
   fmMotherKF(-99), fmMotherErrKF(-99), fEMotherKF(-99), fEMotherErrKF(-99), fpxMotherKF(-99), fpxMotherErrKF(-99), fpyMotherKF(-99), fpyMotherErrKF(-99), fpzMotherKF(-99), fpzMotherErrKF(-99), fptMotherKF(-99), fptMotherErrKF(-99), fpMotherKF(-99), fpMotherErrKF(-99), fyMotherKF(-99), fctMotherKF(-99), fctMotherErrKF(-99),
   //
   fmSubMotherKF(-99), fmSubMotherErrKF(-99), fESubMotherKF(-99), fESubMotherErrKF(-99), fpxSubMotherKF(-99), fpxSubMotherErrKF(-99), fpySubMotherKF(-99), fpySubMotherErrKF(-99), fpzSubMotherKF(-99), fpzSubMotherErrKF(-99), fptSubMotherKF(-99), fptSubMotherErrKF(-99), fpSubMotherKF(-99), fpSubMotherErrKF(-99), fySubMotherKF(-99), fctSubMotherKF(-99), fctSubMotherErrKF(-99),
   //
   fEDaughterKF(-99), fLabelDaughterKF(-99), fptDaughterKF(-99), fpxDaughterKF(-99), fpyDaughterKF(-99), fpzDaughterKF(-99), fyDaughterKF(-99),
   fEDaughter1KF(-99), fLabelDaughter1KF(-99), fptDaughter1KF(-99), fpxDaughter1KF(-99), fpyDaughter1KF(-99), fpzDaughter1KF(-99), fyDaughter1KF(-99),
   fEDaughter2KF(-99), fLabelDaughter2KF(-99), fptDaughter2KF(-99), fpxDaughter2KF(-99), fpyDaughter2KF(-99), fpzDaughter2KF(-99), fyDaughter2KF(-99),
   fEDaughter3KF(-99), fLabelDaughter3KF(-99), fptDaughter3KF(-99), fpxDaughter3KF(-99), fpyDaughter3KF(-99), fpzDaughter3KF(-99), fyDaughter3KF(-99),
   //  
   fptDaughterUnProp(-99), fpxDaughterUnProp(-99), fpyDaughterUnProp(-99), fpzDaughterUnProp(-99), 
   fptDaughter1UnProp(-99), fpxDaughter1UnProp(-99), fpyDaughter1UnProp(-99), fpzDaughter1UnProp(-99), 
   fptDaughter2UnProp(-99), fpxDaughter2UnProp(-99), fpyDaughter2UnProp(-99), fpzDaughter2UnProp(-99), 
   fptDaughter3UnProp(-99), fpxDaughter3UnProp(-99), fpyDaughter3UnProp(-99), fpzDaughter3UnProp(-99), 
   //
   fDcaDaughterXYKF(-99), fDcaDaughter1XYKF(-99), fDcaDaughter2XYKF(-99), fDcaDaughter3XYKF(-99), fDcaDaughter4XYKF(-99), fDcaSecDaughterXYKF(-99), fDcaSecDaughter1XYKF(-99), fDcaSecDaughter2XYKF(-99), fDcaSecDaughter3XYKF(-99), fDcaSecDaughter4XYKF(-99), fDcaDaughterZKF(-99), fDcaDaughter1ZKF(-99), fDcaDaughter2ZKF(-99), fDcaDaughter3ZKF(-99), fDcaDaughter4ZKF(-99), fDcaSecDaughterZKF(-99), fDcaSecDaughter1ZKF(-99), fDcaSecDaughter2ZKF(-99), fDcaSecDaughter3ZKF(-99), fDcaSecDaughter4ZKF(-99),
   //
   fDCA2BXYKF(-99), fDCA3B1XYKF(-99), fDCA3B2XYKF(-99), fDCA3B3XYKF(-99), fDCA2BZKF(-99), fDCA3B1ZKF(-99), fDCA3B2ZKF(-99), fDCA3B3ZKF(-99)
{
  DefineInput(0,   TChain::Class());
  DefineOutput(1,  TList::Class());
  DefineOutput(2,  TTree::Class());
  DefineOutput(3,  TTree::Class());
  DefineOutput(4,  TTree::Class());
  DefineOutput(5,  TTree::Class());
  DefineOutput(6,  TTree::Class());
  DefineOutput(7,  TTree::Class());
  DefineOutput(8,  TTree::Class());
  DefineOutput(9,  TTree::Class());
  DefineOutput(10, TTree::Class());
  DefineOutput(11, TTree::Class());

}
// _________________________________________________ //
// __ Destructor __ //
AliAnalysisTaskDoubleHypNucTreeLS::~AliAnalysisTaskDoubleHypNucTreeLS() {

  ResetVals("Event");

}
// _________________________________________________ //
const Int_t AliAnalysisTaskDoubleHypNucTreeLS::fgkPdgCode[] = {
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
								 1010010030,          //HyperHydrogen 3
								 -1010010030,          //Anti-HyperHydrogen3
								 1010010040,          //HyperHydrogen 4
								 -1010010040,          //Anti-HyperHydrogen 4
								 1010020040,           //HyperHelium 4
								 -1010020040,           //Anti-HyperHelium 4
								 1010010041,          //HyperHydrogen 4*
								 -1010010041,          //Anti-HyperHydrogen 4*
								 1010020041,           //HyperHelium 4*
								 -1010020041,           //Anti-HyperHelium 4*
								 1010020050,          //HyperHelium 5
								 -1010020050,          //Anti-HyperHelium 5
								 1020010040,          //DoubleHyperHydrogen 4
								 -1020010040          //Anti-DoubleHyperHydrogen 4
};
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTreeLS::UserCreateOutputObjects() {
  fInputHandler = dynamic_cast<AliInputEventHandler*>
    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!fInputHandler) {
    AliError("Could not get ESD InputHandler.\n");
    return;
  }
  fPID = fInputHandler->GetPIDResponse();
  if (!fPID) {
    AliError("Could not get PID response.\n");
    return;
  }
  // _________________________________________________ //
  // __ Histograms __ //
  fHistdEdx = new TH2F("fHistdEdX", "dE/dx;#frac{#it{p}}{z} (GeV/#it{c});TPC Signal (a.u.)", 1000, -5.0, 5.0, 1000, 0.0, 2000);

  fHistNumEvents = new TH1F("fHistNumEvents", "Number of Events", 4, 0, 4);
  fHistNumEvents->GetXaxis()->SetBinLabel(1, "before PhysSel");
  fHistNumEvents->GetXaxis()->SetBinLabel(2, "after PhysSel");
  fHistNumEvents->GetXaxis()->SetBinLabel(3, "after PhysSel w TrackVert");
  fHistNumEvents->GetXaxis()->SetBinLabel(4, "after PhysSel w/o TrackVert");

  fHistTrigger = new TH1F("fHistTrigger", "Trigger", 8, 0, 8);
  fHistTrigger->GetXaxis()->SetBinLabel(1, "other");
  fHistTrigger->GetXaxis()->SetBinLabel(2, "kINT7");
  fHistTrigger->GetXaxis()->SetBinLabel(3, "kHighMultV0");
  fHistTrigger->GetXaxis()->SetBinLabel(4, "kHighMultSPD");
  fHistTrigger->GetXaxis()->SetBinLabel(5, "HNU");
  fHistTrigger->GetXaxis()->SetBinLabel(6, "HQU");
  fHistTrigger->GetXaxis()->SetBinLabel(7, "kCentral");
  fHistTrigger->GetXaxis()->SetBinLabel(8, "kSemiCentral");

  fHistTrigger1 = new TH1F("fHistTrigger1", "Trigger", 8, 0, 8);
  fHistTrigger1->GetXaxis()->SetBinLabel(1, "other");
  fHistTrigger1->GetXaxis()->SetBinLabel(2, "kINT7");
  fHistTrigger1->GetXaxis()->SetBinLabel(3, "kHighMultV0");
  fHistTrigger1->GetXaxis()->SetBinLabel(4, "kHighMultSPD");
  fHistTrigger1->GetXaxis()->SetBinLabel(5, "HNU");
  fHistTrigger1->GetXaxis()->SetBinLabel(6, "HQU");
  fHistTrigger1->GetXaxis()->SetBinLabel(7, "kCentral");
  fHistTrigger1->GetXaxis()->SetBinLabel(8, "kSemiCentral");

  fHistTrigger2 = new TH1F("fHistTrigger2", "Trigger", 8, 0, 8);
  fHistTrigger2->GetXaxis()->SetBinLabel(1, "other");
  fHistTrigger2->GetXaxis()->SetBinLabel(2, "kINT7");
  fHistTrigger2->GetXaxis()->SetBinLabel(3, "kHighMultV0");
  fHistTrigger2->GetXaxis()->SetBinLabel(4, "kHighMultSPD");
  fHistTrigger2->GetXaxis()->SetBinLabel(5, "HNU");
  fHistTrigger2->GetXaxis()->SetBinLabel(6, "HQU");
  fHistTrigger2->GetXaxis()->SetBinLabel(7, "kCentral");
  fHistTrigger2->GetXaxis()->SetBinLabel(8, "kSemiCentral");

  fHistCentrality1 = new TH1F("fHistCentrality1", "", 100, 0, 100);
  fHistCentrality2 = new TH1F("fHistCentrality2", "", 100, 0, 100);

  fHistNtracks1 = new TH1F("fHistNtracks1", "", 20000, 0, 20000);
  fHistNtracks2 = new TH1F("fHistNtracks2", "", 20000, 0, 20000);

  fHistogramList = new TList();
  fHistogramList->SetOwner(kTRUE);
  fHistogramList->SetName(GetName());
  fHistogramList->Add(fHistdEdx);
  fHistogramList->Add(fHistNumEvents);
  fHistogramList->Add(fHistTrigger);
  fHistogramList->Add(fHistTrigger1);
  fHistogramList->Add(fHistTrigger2);
  fHistogramList->Add(fHistCentrality1);
  fHistogramList->Add(fHistCentrality2);
  fHistogramList->Add(fHistNtracks1);
  fHistogramList->Add(fHistNtracks2);
  fEventCuts.AddQAplotsToList(fHistogramList);
  // _________________________________________________ //
  // __ associated 4LLH Tree __ //
  fTree = new TTree("fTree", "fTree");
  fTree->Branch("fPeriod", &fPeriod, "fPeriod/I");
  fTree->Branch("frunnumber", &frunnumber, "frunnumber/I");
  fTree->Branch("fEventID", &fEventID, "fEventID/I");
  fTree->Branch("feventclass", &feventclass, "feventclass/I");
  fTree->Branch("fCentrality", &fCentrality, "fCentrality/I");
  fTree->Branch("fMagneticField", &fMagneticField, "fMagneticField/I");
  fTree->Branch("fisWOTrackVertex", &fisWOTrackVertex, "fisWOTrackVertex/I");
  fTree->Branch("fmctruth", &fmctruth, "fmctruth/I");
  fTree->Branch("fTrigMB", &fTrigMB, "fTrigMB/I");
  fTree->Branch("fTrigHMV0", &fTrigHMV0, "fTrigHMV0/I");
  fTree->Branch("fTrigHMSPD", &fTrigHMSPD, "fTrigHMSPD/I");
  fTree->Branch("fTrigHNU", &fTrigHNU, "fTrigHNU/I");
  fTree->Branch("fTrigHQU", &fTrigHQU, "fTrigHQU/I");
  fTree->Branch("fTrigkCentral", &fTrigkCentral, "fTrigkCentral/I");
  fTree->Branch("fTrigkSemiCentral", &fTrigkSemiCentral, "fTrigkSemiCentral/I");

  fTree->Branch("fParticleID", &fParticleID, "fParticleID/I");     
  fTree->Branch("fPDGMother", &fPDGMother, "fPDGMother/I");
  fTree->Branch("fDecayChannel", &fDecayChannel, "fDecayChannel/I");
  fTree->Branch("fRecoMethod", &fRecoMethod, "fRecoMethod/I");
  fTree->Branch("fChargeMother", &fChargeMother, "fChargeMother/I");
  fTree->Branch("fPrimary4LHe", &fPrimary4LHe, "fPrimary4LHe/I");

  fTree->Branch("fPrimVertexX", &fPrimVertexX, "fPrimVertexX/D");
  fTree->Branch("fPrimVertexY", &fPrimVertexY, "fPrimVertexY/D");
  fTree->Branch("fPrimVertexZ", &fPrimVertexZ, "fPrimVertexZ/D");
  fTree->Branch("fSecVertexX", &fSecVertexX, "fSecVertexX/D");
  fTree->Branch("fSecVertexY", &fSecVertexY, "fSecVertexY/D");
  fTree->Branch("fSecVertexZ", &fSecVertexZ, "fSecVertexZ/D");
  fTree->Branch("fTertVertexX", &fTertVertexX, "fTertVertexX/D");
  fTree->Branch("fTertVertexY", &fTertVertexY, "fTertVertexY/D");
  fTree->Branch("fTertVertexZ", &fTertVertexZ, "fTertVertexZ/D");
  fTree->Branch("fPrimVertChi2", &fPrimVertChi2, "fPrimVertChi2/D");
  fTree->Branch("fSecVertChi2", &fSecVertChi2, "fSecVertChi2/D");
  fTree->Branch("fTertVertChi2", &fTertVertChi2, "fTertVertChi2/D");
  fTree->Branch("fPrimVertNDF", &fPrimVertNDF, "fPrimVertNDF/I");
  fTree->Branch("fSecVertNDF", &fSecVertNDF, "fSecVertNDF/I");
  fTree->Branch("fTertVertNDF", &fTertVertNDF, "fTertVertNDF/I");
  /*fTree->Branch("fisOnlineV0_13", &fisOnlineV0_13, "fisOnlineV0_13/I");
    fTree->Branch("fV0VertexX_13", &fV0VertexX_13, "fV0VertexX_13/D");
    fTree->Branch("fV0VertexY_13", &fV0VertexY_13, "fV0VertexY_13/D");
    fTree->Branch("fV0VertexZ_13", &fV0VertexZ_13, "fV0VertexZ_13/D");
    fTree->Branch("fisOnlineV0_14", &fisOnlineV0_14, "fisOnlineV0_14/I");
    fTree->Branch("fV0VertexX_14", &fV0VertexX_14, "fV0VertexX_14/D");
    fTree->Branch("fV0VertexY_14", &fV0VertexY_14, "fV0VertexY_14/D");
    fTree->Branch("fV0VertexZ_14", &fV0VertexZ_14, "fV0VertexZ_14/D");
    fTree->Branch("fisOnlineV0_23", &fisOnlineV0_23, "fisOnlineV0_23/I");
    fTree->Branch("fV0VertexX_23", &fV0VertexX_23, "fV0VertexX_23/D");
    fTree->Branch("fV0VertexY_23", &fV0VertexY_23, "fV0VertexY_23/D");
    fTree->Branch("fV0VertexZ_23", &fV0VertexZ_23, "fV0VertexZ_23/D");
    fTree->Branch("fisOnlineV0_24", &fisOnlineV0_24, "fisOnlineV0_24/I");
    fTree->Branch("fV0VertexX_24", &fV0VertexX_24, "fV0VertexX_24/D");
    fTree->Branch("fV0VertexY_24", &fV0VertexY_24, "fV0VertexY_24/D");
    fTree->Branch("fV0VertexZ_24", &fV0VertexZ_24, "fV0VertexZ_24/D");*/
  fTree->Branch("fDCA2B", &fDCA2B, "fDCA2B/D");
  fTree->Branch("fDCA2Bo", &fDCA2Bo, "fDCA2Bo/D");
  fTree->Branch("fDCA3B1", &fDCA3B1, "fDCA3B1/D");
  fTree->Branch("fDCA3B2", &fDCA3B2, "fDCA3B2/D");
  fTree->Branch("fDCA3B3", &fDCA3B3, "fDCA3B3/D");
  fTree->Branch("fPA", &fPA, "fPA/D");
  fTree->Branch("fSubPA", &fSubPA, "fSubPA/D");
  fTree->Branch("fSubPA2", &fSubPA2, "fSubPA2/D");
  fTree->Branch("fDecAngle", &fDecAngle, "fDecAngle/D");
  fTree->Branch("farmalpha", &farmalpha, "farmalpha/D");
  fTree->Branch("farmpt", &farmpt, "farmpt/D");
   
  fTree->Branch("fmMother", &fmMother, "fmMother/D");
  fTree->Branch("fmMother2", &fmMother2, "fmMother2/D");
  fTree->Branch("fEMother", &fEMother, "fEMother/D");
  fTree->Branch("fpxMother", &fpxMother, "fpxMother/D");
  fTree->Branch("fpyMother", &fpyMother, "fpyMother/D");
  fTree->Branch("fpzMother", &fpzMother, "fpzMother/D");
  fTree->Branch("fptMother", &fptMother, "fptMother/D");
  fTree->Branch("fpMother", &fpMother, "fpMother/D");
  fTree->Branch("fyMother", &fyMother, "fyMother/D");
  fTree->Branch("fctMother", &fctMother, "fctMother/D");
  fTree->Branch("fmSubMother", &fmSubMother, "fmSubMother/D");
  fTree->Branch("fESubMother", &fESubMother, "fESubMother/D");
  fTree->Branch("fpxSubMother", &fpxSubMother, "fpxSubMother/D");
  fTree->Branch("fpySubMother", &fpySubMother, "fpySubMother/D");
  fTree->Branch("fpzSubMother", &fpzSubMother, "fpzSubMother/D");
  fTree->Branch("fptSubMother", &fptSubMother, "fptSubMother/D");
  fTree->Branch("fpSubMother", &fpSubMother, "fpSubMother/D");
  fTree->Branch("fySubMother", &fySubMother, "fySubMother/D");
  fTree->Branch("fctSubMother", &fctSubMother, "fctSubMother/D");

  fTree->Branch("fEDaughter", &fEDaughter, "fEDaughter/D");
  fTree->Branch("fptDaughter", &fptDaughter, "fptDaughter/D");
  fTree->Branch("fpxDaughter", &fpxDaughter, "fpxDaughter/D");
  fTree->Branch("fpyDaughter", &fpyDaughter, "fpyDaughter/D");
  fTree->Branch("fpzDaughter", &fpzDaughter, "fpzDaughter/D");
  fTree->Branch("fyDaughter", &fyDaughter, "fyDaughter/D");
  fTree->Branch("fEDaughter1", &fEDaughter1, "fEDaughter1/D");
  fTree->Branch("fptDaughter1", &fptDaughter1, "fptDaughter1/D");
  fTree->Branch("fpxDaughter1", &fpxDaughter1, "fpxDaughter1/D");
  fTree->Branch("fpyDaughter1", &fpyDaughter1, "fpyDaughter1/D");
  fTree->Branch("fpzDaughter1", &fpzDaughter1, "fpzDaughter1/D");
  fTree->Branch("fyDaughter1", &fyDaughter1, "fyDaughter1/D");
  fTree->Branch("fEDaughter2", &fEDaughter2, "fEDaughter2/D");
  fTree->Branch("fptDaughter2", &fptDaughter2, "fptDaughter2/D");
  fTree->Branch("fpxDaughter2", &fpxDaughter2, "fpxDaughter2/D");
  fTree->Branch("fpyDaughter2", &fpyDaughter2, "fpyDaughter2/D");
  fTree->Branch("fpzDaughter2", &fpzDaughter2, "fpzDaughter2/D");
  fTree->Branch("fyDaughter2", &fyDaughter2, "fyDaughter2/D");
  fTree->Branch("fEDaughter3", &fEDaughter3, "fEDaughter3/D");
  fTree->Branch("fptDaughter3", &fptDaughter3, "fptDaughter3/D");
  fTree->Branch("fpxDaughter3", &fpxDaughter3, "fpxDaughter3/D");
  fTree->Branch("fpyDaughter3", &fpyDaughter3, "fpyDaughter3/D");
  fTree->Branch("fpzDaughter3", &fpzDaughter3, "fpzDaughter3/D");
  fTree->Branch("fyDaughter3", &fyDaughter3, "fyDaughter3/D");

  fTree->Branch("fPropDCADaughter", &fPropDCADaughter, "fPropDCADaughter/I");
  fTree->Branch("fImParDaughter", &fImParDaughter, "fImParDaughter/D");
  fTree->Branch("fImParzDaughter", &fImParzDaughter, "fImParzDaughter/D");
  fTree->Branch("fPropDCADaughter1", &fPropDCADaughter1, "fPropDCADaughter1/I");
  fTree->Branch("fImParDaughter1", &fImParDaughter1, "fImParDaughter1/D");
  fTree->Branch("fImParzDaughter1", &fImParzDaughter1, "fImParzDaughter1/D");
  fTree->Branch("fPropDCADaughter2", &fPropDCADaughter2, "fPropDCADaughter2/I");
  fTree->Branch("fImParDaughter2", &fImParDaughter2, "fImParDaughter2/D");
  fTree->Branch("fImParzDaughter2", &fImParzDaughter2, "fImParzDaughter2/D");
  fTree->Branch("fPropDCADaughter3", &fPropDCADaughter3, "fPropDCADaughter3/I");
  fTree->Branch("fImParDaughter3", &fImParDaughter3, "fImParDaughter3/D");
  fTree->Branch("fImParzDaughter3", &fImParzDaughter3, "fImParzDaughter3/D");
  fTree->Branch("fPropDCADaughter4", &fPropDCADaughter4, "fPropDCADaughter4/I");
  fTree->Branch("fImParDaughter4", &fImParDaughter4, "fImParDaughter4/D");
  fTree->Branch("fImParzDaughter4", &fImParzDaughter4, "fImParzDaughter4/D");
  fTree->Branch("fDcaSecDaughter", &fDcaSecDaughter, "fDcaSecDaughter/D");
  fTree->Branch("fDcaSecDaughter1", &fDcaSecDaughter1, "fDcaSecDaughter1/D");
  fTree->Branch("fDcaSecDaughter2", &fDcaSecDaughter2, "fDcaSecDaughter2/D");
  fTree->Branch("fDcaSecDaughter3", &fDcaSecDaughter3, "fDcaSecDaughter3/D");
  fTree->Branch("fDcaSecDaughter4", &fDcaSecDaughter4, "fDcaSecDaughter4/D");

  fTree->Branch("fTrackPIDDaughter", &fTrackPIDDaughter, "fTrackPIDDaughter/I");
  fTree->Branch("fLabelDaughterKF", &fLabelDaughterKF, "fLabelDaughterKF/D");
  fTree->Branch("fpDaughter", &fpDaughter, "fpDaughter/D");
  fTree->Branch("fptDaughterUnProp", &fptDaughterUnProp, "fptDaughterUnProp/D");
  fTree->Branch("fpxDaughterUnProp", &fpxDaughterUnProp, "fpxDaughterUnProp/D");
  fTree->Branch("fpyDaughterUnProp", &fpyDaughterUnProp, "fpyDaughterUnProp/D");
  fTree->Branch("fpzDaughterUnProp", &fpzDaughterUnProp, "fpzDaughterUnProp/D");
  fTree->Branch("fdEdxDaughter", &fdEdxDaughter, "fdEdxDaughter/D");
  fTree->Branch("fdEdxSigmaDaughter", &fdEdxSigmaDaughter, "fdEdxSigmaDaughter/D");
  fTree->Branch("fDcaDaughter", &fDcaDaughter, "fDcaDaughter/D");
  fTree->Branch("fDcaDaughtero", &fDcaDaughtero, "fDcaDaughtero/D");
  fTree->Branch("fDcazDaughter", &fDcazDaughter, "fDcazDaughter/D");
  fTree->Branch("fNclsDaughter", &fNclsDaughter, "fNclsDaughter/I");
  fTree->Branch("fChi2Daughter", &fChi2Daughter, "fChi2Daughter/D");
  fTree->Branch("fNclsITSDaughter", &fNclsITSDaughter, "fNclsITSDaughter/I");
  fTree->Branch("fEtaDaughter", &fEtaDaughter, "fEtaDaughter/D");
  fTree->Branch("fPhiDaughter", &fPhiDaughter, "fPhiDaughter/D");
  fTree->Branch("fGeoLengthDaughter", &fGeoLengthDaughter, "fGeoLengthDaughter/D");
  fTree->Branch("fTOFSignalDaughter", &fTOFSignalDaughter, "fTOFSignalDaughter/D");
  fTree->Branch("fSigmaYXDaughter", &fSigmaYXDaughter, "fSigmaYXDaughter/D");
  fTree->Branch("fSigmaXYZDaughter", &fSigmaXYZDaughter, "fSigmaXYZDaughter/D");
  fTree->Branch("fSigmaZDaughter", &fSigmaZDaughter, "fSigmaZDaughter/D");
  fTree->Branch("fPtUncertDaughter", &fPtUncertDaughter, "fPtUncertDaughter/D");
  fTree->Branch("fTPCRefitDaughter", &fTPCRefitDaughter, "fTPCRefitDaughter/I");
  fTree->Branch("fITSRefitDaughter", &fITSRefitDaughter, "fITSRefitDaughter/I");
  fTree->Branch("fITSLayer1Daughter", &fITSLayer1Daughter, "fITSLayer1Daughter/I");
  fTree->Branch("fITSLayer2Daughter", &fITSLayer2Daughter, "fITSLayer2Daughter/I");
  fTree->Branch("fITSLayer3Daughter", &fITSLayer3Daughter, "fITSLayer3Daughter/I");
  fTree->Branch("fITSLayer4Daughter", &fITSLayer4Daughter, "fITSLayer4Daughter/I");
  fTree->Branch("fITSLayer5Daughter", &fITSLayer5Daughter, "fITSLayer5Daughter/I");
  fTree->Branch("fITSLayer6Daughter", &fITSLayer6Daughter, "fITSLayer6Daughter/I");

  fTree->Branch("fTrackPIDDaughter1", &fTrackPIDDaughter1, "fTrackPIDDaughter1/I");
  fTree->Branch("fLabelDaughter1KF", &fLabelDaughter1KF, "fLabelDaughter1KF/D");
  fTree->Branch("fpDaughter1", &fpDaughter1, "fpDaughter1/D");
  fTree->Branch("fptDaughter1UnProp", &fptDaughter1UnProp, "fptDaughter1UnProp/D");
  fTree->Branch("fpxDaughter1UnProp", &fpxDaughter1UnProp, "fpxDaughter1UnProp/D");
  fTree->Branch("fpyDaughter1UnProp", &fpyDaughter1UnProp, "fpyDaughter1UnProp/D");
  fTree->Branch("fpzDaughter1UnProp", &fpzDaughter1UnProp, "fpzDaughter1UnProp/D");
  fTree->Branch("fdEdxDaughter1", &fdEdxDaughter1, "fdEdxDaughter1/D");
  fTree->Branch("fdEdxSigmaDaughter1", &fdEdxSigmaDaughter1, "fdEdxSigmaDaughter1/D");
  fTree->Branch("fDcaDaughter1", &fDcaDaughter1, "fDcaDaughter1/D");
  fTree->Branch("fDcaDaughter1o", &fDcaDaughter1o, "fDcaDaughter1o/D");
  fTree->Branch("fDcazDaughter1", &fDcazDaughter1, "fDcazDaughter1/D");
  fTree->Branch("fNclsDaughter1", &fNclsDaughter1, "fNclsDaughter1/I");
  fTree->Branch("fChi2Daughter1", &fChi2Daughter1, "fChi2Daughter1/D");
  fTree->Branch("fNclsITSDaughter1", &fNclsITSDaughter1, "fNclsITSDaughter1/I");
  fTree->Branch("fEtaDaughter1", &fEtaDaughter1, "fEtaDaughter1/D");
  fTree->Branch("fPhiDaughter1", &fPhiDaughter1, "fPhiDaughter1/D");
  fTree->Branch("fGeoLengthDaughter1", &fGeoLengthDaughter1, "fGeoLengthDaughter1/D");
  fTree->Branch("fTOFSignalDaughter1", &fTOFSignalDaughter1, "fTOFSignalDaughter1/D");
  fTree->Branch("fSigmaYXDaughter1", &fSigmaYXDaughter1, "fSigmaYXDaughter1/D");
  fTree->Branch("fSigmaXYZDaughter1", &fSigmaXYZDaughter1, "fSigmaXYZDaughter1/D");
  fTree->Branch("fSigmaZDaughter1", &fSigmaZDaughter1, "fSigmaZDaughter1/D");
  fTree->Branch("fPtUncertDaughter1", &fPtUncertDaughter1, "fPtUncertDaughter1/D");
  fTree->Branch("fTPCRefitDaughter1", &fTPCRefitDaughter1, "fTPCRefitDaughter1/I");
  fTree->Branch("fITSRefitDaughter1", &fITSRefitDaughter1, "fITSRefitDaughter1/I");
  fTree->Branch("fITSLayer1Daughter1", &fITSLayer1Daughter1, "fITSLayer1Daughter1/I");
  fTree->Branch("fITSLayer2Daughter1", &fITSLayer2Daughter1, "fITSLayer2Daughter1/I");
  fTree->Branch("fITSLayer3Daughter1", &fITSLayer3Daughter1, "fITSLayer3Daughter1/I");
  fTree->Branch("fITSLayer4Daughter1", &fITSLayer4Daughter1, "fITSLayer4Daughter1/I");
  fTree->Branch("fITSLayer5Daughter1", &fITSLayer5Daughter1, "fITSLayer5Daughter1/I");
  fTree->Branch("fITSLayer6Daughter1", &fITSLayer6Daughter1, "fITSLayer6Daughter1/I");

  fTree->Branch("fTrackPIDDaughter2", &fTrackPIDDaughter2, "fTrackPIDDaughter2/I");
  fTree->Branch("fLabelDaughter2KF", &fLabelDaughter2KF, "fLabelDaughter2KF/D");
  fTree->Branch("fpDaughter2", &fpDaughter2, "fpDaughter2/D");
  fTree->Branch("fptDaughter2UnProp", &fptDaughter2UnProp, "fptDaughter2UnProp/D");
  fTree->Branch("fpxDaughter2UnProp", &fpxDaughter2UnProp, "fpxDaughter2UnProp/D");
  fTree->Branch("fpyDaughter2UnProp", &fpyDaughter2UnProp, "fpyDaughter2UnProp/D");
  fTree->Branch("fpzDaughter2UnProp", &fpzDaughter2UnProp, "fpzDaughter2UnProp/D");
  fTree->Branch("fdEdxDaughter2", &fdEdxDaughter2, "fdEdxDaughter2/D");
  fTree->Branch("fdEdxSigmaDaughter2", &fdEdxSigmaDaughter2, "fdEdxSigmaDaughter2/D");
  fTree->Branch("fDcaDaughter2", &fDcaDaughter2, "fDcaDaughter2/D");
  fTree->Branch("fDcaDaughter2o", &fDcaDaughter2o, "fDcaDaughter2o/D");
  fTree->Branch("fDcazDaughter2", &fDcazDaughter2, "fDcazDaughter2/D");
  fTree->Branch("fNclsDaughter2", &fNclsDaughter2, "fNclsDaughter2/I");
  fTree->Branch("fChi2Daughter2", &fChi2Daughter2, "fChi2Daughter2/D");
  fTree->Branch("fNclsITSDaughter2", &fNclsITSDaughter2, "fNclsITSDaughter2/I");
  fTree->Branch("fEtaDaughter2", &fEtaDaughter2, "fEtaDaughter2/D");
  fTree->Branch("fPhiDaughter2", &fPhiDaughter2, "fPhiDaughter2/D");
  fTree->Branch("fGeoLengthDaughter2", &fGeoLengthDaughter2, "fGeoLengthDaughter2/D");
  fTree->Branch("fTOFSignalDaughter2", &fTOFSignalDaughter2, "fTOFSignalDaughter2/D");
  fTree->Branch("fSigmaYXDaughter2", &fSigmaYXDaughter2, "fSigmaYXDaughter2/D");
  fTree->Branch("fSigmaXYZDaughter2", &fSigmaXYZDaughter2, "fSigmaXYZDaughter2/D");
  fTree->Branch("fSigmaZDaughter2", &fSigmaZDaughter2, "fSigmaZDaughter2/D");
  fTree->Branch("fPtUncertDaughter2", &fPtUncertDaughter2, "fPtUncertDaughter2/D");
  fTree->Branch("fTPCRefitDaughter2", &fTPCRefitDaughter2, "fTPCRefitDaughter2/I");
  fTree->Branch("fITSRefitDaughter2", &fITSRefitDaughter2, "fITSRefitDaughter2/I");
  fTree->Branch("fITSLayer1Daughter2", &fITSLayer1Daughter2, "fITSLayer1Daughter2/I");
  fTree->Branch("fITSLayer2Daughter2", &fITSLayer2Daughter2, "fITSLayer2Daughter2/I");
  fTree->Branch("fITSLayer3Daughter2", &fITSLayer3Daughter2, "fITSLayer3Daughter2/I");
  fTree->Branch("fITSLayer4Daughter2", &fITSLayer4Daughter2, "fITSLayer4Daughter2/I");
  fTree->Branch("fITSLayer5Daughter2", &fITSLayer5Daughter2, "fITSLayer5Daughter2/I");
  fTree->Branch("fITSLayer6Daughter2", &fITSLayer6Daughter2, "fITSLayer6Daughter2/I");

  fTree->Branch("fTrackPIDDaughter3", &fTrackPIDDaughter3, "fTrackPIDDaughter3/I");
  fTree->Branch("fLabelDaughter3KF", &fLabelDaughter3KF, "fLabelDaughter3KF/D");
  fTree->Branch("fpDaughter3", &fpDaughter3, "fpDaughter3/D");
  fTree->Branch("fptDaughter3UnProp", &fptDaughter3UnProp, "fptDaughter3UnProp/D");
  fTree->Branch("fpxDaughter3UnProp", &fpxDaughter3UnProp, "fpxDaughter3UnProp/D");
  fTree->Branch("fpyDaughter3UnProp", &fpyDaughter3UnProp, "fpyDaughter3UnProp/D");
  fTree->Branch("fpzDaughter3UnProp", &fpzDaughter3UnProp, "fpzDaughter3UnProp/D");
  fTree->Branch("fdEdxDaughter3", &fdEdxDaughter3, "fdEdxDaughter3/D");
  fTree->Branch("fdEdxSigmaDaughter3", &fdEdxSigmaDaughter3, "fdEdxSigmaDaughter3/D");
  fTree->Branch("fDcaDaughter3", &fDcaDaughter3, "fDcaDaughter3/D");
  fTree->Branch("fDcaDaughter3o", &fDcaDaughter3o, "fDcaDaughter3o/D");
  fTree->Branch("fDcazDaughter3", &fDcazDaughter3, "fDcazDaughter3/D");
  fTree->Branch("fNclsDaughter3", &fNclsDaughter3, "fNclsDaughter3/I");
  fTree->Branch("fChi2Daughter3", &fChi2Daughter3, "fChi2Daughter3/D");
  fTree->Branch("fNclsITSDaughter3", &fNclsITSDaughter3, "fNclsITSDaughter3/I");
  fTree->Branch("fEtaDaughter3", &fEtaDaughter3, "fEtaDaughter3/D");
  fTree->Branch("fPhiDaughter3", &fPhiDaughter3, "fPhiDaughter3/D");
  fTree->Branch("fGeoLengthDaughter3", &fGeoLengthDaughter3, "fGeoLengthDaughter3/D");
  fTree->Branch("fTOFSignalDaughter3", &fTOFSignalDaughter3, "fTOFSignalDaughter3/D");
  fTree->Branch("fSigmaYXDaughter3", &fSigmaYXDaughter3, "fSigmaYXDaughter3/D");
  fTree->Branch("fSigmaXYZDaughter3", &fSigmaXYZDaughter3, "fSigmaXYZDaughter3/D");
  fTree->Branch("fSigmaZDaughter3", &fSigmaZDaughter3, "fSigmaZDaughter3/D");
  fTree->Branch("fPtUncertDaughter3", &fPtUncertDaughter3, "fPtUncertDaughter3/D");
  fTree->Branch("fTPCRefitDaughter3", &fTPCRefitDaughter3, "fTPCRefitDaughter3/I");
  fTree->Branch("fITSRefitDaughter3", &fITSRefitDaughter3, "fITSRefitDaughter3/I");
  fTree->Branch("fITSLayer1Daughter3", &fITSLayer1Daughter3, "fITSLayer1Daughter3/I");
  fTree->Branch("fITSLayer2Daughter3", &fITSLayer2Daughter3, "fITSLayer2Daughter3/I");
  fTree->Branch("fITSLayer3Daughter3", &fITSLayer3Daughter3, "fITSLayer3Daughter3/I");
  fTree->Branch("fITSLayer4Daughter3", &fITSLayer4Daughter3, "fITSLayer4Daughter3/I");
  fTree->Branch("fITSLayer5Daughter3", &fITSLayer5Daughter3, "fITSLayer5Daughter3/I");
  fTree->Branch("fITSLayer6Daughter3", &fITSLayer6Daughter3, "fITSLayer6Daughter3/I");

  fTree->Branch("fEDaughter4", &fEDaughter4, "fEDaughter4/D");
  fTree->Branch("fpDaughter4", &fpDaughter4, "fpDaughter4/D");
  fTree->Branch("fptDaughter4", &fptDaughter4, "fptDaughter4/D");
  fTree->Branch("fpxDaughter4", &fpxDaughter4, "fpxDaughter4/D");
  fTree->Branch("fpyDaughter4", &fpyDaughter4, "fpyDaughter4/D");
  fTree->Branch("fpzDaughter4", &fpzDaughter4, "fpzDaughter4/D");
  fTree->Branch("fyDaughter4", &fyDaughter4, "fyDaughter4/D");
  fTree->Branch("fDcaDaughter4", &fDcaDaughter4, "fDcaDaughter4/D");
  fTree->Branch("fDcazDaughter4", &fDcazDaughter4, "fDcazDaughter4/D");
  fTree->Branch("fSigmaYXDaughter4", &fSigmaYXDaughter4, "fSigmaYXDaughter4/D");
  fTree->Branch("fSigmaXYZDaughter4", &fSigmaXYZDaughter4, "fSigmaXYZDaughter4/D");
  fTree->Branch("fSigmaZDaughter4", &fSigmaZDaughter4, "fSigmaZDaughter4/D");
  fTree->Branch("fPtUncertDaughter4", &fPtUncertDaughter4, "fPtUncertDaughter4/D");

  fTree->Branch("fdEdxSigmaPion", &fdEdxSigmaPion, "fdEdxSigmaPion/D");
  fTree->Branch("fdEdxSigmaDeuteron", &fdEdxSigmaDeuteron, "fdEdxSigmaDeuteron/D");
  fTree->Branch("fdEdxSigmaTriton", &fdEdxSigmaTriton, "fdEdxSigmaTriton/D");
  fTree->Branch("fdEdxSigmaAlpha", &fdEdxSigmaAlpha, "fdEdxSigmaAlpha/D");

  fTree->Branch("fCovMatrixTrack", fCovMatrixTrack, "fCovMatrixTrack[21]/D");
  fTree->Branch("fCovMatrixTrack1", fCovMatrixTrack1, "fCovMatrixTrack1[21]/D");
  fTree->Branch("fCovMatrixTrack2", fCovMatrixTrack2, "fCovMatrixTrack2[21]/D");
  fTree->Branch("fCovMatrixTrack3", fCovMatrixTrack3, "fCovMatrixTrack3[21]/D");
  fTree->Branch("fCovMatrixTrack4", fCovMatrixTrack4, "fCovMatrixTrack4[21]/D");
  fTree->Branch("fTrackPar", fTrackPar, "fTrackPar[7]/D");
  fTree->Branch("fTrackPar1", fTrackPar1, "fTrackPar1[7]/D");
  fTree->Branch("fTrackPar2", fTrackPar2, "fTrackPar2[7]/D");
  fTree->Branch("fTrackPar3", fTrackPar3, "fTrackPar3[7]/D");
  fTree->Branch("fTrackPar4", fTrackPar4, "fTrackPar4[7]/D");

  // _________________________________________________ //
  // __ associated 4LLH KF tree __ //
  fTreeKF = new TTree("fTreeKF", "fTreeKF");
  fTreeKF->Branch("fPeriod", &fPeriod, "fPeriod/I");
  fTreeKF->Branch("frunnumber", &frunnumber, "frunnumber/I");
  fTreeKF->Branch("fEventID", &fEventID, "fEventID/I");
  fTreeKF->Branch("feventclass", &feventclass, "feventclass/I");
  fTreeKF->Branch("fCentrality", &fCentrality, "fCentrality/I");
  fTreeKF->Branch("fMagneticField", &fMagneticField, "fMagneticField/I");
  fTreeKF->Branch("fisWOTrackVertex", &fisWOTrackVertex, "fisWOTrackVertex/I");
  fTreeKF->Branch("fmctruth", &fmctruth, "fmctruth/I");
  fTreeKF->Branch("fTrigMB", &fTrigMB, "fTrigMB/I");
  fTreeKF->Branch("fTrigHMV0", &fTrigHMV0, "fTrigHMV0/I");
  fTreeKF->Branch("fTrigHMSPD", &fTrigHMSPD, "fTrigHMSPD/I");
  fTreeKF->Branch("fTrigHNU", &fTrigHNU, "fTrigHNU/I");
  fTreeKF->Branch("fTrigHQU", &fTrigHQU, "fTrigHQU/I");
  fTreeKF->Branch("fTrigkCentral", &fTrigkCentral, "fTrigkCentral/I");
  fTreeKF->Branch("fTrigkSemiCentral", &fTrigkSemiCentral, "fTrigkSemiCentral/I");

  fTreeKF->Branch("fParticleID", &fParticleID, "fParticleID/I");     
  fTreeKF->Branch("fPDGMother", &fPDGMother, "fPDGMother/I");
  fTreeKF->Branch("fDecayChannel", &fDecayChannel, "fDecayChannel/I");
  fTreeKF->Branch("fRecoMethod", &fRecoMethod, "fRecoMethod/I");
  fTreeKF->Branch("fChargeMother", &fChargeMother, "fChargeMother/I");
  fTreeKF->Branch("fPrimary4LHe", &fPrimary4LHe, "fPrimary4LHe/I");

  fTreeKF->Branch("fPrimVertexXKF", &fPrimVertexXKF, "fPrimVertexXKF/D");
  fTreeKF->Branch("fPrimVertexYKF", &fPrimVertexYKF, "fPrimVertexYKF/D");
  fTreeKF->Branch("fPrimVertexZKF", &fPrimVertexZKF, "fPrimVertexZKF/D");
  fTreeKF->Branch("fPrimVertexXErrKF", &fPrimVertexXErrKF, "fPrimVertexXErrKF/D");
  fTreeKF->Branch("fPrimVertexYErrKF", &fPrimVertexYErrKF, "fPrimVertexYErrKF/D");
  fTreeKF->Branch("fPrimVertexZErrKF", &fPrimVertexZErrKF, "fPrimVertexZErrKF/D");
  fTreeKF->Branch("fPrimVertChi2KF", &fPrimVertChi2KF, "fPrimVertChi2KF/D");
  fTreeKF->Branch("fPrimVertNDFKF", &fPrimVertNDFKF, "fPrimVertNDFKF/I");
  fTreeKF->Branch("fSecVertexXKF", &fSecVertexXKF, "fSecVertexXKF/D");
  fTreeKF->Branch("fSecVertexYKF", &fSecVertexYKF, "fSecVertexYKF/D");
  fTreeKF->Branch("fSecVertexZKF", &fSecVertexZKF, "fSecVertexZKF/D");
  fTreeKF->Branch("fSecVertexXErrKF", &fSecVertexXErrKF, "fSecVertexXErrKF/D");
  fTreeKF->Branch("fSecVertexYErrKF", &fSecVertexYErrKF, "fSecVertexYErrKF/D");
  fTreeKF->Branch("fSecVertexZErrKF", &fSecVertexZErrKF, "fSecVertexZErrKF/D");
  fTreeKF->Branch("fSecVertChi2KF", &fSecVertChi2KF, "fSecVertChi2KF/D");
  fTreeKF->Branch("fSecVertNDFKF", &fSecVertNDFKF, "fSecVertNDFKF/I");
  fTreeKF->Branch("fTertVertexXKF", &fTertVertexXKF, "fTertVertexXKF/D");
  fTreeKF->Branch("fTertVertexYKF", &fTertVertexYKF, "fTertVertexYKF/D");
  fTreeKF->Branch("fTertVertexZKF", &fTertVertexZKF, "fTertVertexZKF/D");
  fTreeKF->Branch("fTertVertexXErrKF", &fTertVertexXErrKF, "fTertVertexXErrKF/D");
  fTreeKF->Branch("fTertVertexYErrKF", &fTertVertexYErrKF, "fTertVertexYErrKF/D");
  fTreeKF->Branch("fTertVertexZErrKF", &fTertVertexZErrKF, "fTertVertexZErrKF/D");
  fTreeKF->Branch("fTertVertChi2KF", &fTertVertChi2KF, "fTertVertChi2KF/D");
  fTreeKF->Branch("fTertVertNDFKF", &fTertVertNDFKF, "fTertVertNDFKF/I");
  /*fTreeKF->Branch("fisOnlineV0_13", &fisOnlineV0_13, "fisOnlineV0_13/I");
    fTreeKF->Branch("fV0VertexX_13", &fV0VertexX_13, "fV0VertexX_13/D");
    fTreeKF->Branch("fV0VertexY_13", &fV0VertexY_13, "fV0VertexY_13/D");
    fTreeKF->Branch("fV0VertexZ_13", &fV0VertexZ_13, "fV0VertexZ_13/D");
    fTreeKF->Branch("fisOnlineV0_14", &fisOnlineV0_14, "fisOnlineV0_14/I");
    fTreeKF->Branch("fV0VertexX_14", &fV0VertexX_14, "fV0VertexX_14/D");
    fTreeKF->Branch("fV0VertexY_14", &fV0VertexY_14, "fV0VertexY_14/D");
    fTreeKF->Branch("fV0VertexZ_14", &fV0VertexZ_14, "fV0VertexZ_14/D");
    fTreeKF->Branch("fisOnlineV0_23", &fisOnlineV0_23, "fisOnlineV0_23/I");
    fTreeKF->Branch("fV0VertexX_23", &fV0VertexX_23, "fV0VertexX_23/D");
    fTreeKF->Branch("fV0VertexY_23", &fV0VertexY_23, "fV0VertexY_23/D");
    fTreeKF->Branch("fV0VertexZ_23", &fV0VertexZ_23, "fV0VertexZ_23/D");
    fTreeKF->Branch("fisOnlineV0_24", &fisOnlineV0_24, "fisOnlineV0_24/I");
    fTreeKF->Branch("fV0VertexX_24", &fV0VertexX_24, "fV0VertexX_24/D");
    fTreeKF->Branch("fV0VertexY_24", &fV0VertexY_24, "fV0VertexY_24/D");
    fTreeKF->Branch("fV0VertexZ_24", &fV0VertexZ_24, "fV0VertexZ_24/D");*/
  fTreeKF->Branch("fDCA2BXYKF", &fDCA2BXYKF, "fDCA2BXYKF/D");
  fTreeKF->Branch("fDCA3B1XYKF", &fDCA3B1XYKF, "fDCA3B1XYKF/D");
  fTreeKF->Branch("fDCA3B2XYKF", &fDCA3B2XYKF, "fDCA3B2XYKF/D");
  fTreeKF->Branch("fDCA3B3XYKF", &fDCA3B3XYKF, "fDCA3B3XYKF/D");
  fTreeKF->Branch("fDCA2BZKF", &fDCA2BZKF, "fDCA2BZKF/D");
  fTreeKF->Branch("fDCA3B1ZKF", &fDCA3B1ZKF, "fDCA3B1ZKF/D");
  fTreeKF->Branch("fDCA3B2ZKF", &fDCA3B2ZKF, "fDCA3B2ZKF/D");
  fTreeKF->Branch("fDCA3B3ZKF", &fDCA3B3ZKF, "fDCA3B3ZKF/D");
  fTreeKF->Branch("fPAKF", &fPAKF, "fPAKF/D");
  fTreeKF->Branch("fSubPAKF", &fSubPAKF, "fSubPAKF/D");
  fTreeKF->Branch("fSubPA2KF", &fSubPA2KF, "fSubPA2KF/D");
  fTreeKF->Branch("farmalpha", &farmalpha, "farmalpha/D");
  fTreeKF->Branch("farmpt", &farmpt, "farmpt/D");

  fTreeKF->Branch("fmMother", &fmMother, "fmMother/D");
  fTreeKF->Branch("fmMotherKF", &fmMotherKF, "fmMotherKF/D");
  fTreeKF->Branch("fmMotherErrKF", &fmMotherErrKF, "fmMotherErrKF/D");
  fTreeKF->Branch("fEMotherKF", &fEMotherKF, "fEMotherKF/D");
  fTreeKF->Branch("fEMotherErrKF", &fEMotherErrKF, "fEMotherErrKF/D");
  fTreeKF->Branch("fpxMotherKF", &fpxMotherKF, "fpxMotherKF/D");
  fTreeKF->Branch("fpxMotherErrKF", &fpxMotherErrKF, "fpxMotherErrKF/D");
  fTreeKF->Branch("fpyMotherKF", &fpyMotherKF, "fpyMotherKF/D");
  fTreeKF->Branch("fpyMotherErrKF", &fpyMotherErrKF, "fpyMotherErrKF/D");
  fTreeKF->Branch("fpzMotherKF", &fpzMotherKF, "fpzMotherKF/D");
  fTreeKF->Branch("fpzMotherErrKF", &fpzMotherErrKF, "fpzMotherErrKF/D");
  fTreeKF->Branch("fptMotherKF", &fptMotherKF, "fptMotherKF/D");
  fTreeKF->Branch("fptMotherErrKF", &fptMotherErrKF, "fptMotherErrKF/D");
  fTreeKF->Branch("fpMotherKF", &fpMotherKF, "fpMotherKF/D");
  fTreeKF->Branch("fpMotherErrKF", &fpMotherErrKF, "fpMotherErrKF/D");
  fTreeKF->Branch("fyMotherKF", &fyMotherKF, "fyMotherKF/D");
  fTreeKF->Branch("fctMotherKF", &fctMotherKF, "fctMotherKF/D");
  fTreeKF->Branch("fctMotherErrKF", &fctMotherErrKF, "fctMotherErrKF/D");
  fTreeKF->Branch("fmSubMother", &fmSubMother, "fmSubMother/D");
  fTreeKF->Branch("fmSubMotherKF", &fmSubMotherKF, "fmSubMotherKF/D");
  fTreeKF->Branch("fmSubMotherErrKF", &fmSubMotherErrKF, "fmSubMotherErrKF/D");
  fTreeKF->Branch("fESubMotherKF", &fESubMotherKF, "fESubMotherKF/D");
  fTreeKF->Branch("fESubMotherErrKF", &fESubMotherErrKF, "fESubMotherErrKF/D");
  fTreeKF->Branch("fpxSubMotherKF", &fpxSubMotherKF, "fpxSubMotherKF/D");
  fTreeKF->Branch("fpxSubMotherErrKF", &fpxSubMotherErrKF, "fpxSubMotherErrKF/D");
  fTreeKF->Branch("fpySubMotherKF", &fpySubMotherKF, "fpySubMotherKF/D");
  fTreeKF->Branch("fpySubMotherErrKF", &fpySubMotherErrKF, "fpySubMotherErrKF/D");
  fTreeKF->Branch("fpzSubMotherKF", &fpzSubMotherKF, "fpzSubMotherKF/D");
  fTreeKF->Branch("fpzSubMotherErrKF", &fpzSubMotherErrKF, "fpzSubMotherErrKF/D");
  fTreeKF->Branch("fptSubMotherKF", &fptSubMotherKF, "fptSubMotherKF/D");
  fTreeKF->Branch("fptSubMotherErrKF", &fptSubMotherErrKF, "fptSubMotherErrKF/D");
  fTreeKF->Branch("fpSubMotherKF", &fpSubMotherKF, "fpSubMotherKF/D");
  fTreeKF->Branch("fpSubMotherErrKF", &fpSubMotherErrKF, "fpSubMotherErrKF/D");
  fTreeKF->Branch("fySubMotherKF", &fySubMotherKF, "fySubMotherKF/D");
  fTreeKF->Branch("fctSubMotherKF", &fctSubMotherKF, "fctSubMotherKF/D");
  fTreeKF->Branch("fctSubMotherErrKF", &fctSubMotherErrKF, "fctSubMotherErrKF/D");

  fTreeKF->Branch("fEDaughter", &fEDaughter, "fEDaughter/D");
  fTreeKF->Branch("fptDaughter", &fptDaughter, "fptDaughter/D");
  fTreeKF->Branch("fpxDaughter", &fpxDaughter, "fpxDaughter/D");
  fTreeKF->Branch("fpyDaughter", &fpyDaughter, "fpyDaughter/D");
  fTreeKF->Branch("fpzDaughter", &fpzDaughter, "fpzDaughter/D");
  fTreeKF->Branch("fyDaughter", &fyDaughter, "fyDaughter/D");
  fTreeKF->Branch("fEDaughterKF", &fEDaughterKF, "fEDaughterKF/D");
  fTreeKF->Branch("fptDaughterKF", &fptDaughterKF, "fptDaughterKF/D");
  fTreeKF->Branch("fpxDaughterKF", &fpxDaughterKF, "fpxDaughterKF/D");
  fTreeKF->Branch("fpyDaughterKF", &fpyDaughterKF, "fpyDaughterKF/D");
  fTreeKF->Branch("fpzDaughterKF", &fpzDaughterKF, "fpzDaughterKF/D");
  fTreeKF->Branch("fyDaughterKF", &fyDaughterKF, "fyDaughterKF/D");
  fTreeKF->Branch("fEDaughter1", &fEDaughter1, "fEDaughter1/D");
  fTreeKF->Branch("fptDaughter1", &fptDaughter1, "fptDaughter1/D");
  fTreeKF->Branch("fpxDaughter1", &fpxDaughter1, "fpxDaughter1/D");
  fTreeKF->Branch("fpyDaughter1", &fpyDaughter1, "fpyDaughter1/D");
  fTreeKF->Branch("fpzDaughter1", &fpzDaughter1, "fpzDaughter1/D");
  fTreeKF->Branch("fyDaughter1", &fyDaughter1, "fyDaughter1/D");
  fTreeKF->Branch("fEDaughter1KF", &fEDaughter1KF, "fEDaughter1KF/D");
  fTreeKF->Branch("fptDaughter1KF", &fptDaughter1KF, "fptDaughter1KF/D");
  fTreeKF->Branch("fpxDaughter1KF", &fpxDaughter1KF, "fpxDaughter1KF/D");
  fTreeKF->Branch("fpyDaughter1KF", &fpyDaughter1KF, "fpyDaughter1KF/D");
  fTreeKF->Branch("fpzDaughter1KF", &fpzDaughter1KF, "fpzDaughter1KF/D");
  fTreeKF->Branch("fyDaughter1KF", &fyDaughter1KF, "fyDaughter1KF/D");
  fTreeKF->Branch("fEDaughter2", &fEDaughter2, "fEDaughter2/D");
  fTreeKF->Branch("fptDaughter2", &fptDaughter2, "fptDaughter2/D");
  fTreeKF->Branch("fpxDaughter2", &fpxDaughter2, "fpxDaughter2/D");
  fTreeKF->Branch("fpyDaughter2", &fpyDaughter2, "fpyDaughter2/D");
  fTreeKF->Branch("fpzDaughter2", &fpzDaughter2, "fpzDaughter2/D");
  fTreeKF->Branch("fyDaughter2", &fyDaughter2, "fyDaughter2/D");
  fTreeKF->Branch("fEDaughter2KF", &fEDaughter2KF, "fEDaughter2KF/D");
  fTreeKF->Branch("fptDaughter2KF", &fptDaughter2KF, "fptDaughter2KF/D");
  fTreeKF->Branch("fpxDaughter2KF", &fpxDaughter2KF, "fpxDaughter2KF/D");
  fTreeKF->Branch("fpyDaughter2KF", &fpyDaughter2KF, "fpyDaughter2KF/D");
  fTreeKF->Branch("fpzDaughter2KF", &fpzDaughter2KF, "fpzDaughter2KF/D");
  fTreeKF->Branch("fyDaughter2KF", &fyDaughter2KF, "fyDaughter2KF/D");
  fTreeKF->Branch("fEDaughter3", &fEDaughter3, "fEDaughter3/D");
  fTreeKF->Branch("fptDaughter3", &fptDaughter3, "fptDaughter3/D");
  fTreeKF->Branch("fpxDaughter3", &fpxDaughter3, "fpxDaughter3/D");
  fTreeKF->Branch("fpyDaughter3", &fpyDaughter3, "fpyDaughter3/D");
  fTreeKF->Branch("fpzDaughter3", &fpzDaughter3, "fpzDaughter3/D");
  fTreeKF->Branch("fyDaughter3", &fyDaughter3, "fyDaughter3/D");
  fTreeKF->Branch("fEDaughter3KF", &fEDaughter3KF, "fEDaughter3KF/D");
  fTreeKF->Branch("fptDaughter3KF", &fptDaughter3KF, "fptDaughter3KF/D");
  fTreeKF->Branch("fpxDaughter3KF", &fpxDaughter3KF, "fpxDaughter3KF/D");
  fTreeKF->Branch("fpyDaughter3KF", &fpyDaughter3KF, "fpyDaughter3KF/D");
  fTreeKF->Branch("fpzDaughter3KF", &fpzDaughter3KF, "fpzDaughter3KF/D");
  fTreeKF->Branch("fyDaughter3KF", &fyDaughter3KF, "fyDaughter3KF/D");  

  fTreeKF->Branch("fDcaDaughterXYKF", &fDcaDaughterXYKF, "fDcaDaughterXYKF/D");
  fTreeKF->Branch("fDcaDaughter1XYKF", &fDcaDaughter1XYKF, "fDcaDaughter1XYKF/D");
  fTreeKF->Branch("fDcaDaughter2XYKF", &fDcaDaughter2XYKF, "fDcaDaughter2XYKF/D");
  fTreeKF->Branch("fDcaDaughter3XYKF", &fDcaDaughter3XYKF, "fDcaDaughter3XYKF/D");
  fTreeKF->Branch("fDcaDaughter4XYKF", &fDcaDaughter4XYKF, "fDcaDaughter4XYKF/D");
  fTreeKF->Branch("fDcaDaughterZKF", &fDcaDaughterZKF, "fDcaDaughterZKF/D");
  fTreeKF->Branch("fDcaDaughter1ZKF", &fDcaDaughter1ZKF, "fDcaDaughter1ZKF/D");
  fTreeKF->Branch("fDcaDaughter2ZKF", &fDcaDaughter2ZKF, "fDcaDaughter2ZKF/D");
  fTreeKF->Branch("fDcaDaughter3ZKF", &fDcaDaughter3ZKF, "fDcaDaughter3ZKF/D");
  fTreeKF->Branch("fDcaDaughter4ZKF", &fDcaDaughter4ZKF, "fDcaDaughter4ZKF/D");
  fTreeKF->Branch("fDcaSecDaughterXYKF", &fDcaSecDaughterXYKF, "fDcaSecDaughterXYKF/D");
  fTreeKF->Branch("fDcaSecDaughterZKF", &fDcaSecDaughterZKF, "fDcaSecDaughterZKF/D");
  fTreeKF->Branch("fDcaSecDaughter1XYKF", &fDcaSecDaughter1XYKF, "fDcaSecDaughter1XYKF/D");
  fTreeKF->Branch("fDcaSecDaughter1ZKF", &fDcaSecDaughter1ZKF, "fDcaSecDaughter1ZKF/D");
  fTreeKF->Branch("fDcaSecDaughter2XYKF", &fDcaSecDaughter2XYKF, "fDcaSecDaughter2XYKF/D");
  fTreeKF->Branch("fDcaSecDaughter2ZKF", &fDcaSecDaughter2ZKF, "fDcaSecDaughter2ZKF/D");
  fTreeKF->Branch("fDcaSecDaughter3XYKF", &fDcaSecDaughter3XYKF, "fDcaSecDaughter3XYKF/D");
  fTreeKF->Branch("fDcaSecDaughter3ZKF", &fDcaSecDaughter3ZKF, "fDcaSecDaughter3ZKF/D");
  fTreeKF->Branch("fDcaSecDaughter4XYKF", &fDcaSecDaughter4XYKF, "fDcaSecDaughter4XYKF/D");
  fTreeKF->Branch("fDcaSecDaughter4ZKF", &fDcaSecDaughter4ZKF, "fDcaSecDaughter4ZKF/D");

  fTreeKF->Branch("fPropDCADaughter", &fPropDCADaughter, "fPropDCADaughter/I");
  fTreeKF->Branch("fImParDaughter", &fImParDaughter, "fImParDaughter/D");
  fTreeKF->Branch("fImParzDaughter", &fImParzDaughter, "fImParzDaughter/D");
  fTreeKF->Branch("fPropDCADaughter1", &fPropDCADaughter1, "fPropDCADaughter1/I");
  fTreeKF->Branch("fImParDaughter1", &fImParDaughter1, "fImParDaughter1/D");
  fTreeKF->Branch("fImParzDaughter1", &fImParzDaughter1, "fImParzDaughter1/D");
  fTreeKF->Branch("fPropDCADaughter2", &fPropDCADaughter2, "fPropDCADaughter2/I");
  fTreeKF->Branch("fImParDaughter2", &fImParDaughter2, "fImParDaughter2/D");
  fTreeKF->Branch("fImParzDaughter2", &fImParzDaughter2, "fImParzDaughter2/D");
  fTreeKF->Branch("fPropDCADaughter3", &fPropDCADaughter3, "fPropDCADaughter3/I");
  fTreeKF->Branch("fImParDaughter3", &fImParDaughter3, "fImParDaughter3/D");
  fTreeKF->Branch("fImParzDaughter3", &fImParzDaughter3, "fImParzDaughter3/D");
  fTreeKF->Branch("fPropDCADaughter4", &fPropDCADaughter4, "fPropDCADaughter4/I");
  fTreeKF->Branch("fImParDaughter4", &fImParDaughter4, "fImParDaughter4/D");
  fTreeKF->Branch("fImParzDaughter4", &fImParzDaughter4, "fImParzDaughter4/D");
  fTreeKF->Branch("fDcaSecDaughter", &fDcaSecDaughter, "fDcaSecDaughter/D");
  fTreeKF->Branch("fDcaSecDaughter1", &fDcaSecDaughter1, "fDcaSecDaughter1/D");
  fTreeKF->Branch("fDcaSecDaughter2", &fDcaSecDaughter2, "fDcaSecDaughter2/D");
  fTreeKF->Branch("fDcaSecDaughter3", &fDcaSecDaughter3, "fDcaSecDaughter3/D");
  fTreeKF->Branch("fDcaSecDaughter4", &fDcaSecDaughter4, "fDcaSecDaughter4/D");

  fTreeKF->Branch("fTrackPIDDaughter", &fTrackPIDDaughter, "fTrackPIDDaughter/I");
  fTreeKF->Branch("fLabelDaughterKF", &fLabelDaughterKF, "fLabelDaughterKF/D");
  fTreeKF->Branch("fpDaughter", &fpDaughter, "fpDaughter/D");
  fTreeKF->Branch("fptDaughterUnProp", &fptDaughterUnProp, "fptDaughterUnProp/D");
  fTreeKF->Branch("fpxDaughterUnProp", &fpxDaughterUnProp, "fpxDaughterUnProp/D");
  fTreeKF->Branch("fpyDaughterUnProp", &fpyDaughterUnProp, "fpyDaughterUnProp/D");
  fTreeKF->Branch("fpzDaughterUnProp", &fpzDaughterUnProp, "fpzDaughterUnProp/D");
  fTreeKF->Branch("fdEdxDaughter", &fdEdxDaughter, "fdEdxDaughter/D");
  fTreeKF->Branch("fdEdxSigmaDaughter", &fdEdxSigmaDaughter, "fdEdxSigmaDaughter/D");
  fTreeKF->Branch("fDcaDaughter", &fDcaDaughter, "fDcaDaughter/D");
  fTreeKF->Branch("fDcaDaughtero", &fDcaDaughtero, "fDcaDaughtero/D");
  fTreeKF->Branch("fDcazDaughter", &fDcazDaughter, "fDcazDaughter/D");
  fTreeKF->Branch("fNclsDaughter", &fNclsDaughter, "fNclsDaughter/I");
  fTreeKF->Branch("fChi2Daughter", &fChi2Daughter, "fChi2Daughter/D");
  fTreeKF->Branch("fNclsITSDaughter", &fNclsITSDaughter, "fNclsITSDaughter/I");
  fTreeKF->Branch("fEtaDaughter", &fEtaDaughter, "fEtaDaughter/D");
  fTreeKF->Branch("fPhiDaughter", &fPhiDaughter, "fPhiDaughter/D");
  fTreeKF->Branch("fGeoLengthDaughter", &fGeoLengthDaughter, "fGeoLengthDaughter/D");
  fTreeKF->Branch("fTOFSignalDaughter", &fTOFSignalDaughter, "fTOFSignalDaughter/D");
  fTreeKF->Branch("fSigmaYXDaughter", &fSigmaYXDaughter, "fSigmaYXDaughter/D");
  fTreeKF->Branch("fSigmaXYZDaughter", &fSigmaXYZDaughter, "fSigmaXYZDaughter/D");
  fTreeKF->Branch("fSigmaZDaughter", &fSigmaZDaughter, "fSigmaZDaughter/D");
  fTreeKF->Branch("fPtUncertDaughter", &fPtUncertDaughter, "fPtUncertDaughter/D");
  fTreeKF->Branch("fTPCRefitDaughter", &fTPCRefitDaughter, "fTPCRefitDaughter/I");
  fTreeKF->Branch("fITSRefitDaughter", &fITSRefitDaughter, "fITSRefitDaughter/I");
  fTreeKF->Branch("fITSLayer1Daughter", &fITSLayer1Daughter, "fITSLayer1Daughter/I");
  fTreeKF->Branch("fITSLayer2Daughter", &fITSLayer2Daughter, "fITSLayer2Daughter/I");
  fTreeKF->Branch("fITSLayer3Daughter", &fITSLayer3Daughter, "fITSLayer3Daughter/I");
  fTreeKF->Branch("fITSLayer4Daughter", &fITSLayer4Daughter, "fITSLayer4Daughter/I");
  fTreeKF->Branch("fITSLayer5Daughter", &fITSLayer5Daughter, "fITSLayer5Daughter/I");
  fTreeKF->Branch("fITSLayer6Daughter", &fITSLayer6Daughter, "fITSLayer6Daughter/I");

  fTreeKF->Branch("fTrackPIDDaughter1", &fTrackPIDDaughter1, "fTrackPIDDaughter1/I");
  fTreeKF->Branch("fLabelDaughter1KF", &fLabelDaughter1KF, "fLabelDaughter1KF/D");
  fTreeKF->Branch("fpDaughter1", &fpDaughter1, "fpDaughter1/D");
  fTreeKF->Branch("fptDaughter1UnProp", &fptDaughter1UnProp, "fptDaughter1UnProp/D");
  fTreeKF->Branch("fpxDaughter1UnProp", &fpxDaughter1UnProp, "fpxDaughter1UnProp/D");
  fTreeKF->Branch("fpyDaughter1UnProp", &fpyDaughter1UnProp, "fpyDaughter1UnProp/D");
  fTreeKF->Branch("fpzDaughter1UnProp", &fpzDaughter1UnProp, "fpzDaughter1UnProp/D");
  fTreeKF->Branch("fdEdxDaughter1", &fdEdxDaughter1, "fdEdxDaughter1/D");
  fTreeKF->Branch("fdEdxSigmaDaughter1", &fdEdxSigmaDaughter1, "fdEdxSigmaDaughter1/D");
  fTreeKF->Branch("fDcaDaughter1", &fDcaDaughter1, "fDcaDaughter1/D");
  fTreeKF->Branch("fDcaDaughter1o", &fDcaDaughter1o, "fDcaDaughter1o/D");
  fTreeKF->Branch("fDcazDaughter1", &fDcazDaughter1, "fDcazDaughter1/D");
  fTreeKF->Branch("fNclsDaughter1", &fNclsDaughter1, "fNclsDaughter1/I");
  fTreeKF->Branch("fChi2Daughter1", &fChi2Daughter1, "fChi2Daughter1/D");
  fTreeKF->Branch("fNclsITSDaughter1", &fNclsITSDaughter1, "fNclsITSDaughter1/I");
  fTreeKF->Branch("fEtaDaughter1", &fEtaDaughter1, "fEtaDaughter1/D");
  fTreeKF->Branch("fPhiDaughter1", &fPhiDaughter1, "fPhiDaughter1/D");
  fTreeKF->Branch("fGeoLengthDaughter1", &fGeoLengthDaughter1, "fGeoLengthDaughter1/D");
  fTreeKF->Branch("fTOFSignalDaughter1", &fTOFSignalDaughter1, "fTOFSignalDaughter1/D");
  fTreeKF->Branch("fSigmaYXDaughter1", &fSigmaYXDaughter1, "fSigmaYXDaughter1/D");
  fTreeKF->Branch("fSigmaXYZDaughter1", &fSigmaXYZDaughter1, "fSigmaXYZDaughter1/D");
  fTreeKF->Branch("fSigmaZDaughter1", &fSigmaZDaughter1, "fSigmaZDaughter1/D");
  fTreeKF->Branch("fPtUncertDaughter1", &fPtUncertDaughter1, "fPtUncertDaughter1/D");
  fTreeKF->Branch("fTPCRefitDaughter1", &fTPCRefitDaughter1, "fTPCRefitDaughter1/I");
  fTreeKF->Branch("fITSRefitDaughter1", &fITSRefitDaughter1, "fITSRefitDaughter1/I");
  fTreeKF->Branch("fITSLayer1Daughter1", &fITSLayer1Daughter1, "fITSLayer1Daughter1/I");
  fTreeKF->Branch("fITSLayer2Daughter1", &fITSLayer2Daughter1, "fITSLayer2Daughter1/I");
  fTreeKF->Branch("fITSLayer3Daughter1", &fITSLayer3Daughter1, "fITSLayer3Daughter1/I");
  fTreeKF->Branch("fITSLayer4Daughter1", &fITSLayer4Daughter1, "fITSLayer4Daughter1/I");
  fTreeKF->Branch("fITSLayer5Daughter1", &fITSLayer5Daughter1, "fITSLayer5Daughter1/I");
  fTreeKF->Branch("fITSLayer6Daughter1", &fITSLayer6Daughter1, "fITSLayer6Daughter1/I");

  fTreeKF->Branch("fTrackPIDDaughter2", &fTrackPIDDaughter2, "fTrackPIDDaughter2/I");
  fTreeKF->Branch("fLabelDaughter2KF", &fLabelDaughter2KF, "fLabelDaughter2KF/D");
  fTreeKF->Branch("fpDaughter2", &fpDaughter2, "fpDaughter2/D");
  fTreeKF->Branch("fptDaughter2UnProp", &fptDaughter2UnProp, "fptDaughter2UnProp/D");
  fTreeKF->Branch("fpxDaughter2UnProp", &fpxDaughter2UnProp, "fpxDaughter2UnProp/D");
  fTreeKF->Branch("fpyDaughter2UnProp", &fpyDaughter2UnProp, "fpyDaughter2UnProp/D");
  fTreeKF->Branch("fpzDaughter2UnProp", &fpzDaughter2UnProp, "fpzDaughter2UnProp/D");
  fTreeKF->Branch("fdEdxDaughter2", &fdEdxDaughter2, "fdEdxDaughter2/D");
  fTreeKF->Branch("fdEdxSigmaDaughter2", &fdEdxSigmaDaughter2, "fdEdxSigmaDaughter2/D");
  fTreeKF->Branch("fDcaDaughter2", &fDcaDaughter2, "fDcaDaughter2/D");
  fTreeKF->Branch("fDcaDaughter2o", &fDcaDaughter2o, "fDcaDaughter2o/D");
  fTreeKF->Branch("fDcazDaughter2", &fDcazDaughter2, "fDcazDaughter2/D");
  fTreeKF->Branch("fNclsDaughter2", &fNclsDaughter2, "fNclsDaughter2/I");
  fTreeKF->Branch("fChi2Daughter2", &fChi2Daughter2, "fChi2Daughter2/D");
  fTreeKF->Branch("fNclsITSDaughter2", &fNclsITSDaughter2, "fNclsITSDaughter2/I");
  fTreeKF->Branch("fEtaDaughter2", &fEtaDaughter2, "fEtaDaughter2/D");
  fTreeKF->Branch("fPhiDaughter2", &fPhiDaughter2, "fPhiDaughter2/D");
  fTreeKF->Branch("fGeoLengthDaughter2", &fGeoLengthDaughter2, "fGeoLengthDaughter2/D");
  fTreeKF->Branch("fTOFSignalDaughter2", &fTOFSignalDaughter2, "fTOFSignalDaughter2/D");
  fTreeKF->Branch("fSigmaYXDaughter2", &fSigmaYXDaughter2, "fSigmaYXDaughter2/D");
  fTreeKF->Branch("fSigmaXYZDaughter2", &fSigmaXYZDaughter2, "fSigmaXYZDaughter2/D");
  fTreeKF->Branch("fSigmaZDaughter2", &fSigmaZDaughter2, "fSigmaZDaughter2/D");
  fTreeKF->Branch("fPtUncertDaughter2", &fPtUncertDaughter2, "fPtUncertDaughter2/D");
  fTreeKF->Branch("fTPCRefitDaughter2", &fTPCRefitDaughter2, "fTPCRefitDaughter2/I");
  fTreeKF->Branch("fITSRefitDaughter2", &fITSRefitDaughter2, "fITSRefitDaughter2/I");
  fTreeKF->Branch("fITSLayer1Daughter2", &fITSLayer1Daughter2, "fITSLayer1Daughter2/I");
  fTreeKF->Branch("fITSLayer2Daughter2", &fITSLayer2Daughter2, "fITSLayer2Daughter2/I");
  fTreeKF->Branch("fITSLayer3Daughter2", &fITSLayer3Daughter2, "fITSLayer3Daughter2/I");
  fTreeKF->Branch("fITSLayer4Daughter2", &fITSLayer4Daughter2, "fITSLayer4Daughter2/I");
  fTreeKF->Branch("fITSLayer5Daughter2", &fITSLayer5Daughter2, "fITSLayer5Daughter2/I");
  fTreeKF->Branch("fITSLayer6Daughter2", &fITSLayer6Daughter2, "fITSLayer6Daughter2/I");

  fTreeKF->Branch("fTrackPIDDaughter3", &fTrackPIDDaughter3, "fTrackPIDDaughter3/I");
  fTreeKF->Branch("fLabelDaughter3KF", &fLabelDaughter3KF, "fLabelDaughter3KF/D");
  fTreeKF->Branch("fpDaughter3", &fpDaughter3, "fpDaughter3/D");
  fTreeKF->Branch("fptDaughter3UnProp", &fptDaughter3UnProp, "fptDaughter3UnProp/D");
  fTreeKF->Branch("fpxDaughter3UnProp", &fpxDaughter3UnProp, "fpxDaughter3UnProp/D");
  fTreeKF->Branch("fpyDaughter3UnProp", &fpyDaughter3UnProp, "fpyDaughter3UnProp/D");
  fTreeKF->Branch("fpzDaughter3UnProp", &fpzDaughter3UnProp, "fpzDaughter3UnProp/D");
  fTreeKF->Branch("fdEdxDaughter3", &fdEdxDaughter3, "fdEdxDaughter3/D");
  fTreeKF->Branch("fdEdxSigmaDaughter3", &fdEdxSigmaDaughter3, "fdEdxSigmaDaughter3/D");
  fTreeKF->Branch("fDcaDaughter3", &fDcaDaughter3, "fDcaDaughter3/D");
  fTreeKF->Branch("fDcaDaughter3o", &fDcaDaughter3o, "fDcaDaughter3o/D");
  fTreeKF->Branch("fDcazDaughter3", &fDcazDaughter3, "fDcazDaughter3/D");
  fTreeKF->Branch("fNclsDaughter3", &fNclsDaughter3, "fNclsDaughter3/I");
  fTreeKF->Branch("fChi2Daughter3", &fChi2Daughter3, "fChi2Daughter3/D");
  fTreeKF->Branch("fNclsITSDaughter3", &fNclsITSDaughter3, "fNclsITSDaughter3/I");
  fTreeKF->Branch("fEtaDaughter3", &fEtaDaughter3, "fEtaDaughter3/D");
  fTreeKF->Branch("fPhiDaughter3", &fPhiDaughter3, "fPhiDaughter3/D");
  fTreeKF->Branch("fGeoLengthDaughter3", &fGeoLengthDaughter3, "fGeoLengthDaughter3/D");
  fTreeKF->Branch("fTOFSignalDaughter3", &fTOFSignalDaughter3, "fTOFSignalDaughter3/D");
  fTreeKF->Branch("fSigmaYXDaughter3", &fSigmaYXDaughter3, "fSigmaYXDaughter3/D");
  fTreeKF->Branch("fSigmaXYZDaughter3", &fSigmaXYZDaughter3, "fSigmaXYZDaughter3/D");
  fTreeKF->Branch("fSigmaZDaughter3", &fSigmaZDaughter3, "fSigmaZDaughter3/D");
  fTreeKF->Branch("fPtUncertDaughter3", &fPtUncertDaughter3, "fPtUncertDaughter3/D");
  fTreeKF->Branch("fTPCRefitDaughter3", &fTPCRefitDaughter3, "fTPCRefitDaughter3/I");
  fTreeKF->Branch("fITSRefitDaughter3", &fITSRefitDaughter3, "fITSRefitDaughter3/I");
  fTreeKF->Branch("fITSLayer1Daughter3", &fITSLayer1Daughter3, "fITSLayer1Daughter3/I");
  fTreeKF->Branch("fITSLayer2Daughter3", &fITSLayer2Daughter3, "fITSLayer2Daughter3/I");
  fTreeKF->Branch("fITSLayer3Daughter3", &fITSLayer3Daughter3, "fITSLayer3Daughter3/I");
  fTreeKF->Branch("fITSLayer4Daughter3", &fITSLayer4Daughter3, "fITSLayer4Daughter3/I");
  fTreeKF->Branch("fITSLayer5Daughter3", &fITSLayer5Daughter3, "fITSLayer5Daughter3/I");
  fTreeKF->Branch("fITSLayer6Daughter3", &fITSLayer6Daughter3, "fITSLayer6Daughter3/I");

  fTreeKF->Branch("fEDaughter4", &fEDaughter4, "fEDaughter4/D");
  fTreeKF->Branch("fpDaughter4", &fpDaughter4, "fpDaughter4/D");
  fTreeKF->Branch("fptDaughter4", &fptDaughter4, "fptDaughter4/D");
  fTreeKF->Branch("fpxDaughter4", &fpxDaughter4, "fpxDaughter4/D");
  fTreeKF->Branch("fpyDaughter4", &fpyDaughter4, "fpyDaughter4/D");
  fTreeKF->Branch("fpzDaughter4", &fpzDaughter4, "fpzDaughter4/D");
  fTreeKF->Branch("fyDaughter4", &fyDaughter4, "fyDaughter4/D");
  fTreeKF->Branch("fDcaDaughter4", &fDcaDaughter4, "fDcaDaughter4/D");
  fTreeKF->Branch("fDcazDaughter4", &fDcazDaughter4, "fDcazDaughter4/D");
  fTreeKF->Branch("fSigmaYXDaughter4", &fSigmaYXDaughter4, "fSigmaYXDaughter4/D");
  fTreeKF->Branch("fSigmaXYZDaughter4", &fSigmaXYZDaughter4, "fSigmaXYZDaughter4/D");
  fTreeKF->Branch("fSigmaZDaughter4", &fSigmaZDaughter4, "fSigmaZDaughter4/D");
  fTreeKF->Branch("fPtUncertDaughter4", &fPtUncertDaughter4, "fPtUncertDaughter4/D");

  fTreeKF->Branch("fdEdxSigmaPion", &fdEdxSigmaPion, "fdEdxSigmaPion/D");
  fTreeKF->Branch("fdEdxSigmaDeuteron", &fdEdxSigmaDeuteron, "fdEdxSigmaDeuteron/D");
  fTreeKF->Branch("fdEdxSigmaTriton", &fdEdxSigmaTriton, "fdEdxSigmaTriton/D");
  fTreeKF->Branch("fdEdxSigmaAlpha", &fdEdxSigmaAlpha, "fdEdxSigmaAlpha/D");

  fTreeKF->Branch("fCovMatrixTrack", fCovMatrixTrack, "fCovMatrixTrack[21]/D");
  fTreeKF->Branch("fCovMatrixTrack1", fCovMatrixTrack1, "fCovMatrixTrack1[21]/D");
  fTreeKF->Branch("fCovMatrixTrack2", fCovMatrixTrack2, "fCovMatrixTrack2[21]/D");
  fTreeKF->Branch("fCovMatrixTrack3", fCovMatrixTrack3, "fCovMatrixTrack3[21]/D");
  fTreeKF->Branch("fCovMatrixTrack4", fCovMatrixTrack4, "fCovMatrixTrack4[21]/D");
  fTreeKF->Branch("fTrackPar", fTrackPar, "fTrackPar[7]/D");
  fTreeKF->Branch("fTrackPar1", fTrackPar1, "fTrackPar1[7]/D");
  fTreeKF->Branch("fTrackPar2", fTrackPar2, "fTrackPar2[7]/D");
  fTreeKF->Branch("fTrackPar3", fTrackPar3, "fTrackPar3[7]/D");
  fTreeKF->Branch("fTrackPar4", fTrackPar4, "fTrackPar4[7]/D");

  // _________________________________________________ //
  // __ associated 4LHe tree __ //
  gTree = new TTree("gTree", "gTree");
  gTree->Branch("fEventID", &fEventID, "fEventID/I");  
  gTree->Branch("fPeriod", &fPeriod, "fPeriod/I");
  gTree->Branch("frunnumber", &frunnumber, "frunnumber/I");
  gTree->Branch("feventclass", &feventclass, "feventclass/I");
  gTree->Branch("fCentrality", &fCentrality, "fCentrality/I");
  gTree->Branch("fMagneticField", &fMagneticField, "fMagneticField/I");
  gTree->Branch("fisWOTrackVertex", &fisWOTrackVertex, "fisWOTrackVertex/I");
  gTree->Branch("fmctruth", &fmctruth, "fmctruth/I");
  gTree->Branch("fTrigMB", &fTrigMB, "fTrigMB/I");
  gTree->Branch("fTrigHMV0", &fTrigHMV0, "fTrigHMV0/I");
  gTree->Branch("fTrigHMSPD", &fTrigHMSPD, "fTrigHMSPD/I");
  gTree->Branch("fTrigHNU", &fTrigHNU, "fTrigHNU/I");
  gTree->Branch("fTrigHQU", &fTrigHQU, "fTrigHQU/I");
  gTree->Branch("fTrigkCentral", &fTrigkCentral, "fTrigkCentral/I");
  gTree->Branch("fTrigkSemiCentral", &fTrigkSemiCentral, "fTrigkSemiCentral/I");

  gTree->Branch("fParticleID", &fParticleID, "fParticleID/I");
  gTree->Branch("fisExcited", &fisExcited, "fisExcited/I");
  gTree->Branch("fPDGMother", &fPDGMother, "fPDGMother/I");
  gTree->Branch("fChargeMother", &fChargeMother, "fChargeMother/I");
  gTree->Branch("fDecayChannel", &fDecayChannel, "fDecayChannel/I");
  gTree->Branch("fRecoMethod", &fRecoMethod, "fRecoMethod/I");
  gTree->Branch("fPrimary4LHe", &fPrimary4LHe, "fPrimary4LHe/I");  
  
  gTree->Branch("fPrimVertexX", &fPrimVertexX, "fPrimVertexX/D");
  gTree->Branch("fPrimVertexY", &fPrimVertexY, "fPrimVertexY/D");
  gTree->Branch("fPrimVertexZ", &fPrimVertexZ, "fPrimVertexZ/D");
  gTree->Branch("fTertVertexX", &fTertVertexX, "fTertVertexX/D");
  gTree->Branch("fTertVertexY", &fTertVertexY, "fTertVertexY/D");
  gTree->Branch("fTertVertexZ", &fTertVertexZ, "fTertVertexZ/D");
  gTree->Branch("fDCA3B1", &fDCA3B1, "fDCA3B1/D");
  gTree->Branch("fDCA3B2", &fDCA3B2, "fDCA3B2/D");
  gTree->Branch("fDCA3B3", &fDCA3B3, "fDCA3B3/D");
  gTree->Branch("fSubPA", &fSubPA, "fSubPA/D");
 
  gTree->Branch("fmSubMother", &fmSubMother, "fmSubMother/D");
  gTree->Branch("fESubMother", &fESubMother, "fESubMother/D");
  gTree->Branch("fpxSubMother", &fpxSubMother, "fpxSubMother/D");
  gTree->Branch("fpySubMother", &fpySubMother, "fpySubMother/D");
  gTree->Branch("fpzSubMother", &fpzSubMother, "fpzSubMother/D");
  gTree->Branch("fptSubMother", &fptSubMother, "fptSubMother/D");
  gTree->Branch("fpSubMother", &fpSubMother, "fpSubMother/D");
  gTree->Branch("fySubMother", &fySubMother, "fySubMother/D");
  gTree->Branch("fctSubMother", &fctSubMother, "fctSubMother/D");

  gTree->Branch("fEDaughter", &fEDaughter, "fEDaughter/D");
  gTree->Branch("fptDaughter", &fptDaughter, "fptDaughter/D");
  gTree->Branch("fpxDaughter", &fpxDaughter, "fpxDaughter/D");
  gTree->Branch("fpyDaughter", &fpyDaughter, "fpyDaughter/D");
  gTree->Branch("fpzDaughter", &fpzDaughter, "fpzDaughter/D");
  gTree->Branch("fyDaughter", &fyDaughter, "fyDaughter/D");
  gTree->Branch("fEDaughter1", &fEDaughter1, "fEDaughter1/D");
  gTree->Branch("fptDaughter1", &fptDaughter1, "fptDaughter1/D");
  gTree->Branch("fpxDaughter1", &fpxDaughter1, "fpxDaughter1/D");
  gTree->Branch("fpyDaughter1", &fpyDaughter1, "fpyDaughter1/D");
  gTree->Branch("fpzDaughter1", &fpzDaughter1, "fpzDaughter1/D");
  gTree->Branch("fyDaughter1", &fyDaughter1, "fyDaughter1/D");
  gTree->Branch("fEDaughter2", &fEDaughter2, "fEDaughter2/D");
  gTree->Branch("fptDaughter2", &fptDaughter2, "fptDaughter2/D");
  gTree->Branch("fpxDaughter2", &fpxDaughter2, "fpxDaughter2/D");
  gTree->Branch("fpyDaughter2", &fpyDaughter2, "fpyDaughter2/D");
  gTree->Branch("fpzDaughter2", &fpzDaughter2, "fpzDaughter2/D");
  gTree->Branch("fyDaughter2", &fyDaughter2, "fyDaughter2/D");

  gTree->Branch("fPropDCADaughter", &fPropDCADaughter, "fPropDCADaughter/I");
  gTree->Branch("fImParDaughter", &fImParDaughter, "fImParDaughter/D");
  gTree->Branch("fImParzDaughter", &fImParzDaughter, "fImParzDaughter/D");
  gTree->Branch("fPropDCADaughter1", &fPropDCADaughter1, "fPropDCADaughter1/I");
  gTree->Branch("fImParDaughter1", &fImParDaughter1, "fImParDaughter1/D");
  gTree->Branch("fImParzDaughter1", &fImParzDaughter1, "fImParzDaughter1/D");
  gTree->Branch("fPropDCADaughter2", &fPropDCADaughter2, "fPropDCADaughter2/I");
  gTree->Branch("fImParDaughter2", &fImParDaughter2, "fImParDaughter2/D");
  gTree->Branch("fImParzDaughter2", &fImParzDaughter2, "fImParzDaughter2/D");
  gTree->Branch("fDcaSecDaughter", &fDcaSecDaughter, "fDcaSecDaughter/D");
  gTree->Branch("fDcaSecDaughter1", &fDcaSecDaughter1, "fDcaSecDaughter1/D");
  gTree->Branch("fDcaSecDaughter2", &fDcaSecDaughter2, "fDcaSecDaughter2/D");  

  gTree->Branch("fTrackPIDDaughter", &fTrackPIDDaughter, "fTrackPIDDaughter/I");
  gTree->Branch("fLabelDaughterKF", &fLabelDaughterKF, "fLabelDaughterKF/D");
  gTree->Branch("fpDaughter", &fpDaughter, "fpDaughter/D");
  gTree->Branch("fptDaughterUnProp", &fptDaughterUnProp, "fptDaughterUnProp/D");
  gTree->Branch("fpxDaughterUnProp", &fpxDaughterUnProp, "fpxDaughterUnProp/D");
  gTree->Branch("fpyDaughterUnProp", &fpyDaughterUnProp, "fpyDaughterUnProp/D");
  gTree->Branch("fpzDaughterUnProp", &fpzDaughterUnProp, "fpzDaughterUnProp/D");
  gTree->Branch("fdEdxDaughter", &fdEdxDaughter, "fdEdxDaughter/D");
  gTree->Branch("fdEdxSigmaDaughter", &fdEdxSigmaDaughter, "fdEdxSigmaDaughter/D");
  gTree->Branch("fDcaDaughter", &fDcaDaughter, "fDcaDaughter/D");
  gTree->Branch("fDcaDaughtero", &fDcaDaughtero, "fDcaDaughtero/D");
  gTree->Branch("fDcazDaughter", &fDcazDaughter, "fDcazDaughter/D");
  gTree->Branch("fNclsDaughter", &fNclsDaughter, "fNclsDaughter/I");
  gTree->Branch("fChi2Daughter", &fChi2Daughter, "fChi2Daughter/D");
  gTree->Branch("fNclsITSDaughter", &fNclsITSDaughter, "fNclsITSDaughter/I");
  gTree->Branch("fEtaDaughter", &fEtaDaughter, "fEtaDaughter/D");
  gTree->Branch("fPhiDaughter", &fPhiDaughter, "fPhiDaughter/D");
  gTree->Branch("fGeoLengthDaughter", &fGeoLengthDaughter, "fGeoLengthDaughter/D");
  gTree->Branch("fTOFSignalDaughter", &fTOFSignalDaughter, "fTOFSignalDaughter/D");
  gTree->Branch("fSigmaYXDaughter", &fSigmaYXDaughter, "fSigmaYXDaughter/D");
  gTree->Branch("fSigmaXYZDaughter", &fSigmaXYZDaughter, "fSigmaXYZDaughter/D");
  gTree->Branch("fSigmaZDaughter", &fSigmaZDaughter, "fSigmaZDaughter/D");
  gTree->Branch("fPtUncertDaughter", &fPtUncertDaughter, "fPtUncertDaughter/D");
  gTree->Branch("fTPCRefitDaughter", &fTPCRefitDaughter, "fTPCRefitDaughter/I");
  gTree->Branch("fITSRefitDaughter", &fITSRefitDaughter, "fITSRefitDaughter/I");
  gTree->Branch("fITSLayer1Daughter", &fITSLayer1Daughter, "fITSLayer1Daughter/I");
  gTree->Branch("fITSLayer2Daughter", &fITSLayer2Daughter, "fITSLayer2Daughter/I");
  gTree->Branch("fITSLayer3Daughter", &fITSLayer3Daughter, "fITSLayer3Daughter/I");
  gTree->Branch("fITSLayer4Daughter", &fITSLayer4Daughter, "fITSLayer4Daughter/I");
  gTree->Branch("fITSLayer5Daughter", &fITSLayer5Daughter, "fITSLayer5Daughter/I");
  gTree->Branch("fITSLayer6Daughter", &fITSLayer6Daughter, "fITSLayer6Daughter/I");

  gTree->Branch("fTrackPIDDaughter1", &fTrackPIDDaughter1, "fTrackPIDDaughter1/I");
  gTree->Branch("fLabelDaughter1KF", &fLabelDaughter1KF, "fLabelDaughter1KF/D");
  gTree->Branch("fpDaughter1", &fpDaughter1, "fpDaughter1/D");
  gTree->Branch("fptDaughter1UnProp", &fptDaughter1UnProp, "fptDaughter1UnProp/D");
  gTree->Branch("fpxDaughter1UnProp", &fpxDaughter1UnProp, "fpxDaughter1UnProp/D");
  gTree->Branch("fpyDaughter1UnProp", &fpyDaughter1UnProp, "fpyDaughter1UnProp/D");
  gTree->Branch("fpzDaughter1UnProp", &fpzDaughter1UnProp, "fpzDaughter1UnProp/D");
  gTree->Branch("fdEdxDaughter1", &fdEdxDaughter1, "fdEdxDaughter1/D");
  gTree->Branch("fdEdxSigmaDaughter1", &fdEdxSigmaDaughter1, "fdEdxSigmaDaughter1/D");
  gTree->Branch("fDcaDaughter1", &fDcaDaughter1, "fDcaDaughter1/D");
  gTree->Branch("fDcaDaughter1o", &fDcaDaughter1o, "fDcaDaughter1o/D");
  gTree->Branch("fDcazDaughter1", &fDcazDaughter1, "fDcazDaughter1/D");
  gTree->Branch("fNclsDaughter1", &fNclsDaughter1, "fNclsDaughter1/I");
  gTree->Branch("fChi2Daughter1", &fChi2Daughter1, "fChi2Daughter1/D");
  gTree->Branch("fNclsITSDaughter1", &fNclsITSDaughter1, "fNclsITSDaughter1/I");
  gTree->Branch("fEtaDaughter1", &fEtaDaughter1, "fEtaDaughter1/D");
  gTree->Branch("fPhiDaughter1", &fPhiDaughter1, "fPhiDaughter1/D");
  gTree->Branch("fGeoLengthDaughter1", &fGeoLengthDaughter1, "fGeoLengthDaughter1/D");
  gTree->Branch("fTOFSignalDaughter1", &fTOFSignalDaughter1, "fTOFSignalDaughter1/D");
  gTree->Branch("fSigmaYXDaughter1", &fSigmaYXDaughter1, "fSigmaYXDaughter1/D");
  gTree->Branch("fSigmaXYZDaughter1", &fSigmaXYZDaughter1, "fSigmaXYZDaughter1/D");
  gTree->Branch("fSigmaZDaughter1", &fSigmaZDaughter1, "fSigmaZDaughter1/D");
  gTree->Branch("fPtUncertDaughter1", &fPtUncertDaughter1, "fPtUncertDaughter1/D");
  gTree->Branch("fTPCRefitDaughter1", &fTPCRefitDaughter1, "fTPCRefitDaughter1/I");
  gTree->Branch("fITSRefitDaughter1", &fITSRefitDaughter1, "fITSRefitDaughter1/I");
  gTree->Branch("fITSLayer1Daughter1", &fITSLayer1Daughter1, "fITSLayer1Daughter1/I");
  gTree->Branch("fITSLayer2Daughter1", &fITSLayer2Daughter1, "fITSLayer2Daughter1/I");
  gTree->Branch("fITSLayer3Daughter1", &fITSLayer3Daughter1, "fITSLayer3Daughter1/I");
  gTree->Branch("fITSLayer4Daughter1", &fITSLayer4Daughter1, "fITSLayer4Daughter1/I");
  gTree->Branch("fITSLayer5Daughter1", &fITSLayer5Daughter1, "fITSLayer5Daughter1/I");
  gTree->Branch("fITSLayer6Daughter1", &fITSLayer6Daughter1, "fITSLayer6Daughter1/I");

  gTree->Branch("fTrackPIDDaughter2", &fTrackPIDDaughter2, "fTrackPIDDaughter2/I");
  gTree->Branch("fLabelDaughter2KF", &fLabelDaughter2KF, "fLabelDaughter2KF/D");
  gTree->Branch("fpDaughter2", &fpDaughter2, "fpDaughter2/D");
  gTree->Branch("fptDaughter2UnProp", &fptDaughter2UnProp, "fptDaughter2UnProp/D");
  gTree->Branch("fpxDaughter2UnProp", &fpxDaughter2UnProp, "fpxDaughter2UnProp/D");
  gTree->Branch("fpyDaughter2UnProp", &fpyDaughter2UnProp, "fpyDaughter2UnProp/D");
  gTree->Branch("fpzDaughter2UnProp", &fpzDaughter2UnProp, "fpzDaughter2UnProp/D");
  gTree->Branch("fdEdxDaughter2", &fdEdxDaughter2, "fdEdxDaughter2/D");
  gTree->Branch("fdEdxSigmaDaughter2", &fdEdxSigmaDaughter2, "fdEdxSigmaDaughter2/D");
  gTree->Branch("fDcaDaughter2", &fDcaDaughter2, "fDcaDaughter2/D");
  gTree->Branch("fDcaDaughter2o", &fDcaDaughter2o, "fDcaDaughter2o/D");
  gTree->Branch("fDcazDaughter2", &fDcazDaughter2, "fDcazDaughter2/D");
  gTree->Branch("fNclsDaughter2", &fNclsDaughter2, "fNclsDaughter2/I");
  gTree->Branch("fChi2Daughter2", &fChi2Daughter2, "fChi2Daughter2/D");
  gTree->Branch("fNclsITSDaughter2", &fNclsITSDaughter2, "fNclsITSDaughter2/I");
  gTree->Branch("fEtaDaughter2", &fEtaDaughter2, "fEtaDaughter2/D");
  gTree->Branch("fPhiDaughter2", &fPhiDaughter2, "fPhiDaughter2/D");
  gTree->Branch("fGeoLengthDaughter2", &fGeoLengthDaughter2, "fGeoLengthDaughter2/D");
  gTree->Branch("fTOFSignalDaughter2", &fTOFSignalDaughter2, "fTOFSignalDaughter2/D");
  gTree->Branch("fSigmaYXDaughter2", &fSigmaYXDaughter2, "fSigmaYXDaughter2/D");
  gTree->Branch("fSigmaXYZDaughter2", &fSigmaXYZDaughter2, "fSigmaXYZDaughter2/D");
  gTree->Branch("fSigmaZDaughter2", &fSigmaZDaughter2, "fSigmaZDaughter2/D");
  gTree->Branch("fPtUncertDaughter2", &fPtUncertDaughter2, "fPtUncertDaughter2/D");
  gTree->Branch("fTPCRefitDaughter2", &fTPCRefitDaughter2, "fTPCRefitDaughter2/I");
  gTree->Branch("fITSRefitDaughter2", &fITSRefitDaughter2, "fITSRefitDaughter2/I");
  gTree->Branch("fITSLayer1Daughter2", &fITSLayer1Daughter2, "fITSLayer1Daughter2/I");
  gTree->Branch("fITSLayer2Daughter2", &fITSLayer2Daughter2, "fITSLayer2Daughter2/I");
  gTree->Branch("fITSLayer3Daughter2", &fITSLayer3Daughter2, "fITSLayer3Daughter2/I");
  gTree->Branch("fITSLayer4Daughter2", &fITSLayer4Daughter2, "fITSLayer4Daughter2/I");
  gTree->Branch("fITSLayer5Daughter2", &fITSLayer5Daughter2, "fITSLayer5Daughter2/I");
  gTree->Branch("fITSLayer6Daughter2", &fITSLayer6Daughter2, "fITSLayer6Daughter2/I");

  gTree->Branch("fdEdxSigmaPion", &fdEdxSigmaPion, "fdEdxSigmaPion/D");
  gTree->Branch("fdEdxSigmaDeuteron", &fdEdxSigmaDeuteron, "fdEdxSigmaDeuteron/D");
  gTree->Branch("fdEdxSigmaTriton", &fdEdxSigmaTriton, "fdEdxSigmaTriton/D");
  gTree->Branch("fdEdxSigmaAlpha", &fdEdxSigmaAlpha, "fdEdxSigmaAlpha/D");
  gTree->Branch("fTrackPar", fTrackPar, "fTrackPar[7]/D");
  gTree->Branch("fTrackPar1", fTrackPar1, "fTrackPar1[7]/D");
  gTree->Branch("fTrackPar2", fTrackPar2, "fTrackPar2[7]/D");
  gTree->Branch("fCovMatrixTrack", fCovMatrixTrack, "fCovMatrixTrack[21]/D");
  gTree->Branch("fCovMatrixTrack1", fCovMatrixTrack1, "fCovMatrixTrack1[21]/D");
  gTree->Branch("fCovMatrixTrack2", fCovMatrixTrack2, "fCovMatrixTrack2[21]/D");

  // _________________________________________________ //
  // __ associated 4LHe KF tree __ //
  gTreeKF = new TTree("gTreeKF", "gTreeKF");
  gTreeKF->Branch("fEventID", &fEventID, "fEventID/I");  
  gTreeKF->Branch("fPeriod", &fPeriod, "fPeriod/I");
  gTreeKF->Branch("frunnumber", &frunnumber, "frunnumber/I");
  gTreeKF->Branch("feventclass", &feventclass, "feventclass/I");
  gTreeKF->Branch("fCentrality", &fCentrality, "fCentrality/I");
  gTreeKF->Branch("fMagneticField", &fMagneticField, "fMagneticField/I");
  gTreeKF->Branch("fTrigMB", &fTrigMB, "fTrigMB/I");
  gTreeKF->Branch("fTrigHMV0", &fTrigHMV0, "fTrigHMV0/I");
  gTreeKF->Branch("fTrigHMSPD", &fTrigHMSPD, "fTrigHMSPD/I");
  gTreeKF->Branch("fTrigHNU", &fTrigHNU, "fTrigHNU/I");
  gTreeKF->Branch("fTrigHQU", &fTrigHQU, "fTrigHQU/I");
  gTreeKF->Branch("fTrigkCentral", &fTrigkCentral, "fTrigkCentral/I");
  gTreeKF->Branch("fTrigkSemiCentral", &fTrigkSemiCentral, "fTrigkSemiCentral/I");
  gTreeKF->Branch("fisWOTrackVertex", &fisWOTrackVertex, "fisWOTrackVertex/I");
  gTreeKF->Branch("fTriggerString", (void*)fTriggerString, "fTriggerString/C");

  gTreeKF->Branch("fPrimVertexXKF", &fPrimVertexXKF, "fPrimVertexXKF/D");
  gTreeKF->Branch("fPrimVertexYKF", &fPrimVertexYKF, "fPrimVertexYKF/D");
  gTreeKF->Branch("fPrimVertexZKF", &fPrimVertexZKF, "fPrimVertexZKF/D");
  gTreeKF->Branch("fPrimVertexXErrKF", &fPrimVertexXErrKF, "fPrimVertexXErrKF/D");
  gTreeKF->Branch("fPrimVertexYErrKF", &fPrimVertexYErrKF, "fPrimVertexYErrKF/D");
  gTreeKF->Branch("fPrimVertexZErrKF", &fPrimVertexZErrKF, "fPrimVertexZErrKF/D");
  gTreeKF->Branch("fPrimVertChi2KF", &fPrimVertChi2KF, "fPrimVertChi2KF/D");
  gTreeKF->Branch("fPrimVertNDFKF", &fPrimVertNDFKF, "fPrimVertNDFKF/I");
  gTreeKF->Branch("fTertVertexXKF", &fTertVertexXKF, "fTertVertexXKF/D");
  gTreeKF->Branch("fTertVertexYKF", &fTertVertexYKF, "fTertVertexYKF/D");
  gTreeKF->Branch("fTertVertexZKF", &fTertVertexZKF, "fTertVertexZKF/D");
  gTreeKF->Branch("fTertVertexXErrKF", &fTertVertexXErrKF, "fTertVertexXErrKF/D");
  gTreeKF->Branch("fTertVertexYErrKF", &fTertVertexYErrKF, "fTertVertexYErrKF/D");
  gTreeKF->Branch("fTertVertexZErrKF", &fTertVertexZErrKF, "fTertVertexZErrKF/D");
  gTreeKF->Branch("fTertVertChi2KF", &fTertVertChi2KF, "fTertVertChi2KF/D");
  gTreeKF->Branch("fTertVertNDFKF", &fTertVertNDFKF, "fTertVertNDFKF/I");
  gTreeKF->Branch("fDCA3B1XYKF", &fDCA3B1XYKF, "fDCA3B1XYKF/D");
  gTreeKF->Branch("fDCA3B2XYKF", &fDCA3B2XYKF, "fDCA3B2XYKF/D");
  gTreeKF->Branch("fDCA3B3XYKF", &fDCA3B3XYKF, "fDCA3B3XYKF/D");
  gTreeKF->Branch("fDCA3B1ZKF", &fDCA3B1ZKF, "fDCA3B1ZKF/D");
  gTreeKF->Branch("fDCA3B2ZKF", &fDCA3B2ZKF, "fDCA3B2ZKF/D");
  gTreeKF->Branch("fDCA3B3ZKF", &fDCA3B3ZKF, "fDCA3B3ZKF/D");
  gTreeKF->Branch("fSubPAKF", &fSubPAKF, "fSubPAKF/D");

  gTreeKF->Branch("fPDGMother", &fPDGMother, "fPDGMother/I");
  gTreeKF->Branch("fChargeMother", &fChargeMother, "fChargeMother/I");
  gTreeKF->Branch("fmctruth", &fmctruth, "fmctruth/I");
  gTreeKF->Branch("fPrimary4LHe", &fPrimary4LHe, "fPrimary4LHe/I");
  gTreeKF->Branch("fDecayChannel", &fDecayChannel, "fDecayChannel/I");
  gTreeKF->Branch("fRecoMethod", &fRecoMethod, "fRecoMethod/I");
  gTreeKF->Branch("fParticleID", &fParticleID, "fParticleID/I");
  gTreeKF->Branch("fisExcited", &fisExcited, "fisExcited/I");
  
  gTreeKF->Branch("fmSubMother", &fmSubMother, "fmSubMother/D");
  gTreeKF->Branch("fmSubMotherKF", &fmSubMotherKF, "fmSubMotherKF/D");
  gTreeKF->Branch("fmSubMotherErrKF", &fmSubMotherErrKF, "fmSubMotherErrKF/D");
  gTreeKF->Branch("fESubMotherKF", &fESubMotherKF, "fESubMotherKF/D");
  gTreeKF->Branch("fESubMotherErrKF", &fESubMotherErrKF, "fESubMotherErrKF/D");
  gTreeKF->Branch("fpxSubMotherKF", &fpxSubMotherKF, "fpxSubMotherKF/D");
  gTreeKF->Branch("fpxSubMotherErrKF", &fpxSubMotherErrKF, "fpxSubMotherErrKF/D");
  gTreeKF->Branch("fpySubMotherKF", &fpySubMotherKF, "fpySubMotherKF/D");
  gTreeKF->Branch("fpySubMotherErrKF", &fpySubMotherErrKF, "fpySubMotherErrKF/D");
  gTreeKF->Branch("fpzSubMotherKF", &fpzSubMotherKF, "fpzSubMotherKF/D");
  gTreeKF->Branch("fpzSubMotherErrKF", &fpzSubMotherErrKF, "fpzSubMotherErrKF/D");
  gTreeKF->Branch("fptSubMotherKF", &fptSubMotherKF, "fptSubMotherKF/D");
  gTreeKF->Branch("fptSubMotherErrKF", &fptSubMotherErrKF, "fptSubMotherErrKF/D");
  gTreeKF->Branch("fpSubMotherKF", &fpSubMotherKF, "fpSubMotherKF/D");
  gTreeKF->Branch("fpSubMotherErrKF", &fpSubMotherErrKF, "fpSubMotherErrKF/D");
  gTreeKF->Branch("fySubMotherKF", &fySubMotherKF, "fySubMotherKF/D");
  gTreeKF->Branch("fctSubMotherKF", &fctSubMotherKF, "fctSubMotherKF/D");
  gTreeKF->Branch("fctSubMotherErrKF", &fctSubMotherErrKF, "fctSubMotherErrKF/D");

  gTreeKF->Branch("fEDaughter", &fEDaughter, "fEDaughter/D");
  gTreeKF->Branch("fptDaughter", &fptDaughter, "fptDaughter/D");
  gTreeKF->Branch("fpxDaughter", &fpxDaughter, "fpxDaughter/D");
  gTreeKF->Branch("fpyDaughter", &fpyDaughter, "fpyDaughter/D");
  gTreeKF->Branch("fpzDaughter", &fpzDaughter, "fpzDaughter/D");
  gTreeKF->Branch("fyDaughter", &fyDaughter, "fyDaughter/D");
  gTreeKF->Branch("fEDaughterKF", &fEDaughterKF, "fEDaughterKF/D");
  gTreeKF->Branch("fptDaughterKF", &fptDaughterKF, "fptDaughterKF/D");
  gTreeKF->Branch("fpxDaughterKF", &fpxDaughterKF, "fpxDaughterKF/D");
  gTreeKF->Branch("fpyDaughterKF", &fpyDaughterKF, "fpyDaughterKF/D");
  gTreeKF->Branch("fpzDaughterKF", &fpzDaughterKF, "fpzDaughterKF/D");
  gTreeKF->Branch("fyDaughterKF", &fyDaughterKF, "fyDaughterKF/D");
  gTreeKF->Branch("fEDaughter1", &fEDaughter1, "fEDaughter1/D");
  gTreeKF->Branch("fptDaughter1", &fptDaughter1, "fptDaughter1/D");
  gTreeKF->Branch("fpxDaughter1", &fpxDaughter1, "fpxDaughter1/D");
  gTreeKF->Branch("fpyDaughter1", &fpyDaughter1, "fpyDaughter1/D");
  gTreeKF->Branch("fpzDaughter1", &fpzDaughter1, "fpzDaughter1/D");
  gTreeKF->Branch("fyDaughter1", &fyDaughter1, "fyDaughter1/D");
  gTreeKF->Branch("fEDaughter1KF", &fEDaughter1KF, "fEDaughter1KF/D");
  gTreeKF->Branch("fptDaughter1KF", &fptDaughter1KF, "fptDaughter1KF/D");
  gTreeKF->Branch("fpxDaughter1KF", &fpxDaughter1KF, "fpxDaughter1KF/D");
  gTreeKF->Branch("fpyDaughter1KF", &fpyDaughter1KF, "fpyDaughter1KF/D");
  gTreeKF->Branch("fpzDaughter1KF", &fpzDaughter1KF, "fpzDaughter1KF/D");
  gTreeKF->Branch("fyDaughter1KF", &fyDaughter1KF, "fyDaughter1KF/D");
  gTreeKF->Branch("fEDaughter2", &fEDaughter2, "fEDaughter2/D");
  gTreeKF->Branch("fptDaughter2", &fptDaughter2, "fptDaughter2/D");
  gTreeKF->Branch("fpxDaughter2", &fpxDaughter2, "fpxDaughter2/D");
  gTreeKF->Branch("fpyDaughter2", &fpyDaughter2, "fpyDaughter2/D");
  gTreeKF->Branch("fpzDaughter2", &fpzDaughter2, "fpzDaughter2/D");
  gTreeKF->Branch("fyDaughter2", &fyDaughter2, "fyDaughter2/D");
  gTreeKF->Branch("fEDaughter2KF", &fEDaughter2KF, "fEDaughter2KF/D");
  gTreeKF->Branch("fptDaughter2KF", &fptDaughter2KF, "fptDaughter2KF/D");
  gTreeKF->Branch("fpxDaughter2KF", &fpxDaughter2KF, "fpxDaughter2KF/D");
  gTreeKF->Branch("fpyDaughter2KF", &fpyDaughter2KF, "fpyDaughter2KF/D");
  gTreeKF->Branch("fpzDaughter2KF", &fpzDaughter2KF, "fpzDaughter2KF/D");
  gTreeKF->Branch("fyDaughter2KF", &fyDaughter2KF, "fyDaughter2KF/D");

  gTreeKF->Branch("fDcaDaughterXYKF", &fDcaDaughterXYKF, "fDcaDaughterXYKF/D");
  gTreeKF->Branch("fDcaDaughter1XYKF", &fDcaDaughter1XYKF, "fDcaDaughter1XYKF/D");
  gTreeKF->Branch("fDcaDaughter2XYKF", &fDcaDaughter2XYKF, "fDcaDaughter2XYKF/D");
  gTreeKF->Branch("fDcaDaughterZKF", &fDcaDaughterZKF, "fDcaDaughterZKF/D");
  gTreeKF->Branch("fDcaDaughter1ZKF", &fDcaDaughter1ZKF, "fDcaDaughter1ZKF/D");
  gTreeKF->Branch("fDcaDaughter2ZKF", &fDcaDaughter2ZKF, "fDcaDaughter2ZKF/D");
  gTreeKF->Branch("fDcaSecDaughterXYKF", &fDcaSecDaughterXYKF, "fDcaSecDaughterXYKF/D");
  gTreeKF->Branch("fDcaSecDaughterZKF", &fDcaSecDaughterZKF, "fDcaSecDaughterZKF/D");
  gTreeKF->Branch("fDcaSecDaughter1XYKF", &fDcaSecDaughter1XYKF, "fDcaSecDaughter1XYKF/D");
  gTreeKF->Branch("fDcaSecDaughter1ZKF", &fDcaSecDaughter1ZKF, "fDcaSecDaughter1ZKF/D");
  gTreeKF->Branch("fDcaSecDaughter2XYKF", &fDcaSecDaughter2XYKF, "fDcaSecDaughter2XYKF/D");
  gTreeKF->Branch("fDcaSecDaughter2ZKF", &fDcaSecDaughter2ZKF, "fDcaSecDaughter2ZKF/D");
  gTreeKF->Branch("fPropDCADaughter", &fPropDCADaughter, "fPropDCADaughter/I");
  gTreeKF->Branch("fImParDaughter", &fImParDaughter, "fImParDaughter/D");
  gTreeKF->Branch("fImParzDaughter", &fImParzDaughter, "fImParzDaughter/D");
  gTreeKF->Branch("fPropDCADaughter1", &fPropDCADaughter1, "fPropDCADaughter1/I");
  gTreeKF->Branch("fImParDaughter1", &fImParDaughter1, "fImParDaughter1/D");
  gTreeKF->Branch("fImParzDaughter1", &fImParzDaughter1, "fImParzDaughter1/D");
  gTreeKF->Branch("fPropDCADaughter2", &fPropDCADaughter2, "fPropDCADaughter2/I");
  gTreeKF->Branch("fImParDaughter2", &fImParDaughter2, "fImParDaughter2/D");
  gTreeKF->Branch("fImParzDaughter2", &fImParzDaughter2, "fImParzDaughter2/D");
  gTreeKF->Branch("fDcaSecDaughter", &fDcaSecDaughter, "fDcaSecDaughter/D");
  gTreeKF->Branch("fDcaSecDaughter1", &fDcaSecDaughter1, "fDcaSecDaughter1/D");
  gTreeKF->Branch("fDcaSecDaughter2", &fDcaSecDaughter2, "fDcaSecDaughter2/D");

  gTreeKF->Branch("fTrackPIDDaughter", &fTrackPIDDaughter, "fTrackPIDDaughter/I");
  gTreeKF->Branch("fLabelDaughterKF", &fLabelDaughterKF, "fLabelDaughterKF/D");
  gTreeKF->Branch("fpDaughter", &fpDaughter, "fpDaughter/D");
  gTreeKF->Branch("fptDaughterUnProp", &fptDaughterUnProp, "fptDaughterUnProp/D");
  gTreeKF->Branch("fpxDaughterUnProp", &fpxDaughterUnProp, "fpxDaughterUnProp/D");
  gTreeKF->Branch("fpyDaughterUnProp", &fpyDaughterUnProp, "fpyDaughterUnProp/D");
  gTreeKF->Branch("fpzDaughterUnProp", &fpzDaughterUnProp, "fpzDaughterUnProp/D");
  gTreeKF->Branch("fdEdxDaughter", &fdEdxDaughter, "fdEdxDaughter/D");
  gTreeKF->Branch("fdEdxSigmaDaughter", &fdEdxSigmaDaughter, "fdEdxSigmaDaughter/D");
  gTreeKF->Branch("fDcaDaughter", &fDcaDaughter, "fDcaDaughter/D");
  gTreeKF->Branch("fDcaDaughtero", &fDcaDaughtero, "fDcaDaughtero/D");
  gTreeKF->Branch("fDcazDaughter", &fDcazDaughter, "fDcazDaughter/D");
  gTreeKF->Branch("fNclsDaughter", &fNclsDaughter, "fNclsDaughter/I");
  gTreeKF->Branch("fChi2Daughter", &fChi2Daughter, "fChi2Daughter/D");
  gTreeKF->Branch("fNclsITSDaughter", &fNclsITSDaughter, "fNclsITSDaughter/I");
  gTreeKF->Branch("fEtaDaughter", &fEtaDaughter, "fEtaDaughter/D");
  gTreeKF->Branch("fPhiDaughter", &fPhiDaughter, "fPhiDaughter/D");
  gTreeKF->Branch("fGeoLengthDaughter", &fGeoLengthDaughter, "fGeoLengthDaughter/D");
  gTreeKF->Branch("fTOFSignalDaughter", &fTOFSignalDaughter, "fTOFSignalDaughter/D");
  gTreeKF->Branch("fSigmaYXDaughter", &fSigmaYXDaughter, "fSigmaYXDaughter/D");
  gTreeKF->Branch("fSigmaXYZDaughter", &fSigmaXYZDaughter, "fSigmaXYZDaughter/D");
  gTreeKF->Branch("fSigmaZDaughter", &fSigmaZDaughter, "fSigmaZDaughter/D");
  gTreeKF->Branch("fPtUncertDaughter", &fPtUncertDaughter, "fPtUncertDaughter/D");
  gTreeKF->Branch("fTPCRefitDaughter", &fTPCRefitDaughter, "fTPCRefitDaughter/I");
  gTreeKF->Branch("fITSRefitDaughter", &fITSRefitDaughter, "fITSRefitDaughter/I");
  gTreeKF->Branch("fITSLayer1Daughter", &fITSLayer1Daughter, "fITSLayer1Daughter/I");
  gTreeKF->Branch("fITSLayer2Daughter", &fITSLayer2Daughter, "fITSLayer2Daughter/I");
  gTreeKF->Branch("fITSLayer3Daughter", &fITSLayer3Daughter, "fITSLayer3Daughter/I");
  gTreeKF->Branch("fITSLayer4Daughter", &fITSLayer4Daughter, "fITSLayer4Daughter/I");
  gTreeKF->Branch("fITSLayer5Daughter", &fITSLayer5Daughter, "fITSLayer5Daughter/I");
  gTreeKF->Branch("fITSLayer6Daughter", &fITSLayer6Daughter, "fITSLayer6Daughter/I");

  gTreeKF->Branch("fTrackPIDDaughter1", &fTrackPIDDaughter1, "fTrackPIDDaughter1/I");
  gTreeKF->Branch("fLabelDaughter1KF", &fLabelDaughter1KF, "fLabelDaughter1KF/D");
  gTreeKF->Branch("fpDaughter1", &fpDaughter1, "fpDaughter1/D");
  gTreeKF->Branch("fptDaughter1UnProp", &fptDaughter1UnProp, "fptDaughter1UnProp/D");
  gTreeKF->Branch("fpxDaughter1UnProp", &fpxDaughter1UnProp, "fpxDaughter1UnProp/D");
  gTreeKF->Branch("fpyDaughter1UnProp", &fpyDaughter1UnProp, "fpyDaughter1UnProp/D");
  gTreeKF->Branch("fpzDaughter1UnProp", &fpzDaughter1UnProp, "fpzDaughter1UnProp/D");
  gTreeKF->Branch("fdEdxDaughter1", &fdEdxDaughter1, "fdEdxDaughter1/D");
  gTreeKF->Branch("fdEdxSigmaDaughter1", &fdEdxSigmaDaughter1, "fdEdxSigmaDaughter1/D");
  gTreeKF->Branch("fDcaDaughter1", &fDcaDaughter1, "fDcaDaughter1/D");
  gTreeKF->Branch("fDcaDaughter1o", &fDcaDaughter1o, "fDcaDaughter1o/D");
  gTreeKF->Branch("fDcazDaughter1", &fDcazDaughter1, "fDcazDaughter1/D");
  gTreeKF->Branch("fNclsDaughter1", &fNclsDaughter1, "fNclsDaughter1/I");
  gTreeKF->Branch("fChi2Daughter1", &fChi2Daughter1, "fChi2Daughter1/D");
  gTreeKF->Branch("fNclsITSDaughter1", &fNclsITSDaughter1, "fNclsITSDaughter1/I");
  gTreeKF->Branch("fEtaDaughter1", &fEtaDaughter1, "fEtaDaughter1/D");
  gTreeKF->Branch("fPhiDaughter1", &fPhiDaughter1, "fPhiDaughter1/D");
  gTreeKF->Branch("fGeoLengthDaughter1", &fGeoLengthDaughter1, "fGeoLengthDaughter1/D");
  gTreeKF->Branch("fTOFSignalDaughter1", &fTOFSignalDaughter1, "fTOFSignalDaughter1/D");
  gTreeKF->Branch("fSigmaYXDaughter1", &fSigmaYXDaughter1, "fSigmaYXDaughter1/D");
  gTreeKF->Branch("fSigmaXYZDaughter1", &fSigmaXYZDaughter1, "fSigmaXYZDaughter1/D");
  gTreeKF->Branch("fSigmaZDaughter1", &fSigmaZDaughter1, "fSigmaZDaughter1/D");
  gTreeKF->Branch("fPtUncertDaughter1", &fPtUncertDaughter1, "fPtUncertDaughter1/D");
  gTreeKF->Branch("fTPCRefitDaughter1", &fTPCRefitDaughter1, "fTPCRefitDaughter1/I");
  gTreeKF->Branch("fITSRefitDaughter1", &fITSRefitDaughter1, "fITSRefitDaughter1/I");
  gTreeKF->Branch("fITSLayer1Daughter1", &fITSLayer1Daughter1, "fITSLayer1Daughter1/I");
  gTreeKF->Branch("fITSLayer2Daughter1", &fITSLayer2Daughter1, "fITSLayer2Daughter1/I");
  gTreeKF->Branch("fITSLayer3Daughter1", &fITSLayer3Daughter1, "fITSLayer3Daughter1/I");
  gTreeKF->Branch("fITSLayer4Daughter1", &fITSLayer4Daughter1, "fITSLayer4Daughter1/I");
  gTreeKF->Branch("fITSLayer5Daughter1", &fITSLayer5Daughter1, "fITSLayer5Daughter1/I");
  gTreeKF->Branch("fITSLayer6Daughter1", &fITSLayer6Daughter1, "fITSLayer6Daughter1/I");

  gTreeKF->Branch("fTrackPIDDaughter2", &fTrackPIDDaughter2, "fTrackPIDDaughter2/I");
  gTreeKF->Branch("fLabelDaughter2KF", &fLabelDaughter2KF, "fLabelDaughter2KF/D");
  gTreeKF->Branch("fpDaughter2", &fpDaughter2, "fpDaughter2/D");
  gTreeKF->Branch("fptDaughter2UnProp", &fptDaughter2UnProp, "fptDaughter2UnProp/D");
  gTreeKF->Branch("fpxDaughter2UnProp", &fpxDaughter2UnProp, "fpxDaughter2UnProp/D");
  gTreeKF->Branch("fpyDaughter2UnProp", &fpyDaughter2UnProp, "fpyDaughter2UnProp/D");
  gTreeKF->Branch("fpzDaughter2UnProp", &fpzDaughter2UnProp, "fpzDaughter2UnProp/D");
  gTreeKF->Branch("fdEdxDaughter2", &fdEdxDaughter2, "fdEdxDaughter2/D");
  gTreeKF->Branch("fdEdxSigmaDaughter2", &fdEdxSigmaDaughter2, "fdEdxSigmaDaughter2/D");
  gTreeKF->Branch("fDcaDaughter2", &fDcaDaughter2, "fDcaDaughter2/D");
  gTreeKF->Branch("fDcaDaughter2o", &fDcaDaughter2o, "fDcaDaughter2o/D");
  gTreeKF->Branch("fDcazDaughter2", &fDcazDaughter2, "fDcazDaughter2/D");
  gTreeKF->Branch("fNclsDaughter2", &fNclsDaughter2, "fNclsDaughter2/I");
  gTreeKF->Branch("fChi2Daughter2", &fChi2Daughter2, "fChi2Daughter2/D");
  gTreeKF->Branch("fNclsITSDaughter2", &fNclsITSDaughter2, "fNclsITSDaughter2/I");
  gTreeKF->Branch("fEtaDaughter2", &fEtaDaughter2, "fEtaDaughter2/D");
  gTreeKF->Branch("fPhiDaughter2", &fPhiDaughter2, "fPhiDaughter2/D");
  gTreeKF->Branch("fGeoLengthDaughter2", &fGeoLengthDaughter2, "fGeoLengthDaughter2/D");
  gTreeKF->Branch("fTOFSignalDaughter2", &fTOFSignalDaughter2, "fTOFSignalDaughter2/D");
  gTreeKF->Branch("fSigmaYXDaughter2", &fSigmaYXDaughter2, "fSigmaYXDaughter2/D");
  gTreeKF->Branch("fSigmaXYZDaughter2", &fSigmaXYZDaughter2, "fSigmaXYZDaughter2/D");
  gTreeKF->Branch("fSigmaZDaughter2", &fSigmaZDaughter2, "fSigmaZDaughter2/D");
  gTreeKF->Branch("fPtUncertDaughter2", &fPtUncertDaughter2, "fPtUncertDaughter2/D");
  gTreeKF->Branch("fTPCRefitDaughter2", &fTPCRefitDaughter2, "fTPCRefitDaughter2/I");
  gTreeKF->Branch("fITSRefitDaughter2", &fITSRefitDaughter2, "fITSRefitDaughter2/I");
  gTreeKF->Branch("fITSLayer1Daughter2", &fITSLayer1Daughter2, "fITSLayer1Daughter2/I");
  gTreeKF->Branch("fITSLayer2Daughter2", &fITSLayer2Daughter2, "fITSLayer2Daughter2/I");
  gTreeKF->Branch("fITSLayer3Daughter2", &fITSLayer3Daughter2, "fITSLayer3Daughter2/I");
  gTreeKF->Branch("fITSLayer4Daughter2", &fITSLayer4Daughter2, "fITSLayer4Daughter2/I");
  gTreeKF->Branch("fITSLayer5Daughter2", &fITSLayer5Daughter2, "fITSLayer5Daughter2/I");
  gTreeKF->Branch("fITSLayer6Daughter2", &fITSLayer6Daughter2, "fITSLayer6Daughter2/I");

  gTreeKF->Branch("fdEdxSigmaPion", &fdEdxSigmaPion, "fdEdxSigmaPion/D");
  gTreeKF->Branch("fdEdxSigmaDeuteron", &fdEdxSigmaDeuteron, "fdEdxSigmaDeuteron/D");
  gTreeKF->Branch("fdEdxSigmaTriton", &fdEdxSigmaTriton, "fdEdxSigmaTriton/D");
  gTreeKF->Branch("fdEdxSigmaAlpha", &fdEdxSigmaAlpha, "fdEdxSigmaAlpha/D");
  gTreeKF->Branch("fTrackPar", fTrackPar, "fTrackPar[7]/D");
  gTreeKF->Branch("fTrackPar1", fTrackPar1, "fTrackPar1[7]/D");
  gTreeKF->Branch("fTrackPar2", fTrackPar2, "fTrackPar2[7]/D");
  gTreeKF->Branch("fCovMatrixTrack", fCovMatrixTrack, "fCovMatrixTrack[21]/D");
  gTreeKF->Branch("fCovMatrixTrack1", fCovMatrixTrack1, "fCovMatrixTrack1[21]/D");
  gTreeKF->Branch("fCovMatrixTrack2", fCovMatrixTrack2, "fCovMatrixTrack2[21]/D");

  // _________________________________________________ //
  // __ associated 4LH tree __ //
  hTree = new TTree("hTree", "hTree");
  hTree->Branch("fEventID", &fEventID, "fEventID/I");  
  hTree->Branch("fPeriod", &fPeriod, "fPeriod/I");
  hTree->Branch("frunnumber", &frunnumber, "frunnumber/I");
  hTree->Branch("feventclass", &feventclass, "feventclass/I");
  hTree->Branch("fCentrality", &fCentrality, "fCentrality/I");
  hTree->Branch("fMagneticField", &fMagneticField, "fMagneticField/I");
  hTree->Branch("fisWOTrackVertex", &fisWOTrackVertex, "fisWOTrackVertex/I");
  hTree->Branch("fmctruth", &fmctruth, "fmctruth/I");
  hTree->Branch("fTrigMB", &fTrigMB, "fTrigMB/I");
  hTree->Branch("fTrigHMV0", &fTrigHMV0, "fTrigHMV0/I");
  hTree->Branch("fTrigHMSPD", &fTrigHMSPD, "fTrigHMSPD/I");
  hTree->Branch("fTrigHNU", &fTrigHNU, "fTrigHNU/I");
  hTree->Branch("fTrigHQU", &fTrigHQU, "fTrigHQU/I");
  hTree->Branch("fTrigkCentral", &fTrigkCentral, "fTrigkCentral/I");
  hTree->Branch("fTrigkSemiCentral", &fTrigkSemiCentral, "fTrigkSemiCentral/I");

  hTree->Branch("fParticleID", &fParticleID, "fParticleID/I");
  hTree->Branch("fPDGMother", &fPDGMother, "fPDGMother/I");
  hTree->Branch("fChargeMother", &fChargeMother, "fChargeMother/I");
  hTree->Branch("fDecayChannel", &fDecayChannel, "fDecayChannel/I");
  hTree->Branch("fRecoMethod", &fRecoMethod, "fRecoMethod/I");
  hTree->Branch("fPrimary4LHe", &fPrimary4LHe, "fPrimary4LHe/I");
  hTree->Branch("fisExcited", &fisExcited, "fisExcited/I");
  
  hTree->Branch("fPrimVertexX", &fPrimVertexX, "fPrimVertexX/D");
  hTree->Branch("fPrimVertexY", &fPrimVertexY, "fPrimVertexY/D");
  hTree->Branch("fPrimVertexZ", &fPrimVertexZ, "fPrimVertexZ/D");
  hTree->Branch("fTertVertexX", &fTertVertexX, "fTertVertexX/D");
  hTree->Branch("fTertVertexY", &fTertVertexY, "fTertVertexY/D");
  hTree->Branch("fTertVertexZ", &fTertVertexZ, "fTertVertexZ/D");
  hTree->Branch("fDCA3B1", &fDCA3B1, "fDCA3B1/D");
  hTree->Branch("fSubPA", &fSubPA, "fSubPA/D");
 
  hTree->Branch("fmSubMother", &fmSubMother, "fmSubMother/D");
  hTree->Branch("fESubMother", &fESubMother, "fESubMother/D");
  hTree->Branch("fpxSubMother", &fpxSubMother, "fpxSubMother/D");
  hTree->Branch("fpySubMother", &fpySubMother, "fpySubMother/D");
  hTree->Branch("fpzSubMother", &fpzSubMother, "fpzSubMother/D");
  hTree->Branch("fptSubMother", &fptSubMother, "fptSubMother/D");
  hTree->Branch("fpSubMother", &fpSubMother, "fpSubMother/D");
  hTree->Branch("fySubMother", &fySubMother, "fySubMother/D");
  hTree->Branch("fctSubMother", &fctSubMother, "fctSubMother/D");

  hTree->Branch("fEDaughter", &fEDaughter, "fEDaughter/D");
  hTree->Branch("fptDaughter", &fptDaughter, "fptDaughter/D");
  hTree->Branch("fpxDaughter", &fpxDaughter, "fpxDaughter/D");
  hTree->Branch("fpyDaughter", &fpyDaughter, "fpyDaughter/D");
  hTree->Branch("fpzDaughter", &fpzDaughter, "fpzDaughter/D");
  hTree->Branch("fyDaughter", &fyDaughter, "fyDaughter/D");
  hTree->Branch("fEDaughter2", &fEDaughter2, "fEDaughter2/D");
  hTree->Branch("fptDaughter2", &fptDaughter2, "fptDaughter2/D");
  hTree->Branch("fpxDaughter2", &fpxDaughter2, "fpxDaughter2/D");
  hTree->Branch("fpyDaughter2", &fpyDaughter2, "fpyDaughter2/D");
  hTree->Branch("fpzDaughter2", &fpzDaughter2, "fpzDaughter2/D");
  hTree->Branch("fyDaughter2", &fyDaughter2, "fyDaughter2/D");

  hTree->Branch("fPropDCADaughter", &fPropDCADaughter, "fPropDCADaughter/I");
  hTree->Branch("fImParDaughter", &fImParDaughter, "fImParDaughter/D");
  hTree->Branch("fImParzDaughter", &fImParzDaughter, "fImParzDaughter/D");
  hTree->Branch("fPropDCADaughter2", &fPropDCADaughter2, "fPropDCADaughter2/I");
  hTree->Branch("fImParDaughter2", &fImParDaughter2, "fImParDaughter2/D");
  hTree->Branch("fImParzDaughter2", &fImParzDaughter2, "fImParzDaughter2/D");
  hTree->Branch("fDcaSecDaughter", &fDcaSecDaughter, "fDcaSecDaughter/D");
  hTree->Branch("fDcaSecDaughter2", &fDcaSecDaughter2, "fDcaSecDaughter2/D");  

  hTree->Branch("fTrackPIDDaughter", &fTrackPIDDaughter, "fTrackPIDDaughter/I");
  hTree->Branch("fLabelDaughterKF", &fLabelDaughterKF, "fLabelDaughterKF/D");
  hTree->Branch("fpDaughter", &fpDaughter, "fpDaughter/D");
  hTree->Branch("fptDaughterUnProp", &fptDaughterUnProp, "fptDaughterUnProp/D");
  hTree->Branch("fpxDaughterUnProp", &fpxDaughterUnProp, "fpxDaughterUnProp/D");
  hTree->Branch("fpyDaughterUnProp", &fpyDaughterUnProp, "fpyDaughterUnProp/D");
  hTree->Branch("fpzDaughterUnProp", &fpzDaughterUnProp, "fpzDaughterUnProp/D");
  hTree->Branch("fdEdxDaughter", &fdEdxDaughter, "fdEdxDaughter/D");
  hTree->Branch("fdEdxSigmaDaughter", &fdEdxSigmaDaughter, "fdEdxSigmaDaughter/D");
  hTree->Branch("fDcaDaughter", &fDcaDaughter, "fDcaDaughter/D");
  hTree->Branch("fDcaDaughtero", &fDcaDaughtero, "fDcaDaughtero/D");
  hTree->Branch("fDcazDaughter", &fDcazDaughter, "fDcazDaughter/D");
  hTree->Branch("fNclsDaughter", &fNclsDaughter, "fNclsDaughter/I");
  hTree->Branch("fChi2Daughter", &fChi2Daughter, "fChi2Daughter/D");
  hTree->Branch("fNclsITSDaughter", &fNclsITSDaughter, "fNclsITSDaughter/I");
  hTree->Branch("fEtaDaughter", &fEtaDaughter, "fEtaDaughter/D");
  hTree->Branch("fPhiDaughter", &fPhiDaughter, "fPhiDaughter/D");
  hTree->Branch("fGeoLengthDaughter", &fGeoLengthDaughter, "fGeoLengthDaughter/D");
  hTree->Branch("fTOFSignalDaughter", &fTOFSignalDaughter, "fTOFSignalDaughter/D");
  hTree->Branch("fSigmaYXDaughter", &fSigmaYXDaughter, "fSigmaYXDaughter/D");
  hTree->Branch("fSigmaXYZDaughter", &fSigmaXYZDaughter, "fSigmaXYZDaughter/D");
  hTree->Branch("fSigmaZDaughter", &fSigmaZDaughter, "fSigmaZDaughter/D");
  hTree->Branch("fPtUncertDaughter", &fPtUncertDaughter, "fPtUncertDaughter/D");
  hTree->Branch("fTPCRefitDaughter", &fTPCRefitDaughter, "fTPCRefitDaughter/I");
  hTree->Branch("fITSRefitDaughter", &fITSRefitDaughter, "fITSRefitDaughter/I");
  hTree->Branch("fITSLayer1Daughter", &fITSLayer1Daughter, "fITSLayer1Daughter/I");
  hTree->Branch("fITSLayer2Daughter", &fITSLayer2Daughter, "fITSLayer2Daughter/I");
  hTree->Branch("fITSLayer3Daughter", &fITSLayer3Daughter, "fITSLayer3Daughter/I");
  hTree->Branch("fITSLayer4Daughter", &fITSLayer4Daughter, "fITSLayer4Daughter/I");
  hTree->Branch("fITSLayer5Daughter", &fITSLayer5Daughter, "fITSLayer5Daughter/I");
  hTree->Branch("fITSLayer6Daughter", &fITSLayer6Daughter, "fITSLayer6Daughter/I");

  hTree->Branch("fTrackPIDDaughter2", &fTrackPIDDaughter2, "fTrackPIDDaughter2/I");
  hTree->Branch("fLabelDaughter2KF", &fLabelDaughter2KF, "fLabelDaughter2KF/D");
  hTree->Branch("fpDaughter2", &fpDaughter2, "fpDaughter2/D");
  hTree->Branch("fptDaughter2UnProp", &fptDaughter2UnProp, "fptDaughter2UnProp/D");
  hTree->Branch("fpxDaughter2UnProp", &fpxDaughter2UnProp, "fpxDaughter2UnProp/D");
  hTree->Branch("fpyDaughter2UnProp", &fpyDaughter2UnProp, "fpyDaughter2UnProp/D");
  hTree->Branch("fpzDaughter2UnProp", &fpzDaughter2UnProp, "fpzDaughter2UnProp/D");
  hTree->Branch("fdEdxDaughter2", &fdEdxDaughter2, "fdEdxDaughter2/D");
  hTree->Branch("fdEdxSigmaDaughter2", &fdEdxSigmaDaughter2, "fdEdxSigmaDaughter2/D");
  hTree->Branch("fDcaDaughter2", &fDcaDaughter2, "fDcaDaughter2/D");
  hTree->Branch("fDcaDaughter2o", &fDcaDaughter2o, "fDcaDaughter2o/D");
  hTree->Branch("fDcazDaughter2", &fDcazDaughter2, "fDcazDaughter2/D");
  hTree->Branch("fNclsDaughter2", &fNclsDaughter2, "fNclsDaughter2/I");
  hTree->Branch("fChi2Daughter2", &fChi2Daughter2, "fChi2Daughter2/D");
  hTree->Branch("fNclsITSDaughter2", &fNclsITSDaughter2, "fNclsITSDaughter2/I");
  hTree->Branch("fEtaDaughter2", &fEtaDaughter2, "fEtaDaughter2/D");
  hTree->Branch("fPhiDaughter2", &fPhiDaughter2, "fPhiDaughter2/D");
  hTree->Branch("fGeoLengthDaughter2", &fGeoLengthDaughter2, "fGeoLengthDaughter2/D");
  hTree->Branch("fTOFSignalDaughter2", &fTOFSignalDaughter2, "fTOFSignalDaughter2/D");
  hTree->Branch("fSigmaYXDaughter2", &fSigmaYXDaughter2, "fSigmaYXDaughter2/D");
  hTree->Branch("fSigmaXYZDaughter2", &fSigmaXYZDaughter2, "fSigmaXYZDaughter2/D");
  hTree->Branch("fSigmaZDaughter2", &fSigmaZDaughter2, "fSigmaZDaughter2/D");
  hTree->Branch("fPtUncertDaughter2", &fPtUncertDaughter2, "fPtUncertDaughter2/D");
  hTree->Branch("fTPCRefitDaughter2", &fTPCRefitDaughter2, "fTPCRefitDaughter2/I");
  hTree->Branch("fITSRefitDaughter2", &fITSRefitDaughter2, "fITSRefitDaughter2/I");
  hTree->Branch("fITSLayer1Daughter2", &fITSLayer1Daughter2, "fITSLayer1Daughter2/I");
  hTree->Branch("fITSLayer2Daughter2", &fITSLayer2Daughter2, "fITSLayer2Daughter2/I");
  hTree->Branch("fITSLayer3Daughter2", &fITSLayer3Daughter2, "fITSLayer3Daughter2/I");
  hTree->Branch("fITSLayer4Daughter2", &fITSLayer4Daughter2, "fITSLayer4Daughter2/I");
  hTree->Branch("fITSLayer5Daughter2", &fITSLayer5Daughter2, "fITSLayer5Daughter2/I");
  hTree->Branch("fITSLayer6Daughter2", &fITSLayer6Daughter2, "fITSLayer6Daughter2/I");

  hTree->Branch("fdEdxSigmaPion", &fdEdxSigmaPion, "fdEdxSigmaPion/D");
  hTree->Branch("fdEdxSigmaDeuteron", &fdEdxSigmaDeuteron, "fdEdxSigmaDeuteron/D");
  hTree->Branch("fdEdxSigmaTriton", &fdEdxSigmaTriton, "fdEdxSigmaTriton/D");
  hTree->Branch("fdEdxSigmaAlpha", &fdEdxSigmaAlpha, "fdEdxSigmaAlpha/D");
  hTree->Branch("fTrackPar", fTrackPar, "fTrackPar[7]/D");
  hTree->Branch("fTrackPar2", fTrackPar2, "fTrackPar2[7]/D");
  hTree->Branch("fCovMatrixTrack", fCovMatrixTrack, "fCovMatrixTrack[21]/D");
  hTree->Branch("fCovMatrixTrack2", fCovMatrixTrack2, "fCovMatrixTrack2[21]/D");

  // _________________________________________________ //
  // __ associated 4LH KF tree __ //
  hTreeKF = new TTree("hTreeKF", "hTreeKF");    
  hTreeKF->Branch("fPeriod", &fPeriod, "fPeriod/I");
  hTreeKF->Branch("frunnumber", &frunnumber, "frunnumber/I");
  hTreeKF->Branch("fEventID", &fEventID, "fEventID/I");
  hTreeKF->Branch("feventclass", &feventclass, "feventclass/I");
  hTreeKF->Branch("fCentrality", &fCentrality, "fCentrality/I");
  hTreeKF->Branch("fMagneticField", &fMagneticField, "fMagneticField/I");
  hTreeKF->Branch("fisWOTrackVertex", &fisWOTrackVertex, "fisWOTrackVertex/I");
  hTreeKF->Branch("fTrigMB", &fTrigMB, "fTrigMB/I");
  hTreeKF->Branch("fTrigHMV0", &fTrigHMV0, "fTrigHMV0/I");
  hTreeKF->Branch("fTrigHMSPD", &fTrigHMSPD, "fTrigHMSPD/I");
  hTreeKF->Branch("fTrigHNU", &fTrigHNU, "fTrigHNU/I");
  hTreeKF->Branch("fTrigHQU", &fTrigHQU, "fTrigHQU/I");
  hTreeKF->Branch("fTrigkCentral", &fTrigkCentral, "fTrigkCentral/I");
  hTreeKF->Branch("fTrigkSemiCentral", &fTrigkSemiCentral, "fTrigkSemiCentral/I");

  hTreeKF->Branch("fPDGMother", &fPDGMother, "fPDGMother/I");
  hTreeKF->Branch("fChargeMother", &fChargeMother, "fChargeMother/I");
  hTreeKF->Branch("fmctruth", &fmctruth, "fmctruth/I");
  hTreeKF->Branch("fPrimary4LHe", &fPrimary4LHe, "fPrimary4LHe/I");
  hTreeKF->Branch("fDecayChannel", &fDecayChannel, "fDecayChannel/I");
  hTreeKF->Branch("fRecoMethod", &fRecoMethod, "fRecoMethod/I");
  hTreeKF->Branch("fParticleID", &fParticleID, "fParticleID/I");
  hTreeKF->Branch("fisExcited", &fisExcited, "fisExcited/I");

  hTreeKF->Branch("fPrimVertexXKF", &fPrimVertexXKF, "fPrimVertexXKF/D");
  hTreeKF->Branch("fPrimVertexYKF", &fPrimVertexYKF, "fPrimVertexYKF/D");
  hTreeKF->Branch("fPrimVertexZKF", &fPrimVertexZKF, "fPrimVertexZKF/D");
  hTreeKF->Branch("fPrimVertexXErrKF", &fPrimVertexXErrKF, "fPrimVertexXErrKF/D");
  hTreeKF->Branch("fPrimVertexYErrKF", &fPrimVertexYErrKF, "fPrimVertexYErrKF/D");
  hTreeKF->Branch("fPrimVertexZErrKF", &fPrimVertexZErrKF, "fPrimVertexZErrKF/D");
  hTreeKF->Branch("fPrimVertChi2KF", &fPrimVertChi2KF, "fPrimVertChi2KF/D");
  hTreeKF->Branch("fPrimVertNDFKF", &fPrimVertNDFKF, "fPrimVertNDFKF/I");
  hTreeKF->Branch("fTertVertexXKF", &fTertVertexXKF, "fTertVertexXKF/D");
  hTreeKF->Branch("fTertVertexYKF", &fTertVertexYKF, "fTertVertexYKF/D");
  hTreeKF->Branch("fTertVertexZKF", &fTertVertexZKF, "fTertVertexZKF/D");
  hTreeKF->Branch("fTertVertexXErrKF", &fTertVertexXErrKF, "fTertVertexXErrKF/D");
  hTreeKF->Branch("fTertVertexYErrKF", &fTertVertexYErrKF, "fTertVertexYErrKF/D");
  hTreeKF->Branch("fTertVertexZErrKF", &fTertVertexZErrKF, "fTertVertexZErrKF/D");
  hTreeKF->Branch("fTertVertChi2KF", &fTertVertChi2KF, "fTertVertChi2KF/D");
  hTreeKF->Branch("fTertVertNDFKF", &fTertVertNDFKF, "fTertVertNDFKF/I");
  hTreeKF->Branch("fDCA3B1XYKF", &fDCA3B1XYKF, "fDCA3B1XYKF/D");
  hTreeKF->Branch("fDCA3B1ZKF", &fDCA3B1ZKF, "fDCA3B1ZKF/D");  
  hTreeKF->Branch("fSubPAKF", &fSubPAKF, "fSubPAKF/D");
 
  hTreeKF->Branch("fmSubMother", &fmSubMother, "fmSubMother/D");
  hTreeKF->Branch("fmSubMotherKF", &fmSubMotherKF, "fmSubMotherKF/D");
  hTreeKF->Branch("fmSubMotherErrKF", &fmSubMotherErrKF, "fmSubMotherErrKF/D");
  hTreeKF->Branch("fESubMotherKF", &fESubMotherKF, "fESubMotherKF/D");
  hTreeKF->Branch("fESubMotherErrKF", &fESubMotherErrKF, "fESubMotherErrKF/D");
  hTreeKF->Branch("fpxSubMotherKF", &fpxSubMotherKF, "fpxSubMotherKF/D");
  hTreeKF->Branch("fpxSubMotherErrKF", &fpxSubMotherErrKF, "fpxSubMotherErrKF/D");
  hTreeKF->Branch("fpySubMotherKF", &fpySubMotherKF, "fpySubMotherKF/D");
  hTreeKF->Branch("fpySubMotherErrKF", &fpySubMotherErrKF, "fpySubMotherErrKF/D");
  hTreeKF->Branch("fpzSubMotherKF", &fpzSubMotherKF, "fpzSubMotherKF/D");
  hTreeKF->Branch("fpzSubMotherErrKF", &fpzSubMotherErrKF, "fpzSubMotherErrKF/D");
  hTreeKF->Branch("fptSubMotherKF", &fptSubMotherKF, "fptSubMotherKF/D");
  hTreeKF->Branch("fptSubMotherErrKF", &fptSubMotherErrKF, "fptSubMotherErrKF/D");
  hTreeKF->Branch("fpSubMotherKF", &fpSubMotherKF, "fpSubMotherKF/D");
  hTreeKF->Branch("fpSubMotherErrKF", &fpSubMotherErrKF, "fpSubMotherErrKF/D");
  hTreeKF->Branch("fySubMotherKF", &fySubMotherKF, "fySubMotherKF/D");
  hTreeKF->Branch("fctSubMotherKF", &fctSubMotherKF, "fctSubMotherKF/D");
  hTreeKF->Branch("fctSubMotherErrKF", &fctSubMotherErrKF, "fctSubMotherErrKF/D");

  hTreeKF->Branch("fEDaughter", &fEDaughter, "fEDaughter/D");
  hTreeKF->Branch("fptDaughter", &fptDaughter, "fptDaughter/D");
  hTreeKF->Branch("fpxDaughter", &fpxDaughter, "fpxDaughter/D");
  hTreeKF->Branch("fpyDaughter", &fpyDaughter, "fpyDaughter/D");
  hTreeKF->Branch("fpzDaughter", &fpzDaughter, "fpzDaughter/D");
  hTreeKF->Branch("fyDaughter", &fyDaughter, "fyDaughter/D");
  hTreeKF->Branch("fEDaughterKF", &fEDaughterKF, "fEDaughterKF/D");
  hTreeKF->Branch("fptDaughterKF", &fptDaughterKF, "fptDaughterKF/D");
  hTreeKF->Branch("fpxDaughterKF", &fpxDaughterKF, "fpxDaughterKF/D");
  hTreeKF->Branch("fpyDaughterKF", &fpyDaughterKF, "fpyDaughterKF/D");
  hTreeKF->Branch("fpzDaughterKF", &fpzDaughterKF, "fpzDaughterKF/D");
  hTreeKF->Branch("fyDaughterKF", &fyDaughterKF, "fyDaughterKF/D");
  hTreeKF->Branch("fEDaughter2", &fEDaughter2, "fEDaughter2/D");
  hTreeKF->Branch("fptDaughter2", &fptDaughter2, "fptDaughter2/D");
  hTreeKF->Branch("fpxDaughter2", &fpxDaughter2, "fpxDaughter2/D");
  hTreeKF->Branch("fpyDaughter2", &fpyDaughter2, "fpyDaughter2/D");
  hTreeKF->Branch("fpzDaughter2", &fpzDaughter2, "fpzDaughter2/D");
  hTreeKF->Branch("fyDaughter2", &fyDaughter2, "fyDaughter2/D");
  hTreeKF->Branch("fEDaughter2KF", &fEDaughter2KF, "fEDaughter2KF/D");
  hTreeKF->Branch("fptDaughter2KF", &fptDaughter2KF, "fptDaughter2KF/D");
  hTreeKF->Branch("fpxDaughter2KF", &fpxDaughter2KF, "fpxDaughter2KF/D");
  hTreeKF->Branch("fpyDaughter2KF", &fpyDaughter2KF, "fpyDaughter2KF/D");
  hTreeKF->Branch("fpzDaughter2KF", &fpzDaughter2KF, "fpzDaughter2KF/D");
  hTreeKF->Branch("fyDaughter2KF", &fyDaughter2KF, "fyDaughter2KF/D");

  hTreeKF->Branch("fDcaDaughterXYKF", &fDcaDaughterXYKF, "fDcaDaughterXYKF/D");
  hTreeKF->Branch("fDcaDaughter2XYKF", &fDcaDaughter2XYKF, "fDcaDaughter2XYKF/D");
  hTreeKF->Branch("fDcaDaughterZKF", &fDcaDaughterZKF, "fDcaDaughterZKF/D");
  hTreeKF->Branch("fDcaDaughter2ZKF", &fDcaDaughter2ZKF, "fDcaDaughter2ZKF/D");
  hTreeKF->Branch("fDcaSecDaughterXYKF", &fDcaSecDaughterXYKF, "fDcaSecDaughterXYKF/D");
  hTreeKF->Branch("fDcaSecDaughterZKF", &fDcaSecDaughterZKF, "fDcaSecDaughterZKF/D");
  hTreeKF->Branch("fDcaSecDaughter2XYKF", &fDcaSecDaughter2XYKF, "fDcaSecDaughter2XYKF/D");
  hTreeKF->Branch("fDcaSecDaughter2ZKF", &fDcaSecDaughter2ZKF, "fDcaSecDaughter2ZKF/D");
  hTreeKF->Branch("fPropDCADaughter", &fPropDCADaughter, "fPropDCADaughter/I");
  hTreeKF->Branch("fImParDaughter", &fImParDaughter, "fImParDaughter/D");
  hTreeKF->Branch("fImParzDaughter", &fImParzDaughter, "fImParzDaughter/D");
  hTreeKF->Branch("fPropDCADaughter2", &fPropDCADaughter2, "fPropDCADaughter2/I");
  hTreeKF->Branch("fImParDaughter2", &fImParDaughter2, "fImParDaughter2/D");
  hTreeKF->Branch("fImParzDaughter2", &fImParzDaughter2, "fImParzDaughter2/D");
  hTreeKF->Branch("fDcaSecDaughter", &fDcaSecDaughter, "fDcaSecDaughter/D");
  hTreeKF->Branch("fDcaSecDaughter2", &fDcaSecDaughter2, "fDcaSecDaughter2/D");

  hTreeKF->Branch("fTrackPIDDaughter", &fTrackPIDDaughter, "fTrackPIDDaughter/I");
  hTreeKF->Branch("fLabelDaughterKF", &fLabelDaughterKF, "fLabelDaughterKF/D");
  hTreeKF->Branch("fpDaughter", &fpDaughter, "fpDaughter/D");
  hTreeKF->Branch("fptDaughterUnProp", &fptDaughterUnProp, "fptDaughterUnProp/D");
  hTreeKF->Branch("fpxDaughterUnProp", &fpxDaughterUnProp, "fpxDaughterUnProp/D");
  hTreeKF->Branch("fpyDaughterUnProp", &fpyDaughterUnProp, "fpyDaughterUnProp/D");
  hTreeKF->Branch("fpzDaughterUnProp", &fpzDaughterUnProp, "fpzDaughterUnProp/D");
  hTreeKF->Branch("fdEdxDaughter", &fdEdxDaughter, "fdEdxDaughter/D");
  hTreeKF->Branch("fdEdxSigmaDaughter", &fdEdxSigmaDaughter, "fdEdxSigmaDaughter/D");
  hTreeKF->Branch("fDcaDaughter", &fDcaDaughter, "fDcaDaughter/D");
  hTreeKF->Branch("fDcaDaughtero", &fDcaDaughtero, "fDcaDaughtero/D");
  hTreeKF->Branch("fDcazDaughter", &fDcazDaughter, "fDcazDaughter/D");
  hTreeKF->Branch("fNclsDaughter", &fNclsDaughter, "fNclsDaughter/I");
  hTreeKF->Branch("fChi2Daughter", &fChi2Daughter, "fChi2Daughter/D");
  hTreeKF->Branch("fNclsITSDaughter", &fNclsITSDaughter, "fNclsITSDaughter/I");
  hTreeKF->Branch("fEtaDaughter", &fEtaDaughter, "fEtaDaughter/D");
  hTreeKF->Branch("fPhiDaughter", &fPhiDaughter, "fPhiDaughter/D");
  hTreeKF->Branch("fGeoLengthDaughter", &fGeoLengthDaughter, "fGeoLengthDaughter/D");
  hTreeKF->Branch("fTOFSignalDaughter", &fTOFSignalDaughter, "fTOFSignalDaughter/D");
  hTreeKF->Branch("fSigmaYXDaughter", &fSigmaYXDaughter, "fSigmaYXDaughter/D");
  hTreeKF->Branch("fSigmaXYZDaughter", &fSigmaXYZDaughter, "fSigmaXYZDaughter/D");
  hTreeKF->Branch("fSigmaZDaughter", &fSigmaZDaughter, "fSigmaZDaughter/D");
  hTreeKF->Branch("fPtUncertDaughter", &fPtUncertDaughter, "fPtUncertDaughter/D");
  hTreeKF->Branch("fTPCRefitDaughter", &fTPCRefitDaughter, "fTPCRefitDaughter/I");
  hTreeKF->Branch("fITSRefitDaughter", &fITSRefitDaughter, "fITSRefitDaughter/I");
  hTreeKF->Branch("fITSLayer1Daughter", &fITSLayer1Daughter, "fITSLayer1Daughter/I");
  hTreeKF->Branch("fITSLayer2Daughter", &fITSLayer2Daughter, "fITSLayer2Daughter/I");
  hTreeKF->Branch("fITSLayer3Daughter", &fITSLayer3Daughter, "fITSLayer3Daughter/I");
  hTreeKF->Branch("fITSLayer4Daughter", &fITSLayer4Daughter, "fITSLayer4Daughter/I");
  hTreeKF->Branch("fITSLayer5Daughter", &fITSLayer5Daughter, "fITSLayer5Daughter/I");
  hTreeKF->Branch("fITSLayer6Daughter", &fITSLayer6Daughter, "fITSLayer6Daughter/I");

  hTreeKF->Branch("fTrackPIDDaughter2", &fTrackPIDDaughter2, "fTrackPIDDaughter2/I");
  hTreeKF->Branch("fLabelDaughter2KF", &fLabelDaughter2KF, "fLabelDaughter2KF/D");
  hTreeKF->Branch("fpDaughter2", &fpDaughter2, "fpDaughter2/D");
  hTreeKF->Branch("fptDaughter2UnProp", &fptDaughter2UnProp, "fptDaughter2UnProp/D");
  hTreeKF->Branch("fpxDaughter2UnProp", &fpxDaughter2UnProp, "fpxDaughter2UnProp/D");
  hTreeKF->Branch("fpyDaughter2UnProp", &fpyDaughter2UnProp, "fpyDaughter2UnProp/D");
  hTreeKF->Branch("fpzDaughter2UnProp", &fpzDaughter2UnProp, "fpzDaughter2UnProp/D");
  hTreeKF->Branch("fdEdxDaughter2", &fdEdxDaughter2, "fdEdxDaughter2/D");
  hTreeKF->Branch("fdEdxSigmaDaughter2", &fdEdxSigmaDaughter2, "fdEdxSigmaDaughter2/D");
  hTreeKF->Branch("fDcaDaughter2", &fDcaDaughter2, "fDcaDaughter2/D");
  hTreeKF->Branch("fDcaDaughter2o", &fDcaDaughter2o, "fDcaDaughter2o/D");
  hTreeKF->Branch("fDcazDaughter2", &fDcazDaughter2, "fDcazDaughter2/D");
  hTreeKF->Branch("fNclsDaughter2", &fNclsDaughter2, "fNclsDaughter2/I");
  hTreeKF->Branch("fChi2Daughter2", &fChi2Daughter2, "fChi2Daughter2/D");
  hTreeKF->Branch("fNclsITSDaughter2", &fNclsITSDaughter2, "fNclsITSDaughter2/I");
  hTreeKF->Branch("fEtaDaughter2", &fEtaDaughter2, "fEtaDaughter2/D");
  hTreeKF->Branch("fPhiDaughter2", &fPhiDaughter2, "fPhiDaughter2/D");
  hTreeKF->Branch("fGeoLengthDaughter2", &fGeoLengthDaughter2, "fGeoLengthDaughter2/D");
  hTreeKF->Branch("fTOFSignalDaughter2", &fTOFSignalDaughter2, "fTOFSignalDaughter2/D");
  hTreeKF->Branch("fSigmaYXDaughter2", &fSigmaYXDaughter2, "fSigmaYXDaughter2/D");
  hTreeKF->Branch("fSigmaXYZDaughter2", &fSigmaXYZDaughter2, "fSigmaXYZDaughter2/D");
  hTreeKF->Branch("fSigmaZDaughter2", &fSigmaZDaughter2, "fSigmaZDaughter2/D");
  hTreeKF->Branch("fPtUncertDaughter2", &fPtUncertDaughter2, "fPtUncertDaughter2/D");
  hTreeKF->Branch("fTPCRefitDaughter2", &fTPCRefitDaughter2, "fTPCRefitDaughter2/I");
  hTreeKF->Branch("fITSRefitDaughter2", &fITSRefitDaughter2, "fITSRefitDaughter2/I");
  hTreeKF->Branch("fITSLayer1Daughter2", &fITSLayer1Daughter2, "fITSLayer1Daughter2/I");
  hTreeKF->Branch("fITSLayer2Daughter2", &fITSLayer2Daughter2, "fITSLayer2Daughter2/I");
  hTreeKF->Branch("fITSLayer3Daughter2", &fITSLayer3Daughter2, "fITSLayer3Daughter2/I");
  hTreeKF->Branch("fITSLayer4Daughter2", &fITSLayer4Daughter2, "fITSLayer4Daughter2/I");
  hTreeKF->Branch("fITSLayer5Daughter2", &fITSLayer5Daughter2, "fITSLayer5Daughter2/I");
  hTreeKF->Branch("fITSLayer6Daughter2", &fITSLayer6Daughter2, "fITSLayer6Daughter2/I");

  hTreeKF->Branch("fdEdxSigmaPion", &fdEdxSigmaPion, "fdEdxSigmaPion/D");
  hTreeKF->Branch("fdEdxSigmaDeuteron", &fdEdxSigmaDeuteron, "fdEdxSigmaDeuteron/D");
  hTreeKF->Branch("fdEdxSigmaTriton", &fdEdxSigmaTriton, "fdEdxSigmaTriton/D");
  hTreeKF->Branch("fdEdxSigmaAlpha", &fdEdxSigmaAlpha, "fdEdxSigmaAlpha/D");
  hTreeKF->Branch("fTrackPar", fTrackPar, "fTrackPar[7]/D");
  hTreeKF->Branch("fTrackPar2", fTrackPar2, "fTrackPar2[7]/D");
  hTreeKF->Branch("fCovMatrixTrack", fCovMatrixTrack, "fCovMatrixTrack[21]/D");
  hTreeKF->Branch("fCovMatrixTrack2", fCovMatrixTrack2, "fCovMatrixTrack2[21]/D");
  
  // _________________________________________________ //
  // __ special daughter checks tree __ //
  iTree = new TTree("iTree", "iTree");
  iTree->Branch("fEventID", &fEventID, "fEventID/I");
  iTree->Branch("fParticleID", &fParticleID, "fParticleID/I");
  iTree->Branch("fPeriod", &fPeriod, "fPeriod/I");
  iTree->Branch("frunnumber", &frunnumber, "frunnumber/I");
  iTree->Branch("feventclass", &feventclass, "feventclass/I");
  iTree->Branch("fCentrality", &fCentrality, "fCentrality/I");
  iTree->Branch("fMagneticField", &fMagneticField, "fMagneticField/I");
  iTree->Branch("fisWOTrackVertex", &fisWOTrackVertex, "fisWOTrackVertex/I");
  iTree->Branch("fTriggerString", (void*)fTriggerString, "fTriggerString/C");
  iTree->Branch("fTrigMB", &fTrigMB, "fTrigMB/I");
  iTree->Branch("fTrigHMV0", &fTrigHMV0, "fTrigHMV0/I");
  iTree->Branch("fTrigHMSPD", &fTrigHMSPD, "fTrigHMSPD/I");
  iTree->Branch("fTrigHNU", &fTrigHNU, "fTrigHNU/I");
  iTree->Branch("fTrigHQU", &fTrigHQU, "fTrigHQU/I");
  iTree->Branch("fTrigkCentral", &fTrigkCentral, "fTrigkCentral/I");
  iTree->Branch("fTrigkSemiCentral", &fTrigkSemiCentral, "fTrigkSemiCentral/I");

  iTree->Branch("fPrimVertexX", &fPrimVertexX, "fPrimVertexX/D");
  iTree->Branch("fPrimVertexY", &fPrimVertexY, "fPrimVertexY/D");
  iTree->Branch("fPrimVertexZ", &fPrimVertexZ, "fPrimVertexZ/D");
  iTree->Branch("fTertVertexX", &fTertVertexX, "fTertVertexX/D");
  iTree->Branch("fTertVertexY", &fTertVertexY, "fTertVertexY/D");
  iTree->Branch("fTertVertexZ", &fTertVertexZ, "fTertVertexZ/D");
  iTree->Branch("fDCA3B1", &fDCA3B1, "fDCA3B1/D");
  iTree->Branch("fSubPA", &fSubPA, "fSubPA/D");

  iTree->Branch("fPDGMother", &fPDGMother, "fPDGMother/I");
  iTree->Branch("fChargeMother", &fChargeMother, "fChargeMother/I");
  iTree->Branch("fDecayChannel", &fDecayChannel, "fDecayChannel/I");
  iTree->Branch("fRecoMethod", &fRecoMethod, "fRecoMethod/I");

  iTree->Branch("fmMotherCheck", &fmMotherCheck, "fmMotherCheck/D");
  iTree->Branch("fmSubMotherCheck", &fmSubMotherCheck, "fmSubMotherCheck/D");
  iTree->Branch("fmSubMother", &fmSubMother, "fmSubMother/D");
  iTree->Branch("fESubMother", &fESubMother, "fESubMother/D");
  iTree->Branch("fpxSubMother", &fpxSubMother, "fpxSubMother/D");
  iTree->Branch("fpySubMother", &fpySubMother, "fpySubMother/D");
  iTree->Branch("fpzSubMother", &fpzSubMother, "fpzSubMother/D");
  iTree->Branch("fptSubMother", &fptSubMother, "fptSubMother/D");
  iTree->Branch("fpSubMother", &fpSubMother, "fpSubMother/D");
  iTree->Branch("fySubMother", &fySubMother, "fySubMother/D");
  iTree->Branch("fctSubMother", &fctSubMother, "fctSubMother/D");

  iTree->Branch("fptDaughter", &fptDaughter, "fptDaughter/D");
  iTree->Branch("fpxDaughter", &fpxDaughter, "fpxDaughter/D");
  iTree->Branch("fpyDaughter", &fpyDaughter, "fpyDaughter/D");
  iTree->Branch("fpzDaughter", &fpzDaughter, "fpzDaughter/D");
  iTree->Branch("fyDaughter", &fyDaughter, "fyDaughter/D");
  iTree->Branch("fptDaughter2", &fptDaughter2, "fptDaughter2/D");
  iTree->Branch("fpxDaughter2", &fpxDaughter2, "fpxDaughter2/D");
  iTree->Branch("fpyDaughter2", &fpyDaughter2, "fpyDaughter2/D");
  iTree->Branch("fpzDaughter2", &fpzDaughter2, "fpzDaughter2/D");
  iTree->Branch("fyDaughter2", &fyDaughter2, "fyDaughter2/D");
  iTree->Branch("fPropDCADaughter", &fPropDCADaughter, "fPropDCADaughter/I");
  iTree->Branch("fImParDaughter", &fImParDaughter, "fImParDaughter/D");
  iTree->Branch("fImParzDaughter", &fImParzDaughter, "fImParzDaughter/D");
  iTree->Branch("fPropDCADaughter2", &fPropDCADaughter2, "fPropDCADaughter2/I");
  iTree->Branch("fImParDaughter2", &fImParDaughter2, "fImParDaughter2/D");
  iTree->Branch("fImParzDaughter2", &fImParzDaughter2, "fImParzDaughter2/D");
  iTree->Branch("fDcaSecDaughter", &fDcaSecDaughter, "fDcaSecDaughter/D");
  iTree->Branch("fDcaSecDaughter2", &fDcaSecDaughter2, "fDcaSecDaughter2/D");

  // _________________________________________________ //
  // __ special daughter checks KF tree __ //
  iTreeKF = new TTree("iTreeKF", "iTreeKF");
  iTreeKF->Branch("fEventID", &fEventID, "fEventID/I");
  iTreeKF->Branch("fParticleID", &fParticleID, "fParticleID/I");
  iTreeKF->Branch("fPeriod", &fPeriod, "fPeriod/I");
  iTreeKF->Branch("frunnumber", &frunnumber, "frunnumber/I");
  iTreeKF->Branch("feventclass", &feventclass, "feventclass/I");
  iTreeKF->Branch("fCentrality", &fCentrality, "fCentrality/I");
  iTreeKF->Branch("fMagneticField", &fMagneticField, "fMagneticField/I");
  iTreeKF->Branch("fisWOTrackVertex", &fisWOTrackVertex, "fisWOTrackVertex/I");
  iTreeKF->Branch("fTriggerString", (void*)fTriggerString, "fTriggerString/C");
  iTreeKF->Branch("fTrigMB", &fTrigMB, "fTrigMB/I");
  iTreeKF->Branch("fTrigHMV0", &fTrigHMV0, "fTrigHMV0/I");
  iTreeKF->Branch("fTrigHMSPD", &fTrigHMSPD, "fTrigHMSPD/I");
  iTreeKF->Branch("fTrigHNU", &fTrigHNU, "fTrigHNU/I");
  iTreeKF->Branch("fTrigHQU", &fTrigHQU, "fTrigHQU/I");
  iTreeKF->Branch("fTrigkCentral", &fTrigkCentral, "fTrigkCentral/I");
  iTreeKF->Branch("fTrigkSemiCentral", &fTrigkSemiCentral, "fTrigkSemiCentral/I");

  iTreeKF->Branch("fPrimVertexXKF", &fPrimVertexXKF, "fPrimVertexXKF/D");
  iTreeKF->Branch("fPrimVertexYKF", &fPrimVertexYKF, "fPrimVertexYKF/D");
  iTreeKF->Branch("fPrimVertexZKF", &fPrimVertexZKF, "fPrimVertexZKF/D");
  iTreeKF->Branch("fPrimVertexXErrKF", &fPrimVertexXErrKF, "fPrimVertexXErrKF/D");
  iTreeKF->Branch("fPrimVertexYErrKF", &fPrimVertexYErrKF, "fPrimVertexYErrKF/D");
  iTreeKF->Branch("fPrimVertexZErrKF", &fPrimVertexZErrKF, "fPrimVertexZErrKF/D");
  iTreeKF->Branch("fPrimVertChi2KF", &fPrimVertChi2KF, "fPrimVertChi2KF/D");
  iTreeKF->Branch("fPrimVertNDFKF", &fPrimVertNDFKF, "fPrimVertNDFKF/I");
  iTreeKF->Branch("fTertVertexXKF", &fTertVertexXKF, "fTertVertexXKF/D");
  iTreeKF->Branch("fTertVertexYKF", &fTertVertexYKF, "fTertVertexYKF/D");
  iTreeKF->Branch("fTertVertexZKF", &fTertVertexZKF, "fTertVertexZKF/D");
  iTreeKF->Branch("fTertVertexXErrKF", &fTertVertexXErrKF, "fTertVertexXErrKF/D");
  iTreeKF->Branch("fTertVertexYErrKF", &fTertVertexYErrKF, "fTertVertexYErrKF/D");
  iTreeKF->Branch("fTertVertexZErrKF", &fTertVertexZErrKF, "fTertVertexZErrKF/D");
  iTreeKF->Branch("fTertVertChi2KF", &fTertVertChi2KF, "fTertVertChi2KF/D");
  iTreeKF->Branch("fTertVertNDFKF", &fTertVertNDFKF, "fTertVertNDFKF/I");
  iTreeKF->Branch("fDCA3B1XYKF", &fDCA3B1XYKF, "fDCA3B1XYKF/D");
  iTreeKF->Branch("fDCA3B1ZKF", &fDCA3B1ZKF, "fDCA3B1ZKF/D");
  iTreeKF->Branch("fSubPAKF", &fSubPAKF, "fSubPAKF/D");

  iTreeKF->Branch("fPDGMother", &fPDGMother, "fPDGMother/I");
  iTreeKF->Branch("fChargeMother", &fChargeMother, "fChargeMother/I");
  iTreeKF->Branch("fDecayChannel", &fDecayChannel, "fDecayChannel/I");
  iTreeKF->Branch("fRecoMethod", &fRecoMethod, "fRecoMethod/I");

  iTreeKF->Branch("fmMotherCheck", &fmMotherCheck, "fmMotherCheck/D");
  iTreeKF->Branch("fmSubMotherCheck", &fmSubMotherCheck, "fmSubMotherCheck/D");
  iTreeKF->Branch("fmSubMother", &fmSubMother, "fmSubMother/D");
  iTreeKF->Branch("fmSubMotherKF", &fmSubMotherKF, "fmSubMotherKF/D");
  iTreeKF->Branch("fmSubMotherErrKF", &fmSubMotherErrKF, "fmSubMotherErrKF/D");
  iTreeKF->Branch("fESubMotherKF", &fESubMotherKF, "fESubMotherKF/D");
  iTreeKF->Branch("fESubMotherErrKF", &fESubMotherErrKF, "fESubMotherErrKF/D");
  iTreeKF->Branch("fpxSubMotherKF", &fpxSubMotherKF, "fpxSubMotherKF/D");
  iTreeKF->Branch("fpxSubMotherErrKF", &fpxSubMotherErrKF, "fpxSubMotherErrKF/D");
  iTreeKF->Branch("fpySubMotherKF", &fpySubMotherKF, "fpySubMotherKF/D");
  iTreeKF->Branch("fpySubMotherErrKF", &fpySubMotherErrKF, "fpySubMotherErrKF/D");
  iTreeKF->Branch("fpzSubMotherKF", &fpzSubMotherKF, "fpzSubMotherKF/D");
  iTreeKF->Branch("fpzSubMotherErrKF", &fpzSubMotherErrKF, "fpzSubMotherErrKF/D");
  iTreeKF->Branch("fptSubMotherKF", &fptSubMotherKF, "fptSubMotherKF/D");
  iTreeKF->Branch("fptSubMotherErrKF", &fptSubMotherErrKF, "fptSubMotherErrKF/D");
  iTreeKF->Branch("fpSubMotherKF", &fpSubMotherKF, "fpSubMotherKF/D");
  iTreeKF->Branch("fpSubMotherErrKF", &fpSubMotherErrKF, "fpSubMotherErrKF/D");
  iTreeKF->Branch("fySubMotherKF", &fySubMotherKF, "fySubMotherKF/D");
  iTreeKF->Branch("fctSubMotherKF", &fctSubMotherKF, "fctSubMotherKF/D");
  iTreeKF->Branch("fctSubMotherErrKF", &fctSubMotherErrKF, "fctSubMotherErrKF/D");

  iTreeKF->Branch("fptDaughter", &fptDaughter, "fptDaughter/D");
  iTreeKF->Branch("fpxDaughter", &fpxDaughter, "fpxDaughter/D");
  iTreeKF->Branch("fpyDaughter", &fpyDaughter, "fpyDaughter/D");
  iTreeKF->Branch("fpzDaughter", &fpzDaughter, "fpzDaughter/D");
  iTreeKF->Branch("fyDaughter", &fyDaughter, "fyDaughter/D");
  iTreeKF->Branch("fptDaughter2", &fptDaughter2, "fptDaughter2/D");
  iTreeKF->Branch("fpxDaughter2", &fpxDaughter2, "fpxDaughter2/D");
  iTreeKF->Branch("fpyDaughter2", &fpyDaughter2, "fpyDaughter2/D");
  iTreeKF->Branch("fpzDaughter2", &fpzDaughter2, "fpzDaughter2/D");
  iTreeKF->Branch("fyDaughter2", &fyDaughter2, "fyDaughter2/D");
  iTreeKF->Branch("fDcaDaughterXYKF", &fDcaDaughterXYKF, "fDcaDaughterXYKF/D");
  iTreeKF->Branch("fDcaDaughter2XYKF", &fDcaDaughter2XYKF, "fDcaDaughter2XYKF/D");
  iTreeKF->Branch("fDcaDaughterZKF", &fDcaDaughterZKF, "fDcaDaughterZKF/D");
  iTreeKF->Branch("fDcaDaughter2ZKF", &fDcaDaughter2ZKF, "fDcaDaughter2ZKF/D");
  iTreeKF->Branch("fDcaSecDaughterXYKF", &fDcaSecDaughterXYKF, "fDcaSecDaughterXYKF/D");
  iTreeKF->Branch("fDcaSecDaughterZKF", &fDcaSecDaughterZKF, "fDcaSecDaughterZKF/D");
  iTreeKF->Branch("fDcaSecDaughter2XYKF", &fDcaSecDaughter2XYKF, "fDcaSecDaughter2XYKF/D");
  iTreeKF->Branch("fDcaSecDaughter2ZKF", &fDcaSecDaughter2ZKF, "fDcaSecDaughter2ZKF/D");
  iTreeKF->Branch("fPropDCADaughter", &fPropDCADaughter, "fPropDCADaughter/I");
  iTreeKF->Branch("fImParDaughter", &fImParDaughter, "fImParDaughter/D");
  iTreeKF->Branch("fImParzDaughter", &fImParzDaughter, "fImParzDaughter/D");
  iTreeKF->Branch("fPropDCADaughter2", &fPropDCADaughter2, "fPropDCADaughter2/I");
  iTreeKF->Branch("fImParDaughter2", &fImParDaughter2, "fImParDaughter2/D");
  iTreeKF->Branch("fImParzDaughter2", &fImParzDaughter2, "fImParzDaughter2/D");
  iTreeKF->Branch("fDcaSecDaughter", &fDcaSecDaughter, "fDcaSecDaughter/D");
  iTreeKF->Branch("fDcaSecDaughter2", &fDcaSecDaughter2, "fDcaSecDaughter2/D");

  // _________________________________________________ //
  // __ generated tree 4LLH __ //
  fTreeGen = new TTree("fTreeGen", "fTreeGen");
  fTreeGen->Branch("fPeriod", &fPeriod, "fPeriod/I");
  fTreeGen->Branch("frunnumber", &frunnumber, "frunnumber/I");
  fTreeGen->Branch("fCentrality", &fCentrality, "fCentrality/I");
  fTreeGen->Branch("fMagneticField", &fMagneticField, "fMagneticField/I");
  fTreeGen->Branch("fTrigMB", &fTrigMB, "fTrigMB/I");
  fTreeGen->Branch("fTrigHMV0", &fTrigHMV0, "fTrigHMV0/I");
  fTreeGen->Branch("fTrigHMSPD", &fTrigHMSPD, "fTrigHMSPD/I");
  fTreeGen->Branch("fTrigHNU", &fTrigHNU, "fTrigHNU/I");
  fTreeGen->Branch("fTrigHQU", &fTrigHQU, "fTrigHQU/I");
  fTreeGen->Branch("fTrigkCentral", &fTrigkCentral, "fTrigkCentral/I");
  fTreeGen->Branch("fTrigkSemiCentral", &fTrigkSemiCentral, "fTrigkSemiCentral/I");
  fTreeGen->Branch("fisWOTrackVertex", &fisWOTrackVertex, "fisWOTrackVertex/I");
  fTreeGen->Branch("fTriggerString", (void*)fTriggerString, "fTriggerString[]/C");
  fTreeGen->Branch("fPDGMother", &fPDGMother, "fPDGMother/I");
  fTreeGen->Branch("fDecayChannel", &fDecayChannel, "fDecayChannel/I");
  fTreeGen->Branch("fChargeMother", &fChargeMother, "fChargeMother/I");
  fTreeGen->Branch("fmMother", &fmMother, "fmMother/D");
  fTreeGen->Branch("fptMother", &fptMother, "fptMother/D");
  fTreeGen->Branch("fyMother", &fyMother, "fyMother/D");
  fTreeGen->Branch("fctMother", &fctMother, "fctMother/D");
  fTreeGen->Branch("fmSubMother", &fmSubMother, "fmSubMother/D");
  fTreeGen->Branch("fptSubMother", &fptSubMother, "fptSubMother/D");
  fTreeGen->Branch("fySubMother", &fySubMother, "fySubMother/D");
  fTreeGen->Branch("fctSubMother", &fctSubMother, "fctSubMother/D");	
  fTreeGen->Branch("fptDaughter", &fptDaughter, "fptDaughter/D");
  fTreeGen->Branch("fyDaughter", &fyDaughter, "fyDaughter/D");
  fTreeGen->Branch("fptDaughter1", &fptDaughter1, "fptDaughter1/D");
  fTreeGen->Branch("fyDaughter1", &fyDaughter1, "fyDaughter1/D");
  fTreeGen->Branch("fptDaughter2", &fptDaughter2, "fptDaughter2/D");
  fTreeGen->Branch("fyDaughter2", &fyDaughter2, "fyDaughter2/D");
  fTreeGen->Branch("fptDaughter3", &fptDaughter3, "fptDaughter3/D");
  fTreeGen->Branch("fyDaughter3", &fyDaughter3, "fyDaughter3/D");
  fTreeGen->Branch("fPrimary4LHe", &fPrimary4LHe, "fPrimary4LHe/I");

  // _________________________________________________ //
  // __ generated tree 4LHe + 4LH __ //
  gTreeGen = new TTree("gTreeGen", "gTreeGen");
  gTreeGen->Branch("fPeriod", &fPeriod, "fPeriod/I");
  gTreeGen->Branch("frunnumber", &frunnumber, "frunnumber/I");
  gTreeGen->Branch("fCentrality", &fCentrality, "fCentrality/I");
  gTreeGen->Branch("fMagneticField", &fMagneticField, "fMagneticField/I");
  gTreeGen->Branch("fTrigMB", &fTrigMB, "fTrigMB/I");
  gTreeGen->Branch("fTrigHMV0", &fTrigHMV0, "fTrigHMV0/I");
  gTreeGen->Branch("fTrigHMSPD", &fTrigHMSPD, "fTrigHMSPD/I");
  gTreeGen->Branch("fTrigHNU", &fTrigHNU, "fTrigHNU/I");
  gTreeGen->Branch("fTrigHQU", &fTrigHQU, "fTrigHQU/I");
  gTreeGen->Branch("fTrigkCentral", &fTrigkCentral, "fTrigkCentral/I");
  gTreeGen->Branch("fTrigkSemiCentral", &fTrigkSemiCentral, "fTrigkSemiCentral/I");
  gTreeGen->Branch("fisWOTrackVertex", &fisWOTrackVertex, "fisWOTrackVertex/I");
  gTreeGen->Branch("fTriggerString", (void*)fTriggerString, "fTriggerString[]/C");
  gTreeGen->Branch("fPDGMother", &fPDGMother, "fPDGMother/I");
  gTreeGen->Branch("fisExcited", &fisExcited, "fisExcited/I");
  gTreeGen->Branch("fChargeMother", &fChargeMother, "fChargeMother/I");
  gTreeGen->Branch("fDecayChannel", &fDecayChannel, "fDecayChannel/I");
  gTreeGen->Branch("fmSubMother", &fmSubMother, "fmSubMother/D");
  gTreeGen->Branch("fptSubMother", &fptSubMother, "fptSubMother/D");
  gTreeGen->Branch("fySubMother", &fySubMother, "fySubMother/D");
  gTreeGen->Branch("fctSubMother", &fctSubMother, "fctSubMother/D");
  gTreeGen->Branch("fptDaughter", &fptDaughter, "fptDaughter/D");
  gTreeGen->Branch("fyDaughter", &fyDaughter, "fyDaughter/D");
  gTreeGen->Branch("fptDaughter1", &fptDaughter1, "fptDaughter1/D");
  gTreeGen->Branch("fyDaughter1", &fyDaughter1, "fyDaughter1/D");
  gTreeGen->Branch("fptDaughter2", &fptDaughter2, "fptDaughter2/D");
  gTreeGen->Branch("fyDaughter2", &fyDaughter2, "fyDaughter2/D");
  gTreeGen->Branch("fPrimary4LHe", &fPrimary4LHe, "fPrimary4LHe/I");

  PostData(1,  fHistogramList);
  PostData(2,  fTree);
  PostData(3,  fTreeKF);
  PostData(4,  gTree);
  PostData(5,  gTreeKF); 
  PostData(6,  hTree);
  PostData(7,  hTreeKF);
  PostData(8,  iTree);
  PostData(9,  iTreeKF);
  PostData(10, fTreeGen);
  PostData(11, gTreeGen);


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
void AliAnalysisTaskDoubleHypNucTreeLS::UserExec(Option_t*) {
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
  fAODevent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAODevent) {
    return;
  }
  fevent      = dynamic_cast<AliVEvent*>(InputEvent());
  eventHeader = fevent->GetHeader();
  fEventID    = eventHeader->GetEventIdAsLong();

  // __ runnumber & period __ //
  fPeriod     = fAODevent->GetPeriodNumber();
  frunnumber  = fAODevent->GetRunNumber();

  fHistNumEvents->Fill(0);
  // __ EventCuts & TimeRangeCut __ //
  if (!fMCtrue && (frunnumber == 297219 || frunnumber == 297194 || frunnumber == 297029
		   || frunnumber == 296890 || frunnumber == 296849 || frunnumber == 296750
		   || frunnumber == 296749 || frunnumber == 297481)) fEventCuts.UseTimeRangeCut();

  // __ classify events __ //
  fEventCuts.fRequireTrackVertex = kTRUE; /// switch off track vertex requirement ///
  // __ Trigger Selection __ //
  fEventCuts.OverrideAutomaticTriggerSelection(fTriggerMask);

  // __ Number of acc events __ //
  if (fEventCuts.AcceptEvent(fevent)) feventclass = 1;
  else return;//feventclass = 0;
  // All accepted events
  if (feventclass)fHistNumEvents->Fill(1);
  //  accepted events depending on Track Vertex reco
  if (feventclass)fHistNumEvents->Fill(fEventCuts.PassedCut(AliEventCuts::kVertexTracks) == 1 ? 2 : 3);

  fisWOTrackVertex = fEventCuts.PassedCut(AliEventCuts::kVertexTracks) == 1 ? 0 : 1;	

  // __ Trigger Selection __ //
  TString tTriggerString = fAODevent->GetFiredTriggerClasses();
  fTriggerString = tTriggerString;
  AliAnalysisTaskDoubleHypNucTreeLS::TriggerSelection();

  // __ centrality, 0 = V0M __ //
  Float_t centrality = -1;
  AliMultSelection* fMultSelection = (AliMultSelection*)fAODevent->FindListObject("MultSelection");
  centrality = fMultSelection->GetMultiplicityPercentile("V0M");
  fCentrality = centrality;

  if (!fisWOTrackVertex) fHistCentrality1->Fill(fCentrality);
  if (fisWOTrackVertex)  fHistCentrality2->Fill(fCentrality);

  if (!fisWOTrackVertex) fHistNtracks1->Fill(fAODevent->GetNumberOfTracks());
  if (fisWOTrackVertex)  fHistNtracks2->Fill(fAODevent->GetNumberOfTracks());

	
  // __ trigger __ //
  fTrigMB = MB;
  fTrigHMV0 = HMV0;
  fTrigHMSPD = HMSPD;
  fTrigHNU = HNU;
  fTrigHQU = HQU;
  fTrigkCentral = Central;
  fTrigkSemiCentral = SemiCentral;

  //__ SetBetheBlochParams __ //
  AliAnalysisTaskDoubleHypNucTreeLS::SetBetheBlochParams(frunnumber);

  // __ MagneticField __ //
  kMagF = fAODevent->GetMagneticField(); //double
  fMagneticField = kMagF;               // saved as int
  KFParticle::SetField(kMagF);

  // __ setting up vertex reconstruction __ //
  kStandardReco = 0;
  kKFReco = 1;
  // __ 1 = reconstruction with coord of prim vtx __ //
  // __ 2 = reconstruction with coord of tert vtx __ //
  kVariante = 1;
  // __ skip the 3LH reco part __ //
  kSkip3LHAnalysis = 1;
  kSkip4LHAnalysis = 0;
  // __ copy prim vtx __ //
  vertex = fAODevent->GetPrimaryVertexSPD();

  // __ MC Settings __ //
  if (fMCtrue) {
    konlyBG = kFALSE;
    konlySig = kFALSE;
    if (konlySig) konlyBG = kFALSE;
    kMCPIDCheck = kTRUE;
    fEventID = TMath::Abs(gRandom->Uniform(100000000, 987654321) - gRandom->Uniform(0.2, 0.8) * frunnumber);
  }

  // __ initialize some utils __ //
  sublorentzsum  = new TLorentzVector(0., 0., 0., 0.);
  sublorentzsum2 = new TLorentzVector(0., 0., 0., 0.);
  lorentzsum     = new TLorentzVector(0., 0., 0., 0.);
  lorentzsum2    = new TLorentzVector(0., 0., 0., 0.);
  particle1      = new TLorentzVector(0., 0., 0., 0.);
  particle2      = new TLorentzVector(0., 0., 0., 0.);
  particle3      = new TLorentzVector(0., 0., 0., 0.);
  particle4      = new TLorentzVector(0., 0., 0., 0.);
  exTrack        = new AliExternalTrackParam();
  exTrack1       = new AliExternalTrackParam();
  exTrack2       = new AliExternalTrackParam();
  exTrack3       = new AliExternalTrackParam();
  exTrack4       = new AliExternalTrackParam();
  h              = new TVector3(0., 0., 0.);
  // __ only PID __ //
  if (fPIDCheckOnly) {
    AliAnalysisTaskDoubleHypNucTreeLS::dEdxCheck();
  }// __ Analysis __ //
  else {
    if (fMCtrue) AliAnalysisTaskDoubleHypNucTreeLS::MCGenerated();
    AliAnalysisTaskDoubleHypNucTreeLS::FindTracks();
  }
  // __ delete the utils __ //
  if (sublorentzsum)  delete sublorentzsum;
  if (sublorentzsum2) delete sublorentzsum2;
  if (lorentzsum)     delete lorentzsum;
  if (lorentzsum2)    delete lorentzsum2;
  if (particle1)      delete particle1;
  if (particle2)      delete particle2;
  if (particle3)      delete particle3;
  if (particle4)      delete particle4;
  if (h)              delete h;
  if (exTrack)        delete exTrack;
  if (exTrack1)       delete exTrack1;
  if (exTrack2)       delete exTrack2;
  if (exTrack3)       delete exTrack3;
  if (exTrack4)       delete exTrack4;

  PostData(1,  fHistogramList);
  PostData(2,  fTree);
  PostData(3,  fTreeKF);
  PostData(4,  gTree);
  PostData(5,  gTreeKF); 
  PostData(6,  hTree);
  PostData(7,  hTreeKF);
  PostData(8,  iTree);
  PostData(9,  iTreeKF);
  PostData(10, fTreeGen);
  PostData(11, gTreeGen);
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTreeLS::dEdxCheck() {
  for (Int_t itrack = 0; itrack < fAODevent->GetNumberOfTracks(); itrack++) {

    AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(itrack));
    Double_t momentum = track->GetTPCmomentum();
    fHistdEdx->Fill(momentum * track->GetSign(), track->GetTPCsignal());
  }
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTreeLS::InitArrays() {

  const Int_t nTracksEvent = fAODevent->GetNumberOfTracks();
  // __ Track Array for identified particles __ //
  He4PosCounter   = 0;
  He4NegCounter   = 0;
  He3PosCounter   = 0;
  He3NegCounter   = 0;
  PPosCounter     = 0;
  PNegCounter     = 0;
  PiPosCounter    = 0;
  PiNegCounter    = 0;
  PiPosSecCounter = 0;
  PiNegSecCounter = 0;

  He3PosArray.clear();
  He4PosArray.clear();
  PPosArray.clear();
  PiPosArray.clear();
  PiPosSecArray.clear();
  He3NegArray.clear();
  He4NegArray.clear();
  PNegArray.clear();
  PiNegArray.clear();
  PiNegSecArray.clear();

  // __ cuts __ //
  kDCATracksCut     = 2.;
  kTrackPtUncertCut = 1.0;
  kPointingAngleCut = 0.9;
  kPIDSelecCut      = 4.0;
  int nsigma        = 5;
  kMin4LLHMass      = 4.04;//4.106 - nsigma * 0.004;
  kMax4LLHMass      = 4.18;//4.106 + nsigma * 0.004;
  kMin4LHeMass      = 3.86;//3.924 - nsigma * 0.003;
  kMax4LHeMass      = 4.00;//3.924 + nsigma * 0.003;
  kMin4LHMass       = 3.86;//3.925 - nsigma * 0.003;
  kMax4LHMass       = 4.00;//3.925 + nsigma * 0.003;
  kMin3LHMass       = 2.92;//2.991 - nsigma * 0.002;
  kMax3LHMass       = 3.06;//2.991 + nsigma * 0.002; 
}
// _________________________________________________ //
Int_t AliAnalysisTaskDoubleHypNucTreeLS::CustomTrackCut(const AliAODTrack& track, Int_t particle) {//1 = 3He, 2 = p, 3 = pi
  Float_t parxy, parz;
  track.GetImpactParameters(parxy, parz);
  Double_t mass = 0.0, charge = 1.0;
  if (particle == 1) { mass = AliPID::ParticleMass(AliPID::kHe3);    charge = 2.; }
  if (particle == 2) { mass = AliPID::ParticleMass(AliPID::kProton); charge = 1.; }
  if (particle == 3) { mass = AliPID::ParticleMass(AliPID::kPion);   charge = 1.; }
  TLorentzVector part(0., 0., 0., 0.);
  part.SetXYZM(charge * track.Px(), charge * track.Py(), charge * track.Pz(), mass);
  Int_t rval = 1;
  // __________ //
  // __ reject kinks __ //
  AliAODVertex *vtx1 = (AliAODVertex*)track.GetProdVertex();
  if (Int_t(vtx1->GetType()) == AliAODVertex::kKink)                       rval = 0;
  // __ restrict eta __ //
  if (TMath::Abs(track.Eta()) > 0.9)                                       rval = 0;
  // __ Chi2perTPCCluster __ //
  if ((track.GetTPCchi2() / (float)track.GetTPCNcls()) > 7.0)              rval = 0;
  // __ require TPC refit __ //
  if (!(track.GetStatus() & AliAODTrack::kTPCrefit))                       rval = 0;
  // __ pion momentum __ //
  if (particle == 3 && part.Pt() > 2.)                                     rval = 0;
  // __ proton momentum __ //
  //if (particle == 2 && part.Pt() > 5.)                                     rval = 0;
  // __ he momentum __ //
  //if (particle == 1 && (part.Pt() > 7. || part.Pt() < 0.8))                rval = 0;
  // __ ncls TPC pion __ //
  if (particle == 3 && track.GetTPCNcls() < 50)                            rval = 0;
  // __ ncls TPC proton __ //
  if (particle == 2 && track.GetTPCNcls() < 50)                            rval = 0;
  // __ ncls TPC He3 __ //
  if (particle == 1 && track.GetTPCNcls() < 50)                            rval = 0;
  // __ pion dca to prim vertex __ //
  //if (particle == 3 && TMath::Abs(parxy) < 0.005)                        rval = 0;
  // __ rel 1 pt uncert __ //
  //AliExternalTrackParam* p = new AliExternalTrackParam();
  //p->CopyFromVTrack(&track);
  //if ((TMath::Sqrt(p->GetSigma1Pt2()) * part.Pt()) > kTrackPtUncertCut)  rval = 0;
  //
  return rval;
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTreeLS::FindTracks() {

  const Int_t nTracksEvent = fAODevent->GetNumberOfTracks();

  // __ init track arrays for identified particles __ //
  AliAnalysisTaskDoubleHypNucTreeLS::InitArrays();

  // __ track loop for PID __ //
  for (Int_t qTracks = 0; qTracks < nTracksEvent; qTracks++) {

    AliAODTrack *trackq=dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(qTracks));

    if (!trackq) continue;//

    // __ fill energy loss plot __ //   
    fHistdEdx->Fill(trackq->GetTPCmomentum() * trackq->GetSign(), trackq->GetTPCsignal());

    // __ DATA __ //
    if (!fMCtrue) {
      // __ positive tracks __ //
      if (trackq->GetSign() > 0) {
	if (TMath::Abs(AliAnalysisTaskDoubleHypNucTreeLS::Bethe(*trackq, AliPID::ParticleMass(AliPID::kAlpha), 2, fBetheParamsHe)) <= kPIDSelecCut) {
	  if (!CustomTrackCut(*trackq, 1)) continue;
	  He4PosArray.push_back(qTracks);
	  He4PosCounter++;
	}
	if (TMath::Abs(AliAnalysisTaskDoubleHypNucTreeLS::Bethe(*trackq, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) <= kPIDSelecCut) {
	  if (!CustomTrackCut(*trackq, 1)) continue;
	  He3PosArray.push_back(qTracks);
	  He3PosCounter++;
	}
	if (TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kProton)) <= kPIDSelecCut) {
	  if (!CustomTrackCut(*trackq, 2)) continue;
	  PPosArray.push_back(qTracks);
	  PPosCounter++;
	}
	// -- like-sign flip is here -- //
	if (TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kPion)) <= kPIDSelecCut) {
	  if (!CustomTrackCut(*trackq, 3)) continue;
	  PiNegArray.push_back(qTracks);
	  PiNegCounter++;
	  PiNegSecArray.push_back(qTracks);
	  PiNegSecCounter++;
	}
      }
      // __ negative tracks __ //
      if (trackq->GetSign() < 0) {
	if (TMath::Abs(AliAnalysisTaskDoubleHypNucTreeLS::Bethe(*trackq, AliPID::ParticleMass(AliPID::kAlpha), 2, fBetheParamsHe)) <= kPIDSelecCut) {
	  if (!CustomTrackCut(*trackq, 1)) continue;
	  He4NegArray.push_back(qTracks);
	  He4NegCounter++;
	}
	if (TMath::Abs(AliAnalysisTaskDoubleHypNucTreeLS::Bethe(*trackq, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) <= kPIDSelecCut) {
	  if (!CustomTrackCut(*trackq, 1)) continue;
	  He3NegArray.push_back(qTracks);
	  He3NegCounter++;
	}
	if (TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kProton)) <= kPIDSelecCut) {
	  if (!CustomTrackCut(*trackq, 2)) continue;
	  PNegArray.push_back(qTracks);
	  PNegCounter++;
	}
	// -- like-sign flip is here -- //
	if (TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kPion)) <= kPIDSelecCut) {
	  if (!CustomTrackCut(*trackq, 3)) continue;
	  PiPosArray.push_back(qTracks);
	  PiPosCounter++;
	  PiPosSecArray.push_back(qTracks);
	  PiPosSecCounter++;
	}
      }
    }
    // ___ MC PART __ //
    if (fMCtrue) {
      label1 = trackq->GetLabel();
      AliMCParticle* tparticle = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(label1))->Particle());
      // __ positive tracks __ //
      if (trackq->GetSign() > 0) {
	if ((kMCPIDCheck && TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kAlpha)) <= kPIDSelecCut)
	    || (!kMCPIDCheck && tparticle->PdgCode() == fgkPdgCode[kPDGHelium4])) {
	  if (!CustomTrackCut(*trackq, 1)) continue;
	  He4PosArray.push_back(qTracks);
	  He4PosCounter++;
	}
	if ((kMCPIDCheck && TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kHe3)) <= kPIDSelecCut)
	    || (!kMCPIDCheck && tparticle->PdgCode() == fgkPdgCode[kPDGHelium3])) {
	  if (!CustomTrackCut(*trackq, 1)) continue;
	  He3PosArray.push_back(qTracks);
	  He3PosCounter++;
	}	
	if ((kMCPIDCheck && TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kProton)) <= kPIDSelecCut)
	    || (!kMCPIDCheck && tparticle->PdgCode() == fgkPdgCode[kPDGProton])) {
	  if (!CustomTrackCut(*trackq, 2)) continue;
	  PPosArray.push_back(qTracks);
	  PPosCounter++;
	}
	if ((kMCPIDCheck && TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kPion)) <= kPIDSelecCut)
	    || (!kMCPIDCheck && tparticle->PdgCode() == fgkPdgCode[kPDGPionPlus])) {
	  if (!CustomTrackCut(*trackq, 3)) continue;
	  PiPosArray.push_back(qTracks);
	  PiPosCounter++;
	  PiPosSecArray.push_back(qTracks);
	  PiPosSecCounter++;
	}
      }
      // __ negative tracks __ //
      if (trackq->GetSign() < 0) {
	if ((kMCPIDCheck && TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kAlpha)) <= kPIDSelecCut)
	    || (!kMCPIDCheck && tparticle->PdgCode() == fgkPdgCode[kPDGAntiHelium4])) {
	  if (!CustomTrackCut(*trackq, 1)) continue;
	  He4NegArray.push_back(qTracks);
	  He4NegCounter++;
	}
	if ((kMCPIDCheck && TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kHe3)) <= kPIDSelecCut)
	    || (!kMCPIDCheck && tparticle->PdgCode() == fgkPdgCode[kPDGAntiHelium3])) {
	  if (!CustomTrackCut(*trackq, 1)) continue;
	  He3NegArray.push_back(qTracks);
	  He3NegCounter++;
	}
	if ((kMCPIDCheck && TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kProton)) <= kPIDSelecCut)
	    || (!kMCPIDCheck && tparticle->PdgCode() == fgkPdgCode[kPDGAntiProton])) {
	  if (!CustomTrackCut(*trackq, 2)) continue;
	  PNegArray.push_back(qTracks);
	  PNegCounter++;
	}	
	if ((kMCPIDCheck && TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kPion)) <= kPIDSelecCut)
	    || (!kMCPIDCheck && tparticle->PdgCode() == fgkPdgCode[kPDGPionMinus])) {
	  if (!CustomTrackCut(*trackq, 3)) continue;
	  PiNegArray.push_back(qTracks);
	  PiNegCounter++;
	  PiNegSecArray.push_back(qTracks);
	  PiNegSecCounter++;
	}
      }
    }
  }
  if (!He4PosCounter && !He4NegCounter && !He3PosCounter && !He3NegCounter) return;
  if (!PPosCounter   && !PNegCounter)   return;
  if (!PiPosCounter  && !PiNegCounter)  return;
  // _________________________________________________ // 
  cout << "...finished track finding..." << endl;
  cout << "...found " << He4PosCounter << " helium4 tracks and " << He4NegCounter << " anti-helium4 tracks..." << endl;
  cout << "...found " << He3PosCounter << " helium tracks and " << He3NegCounter << " anti-helium tracks..." << endl;
  cout << "...found " << PPosCounter << " proton tracks and " << PNegCounter << " anti-proton tracks..." << endl;
  cout << "...found " << PiNegCounter << " pion minus tracks and " << PiPosCounter << " pion plus tracks..." << endl;
  cout << "...going to positive hypernuclei..." << endl;
  //__ positive __ //
  AliAnalysisTaskDoubleHypNucTreeLS::CalcPosDaughterHypNuc();
  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
  //__ negative __ //
  cout << "...going to negative hypernuclei..." << endl;
  AliAnalysisTaskDoubleHypNucTreeLS::CalcNegDaughterHypNuc();
  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
  // __ reset __ //
  ResetVals("Event");
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTreeLS::CalcPosDaughterHypNuc() {
  particle1->SetXYZM(0., 0., 0., 0.);
  particle2->SetXYZM(0., 0., 0., 0.);
  particle3->SetXYZM(0., 0., 0., 0.);
  particle4->SetXYZM(0., 0., 0., 0.);
  lorentzsum->SetXYZM(0., 0., 0., 0.);
  sublorentzsum->SetXYZM(0., 0., 0., 0.);
  cout << "...Analyzing pos 4LHe..." << endl;
  // ____________________________ first track _______________________________ //
  for (const Int_t& He3PosTracks : He3PosArray) {

    track1 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(He3PosTracks)); 
    exTrack1->CopyFromVTrack(track1);
    // __ MC part __ //
    if (fMCtrue && konlySig) {
            
      label1 = track1->GetLabel();
      labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));
            
      ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());
      if (ParticleMother1->PdgCode() != fgkPdgCode[kPDGHyperHelium4] && ParticleMother1->PdgCode() != fgkPdgCode[kPDGHyperHelium4Star]) continue;
      if (ParticleMother1) delete ParticleMother1;
    }
    // __ set pxpypz __ //
    particle1->SetXYZM(0., 0., 0., 0.);
    particle1->SetXYZM(2. * track1->Px(), 2. * track1->Py(), 2. * track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
    // ____________________________ second track _______________________________ //
    for (const Int_t& PPosTracks : PPosArray) {

      // __ reject using same track twice __ //
      if (He3PosTracks == PPosTracks) continue;

      track2 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(PPosTracks));
      exTrack2->CopyFromVTrack(track2);
      // __ MC part __ //
      if (fMCtrue && konlySig) {
                
	label2 = track2->GetLabel();
	labelMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label2));
                
	if (labelMother1 != labelMother2) continue;
                
	ParticleMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother2))->Particle());
	if (ParticleMother2->PdgCode() != fgkPdgCode[kPDGHyperHelium4] && ParticleMother2->PdgCode() != fgkPdgCode[kPDGHyperHelium4Star]) continue;
	if (ParticleMother2) delete ParticleMother2;
      }
      // __ dca rejection __ //     
      if (exTrack1->GetDCA(exTrack2, kMagF, xthiss, xpp) > kDCATracksCut) continue;
      // __ set pxpypz __ //
      particle2->SetXYZM(0., 0., 0., 0.);
      particle2->SetXYZM(track2->Px(), track2->Py(), track2->Pz(), AliPID::ParticleMass(AliPID::kProton));
      // ____________________________ third track _______________________________ //
      for (const Int_t& PiNegTracks : PiNegArray) {

	// __ reject using same track twice __ //
	if (PPosTracks == PiNegTracks || He3PosTracks == PiNegTracks) continue;

	track3 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(PiNegTracks));
	exTrack3->CopyFromVTrack(track3);
	// __ MC part __ //
	if (fMCtrue && konlySig) {
                    
	  label3 = track3->GetLabel();
	  labelMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label3));
                    
	  if (labelMother3 != labelMother2) continue;
                    
	  ParticleMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	  if (ParticleMother3->PdgCode() != fgkPdgCode[kPDGHyperHelium4] && ParticleMother3->PdgCode() != fgkPdgCode[kPDGHyperHelium4Star]) continue;
	  if (ParticleMother3) delete ParticleMother3;
	}
	// __ dca rejection  __ //
	if (exTrack1->GetDCA(exTrack2, kMagF, xthiss, xpp) > kDCATracksCut) continue;
	if (exTrack1->GetDCA(exTrack3, kMagF, xthiss, xpp) > kDCATracksCut) continue;
	if (exTrack2->GetDCA(exTrack3, kMagF, xthiss, xpp) > kDCATracksCut) continue;
	// __ set pxpypz __ //
	particle3->SetXYZM(0., 0., 0., 0.);
	particle3->SetXYZM(track3->Px(), track3->Py(), track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
	//______ daughter hypernucleus mass _______ //				
	sublorentzsum->SetXYZM(0., 0., 0., 0.);
	*sublorentzsum = *particle1 + *particle2 + *particle3;
	if (sublorentzsum->M() > kMax4LHeMass || sublorentzsum->M() < kMin4LHeMass) {
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  continue;
	}
	// _________________________________________________ //
	// __ MC part __ //
	if (fMCtrue) {
	  ParticleMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	  ParticleMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother2))->Particle());
	  ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());
                    
	  if (((ParticleMother1->PdgCode() == fgkPdgCode[kPDGHyperHelium4]
		&& ParticleMother2->PdgCode() == fgkPdgCode[kPDGHyperHelium4]
		&& ParticleMother3->PdgCode() == fgkPdgCode[kPDGHyperHelium4])
	       || (ParticleMother1->PdgCode() == fgkPdgCode[kPDGHyperHelium4Star]
		   && ParticleMother2->PdgCode() == fgkPdgCode[kPDGHyperHelium4Star]
		   && ParticleMother3->PdgCode() == fgkPdgCode[kPDGHyperHelium4Star]))
	      && labelMother1 == labelMother2
	      && labelMother2 == labelMother3) {
	    if (TMath::Abs(label1) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGHelium3]))
		&& TMath::Abs(label2) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGProton]))
		&& TMath::Abs(label3) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGPionMinus]))) {
	      fmctruth = 1;
	      fPrimary4LHe = mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(labelMother1)) == 0;
	      fisExcited = ParticleMother1->PdgCode() == fgkPdgCode[kPDGHyperHelium4Star];
	    }
	  }
	  if (ParticleMother1) delete ParticleMother1;
	  if (ParticleMother2) delete ParticleMother2;
	  if (ParticleMother3) delete ParticleMother3;
	}
	// _________________________________________________ //
	// __ vertex reconstruction with KF method __ //
	int status_KFreco = 0;
	if (kKFReco) status_KFreco = AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("4LHe", 1, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe), 8 = 4LH 2B

	if (status_KFreco) {
	  fParticleID = TMath::Abs(fEventID) + He3PosTracks + PPosTracks + PiNegTracks;
	  if (fMCtrue) fParticleID = fParticleID * (labelMother1 - label1);
	  fmSubMotherCheck = fmSubMother;
	  gTreeKF->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  // __ check daughter combinations __ //
	  AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("3LH", 3, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe), 8 = 4LH 2B
	  iTreeKF->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("3LH", 7, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe), 8 = 4LH 2B
	  iTreeKF->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	}
	// _________________________________________________ //
	// __ standard vertex reconstruction __ //
	int status_standardreco = 0;
	if (kStandardReco) status_standardreco = AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("4LHe", 1, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe), 8 = 4LH 2B

	if (status_standardreco) {
	  fParticleID = TMath::Abs(fEventID) + He3PosTracks + PPosTracks + PiNegTracks;
	  if (fMCtrue) fParticleID = fParticleID * (labelMother1 - label1);
	  fmSubMotherCheck = fmSubMother;
	  gTree->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  // __ check daughter combinations __ //
	  AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("3LH", 3, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe), 8 = 4LH 2B
	  iTree->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("3LH", 7, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe), 8 = 4LH 2B
	  iTree->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	}
	// _________________________________________________ //
	// __ reset __ //
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	fmSubMotherCheck = -99;
	fmMotherCheck = -99;
	fParticleID = -1;
	// _________________________________________________ //
	// __ go to 4LLH __ //
	cout << "...positive 4LHe..." << endl;
	cout << "...going to positive 4LLH..." << endl;
	AliAnalysisTaskDoubleHypNucTreeLS::CalcPosMotherHypNuc(He3PosTracks, PPosTracks, PiNegTracks, 1);
      } //track3
    } //track2
  } //track1
  // _____________________________________________________________________________________________________________________ //
  if (!kSkip3LHAnalysis) {
    cout << "...Analyzing pos 3LH..." << endl;
    particle1->SetXYZM(0., 0., 0., 0.);
    particle2->SetXYZM(0., 0., 0., 0.);
    particle3->SetXYZM(0., 0., 0., 0.);
    particle4->SetXYZM(0., 0., 0., 0.);
    lorentzsum->SetXYZM(0., 0., 0., 0.);
    sublorentzsum->SetXYZM(0., 0., 0., 0.);
    // ____________________________ first track _______________________________ //
    for (const Int_t& He3PosTracks : He3PosArray) {

      track1 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(He3PosTracks));
      exTrack1->CopyFromVTrack(track1);
      // __ MC part __ //
      if (fMCtrue && konlySig) {

	label1 = track1->GetLabel();
	labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));

	ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());
	if (ParticleMother1->PdgCode() != fgkPdgCode[kPDGHyperHydrogen3]) continue;
	if (ParticleMother1) delete ParticleMother1;
      }
      // __ set pxpypz __ //
      particle1->SetXYZM(0., 0., 0., 0.);
      particle1->SetXYZM(2. * track1->Px(), 2. * track1->Py(), 2. * track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
      // ____________________________ third track _______________________________ //
      for (const Int_t& PiNegTracks : PiNegArray) {

	// __ reject using same track twice __ //
	if (He3PosTracks == PiNegTracks) continue;

	track3 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(PiNegTracks));
	exTrack3->CopyFromVTrack(track3);
	// __ MC part __ //
	if (fMCtrue && konlySig) {

	  label3 = track3->GetLabel();
	  labelMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label3));

	  if (labelMother1 != labelMother3) continue;

	  ParticleMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	  if (ParticleMother3->PdgCode() != fgkPdgCode[kPDGHyperHydrogen3]) continue;
	  if (ParticleMother3) delete ParticleMother3;
	}
	// __ dca rejection __ //
	if (exTrack1->GetDCA(exTrack3, kMagF, xthiss, xpp) > kDCATracksCut) continue;
	// __ set pxpypz __ //
	particle3->SetXYZM(0., 0., 0., 0.);
	particle3->SetXYZM(track3->Px(), track3->Py(), track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
	//______ daughter hypernucleus mass _______ //			
	sublorentzsum->SetXYZM(0., 0., 0., 0.);
	*sublorentzsum = *particle1 + *particle3;
	if (sublorentzsum->M() > kMax3LHMass || sublorentzsum->M() < kMin3LHMass) {
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  continue;
	}
	// _________________________________________________ //
	// __ MC part __ //
	if (fMCtrue) {
	  ParticleMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	  ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());
	  if (ParticleMother1->PdgCode() == fgkPdgCode[kPDGHyperHydrogen3]
	      && ParticleMother3->PdgCode() == fgkPdgCode[kPDGHyperHydrogen3]
	      && (labelMother1 == labelMother3)) {
	    if (TMath::Abs(label1) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGHelium3]))
		&& TMath::Abs(label3) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGPionMinus]))) {
	      fmctruth = 1;
	      fPrimary4LHe = mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(labelMother1)) == 0;
	    }
	  }
	  if (ParticleMother1) delete ParticleMother1;
	  if (ParticleMother3) delete ParticleMother3;
	}
	// _________________________________________________ //
	// __ vertex reconstruction with KF method __ //
	int status_KFreco = 0;
	if (kKFReco) status_KFreco = AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("3LH", 2, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe), 8 = 4LH 2B
	if (status_KFreco) {
	  fParticleID = TMath::Abs(fEventID) + He3PosTracks + PiNegTracks;
	  if (fMCtrue) fParticleID = fParticleID * (labelMother1 - label1);
	  hTreeKF->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	}
	// _________________________________________________ //
	// __ standard vertex reconstruction __ //
	int status_standardreco = 0;
	if (kStandardReco) status_standardreco = AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("3LH", 2, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe), 8 = 4LH 2B
	if (status_standardreco) {
	  fParticleID = TMath::Abs(fEventID) + He3PosTracks + PiNegTracks;
	  if (fMCtrue) fParticleID = fParticleID * (labelMother1 - label1);
	  hTree->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	}
	// _________________________________________________ //
	// __ reset __ //
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	fmSubMotherCheck = -99;
	fmMotherCheck = -99;
	fParticleID = -1;
	// _________________________________________________ //
	// __ go to 4LLH __ //
	cout << "...positive 3LH..." << endl;
	cout << "...going to positive 4LLH..." << endl;
	AliAnalysisTaskDoubleHypNucTreeLS::CalcPosMotherHypNuc(He3PosTracks, PiNegTracks, 0, 2);
      }//track3
    }//track1
  }//kSkip3LH
  // _____________________________________________________________________________________________________________________ //
  if (!kSkip4LHAnalysis) {
    cout << "...Analyzing pos 4LH..." << endl;
    particle1->SetXYZM(0., 0., 0., 0.);
    particle2->SetXYZM(0., 0., 0., 0.);
    particle3->SetXYZM(0., 0., 0., 0.);
    particle4->SetXYZM(0., 0., 0., 0.);
    lorentzsum->SetXYZM(0., 0., 0., 0.);
    sublorentzsum->SetXYZM(0., 0., 0., 0.);
    // ____________________________ first track _______________________________ //
    for (const Int_t& He4PosTracks : He4PosArray) {

      track1 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(He4PosTracks));
      exTrack1->CopyFromVTrack(track1);
      // __ MC part __ //
      if (fMCtrue && konlySig) {
                
	label1 = track1->GetLabel();
	labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));
                
	ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());
	if (ParticleMother1->PdgCode() != fgkPdgCode[kPDGHyperHydrogen4] && ParticleMother1->PdgCode() != fgkPdgCode[kPDGHyperHydrogen4Star]) continue;
	if (ParticleMother1) delete ParticleMother1;
      }
      // __ set pxpypz __ //
      particle1->SetXYZM(0., 0., 0., 0.);
      particle1->SetXYZM(2. * track1->Px(), 2. * track1->Py(), 2. * track1->Pz(), AliPID::ParticleMass(AliPID::kAlpha));
      // ____________________________ third track _______________________________ //
      for (const Int_t& PiNegTracks : PiNegArray) {

	// __ reject using same track twice __ //
	if (He4PosTracks == PiNegTracks) continue;

	track3 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(PiNegTracks));
	exTrack3->CopyFromVTrack(track3);
	// __ MC part __ //
	if (fMCtrue && konlySig) {
                    
	  label3 = track3->GetLabel();
	  labelMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label3));
                    
	  if (labelMother1 != labelMother3) continue;
                    
	  ParticleMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	  if (ParticleMother3->PdgCode() != fgkPdgCode[kPDGHyperHydrogen4] && ParticleMother3->PdgCode() != fgkPdgCode[kPDGHyperHydrogen4Star]) continue;
	  if (ParticleMother3) delete ParticleMother3;
	}
	// __ dca rejection __ //
	if (exTrack1->GetDCA(exTrack3, kMagF, xthiss, xpp) > kDCATracksCut) continue;
	// __ set pxpypz __ //
	particle3->SetXYZM(0., 0., 0., 0.);
	particle3->SetXYZM(track3->Px(), track3->Py(), track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
	//______ daughter hypernucleus mass _______ //			
	sublorentzsum->SetXYZM(0., 0., 0., 0.);
	*sublorentzsum = *particle1 + *particle3;
	if (sublorentzsum->M() > kMax4LHMass || sublorentzsum->M() < kMin4LHMass) {
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  continue;
	}
	// _________________________________________________ //
	// __ MC part __ //
	if (fMCtrue) {
	  ParticleMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	  ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());
	  if (((ParticleMother1->PdgCode() == fgkPdgCode[kPDGHyperHydrogen4]
		&& ParticleMother3->PdgCode() == fgkPdgCode[kPDGHyperHydrogen4])
	       || (ParticleMother1->PdgCode() == fgkPdgCode[kPDGHyperHydrogen4Star]
		   && ParticleMother3->PdgCode() == fgkPdgCode[kPDGHyperHydrogen4Star]))
	      && labelMother1 == labelMother3) {
	    if (TMath::Abs(label1) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGHelium4]))
		&& TMath::Abs(label3) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGPionMinus]))) {
	      fmctruth = 1;
	      fPrimary4LHe = mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(labelMother1)) == 0;
	      fisExcited = ParticleMother1->PdgCode() == fgkPdgCode[kPDGHyperHydrogen4Star];
	    }
	  }
	  if (ParticleMother1) delete ParticleMother1;
	  if (ParticleMother3) delete ParticleMother3;
	}
	// _________________________________________________ //
	// __ vertex reconstruction with KF method __ //
	int status_KFreco = 0;
	if (kKFReco) status_KFreco = AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("4LH", 8, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe), 8 = 4LH 2B
	if (status_KFreco) {
	  fParticleID = TMath::Abs(fEventID) + He4PosTracks + PiNegTracks;
	  if (fMCtrue) fParticleID = fParticleID * (labelMother1 - label1);
	  hTreeKF->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	}
	// _________________________________________________ //
	// __ standard vertex reconstruction __ //
	int status_standardreco = 0;
	if (kStandardReco) status_standardreco = AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("4LH", 8, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe), 8 = 4LH 2B 
	if (status_standardreco) {
	  fParticleID = TMath::Abs(fEventID) + He4PosTracks + PiNegTracks;
	  if (fMCtrue) fParticleID = fParticleID * (labelMother1 - label1);
	  hTree->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	}
	// _________________________________________________ //
	// __ reset __ //
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	fmSubMotherCheck = -99;
	fmMotherCheck = -99;
	fParticleID = -1;
	// _________________________________________________ //				
      }//track3
    }//track1   
  }//kSkip4LH
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTreeLS::CalcNegDaughterHypNuc() {
  particle1->SetXYZM(0., 0., 0., 0.);
  particle2->SetXYZM(0., 0., 0., 0.);
  particle3->SetXYZM(0., 0., 0., 0.);
  particle4->SetXYZM(0., 0., 0., 0.);
  lorentzsum->SetXYZM(0., 0., 0., 0.);
  sublorentzsum->SetXYZM(0., 0., 0., 0.);
  cout << "...Analyzing neg 4LHe..." << endl;
  // ____________________________ first track _______________________________ //
  for (const Int_t& He3NegTracks : He3NegArray) {

    track1 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(He3NegTracks));
    exTrack1->CopyFromVTrack(track1);
    // __ MC part __ //
    if (fMCtrue && konlySig) {
            
      label1 = track1->GetLabel();
      labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));
            
      ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());
      if (ParticleMother1->PdgCode() != fgkPdgCode[kPDGAntiHyperHelium4] && ParticleMother1->PdgCode() != fgkPdgCode[kPDGAntiHyperHelium4Star]) continue;
      if (ParticleMother1) delete ParticleMother1;
    }
    // __ set pxpypz __ //
    particle1->SetXYZM(0., 0., 0., 0.);
    particle1->SetXYZM(2. * track1->Px(), 2. * track1->Py(), 2. * track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
    // ____________________________ second track _______________________________ //
    for (const Int_t& PNegTracks : PNegArray) {

      // __ reject using same track twice __ //
      if (He3NegTracks == PNegTracks) continue;

      track2 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(PNegTracks));
      exTrack2->CopyFromVTrack(track2);
      // __ MC part __ //
      if (fMCtrue && konlySig) {
                
	label2 = track2->GetLabel();
	labelMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label2));
                
	if (labelMother1 != labelMother2) continue;
                
	ParticleMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother2))->Particle());
	if (ParticleMother2->PdgCode() != fgkPdgCode[kPDGAntiHyperHelium4]  && ParticleMother2->PdgCode() != fgkPdgCode[kPDGAntiHyperHelium4Star]) continue;
	if (ParticleMother2) delete ParticleMother2;
      }
      // __ dca rejection __ //      
      if (exTrack1->GetDCA(exTrack2, kMagF, xthiss, xpp) > kDCATracksCut) continue;
      // __ set pxpypz __ //
      particle2->SetXYZM(0., 0., 0., 0.);
      particle2->SetXYZM(track2->Px(), track2->Py(), track2->Pz(), AliPID::ParticleMass(AliPID::kProton));
      // ____________________________ third track _______________________________ //
      for (const Int_t& PiPosTracks : PiPosArray) {

	// __ reject using same track twice __ //
	if (PNegTracks == PiPosTracks || He3NegTracks == PiPosTracks) continue;

	track3 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(PiPosTracks));
	exTrack3->CopyFromVTrack(track3);
	// __ MC part __ //
	if (fMCtrue && konlySig) {
                    
	  label3 = track3->GetLabel();
	  labelMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label3));
                    
	  if (labelMother3 != labelMother2) continue;
                    
	  ParticleMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	  if (ParticleMother3->PdgCode() != fgkPdgCode[kPDGAntiHyperHelium4]  && ParticleMother3->PdgCode() != fgkPdgCode[kPDGAntiHyperHelium4Star]) continue;
	  if (ParticleMother3) delete ParticleMother3;
	}
	// __ dca rejection __ //
	if (exTrack1->GetDCA(exTrack2, kMagF, xthiss, xpp) > kDCATracksCut) continue;
	if (exTrack1->GetDCA(exTrack3, kMagF, xthiss, xpp) > kDCATracksCut) continue;
	if (exTrack2->GetDCA(exTrack3, kMagF, xthiss, xpp) > kDCATracksCut) continue;
	// __ set pxpypz __ //
	particle3->SetXYZM(0., 0., 0., 0.);
	particle3->SetXYZM(track3->Px(), track3->Py(), track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
	//______ daughter hypernucleus mass _______ //				
	sublorentzsum->SetXYZM(0., 0., 0., 0.);
	*sublorentzsum = *particle1 + *particle2 + *particle3;
	if (sublorentzsum->M() > kMax4LHeMass || sublorentzsum->M() < kMin4LHeMass) {
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  continue;
	}
	// _________________________________________________ //	
	// __ MC part __ //
	if (fMCtrue) {
	  ParticleMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	  ParticleMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother2))->Particle());
	  ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());
                    
	  if (((ParticleMother1->PdgCode() == fgkPdgCode[kPDGAntiHyperHelium4]
		&& ParticleMother2->PdgCode() == fgkPdgCode[kPDGAntiHyperHelium4]
		&& ParticleMother3->PdgCode() == fgkPdgCode[kPDGAntiHyperHelium4])
	       || (ParticleMother1->PdgCode() == fgkPdgCode[kPDGAntiHyperHelium4Star]
		   && ParticleMother2->PdgCode() == fgkPdgCode[kPDGAntiHyperHelium4Star]
		   && ParticleMother3->PdgCode() == fgkPdgCode[kPDGAntiHyperHelium4Star]))
	      && labelMother1 == labelMother2
	      && labelMother2 == labelMother3) {
	    if (TMath::Abs(label1) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGAntiHelium3]))
		&& TMath::Abs(label2) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGAntiProton]))
		&& TMath::Abs(label3) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGPionPlus]))) {
	      fmctruth = 1;
	      fPrimary4LHe = mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(labelMother1)) == 0;
	      fisExcited = ParticleMother1->PdgCode() == fgkPdgCode[kPDGAntiHyperHelium4Star];
	    }
	  }
	  if (ParticleMother1) delete ParticleMother1;
	  if (ParticleMother2) delete ParticleMother2;
	  if (ParticleMother3) delete ParticleMother3;
	}
	// _________________________________________________ //
	// __ vertex reconstruction with KF method __ //
	int status_KFreco = 0;
	if (kKFReco) status_KFreco = AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("4LHe", 1, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe), 8 = 4LH 2B
	if (status_KFreco) {
	  fParticleID = TMath::Abs(fEventID) + He3NegTracks + PNegTracks + PiPosTracks;
	  if (fMCtrue) fParticleID = fParticleID * (labelMother1 - label1);
	  fmSubMotherCheck = fmSubMother;
	  gTreeKF->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  // __ check for combinations __ //
	  AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("3LH", 3, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe), 8 = 4LH 2B
	  iTreeKF->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("3LH", 7, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe), 8 = 4LH 2B
	  iTreeKF->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	}
	// _________________________________________________ //
	// __ standard vertex reconstruction __ //
	int status_standardreco = 0;
	if (kStandardReco) status_standardreco = AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("4LHe", 1, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe), 8 = 4LH 2B
	if (status_standardreco) {
	  fParticleID = TMath::Abs(fEventID) + He3NegTracks + PNegTracks + PiPosTracks;
	  if (fMCtrue) fParticleID = fParticleID * (labelMother1 - label1);
	  fmSubMotherCheck = fmSubMother;
	  gTree->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  // __ check for combinations __ //
	  AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("3LH", 3, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe), 8 = 4LH 2B
	  iTree->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("3LH", 7, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe), 8 = 4LH 2B
	  iTree->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	}
	// _________________________________________________ //
	// __ reset __ //
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	fmSubMotherCheck = -99;
	fmMotherCheck = -99;
	fParticleID = -1;
	// _________________________________________________ //
	// __ go to 4LLH __ //
	cout << "...negative 4LHe..." << endl;
	cout << "...going to negative 4LLH..." << endl;
	AliAnalysisTaskDoubleHypNucTreeLS::CalcNegMotherHypNuc(He3NegTracks, PNegTracks, PiPosTracks, 1);
      } //track3
    } //track2
  } //track1
  // _________________________________________________________________________ //
  if (!kSkip3LHAnalysis) {
    cout << "...Analyzing neg 3LH..." << endl;
    particle1->SetXYZM(0., 0., 0., 0.);
    particle2->SetXYZM(0., 0., 0., 0.);
    particle3->SetXYZM(0., 0., 0., 0.);
    particle4->SetXYZM(0., 0., 0., 0.);
    lorentzsum->SetXYZM(0., 0., 0., 0.);
    sublorentzsum->SetXYZM(0., 0., 0., 0.);
    // ____________________________ first track _______________________________ //
    for (const Int_t& He3NegTracks : He3NegArray) {

      track1 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(He3NegTracks));
      exTrack1->CopyFromVTrack(track1);
      // __ MC part __ //
      if (fMCtrue && konlySig) {

	label1 = track1->GetLabel();
	labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));

	ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());
	if (ParticleMother1->PdgCode() != fgkPdgCode[kPDGAntiHyperHydrogen3]) continue;
	if (ParticleMother1) delete ParticleMother1;
      }
      // __ set pxpypz __ //
      particle1->SetXYZM(0., 0., 0., 0.);
      particle1->SetXYZM(2. * track1->Px(), 2. * track1->Py(), 2. * track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
      // ____________________________ third track _______________________________ //
      for (const Int_t& PiPosTracks : PiPosArray) {

	// __ reject using same track twice __ //
	if (He3NegTracks == PiPosTracks) continue;

	track3 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(PiPosTracks));
	exTrack3->CopyFromVTrack(track3);
	// __ MC part __ //
	if (fMCtrue && konlySig) {

	  label3 = track3->GetLabel();
	  labelMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label3));

	  if (labelMother3 != labelMother1) continue;

	  ParticleMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	  if (ParticleMother3->PdgCode() != fgkPdgCode[kPDGAntiHyperHydrogen3]) continue;
	  if (ParticleMother3) delete ParticleMother3;
	}
	// __ dca rejection __ //
	if (exTrack1->GetDCA(exTrack3, kMagF, xthiss, xpp) > kDCATracksCut) continue;
	// __ set pxpypz __ //
	particle3->SetXYZM(0., 0., 0., 0.);
	particle3->SetXYZM(track3->Px(), track3->Py(), track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
	//______ daughter hypernucleus mass _______ //			
	sublorentzsum->SetXYZM(0., 0., 0., 0.);
	*sublorentzsum = *particle1 + *particle3;
	if ((sublorentzsum->M() > kMax3LHMass || sublorentzsum->M() < kMin3LHMass)) {
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  continue;
	}
	// _________________________________________________ //     
	// __ MC part __ //
	if (fMCtrue) {
	  ParticleMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	  ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());
	  if (ParticleMother1->PdgCode() == fgkPdgCode[kPDGAntiHyperHydrogen3]
	      && ParticleMother3->PdgCode() == fgkPdgCode[kPDGAntiHyperHydrogen3]
	      && labelMother1 == labelMother3) {
	    if (TMath::Abs(label1) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGAntiHelium3]))
		&& TMath::Abs(label3) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGPionPlus]))) {
	      fmctruth = 1;
	      fPrimary4LHe = mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(labelMother1)) == 0;
	    }
	  }
	  if (ParticleMother1) delete ParticleMother1;
	  if (ParticleMother3) delete ParticleMother3;
	}
	// _________________________________________________ //
	// __ vertex reconstruction with KF method __ //
	int status_KFreco = 0;
	if (kKFReco) status_KFreco = AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("3LH", 2, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe), 8 = 4LH 2B
	if (status_KFreco) {
	  fParticleID = TMath::Abs(fEventID) + He3NegTracks + PiPosTracks;
	  if (fMCtrue) fParticleID = fParticleID * (labelMother1 - label1);
	  hTreeKF->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	}
	// _________________________________________________ //
	// __ standard vertex reconstruction __ //
	int status_standardreco = 0;
	if (kStandardReco) status_standardreco = AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("3LH", 2, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe), 8 = 4LH 2B
	if (status_standardreco) {
	  fParticleID = TMath::Abs(fEventID) + He3NegTracks + PiPosTracks;
	  if (fMCtrue) fParticleID = fParticleID * (labelMother1 - label1);
	  hTree->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	}
	// _________________________________________________ //
	// __ reset __ //
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	fmSubMotherCheck = -99;
	fmMotherCheck = -99;
	fParticleID = -1;
	// _________________________________________________ //
	// __ go to 4LLH __ //
	cout << "...negative 3LH..." << endl;
	cout << "...going to negative 4LLH..." << endl;
	AliAnalysisTaskDoubleHypNucTreeLS::CalcNegMotherHypNuc(He3NegTracks, PiPosTracks, 0, 2);
      }//track3
    }//track1
  }//kSkip3LH
  // _________________________________________________________________________ //
  if (!kSkip4LHAnalysis) {
    cout << "...Analyzing neg 4LH..." << endl;
    particle1->SetXYZM(0., 0., 0., 0.);
    particle2->SetXYZM(0., 0., 0., 0.);
    particle3->SetXYZM(0., 0., 0., 0.);
    particle4->SetXYZM(0., 0., 0., 0.);
    lorentzsum->SetXYZM(0., 0., 0., 0.);
    sublorentzsum->SetXYZM(0., 0., 0., 0.);
    // ____________________________ first track _______________________________ //
    for (const Int_t& He4NegTracks : He4NegArray) {

      track1 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(He4NegTracks));
      exTrack1->CopyFromVTrack(track1);
      // __ MC part __ //
      if (fMCtrue && konlySig) {
                
	label1 = track1->GetLabel();
	labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));
                
	ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());
	if (ParticleMother1->PdgCode() != fgkPdgCode[kPDGAntiHyperHydrogen4] && ParticleMother1->PdgCode() != fgkPdgCode[kPDGAntiHyperHydrogen4Star]) continue;
	if (ParticleMother1) delete ParticleMother1;
      }
      // __ set pxpypz __ //
      particle1->SetXYZM(0., 0., 0., 0.);
      particle1->SetXYZM(2. * track1->Px(), 2. * track1->Py(), 2. * track1->Pz(), AliPID::ParticleMass(AliPID::kAlpha));
      // ____________________________ third track _______________________________ //
      for (const Int_t& PiPosTracks : PiPosArray) {

	// __ reject using same track twice __ //
	if (He4NegTracks == PiPosTracks) continue;

	track3 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(PiPosTracks));
	exTrack3->CopyFromVTrack(track3);
	// __ MC part __ //
	if (fMCtrue && konlySig) {
                    
	  label3 = track3->GetLabel();
	  labelMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label3));
                    
	  if (labelMother3 != labelMother1) continue;
                    
	  ParticleMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	  if (ParticleMother3->PdgCode() != fgkPdgCode[kPDGAntiHyperHydrogen4] && ParticleMother3->PdgCode() != fgkPdgCode[kPDGAntiHyperHydrogen4Star]) continue;
	  if (ParticleMother3) delete ParticleMother3;
	}
	// __ dca rejection __ //
	if (exTrack1->GetDCA(exTrack3, kMagF, xthiss, xpp) > kDCATracksCut) continue;
	// __ set pxpypz __ //
	particle3->SetXYZM(0., 0., 0., 0.);
	particle3->SetXYZM(track3->Px(), track3->Py(), track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
	//______ daughter hypernucleus mass _______ //			
	sublorentzsum->SetXYZM(0., 0., 0., 0.);
	*sublorentzsum = *particle1 + *particle3;
	if ((sublorentzsum->M() > kMax4LHMass || sublorentzsum->M() < kMin4LHMass)) {
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  continue;
	}
	// _________________________________________________ //     
	// __ MC part __ //
	if (fMCtrue) {
	  ParticleMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	  ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());
	  if (((ParticleMother1->PdgCode() == fgkPdgCode[kPDGAntiHyperHydrogen4]
		&& ParticleMother3->PdgCode() == fgkPdgCode[kPDGAntiHyperHydrogen4])
	       || (ParticleMother1->PdgCode() == fgkPdgCode[kPDGAntiHyperHydrogen4Star]
		   && ParticleMother3->PdgCode() == fgkPdgCode[kPDGAntiHyperHydrogen4Star]))
	      && labelMother1 == labelMother3) {
	    if (TMath::Abs(label1) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGAntiHelium4]))
		&& TMath::Abs(label3) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGPionPlus]))) {
	      fmctruth = 1;
	      fPrimary4LHe = mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(labelMother1)) == 0;
	      fisExcited = ParticleMother1->PdgCode() == fgkPdgCode[kPDGAntiHyperHydrogen4Star];
	    }
	  }
	  if (ParticleMother1) delete ParticleMother1;
	  if (ParticleMother3) delete ParticleMother3;
	}
	// _________________________________________________ //
	// __ vertex reconstruction with KF method __ //
	int status_KFreco = 0;
	if (kKFReco) status_KFreco = AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("4LH", 8, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe), 8 = 4LH 2B
	if (status_KFreco) {
	  fParticleID = TMath::Abs(fEventID) + He4NegTracks + PiPosTracks;
	  if (fMCtrue) fParticleID = fParticleID * (labelMother1 - label1);
	  hTreeKF->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	}
	// _________________________________________________ //
	// __ standard vertex reconstruction __ //
	int status_standardreco = 0;
	if (kStandardReco) status_standardreco = AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("4LH", 8, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe), 8 = 4LH 2B
	if (status_standardreco) {
	  fParticleID = TMath::Abs(fEventID) + He4NegTracks + PiPosTracks;
	  if (fMCtrue) fParticleID = fParticleID * (labelMother1 - label1);
	  hTree->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	}
	// _________________________________________________ //
	// __ reset __ //
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	fmSubMotherCheck = -99;
	fmMotherCheck = -99;
	fParticleID = -1;
      }//track3
    }//track1 
  }//kSkip4LH
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTreeLS::CalcPosMotherHypNuc(Int_t Track1Entry, Int_t Track2Entry, Int_t Track3Entry, Int_t kDecayChannel) {
  Int_t kCharge = 1;
  particle1->SetXYZM(0., 0., 0., 0.);
  particle2->SetXYZM(0., 0., 0., 0.);
  particle3->SetXYZM(0., 0., 0., 0.);
  particle4->SetXYZM(0., 0., 0., 0.);
  lorentzsum->SetXYZM(0., 0., 0., 0.);
  sublorentzsum->SetXYZM(0., 0., 0., 0.);
  // _________________________________________________ //     
  if (kDecayChannel == 1) {
    if (!Track1Entry) return; if (!Track2Entry) return; if (!Track3Entry) return;
    track1 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(Track1Entry));
    exTrack1->CopyFromVTrack(track1);
    track2 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(Track2Entry));
    exTrack2->CopyFromVTrack(track2);
    track3 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(Track3Entry));
    exTrack3->CopyFromVTrack(track3);
    // __ set pxpypz __ //
    particle1->SetXYZM(2. * track1->Px(), 2. * track1->Py(), 2. * track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
    particle2->SetXYZM(track2->Px(), track2->Py(), track2->Pz(), AliPID::ParticleMass(AliPID::kProton));
    particle3->SetXYZM(track3->Px(), track3->Py(), track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
    // __ MC part __ //
    if (fMCtrue && konlySig) {

      label1 = track1->GetLabel();
      labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));
      labelGrandMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother1));

      ParticleGrandMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother1))->Particle());
      if (ParticleGrandMother1->PdgCode() != fgkPdgCode[kPDGDoubleHyperHydrogen4]) return;
      if (ParticleGrandMother1) delete ParticleGrandMother1;

      label2 = track2->GetLabel();
      labelMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label2));
      labelGrandMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother2));

      ParticleGrandMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother2))->Particle());
      if (ParticleGrandMother2->PdgCode() != fgkPdgCode[kPDGDoubleHyperHydrogen4]) return;
      if (ParticleGrandMother2) delete ParticleGrandMother2;

      label3 = track3->GetLabel();
      labelMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label3));
      labelGrandMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother3));

      ParticleGrandMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother3))->Particle());
      if (ParticleGrandMother3->PdgCode() != fgkPdgCode[kPDGDoubleHyperHydrogen4]) return;
      if (ParticleGrandMother3) delete ParticleGrandMother3;
    }
    // ____________________________ fourth track _______________________________ //
    for (const Int_t& PiNegSecTracks : PiNegSecArray) {

      // __ reject using same track twice __ //
      if (Track2Entry == PiNegSecTracks || Track1Entry == PiNegSecTracks || Track3Entry == PiNegSecTracks) continue;

      track4 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(PiNegSecTracks));
      exTrack4->CopyFromVTrack(track4);
      // __ MC part __ //
      if (fMCtrue && konlySig) {

	label4 = track4->GetLabel();
	labelMother4 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label4));

	if (labelGrandMother3 != labelMother4) continue;

	ParticleMother4 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother4))->Particle());
	if (ParticleMother4->PdgCode() != fgkPdgCode[kPDGDoubleHyperHydrogen4]) continue;
	if (ParticleMother4) delete ParticleMother4;
      }
      // __ set pxpypz __ //
      particle4->SetXYZM(0., 0., 0., 0.);
      particle4->SetXYZM(track4->Px(), track4->Py(), track4->Pz(), AliPID::ParticleMass(AliPID::kPion));
      // _________________________________________________ //      
      lorentzsum->SetXYZM(0., 0., 0., 0.);
      *lorentzsum = *particle1 + *particle2 + *particle3 + *particle4;
      // _________________________________________________ //
      // __ select mother hypernucleus by inv mass __ //
      if (lorentzsum->M() > kMax4LLHMass || lorentzsum->M() < kMin4LLHMass) {
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	continue;
      }
      // _________________________________________________ //
      // __ MC part __ //
      if (fMCtrue) {

	ParticleMother4 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother4))->Particle());
	ParticleMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	ParticleMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother2))->Particle());
	ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());

	ParticleGrandMother4 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother4))->Particle());
	ParticleGrandMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother3))->Particle());
	ParticleGrandMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother2))->Particle());
	ParticleGrandMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother1))->Particle());

	if (ParticleGrandMother1->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	    && ParticleGrandMother2->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	    && ParticleGrandMother3->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	    && ParticleMother4->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]) {

	  if (labelGrandMother1 == labelGrandMother2 && labelGrandMother2 == labelGrandMother3 && labelGrandMother3 == labelMother4) {

	    if (TMath::Abs(label1) == TMath::Abs(GetLabel(labelGrandMother1, kCharge * fgkPdgCode[kPDGHelium3], kCharge * fgkPdgCode[kPDGHyperHelium4]))
		&& TMath::Abs(label2) == TMath::Abs(GetLabel(labelGrandMother1, kCharge * fgkPdgCode[kPDGProton], kCharge * fgkPdgCode[kPDGHyperHelium4]))
		&& TMath::Abs(label3) == TMath::Abs(GetLabel(labelGrandMother1, kCharge * fgkPdgCode[kPDGPionMinus], kCharge * fgkPdgCode[kPDGHyperHelium4]))
		&& TMath::Abs(label4) == TMath::Abs(GetLabel(labelGrandMother1, kCharge * fgkPdgCode[kPDGPionMinus]))) {
	      fmctruth = 1;
	      fPrimary4LHe = mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(labelMother1)) == 0;
	    }
	  }
	}
	if (ParticleGrandMother1) delete ParticleGrandMother1;
	if (ParticleGrandMother2) delete ParticleGrandMother2;
	if (ParticleGrandMother3) delete ParticleGrandMother3;
	if (ParticleGrandMother4) delete ParticleGrandMother4;
	if (ParticleMother1) delete ParticleMother1;
	if (ParticleMother2) delete ParticleMother2;
	if (ParticleMother3) delete ParticleMother3;
	if (ParticleMother4) delete ParticleMother4;
      }
      // _________________________________________________ //
      // __ vertex reconstruction with KF method __ //
      int status_KFreco = 0;
      if (kKFReco) status_KFreco = AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("4LLH", kDecayChannel, kCharge);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
      if (status_KFreco) {
	fParticleID = TMath::Abs(fEventID) + Track1Entry + Track2Entry + Track3Entry + PiNegSecTracks;
	if (fMCtrue) fParticleID = fParticleID * (labelMother1 - label1);
	fmMotherCheck = fmMother;
	fmSubMotherCheck = fmSubMother;
	fTreeKF->Fill();
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("3LH", 3, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	iTreeKF->Fill();
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("3LH", 4, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	iTreeKF->Fill();
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("3LH", 5, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	iTreeKF->Fill();
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("3LH", 6, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	iTreeKF->Fill();
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
      }
      // _________________________________________________ //
      // __ standard vertex reconstruction __ //
      int status_standardreco = 0;
      if (kStandardReco) status_standardreco = AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("4LLH", kDecayChannel, kCharge);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
      if (status_standardreco) {
	fParticleID = TMath::Abs(fEventID) + Track1Entry + Track2Entry + Track3Entry + PiNegSecTracks;
	if (fMCtrue) fParticleID = fParticleID * (labelMother1 - label1);
	fmMotherCheck = fmMother;
	fmSubMotherCheck = fmSubMother;
	fTree->Fill();
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("3LH", 3, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	iTree->Fill();
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("3LH", 4, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	iTree->Fill();
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("3LH", 5, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	iTree->Fill();
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("3LH", 6, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	iTree->Fill();
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
      }
      // __ reset __ //
      AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
      fmMotherCheck = -99;
      fmSubMotherCheck = -99;
      fParticleID = -1;
      // _________________________________________________ //
    }//track4
  }//kDecayChannel
  // ____________________________________________________________________________________________________________________________________ //
  if (kDecayChannel == 2) {
    if (!Track1Entry) return; if (!Track2Entry) return;
    track1 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(Track1Entry));
    exTrack1->CopyFromVTrack(track1);
    track3 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(Track2Entry));
    exTrack3->CopyFromVTrack(track3);
    // __ set pxpypz __ //
    particle1->SetXYZM(2. * track1->Px(), 2. * track1->Py(), 2. * track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
    particle3->SetXYZM(track3->Px(), track3->Py(), track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
    // __ MC part __ //
    if (fMCtrue && konlySig) {

      label1 = track1->GetLabel();
      labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));
      labelGrandMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother1));

      ParticleGrandMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother1))->Particle());
      if (ParticleGrandMother1->PdgCode() != fgkPdgCode[kPDGDoubleHyperHydrogen4]) return;
      if (ParticleGrandMother1) delete ParticleGrandMother1;

      label3 = track3->GetLabel();
      labelMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label3));
      labelGrandMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother3));

      ParticleGrandMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother3))->Particle());
      if (ParticleGrandMother3->PdgCode() != fgkPdgCode[kPDGDoubleHyperHydrogen4]) return;
      if (ParticleGrandMother3) delete ParticleGrandMother3;
    }
    // _________________________________________________ //
    for (const Int_t& PPosTracks : PPosArray) {
      // __ reject using same track twice __ //
      if (Track1Entry == PPosTracks || Track2Entry == PPosTracks) continue;

      track2 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(PPosTracks));
      exTrack2->CopyFromVTrack(track2);
      // __ MC part __ //
      if (fMCtrue && konlySig) {

	label2 = track2->GetLabel();
	labelMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label2));
	labelGrandMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother2));

	if (labelGrandMother3 != labelMother2) continue;

	ParticleMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother2))->Particle());
	if (ParticleMother2->PdgCode() != fgkPdgCode[kPDGDoubleHyperHydrogen4]) continue;
	if (ParticleMother2) delete ParticleMother2;
      }
      // __ set pxpypz __ //
      particle2->SetXYZM(0., 0., 0., 0.);
      particle2->SetXYZM(track2->Px(), track2->Py(), track2->Pz(), AliPID::ParticleMass(AliPID::kProton));
      // ____________________________ fourth track _______________________________ //
      for (const Int_t& PiNegSecTracks : PiNegSecArray) {

	// __ reject using same track twice __ //
	if (PPosTracks == PiNegSecTracks || Track1Entry == PiNegSecTracks || Track3Entry == PiNegSecTracks) continue;

	track4 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(PiNegSecTracks));
	exTrack4->CopyFromVTrack(track4);
	// __ MC part __ //
	if (fMCtrue && konlySig) {

	  label4 = track4->GetLabel();
	  labelMother4 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label4));
	  labelGrandMother4 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother4));

	  if (labelMother2 != labelMother4) continue;

	  ParticleMother4 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother4))->Particle());
	  if (ParticleMother4->PdgCode() != fgkPdgCode[kPDGDoubleHyperHydrogen4]) continue;
	  if (ParticleMother4) delete ParticleMother4;
	}
	// __ dca rejection __ //
	if (exTrack2->GetDCA(exTrack4, kMagF, xthiss, xpp) > kDCATracksCut) continue;
	// __ set pxpypz __ //
	particle4->SetXYZM(0., 0., 0., 0.);
	particle4->SetXYZM(track4->Px(), track4->Py(), track4->Pz(), AliPID::ParticleMass(AliPID::kPion));
	// _________________________________________________ //	  	 	
	lorentzsum->SetXYZM(0., 0., 0., 0.);
	*lorentzsum = *particle1 + *particle2 + *particle3 + *particle4;
	// _________________________________________________ //
	// __ select mother hypernucleus by inv mass __ //
	if (lorentzsum->M() > kMax4LLHMass || lorentzsum->M() < kMin4LLHMass) {
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  continue;
	}
	// _________________________________________________ //
	// __ MC part __ //
	if (fMCtrue) {

	  ParticleMother4 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother4))->Particle());
	  ParticleMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	  ParticleMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother2))->Particle());
	  ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());

	  ParticleGrandMother4 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother4))->Particle());
	  ParticleGrandMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother3))->Particle());
	  ParticleGrandMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother2))->Particle());
	  ParticleGrandMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother1))->Particle());

	  if (ParticleGrandMother1->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	      && ParticleMother2->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	      && ParticleGrandMother3->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	      && ParticleMother4->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]) {

	    if (labelGrandMother1 == labelMother2 && labelMother2 == labelGrandMother3 && labelGrandMother3 == labelMother4) {

	      if (TMath::Abs(label1) == TMath::Abs(GetLabel(labelGrandMother1, kCharge * fgkPdgCode[kPDGHelium3], kCharge * fgkPdgCode[kPDGHyperHydrogen3]))
		  && TMath::Abs(label3) == TMath::Abs(GetLabel(labelGrandMother1, kCharge * fgkPdgCode[kPDGPionMinus], kCharge * fgkPdgCode[kPDGHyperHydrogen3]))
		  && TMath::Abs(label2) == TMath::Abs(GetLabel(labelGrandMother1, kCharge * fgkPdgCode[kPDGProton]))
		  && TMath::Abs(label4) == TMath::Abs(GetLabel(labelGrandMother1, kCharge * fgkPdgCode[kPDGPionMinus]))) {
		fmctruth = 1;
		fPrimary4LHe = mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(labelMother1)) == 0;
	      }
	    }
	  }
	  if (ParticleGrandMother1) delete ParticleGrandMother1;
	  if (ParticleGrandMother2) delete ParticleGrandMother2;
	  if (ParticleGrandMother3) delete ParticleGrandMother3;
	  if (ParticleGrandMother4) delete ParticleGrandMother4;
	  if (ParticleMother1) delete ParticleMother1;
	  if (ParticleMother2) delete ParticleMother2;
	  if (ParticleMother3) delete ParticleMother3;
	  if (ParticleMother4) delete ParticleMother4;
	}
	// _________________________________________________ //
	// __ vertex reconstruction with KF method __ //
	int status_KFreco = 0;
	if (kKFReco) status_KFreco = AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("4LLH", kDecayChannel, kCharge);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	if (status_KFreco) {
	  fParticleID = TMath::Abs(fEventID) + Track1Entry + PPosTracks + Track3Entry + PiNegSecTracks;
	  if (fMCtrue) fParticleID = fParticleID * (labelMother1 - label1);
	  fmMotherCheck = fmMother;
	  fmSubMotherCheck = fmSubMother;
	  fTreeKF->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("3LH", 4, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	  iTreeKF->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("3LH", 5, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	  iTreeKF->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("3LH", 6, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	  iTreeKF->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	}
	// _________________________________________________ //
	// __ standard vertex reconstruction __ //
	int status_standardreco = 0;
	if (kStandardReco) status_standardreco = AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("4LLH", kDecayChannel, kCharge);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	if (status_standardreco) {
	  fParticleID = TMath::Abs(fEventID) + Track1Entry + PPosTracks + Track3Entry + PiNegSecTracks;
	  if (fMCtrue) fParticleID = fParticleID * (labelMother1 - label1);
	  fmMotherCheck = fmMother;
	  fmSubMotherCheck = fmSubMother;
	  fTree->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("3LH", 4, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	  iTree->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("3LH", 5, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	  iTree->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("3LH", 6, 1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	  iTree->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	}
	// __ reset __ //
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	fmMotherCheck = -99;
	fmSubMotherCheck = -99;
	fParticleID = -1;
	// _________________________________________________ //
      }//track4
    }//track2
  }//kDecayChannel
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTreeLS::CalcNegMotherHypNuc(Int_t Track1Entry, Int_t Track2Entry, Int_t Track3Entry, Int_t kDecayChannel) {
  Int_t kCharge = -1;
  particle1->SetXYZM(0., 0., 0., 0.);
  particle2->SetXYZM(0., 0., 0., 0.);
  particle3->SetXYZM(0., 0., 0., 0.);
  particle4->SetXYZM(0., 0., 0., 0.);
  lorentzsum->SetXYZM(0., 0., 0., 0.);
  sublorentzsum->SetXYZM(0., 0., 0., 0.);
  // _________________________________________________ //
  if (kDecayChannel == 1) {
    if (!Track1Entry) return; if (!Track2Entry) return; if (!Track3Entry) return;
    track1 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(Track1Entry));
    exTrack1->CopyFromVTrack(track1);
    track2 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(Track2Entry));
    exTrack2->CopyFromVTrack(track2);
    track3 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(Track3Entry));
    exTrack3->CopyFromVTrack(track3);
    // __ set pxpypz __ //
    particle1->SetXYZM(2. * track1->Px(), 2. * track1->Py(), 2. * track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
    particle2->SetXYZM(track2->Px(), track2->Py(), track2->Pz(), AliPID::ParticleMass(AliPID::kProton));
    particle3->SetXYZM(track3->Px(), track3->Py(), track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
    // __ MC part __ //
    if (fMCtrue && konlySig) {

      label1 = track1->GetLabel();
      labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));
      labelGrandMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother1));

      ParticleGrandMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother1))->Particle());
      if (ParticleGrandMother1->PdgCode() != kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]) return;
      if (ParticleGrandMother1) delete ParticleGrandMother1;

      label2 = track2->GetLabel();
      labelMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label2));
      labelGrandMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother2));

      ParticleGrandMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother2))->Particle());
      if (ParticleGrandMother2->PdgCode() != kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]) return;
      if (ParticleGrandMother2) delete ParticleGrandMother2;

      label3 = track3->GetLabel();
      labelMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label3));
      labelGrandMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother3));

      ParticleGrandMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother3))->Particle());
      if (ParticleGrandMother3->PdgCode() != kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]) return;
      if (ParticleGrandMother3) delete ParticleGrandMother3;
    }
    // ____________________________ fourth track _______________________________ //
    for (const Int_t& PiPosSecTracks : PiPosSecArray) {

      // __ reject using same track twice __ //
      if (Track2Entry == PiPosSecTracks || Track1Entry == PiPosSecTracks || Track3Entry == PiPosSecTracks) continue;

      track4 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(PiPosSecTracks));
      exTrack4->CopyFromVTrack(track4);
      // __ MC part __ //
      if (fMCtrue && konlySig) {

	label4 = track4->GetLabel();
	labelMother4 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label4));

	if (labelGrandMother3 != labelMother4) continue;

	ParticleMother4 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother4))->Particle());
	if (ParticleMother4->PdgCode() != kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]) continue;
	if (ParticleMother4) delete ParticleMother4;
      }
      // __ set pxpypz __ //
      particle4->SetXYZM(0., 0., 0., 0.);
      particle4->SetXYZM(track4->Px(), track4->Py(), track4->Pz(), AliPID::ParticleMass(AliPID::kPion));
      // _________________________________________________ //	  	       
      lorentzsum->SetXYZM(0., 0., 0., 0.);
      *lorentzsum = *particle1 + *particle2 + *particle3 + *particle4;
      // _________________________________________________ //
      // __ select mother hypernucleus by inv mass __ //
      if (lorentzsum->M() > kMax4LLHMass || lorentzsum->M() < kMin4LLHMass) {
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	continue;
      }
      // _________________________________________________ //
      // __ MC part __ //
      if (fMCtrue) {

	ParticleMother4 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother4))->Particle());
	ParticleMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	ParticleMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother2))->Particle());
	ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());

	ParticleGrandMother4 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother4))->Particle());
	ParticleGrandMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother3))->Particle());
	ParticleGrandMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother2))->Particle());
	ParticleGrandMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother1))->Particle());

	if (ParticleGrandMother1->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	    && ParticleGrandMother2->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	    && ParticleGrandMother3->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	    && ParticleMother4->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]) {

	  if (labelGrandMother1 == labelGrandMother2 && labelGrandMother2 == labelGrandMother3 && labelGrandMother3 == labelMother4) {

	    if (TMath::Abs(label1) == TMath::Abs(GetLabel(labelGrandMother1, kCharge * fgkPdgCode[kPDGHelium3], kCharge * fgkPdgCode[kPDGHyperHelium4]))
		&& TMath::Abs(label2) == TMath::Abs(GetLabel(labelGrandMother1, kCharge * fgkPdgCode[kPDGProton], kCharge * fgkPdgCode[kPDGHyperHelium4]))
		&& TMath::Abs(label3) == TMath::Abs(GetLabel(labelGrandMother1, kCharge * fgkPdgCode[kPDGPionMinus], kCharge * fgkPdgCode[kPDGHyperHelium4]))
		&& TMath::Abs(label4) == TMath::Abs(GetLabel(labelGrandMother1, kCharge * fgkPdgCode[kPDGPionMinus]))) {
	      fmctruth = 1;
	      fPrimary4LHe = mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(labelMother1)) == 0;
	    }
	  }
	}
	if (ParticleGrandMother1) delete ParticleGrandMother1;
	if (ParticleGrandMother2) delete ParticleGrandMother2;
	if (ParticleGrandMother3) delete ParticleGrandMother3;
	if (ParticleGrandMother4) delete ParticleGrandMother4;
	if (ParticleMother1) delete ParticleMother1;
	if (ParticleMother2) delete ParticleMother2;
	if (ParticleMother3) delete ParticleMother3;
	if (ParticleMother4) delete ParticleMother4;
      }
      // _________________________________________________ //
      // __ vertex reconstruction with KF method __ //
      int status_KFreco = 0;
      if (kKFReco) status_KFreco = AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("4LLH", kDecayChannel, kCharge);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
      if (status_KFreco) {
	fParticleID = TMath::Abs(fEventID) + Track1Entry + Track2Entry + Track3Entry + PiPosSecTracks;
	if (fMCtrue) fParticleID = fParticleID * (labelMother1 - label1);
	fmMotherCheck = fmMother;
	fmSubMotherCheck = fmSubMother;
	fTreeKF->Fill();
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("3LH", 3, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	iTreeKF->Fill();
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("3LH", 4, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	iTreeKF->Fill();
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("3LH", 5, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	iTreeKF->Fill();
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("3LH", 6, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	iTreeKF->Fill();
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
      }
      // _________________________________________________ //
      // __ standard vertex reconstruction __ //
      int status_standardreco = 0;
      if (kStandardReco) status_standardreco = AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("4LLH", kDecayChannel, kCharge);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
      if (status_standardreco) {
	fParticleID = TMath::Abs(fEventID) + Track1Entry + Track2Entry + Track3Entry + PiPosSecTracks;
	if (fMCtrue) fParticleID = fParticleID * (labelMother1 - label1);
	fmMotherCheck = fmMother;
	fmSubMotherCheck = fmSubMother;
	fTree->Fill();
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("3LH", 3, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	iTree->Fill();
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("3LH", 4, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	iTree->Fill();
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("3LH", 5, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	iTree->Fill();
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("3LH", 6, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	iTree->Fill();
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
      }
      // __ reset __ //
      AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
      fmMotherCheck = -99;
      fmSubMotherCheck = -99;
      fParticleID = -1;
      // _________________________________________________ //
    }//track4
  }//kDecayChannel
  // _____________________________________________________________________ //
  if (kDecayChannel == 2) {
    if (!Track1Entry) return; if (!Track2Entry) return;
    track1 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(Track1Entry));
    exTrack1->CopyFromVTrack(track1);
    track3 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(Track2Entry));
    exTrack3->CopyFromVTrack(track3);
    // __ set pxpypz __ //
    particle1->SetXYZM(2. * track1->Px(), 2. * track1->Py(), 2. * track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
    particle3->SetXYZM(track3->Px(), track3->Py(), track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
    // __ MC part __ //
    if (fMCtrue && konlySig) {

      label1 = track1->GetLabel();
      labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));
      labelGrandMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother1));

      ParticleGrandMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother1))->Particle());
      if (ParticleGrandMother1->PdgCode() != kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]) return;
      if (ParticleGrandMother1) delete ParticleGrandMother1;

      label3 = track3->GetLabel();
      labelMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label3));
      labelGrandMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother3));

      ParticleGrandMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother3))->Particle());
      if (ParticleGrandMother3->PdgCode() != kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]) return;
      if (ParticleGrandMother3) delete ParticleGrandMother3;
    }
    for (const Int_t& PNegTracks : PNegArray) {
      // __ reject using same track twice __ //
      if (Track1Entry == PNegTracks || Track2Entry == PNegTracks) continue;

      track2 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(PNegTracks));
      exTrack2->CopyFromVTrack(track2);
      // __ MC part __ //
      if (fMCtrue && konlySig) {

	label2 = track2->GetLabel();
	labelMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label2));

	if (labelGrandMother3 != labelMother2) continue;

	ParticleMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother2))->Particle());
	if (ParticleMother2->PdgCode() != kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]) continue;
	if (ParticleMother2) delete ParticleMother2;
      }
      // __ set pxpypz __ //
      particle2->SetXYZM(0., 0., 0., 0.);
      particle2->SetXYZM(track2->Px(), track2->Py(), track2->Pz(), AliPID::ParticleMass(AliPID::kProton));
      // ____________________________ fourth track _______________________________ //
      for (const Int_t& PiPosSecTracks : PiPosSecArray) {

	// __ reject using same track twice __ //
	if (PNegTracks == PiPosSecTracks || Track1Entry == PiPosSecTracks || Track2Entry == PiPosSecTracks) continue;

	track4 = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(PiPosSecTracks));
	exTrack4->CopyFromVTrack(track4);
	// __ MC part __ //
	if (fMCtrue && konlySig) {

	  label4 = track4->GetLabel();
	  labelMother4 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label4));

	  if (labelMother2 != labelMother4) continue;

	  ParticleMother4 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother4))->Particle());
	  if (ParticleMother4->PdgCode() != kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]) continue;
	  if (ParticleMother4) delete ParticleMother4;
	}
	// __ dca rejection __ //
	if (exTrack2->GetDCA(exTrack4, kMagF, xthiss, xpp) > kDCATracksCut) continue;
	// __ set pxpypz __ //
	particle4->SetXYZM(0., 0., 0., 0.);
	particle4->SetXYZM(track4->Px(), track4->Py(), track4->Pz(), AliPID::ParticleMass(AliPID::kPion));
	// _________________________________________________ //	  	 	
	lorentzsum->SetXYZM(0., 0., 0., 0.);
	*lorentzsum = *particle1 + *particle2 + *particle3 + *particle4;
	// _________________________________________________ //
	// __ select mother hypernucleus by inv mass __ //
	if (lorentzsum->M() > kMax4LLHMass || lorentzsum->M() < kMin4LLHMass) {
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  continue;
	}
	// _________________________________________________ //
	// __ MC part __ //
	if (fMCtrue) {

	  ParticleMother4 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother4))->Particle());
	  ParticleMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	  ParticleMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother2))->Particle());
	  ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());

	  ParticleGrandMother4 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother4))->Particle());
	  ParticleGrandMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother3))->Particle());
	  ParticleGrandMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother2))->Particle());
	  ParticleGrandMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelGrandMother1))->Particle());

	  if (ParticleGrandMother1->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	      && ParticleMother2->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	      && ParticleGrandMother3->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	      && ParticleMother4->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]) {

	    if (labelGrandMother1 == labelMother2 && labelMother2 == labelGrandMother3 && labelGrandMother3 == labelMother4) {

	      if (TMath::Abs(label1) == TMath::Abs(GetLabel(labelGrandMother1, kCharge * fgkPdgCode[kPDGHelium3], kCharge * fgkPdgCode[kPDGHyperHydrogen3]))
		  && TMath::Abs(label3) == TMath::Abs(GetLabel(labelGrandMother1, kCharge * fgkPdgCode[kPDGPionMinus], kCharge * fgkPdgCode[kPDGHyperHydrogen3]))
		  && TMath::Abs(label2) == TMath::Abs(GetLabel(labelGrandMother1, kCharge * fgkPdgCode[kPDGProton]))
		  && TMath::Abs(label4) == TMath::Abs(GetLabel(labelGrandMother1, kCharge * fgkPdgCode[kPDGPionMinus]))) {
		fmctruth = 1;
		fPrimary4LHe = mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(labelMother1)) == 0;
	      }
	    }
	  }
	  if (ParticleGrandMother1) delete ParticleGrandMother1;
	  if (ParticleGrandMother2) delete ParticleGrandMother2;
	  if (ParticleGrandMother3) delete ParticleGrandMother3;
	  if (ParticleGrandMother4) delete ParticleGrandMother4;
	  if (ParticleMother1) delete ParticleMother1;
	  if (ParticleMother2) delete ParticleMother2;
	  if (ParticleMother3) delete ParticleMother3;
	  if (ParticleMother4) delete ParticleMother4;
	}
	// _________________________________________________ //
	// __ vertex reconstruction with KF method __ //
	int status_KFreco = 0;
	if (kKFReco) status_KFreco = AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("4LLH", kDecayChannel, kCharge);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	if (status_KFreco) {
	  fParticleID = TMath::Abs(fEventID) + Track1Entry + PNegTracks + Track3Entry + PiPosSecTracks;
	  if (fMCtrue) fParticleID = fParticleID * (labelMother1 - label1);
	  fmMotherCheck = fmMother;
	  fmSubMotherCheck = fmSubMother;
	  fTreeKF->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("3LH", 4, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	  iTreeKF->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("3LH", 5, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	  iTreeKF->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction("3LH", 6, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	  iTreeKF->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	}
	// _________________________________________________ //
	// __ standard vertex reconstruction __ //
	int status_standardreco = 0;
	if (kStandardReco) status_standardreco = AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("4LLH", kDecayChannel, kCharge);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	if (status_standardreco) {
	  fParticleID = TMath::Abs(fEventID) + Track1Entry + PNegTracks + Track3Entry + PiPosSecTracks;
	  if (fMCtrue) fParticleID = fParticleID * (labelMother1 - label1);
	  fmMotherCheck = fmMother;
	  fmSubMotherCheck = fmSubMother;
	  fTree->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("3LH", 4, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	  iTree->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("3LH", 5, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	  iTree->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	  AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction("3LH", 6, -1);//1 = 4LHe, 2 = 3LH, 3 = 3LH (tert check, 4LHe), 4 = 3LH (sec check), 5 = Lambda (tert check), 6 = Lambda (sec check), 7 = Lambda (4LHe)
	  iTree->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	}
	// __ reset __ //
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	fmMotherCheck = -99;
	fmSubMotherCheck = -99;
	fParticleID = -1;
	// _________________________________________________ //
      }//track4
    }//track2
  }//kDecayChannel
}
// _________________________________________________ //
Int_t AliAnalysisTaskDoubleHypNucTreeLS::StandardReconstruction(TString Mother, Int_t kDecayChannel, Int_t ksign) {
  Double_t impar[2];
  UShort_t idThree[3] = { 0, 1, 2 };
  UShort_t idTwo[2] = { 0, 1 };
  Int_t    kFinderAlgo = 1; //1: StrLinVertexFinderMinDist(UseWeights=1) (default), 2: StrLinVertexFinderMinDist(UseWeights=0), 3: HelixVertexFinder(), 4: VertexFinder(1), 5: VertexFinder(0)
  // _________________________________ //
  // __ get prim vtx __ //
  primVertex = AliAnalysisTaskDoubleHypNucTreeLS::AODToESDVertex(*vertex);
  PrimVertex[0] = primVertex->GetX();
  PrimVertex[1] = primVertex->GetY();
  PrimVertex[2] = primVertex->GetZ();
  // __ coord of prim vtx __ //
  fPrimVertexX = PrimVertex[0];
  fPrimVertexY = PrimVertex[1];
  fPrimVertexZ = PrimVertex[2];
  fPrimVertChi2 = primVertex->GetChi2();
  fPrimVertNDF = primVertex->GetNDF();
  // _________________________________ //      
  if (Mother == "3LH") {
    exTrack->Reset(); exTrack1->Reset(); exTrack2->Reset(); exTrack3->Reset(); exTrack4->Reset();
    if (kDecayChannel == 2 || kDecayChannel == 3) {
      exTrack1->CopyFromVTrack(track1);
      exTrack3->CopyFromVTrack(track3);
    }
    if (kDecayChannel == 4) {
      exTrack1->CopyFromVTrack(track1);
      exTrack3->CopyFromVTrack(track4);
    }
    if (kDecayChannel == 5 || kDecayChannel == 7) {
      exTrack1->CopyFromVTrack(track2);
      exTrack3->CopyFromVTrack(track3);
    }
    if (kDecayChannel == 6) {
      exTrack1->CopyFromVTrack(track2);
      exTrack3->CopyFromVTrack(track4);
    }
    // _________________________________ //
    // __ init AliVertexerTracks __ //
    tertvertexer = new AliVertexerTracks(fAODevent->GetMagneticField());
    tertvertexer->SetFinderAlgorithm(kFinderAlgo);
    // __start at prim vtx __ //
    tertvertexer->SetVtxStart(primVertex);
    if (primVertex)   delete primVertex;
    // __ tert vtx __ //
    trkArray = new TObjArray(2);
    trkArray->AddAt(exTrack1, 0);
    trkArray->AddAt(exTrack3, 1);
    tertVertex = (AliESDVertex*)tertvertexer->VertexForSelectedTracks(trkArray, idTwo, kTRUE, kTRUE, kFALSE);
    if (tertvertexer) delete tertvertexer;
    if (trkArray) delete trkArray;
    // _________________________________ //
    // __ vertex quality __ //
    fTertVertChi2     = tertVertex->GetChi2();
    fTertVertNDF      = tertVertex->GetNDF();
    // _________________________________________________ //
    // __ coord of tert vtx __ //
    TertVertex[0]     = tertVertex->GetX();
    TertVertex[1]     = tertVertex->GetY();
    TertVertex[2]     = tertVertex->GetZ();
    fTertVertexX      = TertVertex[0];
    fTertVertexY      = TertVertex[1];
    fTertVertexZ      = TertVertex[2];
    // _________________________________ //
    // __ propagate tracks and get impact parameters __ //    
    fPropDCADaughter  = (exTrack1->PropagateToDCA(tertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter    = impar[0];
    fImParzDaughter   = impar[1];
    fDcaSecDaughter   = TMath::Abs(exTrack1->GetD(TertVertex[0], TertVertex[1], kMagF));
    fPropDCADaughter2 = (exTrack3->PropagateToDCA(tertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter2   = impar[0];
    fImParzDaughter2  = impar[1];
    fDcaSecDaughter2  = TMath::Abs(exTrack3->GetD(TertVertex[0], TertVertex[1], kMagF));
    if (tertVertex) delete tertVertex;
    // _________________________________________________ //
    // __ dca between tracks __ //
    fDCA3B1 = TMath::Abs(exTrack3->GetDCA(exTrack1, kMagF, xthiss, xpp));
    // _________________________________________________ //
    if (kDecayChannel < 5) particle1->SetXYZM(2. * exTrack1->Px(), 2. * exTrack1->Py(), 2. * exTrack1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
    else                   particle1->SetXYZM(exTrack1->Px(), exTrack1->Py(), exTrack1->Pz(), AliPID::ParticleMass(AliPID::kProton));
    particle3->SetXYZM(exTrack3->Px(), exTrack3->Py(), exTrack3->Pz(), AliPID::ParticleMass(AliPID::kPion));
    sublorentzsum->SetXYZM(0., 0., 0., 0.);
    *sublorentzsum = *particle1 + *particle3;
    // _________________________________________________ //
    // __ daughter hypernucleus cos(PA) to prim vtx __ //
    dd[0] = PrimVertex[0] - TertVertex[0];
    dd[1] = PrimVertex[1] - TertVertex[1];
    dd[2] = PrimVertex[2] - TertVertex[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    fSubPA = TMath::Cos(sublorentzsum->Angle(*h));
    fctSubMother = (sublorentzsum->M() * TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2))) / sublorentzsum->P();
    // _________________________________________________ //
    if (kDecayChannel == 2 && fSubPA < kPointingAngleCut) return 0;
    // _________________________________________________ // 
    // __ cuts from 04/2020 __ //
    fDcaDaughtero = exTrack1->GetD(fPrimVertexX, fPrimVertexY, kMagF);
    fDcaDaughter2o = exTrack3->GetD(fPrimVertexX, fPrimVertexY, kMagF);
    // _________________________________________________ //
    // __ track information __ //
    if (kDecayChannel == 2) AliAnalysisTaskDoubleHypNucTreeLS::SetDaughterInformation("3LH", kDecayChannel);
    // _________________________________________________ //	  	
    // __ daughter hypernucleus information __ //
    if (kDecayChannel < 5) fPDGMother = ksign * fgkPdgCode[kPDGHyperHydrogen3];
    else                   fPDGMother = ksign * 3122;
    fChargeMother = ksign;
    fRecoMethod   = 1;
    fDecayChannel = kDecayChannel;
    fmSubMother   = sublorentzsum->M();
    fESubMother   = sublorentzsum->E();
    fpxSubMother  = sublorentzsum->Px();
    fpySubMother  = sublorentzsum->Py();
    fpzSubMother  = sublorentzsum->Pz();
    fptSubMother  = sublorentzsum->Pt();
    fpSubMother   = sublorentzsum->P();
    fySubMother   = sublorentzsum->Rapidity();

    fpxDaughter   = particle1->Px();
    fpyDaughter   = particle1->Py();
    fpzDaughter   = particle1->Pz();
    fptDaughter   = particle1->Pt();
    fEDaughter    = particle1->E();
    fyDaughter    = particle1->Rapidity();
    fpxDaughter2  = particle3->Px();
    fpyDaughter2  = particle3->Py();
    fpzDaughter2  = particle3->Pz();
    fptDaughter2  = particle3->Pt();
    fEDaughter2   = particle3->E();
    fyDaughter2   = particle3->Rapidity();
    return 1;
  }
  // _________________________________ //      
  else if (Mother == "4LH") {
    exTrack->Reset(); exTrack1->Reset(); exTrack2->Reset(); exTrack3->Reset(); exTrack4->Reset();
    exTrack1->CopyFromVTrack(track1);
    exTrack3->CopyFromVTrack(track3);
    // _________________________________ //
    // __ init AliVertexerTracks __ //
    tertvertexer = new AliVertexerTracks(fAODevent->GetMagneticField());
    tertvertexer->SetFinderAlgorithm(kFinderAlgo);
    // __start at prim vtx __ //
    tertvertexer->SetVtxStart(primVertex);
    if (primVertex)   delete primVertex;
    // __ tert vtx __ //
    trkArray = new TObjArray(2);
    trkArray->AddAt(exTrack1, 0);
    trkArray->AddAt(exTrack3, 1);
    tertVertex = (AliESDVertex*)tertvertexer->VertexForSelectedTracks(trkArray, idTwo, kTRUE, kTRUE, kFALSE);
    if (tertvertexer) delete tertvertexer;
    if (trkArray) delete trkArray;
    // _________________________________ //
    // __ vertex quality __ //
    fTertVertChi2    = tertVertex->GetChi2();
    fTertVertNDF     = tertVertex->GetNDF();
    // _________________________________________________ //
    // __ coord of tert vtx __ //
    TertVertex[0]    = tertVertex->GetX();
    TertVertex[1]    = tertVertex->GetY();
    TertVertex[2]    = tertVertex->GetZ();
    fTertVertexX     = TertVertex[0];
    fTertVertexY     = TertVertex[1];
    fTertVertexZ     = TertVertex[2];
    // _________________________________ //
    // __ propagate tracks and get impact parameters __ //    
    fPropDCADaughter = (exTrack1->PropagateToDCA(tertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter   = impar[0];
    fImParzDaughter  = impar[1];
    fDcaSecDaughter  = TMath::Abs(exTrack1->GetD(TertVertex[0], TertVertex[1], kMagF));
    fPropDCADaughter2= (exTrack3->PropagateToDCA(tertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter2  = impar[0];
    fImParzDaughter2 = impar[1];
    fDcaSecDaughter2 = TMath::Abs(exTrack3->GetD(TertVertex[0], TertVertex[1], kMagF));
    if (tertVertex) delete tertVertex;
    // _________________________________________________ //
    // __ dca between tracks __ //
    fDCA3B1 = TMath::Abs(exTrack3->GetDCA(exTrack1, kMagF, xthiss, xpp));
    // _________________________________________________ //
    particle1->SetXYZM(2. * exTrack1->Px(), 2. * exTrack1->Py(), 2. * exTrack1->Pz(), AliPID::ParticleMass(AliPID::kAlpha));
    particle3->SetXYZM(exTrack3->Px(), exTrack3->Py(), exTrack3->Pz(), AliPID::ParticleMass(AliPID::kPion));
    sublorentzsum->SetXYZM(0., 0., 0., 0.);
    *sublorentzsum = *particle1 + *particle3;
    // _________________________________________________ //
    // __ daughter hypernucleus cos(PA) to prim vtx __ //
    dd[0] = PrimVertex[0] - TertVertex[0];
    dd[1] = PrimVertex[1] - TertVertex[1];
    dd[2] = PrimVertex[2] - TertVertex[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    fSubPA = TMath::Cos(sublorentzsum->Angle(*h));
    fctSubMother = (sublorentzsum->M() * TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2))) / sublorentzsum->P();
    // _________________________________________________ //
    if (fSubPA < kPointingAngleCut) return 0;
    // _________________________________________________ // 
    // __ cuts from 04/2020 __ //
    fDcaDaughtero = exTrack1->GetD(fPrimVertexX, fPrimVertexY, kMagF);
    fDcaDaughter2o = exTrack3->GetD(fPrimVertexX, fPrimVertexY, kMagF);
    // _________________________________________________ //
    // __ track information __ //
    AliAnalysisTaskDoubleHypNucTreeLS::SetDaughterInformation("4LH", kDecayChannel);
    // _________________________________________________ //	  	
    // __ daughter hypernucleus information __ //
    fPDGMother    = ksign * fgkPdgCode[kPDGHyperHydrogen4];
    fChargeMother = ksign;
    fRecoMethod   = 1;
    fDecayChannel = kDecayChannel;
    fmSubMother   = sublorentzsum->M();
    fESubMother   = sublorentzsum->E();
    fpxSubMother  = sublorentzsum->Px();
    fpySubMother  = sublorentzsum->Py();
    fpzSubMother  = sublorentzsum->Pz();
    fptSubMother  = sublorentzsum->Pt();
    fpSubMother   = sublorentzsum->P();
    fySubMother   = sublorentzsum->Rapidity();

    fpxDaughter   = particle1->Px();
    fpyDaughter   = particle1->Py();
    fpzDaughter   = particle1->Pz();
    fptDaughter   = particle1->Pt();
    fEDaughter    = particle1->E();
    fyDaughter    = particle1->Rapidity();
    fpxDaughter2  = particle3->Px();
    fpyDaughter2  = particle3->Py();
    fpzDaughter2  = particle3->Pz();
    fptDaughter2  = particle3->Pt();
    fEDaughter2   = particle3->E();
    fyDaughter2   = particle3->Rapidity();
    return 1;
  }
  // _________________________________ //      
  else if (Mother == "4LHe") {
    exTrack->Reset(); exTrack1->Reset(); exTrack2->Reset(); exTrack3->Reset(); exTrack4->Reset();
    exTrack1->CopyFromVTrack(track1);
    exTrack2->CopyFromVTrack(track2);
    exTrack3->CopyFromVTrack(track3);
    // _________________________________ //
    // __ init AliVertexerTracks __ //
    tertvertexer = new AliVertexerTracks(fAODevent->GetMagneticField());
    tertvertexer->SetFinderAlgorithm(kFinderAlgo);
    // __start at prim vtx __ //
    tertvertexer->SetVtxStart(primVertex);
    if (primVertex)   delete primVertex;
    // __ tert vtx __ //
    trkArray = new TObjArray(3);
    trkArray->AddAt(exTrack1, 0);
    trkArray->AddAt(exTrack2, 1);
    trkArray->AddAt(exTrack3, 2);
    tertVertex = (AliESDVertex*)tertvertexer->VertexForSelectedTracks(trkArray, idThree, kTRUE, kTRUE, kFALSE);
    if (tertvertexer) delete tertvertexer;
    if (trkArray) delete trkArray;
    // _________________________________ //
    // __ vertex quality __ //
    fTertVertChi2    = tertVertex->GetChi2();
    fTertVertNDF     = tertVertex->GetNDF();
    // _________________________________________________ //
    // __ coord of tert vtx __ //
    TertVertex[0]    = tertVertex->GetX();
    TertVertex[1]    = tertVertex->GetY();
    TertVertex[2]    = tertVertex->GetZ();
    fTertVertexX     = TertVertex[0];
    fTertVertexY     = TertVertex[1];
    fTertVertexZ     = TertVertex[2];
    // _________________________________ //
    // __ propagate tracks and get impact parameters __ //
    fPropDCADaughter = (exTrack1->PropagateToDCA(tertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter   = impar[0];
    fImParzDaughter  = impar[1];
    fDcaSecDaughter  = TMath::Abs(exTrack1->GetD(TertVertex[0], TertVertex[1], kMagF));
    fPropDCADaughter1= (exTrack2->PropagateToDCA(tertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter1  = impar[0];
    fImParzDaughter1 = impar[1];
    fDcaSecDaughter1 = TMath::Abs(exTrack2->GetD(TertVertex[0], TertVertex[1], kMagF));
    fPropDCADaughter2= (exTrack3->PropagateToDCA(tertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter2  = impar[0];
    fImParzDaughter2 = impar[1];
    fDcaSecDaughter2 = TMath::Abs(exTrack3->GetD(TertVertex[0], TertVertex[1], kMagF));
    if (tertVertex) delete tertVertex;
    // _________________________________________________ //
    // __ dca between tracks __ //
    fDCA3B1 = TMath::Abs(exTrack2->GetDCA(exTrack1, kMagF, xthiss, xpp));
    fDCA3B2 = TMath::Abs(exTrack3->GetDCA(exTrack1, kMagF, xthiss, xpp));
    fDCA3B3 = TMath::Abs(exTrack3->GetDCA(exTrack2, kMagF, xthiss, xpp));
    // _________________________________________________ //
    particle1->SetXYZM(2. * exTrack1->Px(), 2. * exTrack1->Py(), 2. * exTrack1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
    particle2->SetXYZM(exTrack2->Px(), exTrack2->Py(), exTrack2->Pz(), AliPID::ParticleMass(AliPID::kProton));
    particle3->SetXYZM(exTrack3->Px(), exTrack3->Py(), exTrack3->Pz(), AliPID::ParticleMass(AliPID::kPion));
    sublorentzsum->SetXYZM(0., 0., 0., 0.);
    *sublorentzsum = *particle1 + *particle2 + *particle3;
    // _________________________________________________ //
    // __ daughter hypernucleus cos(PA) to prim vtx __ //
    dd[0] = PrimVertex[0] - TertVertex[0];
    dd[1] = PrimVertex[1] - TertVertex[1];
    dd[2] = PrimVertex[2] - TertVertex[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    fSubPA = TMath::Cos(sublorentzsum->Angle(*h));
    fctSubMother = (sublorentzsum->M() * TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2))) / sublorentzsum->P();
    // _________________________________________________ //
    if (fSubPA < kPointingAngleCut) return 0;
    // _________________________________________________ // 
    // __ cuts from 04/2020 __ //
    fDcaDaughtero = exTrack1->GetD(fPrimVertexX, fPrimVertexY, kMagF);
    fDcaDaughter1o = exTrack2->GetD(fPrimVertexX, fPrimVertexY, kMagF);
    fDcaDaughter2o = exTrack3->GetD(fPrimVertexX, fPrimVertexY, kMagF);
    // _________________________________________________ //
    // __ track information __ //
    AliAnalysisTaskDoubleHypNucTreeLS::SetDaughterInformation("4LHe", kDecayChannel);
    // _________________________________________________ //	  	
    // __ daughter hypernucleus information __ //
    fPDGMother    = ksign * fgkPdgCode[kPDGHyperHelium4];
    fChargeMother = ksign * 2;
    fRecoMethod   = 1;
    fDecayChannel = kDecayChannel;
    fmSubMother   = sublorentzsum->M();
    fESubMother   = sublorentzsum->E();
    fpxSubMother  = sublorentzsum->Px();
    fpySubMother  = sublorentzsum->Py();
    fpzSubMother  = sublorentzsum->Pz();
    fptSubMother  = sublorentzsum->Pt();
    fpSubMother   = sublorentzsum->P();
    fySubMother   = sublorentzsum->Rapidity();

    fpxDaughter   = particle1->Px();
    fpyDaughter   = particle1->Py();
    fpzDaughter   = particle1->Pz();
    fptDaughter   = particle1->Pt();
    fEDaughter    = particle1->E();
    fyDaughter    = particle1->Rapidity();
    fpxDaughter1  = particle2->Px();
    fpyDaughter1  = particle2->Py();
    fpzDaughter1  = particle2->Pz();
    fptDaughter1  = particle2->Pt();
    fEDaughter1   = particle2->E();
    fyDaughter1   = particle2->Rapidity();
    fpxDaughter2  = particle3->Px();
    fpyDaughter2  = particle3->Py();
    fpzDaughter2  = particle3->Pz();
    fptDaughter2  = particle3->Pt();
    fEDaughter2   = particle3->E();
    fyDaughter2   = particle3->Rapidity();
    return 1;
  }
  // _________________________________ //
  else if (Mother == "4LLH" && kDecayChannel == 1) {
    // ******************************************************************************
    // * Two different methods:                                                     *
    // * kVariante = 2 creates the tertiary vertex first                            *
    // * which coordinates are taken to reconstruct the daughter hypernucleus track *
    // * MC checks showed, that the Pointing resolution is worse in this case.      *
    // * kVariante = 1 uses the coordinates of the primary vertex to reconstruct    *
    // * the daughter hypernucleus track (default).                                 *
    // ******************************************************************************
    // _________________________________ //
    // __ init AliVertexerTracks __ //
    secvertexer = new AliVertexerTracks(fAODevent->GetMagneticField());
    tertvertexer = new AliVertexerTracks(fAODevent->GetMagneticField());
    secvertexer->SetFinderAlgorithm(kFinderAlgo);
    tertvertexer->SetFinderAlgorithm(kFinderAlgo);
    // _________________________________ //
    exTrack->Reset(); exTrack1->Reset(); exTrack2->Reset(); exTrack3->Reset(); exTrack4->Reset();
    exTrack1->CopyFromVTrack(track1);
    exTrack2->CopyFromVTrack(track2);
    exTrack3->CopyFromVTrack(track3);
    exTrack4->CopyFromVTrack(track4);
    // _________________________________ //
    if (kVariante == 2) {
      // __ tert vtx __ //
      trkArray = new TObjArray(3);
      trkArray->AddAt(exTrack1, 0);
      trkArray->AddAt(exTrack2, 1);
      trkArray->AddAt(exTrack3, 2);
      // __start at prim vtx __ //
      tertvertexer->SetVtxStart(primVertex);
      tertVertex = (AliESDVertex*)tertvertexer->VertexForSelectedTracks(trkArray, idThree, kTRUE, kTRUE, kFALSE);
      if (trkArray) delete trkArray;
      if (tertvertexer) delete tertvertexer;
      // __ coord of tert vtx __ //
      TertVertex[0]    = tertVertex->GetX();
      TertVertex[1]    = tertVertex->GetY();
      TertVertex[2]    = tertVertex->GetZ();
      fTertVertexX     = TertVertex[0];
      fTertVertexY     = TertVertex[1];
      fTertVertexZ     = TertVertex[2];
      // __ chi2 of tert vtx __ //
      fTertVertChi2    = tertVertex->GetChi2();
      fTertVertNDF     = tertVertex->GetNDF();
      // _________________________________ //
      // __ propagate tracks and get impact parameters __ //
      fPropDCADaughter = (exTrack1->PropagateToDCA(tertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
      fImParDaughter   = impar[0];
      fImParzDaughter  = impar[1];
      fDcaSecDaughter  = TMath::Abs(exTrack1->GetD(TertVertex[0], TertVertex[1], kMagF));
      fPropDCADaughter1= (exTrack2->PropagateToDCA(tertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
      fImParDaughter1  = impar[0];
      fImParzDaughter1 = impar[1];
      fDcaSecDaughter1 = TMath::Abs(exTrack2->GetD(TertVertex[0], TertVertex[1], kMagF));
      fPropDCADaughter2= (exTrack3->PropagateToDCA(tertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
      fImParDaughter2  = impar[0];
      fImParzDaughter2 = impar[1];
      fDcaSecDaughter2 = TMath::Abs(exTrack3->GetD(TertVertex[0], TertVertex[1], kMagF));
      if (tertVertex) delete tertVertex;
      // _________________________________________________ //
      particle1->SetXYZM(2. * exTrack1->Px(), 2. * exTrack1->Py(), 2. * exTrack1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
      particle2->SetXYZM(exTrack2->Px(), exTrack2->Py(), exTrack2->Pz(), AliPID::ParticleMass(AliPID::kProton));
      particle3->SetXYZM(exTrack3->Px(), exTrack3->Py(), exTrack3->Pz(), AliPID::ParticleMass(AliPID::kPion));
      sublorentzsum->SetXYZM(0., 0., 0., 0.);
      *sublorentzsum = *particle1 + *particle2 + *particle3;
      // _________________________________ //
      // __ Setting up daughter hypernucleus track here __ //
      // __ cov matrix calculated from daughter tracks __ //
      exTrack1->CheckCovariance();            exTrack2->CheckCovariance();            exTrack3->CheckCovariance();
      exTrack1->GetCovarianceXYZPxPyPz(cov0); exTrack2->GetCovarianceXYZPxPyPz(cov1); exTrack3->GetCovarianceXYZPxPyPz(cov2);
      for (int i = 0; i < 21; i++) cov[i] = TMath::Sqrt(cov0[i] * cov0[i] + cov1[i] * cov1[i] + cov2[i] * cov2[i]);
      // __ daughter hypernucleus sign = helium track sign __ //
      sign = ksign;
      // __ use coord from tert vtx for daughter hypernucleus track __ //
      xyz[0] = TertVertex[0]; xyz[1] = TertVertex[1]; xyz[2] = TertVertex[2];
      // __ daughter hypernucleus momentum __ //
      pxpypz[0] = sublorentzsum->Px(); pxpypz[1] = sublorentzsum->Py(); pxpypz[2] = sublorentzsum->Pz();
      // __ Create daughter hypernucleus track __ //
      exTrack->Set(xyz, pxpypz, cov, sign);
      // __ check cov __ //
      exTrack->CheckCovariance(); exTrack4->CheckCovariance();
      // _________________________________ //
      // __ create sec vtx __ //
      trkArray1 = new TObjArray(2);
      trkArray1->AddAt(exTrack, 0);
      trkArray1->AddAt(exTrack4, 1);
      // __ start at prim vtx __ //
      secvertexer->SetVtxStart(primVertex);
      secVertex = (AliESDVertex*)secvertexer->VertexForSelectedTracks(trkArray1, idTwo, kTRUE, kTRUE, kFALSE);
      if (secvertexer) delete secvertexer;
      if (trkArray1) delete trkArray1;
      if (primVertex) delete primVertex;
      // __ coord of sec vtx __ //
      SecVertex[0]      = secVertex->GetX();
      SecVertex[1]      = secVertex->GetY();
      SecVertex[2]      = secVertex->GetZ();
      fSecVertexX       = SecVertex[0];
      fSecVertexY       = SecVertex[1];
      fSecVertexZ       = SecVertex[2];
      // __ chi2 of sec vtx __ //
      fSecVertChi2      = secVertex->GetChi2();
      fSecVertNDF       = secVertex->GetNDF();
      // _________________________________ //
      // __ propagate tracks and get impact parameters __ //
      fPropDCADaughter3 = (exTrack4->PropagateToDCA(secVertex, kMagF, 10, impar) == 1 ? 1 : 0);
      fImParDaughter3   = impar[0];
      fImParzDaughter3  = impar[1];
      fDcaSecDaughter3  = TMath::Abs(exTrack4->GetD(SecVertex[0], SecVertex[1], kMagF));
      fPropDCADaughter4 = (exTrack->PropagateToDCA(secVertex, kMagF, 10, impar) == 1 ? 1 : 0);
      fImParDaughter4   = impar[0];
      fImParzDaughter4  = impar[1];
      fDcaSecDaughter4  = TMath::Abs(exTrack->GetD(SecVertex[0], SecVertex[1], kMagF));
      if (secVertex) delete secVertex;
      // _________________________________ //
      // __ dca between tracks __ //
      fDCA2B  = TMath::Abs(exTrack4->GetDCA(exTrack, kMagF, xthiss, xpp));
      fDCA3B1 = TMath::Abs(exTrack2->GetDCA(exTrack1, kMagF, xthiss, xpp));
      fDCA3B2 = TMath::Abs(exTrack3->GetDCA(exTrack1, kMagF, xthiss, xpp));
      fDCA3B3 = TMath::Abs(exTrack3->GetDCA(exTrack2, kMagF, xthiss, xpp));
    }
    // _________________________________ //
    else {
      // _________________________________________________ //
      // __ tert vtx __ //
      trkArray = new TObjArray(3);
      trkArray->AddAt(exTrack1, 0);
      trkArray->AddAt(exTrack2, 1);
      trkArray->AddAt(exTrack3, 2);
      // __start at prim vtx __ //
      tertvertexer->SetVtxStart(primVertex);
      tertVertex = (AliESDVertex*)tertvertexer->VertexForSelectedTracks(trkArray, idThree, kTRUE, kTRUE, kFALSE);
      if (trkArray)   delete trkArray;
      // __ coord of tert vtx __ //
      TertVertex[0]    = tertVertex->GetX();
      TertVertex[1]    = tertVertex->GetY();
      TertVertex[2]    = tertVertex->GetZ();
      fTertVertexX     = TertVertex[0];
      fTertVertexY     = TertVertex[1];
      fTertVertexZ     = TertVertex[2];
      // __ chi2 of tert vtx __ //
      fTertVertChi2    = tertVertex->GetChi2();
      fTertVertNDF     = tertVertex->GetNDF();
      // _________________________________ //
      // __ propagate tracks and get impact parameters __ //
      fPropDCADaughter = (exTrack1->PropagateToDCA(tertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
      fImParDaughter   = impar[0];
      fImParzDaughter  = impar[1];
      fDcaSecDaughter  = TMath::Abs(exTrack1->GetD(TertVertex[0], TertVertex[1], kMagF));
      fPropDCADaughter1= (exTrack2->PropagateToDCA(tertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
      fImParDaughter1  = impar[0];
      fImParzDaughter1 = impar[1];
      fDcaSecDaughter1 = TMath::Abs(exTrack2->GetD(TertVertex[0], TertVertex[1], kMagF));
      fPropDCADaughter2= (exTrack3->PropagateToDCA(tertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
      fImParDaughter2  = impar[0];
      fImParzDaughter2 = impar[1];
      fDcaSecDaughter2 = TMath::Abs(exTrack3->GetD(TertVertex[0], TertVertex[1], kMagF));
      if (tertVertex) delete tertVertex;
      // _________________________________ //
      particle1->SetXYZM(2. * exTrack1->Px(), 2. * exTrack1->Py(), 2. * exTrack1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
      particle2->SetXYZM(exTrack2->Px(), exTrack2->Py(), exTrack2->Pz(), AliPID::ParticleMass(AliPID::kProton));
      particle3->SetXYZM(exTrack3->Px(), exTrack3->Py(), exTrack3->Pz(), AliPID::ParticleMass(AliPID::kPion));
      sublorentzsum->SetXYZM(0., 0., 0., 0.);
      *sublorentzsum = *particle1 + *particle2 + *particle3;
      // _________________________________ //
      // __ Setting up daughter hypernucleus track here __ //
      // __ cov matrix calculated from daughter tracks __ //
      exTrack1->CheckCovariance();            exTrack2->CheckCovariance();            exTrack3->CheckCovariance();
      exTrack1->GetCovarianceXYZPxPyPz(cov0); exTrack2->GetCovarianceXYZPxPyPz(cov1); exTrack3->GetCovarianceXYZPxPyPz(cov2);
      for (int i = 0; i < 21; i++) cov[i] = TMath::Sqrt(cov0[i] * cov0[i] + cov1[i] * cov1[i] + cov2[i] * cov2[i]);
      // __ daughter hypernucleus sign = helium track sign __ //
      sign = ksign;
      // __ use coord from prim vtx for daughter hypernucleus track __ //
      xyz[0] = PrimVertex[0]; xyz[1] = PrimVertex[1]; xyz[2] = PrimVertex[2];
      // __ daughter hypernucleus momentum __ //
      pxpypz[0] = sublorentzsum->Px(); pxpypz[1] = sublorentzsum->Py(); pxpypz[2] = sublorentzsum->Pz();
      // __ Create daughter hypernucleus track __ //
      exTrack->Set(xyz, pxpypz, cov, sign);
      // __ check cov __ //
      exTrack->CheckCovariance(); exTrack4->CheckCovariance();
      // _________________________________ //
      // __ create sec vtx __ //
      trkArray1 = new TObjArray(2);
      trkArray1->AddAt(exTrack, 0);
      trkArray1->AddAt(exTrack4, 1);
      // __ start at prim vtx __ //
      secvertexer->SetVtxStart(primVertex);
      secVertex = (AliESDVertex*)secvertexer->VertexForSelectedTracks(trkArray1, idTwo, kTRUE, kTRUE, kFALSE);
      if (secvertexer) delete secvertexer;
      if (trkArray1) delete trkArray1;
      if (primVertex) delete primVertex;
      // __ coord of sec vtx __ //
      SecVertex[0]      = secVertex->GetX();
      SecVertex[1]      = secVertex->GetY();
      SecVertex[2]      = secVertex->GetZ();
      fSecVertexX       = SecVertex[0];
      fSecVertexY       = SecVertex[1];
      fSecVertexZ       = SecVertex[2];
      // __ chi2 of sec vtx __ //
      fSecVertChi2      = secVertex->GetChi2();
      fSecVertNDF       = secVertex->GetNDF();
      // _________________________________ //
      // __ propagate tracks and get impact parameters __ //
      fPropDCADaughter3 = (exTrack4->PropagateToDCA(secVertex, kMagF, 10, impar) == 1 ? 1 : 0);
      fImParDaughter3   = impar[0];
      fImParzDaughter3  = impar[1];
      fDcaSecDaughter3  = TMath::Abs(exTrack4->GetD(SecVertex[0], SecVertex[1], kMagF));
      fPropDCADaughter4 = (exTrack->PropagateToDCA(secVertex, kMagF, 10, impar) == 1 ? 1 : 0);
      fImParDaughter4   = impar[0];
      fImParzDaughter4  = impar[1];
      fDcaSecDaughter4  = TMath::Abs(exTrack->GetD(SecVertex[0], SecVertex[1], kMagF));
      if (secVertex) delete secVertex;
      // _________________________________ //
      // __ dca between tracks __ //
      fDCA2B  = TMath::Abs(exTrack4->GetDCA(exTrack, kMagF, xthiss, xpp));
      fDCA3B1 = TMath::Abs(exTrack2->GetDCA(exTrack1, kMagF, xthiss, xpp));
      fDCA3B2 = TMath::Abs(exTrack3->GetDCA(exTrack1, kMagF, xthiss, xpp));
      fDCA3B3 = TMath::Abs(exTrack3->GetDCA(exTrack2, kMagF, xthiss, xpp));
      // _________________________________ //            
    }
    // _________________________________________________ //
    particle1->SetXYZM(2. * exTrack1->Px(), 2. * exTrack1->Py(), 2. * exTrack1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
    particle2->SetXYZM(exTrack2->Px(), exTrack2->Py(), exTrack2->Pz(), AliPID::ParticleMass(AliPID::kProton));
    particle3->SetXYZM(exTrack3->Px(), exTrack3->Py(), exTrack3->Pz(), AliPID::ParticleMass(AliPID::kPion));
    particle4->SetXYZM(exTrack4->Px(), exTrack4->Py(), exTrack4->Pz(), AliPID::ParticleMass(AliPID::kPion));
    sublorentzsum->SetXYZM(0., 0., 0., 0.);
    lorentzsum->SetXYZM(0., 0., 0., 0.);
    *lorentzsum = *particle1 + *particle2 + *particle3 + *particle4;
    *sublorentzsum = *particle1 + *particle2 + *particle3;
    // _________________________________ //
    // __ mother hypernucleus ctau and cos(PA) __ //
    dd[0] = PrimVertex[0] - SecVertex[0];
    dd[1] = PrimVertex[1] - SecVertex[1];
    dd[2] = PrimVertex[2] - SecVertex[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    fctMother = (lorentzsum->M() * TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2))) / lorentzsum->P();
    fPA = TMath::Cos(lorentzsum->Angle(*h));
    // _________________________________________________ //
    // __ daughter hypernucleus ctau and cos(PA) to sec vtx __ //
    dd[0] = SecVertex[0] - TertVertex[0];
    dd[1] = SecVertex[1] - TertVertex[1];
    dd[2] = SecVertex[2] - TertVertex[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    fctSubMother = (sublorentzsum->M() * TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2))) / sublorentzsum->P();
    fSubPA = TMath::Cos(sublorentzsum->Angle(*h));
    // _________________________________________________ //
    // __ daughter hypernucleus cos(PA) to prim vtx __ //
    dd[0] = PrimVertex[0] - TertVertex[0];
    dd[1] = PrimVertex[1] - TertVertex[1];
    dd[2] = PrimVertex[2] - TertVertex[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    fSubPA2 = TMath::Cos(sublorentzsum->Angle(*h));
    fDecAngle = sublorentzsum->Angle(particle4->Vect());
    // _________________________________________________ //
    if (TMath::Abs(fDCA2B) > kDCATracksCut || TMath::Abs(fDCA3B1) > kDCATracksCut
	|| TMath::Abs(fDCA3B2) > kDCATracksCut || TMath::Abs(fDCA3B3) > kDCATracksCut) return 0;
    if (fSubPA2 < kPointingAngleCut) return 0;
    if (fSubPA < kPointingAngleCut) return 0;
    if (fPA < kPointingAngleCut) return 0;
    // _________________________________________________ //
    // __ cuts from 04/2020 __ //
    fDCA2Bo = TMath::Abs(exTrack4->GetD(sublorentzsum->X(), sublorentzsum->Y(), kMagF));
    fDcaDaughtero = exTrack1->GetD(fPrimVertexX, fPrimVertexY, kMagF);
    fDcaDaughter1o = exTrack2->GetD(fPrimVertexX, fPrimVertexY, kMagF);
    fDcaDaughter2o = exTrack3->GetD(fPrimVertexX, fPrimVertexY, kMagF);
    fDcaDaughter3o = exTrack4->GetD(fPrimVertexX, fPrimVertexY, kMagF);
    // _________________________________________________ //
    // __ track information __ //
    AliAnalysisTaskDoubleHypNucTreeLS::SetDaughterInformation("4LLH", kDecayChannel);
    // _________________________________________________ //
    // __ information of hypernuclei track __ //
    Float_t xv[2], yv[3];
    exTrack->GetImpactParameters(xv, yv);
    fpDaughter4 = exTrack->P();
    fptDaughter4 = sublorentzsum->Pt();
    fpxDaughter4 = sublorentzsum->Px();
    fpyDaughter4 = sublorentzsum->Py();
    fpzDaughter4 = sublorentzsum->Pz();
    fEDaughter4 = sublorentzsum->E();
    fyDaughter4 = sublorentzsum->Rapidity();
    fDcaDaughter4 = TMath::Abs(exTrack->GetD(PrimVertex[0], PrimVertex[1], kMagF));//xv[0];
    fDcazDaughter4 = xv[1];
    fSigmaYXDaughter4 = yv[0];
    fSigmaXYZDaughter4 = yv[1];
    fSigmaZDaughter4 = yv[2];
    fPtUncertDaughter4 = TMath::Sqrt(exTrack->GetSigma1Pt2()) * fptDaughter4;
    fDcaSecDaughter4 = TMath::Abs(exTrack->GetD(SecVertex[0], SecVertex[1], kMagF));
    exTrack->GetCovarianceXYZPxPyPz(fCovMatrixTrack4);
    fTrackPar4[0] = exTrack->GetX(); fTrackPar4[1] = exTrack->GetY(); fTrackPar4[2] = exTrack->GetZ(); fTrackPar4[3] = exTrack->GetSnp();
    fTrackPar4[4] = exTrack->GetTgl(); fTrackPar4[5] = exTrack->GetSigned1Pt(); fTrackPar4[6] = exTrack->GetAlpha();
    // _________________________________________________ //
    // __ mother hypernucleus information __ //
    fPDGMother    = ksign * fgkPdgCode[kPDGDoubleHyperHydrogen4];
    fChargeMother = ksign;
    fDecayChannel = kDecayChannel;
    fRecoMethod   = 1;
    fmMother      = lorentzsum->M();
    fmMother2     = lorentzsum2->M();
    fEMother      = lorentzsum->E();
    fpxMother     = lorentzsum->Px();
    fpyMother     = lorentzsum->Py();
    fpzMother     = lorentzsum->Pz();
    fptMother     = lorentzsum->Pt();
    fpMother      = lorentzsum->P();
    fyMother      = lorentzsum->Rapidity();
    // _________________________________________________ //
    // __ daughter hypernucleus information __ //
    fmSubMother   = sublorentzsum->M();
    fESubMother   = sublorentzsum->E();
    fpxSubMother  = sublorentzsum->Px();
    fpySubMother  = sublorentzsum->Py();
    fpzSubMother  = sublorentzsum->Pz();
    fptSubMother  = sublorentzsum->Pt();
    fpSubMother   = sublorentzsum->P();
    fySubMother   = sublorentzsum->Rapidity();
    // __ Check for V0s __ //
    AliAnalysisTaskDoubleHypNucTreeLS::GetV0Status();
    // _________________________________________________ //	  	  
    // daughter momenta
    fpxDaughter   = particle1->Px();
    fpyDaughter   = particle1->Py();
    fpzDaughter   = particle1->Pz();
    fptDaughter   = particle1->Pt();
    fEDaughter    = particle1->E();
    fyDaughter    = particle1->Rapidity();
    fpxDaughter1  = particle2->Px();
    fpyDaughter1  = particle2->Py();
    fpzDaughter1  = particle2->Pz();
    fptDaughter1  = particle2->Pt();
    fEDaughter1   = particle2->E();
    fyDaughter1   = particle2->Rapidity();
    fpxDaughter2  = particle3->Px();
    fpyDaughter2  = particle3->Py();
    fpzDaughter2  = particle3->Pz();
    fptDaughter2  = particle3->Pt();
    fEDaughter2   = particle3->E();
    fyDaughter2   = particle3->Rapidity();
    fpxDaughter3  = particle4->Px();
    fpyDaughter3  = particle4->Py();
    fpzDaughter3  = particle4->Pz();
    fptDaughter3  = particle4->Pt();
    fEDaughter3   = particle4->E();
    fyDaughter3   = particle4->Rapidity();
    // _________________________________________________ //	 
    // __ armenteros podolanski __ //
    TVector3 vecN(0., 0., 0.);
    TVector3 vecP(0., 0., 0.);
    TVector3 vecM(0., 0., 0.);

    if (ksign == 1) {
      vecP = sublorentzsum->Vect();
      vecN = particle4->Vect();
      vecM = lorentzsum->Vect();
    }
    if (ksign == -1) {
      vecN = sublorentzsum->Vect();
      vecP = particle4->Vect();
      vecM = lorentzsum->Vect();
    }

    fthetaP   = TMath::ACos((vecP * vecM) / (vecP.Mag() * vecM.Mag()));
    fthetaN   = TMath::ACos((vecN * vecM) / (vecN.Mag() * vecM.Mag()));
    farmalpha = ((vecP.Mag()) * TMath::Cos(fthetaP) - (vecN.Mag()) * TMath::Cos(fthetaN)) / ((vecP.Mag()) * TMath::Cos(fthetaP) + (vecN.Mag()) * TMath::Cos(fthetaN));
    farmpt    = vecP.Mag() * sin(fthetaP);
    return 1;
  }
  // _________________________________ //
  else if (Mother == "4LLH" && kDecayChannel == 2) {
    // ******************************************************************************
    // * Two different methods:                                                     *
    // * kVariante = 2 creates the tertiary vertex first                            *
    // * which coordinates are taken to reconstruct the daughter hypernucleus track *
    // * MC checks showed, that the Pointing resolution is worse in this case.      *
    // * kVariante = 1 uses the coordinates of the primary vertex to reconstruct    *
    // * the daughter hypernucleus track (default).                                 *
    // ******************************************************************************
    // __ init AliVertexerTracks __ //
    secvertexer = new AliVertexerTracks(fAODevent->GetMagneticField());
    tertvertexer = new AliVertexerTracks(fAODevent->GetMagneticField());
    secvertexer->SetFinderAlgorithm(kFinderAlgo);
    tertvertexer->SetFinderAlgorithm(kFinderAlgo);
    // _________________________________ //
    exTrack->Reset(); exTrack1->Reset(); exTrack2->Reset(); exTrack3->Reset(); exTrack4->Reset();
    exTrack1->CopyFromVTrack(track1);
    exTrack2->CopyFromVTrack(track2);
    exTrack3->CopyFromVTrack(track3);
    exTrack4->CopyFromVTrack(track4);
    // _________________________________ //
    if (kVariante == 2) {
      // __ tert vtx __ //
      trkArray = new TObjArray(2);
      trkArray->AddAt(exTrack1, 0);
      trkArray->AddAt(exTrack3, 1);
      // __start at prim vtx __ //
      tertvertexer->SetVtxStart(primVertex);
      tertVertex = (AliESDVertex*)tertvertexer->VertexForSelectedTracks(trkArray, idTwo, kTRUE, kTRUE, kFALSE);
      if (trkArray) delete trkArray;
      if (tertvertexer) delete tertvertexer;
      // __ coord of tert vtx __ //
      TertVertex[0]    = tertVertex->GetX();
      TertVertex[1]    = tertVertex->GetY();
      TertVertex[2]    = tertVertex->GetZ();
      fTertVertexX     = TertVertex[0];
      fTertVertexY     = TertVertex[1];
      fTertVertexZ     = TertVertex[2];
      // __ chi2 of tert vtx __ //
      fTertVertChi2    = tertVertex->GetChi2();
      fTertVertNDF     = tertVertex->GetNDF();
      // _________________________________ //
      // __ propagate tracks and get impact parameters __ //
      fPropDCADaughter = (exTrack1->PropagateToDCA(tertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
      fImParDaughter   = impar[0];
      fImParzDaughter  = impar[1];
      fDcaSecDaughter  = TMath::Abs(exTrack1->GetD(TertVertex[0], TertVertex[1], kMagF));
      fPropDCADaughter2= (exTrack3->PropagateToDCA(tertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
      fImParDaughter2  = impar[0];
      fImParzDaughter2 = impar[1];
      fDcaSecDaughter2 = TMath::Abs(exTrack3->GetD(TertVertex[0], TertVertex[1], kMagF));
      if (tertVertex) delete tertVertex;
      // _________________________________________________ //
      particle1->SetXYZM(2. * exTrack1->Px(), 2. * exTrack1->Py(), 2. * exTrack1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
      particle3->SetXYZM(exTrack3->Px(), exTrack3->Py(), exTrack3->Pz(), AliPID::ParticleMass(AliPID::kPion));
      sublorentzsum->SetXYZM(0., 0., 0., 0.);
      *sublorentzsum = *particle1 + *particle3;
      // _________________________________ //
      // __ Setting up daughter hypernucleus track here __ //
      // __ cov matrix calculated from daughter tracks __ //
      exTrack1->CheckCovariance();            exTrack3->CheckCovariance();
      exTrack1->GetCovarianceXYZPxPyPz(cov0); exTrack3->GetCovarianceXYZPxPyPz(cov1);
      for (int i = 0; i < 21; i++) cov[i] = TMath::Sqrt(cov0[i] * cov0[i] + cov1[i] * cov1[i]);
      // __ daughter hypernucleus sign = helium track sign __ //
      sign = ksign;
      // __ use coord from tert vtx for daughter hypernucleus track __ //
      xyz[0] = TertVertex[0]; xyz[1] = TertVertex[1]; xyz[2] = TertVertex[2];
      // __ daughter hypernucleus momentum __ //
      pxpypz[0] = sublorentzsum->Px(); pxpypz[1] = sublorentzsum->Py(); pxpypz[2] = sublorentzsum->Pz();
      // __ Create daughter hypernucleus track __ //
      exTrack->Set(xyz, pxpypz, cov, sign);
      // __ check cov __ //
      exTrack->CheckCovariance(); exTrack2->CheckCovariance(); exTrack4->CheckCovariance();
      // _________________________________ //
      // __ create sec vtx __ //
      trkArray1 = new TObjArray(3);
      trkArray1->AddAt(exTrack, 0);
      trkArray1->AddAt(exTrack2, 1);
      trkArray1->AddAt(exTrack4, 2);
      // __ start at prim vtx __ //
      secvertexer->SetVtxStart(primVertex);
      secVertex = (AliESDVertex*)secvertexer->VertexForSelectedTracks(trkArray1, idThree, kTRUE, kTRUE, kFALSE);
      if (trkArray1) delete trkArray1;
      if (secvertexer) delete secvertexer;
      if (primVertex) delete primVertex;
      // __ coord of sec vtx __ //
      SecVertex[0]      = secVertex->GetX();
      SecVertex[1]      = secVertex->GetY();
      SecVertex[2]      = secVertex->GetZ();
      fSecVertexX       = SecVertex[0];
      fSecVertexY       = SecVertex[1];
      fSecVertexZ       = SecVertex[2];
      // __ chi2 of sec vtx __ //
      fSecVertChi2      = secVertex->GetChi2();
      fSecVertNDF       = secVertex->GetNDF();
      // _________________________________ //
      // __ propagate tracks and get impact parameters __ //
      fPropDCADaughter1 = (exTrack2->PropagateToDCA(secVertex, kMagF, 10, impar) == 1 ? 1 : 0);
      fImParDaughter1   = impar[0];
      fImParzDaughter1  = impar[1];
      fDcaSecDaughter1  = TMath::Abs(exTrack2->GetD(SecVertex[0], SecVertex[1], kMagF));
      fPropDCADaughter3 = (exTrack4->PropagateToDCA(secVertex, kMagF, 10, impar) == 1 ? 1 : 0);
      fImParDaughter3   = impar[0];
      fImParzDaughter3  = impar[1];
      fDcaSecDaughter3  = TMath::Abs(exTrack4->GetD(SecVertex[0], SecVertex[1], kMagF));
      fPropDCADaughter4 = (exTrack->PropagateToDCA(secVertex, kMagF, 10, impar) == 1 ? 1 : 0);
      fImParDaughter4   = impar[0];
      fImParzDaughter4  = impar[1];
      fDcaSecDaughter4  = TMath::Abs(exTrack->GetD(SecVertex[0], SecVertex[1], kMagF));
      if (secVertex) delete secVertex;
      // _________________________________ //
      // __ dca between tracks __ //
      fDCA2B  = TMath::Abs(exTrack3->GetDCA(exTrack1, kMagF, xthiss, xpp));
      fDCA3B1 = TMath::Abs(exTrack2->GetDCA(exTrack, kMagF, xthiss, xpp));
      fDCA3B2 = TMath::Abs(exTrack4->GetDCA(exTrack, kMagF, xthiss, xpp));
      fDCA3B3 = TMath::Abs(exTrack4->GetDCA(exTrack2, kMagF, xthiss, xpp));
    }
    // _________________________________ //
    else {
      // _________________________________________________ //
      // __ tert vtx __ //
      trkArray = new TObjArray(2);
      trkArray->AddAt(exTrack1, 0);
      trkArray->AddAt(exTrack3, 1);
      // __start at prim vtx __ //
      tertvertexer->SetVtxStart(primVertex);
      tertVertex = (AliESDVertex*)tertvertexer->VertexForSelectedTracks(trkArray, idTwo, kTRUE, kTRUE, kFALSE);
      if (trkArray) delete trkArray;
      // __ coord of tert vtx __ //
      TertVertex[0]    = tertVertex->GetX();
      TertVertex[1]    = tertVertex->GetY();
      TertVertex[2]    = tertVertex->GetZ();
      fTertVertexX     = TertVertex[0];
      fTertVertexY     = TertVertex[1];
      fTertVertexZ     = TertVertex[2];
      // __ chi2 of tert vtx __ //
      fTertVertChi2    = tertVertex->GetChi2();
      fTertVertNDF     = tertVertex->GetNDF();
      // _________________________________ //
      // __ propagate tracks and get impact parameters __ //
      fPropDCADaughter = (exTrack1->PropagateToDCA(tertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
      fImParDaughter   = impar[0];
      fImParzDaughter  = impar[1];
      fDcaSecDaughter  = TMath::Abs(exTrack1->GetD(TertVertex[0], TertVertex[1], kMagF));
      fPropDCADaughter2= (exTrack3->PropagateToDCA(tertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
      fImParDaughter2  = impar[0];
      fImParzDaughter2 = impar[1];
      fDcaSecDaughter2 = TMath::Abs(exTrack3->GetD(TertVertex[0], TertVertex[1], kMagF));
      if (tertVertex) delete tertVertex;
      // _________________________________ //
      particle1->SetXYZM(2. * exTrack1->Px(), 2. * exTrack1->Py(), 2. * exTrack1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
      particle3->SetXYZM(exTrack3->Px(), exTrack3->Py(), exTrack3->Pz(), AliPID::ParticleMass(AliPID::kPion));
      sublorentzsum->SetXYZM(0., 0., 0., 0.);
      *sublorentzsum = *particle1 + *particle3;
      // __ Setting up daughter hypernucleus track here __ //
      // __ cov matrix calculated from daughter tracks __ //
      exTrack1->CheckCovariance();            exTrack3->CheckCovariance();
      exTrack1->GetCovarianceXYZPxPyPz(cov0); exTrack3->GetCovarianceXYZPxPyPz(cov1);
      for (int i = 0; i < 21; i++) cov[i] = TMath::Sqrt(cov0[i] * cov0[i] + cov1[i] * cov1[i]);
      // __ daughter hypernucleus sign = helium track sign __ //
      sign = ksign;
      // __ use coord from prim vtx for daughter hypernucleus track __ //
      xyz[0] = PrimVertex[0]; xyz[1] = PrimVertex[1]; xyz[2] = PrimVertex[2];
      // __ daughter hypernucleus momentum __ //
      pxpypz[0] = sublorentzsum->Px(); pxpypz[1] = sublorentzsum->Py(); pxpypz[2] = sublorentzsum->Pz();
      // __ Create daughter hypernucleus track __ //
      exTrack->Set(xyz, pxpypz, cov, sign);
      // __ check cov __ //
      exTrack->CheckCovariance(); exTrack2->CheckCovariance(); exTrack4->CheckCovariance();
      // _________________________________ //
      // __ create sec vtx __ //
      trkArray1 = new TObjArray(3);
      trkArray1->AddAt(exTrack, 0);
      trkArray1->AddAt(exTrack2, 1);
      trkArray1->AddAt(exTrack4, 2);
      // __ start at prim vtx __ //
      secvertexer->SetVtxStart(primVertex);
      secVertex = (AliESDVertex*)secvertexer->VertexForSelectedTracks(trkArray1, idThree, kTRUE, kTRUE, kFALSE);
      if (trkArray1) delete trkArray1;
      if (secvertexer) delete secvertexer;
      if (primVertex)   delete primVertex;
      // __ coord of sec vtx __ //
      SecVertex[0]      = secVertex->GetX();
      SecVertex[1]      = secVertex->GetY();
      SecVertex[2]      = secVertex->GetZ();
      fSecVertexX       = SecVertex[0];
      fSecVertexY       = SecVertex[1];
      fSecVertexZ       = SecVertex[2];
      // __ chi2 of sec vtx __ //
      fSecVertChi2      = secVertex->GetChi2();
      fSecVertNDF       = secVertex->GetNDF();
      // _________________________________ //
      // __ propagate tracks and get impact parameters __ //
      fPropDCADaughter1 = (exTrack2->PropagateToDCA(secVertex, kMagF, 10, impar) == 1 ? 1 : 0);
      fImParDaughter1   = impar[0];
      fImParzDaughter1  = impar[1];
      fDcaSecDaughter1  = TMath::Abs(exTrack2->GetD(SecVertex[0], SecVertex[1], kMagF));
      fPropDCADaughter3 = (exTrack4->PropagateToDCA(secVertex, kMagF, 10, impar) == 1 ? 1 : 0);
      fImParDaughter3   = impar[0];
      fImParzDaughter3  = impar[1];
      fDcaSecDaughter3  = TMath::Abs(exTrack4->GetD(SecVertex[0], SecVertex[1], kMagF));
      fPropDCADaughter4 = (exTrack->PropagateToDCA(secVertex, kMagF, 10, impar) == 1 ? 1 : 0);
      fImParDaughter4   = impar[0];
      fImParzDaughter4  = impar[1];
      fDcaSecDaughter4  = TMath::Abs(exTrack->GetD(SecVertex[0], SecVertex[1], kMagF));
      if (secVertex) delete secVertex;
      // _________________________________ //      
      // __ dca between tracks __ //
      fDCA2B  = TMath::Abs(exTrack3->GetDCA(exTrack1, kMagF, xthiss, xpp));
      fDCA3B1 = TMath::Abs(exTrack2->GetDCA(exTrack, kMagF, xthiss, xpp));
      fDCA3B2 = TMath::Abs(exTrack4->GetDCA(exTrack, kMagF, xthiss, xpp));
      fDCA3B3 = TMath::Abs(exTrack4->GetDCA(exTrack2, kMagF, xthiss, xpp));
      // _________________________________ //    
    }
    // _________________________________________________ //
    particle1->SetXYZM(2. * exTrack1->Px(), 2. * exTrack1->Py(), 2. * exTrack1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
    particle2->SetXYZM(exTrack2->Px(), exTrack2->Py(), exTrack2->Pz(), AliPID::ParticleMass(AliPID::kProton));
    particle3->SetXYZM(exTrack3->Px(), exTrack3->Py(), exTrack3->Pz(), AliPID::ParticleMass(AliPID::kPion));
    particle4->SetXYZM(exTrack4->Px(), exTrack4->Py(), exTrack4->Pz(), AliPID::ParticleMass(AliPID::kPion));
    sublorentzsum->SetXYZM(0., 0., 0., 0.);
    lorentzsum->SetXYZM(0., 0., 0., 0.);
    *lorentzsum = *particle1 + *particle2 + *particle3 + *particle4;
    *sublorentzsum = *particle1 + *particle3;
    // _________________________________ //
    // __ mother hypernucleus ctau and cos(PA) __ //
    dd[0] = PrimVertex[0] - SecVertex[0];
    dd[1] = PrimVertex[1] - SecVertex[1];
    dd[2] = PrimVertex[2] - SecVertex[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    fctMother = (lorentzsum->M() * TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2))) / lorentzsum->P();
    fPA = TMath::Cos(lorentzsum->Angle(*h));
    // _________________________________________________ //
    // __ daughter hypernucleus ctau and cos(PA) to sec vtx __ //
    dd[0] = SecVertex[0] - TertVertex[0];
    dd[1] = SecVertex[1] - TertVertex[1];
    dd[2] = SecVertex[2] - TertVertex[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    fctSubMother = (sublorentzsum->M() * TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2))) / sublorentzsum->P();
    fSubPA = TMath::Cos(sublorentzsum->Angle(*h));
    // _________________________________________________ //
    // __ daughter hypernucleus cos(PA) to prim vtx __ //
    dd[0] = PrimVertex[0] - TertVertex[0];
    dd[1] = PrimVertex[1] - TertVertex[1];
    dd[2] = PrimVertex[2] - TertVertex[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    fSubPA2 = TMath::Cos(sublorentzsum->Angle(*h));
    fDecAngle = sublorentzsum->Angle(particle4->Vect());
    // _________________________________________________ //  
    if (TMath::Abs(fDCA2B) > kDCATracksCut || TMath::Abs(fDCA3B1) > kDCATracksCut
	|| TMath::Abs(fDCA3B2) > kDCATracksCut || TMath::Abs(fDCA3B3) > kDCATracksCut) return 0;
    if (fSubPA2 < kPointingAngleCut) return 0;
    if (fSubPA < kPointingAngleCut) return 0;
    if (fPA < kPointingAngleCut) return 0;
    // _________________________________________________ //
    // __ track information __ //
    AliAnalysisTaskDoubleHypNucTreeLS::SetDaughterInformation("4LLH", kDecayChannel);
    // _________________________________________________ //
    // __ information of hypernuclei track __ //
    Float_t xv[2], yv[3];
    exTrack->GetImpactParameters(xv, yv);
    fpDaughter4        = exTrack->P();
    fptDaughter4       = sublorentzsum->Pt();
    fpxDaughter4       = sublorentzsum->Px();
    fpyDaughter4       = sublorentzsum->Py();
    fpzDaughter4       = sublorentzsum->Pz();
    fEDaughter4        = sublorentzsum->E();
    fyDaughter4        = sublorentzsum->Rapidity();
    fDcaDaughter4      = TMath::Abs(exTrack->GetD(PrimVertex[0], PrimVertex[1], kMagF));//xv[0];
    fDcazDaughter4     = xv[1];
    fSigmaYXDaughter4  = yv[0];
    fSigmaXYZDaughter4 = yv[1];
    fSigmaZDaughter4   = yv[2];
    fPtUncertDaughter4 = TMath::Sqrt(exTrack->GetSigma1Pt2()) * fptDaughter4;
    fDcaSecDaughter4   = TMath::Abs(exTrack->GetD(SecVertex[0], SecVertex[1], kMagF));
    exTrack->GetCovarianceXYZPxPyPz(fCovMatrixTrack4);
    fTrackPar4[0]      = exTrack->GetX(); fTrackPar4[1] = exTrack->GetY(); fTrackPar4[2] = exTrack->GetZ(); fTrackPar4[3] = exTrack->GetSnp();
    fTrackPar4[4]      = exTrack->GetTgl(); fTrackPar4[5] = exTrack->GetSigned1Pt(); fTrackPar4[6] = exTrack->GetAlpha();
    // _________________________________________________ //
    // __ mother hypernucleus information __ //
    fPDGMother    = ksign * fgkPdgCode[kPDGDoubleHyperHydrogen4];
    fChargeMother = ksign;
    fDecayChannel = kDecayChannel;
    fRecoMethod   = 1;
    fmMother      = lorentzsum->M();
    fmMother2     = lorentzsum2->M();
    fEMother      = lorentzsum->E();
    fpxMother     = lorentzsum->Px();
    fpyMother     = lorentzsum->Py();
    fpzMother     = lorentzsum->Pz();
    fptMother     = lorentzsum->Pt();
    fpMother      = lorentzsum->P();
    fyMother      = lorentzsum->Rapidity();
    // _________________________________________________ //
    // __ daughter hypernucleus information __ //
    fmSubMother   = sublorentzsum->M();
    fESubMother   = sublorentzsum->E();
    fpxSubMother  = sublorentzsum->Px();
    fpySubMother  = sublorentzsum->Py();
    fpzSubMother  = sublorentzsum->Pz();
    fptSubMother  = sublorentzsum->Pt();
    fpSubMother   = sublorentzsum->P();
    fySubMother   = sublorentzsum->Rapidity();
    // __ Check for V0s __ //
    AliAnalysisTaskDoubleHypNucTreeLS::GetV0Status();
    // _________________________________________________ //	  	
    // daughter momenta
    fpxDaughter  = particle1->Px();
    fpyDaughter  = particle1->Py();
    fpzDaughter  = particle1->Pz();
    fptDaughter  = particle1->Pt();
    fEDaughter   = particle1->E();
    fyDaughter   = particle1->Rapidity();
    fpxDaughter1 = particle2->Px();
    fpyDaughter1 = particle2->Py();
    fpzDaughter1 = particle2->Pz();
    fptDaughter1 = particle2->Pt();
    fEDaughter1  = particle2->E();
    fyDaughter1  = particle2->Rapidity();
    fpxDaughter2 = particle3->Px();
    fpyDaughter2 = particle3->Py();
    fpzDaughter2 = particle3->Pz();
    fptDaughter2 = particle3->Pt();
    fEDaughter2  = particle3->E();
    fyDaughter2  = particle3->Rapidity();
    fpxDaughter3 = particle4->Px();
    fpyDaughter3 = particle4->Py();
    fpzDaughter3 = particle4->Pz();
    fptDaughter3 = particle4->Pt();
    fEDaughter3  = particle4->E();
    fyDaughter3  = particle4->Rapidity();
    // _________________________________________________ //	   
    // __ armenteros podolanski __ //
    TVector3 vecN(0., 0., 0.);
    TVector3 vecP(0., 0., 0.);
    TVector3 vecM(0., 0., 0.);

    if (ksign == 1) {
      vecP = particle1->Vect();
      vecN = particle4->Vect();
      vecM = sublorentzsum->Vect();
    }
    if (ksign == -1) {
      vecN = particle1->Vect();
      vecP = particle4->Vect();
      vecM = sublorentzsum->Vect();
    }

    fthetaP = TMath::ACos((vecP * vecM) / (vecP.Mag() * vecM.Mag()));
    fthetaN = TMath::ACos((vecN * vecM) / (vecN.Mag() * vecM.Mag()));
    farmalpha = ((vecP.Mag()) * TMath::Cos(fthetaP) - (vecN.Mag()) * TMath::Cos(fthetaN)) / ((vecP.Mag()) * TMath::Cos(fthetaP) + (vecN.Mag()) * TMath::Cos(fthetaN));
    farmpt = vecP.Mag() * sin(fthetaP);
    return 1;
  }
  else return 0;
}
// _________________________________________________ //
Int_t AliAnalysisTaskDoubleHypNucTreeLS::KFReconstruction(TString Mother, Int_t kDecayChannel, Int_t ksign) {
  Double_t impar[2];
  Double_t DecLength;
  Double_t DecLengthErr;
  int checkval;
  float l, le;
  // _________________________________________________ //
  // __ get prim vtx __ //
  primVertex        = AliAnalysisTaskDoubleHypNucTreeLS::AODToESDVertex(*vertex);
  primKFVertex      = CreateKFVertex(*primVertex);
  if (primVertex) delete primVertex;
  PrimVertexKF[0]   = primKFVertex.GetX();
  PrimVertexKF[1]   = primKFVertex.GetY();
  PrimVertexKF[2]   = primKFVertex.GetZ();
  fPrimVertexXKF    = PrimVertexKF[0];
  fPrimVertexYKF    = PrimVertexKF[1];
  fPrimVertexZKF    = PrimVertexKF[2];
  fPrimVertexXErrKF = primKFVertex.GetErrX();
  fPrimVertexYErrKF = primKFVertex.GetErrY();
  fPrimVertexZErrKF = primKFVertex.GetErrZ();
  // __ chi2 of prim vtx __ //
  fPrimVertChi2     = primKFVertex.Chi2();
  fPrimVertNDF      = primKFVertex.GetNDF();
  // _________________________________________________ //
  if (Mother == "3LH") {
    // _________________________________________________ //
    // __ initialize KFParticles __ //
    exTrack->Reset(); exTrack1->Reset(); exTrack2->Reset(); exTrack3->Reset(); exTrack4->Reset();
    if (kDecayChannel == 2 || kDecayChannel == 3) {
      exTrack1->CopyFromVTrack(track1);
      exTrack3->CopyFromVTrack(track3);
      KFtrack1 = AliAnalysisTaskDoubleHypNucTreeLS::CreateKFParticle(*exTrack1, AliPID::ParticleMass(AliPID::kHe3),    2);
      KFtrack3 = AliAnalysisTaskDoubleHypNucTreeLS::CreateKFParticle(*exTrack3, AliPID::ParticleMass(AliPID::kPion),   1);
    }
    if (kDecayChannel == 4) {
      exTrack1->CopyFromVTrack(track1);
      exTrack3->CopyFromVTrack(track4);
      KFtrack1 = AliAnalysisTaskDoubleHypNucTreeLS::CreateKFParticle(*exTrack1, AliPID::ParticleMass(AliPID::kHe3),    2);
      KFtrack3 = AliAnalysisTaskDoubleHypNucTreeLS::CreateKFParticle(*exTrack3, AliPID::ParticleMass(AliPID::kPion),   1);
    }
    if (kDecayChannel == 5 || kDecayChannel == 7) {
      exTrack1->CopyFromVTrack(track2);
      exTrack3->CopyFromVTrack(track3);
      KFtrack1 = AliAnalysisTaskDoubleHypNucTreeLS::CreateKFParticle(*exTrack1, AliPID::ParticleMass(AliPID::kProton), 1);
      KFtrack3 = AliAnalysisTaskDoubleHypNucTreeLS::CreateKFParticle(*exTrack3, AliPID::ParticleMass(AliPID::kPion),   1);
    }
    if (kDecayChannel == 6) {
      exTrack1->CopyFromVTrack(track2);
      exTrack3->CopyFromVTrack(track4);
      KFtrack1 = AliAnalysisTaskDoubleHypNucTreeLS::CreateKFParticle(*exTrack1, AliPID::ParticleMass(AliPID::kProton), 1);
      KFtrack3 = AliAnalysisTaskDoubleHypNucTreeLS::CreateKFParticle(*exTrack3, AliPID::ParticleMass(AliPID::kPion),   1);
    }
    // __ construct mother particle __ //
    checkval = 1;
    KFParticleHypNuc KFChecktrack3;
    KFChecktrack3.AddDaughter(KFtrack1);
    if (!KFChecktrack3.CheckDaughter(KFtrack3)) checkval = 0;
    KFParticleHypNuc KFChecktrack1;
    if (!KFChecktrack1.CheckDaughter(KFtrack1)) checkval = 0;
    KFChecktrack1.AddDaughter(KFtrack3);
    KFParticleHypNuc KFSubMother;
    if (checkval) {
      KFSubMother.SetConstructMethod(2);
      KFSubMother.AddDaughter(KFtrack1);
      KFSubMother.AddDaughter(KFtrack3);
      fRecoMethod = 2;
    }
    else {
      KFSubMother.AddDaughter(KFtrack1);
      KFSubMother.AddDaughter(KFtrack3);
      fRecoMethod = 0;
    }
    // __ transport mother particle to its decay vertex __ //
    KFSubMother.TransportToDecayVertex();
    // _________________________________________________ //
    // __ get decay vertex position __ //
    TertVertexKF[0]    = KFSubMother.GetX();
    TertVertexKF[1]    = KFSubMother.GetY();
    TertVertexKF[2]    = KFSubMother.GetZ();
    TertVertexErrKF[0] = KFSubMother.GetErrX();
    TertVertexErrKF[1] = KFSubMother.GetErrY();
    TertVertexErrKF[2] = KFSubMother.GetErrZ();
    // __ coord of tert vtx __ //
    fTertVertexXKF     = TertVertexKF[0];
    fTertVertexYKF     = TertVertexKF[1];
    fTertVertexZKF     = TertVertexKF[2];
    fTertVertexXErrKF  = TertVertexErrKF[0];
    fTertVertexYErrKF  = TertVertexErrKF[1];
    fTertVertexZErrKF  = TertVertexErrKF[2];
    // __ vertex quality KF __ //
    fTertVertChi2KF    = KFSubMother.Chi2();
    fTertVertNDFKF     = KFSubMother.GetNDF();
    // _________________________________________________ //
    particle1->SetXYZM(0., 0., 0., 0.);
    particle3->SetXYZM(0., 0., 0., 0.);
    sublorentzsum->SetXYZM(0., 0., 0., 0.);
    KFtertVertex       = new AliESDVertex(TertVertexKF, TertVertexErrKF);
    fPropDCADaughter   = (exTrack1->PropagateToDCA(KFtertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter     = impar[0];
    fImParzDaughter    = impar[1];
    fDcaSecDaughter    = TMath::Abs(exTrack1->GetD(TertVertexKF[0], TertVertexKF[1], kMagF));
    fPropDCADaughter2  = (exTrack3->PropagateToDCA(KFtertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter2    = impar[0];
    fImParzDaughter2   = impar[1];
    fDcaSecDaughter2   = TMath::Abs(exTrack3->GetD(TertVertexKF[0], TertVertexKF[1], kMagF));
    if (KFtertVertex) delete KFtertVertex;
    if (kDecayChannel < 5) particle1->SetXYZM(2. * exTrack1->Px(), 2. * exTrack1->Py(), 2. * exTrack1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
    else                   particle1->SetXYZM(exTrack1->Px(), exTrack1->Py(), exTrack1->Pz(), AliPID::ParticleMass(AliPID::kProton));
    particle3->SetXYZM(exTrack3->Px(), exTrack3->Py(), exTrack3->Pz(), AliPID::ParticleMass(AliPID::kPion));
    *sublorentzsum = *particle1 + *particle3;
    fmSubMother    = sublorentzsum->M();
    if (kDecayChannel == 2 && (sublorentzsum->M() > kMax3LHMass || sublorentzsum->M() < kMin3LHMass)) return 0;
    // _________________________________________________ //
    // __ KF daughter hypernucleus cos(PA) to prim vtx __ //
    dd[0] = PrimVertexKF[0] - TertVertexKF[0];
    dd[1] = PrimVertexKF[1] - TertVertexKF[1];
    dd[2] = PrimVertexKF[2] - TertVertexKF[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    fSubPAKF = TMath::Cos(sublorentzsum->Angle(*h));
    if (kDecayChannel == 2 && fSubPAKF < kPointingAngleCut) return 0;
    
    DecLength = TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2));
    fctSubMotherKF = (sublorentzsum->M() * DecLength) / sublorentzsum->P();
    DecLengthErr = TMath::Sqrt((fPrimVertexXErrKF * fPrimVertexXErrKF) / (fPrimVertexXKF * fPrimVertexXKF)
			       + (fPrimVertexYErrKF * fPrimVertexYErrKF) / (fPrimVertexYKF * fPrimVertexYKF)
			       + (fPrimVertexZErrKF * fPrimVertexZErrKF) / (fPrimVertexZKF * fPrimVertexZKF)
			       + (fTertVertexXErrKF * fTertVertexXErrKF) / (fTertVertexXKF * fTertVertexXKF)
			       + (fTertVertexYErrKF * fTertVertexYErrKF) / (fTertVertexYKF * fTertVertexYKF)
			       + (fTertVertexZErrKF * fTertVertexZErrKF) / (fTertVertexZKF * fTertVertexZKF)
			       + (KFSubMother.GetErrP() * KFSubMother.GetErrP()) / (KFSubMother.GetP() * KFSubMother.GetP())
			       + (KFSubMother.GetErrMass() * KFSubMother.GetErrMass()) / (KFSubMother.GetMass() * KFSubMother.GetMass())) * fctSubMotherKF;
    fctSubMotherErrKF = DecLengthErr;
    // _________________________________________________ //
    // __ dca tracks KF __ //	
    fDCA3B1XYKF = (Float_t)KFtrack1.GetDistanceFromParticleXY(KFtrack3);
    fDCA3B1ZKF  = (Float_t)KFtrack1.GetDistanceFromParticle(KFtrack3);

    if (kDecayChannel == 2 && TMath::Abs(fDCA3B1XYKF) > kDCATracksCut) return 0;
    // _________________________________________________ //
    // __ daughter hypernucleus KF information __ //
    fChargeMother    = ksign * 1;
    fPDGMother       = ksign * fgkPdgCode[kPDGHyperHydrogen3];
    fDecayChannel    = kDecayChannel;
    fmSubMotherKF    = KFSubMother.GetMass();
    fmSubMotherErrKF = KFSubMother.GetErrMass();
    fESubMotherKF    = KFSubMother.GetE();
    fESubMotherErrKF = KFSubMother.GetErrE();
    fpxSubMotherKF   = KFSubMother.GetPx();
    fpxSubMotherErrKF= KFSubMother.GetErrPx();
    fpySubMotherKF   = KFSubMother.GetPy();
    fpySubMotherErrKF= KFSubMother.GetErrPy();
    fpzSubMotherKF   = KFSubMother.GetPz();
    fpzSubMotherErrKF= KFSubMother.GetErrPz();
    fptSubMotherKF   = KFSubMother.GetPt();
    fptSubMotherErrKF= KFSubMother.GetErrPt();
    fpSubMotherKF    = KFSubMother.GetP();
    fpSubMotherErrKF = KFSubMother.GetErrP();
    fySubMotherKF    = KFSubMother.GetRapidity();
    //
    fpxDaughterKF    = KFtrack1.GetPx();
    fpyDaughterKF    = KFtrack1.GetPy();
    fpzDaughterKF    = KFtrack1.GetPz();
    fptDaughterKF    = KFtrack1.GetPt();
    fEDaughterKF     = KFtrack1.GetE();
    fyDaughterKF     = KFtrack1.GetRapidity();
    fpxDaughter2KF   = KFtrack3.GetPx();
    fpyDaughter2KF   = KFtrack3.GetPy();
    fpzDaughter2KF   = KFtrack3.GetPz();
    fptDaughter2KF   = KFtrack3.GetPt();
    fEDaughter2KF    = KFtrack3.GetE();
    fyDaughter2KF    = KFtrack3.GetRapidity();
    // daughter momenta
    fpxDaughter      = particle1->Px();
    fpyDaughter      = particle1->Py();
    fpzDaughter      = particle1->Pz();
    fptDaughter      = particle1->Pt();
    fEDaughter       = particle1->E();
    fyDaughter       = particle1->Rapidity();
    fpxDaughter2     = particle3->Px();
    fpyDaughter2     = particle3->Py();
    fpzDaughter2     = particle3->Pz();
    fptDaughter2     = particle3->Pt();
    fEDaughter2      = particle3->E();
    fyDaughter2      = particle3->Rapidity();
    // __ dca prim vtx KF __ //
    fDcaDaughterXYKF     = (Float_t)KFtrack1.GetDistanceFromVertexXY(primKFVertex);
    fDcaDaughter2XYKF    = (Float_t)KFtrack3.GetDistanceFromVertexXY(primKFVertex);
    fDcaDaughterZKF      = (Float_t)KFtrack1.GetDistanceFromVertex(primKFVertex);
    fDcaDaughter2ZKF     = (Float_t)KFtrack3.GetDistanceFromVertex(primKFVertex);
    // __ dca sec vtx KF __ //
    fDcaSecDaughterXYKF  = (Float_t)KFtrack1.GetDistanceFromVertexXY(KFSubMother);
    fDcaSecDaughter2XYKF = (Float_t)KFtrack3.GetDistanceFromVertexXY(KFSubMother);
    fDcaSecDaughterZKF   = (Float_t)KFtrack1.GetDistanceFromVertex(KFSubMother);
    fDcaSecDaughter2ZKF  = (Float_t)KFtrack3.GetDistanceFromVertex(KFSubMother);
    // _________________________________________________ //
    if (kDecayChannel == 2) AliAnalysisTaskDoubleHypNucTreeLS::SetDaughterInformation("3LH", kDecayChannel);
    return 1;
  }
  // _________________________________________________ //
  else if (Mother == "4LH") {
    // _________________________________________________ //
    // __ initialize KFParticles __ //
    exTrack->Reset(); exTrack1->Reset(); exTrack2->Reset(); exTrack3->Reset(); exTrack4->Reset();
    exTrack1->CopyFromVTrack(track1);
    exTrack3->CopyFromVTrack(track3);
    KFtrack1 = AliAnalysisTaskDoubleHypNucTreeLS::CreateKFParticle(*exTrack1, AliPID::ParticleMass(AliPID::kAlpha), 2);
    KFtrack3 = AliAnalysisTaskDoubleHypNucTreeLS::CreateKFParticle(*exTrack3, AliPID::ParticleMass(AliPID::kPion),  1);
    // __ construct mother particle __ //
    checkval = 1;
    KFParticleHypNuc KFChecktrack3;
    KFChecktrack3.AddDaughter(KFtrack1);
    if (!KFChecktrack3.CheckDaughter(KFtrack3)) checkval = 0;
    KFParticleHypNuc KFChecktrack1;
    if (!KFChecktrack1.CheckDaughter(KFtrack1)) checkval = 0;
    KFChecktrack1.AddDaughter(KFtrack3);
    KFParticleHypNuc KFSubMother;
    if (checkval) {
      KFSubMother.SetConstructMethod(2);
      KFSubMother.AddDaughter(KFtrack1);
      KFSubMother.AddDaughter(KFtrack3);
      fRecoMethod = 2;
    }
    else {
      KFSubMother.AddDaughter(KFtrack1);
      KFSubMother.AddDaughter(KFtrack3);
      fRecoMethod = 0;
    }
    // __ transport mother particle to its decay vertex __ //
    KFSubMother.TransportToDecayVertex();
    // _________________________________________________ //
    // __ get decay vertex position __ //
    TertVertexKF[0]    = KFSubMother.GetX();
    TertVertexKF[1]    = KFSubMother.GetY();
    TertVertexKF[2]    = KFSubMother.GetZ();
    TertVertexErrKF[0] = KFSubMother.GetErrX();
    TertVertexErrKF[1] = KFSubMother.GetErrY();
    TertVertexErrKF[2] = KFSubMother.GetErrZ();
    // __ coord of tert vtx __ //
    fTertVertexXKF     = TertVertexKF[0];
    fTertVertexYKF     = TertVertexKF[1];
    fTertVertexZKF     = TertVertexKF[2];
    fTertVertexXErrKF  = TertVertexErrKF[0];
    fTertVertexYErrKF  = TertVertexErrKF[1];
    fTertVertexZErrKF  = TertVertexErrKF[2];
    // __ vertex quality KF __ //
    fTertVertChi2KF    = KFSubMother.Chi2();
    fTertVertNDFKF     = KFSubMother.GetNDF();
    // _________________________________________________ //
    particle1->SetXYZM(0., 0., 0., 0.);
    particle3->SetXYZM(0., 0., 0., 0.);
    sublorentzsum->SetXYZM(0., 0., 0., 0.);
    KFtertVertex     = new AliESDVertex(TertVertexKF, TertVertexErrKF);
    fPropDCADaughter = (exTrack1->PropagateToDCA(KFtertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter   = impar[0];
    fImParzDaughter  = impar[1];
    fDcaSecDaughter  = TMath::Abs(exTrack1->GetD(TertVertexKF[0], TertVertexKF[1], kMagF));
    fPropDCADaughter2= (exTrack3->PropagateToDCA(KFtertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter2  = impar[0];
    fImParzDaughter2 = impar[1];
    fDcaSecDaughter2 = TMath::Abs(exTrack3->GetD(TertVertexKF[0], TertVertexKF[1], kMagF));
    if (KFtertVertex) delete KFtertVertex;
    particle1->SetXYZM(2. * exTrack1->Px(), 2. * exTrack1->Py(), 2. * exTrack1->Pz(), AliPID::ParticleMass(AliPID::kAlpha));
    particle3->SetXYZM(exTrack3->Px(), exTrack3->Py(), exTrack3->Pz(), AliPID::ParticleMass(AliPID::kPion));
    *sublorentzsum   = *particle1 + *particle3;
    fmSubMother      = sublorentzsum->M();
    if (sublorentzsum->M() > kMax4LHMass || sublorentzsum->M() < kMin4LHMass) return 0;
    // _________________________________________________ //
    // __ KF daughter hypernucleus cos(PA) to prim vtx __ //
    dd[0] = PrimVertexKF[0] - TertVertexKF[0];
    dd[1] = PrimVertexKF[1] - TertVertexKF[1];
    dd[2] = PrimVertexKF[2] - TertVertexKF[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    fSubPAKF = TMath::Cos(sublorentzsum->Angle(*h));
    if (fSubPAKF < kPointingAngleCut) return 0;

    DecLength = TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2));
    fctSubMotherKF = (sublorentzsum->M() * DecLength) / sublorentzsum->P();
    DecLengthErr = TMath::Sqrt((fPrimVertexXErrKF * fPrimVertexXErrKF) / (fPrimVertexXKF * fPrimVertexXKF)
			       + (fPrimVertexYErrKF * fPrimVertexYErrKF) / (fPrimVertexYKF * fPrimVertexYKF)
			       + (fPrimVertexZErrKF * fPrimVertexZErrKF) / (fPrimVertexZKF * fPrimVertexZKF)
			       + (fTertVertexXErrKF * fTertVertexXErrKF) / (fTertVertexXKF * fTertVertexXKF)
			       + (fTertVertexYErrKF * fTertVertexYErrKF) / (fTertVertexYKF * fTertVertexYKF)
			       + (fTertVertexZErrKF * fTertVertexZErrKF) / (fTertVertexZKF * fTertVertexZKF)
			       + (KFSubMother.GetErrP() * KFSubMother.GetErrP()) / (KFSubMother.GetP() * KFSubMother.GetP())
			       + (KFSubMother.GetErrMass() * KFSubMother.GetErrMass()) / (KFSubMother.GetMass() * KFSubMother.GetMass())) * fctSubMotherKF;
    fctSubMotherErrKF = DecLengthErr;
    // _________________________________________________ //
    // __ dca tracks KF __ //	
    fDCA3B1XYKF = (Float_t)KFtrack1.GetDistanceFromParticleXY(KFtrack3);
    fDCA3B1ZKF  = (Float_t)KFtrack1.GetDistanceFromParticle(KFtrack3);

    if (TMath::Abs(fDCA3B1XYKF) > kDCATracksCut) return 0;
    // _________________________________________________ //
    // __ daughter hypernucleus KF information __ //
    fChargeMother     = ksign * 1;
    fPDGMother        = ksign * fgkPdgCode[kPDGHyperHydrogen4];
    fDecayChannel     = kDecayChannel;
    fmSubMotherKF     = KFSubMother.GetMass();
    fmSubMotherErrKF  = KFSubMother.GetErrMass();
    fESubMotherKF     = KFSubMother.GetE();
    fESubMotherErrKF  = KFSubMother.GetErrE();
    fpxSubMotherKF    = KFSubMother.GetPx();
    fpxSubMotherErrKF = KFSubMother.GetErrPx();
    fpySubMotherKF    = KFSubMother.GetPy();
    fpySubMotherErrKF = KFSubMother.GetErrPy();
    fpzSubMotherKF    = KFSubMother.GetPz();
    fpzSubMotherErrKF = KFSubMother.GetErrPz();
    fptSubMotherKF    = KFSubMother.GetPt();
    fptSubMotherErrKF = KFSubMother.GetErrPt();
    fpSubMotherKF     = KFSubMother.GetP();
    fpSubMotherErrKF  = KFSubMother.GetErrP();
    fySubMotherKF     = KFSubMother.GetRapidity();
    //
    fpxDaughterKF     = KFtrack1.GetPx();
    fpyDaughterKF     = KFtrack1.GetPy();
    fpzDaughterKF     = KFtrack1.GetPz();
    fptDaughterKF     = KFtrack1.GetPt();
    fEDaughterKF      = KFtrack1.GetE();
    fyDaughterKF      = KFtrack1.GetRapidity();
    fpxDaughter2KF    = KFtrack3.GetPx();
    fpyDaughter2KF    = KFtrack3.GetPy();
    fpzDaughter2KF    = KFtrack3.GetPz();
    fptDaughter2KF    = KFtrack3.GetPt();
    fEDaughter2KF     = KFtrack3.GetE();
    fyDaughter2KF     = KFtrack3.GetRapidity();
    // daughter momenta
    fpxDaughter       = particle1->Px();
    fpyDaughter       = particle1->Py();
    fpzDaughter       = particle1->Pz();
    fptDaughter       = particle1->Pt();
    fEDaughter        = particle1->E();
    fyDaughter        = particle1->Rapidity();
    fpxDaughter2      = particle3->Px();
    fpyDaughter2      = particle3->Py();
    fpzDaughter2      = particle3->Pz();
    fptDaughter2      = particle3->Pt();
    fEDaughter2       = particle3->E();
    fyDaughter2       = particle3->Rapidity();
    // __ dca prim vtx KF __ //
    fDcaDaughterXYKF     = (Float_t)KFtrack1.GetDistanceFromVertexXY(primKFVertex);
    fDcaDaughter2XYKF    = (Float_t)KFtrack3.GetDistanceFromVertexXY(primKFVertex);
    fDcaDaughterZKF      = (Float_t)KFtrack1.GetDistanceFromVertex(primKFVertex);
    fDcaDaughter2ZKF     = (Float_t)KFtrack3.GetDistanceFromVertex(primKFVertex);
    // __ dca sec vtx KF __ //
    fDcaSecDaughterXYKF  = (Float_t)KFtrack1.GetDistanceFromVertexXY(KFSubMother);
    fDcaSecDaughter2XYKF = (Float_t)KFtrack3.GetDistanceFromVertexXY(KFSubMother);
    fDcaSecDaughterZKF   = (Float_t)KFtrack1.GetDistanceFromVertex(KFSubMother);
    fDcaSecDaughter2ZKF  = (Float_t)KFtrack3.GetDistanceFromVertex(KFSubMother);
    // _________________________________________________ //
    AliAnalysisTaskDoubleHypNucTreeLS::SetDaughterInformation("4LH", kDecayChannel);
    return 1;
  }
  else if (Mother == "4LHe") {
    // _________________________________________________ //    
    exTrack->Reset(); exTrack1->Reset(); exTrack2->Reset(); exTrack3->Reset(); exTrack4->Reset();
    exTrack1->CopyFromVTrack(track1);
    exTrack2->CopyFromVTrack(track2);
    exTrack3->CopyFromVTrack(track3);
    // __ initialize KFParticles __ //
    KFtrack1 = AliAnalysisTaskDoubleHypNucTreeLS::CreateKFParticle(*exTrack1, AliPID::ParticleMass(AliPID::kHe3),    2);
    KFtrack2 = AliAnalysisTaskDoubleHypNucTreeLS::CreateKFParticle(*exTrack2, AliPID::ParticleMass(AliPID::kProton), 1);
    KFtrack3 = AliAnalysisTaskDoubleHypNucTreeLS::CreateKFParticle(*exTrack3, AliPID::ParticleMass(AliPID::kPion),   1);
    // __ construct mother particle __ //
    checkval = 1;
    KFParticleHypNuc KFChecktrack3;
    KFChecktrack3.AddDaughter(KFtrack1);
    KFChecktrack3.AddDaughter(KFtrack2);
    if (!KFChecktrack3.CheckDaughter(KFtrack3)) checkval = 0;
    KFParticleHypNuc KFChecktrack2;
    KFChecktrack2.AddDaughter(KFtrack1);
    if (!KFChecktrack2.CheckDaughter(KFtrack2)) checkval = 0;
    KFChecktrack2.AddDaughter(KFtrack3);
    KFParticleHypNuc KFChecktrack1;
    if (!KFChecktrack1.CheckDaughter(KFtrack1)) checkval = 0;
    KFChecktrack1.AddDaughter(KFtrack2);
    KFChecktrack1.AddDaughter(KFtrack3);
    KFParticleHypNuc KFSubMother;
    if (checkval) {
      KFSubMother.SetConstructMethod(2);
      KFSubMother.AddDaughter(KFtrack1);
      KFSubMother.AddDaughter(KFtrack2);
      KFSubMother.AddDaughter(KFtrack3);
      fRecoMethod = 2;
    }
    else {
      KFSubMother.AddDaughter(KFtrack1);
      KFSubMother.AddDaughter(KFtrack2);
      KFSubMother.AddDaughter(KFtrack3);
      fRecoMethod = 0;
    }
    // __ transport mother particle to its decay vertex __ //
    KFSubMother.TransportToDecayVertex();
    // _________________________________________________ //
    // __ get decay vertex position __ //
    TertVertexKF[0]    = KFSubMother.GetX();
    TertVertexKF[1]    = KFSubMother.GetY();
    TertVertexKF[2]    = KFSubMother.GetZ();
    TertVertexErrKF[0] = KFSubMother.GetErrX();
    TertVertexErrKF[1] = KFSubMother.GetErrY();
    TertVertexErrKF[2] = KFSubMother.GetErrZ();
    // __ coord of tert vtx __ //
    fTertVertexXKF     = TertVertexKF[0];
    fTertVertexYKF     = TertVertexKF[1];
    fTertVertexZKF     = TertVertexKF[2];
    fTertVertexXErrKF  = TertVertexErrKF[0];
    fTertVertexYErrKF  = TertVertexErrKF[1];
    fTertVertexZErrKF  = TertVertexErrKF[2];
    // __ vertex quality KF __ //
    fTertVertChi2KF    = KFSubMother.Chi2();
    fTertVertNDFKF     = KFSubMother.GetNDF();
    // _________________________________________________ //
    particle1->SetXYZM(0., 0., 0., 0.);
    particle2->SetXYZM(0., 0., 0., 0.);
    particle3->SetXYZM(0., 0., 0., 0.);
    sublorentzsum->SetXYZM(0., 0., 0., 0.);
    KFtertVertex      = new AliESDVertex(TertVertexKF, TertVertexErrKF);
    fPropDCADaughter  = (exTrack1->PropagateToDCA(KFtertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter    = impar[0];
    fImParzDaughter   = impar[1];
    fDcaSecDaughter   = TMath::Abs(exTrack1->GetD(TertVertexKF[0], TertVertexKF[1], kMagF));
    fPropDCADaughter1 = (exTrack2->PropagateToDCA(KFtertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter1   = impar[0];
    fImParzDaughter1  = impar[1];
    fDcaSecDaughter1  = TMath::Abs(exTrack2->GetD(TertVertexKF[0], TertVertexKF[1], kMagF));
    fPropDCADaughter2 = (exTrack3->PropagateToDCA(KFtertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter2   = impar[0];
    fImParzDaughter2  = impar[1];
    fDcaSecDaughter2  = TMath::Abs(exTrack3->GetD(TertVertexKF[0], TertVertexKF[1], kMagF));
    if (KFtertVertex) delete KFtertVertex;
    particle1->SetXYZM(2. * exTrack1->Px(), 2. * exTrack1->Py(), 2. * exTrack1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
    particle2->SetXYZM(exTrack2->Px(), exTrack2->Py(), exTrack2->Pz(), AliPID::ParticleMass(AliPID::kProton));
    particle3->SetXYZM(exTrack3->Px(), exTrack3->Py(), exTrack3->Pz(), AliPID::ParticleMass(AliPID::kPion));
    *sublorentzsum   = *particle1 + *particle2 + *particle3;
    fmSubMother      = sublorentzsum->M();
    if (sublorentzsum->M() > kMax4LHeMass || sublorentzsum->M() < kMin4LHeMass) return 0;
    // _________________________________________________ //
    // __ KF daughter hypernucleus cos(PA) to prim vtx __ //
    dd[0] = PrimVertexKF[0] - TertVertexKF[0];
    dd[1] = PrimVertexKF[1] - TertVertexKF[1];
    dd[2] = PrimVertexKF[2] - TertVertexKF[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    fSubPAKF = TMath::Cos(sublorentzsum->Angle(*h));
    if (fSubPAKF < kPointingAngleCut) return 0;

    DecLength = TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2));
    fctSubMotherKF = (sublorentzsum->M() * DecLength) / sublorentzsum->P();
    DecLengthErr = TMath::Sqrt((fPrimVertexXErrKF * fPrimVertexXErrKF) / (fPrimVertexXKF * fPrimVertexXKF)
			       + (fPrimVertexYErrKF * fPrimVertexYErrKF) / (fPrimVertexYKF * fPrimVertexYKF)
			       + (fPrimVertexZErrKF * fPrimVertexZErrKF) / (fPrimVertexZKF * fPrimVertexZKF)
			       + (fTertVertexXErrKF * fTertVertexXErrKF) / (fTertVertexXKF * fTertVertexXKF)
			       + (fTertVertexYErrKF * fTertVertexYErrKF) / (fTertVertexYKF * fTertVertexYKF)
			       + (fTertVertexZErrKF * fTertVertexZErrKF) / (fTertVertexZKF * fTertVertexZKF)
			       + (KFSubMother.GetErrP() * KFSubMother.GetErrP()) / (KFSubMother.GetP() * KFSubMother.GetP())
			       + (KFSubMother.GetErrMass() * KFSubMother.GetErrMass()) / (KFSubMother.GetMass() * KFSubMother.GetMass())) * fctSubMotherKF;
    fctSubMotherErrKF = DecLengthErr;
    // _________________________________________________ //
    // __ dca tracks KF __ //	
    fDCA3B2XYKF = (Float_t)KFtrack1.GetDistanceFromParticleXY(KFtrack3);
    fDCA3B3XYKF = (Float_t)KFtrack2.GetDistanceFromParticleXY(KFtrack3);
    fDCA3B1XYKF = (Float_t)KFtrack2.GetDistanceFromParticleXY(KFtrack1);
    fDCA3B2ZKF  = (Float_t)KFtrack1.GetDistanceFromParticle(KFtrack3);
    fDCA3B3ZKF  = (Float_t)KFtrack2.GetDistanceFromParticle(KFtrack3);
    fDCA3B1ZKF  = (Float_t)KFtrack2.GetDistanceFromParticle(KFtrack1);

    if (TMath::Abs(fDCA3B1XYKF) > kDCATracksCut || TMath::Abs(fDCA3B2XYKF) > kDCATracksCut
	|| TMath::Abs(fDCA3B3XYKF) > kDCATracksCut) return 0;
    // _________________________________________________ //
    // __ daughter hypernucleus KF information __ //
    fChargeMother     = ksign * 2;
    fDecayChannel     = kDecayChannel;
    fPDGMother        = ksign * fgkPdgCode[kPDGHyperHelium4];
    fmSubMotherKF     = KFSubMother.GetMass();
    fmSubMotherErrKF  = KFSubMother.GetErrMass();
    fESubMotherKF     = KFSubMother.GetE();
    fESubMotherErrKF  = KFSubMother.GetErrE();
    fpxSubMotherKF    = KFSubMother.GetPx();
    fpxSubMotherErrKF = KFSubMother.GetErrPx();
    fpySubMotherKF    = KFSubMother.GetPy();
    fpySubMotherErrKF = KFSubMother.GetErrPy();
    fpzSubMotherKF    = KFSubMother.GetPz();
    fpzSubMotherErrKF = KFSubMother.GetErrPz();
    fptSubMotherKF    = KFSubMother.GetPt();
    fptSubMotherErrKF = KFSubMother.GetErrPt();
    fpSubMotherKF     = KFSubMother.GetP();
    fpSubMotherErrKF  = KFSubMother.GetErrP();
    fySubMotherKF     = KFSubMother.GetRapidity();
    //
    fpxDaughterKF     = KFtrack1.GetPx();
    fpyDaughterKF     = KFtrack1.GetPy();
    fpzDaughterKF     = KFtrack1.GetPz();
    fptDaughterKF     = KFtrack1.GetPt();
    fEDaughterKF      = KFtrack1.GetE();
    fyDaughterKF      = KFtrack1.GetRapidity();
    fpxDaughter1KF    = KFtrack2.GetPx();
    fpyDaughter1KF    = KFtrack2.GetPy();
    fpzDaughter1KF    = KFtrack2.GetPz();
    fptDaughter1KF    = KFtrack2.GetPt();
    fEDaughter1KF     = KFtrack2.GetE();
    fyDaughter1KF     = KFtrack2.GetRapidity();
    fpxDaughter2KF    = KFtrack3.GetPx();
    fpyDaughter2KF    = KFtrack3.GetPy();
    fpzDaughter2KF    = KFtrack3.GetPz();
    fptDaughter2KF    = KFtrack3.GetPt();
    fEDaughter2KF     = KFtrack3.GetE();
    fyDaughter2KF     = KFtrack3.GetRapidity();
    // daughter momenta
    fpxDaughter       = particle1->Px();
    fpyDaughter       = particle1->Py();
    fpzDaughter       = particle1->Pz();
    fptDaughter       = particle1->Pt();
    fEDaughter        = particle1->E();
    fyDaughter        = particle1->Rapidity();
    fpxDaughter1      = particle2->Px();
    fpyDaughter1      = particle2->Py();
    fpzDaughter1      = particle2->Pz();
    fptDaughter1      = particle2->Pt();
    fEDaughter1       = particle2->E();
    fyDaughter1       = particle2->Rapidity();
    fpxDaughter2      = particle3->Px();
    fpyDaughter2      = particle3->Py();
    fpzDaughter2      = particle3->Pz();
    fptDaughter2      = particle3->Pt();
    fEDaughter2       = particle3->E();
    fyDaughter2       = particle3->Rapidity();
    // __ dca prim vtx KF __ //
    fDcaDaughterXYKF     = (Float_t)KFtrack1.GetDistanceFromVertexXY(primKFVertex);
    fDcaDaughter1XYKF    = (Float_t)KFtrack2.GetDistanceFromVertexXY(primKFVertex);
    fDcaDaughter2XYKF    = (Float_t)KFtrack3.GetDistanceFromVertexXY(primKFVertex);
    fDcaDaughterZKF      = (Float_t)KFtrack1.GetDistanceFromVertex(primKFVertex);
    fDcaDaughter1ZKF     = (Float_t)KFtrack2.GetDistanceFromVertex(primKFVertex);
    fDcaDaughter2ZKF     = (Float_t)KFtrack3.GetDistanceFromVertex(primKFVertex);
    // __ dca sec vtx KF __ //
    fDcaSecDaughterXYKF  = (Float_t)KFtrack1.GetDistanceFromVertexXY(KFSubMother);
    fDcaSecDaughter1XYKF = (Float_t)KFtrack2.GetDistanceFromVertexXY(KFSubMother);
    fDcaSecDaughter2XYKF = (Float_t)KFtrack3.GetDistanceFromVertexXY(KFSubMother);
    fDcaSecDaughterZKF   = (Float_t)KFtrack1.GetDistanceFromVertex(KFSubMother);
    fDcaSecDaughter1ZKF  = (Float_t)KFtrack2.GetDistanceFromVertex(KFSubMother);
    fDcaSecDaughter2ZKF  = (Float_t)KFtrack3.GetDistanceFromVertex(KFSubMother);
    // _________________________________________________ //
    AliAnalysisTaskDoubleHypNucTreeLS::SetDaughterInformation("4LHe", kDecayChannel);
    return 1;
  }
  // _________________________________________________ //
  else if (Mother == "4LLH" && kDecayChannel == 1) {
    // _________________________________________________ //
    exTrack->Reset(); exTrack1->Reset(); exTrack2->Reset(); exTrack3->Reset(); exTrack4->Reset();
    exTrack1->CopyFromVTrack(track1);
    exTrack2->CopyFromVTrack(track2);
    exTrack3->CopyFromVTrack(track3);
    exTrack4->CopyFromVTrack(track4);
    // __ initialize KFParticles __ //
    KFtrack1 = AliAnalysisTaskDoubleHypNucTreeLS::CreateKFParticle(*exTrack1, AliPID::ParticleMass(AliPID::kHe3),    2);
    KFtrack2 = AliAnalysisTaskDoubleHypNucTreeLS::CreateKFParticle(*exTrack2, AliPID::ParticleMass(AliPID::kProton), 1);
    KFtrack3 = AliAnalysisTaskDoubleHypNucTreeLS::CreateKFParticle(*exTrack3, AliPID::ParticleMass(AliPID::kPion),   1);
    KFtrack4 = AliAnalysisTaskDoubleHypNucTreeLS::CreateKFParticle(*exTrack4, AliPID::ParticleMass(AliPID::kPion),   1);
    // __ construct mother particle __ //
    checkval = 1;
    KFParticleHypNuc KFChecktrack3;
    KFChecktrack3.AddDaughter(KFtrack1);
    KFChecktrack3.AddDaughter(KFtrack2);
    if (!KFChecktrack3.CheckDaughter(KFtrack3)) checkval = 0;
    KFParticleHypNuc KFChecktrack2;
    KFChecktrack2.AddDaughter(KFtrack1);
    if (!KFChecktrack2.CheckDaughter(KFtrack2)) checkval = 0;
    KFChecktrack2.AddDaughter(KFtrack3);
    KFParticleHypNuc KFChecktrack1;
    if (!KFChecktrack1.CheckDaughter(KFtrack1)) checkval = 0;
    KFChecktrack1.AddDaughter(KFtrack2);
    KFChecktrack1.AddDaughter(KFtrack3);

    KFParticleHypNuc KFSubMother;
    if (checkval) {
      KFSubMother.SetConstructMethod(2);
      KFSubMother.AddDaughter(KFtrack1);
      KFSubMother.AddDaughter(KFtrack2);
      KFSubMother.AddDaughter(KFtrack3);
      fRecoMethod = 2;
    }
    else {
      KFSubMother.AddDaughter(KFtrack1);
      KFSubMother.AddDaughter(KFtrack2);
      KFSubMother.AddDaughter(KFtrack3);
      fRecoMethod = 0;
    }
    checkval = 1;
    KFParticleHypNuc KFChecktrack4;
    if (!KFChecktrack4.CheckDaughter(KFtrack4)) checkval = 0;
    KFChecktrack4.AddDaughter(KFSubMother);
    KFParticleHypNuc KFChecktrack5;
    if (!KFChecktrack5.CheckDaughter(KFSubMother)) checkval = 0;
    KFChecktrack5.AddDaughter(KFtrack4);

    KFParticleHypNuc KFMother;
    if (checkval) {
      KFMother.SetConstructMethod(2);
      KFMother.AddDaughter(KFSubMother);
      KFMother.AddDaughter(KFtrack4);
      fRecoMethod = 2;
    }
    else {
      KFMother.AddDaughter(KFSubMother);
      KFMother.AddDaughter(KFtrack4);
      fRecoMethod = 0;
    }
    // _________________________________________________ //
    // __ transport mother particle to its decay vertex __ //    
    KFSubMother.TransportToDecayVertex();
    KFMother.TransportToDecayVertex();
    // _________________________________________________ //
    // __ get decay vertex position __ //
    SecVertexKF[0]     = KFMother.GetX();
    SecVertexKF[1]     = KFMother.GetY();
    SecVertexKF[2]     = KFMother.GetZ();
    SecVertexErrKF[0]  = KFMother.GetErrX();
    SecVertexErrKF[1]  = KFMother.GetErrY();
    SecVertexErrKF[2]  = KFMother.GetErrZ();
    // __ coord of sec vtx __ //
    fSecVertexXKF      = SecVertexKF[0];
    fSecVertexYKF      = SecVertexKF[1];
    fSecVertexZKF      = SecVertexKF[2];
    fSecVertexXErrKF   = SecVertexErrKF[0];
    fSecVertexYErrKF   = SecVertexErrKF[1];
    fSecVertexZErrKF   = SecVertexErrKF[2];
    // __ vertex quality KF __ //
    fSecVertChi2KF     = KFMother.Chi2();
    fSecVertNDFKF      = KFMother.GetNDF();
    // __ get decay vertex position __ //
    TertVertexKF[0]    = KFSubMother.GetX();
    TertVertexKF[1]    = KFSubMother.GetY();
    TertVertexKF[2]    = KFSubMother.GetZ();
    TertVertexErrKF[0] = KFSubMother.GetErrX();
    TertVertexErrKF[1] = KFSubMother.GetErrY();
    TertVertexErrKF[2] = KFSubMother.GetErrZ();
    // __ coord of tert vtx __ //
    fTertVertexXKF     = TertVertexKF[0];
    fTertVertexYKF     = TertVertexKF[1];
    fTertVertexZKF     = TertVertexKF[2];
    fTertVertexXErrKF  = TertVertexErrKF[0];
    fTertVertexYErrKF  = TertVertexErrKF[1];
    fTertVertexZErrKF  = TertVertexErrKF[2];
    // __ vertex quality KF __ //
    fTertVertChi2KF    = KFSubMother.Chi2();
    fTertVertNDFKF     = KFSubMother.GetNDF();
    // __ create SubMother AliExternalTrackParam __ //    
    Double_t xyzV[3];    for (int indx = 0; indx < 3; indx++)  xyzV[indx]    = KFSubMother.GetParameter(indx);
    Double_t pxpypzV[3]; for (int indx = 3; indx < 6; indx++)  pxpypzV[indx] = KFSubMother.GetParameter(indx);
    Double_t covmV[21];  for (int indx = 0; indx < 21; indx++) covmV[indx]   = KFSubMother.GetCovariance(indx);
    exTrack = new AliExternalTrackParam(xyzV, pxpypzV, covmV, ksign);
    exTrack->GetCovarianceXYZPxPyPz(fCovMatrixTrack4);
    fTrackPar4[0] = exTrack->GetX();   fTrackPar4[1] = exTrack->GetY();         fTrackPar4[2] = exTrack->GetZ();     fTrackPar4[3] = exTrack->GetSnp();
    fTrackPar4[4] = exTrack->GetTgl(); fTrackPar4[5] = exTrack->GetSigned1Pt(); fTrackPar4[6] = exTrack->GetAlpha();
    // _________________________________________________ //
    particle1->SetXYZM(0., 0., 0., 0.);
    particle2->SetXYZM(0., 0., 0., 0.);
    particle3->SetXYZM(0., 0., 0., 0.);
    particle4->SetXYZM(0., 0., 0., 0.);
    sublorentzsum->SetXYZM(0., 0., 0., 0.);
    lorentzsum->SetXYZM(0., 0., 0., 0.);
    KFtertVertex      = new AliESDVertex(TertVertexKF, TertVertexErrKF);
    fPropDCADaughter  = (exTrack1->PropagateToDCA(KFtertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter    = impar[0];
    fImParzDaughter   = impar[1];
    fDcaSecDaughter   = TMath::Abs(exTrack1->GetD(TertVertexKF[0], TertVertexKF[1], kMagF));
    fPropDCADaughter1 = (exTrack2->PropagateToDCA(KFtertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter1   = impar[0];
    fImParzDaughter1  = impar[1];
    fDcaSecDaughter1  = TMath::Abs(exTrack2->GetD(TertVertexKF[0], TertVertexKF[1], kMagF));
    fPropDCADaughter2 = (exTrack3->PropagateToDCA(KFtertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter2   = impar[0];
    fImParzDaughter2  = impar[1];
    fDcaSecDaughter2  = TMath::Abs(exTrack3->GetD(TertVertexKF[0], TertVertexKF[1], kMagF));
    if (KFtertVertex) delete KFtertVertex;
    KFsecVertex       = new AliESDVertex(SecVertexKF, SecVertexErrKF);
    fPropDCADaughter4 = (exTrack->PropagateToDCA(KFsecVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter4   = impar[0];
    fImParzDaughter4  = impar[1];
    fDcaSecDaughter4  = TMath::Abs(exTrack->GetD(SecVertexKF[0], SecVertexKF[1], kMagF));
    fPropDCADaughter3 = (exTrack4->PropagateToDCA(KFsecVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter3   = impar[0];
    fImParzDaughter3  = impar[1];
    fDcaSecDaughter3  = TMath::Abs(exTrack4->GetD(SecVertexKF[0], SecVertexKF[1], kMagF));
    if (KFsecVertex) delete KFsecVertex;
    particle1->SetXYZM(2. * exTrack1->Px(), 2. * exTrack1->Py(), 2. * exTrack1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
    particle2->SetXYZM(exTrack2->Px(), exTrack2->Py(), exTrack2->Pz(), AliPID::ParticleMass(AliPID::kProton));
    particle3->SetXYZM(exTrack3->Px(), exTrack3->Py(), exTrack3->Pz(), AliPID::ParticleMass(AliPID::kPion));
    particle4->SetXYZM(exTrack4->Px(), exTrack4->Py(), exTrack4->Pz(), AliPID::ParticleMass(AliPID::kPion));
    *sublorentzsum   = *particle1 + *particle2 + *particle3;
    *lorentzsum      = *particle1 + *particle2 + *particle3 + *particle4;
    fmSubMother      = sublorentzsum->M();
    fmMother         = lorentzsum->M();
    if (lorentzsum->M() > kMax4LLHMass || lorentzsum->M() < kMin4LLHMass) return 0;
    // _________________________________________________ //
    // __ mother hypernucleus cos(PA) __ //
    dd[0] = PrimVertexKF[0] - SecVertexKF[0];
    dd[1] = PrimVertexKF[1] - SecVertexKF[1];
    dd[2] = PrimVertexKF[2] - SecVertexKF[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    fPAKF = TMath::Cos(lorentzsum->Angle(*h));
    //if (fPAKF < kPointingAngleCut) return 0;

    DecLength = TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2));
    fctMotherKF = (lorentzsum->M() * DecLength) / lorentzsum->P();
    DecLengthErr = TMath::Sqrt((fPrimVertexXErrKF * fPrimVertexXErrKF) / (fPrimVertexXKF * fPrimVertexXKF)
			       + (fPrimVertexYErrKF * fPrimVertexYErrKF) / (fPrimVertexYKF * fPrimVertexYKF)
			       + (fPrimVertexZErrKF * fPrimVertexZErrKF) / (fPrimVertexZKF * fPrimVertexZKF)
			       + (fSecVertexXErrKF * fSecVertexXErrKF) / (fSecVertexXKF * fSecVertexXKF)
			       + (fSecVertexYErrKF * fSecVertexYErrKF) / (fSecVertexYKF * fSecVertexYKF)
			       + (fSecVertexZErrKF * fSecVertexZErrKF) / (fSecVertexZKF * fSecVertexZKF)
			       + (KFMother.GetErrP() * KFMother.GetErrP()) / (KFMother.GetP() * KFMother.GetP())
			       + (KFMother.GetErrMass() * KFMother.GetErrMass()) / (KFMother.GetMass() * KFMother.GetMass())) * fctMotherKF;
    fctMotherErrKF = DecLengthErr;
    // _________________________________________________ //    
    // __ daughter hypernucleus cos(PA) to sec vtx __ //
    dd[0] = SecVertexKF[0] - TertVertexKF[0];
    dd[1] = SecVertexKF[1] - TertVertexKF[1];
    dd[2] = SecVertexKF[2] - TertVertexKF[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    fSubPAKF = TMath::Cos(sublorentzsum->Angle(*h));
    if (fSubPAKF < kPointingAngleCut) return 0;

    DecLength = TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2));
    fctSubMotherKF = (sublorentzsum->M() * DecLength) / sublorentzsum->P();
    DecLengthErr = TMath::Sqrt((fSecVertexXErrKF * fSecVertexXErrKF) / (fSecVertexXKF * fSecVertexXKF)
			       + (fSecVertexYErrKF * fSecVertexYErrKF) / (fSecVertexYKF * fSecVertexYKF)
			       + (fSecVertexZErrKF * fSecVertexZErrKF) / (fSecVertexZKF * fSecVertexZKF)
			       + (fTertVertexXErrKF * fTertVertexXErrKF) / (fTertVertexXKF * fTertVertexXKF)
			       + (fTertVertexYErrKF * fTertVertexYErrKF) / (fTertVertexYKF * fTertVertexYKF)
			       + (fTertVertexZErrKF * fTertVertexZErrKF) / (fTertVertexZKF * fTertVertexZKF)
			       + (KFSubMother.GetErrP() * KFSubMother.GetErrP()) / (KFSubMother.GetP() * KFSubMother.GetP())
			       + (KFSubMother.GetErrMass() * KFSubMother.GetErrMass()) / (KFSubMother.GetMass() * KFSubMother.GetMass())) * fctSubMotherKF;
    fctSubMotherErrKF = DecLengthErr;
    // _________________________________________________ //
    // __ KF daughter hypernucleus cos(PA) to prim vtx __ //
    dd[0] = PrimVertexKF[0] - TertVertexKF[0];
    dd[1] = PrimVertexKF[1] - TertVertexKF[1];
    dd[2] = PrimVertexKF[2] - TertVertexKF[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    fSubPA2KF = TMath::Cos(sublorentzsum->Angle(*h));
    if (fSubPA2KF < kPointingAngleCut) return 0;
    // _________________________________________________ //
    // __ dca tracks KF __ //
    fDCA2BXYKF  = (Float_t)KFSubMother.GetDistanceFromParticleXY(KFtrack4);
    fDCA3B1XYKF = (Float_t)KFtrack1.GetDistanceFromParticleXY(KFtrack2);
    fDCA3B2XYKF = (Float_t)KFtrack1.GetDistanceFromParticleXY(KFtrack3);
    fDCA3B3XYKF = (Float_t)KFtrack2.GetDistanceFromParticleXY(KFtrack3);
    fDCA2BZKF   = (Float_t)KFSubMother.GetDistanceFromParticle(KFtrack4);
    fDCA3B1ZKF  = (Float_t)KFtrack1.GetDistanceFromParticle(KFtrack2);
    fDCA3B2ZKF  = (Float_t)KFtrack1.GetDistanceFromParticle(KFtrack3);
    fDCA3B3ZKF  = (Float_t)KFtrack2.GetDistanceFromParticle(KFtrack3);

    if (TMath::Abs(fDCA3B1XYKF) > kDCATracksCut || TMath::Abs(fDCA3B2XYKF) > kDCATracksCut
	|| TMath::Abs(fDCA3B3XYKF) > kDCATracksCut || TMath::Abs(fDCA2BXYKF) > kDCATracksCut) return 0;
    // _________________________________________________ //
    // __ mother hypernucleus KF information __ //
    fChargeMother    = ksign * 1;
    fDecayChannel    = kDecayChannel;
    fPDGMother       = ksign * fgkPdgCode[kPDGDoubleHyperHydrogen4];
    fmMotherKF       = KFMother.GetMass();
    fmMotherErrKF    = KFMother.GetErrMass();
    fEMotherKF       = KFMother.GetE();
    fEMotherErrKF    = KFMother.GetErrE();
    fpxMotherKF      = KFMother.GetPx();
    fpxMotherErrKF   = KFMother.GetErrPx();
    fpyMotherKF      = KFMother.GetPy();
    fpyMotherErrKF   = KFMother.GetErrPy();
    fpzMotherKF      = KFMother.GetPz();
    fpzMotherErrKF   = KFMother.GetErrPz();
    fptMotherKF      = KFMother.GetPt();
    fptMotherErrKF   = KFMother.GetErrPt();
    fpMotherKF       = KFMother.GetP();
    fpMotherErrKF    = KFMother.GetErrP();
    fyMotherKF       = KFSubMother.GetRapidity();
    // __ daughter hypernucleus KF information __ //
    fmSubMotherKF    = KFSubMother.GetMass();
    fmSubMotherErrKF = KFSubMother.GetErrMass();
    fESubMotherKF    = KFSubMother.GetE();
    fESubMotherErrKF = KFSubMother.GetErrE();
    fpxSubMotherKF   = KFSubMother.GetPx();
    fpxSubMotherErrKF= KFSubMother.GetErrPx();
    fpySubMotherKF   = KFSubMother.GetPy();
    fpySubMotherErrKF= KFSubMother.GetErrPy();
    fpzSubMotherKF   = KFSubMother.GetPz();
    fpzSubMotherErrKF= KFSubMother.GetErrPz();
    fptSubMotherKF   = KFSubMother.GetPt();
    fptSubMotherErrKF= KFSubMother.GetErrPt();
    fpSubMotherKF    = KFSubMother.GetP();
    fpSubMotherErrKF = KFSubMother.GetErrP();
    fySubMotherKF    = KFSubMother.GetRapidity();
    //
    fpxDaughterKF    = KFtrack1.GetPx();
    fpyDaughterKF    = KFtrack1.GetPy();
    fpzDaughterKF    = KFtrack1.GetPz();
    fptDaughterKF    = KFtrack1.GetPt();
    fEDaughterKF     = KFtrack1.GetE();
    fyDaughterKF     = KFtrack1.GetRapidity();
    fpxDaughter1KF   = KFtrack2.GetPx();
    fpyDaughter1KF   = KFtrack2.GetPy();
    fpzDaughter1KF   = KFtrack2.GetPz();
    fptDaughter1KF   = KFtrack2.GetPt();
    fEDaughter1KF    = KFtrack2.GetE();
    fyDaughter1KF    = KFtrack2.GetRapidity();
    fpxDaughter2KF   = KFtrack3.GetPx();
    fpyDaughter2KF   = KFtrack3.GetPy();
    fpzDaughter2KF   = KFtrack3.GetPz();
    fptDaughter2KF   = KFtrack3.GetPt();
    fEDaughter2KF    = KFtrack3.GetE();
    fyDaughter2KF    = KFtrack3.GetRapidity();
    fpxDaughter3KF   = KFtrack4.GetPx();
    fpyDaughter3KF   = KFtrack4.GetPy();
    fpzDaughter3KF   = KFtrack4.GetPz();
    fptDaughter3KF   = KFtrack4.GetPt();
    fEDaughter3KF    = KFtrack4.GetE();
    fyDaughter3KF    = KFtrack4.GetRapidity();
    // daughter momenta
    fpxDaughter      = particle1->Px();
    fpyDaughter      = particle1->Py();
    fpzDaughter      = particle1->Pz();
    fptDaughter      = particle1->Pt();
    fEDaughter       = particle1->E();
    fyDaughter       = particle1->Rapidity();
    fpxDaughter1     = particle2->Px();
    fpyDaughter1     = particle2->Py();
    fpzDaughter1     = particle2->Pz();
    fptDaughter1     = particle2->Pt();
    fEDaughter1      = particle2->E();
    fyDaughter1      = particle2->Rapidity();
    fpxDaughter2     = particle3->Px();
    fpyDaughter2     = particle3->Py();
    fpzDaughter2     = particle3->Pz();
    fptDaughter2     = particle3->Pt();
    fEDaughter2      = particle3->E();
    fyDaughter2      = particle3->Rapidity();
    fpxDaughter3     = particle4->Px();
    fpyDaughter3     = particle4->Py();
    fpzDaughter3     = particle4->Pz();
    fptDaughter3     = particle4->Pt();
    fEDaughter3      = particle4->E();
    fyDaughter3      = particle4->Rapidity();
    // __ Check for V0s __ //
    //AliAnalysisTaskDoubleHypNucTreeLS::GetV0Status();
    // __ dca prim vtx KF __ //
    fDcaDaughterXYKF     = (Float_t)KFtrack1.GetDistanceFromVertexXY(primKFVertex);
    fDcaDaughter1XYKF    = (Float_t)KFtrack2.GetDistanceFromVertexXY(primKFVertex);
    fDcaDaughter2XYKF    = (Float_t)KFtrack3.GetDistanceFromVertexXY(primKFVertex);
    fDcaDaughter3XYKF    = (Float_t)KFtrack4.GetDistanceFromVertexXY(primKFVertex);
    fDcaDaughter4XYKF    = (Float_t)KFSubMother.GetDistanceFromVertexXY(primKFVertex);
    fDcaDaughterZKF      = (Float_t)KFtrack1.GetDistanceFromVertex(primKFVertex);
    fDcaDaughter1ZKF     = (Float_t)KFtrack2.GetDistanceFromVertex(primKFVertex);
    fDcaDaughter2ZKF     = (Float_t)KFtrack3.GetDistanceFromVertex(primKFVertex);
    fDcaDaughter3ZKF     = (Float_t)KFtrack4.GetDistanceFromVertex(primKFVertex);
    fDcaDaughter4ZKF     = (Float_t)KFSubMother.GetDistanceFromVertex(primKFVertex);
    // __ dca sec vtx KF __ //
    fDcaSecDaughterXYKF  = (Float_t)KFtrack1.GetDistanceFromVertexXY(KFSubMother);
    fDcaSecDaughter1XYKF = (Float_t)KFtrack2.GetDistanceFromVertexXY(KFSubMother);
    fDcaSecDaughter2XYKF = (Float_t)KFtrack3.GetDistanceFromVertexXY(KFSubMother);
    fDcaSecDaughter3XYKF = (Float_t)KFtrack4.GetDistanceFromVertexXY(KFMother);
    fDcaSecDaughter4XYKF = (Float_t)KFSubMother.GetDistanceFromVertexXY(KFMother);
    fDcaSecDaughterZKF   = (Float_t)KFtrack1.GetDistanceFromVertex(KFSubMother);
    fDcaSecDaughter1ZKF  = (Float_t)KFtrack2.GetDistanceFromVertex(KFSubMother);
    fDcaSecDaughter2ZKF  = (Float_t)KFtrack3.GetDistanceFromVertex(KFSubMother);
    fDcaSecDaughter3ZKF  = (Float_t)KFtrack4.GetDistanceFromVertex(KFMother);
    fDcaSecDaughter4ZKF  = (Float_t)KFSubMother.GetDistanceFromVertex(KFMother);
    // _________________________________________________ //
    AliAnalysisTaskDoubleHypNucTreeLS::SetDaughterInformation("4LLH", kDecayChannel);
    // _________________________________________________ //
    TVector3 vecN(0., 0., 0.);
    TVector3 vecP(0., 0., 0.);
    TVector3 vecM(0., 0., 0.);

    if (ksign == 1) {
      vecP.SetXYZ(KFSubMother.Px(), KFSubMother.Py(), KFSubMother.Pz());
      vecN.SetXYZ(KFtrack4.Px(), KFtrack4.Py(), KFtrack4.Pz());
      vecM.SetXYZ(KFMother.Px(), KFMother.Py(), KFMother.Pz());
    }
    if (ksign == -1) {
      vecN.SetXYZ(KFSubMother.Px(), KFSubMother.Py(), KFSubMother.Pz());
      vecP.SetXYZ(KFtrack4.Px(), KFtrack4.Py(), KFtrack4.Pz());
      vecM.SetXYZ(KFMother.Px(), KFMother.Py(), KFMother.Pz());
    }

    fthetaP   = TMath::ACos((vecP * vecM) / (vecP.Mag() * vecM.Mag()));
    fthetaN   = TMath::ACos((vecN * vecM) / (vecN.Mag() * vecM.Mag()));
    farmalpha = ((vecP.Mag()) * TMath::Cos(fthetaP) - (vecN.Mag()) * TMath::Cos(fthetaN)) / ((vecP.Mag()) * TMath::Cos(fthetaP) + (vecN.Mag()) * TMath::Cos(fthetaN));
    farmpt    = vecP.Mag() * sin(fthetaP);
    return 1;
  }
  // _________________________________________________ //
  else if (Mother == "4LLH" && kDecayChannel == 2) {
    // _________________________________________________ //
    exTrack->Reset(); exTrack1->Reset(); exTrack2->Reset(); exTrack3->Reset(); exTrack4->Reset();
    exTrack1->CopyFromVTrack(track1);
    exTrack2->CopyFromVTrack(track2);
    exTrack3->CopyFromVTrack(track3);
    exTrack4->CopyFromVTrack(track4);
    // __ initialize KFParticles __ //
    KFtrack1 = AliAnalysisTaskDoubleHypNucTreeLS::CreateKFParticle(*exTrack1, AliPID::ParticleMass(AliPID::kHe3), 2);
    KFtrack2 = AliAnalysisTaskDoubleHypNucTreeLS::CreateKFParticle(*exTrack2, AliPID::ParticleMass(AliPID::kProton), 1);
    KFtrack3 = AliAnalysisTaskDoubleHypNucTreeLS::CreateKFParticle(*exTrack3, AliPID::ParticleMass(AliPID::kPion), 1);
    KFtrack4 = AliAnalysisTaskDoubleHypNucTreeLS::CreateKFParticle(*exTrack4, AliPID::ParticleMass(AliPID::kPion), 1);
    // __ construct mother particle __ //
    checkval = 1;
    KFParticleHypNuc KFChecktrack1;
    if (!KFChecktrack1.CheckDaughter(KFtrack1)) checkval = 0;
    KFChecktrack1.AddDaughter(KFtrack3);
    KFParticleHypNuc KFChecktrack4;
    if (!KFChecktrack4.CheckDaughter(KFtrack3)) checkval = 0;
    KFChecktrack4.AddDaughter(KFtrack1);

    KFParticleHypNuc KFSubMother;
    if (checkval) {
      KFSubMother.SetConstructMethod(2);
      KFSubMother.AddDaughter(KFtrack1);
      KFSubMother.AddDaughter(KFtrack3);
      fRecoMethod = 2;
    }
    else {
      KFSubMother.AddDaughter(KFtrack1);
      KFSubMother.AddDaughter(KFtrack3);
      fRecoMethod = 0;
    }

    checkval = 1;
    KFParticleHypNuc KFChecktrack5;
    if (!KFChecktrack5.CheckDaughter(KFSubMother)) checkval = 0;
    KFChecktrack5.AddDaughter(KFtrack2);
    KFChecktrack5.AddDaughter(KFtrack4);
    KFParticleHypNuc KFChecktrack3;
    KFChecktrack3.AddDaughter(KFSubMother);
    KFChecktrack3.AddDaughter(KFtrack2);
    if (!KFChecktrack3.CheckDaughter(KFtrack4)) checkval = 0;
    KFParticleHypNuc KFChecktrack2;
    KFChecktrack2.AddDaughter(KFSubMother);
    if (!KFChecktrack2.CheckDaughter(KFtrack2)) checkval = 0;
    KFChecktrack2.AddDaughter(KFtrack4);

    KFParticleHypNuc KFMother;
    if (checkval) {
      KFMother.SetConstructMethod(2);
      KFMother.AddDaughter(KFSubMother);
      KFMother.AddDaughter(KFtrack2);
      KFMother.AddDaughter(KFtrack4);
      fRecoMethod = 2;
    }
    else {
      KFMother.AddDaughter(KFSubMother);
      KFMother.AddDaughter(KFtrack2);
      KFMother.AddDaughter(KFtrack4);
      fRecoMethod = 0;
    }
    // _________________________________________________ //
    // __ transport mother particle to its decay vertex __ //
    KFSubMother.TransportToDecayVertex();
    KFMother.TransportToDecayVertex();
    // _________________________________________________ //
    // __ get decay vertex position __ //
    SecVertexKF[0]     = KFMother.GetX();
    SecVertexKF[1]     = KFMother.GetY();
    SecVertexKF[2]     = KFMother.GetZ();
    SecVertexErrKF[0]  = KFMother.GetErrX();
    SecVertexErrKF[1]  = KFMother.GetErrY();
    SecVertexErrKF[2]  = KFMother.GetErrZ();
    // __ coord of sec vtx __ //
    fSecVertexXKF      = SecVertexKF[0];
    fSecVertexYKF      = SecVertexKF[1];
    fSecVertexZKF      = SecVertexKF[2];
    fSecVertexXErrKF   = SecVertexErrKF[0];
    fSecVertexYErrKF   = SecVertexErrKF[1];
    fSecVertexZErrKF   = SecVertexErrKF[2];
    // __ vertex quality KF __ //
    fSecVertChi2KF     = KFMother.Chi2();
    fSecVertNDFKF      = KFMother.GetNDF();
    // __ get decay vertex position __ //
    TertVertexKF[0]    = KFSubMother.GetX();
    TertVertexKF[1]    = KFSubMother.GetY();
    TertVertexKF[2]    = KFSubMother.GetZ();
    TertVertexErrKF[0] = KFSubMother.GetErrX();
    TertVertexErrKF[1] = KFSubMother.GetErrY();
    TertVertexErrKF[2] = KFSubMother.GetErrZ();
    // __ coord of tert vtx __ //
    fTertVertexXKF     = TertVertexKF[0];
    fTertVertexYKF     = TertVertexKF[1];
    fTertVertexZKF     = TertVertexKF[2];
    fTertVertexXErrKF  = TertVertexErrKF[0];
    fTertVertexYErrKF  = TertVertexErrKF[1];
    fTertVertexZErrKF  = TertVertexErrKF[2];
    // __ vertex quality KF __ //
    fTertVertChi2KF    = KFSubMother.Chi2();
    fTertVertNDFKF     = KFSubMother.GetNDF();
    // __ create SubMother AliExternalTrackParam __ //
    Double_t xyzV[3];    for (int indx = 0; indx < 3; indx++)  xyzV[indx]    = KFSubMother.GetParameter(indx);
    Double_t pxpypzV[3]; for (int indx = 3; indx < 6; indx++)  pxpypzV[indx] = KFSubMother.GetParameter(indx);
    Double_t covmV[21];  for (int indx = 0; indx < 21; indx++) covmV[indx]   = KFSubMother.GetCovariance(indx);
    exTrack = new AliExternalTrackParam(xyzV, pxpypzV, covmV, ksign);
    exTrack->GetCovarianceXYZPxPyPz(fCovMatrixTrack4);
    fTrackPar4[0]     = exTrack->GetX(); fTrackPar4[1] = exTrack->GetY(); fTrackPar4[2] = exTrack->GetZ(); fTrackPar4[3] = exTrack->GetSnp();
    fTrackPar4[4]     = exTrack->GetTgl(); fTrackPar4[5] = exTrack->GetSigned1Pt(); fTrackPar4[6] = exTrack->GetAlpha();
    // _________________________________________________ //
    particle1->SetXYZM(0., 0., 0., 0.);
    particle2->SetXYZM(0., 0., 0., 0.);
    particle3->SetXYZM(0., 0., 0., 0.);
    particle4->SetXYZM(0., 0., 0., 0.);
    sublorentzsum->SetXYZM(0., 0., 0., 0.);
    lorentzsum->SetXYZM(0., 0., 0., 0.);
    KFtertVertex      = new AliESDVertex(TertVertexKF, TertVertexErrKF);
    fPropDCADaughter  = (exTrack1->PropagateToDCA(KFtertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter    = impar[0];
    fImParzDaughter   = impar[1];
    fDcaSecDaughter   = TMath::Abs(exTrack1->GetD(TertVertexKF[0], TertVertexKF[1], kMagF));
    fPropDCADaughter2 = (exTrack3->PropagateToDCA(KFtertVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter2   = impar[0];
    fImParzDaughter2  = impar[1];
    fDcaSecDaughter2  = TMath::Abs(exTrack3->GetD(TertVertexKF[0], TertVertexKF[1], kMagF));
    if (KFtertVertex) delete KFtertVertex;
    KFsecVertex       = new AliESDVertex(SecVertexKF, SecVertexErrKF);
    fPropDCADaughter4 = (exTrack->PropagateToDCA(KFsecVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter4   = impar[0];
    fImParzDaughter4  = impar[1];
    fDcaSecDaughter4  = TMath::Abs(exTrack->GetD(SecVertexKF[0], SecVertexKF[1], kMagF));
    fPropDCADaughter1 = (exTrack2->PropagateToDCA(KFsecVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter1   = impar[0];
    fImParzDaughter1  = impar[1];
    fDcaSecDaughter1  = TMath::Abs(exTrack2->GetD(SecVertexKF[0], SecVertexKF[1], kMagF));
    fPropDCADaughter3 = (exTrack4->PropagateToDCA(KFsecVertex, kMagF, 10, impar) == 1 ? 1 : 0);
    fImParDaughter3   = impar[0];
    fImParzDaughter3  = impar[1];
    fDcaSecDaughter3  = TMath::Abs(exTrack4->GetD(SecVertexKF[0], SecVertexKF[1], kMagF));
    if (KFsecVertex) delete KFsecVertex;
    particle1->SetXYZM(2. * exTrack1->Px(), 2. * exTrack1->Py(), 2. * exTrack1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
    particle2->SetXYZM(exTrack2->Px(), exTrack2->Py(), exTrack2->Pz(), AliPID::ParticleMass(AliPID::kProton));
    particle3->SetXYZM(exTrack3->Px(), exTrack3->Py(), exTrack3->Pz(), AliPID::ParticleMass(AliPID::kPion));
    particle4->SetXYZM(exTrack4->Px(), exTrack4->Py(), exTrack4->Pz(), AliPID::ParticleMass(AliPID::kPion));
    *sublorentzsum   = *particle1 + *particle3;
    *lorentzsum      = *sublorentzsum + *particle4 + *particle2;
    fmSubMother      = sublorentzsum->M();
    fmMother         = lorentzsum->M();
    if (lorentzsum->M() > kMax4LLHMass || lorentzsum->M() < kMin4LLHMass) return 0;
    // _________________________________________________ //
    // __ mother hypernucleus cos(PA) __ //
    dd[0] = PrimVertexKF[0] - SecVertexKF[0];
    dd[1] = PrimVertexKF[1] - SecVertexKF[1];
    dd[2] = PrimVertexKF[2] - SecVertexKF[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    fPAKF = TMath::Cos(lorentzsum->Angle(*h));
    //if (fPAKF < kPointingAngleCut) return 0;

    DecLength = TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2));
    fctMotherKF = (lorentzsum->M() * DecLength) / lorentzsum->P();
    DecLengthErr = TMath::Sqrt((fPrimVertexXErrKF * fPrimVertexXErrKF) / (fPrimVertexXKF * fPrimVertexXKF)
			       + (fPrimVertexYErrKF * fPrimVertexYErrKF) / (fPrimVertexYKF * fPrimVertexYKF)
			       + (fPrimVertexZErrKF * fPrimVertexZErrKF) / (fPrimVertexZKF * fPrimVertexZKF)
			       + (fSecVertexXErrKF * fSecVertexXErrKF) / (fSecVertexXKF * fSecVertexXKF)
			       + (fSecVertexYErrKF * fSecVertexYErrKF) / (fSecVertexYKF * fSecVertexYKF)
			       + (fSecVertexZErrKF * fSecVertexZErrKF) / (fSecVertexZKF * fSecVertexZKF)
			       + (KFMother.GetErrP() * KFMother.GetErrP()) / (KFMother.GetP() * KFMother.GetP())
			       + (KFMother.GetErrMass() * KFMother.GetErrMass()) / (KFMother.GetMass() * KFMother.GetMass())) * fctMotherKF;
    fctMotherErrKF = DecLengthErr;
    // _________________________________________________ //
    // __ daughter hypernucleus cos(PA) to sec vtx __ //
    dd[0] = SecVertexKF[0] - TertVertexKF[0];
    dd[1] = SecVertexKF[1] - TertVertexKF[1];
    dd[2] = SecVertexKF[2] - TertVertexKF[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    fSubPAKF = TMath::Cos(sublorentzsum->Angle(*h));
    if (fSubPAKF < kPointingAngleCut) return 0;

    DecLength = TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2));
    fctSubMotherKF = (sublorentzsum->M() * DecLength) / sublorentzsum->P();
    DecLengthErr = TMath::Sqrt((fSecVertexXErrKF * fSecVertexXErrKF) / (fSecVertexXKF * fSecVertexXKF)
			       + (fSecVertexYErrKF * fSecVertexYErrKF) / (fSecVertexYKF * fSecVertexYKF)
			       + (fSecVertexZErrKF * fSecVertexZErrKF) / (fSecVertexZKF * fSecVertexZKF)
			       + (fTertVertexXErrKF * fTertVertexXErrKF) / (fTertVertexXKF * fTertVertexXKF)
			       + (fTertVertexYErrKF * fTertVertexYErrKF) / (fTertVertexYKF * fTertVertexYKF)
			       + (fTertVertexZErrKF * fTertVertexZErrKF) / (fTertVertexZKF * fTertVertexZKF)
			       + (KFSubMother.GetErrP() * KFSubMother.GetErrP()) / (KFSubMother.GetP() * KFSubMother.GetP())
			       + (KFSubMother.GetErrMass() * KFSubMother.GetErrMass()) / (KFSubMother.GetMass() * KFSubMother.GetMass())) * fctSubMotherKF;
    fctSubMotherErrKF = DecLengthErr;
    // _________________________________________________ //
    // __ KF daughter hypernucleus cos(PA) to prim vtx __ //
    dd[0] = PrimVertexKF[0] - TertVertexKF[0];
    dd[1] = PrimVertexKF[1] - TertVertexKF[1];
    dd[2] = PrimVertexKF[2] - TertVertexKF[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    fSubPA2KF = TMath::Cos(sublorentzsum->Angle(*h));
    if (fSubPA2KF < kPointingAngleCut) return 0;
    // _________________________________________________ //
    // __ dca tracks KF __ //
    fDCA2BXYKF  = (Float_t)KFtrack1.GetDistanceFromParticleXY(KFtrack3);
    fDCA3B1XYKF = (Float_t)KFSubMother.GetDistanceFromParticleXY(KFtrack2);
    fDCA3B2XYKF = (Float_t)KFSubMother.GetDistanceFromParticleXY(KFtrack4);
    fDCA3B3XYKF = (Float_t)KFtrack2.GetDistanceFromParticleXY(KFtrack4);
    fDCA2BZKF   = (Float_t)KFtrack1.GetDistanceFromParticle(KFtrack3);
    fDCA3B1ZKF  = (Float_t)KFSubMother.GetDistanceFromParticle(KFtrack2);
    fDCA3B2ZKF  = (Float_t)KFSubMother.GetDistanceFromParticle(KFtrack4);
    fDCA3B3ZKF  = (Float_t)KFtrack2.GetDistanceFromParticle(KFtrack4);
    if (TMath::Abs(fDCA3B1XYKF) > kDCATracksCut || TMath::Abs(fDCA3B2XYKF) > kDCATracksCut
	|| TMath::Abs(fDCA3B3XYKF) > kDCATracksCut || TMath::Abs(fDCA2BXYKF) > kDCATracksCut) return 0;
    // _________________________________________________ //
    // __ mother hypernucleus KF information __ //
    fChargeMother     = ksign * 1;
    fDecayChannel     = kDecayChannel;
    fPDGMother        = ksign * fgkPdgCode[kPDGDoubleHyperHydrogen4];
    fmMotherKF        = KFMother.GetMass();
    fmMotherErrKF     = KFMother.GetErrMass();
    fEMotherKF        = KFMother.GetE();
    fEMotherErrKF     = KFMother.GetErrE();
    fpxMotherKF       = KFMother.GetPx();
    fpxMotherErrKF    = KFMother.GetErrPx();
    fpyMotherKF       = KFMother.GetPy();
    fpyMotherErrKF    = KFMother.GetErrPy();
    fpzMotherKF       = KFMother.GetPz();
    fpzMotherErrKF    = KFMother.GetErrPz();
    fptMotherKF       = KFMother.GetPt();
    fptMotherErrKF    = KFMother.GetErrPt();
    fpMotherKF        = KFMother.GetP();
    fpMotherErrKF     = KFMother.GetErrP();
    fyMotherKF        = KFSubMother.GetRapidity();
    // __ daughter hypernucleus KF information __ //
    fmSubMotherKF     = KFSubMother.GetMass();
    fmSubMotherErrKF  = KFSubMother.GetErrMass();
    fESubMotherKF     = KFSubMother.GetE();
    fESubMotherErrKF  = KFSubMother.GetErrE();
    fpxSubMotherKF    = KFSubMother.GetPx();
    fpxSubMotherErrKF = KFSubMother.GetErrPx();
    fpySubMotherKF    = KFSubMother.GetPy();
    fpySubMotherErrKF = KFSubMother.GetErrPy();
    fpzSubMotherKF    = KFSubMother.GetPz();
    fpzSubMotherErrKF = KFSubMother.GetErrPz();
    fptSubMotherKF    = KFSubMother.GetPt();
    fptSubMotherErrKF = KFSubMother.GetErrPt();
    fpSubMotherKF     = KFSubMother.GetP();
    fpSubMotherErrKF  = KFSubMother.GetErrP();
    fySubMotherKF     = KFSubMother.GetRapidity();
    //
    fpxDaughterKF     = KFtrack1.GetPx();
    fpyDaughterKF     = KFtrack1.GetPy();
    fpzDaughterKF     = KFtrack1.GetPz();
    fptDaughterKF     = KFtrack1.GetPt();
    fEDaughterKF      = KFtrack1.GetE();
    fyDaughterKF      = KFtrack1.GetRapidity();
    fpxDaughter1KF    = KFtrack2.GetPx();
    fpyDaughter1KF    = KFtrack2.GetPy();
    fpzDaughter1KF    = KFtrack2.GetPz();
    fptDaughter1KF    = KFtrack2.GetPt();
    fEDaughter1KF     = KFtrack2.GetE();
    fyDaughter1KF     = KFtrack2.GetRapidity();
    fpxDaughter2KF    = KFtrack3.GetPx();
    fpyDaughter2KF    = KFtrack3.GetPy();
    fpzDaughter2KF    = KFtrack3.GetPz();
    fptDaughter2KF    = KFtrack3.GetPt();
    fEDaughter2KF     = KFtrack3.GetE();
    fyDaughter2KF     = KFtrack3.GetRapidity();
    fpxDaughter3KF    = KFtrack4.GetPx();
    fpyDaughter3KF    = KFtrack4.GetPy();
    fpzDaughter3KF    = KFtrack4.GetPz();
    fptDaughter3KF    = KFtrack4.GetPt();
    fEDaughter3KF     = KFtrack4.GetE();
    fyDaughter3KF     = KFtrack4.GetRapidity();
    // daughter momenta
    fpxDaughter       = particle1->Px();
    fpyDaughter       = particle1->Py();
    fpzDaughter       = particle1->Pz();
    fptDaughter       = particle1->Pt();
    fEDaughter        = particle1->E();
    fyDaughter        = particle1->Rapidity();
    fpxDaughter1      = particle2->Px();
    fpyDaughter1      = particle2->Py();
    fpzDaughter1      = particle2->Pz();
    fptDaughter1      = particle2->Pt();
    fEDaughter1       = particle2->E();
    fyDaughter1       = particle2->Rapidity();
    fpxDaughter2      = particle3->Px();
    fpyDaughter2      = particle3->Py();
    fpzDaughter2      = particle3->Pz();
    fptDaughter2      = particle3->Pt();
    fEDaughter2       = particle3->E();
    fyDaughter2       = particle3->Rapidity();
    fpxDaughter3      = particle4->Px();
    fpyDaughter3      = particle4->Py();
    fpzDaughter3      = particle4->Pz();
    fptDaughter3      = particle4->Pt();
    fEDaughter3       = particle4->E();
    fyDaughter3       = particle4->Rapidity();
    // __ Check for V0s __ //
    //AliAnalysisTaskDoubleHypNucTreeLS::GetV0Status();
    // __ dca prim vtx KF __ //
    fDcaDaughterXYKF  = (Float_t)KFtrack1.GetDistanceFromVertexXY(primKFVertex);
    fDcaDaughter1XYKF = (Float_t)KFtrack2.GetDistanceFromVertexXY(primKFVertex);
    fDcaDaughter2XYKF = (Float_t)KFtrack3.GetDistanceFromVertexXY(primKFVertex);
    fDcaDaughter3XYKF = (Float_t)KFtrack4.GetDistanceFromVertexXY(primKFVertex);
    fDcaDaughter4XYKF = (Float_t)KFSubMother.GetDistanceFromVertexXY(primKFVertex);
    fDcaDaughterZKF   = (Float_t)KFtrack1.GetDistanceFromVertex(primKFVertex);
    fDcaDaughter1ZKF  = (Float_t)KFtrack2.GetDistanceFromVertex(primKFVertex);
    fDcaDaughter2ZKF  = (Float_t)KFtrack3.GetDistanceFromVertex(primKFVertex);
    fDcaDaughter3ZKF  = (Float_t)KFtrack4.GetDistanceFromVertex(primKFVertex);
    fDcaDaughter4ZKF  = (Float_t)KFSubMother.GetDistanceFromVertex(primKFVertex);
    // __ dca sec vtx KF __ //
    fDcaSecDaughterXYKF  = (Float_t)KFtrack1.GetDistanceFromVertexXY(KFSubMother);
    fDcaSecDaughter1XYKF = (Float_t)KFtrack2.GetDistanceFromVertexXY(KFMother);
    fDcaSecDaughter2XYKF = (Float_t)KFtrack3.GetDistanceFromVertexXY(KFSubMother);
    fDcaSecDaughter3XYKF = (Float_t)KFtrack4.GetDistanceFromVertexXY(KFMother);
    fDcaSecDaughter4XYKF = (Float_t)KFSubMother.GetDistanceFromVertexXY(KFMother);
    fDcaSecDaughterZKF   = (Float_t)KFtrack1.GetDistanceFromVertex(KFSubMother);
    fDcaSecDaughter1ZKF  = (Float_t)KFtrack2.GetDistanceFromVertex(KFMother);
    fDcaSecDaughter2ZKF  = (Float_t)KFtrack3.GetDistanceFromVertex(KFSubMother);
    fDcaSecDaughter3ZKF  = (Float_t)KFtrack4.GetDistanceFromVertex(KFMother);
    fDcaSecDaughter4ZKF  = (Float_t)KFSubMother.GetDistanceFromVertex(KFMother);
    // _________________________________________________ //
    AliAnalysisTaskDoubleHypNucTreeLS::SetDaughterInformation("4LLH", kDecayChannel);
    // _________________________________________________ //
    TVector3 vecN(0., 0., 0.);
    TVector3 vecP(0., 0., 0.);
    TVector3 vecM(0., 0., 0.);

    if (ksign == 1) {
      vecP.SetXYZ(KFtrack1.Px(), KFtrack1.Py(), KFtrack1.Pz());
      vecN.SetXYZ(KFtrack3.Px(), KFtrack3.Py(), KFtrack3.Pz());
      vecM.SetXYZ(KFSubMother.Px(), KFSubMother.Py(), KFSubMother.Pz());
    }
    if (ksign == -1) {
      vecN.SetXYZ(KFtrack1.Px(), KFtrack1.Py(), KFtrack1.Pz());
      vecP.SetXYZ(KFtrack3.Px(), KFtrack3.Py(), KFtrack3.Pz());
      vecM.SetXYZ(KFSubMother.Px(), KFSubMother.Py(), KFSubMother.Pz());
    }

    fthetaP   = TMath::ACos((vecP * vecM) / (vecP.Mag() * vecM.Mag()));
    fthetaN   = TMath::ACos((vecN * vecM) / (vecN.Mag() * vecM.Mag()));
    farmalpha = ((vecP.Mag()) * TMath::Cos(fthetaP) - (vecN.Mag()) * TMath::Cos(fthetaN)) / ((vecP.Mag()) * TMath::Cos(fthetaP) + (vecN.Mag()) * TMath::Cos(fthetaN));
    farmpt    = vecP.Mag() * sin(fthetaP);
    return 1;
  }
  else return 0;
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTreeLS::GetV0Status() {
  fV0VertexX_13 = -9999;
  fV0VertexY_13 = -9999;
  fV0VertexZ_13 = -9999;
  fV0VertexX_14 = -9999;
  fV0VertexY_14 = -9999;
  fV0VertexZ_14 = -9999;
  fV0VertexX_23 = -9999;
  fV0VertexY_23 = -9999;
  fV0VertexZ_23 = -9999;
  fV0VertexX_24 = -9999;
  fV0VertexY_24 = -9999;
  fV0VertexZ_24 = -9999;
  fisOnlineV0_13 = 0;
  fisOnlineV0_14 = 0;
  fisOnlineV0_23 = 0;
  fisOnlineV0_24 = 0;
  for (Int_t ivertex = 0; ivertex < fAODevent->GetNumberOfV0s(); ivertex++) {

    //Get V0
    AliAODv0 *fV0 = fAODevent->GetV0(ivertex);

    //check if on the fly
    if(!fV0->GetOnFlyStatus()) continue;

    //ChargeCorrection
    Bool_t v0ChargeCorrect = kTRUE;        
    AliAODTrack* trackN = (AliAODTrack *)(fV0->GetSecondaryVtx()->GetDaughter(0));//pos Track
    AliAODTrack* trackP = (AliAODTrack *)(fV0->GetSecondaryVtx()->GetDaughter(1));//neg Track
    if (trackN->GetSign() > 0 ) {
      trackN = (AliAODTrack *)(fV0->GetSecondaryVtx()->GetDaughter(1));//neg Track
      trackP = (AliAODTrack *)(fV0->GetSecondaryVtx()->GetDaughter(0));//pos Track
      v0ChargeCorrect = kFALSE;
    }
    if (track1->GetSign() > 0) {
      if (trackP->GetLabel() == track1->GetLabel() && trackN->GetLabel() == track3->GetLabel()) {
	fV0VertexX_13 = fV0->Xv();
	fV0VertexY_13 = fV0->Yv();
	fV0VertexZ_13 = fV0->Zv();
	fisOnlineV0_13 = 1;
      }
      if (trackP->GetLabel() == track1->GetLabel() && trackN->GetLabel() == track4->GetLabel()) {
	fV0VertexX_14 = fV0->Xv();
	fV0VertexY_14 = fV0->Yv();
	fV0VertexZ_14 = fV0->Zv();
	fisOnlineV0_14 = 1;
      }
      if (trackP->GetLabel() == track2->GetLabel() && trackN->GetLabel() == track3->GetLabel()) {
	fV0VertexX_23 = fV0->Xv();
	fV0VertexY_23 = fV0->Yv();
	fV0VertexZ_23 = fV0->Zv();
	fisOnlineV0_23 = 1;
      }
      if (trackP->GetLabel() == track2->GetLabel() && trackN->GetLabel() == track4->GetLabel()) {
	fV0VertexX_24 = fV0->Xv();
	fV0VertexY_24 = fV0->Yv();
	fV0VertexZ_24 = fV0->Zv();
	fisOnlineV0_24 = 1;
      }
    }
    if (track1->GetSign() < 0) {
      if (trackN->GetLabel() == track1->GetLabel() && trackP->GetLabel() == track3->GetLabel()) {
	fV0VertexX_13 = fV0->Xv();
	fV0VertexY_13 = fV0->Yv();
	fV0VertexZ_13 = fV0->Zv();
	fisOnlineV0_13 = 1;
      }
      if (trackN->GetLabel() == track1->GetLabel() && trackP->GetLabel() == track4->GetLabel()) {
	fV0VertexX_14 = fV0->Xv();
	fV0VertexY_14 = fV0->Yv();
	fV0VertexZ_14 = fV0->Zv();
	fisOnlineV0_14 = 1;
      }
      if (trackN->GetLabel() == track2->GetLabel() && trackP->GetLabel() == track3->GetLabel()) {
	fV0VertexX_23 = fV0->Xv();
	fV0VertexY_23 = fV0->Yv();
	fV0VertexZ_23 = fV0->Zv();
	fisOnlineV0_23 = 1;
      }
      if (trackN->GetLabel() == track2->GetLabel() && trackP->GetLabel() == track4->GetLabel()) {
	fV0VertexX_24 = fV0->Xv();
	fV0VertexY_24 = fV0->Yv();
	fV0VertexZ_24 = fV0->Zv();
	fisOnlineV0_24 = 1;
      }
    }
  }
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTreeLS::SetDaughterInformation(TString Mother, Int_t kDecayChannel) {

  Float_t xv[2], yv[3];
  AliMCParticle* part;
  // _________________________________ //
  // __ He3 track __ //
  fpDaughter = track1->GetTPCmomentum();
  if (fMCtrue) {
    part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(track1->GetLabel()))->Particle());
    if (Mother != "4LH") fLabelDaughterKF = TMath::Abs(part->PdgCode()) == fgkPdgCode[kPDGHelium3];
    if (Mother == "4LH") fLabelDaughterKF = TMath::Abs(part->PdgCode()) == fgkPdgCode[kPDGHelium4];
  }
  fpxDaughterUnProp = track1->Px();
  fpyDaughterUnProp = track1->Py();
  fpzDaughterUnProp = track1->Pz();
  fptDaughterUnProp = track1->Pt();
  fTrackPIDDaughter = track1->GetPIDForTracking();
  fdEdxDaughter     = track1->GetTPCsignal();
  if (Mother != "4LH") {
    fdEdxSigmaDaughter = AliAnalysisTaskDoubleHypNucTreeLS::Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
    if (fMCtrue) fdEdxSigmaDaughter = fPID->NumberOfSigmasTPC(track1, AliPID::kHe3);
    fdEdxSigmaTriton   = AliAnalysisTaskDoubleHypNucTreeLS::Bethe(*track1, AliPID::ParticleMass(AliPID::kTriton), 1, fBetheParamsT);
    fdEdxSigmaAlpha    = AliAnalysisTaskDoubleHypNucTreeLS::Bethe(*track1, AliPID::ParticleMass(AliPID::kAlpha), 2, fBetheParamsHe);
  }
  else {
    fdEdxSigmaDaughter = AliAnalysisTaskDoubleHypNucTreeLS::Bethe(*track1, AliPID::ParticleMass(AliPID::kAlpha), 2, fBetheParamsHe);
    if (fMCtrue) fdEdxSigmaDaughter = fPID->NumberOfSigmasTPC(track1, AliPID::kAlpha);
    fdEdxSigmaTriton   = AliAnalysisTaskDoubleHypNucTreeLS::Bethe(*track1, AliPID::ParticleMass(AliPID::kTriton), 1, fBetheParamsT);
    fdEdxSigmaAlpha    = AliAnalysisTaskDoubleHypNucTreeLS::Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
  }
  track1->GetImpactParameters(xv, yv);
  fDcaDaughter       = xv[0];
  fDcazDaughter      = xv[1];
  fSigmaYXDaughter   = yv[0];
  fSigmaXYZDaughter  = yv[1];
  fSigmaZDaughter    = yv[2];
  fNclsDaughter      = track1->GetTPCNcls();
  fNclsITSDaughter   = track1->GetITSNcls();
  fChi2Daughter      = track1->GetTPCchi2() / (Float_t)track1->GetTPCNcls();
  fEtaDaughter       = track1->Eta();
  fPhiDaughter       = track1->Phi();
  fGeoLengthDaughter = AliAnalysisTaskDoubleHypNucTreeLS::GeoLength(*track1);
  fTOFSignalDaughter = AliAnalysisTaskDoubleHypNucTreeLS::GetTOFSignal(*track1);
  fPtUncertDaughter  = TMath::Sqrt(exTrack1->GetSigma1Pt2());
  fTPCRefitDaughter  = (track1->GetStatus() & AliAODTrack::kTPCrefit) != 0;
  fITSRefitDaughter  = (track1->GetStatus() & AliAODTrack::kITSrefit) != 0;
  fITSLayer1Daughter = track1->HasPointOnITSLayer(0);
  fITSLayer2Daughter = track1->HasPointOnITSLayer(1);
  fITSLayer3Daughter = track1->HasPointOnITSLayer(2);
  fITSLayer4Daughter = track1->HasPointOnITSLayer(3);
  fITSLayer5Daughter = track1->HasPointOnITSLayer(4);
  fITSLayer6Daughter = track1->HasPointOnITSLayer(5);
  track1->GetCovarianceXYZPxPyPz(fCovMatrixTrack);
  fTrackPar[0] = track1->GetX(); fTrackPar[1] = track1->GetY(); fTrackPar[2] = track1->GetZ(); fTrackPar[3] = track1->GetSnp();
  fTrackPar[4] = track1->GetTgl(); fTrackPar[5] = track1->GetSigned1Pt(); fTrackPar[6] = track1->GetAlpha();
  // _________________________________ //
  // __ p track __ //
  if (Mother != "3LH" && Mother != "4LH") {      
    fpDaughter1 = track2->GetTPCmomentum();
    if (fMCtrue) {
      part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(track2->GetLabel()))->Particle());
      fLabelDaughter1KF = TMath::Abs(part->PdgCode()) == fgkPdgCode[kPDGProton];
    }
    fpxDaughter1UnProp  = track2->Px();
    fpyDaughter1UnProp  = track2->Py();
    fpzDaughter1UnProp  = track2->Pz();
    fptDaughter1UnProp  = track2->Pt();
    fTrackPIDDaughter1  = track2->GetPIDForTracking();
    fdEdxDaughter1      = track2->GetTPCsignal();
    fdEdxSigmaDaughter1 = fPID->NumberOfSigmasTPC(track2, AliPID::kProton);
    track2->GetImpactParameters(xv, yv);
    fDcaDaughter1       = xv[0];
    fDcazDaughter1      = xv[1];
    fSigmaYXDaughter1   = yv[0];
    fSigmaXYZDaughter1  = yv[1];
    fSigmaZDaughter1    = yv[2];
    fNclsDaughter1      = track2->GetTPCNcls();
    fNclsITSDaughter1   = track2->GetITSNcls();
    fChi2Daughter1      = track2->GetTPCchi2() / (Float_t)track2->GetTPCNcls();
    fEtaDaughter1       = track2->Eta();
    fPhiDaughter1       = track2->Phi();
    fGeoLengthDaughter1 = AliAnalysisTaskDoubleHypNucTreeLS::GeoLength(*track2);
    fTOFSignalDaughter1 = AliAnalysisTaskDoubleHypNucTreeLS::GetTOFSignal(*track2);
    fPtUncertDaughter1  = TMath::Sqrt(exTrack2->GetSigma1Pt2());
    fTPCRefitDaughter1  = (track2->GetStatus() & AliAODTrack::kTPCrefit) != 0;
    fITSRefitDaughter1  = (track2->GetStatus() & AliAODTrack::kITSrefit) != 0;
    fdEdxSigmaPion      = fPID->NumberOfSigmasTPC(track2, AliPID::kPion);
    fdEdxSigmaDeuteron  = AliAnalysisTaskDoubleHypNucTreeLS::Bethe(*track2, AliPID::ParticleMass(AliPID::kDeuteron), 1, fBetheParamsT);
    fITSLayer1Daughter1 = track2->HasPointOnITSLayer(0);
    fITSLayer2Daughter1 = track2->HasPointOnITSLayer(1);
    fITSLayer3Daughter1 = track2->HasPointOnITSLayer(2);
    fITSLayer4Daughter1 = track2->HasPointOnITSLayer(3);
    fITSLayer5Daughter1 = track2->HasPointOnITSLayer(4);
    fITSLayer6Daughter1 = track2->HasPointOnITSLayer(5);
    track2->GetCovarianceXYZPxPyPz(fCovMatrixTrack1);
    fTrackPar1[0] = track2->GetX(); fTrackPar1[1] = track2->GetY(); fTrackPar1[2] = track2->GetZ(); fTrackPar1[3] = track2->GetSnp();
    fTrackPar1[4] = track2->GetTgl(); fTrackPar1[5] = track2->GetSigned1Pt(); fTrackPar1[6] = track2->GetAlpha();
  }
  // _________________________________ //
  // __ pion track __ //
  fpDaughter2 = track3->GetTPCmomentum();
  if (fMCtrue) {
    part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(track3->GetLabel()))->Particle());
    fLabelDaughter2KF = TMath::Abs(part->PdgCode()) == fgkPdgCode[kPDGPionPlus];
  }
  fpxDaughter2UnProp  = track3->Px();
  fpyDaughter2UnProp  = track3->Py();
  fpzDaughter2UnProp  = track3->Pz();
  fptDaughter2UnProp  = track3->Pt();
  fTrackPIDDaughter2  = track3->GetPIDForTracking();
  fdEdxDaughter2      = track3->GetTPCsignal();
  fdEdxSigmaDaughter2 = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
  track3->GetImpactParameters(xv, yv);
  fDcaDaughter2       = xv[0];
  fDcazDaughter2      = xv[1];
  fSigmaYXDaughter2   = yv[0];
  fSigmaXYZDaughter2  = yv[1];
  fSigmaZDaughter2    = yv[2];
  fNclsDaughter2      = track3->GetTPCNcls();
  fNclsITSDaughter2   = track3->GetITSNcls();
  fChi2Daughter2      = track3->GetTPCchi2() / (Float_t)track3->GetTPCNcls();
  fEtaDaughter2       = track3->Eta();
  fPhiDaughter2       = track3->Phi();
  fGeoLengthDaughter2 = AliAnalysisTaskDoubleHypNucTreeLS::GeoLength(*track3);
  fTOFSignalDaughter2 = AliAnalysisTaskDoubleHypNucTreeLS::GetTOFSignal(*track3);
  fPtUncertDaughter2  = TMath::Sqrt(exTrack3->GetSigma1Pt2());
  fTPCRefitDaughter2  = (track3->GetStatus() & AliAODTrack::kTPCrefit) != 0;
  fITSRefitDaughter2  = (track3->GetStatus() & AliAODTrack::kITSrefit) != 0;
  fITSLayer1Daughter2 = track3->HasPointOnITSLayer(0);
  fITSLayer2Daughter2 = track3->HasPointOnITSLayer(1);
  fITSLayer3Daughter2 = track3->HasPointOnITSLayer(2);
  fITSLayer4Daughter2 = track3->HasPointOnITSLayer(3);
  fITSLayer5Daughter2 = track3->HasPointOnITSLayer(4);
  fITSLayer6Daughter2 = track3->HasPointOnITSLayer(5);
  track3->GetCovarianceXYZPxPyPz(fCovMatrixTrack2);
  fTrackPar2[0] = track3->GetX(); fTrackPar2[1] = track3->GetY(); fTrackPar2[2] = track3->GetZ(); fTrackPar2[3] = track3->GetSnp();
  fTrackPar2[4] = track3->GetTgl(); fTrackPar2[5] = track3->GetSigned1Pt(); fTrackPar2[6] = track3->GetAlpha();

  if (Mother == "4LLH") {
    // _________________________________ //
    // __ pion track __ //
    fpDaughter3 = track4->GetTPCmomentum();
    if (fMCtrue) {
      part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(track4->GetLabel()))->Particle());
      fLabelDaughter3KF = TMath::Abs(part->PdgCode()) == fgkPdgCode[kPDGPionPlus];
    }
    fpxDaughter3UnProp  = track4->Px();
    fpyDaughter3UnProp  = track4->Py();
    fpzDaughter3UnProp  = track4->Pz();
    fptDaughter3UnProp  = track4->Pt();
    fTrackPIDDaughter3  = track4->GetPIDForTracking();
    fdEdxDaughter3      = track4->GetTPCsignal();
    fdEdxSigmaDaughter3 = fPID->NumberOfSigmasTPC(track4, AliPID::kPion);
    track4->GetImpactParameters(xv, yv);
    fDcaDaughter3       = xv[0];
    fDcazDaughter3      = xv[1];
    fSigmaYXDaughter3   = yv[0];
    fSigmaXYZDaughter3  = yv[1];
    fSigmaZDaughter3    = yv[2];
    fNclsDaughter3      = track4->GetTPCNcls();
    fNclsITSDaughter3   = track4->GetITSNcls();
    fChi2Daughter3      = track4->GetTPCchi2() / (Float_t)track4->GetTPCNcls();
    fEtaDaughter3       = track4->Eta();
    fPhiDaughter3       = track4->Phi();
    fGeoLengthDaughter3 = AliAnalysisTaskDoubleHypNucTreeLS::GeoLength(*track4);
    fTOFSignalDaughter3 = AliAnalysisTaskDoubleHypNucTreeLS::GetTOFSignal(*track4);
    fPtUncertDaughter3  = TMath::Sqrt(exTrack4->GetSigma1Pt2());
    fTPCRefitDaughter3  = (track4->GetStatus() & AliAODTrack::kTPCrefit) != 0;
    fITSRefitDaughter3  = (track4->GetStatus() & AliAODTrack::kITSrefit) != 0;
    fITSLayer1Daughter3 = track4->HasPointOnITSLayer(0);
    fITSLayer2Daughter3 = track4->HasPointOnITSLayer(1);
    fITSLayer3Daughter3 = track4->HasPointOnITSLayer(2);
    fITSLayer4Daughter3 = track4->HasPointOnITSLayer(3);
    fITSLayer5Daughter3 = track4->HasPointOnITSLayer(4);
    fITSLayer6Daughter3 = track4->HasPointOnITSLayer(5);
    track4->GetCovarianceXYZPxPyPz(fCovMatrixTrack3);
    fTrackPar3[0] = track4->GetX(); fTrackPar3[1] = track4->GetY(); fTrackPar3[2] = track4->GetZ(); fTrackPar3[3] = track4->GetSnp();
    fTrackPar3[4] = track4->GetTgl(); fTrackPar3[5] = track4->GetSigned1Pt(); fTrackPar3[6] = track4->GetAlpha();
  }
}
// _________________________________________________ //
KFParticle AliAnalysisTaskDoubleHypNucTreeLS::CreateKFParticle(AliExternalTrackParam& track, Double_t Mass, Int_t Charge) {

  Double_t fP[6];
  track.GetXYZ(fP);
  track.PxPyPz(fP + 3);
  Int_t fQ = track.Charge() * TMath::Abs(Charge);
  fP[3] *= TMath::Abs(Charge);
  fP[4] *= TMath::Abs(Charge);
  fP[5] *= TMath::Abs(Charge);

  Double_t pt = 1. / TMath::Abs(track.GetParameter()[4]) * TMath::Abs(Charge);
  Double_t cs = TMath::Cos(track.GetAlpha()), sn = TMath::Sin(track.GetAlpha());
  Double_t r = TMath::Sqrt((1. - track.GetParameter()[2]) * (1. + track.GetParameter()[2]));

  Double_t m00 = -sn, m10 = cs;
  Double_t m23 = -pt * (sn + track.GetParameter()[2] * cs / r), m43 = -pt * pt * (r * cs - track.GetParameter()[2] * sn);
  Double_t m24 = pt * (cs - track.GetParameter()[2] * sn / r), m44 = -pt * pt * (r * sn + track.GetParameter()[2] * cs);
  Double_t m35 = pt, m45 = -pt * pt * track.GetParameter()[3];

  m43 *= track.GetSign();
  m44 *= track.GetSign();
  m45 *= track.GetSign();

  const Double_t* cTr = track.GetCovariance();
  Double_t fC[21];
  fC[0] = cTr[0] * m00 * m00;
  fC[1] = cTr[0] * m00 * m10;
  fC[2] = cTr[0] * m10 * m10;
  fC[3] = cTr[1] * m00;
  fC[4] = cTr[1] * m10;
  fC[5] = cTr[2];
  fC[6] = m00 * (cTr[3] * m23 + cTr[10] * m43);
  fC[7] = m10 * (cTr[3] * m23 + cTr[10] * m43);
  fC[8] = cTr[4] * m23 + cTr[11] * m43;
  fC[9] = m23 * (cTr[5] * m23 + cTr[12] * m43) + m43 * (cTr[12] * m23 + cTr[14] * m43);
  fC[10] = m00 * (cTr[3] * m24 + cTr[10] * m44);
  fC[11] = m10 * (cTr[3] * m24 + cTr[10] * m44);
  fC[12] = cTr[4] * m24 + cTr[11] * m44;
  fC[13] = m23 * (cTr[5] * m24 + cTr[12] * m44) + m43 * (cTr[12] * m24 + cTr[14] * m44);
  fC[14] = m24 * (cTr[5] * m24 + cTr[12] * m44) + m44 * (cTr[12] * m24 + cTr[14] * m44);
  fC[15] = m00 * (cTr[6] * m35 + cTr[10] * m45);
  fC[16] = m10 * (cTr[6] * m35 + cTr[10] * m45);
  fC[17] = cTr[7] * m35 + cTr[11] * m45;
  fC[18] = m23 * (cTr[8] * m35 + cTr[12] * m45) + m43 * (cTr[13] * m35 + cTr[14] * m45);
  fC[19] = m24 * (cTr[8] * m35 + cTr[12] * m45) + m44 * (cTr[13] * m35 + cTr[14] * m45);
  fC[20] = m35 * (cTr[9] * m35 + cTr[13] * m45) + m45 * (cTr[13] * m35 + cTr[14] * m45);

  KFParticle part;
  part.Create(fP, fC, fQ, Mass);
  return part;
}
// _________________________________________________ //
AliESDVertex* AliAnalysisTaskDoubleHypNucTreeLS::AODToESDVertex(const AliAODVertex &aodVert){  
  Double_t covmatrix[6];  Double_t position[3]; Double_t chi2; Int_t nContributors;
  aodVert.GetCovMatrix(covmatrix); aodVert.GetXYZ(position);
  chi2 = aodVert.GetChi2(); nContributors = aodVert.GetNDF(); nContributors = (nContributors + 3) / 2;
  AliESDVertex *esdVert = new AliESDVertex(position, covmatrix, chi2, nContributors, "primaryVertex");
  return esdVert;
}
// _________________________________________________ //
KFVertex AliAnalysisTaskDoubleHypNucTreeLS::CreateKFVertex(const AliVVertex& vertex) {

  /// GetTrack parameters
  Double_t param[6];
  Double_t cov[6];

  vertex.GetXYZ(param);
  vertex.GetCovarianceMatrix(cov);

  KFPVertex kfpVtx;
  /// Set the values
  Float_t paramF[3] = { (Float_t)param[0],(Float_t)param[1],(Float_t)param[2] };
  kfpVtx.SetXYZ(paramF);
  Float_t covF[6] = { (Float_t)cov[0],(Float_t)cov[1],(Float_t)cov[2],
		      (Float_t)cov[3],(Float_t)cov[4],(Float_t)cov[5] };
  kfpVtx.SetCovarianceMatrix(covF);
  KFVertex KFVtx(kfpVtx);
  return KFVtx;
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTreeLS::Terminate(const Option_t*) {
  if (!GetOutputData(0)) return;
}
/// Set trigger information in reduced event
/// \return returns kTRUE is successful.
// _________________________________________________ //
Bool_t AliAnalysisTaskDoubleHypNucTreeLS::TriggerSelection() {
  //******************************
  //*   get trigger information  *
  //******************************

  MB = 0;
  HMV0 = 0;
  HMSPD = 0;
  HNU = 0;
  HQU = 0;
  Central = 0;
  SemiCentral = 0;

  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7)        MB = kTRUE;
  if (fInputHandler->IsEventSelected() & AliVEvent::kHighMultV0)  HMV0 = kTRUE;
  if (fInputHandler->IsEventSelected() & AliVEvent::kHighMultSPD) HMSPD = kTRUE;
  if (fInputHandler->IsEventSelected() & AliVEvent::kCentral)     Central = kTRUE;
  if (fInputHandler->IsEventSelected() & AliVEvent::kSemiCentral) SemiCentral = kTRUE;

  Int_t nTrdTracks = fAODevent->GetNumberOfTrdTracks();
  // __ Data: get TRD trigger information from trigger classes __ //
  TString classes = fAODevent->GetFiredTriggerClasses();
  if (classes.Contains("HNU")) HNU = 1;
  if (classes.Contains("HQU")) HQU = 1;
  // fill histogram
  fHistTrigger->Fill(0);
  if (MB)          fHistTrigger->Fill(1);
  if (HMV0)        fHistTrigger->Fill(2);
  if (HMSPD)       fHistTrigger->Fill(3);
  if (HNU)         fHistTrigger->Fill(4);
  if (HQU)         fHistTrigger->Fill(5);
  if (Central)     fHistTrigger->Fill(6);
  if (SemiCentral) fHistTrigger->Fill(7);
  if (!fisWOTrackVertex) {
    fHistTrigger1->Fill(0);
    if (MB)          fHistTrigger1->Fill(1);
    if (HMV0)        fHistTrigger1->Fill(2);
    if (HMSPD)       fHistTrigger1->Fill(3);
    if (HNU)         fHistTrigger1->Fill(4);
    if (HQU)         fHistTrigger1->Fill(5);
    if (Central)     fHistTrigger1->Fill(6);
    if (SemiCentral) fHistTrigger1->Fill(7);
  }
  if (fisWOTrackVertex) {
    fHistTrigger2->Fill(0);
    if (MB)          fHistTrigger2->Fill(1);
    if (HMV0)        fHistTrigger2->Fill(2);
    if (HMSPD)       fHistTrigger2->Fill(3);
    if (HNU)         fHistTrigger2->Fill(4);
    if (HQU)         fHistTrigger2->Fill(5);
    if (Central)     fHistTrigger2->Fill(6);
    if (SemiCentral) fHistTrigger2->Fill(7);
  }
  Bool_t isTriggered = kFALSE;
  if (MB || HMV0 || HMSPD || HNU || HQU || Central || SemiCentral) isTriggered = kTRUE;
  return isTriggered;
}
// _________________________________________________ //
// __ calculate stack label of specific daughter __ //
Int_t AliAnalysisTaskDoubleHypNucTreeLS::GetLabel(Int_t labelGrandMother, Int_t particlePdgCode) {
    
  Int_t labelFirstDaughter = mcEvent->GetLabelOfParticleFirstDaughter(TMath::Abs(labelGrandMother));
  Int_t labelLastDaughter = mcEvent->GetLabelOfParticleLastDaughter(TMath::Abs(labelGrandMother));
  Int_t diff = TMath::Abs(labelLastDaughter - labelFirstDaughter) + 1;
  Int_t returnval = -99;
    
  for (Int_t i = 0; i < diff; i++) {
        
    Int_t labelDaughter = labelFirstDaughter + i;
    AliMCParticle* Daughter = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelDaughter))->Particle());
        
    if (Daughter->PdgCode() == particlePdgCode) {
      returnval = labelDaughter;
    }
    if (Daughter) delete Daughter;
  }
    
  return returnval;
}
// _________________________________________________ //
// __ calculate stack label of specific granddaughter __ //
Int_t AliAnalysisTaskDoubleHypNucTreeLS::GetLabel(Int_t labelGrandMother, Int_t particlePdgCode, Int_t motherparticlePdgCode) {
  Int_t labelFirstDaughter = mcEvent->GetLabelOfParticleFirstDaughter(TMath::Abs(labelGrandMother));
  Int_t labelLastDaughter = mcEvent->GetLabelOfParticleLastDaughter(TMath::Abs(labelGrandMother));
  Int_t diff = TMath::Abs(labelLastDaughter - labelFirstDaughter) + 1;
  Int_t returnval = -99;
    
  for (Int_t i = 0; i < diff; i++) {
        
    Int_t labelDaughter = labelFirstDaughter + i;
    AliMCParticle* Daughter = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelDaughter))->Particle());
        
    if (Daughter->PdgCode() == motherparticlePdgCode) {
            
      Int_t labelFirstEnkel = mcEvent->GetLabelOfParticleFirstDaughter(TMath::Abs(labelDaughter));
      Int_t labelLastEnkel = mcEvent->GetLabelOfParticleLastDaughter(TMath::Abs(labelDaughter));
      Int_t diffEnkel = TMath::Abs(labelLastEnkel - labelFirstEnkel) + 1;
            
      for (Int_t j = 0; j < diffEnkel; j++) {
                
	Int_t labelEnkel = labelFirstEnkel + j;
	AliMCParticle* Enkel = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelEnkel))->Particle());
                
	if (Enkel->PdgCode() == particlePdgCode) {
	  returnval = labelEnkel;
	}
	if (Enkel) delete Enkel;
      }
    }
    if (Daughter) delete Daughter;
  }
  return returnval;
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTreeLS::MCGenerated() {
    
  // Monte Carlo for genenerated particles
  stackN = 0;
  for (stackN = 0; stackN < mcEvent->GetNumberOfTracks(); stackN++) //loop over stack
    {
        
      ParticleMother = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(stackN))->Particle());
      PDGCodeMother = ParticleMother->PdgCode();
        
      //DoubleHyperHydrogen4
      if (PDGCodeMother == fgkPdgCode[kPDGDoubleHyperHydrogen4]) //check mother PDG
        {
	  AliAnalysisTaskDoubleHypNucTreeLS::MCFourBodyDecay(stackN, ParticleMother, PDGCodeMother, 1, fgkPdgCode[kPDGHelium3],
							       fgkPdgCode[kPDGProton], fgkPdgCode[kPDGPionMinus], fgkPdgCode[kPDGPionMinus],
							       AliPID::ParticleMass(AliPID::kHe3), AliPID::ParticleMass(AliPID::kProton), AliPID::ParticleMass(AliPID::kPion), AliPID::ParticleMass(AliPID::kPion));
        }
      //AntiDoubleHyperHydrogen4
      if (PDGCodeMother == fgkPdgCode[kPDGAntiDoubleHyperHydrogen4]) //check mother PDG
        {
	  AliAnalysisTaskDoubleHypNucTreeLS::MCFourBodyDecay(stackN, ParticleMother, PDGCodeMother, 1, fgkPdgCode[kPDGAntiHelium3],
							       fgkPdgCode[kPDGAntiProton], fgkPdgCode[kPDGPionPlus], fgkPdgCode[kPDGPionPlus],
							       AliPID::ParticleMass(AliPID::kHe3), AliPID::ParticleMass(AliPID::kProton), AliPID::ParticleMass(AliPID::kPion), AliPID::ParticleMass(AliPID::kPion));
        }
      //HyperHelium4
      if (PDGCodeMother == fgkPdgCode[kPDGHyperHelium4] || PDGCodeMother == fgkPdgCode[kPDGHyperHelium4Star]) //check mother PDG
        {
	  AliAnalysisTaskDoubleHypNucTreeLS::MCThreeBodyDecay(stackN, ParticleMother, PDGCodeMother, fgkPdgCode[kPDGHelium3],
								fgkPdgCode[kPDGProton], fgkPdgCode[kPDGPionMinus],
								AliPID::ParticleMass(AliPID::kHe3), AliPID::ParticleMass(AliPID::kProton), AliPID::ParticleMass(AliPID::kPion));
        }
      //AntiHyperHelium4
      if (PDGCodeMother == fgkPdgCode[kPDGAntiHyperHelium4] || PDGCodeMother == fgkPdgCode[kPDGAntiHyperHelium4Star]) //check mother PDG
        {
            
	  AliAnalysisTaskDoubleHypNucTreeLS::MCThreeBodyDecay(stackN, ParticleMother, PDGCodeMother, fgkPdgCode[kPDGAntiHelium3],
								fgkPdgCode[kPDGAntiProton], fgkPdgCode[kPDGPionPlus],
								AliPID::ParticleMass(AliPID::kHe3), AliPID::ParticleMass(AliPID::kProton), AliPID::ParticleMass(AliPID::kPion));
        }
      if(!kSkip3LHAnalysis){
	//HyperHydrogen3
	if (PDGCodeMother == fgkPdgCode[kPDGHyperHydrogen3]) //check mother PDG
	  {
                
	    AliAnalysisTaskDoubleHypNucTreeLS::MCTwoBodyDecay(stackN, ParticleMother, PDGCodeMother, fgkPdgCode[kPDGHelium3],
                                                                fgkPdgCode[kPDGPionMinus], AliPID::ParticleMass(AliPID::kHe3), AliPID::ParticleMass(AliPID::kPion));
	  }
	//AntiHyperHydrogen3
	if (PDGCodeMother == fgkPdgCode[kPDGAntiHyperHydrogen3]) //check mother PDG
	  {
                
	    AliAnalysisTaskDoubleHypNucTreeLS::MCTwoBodyDecay(stackN, ParticleMother, PDGCodeMother, fgkPdgCode[kPDGAntiHelium3],
                                                                fgkPdgCode[kPDGPionPlus], AliPID::ParticleMass(AliPID::kHe3), AliPID::ParticleMass(AliPID::kPion));
	  }
      }
      if(!kSkip4LHAnalysis){
	//HyperHydrogen4
	if (PDGCodeMother == fgkPdgCode[kPDGHyperHydrogen4] || PDGCodeMother == fgkPdgCode[kPDGHyperHydrogen4Star]) //check mother PDG
	  {
                
	    AliAnalysisTaskDoubleHypNucTreeLS::MCTwoBodyDecay(stackN, ParticleMother, PDGCodeMother, fgkPdgCode[kPDGHelium4],
                                                                fgkPdgCode[kPDGPionMinus], AliPID::ParticleMass(AliPID::kAlpha), AliPID::ParticleMass(AliPID::kPion));
	  }
	//AntiHyperHydrogen4
	if (PDGCodeMother == fgkPdgCode[kPDGAntiHyperHydrogen4] || PDGCodeMother == fgkPdgCode[kPDGAntiHyperHydrogen4Star]) //check mother PDG
	  {
                
	    AliAnalysisTaskDoubleHypNucTreeLS::MCTwoBodyDecay(stackN, ParticleMother, PDGCodeMother, fgkPdgCode[kPDGAntiHelium4],
                                                                fgkPdgCode[kPDGPionPlus], AliPID::ParticleMass(AliPID::kAlpha), AliPID::ParticleMass(AliPID::kPion));
	  }
      }
      if (ParticleMother) delete ParticleMother;
    }//end loop over stack
}
//_________________________________________________
void AliAnalysisTaskDoubleHypNucTreeLS::MCFourBodyDecay(Int_t stackN, AliMCParticle* ParticleMother, Long_t PDGMother, Int_t kDecayChannel,
							  Long_t PDGFirstDaughter, Long_t PDGSecondDaughter, Long_t PDGThirdDaughter, Long_t PDGFourthDaughter,
							  Double_t massFirstDaughter, Double_t massSecondDaughter, Double_t massThirdDaughter, Double_t massFourthDaughter) {
    
  //Check decay channel
  Int_t flabel = mcEvent->GetLabelOfParticleFirstDaughter(stackN);
  Int_t llabel = mcEvent->GetLabelOfParticleLastDaughter(stackN);
  sign = PDGMother / TMath::Abs(PDGMother);
  kDecayChannel = 0;
  for (int i = flabel; i < llabel + 1; i++) {
    AliMCParticle* cparticle = new AliMCParticle(mcEvent->GetTrack(i)->Particle());
    if (cparticle->PdgCode() == sign * fgkPdgCode[kPDGHyperHelium4]) kDecayChannel = 1;
    if (cparticle->PdgCode() == sign * fgkPdgCode[kPDGHyperHydrogen3]) kDecayChannel = 2;
    if (cparticle->PdgCode() == sign * fgkPdgCode[kPDGHyperHydrogen4]) kDecayChannel = 8;
  }
  if (!kDecayChannel) return;
    
  if (kDecayChannel == 1) {
    label1 = TMath::Abs(AliAnalysisTaskDoubleHypNucTreeLS::GetLabel(stackN, PDGFirstDaughter, sign * fgkPdgCode[kPDGHyperHelium4]));
    label2 = TMath::Abs(AliAnalysisTaskDoubleHypNucTreeLS::GetLabel(stackN, PDGSecondDaughter, sign * fgkPdgCode[kPDGHyperHelium4]));
    label3 = TMath::Abs(AliAnalysisTaskDoubleHypNucTreeLS::GetLabel(stackN, PDGThirdDaughter, sign * fgkPdgCode[kPDGHyperHelium4]));
    label4 = TMath::Abs(AliAnalysisTaskDoubleHypNucTreeLS::GetLabel(stackN, PDGFourthDaughter));
        
  }
  else {
    label1 = TMath::Abs(AliAnalysisTaskDoubleHypNucTreeLS::GetLabel(stackN, PDGFirstDaughter, sign * fgkPdgCode[kPDGHyperHydrogen3]));
    label2 = TMath::Abs(AliAnalysisTaskDoubleHypNucTreeLS::GetLabel(stackN, PDGSecondDaughter));
    label3 = TMath::Abs(AliAnalysisTaskDoubleHypNucTreeLS::GetLabel(stackN, PDGThirdDaughter));
    label4 = TMath::Abs(AliAnalysisTaskDoubleHypNucTreeLS::GetLabel(stackN, PDGFourthDaughter, sign * fgkPdgCode[kPDGHyperHydrogen3]));
  }
  FirstDaughter = new AliMCParticle(mcEvent->GetTrack(label1)->Particle());
  SecondDaughter = new AliMCParticle(mcEvent->GetTrack(label2)->Particle());
  ThirdDaughter = new AliMCParticle(mcEvent->GetTrack(label3)->Particle());
  FourthDaughter = new AliMCParticle(mcEvent->GetTrack(label4)->Particle());
    
  particle1->SetXYZM(FirstDaughter->Px(), FirstDaughter->Py(), FirstDaughter->Pz(), massFirstDaughter);
  particle2->SetXYZM(SecondDaughter->Px(), SecondDaughter->Py(), SecondDaughter->Pz(), massSecondDaughter);
  particle3->SetXYZM(ThirdDaughter->Px(), ThirdDaughter->Py(), ThirdDaughter->Pz(), massThirdDaughter);
  particle4->SetXYZM(FourthDaughter->Px(), FourthDaughter->Py(), FourthDaughter->Pz(), massFourthDaughter);
  sublorentzsum->SetXYZM(0., 0., 0., 0.);
  lorentzsum->SetXYZM(0., 0., 0., 0.);
  *lorentzsum = *particle1 + *particle2 + *particle3 + *particle4;
  if (kDecayChannel == 1) *sublorentzsum = *particle1 + *particle2 + *particle3;
  else *sublorentzsum = *particle1 + *particle4;
    
  if (FirstDaughter->PdgCode() == PDGFirstDaughter) {
    if (SecondDaughter->PdgCode() == PDGSecondDaughter) {
      if (ThirdDaughter->PdgCode() == PDGThirdDaughter) {
	if (FourthDaughter->PdgCode() == PDGFourthDaughter) {
                    
	  if (kDecayChannel == 1) {
	    dd[0] = ParticleMother->Xv() - FourthDaughter->Xv();
	    dd[1] = ParticleMother->Yv() - FourthDaughter->Yv();
	    dd[2] = ParticleMother->Zv() - FourthDaughter->Zv();
	  }
	  else {
	    dd[0] = ParticleMother->Xv() - SecondDaughter->Xv();
	    dd[1] = ParticleMother->Yv() - SecondDaughter->Yv();
	    dd[2] = ParticleMother->Zv() - SecondDaughter->Zv();
	  }
                    
	  fPDGMother = (Int_t)PDGMother;
	  fChargeMother = sign;
	  fDecayChannel = kDecayChannel;
	  fmMother = lorentzsum->M();
	  fptMother = lorentzsum->Pt();
	  fyMother = lorentzsum->Rapidity();
	  fctMother = (lorentzsum->M() * TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2))) / lorentzsum->P();
                    
	  if (kDecayChannel == 1) {
	    dd[0] = FourthDaughter->Xv() - FirstDaughter->Xv();
	    dd[1] = FourthDaughter->Yv() - FirstDaughter->Yv();
	    dd[2] = FourthDaughter->Zv() - FirstDaughter->Zv();
	  }
	  else {
	    dd[0] = SecondDaughter->Xv() - FirstDaughter->Xv();
	    dd[1] = SecondDaughter->Yv() - FirstDaughter->Yv();
	    dd[2] = SecondDaughter->Zv() - FirstDaughter->Zv();
	  }
                    
	  fmSubMother = sublorentzsum->M();
	  fptSubMother = sublorentzsum->Pt();
	  fySubMother = sublorentzsum->Rapidity();
	  fctSubMother = (sublorentzsum->M() * TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2))) / sublorentzsum->P();
	  fPrimary4LHe = 0;//mcEvent->IsSecondaryFromWeakDecay(AliAnalysisTaskDoubleHypNucTreeLS::GetLabel(stackN, sign * fgkPdgCode[kPDGHyperHelium4])) == 0;
                    
	  fptDaughter = particle1->Pt();
	  fyDaughter = particle1->Rapidity();
	  fptDaughter1 = particle2->Pt();
	  fyDaughter1 = particle2->Rapidity();
	  fptDaughter2 = particle3->Pt();
	  fyDaughter2 = particle3->Rapidity();
	  fptDaughter3 = particle4->Pt();
	  fyDaughter3 = particle4->Rapidity();
	  fTreeGen->Fill();
	  AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
	}
      }
    }
  }
  if (FirstDaughter)  delete FirstDaughter;
  if (SecondDaughter) delete SecondDaughter;
  if (ThirdDaughter)  delete ThirdDaughter;
  if (FourthDaughter) delete FourthDaughter;
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTreeLS::MCThreeBodyDecay(Int_t stackN, AliMCParticle* ParticleMother, Long_t PDGMother,
							   Long_t PDGFirstDaughter, Long_t PDGSecondDaughter, Long_t PDGThirdDaughter,
							   Double_t massFirstDaughter, Double_t massSecondDaughter, Double_t massThirdDaughter) {
    
  sign = PDGMother / TMath::Abs(PDGMother);
    
  FirstDaughter = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(AliAnalysisTaskDoubleHypNucTreeLS::GetLabel(stackN, PDGFirstDaughter)))->Particle());
  SecondDaughter = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(AliAnalysisTaskDoubleHypNucTreeLS::GetLabel(stackN, PDGSecondDaughter)))->Particle());
  ThirdDaughter = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(AliAnalysisTaskDoubleHypNucTreeLS::GetLabel(stackN, PDGThirdDaughter)))->Particle());
    
  particle1->SetXYZM(FirstDaughter->Px(), FirstDaughter->Py(), FirstDaughter->Pz(), massFirstDaughter);
  particle2->SetXYZM(SecondDaughter->Px(), SecondDaughter->Py(), SecondDaughter->Pz(), massSecondDaughter);
  particle3->SetXYZM(ThirdDaughter->Px(), ThirdDaughter->Py(), ThirdDaughter->Pz(), massThirdDaughter);
  sublorentzsum->SetXYZM(0., 0., 0., 0.);
  *sublorentzsum = *particle1 + *particle2 + *particle3;
    
  if (FirstDaughter->PdgCode() == PDGFirstDaughter) {
    if (SecondDaughter->PdgCode() == PDGSecondDaughter) {
      if (ThirdDaughter->PdgCode() == PDGThirdDaughter) {
                
	dd[0] = ParticleMother->Xv() - FirstDaughter->Xv();
	dd[1] = ParticleMother->Yv() - FirstDaughter->Yv();
	dd[2] = ParticleMother->Zv() - FirstDaughter->Zv();
                
	fPDGMother = (Int_t)PDGMother;
	fChargeMother = sign * 2;
	fPrimary4LHe = mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(stackN)) == 0;
	fisExcited = (TMath::Abs(PDGMother) == fgkPdgCode[kPDGHyperHelium4Star]);
	fDecayChannel = 1;
                
	fmSubMother = sublorentzsum->M();
	fptSubMother = sublorentzsum->Pt();
	fySubMother = sublorentzsum->Rapidity();
	fctSubMother = (sublorentzsum->M() * TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2))) / sublorentzsum->P();
                
	fptDaughter = particle1->Pt();
	fyDaughter = particle1->Rapidity();
	fptDaughter1 = particle2->Pt();
	fyDaughter1 = particle2->Rapidity();
	fptDaughter2 = particle3->Pt();
	fyDaughter2 = particle3->Rapidity();
	gTreeGen->Fill();
	AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
      }
    }
  }
  if (FirstDaughter)  delete FirstDaughter;
  if (SecondDaughter) delete SecondDaughter;
  if (ThirdDaughter)  delete ThirdDaughter;
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTreeLS::MCTwoBodyDecay(Int_t stackN, AliMCParticle* ParticleMother, Long_t PDGMother,
							 Long_t PDGFirstDaughter, Long_t PDGSecondDaughter,
							 Double_t massFirstDaughter, Double_t massSecondDaughter) {
    
  sign = PDGMother / TMath::Abs(PDGMother);
    
  FirstDaughter = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(AliAnalysisTaskDoubleHypNucTreeLS::GetLabel(stackN, PDGFirstDaughter)))->Particle());
  SecondDaughter = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(AliAnalysisTaskDoubleHypNucTreeLS::GetLabel(stackN, PDGSecondDaughter)))->Particle());
    
  particle1->SetXYZM(FirstDaughter->Px(), FirstDaughter->Py(), FirstDaughter->Pz(), massFirstDaughter);
  particle2->SetXYZM(SecondDaughter->Px(), SecondDaughter->Py(), SecondDaughter->Pz(), massSecondDaughter);
  sublorentzsum->SetXYZM(0., 0., 0., 0.);
  *sublorentzsum = *particle1 + *particle2;
    
  if (FirstDaughter->PdgCode() == PDGFirstDaughter) {
    if (SecondDaughter->PdgCode() == PDGSecondDaughter) {
            
      dd[0] = ParticleMother->Xv() - FirstDaughter->Xv();
      dd[1] = ParticleMother->Yv() - FirstDaughter->Yv();
      dd[2] = ParticleMother->Zv() - FirstDaughter->Zv();
            
      fPDGMother = (Int_t)PDGMother;
      fChargeMother = sign * 1;
      fPrimary4LHe = mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(stackN)) == 0;
      fisExcited = (TMath::Abs(PDGMother) == fgkPdgCode[kPDGHyperHydrogen4Star]);
      fDecayChannel = (PDGMother == fgkPdgCode[kPDGHyperHydrogen3]) ? 2 : 8;
            
      fmSubMother = sublorentzsum->M();
      fptSubMother = sublorentzsum->Pt();
      fySubMother = sublorentzsum->Rapidity();
      fctSubMother = (sublorentzsum->M() * TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2))) / sublorentzsum->P();
            
      fptDaughter = particle1->Pt();
      fyDaughter = particle1->Rapidity();
      fptDaughter1 = particle2->Pt();
      fyDaughter1 = particle2->Rapidity();
      gTreeGen->Fill();
      AliAnalysisTaskDoubleHypNucTreeLS::ResetVals("");
    }
  }
  if (FirstDaughter)  delete FirstDaughter;
  if (SecondDaughter) delete SecondDaughter;
}
// _________________________________________________ //
Double_t AliAnalysisTaskDoubleHypNucTreeLS::Bethe(const AliAODTrack& track, Double_t mass, Int_t charge, Double_t* params) {
  Double_t expected = charge * charge * AliExternalTrackParam::BetheBlochAleph(charge * track.GetTPCmomentum() / mass, params[0], params[1], params[2], params[3], params[4]);
  Double_t sigma = expected * params[5];
  if (TMath::IsNaN(expected)) return -999;
  return (track.GetTPCsignal() - expected) / sigma;
}
// _________________________________________________ //
Double_t AliAnalysisTaskDoubleHypNucTreeLS::GetTOFSignal(const AliAODTrack& track) {
  Float_t mass = 0;
  Float_t time = -1;
  Float_t beta = 0;
  Float_t gamma = 0;
  Float_t length = 0;
  Float_t time0 = 0;
  length = track.GetIntegratedLength();
  time0 = fPID->GetTOFResponse().GetStartTime(track.P());//fPID->GetTOFResponse().GetTimeZero();
  time = track.GetTOFsignal() - time0;
  if (time > 0) {
    beta = length / (2.99792457999999984e-02 * time);
    if (beta * beta >= 1) return -99;
    gamma = 1 / TMath::Sqrt(1 - beta * beta);
    if (gamma * gamma <= 1) return -99;
    mass = (track.GetTPCmomentum()) / TMath::Sqrt(gamma * gamma - 1); // using inner TPC mom. as approx.
  }
  return mass;
}
// _________________________________________________ //
Double_t AliAnalysisTaskDoubleHypNucTreeLS::GeoLength(const AliAODTrack& track) {
  Double_t deadZoneWidth = 3.0;
  Double_t lengthInActiveZone;
  AliExternalTrackParam *p = new AliExternalTrackParam();
  p->CopyFromVTrack(&track);
  lengthInActiveZone = GetLengthInActiveZone(p, deadZoneWidth, 220, fAODevent->GetMagneticField(), 0); //track.GetInnerParam()
  return lengthInActiveZone;
}
// _________________________________________________ //
Double_t AliAnalysisTaskDoubleHypNucTreeLS::GetLengthInActiveZone(const AliExternalTrackParam* paramT, Double_t deltaY, Double_t deltaZ, Double_t bz, Double_t exbPhi) {
  //
  // Numerical code to calculate the length of the track in active region of the TPC
  // ( can be speed up if somebody wants to invest time - analysical version shoult be possible) 
  //
  // Input parameters:
  //   paramT - external track parameters 
  //   deltaY - user defined "dead region" in cm
  //   deltaZ - user defined "active region" in cm (250 cm drift lenght - 14 cm L1 delay
  //   bz     - magnetic field 
  //   exbPhi - optional rotation due to the ExB effect
  // return value:
  //   the length of the track in cm in "active volume" of the TPC
  //
  const Double_t rIn = 85;
  const Double_t rOut = 245;
  Double_t xyz[3], pxyz[3];
  if (paramT->GetXYZAt(rIn, bz, xyz)) {
    paramT->GetPxPyPzAt(rIn, bz, pxyz);
  }
  else {
    paramT->GetXYZ(xyz);
    paramT->GetPxPyPz(pxyz);
  }
  //
  Double_t dca = -paramT->GetD(0, 0, bz);  // get impact parameter distance to point (0,0)
  Double_t radius = TMath::Abs(1 / paramT->GetC(bz));  //
  Double_t sign = paramT->GetSign() * TMath::Sign(1., bz) * (-1.);
  Double_t R0 = TMath::Sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1]);   // radius at current point
  Double_t phiR0 = TMath::ATan2(xyz[1], xyz[0]);                // angle of given point
  Double_t dPhiR0 = -TMath::ASin((dca * dca - 2 * dca * radius * sign + R0 * R0) / (2 * R0 * (dca - radius * sign)));
  Double_t phi0 = phiR0 - (dPhiR0);  // global phi offset to be added
  //
  //
  AliExternalTrackParam paramR = (*paramT);
  Double_t length = 0;
  for (Double_t R = rIn; R <= rOut; R++) {
    Double_t sinPhi = (dca * dca - 2 * dca * radius * sign + R * R) / (2 * R * (dca - radius * sign));
    if (TMath::Abs(sinPhi) >= 1) continue;
    Double_t dphi = -TMath::ASin(sinPhi);
    Double_t phi = phi0 + dphi;                           // global phi
    Int_t    sector = TMath::Nint(9 * phi / (TMath::Pi()));
    Double_t dPhiEdge = phi - (sector * TMath::Pi() / 9) + exbPhi;   // distance to sector boundary in rphi
    Double_t dX = R * TMath::Cos(phi) - xyz[0];
    Double_t dY = R * TMath::Sin(phi) - xyz[1];
    Double_t deltaPhi = 2 * TMath::ASin(0.5 * TMath::Sqrt(dX * dX + dY * dY) / radius);
    Double_t z = xyz[2] + deltaPhi * radius * paramT->GetTgl();
    if (TMath::Abs(dPhiEdge * R) > deltaY && TMath::Abs(z) < deltaZ) {
      length++;
    }
  }
  return length;
}
// _________________________________________________ //
Float_t AliAnalysisTaskDoubleHypNucTreeLS::GetInvPtDevFromBC(Int_t b, Int_t c) {
  //returns d(1/Pt) in c/GeV
  //in case of no gtu simulation -> return maximum 0.5
  if (b == 0 && c == 0) return 0.5;
  Int_t tmp = (((b & 0xfff) << 12) ^ 0x800000) - 0x800000;
  tmp += (c & 0xfff);
  Float_t invPtDev = tmp * 0.000001;
  return invPtDev;
}
// _________________________________________________ //
// __ reset variables __ //
void AliAnalysisTaskDoubleHypNucTreeLS::ResetVals(TString mode) {
  if (mode == "Event") {
    fPeriod = -1;
    fMCtrue = 0;
    frunnumber = -1;
    fMagneticField = -1;
    kMagF = -1;
    fCentrality = -1;
    feventclass = -1;
    MB = 0;
    HMV0 = 0;
    HMSPD = 0;
    HNU = 0;
    HQU = 0;
    Central = 0;
    SemiCentral = 0;
    fTrigMB = -99;
    fTrigHMV0 = -99;
    fTrigHMSPD = -99;
    fTrigHNU = -99;
    fTrigHQU = -99;
    fTrigkCentral = -99;
    fTrigkSemiCentral = -99;
    fEventID = -1;
    fisWOTrackVertex = -99;
  }
  if (mode == "") {
    fisExcited = 0;
    PrimVertex[0] = -9999;
    PrimVertex[1] = -9999;
    PrimVertex[2] = -9999;
    SecVertex[0] = -9999;
    SecVertex[1] = -9999;
    SecVertex[2] = -9999;
    TertVertex[0] = -9999;
    TertVertex[1] = -9999;
    TertVertex[2] = -9999;
    fPrimVertexX = -9999;
    fPrimVertexY = -9999;
    fPrimVertexZ = -9999;
    fSecVertexX = -9999;
    fSecVertexY = -9999;
    fSecVertexZ = -9999;
    fTertVertexX = -9999;
    fTertVertexY = -9999;
    fTertVertexZ = -9999;
    fV0VertexX_13 = -9999;
    fV0VertexY_13 = -9999;
    fV0VertexZ_13 = -9999;
    fV0VertexX_14 = -9999;
    fV0VertexY_14 = -9999;
    fV0VertexZ_14 = -9999;
    fV0VertexX_23 = -9999;
    fV0VertexY_23 = -9999;
    fV0VertexZ_23 = -9999;
    fV0VertexX_24 = -9999;
    fV0VertexY_24 = -9999;
    fV0VertexZ_24 = -9999;
    fisOnlineV0_13 = 0;
    fisOnlineV0_14 = 0;
    fisOnlineV0_23 = 0;
    fisOnlineV0_24 = 0;
    fDCA2B = -99;
    fDCA2Bo = -99;
    fDCA3B1 = -99;
    fDCA3B2 = -99;
    fDCA3B3 = -99;
    fPA = -99;
    fSubPA = -99;
    fSubPA2 = -99;
    fDecAngle = -99;
    farmalpha = -99;
    farmpt = -99;
    fPDGMother = -99;
    fRecoMethod = -99;
    fDecayChannel = -99;
    fChargeMother = -99;
    fmctruth = 0;
    fPrimary4LHe = 0;
    fmMother = -99;
    fmMother2 = -99;
    fEMother = -99;
    fpxMother = -99;
    fpyMother = -99;
    fpzMother = -99;
    fptMother = -99;
    fpMother = -99;
    fyMother = -99;
    fctMother = -99;
    fmSubMother = -99;
    fESubMother = -99;
    fpxSubMother = -99;
    fpySubMother = -99;
    fpzSubMother = -99;
    fptSubMother = -99;
    fpSubMother = -99;
    fySubMother = -99;
    fctSubMother = -99;
    fEDaughter = -99;
    fpDaughter = -99;
    fptDaughter = -99;
    fpxDaughter = -99;
    fpyDaughter = -99;
    fpzDaughter = -99;
    fyDaughter = -99;
    fEDaughterKF = -99;
    fLabelDaughterKF = -99;
    fptDaughterKF = -99;
    fpxDaughterKF = -99;
    fpyDaughterKF = -99;
    fpzDaughterKF = -99;
    fyDaughterKF = -99;
    fdEdxDaughter = -99;
    fdEdxSigmaDaughter = -99;
    fDcaDaughter = -99;
    fDcaDaughtero = -99;
    fDcazDaughter = -99;
    fDcaSecDaughter = -99;
    fImParDaughter = -99;
    fImParzDaughter = -99;
    fNclsDaughter = -99;
    fChi2Daughter = -99;
    fNclsITSDaughter = -99;
    fEtaDaughter = -99;
    fPhiDaughter = -99;
    fGeoLengthDaughter = -99;
    fTOFSignalDaughter = -99;
    fSigmaYXDaughter = -99;
    fSigmaXYZDaughter = -99;
    fSigmaZDaughter = -99;
    fPtUncertDaughter = -1;
    fTPCRefitDaughter = -1;
    fITSRefitDaughter = -1;
    fPropDCADaughter = -99;
    fEDaughter1KF = -99;
    fLabelDaughter1KF = -99;
    fptDaughter1KF = -99;
    fpxDaughter1KF = -99;
    fpyDaughter1KF = -99;
    fpzDaughter1KF = -99;
    fyDaughter1KF = -99;
    fEDaughter1 = -99;
    fpDaughter1 = -99;
    fptDaughter1 = -99;
    fpxDaughter1 = -99;
    fpyDaughter1 = -99;
    fpzDaughter1 = -99;
    fyDaughter1 = -99;
    fdEdxDaughter1 = -99;
    fdEdxSigmaDaughter1 = -99;
    fDcaDaughter1 = -99;
    fDcaDaughter1o = -99;
    fDcazDaughter1 = -99;
    fDcaSecDaughter1 = -99;
    fImParDaughter1 = -99;
    fImParzDaughter1 = -99;
    fNclsDaughter1 = -99;
    fChi2Daughter1 = -99;
    fNclsITSDaughter1 = -99;
    fEtaDaughter1 = -99;
    fPhiDaughter1 = -99;
    fGeoLengthDaughter1 = -99;
    fTOFSignalDaughter1 = -99;
    fSigmaYXDaughter1 = -99;
    fSigmaXYZDaughter1 = -99;
    fSigmaZDaughter1 = -99;
    fPtUncertDaughter1 = -1;
    fTPCRefitDaughter1 = -1;
    fITSRefitDaughter1 = -1;
    fPropDCADaughter1 = -99;
    fEDaughter2KF = -99;
    fLabelDaughter2KF = -99;
    fptDaughter2KF = -99;
    fpxDaughter2KF = -99;
    fpyDaughter2KF = -99;
    fpzDaughter2KF = -99;
    fyDaughter2KF = -99;
    fEDaughter2 = -99;
    fpDaughter2 = -99;
    fptDaughter2 = -99;
    fpxDaughter2 = -99;
    fpyDaughter2 = -99;
    fpzDaughter2 = -99;
    fyDaughter2 = -99;
    fdEdxDaughter2 = -99;
    fdEdxSigmaDaughter2 = -99;
    fDcaDaughter2 = -99;
    fDcaDaughter2o = -99;
    fDcazDaughter2 = -99;
    fDcaSecDaughter2 = -99;
    fImParDaughter2 = -99;
    fImParzDaughter2 = -99;
    fNclsDaughter2 = -99;
    fChi2Daughter2 = -99;
    fNclsITSDaughter2 = -99;
    fEtaDaughter2 = -99;
    fPhiDaughter2 = -99;
    fGeoLengthDaughter2 = -99;
    fTOFSignalDaughter2 = -99;
    fSigmaYXDaughter2 = -99;
    fSigmaXYZDaughter2 = -99;
    fSigmaZDaughter2 = -99;
    fPtUncertDaughter2 = -1;
    fTPCRefitDaughter2 = -1;
    fITSRefitDaughter2 = -1;
    fPropDCADaughter2 = -99;
    fEDaughter3KF = -99;
    fLabelDaughter3KF = -99;
    fptDaughter3KF = -99;
    fpxDaughter3KF = -99;
    fpyDaughter3KF = -99;
    fpzDaughter3KF = -99;
    fyDaughter3KF = -99;
    fEDaughter3 = -99;
    fpDaughter3 = -99;
    fptDaughter3 = -99;
    fpxDaughter3 = -99;
    fpyDaughter3 = -99;
    fpzDaughter3 = -99;
    fyDaughter3 = -99;
    fdEdxDaughter3 = -99;
    fdEdxSigmaDaughter3 = -99;
    fDcaDaughter3 = -99;
    fDcaDaughter3o = -99;
    fDcazDaughter3 = -99;
    fDcaSecDaughter3 = -99;
    fImParDaughter3 = -99;
    fImParzDaughter3 = -99;
    fNclsDaughter3 = -99;
    fChi2Daughter3 = -99;
    fNclsITSDaughter3 = -99;
    fEtaDaughter3 = -99;
    fPhiDaughter3 = -99;
    fGeoLengthDaughter3 = -99;
    fTOFSignalDaughter3 = -99;
    fSigmaYXDaughter3 = -99;
    fSigmaXYZDaughter3 = -99;
    fSigmaZDaughter3 = -99;
    fPtUncertDaughter3 = -1;
    fTPCRefitDaughter3 = -1;
    fITSRefitDaughter3 = -1;
    fPropDCADaughter3 = -99;
    fEDaughter4 = -99;
    fpDaughter4 = -99;
    fptDaughter4 = -99;
    fpxDaughter4 = -99;
    fpyDaughter4 = -99;
    fpzDaughter4 = -99;
    fyDaughter4 = -99;
    fDcaDaughter4 = -99;
    fDcazDaughter4 = -99;
    fDcaSecDaughter4 = -99;
    fImParDaughter4 = -99;
    fImParzDaughter4 = -99;
    fSigmaYXDaughter4 = -99;
    fSigmaXYZDaughter4 = -99;
    fSigmaZDaughter4 = -99;
    fPtUncertDaughter4 = -1;
    fPropDCADaughter4 = -99;
    fthetaP = -99;
    fthetaN = -99;
    fdEdxSigmaPion = -99;
    fdEdxSigmaDeuteron = -99;
    fdEdxSigmaTriton = -99;
    fdEdxSigmaAlpha = -99;
    PrimVertexKF[0] = 0.;
    PrimVertexKF[1] = 0.;
    PrimVertexKF[2] = 0.;
    SecVertexKF[0] = 0.;
    SecVertexKF[1] = 0.;
    SecVertexKF[2] = 0.;
    TertVertexKF[0] = 0.;
    TertVertexKF[1] = 0.;
    TertVertexKF[2] = 0.;
    fPrimVertexXKF = -9999;
    fPrimVertexYKF = -9999;
    fPrimVertexZKF = -9999;
    fSecVertexXKF = -9999;
    fSecVertexYKF = -9999;
    fSecVertexZKF = -9999;
    fTertVertexXKF = -9999;
    fTertVertexYKF = -9999;
    fTertVertexZKF = -9999;
    fPrimVertexXErrKF = -9999;
    fPrimVertexYErrKF = -9999;
    fPrimVertexZErrKF = -9999;
    fSecVertexXErrKF = -9999;
    fSecVertexYErrKF = -9999;
    fSecVertexZErrKF = -9999;
    fTertVertexXErrKF = -9999;
    fTertVertexYErrKF = -9999;
    fTertVertexZErrKF = -9999;
    fTertVertChi2KF = -99;
    fTertVertNDFKF = -99;
    fSecVertChi2KF = -99;
    fSecVertNDFKF = -99;
    fPrimVertChi2KF = -99;
    fPrimVertNDFKF = -99;
    fPAKF = -99;
    fSubPAKF = -99;
    fSubPA2KF = -99;
    fmMotherKF = -99;
    fmMotherErrKF = -99;
    fEMotherKF = -99;
    fEMotherErrKF = -99;
    fpxMotherKF = -99;
    fpxMotherErrKF = -99;
    fpyMotherKF = -99;
    fpyMotherErrKF = -99;
    fpzMotherKF = -99;
    fpzMotherErrKF = -99;
    fptMotherKF = -99;
    fptMotherErrKF = -99;
    fpMotherKF = -99;
    fpMotherErrKF = -99;
    fyMotherKF = -99;
    fctMotherKF = -99;
    fctMotherErrKF = -99;
    fmSubMotherKF = -99;
    fmSubMotherErrKF = -99;
    fESubMotherKF = -99;
    fESubMotherErrKF = -99;
    fpxSubMotherKF = -99;
    fpxSubMotherErrKF = -99;
    fpySubMotherKF = -99;
    fpySubMotherErrKF = -99;
    fpzSubMotherKF = -99;
    fpzSubMotherErrKF = -99;
    fptSubMotherKF = -99;
    fptSubMotherErrKF = -99;
    fpSubMotherKF = -99;
    fpSubMotherErrKF = -99;
    fySubMotherKF = -99;
    fctSubMotherKF = -99;
    fctSubMotherErrKF = -99;
    fDcaDaughterXYKF = -99;
    fDcaDaughter1XYKF = -99;
    fDcaDaughter2XYKF = -99;
    fDcaDaughter3XYKF = -99;
    fDcaDaughter4XYKF = -99;
    fDcaSecDaughterXYKF = -99;
    fDcaSecDaughter1XYKF = -99;
    fDcaSecDaughter2XYKF = -99;
    fDcaSecDaughter3XYKF = -99;
    fDcaSecDaughter4XYKF = -99;
    fDcaDaughterZKF = -99;
    fDcaDaughter1ZKF = -99;
    fDcaDaughter2ZKF = -99;
    fDcaDaughter3ZKF = -99;
    fDcaDaughter4ZKF = -99;
    fDcaSecDaughterZKF = -99;
    fDcaSecDaughter1ZKF = -99;
    fDcaSecDaughter2ZKF = -99;
    fDcaSecDaughter3ZKF = -99;
    fDcaSecDaughter4ZKF = -99;
    fDCA2BXYKF = -99;
    fDCA3B1XYKF = -99;
    fDCA3B2XYKF = -99;
    fDCA3B3XYKF = -99;
    fDCA2BZKF = -99;
    fDCA3B1ZKF = -99;
    fDCA3B2ZKF = -99;
    fDCA3B3ZKF = -99;
    fITSLayer1Daughter = 0;
    fITSLayer2Daughter = 0;
    fITSLayer3Daughter = 0;
    fITSLayer4Daughter = 0;
    fITSLayer5Daughter = 0;
    fITSLayer6Daughter = 0;
    fITSLayer1Daughter1 = 0;
    fITSLayer2Daughter1 = 0;
    fITSLayer3Daughter1 = 0;
    fITSLayer4Daughter1 = 0;
    fITSLayer5Daughter1 = 0;
    fITSLayer6Daughter1 = 0;
    fITSLayer1Daughter2 = 0;
    fITSLayer2Daughter2 = 0;
    fITSLayer3Daughter2 = 0;
    fITSLayer4Daughter2 = 0;
    fITSLayer5Daughter2 = 0;
    fITSLayer6Daughter2 = 0;
    fITSLayer1Daughter3 = 0;
    fITSLayer2Daughter3 = 0;
    fITSLayer3Daughter3 = 0;
    fITSLayer4Daughter3 = 0;
    fITSLayer5Daughter3 = 0;
    fITSLayer6Daughter3 = 0;
    fptDaughterUnProp = -99;
    fpxDaughterUnProp = -99;
    fpyDaughterUnProp = -99;
    fpzDaughterUnProp = -99;
    fptDaughter1UnProp = -99;
    fpxDaughter1UnProp = -99;
    fpyDaughter1UnProp = -99;
    fpzDaughter1UnProp = -99;
    fptDaughter2UnProp = -99;
    fpxDaughter2UnProp = -99;
    fpyDaughter2UnProp = -99;
    fpzDaughter2UnProp = -99;
    fptDaughter3UnProp = -99;
    fpxDaughter3UnProp = -99;
    fpyDaughter3UnProp = -99;
    fpzDaughter3UnProp = -99;
    fTrackPIDDaughter = -99;
    fTrackPIDDaughter1 = -99;
    fTrackPIDDaughter2 = -99;
    fTrackPIDDaughter3 = -99;
    std::fill_n(fCovMatrixTrack, 21, -999);
    std::fill_n(fCovMatrixTrack1, 21, -999);
    std::fill_n(fCovMatrixTrack2, 21, -999);
    std::fill_n(fCovMatrixTrack3, 21, -999);
    std::fill_n(fCovMatrixTrack4, 21, -999);
    std::fill_n(fTrackPar, 7, -999);
    std::fill_n(fTrackPar1, 7, -999);
    std::fill_n(fTrackPar2, 7, -999);
    std::fill_n(fTrackPar3, 7, -999);
    std::fill_n(fTrackPar4, 7, -999);
  }
  return;
}
// _________________________________________________ //
// __ set Bethe-Bloch parameter __ //
void AliAnalysisTaskDoubleHypNucTreeLS::SetBetheBlochParams(Int_t runNumber) {
  if (runNumber == 170593) {
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
  else if (runNumber < 246994) return;
  // __ 2015 PbPb __ //
  if (runNumber >= 244917 && runNumber <= 246994) {
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
  else if (runNumber >= 295581 && runNumber <= 297624) {
    fBetheParamsT[0] = 0.648689;
    fBetheParamsT[1] = 56.6706;
    fBetheParamsT[2] = -1.63243e-10;
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
  else if (runNumber >= 252235 && runNumber <= 264347) {
    if (!fMCtrue) { // Data
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
    }
    else { // MC
      if (runNumber >= 262424 || runNumber <= 256418) {
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
      if (runNumber >= 256941 && runNumber <= 258537) {
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
      if (runNumber >= 258962 && runNumber <= 259888) {
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
    if (!fMCtrue) {
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
    }
    else {
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
    if (!fMCtrue) {
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
    }
    else {
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
