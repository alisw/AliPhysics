//--- Task for investigation of the {}^{4}_{#Lambda#Lambda}H ---
//---     Author: Janik Ditzel; janik.ditzel@cern.ch     ---


#include "AliAnalysisTaskDoubleHypNucTree.h"

ClassImp(AliAnalysisTaskDoubleHypNucTree)
// _________________________________________________ //
// Default Constructor
AliAnalysisTaskDoubleHypNucTree::AliAnalysisTaskDoubleHypNucTree()
:AliAnalysisTaskSE("AliAnalysisTaskDoubleHypNucTree"),
  fPIDCheckOnly(kFALSE), fInputHandler(0), fPID(0), fESDevent(0), mcEvent(0), fStack(), fEventCuts(),
//
  fTriggerMask(), MB(0), HMV0(0), HMSPD(0), HNU(0), HQU(0), Central(0), SemiCentral(0),
//
  fTree(0), fTreeD(0), fTreeKF(0), gTree(0), gTreeD(0), gTreeKF(0), fTreeGen(0), gTreeGen(0),
//
  kVariante(1), kStandardReco(0), kKFReco(0),
//
  He3PosArray(), PPosArray(), PiPosArray(), PiPosSecArray(), He3NegArray(), PNegArray(), PiNegArray(), PiNegSecArray(), He3PosCounter(), He3NegCounter(), PPosCounter(), PNegCounter(), PiPosCounter(), PiNegCounter(), PiPosSecCounter(), PiNegSecCounter(),
//
  fBetheParamsHe(), fBetheParamsT(),
//
  fMCtrue(0), konlyBG(0), konlySig(0), kMCPIDCheck(0), fPrimary4LHe(0), stackN(0),
//
  vertex(), primVertex(), secvertexer(), secVertex(), tertvertexer(), tertVertex(), PrimVertex(), SecVertex(), TertVertex(),
//
  track1(), track2(), track3(), track4(), exTrack(), exTrack1(), exTrack2(), exTrack3(), exTrack4(), trkArray(), trkArray1(), 
//
  kDCATracksCut(0), kTrackPtUncertCut(0), kPointingAngleCut(0),  kPIDSelecCut(0), kMin4LLHMass(0), kMax4LLHMass(0), kMin4LHeMass(0), kMax4LHeMass(0), kMin3LHMass(0), kMax3LHMass(0),
//
  kMagF(0),
//
  cov(), cov0(), cov1(), cov2(), pxpypz(), xyz(), sign(0), dd(), xthiss(), xpp(), h(), fthetaP(-99), fthetaN(-99),
//  
  PDGCodeMother(), ParticleMother(NULL), FirstDaughter(NULL), SecondDaughter(NULL), ThirdDaughter(NULL), FourthDaughter(NULL), label1(), labelMother1(), labelGrandMother1(), ParticleMother1(NULL), ParticleGrandMother1(NULL), label2(), labelMother2(), labelGrandMother2(), ParticleMother2(NULL), ParticleGrandMother2(NULL), label3(), labelMother3(), labelGrandMother3(), ParticleMother3(NULL), ParticleGrandMother3(NULL), label4(), labelMother4(), labelGrandMother4(), ParticleMother4(NULL), ParticleGrandMother4(NULL), CheckParticle(NULL),
//
  lorentzsum(NULL), lorentzsum2(NULL), sublorentzsum(NULL), sublorentzsum2(NULL), particle1(NULL), particle2(NULL), particle3(NULL), particle4(NULL),
//Histos
  fHistogramList(NULL), fHistdEdx(0), fHistNumEvents(0), fHistTrigger(0),
//saved in tree
  fCentrality(-99), frunnumber(-99), fPeriod(0), fMagneticField(0), feventclass(0), fmctruth(0), fDecayChannel(0), fRecoMethod(0),
  fTrigMB(-99), fTrigHMV0(-99), fTrigHMSPD(-99), fTrigHNU(0), fTrigHQU(0), fTrigkCentral(0), fTrigkSemiCentral(0),
//
  fisOnlineV0_13(-1), fisOnlineV0_14(-1), fisOnlineV0_23(-1), fisOnlineV0_24(-1),
//
  fDCA2B(-99), fDCA2Bo(-99), fDCA3B1(-99), fDCA3B2(-99), fDCA3B3(-99), fDecAngle(-99), farmalpha(-99), farmpt(-99), fPA(-99), fSubPA(-99), fSubPA2(-99),
//
  fPrimVertexX(-99), fPrimVertexY(-99), fPrimVertexZ(-99), fSecVertexX(-99), fSecVertexY(-99), fSecVertexZ(-99), fTertVertexX(-99), fTertVertexY(-99), fTertVertexZ(-99), fTertVertChi2(-99), fTertVertNDF(-99), fSecVertChi2(-99), fSecVertNDF(-99), fPrimVertChi2(-99), fPrimVertNDF(-99),  fV0VertexX_13(-99), fV0VertexY_13(-99), fV0VertexZ_13(-99), fV0VertexX_14(-99), fV0VertexY_14(-99), fV0VertexZ_14(-99), fV0VertexX_23(-99), fV0VertexY_23(-99), fV0VertexZ_23(-99), fV0VertexX_24(-99), fV0VertexY_24(-99), fV0VertexZ_24(-99),
//
  fPDGMother(-99), fChargeMother(-99), fmMother(-99), fmMother2(-99), fEMother(-99), fpxMother(-99), fpyMother(-99), fpzMother(-99), fptMother(-99), fpMother(-99), fyMother(-99), fctMother(-99),
//
  fmSubMother(-99), fESubMother(-99), fpxSubMother(-99), fpySubMother(-99), fpzSubMother(-99), fptSubMother(-99), fpSubMother(-99), fySubMother(-99), fctSubMother(-99),
//
  fEDaughter(-99), fpDaughter(-99), fptDaughter(-99), fpxDaughter(-99), fpyDaughter(-99), fpzDaughter(-99), fyDaughter(-99), fdEdxDaughter(-99), fdEdxSigmaDaughter(-99), fDcaDaughter(-99), fDcaDaughtero(-99), fDcazDaughter(-99), fDcaSecDaughter(-99), fImParDaughter(-99), fImParzDaughter(-99), fNclsDaughter(-99), fChi2Daughter(-99), fNclsITSDaughter(-99), fEtaDaughter(-99), fPhiDaughter(-99), fGeoLengthDaughter(-99), fTOFSignalDaughter(-99), fSigmaYXDaughter(-99), fSigmaXYZDaughter(-99), fSigmaZDaughter(-99), fPtUncertDaughter(-99), fTPCRefitDaughter(-99), fITSRefitDaughter(-99), fPropDCADaughter(-1),
//
  fEDaughter1(-99), fpDaughter1(-99), fptDaughter1(-99), fpxDaughter1(-99), fpyDaughter1(-99), fpzDaughter1(-99), fyDaughter1(-99), fdEdxDaughter1(-99), fdEdxSigmaDaughter1(-99), fDcaDaughter1(-99), fDcaDaughter1o(-99), fDcazDaughter1(-99), fDcaSecDaughter1(-99), fImParDaughter1(-99), fImParzDaughter1(-99), fNclsDaughter1(-99), fChi2Daughter1(-99), fNclsITSDaughter1(-99), fEtaDaughter1(-99), fPhiDaughter1(-99), fGeoLengthDaughter1(-99), fTOFSignalDaughter1(-99), fSigmaYXDaughter1(-99), fSigmaXYZDaughter1(-99), fSigmaZDaughter1(-99), fPtUncertDaughter1(-99), fTPCRefitDaughter1(-99), fITSRefitDaughter1(-99), fPropDCADaughter1(-1),
//
  fEDaughter2(-99), fpDaughter2(-99), fptDaughter2(-99), fpxDaughter2(-99), fpyDaughter2(-99), fpzDaughter2(-99), fyDaughter2(-99), fdEdxDaughter2(-99), fdEdxSigmaDaughter2(-99), fDcaDaughter2(-99), fDcaDaughter2o(-99), fDcazDaughter2(-99), fDcaSecDaughter2(-99), fImParDaughter2(-99), fImParzDaughter2(-99),  fNclsDaughter2(-99), fChi2Daughter2(-99), fNclsITSDaughter2(-99), fEtaDaughter2(-99), fPhiDaughter2(-99), fGeoLengthDaughter2(-99), fTOFSignalDaughter2(-99), fSigmaYXDaughter2(-99), fSigmaXYZDaughter2(-99), fSigmaZDaughter2(-99), fPtUncertDaughter2(-99), fTPCRefitDaughter2(-99), fITSRefitDaughter2(-99), fPropDCADaughter2(-1),
//
  fEDaughter3(-99), fpDaughter3(-99), fptDaughter3(-99), fpxDaughter3(-99), fpyDaughter3(-99), fpzDaughter3(-99), fyDaughter3(-99), fdEdxDaughter3(-99), fdEdxSigmaDaughter3(-99), fDcaDaughter3(-99), fDcaDaughter3o(-99), fDcazDaughter3(-99), fDcaSecDaughter3(-99), fImParDaughter3(-99), fImParzDaughter3(-99), fNclsDaughter3(-99), fChi2Daughter3(-99), fNclsITSDaughter3(-99), fEtaDaughter3(-99), fPhiDaughter3(-99), fGeoLengthDaughter3(-99), fTOFSignalDaughter3(-99), fSigmaYXDaughter3(-99), fSigmaXYZDaughter3(-99), fSigmaZDaughter3(-99), fPtUncertDaughter3(-99), fTPCRefitDaughter3(-99), fITSRefitDaughter3(-99), fPropDCADaughter3(-1),
//
  fEDaughter4(-99), fpDaughter4(-99), fptDaughter4(-99), fpxDaughter4(-99), fpyDaughter4(-99), fpzDaughter4(-99), fyDaughter4(-99), fDcaDaughter4(-99), fDcazDaughter4(-99), fDcaSecDaughter4(-99), fImParDaughter4(-99),  fImParzDaughter4(-99), fSigmaYXDaughter4(-99), fSigmaXYZDaughter4(-99), fSigmaZDaughter4(-99), fPtUncertDaughter4(-99), fPropDCADaughter4(-1),
//
  fITSLayer1Daughter(-1), fITSLayer2Daughter(-1), fITSLayer3Daughter(-1), fITSLayer4Daughter(-1), fITSLayer5Daughter(-1), fITSLayer6Daughter(-1), fITSLayer1Daughter1(-1), fITSLayer2Daughter1(-1), fITSLayer3Daughter1(-1), fITSLayer4Daughter1(-1), fITSLayer5Daughter1(-1), fITSLayer6Daughter1(-1), fITSLayer1Daughter2(-1), fITSLayer2Daughter2(-1), fITSLayer3Daughter2(-1), fITSLayer4Daughter2(-1), fITSLayer5Daughter2(-1), fITSLayer6Daughter2(-1), fITSLayer1Daughter3(-1), fITSLayer2Daughter3(-1), fITSLayer3Daughter3(-1), fITSLayer4Daughter3(-1), fITSLayer5Daughter3(-1), fITSLayer6Daughter3(-1),
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
  fDcaDaughterXYKF(-99), fDcaDaughter1XYKF(-99), fDcaDaughter2XYKF(-99), fDcaDaughter3XYKF(-99), fDcaDaughter4XYKF(-99), fDcaSecDaughterXYKF(-99), fDcaSecDaughter1XYKF(-99), fDcaSecDaughter2XYKF(-99), fDcaSecDaughter3XYKF(-99), fDcaSecDaughter4XYKF(-99), fDcaDaughterZKF(-99), fDcaDaughter1ZKF(-99), fDcaDaughter2ZKF(-99), fDcaDaughter3ZKF(-99), fDcaDaughter4ZKF(-99), fDcaSecDaughterZKF(-99), fDcaSecDaughter1ZKF(-99), fDcaSecDaughter2ZKF(-99), fDcaSecDaughter3ZKF(-99), fDcaSecDaughter4ZKF(-99),
//
  fDCA2BXYKF(-99), fDCA3B1XYKF(-99), fDCA3B2XYKF(-99), fDCA3B3XYKF(-99), fDCA2BZKF(-99), fDCA3B1ZKF(-99), fDCA3B2ZKF(-99), fDCA3B3ZKF(-99)
{

}
// _________________________________________________ //
// Constructor
AliAnalysisTaskDoubleHypNucTree::AliAnalysisTaskDoubleHypNucTree(const char* name)
  :AliAnalysisTaskSE(name),
   fPIDCheckOnly(kFALSE), fInputHandler(0), fPID(0), fESDevent(0), mcEvent(0), fStack(), fEventCuts(),
   //
   fTriggerMask(), MB(0), HMV0(0), HMSPD(0), HNU(0), HQU(0), Central(0), SemiCentral(0),
   //
   fTree(0), fTreeD(0), fTreeKF(0), gTree(0), gTreeD(0), gTreeKF(0), fTreeGen(0), gTreeGen(0),
   //
   kVariante(1), kStandardReco(0), kKFReco(0),
   //
   He3PosArray(), PPosArray(), PiPosArray(), PiPosSecArray(), He3NegArray(), PNegArray(), PiNegArray(), PiNegSecArray(), He3PosCounter(), He3NegCounter(), PPosCounter(), PNegCounter(), PiPosCounter(), PiNegCounter(), PiPosSecCounter(), PiNegSecCounter(),
   //
   fBetheParamsHe(), fBetheParamsT(),
   //
   fMCtrue(0), konlyBG(0), konlySig(0), kMCPIDCheck(0), fPrimary4LHe(0), stackN(0),
   //
   vertex(), primVertex(), secvertexer(), secVertex(), tertvertexer(), tertVertex(), PrimVertex(), SecVertex(), TertVertex(),
   //
   track1(), track2(), track3(), track4(), exTrack(), exTrack1(), exTrack2(), exTrack3(), exTrack4(), trkArray(), trkArray1(), 
   //
   kDCATracksCut(0), kTrackPtUncertCut(0), kPointingAngleCut(0),  kPIDSelecCut(0), kMin4LLHMass(0), kMax4LLHMass(0), kMin4LHeMass(0), kMax4LHeMass(0), kMin3LHMass(0), kMax3LHMass(0),
   //
   kMagF(0),
   //
   cov(), cov0(), cov1(), cov2(), pxpypz(), xyz(), sign(0), dd(), xthiss(), xpp(), h(), fthetaP(-99), fthetaN(-99),
   //  
   PDGCodeMother(), ParticleMother(NULL), FirstDaughter(NULL), SecondDaughter(NULL), ThirdDaughter(NULL), FourthDaughter(NULL), label1(), labelMother1(), labelGrandMother1(), ParticleMother1(NULL), ParticleGrandMother1(NULL), label2(), labelMother2(), labelGrandMother2(), ParticleMother2(NULL), ParticleGrandMother2(NULL), label3(), labelMother3(), labelGrandMother3(), ParticleMother3(NULL), ParticleGrandMother3(NULL), label4(), labelMother4(), labelGrandMother4(), ParticleMother4(NULL), ParticleGrandMother4(NULL), CheckParticle(NULL),
   //
   lorentzsum(NULL), lorentzsum2(NULL), sublorentzsum(NULL), sublorentzsum2(NULL), particle1(NULL), particle2(NULL), particle3(NULL), particle4(NULL),
   //Histos
   fHistogramList(NULL), fHistdEdx(0), fHistNumEvents(0), fHistTrigger(0),
   //saved in tree
   fCentrality(-99), frunnumber(-99), fPeriod(0), fMagneticField(0), feventclass(0), fmctruth(0), fDecayChannel(0), fRecoMethod(0),
   fTrigMB(-99), fTrigHMV0(-99), fTrigHMSPD(-99), fTrigHNU(0), fTrigHQU(0), fTrigkCentral(0), fTrigkSemiCentral(0), 
   //
   fisOnlineV0_13(-1), fisOnlineV0_14(-1), fisOnlineV0_23(-1), fisOnlineV0_24(-1),
   //
   fDCA2B(-99), fDCA2Bo(-99), fDCA3B1(-99), fDCA3B2(-99), fDCA3B3(-99), fDecAngle(-99), farmalpha(-99), farmpt(-99), fPA(-99), fSubPA(-99), fSubPA2(-99),
   //
   fPrimVertexX(-99), fPrimVertexY(-99), fPrimVertexZ(-99), fSecVertexX(-99), fSecVertexY(-99), fSecVertexZ(-99), fTertVertexX(-99), fTertVertexY(-99), fTertVertexZ(-99), fTertVertChi2(-99), fTertVertNDF(-99), fSecVertChi2(-99), fSecVertNDF(-99), fPrimVertChi2(-99), fPrimVertNDF(-99),  fV0VertexX_13(-99), fV0VertexY_13(-99), fV0VertexZ_13(-99), fV0VertexX_14(-99), fV0VertexY_14(-99), fV0VertexZ_14(-99), fV0VertexX_23(-99), fV0VertexY_23(-99), fV0VertexZ_23(-99), fV0VertexX_24(-99), fV0VertexY_24(-99), fV0VertexZ_24(-99),
   //
   fPDGMother(-99), fChargeMother(-99), fmMother(-99), fmMother2(-99), fEMother(-99), fpxMother(-99), fpyMother(-99), fpzMother(-99), fptMother(-99), fpMother(-99), fyMother(-99), fctMother(-99),
   //
   fmSubMother(-99), fESubMother(-99), fpxSubMother(-99), fpySubMother(-99), fpzSubMother(-99), fptSubMother(-99), fpSubMother(-99), fySubMother(-99), fctSubMother(-99),
   //
   fEDaughter(-99), fpDaughter(-99), fptDaughter(-99), fpxDaughter(-99), fpyDaughter(-99), fpzDaughter(-99), fyDaughter(-99), fdEdxDaughter(-99), fdEdxSigmaDaughter(-99), fDcaDaughter(-99), fDcaDaughtero(-99), fDcazDaughter(-99), fDcaSecDaughter(-99), fImParDaughter(-99), fImParzDaughter(-99), fNclsDaughter(-99), fChi2Daughter(-99), fNclsITSDaughter(-99), fEtaDaughter(-99), fPhiDaughter(-99), fGeoLengthDaughter(-99), fTOFSignalDaughter(-99), fSigmaYXDaughter(-99), fSigmaXYZDaughter(-99), fSigmaZDaughter(-99), fPtUncertDaughter(-99), fTPCRefitDaughter(-99), fITSRefitDaughter(-99), fPropDCADaughter(-1),
   //
   fEDaughter1(-99), fpDaughter1(-99), fptDaughter1(-99), fpxDaughter1(-99), fpyDaughter1(-99), fpzDaughter1(-99), fyDaughter1(-99), fdEdxDaughter1(-99), fdEdxSigmaDaughter1(-99), fDcaDaughter1(-99), fDcaDaughter1o(-99), fDcazDaughter1(-99), fDcaSecDaughter1(-99), fImParDaughter1(-99), fImParzDaughter1(-99), fNclsDaughter1(-99), fChi2Daughter1(-99), fNclsITSDaughter1(-99), fEtaDaughter1(-99), fPhiDaughter1(-99), fGeoLengthDaughter1(-99), fTOFSignalDaughter1(-99), fSigmaYXDaughter1(-99), fSigmaXYZDaughter1(-99), fSigmaZDaughter1(-99), fPtUncertDaughter1(-99), fTPCRefitDaughter1(-99), fITSRefitDaughter1(-99), fPropDCADaughter1(-1),
   //
   fEDaughter2(-99), fpDaughter2(-99), fptDaughter2(-99), fpxDaughter2(-99), fpyDaughter2(-99), fpzDaughter2(-99), fyDaughter2(-99), fdEdxDaughter2(-99), fdEdxSigmaDaughter2(-99), fDcaDaughter2(-99), fDcaDaughter2o(-99), fDcazDaughter2(-99), fDcaSecDaughter2(-99), fImParDaughter2(-99), fImParzDaughter2(-99),  fNclsDaughter2(-99), fChi2Daughter2(-99), fNclsITSDaughter2(-99), fEtaDaughter2(-99), fPhiDaughter2(-99), fGeoLengthDaughter2(-99), fTOFSignalDaughter2(-99), fSigmaYXDaughter2(-99), fSigmaXYZDaughter2(-99), fSigmaZDaughter2(-99), fPtUncertDaughter2(-99), fTPCRefitDaughter2(-99), fITSRefitDaughter2(-99), fPropDCADaughter2(-1),
   //
   fEDaughter3(-99), fpDaughter3(-99), fptDaughter3(-99), fpxDaughter3(-99), fpyDaughter3(-99), fpzDaughter3(-99), fyDaughter3(-99), fdEdxDaughter3(-99), fdEdxSigmaDaughter3(-99), fDcaDaughter3(-99), fDcaDaughter3o(-99), fDcazDaughter3(-99), fDcaSecDaughter3(-99), fImParDaughter3(-99), fImParzDaughter3(-99), fNclsDaughter3(-99), fChi2Daughter3(-99), fNclsITSDaughter3(-99), fEtaDaughter3(-99), fPhiDaughter3(-99), fGeoLengthDaughter3(-99), fTOFSignalDaughter3(-99), fSigmaYXDaughter3(-99), fSigmaXYZDaughter3(-99), fSigmaZDaughter3(-99), fPtUncertDaughter3(-99), fTPCRefitDaughter3(-99), fITSRefitDaughter3(-99), fPropDCADaughter3(-1),
   //
   fEDaughter4(-99), fpDaughter4(-99), fptDaughter4(-99), fpxDaughter4(-99), fpyDaughter4(-99), fpzDaughter4(-99), fyDaughter4(-99), fDcaDaughter4(-99), fDcazDaughter4(-99), fDcaSecDaughter4(-99), fImParDaughter4(-99),  fImParzDaughter4(-99), fSigmaYXDaughter4(-99), fSigmaXYZDaughter4(-99), fSigmaZDaughter4(-99), fPtUncertDaughter4(-99), fPropDCADaughter4(-1),
   //
   fITSLayer1Daughter(-1), fITSLayer2Daughter(-1), fITSLayer3Daughter(-1), fITSLayer4Daughter(-1), fITSLayer5Daughter(-1), fITSLayer6Daughter(-1), fITSLayer1Daughter1(-1), fITSLayer2Daughter1(-1), fITSLayer3Daughter1(-1), fITSLayer4Daughter1(-1), fITSLayer5Daughter1(-1), fITSLayer6Daughter1(-1), fITSLayer1Daughter2(-1), fITSLayer2Daughter2(-1), fITSLayer3Daughter2(-1), fITSLayer4Daughter2(-1), fITSLayer5Daughter2(-1), fITSLayer6Daughter2(-1), fITSLayer1Daughter3(-1), fITSLayer2Daughter3(-1), fITSLayer3Daughter3(-1), fITSLayer4Daughter3(-1), fITSLayer5Daughter3(-1), fITSLayer6Daughter3(-1),
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
   fDcaDaughterXYKF(-99), fDcaDaughter1XYKF(-99), fDcaDaughter2XYKF(-99), fDcaDaughter3XYKF(-99), fDcaDaughter4XYKF(-99), fDcaSecDaughterXYKF(-99), fDcaSecDaughter1XYKF(-99), fDcaSecDaughter2XYKF(-99), fDcaSecDaughter3XYKF(-99), fDcaSecDaughter4XYKF(-99), fDcaDaughterZKF(-99), fDcaDaughter1ZKF(-99), fDcaDaughter2ZKF(-99), fDcaDaughter3ZKF(-99), fDcaDaughter4ZKF(-99), fDcaSecDaughterZKF(-99), fDcaSecDaughter1ZKF(-99), fDcaSecDaughter2ZKF(-99), fDcaSecDaughter3ZKF(-99), fDcaSecDaughter4ZKF(-99),
   //
   fDCA2BXYKF(-99), fDCA3B1XYKF(-99), fDCA3B2XYKF(-99), fDCA3B3XYKF(-99), fDCA2BZKF(-99), fDCA3B1ZKF(-99), fDCA3B2ZKF(-99), fDCA3B3ZKF(-99)
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
  DefineOutput(4, TTree::Class());
  DefineOutput(5, TTree::Class());
  DefineOutput(6, TTree::Class());
  DefineOutput(7, TTree::Class());
  DefineOutput(8, TTree::Class());
  DefineOutput(9, TTree::Class());
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
  if (!fInputHandler) {
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
  fHistdEdx = new TH2F("fHistdEdX", "dE/dx;#frac{#it{p}}{z} (GeV/#it{c});TPC Signal (a.u.)", 1000, -5.0, 5.0, 1000, 0.0, 2000);

  fHistNumEvents = new TH1F("fHistNumEvents", "Number of Events", 2, 0, 2);
  fHistNumEvents->GetXaxis()->SetBinLabel(1, "before PhysSel");
  fHistNumEvents->GetXaxis()->SetBinLabel(2, "after PhysSel");

  fHistTrigger = new TH1F("fHistTrigger", "Trigger", 8, 0, 8);
  fHistTrigger->GetXaxis()->SetBinLabel(1, "other");
  fHistTrigger->GetXaxis()->SetBinLabel(2, "kINT7");
  fHistTrigger->GetXaxis()->SetBinLabel(3, "kHighMultV0");
  fHistTrigger->GetXaxis()->SetBinLabel(4, "kHighMultSPD");
  fHistTrigger->GetXaxis()->SetBinLabel(5, "HNU");
  fHistTrigger->GetXaxis()->SetBinLabel(6, "HQU");
  fHistTrigger->GetXaxis()->SetBinLabel(7, "kCentral");
  fHistTrigger->GetXaxis()->SetBinLabel(8, "kSemiCentral");

  fHistogramList = new TList();
  fHistogramList->SetOwner(kTRUE);
  fHistogramList->SetName(GetName());
  fHistogramList->Add(fHistdEdx);
  fHistogramList->Add(fHistNumEvents);
  fHistogramList->Add(fHistTrigger);
  fEventCuts.AddQAplotsToList(fHistogramList);
  // _________________________________________________ //
  // __ associated Tree __ //
  fTree = new TTree("fTree", "fTree");
  fTree->Branch("fPeriod", &fPeriod, "fPeriod/I");
  fTree->Branch("frunnumber", &frunnumber, "frunnumber/I");
  fTree->Branch("feventclass", &feventclass, "feventclass/I");
  fTree->Branch("fCentrality", &fCentrality, "fCentrality/I");
  fTree->Branch("fMagneticField", &fMagneticField, "fMagneticField/I");
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
  fTree->Branch("fisOnlineV0_13", &fisOnlineV0_13, "fisOnlineV0_13/I");
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
  fTree->Branch("fV0VertexZ_24", &fV0VertexZ_24, "fV0VertexZ_24/D");
  fTree->Branch("fTrigMB", &fTrigMB, "fTrigMB/I");
  fTree->Branch("fTrigHMV0", &fTrigHMV0, "fTrigHMV0/I");
  fTree->Branch("fTrigHMSPD", &fTrigHMSPD, "fTrigHMSPD/I");
  fTree->Branch("fTrigHNU", &fTrigHNU, "fTrigHNU/I");
  fTree->Branch("fTrigHQU", &fTrigHQU, "fTrigHQU/I");
  fTree->Branch("fTrigkCentral", &fTrigkCentral, "fTrigkCentral/I");
  fTree->Branch("fTrigkSemiCentral", &fTrigkSemiCentral, "fTrigkSemiCentral/I");
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
  fTree->Branch("fPDGMother", &fPDGMother, "fPDGMother/I");
  fTree->Branch("fDecayChannel", &fDecayChannel, "fDecayChannel/I");
  fTree->Branch("fRecoMethod", &fRecoMethod, "fRecoMethod/I");
  fTree->Branch("fChargeMother", &fChargeMother, "fChargeMother/I");
  fTree->Branch("fmctruth", &fmctruth, "fmctruth/I");
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
  fTree->Branch("fPrimary4LHe", &fPrimary4LHe, "fPrimary4LHe/I");

  fTreeD = new TTree("fTreeD", "fTreeD");
  fTreeD->Branch("fPDGMother", &fPDGMother, "fPDGMother/I");
  fTreeD->Branch("fDecayChannel", &fDecayChannel, "fDecayChannel/I");
  fTreeD->Branch("fRecoMethod", &fRecoMethod, "fRecoMethod/I");
  fTreeD->Branch("fEDaughter", &fEDaughter, "fEDaughter/D");
  fTreeD->Branch("fpDaughter", &fpDaughter, "fpDaughter/D");
  fTreeD->Branch("fptDaughter", &fptDaughter, "fptDaughter/D");
  fTreeD->Branch("fpxDaughter", &fpxDaughter, "fpxDaughter/D");
  fTreeD->Branch("fpyDaughter", &fpyDaughter, "fpyDaughter/D");
  fTreeD->Branch("fpzDaughter", &fpzDaughter, "fpzDaughter/D");
  fTreeD->Branch("fyDaughter", &fyDaughter, "fyDaughter/D");
  fTreeD->Branch("fEDaughterKF", &fEDaughterKF, "fEDaughterKF/D");
  fTreeD->Branch("fLabelDaughterKF", &fLabelDaughterKF, "fLabelDaughterKF/D");
  fTreeD->Branch("fptDaughterKF", &fptDaughterKF, "fptDaughterKF/D");
  fTreeD->Branch("fpxDaughterKF", &fpxDaughterKF, "fpxDaughterKF/D");
  fTreeD->Branch("fpyDaughterKF", &fpyDaughterKF, "fpyDaughterKF/D");
  fTreeD->Branch("fpzDaughterKF", &fpzDaughterKF, "fpzDaughterKF/D");
  fTreeD->Branch("fyDaughterKF", &fyDaughterKF, "fyDaughterKF/D");
  fTreeD->Branch("fdEdxDaughter", &fdEdxDaughter, "fdEdxDaughter/D");
  fTreeD->Branch("fdEdxSigmaDaughter", &fdEdxSigmaDaughter, "fdEdxSigmaDaughter/D");
  fTreeD->Branch("fDcaDaughter", &fDcaDaughter, "fDcaDaughter/D");
  fTreeD->Branch("fDcaDaughtero", &fDcaDaughtero, "fDcaDaughtero/D");
  fTreeD->Branch("fDcazDaughter", &fDcazDaughter, "fDcazDaughter/D");
  fTreeD->Branch("fDcaSecDaughter", &fDcaSecDaughter, "fDcaSecDaughter/D");
  fTreeD->Branch("fImParDaughter", &fImParDaughter, "fImParDaughter/D");
  fTreeD->Branch("fImParzDaughter", &fImParzDaughter, "fImParzDaughter/D");
  fTreeD->Branch("fNclsDaughter", &fNclsDaughter, "fNclsDaughter/I");
  fTreeD->Branch("fChi2Daughter", &fChi2Daughter, "fChi2Daughter/D");
  fTreeD->Branch("fNclsITSDaughter", &fNclsITSDaughter, "fNclsITSDaughter/I");
  fTreeD->Branch("fEtaDaughter", &fEtaDaughter, "fEtaDaughter/D");
  fTreeD->Branch("fPhiDaughter", &fPhiDaughter, "fPhiDaughter/D");
  fTreeD->Branch("fGeoLengthDaughter", &fGeoLengthDaughter, "fGeoLengthDaughter/D");
  fTreeD->Branch("fTOFSignalDaughter", &fTOFSignalDaughter, "fTOFSignalDaughter/D");
  fTreeD->Branch("fSigmaYXDaughter", &fSigmaYXDaughter, "fSigmaYXDaughter/D");
  fTreeD->Branch("fSigmaXYZDaughter", &fSigmaXYZDaughter, "fSigmaXYZDaughter/D");
  fTreeD->Branch("fSigmaZDaughter", &fSigmaZDaughter, "fSigmaZDaughter/D");
  fTreeD->Branch("fPtUncertDaughter", &fPtUncertDaughter, "fPtUncertDaughter/D");
  fTreeD->Branch("fTPCRefitDaughter", &fTPCRefitDaughter, "fTPCRefitDaughter/I");
  fTreeD->Branch("fITSRefitDaughter", &fITSRefitDaughter, "fITSRefitDaughter/I");
  fTreeD->Branch("fPropDCADaughter", &fPropDCADaughter, "fPropDCADaughter/I");
  fTreeD->Branch("fITSLayer1Daughter", &fITSLayer1Daughter, "fITSLayer1Daughter/I");
  fTreeD->Branch("fITSLayer2Daughter", &fITSLayer2Daughter, "fITSLayer2Daughter/I");
  fTreeD->Branch("fITSLayer3Daughter", &fITSLayer3Daughter, "fITSLayer3Daughter/I");
  fTreeD->Branch("fITSLayer4Daughter", &fITSLayer4Daughter, "fITSLayer4Daughter/I");
  fTreeD->Branch("fITSLayer5Daughter", &fITSLayer5Daughter, "fITSLayer5Daughter/I");
  fTreeD->Branch("fITSLayer6Daughter", &fITSLayer6Daughter, "fITSLayer6Daughter/I");
  fTreeD->Branch("fEDaughter1", &fEDaughter1, "fEDaughter1/D");
  fTreeD->Branch("fpDaughter1", &fpDaughter1, "fpDaughter1/D");
  fTreeD->Branch("fptDaughter1", &fptDaughter1, "fptDaughter1/D");
  fTreeD->Branch("fpxDaughter1", &fpxDaughter1, "fpxDaughter1/D");
  fTreeD->Branch("fpyDaughter1", &fpyDaughter1, "fpyDaughter1/D");
  fTreeD->Branch("fpzDaughter1", &fpzDaughter1, "fpzDaughter1/D");
  fTreeD->Branch("fyDaughter1", &fyDaughter1, "fyDaughter1/D");
  fTreeD->Branch("fEDaughter1KF", &fEDaughter1KF, "fEDaughter1KF/D");
  fTreeD->Branch("fLabelDaughter1KF", &fLabelDaughter1KF, "fLabelDaughter1KF/D");
  fTreeD->Branch("fptDaughter1KF", &fptDaughter1KF, "fptDaughter1KF/D");
  fTreeD->Branch("fpxDaughter1KF", &fpxDaughter1KF, "fpxDaughter1KF/D");
  fTreeD->Branch("fpyDaughter1KF", &fpyDaughter1KF, "fpyDaughter1KF/D");
  fTreeD->Branch("fpzDaughter1KF", &fpzDaughter1KF, "fpzDaughter1KF/D");
  fTreeD->Branch("fyDaughter1KF", &fyDaughter1KF, "fyDaughter1KF/D");
  fTreeD->Branch("fdEdxDaughter1", &fdEdxDaughter1, "fdEdxDaughter1/D");
  fTreeD->Branch("fdEdxSigmaDaughter1", &fdEdxSigmaDaughter1, "fdEdxSigmaDaughter1/D");
  fTreeD->Branch("fDcaDaughter1", &fDcaDaughter1, "fDcaDaughter1/D");
  fTreeD->Branch("fDcaDaughter1o", &fDcaDaughter1o, "fDcaDaughter1o/D");
  fTreeD->Branch("fDcazDaughter1", &fDcazDaughter1, "fDcazDaughter1/D");
  fTreeD->Branch("fDcaSecDaughter1", &fDcaSecDaughter1, "fDcaSecDaughter1/D");
  fTreeD->Branch("fImParDaughter1", &fImParDaughter1, "fImParDaughter1/D");
  fTreeD->Branch("fImParzDaughter1", &fImParzDaughter1, "fImParzDaughter1/D");
  fTreeD->Branch("fNclsDaughter1", &fNclsDaughter1, "fNclsDaughter1/I");
  fTreeD->Branch("fChi2Daughter1", &fChi2Daughter1, "fChi2Daughter1/D");
  fTreeD->Branch("fNclsITSDaughter1", &fNclsITSDaughter1, "fNclsITSDaughter1/I");
  fTreeD->Branch("fEtaDaughter1", &fEtaDaughter1, "fEtaDaughter1/D");
  fTreeD->Branch("fPhiDaughter1", &fPhiDaughter1, "fPhiDaughter1/D");
  fTreeD->Branch("fGeoLengthDaughter1", &fGeoLengthDaughter1, "fGeoLengthDaughter1/D");
  fTreeD->Branch("fTOFSignalDaughter1", &fTOFSignalDaughter1, "fTOFSignalDaughter1/D");
  fTreeD->Branch("fSigmaYXDaughter1", &fSigmaYXDaughter1, "fSigmaYXDaughter1/D");
  fTreeD->Branch("fSigmaXYZDaughter1", &fSigmaXYZDaughter1, "fSigmaXYZDaughter1/D");
  fTreeD->Branch("fSigmaZDaughter1", &fSigmaZDaughter1, "fSigmaZDaughter1/D");
  fTreeD->Branch("fPtUncertDaughter1", &fPtUncertDaughter1, "fPtUncertDaughter1/D");
  fTreeD->Branch("fTPCRefitDaughter1", &fTPCRefitDaughter1, "fTPCRefitDaughter1/I");
  fTreeD->Branch("fITSRefitDaughter1", &fITSRefitDaughter1, "fITSRefitDaughter1/I");
  fTreeD->Branch("fPropDCADaughter1", &fPropDCADaughter1, "fPropDCADaughter1/I");
  fTreeD->Branch("fITSLayer1Daughter1", &fITSLayer1Daughter1, "fITSLayer1Daughter1/I");
  fTreeD->Branch("fITSLayer2Daughter1", &fITSLayer2Daughter1, "fITSLayer2Daughter1/I");
  fTreeD->Branch("fITSLayer3Daughter1", &fITSLayer3Daughter1, "fITSLayer3Daughter1/I");
  fTreeD->Branch("fITSLayer4Daughter1", &fITSLayer4Daughter1, "fITSLayer4Daughter1/I");
  fTreeD->Branch("fITSLayer5Daughter1", &fITSLayer5Daughter1, "fITSLayer5Daughter1/I");
  fTreeD->Branch("fITSLayer6Daughter1", &fITSLayer6Daughter1, "fITSLayer6Daughter1/I");
  fTreeD->Branch("fEDaughter2", &fEDaughter2, "fEDaughter2/D");
  fTreeD->Branch("fpDaughter2", &fpDaughter2, "fpDaughter2/D");
  fTreeD->Branch("fptDaughter2", &fptDaughter2, "fptDaughter2/D");
  fTreeD->Branch("fpxDaughter2", &fpxDaughter2, "fpxDaughter2/D");
  fTreeD->Branch("fpyDaughter2", &fpyDaughter2, "fpyDaughter2/D");
  fTreeD->Branch("fpzDaughter2", &fpzDaughter2, "fpzDaughter2/D");
  fTreeD->Branch("fyDaughter2", &fyDaughter2, "fyDaughter2/D");
  fTreeD->Branch("fEDaughter2KF", &fEDaughter2KF, "fEDaughter2KF/D");
  fTreeD->Branch("fLabelDaughter2KF", &fLabelDaughter2KF, "fLabelDaughter2KF/D");
  fTreeD->Branch("fptDaughter2KF", &fptDaughter2KF, "fptDaughter2KF/D");
  fTreeD->Branch("fpxDaughter2KF", &fpxDaughter2KF, "fpxDaughter2KF/D");
  fTreeD->Branch("fpyDaughter2KF", &fpyDaughter2KF, "fpyDaughter2KF/D");
  fTreeD->Branch("fpzDaughter2KF", &fpzDaughter2KF, "fpzDaughter2KF/D");
  fTreeD->Branch("fyDaughter2KF", &fyDaughter2KF, "fyDaughter2KF/D");
  fTreeD->Branch("fdEdxDaughter2", &fdEdxDaughter2, "fdEdxDaughter2/D");
  fTreeD->Branch("fdEdxSigmaDaughter2", &fdEdxSigmaDaughter2, "fdEdxSigmaDaughter2/D");
  fTreeD->Branch("fDcaDaughter2", &fDcaDaughter2, "fDcaDaughter2/D");
  fTreeD->Branch("fDcaDaughter2o", &fDcaDaughter2o, "fDcaDaughter2o/D");
  fTreeD->Branch("fDcazDaughter2", &fDcazDaughter2, "fDcazDaughter2/D");
  fTreeD->Branch("fDcaSecDaughter2", &fDcaSecDaughter2, "fDcaSecDaughter2/D");
  fTreeD->Branch("fImParDaughter2", &fImParDaughter2, "fImParDaughter2/D");
  fTreeD->Branch("fImParzDaughter2", &fImParzDaughter2, "fImParzDaughter2/D");
  fTreeD->Branch("fNclsDaughter2", &fNclsDaughter2, "fNclsDaughter2/I");
  fTreeD->Branch("fChi2Daughter2", &fChi2Daughter2, "fChi2Daughter2/D");
  fTreeD->Branch("fNclsITSDaughter2", &fNclsITSDaughter2, "fNclsITSDaughter2/I");
  fTreeD->Branch("fEtaDaughter2", &fEtaDaughter2, "fEtaDaughter2/D");
  fTreeD->Branch("fPhiDaughter2", &fPhiDaughter2, "fPhiDaughter2/D");
  fTreeD->Branch("fGeoLengthDaughter2", &fGeoLengthDaughter2, "fGeoLengthDaughter2/D");
  fTreeD->Branch("fTOFSignalDaughter2", &fTOFSignalDaughter2, "fTOFSignalDaughter2/D");
  fTreeD->Branch("fSigmaYXDaughter2", &fSigmaYXDaughter2, "fSigmaYXDaughter2/D");
  fTreeD->Branch("fSigmaXYZDaughter2", &fSigmaXYZDaughter2, "fSigmaXYZDaughter2/D");
  fTreeD->Branch("fSigmaZDaughter2", &fSigmaZDaughter2, "fSigmaZDaughter2/D");
  fTreeD->Branch("fPtUncertDaughter2", &fPtUncertDaughter2, "fPtUncertDaughter2/D");
  fTreeD->Branch("fTPCRefitDaughter2", &fTPCRefitDaughter2, "fTPCRefitDaughter2/I");
  fTreeD->Branch("fITSRefitDaughter2", &fITSRefitDaughter2, "fITSRefitDaughter2/I");
  fTreeD->Branch("fPropDCADaughter2", &fPropDCADaughter2, "fPropDCADaughter2/I");
  fTreeD->Branch("fITSLayer1Daughter2", &fITSLayer1Daughter2, "fITSLayer1Daughter2/I");
  fTreeD->Branch("fITSLayer2Daughter2", &fITSLayer2Daughter2, "fITSLayer2Daughter2/I");
  fTreeD->Branch("fITSLayer3Daughter2", &fITSLayer3Daughter2, "fITSLayer3Daughter2/I");
  fTreeD->Branch("fITSLayer4Daughter2", &fITSLayer4Daughter2, "fITSLayer4Daughter2/I");
  fTreeD->Branch("fITSLayer5Daughter2", &fITSLayer5Daughter2, "fITSLayer5Daughter2/I");
  fTreeD->Branch("fITSLayer6Daughter2", &fITSLayer6Daughter2, "fITSLayer6Daughter2/I");
  fTreeD->Branch("fEDaughter3", &fEDaughter3, "fEDaughter3/D");
  fTreeD->Branch("fpDaughter3", &fpDaughter3, "fpDaughter3/D");
  fTreeD->Branch("fptDaughter3", &fptDaughter3, "fptDaughter3/D");
  fTreeD->Branch("fpxDaughter3", &fpxDaughter3, "fpxDaughter3/D");
  fTreeD->Branch("fpyDaughter3", &fpyDaughter3, "fpyDaughter3/D");
  fTreeD->Branch("fpzDaughter3", &fpzDaughter3, "fpzDaughter3/D");
  fTreeD->Branch("fyDaughter3", &fyDaughter3, "fyDaughter3/D");
  fTreeD->Branch("fEDaughter3KF", &fEDaughter3KF, "fEDaughter3KF/D");
  fTreeD->Branch("fLabelDaughter3KF", &fLabelDaughter3KF, "fLabelDaughter3KF/D");
  fTreeD->Branch("fptDaughter3KF", &fptDaughter3KF, "fptDaughter3KF/D");
  fTreeD->Branch("fpxDaughter3KF", &fpxDaughter3KF, "fpxDaughter3KF/D");
  fTreeD->Branch("fpyDaughter3KF", &fpyDaughter3KF, "fpyDaughter3KF/D");
  fTreeD->Branch("fpzDaughter3KF", &fpzDaughter3KF, "fpzDaughter3KF/D");
  fTreeD->Branch("fyDaughter3KF", &fyDaughter3KF, "fyDaughter3KF/D");
  fTreeD->Branch("fdEdxDaughter3", &fdEdxDaughter3, "fdEdxDaughter3/D");
  fTreeD->Branch("fdEdxSigmaDaughter3", &fdEdxSigmaDaughter3, "fdEdxSigmaDaughter3/D");
  fTreeD->Branch("fDcaDaughter3", &fDcaDaughter3, "fDcaDaughter3/D");
  fTreeD->Branch("fDcaDaughter3o", &fDcaDaughter3o, "fDcaDaughter3o/D");
  fTreeD->Branch("fDcazDaughter3", &fDcazDaughter3, "fDcazDaughter3/D");
  fTreeD->Branch("fDcaSecDaughter3", &fDcaSecDaughter3, "fDcaSecDaughter3/D");
  fTreeD->Branch("fImParDaughter3", &fImParDaughter3, "fImParDaughter3/D");
  fTreeD->Branch("fImParzDaughter3", &fImParzDaughter3, "fImParzDaughter3/D");
  fTreeD->Branch("fNclsDaughter3", &fNclsDaughter3, "fNclsDaughter3/I");
  fTreeD->Branch("fChi2Daughter3", &fChi2Daughter3, "fChi2Daughter3/D");
  fTreeD->Branch("fNclsITSDaughter3", &fNclsITSDaughter3, "fNclsITSDaughter3/I");
  fTreeD->Branch("fEtaDaughter3", &fEtaDaughter3, "fEtaDaughter3/D");
  fTreeD->Branch("fPhiDaughter3", &fPhiDaughter3, "fPhiDaughter3/D");
  fTreeD->Branch("fGeoLengthDaughter3", &fGeoLengthDaughter3, "fGeoLengthDaughter3/D");
  fTreeD->Branch("fTOFSignalDaughter3", &fTOFSignalDaughter3, "fTOFSignalDaughter3/D");
  fTreeD->Branch("fSigmaYXDaughter3", &fSigmaYXDaughter3, "fSigmaYXDaughter3/D");
  fTreeD->Branch("fSigmaXYZDaughter3", &fSigmaXYZDaughter3, "fSigmaXYZDaughter3/D");
  fTreeD->Branch("fSigmaZDaughter3", &fSigmaZDaughter3, "fSigmaZDaughter3/D");
  fTreeD->Branch("fPtUncertDaughter3", &fPtUncertDaughter3, "fPtUncertDaughter3/D");
  fTreeD->Branch("fTPCRefitDaughter3", &fTPCRefitDaughter3, "fTPCRefitDaughter3/I");
  fTreeD->Branch("fITSRefitDaughter3", &fITSRefitDaughter3, "fITSRefitDaughter3/I");
  fTreeD->Branch("fPropDCADaughter3", &fPropDCADaughter3, "fPropDCADaughter3/I");
  fTreeD->Branch("fITSLayer1Daughter3", &fITSLayer1Daughter3, "fITSLayer1Daughter3/I");
  fTreeD->Branch("fITSLayer2Daughter3", &fITSLayer2Daughter3, "fITSLayer2Daughter3/I");
  fTreeD->Branch("fITSLayer3Daughter3", &fITSLayer3Daughter3, "fITSLayer3Daughter3/I");
  fTreeD->Branch("fITSLayer4Daughter3", &fITSLayer4Daughter3, "fITSLayer4Daughter3/I");
  fTreeD->Branch("fITSLayer5Daughter3", &fITSLayer5Daughter3, "fITSLayer5Daughter3/I");
  fTreeD->Branch("fITSLayer6Daughter3", &fITSLayer6Daughter3, "fITSLayer6Daughter3/I");
  fTreeD->Branch("fEDaughter4", &fEDaughter4, "fEDaughter4/D");
  fTreeD->Branch("fpDaughter4", &fpDaughter4, "fpDaughter4/D");
  fTreeD->Branch("fptDaughter4", &fptDaughter4, "fptDaughter4/D");
  fTreeD->Branch("fpxDaughter4", &fpxDaughter4, "fpxDaughter4/D");
  fTreeD->Branch("fpyDaughter4", &fpyDaughter4, "fpyDaughter4/D");
  fTreeD->Branch("fpzDaughter4", &fpzDaughter4, "fpzDaughter4/D");
  fTreeD->Branch("fyDaughter4", &fyDaughter4, "fyDaughter4/D");
  fTreeD->Branch("fDcaDaughter4", &fDcaDaughter4, "fDcaDaughter4/D");
  fTreeD->Branch("fDcazDaughter4", &fDcazDaughter4, "fDcazDaughter4/D");
  fTreeD->Branch("fDcaSecDaughter4", &fDcaSecDaughter4, "fDcaSecDaughter4/D");
  fTreeD->Branch("fImParDaughter4", &fImParDaughter4, "fImParDaughter4/D");
  fTreeD->Branch("fImParzDaughter4", &fImParzDaughter4, "fImParzDaughter4/D");
  fTreeD->Branch("fSigmaYXDaughter4", &fSigmaYXDaughter4, "fSigmaYXDaughter4/D");
  fTreeD->Branch("fSigmaXYZDaughter4", &fSigmaXYZDaughter4, "fSigmaXYZDaughter4/D");
  fTreeD->Branch("fSigmaZDaughter4", &fSigmaZDaughter4, "fSigmaZDaughter4/D");
  fTreeD->Branch("fPtUncertDaughter4", &fPtUncertDaughter4, "fPtUncertDaughter4/D");
  fTreeD->Branch("fPropDCADaughter4", &fPropDCADaughter4, "fPropDCADaughter4/I");
  fTreeD->Branch("fdEdxSigmaPion", &fdEdxSigmaPion, "fdEdxSigmaPion/D");
  fTreeD->Branch("fdEdxSigmaDeuteron", &fdEdxSigmaDeuteron, "fdEdxSigmaDeuteron/D");
  fTreeD->Branch("fdEdxSigmaTriton", &fdEdxSigmaTriton, "fdEdxSigmaTriton/D");
  fTreeD->Branch("fdEdxSigmaAlpha", &fdEdxSigmaAlpha, "fdEdxSigmaAlpha/D");

  fTreeKF = new TTree("fTreeKF", "fTreeKF");
  fTreeKF->Branch("fPeriod", &fPeriod, "fPeriod/I");
  fTreeKF->Branch("frunnumber", &frunnumber, "frunnumber/I");
  fTreeKF->Branch("feventclass", &feventclass, "feventclass/I");
  fTreeKF->Branch("fCentrality", &fCentrality, "fCentrality/I");  
  fTreeKF->Branch("fMagneticField", &fMagneticField, "fMagneticField/I");
  fTreeKF->Branch("fTrigMB", &fTrigMB, "fTrigMB/I");
  fTreeKF->Branch("fTrigHMV0", &fTrigHMV0, "fTrigHMV0/I");
  fTreeKF->Branch("fTrigHMSPD", &fTrigHMSPD, "fTrigHMSPD/I");
  fTreeKF->Branch("fTrigHNU", &fTrigHNU, "fTrigHNU/I");
  fTreeKF->Branch("fTrigHQU", &fTrigHQU, "fTrigHQU/I");
  fTreeKF->Branch("fTrigkCentral", &fTrigkCentral, "fTrigkCentral/I");
  fTreeKF->Branch("fTrigkSemiCentral", &fTrigkSemiCentral, "fTrigkSemiCentral/I");
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
  fTreeKF->Branch("fisOnlineV0_13", &fisOnlineV0_13, "fisOnlineV0_13/I");
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
  fTreeKF->Branch("fV0VertexZ_24", &fV0VertexZ_24, "fV0VertexZ_24/D");
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
  fTreeKF->Branch("fPDGMother", &fPDGMother, "fPDGMother/I");
  fTreeKF->Branch("fDecayChannel", &fDecayChannel, "fDecayChannel/I");
  fTreeKF->Branch("fRecoMethod", &fRecoMethod, "fRecoMethod/I");
  fTreeKF->Branch("fChargeMother", &fChargeMother, "fChargeMother/I");
  fTreeKF->Branch("fmctruth", &fmctruth, "fmctruth/I");
  fTreeKF->Branch("fPrimary4LHe", &fPrimary4LHe, "fPrimary4LHe/I");
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
  // _________________________________________________ //
  // __ associated 4LHe tree __ //
  gTree = new TTree("gTree", "gTree");
  gTree->Branch("fPeriod", &fPeriod, "fPeriod/I");
  gTree->Branch("frunnumber", &frunnumber, "frunnumber/I");
  gTree->Branch("feventclass", &feventclass, "feventclass/I");
  gTree->Branch("fCentrality", &fCentrality, "fCentrality/I");
  gTree->Branch("fMagneticField", &fMagneticField, "fMagneticField/I"); 
  gTree->Branch("fTrigMB", &fTrigMB, "fTrigMB/I");
  gTree->Branch("fTrigHMV0", &fTrigHMV0, "fTrigHMV0/I");
  gTree->Branch("fTrigHMSPD", &fTrigHMSPD, "fTrigHMSPD/I");
  gTree->Branch("fTrigHNU", &fTrigHNU, "fTrigHNU/I");
  gTree->Branch("fTrigHQU", &fTrigHQU, "fTrigHQU/I");
  gTree->Branch("fTrigkCentral", &fTrigkCentral, "fTrigkCentral/I");
  gTree->Branch("fTrigkSemiCentral", &fTrigkSemiCentral, "fTrigkSemiCentral/I");
  gTree->Branch("fPDGMother", &fPDGMother, "fPDGMother/I");
  gTree->Branch("fChargeMother", &fChargeMother, "fChargeMother/I");
  gTree->Branch("fmctruth", &fmctruth, "fmctruth/I");
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

  gTreeD = new TTree("gTreeD", "gTreeD");
  gTreeD->Branch("fPDGMother", &fPDGMother, "fPDGMother/I");
  gTreeD->Branch("fDecayChannel", &fDecayChannel, "fDecayChannel/I");
  gTreeD->Branch("fRecoMethod", &fRecoMethod, "fRecoMethod/I");
  gTreeD->Branch("fEDaughter", &fEDaughter, "fEDaughter/D");
  gTreeD->Branch("fpDaughter", &fpDaughter, "fpDaughter/D");
  gTreeD->Branch("fptDaughter", &fptDaughter, "fptDaughter/D");
  gTreeD->Branch("fpxDaughter", &fpxDaughter, "fpxDaughter/D");
  gTreeD->Branch("fpyDaughter", &fpyDaughter, "fpyDaughter/D");
  gTreeD->Branch("fpzDaughter", &fpzDaughter, "fpzDaughter/D");
  gTreeD->Branch("fyDaughter", &fyDaughter, "fyDaughter/D");
  gTreeD->Branch("fEDaughterKF", &fEDaughterKF, "fEDaughterKF/D");
  gTreeD->Branch("fLabelDaughterKF", &fLabelDaughterKF, "fLabelDaughterKF/D");
  gTreeD->Branch("fptDaughterKF", &fptDaughterKF, "fptDaughterKF/D");
  gTreeD->Branch("fpxDaughterKF", &fpxDaughterKF, "fpxDaughterKF/D");
  gTreeD->Branch("fpyDaughterKF", &fpyDaughterKF, "fpyDaughterKF/D");
  gTreeD->Branch("fpzDaughterKF", &fpzDaughterKF, "fpzDaughterKF/D");
  gTreeD->Branch("fyDaughterKF", &fyDaughterKF, "fyDaughterKF/D");
  gTreeD->Branch("fdEdxDaughter", &fdEdxDaughter, "fdEdxDaughter/D");
  gTreeD->Branch("fdEdxSigmaDaughter", &fdEdxSigmaDaughter, "fdEdxSigmaDaughter/D");
  gTreeD->Branch("fDcaDaughter", &fDcaDaughter, "fDcaDaughter/D");
  gTreeD->Branch("fDcaDaughtero", &fDcaDaughtero, "fDcaDaughtero/D");
  gTreeD->Branch("fDcazDaughter", &fDcazDaughter, "fDcazDaughter/D");
  gTreeD->Branch("fDcaSecDaughter", &fDcaSecDaughter, "fDcaSecDaughter/D");
  gTreeD->Branch("fImParDaughter", &fImParDaughter, "fImParDaughter/D");
  gTreeD->Branch("fImParzDaughter", &fImParzDaughter, "fImParzDaughter/D");
  gTreeD->Branch("fNclsDaughter", &fNclsDaughter, "fNclsDaughter/I");
  gTreeD->Branch("fChi2Daughter", &fChi2Daughter, "fChi2Daughter/D");
  gTreeD->Branch("fNclsITSDaughter", &fNclsITSDaughter, "fNclsITSDaughter/I");
  gTreeD->Branch("fEtaDaughter", &fEtaDaughter, "fEtaDaughter/D");
  gTreeD->Branch("fPhiDaughter", &fPhiDaughter, "fPhiDaughter/D");
  gTreeD->Branch("fGeoLengthDaughter", &fGeoLengthDaughter, "fGeoLengthDaughter/D");
  gTreeD->Branch("fTOFSignalDaughter", &fTOFSignalDaughter, "fTOFSignalDaughter/D");
  gTreeD->Branch("fSigmaYXDaughter", &fSigmaYXDaughter, "fSigmaYXDaughter/D");
  gTreeD->Branch("fSigmaXYZDaughter", &fSigmaXYZDaughter, "fSigmaXYZDaughter/D");
  gTreeD->Branch("fSigmaZDaughter", &fSigmaZDaughter, "fSigmaZDaughter/D");
  gTreeD->Branch("fPtUncertDaughter", &fPtUncertDaughter, "fPtUncertDaughter/D");
  gTreeD->Branch("fTPCRefitDaughter", &fTPCRefitDaughter, "fTPCRefitDaughter/I");
  gTreeD->Branch("fITSRefitDaughter", &fITSRefitDaughter, "fITSRefitDaughter/I");
  gTreeD->Branch("fPropDCADaughter", &fPropDCADaughter, "fPropDCADaughter/I");
  gTreeD->Branch("fITSLayer1Daughter", &fITSLayer1Daughter, "fITSLayer1Daughter/I");
  gTreeD->Branch("fITSLayer2Daughter", &fITSLayer2Daughter, "fITSLayer2Daughter/I");
  gTreeD->Branch("fITSLayer3Daughter", &fITSLayer3Daughter, "fITSLayer3Daughter/I");
  gTreeD->Branch("fITSLayer4Daughter", &fITSLayer4Daughter, "fITSLayer4Daughter/I");
  gTreeD->Branch("fITSLayer5Daughter", &fITSLayer5Daughter, "fITSLayer5Daughter/I");
  gTreeD->Branch("fITSLayer6Daughter", &fITSLayer6Daughter, "fITSLayer6Daughter/I");
  gTreeD->Branch("fEDaughter1", &fEDaughter1, "fEDaughter1/D");
  gTreeD->Branch("fpDaughter1", &fpDaughter1, "fpDaughter1/D");
  gTreeD->Branch("fptDaughter1", &fptDaughter1, "fptDaughter1/D");
  gTreeD->Branch("fpxDaughter1", &fpxDaughter1, "fpxDaughter1/D");
  gTreeD->Branch("fpyDaughter1", &fpyDaughter1, "fpyDaughter1/D");
  gTreeD->Branch("fpzDaughter1", &fpzDaughter1, "fpzDaughter1/D");
  gTreeD->Branch("fyDaughter1", &fyDaughter1, "fyDaughter1/D");
  gTreeD->Branch("fEDaughter1KF", &fEDaughter1KF, "fEDaughter1KF/D");
  gTreeD->Branch("fLabelDaughter1KF", &fLabelDaughter1KF, "fLabelDaughter1KF/D");
  gTreeD->Branch("fptDaughter1KF", &fptDaughter1KF, "fptDaughter1KF/D");
  gTreeD->Branch("fpxDaughter1KF", &fpxDaughter1KF, "fpxDaughter1KF/D");
  gTreeD->Branch("fpyDaughter1KF", &fpyDaughter1KF, "fpyDaughter1KF/D");
  gTreeD->Branch("fpzDaughter1KF", &fpzDaughter1KF, "fpzDaughter1KF/D");
  gTreeD->Branch("fyDaughter1KF", &fyDaughter1KF, "fyDaughter1KF/D");
  gTreeD->Branch("fdEdxDaughter1", &fdEdxDaughter1, "fdEdxDaughter1/D");
  gTreeD->Branch("fdEdxSigmaDaughter1", &fdEdxSigmaDaughter1, "fdEdxSigmaDaughter1/D");
  gTreeD->Branch("fDcaDaughter1", &fDcaDaughter1, "fDcaDaughter1/D");
  gTreeD->Branch("fDcaDaughter1o", &fDcaDaughter1o, "fDcaDaughter1o/D");
  gTreeD->Branch("fDcazDaughter1", &fDcazDaughter1, "fDcazDaughter1/D");
  gTreeD->Branch("fDcaSecDaughter1", &fDcaSecDaughter1, "fDcaSecDaughter1/D");
  gTreeD->Branch("fImParDaughter1", &fImParDaughter1, "fImParDaughter1/D");
  gTreeD->Branch("fImParzDaughter1", &fImParzDaughter1, "fImParzDaughter1/D");
  gTreeD->Branch("fNclsDaughter1", &fNclsDaughter1, "fNclsDaughter1/I");
  gTreeD->Branch("fChi2Daughter1", &fChi2Daughter1, "fChi2Daughter1/D");
  gTreeD->Branch("fNclsITSDaughter1", &fNclsITSDaughter1, "fNclsITSDaughter1/I");
  gTreeD->Branch("fEtaDaughter1", &fEtaDaughter1, "fEtaDaughter1/D");
  gTreeD->Branch("fPhiDaughter1", &fPhiDaughter1, "fPhiDaughter1/D");
  gTreeD->Branch("fGeoLengthDaughter1", &fGeoLengthDaughter1, "fGeoLengthDaughter1/D");
  gTreeD->Branch("fTOFSignalDaughter1", &fTOFSignalDaughter1, "fTOFSignalDaughter1/D");
  gTreeD->Branch("fSigmaYXDaughter1", &fSigmaYXDaughter1, "fSigmaYXDaughter1/D");
  gTreeD->Branch("fSigmaXYZDaughter1", &fSigmaXYZDaughter1, "fSigmaXYZDaughter1/D");
  gTreeD->Branch("fSigmaZDaughter1", &fSigmaZDaughter1, "fSigmaZDaughter1/D");
  gTreeD->Branch("fPtUncertDaughter1", &fPtUncertDaughter1, "fPtUncertDaughter1/D");
  gTreeD->Branch("fTPCRefitDaughter1", &fTPCRefitDaughter1, "fTPCRefitDaughter1/I");
  gTreeD->Branch("fITSRefitDaughter1", &fITSRefitDaughter1, "fITSRefitDaughter1/I");
  gTreeD->Branch("fPropDCADaughter1", &fPropDCADaughter1, "fPropDCADaughter1/I");
  gTreeD->Branch("fITSLayer1Daughter1", &fITSLayer1Daughter1, "fITSLayer1Daughter1/I");
  gTreeD->Branch("fITSLayer2Daughter1", &fITSLayer2Daughter1, "fITSLayer2Daughter1/I");
  gTreeD->Branch("fITSLayer3Daughter1", &fITSLayer3Daughter1, "fITSLayer3Daughter1/I");
  gTreeD->Branch("fITSLayer4Daughter1", &fITSLayer4Daughter1, "fITSLayer4Daughter1/I");
  gTreeD->Branch("fITSLayer5Daughter1", &fITSLayer5Daughter1, "fITSLayer5Daughter1/I");
  gTreeD->Branch("fITSLayer6Daughter1", &fITSLayer6Daughter1, "fITSLayer6Daughter1/I");
  gTreeD->Branch("fEDaughter2", &fEDaughter2, "fEDaughter2/D");
  gTreeD->Branch("fpDaughter2", &fpDaughter2, "fpDaughter2/D");
  gTreeD->Branch("fptDaughter2", &fptDaughter2, "fptDaughter2/D");
  gTreeD->Branch("fpxDaughter2", &fpxDaughter2, "fpxDaughter2/D");
  gTreeD->Branch("fpyDaughter2", &fpyDaughter2, "fpyDaughter2/D");
  gTreeD->Branch("fpzDaughter2", &fpzDaughter2, "fpzDaughter2/D");
  gTreeD->Branch("fyDaughter2", &fyDaughter2, "fyDaughter2/D");
  gTreeD->Branch("fEDaughter2KF", &fEDaughter2KF, "fEDaughter2KF/D");
  gTreeD->Branch("fLabelDaughter2KF", &fLabelDaughter2KF, "fLabelDaughter2KF/D");
  gTreeD->Branch("fptDaughter2KF", &fptDaughter2KF, "fptDaughter2KF/D");
  gTreeD->Branch("fpxDaughter2KF", &fpxDaughter2KF, "fpxDaughter2KF/D");
  gTreeD->Branch("fpyDaughter2KF", &fpyDaughter2KF, "fpyDaughter2KF/D");
  gTreeD->Branch("fpzDaughter2KF", &fpzDaughter2KF, "fpzDaughter2KF/D");
  gTreeD->Branch("fyDaughter2KF", &fyDaughter2KF, "fyDaughter2KF/D");
  gTreeD->Branch("fdEdxDaughter2", &fdEdxDaughter2, "fdEdxDaughter2/D");
  gTreeD->Branch("fdEdxSigmaDaughter2", &fdEdxSigmaDaughter2, "fdEdxSigmaDaughter2/D");
  gTreeD->Branch("fDcaDaughter2", &fDcaDaughter2, "fDcaDaughter2/D");
  gTreeD->Branch("fDcaDaughter2o", &fDcaDaughter2o, "fDcaDaughter2o/D");
  gTreeD->Branch("fDcazDaughter2", &fDcazDaughter2, "fDcazDaughter2/D");
  gTreeD->Branch("fDcaSecDaughter2", &fDcaSecDaughter2, "fDcaSecDaughter2/D");
  gTreeD->Branch("fImParDaughter2", &fImParDaughter2, "fImParDaughter2/D");
  gTreeD->Branch("fImParzDaughter2", &fImParzDaughter2, "fImParzDaughter2/D");
  gTreeD->Branch("fNclsDaughter2", &fNclsDaughter2, "fNclsDaughter2/I");
  gTreeD->Branch("fChi2Daughter2", &fChi2Daughter2, "fChi2Daughter2/D");
  gTreeD->Branch("fNclsITSDaughter2", &fNclsITSDaughter2, "fNclsITSDaughter2/I");
  gTreeD->Branch("fEtaDaughter2", &fEtaDaughter2, "fEtaDaughter2/D");
  gTreeD->Branch("fPhiDaughter2", &fPhiDaughter2, "fPhiDaughter2/D");
  gTreeD->Branch("fGeoLengthDaughter2", &fGeoLengthDaughter2, "fGeoLengthDaughter2/D");
  gTreeD->Branch("fTOFSignalDaughter2", &fTOFSignalDaughter2, "fTOFSignalDaughter2/D");
  gTreeD->Branch("fSigmaYXDaughter2", &fSigmaYXDaughter2, "fSigmaYXDaughter2/D");
  gTreeD->Branch("fSigmaXYZDaughter2", &fSigmaXYZDaughter2, "fSigmaXYZDaughter2/D");
  gTreeD->Branch("fSigmaZDaughter2", &fSigmaZDaughter2, "fSigmaZDaughter2/D");
  gTreeD->Branch("fPtUncertDaughter2", &fPtUncertDaughter2, "fPtUncertDaughter2/D");
  gTreeD->Branch("fTPCRefitDaughter2", &fTPCRefitDaughter2, "fTPCRefitDaughter2/I");
  gTreeD->Branch("fITSRefitDaughter2", &fITSRefitDaughter2, "fITSRefitDaughter2/I");
  gTreeD->Branch("fPropDCADaughter2", &fPropDCADaughter2, "fPropDCADaughter2/I");
  gTreeD->Branch("fITSLayer1Daughter2", &fITSLayer1Daughter2, "fITSLayer1Daughter2/I");
  gTreeD->Branch("fITSLayer2Daughter2", &fITSLayer2Daughter2, "fITSLayer2Daughter2/I");
  gTreeD->Branch("fITSLayer3Daughter2", &fITSLayer3Daughter2, "fITSLayer3Daughter2/I");
  gTreeD->Branch("fITSLayer4Daughter2", &fITSLayer4Daughter2, "fITSLayer4Daughter2/I");
  gTreeD->Branch("fITSLayer5Daughter2", &fITSLayer5Daughter2, "fITSLayer5Daughter2/I");
  gTreeD->Branch("fITSLayer6Daughter2", &fITSLayer6Daughter2, "fITSLayer6Daughter2/I");
  gTreeD->Branch("fdEdxSigmaPion", &fdEdxSigmaPion, "fdEdxSigmaPion/D");
  gTreeD->Branch("fdEdxSigmaDeuteron", &fdEdxSigmaDeuteron, "fdEdxSigmaDeuteron/D");
  gTreeD->Branch("fdEdxSigmaTriton", &fdEdxSigmaTriton, "fdEdxSigmaTriton/D");
  gTreeD->Branch("fdEdxSigmaAlpha", &fdEdxSigmaAlpha, "fdEdxSigmaAlpha/D");

  gTreeKF = new TTree("gTreeKF", "gTreeKF");
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
  // _________________________________________________ //
  // __ generated tree __ //
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
  fTreeGen->Branch("fPDGMother", &fPDGMother, "fPDGMother/I");
  fTreeGen->Branch("fChargeMother", &fChargeMother, "fChargeMother/I");
  fTreeGen->Branch("fptDaughter", &fptDaughter, "fptDaughter/D");
  fTreeGen->Branch("fyDaughter", &fyDaughter, "fyDaughter/D");
  fTreeGen->Branch("fptDaughter1", &fptDaughter1, "fptDaughter1/D");
  fTreeGen->Branch("fyDaughter1", &fyDaughter1, "fyDaughter1/D");
  fTreeGen->Branch("fptDaughter2", &fptDaughter2, "fptDaughter2/D");
  fTreeGen->Branch("fyDaughter2", &fyDaughter2, "fyDaughter2/D");
  fTreeGen->Branch("fptDaughter3", &fptDaughter3, "fptDaughter3/D");
  fTreeGen->Branch("fyDaughter3", &fyDaughter3, "fyDaughter3/D");
  fTreeGen->Branch("fPrimary4LHe", &fPrimary4LHe, "fPrimary4LHe/I");
  // __ generated tree __ //
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
  gTreeGen->Branch("fPDGMother", &fPDGMother, "fPDGMother/I");
  gTreeGen->Branch("fChargeMother", &fChargeMother, "fChargeMother/I");
  gTreeGen->Branch("fDecayChannel", &fDecayChannel, "fDecayChannel/I");
  gTreeGen->Branch("fmSubMother", &fmSubMother, "fmSubMother/D");
  gTreeGen->Branch("fptSubMother", &fptSubMother, "fptSubMother/D");
  gTreeGen->Branch("fySubMother", &fySubMother, "fySubMother/D");
  gTreeGen->Branch("fctSubMother", &fctSubMother, "fctSubMother/D");
  gTreeGen->Branch("fPDGMother", &fPDGMother, "fPDGMother/I");
  gTreeGen->Branch("fChargeMother", &fChargeMother, "fChargeMother/I");
  gTreeGen->Branch("fptDaughter", &fptDaughter, "fptDaughter/D");
  gTreeGen->Branch("fyDaughter", &fyDaughter, "fyDaughter/D");
  gTreeGen->Branch("fptDaughter1", &fptDaughter1, "fptDaughter1/D");
  gTreeGen->Branch("fyDaughter1", &fyDaughter1, "fyDaughter1/D");
  gTreeGen->Branch("fptDaughter2", &fptDaughter2, "fptDaughter2/D");
  gTreeGen->Branch("fyDaughter2", &fyDaughter2, "fyDaughter2/D");
  gTreeGen->Branch("fPrimary4LHe", &fPrimary4LHe, "fPrimary4LHe/I");

  PostData(1, fHistogramList);
  PostData(2, fTree);
  PostData(3, fTreeD);
  PostData(4, fTreeKF);
  PostData(5, gTree);
  PostData(6, gTreeD);
  PostData(7, gTreeKF);
  PostData(8, fTreeGen);
  PostData(9, gTreeGen);


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
void AliAnalysisTaskDoubleHypNucTree::UserExec(Option_t*) {
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

  // __ runnumber & period __ //
  fPeriod = fESDevent->GetPeriodNumber();
  frunnumber = fESDevent->GetRunNumber();

  fHistNumEvents->Fill(0);
  // __ EventCuts & TimeRangeCut __ //
  if (!fMCtrue && (frunnumber == 297219 || frunnumber == 297194 || frunnumber == 297029
		   || frunnumber == 296890 || frunnumber == 296849 || frunnumber == 296750
		   || frunnumber == 296749 || frunnumber == 297481)) fEventCuts.UseTimeRangeCut();

  // __ classify events __ //
  if (!fMCtrue && fEventCuts.AcceptEvent(fESDevent)) feventclass = 1;
  else feventclass = 0;

  if (fMCtrue) feventclass = 1;
  // __ Number of acc events __ //
  if(feventclass)fHistNumEvents->Fill(1);

  // __ centrality, 0 = V0M __ //
  Float_t centrality = -1;
  AliMultSelection* fMultSelection = (AliMultSelection*)fESDevent->FindListObject("MultSelection");
  centrality = fMultSelection->GetMultiplicityPercentile("V0M");
  fCentrality = centrality;

  // __ Trigger Selection __ //
  fEventCuts.OverrideAutomaticTriggerSelection(fTriggerMask);
  AliAnalysisTaskDoubleHypNucTree::TriggerSelection();
  // __ trigger __ //
  fTrigMB = MB;
  fTrigHMV0 = HMV0;
  fTrigHMSPD = HMSPD;
  fTrigHNU = HNU;
  fTrigHQU = HQU;
  fTrigkCentral = Central;
  fTrigkSemiCentral = SemiCentral;

  //__ SetBetheBlochParams __ //
  AliAnalysisTaskDoubleHypNucTree::SetBetheBlochParams(frunnumber);

  // __ MagneticField __ //
  kMagF = fESDevent->GetMagneticField(); //double
  fMagneticField = kMagF;               // saved as int
  KFParticle::SetField(kMagF);

  // __ MC Settings __ //
  if (fMCtrue) {
    konlyBG = kFALSE;
    konlySig = kFALSE;
    if (konlySig) konlyBG = kFALSE;
    kMCPIDCheck = kTRUE;
  }

  // __ setting up vertex reconstruction __ //
  kStandardReco = 0;
  kKFReco       = 1;
  // __ 1 = reconstruction with coord of prim vtx __ //
  // __ 2 = reconstruction with coord of tert vtx __ //
  kVariante = 1;
  // __ copy prim vtx __ //
  vertex = fESDevent->GetPrimaryVertexSPD();

  // __ initialize some utils __ //
  sublorentzsum = new TLorentzVector(0., 0., 0., 0.);
  sublorentzsum2 = new TLorentzVector(0., 0., 0., 0.);
  lorentzsum = new TLorentzVector(0., 0., 0., 0.);
  lorentzsum2 = new TLorentzVector(0., 0., 0., 0.);
  particle1 = new TLorentzVector(0., 0., 0., 0.);
  particle2 = new TLorentzVector(0., 0., 0., 0.);
  particle3 = new TLorentzVector(0., 0., 0., 0.);
  particle4 = new TLorentzVector(0., 0., 0., 0.);    
  exTrack  = new AliExternalTrackParam();
  exTrack1 = new AliExternalTrackParam();
  exTrack2 = new AliExternalTrackParam();
  exTrack3 = new AliExternalTrackParam();
  exTrack4 = new AliExternalTrackParam();
  h = new TVector3(0., 0., 0.);
  // __ only PID __ //
  if (fPIDCheckOnly) {
    AliAnalysisTaskDoubleHypNucTree::dEdxCheck();
  }// __ Analysis __ //
  else {
    if (fMCtrue) AliAnalysisTaskDoubleHypNucTree::MCGenerated();
    AliAnalysisTaskDoubleHypNucTree::FindTracks();
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
  
  PostData(1, fHistogramList);
  PostData(2, fTree);
  PostData(3, fTreeD);
  PostData(4, fTreeKF);
  PostData(5, gTree);
  PostData(6, gTreeD);
  PostData(7, gTreeKF);
  PostData(8, fTreeGen);
  PostData(9, gTreeGen);
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::dEdxCheck() {
  for (Int_t itrack = 0; itrack < fESDevent->GetNumberOfTracks(); itrack++) {
    AliESDtrack* track = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(itrack));
    Double_t momentum = track->GetInnerParam()->GetP();
    fHistdEdx->Fill(momentum * track->GetSign(), track->GetTPCsignal());
  }
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::InitArrays() {

  const Int_t nTracksEvent = fESDevent->GetNumberOfTracks();
  // __ Track Array for identified particles __ //
  He3PosCounter     = 0;
  He3NegCounter     = 0;
  PPosCounter       = 0;
  PNegCounter       = 0;
  PiPosCounter      = 0;
  PiNegCounter      = 0;
  PiPosSecCounter   = 0;
  PiNegSecCounter   = 0;

  He3PosArray.clear();
  PPosArray.clear();
  PiPosArray.clear();
  PiPosSecArray.clear();
  He3NegArray.clear();
  PNegArray.clear();
  PiNegArray.clear();
  PiNegSecArray.clear();

  // __ cuts __ //
  kDCATracksCut     = 1.2;
  kTrackPtUncertCut = 0.9;
  kPointingAngleCut = 0.95;
  kPIDSelecCut      = 4.0;
  int nsigma = 5;
  kMin4LLHMass      = 4.06;//4.106 - nsigma * 0.004;
  kMax4LLHMass      = 4.18;//4.106 + nsigma * 0.004;
  kMin4LHeMass      = 3.88;//3.924 - nsigma * 0.003;
  kMax4LHeMass      = 4.00;//3.924 + nsigma * 0.003;
  kMin3LHMass       = 2.93;//2.991 - nsigma * 0.002;
  kMax3LHMass       = 3.05;//2.991 + nsigma * 0.002; 
}
// _________________________________________________ //
Int_t AliAnalysisTaskDoubleHypNucTree::CustomTrackCut(const AliESDtrack& track, Int_t particle) {//1 = 3He, 2 = p, 3 = pi
  Float_t parxy, parz;
  track.GetImpactParameters(parxy, parz);
  Double_t mass = 0.0, charge = 1.0;
  if(particle == 1) { mass = AliPID::ParticleMass(AliPID::kHe3);    charge = 2.; }
  if(particle == 2) { mass = AliPID::ParticleMass(AliPID::kProton); charge = 1.; }
  if(particle == 3) { mass = AliPID::ParticleMass(AliPID::kPion);   charge = 1.; }
  TLorentzVector part(0.,0.,0.,0.);
  part.SetXYZM(charge*track.Px(), charge*track.Py(), charge*track.Pz(), mass);
  Int_t rval = 1;
  // __________ //
  // __ reject kinks __ //
  if (track.GetKinkIndex(0)>0)                                             rval = 0;
  // __ restrict eta __ //
  if (TMath::Abs(track.Eta()) > 0.9)                                       rval = 0;
  // __ Chi2perTPCCluster __ //
  if ((track.GetTPCchi2() / (float)track.GetTPCclusters(0)) > 5.0)         rval = 0;
  // __ require TPC refit __ //
  if (!(track.GetStatus() & AliVTrack::kTPCrefit))                         rval = 0;
  // __ pion momentum __ //
  if (particle == 3 && part.Pt() > 1.5)                                    rval = 0;
  // __ pion dca to prim vertex __ //
  if (particle == 3 && TMath::Abs(parxy) < 0.005)                          rval = 0;
  // __ ncls TPC pion __ //
  if (particle == 3 && track.GetTPCclusters(0) < 40)                       rval = 0;
  // __ ncls TPC proton __ //
  if (particle == 2 && track.GetTPCclusters(0) < 60)                       rval = 0;
  // __ ncls TPC He3 __ //
  if (particle == 1 && track.GetTPCclusters(0) < 60)                       rval = 0;
  // __ rel 1 pt uncert __ //
  if ((TMath::Sqrt(track.GetSigma1Pt2()) * part.Pt()) > kTrackPtUncertCut) rval = 0;
  //
  return rval;
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::FindTracks() {

  const Int_t nTracksEvent = fESDevent->GetNumberOfTracks();

  // __ init track arrays for identified particles __ //
  AliAnalysisTaskDoubleHypNucTree::InitArrays();

  // __ track loop for PID __ //
  for (Int_t qTracks = 0; qTracks < nTracksEvent; qTracks++) {

    AliESDtrack* trackq = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(qTracks));

    if (!trackq->GetInnerParam()) continue;

    // __ fill energy loss plot __ //   
    fHistdEdx->Fill(trackq->GetInnerParam()->GetP() * trackq->GetSign(), trackq->GetTPCsignal());

    // __ DATA __ //
    if (!fMCtrue) {
      // __ positive tracks __ //
      if (trackq->GetSign() > 0) {
	if (TMath::Abs(AliAnalysisTaskDoubleHypNucTree::Bethe(*trackq, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) <= kPIDSelecCut) {
	  if (!CustomTrackCut(*trackq, 1)) continue;
	  He3PosArray.push_back(qTracks);	  
	  He3PosCounter++;
	}
	//if (TMath::Abs(AliAnalysisTaskDoubleHypNucTree::Bethe(*trackq, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT)) <= kPIDSelecCut) {
	if (TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kProton)) <= kPIDSelecCut) {
	  if (!CustomTrackCut(*trackq, 2)) continue;
	  PPosArray.push_back(qTracks);	  
	  PPosCounter++;
	}
	if (TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kPion)) <= 3.0) {
	  if (!CustomTrackCut(*trackq, 3)) continue;
	  PiPosArray.push_back(qTracks);	  
	  PiPosCounter++;
	  PiPosSecArray.push_back(qTracks);
	  PiPosSecCounter++;
	}
      }
      // __ negative tracks __ //
      if (trackq->GetSign() < 0) {
	if (TMath::Abs(AliAnalysisTaskDoubleHypNucTree::Bethe(*trackq, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) <= kPIDSelecCut) {
	  if (!CustomTrackCut(*trackq, 1)) continue;
	  He3NegArray.push_back(qTracks);
	  He3NegCounter++;
	}
	//if (TMath::Abs(AliAnalysisTaskDoubleHypNucTree::Bethe(*trackq, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT)) <= kPIDSelecCut) {
	if (TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kProton)) <= kPIDSelecCut) {
	  if (!CustomTrackCut(*trackq, 2)) continue;
	  PNegArray.push_back(qTracks);
	  PNegCounter++;
	}
	if (TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kPion)) <= 3.0) {
	  if (!CustomTrackCut(*trackq, 3)) continue;
	  PiNegArray.push_back(qTracks);
	  PiNegCounter++;
	  PiNegSecArray.push_back(qTracks);
	  PiNegSecCounter++;
	}
      }
    }
    // ___ MC PART __ //
    if (fMCtrue) {
      label1 = trackq->GetLabel();
      AliMCParticle* tparticle = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(label1))->Particle());
      // __ positive tracks __ //
      if (trackq->GetSign() > 0) {
	if (tparticle->PdgCode() == fgkPdgCode[kPDGHelium3]) {
	  if ((kMCPIDCheck && TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kHe3)) <= kPIDSelecCut) || !kMCPIDCheck) {
	    if (!CustomTrackCut(*trackq, 1)) continue;
	    He3PosArray.push_back(qTracks);
	    He3PosCounter++;
	  }
	}
	if (tparticle->PdgCode() == fgkPdgCode[kPDGProton]) {
	  if ((kMCPIDCheck && TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kProton)) <= kPIDSelecCut) || !kMCPIDCheck) {
	    if (!CustomTrackCut(*trackq, 2)) continue;
	    PPosArray.push_back(qTracks);
	    PPosCounter++;
	  }
	}
	if (tparticle->PdgCode() == fgkPdgCode[kPDGPionPlus]) {
	  if ((kMCPIDCheck && TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kPion)) <= kPIDSelecCut) || !kMCPIDCheck) {
	    if (!CustomTrackCut(*trackq, 3)) continue;
	    PiPosArray.push_back(qTracks);
	    PiPosCounter++;
	    PiPosSecArray.push_back(qTracks);
	    PiPosSecCounter++;
	  }
	}
      }
      // __ negative tracks __ //
      if (trackq->GetSign() < 0) {
	if (tparticle->PdgCode() == fgkPdgCode[kPDGAntiHelium3]) {
	  if ((kMCPIDCheck && TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kHe3)) <= kPIDSelecCut) || !kMCPIDCheck) {
	    if (!CustomTrackCut(*trackq, 1)) continue;
	    He3NegArray.push_back(qTracks);
	    He3NegCounter++;
	  }
	}
	if (tparticle->PdgCode() == fgkPdgCode[kPDGAntiProton]) {
	  if ((kMCPIDCheck && TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kProton)) <= kPIDSelecCut) || !kMCPIDCheck) {
	    if (!CustomTrackCut(*trackq, 2)) continue;
	    PNegArray.push_back(qTracks);
	    PNegCounter++;
	  }
	}
	if (tparticle->PdgCode() == fgkPdgCode[kPDGPionMinus]) {
	  if ((kMCPIDCheck && TMath::Abs(fPID->NumberOfSigmasTPC(trackq, AliPID::kPion)) <= kPIDSelecCut) || !kMCPIDCheck) {
	    if (!CustomTrackCut(*trackq, 3)) continue;
	    PiNegArray.push_back(qTracks);
	    PiNegCounter++;
	    PiNegSecArray.push_back(qTracks);
	    PiNegSecCounter++;
	  }
	}
      }
    }
  }
  if(!He3PosCounter && !He3NegCounter) return;
  if(!PPosCounter   && !PNegCounter)   return;
  if(!PiPosCounter  && !PiNegCounter)  return;
  // _________________________________________________ // 
  cout << "...finished track finding..." << endl;
  cout << "...found "<<He3PosCounter<<" helium tracks and "    <<He3NegCounter << " anti-helium tracks..." <<endl;
  cout << "...found "<<PPosCounter <<" proton tracks and "     <<PNegCounter   <<" anti-proton tracks..."  <<endl;
  cout << "...found "<<PiNegCounter<<" pion minus tracks and " <<PiPosCounter  << " pion plus tracks..."   <<endl;
  cout << "...going to positive 4LHe..." << endl;
  //__ positive __ //
  AliAnalysisTaskDoubleHypNucTree::CalcPosDaughterHypNuc();
  AliAnalysisTaskDoubleHypNucTree::ResetVals("");
  //__ negative __ //
  cout << "...going to negative 4LHe..." << endl;
  AliAnalysisTaskDoubleHypNucTree::CalcNegDaughterHypNuc();
  AliAnalysisTaskDoubleHypNucTree::ResetVals("");
  // __ reset __ //
  ResetVals("Event");  
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::CalcPosDaughterHypNuc() {
  particle1->SetXYZM(0.,0.,0.,0.);
  particle2->SetXYZM(0.,0.,0.,0.);
  particle3->SetXYZM(0.,0.,0.,0.);
  particle4->SetXYZM(0.,0.,0.,0.);
  lorentzsum->SetXYZM(0.,0.,0.,0.);
  sublorentzsum->SetXYZM(0.,0.,0.,0.);
  cout<<"...Analyzing decay channel 1..."<<endl;
  // ____________________________ first track _______________________________ //
  for (const Int_t& He3PosTracks : He3PosArray) {

    track1 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(He3PosTracks));
    // __ MC part __ //
    if (fMCtrue && konlySig) {

      label1 = track1->GetLabel();
      labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));

      ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());
      if (ParticleMother1->PdgCode() != fgkPdgCode[kPDGHyperHelium4]) continue;
      if (ParticleMother1) delete ParticleMother1;
    }
    // __ set pxpypz __ //
    particle1->SetXYZM(0.,0.,0.,0.);
    particle1->SetXYZM(2.*track1->Px(), 2.*track1->Py(), 2.*track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
    // ____________________________ second track _______________________________ //
    for (const Int_t& PPosTracks : PPosArray) {

      // __ reject using same track twice __ //
      if (He3PosTracks == PPosTracks) continue;

      track2 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(PPosTracks));
      // __ MC part __ //
      if (fMCtrue && konlySig) {

	label2 = track2->GetLabel();
	labelMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label2));

	if (labelMother1 != labelMother2) continue;

	ParticleMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother2))->Particle());
	if (ParticleMother2->PdgCode() != fgkPdgCode[kPDGHyperHelium4]) continue;
	if (ParticleMother2) delete ParticleMother2;
      }     
      // __ dca rejection __ //     
      if (track1->GetDCA(track2, kMagF, xthiss, xpp) > kDCATracksCut) continue;
      // __ set pxpypz __ //
      particle2->SetXYZM(0.,0.,0.,0.);
      particle2->SetXYZM(track2->Px(), track2->Py(), track2->Pz(), AliPID::ParticleMass(AliPID::kProton));
      // ____________________________ third track _______________________________ //
      for (const Int_t& PiNegTracks : PiNegArray) {
	  
	// __ reject using same track twice __ //
	if (PPosTracks == PiNegTracks || He3PosTracks == PiNegTracks) continue;

	track3 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(PiNegTracks));
	// __ MC part __ //
	if (fMCtrue && konlySig) {

	  label3 = track3->GetLabel();
	  labelMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label3));

	  if (labelMother3 != labelMother2) continue;

	  ParticleMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	  if (ParticleMother3->PdgCode() != fgkPdgCode[kPDGHyperHelium4]) continue;
	  if (ParticleMother3) delete ParticleMother3;
	}	
	// __ dca rejection  __ //
	if (track1->GetDCA(track2, kMagF, xthiss, xpp) > kDCATracksCut) continue;
	if (track1->GetDCA(track3, kMagF, xthiss, xpp) > kDCATracksCut) continue;
	if (track2->GetDCA(track3, kMagF, xthiss, xpp) > kDCATracksCut) continue;
	// __ set pxpypz __ //
	particle3->SetXYZM(0.,0.,0.,0.);
	particle3->SetXYZM(track3->Px(), track3->Py(), track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
	//______ daughter hypernucleus mass _______ //				
	sublorentzsum->SetXYZM(0., 0., 0., 0.);
	*sublorentzsum = *particle1 + *particle2 + *particle3;
	if (sublorentzsum->M() > kMax4LHeMass || sublorentzsum->M() < kMin4LHeMass) {	
	  AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	  continue;
	}
	// _________________________________________________ //
	// __ standard vertex reconstruction __ //
	int status_standardreco = 0;
	if(kStandardReco) status_standardreco = AliAnalysisTaskDoubleHypNucTree::StandardReconstruction("4LHe", 1, 1);
	// _________________________________________________ //
	// __ vertex reconstruction with KF method __ //
	int status_KFreco = 0;
	if(kKFReco) status_KFreco = AliAnalysisTaskDoubleHypNucTree::KFReconstruction("4LHe", 1, 1);
	// _________________________________________________ //
	if(!status_standardreco && !status_KFreco) {	     
	  AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	  continue;
	}		
	// _________________________________________________ //
	// __ MC part __ //
	if (fMCtrue) {
	  ParticleMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	  ParticleMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother2))->Particle());
	  ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());
	  
	  if (   ParticleMother1->PdgCode() == fgkPdgCode[kPDGHyperHelium4]
	      && ParticleMother2->PdgCode() == fgkPdgCode[kPDGHyperHelium4]
	      && ParticleMother3->PdgCode() == fgkPdgCode[kPDGHyperHelium4]
	      && labelMother1 == labelMother2
	      && labelMother2 == labelMother3) {
	    if (   TMath::Abs(label1) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGHelium3]))
		&& TMath::Abs(label2) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGProton]))
		&& TMath::Abs(label3) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGPionMinus]))) {
	      fmctruth = 1;
	      fPrimary4LHe = mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(labelMother1)) == 0;
	    }
	  }
	  if (ParticleMother1) delete ParticleMother1;
	  if (ParticleMother2) delete ParticleMother2;
	  if (ParticleMother3) delete ParticleMother3;
	}
	// _________________________________________________ //
	// __ fill tree __ //
	if(kStandardReco) gTree->Fill();
	gTreeD->Fill();
	if(kKFReco) gTreeKF->Fill();
	// __ reset __ //
	AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	// _________________________________________________ //
	// __ go to 4LLH __ //
	cout<<"...positive 4LHe..."<<endl;
	cout<<"...going to positive 4LLH..."<<endl;
	AliAnalysisTaskDoubleHypNucTree::CalcPosMotherHypNuc(He3PosTracks, PPosTracks, PiNegTracks, 1);
      } //track3
    } //track2
  } //track1
  // _____________________________________________________________________________________________________________________ //
  cout<<"...Analyzing decay channel 2..."<<endl;
  particle1->SetXYZM(0.,0.,0.,0.);
  particle2->SetXYZM(0.,0.,0.,0.);
  particle3->SetXYZM(0.,0.,0.,0.);
  particle4->SetXYZM(0.,0.,0.,0.);
  lorentzsum->SetXYZM(0.,0.,0.,0.);
  sublorentzsum->SetXYZM(0.,0.,0.,0.);
  // ____________________________ first track _______________________________ //
  for (const Int_t& He3PosTracks : He3PosArray) {

    track1 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(He3PosTracks));
    // __ MC part __ //
    if (fMCtrue && konlySig) {

      label1 = track1->GetLabel();
      labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));

      ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());
      if (ParticleMother1->PdgCode() != fgkPdgCode[kPDGHyperHydrogen3]) continue;
      if (ParticleMother1) delete ParticleMother1;
    }
    // __ set pxpypz __ //
    particle1->SetXYZM(0.,0.,0.,0.);
    particle1->SetXYZM(2.*track1->Px(), 2.*track1->Py(), 2.*track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));    
    // ____________________________ third track _______________________________ //
    for (const Int_t& PiNegTracks : PiNegArray) {
	
      // __ reject using same track twice __ //
      if (He3PosTracks == PiNegTracks) continue;

      track3 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(PiNegTracks));
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
      if (track1->GetDCA(track3, kMagF, xthiss, xpp) > kDCATracksCut) continue;
      // __ set pxpypz __ //
      particle3->SetXYZM(0.,0.,0.,0.);
      particle3->SetXYZM(track3->Px(), track3->Py(), track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
      //______ daughter hypernucleus mass _______ //			
      sublorentzsum->SetXYZM(0., 0., 0., 0.);
      *sublorentzsum = *particle1 + *particle3;
      if (sublorentzsum->M() > kMax3LHMass || sublorentzsum->M() < kMin3LHMass) {	
	AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	continue;
      }
      // _________________________________________________ //
      // __ standard vertex reconstruction __ //
      int status_standardreco = 0;
      if(kStandardReco) status_standardreco = AliAnalysisTaskDoubleHypNucTree::StandardReconstruction("3LH", 2, 1);
      // _________________________________________________ //
      // __ vertex reconstruction with KF method __ //
      int status_KFreco = 0;
      if(kKFReco) status_KFreco = AliAnalysisTaskDoubleHypNucTree::KFReconstruction("3LH", 2, 1);
      // _________________________________________________ //
      if(!status_standardreco && !status_KFreco) {	     
	AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	continue;
      }           
      // _________________________________________________ //
      // __ MC part __ //
      if (fMCtrue) {
	ParticleMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());
	if (   ParticleMother1->PdgCode() == fgkPdgCode[kPDGHyperHydrogen3]
	    && ParticleMother3->PdgCode() == fgkPdgCode[kPDGHyperHydrogen3]
	    && labelMother1 == labelMother3) {
	  if (   TMath::Abs(label1) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGHelium3]))
	      && TMath::Abs(label3) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGPionMinus]))) {
	    fmctruth = 1;
	    fPrimary4LHe = mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(labelMother1)) == 0;
	  }
	}
	if (ParticleMother1) delete ParticleMother1;
	if (ParticleMother3) delete ParticleMother3;
      }
      // _________________________________________________ //
      // __ fill tree __ //
      if(kStandardReco) gTree->Fill();
      gTreeD->Fill();
      if(kKFReco) gTreeKF->Fill();
      // __ reset __ //
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      // _________________________________________________ //
      // __ go to 4LLH __ //
      cout<<"...positive 3LH..."<<endl;
      cout<<"...going to positive 4LLH..."<<endl;
      AliAnalysisTaskDoubleHypNucTree::CalcPosMotherHypNuc(He3PosTracks, PiNegTracks, 0, 2);
    }//track3
  }//track1
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::CalcNegDaughterHypNuc() {
  particle1->SetXYZM(0.,0.,0.,0.);
  particle2->SetXYZM(0.,0.,0.,0.);
  particle3->SetXYZM(0.,0.,0.,0.);
  particle4->SetXYZM(0.,0.,0.,0.);
  lorentzsum->SetXYZM(0.,0.,0.,0.);
  sublorentzsum->SetXYZM(0.,0.,0.,0.);
  cout<<"...Analyzing decay channel 1..."<<endl;
  // ____________________________ first track _______________________________ //
  for (const Int_t& He3NegTracks : He3NegArray) {

    track1 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(He3NegTracks));
    // __ MC part __ //
    if (fMCtrue && konlySig) {

      label1 = track1->GetLabel();
      labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));

      ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());
      if (ParticleMother1->PdgCode() != fgkPdgCode[kPDGAntiHyperHelium4]) continue;
      if (ParticleMother1) delete ParticleMother1;
    }    
    // __ set pxpypz __ //
    particle1->SetXYZM(0.,0.,0.,0.);
    particle1->SetXYZM(2.*track1->Px(), 2.*track1->Py(), 2.*track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
    // ____________________________ second track _______________________________ //
    for (const Int_t& PNegTracks : PNegArray) {

      // __ reject using same track twice __ //
      if (He3NegTracks == PNegTracks) continue;

      track2 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(PNegTracks));
      // __ MC part __ //
      if (fMCtrue && konlySig) {

	label2 = track2->GetLabel();
	labelMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label2));

	if (labelMother1 != labelMother2) continue;

	ParticleMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother2))->Particle());
	if (ParticleMother2->PdgCode() != fgkPdgCode[kPDGAntiHyperHelium4]) continue;
	if (ParticleMother2) delete ParticleMother2;
      }    
      // __ dca rejection __ //      
      if (track1->GetDCA(track2, kMagF, xthiss, xpp) > kDCATracksCut) continue;
      // __ set pxpypz __ //
      particle2->SetXYZM(0.,0.,0.,0.);
      particle2->SetXYZM(track2->Px(), track2->Py(), track2->Pz(), AliPID::ParticleMass(AliPID::kProton));
      // ____________________________ third track _______________________________ //
      for (const Int_t& PiPosTracks : PiPosArray) {
	
	// __ reject using same track twice __ //
	if (PNegTracks == PiPosTracks || He3NegTracks == PiPosTracks) continue;

	track3 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(PiPosTracks));
	// __ MC part __ //
	if (fMCtrue && konlySig) {

	  label3 = track3->GetLabel();
	  labelMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label3));
	  labelGrandMother3 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother3));

	  if (labelMother3 != labelMother2) continue;

	  ParticleGrandMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	  if (ParticleGrandMother3->PdgCode() != fgkPdgCode[kPDGAntiHyperHelium4]) continue;
	  if (ParticleGrandMother3) delete ParticleGrandMother3;
	}
	// __ dca rejection __ //
	if (track1->GetDCA(track2, kMagF, xthiss, xpp) > kDCATracksCut) continue;
	if (track1->GetDCA(track3, kMagF, xthiss, xpp) > kDCATracksCut) continue;
	if (track2->GetDCA(track3, kMagF, xthiss, xpp) > kDCATracksCut) continue;
	// __ set pxpypz __ //
	particle3->SetXYZM(0.,0.,0.,0.);
	particle3->SetXYZM(track3->Px(), track3->Py(), track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
	//______ daughter hypernucleus mass _______ //				
	sublorentzsum->SetXYZM(0., 0., 0., 0.);
	*sublorentzsum = *particle1 + *particle2 + *particle3;
	if (sublorentzsum->M() > kMax4LHeMass || sublorentzsum->M() < kMin4LHeMass) {	 
	  AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	  continue;
	}
	// _________________________________________________ //
	// __ standard vertex reconstruction __ //
	int status_standardreco = 0;
	if(kStandardReco) status_standardreco = AliAnalysisTaskDoubleHypNucTree::StandardReconstruction("4LHe", 1, -1);
	// _________________________________________________ //
	// __ vertex reconstruction with KF method __ //
	int status_KFreco = 0;
	if(kKFReco) status_KFreco = AliAnalysisTaskDoubleHypNucTree::KFReconstruction("4LHe", 1, -1);
	// _________________________________________________ //
	if(!status_standardreco && !status_KFreco) {	     
	  AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	  continue;
	}
	// _________________________________________________ //	
	// __ MC part __ //
	if (fMCtrue) {
	  ParticleMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	  ParticleMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother2))->Particle());
	  ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());

	  if (   ParticleMother1->PdgCode() == fgkPdgCode[kPDGAntiHyperHelium4]
	      && ParticleMother2->PdgCode() == fgkPdgCode[kPDGAntiHyperHelium4]
	      && ParticleMother3->PdgCode() == fgkPdgCode[kPDGAntiHyperHelium4]
	      && labelMother1 == labelMother2
	      && labelMother2 == labelMother3) {
	    if (   TMath::Abs(label1) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGAntiHelium3]))
		&& TMath::Abs(label2) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGAntiProton]))
		&& TMath::Abs(label3) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGPionPlus]))) {
	      fmctruth = 1;
	      fPrimary4LHe = mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(labelMother1)) == 0;
	    }
	  }
	  if (ParticleMother1) delete ParticleMother1;
	  if (ParticleMother2) delete ParticleMother2;
	  if (ParticleMother3) delete ParticleMother3;
	}
	// _________________________________________________ //
	// __ fill tree __ //
	if(kStandardReco) gTree->Fill();
	gTreeD->Fill();
	if(kKFReco) gTreeKF->Fill();
	// __ reset __ //
	AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	// _________________________________________________ //
	// __ go to 4LLH __ //
	cout<<"...negative 4LHe..."<<endl;
	cout<<"...going to negative 4LLH..."<<endl;
	AliAnalysisTaskDoubleHypNucTree::CalcNegMotherHypNuc(He3NegTracks, PNegTracks, PiPosTracks, 1);
      } //track3
    } //track2
  } //track1
  // _________________________________________________________________________ //  
  cout<<"...Analyzing decay channel 2..."<<endl;
  particle1->SetXYZM(0.,0.,0.,0.);
  particle2->SetXYZM(0.,0.,0.,0.);
  particle3->SetXYZM(0.,0.,0.,0.);
  particle4->SetXYZM(0.,0.,0.,0.);
  lorentzsum->SetXYZM(0.,0.,0.,0.);
  sublorentzsum->SetXYZM(0.,0.,0.,0.);
  // ____________________________ first track _______________________________ //
  for (const Int_t& He3NegTracks : He3NegArray) {

    track1 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(He3NegTracks));
    // __ MC part __ //
    if (fMCtrue && konlySig) {

      label1 = track1->GetLabel();
      labelMother1 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label1));

      ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());
      if (ParticleMother1->PdgCode() != fgkPdgCode[kPDGAntiHyperHydrogen3]) continue;
      if (ParticleMother1) delete ParticleMother1;
    }
    // __ set pxpypz __ //
    particle1->SetXYZM(0.,0.,0.,0.);
    particle1->SetXYZM(2.*track1->Px(), 2.*track1->Py(), 2.*track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
    // ____________________________ third track _______________________________ //
    for (const Int_t& PiPosTracks : PiPosArray) {
	
      // __ reject using same track twice __ //
      if (He3NegTracks == PiPosTracks) continue;

      track3 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(PiPosTracks));
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
      if (track1->GetDCA(track3, kMagF, xthiss, xpp) > kDCATracksCut) continue;
      // __ set pxpypz __ //
      particle3->SetXYZM(0.,0.,0.,0.);
      particle3->SetXYZM(track3->Px(), track3->Py(), track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
      //______ daughter hypernucleus mass _______ //			
      sublorentzsum->SetXYZM(0., 0., 0., 0.);
      *sublorentzsum = *particle1 + *particle3;
      if ((sublorentzsum->M() > kMax3LHMass || sublorentzsum->M() < kMin3LHMass)) {	 
	AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	continue;
      }	
      // _________________________________________________ //
      // __ standard vertex reconstruction __ //
      int status_standardreco = 0;
      if(kStandardReco) status_standardreco = AliAnalysisTaskDoubleHypNucTree::StandardReconstruction("3LH", 2, -1);
      // _________________________________________________ //
      // __ vertex reconstruction with KF method __ //
      int status_KFreco = 0;
      if(kKFReco) status_KFreco = AliAnalysisTaskDoubleHypNucTree::KFReconstruction("3LH", 2, -1);
      // _________________________________________________ //
      if(!status_standardreco && !status_KFreco) {	     
	AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	continue;
      }
      // _________________________________________________ //     
      // __ MC part __ //
      if (fMCtrue) {
	ParticleMother3 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother3))->Particle());
	ParticleMother1 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother1))->Particle());
	if (   ParticleMother1->PdgCode() == fgkPdgCode[kPDGAntiHyperHydrogen3]
	    && ParticleMother3->PdgCode() == fgkPdgCode[kPDGAntiHyperHydrogen3]
	    && labelMother1 == labelMother3) {
	  if (   TMath::Abs(label1) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGAntiHelium3]))
	      && TMath::Abs(label3) == TMath::Abs(GetLabel(labelMother1, fgkPdgCode[kPDGPionPlus]))) {
	    fmctruth = 1;
	    fPrimary4LHe = mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(labelMother1)) == 0;
	  }
	}
	if (ParticleMother1) delete ParticleMother1;
	if (ParticleMother3) delete ParticleMother3;
      }
      // _________________________________________________ //
      // __ fill tree __ //
      if(kStandardReco)gTree->Fill();
      gTreeD->Fill();
      if(kKFReco)gTreeKF->Fill();
      // __ reset __ //
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      // _________________________________________________ //
      // __ go to 4LLH __ //
      cout<<"...negative 3LH..."<<endl;
      cout<<"...going to negative 4LLH..."<<endl;
      AliAnalysisTaskDoubleHypNucTree::CalcNegMotherHypNuc(He3NegTracks, PiPosTracks, 0, 2);
    }//track3
  }//track1
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::CalcPosMotherHypNuc(Int_t Track1Entry, Int_t Track2Entry, Int_t Track3Entry, Int_t kDecayChannel) {
  Int_t kCharge = 1;
  particle1->SetXYZM(0.,0.,0.,0.);
  particle2->SetXYZM(0.,0.,0.,0.);
  particle3->SetXYZM(0.,0.,0.,0.);
  particle4->SetXYZM(0.,0.,0.,0.);
  lorentzsum->SetXYZM(0.,0.,0.,0.);
  sublorentzsum->SetXYZM(0.,0.,0.,0.);
  // _________________________________________________ //     
  if(kDecayChannel == 1){
    if(!Track1Entry) return; if(!Track2Entry) return; if(!Track3Entry) return;
    track1 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(Track1Entry));
    track2 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(Track2Entry));
    track3 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(Track3Entry));
    // __ set pxpypz __ //
    particle1->SetXYZM(2.*track1->Px(), 2.*track1->Py(), 2.*track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
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

      track4 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(PiNegSecTracks));
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
      particle4->SetXYZM(0.,0.,0.,0.);
      particle4->SetXYZM(track4->Px(), track4->Py(), track4->Pz(), AliPID::ParticleMass(AliPID::kPion));
      // _________________________________________________ //      
      lorentzsum->SetXYZM(0., 0., 0., 0.);
      *lorentzsum = *particle1 + *particle2 + *particle3 + *particle4;	 
      // _________________________________________________ //
      // __ select mother hypernucleus by inv mass __ //
      if (lorentzsum->M() > kMax4LLHMass || lorentzsum->M() < kMin4LLHMass) {	   
	AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	continue;
      }	  	      
      // _________________________________________________ //
      // __ standard vertex reconstruction __ //
      int status_standardreco = 0;
      if(kStandardReco) status_standardreco = AliAnalysisTaskDoubleHypNucTree::StandardReconstruction("4LLH", kDecayChannel, kCharge);
      // _________________________________________________ //
      // __ vertex reconstruction with KF method __ //
      int status_KFreco = 0;
      if(kKFReco) status_KFreco = AliAnalysisTaskDoubleHypNucTree::KFReconstruction("4LLH", kDecayChannel, kCharge);
      // _________________________________________________ //
      if(!status_standardreco && !status_KFreco) {	     
	AliAnalysisTaskDoubleHypNucTree::ResetVals("");
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
	
	if (   ParticleGrandMother1->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	    && ParticleGrandMother2->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	    && ParticleGrandMother3->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	    && ParticleMother4->PdgCode()      == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]) {

	  if (labelGrandMother1 == labelGrandMother2 && labelGrandMother2 == labelGrandMother3 && labelGrandMother3 == labelMother4) {

	    if (   TMath::Abs(label1) == TMath::Abs(GetLabel(labelGrandMother1, kCharge * fgkPdgCode[kPDGHelium3],   kCharge * fgkPdgCode[kPDGHyperHelium4]))
		&& TMath::Abs(label2) == TMath::Abs(GetLabel(labelGrandMother1, kCharge * fgkPdgCode[kPDGProton],    kCharge * fgkPdgCode[kPDGHyperHelium4]))
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
      // __ fill tree __ //
      if(kStandardReco)fTree->Fill();
      fTreeD->Fill();
      if(kKFReco)fTreeKF->Fill();		      
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      // _________________________________________________ //
    }//track4
  }//kDecayChannel
  // ____________________________________________________________________________________________________________________________________ //
  if(kDecayChannel == 2){
    if(!Track1Entry) return; if(!Track2Entry) return;
    track1 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(Track1Entry));
    track3 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(Track2Entry));
    // __ set pxpypz __ //
    particle1->SetXYZM(2.*track1->Px(), 2.*track1->Py(), 2.*track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
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

      track2 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(PPosTracks));
      // __ MC part __ //
      if(fMCtrue && konlySig){
	
	label2 = track2->GetLabel();
	labelMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label2));
	labelGrandMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(labelMother2));

	if (labelGrandMother3 != labelMother2) continue;

	ParticleMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother2))->Particle());
	if (ParticleMother2->PdgCode() != fgkPdgCode[kPDGDoubleHyperHydrogen4]) continue;
	if (ParticleMother2) delete ParticleMother2;
      }
      // __ set pxpypz __ //
      particle2->SetXYZM(0.,0.,0.,0.);
      particle2->SetXYZM(track2->Px(), track2->Py(), track2->Pz(), AliPID::ParticleMass(AliPID::kProton));
      // ____________________________ fourth track _______________________________ //
      for (const Int_t& PiNegSecTracks : PiNegSecArray) {

	// __ reject using same track twice __ //
	if (PPosTracks == PiNegSecTracks || Track1Entry == PiNegSecTracks || Track3Entry == PiNegSecTracks) continue;

	track4 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(PiNegSecTracks));
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
	if (track2->GetDCA(track4, kMagF, xthiss, xpp) > kDCATracksCut) continue;
	// __ set pxpypz __ //
	particle4->SetXYZM(0.,0.,0.,0.);
	particle4->SetXYZM(track4->Px(), track4->Py(), track4->Pz(), AliPID::ParticleMass(AliPID::kPion));
	// _________________________________________________ //	  	 	
	lorentzsum->SetXYZM(0., 0., 0., 0.);
	*lorentzsum = *particle1 + *particle2 + *particle3 + *particle4;	 
	// _________________________________________________ //
	// __ select mother hypernucleus by inv mass __ //
	if (lorentzsum->M() > kMax4LLHMass || lorentzsum->M() < kMin4LLHMass) {	   
	  AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	  continue;
	}	  	 	
	// _________________________________________________ //
	// __ standard vertex reconstruction __ //
	int status_standardreco = 0;
	if(kStandardReco) status_standardreco = AliAnalysisTaskDoubleHypNucTree::StandardReconstruction("4LLH", kDecayChannel, kCharge);
	// _________________________________________________ //
	// __ vertex reconstruction with KF method __ //
	int status_KFreco = 0;
	if(kKFReco) status_KFreco = AliAnalysisTaskDoubleHypNucTree::KFReconstruction("4LLH", kDecayChannel, kCharge);
	// _________________________________________________ //
	if(!status_standardreco && !status_KFreco) {	     
	  AliAnalysisTaskDoubleHypNucTree::ResetVals("");
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
	  
	  if (   ParticleGrandMother1->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	      && ParticleMother2->PdgCode()      == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	      && ParticleGrandMother3->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	      && ParticleMother4->PdgCode()      == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]) {

	    if (labelGrandMother1 == labelMother2 && labelMother2 == labelGrandMother3 && labelGrandMother3 == labelMother4) {

	      if (   TMath::Abs(label1) == TMath::Abs(GetLabel(labelGrandMother1, kCharge * fgkPdgCode[kPDGHelium3],   kCharge * fgkPdgCode[kPDGHyperHydrogen3]))
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
	// __ fill tree __ //
	if(kStandardReco)fTree->Fill();
	fTreeD->Fill();
	if(kKFReco)fTreeKF->Fill();		      
	AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	// _________________________________________________ //
      }//track4
    }//track2
  }//kDecayChannel
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::CalcNegMotherHypNuc(Int_t Track1Entry, Int_t Track2Entry, Int_t Track3Entry, Int_t kDecayChannel) {
  Int_t kCharge = -1;
  particle1->SetXYZM(0.,0.,0.,0.);
  particle2->SetXYZM(0.,0.,0.,0.);
  particle3->SetXYZM(0.,0.,0.,0.);
  particle4->SetXYZM(0.,0.,0.,0.);
  lorentzsum->SetXYZM(0.,0.,0.,0.);
  sublorentzsum->SetXYZM(0.,0.,0.,0.);
  // _________________________________________________ //
  if(kDecayChannel == 1){
    if(!Track1Entry) return; if(!Track2Entry) return; if(!Track3Entry) return;
    track1 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(Track1Entry));
    track2 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(Track2Entry));
    track3 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(Track3Entry));
    // __ set pxpypz __ //
    particle1->SetXYZM(2.*track1->Px(), 2.*track1->Py(), 2.*track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
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

      track4 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(PiPosSecTracks));
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
      particle4->SetXYZM(0.,0.,0.,0.);
      particle4->SetXYZM(track4->Px(), track4->Py(), track4->Pz(), AliPID::ParticleMass(AliPID::kPion));  
      // _________________________________________________ //	  	       
      lorentzsum->SetXYZM(0., 0., 0., 0.);
      *lorentzsum = *particle1 + *particle2 + *particle3 + *particle4;	 
      // _________________________________________________ //
      // __ select mother hypernucleus by inv mass __ //
      if (lorentzsum->M() > kMax4LLHMass || lorentzsum->M() < kMin4LLHMass) {	   
	AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	continue;
      }	  	            
      // _________________________________________________ //
      // __ standard vertex reconstruction __ //
      int status_standardreco = 0;
      if(kStandardReco) status_standardreco = AliAnalysisTaskDoubleHypNucTree::StandardReconstruction("4LLH", kDecayChannel, kCharge);
      // _________________________________________________ //
      // __ vertex reconstruction with KF method __ //
      int status_KFreco = 0;
      if(kKFReco) status_KFreco = AliAnalysisTaskDoubleHypNucTree::KFReconstruction("4LLH", kDecayChannel, kCharge);
      // _________________________________________________ //
      if(!status_standardreco && !status_KFreco) {	     
	AliAnalysisTaskDoubleHypNucTree::ResetVals("");
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
	
	if (   ParticleGrandMother1->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	    && ParticleGrandMother2->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	    && ParticleGrandMother3->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	    && ParticleMother4->PdgCode()      == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]) {

	  if (labelGrandMother1 == labelGrandMother2 && labelGrandMother2 == labelGrandMother3 && labelGrandMother3 == labelMother4) {

	    if (   TMath::Abs(label1) == TMath::Abs(GetLabel(labelGrandMother1, kCharge * fgkPdgCode[kPDGHelium3],   kCharge * fgkPdgCode[kPDGHyperHelium4]))
		&& TMath::Abs(label2) == TMath::Abs(GetLabel(labelGrandMother1, kCharge * fgkPdgCode[kPDGProton],    kCharge * fgkPdgCode[kPDGHyperHelium4]))
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
      // __ fill tree __ //
      if(kStandardReco)fTree->Fill();
      fTreeD->Fill();
      if(kKFReco)fTreeKF->Fill();		      
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      // _________________________________________________ //
    }//track4
  }//kDecayChannel
  // _____________________________________________________________________ //
  if(kDecayChannel == 2){
    if(!Track1Entry) return; if(!Track2Entry) return;
    track1 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(Track1Entry));
    track3 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(Track2Entry));
    // __ set pxpypz __ //
    particle1->SetXYZM(2.*track1->Px(), 2.*track1->Py(), 2.*track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
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

      track2 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(PNegTracks));
      // __ MC part __ //
      if(fMCtrue && konlySig){
	
	label2 = track2->GetLabel();
	labelMother2 = mcEvent->GetLabelOfParticleMother(TMath::Abs(label2));

	if (labelGrandMother3 != labelMother2) continue;

	ParticleMother2 = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother2))->Particle());
	if (ParticleMother2->PdgCode() != kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]) continue;
	if (ParticleMother2) delete ParticleMother2;
      }
      // __ set pxpypz __ //
      particle2->SetXYZM(0.,0.,0.,0.);
      particle2->SetXYZM(track2->Px(), track2->Py(), track2->Pz(), AliPID::ParticleMass(AliPID::kProton));
      // ____________________________ fourth track _______________________________ //
      for (const Int_t& PiPosSecTracks : PiPosSecArray) {

	// __ reject using same track twice __ //
	if (PNegTracks == PiPosSecTracks || Track1Entry == PiPosSecTracks || Track2Entry == PiPosSecTracks) continue;

	track4 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(PiPosSecTracks));
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
	if (track2->GetDCA(track4, kMagF, xthiss, xpp) > kDCATracksCut) continue;
	// __ set pxpypz __ //
	particle4->SetXYZM(0.,0.,0.,0.);
	particle4->SetXYZM(track4->Px(), track4->Py(), track4->Pz(), AliPID::ParticleMass(AliPID::kPion));
	// _________________________________________________ //	  	 	
	lorentzsum->SetXYZM(0., 0., 0., 0.);
	*lorentzsum = *particle1 + *particle2 + *particle3 + *particle4;	 
	// _________________________________________________ //
	// __ select mother hypernucleus by inv mass __ //
	if (lorentzsum->M() > kMax4LLHMass || lorentzsum->M() < kMin4LLHMass) {	   
	  AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	  continue;
	}	  	 		
	// _________________________________________________ //
	// __ standard vertex reconstruction __ //
	int status_standardreco = 0;
	if(kStandardReco) status_standardreco = AliAnalysisTaskDoubleHypNucTree::StandardReconstruction("4LLH", kDecayChannel, kCharge);
	// _________________________________________________ //
	// __ vertex reconstruction with KF method __ //
	int status_KFreco = 0;
	if(kKFReco) status_KFreco = AliAnalysisTaskDoubleHypNucTree::KFReconstruction("4LLH", kDecayChannel, kCharge);
	// _________________________________________________ //
	if(!status_standardreco && !status_KFreco) {	     
	  AliAnalysisTaskDoubleHypNucTree::ResetVals("");
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
	  
	  if (   ParticleGrandMother1->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	      && ParticleMother2->PdgCode()      == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	      && ParticleGrandMother3->PdgCode() == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]
	      && ParticleMother4->PdgCode()      == kCharge * fgkPdgCode[kPDGDoubleHyperHydrogen4]) {

	    if (labelGrandMother1 == labelMother2 && labelMother2 == labelGrandMother3 && labelGrandMother3 == labelMother4) {

	      if (   TMath::Abs(label1) == TMath::Abs(GetLabel(labelGrandMother1, kCharge * fgkPdgCode[kPDGHelium3],   kCharge * fgkPdgCode[kPDGHyperHydrogen3]))
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
	// __ fill tree __ //
	if(kStandardReco)fTree->Fill();
	fTreeD->Fill();
	if(kKFReco)fTreeKF->Fill();		      
	AliAnalysisTaskDoubleHypNucTree::ResetVals("");
	// _________________________________________________ //
      }//track4
    }//track2
  }//kDecayChannel
}
// _________________________________________________ //
Int_t AliAnalysisTaskDoubleHypNucTree::StandardReconstruction(TString Mother, Int_t kDecayChannel, Int_t ksign) {
  Double_t impar[2];
  UShort_t idThree[3] = { 0, 1, 2 };
  UShort_t idTwo[2] = { 0, 1 };
  Int_t    kFinderAlgo = 1; //1 StrLinVertexFinderMinDist(UseWeights=1) (default), 2: StrLinVertexFinderMinDist(UseWeights=0), 3: HelixVertexFinder(), 4: VertexFinder(1), 5: VertexFinder(0)
  // _________________________________ //      
  if (Mother == "3LH") {
    exTrack->Reset(); exTrack1->Reset(); exTrack2->Reset(); exTrack3->Reset(); exTrack4->Reset();
    exTrack1->CopyFromVTrack(track1);
    exTrack3->CopyFromVTrack(track3);
    // _________________________________ //
    // __ get prim vtx __ //
    primVertex = new AliESDVertex(*vertex);
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
    // __ init AliVertexerTracks __ //
    tertvertexer = new AliVertexerTracks(fESDevent->GetMagneticField());
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
    fTertVertChi2 = tertVertex->GetChi2();
    fTertVertNDF = tertVertex->GetNDF();
    // _________________________________________________ //
    // __ coord of tert vtx __ //
    TertVertex[0] = tertVertex->GetX();
    TertVertex[1] = tertVertex->GetY();
    TertVertex[2] = tertVertex->GetZ();
    fTertVertexX = TertVertex[0];
    fTertVertexY = TertVertex[1];
    fTertVertexZ = TertVertex[2];
    // _________________________________ //
    // __ propagate tracks and get impact parameters __ //
    if (exTrack1->PropagateToDCA(tertVertex, kMagF, 10, impar)) fPropDCADaughter = 1;
    fImParDaughter = impar[0];
    fImParzDaughter = impar[1];
    if (exTrack3->PropagateToDCA(tertVertex, kMagF, 10, impar)) fPropDCADaughter2 = 1;
    fImParDaughter2 = impar[0];
    fImParzDaughter2 = impar[1];
    // _________________________________ //
    if (tertVertex) delete tertVertex;
    // _________________________________________________ //
    // __ dca between tracks __ //
    fDCA3B1 = TMath::Abs(exTrack3->GetDCA(exTrack1, kMagF, xthiss, xpp));
    // _________________________________________________ //
    particle1->SetXYZM(2. * exTrack1->Px(), 2. * exTrack1->Py(), 2. * exTrack1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
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
    if (fSubPA < kPointingAngleCut) {      
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      return 0;
    }   
    // _________________________________________________ // 
    // __ cuts from 04/2020 __ //
    fDcaDaughtero  = exTrack1->GetD(fPrimVertexX, fPrimVertexY, kMagF);
    fDcaDaughter1o = exTrack2->GetD(fPrimVertexX, fPrimVertexY, kMagF);
    fDcaDaughter2o = exTrack3->GetD(fPrimVertexX, fPrimVertexY, kMagF);
    // _________________________________________________ //
    // __ track information __ //
    AliAnalysisTaskDoubleHypNucTree::SetDaughterInformation("3LH", 2);
    // _________________________________________________ //	  	
    // __ daughter hypernucleus information __ //
    fPDGMother    = ksign*fgkPdgCode[kPDGHyperHydrogen3];
    fChargeMother = ksign;
    fRecoMethod   = 1;
    fmSubMother   = sublorentzsum->M();
    fESubMother   = sublorentzsum->E();
    fpxSubMother  = sublorentzsum->Px();
    fpySubMother  = sublorentzsum->Py();
    fpzSubMother  = sublorentzsum->Pz();
    fptSubMother  = sublorentzsum->Pt();
    fpSubMother   = sublorentzsum->P();
    fySubMother   = sublorentzsum->Rapidity();
    return 1;
  }
  // _________________________________ //      
  else if (Mother == "4LHe") {
    exTrack->Reset(); exTrack1->Reset(); exTrack2->Reset(); exTrack3->Reset(); exTrack4->Reset();
    exTrack1->CopyFromVTrack(track1);
    exTrack2->CopyFromVTrack(track2);
    exTrack3->CopyFromVTrack(track3);
    // _________________________________ //
    // __ get prim vtx __ //
    primVertex = new AliESDVertex(*vertex);
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
    // __ init AliVertexerTracks __ //
    tertvertexer = new AliVertexerTracks(fESDevent->GetMagneticField());
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
    fTertVertChi2 = tertVertex->GetChi2();
    fTertVertNDF  = tertVertex->GetNDF();
    // _________________________________________________ //
    // __ coord of tert vtx __ //
    TertVertex[0] = tertVertex->GetX();
    TertVertex[1] = tertVertex->GetY();
    TertVertex[2] = tertVertex->GetZ();
    fTertVertexX = TertVertex[0];
    fTertVertexY = TertVertex[1];
    fTertVertexZ = TertVertex[2];
    // _________________________________ //
    // __ propagate tracks and get impact parameters __ //
    if (exTrack1->PropagateToDCA(tertVertex, kMagF, 10, impar)) fPropDCADaughter = 1;
    fImParDaughter = impar[0];
    fImParzDaughter = impar[1];
    if (exTrack2->PropagateToDCA(tertVertex, kMagF, 10, impar)) fPropDCADaughter1 = 1;
    fImParDaughter1 = impar[0];
    fImParzDaughter1 = impar[1];
    if (exTrack3->PropagateToDCA(tertVertex, kMagF, 10, impar)) fPropDCADaughter2 = 1;
    fImParDaughter2 = impar[0];
    fImParzDaughter2 = impar[1];
    // _________________________________ //
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
    if (fSubPA < kPointingAngleCut) {      
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      return 0;
    }   
    // _________________________________________________ // 
    // __ cuts from 04/2020 __ //
    fDcaDaughtero  = exTrack1->GetD(fPrimVertexX, fPrimVertexY, kMagF);
    fDcaDaughter1o = exTrack2->GetD(fPrimVertexX, fPrimVertexY, kMagF);
    fDcaDaughter2o = exTrack3->GetD(fPrimVertexX, fPrimVertexY, kMagF);
    // _________________________________________________ //
    // __ track information __ //
    AliAnalysisTaskDoubleHypNucTree::SetDaughterInformation("4LHe", 1);
    // _________________________________________________ //	  	
    // __ daughter hypernucleus information __ //
    fPDGMother    = ksign*fgkPdgCode[kPDGHyperHelium4];
    fChargeMother = ksign*2;
    fRecoMethod   = 1;
    fmSubMother   = sublorentzsum->M();
    fESubMother   = sublorentzsum->E();
    fpxSubMother  = sublorentzsum->Px();
    fpySubMother  = sublorentzsum->Py();
    fpzSubMother  = sublorentzsum->Pz();
    fptSubMother  = sublorentzsum->Pt();
    fpSubMother   = sublorentzsum->P();
    fySubMother   = sublorentzsum->Rapidity();
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
    secvertexer = new AliVertexerTracks(fESDevent->GetMagneticField());
    tertvertexer = new AliVertexerTracks(fESDevent->GetMagneticField());
    secvertexer->SetFinderAlgorithm(kFinderAlgo);
    tertvertexer->SetFinderAlgorithm(kFinderAlgo);
    // _________________________________ //
    // __ get prim vtx __ //
    primVertex = new AliESDVertex(*vertex);
    PrimVertex[0] = primVertex->GetX();
    PrimVertex[1] = primVertex->GetY();
    PrimVertex[2] = primVertex->GetZ();
    // __ coord of prim vtx __ //
    fPrimVertexX = PrimVertex[0];
    fPrimVertexY = PrimVertex[1];
    fPrimVertexZ = PrimVertex[2];
    // __ chi2 of prim vtx __ //
    fPrimVertChi2 = primVertex->GetChi2();
    fPrimVertNDF = primVertex->GetNDF();
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
      TertVertex[0] = tertVertex->GetX();
      TertVertex[1] = tertVertex->GetY();
      TertVertex[2] = tertVertex->GetZ();
      fTertVertexX = TertVertex[0];
      fTertVertexY = TertVertex[1];
      fTertVertexZ = TertVertex[2];
      // __ chi2 of tert vtx __ //
      fTertVertChi2 = tertVertex->GetChi2();
      fTertVertNDF = tertVertex->GetNDF();
      // _________________________________ //
      // __ propagate tracks and get impact parameters __ //
      if (exTrack1->PropagateToDCA(tertVertex, kMagF, 10, impar)) fPropDCADaughter = 1;
      fImParDaughter = impar[0];
      fImParzDaughter = impar[1];
      if (exTrack2->PropagateToDCA(tertVertex, kMagF, 10, impar)) fPropDCADaughter1 = 1;
      fImParDaughter1 = impar[0];
      fImParzDaughter1 = impar[1];
      if (exTrack3->PropagateToDCA(tertVertex, kMagF, 10, impar)) fPropDCADaughter2 = 1;
      fImParDaughter2 = impar[0];
      fImParzDaughter2 = impar[1];
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
      trkArray1->AddAt(exTrack,  0);
      trkArray1->AddAt(exTrack4, 1);
      // __ start at prim vtx __ //
      secvertexer->SetVtxStart(primVertex);
      secVertex = (AliESDVertex*)secvertexer->VertexForSelectedTracks(trkArray1, idTwo, kTRUE, kTRUE, kFALSE);
      if (secvertexer) delete secvertexer;
      if (trkArray1) delete trkArray1;
      // __ coord of sec vtx __ //
      SecVertex[0] = secVertex->GetX();
      SecVertex[1] = secVertex->GetY();
      SecVertex[2] = secVertex->GetZ();
      fSecVertexX = SecVertex[0];
      fSecVertexY = SecVertex[1];
      fSecVertexZ = SecVertex[2];
      // __ chi2 of sec vtx __ //
      fSecVertChi2 = secVertex->GetChi2();
      fSecVertNDF = secVertex->GetNDF();
      // _________________________________ //
      // __ propagate tracks and get impact parameters __ //
      if (exTrack4->PropagateToDCA(secVertex, kMagF, 10, impar)) fPropDCADaughter3 = 1;
      fImParDaughter3 = impar[0];
      fImParzDaughter3 = impar[1];
      if (exTrack->PropagateToDCA(secVertex, kMagF, 10, impar)) fPropDCADaughter4 = 1;
      fImParDaughter4 = impar[0];
      fImParzDaughter4 = impar[1];
      if (secVertex) delete secVertex;
      // _________________________________ //
      // __ dca between tracks __ //
      fDCA2B = TMath::Abs(exTrack4->GetDCA(exTrack, kMagF, xthiss, xpp));
      fDCA3B1 = TMath::Abs(exTrack2->GetDCA(exTrack1, kMagF, xthiss, xpp));
      fDCA3B2 = TMath::Abs(exTrack3->GetDCA(exTrack1, kMagF, xthiss, xpp));
      fDCA3B3 = TMath::Abs(exTrack3->GetDCA(exTrack2, kMagF, xthiss, xpp));
    }
    // _________________________________ //
    else {
      // _________________________________________________ //
      particle1->SetXYZM(2. * exTrack1->Px(), 2. * exTrack1->Py(), 2. * exTrack1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
      particle2->SetXYZM(exTrack2->Px(), exTrack2->Py(), exTrack2->Pz(), AliPID::ParticleMass(AliPID::kProton));
      particle3->SetXYZM(exTrack3->Px(), exTrack3->Py(), exTrack3->Pz(), AliPID::ParticleMass(AliPID::kPion));
      sublorentzsum->SetXYZM(0., 0., 0., 0.);     
      *sublorentzsum = *particle1 + *particle2 + *particle3;
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
      trkArray1->AddAt(exTrack,  0);
      trkArray1->AddAt(exTrack4, 1);
      // __ start at prim vtx __ //
      secvertexer->SetVtxStart(primVertex);
      secVertex = (AliESDVertex*)secvertexer->VertexForSelectedTracks(trkArray1, idTwo, kTRUE, kTRUE, kFALSE);
      if (secvertexer) delete secvertexer;
      if (trkArray1) delete trkArray1;
      // __ coord of sec vtx __ //
      SecVertex[0] = secVertex->GetX();
      SecVertex[1] = secVertex->GetY();
      SecVertex[2] = secVertex->GetZ();
      fSecVertexX = SecVertex[0];
      fSecVertexY = SecVertex[1];
      fSecVertexZ = SecVertex[2];
      // __ chi2 of sec vtx __ //
      fSecVertChi2 = secVertex->GetChi2();
      fSecVertNDF = secVertex->GetNDF();
      // _________________________________ //
      // __ propagate tracks and get impact parameters __ //
      if (exTrack4->PropagateToDCA(secVertex, kMagF, 10, impar)) fPropDCADaughter3 = 1;
      fImParDaughter3 = impar[0];
      fImParzDaughter3 = impar[1];
      if (exTrack->PropagateToDCA(secVertex, kMagF, 10, impar))  fPropDCADaughter4 = 1;
      fImParDaughter4 = impar[0];      
      fImParzDaughter4 = impar[1];
      if (secVertex) delete secVertex;
      // _________________________________ //
      // __ tert vtx __ //
      trkArray = new TObjArray(3);
      trkArray->AddAt(exTrack1, 0);
      trkArray->AddAt(exTrack2, 1);
      trkArray->AddAt(exTrack3, 2);
      // __start at prim vtx __ //
      tertvertexer->SetVtxStart(primVertex);
      tertVertex = (AliESDVertex*)tertvertexer->VertexForSelectedTracks(trkArray, idThree, kTRUE, kTRUE, kFALSE);
      if (trkArray) delete trkArray;
      if (primVertex)   delete primVertex;
      // __ coord of tert vtx __ //
      TertVertex[0] = tertVertex->GetX();
      TertVertex[1] = tertVertex->GetY();
      TertVertex[2] = tertVertex->GetZ();
      fTertVertexX = TertVertex[0];
      fTertVertexY = TertVertex[1];
      fTertVertexZ = TertVertex[2];
      // __ chi2 of tert vtx __ //
      fTertVertChi2 = tertVertex->GetChi2();
      fTertVertNDF = tertVertex->GetNDF();
      // _________________________________ //
      // __ propagate tracks and get impact parameters __ //
      if (exTrack1->PropagateToDCA(tertVertex, kMagF, 10, impar)) fPropDCADaughter = 1;
      fImParDaughter = impar[0];
      fImParzDaughter = impar[1];
      if (exTrack2->PropagateToDCA(tertVertex, kMagF, 10, impar)) fPropDCADaughter1 = 1;
      fImParDaughter1 = impar[0];
      fImParzDaughter1 = impar[1];
      if (exTrack3->PropagateToDCA(tertVertex, kMagF, 10, impar)) fPropDCADaughter2 = 1;
      fImParDaughter2 = impar[0];
      fImParzDaughter2 = impar[1];
      if (tertVertex) delete tertVertex;
      // _________________________________ //
      // __ dca between tracks __ //
      fDCA2B = TMath::Abs(exTrack4->GetDCA(exTrack, kMagF, xthiss, xpp));
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
    if(TMath::Abs(fDCA2B) > kDCATracksCut || TMath::Abs(fDCA3B1) > kDCATracksCut
       || TMath::Abs(fDCA3B2) > kDCATracksCut || TMath::Abs(fDCA3B3) > kDCATracksCut){
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      return 0;
    }
    if (fSubPA2 < kPointingAngleCut) {
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      return 0;
    }
    if (fSubPA  < kPointingAngleCut) {
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      return 0;
    }
    if (fPA     < kPointingAngleCut) {
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      return 0;
    }    
    // _________________________________________________ //
    // __ cuts from 04/2020 __ //
    fDCA2Bo = TMath::Abs(exTrack4->GetD(sublorentzsum->X(), sublorentzsum->Y(), kMagF));
    fDcaDaughtero = exTrack1->GetD(fPrimVertexX, fPrimVertexY, kMagF);
    fDcaDaughter1o = exTrack2->GetD(fPrimVertexX, fPrimVertexY, kMagF);
    fDcaDaughter2o = exTrack3->GetD(fPrimVertexX, fPrimVertexY, kMagF);
    fDcaDaughter3o = exTrack4->GetD(fPrimVertexX, fPrimVertexY, kMagF);
    // _________________________________________________ //
    // __ track information __ //
    AliAnalysisTaskDoubleHypNucTree::SetDaughterInformation("4LLH", kDecayChannel);
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
    // _________________________________________________ //
    // __ mother hypernucleus information __ //
    fPDGMother         = ksign * fgkPdgCode[kPDGDoubleHyperHydrogen4];
    fChargeMother      = ksign;
    fDecayChannel      = kDecayChannel;
    fRecoMethod        = 1;
    fmMother           = lorentzsum->M();
    fmMother2          = lorentzsum2->M();
    fEMother           = lorentzsum->E();
    fpxMother          = lorentzsum->Px();
    fpyMother          = lorentzsum->Py();
    fpzMother          = lorentzsum->Pz();
    fptMother          = lorentzsum->Pt();
    fpMother           = lorentzsum->P();
    fyMother           = lorentzsum->Rapidity();
    // _________________________________________________ //
    // __ daughter hypernucleus information __ //
    fmSubMother        = sublorentzsum->M();
    fESubMother        = sublorentzsum->E();
    fpxSubMother       = sublorentzsum->Px();
    fpySubMother       = sublorentzsum->Py();
    fpzSubMother       = sublorentzsum->Pz();
    fptSubMother       = sublorentzsum->Pt();
    fpSubMother        = sublorentzsum->P();
    fySubMother        = sublorentzsum->Rapidity();
    // __ Check for V0s __ //
    GetV0Status();
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

    fthetaP = TMath::ACos((vecP * vecM) / (vecP.Mag() * vecM.Mag()));
    fthetaN = TMath::ACos((vecN * vecM) / (vecN.Mag() * vecM.Mag()));
    farmalpha = ((vecP.Mag()) * TMath::Cos(fthetaP) - (vecN.Mag()) * TMath::Cos(fthetaN)) / ((vecP.Mag()) * TMath::Cos(fthetaP) + (vecN.Mag()) * TMath::Cos(fthetaN));
    farmpt = vecP.Mag() * sin(fthetaP);
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
    secvertexer = new AliVertexerTracks(fESDevent->GetMagneticField());
    tertvertexer = new AliVertexerTracks(fESDevent->GetMagneticField());
    secvertexer->SetFinderAlgorithm(kFinderAlgo);
    tertvertexer->SetFinderAlgorithm(kFinderAlgo);
    // _________________________________ //
    // __ get prim vtx __ //
    primVertex = new AliESDVertex(*vertex);
    PrimVertex[0] = primVertex->GetX();
    PrimVertex[1] = primVertex->GetY();
    PrimVertex[2] = primVertex->GetZ();
    // __ coord of prim vtx __ //
    fPrimVertexX = PrimVertex[0];
    fPrimVertexY = PrimVertex[1];
    fPrimVertexZ = PrimVertex[2];
    // __ chi2 of prim vtx __ //
    fPrimVertChi2 = primVertex->GetChi2();
    fPrimVertNDF = primVertex->GetNDF();
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
      TertVertex[0] = tertVertex->GetX();
      TertVertex[1] = tertVertex->GetY();
      TertVertex[2] = tertVertex->GetZ();
      fTertVertexX = TertVertex[0];
      fTertVertexY = TertVertex[1];
      fTertVertexZ = TertVertex[2];
      // __ chi2 of tert vtx __ //
      fTertVertChi2 = tertVertex->GetChi2();
      fTertVertNDF = tertVertex->GetNDF();
      // _________________________________ //
      // __ propagate tracks and get impact parameters __ //
      if (exTrack1->PropagateToDCA(tertVertex, kMagF, 10, impar)) fPropDCADaughter = 1;
      fImParDaughter = impar[0];
      fImParzDaughter = impar[1];
      if (exTrack3->PropagateToDCA(tertVertex, kMagF, 10, impar)) fPropDCADaughter2 = 1;
      fImParDaughter2 = impar[0];
      fImParzDaughter2 = impar[1];
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
      SecVertex[0] = secVertex->GetX();
      SecVertex[1] = secVertex->GetY();
      SecVertex[2] = secVertex->GetZ();
      fSecVertexX = SecVertex[0];
      fSecVertexY = SecVertex[1];
      fSecVertexZ = SecVertex[2];
      // __ chi2 of sec vtx __ //
      fSecVertChi2 = secVertex->GetChi2();
      fSecVertNDF = secVertex->GetNDF();
      // _________________________________ //
      // __ propagate tracks and get impact parameters __ //
      if (exTrack2->PropagateToDCA(secVertex, kMagF, 10, impar)) fPropDCADaughter1 = 1;
      fImParDaughter1 = impar[0];
      fImParzDaughter1 = impar[1];
      if (exTrack4->PropagateToDCA(secVertex, kMagF, 10, impar)) fPropDCADaughter3 = 1;
      fImParDaughter3 = impar[0];
      fImParzDaughter3 = impar[1];
      if (exTrack->PropagateToDCA(secVertex, kMagF, 10, impar)) fPropDCADaughter4 = 1;
      fImParDaughter4 = impar[0];
      fImParzDaughter4 = impar[1];
      if (secVertex) delete secVertex;
      // _________________________________ //
      // __ dca between tracks __ //
      fDCA2B  = TMath::Abs(exTrack3->GetDCA(exTrack1, kMagF, xthiss, xpp));
      fDCA3B1 = TMath::Abs(exTrack2->GetDCA(exTrack,  kMagF, xthiss, xpp));
      fDCA3B2 = TMath::Abs(exTrack4->GetDCA(exTrack,  kMagF, xthiss, xpp));
      fDCA3B3 = TMath::Abs(exTrack4->GetDCA(exTrack2, kMagF, xthiss, xpp));
    }
    // _________________________________ //
    else {
      // _________________________________________________ //
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
      trkArray1->AddAt(exTrack,  0);
      trkArray1->AddAt(exTrack2, 1);
      trkArray1->AddAt(exTrack4, 2);
      // __ start at prim vtx __ //
      secvertexer->SetVtxStart(primVertex);
      secVertex = (AliESDVertex*)secvertexer->VertexForSelectedTracks(trkArray1, idThree, kTRUE, kTRUE, kFALSE);
      if (trkArray1) delete trkArray1;
      if (secvertexer) delete secvertexer;
      // __ coord of sec vtx __ //
      SecVertex[0] = secVertex->GetX();
      SecVertex[1] = secVertex->GetY();
      SecVertex[2] = secVertex->GetZ();
      fSecVertexX = SecVertex[0];
      fSecVertexY = SecVertex[1];
      fSecVertexZ = SecVertex[2];
      // __ chi2 of sec vtx __ //
      fSecVertChi2 = secVertex->GetChi2();
      fSecVertNDF = secVertex->GetNDF();
      // _________________________________ //
      // __ propagate tracks and get impact parameters __ //
      if (exTrack2->PropagateToDCA(secVertex, kMagF, 10, impar)) fPropDCADaughter1 = 1;
      fImParDaughter1 = impar[0];
      fImParzDaughter1 = impar[1];
      if (exTrack4->PropagateToDCA(secVertex, kMagF, 10, impar)) fPropDCADaughter3 = 1;
      fImParDaughter3 = impar[0];
      fImParzDaughter3 = impar[1];
      if (exTrack->PropagateToDCA(secVertex, kMagF, 10, impar)) fPropDCADaughter4 = 1;
      fImParDaughter4 = impar[0];
      fImParzDaughter4 = impar[1];
      if(secVertex) delete secVertex;
      // _________________________________ //
      // __ tert vtx __ //
      trkArray = new TObjArray(2);
      trkArray->AddAt(exTrack1, 0);
      trkArray->AddAt(exTrack3, 1);
      // __start at prim vtx __ //
      tertvertexer->SetVtxStart(primVertex);
      tertVertex = (AliESDVertex*)tertvertexer->VertexForSelectedTracks(trkArray, idTwo, kTRUE, kTRUE, kFALSE);
      if (trkArray) delete trkArray;
      if (primVertex)   delete primVertex;
      // __ coord of tert vtx __ //
      TertVertex[0] = tertVertex->GetX();
      TertVertex[1] = tertVertex->GetY();
      TertVertex[2] = tertVertex->GetZ();
      fTertVertexX = TertVertex[0];
      fTertVertexY = TertVertex[1];
      fTertVertexZ = TertVertex[2];
      // __ chi2 of tert vtx __ //
      fTertVertChi2 = tertVertex->GetChi2();
      fTertVertNDF = tertVertex->GetNDF();
      // _________________________________ //
      // __ propagate tracks and get impact parameters __ //
      if (exTrack1->PropagateToDCA(tertVertex, kMagF, 10, impar)) fPropDCADaughter = 1;
      fImParDaughter = impar[0];
      fImParzDaughter = impar[1];
      if (exTrack3->PropagateToDCA(tertVertex, kMagF, 10, impar)) fPropDCADaughter2 = 1;
      fImParDaughter2 = impar[0];
      fImParzDaughter2 = impar[1];
      if (tertVertex) delete tertVertex;
      // _________________________________ //
      // __ dca between tracks __ //
      fDCA2B  = TMath::Abs(exTrack3->GetDCA(exTrack1, kMagF, xthiss, xpp));
      fDCA3B1 = TMath::Abs(exTrack2->GetDCA(exTrack,  kMagF, xthiss, xpp));
      fDCA3B2 = TMath::Abs(exTrack4->GetDCA(exTrack,  kMagF, xthiss, xpp));
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
    if(TMath::Abs(fDCA2B) > kDCATracksCut || TMath::Abs(fDCA3B1) > kDCATracksCut
       || TMath::Abs(fDCA3B2) > kDCATracksCut || TMath::Abs(fDCA3B3) > kDCATracksCut){
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      return 0;
    }
    if (fSubPA2 < kPointingAngleCut) {
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      return 0;
    }
    if (fSubPA  < kPointingAngleCut) {
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      return 0;
    }
    if (fPA     < kPointingAngleCut) {
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      return 0;
    }     
    // _________________________________________________ //
    // __ track information __ //
    AliAnalysisTaskDoubleHypNucTree::SetDaughterInformation("4LLH", kDecayChannel);
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
    // _________________________________________________ //
    // __ mother hypernucleus information __ //
    fPDGMother         = ksign * fgkPdgCode[kPDGDoubleHyperHydrogen4];
    fChargeMother      = ksign;
    fDecayChannel      = kDecayChannel;
    fRecoMethod        = 1;
    fmMother           = lorentzsum->M();
    fmMother2          = lorentzsum2->M();
    fEMother           = lorentzsum->E();
    fpxMother          = lorentzsum->Px();
    fpyMother          = lorentzsum->Py();
    fpzMother          = lorentzsum->Pz();
    fptMother          = lorentzsum->Pt();
    fpMother           = lorentzsum->P();
    fyMother           = lorentzsum->Rapidity();
    // _________________________________________________ //
    // __ daughter hypernucleus information __ //
    fmSubMother        = sublorentzsum->M();
    fESubMother        = sublorentzsum->E();
    fpxSubMother       = sublorentzsum->Px();
    fpySubMother       = sublorentzsum->Py();
    fpzSubMother       = sublorentzsum->Pz();
    fptSubMother       = sublorentzsum->Pt();
    fpSubMother        = sublorentzsum->P();
    fySubMother        = sublorentzsum->Rapidity();
    // __ Check for V0s __ //
    GetV0Status();
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
Int_t AliAnalysisTaskDoubleHypNucTree::KFReconstruction(TString Mother, Int_t kDecayChannel, Int_t ksign) {
  Double_t impar[2];
  Double_t DecLength;
  Double_t DecLengthErr;
  int checkval;
  float l, le;
  // _________________________________________________ //
  // __ get prim vtx __ //
  primVertex = new AliESDVertex(*vertex);
  primKFVertex = CreateKFVertex(*primVertex);
  if(primVertex) delete primVertex;
  PrimVertexKF[0] = primKFVertex.GetX();
  PrimVertexKF[1] = primKFVertex.GetY();
  PrimVertexKF[2] = primKFVertex.GetZ();
  fPrimVertexXKF = PrimVertexKF[0];
  fPrimVertexYKF = PrimVertexKF[1];
  fPrimVertexZKF = PrimVertexKF[2];
  fPrimVertexXErrKF = primKFVertex.GetErrX();
  fPrimVertexYErrKF = primKFVertex.GetErrY();
  fPrimVertexZErrKF = primKFVertex.GetErrZ();
  // __ chi2 of prim vtx __ //
  fPrimVertChi2 = primKFVertex.Chi2();
  fPrimVertNDF = primKFVertex.GetNDF();
  // _________________________________________________ //
  if (Mother == "3LH") {
    // _________________________________________________ //    
    exTrack->Reset(); exTrack1->Reset(); exTrack2->Reset(); exTrack3->Reset(); exTrack4->Reset();
    exTrack1->CopyFromVTrack(track1);
    exTrack3->CopyFromVTrack(track3);
    // __ initialize KFParticles __ //
    KFtrack1 = AliAnalysisTaskDoubleHypNucTree::CreateKFParticle(*exTrack1, AliPID::ParticleMass(AliPID::kHe3),  2);
    KFtrack3 = AliAnalysisTaskDoubleHypNucTree::CreateKFParticle(*exTrack3, AliPID::ParticleMass(AliPID::kPion), 1);
    // __ construct mother particle __ //
    checkval = 1;
    KFParticleHypNuc KFChecktrack3;
    KFChecktrack3.AddDaughter(KFtrack1);
    if (!KFChecktrack3.CheckDaughter(KFtrack3)) checkval = 0;
    KFParticleHypNuc KFChecktrack1;
    if (!KFChecktrack1.CheckDaughter(KFtrack1)) checkval = 0;
    KFChecktrack1.AddDaughter(KFtrack3);
    KFParticleHypNuc KFSubMother;
    if(checkval){
      KFSubMother.SetConstructMethod(2);
      KFSubMother.AddDaughter(KFtrack1);
      KFSubMother.AddDaughter(KFtrack3);
      fRecoMethod = 2;
    } else {
      KFSubMother.AddDaughter(KFtrack1);
      KFSubMother.AddDaughter(KFtrack3);
      fRecoMethod = 0;
    }
    // __ transport mother particle to its decay vertex __ //
    KFSubMother.TransportToDecayVertex();
    // _________________________________________________ //
    // __ get decay vertex position __ //
    TertVertexKF[0] = KFSubMother.GetX();
    TertVertexKF[1] = KFSubMother.GetY();
    TertVertexKF[2] = KFSubMother.GetZ();
    TertVertexErrKF[0] = KFSubMother.GetErrX();
    TertVertexErrKF[1] = KFSubMother.GetErrY();
    TertVertexErrKF[2] = KFSubMother.GetErrZ();
    // __ coord of tert vtx __ //
    fTertVertexXKF = TertVertexKF[0];
    fTertVertexYKF = TertVertexKF[1];
    fTertVertexZKF = TertVertexKF[2];
    fTertVertexXErrKF = TertVertexErrKF[0];
    fTertVertexYErrKF = TertVertexErrKF[1];
    fTertVertexZErrKF = TertVertexErrKF[2];
    // __ vertex quality KF __ //
    fTertVertChi2KF = KFSubMother.Chi2();
    fTertVertNDFKF = KFSubMother.GetNDF();
    // _________________________________________________ //
    particle1->SetXYZM(0.,0.,0.,0.);
    particle3->SetXYZM(0.,0.,0.,0.);
    sublorentzsum->SetXYZM(0.,0.,0.,0.);
    KFtertVertex = new AliESDVertex(TertVertexKF, TertVertexErrKF);
    exTrack1->PropagateToDCA(KFtertVertex, kMagF, 10);
    exTrack3->PropagateToDCA(KFtertVertex, kMagF, 10);
    if(KFtertVertex) delete KFtertVertex;
    particle1->SetXYZM(2.*exTrack1->Px(), 2.*exTrack1->Py(), 2.*exTrack1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
    particle3->SetXYZM(exTrack3->Px(), exTrack3->Py(), exTrack3->Pz(), AliPID::ParticleMass(AliPID::kPion));    
    *sublorentzsum = *particle1 + *particle3;
    fmSubMother = sublorentzsum->M();
    // _________________________________________________ //
    // __ KF daughter hypernucleus cos(PA) to prim vtx __ //
    dd[0] = PrimVertexKF[0] - TertVertexKF[0];
    dd[1] = PrimVertexKF[1] - TertVertexKF[1];
    dd[2] = PrimVertexKF[2] - TertVertexKF[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);  

    if(KFSubMother.GetP() > 0.0){
      sublorentzsum->SetXYZM(0.,0.,0.,0.);
      sublorentzsum->SetXYZM(KFSubMother.GetPx(), KFSubMother.GetPy(), KFSubMother.GetPz(), KFSubMother.GetMass());
    }
    fSubPAKF = TMath::Cos(sublorentzsum->Angle(*h));
    if (fSubPAKF < kPointingAngleCut) {
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      return 0;
    }
    DecLength          = TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2));
    if(KFSubMother.GetMass() > 0.0 && KFSubMother.GetP() > 0.0){
      fctSubMotherKF     = (KFSubMother.GetMass()*DecLength)/KFSubMother.GetP();
      DecLengthErr       = TMath::Sqrt(  (fPrimVertexXErrKF*fPrimVertexXErrKF)/(fPrimVertexXKF*fPrimVertexXKF)
				       + (fPrimVertexYErrKF*fPrimVertexYErrKF)/(fPrimVertexYKF*fPrimVertexYKF)
				       + (fPrimVertexZErrKF*fPrimVertexZErrKF)/(fPrimVertexZKF*fPrimVertexZKF)
				       + (fTertVertexXErrKF*fTertVertexXErrKF)/(fTertVertexXKF*fTertVertexXKF)
				       + (fTertVertexYErrKF*fTertVertexYErrKF)/(fTertVertexYKF*fTertVertexYKF)
				       + (fTertVertexZErrKF*fTertVertexZErrKF)/(fTertVertexZKF*fTertVertexZKF)
				       + (KFSubMother.GetErrP()*KFSubMother.GetErrP())/(KFSubMother.GetP()*KFSubMother.GetP())
				       + (KFSubMother.GetErrMass()*KFSubMother.GetErrMass())/(KFSubMother.GetMass()*KFSubMother.GetMass()))*fctSubMotherKF;
      fctSubMotherErrKF  = DecLengthErr;
    }
    // _________________________________________________ //
    // __ dca tracks KF __ //	
    fDCA3B1XYKF = (Float_t)KFtrack1.GetDistanceFromParticleXY(KFtrack3);
    fDCA3B1ZKF  = (Float_t)KFtrack1.GetDistanceFromParticle(KFtrack3);
    if(TMath::Abs(fDCA3B1XYKF) > kDCATracksCut){
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      return 0;
    }   
    // _________________________________________________ //
    // __ daughter hypernucleus KF information __ //
    fChargeMother = ksign*1;
    fPDGMother = ksign * fgkPdgCode[kPDGHyperHydrogen3];
    fDecayChannel = kDecayChannel;
    fmSubMotherKF = KFSubMother.GetMass();
    fmSubMotherErrKF = KFSubMother.GetErrMass();
    fESubMotherKF = KFSubMother.GetE();
    fESubMotherErrKF = KFSubMother.GetErrE();
    fpxSubMotherKF = KFSubMother.GetPx();
    fpxSubMotherErrKF = KFSubMother.GetErrPx();
    fpySubMotherKF = KFSubMother.GetPy();
    fpySubMotherErrKF = KFSubMother.GetErrPy();
    fpzSubMotherKF = KFSubMother.GetPz();
    fpzSubMotherErrKF = KFSubMother.GetErrPz();
    fptSubMotherKF = KFSubMother.GetPt();
    fptSubMotherErrKF = KFSubMother.GetErrPt();
    fpSubMotherKF = KFSubMother.GetP();
    fpSubMotherErrKF = KFSubMother.GetErrP();
    fySubMotherKF = KFSubMother.GetRapidity();
    // __ dca prim vtx KF __ //
    fDcaDaughterXYKF = (Float_t)KFtrack1.GetDistanceFromVertexXY(primKFVertex);
    fDcaDaughter2XYKF = (Float_t)KFtrack3.GetDistanceFromVertexXY(primKFVertex);
    fDcaDaughterZKF = (Float_t)KFtrack1.GetDistanceFromVertex(primKFVertex);
    fDcaDaughter2ZKF = (Float_t)KFtrack3.GetDistanceFromVertex(primKFVertex);
    // __ dca sec vtx KF __ //
    fDcaSecDaughterXYKF = (Float_t)KFtrack1.GetDistanceFromVertexXY(KFSubMother);
    fDcaSecDaughter2XYKF = (Float_t)KFtrack3.GetDistanceFromVertexXY(KFSubMother);
    fDcaSecDaughterZKF = (Float_t)KFtrack1.GetDistanceFromVertex(KFSubMother);
    fDcaSecDaughter2ZKF = (Float_t)KFtrack3.GetDistanceFromVertex(KFSubMother);    
    // _________________________________________________ //
    AliAnalysisTaskDoubleHypNucTree::SetDaughterInformationKF(KFtrack1, KFtrack2, KFtrack3, kDecayChannel);
    return 1;
  }
  else if (Mother == "4LHe") {
    // _________________________________________________ //    
    exTrack->Reset(); exTrack1->Reset(); exTrack2->Reset(); exTrack3->Reset(); exTrack4->Reset();
    exTrack1->CopyFromVTrack(track1);
    exTrack2->CopyFromVTrack(track2);
    exTrack3->CopyFromVTrack(track3);
    // __ initialize KFParticles __ //
    KFtrack1 = AliAnalysisTaskDoubleHypNucTree::CreateKFParticle(*exTrack1, AliPID::ParticleMass(AliPID::kHe3),    2);
    KFtrack2 = AliAnalysisTaskDoubleHypNucTree::CreateKFParticle(*exTrack2, AliPID::ParticleMass(AliPID::kProton), 1);
    KFtrack3 = AliAnalysisTaskDoubleHypNucTree::CreateKFParticle(*exTrack3, AliPID::ParticleMass(AliPID::kPion),   1);
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
    if(checkval){
      KFSubMother.SetConstructMethod(2);
      KFSubMother.AddDaughter(KFtrack1);
      KFSubMother.AddDaughter(KFtrack2);
      KFSubMother.AddDaughter(KFtrack3);
      fRecoMethod = 2;
    } else {
      KFSubMother.AddDaughter(KFtrack1);
      KFSubMother.AddDaughter(KFtrack2);
      KFSubMother.AddDaughter(KFtrack3);
      fRecoMethod = 0;
    }
    // __ transport mother particle to its decay vertex __ //
    KFSubMother.TransportToDecayVertex();
    // _________________________________________________ //
    // __ get decay vertex position __ //
    TertVertexKF[0] = KFSubMother.GetX();
    TertVertexKF[1] = KFSubMother.GetY();
    TertVertexKF[2] = KFSubMother.GetZ();
    TertVertexErrKF[0] = KFSubMother.GetErrX();
    TertVertexErrKF[1] = KFSubMother.GetErrY();
    TertVertexErrKF[2] = KFSubMother.GetErrZ();
    // __ coord of tert vtx __ //
    fTertVertexXKF = TertVertexKF[0];
    fTertVertexYKF = TertVertexKF[1];
    fTertVertexZKF = TertVertexKF[2];
    fTertVertexXErrKF = TertVertexErrKF[0];
    fTertVertexYErrKF = TertVertexErrKF[1];
    fTertVertexZErrKF = TertVertexErrKF[2];
    // __ vertex quality KF __ //
    fTertVertChi2KF = KFSubMother.Chi2();
    fTertVertNDFKF = KFSubMother.GetNDF();
    // _________________________________________________ //
    particle1->SetXYZM(0.,0.,0.,0.);
    particle2->SetXYZM(0.,0.,0.,0.);
    particle3->SetXYZM(0.,0.,0.,0.);
    sublorentzsum->SetXYZM(0.,0.,0.,0.);
    KFtertVertex = new AliESDVertex(TertVertexKF, TertVertexErrKF);
    exTrack1->PropagateToDCA(KFtertVertex, kMagF, 10);
    exTrack2->PropagateToDCA(KFtertVertex, kMagF, 10);
    exTrack3->PropagateToDCA(KFtertVertex, kMagF, 10);
    if(KFtertVertex) delete KFtertVertex;
    particle1->SetXYZM(2.*exTrack1->Px(), 2.*exTrack1->Py(), 2.*exTrack1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
    particle2->SetXYZM(exTrack2->Px(), exTrack2->Py(), exTrack2->Pz(), AliPID::ParticleMass(AliPID::kProton));
    particle3->SetXYZM(exTrack3->Px(), exTrack3->Py(), exTrack3->Pz(), AliPID::ParticleMass(AliPID::kPion));    
    *sublorentzsum = *particle1 + *particle2 + *particle3;
    fmSubMother = sublorentzsum->M();
    // _________________________________________________ //
    // __ KF daughter hypernucleus cos(PA) to prim vtx __ //
    dd[0] = PrimVertexKF[0] - TertVertexKF[0];
    dd[1] = PrimVertexKF[1] - TertVertexKF[1];
    dd[2] = PrimVertexKF[2] - TertVertexKF[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);   

    if(KFSubMother.GetP() > 0.0){
      sublorentzsum->SetXYZM(0.,0.,0.,0.);
      sublorentzsum->SetXYZM(KFSubMother.GetPx(), KFSubMother.GetPy(), KFSubMother.GetPz(), KFSubMother.GetMass());
    }
    fSubPAKF = TMath::Cos(sublorentzsum->Angle(*h));
    if (fSubPAKF < kPointingAngleCut) {
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      return 0;
    }
    DecLength            = TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2));
    if(KFSubMother.GetMass() > 0.0 && KFSubMother.GetP() > 0.0){
      fctSubMotherKF     = (KFSubMother.GetMass()*DecLength)/KFSubMother.GetP();
      DecLengthErr       = TMath::Sqrt(  (fPrimVertexXErrKF*fPrimVertexXErrKF)/(fPrimVertexXKF*fPrimVertexXKF)
				       + (fPrimVertexYErrKF*fPrimVertexYErrKF)/(fPrimVertexYKF*fPrimVertexYKF)
				       + (fPrimVertexZErrKF*fPrimVertexZErrKF)/(fPrimVertexZKF*fPrimVertexZKF)
				       + (fTertVertexXErrKF*fTertVertexXErrKF)/(fTertVertexXKF*fTertVertexXKF)
				       + (fTertVertexYErrKF*fTertVertexYErrKF)/(fTertVertexYKF*fTertVertexYKF)
				       + (fTertVertexZErrKF*fTertVertexZErrKF)/(fTertVertexZKF*fTertVertexZKF)
				       + (KFSubMother.GetErrP()*KFSubMother.GetErrP())/(KFSubMother.GetP()*KFSubMother.GetP())
				       + (KFSubMother.GetErrMass()*KFSubMother.GetErrMass())/(KFSubMother.GetMass()*KFSubMother.GetMass()))*fctSubMotherKF;
      fctSubMotherErrKF  = DecLengthErr;
    }
    // _________________________________________________ //
    // __ dca tracks KF __ //	
    fDCA3B2XYKF = (Float_t)KFtrack1.GetDistanceFromParticleXY(KFtrack3);
    fDCA3B3XYKF = (Float_t)KFtrack2.GetDistanceFromParticleXY(KFtrack3);
    fDCA3B1XYKF = (Float_t)KFtrack2.GetDistanceFromParticleXY(KFtrack1);
    fDCA3B2ZKF  = (Float_t)KFtrack1.GetDistanceFromParticle(KFtrack3);
    fDCA3B3ZKF  = (Float_t)KFtrack2.GetDistanceFromParticle(KFtrack3);
    fDCA3B1ZKF  = (Float_t)KFtrack2.GetDistanceFromParticle(KFtrack1);
    if(TMath::Abs(fDCA3B1XYKF) > kDCATracksCut || TMath::Abs(fDCA3B2XYKF) > kDCATracksCut
       || TMath::Abs(fDCA3B3XYKF) > kDCATracksCut){
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      return 0;
    }    
    // _________________________________________________ //
    // __ daughter hypernucleus KF information __ //
    fChargeMother = ksign*2;
    fDecayChannel = kDecayChannel;
    fPDGMother = ksign * fgkPdgCode[kPDGHyperHelium4];
    fmSubMotherKF = KFSubMother.GetMass();
    fmSubMotherErrKF = KFSubMother.GetErrMass();
    fESubMotherKF = KFSubMother.GetE();
    fESubMotherErrKF = KFSubMother.GetErrE();
    fpxSubMotherKF = KFSubMother.GetPx();
    fpxSubMotherErrKF = KFSubMother.GetErrPx();
    fpySubMotherKF = KFSubMother.GetPy();
    fpySubMotherErrKF = KFSubMother.GetErrPy();
    fpzSubMotherKF = KFSubMother.GetPz();
    fpzSubMotherErrKF = KFSubMother.GetErrPz();
    fptSubMotherKF = KFSubMother.GetPt();
    fptSubMotherErrKF = KFSubMother.GetErrPt();
    fpSubMotherKF = KFSubMother.GetP();
    fpSubMotherErrKF = KFSubMother.GetErrP();
    fySubMotherKF = KFSubMother.GetRapidity();
    // __ dca prim vtx KF __ //
    fDcaDaughterXYKF = (Float_t)KFtrack1.GetDistanceFromVertexXY(primKFVertex);
    fDcaDaughter1XYKF = (Float_t)KFtrack2.GetDistanceFromVertexXY(primKFVertex);
    fDcaDaughter2XYKF = (Float_t)KFtrack3.GetDistanceFromVertexXY(primKFVertex);
    fDcaDaughterZKF = (Float_t)KFtrack1.GetDistanceFromVertex(primKFVertex);
    fDcaDaughter1ZKF = (Float_t)KFtrack2.GetDistanceFromVertex(primKFVertex);
    fDcaDaughter2ZKF = (Float_t)KFtrack3.GetDistanceFromVertex(primKFVertex);
    // __ dca sec vtx KF __ //
    fDcaSecDaughterXYKF = (Float_t)KFtrack1.GetDistanceFromVertexXY(KFSubMother);
    fDcaSecDaughter1XYKF = (Float_t)KFtrack2.GetDistanceFromVertexXY(KFSubMother);
    fDcaSecDaughter2XYKF = (Float_t)KFtrack3.GetDistanceFromVertexXY(KFSubMother);
    fDcaSecDaughterZKF = (Float_t)KFtrack1.GetDistanceFromVertex(KFSubMother);
    fDcaSecDaughter1ZKF = (Float_t)KFtrack2.GetDistanceFromVertex(KFSubMother);
    fDcaSecDaughter2ZKF = (Float_t)KFtrack3.GetDistanceFromVertex(KFSubMother);    
    // _________________________________________________ //
    AliAnalysisTaskDoubleHypNucTree::SetDaughterInformationKF(KFtrack1, KFtrack2, KFtrack3, kDecayChannel);
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
    KFtrack1 = AliAnalysisTaskDoubleHypNucTree::CreateKFParticle(*exTrack1, AliPID::ParticleMass(AliPID::kHe3),    2);
    KFtrack2 = AliAnalysisTaskDoubleHypNucTree::CreateKFParticle(*exTrack2, AliPID::ParticleMass(AliPID::kProton), 1);
    KFtrack3 = AliAnalysisTaskDoubleHypNucTree::CreateKFParticle(*exTrack3, AliPID::ParticleMass(AliPID::kPion),   1);
    KFtrack4 = AliAnalysisTaskDoubleHypNucTree::CreateKFParticle(*exTrack4, AliPID::ParticleMass(AliPID::kPion),   1);
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
    if(checkval){
      KFSubMother.SetConstructMethod(2);
      KFSubMother.AddDaughter(KFtrack1);
      KFSubMother.AddDaughter(KFtrack2);
      KFSubMother.AddDaughter(KFtrack3);
      fRecoMethod = 2;
    } else {
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
    if(checkval){
      KFMother.SetConstructMethod(2);
      KFMother.AddDaughter(KFSubMother);
      KFMother.AddDaughter(KFtrack4);
      fRecoMethod = 2;
    } else {
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
    SecVertexKF[0] = KFMother.GetX();
    SecVertexKF[1] = KFMother.GetY();
    SecVertexKF[2] = KFMother.GetZ();
    SecVertexErrKF[0] = KFMother.GetErrX();
    SecVertexErrKF[1] = KFMother.GetErrY();
    SecVertexErrKF[2] = KFMother.GetErrZ();
    // __ coord of sec vtx __ //
    fSecVertexXKF = SecVertexKF[0];
    fSecVertexYKF = SecVertexKF[1];
    fSecVertexZKF = SecVertexKF[2];
    fSecVertexXErrKF = SecVertexErrKF[0];
    fSecVertexYErrKF = SecVertexErrKF[1];
    fSecVertexZErrKF = SecVertexErrKF[2];
    // __ vertex quality KF __ //
    fSecVertChi2KF = KFMother.Chi2();
    fSecVertNDFKF = KFMother.GetNDF();
    // __ get decay vertex position __ //
    TertVertexKF[0] = KFSubMother.GetX();
    TertVertexKF[1] = KFSubMother.GetY();
    TertVertexKF[2] = KFSubMother.GetZ();
    TertVertexErrKF[0] = KFSubMother.GetErrX();
    TertVertexErrKF[1] = KFSubMother.GetErrY();
    TertVertexErrKF[2] = KFSubMother.GetErrZ();
    // __ coord of tert vtx __ //
    fTertVertexXKF = TertVertexKF[0];
    fTertVertexYKF = TertVertexKF[1];
    fTertVertexZKF = TertVertexKF[2];
    fTertVertexXErrKF = TertVertexErrKF[0];
    fTertVertexYErrKF = TertVertexErrKF[1];
    fTertVertexZErrKF = TertVertexErrKF[2];
    // __ vertex quality KF __ //
    fTertVertChi2KF = KFSubMother.Chi2();
    fTertVertNDFKF = KFSubMother.GetNDF();
    // _________________________________________________ //
    particle1->SetXYZM(0.,0.,0.,0.);
    particle2->SetXYZM(0.,0.,0.,0.);
    particle3->SetXYZM(0.,0.,0.,0.);
    particle4->SetXYZM(0.,0.,0.,0.);
    sublorentzsum->SetXYZM(0.,0.,0.,0.);
    lorentzsum->SetXYZM(0.,0.,0.,0.);
    KFtertVertex = new AliESDVertex(TertVertexKF, TertVertexErrKF);
    exTrack1->PropagateToDCA(KFtertVertex, kMagF, 10);
    exTrack2->PropagateToDCA(KFtertVertex, kMagF, 10);
    exTrack3->PropagateToDCA(KFtertVertex, kMagF, 10);
    if(KFtertVertex) delete KFtertVertex;
    KFsecVertex = new AliESDVertex(SecVertexKF, SecVertexErrKF);
    exTrack4->PropagateToDCA(KFsecVertex, kMagF, 10);
    if(KFsecVertex) delete KFsecVertex;
    particle1->SetXYZM(2.*exTrack1->Px(), 2.*exTrack1->Py(), 2.*exTrack1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
    particle2->SetXYZM(exTrack2->Px(), exTrack2->Py(), exTrack2->Pz(), AliPID::ParticleMass(AliPID::kProton));
    particle3->SetXYZM(exTrack3->Px(), exTrack3->Py(), exTrack3->Pz(), AliPID::ParticleMass(AliPID::kPion));
    particle4->SetXYZM(exTrack4->Px(), exTrack4->Py(), exTrack4->Pz(), AliPID::ParticleMass(AliPID::kPion));      
    *sublorentzsum = *particle1 + *particle2 + *particle3;    
    *lorentzsum = *sublorentzsum + *particle4;
    fmSubMother = sublorentzsum->M();
    fmMother = lorentzsum->M();
    // _________________________________________________ //
    // __ mother hypernucleus cos(PA) __ //
    dd[0] = PrimVertexKF[0] - SecVertexKF[0];
    dd[1] = PrimVertexKF[1] - SecVertexKF[1];
    dd[2] = PrimVertexKF[2] - SecVertexKF[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    if(KFMother.GetP() > 0.0){
      lorentzsum->SetXYZM(0.,0.,0.,0.);
      lorentzsum->SetXYZM(KFMother.GetPx(), KFMother.GetPy(), KFMother.GetPz(), KFMother.GetMass());
    }
    fPAKF = TMath::Cos(lorentzsum->Angle(*h));
    if(fPAKF < kPointingAngleCut){
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      return 0;
    }
    DecLength            = TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2));
    if(KFMother.GetMass() > 0.0 && KFMother.GetP() > 0.0){
      fctMotherKF        = (KFMother.GetMass()*DecLength)/KFMother.GetP();
      DecLengthErr       = TMath::Sqrt(  (fPrimVertexXErrKF*fPrimVertexXErrKF)/(fPrimVertexXKF*fPrimVertexXKF)
				       + (fPrimVertexYErrKF*fPrimVertexYErrKF)/(fPrimVertexYKF*fPrimVertexYKF)
				       + (fPrimVertexZErrKF*fPrimVertexZErrKF)/(fPrimVertexZKF*fPrimVertexZKF)
				       + (fSecVertexXErrKF*fSecVertexXErrKF)/(fSecVertexXKF*fSecVertexXKF)
				       + (fSecVertexYErrKF*fSecVertexYErrKF)/(fSecVertexYKF*fSecVertexYKF)
				       + (fSecVertexZErrKF*fSecVertexZErrKF)/(fSecVertexZKF*fSecVertexZKF)
				       + (KFMother.GetErrP()*KFMother.GetErrP())/(KFMother.GetP()*KFMother.GetP())
				       + (KFMother.GetErrMass()*KFMother.GetErrMass())/(KFMother.GetMass()*KFMother.GetMass()))*fctMotherKF;
      fctMotherErrKF  = DecLengthErr;
    }
    // _________________________________________________ //    
    // __ daughter hypernucleus cos(PA) to sec vtx __ //
    dd[0] = SecVertexKF[0] - TertVertexKF[0];
    dd[1] = SecVertexKF[1] - TertVertexKF[1];
    dd[2] = SecVertexKF[2] - TertVertexKF[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    if(KFSubMother.GetP() > 0.0){
      sublorentzsum->SetXYZM(0.,0.,0.,0.);
      sublorentzsum->SetXYZM(KFSubMother.GetPx(), KFSubMother.GetPy(), KFSubMother.GetPz(), KFSubMother.GetMass());
    }
    fSubPAKF = TMath::Cos(sublorentzsum->Angle(*h));
    if(fSubPAKF < kPointingAngleCut){
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      return 0;
    }
    DecLength            = TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2));
    if(KFSubMother.GetMass() > 0.0 && KFSubMother.GetP() > 0.0){
      fctSubMotherKF     = (KFSubMother.GetMass()*DecLength)/KFSubMother.GetP();
      DecLengthErr       = TMath::Sqrt(  (fSecVertexXErrKF*fSecVertexXErrKF)/(fSecVertexXKF*fSecVertexXKF)
				       + (fSecVertexYErrKF*fSecVertexYErrKF)/(fSecVertexYKF*fSecVertexYKF)
				       + (fSecVertexZErrKF*fSecVertexZErrKF)/(fSecVertexZKF*fSecVertexZKF)
				       + (fTertVertexXErrKF*fTertVertexXErrKF)/(fTertVertexXKF*fTertVertexXKF)
				       + (fTertVertexYErrKF*fTertVertexYErrKF)/(fTertVertexYKF*fTertVertexYKF)
				       + (fTertVertexZErrKF*fTertVertexZErrKF)/(fTertVertexZKF*fTertVertexZKF)
				       + (KFSubMother.GetErrP()*KFSubMother.GetErrP())/(KFSubMother.GetP()*KFSubMother.GetP())
				       + (KFSubMother.GetErrMass()*KFSubMother.GetErrMass())/(KFSubMother.GetMass()*KFSubMother.GetMass()))*fctSubMotherKF;
      fctSubMotherErrKF  = DecLengthErr;
    }
    // _________________________________________________ //
    // __ KF daughter hypernucleus cos(PA) to prim vtx __ //
    dd[0] = PrimVertexKF[0] - TertVertexKF[0];
    dd[1] = PrimVertexKF[1] - TertVertexKF[1];
    dd[2] = PrimVertexKF[2] - TertVertexKF[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    if(KFSubMother.GetP() > 0.0){
      sublorentzsum->SetXYZM(0.,0.,0.,0.);
      sublorentzsum->SetXYZM(KFSubMother.GetPx(), KFSubMother.GetPy(), KFSubMother.GetPz(), KFSubMother.GetMass());
    }
    fSubPA2KF = TMath::Cos(sublorentzsum->Angle(*h));
    if(fSubPA2KF < kPointingAngleCut){
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      return 0;
    }
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
    if(TMath::Abs(fDCA3B1XYKF) > kDCATracksCut || TMath::Abs(fDCA3B2XYKF) > kDCATracksCut
       || TMath::Abs(fDCA3B3XYKF) > kDCATracksCut || TMath::Abs(fDCA2BXYKF) > kDCATracksCut){
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      return 0;
    }    
    // _________________________________________________ //
    // __ mother hypernucleus KF information __ //
    fChargeMother = ksign*1;
    fDecayChannel = kDecayChannel;
    fPDGMother = ksign * fgkPdgCode[kPDGDoubleHyperHydrogen4];
    fmMotherKF = KFMother.GetMass();
    fmMotherErrKF = KFMother.GetErrMass();
    fEMotherKF = KFMother.GetE();
    fEMotherErrKF = KFMother.GetErrE();
    fpxMotherKF = KFMother.GetPx();
    fpxMotherErrKF = KFMother.GetErrPx();
    fpyMotherKF = KFMother.GetPy();
    fpyMotherErrKF = KFMother.GetErrPy();
    fpzMotherKF = KFMother.GetPz();
    fpzMotherErrKF = KFMother.GetErrPz();
    fptMotherKF = KFMother.GetPt();
    fptMotherErrKF = KFMother.GetErrPt();
    fpMotherKF = KFMother.GetP();
    fpMotherErrKF = KFMother.GetErrP();
    fyMotherKF = KFSubMother.GetRapidity();
    // __ daughter hypernucleus KF information __ //
    fmSubMotherKF = KFSubMother.GetMass();
    fmSubMotherErrKF = KFSubMother.GetErrMass();
    fESubMotherKF = KFSubMother.GetE();
    fESubMotherErrKF = KFSubMother.GetErrE();
    fpxSubMotherKF = KFSubMother.GetPx();
    fpxSubMotherErrKF = KFSubMother.GetErrPx();
    fpySubMotherKF = KFSubMother.GetPy();
    fpySubMotherErrKF = KFSubMother.GetErrPy();
    fpzSubMotherKF = KFSubMother.GetPz();
    fpzSubMotherErrKF = KFSubMother.GetErrPz();
    fptSubMotherKF = KFSubMother.GetPt();
    fptSubMotherErrKF = KFSubMother.GetErrPt();
    fpSubMotherKF = KFSubMother.GetP();
    fpSubMotherErrKF = KFSubMother.GetErrP();
    fySubMotherKF = KFSubMother.GetRapidity();
    // __ Check for V0s __ //
    GetV0Status();
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
    AliAnalysisTaskDoubleHypNucTree::SetDaughterInformationKF(KFtrack1, KFtrack2, KFtrack3, KFtrack4, kDecayChannel);
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

    fthetaP = TMath::ACos((vecP * vecM) / (vecP.Mag() * vecM.Mag()));
    fthetaN = TMath::ACos((vecN * vecM) / (vecN.Mag() * vecM.Mag()));
    farmalpha = ((vecP.Mag()) * TMath::Cos(fthetaP) - (vecN.Mag()) * TMath::Cos(fthetaN)) / ((vecP.Mag()) * TMath::Cos(fthetaP) + (vecN.Mag()) * TMath::Cos(fthetaN));
    farmpt = vecP.Mag() * sin(fthetaP);
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
    KFtrack1 = AliAnalysisTaskDoubleHypNucTree::CreateKFParticle(*exTrack1, AliPID::ParticleMass(AliPID::kHe3),    2);
    KFtrack2 = AliAnalysisTaskDoubleHypNucTree::CreateKFParticle(*exTrack2, AliPID::ParticleMass(AliPID::kProton), 1);
    KFtrack3 = AliAnalysisTaskDoubleHypNucTree::CreateKFParticle(*exTrack3, AliPID::ParticleMass(AliPID::kPion),   1);
    KFtrack4 = AliAnalysisTaskDoubleHypNucTree::CreateKFParticle(*exTrack4, AliPID::ParticleMass(AliPID::kPion),   1);
    // __ construct mother particle __ //
    checkval = 1;
    KFParticleHypNuc KFChecktrack1;
    if (!KFChecktrack1.CheckDaughter(KFtrack1)) checkval = 0;
    KFChecktrack1.AddDaughter(KFtrack3);
    KFParticleHypNuc KFChecktrack4;
    if (!KFChecktrack4.CheckDaughter(KFtrack3)) checkval = 0;
    KFChecktrack4.AddDaughter(KFtrack1);

    KFParticleHypNuc KFSubMother;
    if(checkval){
      KFSubMother.SetConstructMethod(2);
      KFSubMother.AddDaughter(KFtrack1);
      KFSubMother.AddDaughter(KFtrack3);
      fRecoMethod = 2;
    } else {
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
    if(checkval){
      KFMother.SetConstructMethod(2);
      KFMother.AddDaughter(KFSubMother);
      KFMother.AddDaughter(KFtrack2);
      KFMother.AddDaughter(KFtrack4);
      fRecoMethod = 2;
    } else{
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
    SecVertexKF[0] = KFMother.GetX();
    SecVertexKF[1] = KFMother.GetY();
    SecVertexKF[2] = KFMother.GetZ();
    SecVertexErrKF[0] = KFMother.GetErrX();
    SecVertexErrKF[1] = KFMother.GetErrY();
    SecVertexErrKF[2] = KFMother.GetErrZ();
    // __ coord of sec vtx __ //
    fSecVertexXKF = SecVertexKF[0];
    fSecVertexYKF = SecVertexKF[1];
    fSecVertexZKF = SecVertexKF[2];
    fSecVertexXErrKF = SecVertexErrKF[0];
    fSecVertexYErrKF = SecVertexErrKF[1];
    fSecVertexZErrKF = SecVertexErrKF[2];
    // __ vertex quality KF __ //
    fSecVertChi2KF = KFMother.Chi2();
    fSecVertNDFKF = KFMother.GetNDF();
    // __ get decay vertex position __ //
    TertVertexKF[0] = KFSubMother.GetX();
    TertVertexKF[1] = KFSubMother.GetY();
    TertVertexKF[2] = KFSubMother.GetZ();
    TertVertexErrKF[0] = KFSubMother.GetErrX();
    TertVertexErrKF[1] = KFSubMother.GetErrY();
    TertVertexErrKF[2] = KFSubMother.GetErrZ();
    // __ coord of tert vtx __ //
    fTertVertexXKF = TertVertexKF[0];
    fTertVertexYKF = TertVertexKF[1];
    fTertVertexZKF = TertVertexKF[2];
    fTertVertexXErrKF = TertVertexErrKF[0];
    fTertVertexYErrKF = TertVertexErrKF[1];
    fTertVertexZErrKF = TertVertexErrKF[2];
    // __ vertex quality KF __ //
    fTertVertChi2KF = KFSubMother.Chi2();
    fTertVertNDFKF = KFSubMother.GetNDF();
    // _________________________________________________ //
    particle1->SetXYZM(0.,0.,0.,0.);
    particle2->SetXYZM(0.,0.,0.,0.);
    particle3->SetXYZM(0.,0.,0.,0.);
    particle4->SetXYZM(0.,0.,0.,0.);
    sublorentzsum->SetXYZM(0.,0.,0.,0.);
    lorentzsum->SetXYZM(0.,0.,0.,0.);
    KFtertVertex = new AliESDVertex(TertVertexKF, TertVertexErrKF);
    exTrack1->PropagateToDCA(KFtertVertex, kMagF, 10);    
    exTrack3->PropagateToDCA(KFtertVertex, kMagF, 10);
    if(KFtertVertex) delete KFtertVertex;
    KFsecVertex = new AliESDVertex(SecVertexKF, SecVertexErrKF);
    exTrack2->PropagateToDCA(KFsecVertex, kMagF, 10);
    exTrack4->PropagateToDCA(KFsecVertex, kMagF, 10);
    if(KFsecVertex) delete KFsecVertex;
    particle1->SetXYZM(2.*exTrack1->Px(), 2.*exTrack1->Py(), 2.*exTrack1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
    particle2->SetXYZM(exTrack2->Px(), exTrack2->Py(), exTrack2->Pz(), AliPID::ParticleMass(AliPID::kProton));
    particle3->SetXYZM(exTrack3->Px(), exTrack3->Py(), exTrack3->Pz(), AliPID::ParticleMass(AliPID::kPion));
    particle4->SetXYZM(exTrack4->Px(), exTrack4->Py(), exTrack4->Pz(), AliPID::ParticleMass(AliPID::kPion));      
    *sublorentzsum = *particle1  + *particle3;    
    *lorentzsum = *sublorentzsum + *particle4 + *particle2;
    fmSubMother = sublorentzsum->M();
    fmMother = lorentzsum->M();
    // _________________________________________________ //
    // __ mother hypernucleus cos(PA) __ //
    dd[0] = PrimVertexKF[0] - SecVertexKF[0];
    dd[1] = PrimVertexKF[1] - SecVertexKF[1];
    dd[2] = PrimVertexKF[2] - SecVertexKF[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    if(KFMother.GetP()){
      lorentzsum->SetXYZM(0.,0.,0.,0.);
      lorentzsum->SetXYZM(KFMother.GetPx(), KFMother.GetPy(), KFMother.GetPz(), KFMother.GetMass());
    }
    fPAKF = TMath::Cos(lorentzsum->Angle(*h));
    if(fPAKF < kPointingAngleCut){
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      return 0;
    }
    DecLength      = TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2));
    if(KFMother.GetMass() > 0.0 && KFMother.GetP() > 0.0){
      fctMotherKF    = (KFMother.GetMass()*DecLength)/KFMother.GetP();
      DecLengthErr   = TMath::Sqrt(  (fPrimVertexXErrKF*fPrimVertexXErrKF)/(fPrimVertexXKF*fPrimVertexXKF)
				   + (fPrimVertexYErrKF*fPrimVertexYErrKF)/(fPrimVertexYKF*fPrimVertexYKF)
				   + (fPrimVertexZErrKF*fPrimVertexZErrKF)/(fPrimVertexZKF*fPrimVertexZKF)
				   + (fSecVertexXErrKF*fSecVertexXErrKF)/(fSecVertexXKF*fSecVertexXKF)
				   + (fSecVertexYErrKF*fSecVertexYErrKF)/(fSecVertexYKF*fSecVertexYKF)
				   + (fSecVertexZErrKF*fSecVertexZErrKF)/(fSecVertexZKF*fSecVertexZKF)
				   + (KFMother.GetErrP()*KFMother.GetErrP())/(KFMother.GetP()*KFMother.GetP())
				   + (KFMother.GetErrMass()*KFMother.GetErrMass())/(KFMother.GetMass()*KFMother.GetMass()))*fctMotherKF;
      fctMotherErrKF  = DecLengthErr;
    }
    // _________________________________________________ //
    // __ daughter hypernucleus cos(PA) to sec vtx __ //
    dd[0] = SecVertexKF[0] - TertVertexKF[0];
    dd[1] = SecVertexKF[1] - TertVertexKF[1];
    dd[2] = SecVertexKF[2] - TertVertexKF[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    if(KFSubMother.GetP() > 0.0){
      sublorentzsum->SetXYZM(0.,0.,0.,0.);
      sublorentzsum->SetXYZM(KFSubMother.GetPx(), KFSubMother.GetPy(), KFSubMother.GetPz(), KFSubMother.GetMass());
    }
    fSubPAKF = TMath::Cos(sublorentzsum->Angle(*h));
    if(fSubPAKF < kPointingAngleCut){
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      return 0;
    }
    DecLength = TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2));
    if(KFSubMother.GetMass() > 0.0 && KFSubMother.GetP() > 0.0){
      fctSubMotherKF     = (KFSubMother.GetMass()*DecLength)/KFSubMother.GetP();
      DecLengthErr       = TMath::Sqrt(  (fSecVertexXErrKF*fSecVertexXErrKF)/(fSecVertexXKF*fSecVertexXKF)
				       + (fSecVertexYErrKF*fSecVertexYErrKF)/(fSecVertexYKF*fSecVertexYKF)
				       + (fSecVertexZErrKF*fSecVertexZErrKF)/(fSecVertexZKF*fSecVertexZKF)
				       + (fTertVertexXErrKF*fTertVertexXErrKF)/(fTertVertexXKF*fTertVertexXKF)
				       + (fTertVertexYErrKF*fTertVertexYErrKF)/(fTertVertexYKF*fTertVertexYKF)
				       + (fTertVertexZErrKF*fTertVertexZErrKF)/(fTertVertexZKF*fTertVertexZKF)
				       + (KFSubMother.GetErrP()*KFSubMother.GetErrP())/(KFSubMother.GetP()*KFSubMother.GetP())
				       + (KFSubMother.GetErrMass()*KFSubMother.GetErrMass())/(KFSubMother.GetMass()*KFSubMother.GetMass()))*fctSubMotherKF;
      fctSubMotherErrKF  = DecLengthErr;
    }
    // _________________________________________________ //
    // __ KF daughter hypernucleus cos(PA) to prim vtx __ //
    dd[0] = PrimVertexKF[0] - TertVertexKF[0];
    dd[1] = PrimVertexKF[1] - TertVertexKF[1];
    dd[2] = PrimVertexKF[2] - TertVertexKF[2];
    h->SetXYZ(-dd[0], -dd[1], -dd[2]);

    if(KFSubMother.GetP() > 0.0){
      sublorentzsum->SetXYZM(0.,0.,0.,0.);
      sublorentzsum->SetXYZM(KFSubMother.GetPx(), KFSubMother.GetPy(), KFSubMother.GetPz(), KFSubMother.GetMass());
    }
    fSubPA2KF = TMath::Cos(sublorentzsum->Angle(*h));
    if(fSubPA2KF < kPointingAngleCut){
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      return 0;
    }
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
    if(TMath::Abs(fDCA3B1XYKF) > kDCATracksCut || TMath::Abs(fDCA3B2XYKF) > kDCATracksCut
       || TMath::Abs(fDCA3B3XYKF) > kDCATracksCut || TMath::Abs(fDCA2BXYKF) > kDCATracksCut){
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      return 0;
    }    
    // _________________________________________________ //
    // __ mother hypernucleus KF information __ //
    fChargeMother = ksign*1;
    fDecayChannel = kDecayChannel;
    fPDGMother = ksign * fgkPdgCode[kPDGDoubleHyperHydrogen4];
    fmMotherKF = KFMother.GetMass();
    fmMotherErrKF = KFMother.GetErrMass();
    fEMotherKF = KFMother.GetE();
    fEMotherErrKF = KFMother.GetErrE();
    fpxMotherKF = KFMother.GetPx();
    fpxMotherErrKF = KFMother.GetErrPx();
    fpyMotherKF = KFMother.GetPy();
    fpyMotherErrKF = KFMother.GetErrPy();
    fpzMotherKF = KFMother.GetPz();
    fpzMotherErrKF = KFMother.GetErrPz();
    fptMotherKF = KFMother.GetPt();
    fptMotherErrKF = KFMother.GetErrPt();
    fpMotherKF = KFMother.GetP();
    fpMotherErrKF = KFMother.GetErrP();
    fyMotherKF = KFSubMother.GetRapidity();
    // __ daughter hypernucleus KF information __ //
    fmSubMotherKF = KFSubMother.GetMass();
    fmSubMotherErrKF = KFSubMother.GetErrMass();
    fESubMotherKF = KFSubMother.GetE();
    fESubMotherErrKF = KFSubMother.GetErrE();
    fpxSubMotherKF = KFSubMother.GetPx();
    fpxSubMotherErrKF = KFSubMother.GetErrPx();
    fpySubMotherKF = KFSubMother.GetPy();
    fpySubMotherErrKF = KFSubMother.GetErrPy();
    fpzSubMotherKF = KFSubMother.GetPz();
    fpzSubMotherErrKF = KFSubMother.GetErrPz();
    fptSubMotherKF = KFSubMother.GetPt();
    fptSubMotherErrKF = KFSubMother.GetErrPt();
    fpSubMotherKF = KFSubMother.GetP();
    fpSubMotherErrKF = KFSubMother.GetErrP();
    fySubMotherKF = KFSubMother.GetRapidity();
    // __ Check for V0s __ //
    GetV0Status();
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
    AliAnalysisTaskDoubleHypNucTree::SetDaughterInformationKF(KFtrack1, KFtrack2, KFtrack3, KFtrack4, kDecayChannel);
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

    fthetaP = TMath::ACos((vecP * vecM) / (vecP.Mag() * vecM.Mag()));
    fthetaN = TMath::ACos((vecN * vecM) / (vecN.Mag() * vecM.Mag()));
    farmalpha = ((vecP.Mag()) * TMath::Cos(fthetaP) - (vecN.Mag()) * TMath::Cos(fthetaN)) / ((vecP.Mag()) * TMath::Cos(fthetaP) + (vecN.Mag()) * TMath::Cos(fthetaN));
    farmpt = vecP.Mag() * sin(fthetaP);
    return 1;
  }
  else return 0;
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::GetV0Status() {
  for (Int_t ivertex = 0; ivertex < fESDevent->GetNumberOfV0s(); ivertex++) {    

    //Get V0
    AliESDv0 *fV0 = fESDevent->GetV0(ivertex);

    //check if on the fly
    if(!fV0->GetOnFlyStatus()) continue;

    //ChargeCorrection
    Bool_t v0ChargeCorrect = kTRUE;        
    AliESDtrack* trackN = fESDevent->GetTrack(fV0->GetIndex(0));//pos Track
    AliESDtrack* trackP = fESDevent->GetTrack(fV0->GetIndex(1));//neg Track
    if (trackN->GetSign() > 0 ) {
      trackN = fESDevent->GetTrack(fV0->GetIndex(1));//neg Track
      trackP = fESDevent->GetTrack(fV0->GetIndex(0));//pos Track
      v0ChargeCorrect = kFALSE;
    }
    if(track1->GetSign() > 0){
      if(trackP->GetLabel() == track1->GetLabel() && trackN->GetLabel() == track3->GetLabel()){
	fV0VertexX_13 = fV0->Xv();
	fV0VertexY_13 = fV0->Yv();
	fV0VertexZ_13 = fV0->Zv();
	fisOnlineV0_13 = 1;
      }
      if(trackP->GetLabel() == track1->GetLabel() && trackN->GetLabel() == track4->GetLabel()){
	fV0VertexX_14 = fV0->Xv();
	fV0VertexY_14 = fV0->Yv();
	fV0VertexZ_14 = fV0->Zv();
	fisOnlineV0_14 = 1;
      }
      if(trackP->GetLabel() == track2->GetLabel() && trackN->GetLabel() == track3->GetLabel()){
	fV0VertexX_23 = fV0->Xv();
	fV0VertexY_23 = fV0->Yv();
	fV0VertexZ_23 = fV0->Zv();
	fisOnlineV0_23 = 1;
      }
      if(trackP->GetLabel() == track2->GetLabel() && trackN->GetLabel() == track4->GetLabel()){
	fV0VertexX_24 = fV0->Xv();
	fV0VertexY_24 = fV0->Yv();
	fV0VertexZ_24 = fV0->Zv();
	fisOnlineV0_24 = 1;
      }
    }
    if(track1->GetSign() < 0){
      if(trackN->GetLabel() == track1->GetLabel() && trackP->GetLabel() == track3->GetLabel()){
	fV0VertexX_13 = fV0->Xv();
	fV0VertexY_13 = fV0->Yv();
	fV0VertexZ_13 = fV0->Zv();
	fisOnlineV0_13 = 1;
      }
      if(trackN->GetLabel() == track1->GetLabel() && trackP->GetLabel() == track4->GetLabel()){
	fV0VertexX_14 = fV0->Xv();
	fV0VertexY_14 = fV0->Yv();
	fV0VertexZ_14 = fV0->Zv();
	fisOnlineV0_14 = 1;
      }
      if(trackN->GetLabel() == track2->GetLabel() && trackP->GetLabel() == track3->GetLabel()){
	fV0VertexX_23 = fV0->Xv();
	fV0VertexY_23 = fV0->Yv();
	fV0VertexZ_23 = fV0->Zv();
	fisOnlineV0_23 = 1;
      }
      if(trackN->GetLabel() == track2->GetLabel() && trackP->GetLabel() == track4->GetLabel()){
	fV0VertexX_24 = fV0->Xv();
	fV0VertexY_24 = fV0->Yv();
	fV0VertexZ_24 = fV0->Zv();
	fisOnlineV0_24 = 1;
      }
    }
  }
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::SetDaughterInformation(TString Mother, Int_t kDecayChannel) {

  Float_t xv[2], yv[3];
  if (Mother == "4LHe" || Mother == "4LLH" || Mother == "3LH") {
    // _________________________________ //
    // __ He3 track __ //
    fpDaughter = track1->GetInnerParam()->GetP();
    fpxDaughter = particle1->Px();
    fpyDaughter = particle1->Py();
    fpzDaughter = particle1->Pz();
    fptDaughter = particle1->Pt();
    fEDaughter = particle1->E();
    fyDaughter = particle1->Rapidity();
    fdEdxDaughter = track1->GetTPCsignal();
    fdEdxSigmaDaughter = AliAnalysisTaskDoubleHypNucTree::Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
    if (fMCtrue) fdEdxSigmaDaughter = fPID->NumberOfSigmasTPC(track1, AliPID::kHe3);
    track1->GetImpactParameters(xv, yv);
    fDcaDaughter = xv[0];
    fDcazDaughter = xv[1];
    fSigmaYXDaughter = yv[0];
    fSigmaXYZDaughter = yv[1];
    fSigmaZDaughter = yv[2];
    fDcaSecDaughter = TMath::Abs(exTrack1->GetD(TertVertex[0], TertVertex[1], kMagF));
    fNclsDaughter = track1->GetTPCNcls();
    fNclsITSDaughter = track1->GetITSNcls();
    fChi2Daughter = track1->GetTPCchi2() / (Float_t)track1->GetTPCclusters(0);
    fEtaDaughter = track1->Eta();
    fPhiDaughter = track1->Phi();
    fGeoLengthDaughter = AliAnalysisTaskDoubleHypNucTree::GeoLength(*track1);
    fTOFSignalDaughter = AliAnalysisTaskDoubleHypNucTree::GetTOFSignal(*track1);
    fPtUncertDaughter = TMath::Sqrt(exTrack1->GetSigma1Pt2()) * fptDaughter;
    fTPCRefitDaughter = (track1->GetStatus() & AliESDtrack::kTPCrefit) != 0;
    fITSRefitDaughter = (track1->GetStatus() & AliESDtrack::kITSrefit) != 0;
    fdEdxSigmaTriton = AliAnalysisTaskDoubleHypNucTree::Bethe(*track1, AliPID::ParticleMass(AliPID::kTriton), 1, fBetheParamsT);
    fdEdxSigmaAlpha = AliAnalysisTaskDoubleHypNucTree::Bethe(*track1, AliPID::ParticleMass(AliPID::kAlpha), 2, fBetheParamsHe);
    fITSLayer1Daughter = track1->HasPointOnITSLayer(0);
    fITSLayer2Daughter = track1->HasPointOnITSLayer(1);
    fITSLayer3Daughter = track1->HasPointOnITSLayer(2);
    fITSLayer4Daughter = track1->HasPointOnITSLayer(3);
    fITSLayer5Daughter = track1->HasPointOnITSLayer(4);
    fITSLayer6Daughter = track1->HasPointOnITSLayer(5);
    // _________________________________ //
    if(Mother != "3LH"){
      // __ p track __ //
      fpDaughter1 = track2->GetInnerParam()->GetP();
      fpxDaughter1 = particle2->Px();
      fpyDaughter1 = particle2->Py();
      fpzDaughter1 = particle2->Pz();
      fptDaughter1 = particle2->Pt();
      fEDaughter1 = particle2->E();
      fyDaughter1 = particle2->Rapidity();
      fdEdxDaughter1 = track2->GetTPCsignal();
      //fdEdxSigmaDaughter1 = AliAnalysisTaskDoubleHypNucTree::Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
      fdEdxSigmaDaughter1 = fPID->NumberOfSigmasTPC(track2, AliPID::kProton);
      //if (fMCtrue) fdEdxSigmaDaughter1 = fPID->NumberOfSigmasTPC(track2, AliPID::kProton);
      track2->GetImpactParameters(xv, yv);
      fDcaDaughter1 = xv[0];
      fDcazDaughter1 = xv[1];
      fSigmaYXDaughter1 = yv[0];
      fSigmaXYZDaughter1 = yv[1];
      fSigmaZDaughter1 = yv[2];
      if (kDecayChannel == 1) fDcaSecDaughter1 = TMath::Abs(exTrack2->GetD(TertVertex[0], TertVertex[1], kMagF));
      else                   fDcaSecDaughter1 = TMath::Abs(exTrack2->GetD(SecVertex[0], SecVertex[1], kMagF));
      fNclsDaughter1 = track2->GetTPCNcls();
      fNclsITSDaughter1 = track2->GetITSNcls();
      fChi2Daughter1 = track2->GetTPCchi2() / (Float_t)track2->GetTPCclusters(0);
      fEtaDaughter1 = track2->Eta();
      fPhiDaughter1 = track2->Phi();
      fGeoLengthDaughter1 = AliAnalysisTaskDoubleHypNucTree::GeoLength(*track2);
      fTOFSignalDaughter1 = AliAnalysisTaskDoubleHypNucTree::GetTOFSignal(*track2);
      fPtUncertDaughter1 = TMath::Sqrt(exTrack2->GetSigma1Pt2()) * fptDaughter1;
      fTPCRefitDaughter1 = (track2->GetStatus() & AliESDtrack::kTPCrefit) != 0;
      fITSRefitDaughter1 = (track2->GetStatus() & AliESDtrack::kITSrefit) != 0;
      fdEdxSigmaPion = fPID->NumberOfSigmasTPC(track2, AliPID::kPion);
      fdEdxSigmaDeuteron = AliAnalysisTaskDoubleHypNucTree::Bethe(*track2, AliPID::ParticleMass(AliPID::kDeuteron), 1, fBetheParamsT);
      fITSLayer1Daughter1 = track2->HasPointOnITSLayer(0);
      fITSLayer2Daughter1 = track2->HasPointOnITSLayer(1);
      fITSLayer3Daughter1 = track2->HasPointOnITSLayer(2);
      fITSLayer4Daughter1 = track2->HasPointOnITSLayer(3);
      fITSLayer5Daughter1 = track2->HasPointOnITSLayer(4);
      fITSLayer6Daughter1 = track2->HasPointOnITSLayer(5);
    }
    // _________________________________ //
    // __ pion track __ //
    fpDaughter2 = track3->GetInnerParam()->GetP();
    fpxDaughter2 = particle3->Px();
    fpyDaughter2 = particle3->Py();
    fpzDaughter2 = particle3->Pz();
    fptDaughter2 = particle3->Pt();
    fEDaughter2 = particle3->E();
    fyDaughter2 = particle3->Rapidity();
    fdEdxDaughter2 = track3->GetTPCsignal();
    fdEdxSigmaDaughter2 = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
    track3->GetImpactParameters(xv, yv);
    fDcaDaughter2 = xv[0];
    fDcazDaughter2 = xv[1];
    fSigmaYXDaughter2 = yv[0];
    fSigmaXYZDaughter2 = yv[1];
    fSigmaZDaughter2 = yv[2];
    if (kDecayChannel == 1) fDcaSecDaughter2 = TMath::Abs(exTrack3->GetD(TertVertex[0], TertVertex[1], kMagF));
    else                   fDcaSecDaughter2 = TMath::Abs(exTrack3->GetD(SecVertex[0], SecVertex[1], kMagF));
    fNclsDaughter2 = track3->GetTPCNcls();
    fNclsITSDaughter2 = track3->GetITSNcls();
    fChi2Daughter2 = track3->GetTPCchi2() / (Float_t)track3->GetTPCclusters(0);
    fEtaDaughter2 = track3->Eta();
    fPhiDaughter2 = track3->Phi();
    fGeoLengthDaughter2 = AliAnalysisTaskDoubleHypNucTree::GeoLength(*track3);
    fTOFSignalDaughter2 = AliAnalysisTaskDoubleHypNucTree::GetTOFSignal(*track3);
    fPtUncertDaughter2 = TMath::Sqrt(exTrack3->GetSigma1Pt2()) * fptDaughter2;
    fTPCRefitDaughter2 = (track3->GetStatus() & AliESDtrack::kTPCrefit) != 0;
    fITSRefitDaughter2 = (track3->GetStatus() & AliESDtrack::kITSrefit) != 0;
    fITSLayer1Daughter2 = track3->HasPointOnITSLayer(0);
    fITSLayer2Daughter2 = track3->HasPointOnITSLayer(1);
    fITSLayer3Daughter2 = track3->HasPointOnITSLayer(2);
    fITSLayer4Daughter2 = track3->HasPointOnITSLayer(3);
    fITSLayer5Daughter2 = track3->HasPointOnITSLayer(4);
    fITSLayer6Daughter2 = track3->HasPointOnITSLayer(5);
  }
  if (Mother == "4LLH") {
    // _________________________________ //
    // __ pion track __ //
    fpDaughter3 = track4->GetInnerParam()->GetP();
    fpxDaughter3 = particle4->Px();
    fpyDaughter3 = particle4->Py();
    fpzDaughter3 = particle4->Pz();
    fptDaughter3 = particle4->Pt();
    fEDaughter3 = particle4->E();
    fyDaughter3 = particle4->Rapidity();
    fdEdxDaughter3 = track4->GetTPCsignal();
    fdEdxSigmaDaughter3 = fPID->NumberOfSigmasTPC(track4, AliPID::kPion);
    track4->GetImpactParameters(xv, yv);
    fDcaDaughter3 = xv[0];
    fDcazDaughter3 = xv[1];
    fSigmaYXDaughter3 = yv[0];
    fSigmaXYZDaughter3 = yv[1];
    fSigmaZDaughter3 = yv[2];
    if (kDecayChannel == 1) fDcaSecDaughter3 = TMath::Abs(exTrack4->GetD(SecVertex[0], SecVertex[1], kMagF));
    else                   fDcaSecDaughter3 = TMath::Abs(exTrack4->GetD(TertVertex[0], TertVertex[1], kMagF));
    fNclsDaughter3 = track4->GetTPCNcls();
    fNclsITSDaughter3 = track4->GetITSNcls();
    fChi2Daughter3 = track4->GetTPCchi2() / (Float_t)track4->GetTPCclusters(0);
    fEtaDaughter3 = track4->Eta();
    fPhiDaughter3 = track4->Phi();
    fGeoLengthDaughter3 = AliAnalysisTaskDoubleHypNucTree::GeoLength(*track4);
    fTOFSignalDaughter3 = AliAnalysisTaskDoubleHypNucTree::GetTOFSignal(*track4);
    fPtUncertDaughter3 = TMath::Sqrt(exTrack4->GetSigma1Pt2()) * fptDaughter3;
    fTPCRefitDaughter3 = (track4->GetStatus() & AliESDtrack::kTPCrefit) != 0;
    fITSRefitDaughter3 = (track4->GetStatus() & AliESDtrack::kITSrefit) != 0;
    fITSLayer1Daughter3 = track4->HasPointOnITSLayer(0);
    fITSLayer2Daughter3 = track4->HasPointOnITSLayer(1);
    fITSLayer3Daughter3 = track4->HasPointOnITSLayer(2);
    fITSLayer4Daughter3 = track4->HasPointOnITSLayer(3);
    fITSLayer5Daughter3 = track4->HasPointOnITSLayer(4);
    fITSLayer6Daughter3 = track4->HasPointOnITSLayer(5);
  }
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::SetDaughterInformationKF(KFParticle daughter1, KFParticle daughter2, KFParticle daughter3, Int_t kDecayChannel) {

  Float_t xv[2], yv[3];
  AliMCParticle *part;  
  // _________________________________ //
  // __ He3 track __ //
  fpDaughter = track1->GetInnerParam()->GetP();
  if(fMCtrue) {
    part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(track1->GetLabel()))->Particle());
    fLabelDaughterKF = TMath::Abs(part->PdgCode()) == fgkPdgCode[kPDGHelium3];
  }
  fpxDaughterKF = daughter1.GetPx();
  fpyDaughterKF = daughter1.GetPy();
  fpzDaughterKF = daughter1.GetPz();
  fptDaughterKF = daughter1.GetPt();
  fEDaughterKF = daughter1.GetE();
  fyDaughterKF = daughter1.GetRapidity();
  fpxDaughter = particle1->Px();
  fpyDaughter = particle1->Py();
  fpzDaughter = particle1->Pz();
  fptDaughter = particle1->Pt();
  fEDaughter = particle1->E();
  fyDaughter = particle1->Rapidity();
  fdEdxDaughter = track1->GetTPCsignal();
  fdEdxSigmaDaughter = AliAnalysisTaskDoubleHypNucTree::Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
  if (fMCtrue) fdEdxSigmaDaughter = fPID->NumberOfSigmasTPC(track1, AliPID::kHe3);
  track1->GetImpactParameters(xv, yv);
  fDcaDaughter = xv[0];
  fDcazDaughter = xv[1];
  fSigmaYXDaughter = yv[0];
  fSigmaXYZDaughter = yv[1];
  fSigmaZDaughter = yv[2];
  fDcaSecDaughter = TMath::Abs(exTrack1->GetD(TertVertex[0], TertVertex[1], kMagF));
  fNclsDaughter = track1->GetTPCNcls();
  fNclsITSDaughter = track1->GetITSNcls();
  fChi2Daughter = track1->GetTPCchi2() / (Float_t)track1->GetTPCclusters(0);
  fEtaDaughter = track1->Eta();
  fPhiDaughter = track1->Phi();
  fGeoLengthDaughter = AliAnalysisTaskDoubleHypNucTree::GeoLength(*track1);
  fTOFSignalDaughter = AliAnalysisTaskDoubleHypNucTree::GetTOFSignal(*track1);
  fPtUncertDaughter = TMath::Sqrt(exTrack1->GetSigma1Pt2()) * fptDaughter;
  fTPCRefitDaughter = (track1->GetStatus() & AliESDtrack::kTPCrefit) != 0;
  fITSRefitDaughter = (track1->GetStatus() & AliESDtrack::kITSrefit) != 0;
  fdEdxSigmaTriton = AliAnalysisTaskDoubleHypNucTree::Bethe(*track1, AliPID::ParticleMass(AliPID::kTriton), 1, fBetheParamsT);
  fdEdxSigmaAlpha = AliAnalysisTaskDoubleHypNucTree::Bethe(*track1, AliPID::ParticleMass(AliPID::kAlpha), 2, fBetheParamsHe);
  fITSLayer1Daughter = track1->HasPointOnITSLayer(0);
  fITSLayer2Daughter = track1->HasPointOnITSLayer(1);
  fITSLayer3Daughter = track1->HasPointOnITSLayer(2);
  fITSLayer4Daughter = track1->HasPointOnITSLayer(3);
  fITSLayer5Daughter = track1->HasPointOnITSLayer(4);
  fITSLayer6Daughter = track1->HasPointOnITSLayer(5);
  // _________________________________ //
  // __ p track __ //
  if(kDecayChannel == 1){
    fpDaughter1 = track2->GetInnerParam()->GetP();
    if(fMCtrue) {
      part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(track2->GetLabel()))->Particle());
      fLabelDaughter1KF = TMath::Abs(part->PdgCode()) == fgkPdgCode[kPDGProton];
    }
    fpxDaughter1KF = daughter2.GetPx();
    fpyDaughter1KF = daughter2.GetPy();
    fpzDaughter1KF = daughter2.GetPz();
    fptDaughter1KF = daughter2.GetPt();
    fEDaughter1KF = daughter2.GetE();
    fyDaughter1KF = daughter2.GetRapidity();
    fpxDaughter1 = particle2->Px();
    fpyDaughter1 = particle2->Py();
    fpzDaughter1 = particle2->Pz();
    fptDaughter1 = particle2->Pt();
    fEDaughter1 = particle2->E();
    fyDaughter1 = particle2->Rapidity();
    fdEdxDaughter1 = track2->GetTPCsignal();
    //fdEdxSigmaDaughter1 = AliAnalysisTaskDoubleHypNucTree::Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
    fdEdxSigmaDaughter1 = fPID->NumberOfSigmasTPC(track2, AliPID::kProton);
    //if (fMCtrue) fdEdxSigmaDaughter1 = fPID->NumberOfSigmasTPC(track2, AliPID::kProton);
    track2->GetImpactParameters(xv, yv);
    fDcaDaughter1 = xv[0];
    fDcazDaughter1 = xv[1];
    fSigmaYXDaughter1 = yv[0];
    fSigmaXYZDaughter1 = yv[1];
    fSigmaZDaughter1 = yv[2];
    if (kDecayChannel == 1) fDcaSecDaughter1 = TMath::Abs(exTrack2->GetD(TertVertex[0], TertVertex[1], kMagF));
    else                   fDcaSecDaughter1 = TMath::Abs(exTrack2->GetD(SecVertex[0], SecVertex[1], kMagF));
    fNclsDaughter1 = track2->GetTPCNcls();
    fNclsITSDaughter1 = track2->GetITSNcls();
    fChi2Daughter1 = track2->GetTPCchi2() / (Float_t)track2->GetTPCclusters(0);
    fEtaDaughter1 = track2->Eta();
    fPhiDaughter1 = track2->Phi();
    fGeoLengthDaughter1 = AliAnalysisTaskDoubleHypNucTree::GeoLength(*track2);
    fTOFSignalDaughter1 = AliAnalysisTaskDoubleHypNucTree::GetTOFSignal(*track2);
    fPtUncertDaughter1 = TMath::Sqrt(exTrack2->GetSigma1Pt2()) * fptDaughter1;
    fTPCRefitDaughter1 = (track2->GetStatus() & AliESDtrack::kTPCrefit) != 0;
    fITSRefitDaughter1 = (track2->GetStatus() & AliESDtrack::kITSrefit) != 0;
    fdEdxSigmaPion = fPID->NumberOfSigmasTPC(track2, AliPID::kPion);
    fdEdxSigmaDeuteron = AliAnalysisTaskDoubleHypNucTree::Bethe(*track2, AliPID::ParticleMass(AliPID::kDeuteron), 1, fBetheParamsT);
    fITSLayer1Daughter1 = track2->HasPointOnITSLayer(0);
    fITSLayer2Daughter1 = track2->HasPointOnITSLayer(1);
    fITSLayer3Daughter1 = track2->HasPointOnITSLayer(2);
    fITSLayer4Daughter1 = track2->HasPointOnITSLayer(3);
    fITSLayer5Daughter1 = track2->HasPointOnITSLayer(4);
    fITSLayer6Daughter1 = track2->HasPointOnITSLayer(5);
  }
  // _________________________________ //
  // __ pion track __ //
  fpDaughter2 = track3->GetInnerParam()->GetP();
  if(fMCtrue) {
    part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(track3->GetLabel()))->Particle());
    fLabelDaughter2KF = TMath::Abs(part->PdgCode()) == fgkPdgCode[kPDGPionPlus];
  }
  fpxDaughter2KF = daughter3.GetPx();
  fpyDaughter2KF = daughter3.GetPy();
  fpzDaughter2KF = daughter3.GetPz();
  fptDaughter2KF = daughter3.GetPt();
  fEDaughter2KF = daughter3.GetE();
  fyDaughter2KF = daughter3.GetRapidity();
  fpxDaughter2 = particle3->Px();
  fpyDaughter2 = particle3->Py();
  fpzDaughter2 = particle3->Pz();
  fptDaughter2 = particle3->Pt();
  fEDaughter2 = particle3->E();
  fyDaughter2 = particle3->Rapidity();
  fdEdxDaughter2 = track3->GetTPCsignal();
  fdEdxSigmaDaughter2 = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
  track3->GetImpactParameters(xv, yv);
  fDcaDaughter2 = xv[0];
  fDcazDaughter2 = xv[1];
  fSigmaYXDaughter2 = yv[0];
  fSigmaXYZDaughter2 = yv[1];
  fSigmaZDaughter2 = yv[2];
  if (kDecayChannel == 1) fDcaSecDaughter2 = TMath::Abs(exTrack3->GetD(TertVertex[0], TertVertex[1], kMagF));
  else                   fDcaSecDaughter2 = TMath::Abs(exTrack3->GetD(SecVertex[0], SecVertex[1], kMagF));
  fNclsDaughter2 = track3->GetTPCNcls();
  fNclsITSDaughter2 = track3->GetITSNcls();
  fChi2Daughter2 = track3->GetTPCchi2() / (Float_t)track3->GetTPCclusters(0);
  fEtaDaughter2 = track3->Eta();
  fPhiDaughter2 = track3->Phi();
  fGeoLengthDaughter2 = AliAnalysisTaskDoubleHypNucTree::GeoLength(*track3);
  fTOFSignalDaughter2 = AliAnalysisTaskDoubleHypNucTree::GetTOFSignal(*track3);
  fPtUncertDaughter2 = TMath::Sqrt(exTrack3->GetSigma1Pt2()) * fptDaughter2;
  fTPCRefitDaughter2 = (track3->GetStatus() & AliESDtrack::kTPCrefit) != 0;
  fITSRefitDaughter2 = (track3->GetStatus() & AliESDtrack::kITSrefit) != 0;
  fITSLayer1Daughter2 = track3->HasPointOnITSLayer(0);
  fITSLayer2Daughter2 = track3->HasPointOnITSLayer(1);
  fITSLayer3Daughter2 = track3->HasPointOnITSLayer(2);
  fITSLayer4Daughter2 = track3->HasPointOnITSLayer(3);
  fITSLayer5Daughter2 = track3->HasPointOnITSLayer(4);
  fITSLayer6Daughter2 = track3->HasPointOnITSLayer(5);
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::SetDaughterInformationKF(KFParticle daughter1, KFParticle daughter2, KFParticle daughter3, KFParticle daughter4, Int_t kDecayChannel) {

  Float_t xv[2], yv[3];
  AliMCParticle *part;  
  // _________________________________ //
  // __ He3 track __ //
  fpDaughter = track1->GetInnerParam()->GetP();
  if(fMCtrue) {
    part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(track1->GetLabel()))->Particle());
    fLabelDaughterKF = TMath::Abs(part->PdgCode()) == fgkPdgCode[kPDGHelium3];
  }  
  fpxDaughterKF = daughter1.GetPx();
  fpyDaughterKF = daughter1.GetPy();
  fpzDaughterKF = daughter1.GetPz();
  fptDaughterKF = daughter1.GetPt();
  fEDaughterKF = daughter1.GetE();
  fyDaughterKF = daughter1.GetRapidity();
  fpxDaughter = particle1->Px();
  fpyDaughter = particle1->Py();
  fpzDaughter = particle1->Pz();
  fptDaughter = particle1->Pt();
  fEDaughter = particle1->E();
  fyDaughter = particle1->Rapidity();
  fdEdxDaughter = track1->GetTPCsignal();
  fdEdxSigmaDaughter = AliAnalysisTaskDoubleHypNucTree::Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
  if (fMCtrue) fdEdxSigmaDaughter = fPID->NumberOfSigmasTPC(track1, AliPID::kHe3);
  track1->GetImpactParameters(xv, yv);
  fDcaDaughter = xv[0];
  fDcazDaughter = xv[1];
  fSigmaYXDaughter = yv[0];
  fSigmaXYZDaughter = yv[1];
  fSigmaZDaughter = yv[2];
  fDcaSecDaughter = TMath::Abs(exTrack1->GetD(TertVertex[0], TertVertex[1], kMagF));
  fNclsDaughter = track1->GetTPCNcls();
  fNclsITSDaughter = track1->GetITSNcls();
  fChi2Daughter = track1->GetTPCchi2() / (Float_t)track1->GetTPCclusters(0);
  fEtaDaughter = track1->Eta();
  fPhiDaughter = track1->Phi();
  fGeoLengthDaughter = AliAnalysisTaskDoubleHypNucTree::GeoLength(*track1);
  fTOFSignalDaughter = AliAnalysisTaskDoubleHypNucTree::GetTOFSignal(*track1);
  fPtUncertDaughter = TMath::Sqrt(exTrack1->GetSigma1Pt2()) * fptDaughter;
  fTPCRefitDaughter = (track1->GetStatus() & AliESDtrack::kTPCrefit) != 0;
  fITSRefitDaughter = (track1->GetStatus() & AliESDtrack::kITSrefit) != 0;
  fdEdxSigmaTriton = AliAnalysisTaskDoubleHypNucTree::Bethe(*track1, AliPID::ParticleMass(AliPID::kTriton), 1, fBetheParamsT);
  fdEdxSigmaAlpha = AliAnalysisTaskDoubleHypNucTree::Bethe(*track1, AliPID::ParticleMass(AliPID::kAlpha), 2, fBetheParamsHe);
  fITSLayer1Daughter = track1->HasPointOnITSLayer(0);
  fITSLayer2Daughter = track1->HasPointOnITSLayer(1);
  fITSLayer3Daughter = track1->HasPointOnITSLayer(2);
  fITSLayer4Daughter = track1->HasPointOnITSLayer(3);
  fITSLayer5Daughter = track1->HasPointOnITSLayer(4);
  fITSLayer6Daughter = track1->HasPointOnITSLayer(5);
  // _________________________________ //
  // __ p track __ //
  fpDaughter1 = track2->GetInnerParam()->GetP();
  if(fMCtrue) {
    part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(track2->GetLabel()))->Particle());
    fLabelDaughter1KF = TMath::Abs(part->PdgCode()) == fgkPdgCode[kPDGProton];
  }
  fpxDaughter1KF = daughter2.GetPx();
  fpyDaughter1KF = daughter2.GetPy();
  fpzDaughter1KF = daughter2.GetPz();
  fptDaughter1KF = daughter2.GetPt();
  fEDaughter1KF = daughter2.GetE();
  fyDaughter1KF = daughter2.GetRapidity();
  fpxDaughter1 = particle2->Px();
  fpyDaughter1 = particle2->Py();
  fpzDaughter1 = particle2->Pz();
  fptDaughter1 = particle2->Pt();
  fEDaughter1 = particle2->E();
  fyDaughter1 = particle2->Rapidity();
  fdEdxDaughter1 = track2->GetTPCsignal();
  //fdEdxSigmaDaughter1 = AliAnalysisTaskDoubleHypNucTree::Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
  fdEdxSigmaDaughter1 = fPID->NumberOfSigmasTPC(track2, AliPID::kProton);
  //if (fMCtrue) fdEdxSigmaDaughter1 = fPID->NumberOfSigmasTPC(track2, AliPID::kProton);
  track2->GetImpactParameters(xv, yv);
  fDcaDaughter1 = xv[0];
  fDcazDaughter1 = xv[1];
  fSigmaYXDaughter1 = yv[0];
  fSigmaXYZDaughter1 = yv[1];
  fSigmaZDaughter1 = yv[2];
  if (kDecayChannel == 1) fDcaSecDaughter1 = TMath::Abs(exTrack2->GetD(TertVertex[0], TertVertex[1], kMagF));
  else                   fDcaSecDaughter1 = TMath::Abs(exTrack2->GetD(SecVertex[0], SecVertex[1], kMagF));
  fNclsDaughter1 = track2->GetTPCNcls();
  fNclsITSDaughter1 = track2->GetITSNcls();
  fChi2Daughter1 = track2->GetTPCchi2() / (Float_t)track2->GetTPCclusters(0);
  fEtaDaughter1 = track2->Eta();
  fPhiDaughter1 = track2->Phi();
  fGeoLengthDaughter1 = AliAnalysisTaskDoubleHypNucTree::GeoLength(*track2);
  fTOFSignalDaughter1 = AliAnalysisTaskDoubleHypNucTree::GetTOFSignal(*track2);
  fPtUncertDaughter1 = TMath::Sqrt(exTrack2->GetSigma1Pt2()) * fptDaughter1;
  fTPCRefitDaughter1 = (track2->GetStatus() & AliESDtrack::kTPCrefit) != 0;
  fITSRefitDaughter1 = (track2->GetStatus() & AliESDtrack::kITSrefit) != 0;
  fdEdxSigmaPion = fPID->NumberOfSigmasTPC(track2, AliPID::kPion);
  fdEdxSigmaDeuteron = AliAnalysisTaskDoubleHypNucTree::Bethe(*track2, AliPID::ParticleMass(AliPID::kDeuteron), 1, fBetheParamsT);
  fITSLayer1Daughter1 = track2->HasPointOnITSLayer(0);
  fITSLayer2Daughter1 = track2->HasPointOnITSLayer(1);
  fITSLayer3Daughter1 = track2->HasPointOnITSLayer(2);
  fITSLayer4Daughter1 = track2->HasPointOnITSLayer(3);
  fITSLayer5Daughter1 = track2->HasPointOnITSLayer(4);
  fITSLayer6Daughter1 = track2->HasPointOnITSLayer(5);
  // _________________________________ //
  // __ pion track __ //
  fpDaughter2 = track3->GetInnerParam()->GetP();
  if(fMCtrue) {
    part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(track3->GetLabel()))->Particle());
    fLabelDaughter2KF = TMath::Abs(part->PdgCode()) == fgkPdgCode[kPDGPionPlus];
  }
  fpxDaughter2KF = daughter3.GetPx();
  fpyDaughter2KF = daughter3.GetPy();
  fpzDaughter2KF = daughter3.GetPz();
  fptDaughter2KF = daughter3.GetPt();
  fEDaughter2KF = daughter3.GetE();
  fyDaughter2KF = daughter3.GetRapidity();  
  fpxDaughter2 = particle3->Px();
  fpyDaughter2 = particle3->Py();
  fpzDaughter2 = particle3->Pz();
  fptDaughter2 = particle3->Pt();
  fEDaughter2 = particle3->E();
  fyDaughter2 = particle3->Rapidity();
  fdEdxDaughter2 = track3->GetTPCsignal();
  fdEdxSigmaDaughter2 = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
  track3->GetImpactParameters(xv, yv);
  fDcaDaughter2 = xv[0];
  fDcazDaughter2 = xv[1];
  fSigmaYXDaughter2 = yv[0];
  fSigmaXYZDaughter2 = yv[1];
  fSigmaZDaughter2 = yv[2];
  if (kDecayChannel == 1) fDcaSecDaughter2 = TMath::Abs(exTrack3->GetD(TertVertex[0], TertVertex[1], kMagF));
  else                   fDcaSecDaughter2 = TMath::Abs(exTrack3->GetD(SecVertex[0], SecVertex[1], kMagF));
  fNclsDaughter2 = track3->GetTPCNcls();
  fNclsITSDaughter2 = track3->GetITSNcls();
  fChi2Daughter2 = track3->GetTPCchi2() / (Float_t)track3->GetTPCclusters(0);
  fEtaDaughter2 = track3->Eta();
  fPhiDaughter2 = track3->Phi();
  fGeoLengthDaughter2 = AliAnalysisTaskDoubleHypNucTree::GeoLength(*track3);
  fTOFSignalDaughter2 = AliAnalysisTaskDoubleHypNucTree::GetTOFSignal(*track3);
  fPtUncertDaughter2 = TMath::Sqrt(exTrack3->GetSigma1Pt2()) * fptDaughter2;
  fTPCRefitDaughter2 = (track3->GetStatus() & AliESDtrack::kTPCrefit) != 0;
  fITSRefitDaughter2 = (track3->GetStatus() & AliESDtrack::kITSrefit) != 0;
  fITSLayer1Daughter2 = track3->HasPointOnITSLayer(0);
  fITSLayer2Daughter2 = track3->HasPointOnITSLayer(1);
  fITSLayer3Daughter2 = track3->HasPointOnITSLayer(2);
  fITSLayer4Daughter2 = track3->HasPointOnITSLayer(3);
  fITSLayer5Daughter2 = track3->HasPointOnITSLayer(4);
  fITSLayer6Daughter2 = track3->HasPointOnITSLayer(5);
  // _________________________________ //
  // __ pion track __ //
  fpDaughter3 = track4->GetInnerParam()->GetP();
  if(fMCtrue) {
    part = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(track4->GetLabel()))->Particle());
    fLabelDaughter3KF = TMath::Abs(part->PdgCode()) == fgkPdgCode[kPDGPionPlus];
  }
  fpxDaughter3KF = daughter4.GetPx();
  fpyDaughter3KF = daughter4.GetPy();
  fpzDaughter3KF = daughter4.GetPz();
  fptDaughter3KF = daughter4.GetPt();
  fEDaughter3KF = daughter4.GetE();
  fyDaughter3KF = daughter4.GetRapidity();
  fpxDaughter3 = particle4->Px();
  fpyDaughter3 = particle4->Py();
  fpzDaughter3 = particle4->Pz();
  fptDaughter3 = particle4->Pt();
  fEDaughter3 = particle4->E();
  fyDaughter3 = particle4->Rapidity();
  fdEdxDaughter3 = track4->GetTPCsignal();
  fdEdxSigmaDaughter3 = fPID->NumberOfSigmasTPC(track4, AliPID::kPion);
  track4->GetImpactParameters(xv, yv);
  fDcaDaughter3 = xv[0];
  fDcazDaughter3 = xv[1];
  fSigmaYXDaughter3 = yv[0];
  fSigmaXYZDaughter3 = yv[1];
  fSigmaZDaughter3 = yv[2];
  if (kDecayChannel == 1) fDcaSecDaughter3 = TMath::Abs(exTrack4->GetD(SecVertex[0], SecVertex[1], kMagF));
  else                   fDcaSecDaughter3 = TMath::Abs(exTrack4->GetD(TertVertex[0], TertVertex[1], kMagF));
  fNclsDaughter3 = track4->GetTPCNcls();
  fNclsITSDaughter3 = track4->GetITSNcls();
  fChi2Daughter3 = track4->GetTPCchi2() / (Float_t)track4->GetTPCclusters(0);
  fEtaDaughter3 = track4->Eta();
  fPhiDaughter3 = track4->Phi();
  fGeoLengthDaughter3 = AliAnalysisTaskDoubleHypNucTree::GeoLength(*track4);
  fTOFSignalDaughter3 = AliAnalysisTaskDoubleHypNucTree::GetTOFSignal(*track4);
  fPtUncertDaughter3 = TMath::Sqrt(exTrack4->GetSigma1Pt2()) * fptDaughter3;
  fTPCRefitDaughter3 = (track4->GetStatus() & AliESDtrack::kTPCrefit) != 0;
  fITSRefitDaughter3 = (track4->GetStatus() & AliESDtrack::kITSrefit) != 0;
  fITSLayer1Daughter3 = track4->HasPointOnITSLayer(0);
  fITSLayer2Daughter3 = track4->HasPointOnITSLayer(1);
  fITSLayer3Daughter3 = track4->HasPointOnITSLayer(2);
  fITSLayer4Daughter3 = track4->HasPointOnITSLayer(3);
  fITSLayer5Daughter3 = track4->HasPointOnITSLayer(4);
  fITSLayer6Daughter3 = track4->HasPointOnITSLayer(5);
}
// _________________________________________________ //
KFParticle AliAnalysisTaskDoubleHypNucTree::CreateKFParticle(AliExternalTrackParam& track, Double_t Mass, Int_t Charge) {

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
KFVertex AliAnalysisTaskDoubleHypNucTree::CreateKFVertex(const AliVVertex& vertex) {

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

  Int_t nTrdTracks = fESDevent->GetNumberOfTrdTracks();
  // __ Data: get TRD trigger information from trigger classes __ //
  TString classes = fESDevent->GetFiredTriggerClasses();
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
  Bool_t isTriggered = kFALSE;
  if (MB || HMV0 || HMSPD || HNU || HQU || Central || SemiCentral) isTriggered = kTRUE;
  return isTriggered;
}
// _________________________________________________ //
// __ calculate stack label of specific daughter __ //
Int_t AliAnalysisTaskDoubleHypNucTree::GetLabel(Int_t labelGrandMother, Int_t particlePdgCode) {

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
Int_t AliAnalysisTaskDoubleHypNucTree::GetLabel(Int_t labelGrandMother, Int_t particlePdgCode, Int_t motherparticlePdgCode) {
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
void AliAnalysisTaskDoubleHypNucTree::MCGenerated() {

  // Monte Carlo for genenerated particles                                                                                                                                          
  stackN = 0;
  for (stackN = 0; stackN < mcEvent->GetNumberOfTracks(); stackN++) //loop over stack                                                                                                      
    {

      ParticleMother = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(stackN))->Particle());
      PDGCodeMother = ParticleMother->PdgCode();

      //DoubleHyperHydrogen4
      if (PDGCodeMother == fgkPdgCode[kPDGDoubleHyperHydrogen4]) //check mother PDG                                                                                                        
	{
	    AliAnalysisTaskDoubleHypNucTree::MCFourBodyDecay(stackN, ParticleMother, PDGCodeMother, 1, fgkPdgCode[kPDGHelium3],
								fgkPdgCode[kPDGProton], fgkPdgCode[kPDGPionMinus], fgkPdgCode[kPDGPionMinus],								
								AliPID::ParticleMass(AliPID::kHe3), AliPID::ParticleMass(AliPID::kProton), AliPID::ParticleMass(AliPID::kPion), AliPID::ParticleMass(AliPID::kPion));
	}
      //AntiDoubleHyperHydrogen4
      if (PDGCodeMother == fgkPdgCode[kPDGAntiDoubleHyperHydrogen4]) //check mother PDG
	{
	    AliAnalysisTaskDoubleHypNucTree::MCFourBodyDecay(stackN, ParticleMother, PDGCodeMother, 1, fgkPdgCode[kPDGAntiHelium3],
								fgkPdgCode[kPDGAntiProton], fgkPdgCode[kPDGPionPlus], fgkPdgCode[kPDGPionPlus],
								AliPID::ParticleMass(AliPID::kHe3), AliPID::ParticleMass(AliPID::kProton), AliPID::ParticleMass(AliPID::kPion), AliPID::ParticleMass(AliPID::kPion));
	}
      //HyperHelium4
      if (PDGCodeMother == fgkPdgCode[kPDGHyperHelium4]) //check mother PDG                                                                                                        
	{
	  AliAnalysisTaskDoubleHypNucTree::MCThreeBodyDecay(stackN, ParticleMother, PDGCodeMother, fgkPdgCode[kPDGHelium3],
							       fgkPdgCode[kPDGProton], fgkPdgCode[kPDGPionMinus],
							       AliPID::ParticleMass(AliPID::kHe3), AliPID::ParticleMass(AliPID::kProton), AliPID::ParticleMass(AliPID::kPion));
	}
      //AntiHyperHelium4
      if (PDGCodeMother == fgkPdgCode[kPDGAntiHyperHelium4]) //check mother PDG
	{

	  AliAnalysisTaskDoubleHypNucTree::MCThreeBodyDecay(stackN, ParticleMother, PDGCodeMother, fgkPdgCode[kPDGAntiHelium3],
							       fgkPdgCode[kPDGAntiProton], fgkPdgCode[kPDGPionPlus],
							       AliPID::ParticleMass(AliPID::kHe3), AliPID::ParticleMass(AliPID::kProton), AliPID::ParticleMass(AliPID::kPion));
	}      
      //HyperHydrogen3
      if (PDGCodeMother == fgkPdgCode[kPDGHyperHydrogen3]) //check mother PDG                                                                                                        
	{
	  
	  AliAnalysisTaskDoubleHypNucTree::MCTwoBodyDecay(stackN, ParticleMother, PDGCodeMother, fgkPdgCode[kPDGHelium3],
							     fgkPdgCode[kPDGPionMinus], AliPID::ParticleMass(AliPID::kHe3), AliPID::ParticleMass(AliPID::kPion));
	}
      //AntiHyperHydrogen3
      if (PDGCodeMother == fgkPdgCode[kPDGAntiHyperHydrogen3]) //check mother PDG
	{
	  
	  AliAnalysisTaskDoubleHypNucTree::MCTwoBodyDecay(stackN, ParticleMother, PDGCodeMother, fgkPdgCode[kPDGAntiHelium3],
							     fgkPdgCode[kPDGPionPlus], AliPID::ParticleMass(AliPID::kHe3), AliPID::ParticleMass(AliPID::kPion));
	}
      if (ParticleMother) delete ParticleMother;
    }//end loop over stack                                                                                                                                                                              
}
//_________________________________________________ 
void AliAnalysisTaskDoubleHypNucTree::MCFourBodyDecay(Int_t stackN, AliMCParticle* ParticleMother, Long_t PDGMother, Int_t kDecayChannel,
							 Long_t PDGFirstDaughter, Long_t PDGSecondDaughter, Long_t PDGThirdDaughter, Long_t PDGFourthDaughter,
							 Double_t massFirstDaughter, Double_t massSecondDaughter, Double_t massThirdDaughter, Double_t massFourthDaughter) {

  //Check decay channel
  Int_t flabel = mcEvent->GetLabelOfParticleFirstDaughter(stackN);
  Int_t llabel = mcEvent->GetLabelOfParticleLastDaughter(stackN);
  sign = PDGMother / TMath::Abs(PDGMother);
  kDecayChannel = 0; 
  for(int i = flabel; i < llabel +1; i++){
    AliMCParticle *cparticle = new AliMCParticle(mcEvent->GetTrack(i)->Particle());
    if(cparticle->PdgCode() == sign * fgkPdgCode[kPDGHyperHelium4]) kDecayChannel = 1;
    if(cparticle->PdgCode() == sign * fgkPdgCode[kPDGHyperHydrogen3]) kDecayChannel = 2;
  }
  if(!kDecayChannel) return;

  if (kDecayChannel == 1) {
    label1 = TMath::Abs(AliAnalysisTaskDoubleHypNucTree::GetLabel(stackN, PDGFirstDaughter, sign * fgkPdgCode[kPDGHyperHelium4]));
    label2 = TMath::Abs(AliAnalysisTaskDoubleHypNucTree::GetLabel(stackN, PDGSecondDaughter, sign * fgkPdgCode[kPDGHyperHelium4]));
    label3 = TMath::Abs(AliAnalysisTaskDoubleHypNucTree::GetLabel(stackN, PDGThirdDaughter, sign * fgkPdgCode[kPDGHyperHelium4]));
    label4 = TMath::Abs(AliAnalysisTaskDoubleHypNucTree::GetLabel(stackN, PDGFourthDaughter));
    
  }
  else {
    label1 = TMath::Abs(AliAnalysisTaskDoubleHypNucTree::GetLabel(stackN, PDGFirstDaughter, sign * fgkPdgCode[kPDGHyperHydrogen3]));
    label2 = TMath::Abs(AliAnalysisTaskDoubleHypNucTree::GetLabel(stackN, PDGSecondDaughter));
    label3 = TMath::Abs(AliAnalysisTaskDoubleHypNucTree::GetLabel(stackN, PDGThirdDaughter));
    label4 = TMath::Abs(AliAnalysisTaskDoubleHypNucTree::GetLabel(stackN, PDGFourthDaughter, sign * fgkPdgCode[kPDGHyperHydrogen3]));
  }
  FirstDaughter  = new AliMCParticle(mcEvent->GetTrack(label1)->Particle());
  SecondDaughter = new AliMCParticle(mcEvent->GetTrack(label2)->Particle());
  ThirdDaughter  = new AliMCParticle(mcEvent->GetTrack(label3)->Particle());
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
	  fPrimary4LHe = 0;//mcEvent->IsSecondaryFromWeakDecay(AliAnalysisTaskDoubleHypNucTree::GetLabel(stackN, sign * fgkPdgCode[kPDGHyperHelium4])) == 0;	  

	  fptDaughter = particle1->Pt();
	  fyDaughter = particle1->Rapidity();
	  fptDaughter1 = particle2->Pt();
	  fyDaughter1 = particle2->Rapidity();
	  fptDaughter2 = particle3->Pt();
	  fyDaughter2 = particle3->Rapidity();
	  fptDaughter3 = particle4->Pt();
	  fyDaughter3 = particle4->Rapidity();
	  fTreeGen->Fill();
	  AliAnalysisTaskDoubleHypNucTree::ResetVals("");
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
void AliAnalysisTaskDoubleHypNucTree::MCThreeBodyDecay(Int_t stackN, AliMCParticle* ParticleMother, Long_t PDGMother,
							  Long_t PDGFirstDaughter, Long_t PDGSecondDaughter, Long_t PDGThirdDaughter,
							  Double_t massFirstDaughter, Double_t massSecondDaughter, Double_t massThirdDaughter) {

  sign = PDGMother / TMath::Abs(PDGMother);

  FirstDaughter = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(AliAnalysisTaskDoubleHypNucTree::GetLabel(stackN, PDGFirstDaughter)))->Particle());
  SecondDaughter = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(AliAnalysisTaskDoubleHypNucTree::GetLabel(stackN, PDGSecondDaughter)))->Particle());
  ThirdDaughter = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(AliAnalysisTaskDoubleHypNucTree::GetLabel(stackN, PDGThirdDaughter)))->Particle());

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
	AliAnalysisTaskDoubleHypNucTree::ResetVals("");
      }
    }
  }
  if (FirstDaughter)  delete FirstDaughter;
  if (SecondDaughter) delete SecondDaughter;
  if (ThirdDaughter)  delete ThirdDaughter;
}
// _________________________________________________ //
void AliAnalysisTaskDoubleHypNucTree::MCTwoBodyDecay(Int_t stackN, AliMCParticle* ParticleMother, Long_t PDGMother,
							Long_t PDGFirstDaughter, Long_t PDGSecondDaughter,
							Double_t massFirstDaughter, Double_t massSecondDaughter) {

  sign = PDGMother / TMath::Abs(PDGMother);

  FirstDaughter = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(AliAnalysisTaskDoubleHypNucTree::GetLabel(stackN, PDGFirstDaughter)))->Particle());
  SecondDaughter = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(AliAnalysisTaskDoubleHypNucTree::GetLabel(stackN, PDGSecondDaughter)))->Particle());

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
      fDecayChannel = 2;

      fmSubMother = sublorentzsum->M();
      fptSubMother = sublorentzsum->Pt();
      fySubMother = sublorentzsum->Rapidity();
      fctSubMother = (sublorentzsum->M() * TMath::Sqrt(TMath::Power(dd[0], 2) + TMath::Power(dd[1], 2) + TMath::Power(dd[2], 2))) / sublorentzsum->P();

      fptDaughter = particle1->Pt();
      fyDaughter = particle1->Rapidity();
      fptDaughter1 = particle2->Pt();
      fyDaughter1 = particle2->Rapidity();
      gTreeGen->Fill();
      AliAnalysisTaskDoubleHypNucTree::ResetVals("");
    }
  }
  if (FirstDaughter)  delete FirstDaughter;
  if (SecondDaughter) delete SecondDaughter;
}
// _________________________________________________ //
Double_t AliAnalysisTaskDoubleHypNucTree::Bethe(const AliESDtrack& track, Double_t mass, Int_t charge, Double_t* params) {
  Double_t expected = charge * charge * AliExternalTrackParam::BetheBlochAleph(charge*track.GetInnerParam()->GetP()/mass, params[0], params[1], params[2], params[3], params[4]);
  Double_t sigma = expected * params[5];
  if (TMath::IsNaN(expected)) return -999;
  return (track.GetTPCsignal() - expected) / sigma;
}
// _________________________________________________ //
Double_t AliAnalysisTaskDoubleHypNucTree::GetTOFSignal(const AliESDtrack& track) {
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
    if(beta*beta >= 1) return -99;
    gamma = 1 / TMath::Sqrt(1 - beta * beta);
    if(gamma*gamma <= 1) return -99;
    mass = (track.GetInnerParam()->GetP()) / TMath::Sqrt(gamma * gamma - 1); // using inner TPC mom. as approx.
  }
  return mass;
}
// _________________________________________________ //
Double_t AliAnalysisTaskDoubleHypNucTree::GeoLength(const AliESDtrack& track) {
  Double_t deadZoneWidth = 3.0;
  //AliExternalTrackParam *p = new AliExternalTrackParam();
  //p->CopyFromVTrack(&track);
  Double_t lengthInActiveZone;	
  lengthInActiveZone = GetLengthInActiveZone(track.GetInnerParam(), deadZoneWidth, 220, track.GetESDEvent()->GetMagneticField(), 0); 
  return lengthInActiveZone;
}
// _________________________________________________ //
Double_t AliAnalysisTaskDoubleHypNucTree::GetLengthInActiveZone(const AliExternalTrackParam  *paramT, Double_t deltaY, Double_t deltaZ, Double_t bz, Double_t exbPhi) {
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
  const Double_t rIn=85;
  const Double_t rOut=245;
  Double_t xyz[3], pxyz[3];
  if (paramT->GetXYZAt(rIn,bz,xyz)){
    paramT->GetPxPyPzAt(rIn,bz,pxyz);
  }else{
    paramT->GetXYZ(xyz);
    paramT->GetPxPyPz(pxyz);
  }
  //
  Double_t dca   = -paramT->GetD(0,0,bz);  // get impact parameter distance to point (0,0)
  Double_t radius= TMath::Abs(1/paramT->GetC(bz));  //
  Double_t sign  = paramT->GetSign()*TMath::Sign(1.,bz)*(-1.);
  Double_t R0    = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);   // radius at current point
  Double_t phiR0 = TMath::ATan2(xyz[1],xyz[0]);                // angle of given point
  Double_t dPhiR0= -TMath::ASin((dca*dca-2*dca*radius*sign+R0*R0)/(2*R0*(dca-radius*sign)));
  Double_t phi0  = phiR0-(dPhiR0);  // global phi offset to be added
  //
  //
  AliExternalTrackParam paramR=(*paramT);
  Double_t length=0;
  for (Double_t R=rIn; R<=rOut; R++){
    Double_t sinPhi=(dca*dca-2*dca*radius*sign+R*R)/(2*R*(dca-radius*sign));
    if (TMath::Abs(sinPhi)>=1) continue;
    Double_t dphi     = -TMath::ASin(sinPhi);
    Double_t phi      = phi0+dphi;                           // global phi
    Int_t    sector   = TMath::Nint(9*phi/(TMath::Pi()));
    Double_t dPhiEdge = phi-(sector*TMath::Pi()/9)+exbPhi;   // distance to sector boundary in rphi
    Double_t dX   = R*TMath::Cos(phi)-xyz[0];
    Double_t dY   = R*TMath::Sin(phi)-xyz[1];
    Double_t deltaPhi = 2*TMath::ASin(0.5*TMath::Sqrt(dX*dX+dY*dY)/radius);
    Double_t z = xyz[2]+deltaPhi*radius*paramT->GetTgl();
    if (TMath::Abs(dPhiEdge*R)>deltaY && TMath::Abs(z)<deltaZ){
      length++;
    }    
  }
  return length;
}
// _________________________________________________ //
Float_t AliAnalysisTaskDoubleHypNucTree::GetInvPtDevFromBC(Int_t b, Int_t c) {
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
void AliAnalysisTaskDoubleHypNucTree::ResetVals(TString mode) {
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
  }
  if (mode == "") {
    PrimVertex[0] = 0.;
    PrimVertex[1] = 0.;
    PrimVertex[2] = 0.;
    SecVertex[0] = 0.;
    SecVertex[1] = 0.;
    SecVertex[2] = 0.;
    TertVertex[0] = 0.;
    TertVertex[1] = 0.;
    TertVertex[2] = 0.;   
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
    fisOnlineV0_13 = -1;
    fisOnlineV0_14 = -1;
    fisOnlineV0_23 = -1;
    fisOnlineV0_24 = -1;
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
    fPrimary4LHe = -1;
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
    fPropDCADaughter = -1;
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
    fPropDCADaughter1 = -1;
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
    fPropDCADaughter2 = -1;
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
    fPropDCADaughter3 = -1;
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
    fPropDCADaughter4 = -1;
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
    fyMotherKF   = -99;
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
    fySubMotherKF   = -99;
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
    fITSLayer1Daughter = -1;
    fITSLayer2Daughter = -1;
    fITSLayer3Daughter = -1;
    fITSLayer4Daughter = -1;
    fITSLayer5Daughter = -1;
    fITSLayer6Daughter = -1;
    fITSLayer1Daughter1 = -1;
    fITSLayer2Daughter1 = -1;
    fITSLayer3Daughter1 = -1;
    fITSLayer4Daughter1 = -1;
    fITSLayer5Daughter1 = -1;
    fITSLayer6Daughter1 = -1;
    fITSLayer1Daughter2 = -1;
    fITSLayer2Daughter2 = -1;
    fITSLayer3Daughter2 = -1;
    fITSLayer4Daughter2 = -1;
    fITSLayer5Daughter2 = -1;
    fITSLayer6Daughter2 = -1;
    fITSLayer1Daughter3 = -1;
    fITSLayer2Daughter3 = -1;
    fITSLayer3Daughter3 = -1;
    fITSLayer4Daughter3 = -1;
    fITSLayer5Daughter3 = -1;
    fITSLayer6Daughter3 = -1;
  }
  return;
}
// _________________________________________________ //
// __ set Bethe-Bloch parameter __ //
void AliAnalysisTaskDoubleHypNucTree::SetBetheBlochParams(Int_t runNumber) {
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
