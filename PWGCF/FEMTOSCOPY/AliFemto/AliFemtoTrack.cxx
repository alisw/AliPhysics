///
/// \file AliFemtoTrack.cxx
///


#include "AliFemtoTrack.h"

#include "SystemOfUnits.h"   // has "tesla" in it


AliFemtoTrack::AliFemtoTrack():
  fCharge(0),
  fPidProbElectron(0.0f),
  fPidProbPion(0.0f),
  fPidProbKaon(0.0f),
  fPidProbProton(0.0f),
  fPidProbMuon(0.0f),
  fTrackId(0),
  fTofPionTime(-100000.0f),
  fTofKaonTime(-100000.0f),
  fTofProtonTime(-100000.0f),
  fTofDeuteronTime(-100000.0f),
  fTofTritonTime(-100000.0f),
  fTofHe3Time(-100000.0f),
  fTofAlphaTime(-100000.0f),
  fP(0,0,0),
  fPt(0.0f),
  fInnerMomentum(0.0f),
  fHelix(),
  fFlags(0l),
  fLabel(0),
  fImpactD(0.0f),
  fImpactDprim(-10000.0f),
  fImpactDweak(-10000.0f),
  fImpactDmat(-10000.0f),
  fImpactZ(0.0f),
  fCdd(0.0f),
  fCdz(0.0f),
  fCzz(0.0f),
  fITSchi2(0.0f),
  fITSncls(0),
  fTPCchi2(0.0f),
  fTPCncls(0),
  fTPCnclsF(0),
  fTPCsignal(0.0f),
  fTPCsignalN(0),
  fTPCsignalS(0.0f),
  fVTOF(0.0f),
  fTOFsignal(0.0f),
  fNSigmaTPCPi(0.0f),
  fNSigmaTPCK(0.0f),
  fNSigmaTPCP(0.0f),
  fNSigmaTPCE(0.0f),
  fNSigmaTOFPi(0.0f),
  fNSigmaTOFK(0.0f),
  fNSigmaTOFP(0.0f),
  fNSigmaTOFE(0.0f),
  fNSigmaTPCD(0.0f),
  fNSigmaTPCT(0.0f),
  fNSigmaTPCH(0.0f),
  fNSigmaTPCA(0.0f),
  fNSigmaTOFD(0.0f),
  fNSigmaTOFT(0.0f),
  fNSigmaTOFH(0.0f),
  fNSigmaTOFA(0.0f),
  fMassTOF(0.0f),
  fSigmaToVertex(0.0f),
  fClusters(159),
  fShared(159),
  fNominalTpcEntrancePoint(0,0,0),
  fNominalTpcExitPoint(0,0,0),
  fNominalTpcPointShifted(0,0,0),
  fXatDCA(0.0),
  fYatDCA(0.0),
  fZatDCA(0.0),
  fHiddenInfo(nullptr),
  fTrueMomentum(nullptr),
  fEmissionPoint(nullptr),
  fPDGPid(0),
  fMass(0.0),
  fGlobalEmissionPoint(nullptr),
  fCorrPi(0.0),
  fCorrK(0.0),
  fCorrP(0.0),
  fCorrPiMinus(0.0),
  fCorrKMinus(0.0),
  fCorrPMinus(0.0),
  fCorrAll(0.0),
  fCorrD(0.0),
  fCorrT(0.0),
  fCorrH(0.0),
  fCorrA(0.0),
  fCorrDMinus(0.0),
  fCorrTMinus(0.0),
  fCorrHMinus(0.0)
  , fCorrAMinus(0.0)
{
  // Default constructor

  /// fill arrays with a value
  std::fill_n(fKinkIndexes, 3, 0);
  std::fill_n(fVertex, 3, -9999);
  std::fill_n(fHasPointOnITS, 6, kFALSE);
  std::fill_n(fNominalTpcPoints, 9, AliFemtoThreeVector(-9999, -9999, -9999));
}


AliFemtoTrack::AliFemtoTrack(const AliFemtoTrack& t) :
  fCharge(t.fCharge),
  fPidProbElectron(t.fPidProbElectron),
  fPidProbPion(t.fPidProbPion),
  fPidProbKaon(t.fPidProbKaon),
  fPidProbProton(t.fPidProbProton),
  fPidProbMuon(t.fPidProbMuon),
  fTrackId(t.fTrackId),
  fTofPionTime(t.fTofPionTime),
  fTofKaonTime(t.fTofKaonTime),
  fTofProtonTime(t.fTofProtonTime),
  fTofDeuteronTime(t.fTofDeuteronTime),
  fTofTritonTime(t.fTofTritonTime),
  fTofHe3Time(t.fTofHe3Time),
  fTofAlphaTime(t.fTofAlphaTime),
  fP(t.fP),
  fPt(t.fPt),
  fInnerMomentum(t.fInnerMomentum),
  fMultiplicity(t.fMultiplicity),
  fZvtx(t.fZvtx),
  fHelix(t.fHelix),
  fFlags(t.fFlags),
  fLabel(t.fLabel),
  fImpactD(t.fImpactD),
  fImpactDprim(t.fImpactDprim),
  fImpactDweak(t.fImpactDweak),
  fImpactDmat(t.fImpactDmat),
  fImpactZ(t.fImpactZ),
  fCdd(t.fCdd),
  fCdz(t.fCdz),
  fCzz(t.fCzz),
  fITSchi2(t.fITSchi2),
  fITSncls(t.fITSncls),
  fTPCchi2(t.fTPCchi2),
  fTPCncls(t.fTPCncls),
  fTPCnclsF(t.fTPCnclsF),
  fTPCsignal(t.fTPCsignal),
  fTPCsignalN(t.fTPCsignalN),
  fTPCsignalS(t.fTPCsignalS),
  fVTOF(t.fVTOF),
  fTOFsignal(t.fTOFsignal),
  fNSigmaTPCPi(t.fNSigmaTPCPi),
  fNSigmaTPCK(t.fNSigmaTPCK),
  fNSigmaTPCP(t.fNSigmaTPCP),
  fNSigmaTPCE(t.fNSigmaTPCE),
  fNSigmaTOFPi(t.fNSigmaTOFPi),
  fNSigmaTOFK(t.fNSigmaTOFK),
  fNSigmaTOFP(t.fNSigmaTOFP),
  fNSigmaTOFE(t.fNSigmaTOFE),
  fNSigmaTPCD(t.fNSigmaTPCD),
  fNSigmaTPCT(t.fNSigmaTPCT),
  fNSigmaTPCH(t.fNSigmaTPCH),
  fNSigmaTPCA(t.fNSigmaTPCA),
  fNSigmaTOFD(t.fNSigmaTOFD),
  fNSigmaTOFT(t.fNSigmaTOFT),
  fNSigmaTOFH(t.fNSigmaTOFH),
  fNSigmaTOFA(t.fNSigmaTOFA),
  fMassTOF(t.fMassTOF),
  fSigmaToVertex(t.fSigmaToVertex),
  fClusters(t.fClusters),
  fShared(t.fShared),
  fNominalTpcEntrancePoint(t.fNominalTpcEntrancePoint),
  fNominalTpcExitPoint(t.fNominalTpcExitPoint),
  fNominalTpcPointShifted(t.fNominalTpcPointShifted),
  fXatDCA(t.fXatDCA),
  fYatDCA(t.fYatDCA),
  fZatDCA(t.fZatDCA),
  fHiddenInfo(nullptr),
  fTrueMomentum(nullptr),
  fEmissionPoint(nullptr),
  fPDGPid(t.fPDGPid),
  fMass(t.fMass),
  fGlobalEmissionPoint(nullptr),
  fCorrPi(t.fCorrPi),
  fCorrK(t.fCorrK),
  fCorrP(t.fCorrP),
  fCorrPiMinus(t.fCorrPiMinus),
  fCorrKMinus(t.fCorrKMinus),
  fCorrPMinus(t.fCorrPMinus),
  fCorrAll(t.fCorrAll),
  fCorrD(t.fCorrD),
  fCorrT(t.fCorrT),
  fCorrH(t.fCorrH),
  fCorrA(t.fCorrA),
  fCorrDMinus(t.fCorrDMinus),
  fCorrTMinus(t.fCorrTMinus),
  fCorrHMinus(t.fCorrHMinus)
  , fCorrAMinus(t.fCorrAMinus)
 {
   // copy constructor
  fHiddenInfo = t.ValidHiddenInfo() ? t.GetHiddenInfo()->Clone() : nullptr;

  std::copy_n(t.fNominalTpcPoints, 9, fNominalTpcPoints);

  if (t.fTrueMomentum) {
    fTrueMomentum = new AliFemtoThreeVector(*t.fTrueMomentum);
  }

  if (t.fEmissionPoint) {
    fEmissionPoint = new AliFemtoLorentzVector(*t.fEmissionPoint);
  }

  if (t.fGlobalEmissionPoint) {
    fGlobalEmissionPoint = new AliFemtoThreeVector(*t.fGlobalEmissionPoint);
  }

  SetKinkIndexes(t.fKinkIndexes);
  SetPrimaryVertex(t.fVertex);
}

AliFemtoTrack& AliFemtoTrack::operator=(const AliFemtoTrack& aTrack)
{
  // assignment operator
  if (this == &aTrack)
    return *this;
  fCharge = aTrack.fCharge;
  fPidProbElectron = aTrack.fPidProbElectron;
  fPidProbPion = aTrack.fPidProbPion;
  fPidProbKaon = aTrack.fPidProbKaon;
  fPidProbProton = aTrack.fPidProbProton;
  fPidProbMuon=aTrack.fPidProbMuon;
  fTofPionTime=aTrack.fTofPionTime;
  fTofKaonTime=aTrack.fTofKaonTime;
  fTofProtonTime=aTrack.fTofProtonTime;
  fTofDeuteronTime=aTrack.fTofDeuteronTime;
  fTofTritonTime=aTrack.fTofTritonTime;
  fTofHe3Time=aTrack.fTofHe3Time;
  fTofAlphaTime=aTrack.fTofAlphaTime;
  fP=aTrack.fP;
  fPt=aTrack.fPt;
  fInnerMomentum=aTrack.fInnerMomentum;
  fMultiplicity=aTrack.fMultiplicity;
  fZvtx=aTrack.fZvtx;
  fHelix=aTrack.fHelix;
  fTrackId=aTrack.fTrackId;
  fFlags=aTrack.fFlags;
  fLabel=aTrack.fLabel;
  fImpactD=aTrack.fImpactD;
  fImpactDprim=aTrack.fImpactDprim;
  fImpactDweak=aTrack.fImpactDweak;
  fImpactDmat=aTrack.fImpactDmat;
  fImpactZ=aTrack.fImpactZ;
  fCdd=aTrack.fCdd;
  fCdz=aTrack.fCdz;
  fCzz=aTrack.fCzz;
  fITSchi2=aTrack.fITSchi2;
  fITSncls=aTrack.fITSncls;
  fTPCchi2=aTrack.fTPCchi2;
  fTPCncls=aTrack.fTPCncls;
  fTPCnclsF=aTrack.fTPCnclsF;
  fTPCsignal=aTrack.fTPCsignal;
  fTPCsignalN=aTrack.fTPCsignalN;
  fTPCsignalS=aTrack.fTPCsignalS;
  fVTOF=aTrack.fVTOF;
  fTOFsignal=aTrack.fTOFsignal;
  fNSigmaTPCPi=aTrack.fNSigmaTPCPi;
  fNSigmaTPCK=aTrack.fNSigmaTPCK;
  fNSigmaTPCP=aTrack.fNSigmaTPCP;
  fNSigmaTPCE=aTrack.fNSigmaTPCE;
  fNSigmaTPCD=aTrack.fNSigmaTPCD;
  fNSigmaTPCT=aTrack.fNSigmaTPCT;
  fNSigmaTPCH=aTrack.fNSigmaTPCH;
  fNSigmaTPCA=aTrack.fNSigmaTPCA;
  fMassTOF=aTrack.fMassTOF;
  fNSigmaTOFPi=aTrack.fNSigmaTOFPi;
  fNSigmaTOFK=aTrack.fNSigmaTOFK;
  fNSigmaTOFP=aTrack.fNSigmaTOFP;
  fNSigmaTOFE=aTrack.fNSigmaTOFE;
  fNSigmaTOFD=aTrack.fNSigmaTOFD;
  fNSigmaTOFT=aTrack.fNSigmaTOFT;
  fNSigmaTOFH=aTrack.fNSigmaTOFH;
  fNSigmaTOFA=aTrack.fNSigmaTOFA;
  fClusters=aTrack.fClusters;
  fShared=aTrack.fShared;
  fNominalTpcEntrancePoint=aTrack.fNominalTpcEntrancePoint;
  fNominalTpcExitPoint=aTrack.fNominalTpcExitPoint;
  fNominalTpcPointShifted=aTrack.fNominalTpcPointShifted;

  fMass = aTrack.fMass;
  fPDGPid = aTrack.fPDGPid;

  fXatDCA=aTrack.fXatDCA;
  fYatDCA=aTrack.fYatDCA;
  fZatDCA=aTrack.fZatDCA;

  fCorrPi = aTrack.fCorrPi;
  fCorrK = aTrack.fCorrK;
  fCorrP = aTrack.fCorrP;
  fCorrD = aTrack.fCorrD;
  fCorrT = aTrack.fCorrT;
  fCorrH = aTrack.fCorrH;
  fCorrA = aTrack.fCorrA;
  fCorrPiMinus = aTrack.fCorrPiMinus;
  fCorrKMinus = aTrack.fCorrKMinus;
  fCorrPMinus = aTrack.fCorrPMinus;
  fCorrDMinus = aTrack.fCorrDMinus;
  fCorrTMinus = aTrack.fCorrTMinus;
  fCorrHMinus = aTrack.fCorrHMinus;
  fCorrAMinus = aTrack.fCorrAMinus;
  fCorrAll = aTrack.fCorrAll;

  std::copy_n(aTrack.fHasPointOnITS, 6, fHasPointOnITS);
  std::copy_n(aTrack.fNominalTpcPoints, 9, fNominalTpcPoints);

  delete fHiddenInfo;
  fHiddenInfo = aTrack.ValidHiddenInfo()
              ? aTrack.GetHiddenInfo()->Clone()
              : nullptr;

  if (aTrack.fTrueMomentum && fTrueMomentum) {
    *fTrueMomentum = *aTrack.fTrueMomentum;
  } else if (aTrack.fTrueMomentum) {
    fTrueMomentum = new AliFemtoThreeVector(*aTrack.fTrueMomentum);
  } else {
    delete fTrueMomentum;
    fTrueMomentum = nullptr;
  }

  if (aTrack.fEmissionPoint && fEmissionPoint) {
    *fEmissionPoint = *aTrack.fEmissionPoint;
  } else if (aTrack.fEmissionPoint) {
    fEmissionPoint = new AliFemtoLorentzVector(*aTrack.fEmissionPoint);
  } else {
    delete fEmissionPoint;
    fEmissionPoint = nullptr;
  }

  if (aTrack.fGlobalEmissionPoint && fGlobalEmissionPoint) {
    *fGlobalEmissionPoint = *aTrack.fGlobalEmissionPoint;
  } else if (aTrack.fGlobalEmissionPoint) {
    fGlobalEmissionPoint = new AliFemtoThreeVector(*aTrack.fGlobalEmissionPoint);
  } else {
    delete fGlobalEmissionPoint;
    fGlobalEmissionPoint = NULL;
  }

  SetKinkIndexes(aTrack.fKinkIndexes);
  SetPrimaryVertex(aTrack.fVertex);

  return *this;
}


AliFemtoTrack::~AliFemtoTrack()
{
  // destructor
  delete fHiddenInfo;
  delete fTrueMomentum;
  delete fEmissionPoint;
  delete fGlobalEmissionPoint;
}

// void AliFemtoTrack::SetXTPC(const AliFemtoThreeVector& aXTPC)
// {
//   fXTPC = aXTPC;
// }

// void AliFemtoTrack::SetXTPC(double *aXTPC)
// {
//   fXTPC.setX(aXTPC[0]);
//   fXTPC.setY(aXTPC[1]);
//   fXTPC.setZ(aXTPC[2]);
// }

// AliFemtoThreeVector AliFemtoTrack::XTPC() const
// {
//   return fXTPC;
// }
