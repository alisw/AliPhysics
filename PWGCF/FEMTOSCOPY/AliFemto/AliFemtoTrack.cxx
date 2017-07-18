///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoTrack: main class holding all the necessary information       //
// about a track (before the identification) that is required during     //
// femtoscopic analysis. This class is filled with information from the  //
// input stream by the reader. A particle has a link back to the Track   //
// it was created from, so we do not copy the information.               //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "AliFemtoTrack.h"

#include "SystemOfUnits.h"   // has "tesla" in it
//#include "AliFemtoParticleTypes.h"
//#include "AliFemtoTTreeEvent.h"
//#include "AliFemtoTTreeTrack.h"

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
  fNSigmaTPCPi(0.0f),
  fNSigmaTPCK(0.0f),
  fNSigmaTPCP(0.0f),
  fNSigmaTPCE(0.0f),
  fNSigmaTOFPi(0.0f),
  fNSigmaTOFK(0.0f),
  fNSigmaTOFP(0.0f),
  fNSigmaTOFE(0.0f),
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
  fTrueMomentum(nullptr),    // True (simulated) momentum
  fEmissionPoint(nullptr),   // Emission point coordinates
  fPDGPid(0),             // True PID of the particle
  fMass(0.0),               // True particle mass
  fGlobalEmissionPoint(nullptr),
  fCorrPi(0.0),
  fCorrK(0.0),
  fCorrP(0.0),
  fCorrPiMinus(0.0),
  fCorrKMinus(0.0),
  fCorrPMinus(0.0),
  fCorrAll(0.0)
{
  // Default constructor
  fKinkIndexes[0] = 0;
  fKinkIndexes[1] = 0;
  fKinkIndexes[2] = 0;

  for (int i = 0; i < 6; i++) {
    fHasPointOnITS[i] = kFALSE;
  }

  for (int i = 0; i < 9; i++) {
    fNominalTpcPoints[i].SetX(0);
    fNominalTpcPoints[i].SetY(0);
    fNominalTpcPoints[i].SetZ(0);
  }

  fVertex[0] = -9999;
  fVertex[1] = -9999;
  fVertex[2] = -9999;
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
  fNSigmaTPCPi(t.fNSigmaTPCPi),
  fNSigmaTPCK(t.fNSigmaTPCK),
  fNSigmaTPCP(t.fNSigmaTPCP),
  fNSigmaTPCE(t.fNSigmaTPCE),
  fNSigmaTOFPi(t.fNSigmaTOFPi),
  fNSigmaTOFK(t.fNSigmaTOFK),
  fNSigmaTOFP(t.fNSigmaTOFP),
  fNSigmaTOFE(t.fNSigmaTOFE),
  fSigmaToVertex(t.fSigmaToVertex),
  fClusters(t.fClusters),
  fShared(t.fShared),
  fNominalTpcEntrancePoint(t.fNominalTpcEntrancePoint),
  fNominalTpcExitPoint(t.fNominalTpcExitPoint),
  fNominalTpcPointShifted(t.fNominalTpcPointShifted),
  fXatDCA(t.fXatDCA),
  fYatDCA(t.fYatDCA),
  fZatDCA(t.fZatDCA),
  fHiddenInfo(NULL),
  fTrueMomentum(NULL),      // True (simulated) momentum
  fEmissionPoint(NULL),     // Emission point coordinates
  fPDGPid(t.fPDGPid),       // True PID of the particle
  fMass(t.fMass),           // True particle mass
  fGlobalEmissionPoint(0),
  fCorrPi(t.fCorrPi),
  fCorrK(t.fCorrK),
  fCorrP(t.fCorrP),
  fCorrPiMinus(t.fCorrPiMinus),
  fCorrKMinus(t.fCorrKMinus),
  fCorrPMinus(t.fCorrPMinus),
  fCorrAll(t.fCorrAll)
 {
   // copy constructor
  fHiddenInfo = t.ValidHiddenInfo() ? t.GetHiddenInfo()->Clone() : nullptr;

  for (int i = 0; i < 9; i++) {
    fNominalTpcPoints[i] = t.fNominalTpcPoints[i];
  }

  if (t.fTrueMomentum) {
    fTrueMomentum = new AliFemtoThreeVector(*t.fTrueMomentum);
  }

  if (t.fEmissionPoint) {
    fEmissionPoint = new AliFemtoLorentzVector(*t.fEmissionPoint);
  }

  if (t.fGlobalEmissionPoint) {
    fGlobalEmissionPoint = new AliFemtoThreeVector(*t.fGlobalEmissionPoint);
  }

  SetKinkIndexes((int*)t.fKinkIndexes);
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
  fNSigmaTPCPi=aTrack.fNSigmaTPCPi;
  fNSigmaTPCK=aTrack.fNSigmaTPCK;
  fNSigmaTPCP=aTrack.fNSigmaTPCP;
  fNSigmaTPCE=aTrack.fNSigmaTPCE;
  fNSigmaTOFPi=aTrack.fNSigmaTOFPi;
  fNSigmaTOFK=aTrack.fNSigmaTOFK;
  fNSigmaTOFP=aTrack.fNSigmaTOFP;
  fNSigmaTOFE=aTrack.fNSigmaTOFE;
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
  fCorrPiMinus = aTrack.fCorrPiMinus;
  fCorrKMinus = aTrack.fCorrKMinus;
  fCorrPMinus = aTrack.fCorrPMinus;
  fCorrAll = aTrack.fCorrAll;

  for (int i=0; i<6; i++) {
    fHasPointOnITS[i] = aTrack.fHasPointOnITS[i];
  }

  for (int i=0; i<9; i++) {
    fNominalTpcPoints[i] = aTrack.fNominalTpcPoints[i];
  }

  delete fHiddenInfo;
  fHiddenInfo = (aTrack.ValidHiddenInfo()) ? aTrack.GetHiddenInfo()->Clone() : NULL;

  if (aTrack.fTrueMomentum != NULL) {
    if (fTrueMomentum == NULL) {
      fTrueMomentum = new AliFemtoThreeVector(*aTrack.fTrueMomentum);
    } else {
      *fTrueMomentum = *aTrack.fTrueMomentum;
    }
  } else {
    delete fTrueMomentum;
    fTrueMomentum = NULL;
  }

  if (aTrack.fEmissionPoint != NULL) {
    if (fEmissionPoint == NULL) {
      fEmissionPoint = new AliFemtoLorentzVector(*aTrack.fEmissionPoint);
    } else {
      *fEmissionPoint = *aTrack.fEmissionPoint;
    }
  } else {
    delete fEmissionPoint;
    fEmissionPoint = NULL;
  }

  if (aTrack.fGlobalEmissionPoint != NULL) {
    if (fGlobalEmissionPoint == NULL) {
      fGlobalEmissionPoint = new AliFemtoThreeVector(*aTrack.fGlobalEmissionPoint);
    } else {
      *fGlobalEmissionPoint = *aTrack.fGlobalEmissionPoint;
    }
  } else {
    delete fGlobalEmissionPoint;
    fGlobalEmissionPoint = NULL;
  }

  SetKinkIndexes((int*)aTrack.fKinkIndexes);
  SetPrimaryVertex((double*)aTrack.fVertex);

  return *this;
}

void AliFemtoTrack::SetCharge(const short& ch){fCharge=ch;}

void AliFemtoTrack::SetPidProbElectron(const float& x){fPidProbElectron = x;}
void AliFemtoTrack::SetPidProbPion(const float& x){fPidProbPion = x;}
void AliFemtoTrack::SetPidProbKaon(const float& x){fPidProbKaon = x;}
void AliFemtoTrack::SetPidProbProton(const float& x){fPidProbProton = x;}
void AliFemtoTrack::SetPidProbMuon(const float& x){fPidProbMuon = x;}
void AliFemtoTrack::SetTofExpectedTimes(const float& tpi, const float& tkn, const float& tpr){fTofPionTime = tpi; fTofKaonTime = tkn; fTofProtonTime = tpr; }

void AliFemtoTrack::SetP(const AliFemtoThreeVector& p){fP = p;}
void AliFemtoTrack::SetPt(const float& pt){fPt = pt;}
void AliFemtoTrack::SetInnerMomentum(const float& x){fInnerMomentum = x;}
void AliFemtoTrack::SetHelix(const AliFmPhysicalHelixD& h){fHelix = h;}
void AliFemtoTrack::SetTrackId(const int & id) { fTrackId=id;}
void AliFemtoTrack::SetFlags(const long int &flags) {fFlags=flags;}
void AliFemtoTrack::SetLabel(const int &label) {fLabel=label;}
void AliFemtoTrack::SetImpactD(const float& aImpactD){fImpactD=aImpactD;}

void AliFemtoTrack::SetImpactDprim(const float& aImpactDprim){fImpactDprim=aImpactDprim;}
void AliFemtoTrack::SetImpactDweak(const float& aImpactDweak){fImpactDweak=aImpactDweak;}
void AliFemtoTrack::SetImpactDmat(const float& aImpactDmat){fImpactDmat=aImpactDmat;}

void AliFemtoTrack::SetImpactZ(const float& aImpactZ){fImpactZ=aImpactZ;}
void AliFemtoTrack::SetCdd(const float& aCdd){fCdd=aCdd;}
void AliFemtoTrack::SetCdz(const float& aCdz){fCdz=aCdz;}
void AliFemtoTrack::SetCzz(const float& aCzz){fCzz=aCzz;}
void AliFemtoTrack::SetITSchi2(const float& aITSchi2){fITSchi2=aITSchi2;}
void AliFemtoTrack::SetITSncls(const int& aITSncls){fITSncls=aITSncls;}
void AliFemtoTrack::SetTPCchi2(const float& aTPCchi2){fTPCchi2=aTPCchi2;}
void AliFemtoTrack::SetTPCncls(const int& aTPCncls){fTPCncls=aTPCncls;}
void AliFemtoTrack::SetTPCnclsF(const short& aTPCnclsF){fTPCnclsF=aTPCnclsF;}
void AliFemtoTrack::SetTPCsignal(const float& aTPCsig){fTPCsignal=aTPCsig;}
void AliFemtoTrack::SetTPCsignalN(const short& aTPCsignalN){fTPCsignalN=aTPCsignalN;}
void AliFemtoTrack::SetTPCsignalS(const float& aTPCsignalS){fTPCsignalS=aTPCsignalS;}
void AliFemtoTrack::SetVTOF(const float& aVTOF){fVTOF=aVTOF;}
void AliFemtoTrack::SetNSigmaTPCPi(const float& aNSigmaTPCPi){fNSigmaTPCPi=aNSigmaTPCPi;}
void AliFemtoTrack::SetNSigmaTPCK(const float& aNSigmaTPCK){fNSigmaTPCK=aNSigmaTPCK;}
void AliFemtoTrack::SetNSigmaTPCP(const float& aNSigmaTPCP){fNSigmaTPCP=aNSigmaTPCP;}
void AliFemtoTrack::SetNSigmaTPCE(const float& aNSigmaTPCE){fNSigmaTPCE=aNSigmaTPCE;}
void AliFemtoTrack::SetNSigmaTOFPi(const float& aNSigmaTOFPi){fNSigmaTOFPi=aNSigmaTOFPi;}
void AliFemtoTrack::SetNSigmaTOFK(const float& aNSigmaTOFK){fNSigmaTOFK=aNSigmaTOFK;}
void AliFemtoTrack::SetNSigmaTOFP(const float& aNSigmaTOFP){fNSigmaTOFP=aNSigmaTOFP;}
void AliFemtoTrack::SetNSigmaTOFE(const float& aNSigmaTOFE){fNSigmaTOFE=aNSigmaTOFE;}
void AliFemtoTrack::SetSigmaToVertex(const float& aSigma){fSigmaToVertex=aSigma;}

void AliFemtoTrack::SetXatDCA(const double& x) {fXatDCA=x;}
void AliFemtoTrack::SetYatDCA(const double& x) {fYatDCA=x;}
void AliFemtoTrack::SetZatDCA(const double& x) {fZatDCA=x;}


void AliFemtoTrack::SetCorrectionPion(const double& x){fCorrPi=x;}
void AliFemtoTrack::SetCorrectionKaon(const double& x){fCorrK=x;}
void AliFemtoTrack::SetCorrectionProton(const double& x){fCorrP=x;}
void AliFemtoTrack::SetCorrectionPionMinus(const double& x){fCorrPiMinus=x;}
void AliFemtoTrack::SetCorrectionKaonMinus(const double& x){fCorrKMinus=x;}
void AliFemtoTrack::SetCorrectionProtonMinus(const double& x){fCorrPMinus=x;}
void AliFemtoTrack::SetCorrectionAll(const double& x){fCorrAll=x;}

short AliFemtoTrack::Charge() const {return fCharge;}
AliFemtoThreeVector AliFemtoTrack::P() const {return fP;}
float AliFemtoTrack::Pt() const {return fPt;}
float AliFemtoTrack::InnerMomentum() const {return fInnerMomentum;}
const AliFmPhysicalHelixD& AliFemtoTrack::Helix() const {return fHelix;}
int AliFemtoTrack::TrackId() const { return fTrackId; }
long int AliFemtoTrack::Flags() const {return fFlags;}
int AliFemtoTrack::Label()const {return fLabel;}
float AliFemtoTrack::ImpactD()const{return fImpactD;}

float AliFemtoTrack::ImpactDprim()const{return fImpactDprim;}
float AliFemtoTrack::ImpactDweak()const{return fImpactDweak;}
float AliFemtoTrack::ImpactDmat()const{return fImpactDmat;}

float AliFemtoTrack::ImpactZ()const{return fImpactZ;}
float AliFemtoTrack::Cdd() const{return fCdd;}
float AliFemtoTrack::Cdz() const{return fCdz;}
float AliFemtoTrack::Czz() const{return fCzz;}
float AliFemtoTrack::ITSchi2() const{return fITSchi2;}
int   AliFemtoTrack::ITSncls() const{return fITSncls;}
float AliFemtoTrack::TPCchi2() const{return fTPCchi2;}
int   AliFemtoTrack::TPCncls() const{return fTPCncls;}
short AliFemtoTrack::TPCnclsF() const{return fTPCnclsF;}
float AliFemtoTrack::TPCsignal() const{return fTPCsignal;}
short AliFemtoTrack::TPCsignalN() const{return fTPCsignalN;}
float AliFemtoTrack::TPCsignalS() const{return fTPCsignalS;}
float AliFemtoTrack::VTOF() const{return fVTOF;}
float AliFemtoTrack::NSigmaTPCPi() const{return fNSigmaTPCPi;}
float AliFemtoTrack::NSigmaTPCK() const{return fNSigmaTPCK;}
float AliFemtoTrack::NSigmaTPCP() const{return fNSigmaTPCP;}
float AliFemtoTrack::NSigmaTPCE() const{return fNSigmaTPCE;}
float AliFemtoTrack::NSigmaTOFPi() const{return fNSigmaTOFPi;}
float AliFemtoTrack::NSigmaTOFK() const{return fNSigmaTOFK;}
float AliFemtoTrack::NSigmaTOFP() const{return fNSigmaTOFP;}
float AliFemtoTrack::NSigmaTOFE() const{return fNSigmaTOFE;}
float AliFemtoTrack::SigmaToVertex() const{return fSigmaToVertex;}
float AliFemtoTrack::TOFpionTime() const{return fTofPionTime;}
float AliFemtoTrack::TOFkaonTime() const{return fTofKaonTime;}
float AliFemtoTrack::TOFprotonTime() const{return fTofProtonTime;}
//corrections
float AliFemtoTrack::CorrectionPion() const {return fCorrPi;}
float AliFemtoTrack::CorrectionKaon() const {return fCorrK;}
float AliFemtoTrack::CorrectionProton() const {return fCorrP;}
float AliFemtoTrack::CorrectionPionMinus() const {return fCorrPiMinus;}
float AliFemtoTrack::CorrectionKaonMinus() const {return fCorrKMinus;}
float AliFemtoTrack::CorrectionProtonMinus() const {return fCorrPMinus;}
float AliFemtoTrack::CorrectionAll() const {return fCorrAll;}

double AliFemtoTrack::XatDCA() const {return fXatDCA;}
double AliFemtoTrack::YatDCA() const {return fYatDCA;}
double AliFemtoTrack::ZatDCA() const {return fZatDCA;}

void AliFemtoTrack::SetHiddenInfo(AliFemtoHiddenInfo* aHiddenInfo) {fHiddenInfo=aHiddenInfo;}
bool AliFemtoTrack::ValidHiddenInfo() const { return (fHiddenInfo != NULL); }
AliFemtoHiddenInfo* AliFemtoTrack::GetHiddenInfo() const {return fHiddenInfo;}

AliFemtoTrack::~AliFemtoTrack()
{
  // destructor
  delete fHiddenInfo;
  delete fTrueMomentum;
  delete fEmissionPoint;
  delete fGlobalEmissionPoint;
}

const TBits& AliFemtoTrack::TPCclusters() const {return fClusters;}
const TBits& AliFemtoTrack::TPCsharing()  const {return fShared;}

void AliFemtoTrack::SetTPCcluster(const short& aNBit, const Bool_t& aValue)
{
  fClusters.SetBitNumber(aNBit, aValue);
}

void AliFemtoTrack::SetTPCshared(const short& aNBit, const Bool_t& aValue)
{
  fShared.SetBitNumber(aNBit, aValue);
}

void AliFemtoTrack::SetTPCClusterMap(const TBits& aBits)
{
  fClusters = aBits;
}
void AliFemtoTrack::SetTPCSharedMap(const TBits& aBits)
{
  fShared = aBits;
}

void AliFemtoTrack::SetKinkIndexes(int points[3])
{
  // Transfer the Kink indices
  fKinkIndexes[0] = points[0];
  fKinkIndexes[1] = points[1];
  fKinkIndexes[2] = points[2];
}

void AliFemtoTrack::SetITSHitOnLayer(int i, bool val)
{
  // Transfer ITS hit
  fHasPointOnITS[i] = val;
}


int  AliFemtoTrack::KinkIndex(int aIndex) const
{
  // Return Kink index
  if ((aIndex < 3) && (aIndex >= 0))
    return fKinkIndexes[aIndex];
  else
    return 0;
}

bool AliFemtoTrack::HasPointOnITSLayer(int aIndex) const
{
  // Return if i-th ITS layer had a hit for this track
  if ((aIndex < 6) && (aIndex >= 0))
    return fHasPointOnITS[aIndex];
  else
    return false;
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

const AliFemtoThreeVector& AliFemtoTrack::NominalTpcExitPoint() const
{
  return fNominalTpcExitPoint;
}

const AliFemtoThreeVector& AliFemtoTrack::NominalTpcPointShifted() const
{
  return fNominalTpcPointShifted;
}

const AliFemtoThreeVector& AliFemtoTrack::NominalTpcPoint(int i) const
{
  if(i<0)
    return fNominalTpcPoints[0];
  if(i>8)
    return fNominalTpcPoints[8];
  return fNominalTpcPoints[i];
}

const AliFemtoThreeVector& AliFemtoTrack::NominalTpcEntrancePoint() const
{
  return fNominalTpcEntrancePoint;
}

void AliFemtoTrack::SetNominalTPCEntrancePoint(const AliFemtoThreeVector& aXTPC)
{
  fNominalTpcEntrancePoint = aXTPC;
}
void AliFemtoTrack::SetNominalTPCEntrancePoint(double *aXTPC)
{
  // Store the nominal TPC entrance point
  fNominalTpcEntrancePoint.SetX(aXTPC[0]);
  fNominalTpcEntrancePoint.SetY(aXTPC[1]);
  fNominalTpcEntrancePoint.SetZ(aXTPC[2]);
}

void AliFemtoTrack::SetNominalTPCPoints(double **aXTPC)
{
  // Store the nominal TPC points
  for(int i=0;i<9;i++)
    {
      fNominalTpcPoints[i].SetX(aXTPC[i][0]);
      fNominalTpcPoints[i].SetY(aXTPC[i][1]);
      fNominalTpcPoints[i].SetZ(aXTPC[i][2]);
    }
}

void AliFemtoTrack::SetNominalTPCExitPoint(const AliFemtoThreeVector& aXTPC)
{
  fNominalTpcExitPoint = aXTPC;
}
void AliFemtoTrack::SetNominalTPCExitPoint(double *aXTPC)
{
  // Store the nominal TPC exit point
  fNominalTpcExitPoint.SetX(aXTPC[0]);
  fNominalTpcExitPoint.SetY(aXTPC[1]);
  fNominalTpcExitPoint.SetZ(aXTPC[2]);
}

void AliFemtoTrack::SetNominalTPCPointShifted(const AliFemtoThreeVector& aXTPC)
{
  fNominalTpcPointShifted = aXTPC;
}

void AliFemtoTrack::SetNominalTPCPointShifted(double *aXTPC)
{
  fNominalTpcPointShifted.SetX(aXTPC[0]);
  fNominalTpcPointShifted.SetY(aXTPC[1]);
  fNominalTpcPointShifted.SetZ(aXTPC[2]);
}



//_____________________________________________
AliFemtoThreeVector   *AliFemtoTrack::GetTrueMomentum() const
{
  return fTrueMomentum;
}
//_____________________________________________
AliFemtoLorentzVector *AliFemtoTrack::GetEmissionPoint() const
{
  return fEmissionPoint;
}
//_____________________________________________
Int_t                  AliFemtoTrack::GetPDGPid() const
{
  return fPDGPid;
}
//_____________________________________________
Double_t                  AliFemtoTrack::GetMass() const
{
  return fMass;
}
//_____________________________________________
void                   AliFemtoTrack::SetTrueMomentum(AliFemtoThreeVector *aMom)
{
  // Set momentum from vector
  if (fTrueMomentum) {
    *fTrueMomentum = *aMom;
  }
  else {
    fTrueMomentum = new AliFemtoThreeVector(*aMom);
  }
}
//_____________________________________________
void                   AliFemtoTrack::SetTrueMomentum(const AliFemtoThreeVector& aMom)
{
  // Set momentum from vector
  if (fTrueMomentum) {
    *fTrueMomentum = aMom;
  }
  else {
    fTrueMomentum = new AliFemtoThreeVector(aMom);
  }
}
//_____________________________________________
void                   AliFemtoTrack::SetTrueMomentum(Double_t aPx, Double_t aPy, Double_t aPz)
{
  // Set momentum from components
  if (fTrueMomentum) {
    fTrueMomentum->SetX(aPx);
    fTrueMomentum->SetY(aPy);
    fTrueMomentum->SetZ(aPz);
  }
  else {
    fTrueMomentum = new AliFemtoThreeVector(aPx, aPy, aPz);
  }
}
//_____________________________________________
void                   AliFemtoTrack::SetEmissionPoint(AliFemtoLorentzVector *aPos)
{
  // Set position from vector
  if (fEmissionPoint) {
    *fEmissionPoint = *aPos;
  }
  else {
    fEmissionPoint = new AliFemtoLorentzVector(*aPos);
  }
}
//_____________________________________________
void                   AliFemtoTrack::SetEmissionPoint(const AliFemtoLorentzVector& aPos)
{
  // Set position from vector
  if (fEmissionPoint) {
    *fEmissionPoint = aPos;
  }
  else {
    fEmissionPoint = new AliFemtoLorentzVector(aPos);
  }
}
//_____________________________________________
void                   AliFemtoTrack::SetPDGPid(Int_t aPid)
{
  fPDGPid = aPid;
}
//_____________________________________________
void                   AliFemtoTrack::SetMass(Double_t aMass)
{
  fMass = aMass;
}
//_____________________________________________
void                   AliFemtoTrack::SetEmissionPoint(Double_t aRx, Double_t aRy, Double_t aRz, Double_t aT)
{
  // Set position from components
  if (fEmissionPoint) {
    fEmissionPoint->SetX(aRx);
    fEmissionPoint->SetY(aRy);
    fEmissionPoint->SetZ(aRz);
    fEmissionPoint->SetT(aT);
  }
  else {
    fEmissionPoint = new AliFemtoLorentzVector(aRx, aRy, aRz, aT);
  }
}

//_____________________________________________
AliFemtoThreeVector *AliFemtoTrack::GetGlobalEmissionPoint() const
{
  return fGlobalEmissionPoint;
}
//_____________________________________________
void                   AliFemtoTrack::SetGlobalEmissionPoint(const AliFemtoThreeVector& aPos)
{
  // set position from vector
  if (fGlobalEmissionPoint) {
    *fGlobalEmissionPoint = aPos;
  } else {
    fGlobalEmissionPoint = new AliFemtoThreeVector(aPos);
  }
}
//_____________________________________________
void                   AliFemtoTrack::SetGlobalEmissionPoint(Double_t aRx, Double_t aRy, Double_t aRz)
{
  // Set position from components
  if (fGlobalEmissionPoint) {
    fGlobalEmissionPoint->SetX(aRx);
    fGlobalEmissionPoint->SetY(aRy);
    fGlobalEmissionPoint->SetZ(aRz);
  }
  else {
    fGlobalEmissionPoint = new AliFemtoThreeVector(aRx, aRy, aRz);
  }
}
//_______________________


void AliFemtoTrack::SetPrimaryVertex(const double* vertex)
{
  fVertex[0] = vertex[0];
  fVertex[1] = vertex[1];
  fVertex[2] = vertex[2];
}

void AliFemtoTrack::GetPrimaryVertex(double* vertex)
{
  vertex[0] = fVertex[0];
  vertex[1] = fVertex[1];
  vertex[2] = fVertex[2];
}
//______________________
void AliFemtoTrack::SetMultiplicity(int mult)
{
  fMultiplicity=mult;
}
//______________________
void AliFemtoTrack::SetZvtx(double vtx)
{
  fZvtx=vtx;
}
//______________________
