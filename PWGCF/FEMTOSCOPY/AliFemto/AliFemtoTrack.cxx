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
#ifndef __NO_STAR_DEPENDENCE_ALLOWED__
#ifdef __ROOT__
#include "StEvent/StEnumerations.h"
#include "AliFemtoAihongPid.h"
#include "StEventUtilities/StuProbabilityPidAlgorithm.h"
#endif
#endif

#include "SystemOfUnits.h"   // has "tesla" in it
//#include "AliFemtoParticleTypes.h"
//#include "AliFemtoTTreeEvent.h" 
//#include "AliFemtoTTreeTrack.h" 

AliFemtoTrack::AliFemtoTrack():
  fCharge(0),
  fPidProbElectron(0),
  fPidProbPion(0),
  fPidProbKaon(0),
  fPidProbProton(0),
  fPidProbMuon(0),
  fTrackId(0),
  fTofPionTime(-100000.0),
  fTofKaonTime(-100000.0),
  fTofProtonTime(-100000.0),
  fP(0,0,0),
  fPt(0),
  fInnerMomentum(0),
  fHelix(),
  fFlags(0),
  fLabel(0),
  fImpactD(0),
  fImpactDprim(-10000.0),
  fImpactDweak(-10000.0),
  fImpactDmat(-10000.0),
  fImpactZ(0),
  fCdd(0),
  fCdz(0),
  fCzz(0),
  fITSchi2(0),       
  fITSncls(0),        
  fTPCchi2(0),       
  fTPCncls(0),       
  fTPCnclsF(0),      
  fTPCsignal(0),
  fTPCsignalN(0),    
  fTPCsignalS(0),
  fVTOF(0),
  fNSigmaTPCPi(0),
  fNSigmaTPCK(0),
  fNSigmaTPCP(0),
  fNSigmaTOFPi(0),
  fNSigmaTOFK(0),
  fNSigmaTOFP(0),
  fSigmaToVertex(0),
  fClusters(159),
  fShared(159),
  fNominalTpcEntrancePoint(0,0,0),
  fNominalTpcExitPoint(0,0,0),
  fXatDCA(0),
  fYatDCA(0),
  fZatDCA(0),
  fHiddenInfo(0),
  fTrueMomentum(0),  // True (simulated) momentum
  fEmissionPoint(0), // Emission point coordinates
  fPDGPid(0),        // True PID of the particle
  fMass(0),          // True particle mass
  fGlobalEmissionPoint(0)
{
  // Default constructor
  fHiddenInfo = NULL;
  fKinkIndexes[0] = 0;
  fKinkIndexes[1] = 0;
  fKinkIndexes[2] = 0;

  for(int i=0;i<6;i++) {
    fHasPointOnITS[i]=false;
  }

  for(int i=0;i<9;i++)
    {
      fNominalTpcPoints[i].SetX(0);
      fNominalTpcPoints[i].SetY(0);
      fNominalTpcPoints[i].SetZ(0);
    }
  //  cout << "Created track " << this << endl;
}


AliFemtoTrack::AliFemtoTrack(const AliFemtoTrack& t) :
  fCharge(0),
  fPidProbElectron(0),
  fPidProbPion(0),
  fPidProbKaon(0),
  fPidProbProton(0),
  fPidProbMuon(0),
  fTrackId(0),
  fTofPionTime(-100000.0),
  fTofKaonTime(-100000.0),
  fTofProtonTime(-100000.0),
  fP(0,0,0),
  fPt(0),
  fInnerMomentum(0),
  fHelix(),
  fFlags(0),
  fLabel(0),
  fImpactD(0),
  fImpactDprim(-10000.0),
  fImpactDweak(-10000.0),
  fImpactDmat(-10000.0),
  fImpactZ(0),
  fCdd(0),
  fCdz(0),
  fCzz(0),
  fITSchi2(0),       
  fITSncls(0),        
  fTPCchi2(0),       
  fTPCncls(0),       
  fTPCnclsF(0),      
  fTPCsignal(0),
  fTPCsignalN(0),    
  fTPCsignalS(0),
  fVTOF(0),
  fNSigmaTPCPi(0),
  fNSigmaTPCK(0),
  fNSigmaTPCP(0),
  fNSigmaTOFPi(0),
  fNSigmaTOFK(0),
  fNSigmaTOFP(0),
  fSigmaToVertex(0),
  fClusters(159),
  fShared(159),
  fNominalTpcEntrancePoint(0,0,0),
  fNominalTpcExitPoint(0,0,0),
  fXatDCA(0),
  fYatDCA(0),
  fZatDCA(0),
  fHiddenInfo(0),
  fTrueMomentum(0),  // True (simulated) momentum
  fEmissionPoint(0), // Emission point coordinates
  fPDGPid(0),        // True PID of the particle
  fMass(0),          // True particle mass
  fGlobalEmissionPoint(0)
 { 
   // copy constructor
  fCharge = t.fCharge;
  fPidProbElectron = t.fPidProbElectron;
  fPidProbPion = t.fPidProbPion;
  fPidProbKaon = t.fPidProbKaon;
  fPidProbProton = t.fPidProbProton;
  fPidProbMuon=t.fPidProbMuon;
  fTofPionTime=t.fTofPionTime;
  fTofKaonTime=t.fTofKaonTime;
  fTofProtonTime=t.fTofProtonTime;
  fP = t.fP;
  fPt = t.fPt;
  fInnerMomentum = t.fInnerMomentum;
  fHelix = t.fHelix;
  fTrackId = t.fTrackId;
  fFlags=t.fFlags;
  fLabel=t.fLabel;
  fImpactD=t.fImpactD;
  fImpactDprim=t.fImpactDprim;
  fImpactDweak=t.fImpactDweak;
  fImpactDmat=t.fImpactDmat;
  fImpactZ=t.fImpactZ;
  fCdd=t.fCdd;
  fCdz=t.fCdz;
  fCzz=t.fCzz;
  fITSchi2=t.fITSchi2;       
  fITSncls=t.fITSncls;        
  fTPCchi2=t.fTPCchi2;       
  fTPCncls=t.fTPCncls;       
  fTPCnclsF=t.fTPCnclsF;      
  fTPCsignal=t.fTPCsignal;
  fTPCsignalN=t.fTPCsignalN;    
  fTPCsignalS=t.fTPCsignalS;  
  fVTOF=t.fVTOF;
  fNSigmaTPCPi=t.fNSigmaTPCPi;
  fNSigmaTPCK=t.fNSigmaTPCK;
  fNSigmaTPCP=t.fNSigmaTPCP;
  fNSigmaTOFPi=t.fNSigmaTOFPi;
  fNSigmaTOFK=t.fNSigmaTOFK;
  fNSigmaTOFP=t.fNSigmaTOFP;
  fSigmaToVertex=t.fSigmaToVertex;
  fClusters=t.fClusters;
  fShared=t.fShared;
  fNominalTpcEntrancePoint=t.fNominalTpcEntrancePoint;
  fNominalTpcExitPoint=t.fNominalTpcExitPoint;
  if (t.ValidHiddenInfo())
    fHiddenInfo = t.GetHiddenInfo()->Clone();
  else 
    fHiddenInfo = NULL;
  fKinkIndexes[0] = t.fKinkIndexes[0];
  fKinkIndexes[1] = t.fKinkIndexes[1];
  fKinkIndexes[2] = t.fKinkIndexes[2];

  fXatDCA=t.fXatDCA;
  fYatDCA=t.fYatDCA;
  fZatDCA=t.fZatDCA;

  for(int i=0;i<9;i++)
    fNominalTpcPoints[i] = t.fNominalTpcPoints[i];


  fTrueMomentum = new AliFemtoThreeVector();
  if(t.fTrueMomentum){
    fTrueMomentum->SetX(t.fTrueMomentum->x());
    fTrueMomentum->SetY(t.fTrueMomentum->y());
    fTrueMomentum->SetZ(t.fTrueMomentum->z());}

  fEmissionPoint = new AliFemtoLorentzVector();
  if(t.fEmissionPoint){
    fEmissionPoint->SetX(t.fEmissionPoint->x());
    fEmissionPoint->SetY(t.fEmissionPoint->y());
    fEmissionPoint->SetZ(t.fEmissionPoint->z());
    fEmissionPoint->SetT(t.fEmissionPoint->e());
  }

  fPDGPid = t.fPDGPid;
  fMass = t.fMass;
 
  fGlobalEmissionPoint = new AliFemtoThreeVector();
  if(t.fGlobalEmissionPoint){
    fGlobalEmissionPoint->SetX(t.fGlobalEmissionPoint->x());
    fGlobalEmissionPoint->SetY(t.fGlobalEmissionPoint->y());
    fGlobalEmissionPoint->SetZ(t.fGlobalEmissionPoint->z());
    //fGlobalEmissionPoint->SetT(t.fGlobalEmissionPoint->e());
    //  cout << "Created track " << this << endl;
  }
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
  fP = aTrack.fP;
  fPt = aTrack.fPt;
  fInnerMomentum = aTrack.fInnerMomentum;
  fHelix = aTrack.fHelix;
  fTrackId = aTrack.fTrackId;
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
  fNSigmaTOFPi=aTrack.fNSigmaTOFPi;    
  fNSigmaTOFK=aTrack.fNSigmaTOFK;    
  fNSigmaTOFP=aTrack.fNSigmaTOFP;    
  fClusters=aTrack.fClusters;
  fShared=aTrack.fShared;
  fNominalTpcEntrancePoint=aTrack.fNominalTpcEntrancePoint;
  fNominalTpcExitPoint=aTrack.fNominalTpcExitPoint;
  fKinkIndexes[0] = aTrack.fKinkIndexes[0];
  fKinkIndexes[1] = aTrack.fKinkIndexes[1];
  fKinkIndexes[2] = aTrack.fKinkIndexes[2];
  //fTrueMomentum = aTrack.fTrueMomentum;
  //fEmissionPoint(0), // Emission point coordinates
  //fPDGPid(0),        // True PID of the particle
  //  fGlobalEmissionPoint(0)
  fMass = aTrack.fMass;
  fPDGPid = aTrack.fPDGPid;

  for(int i=0;i<6;i++) {
    fHasPointOnITS[i]=false;
  }

  fXatDCA=aTrack.fXatDCA;
  fYatDCA=aTrack.fYatDCA;
  fZatDCA=aTrack.fZatDCA;

  for(int i=0;i<9;i++)
    fNominalTpcPoints[i] = aTrack.fNominalTpcPoints[i];

  if (ValidHiddenInfo())
    delete fHiddenInfo;
  if (aTrack.ValidHiddenInfo())
    fHiddenInfo = aTrack.GetHiddenInfo()->Clone();
  else 
    fHiddenInfo = NULL;


  if(!fTrueMomentum && aTrack.fTrueMomentum)
    fTrueMomentum = new AliFemtoThreeVector();
  if(aTrack.fTrueMomentum){
    fTrueMomentum->SetX(aTrack.fTrueMomentum->x());
    fTrueMomentum->SetY(aTrack.fTrueMomentum->y());
    fTrueMomentum->SetZ(aTrack.fTrueMomentum->z());}

  if(!fEmissionPoint && aTrack.fEmissionPoint) 
    fEmissionPoint = new AliFemtoLorentzVector();
  if(aTrack.fEmissionPoint){
    fEmissionPoint->SetX(aTrack.fEmissionPoint->x());
    fEmissionPoint->SetY(aTrack.fEmissionPoint->y());
    fEmissionPoint->SetZ(aTrack.fEmissionPoint->z());
    fEmissionPoint->SetT(aTrack.fEmissionPoint->e());
  }
 
  if(!fGlobalEmissionPoint && aTrack.fGlobalEmissionPoint)
    fGlobalEmissionPoint = new AliFemtoThreeVector();
  if(aTrack.fGlobalEmissionPoint){
    fGlobalEmissionPoint->SetX(aTrack.fGlobalEmissionPoint->x());
    fGlobalEmissionPoint->SetY(aTrack.fGlobalEmissionPoint->y());
    fGlobalEmissionPoint->SetZ(aTrack.fGlobalEmissionPoint->z());
    //fGlobalEmissionPoint->SetT(t.fGlobalEmissionPoint->e());
    //  cout << "Created track " << this << endl;
  }

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
void AliFemtoTrack::SetNSigmaTOFPi(const float& aNSigmaTOFPi){fNSigmaTOFPi=aNSigmaTOFPi;}
void AliFemtoTrack::SetNSigmaTOFK(const float& aNSigmaTOFK){fNSigmaTOFK=aNSigmaTOFK;}
void AliFemtoTrack::SetNSigmaTOFP(const float& aNSigmaTOFP){fNSigmaTOFP=aNSigmaTOFP;}
void AliFemtoTrack::SetSigmaToVertex(const float& aSigma){fSigmaToVertex=aSigma;} 

void AliFemtoTrack::SetXatDCA(const double& x) {fXatDCA=x;}
void AliFemtoTrack::SetYatDCA(const double& x) {fYatDCA=x;}
void AliFemtoTrack::SetZatDCA(const double& x) {fZatDCA=x;}


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
float AliFemtoTrack::NSigmaTOFPi() const{return fNSigmaTOFPi;}
float AliFemtoTrack::NSigmaTOFK() const{return fNSigmaTOFK;}
float AliFemtoTrack::NSigmaTOFP() const{return fNSigmaTOFP;}
float AliFemtoTrack::SigmaToVertex() const{return fSigmaToVertex;} 
float AliFemtoTrack::TOFpionTime() const{return fTofPionTime;}
float AliFemtoTrack::TOFkaonTime() const{return fTofKaonTime;}
float AliFemtoTrack::TOFprotonTime() const{return fTofProtonTime;}

double AliFemtoTrack::XatDCA() const {return fXatDCA;}
double AliFemtoTrack::YatDCA() const {return fYatDCA;}
double AliFemtoTrack::ZatDCA() const {return fZatDCA;}

void AliFemtoTrack::SetHiddenInfo(AliFemtoHiddenInfo* aHiddenInfo) {fHiddenInfo=aHiddenInfo;}
bool AliFemtoTrack::ValidHiddenInfo() const { if (fHiddenInfo) return true; else return false; }
AliFemtoHiddenInfo* AliFemtoTrack::GetHiddenInfo() const {return fHiddenInfo;}
  
AliFemtoTrack::~AliFemtoTrack()
{
  // destructor
  if (fHiddenInfo)
    delete fHiddenInfo;

  if(fTrueMomentum) delete fTrueMomentum;
  if(fEmissionPoint) delete fEmissionPoint;
  if(fGlobalEmissionPoint) delete fGlobalEmissionPoint;


  //  cout << "Deleted track " << this << endl;
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
  if ((aIndex <3) && (aIndex>=0))
    return fKinkIndexes[aIndex];
  else
    return 0;
}

bool AliFemtoTrack::HasPointOnITSLayer(int aIndex) const
{
  // Return if i-th ITS layer had a hit for this track
  if ((aIndex <6) && (aIndex>=0))
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
    fTrueMomentum->SetX(aMom->x());
    fTrueMomentum->SetY(aMom->y());
    fTrueMomentum->SetZ(aMom->z());
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
    fTrueMomentum->SetX(aMom.x());
    fTrueMomentum->SetY(aMom.y());
    fTrueMomentum->SetZ(aMom.z());
  }
  else {
    fTrueMomentum = new AliFemtoThreeVector();
    *fTrueMomentum = aMom;
  }
}
//_____________________________________________
void                   AliFemtoTrack::SetTrueMomentum(Double_t aPx, Double_t aPy, Double_t aPz)
{
  // Set momentum from components
  if (!fTrueMomentum) fTrueMomentum = new AliFemtoThreeVector();
    fTrueMomentum->SetX(aPx);
    fTrueMomentum->SetY(aPy);
    fTrueMomentum->SetZ(aPz);
}
//_____________________________________________
void                   AliFemtoTrack::SetEmissionPoint(AliFemtoLorentzVector *aPos)
{
  // Set position from vector
  if (fEmissionPoint) {
    fEmissionPoint->SetX(aPos->px());
    fEmissionPoint->SetY(aPos->py());
    fEmissionPoint->SetZ(aPos->pz());
    fEmissionPoint->SetT(aPos->e());
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
    fEmissionPoint->SetX(aPos.px());
    fEmissionPoint->SetY(aPos.py());
    fEmissionPoint->SetZ(aPos.pz());
    fEmissionPoint->SetT(aPos.e());
  }
  else {
    fEmissionPoint = new AliFemtoLorentzVector();
    *fEmissionPoint = aPos;
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
    fGlobalEmissionPoint->SetX(aPos.x());
    fGlobalEmissionPoint->SetY(aPos.y());
    fGlobalEmissionPoint->SetZ(aPos.z());
  }
  else {
    fGlobalEmissionPoint = new AliFemtoThreeVector();
    *fGlobalEmissionPoint = aPos;
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
