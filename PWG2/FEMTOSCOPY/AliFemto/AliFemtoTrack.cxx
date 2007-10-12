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
  fP(0,0,0),
  fPt(0),
  fHelix(),
  fFlags(0),
  fLabel(0),
  fImpactD(0),
  fImpactZ(0),
  fCdd(0),
  fCdz(0),
  fCzz(0),
  fITSchi2(0),       
  fITSncls(0),        
  fTPCchi2(0),       
  fTPCncls(0),       
  fTPCnclsF(0),      
  fTPCsignalN(0),    
  fTPCsignalS(0),
  fClusters(159),
  fShared(159),
  fHiddenInfo(0)
{
  // Default constructor
  fHiddenInfo = NULL;
  fKinkIndexes[0] = 0;
  fKinkIndexes[1] = 0;
  fKinkIndexes[2] = 0;
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
  fP(0,0,0),
  fPt(0),
  fHelix(),
  fFlags(0),
  fLabel(0),
  fImpactD(0),
  fImpactZ(0),
  fCdd(0),
  fCdz(0),
  fCzz(0),
  fITSchi2(0),       
  fITSncls(0),        
  fTPCchi2(0),       
  fTPCncls(0),       
  fTPCnclsF(0),      
  fTPCsignalN(0),    
  fTPCsignalS(0),
  fClusters(159),
  fShared(159),
  fHiddenInfo(0)
 { 
   // copy constructor
  fCharge = t.fCharge;
  fPidProbElectron = t.fPidProbElectron;
  fPidProbPion = t.fPidProbPion;
  fPidProbKaon = t.fPidProbKaon;
  fPidProbProton = t.fPidProbProton;
  fPidProbMuon=t.fPidProbMuon;
  fP = t.fP;
  fPt = t.fPt;
  fHelix = t.fHelix;
  fTrackId = t.fTrackId;
  fFlags=t.fFlags;
  fLabel=t.fLabel;
  fImpactD=t.fImpactD;
  fImpactZ=t.fImpactZ;
  fCdd=t.fCdd;
  fCdz=t.fCdz;
  fCzz=t.fCzz;
  fITSchi2=t.fITSchi2;       
  fITSncls=t.fITSncls;        
  fTPCchi2=t.fTPCchi2;       
  fTPCncls=t.fTPCncls;       
  fTPCnclsF=t.fTPCnclsF;      
  fTPCsignalN=t.fTPCsignalN;    
  fTPCsignalS=t.fTPCsignalS;  
  fClusters=t.fClusters;
  fShared=t.fShared;
  if (t.ValidHiddenInfo())
    fHiddenInfo = t.GetHiddenInfo()->Clone();
  else 
    fHiddenInfo = NULL;
  fKinkIndexes[0] = t.fKinkIndexes[0];
  fKinkIndexes[1] = t.fKinkIndexes[1];
  fKinkIndexes[2] = t.fKinkIndexes[2];
  //  cout << "Created track " << this << endl;
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
  fP = aTrack.fP;
  fPt = aTrack.fPt;
  fHelix = aTrack.fHelix;
  fTrackId = aTrack.fTrackId;
  fFlags=aTrack.fFlags;
  fLabel=aTrack.fLabel;
  fImpactD=aTrack.fImpactD;
  fImpactZ=aTrack.fImpactZ;
  fCdd=aTrack.fCdd;
  fCdz=aTrack.fCdz;
  fCzz=aTrack.fCzz;
  fITSchi2=aTrack.fITSchi2;       
  fITSncls=aTrack.fITSncls;        
  fTPCchi2=aTrack.fTPCchi2;       
  fTPCncls=aTrack.fTPCncls;       
  fTPCnclsF=aTrack.fTPCnclsF;      
  fTPCsignalN=aTrack.fTPCsignalN;    
  fTPCsignalS=aTrack.fTPCsignalS;  
  fClusters=aTrack.fClusters;
  fShared=aTrack.fShared;
  fKinkIndexes[0] = aTrack.fKinkIndexes[0];
  fKinkIndexes[1] = aTrack.fKinkIndexes[1];
  fKinkIndexes[2] = aTrack.fKinkIndexes[2];
  if (ValidHiddenInfo())
    delete fHiddenInfo;
  if (aTrack.ValidHiddenInfo())
    fHiddenInfo = aTrack.GetHiddenInfo()->Clone();
  else 
    fHiddenInfo = NULL;

  return *this;
}

void AliFemtoTrack::SetCharge(const short& ch){fCharge=ch;}

void AliFemtoTrack::SetPidProbElectron(const float& x){fPidProbElectron = x;}
void AliFemtoTrack::SetPidProbPion(const float& x){fPidProbPion = x;}
void AliFemtoTrack::SetPidProbKaon(const float& x){fPidProbKaon = x;}
void AliFemtoTrack::SetPidProbProton(const float& x){fPidProbProton = x;}
void AliFemtoTrack::SetPidProbMuon(const float& x){fPidProbMuon = x;}
 
void AliFemtoTrack::SetP(const AliFemtoThreeVector& p){fP = p;}
void AliFemtoTrack::SetPt(const float& pt){fPt = pt;} 
void AliFemtoTrack::SetHelix(const AliFmPhysicalHelixD& h){fHelix = h;}
void AliFemtoTrack::SetTrackId(const short & id) { fTrackId=id;}
void AliFemtoTrack::SetFlags(const long int &flags) {fFlags=flags;}
void AliFemtoTrack::SetLabel(const int &label) {fLabel=label;}
void AliFemtoTrack::SetImpactD(const float& ImpactD){fImpactD=ImpactD;}
void AliFemtoTrack::SetImpactZ(const float& ImpactZ){fImpactZ=ImpactZ;}
void AliFemtoTrack::SetCdd(const float& Cdd){fCdd=Cdd;}
void AliFemtoTrack::SetCdz(const float& Cdz){fCdz=Cdz;}
void AliFemtoTrack::SetCzz(const float& Czz){fCzz=Czz;}
void AliFemtoTrack::SetITSchi2(const float& ITSchi2){fITSchi2=ITSchi2;}    
void AliFemtoTrack::SetITSncls(const int& ITSncls){fITSncls=ITSncls;}     
void AliFemtoTrack::SetTPCchi2(const float& TPCchi2){fTPCchi2=TPCchi2;}       
void AliFemtoTrack::SetTPCncls(const int& TPCncls){fTPCncls=TPCncls;}       
void AliFemtoTrack::SetTPCnclsF(const short& TPCnclsF){fTPCnclsF=TPCnclsF;}      
void AliFemtoTrack::SetTPCsignalN(const short& TPCsignalN){fTPCsignalN=TPCsignalN;}    
void AliFemtoTrack::SetTPCsignalS(const float& TPCsignalS){fTPCsignalS=TPCsignalS;} 


short AliFemtoTrack::Charge() const {return fCharge;}  
AliFemtoThreeVector AliFemtoTrack::P() const {return fP;}
float AliFemtoTrack::Pt() const {return fPt;}              
const AliFmPhysicalHelixD& AliFemtoTrack::Helix() const {return fHelix;}
short AliFemtoTrack::TrackId() const { return fTrackId; }
long int AliFemtoTrack::Flags() const {return fFlags;}
int AliFemtoTrack::Label()const {return fLabel;}
float AliFemtoTrack::ImpactD()const{return fImpactD;}
float AliFemtoTrack::ImpactZ()const{return fImpactZ;}
float AliFemtoTrack::Cdd() const{return fCdd;}
float AliFemtoTrack::Cdz() const{return fCdz;}
float AliFemtoTrack::Czz() const{return fCzz;} 
float AliFemtoTrack::ITSchi2() const{return fITSchi2;}    
int   AliFemtoTrack::ITSncls() const{return fITSncls;}     
float AliFemtoTrack::TPCchi2() const{return fTPCchi2;}       
int   AliFemtoTrack::TPCncls() const{return fTPCncls;}       
short AliFemtoTrack::TPCnclsF() const{return fTPCnclsF;}      
short AliFemtoTrack::TPCsignalN() const{return fTPCsignalN;}    
float AliFemtoTrack::TPCsignalS() const{return fTPCsignalS;} 

void AliFemtoTrack::SetHiddenInfo(AliFemtoHiddenInfo* aHiddenInfo) {fHiddenInfo=aHiddenInfo;}
bool AliFemtoTrack::ValidHiddenInfo() const { if (fHiddenInfo) return true; else return false; }
AliFemtoHiddenInfo* AliFemtoTrack::GetHiddenInfo() const {return fHiddenInfo;}
  
AliFemtoTrack::~AliFemtoTrack()
{
  if (fHiddenInfo)
    delete fHiddenInfo;
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
  fKinkIndexes[0] = points[0];
  fKinkIndexes[1] = points[1];
  fKinkIndexes[2] = points[2];
}

int  AliFemtoTrack::KinkIndex(int aIndex) const
{
  if ((aIndex <3) && (aIndex>=0))
    return fKinkIndexes[aIndex];
  else
    return 0;
}
