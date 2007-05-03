/***************************************************************************
 *
 * $Id$
 *
 * Author: Frank Laue, Ohio State, laue@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 * Implementation of methods
 *
 ***************************************************************************
 * $Log$
 * Revision 1.3  2007/04/27 07:24:34  akisiel
 * Make revisions needed for compilation from the main AliRoot tree
 *
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
 * Revision 1.3  2007-04-03 16:00:09  mchojnacki
 * Changes to iprove memory managing
 *
 * Revision 1.2  2007/03/08 14:58:04  mchojnacki
 * adding some alice stuff
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.19  2005/07/15 17:41:43  kopytin
 * initialized fHiddenInfo in AliFemtoTrack() to null
 *
 * Revision 1.18  2005/07/10 02:17:21  chajecki
 * Ftpc (Eeat+West) hits included in nHits()
 *
 * Revision 1.17  2005/06/17 21:41:12  chajecki
 * Two bugs fixed:
 * 1. wrong value of mNHits set in one of the constructors
 * AliFemtoTrack::AliFemtoTrack(const StTrack* ST, AliFemtoThreeVector PrimaryVertex)
 *
 * 2. since year 4 the StMuTrack::nHits() returns the total number of hits.
 * We're interested in TPC hits only so this value can be taken by
 *  StMuTrack::topologyMap().numberOfHits(kTpcId);
 * This change is backwards compatible so the code should work also
 * for data <=Y3
 *
 * Revision 1.16  2003/03/18 14:41:48  kisiel
 * Bugfix update for the theoretical part of the code. Reverts the changes to the Lednicky weight calculator, as the previuos one had problems with strong interaction
 *
 * Revision 1.15  2003/01/31 19:57:15  magestro
 * Cleared up simple compiler warnings on i386_linux24
 *
 * Revision 1.14  2003/01/21 17:26:33  magestro
 * Added condition to globalTrack() call so GlobalTracks branch in MuDst can be disabled
 *
 * Revision 1.13  2002/03/21 18:49:31  laue
 * updated for new MuDst reader
 *
 * Revision 1.12  2002/02/04 18:58:33  laue
 * *** empty log message ***
 *
 * Revision 1.11  2001/12/14 23:11:30  fretiere
 * Add class HitMergingCut. Add class fabricesPairCut = HitMerginCut + pair purity cuts. Add TpcLocalTransform function which convert to local tpc coord (not pretty). Modify AliFemtoTrack, AliFemtoParticle, AliFemtoHiddenInfo, AliFemtoPair to handle the hit information and cope with my code
 *
 * Revision 1.4  2001/06/21 19:15:48  laue
 * Modified fiels:
 *   CTH.h : new constructor added
 *   AliFemtoEvent, AliFemtoKink, AliFemtoTrack : constructors from the persistent
 *                                   (TTree) classes added
 *   AliFemtoLikeSignAnalysis : minor changes, for debugging
 *   AliFemtoTypes: split into different files
 * Added files: for the new TTree muDst's
 *   StExceptions.cxx StExceptions.h AliFemtoEnumeration.h
 *   AliFemtoHelix.h AliFemtoHisto.h AliFemtoString.h AliFemtoTFile.h
 *   AliFemtoTTreeEvent.cxx AliFemtoTTreeEvent.h AliFemtoTTreeKink.cxx
 *   AliFemtoTTreeKink.h AliFemtoTTreeTrack.cxx AliFemtoTTreeTrack.h
 *   AliFemtoTTreeV0.cxx AliFemtoTTreeV0.h AliFemtoVector.h
 *
 * Revision 1.3  2001/05/25 23:23:59  lisa
 * Added in AliFemtoKink stuff
 *
 *
 ****************************************************************************/

#include "Infrastructure/AliFemtoTrack.h" 
#ifndef __NO_STAR_DEPENDENCE_ALLOWED__
#ifdef __ROOT__
#include "StEvent/StEnumerations.h"
#include "Infrastructure/AliFemtoAihongPid.h"
#include "StEventUtilities/StuProbabilityPidAlgorithm.h"
#endif
#endif

#include "Base/SystemOfUnits.h"   // has "tesla" in it
//#include "AliFemtoParticleTypes.h"
//#include "Infrastructure/AliFemtoTTreeEvent.h" 
//#include "Infrastructure/AliFemtoTTreeTrack.h" 

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
  fHiddenInfo = NULL;
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
 { // copy constructor
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
    fHiddenInfo = t.getHiddenInfo()->clone();
  else 
    fHiddenInfo = NULL;
  //  cout << "Created track " << this << endl;
};

AliFemtoTrack& AliFemtoTrack::operator=(const AliFemtoTrack& aTrack)
{
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
  if (ValidHiddenInfo())
    delete fHiddenInfo;
  if (aTrack.ValidHiddenInfo())
    fHiddenInfo = aTrack.getHiddenInfo()->clone();
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
AliFemtoHiddenInfo* AliFemtoTrack::getHiddenInfo() const {return fHiddenInfo;}
  
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

