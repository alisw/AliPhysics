////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoEventReaderESDChain - the reader class for the Alice ESD         ///
/// tailored for the Task framework                                 ///
/// Reads in AliESDfriend to create shared hit/quality information           ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoEventReaderESDChain.h"

#include "TFile.h"
#include "TTree.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliMultiplicity.h"
#include "AliCentrality.h"
#include "AliESDVZERO.h"

#include "AliFmPhysicalHelixD.h"
#include "AliFmThreeVectorF.h"

#include "SystemOfUnits.h"

#include "AliFemtoEvent.h"
#include "AliFemtoModelHiddenInfo.h"

ClassImp(AliFemtoEventReaderESDChain)

#if !(ST_NO_NAMESPACES)
  using namespace units;
#endif

using namespace std;
//____________________________
AliFemtoEventReaderESDChain::AliFemtoEventReaderESDChain():
  fFileName(" "),
  fConstrained(true),
  fReadInner(false),
  fUseTPCOnly(false),
  fNumberofEvent(0),
  fCurEvent(0),
  fCurFile(0),
  fEvent(0x0),
  fUsePhysicsSel(kFALSE),
  fSelect(0x0),
  fTrackType(kGlobal),
  fEstEventMult(kITSTPC)
{
  //constructor with 0 parameters , look at default settings 
//   fClusterPerPadrow = (list<Int_t> **) malloc(sizeof(list<Int_t> *) * AliESDfriendTrack::kMaxTPCcluster);
//   for (int tPad=0; tPad<AliESDfriendTrack::kMaxTPCcluster; tPad++) {
//     fClusterPerPadrow[tPad] = new list<Int_t>();
//   }
//   fSharedList = (list<Int_t> **) malloc(sizeof(list<Int_t> *) * AliESDfriendTrack::kMaxTPCcluster);
//   for (int tPad=0; tPad<AliESDfriendTrack::kMaxTPCcluster; tPad++) {
//     fSharedList[tPad] = new list<Int_t>();
//   }
}

//__________________
AliFemtoEventReaderESDChain::AliFemtoEventReaderESDChain(const AliFemtoEventReaderESDChain& aReader):
  AliFemtoEventReader(aReader),
  fFileName(" "),
  fConstrained(true),
  fReadInner(false),
  fUseTPCOnly(false),
  fNumberofEvent(0),
  fCurEvent(0),
  fCurFile(0),
  fEvent(0x0),
  fUsePhysicsSel(kFALSE),
  fSelect(0x0),
  fTrackType(kGlobal),
  fEstEventMult(kITSTPC)
{
  // Copy constructor
  fConstrained = aReader.fConstrained;
  fReadInner = aReader.fReadInner;
  fUseTPCOnly = aReader.fUseTPCOnly;
  fNumberofEvent = aReader.fNumberofEvent;
  fCurEvent = aReader.fCurEvent;
  fCurFile = aReader.fCurFile;
  //  fEvent = new AliESD(*aReader.fEvent);
  fEvent = new AliESDEvent();
  fUsePhysicsSel = aReader.fUsePhysicsSel;
  if (aReader.fUsePhysicsSel)
    fSelect = new AliPhysicsSelection();
  fTrackType = aReader.fTrackType;
  fEstEventMult = aReader.fEstEventMult;
//   fEventFriend = aReader.fEventFriend;
//   fClusterPerPadrow = (list<Int_t> **) malloc(sizeof(list<Int_t> *) * AliESDfriendTrack::kMaxTPCcluster);
//   for (int tPad=0; tPad<AliESDfriendTrack::kMaxTPCcluster; tPad++) {
//     fClusterPerPadrow[tPad] = new list<Int_t>();
//     list<Int_t>::iterator iter;
//     for (iter=aReader.fClusterPerPadrow[tPad]->begin(); iter!=aReader.fClusterPerPadrow[tPad]->end(); iter++) {
//       fClusterPerPadrow[tPad]->push_back(*iter);
//     }
//   }
//   fSharedList = (list<Int_t> **) malloc(sizeof(list<Int_t> *) * AliESDfriendTrack::kMaxTPCcluster);
//   for (int tPad=0; tPad<AliESDfriendTrack::kMaxTPCcluster; tPad++) {
//     fSharedList[tPad] = new list<Int_t>();
//     list<Int_t>::iterator iter;
//     for (iter=aReader.fSharedList[tPad]->begin(); iter!=aReader.fSharedList[tPad]->end(); iter++) {
//       fSharedList[tPad]->push_back(*iter);
//     }
//   }
}
//__________________
AliFemtoEventReaderESDChain::~AliFemtoEventReaderESDChain()
{
  //Destructor
  delete fEvent;

//   for (int tPad=0; tPad<AliESDfriendTrack::kMaxTPCcluster; tPad++) {
//     fClusterPerPadrow[tPad]->clear();
//     delete fClusterPerPadrow[tPad];
//   }
//   delete [] fClusterPerPadrow;
//   for (int tPad=0; tPad<AliESDfriendTrack::kMaxTPCcluster; tPad++) {
//     fSharedList[tPad]->clear();
//     delete fSharedList[tPad];
//   }
//   delete [] fSharedList;
  if (fSelect) delete fSelect;
}

//__________________
AliFemtoEventReaderESDChain& AliFemtoEventReaderESDChain::operator=(const AliFemtoEventReaderESDChain& aReader)
{
  // Assignment operator
  if (this == &aReader)
    return *this;

  fConstrained = aReader.fConstrained;
  fReadInner = aReader.fReadInner;
  fUseTPCOnly = aReader.fUseTPCOnly;
  fNumberofEvent = aReader.fNumberofEvent;
  fCurEvent = aReader.fCurEvent;
  fCurFile = aReader.fCurFile;
  if (fEvent) delete fEvent;
  fEvent = new AliESDEvent();
  fTrackType = aReader.fTrackType;
  fEstEventMult = aReader.fEstEventMult;

  fUsePhysicsSel = aReader.fUsePhysicsSel;
  if (aReader.fUsePhysicsSel)
    fSelect = new AliPhysicsSelection();
  //  fEventFriend = aReader.fEventFriend;
  
//   if (fClusterPerPadrow) {
//     for (int tPad=0; tPad<AliESDfriendTrack::kMaxTPCcluster; tPad++) {
//       fClusterPerPadrow[tPad]->clear();
//       delete fClusterPerPadrow[tPad];
//     }
//     delete [] fClusterPerPadrow;
//   }
  
//   if (fSharedList) {
//     for (int tPad=0; tPad<AliESDfriendTrack::kMaxTPCcluster; tPad++) {
//       fSharedList[tPad]->clear();
//       delete fSharedList[tPad];
//     }
//     delete [] fSharedList;
//   }

//   fClusterPerPadrow = (list<Int_t> **) malloc(sizeof(list<Int_t> *) * AliESDfriendTrack::kMaxTPCcluster);
//   for (int tPad=0; tPad<AliESDfriendTrack::kMaxTPCcluster; tPad++) {
//     fClusterPerPadrow[tPad] = new list<Int_t>();
//     list<Int_t>::iterator iter;
//     for (iter=aReader.fClusterPerPadrow[tPad]->begin(); iter!=aReader.fClusterPerPadrow[tPad]->end(); iter++) {
//       fClusterPerPadrow[tPad]->push_back(*iter);
//     }
//   }
//   fSharedList = (list<Int_t> **) malloc(sizeof(list<Int_t> *) * AliESDfriendTrack::kMaxTPCcluster);
//   for (int tPad=0; tPad<AliESDfriendTrack::kMaxTPCcluster; tPad++) {
//     fSharedList[tPad] = new list<Int_t>();
//     list<Int_t>::iterator iter;
//     for (iter=aReader.fSharedList[tPad]->begin(); iter!=aReader.fSharedList[tPad]->end(); iter++) {
//       fSharedList[tPad]->push_back(*iter);
//     }
//   }
  
  return *this;
}
//__________________
// Simple report
AliFemtoString AliFemtoEventReaderESDChain::Report()
{
  AliFemtoString temp = "\n This is the AliFemtoEventReaderESDChain\n";
  return temp;
}

//__________________
void AliFemtoEventReaderESDChain::SetConstrained(const bool constrained)
{
  // Select whether to read constrained or not constrained momentum
  fConstrained=constrained;
}
//__________________
bool AliFemtoEventReaderESDChain::GetConstrained() const
{
  // Check whether we read constrained or not constrained momentum
  return fConstrained;
}
//__________________
void AliFemtoEventReaderESDChain::SetReadTPCInner(const bool readinner)
{
  fReadInner=readinner;
}

bool AliFemtoEventReaderESDChain::GetReadTPCInner() const
{
  return fReadInner;
}

//__________________
void AliFemtoEventReaderESDChain::SetUseTPCOnly(const bool usetpconly)
{
  fUseTPCOnly=usetpconly;
}

bool AliFemtoEventReaderESDChain::GetUseTPCOnly() const
{
  return fUseTPCOnly;
}

void AliFemtoEventReaderESDChain::SetUsePhysicsSelection(const bool usephysics)
{
  fUsePhysicsSel = usephysics;
  if (!fSelect) fSelect = new AliPhysicsSelection();
}

void AliFemtoEventReaderESDChain::SetUseMultiplicity(EstEventMult aType)
{
  fEstEventMult = aType;
}

AliFemtoEvent* AliFemtoEventReaderESDChain::ReturnHbtEvent()
{
  // Get the event, read all the relevant information
  // and fill the AliFemtoEvent class
  // Returns a valid AliFemtoEvent
  AliFemtoEvent *hbtEvent = 0;
  string tFriendFileName;

  // Get the friend information
  cout<<"starting to read event "<<fCurEvent<<endl;
  //  fEvent->SetESDfriend(fEventFriend);
  if(fEvent->GetAliESDOld())fEvent->CopyFromOldESD();
  
  hbtEvent = new AliFemtoEvent;

  if (fUsePhysicsSel) {
    hbtEvent->SetIsCollisionCandidate(fSelect->IsCollisionCandidate(fEvent));
    if (!(fSelect->IsCollisionCandidate(fEvent)))
      printf("Event not a collision candidate\n");
  }
  else
    hbtEvent->SetIsCollisionCandidate(kTRUE);

  //setting basic things
  //  hbtEvent->SetEventNumber(fEvent->GetEventNumber());
  hbtEvent->SetRunNumber(fEvent->GetRunNumber());
  //hbtEvent->SetNumberOfTracks(fEvent->GetNumberOfTracks());
  hbtEvent->SetMagneticField(fEvent->GetMagneticField()*kilogauss);//to check if here is ok
  hbtEvent->SetZDCN1Energy(fEvent->GetZDCN1Energy());
  hbtEvent->SetZDCP1Energy(fEvent->GetZDCP1Energy());
  hbtEvent->SetZDCN2Energy(fEvent->GetZDCN2Energy());
  hbtEvent->SetZDCP2Energy(fEvent->GetZDCP2Energy());
  hbtEvent->SetZDCEMEnergy(fEvent->GetZDCEMEnergy());
  hbtEvent->SetZDCParticipants(fEvent->GetZDCParticipants());
  hbtEvent->SetTriggerMask(fEvent->GetTriggerMask());
  //  hbtEvent->SetTriggerCluster(fEvent->GetTriggerCluster());

  if ((fEvent->IsTriggerClassFired("CINT1WU-B-NOPF-ALL")) ||
      (fEvent->IsTriggerClassFired("CINT1B-ABCE-NOPF-ALL")))
    hbtEvent->SetTriggerCluster(1);
  else if (fEvent->IsTriggerClassFired("CSH1WU-B-NOPF-ALL"))
    hbtEvent->SetTriggerCluster(2);
  else 
    hbtEvent->SetTriggerCluster(0);
	
  //Vertex
  double fV1[3];
  double fVCov[6];
  if (fUseTPCOnly) {
    fEvent->GetPrimaryVertexTPC()->GetXYZ(fV1);
    fEvent->GetPrimaryVertexTPC()->GetCovMatrix(fVCov);
    if (!fEvent->GetPrimaryVertexTPC()->GetStatus())
      fVCov[4] = -1001.0;
  }
  else {
    fEvent->GetPrimaryVertex()->GetXYZ(fV1);
    fEvent->GetPrimaryVertex()->GetCovMatrix(fVCov);
    if (!fEvent->GetPrimaryVertex()->GetStatus())
      fVCov[4] = -1001.0;
  }

  AliFmThreeVectorF vertex(fV1[0],fV1[1],fV1[2]);
  hbtEvent->SetPrimVertPos(vertex);
  hbtEvent->SetPrimVertCov(fVCov);
  
  Int_t spdetaonecount = 0;
  
  for (int iter=0; iter<fEvent->GetMultiplicity()->GetNumberOfTracklets(); iter++) 
    if (fabs(fEvent->GetMultiplicity()->GetEta(iter)) < 1.0)
      spdetaonecount++;

  //  hbtEvent->SetSPDMult(fEvent->GetMultiplicity()->GetNumberOfTracklets());
  hbtEvent->SetSPDMult(spdetaonecount);

  //starting to reading tracks
  int nofTracks=0;  //number of reconstructed tracks in event
  nofTracks=fEvent->GetNumberOfTracks();
  int realnofTracks=0;//number of track which we use ina analysis

//   // Clear the shared cluster list
//   for (int tPad=0; tPad<AliESDfriendTrack::kMaxTPCcluster; tPad++) {
//     fClusterPerPadrow[tPad]->clear();
//   }
//   for (int tPad=0; tPad<AliESDfriendTrack::kMaxTPCcluster; tPad++) {
//     fSharedList[tPad]->clear();
//   }


//   for (int i=0;i<nofTracks;i++) {
//     const AliESDtrack *esdtrack=fEvent->GetTrack(i);//getting next track

//     list<Int_t>::iterator tClustIter;

//     Int_t tTrackIndices[AliESDfriendTrack::kMaxTPCcluster];
//     Int_t tNClusters = esdtrack->GetTPCclusters(tTrackIndices);
//     for (int tNcl=0; tNcl<AliESDfriendTrack::kMaxTPCcluster; tNcl++) {
//       if (tTrackIndices[tNcl] >= 0) {
// 	tClustIter = find(fClusterPerPadrow[tNcl]->begin(), fClusterPerPadrow[tNcl]->end(), tTrackIndices[tNcl]);
// 	if (tClustIter == fClusterPerPadrow[tNcl]->end()) {
// 	  fClusterPerPadrow[tNcl]->push_back(tTrackIndices[tNcl]);
// 	}
// 	else {
// 	  fSharedList[tNcl]->push_back(tTrackIndices[tNcl]);
// 	}
//       }
//     }
      
//   }

  int tNormMult = 0;
  int tNormMultPos = 0;
  int tNormMultNeg = 0;

  Float_t tTotalPt = 0.0;

  Float_t b[2];
  Float_t bCov[3];

  Int_t tTracklet=0, tITSTPC=0, tITSPure=0;
  
  fEvent->EstimateMultiplicity(tTracklet, tITSTPC, tITSPure, 1.2);
  
  hbtEvent->SetMultiplicityEstimateITSTPC(tITSTPC);
  hbtEvent->SetMultiplicityEstimateTracklets(tTracklet);
  //  hbtEvent->SetMultiplicityEstimateITSPure(tITSPure);
  hbtEvent->SetMultiplicityEstimateITSPure(fEvent->GetMultiplicity()->GetNumberOfITSClusters(1));
  
  for (int i=0;i<nofTracks;i++)
    {
      bool  tGoodMomentum=true; //flaga to chcek if we can read momentum of this track

      const AliESDtrack *esdtrack=fEvent->GetTrack(i);//getting next track
      //      const AliESDfriendTrack *tESDfriendTrack = esdtrack->GetFriendTrack();

      if ((esdtrack->GetStatus() & AliESDtrack::kTPCrefit) &&
	  (esdtrack->GetStatus() & AliESDtrack::kITSrefit)) {
	if (esdtrack->GetTPCNcls() > 70) 
	  if (esdtrack->GetTPCchi2()/esdtrack->GetTPCNcls() < 4.0) {
	    if (TMath::Abs(esdtrack->Eta()) < 1.2) {
	      esdtrack->GetImpactParameters(b,bCov);
	      if ((b[0]<0.2) && (b[1] < 0.25)) {
		tNormMult++;
		tTotalPt += esdtrack->Pt();
	      }
	    }
	  }
      }
      else if (esdtrack->GetStatus() & AliESDtrack::kTPCrefit) {
	if (esdtrack->GetTPCNcls() > 100) 
	  if (esdtrack->GetTPCchi2()/esdtrack->GetTPCNcls() < 4.0) {
	    if (TMath::Abs(esdtrack->Eta()) < 1.2) {
	      esdtrack->GetImpactParameters(b,bCov);
	      if ((b[0]<2.4) && (b[1] < 3.2)) {
		tNormMult++;
		tTotalPt += esdtrack->Pt();
	      }
	    }
	  }
      }
      
      hbtEvent->SetZDCEMEnergy(tTotalPt);
//       if (esdtrack->GetStatus() & AliESDtrack::kTPCrefit)
// 	if (esdtrack->GetTPCNcls() > 80) 
// 	  if (esdtrack->GetTPCchi2()/esdtrack->GetTPCNcls() < 6.0) 
// 	    if (esdtrack->GetConstrainedParam())
// 	      if (fabs(esdtrack->GetConstrainedParam()->Eta()) < 0.5)
// 		if (esdtrack->GetConstrainedParam()->Pt() < 1.0) {
// 		  if (esdtrack->GetSign() > 0)
// 		    tNormMultPos++;
// 		  else if (esdtrack->GetSign() < 0)
// 		    tNormMultNeg--;
// 		}

      // If reading ITS-only tracks, reject all with TPC
      if (fTrackType == kITSOnly) {
	if (esdtrack->GetStatus() & AliESDtrack::kTPCrefit) continue;
	if (!(esdtrack->GetStatus() & AliESDtrack::kITSrefit)) continue;
	if (esdtrack->GetStatus() & AliESDtrack::kTPCin) continue;
	UChar_t iclm = esdtrack->GetITSClusterMap();
	Int_t incls = 0;
	for (int iter=0; iter<6; iter++) if (iclm&(1<<iter)) incls++;
	if (incls<=3) {
	  cout << "Rejecting track with " << incls << " clusters" << endl;
	  continue;
	}
      }

      AliFemtoTrack* trackCopy = new AliFemtoTrack();	
      trackCopy->SetCharge((short)esdtrack->GetSign());

      //in aliroot we have AliPID 
      //0-electron 1-muon 2-pion 3-kaon 4-proton 5-photon 6-pi0 7-neutron 8-kaon0 9-eleCon   
      //we use only 5 first
      double esdpid[5];
//       esdtrack->GetESDpid(esdpid);
      esdtrack->GetTPCpid(esdpid);
      trackCopy->SetPidProbElectron(esdpid[0]);
      trackCopy->SetPidProbMuon(esdpid[1]);
      trackCopy->SetPidProbPion(esdpid[2]);
      trackCopy->SetPidProbKaon(esdpid[3]);
      trackCopy->SetPidProbProton(esdpid[4]);
						
      double pxyz[3];
      double rxyz[3];
      double impact[2];
      double covimpact[3];
      
      if (fUseTPCOnly) {
	if (!esdtrack->GetTPCInnerParam()) {
	  delete trackCopy;
	  continue;
	}


	AliExternalTrackParam *param = new AliExternalTrackParam(*esdtrack->GetTPCInnerParam());
	param->GetXYZ(rxyz);
	param->PropagateToDCA(fEvent->GetPrimaryVertexTPC(), (fEvent->GetMagneticField()), 10000, impact, covimpact);
	param->GetPxPyPz(pxyz);//reading noconstarined momentum

	if (fReadInner == true) {
	  AliFemtoModelHiddenInfo *tInfo = new AliFemtoModelHiddenInfo();
	  tInfo->SetPDGPid(211);
	  tInfo->SetTrueMomentum(pxyz[0], pxyz[1], pxyz[2]);
	  tInfo->SetMass(0.13957);
	  //	  tInfo->SetEmissionPoint(rxyz[0], rxyz[1], rxyz[2], 0.0);
	  //	  tInfo->SetEmissionPoint(fV1[0], fV1[1], fV1[2], 0.0);
	  tInfo->SetEmissionPoint(rxyz[0]-fV1[0], rxyz[1]-fV1[1], rxyz[2]-fV1[2], 0.0);
	  trackCopy->SetHiddenInfo(tInfo);
	}
	
	AliFemtoThreeVector v(pxyz[0],pxyz[1],pxyz[2]);
	if (v.Mag() < 0.0001) {
	  //	cout << "Found 0 momentum ???? " <<endl;
	  delete trackCopy;
	  continue;
	}
	trackCopy->SetP(v);//setting momentum
	trackCopy->SetPt(sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]));

	const AliFmThreeVectorD kP(pxyz[0],pxyz[1],pxyz[2]);
	const AliFmThreeVectorD kOrigin(fV1[0],fV1[1],fV1[2]);
	//setting helix I do not if it is ok
	AliFmPhysicalHelixD helix(kP,kOrigin,(double)(fEvent->GetMagneticField())*kilogauss,(double)(trackCopy->Charge())); 
	trackCopy->SetHelix(helix);

	//some stuff which could be useful 
	trackCopy->SetImpactD(impact[0]);
	trackCopy->SetImpactZ(impact[1]);
	trackCopy->SetCdd(covimpact[0]);
	trackCopy->SetCdz(covimpact[1]);
	trackCopy->SetCzz(covimpact[2]);
	trackCopy->SetSigmaToVertex(GetSigmaToVertex(impact, covimpact));	

	delete param;
      }
      else {
	if (fReadInner == true) {
	  
	  if (esdtrack->GetTPCInnerParam()) {
	    AliExternalTrackParam *param = new AliExternalTrackParam(*esdtrack->GetTPCInnerParam());
	    param->GetXYZ(rxyz);
	    //	    param->PropagateToDCA(fEvent->GetPrimaryVertex(), (fEvent->GetMagneticField()), 10000);
	    param->GetPxPyPz(pxyz);//reading noconstarined momentum
	    delete param;
	    
	    AliFemtoModelHiddenInfo *tInfo = new AliFemtoModelHiddenInfo();
	    tInfo->SetPDGPid(211);
	    tInfo->SetTrueMomentum(pxyz[0], pxyz[1], pxyz[2]);
	    tInfo->SetMass(0.13957);
	    //	    tInfo->SetEmissionPoint(rxyz[0], rxyz[1], rxyz[2], 0.0);
	    //tInfo->SetEmissionPoint(fV1[0], fV1[1], fV1[2], 0.0);
	    tInfo->SetEmissionPoint(rxyz[0]-fV1[0], rxyz[1]-fV1[1], rxyz[2]-fV1[2], 0.0);
	    trackCopy->SetHiddenInfo(tInfo);
	  }
	}

	if (fTrackType == kGlobal) {
	  if (fConstrained==true)		    
	    tGoodMomentum=esdtrack->GetConstrainedPxPyPz(pxyz); //reading constrained momentum
	  else
	    tGoodMomentum=esdtrack->GetPxPyPz(pxyz);//reading noconstarined momentum
	}
	else if (fTrackType == kTPCOnly) {
	  if (esdtrack->GetTPCInnerParam())
	    esdtrack->GetTPCInnerParam()->GetPxPyPz(pxyz);
	  else {
	    delete trackCopy;
	    continue;
	  }
	}
	else if (fTrackType == kITSOnly) {
	  if (fConstrained==true)		    
	    tGoodMomentum=esdtrack->GetConstrainedPxPyPz(pxyz); //reading constrained momentum
	  else
	    tGoodMomentum=esdtrack->GetPxPyPz(pxyz);//reading noconstarined momentum
	}


	AliFemtoThreeVector v(pxyz[0],pxyz[1],pxyz[2]);
	if (v.Mag() < 0.0001) {
	  //	cout << "Found 0 momentum ???? " <<endl;
	  delete trackCopy;
	  continue;
	}
	trackCopy->SetP(v);//setting momentum
	trackCopy->SetPt(sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]));
	const AliFmThreeVectorD kP(pxyz[0],pxyz[1],pxyz[2]);
	const AliFmThreeVectorD kOrigin(fV1[0],fV1[1],fV1[2]);
	//setting helix I do not if it is ok
	AliFmPhysicalHelixD helix(kP,kOrigin,(double)(fEvent->GetMagneticField())*kilogauss,(double)(trackCopy->Charge())); 
	trackCopy->SetHelix(helix);

	//some stuff which could be useful 
	float imp[2];
	float cim[3];
	// if (fTrackType == kTPCOnly) {
	//   esdtrack->GetTPCInnerParam()->GetImpactParameters(imp,cim);
	// }
	// else {
	  esdtrack->GetImpactParameters(imp,cim);
	// }

	impact[0] = imp[0];
	impact[1] = imp[1];
	covimpact[0] = cim[0];
	covimpact[1] = cim[1];
	covimpact[2] = cim[2];

	trackCopy->SetImpactD(impact[0]);
	trackCopy->SetImpactZ(impact[1]);
	trackCopy->SetCdd(covimpact[0]);
	trackCopy->SetCdz(covimpact[1]);
	trackCopy->SetCzz(covimpact[2]);
	trackCopy->SetSigmaToVertex(GetSigmaToVertex(impact,covimpact));
      }

      trackCopy->SetTrackId(esdtrack->GetID());
      trackCopy->SetFlags(esdtrack->GetStatus());
      trackCopy->SetLabel(esdtrack->GetLabel());
		
      trackCopy->SetITSchi2(esdtrack->GetITSchi2());    
      if (esdtrack->GetITSFakeFlag())
	trackCopy->SetITSncls(-esdtrack->GetNcls(0));     
      else
	trackCopy->SetITSncls(esdtrack->GetNcls(0));     
      trackCopy->SetTPCchi2(esdtrack->GetTPCchi2());       
      trackCopy->SetTPCncls(esdtrack->GetTPCNcls());       
      trackCopy->SetTPCnclsF(esdtrack->GetTPCNclsF());      
      trackCopy->SetTPCsignal(esdtrack->GetTPCsignal());
      trackCopy->SetTPCsignalN((short)esdtrack->GetTPCsignalN()); //due to bug in aliesdtrack class   
      trackCopy->SetTPCsignalS(esdtrack->GetTPCsignalSigma()); 

      trackCopy->SetTPCClusterMap(esdtrack->GetTPCClusterMap());
      trackCopy->SetTPCSharedMap(esdtrack->GetTPCSharedMap());

      double xtpc[3];
      esdtrack->GetInnerXYZ(xtpc);
      xtpc[2] -= fV1[2];
      trackCopy->SetNominalTPCEntrancePoint(xtpc);

      esdtrack->GetOuterXYZ(xtpc);
      xtpc[2] -= fV1[2];
      trackCopy->SetNominalTPCExitPoint(xtpc);

      int indexes[3];
      for (int ik=0; ik<3; ik++) {
	indexes[ik] = esdtrack->GetKinkIndex(ik);
      }
      trackCopy->SetKinkIndexes(indexes);
      //decision if we want this track
      //if we using diffrent labels we want that this label was use for first time 
      //if we use hidden info we want to have match between sim data and ESD
      if (tGoodMomentum==true)
	{
	  hbtEvent->TrackCollection()->push_back(trackCopy);//adding track to analysis
	  realnofTracks++;//real number of tracks
	  //	  delete trackCopy;
	}
      else
	{
	  delete  trackCopy;
	}
		
    }

  hbtEvent->SetNumberOfTracks(realnofTracks);//setting number of track which we read in event	

  AliCentrality *cent = fEvent->GetCentrality();
  if (cent) {
    hbtEvent->SetCentralityV0(cent->GetCentralityPercentile("V0M"));
    //    hbtEvent->SetCentralityFMD(cent->GetCentralityPercentile("FMD"));
    hbtEvent->SetCentralitySPD1(cent->GetCentralityPercentile("CL1"));
    //    hbtEvent->SetCentralityTrk(cent->GetCentralityPercentile("TRK"));

    printf("  FemtoReader Got Event with %f %f %f %f\n", cent->GetCentralityPercentile("V0M"), 0.0, cent->GetCentralityPercentile("CL1"), 0.0);
  }

  if (fEstEventMult == kGlobalCount) 
    hbtEvent->SetNormalizedMult(tNormMult);
  else if (fEstEventMult == kTracklet)
    hbtEvent->SetNormalizedMult(tTracklet);
  else if (fEstEventMult == kITSTPC)
    hbtEvent->SetNormalizedMult(tITSTPC);
  else if (fEstEventMult == kITSPure)
    hbtEvent->SetNormalizedMult(tITSPure);
  else if (fEstEventMult == kSPDLayer1)
    hbtEvent->SetNormalizedMult(fEvent->GetMultiplicity()->GetNumberOfITSClusters(1));
  else if (fEstEventMult == kV0Centrality) {
    // centrality between 0 (central) and 1 (very peripheral)

    if (cent) {
      if (cent->GetCentralityPercentile("V0M") < 0.0)
	hbtEvent->SetNormalizedMult(-1);
      else
	hbtEvent->SetNormalizedMult(lrint(10.0*cent->GetCentralityPercentile("V0M")));
      printf ("Set Centrality %i %f %li\n", hbtEvent->UncorrectedNumberOfPrimaries(), 
	      10.0*cent->GetCentralityPercentile("V0M"), lrint(10.0*cent->GetCentralityPercentile("V0M")));
    }
  }

  if (tNormMultPos > tNormMultNeg)
    hbtEvent->SetZDCParticipants(tNormMultPos);
  else
    hbtEvent->SetZDCParticipants(tNormMultNeg);

  fCurEvent++;	
  cout<<"end of reading nt "<<nofTracks<<" real number "<<realnofTracks<<endl;

  return hbtEvent; 
}
//___________________
void AliFemtoEventReaderESDChain::SetESDSource(AliESDEvent *aESD)
{
  // The chain loads the ESD for us
  // You must provide the address where it can be found
  fEvent = aESD;
}
//___________________
// void AliFemtoEventReaderESDChain::SetESDfriendSource(AliESDfriend *aFriend)
// {
//   // We need the ESD tree to obtain
//   // information about the friend file location
//   fEventFriend = aFriend;
// }

//____________________________________________________________________
Float_t AliFemtoEventReaderESDChain::GetSigmaToVertex(double *impact, double *covar)
{
  // Calculates the number of sigma to the vertex.

  Float_t b[2];
  Float_t bRes[2];
  Float_t bCov[3];

  b[0] = impact[0];
  b[1] = impact[1];
  bCov[0] = covar[0];
  bCov[1] = covar[1];
  bCov[2] = covar[2];

  bRes[0] = TMath::Sqrt(bCov[0]);
  bRes[1] = TMath::Sqrt(bCov[2]);

  // -----------------------------------
  // How to get to a n-sigma cut?
  //
  // The accumulated statistics from 0 to d is
  //
  // ->  Erf(d/Sqrt(2)) for a 1-dim gauss (d = n_sigma)
  // ->  1 - Exp(-d**2) for a 2-dim gauss (d*d = dx*dx + dy*dy != n_sigma)
  //
  // It means that for a 2-dim gauss: n_sigma(d) = Sqrt(2)*ErfInv(1 - Exp((-x**2)/2)
  // Can this be expressed in a different way?

  if (bRes[0] == 0 || bRes[1] ==0)
    return -1;

  Float_t d = TMath::Sqrt(TMath::Power(b[0]/bRes[0],2) + TMath::Power(b[1]/bRes[1],2));

  // stupid rounding problem screws up everything:
  // if d is too big, TMath::Exp(...) gets 0, and TMath::ErfInverse(1) that should be infinite, gets 0 :(
  if (TMath::Exp(-d * d / 2) < 1e-10)
    return 1000;

  d = TMath::ErfInverse(1 - TMath::Exp(-d * d / 2)) * TMath::Sqrt(2);
  return d;
}

void AliFemtoEventReaderESDChain::SetReadTrackType(ReadTrackType aType)
{
  fTrackType = aType;
}







