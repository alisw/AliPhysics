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

#include "AliFmPhysicalHelixD.h"
#include "AliFmThreeVectorF.h"

#include "SystemOfUnits.h"

#include "AliFemtoEvent.h"

ClassImp(AliFemtoEventReaderESDChain)

#if !(ST_NO_NAMESPACES)
  using namespace units;
#endif

using namespace std;
//____________________________
AliFemtoEventReaderESDChain::AliFemtoEventReaderESDChain():
  fFileName(" "),
  fConstrained(true),
  fNumberofEvent(0),
  fCurEvent(0),
  fCurFile(0),
  fEvent(0x0)
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
  fFileName(" "),
  fConstrained(true),
  fNumberofEvent(0),
  fCurEvent(0),
  fCurFile(0),
  fEvent(0x0)
{
  // Copy constructor
  fConstrained = aReader.fConstrained;
  fNumberofEvent = aReader.fNumberofEvent;
  fCurEvent = aReader.fCurEvent;
  fCurFile = aReader.fCurFile;
  //  fEvent = new AliESD(*aReader.fEvent);
  fEvent = new AliESDEvent();
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
}

//__________________
AliFemtoEventReaderESDChain& AliFemtoEventReaderESDChain::operator=(const AliFemtoEventReaderESDChain& aReader)
{
  // Assignment operator
  if (this == &aReader)
    return *this;

  fConstrained = aReader.fConstrained;
  fNumberofEvent = aReader.fNumberofEvent;
  fCurEvent = aReader.fCurEvent;
  fCurFile = aReader.fCurFile;
  if (fEvent) delete fEvent;
  fEvent = new AliESDEvent();

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
  vector<int> tLabelTable;//to check labels
	
  hbtEvent = new AliFemtoEvent;
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
  hbtEvent->SetTriggerCluster(fEvent->GetTriggerCluster());
	
  //Vertex
  double fV1[3];
  fEvent->GetVertex()->GetXYZ(fV1);

  AliFmThreeVectorF vertex(fV1[0],fV1[1],fV1[2]);
  hbtEvent->SetPrimVertPos(vertex);
	
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

  for (int i=0;i<nofTracks;i++)
    {
      bool  tGoodMomentum=true; //flaga to chcek if we can read momentum of this track
		
      AliFemtoTrack* trackCopy = new AliFemtoTrack();	
      const AliESDtrack *esdtrack=fEvent->GetTrack(i);//getting next track
      //      const AliESDfriendTrack *tESDfriendTrack = esdtrack->GetFriendTrack();

      trackCopy->SetCharge((short)esdtrack->GetSign());

      //in aliroot we have AliPID 
      //0-electron 1-muon 2-pion 3-kaon 4-proton 5-photon 6-pi0 7-neutron 8-kaon0 9-eleCon   
      //we use only 5 first
      double esdpid[5];
      esdtrack->GetESDpid(esdpid);
      trackCopy->SetPidProbElectron(esdpid[0]);
      trackCopy->SetPidProbMuon(esdpid[1]);
      trackCopy->SetPidProbPion(esdpid[2]);
      trackCopy->SetPidProbKaon(esdpid[3]);
      trackCopy->SetPidProbProton(esdpid[4]);
						
      double pxyz[3];
      if (fConstrained==true)		    
	tGoodMomentum=esdtrack->GetConstrainedPxPyPz(pxyz); //reading constrained momentum
      else
	tGoodMomentum=esdtrack->GetPxPyPz(pxyz);//reading noconstarined momentum

      AliFemtoThreeVector v(pxyz[0],pxyz[1],pxyz[2]);
      if (v.mag() == 0) {
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
	    	
      trackCopy->SetTrackId(esdtrack->GetID());
      trackCopy->SetFlags(esdtrack->GetStatus());
      //trackCopy->SetLabel(esdtrack->GetLabel());
		
      //some stuff which could be useful 
      float impact[2];
      float covimpact[3];
      esdtrack->GetImpactParameters(impact,covimpact);
      trackCopy->SetImpactD(impact[0]);
      trackCopy->SetImpactZ(impact[1]);
      trackCopy->SetCdd(covimpact[0]);
      trackCopy->SetCdz(covimpact[1]);
      trackCopy->SetCzz(covimpact[2]);
      trackCopy->SetITSchi2(esdtrack->GetITSchi2());    
      trackCopy->SetITSncls(esdtrack->GetNcls(0));     
      trackCopy->SetTPCchi2(esdtrack->GetTPCchi2());       
      trackCopy->SetTPCncls(esdtrack->GetTPCNcls());       
      trackCopy->SetTPCnclsF(esdtrack->GetTPCNclsF());      
      trackCopy->SetTPCsignalN((short)esdtrack->GetTPCsignalN()); //due to bug in aliesdtrack class   
      trackCopy->SetTPCsignalS(esdtrack->GetTPCsignalSigma()); 

//       // Fill cluster per padrow information
//       Int_t tTrackIndices[AliESDfriendTrack::kMaxTPCcluster];
//       Int_t tNClusters = esdtrack->GetTPCclusters(tTrackIndices);
//       for (int tNcl=0; tNcl<AliESDfriendTrack::kMaxTPCcluster; tNcl++) {
// 	if (tTrackIndices[tNcl] > 0)
// 	  trackCopy->SetTPCcluster(tNcl, 1);
// 	else
// 	  trackCopy->SetTPCcluster(tNcl, 0);
//       }
      
//       // Fill shared cluster information
//       list<Int_t>::iterator tClustIter;

//       for (int tNcl=0; tNcl<AliESDfriendTrack::kMaxTPCcluster; tNcl++) {
// 	if (tTrackIndices[tNcl] > 0) {
// 	  tClustIter = find(fSharedList[tNcl]->begin(), fSharedList[tNcl]->end(), tTrackIndices[tNcl]);
// 	  if (tClustIter != fSharedList[tNcl]->end()) {
// 	    trackCopy->SetTPCshared(tNcl, 1);
// 	    cout << "Event next" <<  endl;
// 	    cout << "Track: " << i << endl;
// 	    cout << "Shared cluster: " << tNcl << " " << tTrackIndices[tNcl] << endl;
// 	  }
// 	  else {
// 	    trackCopy->SetTPCshared(tNcl, 0);
// 	  }
// 	}
//       }

      trackCopy->SetTPCClusterMap(esdtrack->GetTPCClusterMap());
      trackCopy->SetTPCSharedMap(esdtrack->GetTPCSharedMap());

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









