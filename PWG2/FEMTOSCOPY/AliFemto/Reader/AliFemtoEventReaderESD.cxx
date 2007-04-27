/*
 *$Id$
 *$Log$
 *Revision 1.1.1.1  2007/04/25 15:38:41  panos
 *Importing the HBT code dir
 *
 *Revision 1.5  2007-04-03 16:00:08  mchojnacki
 *Changes to iprove memory managing
 *
 *Revision 1.4  2007/03/13 15:30:03  mchojnacki
 *adding reader for simulated data
 *
 *Revision 1.3  2007/03/08 14:58:03  mchojnacki
 *adding some alice stuff
 *
 *Revision 1.2  2007/03/07 13:36:16  mchojnacki
 *Add some comments
 *
 *Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 *First version on CVS
 *
 */

#include "AliFemtoEventReaderESD.h"

#include "TFile.h"
#include "TTree.h"
#include "AliESD.h"
#include "AliESDtrack.h"

//#include "TSystem.h"

#include "Infrastructure/AliFmPhysicalHelixD.h"
#include "Infrastructure/AliFmThreeVectorF.h"

#include "Base/SystemOfUnits.h"

#include "Infrastructure/AliFemtoEvent.h"

ClassImp(AliFemtoEventReaderESD)

#if !(ST_NO_NAMESPACES)
  using namespace units;
#endif

using namespace std;
//____________________________
//constructor with 0 parameters , look at default settings 
AliFemtoEventReaderESD::AliFemtoEventReaderESD():
  fInputFile(" "),
  fFileName(" "),
  fConstrained(true),
  fNumberofEvent(0),
  fCurEvent(0),
  fCurFile(0),
  fTree(0x0),
  fEvent(0x0),
  fEsdFile(0x0),
  fEventFriend(0)
{
  fClusterPerPadrow = (list<Int_t> **) malloc(sizeof(list<Int_t> *) * AliESDfriendTrack::kMaxTPCcluster);
  for (int tPad=0; tPad<AliESDfriendTrack::kMaxTPCcluster; tPad++) {
    fClusterPerPadrow[tPad] = new list<Int_t>();
  }
  fSharedList = (list<Int_t> **) malloc(sizeof(list<Int_t> *) * AliESDfriendTrack::kMaxTPCcluster);
  for (int tPad=0; tPad<AliESDfriendTrack::kMaxTPCcluster; tPad++) {
    fSharedList[tPad] = new list<Int_t>();
  }
}

//__________________
//Destructor
AliFemtoEventReaderESD::~AliFemtoEventReaderESD()
{
  //delete fListOfFiles;
  delete fTree;
  delete fEvent;
  delete fEsdFile;

  for (int tPad=0; tPad<AliESDfriendTrack::kMaxTPCcluster; tPad++) {
    fClusterPerPadrow[tPad]->clear();
    delete fClusterPerPadrow[tPad];
  }
  delete [] fClusterPerPadrow;
  for (int tPad=0; tPad<AliESDfriendTrack::kMaxTPCcluster; tPad++) {
    fSharedList[tPad]->clear();
    delete fSharedList[tPad];
  }
  delete [] fSharedList;
}

//__________________
AliFemtoString AliFemtoEventReaderESD::Report()
{
  AliFemtoString temp = "\n This is the AliFemtoEventReaderESD\n";
  return temp;
}

//__________________
//setting the name of file where names of ESD file are written 
//it takes only this files which have good trees
void AliFemtoEventReaderESD::SetInputFile(const char* inputFile)
{
  char buffer[256];
  fInputFile=string(inputFile);
  cout<<"Input File set on "<<fInputFile<<endl;
  ifstream infile(inputFile);
  if(infile.good()==true)
    { 
      //checking if all give files have good tree inside
      while (infile.eof()==false)
	{
	  infile.getline(buffer,256);
	  //ifstream test_file(buffer);
	  TFile *esdFile=TFile::Open(buffer,"READ");
	  if (esdFile!=0x0)
	    {	
	      TTree* tree = (TTree*) esdFile->Get("esdTree");
	      if (tree!=0x0)
		{
		  cout<<"putting file  "<<string(buffer)<<" into analysis"<<endl;
		  fListOfFiles.push_back(string(buffer));
		  delete tree;
		}
	      esdFile->Close();	
	    }
	  delete esdFile;
	}
    }
}

//setting the next file to read	
bool AliFemtoEventReaderESD::GetNextFile()
{ 	
  if (fCurFile>=fListOfFiles.size())
    return false;
  fFileName=fListOfFiles.at(fCurFile);	
  cout<<"FileName set on "<<fFileName<<" "<<fCurFile<<endl;

  fCurFile++;
  return true;
}
void AliFemtoEventReaderESD::SetConstrained(const bool constrained)
{
  fConstrained=constrained;
}

bool AliFemtoEventReaderESD::GetConstrained() const
{
  return fConstrained;
}

AliFemtoEvent* AliFemtoEventReaderESD::ReturnHbtEvent()
{
  AliFemtoEvent *hbtEvent = 0;
  string tFriendFileName;

  if (fCurEvent==fNumberofEvent)//open next file  
    {
      cout<<"next file"<<endl;
      if(GetNextFile())	
	{
	  delete fEventFriend;
	  fEventFriend = 0;
	  delete fEvent;//added 1.04.2007
	  fEvent=new AliESD();
	  //	  delete fTree;
	  //fTree=0;
	  delete fEsdFile;
		
	  //ESD data
	  fEsdFile=TFile::Open(fFileName.c_str(),"READ");
	  fTree = (TTree*) fEsdFile->Get("esdTree");			
	  fTree->SetBranchAddress("ESD", &fEvent);			

	  // Attach the friend tree with additional information
 	  tFriendFileName = fFileName;
 	  tFriendFileName.insert(tFriendFileName.find("s.root"),"friend");
 	  cout << "Reading friend " << tFriendFileName.c_str() << endl;;
  	  fTree->AddFriend("esdFriendTree",tFriendFileName.c_str());
  	  fTree->SetBranchAddress("ESDfriend",&fEventFriend);

// 	  chain->SetBranchStatus("*",0);
// 	  chain->SetBranchStatus("fUniqueID",1);
// 	  chain->SetBranchStatus("fTracks",1);
// 	  chain->SetBranchStatus("fTracks.*",1);
// 	  chain->SetBranchStatus("fTracks.fTPCindex[160]",1);
	  fTree->SetBranchStatus("fTracks.fCalibContainer",0);


	  fNumberofEvent=fTree->GetEntries();
	  cout<<"Number of Entries in file "<<fNumberofEvent<<endl;
	  fCurEvent=0;
	  //sim data
	}
      else //no more data to read
	{
	  cout<<"no more files "<<hbtEvent<<endl;
	  fReaderStatus=1;
	  return hbtEvent; 
	}
    }		
  cout<<"starting to read event "<<fCurEvent<<endl;
  fTree->GetEvent(fCurEvent);//getting next event
  fEvent->SetESDfriend(fEventFriend);
  vector<int> label_table;//to check labels
	
  hbtEvent = new AliFemtoEvent;
  //setting basic things
  hbtEvent->SetEventNumber(fEvent->GetEventNumber());
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

  // Clear the shared cluster list
  for (int tPad=0; tPad<AliESDfriendTrack::kMaxTPCcluster; tPad++) {
    fClusterPerPadrow[tPad]->clear();
  }
  for (int tPad=0; tPad<AliESDfriendTrack::kMaxTPCcluster; tPad++) {
    fSharedList[tPad]->clear();
  }


  for (int i=0;i<nofTracks;i++) {
    const AliESDtrack *esdtrack=fEvent->GetTrack(i);//getting next track

    list<Int_t>::iterator tClustIter;

    Int_t tTrackIndices[AliESDfriendTrack::kMaxTPCcluster];
    Int_t tNClusters = esdtrack->GetTPCclusters(tTrackIndices);
    for (int tNcl=0; tNcl<AliESDfriendTrack::kMaxTPCcluster; tNcl++) {
      if (tTrackIndices[tNcl] >= 0) {
	tClustIter = find(fClusterPerPadrow[tNcl]->begin(), fClusterPerPadrow[tNcl]->end(), tTrackIndices[tNcl]);
	if (tClustIter == fClusterPerPadrow[tNcl]->end()) {
	  fClusterPerPadrow[tNcl]->push_back(tTrackIndices[tNcl]);
	}
	else {
	  fSharedList[tNcl]->push_back(tTrackIndices[tNcl]);
	}
      }
    }
      
  }

  for (int i=0;i<nofTracks;i++)
    {
      bool  good_momentum=true; //flaga to chcek if we can read momentum of this track
		
      AliFemtoTrack* trackCopy = new AliFemtoTrack();	
      const AliESDtrack *esdtrack=fEvent->GetTrack(i);//getting next track
      const AliESDfriendTrack *tESDfriendTrack = esdtrack->GetFriendTrack();

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
	good_momentum=esdtrack->GetConstrainedPxPyPz(pxyz); //reading constrained momentum
      else
	good_momentum=esdtrack->GetPxPyPz(pxyz);//reading noconstarined momentum
      AliFemtoThreeVector v(pxyz[0],pxyz[1],pxyz[2]);
      trackCopy->SetP(v);//setting momentum
      trackCopy->SetPt(sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]));
      const AliFmThreeVectorD p(pxyz[0],pxyz[1],pxyz[2]);
      if (p.mag() == 0) {
	delete trackCopy;
	continue;
      }
      const AliFmThreeVectorD origin(fV1[0],fV1[1],fV1[2]);
      //setting helix I do not if it is ok
      AliFmPhysicalHelixD helix(p,origin,(double)(fEvent->GetMagneticField())*kilogauss,(double)(trackCopy->Charge())); 
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

      // Fill cluster per padrow information
      Int_t tTrackIndices[AliESDfriendTrack::kMaxTPCcluster];
      Int_t tNClusters = esdtrack->GetTPCclusters(tTrackIndices);
      for (int tNcl=0; tNcl<AliESDfriendTrack::kMaxTPCcluster; tNcl++) {
	if (tTrackIndices[tNcl] > 0)
	  trackCopy->SetTPCcluster(tNcl, 1);
	else
	  trackCopy->SetTPCcluster(tNcl, 0);
      }
      
      // Fill shared cluster information
      list<Int_t>::iterator tClustIter;

      for (int tNcl=0; tNcl<AliESDfriendTrack::kMaxTPCcluster; tNcl++) {
	if (tTrackIndices[tNcl] > 0) {
	  tClustIter = find(fSharedList[tNcl]->begin(), fSharedList[tNcl]->end(), tTrackIndices[tNcl]);
	  if (tClustIter != fSharedList[tNcl]->end()) {
	    trackCopy->SetTPCshared(tNcl, 1);
	    cout << "Event next" <<  endl;
	    cout << "Track: " << i << endl;
	    cout << "Shared cluster: " << tNcl << " " << tTrackIndices[tNcl] << endl;
	  }
	  else {
	    trackCopy->SetTPCshared(tNcl, 0);
	  }
	}
      }

      //decision if we want this track
      //if we using diffrent labels we want that this label was use for first time 
      //if we use hidden info we want to have match between sim data and ESD
      if (good_momentum==true)
	{
	  hbtEvent->TrackCollection()->push_back(trackCopy);//adding track to analysis
	  realnofTracks++;//real number of tracks
	}
      else
	{
	  delete  trackCopy;
	}
		
    }

  hbtEvent->SetNumberOfTracks(realnofTracks);//setting number of track which we read in event	
  fCurEvent++;	
  cout<<"end of reading nt "<<nofTracks<<" real number "<<realnofTracks<<endl;
  if (fCurEvent== fNumberofEvent)//if end of current file close all
    {   
      fTree->Reset(); 
      delete fTree;
      fEsdFile->Close();
    }
  return hbtEvent; 
}









