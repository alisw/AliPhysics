////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoEventReaderESDKine - the reader class for the Alice ESD              ///
/// Reads in ESD information and converts it into internal AliFemtoEvent     ///
/// Reads in AliESDfriend to create shared hit/quality information           ///
/// Authors: Marek Chojnacki mchojnacki@knf.pw.edu.pl                        ///
///          Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////

/*
 *$Id$
 *$Log$
 *Revision 1.1  2007/05/25 12:42:54  akisiel
 *Adding a reader for the Kine information
 *
 *Revision 1.2  2007/05/22 09:01:42  akisiel
 *Add the possibiloity to save cut settings in the ROOT file
 *
 *Revision 1.1  2007/05/16 10:22:11  akisiel
 *Making the directory structure of AliFemto flat. All files go into one common directory
 *
 *Revision 1.5  2007/05/03 09:45:20  akisiel
 *Fixing Effective C++ warnings
 *
 *Revision 1.4  2007/04/27 07:28:34  akisiel
 *Remove event number reading due to interface changes
 *
 *Revision 1.3  2007/04/27 07:25:16  akisiel
 *Make revisions needed for compilation from the main AliRoot tree
 *
 *Revision 1.1.1.1  2007/04/25 15:38:41  panos
 *Importing the HBT code dir
 *
 */

#include "AliFemtoEventReaderESDKine.h"

#include "TFile.h"
#include "TChain.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliStack.h"
#include "AliAODParticle.h"
#include "TParticle.h"

//#include "TSystem.h"

#include "AliFmPhysicalHelixD.h"
#include "AliFmThreeVectorF.h"

#include "SystemOfUnits.h"

#include "AliFemtoEvent.h"

ClassImp(AliFemtoEventReaderESDKine)

#if !(ST_NO_NAMESPACES)
  using namespace units;
#endif

using namespace std;
//____________________________
//constructor with 0 parameters , look at default settings 
AliFemtoEventReaderESDKine::AliFemtoEventReaderESDKine():
  fInputFile(" "),
  fFileName(" "),
  fConstrained(true),
  fNumberofEvent(0),
  fCurEvent(0),
  fCurRLEvent(0),
  fTree(0x0),
  fEvent(0x0),
  fRunLoader(0x0)
{
}

AliFemtoEventReaderESDKine::AliFemtoEventReaderESDKine(const AliFemtoEventReaderESDKine &aReader) :
  fInputFile(" "),
  fFileName(" "),
  fConstrained(true),
  fNumberofEvent(0),
  fCurEvent(0),
  fCurRLEvent(0),
  fTree(0x0),
  fEvent(0x0),
  fRunLoader(0x0)
{
  // copy constructor
  fInputFile = aReader.fInputFile;
  fFileName  = aReader.fFileName;
  fConstrained = aReader.fConstrained;
  fNumberofEvent = aReader.fNumberofEvent;
  fCurEvent = aReader.fCurEvent;
  fEvent = new AliESDEvent();
}
//__________________
//Destructor
AliFemtoEventReaderESDKine::~AliFemtoEventReaderESDKine()
{
  // destructor
  //delete fListOfFiles;
  delete fTree;
  delete fEvent;
  if (fRunLoader) delete fRunLoader;
}

//__________________
AliFemtoEventReaderESDKine& AliFemtoEventReaderESDKine::operator=(const AliFemtoEventReaderESDKine& aReader)
{
  // assignment operator
  if (this == &aReader)
    return *this;

  fInputFile = aReader.fInputFile;
  fFileName  = aReader.fFileName;
  fConstrained = aReader.fConstrained;
  fNumberofEvent = aReader.fNumberofEvent;
  fCurEvent = aReader.fCurEvent;
  fCurRLEvent = aReader.fCurRLEvent;
  if (fTree) delete fTree;
  //  fTree = aReader.fTree->CloneTree();
  if (fEvent) delete fEvent;
  fEvent = new AliESDEvent();
  if (fRunLoader) delete fRunLoader;
  fRunLoader = new AliRunLoader(*aReader.fRunLoader);

  return *this;
}
//__________________
AliFemtoString AliFemtoEventReaderESDKine::Report() const
{
  // create reader report
  AliFemtoString temp = "\n This is the AliFemtoEventReaderESDKine\n";
  return temp;
}

//__________________
void AliFemtoEventReaderESDKine::SetInputFile(const char* inputFile)
{
  //setting the name of file where names of ESD file are written 
  //it takes only this files which have good trees
  char buffer[256];
  fInputFile=string(inputFile);
  cout<<"Input File set on "<<fInputFile<<endl;
  ifstream infile(inputFile);

  fTree = new TChain("esdTree");

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
		  fTree->AddFile(buffer);
		  delete tree;
		}
	      esdFile->Close();	
	    }
	  delete esdFile;
	}
    }
}

void AliFemtoEventReaderESDKine::SetConstrained(const bool constrained)
{
  fConstrained=constrained;
}

bool AliFemtoEventReaderESDKine::GetConstrained() const
{
  return fConstrained;
}

AliFemtoEvent* AliFemtoEventReaderESDKine::ReturnHbtEvent()
{
  // read in a next hbt event from the chain
  // convert it to AliFemtoEvent and return
  // for further analysis
  AliFemtoEvent *hbtEvent = 0;
  TString tGAliceFilename;

  if (fCurEvent==fNumberofEvent)//open next file  
    {
      if (fNumberofEvent == 0) {
	fEvent=new AliESDEvent();
		
	  //ESD data
// 	  fEsdFile=TFile::Open(fFileName.c_str(),"READ");
// 	  fTree = (TTree*) fEsdFile->Get("esdTree");			

	  fTree->SetBranchStatus("MuonTracks*",0);
	  fTree->SetBranchStatus("PmdTracks*",0);
	  fTree->SetBranchStatus("TrdTracks*",0);
	  fTree->SetBranchStatus("V0s*",0);
	  fTree->SetBranchStatus("Cascades*",0);
	  fTree->SetBranchStatus("Kinks*",0);
	  fTree->SetBranchStatus("CaloClusters*",0);
	  fTree->SetBranchStatus("AliRawDataErrorLogs*",0);
	  fTree->SetBranchStatus("ESDfriend*",0);
	  fEvent->ReadFromTree(fTree);

// 	  chain->SetBranchStatus("*",0);
// 	  chain->SetBranchStatus("fUniqueID",1);
// 	  chain->SetBranchStatus("fTracks",1);
// 	  chain->SetBranchStatus("fTracks.*",1);
// 	  chain->SetBranchStatus("fTracks.fTPCindex[160]",1);
//	  fTree->SetBranchStatus("fTracks.fCalibContainer",0);


	fNumberofEvent=fTree->GetEntries();

	if (fNumberofEvent == 0) {
	  cout<<"no event in input "<<endl;
	  fReaderStatus=1;
	  return hbtEvent; 
	}

	cout<<"Number of Entries in the input "<<fNumberofEvent<<endl;
	fCurEvent=0;
	// simulation data reading setup
	
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
  //  vector<int> tLabelTable;//to check labels
  
  cout << "fFileName is " << fFileName.Data() << endl;
  cout << "Current file is " << fTree->GetCurrentFile()->GetName() << endl;
  if (fFileName.CompareTo(fTree->GetCurrentFile()->GetName())) {
    fFileName = fTree->GetCurrentFile()->GetName();
    tGAliceFilename = fFileName;
    tGAliceFilename.ReplaceAll("AliESDs","galice");
    cout << "Reading RunLoader from " << tGAliceFilename.Data() << endl;
    if (fRunLoader) delete fRunLoader;
    fRunLoader = AliRunLoader::Open(tGAliceFilename.Data());
    if (fRunLoader==0x0)
      {
	cout << "No Kine tree in file " << tGAliceFilename.Data() << endl;
	exit(0);
      }
    if(fRunLoader->LoadHeader())
      {
	cout << "Could not read RunLoader header in file " << tGAliceFilename.Data() << endl;
	exit(0);
      }
    fRunLoader->LoadKinematics();
    fCurRLEvent = 0;
  }

  fRunLoader->GetEvent(fCurRLEvent);
  AliStack* tStack = 0x0;
  tStack = fRunLoader->Stack();
	
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
      trackCopy->SetP(v);//setting momentum
      trackCopy->SetPt(sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]));
      const AliFmThreeVectorD ktP(pxyz[0],pxyz[1],pxyz[2]);
      if (ktP.mag() == 0) {
	delete trackCopy;
	continue;
      }
      const AliFmThreeVectorD kOrigin(fV1[0],fV1[1],fV1[2]);
      //setting helix I do not if it is ok
      AliFmPhysicalHelixD helix(ktP,kOrigin,(double)(fEvent->GetMagneticField())*kilogauss,(double)(trackCopy->Charge())); 
      trackCopy->SetHelix(helix);
	    	
      trackCopy->SetTrackId(esdtrack->GetID());
      trackCopy->SetFlags(esdtrack->GetStatus());
      trackCopy->SetLabel(esdtrack->GetLabel());
		
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

      trackCopy->SetTPCClusterMap(esdtrack->GetTPCClusterMap());
      trackCopy->SetTPCSharedMap(esdtrack->GetTPCSharedMap());

      // Fill the hidden information with the simulated data
      TParticle *tPart = tStack->Particle(TMath::Abs(esdtrack->GetLabel()));
      AliAODParticle* tParticle= new AliAODParticle(*tPart,i);
      AliFemtoModelHiddenInfo *tInfo = new AliFemtoModelHiddenInfo();
      tInfo->SetPDGPid(tParticle->GetMostProbable());
      tInfo->SetTrueMomentum(tParticle->Px(), tParticle->Py(), tParticle->Pz());
      Double_t mass2 = (tParticle->E()*tParticle->E() -
			 tParticle->Px()*tParticle->Px() -
			 tParticle->Py()*tParticle->Py() -
			 tParticle->Pz()*tParticle->Pz());
      if (mass2>0.0)
	tInfo->SetMass(TMath::Sqrt(mass2));
      else 
	tInfo->SetMass(0.0);
      trackCopy->SetHiddenInfo(tInfo);

      //decision if we want this track
      //if we using diffrent labels we want that this label was use for first time 
      //if we use hidden info we want to have match between sim data and ESD
      if (tGoodMomentum==true)
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
  fCurRLEvent++;
  cout<<"end of reading nt "<<nofTracks<<" real number "<<realnofTracks<<endl;
  if (fCurEvent== fNumberofEvent)//if end of current file close all
    {   
      fTree->Reset(); 
      delete fTree;
    }
  return hbtEvent; 
}
