////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoEventReaderESD - the reader class for the Alice ESD              ///
/// Reads in ESD information and converts it into internal AliFemtoEvent     ///
/// Reads in AliESDfriend to create shared hit/quality information           ///
/// Authors: Marek Chojnacki mchojnacki@knf.pw.edu.pl                        ///
///          Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////

/*
 *$Id$
 *$Log$
 *Revision 1.2.2.2  2007/10/04 13:10:52  akisiel
 *Add Kink index storageAliFemtoEventReaderESD.cxx AliFemtoTrack.cxx AliFemtoTrack.h
 *
 *Revision 1.2.2.1  2007/09/30 11:38:59  akisiel
 *Adapt the readers to the new AliESDEvent structure
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

#include "AliFemtoEventReaderESD.h"

#include "TFile.h"
#include "TTree.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"

//#include "TSystem.h"

#include "AliFmPhysicalHelixD.h"
#include "AliFmThreeVectorF.h"

#include "SystemOfUnits.h"

#include "AliFemtoEvent.h"
#include "AliFemtoModelHiddenInfo.h"

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
  fReadInner(false),
  fNumberofEvent(0),
  fCurEvent(0),
  fTree(0x0),
  fEsdFile(0x0),
  fEvent(0x0)
{
  // default constructor
}

AliFemtoEventReaderESD::AliFemtoEventReaderESD(const AliFemtoEventReaderESD &aReader) :
  fInputFile(" "),
  fFileName(" "),
  fConstrained(true),
  fReadInner(false),
  fNumberofEvent(0),
  fCurEvent(0),
  fTree(0x0),
  fEsdFile(0x0),
  fEvent(0x0)
{
  // copy constructor
  fInputFile = aReader.fInputFile;
  fFileName  = aReader.fFileName;
  fConstrained = aReader.fConstrained;
  fReadInner = aReader.fReadInner;
  fNumberofEvent = aReader.fNumberofEvent;
  fCurEvent = aReader.fCurEvent;
  //  fTree = aReader.fTree->CloneTree();
  //  fEvent = new AliESD(*aReader.fEvent);
  fEvent = new AliESDEvent();
  fEsdFile = new TFile(aReader.fEsdFile->GetName());
}
//__________________
//Destructor
AliFemtoEventReaderESD::~AliFemtoEventReaderESD()
{
  // destructor
  //delete fListOfFiles;
  delete fTree;
  delete fEvent;
  delete fEsdFile;
}

//__________________
AliFemtoEventReaderESD& AliFemtoEventReaderESD::operator=(const AliFemtoEventReaderESD& aReader)
{
  // assignment operator
  if (this == &aReader)
    return *this;

  fInputFile = aReader.fInputFile;
  fFileName  = aReader.fFileName;
  fConstrained = aReader.fConstrained;
  fReadInner = aReader.fReadInner;
  fNumberofEvent = aReader.fNumberofEvent;
  fCurEvent = aReader.fCurEvent;
  if (fTree) delete fTree;
  //  fTree = aReader.fTree->CloneTree();
  if (fEvent) delete fEvent;
  fEvent = new AliESDEvent();
  if (fEsdFile) delete fEsdFile;
  fEsdFile = new TFile(aReader.fEsdFile->GetName());

  return *this;
}
//__________________
AliFemtoString AliFemtoEventReaderESD::Report()
{
  // create reader report
  AliFemtoString temp = "\n This is the AliFemtoEventReaderESD\n";
  return temp;
}

//__________________
void AliFemtoEventReaderESD::SetInputFile(const char* inputFile)
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

void AliFemtoEventReaderESD::SetConstrained(const bool constrained)
{
  fConstrained=constrained;
}

bool AliFemtoEventReaderESD::GetConstrained() const
{
  return fConstrained;
}

void AliFemtoEventReaderESD::SetReadTPCInner(const bool readinner)
{
  fReadInner=readinner;
}

bool AliFemtoEventReaderESD::GetReadTPCInner() const
{
  return fReadInner;
}

AliFemtoEvent* AliFemtoEventReaderESD::ReturnHbtEvent()
{
  // read in a next hbt event from the chain
  // convert it to AliFemtoEvent and return
  // for further analysis
  AliFemtoEvent *hbtEvent = 0;

  if (fCurEvent==fNumberofEvent)//open next file  
    {
      if(fNumberofEvent==0)	
	{
	  //	  delete fEvent;//added 1.04.2007
	  fEvent=new AliESDEvent();
	  //	  delete fTree;
	  //fTree=0;
	  //	  delete fEsdFile;
		
	  //ESD data
	  //	  fEsdFile=TFile::Open(fFileName.c_str(),"READ");
	  //	  fTree = (TTree*) fEsdFile->Get("esdTree");			
	  //	  fTree->SetBranchAddress("ESD", &fEvent);			
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
  cout << "Read event " << fEvent << " from file " << fTree << endl;
  //  vector<int> tLabelTable;//to check labels
	
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
  cout << "Event has " << nofTracks << " tracks " << endl;

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
      if (fReadInner == true) {
	
	if (esdtrack->GetTPCInnerParam()) {
	  AliExternalTrackParam *param = new AliExternalTrackParam(*esdtrack->GetTPCInnerParam());
	  param->PropagateToDCA(fEvent->GetPrimaryVertex(), (fEvent->GetMagneticField()), 10000);
	  param->GetPxPyPz(pxyz);//reading noconstarined momentum
	  delete param;

	  AliFemtoModelHiddenInfo *tInfo = new AliFemtoModelHiddenInfo();
	  tInfo->SetPDGPid(211);
	  tInfo->SetTrueMomentum(pxyz[0], pxyz[1], pxyz[2]);
	  tInfo->SetMass(0.13957);
	  trackCopy->SetHiddenInfo(tInfo);
	}
      }
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
      const AliFmThreeVectorD origin(fV1[0],fV1[1],fV1[2]);
      //setting helix I do not if it is ok
      AliFmPhysicalHelixD helix(ktP,origin,(double)(fEvent->GetMagneticField())*kilogauss,(double)(trackCopy->Charge())); 
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
	}
      else
	{
	  delete  trackCopy;
	}
		
    }

  hbtEvent->SetNumberOfTracks(realnofTracks);//setting number of track which we read in event	
  fCurEvent++;	
  cout<<"end of reading nt "<<nofTracks<<" real number "<<realnofTracks<<endl;
//   if (fCurEvent== fNumberofEvent)//if end of current file close all
//     {   
//       fTree->Reset(); 
//       delete fTree;
//       fEsdFile->Close();
//     }
  return hbtEvent; 
}









