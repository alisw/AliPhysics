////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoEventReaderAOD - the reader class for the Alice AOD                //
// Reads in AOD information and converts it into internal AliFemtoEvent       //
// Authors: Marek Chojnacki mchojnacki@knf.pw.edu.pl                          //
//          Adam Kisiel kisiel@mps.ohio-state.edu                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoEventReaderAOD.h"

#include "TFile.h"
#include "TTree.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"

#include "AliFmPhysicalHelixD.h"
#include "AliFmThreeVectorF.h"

#include "SystemOfUnits.h"

#include "AliFemtoEvent.h"
#include "AliFemtoModelHiddenInfo.h"

ClassImp(AliFemtoEventReaderAOD)

#if !(ST_NO_NAMESPACES)
  using namespace units;
#endif

using namespace std;
//____________________________
//constructor with 0 parameters , look at default settings 
AliFemtoEventReaderAOD::AliFemtoEventReaderAOD():
  fInputFile(" "),
  fFileName(" "),
  fNumberofEvent(0),
  fCurEvent(0),
  fTree(0x0),
  fAodFile(0x0),
  fEvent(0x0),
  fAllTrue(160),
  fAllFalse(160),
  fFilterBit(0),
  fPWG2AODTracks(0x0)
{
  // default constructor
  fAllTrue.ResetAllBits(kTRUE);
  fAllFalse.ResetAllBits(kFALSE);
}

AliFemtoEventReaderAOD::AliFemtoEventReaderAOD(const AliFemtoEventReaderAOD &aReader) :
  fInputFile(" "),
  fFileName(" "),
  fNumberofEvent(0),
  fCurEvent(0),
  fTree(0x0),
  fAodFile(0x0),
  fEvent(0x0),
  fAllTrue(160),
  fAllFalse(160),
  fFilterBit(0),
  fPWG2AODTracks(0x0)
{
  // copy constructor
  fInputFile = aReader.fInputFile;
  fFileName  = aReader.fFileName;
  fNumberofEvent = aReader.fNumberofEvent;
  fCurEvent = aReader.fCurEvent;
  fEvent = new AliAODEvent();
  fAodFile = new TFile(aReader.fAodFile->GetName());
  fAllTrue.ResetAllBits(kTRUE);
  fAllFalse.ResetAllBits(kFALSE);
  fFilterBit = aReader.fFilterBit;
  fPWG2AODTracks = aReader.fPWG2AODTracks;
}
//__________________
//Destructor
AliFemtoEventReaderAOD::~AliFemtoEventReaderAOD()
{
  // destructor
  delete fTree;
  delete fEvent;
  delete fAodFile;
  if (fPWG2AODTracks) {
    fPWG2AODTracks->Delete();
    delete fPWG2AODTracks;
  }
}

//__________________
AliFemtoEventReaderAOD& AliFemtoEventReaderAOD::operator=(const AliFemtoEventReaderAOD& aReader)
{
  // assignment operator
  if (this == &aReader)
    return *this;

  fInputFile = aReader.fInputFile;
  fFileName  = aReader.fFileName;
  fNumberofEvent = aReader.fNumberofEvent;
  fCurEvent = aReader.fCurEvent;
  if (fTree) delete fTree;
  if (fEvent) delete fEvent;
  fEvent = new AliAODEvent();
  if (fAodFile) delete fAodFile;
  fAodFile = new TFile(aReader.fAodFile->GetName());
  fAllTrue.ResetAllBits(kTRUE);
  fAllFalse.ResetAllBits(kFALSE);
  fFilterBit = aReader.fFilterBit;
  fPWG2AODTracks = aReader.fPWG2AODTracks;

  return *this;
}
//__________________
AliFemtoString AliFemtoEventReaderAOD::Report()
{
  // create reader report
  AliFemtoString temp = "\n This is the AliFemtoEventReaderAOD\n";
  return temp;
}

//__________________
void AliFemtoEventReaderAOD::SetInputFile(const char* inputFile)
{
  //setting the name of file where names of AOD file are written 
  //it takes only this files which have good trees
  char buffer[256];
  fInputFile=string(inputFile);
  cout<<"Input File set on "<<fInputFile<<endl;
  ifstream infile(inputFile);

  fTree = new TChain("aodTree");

  if(infile.good()==true)
    { 
      //checking if all give files have good tree inside
      while (infile.eof()==false)
	{
	  infile.getline(buffer,256);
	  TFile *aodFile=TFile::Open(buffer,"READ");
	  if (aodFile!=0x0)
	    {	
	      TTree* tree = (TTree*) aodFile->Get("aodTree");
	      if (tree!=0x0)
		{
		  cout<<"putting file  "<<string(buffer)<<" into analysis"<<endl;
		  fTree->AddFile(buffer);
		  delete tree;
		}
	      aodFile->Close();	
	    }
	  delete aodFile;
	}
    }
}

AliFemtoEvent* AliFemtoEventReaderAOD::ReturnHbtEvent()
{
  // read in a next hbt event from the chain
  // convert it to AliFemtoEvent and return
  // for further analysis
  AliFemtoEvent *hbtEvent = 0;

  if (fCurEvent==fNumberofEvent)//open next file  
    {
      if(fNumberofEvent==0)	
	{
	  fEvent=new AliAODEvent();
	  fEvent->ReadFromTree(fTree);

	  // Check for the existence of the additional information
	  fPWG2AODTracks = (TClonesArray *) fEvent->GetList()->FindObject("pwg2aodtracks");

	  if (fPWG2AODTracks) {
	    cout << "Found additional PWG2 specific information in the AOD!" << endl;
	    cout << "Reading only tracks with the additional information" << endl;
	  }

	  fNumberofEvent=fTree->GetEntries();
	  cout<<"Number of Entries in file "<<fNumberofEvent<<endl;
	  fCurEvent=0;
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
	
  hbtEvent = new AliFemtoEvent;

  CopyAODtoFemtoEvent(hbtEvent);

  fCurEvent++;	

  return hbtEvent; 
}

void AliFemtoEventReaderAOD::CopyAODtoFemtoEvent(AliFemtoEvent *tEvent)
{
  // A function that reads in the AOD event
  // and transfers the neccessary information into
  // the internal AliFemtoEvent

  // setting global event characteristics
  tEvent->SetRunNumber(fEvent->GetRunNumber());
  tEvent->SetMagneticField(fEvent->GetMagneticField()*kilogauss);//to check if here is ok
  tEvent->SetZDCN1Energy(fEvent->GetZDCN1Energy());
  tEvent->SetZDCP1Energy(fEvent->GetZDCP1Energy());
  tEvent->SetZDCN2Energy(fEvent->GetZDCN2Energy());
  tEvent->SetZDCP2Energy(fEvent->GetZDCP2Energy());
  tEvent->SetZDCEMEnergy(fEvent->GetZDCEMEnergy(0));
  tEvent->SetZDCParticipants(-1);
  tEvent->SetTriggerMask(fEvent->GetTriggerMask());
  tEvent->SetTriggerCluster(fEvent->GetTriggerCluster());
	
  // Primary Vertex position
  double fV1[3];
  fEvent->GetPrimaryVertex()->GetPosition(fV1);

  AliFmThreeVectorF vertex(fV1[0],fV1[1],fV1[2]);
  tEvent->SetPrimVertPos(vertex);
	
  //starting to reading tracks
  int nofTracks=0;  //number of reconstructed tracks in event

  // Check to see whether the additional info exists
  if (fPWG2AODTracks)
    nofTracks=fPWG2AODTracks->GetEntries();
  else
    nofTracks=fEvent->GetNumberOfTracks();

  int realnofTracks=0;   // number of track which we use in a analysis
  cout << "Event has " << nofTracks << " tracks " << endl;

  for (int i=0;i<nofTracks;i++)
    {
      AliFemtoTrack* trackCopy = new AliFemtoTrack();	

      if (fPWG2AODTracks) {
	// Read tracks from the additional pwg2 specific AOD part
	// if they exist
	// Note that in that case all the AOD tracks without the 
	// additional information will be ignored !
	AliPWG2AODTrack *pwg2aodtrack = (AliPWG2AODTrack *) fPWG2AODTracks->At(i);

	// Getting the AOD track through the ref of the additional info
	AliAODTrack *aodtrack = pwg2aodtrack->GetRefAODTrack();	
	if (!aodtrack->TestFilterBit(fFilterBit))
	  continue;
	
	CopyAODtoFemtoTrack(aodtrack, trackCopy, pwg2aodtrack);
	
	double pxyz[3];
	aodtrack->PxPyPz(pxyz);//reading noconstarined momentum
	const AliFmThreeVectorD ktP(pxyz[0],pxyz[1],pxyz[2]);
	// Check the sanity of the tracks - reject zero momentum tracks
	if (ktP.mag() == 0) {
	  delete trackCopy;
	  continue;
	}
      }
      else {
	// No additional information exists
	// Read in the normal AliAODTracks 
	const AliAODTrack *aodtrack=fEvent->GetTrack(i); // getting the AODtrack directly
	
	if (!aodtrack->TestFilterBit(fFilterBit))
	  continue;
	
	CopyAODtoFemtoTrack(aodtrack, trackCopy, 0);
	
	double pxyz[3];
	aodtrack->PxPyPz(pxyz);//reading noconstarined momentum
	const AliFmThreeVectorD ktP(pxyz[0],pxyz[1],pxyz[2]);
	// Check the sanity of the tracks - reject zero momentum tracks
	if (ktP.mag() == 0) {
	  delete trackCopy;
	  continue;
	}
      }


      tEvent->TrackCollection()->push_back(trackCopy);//adding track to analysis
      realnofTracks++;//real number of tracks		
    }

  tEvent->SetNumberOfTracks(realnofTracks);//setting number of track which we read in event	

  cout<<"end of reading nt "<<nofTracks<<" real number "<<realnofTracks<<endl;
}

void AliFemtoEventReaderAOD::CopyAODtoFemtoTrack(const AliAODTrack *tAodTrack, 
						 AliFemtoTrack *tFemtoTrack, 
						 AliPWG2AODTrack *tPWG2AODTrack)
{
  // Copy the track information from the AOD into the internal AliFemtoTrack
  // If it exists, use the additional information from the PWG2 AOD

  // Primary Vertex position
  double fV1[3];
  fEvent->GetPrimaryVertex()->GetPosition(fV1);

  tFemtoTrack->SetCharge(tAodTrack->Charge());
  
  //in aliroot we have AliPID 
  //0-electron 1-muon 2-pion 3-kaon 4-proton 5-photon 6-pi0 7-neutron 8-kaon0 9-eleCon   
  //we use only 5 first

  // AOD pid has 10 components
  double aodpid[10];
  tAodTrack->GetPID(aodpid);
  tFemtoTrack->SetPidProbElectron(aodpid[0]);
  tFemtoTrack->SetPidProbMuon(aodpid[1]);
  tFemtoTrack->SetPidProbPion(aodpid[2]);
  tFemtoTrack->SetPidProbKaon(aodpid[3]);
  tFemtoTrack->SetPidProbProton(aodpid[4]);
						
  double pxyz[3];
  tAodTrack->PxPyPz(pxyz);//reading noconstrained momentum
  AliFemtoThreeVector v(pxyz[0],pxyz[1],pxyz[2]);
  tFemtoTrack->SetP(v);//setting momentum
  tFemtoTrack->SetPt(sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]));
  const AliFmThreeVectorD kOrigin(fV1[0],fV1[1],fV1[2]);
  //setting track helix 
  const AliFmThreeVectorD ktP(pxyz[0],pxyz[1],pxyz[2]);
  AliFmPhysicalHelixD helix(ktP,kOrigin,(double)(fEvent->GetMagneticField())*kilogauss,(double)(tFemtoTrack->Charge())); 
  tFemtoTrack->SetHelix(helix);
	    	
  // Flags
  tFemtoTrack->SetTrackId(tAodTrack->GetID());
  tFemtoTrack->SetFlags(1);
  tFemtoTrack->SetLabel(tAodTrack->GetLabel());
		
  // Track quality information 
  float covmat[6];
  tAodTrack->GetCovMatrix(covmat);
  tFemtoTrack->SetImpactD(covmat[0]);
  tFemtoTrack->SetImpactZ(covmat[2]);
  tFemtoTrack->SetCdd(covmat[3]);
  tFemtoTrack->SetCdz(covmat[4]);
  tFemtoTrack->SetCzz(covmat[5]);
  // This information is only available in the ESD
  // We put in fake values or reasonable estimates
  tFemtoTrack->SetITSchi2(tAodTrack->Chi2perNDF());    
  tFemtoTrack->SetITSncls(1);     
  tFemtoTrack->SetTPCchi2(tAodTrack->Chi2perNDF());       
  tFemtoTrack->SetTPCncls(1);       
  tFemtoTrack->SetTPCnclsF(-1);      
  tFemtoTrack->SetTPCsignalN(-1); 
  tFemtoTrack->SetTPCsignalS(-1); 

  if (tPWG2AODTrack) {
    // Copy the PWG2 specific information if it exists
    tFemtoTrack->SetTPCClusterMap(tPWG2AODTrack->GetTPCClusterMap());
    tFemtoTrack->SetTPCSharedMap(tPWG2AODTrack->GetTPCSharedMap());
    
    double xtpc[3] = {0,0,0};
    tPWG2AODTrack->GetTPCNominalEntrancePoint(xtpc);
    tFemtoTrack->SetNominalTPCEntrancePoint(xtpc);
    tPWG2AODTrack->GetTPCNominalExitPoint(xtpc);
    tFemtoTrack->SetNominalTPCExitPoint(xtpc);
  }
  else {
    // If not use dummy values
    tFemtoTrack->SetTPCClusterMap(fAllTrue);
    tFemtoTrack->SetTPCSharedMap(fAllFalse);
    
    double xtpc[3] = {0,0,0};
    tFemtoTrack->SetNominalTPCEntrancePoint(xtpc);
    tFemtoTrack->SetNominalTPCExitPoint(xtpc);
  }

  int indexes[3];
  for (int ik=0; ik<3; ik++) {
    indexes[ik] = 0;
  }
  tFemtoTrack->SetKinkIndexes(indexes);
}

void AliFemtoEventReaderAOD::SetFilterBit(UInt_t ibit)
{
  fFilterBit = (1 << (ibit));
}






