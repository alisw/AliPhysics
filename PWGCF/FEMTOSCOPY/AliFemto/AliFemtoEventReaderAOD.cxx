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
#include "AliAODMCHeader.h"
#include "AliESDtrack.h"

#include "AliFmPhysicalHelixD.h"
#include "AliFmThreeVectorF.h"

#include "SystemOfUnits.h"

#include "AliFemtoEvent.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoModelGlobalHiddenInfo.h"
#include "AliPID.h"

#include "AliAODpidUtil.h"

ClassImp(AliFemtoEventReaderAOD)

#if !(ST_NO_NAMESPACES)
  using namespace units;
#endif

using namespace std;

double fV1[3];

//____________________________
//constructor with 0 parameters , look at default settings 
AliFemtoEventReaderAOD::AliFemtoEventReaderAOD():
  fNumberofEvent(0),
  fCurEvent(0),
  fEvent(0x0),
  fAllTrue(160),
  fAllFalse(160),
  fFilterBit(0),
  fFilterMask(0),
  //  fPWG2AODTracks(0x0),
  fReadMC(0),
  fReadV0(0),
  fUsePreCent(0),
  fEstEventMult(kCentrality),
  fAODpidUtil(0),
  fAODheader(0),
  fInputFile(" "),
  fFileName(" "),
  fTree(0x0),
  fAodFile(0x0),
    fMagFieldSign(1),
    fisEPVZ(kTRUE)
{
  // default constructor
  fAllTrue.ResetAllBits(kTRUE);
  fAllFalse.ResetAllBits(kFALSE);
  fCentRange[0] = 0;
  fCentRange[1] = 1000;
}

AliFemtoEventReaderAOD::AliFemtoEventReaderAOD(const AliFemtoEventReaderAOD &aReader) :
  AliFemtoEventReader(),
  fNumberofEvent(0),
  fCurEvent(0),
  fEvent(0x0),
  fAllTrue(160),
  fAllFalse(160),
  fFilterBit(0),
  fFilterMask(0),
  //  fPWG2AODTracks(0x0),
  fReadMC(0),
  fReadV0(0),
  fUsePreCent(0),
  fEstEventMult(kCentrality),
  fAODpidUtil(0),
  fAODheader(0),
  fInputFile(" "),
  fFileName(" "),
  fTree(0x0),
  fAodFile(0x0),
    fMagFieldSign(1),
    fisEPVZ(kTRUE)
{
  // copy constructor
  fReadMC = aReader.fReadMC;
  fReadV0 = aReader.fReadV0;
  fInputFile = aReader.fInputFile;
  fFileName  = aReader.fFileName;
  fNumberofEvent = aReader.fNumberofEvent;
  fCurEvent = aReader.fCurEvent;
  fEvent = new AliAODEvent();
  fAodFile = new TFile(aReader.fAodFile->GetName());
  fAllTrue.ResetAllBits(kTRUE);
  fAllFalse.ResetAllBits(kFALSE);
  fFilterBit = aReader.fFilterBit;
  //  fPWG2AODTracks = aReader.fPWG2AODTracks;
  fAODpidUtil = aReader.fAODpidUtil;
  fAODheader = aReader.fAODheader;
  fCentRange[0] = aReader.fCentRange[0];
  fCentRange[1] = aReader.fCentRange[1];
  fEstEventMult = aReader.fEstEventMult;
  fUsePreCent = aReader.fUsePreCent;
}
//__________________
//Destructor
AliFemtoEventReaderAOD::~AliFemtoEventReaderAOD()
{
  // destructor
  delete fTree;
  delete fEvent;
  delete fAodFile;
//   if (fPWG2AODTracks) {
//     fPWG2AODTracks->Delete();
//     delete fPWG2AODTracks;
//   }
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
  fFilterMask = aReader.fFilterMask;
  //  fPWG2AODTracks = aReader.fPWG2AODTracks;
  fAODpidUtil = aReader.fAODpidUtil;
  fAODheader = aReader.fAODheader;
  fCentRange[0] = aReader.fCentRange[0];
  fCentRange[1] = aReader.fCentRange[1];
  fUsePreCent = aReader.fUsePreCent;
  fEstEventMult = aReader.fEstEventMult;

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
		  // 		  cout<<"putting file  "<<string(buffer)<<" into analysis"<<endl;
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
    // cout<<"reader"<<endl;
  if (fCurEvent==fNumberofEvent)//open next file  
    {
      if(fNumberofEvent==0)	
	{
	  fEvent=new AliAODEvent();
	  fEvent->ReadFromTree(fTree);

	  // Check for the existence of the additional information
// 	  fPWG2AODTracks = (TClonesArray *) fEvent->GetList()->FindObject("pwg2aodtracks");

// 	  if (fPWG2AODTracks) {
// 	    cout << "Found additional PWG2 specific information in the AOD!" << endl;
// 	    cout << "Reading only tracks with the additional information" << endl;
// 	  }

	  fNumberofEvent=fTree->GetEntries();
	  //	  cout<<"Number of Entries in file "<<fNumberofEvent<<endl;
	  fCurEvent=0;
	}
      else //no more data to read
	{
            // cout<<"no more files "<<hbtEvent<<endl;
	  fReaderStatus=1;
	  return hbtEvent; 
	}
    }		

    // cout<<"starting to read event "<<fCurEvent<<endl;
  fTree->GetEvent(fCurEvent);//getting next event
  //  cout << "Read event " << fEvent << " from file " << fTree << endl;
	
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
  tEvent->SetZDCParticipants(0);
  tEvent->SetTriggerMask(fEvent->GetTriggerMask());
  tEvent->SetTriggerCluster(fEvent->GetTriggerCluster());
  
  // Attempt to access MC header
  AliAODMCHeader *mcH;
  TClonesArray *mcP=0;
  if (fReadMC) {
    mcH = (AliAODMCHeader *) fEvent->FindListObject(AliAODMCHeader::StdBranchName());
    if (!mcH) {
      cout << "AOD MC information requested, but no header found!" << endl;
    }

    mcP = (TClonesArray *) fEvent->FindListObject(AliAODMCParticle::StdBranchName());
    if (!mcP) {
      cout << "AOD MC information requested, but no particle array found!" << endl;
    }
  }

  tEvent->SetReactionPlaneAngle(fEvent->GetHeader()->GetQTheta(0)/2.0);

  Int_t *motherids=0;
  if (mcP) {
    motherids = new Int_t[((AliAODMCParticle *) mcP->At(mcP->GetEntries()-1))->GetLabel()];
    for (int ip=0; ip<mcP->GetEntries(); ip++) motherids[ip] = 0;

    // Read in mother ids
    AliAODMCParticle *motherpart;
    for (int ip=0; ip<mcP->GetEntries(); ip++) {
      motherpart = (AliAODMCParticle *) mcP->At(ip);
      if (motherpart->GetDaughter(0) > 0)
	motherids[motherpart->GetDaughter(0)] = ip;
      if (motherpart->GetDaughter(1) > 0)
	motherids[motherpart->GetDaughter(1)] = ip;
    }
  }

  // Primary Vertex position
  //  double fV1[3];
  fEvent->GetPrimaryVertex()->GetPosition(fV1);

  AliFmThreeVectorF vertex(fV1[0],fV1[1],fV1[2]);
  tEvent->SetPrimVertPos(vertex);
	
  //starting to reading tracks
  int nofTracks=0;  //number of reconstructed tracks in event

  // Check to see whether the additional info exists
  //   if (fPWG2AODTracks)
  //     nofTracks=fPWG2AODTracks->GetEntries();
  //   else
  nofTracks=fEvent->GetNumberOfTracks();
  
  AliEventplane *ep = fEvent->GetEventplane();
  if (ep) {
    tEvent->SetEP(ep);
        if (fisEPVZ)
            tEvent->SetReactionPlaneAngle(ep->GetEventplane("V0",fEvent,2));
        else
    tEvent->SetReactionPlaneAngle(ep->GetEventplane("Q"));
  }

  AliCentrality *cent = fEvent->GetCentrality();
  
  if (!fEstEventMult && cent && fUsePreCent) {
    if ((cent->GetCentralityPercentile("V0M")*10 < fCentRange[0]) ||
	(cent->GetCentralityPercentile("V0M")*10 > fCentRange[1]))
      {
            // cout << "Centrality " << cent->GetCentralityPercentile("V0M") << " outside of preselection range " << fCentRange[0] << " - " << fCentRange[1] << endl;
	
	return;
      }
  }

  int realnofTracks=0;   // number of track which we use in a analysis
  int tracksPrim=0;     

  int labels[20000];	
  for (int il=0; il<20000; il++) labels[il] = -1;

  // looking for global tracks and saving their numbers to copy from them PID information to TPC-only tracks in the main loop over tracks
  for (int i=0;i<nofTracks;i++) {
    const AliAODTrack *aodtrack=fEvent->GetTrack(i);
    if (!aodtrack->TestFilterBit(fFilterBit)) {
      if(aodtrack->GetID() < 0) continue;
      labels[aodtrack->GetID()] = i;
    }
  }

  int tNormMult = 0;
  for (int i=0;i<nofTracks;i++)
    {
      AliFemtoTrack* trackCopy = new AliFemtoTrack();	

//       if (fPWG2AODTracks) {
// 	// Read tracks from the additional pwg2 specific AOD part
// 	// if they exist
// 	// Note that in that case all the AOD tracks without the 
// 	// additional information will be ignored !
// 	AliPWG2AODTrack *pwg2aodtrack = (AliPWG2AODTrack *) fPWG2AODTracks->At(i);

// 	// Getting the AOD track through the ref of the additional info
// 	AliAODTrack *aodtrack = pwg2aodtrack->GetRefAODTrack();	
// 	if (!aodtrack->TestFilterBit(fFilterBit)) {
// 	  delete trackCopy;
// 	  continue;
// 	}

 
// 	if (aodtrack->IsOn(AliESDtrack::kTPCrefit))
// 	  if (aodtrack->Chi2perNDF() < 6.0) 
// 	    if (aodtrack->Eta() < 0.9)
// 	      tNormMult++;


// 	CopyAODtoFemtoTrack(aodtrack, trackCopy, pwg2aodtrack);
	
// 	if (mcP) {
// 	  // Fill the hidden information with the simulated data
// 	  //	  Int_t pLabel = aodtrack->GetLabel();
// 	  AliAODMCParticle *tPart = GetParticleWithLabel(mcP, (TMath::Abs(aodtrack->GetLabel())));

// 	  // Check the mother information
	  
// 	  // Using the new way of storing the freeze-out information
// 	  // Final state particle is stored twice on the stack
// 	  // one copy (mother) is stored with original freeze-out information
// 	  //   and is not tracked
// 	  // the other one (daughter) is stored with primary vertex position
// 	  //   and is tracked
	  
// 	  // Freeze-out coordinates
// 	  double fpx=0.0, fpy=0.0, fpz=0.0, fpt=0.0;
// 	  fpx = tPart->Xv() - fV1[0];
// 	  fpy = tPart->Yv() - fV1[1];
// 	  fpz = tPart->Zv() - fV1[2];
// 	  fpt = tPart->T();

// 	  AliFemtoModelGlobalHiddenInfo *tInfo = new AliFemtoModelGlobalHiddenInfo();
// 	  tInfo->SetGlobalEmissionPoint(fpx, fpy, fpz);

// 	  fpx *= 1e13;
// 	  fpy *= 1e13;
// 	  fpz *= 1e13;
// 	  fpt *= 1e13;
	  
// 	  //      cout << "Looking for mother ids " << endl;
// 	  if (motherids[TMath::Abs(aodtrack->GetLabel())]>0) {
// 	    //	cout << "Got mother id" << endl;
// 	    AliAODMCParticle *mother = GetParticleWithLabel(mcP, motherids[TMath::Abs(aodtrack->GetLabel())]);
// 	    // Check if this is the same particle stored twice on the stack
// 	    if ((mother->GetPdgCode() == tPart->GetPdgCode() || (mother->Px() == tPart->Px()))) {
// 	      // It is the same particle
// 	      // Read in the original freeze-out information
// 	      // and convert it from to [fm]
	      
// 	      // EPOS style 
// 	      // 	  fpx = mother->Xv()*1e13*0.197327;
// 	      // 	  fpy = mother->Yv()*1e13*0.197327;
// 	      // 	  fpz = mother->Zv()*1e13*0.197327;
// 	      // 	  fpt = mother->T() *1e13*0.197327*0.5;
	      
	      
// 	      // Therminator style 
// 	      fpx = mother->Xv()*1e13;
// 	      fpy = mother->Yv()*1e13;
// 	      fpz = mother->Zv()*1e13;
// 	      fpt = mother->T() *1e13*3e10;
	      
// 	    }
// 	  }
	  
// 	  //       if (fRotateToEventPlane) {
// 	  // 	double tPhi = TMath::ATan2(fpy, fpx);
// 	  // 	double tRad = TMath::Hypot(fpx, fpy);
	
// 	  // 	fpx = tRad*TMath::Cos(tPhi - tReactionPlane);
// 	  // 	fpy = tRad*TMath::Sin(tPhi - tReactionPlane);
// 	  //       }

// 	  tInfo->SetPDGPid(tPart->GetPdgCode());

// 	  // 	  if (fRotateToEventPlane) {
// 	  // 	    double tPhi = TMath::ATan2(tPart->Py(), tPart->Px());
// 	  // 	    double tRad = TMath::Hypot(tPart->Px(), tPart->Py());
	    
// 	  // 	    tInfo->SetTrueMomentum(tRad*TMath::Cos(tPhi - tReactionPlane),
// 	  // 				   tRad*TMath::Sin(tPhi - tReactionPlane),
// 	  // 				   tPart->Pz());
// 	  // 	  }
// 	  //       else
// 	  tInfo->SetTrueMomentum(tPart->Px(), tPart->Py(), tPart->Pz());
// 	  Double_t mass2 = (tPart->E() *tPart->E() -
// 			    tPart->Px()*tPart->Px() -
// 			    tPart->Py()*tPart->Py() -
// 			    tPart->Pz()*tPart->Pz());
// 	  if (mass2>0.0)
// 	    tInfo->SetMass(TMath::Sqrt(mass2));
// 	  else 
// 	    tInfo->SetMass(0.0);
	  
// 	  tInfo->SetEmissionPoint(fpx, fpy, fpz, fpt);
// 	  trackCopy->SetHiddenInfo(tInfo);

// 	}

// 	double pxyz[3];
// 	aodtrack->PxPyPz(pxyz);//reading noconstarined momentum
// 	const AliFmThreeVectorD ktP(pxyz[0],pxyz[1],pxyz[2]);
// 	// Check the sanity of the tracks - reject zero momentum tracks
// 	if (ktP.Mag() == 0) {
// 	  delete trackCopy;
// 	  continue;
// 	}
//       }
//       else {
	// No additional information exists
	// Read in the normal AliAODTracks 

	//	const AliAODTrack *aodtrack=fEvent->GetTrack(i); // getting the AODtrack directly
	AliAODTrack *aodtrack=fEvent->GetTrack(i); // getting the AODtrack directly
	


	if (aodtrack->IsPrimaryCandidate()) tracksPrim++;
	
	if (fFilterBit && !aodtrack->TestFilterBit(fFilterBit)) {
	  delete trackCopy;
	  continue;
	}

	if (fFilterMask && !aodtrack->TestFilterBit(fFilterMask)) {
	  delete trackCopy;
	  continue;
	}		

	//counting particles to set multiplicity
	double impact[2];
	double covimpact[3];
	if (aodtrack->PropagateToDCA(fEvent->GetPrimaryVertex(),fEvent->GetMagneticField(),10000,impact,covimpact)) {
	  if(impact[0]<0.2 && TMath::Abs(impact[1]+fV1[2])<2.0)
	    //if (aodtrack->IsPrimaryCandidate()) //? instead of kinks?
	      if (aodtrack->Chi2perNDF() < 4.0) 
		if (aodtrack->Pt() > 0.15 && aodtrack->Pt() < 20) 
		  if (aodtrack->GetTPCNcls() > 70)
		    if (aodtrack->Eta() < 0.8)
		      tNormMult++;
	} 

	CopyAODtoFemtoTrack(aodtrack, trackCopy);

	// copying PID information from the correspondent track
	//	const AliAODTrack *aodtrackpid = fEvent->GetTrack(labels[-1-fEvent->GetTrack(i)->GetID()]);


	AliAODTrack *aodtrackpid;
	if((fFilterBit ==  (1 << (7))) || fFilterMask==128) //for TPC Only tracks we have to copy PID information from corresponding global tracks
	  aodtrackpid = fEvent->GetTrack(labels[-1-fEvent->GetTrack(i)->GetID()]);
	else
	  aodtrackpid = fEvent->GetTrack(i);
        CopyPIDtoFemtoTrack(aodtrackpid, trackCopy);
	
	if (mcP) {
	  // Fill the hidden information with the simulated data
	  //	  Int_t pLabel = aodtrack->GetLabel();
	  AliAODMCParticle *tPart = GetParticleWithLabel(mcP, (TMath::Abs(aodtrack->GetLabel())));
	  
	  AliFemtoModelGlobalHiddenInfo *tInfo = new AliFemtoModelGlobalHiddenInfo();
	  double fpx=0.0, fpy=0.0, fpz=0.0, fpt=0.0;
	  if (!tPart) {
	    fpx = fV1[0];
	    fpy = fV1[1];
	    fpz = fV1[2];
	    tInfo->SetGlobalEmissionPoint(fpx, fpy, fpz);
	    tInfo->SetPDGPid(0);
	    tInfo->SetTrueMomentum(0.0, 0.0, 0.0);
	    tInfo->SetEmissionPoint(0.0, 0.0, 0.0, 0.0);
	    tInfo->SetMass(0);
	  }
	  else {
	    // Check the mother information
	  
	    // Using the new way of storing the freeze-out information
	    // Final state particle is stored twice on the stack
	    // one copy (mother) is stored with original freeze-out information
	    //   and is not tracked
	    // the other one (daughter) is stored with primary vertex position
	    //   and is tracked
	    
	    // Freeze-out coordinates
	    fpx = tPart->Xv() - fV1[0];
	    fpy = tPart->Yv() - fV1[1];
	    fpz = tPart->Zv() - fV1[2];
	    //	  fpt = tPart->T();
	    
	    tInfo->SetGlobalEmissionPoint(fpx, fpy, fpz);
	    
	    fpx *= 1e13;
	    fpy *= 1e13;
	    fpz *= 1e13;
	    //	  fpt *= 1e13;
	    
	    //      cout << "Looking for mother ids " << endl;
	    if (motherids[TMath::Abs(aodtrack->GetLabel())]>0) {
	      //	cout << "Got mother id" << endl;
	      AliAODMCParticle *mother = GetParticleWithLabel(mcP, motherids[TMath::Abs(aodtrack->GetLabel())]);
	      // Check if this is the same particle stored twice on the stack
	      if (mother) {
		if ((mother->GetPdgCode() == tPart->GetPdgCode() || (mother->Px() == tPart->Px()))) {
		  // It is the same particle
		  // Read in the original freeze-out information
		  // and convert it from to [fm]
		  
		  // EPOS style 
		  // 	  fpx = mother->Xv()*1e13*0.197327;
		  // 	  fpy = mother->Yv()*1e13*0.197327;
		  // 	  fpz = mother->Zv()*1e13*0.197327;
		  // 	  fpt = mother->T() *1e13*0.197327*0.5;
		  
		  
		  // Therminator style 
		  fpx = mother->Xv()*1e13;
		  fpy = mother->Yv()*1e13;
		  fpz = mother->Zv()*1e13;
		  //	      fpt = mother->T() *1e13*3e10;
		  
		}
	      }
	    }
	    
	    //       if (fRotateToEventPlane) {
	    // 	double tPhi = TMath::ATan2(fpy, fpx);
	    // 	double tRad = TMath::Hypot(fpx, fpy);
	    
	    // 	fpx = tRad*TMath::Cos(tPhi - tReactionPlane);
	    // 	fpy = tRad*TMath::Sin(tPhi - tReactionPlane);
	    //       }
	    
	    tInfo->SetPDGPid(tPart->GetPdgCode());
	    
	    // 	  if (fRotateToEventPlane) {
	    // 	    double tPhi = TMath::ATan2(tPart->Py(), tPart->Px());
	    // 	    double tRad = TMath::Hypot(tPart->Px(), tPart->Py());
	    
	    // 	    tInfo->SetTrueMomentum(tRad*TMath::Cos(tPhi - tReactionPlane),
	    // 				   tRad*TMath::Sin(tPhi - tReactionPlane),
	    // 				   tPart->Pz());
	    // 	  }
	    //       else
	    tInfo->SetTrueMomentum(tPart->Px(), tPart->Py(), tPart->Pz());
	    Double_t mass2 = (tPart->E() *tPart->E() -
			      tPart->Px()*tPart->Px() -
			      tPart->Py()*tPart->Py() -
			      tPart->Pz()*tPart->Pz());
	    if (mass2>0.0)
	      tInfo->SetMass(TMath::Sqrt(mass2));
	    else 
	      tInfo->SetMass(0.0);
	    
	    tInfo->SetEmissionPoint(fpx, fpy, fpz, fpt);
	  }
	  trackCopy->SetHiddenInfo(tInfo);
	}

	double pxyz[3];

	//AliExternalTrackParam *param = new AliExternalTrackParam(*aodtrack->GetInnerParam());
	trackCopy->SetInnerMomentum(aodtrack->GetTPCmomentum());

	aodtrack->PxPyPz(pxyz);//reading noconstarined momentum
	const AliFmThreeVectorD ktP(pxyz[0],pxyz[1],pxyz[2]);
	// Check the sanity of the tracks - reject zero momentum tracks
	if (ktP.Mag() == 0) {
	  delete trackCopy;
	  continue;
	}
	//    }
  
	
	tEvent->TrackCollection()->push_back(trackCopy);//adding track to analysis
	realnofTracks++;//real number of tracks		
    }
  
  tEvent->SetNumberOfTracks(realnofTracks);//setting number of track which we read in event	
  tEvent->SetNormalizedMult(tracksPrim);

  if (fEstEventMult==kCentrality) {
    //  AliCentrality *cent = fEvent->GetCentrality();
    //cout<<"AliFemtoEventReaderAOD:"<<lrint(10*cent->GetCentralityPercentile("V0M"))<<endl;
    if (cent) tEvent->SetNormalizedMult(lrint(10*cent->GetCentralityPercentile("V0M")));
    //  if (cent) tEvent->SetNormalizedMult((int) cent->GetCentralityPercentile("V0M"));
    
    if (cent) {
      tEvent->SetCentralityV0(cent->GetCentralityPercentile("V0M"));
      //    tEvent->SetCentralityFMD(cent->GetCentralityPercentile("FMD"));
      tEvent->SetCentralitySPD1(cent->GetCentralityPercentile("CL1"));
      //    tEvent->SetCentralityTrk(cent->GetCentralityPercentile("TRK"));
    }
  }
  else if(fEstEventMult==kGlobalCount){
    tEvent->SetNormalizedMult(tNormMult); //particles counted in the loop, trying to reproduce GetReferenceMultiplicity. If better (default) method appears it should be changed
  }
  else if(fEstEventMult==kReference)
    {
      tEvent->SetNormalizedMult(fAODheader->GetRefMultiplicity());
    }
  else if(fEstEventMult==kTPCOnlyRef)
    {
      tEvent->SetNormalizedMult(fAODheader->GetTPConlyRefMultiplicity());
    }
  else if(fEstEventMult == kVZERO)
    {
      Float_t multV0 = 0;
      for (Int_t i=0; i<64; i++)
	multV0 += fEvent->GetVZEROData()->GetMultiplicity(i);
      tEvent->SetNormalizedMult(multV0);
    }

  if (mcP) delete [] motherids;

    // cout<<"end of reading nt "<<nofTracks<<" real number "<<realnofTracks<<endl;

  if(fReadV0)
    {
      int count_pass = 0;
      for (Int_t i = 0; i < fEvent->GetNumberOfV0s(); i++) {
	AliAODv0* aodv0 = fEvent->GetV0(i);
	if (!aodv0) continue;
	if(aodv0->GetNDaughters()>2) continue;
	if(aodv0->GetNProngs()>2) continue;
	if(aodv0->GetCharge()!=0) continue;
	if(aodv0->ChargeProng(0)==aodv0->ChargeProng(1)) continue;
	if(aodv0->CosPointingAngle(fV1)<0.998) continue;
	AliFemtoV0* trackCopyV0 = new AliFemtoV0();
	count_pass++;
	CopyAODtoFemtoV0(aodv0, trackCopyV0);
	tEvent->V0Collection()->push_back(trackCopyV0);
	//cout<<"Pushback v0 to v0collection"<<endl;
      }
    }

}

void AliFemtoEventReaderAOD::CopyAODtoFemtoTrack(AliAODTrack *tAodTrack, 
						 AliFemtoTrack *tFemtoTrack 
						 //						 AliPWG2AODTrack *tPWG2AODTrack
						 )
{
  // Copy the track information from the AOD into the internal AliFemtoTrack
  // If it exists, use the additional information from the PWG2 AOD

  // Primary Vertex position
  
  fEvent->GetPrimaryVertex()->GetPosition(fV1);
  //  fEvent->GetPrimaryVertex()->GetXYZ(fV1);

  tFemtoTrack->SetCharge(tAodTrack->Charge());
  
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
  tFemtoTrack->SetFlags(tAodTrack->GetFlags());
  tFemtoTrack->SetLabel(tAodTrack->GetLabel());
		
  // Track quality information 
  float covmat[6];
  tAodTrack->GetCovMatrix(covmat);  

        // ! DCA information is done in CopyPIDtoFemtoTrack()
  
	// double impact[2];
	// double covimpact[3];

	// if (!tAodTrack->PropagateToDCA(fEvent->GetPrimaryVertex(),fEvent->GetMagneticField(),10000,impact,covimpact)) {
	//   //cout << "sth went wrong with dca propagation" << endl;
	//   tFemtoTrack->SetImpactD(-1000.0);
	//   tFemtoTrack->SetImpactZ(-1000.0);

	// }
	// else {
	//   tFemtoTrack->SetImpactD(impact[0]);
	//   tFemtoTrack->SetImpactZ(impact[1]+fV1[2]);
	// }

  //   if (TMath::Abs(tAodTrack->Xv()) > 0.00000000001)
  //     tFemtoTrack->SetImpactD(TMath::Hypot(tAodTrack->Xv(), tAodTrack->Yv())*(tAodTrack->Xv()/TMath::Abs(tAodTrack->Xv())));
  //   else
  //     tFemtoTrack->SetImpactD(0.0);
  //   tFemtoTrack->SetImpactD(tAodTrack->DCA());
    
  //   tFemtoTrack->SetImpactZ(tAodTrack->ZAtDCA());


  //   tFemtoTrack->SetImpactD(TMath::Hypot(tAodTrack->Xv() - fV1[0], tAodTrack->Yv() - fV1[1]));
  //   tFemtoTrack->SetImpactZ(tAodTrack->Zv() - fV1[2]);


  //   cout 
    //    << "dca" << TMath::Hypot(tAodTrack->Xv() - fV1[0], tAodTrack->Yv() - fV1[1]) 
    //    << "xv - fv10 = "<< tAodTrack->Xv() - fV1[0] 
    //    << tAodTrack->Yv() - fV1[1] 
//     << "xv = " << tAodTrack->Xv() << endl 
//     << "fv1[0] = " << fV1[0]  << endl 
//     << "yv = " << tAodTrack->Yv()  << endl 
//     << "fv1[1] = " << fV1[1]  << endl 
//     << "zv = " << tAodTrack->Zv()  << endl 
//     << "fv1[2] = " << fV1[2]  << endl 
//     << "impact[0] = " << impact[0]  << endl 
//     << "impact[1] = " << impact[1]  << endl 
//     << endl << endl ;

  tFemtoTrack->SetCdd(covmat[0]);
  tFemtoTrack->SetCdz(covmat[1]);
  tFemtoTrack->SetCzz(covmat[2]);
  tFemtoTrack->SetITSchi2(tAodTrack->Chi2perNDF());
  tFemtoTrack->SetITSncls(tAodTrack->GetITSNcls());
  tFemtoTrack->SetTPCchi2(tAodTrack->Chi2perNDF());
  tFemtoTrack->SetTPCncls(tAodTrack->GetTPCNcls());
  tFemtoTrack->SetTPCnclsF(tAodTrack->GetTPCNcls());
  tFemtoTrack->SetTPCsignalN(1); 
  tFemtoTrack->SetTPCsignalS(1); 
  tFemtoTrack->SetTPCsignal(tAodTrack->GetTPCsignal());

//   if (tPWG2AODTrack) {
//     // Copy the PWG2 specific information if it exists
//     tFemtoTrack->SetTPCClusterMap(tPWG2AODTrack->GetTPCClusterMap());
//     tFemtoTrack->SetTPCSharedMap(tPWG2AODTrack->GetTPCSharedMap());
    
//     double xtpc[3] = {0,0,0};
//     tPWG2AODTrack->GetTPCNominalEntrancePoint(xtpc);
//     tFemtoTrack->SetNominalTPCEntrancePoint(xtpc);
//     tPWG2AODTrack->GetTPCNominalExitPoint(xtpc);
//     tFemtoTrack->SetNominalTPCExitPoint(xtpc);
//   }
//   else {
    // If not use dummy values
  tFemtoTrack->SetTPCClusterMap(tAodTrack->GetTPCClusterMap());
  tFemtoTrack->SetTPCSharedMap(tAodTrack->GetTPCSharedMap());
  

  float globalPositionsAtRadii[9][3];
  float bfield = 5*fMagFieldSign;
  GetGlobalPositionAtGlobalRadiiThroughTPC(tAodTrack,bfield,globalPositionsAtRadii);
  double tpcEntrance[3]={globalPositionsAtRadii[0][0],globalPositionsAtRadii[0][1],globalPositionsAtRadii[0][2]};
  double **tpcPositions;
  tpcPositions = new double*[9];
  for(int i=0;i<9;i++)
    tpcPositions[i] = new double[3];
  double tpcExit[3]={globalPositionsAtRadii[8][0],globalPositionsAtRadii[8][1],globalPositionsAtRadii[8][2]};
  for(int i=0;i<9;i++)
    {
      tpcPositions[i][0] = globalPositionsAtRadii[i][0];
      tpcPositions[i][1] = globalPositionsAtRadii[i][1];
      tpcPositions[i][2] = globalPositionsAtRadii[i][2];
    }
  tFemtoTrack->SetNominalTPCEntrancePoint(tpcEntrance);
  tFemtoTrack->SetNominalTPCPoints(tpcPositions);
  tFemtoTrack->SetNominalTPCExitPoint(tpcExit);

  //   }
  
  //   //  cout << "Track has " << TMath::Hypot(tAodTrack->Xv(), tAodTrack->Yv()) << "  " << tAodTrack->Zv() << "  " << tAodTrack->GetTPCNcls() << endl;
  
  
  int indexes[3];
  for (int ik=0; ik<3; ik++) {
    indexes[ik] = 0;
  }
  tFemtoTrack->SetKinkIndexes(indexes);


  for (int ii=0; ii<6; ii++){
    tFemtoTrack->SetITSHitOnLayer(ii,tAodTrack->HasPointOnITSLayer(ii));
  }


}

void AliFemtoEventReaderAOD::CopyAODtoFemtoV0(AliAODv0 *tAODv0, AliFemtoV0 *tFemtoV0)
{
  tFemtoV0->SetdecayLengthV0(tAODv0->DecayLength(fV1));
  tFemtoV0->SetdecayVertexV0X(tAODv0->DecayVertexV0X());
  tFemtoV0->SetdecayVertexV0Y(tAODv0->DecayVertexV0Y());
  tFemtoV0->SetdecayVertexV0Z(tAODv0->DecayVertexV0Z());
  AliFemtoThreeVector decayvertex(tAODv0->DecayVertexV0X(),tAODv0->DecayVertexV0Y(),tAODv0->DecayVertexV0Z());
  tFemtoV0->SetdecayVertexV0(decayvertex);
  tFemtoV0->SetdcaV0Daughters(tAODv0->DcaV0Daughters());
  tFemtoV0->SetdcaV0ToPrimVertex(tAODv0->DcaV0ToPrimVertex());
  tFemtoV0->SetdcaPosToPrimVertex(tAODv0->DcaPosToPrimVertex());
  tFemtoV0->SetdcaNegToPrimVertex(tAODv0->DcaNegToPrimVertex());
  tFemtoV0->SetmomPosX(tAODv0->MomPosX());
  tFemtoV0->SetmomPosY(tAODv0->MomPosY());
  tFemtoV0->SetmomPosZ(tAODv0->MomPosZ());
  AliFemtoThreeVector mompos(tAODv0->MomPosX(),tAODv0->MomPosY(),tAODv0->MomPosZ());
  tFemtoV0->SetmomPos(mompos);
  tFemtoV0->SetmomNegX(tAODv0->MomNegX());
  tFemtoV0->SetmomNegY(tAODv0->MomNegY());
  tFemtoV0->SetmomNegZ(tAODv0->MomNegZ());
  AliFemtoThreeVector momneg(tAODv0->MomNegX(),tAODv0->MomNegY(),tAODv0->MomNegZ());
  tFemtoV0->SetmomNeg(momneg);

  //jest cos takiego w AliFemtoV0.h czego nie ma w AliAODv0.h
  //void SettpcHitsPos(const int& i);      
  //void SettpcHitsNeg(const int& i);      

  //void SetTrackTopologyMapPos(unsigned int word, const unsigned long& m);
  //void SetTrackTopologyMapNeg(unsigned int word, const unsigned long& m);

  tFemtoV0->SetmomV0X(tAODv0->MomV0X());
  tFemtoV0->SetmomV0Y(tAODv0->MomV0Y());
  tFemtoV0->SetmomV0Z(tAODv0->MomV0Z());
  AliFemtoThreeVector momv0(tAODv0->MomV0X(),tAODv0->MomV0Y(),tAODv0->MomV0Z());
  tFemtoV0->SetmomV0(momv0);
  tFemtoV0->SetalphaV0(tAODv0->AlphaV0());
  tFemtoV0->SetptArmV0(tAODv0->PtArmV0());
  tFemtoV0->SeteLambda(tAODv0->ELambda());
  tFemtoV0->SeteK0Short(tAODv0->EK0Short());
  tFemtoV0->SetePosProton(tAODv0->EPosProton());
  tFemtoV0->SeteNegProton(tAODv0->ENegProton());
  tFemtoV0->SetmassLambda(tAODv0->MassLambda());
  tFemtoV0->SetmassAntiLambda(tAODv0->MassAntiLambda());
  tFemtoV0->SetmassK0Short(tAODv0->MassK0Short());
  tFemtoV0->SetrapLambda(tAODv0->RapLambda());
  tFemtoV0->SetrapK0Short(tAODv0->RapK0Short());
  
  //void SetcTauLambda( float x);   
  //void SetcTauK0Short( float x); 
  
  //tFemtoV0->SetptV0(::sqrt(tAODv0->Pt2V0())); //!
  tFemtoV0->SetptV0(tAODv0->Pt());
  tFemtoV0->SetptotV0(::sqrt(tAODv0->Ptot2V0()));
  //tFemtoV0->SetptPos(::sqrt(tAODv0->MomPosX()*tAODv0->MomPosX()+tAODv0->MomPosY()*tAODv0->MomPosY()));
  //tFemtoV0->SetptotPos(::sqrt(tAODv0->Ptot2Pos()));
  //tFemtoV0->SetptNeg(::sqrt(tAODv0->MomNegX()*tAODv0->MomNegX()+tAODv0->MomNegY()*tAODv0->MomNegY()));
  //tFemtoV0->SetptotNeg(::sqrt(tAODv0->Ptot2Neg()));
  
  tFemtoV0->SetidNeg(tAODv0->GetNegID());
  //cout<<"tAODv0->GetNegID(): "<<tAODv0->GetNegID()<<endl;
  //cout<<"tFemtoV0->IdNeg(): "<<tFemtoV0->IdNeg()<<endl;
  tFemtoV0->SetidPos(tAODv0->GetPosID());

  tFemtoV0->SetEtaV0(tAODv0->Eta());
  tFemtoV0->SetPhiV0(tAODv0->Phi());
  tFemtoV0->SetCosPointingAngle(tAODv0->CosPointingAngle(fV1));
  //tFemtoV0->SetYV0(tAODv0->Y());

  //void SetdedxNeg(float x);
  //void SeterrdedxNeg(float x);//Gael 04Fev2002
  //void SetlendedxNeg(float x);//Gael 04Fev2002
  //void SetdedxPos(float x);
  //void SeterrdedxPos(float x);//Gael 04Fev2002
  //void SetlendedxPos(float x);//Gael 04Fev2002

  //tFemtoV0->SetEtaPos(tAODv0->PseudoRapPos());
  //tFemtoV0->SetEtaNeg(tAODv0->PseudoRapNeg());

  AliAODTrack *trackpos = (AliAODTrack*)tAODv0->GetDaughter(0);
  AliAODTrack *trackneg = (AliAODTrack*)tAODv0->GetDaughter(1);

  if(trackpos && trackneg)
    {
      tFemtoV0->SetEtaPos(trackpos->Eta());
      tFemtoV0->SetEtaNeg(trackneg->Eta());
      tFemtoV0->SetptotPos(tAODv0->PProng(0));
      tFemtoV0->SetptotNeg(tAODv0->PProng(1));
      tFemtoV0->SetptPos(trackpos->Pt());
      tFemtoV0->SetptNeg(trackneg->Pt());

      //tFemtoV0->SetEtaPos(trackpos->Eta()); //tAODv0->PseudoRapPos()
      //tFemtoV0->SetEtaNeg(trackneg->Eta()); //tAODv0->PseudoRapNeg()
      tFemtoV0->SetTPCNclsPos(trackpos->GetTPCNcls());
      tFemtoV0->SetTPCNclsNeg(trackneg->GetTPCNcls());
      tFemtoV0->SetTPCclustersPos(trackpos->GetTPCClusterMap());
      tFemtoV0->SetTPCclustersNeg(trackneg->GetTPCClusterMap());
      tFemtoV0->SetTPCsharingPos(trackpos->GetTPCSharedMap());
      tFemtoV0->SetTPCsharingNeg(trackneg->GetTPCSharedMap());
      tFemtoV0->SetNdofPos(trackpos->Chi2perNDF());
      tFemtoV0->SetNdofNeg(trackneg->Chi2perNDF());
      tFemtoV0->SetStatusPos(trackpos->GetStatus());
      tFemtoV0->SetStatusNeg(trackneg->GetStatus());

      tFemtoV0->SetPosNSigmaTPCK(fAODpidUtil->NumberOfSigmasTPC(trackpos,AliPID::kKaon));
      tFemtoV0->SetNegNSigmaTPCK(fAODpidUtil->NumberOfSigmasTPC(trackneg,AliPID::kKaon));
      tFemtoV0->SetPosNSigmaTPCP(fAODpidUtil->NumberOfSigmasTPC(trackpos,AliPID::kProton));
      tFemtoV0->SetNegNSigmaTPCP(fAODpidUtil->NumberOfSigmasTPC(trackneg,AliPID::kProton));
      tFemtoV0->SetPosNSigmaTPCPi(fAODpidUtil->NumberOfSigmasTPC(trackpos,AliPID::kPion));
      tFemtoV0->SetNegNSigmaTPCPi(fAODpidUtil->NumberOfSigmasTPC(trackneg,AliPID::kPion));


      float bfield = 5*fMagFieldSign;
      float globalPositionsAtRadiiPos[9][3];
      GetGlobalPositionAtGlobalRadiiThroughTPC(trackpos,bfield,globalPositionsAtRadiiPos);
      double tpcEntrancePos[3]={globalPositionsAtRadiiPos[0][0],globalPositionsAtRadiiPos[0][1],globalPositionsAtRadiiPos[0][2]};
      double tpcExitPos[3]={globalPositionsAtRadiiPos[8][0],globalPositionsAtRadiiPos[8][1],globalPositionsAtRadiiPos[8][2]};

      float globalPositionsAtRadiiNeg[9][3];
      GetGlobalPositionAtGlobalRadiiThroughTPC(trackneg,bfield,globalPositionsAtRadiiNeg);
      double tpcEntranceNeg[3]={globalPositionsAtRadiiNeg[0][0],globalPositionsAtRadiiNeg[0][1],globalPositionsAtRadiiNeg[0][2]};
      double tpcExitNeg[3]={globalPositionsAtRadiiNeg[8][0],globalPositionsAtRadiiNeg[8][1],globalPositionsAtRadiiNeg[8][2]};

      AliFemtoThreeVector tmpVec;
      tmpVec.SetX(tpcEntrancePos[0]); tmpVec.SetX(tpcEntrancePos[1]); tmpVec.SetX(tpcEntrancePos[2]);
      tFemtoV0->SetNominalTpcEntrancePointPos(tmpVec);

      tmpVec.SetX(tpcExitPos[0]); tmpVec.SetX(tpcExitPos[1]); tmpVec.SetX(tpcExitPos[2]);
      tFemtoV0->SetNominalTpcExitPointPos(tmpVec);

      tmpVec.SetX(tpcEntranceNeg[0]); tmpVec.SetX(tpcEntranceNeg[1]); tmpVec.SetX(tpcEntranceNeg[2]);
      tFemtoV0->SetNominalTpcEntrancePointNeg(tmpVec);

      tmpVec.SetX(tpcExitNeg[0]); tmpVec.SetX(tpcExitNeg[1]); tmpVec.SetX(tpcExitNeg[2]);
      tFemtoV0->SetNominalTpcExitPointNeg(tmpVec);

      AliFemtoThreeVector vecTpcPos[9];
      AliFemtoThreeVector vecTpcNeg[9];
      for(int i=0;i<9;i++)
	{
	  vecTpcPos[i].SetX(globalPositionsAtRadiiPos[i][0]); vecTpcPos[i].SetY(globalPositionsAtRadiiPos[i][1]); vecTpcPos[i].SetZ(globalPositionsAtRadiiPos[i][2]);
	  vecTpcNeg[i].SetX(globalPositionsAtRadiiNeg[i][0]); vecTpcNeg[i].SetY(globalPositionsAtRadiiNeg[i][1]); vecTpcNeg[i].SetZ(globalPositionsAtRadiiNeg[i][2]);
	}
      tFemtoV0->SetNominalTpcPointPos(vecTpcPos);
      tFemtoV0->SetNominalTpcPointNeg(vecTpcNeg);

      tFemtoV0->SetTPCMomentumPos(trackpos->GetTPCmomentum());
      tFemtoV0->SetTPCMomentumNeg(trackneg->GetTPCmomentum());

      tFemtoV0->SetdedxPos(trackpos->GetTPCsignal());
      tFemtoV0->SetdedxNeg(trackneg->GetTPCsignal());

        if((tFemtoV0->StatusPos()&AliESDtrack::kTOFpid)==0 || (tFemtoV0->StatusPos()&AliESDtrack::kTIME)==0 || (tFemtoV0->StatusPos()&AliESDtrack::kTOFout)==0 )
	{
            if((tFemtoV0->StatusNeg()&AliESDtrack::kTOFpid)==0 || (tFemtoV0->StatusNeg()&AliESDtrack::kTIME)==0 || (tFemtoV0->StatusNeg()&AliESDtrack::kTOFout)==0 )
	    {
	      tFemtoV0->SetPosNSigmaTOFK(-1000);
	      tFemtoV0->SetNegNSigmaTOFK(-1000);
	      tFemtoV0->SetPosNSigmaTOFP(-1000);
	      tFemtoV0->SetNegNSigmaTOFP(-1000);
	      tFemtoV0->SetPosNSigmaTOFPi(-1000);
	      tFemtoV0->SetNegNSigmaTOFPi(-1000);

	      tFemtoV0->SetTOFProtonTimePos(-1000);
	      tFemtoV0->SetTOFPionTimePos(-1000);
	      tFemtoV0->SetTOFKaonTimePos(-1000);
	      tFemtoV0->SetTOFProtonTimeNeg(-1000);
	      tFemtoV0->SetTOFPionTimeNeg(-1000);
	      tFemtoV0->SetTOFKaonTimeNeg(-1000);
	    }
	}
      else
	{
	  tFemtoV0->SetPosNSigmaTOFK(fAODpidUtil->NumberOfSigmasTOF(trackpos,AliPID::kKaon));
	  tFemtoV0->SetNegNSigmaTOFK(fAODpidUtil->NumberOfSigmasTOF(trackneg,AliPID::kKaon));
	  tFemtoV0->SetPosNSigmaTOFP(fAODpidUtil->NumberOfSigmasTOF(trackpos,AliPID::kProton));
	  tFemtoV0->SetNegNSigmaTOFP(fAODpidUtil->NumberOfSigmasTOF(trackneg,AliPID::kProton));
	  tFemtoV0->SetPosNSigmaTOFPi(fAODpidUtil->NumberOfSigmasTOF(trackpos,AliPID::kPion));
	  tFemtoV0->SetNegNSigmaTOFPi(fAODpidUtil->NumberOfSigmasTOF(trackneg,AliPID::kPion));

	  double TOFSignalPos = trackpos->GetTOFsignal();
	  double TOFSignalNeg = trackneg->GetTOFsignal();
	  double pidPos[5];
	  double pidNeg[5];
	  trackpos->GetIntegratedTimes(pidPos);
	  trackneg->GetIntegratedTimes(pidNeg);

	  tFemtoV0->SetTOFPionTimePos(TOFSignalPos-pidPos[2]);
	  tFemtoV0->SetTOFKaonTimePos(TOFSignalPos-pidPos[3]);
	  tFemtoV0->SetTOFProtonTimePos(TOFSignalPos-pidPos[4]);
	  tFemtoV0->SetTOFPionTimeNeg(TOFSignalNeg-pidNeg[2]);
	  tFemtoV0->SetTOFKaonTimeNeg(TOFSignalNeg-pidNeg[3]);
	  tFemtoV0->SetTOFProtonTimeNeg(TOFSignalNeg-pidNeg[4]);
	}
    }
  else
    {
      tFemtoV0->SetStatusPos(999);
      tFemtoV0->SetStatusNeg(999);
    }
  tFemtoV0->SetOnFlyStatusV0(tAODv0->GetOnFlyStatus());
}

void AliFemtoEventReaderAOD::SetFilterBit(UInt_t ibit)
{
  fFilterBit = (1 << (ibit));
}


void AliFemtoEventReaderAOD::SetFilterMask(int ibit)
{
  fFilterMask = ibit;
}

void AliFemtoEventReaderAOD::SetReadMC(unsigned char a)
{
  fReadMC = a;
}


void AliFemtoEventReaderAOD::SetReadV0(unsigned char a)
{
  fReadV0 = a;
}

void AliFemtoEventReaderAOD::SetUseMultiplicity(EstEventMult aType)
{
  fEstEventMult = aType;
}

AliAODMCParticle* AliFemtoEventReaderAOD::GetParticleWithLabel(TClonesArray *mcP, Int_t aLabel)
{
  if (aLabel < 0) return 0;
  AliAODMCParticle *aodP;
  Int_t posstack = 0;
  if (aLabel > mcP->GetEntries())
    posstack = mcP->GetEntries();
  else
    posstack = aLabel;

  aodP = (AliAODMCParticle *) mcP->At(posstack);
  if (aodP->GetLabel() > posstack) {
    do {
      aodP = (AliAODMCParticle *) mcP->At(posstack);
      if (aodP->GetLabel() == aLabel) return aodP;
      posstack--;
    }
    while (posstack > 0);
  }
  else {
    do {
      aodP = (AliAODMCParticle *) mcP->At(posstack);
      if (aodP->GetLabel() == aLabel) return aodP;
      posstack++;
    }
    while (posstack < mcP->GetEntries());
  }
  
  return 0;
}

void AliFemtoEventReaderAOD::CopyPIDtoFemtoTrack(AliAODTrack *tAodTrack, 
						 AliFemtoTrack *tFemtoTrack)
{

	// copying DCA information (taking it from global tracks gives better resolution than from TPC-only)

	double impact[2];
	double covimpact[3];

	if (!tAodTrack->PropagateToDCA(fEvent->GetPrimaryVertex(),fEvent->GetMagneticField(),10000,impact,covimpact)) {
		//cout << "sth went wrong with dca propagation" << endl;
		tFemtoTrack->SetImpactD(-1000.0);
		tFemtoTrack->SetImpactZ(-1000.0);

	}
	else {
		tFemtoTrack->SetImpactD(impact[0]);
		tFemtoTrack->SetImpactZ(impact[1]);
	}

  double aodpid[10];
  tAodTrack->GetPID(aodpid);
  tFemtoTrack->SetPidProbElectron(aodpid[0]);
  tFemtoTrack->SetPidProbMuon(aodpid[1]);
  tFemtoTrack->SetPidProbPion(aodpid[2]);
  tFemtoTrack->SetPidProbKaon(aodpid[3]);
  tFemtoTrack->SetPidProbProton(aodpid[4]);

  aodpid[0] = -100000.0;
  aodpid[1] = -100000.0;
  aodpid[2] = -100000.0;
  aodpid[3] = -100000.0;
  aodpid[4] = -100000.0;
		
  double tTOF = 0.0;

  //what is that code? for what do we need it? nsigma values are not enough?
   if (tAodTrack->GetStatus() & AliESDtrack::kTOFpid) {  //AliESDtrack::kTOFpid=0x8000
     tTOF = tAodTrack->GetTOFsignal();
     tAodTrack->GetIntegratedTimes(aodpid);

     tTOF -= fAODpidUtil->GetTOFResponse().GetStartTime(tAodTrack->P());
   }

   tFemtoTrack->SetTofExpectedTimes(tTOF-aodpid[2], tTOF-aodpid[3], tTOF-aodpid[4]);
 
  //////  TPC ////////////////////////////////////////////

  float nsigmaTPCK=-1000.;                                                  
  float nsigmaTPCPi=-1000.;                                                 
  float nsigmaTPCP=-1000.;                                                  
          
  //   cout<<"in reader fESDpid"<<fESDpid<<endl;

  if (tAodTrack->IsOn(AliESDtrack::kTPCpid)){ //AliESDtrack::kTPCpid=0x0080
    nsigmaTPCK = fAODpidUtil->NumberOfSigmasTPC(tAodTrack,AliPID::kKaon);
    nsigmaTPCPi = fAODpidUtil->NumberOfSigmasTPC(tAodTrack,AliPID::kPion);
    nsigmaTPCP = fAODpidUtil->NumberOfSigmasTPC(tAodTrack,AliPID::kProton);
  }

  tFemtoTrack->SetNSigmaTPCPi(nsigmaTPCPi);
  tFemtoTrack->SetNSigmaTPCK(nsigmaTPCK);
  tFemtoTrack->SetNSigmaTPCP(nsigmaTPCP);

  tFemtoTrack->SetTPCchi2(tAodTrack->Chi2perNDF());       
  tFemtoTrack->SetTPCncls(tAodTrack->GetTPCNcls());       
  tFemtoTrack->SetTPCnclsF(tAodTrack->GetTPCNcls());      
  
  tFemtoTrack->SetTPCsignalN(1); 
  tFemtoTrack->SetTPCsignalS(1); 
  tFemtoTrack->SetTPCsignal(tAodTrack->GetTPCsignal());
 
  ///////TOF//////////////////////

    float vp=-1000.;
    float nsigmaTOFPi=-1000.;
    float nsigmaTOFK=-1000.;
    float nsigmaTOFP=-1000.;

    if ((tAodTrack->GetStatus() & AliESDtrack::kTOFpid) && //AliESDtrack::kTOFpid=0x8000
	(tAodTrack->GetStatus() & AliESDtrack::kTOFout) && //AliESDtrack::kTOFout=0x2000
        (tAodTrack->GetStatus() & AliESDtrack::kTIME) ) //AliESDtrack::kTIME=0x80000000
      {
	if(tAodTrack->IsOn(AliESDtrack::kTOFpid)) //AliESDtrack::kTOFpid=0x8000
	  {

	    nsigmaTOFPi = fAODpidUtil->NumberOfSigmasTOF(tAodTrack,AliPID::kPion);
	    nsigmaTOFK = fAODpidUtil->NumberOfSigmasTOF(tAodTrack,AliPID::kKaon);
	    nsigmaTOFP = fAODpidUtil->NumberOfSigmasTOF(tAodTrack,AliPID::kProton);

	    Double_t len=200;// esdtrack->GetIntegratedLength(); !!!!!
	    Double_t tof=tAodTrack->GetTOFsignal();
	    if(tof > 0.) vp=len/tof/0.03;
	  }
      }
    tFemtoTrack->SetVTOF(vp);
    tFemtoTrack->SetNSigmaTOFPi(nsigmaTOFPi);
    tFemtoTrack->SetNSigmaTOFK(nsigmaTOFK);
    tFemtoTrack->SetNSigmaTOFP(nsigmaTOFP);

    
    //////////////////////////////////////

}

void AliFemtoEventReaderAOD::SetCentralityPreSelection(double min, double max)
{
  fCentRange[0] = min; fCentRange[1] = max;
  fUsePreCent = 1; 
  fEstEventMult = kCentrality;
}

void AliFemtoEventReaderAOD::SetNoCentrality(bool anocent)
{
  if(anocent==false) {fEstEventMult=kCentrality;}
  else {fEstEventMult=kReference; fUsePreCent = 0; }
}

void AliFemtoEventReaderAOD::SetAODpidUtil(AliAODpidUtil *aAODpidUtil)
{
  fAODpidUtil = aAODpidUtil;
  //  printf("fAODpidUtil: %x\n",fAODpidUtil);
}

void AliFemtoEventReaderAOD::SetAODheader(AliAODHeader *aAODheader)
{
  fAODheader = aAODheader;
}


void AliFemtoEventReaderAOD::SetMagneticFieldSign(int s)
{
  if(s>0)
    fMagFieldSign = 1;
  else if(s<0)
    fMagFieldSign = -1;
  else
    fMagFieldSign = 0;
}

void AliFemtoEventReaderAOD::SetEPVZERO(Bool_t iepvz)
{
    fisEPVZ = iepvz;
}

void AliFemtoEventReaderAOD::GetGlobalPositionAtGlobalRadiiThroughTPC(AliAODTrack *track, Float_t bfield, Float_t globalPositionsAtRadii[9][3])
{
  // Gets the global position of the track at nine different radii in the TPC
  // track is the track you want to propagate
  // bfield is the magnetic field of your event
  // globalPositionsAtRadii is the array of global positions in the radii and xyz

  // Initialize the array to something indicating there was no propagation
  for(Int_t i=0;i<9;i++){
    for(Int_t j=0;j<3;j++){
      globalPositionsAtRadii[i][j]=-9999.;
    }
  }

  // Make a copy of the track to not change parameters of the track
  AliExternalTrackParam etp; etp.CopyFromVTrack(track);
  //printf("\nAfter CopyFromVTrack\n");
  //etp.Print();

  // The global position of the the track
  Double_t xyz[3]={-9999.,-9999.,-9999.};

  // Counter for which radius we want
  Int_t iR=0;
  // The radii at which we get the global positions
  // IROC (OROC) from 84.1 cm to 132.1 cm (134.6 cm to 246.6 cm)
  Float_t Rwanted[9]={85.,105.,125.,145.,165.,185.,205.,225.,245.};
  // The global radius we are at
  Float_t globalRadius=0;

  // Propagation is done in local x of the track
  for (Float_t x = etp.GetX();x<247.;x+=1.){ // GetX returns local coordinates
    // Starts at the tracks fX and goes outwards. x = 245 is the outer radial limit
    // of the TPC when the track is straight, i.e. has inifinite pt and doesn't get bent.
    // If the track's momentum is smaller than infinite, it will develop a y-component, which
    // adds to the global radius

    // Stop if the propagation was not succesful. This can happen for low pt tracks
    // that don't reach outer radii
    if(!etp.PropagateTo(x,bfield))break;
    etp.GetXYZ(xyz); // GetXYZ returns global coordinates
    globalRadius = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]); //Idea to speed up: compare squared radii

    // Roughly reached the radius we want
    if(globalRadius > Rwanted[iR]){

      // Bigger loop has bad precision, we're nearly one centimeter too far, go back in small steps.
      while (globalRadius>Rwanted[iR]){
	x-=.1;
	//      printf("propagating to x %5.2f\n",x);
	if(!etp.PropagateTo(x,bfield))break;
	etp.GetXYZ(xyz); // GetXYZ returns global coordinates
	globalRadius = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]); //Idea to speed up: compare squared radii
      }
      //printf("At Radius:%05.2f (local x %5.2f). Setting position to x %4.1f y %4.1f z %4.1f\n",globalRadius,x,xyz[0],xyz[1],xyz[2]);
      globalPositionsAtRadii[iR][0]=xyz[0];
      globalPositionsAtRadii[iR][1]=xyz[1];
      globalPositionsAtRadii[iR][2]=xyz[2];
      // Indicate we want the next radius
      iR+=1;
    }
    if(iR>=8){
      // TPC edge reached
      return;
    }
  }
}


