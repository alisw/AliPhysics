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

#include "AliFmPhysicalHelixD.h"
#include "AliFmThreeVectorF.h"

#include "SystemOfUnits.h"

#include "AliFemtoEvent.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoModelGlobalHiddenInfo.h"

ClassImp(AliFemtoEventReaderAOD)

#if !(ST_NO_NAMESPACES)
  using namespace units;
#endif

using namespace std;
//____________________________
//constructor with 0 parameters , look at default settings 
AliFemtoEventReaderAOD::AliFemtoEventReaderAOD():
  fNumberofEvent(0),
  fCurEvent(0),
  fEvent(0x0),
  fAllTrue(160),
  fAllFalse(160),
  fFilterBit(0),
  fPWG2AODTracks(0x0),
  fReadMC(0),
  fInputFile(" "),
  fFileName(" "),
  fTree(0x0),
  fAodFile(0x0)
{
  // default constructor
  fAllTrue.ResetAllBits(kTRUE);
  fAllFalse.ResetAllBits(kFALSE);
}

AliFemtoEventReaderAOD::AliFemtoEventReaderAOD(const AliFemtoEventReaderAOD &aReader) :
  AliFemtoEventReader(),
  fNumberofEvent(0),
  fCurEvent(0),
  fEvent(0x0),
  fAllTrue(160),
  fAllFalse(160),
  fFilterBit(0),
  fPWG2AODTracks(0x0),
  fReadMC(0),
  fInputFile(" "),
  fFileName(" "),
  fTree(0x0),
  fAodFile(0x0)
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
	  //	  cout<<"Number of Entries in file "<<fNumberofEvent<<endl;
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
  int tracksPrim=0;     

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
	if (!aodtrack->TestFilterBit(fFilterBit)) {
	  delete trackCopy;
	  continue;
	}

	CopyAODtoFemtoTrack(aodtrack, trackCopy, pwg2aodtrack);
	
	if (mcP) {
	  // Fill the hidden information with the simulated data
	  //	  Int_t pLabel = aodtrack->GetLabel();
	  AliAODMCParticle *tPart = GetParticleWithLabel(mcP, (TMath::Abs(aodtrack->GetLabel())));

	  // Check the mother information
	  
	  // Using the new way of storing the freeze-out information
	  // Final state particle is stored twice on the stack
	  // one copy (mother) is stored with original freeze-out information
	  //   and is not tracked
	  // the other one (daughter) is stored with primary vertex position
	  //   and is tracked
	  
	  // Freeze-out coordinates
	  double fpx=0.0, fpy=0.0, fpz=0.0, fpt=0.0;
	  fpx = tPart->Xv() - fV1[0];
	  fpy = tPart->Yv() - fV1[1];
	  fpz = tPart->Zv() - fV1[2];
	  fpt = tPart->T();

	  AliFemtoModelGlobalHiddenInfo *tInfo = new AliFemtoModelGlobalHiddenInfo();
	  tInfo->SetGlobalEmissionPoint(fpx, fpy, fpz);

	  fpx *= 1e13;
	  fpy *= 1e13;
	  fpz *= 1e13;
	  fpt *= 1e13;
	  
	  //      cout << "Looking for mother ids " << endl;
	  if (motherids[TMath::Abs(aodtrack->GetLabel())]>0) {
	    //	cout << "Got mother id" << endl;
	    AliAODMCParticle *mother = GetParticleWithLabel(mcP, motherids[TMath::Abs(aodtrack->GetLabel())]);
	    // Check if this is the same particle stored twice on the stack
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
	      fpt = mother->T() *1e13*3e10;
	      
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
	  trackCopy->SetHiddenInfo(tInfo);

	}

	double pxyz[3];
	aodtrack->PxPyPz(pxyz);//reading noconstarined momentum
	const AliFmThreeVectorD ktP(pxyz[0],pxyz[1],pxyz[2]);
	// Check the sanity of the tracks - reject zero momentum tracks
	if (ktP.Mag() == 0) {
	  delete trackCopy;
	  continue;
	}
      }
      else {
	// No additional information exists
	// Read in the normal AliAODTracks 
	const AliAODTrack *aodtrack=fEvent->GetTrack(i); // getting the AODtrack directly

	if (aodtrack->IsPrimaryCandidate()) tracksPrim++;
	
// 	if (!aodtrack->TestFilterBit(fFilterBit))
// 	  continue;
	
	CopyAODtoFemtoTrack(aodtrack, trackCopy, 0);
	
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
	aodtrack->PxPyPz(pxyz);//reading noconstarined momentum
	const AliFmThreeVectorD ktP(pxyz[0],pxyz[1],pxyz[2]);
	// Check the sanity of the tracks - reject zero momentum tracks
	if (ktP.Mag() == 0) {
	  delete trackCopy;
	  continue;
	}
      }


      tEvent->TrackCollection()->push_back(trackCopy);//adding track to analysis
      realnofTracks++;//real number of tracks		
    }

  tEvent->SetNumberOfTracks(realnofTracks);//setting number of track which we read in event	
  tEvent->SetNormalizedMult(tracksPrim);

  AliCentrality *cent = fEvent->GetCentrality();
  if (cent) tEvent->SetNormalizedMult((int) cent->GetCentralityPercentile("V0M"));

  if (cent) {
    tEvent->SetCentralityV0(cent->GetCentralityPercentile("V0M"));
    //    tEvent->SetCentralityFMD(cent->GetCentralityPercentile("FMD"));
    tEvent->SetCentralitySPD1(cent->GetCentralityPercentile("CL1"));
    //    tEvent->SetCentralityTrk(cent->GetCentralityPercentile("TRK"));
  }
  

  if (mcP) delete [] motherids;

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
  tFemtoTrack->SetFlags(tAodTrack->GetFlags());
  tFemtoTrack->SetLabel(tAodTrack->GetLabel());
		
  // Track quality information 
  float covmat[6];
  tAodTrack->GetCovMatrix(covmat);  

//   if (TMath::Abs(tAodTrack->Xv()) > 0.00000000001)
//     tFemtoTrack->SetImpactD(TMath::Hypot(tAodTrack->Xv(), tAodTrack->Yv())*(tAodTrack->Xv()/TMath::Abs(tAodTrack->Xv())));
//   else
//     tFemtoTrack->SetImpactD(0.0);
//   tFemtoTrack->SetImpactD(tAodTrack->DCA());
    
//   tFemtoTrack->SetImpactZ(tAodTrack->ZAtDCA());
  tFemtoTrack->SetImpactD(TMath::Hypot(tAodTrack->Xv() - fV1[0], tAodTrack->Yv() - fV1[1]));
  tFemtoTrack->SetImpactZ(tAodTrack->Zv() - fV1[2]);

  //  cout << fV1[0] << " " << tAodTrack->XAtDCA() << " " << tAodTrack->Xv() << endl;

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
    tFemtoTrack->SetTPCClusterMap(tAodTrack->GetTPCClusterMap());
    tFemtoTrack->SetTPCSharedMap(tAodTrack->GetTPCSharedMap());
    
    double xtpc[3] = {0,0,0};
    tFemtoTrack->SetNominalTPCEntrancePoint(xtpc);
    tFemtoTrack->SetNominalTPCExitPoint(xtpc);
  }

  //  cout << "Track has " << TMath::Hypot(tAodTrack->Xv(), tAodTrack->Yv()) << "  " << tAodTrack->Zv() << "  " << tAodTrack->GetTPCNcls() << endl;


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

void AliFemtoEventReaderAOD::SetReadMC(unsigned char a)
{
  fReadMC = a;
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





