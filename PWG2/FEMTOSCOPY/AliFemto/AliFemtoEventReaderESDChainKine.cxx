////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoEventReaderESDChainKine - the reader class for the Alice ESD and   //
// the model Kinematics information tailored for the Task framework and the   //
// Reads in AliESDfriend to create shared hit/quality information             //
// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoEventReaderESDChainKine.h"

#include "TFile.h"
#include "TTree.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"

#include "AliFmPhysicalHelixD.h"
#include "AliFmThreeVectorF.h"

#include "SystemOfUnits.h"

#include "AliFemtoEvent.h"

#include "TParticle.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoModelGlobalHiddenInfo.h"
#include "AliGenHijingEventHeader.h"

#include "AliVertexerTracks.h"

ClassImp(AliFemtoEventReaderESDChainKine)

#if !(ST_NO_NAMESPACES)
  using namespace units;
#endif

using namespace std;
//____________________________
AliFemtoEventReaderESDChainKine::AliFemtoEventReaderESDChainKine():
  fFileName(" "),
  fConstrained(true),
  fUseTPCOnly(false),
  fNumberofEvent(0),
  fCurEvent(0),
  fCurFile(0),
  fEvent(0x0),
  fStack(0x0),
  fGenHeader(0x0)
{
  //constructor with 0 parameters , look at default settings 
}

//__________________
AliFemtoEventReaderESDChainKine::AliFemtoEventReaderESDChainKine(const AliFemtoEventReaderESDChainKine& aReader):
  AliFemtoEventReader(aReader),
  fFileName(" "),
  fConstrained(true),
  fUseTPCOnly(false),
  fNumberofEvent(0),
  fCurEvent(0),
  fCurFile(0),
  fEvent(0x0),
  fStack(0x0),
  fGenHeader(0x0)
{
  // Copy constructor
  fConstrained = aReader.fConstrained;
  fUseTPCOnly = aReader.fUseTPCOnly;
  fNumberofEvent = aReader.fNumberofEvent;
  fCurEvent = aReader.fCurEvent;
  fCurFile = aReader.fCurFile;
  fEvent = new AliESDEvent();
  fStack = aReader.fStack;
}
//__________________
AliFemtoEventReaderESDChainKine::~AliFemtoEventReaderESDChainKine()
{
  //Destructor
  delete fEvent;
}

//__________________
AliFemtoEventReaderESDChainKine& AliFemtoEventReaderESDChainKine::operator=(const AliFemtoEventReaderESDChainKine& aReader)
{
  // Assignment operator
  if (this == &aReader)
    return *this;

  fConstrained = aReader.fConstrained;
  fUseTPCOnly = aReader.fUseTPCOnly;
  fNumberofEvent = aReader.fNumberofEvent;
  fCurEvent = aReader.fCurEvent;
  fCurFile = aReader.fCurFile;
  if (fEvent) delete fEvent;
  fEvent = new AliESDEvent();
  fStack = aReader.fStack;
  fGenHeader = aReader.fGenHeader;

  return *this;
}
//__________________
// Simple report
AliFemtoString AliFemtoEventReaderESDChainKine::Report()
{
  AliFemtoString temp = "\n This is the AliFemtoEventReaderESDChainKine\n";
  return temp;
}

//__________________
void AliFemtoEventReaderESDChainKine::SetConstrained(const bool constrained)
{
  // Select whether to read constrained or not constrained momentum
  fConstrained=constrained;
}
//__________________
bool AliFemtoEventReaderESDChainKine::GetConstrained() const
{
  // Check whether we read constrained or not constrained momentum
  return fConstrained;
}
//__________________
AliFemtoEvent* AliFemtoEventReaderESDChainKine::ReturnHbtEvent()
{
  // Get the event, read all the relevant information
  // and fill the AliFemtoEvent class
  // Returns a valid AliFemtoEvent
  AliFemtoEvent *hbtEvent = 0;
  string tFriendFileName;

  // Get the friend information
  cout<<"starting to read event "<<fCurEvent<<endl;
  //  fEvent->SetESDfriend(fEventFriend);
	
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
  double fVCov[6];
  if (fUseTPCOnly) {
    fEvent->GetPrimaryVertexTPC()->GetXYZ(fV1);
    fEvent->GetPrimaryVertexTPC()->GetCovMatrix(fVCov);
    if (!fEvent->GetPrimaryVertexTPC()->GetStatus())
      fVCov[4] = -1001.0;
  }
  else {
    if (fEvent->GetPrimaryVertex()) {
      fEvent->GetPrimaryVertex()->GetXYZ(fV1);
      fEvent->GetPrimaryVertex()->GetCovMatrix(fVCov);

      if (!fEvent->GetPrimaryVertex()->GetStatus()) {
	// Get the vertex from SPD
	fEvent->GetPrimaryVertexSPD()->GetXYZ(fV1);
	fEvent->GetPrimaryVertexSPD()->GetCovMatrix(fVCov);
	
	
	if (!fEvent->GetPrimaryVertexSPD()->GetStatus())
	  fVCov[4] = -1001.0;
	else {
	  fEvent->GetPrimaryVertexSPD()->GetXYZ(fV1);
	  fEvent->GetPrimaryVertexSPD()->GetCovMatrix(fVCov);
	}
      }
    }
    else {
      if (fEvent->GetPrimaryVertexSPD()) {
	fEvent->GetPrimaryVertexSPD()->GetXYZ(fV1);
	fEvent->GetPrimaryVertexSPD()->GetCovMatrix(fVCov);
      }
    }
    if ((!fEvent->GetPrimaryVertex()) && (!fEvent->GetPrimaryVertexSPD()))
      {
	cout << "No vertex found !!!" << endl;
	fV1[0] = 10000.0;
	fV1[1] = 10000.0;
	fV1[2] = 10000.0;
	fVCov[4] = -1001.0;
      }
  }

  AliFmThreeVectorF vertex(fV1[0],fV1[1],fV1[2]);
  
  hbtEvent->SetPrimVertPos(vertex);
  hbtEvent->SetPrimVertCov(fVCov);

  AliGenHijingEventHeader *hdh = dynamic_cast<AliGenHijingEventHeader *> (fGenHeader);
	
  Double_t tReactionPlane = 0;
  if (hdh)
    {
      tReactionPlane = hdh->ReactionPlaneAngle();
    }
  //starting to reading tracks
  int nofTracks=0;  //number of reconstructed tracks in event
  nofTracks=fEvent->GetNumberOfTracks();
  int realnofTracks=0;//number of track which we use ina analysis

  Int_t *motherids;
  motherids = new Int_t[fStack->GetNtrack()];
  for (int ip=0; ip<fStack->GetNtrack(); ip++) motherids[ip] = 0;

  // Read in mother ids
  TParticle *motherpart;
  for (int ip=0; ip<fStack->GetNtrack(); ip++) {
    motherpart = fStack->Particle(ip);
    if (motherpart->GetDaughter(0) > 0)
      motherids[motherpart->GetDaughter(0)] = ip;
    if (motherpart->GetDaughter(1) > 0)
      motherids[motherpart->GetDaughter(1)] = ip;

//     if (motherpart->GetPdgCode() == 211) {
//       cout << "Mother " << ip << " has daughters " 
// 	   << motherpart->GetDaughter(0) << " " 
// 	   << motherpart->GetDaughter(1) << " " 
// 	   << motherpart->Vx() << " " 
// 	   << motherpart->Vy() << " " 
// 	   << motherpart->Vz() << " " 
// 	   << endl;
      
//     }
  }

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

	AliFemtoThreeVector v(pxyz[0],pxyz[1],pxyz[2]);
	if (v.mag() < 0.0001) {
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
	if (fConstrained==true)		    
	  tGoodMomentum=esdtrack->GetConstrainedPxPyPz(pxyz); //reading constrained momentum
	else
	  tGoodMomentum=esdtrack->GetPxPyPz(pxyz);//reading noconstarined momentum
	
	AliFemtoThreeVector v(pxyz[0],pxyz[1],pxyz[2]);
	if (v.mag() < 0.0001) {
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
	esdtrack->GetImpactParameters(imp,cim);

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
      trackCopy->SetITSncls(esdtrack->GetNcls(0));     
      trackCopy->SetTPCchi2(esdtrack->GetTPCchi2());       
      trackCopy->SetTPCncls(esdtrack->GetTPCNcls());       
      trackCopy->SetTPCnclsF(esdtrack->GetTPCNclsF());      
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

      // Fill the hidden information with the simulated data
      TParticle *tPart = fStack->Particle(TMath::Abs(esdtrack->GetLabel()));

      // Check the mother information

      // Using the new way of storing the freeze-out information
      // Final state particle is stored twice on the stack
      // one copy (mother) is stored with original freeze-out information
      //   and is not tracked
      // the other one (daughter) is stored with primary vertex position
      //   and is tracked
	
      // Freeze-out coordinates
      double fpx=0.0, fpy=0.0, fpz=0.0, fpt=0.0;
      fpx = tPart->Vx() - fV1[0];
      fpy = tPart->Vy() - fV1[1];
      fpz = tPart->Vz() - fV1[2];
      fpt = tPart->T();

      AliFemtoModelGlobalHiddenInfo *tInfo = new AliFemtoModelGlobalHiddenInfo();
      tInfo->SetGlobalEmissionPoint(fpx, fpy, fpz);

      fpx *= 1e13;
      fpy *= 1e13;
      fpz *= 1e13;
      fpt *= 1e13*3e10;

      if (motherids[TMath::Abs(esdtrack->GetLabel())]>0) {
 	TParticle *mother = fStack->Particle(motherids[TMath::Abs(esdtrack->GetLabel())]);
 	// Check if this is the same particle stored twice on the stack
  	if ((mother->GetPdgCode() == tPart->GetPdgCode() || (mother->Px() == tPart->Px()))) {
	  // It is the same particle
	  // Read in the original freeze-out information
	  // and convert it from to [fm]
 	  fpx = mother->Vx()*1e13;
 	  fpy = mother->Vy()*1e13;
 	  fpz = mother->Vz()*1e13;
 	  fpt = mother->T()*1e13*3e10;
	  
	}
      }

      tInfo->SetPDGPid(tPart->GetPdgCode());
      tInfo->SetTrueMomentum(tPart->Px(), tPart->Py(), tPart->Pz());
      Double_t mass2 = (tPart->Energy() *tPart->Energy() -
			tPart->Px()*tPart->Px() -
			tPart->Py()*tPart->Py() -
			tPart->Pz()*tPart->Pz());
      if (mass2>0.0)
	tInfo->SetMass(TMath::Sqrt(mass2));
      else 
	tInfo->SetMass(0.0);

      tInfo->SetEmissionPoint(fpx, fpy, fpz, fpt);
      trackCopy->SetHiddenInfo(tInfo);
      
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
  fCurEvent++;	
  cout<<"end of reading nt "<<nofTracks<<" real number "<<realnofTracks<<endl;
  return hbtEvent; 
}
//___________________
void AliFemtoEventReaderESDChainKine::SetESDSource(AliESDEvent *aESD)
{
  // The chain loads the ESD for us
  // You must provide the address where it can be found
  fEvent = aESD;
}
//___________________
void AliFemtoEventReaderESDChainKine::SetStackSource(AliStack *aStack)
{
  // The chain loads the stack for us
  // You must provide the address where it can be found
  fStack = aStack;
}
//___________________
void AliFemtoEventReaderESDChainKine::SetGenEventHeader(AliGenEventHeader *aGenHeader)
{
  // The chain loads the generator event header for us
  // You must provide the address where it can be found
  fGenHeader = aGenHeader;
}

//__________________
void AliFemtoEventReaderESDChainKine::SetUseTPCOnly(const bool usetpconly)
{
  fUseTPCOnly=usetpconly;
}

bool AliFemtoEventReaderESDChainKine::GetUseTPCOnly() const
{
  return fUseTPCOnly;
}
Float_t AliFemtoEventReaderESDChainKine::GetSigmaToVertex(double *impact, double *covar)
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
