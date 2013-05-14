/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
// AliFemtoEventReaderKinematicsChainESD - the reader class for the Alice ESD and     //
// the model Kinematics information tailored for the Task framework and the        //
// Reads in AliESDfriend to create shared hit/quality information                  //
// Authors: Malgorzata Janik, Warsaw University of Technology, majanik@cern.ch     //
//          Lukasz Graczykowski, Warsaw University of Technology, lgraczyk@cern.ch //
//  	       				                                           //
/////////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoEventReaderKinematicsChainESD.h"

#include "TFile.h"
#include "TTree.h"
#include "TList.h"

#include "AliFmPhysicalHelixD.h"
#include "AliFmThreeVectorF.h"

#include "SystemOfUnits.h"

#include "AliFemtoEvent.h"

#include "TParticle.h"
#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliCentrality.h"
#include "TParticlePDG.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoModelGlobalHiddenInfo.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"

#include "AliVertexerTracks.h"

ClassImp(AliFemtoEventReaderKinematicsChainESD)

#if !(ST_NO_NAMESPACES)
  using namespace units;
#endif

using namespace std;
//____________________________
AliFemtoEventReaderKinematicsChainESD::AliFemtoEventReaderKinematicsChainESD():
  fFileName(" "),
  fConstrained(true),
  fNumberofEvent(0),
  fCurEvent(0),
  fCurFile(0),
  fStack(0x0),
  fEvent(0x0),
  fGenHeader(0x0),
  fEstEventMult(kGlobalCount),
  fRotateToEventPlane(0)
{
  //constructor with 0 parameters , look at default settings 
}

//__________________
AliFemtoEventReaderKinematicsChainESD::AliFemtoEventReaderKinematicsChainESD(const AliFemtoEventReaderKinematicsChainESD& aReader):
  AliFemtoEventReader(aReader),
  fFileName(" "),
  fConstrained(true),
  fNumberofEvent(0),
  fCurEvent(0),
  fCurFile(0),
  fStack(0x0),
  fEvent(0x0),
  fGenHeader(0x0),
  fEstEventMult(kGlobalCount),
  fRotateToEventPlane(0)
{
  // Copy constructor
  fConstrained = aReader.fConstrained;
  fNumberofEvent = aReader.fNumberofEvent;
  fCurEvent = aReader.fCurEvent;
  fCurFile = aReader.fCurFile;
  fStack = aReader.fStack;
  fEvent = aReader.fEvent;
  fEstEventMult = aReader.fEstEventMult;
  fRotateToEventPlane = aReader.fRotateToEventPlane;
}
//__________________
AliFemtoEventReaderKinematicsChainESD::~AliFemtoEventReaderKinematicsChainESD()
{
  //Destructor
  //delete fEvent;
}

//__________________
AliFemtoEventReaderKinematicsChainESD& AliFemtoEventReaderKinematicsChainESD::operator=(const AliFemtoEventReaderKinematicsChainESD& aReader)
{
  // Assignment operator
  if (this == &aReader)
    return *this;

  fConstrained = aReader.fConstrained;
  fNumberofEvent = aReader.fNumberofEvent;
  fCurEvent = aReader.fCurEvent;
  fCurFile = aReader.fCurFile;
  fStack = aReader.fStack;
  fEvent = aReader.fEvent;
  fGenHeader = aReader.fGenHeader;
  fEstEventMult = aReader.fEstEventMult;
  fRotateToEventPlane = aReader.fRotateToEventPlane;
  return *this;
}
//__________________
// Simple report
AliFemtoString AliFemtoEventReaderKinematicsChainESD::Report()
{
  AliFemtoString temp = "\n This is the AliFemtoEventReaderKinematicsChainESD\n";
  return temp;
}

//__________________
void AliFemtoEventReaderKinematicsChainESD::SetConstrained(const bool constrained)
{
  // Select whether to read constrained or not constrained momentum
  fConstrained=constrained;
}
//__________________
bool AliFemtoEventReaderKinematicsChainESD::GetConstrained() const
{
  // Check whether we read constrained or not constrained momentum
  return fConstrained;
}
//__________________
AliFemtoEvent* AliFemtoEventReaderKinematicsChainESD::ReturnHbtEvent()
{
  // Get the event, read all the relevant information from the stack
  // and fill the AliFemtoEvent class
  // Returns a valid AliFemtoEvent
  AliFemtoEvent *hbtEvent = 0;
  string tFriendFileName;

  cout << "AliFemtoEventReaderKinematicsChainESD::Starting to read event: "<<fCurEvent<<endl;
	
  hbtEvent = new AliFemtoEvent;
  //setting basic things
  //  hbtEvent->SetEventNumber(fEvent->GetEventNumber());
  hbtEvent->SetRunNumber(0); //No Run number in Kinematics!
  hbtEvent->SetMagneticField(0*kilogauss);//to check if here is ok
  hbtEvent->SetZDCN1Energy(0);
  hbtEvent->SetZDCP1Energy(0);
  hbtEvent->SetZDCN2Energy(0);
  hbtEvent->SetZDCP2Energy(0);
  hbtEvent->SetZDCEMEnergy(0);
  hbtEvent->SetZDCParticipants(0);
  hbtEvent->SetTriggerMask(0);
  hbtEvent->SetTriggerCluster(0);
	
  //Vertex
  double fV1[3] = {0.0,0.0,0.0};
  double fVCov[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};


  AliFmThreeVectorF vertex(0,0,0);
    
  
  hbtEvent->SetPrimVertPos(vertex);
  hbtEvent->SetPrimVertCov(fVCov);

  Double_t tReactionPlane = 0;

  AliGenHijingEventHeader *hdh = dynamic_cast<AliGenHijingEventHeader *> (fGenHeader);	
  if (!hdh) {
    AliGenCocktailEventHeader *cdh = dynamic_cast<AliGenCocktailEventHeader *> (fGenHeader);
    if (cdh) {
      TList *tGenHeaders = cdh->GetHeaders();
      for (int ihead = 0; ihead<tGenHeaders->GetEntries(); ihead++) {
	hdh = dynamic_cast<AliGenHijingEventHeader *> (fGenHeader);	
	if (hdh) break;
      }
    }
  }

  if (hdh)
    {
      tReactionPlane = hdh->ReactionPlaneAngle();
      cout << "Got reaction plane " << tReactionPlane << endl;
    }

  hbtEvent->SetReactionPlaneAngle(tReactionPlane);

  //starting to reading tracks
  int nofTracks=0;  //number of all tracks in MC event
  nofTracks=fStack->GetNtrack(); 
  int realnofTracks=0;//number of track which we use in analysis


  int tNormMult = 0;
  int tV0direction = 0;
  for (int i=0;i<nofTracks;i++)
    {
      //take only primaries
      if(!fStack->IsPhysicalPrimary(i)) {continue;}
	  	  
      AliFemtoTrack* trackCopy = new AliFemtoTrack();	
	
      	  //getting next track
      TParticle *kinetrack= fStack->Particle(i);

      //setting multiplicity
        realnofTracks++;//real number of tracks (only primary particles)

      //setting normalized multiplicity
      
	if(kinetrack->GetPDG()->Charge()/3!=0)
	  if (kinetrack->Pt() > 0.15 && kinetrack->Pt() < 20) 
	    if (kinetrack->Eta() < 0.8)
	      tNormMult++;
	  
	//counting particles that go into direction of VZERO detector
	if(kinetrack->Eta() > 2.8 && kinetrack->Eta() < 5.1) //VZERO-A
	  tV0direction++;
	if(kinetrack->Eta() > -3.7 && kinetrack->Eta() < -1.7)//VZERO-C
	  tV0direction++;	
 
	  //charge
      trackCopy->SetCharge((short)(fStack->Particle(i)->GetPDG()->Charge()/3));


      //in aliroot we have AliPID 
      //0-electron 1-muon 2-pion 3-kaon 4-proton 5-photon 6-pi0 7-neutron 8-kaon0 9-eleCon   
      //we use only 5 first
      double kinepid[5];
      for(int pid_iter=0;pid_iter<5;pid_iter++)
	  kinepid[pid_iter]=0;

      int pdgcode = kinetrack->GetPdgCode();
      //proton
      if(pdgcode==2212 || pdgcode==-2212)
        kinepid[4]=1000;
      //kaon
      else if(pdgcode==321 || pdgcode==-321 )
        kinepid[3]=1000;
      //pion
      else if( pdgcode==211 || pdgcode==-211)
        kinepid[2]=1000;
      //electron
      else if(pdgcode==11 || pdgcode==-11)
        kinepid[0]=1000;
      //muon
      else if(pdgcode==13 || pdgcode==-13)
        kinepid[1]=1000;
      else
	continue;
      trackCopy->SetPidProbElectron(kinepid[0]);
      trackCopy->SetPidProbMuon(kinepid[1]);
      trackCopy->SetPidProbPion(kinepid[2]);
      trackCopy->SetPidProbKaon(kinepid[3]);
      trackCopy->SetPidProbProton(kinepid[4]);
					
					
	//Momentum
      double pxyz[3];
      // double rxyz[3];
     
	pxyz[0]=kinetrack->Px();
	pxyz[1]=kinetrack->Py();
	pxyz[2]=kinetrack->Pz();

	// rxyz[0]=kinetrack->Vx();
	// rxyz[1]=kinetrack->Vy();
	// rxyz[2]=kinetrack->Vz();

	if (fRotateToEventPlane) {
	  double tPhi = TMath::ATan2(pxyz[1], pxyz[0]);
	  double tRad = TMath::Hypot(pxyz[0], pxyz[1]);
	  
	  pxyz[0] = tRad*TMath::Cos(tPhi - tReactionPlane);
	  pxyz[1] = tRad*TMath::Sin(tPhi - tReactionPlane);
	}

	AliFemtoThreeVector v(pxyz[0],pxyz[1],pxyz[2]);
	if (v.Mag() < 0.0001) {
	  //cout << "Found 0 momentum ???? "  << pxyz[0] << " " << pxyz[1] << " " << pxyz[2] << endl;
	  delete trackCopy;
	  continue;
	}

	trackCopy->SetP(v);//setting momentum
	trackCopy->SetPt(sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]));
	const AliFmThreeVectorD kP(pxyz[0],pxyz[1],pxyz[2]);
	const AliFmThreeVectorD kOrigin(fV1[0],fV1[1],fV1[2]);

	//label
	trackCopy->SetLabel(i);


	hbtEvent->TrackCollection()->push_back(trackCopy);//adding track to analysis
	//cout<<"Track added: "<<i<<endl;
		
    }
  
  hbtEvent->SetNumberOfTracks(realnofTracks);//setting number of track which we read in event

  if (fEstEventMult == kGlobalCount) 
    hbtEvent->SetNormalizedMult(tNormMult);
  else if(fEstEventMult == kVZERO)
    hbtEvent->SetNormalizedMult(tV0direction);
  else if (fEstEventMult == kReferenceITSTPC)
    {
      hbtEvent->SetNormalizedMult(AliESDtrackCuts::GetReferenceMultiplicity(fEvent,AliESDtrackCuts::kTrackletsITSTPC,1.2));
    }
  else if (fEstEventMult == kCentrality) {
    // centrality between 0 (central) and 1 (very peripheral)
    AliCentrality *cent = fEvent->GetCentrality();
    if (cent) {
      if (cent->GetCentralityPercentile("V0M") < 0.00001)
	hbtEvent->SetNormalizedMult(-1);
      else
	hbtEvent->SetNormalizedMult(lrint(10.0*cent->GetCentralityPercentile("V0M")));
      if (Debug()>1) printf ("Set Centrality %i %f %li\n", hbtEvent->UncorrectedNumberOfPrimaries(), 
			     10.0*cent->GetCentralityPercentile("V0M"), lrint(10.0*cent->GetCentralityPercentile("V0M")));
    }
  }

  fCurEvent++;	


  //V0 analysis code - no V0 finder for Kinematics, we can only check if it is primary and if it has at least 2 daughters.

  for (int i=0;i<nofTracks;i++)
    {
      //do not take primaries
      if(!fStack->IsPhysicalPrimary(i)) {continue;}
      //getting next track
      TParticle *kinetrack= fStack->Particle(i);
      if (!kinetrack) continue;
      
      if(kinetrack->GetPDG()->Charge()!=0) continue; //charge - neutral
      //if(kinetrack->GetPDG()->Stable()==1) continue; //particle is not stable
      if(kinetrack->GetDaughter(0)<1) continue; //has 1'st daughter
      if(kinetrack->GetDaughter(1)<1) continue;  //has 2'nd daughter
 
      //we want one positive, one negative particle. Or two neutral.
      // if((fStack->Particle(kinetrack->GetDaughter(0)))->GetPDG()->Charge()>=0)
      // 	if((fStack->Particle(kinetrack->GetDaughter(1)))->GetPDG()->Charge()>0)
      // 	  continue;
      // if((fStack->Particle(kinetrack->GetDaughter(0)))->GetPDG()->Charge()<=0)
      // 	if((fStack->Particle(kinetrack->GetDaughter(0)))->GetPDG()->Charge()<0)
      // 	  continue;
      
      if(kinetrack->Pt()<0.00001)
	continue;
      
      AliFemtoV0* trackCopyV0 = new AliFemtoV0();
      CopyAODtoFemtoV0(kinetrack, trackCopyV0);
      hbtEvent->V0Collection()->push_back(trackCopyV0);
    //cout<<"Pushback v0 to v0collection"<<endl;
    }


  cout<<"Number of tracks: "<<realnofTracks<<endl;

  return hbtEvent; 
}

//___________________
void AliFemtoEventReaderKinematicsChainESD::SetStackSource(AliStack *aStack)
{
  // The chain loads the stack for us
  // You must provide the address where it can be found
  fStack = aStack;
}
void AliFemtoEventReaderKinematicsChainESD::SetESDSource(AliESDEvent *aESD)
{
  // The chain loads the ESD for us
  // You must provide the address where it can be found
  fEvent = aESD;
}
//___________________
void AliFemtoEventReaderKinematicsChainESD::SetGenEventHeader(AliGenEventHeader *aGenHeader)
{
  // The chain loads the generator event header for us
  // You must provide the address where it can be found
  fGenHeader = aGenHeader;
}

//__________________
void AliFemtoEventReaderKinematicsChainESD::SetRotateToEventPlane(short dorotate)
{
  fRotateToEventPlane=dorotate;
}

void AliFemtoEventReaderKinematicsChainESD::SetUseMultiplicity(EstEventMult aType)
{
  fEstEventMult = aType;
}

Float_t AliFemtoEventReaderKinematicsChainESD::GetSigmaToVertex(double *impact, double *covar)
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



 void AliFemtoEventReaderKinematicsChainESD::CopyAODtoFemtoV0(TParticle *tv0, AliFemtoV0 *tFemtoV0 )
{
  tFemtoV0->SetEtaV0(tv0->Eta());
  tFemtoV0->SetEtaV0(tv0->Phi());
  tFemtoV0->SetptV0(tv0->Pt()); 
  tFemtoV0->SetptotV0(tv0->P());

  tFemtoV0->SetmomV0X(tv0->Px());
  tFemtoV0->SetmomV0Y(tv0->Py());
  tFemtoV0->SetmomV0Z(tv0->Pz());
  AliFemtoThreeVector momv0(tv0->Px(),tv0->Py(),tv0->Pz());
  tFemtoV0->SetmomV0(momv0);


  TParticle *trackpos;
  TParticle *trackneg;

  //daughters
  if(fStack->Particle(tv0->GetDaughter(0))->GetPDG()->Charge()>=0) //first positive, second negative
    {
      trackpos = (TParticle*)(fStack->Particle(tv0->GetDaughter(0)));
      trackneg = (TParticle*)(fStack->Particle(tv0->GetDaughter(1)));
      tFemtoV0->SetidPos(tv0->GetDaughter(0));
      tFemtoV0->SetidNeg(tv0->GetDaughter(1));
    }
  else //first negative, second positive
    {
      trackpos = (TParticle*)(fStack->Particle(tv0->GetDaughter(1)));
      trackneg = (TParticle*)(fStack->Particle(tv0->GetDaughter(0)));
      tFemtoV0->SetidPos(tv0->GetDaughter(1));
      tFemtoV0->SetidNeg(tv0->GetDaughter(0));
    }

  tFemtoV0->SetEtaPos(trackpos->Eta());
  tFemtoV0->SetEtaNeg(trackneg->Eta());
  
  tFemtoV0->SetptPos(trackpos->Pt());
  tFemtoV0->SetptNeg(trackneg->Pt());
   
  tFemtoV0->SetptotPos(trackpos->P());
  tFemtoV0->SetptotNeg(trackneg->P());

  tFemtoV0->SetmomPosX(trackpos->Px());
  tFemtoV0->SetmomPosY(trackpos->Py());
  tFemtoV0->SetmomPosZ(trackpos->Pz());
  AliFemtoThreeVector mompos(trackpos->Px(),trackpos->Py(),trackpos->Pz());
  tFemtoV0->SetmomPos(mompos);

  tFemtoV0->SetmomNegX(trackneg->Px());
  tFemtoV0->SetmomNegY(trackneg->Py());
  tFemtoV0->SetmomNegZ(trackneg->Pz());
  AliFemtoThreeVector momneg(trackneg->Px(),trackneg->Py(),trackneg->Pz());
  tFemtoV0->SetmomNeg(momneg);

    
  tFemtoV0->SetmassLambda(tv0->GetMass());
  tFemtoV0->SetmassAntiLambda(tv0->GetMass());
  tFemtoV0->SetmassK0Short(tv0->GetMass());

  tFemtoV0->SetYV0(tv0->Y());

  tFemtoV0->SetdecayVertexV0X(trackpos->Vx()); //vertex of the decay is set as the vertex of creation of daughters
  tFemtoV0->SetdecayVertexV0Y(trackpos->Vy());
  tFemtoV0->SetdecayVertexV0Z(trackpos->Vz());
  AliFemtoThreeVector decayvertex(trackpos->Vx(),trackpos->Vy(),trackpos->Vz());
  tFemtoV0->SetdecayVertexV0(decayvertex);

  tFemtoV0->SetdcaV0Daughters(0);
  tFemtoV0->SetCosPointingAngle(1);


  tFemtoV0->SetStatusPos(1);
  tFemtoV0->SetStatusNeg(1);

  
  if(trackpos->GetPdgCode()==2212) //proton
    {
      tFemtoV0->SetPosNSigmaTPCK(1000);
      tFemtoV0->SetPosNSigmaTPCPi(1000);
      tFemtoV0->SetPosNSigmaTPCP(0);
    }
  if(trackneg->GetPdgCode()==-2212) //antiproton
    {
      tFemtoV0->SetNegNSigmaTPCK(1000);
      tFemtoV0->SetNegNSigmaTPCPi(1000);
      tFemtoV0->SetNegNSigmaTPCP(0);
    }
  if(trackpos->GetPdgCode()==211) //pion plus
    {
      tFemtoV0->SetPosNSigmaTPCK(1000);
      tFemtoV0->SetPosNSigmaTPCPi(0);
      tFemtoV0->SetPosNSigmaTPCP(1000);
    }
  if(trackneg->GetPdgCode()==-211) //pion minus
    {
      tFemtoV0->SetNegNSigmaTPCK(1000);
      tFemtoV0->SetNegNSigmaTPCPi(0);
      tFemtoV0->SetNegNSigmaTPCP(1000);
    }
  if(trackpos->GetPdgCode()==321) //K+
    {
      tFemtoV0->SetPosNSigmaTPCK(0);
      tFemtoV0->SetPosNSigmaTPCPi(1000);
      tFemtoV0->SetPosNSigmaTPCP(1000);
    }
  if(trackneg->GetPdgCode()==-321) //K-
    {
      tFemtoV0->SetNegNSigmaTPCK(0);
      tFemtoV0->SetNegNSigmaTPCPi(1000);
      tFemtoV0->SetNegNSigmaTPCP(1000);
    }

  
}
