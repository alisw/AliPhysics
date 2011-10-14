/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
// AliFemtoEventReaderKinematicsChain - the reader class for the Alice ESD and     //
// the model Kinematics information tailored for the Task framework and the        //
// Reads in AliESDfriend to create shared hit/quality information                  //
// Authors: Malgorzata Janik, Warsaw University of Technology, majanik@cern.ch     //
//          Lukasz Graczykowski, Warsaw University of Technology, lgraczyk@cern.ch //
//  	       				                                           //
/////////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoEventReaderKinematicsChain.h"

#include "TFile.h"
#include "TTree.h"
#include "TList.h"

#include "AliFmPhysicalHelixD.h"
#include "AliFmThreeVectorF.h"

#include "SystemOfUnits.h"

#include "AliFemtoEvent.h"

#include "TParticle.h"
#include "AliStack.h"
#include "TParticlePDG.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoModelGlobalHiddenInfo.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"

#include "AliVertexerTracks.h"

ClassImp(AliFemtoEventReaderKinematicsChain)

#if !(ST_NO_NAMESPACES)
  using namespace units;
#endif

using namespace std;
//____________________________
AliFemtoEventReaderKinematicsChain::AliFemtoEventReaderKinematicsChain():
  fFileName(" "),
  fConstrained(true),
  fNumberofEvent(0),
  fCurEvent(0),
  fCurFile(0),
  fStack(0x0),
  fGenHeader(0x0),
  fRotateToEventPlane(0)
{
  //constructor with 0 parameters , look at default settings 
}

//__________________
AliFemtoEventReaderKinematicsChain::AliFemtoEventReaderKinematicsChain(const AliFemtoEventReaderKinematicsChain& aReader):
  AliFemtoEventReader(aReader),
  fFileName(" "),
  fConstrained(true),
  fNumberofEvent(0),
  fCurEvent(0),
  fCurFile(0),
  fStack(0x0),
  fGenHeader(0x0),
  fRotateToEventPlane(0)
{
  // Copy constructor
  fConstrained = aReader.fConstrained;
  fNumberofEvent = aReader.fNumberofEvent;
  fCurEvent = aReader.fCurEvent;
  fCurFile = aReader.fCurFile;
  fStack = aReader.fStack;
  fRotateToEventPlane = aReader.fRotateToEventPlane;
}
//__________________
AliFemtoEventReaderKinematicsChain::~AliFemtoEventReaderKinematicsChain()
{
  //Destructor
  //delete fEvent;
}

//__________________
AliFemtoEventReaderKinematicsChain& AliFemtoEventReaderKinematicsChain::operator=(const AliFemtoEventReaderKinematicsChain& aReader)
{
  // Assignment operator
  if (this == &aReader)
    return *this;

  fConstrained = aReader.fConstrained;
  fNumberofEvent = aReader.fNumberofEvent;
  fCurEvent = aReader.fCurEvent;
  fCurFile = aReader.fCurFile;
  fStack = aReader.fStack;
  fGenHeader = aReader.fGenHeader;
  fRotateToEventPlane = aReader.fRotateToEventPlane;
  return *this;
}
//__________________
// Simple report
AliFemtoString AliFemtoEventReaderKinematicsChain::Report()
{
  AliFemtoString temp = "\n This is the AliFemtoEventReaderKinematicsChain\n";
  return temp;
}

//__________________
void AliFemtoEventReaderKinematicsChain::SetConstrained(const bool constrained)
{
  // Select whether to read constrained or not constrained momentum
  fConstrained=constrained;
}
//__________________
bool AliFemtoEventReaderKinematicsChain::GetConstrained() const
{
  // Check whether we read constrained or not constrained momentum
  return fConstrained;
}
//__________________
AliFemtoEvent* AliFemtoEventReaderKinematicsChain::ReturnHbtEvent()
{
  // Get the event, read all the relevant information from the stack
  // and fill the AliFemtoEvent class
  // Returns a valid AliFemtoEvent
  AliFemtoEvent *hbtEvent = 0;
  string tFriendFileName;

  cout << "AliFemtoEventReaderKinematlaicsChain::Starting to read event: "<<fCurEvent<<endl;
	
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
  double fV1[3];
  double fVCov[6];


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
      if (kinetrack->Eta() < 0.9)
	if(kinetrack->GetPDG()->Charge()/3!=0)
	  tNormMult++;
	  
	  
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
      if(pdgcode==321 || pdgcode==-321 )
        kinepid[3]=1000;
      //pion
      if( pdgcode==211 || pdgcode==-211)
        kinepid[2]=1000;
      //electron
      if(pdgcode==11 || pdgcode==-11)
        kinepid[0]=1000;
      //muon
      if(pdgcode==13 || pdgcode==-13)
        kinepid[1]=1000;

      trackCopy->SetPidProbElectron(kinepid[0]);
      trackCopy->SetPidProbMuon(kinepid[1]);
      trackCopy->SetPidProbPion(kinepid[2]);
      trackCopy->SetPidProbKaon(kinepid[3]);
      trackCopy->SetPidProbProton(kinepid[4]);
					
					
	//Momentum
      double pxyz[3];
      double rxyz[3];
     
	pxyz[0]=kinetrack->Px();
	pxyz[1]=kinetrack->Py();
	pxyz[2]=kinetrack->Pz();

	rxyz[0]=kinetrack->Vx();
	rxyz[1]=kinetrack->Vy();
	rxyz[2]=kinetrack->Vz();

	if (fRotateToEventPlane) {
	  double tPhi = TMath::ATan2(pxyz[1], pxyz[0]);
	  double tRad = TMath::Hypot(pxyz[0], pxyz[1]);
	  
	  pxyz[0] = tRad*TMath::Cos(tPhi - tReactionPlane);
	  pxyz[1] = tRad*TMath::Sin(tPhi - tReactionPlane);
	}

	AliFemtoThreeVector v(pxyz[0],pxyz[1],pxyz[2]);
	if (v.Mag() < 0.0001) {
	  //	  cout << "Found 0 momentum ???? "  << pxyz[0] << " " << pxyz[1] << " " << pxyz[2] << endl;
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
	
		
    }
  
  hbtEvent->SetNumberOfTracks(realnofTracks);//setting number of track which we read in event
  hbtEvent->SetNormalizedMult(tNormMult);
  fCurEvent++;	

  return hbtEvent; 
}

//___________________
void AliFemtoEventReaderKinematicsChain::SetStackSource(AliStack *aStack)
{
  // The chain loads the stack for us
  // You must provide the address where it can be found
  fStack = aStack;
}
//___________________
void AliFemtoEventReaderKinematicsChain::SetGenEventHeader(AliGenEventHeader *aGenHeader)
{
  // The chain loads the generator event header for us
  // You must provide the address where it can be found
  fGenHeader = aGenHeader;
}

//__________________
void AliFemtoEventReaderKinematicsChain::SetRotateToEventPlane(short dorotate)
{
  fRotateToEventPlane=dorotate;
}

Float_t AliFemtoEventReaderKinematicsChain::GetSigmaToVertex(double *impact, double *covar)
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
