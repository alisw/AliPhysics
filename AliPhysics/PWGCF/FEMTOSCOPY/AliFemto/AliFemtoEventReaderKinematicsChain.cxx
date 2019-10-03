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

//#include "AliVertexerTracks.h"


#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoEventReaderKinematicsChain);
  /// \endcond
#endif


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
  fEstEventMult(kGlobalCount),
  fRotateToEventPlane(0),
  fReadOnlyPrimaries(true),
  fReadOnlyPrimariesV0(true),
  fReadPrimariesSecWeakMaterial(false),
  fReadPrimariesSecWeakMaterialV0(false),
  fIsMisalignment(false),
  fRemoveWeakDecaysInMC(false),
  fDiscardStatusCodeFlag(false),
  fDiscardStatusCode(0)
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
  fEstEventMult(kGlobalCount),
  fRotateToEventPlane(0),
  fReadOnlyPrimaries(true),
  fReadOnlyPrimariesV0(true),
  fReadPrimariesSecWeakMaterial(false),
  fReadPrimariesSecWeakMaterialV0(false),
  fIsMisalignment(false),
  fRemoveWeakDecaysInMC(false),
  fDiscardStatusCodeFlag(false),
  fDiscardStatusCode(0)
{
  // Copy constructor
  fConstrained = aReader.fConstrained;
  fNumberofEvent = aReader.fNumberofEvent;
  fCurEvent = aReader.fCurEvent;
  fCurFile = aReader.fCurFile;
  fStack = aReader.fStack;
  fEstEventMult = aReader.fEstEventMult;
  fRotateToEventPlane = aReader.fRotateToEventPlane;
  fReadOnlyPrimaries = aReader.fReadOnlyPrimaries;
  fReadOnlyPrimariesV0 = aReader.fReadOnlyPrimariesV0;
  fReadPrimariesSecWeakMaterial = aReader.fReadPrimariesSecWeakMaterial;
  fReadPrimariesSecWeakMaterialV0 = aReader.fReadPrimariesSecWeakMaterialV0;
  fIsMisalignment  = aReader.fIsMisalignment;
  fRemoveWeakDecaysInMC = aReader.fRemoveWeakDecaysInMC;
  fDiscardStatusCodeFlag = aReader.fDiscardStatusCodeFlag;
  fDiscardStatusCode = aReader.fDiscardStatusCode;
  
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
  fEstEventMult = aReader.fEstEventMult;
  fRotateToEventPlane = aReader.fRotateToEventPlane;
  fReadOnlyPrimaries = aReader.fReadOnlyPrimaries;
  fReadOnlyPrimariesV0 = aReader.fReadOnlyPrimariesV0;
  fReadPrimariesSecWeakMaterial = aReader.fReadPrimariesSecWeakMaterial;
  fReadPrimariesSecWeakMaterialV0 = aReader.fReadPrimariesSecWeakMaterialV0;
  fIsMisalignment  = aReader.fIsMisalignment;
  fRemoveWeakDecaysInMC = aReader.fRemoveWeakDecaysInMC;
  fDiscardStatusCodeFlag = aReader.fDiscardStatusCodeFlag;
  fDiscardStatusCode = aReader.fDiscardStatusCode;
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

void AliFemtoEventReaderKinematicsChain::ReadOnlyPrimaries(bool primaries)
{
  fReadOnlyPrimaries = primaries;
}

void AliFemtoEventReaderKinematicsChain::ReadOnlyPrimariesV0(bool primaries)
{
  fReadOnlyPrimariesV0 = primaries;
}

void AliFemtoEventReaderKinematicsChain::ReadPrimariesSecWeakMaterial(bool primaries)
{
  fReadPrimariesSecWeakMaterial = primaries;
}

void AliFemtoEventReaderKinematicsChain::ReadPrimariesSecWeakMaterialV0(bool primaries)
{
  fReadPrimariesSecWeakMaterialV0 = primaries;
}


//__________________
AliFemtoEvent* AliFemtoEventReaderKinematicsChain::ReturnHbtEvent()
{
  // Get the event, read all the relevant information from the stack
  // and fill the AliFemtoEvent class
  // Returns a valid AliFemtoEvent
  AliFemtoEvent *hbtEvent = 0;
  string tFriendFileName;

  //cout << "AliFemtoEventReaderKinematlaicsChain::Starting to read event: "<<fCurEvent<<endl;

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
      //cout << "Got reaction plane " << tReactionPlane << endl;
    }
  
  
  hbtEvent->SetReactionPlaneAngle(tReactionPlane);

  //starting to reading tracks
  int nofTracks=0;  //number of all tracks in MC event
  nofTracks=fStack->GetNtrack();
  int realnofTracks=0;//number of track which we use in analysis


  int tNormMult = 0;
  int tV0direction = 0;

  double previousTrackPt = 0;
  for (int i=0;i<nofTracks;i++)	
    {	
     
      //getting next track
      TParticle *kinetrack= fStack->Particle(i);

      if(fReadOnlyPrimaries)
	{
	  //take only primaries
	  if(!fStack->IsPhysicalPrimary(i)) {  continue;}
	}
      else if(fReadPrimariesSecWeakMaterial)
	{
	  //take only primaries
	  if(!(fStack->IsPhysicalPrimary(i) || fStack->IsSecondaryFromWeakDecay(i) || fStack->IsSecondaryFromMaterial(i))) {  continue;}
	}

      AliFemtoTrack* trackCopy = new AliFemtoTrack();
       
      if(fIsMisalignment){
	TParticle *motherParticle1;
	Int_t motherIndex1 = kinetrack->GetFirstMother();
	bool deleted_kinetrack = false;


	if((kinetrack->GetPdgCode()==3122)|| (kinetrack->GetPdgCode()==-3122))
	  {
            if (motherIndex1 != -1){
             
	      motherParticle1 = fStack->Particle(motherIndex1);
	      if(motherParticle1->GetPdgCode()==kinetrack->GetPdgCode())
		{
		  delete trackCopy;
		  deleted_kinetrack = true;
		  continue;
            
		}

	      //   cout << "mother after: "<< motherParticle1->GetPdgCode()<< "kinetrack after: "<< kinetrack->GetPdgCode()<<endl;

	    }
	  }

      }



      if (fRemoveWeakDecaysInMC)
	{
	  // Exclude weak decay products (if not done by IsPhysicalPrimary)
	  // In order to prevent analyzing daughters from weak decays 
	  // - AMPT does not only strong decays, so IsPhysicalPrimary does not catch it

	  //Checking mother
	  /*
	  Int_t motherIndex = kinetrack->GetFirstMother();
	  if (motherIndex != -1){
	    TParticle *motherParticle = fStack->Particle(motherIndex);
		if (motherParticle){
			int pdgcode =  motherParticle->GetPdgCode();
		
			const Int_t kNWeakParticles = 7;
			const Int_t kWeakParticles[kNWeakParticles] = { 3322, 3312, 3222, // Xi0 Xi+- Sigma-+
								  3122, 3112, // Lambda0 Sigma+-
								  130, 310 // K_L0 K_S0
			};

			bool fromWeak = false;

			for (Int_t j=0; j != kNWeakParticles; ++j) {
				if (kWeakParticles[j] == pdgcode) {
					//std::cout<<Form("Removing particle %d (pdg code mother %d)", i, motherParticle->GetPdgCode())<<std::endl;
					fromWeak=true;
				break;
				}
			}
				if(fromWeak) {
					delete trackCopy;
					continue;
				}
		}
		}*/

	  //checking distance from the interaction point
	  double vtx = kinetrack->Vx();
	  double vty = kinetrack->Vy();
	
	  if(sqrt(vtx*vtx+vty*vty)<0.1){
	    delete trackCopy;
	    continue;
	  }

     // if(sqrt(vtx*vtx+vty*vty)<0.1){
     // 	if (motherIndex1 != -1){
     // 	  TParticle *motherParticle2 = fStack->Particle(motherIndex1);
     // 	  cout<<"PASS Particle "<<i<<" PDG: "<<kinetrack->GetPdgCode()<<" Mother: "<<motherParticle2->GetPdgCode()<<" "<<kinetrack->Vx()<<" "<<kinetrack->Vy()<<" "<<kinetrack->Vz()<<endl;
     // 	}
     // 	else
     // 	  cout<<"PASS Particle "<<i<<" PDG: "<<kinetrack->GetPdgCode()<<" Mother: "<<"-1"<<" "<<kinetrack->Vx()<<" "<<kinetrack->Vy()<<" "<<kinetrack->Vz()<<endl;
     //  }
     //  else{
     // 	if (motherIndex1 != -1){
     // 	  TParticle *motherParticle2 = fStack->Particle(motherIndex1);
     // 	  cout<<"FAIL Particle "<<i<<" PDG: "<<kinetrack->GetPdgCode()<<" Mother: "<<motherParticle2->GetPdgCode()<<" "<<kinetrack->Vx()<<" "<<kinetrack->Vy()<<" "<<kinetrack->Vz()<<endl;
     // 	}
     // 	else
     // 	  cout<<"FAIL Particle "<<i<<" PDG: "<<kinetrack->GetPdgCode()<<" Mother: "<<"-1"<<" "<<kinetrack->Vx()<<" "<<kinetrack->Vy()<<" "<<kinetrack->Vz()<<endl;
     //  }	

	  
	}

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
      else if(pdgcode==3122 || pdgcode==-3122 || abs(pdgcode)==310 ) //Lambda, AntiLambda, K0
	{; }
      else if(pdgcode==3312 || pdgcode==-3312) //Xi-, Xi+
	{; }
      else {
		delete trackCopy;
		continue;
      }
  


      trackCopy->SetPidProbElectron(kinepid[0]);
      trackCopy->SetPidProbMuon(kinepid[1]);
      trackCopy->SetPidProbPion(kinepid[2]);
      trackCopy->SetPidProbKaon(kinepid[3]);
      trackCopy->SetPidProbProton(kinepid[4]);


	//Momentum
      double pxyz[3];
      double rxyz[3];

      if(kinetrack->Px()==0 && kinetrack->Py()==0)
        continue;

      pxyz[0]=kinetrack->Px();
      pxyz[1]=kinetrack->Py();
      pxyz[2]=kinetrack->Pz();

      rxyz[0]=kinetrack->Vx();
      rxyz[1]=kinetrack->Vy();
      rxyz[2]=kinetrack->Vz();

      AliFemtoModelHiddenInfo *tInfo = new AliFemtoModelHiddenInfo();
      tInfo->SetPDGPid(pdgcode);
      tInfo->SetTrueMomentum(pxyz[0], pxyz[1], pxyz[2]);
      tInfo->SetMass(kinetrack->GetMass());
      tInfo->SetEmissionPoint(rxyz[0]-fV1[0], rxyz[1]-fV1[1], rxyz[2]-fV1[2], 0.0);
      trackCopy->SetHiddenInfo(tInfo);

      trackCopy->SetTrueMomentum(pxyz[0], pxyz[1], pxyz[2]);
      trackCopy->SetEmissionPoint(rxyz[0]-fV1[0], rxyz[1]-fV1[1], rxyz[2]-fV1[2], 0.0);


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

      if(fIsMisalignment){
	if(previousTrackPt==trackCopy->Pt())
	  continue;
	previousTrackPt=trackCopy->Pt();



      }
      
      //remove particles with specific status code
      if(fDiscardStatusCodeFlag)
	{
	  if(kinetrack->GetStatusCode() == fDiscardStatusCode) continue;
	}
      
      hbtEvent->TrackCollection()->push_back(trackCopy);//adding track to analysis
      //cout<<"Track added: "<<i<<endl;

    }

  if(fEstEventMult == kImpactParameter){
    AliGenEventHeader* eventHeader = 0;
    AliGenCocktailEventHeader* cocktailHeader = dynamic_cast<AliGenCocktailEventHeader*> (fGenHeader);
    if (cocktailHeader)
      eventHeader = dynamic_cast<AliGenEventHeader*> (cocktailHeader->GetHeaders()->First());
    else
      eventHeader = dynamic_cast<AliGenEventHeader *> (fGenHeader);

    if (!eventHeader)
      {
	// We avoid AliFatal here, because the AOD productions sometimes have events where the MC header is missing 
	// (due to unreadable Kinematics) and we don't want to loose the whole job because of a few events
	std::cout<<"Event header not found. Skipping this event."<<std::endl;
	return 0;
      }
      
    AliCollisionGeometry* collGeometry = dynamic_cast<AliCollisionGeometry*> (eventHeader);
    if (!collGeometry)
      {
	eventHeader->Dump();
	std::cout<<"Asking for MC_b centrality, but event header has no collision geometry information"<<std::endl;
      }
      
    tNormMult = 100 * collGeometry->ImpactParameter();
    //cout<<"Impact: "<<collGeometry->ImpactParameter()<<endl;
  }


  hbtEvent->SetNumberOfTracks(realnofTracks);//setting number of track which we read in event
  if (fEstEventMult == kGlobalCount)
    hbtEvent->SetNormalizedMult(tNormMult);
  else if(fEstEventMult == kVZERO)
    hbtEvent->SetNormalizedMult(tV0direction);
  else if(fEstEventMult == kImpactParameter)
    hbtEvent->SetNormalizedMult(tNormMult);

  fCurEvent++;


  //V0 analysis code - no V0 finder for Kinematics, we can only check if it is primary and if it has at least 2 daughters.

  for (int i=0;i<nofTracks;i++)
    {

      if(fReadOnlyPrimariesV0)
	{
	  //take only primaries
	  if(!fStack->IsPhysicalPrimary(i)) {continue;}
	}
      else if(fReadPrimariesSecWeakMaterialV0)
	{
	  //take only primaries
	  if(!(fStack->IsPhysicalPrimary(i) || fStack->IsSecondaryFromWeakDecay(i) || fStack->IsSecondaryFromMaterial(i))) {continue;}
	}

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


  //cout<<"Number of tracks: "<<realnofTracks<<endl;

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

void AliFemtoEventReaderKinematicsChain::SetUseMultiplicity(EstEventMult aType)
{
  fEstEventMult = aType;
}

void AliFemtoEventReaderKinematicsChain::IsMisalignment(bool isMisalignment){
  fIsMisalignment = isMisalignment;
}

void AliFemtoEventReaderKinematicsChain::RemoveWeakDecaysManually(bool removeWeakDecaysInMC){
  fRemoveWeakDecaysInMC = removeWeakDecaysInMC;
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



 void AliFemtoEventReaderKinematicsChain::CopyAODtoFemtoV0(TParticle *tv0, AliFemtoV0 *tFemtoV0 )
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



void AliFemtoEventReaderKinematicsChain::DiscardStatusCode(bool flag, int code)
{
  fDiscardStatusCodeFlag = flag;
  fDiscardStatusCode = code;
}
