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
#include "TList.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"

#include "AliMultiplicity.h"
#include "AliCentrality.h"
#include "AliEventplane.h"

#include "AliFmPhysicalHelixD.h"
#include "AliFmThreeVectorF.h"

#include "SystemOfUnits.h"

#include "AliFemtoEvent.h"

#include "TParticle.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoModelGlobalHiddenInfo.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"

#include "AliVertexerTracks.h"

#include "AliPID.h"

ClassImp(AliFemtoEventReaderESDChainKine)

#if !(ST_NO_NAMESPACES)
  using namespace units;
#endif

using namespace std;
//____________________________
AliFemtoEventReaderESDChainKine::AliFemtoEventReaderESDChainKine():
  fFileName(" "),
  fConstrained(true),
  fReadInner(false),
  fUseTPCOnly(false),
  fNumberofEvent(0),
  fCurEvent(0),
  fCurFile(0),
  fEvent(0x0),
  fStack(0x0),
  fGenHeader(0x0),
  fTrackType(kGlobal),
  fEstEventMult(kV0Centrality),
  fRotateToEventPlane(0),
  fESDpid(0),
  fIsPidOwner(0),
  fMagFieldSign(0),
  fReadV0(0)

{
  //constructor with 0 parameters , look at default settings 
}

//__________________
AliFemtoEventReaderESDChainKine::AliFemtoEventReaderESDChainKine(const AliFemtoEventReaderESDChainKine& aReader):
  AliFemtoEventReader(aReader),
  fFileName(" "),
  fConstrained(true),
  fReadInner(false),
  fUseTPCOnly(false),
  fNumberofEvent(0),
  fCurEvent(0),
  fCurFile(0),
  fEvent(0x0),
  fStack(0x0),
  fGenHeader(0x0),
  fTrackType(kGlobal),
  fEstEventMult(kV0Centrality),
  fRotateToEventPlane(0),
  fESDpid(0),
  fIsPidOwner(0),
  fMagFieldSign(0),
  fReadV0(0)
{
  // Copy constructor
  fConstrained = aReader.fConstrained;
  fReadInner = aReader.fReadInner;
  fUseTPCOnly = aReader.fUseTPCOnly;
  fNumberofEvent = aReader.fNumberofEvent;
  fCurEvent = aReader.fCurEvent;
  fCurFile = aReader.fCurFile;
  fEvent = new AliESDEvent();
  fStack = aReader.fStack;
  fTrackType = aReader.fTrackType;
  fEstEventMult = aReader.fEstEventMult;
  fRotateToEventPlane = aReader.fRotateToEventPlane;
  fESDpid = aReader.fESDpid;
  fIsPidOwner = aReader.fIsPidOwner;
  fMagFieldSign = aReader.fMagFieldSign;
  fReadV0 = aReader.fReadV0;
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
  fRotateToEventPlane = aReader.fRotateToEventPlane;
  fESDpid = aReader.fESDpid;
  fIsPidOwner = aReader.fIsPidOwner;
  fMagFieldSign = aReader.fMagFieldSign;
  fReadV0 = aReader.fReadV0;

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
void AliFemtoEventReaderESDChainKine::SetReadTPCInner(const bool readinner)
{
  fReadInner=readinner;
}

bool AliFemtoEventReaderESDChainKine::GetReadTPCInner() const
{
  return fReadInner;
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
void AliFemtoEventReaderESDChainKine::SetUseMultiplicity(EstEventMult aType)
{
  fEstEventMult = aType;
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
  cout << "AliFemtoEventReaderESDChainKine::Starting to read event: "<<fCurEvent<<endl;
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
  //hbtEvent->SetTriggerCluster(fEvent->GetTriggerCluster());

  if ((fEvent->IsTriggerClassFired("CINT1WU-B-NOPF-ALL")) ||
      (fEvent->IsTriggerClassFired("CINT1B-ABCE-NOPF-ALL")) ||
      (fEvent->IsTriggerClassFired("CINT1-B-NOPF-ALLNOTRD")) ||
      (fEvent->IsTriggerClassFired("CINT1-B-NOPF-FASTNOTRD")))
    hbtEvent->SetTriggerCluster(1);
  else if ((fEvent->IsTriggerClassFired("CSH1WU-B-NOPF-ALL")) ||
	   (fEvent->IsTriggerClassFired("CSH1-B-NOPF-ALLNOTRD")))
    hbtEvent->SetTriggerCluster(2);
  else 
    hbtEvent->SetTriggerCluster(0);

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
    fEvent->GetPrimaryVertex()->GetXYZ(fV1);
    fEvent->GetPrimaryVertex()->GetCovMatrix(fVCov);
    if (!fEvent->GetPrimaryVertex()->GetStatus())
      fVCov[4] = -1001.0;
  }

  AliFmThreeVectorF vertex(fV1[0],fV1[1],fV1[2]);
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

  Int_t spdetaonecount = 0;
  
  for (int iter=0; iter<fEvent->GetMultiplicity()->GetNumberOfTracklets(); iter++) 
    if (fabs(fEvent->GetMultiplicity()->GetEta(iter)) < 1.0)
      spdetaonecount++;

  //  hbtEvent->SetSPDMult(fEvent->GetMultiplicity()->GetNumberOfTracklets());
  hbtEvent->SetSPDMult(spdetaonecount);

  //starting to reading tracks
  int nofTracks=0;  //number of reconstructed tracks in event
  nofTracks=fEvent->GetNumberOfTracks();
  int realnofTracks=0;//number of track which we use ina analysis

  int tNormMult = 0;
  int tNormMultPos = 0;
  int tNormMultNeg = 0;

  Float_t tTotalPt = 0.0;

  Float_t b[2];
  Float_t bCov[3];

  Int_t tTracklet=0, tITSTPC=0, tITSPure=0;
  
  fEvent->EstimateMultiplicity(tTracklet, tITSTPC, tITSPure, 1.2);
  
  hbtEvent->SetMultiplicityEstimateITSTPC(tITSTPC);
  hbtEvent->SetMultiplicityEstimateTracklets(tTracklet);
  //  hbtEvent->SetMultiplicityEstimateITSPure(tITSPure);
  hbtEvent->SetMultiplicityEstimateITSPure(fEvent->GetMultiplicity()->GetNumberOfITSClusters(1));

  Int_t *motherids;
  if (fStack) {
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
  }
  else {
    printf ("No Stack ???\n");
    delete hbtEvent;
    return 0;
  }

  for (int i=0;i<nofTracks;i++)
    {
      //      cout << "Reading track " << i << endl;
      bool  tGoodMomentum=true; //flaga to chcek if we can read momentum of this track
		
	
      const AliESDtrack *esdtrack=fEvent->GetTrack(i);//getting next track
      //      const AliESDfriendTrack *tESDfriendTrack = esdtrack->GetFriendTrack();


      if ((esdtrack->GetStatus() & AliESDtrack::kTPCrefit) &&
	  (esdtrack->GetStatus() & AliESDtrack::kITSrefit)) {
	if (esdtrack->GetTPCNcls() > 70) 
	  if (esdtrack->GetTPCchi2()/esdtrack->GetTPCNcls() < 4.0) {
	    if (TMath::Abs(esdtrack->Eta()) < 1.2) {
	      esdtrack->GetImpactParameters(b,bCov);
	      if ((b[0]<0.2) && (b[1] < 0.25)) {
		tNormMult++;
		tTotalPt += esdtrack->Pt();
	      }
	    }
	  }
      }
      else if (esdtrack->GetStatus() & AliESDtrack::kTPCrefit) {
	if (esdtrack->GetTPCNcls() > 100) 
	  if (esdtrack->GetTPCchi2()/esdtrack->GetTPCNcls() < 4.0) {
	    if (TMath::Abs(esdtrack->Eta()) < 1.2) {
	      esdtrack->GetImpactParameters(b,bCov);
	      if ((b[0]<2.4) && (b[1] < 3.2)) {
		tNormMult++;
		tTotalPt += esdtrack->Pt();
	      }
	    }
	  }
      }
      
      hbtEvent->SetZDCEMEnergy(tTotalPt);
      //       if (esdtrack->GetStatus() & AliESDtrack::kTPCrefit)
      // 	if (esdtrack->GetTPCNcls() > 80) 
      // 	  if (esdtrack->GetTPCchi2()/esdtrack->GetTPCNcls() < 6.0) 
      // 	    if (esdtrack->GetConstrainedParam())
      // 	      if (fabs(esdtrack->GetConstrainedParam()->Eta()) < 0.5)
      // 		if (esdtrack->GetConstrainedParam()->Pt() < 1.0) {
      // 		  if (esdtrack->GetSign() > 0)
      // 		    tNormMultPos++;
      // 		  else if (esdtrack->GetSign() < 0)
      // 		    tNormMultNeg--;
      // 		}

      // If reading ITS-only tracks, reject all with TPC
      if (fTrackType == kITSOnly) {
	if (esdtrack->GetStatus() & AliESDtrack::kTPCrefit) continue;
	if (!(esdtrack->GetStatus() & AliESDtrack::kITSrefit)) continue;
	if (esdtrack->GetStatus() & AliESDtrack::kTPCin) continue;
	UChar_t iclm = esdtrack->GetITSClusterMap();
	Int_t incls = 0;
	for (int iter=0; iter<6; iter++) if (iclm&(1<<iter)) incls++;
	if (incls<=3) {
	  if (Debug()>1) cout << "Rejecting track with " << incls << " clusters" << endl;
	  continue;
	}
      }



      AliFemtoTrack* trackCopy = new AliFemtoTrack();
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

      esdpid[0] = -100000.0;
      esdpid[1] = -100000.0;
      esdpid[2] = -100000.0;
      esdpid[3] = -100000.0;
      esdpid[4] = -100000.0;

      double tTOF = 0.0;

      if (esdtrack->GetStatus()&AliESDtrack::kTOFpid) {
	tTOF = esdtrack->GetTOFsignal();
	esdtrack->GetIntegratedTimes(esdpid);
      }

      trackCopy->SetTofExpectedTimes(tTOF-esdpid[2], tTOF-esdpid[3], tTOF-esdpid[4]);


      //////  TPC ////////////////////////////////////////////
      float nsigmaTPCK=-1000.;                                                        
      float nsigmaTPCPi=-1000.;                                                        
      float nsigmaTPCP=-1000.;  

     if ((fESDpid) && (esdtrack->IsOn(AliESDtrack::kTPCpid))){
        nsigmaTPCK = fESDpid->NumberOfSigmasTPC(esdtrack,AliPID::kKaon);
        nsigmaTPCPi = fESDpid->NumberOfSigmasTPC(esdtrack,AliPID::kPion);
        nsigmaTPCP = fESDpid->NumberOfSigmasTPC(esdtrack,AliPID::kProton);

      }
      trackCopy->SetNSigmaTPCPi(nsigmaTPCPi);
      trackCopy->SetNSigmaTPCK(nsigmaTPCK);
      trackCopy->SetNSigmaTPCP(nsigmaTPCP);

      ///// TOF ///////////////////////////////////////////////

      float vp=-1000.;
      float nsigmaTOFPi=-1000.;
      float nsigmaTOFK=-1000.;
      float nsigmaTOFP=-1000.;

	if ((esdtrack->GetStatus()&AliESDtrack::kTOFpid) &&
	    (esdtrack->GetStatus()&AliESDtrack::kTOFout) &&
	    (esdtrack->GetStatus()&AliESDtrack::kTIME) &&
	    !(esdtrack->GetStatus()&AliESDtrack::kTOFmismatch))
	  {

	    //if ((esdtrack->GetStatus()&AliESDtrack::kTOFpid) &&
	    //(esdtrack->GetStatus()&AliESDtrack::kTOFout) &&
	    //(esdtrack->GetStatus()&AliESDtrack::kTIME)){
	    // collect info from ESDpid class
	  
	      if ((fESDpid) && (esdtrack->IsOn(AliESDtrack::kTOFpid))) {
	      
	      
	      double tZero = fESDpid->GetTOFResponse().GetStartTime(esdtrack->P());
	      
	      nsigmaTOFPi = fESDpid->NumberOfSigmasTOF(esdtrack,AliPID::kPion,tZero);
	      nsigmaTOFK = fESDpid->NumberOfSigmasTOF(esdtrack,AliPID::kKaon,tZero);
	      nsigmaTOFP = fESDpid->NumberOfSigmasTOF(esdtrack,AliPID::kProton,tZero);
	      
	      Double_t len=esdtrack->GetIntegratedLength();
	      Double_t tof=esdtrack->GetTOFsignal();
	      if(tof > 0.) vp=len/tof/0.03;
	    }
	  }

      trackCopy->SetVTOF(vp);
      trackCopy->SetNSigmaTOFPi(nsigmaTOFPi);
      trackCopy->SetNSigmaTOFK(nsigmaTOFK);
      trackCopy->SetNSigmaTOFP(nsigmaTOFP);

						
      double pxyz[3];
      double rxyz[3];
      double impact[2];
      double covimpact[3];
      
      if (fUseTPCOnly) {
	if (!esdtrack->GetTPCInnerParam()) {
	  cout << "No TPC inner param !" << endl;
	  delete trackCopy;
	  continue;
	}

	AliExternalTrackParam *param = new AliExternalTrackParam(*esdtrack->GetTPCInnerParam());
	param->GetXYZ(rxyz);
	param->PropagateToDCA(fEvent->GetPrimaryVertexTPC(), (fEvent->GetMagneticField()), 10000, impact, covimpact);
	param->GetPxPyPz(pxyz);//reading noconstarined momentum

	  if (fReadInner == true) {
	    AliFemtoModelHiddenInfo *tInfo = new AliFemtoModelHiddenInfo();
	    tInfo->SetPDGPid(211);
	    tInfo->SetTrueMomentum(pxyz[0], pxyz[1], pxyz[2]);
	    tInfo->SetMass(0.13957);
	    //	  tInfo->SetEmissionPoint(rxyz[0], rxyz[1], rxyz[2], 0.0);
	    //	  tInfo->SetEmissionPoint(fV1[0], fV1[1], fV1[2], 0.0);
	    tInfo->SetEmissionPoint(rxyz[0]-fV1[0], rxyz[1]-fV1[1], rxyz[2]-fV1[2], 0.0);
	    trackCopy->SetHiddenInfo(tInfo);
	  }

	if (fRotateToEventPlane) {
	  double tPhi = TMath::ATan2(pxyz[1], pxyz[0]);
	  double tRad = TMath::Hypot(pxyz[0], pxyz[1]);
	  
	  pxyz[0] = tRad*TMath::Cos(tPhi - tReactionPlane);
	  pxyz[1] = tRad*TMath::Sin(tPhi - tReactionPlane);
	}

	AliFemtoThreeVector v(pxyz[0],pxyz[1],pxyz[2]);
	if (v.Mag() < 0.0001) {
	  //	  cout << "Found 0 momentum ???? " << pxyz[0] << " " << pxyz[1] << " " << pxyz[2] << endl;
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
	  if (fReadInner == true) {
	  
	    if (esdtrack->GetTPCInnerParam()) {
	      AliExternalTrackParam *param = new AliExternalTrackParam(*esdtrack->GetTPCInnerParam());
	      trackCopy->SetInnerMomentum(esdtrack->GetTPCmomentum());
	      param->GetXYZ(rxyz);
	      //	    param->PropagateToDCA(fEvent->GetPrimaryVertex(), (fEvent->GetMagneticField()), 10000);
	      param->GetPxPyPz(pxyz);//reading noconstarined momentum
	      delete param;
	    
	      AliFemtoModelHiddenInfo *tInfo = new AliFemtoModelHiddenInfo();
	      tInfo->SetPDGPid(211);
	      tInfo->SetTrueMomentum(pxyz[0], pxyz[1], pxyz[2]);
	      tInfo->SetMass(0.13957);
	      //	    tInfo->SetEmissionPoint(rxyz[0], rxyz[1], rxyz[2], 0.0);
	      //tInfo->SetEmissionPoint(fV1[0], fV1[1], fV1[2], 0.0);
	      tInfo->SetEmissionPoint(rxyz[0]-fV1[0], rxyz[1]-fV1[1], rxyz[2]-fV1[2], 0.0);
	      trackCopy->SetHiddenInfo(tInfo);
	    }
	  }


	if(fTrackType == kGlobal) {
	if (fConstrained==true)		    
	  tGoodMomentum=esdtrack->GetConstrainedPxPyPz(pxyz); //reading constrained momentum
	else
	  tGoodMomentum=esdtrack->GetPxPyPz(pxyz);//reading noconstarined momentum
	}
	  else if (fTrackType == kTPCOnly) {
	    if (esdtrack->GetTPCInnerParam())
	      esdtrack->GetTPCInnerParam()->GetPxPyPz(pxyz);
	    else {
	      delete trackCopy;
	      continue;
	    }
	  }
	  else if (fTrackType == kITSOnly) {
	    if (fConstrained==true)		    
	      tGoodMomentum=esdtrack->GetConstrainedPxPyPz(pxyz); //reading constrained momentum
	    else
	      tGoodMomentum=esdtrack->GetPxPyPz(pxyz);//reading noconstarined momentum
	  }

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

	for (int ii=0; ii<6; ii++){
	  trackCopy->SetITSHitOnLayer(ii,esdtrack->HasPointOnITSLayer(ii));
	}


      // Fill the hidden information with the simulated data
      if (TMath::Abs(esdtrack->GetLabel()) < fStack->GetNtrack()) {
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
	fpt *= 1e13;
	
	// fillDCA

	if (TMath::Abs(impact[0]) > 0.001) {
	  if (fStack->IsPhysicalPrimary(TMath::Abs(esdtrack->GetLabel()))){
	    trackCopy->SetImpactDprim(impact[0]);
	    //cout << "prim" << endl;

	  }
	  else if (fStack->IsSecondaryFromWeakDecay(TMath::Abs(esdtrack->GetLabel()))) {
	    trackCopy->SetImpactDweak(impact[0]);
	    //cout << "wea" << endl;
	  }
	  else if (fStack->IsSecondaryFromMaterial(TMath::Abs(esdtrack->GetLabel()))) {
	    trackCopy->SetImpactDmat(impact[0]);
	    //cout << "mat" << endl;
	  }
	}
	//  end fillDCA 

	//      cout << "Looking for mother ids " << endl;
	if (motherids[TMath::Abs(esdtrack->GetLabel())]>0) {
	  //	cout << "Got mother id" << endl;
	  TParticle *mother = fStack->Particle(motherids[TMath::Abs(esdtrack->GetLabel())]);
	  // Check if this is the same particle stored twice on the stack
	  if ((mother->GetPdgCode() == tPart->GetPdgCode() || (mother->Px() == tPart->Px()))) {
	    // It is the same particle
	    // Read in the original freeze-out information
	    // and convert it from to [fm]
	    
	    // EPOS style 
	    fpx = mother->Vx()*1e13*0.197327;
	    fpy = mother->Vy()*1e13*0.197327;
	    fpz = mother->Vz()*1e13*0.197327;
	    fpt = mother->T() *1e13*0.197327;
	    
	    
	    // Therminator style 
// 	    fpx = mother->Vx()*1e13;
// 	    fpy = mother->Vy()*1e13;
// 	    fpz = mother->Vz()*1e13;
// 	    fpt = mother->T() *1e13*3e10;
	    
	  }
	}
	
	if (fRotateToEventPlane) {
	  double tPhi = TMath::ATan2(fpy, fpx);
	  double tRad = TMath::Hypot(fpx, fpy);
	  
	  fpx = tRad*TMath::Cos(tPhi - tReactionPlane);
	  fpy = tRad*TMath::Sin(tPhi - tReactionPlane);
	}
	
	tInfo->SetPDGPid(tPart->GetPdgCode());
	
	if (fRotateToEventPlane) {
	  double tPhi = TMath::ATan2(tPart->Py(), tPart->Px());
	  double tRad = TMath::Hypot(tPart->Px(), tPart->Py());
	  
	  tInfo->SetTrueMomentum(tRad*TMath::Cos(tPhi - tReactionPlane),
				 tRad*TMath::Sin(tPhi - tReactionPlane),
				 tPart->Pz());
	}
	else
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
      }
      else {
	AliFemtoModelGlobalHiddenInfo *tInfo = new AliFemtoModelGlobalHiddenInfo();
	tInfo->SetMass(0.0);
	double fpx=0.0, fpy=0.0, fpz=0.0, fpt=0.0;
	fpx = fV1[0]*1e13;
	fpy = fV1[1]*1e13;
	fpz = fV1[2]*1e13;
	fpt = 0.0;
	tInfo->SetEmissionPoint(fpx, fpy, fpz, fpt);

	tInfo->SetTrueMomentum(pxyz[0],pxyz[1],pxyz[2]);
	
	trackCopy->SetHiddenInfo(tInfo);
      }
      //      cout << "Got freeze-out " << fpx << " " << fpy << " " << fpz << " " << fpt << " " <<  mass2 << " " << tPart->GetPdgCode() << endl; 
      
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

  delete [] motherids;
  
  hbtEvent->SetNumberOfTracks(realnofTracks);//setting number of track which we read in event	

  AliCentrality *cent = fEvent->GetCentrality();
  if (cent) {
    hbtEvent->SetCentralityV0(cent->GetCentralityPercentile("V0M"));
    //    hbtEvent->SetCentralityFMD(cent->GetCentralityPercentile("FMD"));
    hbtEvent->SetCentralitySPD1(cent->GetCentralityPercentile("CL1"));
    //    hbtEvent->SetCentralityTrk(cent->GetCentralityPercentile("TRK"));

    if (Debug()>1) printf("  FemtoReader Got Event with %f %f %f %f\n", cent->GetCentralityPercentile("V0M"), 0.0, cent->GetCentralityPercentile("CL1"), 0.0);
  }

  if (fEstEventMult == kGlobalCount) 
    hbtEvent->SetNormalizedMult(tNormMult);
  else if (fEstEventMult == kReferenceITSTPC)
    hbtEvent->SetNormalizedMult(AliESDtrackCuts::GetReferenceMultiplicity(fEvent,AliESDtrackCuts::kTrackletsITSTPC,1.2));
  else if(fEstEventMult == kReferenceITSSA)
    hbtEvent->SetNormalizedMult(AliESDtrackCuts::GetReferenceMultiplicity(fEvent,AliESDtrackCuts::kTrackletsITSSA,1.2));
  else if(fEstEventMult == kReferenceTracklets)
    hbtEvent->SetNormalizedMult(AliESDtrackCuts::GetReferenceMultiplicity(fEvent,AliESDtrackCuts::kTracklets,1.2));
  else if (fEstEventMult == kTracklet)
    hbtEvent->SetNormalizedMult(tTracklet);
  else if (fEstEventMult == kITSTPC)
    hbtEvent->SetNormalizedMult(tITSTPC);
  else if (fEstEventMult == kITSPure)
    hbtEvent->SetNormalizedMult(tITSPure);
  else if (fEstEventMult == kSPDLayer1)
    hbtEvent->SetNormalizedMult(fEvent->GetMultiplicity()->GetNumberOfITSClusters(1));
  else if (fEstEventMult == kV0Centrality) {
    // centrality between 0 (central) and 1 (very peripheral)

    if (cent) {
      if (cent->GetCentralityPercentile("V0M") < 0.00001)
	hbtEvent->SetNormalizedMult(-1);
      else
	hbtEvent->SetNormalizedMult(lrint(10.0*cent->GetCentralityPercentile("V0M")));
      if (Debug()>1) printf ("Set Centrality %i %f %li\n", hbtEvent->UncorrectedNumberOfPrimaries(), 
			     10.0*cent->GetCentralityPercentile("V0M"), lrint(10.0*cent->GetCentralityPercentile("V0M")));
    }
  }

  if (tNormMultPos > tNormMultNeg)
    hbtEvent->SetZDCParticipants(tNormMultPos);
  else
    hbtEvent->SetZDCParticipants(tNormMultNeg);
  
  AliEventplane* ep = fEvent->GetEventplane();
  if (ep) {
    hbtEvent->SetEP(ep);
    hbtEvent->SetReactionPlaneAngle(ep->GetEventplane("Q"));
  }

  if(fReadV0)
    {
      for (Int_t i = 0; i < fEvent->GetNumberOfV0s(); i++) {
	AliESDv0* esdv0 = fEvent->GetV0(i);
	if (!esdv0) continue;
	//if(esdv0->GetNDaughters()>2) continue;
	//if(esdv0->GetNProngs()>2) continue;
	if(esdv0->Charge()!=0) continue;
	AliESDtrack *trackPos = fEvent->GetTrack(esdv0->GetPindex());
	if(!trackPos) continue;
	AliESDtrack *trackNeg = fEvent->GetTrack(esdv0->GetNindex());
	if(!trackNeg) continue;
	if(trackPos->Charge()==trackNeg->Charge()) continue;
	AliFemtoV0* trackCopyV0 = new AliFemtoV0();
	CopyESDtoFemtoV0(esdv0, trackCopyV0, fEvent);
	hbtEvent->V0Collection()->push_back(trackCopyV0);
	//cout<<"Pushback v0 to v0collection"<<endl;
      }
    }


  fCurEvent++;	
  //  cout<<"end of reading nt "<<nofTracks<<" real number "<<realnofTracks<<endl;
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
void AliFemtoEventReaderESDChainKine::SetRotateToEventPlane(short dorotate)
{
  fRotateToEventPlane=dorotate;
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

void AliFemtoEventReaderESDChainKine::SetReadTrackType(ReadTrackType aType)
{
  fTrackType = aType;
}

void AliFemtoEventReaderESDChainKine::SetReadV0(bool a)
{
  fReadV0 = a;
}

void AliFemtoEventReaderESDChainKine::CopyESDtoFemtoV0(AliESDv0 *tESDv0, AliFemtoV0 *tFemtoV0, AliESDEvent *tESDevent)
{
  double fPrimaryVtxPosition[3];
  tESDevent->GetPrimaryVertex()->GetXYZ(fPrimaryVtxPosition);
  tFemtoV0->SetdcaV0ToPrimVertex(tESDv0->GetD(fPrimaryVtxPosition[0],fPrimaryVtxPosition[1],fPrimaryVtxPosition[2]));

  //tFemtoV0->SetdecayLengthV0(tESDv0->DecayLengthV0()); //wrocic do tego
  //tFemtoV0->SetdecayVertexV0X(tESDv0->DecayVertexV0X());
  //tFemtoV0->SetdecayVertexV0Y(tESDv0->DecayVertexV0Y());
  //tFemtoV0->SetdecayVertexV0Z(tESDv0->DecayVertexV0Z()); //nie ma w AliESDv0
  //AliFemtoThreeVector decayvertex(tESDv0->DecayVertexV0X(),tESDv0->DecayVertexV0Y(),tESDv0->DecayVertexV0Z());
  //tFemtoV0->SetdecayVertexV0(decayvertex);
  tFemtoV0->SetdcaV0Daughters(tESDv0->GetDcaV0Daughters());
  tFemtoV0->SetmomV0X(tESDv0->Px());
  tFemtoV0->SetmomV0Y(tESDv0->Py());
  tFemtoV0->SetmomV0Z(tESDv0->Pz());
  AliFemtoThreeVector momv0(tESDv0->Px(),tESDv0->Py(),tESDv0->Pz());
  tFemtoV0->SetmomV0(momv0);
  tFemtoV0->SetalphaV0(tESDv0->AlphaV0());
  tFemtoV0->SetptArmV0(tESDv0->PtArmV0());
  //tFemtoV0->SeteLambda(tESDv0->ELambda());
  //tFemtoV0->SeteK0Short(tESDv0->EK0Short());
  //tFemtoV0->SetePosProton(tESDv0->EPosProton());
  //tFemtoV0->SeteNegProton(tESDv0->ENegProton());
  tFemtoV0->SetmassLambda(tESDv0->GetEffMass(4,2));
  tFemtoV0->SetmassAntiLambda(tESDv0->GetEffMass(2,4));
  tFemtoV0->SetmassK0Short(tESDv0->GetEffMass(2,2));
  //tFemtoV0->SetrapLambda(tESDv0->RapLambda());
  //tFemtoV0->SetrapK0Short(tESDv0->RapK0Short());
  tFemtoV0->SetptV0(tESDv0->Pt());
  tFemtoV0->SetptotV0(tESDv0->P());
  tFemtoV0->SetEtaV0(tESDv0->Eta());
  tFemtoV0->SetPhiV0(tESDv0->Phi());
  tFemtoV0->SetCosPointingAngle(tESDv0->GetV0CosineOfPointingAngle(fPrimaryVtxPosition[0],fPrimaryVtxPosition[1], fPrimaryVtxPosition[2]));
  tFemtoV0->SetYV0(tESDv0->Y());

  AliESDtrack *trackpos = tESDevent->GetTrack(tESDv0->GetPindex()); //AliAODTrack *trackpos = (AliAODTrack*)tESDv0->GetDaughter(0);
  AliESDtrack *trackneg = tESDevent->GetTrack(tESDv0->GetNindex()); //AliAODTrack *trackneg = (AliAODTrack*)tESDv0->GetDaughter(1);

  if(trackpos && trackneg)
    {
      tFemtoV0->SetdcaPosToPrimVertex(TMath::Abs(trackpos->GetD(fPrimaryVtxPosition[0],fPrimaryVtxPosition[1],tESDevent->GetMagneticField())));
      tFemtoV0->SetdcaNegToPrimVertex(TMath::Abs(trackneg->GetD(fPrimaryVtxPosition[0],fPrimaryVtxPosition[1],tESDevent->GetMagneticField())));
      double MomPos[3];
      trackpos->PxPyPz(MomPos);
      tFemtoV0->SetmomPosX(MomPos[0]);
      tFemtoV0->SetmomPosY(MomPos[1]);
      tFemtoV0->SetmomPosZ(MomPos[2]);
      AliFemtoThreeVector mompos(MomPos[0],MomPos[1],MomPos[2]);
      tFemtoV0->SetmomPos(mompos);

      double MomNeg[3];
      trackneg->PxPyPz(MomNeg);
      tFemtoV0->SetmomNegX(MomNeg[0]);
      tFemtoV0->SetmomNegY(MomNeg[1]);
      tFemtoV0->SetmomNegZ(MomNeg[2]);
      AliFemtoThreeVector momneg(MomNeg[0],MomNeg[1],MomNeg[2]);
      tFemtoV0->SetmomNeg(momneg);

      tFemtoV0->SetptPos(trackpos->Pt());
      tFemtoV0->SetptotPos(trackpos->P());
      tFemtoV0->SetptNeg(trackneg->Pt());
      tFemtoV0->SetptotNeg(trackneg->P());
      
      tFemtoV0->SetidNeg(trackneg->GetID());
      //cout<<"tESDv0->GetNegID(): "<<tESDv0->GetNegID()<<endl;
      //cout<<"tFemtoV0->IdNeg(): "<<tFemtoV0->IdNeg()<<endl;
      tFemtoV0->SetidPos(trackpos->GetID());
      
      tFemtoV0->SetEtaPos(trackpos->Eta());
      tFemtoV0->SetEtaNeg(trackneg->Eta());

      tFemtoV0->SetEtaPos(trackpos->Eta()); //tESDv0->PseudoRapPos()
      tFemtoV0->SetEtaNeg(trackneg->Eta()); //tESDv0->PseudoRapNeg()
      tFemtoV0->SetTPCNclsPos(trackpos->GetTPCNcls());
      tFemtoV0->SetTPCNclsNeg(trackneg->GetTPCNcls());
      tFemtoV0->SetTPCclustersPos(trackpos->GetTPCClusterMap());
      tFemtoV0->SetTPCclustersNeg(trackneg->GetTPCClusterMap());
      tFemtoV0->SetTPCsharingPos(trackpos->GetTPCSharedMap());
      tFemtoV0->SetTPCsharingNeg(trackneg->GetTPCSharedMap());
      tFemtoV0->SetNdofPos(trackpos->GetTPCchi2()/trackpos->GetTPCNcls());
      tFemtoV0->SetNdofNeg(trackneg->GetTPCchi2()/trackneg->GetTPCNcls());
      tFemtoV0->SetStatusPos(trackpos->GetStatus());
      tFemtoV0->SetStatusNeg(trackneg->GetStatus());

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

      tFemtoV0->SetTPCMomentumPos(trackpos->GetTPCInnerParam()->P()); //trackpos->GetTPCmomentum();
      tFemtoV0->SetTPCMomentumNeg(trackneg->GetTPCInnerParam()->P()); //trackneg->GetTPCmomentum();

      tFemtoV0->SetdedxPos(trackpos->GetTPCsignal());
      tFemtoV0->SetdedxNeg(trackneg->GetTPCsignal());


      if (fESDpid) {
	tFemtoV0->SetPosNSigmaTPCK(fESDpid->NumberOfSigmasTPC(trackpos,AliPID::kKaon));
	tFemtoV0->SetNegNSigmaTPCK(fESDpid->NumberOfSigmasTPC(trackneg,AliPID::kKaon));
	tFemtoV0->SetPosNSigmaTPCP(fESDpid->NumberOfSigmasTPC(trackpos,AliPID::kProton));
	tFemtoV0->SetNegNSigmaTPCP(fESDpid->NumberOfSigmasTPC(trackneg,AliPID::kProton));
	tFemtoV0->SetPosNSigmaTPCPi(fESDpid->NumberOfSigmasTPC(trackpos,AliPID::kPion));
	tFemtoV0->SetNegNSigmaTPCPi(fESDpid->NumberOfSigmasTPC(trackneg,AliPID::kPion));
      }
      else {
	tFemtoV0->SetPosNSigmaTPCK(-1000);
	tFemtoV0->SetNegNSigmaTPCK(-1000);
	tFemtoV0->SetPosNSigmaTPCP(-1000);
	tFemtoV0->SetNegNSigmaTPCP(-1000);
	tFemtoV0->SetPosNSigmaTPCPi(-1000);
	tFemtoV0->SetNegNSigmaTPCPi(-1000);
      }

      if((tFemtoV0->StatusPos()&AliESDtrack::kTOFpid)==0 || (tFemtoV0->StatusPos()&AliESDtrack::kTIME)==0 || (tFemtoV0->StatusPos()&AliESDtrack::kTOFout)==0 || (tFemtoV0->StatusPos()&AliESDtrack::kTOFmismatch)>0)
	{
	  if((tFemtoV0->StatusNeg()&AliESDtrack::kTOFpid)==0 || (tFemtoV0->StatusNeg()&AliESDtrack::kTIME)==0 || (tFemtoV0->StatusNeg()&AliESDtrack::kTOFout)==0 || (tFemtoV0->StatusNeg()&AliESDtrack::kTOFmismatch)>0)
	    {
	      tFemtoV0->SetPosNSigmaTOFK(-1000);
	      tFemtoV0->SetNegNSigmaTOFK(-1000);
	      tFemtoV0->SetPosNSigmaTOFP(-1000);
	      tFemtoV0->SetNegNSigmaTOFP(-1000);
	      tFemtoV0->SetPosNSigmaTOFPi(-1000);
	      tFemtoV0->SetNegNSigmaTOFPi(-1000);
	    }
	}
      else
	{
	  if (fESDpid) {
	    tFemtoV0->SetPosNSigmaTOFK(fESDpid->NumberOfSigmasTOF(trackpos,AliPID::kKaon));
	    tFemtoV0->SetNegNSigmaTOFK(fESDpid->NumberOfSigmasTOF(trackneg,AliPID::kKaon));
	    tFemtoV0->SetPosNSigmaTOFP(fESDpid->NumberOfSigmasTOF(trackpos,AliPID::kProton));
	    tFemtoV0->SetNegNSigmaTOFP(fESDpid->NumberOfSigmasTOF(trackneg,AliPID::kProton));
	    tFemtoV0->SetPosNSigmaTOFPi(fESDpid->NumberOfSigmasTOF(trackpos,AliPID::kPion));
	    tFemtoV0->SetNegNSigmaTOFPi(fESDpid->NumberOfSigmasTOF(trackneg,AliPID::kPion));
	  }
	  else {
	    tFemtoV0->SetPosNSigmaTOFK(-1000);
	    tFemtoV0->SetNegNSigmaTOFK(-1000);
	    tFemtoV0->SetPosNSigmaTOFP(-1000);
	    tFemtoV0->SetNegNSigmaTOFP(-1000);
	    tFemtoV0->SetPosNSigmaTOFPi(-1000);
	    tFemtoV0->SetNegNSigmaTOFPi(-1000);
	  }
	}

      if ((TMath::Abs(trackpos->GetLabel()) < fStack->GetNtrack()) && (TMath::Abs(trackneg->GetLabel()) < fStack->GetNtrack())) {
  
	//cout<<"tESDv0->GetPdgCode(): "<<tESDv0->GetPdgCode()<<endl;
	// cout<<"Labels: "<<trackpos->GetLabel()<<" "<<trackneg->GetLabel()<<endl;
	  AliFemtoModelHiddenInfo *tInfo = new AliFemtoModelHiddenInfo();
	  //TParticle *tPart = fStack->Particle(tESDv0->GetLabel()); //zle

	  int labelpos = TMath::Abs(trackpos->GetLabel());
	  int labelneg = TMath::Abs(trackneg->GetLabel());
	  TParticle *tPartPos = fStack->Particle(labelpos);
	  TParticle *tPartNeg = fStack->Particle(labelneg);


	  double impactpos[2];
	  double impactneg[2];
	  double covimpact[3];
      
	  AliExternalTrackParam *parampos = new AliExternalTrackParam(*trackpos->GetTPCInnerParam());
	  parampos->PropagateToDCA(fEvent->GetPrimaryVertexTPC(), (fEvent->GetMagneticField()), 10000, impactpos, covimpact);
	  AliExternalTrackParam *paramneg = new AliExternalTrackParam(*trackneg->GetTPCInnerParam());
	  paramneg->PropagateToDCA(fEvent->GetPrimaryVertexTPC(), (fEvent->GetMagneticField()), 10000, impactneg, covimpact);


	  // fillDCA
	  if (TMath::Abs(impactpos[0]) > 0.001) {
	    if (fStack->IsPhysicalPrimary(labelpos)){
	      tFemtoV0->SetImpactDprimPos(impactpos[0]);
	    }
	    else if (fStack->IsSecondaryFromWeakDecay(labelpos)) {
	      tFemtoV0->SetImpactDweakPos(impactpos[0]);
	      //cout << "wea" << endl;
	    }
	    else if (fStack->IsSecondaryFromMaterial(labelpos)) {
	      tFemtoV0->SetImpactDmatPos(impactpos[0]);
	      //cout << "mat" << endl;
	    }
	  }
	  if (TMath::Abs(impactneg[0]) > 0.001) {
	    if (fStack->IsPhysicalPrimary(labelneg)){
	      tFemtoV0->SetImpactDprimNeg(impactneg[0]);
	      //cout << "prim" << endl;
	    }
	    else if (fStack->IsSecondaryFromWeakDecay(labelneg)) {
	      tFemtoV0->SetImpactDweakNeg(impactneg[0]);
	      //cout << "wea" << endl;
	    }
	    else if (fStack->IsSecondaryFromMaterial(labelneg)) {
	      tFemtoV0->SetImpactDmatNeg(impactneg[0]);
	      //cout << "mat" << endl;
	    }

	  }
	  //  end fillDCA 
      

	  //tInfo->SetPDGPid();
	  //tInfo->SetMass();
	  //tInfo->SetTrueMomentum();
	  //tInfo->SetEmissionPoint();

	  // Freeze-out coordinates
	  double fpx=0.0, fpy=0.0, fpz=0.0, fpt=0.0;

	  fpx = tPartPos->Vx() - fPrimaryVtxPosition[0];
	  fpy = tPartPos->Vy() - fPrimaryVtxPosition[1];
	  fpz = tPartPos->Vz() - fPrimaryVtxPosition[2];
	  fpt = tPartPos->T();
	
	  fpx *= 1e13;
	  fpy *= 1e13;
	  fpz *= 1e13;
	  fpt *= 1e13;

	  tInfo->SetPDGPidPos(tPartPos->GetPdgCode());
	  tInfo->SetMassPos(tPartPos->GetMass());
	  tInfo->SetTrueMomentumPos(tPartPos->Px(),tPartPos->Py(),tPartPos->Pz());
	  tInfo->SetEmissionPointPos(fpx,fpy,fpz,fpt);

	  fpx = tPartNeg->Vx() - fPrimaryVtxPosition[0];
	  fpy = tPartNeg->Vy() - fPrimaryVtxPosition[1];
	  fpz = tPartNeg->Vz() - fPrimaryVtxPosition[2];
	  fpt = tPartNeg->T();
	
	  fpx *= 1e13;
	  fpy *= 1e13;
	  fpz *= 1e13;
	  fpt *= 1e13;

	  tInfo->SetPDGPidNeg(tPartNeg->GetPdgCode());
	  tInfo->SetMassNeg(tPartNeg->GetMass());
	  tInfo->SetTrueMomentumNeg(tPartNeg->Px(),tPartNeg->Py(),tPartNeg->Pz());
	  tInfo->SetEmissionPointNeg(fpx,fpy,fpz,fpt);

	  tFemtoV0->SetHiddenInfo(tInfo);
      }
    }
  else
    {
      tFemtoV0->SetStatusPos(999);
      tFemtoV0->SetStatusNeg(999);
    }
  tFemtoV0->SetOnFlyStatusV0(tESDv0->GetOnFlyStatus());
}

void AliFemtoEventReaderESDChainKine::SetMagneticFieldSign(int s)
{
  if(s>0)
    fMagFieldSign = 1;
  else if(s<0)
    fMagFieldSign = -1;
  else
    fMagFieldSign = 0;
}

void AliFemtoEventReaderESDChainKine::GetGlobalPositionAtGlobalRadiiThroughTPC(AliESDtrack *track, Float_t bfield, Float_t globalPositionsAtRadii[9][3])
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

