////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoEventReaderStandard - the reader class for the Alice ESD, AOD      //
// the model Kinematics information tailored for the Task framework           //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoEventReaderStandard.h"

#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TBits.h"
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
#include "AliGenCocktailEventHeader.h"

#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"

#include "AliVertexerTracks.h"
#include "assert.h"

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoEventReaderStandard);
  /// \endcond
#endif

#if !(ST_NO_NAMESPACES)
  using namespace units;
#endif

using namespace std;
//____________________________
AliFemtoEventReaderStandard::AliFemtoEventReaderStandard():
  AliFemtoEventReader(),
  fFileName(" "),
  fNumberofEvent(0),
  fCurEvent(0),
  fCurFile(0),
  fESDEvent(0x0),
  fAODEvent(0x0),
  fStack(0x0),
  fGenHeader(0x0),
  fInputType(kUnknown),
  fUsePhysicsSel(kFALSE),
  fSelect(0),
  fTrackCuts(0x0),
  fUseTPCOnly(kFALSE)
{
  //constructor with 0 parameters , look at default settings
}

//__________________
AliFemtoEventReaderStandard::AliFemtoEventReaderStandard(const AliFemtoEventReaderStandard& aReader):
  AliFemtoEventReader(aReader),
  fFileName(" "),
  fNumberofEvent(0),
  fCurEvent(0),
  fCurFile(0),
  fESDEvent(0x0),
  fAODEvent(0x0),
  fStack(0x0),
  fGenHeader(0x0),
  fInputType(kUnknown),
  fUsePhysicsSel(kFALSE),
  fSelect(0),
  fTrackCuts(0x0),
  fUseTPCOnly(kFALSE)
{
  // Copy constructor
  fCurEvent = aReader.fCurEvent;
  fCurFile = aReader.fCurFile;
  fESDEvent = new AliESDEvent();
  fAODEvent = new AliAODEvent();
  fStack = aReader.fStack;
  fInputType = aReader.fInputType;
  fUsePhysicsSel = aReader.fUsePhysicsSel;
  if (fUsePhysicsSel) fSelect = new AliPhysicsSelection();
  fTrackCuts = new AliESDtrackCuts(*(aReader.fTrackCuts));
  fUseTPCOnly = aReader.fUseTPCOnly;
}
//__________________
AliFemtoEventReaderStandard::~AliFemtoEventReaderStandard()
{
  //Destructor
  delete fESDEvent;
  delete fAODEvent;
  if (fTrackCuts) delete fTrackCuts;
}

//__________________
AliFemtoEventReaderStandard& AliFemtoEventReaderStandard::operator=(const AliFemtoEventReaderStandard& aReader)
{
  // Assignment operator
  if (this == &aReader)
    return *this;

  fCurEvent = aReader.fCurEvent;
  fCurFile = aReader.fCurFile;
  if (fESDEvent) delete fESDEvent;
  fESDEvent = new AliESDEvent();
  if (fAODEvent) delete fAODEvent;
  fAODEvent = new AliAODEvent();
  fStack = aReader.fStack;
  fGenHeader = aReader.fGenHeader;
  fInputType = aReader.fInputType;
  fUsePhysicsSel = aReader.fUsePhysicsSel;
  if (fUsePhysicsSel) fSelect = new AliPhysicsSelection();
  if (fTrackCuts) delete fTrackCuts;
  fTrackCuts = new AliESDtrackCuts(*(aReader.fTrackCuts));
  fUseTPCOnly = aReader.fUseTPCOnly;

  return *this;
}
//__________________
// Simple report
AliFemtoString AliFemtoEventReaderStandard::Report()
{
  AliFemtoString temp = "\n This is the AliFemtoEventReaderStandard\n";
  return temp;
}
//__________________
AliFemtoEvent* AliFemtoEventReaderStandard::ReturnHbtEvent()
{
  // Get the event, read all the relevant information
  // and fill the AliFemtoEvent class
  // Returns a valid AliFemtoEvent
  AliFemtoEvent *hbtEvent = 0;
  string tFriendFileName;

  hbtEvent = new AliFemtoEvent;

  // Get the friend information
  cout<<"starting to read event "<<fCurEvent<<endl;
  if ((fInputType == kESD) || (fInputType == kESDKine)) {
    if(fESDEvent->GetAliESDOld())fESDEvent->CopyFromOldESD();

    if (fUsePhysicsSel) {
      hbtEvent->SetIsCollisionCandidate(fSelect->IsCollisionCandidate(fESDEvent));
      if (!(fSelect->IsCollisionCandidate(fESDEvent)))
	printf("Event not a collision candidate\n");
    }
    else
      hbtEvent->SetIsCollisionCandidate(kTRUE);
  }
  else {
    hbtEvent->SetIsCollisionCandidate(kTRUE);
  }

  double fV1[3];

  //setting basic things
  if ((fInputType == kESD) || (fInputType == kESDKine)) {
    hbtEvent->SetRunNumber(fESDEvent->GetRunNumber());
    hbtEvent->SetMagneticField(fESDEvent->GetMagneticField()*kilogauss);//to check if here is ok
    hbtEvent->SetZDCN1Energy(fESDEvent->GetZDCN1Energy());
    hbtEvent->SetZDCP1Energy(fESDEvent->GetZDCP1Energy());
    hbtEvent->SetZDCN2Energy(fESDEvent->GetZDCN2Energy());
    hbtEvent->SetZDCP2Energy(fESDEvent->GetZDCP2Energy());
    hbtEvent->SetZDCEMEnergy(fESDEvent->GetZDCEMEnergy());
    hbtEvent->SetZDCParticipants(fESDEvent->GetZDCParticipants());
    hbtEvent->SetTriggerMask(fESDEvent->GetTriggerMask());
    hbtEvent->SetTriggerCluster(fESDEvent->GetTriggerCluster());

    printf("Got event type %i\n", fESDEvent->GetEventType());

    //Vertex
    double fVCov[6];
    //   if (fUseTPCOnly) {
    //     fESDEvent->GetPrimaryVertexTPC()->GetXYZ(fV1);
    //     fESDEvent->GetPrimaryVertexTPC()->GetCovMatrix(fVCov);
    //     if (!fESDEvent->GetPrimaryVertexTPC()->GetStatus())
    //       fVCov[4] = -1001.0;
    //   }
    //   else {
    if (fESDEvent->GetPrimaryVertex()) {
      fESDEvent->GetPrimaryVertex()->GetXYZ(fV1);
      fESDEvent->GetPrimaryVertex()->GetCovMatrix(fVCov);

      if (!fESDEvent->GetPrimaryVertex()->GetStatus()) {
	// Get the vertex from SPD
	fESDEvent->GetPrimaryVertexSPD()->GetXYZ(fV1);
	fESDEvent->GetPrimaryVertexSPD()->GetCovMatrix(fVCov);


	if (!fESDEvent->GetPrimaryVertexSPD()->GetStatus())
	  fVCov[4] = -1001.0;
	else {
	  fESDEvent->GetPrimaryVertexSPD()->GetXYZ(fV1);
	  fESDEvent->GetPrimaryVertexSPD()->GetCovMatrix(fVCov);
	}
      }
    }
    else {
      if (fESDEvent->GetPrimaryVertexSPD()) {
	fESDEvent->GetPrimaryVertexSPD()->GetXYZ(fV1);
	fESDEvent->GetPrimaryVertexSPD()->GetCovMatrix(fVCov);
      }
    }
    if ((!fESDEvent->GetPrimaryVertex()) && (!fESDEvent->GetPrimaryVertexSPD()))
      {
	cout << "No vertex found !!!" << endl;
	fV1[0] = 10000.0;
	fV1[1] = 10000.0;
	fV1[2] = 10000.0;
	fVCov[4] = -1001.0;
      }

    AliFmThreeVectorF vertex(fV1[0],fV1[1],fV1[2]);

    hbtEvent->SetPrimVertPos(vertex);
    hbtEvent->SetPrimVertCov(fVCov);
  }

  if ((fInputType == kAOD) || (fInputType == kAODKine)) {
    hbtEvent->SetRunNumber(fAODEvent->GetRunNumber());
    hbtEvent->SetMagneticField(fAODEvent->GetMagneticField()*kilogauss);//to check if here is ok
    hbtEvent->SetZDCN1Energy(fAODEvent->GetZDCN1Energy());
    hbtEvent->SetZDCP1Energy(fAODEvent->GetZDCP1Energy());
    hbtEvent->SetZDCN2Energy(fAODEvent->GetZDCN2Energy());
    hbtEvent->SetZDCP2Energy(fAODEvent->GetZDCP2Energy());
    hbtEvent->SetZDCEMEnergy(fAODEvent->GetZDCEMEnergy(0));
    hbtEvent->SetZDCParticipants(0);
    hbtEvent->SetTriggerMask(fAODEvent->GetTriggerMask());
    hbtEvent->SetTriggerCluster(fAODEvent->GetTriggerCluster());

    // Primary Vertex position
    fAODEvent->GetPrimaryVertex()->GetPosition(fV1);

    AliFmThreeVectorF vertex(fV1[0],fV1[1],fV1[2]);
    hbtEvent->SetPrimVertPos(vertex);
  }

  if ((fInputType == kESDKine) || (fInputType == kAODKine)) {
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
  }

  //starting to reading tracks
  int nofTracks=0;  //number of reconstructed tracks in event
  if ((fInputType == kESD) || (fInputType == kESDKine))
    nofTracks = fESDEvent->GetNumberOfTracks();
  else if ((fInputType == kAOD) || (fInputType == kAODKine))
    nofTracks = fAODEvent->GetNumberOfTracks();

  int realnofTracks=0;//number of track which we use ina analysis

  TClonesArray *mcP = 0;
  Int_t *motherids=0;

  if (fInputType == kAODKine) {
    // Attempt to access MC header
    AliAODMCHeader *mcH;
    mcH = (AliAODMCHeader *) fAODEvent->FindListObject(AliAODMCHeader::StdBranchName());
    if (!mcH) {
      cout << "AOD MC information requested, but no header found!" << endl;
    }

    mcP = (TClonesArray *) fAODEvent->FindListObject(AliAODMCParticle::StdBranchName());
    if (!mcP) {
      cout << "AOD MC information requested, but no particle array found!" << endl;
    }
    AliAODHeader * header = dynamic_cast<AliAODHeader*>(fAODEvent->GetHeader());
    assert(header&&"Not a standard AOD");

    hbtEvent->SetReactionPlaneAngle(header->GetQTheta(0)/2.0);

    if (mcP) {
      motherids = new Int_t[((AliAODMCParticle *) mcP->At(mcP->GetEntries()-1))->GetLabel()];
      for (int ip=0; ip<mcP->GetEntries(); ip++) motherids[ip] = 0;

      // Read in mother ids
      AliAODMCParticle *motherpart;
      for (int ip=0; ip<mcP->GetEntries(); ip++) {
	motherpart = (AliAODMCParticle *) mcP->At(ip);
	if (motherpart->GetDaughterLabel(0) > 0)
	  motherids[motherpart->GetDaughterLabel(0)] = ip;
	if (motherpart->GetDaughterLabel(1) > 0)
	  motherids[motherpart->GetDaughterLabel(1)] = ip;
      }
    }
  }

  if (fInputType == kESDKine) {
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
    }
  }

  for (int i=0;i<nofTracks;i++)
    {
      //      cout << "Reading track " << i << endl;
      bool  tGoodMomentum=true; //flaga to chcek if we can read momentum of this track

      AliFemtoTrack* trackCopy = new AliFemtoTrack();

      if ((fInputType == kESD) || (fInputType == kESDKine)) {

	AliESDtrack *esdtrack = 0x0;
	if (fUseTPCOnly) {
	  AliESDtrack *mcp = fESDEvent->GetTrack(i);
	  esdtrack = AliESDtrackCuts::GetTPCOnlyTrack(fESDEvent, mcp->GetID());
	  //	  printf("Got %p for track %i | ", esdtrack, mcp->GetID());
	}
	else {
	  esdtrack = fESDEvent->GetTrack(i);//getting next track
	}

	if (esdtrack && (fTrackCuts->AcceptTrack(esdtrack))) {

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
	  double impact[2];
	  double covimpact[3];

	  // 	  if (fUseTPCOnly) {
	  // 	    if (!esdtrack->GetTPCInnerParam()) {
	  // 	      cout << "No TPC inner param !" << endl;
	  // 	      delete trackCopy;
	  // 	      continue;
	  // 	    }

	  // 	    AliExternalTrackParam *param = new AliExternalTrackParam(*esdtrack->GetTPCInnerParam());
	  // 	    param->GetXYZ(rxyz);
	  // 	    param->PropagateToDCA(fESDEvent->GetPrimaryVertexTPC(), (fESDEvent->GetMagneticField()), 10000, impact, covimpact);
	  // 	    param->GetPxPyPz(pxyz);//reading noconstarined momentum

	  // 	    if (fRotateToEventPlane) {
	  // 	      double tPhi = TMath::ATan2(pxyz[1], pxyz[0]);
	  // 	      double tRad = TMath::Hypot(pxyz[0], pxyz[1]);

	  // 	      pxyz[0] = tRad*TMath::Cos(tPhi - tReactionPlane);
	  // 	      pxyz[1] = tRad*TMath::Sin(tPhi - tReactionPlane);
	  // 	    }

	  // 	    AliFemtoThreeVector v(pxyz[0],pxyz[1],pxyz[2]);
	  // 	    if (v.mag() < 0.0001) {
	  // 	      //	  cout << "Found 0 momentum ???? " << pxyz[0] << " " << pxyz[1] << " " << pxyz[2] << endl;
	  // 	      delete trackCopy;
	  // 	      continue;
	  // 	    }

	  // 	    trackCopy->SetP(v);//setting momentum
	  // 	    trackCopy->SetPt(sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]));

	  // 	    const AliFmThreeVectorD kP(pxyz[0],pxyz[1],pxyz[2]);
	  // 	    const AliFmThreeVectorD kOrigin(fV1[0],fV1[1],fV1[2]);
	  // 	    //setting helix I do not if it is ok
	  // 	    AliFmPhysicalHelixD helix(kP,kOrigin,(double)(fESDEvent->GetMagneticField())*kilogauss,(double)(trackCopy->Charge()));
	  // 	    trackCopy->SetHelix(helix);

	  // 	    //some stuff which could be useful
	  // 	    trackCopy->SetImpactD(impact[0]);
	  // 	    trackCopy->SetImpactZ(impact[1]);
	  // 	    trackCopy->SetCdd(covimpact[0]);
	  // 	    trackCopy->SetCdz(covimpact[1]);
	  // 	    trackCopy->SetCzz(covimpact[2]);
	  // 	    trackCopy->SetSigmaToVertex(GetSigmaToVertex(impact, covimpact));

	  // 	    delete param;
	  // 	  }
	  // 	  else {
	  if (fUseTPCOnly)
	    tGoodMomentum=esdtrack->GetPxPyPz(pxyz);
	  else
	    tGoodMomentum=esdtrack->GetConstrainedPxPyPz(pxyz); //reading constrained momentum
	  //	  printf("Got good momentum %i\n", tGoodMomentum);

	  AliFemtoThreeVector v(pxyz[0],pxyz[1],pxyz[2]);
	  if (v.Mag() < 0.0001) {

	    delete trackCopy;
	    continue;
	  }

	  trackCopy->SetP(v);//setting momentum
	  trackCopy->SetPt(sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]));
	  const AliFmThreeVectorD kP(pxyz[0],pxyz[1],pxyz[2]);
	  const AliFmThreeVectorD kOrigin(fV1[0],fV1[1],fV1[2]);
	  //setting helix I do not if it is ok
	  AliFmPhysicalHelixD helix(kP,kOrigin,(double)(fESDEvent->GetMagneticField())*kilogauss,(double)(trackCopy->Charge()));
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
	  trackCopy->SetSigmaToVertex(AliESDtrackCuts::GetSigmaToVertex(esdtrack));

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
	  if (fInputType == kESDKine) {
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

	      if (motherids[TMath::Abs(esdtrack->GetLabel())]>0) {
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
	  }
	  //      cout << "Got freeze-out " << fpx << " " << fpy << " " << fpz << " " << fpt << " " <<  mass2 << " " << tPart->GetPdgCode() << endl;
	}
	else
	  tGoodMomentum = false;

	if (fUseTPCOnly)
	  if (esdtrack) delete esdtrack;
      }

      if ((fInputType == kAOD) || (fInputType == kAODKine)) {
	// Read in the normal AliAODTracks
	const AliAODTrack *aodtrack=dynamic_cast<const AliAODTrack*>(fAODEvent->GetTrack(i));
        assert(aodtrack&&"Not a standard AOD"); // getting the AODtrack directly

	// 	if (!aodtrack->TestFilterBit(fFilterBit))
	// 	  continue;

	CopyAODtoFemtoTrack(aodtrack, trackCopy);

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

  if (motherids)
    delete [] motherids;

  hbtEvent->SetNumberOfTracks(realnofTracks);//setting number of track which we read in event
  fCurEvent++;
  //  cout<<"end of reading nt "<<nofTracks<<" real number "<<realnofTracks<<endl;
  return hbtEvent;
}
//___________________
void AliFemtoEventReaderStandard::SetESDSource(AliESDEvent *aESD)
{
  // The chain loads the ESD for us
  // You must provide the address where it can be found
  fESDEvent = aESD;
}
//___________________
void AliFemtoEventReaderStandard::SetAODSource(AliAODEvent *aAOD)
{
  // The chain loads the ESD for us
  // You must provide the address where it can be found
  fAODEvent = aAOD;
}
//___________________
void AliFemtoEventReaderStandard::SetStackSource(AliStack *aStack)
{
  // The chain loads the stack for us
  // You must provide the address where it can be found
  fStack = aStack;
}
void AliFemtoEventReaderStandard::SetInputType(AliFemtoInputType aInput)
{
  // Set the proper input type
  fInputType = aInput;
}
//___________________
void AliFemtoEventReaderStandard::SetGenEventHeader(AliGenEventHeader *aGenHeader)
{
  // The chain loads the generator event header for us
  // You must provide the address where it can be found
  fGenHeader = aGenHeader;
}

void AliFemtoEventReaderStandard::SetUsePhysicsSelection(const bool usephysics)
{
  fUsePhysicsSel = usephysics;
  if (!fSelect) fSelect = new AliPhysicsSelection();
}

void AliFemtoEventReaderStandard::SetESDTrackCuts(AliESDtrackCuts *esdcuts)
{
  // Set external ESD track cuts
  fTrackCuts = esdcuts;
}

void AliFemtoEventReaderStandard::SetUseTPCOnly(const bool usetpconly)
{
  // Set flag to use TPC only tracks
  fUseTPCOnly = usetpconly;
}

void AliFemtoEventReaderStandard::CopyAODtoFemtoTrack(const AliAODTrack *tAodTrack,
						      AliFemtoTrack *tFemtoTrack)
{
  // Copy the track information from the AOD into the internal AliFemtoTrack
  // If it exists, use the additional information from the PWG2 AOD

  // Primary Vertex position
  double fV1[3];
  fAODEvent->GetPrimaryVertex()->GetPosition(fV1);

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
  AliFmPhysicalHelixD helix(ktP,kOrigin,(double)(fAODEvent->GetMagneticField())*kilogauss,(double)(tFemtoTrack->Charge()));
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
  tFemtoTrack->SetTPCnclsF(1);
  tFemtoTrack->SetTPCsignalN(1);
  tFemtoTrack->SetTPCsignalS(1);

  TBits tAllTrue;
  TBits tAllFalse;
  tAllTrue.ResetAllBits(kTRUE);
  tAllFalse.ResetAllBits(kFALSE);


  // If not use dummy values
  tFemtoTrack->SetTPCClusterMap(tAllTrue);
  tFemtoTrack->SetTPCSharedMap(tAllFalse);

  double xtpc[3] = {0,0,0};
  tFemtoTrack->SetNominalTPCEntrancePoint(xtpc);
  tFemtoTrack->SetNominalTPCExitPoint(xtpc);

  int indexes[3];
  for (int ik=0; ik<3; ik++) {
    indexes[ik] = 0;
  }
  tFemtoTrack->SetKinkIndexes(indexes);
}

AliAODMCParticle* AliFemtoEventReaderStandard::GetParticleWithLabel(TClonesArray *mcP, Int_t aLabel)
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
