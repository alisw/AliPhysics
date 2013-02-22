/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// Debug tree task
// the tree is represented as reduced events
// 
// Authors:
//   M.Fasel <M.Fasel@gsi.de>
//
#include <TArrayI.h>
#include <TBits.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TMath.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisUtils.h"
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODPid.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"
#include "AliInputEventHandler.h"
#include "AliHFEcuts.h"
#include "AliHFEextraCuts.h"
#include "AliHFEpidTPC.h"
#include "AliHFEreducedEvent.h"
#include "AliHFEreducedTrack.h"
#include "AliHFEreducedMCParticle.h"
#include "AliHFEsignalCuts.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliPIDResponse.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliVVZERO.h"
#include "AliVZDC.h"
#include "TTreeStream.h"

#include "AliHFEreducedEventCreatorAOD.h"

ClassImp(AliHFEreducedEventCreatorAOD)

AliHFEreducedEventCreatorAOD::AliHFEreducedEventCreatorAOD():
  AliAnalysisTaskSE(),
  fAODMCHeader(NULL),
  fAODArrayMCInfo(NULL),
  fHFEtree(NULL),
  fAnalysisUtils(NULL),
  fHFEevent(NULL),
  fTrackCuts(NULL),
  fExtraCuts(NULL),
  fSignalCuts(NULL),
  fTPCpid(NULL),
  fEventNumber(0),
  fNclustersTPC(70),
  fNclustersTPCPID(0),
  fNclustersITS(2),
  fNbOfTOFSigma(-1.0),
  fRemoveFirstEvent(kFALSE)
{
  //
  // Default constructor
  //
}

AliHFEreducedEventCreatorAOD::AliHFEreducedEventCreatorAOD(const char *name):
  AliAnalysisTaskSE(name),
  fAODMCHeader(NULL),
  fAODArrayMCInfo(NULL),
  fHFEtree(NULL),
  fAnalysisUtils(NULL),
  fHFEevent(NULL),
  fTrackCuts(NULL),
  fExtraCuts(NULL),
  fSignalCuts(NULL),
  fTPCpid(NULL),
  fEventNumber(0),
  fNclustersTPC(70),
  fNclustersTPCPID(0),
  fNclustersITS(2),
  fNbOfTOFSigma(-1.0),
  fRemoveFirstEvent(kFALSE)
{
  //
  // Default constructor
  //
  fTPCpid = new AliHFEpidTPC("QAtpcPID");
  fAnalysisUtils = new AliAnalysisUtils;
  DefineOutput(1, TTree::Class());
}

AliHFEreducedEventCreatorAOD::~AliHFEreducedEventCreatorAOD(){
  //
  // Default destructor
  //
    if(fTPCpid) delete fTPCpid;
    if(fAnalysisUtils) delete fAnalysisUtils;
    if(fHFEevent) delete fHFEevent;
}

void AliHFEreducedEventCreatorAOD::UserCreateOutputObjects(){
  //
  // Create debug tree, signal cuts and track cuts
  //

  //printf("test\n");

 // printf("testa\n");
  fSignalCuts = new AliHFEsignalCuts("HFEsignalCuts", "HFE MC Signal definition");
  //printf("testb\n");
  
  fTrackCuts = new AliHFEcuts("fTrackCuts", "Basic HFE track cuts");
  fTrackCuts->CreateStandardCuts();
  fTrackCuts->SetAOD();
  // Track cuts
  fTrackCuts->SetMinNClustersTPC(fNclustersTPC);
  fTrackCuts->SetMinRatioTPCclusters(0);
  fTrackCuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable); 
  fTrackCuts->SetMinNClustersTPCPID(fNclustersTPCPID);
  fTrackCuts->SetMinNClustersITS(fNclustersITS);
  // Event cuts
  fTrackCuts->SetUseMixedVertex(kTRUE);
  fTrackCuts->SetVertexRange(10.);
  //printf("testa\n");
  fTrackCuts->Initialize();
  //printf("testb\n");

  fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");

  fHFEevent = new AliHFEreducedEvent;
  OpenFile(1);
  fHFEtree = new TTree("HFEtree", "HFE event tree");
  fHFEtree->Branch("HFEevent", "AliHFEreducedEvent", fHFEevent,128000,0);
  PostData(1, fHFEtree);
}

void AliHFEreducedEventCreatorAOD::UserExec(Option_t *){
  //
  // User Exec: Fill debug Tree
  // 

  // Get PID response
  AliPIDResponse *pid = NULL;
  AliInputEventHandler *handler = dynamic_cast<AliInputEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(handler){
    //    printf("testb\n");
    pid = handler->GetPIDResponse();
  } else {
    AliError("No Handler");
  }
  if(!pid){
    //   printf("testc\n");
    AliError("No PID response");
    return;
  }
  if(!fInputEvent) {
    AliError("No Input event");
    return;
  }

  if(fRemoveFirstEvent){
    if(fAnalysisUtils->IsFirstEventInChunk(fInputEvent)) return;
  }

  AliAODTrack copyTrack;

  // MC info
  Bool_t mcthere = kTRUE;
  AliAODEvent *aodE = dynamic_cast<AliAODEvent *>(fInputEvent);
  if(!aodE){
    //  printf("testd\n");
    AliError("No AOD Event");
    return;
  }
  fAODMCHeader = dynamic_cast<AliAODMCHeader *>(fInputEvent->FindListObject(AliAODMCHeader::StdBranchName()));
  if(!fAODMCHeader){ 
      mcthere = kFALSE;
      //   return;
  }
  fAODArrayMCInfo = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if(!fAODArrayMCInfo){ 
      mcthere = kFALSE;
      //  return;
  }
  else {
    fSignalCuts->SetMCAODInfo(fAODArrayMCInfo);
    fTrackCuts->SetMCEvent(aodE);
  }

  // Set for track cuts
  if(mcthere){
    // Check Monte-Carlo events for AODs
    AliDebug(1, Form("Monte-Carlo Event: %p\n", fMCEvent));
    if(fMCEvent){
      AliDebug(1, Form("Available! Number of particles: %d\n", fMCEvent->GetNumberOfTracks()));
    }
  }
  fTrackCuts->SetRecEvent(fInputEvent);

  if(!fTrackCuts->CheckEventCuts("fCutsEvRec", fInputEvent)){
    AliDebug(1, "Event rejected by the event cuts\n");
    return;
  }
  if(!fExtraCuts){
    fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");
  }
  fExtraCuts->SetRecEventInfo(fInputEvent);

 
  // Make Reduced Event 
  //AliHFEreducedEvent hfeevent;
  fHFEevent->~AliHFEreducedEvent();
  new(fHFEevent)AliHFEreducedEvent();

  // Get run number
  fHFEevent->SetRunNumber(fInputEvent->GetRunNumber());

  // Derive trigger 
  AliDebug(1, "Get triggers\n");
  UInt_t trigger = fInputHandler->IsEventSelected();
  if(trigger & AliVEvent::kMB) fHFEevent->SetMBTrigger();
  if(trigger & AliVEvent::kCentral) fHFEevent->SetCentralTrigger();
  if(trigger & AliVEvent::kSemiCentral) fHFEevent->SetCentralTrigger();
  if(trigger & AliVEvent::kEMCEJE) fHFEevent->SetEMCALTrigger();

  // Get Primary Vertex
  AliDebug(1, "Get Primary Vertex\n");
  const AliVVertex *vertex = fInputEvent->GetPrimaryVertex();
  Double_t vtx[3];
  Double_t vcov[6];
  Int_t ncontrib = -1;
  if(vertex) {
    AliDebug(1, "Found vertex\n");
    vertex->GetXYZ(vtx);
    ncontrib = vertex->GetNContributors();
    vertex->GetCovarianceMatrix(vcov);
  }
  fHFEevent->SetVX(vtx[0]);
  fHFEevent->SetVY(vtx[1]);
  fHFEevent->SetVZ(vtx[2]);
  fHFEevent->SetNContribVertex(ncontrib);
  fHFEevent->SetVertexResolution(TMath::Sqrt(TMath::Abs(vcov[5])));
  // Get Primary Vertex from SPD
  const AliVVertex *vertexSPD = aodE->GetPrimaryVertexSPD();
  if(vertexSPD){
    AliDebug(1, "Found SPD vertex\n");
    memset(vtx, 0, sizeof(Double_t) *3);
    vertexSPD->GetXYZ(vtx);
    fHFEevent->SetVXSPD(vtx[0]);
    fHFEevent->SetVYSPD(vtx[1]);
    fHFEevent->SetVZSPD(vtx[2]);
    fHFEevent->SetNContribVertexSPD(vertexSPD->GetNContributors());
    memset(vcov, 0, sizeof(Double_t)*6);
    vertexSPD->GetCovarianceMatrix(vcov);
    AliDebug(1, Form("Covariance Matrix vcov[5] %f\n",vcov[5]));
    fHFEevent->SetVertexResolutionSPD(TMath::Sqrt(TMath::Abs(vcov[5])));
  }

  // Get centrality
  AliDebug(1, "Centrality\n");
  AliCentrality *hicent = fInputEvent->GetCentrality();
  if(hicent) fHFEevent->SetCentrality(
    hicent->GetCentralityPercentile("V0M"),
    hicent->GetCentralityPercentile("V0A"),
    hicent->GetCentralityPercentile("V0C"),
    hicent->GetCentralityPercentile("TKL"),
    hicent->GetCentralityPercentile("TRK"),
    hicent->GetCentralityPercentile("ZNA")
  );

  // Get VZERO Information
  AliDebug(1, "VZERO info\n");
  AliVVZERO *vzeroinfo = fInputEvent->GetVZEROData();
  if(vzeroinfo) fHFEevent->SetV0Multiplicity(vzeroinfo->GetMTotV0A(), vzeroinfo->GetMTotV0C());

  // Get ZDC Information
  AliDebug(1, "ZDC info\n");
  AliVZDC *zdcinfo = fInputEvent->GetZDCData();
  if(zdcinfo) fHFEevent->SetZDCEnergy(zdcinfo->GetZNAEnergy(), zdcinfo->GetZNCEnergy(), zdcinfo->GetZPAEnergy(), zdcinfo->GetZPCEnergy()); 
  
  // Set SPD multiplicity
  AliDebug(1, "SPD multiplicity\n");
  AliAODTracklets *tls = aodE->GetTracklets();
  if(tls) fHFEevent->SetSPDMultiplicity(tls->GetNumberOfTracklets());

  // Look for kink mother
  AliDebug(1, "Vertices\n");
  Int_t numberofvertices = aodE->GetNumberOfVertices();
  Double_t listofmotherkink[numberofvertices];
  Int_t numberofmotherkink = 0;
  //printf("Number of vertices %d\n",numberofvertices);
  for(Int_t ivertex=0; ivertex < numberofvertices; ivertex++) {
    AliAODVertex *aodvertex = aodE->GetVertex(ivertex);
    if(!aodvertex) continue;
    //printf("Type %d\n",aodvertex->GetType());
    if(aodvertex->GetType()==AliAODVertex::kKink) {
      //printf("Find one kink\n");
      AliAODTrack *mother = (AliAODTrack *) aodvertex->GetParent();
      if(!mother) continue;
      Int_t idmother = mother->GetID();
      listofmotherkink[numberofmotherkink] = idmother;
      //printf("ID %d\n",idmother);
      numberofmotherkink++;
    }
  }
  //printf("Number of kink mother in the events %d\n",numberofmotherkink);

  //
  // Loop on MC tracks only
  //
  AliAODMCParticle *mctrack(NULL);
  // Monte-Carlo info
  Int_t source(5);
  if(mcthere){
    AliDebug(1, "Loop MC tracks\n");
    for(Int_t itrack = 0; itrack < fAODArrayMCInfo->GetEntriesFast(); itrack++) {
      mctrack = (AliAODMCParticle *)(fAODArrayMCInfo->At(itrack));
      if(!mctrack) continue;
      AliHFEreducedMCParticle hfemcpart;
      if(fTrackCuts->CheckParticleCuts(static_cast<UInt_t>(AliHFEcuts::kStepMCGenerated), mctrack)) hfemcpart.SetSignal();
      // Kinematics
      hfemcpart.SetSignedPt(mctrack->Pt(), mctrack->Charge() > 0.);
      hfemcpart.SetP(mctrack->P());
      hfemcpart.SetEta(mctrack->Eta());
      hfemcpart.SetPhi(mctrack->Phi());
      hfemcpart.SetPdg(mctrack->GetPdgCode());
      
      // Get Production Vertex in radial direction
      hfemcpart.SetProductionVertex(mctrack->Xv(),mctrack->Yv(),mctrack->Zv());

      // Get Mother PDG code of the particle
      Int_t motherlabel = TMath::Abs(mctrack->GetMother());
      if(motherlabel >= 0 && motherlabel < fAODArrayMCInfo->GetEntriesFast()){
        AliAODMCParticle *mother = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(motherlabel));
        if(mother) hfemcpart.SetMotherPdg(mother->GetPdgCode());
      }
      
      // derive source
      source = 5;
      if(fSignalCuts->IsCharmElectron(mctrack)) source = 0;
      else if(fSignalCuts->IsBeautyElectron(mctrack)) source = 1;
      else if(fSignalCuts->IsGammaElectron(mctrack)) source = 2;
      else if(fSignalCuts->IsNonHFElectron(mctrack)) source = 3;
      else if(TMath::Abs(mctrack->GetPdgCode()) == 11) source = 4;
      else source = 5;
      hfemcpart.SetSource(source);

      fHFEevent->AddMCParticle(&hfemcpart);
    }
  }
  
  //
  // Loop on reconstructed tracks
  //
  TArrayI arraytrack(fInputEvent->GetNumberOfTracks());
  Int_t counterdc=0;
  
  AliAODTrack *track = 0x0;
  AliDebug(1, "Loop reconstructed tracks\n");
  for(Int_t itrack = 0; itrack < fInputEvent->GetNumberOfTracks(); itrack++){
    // Run track loop
    track = dynamic_cast<AliAODTrack *>(fInputEvent->GetTrack(itrack));
    if(!track) continue;
    // Cut track (Only basic track cuts)
    // printf("testv\n");
    if(!fTrackCuts->CheckParticleCuts(AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepRecKineITSTPC, track)) continue;
    //
    //printf("testu\n");

    // Kinematics
    AliHFEreducedTrack hfetrack;
    hfetrack.SetSignedPt(track->Pt(), track->Charge() > 0);
    hfetrack.SetP(track->P());
    hfetrack.SetEta(track->Eta());
    hfetrack.SetPhi(track->Phi());
    hfetrack.SetTPCmomentum(track->GetDetPid() ? track->GetDetPid()->GetTPCmomentum() : track->P());

    // Track ID
    hfetrack.SetTrackID(track->GetID());

    // status
    ULong_t status = track->GetStatus();
    if((status & AliVTrack::kITSrefit) == AliVTrack::kITSrefit) hfetrack.SetITSrefit();
    if((status & AliVTrack::kTPCrefit) == AliVTrack::kTPCrefit) hfetrack.SetTPCrefit();
    if((status & AliVTrack::kTOFpid) == AliVTrack::kTOFpid) hfetrack.SetTOFpid();
    //if((status & AliVTrack::kTOFmismatch) == AliVTrack::kTOFmismatch) hfetrack.SetTOFmismatch();
    if(IsTOFmismatch(track, pid)) hfetrack.SetTOFmismatch(); // New version suggested by Pietro Antonioli
    if(track->IsEMCAL()) hfetrack.SetEMCALpid();
    // Filter
    for(Int_t k=0; k<20; k++) {
      Int_t u = 1<<k;     
      if((track->TestFilterBit(u))) {
	      hfetrack.SetFilterBit(k);
      } 
    }

    if(mcthere){
      // Fill Monte-Carlo Information
      Int_t label = TMath::Abs(track->GetLabel());
      if(label && label < fAODArrayMCInfo->GetEntriesFast())
        mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(label));
      if(!mctrack) continue;
      if(fTrackCuts->CheckParticleCuts(static_cast<UInt_t>(AliHFEcuts::kStepMCGenerated), mctrack)) hfetrack.SetMCSignal();
      // Kinematics
      hfetrack.SetMCSignedPt(mctrack->Pt(),mctrack->Charge() > 0.);
      hfetrack.SetMCP(mctrack->P());
      hfetrack.SetMCEta(mctrack->Eta());
      hfetrack.SetMCPhi(mctrack->Phi());
      hfetrack.SetMCPDG(mctrack->GetPdgCode());
      
      // Get Production Vertex in radial direction
      hfetrack.SetMCProdVtx(mctrack->Xv(),mctrack->Yv(),mctrack->Zv());
      
      // Get Mother PDG code of the particle
      Int_t motherlabel = TMath::Abs(mctrack->GetMother());
      if(motherlabel >= 0 && motherlabel < fAODArrayMCInfo->GetEntriesFast()){
        AliAODMCParticle *mother = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(motherlabel));
        if(mother) hfetrack.SetMCMotherPdg(mother->GetPdgCode());
      }
      
      // derive source
      source = 5;
      if(fSignalCuts->IsCharmElectron(mctrack)) source = 0;
      else if(fSignalCuts->IsBeautyElectron(mctrack)) source = 1;
      else if(fSignalCuts->IsGammaElectron(mctrack)) source = 2;
      else if(fSignalCuts->IsNonHFElectron(mctrack)) source = 3;
      else if(TMath::Abs(mctrack->GetPdgCode()) == 11) source = 4;
      else source = 5;
      hfetrack.SetMCSource(source); 
    }

    // HFE DCA
    Float_t dcaxy = -999.,
            dcaz = -999.;
    fExtraCuts->GetImpactParameters((AliVTrack *)track,dcaxy,dcaz);
    hfetrack.SetDCA(dcaxy, dcaz);

    // Different number of clusters definitions
    Int_t nclustersITS(track->GetITSNcls()),
          nclustersTPC(track->GetTPCNcls()),
          nclustersTPCall(track->GetTPCClusterMap().CountBits()),
          nclustersTPCshared(0);
    UChar_t nfindableTPC = track->GetTPCNclsF();
    const TBits &sharedTPC = track->GetTPCSharedMap();
    for(Int_t ibit = 0; ibit < 160; ibit++) if(sharedTPC.TestBitNumber(ibit)) nclustersTPCshared++;
    hfetrack.SetITSnclusters(nclustersITS);
    hfetrack.SetTPCnclusters(track->GetTPCNcls());
    hfetrack.SetTRDnclusters(track->GetTRDncls());
    hfetrack.SetTPCnclustersPID(track->GetTPCsignalN());
    hfetrack.SetTPCcrossedRows(track->GetTPCNCrossedRows());
    hfetrack.SetTPCnclustersAll(nclustersTPCall);
    hfetrack.SetTPCsharedClusters(nclustersTPCshared);
    hfetrack.SetTPCclusterRatio(nfindableTPC ? static_cast<Float_t>(nclustersTPC)/static_cast<Float_t>(nfindableTPC) : 0);
    hfetrack.SetTPCclusterRatioAll(nfindableTPC ? static_cast<Float_t>(nclustersTPCall)/static_cast<Float_t>(nfindableTPC) : 0);
    UChar_t itsPixel = track->GetITSClusterMap();
    for(int ily = 0; ily < 6; ily++) 
            if(TESTBIT(itsPixel, ily)) hfetrack.SetITScluster(ily);
   
    // TRD related quantities (Yvonne)
    hfetrack.SetTRDntrackletsPID(track->GetTRDntrackletsPID());
    hfetrack.SetTRDnslices(track->GetDetPid()->GetTRDnSlices());
    hfetrack.SetTRDchi2(track->GetTRDchi2());
    AliAODPid* aodpid= track->GetDetPid();
    Double_t* arraytrdsignals;
    arraytrdsignals=aodpid->GetTRDslices();
    Int_t nslicetemp=0;
    for(Int_t iplane = 0; iplane < 6; iplane++){
	    nslicetemp=0;
	    for(Int_t b=(iplane*8);b<((iplane*8)+8);b++){
	      if(track->GetTRDntrackletsPID()>0){
		      if(arraytrdsignals[b]>0.001) nslicetemp++;
	      }
	    }
	    if(nslicetemp > 0) hfetrack.SetTRDstatus(iplane);
    }


    // Kink
    Int_t kink = 0;
    if(fExtraCuts->IsKinkDaughter(track)) kink = 1;

    // kink mother
    Int_t kinkmotherpass = 0;
    for(Int_t kinkmother = 0; kinkmother < numberofmotherkink; kinkmother++) {
      if(track->GetID() == listofmotherkink[kinkmother]) {
	      kinkmotherpass = 1;
	      continue;
      }
    }
    
    // Double counted
    Int_t id(track->GetID());
    for(Int_t l=0; l < counterdc; l++){
      Int_t iTrack2 = arraytrack.At(l);
      if(iTrack2==id){
         hfetrack.SetDoubleCounted();
         break;
      }
    }
    // Add the id at this place
    arraytrack.AddAt(id,counterdc);
    counterdc++;

    // PID
    hfetrack.SetTPCdEdx(track->GetDetPid() ? track->GetDetPid()->GetTPCsignal() : 0.);
    hfetrack.SetTPCsigmaEl(pid->NumberOfSigmasTPC(track, AliPID::kElectron));
    hfetrack.SetTOFsigmaEl(pid->NumberOfSigmasTOF(track, AliPID::kElectron));
    hfetrack.SetTOFmismatchProbability(pid->GetTOFMismatchProbability(track));
    // Eta correction
    copyTrack.~AliAODTrack();
    new(&copyTrack) AliAODTrack(*track);
    if(fTPCpid->HasCentralityCorrection()) fTPCpid->ApplyCentralityCorrection(&copyTrack, static_cast<Double_t>(ncontrib),AliHFEpidObject::kAODanalysis);
    if(fTPCpid->HasEtaCorrection()) fTPCpid->ApplyEtaCorrection(&copyTrack, AliHFEpidObject::kAODanalysis);
    hfetrack.SetTPCsigmaElCorrected(pid->NumberOfSigmasTPC(&copyTrack, AliPID::kElectron));
    hfetrack.SetTPCdEdxCorrected(copyTrack.GetDetPid() ? copyTrack.GetDetPid()->GetTPCsignal() : 0.);
    if(track->IsEMCAL()){
      // EMCAL cluster
      Double_t emcalEnergyOverP = -1.,
               showershape[4] = {0.,0.,0.,0.};
      hfetrack.SetEMCALSigmaEl(pid->NumberOfSigmasEMCAL(track, AliPID::kElectron, emcalEnergyOverP, &showershape[0]));
      hfetrack.SetEMCALEoverP(emcalEnergyOverP);
      hfetrack.SetEMCALShowerShape(showershape);
    }

    // If TOF cut
    if(fNbOfTOFSigma>0.0){
      AliDebug(1, "TOF cut\n");
      if(!((status & AliVTrack::kTOFpid) == AliVTrack::kTOFpid)) continue; 
      if(TMath::Abs(pid->NumberOfSigmasTOF(track, AliPID::kElectron))> fNbOfTOFSigma) continue;
    }

    // Track finished, add NOW to the Event
    fHFEevent->AddTrack(&hfetrack);
    //printf("after\n");
  }
  
  // Fill the debug tree
  AliInfo(Form("Number of tracks: %d\n", fHFEevent->GetNumberOfTracks()));
  AliInfo(Form("Number of MC particles: %d\n", fHFEevent->GetNumberOfMCParticles()));
  fHFEtree->Fill();

  fEventNumber++;
  PostData(1, fHFEtree);
}

void AliHFEreducedEventCreatorAOD::Terminate(Option_t *){
  //
  // Terminate
  //
  AliInfo("terminating...\n");

}

Bool_t AliHFEreducedEventCreatorAOD::IsTOFmismatch(const AliVTrack *const track, const AliPIDResponse *const pid) const {
  //
  // is TOF mismatch
  //
  Double_t probs[AliPID::kSPECIESC];
  AliPIDResponse::EDetPidStatus status = pid->ComputeTOFProbability(track, AliPID::kSPECIESC, probs);
  return status == AliPIDResponse::kDetMismatch;
}

