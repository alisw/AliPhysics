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
//
#include <TArrayI.h>
#include <TBits.h>
#include <TFile.h>
#include <TParticle.h>
#include <TString.h>
#include <TTree.h>

#include "AliAnalysisManager.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisUtils.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliInputEventHandler.h"
#include "AliHFEcuts.h"
#include "AliHFEextraCuts.h"
#include "AliHFEmcQA.h"
#include "AliHFEpidTPC.h"
#include "AliHFEreducedEvent.h"
#include "AliHFEreducedTrack.h"
#include "AliHFEreducedMCParticle.h"
#include "AliHFEsignalCuts.h"
#include "AliHFEV0taginfo.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliMCEventHandler.h"
#include "AliMultiplicity.h"
#include "AliPIDResponse.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliVVZERO.h"
#include "AliVZDC.h"
#include "AliTRDTriggerAnalysis.h"
#include "TTreeStream.h"

#include "AliHFEreducedEventCreatorESD.h"

ClassImp(AliHFEreducedEventCreatorESD)

AliHFEreducedEventCreatorESD::AliHFEreducedEventCreatorESD():
  AliAnalysisTaskSE(),
  fHFEtree(NULL),
  fAnalysisUtils(NULL),
  fHFEevent(NULL),
  fTrackCuts(NULL),
  fExtraCuts(NULL),
  fSignalCuts(NULL),
  fTPCpid(NULL),
  fV0Tagger(NULL),
  fTRDTriggerAnalysis(NULL),
  fEventNumber(0),
  fNclustersTPC(70),
  fNclustersTPCPID(0),
  fNclustersITS(2),
  fRemoveFirstEvent(kFALSE),
  fFlagPileupEvents(kFALSE),
  fSelectSignalOnly(kFALSE)
{
  //
  // Default constructor
  //
}

AliHFEreducedEventCreatorESD::AliHFEreducedEventCreatorESD(const char *name):
  AliAnalysisTaskSE(name),
  fHFEtree(NULL),
  fAnalysisUtils(NULL),
  fHFEevent(NULL),
  fTrackCuts(NULL),
  fExtraCuts(NULL),
  fSignalCuts(NULL),
  fTPCpid(NULL),
  fV0Tagger(NULL),
  fTRDTriggerAnalysis(NULL),
  fEventNumber(0),
  fNclustersTPC(70),
  fNclustersTPCPID(0),
  fNclustersITS(2),
  fRemoveFirstEvent(kFALSE),
  fFlagPileupEvents(kFALSE),
  fSelectSignalOnly(kFALSE)
{
  //
  // Default constructor
  //
  fTPCpid = new AliHFEpidTPC("QAtpcPID");
  fAnalysisUtils = new AliAnalysisUtils;
  fTRDTriggerAnalysis = new AliTRDTriggerAnalysis();
  fTRDTriggerAnalysis->SetRequireMatchElectron(kTRUE);
  fV0Tagger = new AliHFEV0taginfo("Tagger");
  DefineOutput(1, TTree::Class());
}

AliHFEreducedEventCreatorESD::~AliHFEreducedEventCreatorESD(){
  //
  // Default destructor
  //
  if(fAnalysisUtils) delete fAnalysisUtils;
  if(fTPCpid) delete fTPCpid;
  if(fV0Tagger) delete fV0Tagger;
  if(fTRDTriggerAnalysis) delete fTRDTriggerAnalysis;
  if(fHFEevent) delete fHFEevent;
  if(fSignalCuts) delete fSignalCuts;
  if(fTrackCuts) delete fTrackCuts;

}

void AliHFEreducedEventCreatorESD::UserCreateOutputObjects(){
  //
  // Create debug tree, signal cuts and track cuts
  //

  //printf("test\n");

 // printf("testa\n");
  fSignalCuts = new AliHFEsignalCuts("HFEsignalCuts", "HFE MC Signal definition");
  //printf("testb\n");
  
  fTrackCuts = new AliHFEcuts("fTrackCuts", "Basic HFE track cuts");
  fTrackCuts->CreateStandardCuts();
  // Track cuts
  fTrackCuts->SetMaxChi2perClusterTPC(1000.);   // switch off this cut (for the moment);
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

void AliHFEreducedEventCreatorESD::UserExec(Option_t *){
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

  AliDebug(1, "Event Selected");

  AliESDtrack copyTrack;

  // MC info
  Bool_t mcthere = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()) != NULL;
  if(mcthere){ 
    fTrackCuts->SetMCEvent(fMCEvent);
    fSignalCuts->SetMCEvent(fMCEvent);
  }  

  fTrackCuts->SetRecEvent(fInputEvent);

  if(!fTrackCuts->CheckEventCuts("fCutsEvRec", fInputEvent)){
    AliDebug(1, "Event rejected by the event cuts\n");
    return;
  }
  
  // reject pile up in case of pp
  AliESDEvent *event = dynamic_cast<AliESDEvent *>(fInputEvent);
  if(event) {
    TString beamtype = event->GetBeamType();
    //printf("beamtype %s\n",(const char*)beamtype);
    if (strstr(beamtype,"p-p")) {
      //printf("Reject\n");
      if(fInputEvent->IsPileupFromSPD(3, 0.8, 3., 2., 5)){
	      AliDebug(1, "Event flagged as pileup\n");
	      return;
      }
    }
  }

  if(!fExtraCuts){
    fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");
  }
  fExtraCuts->SetRecEventInfo(fInputEvent);

  // Tag all v0s in current event
  if(fV0Tagger){
    fV0Tagger->Reset();
    fV0Tagger->TagV0Tracks(fInputEvent);
  }

  // Make Reduced Event 
  //AliHFEreducedEvent hfeevent;
  fHFEevent->~AliHFEreducedEvent();
  new(fHFEevent)AliHFEreducedEvent();

  // Get run number
  fHFEevent->SetRunNumber(fInputEvent->GetRunNumber());

  // Derive trigger 
  UInt_t trigger = fInputHandler->IsEventSelected();
  if(trigger & AliVEvent::kMB) fHFEevent->SetMBTrigger();
  if((trigger & AliVEvent::kINT7)||(trigger & AliVEvent::kINT8)) fHFEevent->SetINTTrigger();
  if(trigger & AliVEvent::kCentral) fHFEevent->SetCentralTrigger();
  if(trigger & AliVEvent::kSemiCentral) fHFEevent->SetCentralTrigger();
  if(trigger & AliVEvent::kEMCEJE) fHFEevent->SetEMCALTrigger();

  fTRDTriggerAnalysis->CalcTriggers(event);
  if(fTRDTriggerAnalysis->HasTriggeredConfirmed(AliTRDTriggerAnalysis::kHSE)) fHFEevent->SetTRDSETrigger();
  if(fTRDTriggerAnalysis->HasTriggeredConfirmed(AliTRDTriggerAnalysis::kHQU)) fHFEevent->SetTRDDQTrigger();

  // Get Primary Vertex
  const AliVVertex *vertex = fInputEvent->GetPrimaryVertex();
  Double_t vtx[3];
  vertex->GetXYZ(vtx);
  fHFEevent->SetVX(vtx[0]);
  fHFEevent->SetVY(vtx[1]);
  fHFEevent->SetVZ(vtx[2]);
  Int_t ncontrib(vertex->GetNContributors());
  fHFEevent->SetNContribVertex(ncontrib);
  Double_t vcov[6];
  vertex->GetCovarianceMatrix(vcov);
  fHFEevent->SetVertexResolution(TMath::Sqrt(vcov[5]));
  fHFEevent->SetVertexDispersion(static_cast<const AliESDVertex *>(vertex)->GetDispersion());
  // Get Primary Vertex from SPD
  const AliVVertex *vertexSPD = event->GetPrimaryVertexSPD();
  if(vertexSPD){
    memset(vtx, 0, sizeof(Double_t) *3);
    vertexSPD->GetXYZ(vtx);
    fHFEevent->SetVXSPD(vtx[0]);
    fHFEevent->SetVYSPD(vtx[1]);
    fHFEevent->SetVZSPD(vtx[2]);
    fHFEevent->SetNContribVertexSPD(vertexSPD->GetNContributors());
    memset(vcov, 0, sizeof(Double_t)*6);
    vertex->GetCovarianceMatrix(vcov);
    fHFEevent->SetVertexResolutionSPD(TMath::Sqrt(vcov[5]));
    fHFEevent->SetVertexDispersionSPD(static_cast<const AliESDVertex *>(vertex)->GetDispersion());
  }

  // Get centrality
  AliCentrality *hicent = fInputEvent->GetCentrality();
  fHFEevent->SetCentrality(
    hicent->GetCentralityPercentile("V0M"),
    hicent->GetCentralityPercentile("V0A"),
    hicent->GetCentralityPercentile("V0C"),
    hicent->GetCentralityPercentile("TKL"),
    hicent->GetCentralityPercentile("TRK"),
    hicent->GetCentralityPercentile("ZNA"),
    hicent->GetCentralityPercentile("ZNC"),
    hicent->GetCentralityPercentile("CL0"),
    hicent->GetCentralityPercentile("CL1"),
    hicent->GetCentralityPercentile("CND")
  );
  
  // Get VZERO Information
  AliVVZERO *vzeroinfo = fInputEvent->GetVZEROData();
  if(vzeroinfo) fHFEevent->SetV0Multiplicity(vzeroinfo->GetMTotV0A(), vzeroinfo->GetMTotV0C());

  // Get ZDC Information
  AliVZDC *zdcinfo = fInputEvent->GetZDCData();
  if(zdcinfo) fHFEevent->SetZDCEnergy(zdcinfo->GetZNAEnergy(), zdcinfo->GetZNCEnergy(), zdcinfo->GetZPAEnergy(), zdcinfo->GetZPCEnergy()); 

  // Set SPD multiplicity
  const AliMultiplicity *mult = event->GetMultiplicity();
  if(mult) fHFEevent->SetSPDMultiplicity(mult->GetNumberOfTracklets());

  // Flag Pileup
  if(fFlagPileupEvents){
    if(fAnalysisUtils->IsPileUpEvent(fInputEvent)) fHFEevent->SetPileupFlag();
  }

  //
  // Loop on MC tracks only
  //
  AliMCParticle *mctrack(NULL);
  // Monte-Carlo info
  //Int_t source(5);
  if(mcthere){
    const AliVVertex *mcvtx = fMCEvent->GetPrimaryVertex();
    fHFEevent->SetVMC(mcvtx->GetX(), mcvtx->GetY(), mcvtx->GetX());
    for(Int_t itrack = 0; itrack < fMCEvent->GetNumberOfTracks(); itrack++) {
      mctrack = (AliMCParticle *)(fMCEvent->GetTrack(itrack));
      if(!mctrack) continue;
      if(TMath::Abs(mctrack->PdgCode()) != 11) continue;  // in MC truth list store only electrons
      if(fSelectSignalOnly && !fTrackCuts->CheckParticleCuts(static_cast<UInt_t>(AliHFEcuts::kStepMCGenerated), mctrack))  continue;        
      AliHFEreducedMCParticle hfemcpart;
      if(fTrackCuts->CheckParticleCuts(static_cast<UInt_t>(AliHFEcuts::kStepMCGenerated), mctrack)) hfemcpart.SetSignal();
      // Kinematics
      hfemcpart.SetSignedPt(mctrack->Pt(), mctrack->Charge() > 0.);
      hfemcpart.SetP(mctrack->P());
      hfemcpart.SetEta(mctrack->Eta());
      hfemcpart.SetPhi(mctrack->Phi());
      hfemcpart.SetPdg(mctrack->PdgCode());
      
      // Get Production Vertex in radial direction
      hfemcpart.SetProductionVertex(mctrack->Xv(),mctrack->Yv(),mctrack->Zv());

      // Get Mother PDG code of the particle
      Int_t motherlabel = TMath::Abs(mctrack->GetMother());
      if(motherlabel >= 0 && motherlabel < fMCEvent->GetNumberOfTracks()){
        AliMCParticle *mother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(motherlabel));
        if(mother){
          hfemcpart.SetMotherPdg(mother->PdgCode());
          hfemcpart.SetMotherProductionVertex(mother->Xv(),mother->Yv(),mother->Zv());
        }
      }
      
      // derive source
      /*
      source = 7;	      
      if(fSignalCuts->IsCharmElectron(mctrack)) source = 0;
	    else if(fSignalCuts->IsBeautyElectron(mctrack)) source = 1;
	    else if(fSignalCuts->IsGammaElectron(mctrack)) source = 2;
      else if(fSignalCuts->IsNonHFElectron(mctrack)) source = 3;
      else if(fSignalCuts->IsJpsiElectron(mctrack)) source = 4;
      else if(fSignalCuts->IsB2JpsiElectron(mctrack)) source = 5;
      else if(fSignalCuts->IsKe3Electron(mctrack)) source = 6;
	    else source = 7;
      hfemcpart.SetSource(source);
      */
      hfemcpart.SetSource(static_cast<Int_t>(fSignalCuts->GetSignalSource(mctrack)));
      hfemcpart.SetElectronSource(fSignalCuts->GetMCQAObject()->GetElecSource(mctrack));
      fHFEevent->AddMCParticle(&hfemcpart);
    }
  }
  
  //
  // Loop on reconstructed tracks
  //
  TArrayI arraytrack(fInputEvent->GetNumberOfTracks());
  Int_t counterdc=0;
  
  AliESDtrack *track = 0x0;
  for(Int_t itrack = 0; itrack < fInputEvent->GetNumberOfTracks(); itrack++){
    // Run track loop
    track = dynamic_cast<AliESDtrack *>(fInputEvent->GetTrack(itrack));
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
    hfetrack.SetTPCmomentum(track->GetTPCmomentum());

    // Track ID
    hfetrack.SetTrackID(track->GetID());

    // status
    ULong_t status = track->GetStatus();
    if((status & AliVTrack::kITSrefit) == AliVTrack::kITSrefit) hfetrack.SetITSrefit();
    if((status & AliVTrack::kTPCrefit) == AliVTrack::kTPCrefit) hfetrack.SetTPCrefit();
    if((status & AliVTrack::kTOFpid) == AliVTrack::kTOFpid) hfetrack.SetTOFpid();
    //if((status & AliVTrack::kTOFmismatch) == AliVTrack::kTOFmismatch) hfetrack.SetTOFmismatch();
    if(IsTOFmismatch(track, pid)) hfetrack.SetTOFmismatch(); // New version suggested by Pietro Antonioli
    Bool_t isEMCAL(kFALSE);
    Int_t fClsId = track->GetEMCALcluster();
    if(fClsId >= 0) isEMCAL = kTRUE;
    AliDebug(2, Form("cluster ID: %d, EMCAL: %s", fClsId, isEMCAL ? "yes" : "no"));
    if(isEMCAL) hfetrack.SetEMCALpid();
    // no filter bits available for ESDs

    // fill counts of v0-identified particles
    AliPID::EParticleType myv0pid = fV0Tagger->GetV0Info(track->GetID());
    AliHFEreducedTrack::EV0PID_t v0pid = AliHFEreducedTrack::kV0undef;
    if(myv0pid == AliPID::kElectron) v0pid = AliHFEreducedTrack::kV0electron;
    else if(myv0pid == AliPID::kPion) v0pid = AliHFEreducedTrack::kV0pion;
    else if(myv0pid == AliPID::kProton) v0pid = AliHFEreducedTrack::kV0proton;
    hfetrack.SetV0PID(v0pid);

    Double_t v0prodR = fV0Tagger->GetV0ProdR(track->GetID());
    hfetrack.SetV0prodR(v0prodR);

    if(mcthere){
      // Fill Monte-Carlo Information
      Int_t label = TMath::Abs(track->GetLabel());
      if(label < fMCEvent->GetNumberOfTracks())
        mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(label));
      if(mctrack){
        AliDebug(2, "Associated MC particle found");
        if(fTrackCuts->CheckParticleCuts(static_cast<UInt_t>(AliHFEcuts::kStepMCGenerated), mctrack)) hfetrack.SetMCSignal();
        // Kinematics
        hfetrack.SetMCSignedPt(mctrack->Pt(),mctrack->Charge() > 0.);
        hfetrack.SetMCP(mctrack->P());
        hfetrack.SetMCEta(mctrack->Eta());
        hfetrack.SetMCPhi(mctrack->Phi());
        hfetrack.SetMCPDG(mctrack->PdgCode());
      
        // Get Production Vertex in radial direction
        hfetrack.SetMCProdVtx(mctrack->Xv(),mctrack->Yv(),mctrack->Zv());
      
        // Get Mother PDG code of the particle
        Int_t motherlabel = TMath::Abs(mctrack->GetMother());
        if(motherlabel >= 0 && motherlabel < fMCEvent->GetNumberOfTracks()){
          AliMCParticle *mother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(motherlabel));
          if(mother){ 
            hfetrack.SetMCMotherPdg(mother->PdgCode());
            hfetrack.SetMCMotherProdVtx(mother->Xv(),mother->Yv(),mother->Zv());
          }
        }
        
        // derive source
        /*
        source = 7;
        if(fSignalCuts->IsCharmElectron(track)) source = 0;
	      else if(fSignalCuts->IsBeautyElectron(track)) source = 1;
	      else if(fSignalCuts->IsGammaElectron(track)) source = 2;
        else if(fSignalCuts->IsNonHFElectron(track)) source = 3;
        else if(fSignalCuts->IsJpsiElectron(track)) source = 4;
        else if(fSignalCuts->IsB2JpsiElectron(track)) source = 5;
        else if(fSignalCuts->IsKe3Electron(track)) source = 6;
        else source = 7;
        hfetrack.SetMCSource(source); 
        */
        hfetrack.SetMCSource(static_cast<Int_t>(fSignalCuts->GetSignalSource(track))); 
        hfetrack.SetMCElectronSource(fSignalCuts->GetMCQAObject()->GetElecSource(mctrack));
      } else {
        AliDebug(2, "Associated MC particle not found");
      }
    }

    // HFE DCA
    Float_t dcaxy = -999.,
            dcaz = -999.;
    Double_t dcaErr, dcaxyD;
    fExtraCuts->GetImpactParameters((AliVTrack *)track,dcaxy,dcaz);
    fExtraCuts->GetHFEImpactParameters((AliVTrack *)track,dcaxyD,dcaErr);
    hfetrack.SetDCA(dcaxyD, dcaz);
    hfetrack.SetDCAerr(dcaErr);
    Double_t hfeImpactParam(-999.), hfeImpactParamResol(-999.);
    fExtraCuts->GetHFEImpactParameters((AliVTrack *)track,hfeImpactParam,hfeImpactParamResol);
    hfetrack.SetHFEImpactParam(hfeImpactParam,hfeImpactParamResol);

    // Different number of clusters definitions
    Int_t nclustersITS(track->GetITSclusters(NULL)),
          nclustersTPC(track->GetTPCNcls()),
          nclustersTPCall(track->GetTPCClusterMap().CountBits()),
          nclustersTPCshared(0);
    UChar_t nfindableTPC = track->GetTPCNclsF();
    const TBits &sharedTPC = track->GetTPCSharedMap();
    for(Int_t ibit = 0; ibit < 160; ibit++) if(sharedTPC.TestBitNumber(ibit)) nclustersTPCshared++;
    hfetrack.SetChi2PerTPCcluster(track->GetTPCchi2()/Double_t(nclustersTPC));
    hfetrack.SetITSnclusters(nclustersITS);
    hfetrack.SetTPCnclusters(nclustersTPC);
    hfetrack.SetTRDnclusters(track->GetTRDncls());
    hfetrack.SetTPCnclustersPID(track->GetTPCsignalN());
    hfetrack.SetTPCcrossedRows(track->GetTPCCrossedRows());
    hfetrack.SetTPCnclustersAll(nclustersTPCall);
    hfetrack.SetTPCsharedClusters(nclustersTPCshared);
    hfetrack.SetTPCclusterRatio(nfindableTPC ? static_cast<Float_t>(nclustersTPC)/static_cast<Float_t>(nfindableTPC) : 0);
    hfetrack.SetTPCclusterRatioAll(nfindableTPC ? static_cast<Float_t>(nclustersTPCall)/static_cast<Float_t>(nfindableTPC) : 0);
    UChar_t itsPixel = track->GetITSClusterMap();
    for(int ily = 0; ily < 6; ily++) 
            if(TESTBIT(itsPixel, ily)) hfetrack.SetITScluster(ily);
   
    // TRD related quantities (Yvonne)
    Int_t nslices = track->GetNumberOfTRDslices();
    hfetrack.SetTRDntrackletsPID(track->GetTRDntrackletsPID());
    hfetrack.SetTRDnslices(nslices);
    hfetrack.SetTRDchi2(track->GetTRDchi2());
    Int_t nslicetemp=0;
    for(Int_t iplane = 0; iplane < 6; iplane++){
	    nslicetemp=0;
	    for(Int_t isl = 0; isl < nslices; isl++){
	      if(track->GetTRDntrackletsPID()>0){
		      if(track->GetTRDslice(iplane, isl)>0.001) nslicetemp++;
	      }
	    }
	    if(nslicetemp > 0) hfetrack.SetTRDstatus(iplane);
    }


    //test for kink tracks
    if(fExtraCuts->IsKinkMother(track)) hfetrack.SetIsKinkMother();
    else if(fExtraCuts->IsKinkDaughter(track)) hfetrack.SetIsKinkDaughter();
    
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
    hfetrack.SetTPCdEdx(track->GetTPCsignal());
    hfetrack.SetTPCsigmaEl(pid->NumberOfSigmasTPC(track, AliPID::kElectron));
    hfetrack.SetTOFsigmaEl(pid->NumberOfSigmasTOF(track, AliPID::kElectron));
    if(TMath::Abs(pid->NumberOfSigmasTOF(track, AliPID::kDeuteron)) < 40.) hfetrack.SetTOFsigmaDeuteron(pid->NumberOfSigmasTOF(track, AliPID::kDeuteron));   else hfetrack.SetTOFsigmaDeuteron(100);
    hfetrack.SetTOFmismatchProbability(pid->GetTOFMismatchProbability(track));
    hfetrack.SetITSsigmaEl(pid->NumberOfSigmasITS(track, AliPID::kElectron));
    // Eta correction
    copyTrack.~AliESDtrack();
    new(&copyTrack) AliESDtrack(*track);
    if(fTPCpid->HasCentralityCorrection()) fTPCpid->ApplyCentralityCorrection(&copyTrack, static_cast<Double_t>(ncontrib),AliHFEpidObject::kESDanalysis);
    if(fTPCpid->HasEtaCorrection()) fTPCpid->ApplyEtaCorrection(&copyTrack, AliHFEpidObject::kESDanalysis);
    hfetrack.SetTPCsigmaElCorrected(pid->NumberOfSigmasTPC(&copyTrack, AliPID::kElectron));
    hfetrack.SetTPCdEdxCorrected(copyTrack.GetTPCsignal());
    if(isEMCAL){
      AliDebug(2, "Adding EMCAL PID information");
      // EMCAL cluster
      Double_t emcalEnergyOverP = -1.,
               showershape[4] = {0.,0.,0.,0.};
      hfetrack.SetEMCALSigmaEl(pid->NumberOfSigmasEMCAL(track, AliPID::kElectron, emcalEnergyOverP, &showershape[0]));
      hfetrack.SetEMCALEoverP(emcalEnergyOverP);
      hfetrack.SetEMCALShowerShape(showershape);
    }

    // Track finished, add NOW to the Event
    fHFEevent->AddTrack(&hfetrack);
    //printf("after\n");
  }
  
  // Fill the debug tree
  //AliInfo(Form("Number of tracks: %d\n", fHFEevent->GetNumberOfTracks()));
  //AliInfo(Form("Number of MC particles: %d\n", fHFEevent->GetNumberOfMCParticles()));
  fHFEtree->Fill();

  fEventNumber++;
  PostData(1, fHFEtree);
}

void AliHFEreducedEventCreatorESD::Terminate(Option_t *){
  //
  // Terminate
  //
  AliInfo("terminating...\n");

}

Bool_t AliHFEreducedEventCreatorESD::IsTOFmismatch(const AliVTrack *const track, const AliPIDResponse *const pid) const {
  //
  // Is TOF mismatch
  //
  Double_t probs[AliPID::kSPECIESC];
  AliPIDResponse::EDetPidStatus status = pid->ComputeTOFProbability(track, AliPID::kSPECIESC, probs);
  return status == AliPIDResponse::kDetMismatch;
}

