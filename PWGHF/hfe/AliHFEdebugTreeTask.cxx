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
// 
// Authors:
//   Markus Fasel <M.Fasel@gsi.de>
//
#include <TBits.h>
#include <TString.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliHFEcuts.h"
#include "AliHFEsignalCuts.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliPIDResponse.h"
#include "AliTrackReference.h"
#include "AliVEvent.h"
#include "AliHFEpidTPC.h"
#include "AliHFEpidTRD.h"
#include "AliHFEmcQA.h"
#include "TTreeStream.h"

#include "AliHFEdebugTreeTask.h"

ClassImp(AliHFEdebugTreeTask)

AliHFEdebugTreeTask::AliHFEdebugTreeTask():
  AliAnalysisTaskSE(),
  fTrackCuts(NULL),
  fSignalCuts(NULL),
  fTRDpid(NULL),
  fTPCpid(NULL),
  fExtraCuts(NULL),
  fNclustersTPC(70),
  fNclustersTPCPID(0),
  fNclustersITS(2),
  fFilename("HFEtree.root"),
  fDebugTree(NULL),
  fNparents(-1)
{

}

AliHFEdebugTreeTask::AliHFEdebugTreeTask(const char *name):
  AliAnalysisTaskSE(name),
  fTrackCuts(NULL),
  fSignalCuts(NULL),
  fTRDpid(NULL),
  fTPCpid(NULL),
  fExtraCuts(NULL),
  fNclustersTPC(70),
  fNclustersTPCPID(0),
  fNclustersITS(2),
  fFilename("HFEtree.root"),
  fDebugTree(NULL),
  fNparents(-1)
{
  fTRDpid = new AliHFEpidTRD("QAtrdPID");
  fTPCpid = new AliHFEpidTPC("QAtpcPID");
}

AliHFEdebugTreeTask::~AliHFEdebugTreeTask(){
    if(fDebugTree) delete fDebugTree;
    if(fTRDpid) delete fTRDpid;
    if(fTPCpid) delete fTPCpid;
}

void AliHFEdebugTreeTask::UserCreateOutputObjects(){
  //
  // Create debug tree, signal cuts and track cuts
  //
  fDebugTree = new TTreeSRedirector(fFilename.Data());

  fSignalCuts = new AliHFEsignalCuts("HFEsignalCuts", "HFE MC Signal definition");

  fTrackCuts = new AliHFEcuts("fTrackCuts", "Basic HFE track cuts");
  fTrackCuts->CreateStandardCuts();
  // Track cuts
  fTrackCuts->SetMinNClustersTPC(fNclustersTPC);
  fTrackCuts->SetMinRatioTPCclusters(0);
  fTrackCuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable); 
  fTrackCuts->SetMinNClustersTPCPID(fNclustersTPCPID);
  fTrackCuts->SetMinNClustersITS(fNclustersITS);
  // Event cuts
  fTrackCuts->SetUseMixedVertex(kTRUE);
  fTrackCuts->SetVertexRange(10.);
  fTrackCuts->Initialize();

  fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");

}

void AliHFEdebugTreeTask::UserExec(Option_t *){
  //
  // User Exec: Fill debug Tree
  // 

  // Get PID response
  AliPIDResponse *pid = NULL;
  AliInputEventHandler *handler = dynamic_cast<AliInputEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(handler){
    pid = handler->GetPIDResponse();
  } else {
    AliError("No Handler");
  }
  if(!pid){
    AliError("No PID response");
    return;
  }
  if(!fInputEvent) {
    AliError("No Input event");
    return;
  }

  AliESDtrack copyTrack;
  
  fTrackCuts->SetRecEvent(fInputEvent);

  // Cut event
  if(fInputEvent->IsPileupFromSPD(3, 0.8, 3., 2., 5)){
  AliDebug(1, "Event flagged as pileup\n");
    return;
  }
  if(!fTrackCuts->CheckEventCuts("fCutsEvRec", fInputEvent)){
    AliDebug(1, "Event rejected by the event cuts\n");
    return;
  }

  // Get run number
  Int_t run = fInputEvent->GetRunNumber();

  // Derive trigger 
  UInt_t trigger = fInputHandler->IsEventSelected();
  Bool_t isMBTrigger = trigger & AliVEvent::kMB;
  Bool_t isCentralTrigger = trigger & AliVEvent::kCentral;
  Bool_t isSemicentralTrigger = trigger & AliVEvent::kSemiCentral;
  Bool_t isEMCALTrigger = trigger & AliVEvent::kEMCEJE;

  // Check if MC information is available
  Bool_t mcthere = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()) != NULL;
  if(mcthere){ 
    fTrackCuts->SetMCEvent(fMCEvent);
    fSignalCuts->SetMCEvent(fMCEvent);
  }

  // MC get stack
  AliStack* stack = 0x0;
  if(mcthere){
      stack = fMCEvent->Stack();
      if(!stack) AliError("No Stack");
  }

  // Get Primary Vertex
  const AliVVertex *vertex = fInputEvent->GetPrimaryVertex();
  Double_t vtx[3];
  vertex->GetXYZ(vtx);
  Double_t ncontrib = fInputEvent->GetPrimaryVertex()->GetNContributors();

  // Get centrality
  Float_t centrality = -1.;
  AliESDEvent *event = (dynamic_cast<AliESDEvent *>(fInputEvent));
  if(!event) return;
  TString beamtype = event->GetBeamType();
  //printf("Beam type %s\n",(const char*)beamtype);
  if(!beamtype.CompareTo("Pb-Pb") || !beamtype.CompareTo("A-A")){
    // Heavy ion run
    AliDebug(1, "Heavy-Ion event\n");
    AliCentrality *hicent = fInputEvent->GetCentrality();
    centrality = hicent->GetCentralityPercentile("V0M");
  }

  if(!fExtraCuts){
    fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");
  }
  fExtraCuts->SetRecEventInfo(event);

  // Store event selection variables
  (*fDebugTree) << "EventDebug"
                << "Centrality="              << centrality
                << "VertexZ="                 << vtx[2]
                << "NumberOfContributors="    << ncontrib
                << "run="                     << run
                << "\n";

  // Common variables
  Double_t charge, eta, phi, momentum, transversemomentum;
  Int_t source;
      
  // Monte-Carlo loop
  if(mcthere){
    AliMCParticle * mcpart;
    for(Int_t itrk = 0; itrk < fMCEvent->GetNumberOfTracks(); itrk++){
      mcpart = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(itrk));
      if(!mcpart) continue;
      if(!fTrackCuts->CheckParticleCuts(static_cast<UInt_t>(AliHFEcuts::kStepMCGenerated), mcpart)) continue;
      // Kinematics
      charge = mcpart->Charge() > 0. ? 1. : -1.;
      momentum = mcpart->P() * charge;
      transversemomentum = mcpart->Pt() * charge;
      eta = mcpart->Eta();
      phi = mcpart->Phi();
      Int_t pdg = mcpart->Particle()->GetPdgCode();
      
      // Get Production Vertex in radial direction
      Double_t vx = mcpart->Particle()->Vx(), 
              vy = mcpart->Particle()->Vy(); 
      Double_t productionVertex = TMath::Sqrt(vx*vx+vy*vy);
      
      // Get Mother PDG code of the particle
      Int_t motherPdg = 0;
      Int_t motherlabel = mcpart->Particle()->GetFirstMother();
      if(motherlabel >= 0 && motherlabel < fMCEvent->GetNumberOfTracks()){
        AliMCParticle *mother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(motherlabel));
        if(mother) motherPdg = mother->Particle()->GetPdgCode();
      }
      
      // derive source
      source = 5;
      if(fSignalCuts->IsCharmElectron(mcpart)) source = 0;
      else if(fSignalCuts->IsBeautyElectron(mcpart)) source = 1;
      else if(fSignalCuts->IsGammaElectron(mcpart)) source = 2;
      else if(fSignalCuts->IsNonHFElectron(mcpart)) source = 3;
      else if(TMath::Abs(pdg) == 11) source = 4;
      else source = 5;
      
      // Momemntum at the inner wall of the TPC
      Double_t pTPC = 0., ptTPC = 0.;
      AliTrackReference *ref = FindTrackReference(mcpart, 80, 270, AliTrackReference::kTPC);
      if(ref){
        pTPC = ref->P();
        ptTPC = ref->Pt();
      }



    (*fDebugTree) << "MCDebug"
                    << "centrality="          << centrality
                    << "MBtrigger="           << isMBTrigger 
                    << "CentralTrigger="      << isCentralTrigger
                    << "SemicentralTrigger="  << isSemicentralTrigger
                    << "EMCALtrigger="        << isEMCALTrigger
                    << "run="                 << run
                    << "p="                   << momentum
                    << "pt="                  << transversemomentum
                    << "eta="                 << eta
                    << "phi="                 << phi
                    << "pdg="                 << pdg
                    << "ProductionVertex="    << productionVertex
                    << "motherPdg="           << motherPdg
                    << "source="              << source
                    << "\n";
    }
  }
  
  AliESDtrack *track;
  Double_t mcp, mcpt, mcptTPC, mcpTPC;  // MC Variables added to the debug tree
  for(Int_t itrack = 0; itrack < fInputEvent->GetNumberOfTracks(); itrack++){
    // fill the tree
    track = dynamic_cast<AliESDtrack *>(fInputEvent->GetTrack(itrack));
    if(!track) continue;
    // Cut track (Only basic track cuts)
    if(!fTrackCuts->CheckParticleCuts(AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepRecKineITSTPC, track)) continue;
    // Debug streaming of PID-related quantities
    copyTrack.~AliESDtrack();
    new(&copyTrack) AliESDtrack(*track);
    if(fTPCpid->HasEtaCorrection()) fTPCpid->ApplyEtaCorrection(&copyTrack, AliHFEpidObject::kESDanalysis); // Apply Eta Correction on copy track
    Double_t nSigmaTOF = pid->NumberOfSigmasTOF(track, AliPID::kElectron);
    Double_t nSigmaTPC = pid->NumberOfSigmasTPC(&copyTrack, AliPID::kElectron);
    //if(TMath::Abs(nSigmaTOF) > 5) continue;
    // we are not interested in tracks which are more than 5 sigma away from the electron hypothesis in either TOF or TPC
    Double_t tPCdEdx = copyTrack.GetTPCsignal();
    // Signal, source and MCPID
    Bool_t signal = kTRUE;
    source = 5;
    mcp = mcpt = mcpTPC = mcptTPC = 0.;



    Double_t bgcategory = 0.;
    Int_t mArr = -1;
    Int_t mesonID = -999;
    Double_t xr[3]={-999,-999,-999};
    Double_t eR=-999;
    Double_t eZ=-999;
    Double_t unique=-999;
    Double_t mesonunique=-999;
    Double_t mesonR=-999;
    Double_t mesonZ=-999;
    Double_t mesonMomPdg=-999;
    Double_t mesonMomPt=-999;
    Double_t mesonGMomPdg=-999;
    Double_t mesonGGMomPdg=-999;
    Double_t mesonPt = -999;
    Double_t mceta = -999;
    Double_t mcphi = -999;
    Int_t mcpdg;

    if(mcthere){
      // Signal
      AliMCParticle *mctrack;
      if((mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(track->GetLabel()))))){
        if(!fTrackCuts->CheckParticleCuts(AliHFEcuts::kStepMCGenerated, mctrack)) signal = kFALSE; 
      }
  
      // Source
      if(fSignalCuts->IsCharmElectron(track)) source = 0;
      else if(fSignalCuts->IsBeautyElectron(track)) source = 1;
      else if(fSignalCuts->IsGammaElectron(track)) source = 2;
      else if(fSignalCuts->IsNonHFElectron(track)) source = 3;
      else if(mctrack && (TMath::Abs(mctrack->Particle()->GetPdgCode()) == 11)) source = 4;
      else source = 5;

      if(!mctrack) continue;

      // Kinematics
      mcpt  = mctrack->Pt();
      mcp   = mctrack->P();
      mceta = mctrack->Eta();
      mcphi = mctrack->Phi();
      mcpdg = mctrack->Particle()->GetPdgCode();

      AliTrackReference *ref = FindTrackReference(mctrack, 80, 270, AliTrackReference::kTPC);
      if(ref){
        mcpTPC = ref->P();
        mcptTPC = ref->Pt();
      }

    
      TParticle *mctrack1 = mctrack->Particle();
      mesonID=GetElecSourceMC(mctrack1);
      if(mesonID==AliHFEmcQA::kGammaPi0 || mesonID==AliHFEmcQA::kPi0) mArr=0;                //pion
      else if(mesonID==AliHFEmcQA::kGammaEta || mesonID==AliHFEmcQA::kEta) mArr=1;           //eta
      else if(mesonID==AliHFEmcQA::kGammaOmega || mesonID==AliHFEmcQA::kOmega) mArr=2;       //omega
      else if(mesonID==AliHFEmcQA::kGammaPhi || mesonID==AliHFEmcQA::kPhi) mArr=3;           //phi
      else if(mesonID==AliHFEmcQA::kGammaEtaPrime || mesonID==AliHFEmcQA::kEtaPrime) mArr=4; //etaprime
      else if(mesonID==AliHFEmcQA::kGammaRho0 || mesonID==AliHFEmcQA::kRho0) mArr=5;         //rho
    
      mctrack->XvYvZv(xr);
    
      eR= TMath::Sqrt(xr[0]*xr[0]+xr[1]*xr[1]);
      eZ = xr[2];
      TParticle *mctrackt = mctrack->Particle();
      unique=mctrackt->GetUniqueID();
    
      AliMCParticle *mctrackmother = NULL;
      AliMCParticle *mctrackmother2 = NULL;

      if(!(mArr<0)){
    if(mesonID>=AliHFEmcQA::kGammaPi0) {  // conversion electron, be careful with the enum odering
        Int_t glabel=TMath::Abs(mctrack->GetMother()); // gamma label
        if((mctrackmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(glabel)))){
      glabel=TMath::Abs(mctrackmother->GetMother()); // gamma's mother's label
      if((mctrackmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(glabel)))){
          mesonPt = mctrackmother->Pt(); //meson pt
          bgcategory = 1.;
          mctrackmother->XvYvZv(xr);
          mesonR = TMath::Sqrt(xr[0]*xr[0]+xr[1]*xr[1]);
          mesonZ = xr[2];

          mctrackt = mctrackmother->Particle();
          if(mctrackt){
            mesonunique = mctrackt->GetUniqueID();
          }

        Int_t glabel2=TMath::Abs(mctrackmother->GetMother()); // gamma's mother's mother
        if((mctrackmother2 = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(glabel2)))){
            mesonMomPdg=mctrackmother2->PdgCode();
            mesonMomPt=mctrackmother2->Pt();
        }

          if(glabel>fMCEvent->GetNumberOfPrimaries()) {
        bgcategory = 2.;
        glabel=TMath::Abs(mctrackmother->GetMother()); // gamma's mother's mother
        if((mctrackmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(glabel)))){
            mesonMomPdg=mctrackmother->PdgCode();
            mesonMomPt=mctrackmother->Pt();
            if(TMath::Abs(mctrackmother->PdgCode())==310){
          bgcategory = 3.;
          glabel=TMath::Abs(mctrackmother->GetMother()); // gamma's mother's mother's mother
          if((mctrackmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(glabel)))){
              mesonGMomPdg=mctrackmother->PdgCode();
              glabel=TMath::Abs(mctrackmother->GetMother()); // gamma's mother's mother
              if((mctrackmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(glabel)))){
            mesonGGMomPdg=mctrackmother->PdgCode();
              }
          }
            }
        }
          }
      }
        }
    }
    else{ // nonHFE except for the conversion electron
        Int_t glabel=TMath::Abs(mctrack->GetMother());
        if((mctrackmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(glabel)))){
      mesonPt = mctrackmother->Pt(); //meson pt
      bgcategory = -1.;
      mctrackmother->XvYvZv(xr);
      mesonR = TMath::Sqrt(xr[0]*xr[0]+xr[1]*xr[1]);
      mesonZ = xr[2];

      mctrackt = mctrackmother->Particle();
      if(mctrackt){
          mesonunique = mctrackt->GetUniqueID();
      }
      if(glabel>fMCEvent->GetNumberOfPrimaries()) {
          bgcategory = -2.;
          glabel=TMath::Abs(mctrackmother->GetMother()); // gamma's mother's mother
          if((mctrackmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(glabel)))){
        mesonMomPdg=mctrackmother->PdgCode();
        mesonMomPt=mctrackmother->Pt();
        if(TMath::Abs(mctrackmother->PdgCode())==310){
            bgcategory = -3.;
            glabel=TMath::Abs(mctrackmother->GetMother()); // gamma's mother's mother's mother
            if((mctrackmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(glabel)))){
                                  mesonGMomPdg=mctrackmother->PdgCode();
          glabel=TMath::Abs(mctrackmother->GetMother()); // gamma's mother's mother
          if((mctrackmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(glabel)))){
                                    mesonGGMomPdg=mctrackmother->PdgCode();
          }
            }
        }
          }
      }
        }
    }
      }



    }
    // Get V0 tag (if available)
    Int_t v0pid = -1;
    if(track->TestBit(BIT(14))) v0pid = AliPID::kElectron;
    else if(track->TestBit(BIT(15))) v0pid = AliPID::kPion;
    else if(track->TestBit(BIT(16))) v0pid = AliPID::kProton;
    // Kinematics
    charge = track->Charge() > 0 ? 1. : -1.;
    eta = track->Eta();
    phi = track->Phi();
    momentum = track->P() * charge;
    transversemomentum = track->Pt() * charge;
    Double_t momentumTPC = track->GetTPCInnerParam() ? track->GetTPCInnerParam()->P() : 0.;
    Double_t transversemomentumTPC = track->GetTPCInnerParam() ? track->GetTPCInnerParam()->Pt() : 0.;
    // ITS number of clusters
    UChar_t nclustersITS = track->GetITSclusters(NULL);
    Double_t chi2matching =  track->GetChi2TPCConstrainedVsGlobal(dynamic_cast<const AliESDVertex *>(vertex));
    Double_t chi2PerClusterITS = 0.0;
    if (nclustersITS != 0) chi2PerClusterITS = track->GetITSchi2() / Float_t(nclustersITS);
    // TPC number of clusters (different definitions)
    UChar_t nclustersTPCfit = track->GetTPCNcls();
    UChar_t nclustersTPCall = 0;
    const TBits &clusterTPC = track->GetTPCClusterMap();
    nclustersTPCall = clusterTPC.CountBits();
    UChar_t nclustersTPCPID = track->GetTPCsignalN();
    UChar_t nfindableTPC =  track->GetTPCNclsF();
    Double_t clusterRatioTPCfit = 0.0;
    if((static_cast<Double_t>(nfindableTPC))>0.0) clusterRatioTPCfit = static_cast<Double_t>(nclustersTPCfit)/static_cast<Double_t>(nfindableTPC);
    Double_t clusterRatioTPCall = 0.0;
    if((static_cast<Double_t>(nfindableTPC))>0.0) clusterRatioTPCall = static_cast<Double_t>(nclustersTPCall)/static_cast<Double_t>(nfindableTPC);
    UChar_t nclustersTPCshared = 0;
    Float_t ncrossedRowsTPC = track->GetTPCCrossedRows();
    const TBits &sharedTPC = track->GetTPCSharedMap();
    for(Int_t ibit = 0; ibit < 160; ibit++) if(sharedTPC.TestBitNumber(ibit)) nclustersTPCshared++;
    // TRD clusters and tracklets
    UChar_t nclustersTRD = track->GetTRDncls();
    UChar_t ntrackletsTRDPID = track->GetTRDntrackletsPID();
    // ITS and TRD acceptance maps
    UChar_t hasClusterITS[6], statusITS[6], hasTrackletTRD[6];
    UChar_t itsPixel = track->GetITSClusterMap();
    for(Int_t icl = 0; icl < 6; icl++){ 
      hasClusterITS[icl] = TESTBIT(itsPixel, icl) ? 1 : 0;
      if(CheckITSstatus(track, icl)) statusITS[icl] = 1;
      else statusITS[icl] = 0;
    }
    Double_t trddEdxSum[6];
    for(Int_t a=0;a<6;a++) { trddEdxSum[a]= 0.;}
    for(Int_t itl = 0; itl < 6; itl++){
      Int_t nSliceNonZero = 0;
      trddEdxSum[itl] = track->GetTRDslice(itl, 0); // in new reconstruction slice 0 contains the total charge
      for(Int_t islice = 0; islice < 8; islice++){
        if(track->GetTRDslice(itl, islice) > 0.001) nSliceNonZero++;
      }
      hasTrackletTRD[itl] = nSliceNonZero ? 1 : 0;
    }
    // TRD PID
    Double_t pidprobs[5];
    track->GetTRDpid(pidprobs);
    Double_t likeEleTRD = pidprobs[0];
    Double_t likeEleTRDn = likeEleTRD/(likeEleTRD + pidprobs[2]);
    Double_t trdtruncmean1 = fTRDpid->GetTRDSignalV1(track, 0.6);
    Double_t trdtruncmean2 = fTRDpid->GetTRDSignalV2(track, 0.6);

    // DCA
    Float_t b[2] = {0.,0.};
    Float_t bCov[3] = {0.,0.,0.};
    track->GetImpactParameters(b,bCov);
    Double_t dca = TMath::Sqrt(b[0]*b[0]+b[1]*b[1]); // impact parameter space
    Double_t dcaSR=0, dcaSZ=0;
    if(bCov[0]>0) dcaSR = b[0]/TMath::Sqrt(bCov[0]); // normalised impact parameter xy
    if(bCov[2]>0) dcaSZ = b[1]/TMath::Sqrt(bCov[2]); // normalised impact parameter z
    Double_t dcaS = AliESDtrackCuts::GetSigmaToVertex(track); // n_sigma

    // HFE DCA
    Double_t hfeb[2] = {-99.,-99.};
    Double_t hfebCov[3] = {-999.,-999.,-999.};
    fExtraCuts->GetHFEImpactParameters(track, hfeb, hfebCov);

    Double_t tofdx= -999.0;
    Double_t tofdz= -999.0;
    tofdx=track->GetTOFsignalDx();
    tofdz=track->GetTOFsignalDz();

    // TOF track status
    UInt_t status = 0;
    status = track->GetStatus();
    Bool_t hasTOFout  = status&AliESDtrack::kTOFout; 
    Bool_t hasTOFtime = status&AliESDtrack::kTIME;
    Bool_t hasTOFpid  = status&AliESDtrack::kTOFpid;
    Bool_t hasgoodTOF     = kFALSE;
    if (hasTOFout && hasTOFtime && hasTOFpid) hasgoodTOF = kTRUE;

    // TRD track status
    Bool_t hasTRDin = status&AliESDtrack::kTRDin;


    // TOF mismatch (particle spectra group)
    Int_t mismatchlevel=0; // default value; in data always 0
    if(mcthere){
	Int_t tofLabel[3];
	track->GetTOFLabel(tofLabel);
        if(TMath::Abs(track->GetLabel()) != TMath::Abs(tofLabel[0]) || tofLabel[1] > 0) mismatchlevel=1;
	TParticle *matchedTrack = stack->Particle(TMath::Abs(tofLabel[0]));
	if(TMath::Abs(matchedTrack->GetFirstMother()) == TMath::Abs(track->GetLabel()))
	{
	    if(mismatchlevel==1) mismatchlevel=3;
            else mismatchlevel=2;
	}
    }


    // Fill Tree
    (*fDebugTree) << "PIDdebug"
                  << "centrality="          << centrality
                  << "MBtrigger="           << isMBTrigger 
                  << "CentralTrigger="      << isCentralTrigger
                  << "SemicentralTrigger="  << isSemicentralTrigger
                  << "EMCALtrigger="        << isEMCALTrigger
                  << "signal="              << signal
                  << "source="              << source
                  << "v0pid="               << v0pid
                  << "run="                 << run
                  << "p="                   << momentum
                  << "ptpc="                << momentumTPC
                  << "pt="                  << transversemomentum
                  << "pttpc="               << transversemomentumTPC
                  << "mcp="                 << mcp
                  << "mcpt="                << mcpt
                  << "mcpTPC="              << mcpTPC
                  << "mcptTPC="             << mcptTPC
                  << "mceta="               << mceta
                  << "mcphi="               << mcphi
                  << "mcpdg="               << mcpdg
                  << "eta="                 << eta
                  << "phi="                 << phi
                  << "ntracklets="          << ntrackletsTRDPID
                  << "nclustersTPC="        << nclustersTPCfit
		  << "nclustersTPCall="     << nclustersTPCall
                  << "nclustersTPCPID="     << nclustersTPCPID
                  << "nclustersTPCshared="  << nclustersTPCshared
                  << "ncrossedRowsTPC="     << ncrossedRowsTPC
                  << "clusterRatioTPC="     << clusterRatioTPCfit
                  << "clusterRatioTPCall="  << clusterRatioTPCall
                  << "nclustersITS="        << nclustersITS
                  << "nclusters="           << nclustersTRD
                  << "chi2matching="        << chi2matching
                  << "chi2PerClusterITS="   << chi2PerClusterITS
                  << "its0="                << hasClusterITS[0]
                  << "its1="                << hasClusterITS[1]
                  << "its2="                << hasClusterITS[2]
                  << "its3="                << hasClusterITS[3]
                  << "its4="                << hasClusterITS[4]
                  << "its5="                << hasClusterITS[5]
                  << "statusITS0="          << statusITS[0]
                  << "statusITS1="          << statusITS[1]
                  << "statusITS2="          << statusITS[2]
                  << "statusITS3="          << statusITS[3]
                  << "statusITS4="          << statusITS[4]
                  << "statusITS5="          << statusITS[5]
                  << "trd0="                << hasTrackletTRD[0]
                  << "trd1="                << hasTrackletTRD[1]
                  << "trd2="                << hasTrackletTRD[2]
                  << "trd3="                << hasTrackletTRD[3]
                  << "trd4="                << hasTrackletTRD[4]
                  << "trd5="                << hasTrackletTRD[5]
                  << "TRDdEdxl0="           << trddEdxSum[0]
                  << "TRDdEdxl1="           << trddEdxSum[1]
                  << "TRDdEdxl2="           << trddEdxSum[2]
                  << "TRDdEdxl3="           << trddEdxSum[3]
                  << "TRDdEdxl4="           << trddEdxSum[4]
                  << "TRDdEdxl5="           << trddEdxSum[5]
                  << "TOFsigmaEl="          << nSigmaTOF
                  << "TPCsigmaEl="          << nSigmaTPC
                  << "TPCdEdx="             << tPCdEdx
                  << "TRDlikeEl="           << likeEleTRD
                  << "TRDlikeEln="          << likeEleTRDn
                  << "trdtruncmean1="       << trdtruncmean1
                  << "trdtruncmean2="       << trdtruncmean2
                  << "dcaR="                << b[0]
                  << "dcaZ="                << b[1]
                  << "dca="                 << dca
                  << "dcaSR="               << dcaSR
                  << "dcaSZ="               << dcaSZ
                  << "dcaS="                << dcaS
                  << "hfedcaR="             << hfeb[0]
                  << "hfedcaZ="             << hfeb[1]
                  << "hfedcacovR="          << hfebCov[0]
                  << "hfedcacovZ="          << hfebCov[2]
                  << "vx="                  << vtx[0]
                  << "vy="                  << vtx[1]
                  << "vz="                  << vtx[2]
                  << "tofdx="                << tofdx
	          << "tofdz="                << tofdz
                  << "statusTOFtracking="   << hasgoodTOF
                  << "TOFmismatchlevel="    << mismatchlevel
                  << "statusTRDtracking="   << hasTRDin
                  << "ncontrib="            << ncontrib
                  << "mesonID="             << mesonID
                  << "eR="                  << eR
                  << "mesonR="              << mesonR
                  << "eZ="                  << eZ
                  << "mesonZ="              << mesonZ
                  << "unique="              << unique
                  << "mesonunique="         << mesonunique
                  << "bgcategory="          << bgcategory
                  << "mesonpt="             << mesonPt
                  << "mesonMomPdg="         << mesonMomPdg
                  << "mesonGMomPdg="         << mesonGMomPdg
                  << "mesonGGMomPdg="         << mesonGGMomPdg
                  << "mesonMomPt="         << mesonMomPt
                  << "\n";
  }
}


void AliHFEdebugTreeTask::SetFileName(const char *filename){ fFilename = filename; }
//___________________________________________________________
AliTrackReference *AliHFEdebugTreeTask::FindTrackReference(AliMCParticle *track, Float_t minRadius, Float_t maxRadius, Int_t detectorID)
{
  //
  // Find the track reference
  //
  AliTrackReference *ref = NULL, *reftmp;
  Float_t radius;
  for(Int_t iref = 0; iref < track->GetNumberOfTrackReferences(); iref++){
    reftmp = track->GetTrackReference(iref);
    if(reftmp->DetectorId() != detectorID) continue;
    radius = reftmp->R();
    if(radius >= minRadius && radius < maxRadius){
      ref = reftmp;
      break;
    } 
    if(radius > maxRadius) break;
  }
  return ref;
}

//__________________________________________
Int_t AliHFEdebugTreeTask::GetElecSourceMC(TParticle * const mcpart)
{
  // decay particle's origin 

  if(!mcpart){
    AliDebug(1, "no mcparticle, return\n");
    return -1;
  }

  if ( abs(mcpart->GetPdgCode()) != AliHFEmcQA::kElectronPDG ) return AliHFEmcQA::kMisID;

  Int_t origin = -1;
  Bool_t isFinalOpenCharm = kFALSE;

  Int_t iLabel = mcpart->GetFirstMother();
  if (iLabel<0){
    AliDebug(1, "Stack label is negative, return\n");
    return -1;
  }

  AliMCParticle *mctrack = NULL;
  Int_t tmpMomLabel=0;
  if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(iLabel))))) return -1; 
  TParticle *partMother = mctrack->Particle();
  TParticle *partMotherCopy = mctrack->Particle();
  Int_t maPdgcode = partMother->GetPdgCode();

  // if the mother is charmed hadron  
  if ( (int(abs(maPdgcode)/100.)%10) == AliHFEmcQA::kCharm || (int(abs(maPdgcode)/1000.)%10) == AliHFEmcQA::kCharm ) {

    for (Int_t i=0; i<fNparents; i++){
        if (abs(maPdgcode)==fParentSelect[0][i]){
          isFinalOpenCharm = kTRUE;
        }
    }
    if (!isFinalOpenCharm) return -1;

    // iterate until you find B hadron as a mother or become top ancester 
    for (Int_t i=1; i<fgkMaxIter; i++){

        Int_t jLabel = partMother->GetFirstMother();
        if (jLabel == -1){
          origin = AliHFEmcQA::kDirectCharm;
          return origin;
        }
        if (jLabel < 0){ // safety protection
          AliDebug(1, "Stack label is negative, return\n");
          return -1;
        }

        // if there is an ancester
        if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(jLabel))))) return -1; 
        TParticle* grandMa = mctrack->Particle();
        Int_t grandMaPDG = grandMa->GetPdgCode();

        for (Int_t j=0; j<fNparents; j++){
          if (abs(grandMaPDG)==fParentSelect[1][j]){
            origin = AliHFEmcQA::kBeautyCharm;
            return origin;
          }
        }

        partMother = grandMa;
    } // end of iteration 
  } // end of if
  else if ( (int(abs(maPdgcode)/100.)%10) == AliHFEmcQA::kBeauty || (int(abs(maPdgcode)/1000.)%10) == AliHFEmcQA::kBeauty ) {
    for (Int_t i=0; i<fNparents; i++){
        if (abs(maPdgcode)==fParentSelect[1][i]){
          origin = AliHFEmcQA::kDirectBeauty;
          return origin;
        }
    }
  } // end of if
  else if ( abs(maPdgcode) == 22 ) { //conversion

    tmpMomLabel = partMotherCopy->GetFirstMother();
    if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(tmpMomLabel))))) return -1;
    partMother = mctrack->Particle();
    maPdgcode = partMother->GetPdgCode();
    if ( abs(maPdgcode) == 111 ) {
      origin = AliHFEmcQA::kGammaPi0;
      return origin;
    } 
    else if ( abs(maPdgcode) == 221 ) {
      origin = AliHFEmcQA::kGammaEta;
      return origin;
    } 
    else if ( abs(maPdgcode) == 223 ) {
      origin = AliHFEmcQA::kGammaOmega;
      return origin;
    } 
    else if ( abs(maPdgcode) == 333 ) {
      origin = AliHFEmcQA::kGammaPhi;
      return origin;
    }
    else if ( abs(maPdgcode) == 331 ) {
      origin = AliHFEmcQA::kGammaEtaPrime;
      return origin; 
    }
    else if ( abs(maPdgcode) == 113 ) {
      origin = AliHFEmcQA::kGammaRho0;
      return origin;
    }
    else origin = AliHFEmcQA::kElse;
    //origin = kGamma; // finer category above
    return origin;

  } // end of if
  else if ( abs(maPdgcode) == 111 ) {
    origin = AliHFEmcQA::kPi0;
    return origin;
  } // end of if
  else if ( abs(maPdgcode) == 221 ) {
    origin = AliHFEmcQA::kEta;
    return origin;
  } // end of if
  else if ( abs(maPdgcode) == 223 ) {
    origin = AliHFEmcQA::kOmega;
    return origin;
  } // end of if
  else if ( abs(maPdgcode) == 333 ) {
    origin = AliHFEmcQA::kPhi;
    return origin;
  } // end of if
  else if ( abs(maPdgcode) == 331 ) {
    origin = AliHFEmcQA::kEtaPrime;
    return origin;
  } // end of if
  else if ( abs(maPdgcode) == 113 ) {
    origin = AliHFEmcQA::kRho0;
    return origin;
  } // end of if
  else{ 
    origin = AliHFEmcQA::kElse;
  }
  return origin;
}

//______________________________________________________
Bool_t AliHFEdebugTreeTask::CheckITSstatus( const AliESDtrack * const esdtrack, Int_t layer) const {
  //
  // Check whether ITS area is dead
  //
  Int_t itsStatus = 0;
  Int_t det;
  Float_t xloc, zloc;
  esdtrack->GetITSModuleIndexInfo(layer, det, itsStatus, xloc, zloc);
  Bool_t status;
  switch(itsStatus){
    case 2: status = kFALSE; break;
    case 3: status = kFALSE; break;
    case 7: status = kFALSE; break;
    default: status = kTRUE;
  }
  return status;
}

