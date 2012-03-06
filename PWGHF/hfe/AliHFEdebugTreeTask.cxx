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
#include "AliPIDResponse.h"
#include "AliVEvent.h"
#include "AliHFEpidTRD.h"
#include "TTreeStream.h"

#include "AliHFEdebugTreeTask.h"

ClassImp(AliHFEdebugTreeTask)

AliHFEdebugTreeTask::AliHFEdebugTreeTask():
  AliAnalysisTaskSE(),
  fTrackCuts(NULL),
  fSignalCuts(NULL),
  fTRDpid(NULL),
  fNclustersTPC(70),
  fNclustersTPCPID(0),
  fNclustersITS(2),
  fFilename("HFEtree.root"),
  fDebugTree(NULL)
{

}

AliHFEdebugTreeTask::AliHFEdebugTreeTask(const char *name):
  AliAnalysisTaskSE(name),
  fTrackCuts(NULL),
  fSignalCuts(NULL),
  fTRDpid(NULL),
  fNclustersTPC(70),
  fNclustersTPCPID(0),
  fNclustersITS(2),
  fFilename("HFEtree.root"),
  fDebugTree(NULL)
{

}

AliHFEdebugTreeTask::~AliHFEdebugTreeTask(){
    if(fDebugTree) delete fDebugTree;
    if(fTRDpid) delete fTRDpid;
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

  fTRDpid = new AliHFEpidTRD("QAtrdPID");

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

  // Get Primary Vertex
  const AliVVertex *vertex = fInputEvent->GetPrimaryVertex();

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
  for(Int_t itrack = 0; itrack < fInputEvent->GetNumberOfTracks(); itrack++){
    // fill the tree
    track = dynamic_cast<AliESDtrack *>(fInputEvent->GetTrack(itrack));
    if(!track) continue;
    // Cut track (Only basic track cuts)
    if(!fTrackCuts->CheckParticleCuts(AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepRecKineITSTPC, track)) continue;
    // Debug streaming of PID-related quantities
    Double_t nSigmaTOF = pid->NumberOfSigmasTOF(track, AliPID::kElectron);
    Double_t nSigmaTPC = pid->NumberOfSigmasTPC(track, AliPID::kElectron);
    if(TMath::Abs(nSigmaTOF) > 5) continue;
    // we are not interested in tracks which are more than 5 sigma away from the electron hypothesis in either TOF or TPC
    Double_t tPCdEdx = track->GetTPCsignal();
    // Signal, source and MCPID
    Bool_t signal = kTRUE;
    source = 5;
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
    // ITS number of clusters
    UChar_t nclustersITS = track->GetITSclusters(NULL);
    Double_t chi2matching =  track->GetChi2TPCConstrainedVsGlobal(dynamic_cast<const AliESDVertex *>(vertex));
    Double_t chi2PerClusterITS = 0.0;
    if (nclustersITS != 0) chi2PerClusterITS = track->GetITSchi2() / Float_t(nclustersITS);
    // TPC number of clusters (different definitions)
    UChar_t nclustersTPC = track->GetTPCncls();
    UChar_t nclustersTPCPID = track->GetTPCsignalN();
    UChar_t nfindableTPC =  track->GetTPCNclsF();
    Double_t clusterRatioTPC = 0.0;
    if((static_cast<Double_t>(nfindableTPC))>0.0) clusterRatioTPC = static_cast<Double_t>(nclustersTPC)/static_cast<Double_t>(nfindableTPC);
    UChar_t nclustersTPCshared = 0;
    Float_t ncrossedRowsTPC = track->GetTPCCrossedRows();
    const TBits &sharedTPC = track->GetTPCSharedMap();
    for(Int_t ibit = 0; ibit < 160; ibit++) if(sharedTPC.TestBitNumber(ibit)) nclustersTPCshared++;
    // TRD clusters and tracklets
    UChar_t nclustersTRD = track->GetTRDncls();
    UChar_t ntrackletsTRDPID = track->GetTRDntrackletsPID();
    // ITS and TRD acceptance maps
    UChar_t hasClusterITS[6], hasTrackletTRD[6];
    UChar_t itsPixel = track->GetITSClusterMap();
    for(Int_t icl = 0; icl < 6; icl++) hasClusterITS[icl] = TESTBIT(itsPixel, icl) ? 1 : 0;
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
    // Vertex
    Double_t vtx[3];
    vertex->GetXYZ(vtx);
    Double_t ncontrib = fInputEvent->GetPrimaryVertex()->GetNContributors();
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
                  << "pt="                  << transversemomentum
                  << "eta="                 << eta
                  << "phi="                 << phi
                  << "ntracklets="          << ntrackletsTRDPID
                  << "nclustersTPC="        << nclustersTPC
                  << "nclustersTPCPID="     << nclustersTPCPID
                  << "nclustersTPCshared="  << nclustersTPCshared
                  << "ncrossedRowsTPC="     << ncrossedRowsTPC
                  << "clusterRatioTPC="     << clusterRatioTPC
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
                  << "vx="                  << vtx[0]
                  << "vy="                  << vtx[1]
                  << "vz="                  << vtx[2] 
                  << "ncontrib="            << ncontrib
                  << "\n";
  }
}


void AliHFEdebugTreeTask::SetFileName(const char *filename){ fFilename = filename; }
