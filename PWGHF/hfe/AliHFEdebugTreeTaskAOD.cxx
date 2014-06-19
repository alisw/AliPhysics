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
//   R.Bailhache <R.Bailhache@gsi.de>
//
#include <TBits.h>
#include <TString.h>
#include <TArrayI.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODPid.h"
#include "AliHFEcuts.h"
#include "AliHFEextraCuts.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliPIDResponse.h"
#include "AliVEvent.h"
#include "AliHFEpidTPC.h"
#include "TTreeStream.h"
#include "AliESDtrack.h"
#include "TClonesArray.h"
#include "AliAODMCHeader.h"
#include "AliHFEsignalCuts.h"
#include "AliAODMCParticle.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "AliAODVertex.h"

#include "AliHFEdebugTreeTaskAOD.h"

ClassImp(AliHFEdebugTreeTaskAOD)

AliHFEdebugTreeTaskAOD::AliHFEdebugTreeTaskAOD():
  AliAnalysisTaskSE(),
  fAODMCHeader(NULL),
  fAODArrayMCInfo(NULL),
  fTrackCuts(NULL),
  fExtraCuts(NULL),
  fSignalCuts(NULL),
  fTPCpid(NULL),
  fEventNumber(0),
  fNclustersTPC(70),
  fNclustersTPCPID(0),
  fNclustersITS(2),
  fFilename("HFEtree.root"),
  fDebugTree(NULL)
{

}

AliHFEdebugTreeTaskAOD::AliHFEdebugTreeTaskAOD(const char *name):
  AliAnalysisTaskSE(name),
  fAODMCHeader(NULL),
  fAODArrayMCInfo(NULL),
  fTrackCuts(NULL),
  fExtraCuts(NULL),
  fSignalCuts(NULL),
  fTPCpid(NULL),
  fEventNumber(0),
  fNclustersTPC(70),
  fNclustersTPCPID(0),
  fNclustersITS(2),
  fFilename("HFEtree.root"),
  fDebugTree(NULL)
{
  fTPCpid = new AliHFEpidTPC("QAtpcPID");
}

AliHFEdebugTreeTaskAOD::~AliHFEdebugTreeTaskAOD(){

    if(fDebugTree) delete fDebugTree;
    if(fTPCpid) delete fTPCpid;
}

void AliHFEdebugTreeTaskAOD::UserCreateOutputObjects(){
  //
  // Create debug tree, signal cuts and track cuts
  //

  //printf("test\n");
  fDebugTree = new TTreeSRedirector(fFilename.Data());

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

}

void AliHFEdebugTreeTaskAOD::UserExec(Option_t *){
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

  AliAODTrack copyTrack;

  // MC info
  Bool_t mcthere = kTRUE;
  AliAODEvent *aodE = dynamic_cast<AliAODEvent *>(fInputEvent);
  if(!aodE){
 //        printf("testd\n");
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
  fTrackCuts->SetRecEvent(fInputEvent);

  if(!fTrackCuts->CheckEventCuts("fCutsEvRec", fInputEvent)){
    AliDebug(1, "Event rejected by the event cuts\n");
    return;
  }
  if(!fExtraCuts){
    fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");
  }
  fExtraCuts->SetRecEventInfo(fInputEvent);

 
 

  // Get run number
  Int_t run = fInputEvent->GetRunNumber();

  // Derive trigger 
  UInt_t trigger = fInputHandler->IsEventSelected();
  Bool_t isMBTrigger = trigger & AliVEvent::kMB;
  Bool_t isCentralTrigger = trigger & AliVEvent::kCentral;
  Bool_t isSemicentralTrigger = trigger & AliVEvent::kSemiCentral;
  Bool_t isEMCALTrigger = trigger & AliVEvent::kEMCEJE;

  // Get Primary Vertex
  const AliVVertex *vertex = fInputEvent->GetPrimaryVertex();
  Double_t vtx[3];
  vertex->GetXYZ(vtx);
  Double_t ncontrib = fInputEvent->GetPrimaryVertex()->GetNContributors();

  // Get centrality
  Float_t centrality = -1.;
  AliCentrality *hicent = fInputEvent->GetCentrality();
  centrality = hicent->GetCentralityPercentile("V0M");

  // Store event selection variables
  (*fDebugTree) << "EventDebug"
                << "Centrality="              << centrality
                << "VertexZ="                 << vtx[2]
                << "NumberOfContributors="    << ncontrib
                << "run="                 << run
                << "\n";



  // Look for kink mother
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
  AliAODMCParticle *mctrack = NULL;
  // Monte-Carlo info
  Double_t eR,vx,vy,vz;
  Double_t chargemc, etamc, phimc, momentummc, transversemomentummc;
  Int_t source,pdg,signal;
  if(mcthere) {
    for(Int_t itrack = 0; itrack < fAODArrayMCInfo->GetEntriesFast(); itrack++) {
      mctrack = (AliAODMCParticle *)(fAODArrayMCInfo->At(itrack));
      if(!mctrack) continue;
      signal = 0;
      if(fTrackCuts->CheckParticleCuts(static_cast<UInt_t>(AliHFEcuts::kStepMCGenerated), mctrack)) signal = 1;
      // Kinematics
      chargemc = mctrack->Charge() > 0. ? 1. : -1.;
      momentummc = mctrack->P() * chargemc;
      transversemomentummc = mctrack->Pt() * chargemc;
      etamc = mctrack->Eta();
      phimc = mctrack->Phi();
      pdg = mctrack->GetPdgCode();
      
      // Get Production Vertex in radial direction
      vx = mctrack->Xv();
      vy = mctrack->Yv(); 
      vz = mctrack->Zv(); 
      eR = TMath::Sqrt(vx*vx+vy*vy);

      // Get Mother PDG code of the particle
      Int_t motherPdg = 0;
      Int_t motherlabel = TMath::Abs(mctrack->GetMother());
      if(motherlabel >= 0 && motherlabel < fAODArrayMCInfo->GetEntriesFast()){
        AliAODMCParticle *mother = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(motherlabel));
        if(mother) motherPdg = mother->GetPdgCode();
      }
      
      // derive source
      source = 5;
      if(fSignalCuts->IsCharmElectron(mctrack)) source = 0;
      else if(fSignalCuts->IsBeautyElectron(mctrack)) source = 1;
      else if(fSignalCuts->IsGammaElectron(mctrack)) source = 2;
      else if(fSignalCuts->IsNonHFElectron(mctrack)) source = 3;
      else if(TMath::Abs(pdg) == 11) source = 4;
      else source = 5;

      (*fDebugTree) << "MCDebug"
                    << "centrality="          << centrality
                    << "MBtrigger="           << isMBTrigger 
                    << "CentralTrigger="      << isCentralTrigger
                    << "SemicentralTrigger="  << isSemicentralTrigger
                    << "EMCALtrigger="        << isEMCALTrigger
                    << "run="                 << run
                    << "p="                   << momentummc
                    << "pt="                  << transversemomentummc
                    << "eta="                 << etamc
                    << "phi="                 << phimc
                    << "pdg="                 << pdg
		                << "px="                  << vx
		                << "py="                  << vy
		                << "pz="                  << vz
                    << "ProductionVertex="    << eR
                    << "motherPdg="           << motherPdg
                    << "source="              << source
                    << "\n";
    }
  }
  
  // Common variables
  Double_t charge, eta, phi, momentum, momentumTPC, transversemomentum;
  
  //
  // Loop on reconstructed tracks
  //

  TArrayI arraytrack(fInputEvent->GetNumberOfTracks());
  Int_t counterdc=0;
  
  AliAODTrack *track = 0x0;
  for(Int_t itrack = 0; itrack < fInputEvent->GetNumberOfTracks(); itrack++){
    // fill the tree
    track = dynamic_cast<AliAODTrack *>(fInputEvent->GetTrack(itrack));
    if(!track) continue;
    // Cut track (Only basic track cuts)
    // printf("testv\n");
    if(!fTrackCuts->CheckParticleCuts(AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepRecKineITSTPC, track)) continue;
    //
    //printf("testu\n");
    Double_t nSigmaTOF = pid->NumberOfSigmasTOF(track, AliPID::kElectron);
    Double_t nSigmaTPC = pid->NumberOfSigmasTPC(track, AliPID::kElectron);
    Double_t tPCdEdx = track->GetDetPid() ? track->GetDetPid()->GetTPCsignal() : 0.;
    // Eta correction
    copyTrack.~AliAODTrack();
    new(&copyTrack) AliAODTrack(*track);
    if(fTPCpid->HasCentralityCorrection()) fTPCpid->ApplyCentralityCorrection(&copyTrack, static_cast<Double_t>(ncontrib),AliHFEpidObject::kAODanalysis);
    if(fTPCpid->HasEtaCorrection()) fTPCpid->ApplyEtaCorrection(&copyTrack, AliHFEpidObject::kAODanalysis);
    Double_t nSigmaTPCcorr = pid->NumberOfSigmasTPC(&copyTrack, AliPID::kElectron);
    Double_t tPCdEdxcorr = copyTrack.GetDetPid() ? copyTrack.GetDetPid()->GetTPCsignal() : 0.;
   
    // Kinematics
    charge = track->Charge() > 0 ? 1. : -1.;
    eta = track->Eta();
    phi = track->Phi();
    momentum = track->P() * charge;
    transversemomentum = track->Pt() * charge;
    momentumTPC = track->GetDetPid() ? track->GetDetPid()->GetTPCmomentum() : track->P();
  
    // status
    ULong_t status = track->GetStatus();
    Int_t itsrefit=0;
    if((status & AliESDtrack::kITSrefit) == AliESDtrack::kITSrefit) itsrefit = 1;
    Int_t tpcrefit=0;
    if((status & AliESDtrack::kTPCrefit) == AliESDtrack::kTPCrefit) tpcrefit = 1;
    Int_t tofpid=0;
    if((status & AliESDtrack::kTOFpid) == AliESDtrack::kTOFpid) tofpid = 1;

    // ITS number of clusters
    UChar_t nclustersITS = track->GetITSNcls();
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
    Float_t ncrossedRowsTPC = track->GetTPCNCrossedRows();
    const TBits &sharedTPC = track->GetTPCSharedMap();
    for(Int_t ibit = 0; ibit < 160; ibit++) if(sharedTPC.TestBitNumber(ibit)) nclustersTPCshared++;
    // TRD clusters and tracklets
    UChar_t nclustersTRD = track->GetTRDncls();
    UChar_t ntrackletsTRDPID = track->GetTRDntrackletsPID();
    Int_t   nslicesTRD = track->GetDetPid()->GetTRDnSlices();
    Int_t   chi2TRD = track->GetTRDchi2();
    AliAODPid* aodpid= track->GetDetPid();
    Double_t* arraytrdsignals;
    arraytrdsignals=aodpid->GetTRDslices();
    Int_t nslicetemp=0;
    Int_t trdlayer[6];
    for(Int_t iplane = 0; iplane < 6; iplane++){
	    trdlayer[iplane]=0;
    }

    for(Int_t iplane = 0; iplane < 6; iplane++){
	    nslicetemp=0;
	    for(Int_t b=(iplane*8);b<((iplane*8)+8);b++){
	      if(ntrackletsTRDPID>0){
		      if(arraytrdsignals[b]>0.001) nslicetemp++;
	      }
	    }
	    if(nslicetemp > 0) trdlayer[iplane]=1;
    }

    // ITS and TRD acceptance maps
    UChar_t itsPixel = track->GetITSClusterMap();
    Bool_t statusL0 = kFALSE;
    if(TESTBIT(itsPixel, 0)) statusL0 = kTRUE; 
    Bool_t statusL1 = kFALSE;
    if(TESTBIT(itsPixel, 1)) statusL1 = kTRUE; 

    // HFE DCA
    Float_t dcaxy = -999.;
    Float_t dcaz = -999.;
    fExtraCuts->GetImpactParameters((AliVTrack *)track,dcaxy,dcaz);

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
    
    // ID track to see if negative
    Int_t id = track->GetID();

    // Double counted
    Int_t doublec = 0;
    for(Int_t l=0; l < counterdc; l++){
      Int_t iTrack2 = arraytrack.At(l);
      if(iTrack2==id) doublec=1;
    }
    //printf("Doublec %d\n",doublec);

    // Add the id at this place
    arraytrack.AddAt(id,counterdc);
    counterdc++;

    // Filter
    Int_t filter[20];
    for(Int_t k=0; k<20; k++) {
      filter[k]=0;
      Int_t u = 1<<k;     
      if((track->TestFilterBit(u))) {
	      filter[k] = 1;
      } 
    }
    Int_t filter0 = filter[0];
    Int_t filter1 = filter[1];
    Int_t filter2 = filter[2];
    Int_t filter3 = filter[3];
    Int_t filter4 = filter[4];
    Int_t filter5 = filter[5];
    Int_t filter6 = filter[6];
    Int_t filter7 = filter[7];
    Int_t filter8 = filter[8];
    Int_t filter9 = filter[9];
    Int_t filter10 = filter[10];
    Int_t filter11 = filter[11];
    Int_t filter12 = filter[12];
    Int_t filter13 = filter[13];
    Int_t filter14 = filter[14];
    Int_t filter15 = filter[15];
    Int_t filter16 = filter[16];
    Int_t filter17 = filter[17];
    Int_t filter18 = filter[18];
    Int_t filter19 = filter[19];

    // EMCAL cluster
    Double_t emcalEnergyOverP = -1.,
             emcalNSigma = -1.,
             showershape[4] = {0.,0.,0.,0.};
    if(track->IsEMCAL()){
      emcalNSigma = pid->NumberOfSigmasEMCAL(track, AliPID::kElectron, emcalEnergyOverP, &showershape[0]);
    }

    Int_t eventnb = fEventNumber;

    //printf("track\n");

    signal = 0;
    if(mcthere){
      Int_t label = TMath::Abs(track->GetLabel());
      if(label && label < fAODArrayMCInfo->GetEntriesFast())
        mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(label));
      if(!mctrack) continue;
      if(fTrackCuts->CheckParticleCuts(static_cast<UInt_t>(AliHFEcuts::kStepMCGenerated), mctrack)) signal = 1;
      // Kinematics
      chargemc = mctrack->Charge() > 0. ? 1. : -1.;
      momentummc = mctrack->P() * chargemc;
      transversemomentummc = mctrack->Pt() * chargemc;
      etamc = mctrack->Eta();
      phimc = mctrack->Phi();
      pdg = mctrack->GetPdgCode();
      
      // Get Production Vertex in radial direction
      vx = mctrack->Xv();
      vy = mctrack->Yv(); 
      vz = mctrack->Zv(); 
      eR = TMath::Sqrt(vx*vx+vy*vy);
      
      // Get Mother PDG code of the particle
      /*
      Int_t motherPdg = 0;
      Int_t motherlabel = TMath::Abs(mctrack->GetMother());
      if(motherlabel >= 0 && motherlabel < fAODArrayMCInfo->GetEntriesFast()){
        AliAODMCParticle *mother = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(motherlabel));
        if(mother) motherPdg = mother->GetPdgCode();
      }
      */
      
      // derive source
      source = 5;
      if(fSignalCuts->IsCharmElectron(mctrack)) source = 0;
      else if(fSignalCuts->IsBeautyElectron(mctrack)) source = 1;
      else if(fSignalCuts->IsGammaElectron(mctrack)) source = 2;
      else if(fSignalCuts->IsNonHFElectron(mctrack)) source = 3;
      else if(TMath::Abs(pdg) == 11) source = 4;
      else source = 5;
      
    }

    
    // Fill Tree
    //printf("Fill\n");
    (*fDebugTree) << "PIDdebug"
                  << "centrality="          << centrality
                  << "MBtrigger="           << isMBTrigger 
                  << "CentralTrigger="      << isCentralTrigger
                  << "SemicentralTrigger="  << isSemicentralTrigger
                  << "EMCALtrigger="        << isEMCALTrigger
		              << "run="                 << run
		              << "eventnb="             << eventnb
		              << "vx="                  << vtx[0]
                  << "vy="                  << vtx[1]
                  << "vz="                  << vtx[2]
		              << "ncontrib="            << ncontrib
		              << "id="                  << id
		              << "dc="                  << doublec
                  << "p="                   << momentum
		              << "ptpc="                << momentumTPC
		              << "pt="                  << transversemomentum
		              << "eta="                 << eta
                  << "phi="                 << phi
		              << "itsrefit="            << itsrefit
		              << "tpcrefit="            << tpcrefit
		              << "tofpid="              << tofpid
		              << "nclustersTPC="        << nclustersTPCfit
		              << "nclustersTPCall="     << nclustersTPCall
                  << "nclustersTPCPID="     << nclustersTPCPID
                  << "nclustersTPCshared="  << nclustersTPCshared
                  << "ncrossedRowsTPC="     << ncrossedRowsTPC
                  << "clusterRatioTPC="     << clusterRatioTPCfit
		              << "clusterRatioTPCall="  << clusterRatioTPCall
                  << "nclustersITS="        << nclustersITS
	                << "nclustersTRD="        << nclustersTRD
      	          << "ntrackletsTRD="       << ntrackletsTRDPID
	                << "nslicesTRD="          << nslicesTRD
                  << "trd0="                << trdlayer[0]
                  << "trd1="                << trdlayer[1]
                  << "trd2="                << trdlayer[2]
                  << "trd3="                << trdlayer[3]
                  << "trd4="                << trdlayer[4]
                  << "trd5="                << trdlayer[5]
                  << "chi2TRD="             << chi2TRD
		              << "statusITS0="          << statusL0
                  << "statusITS1="          << statusL1
 		              << "TOFsigmaEl="          << nSigmaTOF
                  << "TPCsigmaEl="          << nSigmaTPC
		              << "TPCsigmaElcorr="      << nSigmaTPCcorr
                  << "TPCdEdx="             << tPCdEdx
                  << "TPCdEdxcorr="         << tPCdEdxcorr
		              << "dcaR="                << dcaxy
                  << "EMCALEOP="            << emcalEnergyOverP
                  << "EMCALsigmaEl="        << emcalNSigma
                  << "showershape0="        << showershape[0]
                  << "showershape1="        << showershape[1]
                  << "showershape2="        << showershape[2]
                  << "showershape3="        << showershape[3]
                  << "dcaZ="                << dcaz
		              << "kinkdaughter="        << kink
		              << "kinkmother="          << kinkmotherpass
		              << "nbofmotherkink="      << numberofmotherkink
		              << "filter0="             << filter0
		              << "filter1="             << filter1
		              << "filter2="             << filter2
		              << "filter3="             << filter3
		              << "filter4="             << filter4
		              << "filter5="             << filter5
		              << "filter6="             << filter6
		              << "filter7="             << filter7
		              << "filter8="             << filter8
		              << "filter9="             << filter9
		              << "filter10="            << filter10
		              << "filter11="            << filter11
		              << "filter12="            << filter12
		              << "filter13="            << filter13
		              << "filter14="            << filter14
		              << "filter15="            << filter15
		              << "filter16="            << filter16
		              << "filter17="            << filter17
		              << "filter18="            << filter18
		              << "filter19="            << filter19
		              << "mcp="                 << momentummc
                  << "mcpt="                << transversemomentummc
		              << "mceta="               << etamc
                  << "mcphi="               << phimc
                  << "mcpdg="               << pdg
		              << "source="              << source
		              << "px="                  << vx
                  << "py="                  << vy
                  << "pz="                  << vz
		              << "eR="                  << eR
		              << "mccharge="            << chargemc
		              << "signal="              << signal
		              << "\n";

 
    //printf("after\n");

  }
  
  fEventNumber++;


}


void AliHFEdebugTreeTaskAOD::SetFileName(const char *filename){ fFilename = filename; }
void AliHFEdebugTreeTaskAOD::Terminate(Option_t *){

  if(fDebugTree) delete fDebugTree;
  fDebugTree=0x0;

}

