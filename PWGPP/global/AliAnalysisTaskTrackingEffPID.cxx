
// ROOT includes
#include <TAxis.h>
#include <TChain.h>
#include <TF1.h>
#include <TList.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TRandom3.h>

// ALIROOT includes
#include "AliAnalysisManager.h"
#include "AliPIDResponse.h"
#include "AliPDG.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVParticle.h"
#include "AliMCEvent.h"
#include "AliInputEventHandler.h"
#include "AliVEventHandler.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"
#include "AliAODMCParticle.h"
#include "AliVParticle.h"
#include "AliAODVertex.h"
#include "AliAODMCHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenPythiaEventHeader.h"

#include "AliAnalysisTaskTrackingEffPID.h"

/**************************************************************************
 * Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
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

//*************************************************************************
// Implementation of class AliAnalysisTaskTrackingEffPID
// AliAnalysisTaskSE to compute tracking and PID efficiencies for 
//  different particle species
//
// Authors: 
//          M. Puccio
//          F. Prino
//          
//*************************************************************************


///\cond CLASSIMP
ClassImp(AliAnalysisTaskTrackingEffPID);
///\endcond

AliAnalysisTaskTrackingEffPID::AliAnalysisTaskTrackingEffPID() :
  AliAnalysisTaskSE("TrackingEffPID"),
  fEventCut{false},
  fUseTrackCutsForAOD{false},
  fUseGeneratedKine{true},
  fRejectPileupParticles{true},
  fRejectTracksOfPileupPart{false},
  fPrimarySelectionOpt{1},
  fMultEstimator{0},
  fIsAA{false},
  fFilterBit{4},
  fTrackCuts{0x0},
  fSelectOnGenerator{false},
  fGenerToKeep{""},
  fGenerToExclude{""},
  fKeepOnlyInjected{false},
  fKeepOnlyUE{false},
  fUseImpPar{false},
  fUseLocDen{false},
  fDeltaRcut{0.2},
  fMaxTracksInCone{-1.},
  fSelectPtHardRange{false},
  fMinPtHard{0.},
  fMaxPtHard{99999.},
  fOutputList{0x0},
  fListCuts{0x0},
  fHistNEvents{0x0},
  fHistNParticles{0x0},
  fHistNTracks{0x0},
  fHistPileupTagAOD{0x0},
  hHistXsecVsPtHard{0x0}
{
  // default: use the filter bit 4 cuts
  fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
  fTrackCuts->SetMaxDCAToVertexXY(2.4);
  fTrackCuts->SetMaxDCAToVertexZ(3.2);
  fTrackCuts->SetDCAToVertex2D(kTRUE);
  // default: physics selection with loose trigger mask request, pileup disabled
  fEventCut.OverrideAutomaticTriggerSelection(AliVEvent::kMB|AliVEvent::kINT7);
  fEventCut.OverridePileUpCuts(999999,999.,999.,3.,3.,kTRUE); 
  fEventCut.fPileUpCutMV=kFALSE;

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

/// Standard destructor
///
AliAnalysisTaskTrackingEffPID::~AliAnalysisTaskTrackingEffPID(){
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  if (fOutputList) delete fOutputList;
  if (fListCuts) delete fListCuts;
  if (fTrackCuts) delete fTrackCuts;
}

/// This function creates all the histograms and all the objects in general used during the analysis
/// \return void
///
void AliAnalysisTaskTrackingEffPID::UserCreateOutputObjects() {

  fOutputList = new TList();
  fOutputList->SetOwner(true);

  fHistNEvents = new TH1F("hNEvents", "Number of processed events",4,-0.5,3.5);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"All events");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"Generator name selected");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"MC selected");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"Reco Selected");
  fOutputList->Add(fHistNEvents);

  fHistNParticles = new TH1D("hNParticles", "Number of particles",5,-0.5,4.5);
  fHistNParticles->GetXaxis()->SetBinLabel(1,"All particles");
  fHistNParticles->GetXaxis()->SetBinLabel(2,"Pileup events");
  fHistNParticles->GetXaxis()->SetBinLabel(3,"Trigger event");
  fHistNParticles->GetXaxis()->SetBinLabel(4,"Phys. Primary");
  fHistNParticles->GetXaxis()->SetBinLabel(5,"Injected/UE sel.");
  fOutputList->Add(fHistNParticles);

  fHistNTracks = new TH1D("hNTracks", "Number of tracks",7,-0.5,6.5);
  fHistNTracks->GetXaxis()->SetBinLabel(1,"All tracks");
  fHistNTracks->GetXaxis()->SetBinLabel(2,"After track sel.");
  fHistNTracks->GetXaxis()->SetBinLabel(3,"Pileup events");
  fHistNTracks->GetXaxis()->SetBinLabel(4,"Trigger event");
  fHistNTracks->GetXaxis()->SetBinLabel(5,"Phys. Primary");
  fHistNTracks->GetXaxis()->SetBinLabel(6,"Injected/UE sel.");
  fHistNTracks->GetXaxis()->SetBinLabel(7,"Species sel.");
  fOutputList->Add(fHistNTracks);

  fHistPileupTagAOD = new TH2D("hPileupTagAOD"," OOB pileup particles ; tag with AliMCEvent ; tag with AliAODMCHeader",2,-0.5,1.5,2,-0.5,1.5);
  fOutputList->Add(fHistPileupTagAOD);
  
  hHistXsecVsPtHard = new TH1D("hXsecVsPtHard", " ; pthard (GeV/c) ; Xsec", 200,0.,100.);
  fOutputList->Add(hHistXsecVsPtHard);
  
  TString axTit[5]={"#eta","#varphi","#it{p}_{T} (GeV/#it{c})","Multiplicity","z_{vertex} (cm)"};
  if(fMultEstimator==0)  axTit[3]="N_{tracklets}";
  else if(fMultEstimator==1) axTit[3]="N_{contributors}";
  else if(fMultEstimator==2) axTit[3]="N_{TPCITStracks}";
  else if(fMultEstimator==3) axTit[3]="N_{TPCtracks}";
  else if(fMultEstimator==4) axTit[3]="N_{TPCclusters}/1000";
  const int nPtBins=32;
  const int nMultBins=10;
  int nbins[5]={10,18,nPtBins,nMultBins,4};
  double xmin[5]={-1.,0.,0.,0,-10.};
  double xmax[5]={1.,2*TMath::Pi(),30.,200.,10.};
  TString charge[2] = {"pos","neg"};
  double ptBins[nPtBins+1] = {0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.50,
                              0.60,0.70,0.80,0.90,1.00,1.25,1.50,1.75,2.00,2.50,
                              3.00,3.50,4.00,4.50,5.00,6.00,7.00,8.00,10.0,12.0,
                              16.0,20.0,30.0};
  double multBins[nMultBins+1] = {0.,5.,10.,20.,30.,40.,50.,60.,80.,100.,200.};
  if(fIsAA){
    multBins[0]=0.;
    multBins[1]=100.;
    multBins[2]=500.;
    multBins[3]=1000.;
    multBins[4]=1500.;
    multBins[5]=2000.;
    multBins[6]=2500.;
    multBins[7]=3000.;
    multBins[8]=4000.;
    multBins[9]=5000.;
    multBins[10]=7500.;
    if(fMultEstimator==3) for(Int_t jm=0; jm<=nMultBins; jm++) multBins[jm]*=2.;
  }
  xmax[3]=multBins[nMultBins];
  if(fUseImpPar && !fUseLocDen){
    // use impact parameter instead of multiplicity
    axTit[3]="b (fm)";
    nbins[3]=15;
    xmin[3]=0.;
    xmax[3]=15.;
  }
  if(fUseLocDen){
    // use local track density instead of multiplicity
    axTit[3]="tracks in cone";
    if(fIsAA){
      nbins[3]=25;
      xmin[3]=0.;
      if(fMaxTracksInCone>0) xmax[3]=fMaxTracksInCone;
      else xmax[3]=125.;
    }else{
      nbins[3]=15;
      xmin[3]=0.;
      if(fMaxTracksInCone>0) xmax[3]=fMaxTracksInCone;
      else xmax[3]=45.;
    }
  }
  
  for (int iSpecies = 0; iSpecies < AliPID::kSPECIESC; iSpecies++) {
    for (int iCharge = 0; iCharge < 2; ++iCharge) {
      fGenerated[iSpecies][iCharge] = new THnSparseF(Form("hGen_%s_%s",AliPID::ParticleShortName(iSpecies),charge[iCharge].Data()),
                                                     "Generated particles",5,nbins,xmin,xmax);
      fGeneratedEvSel[iSpecies][iCharge] = new THnSparseF(Form("hGenEvSel_%s_%s",AliPID::ParticleShortName(iSpecies),charge[iCharge].Data()),
                                                          "Generated particles in selected events",5,nbins,xmin,xmax);
      fReconstructed[iSpecies][iCharge] = new THnSparseF(Form("hReconstructed_%s_%s",AliPID::ParticleShortName(iSpecies),charge[iCharge].Data()),
                                                          "Reconstructed tracks",5,nbins,xmin,xmax);
      fReconstructedTOF[iSpecies][iCharge] = new THnSparseF(Form("hReconstructedTOF_%s_%s",AliPID::ParticleShortName(iSpecies),charge[iCharge].Data()),
                                                            "Reconstructed tracks with TOF",5,nbins,xmin,xmax);
      fReconstructedPID[iSpecies][iCharge] = new THnSparseF(Form("hReconstructedPID_%s_%s",AliPID::ParticleShortName(iSpecies),charge[iCharge].Data()),
                                                            "Reconstructed tracks with PID",5,nbins,xmin,xmax);
      for(int iax=0; iax<5; iax++){
        fGenerated[iSpecies][iCharge]->GetAxis(iax)->SetTitle(axTit[iax].Data());
        fGeneratedEvSel[iSpecies][iCharge]->GetAxis(iax)->SetTitle(axTit[iax].Data());
        fReconstructed[iSpecies][iCharge]->GetAxis(iax)->SetTitle(axTit[iax].Data());
        fReconstructedTOF[iSpecies][iCharge]->GetAxis(iax)->SetTitle(axTit[iax].Data());
        fReconstructedPID[iSpecies][iCharge]->GetAxis(iax)->SetTitle(axTit[iax].Data());
      }
      fGenerated[iSpecies][iCharge]->GetAxis(2)->Set(nPtBins,ptBins);
      fGeneratedEvSel[iSpecies][iCharge]->GetAxis(2)->Set(nPtBins,ptBins);
      fReconstructed[iSpecies][iCharge]->GetAxis(2)->Set(nPtBins,ptBins);
      fReconstructedTOF[iSpecies][iCharge]->GetAxis(2)->Set(nPtBins,ptBins);
      fReconstructedPID[iSpecies][iCharge]->GetAxis(2)->Set(nPtBins,ptBins);

      if(!fUseImpPar && !fUseLocDen){
        fGenerated[iSpecies][iCharge]->GetAxis(3)->Set(nMultBins,multBins);
        fGeneratedEvSel[iSpecies][iCharge]->GetAxis(3)->Set(nMultBins,multBins);
        fReconstructed[iSpecies][iCharge]->GetAxis(3)->Set(nMultBins,multBins);
        fReconstructedTOF[iSpecies][iCharge]->GetAxis(3)->Set(nMultBins,multBins);
        fReconstructedPID[iSpecies][iCharge]->GetAxis(3)->Set(nMultBins,multBins);
      }
      
      fOutputList->Add(fGenerated[iSpecies][iCharge]);
      fOutputList->Add(fGeneratedEvSel[iSpecies][iCharge]);
      fOutputList->Add(fReconstructed[iSpecies][iCharge]);
      fOutputList->Add(fReconstructedTOF[iSpecies][iCharge]);
      fOutputList->Add(fReconstructedPID[iSpecies][iCharge]);
    }
  }
  fEventCut.AddQAplotsToList(fOutputList);

  fListCuts = new TList();
  fListCuts->SetOwner();
  if(fTrackCuts){
    AliESDtrackCuts* ttosave=new AliESDtrackCuts(*fTrackCuts);
    fListCuts->Add(ttosave);
    TH1F* hAODCuts=new TH1F("hAODCuts","",2,0.,2.);
    hAODCuts->GetXaxis()->SetBinLabel(1,"filter bit");
    hAODCuts->SetBinContent(1,fFilterBit);
    hAODCuts->GetXaxis()->SetBinLabel(2,"use track cuts for AOD");
    hAODCuts->SetBinContent(2,fUseTrackCutsForAOD);
    fListCuts->Add(hAODCuts);
  }
  PostData(2, fListCuts);

  PostData(1,fOutputList);
}

/// This is the function that is evaluated for each event. The analysis code stays here.
///
/// \param options Deprecated parameter
/// \return void
///
void AliAnalysisTaskTrackingEffPID::UserExec(Option_t *){

  AliVEvent *ev = fInputEvent;
  if(!ev){
    AliFatal("NO EVENT FOUND!");
    return;    
  }
  if(!fMCEvent) {
    AliFatal("NO MC INFO FOUND");
    return;
  }

  bool isAOD = fInputEvent->IsA()->InheritsFrom("AliAODEvent");
  bool eventAccepted = fEventCut.AcceptEvent(fInputEvent);

  fHistNEvents->Fill(0);

  // check the generator name
  TList *lh=0x0;
  double imppar=-999.;
  AliAODMCHeader *aodMcHeader = 0x0;
  if(fSelectOnGenerator || fKeepOnlyInjected || fKeepOnlyUE || fUseImpPar || fSelectPtHardRange){
    if(isAOD){
      aodMcHeader = dynamic_cast<AliAODMCHeader*>(fInputEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
      lh=aodMcHeader->GetCocktailHeaders();
    }else{
      TString genname=fMCEvent->GenEventHeader()->ClassName();
      if(genname.Contains("CocktailEventHeader")){
        AliGenCocktailEventHeader *cockhead=(AliGenCocktailEventHeader*)fMCEvent->GenEventHeader();
        lh=cockhead->GetHeaders();
      }
    }
    if(fSelectPtHardRange && lh){
      Int_t nh=lh->GetEntries();
      for(Int_t i=0;i<nh;i++){
        AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(i);
        TString genname=gh->GetName();
        if(genname.Contains("ythia") || genname.Contains("YTHIA")){
          AliGenPythiaEventHeader* pyth=(AliGenPythiaEventHeader*)lh->At(i);
          double ptha=pyth->GetPtHard();
          double xsec=pyth->GetXsection();
          if(ptha<fMinPtHard || ptha>fMaxPtHard) return;
          hHistXsecVsPtHard->SetBinContent(hHistXsecVsPtHard->GetXaxis()->FindBin(ptha),xsec);
        }
      }
    }
    if(fUseImpPar && lh){
      Int_t nh=lh->GetEntries();
      for(Int_t i=0;i<nh;i++){
        AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(i);
        TString genname=gh->GetName();
        if(genname.Contains("hijing") || genname.Contains("Hijing")){
          AliGenHijingEventHeader* hijh=(AliGenHijingEventHeader*)lh->At(i);
          imppar=hijh->ImpactParameter();
        }
      }
    }
    if(fSelectOnGenerator && lh){
      Bool_t keep=kTRUE;
      if(fGenerToExclude.Length()==0) keep=kFALSE;
      Int_t nh=lh->GetEntries();
      for(Int_t i=0;i<nh;i++){
        AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(i);
        TString genname=gh->GetName();
        if(fGenerToKeep.Length()>0 && genname.Contains(fGenerToKeep.Data())) keep=kTRUE;
        if(fGenerToExclude.Length()>0 && genname.Contains(fGenerToExclude.Data())) keep=kFALSE;
      }
      if(!keep) return;
    }
  }

  fHistNEvents->Fill(1);

  double zMCVertex =99999;
  int nTracklets = 0;
  TClonesArray *aodMcArray = 0x0;
  int nTPCclusters=0;
  if(isAOD){
    aodMcArray = dynamic_cast<TClonesArray*>(fInputEvent->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
    aodMcHeader = dynamic_cast<AliAODMCHeader*>(fInputEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if (!aodMcHeader) {
      AliError("Could not find MC Header in AOD");
      return;
    }
    zMCVertex = aodMcHeader->GetVtxZ();
    AliAODTracklets *mult=((AliAODEvent*)fInputEvent)->GetTracklets();
    if(mult) nTracklets=mult->GetNumberOfTracklets();
    nTPCclusters=((AliAODEvent*)fInputEvent)->GetNumberOfTPCClusters();
  }else{
    const AliVVertex* mcVert=fMCEvent->GetPrimaryVertex();
    if(!mcVert){
      AliError("Generated vertex not available");
      return;
    }
    zMCVertex=mcVert->GetZ();
    const AliMultiplicity *mult = ((AliESDEvent*)fInputEvent)->GetMultiplicity();
    if(mult) nTracklets=mult->GetNumberOfTracklets();
    nTPCclusters=((AliESDEvent*)fInputEvent)->GetNumberOfTPCClusters();
  }

  if(zMCVertex<fEventCut.fMinVtz || zMCVertex>fEventCut.fMaxVtz){
    PostData(1, fOutputList);
    return;
  }
  fHistNEvents->Fill(2);

  const AliVVertex* vtTrc = fInputEvent->GetPrimaryVertex();
  double pos[3],cov[6];
  vtTrc->GetXYZ(pos);
  vtTrc->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos,cov,100.,100);
  double magField = fInputEvent->GetMagneticField();
  int nContrib = vtTrc->GetNContributors();;

  int nTracksTPCITS = 0;
  int nTracksTPC = 0;
  TNtuple* trEtaPhiMap = new TNtuple("trEtaPhiMap", "tracks", "eta:phi");
  for (int iT = 0; iT < (int)fInputEvent->GetNumberOfTracks(); ++iT) {
    /// count tracks passing ITS TPC selections
    AliVTrack *track = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(iT));
    int cluITS=track->GetNcls(0);
    int cluTPC=track->GetNcls(1);
    if(track->GetStatus()&AliESDtrack::kTPCin && cluTPC>70) nTracksTPC++;
    if(track->GetStatus()&AliESDtrack::kITSrefit && cluITS>3 && cluTPC>70) nTracksTPCITS++;
    if(fUseLocDen && track->GetStatus()&AliESDtrack::kTPCin && track->GetID()>=0) trEtaPhiMap->Fill(track->Eta(),track->Phi());
  }

  double multEstim=nTracklets;
  if(fMultEstimator==1) multEstim=nContrib;
  else if(fMultEstimator==2) multEstim=nTracksTPCITS;
  else if(fMultEstimator==3) multEstim=nTracksTPC;
  else if(fMultEstimator==4) multEstim=nTPCclusters/1000.;

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* handl = (AliInputEventHandler*)mgr->GetInputEventHandler();
  AliPIDResponse* pid = handl->GetPIDResponse();
  if (!pid) {
    AliFatal("Missing PID response. Did you attach the AliPIDresponseTask to your analysis?");
  }

  for (int iMC = 0; iMC < fMCEvent->GetNumberOfTracks(); ++iMC) {
    AliVParticle *part = (AliVParticle*)fMCEvent->GetTrack(iMC);
    fHistNParticles->Fill(0);
    if(aodMcHeader && aodMcArray){
      Bool_t isPil1=AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iMC,fMCEvent);
      Bool_t isPil2=AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iMC,aodMcHeader,aodMcArray);
      fHistPileupTagAOD->Fill(isPil1,isPil2);
    }
    if(AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iMC,fMCEvent)){
      fHistNParticles->Fill(1);
      if(fRejectPileupParticles) continue;
    }else{
      fHistNParticles->Fill(2);
    }
    if(fPrimarySelectionOpt==1 && !part->IsPhysicalPrimary()) continue;
    if(fPrimarySelectionOpt==2){
      // primary particle selection based on origin of particle
      double pRad2=part->Xv()*part->Xv()+part->Yv()*part->Yv();
      double distz=TMath::Abs(part->Zv()-zMCVertex);
      if(pRad2>8 || distz>1) continue;
    }
    fHistNParticles->Fill(3);
    if(lh && (fKeepOnlyInjected || fKeepOnlyUE)){
      bool isInjected=IsInjectedParticle(iMC,lh);
      if(fKeepOnlyInjected && !isInjected) continue;
      if(fKeepOnlyUE && isInjected) continue;
    }
    fHistNParticles->Fill(4);
    double arrayForSparse[5]={part->Eta(),part->Phi(),part->Pt(),multEstim,zMCVertex};
    if(fUseImpPar) arrayForSparse[3]=imppar;
    const int pdg = std::abs(part->PdgCode());
    Int_t jPDG=-1;
    for (int iSpecies = 0; iSpecies < AliPID::kSPECIESC; ++iSpecies) {
      if (pdg == AliPID::ParticleCode(iSpecies)) {
        jPDG=iSpecies;
        break;
      }
    }
    if(jPDG>0){
      if(fUseLocDen) arrayForSparse[3]=GetLocalTrackDens(trEtaPhiMap,part->Eta(),part->Phi());
      const int iCharge = part->Charge() > 0 ? 1 : 0;
      fGenerated[jPDG][iCharge]->Fill(arrayForSparse);
      if(eventAccepted) fGeneratedEvSel[jPDG][iCharge]->Fill(arrayForSparse);
    }
  }

  if (!eventAccepted) {
    delete trEtaPhiMap;
    PostData(1, fOutputList);
    return;
  }
  fHistNEvents->Fill(3);


  for (int iT = 0; iT < (int)fInputEvent->GetNumberOfTracks(); ++iT) {
    /// Get the track and do the minimal cuts
    AliVTrack *track = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(iT));
    fHistNTracks->Fill(0);
    
    if(!isAOD){
      AliESDtrack *esdtrack = dynamic_cast<AliESDtrack*>(track); 
      if(fTrackCuts && !fTrackCuts->AcceptTrack(esdtrack)) continue;
    }else{
      AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(track); 
      if(fFilterBit<0 && aodtrack->GetID() < 0) continue;
      if(fFilterBit>=0 && !aodtrack->TestFilterBit(BIT(fFilterBit))) continue;
      if(fTrackCuts && fUseTrackCutsForAOD){
        bool accept=ConvertAndSelectAODTrack(aodtrack,vESD,magField);
        if(!accept) continue;
      }
    }
    fHistNTracks->Fill(1);
    
    int lab=TMath::Abs(track->GetLabel());
    AliVParticle *mcPart  = (AliVParticle*)fMCEvent->GetTrack(lab);
    if(!mcPart) continue;
    if(AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(lab,fMCEvent)){
      fHistNTracks->Fill(2);
      if(fRejectTracksOfPileupPart) continue;
    }else{
      fHistNTracks->Fill(3);
    }
    
    if (fPrimarySelectionOpt==1 && !mcPart->IsPhysicalPrimary()) continue;
    if(fPrimarySelectionOpt==2){
      // primary particle selection based on origin of particle
      double pRad2=mcPart->Xv()*mcPart->Xv()+mcPart->Yv()*mcPart->Yv();
      double distz=TMath::Abs(mcPart->Zv()-zMCVertex);
      if(pRad2>8 || distz>1) continue;
    }
    fHistNTracks->Fill(4);
    if(lh && (fKeepOnlyInjected || fKeepOnlyUE)){
      bool isInjected=IsInjectedParticle(lab,lh);
      if(fKeepOnlyInjected && !isInjected) continue;
      if(fKeepOnlyUE && isInjected) continue;
    }
    fHistNTracks->Fill(5);

    const int iCharge = mcPart->Charge() > 0 ? 1 : 0;
    int iSpecies = -1;
    for (int iS = 0; iS < AliPID::kSPECIESC; ++iS) {
      if (std::abs(mcPart->PdgCode()) == AliPID::ParticleCode(iS)) {
        iSpecies = iS;
        break;
      }
    }
    if (iSpecies < 0) continue;
    fHistNTracks->Fill(6);

    const double pt = fUseGeneratedKine ? mcPart->Pt() : track->Pt() * AliPID::ParticleCharge(iSpecies);
    const double eta = fUseGeneratedKine ? mcPart->Eta() : track->Eta();
    const double phi = fUseGeneratedKine ? mcPart->Phi() : track->Phi();
    double arrayForSparseData[5]={eta,phi,pt,multEstim,zMCVertex};
    if(fUseImpPar) arrayForSparseData[3]=imppar;
    if(fUseLocDen) arrayForSparseData[3]=GetLocalTrackDens(trEtaPhiMap,eta,phi);
    bool TPCpid = std::abs(pid->NumberOfSigmasTPC(track, static_cast<AliPID::EParticleType>(iSpecies))) < 3;
    bool hasTOF = HasTOF(track);
    bool TOFpid = std::abs(pid->NumberOfSigmasTOF(track, static_cast<AliPID::EParticleType>(iSpecies))) < 3;

    fReconstructed[iSpecies][iCharge]->Fill(arrayForSparseData);
    if(hasTOF) fReconstructedTOF[iSpecies][iCharge]->Fill(arrayForSparseData);
    if(TPCpid && TOFpid) fReconstructedPID[iSpecies][iCharge]->Fill(arrayForSparseData);

  } // End track loop
  delete trEtaPhiMap;

  //  Post output data.
  PostData(1,fOutputList);
}

//______________________________________________________________________________
/// Merge the output. Called once at the end of the query.
///
/// \return void
///
void AliAnalysisTaskTrackingEffPID::Terminate(Option_t *) {
  return;
}

//______________________________________________________________________________
/// This function checks whether a track has or has not a prolongation in TOF.
///
/// \param track Track that has to be checked
///
bool AliAnalysisTaskTrackingEffPID::HasTOF(AliVTrack *track) {
  const bool hasTOFout  = track->GetStatus() & AliVTrack::kTOFout;
  const bool hasTOFtime = track->GetStatus() & AliVTrack::kTIME;
  const bool hasGoodLength = track->GetIntegratedLength() > 350.;
  return hasTOFout && hasTOFtime && hasGoodLength;
}

//______________________________________________________________________________
/// This function convertes and AOD track into and ESD track and applires the cuts
///
/// \param aTrack Track that has to be selected
/// \param vESD primary vertex position
/// \param magField Magnetic Field value
bool AliAnalysisTaskTrackingEffPID::ConvertAndSelectAODTrack(AliAODTrack* aTrack, const AliESDVertex vESD, Double_t magField)
{

  AliESDtrack esdTrack(aTrack);
  esdTrack.SetTPCClusterMap(aTrack->GetTPCClusterMap());
  esdTrack.SetTPCSharedMap(aTrack->GetTPCSharedMap());
  esdTrack.SetTPCPointsF(aTrack->GetTPCNclsF());
  esdTrack.SetTPCNcls(aTrack->GetTPCNcls());
  Int_t nTPCclus=aTrack->GetNcls(1);
  Double_t chi2ndf=aTrack->Chi2perNDF();
  Double_t chi2tpc=999.;
  if(chi2ndf>0. && nTPCclus > 5){
    chi2tpc=Float_t(nTPCclus-5)*chi2ndf;
  }
  esdTrack.SetTPCchi2(chi2tpc);
  // needed to calculate the impact parameters
  Bool_t okDCA=esdTrack.RelateToVertex(&vESD,magField,99999.);
  if(!okDCA) return kFALSE;
  AliAODVertex* av=aTrack->GetProdVertex();
  if(av->GetType()==AliAODVertex::kKink) return kFALSE;
  return fTrackCuts->AcceptTrack(&esdTrack);
}

//______________________________________________________________________________
TString AliAnalysisTaskTrackingEffPID::GetGenerator(int label, TList *lh){
  /// get the name of the generator that produced a given particle
  
  Int_t nsumpart=0;
  Int_t nh=lh->GetEntries();
  for(Int_t i=0;i<nh;i++){
    AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(i);
    TString genname=gh->GetName();
    Int_t npart=gh->NProduced();
    if(label>=nsumpart && label<(nsumpart+npart)) return genname;
    nsumpart+=npart;
  }
  TString empty="";
  return empty;
}

//______________________________________________________________________________
bool AliAnalysisTaskTrackingEffPID::IsInjectedParticle(int lab, TList *lh){
  /// check if a particle is injected signal or hijing UE
  TString nameGen=GetGenerator(lab,lh);
  while(nameGen.IsWhitespace()){
    AliVParticle* mcpart=(AliVParticle*)fMCEvent->GetTrack(lab);
    if(!mcpart) break;
    int mother=fMCEvent->GetLabelOfParticleMother(lab);
    if(mother<0) break;
    lab=mother;
    nameGen=GetGenerator(lab,lh);
  }
  if(nameGen.IsWhitespace() || nameGen.Contains("ijing")) return kFALSE;
  else return kTRUE;
}
//______________________________________________________________________________
double AliAnalysisTaskTrackingEffPID::GetLocalTrackDens(TNtuple* trEtaPhiMap, double eta, double phi) const {
  /// count tracks in a cone around selected particle

  if(TMath::Abs(eta)>1) return -1;
  double nTracksInCone=0.;
  float etatr,phitr;
  trEtaPhiMap->SetBranchAddress("eta",&etatr);
  trEtaPhiMap->SetBranchAddress("phi",&phitr);
  double etamin=eta-fDeltaRcut;
  double etamax=eta+fDeltaRcut;
  double scalFac=1;
  if(fDeltaRcut<0.8){
    if(etamax>0.8){
      etamax=eta;
      scalFac=2.;
    }
    if(etamin<-0.8){
      etamin=eta;
      scalFac=2.;
    }
  }
  for (int iT = 0; iT < trEtaPhiMap->GetEntriesFast(); ++iT) {
    trEtaPhiMap->GetEvent(iT);
    if(etatr>etamin && etatr<etamax){
      double deltaEta=etatr-eta;
      double deltaPhi=phitr-phi;
      if(deltaPhi<-TMath::Pi()) deltaPhi+=2*TMath::Pi();
      else if(deltaPhi>TMath::Pi()) deltaPhi-=2*TMath::Pi();
      double deltaR2=deltaEta*deltaEta+deltaPhi*deltaPhi;
      if(deltaR2<fDeltaRcut*fDeltaRcut) nTracksInCone+=1.;
    }
  }
  return nTracksInCone*scalFac;
}
