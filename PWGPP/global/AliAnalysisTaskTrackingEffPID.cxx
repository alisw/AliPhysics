
// ROOT includes
#include <TAxis.h>
#include <TChain.h>
#include <TF1.h>
#include <TList.h>
#include <TMath.h>
#include <TParticle.h>
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
  fIsAA{false},
  fFilterBit{4},
  fTrackCuts{0x0},
  fOutputList{0x0},
  fHistNEvents{0x0}
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
}

/// Standard destructor
///
AliAnalysisTaskTrackingEffPID::~AliAnalysisTaskTrackingEffPID(){
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  if (fOutputList) delete fOutputList;
  if (fTrackCuts) delete fTrackCuts;
}

/// This function creates all the histograms and all the objects in general used during the analysis
/// \return void
///
void AliAnalysisTaskTrackingEffPID::UserCreateOutputObjects() {

  fOutputList = new TList();
  fOutputList->SetOwner(true);

  fHistNEvents = new TH1F("hNEvents", "Number of processed events",3,-0.5,2.5);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"All events");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"MC selected");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"Reco Selected");
  fOutputList->Add(fHistNEvents);


  TString axTit[5]={"#eta","#varphi","#it{p}_{T} (GeV/#it{c})","N_{tracklets}","z_{vertex} (cm)"};
  const int nPtBins=32;
  const int nMultBins=8;
  int nbins[5]={10,18,nPtBins,nMultBins,4};
  double xmin[5]={-1.,0.,0.,0,-10.};
  double xmax[5]={1.,2*TMath::Pi(),30.,200.,10.};
  TString charge[2] = {"pos","neg"};
  double ptBins[nPtBins+1] = {0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.50,
			      0.60,0.70,0.80,0.90,1.00,1.25,1.50,1.75,2.00,2.50,
			      3.00,3.50,4.00,4.50,5.00,6.00,7.00,8.00,10.0,12.0,
			      16.0,20.0,30.0};
  double multBins[nMultBins+1] = {0.,5.,10.,20.,30.,40.,50.,80.,200.};
  if(fIsAA){
    multBins[0]=0.;
    multBins[1]=100.;
    multBins[2]=500.;
    multBins[3]=1000.;
    multBins[4]=2000.;
    multBins[5]=3000.;
    multBins[6]=4000.;
    multBins[7]=5000.;
    multBins[8]=10000.;
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

      fGenerated[iSpecies][iCharge]->GetAxis(3)->Set(nMultBins,multBins);
      fGeneratedEvSel[iSpecies][iCharge]->GetAxis(3)->Set(nMultBins,multBins);
      fReconstructed[iSpecies][iCharge]->GetAxis(3)->Set(nMultBins,multBins);
      fReconstructedTOF[iSpecies][iCharge]->GetAxis(3)->Set(nMultBins,multBins);
      fReconstructedPID[iSpecies][iCharge]->GetAxis(3)->Set(nMultBins,multBins);

      fOutputList->Add(fGenerated[iSpecies][iCharge]);
      fOutputList->Add(fGeneratedEvSel[iSpecies][iCharge]);
      fOutputList->Add(fReconstructed[iSpecies][iCharge]);
      fOutputList->Add(fReconstructedTOF[iSpecies][iCharge]);
      fOutputList->Add(fReconstructedPID[iSpecies][iCharge]);
    }
  }
  fEventCut.AddQAplotsToList(fOutputList);

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

  double zMCVertex =99999;
  int nTracklets = 0;
  if(isAOD){
    AliAODMCHeader *mcHeader = dynamic_cast<AliAODMCHeader*>(fInputEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if (!mcHeader) {
      AliError("Could not find MC Header in AOD");
      return;
    }
    zMCVertex = mcHeader->GetVtxZ();
    AliAODTracklets *mult=((AliAODEvent*)fInputEvent)->GetTracklets();
    if(mult) nTracklets=mult->GetNumberOfTracklets();
  }else{
    const AliVVertex* mcVert=fMCEvent->GetPrimaryVertex();
    if(!mcVert){
      AliError("Generated vertex not available");
      return;
    }
    zMCVertex=mcVert->GetZ();
    const AliMultiplicity *mult = ((AliESDEvent*)fInputEvent)->GetMultiplicity();
    if(mult) nTracklets=mult->GetNumberOfTracklets();
  }

  if(zMCVertex<fEventCut.fMinVtz || zMCVertex>fEventCut.fMaxVtz){
    PostData(1, fOutputList);
    return;
  }
  fHistNEvents->Fill(1);

  const AliVVertex* vtTrc = fInputEvent->GetPrimaryVertex();
  double pos[3],cov[6];
  vtTrc->GetXYZ(pos);
  vtTrc->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos,cov,100.,100);
  double magField = fInputEvent->GetMagneticField();

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* handl = (AliInputEventHandler*)mgr->GetInputEventHandler();
  AliPIDResponse* pid = handl->GetPIDResponse();
  if (!pid) {
    AliFatal("Missing PID response. Did you attach the AliPIDresponseTask to your analysis?");
  }

  for (int iMC = 0; iMC < fMCEvent->GetNumberOfTracks(); ++iMC) {
    AliVParticle *part = (AliVParticle*)fMCEvent->GetTrack(iMC);
    if (!part->IsPhysicalPrimary()) continue;
    double arrayForSparse[5]={part->Eta(),part->Phi(),part->Pt(),(double)nTracklets,zMCVertex};
    const int pdg = std::abs(part->PdgCode());
    const int iCharge = part->Charge() > 0 ? 1 : 0;
    for (int iSpecies = 0; iSpecies < AliPID::kSPECIESC; ++iSpecies) {
      if (pdg == AliPID::ParticleCode(iSpecies)) {
        fGenerated[iSpecies][iCharge]->Fill(arrayForSparse);
	if(eventAccepted) fGeneratedEvSel[iSpecies][iCharge]->Fill(arrayForSparse);
        break;
      }
    }
  }

  if (!eventAccepted) {
    PostData(1, fOutputList);
    return;
  }
  fHistNEvents->Fill(2);


  for (int iT = 0; iT < (int)fInputEvent->GetNumberOfTracks(); ++iT) {
    /// Get the track and do the minimal cuts
    AliVTrack *track = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(iT));

    if(!isAOD){
      AliESDtrack *esdtrack = dynamic_cast<AliESDtrack*>(track); 
      if(fTrackCuts && !fTrackCuts->AcceptTrack(esdtrack)) continue;
    }else{
      if(track->GetID() < 0) continue;
      AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(track); 
      if(fFilterBit>=0 && !aodtrack->TestFilterBit(BIT(fFilterBit))) continue;
      if(fTrackCuts && fUseTrackCutsForAOD){
	bool accept=ConvertAndSelectAODTrack(aodtrack,vESD,magField);
	if(!accept) continue;
      }
    }

    AliVParticle *mcPart  = (AliVParticle*)fMCEvent->GetTrack(TMath::Abs(track->GetLabel()));
    if(!mcPart) continue;
    if (!mcPart->IsPhysicalPrimary()) continue;
    const int iCharge = mcPart->Charge() > 0 ? 1 : 0;
    int iSpecies = -1;
    for (int iS = 0; iS < AliPID::kSPECIESC; ++iS) {
      if (std::abs(mcPart->PdgCode()) == AliPID::ParticleCode(iS)) {
        iSpecies = iS;
        break;
      }
    }
    if (iSpecies < 0) continue;

    const double pt = fUseGeneratedKine ? mcPart->Pt() : track->Pt() * AliPID::ParticleCharge(iSpecies);
    const double eta = fUseGeneratedKine ? mcPart->Eta() : track->Eta();
    const double phi = fUseGeneratedKine ? mcPart->Phi() : track->Phi();
    double arrayForSparseData[5]={eta,phi,pt,(double)nTracklets,zMCVertex};
    bool TPCpid = std::abs(pid->NumberOfSigmasTPC(track, static_cast<AliPID::EParticleType>(iSpecies))) < 3;
    bool hasTOF = HasTOF(track);
    bool TOFpid = std::abs(pid->NumberOfSigmasTOF(track, static_cast<AliPID::EParticleType>(iSpecies))) < 3;

    fReconstructed[iSpecies][iCharge]->Fill(arrayForSparseData);
    if(hasTOF) fReconstructedTOF[iSpecies][iCharge]->Fill(arrayForSparseData);
    if(TPCpid && TOFpid) fReconstructedPID[iSpecies][iCharge]->Fill(arrayForSparseData);

  } // End track loop

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
