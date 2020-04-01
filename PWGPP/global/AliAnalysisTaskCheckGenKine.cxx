#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliHeader.h"
#include "AliMultiplicity.h"
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TChain.h>
#include "AliGenerator.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliAnalysisTaskCheckGenKine.h"

/**************************************************************************
 * Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */ 

//*************************************************************************
// Implementation of class AliAnalysisTaskCheckGenKine
// AliAnalysisTask to check MC production at ESD+Kine level
// 
//
// Authors: F. Prino, prino@to.infn.it
//          
//*************************************************************************

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskCheckGenKine);
/// \endcond

//______________________________________________________________________________
AliAnalysisTaskCheckGenKine::AliAnalysisTaskCheckGenKine() : 
  AliAnalysisTaskSE("CheckGenKine"), 
  fOutput(0),
  fHistoNEvents(0),
  fHistoGenMult(0),
  fHistoVtxContrib(0),
  fHistoSPD3DVtxX(0),
  fHistoSPD3DVtxY(0),
  fHistoSPD3DVtxZ(0),
  fHistoSPDZVtxZ(0),
  fHistoTrkVtxX(0),
  fHistoTrkVtxY(0),
  fHistoTrkVtxZ(0),
  fHistoSPD3DVtxResidX(0),
  fHistoSPD3DVtxResidY(0),
  fHistoSPD3DVtxResidZ(0),
  fHistoSPDZVtxResidZ(0),
  fHistoTrkVtxResidX(0),
  fHistoTrkVtxResidY(0),
  fHistoTrkVtxResidZ(0),
  fHistoTracklets(0),
  fHistoSelTracks(0),
  fHistoGenMultVsb(0),
  fHistoTrackletsVsb(0),
  fHistoSelTracksVsb(0),
  fSpeciesAbundance(0),
  fIsAA(kFALSE),
  fNumOfSpeciesToCheck(0)
{
  //
  fNumOfSpeciesToCheck=0;
  fPdgCodes[fNumOfSpeciesToCheck++]=11;
  fPdgCodes[fNumOfSpeciesToCheck++]=-11;
  fPdgCodes[fNumOfSpeciesToCheck++]=211;
  fPdgCodes[fNumOfSpeciesToCheck++]=-211;
  fPdgCodes[fNumOfSpeciesToCheck++]=111;
  fPdgCodes[fNumOfSpeciesToCheck++]=321;
  fPdgCodes[fNumOfSpeciesToCheck++]=-321;
  fPdgCodes[fNumOfSpeciesToCheck++]=310;
  fPdgCodes[fNumOfSpeciesToCheck++]=130;
  fPdgCodes[fNumOfSpeciesToCheck++]=2212;
  fPdgCodes[fNumOfSpeciesToCheck++]=-2212;
  fPdgCodes[fNumOfSpeciesToCheck++]=3122;
  fPdgCodes[fNumOfSpeciesToCheck++]=-3122;
  fPdgCodes[fNumOfSpeciesToCheck++]=3312;
  fPdgCodes[fNumOfSpeciesToCheck++]=-3312;
  fPdgCodes[fNumOfSpeciesToCheck++]=3334;
  fPdgCodes[fNumOfSpeciesToCheck++]=-3334;
  fPdgCodes[fNumOfSpeciesToCheck++]=421;
  fPdgCodes[fNumOfSpeciesToCheck++]=-421;
  fPdgCodes[fNumOfSpeciesToCheck++]=511;
  fPdgCodes[fNumOfSpeciesToCheck++]=-511;
  fPdgCodes[fNumOfSpeciesToCheck++]=1000010020;
  fPdgCodes[fNumOfSpeciesToCheck++]=-1000010020;
  
  for(Int_t j=0; j<kMaxNumOfSpeciesToCheck; j++){
    fEtaPt[j]=0x0;
    fPrimSec[j]=0x0;
    fNumOfDau[j]=0x0;
    fDecLen[j]=0x0;
    fCt[j]=0x0;
    fMassDiff[j]=0x0;
    fMomDiff[j]=0x0;
    fPrimSecb[j]=0x0;
  }

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}


//___________________________________________________________________________
AliAnalysisTaskCheckGenKine::~AliAnalysisTaskCheckGenKine(){
  //
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
}
   
//___________________________________________________________________________
void AliAnalysisTaskCheckGenKine::AddParticleToCheck(Int_t pdg){
  // add particle PDG code in the list of species to be checked

  for(Int_t kp=0; kp<fNumOfSpeciesToCheck; kp++){
    if(fPdgCodes[kp]==pdg){
      printf("Particle %d already in list of species to be checked\n",pdg);
      return;
    }
  }
  if(fNumOfSpeciesToCheck<kMaxNumOfSpeciesToCheck) fPdgCodes[fNumOfSpeciesToCheck++]=pdg;
  else printf("Maximum number of particles species (%d) already reached, won't add %d\n",kMaxNumOfSpeciesToCheck,pdg);
  return;
}
//___________________________________________________________________________
void AliAnalysisTaskCheckGenKine::PrintSpeciesToCheck(){
  // print the list of particle species that will be checked
  printf("--- Number of particle species that will be checked = %d\n",fNumOfSpeciesToCheck);
  for(Int_t kp=0; kp<fNumOfSpeciesToCheck; kp++) printf(" %d)  PDGcode=%d\n",kp,fPdgCodes[kp]);
  return;
}
//___________________________________________________________________________
void AliAnalysisTaskCheckGenKine::UserCreateOutputObjects() {
  /// create output histos

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  Double_t minMult=-0.5;
  Double_t maxMult=99.5;
  if(fIsAA){
    maxMult=5000.;
    minMult=0.;
  }
  
  fHistoNEvents = new TH1F("hNEvents", "Number of processed events",4,-0.5,3.5);
  fHistoNEvents->SetMinimum(0);
  fHistoNEvents->GetXaxis()->SetBinLabel(1,"Analyzed");
  fHistoNEvents->GetXaxis()->SetBinLabel(2,"SPD vertex 3D");
  fHistoNEvents->GetXaxis()->SetBinLabel(3,"SPD vertex Z");
  fHistoNEvents->GetXaxis()->SetBinLabel(4,"Track vertex");
  fOutput->Add(fHistoNEvents);

  fHistoGenMult = new TH1F("hGenMult"," ;  N_{gen} (charged, |#eta|<0.9)",100,minMult,5.*maxMult);
  fOutput->Add(fHistoGenMult);

  // vertex histos
  fHistoVtxContrib = new TH3F("kVtxContrib"," ;  Nch (pi,K,p |#eta|<0.9) ; SPDVert contrib. ; TrackVert contrib.",100,-0.5,99.5,102,-2.5,99.5,102,-2.5,99.5);
  fHistoSPD3DVtxX = new TH1F("hSPD3DVtxX"," ; SPD 3Dvertex X (cm)",200,-1.,1.);
  fHistoSPD3DVtxY = new TH1F("hSPD3DVtxY"," ; SPD 3Dvertex Y (cm)",200,-1.,1.);
  fHistoSPD3DVtxZ = new TH1F("hSPD3DVtxZ"," ; SPD 3Dvertex Z (cm)",200,-20.,20.);
  fHistoSPDZVtxZ = new TH1F("hSPDZVtxZ"," ; SPD Zvertex Z (cm)",200,-20.,20.);
  fHistoTrkVtxX = new TH1F("hTrkVtxX"," ; Track vertex X (cm)",200,-1.,1.);
  fHistoTrkVtxY = new TH1F("hTrkVtxY"," ; Track vertex Y (cm)",200,-1.,1.);
  fHistoTrkVtxZ = new TH1F("hTrkVtxZ"," ; Track vertex Z (cm)",200,-20.,20.);
  fHistoSPD3DVtxResidX = new TH3F("hSPD3DVtxResidX"," ; SPDVert contrib. ; z_{vertex} (cm) ; x_{reco}-x_{gen} (cm)",100,minMult,maxMult,40,-20.,20.,100,-0.1,0.1);
  fHistoSPD3DVtxResidY = new TH3F("hSPD3DVtxResidY"," ; SPDVert contrib. ; z_{vertex} (cm) ; y_{reco}-y_{gen} (cm)",100,minMult,maxMult,40,-20.,20.,100,-0.1,0.1);
  fHistoSPD3DVtxResidZ = new TH3F("hSPD3DVtxResidZ"," ; SPDVert contrib. ; z_{vertex} (cm) ; z_{reco}-z_{gen} (cm)",100,minMult,maxMult,40,-20.,20.,100,-0.1,0.1);
  fHistoSPDZVtxResidZ = new TH3F("hSPDZVtxResidZ"," ; SPDVert contrib. ; z_{vertex} (cm) ; z_{reco}-z_{gen} (cm)",100,minMult,maxMult,40,-20.,20.,100,-0.1,0.1);
  fHistoTrkVtxResidX = new TH3F("hTrkVtxResidX"," ; SPDVert contrib. ; z_{vertex} (cm) ; x_{reco}-x_{gen} (cm)",100,minMult,maxMult,40,-20.,20.,100,-0.1,0.1);
  fHistoTrkVtxResidY = new TH3F("hTrkVtxResidY"," ; SPDVert contrib. ; z_{vertex} (cm) ; y_{reco}-y_{gen} (cm)",100,minMult,maxMult,40,-20.,20.,100,-0.1,0.1);
  fHistoTrkVtxResidZ = new TH3F("hTrkVtxResidZ"," ; SPDVert contrib. ; z_{vertex} (cm) ; z_{reco}-z_{gen} (cm)",100,minMult,maxMult,40,-20.,20.,100,-0.1,0.1);
  fOutput->Add(fHistoVtxContrib);
  fOutput->Add(fHistoSPD3DVtxX);  
  fOutput->Add(fHistoSPD3DVtxY);
  fOutput->Add(fHistoSPD3DVtxZ);
  fOutput->Add(fHistoSPDZVtxZ);
  fOutput->Add(fHistoTrkVtxX);
  fOutput->Add(fHistoTrkVtxY);
  fOutput->Add(fHistoTrkVtxZ);
  fOutput->Add(fHistoSPD3DVtxResidX);
  fOutput->Add(fHistoSPD3DVtxResidY);
  fOutput->Add(fHistoSPD3DVtxResidZ);
  fOutput->Add(fHistoSPDZVtxResidZ);
  fOutput->Add(fHistoTrkVtxResidX);
  fOutput->Add(fHistoTrkVtxResidY);
  fOutput->Add(fHistoTrkVtxResidZ);

  // multiplicity histos
  fHistoTracklets = new TH2F("hTracklets"," ; N_{gen} (phys prim, |#eta|<0.9) ; N_{tracklets}",100,minMult,maxMult,100,minMult,maxMult);
  fHistoSelTracks = new TH2F("hSelTracks"," ; N_{gen} (phys prim, |#eta|<0.9) ; N_{TPC+ITS tracks}",100,minMult,maxMult,100,minMult,maxMult);
  fOutput->Add(fHistoTracklets);
  fOutput->Add(fHistoSelTracks);
  if(fIsAA){
    fHistoGenMultVsb = new TH2F("hGenMultVsb"," ;  impact parameter (fm) ; N_{gen} (phys prim, |#eta|<0.9)",150,0.,15.,100,minMult,maxMult);
    fHistoTrackletsVsb = new TH2F("hTrackletsVsb"," ; impact parameter (fm) ; N_{tracklets}",150,0.,15.,100,minMult,maxMult);
    fHistoSelTracksVsb = new TH2F("hSelTracksVsb"," ; impact parameter (fm) ; N_{TPC+ITS tracks}",150,0.,15.,100,minMult,maxMult);
    fOutput->Add(fHistoGenMultVsb);
    fOutput->Add(fHistoTrackletsVsb);
    fOutput->Add(fHistoSelTracksVsb);
  }
  
  // per-particle histos
  fSpeciesAbundance = new TH2F("hSpeciesAbundance","",fNumOfSpeciesToCheck,-0.5,fNumOfSpeciesToCheck-0.5,2,-0.5,1.5);
  fSpeciesAbundance->GetYaxis()->SetBinLabel(1,"From generator");
  fSpeciesAbundance->GetYaxis()->SetBinLabel(2,"From transport");
  fOutput->Add(fSpeciesAbundance);
  for(Int_t j=0; j<fNumOfSpeciesToCheck; j++){
    TString pname=TDatabasePDG::Instance()->GetParticle(fPdgCodes[j])->GetName();
    fSpeciesAbundance->GetXaxis()->SetBinLabel(j+1,pname.Data());
    fEtaPt[j] = new TH2F(TString::Format("hEtaPt%s",pname.Data())," ; #eta ; p_{T} (GeV/c)",20,-10,10,100,0.,20.);
    fPrimSec[j] = new TH3F(TString::Format("hPrimSec%s",pname.Data())," ; ; p_{T} (GeV/c) ; dist from vert (cm)",4,-0.5,3.5,100,0.,20.,100,0.,100.);
    fPrimSec[j]->GetXaxis()->SetBinLabel(1,"Primary");
    fPrimSec[j]->GetXaxis()->SetBinLabel(2,"Secondary from weak");
    fPrimSec[j]->GetXaxis()->SetBinLabel(3,"Secondary from material");
    fPrimSec[j]->GetXaxis()->SetBinLabel(4,"Other");
    fNumOfDau[j] = new TH2F(TString::Format("hNumOfDau%s",pname.Data())," ; Number of daughters",11,-0.5,10.5,6,-0.5,5.5);
    fNumOfDau[j]->GetYaxis()->SetBinLabel(1,"From generator, decay generator");
    fNumOfDau[j]->GetYaxis()->SetBinLabel(2,"From generator, decay transport");
    fNumOfDau[j]->GetYaxis()->SetBinLabel(3,"From generator, interaction transport");
    fNumOfDau[j]->GetYaxis()->SetBinLabel(4,"From transport, decay transport");
    fNumOfDau[j]->GetYaxis()->SetBinLabel(5,"From transport, interaction transport");
    fNumOfDau[j]->GetYaxis()->SetBinLabel(6,"Ohter");
    Double_t maxDL=20.;
    if((fPdgCodes[j]>400 && fPdgCodes[j]<600) || (fPdgCodes[j]>4000 && fPdgCodes[j]<6000)) maxDL=2.;
    fDecLen[j] = new TH2F(TString::Format("hDecLen%s",pname.Data())," ; p (GeV/c) ; decay length (cm)",100,0.,20.,200,0.,maxDL);
    fCt[j] = new TH2F(TString::Format("hCt%s",pname.Data())," ; p (GeV/c) ; ct (cm)",100,0.,20.,200,0.,maxDL);
    fMassDiff[j] = new TH2F(TString::Format("hMassDiff%s",pname.Data())," ; (M_{mother} - M_{daughters})/ M_{mother}",101,-0.0505,0.0505,3,-0.5,2.5);
    fMomDiff[j] = new TH2F(TString::Format("hMomDiff%s",pname.Data())," ; (p_{mother} - p_{daughters})/ p_{mother}",101,-0.0505,0.0505,3,-0.5,2.5);

    fMassDiff[j]->GetYaxis()->SetBinLabel(1,"Decay generator");
    fMassDiff[j]->GetYaxis()->SetBinLabel(2,"Decay transport");
    fMassDiff[j]->GetYaxis()->SetBinLabel(3,"Interaction transport");
    fMomDiff[j]->GetYaxis()->SetBinLabel(1,"Decay generator");
    fMomDiff[j]->GetYaxis()->SetBinLabel(2,"Decay transport");
    fMomDiff[j]->GetYaxis()->SetBinLabel(3,"Interaction transport");

    fOutput->Add(fEtaPt[j]);
    fOutput->Add(fPrimSec[j]);
    fOutput->Add(fNumOfDau[j]);
    fOutput->Add(fDecLen[j]);
    fOutput->Add(fCt[j]);
    fOutput->Add(fMassDiff[j]);
    fOutput->Add(fMomDiff[j]);
    if(fIsAA){
      fPrimSecb[j] = new TH3F(TString::Format("hPrimSecb%s",pname.Data())," ; ; p_{T} (GeV/c) ; impact parameter (fm)",4,-0.5,3.5,100,0.,20.,150,0.,15.);
      fPrimSecb[j]->GetXaxis()->SetBinLabel(1,"Primary");
      fPrimSecb[j]->GetXaxis()->SetBinLabel(2,"Secondary from weak");
      fPrimSecb[j]->GetXaxis()->SetBinLabel(3,"Secondary from material");
      fPrimSecb[j]->GetXaxis()->SetBinLabel(4,"Other");
      fOutput->Add(fPrimSecb[j]);
    }
  }

  PostData(1,fOutput);

}
//______________________________________________________________________________
void AliAnalysisTaskCheckGenKine::UserExec(Option_t *)
{
  //

  AliESDEvent *esd = (AliESDEvent*) (InputEvent());

  if(!esd) {
    printf("AliAnalysisTaskCheckGenKine::Exec(): bad ESD\n");
    return;
  } 

  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler) {
    printf("AliAnalysisTaskCheckGenKine::Exec(): could not retrieve MC event handler");
    return;
  }
  AliMCEvent* mcEvent = eventHandler->MCEvent();
  if (!mcEvent) {
    Printf("AliAnalysisTaskCheckGenKine::Exec(): could not retrieve MC event");
    return;
  }
  const AliVVertex* mcVert=mcEvent->GetPrimaryVertex();
  if(!mcVert){
    Printf("AliAnalysisTaskCheckGenKine::Exec(): generated vertex not available");
    return;
  }
  Double_t imppar=-999.;
  TString genname=mcEvent->GenEventHeader()->ClassName();
  if(genname.Contains("CocktailEventHeader")){
    AliGenCocktailEventHeader *cockhead=(AliGenCocktailEventHeader*)mcEvent->GenEventHeader();
    TList* lgen=cockhead->GetHeaders();
    for(Int_t ig=0; ig<lgen->GetEntries(); ig++){
      AliGenerator* gen=(AliGenerator*)lgen->At(ig);
      TString title=gen->GetName();
      if(title.Contains("hijing") || title.Contains("Hijing")){
	AliGenHijingEventHeader* hijh=(AliGenHijingEventHeader*)lgen->At(ig);
	imppar=hijh->ImpactParameter();
      }
    }
  }else if(genname.Contains("hijing") || genname.Contains("Hijing")){
    AliGenHijingEventHeader* hijh=(AliGenHijingEventHeader*)mcEvent->GenEventHeader();
    imppar=hijh->ImpactParameter();
  }
  fHistoNEvents->Fill(0);

  //vertices
  const AliESDVertex *spdv=esd->GetPrimaryVertexSPD();
  const AliESDVertex *trkv=esd->GetPrimaryVertexTracks();
  Int_t nContribSPD=-3;
  if(spdv) nContribSPD=spdv->GetNContributors();
  Int_t nContribTRK=-1;
  if(trkv) nContribTRK=trkv->GetNContributors();
  
  // generated multiplicities
  Int_t nParticles=mcEvent->GetNumberOfTracks();
  Int_t nChEta09 = 0.;
  Int_t nPhysPrimEta09=0;
  Int_t nChPhysPrimEta09=0;
  Int_t nPiKPEta09=0;
  for (Int_t i=0;i<nParticles;i++){
    TParticle* part = (TParticle*)mcEvent->Particle(i);
    TParticlePDG* ppdg = part->GetPDG();
    Int_t absPdg=TMath::Abs(part->GetPdgCode());
    Double_t eta=part->Eta();
    if(TMath::Abs(eta)<0.9){
      if(mcEvent->IsPhysicalPrimary(i)){
	nPhysPrimEta09++;
	if(ppdg && TMath::Abs(ppdg->Charge())>0.) nChPhysPrimEta09++;
      }
      if(absPdg==211 || absPdg==321 || absPdg==2212) nPiKPEta09++;
      if(ppdg && TMath::Abs(ppdg->Charge())>0.) nChEta09++;
    }
  }
  fHistoGenMult->Fill(nChEta09);
  fHistoVtxContrib->Fill(nPiKPEta09,nContribSPD,nContribTRK);

  // primary vertex plots
  if(nContribSPD>=1){
    Double_t dx=spdv->GetX()-mcVert->GetX();
    Double_t dy=spdv->GetY()-mcVert->GetY();
    Double_t dz=spdv->GetZ()-mcVert->GetZ();
    if(spdv->IsFromVertexer3D()){
      fHistoNEvents->Fill(1);
      fHistoSPD3DVtxX->Fill(spdv->GetX());
      fHistoSPD3DVtxY->Fill(spdv->GetY());
      fHistoSPD3DVtxZ->Fill(spdv->GetZ());
      fHistoSPD3DVtxResidX->Fill(nContribSPD,mcVert->GetZ(),dx);
      fHistoSPD3DVtxResidY->Fill(nContribSPD,mcVert->GetZ(),dy);
      fHistoSPD3DVtxResidZ->Fill(nContribSPD,mcVert->GetZ(),dz);
    }else if(spdv->IsFromVertexerZ()){
      fHistoNEvents->Fill(2);
      fHistoSPDZVtxZ->Fill(spdv->GetZ());
      fHistoSPDZVtxResidZ->Fill(nContribSPD,mcVert->GetZ(),dz);
    }
  }
  if(nContribTRK>=1){
    Double_t dx=trkv->GetX()-mcVert->GetX();
    Double_t dy=trkv->GetY()-mcVert->GetY();
    Double_t dz=trkv->GetZ()-mcVert->GetZ();
    fHistoNEvents->Fill(3);
    fHistoTrkVtxX->Fill(trkv->GetX());
    fHistoTrkVtxY->Fill(trkv->GetY());
    fHistoTrkVtxZ->Fill(trkv->GetZ());
    fHistoTrkVtxResidX->Fill(nContribTRK,mcVert->GetZ(),dx);
    fHistoTrkVtxResidY->Fill(nContribTRK,mcVert->GetZ(),dy);
    fHistoTrkVtxResidZ->Fill(nContribTRK,mcVert->GetZ(),dz);
   }



  const AliMultiplicity* mult=esd->GetMultiplicity();
  Int_t nTracklets=mult->GetNumberOfTracklets();
  Int_t nTrackletsEta09=0;
  for(Int_t it=0; it<nTracklets; it++){
    Double_t eta=TMath::Abs(mult->GetEta(it));
    if(eta<0.9) nTrackletsEta09++;
  }
  fHistoTracklets->Fill(nChPhysPrimEta09,nTrackletsEta09);

  Int_t nTracks=esd->GetNumberOfTracks();
  Int_t nSelTracks=0;
  for(Int_t it=0; it<nTracks; it++){
    AliESDtrack* tr=esd->GetTrack(it);
    UInt_t status=tr->GetStatus();
    if(!(status&AliESDtrack::kITSrefit)) continue;
    if(!(status&AliESDtrack::kTPCin)) continue;
    nSelTracks++;
  }
  fHistoSelTracks->Fill(nChPhysPrimEta09,nSelTracks);

  if(fIsAA){
    fHistoGenMultVsb->Fill(imppar,nChPhysPrimEta09);
    fHistoTrackletsVsb->Fill(imppar,nTrackletsEta09);
    fHistoSelTracksVsb->Fill(imppar,nSelTracks);
  }
  
  for (Int_t i=0;i<nParticles;i++){
    AliMCParticle* mcPart=(AliMCParticle*)mcEvent->GetTrack(i);
    TParticle* part = (TParticle*)mcEvent->Particle(i);
    if(!mcPart || !part) continue;
    Int_t pdg=part->GetPdgCode();
    Int_t spId=GetSpeciesIndex(pdg);
    if(spId<0) continue;
    Double_t pt=part->Pt();
    Double_t mom=part->P();
    Double_t mass=part->GetMass();
    Double_t eta=part->Eta();
    Int_t fromGener=1;
    if(i<mcEvent->GetNumberOfPrimaries()) fromGener=0;
    fSpeciesAbundance->Fill(spId,fromGener);
    if(fromGener==0) fEtaPt[spId]->Fill(eta,pt);
    if(TMath::Abs(eta)>0.8) continue;
    Double_t distx=part->Vx()-mcVert->GetX();
    Double_t disty=part->Vy()-mcVert->GetY();
    Double_t distz=part->Vz()-mcVert->GetZ();
    Double_t distToVert=TMath::Sqrt(distx*distx+disty*disty+distz*distz);
    Int_t primSec=3;
    if(mcEvent->IsPhysicalPrimary(i)) primSec=0;
    else{
      if(mcEvent->IsSecondaryFromWeakDecay(i)) primSec=1;
      else if(mcEvent->IsSecondaryFromMaterial(i)) primSec=2;
    }
    fPrimSec[spId]->Fill(primSec,pt,distToVert);
    if(fIsAA) fPrimSecb[spId]->Fill(primSec,pt,imppar);
    Int_t nDau=mcPart->GetNDaughters();
    Int_t iDau=mcPart->GetDaughterFirst();
    if(iDau>=0){

      TParticle* firstDau = (TParticle*)mcEvent->Particle(iDau);
      Int_t yVal=5;
      Int_t dauOrig=-1;
      if(iDau<mcEvent->GetNumberOfPrimaries()){
	// produced by generator, decayed by generator
	dauOrig=0;
      }
      Int_t process=firstDau->GetUniqueID();
      if(process==4) dauOrig=1; //decay
      else if(process>=13) dauOrig=2; //hadronic interaction
      if(dauOrig>=0){
	if(fromGener==0) yVal=dauOrig;
	else if(fromGener==1 && dauOrig>=1) yVal=dauOrig+2;
      } 
      fNumOfDau[spId]->Fill(nDau,yVal);
      Double_t dxDau=firstDau->Vx()-part->Vx();
      Double_t dyDau=firstDau->Vy()-part->Vy();
      Double_t dzDau=firstDau->Vz()-part->Vz();
      Double_t decLen=TMath::Sqrt(dxDau*dxDau+dyDau*dyDau+dzDau*dzDau);
      Double_t sumPxDau=0.;
      Double_t sumPyDau=0.;
      Double_t sumPzDau=0.;
      Double_t sumEDau=0.;
      for(Int_t j=0; j<nDau; j++){
	TParticle* partDau = (TParticle*)mcEvent->Particle(iDau+j);
	if(partDau){
	  sumPxDau+=partDau->Px();
	  sumPyDau+=partDau->Py();
	  sumPzDau+=partDau->Pz();
	  sumEDau+=partDau->Energy();
	}
      }
      Double_t pSquareDau=sumPxDau*sumPxDau+sumPyDau*sumPyDau+sumPzDau*sumPzDau;
      Double_t pDau=TMath::Sqrt(pSquareDau);
      Double_t pDiff=(mom-pDau)/mom;
      fMomDiff[spId]->Fill(pDiff,dauOrig);
      Double_t invMass2=sumEDau*sumEDau-pSquareDau;
      if(invMass2>=0){
	Double_t invMass=TMath::Sqrt(sumEDau*sumEDau-pSquareDau);
	Double_t mDiff=(mass-invMass)/mass;
	fMassDiff[spId]->Fill(mDiff,dauOrig);
      }
      fDecLen[spId]->Fill(mom,decLen);
      fCt[spId]->Fill(mom,decLen*mass/mom);
    }
  }

  PostData(1,fOutput);
  
}

//______________________________________________________________________________
void AliAnalysisTaskCheckGenKine::Terminate(Option_t */*option*/)
{
  /// Terminate analysis
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }

  return;
}




