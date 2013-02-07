#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliESDEvent.h"
#include "AliStack.h"
#include "AliCentrality.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMultiplicity.h"
#include <TParticle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TChain.h>
#include "AliESDInputHandlerRP.h"
#include "AliAnalysisTaskCheckHFMCProd.h"

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
// Implementation of class AliAnalysisTaskCheckHFMCProd
// AliAnalysisTask to check MC production at ESD+Kine level
// 
//
// Authors: F. Prino, prino@to.infn.it
//          
//*************************************************************************

ClassImp(AliAnalysisTaskCheckHFMCProd)
//______________________________________________________________________________
AliAnalysisTaskCheckHFMCProd::AliAnalysisTaskCheckHFMCProd() : AliAnalysisTaskSE("HFMCChecks"), 
  fOutput(0),
  fHistoNEvents(0),
  fHistoPhysPrim(0),
  fHistoTracks(0),
  fHistoSelTracks(0),
  fHistoTracklets(0),
  fHistoSPD3DVtxX(0),
  fHistoSPD3DVtxY(0),
  fHistoSPD3DVtxZ(0),
  fHistoSPDZVtxX(0),
  fHistoSPDZVtxY(0),
  fHistoSPDZVtxZ(0),
  fHistoTRKVtxX(0),
  fHistoTRKVtxY(0),
  fHistoTRKVtxZ(0),
  fHistoNcharmed(0),
  fHistoNbVsNc(0),
  fSystem(0),
  fReadMC(kTRUE)
{
  //
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}


//___________________________________________________________________________
AliAnalysisTaskCheckHFMCProd::~AliAnalysisTaskCheckHFMCProd(){
  //
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
}
   
//___________________________________________________________________________
void AliAnalysisTaskCheckHFMCProd::UserCreateOutputObjects() {
  // create output histos

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  fHistoNEvents = new TH1F("hNEvents", "Number of processed events",3,-0.5,2.5);
  fHistoNEvents->Sumw2();
  fHistoNEvents->SetMinimum(0);
  fOutput->Add(fHistoNEvents);

  Double_t maxMult=100.;
  if(fSystem==1) maxMult=10000.;
  if(fSystem==2) maxMult=500.;
  fHistoPhysPrim = new TH1F("hPhysPrim","",100,0.,maxMult);
  fHistoPhysPrim->Sumw2();
  fOutput->Add(fHistoPhysPrim);
  fHistoTracks = new TH1F("hTracks","",100,0.,maxMult*2);
  fHistoTracks->Sumw2();
  fOutput->Add(fHistoTracks);
  fHistoSelTracks = new TH1F("hSelTracks","",100,0.,maxMult);
  fHistoSelTracks->Sumw2();
  fOutput->Add(fHistoSelTracks);
  fHistoTracklets = new TH1F("hTracklets","",100,0.,maxMult);
  fHistoTracklets->Sumw2();
  fOutput->Add(fHistoTracklets);

  fHistoSPD3DVtxX = new TH1F("hSPD3DvX","",100,-1.,1.);
  fHistoSPD3DVtxX->Sumw2();
  fOutput->Add(fHistoSPD3DVtxX);
  fHistoSPD3DVtxY = new TH1F("hSPD3DvY","",100,-1.,1.);
  fHistoSPD3DVtxY->Sumw2();
  fOutput->Add(fHistoSPD3DVtxY);
  fHistoSPD3DVtxZ = new TH1F("hSPD3DvZ","",100,-15.,15.);
  fHistoSPD3DVtxZ->Sumw2();
  fOutput->Add(fHistoSPD3DVtxZ);

  fHistoSPDZVtxX = new TH1F("hSPDZvX","",100,-1.,1.);
  fHistoSPDZVtxX->Sumw2();
  fOutput->Add(fHistoSPDZVtxX);
  fHistoSPDZVtxY = new TH1F("hSPDZvY","",100,-1.,1.);
  fHistoSPDZVtxY->Sumw2();
  fOutput->Add(fHistoSPDZVtxY);
  fHistoSPDZVtxZ = new TH1F("hSPDZvZ","",100,-15.,15.);
  fHistoSPDZVtxZ->Sumw2();
  fOutput->Add(fHistoSPDZVtxZ);


  fHistoTRKVtxX = new TH1F("hTRKvX","",100,-1.,1.);
  fHistoTRKVtxX->Sumw2();
  fOutput->Add(fHistoTRKVtxX);
  fHistoTRKVtxY = new TH1F("hTRKvY","",100,-1.,1.);
  fHistoTRKVtxY->Sumw2();
  fOutput->Add(fHistoTRKVtxY);
  fHistoTRKVtxZ = new TH1F("hTRKvZ","",100,-15.,15.);
  fHistoTRKVtxZ->Sumw2();
  fOutput->Add(fHistoTRKVtxZ);

  Int_t nBinscb=11;
  if(fSystem==1) nBinscb=200;
  if(fSystem==2) nBinscb=21;
  Double_t maxncn=nBinscb-0.5;
  fHistoNcharmed = new TH2F("hncharmed","",100,0.,maxMult,nBinscb,-0.5,maxncn);
  fHistoNcharmed->Sumw2();
  fOutput->Add(fHistoNcharmed);
  fHistoNbVsNc = new TH2F("hnbvsnc","",nBinscb,-0.5,maxncn,nBinscb,-0.5,maxncn);
  fHistoNbVsNc->Sumw2();
  fOutput->Add(fHistoNbVsNc);

  fHistYPtPrompt[0] = new TH2F("hyptD0prompt","D0 - Prompt",40,0.,40.,20,-2.,2.);
  fHistYPtPrompt[1] = new TH2F("hyptDplusprompt","Dplus - Prompt",40,0.,40.,20,-2.,2.);
  fHistYPtPrompt[2] = new TH2F("hyptDstarprompt","Dstar - Prompt",40,0.,40.,20,-2.,2.);
  fHistYPtPrompt[3] = new TH2F("hyptDsprompt","Ds - Prompt",40,0.,40.,20,-2.,2.);
  fHistYPtPrompt[4] = new TH2F("hyptLcprompt","Lc - Prompt",40,0.,40.,20,-2.,2.);

  fHistBYPtAllDecay[0] = new TH2F("hyptB0AllDecay","B0 - All",40,0.,40.,40,-2.,2.);
  fHistBYPtAllDecay[1] = new TH2F("hyptBplusAllDecay","Bplus - All",40,0.,40.,40,-2.,2.);
  fHistBYPtAllDecay[2] = new TH2F("hyptBstarAllDecay","Bstar - All",40,0.,40.,40,-2.,2.);
  fHistBYPtAllDecay[3] = new TH2F("hyptBsAllDecay","Bs - All",40,0.,40.,40,-2.,2.);
  fHistBYPtAllDecay[4] = new TH2F("hyptLbAllDecay","LB - All",40,0.,40.,40,-2.,2.);

  fHistYPtAllDecay[0] = new TH2F("hyptD0AllDecay","D0 - All",40,0.,40.,40,-2.,2.);
  fHistYPtAllDecay[1] = new TH2F("hyptDplusAllDecay","Dplus - All",40,0.,40.,40,-2.,2.);
  fHistYPtAllDecay[2] = new TH2F("hyptDstarAllDecay","Dstar - All",40,0.,40.,40,-2.,2.);
  fHistYPtAllDecay[3] = new TH2F("hyptDsAllDecay","Ds - All",40,0.,40.,40,-2.,2.);
  fHistYPtAllDecay[4] = new TH2F("hyptLcAllDecay","Lc - All",40,0.,40.,40,-2.,2.);

  fHistYPtPromptAllDecay[0] = new TH2F("hyptD0promptAllDecay","D0 - Prompt",40,0.,40.,40,-2.,2.);
  fHistYPtPromptAllDecay[1] = new TH2F("hyptDpluspromptAllDecay","Dplus - Prompt",40,0.,40.,40,-2.,2.);
  fHistYPtPromptAllDecay[2] = new TH2F("hyptDstarpromptAllDecay","Dstar - Prompt",40,0.,40.,40,-2.,2.);
  fHistYPtPromptAllDecay[3] = new TH2F("hyptDspromptAllDecay","Ds - Prompt",40,0.,40.,40,-2.,2.);
  fHistYPtPromptAllDecay[4] = new TH2F("hyptLcpromptAllDecay","Lc - Prompt",40,0.,40.,40,-2.,2.);

  fHistYPtFeeddownAllDecay[0] = new TH2F("hyptD0feeddownAllDecay","D0 - FromB",40,0.,40.,40,-2.,2.);
  fHistYPtFeeddownAllDecay[1] = new TH2F("hyptDplusfeeddownAllDecay","Dplus - FromB",40,0.,40.,40,-2.,2.);
  fHistYPtFeeddownAllDecay[2] = new TH2F("hyptDstarfeeddownAllDecay","Dstar - FromB",40,0.,40.,40,-2.,2.);
  fHistYPtFeeddownAllDecay[3] = new TH2F("hyptDsfeeddownAllDecay","Ds - FromB",40,0.,40.,40,-2.,2.);
  fHistYPtFeeddownAllDecay[4] = new TH2F("hyptLcfeeddownAllDecay","Lc - FromB",40,0.,40.,40,-2.,2.);


  fHistYPtFeeddown[0] = new TH2F("hyptD0feeddown","D0 - Feeddown",40,0.,40.,20,-2.,2.);
  fHistYPtFeeddown[1] = new TH2F("hyptDplusfeeddown","Dplus - Feeddown",40,0.,40.,20,-2.,2.);
  fHistYPtFeeddown[2] = new TH2F("hyptDstarfeedown","Dstar - Feeddown",40,0.,40.,20,-2.,2.);
  fHistYPtFeeddown[3] = new TH2F("hyptDsfeedown","Ds - Feeddown",40,0.,40.,20,-2.,2.);
  fHistYPtFeeddown[4] = new TH2F("hyptLcfeedown","Lc - Feeddown",40,0.,40.,20,-2.,2.);

  for(Int_t ih=0; ih<5; ih++){
    fHistBYPtAllDecay[ih]->Sumw2();
    fHistBYPtAllDecay[ih]->SetMinimum(0);
    fOutput->Add(fHistBYPtAllDecay[ih]);
    fHistYPtAllDecay[ih]->Sumw2();
    fHistYPtAllDecay[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtAllDecay[ih]);
    fHistYPtPromptAllDecay[ih]->Sumw2();
    fHistYPtPromptAllDecay[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtPromptAllDecay[ih]);
    fHistYPtFeeddownAllDecay[ih]->Sumw2();
    fHistYPtFeeddownAllDecay[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtFeeddownAllDecay[ih]);
    fHistYPtPrompt[ih]->Sumw2();
    fHistYPtPrompt[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtPrompt[ih]);
    fHistYPtFeeddown[ih]->Sumw2();
    fHistYPtFeeddown[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtFeeddown[ih]);
  }

  fHistYPtD0byDecChannel[0] = new TH2F("hyptD02","D0 - 2prong",40,0.,40.,20,-2.,2.);
  fHistYPtD0byDecChannel[1] = new TH2F("hyptD04","D0 - 4prong",40,0.,40.,20,-2.,2.);
  fHistYPtDplusbyDecChannel[0] = new TH2F("hyptDplusnonreson","Dplus - non reson",40,0.,40.,20,-2.,2.);
  fHistYPtDplusbyDecChannel[1] = new TH2F("hyptDplusreson","Dplus - reson via K0*",40,0.,40.,20,-2.,2.);
  fHistYPtDsbyDecChannel[0] = new TH2F("hyptDsphi","Ds - vis Phi",40,0.,40.,20,-2.,2.);
  fHistYPtDsbyDecChannel[1] = new TH2F("hyptDsk0st","Ds - via k0*",40,0.,40.,20,-2.,2.);

  for(Int_t ih=0; ih<2; ih++){

    fHistYPtD0byDecChannel[ih]->Sumw2();
    fHistYPtD0byDecChannel[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtD0byDecChannel[ih]);
    fHistYPtDplusbyDecChannel[ih]->Sumw2();
    fHistYPtDplusbyDecChannel[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtDplusbyDecChannel[ih]);
    fHistYPtDsbyDecChannel[ih]->Sumw2();
    fHistYPtDsbyDecChannel[ih]->SetMinimum(0);
    fOutput->Add(fHistYPtDsbyDecChannel[ih]);
  }
    

  PostData(1,fOutput);

}
//______________________________________________________________________________
void AliAnalysisTaskCheckHFMCProd::UserExec(Option_t *)
{
  //

  AliESDEvent *esd = (AliESDEvent*) (InputEvent());


  if(!esd) {
    printf("AliAnalysisTaskSDDRP::Exec(): bad ESD\n");
    return;
  } 

  fHistoNEvents->Fill(0);

  Int_t nTracks=esd->GetNumberOfTracks();
  fHistoTracks->Fill(nTracks);
  Int_t nSelTracks=0;
  for(Int_t it=0; it<nTracks; it++){
    AliESDtrack* tr=esd->GetTrack(it);
    UInt_t status=tr->GetStatus();
    if(!(status&AliESDtrack::kITSrefit)) continue;
    if(!(status&AliESDtrack::kTPCin)) continue;
    nSelTracks++;
  }
  fHistoSelTracks->Fill(nSelTracks);

  const AliMultiplicity* mult=esd->GetMultiplicity();
  Int_t nTracklets=mult->GetNumberOfTracklets();
  fHistoTracklets->Fill(nTracklets);

  const AliESDVertex *spdv=esd->GetVertex();
  if(spdv && spdv->IsFromVertexer3D()){
    fHistoSPD3DVtxX->Fill(spdv->GetXv());
    fHistoSPD3DVtxY->Fill(spdv->GetYv());
    fHistoSPD3DVtxZ->Fill(spdv->GetZv());
  }
  if(spdv && spdv->IsFromVertexerZ()){
    fHistoSPDZVtxX->Fill(spdv->GetXv());
    fHistoSPDZVtxY->Fill(spdv->GetYv());
    fHistoSPDZVtxZ->Fill(spdv->GetZv());
  }
  const AliESDVertex *trkv=esd->GetPrimaryVertex();
  if(trkv && trkv->GetNContributors()>1){
    fHistoTRKVtxX->Fill(trkv->GetXv());
    fHistoTRKVtxY->Fill(trkv->GetYv());
    fHistoTRKVtxZ->Fill(trkv->GetZv());
  }

  AliStack* stack=0;

  if(fReadMC){
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }
    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }
    stack = mcEvent->Stack();
    if (!stack) {
      Printf("ERROR: stack not available");
      return;
    }
  

    Int_t nParticles=stack->GetNtrack();
    Double_t dNchdy = 0.;
    Int_t nb = 0, nc=0;
    Int_t nCharmed=0;
    Int_t nPhysPrim=0;
    for (Int_t i=0;i<nParticles;i++){
      TParticle* part = (TParticle*)stack->Particle(i);
      Int_t absPdg=TMath::Abs(part->GetPdgCode());
      if(absPdg==4) nc++;
      if(absPdg==5) nb++;
      if(stack->IsPhysicalPrimary(i)){
	Double_t eta=part->Eta();
	if(TMath::Abs(eta)<0.5){
	  dNchdy+=0.6666;   // 2/3 for the ratio charged/all
	  nPhysPrim++;
	}
      }
      Float_t rapid=-999.;
      if (part->Energy() != TMath::Abs(part->Pz())){
	rapid=0.5*TMath::Log((part->Energy()+part->Pz())/(part->Energy()-part->Pz()));
      }

       
      Int_t iPart=-1;
      Int_t iType=0;
      Int_t iSpecies=-1;
      if(absPdg==421){
	iSpecies=0;
	iType=CheckD0Decay(i,stack);
	if(iType>=0) iPart=0;	
      }
      else if(absPdg==411){
	iSpecies=1;
	iType=CheckDplusDecay(i,stack);
	if(iType>=0) iPart=1;
      }
      else if(absPdg==413){
	iSpecies=2;
	iType=CheckDstarDecay(i,stack);
	if(iType>=0) iPart=2;
      }
      else if(absPdg==431){
	iSpecies=3;
	iType=CheckDsDecay(i,stack);
	if(iType==0 || iType==1) iPart=3;
      }
      else if(absPdg==4122){
	iSpecies=4;
	iType=CheckLcDecay(i,stack);
	if(iType>=0) iPart=4;
      }
      if(iSpecies>=0) fHistYPtAllDecay[iSpecies]->Fill(part->Pt(),rapid);

      // check beauty mesons
      if(absPdg==511) fHistBYPtAllDecay[0]->Fill(part->Pt(),rapid);
      else if(absPdg==521) fHistBYPtAllDecay[1]->Fill(part->Pt(),rapid);
      else if(absPdg==513) fHistBYPtAllDecay[2]->Fill(part->Pt(),rapid);
      else if(absPdg==531) fHistBYPtAllDecay[3]->Fill(part->Pt(),rapid);
      else if(absPdg==5122) fHistBYPtAllDecay[4]->Fill(part->Pt(),rapid);

      if(iSpecies<0) continue; // not a charm meson

      TParticle* runningpart=part;
      Int_t iFromB=-1;
      Int_t pdgmoth=-1;
      while(1){
	Int_t labmoth=runningpart->GetFirstMother();
	if(labmoth==-1) break;
	TParticle *mot=(TParticle*)stack->Particle(labmoth);
	pdgmoth=TMath::Abs(mot->GetPdgCode());
	if(pdgmoth==5){ 
	  iFromB=1;
	  break;
	}else if(pdgmoth==4){
	  iFromB=0;
	  break;
	}
	runningpart=mot;
      }
      if(iFromB<0) continue;
      else if(iFromB==0) fHistYPtPromptAllDecay[iSpecies]->Fill(part->Pt(),rapid);
      else if(iFromB==1) fHistYPtFeeddownAllDecay[iSpecies]->Fill(part->Pt(),rapid);

      if(iPart<0) continue;
      if(iType<0) continue;
      nCharmed++;
      if(iPart==0 && iType<=1){
	fHistYPtD0byDecChannel[iType]->Fill(part->Pt(),rapid);
      }else if(iPart==1 && iType<=1){
	fHistYPtDplusbyDecChannel[iType]->Fill(part->Pt(),rapid);
      }else if(iPart==3 && iType<=1){
	fHistYPtDsbyDecChannel[iType]->Fill(part->Pt(),rapid);
      }
      
      if(iFromB==0 && iPart>=0 && iPart<5) fHistYPtPrompt[iPart]->Fill(part->Pt(),rapid);
      else if(iFromB==1 && iPart>=0 && iPart<5) fHistYPtFeeddown[iPart]->Fill(part->Pt(),rapid);
      
    }
    fHistoNcharmed->Fill(dNchdy,nCharmed);
    fHistoNbVsNc->Fill(nc,nb);
    fHistoPhysPrim->Fill(nPhysPrim);
  }

  PostData(1,fOutput);
  
}
//______________________________________________________________________________
void AliAnalysisTaskCheckHFMCProd::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }

  return;
}




//______________________________________________________________________________
Int_t AliAnalysisTaskCheckHFMCProd::CheckD0Decay(Int_t labD0, AliStack* stack) const{

  if(labD0<0) return -1;
  TParticle* dp = (TParticle*)stack->Particle(labD0);
  Int_t pdgdp=dp->GetPdgCode();
  Int_t nDau=dp->GetNDaughters();

  if(nDau==2){
    Int_t nKaons=0;
    Int_t nPions=0;
    for(Int_t iDau=dp->GetFirstDaughter(); iDau<=dp->GetLastDaughter(); iDau++){
      if(iDau<0) return -1;
      TParticle* dau=(TParticle*)stack->Particle(iDau);
      Int_t pdgdau=dau->GetPdgCode();
      if(TMath::Abs(pdgdau)==321){
	if(pdgdp>0 && pdgdau>0) return -1;
	if(pdgdp<0 && pdgdau<0) return -1;
	nKaons++;
      }else if(TMath::Abs(pdgdau)==211){
	if(pdgdp<0 && pdgdau>0) return -1;
	if(pdgdp>0 && pdgdau<0) return -1;
	nPions++;
      }else{
	return -1;
      }
    }
    if(nPions!=1) return -1;
    if(nKaons!=1) return -1;
    for(Int_t iDau=dp->GetFirstDaughter(); iDau<=dp->GetLastDaughter(); iDau++){
      if(iDau<0) return -1;
    }
    return 0;
  }

  if(nDau==3 || nDau==4){
    Int_t nKaons=0;
    Int_t nPions=0;
    for(Int_t iDau=dp->GetFirstDaughter(); iDau<=dp->GetLastDaughter(); iDau++){
      if(iDau<0) return -1;
      TParticle* dau=(TParticle*)stack->Particle(iDau);
      Int_t pdgdau=dau->GetPdgCode();
      if(TMath::Abs(pdgdau)==321){
	if(pdgdp>0 && pdgdau>0) return -1;
	if(pdgdp<0 && pdgdau<0) return -1;
	nKaons++;
      }else if(TMath::Abs(pdgdau)==113 || TMath::Abs(pdgdau)==313){
	for(Int_t resDau=dau->GetFirstDaughter(); resDau<=dau->GetLastDaughter(); resDau++){
	  if(resDau<0) return -1;
	  TParticle* resdau=(TParticle*)stack->Particle(resDau);
	  Int_t pdgresdau=resdau->GetPdgCode();
	  if(TMath::Abs(pdgresdau)==321){
	    if(pdgdp>0 && pdgresdau>0) return -1;
	    if(pdgdp<0 && pdgresdau<0) return -1;
	    nKaons++;
	  }
	  if(TMath::Abs(pdgresdau)==211){
	    nPions++;
	  }
	}
      }else if(TMath::Abs(pdgdau)==211){
	  nPions++;
      }else{
	return -1;
      }
    }
    if(nPions!=3) return -1;
    if(nKaons!=1) return -1;
    for(Int_t iDau=dp->GetFirstDaughter(); iDau<=dp->GetLastDaughter(); iDau++){
      if(iDau<0) return -1;
    }
    return 1;
  }
  
  return -1;
}


//______________________________________________________________________________
Int_t AliAnalysisTaskCheckHFMCProd::CheckDplusDecay(Int_t labDplus, AliStack* stack) const{

  if(labDplus<0) return -1;
  TParticle* dp = (TParticle*)stack->Particle(labDplus);
  Int_t pdgdp=dp->GetPdgCode();
  Int_t nDau=dp->GetNDaughters();

  if(nDau==3){
    Int_t nKaons=0;
    Int_t nPions=0;
    for(Int_t iDau=dp->GetFirstDaughter(); iDau<=dp->GetLastDaughter(); iDau++){
      if(iDau<0) return -1;
      TParticle* dau=(TParticle*)stack->Particle(iDau);
      Int_t pdgdau=dau->GetPdgCode();
      if(TMath::Abs(pdgdau)==321){
	if(pdgdp>0 && pdgdau>0) return -1;
	if(pdgdp<0 && pdgdau<0) return -1;
	nKaons++;
      }else if(TMath::Abs(pdgdau)==211){
	if(pdgdp<0 && pdgdau>0) return -1;
	if(pdgdp>0 && pdgdau<0) return -1;
	nPions++;
      }else{
	return -1;
      }
    }
    if(nPions!=2) return -1;
    if(nKaons!=1) return -1;
    for(Int_t iDau=dp->GetFirstDaughter(); iDau<=dp->GetLastDaughter(); iDau++){
      if(iDau<0) return -1;
    }
    return 0;
  }

  if(nDau==2){
    Int_t nKaons=0;
    Int_t nPions=0;
    for(Int_t iDau=dp->GetFirstDaughter(); iDau<=dp->GetLastDaughter(); iDau++){
      if(iDau<0) return -1;
      TParticle* dau=(TParticle*)stack->Particle(iDau);
      Int_t pdgdau=dau->GetPdgCode();
      if(TMath::Abs(pdgdau)==313){
	for(Int_t resDau=dau->GetFirstDaughter(); resDau<=dau->GetLastDaughter(); resDau++){
	  if(resDau<0) return -1;
	  TParticle* resdau=(TParticle*)stack->Particle(resDau);
	  Int_t pdgresdau=resdau->GetPdgCode();
	  if(TMath::Abs(pdgresdau)==321){
	    if(pdgdp>0 && pdgresdau>0) return -1;
	    if(pdgdp<0 && pdgresdau<0) return -1;
	    nKaons++;
	  }
	  if(TMath::Abs(pdgresdau)==211){
	    if(pdgdp<0 && pdgresdau>0) return -1;
	    if(pdgdp>0 && pdgresdau<0) return -1;
	    nPions++;
	  }
	}
      }else if(TMath::Abs(pdgdau)==211){
	  if(pdgdp<0 && pdgdau>0) return -1;
	  if(pdgdp>0 && pdgdau<0) return -1;
	  nPions++;
      }else{
	return -1;
      }
    }
    if(nPions!=2) return -1;
    if(nKaons!=1) return -1;
    for(Int_t iDau=dp->GetFirstDaughter(); iDau<=dp->GetLastDaughter(); iDau++){
      if(iDau<0) return -1;
    }
    return 1;
  }
  return -1;
}

//______________________________________________________________________________
Int_t AliAnalysisTaskCheckHFMCProd::CheckDsDecay(Int_t labDs, AliStack* stack) const{
  // Ds decay
  if(labDs<0) return -1;
  TParticle* dp = (TParticle*)stack->Particle(labDs);
  Int_t pdgdp=dp->GetPdgCode();
  Int_t nDau=dp->GetNDaughters();

  if(nDau==3){
    Int_t nKaons=0;
    Int_t nPions=0;
    for(Int_t iDau=dp->GetFirstDaughter(); iDau<=dp->GetLastDaughter(); iDau++){
      if(iDau<0) return -1;
      TParticle* dau=(TParticle*)stack->Particle(iDau);
      Int_t pdgdau=dau->GetPdgCode();
      if(TMath::Abs(pdgdau)==321){
	nKaons++;
      }else if(TMath::Abs(pdgdau)==211){
	if(pdgdp<0 && pdgdau>0) return -1;
	if(pdgdp>0 && pdgdau<0) return -1;
	nPions++;
      }else{
	return -1;
      }
    }
    if(nPions!=1) return -1;
    if(nKaons!=2) return -1;
    for(Int_t iDau=dp->GetFirstDaughter(); iDau<=dp->GetLastDaughter(); iDau++){
      if(iDau<0) return -1;
    }
    return 2;
  }

  if(nDau==2){
    Int_t nKaons=0;
    Int_t nPions=0;
    Bool_t isPhi=kFALSE;
    Bool_t isk0st=kFALSE;
    for(Int_t iDau=dp->GetFirstDaughter(); iDau<=dp->GetLastDaughter(); iDau++){
      if(iDau<0) return -1;
      TParticle* dau=(TParticle*)stack->Particle(iDau);
      Int_t pdgdau=dau->GetPdgCode();
      if(TMath::Abs(pdgdau)==313){
	isk0st=kTRUE;
	for(Int_t resDau=dau->GetFirstDaughter(); resDau<=dau->GetLastDaughter(); resDau++){
	  if(resDau<0) return -1;
	  TParticle* resdau=(TParticle*)stack->Particle(resDau);
	  Int_t pdgresdau=resdau->GetPdgCode();
	  if(TMath::Abs(pdgresdau)==321){
	    nKaons++;
	  }
	  if(TMath::Abs(pdgresdau)==211){
	    if(pdgdp<0 && pdgresdau>0) return -1;
	    if(pdgdp>0 && pdgresdau<0) return -1;
	    nPions++;
	  }
	}
      }else if(TMath::Abs(pdgdau)==333){
	isPhi=kTRUE;
	for(Int_t resDau=dau->GetFirstDaughter(); resDau<=dau->GetLastDaughter(); resDau++){
	  if(resDau<0) return -1;	  
	  TParticle* resdau=(TParticle*)stack->Particle(resDau);
	  if(!resdau) return -1;
	  Int_t pdgresdau=resdau->GetPdgCode();	  
	  if(TMath::Abs(pdgresdau)==321){
	    nKaons++;
	  }else{
	    return -1;
	  }
	}
      }else if(TMath::Abs(pdgdau)==211){
	if(pdgdp<0 && pdgdau>0) return -1;
	if(pdgdp>0 && pdgdau<0) return -1;
	nPions++;
      }else if(TMath::Abs(pdgdau)==321){
	nKaons++;
      }else{
	return -1;
      }
    }
    if(nPions!=1) return -1;
    if(nKaons!=2) return -1;
    for(Int_t iDau=dp->GetFirstDaughter(); iDau<=dp->GetLastDaughter(); iDau++){
      if(iDau<0) return -1;
    }
    if(isk0st) return 1;
    else if(isPhi) return 0;
    else return 3;
  } 
  return -1;
}

//______________________________________________________________________________
Int_t AliAnalysisTaskCheckHFMCProd::CheckDstarDecay(Int_t labDstar, AliStack* stack) const{

  if(labDstar<0) return -1;
  TParticle* dp = (TParticle*)stack->Particle(labDstar);
  Int_t pdgdp=dp->GetPdgCode();
  Int_t nDau=dp->GetNDaughters();
  if(nDau!=2) return -1;

  Int_t nKaons=0;
  Int_t nPions=0;
  for(Int_t iDau=dp->GetFirstDaughter(); iDau<=dp->GetLastDaughter(); iDau++){
    if(iDau<0) return -1;
    TParticle* dau=(TParticle*)stack->Particle(iDau);
    Int_t pdgdau=dau->GetPdgCode();
    if(TMath::Abs(pdgdau)==421){
      for(Int_t resDau=dau->GetFirstDaughter(); resDau<=dau->GetLastDaughter(); resDau++){
	if(resDau<0) return -1;
	TParticle* resdau=(TParticle*)stack->Particle(resDau);
	Int_t pdgresdau=resdau->GetPdgCode();
	if(TMath::Abs(pdgresdau)==321){
	  if(pdgdp>0 && pdgresdau>0) return -1;
	  if(pdgdp<0 && pdgresdau<0) return -1;
	  nKaons++;
	}
	if(TMath::Abs(pdgresdau)==211){
	  if(pdgdp<0 && pdgresdau>0) return -1;
	  if(pdgdp>0 && pdgresdau<0) return -1;
	  nPions++;
	}
      }
    }else if(TMath::Abs(pdgdau)==211){
      if(pdgdp<0 && pdgdau>0) return -1;
      if(pdgdp>0 && pdgdau<0) return -1;
      nPions++;
    }else{
      return -1;
    }
  }
  if(nPions!=2) return -1;
  if(nKaons!=1) return -1;
  for(Int_t iDau=dp->GetFirstDaughter(); iDau<=dp->GetLastDaughter(); iDau++){
    if(iDau<0) return -1;
  }
  return 0;  
  
}

//______________________________________________________________________________
Int_t AliAnalysisTaskCheckHFMCProd::CheckLcDecay(Int_t labLc, AliStack* stack) const{
  if(labLc<0) return -1;
  TParticle* dp = (TParticle*)stack->Particle(labLc);
  Int_t pdgdp=dp->GetPdgCode();
  Int_t nDau=dp->GetNDaughters();

  if(nDau==3){
    Int_t nKaons=0;
    Int_t nPions=0;
    Int_t nProtons=0;
    for(Int_t iDau=dp->GetFirstDaughter(); iDau<=dp->GetLastDaughter(); iDau++){
      if(iDau<0) return -1;
      TParticle* dau=(TParticle*)stack->Particle(iDau);
      Int_t pdgdau=dau->GetPdgCode();
      if(TMath::Abs(pdgdau)==321){
	if(pdgdp>0 && pdgdau>0) return -1;
	if(pdgdp<0 && pdgdau<0) return -1;
	nKaons++;
      }else if(TMath::Abs(pdgdau)==211){
	if(pdgdp<0 && pdgdau>0) return -1;
	if(pdgdp>0 && pdgdau<0) return -1;
	nPions++;
      }else if(TMath::Abs(pdgdau)==2212){
	if(pdgdp<0 && pdgdau>0) return -1;
	if(pdgdp>0 && pdgdau<0) return -1;
	nProtons++;
      }else{
	return -1;
      }
    }
    if(nPions!=1) return -1;
    if(nKaons!=1) return -1;
    if(nProtons!=1) return -1;
    for(Int_t iDau=dp->GetFirstDaughter(); iDau<=dp->GetLastDaughter(); iDau++){
      if(iDau<0) return -1;
    }
    return 0;
  }

  if(nDau==2){
    Int_t nKaons=0;
    Int_t nPions=0;
    Int_t nProtons=0;
    for(Int_t iDau=dp->GetFirstDaughter(); iDau<=dp->GetLastDaughter(); iDau++){
      if(iDau<0) return -1;
      TParticle* dau=(TParticle*)stack->Particle(iDau);
      Int_t pdgdau=dau->GetPdgCode();
      if(TMath::Abs(pdgdau)==313 || TMath::Abs(pdgdau)==3124 || TMath::Abs(pdgdau)==2224
	 || TMath::Abs(pdgdau)==3122 || TMath::Abs(pdgdau)==311){
	Int_t nDauRes=dau->GetNDaughters();
	if(nDauRes!=2)  return -1;
	for(Int_t resDau=dau->GetFirstDaughter(); resDau<=dau->GetLastDaughter(); resDau++){
	  if(resDau<0) return -1;
	  TParticle* resdau=(TParticle*)stack->Particle(resDau);
	  Int_t pdgresdau=resdau->GetPdgCode();
	  if(TMath::Abs(pdgresdau)==321){
	    if(pdgdp>0 && pdgresdau>0) return -1;
	    if(pdgdp<0 && pdgresdau<0) return -1;
	    nKaons++;
	  }
	  if(TMath::Abs(pdgresdau)==211){
	    if(pdgdp<0 && pdgresdau>0) return -1;
	    if(pdgdp>0 && pdgresdau<0) return -1;
	    nPions++;
	  }
	  if(TMath::Abs(pdgresdau)==2212){
	    if(pdgdp<0 && pdgresdau>0) return -1;
	    if(pdgdp>0 && pdgresdau<0) return -1;
	    nProtons++;
	  }
	}
      }else if(TMath::Abs(pdgdau)==321){
	if(pdgdp>0 && pdgdau>0) return -1;
	if(pdgdp<0 && pdgdau<0) return -1;
	nKaons++;
      }else if(TMath::Abs(pdgdau)==211){
	if(pdgdp<0 && pdgdau>0) return -1;
	if(pdgdp>0 && pdgdau<0) return -1;
	nPions++;
      }else if(TMath::Abs(pdgdau)==2212){
	if(pdgdp<0 && pdgdau>0) return -1;
	if(pdgdp>0 && pdgdau<0) return -1;
	nProtons++;
      }else{
	return -1;
      }
    }
    if(nPions!=1) return -1;
    if(nKaons!=1) return -1;
    if(nProtons!=1) return -1;
    for(Int_t iDau=dp->GetFirstDaughter(); iDau<=dp->GetLastDaughter(); iDau++){
      if(iDau<0) return -1;
    }
    return 1;
  } 
  return -1;
}

