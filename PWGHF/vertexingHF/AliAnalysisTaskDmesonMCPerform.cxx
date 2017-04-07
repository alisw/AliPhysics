/**************************************************************************
 * Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
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
// Class AliAnalysisTaskDmesonMCPerform
// AliAnalysisTaskSE for performance studies of 
// D meson hadronic decays in MC simulations
// F. Prino, prino@to.infn.it
//*************************************************************************

#include <TList.h>

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliAODVertex.h"
#include "AliAnalysisVertexingHF.h"
#include "AliVertexingHFUtils.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAnalysisTaskDmesonMCPerform.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonMCPerform);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskDmesonMCPerform::AliAnalysisTaskDmesonMCPerform():
  AliAnalysisTaskSE("taskDperf"),
  fOutput(0x0),
  fHistNEvents(0x0),
  fHistNGenD(0x0),
  fHistNCand(0x0),
  fAODProtection(1),
  fRDHFCuts(0x0),
  fRDHFCutsDplus(0x0)
{
  //
  /// Default constructor
  //

  fRDHFCuts=new AliRDHFCutsD0toKpi("EvSelCuts");
  fRDHFCuts->SetUsePhysicsSelection(kTRUE);
  fRDHFCuts->SetUseAnyTrigger();
  fRDHFCuts->SetTriggerClass("");
  fRDHFCuts->SetMaxVtxZ(10.);
  fRDHFCuts->SetOptPileup(AliRDHFCuts::kNoPileupSelection);
  fRDHFCuts->SetUseCentrality(AliRDHFCuts::kCentV0M);

  fPartName[0]="Dzero";
  fPartName[1]="Dplus";
  fPartName[2]="Dstar";
  fPartName[3]="Ds";
  fPartName[4]="Lc2pkpi";

  // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output
}

//________________________________________________________________________
AliAnalysisTaskDmesonMCPerform::~AliAnalysisTaskDmesonMCPerform()
{
  //
  /// Destructor
  //
  if(fOutput && !fOutput->IsOwner()){
    // delete histograms stored in fOutput
    delete fHistNEvents;
    delete fHistNGenD;
    delete fHistNCand;
    for(Int_t j=0; j<2*kDecays; j++){
      delete fHistPtYMultGen[j];
      delete fHistPtYMultGenDauInAcc[j];
      delete fHistPtYMultReco[j];
      delete fHistPtYMultRecoFilt[j];
      delete fHistPtYMultRecoSel[j];
      delete fHistXvtxResVsPt[j];
      delete fHistYvtxResVsPt[j];
      delete fHistZvtxResVsPt[j];
      delete fHistInvMassVsPt[j];
      delete fHistDecLenVsPt[j];
      delete fHistNormDLxyVsPt[j];
      delete fHistCosPointVsPt[j];
    }
  }
  delete fOutput;
  delete fRDHFCuts;
}
//________________________________________________________________________
void AliAnalysisTaskDmesonMCPerform::UserCreateOutputObjects()
{
  //
  /// Create the output container
  //
  if(fDebug > 1) printf("AliAnalysisTaskDmesonMCPerform::UserCreateOutputObjects() \n");
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  fHistNEvents = new TH1F("fHistNEvents", "number of events ",10,-0.5,9.5);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"nEvents read");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"Rejected due to mismatch in trees");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"nEvents with good AOD");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"Rejected due to trigger");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"Rejected due to vertex reco");
  fHistNEvents->GetXaxis()->SetBinLabel(6,"Rejected due to pileup");
  fHistNEvents->GetXaxis()->SetBinLabel(7,"Rejected due to centrality");
  fHistNEvents->GetXaxis()->SetBinLabel(8,"Rejected due to vtxz");
  fHistNEvents->GetXaxis()->SetBinLabel(9,"Rejected due to Physics Sel");
  fHistNEvents->GetXaxis()->SetBinLabel(10,"nEvents accepted");
  fHistNEvents->GetXaxis()->SetNdivisions(1,kFALSE);
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);

  fHistNGenD=new TH1F("fHistNGenD","number of D mesons",10,-0.5,9.5);
  TString names[5]={"Dzero","Dplus","Dstar","Ds","Lc"};
  TString type[2]={"Prompt","Feeddown"};
  for(Int_t j=0; j<5; j++){
    for(Int_t i=0; i<2; i++){
      Int_t index=j*2+i; 
      fHistNGenD->GetXaxis()->SetBinLabel(index+1,Form("%s %s",type[i].Data(),names[j].Data()));
    }
  }
  fHistNGenD->GetXaxis()->SetNdivisions(1,kFALSE);
  fHistNGenD->SetMinimum(0);
  fOutput->Add(fHistNGenD);

  fHistNCand=new TH2F("fHistNCand","number of D meson candidates",5,-0.5,4.5,16,0.,16.);
  fHistNCand->GetXaxis()->SetBinLabel(1,"D0#rightarrowK#pi");
  fHistNCand->GetXaxis()->SetBinLabel(2,"D+#rightarrowK#pi#pi");
  fHistNCand->GetXaxis()->SetBinLabel(3,"D*+#rightarrowD0#pi");
  fHistNCand->GetXaxis()->SetBinLabel(4,"D_s#rightarrowKK#pi");
  fHistNCand->GetXaxis()->SetBinLabel(5,"Lc#rightarrowpK#pi");
  fHistNCand->SetMinimum(0);
  fOutput->Add(fHistNCand);

  for(Int_t j=0; j<kDecays; j++){
    for(Int_t i=0; i<2; i++){
      Int_t index=j*2+i;
      fHistPtYMultGen[index]=new TH3F(Form("hPtYMult%sGen%s",type[i].Data(),fPartName[j].Data())," ; p_{T} (GeV/c) ; y; multPercentile",16,0.,16.,24,-6.,6.,40,0.,100.);
      fHistPtYMultGenDauInAcc[index]=new TH3F(Form("hPtYMult%sGenDauInAcc%s",type[i].Data(),fPartName[j].Data())," ; p_{T} (GeV/c) ; y; multPercentile",16,0.,16.,24,-6.,6.,40,0.,100.);
      fHistPtYMultReco[index]=new TH3F(Form("hPtYMult%sReco%s",type[i].Data(),fPartName[j].Data())," ; p_{T} (GeV/c) ; y; multPercentile",16,0.,16.,24,-6.,6.,40,0.,100.);
      fHistPtYMultRecoFilt[index]=new TH3F(Form("hPtYMult%sRecoFilt%s",type[i].Data(),fPartName[j].Data())," ; p_{T} (GeV/c) ; y; multPercentile",16,0.,16.,24,-6.,6.,40,0.,100.);
      fHistPtYMultRecoSel[index]=new TH3F(Form("hPtYMult%sRecoSel%s",type[i].Data(),fPartName[j].Data())," ; p_{T} (GeV/c) ; y; multPercentile",16,0.,16.,24,-6.,6.,40,0.,100.);
      fOutput->Add(fHistPtYMultGen[index]);
      fOutput->Add(fHistPtYMultGenDauInAcc[index]);
      fOutput->Add(fHistPtYMultReco[index]);
      fOutput->Add(fHistPtYMultRecoFilt[index]);
      fOutput->Add(fHistPtYMultRecoSel[index]);
      fHistXvtxResVsPt[index]=new TH2F(Form("hXvtxResVsPt%s%s",type[i].Data(),fPartName[j].Data())," ;  p_{T} (GeV/c) ; x_{v}(rec)-x_{v}(gen) (#mum)",16,0.,16.,100,-500.,500.);
      fHistYvtxResVsPt[index]=new TH2F(Form("hYvtxResVsPt%s%s",type[i].Data(),fPartName[j].Data())," ;  p_{T} (GeV/c) ; y_{v}(rec)-y_{v}(gen) (#mum)",16,0.,16.,100,-500.,500.);
      fHistZvtxResVsPt[index]=new TH2F(Form("hZvtxResVsPt%s%s",type[i].Data(),fPartName[j].Data())," ;  p_{T} (GeV/c) ; z_{v}(rec)-z_{v}(gen) (#mum)",16,0.,16.,100,-500.,500.);
      fHistInvMassVsPt[index]=new TH2F(Form("hInvMassVsPt%s%s",type[i].Data(),fPartName[j].Data())," ;  p_{T} (GeV/c) ; Inv. Mass (GeV/c^{2})",16,0.,16.,300,1.75,2.35);
      fHistDecLenVsPt[index]=new TH2F(Form("hDecLenVsPt%s%s",type[i].Data(),fPartName[j].Data())," ;  p_{T} (GeV/c) ; dec. len. (#mum)",16,0.,16.,100,0.,5000.);
      fHistNormDLxyVsPt[index]=new TH2F(Form("hNormDLxyVsPt%s%s",type[i].Data(),fPartName[j].Data())," ;  p_{T} (GeV/c) ; dec. len. (#mum)",16,0.,16.,100,0.,20.);
      fHistCosPointVsPt[index]=new TH2F(Form("hCosPointVsPt%s%s",type[i].Data(),fPartName[j].Data())," ;  p_{T} (GeV/c) ; cos(#vartheta_{point})",16,0.,16.,100,-1.,1.);
      fOutput->Add(fHistXvtxResVsPt[index]);
      fOutput->Add(fHistYvtxResVsPt[index]);
      fOutput->Add(fHistZvtxResVsPt[index]);
      fOutput->Add(fHistInvMassVsPt[index]);
      fOutput->Add(fHistDecLenVsPt[index]);
      fOutput->Add(fHistNormDLxyVsPt[index]);
      fOutput->Add(fHistCosPointVsPt[index]);
    }
  }


  PostData(1,fOutput);
  return;
}
//________________________________________________________________________
void AliAnalysisTaskDmesonMCPerform::UserExec(Option_t */*option*/)
{
  /// Execute analysis for current event:
  /// heavy flavor candidates association to MC truth

  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());

  fHistNEvents->Fill(0); // count event

  if(fAODProtection>=0){
    //   Protection against different number of events in the AOD and deltaAOD
    //   In case of discrepancy the event is rejected.
    Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      fHistNEvents->Fill(1);
      return;
    }
  }

  // Load all the branches of the DeltaAOD - needed for SelectionBit counting
  TClonesArray *arrayD0toKpi  =0;
  TClonesArray *array3Prong   =0;
  TClonesArray *arrayDstar    =0;
  TClonesArray *arrayCascades =0;
  if(!aod && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD
    // event in memory rather than the input (ESD) event.
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*)
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {

      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();

      array3Prong  =(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
      arrayD0toKpi =(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
      arrayDstar   =(TClonesArray*)aodFromExt->GetList()->FindObject("Dstar");
      arrayCascades=(TClonesArray*)aodFromExt->GetList()->FindObject("CascadesHF");
    }
  } else if(aod) {
    array3Prong  =(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
    arrayD0toKpi =(TClonesArray*)aod->GetList()->FindObject("D0toKpi");
    arrayDstar   =(TClonesArray*)aod->GetList()->FindObject("Dstar");
    arrayCascades=(TClonesArray*)aod->GetList()->FindObject("CascadesHF");
  }

  if(!aod || !array3Prong || !arrayD0toKpi) {
    printf("AliAnalysisTaskSEDplus::UserExec: Charm3Prong branch not found!\n");
    return;
  }
  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex()||TMath::Abs(aod->GetMagneticField())<0.001) return;
  fHistNEvents->Fill(2);

  //Event selection
  Bool_t isEvSel=fRDHFCuts->IsEventSelected(aod);
  if(fRDHFCuts->GetWhyRejection()==5) fHistNEvents->Fill(3);
  if(!isEvSel && fRDHFCuts->GetWhyRejection()==0) fHistNEvents->Fill(4);
  if(fRDHFCuts->GetWhyRejection()==1) fHistNEvents->Fill(5);
  if(fRDHFCuts->GetWhyRejection()==2) fHistNEvents->Fill(6);
  if(fRDHFCuts->GetWhyRejection()==6)fHistNEvents->Fill(7);
  if(fRDHFCuts->GetWhyRejection()==7)fHistNEvents->Fill(8);
 
  // load MC header
  AliAODMCHeader *mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  if(!mcHeader) {
    printf("AliAnalysisTaskDmesonMCPerform:UserExec: MC header branch not found!\n");
    return;
  }

  // load MC particles
  TClonesArray *arrayMC=  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if(!arrayMC) {
    printf("AliAnalysisTaskDmesonMCPerform::UserExec: MC particles branch not found!\n");
    return;
  }

  Int_t runNumber=aod->GetRunNumber();

  if(aod->GetTriggerMask()==0 &&
     (runNumber>=195344 && runNumber<=195677)){
    // protection for events with empty trigger mask in p-Pb
    return;
  }
  if(fRDHFCuts->GetUseCentrality()>0 && fRDHFCuts->IsEventSelectedInCentrality(aod)!=0) return;

  PostData(1,fOutput);

  for(Int_t i=0; i<kMaxLabel; i++) fMapTrLabel[i]=-999;
  MapTrackLabels(aod);
  FillGenLevelHistos(aod,arrayMC,mcHeader);
  if(!isEvSel)return;
  fHistNEvents->Fill(9);

  FillCandLevelHistos(0,aod,arrayD0toKpi,arrayMC);
  FillCandLevelHistos(1,aod,array3Prong,arrayMC);
  FillCandLevelHistos(2,aod,arrayDstar,arrayMC);

  return;
}
//________________________________________________________________________
void AliAnalysisTaskDmesonMCPerform::FillGenLevelHistos(AliAODEvent* aod, TClonesArray *arrayMC, AliAODMCHeader *mcHeader){
  //
  /// fill histograms starting from generated particles

  Double_t zMCVertex = mcHeader->GetVtxZ(); //vertex MC
  if(TMath::Abs(zMCVertex)>fRDHFCuts->GetMaxVtxZ()) return;
  Double_t evCentr=1.;
  if(fRDHFCuts->GetUseCentrality()>0) evCentr=fRDHFCuts->GetCentrality(aod);

  for(Int_t iPart=0; iPart<arrayMC->GetEntriesFast(); iPart++){
    AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(arrayMC->At(iPart));
    Int_t pdgCode=TMath::Abs(mcPart->GetPdgCode());
    Int_t iSpec=-1;
    if(pdgCode == 411) iSpec=1;
    else if(pdgCode == 413) iSpec=2;
    else if(pdgCode == 421) iSpec=0;
    else if(pdgCode == 431) iSpec=3;
    else if(pdgCode == 4122) iSpec=4;
    if (iSpec>=0){
      Int_t orig=AliVertexingHFUtils::CheckOrigin(arrayMC,mcPart,kTRUE);//Prompt = 4, FeedDown = 5
      if(orig<4 || orig>5) continue;
      Int_t indexh=iSpec*2+(orig-4);
      fHistNGenD->Fill(indexh);
      Int_t deca=0;
      Bool_t isGoodDecay=kFALSE;
      Int_t labDau[4]={-1,-1,-1,-1};
      Int_t nProng=3;
      if (pdgCode == 411){
	deca=AliVertexingHFUtils::CheckDplusDecay(arrayMC,mcPart,labDau);
        if(deca>0) isGoodDecay=kTRUE;
      }else if(pdgCode == 421){
	deca=AliVertexingHFUtils::CheckD0Decay(arrayMC,mcPart,labDau);
        if(mcPart->GetNDaughters()!=2) continue;
        if(deca==1) isGoodDecay=kTRUE;
	nProng=2;
      }else if(pdgCode == 431){
        deca=AliVertexingHFUtils::CheckDsDecay(arrayMC,mcPart,labDau);
        if(deca==1) isGoodDecay=kTRUE;	
      }else if(pdgCode==413){
	deca=AliVertexingHFUtils::CheckDstarDecay(arrayMC,mcPart,labDau);
	if(deca==1) isGoodDecay=kTRUE;
      }else if(pdgCode==4122){
 	deca=AliVertexingHFUtils::CheckLcpKpiDecay(arrayMC,mcPart,labDau);
	if(deca>=1) isGoodDecay=kTRUE;
     }
      if(labDau[0]==-1){
	continue; //protection against unfilled array of labels
      }
      if(isGoodDecay && indexh>=0){
        Double_t ptgen=mcPart->Pt();
        Double_t ygen=mcPart->Y();
	fHistPtYMultGen[indexh]->Fill(ptgen,ygen,evCentr);
	if(fRDHFCuts->IsInFiducialAcceptance(ptgen,ygen)){
	  Bool_t dauInAcc=CheckAcceptance(arrayMC,nProng,labDau);
	  if(dauInAcc){
	    fHistPtYMultGenDauInAcc[indexh]->Fill(ptgen,ygen,evCentr);
	    AliAODRecoDecayHF* dmes=GetRecoDecay(aod,nProng,labDau);
	    if(dmes){
	      fHistPtYMultReco[indexh]->Fill(ptgen,ygen,evCentr);
	    }
	    delete dmes;
	  }
	}
      }
    }
  }
}
//________________________________________________________________________
void AliAnalysisTaskDmesonMCPerform::FillCandLevelHistos(Int_t idCase, AliAODEvent* aod, TClonesArray *arrayDcand, TClonesArray *arrayMC){
  //
  /// fill histograms starting from candidates in deltaAODs.

  Double_t evCentr=1.;
  if(fRDHFCuts->GetUseCentrality()>0) evCentr=fRDHFCuts->GetCentrality(aod);

  Int_t pdg0[2]={321,211};
  Int_t pdgp[3]={321,211,211};
  Int_t pdgs[3]={321,211,321};
  Int_t pdgst[2]={421,211};
  Int_t pdgl[3]={2212,321,211};
  // vHF object is needed to call the method that refills the missing info of the candidates
  // if they have been deleted in dAOD reconstruction phase
  // in order to reduce the size of the file
  AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();

  Int_t nCand=arrayDcand->GetEntriesFast();
  for (Int_t iCand = 0; iCand < nCand; iCand++) {
    AliAODRecoDecayHF *d=(AliAODRecoDecayHF*)arrayDcand->UncheckedAt(iCand);
    Double_t ptCand=-999.;
    Int_t iSpec=-1;
    Int_t labD=-1;
    Double_t invMass=0.;
    if(idCase==0){
      if(!vHF->FillRecoCand(aod,(AliAODRecoDecayHF2Prong*)d))continue;
      ptCand=d->Pt();
      fHistNCand->Fill(0.,ptCand);
      labD = d->MatchToMC(421,arrayMC,2,pdg0);
      if(labD>=0){ 
	iSpec=0;
	AliAODMCParticle *pd0 = (AliAODMCParticle*)arrayMC->At(labD);
	if(pd0->GetPdgCode()==421) invMass=((AliAODRecoDecayHF2Prong*)d)->InvMassD0();
	else invMass=((AliAODRecoDecayHF2Prong*)d)->InvMassD0bar();
      }
    }else if(idCase==1){
      if(!vHF->FillRecoCand(aod,(AliAODRecoDecayHF3Prong*)d))continue;
      ptCand=d->Pt();
      if(d->HasSelectionBit(AliRDHFCuts::kDplusCuts)) fHistNCand->Fill(1.,ptCand);
      if(d->HasSelectionBit(AliRDHFCuts::kDsCuts)) fHistNCand->Fill(3.,ptCand);
      if(d->HasSelectionBit(AliRDHFCuts::kLcCuts)) fHistNCand->Fill(4.,ptCand);
      labD = d->MatchToMC(411,arrayMC,3,pdgp);
      if(labD>=0){
	if(d->HasSelectionBit(AliRDHFCuts::kDplusCuts)) iSpec=1;
	invMass=((AliAODRecoDecayHF3Prong*)d)->InvMassDplus();
      }else{
	Int_t labDau0=((AliAODTrack*)d->GetDaughter(0))->GetLabel();
	AliAODMCParticle* pdau0=(AliAODMCParticle*)arrayMC->UncheckedAt(TMath::Abs(labDau0));
	Int_t pdgCode0=TMath::Abs(pdau0->GetPdgCode());
	labD = d->MatchToMC(431,arrayMC,3,pdgs);
	if(labD>=0){
	  if(d->HasSelectionBit(AliRDHFCuts::kDsCuts)) iSpec=3;
	  if(pdgCode0==321) invMass=((AliAODRecoDecayHF3Prong*)d)->InvMassDsKKpi();
	  else invMass=((AliAODRecoDecayHF3Prong*)d)->InvMassDspiKK();
	}else{
	  labD = d->MatchToMC(4122,arrayMC,3,pdgl);
	  if(labD>=0){
	    if(d->HasSelectionBit(AliRDHFCuts::kLcCuts)) iSpec=4;
	    if(pdgCode0==2212) invMass=((AliAODRecoDecayHF3Prong*)d)->InvMassLcpKpi();
	    else invMass=((AliAODRecoDecayHF3Prong*)d)->InvMassLcpiKp();
	  }
	}
      }
    }else if(idCase==2){
      if(!vHF->FillRecoCasc(aod,((AliAODRecoCascadeHF*)d),kTRUE))continue;
      ptCand=d->Pt();
      fHistNCand->Fill(2.,ptCand);
      labD = ((AliAODRecoCascadeHF*)d)->MatchToMC(413,421,pdgst,pdg0,arrayMC);
      if(labD>=0){
	iSpec=2;
	invMass=((AliAODRecoCascadeHF*)d)->InvMassDstarKpipi();
      }
    }
    
    if(labD>=0 && iSpec>=0){
      AliAODMCParticle *partD = (AliAODMCParticle*)arrayMC->At(labD);
      Double_t ptgen=partD->Pt();
      Double_t ygen=partD->Y();
      Double_t ptreco=d->Pt();
      Double_t dlen=d->DecayLength()*10000.; //um
      Double_t ndlenxy=d->NormalizedDecayLengthXY();
      Double_t cosp=d->CosPointingAngle();
      Double_t dx=(d->Xv()-partD->Xv())*10000.;
      Double_t dy=(d->Yv()-partD->Yv())*10000.;
      Double_t dz=(d->Zv()-partD->Zv())*10000.;
      Int_t orig=AliVertexingHFUtils::CheckOrigin(arrayMC,partD,kTRUE);//Prompt = 4, FeedDown = 5
      if(orig<4 || orig>5) continue;
      Int_t indexh=iSpec*2+(orig-4);
      fHistPtYMultRecoFilt[indexh]->Fill(ptgen,ygen,evCentr);
      if(iSpec==1){
	if(fRDHFCutsDplus){
	  if(fRDHFCutsDplus->IsEventSelected(aod) &&
	     fRDHFCutsDplus->IsSelected(d,AliRDHFCuts::kAll,aod)){
	    fHistPtYMultRecoSel[indexh]->Fill(ptgen,ygen,evCentr);
	  }
	}
      }
      fHistXvtxResVsPt[indexh]->Fill(ptreco,dx);
      fHistYvtxResVsPt[indexh]->Fill(ptreco,dy);
      fHistZvtxResVsPt[indexh]->Fill(ptreco,dz);
      fHistInvMassVsPt[indexh]->Fill(ptreco,invMass);
      fHistDecLenVsPt[indexh]->Fill(ptreco,dlen);
      fHistNormDLxyVsPt[indexh]->Fill(ptreco,ndlenxy);
      fHistCosPointVsPt[indexh]->Fill(ptreco,cosp);
    }
  }
  delete vHF;
  return;
}
//_________________________________________________________________
Bool_t AliAnalysisTaskDmesonMCPerform::CheckAcceptance(TClonesArray* arrayMC,Int_t nProng, Int_t *labDau){
  /// check if the decay products are in the good eta and pt range
  for (Int_t iProng = 0; iProng<nProng; iProng++){
    AliAODMCParticle* mcPartDaughter=dynamic_cast<AliAODMCParticle*>(arrayMC->At(labDau[iProng]));
    if(!mcPartDaughter) return kFALSE;
    Double_t eta = mcPartDaughter->Eta();
    Double_t pt = mcPartDaughter->Pt();
    if (TMath::Abs(eta) > 0.9 || pt < 0.1) return kFALSE;
  }
  return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskDmesonMCPerform::MapTrackLabels(AliAODEvent* aod){
  /// Fill array of correspondence track lables <-> id
  //
  Int_t nTracks=aod->GetNumberOfTracks();

  for(Int_t it=0; it<nTracks; it++) {
    AliAODTrack *tr=dynamic_cast<AliAODTrack*>(aod->GetTrack(it));
    if(!tr) continue;
    if(tr->GetID()<0) continue;
    Int_t lab=TMath::Abs(tr->GetLabel());
    if(lab<kMaxLabel) fMapTrLabel[lab]=it;
    else printf("Label %d exceeds upper limit\n",lab);
  }
  return;
}
//________________________________________________________________________
AliAODRecoDecayHF* AliAnalysisTaskDmesonMCPerform::GetRecoDecay(AliAODEvent* aod, Int_t nProng, Int_t *labDau){
  //
  /// create the AliAODRecoDecayHF object fromt he tracks
  //
  Double_t px[3],py[3],pz[3],d0[3],d0err[3];
  Int_t nFound=0;
  for(Int_t jd=0; jd<nProng; jd++){
    Int_t it=fMapTrLabel[labDau[jd]];
    if(it>0){
      AliAODTrack* tr=dynamic_cast<AliAODTrack*>(aod->GetTrack(it));
      if(tr){
	if(TMath::Abs(tr->GetLabel())!=labDau[jd]){
	  printf("Mismatched track label\n");
	  continue;
	}
	px[jd]=tr->Px();
	py[jd]=tr->Py();
	px[jd]=tr->Pz();
	d0[jd]=0.; // temporary
	d0err[jd]=0.; // temporary
	nFound++;
      }
    }
  }
  if(nFound!=nProng) return 0x0;
  AliAODRecoDecayHF* dmes=new AliAODRecoDecayHF(0x0,nProng,nProng-2,px,py,pz,d0,d0err);
  return dmes;
}
//________________________________________________________________________
void AliAnalysisTaskDmesonMCPerform::Terminate(Option_t */*option*/)
{
  /// Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskDmesonMCPerform: Terminate() \n");

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {
    printf("ERROR: fOutput not available\n");
    return;
  }

  fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEvents"));
  if(fHistNEvents){
    printf("Number of analyzed events = %d\n",(Int_t)fHistNEvents->GetBinContent(10));
  }else{
    printf("ERROR: fHistNEvents not available\n");
    return;
  }

  return;
}
