/*************************************************************************
 * Copyright(c) 1998-2011, ALICE Experiment at CERN, All rights reserved. *
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

#include <iostream>
#include <TList.h>
#include <TString.h>
#include <TDatabasePDG.h>

#include <TObjArray.h>
#include <TRandom.h>
#include <TClonesArray.h>
#include <TGraph.h>
#include <TFile.h>
#include <TList.h>
#include <TNtuple.h>


#include "AliVertex.h"
#include "AliVVertex.h"
#include "AliESDVertex.h"
#include "AliVertexerTracks.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisVertexingHF.h"
#include "AliExternalTrackParam.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliVertexingHFUtils.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandlerRP.h"
#include "AliAnalysisTaskSELbtoLcpi4.h"

using namespace std;

ClassImp(AliAnalysisTaskSELbtoLcpi4);
//
//  Analysis Task for Lb in Lc and pion
//

AliAnalysisTaskSELbtoLcpi4::AliAnalysisTaskSELbtoLcpi4()
:AliAnalysisTaskSE(),
  fPIDResponse(0),
  fOutput(0),
  fHistNEvents(0),
  fRDCutsAnalysisLc(0),
  fRDCutsProductionLb(0),
  fListCuts(0),
  fBzkG(0.),
  fvtx1(0x0),
  fNtupleLambdabUPG(0),
  fFillNtupleSignal(kFALSE),
  fFillNtupleBackgroundRotated(kFALSE),
  fFillNtupleBackgroundNonRotated(kFALSE),
  fCutsond0Lcdaughters(kFALSE),
  fIsPromptLc(kFALSE),
  fApplyFixesITS3AnalysisBit(kFALSE),
  fApplyFixesITS3AnalysiskAll(kFALSE),
  fApplyFixesITS3AnalysisHijing(kFALSE),
  fPreSelectLctopKpi(kFALSE)
{
  //
  // Default constructor.
  //
  for(Int_t icut=0; icut<2; icut++) fCutD0Daughter[icut]=0.;
  for(Int_t icut=0; icut<7; icut++) fCutsPerPt[icut]=0.;
  fNRotations=13.;

}

AliAnalysisTaskSELbtoLcpi4::AliAnalysisTaskSELbtoLcpi4(const char *name,
    Bool_t fillNtuple,
    AliRDHFCutsLctopKpi *lccutsana,
    AliRDHFCutsLctopKpi *lccutsprod,
    Int_t ndebug)
:AliAnalysisTaskSE(name),
  fPIDResponse(0),
  fOutput(0),
  fHistNEvents(0),
  fRDCutsAnalysisLc(lccutsana),
  fRDCutsProductionLb(lccutsprod),
  fListCuts(0),
  fBzkG(0.),
  fvtx1(0x0),
  fNtupleLambdabUPG(0),
  fFillNtupleSignal(kFALSE),
  fFillNtupleBackgroundRotated(kFALSE),
  fFillNtupleBackgroundNonRotated(kFALSE),
  fCutsond0Lcdaughters(kFALSE),
  fIsPromptLc(kFALSE),
  fApplyFixesITS3AnalysisBit(kFALSE),
  fApplyFixesITS3AnalysiskAll(kFALSE),
  fApplyFixesITS3AnalysisHijing(kFALSE),
  fPreSelectLctopKpi(kFALSE)
{
  //SetPtBinLimit(fRDCutsAnalysisLc->GetNPtBins()+1,fRDCutsAnalysisLc->GetPtBinLimits());
    for(Int_t icut=0; icut<2; icut++) fCutD0Daughter[icut]=0.;
    for(Int_t icut=0; icut<7; icut++) fCutsPerPt[icut]=0.;
    fNRotations=13.;
  //
  // Constructor to be used to create the task.
  // The the URIs specify the resolution files to be used.
  // They are expected to contain TGraphs with the resolutions
  // for the current and the upgraded ITS (see code for details).
  // One may also specify for how many tracks debug information
  // is written to the output.
  //
  DefineOutput(1,TList::Class());
  DefineOutput(2,TNtuple::Class());
}

AliAnalysisTaskSELbtoLcpi4::~AliAnalysisTaskSELbtoLcpi4() {
  //
  // Destructor.
  //

  delete fPIDResponse;
  delete fOutput;
  delete fHistNEvents;
  delete fRDCutsAnalysisLc;
  delete fRDCutsProductionLb;
  delete fListCuts;
  delete fvtx1;
  delete fNtupleLambdabUPG;

}
//-----------------------------------------------------
void AliAnalysisTaskSELbtoLcpi4::Init()
{
  // Initialization

  fListCuts=new TList();
  fListCuts->Add(new AliRDHFCutsLctopKpi(*fRDCutsAnalysisLc));
  fListCuts->Add(new AliRDHFCutsLctopKpi(*fRDCutsProductionLb));

  return;
}

//-----------------------------------------

void AliAnalysisTaskSELbtoLcpi4::UserExec(Option_t*) {
  //
  // The event loop
  //

  AliAODEvent *aod=dynamic_cast<AliAODEvent*>(InputEvent());
  fHistNEvents->Fill(0);

  TClonesArray *array3Prong = 0;
  TClonesArray *arrayLikeSign =0;
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
      array3Prong=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
      arrayLikeSign=(TClonesArray*)aodFromExt->GetList()->FindObject("LikeSign3Prong");
    }
  } else if(aod) {
    array3Prong=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
    arrayLikeSign=(TClonesArray*)aod->GetList()->FindObject("LikeSign3Prong");
  }
  if(!aod || !array3Prong) return;

  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  fBzkG = aod->GetMagneticField();
  fvtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  if(!fvtx1 || TMath::Abs(fBzkG)<0.001) return;

  AliAODMCHeader *mcHeader = 0;
  mcHeader = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  if(!mcHeader) {
    printf(" MC header branch not found!\n");
    return;
  }
  TClonesArray *mcs=static_cast<TClonesArray*>(aod->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
  if (!mcs) return;

  AliAnalysisVertexingHF *vHF=new AliAnalysisVertexingHF();
  //  loop 3 prongs
  Int_t n3Prong = array3Prong->GetEntriesFast();
  for (Int_t icand = 0; icand < n3Prong; icand++) {
    AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(icand);

    if(fApplyFixesITS3AnalysisBit){
      if(!(d->HasSelectionBit(AliRDHFCuts::kLcCuts))) continue;
    }
      
    if(fPreSelectLctopKpi){
       TObjArray arrTracks(3);
       for(Int_t ipr=0;ipr<3;ipr++){
           AliAODTrack *tr=vHF->GetProng(aod,d,ipr);
           arrTracks.AddAt(tr,ipr);
       }
       Bool_t recoLc=kTRUE;
       recoLc=fRDCutsAnalysisLc->PreSelectMass(arrTracks);
       if (!recoLc) continue;
    }

    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()){
      d->SetOwnPrimaryVtx(fvtx1);
      unsetvtx=kTRUE;
    }

    if((vHF->FillRecoCand(aod,d))) {
    Int_t selection;
    //is selected for LambdaC
    if(fApplyFixesITS3AnalysiskAll)selection=fRDCutsAnalysisLc->IsSelected(d,AliRDHFCuts::kAll,aod);
    else                           selection=fRDCutsAnalysisLc->IsSelected(d,AliRDHFCuts::kCandidate,aod);
    if(selection==0){
      if(unsetvtx) d->UnsetOwnPrimaryVtx();
      continue;
    }

    //lc preliminary large pt cuts:
    if(d->Pt()>fCutsPerPt[1] || d->Pt()<fCutsPerPt[2]){
      if(unsetvtx) d->UnsetOwnPrimaryVtx();
      continue;
    }

    //Additional Cut on Lc
    // d0p and d0pi of Lc
    //large pt cuts on the d0's of the Lc daughters //keep them out for the moment
    if (fCutsond0Lcdaughters){
      if(TMath::Abs(d->Getd0Prong(0))<fCutD0Daughter[0] || (TMath::Abs(d->Getd0Prong(2))<fCutD0Daughter[1])){
        if(unsetvtx) d->UnsetOwnPrimaryVtx();
        continue;
      }
    }

    FillHistos(d,mcs,aod,mcHeader);
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  }
}
  delete vHF;
  return;
}

void AliAnalysisTaskSELbtoLcpi4::FillHistos(AliAODRecoDecayHF3Prong* d,TClonesArray* arrayMC,AliAODEvent *ev,AliAODMCHeader *mcHeader){
  Int_t countLc=0;

  //check ID d prongs
  Int_t idProng1 = d->GetProngID(0);
  Int_t idProng2 = d->GetProngID(1);
  Int_t idProng3 = d->GetProngID(2);

  Int_t lc=0;
  fIsPromptLc=kFALSE;
  Int_t labLb=CheckMCLc(d,arrayMC);//after track cuts
  //here we know if it is a prompt Lc
  AliExternalTrackParam *LcCand = new AliExternalTrackParam;
  LcCand->CopyFromVTrack(d);

  for(Int_t k = 0 ; k < ev->GetNumberOfTracks() ; k++) {
    Double_t xdummy,ydummy;
    Double_t dzdummy[2];
    Double_t covardummy[3];

    AliAODTrack * HPiAODtrk = dynamic_cast<AliAODTrack*>(ev->GetTrack(k));
    if(HPiAODtrk->GetID()==idProng1 || HPiAODtrk->GetID()==idProng2 || HPiAODtrk->GetID()==idProng3){
      HPiAODtrk=0;
      continue;
    }
    if(HPiAODtrk->Charge()==d->Charge()){
      HPiAODtrk=0;
      continue;
    }
    if(HPiAODtrk->Charge()==0){
      HPiAODtrk=0;
      continue;
    }
    if(HPiAODtrk->Pt()>fCutsPerPt[3] || HPiAODtrk->Pt()<fCutsPerPt[4]){
      HPiAODtrk=0;
      continue;
    }

    //basic PID pion
    if(fRDCutsAnalysisLc->GetIsUsePID()){
      Double_t nsigmatofPi= fPIDResponse->NumberOfSigmasTOF(HPiAODtrk,AliPID::kPion);
      if(nsigmatofPi>-990. && (nsigmatofPi<-3 || nsigmatofPi>3)){HPiAODtrk=0;continue;}
      Double_t nsigmatpcPi= fPIDResponse->NumberOfSigmasTPC(HPiAODtrk,AliPID::kPion);
      if(nsigmatpcPi>-990. && (nsigmatpcPi<-3 || nsigmatpcPi>3)){HPiAODtrk=0;continue;}
    }

    // Track cuts
    AliESDtrackCuts *trackCutsHPi = fRDCutsAnalysisLc->GetTrackCuts();//
    if(!fRDCutsAnalysisLc->IsDaughterSelected(HPiAODtrk,(AliESDVertex*)fvtx1,trackCutsHPi)){
      HPiAODtrk=0;
      continue;
    }

    AliExternalTrackParam *chargedHPi = new AliExternalTrackParam;
    chargedHPi->CopyFromVTrack(HPiAODtrk);

    //further cuts on candidate charged track

    const Double_t max = 1;
    Double_t d0cut[2],covd0cut[3];
    Double_t d0cutLc[2],covd0cutLc[3];
    chargedHPi->PropagateToDCA(fvtx1,fBzkG,max,d0cut,covd0cut);
    LcCand->PropagateToDCA(fvtx1,fBzkG,max,d0cutLc,covd0cutLc);

    // Construction of Secondary vertex
    TObjArray *recoArray = new TObjArray(2);
    Double_t dispersion;

    //cri mod
    recoArray->AddAt(chargedHPi,0);//pion is the second Prong
    recoArray->AddAt(LcCand,1);//Lc is the first Prong

    AliAODVertex *vtxAODNew = ReconstructSecondaryVertex(recoArray,dispersion,kFALSE);
    if(!vtxAODNew) {
      recoArray->Clear();
      HPiAODtrk=0;
      delete chargedHPi;
      delete recoArray;
      continue;
    }
    // Add daughter information
    AddDaughterRefs(vtxAODNew,(AliAODEvent*)ev,recoArray);

    //______________________________________________________________________________________
    // construction of lb (with secondary vertex)

    const Double_t maxd = 1;
    //___
    // Propagate candidates to secondary vertex
    chargedHPi->PropagateToDCA(vtxAODNew,fBzkG,maxd,dzdummy,covardummy);
    LcCand->PropagateToDCA(vtxAODNew,fBzkG,maxd,dzdummy,covardummy);

    // Calculate momenta
    Double_t momentum[3];
    chargedHPi->GetPxPyPz(momentum);

    Double_t px[2],py[2],pz[2],d0[2],d0err[2],dcaCand;
    px[0] = momentum[0]; py[0] = momentum[1]; pz[0] = momentum[2];
    LcCand->GetPxPyPz(momentum);
    px[1] = momentum[0]; py[1] = momentum[1]; pz[1] = momentum[2];

    // Calculate impact parameters
    Double_t d0z0[2],covd0z0[3];
    LcCand->PropagateToDCA(fvtx1,fBzkG,maxd,d0z0,covd0z0);
    d0[1] = d0z0[0];
    d0err[1] = TMath::Sqrt(covd0z0[0]);
    chargedHPi->PropagateToDCA(fvtx1,fBzkG,maxd,d0z0,covd0z0);
    d0[0] = d0z0[0];
    d0err[0] = TMath::Sqrt(covd0z0[0]);

    dcaCand = chargedHPi->GetDCA(LcCand,fBzkG,xdummy,ydummy);

    // Create lbcandidate as AliAODRecoDecayHF2Prong
    AliAODRecoDecayHF2Prong *lbcandProng = new AliAODRecoDecayHF2Prong(vtxAODNew,px,py,pz,d0,d0err,dcaCand);
    if(lbcandProng->Pt()>fCutsPerPt[5]||lbcandProng->Pt()<fCutsPerPt[6]){
      HPiAODtrk=0;
      delete chargedHPi;
      recoArray->Clear();
      delete recoArray;
      if(vtxAODNew){delete vtxAODNew;vtxAODNew=NULL;}
      delete lbcandProng;
      continue;
    }
    lbcandProng->SetOwnPrimaryVtx(fvtx1);

    UShort_t id[2];
    id[0]=(UShort_t)HPiAODtrk->GetID();
    id[1]=(UShort_t)LcCand->GetID();
    lbcandProng->SetProngIDs(2,id);
    AddDaughterRefs(vtxAODNew,(AliAODEvent*)ev,recoArray);
    Int_t idHardPi=(Int_t)HPiAODtrk->GetID();
    if (idHardPi > -1) {
      vtxAODNew->AddDaughter(HPiAODtrk);
    }
    vtxAODNew->AddDaughter(LcCand);

    Int_t lb=0;
    Int_t labPi2=CheckMCpartPIONaf(HPiAODtrk,arrayMC);
    if(labPi2==labLb){//signal
      lb=1;
      lc=1;
    }

    Bool_t isHijing = CheckGenerator(HPiAODtrk,d,mcHeader,arrayMC);
    // JJJ - check whether bkg from hijing
    if(lb==0 && !isHijing){
      HPiAODtrk=0;
      delete chargedHPi;
      recoArray->Clear();
      delete recoArray;
      if(vtxAODNew){delete vtxAODNew;vtxAODNew=NULL;}
      delete lbcandProng;
      continue;
    }

    //fill not rotated

    Int_t selectionlb=IsSelectedLbMY(lbcandProng,AliRDHFCuts::kCandidate,lb,0,isHijing);
    if(selectionlb!=0){
      Int_t dgLabelsnr[3];
      for(Int_t i=0; i<3; i++) {
        AliAODTrack *trknr = (AliAODTrack*)d->GetDaughter(i);
        dgLabelsnr[i] = trknr->GetLabel();
      }
      Int_t LabelPionnr= HPiAODtrk->GetLabel();
      FillLbHistsnr(lbcandProng,lb,mcHeader,arrayMC,HPiAODtrk,d, lc, ev, fIsPromptLc);
    }
    lbcandProng->UnsetOwnPrimaryVtx();

    //Add possibility to skip this heavy operation for checks
    if(fNRotations>0) DoRotations(ev,lbcandProng,d,HPiAODtrk,fNRotations,isHijing,lb,arrayMC,mcHeader);

    HPiAODtrk=0;
    delete chargedHPi;
    recoArray->Clear();
    delete recoArray;
    if(vtxAODNew){delete vtxAODNew;vtxAODNew=NULL;}
    if(lbcandProng)delete lbcandProng;
  }//loop on pion tracks

  delete LcCand;
  PostData(1,fOutput);

  return;
}

void AliAnalysisTaskSELbtoLcpi4::UserCreateOutputObjects()
{
  // Create the output container
  // Load PID Response
  AliAnalysisManager *man    = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  if (!inputHandler) return;

  fPIDResponse      = inputHandler->GetPIDResponse();
  if(!fPIDResponse) AliFatal("PID response not found.");

  if(fRDCutsProductionLb->GetIsUsePID()){
    fRDCutsProductionLb->GetPidHF()->SetPidResponse(fPIDResponse);
    fRDCutsProductionLb->GetPidpion()->SetPidResponse(fPIDResponse);
    fRDCutsProductionLb->GetPidprot()->SetPidResponse(fPIDResponse);
    fRDCutsProductionLb->GetPidHF()->SetOldPid(kFALSE);
    fRDCutsProductionLb->GetPidpion()->SetOldPid(kFALSE);
    fRDCutsProductionLb->GetPidprot()->SetOldPid(kFALSE);
  }
  if(fRDCutsAnalysisLc->GetIsUsePID()){
    fRDCutsAnalysisLc->GetPidHF()->SetPidResponse(fPIDResponse);
    fRDCutsAnalysisLc->GetPidpion()->SetPidResponse(fPIDResponse);
    fRDCutsAnalysisLc->GetPidprot()->SetPidResponse(fPIDResponse);
    fRDCutsAnalysisLc->GetPidHF()->SetOldPid(kFALSE);
    fRDCutsAnalysisLc->GetPidpion()->SetOldPid(kFALSE);
    fRDCutsAnalysisLc->GetPidprot()->SetOldPid(kFALSE);
  }

  //
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("Histos");

  OpenFile(2);
  fNtupleLambdabUPG = new TNtuple("fNtupleLambdabUPG"," Lb ","massCand:ptLb:pt_Prong0:pt_Prong1:d0_Prong1:d0_Prong0:cosThetaStar:Ct:Prodd0:cosp:cospXY:NormDL:ImpPar:dca:signal:rotated:ptLc:d0_Prong0Lc:d0_Prong1Lc:d0_Prong2Lc:pt_Prong0Lc:pt_Prong1Lc:pt_Prong2Lc:dist12Lc:sigmavertLc:distprimsecLc:costhetapointLc:dcaLc:signalLc:promptLc");
  fNtupleLambdabUPG->SetMaxVirtualSize(1.e+8);
  PostData(2,fNtupleLambdabUPG);

 
  fHistNEvents = new TH1F("fHistNEvents", "number of events ",6,-0.5,5.5);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"n lc");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"n lb");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"n lc mother");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"lb lc3ok pion ");
  fHistNEvents->GetXaxis()->SetBinLabel(6,"lb lc3okacc pion acc");
  fHistNEvents->GetXaxis()->SetNdivisions(1,kFALSE);
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);

  PostData(1,fOutput);

  return;
}

//_____________________________________________
void AliAnalysisTaskSELbtoLcpi4::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEvents"));
  fNtupleLambdabUPG = dynamic_cast<TNtuple*>(GetOutputData(2));
  if (!fOutput) {
    printf("ERROR: fOutput not available\n");
    return;
  }

  return;
}

//------------------------------------------------------------------------
void AliAnalysisTaskSELbtoLcpi4::FillLbHists(AliAODRecoDecayHF2Prong *part,Int_t lb,AliAODMCHeader *mcHeader,TClonesArray* arrayMC, AliAODTrack *pion,AliAODRecoDecayHF3Prong *d, Int_t lc, AliAODEvent *ev, Bool_t IsPromptLc){

  Int_t promptLc=0;
  if (IsPromptLc) promptLc=1;
  Bool_t gen = CheckGenerator(pion,d,mcHeader,arrayMC);

  Double_t massTrueLB = 5.641;
  UInt_t pdgLb[2]={211,4122};
  Double_t massCandLb = part->InvMass(2,pdgLb);
  if(TMath::Abs(massCandLb - massTrueLB)>1.) return;

  Double_t ptCandlb = part->Pt();
  
  //fill ntuple
  Float_t lbVarC[30] = {0};
  lbVarC[0] = massCandLb;
  lbVarC[1] = ptCandlb;
  lbVarC[2] = part->PtProng(0);
  lbVarC[3] = part->PtProng(1);
  lbVarC[4] = part->Getd0Prong(1);
  lbVarC[5] = part->Getd0Prong(0);
  lbVarC[6] = part->CosThetaStar(0,5122,4122,211);
  lbVarC[7] = part->Ct(5122);
  lbVarC[8] = part->Prodd0d0();
  lbVarC[9] = part->CosPointingAngle();
  lbVarC[10] = part->CosPointingAngleXY();
  lbVarC[11] = part->NormalizedDecayLengthXY();
  lbVarC[12] = part->ImpParXY()*10000;
  lbVarC[13] = part->GetDCA();
  lbVarC[14] = lb;
  lbVarC[15] = 1; // rotated
  lbVarC[16] = d->Pt();
  lbVarC[17] = d->Getd0Prong(0);
  lbVarC[18] = d->Getd0Prong(1);
  lbVarC[19] = d->Getd0Prong(2);
  lbVarC[20] = d->PtProng(0);
  lbVarC[21] = d->PtProng(1);
  lbVarC[22] = d->PtProng(2);
  lbVarC[23] = d->GetDist12toPrim();
  lbVarC[24] = d->GetSigmaVert(ev);
  lbVarC[25] = d->GetDist23toPrim();
  lbVarC[26] = d->CosPointingAngle();
  lbVarC[27] = d->GetDCA();
  lbVarC[28] = lc;
  lbVarC[29] = promptLc;

  
  if(lb!=1){
    if(gen){
      if(fFillNtupleBackgroundRotated) {
        fNtupleLambdabUPG->Fill(lbVarC);
        PostData(2,fNtupleLambdabUPG);
      }
    }
  }

  return;

}

//-----------------------------------------------------------------------------
AliAODVertex* AliAnalysisTaskSELbtoLcpi4::ReconstructSecondaryVertex(TObjArray *trkArray,Double_t &dispersion,Bool_t useTRefArray) const {
  // Secondary vertex reconstruction with AliVertexerTracks or AliKFParticle
  //        cout<<" Secondary vertex reconstruction with AliVertexerTracks or AliKFParticle"<<endl;
  AliESDVertex *vertexESD = 0;
  AliAODVertex *vertexAOD = 0;
  //Double_t covmatrix[6];
  // AliVertexerTracks
  AliVertexerTracks *vertexer = new AliVertexerTracks(fBzkG);
  vertexer->SetVtxStart((AliESDVertex*)fvtx1);//primary vertex
  vertexESD = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(trkArray);
  delete vertexer; vertexer=NULL;

  if(!vertexESD) return vertexAOD;

  if(vertexESD->GetNContributors()!=trkArray->GetEntriesFast()) {
    //AliDebug(2,"vertexing failed"); i
    //cout<< " vertex failed " << endl;
    delete vertexESD; vertexESD=NULL;
    return vertexAOD;
  }

  Double_t vertRadius2=vertexESD->GetX()*vertexESD->GetX()+vertexESD->GetY()*vertexESD->GetY();
  if(vertRadius2>8.){//(2.82)^2 radius beam pipe
    //clout<<"  // vertex outside beam pipe, reject candidate to avoid propagation through material "<< endl;
    delete vertexESD; vertexESD=NULL;
    return vertexAOD;
  }
  // convert to AliAODVertex
  //
  Double_t pos[3],cov[6],chi2perNDF;
  for(Int_t a=0;a<3;a++)pos[a]=0.;
  for(Int_t b=0;b<6;b++)cov[b]=0.;
  chi2perNDF=0;
  //
  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix
  chi2perNDF = vertexESD->GetChi2toNDF();
  dispersion = vertexESD->GetDispersion();
  delete vertexESD; vertexESD=NULL;
  Int_t nprongs= (useTRefArray ? 0 : trkArray->GetEntriesFast());
  vertexAOD = new AliAODVertex(pos,cov,chi2perNDF,0x0,-1,AliAODVertex::kUndef,nprongs);
  return vertexAOD;
}

//----------------------------------------------------------------------------
void AliAnalysisTaskSELbtoLcpi4::AddDaughterRefs(AliAODVertex *v, const AliVEvent *event,const TObjArray *trkArray) const
{
  // Add the AOD tracks as daughters of the vertex (TRef)
  Int_t nDg = v->GetNDaughters();
  TObject *dg = 0;
  if(nDg) dg = v->GetDaughter(0);
  if(dg) return; // daughters already added

  Int_t nTrks = trkArray->GetEntriesFast();

  AliExternalTrackParam *track = 0;
  AliAODTrack *aodTrack = 0;
  Int_t id;

  for(Int_t i=0; i<nTrks; i++) {
    track = (AliExternalTrackParam*)trkArray->UncheckedAt(i);
    id = (Int_t)track->GetID();
    // printf("---> %d\n",id);
    if(id<0) continue; // this track is a AliAODRecoDecay
    aodTrack = (AliAODTrack*)event->GetTrack(id); // Johannes check this
    v->AddDaughter(aodTrack);
  }
  return;
}
//-------------------------------------------------------------
Int_t AliAnalysisTaskSELbtoLcpi4::CheckMCpartPIONaf(AliAODTrack *HPiAODtrk, TClonesArray* arrayMC){
  Int_t lab=HPiAODtrk->GetLabel();
  if(lab<0)return -2112;
  AliAODMCParticle *partPion= (AliAODMCParticle*)arrayMC->At(HPiAODtrk->GetLabel());
  if(!partPion)return -2112;
  Int_t pdgsPion=0;
  pdgsPion=partPion->GetPdgCode();
  if(TMath::Abs(pdgsPion)==211){
    Int_t labMPion = partPion->GetMother();
    if(labMPion==-1){
      return -2112;
    }
    if(labMPion>=0){
      AliAODMCParticle *partMPion= (AliAODMCParticle*)arrayMC->At(labMPion);
      if(!partMPion) return -2112;
      if(TMath::Abs(partMPion->GetPdgCode())==5122){
        if(partMPion->GetNDaughters()==2){
          AliAODMCParticle *part0dLb=(AliAODMCParticle*)arrayMC->At(partMPion->GetDaughterLabel(0));
          AliAODMCParticle *part1dLb=(AliAODMCParticle*)arrayMC->At(partMPion->GetDaughterLabel(1));
          if(part0dLb && part1dLb){
            Int_t pdgcode0=TMath::Abs(part0dLb->GetPdgCode());
            Int_t pdgcode1=TMath::Abs(part1dLb->GetPdgCode());
            if((pdgcode0==4122 && pdgcode1==211) || (pdgcode1==4122 && pdgcode0==211)){
              return partMPion->GetLabel();
            }else return -2112;
          }else return -2112;
        }else return -2122;//lb in 2 daugh
      }else return -2122;//lb
    }else return -2122;
  }else return -2112;
}
///_____________________________________________________________
Int_t AliAnalysisTaskSELbtoLcpi4::CheckMCLc(AliAODRecoDecayHF3Prong *d,TClonesArray* arrayMC){
  //return mother labe
  const Int_t pdgdaughtersLc[3]={2212,321,211};
  const Int_t labDpL =  d->MatchToMC(4122,arrayMC,3,&pdgdaughtersLc[0]);
  Int_t pdgsl[3]={0,0,0};
  for(Int_t i=0;i<3;i++){
    AliAODTrack *daughl=(AliAODTrack*)d->GetDaughter(i);
    Int_t labl=daughl->GetLabel();
    if(labl<0) continue;
    AliAODMCParticle *partl= (AliAODMCParticle*)arrayMC->At(labl);
    if(!partl) continue;
    pdgsl[i]=TMath::Abs(partl->GetPdgCode());
  }
  
  //in match to mc lc return 0 is not  lambdac
  if(labDpL<0){// from MC if label particle is positive
    return 999;
  }
  AliAODMCParticle *partLc= (AliAODMCParticle*)arrayMC->At(labDpL);
  if(!partLc)return 999;
  Int_t labMLc = partLc->GetMother();//mother Lc
  if(labMLc<0) {fIsPromptLc=kTRUE; return 999;}
  AliAODMCParticle *partMLc= (AliAODMCParticle*)arrayMC->At(labMLc);//MC mother Lc
  if(!partMLc) return 999;
  Int_t pdgsLb=partMLc->GetPdgCode();
  if(TMath::Abs(pdgsLb)==5122){
    //is mother Lc and is Lb
    //check if the lb decays in lc and pion
    if(partMLc->GetNDaughters()==2){
      AliAODMCParticle *part0dLb=(AliAODMCParticle*)arrayMC->At(partMLc->GetDaughterLabel(0));
      AliAODMCParticle *part1dLb=(AliAODMCParticle*)arrayMC->At(partMLc->GetDaughterLabel(1));
      if(part0dLb && part1dLb){
        Int_t pdgcode0=TMath::Abs(part0dLb->GetPdgCode());
        Int_t pdgcode1=TMath::Abs(part1dLb->GetPdgCode());
        //if lb has 2 daught
        if((pdgcode0==4122 && pdgcode1==211) || (pdgcode1==4122 && pdgcode0==211)){
          if(labDpL<0){
            // from MC if label particle is positive
            return 999;//
          }else{
            return labMLc;
          }
        }
      }else return 999;//daughters exist
    }else return 999;
  }else return 999;//mother exist
  return 999;
}
//----------------------------------------------------------------
Bool_t AliAnalysisTaskSELbtoLcpi4::CheckGenerator(AliAODTrack *p, AliAODRecoDecayHF3Prong *d, AliAODMCHeader *mcHeader,TClonesArray* arrayMC){
  Bool_t LcNotHijing=IsCandidateInjected(d, mcHeader,arrayMC);
  Bool_t pionNotHijing=IsTrackInjected(p,mcHeader,arrayMC);
  if(!LcNotHijing && !pionNotHijing) return kTRUE;
  else return kFALSE;
}
//__________________________________________________________________________
Int_t AliAnalysisTaskSELbtoLcpi4::IsSelectedLbMY(TObject* obj,Int_t selectionLevel,Int_t lb,Int_t isRot, Bool_t isHijing) const{
  //
  // Apply selection for Lb. The selection in on Lb prongs
  //
  //called for each Lb. Cut Applied for pt bin
  AliAODRecoDecayHF2Prong* dd = (AliAODRecoDecayHF2Prong*)obj;
  if(!dd){
    cout<<"AliAODRecoDecayHF2Prong null"<<endl;
    return 0;
  }

  Double_t massTrueLB = 5.641;
  Double_t massCandLb = 0;
  UInt_t pdgLb[2]={0,0};
  pdgLb[0] = 211;//lambdac
  pdgLb[1] = 4122;//pion


  Int_t iPtBinlb = -1;
  Float_t lbVar[14];
  Float_t lbVarbg[14];
  Double_t ptCandlb = dd->Pt();
  if(ptCandlb>0. && ptCandlb<2.) iPtBinlb=0;
  if(ptCandlb>=2. && ptCandlb<4.) iPtBinlb=1;
  if(ptCandlb>=4. && ptCandlb<7.) iPtBinlb=2;
  if(ptCandlb>=7. && ptCandlb<10.) iPtBinlb=3;
  if(ptCandlb>=10. && ptCandlb<14.) iPtBinlb=4;
  if(ptCandlb>=14.) iPtBinlb=5;

  Float_t cutV[6][22]=
  //0,   1     2      3       4      5    6     7      8     9       10    11    12   13        14      15     16     17         18       19     20     21
  //M, ptp0<, ptp1<, ptp0>, ptp1>, d0p1>,d0p1<,d0p0>,d0p0<,cosTSt<,cosTSt>,ct<, ct>,prodd0d0>,prodd0d0<,cosp>,cospXY<,normDLXY>,normDLXY<,dca>, dca<, ImpParXY>
  /*  {{1., 0.,  0.,     0.,     0.,   0.,     0.,   0.,   0.,   0.,     0.,     0.,  0.,   0.,        0.,   0.,    0.,     0.,       0.,        0., 0.,    0.},
   {1., 0.,  0.,     0.,     0.,   0.,     0.,   0.,   0.,   0.,     0.,     0.,  0.,   0.,        0.,   0.,    0.,     0.,       0.,        0., 0.,    0.},
   {1.,1.1,  3.,     5.,    7.5, 0.04,  0.004, 999.,0.002,  -1.,    0.6,   0.02, 0.175, 0.,      -999, 0.9967,  0.998,  13.6,    1.3, 0.0045,0.00009,  26.},//4-7
   {1.,2. ,  3.,     6.,     8.,  999, 0.0045, 999.,0.0055,-0.6,     1.,  0.015, 9999.,-0.00004, -999,  0.998,  0.999,  13.5,     1.,0.00333,0.00005,18.},//7-10
   {1.,1. , 2.5,    14.,    22., 0.04,  0.002, 0.08,0.003, -0.7,   1.4,  -999,   0.08,-0.03*1e-3,-999, 0.998,  0.998,  14.,      0.2,  9999, -9999,   22.},//10-14
   {1.,2. ,  6.,    14.,    22., 0.04,  0.0005, 0.08,0.001, -0.6,   1. ,  0.0018, 9999., 0.,      -0.003,0.998, 0.9993,  14.,      0.3, 0.015,-9999.,  20.}};//>14
   // {1.,2. ,  6.,    14.,    22., 0.04,  0.0005, 0.08,0.001, -0.6,   1. ,  0.002, 0.08, 0.,      -0.003,0.998, 0.9993,  14.,      0.3, 0.015,-9999.,  20.}};//>14*/
  {{1., 0.,  0.,   999., 999., 999.,-999.,999.,-999.,         -10., 10.,  0.005, 999.,   0.     -999, 0.96,    0.96,        999., 0.5, 0.008,  0., 50.},
    {1., 0.,  0.,  999., 999., 999.,-999.,999.,-999.,         -10., 10.,  0.005, 999.,   0.,    -999, 0.96,    0.96,        999., 0.5, 0.008,  0., 50.},
    {1.,0.,  0.,   999., 999., 999.,-999.,999.,-999.,        -10.,  10.,  0.005, 999.,   0.,    -999, 0.96,    0.96,        999., 0.5, 0.008,  0., 50.},//4-7
    {1.,0. ,  0.,  999., 999., 999.,-999.,999.,-999.,       -10,    10.,  0.005, 999.,   0.,    -999, 0.96,    0.96,        999., 0.5, 0.008,  0., 50.},//7-10
    {1.,0. , 0.,   999., 999., 999.,-999.,999.,-999.,       -10,    10.,  0.005, 999.,   0.,    -999, 0.96,    0.96,        999., 0.5, 0.008,  0., 50.},//10-14
    {1.,0. ,  0.,  999., 999., 999.,-999.,999.,-999.,        -10,   10.,  0.005, 999.,   0.,    -999, 0.96,    0.96,        999., 0.5, 0.008,  0., 50.}};//>14

  massCandLb = dd->InvMass(2,pdgLb);
  if(TMath::Abs(massCandLb - massTrueLB)>cutV[iPtBinlb][0]) return 0;
  else if(dd->PtProng(0) < cutV[iPtBinlb][1]) return 0;
  else if(dd->PtProng(1) < cutV[iPtBinlb][2]) return 0;
  else if(dd->PtProng(0) > cutV[iPtBinlb][3]) return 0;
  else if(dd->PtProng(1) > cutV[iPtBinlb][4]) return 0;
  else if(TMath::Abs(dd->Getd0Prong(1)) > cutV[iPtBinlb][5] || TMath::Abs(dd->Getd0Prong(1)) < cutV[iPtBinlb][6]) return 0;
  else if(TMath::Abs(dd->Getd0Prong(0)) > cutV[iPtBinlb][7] || TMath::Abs(dd->Getd0Prong(0)) < cutV[iPtBinlb][8]) return 0;
  //NOT checked for the moment if(dd->CosThetaStar(0,5122,4122,211)<cutV[iPtBinlb][9] || dd->CosThetaStar(0,5122,4122,211)>cutV[iPtBinlb][10])
  else if(dd->Ct(5122)<cutV[iPtBinlb][11] || dd->Ct(5122)>cutV[iPtBinlb][12]) return 0;
  else if((dd->Prodd0d0()) >cutV[iPtBinlb][13] || dd->Prodd0d0() <cutV[iPtBinlb][14]) return 0;
  else if(dd->CosPointingAngle() <cutV[iPtBinlb][15]) return 0;
  else if(dd->CosPointingAngleXY() < cutV[iPtBinlb][16]) return 0;
  else if(dd->NormalizedDecayLengthXY()>cutV[iPtBinlb][17] || dd->NormalizedDecayLengthXY()<cutV[iPtBinlb][18]) return 0;
  else if(dd->GetDCA() > cutV[iPtBinlb][19] || dd->GetDCA() <cutV[iPtBinlb][20]) return 0;
  else if(TMath::Abs(dd->ImpParXY())*10000.>cutV[iPtBinlb][21]) return 0;
  else return 1;
}

//__________________________________________________________________________________
Bool_t AliAnalysisTaskSELbtoLcpi4::CountLc(AliAODRecoDecayHF3Prong* Lc,AliAODTrack* pion, TClonesArray* arrayMC,Int_t motherLabelLc,Int_t motherLabelpione)
{//if Lc is signal how many time is attached to pion to form background
  const Int_t pdgdaughtersLc[3]={2212,321,211};
  const Int_t labLc =  Lc->MatchToMC(4122,arrayMC,3,&pdgdaughtersLc[0]);
  if(labLc<0) return kFALSE;
  if(motherLabelLc != motherLabelpione)return kTRUE;
  else return kFALSE;
}
//_____________________________________________________________
void AliAnalysisTaskSELbtoLcpi4::DoRotations(AliAODEvent* ev, AliAODRecoDecayHF2Prong *decay, AliAODRecoDecayHF3Prong* lc3prong, AliAODTrack* piontrack, Int_t nRot, Bool_t isHijing, Int_t lb, TClonesArray* arrayMC, AliAODMCHeader *mcHeader) {
  //
  // Do full rotation procedure in one function and delete everything at the end to fix memory leak
  //

  // magnetic field
  Double_t bz=ev->GetMagneticField();

  //primary vertex
  AliVVertex *primaryVertex=ev->GetPrimaryVertex();
  AliAODVertex *primaryVertexAOD=ev->GetPrimaryVertex();
  if(!primaryVertex || !primaryVertexAOD) return;

  Double_t pseudoX2[3], pseudoP2[3];
  Double_t CovPseudo2[21],CovPseudo1[21];
  Double_t pseudoX1[3],pseudoP1[3];

  // positive track to be rotated
  AliAODTrack* positive = (AliAODTrack*)decay->GetDaughter(0);//pion
  positive->GetCovarianceXYZPxPyPz(CovPseudo2);
  positive->GetXYZ(pseudoX2);
  positive->GetPxPyPz(pseudoP2);
  Short_t sign = positive->Charge();
  //AliExternalTrackParam* et1 = new AliExternalTrackParam(pseudoX2,pseudoP2,CovPseudo2,sign); //made below

  // negative track
  AliAODTrack* negative = (AliAODTrack*)decay->GetDaughter(1);//Lc
  negative->GetCovarianceXYZPxPyPz(CovPseudo1);
  negative->GetXYZ(pseudoX1);
  negative->GetPxPyPz(pseudoP1);
  Short_t sign1 = negative->Charge();
  AliExternalTrackParam* et2 = new AliExternalTrackParam(pseudoX1,pseudoP1,CovPseudo1,sign1);
  Double_t Prot[3];

  //Double_t fRot=13.;//20
  Double_t fRot=nRot;
  Double_t fAngle=0.0872;//0.1047
  Double_t fAngleFirst=3.14;
  if(nRot==20){
    fAngle=0.1047;//0.0872;//5 gradi...
  }

  Double_t d0z0[2],covd0z0[3],d0[2],d0err[2];
  Double_t d0z02[2],covd0z02[3];
  Double_t xdummy=0.,ydummy=0.,dca;
  Double_t px[2],py[2],pz[2];
  UShort_t id[2];

  // parameters of the actual 2Prong
  for (Int_t i=0;i<2;++i) {
    const AliAODTrack *t=static_cast<AliAODTrack*>(decay->GetDaughter(i));
    px[i]=t->Px();
    py[i]=t->Py();
    pz[i]=t->Pz();
    d0[i]= decay->Getd0Prong(i);           // momenta of the pion change. It need to be recalculated
    d0err[i]= decay->Getd0errProng(i);     // same as before
    id[i]= decay->GetProngID(i);
  }
  dca = decay->GetDCA(); // will be recalculated

  // rotate the first track in the xy plane
  for (Int_t r = 0; r < fRot; r++)
  {
    AliExternalTrackParam* et1;
    if(nRot==20){
      Prot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r*fAngle)*pseudoP2[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r*fAngle)*pseudoP2[1];
      Prot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r*fAngle)*pseudoP2[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r*fAngle)*pseudoP2[1];
      Prot[2] = pseudoP2[2];
      et1 = new AliExternalTrackParam(pseudoX2,Prot,CovPseudo2,sign);
    } else {
      Prot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r*fAngle)*pseudoP2[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r*fAngle)*pseudoP2[1];
      Prot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r*fAngle)*pseudoP2[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r*fAngle)*pseudoP2[1];
      Prot[2] = pseudoP2[2];
      et1 = new AliExternalTrackParam(pseudoX2,Prot,CovPseudo2,sign);
    }
    TObjArray ta12;
    ta12.Add(et1); ta12.Add(et2);

    //recalculate the secondary vertex
    //not crucial now but needed if you rotate the momenta
    AliAODVertex *vtxt=RecalculateVertex(primaryVertex,&ta12 ,bz);
    if(!vtxt){
      et1->Reset(); delete et1; et1=0x0;
      ta12.Clear(); ta12=0x0;
      continue;
    }

    // this are the new impact prameters
    // with relative errors
    et1->PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
    d0[0]=d0z0[0];
    d0err[0] = TMath::Sqrt(covd0z0[0]);
    et2->PropagateToDCA(primaryVertex,bz,100.,d0z02,covd0z02);
    d0[1]=d0z02[0];
    d0err[1] = TMath::Sqrt(covd0z02[0]);

    //this is the new DCA
    dca=et1->GetDCA(et2,bz,xdummy,ydummy);
    //daughters momenta
    Double_t px1[2],py1[2],pz1[2];
    px1[1]= pseudoP1[0];
    py1[1]= pseudoP1[1];
    pz1[1]= pseudoP1[2];
    px1[0]= Prot[0];
    py1[0]= Prot[1];
    pz1[0]= Prot[2];

    // make the rotated Lb candidate and add to the candidate vector
    // this is the new rotated candidate
    AliAODRecoDecayHF2Prong *the2Prong = new AliAODRecoDecayHF2Prong(vtxt,px1,py1,pz1,d0,d0err,dca);
    the2Prong->SetCharge(decay->Charge());
    the2Prong->GetSecondaryVtx()->AddDaughter(positive);
    the2Prong->GetSecondaryVtx()->AddDaughter(negative);
    the2Prong->SetProngIDs(2,id);

    /*Update 20/09/19: From here on copied from FillHistos, instead of using array that has big memory leak*/

    the2Prong->SetOwnPrimaryVtx(primaryVertexAOD);
   
    Int_t selectionlbR=IsSelectedLbMY(the2Prong,AliRDHFCuts::kCandidate,lb,1,isHijing);//analysis cut Lb --> to be improved
    if(selectionlbR==0){
      et1->Reset(); delete et1; et1=0x0;
      ta12.Clear(); ta12=0x0;
      the2Prong->UnsetOwnPrimaryVtx();
      delete the2Prong;
      delete vtxt;
      continue;
    }

    FillLbHists(the2Prong, lb, mcHeader, arrayMC, piontrack, lc3prong, lb /*is same as lc*/, ev, fIsPromptLc);

    et1->Reset(); delete et1; et1=0x0;
    ta12.Clear(); ta12=0x0;
    the2Prong->UnsetOwnPrimaryVtx();
    delete the2Prong;
    delete vtxt;
  }
  et2->Reset(); delete et2; et2=0x0;
}

//_____________________________________________________________
AliAODVertex* AliAnalysisTaskSELbtoLcpi4::RecalculateVertex(const AliVVertex *primary,TObjArray *tracks,Double_t bField) {
  //
  // Helper function to recalculate a vertex.
  //

  AliESDVertex *vertexESD = 0;
  AliAODVertex *vertexAOD = 0;
  //Double_t covmatrix[6];
  // AliVertexerTracks
  AliVertexerTracks *vertexer = new AliVertexerTracks(bField);
  vertexer->SetVtxStart((AliESDVertex*)primary);//primary vertex
  vertexESD = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(tracks);
  delete vertexer; vertexer=NULL;

  if(!vertexESD) return vertexAOD;

  if(vertexESD->GetNContributors()!=tracks->GetEntriesFast()) {
    delete vertexESD; vertexESD=NULL;
    return vertexAOD;
  }

  Double_t vertRadius2=vertexESD->GetX()*vertexESD->GetX()+vertexESD->GetY()*vertexESD->GetY();
  if(vertRadius2>8.){//(2.82)^2 radius beam pipe
    delete vertexESD; vertexESD=NULL;
    return vertexAOD;
  }
  // convert to AliAODVertex
  //
  Double_t dispersion;
  Double_t pos[3],cov[6],chi2perNDF;
  for(Int_t a=0;a<3;a++)pos[a]=0.;
  for(Int_t b=0;b<6;b++)cov[b]=0.;
  chi2perNDF=0;
  //

  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix
  chi2perNDF = vertexESD->GetChi2toNDF();
  dispersion = vertexESD->GetDispersion();
  delete vertexESD; vertexESD=NULL;
  Int_t nprongs= tracks->GetEntriesFast();
  vertexAOD = new AliAODVertex(pos,cov,chi2perNDF,0x0,-1,AliAODVertex::kUndef,nprongs);
  return vertexAOD;

}
//___________________________
void AliAnalysisTaskSELbtoLcpi4::FillLbHistsnr(AliAODRecoDecayHF2Prong *part,Int_t lb,AliAODMCHeader *mcHeader,TClonesArray* arrayMC, AliAODTrack *pion,AliAODRecoDecayHF3Prong *d, Int_t lc,AliAODEvent *ev, Bool_t IsPromptLc){
  //ptlb cut

  Int_t promptLc=0;
  if (IsPromptLc) promptLc=1;
  Bool_t gen = CheckGenerator(pion,d,mcHeader,arrayMC);

  Double_t massTrueLB = 5.641;
  UInt_t pdgLb[2]={211,4122};
  Double_t massCandLb = part->InvMass(2,pdgLb);
  if(TMath::Abs(massCandLb - massTrueLB)>1.) return;

  Float_t lbVar[14];
  Float_t lbVarbg[14];
  Double_t ptCandlb = part->Pt();
 
  //fill ntuple
  Float_t lbVarC[30] = {0};
  lbVarC[0] = massCandLb;
  lbVarC[1] = ptCandlb;
  lbVarC[2] = part->PtProng(0);
  lbVarC[3] = part->PtProng(1);
  lbVarC[4] = part->Getd0Prong(1);
  lbVarC[5] = part->Getd0Prong(0);
  lbVarC[6] = part->CosThetaStar(0,5122,4122,211);
  lbVarC[7] = part->Ct(5122);
  lbVarC[8] = part->Prodd0d0();
  lbVarC[9] = part->CosPointingAngle();
  lbVarC[10] = part->CosPointingAngleXY();
  lbVarC[11] = part->NormalizedDecayLengthXY();
  lbVarC[12] = part->ImpParXY()*10000;
  lbVarC[13] = part->GetDCA();
  lbVarC[14] = lb;
  lbVarC[15] = 0; // not rotated
  lbVarC[16] = d->Pt();
  lbVarC[17] = d->Getd0Prong(0);
  lbVarC[18] = d->Getd0Prong(1);
  lbVarC[19] = d->Getd0Prong(2);
  lbVarC[20] = d->PtProng(0);
  lbVarC[21] = d->PtProng(1);
  lbVarC[22] = d->PtProng(2);
  lbVarC[23] = d->GetDist12toPrim();
  lbVarC[24] = d->GetSigmaVert(ev);
  lbVarC[25] = d->GetDist23toPrim();
  lbVarC[26] = d->CosPointingAngle();
  lbVarC[27] = d->GetDCA();
  lbVarC[28] = lc;
  lbVarC[29] = promptLc;

  if(lb==1){
    if(fFillNtupleSignal) {
      fNtupleLambdabUPG->Fill(lbVarC);
      PostData(2,fNtupleLambdabUPG);
    }
  }
  if(lb!=1){
    if(gen){
      if(fFillNtupleBackgroundNonRotated) {
        fNtupleLambdabUPG->Fill(lbVarC);
        PostData(2,fNtupleLambdabUPG);
      }
    }
  }
  return;
}
//___________________________________________________________
Int_t AliAnalysisTaskSELbtoLcpi4::IsTrackInjected(AliAODTrack *part,AliAODMCHeader *header,TClonesArray *arrayMC){

  AliVertexingHFUtils* ggg=new  AliVertexingHFUtils();

  Int_t lab=-999.;
  if(fApplyFixesITS3AnalysisHijing)lab=TMath::Abs(part->GetLabel());
  else                             lab=part->GetLabel();
  if(lab<0) {delete ggg;return 1;} //
  TString nameGen=ggg->GetGenerator(lab,header);
  TString empty="";
  Int_t countControl =0;
  while(nameGen.IsWhitespace()){
    AliAODMCParticle *mcpart= (AliAODMCParticle*)arrayMC->At(lab);
    if(!mcpart){
      printf("AliVertexingHFUtils::IsTrackInjected - BREAK: No valid AliAODMCParticle at label %i\n",lab);
      break;
    }
    Int_t mother = mcpart->GetMother();
    if(mother<0){
      printf("AliVertexingHFUtils::IsTrackInjected - BREAK: Reached primary particle without valid mother\n");
      break;
    }
    lab=mother;
    nameGen=ggg->GetGenerator(mother,header);
    countControl++;
    if(countControl>=10){ // 10 = arbitrary number; protection from infinite loops
      printf("AliVertexingHFUtils::IsTrackInjected - BREAK: Protection from infinite loop active\n");
      break;
    }
  }
  if(nameGen.IsWhitespace() || nameGen.Contains("ijing")){delete ggg; return 0;}

  delete ggg;
  return 1;
}

//_____________________________________________________________
Bool_t AliAnalysisTaskSELbtoLcpi4::IsCandidateInjected(AliAODRecoDecayHF *part, AliAODMCHeader *header,TClonesArray *arrayMC){

  Int_t nprongs=part->GetNProngs();
  for(Int_t i=0;i<nprongs;i++){
    AliAODTrack *daugh=(AliAODTrack*)part->GetDaughter(i);
    Int_t lab=-999.;
    if(fApplyFixesITS3AnalysisHijing)lab=TMath::Abs(daugh->GetLabel());
    else                             lab=daugh->GetLabel();
    if(lab<0) return 0;
    if(IsTrackInjected(daugh,header,arrayMC)) return kTRUE;
  }
  return kFALSE;
}
        
