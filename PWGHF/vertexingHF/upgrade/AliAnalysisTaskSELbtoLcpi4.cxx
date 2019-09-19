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
  fHistNEventsCuts(0),
  fHistNEventsCutsLb(0),
  fRDCutsAnalysisLc(0),
  fRDCutsProductionLb(0),
  fListCuts(0),
  fBzkG(0.),
  fvtx1(0x0),
  fInvMassLbSign0(0),
  fInvMassLbSign1(0),
  fInvMassLbSign2(0),
  fInvMassLbSign3(0),
  fInvMassLbSign4(0),
  fInvMassLbSign5(0),
  fSelMC(0),
  fCountLc(0),
  fNtupleLambdabUPG(0),
  //fNtupleDiffD0rot(0)
  fFillNtupleSignal(kFALSE), 
  fFillNtupleBackgroundRotated(kFALSE),
  fFillNtupleBackgroundNonRotated(kFALSE),
  fCutsond0Lcdaughters(kFALSE),
  fIsPromptLc(kFALSE),
  fApplyFixesITS3AnalysisBit(kFALSE),
  fApplyFixesITS3AnalysiskAll(kFALSE),
  fApplyFixesITS3AnalysisHijing(kFALSE)
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
  fHistNEventsCuts(0),
  fHistNEventsCutsLb(0),
  fRDCutsAnalysisLc(lccutsana),
  fRDCutsProductionLb(lccutsprod),
  fListCuts(0),
  fBzkG(0.),
  fvtx1(0x0),
  fInvMassLbSign0(0),
  fInvMassLbSign1(0),
  fInvMassLbSign2(0),
  fInvMassLbSign3(0),
  fInvMassLbSign4(0),
  fInvMassLbSign5(0),
  fSelMC(0),
  fCountLc(0),
  fNtupleLambdabUPG(0),
  //fNtupleDiffD0rot(0)
  fFillNtupleSignal(kFALSE), 
  fFillNtupleBackgroundRotated(kFALSE),
  fFillNtupleBackgroundNonRotated(kFALSE),
  fCutsond0Lcdaughters(kFALSE),
  fIsPromptLc(kFALSE),
  fApplyFixesITS3AnalysisBit(kFALSE),
  fApplyFixesITS3AnalysiskAll(kFALSE),
  fApplyFixesITS3AnalysisHijing(kFALSE)
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
  //DefineOutput(3,TNtuple::Class());
  //DefineOutput(4,TNtuple::Class());
  //DefineOutput(5,TNtuple::Class());
  //DefineOutput(6,TNtuple::Class());
  //DefineOutput(7,TNtuple::Class());


}

AliAnalysisTaskSELbtoLcpi4::~AliAnalysisTaskSELbtoLcpi4() {
  //
  // Destructor.
  //

  if (fPIDResponse) {
    delete  fPIDResponse;
  }
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
  if(fHistNEvents){
    delete fHistNEvents;
    fHistNEvents=0;
  }
  if(fHistNEventsCuts){
    delete fHistNEventsCuts;
    fHistNEventsCuts=0;
  }
  if(fHistNEventsCutsLb){
    delete fHistNEventsCutsLb;
    fHistNEventsCutsLb=0;
  }


  if(fRDCutsAnalysisLc){
    delete fRDCutsAnalysisLc;
    fRDCutsAnalysisLc = 0;
  }

  if(fRDCutsProductionLb){
    delete fRDCutsProductionLb;
    fRDCutsProductionLb = 0;
  }
  if (fListCuts) {
    delete fListCuts;
    fListCuts = 0;
  }
  if(fvtx1){
    delete fvtx1; 
    fvtx1 = 0x0;
  }
 
 delete fNtupleLambdabUPG;

}
//-----------------------------------------------------
void AliAnalysisTaskSELbtoLcpi4::Init()
{
  // Initialization


  fListCuts=new TList();

  fListCuts->Add(new AliRDHFCutsLctopKpi(*fRDCutsAnalysisLc));
  fListCuts->Add(new AliRDHFCutsLctopKpi(*fRDCutsProductionLb));


  //PostData(7,fListCuts);
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

  if(!aod) return;
  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField())<0.001) return;

  fBzkG = aod->GetMagneticField();
  // AOD primary vertex
  fvtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  if(!fvtx1) return;
  if(array3Prong==NULL) return;
  Int_t n3Prong = array3Prong->GetEntriesFast();
  //if (array3Prong)std::cout<<"the array is there"<<std::endl;
  //if(!array3Prong==NULL)
  //Int_t n3Prong = array3Prong->GetEntriesFast();
 // else{return;}

  AliAODMCHeader *mcHeader = 0;
  mcHeader = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  if(!mcHeader) {
    printf(" MC header branch not found!\n");
    return;
  }
  TClonesArray *mcs=static_cast<TClonesArray*>(aod->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
  if (!mcs) return;

  //CheckMCKine(mcs);

  //  loop 3 prongs 
  for (Int_t icand = 0; icand < n3Prong; icand++) {
    AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(icand);

    if(fApplyFixesITS3AnalysisBit){
      if(!(d->HasSelectionBit(AliRDHFCuts::kLcCuts))) continue;
    }
    
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()){
      d->SetOwnPrimaryVtx(fvtx1);
      unsetvtx=kTRUE;
    }
    Int_t selection;
    if(fApplyFixesITS3AnalysiskAll)selection=fRDCutsAnalysisLc->IsSelected(d,AliRDHFCuts::kAll,aod);//is selected for LambdaC
    else                           selection=fRDCutsAnalysisLc->IsSelected(d,AliRDHFCuts::kCandidate,aod);//is selected for LambdaC
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
   if (fCutsond0Lcdaughters)
   {   
      if(TMath::Abs(d->Getd0Prong(0))<fCutD0Daughter[0] || (TMath::Abs(d->Getd0Prong(2))<fCutD0Daughter[1])){
        if(unsetvtx) d->UnsetOwnPrimaryVtx();
        continue;
      }
   }
     // if(TMath::Abs(d->Getd0Prong(0))<0.002 || (TMath::Abs(d->Getd0Prong(2))<0.002))continue;
    //Dist12 and Dist23 and DecayLength on Lc are hardcoded at 0.5
    //     if(d->GetDist12toPrim()>1.) continue;
    //     if(d->GetDist23toPrim()>1.) continue;
    //     if(d->DecayLength()>0.6) continue;
    FillHistos(d,mcs,aod,mcHeader);
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  }
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
    //______________________________________________________________________________________
    // Track cuts
    //cout << "Applying Track Cuts" << endl;
    AliESDtrackCuts *trackCutsHPi = fRDCutsAnalysisLc->GetTrackCuts();//
    if(!fRDCutsAnalysisLc->IsDaughterSelected(HPiAODtrk,(AliESDVertex*)fvtx1,trackCutsHPi)){
      //cout << " *** Track Rejected *** " << endl;
      HPiAODtrk=0;
      continue;
    }
    //AliExternalTrackParam *chargedHPi_old = new AliExternalTrackParam(HPiAODtrk);
    AliExternalTrackParam *chargedHPi = new AliExternalTrackParam;
    chargedHPi->CopyFromVTrack(HPiAODtrk);
    //______________________________________________________________________________________
    // implement check for pion here later
    // Cut on impact parameter of Pion track candidate w.r.t. the primary vertex

    //              Double_t dAtDCA = chargedHPi->GetDCA();
    //               ((TH1F*)fOutput->FindObject("fDCApionBg"))->Fill(dAtDCA);
    Double_t dAtDCALc = d->GetDCA();
      //keep this since we know we have bg above this number
    if(dAtDCALc>0.05){
      HPiAODtrk=0;
      delete chargedHPi;
      continue;
    }
      
       //out for the large pt cuts
    ((TH1F*)fOutput->FindObject("fDCALcBg"))->Fill(dAtDCALc);
    //further cuts on candidate charged track
    const Double_t max = 1;
    Double_t d0cut[2],covd0cut[3];
    Double_t d0cutLc[2],covd0cutLc[3];
    chargedHPi->PropagateToDCA(fvtx1,fBzkG,max,d0cut,covd0cut);
    LcCand->PropagateToDCA(fvtx1,fBzkG,max,d0cutLc,covd0cutLc);
  /*    if(TMath::Abs(d0cutLc[0])<0.002||TMath::Abs(d0cutLc[0])>0.04||TMath::Abs(d0cut[0])<0.003||TMath::Abs(d0cut[0])>0.08) {
      HPiAODtrk=0;
      delete chargedHPi;
      continue;
   }
      */
 
    // Construction of Secondary vertex
    //cout << "Building Secondary Vertex" << endl;

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
    if(vtxAODNew) {
      AddDaughterRefs(vtxAODNew,(AliAODEvent*)ev,recoArray);
    }
    //______________________________________________________________________________________
    // construction of lb (with secondary vertex)
    //cout << "Constructing lb Candidate" << endl;

    const Double_t maxd = 1;
    //___
    // Propagate candidates to secondary vertex
    //cout << "Propagating Daughter Tracks to Secondary Vertex" << endl;

    Double_t px[2],py[2],pz[2],d0[2],d0err[2],dcaCand;

    chargedHPi->PropagateToDCA(vtxAODNew,fBzkG,maxd,dzdummy,covardummy);
    LcCand->PropagateToDCA(vtxAODNew,fBzkG,maxd,dzdummy,covardummy);
    // Calculate momenta
    //cout << "Calculating Momenta" << endl;
    Double_t momentum[3];
    chargedHPi->GetPxPyPz(momentum);
    px[0] = momentum[0]; py[0] = momentum[1]; pz[0] = momentum[2];
    LcCand->GetPxPyPz(momentum);
    px[1] = momentum[0]; py[1] = momentum[1]; pz[1] = momentum[2];
    // Calculate impact parameters
    //cout << "Calculating Impact Parameters" << endl;
    Double_t d0z0[2],covd0z0[3];
    LcCand->PropagateToDCA(fvtx1,fBzkG,maxd,d0z0,covd0z0);
    d0[1] = d0z0[0];
    d0err[1] = TMath::Sqrt(covd0z0[0]);
    chargedHPi->PropagateToDCA(fvtx1,fBzkG,maxd,d0z0,covd0z0);
    d0[0] = d0z0[0];
    d0err[0] = TMath::Sqrt(covd0z0[0]);

    // Create AliExternalTrackParameter to calculate DCA between both tracks

    dcaCand = chargedHPi->GetDCA(LcCand,fBzkG,xdummy,ydummy);
    //cout << "Calculating DCA between the Tracks dca: " <<dcaCand<< endl;

    // Create lbcandidate as AliAODRecoDecayHF2Prong
    AliAODRecoDecayHF2Prong *lbcandProng = new AliAODRecoDecayHF2Prong(vtxAODNew,px,py,pz,d0,d0err,dcaCand);
    //if(lbcandProng->Pt()<4.){//FOR THE MOMENT ONLY CANDIDATES WITH pt>4GeV/c //old configuration
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
      Float_t DCALc=d->GetDCA();
      ((TH1F*)fOutput->FindObject("fDCALc"))->Fill(DCALc);
      fSelMC->Fill(8);
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
    Double_t pionPt=lbcandProng->PtProng(0);
    Double_t pionP=lbcandProng->PProng(0);
    Double_t LcPt=lbcandProng->PtProng(1);
    ((TH1F*)fOutput->FindObject("fpionPt"))->Fill(pionPt);
    ((TH1F*)fOutput->FindObject("fLcPt"))->Fill(LcPt);
    ((TH1F*)fOutput->FindObject("fpionP"))->Fill(pionP);
    //fill not rotated
    //CRI
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
   
    UInt_t pdgLb[2]={0,0};
    pdgLb[1] = 4122;
    pdgLb[0] = 211;
    Double_t massTrueLB = 5.641;
    UInt_t pdgLb2[2]={0,0};
    pdgLb2[1] = 4122;
    pdgLb2[0] = 211;
    Double_t    massLb = lbcandProng->InvMass(2,pdgLb);
    if(TMath::Abs(massLb - massTrueLB)<1.) ((TH1F*)fOutput->FindObject("fMassLbbkg"))->Fill(massLb);
    
    //Add possibility to skip this heavy operation for checks
    if(fNRotations>0){
      //Int_t nRot=13.;
      //if(lbcandProng->Pt()>10.)nRot=20.;
      TObjArray *tob=GetArrayCandRotated(ev,lbcandProng,arrayMC,fNRotations);

      Int_t candidates = tob->GetEntriesFast();

      for(Int_t ic = 0; ic < candidates; ic++) {
        AliAODRecoDecayHF2Prong *lb2 =(AliAODRecoDecayHF2Prong*)tob->At(ic);
        if(!lb2){
          delete lb2;
          continue;
        }
        lb2->SetOwnPrimaryVtx(fvtx1);

        Double_t pionPt2=lb2->PtProng(0);
        Double_t pionP2=lb2->PProng(0);
        Double_t LcPt2=lb2->PtProng(1);

        ((TH1F*)fOutput->FindObject("fpionPt2"))->Fill(pionPt2);
        ((TH1F*)fOutput->FindObject("fLcPt2"))->Fill(LcPt2);
        ((TH1F*)fOutput->FindObject("fpionP2"))->Fill(pionP2);
        //   if(lb2->PtProng(1)<0.3){
        //cout << " *** Track Rejected *** " << endl;
        //   lb2->UnsetOwnPrimaryVtx();
        //   delete lb2;
        //   continue;
        // }
        Double_t    massLb2 = lb2->InvMass(2,pdgLb2);
        if(TMath::Abs(massLb2 - massTrueLB)<1.)((TH1F*)fOutput->FindObject("fMassLb2bkg"))->Fill(massLb2);
        Int_t selectionlbR=IsSelectedLbMY(lb2,AliRDHFCuts::kCandidate,lb,1,isHijing);//analysis cut Lb --> to be improved
        if(selectionlbR==0){
          lb2->UnsetOwnPrimaryVtx();
          delete lb2;
          continue;
        }


        Int_t dgLabels[3];
        for(Int_t i=0; i<3; i++) {
          AliAODTrack *trk = (AliAODTrack*)d->GetDaughter(i);
           dgLabels[i] = trk->GetLabel();
        }
        Int_t LabelPion= HPiAODtrk->GetLabel();
        // if(CountLc(d,HPiAODtrk,arrayMC,labPi2,labLb))countLc++;
        FillLbHists(lb2,lb,mcHeader,arrayMC,HPiAODtrk,d, lc, ev, fIsPromptLc);

        //cout << "__________________________________________*Done*__________________________________________" << endl;
        lb2->UnsetOwnPrimaryVtx();
        delete lb2;
      }//loop on candidates lb2
    
      if(lbcandProng){
        lbcandProng->UnsetOwnPrimaryVtx();
      }
    }
    
    HPiAODtrk=0;
    delete chargedHPi;
    recoArray->Clear();
    delete recoArray;
    if(vtxAODNew){delete vtxAODNew;vtxAODNew=NULL;}
    if(lbcandProng)delete lbcandProng;
  }//loop on pion tracks
  // if(countLc>0) fCountLc->Fill(countLc); 
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

  fInvMassLbSign0 = new TH1F("fMassLbSign0", "Lb signal invariant mass 0 ; M [GeV]; Entries",20000,5.641-1.,5.641+1.);
  fInvMassLbSign1 = new TH1F("fMassLbSign1", "Lb signal invariant mass 1; M [GeV]; Entries",20000,5.641-1.,5.641+1.);
  fInvMassLbSign2 = new TH1F("fMassLbSign2", "Lb signal invariant mass 2; M [GeV]; Entries",20000,5.641-1.,5.641+1.);
  fInvMassLbSign3 = new TH1F("fMassLbSign3", "Lb signal invariant mass 3; M [GeV]; Entries",20000,5.641-1.,5.641+1.);
  fInvMassLbSign4 = new TH1F("fMassLbSign4", "Lb signal invariant mass 4; M [GeV]; Entries",20000,5.641-1.,5.641+1.);
  fInvMassLbSign5 = new TH1F("fMassLbSign5", "Lb signal invariant mass 5; M [GeV]; Entries",20000,5.641-1.,5.641+1.);
  fOutput->Add(fInvMassLbSign0);
  fOutput->Add(fInvMassLbSign1);
  fOutput->Add(fInvMassLbSign2);
  fOutput->Add(fInvMassLbSign3);
  fOutput->Add(fInvMassLbSign4);
  fOutput->Add(fInvMassLbSign5);

    fNtupleLambdabUPG = new TNtuple("fNtupleLambdabUPG"," Lb ","massCand:ptLb:pt_Prong0:pt_Prong1:d0_Prong1:d0_Prong0:cosThetaStar:Ct:Prodd0:cosp:cospXY:NormDL:ImpPar:dca:signal:rotated:ptLc:d0_Prong0Lc:d0_Prong1Lc:d0_Prong2Lc:pt_Prong0Lc:pt_Prong1Lc:pt_Prong2Lc:dist12Lc:sigmavertLc:distprimsecLc:costhetapointLc:dcaLc:signalLc:promptLc");
  PostData(2,fNtupleLambdabUPG);

  /*fNtupleLambdacUPG = new TNtuple("fNtupleLambdacUPG"," Lc ","ptLc:d0_Prong0:d0_Prong1:d0_Prong2:pt_Prong0:pt_Prong1:pt_Prong2:dist12:sigmavert:distprimsec:costhetapoint:dca:signal");
  PostData(3,fNtupleLambdacUPG);
   */
//  fNtupleDiffD0rot = new TNtuple("fNtupleDiffD0rot"," diff d0 rot vs angle ","d0_1:d0_2:d0_3:d0_4:d0_5:d0_6:d0_7:d0_8:d0_9:d0_10:d0_11:d0_12:d0_13");
//  PostData(7,fNtupleDiffD0rot);


  fSelMC = new TH1I("trackSelMC", "SelMC",9,-0.5,8.5);
  fSelMC->GetXaxis()->SetBinLabel(1,"is lc from MTMC");
  fSelMC->GetXaxis()->SetBinLabel(2,"sel lab is negative ") ;
  fSelMC->GetXaxis()->SetBinLabel(3,"sel on daugh exis");
  fSelMC->GetXaxis()->SetBinLabel(4,"remaining lc");
  fSelMC->GetXaxis()->SetBinLabel(5,"is secondary Lc ");
  fSelMC->GetXaxis()->SetBinLabel(6,"mother is lb");
  fSelMC->GetXaxis()->SetBinLabel(7,"lb in 2 daughters");
  fSelMC->GetXaxis()->SetBinLabel(8,"lb in lc e pions");
  fSelMC->GetXaxis()->SetBinLabel(9," lc e pions same label");
  fOutput->Add(fSelMC);



  fCountLc = new TH1I("countLc for background","count Lc for background",500,0,499);
  fOutput->Add(fCountLc);

  TH1F *fMassUpg_pt0=new TH1F("fMassUpg_pt0","fMassUpg_pt0",100,2.086,2.486);
  TH1F *fMassUpg_pt1=new TH1F("fMassUpg_pt1","fMassUpg_pt1",100,2.086,2.486);
  TH1F *fMassUpg_pt2=new TH1F("fMassUpg_pt2","fMassUpg_pt2",100,2.086,2.486);
  TH1F *fMassUpg_pt3=new TH1F("fMassUpg_pt3","fMassUpg_pt3",100,2.086,2.486);
  TH1F *fMassUpg_pt4=new TH1F("fMassUpg_pt4","fMassUpg_pt4",100,2.086,2.486);
  TH1F *fMassUpg_pt5=new TH1F("fMassUpg_pt5","fMassUpg_pt5",100,2.086,2.486);

  fOutput->Add(fMassUpg_pt0);
  fOutput->Add(fMassUpg_pt1);
  fOutput->Add(fMassUpg_pt2);
  fOutput->Add(fMassUpg_pt3);
  fOutput->Add(fMassUpg_pt4);
  fOutput->Add(fMassUpg_pt5);

  TH1F *fMassUpg_pt0lcb=new TH1F("fMassUpg_pt0lbNR","fMassUpg_pt0lbNR",20000,5.641-1.,5.641+1.);
  TH1F *fMassUpg_pt1lcb=new TH1F("fMassUpg_pt1lbNR","fMassUpg_pt1lbNR",20000,5.641-1.,5.641+1.);
  TH1F *fMassUpg_pt2lcb=new TH1F("fMassUpg_pt2lbNR","fMassUpg_pt2lbNR",20000,5.641-1.,5.641+1.);
  TH1F *fMassUpg_pt3lcb=new TH1F("fMassUpg_pt3lbNR","fMassUpg_pt3lbNR",20000,5.641-1.,5.641+1.);
  TH1F *fMassUpg_pt4lcb=new TH1F("fMassUpg_pt4lbNR","fMassUpg_pt4lbNR",20000,5.641-1.,5.641+1.);
  TH1F *fMassUpg_pt5lcb=new TH1F("fMassUpg_pt5lbNR","fMassUpg_pt5lbNR",20000,5.641-1.,5.641+1.);

  fOutput->Add(fMassUpg_pt0lcb);
  fOutput->Add(fMassUpg_pt1lcb);
  fOutput->Add(fMassUpg_pt2lcb);
  fOutput->Add(fMassUpg_pt3lcb);
  fOutput->Add(fMassUpg_pt4lcb);
  fOutput->Add(fMassUpg_pt5lcb);


  TH1F *fMassUpg_pt0lb=new TH1F("fMassUpg_pt0lb","fMassUpg_pt0lb",20000,5.641-1.,5.641+1.);
  TH1F *fMassUpg_pt1lb=new TH1F("fMassUpg_pt1lb","fMassUpg_pt1lb",20000,5.641-1.,5.641+1.);
  TH1F *fMassUpg_pt2lb=new TH1F("fMassUpg_pt2lb","fMassUpg_pt2lb",20000,5.641-1.,5.641+1.);
  TH1F *fMassUpg_pt3lb=new TH1F("fMassUpg_pt3lb","fMassUpg_pt3lb",20000,5.641-1.,5.641+1.);
  TH1F *fMassUpg_pt4lb=new TH1F("fMassUpg_pt4lb","fMassUpg_pt4lb",20000,5.641-1.,5.641+1.);
  TH1F *fMassUpg_pt5lb=new TH1F("fMassUpg_pt5lb","fMassUpg_pt5lb",20000,5.641-1.,5.641+1.);
  fOutput->Add(fMassUpg_pt0lb);
  fOutput->Add(fMassUpg_pt1lb);
  fOutput->Add(fMassUpg_pt2lb);
  fOutput->Add(fMassUpg_pt3lb);
  fOutput->Add(fMassUpg_pt4lb);
  fOutput->Add(fMassUpg_pt5lb);


  TH1F *fMassUpg_pt0lbbgOnly=new TH1F("fMassUpg_pt0lbbgOnly","fMassUpg_pt0lbbgOnly",20000,4.641,6.641);
  TH1F *fMassUpg_pt1lbbgOnly=new TH1F("fMassUpg_pt1lbbgOnly","fMassUpg_pt1lbbgOnly",20000,4.641,6.641);
  TH1F *fMassUpg_pt2lbbgOnly=new TH1F("fMassUpg_pt2lbbgOnly","fMassUpg_pt2lbbgOnly",20000,4.641,6.641);
  TH1F *fMassUpg_pt3lbbgOnly=new TH1F("fMassUpg_pt3lbbgOnly","fMassUpg_pt3lbbgOnly",20000,4.641,6.641);
  TH1F *fMassUpg_pt4lbbgOnly=new TH1F("fMassUpg_pt4lbbgOnly","fMassUpg_pt4lbbgOnly",20000,4.641,6.641);
  TH1F *fMassUpg_pt5lbbgOnly=new TH1F("fMassUpg_pt5lbbgOnly","fMassUpg_pt5lbbgOnly",20000,4.641,6.641);

  fOutput->Add(fMassUpg_pt0lbbgOnly);
  fOutput->Add(fMassUpg_pt1lbbgOnly);
  fOutput->Add(fMassUpg_pt2lbbgOnly);
  fOutput->Add(fMassUpg_pt3lbbgOnly);
  fOutput->Add(fMassUpg_pt4lbbgOnly);
  fOutput->Add(fMassUpg_pt5lbbgOnly);

  TH1F *fMassUpg_pt0lbbg=new TH1F("fMassUpg_pt0lbbgNR","fMassUpg_pt0lbbgNR",20000,5.641-1.,5.641+1.);
  TH1F *fMassUpg_pt1lbbg=new TH1F("fMassUpg_pt1lbbgNR","fMassUpg_pt1lbbgNR",20000,5.641-1.,5.641+1.);
  TH1F *fMassUpg_pt2lbbg=new TH1F("fMassUpg_pt2lbbgNR","fMassUpg_pt2lbbgNR",20000,5.641-1.,5.641+1.);
  TH1F *fMassUpg_pt3lbbg=new TH1F("fMassUpg_pt3lbbgNR","fMassUpg_pt3lbbgNR",20000,5.641-1.,5.641+1.);
  TH1F *fMassUpg_pt4lbbg=new TH1F("fMassUpg_pt4lbbgNR","fMassUpg_pt4lbbgNR",20000,5.641-1.,5.641+1.);
  TH1F *fMassUpg_pt5lbbg=new TH1F("fMassUpg_pt5lbbgNR","fMassUpg_pt5lbbgNR",20000,5.641-1.,5.641+1.);

  fOutput->Add(fMassUpg_pt0lbbg);
  fOutput->Add(fMassUpg_pt1lbbg);
  fOutput->Add(fMassUpg_pt2lbbg);
  fOutput->Add(fMassUpg_pt3lbbg);
  fOutput->Add(fMassUpg_pt4lbbg);
  fOutput->Add(fMassUpg_pt5lbbg);




  //
  // JJJ histograms for signal and background
  //
  TH1F *fMassBkg_pt0lb_NoCuts=new TH1F("fMassBkg_pt0lb_NoCuts","fMassBkg_pt0lb_NoCuts",2000,5.641-1.,5.641+1.);
  TH1F *fMassBkg_pt1lb_NoCuts=new TH1F("fMassBkg_pt1lb_NoCuts","fMassBkg_pt1lb_NoCuts",2000,5.641-1.,5.641+1.);
  TH1F *fMassBkg_pt2lb_NoCuts=new TH1F("fMassBkg_pt2lb_NoCuts","fMassBkg_pt2lb_NoCuts",2000,5.641-1.,5.641+1.);
  TH1F *fMassBkg_pt3lb_NoCuts=new TH1F("fMassBkg_pt3lb_NoCuts","fMassBkg_pt3lb_NoCuts",2000,5.641-1.,5.641+1.);
  TH1F *fMassBkg_pt4lb_NoCuts=new TH1F("fMassBkg_pt4lb_NoCuts","fMassBkg_pt4lb_NoCuts",2000,5.641-1.,5.641+1.);
  TH1F *fMassBkg_pt5lb_NoCuts=new TH1F("fMassBkg_pt5lb_NoCuts","fMassBkg_pt5lb_NoCuts",2000,5.641-1.,5.641+1.);
  fOutput->Add(fMassBkg_pt0lb_NoCuts);
  fOutput->Add(fMassBkg_pt1lb_NoCuts);
  fOutput->Add(fMassBkg_pt2lb_NoCuts);
  fOutput->Add(fMassBkg_pt3lb_NoCuts);
  fOutput->Add(fMassBkg_pt4lb_NoCuts);
  fOutput->Add(fMassBkg_pt5lb_NoCuts);
  TH1F *fMassBkgRot_pt0lb_NoCuts=new TH1F("fMassBkgRot_pt0lb_NoCuts","fMassBkgRot_pt0lb_NoCuts",2000,5.641-1.,5.641+1.);
  TH1F *fMassBkgRot_pt1lb_NoCuts=new TH1F("fMassBkgRot_pt1lb_NoCuts","fMassBkgRot_pt1lb_NoCuts",2000,5.641-1.,5.641+1.);
  TH1F *fMassBkgRot_pt2lb_NoCuts=new TH1F("fMassBkgRot_pt2lb_NoCuts","fMassBkgRot_pt2lb_NoCuts",2000,5.641-1.,5.641+1.);
  TH1F *fMassBkgRot_pt3lb_NoCuts=new TH1F("fMassBkgRot_pt3lb_NoCuts","fMassBkgRot_pt3lb_NoCuts",2000,5.641-1.,5.641+1.);
  TH1F *fMassBkgRot_pt4lb_NoCuts=new TH1F("fMassBkgRot_pt4lb_NoCuts","fMassBkgRot_pt4lb_NoCuts",2000,5.641-1.,5.641+1.);
  TH1F *fMassBkgRot_pt5lb_NoCuts=new TH1F("fMassBkgRot_pt5lb_NoCuts","fMassBkgRot_pt5lb_NoCuts",2000,5.641-1.,5.641+1.);
  fOutput->Add(fMassBkgRot_pt0lb_NoCuts);
  fOutput->Add(fMassBkgRot_pt1lb_NoCuts);
  fOutput->Add(fMassBkgRot_pt2lb_NoCuts);
  fOutput->Add(fMassBkgRot_pt3lb_NoCuts);
  fOutput->Add(fMassBkgRot_pt4lb_NoCuts);
  fOutput->Add(fMassBkgRot_pt5lb_NoCuts);
  TH1F *fMassSig_pt0lb_NoCuts=new TH1F("fMassSig_pt0lb_NoCuts","fMassSig_pt0lb_NoCuts",2000,5.641-1.,5.641+1.);
  TH1F *fMassSig_pt1lb_NoCuts=new TH1F("fMassSig_pt1lb_NoCuts","fMassSig_pt1lb_NoCuts",2000,5.641-1.,5.641+1.);
  TH1F *fMassSig_pt2lb_NoCuts=new TH1F("fMassSig_pt2lb_NoCuts","fMassSig_pt2lb_NoCuts",2000,5.641-1.,5.641+1.);
  TH1F *fMassSig_pt3lb_NoCuts=new TH1F("fMassSig_pt3lb_NoCuts","fMassSig_pt3lb_NoCuts",2000,5.641-1.,5.641+1.);
  TH1F *fMassSig_pt4lb_NoCuts=new TH1F("fMassSig_pt4lb_NoCuts","fMassSig_pt4lb_NoCuts",2000,5.641-1.,5.641+1.);
  TH1F *fMassSig_pt5lb_NoCuts=new TH1F("fMassSig_pt5lb_NoCuts","fMassSig_pt5lb_NoCuts",2000,5.641-1.,5.641+1.);
  fOutput->Add(fMassSig_pt0lb_NoCuts);
  fOutput->Add(fMassSig_pt1lb_NoCuts);
  fOutput->Add(fMassSig_pt2lb_NoCuts);
  fOutput->Add(fMassSig_pt3lb_NoCuts);
  fOutput->Add(fMassSig_pt4lb_NoCuts);
  fOutput->Add(fMassSig_pt5lb_NoCuts);
//  //BDT histograms
//  // comment out for now
//  TH2F *fMassBkg_pt0lb_BDT=new TH2F("fMassBkg_pt0lb_BDT","fMassBkg_pt0lb_BDT",2000,-1,1,2000,5.641-1.,5.641+1.);
//  TH2F *fMassBkg_pt1lb_BDT=new TH2F("fMassBkg_pt1lb_BDT","fMassBkg_pt1lb_BDT",2000,-1,1,2000,5.641-1.,5.641+1.);
//  TH2F *fMassBkg_pt2lb_BDT=new TH2F("fMassBkg_pt2lb_BDT","fMassBkg_pt2lb_BDT",2000,-1,1,2000,5.641-1.,5.641+1.);
//  TH2F *fMassBkg_pt3lb_BDT=new TH2F("fMassBkg_pt3lb_BDT","fMassBkg_pt3lb_BDT",2000,-1,1,2000,5.641-1.,5.641+1.);
//  TH2F *fMassBkg_pt4lb_BDT=new TH2F("fMassBkg_pt4lb_BDT","fMassBkg_pt4lb_BDT",2000,-1,1,2000,5.641-1.,5.641+1.);
//  TH2F *fMassBkg_pt5lb_BDT=new TH2F("fMassBkg_pt5lb_BDT","fMassBkg_pt5lb_BDT",2000,-1,1,2000,5.641-1.,5.641+1.);
//  fOutput->Add(fMassBkg_pt0lb_BDT);
//  fOutput->Add(fMassBkg_pt1lb_BDT);
//  fOutput->Add(fMassBkg_pt2lb_BDT);
//  fOutput->Add(fMassBkg_pt3lb_BDT);
//  fOutput->Add(fMassBkg_pt4lb_BDT);
//  fOutput->Add(fMassBkg_pt5lb_BDT);
//  TH2F *fMassBkgRot_pt0lb_BDT=new TH2F("fMassBkgRot_pt0lb_BDT","fMassBkgRot_pt0lb_BDT",2000,-1,1,2000,5.641-1.,5.641+1.);
//  TH2F *fMassBkgRot_pt1lb_BDT=new TH2F("fMassBkgRot_pt1lb_BDT","fMassBkgRot_pt1lb_BDT",2000,-1,1,2000,5.641-1.,5.641+1.);
//  TH2F *fMassBkgRot_pt2lb_BDT=new TH2F("fMassBkgRot_pt2lb_BDT","fMassBkgRot_pt2lb_BDT",2000,-1,1,2000,5.641-1.,5.641+1.);
//  TH2F *fMassBkgRot_pt3lb_BDT=new TH2F("fMassBkgRot_pt3lb_BDT","fMassBkgRot_pt3lb_BDT",2000,-1,1,2000,5.641-1.,5.641+1.);
//  TH2F *fMassBkgRot_pt4lb_BDT=new TH2F("fMassBkgRot_pt4lb_BDT","fMassBkgRot_pt4lb_BDT",2000,-1,1,2000,5.641-1.,5.641+1.);
//  TH2F *fMassBkgRot_pt5lb_BDT=new TH2F("fMassBkgRot_pt5lb_BDT","fMassBkgRot_pt5lb_BDT",2000,-1,1,2000,5.641-1.,5.641+1.);
//  fOutput->Add(fMassBkgRot_pt0lb_BDT);
//  fOutput->Add(fMassBkgRot_pt1lb_BDT);
//  fOutput->Add(fMassBkgRot_pt2lb_BDT);
//  fOutput->Add(fMassBkgRot_pt3lb_BDT);
//  fOutput->Add(fMassBkgRot_pt4lb_BDT);
//  fOutput->Add(fMassBkgRot_pt5lb_BDT);
//  TH2F *fMassSig_pt0lb_BDT=new TH2F("fMassSig_pt0lb_BDT","fMassSig_pt0lb_BDT",2000,-1,1,2000,5.641-1.,5.641+1.);
//  TH2F *fMassSig_pt1lb_BDT=new TH2F("fMassSig_pt1lb_BDT","fMassSig_pt1lb_BDT",2000,-1,1,2000,5.641-1.,5.641+1.);
//  TH2F *fMassSig_pt2lb_BDT=new TH2F("fMassSig_pt2lb_BDT","fMassSig_pt2lb_BDT",2000,-1,1,2000,5.641-1.,5.641+1.);
//  TH2F *fMassSig_pt3lb_BDT=new TH2F("fMassSig_pt3lb_BDT","fMassSig_pt3lb_BDT",2000,-1,1,2000,5.641-1.,5.641+1.);
//  TH2F *fMassSig_pt4lb_BDT=new TH2F("fMassSig_pt4lb_BDT","fMassSig_pt4lb_BDT",2000,-1,1,2000,5.641-1.,5.641+1.);
//  TH2F *fMassSig_pt5lb_BDT=new TH2F("fMassSig_pt5lb_BDT","fMassSig_pt5lb_BDT",2000,-1,1,2000,5.641-1.,5.641+1.);
//  fOutput->Add(fMassSig_pt0lb_BDT);
//  fOutput->Add(fMassSig_pt1lb_BDT);
//  fOutput->Add(fMassSig_pt2lb_BDT);
//  fOutput->Add(fMassSig_pt3lb_BDT);
//  fOutput->Add(fMassSig_pt4lb_BDT);
//  fOutput->Add(fMassSig_pt5lb_BDT);



  TH1F *fPtBkg = new TH1F("fPtBkg","fPtBkg",100,0.,100.);
  TH1F *fPtBkg_TC = new TH1F("fPtBkg_TC","fPtBkg_TC",100,0.,100.);
  TH1F *fPtBkgRot = new TH1F("fPtBkgRot","fPtBkgRot",100,0.,100.);
  TH1F *fPtBkgRot_TC = new TH1F("fPtBkgRot_TC","fPtBkgRot_TC",100,0.,100.);
  TH1F *fPtSig = new TH1F("fPtSig","fPtSig",100,0.,100.);
  TH1F *fPtSig_TC = new TH1F("fPtSig_TC","fPtSigTC",100,0.,100.);
  fOutput->Add(fPtBkg);
  fOutput->Add(fPtBkg_TC);
  fOutput->Add(fPtBkgRot);
  fOutput->Add(fPtBkgRot_TC);
  fOutput->Add(fPtSig);
  fOutput->Add(fPtSig_TC);




  TH1F *fMassLbbkg=new TH1F("fMassLbbkg","fMassLbbkg",500,5.641-1.,5.641+1.);
  fOutput->Add(fMassLbbkg);
  TH1F *fMassLb2bkg=new TH1F("fMassLb2bkg","fMassLb2bkg",500,5.641-1.,5.641+1.);
  fOutput->Add(fMassLb2bkg);
  TH1F *fpionPt=new TH1F("fpionPt","fpionPt",100,0.,25.);
  fOutput->Add(fpionPt);
  TH1F *fLcPt=new TH1F("fLcPt","fLcPt",100,0.,25.);
  fOutput->Add(fLcPt);
  TH1F *fpionP=new TH1F("fpionP","fpionP",100,0.,25.);
  fOutput->Add(fpionP);
  TH1F *fpionPt2=new TH1F("fpionPt2","fpionPt2",100,0.,25.);
  fOutput->Add(fpionPt2);
  TH1F *fLcPt2=new TH1F("fLcPt2","fLcPt2",100,0.,25.);
  fOutput->Add(fLcPt2);
  TH1F *fpionP2=new TH1F("fpionP2","fpionP2",100,0.,25.);
  fOutput->Add(fpionP2);

  TH1F *fd0Pion=new TH1F("fd0Pion","fd0Pion",100,-1.,1.);
  TH1F *fd0Lc=new TH1F("fd0Lc","fd0Lc",100,-1.,1.);
  TH1F *fDCApion=new TH1F("fDCApion","fDCApion",100,-1.,1.);
  TH1F *fDCApionBg=new TH1F("fDCApionBg","fDCApionBg",100,-1.,1.);
  TH1F *fDCALc=new TH1F("fDCALc","fDCALc",100,-0.09,0.09);
  TH1F *fDCALcBg=new TH1F("fDCALcBg","fDCALcBg",100,-0.09,0.09);
  TH1F *fd0Lcprong0=new TH1F("fd0Lcprong0","fd0Lcprong0",100,-0.09,0.09);
  TH1F *fd0Lcprong1=new TH1F("fd0Lcprong1","fd0Lcprong1",100,-0.09,0.09);
  TH1F *fd0Lcprong0Bg=new TH1F("fd0Lcprong0Bg","fd0Lcprong0Bg",100,-0.09,0.09);
  TH1F *fd0Lcprong1Bg=new TH1F("fd0Lcprong1Bg","fd0Lcprong1Bg",100,-0.09,0.09);
    TH1F *fd0Lcprong0nr=new TH1F("fd0Lcprong0nr","fd0Lcprong0nr",100,-0.09,0.09);
    TH1F *fd0Lcprong1nr=new TH1F("fd0Lcprong1nr","fd0Lcprong1nr",100,-0.09,0.09);
    TH1F *fd0Lcprong0Bgnr=new TH1F("fd0Lcprong0Bgnr","fd0Lcprong0Bgnr",100,-0.09,0.09);
    TH1F *fd0Lcprong1Bgnr=new TH1F("fd0Lcprong1Bgnr","fd0Lcprong1Bgnr",100,-0.09,0.09);
    
  fOutput->Add(fd0Pion);
  fOutput->Add(fd0Lc);
  fOutput->Add(fDCApion);
  fOutput->Add(fDCApionBg);
  fOutput->Add(fDCALc);
  fOutput->Add(fDCALcBg);
  fOutput->Add(fd0Lcprong0);
  fOutput->Add(fd0Lcprong1);
  fOutput->Add(fd0Lcprong0Bg);
  fOutput->Add(fd0Lcprong1Bg);
    fOutput->Add(fd0Lcprong0nr);
    fOutput->Add(fd0Lcprong1nr);
    fOutput->Add(fd0Lcprong0Bgnr);
    fOutput->Add(fd0Lcprong1Bgnr);

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


  fHistNEventsCuts = new TH1F("fHistNEventsCuts", "pass cuts ",14,-0.5,13.5);
  fHistNEventsCuts->GetXaxis()->SetBinLabel(1,"pt inf prong 0");
  fHistNEventsCuts->GetXaxis()->SetBinLabel(2,"pt inf prong 1");
  fHistNEventsCuts->GetXaxis()->SetBinLabel(3,"pt sup prong 0");
  fHistNEventsCuts->GetXaxis()->SetBinLabel(4,"pt sup prong 1");
  fHistNEventsCuts->GetXaxis()->SetBinLabel(5,"|d0|>value prong 0");
  fHistNEventsCuts->GetXaxis()->SetBinLabel(6,"|d0|>value prong 1");
  fHistNEventsCuts->GetXaxis()->SetBinLabel(7,"cosStar ");
  fHistNEventsCuts->GetXaxis()->SetBinLabel(8,"ct lb ");
  fHistNEventsCuts->GetXaxis()->SetBinLabel(9,"Prodd0d0");
  fHistNEventsCuts->GetXaxis()->SetBinLabel(10,"cosp<0");
  fHistNEventsCuts->GetXaxis()->SetBinLabel(11,"cospXY");
  fHistNEventsCuts->GetXaxis()->SetBinLabel(12,"dca");
  fHistNEventsCuts->GetXaxis()->SetBinLabel(13,"ImpPar");
  fHistNEventsCuts->GetXaxis()->SetBinLabel(14,"normDecL");
  fHistNEventsCuts->GetXaxis()->SetNdivisions(1,kFALSE);
  fHistNEventsCuts->Sumw2();
  fHistNEventsCuts->SetMinimum(0);
  fOutput->Add(fHistNEventsCuts);

  fHistNEventsCutsLb= new TH1F("fHistNEventsCutsLb", "pass cuts Lb ",14,-0.5,13.5);
  fHistNEventsCutsLb->GetXaxis()->SetBinLabel(1,"pt inf prong 0");
  fHistNEventsCutsLb->GetXaxis()->SetBinLabel(2,"pt inf prong 1");
  fHistNEventsCutsLb->GetXaxis()->SetBinLabel(3,"pt sup prong 0");
  fHistNEventsCutsLb->GetXaxis()->SetBinLabel(4,"pt sup prong 1");
  fHistNEventsCutsLb->GetXaxis()->SetBinLabel(5,"|d0|>value prong 0");
  fHistNEventsCutsLb->GetXaxis()->SetBinLabel(6,"|d0|>value prong 1");
  fHistNEventsCutsLb->GetXaxis()->SetBinLabel(7,"cosStar ");
  fHistNEventsCutsLb->GetXaxis()->SetBinLabel(8,"ct lb ");
  fHistNEventsCutsLb->GetXaxis()->SetBinLabel(9,"Prodd0d0");
  fHistNEventsCutsLb->GetXaxis()->SetBinLabel(10,"cosp");
  fHistNEventsCutsLb->GetXaxis()->SetBinLabel(11,"cospXY");
  fHistNEventsCutsLb->GetXaxis()->SetBinLabel(12,"dca");
  fHistNEventsCutsLb->GetXaxis()->SetBinLabel(13,"ImpPar ");
  fHistNEventsCutsLb->GetXaxis()->SetBinLabel(14,"normDecL ");
  fHistNEventsCutsLb->GetXaxis()->SetNdivisions(1,kFALSE);
  fHistNEventsCutsLb->Sumw2();
  fHistNEventsCutsLb->SetMinimum(0);
  fOutput->Add(fHistNEventsCutsLb);
  PostData(1,fOutput);

//  TString inputVariables = "Ptp,PtK,Ptpi,CosP,DecayL,Dist12Min,SigVert,DCAMax";
//  AliInfo(Form("adding variables %s to reader",inputVariables.Data()));
//  TObjArray *tokens = inputVariables.Tokenize(",");
//  tokens->Print();
//  std::vector<std::string> inputNamesVec;
//  for(Int_t i=0; i<tokens->GetEntries(); i++) {
//    TString variable = ((TObjString*)(tokens->At(i)))->String();
//    AliInfo(Form("* * * added %s to vector",variable.Data()));
//    string tmpvar = variable.Data();
//    inputNamesVec.push_back(tmpvar);
//  }
//  fBDTReader[0]= new ReadBDT_pt4to7( inputNamesVec );
//  fBDTReader[1]= new ReadBDT_pt7to10( inputNamesVec );
//  fBDTReader[2]= new ReadBDT_pt10to14( inputNamesVec );
//  fBDTReader[3]= new ReadBDT_pt14to9999( inputNamesVec );
//  AliInfo("* * * Created BDT reader");



  return;

}

//_____________________________________________
void AliAnalysisTaskSELbtoLcpi4::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEvents"));
  fHistNEventsCuts = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEventsCuts"));
  fHistNEventsCutsLb = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEventsCutsLb"));
  fNtupleLambdabUPG = dynamic_cast<TNtuple*>(GetOutputData(2));
  if (!fOutput) {
    printf("ERROR: fOutput not available\n");
    return;
  }

  return;
}
//------------------------------------------------------------------------
void AliAnalysisTaskSELbtoLcpi4::FillLbHists(AliAODRecoDecayHF2Prong *part,Int_t lb,AliAODMCHeader *mcHeader,TClonesArray* arrayMC, AliAODTrack *pion,AliAODRecoDecayHF3Prong *d, Int_t lc, AliAODEvent *ev, Bool_t IsPromptLc){
  //ptlb cut
  Int_t promptLc=0;
  if (IsPromptLc) promptLc=1;
  Bool_t gen = CheckGenerator(pion,d,mcHeader,arrayMC);
  Double_t massTrueLB = 5.641;
  Double_t massCandLb = 0;
  UInt_t pdgLb[2]={0,0};
  //lc always at first place
  pdgLb[1] = 4122;//lambdac
  pdgLb[0] = 211;//pion
  massCandLb = part->InvMass(2,pdgLb);
  if(TMath::Abs(massCandLb - massTrueLB)>1.) return;
  Int_t iPtBinlb = -1;
  Double_t ptCandlb = part->Pt();
  //if(ptCandlb<2.) return; // dont save anything with pt<2 //old configuration
  if(ptCandlb>0. && ptCandlb<2.) iPtBinlb=0;
  if(ptCandlb>=2. && ptCandlb<4.) iPtBinlb=1;
  if(ptCandlb>=4. && ptCandlb<7.) iPtBinlb=2;
  if(ptCandlb>=7. && ptCandlb<10.) iPtBinlb=3;
  if(ptCandlb>10. && ptCandlb<14.) iPtBinlb=4;
  if(ptCandlb>=14.) iPtBinlb=5;

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

  if(lb==1){ //
    if(iPtBinlb==0)((TH1F*)fOutput->FindObject("fMassUpg_pt0lb"))->Fill(massCandLb);
    if(iPtBinlb==1)((TH1F*)fOutput->FindObject("fMassUpg_pt1lb"))->Fill(massCandLb);
    if(iPtBinlb==2)((TH1F*)fOutput->FindObject("fMassUpg_pt2lb"))->Fill(massCandLb);
    if(iPtBinlb==3)((TH1F*)fOutput->FindObject("fMassUpg_pt3lb"))->Fill(massCandLb);
    if(iPtBinlb==4)((TH1F*)fOutput->FindObject("fMassUpg_pt4lb"))->Fill(massCandLb);
    if(iPtBinlb==5)((TH1F*)fOutput->FindObject("fMassUpg_pt5lb"))->Fill(massCandLb);
      
      
    // note - don't fill rotated signal (not needed)
    }
  if(lb!=1){
    if(gen){
      cout << " gen " << gen << endl;
      if(iPtBinlb==0)((TH1F*)fOutput->FindObject("fMassUpg_pt0lbbgOnly"))->Fill(massCandLb);
      if(iPtBinlb==1)((TH1F*)fOutput->FindObject("fMassUpg_pt1lbbgOnly"))->Fill(massCandLb);
      if(iPtBinlb==2)((TH1F*)fOutput->FindObject("fMassUpg_pt2lbbgOnly"))->Fill(massCandLb);
      if(iPtBinlb==3)((TH1F*)fOutput->FindObject("fMassUpg_pt3lbbgOnly"))->Fill(massCandLb);
      if(iPtBinlb==4)((TH1F*)fOutput->FindObject("fMassUpg_pt4lbbgOnly"))->Fill(massCandLb);
      if(iPtBinlb==5)((TH1F*)fOutput->FindObject("fMassUpg_pt5lbbgOnly"))->Fill(massCandLb);
        
        if(fFillNtupleBackgroundRotated) {
        fNtupleLambdabUPG->Fill(lbVarC);
        PostData(2,fNtupleLambdabUPG); 
      }
    }
  }

  return;

}
//
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
    if(labl<0){
      continue;
    }
    AliAODMCParticle *partl= (AliAODMCParticle*)arrayMC->At(labl);
    if(!partl){
      fSelMC->Fill(2);// from MC if label particle is positive
      continue;
    }
    pdgsl[i]=TMath::Abs(partl->GetPdgCode());
  }
  Double_t massl=0.;
  Int_t iPtBinl = -1;
  if(pdgsl[0]==211 && pdgsl[1]==321 && pdgsl[2]==2212) massl=d->InvMassLcpiKp();
  if(pdgsl[0]==2212 && pdgsl[1]==321 && pdgsl[2]==211) massl=d->InvMassLcpKpi();

  Double_t ptCandl = d->Pt();
  
  if(ptCandl<2.) iPtBinl=0;
  if(ptCandl>=2. && ptCandl<4.) iPtBinl=1;
  if(ptCandl>=4. && ptCandl<6.) iPtBinl=2;
  if(ptCandl>6. && ptCandl<14.) iPtBinl=3;
  if(ptCandl>6. && ptCandl<14.) iPtBinl=3;
  if(ptCandl>14.)iPtBinl=4;



  fSelMC->Fill(0);// from MC if label particle is positive
  //in match to mc lc return 0 is not  lambdac
  if(labDpL<0){
    fSelMC->Fill(1);// from MC if label particle is positive
    return 999;
  }
  AliAODMCParticle *partLc= (AliAODMCParticle*)arrayMC->At(labDpL);
  if(!partLc)return 999;
  Int_t labMLc = partLc->GetMother();//mother Lc
  if(labMLc<0) {fIsPromptLc=kTRUE; return 999;}
  fSelMC->Fill(4);
  AliAODMCParticle *partMLc= (AliAODMCParticle*)arrayMC->At(labMLc);//MC mother Lc
  if(!partMLc) return 999;
  Int_t pdgsLb=partMLc->GetPdgCode();
  if(TMath::Abs(pdgsLb)==5122){
    fSelMC->Fill(5);//is mother Lc and is Lb
    //check if the lb decays in lc and pion
    if(partMLc->GetNDaughters()==2){
      AliAODMCParticle *part0dLb=(AliAODMCParticle*)arrayMC->At(partMLc->GetDaughterLabel(0));
      AliAODMCParticle *part1dLb=(AliAODMCParticle*)arrayMC->At(partMLc->GetDaughterLabel(1));
      if(part0dLb && part1dLb){
        Int_t pdgcode0=TMath::Abs(part0dLb->GetPdgCode());
        Int_t pdgcode1=TMath::Abs(part1dLb->GetPdgCode());
        fSelMC->Fill(6);//if lb has 2 daught
        if((pdgcode0==4122 && pdgcode1==211) || (pdgcode1==4122 && pdgcode0==211)){ 
          if(labDpL<0){
            fSelMC->Fill(1);// from MC if label particle is positive
            return 999;//
          }else{
            fSelMC->Fill(7);  
            if(iPtBinl==0)((TH1F*)fOutput->FindObject("fMassUpg_pt0"))->Fill(massl);
            if(iPtBinl==1)((TH1F*)fOutput->FindObject("fMassUpg_pt1"))->Fill(massl);
            if(iPtBinl==2)((TH1F*)fOutput->FindObject("fMassUpg_pt2"))->Fill(massl);
            if(iPtBinl==3)((TH1F*)fOutput->FindObject("fMassUpg_pt3"))->Fill(massl);
            if(iPtBinl==4)((TH1F*)fOutput->FindObject("fMassUpg_pt4"))->Fill(massl);
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
  Bool_t LcNotHijing;
  Bool_t pionNotHijing;
  LcNotHijing=IsCandidateInjected(d, mcHeader,arrayMC); 
  pionNotHijing=IsTrackInjected(p,mcHeader,arrayMC);
  //cout << " LcNotHijing "<< LcNotHijing << " pionNotHijing " << pionNotHijing << endl;
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
  Int_t cut=1;
  //cut value pt dependent
  Double_t massTrueLB = 5.641;
  Double_t massCandLb = 0;
  UInt_t pdgLb[2]={0,0};
  pdgLb[0] = 211;//lambdac
  pdgLb[1] = 4122;//pion

  //
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
  if(TMath::Abs(massCandLb - massTrueLB)>cutV[iPtBinlb][0]){
    cut=0;
  }
  if(dd->PtProng(0) < cutV[iPtBinlb][1]){
    fHistNEventsCuts->Fill(0);
    if(lb==1)fHistNEventsCutsLb->Fill(0);
    cut = 0;
  }
  if(dd->PtProng(1) < cutV[iPtBinlb][2]) {
    fHistNEventsCuts->Fill(1);
    if(lb==1)fHistNEventsCutsLb->Fill(1);
    cut = 0;
  }
  if(dd->PtProng(0) > cutV[iPtBinlb][3]) {
    fHistNEventsCuts->Fill(2);
    if(lb==1)fHistNEventsCutsLb->Fill(2);
    cut = 0;
  }
  if(dd->PtProng(1) > cutV[iPtBinlb][4]) {
    fHistNEventsCuts->Fill(3);
    if(lb==1)fHistNEventsCutsLb->Fill(3);
    cut = 0;}


    if(TMath::Abs(dd->Getd0Prong(1)) > cutV[iPtBinlb][5] || TMath::Abs(dd->Getd0Prong(1)) < cutV[iPtBinlb][6]){ 
      fHistNEventsCuts->Fill(5);
      if(lb==1)fHistNEventsCutsLb->Fill(5);
      cut = 0;
    } 

    if(TMath::Abs(dd->Getd0Prong(0)) > cutV[iPtBinlb][7] || TMath::Abs(dd->Getd0Prong(0)) < cutV[iPtBinlb][8]){
      fHistNEventsCuts->Fill(4);
      if(lb==1)fHistNEventsCutsLb->Fill(4);
      cut = 0;
    }

    if(dd->CosThetaStar(0,5122,4122,211)<cutV[iPtBinlb][9] || dd->CosThetaStar(0,5122,4122,211)>cutV[iPtBinlb][10]){//era -0.6
      //    cout << " OUT OUT ------------------- THETA STAR " << dd->CosThetaStar(0,5122,211,4122)<<endl;
      fHistNEventsCuts->Fill(6);
      if(lb==1)fHistNEventsCutsLb->Fill(6);
      //    cut = 0;
    }
    if(dd->Ct(5122)<cutV[iPtBinlb][11] || dd->Ct(5122)>cutV[iPtBinlb][12]) {//
      fHistNEventsCuts->Fill(7);
      if(lb==1)fHistNEventsCutsLb->Fill(7);
      //    //cout<< " OUTOUT ------------------------ CTAU "<< dd->Ct(5122)<< endl;
      cut = 0;}


      if((dd->Prodd0d0()) >cutV[iPtBinlb][13] || dd->Prodd0d0() <cutV[iPtBinlb][14]){
        fHistNEventsCuts->Fill(8);
        if(lb==1)fHistNEventsCutsLb->Fill(8);
        //cout << " OUT OUT ------------ PRODDO " << dd->Prodd0d0()<< endl;
        cut = 0;
      }
      if(dd->CosPointingAngle() <cutV[iPtBinlb][15]) {
        fHistNEventsCuts->Fill(9);
        if(lb==1)fHistNEventsCutsLb->Fill(9);
        //cout << " OUT OUT  ---------------- COSP " << dd->CosPointingAngle()<<endl;
        cut = 0;}

        if(dd->CosPointingAngleXY() < cutV[iPtBinlb][16]) {
          fHistNEventsCuts->Fill(10);
          if(lb==1)fHistNEventsCutsLb->Fill(10);
          //cout << " OUT OUT ---------------------cosp XY  " <<dd->CosPointingAngleXY()<< endl; 
          cut = 0;}


          if(dd->NormalizedDecayLengthXY()>cutV[iPtBinlb][17] || dd->NormalizedDecayLengthXY()<cutV[iPtBinlb][18]){ 
            fHistNEventsCuts->Fill(13);
            if(lb==1)fHistNEventsCuts->Fill(13);
            cut =0;//era 30
          }


          if(dd->GetDCA() > cutV[iPtBinlb][19] || dd->GetDCA() <cutV[iPtBinlb][20]) {//era 0.01 puo' essere 0.005  0.004
            fHistNEventsCuts->Fill(11);
            if(lb==1)fHistNEventsCutsLb->Fill(11);
            //    cout << " OUT OUT ---------------------cosp XY  " <<dd->CosPointingAngleXY()<< endl;
            cut= 0;}

            if(TMath::Abs(dd->ImpParXY())*10000.>cutV[iPtBinlb][21]) {//era 60
              fHistNEventsCuts->Fill(12);
              if(lb==1)fHistNEventsCutsLb->Fill(12);
              //cout << " OUT OUT ---------------------cosp XY  " <<dd->CosPointingAngleXY()<< endl; 
              cut = 0;
            }


            // fill histograms

            if(lb==0 && isHijing){
              if(isRot==0) { //not rotated bkg
                ((TH1F*)fOutput->FindObject("fPtBkg"))->Fill(ptCandlb);
                if(cut==1)((TH1F*)fOutput->FindObject("fPtBkg_TC"))->Fill(ptCandlb);
                // these histograms are for the background with no cut
                if(iPtBinlb==0) ((TH1F*)fOutput->FindObject("fMassBkg_pt0lb_NoCuts"))->Fill(massCandLb);
                else if(iPtBinlb==1) ((TH1F*)fOutput->FindObject("fMassBkg_pt1lb_NoCuts"))->Fill(massCandLb);
                else if(iPtBinlb==2) ((TH1F*)fOutput->FindObject("fMassBkg_pt2lb_NoCuts"))->Fill(massCandLb);
                else if(iPtBinlb==3) ((TH1F*)fOutput->FindObject("fMassBkg_pt3lb_NoCuts"))->Fill(massCandLb);
                else if(iPtBinlb==4) ((TH1F*)fOutput->FindObject("fMassBkg_pt4lb_NoCuts"))->Fill(massCandLb);
                else if(iPtBinlb==5) ((TH1F*)fOutput->FindObject("fMassBkg_pt5lb_NoCuts"))->Fill(massCandLb);
              }
              else if(isRot==1) { //rotated bkg
                ((TH1F*)fOutput->FindObject("fPtBkgRot"))->Fill(ptCandlb);
                if(cut==1)((TH1F*)fOutput->FindObject("fPtBkgRot_TC"))->Fill(ptCandlb);
                // these histograms are for the rotated background with no cut
                if(iPtBinlb==0) ((TH1F*)fOutput->FindObject("fMassBkgRot_pt0lb_NoCuts"))->Fill(massCandLb);
                else if(iPtBinlb==1) ((TH1F*)fOutput->FindObject("fMassBkgRot_pt1lb_NoCuts"))->Fill(massCandLb);
                else if(iPtBinlb==2) ((TH1F*)fOutput->FindObject("fMassBkgRot_pt2lb_NoCuts"))->Fill(massCandLb);
                else if(iPtBinlb==3) ((TH1F*)fOutput->FindObject("fMassBkgRot_pt3lb_NoCuts"))->Fill(massCandLb);
                else if(iPtBinlb==4) ((TH1F*)fOutput->FindObject("fMassBkgRot_pt4lb_NoCuts"))->Fill(massCandLb);
                else if(iPtBinlb==5) ((TH1F*)fOutput->FindObject("fMassBkgRot_pt5lb_NoCuts"))->Fill(massCandLb);
              }
            }
            else if(lb==1 && isRot==0){
              ((TH1F*)fOutput->FindObject("fPtSig"))->Fill(ptCandlb);
              if(cut==1)((TH1F*)fOutput->FindObject("fPtSig_TC"))->Fill(ptCandlb);
                // these histograms are for the signal with no cut
              if(iPtBinlb==0) ((TH1F*)fOutput->FindObject("fMassSig_pt0lb_NoCuts"))->Fill(massCandLb);
              else if(iPtBinlb==1) ((TH1F*)fOutput->FindObject("fMassSig_pt1lb_NoCuts"))->Fill(massCandLb);
              else if(iPtBinlb==2) ((TH1F*)fOutput->FindObject("fMassSig_pt2lb_NoCuts"))->Fill(massCandLb);
              else if(iPtBinlb==3) ((TH1F*)fOutput->FindObject("fMassSig_pt3lb_NoCuts"))->Fill(massCandLb);
              else if(iPtBinlb==4) ((TH1F*)fOutput->FindObject("fMassSig_pt4lb_NoCuts"))->Fill(massCandLb);
              else if(iPtBinlb==5) ((TH1F*)fOutput->FindObject("fMassSig_pt5lb_NoCuts"))->Fill(massCandLb);
            }



            if(cut==0)return 0;
            else return 1;//returnvalue;
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
TObjArray* AliAnalysisTaskSELbtoLcpi4::GetArrayCandRotated(AliAODEvent* ev,AliAODRecoDecayHF2Prong *decay,TClonesArray* arrayMC,Int_t nRot) {
  //
  // fill a TObjArray with the rotated candidates
  //
  Bool_t rotateFirst = kTRUE;//pion
  Bool_t rotateSecond = kFALSE;//Lc
  // magnetic field
  Double_t bz=ev->GetMagneticField();
  //primary vertex
  AliVVertex *primaryVertex=ev->GetPrimaryVertex();
  if(!primaryVertex) return 0x0;
  Double_t pseudoX2[3], pseudoP2[3];
  Double_t CovPseudo2[21],CovPseudo1[21];
  Double_t pseudoX1[3],pseudoP1[3];
  // for rotations
  AliExternalTrackParam * et1;//
  AliExternalTrackParam * et2;//
  // positive track to be rotated 
  AliAODTrack* positive = (AliAODTrack*)decay->GetDaughter(0);//pion
  positive->GetCovarianceXYZPxPyPz(CovPseudo2);
  positive->GetXYZ(pseudoX2);
  /*  TRandom *xg = new TRandom(0);
      TRandom *yg = new TRandom(0);
      TRandom *zg = new TRandom(0);
      Double_t x[3];
      x[0]=xg->Gaus(pseudoX2[0],TMath::Sqrt(TMath::Abs(CovPseudo2[3])));
      x[1]=yg->Gaus(pseudoX2[1],TMath::Sqrt(TMath::Abs(CovPseudo2[4])));
      x[2]=zg->Gaus(pseudoX2[2],TMath::Sqrt(TMath::Abs(CovPseudo2[5])));
      positive->GetXYZ(x);
  //cout << " x   " << x << " y   " <<y << " z   " << z << endl;
  //cout << " xor " << pseudoX2[0] << " y or " <<pseudoX2[1] << " z or " << pseudoX2[2] << endl;
  */
  positive->GetPxPyPz(pseudoP2);
  Short_t sign = positive->Charge();
  if(rotateSecond)  et1 = new AliExternalTrackParam(pseudoX2,pseudoP2,CovPseudo2,sign);
  // negative track 
  AliAODTrack* negative = (AliAODTrack*)decay->GetDaughter(1);//Lc
  negative->GetCovarianceXYZPxPyPz(CovPseudo1);
  negative->GetXYZ(pseudoX1);
  negative->GetPxPyPz(pseudoP1);
  Short_t sign1 = negative->Charge();
  if(rotateFirst) et2 = new AliExternalTrackParam(pseudoX1,pseudoP1,CovPseudo1,sign1);
  Double_t Prot[3];

  Double_t fRot=13.;//20
  Double_t fAngle=0.0872;//0.1047
  Double_t fAngleFirst=3.14;
  if(nRot==20){
    fRot=20.;
    fAngle=0.1047;//0.0872;//5 gradi... 
  }

  // vector with rotated candidates
  TObjArray *charmArray = new TObjArray(fRot);

  Double_t d0z0[2],covd0z0[3],d0[2],d0err[2];
  Double_t d0z02[2],covd0z02[3];
  Double_t xdummy=0.,ydummy=0.,dca;
  // Double_t Angle = (TMath::Pi() - fAngle*(fRot - 1.)/2.);
  //  Double_t Angle = /*TMath::Pi()*/-6*fAngle;
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
    TObjArray ta12;
    if(rotateFirst){
      Prot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r*fAngle)*pseudoP2[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r*fAngle)*pseudoP2[1];
      Prot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r*fAngle)*pseudoP2[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r*fAngle)*pseudoP2[1];
      Prot[2] = pseudoP2[2];
      et1 = new AliExternalTrackParam(pseudoX2,Prot,CovPseudo2,sign);
    }else{
      Prot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r*fAngle)*pseudoP1[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r*fAngle)*pseudoP1[1];
      Prot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r*fAngle)*pseudoP1[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(5.*(fRot-1))/2.+ r*fAngle)*pseudoP1[1];
      Prot[2] = pseudoP1[2];
      et2 = new AliExternalTrackParam(pseudoX1,Prot,CovPseudo1,sign1);
    }
    if(nRot==20){
      if(rotateFirst){
        Prot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r*fAngle)*pseudoP2[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r*fAngle)*pseudoP2[1];
        Prot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r*fAngle)*pseudoP2[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r*fAngle)*pseudoP2[1];
        Prot[2] = pseudoP2[2];
        et1 = new AliExternalTrackParam(pseudoX2,Prot,CovPseudo2,sign);
      }else{
        Prot[0] = TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r*fAngle)*pseudoP1[0] - TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r*fAngle)*pseudoP1[1];
        Prot[1] = TMath::Sin(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r*fAngle)*pseudoP1[0] + TMath::Cos(2.*TMath::Pi()-(TMath::Pi()/180.)*(6.*(fRot-1))/2.+ r*fAngle)*pseudoP1[1];
        Prot[2] = pseudoP1[2];
        et2 = new AliExternalTrackParam(pseudoX1,Prot,CovPseudo1,sign1);
      }
    }
    ta12.Add(et1); ta12.Add(et2);
    //recalculate the secondary vertex 
    // not crucial now but needed if you rotate the momenta
    //AliAODVertex *vtxt=(AliAODVertex*)primaryVertex; //=RecalculateVertex(primaryVertex,&ta12 ,bz);
    AliAODVertex *vtxt=RecalculateVertex(primaryVertex,&ta12 ,bz);
    if(!vtxt) 
    {//cout<<"no vertex---------------------------------******************* " << endl;
      charmArray->AddAt(0,r);
      if(rotateFirst) et1->Reset();
      if(rotateSecond) et2->Reset();
      ta12.Clear();ta12=0x0;
      continue;}	
      // this are the new impact prameters
      // with relative errors 
      et1->PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
      d0[0]=d0z0[0];
      d0err[0] = TMath::Sqrt(covd0z0[0]);
      et2->PropagateToDCA(primaryVertex,bz,100.,d0z02,covd0z02);
      d0[1]=d0z02[0];
      d0err[1] = TMath::Sqrt(covd0z02[0]);

      /*   if(r==6){ 
           cout<< " propagate to DCA et1 vs et2 d00 " << "  d0z0[0]     " <<d0z0[0] << "    d0z02[0]    "<< d0z02[0]<< endl;
           cout<< " propagate to DCA et1 vs et2 d01 " << "  d0z0[1]     " <<d0z0[1] << "    d0z02[1]    "<< d0z02[1]<< endl;
           cout<< "  covd0 0 1 vs 2		 " <<covd0z0[0] << "		  "<< covd0z02[0]<<endl;
           cout<< "  covd0 1        		 " << covd0z0[1] << "  		  " << covd0z02[1]<<endl;
           cout<< "  covd0 2        		 " << covd0z0[2] << "  		  " << covd0z02[2]<<endl;
           cout<< "  covd0 3       		 " << covd0z0[3] << "  		  " << covd0z02[3]<<endl;
           cout<< "  covd0 4        		 " << covd0z0[4] << "  		  " << covd0z02[4]<<endl;
           cout<< "  covd0 5       		 " << covd0z0[5] << "  		  " << covd0z02[5]<<endl;
           cout<< "  covd0 6       		 " << covd0z0[6] << "  		  " << covd0z02[6]<<endl;
           cout<< "  covd0 7        		 " << covd0z0[7] << "  		  " << covd0z02[7]<<endl;
           cout<< "  covd0 8        		 " << covd0z0[8] << " 		  " << covd0z02[8]<<endl;
           cout<< "  covd0 9        		 " << covd0z0[9] << "  		  " << covd0z02[9]<<endl;
           cout<< "  covd0 10      		 " << covd0z0[10] << "  	  " << covd0z02[10]<<endl;
           cout<< "  covd0 11      		 " << covd0z0[11] << "  	  " << covd0z02[11]<<endl;
           cout<< "  covd0 12      		 " << covd0z0[12] << "  	  " << covd0z02[12]<<endl;
           cout<< "  covd0 13       		 " << covd0z0[13] << "  	  " << covd0z02[13]<<endl;
           cout<< "  covd0 14      		 " << covd0z0[14] << "  	  " << covd0z02[14]<<endl;
           }
           */
      //this is the new DCA      
      dca=et1->GetDCA(et2,bz,xdummy,ydummy);
      //daughters momenta
      Double_t px1[2],py1[2],pz1[2];
      if(rotateFirst){
        px1[1]= pseudoP1[0];
        py1[1]= pseudoP1[1];
        pz1[1]= pseudoP1[2];
        px1[0]= Prot[0];
        py1[0]= Prot[1];
        pz1[0]= Prot[2];
      }else{
        px1[0]= pseudoP2[0];
        py1[0]= pseudoP2[1];
        pz1[0]= pseudoP2[2];
        px1[1]= Prot[0];
        py1[1]= Prot[1];
        pz1[1]= Prot[2];
      }
      // make the rotated D0 candidate and add to the candidate vector
      // this is the new rotated candidate
      AliAODRecoDecayHF2Prong *the2Prong = new AliAODRecoDecayHF2Prong(vtxt,px1,py1,pz1,d0,d0err,dca);
      //     AliAODRecoDecayHF2Prong *the2Prong = new AliAODRecoDecayHF2Prong(primaryVertex,px1,py1,pz1,d0,d0err,dca);
      the2Prong->SetCharge(decay->Charge());
      //   the2Prong->SetPrimaryVtxRef(primaryVertex);
      the2Prong->GetSecondaryVtx()->AddDaughter(positive);
      the2Prong->GetSecondaryVtx()->AddDaughter(negative);
      the2Prong->SetProngIDs(2,id);//
      charmArray->AddAt(the2Prong,r);
      //PostData(7,fNtupleDiffD0rot);
      // delete vtxt; vtxt=NULL;
      ta12.Clear();
      ta12.Delete();
      if(rotateFirst) et1->Reset();
      if(rotateSecond) et2->Reset();
  }

  et1->Delete(); et1=0x0;
  et2->Delete(); et2=0x0;

  return charmArray;
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
  Double_t massCandLb = 0;
  UInt_t pdgLb[2]={0,0};
  //lc always at first place
  pdgLb[1] = 4122;//lambdac
  pdgLb[0] = 211;//pion
  massCandLb = part->InvMass(2,pdgLb);
  if(TMath::Abs(massCandLb - massTrueLB)>1.) return;
  Int_t iPtBinlb = -1;
  Float_t lbVar[14];
  Float_t lbVarbg[14];
  Double_t ptCandlb = part->Pt();
  // if(ptCandlb<2.) return; // dont save anything with pt<2
  if(ptCandlb>0. && ptCandlb<2.) iPtBinlb=0;
  if(ptCandlb>=2. && ptCandlb<4.) iPtBinlb=1;
  if(ptCandlb>=4. && ptCandlb<7.) iPtBinlb=2;
  if(ptCandlb>=7. && ptCandlb<10.) iPtBinlb=3;
  if(ptCandlb>=10. && ptCandlb<14.) iPtBinlb=4;
  if(ptCandlb>=14.) iPtBinlb=5;
  if(ptCandlb<2. || ptCandlb>999.) return;

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
    if(iPtBinlb==0)fInvMassLbSign0->Fill(massCandLb);
    if(iPtBinlb==1)fInvMassLbSign1->Fill(massCandLb);
    if(iPtBinlb==2)fInvMassLbSign2->Fill(massCandLb);
    if(iPtBinlb==3)fInvMassLbSign3->Fill(massCandLb);
    if(iPtBinlb==4)fInvMassLbSign4->Fill(massCandLb);
    if(iPtBinlb==5)fInvMassLbSign5->Fill(massCandLb);
    if(iPtBinlb==0)((TH1F*)fOutput->FindObject("fMassUpg_pt0lbNR"))->Fill(massCandLb);
    if(iPtBinlb==1)((TH1F*)fOutput->FindObject("fMassUpg_pt1lbNR"))->Fill(massCandLb);
    if(iPtBinlb==2)((TH1F*)fOutput->FindObject("fMassUpg_pt2lbNR"))->Fill(massCandLb);
    if(iPtBinlb==3)((TH1F*)fOutput->FindObject("fMassUpg_pt3lbNR"))->Fill(massCandLb);
    if(iPtBinlb==4)((TH1F*)fOutput->FindObject("fMassUpg_pt4lbNR"))->Fill(massCandLb);
    if(iPtBinlb==5)((TH1F*)fOutput->FindObject("fMassUpg_pt5lbNR"))->Fill(massCandLb);
      
      if(fFillNtupleSignal) {
      fNtupleLambdabUPG->Fill(lbVarC);
      PostData(2,fNtupleLambdabUPG);
    }
  }
  if(lb!=1){
    if(gen){
//      cout << " gen " << gen << endl;
      if(iPtBinlb==0)((TH1F*)fOutput->FindObject("fMassUpg_pt0lbbgNR"))->Fill(massCandLb);
      if(iPtBinlb==1)((TH1F*)fOutput->FindObject("fMassUpg_pt1lbbgNR"))->Fill(massCandLb);
      if(iPtBinlb==2)((TH1F*)fOutput->FindObject("fMassUpg_pt2lbbgNR"))->Fill(massCandLb);
      if(iPtBinlb==3)((TH1F*)fOutput->FindObject("fMassUpg_pt3lbbgNR"))->Fill(massCandLb);
      if(iPtBinlb==4)((TH1F*)fOutput->FindObject("fMassUpg_pt4lbbgNR"))->Fill(massCandLb);
      if(iPtBinlb==5)((TH1F*)fOutput->FindObject("fMassUpg_pt5lbbgNR"))->Fill(massCandLb);
        
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

  Int_t lab;
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
    Int_t lab;
    if(fApplyFixesITS3AnalysisHijing)lab=TMath::Abs(daugh->GetLabel());
    else                             lab=daugh->GetLabel();
    if(lab<0) return 0;
    if(IsTrackInjected(daugh,header,arrayMC)) return kTRUE;
  }
  return kFALSE;
}

//____________________________________________________________


