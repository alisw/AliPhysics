/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSE for the extraction of signal(e.g D+) of heavy flavor
// decay candidates with the MC truth.
// Authors: Renu Bala, bala@to.infn.it
// F. Prino, prino@to.infn.it
// G. Ortona, ortona@to.infn.it
/////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TNtuple.h>
#include <TList.h>
#include <TString.h>
#include <TH1F.h>
#include <TDatabasePDG.h>

#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEDplus.h"

ClassImp(AliAnalysisTaskSEDplus)


//________________________________________________________________________
AliAnalysisTaskSEDplus::AliAnalysisTaskSEDplus():
AliAnalysisTaskSE(),
fOutput(0), 
fHistNEvents(0),
fNtupleDplus(0),
fHistLS(0),
fHistLStrip(0),
fHistLStripcut(0),
fHistOS(0),
fHistOSbkg(0),
fHistLSnoweight(0),
fHistLSsinglecut(0),
fFillNtuple(kFALSE),
fReadMC(kFALSE),
fDoLS(kFALSE),
fVHF(0)
{
  // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSEDplus::AliAnalysisTaskSEDplus(const char *name,Bool_t fillNtuple):
AliAnalysisTaskSE(name),
fOutput(0),
fHistNEvents(0),
fNtupleDplus(0),
fHistLS(0),
fHistLStrip(0),
fHistLStripcut(0),
fHistOS(0),
fHistOSbkg(0),
fHistLSnoweight(0),
fHistLSsinglecut(0),
fFillNtuple(fillNtuple),
fReadMC(kFALSE),
fDoLS(kFALSE),
fVHF(0)
{
  // Default constructor

  // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output

  if(fFillNtuple){
    // Output slot #2 writes into a TNtuple container
    DefineOutput(2,TNtuple::Class());  //My private output
  }
}

//________________________________________________________________________
AliAnalysisTaskSEDplus::~AliAnalysisTaskSEDplus()
{
  // Destructor
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
  if(fNtupleDplus){
    delete fNtupleDplus;
    fNtupleDplus=0;
  }
  if (fVHF) {
    delete fVHF;
    fVHF = 0;
  }

}  

//________________________________________________________________________
void AliAnalysisTaskSEDplus::Init()
{
  // Initialization

  if(fDebug > 1) printf("AnalysisTaskSEDplus::Init() \n");

  gROOT->LoadMacro("ConfigVertexingHF.C");

  fVHF = (AliAnalysisVertexingHF*)gROOT->ProcessLine("ConfigVertexingHF()");  
  fVHF->PrintStatus();

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDplus::UserCreateOutputObjects()
{
  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSEDplus::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  Int_t nPtBins=4;
  TString hisname;
  for(Int_t i=0;i<nPtBins;i++){
    hisname.Form("hMassPt%d",i);
    TH1F* hm=new TH1F(hisname.Data(),hisname.Data(),100,1.765,1.965);
    hm->Sumw2();
    hisname.Form("hSigPt%d",i);
    TH1F* hs=new TH1F(hisname.Data(),hisname.Data(),100,1.765,1.965);
    hs->Sumw2();
    hisname.Form("hBkgPt%d",i);
    TH1F* hb=new TH1F(hisname.Data(),hisname.Data(),100,1.765,1.965);
    hb->Sumw2();

    hisname.Form("hMassPt%dTC",i);
    TH1F* hmtc=new TH1F(hisname.Data(),hisname.Data(),100,1.765,1.965);
    hmtc->Sumw2();
    hisname.Form("hSigPt%dTC",i);
    TH1F* hstc=new TH1F(hisname.Data(),hisname.Data(),100,1.765,1.965);
    hstc->Sumw2();
    hisname.Form("hBkgPt%dTC",i);
    TH1F* hbtc=new TH1F(hisname.Data(),hisname.Data(),100,1.765,1.965);
    hbtc->Sumw2();


    fOutput->Add(hm);
    fOutput->Add(hs);
    fOutput->Add(hb);
    fOutput->Add(hmtc);
    fOutput->Add(hstc);
    fOutput->Add(hbtc);
  }
  fHistNEvents = new TH1F("fHistNEvents", "Number of processed events; ; Events",3,-1.5,1.5);
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);
  if(fDoLS){
    Double_t massDplus=TDatabasePDG::Instance()->GetParticle(411)->Mass();
  
    fHistLS = new TH1F("fHistLS","fHistLS",100,massDplus-0.2,massDplus+0.2);
    fHistLS->Sumw2();
    fOutput->Add(fHistLS);

    fHistLStrip = new TH1F("fHistLStrip","fHistLStrip",100,massDplus-0.2,massDplus+0.2);
    fHistLStrip->Sumw2();
    fOutput->Add(fHistLStrip);

    fHistLStripcut = new TH1F("fHistLStripcut","fHistLStripcut",100,massDplus-0.2,massDplus+0.2);
    fHistLStripcut->Sumw2();
    fOutput->Add(fHistLStripcut);

    fHistOS = new TH1F("fHistOS","fHistOS",100,massDplus-0.2,massDplus+0.2);
    fHistOS->Sumw2();
    fOutput->Add(fHistOS);

    fHistOSbkg = new TH1F("fHistOSbkg","fHistOSbkg",100,massDplus-0.2,massDplus+0.2);
    fHistOSbkg->Sumw2();
    fOutput->Add(fHistOSbkg);

    fHistLSnoweight = new TH1F("fHistLSnoweight","fHistLSnoweight",100,massDplus-0.2,massDplus+0.2);
    fHistLSnoweight->Sumw2();
    fOutput->Add(fHistLSnoweight);

    fHistLSsinglecut = new TH1F("fHistLSsinglecut","fHistLSsinglecut",100,massDplus-0.2,massDplus+0.2);
    fHistLSsinglecut->Sumw2();
    fOutput->Add(fHistLSsinglecut);
  }

  if(fFillNtuple){
    OpenFile(2); // 2 is the slot number of the ntuple
    fNtupleDplus = new TNtuple("fNtupleDplus","D +","pdg:Px:Py:Pz:PtTrue:VxTrue:VyTrue:VzTrue:Ptpi:PtK:Ptpi2:PtRec:PointingAngle:DecLeng:VxRec:VyRec:VzRec:InvMass:sigvert:d0Pi:d0K:d0Pi2");
  }

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDplus::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth

  
  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  fHistNEvents->Fill(0); // count event
  // Post the data already here
  PostData(1,fOutput);
  

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
  } else {
    array3Prong=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
    arrayLikeSign=(TClonesArray*)aod->GetList()->FindObject("LikeSign3Prong");
  }

  if(!array3Prong) {
    printf("AliAnalysisTaskSEDplus::UserExec: Charm3Prong branch not found!\n");
    return;
  }
  if(!arrayLikeSign) {
    printf("AliAnalysisTaskSEDplus::UserExec: LikeSign3Prong branch not found!\n");
    // do in any case the OS analysis.
  }

 
  TClonesArray *arrayMC=0;
  AliAODMCHeader *mcHeader=0;

  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  //    vtx1->Print();
  
  // load MC particles
  if(fReadMC){
    
    arrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC) {
      printf("AliAnalysisTaskSEDplus::UserExec: MC particles branch not found!\n");
      return;
    }
    
  // load MC header
    mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSEDplus::UserExec: MC header branch not found!\n");
      return;
    }
  }
  
  Int_t n3Prong = array3Prong->GetEntriesFast();
  if(fDebug>1) printf("Number of D+->Kpipi: %d\n",n3Prong);
  
  
  Int_t nOS=0;
  Int_t pdgDgDplustoKpipi[3]={321,211,211};
  Double_t cutsDplus[12]={0.2,0.4,0.4,0.,0.,0.01,0.06,0.02,0.,0.85,0.,10000000000.};
  for (Int_t i3Prong = 0; i3Prong < n3Prong; i3Prong++) {
    AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(i3Prong);
    
    
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()){
      d->SetOwnPrimaryVtx(vtx1);
      unsetvtx=kTRUE;
    }

    if(d->SelectDplus(fVHF->GetDplusCuts())) {
      Int_t iPtBin=0;
      Double_t ptCand = d->Pt();
      if(ptCand<2.){
	iPtBin=0;
	cutsDplus[7]=0.08;
 	cutsDplus[8]=0.5;
 	cutsDplus[9]=0.979;
	cutsDplus[10]=0.0000055;
      }
      else if(ptCand>2. && ptCand<3){ 
	iPtBin=1;
	cutsDplus[7]=0.08;
 	cutsDplus[8]=0.5;
 	cutsDplus[9]=0.991;
	cutsDplus[10]=0.000005;
      }else if(ptCand>3. && ptCand<5){ 
	iPtBin=2;
	cutsDplus[7]=0.1;
 	cutsDplus[8]=0.5;
 	cutsDplus[9]=0.9955;
	cutsDplus[10]=0.0000035;
      }else{
	iPtBin=3;
	cutsDplus[7]=0.1;
 	cutsDplus[8]=0.5;
  	cutsDplus[9]=0.997;
	cutsDplus[10]=0.000001;
      }
      Bool_t passTightCuts=d->SelectDplus(cutsDplus);
      Int_t labDp=-1;
      Float_t deltaPx=0.;
      Float_t deltaPy=0.;
      Float_t deltaPz=0.;
      Float_t truePt=0.;
      Float_t xDecay=0.;
      Float_t yDecay=0.;
      Float_t zDecay=0.;
      Float_t pdgCode=-2;
      if(fReadMC){
	labDp = d->MatchToMC(411,arrayMC,3,pdgDgDplustoKpipi);
	if(labDp>=0){
	  AliAODMCParticle *partDp = (AliAODMCParticle*)arrayMC->At(labDp);
	  AliAODMCParticle *dg0 = (AliAODMCParticle*)arrayMC->At(partDp->GetDaughter(0));
	  deltaPx=partDp->Px()-d->Px();
	  deltaPy=partDp->Py()-d->Py();
	  deltaPz=partDp->Pz()-d->Pz();
	  truePt=partDp->Pt();
	  xDecay=dg0->Xv();	  
	  yDecay=dg0->Yv();	  
	  zDecay=dg0->Zv();
	  pdgCode=TMath::Abs(partDp->GetPdgCode());
	}else{
	  pdgCode=-1;
	}
      }
      Double_t invMass=d->InvMassDplus();

      TString hisNameA(Form("hMassPt%d",iPtBin));
      TString hisNameS(Form("hSigPt%d",iPtBin));
      TString hisNameB(Form("hBkgPt%d",iPtBin));
      TString hisNameATC(Form("hMassPt%dTC",iPtBin));
      TString hisNameSTC(Form("hSigPt%dTC",iPtBin));
      TString hisNameBTC(Form("hBkgPt%dTC",iPtBin));
      
      ((TH1F*)(fOutput->FindObject(hisNameA)))->Fill(invMass);
      if(passTightCuts){
	((TH1F*)(fOutput->FindObject(hisNameATC)))->Fill(invMass);
      }

      Float_t tmp[22];
      if(fFillNtuple){  	  
	tmp[0]=pdgCode;
	tmp[1]=deltaPx;
	tmp[2]=deltaPy;
	tmp[3]=deltaPz;
	tmp[4]=truePt;
	tmp[5]=xDecay;	  
	tmp[6]=yDecay;	  
	tmp[7]=zDecay;	  
	tmp[8]=d->PtProng(0);
	tmp[9]=d->PtProng(1);
	tmp[10]=d->PtProng(2);
	tmp[11]=d->Pt();
	tmp[12]=d->CosPointingAngle();
	tmp[13]=d->DecayLength();
	tmp[14]=d->Xv();
	tmp[15]=d->Yv();
	tmp[16]=d->Zv();
	tmp[17]=d->InvMassDplus();
	tmp[18]=d->GetSigmaVert();
	tmp[19]=d->Getd0Prong(0);
	tmp[20]=d->Getd0Prong(1);
	tmp[21]=d->Getd0Prong(2);	  
	fNtupleDplus->Fill(tmp);
	PostData(2,fNtupleDplus);
      }

      if(fReadMC){
	if(labDp>=0) {
	  ((TH1F*)(fOutput->FindObject(hisNameS)))->Fill(invMass);
	  if(passTightCuts){
	    ((TH1F*)(fOutput->FindObject(hisNameSTC)))->Fill(invMass);
	  }
	  
	}else{
	  ((TH1F*)(fOutput->FindObject(hisNameB)))->Fill(invMass);
	  if(passTightCuts){
	    ((TH1F*)(fOutput->FindObject(hisNameBTC)))->Fill(invMass);
	  }
 	  fHistOSbkg->Fill(invMass);
	}
      }

      fHistOS->Fill(invMass);
      nOS++;
    }
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  }
 
  //start LS analysis
  if(fDoLS && arrayLikeSign) LSAnalysis(array3Prong,arrayLikeSign,aod,vtx1,nOS);
   
  PostData(1,fOutput);    
  return;
}



//_________________________________________________________________
void AliAnalysisTaskSEDplus::LSAnalysis(TClonesArray *arrayOppositeSign,TClonesArray *arrayLikeSign,AliAODEvent *aod,AliAODVertex *vtx1, Int_t nDplusOS){

/*
 * Fill the Like Sign histograms
 */




  //count pos/neg tracks
  Int_t nPosTrks=0,nNegTrks=0;
  //counter for particles passing single particle cuts
  Int_t nspcplus=0;
  Int_t nspcminus=0;

  for(Int_t it=0;it<aod->GetNumberOfTracks();it++) {
    AliAODTrack *track = aod->GetTrack(it);
    if(track->Charge()>0){
      nPosTrks++;
      if(track->Pt()>=0.4){
	nspcplus++;
      }
    }else{      
      nNegTrks++;
      if(track->Pt()>=0.4){
	nspcminus++;
      }
    }
  }

  Int_t nOStriplets = arrayOppositeSign->GetEntriesFast();

  Int_t nDplusLS=0;
  Int_t nLikeSign = arrayLikeSign->GetEntriesFast();
 
  for(Int_t iLikeSign = 0; iLikeSign < nLikeSign; iLikeSign++) {
    AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)arrayLikeSign->UncheckedAt(iLikeSign);
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()) {
      d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
      unsetvtx=kTRUE;
    }
    if(d->SelectDplus(fVHF->GetDplusCuts()))nDplusLS++;
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  }

  Float_t wei2=0;
  if(nLikeSign!=0)wei2 = (Float_t)nOStriplets/(Float_t)nLikeSign;
  Float_t wei3=0;
  if(nDplusLS!=0)wei3 = (Float_t)nDplusOS/(Float_t)nDplusLS;

 // loop over  sign candidates
  for(Int_t iLikeSign = 0; iLikeSign < nLikeSign; iLikeSign++) {
    AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)arrayLikeSign->UncheckedAt(iLikeSign);
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()) {
      d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
      unsetvtx=kTRUE;
    }
 
    if(d->SelectDplus(fVHF->GetDplusCuts())){
      Int_t sign= d->GetCharge();
      Float_t wei1=1;
      Float_t wei4=1;
      if(sign>0&&nPosTrks>2&&nspcplus>2){//wei* should be automatically protected, since to get a triplet there must be at least 3 good tracks in the event
	wei1=3.*(Float_t)nNegTrks/((Float_t)nPosTrks-2.);
	wei4=3.*(Float_t)nspcminus/((Float_t)nspcplus-2.);
      }
      if(sign<0&&nNegTrks>2&&nspcminus>2){
	wei1=3.*(Float_t)nPosTrks/((Float_t)nNegTrks-2.);
	wei4=3.*(Float_t)nspcplus/((Float_t)nspcminus-2.);
      }
      Double_t invMass=d->InvMassDplus();
      fHistLS->Fill(invMass,wei1);
      fHistLSnoweight->Fill(invMass);
      fHistLStrip->Fill(invMass,wei2);
      fHistLStripcut->Fill(invMass,wei3);
      fHistLSsinglecut->Fill(invMass,wei4);
     }
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  }
  
  //printf("------------ N. of positive tracks in Event ----- %d \n", nPosTrks);
  //printf("------------ N. of negative tracks in Event ----- %d \n", nNegTrks);

  //  printf("LS analysis...done\n");

}

//_________________________________________________________________


void AliAnalysisTaskSEDplus::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSEDplus: Terminate() \n");
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
  fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEvents"));
  if(fFillNtuple){
    fNtupleDplus = dynamic_cast<TNtuple*>(GetOutputData(2));
  }
  fHistOS =  dynamic_cast<TH1F*>(fOutput->FindObject("fHistOS"));
  fHistOSbkg =  dynamic_cast<TH1F*>(fOutput->FindObject("fHistOSbkg"));
  if(fDoLS){
    fHistLS =  dynamic_cast<TH1F*>(fOutput->FindObject("fHistLS"));
    fHistLStrip =  dynamic_cast<TH1F*>(fOutput->FindObject("fHistLStrip"));
    fHistLStripcut =  dynamic_cast<TH1F*>(fOutput->FindObject("fHistLStripcut"));
    fHistLSnoweight =  dynamic_cast<TH1F*>(fOutput->FindObject("fHistLSnoweight"));
    fHistLSsinglecut =  dynamic_cast<TH1F*>(fOutput->FindObject("fHistLSsinglecut"));
  }


  return;
}
