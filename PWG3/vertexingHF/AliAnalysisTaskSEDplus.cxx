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
// Author: Renu Bala, bala@to.infn.it
/////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TNtuple.h>
#include <TList.h>
#include <TString.h>
#include <TH1F.h>
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
fNtupleDplus(0),
fNtupleDplusbackg(0),
fFillNtuple(kFALSE),
fVHF(0)
{
  // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSEDplus::AliAnalysisTaskSEDplus(const char *name,Bool_t fillNtuple):
AliAnalysisTaskSE(name),
fOutput(0), 
fNtupleDplus(0),
fNtupleDplusbackg(0),
fFillNtuple(fillNtuple),
fVHF(0)
{
  // Default constructor

  // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output

  if(fFillNtuple){
    // Output slot #2 writes into a TNtuple container
    DefineOutput(2,TNtuple::Class());  //My private output
    // Output slot #3 writes into a TNtuple container
    DefineOutput(3,TNtuple::Class());  //My private output
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

  
  if(fFillNtuple){
    OpenFile(2); // 2 is the slot number of the ntuple
    fNtupleDplus = new TNtuple("fNtupleDplus","D +","pdg:Px:Py:Pz:Ptpi:PtK:Ptpi2:PtRec:PtTrue:PointingAngle:DecLeng:VxTrue:VxRec:VyRec:VzRec:InvMass:sigvert:d0Pi:d0K:d0Pi2");
    OpenFile(3); // 3 is the slot number of the ntuple
    fNtupleDplusbackg = new TNtuple("fNtupleDplusbackg","D + backg","Ptpibkg:Ptkbkg:Ptpi2bkg:PtRecbkg:PointingAnglebkg:DLbkg:VxRecbkg:VyRecbkg:VzRecbkg:InvMassbkg:sigvertbkg:d0Pibkg:d0Kbkg:d0Pi2bkg");
  }

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEDplus::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth

  
  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());

  // In case there is an AOD handler writing a standard AOD, use the AOD 
  // event in memory rather than the input (ESD) event.
  if (!aod && AODEvent() && IsStandardAOD()) aod = dynamic_cast<AliAODEvent*> (AODEvent());

  // load Dplus->Kpipi candidates                                                   
  TClonesArray *array3Prong =
    (TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
  if(!array3Prong) {
    printf("AliAnalysisTaskSEDplus::UserExec: Charm3Prong branch not found!\n");
    return;
  }

  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  //    vtx1->Print();
  
  // load MC particles
  TClonesArray *arrayMC = 
    (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if(!arrayMC) {
    printf("AliAnalysisTaskSEDplus::UserExec: MC particles branch not found!\n");
    return;
  }
    
  // load MC header
  AliAODMCHeader *mcHeader = 
    (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  if(!mcHeader) {
    printf("AliAnalysisTaskSEDplus::UserExec: MC header branch not found!\n");
    return;
  }
  Int_t n3Prong = array3Prong->GetEntriesFast();
  if(fDebug>1) printf("Number of D+->Kpipi: %d\n",n3Prong);


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
 	cutsDplus[10]=0.979;
      }
      else if(ptCand>2. && ptCand<3){ 
	iPtBin=1;
	cutsDplus[7]=0.08;
 	cutsDplus[8]=0.5;
 	cutsDplus[9]=0.991;
      }else if(ptCand>3. && ptCand<5){ 
	iPtBin=2;
	cutsDplus[7]=0.1;
 	cutsDplus[8]=0.5;
 	cutsDplus[9]=0.9955;
      }else{
	iPtBin=3;
	cutsDplus[7]=0.1;
 	cutsDplus[8]=0.5;
  	cutsDplus[9]=0.997;
      }
      Bool_t passTightCuts=d->SelectDplus(cutsDplus);
      Int_t labDp = d->MatchToMC(411,arrayMC,3,pdgDgDplustoKpipi);
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


      if(labDp>=0) {
	((TH1F*)(fOutput->FindObject(hisNameS)))->Fill(invMass);
	if(passTightCuts){
	  ((TH1F*)(fOutput->FindObject(hisNameSTC)))->Fill(invMass);
	}
	if(fFillNtuple){
	  AliAODMCParticle *partDp = (AliAODMCParticle*)arrayMC->At(labDp);
	  AliAODMCParticle *dg0 = (AliAODMCParticle*)arrayMC->At(partDp->GetDaughter(0));	  
	  Float_t tmp[20];
	  tmp[0]=TMath::Abs(partDp->GetPdgCode());
	  tmp[1]=partDp->Px()-d->Px();
	  tmp[2]=partDp->Py()-d->Py();
	  tmp[3]=partDp->Pz()-d->Pz();
	  tmp[4]=d->PtProng(0);
	  tmp[5]=d->PtProng(1);
	  tmp[6]=d->PtProng(2);
	  tmp[7]=d->Pt();
	  tmp[8]=partDp->Pt();
	  tmp[9]=d->CosPointingAngle();
	  tmp[10]=d->DecayLength();
	  tmp[11]=dg0->Xv();
	  tmp[12]=d->Xv();
	  tmp[13]=d->Yv();
	  tmp[14]=d->Zv();
	  tmp[15]=d->InvMassDplus();
	  tmp[16]=d->GetSigmaVert();
	  tmp[17]=d->Getd0Prong(0);
	  tmp[18]=d->Getd0Prong(1);
	  tmp[19]=d->Getd0Prong(2);	  
	  fNtupleDplus->Fill(tmp);
	  PostData(2,fNtupleDplus);
	}
      }else{
	((TH1F*)(fOutput->FindObject(hisNameB)))->Fill(invMass);
	if(passTightCuts){
	  ((TH1F*)(fOutput->FindObject(hisNameBTC)))->Fill(invMass);
	}
	if(fFillNtuple){
	  Float_t tmpbkg[14];
	  tmpbkg[0]=d->PtProng(0);
	  tmpbkg[1]=d->PtProng(1);
	  tmpbkg[2]=d->PtProng(2);
	  tmpbkg[3]=d->Pt();
	  tmpbkg[4]=d->CosPointingAngle();
	  tmpbkg[5]=d->DecayLength();
	  tmpbkg[6]=d->Xv();
	  tmpbkg[7]=d->Yv();
	  tmpbkg[8]=d->Zv();
	  tmpbkg[9]=d->InvMassDplus();
	  tmpbkg[10]=d->GetSigmaVert();
	  tmpbkg[11]=d->Getd0Prong(0);
	  tmpbkg[12]=d->Getd0Prong(1);
	  tmpbkg[13]=d->Getd0Prong(2);
	  fNtupleDplusbackg->Fill(tmpbkg);	    
	  PostData(3,fNtupleDplusbackg);
	}
      }
      PostData(1,fOutput);    
    }
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  }
  
  return;
}

//________________________________________________________________________
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

  if(fFillNtuple){
    fNtupleDplus = dynamic_cast<TNtuple*>(GetOutputData(2));
    fNtupleDplusbackg = dynamic_cast<TNtuple*>(GetOutputData(3));
  }

  return;
}

