/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
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

//#####################################################
//#                                                   # 
//#          Analysis Task for Lc analysis on ESD     #
//#Authors: C. Bianchin (Utrecht University)	      #
//#         and R. Romita (Univ of Liverpool,         # 
//#         Daresbury Lab),                           #
//#         based on a class                          #
//#         by MinJung Kweon, Universitaet Heidelberg #
//#                                                   #
//#####################################################

#include "TChain.h"
#include "TTree.h"
#include "TList.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TNtuple.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include <exception>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliLog.h"
#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"

#include "AliAnalysisTaskCountLcEta.h"


ClassImp(AliAnalysisTaskCountLcEta) // adding the class to ROOT

//__________________________________________________________________
AliAnalysisTaskCountLcEta::AliAnalysisTaskCountLcEta(const char *name,const Int_t ncuts,Double_t *cuts)
: AliAnalysisTaskSE(name)
  , fESD(0)
  , fAOD(0)
  , fAnalysisType("ESD") 
  , fEvt(0) 
  , fOutList(0)
  , fEnableMCQA(kTRUE)
  , fhNevt(0)
  , fEtaAbs(0.9)
  , fEtaAbsMax(1.5)
  , fFillBkg(0)
  , fNcuts(ncuts)
  , fCuts(cuts)
  , fCutNames(0)
  , fLooserPtTrack(0)
  , fInvMassCut(0.024)
{
  // Default constructor 
	
  // Define input and output slots here
  DefineInput(0, TChain::Class()); // Input slot #0 works with a TChain
  DefineOutput(1, TList::Class()); // Output slot #0 writes into a TList container for mcQA
  const Int_t ncutsrightnow=3;
  fCutNames=new TString[ncutsrightnow];
  if(fNcuts!=ncutsrightnow) {
    Printf("ERROR!! Given %d cuts while %d are expected! Fix the class",fNcuts,ncutsrightnow);
    return;
  }
  fCutNames[0]="ptpi";
  fCutNames[1]="ptK";
  fCutNames[2]="ptp";
  Double_t pt=9999999.;
  for (Int_t i=0; i<fNcuts; i++){
    if(fCutNames[i].Contains("pt") && pt > fCuts[i]) pt=fCuts[i];
  }
  fLooserPtTrack=pt;
  Printf("INFO: Looser pt cuts = %f",fLooserPtTrack);
}


//__________________________________________________________________
AliAnalysisTaskCountLcEta::AliAnalysisTaskCountLcEta(): AliAnalysisTaskSE()
  , fESD(0)
  , fAOD(0)
  , fAnalysisType("ESD") 
  , fEvt(0) 
  , fOutList(0)
  , fEnableMCQA(kTRUE)
  , fhNevt(0)
  , fEtaAbs(0.9)
  , fEtaAbsMax(1.5)
  , fFillBkg(0)
  , fNcuts(3)
  , fCuts(0)
  , fCutNames(0)
  , fLooserPtTrack(0)
  , fInvMassCut(0.024)
{

}

//________________________________________________________________________
AliAnalysisTaskCountLcEta::~AliAnalysisTaskCountLcEta()
{
  // Destructor

  if(fESD) { delete fESD; fESD=0; }
  if(fAOD) { delete fAOD; fAOD=0; }
  if(fhNevt) { delete fhNevt; fhNevt=0; }
  if(fCuts) { delete [] fCuts; fCuts=0; }
  if(fCutNames) { delete [] fCutNames; fCutNames=0; }
  if(fOutList) { delete fOutList; fOutList=0; }


}
//__________________________________________________________________
void AliAnalysisTaskCountLcEta::UserCreateOutputObjects() {
  // Create histograms. Called once
	
  //printf("CreateOutputObjects of task %s\n", GetName());
 
  fOutList = new TList();
  fOutList->SetOwner();	
  // Open a output file
  //OpenFile(0);

  fhNevt=new TH1F("fhNevt","Number of events",1,-0.5,0.5);
  fOutList->Add(fhNevt);
  //create histograms for Lambdac

  TString hname="hLc3Prongs";
  TH1F *hLc3Prongs=new TH1F(hname,"Pt of generated Lambdac in 3 prongs;p_{T} (GeV/c)",100,0.,20.);
  hLc3Prongs->Sumw2();
  fOutList->Add(hLc3Prongs);

  hname="hPLc3Prongs";
  TH1F *hPLc3Prongs=new TH1F(hname,"P of generated Lambdac in 3 prongs;p (GeV/c)",100,0.,20.);
  hPLc3Prongs->Sumw2();
  fOutList->Add(hPLc3Prongs);

  hname="hPtpi";
  TH1F *hPtpi=new TH1F(hname,"Pt of Lc pi;p_{T} (GeV/c)",100,0.,20.);
  hPtpi->Sumw2();
  fOutList->Add(hPtpi);
  hname="hPpi";
  TH1F *hPpi=new TH1F(hname,"P of Lc pi;p (GeV/c)",100,0.,20.);
  hPpi->Sumw2();
  fOutList->Add(hPpi);

  hname="hPtp";
  TH1F *hPtp=new TH1F(hname,"Pt of Lc p;p_{T} (GeV/c)",100,0.,20.);
  hPtp->Sumw2();
  fOutList->Add(hPtp);
  hname="hPp";
  TH1F *hPp=new TH1F(hname,"P of Lc p;p (GeV/c)",100,0.,20.);
  hPp->Sumw2();
  fOutList->Add(hPp);

  hname="hPtK";
  TH1F *hPtK=new TH1F(hname,"Pt of Lc K;p_{T} (GeV/c)",100,0.,20.);
  hPtK->Sumw2();
  fOutList->Add(hPtK);
  hname="hPK";
  TH1F *hPK=new TH1F(hname,"P of Lc K;p (GeV/c)",100,0.,20.);
  hPK->Sumw2();
  fOutList->Add(hPK);

  hname="hEtaLc";
  TH1F *hEtaLc=new TH1F(hname,"Eta of Lc",400,-20.,20.);
  hEtaLc->Sumw2();
  fOutList->Add(hEtaLc);

  hname="hYLc";
  TH1F* hYLc=new TH1F(hname,"Rapidity (y) of Lc",400,-20.,20.);
  hYLc->Sumw2();
  fOutList->Add(hYLc);

  
  hname="hLc3ProngsInEta";
  TH1F *hLc3ProngsInEta=new TH1F(hname,Form("Pt of generated Lambdac in 3 prongs within %.1f ;p_{T} (GeV/c)",fEtaAbs),100,0.,20.);
  hLc3ProngsInEta->Sumw2();
  fOutList->Add(hLc3ProngsInEta);

  hname="hLc3DaughInEta";
  TH1F *hLc3DaughInEta=new TH1F(hname,Form("Pt of generated Lambdac with 3 daughters within %.1f;p_{T} (GeV/c)",fEtaAbs),100,0.,20.);
  hLc3DaughInEta->Sumw2();
  fOutList->Add(hLc3DaughInEta);

  hname="hLc2DaughInEta";
  TH1F *hLc2DaughInEta=new TH1F(hname,Form("Pt of generated Lambdac with 2 daughters within %.1f;p_{T} (GeV/c)",fEtaAbs),100,0.,20.);
  hLc2DaughInEta->Sumw2();
  fOutList->Add(hLc2DaughInEta);

  hname="hLc1DaughInEta";
  TH1F *hLc1DaughInEta=new TH1F(hname,Form("Pt of generated Lambdac with 1 daughter within %.1f;p_{T} (GeV/c)",fEtaAbs),100,0.,20.);
  hLc1DaughInEta->Sumw2();
  fOutList->Add(hLc1DaughInEta);

  TH1F* hLcpiInEta=new TH1F("hLcpiInEta", Form("Pt of generated Lambdac with #pi%.1f < |#eta| < %.1f;p_{T} (GeV/c)",fEtaAbs,fEtaAbsMax),100,0.,20.);
  hLcpiInEta->Sumw2();
  fOutList->Add(hLcpiInEta);

  TH1F* hLcpInEta=new TH1F("hLcpInEta", Form("Pt of generated Lambdac with p %.1f < |#eta| < %.1f;p_{T} (GeV/c)",fEtaAbs,fEtaAbsMax),100,0.,20.);
  hLcpInEta->Sumw2();
  fOutList->Add(hLcpInEta);

  TH1F* hLcKInEta=new TH1F("hLcKInEta", Form("Pt of generated Lambdac with K %.1f < |#eta| < %.1f;p_{T} (GeV/c)",fEtaAbs,fEtaAbsMax),100,0.,20.);
  hLcKInEta->Sumw2();
  fOutList->Add(hLcKInEta);

  TH1F* hLcpiKInEta=new TH1F("hLcpiKInEta", Form("Pt of generated Lambdac with #pi and K %.1f < |#eta| < %.1f;p_{T} (GeV/c)",fEtaAbs,fEtaAbsMax),100,0.,20.);
  hLcpiKInEta->Sumw2();
  fOutList->Add(hLcpiKInEta);

  TH1F* hLcpipInEta=new TH1F("hLcpipInEta", Form("Pt of generated Lambdac with #pi and p %.1f < |#eta| < %.1f;p_{T} (GeV/c)",fEtaAbs,fEtaAbsMax),100,0.,20.);
  hLcpipInEta->Sumw2();
  fOutList->Add(hLcpipInEta);

  TH1F* hLcKpInEta=new TH1F("hLcKpInEta", Form("Pt of generated Lambdac with K and p %.1f < |#eta| < %.1f;p_{T} (GeV/c)",fEtaAbs,fEtaAbsMax),100,0.,20.);
  hLcKpInEta->Sumw2();
  fOutList->Add(hLcKpInEta);

  TH1F* hLcpKpiInEta=new TH1F("hLcpKpiInEta", Form("Pt of generated Lambdac with #pi K and p %.1f < |#eta| < %.1f;p_{T} (GeV/c)",fEtaAbs,fEtaAbsMax),100,0.,20.);
  hLcpKpiInEta->Sumw2();
  fOutList->Add(hLcpKpiInEta);

  TH1F* hRejection=new TH1F("hRejection","Reason of track rejection",4,-0.5,3.5);
  hRejection->GetXaxis()->SetBinLabel(1,"not primary");
  hRejection->GetXaxis()->SetBinLabel(2,"out of beam pipe");
  hRejection->GetXaxis()->SetBinLabel(3,"p_{T} cut");
  hRejection->GetXaxis()->SetBinLabel(4,Form("|#eta|>%.1f",fEtaAbs));
  hRejection->GetXaxis()->SetNdivisions(1,kFALSE);
  fOutList->Add(hRejection);
  //background
  if(fFillBkg){
  	TH2F* hPtEtaBkg=new TH2F("hPtEtaBkg","P_{T} distribution of background tracks;p_{T} (GeV/c);#eta",100,0.,20.,40,-10.,10.);
  	fOutList->Add(hPtEtaBkg);
	TH2F* hPtEtaMCCandB=new TH2F("hPtEtaMCCandB","p_{T} distribution of background candidates;p_{T} (GeV/c);#eta_{acc}",100,0.,20.,60,-1.5,1.5);
	hPtEtaMCCandB->Sumw2();
	fOutList->Add(hPtEtaMCCandB);
	TH2F* hMassEtaMCCandB=new TH2F("hMassEtaMCCandB","Invariant mass distribution of background candidates;inv mass (GeV);#eta_{acc}",400,2.261,2.309,60,-1.5,1.5);
	hMassEtaMCCandB->Sumw2();
	fOutList->Add(hMassEtaMCCandB);


 }

PostData(1,fOutList);

}


//__________________________________________________________________
void AliAnalysisTaskCountLcEta::UserExec(Option_t *) {
  // Main loop. Called for each event
  fESD=dynamic_cast<AliESDEvent*> (InputEvent());
  if (!fESD) { Printf("ERROR: fESD not available"); return;}
  // MC Event
  AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!mcHandler) { printf("ERROR: Could not retrieve MC event handler"); return;}
  
  AliMCEvent* mcEvent = mcHandler->MCEvent();
  if (!mcEvent) { Printf("ERROR: Could not retrieve MC event"); return;}
  printf("MC Event %p \n",mcEvent);
  AliStack* stack = mcEvent->Stack();
  if (!stack) { printf( "Stack not available"); return;}
  const AliVVertex *vtx=mcEvent->GetPrimaryVertex();
  Double_t position[3]={0.,0.,0.};
  vtx->GetXYZ(position);

  // fill to count total number of event
  fhNevt->Fill(0);

  // run MC QA 
  if (fEnableMCQA) {

    //Int_t nPrims = stack->GetNprimary();
    Int_t nMCTracks = stack->GetNtrack();
    //histograms
    TH2F* hPtEtaBkg=(TH2F*)fOutList->FindObject("hPtEtaBkg");
    TH2F* hPtEtaCand=(TH2F*)fOutList->FindObject("hPtEtaMCCandB");
    TH2F* hMassEtaCand=(TH2F*)fOutList->FindObject("hMassEtaMCCandB");
    
    // loop over primary particles for quark and heavy hadrons
    for (Int_t igen = 0; igen < nMCTracks; igen++){
      TParticle* particle=(TParticle*)stack->Particle(igen);
      Int_t ch1=particle->GetPdgCode();
      Int_t pdgcode=TMath::Abs(ch1);
      if(ch1>0) ch1=+1;
      if(ch1<0) ch1=-1;
      
      if(pdgcode==4122) FillHistosL(particle,stack);

      if(fFillBkg && particle->IsPrimary() && (pdgcode==2212 || pdgcode==321 || pdgcode==211)) hPtEtaBkg->Fill(particle->Pt(),particle->Eta());
  
      //Printf("First loop");
      Bool_t selected=SelectTrack(particle,kTRUE);
      if(!selected) continue;
      for(Int_t j=0;fFillBkg && j<nMCTracks;j++){//second loop
	if(igen==j) continue;
	TParticle* part2=(TParticle*)stack->Particle(j);
	Int_t ch2=part2->GetPdgCode();
	pdgcode=TMath::Abs(ch2);
	if(ch2>0) ch2=+1;
	if(ch2<0) ch2=-1;

	if(!part2->IsPrimary()) continue;
	//check charge 
	if(ch1==ch2) continue;
	selected=SelectTrack(part2);
	if(!selected) continue;


	for(Int_t k=0;fFillBkg && k<nMCTracks;k++){//third loop
	  if(igen==k || j==k) continue;
	  TParticle* part3=(TParticle*)stack->Particle(k);
	  Int_t ch3=part3->GetPdgCode();
	  pdgcode=TMath::Abs(ch3);
	  if(ch3>0) ch3=+1;
	  if(ch3<0) ch3=-1;

	  if(!part3->IsPrimary()) continue;
	  //check charge 
	  if(ch1==(-1)*ch3) continue;
	  selected=SelectTrack(part3);
	  if(!selected) continue;
	  Float_t eta1=particle->Eta();
	  Float_t eta2=part2->Eta();
	  Float_t eta3=part3->Eta();
	  Float_t etamax=TMath::Abs(eta1); if(TMath::Abs(eta2)>etamax) etamax=TMath::Abs(eta2); if(TMath::Abs(eta3)>etamax) etamax=TMath::Abs(eta3);

	  //candidates
	  Bool_t hyppKpi=kTRUE;
	  hyppKpi=SelectTracksForCandidate(particle,part2,part3);
	  Bool_t hyppiKp=kTRUE;
	  hyppiKp=SelectTracksForCandidate(part3,part2,particle);
	  Double_t invmasspKpi=0,invmasspiKp=0;

	  TLorentzVector *candpKpi;
	  invmasspKpi=InvMass(particle,part2,part3,candpKpi);
	  
	  TLorentzVector *candpiKp;
	  invmasspiKp=InvMass(part3,part2,particle,candpiKp);
	  
	  Double_t mLc=TDatabasePDG::Instance()->GetParticle(4122)->Mass();
	  Double_t threesigmas[2]={mLc-fInvMassCut,mLc+fInvMassCut};
	  
	  if(invmasspKpi < threesigmas[0] || invmasspKpi > threesigmas[1] ) hyppKpi=kFALSE;
	  if(invmasspiKp < threesigmas[0] || invmasspiKp > threesigmas[1] ) hyppiKp=kFALSE;
	  if(!hyppKpi && !hyppiKp) {
	    delete candpiKp;
	    delete candpKpi;

	    continue;
	  }

	  //Double_t pcandpKpi=candpKpi->P();
	  //Double_t pcandpiKp=candpiKp->P();
	  Double_t ptcandpKpi=candpKpi->Pt();
	  Double_t ptcandpiKp=candpiKp->Pt();

	  //Fill histograms

	  Int_t binetamaxpos=hPtEtaCand->GetYaxis()->FindBin(etamax);
	  Int_t binetamaxneg=hPtEtaCand->GetYaxis()->FindBin((-1)*etamax);
	  Double_t width=hPtEtaCand->GetYaxis()->GetBinWidth(binetamaxpos);
	  //Double_t etaminhist=hPtEtaCand->GetYaxis()->GetBinLowEdge(binetamaxneg);

	  for(Int_t i=binetamaxneg;i<=binetamaxpos;i++){
	    Double_t content=hMassEtaCand->GetYaxis()->GetBinLowEdge(i)+0.5*width;

	    if(hyppKpi) {
	      hPtEtaCand->Fill(ptcandpKpi,content);
	      hMassEtaCand->Fill(candpKpi->M(),content);
	    }
	    if(hyppiKp) {
	      hPtEtaCand->Fill(ptcandpiKp,content);
	      hMassEtaCand->Fill(candpiKp->M(),content);
	      
	    }

	  }

	  delete candpiKp;
	  delete candpKpi;
	}//end third loop
      }//end second loop

    }//end loop on generated tracks

  } // end of MC QA loop


  fEvt++;		// event number

  PostData(1, fOutList);
}

Bool_t AliAnalysisTaskCountLcEta::SelectTrack(TParticle *p,Bool_t fillh){

  TH1F* hrej=(TH1F*)fOutList->FindObject("hRejection");

  if(!p->IsPrimary()) {
    if(fillh) hrej->Fill(0);
    return kFALSE;
  }
  Float_t bpipe=2.8; //cm
  Float_t zvtxdist=1.; //cm


  Double_t vx=p->Vx(),vy=p->Vy(),vz=p->Vz();
  //Printf("Production point %f, %f", vx,vy);
  if((vx*vx + vy*vy) > bpipe*bpipe || vz > zvtxdist) {
    if(fillh) hrej->Fill(1);
    return kFALSE;
  }
  Double_t pt=p->Pt();
  if(pt<fLooserPtTrack){
    if(fillh) hrej->Fill(2);
    return kFALSE;
  }
  Double_t eta=TMath::Abs(p->Eta());
  if(eta>fEtaAbs) {
    if(fillh) hrej->Fill(3);
    if(eta>1.5) return kFALSE;
  }
  return kTRUE;
}

Bool_t AliAnalysisTaskCountLcEta::SelectTracksForCandidate(TParticle* pion, TParticle* kaon, TParticle* proton){

  if(pion->Pt()<fCuts[0]) return kFALSE;
  if(kaon->Pt()<fCuts[1]) return kFALSE;
  if(proton->Pt()<fCuts[2]) return kFALSE;

  return kTRUE;
}

Double_t AliAnalysisTaskCountLcEta::InvMass(TParticle *p1p,TParticle *pn,TParticle *p2p,TLorentzVector *&candp1ppnp2p){

  Double_t pxp1=p1p->Px(),pyp1=p1p->Py(),pzp1=p1p->Pz();
  Double_t pxp2=p2p->Px(),pyp2=p2p->Py(),pzp2=p2p->Pz();
  Double_t pxpn=pn->Px(),pypn=pn->Py(),pzpn=pn->Pz();
  //Printf("pxp1 = %f, pyp1 = %f,pzp1=%f ",pxp1,pyp1,pzp1);
  Double_t p[3];
  p[0]=pxp1+pxpn+pxp2;
  p[1]=pyp1+pypn+pyp2;
  p[2]=pzp1+pzpn+pzp2;
  Double_t masspi=TDatabasePDG::Instance()->GetParticle(211)->Mass();
  Double_t massp= TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  Double_t massK= TDatabasePDG::Instance()->GetParticle(321)->Mass();
  Double_t energyP1=TMath::Sqrt(massp*massp+(p1p->P()*p1p->P()));
  Double_t energyPn=TMath::Sqrt(massK*massK+(pn->P()*pn->P()));
  Double_t energyP2=TMath::Sqrt(masspi*masspi+(p2p->P()*p2p->P()));
 

  Double_t energy=energyP1 + energyPn +energyP2;
  Double_t p2=p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
  
  Double_t invmass=TMath::Sqrt(energy*energy-p2);
  candp1ppnp2p=new TLorentzVector(p[0],p[1],p[2],energy);

  return invmass;

}
//__________________________________________________________________
void AliAnalysisTaskCountLcEta::Terminate(Option_t *) {
		
  
}

void AliAnalysisTaskCountLcEta::FillHistosL(TParticle *part, AliStack* stack){

  //Histograms
  TH1F* hPtLc3Prongs=(TH1F*)fOutList->FindObject("hLc3Prongs");
  TH1F* hPLc3Prongs=(TH1F*)fOutList->FindObject("hPLc3Prongs");
  
  TH1F* hPtpi=(TH1F*)fOutList->FindObject("hPtpi");
  TH1F* hPpi=(TH1F*)fOutList->FindObject("hPpi");
  TH1F* hPtK=(TH1F*)fOutList->FindObject("hPtK");
  TH1F* hPK=(TH1F*)fOutList->FindObject("hPK");
  TH1F* hPtp=(TH1F*)fOutList->FindObject("hPtp");
  TH1F* hPp=(TH1F*)fOutList->FindObject("hPp");
  TH1F* hLc3ProngsInEta=(TH1F*)fOutList->FindObject("hLc3ProngsInEta");
  TH1F* hLc3DaughInEta=(TH1F*)fOutList->FindObject("hLc3DaughInEta");
  TH1F* hLc2DaughInEta=(TH1F*)fOutList->FindObject("hLc2DaughInEta");
  TH1F* hLc1DaughInEta=(TH1F*)fOutList->FindObject("hLc1DaughInEta");

  TH1F* hEtaLc=(TH1F*)fOutList->FindObject("hEtaLc");
  TH1F* hYLc=(TH1F*)fOutList->FindObject("hYLc");
  
  TH1F* hLcpiInEta=(TH1F*)fOutList->FindObject("hLcpiInEta");
  TH1F* hLcpInEta=(TH1F*)fOutList->FindObject("hLcpInEta");
  TH1F* hLcKInEta=(TH1F*)fOutList->FindObject("hLcKInEta");
  TH1F* hLcpiKInEta=(TH1F*)fOutList->FindObject("hLcpiKInEta");
  TH1F* hLcpipInEta=(TH1F*)fOutList->FindObject("hLcpipInEta");
  TH1F* hLcKpInEta=(TH1F*)fOutList->FindObject("hLcKpInEta");
  TH1F* hLcpKpiInEta=(TH1F*)fOutList->FindObject("hLcpKpiInEta");
  
 
	
  Double_t pt_part=part->Pt();
  Double_t p_part=part->P();
  Double_t eta_part=part->Eta();
  Double_t y_part=part->Y();
  Int_t nDaugh = part->GetNDaughters();
  //printf("Fillhistos L\n");
  if(nDaugh<2) return;
  if(nDaugh>3) return;
  //Printf("Lc in 3 prongs");
  TParticle* pdaugh1 = stack->Particle(part->GetFirstDaughter());
  if(!pdaugh1) return;
  Int_t number1 = TMath::Abs(pdaugh1->GetPdgCode());
  TParticle* pdaugh2 = stack->Particle(part->GetLastDaughter());
  if(!pdaugh2) return;
  Int_t number2 = TMath::Abs(pdaugh2->GetPdgCode());

  Double_t eta1=0.;
  Double_t eta2=0.;
  Double_t eta3=0.;
  Bool_t okLambdac=kFALSE;
  Int_t pdgs[3]={0,0,0};

  //TMath::Abs(pdaugh2->Eta());

  //Lambda in 3 prongs
	
  Double_t mom[3];
  Double_t mom_t[3];
  if(nDaugh==3){
    //Printf("Pt part %f",pt_part);
    Int_t thirdDaugh=part->GetLastDaughter()-1;
    TParticle* pdaugh3 = stack->Particle(thirdDaugh);
    //printf("Fillhistos L 3 daugh\n");
    if(!pdaugh3) return;
    mom[0]=pdaugh1->P();
    mom[1]=pdaugh2->P();
    mom[2]=pdaugh3->P();
    mom_t[0]=pdaugh1->Pt();
    mom_t[1]=pdaugh2->Pt();
    mom_t[2]=pdaugh3->Pt();
    eta1=TMath::Abs(pdaugh1->Eta());
    eta2=TMath::Abs(pdaugh2->Eta());
    eta3=TMath::Abs(pdaugh3->Eta());
    Int_t number3 = TMath::Abs(pdaugh3->GetPdgCode());
    pdgs[0]=number1; pdgs[1]=number2; pdgs[2]=number3;
    okLambdac=kFALSE;
    if(number1==211 && number2==321 && number3==2212) okLambdac=kTRUE;
    if(number1==321 && number2==211 && number3==2212) okLambdac=kTRUE;
    if(number1==321 && number2==2212 && number3==211) okLambdac=kTRUE;
    if(number1==211 && number2==2212 && number3==321) okLambdac=kTRUE;
    if(number1==2212 && number2==211 && number3==321) okLambdac=kTRUE;
    if(number1==2212 && number2==321 && number3==211) okLambdac=kTRUE;
  }
  if(nDaugh==2){

    //Lambda resonant
    //Lambda -> p K*0
    Int_t nfiglieK=0;
    if(number1==2212 && number2==313){
      nfiglieK=pdaugh2->GetNDaughters();
      if(nfiglieK!=2) return;
      TParticle* pdaughK1 = stack->Particle(pdaugh2->GetFirstDaughter());
      TParticle* pdaughK2 = stack->Particle(pdaugh2->GetLastDaughter());
      if(!pdaughK1) return;
      if(!pdaughK2) return;
      Int_t number2K=TMath::Abs(pdaughK1->GetPdgCode());
      Int_t number3K=TMath::Abs(pdaughK2->GetPdgCode());
      mom[0]=pdaugh1->P();
      mom[1]=pdaughK1->P();
      mom[2]=pdaughK2->P();
      mom_t[0]=pdaugh1->Pt();
      mom_t[1]=pdaughK1->Pt();
      mom_t[2]=pdaughK2->Pt();
      pdgs[0]=number1; pdgs[1]=number2K; pdgs[2]=number3K;
      okLambdac=kFALSE;
      if(number1==211 && number2K==321 && number3K==2212) okLambdac=kTRUE;
      if(number1==321 && number2K==211 && number3K==2212) okLambdac=kTRUE;
      if(number1==321 && number2K==2212 && number3K==211) okLambdac=kTRUE;
      if(number1==211 && number2K==2212 && number3K==321) okLambdac=kTRUE;
      if(number1==2212 && number2K==211 && number3K==321) okLambdac=kTRUE;
      if(number1==2212 && number2K==321 && number3K==211) okLambdac=kTRUE;
      eta1=TMath::Abs(pdaugh1->Eta());
      eta2=TMath::Abs(pdaughK1->Eta());
      eta3=TMath::Abs(pdaughK2->Eta());
    }

    if(number1==313 && number2==2212){
      nfiglieK=pdaugh1->GetNDaughters();
      if(nfiglieK!=2) return;
      TParticle* pdaughK1 = stack->Particle(pdaugh1->GetFirstDaughter());
      TParticle* pdaughK2 = stack->Particle(pdaugh1->GetLastDaughter());
      if(!pdaughK1) return;
      if(!pdaughK2) return;
      Int_t number2K=TMath::Abs(pdaughK1->GetPdgCode());
      Int_t number3K=TMath::Abs(pdaughK2->GetPdgCode());
      mom[0]=pdaugh2->P();
      mom[1]=pdaughK1->P();
      mom[2]=pdaughK2->P();
      mom_t[0]=pdaugh2->Pt();
      mom_t[1]=pdaughK1->Pt();
      mom_t[2]=pdaughK2->Pt();
      pdgs[0]=number2; pdgs[1]=number2K; pdgs[2]=number3K;
      okLambdac=kFALSE;
      if(number2==211 && number2K==321 && number3K==2212) okLambdac=kTRUE;
      if(number2==321 && number2K==211 && number3K==2212) okLambdac=kTRUE;
      if(number2==321 && number2K==2212 && number3K==211) okLambdac=kTRUE;
      if(number2==211 && number2K==2212 && number3K==321) okLambdac=kTRUE;
      if(number2==2212 && number2K==211 && number3K==321) okLambdac=kTRUE;
      if(number2==2212 && number2K==321 && number3K==211) okLambdac=kTRUE;
      eta1=TMath::Abs(pdaugh2->Eta());
      eta2=TMath::Abs(pdaughK1->Eta());
      eta3=TMath::Abs(pdaughK2->Eta());
    }
    //Lambda -> Delta++ k
    Int_t nfiglieDelta=0;
    if(number1==321 && number2==2224){
      nfiglieDelta=pdaugh2->GetNDaughters();
      if(nfiglieDelta!=2) return;
      TParticle *pdaughD1=stack->Particle(pdaugh2->GetFirstDaughter());
      TParticle *pdaughD2=stack->Particle(pdaugh2->GetLastDaughter());
      if(!pdaughD1) return;
      if(!pdaughD2) return;
      Int_t number2D=TMath::Abs(pdaughD1->GetPdgCode());
      Int_t number3D=TMath::Abs(pdaughD2->GetPdgCode());
      mom[0]=pdaugh1->P();
      mom[1]=pdaughD1->P();
      mom[2]=pdaughD2->P();
      mom_t[0]=pdaugh1->Pt();
      mom_t[1]=pdaughD1->Pt();
      mom_t[2]=pdaughD2->Pt();
      okLambdac=kFALSE;
      pdgs[0]=number1; pdgs[1]=number2D; pdgs[2]=number3D;
      if(number1==211 && number2D==321 && number3D==2212) okLambdac=kTRUE;
      if(number1==321 && number2D==211 && number3D==2212) okLambdac=kTRUE;
      if(number1==321 && number2D==2212 && number3D==211) okLambdac=kTRUE;
      if(number1==211 && number2D==2212 && number3D==321) okLambdac=kTRUE;
      if(number1==2212 && number2D==211 && number3D==321) okLambdac=kTRUE;
      if(number1==2212 && number2D==321 && number3D==211) okLambdac=kTRUE;
      eta1=TMath::Abs(pdaugh1->Eta());
      eta2=TMath::Abs(pdaughD1->Eta());
      eta3=TMath::Abs(pdaughD2->Eta());
    }
    if(number1==2224 && number2==321){
      nfiglieDelta=pdaugh1->GetNDaughters();
      if(nfiglieDelta!=2) return;
      TParticle* pdaughD1 = stack->Particle(pdaugh1->GetFirstDaughter());
      TParticle* pdaughD2 = stack->Particle(pdaugh1->GetLastDaughter());
      if(!pdaughD1) return;
      if(!pdaughD2) return;
      Int_t number2D=TMath::Abs(pdaughD1->GetPdgCode());
      Int_t number3D=TMath::Abs(pdaughD2->GetPdgCode());
      mom[0]=pdaugh2->P();
      mom[1]=pdaughD1->P();
      mom[2]=pdaughD2->P();
      mom_t[0]=pdaugh2->Pt();
      mom_t[1]=pdaughD1->Pt();
      mom_t[2]=pdaughD2->Pt();
      pdgs[0]=number2; pdgs[1]=number2D; pdgs[2]=number3D;
      okLambdac=kFALSE;
      if(number2==211 && number2D==321 && number3D==2212) okLambdac=kTRUE;
      if(number2==321 && number2D==211 && number3D==2212) okLambdac=kTRUE;
      if(number2==321 && number2D==2212 && number3D==211) okLambdac=kTRUE;
      if(number2==211 && number2D==2212 && number3D==321) okLambdac=kTRUE;
      if(number2==2212 && number2D==211 && number3D==321) okLambdac=kTRUE;
      if(number2==2212 && number2D==321 && number3D==211) okLambdac=kTRUE;
      eta1=TMath::Abs(pdaugh2->Eta());
      eta2=TMath::Abs(pdaughD1->Eta());
      eta3=TMath::Abs(pdaughD2->Eta());
    }

    //Lambdac -> Lambda(1520) pi
    Int_t nfiglieLa=0;
    if(number1==3124 && number2==211){
      nfiglieLa=pdaugh1->GetNDaughters();
      if(nfiglieLa!=2) return;
      TParticle *pdaughL1=stack->Particle(pdaugh1->GetFirstDaughter());
      TParticle *pdaughL2=stack->Particle(pdaugh1->GetLastDaughter());
      if(!pdaughL1) return;
      if(!pdaughL2) return;
      Int_t number2L=TMath::Abs(pdaughL1->GetPdgCode());
      Int_t number3L=TMath::Abs(pdaughL2->GetPdgCode());
      mom[0]=pdaugh2->P();
      mom[1]=pdaughL1->P();
      mom[2]=pdaughL2->P();
      mom_t[0]=pdaugh2->Pt();
      mom_t[1]=pdaughL1->Pt();
      mom_t[2]=pdaughL2->Pt();
      pdgs[0]=number2; pdgs[1]=number2L; pdgs[2]=number3L;
      okLambdac=kFALSE;
      if(number2==211 && number2L==321 && number3L==2212) okLambdac=kTRUE;
      if(number2==321 && number2L==211 && number3L==2212) okLambdac=kTRUE;
      if(number2==321 && number2L==2212 && number3L==211) okLambdac=kTRUE;
      if(number2==211 && number2L==2212 && number3L==321) okLambdac=kTRUE;
      if(number2==2212 && number2L==211 && number3L==321) okLambdac=kTRUE;
      if(number2==2212 && number2L==321 && number3L==211) okLambdac=kTRUE;
      eta1=TMath::Abs(pdaugh2->Eta());
      eta2=TMath::Abs(pdaughL1->Eta());
      eta3=TMath::Abs(pdaughL2->Eta());
    }
    if(number1==211 && number2==3124){
      nfiglieLa=pdaugh2->GetNDaughters();
      if(nfiglieLa!=2) return;
      TParticle *pdaughL1=stack->Particle(pdaugh2->GetFirstDaughter());
      TParticle *pdaughL2=stack->Particle(pdaugh2->GetLastDaughter());
      if(!pdaughL1) return;
      if(!pdaughL2) return;
      Int_t number2L=TMath::Abs(pdaughL1->GetPdgCode());
      Int_t number3L=TMath::Abs(pdaughL2->GetPdgCode());
      mom[0]=pdaugh1->P();
      mom[1]=pdaughL1->P();
      mom[2]=pdaughL2->P();
      mom_t[0]=pdaugh1->Pt();
      mom_t[1]=pdaughL1->Pt();
      mom_t[2]=pdaughL2->Pt();
      pdgs[0]=number1; pdgs[1]=number2L; pdgs[2]=number3L;
      okLambdac=kFALSE;
      if(number1==211 && number2L==321 && number3L==2212) okLambdac=kTRUE;
      if(number1==321 && number2L==211 && number3L==2212) okLambdac=kTRUE;
      if(number1==321 && number2L==2212 && number3L==211) okLambdac=kTRUE;
      if(number1==211 && number2L==2212 && number3L==321) okLambdac=kTRUE;
      if(number1==2212 && number2L==211 && number3L==321) okLambdac=kTRUE;
      if(number1==2212 && number2L==321 && number3L==211) okLambdac=kTRUE;
      eta1=TMath::Abs(pdaugh1->Eta());
      eta2=TMath::Abs(pdaughL1->Eta());
      eta3=TMath::Abs(pdaughL2->Eta());
    }

  }
  if(!okLambdac) return;
  //Printf("Filling");
  hPtLc3Prongs->Fill(pt_part);
  hPLc3Prongs->Fill(p_part);
  hEtaLc->Fill(eta_part);
  hYLc->Fill(y_part);
  if(eta_part<fEtaAbs) hLc3ProngsInEta->Fill(pt_part);
  for(Int_t iind=0;iind<3;iind++){
    if(pdgs[iind]==211) {
      hPtpi->Fill(mom_t[iind]);
      hPpi->Fill(mom[iind]);
    }
    if(pdgs[iind]==321) {
      hPtK->Fill(mom_t[iind]);
      hPK->Fill(mom[iind]);
    }
    if(pdgs[iind]==2212) {
      hPtp->Fill(mom_t[iind]);
      hPp->Fill(mom[iind]);
    }
  }

  if(eta1<fEtaAbs && eta2<fEtaAbs && eta3<fEtaAbs){
    hLc3DaughInEta->Fill(pt_part);
  }

  if((eta1<fEtaAbs && eta2<fEtaAbs &&  eta3>fEtaAbs) || (eta1<fEtaAbs && eta3<fEtaAbs &&  eta2>fEtaAbs) || (eta3<fEtaAbs && eta2<fEtaAbs &&  eta1>fEtaAbs)){
    hLc2DaughInEta->Fill(pt_part);

    if(eta1<fEtaAbsMax && eta2<fEtaAbsMax && eta3<fEtaAbsMax){ //upper limit in eta (default 1.5)
      //PID : 1 id particle out of acceptance
      Int_t id=-1;
      if(eta1>fEtaAbs) id=0;
      if(eta2>fEtaAbs) id=1;
      if(eta3>fEtaAbs) id=2;
      if(pdgs[id]==211) hLcpiInEta->Fill(pt_part);
      if(pdgs[id]==321) hLcKInEta->Fill(pt_part);
      if(pdgs[id]==2212) hLcpInEta->Fill(pt_part);
    }
  }

  if( (eta1<fEtaAbs && eta2>fEtaAbs && eta3>fEtaAbs) || (eta2<fEtaAbs && eta1>fEtaAbs && eta3>fEtaAbs) || (eta3<fEtaAbs && eta2>fEtaAbs &&  eta1>fEtaAbs)){
    hLc1DaughInEta->Fill(pt_part);

    if(eta1<fEtaAbsMax && eta2<fEtaAbsMax && eta3<fEtaAbsMax){ //upper limit in eta (default 1.5)
      //PID : 2 id particles out of acceptance
      Int_t id1=-1,id2=-1;
      if(eta1>fEtaAbs && eta2>fEtaAbs) {
	id1=0;
	id2=1;
      }
      if(eta2>fEtaAbs && eta3>fEtaAbs) {
	id1=1;
	id2=2;
      }
      if(eta1>fEtaAbs && eta3>fEtaAbs) {
	id1=0;
	id2=2;
      }
      if((pdgs[id1]==211 && pdgs[id2]==321) || (pdgs[id1]==321 && pdgs[id2]==211))
	hLcpiKInEta->Fill(pt_part);
      if((pdgs[id1]==321 && pdgs[id2]==2212) || (pdgs[id2]==321 && pdgs[id1]==2212)) 
	hLcKpInEta->Fill(pt_part);
      if((pdgs[id1]==211 && pdgs[id2]==2212) || (pdgs[id2]==211 && pdgs[id1]==2212) ) 
	hLcpipInEta->Fill(pt_part);
    }
  }

  if(eta1<fEtaAbsMax && eta2<fEtaAbsMax && eta3<fEtaAbsMax){ //upper limit in eta (default 1.5)
    //PID : 3 id particles out of acceptance
    if (eta1>fEtaAbs && eta2>fEtaAbs && eta3>fEtaAbs) {
      hLcpKpiInEta->Fill(pt_part);
    }
  }

	 

  return;
}


