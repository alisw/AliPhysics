/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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
// Class AliAnalysisTaskCharmFraction
// AliAnalysisTask for the extraction of the fraction of prompt charm
// using the charm hadron impact parameter to the primary vertex
//
// Author: Andrea Rossi, andrea.rossi@ts.infn.it
/////////////////////////////////////////////////////////////

#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TMath.h>

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecay.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisTaskCharmFraction.h"

ClassImp(AliAnalysisTaskCharmFraction)
 
  //________________________________________________________________________
  AliAnalysisTaskCharmFraction::AliAnalysisTaskCharmFraction(const char *name) 
    : AliAnalysisTask(name, ""), fAOD(0), fVHF(0),
      fArrayD0toKpi(0),
      fArrayMC(0),
      fAODmcHeader(0),
      fhMass(0),
      fhMassTrue(0),
      fhCPtaVSd0d0(0),
      fhd0xd0(0),
      fhCPta(0),
      fhSecVtxZ(0),
      fhSecVtxX(0),
      fhSecVtxY(0),
      fhSecVtxXY(0),
      fhSecVtxPhi(0),
      fhd0D0(0),
      fhd0D0pt(0x0),
      fhd0D0VtxTrue(0),
      fhd0D0VtxTruept(0x0),
      fhMCd0D0(0),
      fhMCd0D0pt(0x0),
      fnbins(),
      fD0usecuts(kTRUE),
      fcheckMC(kFALSE),
      fcheckMCD0(kFALSE),
      fcheckMC2prongs(kFALSE),
      fcheckMCprompt(kFALSE),
      fcheckMCfromB(kFALSE),
      fSkipD0star(kFALSE),
      fStudyd0fromBTrue(kTRUE),
      fStudyPureBackground(kFALSE),
      fSideBands(0)
{
  // Constructor
  //Side bands selection (see cxx): 999 -> no inv mass selection at all,
  //                                >0 ->  both D0 and D0bar hypothesis are required to be kFALSE, 
  //                                <0 -> only 1 inv mass hypothesis is required to be KFALSE
 
  SetNPtBins();
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  //  DefineInput(1, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0, TH2F::Class());// hCPtaVsd0d0
  DefineOutput(1, TH2F::Class());// hVtxXY
  for(Int_t j=2;j<3*fnbins+11;j++){ // hCPta, hd0xd0, hVtZ, hVtX hVtY hVtPhi + nptbin hd0D0 +hd0D0 integreated +nptbin hd0VtxTrueD0 +hd0VtxTrueD0 integrated +nptbin hMCd0D0 +hMCd0D0 integrated
    DefineOutput(j, TH1F::Class());
  }
  //  DefineOutput(fnbins+8, TH1F::Class());// hd0D0 pt integrated
}

//________________________________________________________________________
AliAnalysisTaskCharmFraction::AliAnalysisTaskCharmFraction(const char *name,Int_t nptbins) 
  : AliAnalysisTask(name, ""), fAOD(0), fVHF(0),
    fArrayD0toKpi(0),
    fArrayMC(0),
    fAODmcHeader(0),
    fhMass(0),
    fhMassTrue(0),
    fhCPtaVSd0d0(0),
    fhd0xd0(0),
    fhCPta(0),
    fhSecVtxZ(0),
    fhSecVtxX(0),
    fhSecVtxY(0),
    fhSecVtxXY(0),
    fhSecVtxPhi(0),
    fhd0D0(0),
    fhd0D0pt(0x0),
    fhd0D0VtxTrue(0),
    fhd0D0VtxTruept(0x0),
    fhMCd0D0(0),
    fhMCd0D0pt(0x0),
    fnbins(),
    fD0usecuts(kTRUE),
    fcheckMC(kFALSE),
    fcheckMCD0(kFALSE),
    fcheckMC2prongs(kFALSE),
    fcheckMCprompt(kFALSE),
    fcheckMCfromB(kFALSE),
    fSkipD0star(kFALSE),
    fStudyd0fromBTrue(kTRUE),
    fStudyPureBackground(kFALSE),
    fSideBands(0)        
{
  // Constructor
  //Side bands selection (see cxx): 999 -> no inv mass selection at all,
  //                                >0 ->  both D0 and D0bar hypothesis are required to be kFALSE, 
  //                                <0 -> only 1 inv mass hypothesis is required to be KFALSE
  SetNPtBins(nptbins);
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  //  DefineInput(1, TChain::Class());
  // Output slot #0 writes into a TH1 container
  
  DefineOutput(0, TH2F::Class());// hCPtaVsd0d0
  DefineOutput(1, TH2F::Class());// hVtxXY
  for(Int_t j=2;j<3*fnbins+13;j++){ // hCPta, hd0xd0, hVtZ, hVtX hVtY hVtPhi + nptbin hd0D0 +hd0D0 integreated +nptbin hd0VtxTrueD0 +hd0VtxTrueD0 integrated +nptbin hMCd0D0 +hMCd0D0 integrated
    DefineOutput(j, TH1F::Class());
  }
  //  DefineOutput(fnbins+8, TH1F::Class());// hd0D0 pt integrated 
}

//________________________________________________________________________
void AliAnalysisTaskCharmFraction::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  //  TTree* treeFriend = dynamic_cast<TTree*> (GetInputData(1));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    //  tree->AddFriend(treeFriend);
    
    // Disable all branches and enable only the needed ones
    // The next two lines are different when data produced as AliESDEvent is read
    
    AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    fAOD=new AliAODEvent();
    fAOD->ReadFromTree(tree);
    
    fArrayD0toKpi = (TClonesArray*)fAOD->GetList()->FindObject("D0toKpi"); 
    fArrayMC =      (TClonesArray*)fAOD->GetList()->FindObject("mcparticles"); 
    
    fAODmcHeader=(AliAODMCHeader*)fAOD->GetList()->FindObject("mcHeader");
    
    if (!aodH) {
      Printf("ERROR: Could not get AODInputHandler");
    } else
      fAOD = aodH->GetEvent();
    
  }
}

//________________________________________________________________________
void AliAnalysisTaskCharmFraction::CreateOutputObjects()
{
  // Create histograms
  // Called once
  fhd0D0 = new TH1F("hd0D0","D^{0} impact par. plot  (All momenta)",1000,-1000.,1000.);
  fhd0D0->SetXTitle("Impact parameter [#mum]");
  fhd0D0->SetYTitle("Entries");

  fhd0D0VtxTrue = new TH1F("hd0D0VtxTrue","D^{0} impact par. w.r.t. True Vtx (All momenta)",1000,-1000.,1000.);
  fhd0D0VtxTrue->SetXTitle("Impact parameter [#mum]");
  fhd0D0VtxTrue->SetYTitle("Entries");

  fhMCd0D0 = new TH1F("hMCd0D0","D^{0} impact par. plot  (All momenta)",1000,-1000.,1000.);
  fhMCd0D0->SetXTitle("MC Impact parameter [#mum]");
  fhMCd0D0->SetYTitle("Entries");

  fhCPtaVSd0d0=new TH2F("hCPtaVSd0d0","hCPtaVSd0d0",1000,-100000.,100000.,100,0.,1.);
  fhSecVtxZ=new TH1F("hSecVtxZ","hSecVtxZ",1000,-8.,8.);
  fhSecVtxX=new TH1F("hSecVtxX","hSecVtxX",1000,-3000.,3000.);
  fhSecVtxY=new TH1F("hSecVtxY","hSecVtxY",1000,-3000.,3000.);
  fhSecVtxXY=new TH2F("hSecVtxXY","hSecVtxXY",1000,-3000.,3000.,1000,-3000.,3000.);
  fhSecVtxPhi=new TH1F("hSecVtxPhi","hSecVtxPhi",180,-180.1,180.1);
  fhCPta=new TH1F("hCPta","hCPta",100,0.,1.);
  fhd0xd0=new TH1F("hd0xd0","hd0xd0",1000,-100000.,100000.);
  fhMassTrue=new TH1F("hMassTrue","D^{0} MC inv. Mass (All momenta)",600,1.600,2.200);
  fhMass=new TH1F("hMass","D^{0} inv. Mass (All momenta)",600,1.600,2.200);
  Double_t standardptbins[14]={0.,0.5,1.,2.,3.,5.,7.,10.,15.,20.,30.,40.,50.,100.};
  
  fhd0D0pt=new TH1F*[fnbins];
  fhMCd0D0pt=new TH1F*[fnbins];
  fhd0D0VtxTruept=new TH1F*[fnbins];
  TString namehist="hd0D0pt";
  TString titlehist="D^{0} impact par. plot ";
  TString strnamept,strtitlept;
  for(Int_t i=0;i<fnbins;i++){
    strnamept=namehist;
    strnamept+=i;// strnamept+=standardptbins[i+1];
    //    strnamept.Append("_");
    //strnamept+=standardptbins[i+1];

    strtitlept=namehist;
    strtitlept+=standardptbins[i];
    strtitlept.Append("<= pt <");
    strtitlept+=standardptbins[i+1];
    strtitlept.Append(" [GeV/c]");
    
    fhd0D0pt[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    fhd0D0pt[i]->SetXTitle("Impact parameter [#mum] ");
    fhd0D0pt[i]->SetYTitle("Entries");

    strnamept.ReplaceAll("hd0D0","hMCd0D0");
    fhMCd0D0pt[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    fhMCd0D0pt[i]->SetXTitle("MC Impact parameter [#mum] ");
    fhMCd0D0pt[i]->SetYTitle("Entries");

    strnamept.ReplaceAll("hMCd0D0","hd0D0VtxTrue");
    fhd0D0VtxTruept[i] = new TH1F(strnamept.Data(),strtitlept.Data(),1000,-1000.,1000.);
    fhd0D0VtxTruept[i]->SetXTitle("Impact parameter w.r.t. True Vtx [#mum] ");
    fhd0D0VtxTruept[i]->SetYTitle("Entries");
    
  }


  /*  fHistPt = new TH1F("fHistPt", "P_{T} distribution", 15, 0.1, 3.1);
  fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPt->SetMarkerStyle(kFullCircle);*/
}

//________________________________________________________________________
void AliAnalysisTaskCharmFraction::Exec(Option_t *) 
{
  // Main loop
  // Called for each event

  if (!fAOD) {
    Printf("ERROR: fAOD not available");
    return;
  }

  //Printf("There are %d tracks in this event", fAOD->GetNumberOfTracks());
  //  Int_t nTotHF=0,nTotDstar=0,nTot3Prong=0;
  Int_t nTotD0toKpi=0;
  //  const  Int_t nptbins=10;
  Double_t standardptbins[14]={0.,0.5,1.,2.,3.,5.,7.,10.,15.,20.,30.,40.,50.,100.};
  Double_t pD[3],xD[3],pXtrTrue[2],pYtrTrue[2],pZtrTrue[2],xtr1[3],xtr2[3];
  Double_t cutsD0[9]=
	// cutsD0[0] = inv. mass half width [GeV]
	// cutsD0[1] = dca [cm]
	// cutsD0[2] = cosThetaStar
	// cutsD0[3] = pTK [GeV/c]
	// cutsD0[4] = pTPi [GeV/c]
	// cutsD0[5] = d0K [cm]   upper limit!
	// cutsD0[6] = d0Pi [cm]  upper limit!
	// cutsD0[7] = d0d0 [cm^2]
	// cutsD0[8] = cosThetaPoint
    {1000.,
     100000.,
     1.1,
     0.,
     0.,
     100000.,
     100000.,
     100000000.,
     -1.1};   
  if(fD0usecuts){
    cutsD0[0]=0.033;
    cutsD0[1]=0.03;
    cutsD0[2]=0.8;
    cutsD0[3]=0.;
    cutsD0[4]=0.;
    cutsD0[5]=100000.;
    cutsD0[6]=100000.;
    cutsD0[7]=-0.00005;
    cutsD0[8]=0.8; 
  }
  const Double_t fixedInvMassCut=cutsD0[0];
  AliAODVertex *vtx1=0;
  // primary vertex
  vtx1 = (AliAODVertex*)fAOD->GetPrimaryVertex();
  Double_t vtxTrue[3];
  fAODmcHeader->GetVertex(vtxTrue);
  //	vtx1->Print();
  AliAODRecoDecayHF *aodDMC=0x0;
  Bool_t nomum=kFALSE,background=kFALSE;
  AliAODMCParticle *mum1=0x0,*mum2=0x0;
  AliAODMCParticle *b1=0x0,*b2=0x0;
  AliAODMCParticle *sonpart=0x0,*grandmoth1=0x0;
  
  // make trkIDtoEntry register (temporary)
  Int_t trkIDtoEntry[100000];
  for(Int_t it=0;it<fAOD->GetNumberOfTracks();it++) {
    AliAODTrack *track = fAOD->GetTrack(it);
    if(track->GetID()<0) {
      //printf("Track ID <0, id= %d\n",track->GetID());
      return;
    }
    trkIDtoEntry[track->GetID()]=it;
  }
  
  
  // loop over D0->Kpi candidates
  Int_t nD0toKpi = fArrayD0toKpi->GetEntriesFast();
  nTotD0toKpi += nD0toKpi;
  //	cout<<"Number of D0->Kpi: "<<nD0toKpi<<endl;
  
  for (Int_t iD0toKpi = 0; iD0toKpi < nD0toKpi; iD0toKpi++) {
    if(aodDMC!=0x0)delete aodDMC;
    mum1=0x0;
    mum2=0x0;
    b1=0x0;
    b2=0x0;
    sonpart=0x0;
    grandmoth1=0x0;
    
    AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)fArrayD0toKpi->UncheckedAt(iD0toKpi);
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()) {
      d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
      unsetvtx=kTRUE;
    }
    
    //    d->SetOwnPrimaryVtx(vtx1); 
    //unsetvtx=kTRUE;
    
    Int_t okD0=0,okD0bar=0; 
    Double_t invMassD0,invMassD0bar,cutmassvalue=-1.;
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();
    nomum=kFALSE;
    background=kFALSE;

    cutsD0[0]=fixedInvMassCut;
    if(fSideBands){
      cutmassvalue=cutsD0[0];
      if(fSideBands==999.)cutsD0[0]=1000.;
      else cutsD0[0]=2.*TMath::Abs(fSideBands)*cutsD0[0];
    }
    if(d->SelectD0(cutsD0,okD0,okD0bar)) {
      if(fSideBands){
	if(fSideBands==999.)cutmassvalue=-1.;
	d->InvMassD0(invMassD0,invMassD0bar);
	
	if(TMath::Abs(invMassD0-mD0PDG)    > TMath::Abs(fSideBands)*cutmassvalue) okD0 = okD0 && kTRUE;
	else okD0=kFALSE;
	if(TMath::Abs(invMassD0bar-mD0PDG) > TMath::Abs(fSideBands)*cutmassvalue) okD0bar = okD0bar && kTRUE;    // SIDE BANDS SELECTION
	else okD0bar=kFALSE;

	if(fSideBands<0.){
	  if(!okD0 && !okD0bar) continue;
	}
	else if (!okD0 || !okD0bar)continue;
      }
      
      
      // get daughter AOD tracks
      AliAODTrack *trk0 = (AliAODTrack*)d->GetDaughter(0);
      AliAODTrack *trk1 = (AliAODTrack*)d->GetDaughter(1);
      
      //Check if it's a True D0
      if(fcheckMC){
	while(!background){
	  
	  if(trk0->GetLabel()==-1||trk1->GetLabel()==-1){
	    //fake tracks!!
	    background=kTRUE;  //FAKE TRACK
	    break;
	  }
	  //      printf("Before entering the MC checks \n");
	  
	  b1=(AliAODMCParticle*)fArrayMC->At(trk0->GetLabel());
	  b2=(AliAODMCParticle*)fArrayMC->At(trk1->GetLabel());
	  
	  
	  if(b1->GetMother()==-1||b2->GetMother()==-1){
	    nomum=kTRUE; //Tracks with no mother 
	    //printf("Inside problem mothering \n");// FAKE DECAY VERTEX
	    background=kTRUE;
	    break;
	    // if(!fStudyPureBackground)continue;
	  }
       
	  if(b1->GetMother()!=b2->GetMother()){
	    //Check the label of the mother is the same
	    background=kTRUE; // NOT SAME MOTHER
	    break;
	  }
	  //	  printf("After first mothering \n");
	  mum1=(AliAODMCParticle*)fArrayMC->At(b1->GetMother());
	  mum2=(AliAODMCParticle*)fArrayMC->At(b2->GetMother());
	  //	if((mum1->GetPdgCode()!=mum2->GetPdgCode()))continue; //Check the mother is the same particle
	  // printf("After second mothering \n");
	  //      printf("Particle codes: tr1: %d, tr2: %d, mum1: %d, mum 2: %d \n",b1->GetPdgCode(),b2->GetPdgCode(),mum1->GetPdgCode(),mum2->GetPdgCode());
	  if(fcheckMCD0)if(TMath::Abs(mum1->GetPdgCode())!=421){
	    //^fStudyPureBackground)continue;//NPRONG
	    background=kTRUE;// NOT A D0
	    break;
	  }
	  if(fcheckMC2prongs)if(mum1->GetDaughter(1)-mum1->GetDaughter(0)+1!=2){
	    //^fStudyPureBackground)continue;// NPRONG
	    background=kTRUE; // NOT A 2 PRONG DECAY
	    break;
	  }
	  
	  if(mum1->GetMother()==-1){
	    background=kTRUE; // A particle coming from nothing
	    break;
	  }
	      
	  // ^fStudyPureBackground)continue;
	  //	  printf("After second mothering \n");
	  grandmoth1=(AliAODMCParticle*)fArrayMC->At(mum1->GetMother());
	  if(fSkipD0star)if(TMath::Abs(grandmoth1->GetPdgCode())==413||TMath::Abs(grandmoth1->GetPdgCode())==423){
	    //^fStudyPureBackground)continue;
	    background=kTRUE;// D0 COMING FROM A D0*
	    break;
	  }
	  if(fcheckMCD0){//THIS CHECK IS NEEDE TO AVOID POSSIBLE FAILURE IN THE SECOND WHILE
	    while(TMath::Abs(grandmoth1->GetPdgCode())!=4&&TMath::Abs(grandmoth1->GetPdgCode())!=5){
	      if(grandmoth1->GetMother()==-1){
		//	      printf("Inside grandmothering \n");
		//printf("mother=-1, pdgcode: %d \n",grandmoth1->GetPdgCode());
		Int_t son=grandmoth1->GetDaughter(0);
		sonpart=(AliAODMCParticle*)fArrayMC->At(son);
		while(TMath::Abs(sonpart->GetPdgCode())!=421){
		  //printf("mother=-1, pdgcode: %d \n",sonpart->GetPdgCode());
		  son++;
		  sonpart=(AliAODMCParticle*)fArrayMC->At(son);
		}
		break;
	      }
	      grandmoth1=(AliAODMCParticle*)fArrayMC->At(grandmoth1->GetMother());
	    }
	  }
	  if(fcheckMCprompt)if(TMath::Abs(grandmoth1->GetPdgCode())!=4){
	    //^fStudyPureBackground)continue;
	    background=kTRUE;
	    break;
	  }
	  if(fcheckMCfromB){
	    if(TMath::Abs(grandmoth1->GetPdgCode())!=5){
	      //^fStudyPureBackground)continue;
	      background=kTRUE;
	      break;
	    }
	  }
	  //else if(!fStudyPureBackground)continue;
	  break;
	}
	if(background^fStudyPureBackground)continue;
	if(fStudyd0fromBTrue){
	  if(b1!=0x0&&b2!=0x0){
	    Int_t charge[2]={0,0};
	    if(b1->Charge()==-1)charge[0]=1;
	    else {
	      if(b2->Charge()==-1){
		//printf("Same charges for prongs \n");
		continue;
	      }
	      charge[1]=1;
	    }
	    /* if(!pXtrTrue[charge[0]]=b1->Px()){printf("!first MCtrack:Get momentum X failed \n");continue;}
	       if(!b1->Py(pYtrTrue[charge[0]])){printf("!first MCtrack:Get momentum Y failed \n");continue;}
	       if(!b1->Pz(pZtrTrue[charge[0]])){printf("!first MCtrack:Get momentum Z failed \n");continue;}
	    */
	    
	    pXtrTrue[charge[0]]=b1->Px();
	    pYtrTrue[charge[0]]=b1->Py();
	    pZtrTrue[charge[0]]=b1->Pz();
	    if(!b1->XvYvZv(xtr1)){
	      //printf("!first MCtrack:Get position failed \n");
	      continue;
	    }
	    /*if(!b2->Px(pXtrTrue[charge[1]])){printf("!second MCtrack:Get momentum X failed \n");continue;}
	      if(!b2->Py(pYtrTrue[charge[1]])){printf("!second MCtrack:Get momentum Y failed \n");continue;}
	      if(!b2->Pz(pZtrTrue[charge[1]])){printf("!second MCtrack:Get momentum Z failed \n");continue;}*/
	    pXtrTrue[charge[1]]=b2->Px();
	    pYtrTrue[charge[1]]=b2->Py();
	    pZtrTrue[charge[1]]=b2->Pz();
	    if(!nomum){
	      if(!(mum1==0x0||mum2==0x0)){
		if(!mum1->PxPyPz(pD)){
		  //printf("!D from B:Get momentum failed \n");
		  continue;
		}
		if(!mum1->XvYvZv(xD)){
		  //printf("!D from B:Get position failed \n");
		  continue;
		}
		if(pXtrTrue[0]+pXtrTrue[1]!=pD[0]){
		  //printf("Warning! MCinfo: Pt of the mother doesn't match the sum of the daughters pt: X component  %f + %f = %f Vs %f \n",pXtrTrue[0],pXtrTrue[1],pXtrTrue[0]+pXtrTrue[1],pD[0]);
		  //printf("Warning! MCinfo: Pt of the mother doesn't match the sum of the daughters pt: Y componente %f + %f = %f Vs %f \n",pYtrTrue[0],pYtrTrue[1],pYtrTrue[0]+pYtrTrue[1],pD[1]);
		  //printf("Warning!: Pt comparison from the sum of the daughters:%f Vs mother pT: %f \n",TMath::Sqrt((pXtrTrue[0]+pXtrTrue[1])*(pXtrTrue[0]+pXtrTrue[1])+(pYtrTrue[0]+pYtrTrue[1])*(pYtrTrue[0]+pYtrTrue[1])),TMath::Sqrt(pD[0]*pD[0]+pD[1]*pD[1]));
		  //	      return;
		}
		
		//	    if(!b2->PxPyPz(ptr2True)){printf("!second MCtrack:Get momentum failed \n");continue;}
		if(!b2->XvYvZv(xtr2)){
		  //printf("!second MCtrack:Get position failed \n");
		  continue;
		}
		Double_t d0dummy[2]={0.,0.};//TEMPORARY : dummy d0 for AliAODRecoDeay constructor
		aodDMC=new AliAODRecoDecayHF(vtxTrue,xD,2,0,pXtrTrue,pYtrTrue,pZtrTrue,d0dummy);
		fhMassTrue->Fill(mum1->GetCalcMass());
	      }
	      //else printf("Warning! The candidate comes from two track with a different mother, \n -> the true mother d0 doesn't exist ! -> There will be different entries between observed and MC d0distributions \n");
	    }
	  }
	}
	//printf("End of MC checks \n");
      }
     
      //printf("A real D0 was found!! \n");
      
    
      //cout<<1e8*d->Prodd0d0()<<endl;
      /*      if(fSideBands){	// The integral of the Hists will in general differ 
	                // from the number of candidates which passed the selection
	if(okD0)fhMass->Fill(d->InvMassD0(),1.);
	if(okD0bar)fhMass->Fill(d->InvMassD0bar(),1.);
      }
      else { // The integral of the Hists will in general differ 
	// from the number of candidates which passed the selection
	if(okD0)fhMass->Fill(d->InvMassD0(),1.);
	if(okD0bar)fhMass->Fill(d->InvMassD0bar(),1.);
      }
      */
      
      if(okD0)fhMass->Fill(d->InvMassD0(),1.);
      if(okD0bar)fhMass->Fill(d->InvMassD0bar(),1.);
   
      fhCPtaVSd0d0->Fill(1e8*d->Prodd0d0(),d->CosPointingAngle());
      fhSecVtxZ->Fill(d->GetSecVtxZ());
      fhSecVtxX->Fill(d->GetSecVtxX()*10000.);
      fhSecVtxY->Fill(d->GetSecVtxY()*10000.);
      fhSecVtxXY->Fill(d->GetSecVtxX()*10000.,d->GetSecVtxY()*10000.);
      fhSecVtxPhi->Fill(TMath::ATan2(d->GetSecVtxY()*10000.,d->GetSecVtxX()*10000.)*TMath::RadToDeg());
      fhd0xd0->Fill(1e8*d->Prodd0d0());
      fhCPta->Fill(d->CosPointingAngle());
      
    
      //cout<<d->GetSecVtxX() <<endl;
      
      if(!trk0 || !trk1) {
	/*	if(d->GetProngID(0)<0||d->GetProngID(1)<0){
		printf("Wrong ProngID,%d or %d\n",d->GetProngID(0),d->GetProngID(1));continue;
		}
	*/
	trk0=fAOD->GetTrack(trkIDtoEntry[d->GetProngID(0)]);
	trk1=fAOD->GetTrack(trkIDtoEntry[d->GetProngID(1)]);
	//	      cout<<"references to standard AOD not available"<<endl;
      }
      //	    cout<<"pt of positive track: "<<trk0->Pt()<<endl;
      
      // make a AliExternalTrackParam from the D0 
      // and calculate impact parameters to primary vertex
      
      
      //	AliExternalTrackParam trackD0(d);
      //cout<<"pt of D0: "<<d->Pt()<<" (AliAODRecoDecay); "<<trackD0.Pt()<<" (track)"<<endl;
      //trackD0.Print();
      Double_t ptD0=d->Pt();
      for(Int_t j=0;j<fnbins;j++){
	if(ptD0<standardptbins[j+1]){
	  fhd0D0VtxTruept[j]->Fill(d->AliAODRecoDecay::ImpParXY(vtxTrue)*10000.);
	  fhd0D0pt[j]->Fill(d->ImpParXY()*10000.);
	  // printf("HERE in tha class \n");
	  if(aodDMC!=0x0)fhMCd0D0pt[j]->Fill(aodDMC->ImpParXY()*10000.);
	  // printf("HERE2 in tha class \n");
	  break;
	}
      }
      fhd0D0->Fill(d->ImpParXY()*10000.);	   
      fhd0D0VtxTrue->Fill(d->AliAODRecoDecay::ImpParXY(vtxTrue)*10000.);	   
      if(aodDMC!=0x0)fhMCd0D0->Fill(aodDMC->ImpParXY()*10000.);
      //  if((1.864-3.*0.011.<d->InvMassD0()&&d->InvMassD0()<1.864+3.*0.011)||(1.864-3.*0.011<d->InvMassD0bar()&&d->InvMassD0bar()<1.864+3.*0.011)){
      if(aodDMC!=0x0){
	delete aodDMC;
	aodDMC=0x0;
      }
      mum1=0x0;
      mum2=0x0;
      b1=0x0;
      b2=0x0;
      sonpart=0x0;
      grandmoth1=0x0;
    
      //	    }
      //	Double_t dz[2],covdz[3];
      //trackD0.PropagateToDCA(vtx1,fAOD->GetMagneticField(),1000.,dz,covdz);
      //	cout<<"D0 impact parameter rphi: "<<dz[0]<<" +- "<<TMath::Sqrt(covdz[0])<<endl;
    }
    // printf("Before end of the event \n");
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
	  
  }

  
  
  // Post output data.
  
  PostData(0,fhCPtaVSd0d0);
  PostData(1,fhSecVtxXY);
  PostData(2,fhd0xd0);
  PostData(3,fhCPta);
  PostData(4,fhSecVtxZ);
  PostData(5,fhSecVtxX);
  PostData(6,fhSecVtxY);
  PostData(7,fhSecVtxPhi);

  
  for(Int_t j=0;j<fnbins;j++){
    
    PostData(j+8, fhd0D0pt[j]);
  }
  
  PostData(fnbins+8, fhd0D0);

  for(Int_t j=0;j<fnbins;j++){
    
    PostData(j+fnbins+9,fhMCd0D0pt[j]);
  }
 
  PostData(2*fnbins+9, fhMCd0D0);

  for(Int_t j=0;j<fnbins;j++){
    
    PostData(j+2*fnbins+10,fhd0D0VtxTruept[j]);
  }
 
  PostData(3*fnbins+10,fhd0D0VtxTrue);
  PostData(3*fnbins+11,fhMassTrue);
  PostData(3*fnbins+12,fhMass);
}

//________________________________________________________________________
void AliAnalysisTaskCharmFraction::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  fhd0D0 = dynamic_cast<TH1F*> (GetOutputData(fnbins+8));
  if (!fhd0D0) {
    Printf("ERROR: fhd0D0 not available");
    return;
  }
   
  TCanvas *c1 = new TCanvas("AliAnalysisTaskCharmFraction","D0 impact par",10,10,510,510);
  c1->cd(1)->SetLogy();
  fhd0D0->DrawCopy("E");

  TCanvas *cd0D0pt=new TCanvas("cd0D0pt","cd0D0pt",0,0,800,700);
  cd0D0pt->Divide(fnbins/2+fnbins%2,2);
  for(Int_t j=0;j<fnbins;j++){
    fhd0D0pt[j] = dynamic_cast<TH1F*> (GetOutputData(8+j));
    if (!fhd0D0pt[j]) {
      Printf("ERROR: fhd0D0pt not available");
      return;
    }
    cd0D0pt->cd(j+1);
    fhd0D0pt[j]->Draw();
  }
  
}
