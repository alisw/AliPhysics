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

///////////////////////////////////////////////////////////////////////////
// AliAnalysisTaskSE for the extraction of the various histograms to
// study the pt spectra of identified hadrons:
// - log(dEdx)-log(dEdxBB) distributions for pions, kaons and protons in pt bins
// - Pt distributions of pions, kaons and protons with nSigma PID
// Authors: 
// E. Biolcati, biolcati@to.infn.it
// L. Milano, milano@to.infn.it
// F. Prino, prino@to.infn.it
///////////////////////////////////////////////////////////////////////////

#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TH1F.h>
#include <TRandom3.h>
#include <TH2F.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TParticle.h>
#include <Rtypes.h>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliPhysicsSelection.h"
#include "AliAnalysisTaskSEITSsaSpectra.h"

ClassImp(AliAnalysisTaskSEITSsaSpectra)

//________________________________________________________________________
AliAnalysisTaskSEITSsaSpectra::AliAnalysisTaskSEITSsaSpectra():
AliAnalysisTaskSE("Task CFits"),
  fESD(0),
  fOutput(0),
  fHistNEvents(0),
  fHistNTracks(0),
  fHistDEDX(0),
  fHistDEDXdouble(0),
  fHistBeforeEvSel(0),
  fHistAfterEvSel(0),
  fMinSPDPts(1),
  fMinNdEdxSamples(3),
  fMindEdx(0.),
  fMinNSigma(3.),
  fMaxY(0.5),
  fMaxChi2Clu(1.),
  fNSigmaDCAxy(7.),
  fNSigmaDCAz(7.),
  fMC(kFALSE), 
  fSmearMC(kFALSE),
  fSmearP(0.),
  fSmeardEdx(0.),
  fRandGener(0),
  fFillNtuple(kFALSE),
  fNtupleNSigma(0),
  fNtupleMC(0)
{
  // Constructor
  Double_t xbins[kNbins+1]={0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.0};
  for(Int_t iBin=0; iBin<kNbins+1; iBin++) fPtBinLimits[iBin]=xbins[iBin];
  fRandGener=new TRandom3(0);

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  Printf("end of AliAnalysisTaskSEITSsaSpectra");
}

//___________________________________________________________________________
AliAnalysisTaskSEITSsaSpectra::~AliAnalysisTaskSEITSsaSpectra(){
  // Destructor
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
  if(fRandGener) delete fRandGener;
}

//________________________________________________________________________
Double_t AliAnalysisTaskSEITSsaSpectra::CookdEdx(Double_t *s){
  // truncated mean for the dEdx
  Int_t nc=0; 
  Double_t dedx[4]={0.,0.,0.,0.};
  for (Int_t il=0; il<4; il++) { // count good (>0) dE/dx values
    if(s[il]>fMindEdx){
      dedx[nc]= s[il];
      nc++;
    }
  }
  if(nc<fMinNdEdxSamples) return -1.;
  
  Double_t tmp;
  Int_t swap; // sort in ascending order 
  do {
    swap=0;
    for (Int_t i=0; i<nc-1; i++) {
      if (dedx[i]<=dedx[i+1]) continue;
      tmp=dedx[i];
      dedx[i]=dedx[i+1];
      dedx[i+1]=tmp;
      swap++;
    } 
  } while (swap);
  
  Double_t sumamp=0,sumweight=0;
  Double_t weight[4]={1.,1.,0.,0.};
  if(nc==3) weight[1]=0.5;
  else if(nc<3) weight[1]=0.;
  for (Int_t i=0; i<nc; i++) {
    sumamp+= dedx[i]*weight[i];
    sumweight+=weight[i];
  }
  return sumamp/sumweight;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskSEITSsaSpectra::DCAcut(Double_t impactXY, Double_t impactZ, Double_t pt, Bool_t optMC){
  // cut on transverse impact parameter updaated on 20-5-2010
  // from the study of L. Milano, F. Prino on the ITS standalone tracks
  // using the common binning of the TPC tracks
  Double_t xyP[3];
  Double_t zP[3];
  if(optMC){
    xyP[0]=88.63;//SIMpass4a12
    xyP[1]=19.57;
    xyP[2]=1.65;
    zP[0]=140.98;
    zP[1]=62.33;
    zP[2]=1.15;
  }
  else{
    xyP[0]=85.28;//DATApass6
    xyP[1]=25.78;
    xyP[2]=1.55;
    zP[0]=146.80;
    zP[1]=70.07;
    zP[2]=1.11;
  }
  Double_t xySigma = xyP[0] + xyP[1]/TMath::Power(TMath::Abs(pt),xyP[2]);
  Double_t xyMax = fNSigmaDCAxy*xySigma; //in micron
  if((TMath::Abs(impactXY)*10000)>xyMax) return kFALSE;
  
  Double_t zSigma = zP[0] + zP[1]/TMath::Power(TMath::Abs(pt),zP[2]);
  Double_t zMax = fNSigmaDCAz*zSigma; //in micron
  if((TMath::Abs(impactZ)*10000)>zMax) return kFALSE;
  
  return kTRUE;
}

//________________________________________________________________________
Double_t AliAnalysisTaskSEITSsaSpectra::Eta2y(Double_t pt, Double_t m, Double_t eta){
  Double_t mt = TMath::Sqrt(m*m + pt*pt);
  return TMath::ASinH(pt/mt*TMath::SinH(eta));
}


//________________________________________________________________________
Double_t AliAnalysisTaskSEITSsaSpectra::BetheBloch(Double_t bg,Bool_t optMC) {
  // BB PHOBOS parametrization tuned by G. Ortona on 900 GeV pp data of 2009
  Double_t par[5];
  if(optMC){//Double_t par[5]={139.1,23.36,0.06052,0.2043,-0.0004999};
    par[0]=139.1;
    par[1]=23.36;
    par[2]=0.06052;
    par[3]=0.2043;
    par[4]=-0.0004999; 
  }else {
    //Double_t par[5]={5.33458e+04,1.65303e+01,2.60065e-03,3.59533e-04,7.51168e-05};}
  par[0]=5.33458e+04;
  par[1]=1.65303e+01;
  par[2]=2.60065e-03;
  par[3]=3.59533e-04;
  par[4]=7.51168e-05;
  }
  Double_t beta = bg/TMath::Sqrt(1.+ bg*bg);
  Double_t gamma=bg/beta;
  Double_t eff=1.0;
  if(bg<par[2]) eff=(bg-par[3])*(bg-par[3])+par[4];
  else eff=(par[2]-par[3])*(par[2]-par[3])+par[4];
  return (par[1]+2.0*TMath::Log(gamma)-beta*beta)*(par[0]/(beta*beta))*eff;
}


//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::UserCreateOutputObjects(){
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("Spiderman");
  
  fHistNEvents = new TH1F("fHistNEvents", "Number of processed events",6,0.5,6.5);
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);
  
  fHistNTracks = new TH1F("fHistNTracks", "Number of ITSsa tracks",20,0.5,20.5);
  fHistNTracks->Sumw2();
  fHistNTracks->SetMinimum(0);
  fOutput->Add(fHistNTracks);
  
  //binning for the histogram
  const Int_t hnbins=400;
  Double_t hxmin = 0.01;
  Double_t hxmax = 10;
  Double_t hlogxmin = TMath::Log10(hxmin);
  Double_t hlogxmax = TMath::Log10(hxmax);
  Double_t hbinwidth = (hlogxmax-hlogxmin)/hnbins;
  Double_t hxbins[hnbins+1];
  hxbins[0] = 0.01; 
  for (Int_t i=1;i<=hnbins;i++) {
    hxbins[i] = hxmin + TMath::Power(10,hlogxmin+i*hbinwidth);
  }
  
  fHistDEDX = new TH2F("fHistDEDX","",hnbins,hxbins,900,0,1000);
  fOutput->Add(fHistDEDX);
  
  fHistDEDXdouble = new TH2F("fHistDEDXdouble","",500,-5,5,900,0,1000);
  fOutput->Add(fHistDEDXdouble);
  

  fHistBeforeEvSel = new TH1F("fHistBeforeEvSel","fHistBeforeEvSel",kNbins,fPtBinLimits);
  fHistAfterEvSel = new TH1F("fHistAfterEvSel","fHistAfterEvSel",kNbins,fPtBinLimits);
  fOutput->Add(fHistBeforeEvSel);
  fOutput->Add(fHistAfterEvSel);
  
  
  
  for(Int_t j=0;j<3;j++){
    fHistMCpos[j] = new TH1F(Form("fHistMCpos%d",j),Form("fHistMCpos%d",j),kNbins,fPtBinLimits);
    fHistMCneg[j] = new TH1F(Form("fHistMCneg%d",j),Form("fHistMCneg%d",j),kNbins,fPtBinLimits);
    fOutput->Add(fHistMCneg[j]);
		fOutput->Add(fHistMCpos[j]);
  }
  

  for(Int_t j=0;j<3;j++){
    fHistMCposBefEvSel[j] = new TH1F(Form("fHistMCposBefEvSel%d",j),Form("fHistMCposBefEvSel%d",j),kNbins,fPtBinLimits);
    fHistMCnegBefEvSel[j] = new TH1F(Form("fHistMCnegBefEvSel%d",j),Form("fHistMCnegBefEvSel%d",j),kNbins,fPtBinLimits);
    fOutput->Add(fHistMCnegBefEvSel[j]);
    fOutput->Add(fHistMCposBefEvSel[j]);
  }
  
  for(Int_t i=0; i<4; i++){
    fHistCharge[i] = new TH1F(Form("fHistChargeLay%d",i),Form("fHistChargeLay%d",i),100,0,300);
    fOutput->Add(fHistCharge[i]);
  }
  
  for(Int_t i=0; i<kNbins; i++){
    fHistPosPi[i] = new TH1F(Form("fHistPosPi%d",i),Form("fHistPosPi%d",i),175,-3.5,3.5);	
    fHistPosK[i]  = new TH1F(Form("fHistPosK%d",i),Form("fHistPosK%d",i),175,-3.5,3.5);	
    fHistPosP[i]  = new TH1F(Form("fHistPosP%d",i),Form("fHistPosP%d",i),175,-3.5,3.5);	
    fHistNegPi[i] = new TH1F(Form("fHistNegPi%d",i),Form("fHistNegPi%d",i),175,-3.5,3.5);	
    fHistNegK[i]  = new TH1F(Form("fHistNegK%d",i),Form("fHistNegK%d",i),175,-3.5,3.5);	
    fHistNegP[i]  = new TH1F(Form("fHistNegP%d",i),Form("fHistNegP%d",i),175,-3.5,3.5);	
    
    fHistDCAPosPi[i] = new TH1F(Form("fHistDCAPosPi%d",i),Form("fHistDCAPosPi%d",i),2000,-1,1);  //DCA distr.	
    fHistDCAPosK[i]  = new TH1F(Form("fHistDCAPosK%d",i),Form("fHistDCAPosK%d",i),2000,-1,1);	
    fHistDCAPosP[i]  = new TH1F(Form("fHistDCAPosP%d",i),Form("fHistDCAPosP%d",i),2000,-1,1);	
    fHistDCANegPi[i] = new TH1F(Form("fHistDCANegPi%d",i),Form("fHistDCANegPi%d",i),2000,-1,1);	
    fHistDCANegK[i]  = new TH1F(Form("fHistDCANegK%d",i),Form("fHistDCANegK%d",i),2000,-1,1);	
    fHistDCANegP[i]  = new TH1F(Form("fHistDCANegP%d",i),Form("fHistDCANegP%d",i),2000,-1,1);	
    
    fHistMCPosPi[i] = new TH1F(Form("fHistMCPosPi%d",i),Form("fHistMCPosPi%d",i),175,-3.5,3.5);	//MC truth
    fHistMCPosK[i]  = new TH1F(Form("fHistMCPosK%d",i),Form("fHistMCPosK%d",i),175,-3.5,3.5);	
    fHistMCPosP[i]  = new TH1F(Form("fHistMCPosP%d",i),Form("fHistMCPosP%d",i),175,-3.5,3.5);	
    fHistMCNegPi[i] = new TH1F(Form("fHistMCNegPi%d",i),Form("fHistMCNegPi%d",i),175,-3.5,3.5);	
    fHistMCNegK[i]  = new TH1F(Form("fHistMCNegK%d",i),Form("fHistMCNegK%d",i),175,-3.5,3.5);	
    fHistMCNegP[i]  = new TH1F(Form("fHistMCNegP%d",i),Form("fHistMCNegP%d",i),175,-3.5,3.5);	
    fOutput->Add(fHistPosPi[i]);
    fOutput->Add(fHistPosK[i]);
    fOutput->Add(fHistPosP[i]);
    fOutput->Add(fHistNegPi[i]);
    fOutput->Add(fHistNegK[i]);
    fOutput->Add(fHistNegP[i]);
    
    fOutput->Add(fHistDCAPosPi[i]);//DCA distr.
    fOutput->Add(fHistDCAPosK[i]);
    fOutput->Add(fHistDCAPosP[i]);
    fOutput->Add(fHistDCANegPi[i]);
    fOutput->Add(fHistDCANegK[i]);
    fOutput->Add(fHistDCANegP[i]);
    
    fOutput->Add(fHistMCPosPi[i]);//MC truth
    fOutput->Add(fHistMCPosK[i]);
    fOutput->Add(fHistMCPosP[i]);
    fOutput->Add(fHistMCNegPi[i]);
    fOutput->Add(fHistMCNegK[i]);
    fOutput->Add(fHistMCNegP[i]);
  }
  
  //NSigma Histos
  for(Int_t j=0;j<3;j++){
    fHistPosNSigma[j] = new TH1F(Form("hHistPosNSigma%d",j),Form("hHistPosNSigma%d",j),kNbins,fPtBinLimits);
    fHistNegNSigma[j] = new TH1F(Form("hHistNegNSigma%d",j),Form("hHistNegNSigma%d",j),kNbins,fPtBinLimits);
    fHistPosNSigmaPrim[j] = new TH1F(Form("hHistPosNSigmaPrim%d",j),Form("hHistPosNSigmaPrim%d",j),kNbins,fPtBinLimits);
    fHistNegNSigmaPrim[j] = new TH1F(Form("hHistNegNSigmaPrim%d",j),Form("hHistNegNSigmaPrim%d",j),kNbins,fPtBinLimits);
    fHistPosNSigmaPrimMC[j] = new TH1F(Form("hHistPosNSigmaPrimMC%d",j),Form("hHistPosNSigmaPrimMC%d",j),kNbins,fPtBinLimits);
    fHistNegNSigmaPrimMC[j] = new TH1F(Form("hHistNegNSigmaPrimMC%d",j),Form("hHistNegNSigmaPrimMC%d",j),kNbins,fPtBinLimits);
    fOutput->Add(fHistPosNSigma[j]);
    fOutput->Add(fHistNegNSigma[j]);
    fOutput->Add(fHistPosNSigmaPrim[j]);
    fOutput->Add(fHistNegNSigmaPrim[j]);
    fOutput->Add(fHistPosNSigmaPrimMC[j]);
    fOutput->Add(fHistNegNSigmaPrimMC[j]);
  }
  
  fNtupleNSigma = new TNtuple("fNtupleNSigma","fNtupleNSigma","p:pt:dedx:ncls:sign:run:eta:impactXY:impactZ:isph:pdgcode:mfl");
  fOutput->Add(fNtupleNSigma);
  fNtupleMC = new TNtuple("fNtupleMC","fNtupleMC","ptMC:pdgcode:signMC:etaMC:yMC:isph:evSel:run");
  fOutput->Add(fNtupleMC);
  
  Printf("end of CreateOutputObjects");
}

//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::UserExec(Option_t *){
  
  fESD=(AliESDEvent*)InputEvent();
  if(!fESD) {
    printf("AliAnalysisTaskSDDRP::Exec(): bad ESD\n");
    return;
  } 
  
  //binning for the dEdx distributions

  //variables
  Float_t pdgmass[3]={0.13957,0.493677,0.938272}; //mass for pi, K, P (Gev/c^2)
  Int_t listcode[3]={211,321,2212};//code for pi, K, P (Gev/c^2)
  Double_t s[4];
  Float_t ptMC=-999,yMC=-999,code=-999, signMC=-999,isph=-999,mfl=-999.;
  Float_t impactXY=-999, impactZ=-999;
  Int_t evSel=1;
  AliESDtrack* track;
  UInt_t status; 
  AliStack* stack=0;
  TParticle *part=0;
  TParticlePDG *pdgPart=0;
	
  //Nsigma Method
  Float_t resodedx[4];
  if(fMC){
    resodedx[0]=0.13;
    resodedx[1]=0.13;
    resodedx[2]=0.134;
    resodedx[3]=0.127;
  }else{
    resodedx[0]=0.23;
    resodedx[1]=0.18;
    resodedx[2]=0.16;
    resodedx[3]=0.14;
  }
  /////////////////////
  if(fMC){
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
      printf("ERROR: stack not available\n");
      return;
    }
  }
  //flags for MC
  Int_t nTrackMC=0; 
  if(stack) nTrackMC = stack->GetNtrack();	
  const AliESDVertex *vtx =  fESD->GetPrimaryVertexSPD();
  
  //event selection
  fHistNEvents->Fill(1);
  if(!vtx)evSel=0;
  else{
    fHistNEvents->Fill(2);
    if(vtx->GetNContributors()<0) evSel=0;
    else{
      fHistNEvents->Fill(3);
      if(TMath::Abs(vtx->GetZv())>10) evSel=0;
      else{
	fHistNEvents->Fill(4);
	if(vtx->GetZRes()>0.5) evSel=0;
	else{
	  fHistNEvents->Fill(5);
	  if(vtx->IsFromVertexerZ() && vtx->GetDispersion()>0.03) evSel=0;
	  else fHistNEvents->Fill(6);
	}
      }
    }
  }
	
  /////first loop on stack, before event selection, filling MC ntuple
	
  for(Int_t imc=0; imc<nTrackMC; imc++){
    part = stack->Particle(imc);
    if(!stack->IsPhysicalPrimary(imc))continue;//no secondary in the MC sample 
    isph=1.;
    pdgPart = part->GetPDG();
    if(pdgPart->Charge()==0) continue; //no neutral particles
    if(TMath::Abs(part->Eta()) > 0.9) continue; //pseudorapidity-acceptance cut
    if(part->Energy() != TMath::Abs(part->Pz())) yMC = 0.5*TMath::Log((part->Energy()+part->Pz())/(part->Energy()-part->Pz()));
    if(TMath::Abs(yMC) > fMaxY) continue; //rapidity cut
    
    if(pdgPart->Charge()>0) signMC=1;
    else signMC=-1;
    ptMC=part->Pt();
    code=pdgPart->PdgCode();
    
	  
    //filling MC ntuple
    if(TMath::Abs(code)==211 || TMath::Abs(code)==321 || TMath::Abs(code)==2212){
      Float_t xntMC[8];
      Int_t indexMC=0;
      xntMC[indexMC++]=(Float_t)ptMC;
      xntMC[indexMC++]=(Float_t)code;
      xntMC[indexMC++]=(Float_t)signMC;
      xntMC[indexMC++]=(Float_t)part->Eta();
      xntMC[indexMC++]=(Float_t)yMC;
      xntMC[indexMC++]=(Float_t)isph;
      xntMC[indexMC++]=(Float_t)evSel;
      xntMC[indexMC++]=(Float_t)fESD->GetRunNumber();

      if(fFillNtuple) fNtupleMC->Fill(xntMC);
    }
    
    for(Int_t j=0; j<3; j++){
      if(TMath::Abs(code)==listcode[j]){
	if(signMC>0) fHistMCposBefEvSel[j]->Fill(TMath::Abs(ptMC));
	else  fHistMCnegBefEvSel[j]->Fill(TMath::Abs(ptMC));
      }
    }
    if(evSel==1){
      for(Int_t j=0; j<3; j++){
	if(TMath::Abs(code)==listcode[j]){
	  if(signMC>0) fHistMCpos[j]->Fill(TMath::Abs(ptMC));
	  else  fHistMCneg[j]->Fill(TMath::Abs(ptMC));
	}
      }
    }	
  }
  
  if(evSel==0)return;
	
  //loop on tracks
  for (Int_t iTrack=0; iTrack<fESD->GetNumberOfTracks(); iTrack++) {  
    isph=-999.;
    code=-999;
    mfl=-999;
	  
    track = (AliESDtrack*)fESD->GetTrack(iTrack);      
    if (!track) continue;
    
    //track selection
    fHistNTracks->Fill(1);
    status=track->GetStatus();
    if((status&AliESDtrack::kITSpureSA)==0) continue; //its standalone
    fHistNTracks->Fill(2);
    if((status&AliESDtrack::kITSrefit)==0) continue; //its refit
    fHistNTracks->Fill(3);
    if(track->GetSign()==0.) continue; //no neutral particles
    fHistNTracks->Fill(4);

	  
    //cluster in ITS
    UInt_t clumap = track->GetITSClusterMap();
    Int_t nSPD=0;
    for(Int_t il=0; il<2; il++) if(TESTBIT(clumap,il)) nSPD++;
    if(nSPD<fMinSPDPts) continue;
    fHistNTracks->Fill(5);
    Int_t count=0;
    for(Int_t j=2;j<6;j++) if(TESTBIT(clumap,j)) count++;
    if(count<fMinNdEdxSamples) continue; //at least 3 points on SSD/SDD
    fHistNTracks->Fill(6);
    //chisquare/nclusters	
    Int_t nclu=nSPD+count;
    if(track->GetITSchi2()/nclu > fMaxChi2Clu) continue; 
    fHistNTracks->Fill(7);
    //pseudorapidity and rapidity
    if(TMath::Abs(track->Eta()) > 0.9) continue;
    fHistNTracks->Fill(8);
    //truncated mean
    //if(fMC) for(Int_t j=0;j<2;j++) s[j]*=3.34/5.43;//correction for SDD miscalibration of the MCpass4
    track->GetITSdEdxSamples(s);
    Double_t dedx = CookdEdx(s);
    if(dedx<0) continue;
    fHistNTracks->Fill(9);


    Float_t pt = track->Pt();
    Int_t theBin=-1;
    for(Int_t m=0; m<kNbins; m++){
      if(TMath::Abs(pt) > fPtBinLimits[m] && TMath::Abs(pt) < fPtBinLimits[m+1]){
	theBin=m;
	break;
      }
    }
    track->GetImpactParameters(impactXY, impactZ);
	  
    //Filling Ntuple
    //information from the MC kinematics
    if(fMC){
      if(track->GetLabel()<0)isph=-1.;
      if(track->GetLabel()>=0){
	part = (TParticle*)stack->Particle(track->GetLabel());
	pdgPart = part->GetPDG();
	code = pdgPart->PdgCode();
	if(stack->IsPhysicalPrimary(track->GetLabel()))isph=1.;
	else{ 
	  isph=0.;
	  TParticle* moth = stack->Particle(part->GetFirstMother());
	  Float_t codemoth = TMath::Abs(moth->GetPdgCode());
	  mfl = Int_t (codemoth/ TMath::Power(10, Int_t(TMath::Log10(codemoth))));
	}
      }
    }
    Float_t xnt[12];
    Int_t index=0;
    xnt[index++]=(Float_t)track->GetP();
    xnt[index++]=(Float_t)track->Pt();
    xnt[index++]=(Float_t)dedx;
    xnt[index++]=(Float_t)count;
    xnt[index++]=(Float_t)track->GetSign();
    xnt[index++]=(Float_t)fESD->GetRunNumber();
    xnt[index++]=(Float_t)track->Eta();
    xnt[index++]=(Float_t)impactXY;
    xnt[index++]=(Float_t)impactZ;
    xnt[index++]=(Float_t)isph;
    xnt[index++]=(Float_t)code;
    xnt[index]=(Float_t)mfl;
	  
    if(fFillNtuple) fNtupleNSigma->Fill(xnt);
    
    
	
    //Compute y and bb
    Double_t y[3],bbtheo[3],logdiff[3];
    Float_t p=track->GetP();
    if(fMC && fSmearMC){
      dedx=fRandGener->Gaus(dedx,fSmeardEdx*dedx);
      p=fRandGener->Gaus(p,fSmearP*p);     
    }

    for(Int_t i=0;i<3;i++){
      y[i] = Eta2y(pt,pdgmass[i],track->Eta());
      bbtheo[i]=BetheBloch(p/pdgmass[i],fMC);
      logdiff[i]=TMath::Log(dedx) - TMath::Log(bbtheo[i]);
    }

    //NSigma Method
    Int_t resocls=(Int_t)count-1;

    Double_t nsigmas[3];
    Double_t min=999999.;
    Int_t minPos=-1;
    for(Int_t isp=0; isp<3; isp++){
      Double_t bb=bbtheo[isp];
      nsigmas[isp]=TMath::Abs((dedx-bb)/(resodedx[resocls]*bb));
      if(nsigmas[isp]<min){
	min=nsigmas[isp];
	minPos=isp;
      }
    }
    Double_t yPart=y[minPos];
    
    if(min<fMinNSigma && yPart<fMaxY){     
      //DCA distributions, before the DCA cuts
      if(theBin>=0 && theBin<kNbins){
	if(track->GetSign()>0){
	  if(minPos==0) fHistDCAPosPi[theBin]->Fill(impactXY);
	  else if(minPos==1) fHistDCAPosK[theBin]->Fill(impactXY);
	  else if(minPos==2) fHistDCAPosP[theBin]->Fill(impactXY);
	}else{
	  if(minPos==0) fHistDCANegPi[theBin]->Fill(impactXY);
	  else if(minPos==1) fHistDCANegK[theBin]->Fill(impactXY);
	  else if(minPos==2) fHistDCANegP[theBin]->Fill(impactXY);
	}
      } 
    }
    
    //DCA cut on xy and z
    if(!DCAcut(impactXY,impactZ,pt,fMC)) continue;
    fHistNTracks->Fill(10);
    
    if(min<fMinNSigma && yPart<fMaxY){
      if(track->GetSign()>0) fHistPosNSigma[minPos]->Fill(pt);
      else fHistNegNSigma[minPos]->Fill(pt);
      if(fMC){
	if(isph==1){
	  if(track->GetSign()>0) fHistPosNSigmaPrim[minPos]->Fill(pt);
	  else fHistNegNSigmaPrim[minPos]->Fill(pt);
	  if(TMath::Abs(code)==listcode[minPos]){
	    if(track->GetSign()>0) fHistPosNSigmaPrimMC[minPos]->Fill(pt);
	    else fHistNegNSigmaPrimMC[minPos]->Fill(pt);
	  }
	}
      }
    }
    
    if(theBin>=0 && theBin<kNbins){
      if(track->GetSign()>0){
	if(TMath::Abs(y[0]) < fMaxY)fHistPosPi[theBin]->Fill(logdiff[0]);
	if(TMath::Abs(y[1]) < fMaxY)fHistPosK[theBin]->Fill(logdiff[1]);
	if(TMath::Abs(y[2]) < fMaxY)fHistPosP[theBin]->Fill(logdiff[2]);
	if(fMC){
	  if((TMath::Abs(y[0])<fMaxY) && (TMath::Abs(code)==211))  fHistMCPosPi[theBin]->Fill(logdiff[0]);
	  if((TMath::Abs(y[1])<fMaxY) && (TMath::Abs(code)==321))  fHistMCPosK[theBin]->Fill(logdiff[1]);
	  if((TMath::Abs(y[2])<fMaxY) && (TMath::Abs(code)==2212)) fHistMCPosP[theBin]->Fill(logdiff[2]);
	}
      }else{
	if(TMath::Abs(y[0]) < fMaxY)fHistNegPi[theBin]->Fill(logdiff[0]);
	if(TMath::Abs(y[1]) < fMaxY)fHistNegK[theBin]->Fill(logdiff[1]);
	if(TMath::Abs(y[2]) < fMaxY)fHistNegP[theBin]->Fill(logdiff[2]);
	if(fMC){
	  if((TMath::Abs(y[0])<fMaxY) && (TMath::Abs(code)==211))  fHistMCNegPi[theBin]->Fill(logdiff[0]);
	  if((TMath::Abs(y[1])<fMaxY) && (TMath::Abs(code)==321))  fHistMCNegK[theBin]->Fill(logdiff[1]);
	  if((TMath::Abs(y[2])<fMaxY) && (TMath::Abs(code)==2212)) fHistMCNegP[theBin]->Fill(logdiff[2]);
	}
      }
    }							     
  
	  
    //fill propaganda plot with dedx
    fHistDEDX->Fill(track->GetP(),dedx);
    fHistDEDXdouble->Fill(track->GetP()*track->GetSign(),dedx);
	  
    //fill charge distribution histo to check the calibration
    for(Int_t j=0;j<4;j++){
      if(s[j]<5) continue;
      fHistCharge[j]->Fill(s[j]);
    }
  }
		
  // Post output data.
  PostData(1,fOutput);
  Printf("............. end of Exec");
}      

//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectra::Terminate(Option_t *) {
  
  fOutput = dynamic_cast<TList*>(GetOutputData(1));
  if (!fOutput) {
    printf("ERROR: fOutput not available\n");
    return;
  } 
  fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEvents"));
  fHistNTracks = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNTracks"));
  fHistDEDX = dynamic_cast<TH2F*>(fOutput->FindObject("fHistDEDX"));
  fHistDEDXdouble = dynamic_cast<TH2F*>(fOutput->FindObject("fHistDEDXdouble"));

  fHistBeforeEvSel = dynamic_cast<TH1F*>(fOutput->FindObject("fHistBeforeEvSel"));
  fHistAfterEvSel = dynamic_cast<TH1F*>(fOutput->FindObject("fHistAfterEvSel"));

	
  for(Int_t j=0;j<3;j++){
    fHistMCpos[j] = dynamic_cast<TH1F*>(fOutput->FindObject(Form("fHistMCpos%d",j)));
    fHistMCneg[j] = dynamic_cast<TH1F*>(fOutput->FindObject(Form("fHistMCneg%d",j)));
  }

  
  for(Int_t j=0;j<3;j++){
    fHistMCposBefEvSel[j] = dynamic_cast<TH1F*>(fOutput->FindObject(Form("fHistMCposBefEvSel%d",j)));
    fHistMCnegBefEvSel[j] = dynamic_cast<TH1F*>(fOutput->FindObject(Form("fHistMCnegBefEvSel%d",j)));
  }

  for(Int_t i=0; i<4; i++){
    fHistCharge[i] = dynamic_cast<TH1F*>(fOutput->FindObject(Form("fHistChargeLay%d",i)));
  }

  for(Int_t i=0; i<kNbins; i++){
    fHistPosPi[i] = dynamic_cast<TH1F*>(fOutput->FindObject(Form("fHistPosPi%d",i)));
    fHistPosK[i]  = dynamic_cast<TH1F*>(fOutput->FindObject(Form("fHistPosK%d",i)));
    fHistPosP[i]  = dynamic_cast<TH1F*>(fOutput->FindObject(Form("fHistPosP%d",i)));
    fHistNegPi[i] = dynamic_cast<TH1F*>(fOutput->FindObject(Form("fHistNegPi%d",i)));
    fHistNegK[i]  = dynamic_cast<TH1F*>(fOutput->FindObject(Form("fHistNegK%d",i)));
    fHistNegP[i]  = dynamic_cast<TH1F*>(fOutput->FindObject(Form("fHistNegP%d",i)));
    
    fHistDCAPosPi[i] = dynamic_cast<TH1F*>(fOutput->FindObject(Form("fHistDCAPosPi%d",i)));
    fHistDCAPosK[i]  = dynamic_cast<TH1F*>(fOutput->FindObject(Form("fHistDCAPosK%d",i)));
    fHistDCAPosP[i]  = dynamic_cast<TH1F*>(fOutput->FindObject(Form("fHistDCAPosP%d",i)));
    fHistDCANegPi[i] = dynamic_cast<TH1F*>(fOutput->FindObject(Form("fHistDCANegPi%d",i)));
    fHistDCANegK[i]  = dynamic_cast<TH1F*>(fOutput->FindObject(Form("fHistDCANegK%d",i)));
    fHistDCANegP[i]  = dynamic_cast<TH1F*>(fOutput->FindObject(Form("fHistDCANegP%d",i)));
    
    
    if(fMC){	
      fHistMCPosPi[i] = dynamic_cast<TH1F*>(fOutput->FindObject(Form("fHistMCPosPi%d",i)));
      fHistMCPosK[i]  = dynamic_cast<TH1F*>(fOutput->FindObject(Form("fHistMCPosK%d",i)));
      fHistMCPosP[i]  = dynamic_cast<TH1F*>(fOutput->FindObject(Form("fHistMCPosP%d",i)));
      fHistMCNegPi[i] = dynamic_cast<TH1F*>(fOutput->FindObject(Form("fHistMCNegPi%d",i)));
      fHistMCNegK[i]  = dynamic_cast<TH1F*>(fOutput->FindObject(Form("fHistMCNegK%d",i)));
      fHistMCNegP[i]  = dynamic_cast<TH1F*>(fOutput->FindObject(Form("fHistMCNegP%d",i)));
    }
  }

  for(Int_t j=0;j<3;j++){
    fHistPosNSigma[j] = dynamic_cast<TH1F*>(fOutput->FindObject(Form("hHistPosNSigma%d",j)));
    fHistNegNSigma[j] = dynamic_cast<TH1F*>(fOutput->FindObject(Form("hHistNegNSigma%d",j)));
    fHistPosNSigmaPrim[j] = dynamic_cast<TH1F*>(fOutput->FindObject(Form("hHistPosNSigmaPrim%d",j)));
    fHistNegNSigmaPrim[j] = dynamic_cast<TH1F*>(fOutput->FindObject(Form("hHistNegNSigmaPrim%d",j)));
    fHistPosNSigmaPrimMC[j] = dynamic_cast<TH1F*>(fOutput->FindObject(Form("hHistPosNSigmaPrimMC%d",j)));
    fHistNegNSigmaPrimMC[j] = dynamic_cast<TH1F*>(fOutput->FindObject(Form("hHistNegNSigmaPrimMC%d",j)));
  }
  
  fNtupleNSigma = dynamic_cast<TNtuple*>(fOutput->FindObject("fNtupleNSigma"));
  fNtupleMC = dynamic_cast<TNtuple*>(fOutput->FindObject("fNtupleMC"));
  
  Printf("end of Terminate");
  return;
}
