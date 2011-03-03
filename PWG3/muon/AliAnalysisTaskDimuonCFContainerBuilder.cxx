#ifndef ALIANALYSISTASKDIMUONCFCONTAINERBUILDER_CXX
#define ALIANALYSISTASKDIMUONCFCONTAINERBUILDER_CXX

/* $Id$ */

#include "AliAnalysisTaskDimuonCFContainerBuilder.h"
#include "AliStack.h"
#include "TParticle.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TH2D.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliCFManager.h"
#include "AliCFContainer.h"
#include "AliESDMuonTrack.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"

//	Analysis task for the building of a dimuon CF container
//	Also some single-muon variables are stored
//	All the variable ranges and binning are set in the AddTask macro
//	
//	author: L. Bianchi - Universita' & INFN Torino


ClassImp(AliAnalysisTaskDimuonCFContainerBuilder)

//__________________________________________________________________________
AliAnalysisTaskDimuonCFContainerBuilder::AliAnalysisTaskDimuonCFContainerBuilder() :
  fReadAODData(0),
  fReadMCInfo(kTRUE),
  fIsAccProduction(kTRUE),
  fCFManager(0x0),
  fQAHistList(0x0),
  fNevt(0),
  fBeamEnergy(3500.),
  fOutput(0x0),
  fCutOnzVtxSPD(kFALSE),
  fCutOnNContributors(kFALSE),
  fTrigClassMuon(""),
  fTrigClassInteraction(""),
  fDistinguishTrigClass(kFALSE)
{
  
  //Default constructor
  
  Double_t chilims[2]={0.,10000.};
  Double_t ptlims[2]={0.,100.};
  Double_t thetalims[2]={0.,180.};
  Double_t vtxlims[2]={-1000.,1000.};
  SetChi2Limits(chilims);
  SetChi2MatchLimits(chilims);
  SetPtSingMuLimits(ptlims);
  SetThetaSingMuLimits(thetalims);
  SetZprimVertLimits(vtxlims);
  SetTrigClassMuonName();
  SetTrigClassInteracName();
  TString namemuonside[4]={"CMUS1A-","CMUS1B-","CMUS1C-","CMUS1-E-"};
  TString nameinteractionside[4]={"CINT1A-","CINT1B-","CINT1C-","CINT1-E"};
  SetTrigClassMuonSideName(namemuonside);
  SetTrigClassInteracSideName(nameinteractionside);
}
//___________________________________________________________________________
AliAnalysisTaskDimuonCFContainerBuilder::AliAnalysisTaskDimuonCFContainerBuilder(const Char_t* name, Bool_t readaod, Bool_t readMC, Bool_t isaccept, Double_t beamEn) :
  AliAnalysisTaskSE(name),
  fReadAODData(0),
  fReadMCInfo(kTRUE),
  fIsAccProduction(kTRUE),
  fCFManager(0x0),
  fQAHistList(0x0),
  fNevt(0),
  fBeamEnergy(3500.),
  fOutput(0x0),
  fCutOnzVtxSPD(kFALSE),
  fCutOnNContributors(kFALSE),
  fTrigClassMuon(""),
  fTrigClassInteraction(""),
  fDistinguishTrigClass(kFALSE)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskDimuonCFContainerBuilder","Calling Constructor");

  SetReadAODData(readaod);
  SetReadMCinfo(readMC);
  SetIsAccProd(isaccept);
  SetBeamEnergy(beamEn);

  Double_t chilims[2]={0.,10000.};
  Double_t ptlims[2]={0.,100.};
  Double_t thetalims[2]={0.,180.};
  Double_t vtxlims[2]={-1000.,1000.};
  SetChi2Limits(chilims);
  SetChi2MatchLimits(chilims);
  SetPtSingMuLimits(ptlims);
  SetThetaSingMuLimits(thetalims);
  SetZprimVertLimits(vtxlims);
  SetTrigClassMuonName();
  SetTrigClassInteracName();
  TString namemuonside[4]={"CMUS1A-","CMUS1B-","CMUS1C-","CMUS1-E-"};
  TString nameinteractionside[4]={"CINT1A-","CINT1B-","CINT1C-","CINT1-E"};
  SetTrigClassMuonSideName(namemuonside);
  SetTrigClassInteracSideName(nameinteractionside);
  
  DefineOutput(1,TList::Class());
  DefineOutput(2,AliCFContainer::Class());

}

//___________________________________________________________________________
AliAnalysisTaskDimuonCFContainerBuilder& AliAnalysisTaskDimuonCFContainerBuilder::operator=(const AliAnalysisTaskDimuonCFContainerBuilder& c) 
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c) ;
    fReadAODData = c.fReadAODData ;
    fCFManager  = c.fCFManager;
    fQAHistList = c.fQAHistList ;
    fNevt = c.fNevt ;
  }
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskDimuonCFContainerBuilder::AliAnalysisTaskDimuonCFContainerBuilder(const AliAnalysisTaskDimuonCFContainerBuilder& c) :
  AliAnalysisTaskSE(c),
  fReadAODData(c.fReadAODData),
  fReadMCInfo(c.fReadMCInfo),
  fIsAccProduction(c.fIsAccProduction),
  fCFManager(c.fCFManager),
  fQAHistList(c.fQAHistList),
  fNevt(c.fNevt),
  fBeamEnergy(c.fBeamEnergy),
  fOutput(c.fOutput),
  fCutOnzVtxSPD(c.fCutOnzVtxSPD),
  fCutOnNContributors(c.fCutOnNContributors),
  fTrigClassMuon(c.fTrigClassMuon),
  fTrigClassInteraction(c.fTrigClassInteraction),
  fDistinguishTrigClass(c.fDistinguishTrigClass)
{
  
  // Copy Constructor
  
}

//___________________________________________________________________________
AliAnalysisTaskDimuonCFContainerBuilder::~AliAnalysisTaskDimuonCFContainerBuilder() {
  //
  //destructor
  //
  Info("~AliAnalysisTaskDimuonCFContainerBuilder","Calling Destructor");
  if (fCFManager)           delete fCFManager ;
  if (fQAHistList) {fQAHistList->Clear(); delete fQAHistList;}
}

//___________________________________________________________________________
void AliAnalysisTaskDimuonCFContainerBuilder::UserCreateOutputObjects(){
 //	UserCreateOutputObjects
 fOutput = new TList();
 fOutput->SetOwner(); 
 
 TH1D *hnevts	= new TH1D("hnevts","hnevts",1,0,1);					// Stat check histos
 TH1D *jpsiMult = new TH1D("jpsiMult","jpsiMult",20,0,20);	

 TH1D *ptdimuREC    = new TH1D("ptdimuREC","ptdimuREC",200,0,20);			// Dimu check histos
 TH1D *ydimuREC     = new TH1D("ydimuREC","ydimuREC",200,-10.,10.);	
 TH1D *ydimuFail     = new TH1D("ydimuFail","ydimuFail",10,660,670);
 TH1D *costHEdimuREC= new TH1D("costHEdimuREC","costHEdimuREC",200,-1.,1.);	
 TH1D *costCSdimuREC= new TH1D("costCSdimuREC","costCSdimuREC",200,-1.,1.);	
 TH1D *costHEdimuFail= new TH1D("costHEdimuFail","costHEdimuFail",10,660.,670.);	
 TH1D *costCSdimuFail= new TH1D("costCSdimuFail","costCSdimuFail",10,660.,670.);	
 TH1D *phiHEdimuREC = new TH1D("phiHEdimuREC","phiHEdimuREC",100,0.,TMath::Pi());	
 TH1D *phiCSdimuREC = new TH1D("phiCSdimuREC","phiCSdimuREC",100,0.,TMath::Pi());	
 TH1D *phiHEdimuFail = new TH1D("phiHEdimuFail","phiHEdimuFail",10,660.,670.);	
 TH1D *phiCSdimuFail = new TH1D("phiCSdimuFail","phiCSdimuFail",10,660.,670.);	
 TH1D *imassTot     = new TH1D("imassTot","imassTot",100,0,8);	
 TH1D *trigCond     = new TH1D("trigCond","trigCond",40,0,4);	

 TH1D *emuonREC	= new TH1D("emuonREC","emuonREC",200,0.,20.);				// Mu check histos
 TH1D *ptmuonREC= new TH1D("ptmuonREC","ptmuonREC",200,0.,20.);
 TH1D *ymuonREC	= new TH1D("ymuonREC","ymuonREC",200,-10.,10.);	
 TH1D *hdca	= new TH1D("hdca","hdca",20,0,200);
 TH2D *hdcay	= new TH2D("hdcay","hdcay",200,0,200,20,-5,0);

 TH1D *zvSPD 	= new TH1D("zvSPD","zvSPD",100,-50,50.);				// Event check histos
 TH1D *zvSPDcut	= new TH1D("zvSPDcut","zvSPDcut",100,-50,50.);
 TH1D *nContrib	= new TH1D("nContrib","nContrib",100,0,100);


 fOutput->Add(hnevts);
 fOutput->Add(jpsiMult); 

 fOutput->Add(ptdimuREC); 
 fOutput->Add(ydimuREC); 
 fOutput->Add(ydimuFail); 
 fOutput->Add(costHEdimuREC);
 fOutput->Add(costCSdimuREC);
 fOutput->Add(costHEdimuFail);
 fOutput->Add(costCSdimuFail);
 fOutput->Add(phiHEdimuREC);
 fOutput->Add(phiCSdimuREC);
 fOutput->Add(phiHEdimuFail);
 fOutput->Add(phiCSdimuFail);
 fOutput->Add(imassTot); 
 fOutput->Add(trigCond); 

 fOutput->Add(emuonREC); 
 fOutput->Add(ptmuonREC);
 fOutput->Add(ymuonREC);
 fOutput->Add(hdca);
 fOutput->Add(hdcay);
 
 fOutput->Add(zvSPD); 
 fOutput->Add(zvSPDcut); 
 fOutput->Add(nContrib); 

	
 fOutput->ls();
 

} 


//_________________________________________________
void AliAnalysisTaskDimuonCFContainerBuilder::UserExec(Option_t *)
{
  
  
  //Info("UserExec","");
  fNevt++;
  ((TH1D*)(fOutput->FindObject("hnevts")))->Fill(0.5);

  Double_t containerInput[15]={666,666,666,666,666,666,666,666,666,666,666,666,666,666,666};   //y, pT, costHE, phiHE, costCS, phiCS, mass, TrigCond
  Int_t cutAccept =1;
  Int_t numbJpsis =0;
   
  if (!fReadAODData){     // ESD-based ANALYSIS

    if (fReadMCInfo){

      if (!fMCEvent) {
	Error("UserExec","NO MC EVENT FOUND!");
	//return;
      }
 
      // MC part  ---------------------------------------------------------------------------------------------------

      //fCFManager->SetEventInfo(fMCEvent);					CHANGES IN NEW VERSIONS - MANUAL CUT ON PDG
      AliStack* stack = fMCEvent->Stack();
 
      // loop on the MC event
      for (Int_t ipart=0; ipart<fMCEvent->GetNumberOfTracks(); ipart++) { 

	AliMCParticle *mcPart  = (AliMCParticle*) fMCEvent->GetTrack(ipart);
  
	TParticle *part = mcPart->Particle(); 

	// Mother kinematics
	Double_t e = part->Energy();
	Double_t pz = part->Pz();		     
	Double_t rapmc = Rap(e,pz);

	// Selection of the resonance
	//if (!fCFManager->CheckParticleCuts(0,mcPart)) continue;
	if (part->GetPdgCode()!=443) continue;					// MANUAL CUT ON PDG
	numbJpsis++;
	
	// Decays kinematics
	Int_t p0 = part->GetDaughter(0);
	TParticle *part0 = stack->Particle(p0); 
	Int_t pdg0 = part0->GetPdgCode();
 
	Int_t p1 = part->GetDaughter(1);
	TParticle *part1 = stack->Particle(p1);
	Int_t pdg1 = part1->GetPdgCode();
  
	Double_t e0 = part0->Energy();
	Double_t pz0 = part0->Pz();
	Double_t py0 = part0->Py();
	Double_t px0 = part0->Px();
	Double_t charge0 = part0->GetPDG()->Charge()/3;
	Double_t theta0 = (180./TMath::Pi())*TMath::ATan2(TMath::Sqrt(px0*px0+py0*py0),pz0);
	Double_t pt0 = TMath::Sqrt(px0*px0+py0*py0);
	Double_t mom0 = part0->P();

	Double_t e1 = part1->Energy();
	Double_t pz1 = part1->Pz();
	Double_t py1 = part1->Py();
	Double_t px1 = part1->Px();
	//Double_t charge1 = part1->GetPDG()->Charge()/3;
	Double_t theta1 = (180./TMath::Pi())*TMath::ATan2(TMath::Sqrt(px1*px1+py1*py1),pz1);
	Double_t pt1 = TMath::Sqrt(px1*px1+py1*py1);
	Double_t mom1 = part1->P();

     
	if(pdg0==13 || pdg1==13) { 
 
	  Double_t ptmc = TMath::Sqrt((px0+px1)*(px0+px1)+(py0+py1)*(py0+py1));
	  Double_t imassmc = Imass(e0,px0,py0,pz0,e1,px1,py1,pz1);
      
	  Double_t costCS = CostCS(px0,py0,pz0,e0,charge0,px1,py1,pz1,e1);
	  Double_t costHE = CostHE(px0,py0,pz0,e0,charge0,px1,py1,pz1,e1);
	  Double_t phiCS  = PhiCS(px0,py0,pz0,e0,charge0,px1,py1,pz1,e1);
	  Double_t phiHE  = PhiHE(px0,py0,pz0,e0,charge0,px1,py1,pz1,e1);
 
	  containerInput[0]  = rapmc ;   
	  containerInput[1]  = ptmc;
	  containerInput[2]  = costHE;    
	  containerInput[3]  = TMath::Abs(phiHE);     
	  containerInput[4]  = costCS;	
	  containerInput[5]  = TMath::Abs(phiCS);
	  containerInput[6]  = imassmc;      
	  containerInput[7]  = 1.;		//for generated no trigger condition
	  if (pt0<pt1) {
	    containerInput[8]=pt0; 
	    containerInput[9]=pt1;
	  } else {
	    containerInput[8]=pt1; 
	    containerInput[9]=pt0;
	  }
	  if (theta0<theta1) {
	    containerInput[10]=theta0; 
	    containerInput[11]=theta1;
	  } else {
	    containerInput[10]=theta1; 
	    containerInput[11]=theta0;
	  }
	  if (mom0<mom1) {
	    containerInput[12]=mom0; 
	    containerInput[13]=mom1;
	  } else {
	    containerInput[12]=mom1; 
	    containerInput[13]=mom0;
	  }
	  containerInput[14]=1.;
 
	  // fill the container at the first step
	  fCFManager->GetParticleContainer()->Fill(containerInput,0);
	  	
	  if(fIsAccProduction){
	    // acceptance cuts on single mu
	    if(theta0<fThetaSingMuCut[0]||theta0>fThetaSingMuCut[1]||theta1<fThetaSingMuCut[0]||theta1>fThetaSingMuCut[1]) cutAccept=0;
	    if(pt0<fPtSingMuCut[0] || pt0>fPtSingMuCut[1] || pt1<fPtSingMuCut[0] || pt1>fPtSingMuCut[1]) cutAccept=0;
	    if (cutAccept ==1) fCFManager->GetParticleContainer()->Fill(containerInput,2);
	  }
	}
      } 
    } //fReadMCInfo


      // ESD part  ---------------------------------------------------------------------------------------------------
  
      AliESDEvent *fESD; 
      AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
      if ( ! esdH ) {
        AliError("Cannot get input event handler");
        return;
      }  
      fESD = esdH->GetEvent();
      Int_t mult1 = fESD->GetNumberOfMuonTracks() ;

      Int_t trigfired=-1;
      Int_t trigside=-1;
      if(fDistinguishTrigClass){
	//aod->GetHeader()->SetFiredTriggerClasses("CMU");
	TString trigclass = fESD->GetFiredTriggerClasses();
	printf("TrigClass: %s\n",trigclass.Data());
        if(trigclass.Contains(fTrigClassMuon)) trigfired = 1;
        else if (trigclass.Contains(fTrigClassInteraction)) trigfired = 0;
	printf("Muon container: %d\n",trigfired);
	for(Int_t i=0;i<4;i++){
	  if(trigfired==1 && trigclass.Contains(fTrigClassMuonSide[i])) {trigside=i;printf("Muon Side: %d\n",trigside);}
	  if(trigfired==0 && trigclass.Contains(fTrigClassInteractionSide[i])) {trigside=i;printf("Interaction Side: %d\n",trigside);}
	}
      }

      ((TH1D*)(fOutput->FindObject("zvSPD")))->Fill(fESD->GetPrimaryVertexSPD()->GetZ());
      ((TH1D*)(fOutput->FindObject("nContrib")))->Fill(fESD->GetPrimaryVertexSPD()->GetNContributors());
      if (!fCutOnzVtxSPD || (fCutOnzVtxSPD && (fESD->GetPrimaryVertexSPD()->GetZ()>fzPrimVertexSPD[0] && fESD->GetPrimaryVertexSPD()->GetZ()<fzPrimVertexSPD[1]))){
      ((TH1D*)(fOutput->FindObject("zvSPDcut")))->Fill(fESD->GetPrimaryVertexSPD()->GetZ());
      if (!fCutOnNContributors || (fCutOnNContributors && (fESD->GetPrimaryVertexSPD()->GetNContributors()>fNContributors[0] && fESD->GetPrimaryVertexSPD()->GetNContributors()<fNContributors[1]))){
      for (Int_t j = 0; j<mult1; j++) { 

	AliESDMuonTrack* mu1 = new AliESDMuonTrack(*(fESD->GetMuonTrack(j)));
	if (!mu1->ContainTrackerData()) continue;
	if (mu1->GetChi2()<fChi2Track[0] || mu1->GetChi2()>fChi2Track[1]) continue;
	if (mu1->GetChi2MatchTrigger()<fChi2MatchTrig[0] || mu1->GetChi2MatchTrigger()>fChi2MatchTrig[1]) continue;
	Double_t zr1 = mu1->Charge();
	Double_t pxr1  = mu1->Px();
	Double_t pyr1  = mu1->Py();
	Double_t pzr1  = mu1->Pz();
        Double_t ptmu1 = TMath::Sqrt(pxr1*pxr1+pyr1*pyr1);
 	Double_t er1 = mu1->E();
	Double_t mom1 = mu1->P();
	Double_t theta1 = (180./TMath::Pi())*mu1->Theta();
	Double_t rapiditymu1 = Rap(er1,pzr1); 
	((TH1D*)(fOutput->FindObject("emuonREC")))->Fill(er1);
	((TH1D*)(fOutput->FindObject("ptmuonREC")))->Fill(ptmu1);
	((TH1D*)(fOutput->FindObject("ymuonREC")))->Fill(rapiditymu1);

	if(zr1<0){

	    for (Int_t jj = 0; jj<mult1; jj++) {

	      AliESDMuonTrack* mu2 = new AliESDMuonTrack(*(fESD->GetMuonTrack(jj)));
	      if (!mu2->ContainTrackerData()) continue;
	      if (mu2->GetChi2()<fChi2Track[0] || mu2->GetChi2()>fChi2Track[1]) continue;
	      if (mu2->GetChi2MatchTrigger()<fChi2MatchTrig[0] || mu2->GetChi2MatchTrigger()>fChi2MatchTrig[1]) continue;
	      Double_t zr2 = mu2->Charge();

	      if(zr2>0){
		Double_t trigCondition=0;
		if (mu1->GetMatchTrigger()>=mu2->GetMatchTrigger()) trigCondition = mu1->GetMatchTrigger()+0.1*mu2->GetMatchTrigger();
		else trigCondition = mu2->GetMatchTrigger()+0.1*mu1->GetMatchTrigger();
		((TH1D*)(fOutput->FindObject("trigCond")))->Fill(trigCondition);
		Double_t pxr2 = mu2->Px();
		Double_t pyr2 = mu2->Py();
		Double_t pzr2 = mu2->Pz();
                Double_t ptmu2 = TMath::Sqrt(pxr2*pxr2+pyr2*pyr2);
		Double_t er2 = mu2->E();
		Double_t mom2 = mu2->P();
		Double_t theta2 = (180./TMath::Pi())*mu2->Theta();

		Double_t ptrec = TMath::Sqrt((pxr1+pxr2)*(pxr1+pxr2)+(pyr1+pyr2)*(pyr1+pyr2));
		((TH1D*)(fOutput->FindObject("ptdimuREC")))->Fill(ptrec);
		Double_t raprec= Rap((er1+er2),(pzr1+pzr2));
		((TH1D*)(fOutput->FindObject("ydimuREC")))->Fill(raprec);
		Double_t imassrec = Imass(er1,pxr1,pyr1,pzr1,er2,pxr2,pyr2,pzr2);
		((TH1D*)(fOutput->FindObject("imassTot")))->Fill(imassrec);

		if (imassrec>12.) continue;
		    
		Double_t costCSrec = CostCS(pxr1,pyr1,pzr1,er1,zr1,pxr2,pyr2,pzr2,er2);
		((TH1D*)(fOutput->FindObject("costCSdimuREC")))->Fill(costCSrec);
		if((Int_t)costCSrec==666) ((TH1D*)(fOutput->FindObject("costCSdimuFail")))->Fill(costCSrec);
		Double_t costHErec = CostHE(pxr1,pyr1,pzr1,er1,zr1,pxr2,pyr2,pzr2,er2);
		((TH1D*)(fOutput->FindObject("costHEdimuREC")))->Fill(costHErec);
		if((Int_t)costHErec==666) ((TH1D*)(fOutput->FindObject("costHEdimuFail")))->Fill(costHErec);
		Double_t phiCSrec  = PhiCS(pxr1,pyr1,pzr1,er1,zr1,pxr2,pyr2,pzr2,er2);
		((TH1D*)(fOutput->FindObject("phiCSdimuREC")))->Fill(phiCSrec);
		if((Int_t)phiCSrec==666) ((TH1D*)(fOutput->FindObject("phiCSdimuFail")))->Fill(phiCSrec);
		Double_t phiHErec  = PhiHE(pxr1,pyr1,pzr1,er1,zr1,pxr2,pyr2,pzr2,er2);
		((TH1D*)(fOutput->FindObject("phiHEdimuREC")))->Fill(phiHErec);
		if((Int_t)phiHErec==666) ((TH1D*)(fOutput->FindObject("phiHEdimuFail")))->Fill(phiHErec);

		containerInput[0] = raprec ;   
		containerInput[1] = ptrec;
		containerInput[2] = costHErec;	
		containerInput[3] = TMath::Abs(phiHErec);	
		containerInput[4] = costCSrec;	
		containerInput[5] = TMath::Abs(phiCSrec);
		containerInput[6] = imassrec;  
		containerInput[7] = trigCondition+0.05;
		if (ptmu1<ptmu2) {
		  containerInput[8]=ptmu1; 
		  containerInput[9]=ptmu2;
		} else {
		  containerInput[8]=ptmu2; 
		  containerInput[9]=ptmu1;
		}
		if (theta1<theta2) {
		  containerInput[10]=theta1; 
		  containerInput[11]=theta2;
		} else {
		  containerInput[10]=theta2; 
		  containerInput[11]=theta1;
		}
		if (mom1<mom2) {
		  containerInput[12]=mom1; 
		  containerInput[13]=mom2;
		} else {
		  containerInput[12]=mom2; 
		  containerInput[13]=mom1;
		}
		if (fDistinguishTrigClass && trigside==0) containerInput[14]=0.;
		else if (fDistinguishTrigClass && trigside==1) containerInput[14]=1.;
		else if (fDistinguishTrigClass && trigside==2) containerInput[14]=2.;
		else if (fDistinguishTrigClass && trigside==3) containerInput[14]=3.;
		else containerInput[14]=0.;
		    
		if(fDistinguishTrigClass){    
		  if (trigfired==1) fCFManager->GetParticleContainer()->Fill(containerInput,5);
		  else if (trigfired==0) fCFManager->GetParticleContainer()->Fill(containerInput,1);
		} else {
		  fCFManager->GetParticleContainer()->Fill(containerInput,1);
		  if(fIsAccProduction){
		    if(cutAccept==1) fCFManager->GetParticleContainer()->Fill(containerInput,3);
		  }
		}
	      }  // mu+ Selection

	    }      // second mu Loop
	 }          // mu- Selection
     }  
	}
      }
  } else {     // AOD-based ANALYSIS

    AliAODEvent *aod;
    AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    aod = aodH->GetEvent();
    Int_t ntracks=aod->GetNumberOfTracks(); 

    // MC part  ---------------------------------------------------------------------------------------------------

    if (fReadMCInfo){
      TClonesArray *mcarray = dynamic_cast<TClonesArray*> (aod->FindListObject(AliAODMCParticle::StdBranchName()));  //array of MC particles in this event
      for(int ii=0;ii<mcarray->GetEntries();ii++){
	AliAODMCParticle *mctrack = (AliAODMCParticle*) mcarray->At(ii);
	if(mctrack->GetPdgCode()!=13) continue;
	Int_t numbMother = mctrack->GetMother();
	if (numbMother==-1) continue;
	AliAODMCParticle *mother = (AliAODMCParticle*) mcarray->At(numbMother);
	if (mother->GetPdgCode()!=443) continue;
	numbJpsis++;
	Int_t daught0 = mother->GetDaughter(0);
	Int_t daught1 = mother->GetDaughter(1);
	AliAODMCParticle *mcDaughter0 = (AliAODMCParticle*) mcarray->At(daught0);
	Double_t pxmc0 = mcDaughter0->Px();
	Double_t pymc0 = mcDaughter0->Py();
	Double_t pzmc0 = mcDaughter0->Pz();
	Double_t ptmc0 = mcDaughter0->Pt();
	Double_t mommc0 = mcDaughter0->P();
	Double_t emc0  = mcDaughter0->E();
	Double_t thetamc0 = (180./TMath::Pi())*mcDaughter0->Theta();
	Int_t charge0 = (Int_t) mcDaughter0->Charge()/3;
	AliAODMCParticle *mcDaughter1 = (AliAODMCParticle*) mcarray->At(daught1);
	Double_t pxmc1 = mcDaughter1->Px();
	Double_t pymc1 = mcDaughter1->Py();
	Double_t pzmc1 = mcDaughter1->Pz();
	Double_t ptmc1 = mcDaughter1->Pt();
	Double_t mommc1 = mcDaughter1->P();
	Double_t emc1  = mcDaughter1->E();
	Double_t thetamc1 = (180./TMath::Pi())*mcDaughter1->Theta();
	Int_t charge1 = (Int_t) mcDaughter1->Charge()/3;
      
	if (charge0==charge1 || TMath::Abs(mcDaughter0->GetPdgCode())!=13 || TMath::Abs(mcDaughter1->GetPdgCode())!=13) continue;
	Double_t rapJpsi = Rap(mother->E(),mother->Pz());
	Double_t ptJpsi = mother->Pt();
	Double_t costHE = CostHE(pxmc0,pymc0,pzmc0,emc0,(Double_t)charge0,pxmc1,pymc1,pzmc1,emc1);
	Double_t phiHE  = PhiHE(pxmc0,pymc0,pzmc0,emc0,(Double_t)charge0,pxmc1,pymc1,pzmc1,emc1);
	Double_t costCS = CostCS(pxmc0,pymc0,pzmc0,emc0,(Double_t)charge0,pxmc1,pymc1,pzmc1,emc1);
	Double_t phiCS  = PhiCS(pxmc0,pymc0,pzmc0,emc0,(Double_t)charge0,pxmc1,pymc1,pzmc1,emc1);
	Double_t massJpsi = mother->M();
	containerInput[0] = rapJpsi;
	containerInput[1] = ptJpsi;
	containerInput[2] = costHE;
	containerInput[3] = TMath::Abs(phiHE);
	containerInput[4] = costCS;
	containerInput[5] = TMath::Abs(phiCS);
	containerInput[6] = massJpsi;
	containerInput[7] = 1.;		//for generated no trigger condition
	if (ptmc0<ptmc1) {
	  containerInput[8]=ptmc0; 
	  containerInput[9]=ptmc1;
	} else {
	  containerInput[8]=ptmc1; 
	  containerInput[9]=ptmc0;
	}
	if (thetamc0<thetamc1) {
	  containerInput[10]=thetamc0; 
	  containerInput[11]=thetamc1;
	} else {
	  containerInput[10]=thetamc1; 
	  containerInput[11]=thetamc0;
	}
	if (mommc0<mommc1) {
	  containerInput[12]=mommc0; 
	  containerInput[13]=mommc1;
	} else {
	  containerInput[12]=mommc1; 
	  containerInput[13]=mommc0;
	}
	containerInput[14]=1.;
      
	if((Int_t)rapJpsi!=666 && (Int_t)costHE!=666 && (Int_t)phiHE!=666 && (Int_t)costCS!=666 && (Int_t)phiCS!=666){
	  fCFManager->GetParticleContainer()->Fill(containerInput,0);
	  if(fIsAccProduction){
	    // acceptance cuts on single mu
	    if(thetamc0<fThetaSingMuCut[0]||thetamc0>fThetaSingMuCut[1]||thetamc1<fThetaSingMuCut[0]||thetamc1>fThetaSingMuCut[1]) cutAccept=0;
	    if(ptmc0<fPtSingMuCut[0] || ptmc0>fPtSingMuCut[1] || ptmc1<fPtSingMuCut[0] || ptmc1>fPtSingMuCut[1]) cutAccept=0;
	    if (cutAccept ==1) fCFManager->GetParticleContainer()->Fill(containerInput,2);
	  }
	  
	}

	}
      }
      
      
      // AOD part  ---------------------------------------------------------------------------------------------------
      
      Int_t trigfired=-1;
      Int_t trigside=-1;
      if(fDistinguishTrigClass){
	TString trigclass = aod->GetFiredTriggerClasses();
	printf("TrigClass: %s\n",trigclass.Data());
        if(trigclass.Contains(fTrigClassMuon)) trigfired = 1;
        else if (trigclass.Contains(fTrigClassInteraction)) trigfired = 0;
	printf("Muon container: %d\n",trigfired);
	for(Int_t i=0;i<4;i++){
	  if(trigfired==1 && trigclass.Contains(fTrigClassMuonSide[i])) {trigside=i;printf("Muon Side: %d\n",trigside);}
	  if(trigfired==0 && trigclass.Contains(fTrigClassInteractionSide[i])) {trigside=i;printf("Interaction Side: %d\n",trigside);}
	}
      }
      
      ((TH1D*)(fOutput->FindObject("zvSPD")))->Fill(aod->GetPrimaryVertex()->GetZ());
      Int_t ncontr = aod->GetPrimaryVertex()->GetNContributors();
      ((TH1D*)(fOutput->FindObject("nContrib")))->Fill(ncontr);
      if (!fCutOnzVtxSPD || (fCutOnzVtxSPD && (aod->GetPrimaryVertex()->GetZ()>fzPrimVertexSPD[0] && aod->GetPrimaryVertex()->GetZ()<fzPrimVertexSPD[1]))){	//NOT REALLY SPD VERTEX
      ((TH1D*)(fOutput->FindObject("zvSPDcut")))->Fill(aod->GetPrimaryVertex()->GetZ());
      if (!fCutOnNContributors || (fCutOnNContributors && (aod->GetPrimaryVertex()->GetNContributors()>fNContributors[0] && aod->GetPrimaryVertex()->GetNContributors()<fNContributors[1]))){
      for (Int_t j = 0; j<ntracks; j++) {
	AliAODTrack *mu1 = aod->GetTrack(j);
	if(!mu1->IsMuonTrack()) continue;
	if (mu1->Chi2perNDF()<fChi2Track[0] || mu1->Chi2perNDF()>fChi2Track[1]) continue;
	if (mu1->GetChi2MatchTrigger()<fChi2MatchTrig[0] || mu1->GetChi2MatchTrigger()>fChi2MatchTrig[1]) continue;
	Double_t chargemu1 = mu1->Charge();
	Double_t pxmu1 = mu1->Px();
	Double_t pymu1 = mu1->Py();
	Double_t pzmu1 = mu1->Pz();
	Double_t emu1 = mu1->E();
	Double_t pmu1 = mu1->P();
	Double_t ptmu1 = mu1->Pt();
	Double_t rapiditymu1 = Rap(emu1,pzmu1); 
	Double_t thetamu1 = (180./TMath::Pi())*mu1->Theta();
	Double_t rdcamu1 = mu1->DCA();
	((TH1D*)(fOutput->FindObject("emuonREC")))->Fill(emu1);
	((TH1D*)(fOutput->FindObject("ptmuonREC")))->Fill(ptmu1);
	((TH1D*)(fOutput->FindObject("ymuonREC")))->Fill(rapiditymu1);
	if(chargemu1<0){
	  for (Int_t jj = 0; jj<ntracks; jj++) { 
	    AliAODTrack *mu2 = aod->GetTrack(jj);
	    if(!mu2->IsMuonTrack()) continue;
	    if (mu2->Chi2perNDF()<fChi2Track[0] || mu2->Chi2perNDF()>fChi2Track[1]) continue;
	    if (mu2->GetChi2MatchTrigger()<fChi2MatchTrig[0] || mu2->GetChi2MatchTrigger()>fChi2MatchTrig[1]) continue;
	    Double_t chargemu2 = mu2->Charge();
	    Double_t pxmu2 = mu2->Px();
	    Double_t pymu2 = mu2->Py();
	    Double_t pzmu2 = mu2->Pz();
	    Double_t emu2 = mu2->E();
	    Double_t pmu2 = mu2->P();
	    Double_t ptmu2 = mu2->Pt();
	    //Double_t rapiditymu2 = Rap(emu2,pzmu2); 
	    Double_t thetamu2 = (180./TMath::Pi())*mu2->Theta();
	    Double_t rdcamu2 = mu2->DCA();
	    if(chargemu2>0){
	    if (mu1->GetMatchTrigger()>0) printf("Mu1: charge: %f, match: %d\n",chargemu1,1);
	    else  printf("Mu1: charge: %f, match: %d\n",chargemu1,0);
	    if (mu2->GetMatchTrigger()>0) printf("Mu2: charge: %f, match: %d\n",chargemu2,1);
	    else printf("Mu2: charge: %f, match: %d\n",chargemu2,0);
	      ((TH1D*)(fOutput->FindObject("hdca")))->Fill(rdcamu2);
	      ((TH1D*)(fOutput->FindObject("hdca")))->Fill(rdcamu1);
	      Double_t trigCondition=0;
	      if (mu1->GetMatchTrigger()>=mu2->GetMatchTrigger()) trigCondition = mu1->GetMatchTrigger()+0.1*mu2->GetMatchTrigger();
		else trigCondition = mu2->GetMatchTrigger()+0.1*mu1->GetMatchTrigger();	    
	      containerInput[0] = (Double_t) Rap((emu1+emu2),(pzmu1+pzmu2));  
	      ((TH1D*)(fOutput->FindObject("ydimuREC")))->Fill(containerInput[0]);
	      ((TH2D*)(fOutput->FindObject("hdcay")))->Fill(TMath::Max(rdcamu1,rdcamu2),containerInput[0]);
	      //((TH2D*)(fOutput->FindObject("hdcay")))->Fill(Rdcamu2,containerInput[0]);
	      containerInput[1] = (Double_t) TMath::Sqrt((pxmu1+pxmu2)*(pxmu1+pxmu2)+(pymu1+pymu2)*(pymu1+pymu2));
	      ((TH1D*)(fOutput->FindObject("ptdimuREC")))->Fill(containerInput[1]);
	      containerInput[2] = CostHE(pxmu1,pymu1,pzmu1,emu1,chargemu1,pxmu2,pymu2,pzmu2,emu2);
	      ((TH1D*)(fOutput->FindObject("costHEdimuREC")))->Fill(containerInput[2]);
	      if(containerInput[2]==666) ((TH1D*)(fOutput->FindObject("costHEdimuFail")))->Fill(containerInput[2]);
	      containerInput[3] = TMath::Abs(PhiHE(pxmu1,pymu1,pzmu1,emu1,chargemu1,pxmu2,pymu2,pzmu2,emu2));
	      ((TH1D*)(fOutput->FindObject("phiHEdimuREC")))->Fill(containerInput[3]);
	      if(containerInput[3]==666) ((TH1D*)(fOutput->FindObject("phiHEdimuFail")))->Fill(containerInput[3]);
	      containerInput[4] = CostCS(pxmu1,pymu1,pzmu1,emu1,chargemu1,pxmu2,pymu2,pzmu2,emu2);
	      ((TH1D*)(fOutput->FindObject("costCSdimuREC")))->Fill(containerInput[4]);
	      if(containerInput[4]==666) ((TH1D*)(fOutput->FindObject("costCSdimuFail")))->Fill(containerInput[4]);
	      containerInput[5] = TMath::Abs(PhiCS(pxmu1,pymu1,pzmu1,emu1,chargemu1,pxmu2,pymu2,pzmu2,emu2));
	      ((TH1D*)(fOutput->FindObject("phiCSdimuREC")))->Fill(containerInput[5]);
	      if(containerInput[5]==666) ((TH1D*)(fOutput->FindObject("phiCSdimuFail")))->Fill(containerInput[5]);
	      containerInput[6] = Imass(emu1,pxmu1,pymu1,pzmu1,emu2,pxmu2,pymu2,pzmu2);
	      ((TH1D*)(fOutput->FindObject("imassTot")))->Fill(containerInput[6]);
	      containerInput[7] = trigCondition+0.05;
	      ((TH1D*)(fOutput->FindObject("trigCond")))->Fill(containerInput[7]);
	      if (ptmu1<ptmu2) {
		containerInput[8]=ptmu1; 
		containerInput[9]=ptmu2;
	      } else {
		containerInput[8]=ptmu2; 
		containerInput[9]=ptmu1;
	      }
	      if (thetamu1<thetamu2) {
		containerInput[10]=thetamu1; 
		containerInput[11]=thetamu2;
	      } else {
		containerInput[10]=thetamu2; 
		containerInput[11]=thetamu1;
	      }
	      if (pmu1<pmu2) {
		containerInput[12]=pmu1; 
		containerInput[13]=pmu2;
	      } else {
		containerInput[12]=pmu2; 
		containerInput[13]=pmu1;
	      }
	      if (fDistinguishTrigClass && trigside==0) containerInput[14]=0.;
	      else if (fDistinguishTrigClass && trigside==1) containerInput[14]=1.;
	      else if (fDistinguishTrigClass && trigside==2) containerInput[14]=2.;
	      else if (fDistinguishTrigClass && trigside==3) containerInput[14]=3.;
	      else containerInput[14]=0.;
 	      
	      if(containerInput[2]!=666 && containerInput[3]!=666 && containerInput[4]!=666 && containerInput[5]!=666){
		if(fDistinguishTrigClass){
		  if (trigfired==1) fCFManager->GetParticleContainer()->Fill(containerInput,5);
		  else if (trigfired==0) fCFManager->GetParticleContainer()->Fill(containerInput,1);
		} else {
		  fCFManager->GetParticleContainer()->Fill(containerInput,1);
		  if(fIsAccProduction){
		    if (cutAccept ==1) fCFManager->GetParticleContainer()->Fill(containerInput,3);
		  }
		}
	      } else if (containerInput[0]==666) ((TH1D*)(fOutput->FindObject("ydimuFail")))->Fill(containerInput[0]);
	    
	  }

	}

      }

    } 
    }
    }
   

  }


//  ----------
  if (fReadMCInfo) ((TH1D*)(fOutput->FindObject("jpsiMult")))->Fill(numbJpsis);

  PostData(1,fOutput) ;
  PostData(2,fCFManager->GetParticleContainer()) ;
  
  
}

//________________________________________________________________________
Double_t AliAnalysisTaskDimuonCFContainerBuilder::Imass(Double_t e1, Double_t px1, Double_t py1, Double_t pz1,
				   Double_t e2, Double_t px2, Double_t py2, Double_t pz2) const
{
// invariant mass calculation
    Double_t imassrec = TMath::Sqrt((e1+e2)*(e1+e2)-((px1+px2)*(px1+px2)+
                                    (py1+py2)*(py1+py2)+(pz1+pz2)*(pz1+pz2)));
    return imassrec;
}

//________________________________________________________________________
Double_t AliAnalysisTaskDimuonCFContainerBuilder::Rap(Double_t e, Double_t pz) const
{
// calculate rapidity
    Double_t rap;
    if(e>TMath::Abs(pz)){                         // in order to avoid problems with AODs
	rap = 0.5*TMath::Log((e+pz)/(e-pz));
	return rap;
    }
    else{
	rap = 666.;
	return rap;
    }
}

//________________________________________________________________________
Double_t AliAnalysisTaskDimuonCFContainerBuilder::CostCS(Double_t px1, Double_t py1, Double_t pz1, Double_t e1,
Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2)
{
// Cosine of the theta decay angle (mu+) in the Collins-Soper frame

  TLorentzVector pMu1CM, pMu2CM, pProjCM, pTargCM, pDimuCM; // In the CM. frame
  TLorentzVector pMu1Dimu, pMu2Dimu, pProjDimu, pTargDimu; // In the dimuon rest frame
  TVector3 beta,zaxisCS;
  Double_t mp=0.93827231;
  //
  // --- Fill the Lorentz vector for projectile and target in the CM frame
  //
  pProjCM.SetPxPyPzE(0.,0.,-fBeamEnergy,TMath::Sqrt(fBeamEnergy*fBeamEnergy+mp*mp)); 
  pTargCM.SetPxPyPzE(0.,0.,fBeamEnergy,TMath::Sqrt(fBeamEnergy*fBeamEnergy+mp*mp)); 
  //
  // --- Get the muons parameters in the CM frame 
  //
  pMu1CM.SetPxPyPzE(px1,py1,pz1,e1);
  pMu2CM.SetPxPyPzE(px2,py2,pz2,e2);
  //
  // --- Obtain the dimuon parameters in the CM frame
  //
  pDimuCM=pMu1CM+pMu2CM;
  //
  // --- Translate the dimuon parameters in the dimuon rest frame
  //
  beta=(-1./pDimuCM.E())*pDimuCM.Vect();
  if(beta.Mag()>=1) return 666.;
  pMu1Dimu=pMu1CM;
  pMu2Dimu=pMu2CM;
  pProjDimu=pProjCM;
  pTargDimu=pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  
  //Debugging part -------------------------------------
  Double_t debugProj[4]={0.,0.,0.,0.};
  Double_t debugTarg[4]={0.,0.,0.,0.};
  Double_t debugMu1[4]={0.,0.,0.,0.};
  Double_t debugMu2[4]={0.,0.,0.,0.};
  pMu1Dimu.GetXYZT(debugMu1);
  pMu2Dimu.GetXYZT(debugMu2);
  pProjDimu.GetXYZT(debugProj);
  pTargDimu.GetXYZT(debugTarg);
  if (debugProj[0]!=debugProj[0] ||debugProj[1]!=debugProj[1] || debugProj[2]!=debugProj[2] ||debugProj[3]!=debugProj[3]) return 666; 
  if (debugTarg[0]!=debugTarg[0] ||debugTarg[1]!=debugTarg[1] || debugTarg[2]!=debugTarg[2] ||debugTarg[3]!=debugTarg[3]) return 666; 
  if (debugMu1[0]!=debugMu1[0] ||debugMu1[1]!=debugMu1[1] || debugMu1[2]!=debugMu1[2] ||debugMu1[3]!=debugMu1[3]) return 666; 
  if (debugMu2[0]!=debugMu2[0] ||debugMu2[1]!=debugMu2[1] || debugMu2[2]!=debugMu2[2] ||debugMu2[3]!=debugMu2[3]) return 666; 
  //----------------------------------------------------

  // --- Determine the z axis for the CS angle 
  zaxisCS=(((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
  				     
  // --- Determine the CS angle (angle between mu+ and the z axis defined above)
  Double_t cost;
  
  if(charge1>0) {cost = zaxisCS.Dot((pMu1Dimu.Vect()).Unit());}
  else {cost = zaxisCS.Dot((pMu2Dimu.Vect()).Unit());}
  
  return cost;
}

//________________________________________________________________________
Double_t AliAnalysisTaskDimuonCFContainerBuilder::CostHE(Double_t px1, Double_t py1, Double_t pz1, Double_t e1,
Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2)
{
// Cosine of the theta decay angle (mu+) in the Helicity frame
  
  TLorentzVector pMu1CM, pMu2CM, pDimuCM; // In the CM frame 
  TLorentzVector pMu1Dimu, pMu2Dimu; // In the dimuon rest frame
  TVector3 beta,zaxisCS;
  //
  // --- Get the muons parameters in the CM frame
  //
  pMu1CM.SetPxPyPzE(px1,py1,pz1,e1);
  pMu2CM.SetPxPyPzE(px2,py2,pz2,e2);
  //
  // --- Obtain the dimuon parameters in the CM frame
  //
  pDimuCM=pMu1CM+pMu2CM;
  //
  // --- Translate the muon parameters in the dimuon rest frame
  //
  beta=(-1./pDimuCM.E())*pDimuCM.Vect();
  if(beta.Mag()>=1) return 666.;
  pMu1Dimu=pMu1CM;
  pMu2Dimu=pMu2CM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  
  //Debugging part -------------------------------------
  Double_t debugMu1[4]={0.,0.,0.,0.};
  Double_t debugMu2[4]={0.,0.,0.,0.};
  pMu1Dimu.GetXYZT(debugMu1);
  pMu2Dimu.GetXYZT(debugMu2);
  if (debugMu1[0]!=debugMu1[0] ||debugMu1[1]!=debugMu1[1] || debugMu1[2]!=debugMu1[2] ||debugMu1[3]!=debugMu1[3]) return 666; 
  if (debugMu2[0]!=debugMu2[0] ||debugMu2[1]!=debugMu2[1] || debugMu2[2]!=debugMu2[2] ||debugMu2[3]!=debugMu2[3]) return 666; 
  //----------------------------------------------------
 
  // --- Determine the z axis for the calculation of the polarization angle (i.e. the direction of the dimuon in the CM system)
  TVector3 zaxis;
  zaxis=(pDimuCM.Vect()).Unit();
  
  // --- Calculation of the polarization angle (angle between mu+ and the z axis defined above)
  Double_t cost;
  if(charge1>0) {cost = zaxis.Dot((pMu1Dimu.Vect()).Unit());} 
  else {cost = zaxis.Dot((pMu2Dimu.Vect()).Unit());} 
  return cost;
}

//________________________________________________________________________
Double_t AliAnalysisTaskDimuonCFContainerBuilder::PhiCS(Double_t px1, Double_t py1, Double_t pz1, Double_t e1,
Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2)
{
// Phi decay angle (mu+) in the Collins-Soper frame

   TLorentzVector pMu1CM, pMu2CM, pProjCM, pTargCM, pDimuCM; // In the CM frame
   TLorentzVector pMu1Dimu, pMu2Dimu, pProjDimu, pTargDimu; // In the dimuon rest frame
   TVector3 beta,yaxisCS, xaxisCS, zaxisCS;
   Double_t mp=0.93827231;
   
   // --- Fill the Lorentz vector for projectile and target in the CM frame
   pProjCM.SetPxPyPzE(0.,0.,-fBeamEnergy,TMath::Sqrt(fBeamEnergy*fBeamEnergy+mp*mp)); 
   pTargCM.SetPxPyPzE(0.,0.,fBeamEnergy,TMath::Sqrt(fBeamEnergy*fBeamEnergy+mp*mp)); 
   
   // --- Get the muons parameters in the CM frame 
   pMu1CM.SetPxPyPzE(px1,py1,pz1,e1);
   pMu2CM.SetPxPyPzE(px2,py2,pz2,e2);
   
   // --- Obtain the dimuon parameters in the CM frame
   pDimuCM=pMu1CM+pMu2CM;
   
   // --- Translate the dimuon parameters in the dimuon rest frame
   beta=(-1./pDimuCM.E())*pDimuCM.Vect();
   if(beta.Mag()>=1) return 666.;
   pMu1Dimu=pMu1CM;
   pMu2Dimu=pMu2CM;
   pProjDimu=pProjCM;
   pTargDimu=pTargCM;
   pMu1Dimu.Boost(beta);
   pMu2Dimu.Boost(beta);
   pProjDimu.Boost(beta);
   pTargDimu.Boost(beta);

   //Debugging part -------------------------------------
   Double_t debugProj[4]={0.,0.,0.,0.};
   Double_t debugTarg[4]={0.,0.,0.,0.};
   Double_t debugMu1[4]={0.,0.,0.,0.};
   Double_t debugMu2[4]={0.,0.,0.,0.};
   pMu1Dimu.GetXYZT(debugMu1);
   pMu2Dimu.GetXYZT(debugMu2);
   pProjDimu.GetXYZT(debugProj);
   pTargDimu.GetXYZT(debugTarg);
   if (debugProj[0]!=debugProj[0] ||debugProj[1]!=debugProj[1] || debugProj[2]!=debugProj[2] ||debugProj[3]!=debugProj[3]) return 666; 
   if (debugTarg[0]!=debugTarg[0] ||debugTarg[1]!=debugTarg[1] || debugTarg[2]!=debugTarg[2] ||debugTarg[3]!=debugTarg[3]) return 666; 
   if (debugMu1[0]!=debugMu1[0] ||debugMu1[1]!=debugMu1[1] || debugMu1[2]!=debugMu1[2] ||debugMu1[3]!=debugMu1[3]) return 666; 
   if (debugMu2[0]!=debugMu2[0] ||debugMu2[1]!=debugMu2[1] || debugMu2[2]!=debugMu2[2] ||debugMu2[3]!=debugMu2[3]) return 666; 
   //----------------------------------------------------

   // --- Determine the z axis for the CS angle 
   zaxisCS=(((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
   yaxisCS=(((pProjDimu.Vect()).Unit()).Cross((pTargDimu.Vect()).Unit())).Unit();
   xaxisCS=(yaxisCS.Cross(zaxisCS)).Unit();
 
   Double_t phi=0.;
   if(charge1>0) {
       phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxisCS),((pMu1Dimu.Vect()).Dot(xaxisCS)));
   } else {
       phi = TMath::ATan2((pMu2Dimu.Vect()).Dot(yaxisCS),((pMu2Dimu.Vect()).Dot(xaxisCS)));
   }
   if (phi>TMath::Pi()) phi=phi-TMath::Pi();
   
   return phi;
}

//________________________________________________________________________
Double_t AliAnalysisTaskDimuonCFContainerBuilder::PhiHE(Double_t px1, Double_t py1, Double_t pz1, Double_t e1,
Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2)
{
// Phi decay angle (mu+) in the Helicity frame
  TLorentzVector pMu1Lab, pMu2Lab, pProjLab, pTargLab, pDimuLab; // In the lab. frame 
  TLorentzVector pMu1Dimu, pMu2Dimu, pProjDimu, pTargDimu; // In the dimuon rest frame
  TVector3 beta,xaxis, yaxis,zaxis;
  Double_t mp=0.93827231;

  // --- Get the muons parameters in the LAB frame
  pMu1Lab.SetPxPyPzE(px1,py1,pz1,e1);
  pMu2Lab.SetPxPyPzE(px2,py2,pz2,e2);
  
  // --- Obtain the dimuon parameters in the LAB frame
  pDimuLab=pMu1Lab+pMu2Lab;
  zaxis=(pDimuLab.Vect()).Unit();
  
  // --- Translate the muon parameters in the dimuon rest frame
  beta=(-1./pDimuLab.E())*pDimuLab.Vect();
  if(beta.Mag()>=1.) return 666.;

  pProjLab.SetPxPyPzE(0.,0.,-fBeamEnergy,TMath::Sqrt(fBeamEnergy*fBeamEnergy+mp*mp)); // proiettile
  pTargLab.SetPxPyPzE(0.,0.,fBeamEnergy,TMath::Sqrt(fBeamEnergy*fBeamEnergy+mp*mp)); // bersaglio

  pProjDimu=pProjLab;
  pTargDimu=pTargLab;

  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  
  yaxis=((pProjDimu.Vect()).Cross(pTargDimu.Vect())).Unit();
  xaxis=(yaxis.Cross(zaxis)).Unit();
  
  pMu1Dimu=pMu1Lab;
  pMu2Dimu=pMu2Lab;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  
  //Debugging part -------------------------------------
  Double_t debugProj[4]={0.,0.,0.,0.};
  Double_t debugTarg[4]={0.,0.,0.,0.};
  Double_t debugMu1[4]={0.,0.,0.,0.};
  Double_t debugMu2[4]={0.,0.,0.,0.};
  pMu1Dimu.GetXYZT(debugMu1);
  pMu2Dimu.GetXYZT(debugMu2);
  pProjDimu.GetXYZT(debugProj);
  pTargDimu.GetXYZT(debugTarg);
  if (debugProj[0]!=debugProj[0] ||debugProj[1]!=debugProj[1] || debugProj[2]!=debugProj[2] ||debugProj[3]!=debugProj[3]) return 666; 
  if (debugTarg[0]!=debugTarg[0] ||debugTarg[1]!=debugTarg[1] || debugTarg[2]!=debugTarg[2] ||debugTarg[3]!=debugTarg[3]) return 666; 
  if (debugMu1[0]!=debugMu1[0] ||debugMu1[1]!=debugMu1[1] || debugMu1[2]!=debugMu1[2] ||debugMu1[3]!=debugMu1[3]) return 666; 
  if (debugMu2[0]!=debugMu2[0] ||debugMu2[1]!=debugMu2[1] || debugMu2[2]!=debugMu2[2] ||debugMu2[3]!=debugMu2[3]) return 666; 
  //----------------------------------------------------
  
  Double_t phi=0.;
   if(charge1 > 0) {
      phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxis),(pMu1Dimu.Vect()).Dot(xaxis));
     } else { 
      phi = TMath::ATan2((pMu2Dimu.Vect()).Dot(yaxis),(pMu2Dimu.Vect()).Dot(xaxis));
   }  
   return phi;
}

//________________________________________________________________________
void AliAnalysisTaskDimuonCFContainerBuilder::Terminate(Option_t *) 
{
 
}

#endif
