/*************************************************************************
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

//  ------------------------------------------
//  analysis task for azimuthal isotropic
//  expansion in highly central collisions
//  author: Cristian Andrei
//          acristian@niham.nipne.ro
//  ------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TList.h"
#include "TObjArray.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "AliAnalysisManager.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliCFContainer.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisDataSlot.h"

#include "AliAnalysisCentralCutMC.h"
#include "AliAnalysisCentralCutESD.h"
#include "AliAnalysisCentralCutEvtMC.h"
#include "AliAnalysisCentralCutEvtESD.h"
#include "AliAnalysisCentralExtrapolate.h"
#include "AliAnalysisTaskCentral.h"


ClassImp(AliAnalysisTaskCentral)


//________________________________________________________________________
AliAnalysisTaskCentral::AliAnalysisTaskCentral(const char *name) 
  : AliAnalysisTask(name, "")
  ,fESD(0), fMC(0)
  ,fNoEvt(0)
  ,fCFContainerPi(0)
  ,fCFContainerK(0)
  ,fCFContainerP(0)
  ,fSim(kFALSE)
  ,fOutList(NULL)

{
  // Constructor
	printf("AliAnalysisTaskCentral::AliAnalysisTaskCentral(const char *name)\n");

  // Define input and output slots here
	DefineInput(0, TChain::Class());

	DefineOutput(0, TList::Class()); 


	for(Int_t i=0; i<10; i++){
		fCutsList[i] = 0;
    }

	InitCuts(); //initialize the analysis specific cuts	

}


//________________________________________________________________________
AliAnalysisTaskCentral::~AliAnalysisTaskCentral() 
{
// Destructor
// Delete the created objects

	if(fESD) delete fESD;
	if(fMC) delete fMC;

	for (Int_t i=0; i<10; i++)
	  delete fCutsList[i];


	if(fNoEvt) delete fNoEvt;

	if(fCFContainerPi) delete fCFContainerPi;
	if(fCFContainerPi) delete fCFContainerK;
	if(fCFContainerPi) delete fCFContainerP;

	if(fOutList) delete fOutList;


}

//________________________________________________________________________
void AliAnalysisTaskCentral::ConnectInputData(Option_t *) {
// get the event from the input chain

    TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
	printf("tree(%p)\n", (void*)tree);

    if (!tree) {
	Printf("ERROR: Could not read chain from input slot 0");
    } else {
	AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
	if (!esdH) {
    	    Printf("ERROR: Could not get ESDInputHandler");
	}
	else fESD = esdH->GetEvent();

	AliMCEventHandler *mcH = (AliMCEventHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
	if (!mcH) {
	    Printf("ERROR: Could not get MCInputHandler");
	    fSim = kFALSE;
	}
	else{
	    fMC = mcH->MCEvent();
	    fSim = kTRUE;
	}
    }
}

//________________________________________________________________________
void AliAnalysisTaskCentral::CreateOutputObjects(){
//Creates the output list

// 	fNoEvt = new TH1F(Form("%s_NoEvt",GetName()), "Number of processed events", 1, 0.0, 1.0);

	fNoEvt = new TH1D("TaskCentral_NoEvt", "Number of processed events", 1, 0.0, 1.0);
//set here the min and max variable values (to do: move these to a separate cuts class)
    const Double_t ptmin  = 0.0 ; //bins in pt
    const Double_t ptmax  = 5.0 ; //bins in y
    const Double_t etamin  = -0.9 ; //bins in pt
    const Double_t etamax  = 0.9 ; //bins in y
    //apropiate number of bins for histos and CF Container
    const Int_t nbinpt  = 25 ; //bins in pt
    const Int_t nbineta  = 25 ; //bins in y

//Correction Framework CONTAINER DEFINITION
    //the sensitive variables, their indices
	UInt_t ipt = 0;
	UInt_t ieta  = 1;
    //Setting up the container grid...
    const UInt_t nstep = 2 ; //number of selection steps MC & ESD 
    const Int_t nvar = 2 ; //number of variables on the grid:pt,eta

    //arrays for the number of bins in each dimension
    //nbin is defined above
    Int_t iBin[nvar];
    iBin[0]=nbinpt;
    iBin[1]=nbineta;

    //arrays for lower bounds :
    Double_t *binLim1=new Double_t[nbinpt+1];
    Double_t *binLim2=new Double_t[nbineta+1];

    //values for bin lower bounds
    for(Int_t i=0; i<=nbinpt; i++) binLim1[i]=(Double_t)ptmin + (ptmax-ptmin)/nbinpt*(Double_t)i ;
    for(Int_t i=0; i<=nbineta; i++) binLim2[i]=(Double_t)etamin  + (etamax-etamin)  /nbineta*(Double_t)i ;

	//container creation
	fCFContainerPi = new AliCFContainer("TaskCentral_CFCont_Pi","container for tracks",nstep,nvar,iBin);
	fCFContainerK = new AliCFContainer("TaskCentral_CFCont_K","container for tracks",nstep,nvar,iBin);
	fCFContainerP = new AliCFContainer("TaskCentral_CFCont_P","container for tracks",nstep,nvar,iBin);

    //setting the bin limits
    fCFContainerPi -> SetBinLimits(ipt,binLim1);
    fCFContainerPi -> SetBinLimits(ieta,binLim2);

	fCFContainerK -> SetBinLimits(ipt,binLim1);
	fCFContainerK -> SetBinLimits(ieta,binLim2);

	fCFContainerP -> SetBinLimits(ipt,binLim1);
	fCFContainerP -> SetBinLimits(ieta,binLim2);

//outout list creation
	fOutList = new TList();
	fOutList->Add(fNoEvt);
	fOutList->Add(fCFContainerPi);
	fOutList->Add(fCFContainerK);
	fOutList->Add(fCFContainerP);
}


//________________________________________________________________
void AliAnalysisTaskCentral::InitCuts(){

//Create and set cuts

//------------------EVENT LEVEL CUTS--------------------
//-----------------------MC-----------------------------
	AliAnalysisCentralCutEvtMC *evMCCuts = new AliAnalysisCentralCutEvtMC();
	evMCCuts->SetMultiplicityRange(1,100);
	evMCCuts->SetDirectivityRange(0.0, 1.0);
    evMCCuts->SetDirUnitRange(0.0, 1.0);

//-----------------------ESD----------------------------
	AliAnalysisCentralCutEvtESD *evESDCuts = new AliAnalysisCentralCutEvtESD();
	evESDCuts->SetMultiplicityRange(1,100);
	evESDCuts->SetDirectivityRange(0.0, 1.0);
    evESDCuts->SetSPDMultiplicityRange(1,20000);
    evESDCuts->SetSPDDirectivityRange(0.0, 1.0);

	TObjArray* mcEventCuts = new TObjArray();
	mcEventCuts->AddLast(evMCCuts);

	TObjArray* esdEventCuts = new TObjArray();
	esdEventCuts->AddLast(evESDCuts);


//------------------PARTICLE LEVEL CUTS-----------------
//-------------------General MC Cuts--------------------
	AliAnalysisCentralCutMC *mcCutsGen = new AliAnalysisCentralCutMC();
	mcCutsGen->SetOnlyPrimaries(kTRUE);
	mcCutsGen->SetPtRange(0.2,4.0);
	mcCutsGen->SetEtaRange(-0.2,0.2);

//-------------------Specific MC Cuts-------------------
	AliAnalysisCentralCutMC *mcCutsPi = new AliAnalysisCentralCutMC();
	mcCutsPi->SetPDGCode(211); //211 pion;  321 kaon; 2212 proton

	AliAnalysisCentralCutMC *mcCutsK = new AliAnalysisCentralCutMC();
	mcCutsK->SetPDGCode(321); //211 pion;  321 kaon; 2212 proton

	AliAnalysisCentralCutMC *mcCutsP = new AliAnalysisCentralCutMC();
	mcCutsP->SetPDGCode(2212); //211 pion;  321 kaon; 2212 proton

//--------each task has its own mcList of cuts----------
	TObjArray *mcListGen = new TObjArray(); //general MC cuts
	mcListGen->AddLast(mcCutsGen);
	
	TObjArray *mcListPi = new TObjArray(); //task pt pions
	mcListPi->AddLast(mcCutsPi);

	TObjArray *mcListK = new TObjArray(); //task pt kaons
	mcListK->AddLast(mcCutsK);

	TObjArray *mcListP = new TObjArray(); //task pt protons
	mcListP->AddLast(mcCutsP);


//-------------------General ESD Cuts-------------------
	AliESDtrackCuts *esdCutsGen = new AliESDtrackCuts("AliESDtrackCuts", "Loose");
	esdCutsGen->SetMinNClustersTPC(50);
	esdCutsGen->SetMaxChi2PerClusterTPC(2.2);
	esdCutsGen->SetMaxCovDiagonalElements(0.5,0.5,0.5,0.5,0.5);
	esdCutsGen->SetRequireTPCRefit(kTRUE);
	esdCutsGen->SetAcceptKinkDaughters(kFALSE);
	esdCutsGen->SetMaxNsigmaToVertex(2.0);
	esdCutsGen->SetRequireSigmaToVertex(kTRUE);
	esdCutsGen->SetPtRange(0.2,4.0);
	esdCutsGen->SetEtaRange(-0.2,0.2);

//-------------------Specific ESD Cuts------------------
	AliAnalysisCentralCutESD *esdCutsPi = new AliAnalysisCentralCutESD("AliAnalysisCentralCutESD","NIHAM");
	esdCutsPi->SetPIDtype("Custom");
	esdCutsPi->SetPriorFunctions(kFALSE);
	esdCutsPi->SetPartType(kPiPlus);

	AliAnalysisCentralCutESD *esdCutsK = new AliAnalysisCentralCutESD("AliAnalysisCentralCutESD","NIHAM");
	esdCutsK->SetPIDtype("Custom");
	esdCutsK->SetPriorFunctions(kFALSE);
	esdCutsK->SetPartType(kKPlus);

	AliAnalysisCentralCutESD *esdCutsP = new AliAnalysisCentralCutESD("AliAnalysisCentralCutESD","NIHAM");
	esdCutsP->SetPIDtype("Custom");
	esdCutsP->SetPriorFunctions(kFALSE);
	esdCutsP->SetPartType(kProton);

//--------each task has its own esdList of cuts---------
	TObjArray* esdListGen = new TObjArray(); //general ESD track cuts
	esdListGen->AddLast(esdCutsGen);
	
	TObjArray* esdListPi = new TObjArray();
	esdListPi->AddLast(esdCutsPi);

	TObjArray* esdListK = new TObjArray();
	esdListK->AddLast(esdCutsK);

	TObjArray* esdListP = new TObjArray();
	esdListP->AddLast(esdCutsP);

//------set the cuts to the RIGHT! fCutsList slots-------
// event level cuts
	SetCuts(0, mcEventCuts);
	SetCuts(1, esdEventCuts);

// particle level cuts
	SetCuts(2, mcListGen);
	SetCuts(3, mcListPi);
	SetCuts(4, mcListK);
	SetCuts(5, mcListP);

	SetCuts(6, esdListGen);
	SetCuts(7, esdListPi);
	SetCuts(8, esdListK);
	SetCuts(9, esdListP);

}


//_______________________________________________________________________
void AliAnalysisTaskCentral::SendEvent(TObject *obj) const{

// Some cuts (ie MC IsPrimary) need the MC Event info

		for(Int_t isel=0;isel< 10; isel++){
			if(!fCutsList[isel]) continue;
			TObjArrayIter iter(fCutsList[isel]);
			AliAnalysisCuts *cut = 0;
	
			while ((cut = (AliAnalysisCuts*)iter.Next())) {
				
				TString cutName=cut->GetName();
				if(!cutName){
					printf("No cutname!\n");
					return;
				}
		
				Bool_t checkCut=cutName.Contains("AliAnalysisCentralCutMC");
					
				if(checkCut){
					AliAnalysisCentralCutMC *newcut = dynamic_cast<AliAnalysisCentralCutMC *>(cut);
					if (newcut) newcut->ReceiveEvt(obj);
				}
			}
		}

}


//________________________________________________________________________
Bool_t AliAnalysisTaskCentral::CheckCuts(Int_t no, TObject *obj) const{ 

// For each cut run IsSelected();
// 	printf("AliAnalysisTaskCentral::CheckCuts IN\n");

    if(no > 9){
		printf("\nAliAnalysisTaskCentral::CheckCuts -> Cut number is not ok! \n");
		return kFALSE;
    }

    if(!fCutsList[no]){
		printf("AliAnalysisTaskCentral::CheckCuts -> cuts list problem! \n");
		return kFALSE;
    }

    TObjArrayIter iter(fCutsList[no]);
    AliAnalysisCuts *cut = 0;

    while((cut = (AliAnalysisCuts*)iter.Next())){
		if(!cut->IsSelected(obj)){
// 			printf("AliAnalysisTaskCentral::CheckCuts OUT\n");
			return kFALSE;
		}
    }

// 	printf("AliAnalysisTaskCentral::CheckCuts OUT\n");
    return kTRUE;
}


//________________________________________________________________________
void AliAnalysisTaskCentral::Exec(Option_t *) {

// Main loop
// Called for each event

	Double_t pt, eta;
	Double_t w = 1.0;
	const Int_t nvar=2; //number of variables on the grid:pt,vtx
	Double_t value[nvar];//used to fill the CF container

	if(fSim){  // if running on simulations -> look at MC Truth
		
		if (!fMC) {
			Printf("ERROR: fMC not available");
			return;
		}

		if(CheckCuts(0, fMC)){ //check event level cuts
	
			SendEvent(fMC);
		
			// MC loop
			for (Int_t ipart=0; ipart<fMC->GetNumberOfTracks(); ipart++) {
		
				AliMCParticle *particle  = dynamic_cast<AliMCParticle*>(fMC->GetTrack(ipart));
						
				if(!particle){
						printf("\nMCParticle pointer is null!!!\n");
						continue;
				}
				
				if (!CheckCuts(2, particle)) continue; //check the MC general particle cuts 
				
				pt = particle->Pt();
				eta = particle->Eta();
				
				if(pt>0) w = 1.0/pt; //invariant distribution
				
				value[0]=pt;
				value[1]=eta;
				
				if(CheckCuts(3, particle)){  //fill the right container for each particle
					fCFContainerPi->Fill(value,0,w);
				}
				else if(CheckCuts(4, particle)){
					fCFContainerK->Fill(value,0,w);
				}
				else if(CheckCuts(5, particle)){
					fCFContainerP->Fill(value,0,w);
				}
			} //end MC particle loop 
		}
	}
	else{  // if we DONT run in simulated data we fill the MC step of the CFCont with 0
		value[0]=0;
		value[1]=0;
			
		fCFContainerPi->Fill(value,0);
		fCFContainerK->Fill(value,0);
		fCFContainerP->Fill(value,0);
	}
	
	
    if (!fESD) {
	Printf("ERROR: fESD not available");
	return;
    }

    if(CheckCuts(1, fESD)){

		Printf("There are %d ESD tracks in this event\n", fESD->GetNumberOfTracks());
	
		// ESD loop 
		for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
	
			AliESDtrack* track = fESD->GetTrack(iTracks);
			
			if (!track) {
					Printf("ERROR: Could not receive track %d", iTracks);
				continue;
			}
		
			if(!CheckCuts(6, track)) continue; //check general ESD track cuts
		
			pt = track->Pt();
			eta = track->Eta();
		
			if(pt>0) w = 1.0/pt; // invariant
			
			value[0]=pt;
			value[1]=eta;
		
			if(CheckCuts(7, track)){
				fCFContainerPi->Fill(value,1,w);
			}
			else if(CheckCuts(8, track)){
				fCFContainerK->Fill(value,1,w);
			}
			else if(CheckCuts(9, track)){
				fCFContainerP->Fill(value,1,w);
			}
			
		} //end ESD track loop 
		
		fNoEvt->Fill(0); //get the number of analyzed events
	}

  // Post output data.
	PostData(0, fOutList);

}

//________________________________________________________________________
void AliAnalysisTaskCentral::Terminate(Option_t *) {

// Called once at the end of the query
	printf("\n\n****************************************\n");
    printf("\tAliAnalysisCentralExtrapolate Terminate... \n");

	TList *outList = dynamic_cast<TList*>(GetOutputData(0));
	if(!outList){
		printf("Unable to get output list! \n");
		return;
	}

	AliAnalysisCentralExtrapolate *extPi = new AliAnalysisCentralExtrapolate("extrapolationpi");
	extPi->SetInputList(outList);  //send the outlist to the extrapolation class
	extPi->SetParticle("kPiPlus"); //set the particle type
	extPi->ApplyEff();             //correct the pt distribution !!HAS TO RUN BEFORE extrapolation!!
	extPi->BoltzmannFit();         //fit and extrapolate using Boltzmann-Gibbs Blast wave model
	extPi->TsallisFit();           //fit and extrapolate using Tsallis Blast wave model
	TList *extOutListPi = extPi->GetOutputList();

	AliAnalysisCentralExtrapolate *extK = new AliAnalysisCentralExtrapolate("extrapolationK");
	extK->SetInputList(outList);
	extK->SetParticle("kKPlus");
	extK->ApplyEff();
	extK->BoltzmannFit();
	extK->TsallisFit();
	TList *extOutListK = extK->GetOutputList();

	AliAnalysisCentralExtrapolate *extP = new AliAnalysisCentralExtrapolate("extrapolationP");
	extP->SetInputList(outList);
	extP->SetParticle("kProton");
	extP->ApplyEff();
	extP->BoltzmannFit();
	extP->TsallisFit();
	TList *extOutListP = extP->GetOutputList();


//----------- Plot the extrapolated spectra -----------------
	TCanvas *ccorrdata = new TCanvas();
	ccorrdata->Divide(3,2);
	
	ccorrdata->cd(1);
	ccorrdata->cd(1)->SetLogy();
	TH1D *extBoltzPi = dynamic_cast<TH1D*>(extOutListPi->FindObject("PtExtBoltzmann"));
	if(extBoltzPi){
		extBoltzPi->Draw("p e1");
	}

	ccorrdata->cd(4);
	ccorrdata->cd(4)->SetLogy();
	TH1D *extTsallisPi = dynamic_cast<TH1D*>(extOutListPi->FindObject("PtExtTsallis"));
	if(extTsallisPi){
		extTsallisPi->Draw("p e1");
	}


	ccorrdata->cd(2);
	ccorrdata->cd(2)->SetLogy();
	TH1D *extBoltzK = dynamic_cast<TH1D*>(extOutListK->FindObject("PtExtBoltzmann"));
	if(extBoltzK){
		extBoltzK->Draw("p e1");
	}

	ccorrdata->cd(5);
	ccorrdata->cd(5)->SetLogy();
	TH1D *extTsallisK = dynamic_cast<TH1D*>(extOutListK->FindObject("PtExtTsallis"));
	if(extTsallisK){
		extTsallisK->Draw("p e1");
	}


	ccorrdata->cd(3);
	ccorrdata->cd(3)->SetLogy();
	TH1D *extBoltzP = dynamic_cast<TH1D*>(extOutListP->FindObject("PtExtBoltzmann"));
	if(extBoltzP){
		extBoltzP->Draw("p e1");
	}

	ccorrdata->cd(6);
	ccorrdata->cd(6)->SetLogy();
	TH1D *extTsallisP = dynamic_cast<TH1D*>(extOutListP->FindObject("PtExtTsallis"));
	if(extTsallisP){
		extTsallisP->Draw("p e1");
	}

// ------------------ Save the results -----------------------
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if(!mgr){
		printf("Unable to get AnalysisManager! \n");
        return;
	}

	TString anaType;
	mgr->GetAnalysisTypeString(anaType);

	if(anaType.Contains("local")){
		fOutList->Add(extOutListPi);
		fOutList->Add(extOutListK);
		fOutList->Add(extOutListP);
	}
	else{
		AliAnalysisDataContainer *cont = GetOutputSlot(0)->GetContainer();
		if(!cont){
			printf("Unable to get DataContainer! \n");
			return;
		}

		printf("file name = %s\n", cont->GetFileName());
		TFile file(cont->GetFileName(),"update");

		file.cd("PWG2Central");

		gFile->WriteObject(extOutListPi,"pion_list","SingleKey");
		gFile->WriteObject(extOutListK,"kaon_list","SingleKey");
		gFile->WriteObject(extOutListP,"proton_list","SingleKey");
		file.Close();
	}

}
