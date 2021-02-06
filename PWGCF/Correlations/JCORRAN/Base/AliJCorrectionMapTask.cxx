/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

//______________________________________________________________________________
// Analysis task for high pt particle correlations
// author: Jasper Parkila, D.J. Kim
// ALICE Group University of Jyvaskyla
// Finland
//////////////////////////////////////////////////////////////////////////////

#include <TGrid.h>
#include <TRandom.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <AliAnalysisTaskSE.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisDataContainer.h>
#include "AliJCorrectionMapTask.h"
#include "AliLog.h"

ClassImp(AliJCorrectionMapTask)

//#pragma GCC diagnostic warning "-Wall"
//______________________________________________________________________________
AliJCorrectionMapTask::AliJCorrectionMapTask():
	AliAnalysisTaskSE(),
	//fOutputList(NULL),
	inputIndex(1)
{
	//
}

//______________________________________________________________________________
AliJCorrectionMapTask::AliJCorrectionMapTask(const char *name):
	AliAnalysisTaskSE(name),
	//fOutputList(NULL),
	inputIndex(1)
{
	//DefineOutput(1, TDirectory::Class());
}

//____________________________________________________________________________
AliJCorrectionMapTask::AliJCorrectionMapTask(const AliJCorrectionMapTask& ap) :
	AliAnalysisTaskSE(ap.GetName()),
	//fOutputList(ap.fOutputList),
	inputIndex(ap.inputIndex)
	//fOutput(ap.fOutput)
{
	AliInfo("----DEBUG AliJCorrectionMapTask COPY ----");
}

//_____________________________________________________________________________
AliJCorrectionMapTask& AliJCorrectionMapTask::operator = (const AliJCorrectionMapTask& ap)
{
	// assignment operator
	AliInfo("----DEBUG AliJCorrectionMapTask operator= ----");
	this->~AliJCorrectionMapTask();
	new(this) AliJCorrectionMapTask(ap);
	return *this;
}

//______________________________________________________________________________
AliJCorrectionMapTask::~AliJCorrectionMapTask()
{
	//delete fOutputList; 
}

//________________________________________________________________________
void AliJCorrectionMapTask::UserCreateOutputObjects()
{
	//OpenFile(1);
	//if(!fOutputList) fOutputList = new TList();
	//fOutputList->SetOwner(kTRUE);

	//Bool_t oldStatus = TH1::AddDirectoryStatus();
	//TH1::AddDirectory(kFALSE);

	//PostData(1, fOutputList);
}

//______________________________________________________________________________
void AliJCorrectionMapTask::UserExec(Option_t* /*option*/)
{
	// Processing of one event
	//if(!((Entry()-1)%100))  AliInfo(Form(" Processing event # %lld",  Entry()));
	//PostData(1, fOutputList);

}

void AliJCorrectionMapTask::Init()
{
	AliInfo("Doing initialization") ;
}

//______________________________________________________________________________
void AliJCorrectionMapTask::Terminate(Option_t *)
{
	//fFFlucAna->Terminate("");
	cout<<"AliJCorrectionMapTask Analysis DONE !!"<<endl;
	if(!GetOutputData(1)) return;
}

UInt_t AliJCorrectionMapTask::ConnectInputContainer(const TString fname, const TString listName){
	DefineInput(inputIndex,TList::Class());

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
		
	TString containerName = Form("CorrectionMap-%u",fname.Hash());
	cout<<"Container: "<<containerName<<endl;

	TObjArray *ptaskContainer = mgr->GetContainers();
	AliAnalysisDataContainer *pCorrMapCont = (AliAnalysisDataContainer*)ptaskContainer->FindObject(containerName);
	if(!pCorrMapCont){
		TGrid::Connect("alien:");
		TFile *pfile = TFile::Open(fname);
		TList *plist = (TList*)pfile->Get(listName);
		pCorrMapCont = mgr->CreateContainer(containerName,
			TList::Class(),AliAnalysisManager::kInputContainer);
		pCorrMapCont->SetData(plist);
	}
	
	mgr->ConnectInput(this,inputIndex,pCorrMapCont);

	return inputIndex++; // need this!!!
}

void AliJCorrectionMapTask::EnablePhiCorrection(UInt_t id, const TString fname){
	phiInputIndex[id] = ConnectInputContainer(fname,Form("PhiWeights_%u",id));
	cout<<"Phi correction enabled: "<<fname.Data()<<" (index "<<phiInputIndex[id]<<")"<<endl;
}

void AliJCorrectionMapTask::EnableCentFlattening(const TString fname){
	centInputIndex = ConnectInputContainer(fname,"CentralityWeights");//inputIndex++;
	cout<<"Centrality flattening enabled: "<<fname.Data()<<" (index "<<centInputIndex<<")"<<endl;
}

void AliJCorrectionMapTask::EnableEffCorrection(const TString fname){
	effInputIndex = ConnectInputContainer(fname,"EffCorrections");//inputIndex++;
	cout<<"Efficiency Correction enabled: "<<fname.Data()<<" (index "<<effInputIndex<<")"<<endl;
}

TH1 * AliJCorrectionMapTask::GetCorrectionMap(UInt_t it, UInt_t run, UInt_t cbin){ // mapID, run#, cent
	if(cbin>8) cbin = 8; // Temp fix since >60% centrality is missing in the map
	auto m = PhiWeightMap[it][cbin].find(run);
	if(m == PhiWeightMap[it][cbin].end()){
		TList *plist = (TList*)GetInputData(phiInputIndex[it]);
		if(!plist)
			return 0;
		TH1 *pmap = (TH1*)plist->FindObject(Form("PhiWeights_%u_%02u",run,cbin));
		if(!pmap)
			return 0;
		PhiWeightMap[it][cbin][run] = pmap;
		return pmap;
	}

	return (*m).second;
}

TH1 * AliJCorrectionMapTask::GetCentCorrection(){
	TList *plist = (TList*)GetInputData(centInputIndex);
	if(!plist)
		return 0;
	TH1 *pmap = (TH1*)plist->FindObject("CentCorrection");
	return pmap;
}

TAxis * AliJCorrectionMapTask::GetCentBinEff() {
	TList *plist = (TList*)GetInputData(effInputIndex); 
	if(!plist) {
		AliError("ERROR: effInputIndex not available");
	}
	fCentBinEff = (TAxis*)plist->FindObject("CentralityBin");
	return fCentBinEff;
}

// Pt dependent efficiency loaded here to improve the CPU time.
TGraphErrors * AliJCorrectionMapTask::GetEffCorrectionMap(UInt_t run, Double_t cent, UInt_t fEffFilterBit){ //need to be centrality itself since efficiency is merged in bins
	int centBin = fCentBinEff->FindBin( cent ) -1 ;
	auto m = EffWeightMap[centBin].find(run);
	if(m == EffWeightMap[centBin].end()){
		TList *plist = (TList*)GetInputData(effInputIndex); 
		if(!plist)
			return 0;
		  
		if( centBin < 0 || centBin > fCentBinEff->GetNbins()-1 ) {
			cout<<"J_WARNING : Centrality "<<cent<<" is out of CentBinBorder"<<endl;
			return 0;
		}
		
		TGraphErrors *grmap = (TGraphErrors*)plist->FindObject(Form("gCor%02d%02d%02d", 0,centBin,fEffFilterBit));
		
		EffWeightMap[centBin][run] = grmap;
		return grmap;
	}
	return (*m).second;
}

// Pt dependent efficiency loaded here to improve the CPU time.
double AliJCorrectionMapTask::GetEffCorrection( TGraphErrors * gr, double pt ) const {
	// TODO : Function mode
	//=== TEMPERORY SETTING. IT will be removed soon.
	if( pt > 30 ) pt = 30; // Getting eff of 30GeV for lager pt
	double cor = gr->Eval(pt);
	if ( cor < 0.2 ) cor = 0.2;
	return cor;
}
