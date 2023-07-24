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

UInt_t AliJCorrectionMapTask::ConnectInputContainer(const TString &fname, const TString &listName){
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
		if(!plist)
			cout<<"ERROR: Correction map TList does not exists: "<<listName.Data()<<endl;
		pCorrMapCont = mgr->CreateContainer(containerName,
			TList::Class(),AliAnalysisManager::kInputContainer);
		pCorrMapCont->SetData(plist);
	}
	
	mgr->ConnectInput(this,inputIndex,pCorrMapCont);

	return inputIndex++; // need this!!!
}

void AliJCorrectionMapTask::EnablePhiCorrection(UInt_t id, const TString &fname){
	phiInputIndex[id] = ConnectInputContainer(fname,"PhiWeights");
	cout<<"Phi correction enabled: "<<fname.Data()<<" (id "<<id<<", index "<<phiInputIndex[id]<<")"<<endl;
}

void AliJCorrectionMapTask::EnableCentFlattening(const TString &fname){
	centInputIndex = ConnectInputContainer(fname,"CentralityWeights");//inputIndex++;
	cout<<"Centrality flattening enabled: "<<fname.Data()<<" (index "<<centInputIndex<<")"<<endl;
}

void AliJCorrectionMapTask::EnableEffCorrection(const TString &fname){
	effInputIndex = ConnectInputContainer(fname,"EffCorrections");//inputIndex++;
	cout<<"Efficiency correction enabled: "<<fname.Data()<<" (index "<<effInputIndex<<")"<<endl;
}

void AliJCorrectionMapTask::EnableEffCorrection2(const TString &fname){
	effInputIndex2 = ConnectInputContainer(fname,"EffCorrections");//inputIndex++;
	cout<<"Efficiency correction2 enabled: "<<fname.Data()<<" (index "<<effInputIndex2<<")"<<endl;
	runPeriods = {
		{252235,252330,"LHC16d",120.0},
		{253437,253591,"LHC16e",120.0},
		{253659,253978,"LHC16f",120.0},
		{254128,254332,"LHC16g",120.0},
		{254604,255467,"LHC16h",120.0},
		{255539,255618,"LHC16i",120.0},
		{256219,256418,"LHC16j",120.0},
		{256941,256219,"LHC16k",120.0},
		{258962,259888,"LHC16l",89.9003},// V0M_mean=89.9003; }
		{262424,264035,"LHC16o",86.3912},// V0M_mean=86.3912; }
		{264076,264347,"LHC16p",138.814},// V0M_mean=138.814; }

		{270581,270667,"LHC17c",120.0},
		{270822,270830,"LHC17e",120.0},
		{270854,270865,"LHC17f",120.0},
		{270882,271777,"LHC17g",120.0},
		{271870,273103,"LHC17h",127.895},// V0M_mean=127.895; }
		{273591,274442,"LHC17i",124.276},// V0M_mean=124.276; }
		{274593,274671,"LHC17j",120.0},
		{274690,276508,"LHC17k",121.31},// V0M_mean=121.31; }
		{276551,278216,"LHC17l",119.144},// V0M_mean=119.144; }
		{278914,280140,"LHC17m",117.165},// V0M_mean=117.165; }
		{280282,281961,"LHC17o",113.45},// V0M_mean=113.45; }
		{282528,282704,"LHC17r",111.462},// V0M_mean=111.462; }

		{285009,285396,"LHC18b",120.0},
		{285978,286350,"LHC18d",131.868},// V0M_mean=131.868; }
		{286380,286937,"LHC18e",131.397},// V0M_mean=131.397; }
		{287000,287658,"LHC18f",130.591},// V0M_mean=130.591; }
		{288750,288619,"LHC18g",120.0},
		{288806,288804,"LHC18h",130.86},//V0M_mean=130.86; }
		{288909,288861,"LHC18i",120.0},
		{288943,288943,"LHC18j",131.17},//V0M_mean=131.17; }
		{289240,289971,"LHC18l",131.59},//V0M_mean=131.59; }
		{290323,292839,"LHC18m",130.467},//V0M_mean=130.467; }
		{293359,293357,"LHC18n",120.0},
		{289201,289165,"LHC18k",127.642},//V0M_mean=127.642; }
		{293475,293898,"LHC18o",124.973},//V0M_mean=124.973; }
		{294009,294925,"LHC18p",120.0}
	};
}

TH3D * AliJCorrectionMapTask::GetCorrectionMap(UInt_t it, UInt_t run, UInt_t cbin){ // mapID, run#, cent
	if(cbin>8) cbin = 8; // Temp fix since >60% centrality is missing in the map
	auto m = PhiWeightMap[it][cbin].find(run);
	if(m == PhiWeightMap[it][cbin].end()){
		TList *plist = (TList*)GetInputData(phiInputIndex[it]);
		if(!plist)
			return 0;
		TH3D *pmap = (TH3D*)plist->FindObject(Form("PhiWeights_%u_%02u",run,cbin));
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
		
		TGraphErrors *grmap = (TGraphErrors*)plist->FindObject(Form("gCor%02d%02d%02d",0,centBin,fEffFilterBit));
		
		EffWeightMap[centBin][run] = grmap;
		return grmap;
	}
	return (*m).second;
}

// Pt dependent efficiency loaded here to improve the CPU time.
double AliJCorrectionMapTask::GetEffCorrection(TGraphErrors *gr, double pt) const{
	// TODO : Function mode
	//=== TEMPERORY SETTING. IT will be removed soon.
	if( pt > 30 ) pt = 30; // Getting eff of 30GeV for lager pt
	double cor = gr->Eval(pt);
	if ( cor < 0.2 ) cor = 0.2;
	return cor;
}

//std::tuple<TH1 *, double> AliJCorrectionMapTask::GetEffCorrectionMap2(UInt_t run, EFF2_LABEL effLabel){
TH1 * AliJCorrectionMapTask::GetEffCorrectionMap2(UInt_t run, EFF2_LABEL effLabel, double &V0mean){
	auto m = std::find_if(runPeriods.begin(),runPeriods.end(),[&](/*auto &t*/const RunPeriod &t)->bool{
		//return std::get<0>(t) <= run && run <= std::get<1>(t);
		return t.runStart <= run && run <= t.runEnd;
	});
	const char *pPeriod;
	if(m != runPeriods.end()){
		pPeriod = (*m).pPeriod;
		V0mean = (*m).V0mean;
	}else{
		pPeriod = "LHC16q";
		V0mean = 120.0;
	}
	//const char *pPeriod = m != runPeriods.end()?
	//	std::get<2>(*m):"LHC16q"; //assume pPb if not on the pp list
	
	TList *plist = (TList*)GetInputData(effInputIndex2);
	if(!plist)
		return 0;//std::make_tuple((TH1*)0,V0mean);
	
	TH1 *pmap = (TH1*)plist->FindObject(Form("%s_%s",pPeriod,peff2Labels[effLabel]));
	//return std::make_tuple(pmap,V0mean);
	return pmap;
}

TH1 * AliJCorrectionMapTask::GetEffCorrectionMap2(const TString &tag, EFF2_LABEL effLabel){
	TList *plist = (TList*)GetInputData(effInputIndex2);
	if(!plist)
		return 0;
	
	TH1 *pmap = (TH1*)plist->FindObject(Form("%s_%s",tag.Data(),peff2Labels[effLabel]));
	return pmap;
}

//double AliJCorrectionMapTask::GetEffCorrection2(TH1 *ph, double pt, 

const char * AliJCorrectionMapTask::peff2Labels[AliJCorrectionMapTask::EFF2_LABEL_COUNT] = {
	"Glb8cm",
	"GlbSDD8cm",
	"Hyb6cm",
	"Hyb8cm",
	"Hyb10cm"
};

