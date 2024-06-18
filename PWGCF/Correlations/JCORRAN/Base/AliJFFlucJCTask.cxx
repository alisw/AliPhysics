
#include <AliAnalysisManager.h>
//#include <TRandom.h>
//#include <TFile.h>
#include <TGrid.h>
#include <TFile.h>
#include <TRandom.h>
#include <AliAnalysisTaskSE.h>
#include <AliAODHandler.h>
#include <AliAODMCParticle.h>
#include <AliMCEvent.h>
#include <AliGenHijingEventHeader.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisDataContainer.h>
#include <AliAODEvent.h>
#include <AliMultSelection.h>
#include "AliJCatalystTask.h"
#include "AliJFFlucAnalysis.h"
#include "AliJFFlucJCTask.h"

AliJFFlucJCTask::AliJFFlucJCTask(): fJCatalystTask(0), fTaskName("AliJFFlucJCTask"), fJCatalystTaskName("JCatalystTask"),
	pfa(NULL), subeventMask(SUBEVENT_A|SUBEVENT_B),
	flags(0){
	//
}

AliJFFlucJCTask::AliJFFlucJCTask(const char *name):
	AliAnalysisTaskSE(name),
	fJCatalystTask(0),
	fTaskName(name),
	fJCatalystTaskName("JCatalystTask"),
	pfa(NULL),
	subeventMask(SUBEVENT_A|SUBEVENT_B),
	fEtaMin(0.0),
	fEtaMax(0.0),
	flags(0){
	//DefineOutput(1,TList::Class());
	DefineOutput(1,TDirectory::Class());
}

// ------------------------------------------------------------------------- //
AliJFFlucJCTask::AliJFFlucJCTask(const AliJFFlucJCTask &a):
	AliAnalysisTaskSE(a.GetName()), 
	fJCatalystTask(a.fJCatalystTask),
	fTaskName(a.fTaskName),
	fJCatalystTaskName(a.fJCatalystTaskName),
	pfa(a.pfa),
	fEtaMin(a.fEtaMin),
	fEtaMax(a.fEtaMax),
	subeventMask(SUBEVENT_A|SUBEVENT_B),
	flags(0){ 
	//
}

/*AliJFFlucJCTask& AliJFFlucJCTask::operator = (const AliJFFlucJCTask& ap)
{
	this->~AliJFFlucJCTask();
	new(this) AliJFFlucJCTask(ap);
	return *this;
}*/

AliJFFlucJCTask::~AliJFFlucJCTask(){
	//
	delete pfa;
}

void AliJFFlucJCTask::UserCreateOutputObjects(){
	//
	pfa = new AliJFFlucAnalysis(fTaskName);
	pfa->SetBinning(AliJFFlucAnalysis::BINNING_CENT_PbPb); //standard 0-1, 1-2, 2-5, 5-10, ... centrality binning (not exclusive to PbPb)
	pfa->SelectSubevents(subeventMask);

	if(flags & FLUC_EBE_WEIGHTING)
		pfa->AddFlags(AliJFFlucAnalysis::FLUC_EBE_WEIGHTING);
	if(flags & FLUC_PHI_CORRECTION)
		pfa->AddFlags(AliJFFlucAnalysis::FLUC_PHI_CORRECTION);
	
	gRandom->SetSeed(151222);

	AliAnalysisManager *pman = AliAnalysisManager::GetAnalysisManager();
	fJCatalystTask = (AliJCatalystTask*)pman->GetTask(fJCatalystTaskName);
	if(!fJCatalystTask)
		AliError("AliJCatalystTask unavailable.");

	OpenFile(1);
	TDirectory *fOutput = gDirectory;
	fOutput->cd();
	pfa->UserCreateOutputObjects();
	

	PostData(1,fOutput);
}

void AliJFFlucJCTask::UserExec(Option_t*){
	const double fVertex[] = {0.0,0.0,fJCatalystTask->GetZVertex()};

	if( fJCatalystTask->GetJCatalystEntry() != fEntry ) return;

	if(fJCatalystTask->GetCentrality()>80. || fJCatalystTask->GetCentrality()<0.) return;
	pfa->Init();
	pfa->SetInputList(fJCatalystTask->GetInputList());
	pfa->SetEventCentrality(fJCatalystTask->GetCentrality());
	pfa->SetEventVertex(fVertex);
	pfa->SetEtaRange(fEtaMin, fEtaMax);

	pfa->UserExec("");
}

