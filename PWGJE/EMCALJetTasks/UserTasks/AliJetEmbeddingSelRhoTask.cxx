#include <AliRhoParameter.h>
#include <AliLog.h>

#include "AliJetEmbeddingSelRhoTask.h"

ClassImp(AliJetEmbeddingSelRhoTask)

//____________________________________________________________________________

AliJetEmbeddingSelRhoTask::AliJetEmbeddingSelRhoTask() : AliJetEmbeddingTask(),
fRhoMin(0),
fRhoMax(1000),
fRhoName(""),
fhQARhoEventRejection(),
fhQARho()
{
	// Default constructor
}

//____________________________________________________________________________

AliJetEmbeddingSelRhoTask::AliJetEmbeddingSelRhoTask(const char *name) : AliJetEmbeddingTask(name),
fRhoMin(0),
fRhoMax(1000),
fRhoName(""),
fhQARhoEventRejection(),
fhQARho()
{
	// Standard constructor
}

//____________________________________________________________________________

void AliJetEmbeddingSelRhoTask::UserCreateOutputObjects(){
	
	AliJetEmbeddingTask::UserCreateOutputObjects();
	
	Int_t    nBinsRho = 50;
	Double_t minRho = 0.;
	Double_t maxRho = 20.;

	fhQARhoEventRejection = new TH1F("fhQARhoEventRejection", "Event selection: #rho", 2, 0, 2.);
	fOutput->Add(fhQARhoEventRejection);
	
	fhQARho = new TH1F("fhQARho", "#rho of selected events; #rho;", nBinsRho, minRho, maxRho);
	fOutput->Add(fhQARho);
	
}
//____________________________________________________________________________

void AliJetEmbeddingSelRhoTask::Run(){
	
	/// Select the event with rho in the desired range
	AliRhoParameter* rhomParam = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(fRhoName));
	if(!rhomParam){
		AliFatal(Form("Cannot find %s in the event", fRhoName.Data()));
	}
	Double_t rhoevent = rhomParam->GetVal();

	Bool_t outofrhorange = kFALSE;

	if(rhoevent < fRhoMin ||  rhoevent >= fRhoMax) outofrhorange = kTRUE;

	if(!outofrhorange){
		fhQARhoEventRejection->Fill(0);
		fhQARho->Fill(rhoevent);
	}
	else {
		fhQARhoEventRejection->Fill(1);
		//don't run the embedding in this event
		return;
	}

	/// now continue with the embedding
	AliJetEmbeddingTask::Run();
}
//____________________________________________________________________________

void AliJetEmbeddingSelRhoTask::Terminate(Option_t *option){
   Printf("Task terminated");
   AliJetEmbeddingTask::Terminate(option);
}

