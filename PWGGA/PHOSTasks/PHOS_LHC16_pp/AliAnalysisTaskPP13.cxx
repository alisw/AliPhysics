#include "iostream"

// --- Custom header files ---
#include "AliAnalysisTaskPP13.h"
#include <AliPP13AnalysisCluster.h>
#include <AliPP13TriggerProperties.h>


// --- ROOT header files ---
#include <TFile.h>
#include <TObjArray.h>
#include <TROOT.h>


// --- AliRoot header files ---
#include "AliAnalysisManager.h"
#include <AliVEvent.h>
#include <AliVCaloCells.h>
#include <AliVCluster.h>
#include <AliAODCaloCluster.h>
#include <AliVVertex.h>
#include <AliPHOSGeometry.h>
#include <AliLog.h>
#include <AliAODMCParticle.h>
#include <AliAODEvent.h>
#include "AliPHOSTriggerUtils.h"


// --- AliRoot MC headers ---
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>


ClassImp(AliAnalysisTaskPP13)

//________________________________________________________________
AliAnalysisTaskPP13::AliAnalysisTaskPP13() : AliAnalysisTaskSE(),
	fPreviousEvents(0),
	fSelections(0),
	fNMixedEvents(0)
{
	// Constructor for root I/O, do not use it
}

//________________________________________________________________
AliAnalysisTaskPP13::AliAnalysisTaskPP13(const char * name, TList * selections, Int_t nmix):
	AliAnalysisTaskSE(name),
	fPreviousEvents(0),
	fSelections(selections),
	fNMixedEvents(nmix)
{
	fSelections->SetOwner(kTRUE);

	for (int i = 0; i < fSelections->GetEntries(); ++i)
		DefineOutput(i + 1, TList::Class()); // Output starts from 1
}

//________________________________________________________________
AliAnalysisTaskPP13::~AliAnalysisTaskPP13()
{
	if (fSelections) delete fSelections;
	if (fPreviousEvents) delete fPreviousEvents;
}

//________________________________________________________________
void AliAnalysisTaskPP13::UserCreateOutputObjects()
{
	// Initialization of all outputs
	for (int i = 0; i < fSelections->GetEntries(); ++i)
	{
		AliPP13PhysicsSelection * selection = dynamic_cast<AliPP13PhysicsSelection *> (fSelections->At(i));
		selection->InitSummaryHistograms();
		PostData(i + 1, selection->GetListOfHistos()); // Output starts from 1
	}

	fPreviousEvents = new AliPP13MixingSample(fNMixedEvents);
}

//________________________________________________________________
void AliAnalysisTaskPP13::UserExec(Option_t *)
{
	// Does the job for one event

	// event
	AliVEvent * event = InputEvent();
	if (!event)
	{
		AliWarning("Can't get event");
		return;
	}

	// AliPHOSTriggerUtils * triggerUtils = new AliPHOSTriggerUtils("PHOSTrig"); 
	// triggerUtils->SetEvent(event);

	AliPP13TriggerProperties triggerProperties(
		event->GetCaloTrigger("PHOS")
	);

	// Count MB event before event cuts for every selection
	for (int i = 0; i < fSelections->GetEntries(); ++i)
	{
		AliPP13PhysicsSelection * selection = dynamic_cast<AliPP13PhysicsSelection *> (fSelections->At(i));
		selection->CountMBEvent();
	}


	// check geometry
	if (!AliPHOSGeometry::GetInstance())
	{
		AliInfo("PHOS geometry not initialized, initializing it for you");
		// Don't instantinate geometry: Use tender
		AliPHOSGeometry::GetInstance();
	}

	// Select Event
	EventFlags evtProperties;
	if (!EventSelected(event, evtProperties))
		return;

	// Set the Event handler if it's available
	AliAnalysisManager * manager = AliAnalysisManager::GetAnalysisManager();
	if (manager) 
	{
		AliInputEventHandler * inputHandler = dynamic_cast<AliInputEventHandler * >(manager->GetInputEventHandler());
		if (inputHandler)
			evtProperties.fPIDResponse = inputHandler->GetPIDResponse();
	}

	// NB: Use don't use TClonesArray as you don't want to copy the clusters
	// just use pointers
	//
	Int_t nTriggered = 0;
	TObjArray clusArray;
	clusArray.SetOwner(kTRUE);
	for (Int_t i = 0; i < event->GetNumberOfCaloClusters(); i++)
	{
		AliVCluster * clus = event->GetCaloCluster(i);

		// AliAODCaloCluster * c = dynamic_cast<AliAODCaloCluster *> (event->GetCaloCluster(i));
		if (!clus)
		{
			AliWarning("Can't get cluster");
			return;
		}

		// AliPP13AnalysisCluster * clus = new AliPP13AnalysisCluster(*c);

		if (!clus->IsPHOS())
			continue;

		// triggerProperties.FillTriggerInformation(clus);
		// Use only with triggerUtils
		// clus->SetTrigger(triggerUtils->IsFiredTrigger(c));
		// nTriggered += clus->IsTrigger();
		clusArray.Add(
			new AliAODCaloCluster(*dynamic_cast<AliAODCaloCluster *>(clus))
		);
	}
	// delete triggerUtils;

	evtProperties.fTriggerEvent = nTriggered > 0;
	evtProperties.fMcParticles = GetMCParticles(event);

	TList * pool = fPreviousEvents->GetPool(evtProperties);
	for (int i = 0; i < fSelections->GetEntries(); ++i) // Fill and Post Data to outputs
	{
		AliPP13PhysicsSelection * selection = dynamic_cast<AliPP13PhysicsSelection *> (fSelections->At(i));

		if (!selection->SelectEvent(evtProperties))
			continue;

		selection->FillHistograms(&clusArray, pool, evtProperties);
		selection->ConsiderGeneratedParticles(evtProperties);

		PostData(i + 1, selection->GetListOfHistos()); // Output starts from 1
	}
	fPreviousEvents->UpdatePool(clusArray, evtProperties);
}

//________________________________________________________________
TClonesArray * AliAnalysisTaskPP13::GetMCParticles(const AliVEvent * event) const
{
	const AliAODEvent * aodevent = dynamic_cast<const AliAODEvent*>(event);

	if (!aodevent)
		return 0;

	TClonesArray * mc = (TClonesArray*)aodevent->FindListObject(AliAODMCParticle::StdBranchName());
	return mc;
}

//________________________________________________________________
void AliAnalysisTaskPP13::Terminate(Option_t *)
{
}

//________________________________________________________________
Bool_t AliAnalysisTaskPP13::EventSelected(const AliVEvent * event, EventFlags & eprops) const
{
	// pileup
	if (event->IsPileupFromSPD(3, 0.8, 3., 2., 5.))
		return kFALSE;

	// TODO: DO we need this check?
	// cells
	AliVCaloCells * cells = event->GetPHOSCells();

	if (!cells)
	{
		AliWarning("Can't get cells");
		return kFALSE;
	}

	// primary vertex
	AliVVertex * vertex = (AliVVertex *) event->GetPrimaryVertex();
	if (!vertex)
	{
		AliWarning("Can't get primary vertex");
		return kFALSE;
	}

	vertex->GetXYZ(eprops.vtxBest);
	eprops.BC = event->GetBunchCrossNumber();

	eprops.ncontributors = vertex->GetNContributors();

	return kTRUE;
}
