#include "iostream"

// --- ROOT header files ---
#include <TFile.h>
#include <TObjArray.h>
#include <TROOT.h>

// --- AliRoot header files ---
#include "AliAnalysisTaskPP.h"
#include "AliAnalysisManager.h"
#include <AliVEvent.h>
#include <AliVCaloCells.h>
#include <AliVCluster.h>
#include <AliVVertex.h>
#include <AliPHOSGeometry.h>
#include <AliLog.h>
#include <AliAODMCParticle.h>
#include <AliAODEvent.h>


// --- AliRoot MC headers ---
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>


ClassImp(AliAnalysisTaskPP)

//________________________________________________________________
AliAnalysisTaskPP::AliAnalysisTaskPP() : AliAnalysisTaskSE(),
	fPreviousEvents(0),
	fSelections(0),
	fPHOSBadMap(),
	fNMixedEvents(0),
	fNBad(0),
	fBadCells(0)
{
	// Constructor for root I/O, do not use it
}

//________________________________________________________________
AliAnalysisTaskPP::AliAnalysisTaskPP(const char * name, TList * selections, Int_t nmix):
	AliAnalysisTaskSE(name),
	fPreviousEvents(0),
	fSelections(selections),
	fPHOSBadMap(),
	fNMixedEvents(nmix),
	fNBad(0),
	fBadCells(0)
{
	fSelections->SetOwner(kTRUE);

	for (int i = 0; i < fSelections->GetEntries(); ++i)
		DefineOutput(i + 1, TList::Class()); // Output starts from 1
}

//________________________________________________________________
AliAnalysisTaskPP::~AliAnalysisTaskPP()
{
	if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) delete fSelections;
	if (fBadCells) delete [] fBadCells;
	if (fPreviousEvents) delete fPreviousEvents;
}

//________________________________________________________________
void AliAnalysisTaskPP::UserCreateOutputObjects()
{
	// Initialization of all outputs
	for (int i = 0; i < fSelections->GetEntries(); ++i)
	{
		PhotonSelection * selection = dynamic_cast<PhotonSelection *> (fSelections->At(i));
		selection->InitSummaryHistograms();
		PostData(i + 1, selection->GetListOfHistos()); // Output starts from 1
	}

	fPreviousEvents = new MixingSample(fNMixedEvents);
}

//________________________________________________________________
void AliAnalysisTaskPP::UserExec(Option_t *)
{
	// Does the job for one event

	// event
	AliVEvent * event = InputEvent();
	if (!event)
	{
		AliWarning("Can't get event");
		return;
	}

	// Count MB event before event cuts for every selection 
	for (int i = 0; i < fSelections->GetEntries(); ++i) 
	{
		PhotonSelection * selection = dynamic_cast<PhotonSelection *> (fSelections->At(i));
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


	// NB: Use don't use TClonesArray as you don't want to copy the clusters
	// just use pointers
	//
	TObjArray clusArray;
	for (Int_t i = 0; i < event->GetNumberOfCaloClusters(); i++)
	{
		AliVCluster * clus = event->GetCaloCluster(i);
		if (!clus)
		{
			AliWarning("Can't get cluster");
			return;
		}

		// only basic filtering
		if (!clus->IsPHOS()) continue;
		if (IsClusterBad(clus)) continue;

		clusArray.Add(clus);
	}

	evtProperties.fMcParticles = GetMCParticles(event);

	TList * pool = fPreviousEvents->GetPool(evtProperties);
	for (int i = 0; i < fSelections->GetEntries(); ++i) // Fill and Post Data to outputs
	{
		PhotonSelection * selection = dynamic_cast<PhotonSelection *> (fSelections->At(i));

		if (!selection->SelectEvent(evtProperties))
			continue;

		selection->FillPi0Mass(&clusArray, pool, evtProperties);
		selection->ConsiderGeneratedParticles(evtProperties);

		PostData(i + 1, selection->GetListOfHistos()); // Output starts from 1
	}
	fPreviousEvents->UpdatePool(clusArray, evtProperties);
}

//________________________________________________________________
TClonesArray * AliAnalysisTaskPP::GetMCParticles(const AliVEvent * event) const
{
	// TODO: Handle the ESD case here
	const AliAODEvent * aodevent = dynamic_cast<const AliAODEvent*>(event);

	if (!aodevent)
		return 0;

	TClonesArray * mc = (TClonesArray*)aodevent->FindListObject(AliAODMCParticle::StdBranchName());
	return mc;
}

//________________________________________________________________
void AliAnalysisTaskPP::Terminate(Option_t *)
{
}

//________________________________________________________________
Bool_t AliAnalysisTaskPP::EventSelected(const AliVEvent * event, EventFlags & eprops) const
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

//____________________________________________________________
void AliAnalysisTaskPP::SetBadCells(Int_t badcells[], Int_t nbad)
{
	// Set absId numbers for bad cells;
	// clusters which contain a bad cell will be rejected.

	if (fBadCells) delete [] fBadCells;

	// switch off bad cells, if asked
	if (nbad <= 0)
	{
		fNBad = 0;
		return;
	}

	fNBad = nbad;
	fBadCells = new Int_t[nbad];

	for (Int_t i = 0; i < nbad; i++)
		fBadCells[i] = badcells[i];
}

//________________________________________________________________
Bool_t AliAnalysisTaskPP::CellInPhos(Int_t absId, Int_t & sm, Int_t & ix, Int_t & iz) const
{
	// Converts cell absId --> (sm,ix,iz);
	AliPHOSGeometry * geomPHOS = AliPHOSGeometry::GetInstance();
	// AliPHOSGeometry * geomPHOS = AliPHOSGeometry::GetInstance("IHEP");
	if (!geomPHOS)
	{
		AliWarning("Something is wrong with PHOS Geometry. Check if you initialize it in UserExec!");
		return kTRUE;
	}

	Int_t relid[4];
	geomPHOS->AbsToRelNumbering(absId, relid);
	sm = relid[0];
	ix = relid[2];
	iz = relid[3];
	return (sm >= kMinModule) && (sm <= kMaxModule);
}

//________________________________________________________________
Bool_t AliAnalysisTaskPP::IsClusterBad(AliVCluster * clus) const
{
	// Returns true if cluster contains a bad cell
	for (Int_t b = 0; b < fNBad; b++)
		for (Int_t c = 0; c < clus->GetNCells(); c++)
			if (clus->GetCellAbsId(c) == fBadCells[b])
				return kTRUE;

	// If fBadCells array is empty then use BadMap

	if ( !fPHOSBadMap[0] ) return kFALSE;

	Int_t sm, ix, iz;
	for (Int_t c = 0; c < clus->GetNCells(); c++) // Loop over all cells in cluster
	{
		if (!CellInPhos(clus->GetCellAbsId(c), sm, ix, iz)) // Reject cells outside PHOS
			return kTRUE;

		if (!fPHOSBadMap[sm - kMinModule]) // Warn if something is wrong
			AliError(Form("No Bad map for PHOS module %d", sm));

		if (fPHOSBadMap[sm - kMinModule]->GetBinContent(ix, iz) > 0) // Check if cell is bad
			return kTRUE;

	}

	return kFALSE;
}

//________________________________________________________________
void AliAnalysisTaskPP::SetBadMap(const char * filename)
{
	TFile * fBadMap = TFile::Open(filename);
	if (!fBadMap->IsOpen())
		AliFatal(Form("Cannot set BadMap %s doesn't exist", filename));

	std::cout << "\n\n...Adding PHOS bad channel map \n"  << std::endl;
	gROOT->cd();

	for (Int_t module = kMinModule; module <= kMaxModule; ++module)
	{
		TH2I * h = (TH2I *) fBadMap->Get(Form("PHOS_BadMap_mod%d", module));
		if (!h) AliFatal( Form("PHOS_BadMap_mod%d doesn't exist", module));

		fPHOSBadMap[module - kMinModule] = new TH2I(*h);
		std::cout << "Set " <<  fPHOSBadMap[module - kMinModule]->GetName() << std::endl;
	}
	fBadMap->Close();
	std::cout << "\n\n...PHOS BadMap is set now." << std::endl;
}
