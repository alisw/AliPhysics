// --- Custom header files ---
#include "AliPP13PhysicsSelection.h"

// --- ROOT system ---
#include <TH2F.h>
#include <TH3F.h>


// --- AliRoot header files ---
#include <AliLog.h>

#include <iostream>
using namespace std;



ClassImp(AliPP13PhysicsSelection);

//________________________________________________________________
void AliPP13PhysicsSelection::FillHistograms(TObjArray * clusArray, TList * pool, const EventFlags & eflags)
{
	// Ensure that we are not doing mixing
	EventFlags flags = eflags;
	flags.isMixing = kFALSE;

	// Select photons
	TObjArray photonCandidates;

	// Select particles and fill single-particle distributions
	SelectPhotonCandidates(clusArray, &photonCandidates, flags);

	// Fill two-particle distributions (e.g. invariant masses)
	SelectTwoParticleCombinations(photonCandidates, flags);

	// Fill all mixing combinations
	MixPhotons(photonCandidates, pool, flags);
}

//________________________________________________________________
void AliPP13PhysicsSelection::SelectTwoParticleCombinations(const TObjArray & photonCandidates, const EventFlags & eflags)
{
	// All possible combinations on photon candadates
	// Int_t counter = 0;
	for (Int_t i = 0; i < photonCandidates.GetEntriesFast(); i++)
	{
		AliVCluster * clus1 = (AliVCluster *) photonCandidates.At(i);

		// second cluster loop
		for (Int_t j = i + 1; j < photonCandidates.GetEntriesFast(); j++)
		{
			// ++counter;
			AliVCluster * clus2 = (AliVCluster *) photonCandidates.At(j);

			if(!fCuts.AcceptPair(clus1, clus2, eflags))
				continue;

			ConsiderPair(clus1, clus2, eflags);
		} // second cluster loop
	} // cluster loop

	// Int_t Nn = photonCandidates.GetEntriesFast();
	// std::cout << "Number of combinations: " << counter << " should be " << Nn * (Nn - 1.) / 2. << std::endl;
}

//________________________________________________________________
void AliPP13PhysicsSelection::MixPhotons(TObjArray & photonCandidates, TList * pool, const EventFlags & eflags)
{
	// Notify all selections that this is mixing
	EventFlags mflags = eflags;
	mflags.isMixing = kTRUE;

	// Apply user selection
	TObjArray previousPhotons;
	for (Int_t ev = 0; ev < pool->GetEntries(); ++ev)
	{
		TObjArray * previousClusters = dynamic_cast<TObjArray *>(pool->At(ev));
		if (!previousClusters) continue;
		SelectPhotonCandidates(previousClusters, &previousPhotons, mflags);
	}

	// old cluster loop
	for (Int_t j = 0; j < previousPhotons.GetEntriesFast(); ++j)
	{
		AliVCluster * clus1 = (AliVCluster *) previousPhotons.At(j);

		// Check all possible combinations
		for (Int_t i = 0; i < photonCandidates.GetEntriesFast(); ++i)
		{
			AliVCluster * clus2 = (AliVCluster *) photonCandidates.At(i);

			if(!fCuts.AcceptPair(clus1, clus2, eflags))
				continue;

			ConsiderPair(clus1, clus2, mflags);
		}
	} // old cluster loop
}

//________________________________________________________________
Int_t AliPP13PhysicsSelection::CheckClusterGetSM(const AliVCluster * clus, Int_t & x, Int_t & z) const
{
	// Apply common cluster cuts and return supermodule number on success.
	// Return -1 if cuts not passed or an error occured.

	if (!clus->IsPHOS()) return -1;
	if (clus->GetType() != AliVCluster::kPHOSNeutral) return -1; // don't use CPV
	if (clus->GetNCells() < 1) return -1;

	Float_t  position[3];
	clus->GetPosition(position);
	TVector3 global(position);

	Int_t relId[4];
	AliPHOSGeometry * phosGeometry = AliPHOSGeometry::GetInstance();
	phosGeometry->GlobalPos2RelId(global, relId) ;

	Int_t sm = relId[0];
	x = relId[2];
	z = relId[3];

	// check for data corruption to avoid segfaults
	if (sm < 1 || sm > 5)
	{
		AliError("Data is corrupted!!!");
		return -1;
	}

	return sm;
}

//________________________________________________________________
TLorentzVector AliPP13PhysicsSelection::ClusterMomentum(const AliVCluster * c1, const EventFlags & eflags) const
{
	TLorentzVector p;
	c1->GetMomentum(p, eflags.vtxBest);
	return p;
}

//________________________________________________________________
void AliPP13PhysicsSelection::InitSummaryHistograms()
{
	// Find better place to apply this
	fListOfHistos = new TList();
	fListOfHistos->SetOwner(kTRUE);
	InitSelectionHistograms();

	TString cuts = Form(
	                   ";Cuts: |Z_{vtx}| < 10 cm, no pileup spd, E_{min}^{clu} = %.2g GeV, A =  %.2g, N_{min}^{cell} = %d, t_{clus} = %0.3g ns, dist = %.2g",
	                   fCuts.fClusterMinE,
	                   fCuts.fAsymmetryCut,
	                   fCuts.fNCellsCut,
	                   fCuts.fTimingCut * 1e+9,
	                   fCuts.fMinimalDistance
	               );

	cout << "Adding " << this->GetName() << ": " << this->GetTitle() << endl;
	this->SetTitle(this->GetTitle() + cuts);
	fWeights->Report(fListOfHistos);
	fListOfHistos->AddFirst(
	    new TH1C("selection_description", this->GetTitle(), 1, 0, 1)
	); // Very important!!! Description, dummy way

	// The true event counter
	fEventCounter = new TH1F("EventCounter", "Event cuts", 5, 0, 5);
	fEventCounter->GetXaxis()->SetBinLabel(1, "MB");
	fEventCounter->GetXaxis()->SetBinLabel(2, "all good");
	fEventCounter->GetXaxis()->SetBinLabel(3, "|Z_{vtx}| < 10");
	fEventCounter->GetXaxis()->SetBinLabel(4, Form("N_{vtx contrib} > %d", fCuts.fNContributors - 1));
	fEventCounter->GetXaxis()->SetBinLabel(5, "N_{#gamma} > 2");
	fListOfHistos->AddFirst(fEventCounter);

}

//________________________________________________________________
void AliPP13PhysicsSelection::CountMBEvent()
{
	fEventCounter->Fill(EventFlags::kMB);
}


//________________________________________________________________
AliPP13PhysicsSelection::~AliPP13PhysicsSelection()
{
	// if (fWeights)
	// delete fWeights;

	// Don't delete fEventCounter and other ROOT objects
	// root has it's own memory management.
	//
	delete fListOfHistos;
}

//________________________________________________________________
Bool_t AliPP13PhysicsSelection::SelectEvent(const EventFlags & flgs)
{
	// All events
	fEventCounter->Fill(EventFlags::kGood);


	if (TMath::Abs(flgs.vtxBest[2]) > 10)
		return kFALSE;

	// Z-vtx events
	fEventCounter->Fill(EventFlags::kZvertex);


	// Number of contributors > 0
	if (flgs.ncontributors < fCuts.fNContributors)
		return kFALSE;

	fEventCounter->Fill(EventFlags::kNcontributors);

	// Physical Events
	return kTRUE;
}

//________________________________________________________________
void AliPP13PhysicsSelection::SelectPhotonCandidates(const TObjArray * clusArray, TObjArray * candidates, const EventFlags & eflags)
{
	// Don't return TObjArray: force user to handle candidates lifetime
	for (Int_t i = 0; i < clusArray->GetEntriesFast(); i++)
	{
		AliVCluster * clus = (AliVCluster *) clusArray->At(i);

		if(!AcceptPhotonCandidate(clus, eflags))
			continue;

		candidates->Add(clus);

		// Fill histograms only for real events
		if (eflags.isMixing)
			continue;

		FillClusterHistograms(clus, eflags);
	}

	if (candidates->GetEntriesFast() > 1 && !eflags.isMixing)
		fEventCounter->Fill(EventFlags::kTwoPhotons);
}

//________________________________________________________________
Bool_t AliPP13PhysicsSelection::AcceptPhotonCandidate(AliVCluster * clus, const EventFlags & eflags)
{
	(void) eflags;

	Int_t sm, x, z;
	if ((sm = CheckClusterGetSM(clus, x, z)) < 0)
		return kFALSE;

	if (!fCuts.AcceptCluster(clus))
		return kFALSE;

	return kTRUE;
}
