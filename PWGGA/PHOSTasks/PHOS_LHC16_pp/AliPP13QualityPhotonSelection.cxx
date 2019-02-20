// --- Custom header files ---
#include "AliPP13QualityPhotonSelection.h"

// --- AliRoot header files ---
#include <AliPHOSGeometry.h>

#include <iostream>
using namespace std;


ClassImp(AliPP13QualityPhotonSelection);

//________________________________________________________________
void AliPP13QualityPhotonSelection::InitSelectionHistograms()
{
	// Z-vertex
	fZvertex = new TH1F("hZvertex", "Reconstructed vertex Z-coordinate; z_{vtx}, cm; counts", 200, -12, 12);

	// Info about selected clusters
	fNcellsE = new TH2F("hNcellsE", "Cell multiplicity; E, GeV; N_{cell}", 41, 0, 40, 81, 0, 80);
	fShapeE  = new TH2F("hShapeE", "Cluster shape; E, GeV; M20, cm", 41, 0, 40, 41, 0, 40);

	fListOfHistos->Add(fZvertex);
	fListOfHistos->Add(fNcellsE);
	fListOfHistos->Add(fShapeE);

	// Cluster occupancy
	for(Int_t i = 0; i < 2; ++i)
	{
		TString tenergy = Form(", E %s 1,", i == 0 ? "<" : ">");
		fClusterNXZ[i] = new AliPP13DetectorHistogram(new TH2F(Form("hCluNXZM_%d_", i), "Cluster N(X,Z)" + tenergy + " ; x; z", 64, 0.5, 64.5, 56, 0.5, 56.5), fListOfHistos, AliPP13DetectorHistogram::kModules);
		fClusterEXZ[i] = new AliPP13DetectorHistogram(new TH2F(Form("hCluEXZM_%d_", i), "Cluster E(X,Z)" + tenergy + " ; x; z", 64, 0.5, 64.5, 56, 0.5, 56.5), fListOfHistos, AliPP13DetectorHistogram::kModules);
	}

	// Time maps
	fClusterTime    = new AliPP13DetectorHistogram(new TH1F("hClusterTime", "Cluster Time scaled by E, ;t, s", 4800, -0.25 * 1e-6, 0.25 * 1e-6), fListOfHistos, AliPP13DetectorHistogram::kModules);
	fClusterEvsT    = new AliPP13DetectorHistogram(new TH2F("hClusterEvsT", "Cluster energy vs time, ; cluster energy, GeV; time, s", 100, 0., 12., 1200, -0.25 * 1e-6, 0.25 * 1e-6), fListOfHistos, AliPP13DetectorHistogram::kModules);
	fClusterTimeWide = new AliPP13DetectorHistogram(new TH1F("hClusterTimeWide", "Cluster Time scaled by E, ;t, s", 4800, -0.25 * 1e-3, 0.25 * 1e-3), fListOfHistos, AliPP13DetectorHistogram::kModules);
	fClusterEvsTWide = new AliPP13DetectorHistogram(new TH2F("hClusterEvsTWide", "Cluster energy vs time, ; cluster energy, GeV; time, s", 100, 0., 12., 1200, -0.25 * 1e-2, 0.25 * 1e-2), fListOfHistos, AliPP13DetectorHistogram::kModules);
	fClusterTimeMap = new AliPP13DetectorHistogram(new TH2F("hClusterTimeMap", "Cluster time map, ; X; Z", 64, 0.5, 64.5, 56, 0.5, 56.5), fListOfHistos, AliPP13DetectorHistogram::kModules);
	fAsymmetry      = new TH1F("hAsymmetry", "Asymmetry between clusters; asymmetry A = (E_{1} - E_{2})/(E_{1} + E_{2})", 500, -0.5, 1.5);

	for(Int_t i = 0; i < 2; ++i)
	{
		const char * s = (i == 0) ? "E < 1 GeV": "E > 1 GeV";
		fClusterIdN[i] = new TH1F(Form("hClusterIdN_%d", i), Form("Cluster N(Id), %s; id ", s), 4 * 3584, 0.5, 4 * 3584 + 0.5);
		fClusterIdE[i] = new TH1F(Form("hClusterIdE_%d", i), Form("Cluster E(Id), %s; id ", s), 4 * 3584, 0.5, 4 * 3584 + 0.5);

		fListOfHistos->Add(fClusterIdN[i]);
		fListOfHistos->Add(fClusterIdE[i]);
	}
	fListOfHistos->Add(fAsymmetry);
}


//________________________________________________________________
void AliPP13QualityPhotonSelection::ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags)
{
	TLorentzVector p1, p2, psum;
	c1->GetMomentum(p1, eflags.vtxBest);
	c2->GetMomentum(p2, eflags.vtxBest);
	psum = p1 + p2;

	// Pair cuts can be applied here
	if (psum.M2() < 0)  return;

	Int_t sm1, sm2, x1, z1, x2, z2;
	if ((sm1 = CheckClusterGetSM(c1, x1, z1)) < 0) return; //  To be sure that everything is Ok
	if ((sm2 = CheckClusterGetSM(c2, x2, z2)) < 0) return; //  To be sure that everything is Ok

	Double_t asym = TMath::Abs( (p1.E() - p2.E()) / (p1.E() + p2.E()) );

	if(!eflags.isMixing)
		fAsymmetry->Fill(asym);
}


//________________________________________________________________
void AliPP13QualityPhotonSelection::SelectPhotonCandidates(const TObjArray * clusArray, TObjArray * candidates, const EventFlags & eflags)
{
	// Don't return TObjArray: force user to handle candidates lifetime
	Int_t sm, x, z;
	for (Int_t i = 0; i < clusArray->GetEntriesFast(); i++)
	{

		// cout << "Start" << endl;
		AliVCluster * clus = (AliVCluster *) clusArray->At(i);
		TLorentzVector p;
		clus->GetMomentum(p, eflags.vtxBest);

		if ((sm = CheckClusterGetSM(clus, x, z)) < 0) 
			continue;

		if (clus->E() < fCuts.fClusterMinE) continue;

		Float_t tof = clus->GetTOF();

		if (!eflags.isMixing)
		{
			fClusterTime->FillAll(sm, sm, clus->E());
			fClusterTimeWide->FillAll(sm, sm, clus->E());
			fClusterEvsT->FillAll(sm, sm, tof, clus->E());
			fClusterEvsTWide->FillAll(sm, sm, tof, clus->E());
			fClusterTimeMap->FillAll(sm, sm, x, z, tof);
		}
		if (TMath::Abs(clus->GetTOF()) > fCuts.fTimingCut) continue;

		if (!eflags.isMixing) fNcellsE->Fill(p.E(), clus->GetNCells());
		if (!eflags.isMixing) fShapeE->Fill(clus->GetM20(), clus->GetNCells());

		if (clus->GetNCells() < fCuts.fNCellsCut) continue;
		candidates->Add(clus);


		// Fill histograms only for real events
		if (eflags.isMixing)
			continue;

		Float_t energy = clus->E();
		Int_t isHighECluster = Int_t(energy > 1.);

		fClusterNXZ[isHighECluster]->FillAll(sm, sm, x, z, 1.);
		fClusterEXZ[isHighECluster]->FillAll(sm, sm, x, z, energy);

		// cout << "Stop" << endl;
		// Use sm number as weight in order to distinguish different modules
		//
		Int_t absId = AbsId(x, z, sm);
		fClusterIdN[isHighECluster]->Fill(absId, sm);
		fClusterIdE[isHighECluster]->Fill(absId, sm * energy);
	}

	if (candidates->GetEntriesFast() > 1 && !eflags.isMixing)
		fEventCounter->Fill(EventFlags::kTwoPhotons);
}

// Default version defined in PHOSutils uses ideal geometry
// use this instead
Int_t AliPP13QualityPhotonSelection::AbsId(Int_t x, Int_t z, Int_t sm) const
{
	// Converts cell absId --> (sm,eta,phi);
	// AliPHOSGeometry * geomPHOS = AliPHOSGeometry::GetInstance("Run2");
	// AliPHOSGeometry * geomPHOS = AliPHOSGeometry::GetInstance("IHEP");
	AliPHOSGeometry * geomPHOS = AliPHOSGeometry::GetInstance();

	if (!geomPHOS)
		AliFatal("Geometry is not defined");

	for (Int_t id = (3584 * (sm - 1) + 1); id <= (3584 * (sm)); ++id)
	{
		Int_t rel[4];
		geomPHOS->AbsToRelNumbering(id, rel);

		if (rel[2] == x && rel[3] == z)
			return id;
	}

	// There is no such a cell
	return -1;
}

//________________________________________________________________
Bool_t AliPP13QualityPhotonSelection::SelectEvent(const EventFlags & flgs)
{
	// Keep it this way if you decide to switch Bool_t -> Some_Other_type

	Bool_t accepted = AliPP13PhysicsSelection::SelectEvent(flgs);
	if (accepted)
		fZvertex->Fill(flgs.vtxBest[2]);

	return accepted;
}
