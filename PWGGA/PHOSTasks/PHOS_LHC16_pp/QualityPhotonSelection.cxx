// --- Custom header files ---
#include "QualityPhotonSelection.h"
// #include "AliAnalysisTaskPP.h"

// --- AliRoot header files ---
#include <AliPHOSGeometry.h>

#include <iostream>
using namespace std;


ClassImp(QualityPhotonSelection);

//________________________________________________________________
void QualityPhotonSelection::InitSelectionHistograms()
{
	// pi0 mass spectrum
	Int_t nM       = 750;
	Double_t mMin  = 0.0;
	Double_t mMax  = 1.5;
	Int_t nPt      = 400;
	Double_t ptMin = 0;
	Double_t ptMax = 20;


	// Z-vertex
	fZvertex = new TH1F("hZvertex", "Reconstructed vertex Z-coordinate; z_{vtx}, cm; counts", 200, -12, 12);

	// Info about selected clusters
	fNcellsE = new TH2F("hNcellsE", "Cell multiplicity; E, GeV; N_{cell}", 41, 0, 40, 81, 0, 80);
	fShapeE  = new TH2F("hShapeE", "Cluster shape; E, GeV; M20, cm", 41, 0, 40, 41, 0, 40);

	fListOfHistos->Add(fZvertex);
	fListOfHistos->Add(fNcellsE);
	fListOfHistos->Add(fShapeE);

	// Test Assymetry cut
	for(Int_t i = 0; i < 2; ++i)
	{
		const char * s = (i == 0) ? "": "Mix";
		fMassPtA[i] = new TH3F(Form("h%sMassPtA", s), "(M,p_{T}, A)_{#gamma#gamma}, N_{cell}>2; M_{#gamma#gamma}, GeV; p_{T}, GeV/c", nM, mMin, mMax, nPt, ptMin, ptMax, 20, 0., 1.);
		fListOfHistos->Add(fMassPtA[i]);
	}

	// Cluster occupancy
	for(Int_t i = 0; i < 2; ++i)
	{
		TString tenergy = Form(", E %s 1,", i == 0 ? "<" : ">");
		fClusterNXZ[i] = new DetectorHistogram(new TH2F(Form("hCluNXZM_%d_", i), "Cluster N(X,Z)" + tenergy + " ; x; z", 64, 0.5, 64.5, 56, 0.5, 56.5), fListOfHistos, DetectorHistogram::kModules);
		fClusterEXZ[i] = new DetectorHistogram(new TH2F(Form("hCluEXZM_%d_", i), "Cluster E(X,Z)" + tenergy + " ; x; z", 64, 0.5, 64.5, 56, 0.5, 56.5), fListOfHistos, DetectorHistogram::kModules);
	}

	// Time maps
	fClusterTime    = new DetectorHistogram(new TH1F("hClusterTime", "Cluster Time scaled by E, ;t, s", 4800, -0.25 * 1e-6, 0.25 * 1e-6), fListOfHistos, DetectorHistogram::kModules);
	fClusterEvsT    = new DetectorHistogram(new TH2F("hClusterEvsT", "Cluster energy vs time, ; cluster energy, GeV; time, s", 100, 0., 12., 1200, -0.25 * 1e-6, 0.25 * 1e-6), fListOfHistos, DetectorHistogram::kModules);
	fClusterTimeMap = new DetectorHistogram(new TH2F("hClusterTimeMap", "Cluster time map, ; X; Z", 64, 0.5, 64.5, 56, 0.5, 56.5), fListOfHistos, DetectorHistogram::kModules);

	for(Int_t i = 0; i < 2; ++i)
	{
		const char * s = (i == 0) ? "E < 1 GeV": "E > 1 GeV";
		fClusterIdN[i] = new TH1F(Form("hClusterIdN_%d", i), Form("Cluster N(Id), %s; id ", s), 4 * 3584, 0.5, 4 * 3584 + 0.5);
		fClusterIdE[i] = new TH1F(Form("hClusterIdE_%d", i), Form("Cluster E(Id), %s; id ", s), 4 * 3584, 0.5, 4 * 3584 + 0.5);

		fListOfHistos->Add(fClusterIdN[i]);
		fListOfHistos->Add(fClusterIdE[i]);
	}
}


//________________________________________________________________
void QualityPhotonSelection::ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags)
{
	TLorentzVector p1, p2, psum;
	c1->GetMomentum(p1, eflags.vtxBest);
	c2->GetMomentum(p2, eflags.vtxBest);
	psum = p1 + p2;

	// Pair cuts can be applied here
	if (psum.M2() < 0)  return;

	// Appply asymmetry cut for pair
	Double_t asym = TMath::Abs( (p1.E() - p2.E()) / (p1.E() + p2.E()) );
	if (asym > fCuts.fAsymmetryCut) return;


	Int_t sm1, sm2, x1, z1, x2, z2;
	if ((sm1 = CheckClusterGetSM(c1, x1, z1)) < 0) return; //  To be sure that everything is Ok
	if ((sm2 = CheckClusterGetSM(c2, x2, z2)) < 0) return; //  To be sure that everything is Ok

	Double_t ma12 = psum.M();
	Double_t pt12 = psum.Pt();
	
	fMassPtA[Int_t(eflags.isMixing)]->Fill(ma12, pt12, asym);
}


//________________________________________________________________
void QualityPhotonSelection::SelectPhotonCandidates(const TObjArray * clusArray, TObjArray * candidates, const EventFlags & eflags)
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
			fClusterEvsT->FillAll(sm, sm, tof, clus->E());
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
Int_t QualityPhotonSelection::AbsId(Int_t x, Int_t z, Int_t sm) const
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
Bool_t QualityPhotonSelection::SelectEvent(const EventFlags & flgs)
{
	// Keep it this way if you decide to switch Bool_t -> Some_Other_type

	Bool_t accepted = PhotonSelection::SelectEvent(flgs);
	if (accepted)
		fZvertex->Fill(flgs.vtxBest[2]);

	return accepted;
}