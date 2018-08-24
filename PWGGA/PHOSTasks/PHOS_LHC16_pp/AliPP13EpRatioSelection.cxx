// --- Custom header files ---
#include "AliPP13EpRatioSelection.h"

// --- ROOT system ---

// --- AliRoot header files ---
#include <AliPP13EpRatioSelection.h>
#include <AliAODTrack.h>


ClassImp(AliPP13EpRatioSelection);

//________________________________________________________________
void AliPP13EpRatioSelection::FillHistograms(TObjArray * clusArray, TList * pool, const EventFlags & eflags)
{
	// NB: This selection doesn't fill 2-particle distributions
	(void) pool;
	
	// Select photons
	TObjArray photonCandidates;
	SelectPhotonCandidates(clusArray, &photonCandidates, eflags);
}

//________________________________________________________________
void AliPP13EpRatioSelection::InitSelectionHistograms()
{
	// pi0 mass spectrum
	Int_t nM       = 50;
	Double_t mMin  = 0.0;
	Double_t mMax  = 2.;
	Int_t nPt      = 400;
	Double_t ptMin = 0;
	Double_t ptMax = 20;

	for (Int_t i = 0; i < 2; ++i)
	{
		const char * species = (i == 0) ? "Electrons" : "Non-Electron";

		TH2 * patternP = new TH2F(
		    Form("hEp%sE", species),
		    "E/p ratio vs. E^{cluster}, ; E/p; E^{track} ,GeV",
		    nM, mMin, mMax, nPt, ptMin, ptMax
		);

		fEpE[i] = new AliPP13DetectorHistogram(
		    patternP,
		    fListOfHistos,
		    AliPP13DetectorHistogram::kModules
		);

		TH2 * patternPt = new TH2F(
		    Form("hEp%sP", species),
		    "E/p ratio vs. p_{T}^{cluster}, ; E/p; p_{T}^{track} ,GeV/c",
		    nM, mMin, mMax, nPt, ptMin, ptMax
		);

		fEpPt[i] = new AliPP13DetectorHistogram(
		    patternPt,
		    fListOfHistos,
		    AliPP13DetectorHistogram::kModules
		);
	}

	fTPCSignal[0] = new TH2F("hEpRatioNSigmaElectron", "E/p ratio vs. N_{#sigma}^{electron}; E/p; n#sigma^{electron}", nM, mMin, mMax, 40, -5, 5);
	fTPCSignal[1] = new TH2F("hTPCSignal_Electron", "TPC dE/dx vs. electron momentum; p^{track}, GeV/c; dE/dx, a.u.", 40, 0, 20, 200, 0, 200);
	fTPCSignal[2] = new TH2F("hTPCSignal_Non-Electron", "TPC dE/dx vs. non-electron momentum; p^{track}, GeV/c; dE/dx, a.u.", 40, 0, 20, 200, 0, 200);
	fTPCSignal[3] = new TH2F("hTPCSignal_All", "TPC dE/dx vs. all particles momentum; p^{track}, GeV/c; dE/dx, a.u.", 40, 0, 20, 200, 0, 200);

	fPosition[0] = new AliPP13DetectorHistogram(
	    new TH3F("hdXvsXvsPt_plus", "dX vs. X positive, ; X (cm); dX (cm); p_{T}^{track +} (GeV/c)", 160, -80, 80, 80, -20, 20, 40, 0, 20),
	    fListOfHistos,
	    AliPP13DetectorHistogram::kModules
	);

	fPosition[1] = new AliPP13DetectorHistogram(
	    new TH3F("hdXvsXvsPt_minus", "dX vs. X negative, ; X (cm); dX (cm); p_{T}^{track -} (GeV/c)", 160, -80, 80, 80, -20, 20, 40, 0, 20),
	    fListOfHistos,
	    AliPP13DetectorHistogram::kModules
	);

	fPosition[2] = new AliPP13DetectorHistogram(
	    new TH3F("hdZvsZvsPt", "dZ vs. Z; Z (cm); dZ (cm), ; p_{T}^{track} (GeV/c)", 160, -80, 80, 80, -20, 20, 40, 0, 20),
	    fListOfHistos,
	    AliPP13DetectorHistogram::kModules
	);


	fPosition[3] = new AliPP13DetectorHistogram(
	    new TH3F("hdZvsZvsPtElectron", "dZ vs. Z of e^{#pm}, ; Z (cm); dZ (cm); p_{T}^{track} (GeV/c)", 160, -80, 80, 80, -20, 20, 40, 0, 20),
	    fListOfHistos,
	    AliPP13DetectorHistogram::kModules
	);//for radial displacement


	for (int i = 0; i < 4; ++i)
		fListOfHistos->Add(fTPCSignal[i]);

	for (Int_t i = 0; i < fListOfHistos->GetEntries(); ++i)
	{
		TH1 * hist = dynamic_cast<TH1 *>(fListOfHistos->At(i));
		if (!hist) continue;
		hist->Sumw2();
	}
}

//________________________________________________________________
void AliPP13EpRatioSelection::FillClusterHistograms(const AliVCluster * cluster, const EventFlags & eflags)
{
	Float_t nsigma_min = -1.5;
	Float_t nsigma_max = 3;

	// Don't do anything if pidresponse wasn't defined
	if (!eflags.fPIDResponse)
		return;

	// No tracks
	if ( !(cluster->GetNTracksMatched() > 0) )
		return;

	AliVTrack * track = dynamic_cast<AliVTrack*>(cluster->GetTrackMatched(0));

	// The track wasn't found
	if (!track)
		return;

	Bool_t isHybridTrack = dynamic_cast<AliAODTrack*>(track)->IsHybridGlobalConstrainedGlobal();//hybrid track

	// Take only hybrid tracks
	if (!isHybridTrack)
		return;

	Int_t sm, x1, z1;
	if ((sm = CheckClusterGetSM(cluster, x1, z1)) < 0)
		return; //  To be sure that everything is Ok	

	Double_t trackPt = track->Pt();

	Int_t charge = track->Charge();
	Double_t dx = cluster->GetTrackDx();
	Double_t dz = cluster->GetTrackDz();
	TVector3 local = LocalPosition(cluster);

	fPosition[Int_t(charge > 0)]->FillAll(sm, sm, local.X(), dx, trackPt);
	fPosition[1]->FillAll(sm, sm, local.Z(), dz, trackPt);

	Double_t energy = cluster->E();
	Double_t trackP = track->P();
	Double_t dEdx = track->GetTPCsignal();
	Double_t EpRatio = energy / trackP;

	fTPCSignal[3]->Fill(trackP, dEdx);

	// TODO: Accept electron cuts
	// TODO: Ensure not neutra particles ?
	//

	Double_t nSigma = eflags.fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);

	fTPCSignal[0]->Fill(EpRatio, nSigma);
	Bool_t isElectron = (nsigma_min < nSigma && nSigma < nsigma_max) ;



	if (isElectron)
	{
		fEpE[0]->FillAll(sm, sm, EpRatio, energy);
		fEpPt[0]->FillAll(sm, sm, EpRatio, trackPt);
		fTPCSignal[1]->Fill(trackP, dEdx);
        if(0.8 < energy / trackP && energy / trackP < 1.2) 
			fPosition[3]->FillAll(sm, sm, local.Z(), dz, trackPt);
	}
	if (!isElectron)
	{
		fEpE[1]->FillAll(sm, sm, EpRatio, energy);
		fEpPt[1]->FillAll(sm, sm, EpRatio, trackPt);
		fTPCSignal[2]->Fill(trackP, dEdx);
	}
}

//________________________________________________________________
TVector3 AliPP13EpRatioSelection::LocalPosition(const AliVCluster * cluster) const
{
	// Fill the global cluster position
	Float_t  position[3];
	cluster->GetPosition(position);
	TVector3 global(position);

	// Find the module
	Int_t x, z;
	Int_t sm = this->CheckClusterGetSM(cluster, x, z);

	TVector3 local;

	// Allow errors, don't check this pointer
	AliPHOSGeometry * phosGeometry = AliPHOSGeometry::GetInstance();
	phosGeometry->Global2Local(local, global, sm);

	return local;
}
