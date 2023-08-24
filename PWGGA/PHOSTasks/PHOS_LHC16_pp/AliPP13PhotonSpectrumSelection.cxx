// --- Custom header files ---
#include "AliPP13PhotonSpectrumSelection.h"

// --- ROOT system ---
#include <TH2F.h>
#include <TH3F.h>

// --- AliRoot header files ---
#include <AliPHOSGeometry.h>

#include <iostream>
using namespace std;

Bool_t TestLambda(Double_t l1, Double_t l2, Double_t R)
{
	// TODO: Check these parameters
	Double_t l1Mean = 1.22 ;
	Double_t l2Mean = 2.0 ;
	Double_t l1Sigma = 0.42 ;
	Double_t l2Sigma = 0.71 ;
	Double_t c = -0.59 ;
	Double_t R2 = (l1 - l1Mean) * (l1 - l1Mean) / l1Sigma / l1Sigma + (l2 - l2Mean) * (l2 - l2Mean) / l2Sigma / l2Sigma - c * (l1 - l1Mean) * (l2 - l2Mean) / l1Sigma / l2Sigma;
	return (R2 < R * R) ;
}


ClassImp(AliPP13PhotonSpectrumSelection);

//________________________________________________________________
void AliPP13PhotonSpectrumSelection::InitSelectionHistograms()
{

	Int_t nPt      = 400;
	Double_t ptMin = 0;
	Double_t ptMax = 40;

	this->SetTitle(Form("%s ## CPV = %.1f cm, Disp = %1.f cm", this->GetTitle(), fDistanceCPV, fDispersionCut));

	fSpectrum     = new AliPP13DetectorHistogram(new TH1F("hClusterPt_",               "Cluster p_{T}, ; cluster p_{T} (GeV/#it{c}); counts", nPt, ptMin, ptMax), fListOfHistos, AliPP13DetectorHistogram::kModules);
	fSpectrumCPV  = new AliPP13DetectorHistogram(new TH1F("hClusterPt_cpv_",      Form("Cluster p_{T} with CPV cut %.1f cm, ; cluster p_{T} (GeV/#it{c}); counts", fDistanceCPV) , nPt, ptMin, ptMax), fListOfHistos, AliPP13DetectorHistogram::kModules);
	fSpectrumDisp = new AliPP13DetectorHistogram(new TH1F("hClusterPt_disp_",     Form("Cluster p_{T} with dispersion cut %.1f cm, ; cluster p_{T} (GeV/#it{c}); counts", fDispersionCut), nPt, ptMin, ptMax), fListOfHistos, AliPP13DetectorHistogram::kModules);
	fSpectrumBoth = new AliPP13DetectorHistogram(new TH1F("hClusterPt_cpv_disp_", Form("Cluster p_{T} with CPV %.1f cm and dispersion %.1f cm cuts, ; cluster p_{T} (GeV/#it{c}); counts", fDistanceCPV, fDispersionCut), nPt, ptMin, ptMax), fListOfHistos, AliPP13DetectorHistogram::kModules);


	for(Int_t i = 0; i < fListOfHistos->GetEntries(); ++i)
	{
		TH1 * hist = dynamic_cast<TH1 *>(fListOfHistos->At(i));
		if(!hist) continue;
		hist->Sumw2();
	}	
}

//________________________________________________________________
void AliPP13PhotonSpectrumSelection::FillClusterHistograms(const AliVCluster * clus, const EventFlags & eflags)
{
	Int_t x, z;
	Int_t sm = CheckClusterGetSM(clus, x, z);

	TLorentzVector p;
	clus->GetMomentum(p, eflags.vtxBest);

	Float_t eff = 1. / fWeights->TofEfficiency(p.E());
	fSpectrum->FillAll(sm, sm, p.Pt(), eff);

	Bool_t cpv = clus->GetEmcCpvDistance() > fDistanceCPV;
	if (cpv)
		fSpectrumCPV->FillAll(sm, sm, p.Pt(), eff);

	Bool_t disp = TestLambda(clus->GetM20(), clus->GetM02(), fDispersionCut);

	if (disp)
		fSpectrumDisp->FillAll(sm, sm, p.Pt(), eff);

	if (cpv && disp)
		fSpectrumBoth->FillAll(sm, sm, p.Pt(), eff);

}
