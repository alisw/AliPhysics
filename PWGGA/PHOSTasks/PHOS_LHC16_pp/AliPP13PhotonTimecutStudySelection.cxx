// --- Custom header files ---
#include "AliPP13PhotonTimecutStudySelection.h"

// --- ROOT system ---
#include <TH2F.h>

// --- AliRoot header files ---
#include <AliVCluster.h>
#include <AliLog.h>

#include <iostream>
using namespace std;

ClassImp(AliPP13PhotonTimecutStudySelection);


//________________________________________________________________
Bool_t AliPP13PhotonTimecutStudySelection::IsMainBC(const AliVCluster * clus) const
{
	return TMath::Abs(clus->GetTOF()) < fTimingCutPair;
}

//________________________________________________________________
void AliPP13PhotonTimecutStudySelection::InitSelectionHistograms()
{

	// pi0 mass spectrum
	Int_t nM       = fLimits.nM;
	Double_t mMin  = fLimits.mMin;
	Double_t mMax  = fLimits.mMax;
	Int_t nPt      = fLimits.nPt;
	Double_t ptMin = fLimits.ptMin;
	Double_t ptMax = fLimits.ptMax;


	// Timecut study
	for (Int_t i = 0; i < 2; ++i)
	{
		const char * s = (i == 0) ? "" : "Mix";
		AliPP13DetectorHistogram::Mode m = AliPP13DetectorHistogram::kSingleHist;
		fMassPt[i]             = new AliPP13DetectorHistogram(new TH2F(Form("h%sMassPt", s), "(M,p_{T})_{#gamma#gamma}, N_{cell}>2 ; M_{#gamma#gamma} (GeV/#it{c}^{2}); p_{T} (GeV/#it{c})", nM, mMin, mMax, nPt, ptMin, ptMax), fListOfHistos, m);
		fMassPtMainMain[i]     = new AliPP13DetectorHistogram(new TH2F(Form("h%sMassPtMainMain", s), "(M,p_{T})_{#gamma#gamma}, main-main ; M_{#gamma#gamma} (GeV/#it{c}^{2}); p_{T} (GeV/#it{c})", nM, mMin, mMax, nPt, ptMin, ptMax), fListOfHistos, m);
		fMassPtMainPileup[i]   = new AliPP13DetectorHistogram(new TH2F(Form("h%sMassPtMainPileup", s), "(M,p_{T})_{#gamma#gamma}, main-pileup ; M_{#gamma#gamma} (GeV/#it{c}^{2}); p_{T} (GeV/#it{c})", nM, mMin, mMax, nPt, ptMin, ptMax), fListOfHistos, m);
		fMassPtPileupPileup[i] = new AliPP13DetectorHistogram(new TH2F(Form("h%sMassPtPileupPileup", s), "(M,p_{T})_{#gamma#gamma}, pileup-pileup ; M_{#gamma#gamma} (GeV/#it{c}^{2}); p_{T} (GeV/#it{c})", nM, mMin, mMax, nPt, ptMin, ptMax), fListOfHistos, m);
	}



	// Use these histograms for analysis
	for (Int_t i = 0; i < fListOfHistos->GetEntries(); ++i)
	{
		TH1 * hist = dynamic_cast<TH1 *>(fListOfHistos->At(i));
		if (!hist) continue;
		hist->Sumw2();
	}
}

//________________________________________________________________
void AliPP13PhotonTimecutStudySelection::ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags)
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

	Double_t ma12 = psum.M();
	Double_t pt12 = psum.Pt();

	Int_t mix = eflags.isMixing;
	fMassPt[mix]->FillAll(sm1, sm2, ma12, pt12);

	Bool_t bc1 = IsMainBC(c1);
	Bool_t bc2 = IsMainBC(c2);

	if (bc1 && bc1)
		fMassPtMainMain[mix]->FillAll(sm1, sm2, ma12, pt12);

	if (bc1 ^ bc2)
		fMassPtMainPileup[mix]->FillAll(sm1, sm2, ma12, pt12);

	if ((!bc1) && (!bc2))
		fMassPtPileupPileup[mix]->FillAll(sm1, sm2, ma12, pt12);
}
