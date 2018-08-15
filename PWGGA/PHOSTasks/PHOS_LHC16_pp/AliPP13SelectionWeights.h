#ifndef ALIPP13SELECTIONWEIGHTS_H
#define ALIPP13SELECTIONWEIGHTS_H


// --- Custom header files ---

// --- ROOT system ---
#include <TLorentzVector.h>
#include <TObject.h>
#include <TF1.h>
#include <TH1.h>
#include <TList.h>

// --- AliRoot header files ---
#include <AliAODMCParticle.h>
#include <AliVCluster.h>
#include <AliStack.h>
#include <AliLog.h>

#include <AliPIDResponse.h>
#include <AliPHOSGeometry.h>
#include <AliVCluster.h>
#include <AliLog.h>


struct EventFlags
{
	enum EventType {kMB = 0, kGood = 1, kZvertex = 2, kNcontributors = 3, kTwoPhotons = 4};

	EventFlags(Int_t c = 0, Int_t z = 0, Bool_t m = kFALSE, Bool_t p = kFALSE, Bool_t vtx = kFALSE, UShort_t bc = 0. /*, Bool_t v0 = kFalse*/):
		centr(c),
		zvtx(z),
		BC(bc),
		isMixing(m),
		eventPileup(p),
		eventVtxExists(vtx),
		ncontributors(0),
		fMcParticles(0),
		fPIDResponse(0)
		//, eventV0AND(v0)
	{}

	Double_t vtxBest[3];   // Calculated vertex position
	Int_t  centr;
	Int_t  zvtx;
	UShort_t BC;
	Bool_t isMixing;
	Bool_t eventPileup;
	Bool_t eventVtxExists;
	Int_t ncontributors;
	TClonesArray * fMcParticles;
	AliPIDResponse * fPIDResponse;
	// Bool_t eventV0AND;
};


struct AliPP13SelectionWeights: TObject
{
	enum Mode {kData, kMC, kFeeddown, kSinglePi0MC, kSingleEtaMC, kPlain};

	// NB: One needs default constructor for IO readsons
	AliPP13SelectionWeights(): TObject() {}

	virtual Double_t TofEfficiency(Double_t x) const
	{
		(void) x;
		return 1.0;
	}

	virtual Double_t Weights(Double_t x, const EventFlags & eflags) const
	{
		(void) x;
		(void) eflags;
		return 1.0;
	}

	virtual Double_t Nonlinearity(Double_t x) const
	{
		(void) x;
		return 1.0;
	}

	virtual void Report(TList * listOfHistos) const
	{
		TString weights = "No weights";
		listOfHistos->AddFirst(
			new TH1C("selection_weighs", weights, 1, 0, 1)
		); 
	}

	// TODO: Use Shared_ptr
	static AliPP13SelectionWeights & Init(Mode m);
protected:
	ClassDef(AliPP13SelectionWeights, 2)
};

struct AliPP13SelectionWeightsTOF: public AliPP13SelectionWeights
{
	// NB: One needs default constructor for IO readsons
	AliPP13SelectionWeightsTOF(
	    Double_t la = -1.02802e+00,
	    Double_t lb = 8.51857e+00,
	    Double_t ls = 6.19420e-01,
	    Double_t ea = 1.30498e+00,
	    Double_t eal = -1.78716e+00

	):
		AliPP13SelectionWeights(),
		fLogA(la),
		fLogB(lb),
		fLogScale(ls),
		fExpA(ea),
		fExpAlpha(eal)
	{
	}


	virtual Double_t TofEfficiency(Double_t energy) const;
	virtual void Report(TList * listOfHistos) const
	{
		const char * lab = "TOF logA = %.4g, logB = %.4g, logS = %.4g, ExpA = %.4g, ExpAlpha = %.4g";
		TString weights = Form(lab, fLogA, fLogB, fLogScale, fExpA, fExpAlpha);
		listOfHistos->AddFirst(
			new TH1C("selection_tof", weights, 1, 0, 1)
		);

	}
	// Parameters for TOF cut efficiency
	Double_t fLogA;
	Double_t fLogB;
	Double_t fLogScale;
	Double_t fExpA;
	Double_t fExpAlpha;

protected:
	ClassDef(AliPP13SelectionWeightsTOF, 2)

};

struct AliPP13SelectionWeightsMC: public AliPP13SelectionWeights
{
	// NB: One needs default constructor for IO readsons
	AliPP13SelectionWeightsMC(Double_t a = -0.035, Double_t s = 0.95, Double_t g = 1.02):
		AliPP13SelectionWeights(),
		fNonGlobal(g),
		fNonA(a),
		fNonSigma(s)
	{
	}

	virtual Double_t Nonlinearity(Double_t x) const;
	virtual void Report(TList * listOfHistos) const
	{
		const char * lab = "Nonlinearity parameters; NonGlobal = %.6g, NonA = %.6g, NonSigma = %.6g";
		TString weights = Form(lab, fNonGlobal, fNonA, fNonSigma);
		listOfHistos->AddFirst(
			new TH1C("selection_nonlinearity", weights, 1, 0, 1)
		);
	}

	// Parameters for Nonlinearity
	Double_t fNonGlobal;
	Double_t fNonA;
	Double_t fNonSigma;

protected:
	ClassDef(AliPP13SelectionWeightsMC, 2)

};

struct AliPP13SelectionWeightsFeeddown: public AliPP13SelectionWeightsMC
{
	// NB: One needs default constructor for IO readsons
	AliPP13SelectionWeightsFeeddown(Double_t a = -0.035, Double_t s = 0.95, Double_t g = 1.02):
		AliPP13SelectionWeightsMC(a, s, g),
		fDataMCRatio(0)
	{
		fDataMCRatio =  new TF1(
			"feeddown_ratio",
            "[3] * x *(([4] + [5]) * x - [5]) + [0] * (1 + [1] * TMath::Exp(-x * x/ [2]))", 0, 100);
        fDataMCRatio->SetParameter(0, 1.53561e+00);
        fDataMCRatio->SetParameter(1, -4.69350e-01);
        fDataMCRatio->SetParameter(2, 2.38042e-01);
        fDataMCRatio->SetParameter(3, -8.01155e-02);
        fDataMCRatio->SetParameter(4, 6.30860e-01);
        fDataMCRatio->SetParameter(5, -7.21683e-01);
	}

	virtual Double_t Weights(Double_t x, const EventFlags & eflags) const;
	virtual void Report(TList * listOfHistos) const
	{
		const char * lab = "Nonlinearity parameters; NonGlobal = %.6g, NonA = %.6g, NonSigma = %.6g";
		TString weights = Form(lab, fNonGlobal, fNonA, fNonSigma);
		listOfHistos->AddFirst(
			new TH1C("selection_nonlinearity", weights, 1, 0, 1)
		);
		fDataMCRatio->SetNpx(1000);
		TH1 * mcratio = fDataMCRatio->GetHistogram();
		mcratio->SetName("selection_weights");
		listOfHistos->AddFirst(mcratio);
	}

	// Parameters for Nonlinearity
	TF1 * fDataMCRatio;

protected:
	ClassDef(AliPP13SelectionWeightsFeeddown, 2)

};

struct AliPP13SelectionWeightsSPMC: public AliPP13SelectionWeightsMC
{
	// NB: One needs default constructor for IO readsons
	AliPP13SelectionWeightsSPMC(Double_t a = -0.035, Double_t s = 0.95, Double_t g = 1.02):
		AliPP13SelectionWeightsMC(g, a, s),
		fW0(0.014875782846110793),
		fW1(0.28727403800708634),
		fW2(9.9198075195331),
		fW3(0.135),
		fW4(0.135)
	{
	}
	virtual Double_t Weights(Double_t x, const EventFlags & eflags) const;
	virtual void Report(TList * listOfHistos) const
	{
		const char * lab = "Nonlinearity parameters; NonGlobal = %.6g, NonA = %.6g, NonSigma = %.6g";
		TString weights = Form(lab, fNonGlobal, fNonA, fNonSigma);
		listOfHistos->AddFirst(
			new TH1C("selection_nonlinearity", weights, 1, 0, 1)
		);

		const char * labt = "Tsallis parameters for Single Particle MC; fW0 = %.4g, fW1 = %.4g, fW2 = %.4g, fW3 = %.4g, fW4 = %.4g";
		TString weightst = Form(labt, fW0, fW1, fW2, fW3, fW4);
		listOfHistos->AddFirst(
			new TH1C("selection_weighs", weightst, 1, 0, 1)
		);

	}
	// TODO: Use Shared_ptr
	static AliPP13SelectionWeights & SinglePi0();
	static AliPP13SelectionWeights & SingleEta();

	// Parameters for Nonlinearity function
	Double_t fW0;
	Double_t fW1;
	Double_t fW2;
	Double_t fW3;
	Double_t fW4;

protected:
	ClassDef(AliPP13SelectionWeightsSPMC, 2)
};

#endif