#ifndef ALIPP13MESONSELECTIONMC_H
#define ALIPP13MESONSELECTIONMC_H


#include <map>

// --- Custom header files ---
#include "AliPP13SpectrumSelectionMC.h"
#include "AliPP13ParticlesHistogram.h"
#include "AliPP13SelectionWeights.h"

// --- ROOT system ---
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TList.h>
#include <TF1.h>

// --- AliRoot header files ---
#include <AliAODMCParticle.h>
#include <AliVCluster.h>
#include <AliStack.h>
#include <AliLog.h>


// NB: This will simplify the code
//

// TODO: Split this class to have separate efficiency and contamination estimators
// 

struct ParticleSpectrum
{
	ParticleSpectrum(const char * n, TList * fListOfHistos, Int_t ptsize, Float_t * ptbins, Bool_t full = kTRUE):
		fPtAllRange(0),
		fPtRadius(0),
		fEtaPhi(0),
		fPtLong(0),
		fPt(0),
		fPtPrimaries(),
		fPtPrimariesStandard()
	{
		fPtAllRange = new TH1F(Form("hPt_allrange_%s", n), Form("Generated p_{T} spectrum of %ss in 4 #pi ; p_{T} (GeV/#it{c})", n), ptsize, ptbins);
		fPtRadius   = new TH2F(Form("hPt_%s_radius", n), Form("Generated radius, p_{T} spectrum of all %ss; r, cm; p_{T} (GeV/#it{c})", n), 500, 0., 500., 400, 0, 20);
		fEtaPhi     = new TH2F(Form("hEtaPhi_%s", n), Form("Generated %ss y vs #phi plot; #phi (rad); y", n), 100, 0, TMath::Pi() * 2, 100, -1, 1);
		fPtLong     = new TH1F(Form("hPtLong_%s", n), Form("Generated p_{T} spectrum of %ss; p_{T} (GeV/#it{c})", n), 1000, 0, 100);
		fPt         = new TH1F(Form("hPt_%s", n), Form("Generated p_{T} spectrum of %ss; p_{T} (GeV/#it{c})", n), ptsize, ptbins);

		fListOfHistos->Add(fPtAllRange);
		fListOfHistos->Add(fPtRadius);
		fListOfHistos->Add(fEtaPhi);
		fListOfHistos->Add(fPtLong);
		fListOfHistos->Add(fPt);

		if (!full)
			return;

		for(Int_t i = 0; i < 2; ++i)
		{
			const char * s = (i == 0) ? "secondary": "primary";
			fPtPrimaries[i] = new TH1F(Form("hPt_%s_%s_", n, s), Form("Generated p_{T} spectrum of %s %ss; p_{T} (GeV/#it{c})", s, n), ptsize, ptbins);
			fListOfHistos->Add(fPtPrimaries[i]);

			fPtPrimariesStandard[i] = new TH1F(Form("hPt_%s_%s_standard", n, s), Form("Generated p_{T} spectrum of %s %ss; p_{T} (GeV/#it{c})", s, n), 200, 0, 20);
			fListOfHistos->Add(fPtPrimariesStandard[i]);

		}
	}

// private:

	TH1F * fPtAllRange; //!
	TH2F * fPtRadius;   //!
	TH2F * fEtaPhi;     //!
	TH1F * fPtLong;     //!
	TH1F * fPt;         //!
	TH1F * fPtPrimaries[2]; //!
	TH1F * fPtPrimariesStandard[2]; //!

};


class AliPP13MesonSelectionMC: public AliPP13SpectrumSelectionMC
{
public:
	enum Modes {kGenerated = 0, kReconstructed = 1, kNhists = 2};
	enum Particles
	{
		kGamma = 22, kPi0 = 111, kEta = 221, kK0s = 310,
		kOmega = 223, kLambda = 3122, kPPion = 211, kNPion = -211,
		kPRho = 213, kNRho = -213,
		kKStarP = 323, kKStarN = -323, kKStar0 = 313, kBarKstar0 = -313,
		kKplus = 321, kKminus = -321, kSigmaZero = 3212
	};

	AliPP13MesonSelectionMC():
		AliPP13SpectrumSelectionMC(),
		fPrimaryPi0(),
		fSecondaryPi0(),
		fFeedDownPi0(),
		fInvMass(),
		fPi0Sources()
	{
		fPartNames[kGamma] = "#gamma";
		fPartNames[kPi0] = "#pi^{0}";
		fPartNames[kEta] = "#eta";

		// Define sources of pi0s
		fPi0SourcesNames[kPRho] = "#rho^{+}";
		fPi0SourcesNames[kNRho] = "#rho^{-}";
		fPi0SourcesNames[kK0s] = "K^{s}_{0}";
		fPi0SourcesNames[kLambda] = "#Lambda";
		fPi0SourcesNames[kPPion] = "#pi^{+}";
		fPi0SourcesNames[kNPion] = "#pi^{-}";
		fPi0SourcesNames[kEta] = "#eta";
		fPi0SourcesNames[kOmega] = "#omega";
		fPi0SourcesNames[kKStarP] = "K^{*+}";
		fPi0SourcesNames[kKStarN] = "K^{*-}";
		fPi0SourcesNames[kKStar0] = "K^{*0}";
		fPi0SourcesNames[kBarKstar0] = "#barK^{*0}";
		fPi0SourcesNames[kKplus] = "K^{+}";
		fPi0SourcesNames[kKminus] = "K^{-}";
		fPi0SourcesNames[kSigmaZero] = "#Sigma^{0}";
	}

	AliPP13MesonSelectionMC(const char * name, const char * title, 
			AliPP13ClusterCuts cuts, AliPP13SelectionWeights * w):
		AliPP13SpectrumSelectionMC(name, title, cuts, w),
		fPrimaryPi0(),
		fSecondaryPi0(),
		fFeedDownPi0(),
		fInvMass(),
		fPi0Sources()
	{
		// Force no timing cut for MC,
		// as there is no photons from different bunches
		fCuts.fTimingCut = 9999; 

		// Don't use c++11 here, as it might fail at some nodes
		fPartNames[kGamma] = "#gamma";
		fPartNames[kPi0] = "#pi^{0}";
		fPartNames[kEta] = "#eta";

		// Define sources of pi0s
		fPi0SourcesNames[kPRho] = "#rho^{+}";
		fPi0SourcesNames[kNRho] = "#rho^{-}";
		fPi0SourcesNames[kK0s] = "K^{s}_{0}";
		fPi0SourcesNames[kLambda] = "#Lambda";
		fPi0SourcesNames[kPPion] = "#pi^{+}";
		fPi0SourcesNames[kNPion] = "#pi^{-}";
		fPi0SourcesNames[kEta] = "#eta";
		fPi0SourcesNames[kOmega] = "#omega";
		fPi0SourcesNames[kKStarP] = "K^{*+}";
		fPi0SourcesNames[kKStarN] = "K^{*-}";
		fPi0SourcesNames[kKStar0] = "K^{*0}";
		fPi0SourcesNames[kBarKstar0] = "#barK^{*0}";
		fPi0SourcesNames[kKplus] = "K^{+}";
		fPi0SourcesNames[kKminus] = "K^{-}";
		fPi0SourcesNames[kSigmaZero] = "#Sigma^{0}";

	}

	virtual void InitSelectionHistograms();
	virtual void ConsiderGeneratedParticles(const EventFlags & eflags);

	virtual ~AliPP13MesonSelectionMC()
	{
		for (int i = 0; i < kNhists; ++i)
		{
			if (fPrimaryPi0[i]) delete fPrimaryPi0[i];
			if (fSecondaryPi0[i]) delete fSecondaryPi0[i];
			if (fFeedDownPi0[i]) delete fFeedDownPi0[i];
		}

		for (ParticleSpectrums::iterator i = fSpectrums.begin(); i != fSpectrums.end(); ++i)
			delete i->second;
	}

protected:
	virtual void ConsiderPair(const AliVCluster * c1, const AliVCluster * c2, const EventFlags & eflags);

	virtual AliAODMCParticle * GetParent(Int_t label, Int_t & plabel, TClonesArray * particles) const;
	virtual AliAODMCParticle * GetParent(Int_t label, TClonesArray * particles) const
	{
		Int_t plabel;
		return GetParent(label, plabel, particles);
	}
	void ConsiderGeneratedPi0(Int_t i, Double_t pt, Bool_t primary, const EventFlags & flags);

	AliPP13MesonSelectionMC(const AliPP13MesonSelectionMC &);
	AliPP13MesonSelectionMC & operator = (const AliPP13MesonSelectionMC &);

	AliPP13ParticlesHistogram * fPrimaryPi0[kNhists];
	AliPP13ParticlesHistogram * fSecondaryPi0[kNhists];
	AliPP13ParticlesHistogram * fFeedDownPi0[kNhists];


	// This data structure contains all necesary histograms
	// for the particles we want to get
	typedef std::map<Int_t, ParticleSpectrum * > ParticleSpectrums;
	ParticleSpectrums fSpectrums;

	typedef std::map<Int_t, TString> EnumNames;
	EnumNames fPartNames;
	EnumNames fPi0SourcesNames;

	TH1 * fInvMass[2];     //!
	TH1 * fPi0Sources[2];  //!

	// Parameters of weighed MC parametrization
	ClassDef(AliPP13MesonSelectionMC, 2)
};
#endif
