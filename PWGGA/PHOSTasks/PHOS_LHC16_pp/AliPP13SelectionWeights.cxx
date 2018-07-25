// --- Custom header files ---
#include <AliPP13SelectionWeights.h>

// --- ROOT system ---
#include <TMath.h>

// --- AliRoot header files ---


ClassImp(AliPP13SelectionWeights);
ClassImp(AliPP13SelectionWeightsTOF);
ClassImp(AliPP13SelectionWeightsMC);
ClassImp(AliPP13SelectionWeightsSPMC);

//________________________________________________________________
Double_t AliPP13SelectionWeightsTOF::TofEfficiency(Double_t energy) const
{
    // TOF efficiency was parametrized as photon energy
    //
    Double_t logisitc = fLogScale / (1. + TMath::Exp(energy * fLogA + fLogB));
    Double_t expo = fExpA * TMath::Exp(energy * fExpAlpha);

    return 1. - logisitc - expo;
}


//________________________________________________________________
Double_t AliPP13SelectionWeightsMC::Nonlinearity(Double_t x) const
{
    return fNonGlobal * (1. + fNonA / (1 + TMath::Power(x / fNonSigma, 2)));
}



//________________________________________________________________
Double_t AliPP13SelectionWeightsSPMC::Weights(Double_t pT, const EventFlags & eflags) const
{
    // NB: Don't use origin pT
    (void) pT;
    AliAODMCParticle * origin = (AliAODMCParticle*)eflags.fMcParticles->At(0);//0 is always generated particle by AliGenBox.
    Double_t opT = origin->Pt();

    // NB: Try generating the yield instead of invariant yield
    Double_t w = /* opT * */ opT * fW0 / 2. / TMath::Pi();
    Double_t fraction = (fW2 - 1.) * (fW2 - 2.) / (fW2 * fW1 * (fW2 * fW1 + fW4 * (fW2 - 2.)));
    Double_t power = TMath::Power(1. + (TMath::Sqrt(opT * opT + fW3 * fW3) - fW4) / (fW2 * fW1), -fW2);
    return w * fraction * power;
}

//________________________________________________________________
AliPP13SelectionWeights & AliPP13SelectionWeightsSPMC::SinglePi0()
{
    AliPP13SelectionWeightsSPMC & ws = * new AliPP13SelectionWeightsSPMC();

    // Weights 3
    // ws.fW0 = 0.014875782846110793;
    // ws.fW1 = 0.28727403800708634;
    // ws.fW2 = 9.9198075195331;

    // Weights 0 (new efficiency)
    // ws.fW0 = 0.10325998438001027;
    // ws.fW1 = 0.1710556728057399;
    // ws.fW2 = 8.613628140871766;

    // Debug
    // ws.fW0 = 21.339890553914014;
    // ws.fW1 = 0.08359755308503322;
    // ws.fW2 = 7.334946541612603;

    // The latest iteration
    //
    ws.fW0 = 0.2622666606436988 / 0.0119143016137;
    ws.fW1 = 0.08435275173194286;
    ws.fW2 = 7.356520553419461;
    ws.fW3 = 0.135;
    ws.fW4 = 0.135;

    ws.fNonA = -0.06;
    ws.fNonSigma = 0.7;
    ws.fNonGlobal = 1.015;
    return ws;
}


//________________________________________________________________
AliPP13SelectionWeights & AliPP13SelectionWeightsSPMC::SingleEta()
{
    AliPP13SelectionWeightsSPMC & ws = * new AliPP13SelectionWeightsSPMC();

    // NB: Note Different Parameters
    // Weights Initial (form 7 TeV paper)
    ws.fW0 = 0.201;
    ws.fW1 = 0.229;
    ws.fW2 = 7.0;
    ws.fW3 = 0.547;
    ws.fW4 = 0.547;

    // The latest nonlinarity tested on the simples data
    ws.fNonA = -0.06;
    ws.fNonSigma = 0.7;
    ws.fNonGlobal = 1.015;
    return ws;
}


//________________________________________________________________
AliPP13SelectionWeights & AliPP13SelectionWeights::Init(Mode m)
{
    if (m == kSinglePi0MC)
        return AliPP13SelectionWeightsSPMC::SinglePi0();

    if (m == kSingleEtaMC)
        return AliPP13SelectionWeightsSPMC::SingleEta();

    if (m == kMC)
        return * new AliPP13SelectionWeightsMC();

    if (m == kData)
        return * new AliPP13SelectionWeightsTOF();

    return * new AliPP13SelectionWeights();
}
