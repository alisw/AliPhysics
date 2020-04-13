// --- Custom header files ---
#include <AliPP13SelectionWeights.h>

// --- ROOT system ---
#include <TMath.h>

// --- AliRoot header files ---


ClassImp(AliPP13SelectionWeights);
ClassImp(AliPP13SelectionWeightsTOF);
ClassImp(AliPP13SelectionWeightsMC);
ClassImp(AliPP13SelectionWeightsFeeddown);
ClassImp(AliPP13SelectionWeightsSPMC);
ClassImp(AliPP13SelectionWeightsScan);

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
Double_t AliPP13SelectionWeightsScan::Nonlinearity(Double_t x) const
{
    // These magic numbers are taken from the official tender configuration
    // https://github.com/alisw/AliPhysics/blob/master/TENDER/TenderSupplies/AliPHOSTenderSupply.cxx

    // When changing these parameters chagne the fE and fD as well
    Double_t p0 = 1.04397;
    Double_t p1 = 0.512307;
    Double_t p2 = 0.133812;
    Double_t p3 = -0.150093;
    Double_t p4 = -0.455062;

    // Override the parameters for scan
    p2 = fE;
    p4 = fD;

    // Correct in the following way: p *= nonlin(E);
    return p0 + p1 / x + p2 / x / x + p3 / TMath::Sqrt(x) + p4 / x / TMath::Sqrt(x);
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
Double_t AliPP13SelectionWeightsFeeddown::Weights(Double_t pT, const EventFlags & eflags) const
{
    (void) eflags;
    return fDataMCRatio->Eval(pT);
}

//________________________________________________________________
AliPP13SelectionWeights & AliPP13SelectionWeightsSPMC::SinglePi0()
{
    AliPP13SelectionWeightsSPMC & ws = * new AliPP13SelectionWeightsSPMC();
    // The latest iteration
    //
    ws.fW0 = 0.2622666606436988 / 0.0119143016137;
    ws.fW1 = 0.08435275173194286;
    ws.fW2 = 7.356520553419461;
    ws.fW3 = 0.135;
    ws.fW4 = 0.135;

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

    return ws;
}


//________________________________________________________________
AliPP13SelectionWeights & AliPP13SelectionWeights::Init(Mode m)
{
    if (m == kSinglePi0MC)
        return AliPP13SelectionWeightsSPMC::SinglePi0();

    if (m == kSingleEtaMC)
        return AliPP13SelectionWeightsSPMC::SingleEta();

    if (m == kScan)
        return * new AliPP13SelectionWeightsScan();

    if (m == kFeeddown)
        return * new AliPP13SelectionWeightsFeeddown();

    if (m == kMC)
        return * new AliPP13SelectionWeightsMC();

    if (m == kData)
        return * new AliPP13SelectionWeightsTOF();

    return * new AliPP13SelectionWeights();
}
