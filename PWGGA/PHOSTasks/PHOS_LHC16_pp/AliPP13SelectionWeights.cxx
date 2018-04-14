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
Double_t AliPP13SelectionWeightsTOF::Weight(Double_t energy) const
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
    return fNonGlobal * (1. + fNonA * TMath::Exp(-x * x / 2. / fNonSigma / fNonSigma));
}



//________________________________________________________________
Double_t AliPP13SelectionWeightsSPMC::Weight(Double_t pT) const
{
	Double_t w = pT * pT * fW0 / 2. / TMath::Pi();
	Double_t fraction = (fW2 - 1.) * (fW2 - 2.) / (fW2 * fW1 * (fW2 * fW1 + fW4 * (fW2 - 2.)));
	Double_t power = TMath::Power(1. + (TMath::Sqrt(pT * pT + fW3 * fW3) - fW4) / (fW2 * fW1), -fW2);
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


    // Debug
    ws.fW0 = 21.339890553914014;
    ws.fW1 = 0.08359755308503322;
    ws.fW2 = 7.334946541612603;

    // Weights 0 (new efficiency)
    // ws.fW0 = 0.10325998438001027;
    // ws.fW1 = 0.1710556728057399;
    // ws.fW2 = 8.613628140871766;

    ws.fW3 = 0.135;
    ws.fW4 = 0.135;



    // Nonlinearity Naive estimation
    //
    // ws.fNonA = -0.022934923767457753;
    // ws.fNonSigma = 1.4188237289034245;
    // ws.fNonGlobal = 1.0579663356860527;

     ws.fNonA = -0.014719244288611932;
     ws.fNonSigma = 2 * 0.8017501954719543;
     ws.fNonGlobal = 1.050000000000015;

    // ws.fNonGlobal = 1.0;
    // ws.fNonA = 0;
    // ws.fNonSigma = 1.0579663356860527;

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

     ws.fNonA = -0.014719244288611932;
     ws.fNonSigma = 2 * 0.8017501954719543;
     ws.fNonGlobal = 1.050000000000015;
     return ws;
}


//________________________________________________________________
AliPP13SelectionWeights & AliPP13SelectionWeights::Init(Mode m)
{
    if(m == kSinglePi0MC)
        return AliPP13SelectionWeightsSPMC::SinglePi0();

    if(m == kSingleEtaMC)
        return AliPP13SelectionWeightsSPMC::SingleEta();

    if(m == kMC)
        return * new AliPP13SelectionWeightsMC();

    if(m == kData)
        return * new AliPP13SelectionWeightsTOF();
      
    return * new AliPP13SelectionWeights();
}
