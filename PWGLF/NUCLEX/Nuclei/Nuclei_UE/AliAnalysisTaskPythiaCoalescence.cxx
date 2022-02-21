#include "AliAnalysisTaskPythiaCoalescence.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTask.h"
#include "TLorentzVector.h"
#include "AliMCParticle.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;
ClassImp(AliAnalysisTaskPythiaCoalescence)

    //____________________________________________________________________________________________________________________________________________________
    AliAnalysisTaskPythiaCoalescence::AliAnalysisTaskPythiaCoalescence() : AliAnalysisTaskSE(),
fAODevent(nullptr),
fMCEvent(nullptr),
fOutputList(nullptr),
fQAList(nullptr),
fAverage_Nch_Transv(7.2),
hProtonWeights(nullptr),
fProtonWeights(nullptr),
hWeightToward(nullptr),
hWeightTransv(nullptr),
hWeightAway(nullptr),
fDeuteronWF(nullptr),
hSourceSize(nullptr),
hNumberOfEvents(nullptr),
hTransverseMult(nullptr),
hRtDistribution(nullptr),
hProtonsINELgtZERO(nullptr),
hProtonsINELgtZERO_reshaped(nullptr),
hDeuteronsINELgtZERO_TruewignerGaus(nullptr),
hDeuteronsINELgtZERO_TruewignerDoubleGaus(nullptr),
hProtons_Toward(nullptr),
hProtons_Transv(nullptr),
hProtons_Away(nullptr),
hProtons_Toward_reshaped(nullptr),
hProtons_Transv_reshaped(nullptr),
hProtons_Away_reshaped(nullptr),
hNeutronsINELgtZERO(nullptr),
hNeutronsINELgtZERO_reshaped(nullptr),
hNeutrons_Toward(nullptr),
hNeutrons_Transv(nullptr),
hNeutrons_Away(nullptr),
hNeutrons_Toward_reshaped(nullptr),
hNeutrons_Transv_reshaped(nullptr),
hNeutrons_Away_reshaped(nullptr),
hRapidityProtons(nullptr),
hRapidityNeutrons(nullptr),
hDeltaP(nullptr),
hSourceRadius_Prot(nullptr),
hSourceRadius_Neut(nullptr),
hDistanceLab(nullptr),
hDistanceDeut(nullptr),
hDistanceDiff(nullptr),
hPtProtonsFirstBinDeut(nullptr),
hDeltaPhi_Toward(nullptr),
hDeltaPhi_Transv(nullptr),
hDeltaPhi_Away(nullptr),
hDeltaPhi_INELgtZERO(nullptr)
{
}
//____________________________________________________________________________________________________________________________________________________
AliAnalysisTaskPythiaCoalescence::AliAnalysisTaskPythiaCoalescence(const char *name) : AliAnalysisTaskSE(name),
fAODevent(nullptr),
fMCEvent(nullptr),
fOutputList(nullptr),
fQAList(nullptr),
fAverage_Nch_Transv(7.2),
hProtonWeights(nullptr),
fProtonWeights(nullptr),
hWeightToward(nullptr),
hWeightTransv(nullptr),
hWeightAway(nullptr),
fDeuteronWF(nullptr),
hSourceSize(nullptr),
hNumberOfEvents(nullptr),
hTransverseMult(nullptr),
hRtDistribution(nullptr),
hProtonsINELgtZERO(nullptr),
hProtonsINELgtZERO_reshaped(nullptr),
hDeuteronsINELgtZERO_TruewignerGaus(nullptr),
hDeuteronsINELgtZERO_TruewignerDoubleGaus(nullptr),
hProtons_Toward(nullptr),
hProtons_Transv(nullptr),
hProtons_Away(nullptr),
hProtons_Toward_reshaped(nullptr),
hProtons_Transv_reshaped(nullptr),
hProtons_Away_reshaped(nullptr),
hNeutronsINELgtZERO(nullptr),
hNeutronsINELgtZERO_reshaped(nullptr),
hNeutrons_Toward(nullptr),
hNeutrons_Transv(nullptr),
hNeutrons_Away(nullptr),
hNeutrons_Toward_reshaped(nullptr),
hNeutrons_Transv_reshaped(nullptr),
hNeutrons_Away_reshaped(nullptr),
hRapidityProtons(nullptr),
hRapidityNeutrons(nullptr),
hDeltaP(nullptr),
hSourceRadius_Prot(nullptr),
hSourceRadius_Neut(nullptr),
hDistanceLab(nullptr),
hDistanceDeut(nullptr),
hDistanceDiff(nullptr),
hPtProtonsFirstBinDeut(nullptr),
hDeltaPhi_Toward(nullptr),
hDeltaPhi_Transv(nullptr),
hDeltaPhi_Away(nullptr),
hDeltaPhi_INELgtZERO(nullptr)
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}
//____________________________________________________________________________________________________________________________________________________
AliAnalysisTaskPythiaCoalescence::~AliAnalysisTaskPythiaCoalescence()
{

    fOutputList->Clear();
    delete fAODevent;
    delete fMCEvent;
    delete fOutputList;
    delete fQAList;
    delete hProtonWeights;
    delete fProtonWeights;
    delete fDeuteronWF;
    delete hWeightToward;
    delete hWeightTransv;
    delete hWeightAway;
    delete hSourceSize;
    delete hNumberOfEvents;
    delete hTransverseMult;
    delete hRtDistribution;
    delete hProtonsINELgtZERO;
    delete hProtonsINELgtZERO_reshaped;
    delete hDeuteronsINELgtZERO_TruewignerGaus;
    delete hDeuteronsINELgtZERO_TruewignerDoubleGaus;
    delete hProtons_Toward;
    delete hProtons_Transv;
    delete hProtons_Away;
    delete hProtons_Toward_reshaped;
    delete hProtons_Transv_reshaped;
    delete hProtons_Away_reshaped;
    delete hNeutronsINELgtZERO;
    delete hNeutronsINELgtZERO_reshaped;
    delete hNeutrons_Toward;
    delete hNeutrons_Transv;
    delete hNeutrons_Away;
    delete hNeutrons_Toward_reshaped;
    delete hNeutrons_Transv_reshaped;
    delete hNeutrons_Away_reshaped;
    delete hRapidityProtons;
    delete hRapidityNeutrons;
    delete hDeltaP;
    delete hSourceRadius_Prot;
    delete hSourceRadius_Neut;
    delete hDistanceLab;
    delete hDistanceDeut;
    delete hDistanceDiff;
    delete hPtProtonsFirstBinDeut;
    delete hDeltaPhi_Toward;
    delete hDeltaPhi_Transv;
    delete hDeltaPhi_Away;
    delete hDeltaPhi_INELgtZERO;
}
//____________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskPythiaCoalescence::UserCreateOutputObjects()
{

    // Create Output List
    fOutputList = new TList();
    fQAList = new TList();
    fOutputList->SetOwner();
    fQAList->SetOwner();

    // Event Counter
    hNumberOfEvents = new TH1D("hNumberOfEvents", "", 20, 0, 20);
    hNumberOfEvents->Sumw2();
    fOutputList->Add(hNumberOfEvents);

    // Transverse Multiplicity
    hTransverseMult = new TH1D("hTransverseMult", "", 200, 0, 200);
    hTransverseMult->Sumw2();
    fOutputList->Add(hTransverseMult);

    // R_{T} Distribution
    hRtDistribution = new TH1D("hRtDistribution", "", 1000, 0, 10);
    hRtDistribution->Sumw2();
    fOutputList->Add(hRtDistribution);

    // p_{T} Intervals
    Double_t pt_proton[] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5,
                            1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 10.0};
    Double_t pt_deuteron[] = {0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.6, 3.0, 3.4, 3.8};
    const Int_t nPtProton = sizeof(pt_proton) / sizeof(Double_t) - 1;
    const Int_t nPtDeuteron = sizeof(pt_deuteron) / sizeof(Double_t) - 1;

    // p_{T} Spectra in Events INEL>0
    hProtonsINELgtZERO = new TH1D("hProtonsINELgtZERO", "", nPtProton, pt_proton);
    hProtonsINELgtZERO_reshaped = new TH1D("hProtonsINELgtZERO_reshaped", "", nPtProton, pt_proton);
    hNeutronsINELgtZERO = new TH1D("hNeutronsINELgtZERO", "", nPtProton, pt_proton);
    hNeutronsINELgtZERO_reshaped = new TH1D("hNeutronsINELgtZERO_reshaped", "", nPtProton, pt_proton);
    hProtonsINELgtZERO->Sumw2();
    hNeutronsINELgtZERO->Sumw2();
    hProtonsINELgtZERO_reshaped->Sumw2();
    hNeutronsINELgtZERO_reshaped->Sumw2();
    fOutputList->Add(hProtonsINELgtZERO);
    fOutputList->Add(hNeutronsINELgtZERO);
    fOutputList->Add(hProtonsINELgtZERO_reshaped);
    fOutputList->Add(hNeutronsINELgtZERO_reshaped);

    // p_{T} Spectra Deuterons
    hDeuteronsINELgtZERO_TruewignerGaus = new TH1D("hDeuteronsINELgtZERO_TruewignerGaus", "Distance: 1.2 fm", nPtDeuteron, pt_deuteron);
    hDeuteronsINELgtZERO_TruewignerGaus->Sumw2();
    fOutputList->Add(hDeuteronsINELgtZERO_TruewignerGaus);

    hDeuteronsINELgtZERO_TruewignerDoubleGaus = new TH1D("hDeuteronsINELgtZERO_TruewignerDoubleGaus", "Distance: 1.2 fm", nPtDeuteron, pt_deuteron); // min bias
    hDeuteronsINELgtZERO_TruewignerDoubleGaus->Sumw2();
    fOutputList->Add(hDeuteronsINELgtZERO_TruewignerDoubleGaus);

    for (Int_t i = 0; i < 50; i++)
    {

        // p_{T} Spectra: Deuterons (Simple Coalescence)
        hDeuteronsINELgtZERO_simpleCoal[i] = new TH1D(Form("hDeuteronsINELgtZERO_simpleCoal[%d]", i), "", nPtDeuteron, pt_deuteron);
        hDeuterons_Toward_simpleCoal[i] = new TH1D(Form("hDeuterons_Toward_simpleCoal[%d]", i), "", 500, 0, 5);
        hDeuterons_Transv_simpleCoal[i] = new TH1D(Form("hDeuterons_Transv_simpleCoal[%d]", i), "", 500, 0, 5);
        hDeuteronsINELgtZERO_simpleCoal[i]->Sumw2();
        hDeuterons_Toward_simpleCoal[i]->Sumw2();
        hDeuterons_Transv_simpleCoal[i]->Sumw2();
        fOutputList->Add(hDeuteronsINELgtZERO_simpleCoal[i]);
        fOutputList->Add(hDeuterons_Toward_simpleCoal[i]);
        fOutputList->Add(hDeuterons_Transv_simpleCoal[i]);

        // p_{T} Spectra: Deuterons (Wigner Gaussian)
        hDeuteronsINELgtZERO_wignerGaus[i] = new TH1D(Form("hDeuteronsINELgtZERO_wignerGaus[%d]", i), "", nPtDeuteron, pt_deuteron);
        hDeuterons_Toward_wignerGaus[i] = new TH1D(Form("hDeuterons_Toward_wignerGaus[%d]", i), "", 500, 0, 5);
        hDeuterons_Transv_wignerGaus[i] = new TH1D(Form("hDeuterons_Transv_wignerGaus[%d]", i), "", 500, 0, 5);
        hDeuteronsINELgtZERO_wignerGaus[i]->Sumw2();
        hDeuterons_Toward_wignerGaus[i]->Sumw2();
        hDeuterons_Transv_wignerGaus[i]->Sumw2();
        fOutputList->Add(hDeuteronsINELgtZERO_wignerGaus[i]);
        fOutputList->Add(hDeuterons_Toward_wignerGaus[i]);
        fOutputList->Add(hDeuterons_Transv_wignerGaus[i]);

        // p_{T} Spectra: Deuterons (Wigner Argonne)
        hDeuteronsINELgtZERO_wignerArg[i] = new TH1D(Form("hDeuteronsINELgtZERO_wignerArg[%d]", i), "", nPtDeuteron, pt_deuteron);
        hDeuterons_Toward_wignerArg[i] = new TH1D(Form("hDeuterons_Toward_wignerArg[%d]", i), "", 500, 0, 5);
        hDeuterons_Transv_wignerArg[i] = new TH1D(Form("hDeuterons_Transv_wignerArg[%d]", i), "", 500, 0, 5);
        hDeuteronsINELgtZERO_wignerArg[i]->Sumw2();
        hDeuterons_Toward_wignerArg[i]->Sumw2();
        hDeuterons_Transv_wignerArg[i]->Sumw2();
        fOutputList->Add(hDeuteronsINELgtZERO_wignerArg[i]);
        fOutputList->Add(hDeuterons_Toward_wignerArg[i]);
        fOutputList->Add(hDeuterons_Transv_wignerArg[i]);

        // p_{T} Spectra: Deuterons (True Wigner Gaus)
        hDeuterons_Toward_TruewignerGaus[i] = new TH1D(Form("hDeuterons_Toward_TruewignerGaus[%d]", i), Form("Distance: %d*0.05 fm", i ), 500, 0, 5);
        hDeuterons_Transv_TruewignerGaus[i] = new TH1D(Form("hDeuterons_Transv_TruewignerGaus[%d]", i), Form("Distance: %d*0.05 fm", i ), 500, 0, 5);
        hDeuterons_Toward_TruewignerGaus[i]->Sumw2();
        hDeuterons_Transv_TruewignerGaus[i]->Sumw2();
        fOutputList->Add(hDeuterons_Toward_TruewignerGaus[i]);
        fOutputList->Add(hDeuterons_Transv_TruewignerGaus[i]);

        // p_{T} Spectra: Deuterons (True Wigner DoubleGaus)
        hDeuterons_Toward_TruewignerDoubleGaus[i] = new TH1D(Form("hDeuterons_Toward_TruewignerDoubleGaus[%d]", i), Form("Distance: %d*0.05 fm", i ), 500, 0, 5);                      // Toward region
        hDeuterons_Transv_TruewignerDoubleGaus[i] = new TH1D(Form("hDeuterons_Transv_TruewignerDoubleGaus[%d]", i), Form("Distance: %d*0.05 fm", i ), 500, 0, 5);                      // away region
        hDeuterons_Toward_TruewignerDoubleGaus[i]->Sumw2();
        hDeuterons_Transv_TruewignerDoubleGaus[i]->Sumw2();
        fOutputList->Add(hDeuterons_Toward_TruewignerDoubleGaus[i]);
        fOutputList->Add(hDeuterons_Transv_TruewignerDoubleGaus[i]);
    }

    // R_{T} Intervals
    Double_t Rt_Intervals[] = {0.0, 0.5, 1.5, 5.0};
    Int_t nRtIntervals = sizeof(Rt_Intervals) / sizeof(Double_t) - 1;

    // Proton Spectra in the Azimuthal Regions
    hProtons_Toward = new TH1D("hProtons_Toward", "", 500, 0, 5);
    hProtons_Transv = new TH1D("hProtons_Transv", "", 500, 0, 5);
    hProtons_Away = new TH1D("hProtons_Away", "", 500, 0, 5);
    hProtons_Toward_reshaped = new TH1D("hProtons_Toward_reshaped", "", 500, 0, 5);
    hProtons_Transv_reshaped = new TH1D("hProtons_Transv_reshaped", "", 500, 0, 5);
    hProtons_Away_reshaped = new TH1D("hProtons_Away_reshaped", "", 500, 0, 5);
    hProtons_Toward->Sumw2();
    hProtons_Transv->Sumw2();
    hProtons_Away->Sumw2();
    hProtons_Toward_reshaped->Sumw2();
    hProtons_Transv_reshaped->Sumw2();
    hProtons_Away_reshaped->Sumw2();
    fOutputList->Add(hProtons_Toward);
    fOutputList->Add(hProtons_Transv);
    fOutputList->Add(hProtons_Away);
    fOutputList->Add(hProtons_Toward_reshaped);
    fOutputList->Add(hProtons_Transv_reshaped);
    fOutputList->Add(hProtons_Away_reshaped);

    // Neutron Spectra in the Azimuthal Regions
    hNeutrons_Toward = new TH1D("hNeutrons_Toward", "", 500, 0, 5);
    hNeutrons_Transv = new TH1D("hNeutrons_Transv", "", 500, 0, 5);
    hNeutrons_Away = new TH1D("hNeutrons_Away", "", 500, 0, 5);
    hNeutrons_Toward_reshaped = new TH1D("hNeutrons_Toward_reshaped", "", 500, 0, 5);
    hNeutrons_Transv_reshaped = new TH1D("hNeutrons_Transv_reshaped", "", 500, 0, 5);
    hNeutrons_Away_reshaped = new TH1D("hNeutrons_Away_reshaped", "", 500, 0, 5);
    hNeutrons_Toward->Sumw2();
    hNeutrons_Transv->Sumw2();
    hNeutrons_Away->Sumw2();
    hNeutrons_Toward_reshaped->Sumw2();
    hNeutrons_Transv_reshaped->Sumw2();
    hNeutrons_Away_reshaped->Sumw2();
    fOutputList->Add(hNeutrons_Toward);
    fOutputList->Add(hNeutrons_Transv);
    fOutputList->Add(hNeutrons_Away);
    fOutputList->Add(hNeutrons_Toward_reshaped);
    fOutputList->Add(hNeutrons_Transv_reshaped);
    fOutputList->Add(hNeutrons_Away_reshaped);

    // Rapidity Distributions of Protons & Neutrons that form Deuteron in |y|<0.5
    hRapidityProtons = new TH1D("hRapidityProtons", "", 200, -2, 2);
    hRapidityNeutrons = new TH1D("hRapidityNeutrons", "", 200, -2, 2);
    hRapidityProtons->Sumw2();
    hRapidityNeutrons->Sumw2();
    fOutputList->Add(hRapidityProtons);
    fOutputList->Add(hRapidityNeutrons);

    // DeltaP Distribution
    hDeltaP = new TH1D("hDeltaP", "", 1000, 0, 1);
    hDeltaP->Sumw2();
    fOutputList->Add(hDeltaP);

    // Source radii
    hSourceRadius_Prot = new TH1D("hSourceRadius_Prot", "", 500, 0, 5);
    hSourceRadius_Neut = new TH1D("hSourceRadius_Neut", "", 500, 0, 5);
    hSourceRadius_Prot->Sumw2();
    hSourceRadius_Neut->Sumw2();
    fOutputList->Add(hSourceRadius_Prot);
    fOutputList->Add(hSourceRadius_Neut);

    // Control Histograms
    hDistanceLab = new TH1D("hDistanceLab", "", 500, 0, 5);
    hDistanceDeut = new TH1D("hDistanceDeut", "", 500, 0, 5);
    hDistanceDiff = new TH1D("hDistanceDiff", "", 500, -5, 5);
    hDistanceLab->Sumw2();
    hDistanceDeut->Sumw2();
    hDistanceDiff->Sumw2();
    fOutputList->Add(hDistanceLab);
    fOutputList->Add(hDistanceDeut);
    fOutputList->Add(hDistanceDiff);

    // pT Distribution Protons contributing to first bin of Deuterons
    hPtProtonsFirstBinDeut = new TH1D("hPtProtonsFirstBinDeut", "", 100, 0, 1);
    hPtProtonsFirstBinDeut->Sumw2();
    fOutputList->Add(hPtProtonsFirstBinDeut);

    // Angular Distributions
    hDeltaPhi_Toward = new TH2D("hDeltaPhi_Toward", "", 50, 0, 5, 200, 0, 100);
    hDeltaPhi_Transv = new TH2D("hDeltaPhi_Transv", "", 50, 0, 5, 200, 0, 100);
    hDeltaPhi_Away = new TH2D("hDeltaPhi_Away", "", 50, 0, 5, 200, 0, 100);
    hDeltaPhi_INELgtZERO = new TH2D("hDeltaPhi_INELgtZERO", "", 50, 0, 5, 200, 0, 100);
    hDeltaPhi_Toward->Sumw2();
    hDeltaPhi_Transv->Sumw2();
    hDeltaPhi_Away->Sumw2();
    hDeltaPhi_INELgtZERO->Sumw2();
    fOutputList->Add(hDeltaPhi_Toward);
    fOutputList->Add(hDeltaPhi_Transv);
    fOutputList->Add(hDeltaPhi_Away);
    fOutputList->Add(hDeltaPhi_INELgtZERO);

    PostData(1, fOutputList);
    PostData(2, fQAList);
}
//____________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskPythiaCoalescence::UserExec(Option_t *)
{

    // Seed Random Number
    gRandom->SetSeed(0);

    // Coalescence Momentum
    Double_t p0[] = {0.200, 0.210, 0.220, 0.230, 0.240, 0.245, 0.250, 0.255, 0.260, 0.265, 0.270, 0.275, 0.280, 0.284, 0.285, 0.286, 0.290, 0.295, 0.300, 0.305, 0.310, 0.315, 0.320, 0.325, 0.330, 0.335, 0.340, 0.345, 0.350, 0.360, 0.370, 0.380, 0.390, 0.400};
    const Int_t nTrials = sizeof(p0) / sizeof(Double_t);

    // Particle Masses in GeV
    Double_t mp = 0.93827208816; // Proton
    Double_t mn = 0.93956542052; // Neutron
    Double_t md = 1.87561294257; // Deuteron

    // Get Input Event (INEL>0 Selection)
    if (!GetEvent())
        return;

    // Protons and Neutrons IDs
    vector<Int_t> proton_ID;
    vector<Int_t> neutron_ID;
    vector<Int_t> neutron_status;

    // Loop over Generated Particles
    for (Int_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++)
    {

        // Get Primary Particle
        AliMCParticle *particle = (AliMCParticle *)fMCEvent->GetTrack(i);
        if (!particle)
            continue;
        if (!particle->IsPhysicalPrimary())
            continue;
        if (IsInjectedParticle(particle))
            continue;
        if (TMath::Abs(particle->Y()) > 1.0)
            continue;

        // Variables
        Double_t pt = particle->Pt();
        Double_t protWeight = GetProtonWeight(pt);

        // Store Protons and Neutrons (no rapidity cut)
        if (particle->PdgCode() == -2212)
        {
            proton_ID.push_back(i);
        }
        if (particle->PdgCode() == -2112)
        {
            neutron_ID.push_back(i);
            neutron_status.push_back(0);
        }

        // Protons and Neutrons p_{T} Spectra from Pythia
        if (TMath::Abs(particle->Y()) > 0.5)
            continue;
        if (particle->PdgCode() == -2212)
        {
            hProtonsINELgtZERO->Fill(pt);
            hProtonsINELgtZERO_reshaped->Fill(pt, protWeight);
        }
        if (particle->PdgCode() == -2112)
        {
            hNeutronsINELgtZERO->Fill(pt);
            hNeutronsINELgtZERO_reshaped->Fill(pt, protWeight);
        }
    }

    // Save Deuteron Properties (Simple Coalescence)
    Double_t pt_deuteron_simpleCoal[nTrials][20];
    Double_t phi_deuteron_simpleCoal[nTrials][20];
    Double_t weight_deuteron_simpleCoal[nTrials][20];
    Int_t nDeuterons_simpleCoal[nTrials];

    // Save Deuteron Properties (Wigner Gaus)
    Double_t pt_deuteron_wignerGaus[nTrials][20];
    Double_t phi_deuteron_wignerGaus[nTrials][20];
    Double_t weight_deuteron_wignerGaus[nTrials][20];
    Int_t nDeuterons_wignerGaus[nTrials];

    // Save Deuteron Properties (Wigner Argonne18)
    Double_t pt_deuteron_wignerArg[nTrials][20];
    Double_t phi_deuteron_wignerArg[nTrials][20];
    Double_t weight_deuteron_wignerArg[nTrials][20];
    Int_t nDeuterons_wignerArg[nTrials];

    // Save Deuteron Properties (Wigner Argonne18)
    Double_t pt_deuteron_TruewignerGaus[nTrials][20];
    Double_t phi_deuteron_TruewignerGaus[nTrials][20];
    Double_t weight_deuteron_TruewignerGaus[nTrials][20];
    Int_t nDeuterons_TruewignerGaus[nTrials];

    // Save Deuteron Properties (Wigner Argonne18)
    Double_t pt_deuteron_TruewignerDoubleGaus[nTrials][20];
    Double_t phi_deuteron_TruewignerDoubleGaus[nTrials][20];
    Double_t weight_deuteron_TruewignerDoubleGaus[nTrials][20];
    Int_t nDeuterons_TruewignerDoubleGaus[nTrials];

    // Initialize Deuteron Counters
    for (Int_t iTrial = 0; iTrial < nTrials; iTrial++)
    {

        nDeuterons_simpleCoal[iTrial] = 0;
        nDeuterons_wignerGaus[iTrial] = 0;
        nDeuterons_wignerArg[iTrial] = 0;
        nDeuterons_TruewignerGaus[iTrial] = 0;
        nDeuterons_TruewignerDoubleGaus[iTrial] = 0;

    }
    //************************************************** SIMPLE COALESCENCE **************************************************//
    for (Int_t iTrial = 0; iTrial < nTrials; iTrial++)
    {

        // Reset Neutron Status
        for (Int_t in = 0; in < (Int_t)neutron_ID.size(); in++)
        {
            neutron_status[in] = 0;
        }

        for (Int_t ip = 0; ip < (Int_t)proton_ID.size(); ip++)
        {

            // Proton 4-Momentum
            AliMCParticle *proton = (AliMCParticle *)fMCEvent->GetTrack(proton_ID[ip]);
            TLorentzVector p_proton;
            p_proton.SetXYZM(proton->Px(), proton->Py(), proton->Pz(), mp);

            for (Int_t in = 0; in < (Int_t)neutron_ID.size(); in++)
            {

                // Neutron 4-Momentum
                AliMCParticle *neutron = (AliMCParticle *)fMCEvent->GetTrack(neutron_ID[in]);
                TLorentzVector p_neutron;
                p_neutron.SetXYZM(neutron->Px(), neutron->Py(), neutron->Pz(), mn);

                // Deuteron 4-Momentum
                TLorentzVector p_deuteron;
                p_deuteron.SetXYZM(proton->Px() + neutron->Px(), proton->Py() + neutron->Py(), proton->Pz() + neutron->Pz(), md);
                Double_t beta_x = p_deuteron.Px() / p_deuteron.E();
                Double_t beta_y = p_deuteron.Py() / p_deuteron.E();
                Double_t beta_z = p_deuteron.Pz() / p_deuteron.E();
                TVector3 beta(beta_x, beta_y, beta_z);

                // Lorentz Transformations (from Lab to Deuteron Frame)
                TLorentzVector p_proton_prime = LorentzTransform(p_proton, beta);
                TLorentzVector p_neutron_prime = LorentzTransform(p_neutron, beta);
                Double_t deltaP = (p_proton_prime - p_neutron_prime).P();
                Double_t deutWeight = GetDeuteronWeight(p_proton.Pt(), p_neutron.Pt());

                // Fill DeltaP Distribution
                if (iTrial == 0)
                    hDeltaP->Fill(deltaP);

                if (neutron_status[in] == 1)
                    continue; // Skip already used neutrons

                // Simple Coalescence Condition
                if (deltaP < p0[iTrial])
                {

                    neutron_status[in] = 1;
                    Double_t y = p_deuteron.Rapidity();
                    if (TMath::Abs(y) < 0.5)
                    {

                        hDeuteronsINELgtZERO_simpleCoal[iTrial]->Fill(p_deuteron.Pt(), deutWeight);

                        // Store Deuteron Properties
                        phi_deuteron_simpleCoal[iTrial][nDeuterons_simpleCoal[iTrial]] = TVector2::Phi_0_2pi(p_deuteron.Phi());
                        pt_deuteron_simpleCoal[iTrial][nDeuterons_simpleCoal[iTrial]] = p_deuteron.Pt();
                        weight_deuteron_simpleCoal[iTrial][nDeuterons_simpleCoal[iTrial]] = deutWeight;
                        nDeuterons_simpleCoal[iTrial]++;

                        // Rapidity Distributions of Protons and Neutrons
                        hRapidityProtons->Fill(proton->Y());
                        hRapidityNeutrons->Fill(neutron->Y());

                        // DeltaPhi
                        TVector3 mom_proton = p_proton.Vect();
                        TVector3 mom_neutron = p_neutron.Vect();
                        Double_t deltaPhi = (180.0 / TMath::Pi()) * mom_proton.Angle(mom_neutron);
                        if (iTrial == 13)
                            hDeltaPhi_INELgtZERO->Fill(p_deuteron.Pt(), deltaPhi);

                        // QA Histogram
                        if (p_deuteron.Pt() > 0.7 && p_deuteron.Pt() < 0.8)
                            hPtProtonsFirstBinDeut->Fill(proton->Pt());
                    }
                    break;
                }
            }
        }
    }
    //************************************************** COALESCENCE: WIGNER GAUS **************************************************//
    for (Int_t iTrial = 0; iTrial < nTrials; iTrial++)
    {
        // std::cout << "Trial " << iTrial << std::endl;

        // Reset Neutron Status
        for (Int_t in = 0; in < (Int_t)neutron_ID.size(); in++)
        {
            neutron_status[in] = 0;
        }

        for (Int_t ip = 0; ip < (Int_t)proton_ID.size(); ip++)
        {
            // std::cout << "Proton " << ip  << std::endl;
            // Proton 4-Momentum
            AliMCParticle *proton = (AliMCParticle *)fMCEvent->GetTrack(proton_ID[ip]);
            TLorentzVector p_proton;
            p_proton.SetXYZM(proton->Px(), proton->Py(), proton->Pz(), mp);

            for (Int_t in = 0; in < (Int_t)neutron_ID.size(); in++)
            {
                // std::cout << "HOLY CANOLY" << std::endl;
                //  Neutron 4-Momentum
                AliMCParticle *neutron = (AliMCParticle *)fMCEvent->GetTrack(neutron_ID[in]);
                TLorentzVector p_neutron;
                p_neutron.SetXYZM(neutron->Px(), neutron->Py(), neutron->Pz(), mn);

                // Deuteron 4-Momentum
                TLorentzVector p_deuteron;
                p_deuteron.SetXYZM(proton->Px() + neutron->Px(), proton->Py() + neutron->Py(), proton->Pz() + neutron->Pz(), md);
                Double_t beta_x = p_deuteron.Px() / p_deuteron.E();
                Double_t beta_y = p_deuteron.Py() / p_deuteron.E();
                Double_t beta_z = p_deuteron.Pz() / p_deuteron.E();
                TVector3 beta(beta_x, beta_y, beta_z);

                // Lorentz Transformation
                TLorentzVector p_proton_prime = LorentzTransform(p_proton, beta);
                TLorentzVector p_neutron_prime = LorentzTransform(p_neutron, beta);

                Double_t deltaX = GetSpatialDistance(p_proton, p_neutron, beta);
                Double_t deltaP = (p_proton_prime - p_neutron_prime).P();
                Double_t deutWeight = GetDeuteronWeight(p_proton.Pt(), p_neutron.Pt());

                if (neutron_status[in] == 1)
                    continue; // Skip already used neutrons

                // Coalescence Condition
                if (DoCoalescence(deltaX, deltaP, p0[iTrial], "Gaus"))
                {

                    neutron_status[in] = 1;
                    Double_t y = p_deuteron.Rapidity();
                    if (TMath::Abs(y) < 0.5)
                    {

                        hDeuteronsINELgtZERO_wignerGaus[iTrial]->Fill(p_deuteron.Pt(), deutWeight);

                        // Store Deuteron Properties
                        phi_deuteron_wignerGaus[iTrial][nDeuterons_wignerGaus[iTrial]] = TVector2::Phi_0_2pi(p_deuteron.Phi());
                        pt_deuteron_wignerGaus[iTrial][nDeuterons_wignerGaus[iTrial]] = p_deuteron.Pt();
                        weight_deuteron_wignerGaus[iTrial][nDeuterons_wignerGaus[iTrial]] = deutWeight;
                        nDeuterons_wignerGaus[iTrial]++;
                    }
                    break;
                }
            }
        }
    }
    //******************************************************COALESCENCE: TRUE WIGNER GAUS***************************************************//

    // std::cout << "Trial " << iTrial << std::endl;

    // Reset Neutron Status
    for (Int_t in = 0; in < (Int_t)neutron_ID.size(); in++)
    {
        neutron_status[in] = 0;
    }

    for (Int_t ip = 0; ip < (Int_t)proton_ID.size(); ip++)
    {
        // std::cout << "Proton " << ip  << std::endl;
        // Proton 4-Momentum
        AliMCParticle *proton = (AliMCParticle *)fMCEvent->GetTrack(proton_ID[ip]);
        TLorentzVector p_proton;
        p_proton.SetXYZM(proton->Px(), proton->Py(), proton->Pz(), mp);

        for (Int_t in = 0; in < (Int_t)neutron_ID.size(); in++)
        {
            // std::cout << "HOLY CANOLY" << std::endl;
            //  Neutron 4-Momentum
            AliMCParticle *neutron = (AliMCParticle *)fMCEvent->GetTrack(neutron_ID[in]);
            TLorentzVector p_neutron;
            p_neutron.SetXYZM(neutron->Px(), neutron->Py(), neutron->Pz(), mn);

            // Deuteron 4-Momentum
            TLorentzVector p_deuteron;
            p_deuteron.SetXYZM(proton->Px() + neutron->Px(), proton->Py() + neutron->Py(), proton->Pz() + neutron->Pz(), md);
            Double_t beta_x = p_deuteron.Px() / p_deuteron.E();
            Double_t beta_y = p_deuteron.Py() / p_deuteron.E();
            Double_t beta_z = p_deuteron.Pz() / p_deuteron.E();
            TVector3 beta(beta_x, beta_y, beta_z);

            // Lorentz Transformation
            TLorentzVector p_proton_prime = LorentzTransform(p_proton, beta);
            TLorentzVector p_neutron_prime = LorentzTransform(p_neutron, beta);

            Double_t deltaX = GetSpatialDistance(p_proton, p_neutron, beta);
            Double_t deltaP = (p_proton_prime - p_neutron_prime).P();
            Double_t deutWeight = GetDeuteronWeight(p_proton.Pt(), p_neutron.Pt());

            if (neutron_status[in] == 1)
                continue; // Skip already used neutrons

            //  True Wigner function Coalescence
            Double_t rndmG = gRandom->Uniform(0.0, 1.0); // random number for Gaus Wavefunction
            // Gaussian Wigner function Coalescence Condition
            Double_t SourceSize = 1.2; // fm. Bambi This should be done as a function of mT
            Double_t zeta = pow(pow(3.2, 2) / (pow(3.2, 2) + 4 * SourceSize), 3 / 2);
            if (rndmG < 3 * zeta * exp(-deltaP / 2 * deltaP / 2 * 3.2 * 3.2 * 5.08 * 5.08))
            {

                neutron_status[in] = 1;
                Double_t y = p_deuteron.Rapidity();
                if (TMath::Abs(y) < 0.5)
                {

                    hDeuteronsINELgtZERO_TruewignerGaus->Fill(p_deuteron.Pt(), deutWeight);

                    // Store Deuteron Properties
                    /*phi_deuteron_TruewignerGaus[iTrial][nDeuterons_TruewignerGaus[iTrial]] = TVector2::Phi_0_2pi(p_deuteron.Phi());
                    pt_deuteron_TruewignerGaus[iTrial][nDeuterons_TruewignerGaus[iTrial]] = p_deuteron.Pt();
                    weight_deuteron_TruewignerGaus[iTrial][nDeuterons_TruewignerGaus[iTrial]] = deutWeight;
                    nDeuterons_TruewignerGaus[iTrial]++;*/
                }
                break;
            }
        }
    }

    //******************************************************COALESCENCE: TRUE WIGNER DOUBLE GAUS***************************************************//

    // Reset Neutron Status
    for (Int_t in = 0; in < (Int_t)neutron_ID.size(); in++)
    {
        neutron_status[in] = 0;
    }

    for (Int_t ip = 0; ip < (Int_t)proton_ID.size(); ip++)
    {
        // std::cout << "Proton " << ip  << std::endl;
        // Proton 4-Momentum
        AliMCParticle *proton = (AliMCParticle *)fMCEvent->GetTrack(proton_ID[ip]);
        TLorentzVector p_proton;
        p_proton.SetXYZM(proton->Px(), proton->Py(), proton->Pz(), mp);

        for (Int_t in = 0; in < (Int_t)neutron_ID.size(); in++)
        {
            // std::cout << "HOLY CANOLY" << std::endl;
            //  Neutron 4-Momentum
            AliMCParticle *neutron = (AliMCParticle *)fMCEvent->GetTrack(neutron_ID[in]);
            TLorentzVector p_neutron;
            p_neutron.SetXYZM(neutron->Px(), neutron->Py(), neutron->Pz(), mn);

            // Deuteron 4-Momentum
            TLorentzVector p_deuteron;
            p_deuteron.SetXYZM(proton->Px() + neutron->Px(), proton->Py() + neutron->Py(), proton->Pz() + neutron->Pz(), md);
            Double_t beta_x = p_deuteron.Px() / p_deuteron.E();
            Double_t beta_y = p_deuteron.Py() / p_deuteron.E();
            Double_t beta_z = p_deuteron.Pz() / p_deuteron.E();
            TVector3 beta(beta_x, beta_y, beta_z);

            // Lorentz Transformation
            TLorentzVector p_proton_prime = LorentzTransform(p_proton, beta);
            TLorentzVector p_neutron_prime = LorentzTransform(p_neutron, beta);

            Double_t deltaX = GetSpatialDistance(p_proton, p_neutron, beta);
            Double_t deltaP = (p_proton_prime - p_neutron_prime).P();
            Double_t deutWeight = GetDeuteronWeight(p_proton.Pt(), p_neutron.Pt());

            if (neutron_status[in] == 1)
                continue; // Skip already used neutrons

            //  True Wigner function Coalescence
            Double_t rndmDG = gRandom->Uniform(0.0, 1.0); // random number for Gaus Wavefunction
            // Gaussian Wigner function Coalescence Condition
            Double_t SourceSize = 1.2; // fm. Bambi This should be done as a function of mT
            Double_t zeta1 = pow((pow(3.979, 2)) / (pow(3.979, 2) + 4 * SourceSize), 3 / 2);
            Double_t zeta2 = pow((pow(0.89, 2)) / (pow(0.89, 2) + 4 * SourceSize), 3 / 2);
            Double_t Delta = 0.581;
            if (rndmDG < 3 * (Delta * zeta1 * exp(-deltaP / 2 * deltaP / 2 * 3.979 * 5.08 * 3.979 * 5.08) + (1 - Delta) * zeta2 * exp(-deltaP / 2 * deltaP / 2 * 0.89 * 5.08 * 0.89 * 5.08)))
            {

                neutron_status[in] = 1;
                Double_t y = p_deuteron.Rapidity();
                if (TMath::Abs(y) < 0.5)
                {

                    hDeuteronsINELgtZERO_TruewignerDoubleGaus->Fill(p_deuteron.Pt(), deutWeight);

                    // Store Deuteron Properties
                    /*phi_deuteron_TruewignerDoubleGaus[iTrial][nDeuterons_TruewignerDoubleGaus[iTrial]] = TVector2::Phi_0_2pi(p_deuteron.Phi());
                    pt_deuteron_TruewignerDoubleGaus[iTrial][nDeuterons_TruewignerDoubleGaus[iTrial]] = p_deuteron.Pt();
                    weight_deuteron_TruewignerDoubleGaus[iTrial][nDeuterons_TruewignerDoubleGaus[iTrial]] = deutWeight;
                    nDeuterons_TruewignerDoubleGaus[iTrial]++;*/
                }
                break;
            }
        }
    }

    //************************************************** COALESCENCE: WIGNER ARGONNE **************************************************//
    for (Int_t iTrial = 0; iTrial < nTrials; iTrial++)
    {

        // Reset Neutron Status
        for (Int_t in = 0; in < (Int_t)neutron_ID.size(); in++)
        {
            neutron_status[in] = 0;
        }

        for (Int_t ip = 0; ip < (Int_t)proton_ID.size(); ip++)
        {

            // Proton 4-Momentum
            AliMCParticle *proton = (AliMCParticle *)fMCEvent->GetTrack(proton_ID[ip]);
            TLorentzVector p_proton;
            p_proton.SetXYZM(proton->Px(), proton->Py(), proton->Pz(), mp);

            for (Int_t in = 0; in < (Int_t)neutron_ID.size(); in++)
            {

                // Neutron 4-Momentum
                AliMCParticle *neutron = (AliMCParticle *)fMCEvent->GetTrack(neutron_ID[in]);
                TLorentzVector p_neutron;
                p_neutron.SetXYZM(neutron->Px(), neutron->Py(), neutron->Pz(), mn);

                // Deuteron 4-Momentum
                TLorentzVector p_deuteron;
                p_deuteron.SetXYZM(proton->Px() + neutron->Px(), proton->Py() + neutron->Py(), proton->Pz() + neutron->Pz(), md);
                Double_t beta_x = p_deuteron.Px() / p_deuteron.E();
                Double_t beta_y = p_deuteron.Py() / p_deuteron.E();
                Double_t beta_z = p_deuteron.Pz() / p_deuteron.E();
                TVector3 beta(beta_x, beta_y, beta_z);

                // Lorentz Transformation
                TLorentzVector p_proton_prime = LorentzTransform(p_proton, beta);
                TLorentzVector p_neutron_prime = LorentzTransform(p_neutron, beta);

                Double_t deltaX = GetSpatialDistance(p_proton, p_neutron, beta);
                Double_t deltaP = (p_proton_prime - p_neutron_prime).P();
                Double_t deutWeight = GetDeuteronWeight(p_proton.Pt(), p_neutron.Pt());

                if (neutron_status[in] == 1)
                    continue; // Skip already used neutrons

                // Coalescence Condition
                if (DoCoalescence(deltaX, deltaP, p0[iTrial], "Arg"))
                {

                    neutron_status[in] = 1;
                    Double_t y = p_deuteron.Rapidity();
                    if (TMath::Abs(y) < 0.5)
                    {

                        hDeuteronsINELgtZERO_wignerArg[iTrial]->Fill(p_deuteron.Pt(), deutWeight);

                        // Store Deuteron Properties
                        phi_deuteron_wignerArg[iTrial][nDeuterons_wignerArg[iTrial]] = TVector2::Phi_0_2pi(p_deuteron.Phi());
                        pt_deuteron_wignerArg[iTrial][nDeuterons_wignerArg[iTrial]] = p_deuteron.Pt();
                        weight_deuteron_wignerArg[iTrial][nDeuterons_wignerArg[iTrial]] = deutWeight;
                        nDeuterons_wignerArg[iTrial]++;
                    }
                    break;
                }
            }
        }
    }
    //*********************************************************************************************************************************//

    // Selection of Leading Particle
    Int_t lp = GetLeadingParticle();
    AliMCParticle *leading_particle = (AliMCParticle *)fMCEvent->GetTrack(lp);
    Double_t pt_leading = leading_particle->Pt();
    Double_t phi_leading = TVector2::Phi_0_2pi(leading_particle->Phi());
    if (pt_leading < 5.0)
        return;
    hNumberOfEvents->Fill(17.5);

    // Transverse Multiplicity
    Int_t nParticles_Transv = GetTransverseMultiplicity(lp);
    hTransverseMult->Fill(nParticles_Transv);

    // R_{T}
    Double_t Rt = static_cast<Double_t>(nParticles_Transv) / fAverage_Nch_Transv;
    hRtDistribution->Fill(Rt);

    // Fill Proton Spectra in Azimuthal Regions
    vector<Int_t> proton_region;

    for (Int_t i = 0; i < (Int_t)proton_ID.size(); i++)
    {

        // Get Particle
        AliMCParticle *particle = (AliMCParticle *)fMCEvent->GetTrack(proton_ID[i]);

        // Variables
        Double_t pt = particle->Pt();
        Double_t phi_particle = TVector2::Phi_0_2pi(particle->Phi());

        // Weights
        Double_t wp_toward = GetProtonWeight(pt, 0);
        Double_t wp_transv = GetProtonWeight(pt, 1);
        Double_t wp_away = GetProtonWeight(pt, 2);

        // Save Information on Azimuthal Region
        if (IsParticleInTowardRegion(phi_particle, phi_leading))
            proton_region.push_back(0);
        if (IsParticleInTransverseRegion(phi_particle, phi_leading))
            proton_region.push_back(1);
        if (IsParticleInAwayRegion(phi_particle, phi_leading))
            proton_region.push_back(2);

        // Fill p_{T} Spectra of Protons
        if (TMath::Abs(particle->Y()) > 0.5)
            continue;
        if (IsParticleInTowardRegion(phi_particle, phi_leading))
        {
            hProtons_Toward->Fill(pt);
            hProtons_Toward_reshaped->Fill(pt, wp_toward);
        }
        if (IsParticleInTransverseRegion(phi_particle, phi_leading))
        {
            hProtons_Transv->Fill(pt);
            hProtons_Transv_reshaped->Fill(pt, wp_transv);
        }
        if (IsParticleInAwayRegion(phi_particle, phi_leading))
        {
            hProtons_Away->Fill(pt);
            hProtons_Away_reshaped->Fill(pt, wp_away);
        }
    }

    // Fill Neutron Spectra in Azimuthal Regions
    vector<Int_t> neutron_region;

    for (Int_t i = 0; i < (Int_t)neutron_ID.size(); i++)
    {

        // Get Particle
        AliMCParticle *particle = (AliMCParticle *)fMCEvent->GetTrack(neutron_ID[i]);

        // Variables
        Double_t pt = particle->Pt();
        Double_t phi_particle = TVector2::Phi_0_2pi(particle->Phi());

        // Weights
        Double_t wn_toward = GetProtonWeight(pt, 0);
        Double_t wn_transv = GetProtonWeight(pt, 1);
        Double_t wn_away = GetProtonWeight(pt, 2);

        // Save Information on Azimuthal Region
        if (IsParticleInTowardRegion(phi_particle, phi_leading))
            neutron_region.push_back(0);
        if (IsParticleInTransverseRegion(phi_particle, phi_leading))
            neutron_region.push_back(1);
        if (IsParticleInAwayRegion(phi_particle, phi_leading))
            neutron_region.push_back(2);

        // Fill p_{T} Spectra of Neutrons
        if (TMath::Abs(particle->Y()) > 0.5)
            continue;
        if (IsParticleInTowardRegion(phi_particle, phi_leading))
        {
            hNeutrons_Toward->Fill(pt);
            hNeutrons_Toward_reshaped->Fill(pt, wn_toward);
        }
        if (IsParticleInTransverseRegion(phi_particle, phi_leading))
        {
            hNeutrons_Transv->Fill(pt);
            hNeutrons_Transv_reshaped->Fill(pt, wn_transv);
        }
        if (IsParticleInAwayRegion(phi_particle, phi_leading))
        {
            hNeutrons_Away->Fill(pt);
            hNeutrons_Away_reshaped->Fill(pt, wn_away);
        }
    }

    // Simple Coalescence
    for (Int_t iTrial = 0; iTrial < nTrials; iTrial++)
    {

        // Reset Neutron Status
        for (Int_t in = 0; in < (Int_t)neutron_ID.size(); in++)
        {
            neutron_status[in] = 0;
        }

        for (Int_t ip = 0; ip < (Int_t)proton_ID.size(); ip++)
        {

            // Proton 4-Momentum
            AliMCParticle *proton = (AliMCParticle *)fMCEvent->GetTrack(proton_ID[ip]);
            TLorentzVector p_proton;
            p_proton.SetXYZM(proton->Px(), proton->Py(), proton->Pz(), mp);

            for (Int_t in = 0; in < (Int_t)neutron_ID.size(); in++)
            {

                if (neutron_status[in] == 1)
                    continue; // Skip already used neutrons

                // Neutron 4-Momentum
                AliMCParticle *neutron = (AliMCParticle *)fMCEvent->GetTrack(neutron_ID[in]);
                TLorentzVector p_neutron;
                p_neutron.SetXYZM(neutron->Px(), neutron->Py(), neutron->Pz(), mn);

                // Deuteron 4-Momentum
                TLorentzVector p_deuteron;
                p_deuteron.SetXYZM(proton->Px() + neutron->Px(), proton->Py() + neutron->Py(), proton->Pz() + neutron->Pz(), md);
                Double_t beta_x = p_deuteron.Px() / p_deuteron.E();
                Double_t beta_y = p_deuteron.Py() / p_deuteron.E();
                Double_t beta_z = p_deuteron.Pz() / p_deuteron.E();
                TVector3 beta(beta_x, beta_y, beta_z);

                // Lorentz Transformations (from Lab to Deuteron Frame)
                TLorentzVector p_proton_prime = LorentzTransform(p_proton, beta);
                TLorentzVector p_neutron_prime = LorentzTransform(p_neutron, beta);
                Double_t deltaP = (p_proton_prime - p_neutron_prime).P();
                Double_t deutWeight = GetDeuteronWeight(p_proton.Pt(), p_neutron.Pt(), proton_region[ip], neutron_region[in]);

                // Simple Coalescence Condition
                if (deltaP < p0[iTrial])
                {

                    neutron_status[in] = 1;
                    Double_t y = p_deuteron.Rapidity();
                    if (TMath::Abs(y) < 0.5)
                    {

                        Double_t pt = p_deuteron.Pt();
                        Double_t phi_particle = TVector2::Phi_0_2pi(p_deuteron.Phi());

                        // DeltaPhi
                        TVector3 mom_proton = p_proton.Vect();
                        TVector3 mom_neutron = p_neutron.Vect();
                        Double_t deltaPhi = (180.0 / TMath::Pi()) * mom_proton.Angle(mom_neutron);

                        if (IsParticleInTowardRegion(phi_particle, phi_leading))
                            hDeuterons_Toward_simpleCoal[iTrial]->Fill(pt, deutWeight);
                        if (IsParticleInTransverseRegion(phi_particle, phi_leading))
                            hDeuterons_Transv_simpleCoal[iTrial]->Fill(pt, deutWeight);

                        // Fill Angular Correlations
                        if (iTrial == 13)
                        {

                            if (IsParticleInTowardRegion(phi_particle, phi_leading))
                                hDeltaPhi_Toward->Fill(pt, deltaPhi);
                            if (IsParticleInTransverseRegion(phi_particle, phi_leading))
                                hDeltaPhi_Transv->Fill(pt, deltaPhi);
                            if (IsParticleInAwayRegion(phi_particle, phi_leading))
                                hDeltaPhi_Away->Fill(pt, deltaPhi);
                        }
                    }
                    break;
                }
            }
        }
    }
    for (Int_t iTrial = 0; iTrial < nTrials; iTrial++)
    {
        for (Int_t ip = 0; ip < (Int_t)proton_ID.size(); ip++)
        {

            // Proton 4-Momentum
            AliMCParticle *proton = (AliMCParticle *)fMCEvent->GetTrack(proton_ID[ip]);
            TLorentzVector p_proton;
            p_proton.SetXYZM(proton->Px(), proton->Py(), proton->Pz(), mp);

            for (Int_t in = 0; in < (Int_t)neutron_ID.size(); in++)
            {

                if (neutron_status[in] == 1)
                    continue; // Skip already used neutrons

                // Neutron 4-Momentum
                AliMCParticle *neutron = (AliMCParticle *)fMCEvent->GetTrack(neutron_ID[in]);
                TLorentzVector p_neutron;
                p_neutron.SetXYZM(neutron->Px(), neutron->Py(), neutron->Pz(), mn);

                // Deuteron 4-Momentum
                TLorentzVector p_deuteron;
                p_deuteron.SetXYZM(proton->Px() + neutron->Px(), proton->Py() + neutron->Py(), proton->Pz() + neutron->Pz(), md);
                Double_t beta_x = p_deuteron.Px() / p_deuteron.E();
                Double_t beta_y = p_deuteron.Py() / p_deuteron.E();
                Double_t beta_z = p_deuteron.Pz() / p_deuteron.E();
                TVector3 beta(beta_x, beta_y, beta_z);

                // Lorentz Transformations (from Lab to Deuteron Frame)
                TLorentzVector p_proton_prime = LorentzTransform(p_proton, beta);
                TLorentzVector p_neutron_prime = LorentzTransform(p_neutron, beta);
                Double_t deltaP = (p_proton_prime - p_neutron_prime).P();
                Double_t deutWeight = GetDeuteronWeight(p_proton.Pt(), p_neutron.Pt(), proton_region[ip], neutron_region[in]);

                // Simple Coalescence Condition
                Double_t rndmG = gRandom->Uniform(0.0, 1.0); // random number for Gaus Wavefunction
                // Gaussian Wigner function Coalescence Condition
                Double_t SourceSize = iTrial * 0.05; // fm. Bambi This should be done as a function of mT
                Double_t zeta = pow((pow(3.2, 2)) / (pow(3.2, 2) + 4 * SourceSize), 3 / 2);
                if (rndmG < 3 * zeta * exp(-deltaP / 2 * deltaP / 2 * 3.2 * 5.08 * 3.2 * 5.08))
                {

                    neutron_status[in] = 1;
                    Double_t y = p_deuteron.Rapidity();
                    if (TMath::Abs(y) < 0.5)
                    {
                        Double_t pt = p_deuteron.Pt();
                        Double_t phi_particle = TVector2::Phi_0_2pi(p_deuteron.Phi());

                        // DeltaPhi
                        TVector3 mom_proton = p_proton.Vect();
                        TVector3 mom_neutron = p_neutron.Vect();
                        Double_t deltaPhi = (180.0 / TMath::Pi()) * mom_proton.Angle(mom_neutron);

                        if (IsParticleInTowardRegion(phi_particle, phi_leading))
                            hDeuterons_Toward_TruewignerGaus[iTrial]->Fill(pt, deutWeight);
                        if (IsParticleInTransverseRegion(phi_particle, phi_leading))
                            hDeuterons_Transv_TruewignerGaus[iTrial]->Fill(pt, deutWeight);
                    }
                    break;
                }
            }
        }
    }
    for (Int_t iTrial = 0; iTrial < nTrials; iTrial++)
    {
        for (Int_t ip = 0; ip < (Int_t)proton_ID.size(); ip++)
        {

            // Proton 4-Momentum
            AliMCParticle *proton = (AliMCParticle *)fMCEvent->GetTrack(proton_ID[ip]);
            TLorentzVector p_proton;
            p_proton.SetXYZM(proton->Px(), proton->Py(), proton->Pz(), mp);

            for (Int_t in = 0; in < (Int_t)neutron_ID.size(); in++)
            {

                if (neutron_status[in] == 1)
                    continue; // Skip already used neutrons

                // Neutron 4-Momentum
                AliMCParticle *neutron = (AliMCParticle *)fMCEvent->GetTrack(neutron_ID[in]);
                TLorentzVector p_neutron;
                p_neutron.SetXYZM(neutron->Px(), neutron->Py(), neutron->Pz(), mn);

                // Deuteron 4-Momentum
                TLorentzVector p_deuteron;
                p_deuteron.SetXYZM(proton->Px() + neutron->Px(), proton->Py() + neutron->Py(), proton->Pz() + neutron->Pz(), md);
                Double_t beta_x = p_deuteron.Px() / p_deuteron.E();
                Double_t beta_y = p_deuteron.Py() / p_deuteron.E();
                Double_t beta_z = p_deuteron.Pz() / p_deuteron.E();
                TVector3 beta(beta_x, beta_y, beta_z);

                // Lorentz Transformations (from Lab to Deuteron Frame)
                TLorentzVector p_proton_prime = LorentzTransform(p_proton, beta);
                TLorentzVector p_neutron_prime = LorentzTransform(p_neutron, beta);
                Double_t deltaP = (p_proton_prime - p_neutron_prime).P();
                Double_t deutWeight = GetDeuteronWeight(p_proton.Pt(), p_neutron.Pt(), proton_region[ip], neutron_region[in]);

                // Simple Coalescence Condition
                // Gaussian Wigner function Coalescence Condition
                Double_t SourceSize = iTrial * 0.05;          // fm. Bambi This should be done as a function of mT
                Double_t rndmDG = gRandom->Uniform(0.0, 1.0); // random number for double Gaus Wigner function
                Double_t zeta1 = pow((pow(3.979, 2)) / (pow(3.979, 2) + 4 * SourceSize), 3 / 2);
                Double_t zeta2 = pow((pow(0.89, 2)) / (pow(0.89, 2) + 4 * SourceSize), 3 / 2);
                Double_t Delta = 0.581;
                if (rndmDG < 3 * (Delta * zeta1 * exp(-deltaP / 2 * deltaP / 2 * 3.979 * 5.08 * 3.979 * 5.08) + (1 - Delta) * zeta2 * exp(-deltaP / 2 * deltaP / 2 * 0.89 * 5.08 * 0.89 * 5.08)))
                {
                    neutron_status[in] = 1;
                    Double_t y = p_deuteron.Rapidity();
                    if (TMath::Abs(y) < 0.5)
                    {
                        Double_t pt = p_deuteron.Pt();
                        Double_t phi_particle = TVector2::Phi_0_2pi(p_deuteron.Phi());

                        // DeltaPhi
                        TVector3 mom_proton = p_proton.Vect();
                        TVector3 mom_neutron = p_neutron.Vect();
                        Double_t deltaPhi = (180.0 / TMath::Pi()) * mom_proton.Angle(mom_neutron);

                        if (IsParticleInTowardRegion(phi_particle, phi_leading))
                            hDeuterons_Toward_TruewignerDoubleGaus[iTrial]->Fill(pt, deutWeight);
                        if (IsParticleInTransverseRegion(phi_particle, phi_leading))
                            hDeuterons_Transv_TruewignerDoubleGaus[iTrial]->Fill(pt, deutWeight);
                    }
                    break;
                }
            }
        }
    }

    /*

     //Fill Deuteron Spectra in Azimuthal Regions: Simple Coalescence
     for (Int_t iTrial=0 ; iTrial<nTrials ; iTrial++)  {
         for (Int_t i=0 ; i<nDeuterons_simpleCoal[iTrial] ; i++)  {

             //Variables
             Double_t pt           = pt_deuteron_simpleCoal[iTrial][i];
             Double_t phi_particle = phi_deuteron_simpleCoal[iTrial][i];
             Double_t deutWeight   = weight_deuteron_simpleCoal[iTrial][i];

             //Fill p_{T} Spectra Deuterons
             if (IsParticleInTowardRegion(phi_particle,phi_leading))     hDeuterons_Toward_simpleCoal[iTrial]->Fill(pt,deutWeight);
             if (IsParticleInTransverseRegion(phi_particle,phi_leading)) hDeuterons_Transv_simpleCoal[iTrial]->Fill(pt,deutWeight);
         }
     }


    //Fill Deuteron Spectra in Azimuthal Regions: Wigner Gaus
    for (Int_t iTrial=0 ; iTrial<nTrials ; iTrial++)  {
        for (Int_t i=0 ; i<nDeuterons_wignerGaus[iTrial] ; i++)  {

            //Variables
            Double_t pt           = pt_deuteron_wignerGaus[iTrial][i];
            Double_t phi_particle = phi_deuteron_wignerGaus[iTrial][i];
            Double_t deutWeight   = weight_deuteron_wignerGaus[iTrial][i];

            //Fill p_{T} Spectra Deuterons
            if (IsParticleInTowardRegion(phi_particle,phi_leading))     hDeuterons_Toward_wignerGaus[iTrial]->Fill(pt,deutWeight);
            if (IsParticleInTransverseRegion(phi_particle,phi_leading)) hDeuterons_Transv_wignerGaus[iTrial]->Fill(pt,deutWeight);
        }
    }

    //Fill Deuteron Spectra in Azimuthal Regions: Wigner Argonne
    for (Int_t iTrial=0 ; iTrial<nTrials ; iTrial++)  {
        for (Int_t i=0 ; i<nDeuterons_wignerArg[iTrial] ; i++)  {

            //Variables
            Double_t pt           = pt_deuteron_wignerArg[iTrial][i];
            Double_t phi_particle = phi_deuteron_wignerArg[iTrial][i];
            Double_t deutWeight   = weight_deuteron_wignerArg[iTrial][i];

            //Fill p_{T} Spectra Deuterons
            if (IsParticleInTowardRegion(phi_particle,phi_leading))     hDeuterons_Toward_wignerArg[iTrial]->Fill(pt,deutWeight);
            if (IsParticleInTransverseRegion(phi_particle,phi_leading)) hDeuterons_Transv_wignerArg[iTrial]->Fill(pt,deutWeight);
        }
    }
    */

    // Post Output Data
    PostData(1, fOutputList);
}
//____________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskPythiaCoalescence::GetEvent()
{

    // Get AOD Event
    fAODevent = dynamic_cast<AliAODEvent *>(InputEvent());
    if (!fAODevent)
        return kFALSE;
    hNumberOfEvents->Fill(0.5);

    // Get MC Event
    fMCEvent = MCEvent();
    if (!fMCEvent)
        return kFALSE;
    hNumberOfEvents->Fill(1.5);
    hNumberOfEvents->Fill(2.5);

    // Selection of INEL>0 Events
    if (!IsINELgtZERO())
        return kFALSE;
    hNumberOfEvents->Fill(3.5);

    // No Vertex & Multiplicity Selection
    hNumberOfEvents->Fill(4.5);
    hNumberOfEvents->Fill(5.5);
    hNumberOfEvents->Fill(6.5);
    hNumberOfEvents->Fill(7.5);
    hNumberOfEvents->Fill(8.5);
    hNumberOfEvents->Fill(9.5);
    hNumberOfEvents->Fill(10.5);
    hNumberOfEvents->Fill(11.5);
    hNumberOfEvents->Fill(12.5);
    hNumberOfEvents->Fill(13.5);
    hNumberOfEvents->Fill(14.5);
    hNumberOfEvents->Fill(15.5);
    hNumberOfEvents->Fill(16.5);

    return kTRUE;
}
//____________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskPythiaCoalescence::IsINELgtZERO()
{

    // Initialization
    Bool_t isEventINELgtZERO = (kFALSE);

    // Loop over Generated Particles
    for (Int_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++)
    {

        AliMCParticle *particle = (AliMCParticle *)fMCEvent->GetTrack(i);
        if (!particle)
            continue;
        if (particle->Charge() == 0)
            continue;
        if (IsInjectedParticle(particle))
            continue;
        if (particle->MCStatusCode() <= 0)
            continue;

        Bool_t isFSParticle = (kFALSE);
        if (TMath::Abs(particle->PdgCode()) == 211)
            isFSParticle = kTRUE;
        if (TMath::Abs(particle->PdgCode()) == 321)
            isFSParticle = kTRUE;
        if (TMath::Abs(particle->PdgCode()) == 2212)
            isFSParticle = kTRUE;
        if (TMath::Abs(particle->PdgCode()) == 11)
            isFSParticle = kTRUE;
        if (isFSParticle && TMath::Abs(particle->Eta()) < 1.0)
        {
            isEventINELgtZERO = kTRUE;
            break;
        }
    }

    return isEventINELgtZERO;
}
//____________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskPythiaCoalescence::GetProtonWeight(Double_t pt)
{

    // Initialization
    Double_t w(1);

    /*
    Int_t ibin = hProtonWeights->FindBin(pt);
    w = hProtonWeights->GetBinContent(ibin);*/
    w = fProtonWeights->Eval(pt);

    return w;
}
//____________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskPythiaCoalescence::GetDeuteronWeight(Double_t pt_prot, Double_t pt_neut)
{

    // Initialization
    Double_t w(1);

    // Double_t wp = hProtonWeights->GetBinContent(hProtonWeights->FindBin(pt_prot));
    // Double_t wn = hProtonWeights->GetBinContent(hProtonWeights->FindBin(pt_neut));

    Double_t wp = fProtonWeights->Eval(pt_prot);
    Double_t wn = fProtonWeights->Eval(pt_neut);

    Double_t spin_factor = 3.0 / 4.0;
    Double_t isospin_factor = 1.0 / 2.0;
    w = spin_factor * isospin_factor * wp * wn;

    return w;
}
//____________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskPythiaCoalescence::GetDeuteronWeight(Double_t pt_prot, Double_t pt_neut, Int_t prot_Reg, Int_t neut_Reg)
{

    // Initialization
    Double_t w(1), wp(1), wn(1);

    // Proton
    if (prot_Reg == 0)
        wp = hWeightToward->GetBinContent(hWeightToward->FindBin(pt_prot));
    if (prot_Reg == 1)
        wp = hWeightTransv->GetBinContent(hWeightTransv->FindBin(pt_prot));
    if (prot_Reg == 2)
        wp = hWeightAway->GetBinContent(hWeightAway->FindBin(pt_prot));

    // Neutron
    if (neut_Reg == 0)
        wn = hWeightToward->GetBinContent(hWeightToward->FindBin(pt_neut));
    if (neut_Reg == 1)
        wn = hWeightTransv->GetBinContent(hWeightTransv->FindBin(pt_neut));
    if (neut_Reg == 2)
        wn = hWeightAway->GetBinContent(hWeightAway->FindBin(pt_neut));

    // Statistical Factor
    Double_t spin_factor = 3.0 / 4.0;
    Double_t isospin_factor = 1.0 / 2.0;
    w = spin_factor * isospin_factor * wp * wn;

    return w;
}
//____________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskPythiaCoalescence::GetProtonWeight(Double_t pt, Int_t region)
{

    // Initialization
    Double_t w(1);

    // Proton
    if (region == 0)
        w = hWeightToward->GetBinContent(hWeightToward->FindBin(pt));
    if (region == 1)
        w = hWeightTransv->GetBinContent(hWeightTransv->FindBin(pt));
    if (region == 2)
        w = hWeightAway->GetBinContent(hWeightAway->FindBin(pt));

    return w;
}
//____________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskPythiaCoalescence::GetSpatialDistance(TLorentzVector p_proton, TLorentzVector p_neutron, TVector3 beta_vect)
{

    // Initialization
    Double_t deltaX(0);

    // Source Radius: m_{T} dependence
    Double_t r_source_proton = hSourceSize->GetBinContent(hSourceSize->FindBin(p_proton.Mt()));
    Double_t r_source_neutron = hSourceSize->GetBinContent(hSourceSize->FindBin(p_neutron.Mt()));

    // Proton Coordinates in the Lab
    Double_t rp = TMath::Abs(gRandom->Gaus(0.0, r_source_proton));
    Double_t theta_p = p_proton.Theta();
    Double_t phi_p = p_proton.Phi();
    Double_t xp = rp * TMath::Sin(theta_p) * TMath::Cos(phi_p);
    Double_t yp = rp * TMath::Sin(theta_p) * TMath::Sin(phi_p);
    Double_t zp = rp * TMath::Cos(theta_p);
    TLorentzVector r_proton_lab;
    r_proton_lab.SetXYZT(xp, yp, zp, 0);

    // Neutron Coordinates in the Lab
    Double_t rn = TMath::Abs(gRandom->Gaus(0.0, r_source_neutron));
    Double_t theta_n = p_neutron.Theta();
    Double_t phi_n = p_neutron.Phi();
    Double_t xn = rn * TMath::Sin(theta_n) * TMath::Cos(phi_n);
    Double_t yn = rn * TMath::Sin(theta_n) * TMath::Sin(phi_n);
    Double_t zn = rn * TMath::Cos(theta_n);
    TLorentzVector r_neutron_lab;
    r_neutron_lab.SetXYZT(xn, yn, zn, 0);

    // Fill Control Histograms
    hSourceRadius_Prot->Fill(rp);
    hSourceRadius_Neut->Fill(rn);

    // Distance in the Deuteron Frame
    TLorentzVector r_proton_prime = LorentzTransform(r_proton_lab, beta_vect);
    TLorentzVector r_neutron_prime = LorentzTransform(r_neutron_lab, beta_vect);
    deltaX = (r_proton_prime - r_neutron_prime).P();

    // Control Histograms
    Double_t distance_lab = (r_proton_lab - r_neutron_lab).P();
    Double_t distance_deut = (r_proton_prime - r_neutron_prime).P();
    hDistanceLab->Fill(distance_lab);
    hDistanceDeut->Fill(distance_deut);
    hDistanceDiff->Fill(distance_lab - distance_deut);

    return deltaX;
}
//____________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskPythiaCoalescence::DoCoalescence(Double_t deltaX, Double_t deltaP, Double_t sigma_p, const char *func)
{

    // Initialization
    Bool_t doCoalescence = (kFALSE);

    // Constants
    Double_t sigma_x = 3.2; //(fm) Parameter related to deuteron radius (arxiv.org/abs/1812.05175)
    Double_t max = fDeuteronWF->GetMaximum();

    // Wave Function in Momentum Space
    Double_t Psi_momentum = TMath::Exp(-(deltaP * deltaP) / (sigma_p * sigma_p));

    // Wave Function in Coordinate Space
    Double_t Psi_space(0);
    if (strcmp(func, "Gaus") == 0)
        Psi_space = TMath::Exp(-(deltaX * deltaX) / (sigma_x * sigma_x));
    if (strcmp(func, "Arg") == 0)
        Psi_space = (1.0 / max) * fDeuteronWF->Eval(deltaX);

    // Wigner Density
    Double_t wigner = Psi_space * Psi_momentum;

    // Coalescence Probability
    Double_t rndm = gRandom->Uniform(0.0, 1.0);
    if (rndm < wigner)
        doCoalescence = kTRUE;

    return doCoalescence;
}
//____________________________________________________________________________________________________________________________________________________
Int_t AliAnalysisTaskPythiaCoalescence::GetLeadingParticle()
{

    // Initialization
    Int_t ID_leading_particle(0);
    Double_t pt_max(0);

    // Loop over Generated Particles
    for (Int_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++)
    {

        // Get Particle
        AliMCParticle *particle = (AliMCParticle *)fMCEvent->GetTrack(i);
        if (!particle)
            continue;
        if (particle->Charge() == 0)
            continue;
        if (IsInjectedParticle(particle))
            continue;
        if (TMath::Abs(particle->Eta()) > 0.8)
            continue;
        if (particle->MCStatusCode() <= 0)
            continue;
        if (particle->Pt() < 0.15)
            continue;

        // Primary Particle
        Double_t x = particle->Xv();
        Double_t y = particle->Yv();
        Double_t z = particle->Zv();
        Double_t r = TMath::Sqrt(x * x + y * y + z * z);

        Bool_t isFSParticle = (kFALSE);
        if (TMath::Abs(particle->PdgCode()) == 211)
            isFSParticle = kTRUE;
        if (TMath::Abs(particle->PdgCode()) == 321)
            isFSParticle = kTRUE;
        if (TMath::Abs(particle->PdgCode()) == 2212)
            isFSParticle = kTRUE;
        if (TMath::Abs(particle->PdgCode()) == 11)
            isFSParticle = kTRUE;
        if (!isFSParticle)
            continue;

        if (particle->Pt() > pt_max)
        {

            pt_max = particle->Pt();
            ID_leading_particle = i;
        }
    }

    return ID_leading_particle;
}
//____________________________________________________________________________________________________________________________________________________
Int_t AliAnalysisTaskPythiaCoalescence::GetTransverseMultiplicity(Int_t leading_particle_ID)
{

    // Initialization
    Int_t mult_Transverse(0);

    // Leading Particle
    AliMCParticle *leading_particle = (AliMCParticle *)fMCEvent->GetTrack(leading_particle_ID);

    // Loop over Generated Particles
    for (Int_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++)
    {

        // Get Particle
        AliMCParticle *particle = (AliMCParticle *)fMCEvent->GetTrack(i);
        if (!particle)
            continue;
        if (particle->Charge() == 0)
            continue;
        if (IsInjectedParticle(particle))
            continue;
        if (TMath::Abs(particle->Eta()) > 0.8)
            continue;
        if (particle->MCStatusCode() <= 0)
            continue;
        if (particle->Pt() < 0.15)
            continue;

        // Primary Particle
        Double_t x = particle->Xv();
        Double_t y = particle->Yv();
        Double_t z = particle->Zv();
        Double_t r = TMath::Sqrt(x * x + y * y + z * z);
        // if (r>1e-10) continue;

        Bool_t isFSParticle = (kFALSE);
        if (TMath::Abs(particle->PdgCode()) == 211)
            isFSParticle = kTRUE;
        if (TMath::Abs(particle->PdgCode()) == 321)
            isFSParticle = kTRUE;
        if (TMath::Abs(particle->PdgCode()) == 2212)
            isFSParticle = kTRUE;
        if (TMath::Abs(particle->PdgCode()) == 11)
            isFSParticle = kTRUE;
        if (!isFSParticle)
            continue;

        Double_t phi_particle = TVector2::Phi_0_2pi(particle->Phi());
        Double_t phi_leading = TVector2::Phi_0_2pi(leading_particle->Phi());

        if (IsParticleInTransverseRegion(phi_particle, phi_leading))
            mult_Transverse++;
    }

    return mult_Transverse;
}
//____________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskPythiaCoalescence::IsParticleInTowardRegion(Double_t phi, Double_t phi_leading)
{

    // Initialization
    Bool_t isInTowardRegion = kFALSE;

    // DeltaPhi
    Double_t delta_phi = (180.0 / TMath::Pi()) * TVector2::Phi_0_2pi(phi - phi_leading);

    if (delta_phi >= 0.0 && delta_phi < 60.0)
        isInTowardRegion = kTRUE;
    if (delta_phi >= 300.0 && delta_phi <= 360.0)
        isInTowardRegion = kTRUE;

    return isInTowardRegion;
}
//____________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskPythiaCoalescence::IsParticleInTransverseRegion(Double_t phi, Double_t phi_leading)
{

    // Initialization
    Bool_t isInTransverseRegion = kFALSE;

    // DeltaPhi
    Double_t delta_phi = (180.0 / TMath::Pi()) * TVector2::Phi_0_2pi(phi - phi_leading);

    if (delta_phi >= 60.0 && delta_phi < 120.0)
        isInTransverseRegion = kTRUE;
    if (delta_phi >= 240.0 && delta_phi < 300.0)
        isInTransverseRegion = kTRUE;

    return isInTransverseRegion;
}
//____________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskPythiaCoalescence::IsParticleInAwayRegion(Double_t phi, Double_t phi_leading)
{

    // Initialization
    Bool_t isInAwayRegion = kFALSE;

    // DeltaPhi
    Double_t delta_phi = (180.0 / TMath::Pi()) * TVector2::Phi_0_2pi(phi - phi_leading);
    if (delta_phi >= 120.0 && delta_phi < 240.0)
        isInAwayRegion = kTRUE;

    return isInAwayRegion;
}
//____________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskPythiaCoalescence::IsInjectedParticle(AliMCParticle *particle)
{

    // Initialization
    Bool_t isInjected = kFALSE;

    if (TMath::Abs(particle->PdgCode()) == 1000010020)
        isInjected = kTRUE; // Deuteron
    if (TMath::Abs(particle->PdgCode()) == 1000020030)
        isInjected = kTRUE; // Helium-3
    if (TMath::Abs(particle->PdgCode()) == 1000010030)
        isInjected = kTRUE; // Triton
    if (TMath::Abs(particle->PdgCode()) == 1000020040)
        isInjected = kTRUE; // Helium-4
    if (TMath::Abs(particle->PdgCode()) == 1010010030)
        isInjected = kTRUE; // Hypertriton

    return isInjected;
}
//____________________________________________________________________________________________________________________________________________________
TLorentzVector AliAnalysisTaskPythiaCoalescence::LorentzTransform(TLorentzVector R, TVector3 beta_vect)
{

    // Inizialization
    TLorentzVector R_prime(0.0, 0.0, 0.0, 0.0);

    // Beta Components
    Double_t Bx = beta_vect.X();
    Double_t By = beta_vect.Y();
    Double_t Bz = beta_vect.Z();

    // Beta & Gamma
    Double_t beta = TMath::Sqrt(Bx * Bx + By * By + Bz * Bz);
    if (beta >= 1.0)
    {
        return R_prime;
    }
    Double_t gamma = 1.0 / TMath::Sqrt(1.0 - (beta * beta));

    // Coordinates in the Lab System
    Double_t t = R.T();
    Double_t x = R.X();
    Double_t y = R.Y();
    Double_t z = R.Z();

    // Coordinates in the Deuteron Frame
    Double_t t_prime = gamma * t - gamma * Bx * x - gamma * By * y - gamma * Bz * z;
    Double_t x_prime = -gamma * Bx * t + (1.0 + (gamma - 1.0) * Bx * Bx / (beta * beta)) * x + (gamma - 1.0) * (Bx * By / (beta * beta)) * y + (gamma - 1.0) * (Bx * Bz / (beta * beta)) * z;
    Double_t y_prime = -gamma * By * t + (gamma - 1.0) * (Bx * By / (beta * beta)) * x + (1.0 + (gamma - 1.0) * By * By / (beta * beta)) * y + (gamma - 1.0) * (By * Bz / (beta * beta)) * z;
    Double_t z_prime = -gamma * Bz * t + (gamma - 1.0) * (Bx * Bz / (beta * beta)) * x + (gamma - 1.0) * (By * Bz / (beta * beta)) * y + (1.0 + (gamma - 1.0) * Bz * Bz / (beta * beta)) * z;

    // Set Coordinates
    R_prime.SetXYZT(x_prime, y_prime, z_prime, t_prime);

    return R_prime;
}
//____________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskPythiaCoalescence::Terminate(Option_t *)
{

    fOutputList = dynamic_cast<TList *>(GetOutputData(1));
    if (!fOutputList)
        return;
}
//____________________________________________________________________________________________________________________________________________________
