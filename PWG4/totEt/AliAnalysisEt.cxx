#include "AliAnalysisEt.h"
#include "TMath.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include "AliAnalysisEtCuts.h"
ClassImp(AliAnalysisEt);


AliAnalysisEt::AliAnalysisEt() :
        fHistogramNameSuffix("")
        ,fPdgDB(0)
        ,fSumEt(0)
        ,fSumEtAcc(0)
        ,fTotEt(0)
        ,fTotEtAcc(0)
        ,fTotNeutralEt(0)
        ,fTotNeutralEtAcc(0)
        ,fTotChargedEt(0)
        ,fTotChargedEtAcc(0)
        ,fMultiplicity(0)
        ,fChargedMultiplicity(0)
        ,fNeutralMultiplicity(0)
        ,fEtaCut(EtCommonCuts::kEtaCut)
        ,fEtaCutAcc(0)
        ,fPhiCutAccMin(0)
        ,fPhiCutAccMax(360.)
        ,fDetectorRadius(460.)
        ,fVertexXCut(0)
        ,fVertexYCut(0)
        ,fVertexZCut(0)
        ,fIPxyCut(0)
        ,fIPzCut(0)
        ,fClusterEnergyCut(EtCommonCuts::kClusterEnergyCut)
        ,fTrackPtCut(EtCommonCuts::kTrackPtCut)
        ,fSingleCellEnergyCut(0)
        ,fHistEt(0)
        ,fHistChargedEt(0)
        ,fHistNeutralEt(0)
        ,fHistEtAcc(0)
        ,fHistChargedEtAcc(0)
        ,fHistNeutralEtAcc(0)
        ,fHistMult(0)
        ,fHistChargedMult(0)
        ,fHistNeutralMult(0)
        ,fHistPhivsPtPos(0)
        ,fHistPhivsPtNeg(0)
        ,fHistBaryonEt(0)
        ,fHistAntiBaryonEt(0)
        ,fHistMesonEt(0)
        ,fHistBaryonEtAcc(0)
        ,fHistAntiBaryonEtAcc(0)
        ,fHistMesonEtAcc(0)
        ,fHistEtRecvsEtMC(0)
        ,fHistTMDeltaR(0)
{

}

AliAnalysisEt::~AliAnalysisEt()
{

}

void AliAnalysisEt::FillOutputList(TList *list)
{
    list->Add(fHistEt);
    list->Add(fHistChargedEt);
    list->Add(fHistNeutralEt);

    list->Add(fHistEtAcc);
    list->Add(fHistChargedEtAcc);
    list->Add(fHistNeutralEtAcc);

    list->Add(fHistMult);
    list->Add(fHistChargedMult);
    list->Add(fHistNeutralMult);

    list->Add(fHistPhivsPtPos);
    list->Add(fHistPhivsPtNeg);

    list->Add(fHistBaryonEt);
    list->Add(fHistAntiBaryonEt);
    list->Add(fHistMesonEt);

    list->Add(fHistBaryonEtAcc);
    list->Add(fHistAntiBaryonEtAcc);
    list->Add(fHistMesonEtAcc);

    list->Add(fHistEtRecvsEtMC);

    list->Add(fHistTMDeltaR);
}

void AliAnalysisEt::Init()
{
    fPdgDB = new TDatabasePDG();
}

void AliAnalysisEt::CreateHistograms()
{

    TString histname = "fHistEt" + fHistogramNameSuffix;

    fHistEt = new TH1F(histname.Data(), "Total E_{T} Distribution", 1000, 0.00, 99);
    fHistEt->GetXaxis()->SetTitle("E_{T} (GeV/c^{2})");
    fHistEt->GetYaxis()->SetTitle("dN/dE_{T} (c^{2}/GeV)");

    histname = "fHistChargedEt" + fHistogramNameSuffix;
    fHistChargedEt = new TH1F(histname.Data(), "Total Charged E_{T} Distribution", 1000, 0.00, 99);
    fHistChargedEt->GetXaxis()->SetTitle("E_{T} (GeV/c^{2})");
    fHistChargedEt->GetYaxis()->SetTitle("dN/dE_{T} (c^{2}/GeV)");

    histname = "fHistNeutralEt" + fHistogramNameSuffix;
    fHistNeutralEt = new TH1F(histname.Data(), "Total Neutral E_{T} Distribution", 1000, 0.00, 99);
    fHistNeutralEt->GetXaxis()->SetTitle("E_{T} (GeV/c^{2})");
    fHistNeutralEt->GetYaxis()->SetTitle("dN/dE_{T} (c^{2}/GeV)");

    histname = "fHistEtAcc" + fHistogramNameSuffix;
    fHistEtAcc = new TH1F(histname.Data(), "Total E_{T} Distribution in Acceptance", 1000, 0.00, 99);
    fHistEtAcc->GetXaxis()->SetTitle("E_{T} (GeV/c^{2})");
    fHistEtAcc->GetYaxis()->SetTitle("dN/dE_{T} (c^{2}/GeV)");

    histname = "fHistChargedEtAcc" + fHistogramNameSuffix;
    fHistChargedEtAcc = new TH1F(histname.Data(), "Total Charged E_{T} Distribution in Acceptance", 1000, 0.00, 99);
    fHistChargedEtAcc->GetXaxis()->SetTitle("E_{T} (GeV/c^{2})");
    fHistChargedEtAcc->GetYaxis()->SetTitle("dN/dE_{T} (c^{2}/GeV)");

    histname = "fHistNeutralEtAcc" + fHistogramNameSuffix;
    fHistNeutralEtAcc = new TH1F(histname.Data(), "Total Neutral E_{T} Distribution in Acceptance", 1000, 0.00, 99);
    fHistNeutralEtAcc->GetXaxis()->SetTitle("E_{T} (GeV/c^{2})");
    fHistNeutralEtAcc->GetYaxis()->SetTitle("dN/dE_{T} (c^{2}/GeV)");
    std::cout << histname << std::endl;
    histname = "fHistMult" + fHistogramNameSuffix;
    fHistMult = new TH1F(histname.Data(), "Total Multiplicity", 200, 0, 199);
    fHistMult->GetXaxis()->SetTitle("N");
    fHistMult->GetYaxis()->SetTitle("Multiplicity");

    histname = "fHistChargedMult" + fHistogramNameSuffix;
    fHistChargedMult = new TH1F(histname.Data(), "Charged Multiplicity", 200, 0, 199);
    fHistChargedMult->GetXaxis()->SetTitle("N");
    fHistChargedMult->GetYaxis()->SetTitle("Multiplicity");

    histname = "fHistNeutralMult" + fHistogramNameSuffix;
    fHistNeutralMult = new TH1F(histname.Data(), "Charged Multiplicity", 200, 0, 199);
    fHistNeutralMult->GetXaxis()->SetTitle("N");
    fHistNeutralMult->GetYaxis()->SetTitle("Multiplicity");

    histname = "fHistPhivsPtPos" + fHistogramNameSuffix;
    fHistPhivsPtPos = new TH2F(histname.Data(), "Phi vs pT of positively charged tracks hitting the calorimeter", 	200, 0, 2*TMath::Pi(), 2000, 0, 100);

    histname = "fHistPhivsPtNeg" + fHistogramNameSuffix;
    fHistPhivsPtNeg = new TH2F(histname.Data(), "Phi vs pT of negatively charged tracks hitting the calorimeter", 	200, 0, 2*TMath::Pi(), 2000, 0, 100);

    histname = "fHistBaryonEt" + fHistogramNameSuffix;
    fHistBaryonEt = new TH1F(histname.Data(), "E_{T} for baryons",  1000, 0.0001, 100);

    histname = "fHistAntiBaryonEt" + fHistogramNameSuffix;
    fHistAntiBaryonEt = new TH1F(histname.Data(), "E_{T} for anti baryons",  1000, 0.0001, 100);

    histname = "fHistMesonEt" + fHistogramNameSuffix;
    fHistMesonEt = new TH1F(histname.Data(), "E_{T} for mesons",  1000, 0.0001, 100);

    histname = "fHistBaryonEtAcc" + fHistogramNameSuffix;
    fHistBaryonEtAcc = new TH1F(histname.Data(), "E_{T} for baryons in calorimeter acceptance",  1000, 0.0001, 100);

    histname = "fHistAntiBaryonEtAcc" + fHistogramNameSuffix;
    fHistAntiBaryonEtAcc = new TH1F(histname.Data(), "E_{T} for anti baryons in calorimeter acceptance",  1000, 0.0001, 100);

    histname = "fHistMesonEtAcc" + fHistogramNameSuffix;
    fHistMesonEtAcc = new TH1F(histname.Data(), "E_{T} for mesons in calorimeter acceptance",  1000, 0.0001, 100);

    histname = "fHistEtRecvsEtMC" + fHistogramNameSuffix;
    fHistEtRecvsEtMC = new TH2F(histname.Data(), "Reconstructed E_{t} vs MC E_{t}", 1000, 0.000, 100, 1000, 0.0001, 100);

    histname = "fHistTMDeltaR" + fHistogramNameSuffix;
    fHistTMDeltaR = new TH1F(histname.Data(), "#Delta R for calorimeter clusters", 200, 0, 50);

}

void AliAnalysisEt::FillHistograms()
{
    fHistEt->Fill(fTotEt);
    fHistChargedEt->Fill(fTotChargedEt);
    fHistNeutralEt->Fill(fTotNeutralEt);

    fHistEtAcc->Fill(fTotEtAcc);
    fHistChargedEtAcc->Fill(fTotChargedEtAcc);
    fHistNeutralEtAcc->Fill(fTotNeutralEtAcc);

    fHistMult->Fill(fMultiplicity);
    fHistChargedMult->Fill(fChargedMultiplicity);
    fHistNeutralMult->Fill(fNeutralMultiplicity);

    /* // DS commented out non-fills to prevent compilation warnings
    fHistPhivsPtPos;
    fHistPhivsPtNeg;

    fHistBaryonEt;
    fHistAntiBaryonEt;
    fHistMesonEt;

    fHistBaryonEtAcc;
    fHistAntiBaryonEtAcc;
    fHistMesonEtAcc;

    fHistTMDeltaR;
    */
}

void AliAnalysisEt::ResetEventValues()
{
    fTotEt = 0;
    fTotEtAcc = 0;
    fTotNeutralEt = 0;
    fTotNeutralEtAcc = 0;
    fTotChargedEt  = 0;
    fTotChargedEtAcc = 0;
    fMultiplicity = 0;
    fChargedMultiplicity = 0;
    fNeutralMultiplicity = 0;
}
