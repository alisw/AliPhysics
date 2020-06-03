#if !defined(__CINT__) || defined(__CLING__)

#include <TROOT.h>
#include <Riostream.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>

#include "AliOADBContainer.h"

#endif

enum
{
    kLHC18q,
    kLHC18r
};

int ComputeTPCCalibrations(int period = kLHC18q, TString infilename = "$HOME/Downloads/AnalysisResults.root", TString dirname = "PWGHF_D2H_QnVectorTenderPhiDistr")
{
    if(period != kLHC18q && period != kLHC18r)
    {
        std::cerr << "ERROR: period not implemented! Exit" << std::endl;
        return -1;
    }

    gStyle->SetPadRightMargin(0.035);
    gStyle->SetPadLeftMargin(0.16);
    gStyle->SetPadBottomMargin(0.1);
    gStyle->SetPadTopMargin(0.035);
    gStyle->SetTitleSize(0.045, "xy");
    gStyle->SetLabelSize(0.04, "xy");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetHistLineWidth(2);
    gStyle->SetOptStat(0);
    gROOT->SetBatch(kTRUE);

    int runs[150];
    int nRuns = 0;
    TString periodSuffix = "";
    if (period == kLHC18q)
    {
        nRuns = 126;
        std::array<int, 126> runlist = {296623, 296622, 296621, 296619, 296618, 296616, 296615, 296594, 296553, 296552, 296551,
                                        296550, 296549, 296548, 296547, 296516, 296512, 296511, 296510, 296509, 296472, 296433, 296424, 296423, 296420,
                                        296419, 296415, 296414, 296383, 296381, 296380, 296379, 296378, 296377, 296376, 296375, 296312, 296309, 296304,
                                        296303, 296280, 296279, 296273, 296270, 296269, 296247, 296246, 296244, 296243, 296242, 296241, 296240, 296198,
                                        296197, 296196, 296195, 296194, 296192, 296191, 296143, 296142, 296135, 296134, 296133, 296132, 296123, 296074,
                                        296066, 296065, 296063, 296062, 296060, 296016, 295942, 295941, 295937, 295936, 295913, 295910, 295909, 295861,
                                        295860, 295859, 295856, 295855, 295854, 295853, 295831, 295829, 295826, 295825, 295822, 295819, 295818, 295816,
                                        295791, 295788, 295786, 295763, 295762, 295759, 295758, 295755, 295754, 295725, 295723, 295721, 295719, 295718,
                                        295717, 295714, 295712, 295676, 295675, 295673, 295668, 295667, 295666, 295615, 295612, 295611, 295610, 295589,
                                        295588, 295586, 295585};
        std::copy(runlist.begin(), runlist.end(), runs);
        periodSuffix = "18q";
    }
    else if (period == kLHC18r)
    {
        nRuns = 97;
        std::array<int, 97> runlist = {297595, 297590, 297588, 297558, 297544, 297542, 297541, 297540, 297537, 297512, 297483,
                                       297481, 297479, 297452, 297451, 297450, 297446, 297442, 297441, 297415, 297414, 297413, 297406, 297405, 297380,
                                       297379, 297372, 297367, 297366, 297363, 297336, 297335, 297333, 297332, 297317, 297311, 297310, 297278, 297222,
                                       297221, 297218, 297196, 297195, 297193, 297133, 297132, 297129, 297128, 297124, 297123, 297119, 297118, 297117,
                                       297085, 297035, 297031, 296966, 296941, 296938, 296935, 296934, 296932, 296931, 296930, 296903, 296900, 296899,
                                       296894, 296852, 296851, 296850, 296848, 296839, 296838, 296836, 296835, 296799, 296794, 296793, 296790, 296787,
                                       296786, 296785, 296784, 296781, 296752, 296694, 296693, 296691, 296690, 296749, 296750, 296849, 296890, 297029,
                                       297194, 297219};
        std::copy(runlist.begin(), runlist.end(), runs);
        periodSuffix = "18r";
    }

    TFile *infile = TFile::Open(infilename.Data());
    if (!infile || !infile->IsOpen())
        return -1;

    AliOADBContainer *OADBcontPosEta[9];
    AliOADBContainer *OADBcontNegEta[9];
    for (int iCent = 0; iCent < 9; iCent++)
    {
        OADBcontPosEta[iCent] = new AliOADBContainer(Form("fphidistr_poseta_%d_%d", iCent * 10, (iCent + 1) * 10));
        OADBcontNegEta[iCent] = new AliOADBContainer(Form("fphidistr_negeta_%d_%d", iCent * 10, (iCent + 1) * 10));
    }

    TCanvas *cExample = new TCanvas("cExample", "", 800, 800);
    for (int iRun = 0; iRun < nRuns; iRun++)
    {
        TH2F *hPhiDistrPosEtaVsCentr = (TH2F *)infile->Get(Form("%s/fPhiVsCentrTPCTPCPosEta_%d", dirname.Data(), runs[iRun]));
        TH2F *hPhiDistrNegEtaVsCentr = (TH2F *)infile->Get(Form("%s/fPhiVsCentrTPCTPCNegEta_%d", dirname.Data(), runs[iRun]));
        if (!hPhiDistrPosEtaVsCentr || !hPhiDistrNegEtaVsCentr)
            return runs[iRun];
        hPhiDistrPosEtaVsCentr->SetDirectory(0);
        hPhiDistrNegEtaVsCentr->SetDirectory(0);
        for (int iCent = 0; iCent < 9; iCent++)
        {
            TH1D *hPhiTPCPosEta = static_cast<TH1D *>(hPhiDistrPosEtaVsCentr->ProjectionY("fPhiTPCPosEta", iCent + 1, iCent + 1));
            TH1D *hPhiTPCNegEta = static_cast<TH1D *>(hPhiDistrNegEtaVsCentr->ProjectionY("fPhiTPCNegEta", iCent + 1, iCent + 1));
            hPhiTPCPosEta->Scale(hPhiTPCPosEta->GetNbinsX() / hPhiTPCPosEta->Integral());
            hPhiTPCNegEta->Scale(hPhiTPCNegEta->GetNbinsX() / hPhiTPCNegEta->Integral());
            TLegend *leg = new TLegend(0.4, 0.3, 0.8, 0.6);
            leg->SetTextSize(0.045);
            leg->AddEntry("", Form("Run %d, %02d-%02d%%", runs[iRun], iCent*10, (iCent+1)*10), "");
            leg->AddEntry(hPhiTPCPosEta, "0 < #eta < 0.8", "l");
            leg->AddEntry(hPhiTPCNegEta, "-0.8 < #eta < 0", "l");
            hPhiTPCPosEta->SetLineColor(kRed + 1);
            hPhiTPCNegEta->SetLineColor(kBlue + 1);
            hPhiTPCPosEta->SetLineWidth(2);
            hPhiTPCNegEta->SetLineWidth(2);
            hPhiTPCPosEta->SetDirectory(0);
            hPhiTPCNegEta->SetDirectory(0);
            cExample->cd()->DrawFrame(0, 0., 2 * TMath::Pi(), 1.2, ";#varphi;N(#varphi) / N_{tot}");
            hPhiTPCPosEta->DrawCopy("histsame");
            hPhiTPCNegEta->DrawCopy("histsame");
            leg->Draw();
            if(iCent == 0 && iRun == 0)
                cExample->SaveAs(Form("PhiDistr_%s.pdf[", periodSuffix.Data()));
            cExample->SaveAs(Form("PhiDistr_%s.pdf", periodSuffix.Data()));
            if(iCent == 8 && iRun == nRuns-1)
                cExample->SaveAs(Form("PhiDistr_%s.pdf]", periodSuffix.Data()));

            OADBcontPosEta[iCent]->AppendObject(hPhiTPCPosEta->Clone(), runs[iRun], runs[iRun]);
            OADBcontNegEta[iCent]->AppendObject(hPhiTPCNegEta->Clone(), runs[iRun], runs[iRun]);
        }
    }

    TFile outfile(Form("calibTPCRun2%sVtx10", periodSuffix.Data()), "recreate");
    for (int iCent = 0; iCent < 9; iCent++)
    {
        OADBcontPosEta[iCent]->Write();
        OADBcontNegEta[iCent]->Write();
    }
    outfile.Close();

    return 0;
}
