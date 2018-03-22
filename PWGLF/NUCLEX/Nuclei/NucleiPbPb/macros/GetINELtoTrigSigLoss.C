#if !defined(__CINT__) || defined(__MAKECINT__) || defined(__CLING__)
#include <map>
#include <string>
#include <sstream>
#include <vector>

#include <TStyle.h>
#include <TFile.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TF1.h>

#include "MyUtils.h"
#endif

using std::map;
using std::string;
using std::vector;

typedef enum { kpp13_LHC17_Pythia8_EvtCut, kpp13_LHC17_Pythia8, kNDataSet } EDataSet_Type;

const char *base_dir = "/Users/ycorrales/Alice/PWGLF_SPECTRA/INELtokINT7_SigLoss";

const char *prtName_ltx[] = { "#pi^{+} + #pi^{-}", "K^{+} + K^{-}", "p + #bar{p}", "#phi",  "K^{*}",
                              "K^{0}_{s}",         "#Lambda",       "#Xi",         "#Omega" };
const vector<string> v_prtName_ltx(prtName_ltx, prtName_ltx + 9);

const char *prtName_str[] = { "PiPos", "PiNeg", "KaPos", "KaNeg",  "PrPos", "PrNeg",
                              "Phi",   "Ks",    "K0",    "Lambda", "Xi",    "Om" };
const unsigned int kNspc = 12;
const vector<string> v_prtName_str(prtName_str, prtName_str + kNspc);
//
const char *buf = "_EvCut";
//
// const string rat_lab[3] = {"", "(K^{+} + K^{-}) / (#pi^{+} + #pi^{-})", "(p + #bar{p}) / (#pi^{+} + #pi^{-})"};

//
// Main Function
//________________________________________________________________________________
void GetINELtoTrigSigLoss()
{
  gStyle->SetOptStat(0);
  // gStyle->SetOptTitle(0);

  TH1F *frame = new TH1F("frame", "frame", 10000, 0., 20);
  frame->GetYaxis()->SetRangeUser(0.98, 1.1);
  frame->GetYaxis()->SetTitle("Signal Loss");
  frame->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  frame->GetYaxis()->SetTitleOffset(1.5);

  vector<EDataSet_Type> v_dataset{ kpp13_LHC17_Pythia8_EvtCut, kpp13_LHC17_Pythia8 };

  string flName = base_dir + string("/AnTask_output/AR_SigLoss_NEL-kINT7_LHC17d20a2-Pythia8.root");
  TFile *flin = TFile::Open(flName.data());
  for (EDataSet_Type &dtSet : v_dataset) {

    // output file
    string ofl_name("SignalLoss");
    string lst_name("cListSpectraINEL");
    if (dtSet == kpp13_LHC17_Pythia8_EvtCut) {
      ofl_name += string("_EvCut");
      lst_name += string("_EvCut");
    }
    ofl_name += string(".root");

    TFile *ofl = TFile::Open(ofl_name.data(), "RECREATE");

    TList *lst = (TList *)flin->Get(lst_name.data());
    if (!lst)
      gMyUtils->Error("__FUNCTION__", "Wrong list name");

    frame->SetName(lst_name.c_str());
    frame->SetTitle(lst_name.c_str());

    TH2F *fHistRecEvtSelected = (TH2F *)lst->FindObject("fHistRecEvtSelected");
    const int lastBinSelected = fHistRecEvtSelected->GetXaxis()->GetNbins();
    TH1F *fHistAccepted = (TH1F *)fHistRecEvtSelected->ProjectionX("fHistAcceptedEvents", lastBinSelected,
                                                                   lastBinSelected); // Accepted Events
    TH2F *fHistGenEvtSelected = (TH2F *)lst->FindObject("fHistGenEvtSelected");
    TH1F *fHistTrueINELgtZero = (TH1F *)fHistGenEvtSelected->ProjectionX("fHistTrueINELgtZero", 2, 2); // TrueINEL>0
                                                                                                       // gVtxZ<10cm
    TH1F *fHistEvtLoss = (TH1F *)fHistAccepted->Clone("fHistEvtLoss");
    fHistEvtLoss->Divide(fHistAccepted, fHistTrueINELgtZero, 1, 1);
    fHistEvtLoss->SetTitle("fHistEvtLoss;Cent[Percentile];#epsilon_{c}[Counts]");
    fHistEvtLoss->GetXaxis()->CenterTitle();
    fHistEvtLoss->GetYaxis()->CenterTitle();

    ofl->cd();
    fHistAccepted->Write();
    fHistTrueINELgtZero->Write();
    fHistEvtLoss->Write();

    vector<TH2D *> fHistSigLoss(kNspc, nullptr);
    vector<TH2D *> fHistGenINEL(kNspc, nullptr);
    vector<TH2D *> fHistGenAccepted(kNspc, nullptr);
    vector<TH2D *> fHistGenZvtxTrue(kNspc, nullptr);

    //        0-1%, 1-5%, 5-10%, 10-20%,         20-30%, 30-40%, 40-60%,         60-80%,         80-100%
    //-.5-0%, 0-1%, 1-5%, 5-10%, 10-15%, 15-20%, 20-30%, 30-40%, 40-50%, 50-60%, 60-70%, 70-80%, 80-90%, 90-100%
    const unsigned int nCentBins = 10;
    const unsigned int lowCentBin[nCentBins] = { 2, 3, 4, 5, 7, 8, 9, 10, 12, 2 };
    const unsigned int upCentBin[nCentBins] = { 2, 3, 4, 6, 7, 8, 9, 11, 14, 14 };

    bool showOne = true;
    for (unsigned int ispc(0); ispc < kNspc; ++ispc) {
      TH3F *fHistCentPtAccepted = (TH3F *)lst->FindObject(Form("fHistPtAftAllEventCuts%s", v_prtName_str[ispc].data()));
      TH3F *fHistCentPtTrueINELgtZ =
        (TH3F *)lst->FindObject(Form("fHistPtGoodTrueINELgtZ%s", v_prtName_str[ispc].data()));

      std::ostringstream fHistName;
      TList *lst_part = new TList();
      lst_part->SetName(v_prtName_str[ispc].data());
      lst_part->SetOwner(true);
      for (unsigned int i_cent = 0; i_cent < nCentBins; ++i_cent) {
        float lowEdge = fHistCentPtAccepted->GetXaxis()->GetBinLowEdge(lowCentBin[i_cent]);
        float upEdge = fHistCentPtAccepted->GetXaxis()->GetBinUpEdge(upCentBin[i_cent]);
        std::ostringstream os_suffix;
        os_suffix << "_" << lowEdge << "-" << upEdge;

        fHistName.str(string());
        fHistName << "GenSpectraAccepted"
                  << "_" << v_prtName_str[ispc] << os_suffix.str();
        TH1D *fHistPtAccepted = fHistCentPtAccepted->ProjectionZ(fHistName.str().data(), lowCentBin[i_cent],
                                                                 upCentBin[i_cent], lastBinSelected, lastBinSelected);
        fHistPtAccepted->SetTitle(fHistName.str().data());

        fHistName.str(string());
        fHistName << "GenSpectraTrueINELgtZ"
                  << "_" << v_prtName_str[ispc] << os_suffix.str();

        TH1D *fHistPtTrueINELgtZ =
          fHistCentPtTrueINELgtZ->ProjectionZ(fHistName.str().data(), lowCentBin[i_cent], upCentBin[i_cent], 2, 2);
        fHistPtTrueINELgtZ->SetTitle(fHistName.str().data());

        fHistPtAccepted->Sumw2();
        fHistPtTrueINELgtZ->Sumw2();

        //        double evtLoss = fHistEvtLoss->GetBinContent(i_cent);
        //        if (showOne) {
        //          Printf("Evt Loss for cent bin %d and part %d, %f", i_cent, ispc, evtLoss);
        //        }
        fHistName.str(string());
        fHistName << "sgnLoss"
                  << "_" << v_prtName_str[ispc] << os_suffix.str();
        TH1D *fHistSigLoss = (TH1D *)fHistPtTrueINELgtZ->Clone(fHistName.str().data());
        fHistSigLoss->SetTitle(Form("%s;#it{p}_{T};SigLoss * #epsilon_{c}[a.u]", fHistName.str().data()));
        fHistSigLoss->Divide(fHistPtTrueINELgtZ, fHistPtAccepted, 1, 1, "B");
        fHistSigLoss->SetLineColor(ispc + 1);

        lst_part->Add(fHistSigLoss);
      }
      showOne = false;
      ofl->cd();
      lst_part->Write("", TObject::kSingleKey);
    }
  }
  flin->Close();
  // delete frame;
}
