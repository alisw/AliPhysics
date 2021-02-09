#include "src/Common.h"
#include "src/Utils.h"
#include "src/Plotting.h"
using namespace utils;
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TH1F.h>
#include <TDirectory.h>
#include <TLegend.h>

#include <string>
using std::string;

string output_name = "/Users/lbariogl/cernbox/Deuterons13TeV/macros/results/sec_ratio.root";

void Divider(bool short_mode = true){
  TFile file_TFF_TOF(kSecondariesOutput.data());
  TFile file_TFF_TPC(kSecondariesTPCoutput.data());
  TFile file_RooFit_TOF(kSecondariesOutputRooFit.data());
  TFile file_RooFit_TPC(kSecondariesTPCoutputRooFit.data());
  TFile file_out(output_name.data(),"RECREATE");

  const int kNCustomPtBins = 8;
  const double kCustomPtBins[9] = {0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6};

  int counter = 0;

  for (auto list_key : *file_TFF_TOF.GetListOfKeys()) {
    /// Preliminary operation to read the list and create an output dir
    if(short_mode)
      if(counter>0) return;
    counter++;
    if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos) continue;

    TDirectory* base_dir = file_out.mkdir(list_key->GetName());
    file_out.cd(list_key->GetName());

    TH1F* TofFraction[kCentLength];
    TH1F* TpcFraction[kCentLength];
    TH1F* TofTFFCorrection[kCentLength];
    TH1F* TpcTFFCorrection[kCentLength];


    for (int iC = 0; iC < kCentLength; ++iC) {
      TH1* hResTFFtof = (TH1*)file_TFF_TOF.Get(Form("%s/Results/hResTFF_%i",list_key->GetName(),iC));
      Requires(hResTFFtof,"hResTFFtof");
      TF1* prim_frac_TFF_TOF = hResTFFtof->GetFunction("fitFrac");
      Requires(prim_frac_TFF_TOF,"prim_frac_TFF_TOF");
      TofTFFCorrection[iC] = (TH1F*) hResTFFtof->Rebin(kNCustomPtBins,Form("TofTFFCorrection_%d",iC),kCustomPtBins);
      TofTFFCorrection[iC]->GetListOfFunctions()->Clear();


      TH1* hResTFFtpc = (TH1*)file_TFF_TPC.Get(Form("%s/Results/hResTFF_%i",list_key->GetName(),iC));
      Requires(hResTFFtpc,"hResTFFtpc");
      TF1* prim_frac_TFF_TPC = hResTFFtpc->GetFunction("fitFrac");
      Requires(prim_frac_TFF_TPC,"prim_frac_TFF_TPC");
      TpcTFFCorrection[iC] = (TH1F*) hResTFFtpc->Rebin(kNCustomPtBins,Form("TpcTFFCorrection_%d",iC),kCustomPtBins);
      TpcTFFCorrection[iC]->GetListOfFunctions()->Clear();

      TH1* hResRooFitTof = (TH1*)file_RooFit_TOF.Get(Form("%s/Results/hResTFF_%i",list_key->GetName(),iC));
      Requires(hResRooFitTof,"hResRooFitTof");
      TF1* prim_frac_RooFit_TOF = hResRooFitTof->GetFunction("fitFrac");
      Requires(prim_frac_RooFit_TOF,"prim_frac_RooFit__TOF");

      TH1* hResRooFitTpc = (TH1*)file_RooFit_TPC.Get(Form("%s/Results/hResTFF_%i",list_key->GetName(),iC));
      Requires(hResTFFtpc,"hResTFFtpc");
      TF1* prim_frac_RooFit_TPC = hResRooFitTpc->GetFunction("fitFrac");
      Requires(prim_frac_RooFit_TPC,"prim_frac_RooFit_TPC");

      TofFraction[iC] = new TH1F(Form("FractionTof_%i",iC),";p_{T} (GeV/c);RooFit / TFF",kNCustomPtBins,kCustomPtBins);
      TpcFraction[iC] = new TH1F(Form("FractionTpc_%i",iC),";p_{T} (GeV/c);RooFit / TFF",kNCustomPtBins,kCustomPtBins);

      for(int iB=1; iB<=kNCustomPtBins; iB++){
        TofFraction[iC]->SetBinContent(iB,prim_frac_RooFit_TOF->Eval(TofFraction[iC]->GetBinCenter(iB))/prim_frac_TFF_TOF->Eval(TofFraction[iC]->GetBinCenter(iB)));
        TofFraction[iC]->SetBinError(iB,0.);

        TpcFraction[iC]->SetBinContent(iB,prim_frac_RooFit_TPC->Eval(TpcFraction[iC]->GetBinCenter(iB))/prim_frac_TFF_TPC->Eval(TpcFraction[iC]->GetBinCenter(iB)));
        TpcFraction[iC]->SetBinError(iB,0.);
      }
      TofFraction[iC]->Write();
      TpcFraction[iC]->Write();
    }

    TCanvas CorrTof("corr_tof","corr_tof",3200,2400);
    CorrTof.DrawFrame(0.6,0.6,2.,1.,";#it{p}_{T} (GeV/#it{c});");
    TLegend legcorrtof(0.77,0.56,0.87,0.84);
    legcorrtof.SetBorderSize(0);
    legcorrtof.SetHeader("TOF (TFF)");
    for (int iC = 0; iC < kCentLength; ++iC) {
      plotting::SetHistStyle(TofTFFCorrection[iC],plotting::kSpectraColors[iC]);
      TofTFFCorrection[iC]->Draw("pex0same");
      legcorrtof.AddEntry(TofTFFCorrection[iC],Form("%4.0f - %2.0f %%",kCentLabels[iC][0],kCentLabels[iC][1]),"pe");
    }
    legcorrtof.Draw();
    CorrTof.Write();

    TCanvas CorrTpc("corr_tpc","corr_tpc",3200,2400);
    CorrTpc.DrawFrame(0.6,0.6,2.,1.,";#it{p}_{T} (GeV/#it{c});");
    TLegend legcorrtpc(0.77,0.56,0.87,0.84);
    legcorrtpc.SetBorderSize(0);
    legcorrtpc.SetHeader("TPC (TFF)");
    for (int iC = 0; iC < kCentLength; ++iC) {
      plotting::SetHistStyle(TpcTFFCorrection[iC],plotting::kSpectraColors[iC]);
      TpcTFFCorrection[iC]->Draw("pex0same");
      legcorrtpc.AddEntry(TpcTFFCorrection[iC],Form("%4.0f - %2.0f %%",kCentLabels[iC][0],kCentLabels[iC][1]),"pe");
    }
    legcorrtpc.Draw();
    CorrTpc.Write();

    TCanvas RatioTof("ratio_tof","ratio_tof",3200,2400);
    RatioTof.DrawFrame(0.6,0.95,1.6,1.15,";#it{p}_{T} (GeV/#it{c});RooFit/TFF");
    TLegend legtof(0.77,0.56,0.87,0.84);
    legtof.SetBorderSize(0);
    legtof.SetHeader("TOF");
    for (int iC = 0; iC < kCentLength; ++iC) {
      plotting::SetHistStyle(TofFraction[iC],plotting::kSpectraColors[iC]);
      TofFraction[iC]->Draw("psame");
      legtof.AddEntry(TofFraction[iC],Form("%4.0f - %2.0f %%",kCentLabels[iC][0],kCentLabels[iC][1]),"p");
    }
    legtof.Draw();
    RatioTof.Write();

    TCanvas RatioTpc("ratio_tpc","ratio_tpc",3200,2400);
    RatioTpc.DrawFrame(0.6,0.95,1.4,1.15,";#it{p}_{T} (GeV/#it{c});RooFit/TFF");
    TLegend legtpc(0.77,0.56,0.87,0.84);
    legtpc.SetBorderSize(0);
    legtpc.SetHeader("TPC");
    for (int iC = 0; iC < kCentLength; ++iC) {
      plotting::SetHistStyle(TpcFraction[iC],plotting::kSpectraColors[iC]);
      TpcFraction[iC]->Draw("psame");
      legtpc.AddEntry(TpcFraction[iC],Form("%4.0f - %2.0f %%",kCentLabels[iC][0],kCentLabels[iC][1]),"p");
    }
    legtpc.Draw();
    RatioTpc.Write();


  }
}
