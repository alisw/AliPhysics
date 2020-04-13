#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include "src/Plotting.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TPaveText.h"

constexpr int nCol = 3;
constexpr int nRow = 3;
constexpr double sx[2]{0.36,1.-2*sx[0]};
constexpr double sy[2]{0.38,1.-2*sy[0]};
constexpr double fx = sx[1];
constexpr double fy = sy[1];
double global_y = 1.;
double global_x = 0.;
double xAxisEdges[9] = {4.2,3.7,3.2,4.2,3.7,3.2,4.2,3.7,3.2};
int place_holder[9] = {0,3,6,1,4,7,2,5,8};

std::array<TPad*,9> CreatePads(TCanvas* &cv);

void RatioErrorMakerFromFile(TH1F* hRatio, const char* file_name, const char* path_num, const char* path_den);

void Ratio(){

  gStyle->SetOptStat(0);

  TFile input_file(kSpectraOutput.data());
  TFile output_file(kRatioOutput.data(),"recreate");

  TCanvas *cRatioAll = new TCanvas("cv","cv",800,600);
  auto pads = CreatePads(cRatioAll);  

  TLatex text;
  text.SetTextFont(63);
  text.SetTextSize(18);
  pads[0]->cd();
  //text.DrawText(1.,2.20,"This work");
  pads[1]->cd();
  text.DrawLatex(1.2,2.20,"#bf{pp, #sqrt{#it{s}} = 13 TeV}");
  pads[2]->cd();
  text.DrawLatex(0.6,2.20,"#bf{V0M Multiplicity Classes}");
  TLine l;
  l.SetLineStyle(kDashed);
  l.SetLineColor(kBlack);

  for (int iC = 0; iC < kCentLength; ++iC) {
    TDirectory* cent_dir = output_file.mkdir(std::to_string(iC).data());
    TH1F* hSpectraRef[2] = {nullptr};
    TDirectory* ref_dir = cent_dir->mkdir("ref");
    string basepath = kFilterListNames + "/%s/" + std::to_string(iC) + "/Joined/JoinedSpectra%c" + std::to_string(iC);
    hSpectraRef[0] = (TH1F*)input_file.Get(Form(basepath.data(),kNames[0].data(),kLetter[0]));
    utils::Requires(hSpectraRef[0],"hSpectraRef[0]");
    hSpectraRef[1] = (TH1F*)input_file.Get(Form(basepath.data(),kNames[1].data(),kLetter[1]));
    utils::Requires(hSpectraRef[1],"hSpectraRef[1]");
    TH1F* hRatioRef = (TH1F*)hSpectraRef[1]->Clone("hRatio_ref");
    hRatioRef->Divide(hSpectraRef[0]);
    hRatioRef->SetTitle("ref;#it{p}_{T} (GeV/#it{c});#bar{d}/d");
    TH1F* hRatioSyst = (TH1F*)hRatioRef->Clone("hRatioSyst");
    hRatioSyst->Reset();
    ref_dir->cd();
    plotting::SetHistStyle(hSpectraRef[0],kBlack);
    plotting::SetHistStyle(hSpectraRef[1],kBlack);
    plotting::SetHistStyle(hRatioRef,kBlack);
    hSpectraRef[0]->Write();
    hSpectraRef[1]->Write();
    hRatioRef->Write();

    // Bin-counting systematics   
    std::string width_sys_path = kFilterListNames + "/%s" + Form("/Systematic/C_%i/hWidenRangeSystJoined",iC) + "%c" + std::to_string(iC);
    TH1F* hRatioSyst_width = new TH1F("hRatioSyst_width",";#it{p}_{T} (GeV/#it{c});",kNPtBins,kPtBins);
    RatioErrorMakerFromFile(hRatioSyst_width,kSignalOutput.data(),Form(width_sys_path.data(),kNames[0].data(),kLetter[0]),Form(width_sys_path.data(),kNames[1].data(),kLetter[1]));
    utils::SmoothInRange(hRatioSyst_width,kPtRange[0]+0.01,kCentPtLimits[iC]-0.01);
    cent_dir->cd();
    hRatioSyst_width->SetLineColor(plotting::kHighContrastColors[1]);
    hRatioSyst_width->Write();

    std::string shift_sys_path = kFilterListNames + "/%s" + Form("/Systematic/C_%i/hShiftRangeSystJoined",iC) + "%c" + std::to_string(iC);
    TH1F* hRatioSyst_shift = new TH1F("hRatioSyst_shift",";#it{p}_{T} (GeV/#it{c});",kNPtBins,kPtBins);
    RatioErrorMakerFromFile(hRatioSyst_shift,kSignalOutput.data(),Form(shift_sys_path.data(),kNames[0].data(),kLetter[0]),Form(shift_sys_path.data(),kNames[1].data(),kLetter[1]));
    utils::SmoothInRange(hRatioSyst_shift,kPtRange[0]+0.01,kCentPtLimits[iC]-0.01);
    cent_dir->cd();
    hRatioSyst_shift->SetLineColor(plotting::kHighContrastColors[5]);
    hRatioSyst_shift->Write();

    // Material Budget
    string mat_sys_path = "deuterons%ctpc";
    TH1F* hRatioSyst_mat = new TH1F("hRatioSyst_mat",";#it{p}_{T} (GeV/#it{c});",kNPtBins,kPtBins);
    RatioErrorMakerFromFile(hRatioSyst_mat,kMaterialOutput.data(),Form(mat_sys_path.data(),kLetter[0]),Form(mat_sys_path.data(),kLetter[1]));
    utils::SmoothInRange(hRatioSyst_mat,kPtRange[0]+0.01,kCentPtLimits[iC]-0.01);
    cent_dir->cd();
    hRatioSyst_mat->SetLineColor(plotting::kHighContrastColors[2]);
    hRatioSyst_mat->Write();

    // Hadronic interaction
    TH1F* hHadSyst[2] = {nullptr};
    for(int iS=0; iS<2; iS++){
      hHadSyst[iS] = new TH1F(Form("hHadSyst_%c",kLetter[iS]),";#it{p}_{T} (GeV/#it{c});",kNPtBins,kPtBins);
      for(int iBin=1; iBin<=kNPtBins; iBin++){
        if(hHadSyst[iS]->GetBinCenter(iBin)>kCentPtLimits[iC]) break;
        hHadSyst[iS]->SetBinContent(iBin,kAbsSyst[iS]);
      }
    }
    TH1F* hRatioSyst_had = new TH1F("hRatioSyst_had",";#it{p}_{T} (GeV/#it{c});",kNPtBins,kPtBins);
    for(int iBin=1; iBin<=kNPtBins; iBin++){
      auto value = TMath::Sqrt(Sq(hHadSyst[0]->GetBinContent(iBin))+ Sq(hHadSyst[1]->GetBinContent(iBin)));
      hRatioSyst_had->SetBinContent(iBin,value);
    }
    utils::SmoothInRange(hRatioSyst_had,kPtRange[0]+0.01,kCentPtLimits[iC]-0.01);
    cent_dir->cd();
    hRatioSyst_had->SetLineColor(plotting::kHighContrastColors[3]);
    hRatioSyst_had->Write();


    // Selection systematics
    TCanvas* cRatioSyst = new TCanvas(Form("cRatioSyst%d",iC),Form("cRatioSyst%d",iC));
    TH1F* hRatioSyst_dcaz = (TH1F*)hRatioRef->Clone("hRatioSyst_dcaz");
    hRatioSyst_dcaz->Reset();
    hRatioSyst_dcaz->SetLineColor(plotting::kHighContrastColors[0]);
    TH1F* hRatioSyst_dcaxy = (TH1F*)hRatioRef->Clone("hRatioSyst_dcaxy");
    hRatioSyst_dcaxy->Reset();
    hRatioSyst_dcaxy->SetLineColor(plotting::kHighContrastColors[8]);
    TH1F* hRatioSyst_pid = (TH1F*)hRatioRef->Clone("hRatioSyst_pid");
    hRatioSyst_pid->SetLineColor(plotting::kHighContrastColors[6]);
    hRatioSyst_pid->Reset();
    TH1F* hRatioSyst_tpc = (TH1F*)hRatioRef->Clone("hRatioSyst_tpc");
    hRatioSyst_tpc->SetLineColor(plotting::kHighContrastColors[7]);
    hRatioSyst_tpc->Reset();
    TH1F* hRatioSyst_cuts = (TH1F*)hRatioRef->Clone("hRatioSyst_cuts");
    hRatioSyst_cuts->Reset();
    hRatioSyst_cuts->SetLineColor(plotting::kHighContrastColors[0]);
    for (auto& syst : kCutNames) {
      const char* kCutName = syst.first.data();
      auto& kCutString = syst.first;
      basepath = kFilterListNames + kCutName + "%i" + "/%s/" + std::to_string(iC) + "/Joined/JoinedSpectra%c" + std::to_string(iC);
      TDirectory* cut_dir = cent_dir->mkdir(kCutName);
      const int kNvariations = static_cast<int>(syst.second.size());
      TH1F* hSpectraVar[2][kNvariations];
      TH1F* hRatioVar[kNvariations];
      TCanvas* cRatioCut = new TCanvas(Form("ratio_%s_c%d",kCutName,iC));
      cRatioCut->DrawFrame(0.4,.1,4.0,2.1,";#it{p}_{T} (GeV/#it{c});#bar{d}/d");
      TLine *line_one = new TLine(0.4,1.,4.0,1.);
      line_one->SetLineColor(kBlack);
      line_one->SetLineStyle(2);
      line_one->Draw();
      TLegend* leg = new TLegend(0.16,0.65,0.36,0.85);
      hRatioRef->Draw("PESAME");
      leg->AddEntry(hRatioRef);
      /// Filling the histograms for each cuts.
      for (int iCut = 0; iCut < kNvariations; iCut++) {
        hSpectraVar[0][iCut] = (TH1F*)input_file.Get(Form(basepath.data(),iCut,kNames[0].data(),kLetter[0]));
        utils::Requires(hSpectraVar[0][iCut],"hSpectraVar[0]");
        hSpectraVar[1][iCut] = (TH1F*)input_file.Get(Form(basepath.data(),iCut,kNames[1].data(),kLetter[1]));
        utils::Requires(hSpectraVar[1][iCut],"hSpectraVar[1]");
        hRatioVar[iCut] = (TH1F*)hSpectraVar[1][iCut]->Clone(Form("hRatio_%s_%d",kCutName,iCut));
        hRatioVar[iCut]->Divide(hSpectraVar[0][iCut]);
        hRatioVar[iCut]->SetTitle(Form("%s_%d;#it{p}_{T} (GeV/#it{c});#bar{d}/d",kCutName,iCut));
        TDirectory* var_dir = cut_dir->mkdir(Form("%s_%d",kCutName,iCut));
        var_dir->cd();
        plotting::SetHistStyle(hSpectraVar[0][iCut],plotting::kSpectraColors[iCut]);
        plotting::SetHistStyle(hSpectraVar[1][iCut],plotting::kSpectraColors[iCut]);
        plotting::SetHistStyle(hRatioVar[iCut],plotting::kSpectraColors[iCut]);
        hSpectraVar[0][iCut]->Write();
        hSpectraVar[1][iCut]->Write();
        hRatioVar[iCut]->Write();
        cRatioCut->cd();
        hRatioVar[iCut]->Draw("PESAME");
        leg->AddEntry(hRatioVar[iCut]);
      }
      cut_dir->cd();
      cRatioCut->cd();
      leg->SetBorderSize(0);
      leg->Draw();
      cRatioCut->Write();
      /// computing the systematic uncertainty as RMS
      for(int iBin = 1; iBin <= kNPtBins; iBin++){
        if (hRatioRef->GetXaxis()->GetBinCenter(iBin) < kPtRange[0]||
              hRatioRef->GetXaxis()->GetBinCenter(iBin) > kPtRange[1])
            continue;
        if (hRatioRef->GetXaxis()->GetBinCenter(iBin)>kCentPtLimits[iC]) continue;
        const float value_ref = hRatioRef->GetBinContent(iBin);
        const float err_ref = hRatioRef->GetBinError(iBin);
        vector<float> value_vec = {value_ref};
        vector<float> error_vec = {err_ref};
        for(int iCut = 0; iCut < kNvariations; iCut++){
          const float value_var = hRatioVar[iCut]->GetBinContent(iBin);
          const float err_var = hRatioVar[iCut]->GetBinError(iBin);
          const float z_test = utils::zTest(value_ref,err_ref,value_var,err_var);
          if(TMath::Abs(z_test)>2){
            value_vec.push_back(value_var);
            error_vec.push_back(err_var);
          }
        }
        float bin_syst = TMath::RMS(value_vec.begin(),value_vec.end())/value_ref;
        if(kCutString=="dcaz"){
          hRatioSyst_dcaz->SetBinContent(iBin,bin_syst);
        }
        if(kCutString=="dcaxy"){
          hRatioSyst_dcaxy->SetBinContent(iBin,bin_syst);
        }
        else if(kCutString=="pid"){
          hRatioSyst_pid->SetBinContent(iBin,bin_syst);
        }
        else if(kCutString=="tpc"){
          hRatioSyst_tpc->SetBinContent(iBin,bin_syst);
        }
      }
      if(kCutString=="dcaz"){
        utils::SmoothInRange(hRatioSyst_dcaz,kPtRange[0]+0.01,kCentPtLimits[iC]-0.01);
        hRatioSyst_dcaz->Write();
      }
      if(kCutString=="dcaxy"){
        utils::SmoothInRange(hRatioSyst_dcaxy,kPtRange[0]+0.01,kCentPtLimits[iC]-0.01);
        hRatioSyst_dcaxy->Write();
      }
      else if(kCutString=="pid"){
        utils::SmoothInRange(hRatioSyst_pid,kPtRange[0]+0.01,kCentPtLimits[iC]-0.01);
        hRatioSyst_pid->Write();
      }
      else if(kCutString=="tpc"){
        utils::SmoothInRange(hRatioSyst_pid,kPtRange[0]+0.01,kCentPtLimits[iC]-0.01);
        hRatioSyst_tpc->Write();
      }
    }
    for(int iBin=1; iBin<=kNPtBins; iBin++){
      if (hRatioRef->GetXaxis()->GetBinCenter(iBin) < kPtRange[0]||
              hRatioRef->GetXaxis()->GetBinCenter(iBin) > kPtRange[1])
            continue;
      if (hRatioRef->GetXaxis()->GetBinCenter(iBin)>kCentPtLimits[iC]) continue;
      float dcaz_err = hRatioSyst_dcaz->GetBinContent(iBin);
      float dcaxy_err = hRatioSyst_dcaxy->GetBinContent(iBin);
      float pid_err = hRatioSyst_pid->GetBinContent(iBin);
      float tpc_err = hRatioSyst_tpc->GetBinContent(iBin);
      hRatioSyst_cuts->SetBinContent(iBin,TMath::Sqrt(dcaz_err*dcaz_err+dcaxy_err*dcaxy_err+pid_err*pid_err+tpc_err*tpc_err));
    }
    cent_dir->cd();
    hRatioSyst_cuts->Write();
    
    // Total systematic uncertainty
    for(int iBin=1; iBin<=kNPtBins; iBin++){
      float tot = sqrt(
              hRatioSyst_cuts->GetBinContent(iBin) * hRatioSyst_cuts->GetBinContent(iBin) +
              hRatioSyst_mat->GetBinContent(iBin) * hRatioSyst_mat->GetBinContent(iBin) +
              hRatioSyst_had->GetBinContent(iBin) * hRatioSyst_had->GetBinContent(iBin) +
              hRatioSyst_width->GetBinContent(iBin) * hRatioSyst_width->GetBinContent(iBin)+
              hRatioSyst_shift->GetBinContent(iBin) * hRatioSyst_shift->GetBinContent(iBin)
              );
      hRatioSyst->SetBinContent(iBin,tot);
    }
    hRatioSyst->SetLineColor(plotting::kHighContrastColors[4]);

    cRatioSyst->cd();
    cRatioSyst->DrawFrame(0.5,0.,kCentPtLimits[iC]+0.1,0.3,";#it{p}_{T} (GeV/#it{c}); Systematics uncertainties");
    hRatioSyst_cuts->Draw("SAME");
    hRatioSyst_width->Draw("SAME");
    hRatioSyst_shift->Draw("SAME");
    hRatioSyst_mat->Draw("SAME");
    hRatioSyst_had->Draw("SAME");
    hRatioSyst->Draw("same");
    TLegend* leg = new TLegend(0.6,0.56,0.89,0.84);
    leg->AddEntry(hRatioSyst_cuts,"PID and cuts","l");
    leg->AddEntry(hRatioSyst_width,"Range widening","l");
    leg->AddEntry(hRatioSyst_shift,"Range shifting","l");
    leg->AddEntry(hRatioSyst_mat,"Material budget","l");
    leg->AddEntry(hRatioSyst_had,"Hadronic interaction","l");
    leg->AddEntry(hRatioSyst,"total","l");
    leg->Draw();
    cRatioSyst->Write();

    TCanvas* cRatio = new TCanvas(Form("ratio_%d",iC));
    cRatio->DrawFrame(0.4,0.1,4.0,2.1,";#it{p}_{T} (GeV/#it{c});#bar{d}/d");
    plotting::SetHistStyle(hRatioRef,plotting::kSpectraColors[iC]);
    hRatioRef->SetTitle("");
    TH1F* hRatioSystFinal = (TH1F*)hRatioRef->Clone(";#it{p}_{T} (GeV/#it{c});#bar{d}/d");
    for(int iBin=1; iBin<=kNPtBins; iBin++){
      hRatioSystFinal->SetBinError(iBin,hRatioRef->GetBinContent(iBin) * hRatioSyst->GetBinContent(iBin));
    }
    plotting::SetHistStyle(hRatioSystFinal,plotting::kSpectraColors[iC]);
    TH1F* hRatioTot = (TH1F*)hRatioSystFinal->Clone("hRatioTot");
    for(int iBin=1; iBin<=kNPtBins; iBin++){
      hRatioTot->SetBinError(iBin,TMath::Sqrt(hRatioSystFinal->GetBinError(iBin)*hRatioSystFinal->GetBinError(iBin)+hRatioRef->GetBinError(iBin)*hRatioRef->GetBinError(iBin)));
    }
    hRatioTot->SetFillStyle(3004);
    hRatioTot->SetFillColor(plotting::kSpectraColors[iC]);
    hRatioTot->SetLineColor(plotting::kSpectraColors[iC]);
    TF1* func = new TF1(Form("ratiopol_%d",iC),"pol0",0.6,kCentPtLimits[iC]);
    func->SetLineWidth(2);
    func->SetLineColor(plotting::kSpectraColors[iC]);
    hRatioTot->Fit(Form("ratiopol_%d",iC),"QR0");
    hRatioTot->Draw("e2same");
    hRatioRef->Draw("esamex0");
    hRatioSystFinal->Draw("e2same");
    func->Draw("same");
    TLatex lat;
    lat.DrawLatex(1.5,1.7,Form("p0: %f #pm %f",func->GetParameter(0),func->GetParError(0)));
    lat.DrawLatex(1.5,1.6,Form("#chi^{2} / NDF: %.2f / %d",func->GetChisquare(),func->GetNDF()));
    TLine* line = new TLine(0.4,1.,4.0,1.);
    line->SetLineColor(kBlack);
    line->SetLineStyle(kDashed);
    line->Draw("same");
    cRatio->Write();
    if(iC==9) {
      TCanvas* cRatioMB = new TCanvas("RatioMB","RatioMB");
      cRatioMB->DrawFrame(0.4,0.1,4.0,2.1,";#it{p}_{T} (GeV/#it{c});#bar{d}/d");
      hRatioSystFinal->Draw("e2same");
      hRatioRef->Draw("esamex0");
      line->Draw("same");
      TPaveText paveMB(0.16,0.73,0.40,0.84,"blNDC");
      paveMB.SetBorderSize(0);
      paveMB.SetFillColor(0);
      paveMB.SetTextSize(0.03);
      paveMB.AddText("This work");
      paveMB.AddText("#bf{pp, #sqrt{#it{s}} = 13 TeV}");
      paveMB.Draw();
      cRatioMB->Write();
      break;
    }
    pads[place_holder[iC]]->cd();
    text.DrawLatex(0.5,1.7,Form("#bf{%s}",kRomanLabels[iC]));
    hRatioRef->Draw("esamex0");
    hRatioSystFinal->Draw("e2same");
    l.DrawLine(0.3,1.,xAxisEdges[place_holder[iC]],1.);

  }
  output_file.cd();
  cRatioAll->Write();
  cRatioAll->SaveAs(Form("%sdeuteron_ratio.pdf",kFiguresFolder.data()));

}

void RatioErrorMakerFromFile(TH1F* hRatio, const char* file_name, const char* path_num, const char* path_den){
  TFile file(file_name);
  TH1F* hHist[2] = {nullptr};
  hHist[0] = (TH1F*)file.Get(path_num);
  utils::Requires(path_num);
  hHist[1] = (TH1F*)file.Get(path_den);
  for(int iBin=1; iBin<=hRatio->GetXaxis()->GetNbins(); iBin++){
    float bin_center = hRatio->GetXaxis()->GetBinCenter(iBin);
    int bin_num = hHist[0]->FindBin(bin_center);
    int bin_den = hHist[0]->FindBin(bin_center);
    auto value = TMath::Sqrt(Sq(hHist[0]->GetBinContent(bin_num)) + Sq(hHist[1]->GetBinContent(bin_den)));
    hRatio->SetBinContent(iBin,value);
  }
}

std::array<TPad*,9> CreatePads(TCanvas* &cv)
{
  if (!cv) cv = new TCanvas;
  std::array<TPad*,9> pads{nullptr};

  for (int iP = 0; iP < 9; ++iP) {
    const int col = iP % nRow;
    const int row = iP / nRow;

    const int top = (row == 0);
    const int bot = (row == (nRow - 1));
    const int mid = !(top || bot);
    const int left = col==0;
    const int right = col==(nCol-1);
    const int center = !(left || right);

    pads[iP] = new TPad(Form("ratio%i",iP),"",global_x, std::abs(global_y - sy[mid]), global_x + sx[center], global_y);
    if (col == nCol - 1) global_y -= sy[mid];
    global_x = (right) ? 0. : global_x + sx[center];

    pads[iP]->SetRightMargin(right * (1 - fx / sx[center]));
    pads[iP]->SetLeftMargin(left * (1 - fx / sx[center]));
    pads[iP]->SetTopMargin(top * (1 - fy / sy[mid]));
    pads[iP]->SetBottomMargin(bot * (1 - fy / sy[mid]));

    cv->cd();
    pads[iP]->Draw();
    pads[iP]->cd();
    // if (iP == nCol - 1) {
    //   TLine border;
    //   border.SetLineColor(kBlack);
    //   border.DrawLineNDC(0.,0.,0.,fy / sy[mid]);
    //   continue;
    // }
    TH2F *rframe = new TH2F(Form("rframe%i",iP),";#it{p}_{T} (GeV/#it{c});#bar{d}/d;",100,0.3,xAxisEdges[iP],100,0.1,2.1);

    rframe->GetYaxis()->CenterTitle();
    rframe->GetYaxis()->SetTickLength(0.012 / sx[center]);
    rframe->GetYaxis()->SetTitleSize(20);
    rframe->GetYaxis()->SetTitleFont((!col) * 43);
    rframe->GetYaxis()->SetTitleOffset(2.2);
    rframe->GetYaxis()->SetLabelOffset(0.01);
    rframe->GetYaxis()->SetNdivisions(505);
    rframe->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    rframe->GetYaxis()->SetLabelSize((!col) * 15);

    rframe->GetXaxis()->CenterTitle();
    rframe->GetXaxis()->SetTickLength(0.012 / sy[mid]);
    rframe->GetXaxis()->SetTitleSize((bot && center) * 20);
    rframe->GetXaxis()->SetTitleFont(43);
    rframe->GetXaxis()->SetTitleOffset(4);
    if(iP%3==2)rframe->GetXaxis()->SetNdivisions(505);
    else rframe->GetXaxis()->SetNdivisions(510);
    rframe->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    rframe->GetXaxis()->SetLabelSize(bot * 15);

    rframe->Draw("col");
  }
  return pads;
}

