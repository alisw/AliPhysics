#include "src/Common.h"
#include "TH2F.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"

constexpr int nCol = 3;
constexpr int nRow = 3;
constexpr double sx[2]{0.36,1.-2*sx[0]};
constexpr double sy[2]{0.38,1.-2*sy[0]};
constexpr double fx = sx[1];
constexpr double fy = sy[1];
double global_y = 1.;
double global_x = 0.;
double xAxisEdges[9] = {3.9,3.6,2.4,3.9,3.6,2.4,3.9,3.6,2.4};

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

    std::cout << global_x << "\t" << global_y - sy[mid] << "\t" << global_x + sx[center] << "\t" << global_y << std::endl;
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


void MakeRatioPlot() {
  gStyle->SetOptStat(0);

  TLine l;
  l.SetLineStyle(kDashed);
  l.SetLineColor(kBlack);

  TCanvas *cv = new TCanvas("cv","cv",800,600);
  auto pads = CreatePads(cv);

  TLatex text;
  text.SetTextFont(63);
  text.SetTextSize(18);
  pads[0]->cd();
  //text.DrawText(1.,2.20,"ALICE Preliminary");
  pads[1]->cd();
  text.DrawLatex(1.2,2.20,"#bf{pp, #sqrt{#it{s}} = 13 TeV}");
  pads[2]->cd();
  text.DrawLatex(0.6,2.20,"#bf{V0M Multiplicity Classes}");

  const string labels[9]{"0-1%","1-5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-70%","70-100%"};
  TFile input(kFinalOutput.data());
  int place_holder[9] = {0,3,6,1,4,7,2,5,8};
  TH1* stat[9] = {nullptr};
  TH1* syst[9] = {nullptr};
  for (int iC = 0; iC < 9; ++iC) {
    pads[place_holder[iC]]->cd();
    text.DrawLatex(0.5,1.7,Form("#bf{%s}",kRomanLabels[iC]));
    stat[iC] = (TH1*)input.Get(Form("ratio/%i/ratio_stat",iC));
    stat[iC]->SetDirectory(0);
    syst[iC] = (TH1*)input.Get(Form("ratio/%i/ratio_syst",iC));
    syst[iC]->SetDirectory(0);
    cout << iC  << "\t" << stat << "\t" << syst << endl;
    stat[iC]->Draw("esamex0");
    syst[iC]->Draw("e2same");
    l.DrawLine(0.3,1.,xAxisEdges[place_holder[iC]],1.);
  }

  cv->SaveAs(Form("%sdeuteron_ratio.pdf",kFiguresFolder.data()));

}
