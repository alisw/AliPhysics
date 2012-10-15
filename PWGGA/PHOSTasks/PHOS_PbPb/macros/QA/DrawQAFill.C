#include "TCanvas.h" // needed for some reason.

TFile* file;
const char* prefixToName = "imgsFill/";
const char* appendToName = ".pdf";
const int kNCents = 1;

void Draw(const char* name, const char* options = "", double yFrom=0., double yTo=-1.)
{
  TH1* hist = ((TH1*)file->Get(name))->Clone();
  hist->GetXaxis()->SetTitle("Fill");

  if( yFrom < yTo )
    hist->GetYaxis()->SetRangeUser(yFrom, yTo);

  TCanvas* canv = new TCanvas;
  //canv->SetGrid();
  //hist->GetYaxis()->SetNdivisions(16);

  if( TString(options).Contains("LINFIT") )
    hist->Fit("pol0", "Q");

  hist->DrawCopy(options);

  canv->SaveAs(Form("%s%s%s", prefixToName, hist->GetName(), appendToName ));
  delete hist;
}

void DrawQAFill()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  file = TFile::Open("outputQAFill.root", "read");

  Draw("grVtxZ10Cent", "", 0.7, 1.);
  // Draw("grNCellsM1", "E");
  // Draw("grNCellsM2");
  // Draw("grNCellsM3");
  // Draw("grECluster", "", 0.5, 0.7);
  Draw("grNCluster", "", 0, 40);
  Draw("grNTracks0", "", 0 , 12000);
  // Draw("grNPhotAll_cen0", "", 0, 40);
  // Draw("grNPhotAllcore_cen0", "", 0, 40);
  // Draw("grNPhotAllwou_cen0", "", 0, 40);
  // Draw("grNPhotDisp_cen0", "", 0, 40);
  // Draw("grNPhotDisp2_cen0", "", 0, 40);
  // Draw("grNPhotDispwou_cen0", "", 0, 40);
  // Draw("grNPhotCPV_cen0", "", 0, 40);
  // Draw("grNPhotCPV2_cen0", "", 0, 40);
  // Draw("grNPhotBoth_cen0", "", 0, 40);
  // Draw("grEnAll_cen0", "", 0.4, 0.7);
  // Draw("grEnAllcore_cen0", "", 0.4, 0.7);
  // Draw("grEnAllwou_cen0", "", 0.4, 0.7);
  // Draw("grEnDisp_cen0", "", 0.4, 0.7);
  // Draw("grEnDisp2_cen0", "", 0.4, 0.7);
  // Draw("grEnDispcore_cen0", "", 0.4, 0.7);
  // Draw("grEnDispwou_cen0", "", 0.4, 0.7);
  // Draw("grEnCPV_cen0", "", 0.4, 0.7);
  // Draw("grEnCPVcore_cen0", "", 0.4, 0.7);
  // Draw("grEnCPV2_cen0", "", 0.4, 0.7);
  // Draw("grEnBoth_cen0", "", 0.4, 0.7);
  // Draw("grEnBothcore_cen0", "", 0.4, 0.7);


  Draw("grMPi0", "LINFIT", 0.13, 0.15);
  Draw("grWPi0", "LINFIT");
  Draw("grNPi0", "LINFIT");

  file->Close();
}
