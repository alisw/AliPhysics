////////////////////////////////////////////////////////////
//
//  Set Style for Flow papers plots
//
////////////////////////////////////////////////////////////
void setFlowStyle()
{
  gStyle->Reset();
  //  gStyle->SetOptTitle(0);
  //  gStyle->SetOptStat(0);
  // gStyle->SetOptDate(1);
  //  gStyle->SetPalette(8,0);  // (1,0)

  gStyle->SetPalette(1);  // (1,0)
  gStyle->SetDrawBorder(0);

  //  gStyle->SetFillColor(0);  // kills palette ???
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  //  gStyle->SetFillColor(0);  othewize it affects Fill colors later
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(2);
  //  gStyle->SetFrameFillStyle(4000);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(2);
  gStyle->SetFuncWidth(2);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.06,"xyz");
  gStyle->SetTitleOffset(1.3,"y");
  gStyle->SetTitleFillColor(0);
  gStyle->SetLineWidth(2);
  
  gStyle->SetHistLineColor(1);
  gStyle->SetTextColor(1);
  gStyle->SetTitleTextColor(1);
  TGaxis::SetMaxDigits(4);

  gROOT->ForceStyle();

}

