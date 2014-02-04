static  int      myDarkRed  = TColor::GetColor(128,0,0);
static  int      myLightRed  = TColor::GetColor(128,0,0);
static  int      myBlue     = 9;
static  int      myGreen     = kGreen+3;
static  int      mySysErrColorOpp   = 17;
static  int      mySysErrColorSame    = 11;
// static  int      myv1FluctColorSame    = 29;
static  int      myv1FluctColorSame    = TColor::GetColor(207,206,232);
// static  int      myv1FluctColorSame    = TColor::GetColor(212,212,212);
static  int      myToneevColor    = myDarkRed;
//static  int      myMarkerSize    = 2.0;

// =======================================================================//
void SetPlotStyle() {
  // Set style which will affect all plots.
  
  gStyle->Reset();
  gStyle->SetStatColor(0);
  //gStyle->SetHistFillColor(17);
  // gStyle->SetOptitle(0);
  gStyle->SetOptStat("");
  // gStyle->SetOptStat(0);
  //gStyle->SetOptStat(111111);
  //gStyle->SetOptDate(1);
  // gStyle->SetPalette(8,0);  // (1,0)
  gStyle->SetPalette(1);  // (1,0)
  gStyle->SetDrawBorder(0);
  // gStyle->SetFillColor(0);  // kills palete ???
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  // gStyle->SetFillColor(0); // otherwize it affects Fill colors later
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameLineWidth(2);
  // gStyle->SetFrameFillStyle(4000);
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
  gStyle->SetTitleOffset(1.1,"y");
  gStyle->SetTitleFillColor(0);
  gStyle->SetLineWidth(2);
  gStyle->SetHistLineColor(1);
  gStyle->SetTextColor(1);
  gStyle->SetTitleTextColor(1);
  TGaxis::SetMaxDigits(4);
  gROOT->ForceStyle();
}
