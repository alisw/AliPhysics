// =====================================================================================
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
//static  int      myMarkerSize    = 2.5;
static  int      myMarkerSize    = 2.0;
void SetFlowStyle()
{
  // Set style which will affect all plots.

  // Set gStyle
  int font = 42;
  // From plain
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(1);
  gStyle->SetLegendBorderSize(1);
  //
  gStyle->SetDrawBorder(0);
  gStyle->SetTextFont(font);
  gStyle->SetStatFont(font);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);
  gStyle->SetTickLength(0.02,"y");
  gStyle->SetEndErrorSize(3);
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelFont(font,"xyz"); 
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetTitleFont(font,"xyz");  
  gStyle->SetTitleOffset(1.0,"xyz");  
  gStyle->SetTitleSize(0.06,"xyz");  
  gStyle->SetMarkerSize(1); 
  gStyle->SetPalette(1,0); 
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gROOT->ForceStyle();
  gStyle->SetFillStyle(1001);
}

//_______________________________________________________________//
void myPadSetUp(TPad *currentPad, float currentLeft=0.11, 
		float currentTop=0.04, float currentRight=0.04, 
		float currentBottom=0.15,
		Int_t gColor = 10){
  currentPad->SetFillColor(gColor);
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}

//_______________________________________________________________//
void DrawMarker(Double_t x, Double_t y, Int_t style, 
		Double_t size, Int_t color) {
  TMarker *m = new TMarker(x,y,style);
  m->SetMarkerSize(size);
  m->SetMarkerColor(color);
  m->Draw();
}

void myLegendSetUp(TLegend *currentLegend=0,float currentTextSize=0.07){
  currentLegend->SetTextFont(42);
  currentLegend->SetBorderSize(0);
  currentLegend->SetFillStyle(0);
  currentLegend->SetFillColor(0);
  currentLegend->SetMargin(0.25);
  currentLegend->SetTextSize(currentTextSize);
  currentLegend->SetEntrySeparation(0.5);
  return;
}

void myTGraphSetUp
(
  TGraphErrors *currentGraph=0,
  int myMarkerStyle=8,
  int myMarkerColor=1,
  float myMarkerSize=1,
  int myLineStyle=1,
  int myLineColor=1,
  float myLineWidth=1,
  int myFillStyle =1001,
  int myFillColor =1 
)
{
  currentGraph->SetMarkerStyle(myMarkerStyle);
  currentGraph->SetMarkerColor(myMarkerColor);
  currentGraph->SetMarkerSize(myMarkerSize);
  currentGraph->SetLineColor(myLineColor);
  currentGraph->SetLineStyle(myLineStyle);
  currentGraph->SetLineWidth(myLineWidth);
  currentGraph->SetFillStyle(myFillStyle);
  currentGraph->SetFillColor(myFillColor);
//   currentGraph->Set();
}

void myTGraphSetUp_Asym
(
  TGraphAsymmErrors *currentGraph=0,
  int myMarkerStyle=8,
  int myMarkerColor=1,
  float myMarkerSize=1,
  int myLineStyle=1,
  int myLineColor=1,
  float myLineWidth=1,
  int myFillStyle =1001,
  int myFillColor =1 
)
{
  currentGraph->SetMarkerStyle(myMarkerStyle);
  currentGraph->SetMarkerColor(myMarkerColor);
  currentGraph->SetMarkerSize(myMarkerSize);
  currentGraph->SetLineColor(myLineColor);
  currentGraph->SetLineStyle(myLineStyle);
  currentGraph->SetLineWidth(myLineWidth);
  currentGraph->SetFillStyle(myFillStyle);
  currentGraph->SetFillColor(myFillColor);
//   currentGraph->Set();
}

void myTGraphSetUpPr
(
  TProfile *currentGraph=0,
  int myMarkerStyle=8,
  int myMarkerColor=1,
  float myMarkerSize=1,
  int myLineStyle=1,
  int myLineColor=1,
  float myLineWidth=1,
  int myFillStyle =1001,
  int myFillColor =1 
)
{
  currentGraph->SetMarkerStyle(myMarkerStyle);
  currentGraph->SetMarkerColor(myMarkerColor);
  currentGraph->SetMarkerSize(myMarkerSize);
  currentGraph->SetLineColor(myLineColor);
  currentGraph->SetLineStyle(myLineStyle);
  currentGraph->SetLineWidth(myLineWidth);
  currentGraph->SetFillStyle(myFillStyle);
  currentGraph->SetFillColor(myFillColor);
//   currentGraph->Set();
}

void myTGraphSetUpTH
(
  TH1F *currentGraph=0,
  int myMarkerStyle=8,
  int myMarkerColor=1,
  float myMarkerSize=1,
  int myLineStyle=1,
  int myLineColor=1,
  float myLineWidth=1,
  int myFillStyle =1001,
  int myFillColor =1 
)
{
  currentGraph->SetMarkerStyle(myMarkerStyle);
  currentGraph->SetMarkerColor(myMarkerColor);
  currentGraph->SetMarkerSize(myMarkerSize);
  currentGraph->SetLineColor(myLineColor);
  currentGraph->SetLineStyle(myLineStyle);
  currentGraph->SetLineWidth(myLineWidth);
  currentGraph->SetFillStyle(myFillStyle);
  currentGraph->SetFillColor(myFillColor);
//   currentGraph->Set();
}

void zeroTGraphHorisontalErrors(TGraphErrors *currentGraph=0)
{
  for (int i)
  return;
}

void ShiftAlongXaxis_TGraphErrors(TGraphErrors *ge, Double_t shift)
{
 if(!ge)
 {
  printf("\n WARNING: ge is NULL in ShiftAlongXaxis() !!!! \n\n");
  return;
 }
 Int_t nPoints = ge->GetN();
 Double_t x = 0.;
 Double_t y = 0.;
 for(Int_t p=0;p<nPoints;p++)
 {
  ge->GetPoint(p,x,y);
  x+=shift;
  ge->SetPoint(p,x,y);
 }
}

void ShiftAlongXaxis_TGraphAsymmErrors(TGraphAsymmErrors *ge, Double_t shift)
{
 if(!ge)
 {
  printf("\n WARNING: ge is NULL in ShiftAlongXaxis() !!!! \n\n");
  return;
 }
 Int_t nPoints = ge->GetN();
 Double_t x = 0.;
 Double_t y = 0.;
 for(Int_t p=0;p<nPoints;p++)
 {
  ge->GetPoint(p,x,y);
  x+=shift;
  ge->SetPoint(p,x,y);
 }
}
