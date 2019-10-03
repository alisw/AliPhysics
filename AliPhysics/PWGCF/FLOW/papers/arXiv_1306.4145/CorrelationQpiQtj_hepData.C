const int nTerms = 4;

TGraphErrors*ResolutionXY[nTerms];
TString title[nTerms]={"#LTQ_{x}^{p}Q_{x}^{t}#GT","#LTQ_{y}^{p}Q_{y}^{t}#GT","#LTQ_{x}^{p}Q_{y}^{t}#GT","#LTQ_{y}^{p}Q_{x}^{t}#GT"};
TString name[nTerms]={"Qpx_Qtx","Qpy_Qty","Qpx_Qty","Qpy_Qtx"};
int color[nTerms]={kBlue,kMagenta+2,kBlue+2,kGreen+2};
int symbol[nTerms]={kFullSquare,kFullDiamond,kOpenTriangleUp,kOpenTriangleDown};
float myMarkerSize=1.5;
float markerSize[nTerms]={myMarkerSize,1.5*myMarkerSize,myMarkerSize,myMarkerSize};
float mainFont = 20;
int fillStyle = 1001;
float myShift =0.;
float shift[nTerms]={-myShift,myShift,-myShift,myShift};
TString centrality_name[9]={"0-5","5-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80"};
void CorrelationQpiQtj_hepData(bool rWrite = false, bool dumpHepData=false,TString dataFileName="ALICE_v1_arxiv_1306_4145.root")
{
  TGaxis::SetMaxDigits(3);
  myOptions();
  gROOT->ForceStyle();
  TCanvas *myCan = new TCanvas("myCan","Resolution",800,600);
  TPad *myPad1 = new TPad("myPad1","myPad1",0.,0,1,1);
  myPadSetUp(myPad1,0.08, 0.05, 0.01, 0.1);
  myPad1->Draw();
  myPad1->cd();
  TFile *dataFile = TFile::Open(dataFileName,"READ");
  ResolutionXY[0] = (TGraphErrors*)dataFile->Get("Qpx_Qtx");
  ResolutionXY[1] = (TGraphErrors*)dataFile->Get("Qpy_Qty");
  ResolutionXY[2] = (TGraphErrors*)dataFile->Get("Qpx_Qty");
  ResolutionXY[3] = (TGraphErrors*)dataFile->Get("Qpy_Qtx");
  
  TH1F* myBlankHisto = new TH1F("myBlankHisto","#LTQ_{i}^{p}Q_{j}^{t}#GT",10,0,82);
  myHistoSetUp(myBlankHisto,"centrality (%)","#LTQ_{r}^{p}Q_{r'}^{t}#GT",-1.05e-2,2.5e-3, 1, 2, 315, 505);
  TLegend* myLegend = new TLegend(0.35,0.45,0.6,0.7);
  myLegendSetUp(myLegend,mainFont);
  myBlankHisto->Draw();
  for (int i=0; i<nTerms; i++)
  {
    myTGraphSetUp(ResolutionXY[i],symbol[i],color[i],markerSize[i],1,color[i],2,fillStyle,color[i]);
    TString draw_opt = "p,eZ,same";
    myLegend->AddEntry(ResolutionXY[i],title[i],"P");
    ResolutionXY[i]->Draw(draw_opt);
    if (dumpHepData)
    {
      cout<< name[i] <<endl; cout <<"x" << "\t" <<"xlow" << "\t"<< "xhigh"<< "\t"<< "y"<< "\t\t"<<"dy+"<< "\t\t"<<"dy-"<< "\t\t"<<"dy+"<< "\t"<<"dy-"<<endl;
      Double_t x = 0.;
      Double_t y = 0.;
      for (int j=0; j< ResolutionXY[i]->GetN(); j++)
      {
	ResolutionXY[i]->GetPoint(j,x,y);
	cout<< x<< "\t" << x-ResolutionXY[i]->GetErrorXlow(j)<< "\t" << x+ResolutionXY[i]->GetErrorXhigh(j)<< "\t" << y << "\t"<< ResolutionXY[i]->GetErrorY(j)<< "\t"<< ResolutionXY[i]->GetErrorY(j)<<"\tN/A\tN/A" <<endl;
      }
    }	
    ShiftAlongXaxis(ResolutionXY[i], shift[i]);
  }
  myLegend->Draw();
  TLatex *myText = new TLatex();
  myText->SetNDC();
  myText->SetTextSize(mainFont);
  myText->DrawLatex(0.12,0.9,"ALICE Pb-Pb #sqrt{s_{NN}}=2.76TeV");
  myText->DrawLatex(0.57,0.9,"Q^{p,t}_{x,y} from Neutron ZDCs (|#eta|>8.78)");
  
  TString fileName="Correlation_QpQt";
  myCan->Update();
  if (rWrite)  
  {
    myCan->SaveAs(fileName+".png");
    myCan->SaveAs(fileName+".eps");
    myCan->SaveAs(fileName+".pdf");
  }
}

void myLegendSetUp(TLegend *currentLegend=0,float currentTextSize=0.07){
  currentLegend->SetBorderSize(0);
  currentLegend->SetFillStyle(0);
  currentLegend->SetFillColor(0);
  currentLegend->SetMargin(0.25);
  currentLegend->SetTextSize(currentTextSize);
  currentLegend->SetEntrySeparation(0.5);
  return;
}
//   myTGraphSetUp(v1etaStat_cent10_20,kFullTriangleUp,kRed+2,myMarkerSize,1,kRed+2,2,fillStyle,kRed+2);

void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}

void myGraphSetUp(TGraphErrors *currentGraph=0, Float_t currentMarkerSize = 1.0,
		  int currentMarkerStyle=20, int currentMarkerColor=0,
		  int currentLineStyle=1, int currentLineColor=0)
{
  currentGraph->SetMarkerSize(currentMarkerSize);
  currentGraph->SetMarkerStyle(currentMarkerStyle);
  currentGraph->SetMarkerColor(currentMarkerColor);
  currentGraph->SetLineStyle(currentLineStyle);
  currentGraph->SetLineColor(currentLineColor);
  return;
}

void myOptions(Int_t lStat=0){
  // Set gStyle
  int font = 43;
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
  gStyle->SetLegendFont(font);
  gStyle->SetStatFontSize(20);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);
  gStyle->SetTickLength(0.02,"y");
  gStyle->SetEndErrorSize(3);
  gStyle->SetLabelSize(20,"xyz");
  gStyle->SetLabelFont(font,"xyz"); 
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetTitleFont(font,"xyz");
  gStyle->SetTitleOffset(.8,"y");  
  gStyle->SetTitleOffset(1.,"xz");  
  gStyle->SetTitleSize(24,"x");  
  gStyle->SetTitleSize(24,"y");  
  gStyle->SetMarkerSize(1); 
  gStyle->SetPalette(1,0); 
  
  if (lStat){
    gStyle->SetOptTitle(1);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
  }
  else {
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
  }
}

TGraphErrors* makeGraphH1(TH1* hist, TString name="")
{
  name.ReplaceAll(" ","");
  Int_t nbins = hist->GetNbinsX();
  Double_t* x = new Double_t[nbins];
  Double_t* y = new Double_t[nbins];
  Double_t* xerr = new Double_t[nbins]; 
  Double_t* yerr = new Double_t[nbins];
  Int_t n=0;
  for (Int_t i=0; i<nbins; i++)
  {
    //if (hist->GetBinContent(i+1)==0.0) continue;
    
    x[n] = hist->GetXaxis()->GetBinCenter(i+1);
    y[n] = hist->GetBinContent(i+1);
    xerr[n] = 0.0;
    yerr[n] = hist->GetBinError(i+1);
    n++;
  }
  TGraphErrors* gr = new TGraphErrors(n,x,y,xerr,yerr);
  delete [] x;
  delete [] y;
  delete [] xerr;
  delete [] yerr;
  return gr;
}

TGraphErrors* makeGraphPr(TProfile* hist, TString name="")
{
  name.ReplaceAll(" ","");
  Int_t nbins = hist->GetNbinsX();
  Double_t* x = new Double_t[nbins];
  Double_t* y = new Double_t[nbins];
  Double_t* xerr = new Double_t[nbins]; 
  Double_t* yerr = new Double_t[nbins];
  Int_t n=0;
  for (Int_t i=0; i<nbins; i++)
  {
    //if (hist->GetBinContent(i+1)==0.0) continue;
    
    x[n] = hist->GetXaxis()->GetBinCenter(i+1);
    y[n] = hist->GetBinContent(i+1);
    xerr[n] = 0.0;
    yerr[n] = hist->GetBinError(i+1);
    n++;
  }
  TGraphErrors* gr = new TGraphErrors(n,x,y,xerr,yerr);
  delete [] x;
  delete [] y;
  delete [] xerr;
  delete [] yerr;
  return gr;
}

void myHistoSetUp
(
  TH1F *hist=0, TString xTitle="xTitle", TString yTitle="yTitle", float minY=-1, float maxY=1, int lineColor=1, int lineStyle=2,
 int nDivisionsX=305,int nDivisionsY=305
)
{
  hist->GetXaxis()->SetTitle(xTitle);
  hist->GetYaxis()->SetTitle(yTitle);
  hist->SetMinimum(minY);
  hist->SetMaximum(maxY);
  hist->SetLineColor(lineColor);
  hist->SetLineStyle(lineStyle);
  hist->SetNdivisions(nDivisionsX,"x");
  hist->SetNdivisions(nDivisionsY,"y");
  return;
}

// Shift original TGraphErrors along x-axis by amount determined by 'shift'.
void ShiftAlongXaxis(TGraphErrors *ge, Double_t shift)
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

void ScaleXaxis(TGraphErrors *ge, Double_t scale)
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
    x*=scale;
    ge->SetPoint(p,x,y);
  }
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
 int myFillStyle =fillStyle,
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
