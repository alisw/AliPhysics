#include "periodLevelQA.C"
#include "map"
using namespace std;

void trending_INT7_vs_INT8(TString inputFileName = "trending.root"){
  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(1.5);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(1);

  TFile* f = new TFile(inputFileName.Data());
  TTree* t = (TTree*) f->Get("trending");
  Int_t run                = 0;
  Int_t fill               = 0;
  Double_t mu              = 0;
  TObjString* lhcPeriod = new TObjString();
  ULong64_t alias_accepted[NBITS]={0};
  
  t->SetBranchAddress("run",&run);
  t->SetBranchAddress("alias_accepted",&alias_accepted);
  t->SetBranchAddress("fill",&fill);
  t->SetBranchAddress("lhcPeriod",&lhcPeriod);
  t->SetBranchAddress("mu",&mu);

  TH1D* hINT7  = new TH1D("hINT7","",1,0,1);
  TH1D* hINT8  = new TH1D("hINT8","",1,0,1);
  Int_t nRuns = t->GetEntries();
  map<Int_t,Int_t> fills;
  map<Int_t,TString> periods;

  for (Int_t r=0;r<nRuns;r++){
    t->GetEntry(r);
    char* srun = Form("%i",run);
    UInt_t int7 = alias_accepted[1];
    UInt_t int8 = alias_accepted[21];
    if (int8<=100000) continue;
//    if (mu>0.04) continue;
    fills[run]=fill;
    periods[run]=lhcPeriod->String();
    hINT7->Fill(srun,int7);
    hINT8->Fill(srun,int8);
  }
  hINT7->LabelsDeflate("X");
  hINT8->LabelsDeflate("X");
  hINT7->Sumw2();
  hINT8->Sumw2();
  TH1D* hRatio = (TH1D*) hINT7->Clone("hRatio");
  hRatio->Divide(hINT7,hINT8,1,1,"B");
  TCanvas* c = new TCanvas("cRatio","cRatio",1800,500);
  c->SetMargin(0.03,0.003,0.15,0.06);
  SetHisto(hRatio);
  hRatio->SetTitle("CINT7/CINT8 ratio after physics selection");
  hRatio->GetYaxis()->SetRangeUser(1.97,2.03);
  hRatio->Draw();
  AddFillSeparationLines(hRatio,fills);
  AddPeriodSeparationLines(hRatio,periods);
  gPad->Print("CINT7to8ratio.png");
}
