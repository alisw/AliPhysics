#define NMAXCLASSES 100

void SetFrame(TH1F* frame){
  frame->GetXaxis()->SetTitleOffset(1.);
  frame->GetYaxis()->SetTitleOffset(1.3);
  frame->GetXaxis()->SetTitleSize(0.05);
  frame->GetYaxis()->SetTitleSize(0.05);
  frame->GetXaxis()->SetLabelSize(0.05);
  frame->GetYaxis()->SetLabelSize(0.05);
  frame->GetXaxis()->SetTitleFont(42);
  frame->GetYaxis()->SetTitleFont(42);
  frame->GetXaxis()->SetLabelFont(42);
  frame->GetYaxis()->SetLabelFont(42);
  frame->GetXaxis()->SetNdivisions(408);
}

TChain* t;
ULong64_t class_l2a[NMAXCLASSES] = {0};
Double_t  class_lumi[NMAXCLASSES] = {0};
TObjArray* classes;
TObjString* partition;
TObjString* lhcState;
TObjString* lhcPeriod;
TObjString* activeDetectors;
Int_t run;
Int_t fill;
Int_t timeStart;
Int_t timeEnd;



TGraph* GetStat(TString className="CINT7-B-NOPF-CENT",Bool_t lumi=1, Bool_t goodOnly=0, TString part="PHYSICS_1"){
  Double_t stat[100000];
  Double_t time_stamp[100000];
  TH1D* hStat      = new TH1D("hStat","",1,0,1);
  TH1D* hTimeStart = new TH1D("hTimeStart","",1,0,1);
  TH1D* hTimeEnd   = new TH1D("hTimeEnd","",1,0,1);
  TString classNameRun;

  for (Int_t r=0;r<t->GetEntries();r++){
    t->GetEntry(r);
    if (!partition->String().Contains(part.Data())) continue;
    if (!lhcState->String().Contains("STABLE")) continue;
    if (!lhcPeriod->String().Contains("LHC15m")) continue;
    hTimeStart->Fill(Form("%i",run),timeStart);
    hTimeEnd->Fill(Form("%i",run),timeEnd);
    classNameRun = TString(className);

    AliTriggerClass* cl = (AliTriggerClass*) classes->FindObject(classNameRun.Data());
    if (!cl) { 
      hStat->Fill(Form("%i",run),0.); 
      continue; 
    }
    
    if (goodOnly){
      if (classNameRun.Contains("CMUL")){
        if (!activeDetectors->String().Contains("MUONTRK") ||
            !activeDetectors->String().Contains("MUONTRG") ||
            !activeDetectors->String().Contains("T0") ||
            !activeDetectors->String().Contains("VZERO") ||
            !activeDetectors->String().Contains("ITSSPD")
        ){
          hStat->Fill(Form("%i",run),0.);
          continue;
        }
      }
    }
    Double_t x = lumi ? class_lumi[classes->IndexOf(cl)] : class_l2a[classes->IndexOf(cl)];
    hStat->Fill(Form("%i",run),x);
  }
  hStat->GetXaxis()->LabelsOption("a");
  hTimeEnd->GetXaxis()->LabelsOption("a");
  hTimeStart->GetXaxis()->LabelsOption("a");
  hStat->LabelsDeflate("x");
  hTimeStart->LabelsDeflate("x");
  hTimeEnd->LabelsDeflate("x");
  
  Int_t n = 2*hStat->GetNbinsX();

  for (Int_t i=0;i<hStat->GetNbinsX();i++){
    time_stamp[2*i  ] = hTimeStart->GetBinContent(i+1);
    time_stamp[2*i+1] = hTimeEnd->GetBinContent(i+1);
    stat[2*i]         = i>0?stat[2*i-1]:0;
    stat[2*i+1]       = stat[2*i]+hStat->GetBinContent(i+1)/(lumi ? 1000.: 1000000.);
  }

  delete hStat;
  delete hTimeStart;
  delete hTimeEnd;
  TGraph* gStat = new TGraph(n,time_stamp,stat);
  gStat->SetLineColor(kBlue);
  gStat->SetLineWidth(2);
  return gStat;
}

void integrated_lumi_ppref(Bool_t goodOnly=0){
  gStyle->SetPadTopMargin(0.01);
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadBottomMargin(0.06);
  gStyle->SetPadLeftMargin(0.13);
  TGaxis::SetMaxDigits(3);
  
  gStyle->SetOptTitle(1);
  gStyle->SetTitleOffset(1.6,"Y");
  gStyle->SetOptStat(0);
  TLatex* latex = new TLatex();
  latex->SetTextSize(0.05);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  
  t = new TChain("trending");
  t->AddFile("trending.root");
  classes = new TObjArray();
  partition = new TObjString();
  lhcState = new TObjString();
  lhcPeriod = new TObjString();
  activeDetectors = new TObjString();
  
  t->SetBranchAddress("run",&run);
  t->SetBranchAddress("fill",&fill);
  t->SetBranchAddress("classes",&classes);
  t->SetBranchAddress("class_l2a",&class_l2a);
  t->SetBranchAddress("class_lumi",&class_lumi);
  t->SetBranchAddress("timeStart",&timeStart);
  t->SetBranchAddress("timeEnd",&timeEnd);
  t->SetBranchAddress("partition",&partition);
  t->SetBranchAddress("lhcState",&lhcState);
  t->SetBranchAddress("lhcPeriod",&lhcPeriod);
  t->SetBranchAddress("activeDetectors",&activeDetectors);
  TGraph* gINT = GetStat("CINT7-B-NOPF-CENT",1,goodOnly);
  TGraph* gMUL = GetStat("CMUL7-B-NOPF-MUFAST",1,goodOnly);
  TGraph* gStatINT = GetStat("CINT7-B-NOPF-CENT",0,goodOnly);
  gMUL->SetLineColor(kBlack);
  gINT->SetLineColor(kBlue);
  gStatINT->SetLineColor(kBlue);

  TCanvas* c1 = new TCanvas("c1","",800,700);
  Double_t xminLumi   = gMUL->GetXaxis()->GetXmin();
  Double_t xmaxLumi   = gMUL->GetXaxis()->GetXmax();
  Double_t ymaxLumi   = gMUL->GetYaxis()->GetXmax()*1.1;
  TH1F* f1 = gPad->DrawFrame(xminLumi,0,xmaxLumi,ymaxLumi);
  SetFrame(f1);
  f1->GetYaxis()->SetTitle("Integrated luminosity, nb^{-1}");
  f1->GetXaxis()->SetTimeDisplay(1);
  f1->GetXaxis()->SetTimeFormat("%d %b");
  Double_t y = 0.94;
  latex->DrawLatex(0.18,0.94,"ALICE Performance, pp #sqrt{s} = 5.02 TeV");
  latex->SetTextColor(gINT->GetLineColor());
  latex->DrawLatex(0.18,y-=0.07,Form("MB triggers: L = %.3f nb^{-1}",gINT->GetY()[gINT->GetN()-1]));
  latex->SetTextColor(gMUL->GetLineColor());
  latex->DrawLatex(0.18,y-=0.07,Form("Dimuon triggers: L = %.3f nb^{-1}",gMUL->GetY()[gMUL->GetN()-1]));
  gMUL->Draw();
  gINT->Draw();
  gPad->Print("lumi_ppref.png");

  TCanvas* c2 = new TCanvas("c2","",800,700);
  Double_t xminEvents = gStatINT->GetXaxis()->GetXmin();
  Double_t xmaxEvents = gStatINT->GetXaxis()->GetXmax();
  Double_t ymaxEvents = gStatINT->GetYaxis()->GetXmax()*1.1;

  TH1F* f2 = gPad->DrawFrame(xminEvents,0,xmaxEvents,ymaxEvents);
  SetFrame(f2);
  f2->GetYaxis()->SetTitle("Recorded triggers, 10^{6}");
  f2->GetXaxis()->SetTimeDisplay(1);
  f2->GetXaxis()->SetTimeFormat("%d %b");
  gStatINT->Draw();
  latex->SetTextColor(1);
  latex->DrawLatex(0.18,0.94,"ALICE Performance, pp #sqrt{s} = 5.02 TeV");
  y = 0.94;
  latex->SetTextColor(gStatINT->GetLineColor());
  latex->DrawLatex(0.18,y-=0.07,Form("MB triggers: %.0fM",gStatINT->GetY()[gStatINT->GetN()-1]));
  gPad->Print("stat_ppref.png");
}
