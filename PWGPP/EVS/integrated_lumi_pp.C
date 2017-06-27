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
  frame->GetXaxis()->SetNdivisions(306);
}

TChain* t;
ULong64_t class_l2a[NMAXCLASSES] = {0};
Double_t  class_lumi[NMAXCLASSES] = {0};
Double_t lumi_seen;
Double_t mu;
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
//    if (lhcPeriod->String().Contains("LHC16c")) continue;
//    if (run>255469) continue;
//    if (!lhcPeriod->String().Contains("LHC16f")) continue;
//    if (run>=265304) continue;
//    if (run<252238) continue;
    if (run<270565) continue;
    if (mu>0.1) continue;
    hTimeStart->Fill(Form("%i",run),timeStart);
    hTimeEnd->Fill(Form("%i",run),timeEnd);
    classNameRun = TString(className);
    if (className.Contains("CVHMV0M-B-NOPF-CENT")    && run>=254124 && run<=256503) classNameRun="CVHMV0M-B-NOPF-CENTNOTRD";
    if (className.Contains("CVHMV0M-B-NOPF-CENT")    && run>=256504               ) classNameRun="CVHMV0M-B-SPD2-CENT";
    if (className.Contains("CPHI7-B-NOPF-CENTNOTRD") && run>=257606 && run<=257684) classNameRun="CPHI7-B-NOPF-CENT";
    AliTriggerClass* cl = (AliTriggerClass*) classes->FindObject(classNameRun.Data());
    if (!cl && !className.Contains("Seen")) { 
      hStat->Fill(Form("%i",run),0.); 
      continue; 
    }

    if (goodOnly){
      if (classNameRun.Contains("CINT7")){
        if (!activeDetectors->String().Contains("TPC") ||
            !activeDetectors->String().Contains("T0") ||
            !activeDetectors->String().Contains("VZERO") ||
            !activeDetectors->String().Contains("ITSSPD")
        ){
          hStat->Fill(Form("%i",run),0.);
          continue;
        }
      }
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
    Double_t x;
    if (classNameRun.Contains("Seen")){
      x = lumi_seen;
    } else {
      x  = lumi ? class_lumi[classes->IndexOf(cl)] : class_l2a[classes->IndexOf(cl)];
    }

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
    stat[2*i+1]       = stat[2*i]+hStat->GetBinContent(i+1)/(lumi ? 1000000.: 1000000.);
  }

  delete hStat;
  delete hTimeStart;
  delete hTimeEnd;
  TGraph* gStat = new TGraph(n,time_stamp,stat);
  gStat->SetLineColor(kBlue);
  gStat->SetLineWidth(2);
  return gStat;
}

void integrated_lumi_pp(Bool_t goodOnly=0){
  gStyle->SetPadTopMargin(0.01);
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadBottomMargin(0.06);
  gStyle->SetPadLeftMargin(0.13);
  TGaxis::SetMaxDigits(3);
  
  gStyle->SetOptTitle(1);
  gStyle->SetTitleOffset(1.6,"Y");
  gStyle->SetOptStat(0);
  TLatex* latex = new TLatex();
  latex->SetTextSize(0.045);
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
  t->SetBranchAddress("lumi_seen",&lumi_seen);
  t->SetBranchAddress("mu",&mu);
  t->GetEntry(t->GetEntries()-1);
  TTimeStamp stamp = TTimeStamp(timeEnd,0,3600);
  printf("%s\n",stamp.AsString("s"));
//  TFile* fDelivered = new TFile("Fill_Summaries_2016/lumi.root");
//  TGraph* gDelivered = (TGraph*) fDelivered->Get("gIntegrated");
  TGraph* gSeen    = GetStat("Seen",1,0);

  TGraph* gINT = GetStat("CINT7-B-NOPF-CENT",1,goodOnly);
  TGraph* gV0M = GetStat("CVHMV0M-B-NOPF-CENT",1,goodOnly);
  TGraph* gMUL = GetStat("CMUL7-B-NOPF-MUFAST",1,goodOnly);
  TGraph* gMSL = GetStat("CMSL7-B-NOPF-MUFAST",1,goodOnly);
  TGraph* gEMC1 = GetStat("CEMC7EJ1-B-NOPF-CENTNOTRD",1,goodOnly);
  TGraph* gEMC2 = GetStat("CEMC7EJ2-B-NOPF-CENT",1,goodOnly);
  TGraph* gPHI  = GetStat("CPHI7-B-NOPF-CENTNOTRD",1,goodOnly);
  TGraph* gDG   = GetStat("CCUP25-B-SPD1-CENTNOTRD",1,goodOnly);
  TGraph* gStatV0M   = GetStat("CVHMV0M-B-NOPF-CENT",0,goodOnly);
  TGraph* gStatINT   = GetStat("CINT7-B-NOPF-CENT",0,goodOnly);
  TGraph* gStatINT10 = GetStat("CINT10-B-NOPF-CENTNOTRD",0,goodOnly);
  TGraph* gStatINT11 = GetStat("CINT11-B-NOPF-CENTNOTRD",0,goodOnly);

//  gDelivered->SetLineWidth(2);
//  gDelivered->SetLineColor(kBlack);
  gSeen->SetLineColor(kRed+2);
  gMUL->SetLineColor(kBlack);
  gMSL->SetLineColor(kGray+2);
  gV0M->SetLineColor(kGreen+2);
  gINT->SetLineColor(kBlue);
  gEMC1->SetLineColor(kMagenta);
  gEMC2->SetLineColor(kMagenta+2);
  gPHI->SetLineColor(kBlue+2);
  gDG->SetLineColor(kYellow+2);
  gStatV0M->SetLineColor(kGreen+2);
  gStatINT->SetLineColor(kBlue);
  gStatINT10->SetLineColor(kRed);
  gStatINT11->SetLineColor(kMagenta);

  TCanvas* c1 = new TCanvas("c1","",800,700);
  Double_t xminLumi = gSeen->GetXaxis()->GetXmin();
  Double_t xmaxLumi = gSeen->GetXaxis()->GetXmax();
  Double_t ymaxLumi = gSeen->GetYaxis()->GetXmax()*1.1;
  TH1F* f1 = gPad->DrawFrame(xminLumi,0,xmaxLumi,ymaxLumi);
  SetFrame(f1);
  f1->GetYaxis()->SetTitle("Integrated luminosity, pb^{-1}");
  f1->GetXaxis()->SetTimeDisplay(1);
  f1->GetXaxis()->SetTimeFormat("%d %b");
  Double_t y = 0.94;
  Double_t dy = 0.055;
  latex->DrawLatex(0.18,0.94,"ALICE Performance 2017, pp #sqrt{s} = 13 TeV");
  latex->DrawLatex(0.18,y-=dy,Form("%s",stamp.AsString("s")));
//  latex->SetTextColor(gDelivered->GetLineColor());
//  latex->DrawLatex(0.18,y-=0.06,Form("Delivered: %.1f ub^{-1}",gDelivered->GetY()[gDelivered->GetN()-1]));
  latex->SetTextColor(gSeen->GetLineColor());
  latex->DrawLatex(0.18,y-=dy,Form("Seen: %.3f pb^{-1}",gSeen->GetY()[gSeen->GetN()-1]));
  latex->SetTextColor(gMUL->GetLineColor());
  latex->DrawLatex(0.18,y-=dy,Form("Dimuon: %.3f pb^{-1}",gMUL->GetY()[gMUL->GetN()-1]));
  latex->SetTextColor(gMSL->GetLineColor());
  latex->DrawLatex(0.18,y-=dy,Form("MSL: %.3f pb^{-1}",gMSL->GetY()[gMSL->GetN()-1]));
  latex->SetTextColor(gV0M->GetLineColor());
  latex->DrawLatex(0.18,y-=dy,Form("V0 HM: %.3f pb^{-1}",gV0M->GetY()[gV0M->GetN()-1]));
  latex->SetTextColor(gEMC1->GetLineColor());
  latex->DrawLatex(0.18,y-=dy,Form("CALO high: %.3f pb^{-1}",gEMC1->GetY()[gEMC1->GetN()-1]));
  latex->SetTextColor(gEMC2->GetLineColor());
  latex->DrawLatex(0.18,y-=dy,Form("CALO low: %.3f pb^{-1}",gEMC2->GetY()[gEMC2->GetN()-1]));
  latex->SetTextColor(gPHI->GetLineColor());
  latex->DrawLatex(0.18,y-=dy,Form("PHOS: %.3f pb^{-1}",gPHI->GetY()[gPHI->GetN()-1]));
  latex->SetTextColor(gDG->GetLineColor());
  latex->DrawLatex(0.18,y-=dy,Form("DG: %.3f pb^{-1}",gDG->GetY()[gDG->GetN()-1]));
  latex->SetTextColor(gINT->GetLineColor());
  latex->DrawLatex(0.18,y-=dy,Form("INT7: %.3f pb^{-1}",gINT->GetY()[gINT->GetN()-1]));
  gSeen->Draw();
  gMSL->Draw();
  gMUL->Draw();
  gINT->Draw();
  gV0M->Draw();
  gEMC1->Draw();
  gEMC2->Draw();
  gPHI->Draw();
  gDG->Draw();

  gPad->Print("lumi.png");

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
  gStatV0M->Draw();
  gStatINT10->Draw();
  gStatINT11->Draw();
  latex->SetTextColor(1);
  latex->DrawLatex(0.18,0.94,"ALICE Performance 2017, pp #sqrt{s} = 13 TeV");
  y = 0.94;
//  y -= dy;
  latex->DrawLatex(0.18,y-=dy,Form("%s",stamp.AsString("s")));
  latex->SetTextColor(gStatINT->GetLineColor());
  latex->DrawLatex(0.18,y-=dy,Form("INT7: %.0fM",gStatINT->GetY()[gStatINT->GetN()-1]));
  latex->SetTextColor(gStatINT11->GetLineColor());
  latex->DrawLatex(0.18,y-=dy,Form("INT11: %.0fM",gStatINT11->GetY()[gStatINT11->GetN()-1]));
  latex->SetTextColor(gStatV0M->GetLineColor());
  latex->DrawLatex(0.18,y-=dy,Form("High Multiplicity: %.0fM",gStatV0M->GetY()[gStatV0M->GetN()-1]));
//  latex->SetTextColor(gStatINT10->GetLineColor());
//  latex->DrawLatex(0.18,y-=dy,Form("CINT10: %.0fM",gStatINT10->GetY()[gStatINT10->GetN()-1]));
  gPad->Print("stat.png");
}
