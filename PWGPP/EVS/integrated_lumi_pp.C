#ifndef __CINT__
#include "Includes.h"
#endif

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
Int_t timeStartSB = 1497139200;// June 12,2017, 00:00:00 (TTimeStamp)
Int_t timeDay = 84600;// 1 day in seconds
Int_t timeCorr = 8*timeDay;
Int_t timeSBday[366]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,29,30,31,32,33,34,35,36,37,38,39,40,41,42,47,48,49,50,// 37 days in this line
                  51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,
                  84,85,86,87,88,89,90,91,92,93,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,
                  123,124,125,126,127,128,129,130,131,132,133,139,140,141,142,143,144,145,146,147,148,149,150,151,152,
                  153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169};// Rank of Stable Beam Day since the first Stable Beam Day



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

TGraph* GetProjection(Double_t lumiProjectedPerDay,Double_t timeLimit,Int_t startDay = 0){
  Double_t stat[100000];
  Double_t time_stamp[100000];
  Int_t n = 0;

//  // Loop to prepare points for TGraph - linear version
//  for (Int_t iDay=0;iDay<366;iDay++){
//    if (timeStartSB+timeDay*timeSBday[iDay]>timeLimit) break;
//    n++;n++;
//    time_stamp[2*iDay  ] = timeStartSB+timeDay*timeSBday[iDay];
//    time_stamp[2*iDay+1] = timeStartSB+timeDay*timeSBday[iDay+1];
//    stat[2*iDay]         = lumiProjectedPerDay*iDay;
//    stat[2*iDay+1]       = lumiProjectedPerDay*(iDay+1);
////    printf("Iteration number %i, n: %i, time: %f, stat: %f\n",iDay,n,time_stamp[2*iDay],stat[2*iDay]);
//  }

  // Loop to prepare points for TGraph - step version
  time_stamp[0] = timeStartSB+timeDay*timeSBday[startDay];
  stat[0]       = 0;
  Int_t iDay(0);
  for (Int_t iSBday=1;iSBday<366-startDay;iSBday++){
    if (timeStartSB+timeDay*timeSBday[iSBday+startDay]>timeLimit) break;
    n++;
    if (timeSBday[iSBday+startDay]-timeSBday[iSBday-1+startDay]>1){
      time_stamp[iSBday] = timeStartSB+timeDay*timeSBday[iSBday+startDay];
      stat[iSBday]       = lumiProjectedPerDay*iDay;
      }
    else {
      iDay++;
      time_stamp[iSBday] = timeStartSB+timeDay*timeSBday[iSBday+startDay];
      stat[iSBday]       = lumiProjectedPerDay*iDay;
    }
//    printf("Iteration number %i, n: %i, time: %f, stat: %f\n",iDay,n,time_stamp[iDay],stat[iDay]);
  }

  TGraph* gProj = new TGraph(n,time_stamp,stat);
  gProj->SetLineColor(kBlue);
  gProj->SetLineWidth(2);

  return gProj;
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
  latex->SetTextSize(0.035);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  
  t = new TChain("trending");
  t->AddFile("trending_merged.root");
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

  TGraph* gINT          = GetStat("CINT7-B-NOPF-CENT",1,goodOnly);
  TGraph* gINTFAST      = GetStat("CINT7-B-NOPF-MUFAST",1,goodOnly);
  TGraph* gV0M          = GetStat("CVHMV0M-B-NOPF-CENT",1,goodOnly);
  TGraph* gMUL          = GetStat("CMUL7-B-NOPF-MUFAST",1,goodOnly);
  TGraph* gMSL          = GetStat("CMSL7-B-NOPF-MUFAST",1,goodOnly);
  TGraph* gEMC1         = GetStat("CEMC7EJ1-B-NOPF-CENTNOTRD",1,goodOnly);
  TGraph* gEMC2         = GetStat("CEMC7EJ2-B-NOPF-CENT",1,goodOnly);
  TGraph* gPHI          = GetStat("CPHI7-B-NOPF-CENTNOTRD",1,goodOnly);
  TGraph* gDG           = GetStat("CCUP25-B-SPD1-CENTNOTRD",1,goodOnly);
  TGraph* gStatV0M      = GetStat("CVHMV0M-B-NOPF-CENT",0,goodOnly);
  TGraph* gStatINT      = GetStat("CINT7-B-NOPF-CENT",0,goodOnly);
  TGraph* gStatINT10    = GetStat("CINT10-B-NOPF-CENTNOTRD",0,goodOnly);
  TGraph* gStatINT11    = GetStat("CINT11-B-NOPF-CENTNOTRD",0,goodOnly);
  TGraph* gMuonCaloHigh = GetStat("CEMC7MUL-B-NOPF-ALLNOTRD",1,goodOnly);
  TGraph* gMuonCaloLow  = GetStat("CEMC7MSL-B-NOPF-ALLNOTRD",1,goodOnly);
  TGraph* gTRDqrkNuc    = GetStat("CINT7HNU-T-NOPF-CENT",1,goodOnly);
  TGraph* gTRDjet       = GetStat("CINT7HJT-T-NOPF-CENT",1,goodOnly);
  // Projection graphs
  TGraph* gProjMUL      = GetProjection(0.115,gSeen->GetXaxis()->GetXmax()-timeCorr);
  TGraph* gProjV0M      = GetProjection(0.053,gSeen->GetXaxis()->GetXmax()-timeCorr);
  TGraph* gProjStatINT  = GetProjection(7.390,gSeen->GetXaxis()->GetXmax()-timeCorr);
  TGraph* gProjTRDjet   = GetProjection(0.0095,gSeen->GetXaxis()->GetXmax()-timeCorr,53);// Start 17/8/2017, which is 53. day of stable beams


//  gDelivered->SetLineWidth(2);
//  gDelivered->SetLineColor(kBlack);
  gSeen->SetLineColor(kRed+2);
  gMUL->SetLineColor(kBlack);
  gMSL->SetLineColor(kGray+2);
  gV0M->SetLineColor(kGreen+2);
  gINT->SetLineColor(kBlue);
  gINTFAST->SetLineColor(kBlue);
  gEMC1->SetLineColor(kMagenta);
  gEMC2->SetLineColor(kMagenta+2);
  gPHI->SetLineColor(kBlue+2);
  gDG->SetLineColor(kYellow+2);
  gStatV0M->SetLineColor(kGreen+2);
  gStatINT->SetLineColor(kBlue);
  gStatINT10->SetLineColor(kRed);
  gStatINT11->SetLineColor(kMagenta);
  gMuonCaloHigh->SetLineColor(kRed);
  gMuonCaloLow->SetLineColor(kYellow+4);
  gTRDqrkNuc->SetLineColor(kCyan+2);
  gTRDjet->SetLineColor(kCyan+4);
  gProjMUL->SetLineColor(kGray+2);
  gProjMUL->SetLineStyle(9);
  gProjV0M->SetLineColor(kRed+1);
  gProjV0M->SetLineStyle(9);
  gProjStatINT->SetLineColor(kRed+1);
  gProjStatINT->SetLineStyle(9);
  gProjTRDjet->SetLineColor(kCyan-5);
  gProjTRDjet->SetLineStyle(9);

  TCanvas* c1 = new TCanvas("c1","",800,700);
  Double_t xminLumi = gSeen->GetXaxis()->GetXmin();
  Double_t xmaxLumi = gSeen->GetXaxis()->GetXmax();
  Double_t ymaxLumi = gSeen->GetYaxis()->GetXmax()*1.1;
  TH1F* f1 = gPad->DrawFrame(xminLumi,0,xmaxLumi,ymaxLumi);
  SetFrame(f1);
  f1->GetYaxis()->SetTitle("Integrated luminosity, pb^{-1}");
  f1->GetXaxis()->SetTimeDisplay(1);
  f1->GetXaxis()->SetTimeFormat("%d %b");
  Double_t x = 0.17;
  Double_t y = 0.94;
  Double_t dy = 0.045;
  latex->DrawLatex(x,0.94,"ALICE Performance 2017, pp #sqrt{s} = 13 TeV");
  latex->DrawLatex(x,y-=dy,Form("%s",stamp.AsString("s")));
//  latex->SetTextColor(gDelivered->GetLineColor());
//  latex->DrawLatex(x,y-=0.06,Form("Delivered: %.1f ub^{-1}",gDelivered->GetY()[gDelivered->GetN()-1]));
  latex->SetTextColor(gSeen->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("Seen: %.3f pb^{-1}",gSeen->GetY()[gSeen->GetN()-1]));
  latex->SetTextColor(gMUL->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("Di-#mu, single #mu high-p_{T}: %.3f pb^{-1} / 15 pb^{-1}",gMUL->GetY()[gMUL->GetN()-1]));
  latex->SetTextColor(gMSL->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("Single #mu low-p_{T}: %.3f pb^{-1} / 0.9 pb^{-1}",gMSL->GetY()[gMSL->GetN()-1]));
  latex->SetTextColor(gTRDqrkNuc->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("TRD quarkonia and nuclei: %.3f pb^{-1} / 0.85 pb^{-1}",gTRDqrkNuc->GetY()[gTRDqrkNuc->GetN()-1]));
  latex->SetTextColor(gTRDjet->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("TRD jet: %.3f pb^{-1} / 0.85 pb^{-1}",gTRDjet->GetY()[gTRDjet->GetN()-1]));
  latex->SetTextColor(gMuonCaloHigh->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("#mu-CALO high: %.3f pb^{-1}",gMuonCaloHigh->GetY()[gMuonCaloHigh->GetN()-1]));
  latex->SetTextColor(gMuonCaloLow->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("#mu-CALO low: %.3f pb^{-1}",gMuonCaloLow->GetY()[gMuonCaloLow->GetN()-1]));
  latex->SetTextColor(gEMC1->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("E/DCal high: %.3f pb^{-1} / 7 pb^{-1}",gEMC1->GetY()[gEMC1->GetN()-1]));
  latex->SetTextColor(gEMC2->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("E/DCal low: %.3f pb^{-1} / 0.8 pb^{-1}",gEMC2->GetY()[gEMC2->GetN()-1]));
  latex->SetTextColor(gPHI->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("PHOS: %.3f pb^{-1} / 7 pb^{-1}",gPHI->GetY()[gPHI->GetN()-1]));
  latex->SetTextColor(gDG->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("DG: %.3f pb^{-1} / 7 pb^{-1}",gDG->GetY()[gDG->GetN()-1]));
  latex->SetTextColor(gV0M->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("V0 HM: %.3f pb^{-1} / 7 pb^{-1}",gV0M->GetY()[gV0M->GetN()-1]));
  latex->SetTextColor(gINT->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("INT7: %.4f pb^{-1}",gINT->GetY()[gINT->GetN()-1]));

  gSeen->Draw();
  gMSL->Draw();
  gMUL->Draw();
  gINT->Draw();
  gV0M->Draw();
  gEMC1->Draw();
  gEMC2->Draw();
  gPHI->Draw();
  gDG->Draw();
  gMuonCaloHigh->Draw();
  gMuonCaloLow->Draw();
  gTRDqrkNuc->Draw();
  gTRDjet->Draw();

  gPad->Print("lumi.png");
  delete f1;

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
  latex->DrawLatex(x,0.94,"ALICE Performance 2017, pp #sqrt{s} = 13 TeV");
  y = 0.94;
//  y -= dy;
  latex->DrawLatex(x,y-=dy,Form("%s",stamp.AsString("s")));
  latex->SetTextColor(gStatINT->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("INT7: %.0fM / 845M full-B, 145M / 150M low-B",gStatINT->GetY()[gStatINT->GetN()-1]-145));
  latex->SetTextColor(gStatINT11->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("INT11: %.0fM / 100M",gStatINT11->GetY()[gStatINT11->GetN()-1]));
  latex->SetTextColor(gStatV0M->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("High Multiplicity: %.0fM / 750M",gStatV0M->GetY()[gStatV0M->GetN()-1]));
//  latex->SetTextColor(gStatINT10->GetLineColor());
//  latex->DrawLatex(x,y-=dy,Form("CINT10: %.0fM",gStatINT10->GetY()[gStatINT10->GetN()-1]));
  gPad->Print("stat.png");
  delete f2;

  TCanvas* cCentralBarrel = new TCanvas("cCentralBarrel","",800,700);
  xminLumi = gPHI->GetXaxis()->GetXmin();
  xmaxLumi = gPHI->GetXaxis()->GetXmax();
  ymaxLumi = gPHI->GetYaxis()->GetXmax()*1.1;
  TH1F* f1 = gPad->DrawFrame(xminLumi,0,xmaxLumi,ymaxLumi);
  SetFrame(f1);
  f1->GetYaxis()->SetTitle("Integrated luminosity, pb^{-1}");
  f1->GetXaxis()->SetTimeDisplay(1);
  f1->GetXaxis()->SetTimeFormat("%d %b");
  y = 0.94;
  latex->SetTextColor(1);
  latex->DrawLatex(x,0.94,"ALICE Performance 2017, pp #sqrt{s} = 13 TeV");
  latex->DrawLatex(x,y-=dy,"Central Barrel");
  latex->DrawLatex(x,y-=dy,Form("%s",stamp.AsString("s")));
  latex->SetTextColor(gMuonCaloHigh->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("#mu-CALO high: %.3f pb^{-1}",gMuonCaloHigh->GetY()[gMuonCaloHigh->GetN()-1]));
  latex->SetTextColor(gMuonCaloLow->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("#mu-CALO low: %.3f pb^{-1}",gMuonCaloLow->GetY()[gMuonCaloLow->GetN()-1]));
  latex->SetTextColor(gEMC1->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("E/DCal high: %.3f pb^{-1} / 7 pb^{-1}",gEMC1->GetY()[gEMC1->GetN()-1]));
  latex->SetTextColor(gEMC2->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("E/DCal low: %.3f pb^{-1} / 0.8 pb^{-1}",gEMC2->GetY()[gEMC2->GetN()-1]));
  latex->SetTextColor(gPHI->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("PHOS: %.3f pb^{-1} / 7 pb^{-1}",gPHI->GetY()[gPHI->GetN()-1]));
  latex->SetTextColor(gDG->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("DG: %.3f pb^{-1} / 7 pb^{-1}",gDG->GetY()[gDG->GetN()-1]));
  latex->SetTextColor(gV0M->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("V0-HM: %.3f pb^{-1} / 7 pb^{-1}",gV0M->GetY()[gV0M->GetN()-1]));
  latex->SetTextColor(gINT->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("INT7: %.4f pb^{-1}",gINT->GetY()[gINT->GetN()-1]));

  gMuonCaloHigh->Draw();
  gMuonCaloLow->Draw();;
  gEMC1->Draw();
  gEMC2->Draw();
  gPHI->Draw();
  gDG->Draw();
  gV0M->Draw();
  gINT->Draw();

  gPad->Print("lumi_central-barrel.png");
  delete f1;

  TCanvas* cMuons = new TCanvas("cMuons","",800,700);
  xminLumi = gMUL->GetXaxis()->GetXmin();
  xmaxLumi = gMUL->GetXaxis()->GetXmax();
  ymaxLumi = gMUL->GetYaxis()->GetXmax()*1.1;
  TH1F* f1 = gPad->DrawFrame(xminLumi,0,xmaxLumi,ymaxLumi);
  SetFrame(f1);
  f1->GetYaxis()->SetTitle("Integrated luminosity, pb^{-1}");
  f1->GetXaxis()->SetTimeDisplay(1);
  f1->GetXaxis()->SetTimeFormat("%d %b");
  y = 0.94;
  latex->SetTextColor(1);
  latex->DrawLatex(x,0.94,"ALICE Performance 2017, pp #sqrt{s} = 13 TeV");
  latex->DrawLatex(x,y-=dy,"Muon Arm");
  latex->DrawLatex(x,y-=dy,Form("%s",stamp.AsString("s")));
  latex->SetTextColor(gMUL->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("Di-#mu, single #mu high-p_{T}: %.3f pb^{-1} / 15 pb^{-1}",gMUL->GetY()[gMUL->GetN()-1]));
  latex->SetTextColor(gMSL->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("Single #mu low-p_{T}: %.3f pb^{-1} / 0.9 pb^{-1}",gMSL->GetY()[gMSL->GetN()-1]));
  latex->SetTextColor(gMuonCaloHigh->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("#mu-CALO high: %.3f pb^{-1}",gMuonCaloHigh->GetY()[gMuonCaloHigh->GetN()-1]));
  latex->SetTextColor(gMuonCaloLow->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("#mu-CALO low: %.3f pb^{-1}",gMuonCaloLow->GetY()[gMuonCaloLow->GetN()-1]));
  latex->SetTextColor(gINTFAST->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("INT7-MUFAST: %.4f pb^{-1}",gINTFAST->GetY()[gINTFAST->GetN()-1]));

  gMSL->Draw();
  gMUL->Draw();
  gMuonCaloHigh->Draw();
  gMuonCaloLow->Draw();;
  gINT->Draw();

  gPad->Print("lumi_muon-arm.png");
  delete f1;

  TCanvas* cProj = new TCanvas("cProj","",800,700);
  xminLumi = gProjMUL->GetXaxis()->GetXmin();
  xmaxLumi = gProjMUL->GetXaxis()->GetXmax();
  ymaxLumi = gProjMUL->GetYaxis()->GetXmax()*1.1;
  TH1F* f1 = gPad->DrawFrame(xminLumi,0,xmaxLumi,ymaxLumi);
  SetFrame(f1);
  f1->GetYaxis()->SetTitle("Integrated luminosity, pb^{-1}");
  f1->GetXaxis()->SetTimeDisplay(1);
  f1->GetXaxis()->SetTimeFormat("%d %b");
  y = 0.94;
  latex->SetTextColor(1);
  latex->DrawLatex(x,0.94,"ALICE Performance 2017, pp #sqrt{s} = 13 TeV, Muon Triggers");
  latex->DrawLatex(x,y-=dy,Form("%s",stamp.AsString("s")));
  latex->SetTextColor(gProjMUL->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("Lumi expected (Di-#mu): %.3f pb^{-1} / 15 pb^{-1}",gProjMUL->GetY()[gProjMUL->GetN()-1]));
  latex->SetTextColor(gProjV0M->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("Lumi expected (E/DCal, PHOS, V0-HM): %.3f pb^{-1} / 7 pb^{-1}",gProjV0M->GetY()[gProjV0M->GetN()-1]));
  latex->SetTextColor(gProjTRDjet->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("Lumi expected (TRD): %.3f pb^{-1} / 0.85 pb^{-1}",gProjTRDjet->GetY()[gProjTRDjet->GetN()-1]));
  latex->SetTextColor(gSeen->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("Seen: %.3f pb^{-1}",gSeen->GetY()[gSeen->GetN()-1]));
  latex->SetTextColor(gMUL->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("Di-#mu, single #mu high-p_{T}: %.3f pb^{-1} / 15 pb^{-1}",gMUL->GetY()[gMUL->GetN()-1]));
  latex->SetTextColor(gTRDqrkNuc->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("TRD quarkonia and nuclei: %.3f pb^{-1} / 0.85 pb^{-1}",gTRDqrkNuc->GetY()[gTRDqrkNuc->GetN()-1]));
  latex->SetTextColor(gTRDjet->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("TRD jet: %.3f pb^{-1} / 0.85 pb^{-1}",gTRDjet->GetY()[gTRDjet->GetN()-1]));
  latex->SetTextColor(gEMC1->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("E/DCal high: %.3f pb^{-1} / 7 pb^{-1}",gEMC1->GetY()[gEMC1->GetN()-1]));
  latex->SetTextColor(gPHI->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("PHOS: %.3f pb^{-1} / 7 pb^{-1}",gPHI->GetY()[gPHI->GetN()-1]));
  latex->SetTextColor(gV0M->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("V0-HM: %.3f pb^{-1} / 7 pb^{-1}",gV0M->GetY()[gV0M->GetN()-1]));

  gSeen->Draw();
  gProjMUL->Draw();
  gMUL->Draw();
  gEMC1->Draw();
  gPHI->Draw();
  gV0M->Draw();
  gTRDqrkNuc->Draw();
  gTRDjet->Draw();
  gProjV0M->Draw();
  gProjTRDjet->Draw();

  gPad->Print("lumi_projection.png");
  delete f1;

  TCanvas* cStatProj = new TCanvas("cStatProj","",800,700);
  xminEvents = gStatINT->GetXaxis()->GetXmin();
  xmaxEvents = gStatINT->GetXaxis()->GetXmax();
  ymaxEvents = gStatINT->GetYaxis()->GetXmax()*1.1;

  TH1F* f2 = gPad->DrawFrame(xminEvents,0,xmaxEvents,ymaxEvents);
  SetFrame(f2);
  f2->GetYaxis()->SetTitle("Recorded triggers, 10^{6}");
  f2->GetXaxis()->SetTimeDisplay(1);
  f2->GetXaxis()->SetTimeFormat("%d %b");
  gStatINT->Draw();
  gProjStatINT->Draw();
  gStatV0M->Draw();
  gStatINT10->Draw();
  gStatINT11->Draw();
  latex->SetTextColor(1);
  latex->DrawLatex(x,0.94,"ALICE Performance 2017, pp #sqrt{s} = 13 TeV");
  y = 0.94;
//  y -= dy;
  latex->DrawLatex(x,y-=dy,Form("%s",stamp.AsString("s")));
  latex->SetTextColor(gProjStatINT->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("Projected INT7: %.0fM / 845M full-B, 145M / 150M low-B",gProjStatINT->GetY()[gProjStatINT->GetN()-1]-145));
  latex->SetTextColor(gStatINT->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("INT7: %.0fM / 845M full-B, 145M / 150M low-B",gStatINT->GetY()[gStatINT->GetN()-1]-145));
  latex->SetTextColor(gStatINT11->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("INT11: %.0fM / 100M",gStatINT11->GetY()[gStatINT11->GetN()-1]));
  latex->SetTextColor(gStatV0M->GetLineColor());
  latex->DrawLatex(x,y-=dy,Form("High Multiplicity: %.0fM / 750M",gStatV0M->GetY()[gStatV0M->GetN()-1]));
//  latex->SetTextColor(gStatINT10->GetLineColor());
//  latex->DrawLatex(x,y-=dy,Form("CINT10: %.0fM",gStatINT10->GetY()[gStatINT10->GetN()-1]));
  gPad->Print("stat_projection.png");
  delete f2;
//  TTimeStamp* ts = new TTimeStamp();
//  ts->SetSec(gSeen->GetXaxis()->GetXmax()-timeCorr);
//  Printf("%f = %s \n",gSeen->GetXaxis()->GetXmax()-timeCorr,ts->AsString("s"));

}
