#define NBITS 29
#define NMAXCLASSES 100
//#include "TObjArray.h"
//#include "AliTriggerClass.h"

//using namespace std;


//vector<Int_t> goodRuns;

Int_t goodRuns[100000];
Int_t nGoodRuns=0;
Int_t nGoodRuns15f=21;
Int_t nGoodRuns15h=48;
Int_t nGoodRuns15i_badIR = 18;
Int_t nGoodRuns15i_badSPD = 10;
Int_t nGoodRuns15i=107;
Int_t nGoodRuns15j=145;
Int_t nGoodRuns15l=90;

Int_t goodRuns15f[] = {
225305, 225307, 225309, 225310, 225313, 225314, 225315, 225322, 225576, 225578,
225579, 225580, 225582, 225586, 225587, 225705, 225707, 225708, 225709, 225710,
225716
};

Int_t goodRuns15h[] = {
232914, 232916, 232977, 232986, 232993, 232995, 233020, 233059, 233061, 233093,
233116, 233120, 233144, 233169, 233172, 233217, 233232, 233233, 233239, 233242,
233361, 233465, 233472, 233614, 233621, 233623, 233627, 233743, 233799, 233828,
233830, 233837, 233858, 233912, 233969, 233971, 233972, 233973, 233974, 233975,
233976, 233977, 233978, 234031, 234039, 234048, 234049, 234050
};

Int_t goodRuns15i_badIR[] = {
    235345, 235346, 235347, 235356, 235362, 235364, 235380, 235383, 235423, 235432,
    235435, 235436, 235443, 235450, 235454, 235457, 235459, 235462
};

Int_t goodRuns15i_badSPD[] = {
    235196, 235201, 235203, 235204, 235226, 235242, 235245, 235344, 235547, 235573
};

Int_t goodRuns15i[] = {
235196, 235201, 235203, 235204, 235226, 235242, 235245, 235344, 235345, 235346,
235347, 235356, 235362, 235364, 235380, 235383, 235423, 235432, 235435, 235436,
235443, 235450, 235454, 235457, 235459, 235462, 235547, 235573, 235683, 235684,
235685, 235687, 235694, 235710, 235714, 235721, 235759, 235811, 235813, 235839,
235841, 235892, 235893, 235895, 235896, 235897, 235898, 236062, 236137, 236138,
236150, 236151, 236153, 236158, 236159, 236161, 236163, 236164, 236203, 236204,
236222, 236227, 236234, 236238, 236240, 236242, 236244, 236246, 236248, 236281,
236284, 236285, 236331, 236334, 236337, 236348, 236349, 236352, 236353, 236354,
236356, 236357, 236359, 236360, 236386, 236389, 236393, 236395, 236397, 236441,
236443, 236444, 236446, 236453, 236459, 236462, 236541, 236554, 236556, 236558,
236562, 236563, 236564, 236565, 236569, 236575, 236866
};

Int_t goodRuns15j[] = {
236892, 236893, 236965, 236968, 236969, 236972, 236973, 237001, 237029, 237049,
237051, 237056, 237061, 237104, 237111, 237112, 237115, 237174, 237180, 237244,
237245, 237251, 237253, 237255, 237259, 237286, 237288, 237289, 237291, 237330,
237335, 237342, 237350, 237356, 237364, 237368, 237372, 237391, 237394, 237396,
237397, 237400, 237401, 237406, 237408, 237409, 237507, 237512, 237515, 237645,
237670, 237671, 237675, 237676, 237678, 237681, 237684, 237691, 237698, 237699,
237705, 237706, 237707, 237708, 237710, 237711, 237713, 237765, 237768, 237777,
237779, 237780, 237782, 237787, 237789, 237790, 237791, 237793, 237795, 237796,
237806, 237842, 237844, 237845, 237847, 237945, 237948, 237969, 237978, 237982,
237983, 238073, 238091, 238097, 238129, 238131, 238132, 238133, 238136, 238139,
238140, 238142, 238144, 238145, 238147, 238148, 238154, 238159, 238160, 238164,
238170, 238176, 238179, 238184, 238185, 238187, 238395, 238451, 238454, 238455,
238456, 238457, 238458, 238459, 238460, 238466, 238468, 238469, 238470, 238472,
238474, 238570, 238576, 238578, 238579, 238583, 238587, 238594, 238598, 238604,
238606, 238607, 238610, 238614, 238621 
};

Int_t goodRuns15l[] = {
      239319, 239324, 239518, 239519, 240069, 240165, 240183, 240194, 240196, 240201, 
      240204, 240212, 240220, 240241, 240250, 240256, 240262, 240263, 240265, 240271, 
      240274, 240293, 240303, 240312, 240376, 240380, 240381, 240382, 240385, 240392, 
      240394, 240404, 240411, 240427, 240443, 240444, 240447, 240450, 240452, 240610, 
      240612, 240845, 240854, 240860, 240864, 240872, 240874, 240875, 240880, 241001, 
      241010, 241014, 241021, 241032, 241043, 241047, 241050, 241054, 241055, 241056, 
      241057, 241062, 241069, 241075, 241141, 241144, 241152, 241176, 241187, 241196, 
      241222, 241245, 241257, 241261, 241263, 241267, 241268, 241269, 241281, 241288, 
      241295, 241296, 241354, 241360, 241361, 241393, 241396, 241407, 241412, 241526
};

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



TGraph* GetStat(TString className="CVHMV0M-B-NOPF-CENTNOTRD",Bool_t lumi=1, Bool_t goodOnly=0, TString part="PHYSICS_1"){
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
    if (run<225000) continue;
    hTimeStart->Fill(Form("%i",run),timeStart);
    hTimeEnd->Fill(Form("%i",run),timeEnd);
    classNameRun = TString(className);
    if (className.EqualTo("CEMC7-B-NOPF-CENTNOTRD")                  && run<=236224) classNameRun="CEMC7-ABCE-NOPF-ALLNOTRD";
    if (className.EqualTo("CDMC7-B-NOPF-CENTNOTRD")                  && run<=236224) classNameRun="CDMC7-ABCE-NOPF-ALLNOTRD";
    if (className.EqualTo("CMUL7-B-NOPF-MUFAST")                     && run<=226606) classNameRun="CMUL7-B-NOPF-ALLNOTRD";
    if (className.EqualTo("CVHMV0M-B-NOPF-CENTNOTRD")                && run<=236221) classNameRun="CVHMV0M-B-NOPF-CENT";
    if (className.EqualTo("CVHMSH2-B-NOPF-CENTNOTRD")                && run<=236226) classNameRun="CVHMSH2-B-NOPF-CENT";
    if (className.EqualTo("CVHMV0M-B-NOPF-CENTNOTRD") && run>=238432 && run<=239144) classNameRun="CVHMV0M-B-SPD1-CENTNOTRD"; 
    if (className.EqualTo("CVHMSH2-B-NOPF-CENTNOTRD") && run>=238432               ) classNameRun="CVHMSH2-B-SPD1-CENTNOTRD";
    if (className.EqualTo("CINT7-B-NOPF-CENT")        && run>=225000 && run<=228935) classNameRun="CINT7-B-NOPF-ALLNOTRD";
    if (className.EqualTo("CINT7-B-NOPF-CENT")        && run>=228936 && run<=229893) classNameRun="CINT7-B-NOPF-CENTNOTRD";
    if (className.EqualTo("CINT7-B-NOPF-CENT")        && run>=229894 && run<=229899) classNameRun="CINT7-B-NOPF-ALLNOTRD";
    if (className.EqualTo("CINT7-B-NOPF-CENT")        && run>=229900 && run<=233911) classNameRun="CINT7-B-NOPF-CENT";
    if (className.EqualTo("CINT7-B-NOPF-CENT")        && run>=233912 && run<=234050) classNameRun="CINT7-B-NOPF-ALLNOTRD";
    if (className.EqualTo("CINT7-B-NOPF-CENT")        && run>=238890 && run<=239144) classNameRun="CINT7-I-NOPF-CENTNOTRD";
    if (className.EqualTo("CEMC7-B-NOPF-CENTNOTRD")                  && run<=236224) classNameRun="CEMC7-ABCE-NOPF-ALLNOTRD";
    if (className.EqualTo("CDMC7-B-NOPF-CENTNOTRD")                  && run<=236224) classNameRun="CDMC7-ABCE-NOPF-ALLNOTRD";

    AliTriggerClass* cl = (AliTriggerClass*) classes->FindObject(classNameRun.Data());
    if (!cl) { 
      hStat->Fill(Form("%i",run),0.); 
      continue; 
    }
    
    if (goodOnly){
      if (classNameRun.Contains("CINT7") || classNameRun.Contains("CVHM")){
        Int_t good=0;
        for (Int_t i=0;i<nGoodRuns;i++) good|= (goodRuns[i]==run);
        if (!good) {
          hStat->Fill(Form("%i",run),0.);
          continue;
        }
      }

      if (className.Contains("MC7")){
        if (run<235709 || run==236855 || run==236858 || run==236861 || fill==4440) {
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
//    printf("%i\n",classes->IndexOf(cl));
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
    stat[2*i+1]       = stat[2*i]+hStat->GetBinContent(i+1)/1000000.;
  }

  delete hStat;
  delete hTimeStart;
  delete hTimeEnd;
  TGraph* gStat = new TGraph(n,time_stamp,stat);
  gStat->SetLineColor(kBlue);
  gStat->SetLineWidth(2);
  return gStat;
}

void integrated_lumi(Bool_t goodOnly=0){
  nGoodRuns+= nGoodRuns15f;
  nGoodRuns+= nGoodRuns15h;
  nGoodRuns+= nGoodRuns15i;
  nGoodRuns+= nGoodRuns15j;
  nGoodRuns+= nGoodRuns15l;
  Int_t j=0;
  for (Int_t i=0;i<nGoodRuns15f;i++) goodRuns[j++]=goodRuns15f[i];
  for (Int_t i=0;i<nGoodRuns15h;i++) goodRuns[j++]=goodRuns15h[i];
//  for (Int_t i=0;i<nGoodRuns15i_badIR;i++) goodRuns[j++]=goodRuns15i_badIR[i];
//  for (Int_t i=0;i<nGoodRuns15i_badSPD;i++) goodRuns[j++]=goodRuns15i_badSPD[i];
  for (Int_t i=0;i<nGoodRuns15i;i++) goodRuns[j++]=goodRuns15i[i];
  for (Int_t i=0;i<nGoodRuns15j;i++) goodRuns[j++]=goodRuns15j[i];
  for (Int_t i=0;i<nGoodRuns15l;i++) goodRuns[j++]=goodRuns15l[i];


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
  TGraph* gV0M = GetStat("CVHMV0M-B-NOPF-CENTNOTRD",1,goodOnly);
  TGraph* gSH2 = GetStat("CVHMSH2-B-NOPF-CENTNOTRD",1,goodOnly);
  TGraph* gStatINT = GetStat("CINT7-B-NOPF-CENT",0,goodOnly);
  TGraph* gStatV0M = GetStat("CVHMV0M-B-NOPF-CENTNOTRD",0,goodOnly);
  TGraph* gStatSH2 = GetStat("CVHMSH2-B-NOPF-CENTNOTRD",0,goodOnly);
  TGraph* gStatEMC = GetStat("CEMC7-B-NOPF-CENTNOTRD",0,goodOnly,"PHYSICS_2");
  TGraph* gStatDMC = GetStat("CDMC7-B-NOPF-CENTNOTRD",0,goodOnly,"PHYSICS_2");
  gMUL->SetLineColor(kBlack);
  gINT->SetLineColor(kBlue);
  gV0M->SetLineColor(kMagenta);
  gSH2->SetLineColor(kRed);
  gStatV0M->SetLineColor(kMagenta);
  gStatSH2->SetLineColor(kRed);
  gStatEMC->SetLineColor(kGray);
  gStatDMC->SetLineColor(kGreen-2);

  TCanvas* c1 = new TCanvas("c1","",800,700);
  Double_t xminLumi   = gMUL->GetXaxis()->GetXmin();
  Double_t xmaxLumi   = gMUL->GetXaxis()->GetXmax();
  Double_t ymaxLumi   = gMUL->GetYaxis()->GetXmax()*1.1;
  TH1F* f1 = gPad->DrawFrame(xminLumi,0,xmaxLumi,ymaxLumi);
  SetFrame(f1);
  f1->GetYaxis()->SetTitle("Integrated luminosity, pb^{-1}");
  f1->GetXaxis()->SetTimeDisplay(1);
  f1->GetXaxis()->SetTimeFormat("%d %b");
  Double_t y = 0.94;
  latex->DrawLatex(0.18,0.94,"ALICE Performance, pp #sqrt{s} = 13 TeV");
  latex->SetTextColor(gINT->GetLineColor());
  latex->DrawLatex(0.18,y-=0.07,Form("MB triggers: L = %.3f pb^{-1}",gINT->GetY()[gINT->GetN()-1]));
  latex->SetTextColor(gV0M->GetLineColor());
  latex->DrawLatex(0.18,y-=0.07,Form("V0 HM triggers:  L = %.3f pb^{-1}",gV0M->GetY()[gV0M->GetN()-1]));
  latex->SetTextColor(gSH2->GetLineColor());
  latex->DrawLatex(0.18,y-=0.07,Form("SPD HM triggers: L = %.3f pb^{-1}",gSH2->GetY()[gSH2->GetN()-1]));
  latex->SetTextColor(gMUL->GetLineColor());
  latex->DrawLatex(0.18,y-=0.07,Form("Dimuon triggers: L = %.3f pb^{-1}",gMUL->GetY()[gMUL->GetN()-1]));
  gMUL->Draw();
  gINT->Draw();
  gV0M->Draw();
  gSH2->Draw();
  gPad->Print("lumi_dimuon_triggers.png");

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
  gStatSH2->Draw();
//  gStatDMC->Draw();
//  gStatEMC->Draw();
  latex->SetTextColor(1);
  latex->DrawLatex(0.18,0.94,"ALICE Performance, pp #sqrt{s} = 13 TeV");
  y = 0.94;
  latex->SetTextColor(gStatINT->GetLineColor());
  latex->DrawLatex(0.18,y-=0.07,Form("MB triggers: %.0fM",gStatINT->GetY()[gStatINT->GetN()-1]));
  latex->SetTextColor(gStatV0M->GetLineColor());
  latex->DrawLatex(0.18,y-=0.07,Form("V0 HM triggers:  %.0fM",gStatV0M->GetY()[gStatV0M->GetN()-1]));
  latex->SetTextColor(gStatSH2->GetLineColor());
  latex->DrawLatex(0.18,y-=0.07,Form("SPD HM triggers: %.0fM",gStatSH2->GetY()[gStatSH2->GetN()-1]));
//  latex->SetTextColor(gStatEMC->GetLineColor());
//  latex->DrawLatex(0.18,y-=0.07,Form("EMCAL triggers: %.0fM",gStatEMC->GetY()[gStatDMC->GetN()-1]));
//  latex->SetTextColor(gStatDMC->GetLineColor());
//  latex->DrawLatex(0.18,y-=0.07,Form("DCAL triggers: %.0fM",gStatDMC->GetY()[gStatDMC->GetN()-1]));
  gPad->Print("stat_mb_triggers.png");
}
