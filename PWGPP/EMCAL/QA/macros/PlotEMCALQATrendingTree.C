#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TStyle.h>
#include <TROOT.h>
#include <Riostream.h>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TTree.h>
#include <TTreeFormula.h>
#include <TEventList.h>
#include <TObjString.h>
#include <TList.h>
#include <TCut.h>
#include <TError.h>
#include <TDatime.h>
#include <TPRegexp.h>
#include <TInterpreter.h>
#include <TGaxis.h>
#endif

///
/// \file PlotEMCALQATrendingTree.C
/// \ingroup EMCALOfflineMacros
/// \brief QA at run level
///
/// This macro produces periodLevelQA for EMCAL from a trending tree.
/// Re-factored for automatic QA processing and trending by Alexandre Shabetai.
///
/// \author Alexis Mas, <Alexis.Mas@subatech.in2p3.fr>, SUBATECH
/// \author Alexandre Shabetai, <Alexandre.Shabetai@subatech.in2p3.fr>  SUBATECH
/// \author Marie Germain, <Marie.Germain@subatech.in2p3.fr>, SUBATECH
///


int PlotEMCALQATrendingTree(TTree* tree, const char* trig,TFile* fout,Bool_t SavePlots,TString expr);
TH1F* ZoomFromTree(TH1F* h, TTree* atree, Int_t n, const char* aVar, UShort_t aScaleFactor=1);
Double_t GetTreeMinimum(TTree* aTree, Int_t n, const char* columname);
Double_t GetTreeMaximum(TTree* aTree, Int_t n, const char* columname);

TString QAPATH;
TString QAPATHF= "./";

//--------------------------------------------------------------------------------------------------------------------
int NextInt(int newMax=0)
{

  static int N=0;
  static int nMax=1;

  if(newMax)  {nMax = newMax; N=0; return 0;}
  if (N == nMax) N=0;
  N = N +1;
  return N;

}

//--------------------------------------------------------------------------------------------------------------------
int PlotEMCALQATrendingTree(const char* filename="trending.root",Bool_t SavePlots=0, TString expr ="", TString fTrigger="")
{

  QAPATH = TString(gSystem->Getenv("QAPATH"));
  if(QAPATH.IsNull()) QAPATH = QAPATHF;

  Int_t ret=0;
  TFile* f = TFile::Open(filename);
  if(!f) { return -1;}
  TTree* tree = (TTree*)f->Get("trending");
  if (! tree) {Error("PlotEMCALQATrendingTree","No Tree found!"); return -1;}
  TFile* fout = new TFile(Form("%s/trendingPlots.root",QAPATH.Data()),"RECREATE");

  TList* TriggersList = new TList();
  if (fTrigger=="")
    {
      tree->Draw("fTrigger","","goff");
      const char* obj;
      for(Int_t i = 0 ; i < tree->GetSelectedRows() ; i++){
	tree->GetEntry(i);
	obj = tree->GetVar1()->PrintValue(0);
	if(! TriggersList->FindObject(obj)) {TriggersList->Add(new TObjString(obj));}
      }
    }
  else
    {
      if(!fTrigger.Contains("QA")) {fTrigger = "CaloQA_" + fTrigger;}
      TriggersList->Add(new TObjString(fTrigger.Data()));
    }
  TIter next(TriggersList);
  TObject* obj1;
  while ((obj1 = next()))
    {
      ret = PlotEMCALQATrendingTree(tree,obj1->GetName(),fout,SavePlots,expr);
    }

  f->Close();
  return ret;

}

//--------------------------------------------------------------------------------------------------------------------
int PlotEMCALQATrendingTree(TTree* tree, const char* Trig, TFile* fout, Bool_t SavePlots,TString Expr)
{

  TCut trig = Form("fTrigger==\"%s\"",Trig);
  TCut NotZero = TCut("Nevent>0.&&Npi0SM>0.");
  TCut select = trig;

  if (Expr.Contains(".C"))
    {
      Info("PlotEMCALQATrendingTree",Form("Additional selections from %s: ", Expr.Data()));
      gInterpreter->ExecuteMacro(Expr.Data());
      select = trig + expr;
   }

  if (! tree) {Error("PlotEMCALQATrendingTree","No Tree found!"); return -1;}
  select.Print();
  int CurN=0;
  TString* fCalorimeter;
  TString* period;
  TString* pass;
  TString* fTrigger;
  TString* system;
  TDatime* dtime;

  tree->SetBranchAddress("fDate",&dtime);
  tree->SetBranchAddress("nSM",&CurN);
  tree->SetBranchAddress("fCalorimeter",&fCalorimeter);
  tree->SetBranchAddress("system",&system);
  tree->SetBranchAddress("period",&period);
  tree->SetBranchAddress("pass",&pass);
  tree->SetBranchAddress("fTrigger",&fTrigger);


  tree->SetEventList(0);
  tree->Draw(">>elist",select);
  tree->Draw(">>listNotZero",select+NotZero);
  TEventList* listNotZero = (TEventList*)gDirectory->Get("listNotZero");
  TEventList* elist = (TEventList*)gDirectory->Get("elist");
  tree->SetEventList(elist);
  if(! elist->GetN()) { Error("PlotEMCALQATrendingTree","The current selection doess not match any entry!"); return -2; }
  CurN = tree->GetMinimum("nSM");
  const Int_t n = CurN;
  if(n<=12) const Int_t nEMCAL = n;
  if(n>12) const Int_t nEMCAL = 12;
  tree->GetEntry(elist->GetEntry(0));

  TGraphErrors* AverNclustersSM[n];
  TGraph* TotNclustersSM[n];
  TGraphErrors* AverNcellsPerClusterSM[n];
  TGraphErrors* AverESM[n];
  TGraphErrors* AverMeanSM[n];
  TGraphErrors* AverWidthSM[n];
  TGraphErrors* AverNpi0SM[n];

  // --------------------------------- plots ------------------------------

  TString base = QAPATH + period->Data() + "_" + pass->Data() + "_";
  TPRegexp r("_\\w+");

  TString ClusterAverages;          ClusterAverages         = base + "ClAv"        + (*fTrigger)(r) + ".png";
  TString Entries;                  Entries                 = base + "Nentries"    + (*fTrigger)(r) + ".png";
  TString ClusterAveragesEnergy;    ClusterAveragesEnergy   = base + "ClAvEne"     + (*fTrigger)(r) + ".png";
  TString ClusterAveragesEnergyD;   ClusterAveragesEnergyD   = base + "ClAvEne"     + (*fTrigger)(r) + "DCAL.png";
  TString ClusterAveragesEnergy2;   ClusterAveragesEnergy2  = base + "ClAvEne"     + (*fTrigger)(r) + ".pdf";
  TString ClusterAveragesEnergy2D;  ClusterAveragesEnergy2D  = base + "ClAvEne"     + (*fTrigger)(r) + "DCAL.pdf";
  TString ClusterAveragesEntries;   ClusterAveragesEntries  = base + "ClAvEnt"     + (*fTrigger)(r) + ".png";
  TString ClusterTotEntries;        ClusterTotEntries  = base + "ClTotEnt"     + (*fTrigger)(r) + ".png";
  TString ClusterTotEntries2;       ClusterTotEntries2  = base + "ClTotEnt"     + (*fTrigger)(r) + ".png";
  TString ClusterTotEntriesD;        ClusterTotEntries  = base + "ClTotEnt"     + (*fTrigger)(r) + ".png";
  TString ClusterTotEntriesD2;       ClusterTotEntries2  = base + "ClTotEnt"     + (*fTrigger)(r) + ".png";
  TString ClusterAveragesEntriesD;  ClusterAveragesEntriesD  = base + "ClAvEnt"     + (*fTrigger)(r) + "DCAL.png";
  TString ClusterAveragesEntries2;  ClusterAveragesEntries2 = base + "ClAvEnt"     + (*fTrigger)(r) + ".pdf";
  TString ClusterAveragesEntries2D; ClusterAveragesEntries2D = base + "ClAvEnt"     + (*fTrigger)(r) + "DCAL.pdf";
  TString ClusterAveragesCells;     ClusterAveragesCells    = base + "ClAvCells"   + (*fTrigger)(r) + ".png";
  TString ClusterAveragesCells2;    ClusterAveragesCells2   = base + "ClAvCells"   + (*fTrigger)(r) + ".pdf";
  TString ClusterChargedvsTot;      ClusterChargedvsTot     = base + "ClCharged"   + (*fTrigger)(r) + ".png";
  TString ClusterChargedvsTot2;     ClusterChargedvsTot2    = base + "ClCharged"   + (*fTrigger)(r) + ".pdf";
  TString Pi0Entries;               Pi0Entries              = base + "Pi0Entries"  + (*fTrigger)(r) + ".png";
  TString Pi0Entries2;              Pi0Entries2             = base + "Pi0Entries"  + (*fTrigger)(r) + ".pdf";
  TString Pi0Mass;                  Pi0Mass                 = base + "Pi0Mass"     + (*fTrigger)(r) + ".png";
  TString Pi0Mass2;                 Pi0Mass2                = base + "Pi0Mass"     + (*fTrigger)(r) + ".pdf";
  TString Pi0Width;                 Pi0Width                = base + "Pi0Width"    + (*fTrigger)(r) + ".png";
  TString Pi0Width2;                Pi0Width2               = base + "Pi0Width"    + (*fTrigger)(r) + ".pdf";


  int nEmptyRuns = tree->Draw("run","Nevent==0","goff");
  if (nEmptyRuns && (nEmptyRuns != -1)) {
   Info("PlotEMCALQATrendingTree",Form("The following %i runs are empty for trigger %s:",nEmptyRuns,Trig));
   for(Int_t i = 0 ; i < nEmptyRuns ; i++){
     cout<<tree->GetV1()[i]<<endl;
   }
  }

  int nNoEMCALRuns = tree->Draw("run","Nevent!=0&&EtotalMean==0","goff");
  if (nNoEMCALRuns && (nNoEMCALRuns != -1)) {
   Info("PlotEMCALQATrendingTree",Form("The following %i runs are without EMCAL for trigger %s:",nNoEMCALRuns,Trig));
   for(Int_t i = 0 ; i < nNoEMCALRuns ; i++){
     cout<<tree->GetV1()[i]<<endl;
   }
  }

  int nRun =  tree->Draw("run","","goff");
  NextInt(nRun);
  TH1F* h1 = new TH1F("h1", "dummy", nRun, 0., nRun+0.5);
  TGaxis::SetMaxDigits(3);
  h1->SetTitle("") ;
  h1->SetStats(kFALSE) ;
  h1->SetAxisRange(0, nRun, "X") ;
  h1->GetXaxis()->SetTitle("RUN Index");
  h1->GetXaxis()->SetTitleOffset(1.86);
  h1->GetXaxis()->SetTitleSize(0.03);

  for(Int_t i = 0 ; i < nRun ; i++){
    TString label = " ";
    label+=tree->GetV1()[i];
    h1->GetXaxis()->SetBinLabel(i+1,label.Data());
    h1->GetXaxis()->LabelsOption("v");
  }

  //number of events
  TCanvas* c1 = new TCanvas("Nevents","Nb of events", 1000, 500);
  c1->SetFillColor(0);
  c1->SetBorderSize(0);
  c1->SetFrameBorderMode(0);
  gStyle->SetOptStat(0);
  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.02);
  c1->SetGrid();
  tree->Draw("NextInt():Nevent","","goff");
  h1->GetYaxis()->SetTitle("N_{events}");
  ZoomFromTree(h1,tree,n,"Nevent",2);
  if (h1->GetMinimum() > 0.) {c1->SetLogy();}
  h1->Draw();

  TGraph* Nevents = new TGraph(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2());
  Nevents->SetMarkerStyle(20);
  Nevents->SetMarkerColor(1);
  Nevents->SetLineColor(2);
  Nevents->Draw("same lp") ;

  c1->Update();
  if(SavePlots) c1->SaveAs(Entries);

  TCanvas* c2 = new TCanvas("ClusterAveragesEvents", "Mean Nb of Cluster per Event", 1000, 500);
  c2->SetFillColor(0);
  c2->SetBorderSize(0);
  c2->SetFrameBorderMode(0);
  c2->SetGrid();

  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.02);
  gPad->SetGrid();

  TH1F* h2 = (TH1F*)h1->Clone("");
  h2->GetYaxis()->SetTitle("<N_{clusters}>/event");
  ZoomFromTree(h2,tree,n,"ClusterMeanSM");
  h2->GetXaxis()->SetTitle("RUN Index");
  h2->GetXaxis()->SetTitleOffset(1.86);
  h2->GetXaxis()->SetTitleSize(0.03);
  h2->Draw();

  tree->Draw("NextInt():ClusterMean:xe:ClusterRMS","","goff");
  TGraphErrors * AverNclusters = new TGraphErrors(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2(),tree->GetV3(), tree->GetV4());
  AverNclusters->SetMarkerStyle(20);
  AverNclusters->SetMarkerColor(1);
  AverNclusters->Draw("same P") ;

  for(Int_t ism = 0 ; ism < nEMCAL ; ism++){
    tree->Draw(Form("NextInt():ClusterMeanSM[%i]:xe:ClusterRMSSM[%i]",ism,ism),"","goff");
    AverNclustersSM[ism] = new  TGraphErrors(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2(),tree->GetV3(),tree->GetV4());
    if (ism !=8)AverNclustersSM[ism]->SetMarkerColor(ism<10?ism+2:ism+1);else AverNclustersSM[ism]->SetMarkerColor(7);
    AverNclustersSM[ism]->SetMarkerStyle(21+(ism<10 ? ism: ism-10));

    AverNclustersSM[ism]->Draw("same P");
  }
 TLegend* l2 = new TLegend(0.123, 0.744, 0.933, 0.894);
  l2->SetNColumns((n+1)/2.);
  l2->SetFillColor(0);
  l2->SetBorderSize(0);
  l2->SetTextSize(0.04);
  l2->SetHeader(Form("<# of clusters> in %s (period %s trigger %s)",fCalorimeter->Data(),period->Data(),((*fTrigger)(r)).Data()));
  l2->AddEntry(AverNclusters,"average", "p");
  for(Int_t ism = 0; ism < nEMCAL ; ism++){
    TString projname = Form("SM %d",ism);
    l2->AddEntry(AverNclustersSM[ism],projname.Data(), "p");
  }
  l2->Draw("same");
  c2->Update();
  if(SavePlots)     c2->SaveAs(ClusterAveragesEntries);
  if(SavePlots==2)  c2->SaveAs(ClusterAveragesEntries2);

  if(n>12)
    {
  TCanvas* c2D = new TCanvas("ClusterAveragesEventsD", "Mean Nb of Cluster per Event", 1000, 500);
  c2D->SetFillColor(0);
  c2D->SetBorderSize(0);
  c2D->SetFrameBorderMode(0);
  c2D->SetGrid();

  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.02);
  gPad->SetGrid();

  TH1F* h2D = (TH1F*)h1->Clone("");
  h2D->GetYaxis()->SetTitle("<N_{clusters}>/event");
  ZoomFromTree(h2D,tree,n,"ClusterMeanSM");
  h2D->GetXaxis()->SetTitle("RUN Index");
  h2D->GetXaxis()->SetTitleOffset(1.86);
  h2D->GetXaxis()->SetTitleSize(0.03);
  h2D->Draw();

  tree->Draw("NextInt():ClusterMean:xe:ClusterRMS","","goff");
  TGraphErrors * AverNclustersD = new TGraphErrors(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2(),tree->GetV3(), tree->GetV4());
  AverNclustersD->SetMarkerStyle(20);
  AverNclustersD->SetMarkerColor(1);
  AverNclustersD->Draw("same P") ;

  for(Int_t ism = 12 ; ism < n ; ism++){
    tree->Draw(Form("NextInt():ClusterMeanSM[%i]:xe:ClusterRMSSM[%i]",ism,ism),"","goff");
    AverNclustersSM[ism] = new  TGraphErrors(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2(),tree->GetV3(),tree->GetV4());
    if (ism !=8)AverNclustersSM[ism]->SetMarkerColor(ism<10?ism+2:ism+1);else AverNclustersSM[ism]->SetMarkerColor(7);
    AverNclustersSM[ism]->SetMarkerStyle(21+(ism<10 ? ism: ism-10));

    AverNclustersSM[ism]->Draw("same P");
  }

  TLegend* l2D = new TLegend(0.123, 0.744, 0.933, 0.894);
  l2D->SetNColumns((n+1)/2.);
  l2D->SetFillColor(0);
  l2D->SetBorderSize(0);
  l2D->SetTextSize(0.04);
  l2D->SetHeader(Form("<# of clusters> in %s (period %s trigger %s)",fCalorimeter->Data(),period->Data(),((*fTrigger)(r)).Data()));
  l2D->AddEntry(AverNclusters,"average", "p");
  for(Int_t ism = 12 ; ism < n ; ism++){
    TString projname = Form("SM %d",ism);
    l2D->AddEntry(AverNclustersSM[ism],projname.Data(), "p");
  }
  l2D->Draw("same");
  c2D->Update();
  if(SavePlots)     c2D->SaveAs(ClusterAveragesEntriesD);
  if(SavePlots==2)  c2D->SaveAs(ClusterAveragesEntries2D);
    }

  TCanvas* c2Tot = new TCanvas("ClusterTotEvents", "Tot Nb of Cluster per Event", 1000, 500);
  c2Tot->SetFillColor(0);
  c2Tot->SetBorderSize(0);
  c2Tot->SetFrameBorderMode(0);
  c2Tot->SetGrid();

  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.02);
  gPad->SetGrid();

  TH1F* h2Tot = (TH1F*)h1->Clone("");
  h2Tot->GetYaxis()->SetTitle("N_{clusters}/event");
  ZoomFromTree(h2Tot,tree,n,"ClusterTotSM");
  h2Tot->GetXaxis()->SetTitle("RUN Index");
  h2Tot->GetXaxis()->SetTitleOffset(1.86);
  h2Tot->GetXaxis()->SetTitleSize(0.03);
  h2Tot->Draw();



  for(Int_t ism = 0 ; ism < nEMCAL ; ism++){
    tree->Draw(Form("NextInt():ClusterTotSM[%i]",ism),"","goff");
    TotNclustersSM[ism] = new  TGraph(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2());
    if (ism !=8)TotNclustersSM[ism]->SetMarkerColor(ism<10?ism+2:ism+1);else TotNclustersSM[ism]->SetMarkerColor(7);
    TotNclustersSM[ism]->SetMarkerStyle(21+(ism<10 ? ism: ism-10));

    TotNclustersSM[ism]->Draw("same P");
  }
 TLegend* l2Tot = new TLegend(0.123, 0.744, 0.933, 0.894);
  l2Tot->SetNColumns((n+1)/2.);
  l2Tot->SetFillColor(0);
  l2Tot->SetBorderSize(0);
  l2Tot->SetTextSize(0.04);
  l2Tot->SetHeader(Form("# of clusters in %s (period %s trigger %s)",fCalorimeter->Data(),period->Data(),((*fTrigger)(r)).Data()));
  //l2Tot->AddEntry(TotNclusters,"total", "p");
  for(Int_t ism = 0; ism < nEMCAL ; ism++){
    TString projname = Form("SM %d",ism);
    l2Tot->AddEntry(TotNclustersSM[ism],projname.Data(), "p");
  }
  l2Tot->Draw("same");
  c2Tot->Update();
  if(SavePlots)     c2Tot->SaveAs(ClusterTotEntries);
  if(SavePlots==2)  c2Tot->SaveAs(ClusterTotEntries2);

  if(n>12)
    {
  TCanvas* c2TotD = new TCanvas("ClusterTotEventsD", "Tot Nb of Cluster per Event in DCAL", 1000, 500);
  c2TotD->SetFillColor(0);
  c2TotD->SetBorderSize(0);
  c2TotD->SetFrameBorderMode(0);
  c2TotD->SetGrid();

  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.02);
  gPad->SetGrid();

  TH1F* h2TotD = (TH1F*)h1->Clone("");
  h2TotD->GetYaxis()->SetTitle("N_{clusters}/event");
  ZoomFromTree(h2TotD,tree,n,"ClusterTotSM");
  h2TotD->GetXaxis()->SetTitle("RUN Index");
  h2TotD->GetXaxis()->SetTitleOffset(1.86);
  h2TotD->GetXaxis()->SetTitleSize(0.03);
  h2TotD->Draw();



  for(Int_t ism = 12 ; ism < n ; ism++){
    tree->Draw(Form("NextInt():ClusterTotSM[%i]",ism),"","goff");
    TotNclustersSM[ism] = new  TGraph(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2());
    if (ism !=8)TotNclustersSM[ism]->SetMarkerColor(ism<10?ism+2:ism+1);else TotNclustersSM[ism]->SetMarkerColor(7);
    TotNclustersSM[ism]->SetMarkerStyle(21+(ism<10 ? ism: ism-10));

    TotNclustersSM[ism]->Draw("same P");
  }
 TLegend* l2TotD = new TLegend(0.123, 0.744, 0.933, 0.894);
  l2TotD->SetNColumns((n+1)/2.);
  l2TotD->SetFillColor(0);
  l2TotD->SetBorderSize(0);
  l2TotD->SetTextSize(0.04);
  l2TotD->SetHeader(Form("# of clusters in %s (period %s trigger %s)",fCalorimeter->Data(),period->Data(),((*fTrigger)(r)).Data()));
  //l2Tot->AddEntry(TotNclusters,"total", "p");
  for(Int_t ism = 12; ism < n ; ism++){
    TString projname = Form("SM %d",ism);
    l2TotD->AddEntry(TotNclustersSM[ism],projname.Data(), "p");
  }
  l2TotD->Draw("same");
  c2TotD->Update();
  if(SavePlots)     c2TotD->SaveAs(ClusterTotEntriesD);
  if(SavePlots==2)  c2TotD->SaveAs(ClusterTotEntriesD2);

    }



  TCanvas* c3 = new TCanvas("ClusterAveragesEnergy", "Mean Cluster Energy", 1000, 500);
  c3->SetFillColor(0);
  c3->SetBorderSize(0);
  c3->SetFrameBorderMode(0);
  c3->SetGrid();

  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.02);
  gPad->SetGrid();

  TH1F* h3 = (TH1F*)h1->Clone("");
  h3->GetYaxis()->SetTitle("<E> (GeV)");
  ZoomFromTree(h3,tree,n,"EtotalMeanSM");
  h3->GetXaxis()->SetTitle("RUN Index");
  h3->GetXaxis()->SetTitleOffset(1.86);
  h3->GetXaxis()->SetTitleSize(0.03);
  h3->Draw();

  tree->Draw("NextInt():EtotalMean:xe:EtotalRMS","","goff");
  TGraphErrors * AverE = new TGraphErrors(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2(),tree->GetV3(),tree->GetV4());
  AverE->SetMarkerStyle(20);
  AverE->SetMarkerColor(1);
  AverE->Draw("same P");

  for(Int_t ism = 0 ; ism < n ; ism++){

    tree->Draw(Form("NextInt():EtotalMeanSM[%i]:xe:EtotalRMSSM[%i]",ism,ism),"","goff");
    AverESM[ism] = new  TGraphErrors(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2(),tree->GetV3(),tree->GetV4());
    if (ism !=8)AverESM[ism]->SetMarkerColor(ism<10?ism+2:ism+1);else AverESM[ism]->SetMarkerColor(7);
    AverESM[ism]->SetMarkerStyle(21+(ism<10 ? ism: ism-10));
    AverESM[ism]->Draw("same P");

  }

  TLegend* l3 = new TLegend(0.123, 0.744, 0.933, 0.894);
  l3->SetNColumns((n+1)/2.);
  l3->SetFillColor(0);
  l3->SetBorderSize(0);
  l3->SetTextSize(0.04);
  l3->SetHeader(Form("<E> in %s (period %s trigger %s)",fCalorimeter->Data(),period->Data(),((*fTrigger)(r)).Data()));
  l3->AddEntry(AverE,"average", "p");
  for(Int_t ism = 0 ; ism < n ; ism++){
    TString projname = Form("SM %d",ism);
    l3->AddEntry(AverESM[ism],projname.Data(), "p");
  }
  l3->Draw("same");

  if(SavePlots) c3->SaveAs(ClusterAveragesEnergy);
  if(SavePlots==2) c3->SaveAs(ClusterAveragesEnergy2);

  TCanvas* c4 = new TCanvas("ClusterAveragesCells", "Mean Nb of Cells per Cluster", 1000, 500);
  c4->SetFillColor(0);
  c4->SetBorderSize(0);
  c4->SetFrameBorderMode(0);
  c4->SetGrid();

  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.02);
  gPad->SetGrid();

  TH1F* h4 = (TH1F*)h1->Clone("");
  h4->GetYaxis()->SetTitle("<N_{CellsPerCluster}>");
  ZoomFromTree(h4,tree,n,"CellPerClusterMeanSM");
  h4->GetXaxis()->SetTitle("RUN Index");
  h4->GetXaxis()->SetTitleOffset(1.86);
  h4->GetXaxis()->SetTitleSize(0.03);
  h4->Draw();

  //
  tree->Draw("NextInt():CellPerClusterMean:xe:CellPerClusterRMS","","goff");
  TGraphErrors * AverCellPerCluster = new TGraphErrors(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2(),tree->GetV3(),tree->GetV4());
  AverCellPerCluster->SetMarkerStyle(20);
  AverCellPerCluster->SetMarkerColor(1);

  for(Int_t ism = 0 ; ism < n ; ism++){
    tree->Draw(Form("NextInt():CellPerClusterMeanSM[%i]:xe:CellPerClusterRMSSM[%i]",ism,ism),"","goff");
    AverNcellsPerClusterSM[ism] =  new TGraphErrors(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2(),tree->GetV3(),tree->GetV4());
    if (ism !=8)AverNcellsPerClusterSM[ism]->SetMarkerColor(ism<10?ism+2:ism+1);else AverNcellsPerClusterSM[ism]->SetMarkerColor(7);
    AverNcellsPerClusterSM[ism]->SetMarkerStyle(21+(ism<10 ? ism: ism-10));
    AverNcellsPerClusterSM[ism]->Draw("same P");

  }

  TLegend* l4 = new TLegend(0.123, 0.744, 0.933, 0.894);
  l4->SetNColumns((n+1)/2.);
  l4->SetFillColor(0);
  l4->SetBorderSize(0);
  l4->SetTextSize(0.04);
  l4->SetHeader(Form("<# of cells per cluster> in %s (period %s trigger %s)",fCalorimeter->Data(),period->Data(),((*fTrigger)(r)).Data()));
  l4->AddEntry(AverCellPerCluster,"average", "p");
  for(Int_t ism = 0 ; ism < n ; ism++){
    TString projname = Form("SM %d",ism);
    l4->AddEntry(AverNcellsPerClusterSM[ism],projname.Data(), "p");
  }
  l4->Draw("same");

  if(SavePlots) c4->SaveAs(ClusterAveragesCells);

    TCanvas* c8 = new TCanvas("NMatchClusters","x100 % of matched clusters", 1000, 500);
  c8->SetFillColor(0);
  c8->SetBorderSize(0);
  c8->SetFrameBorderMode(0);
  gStyle->SetOptStat(0);
  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.02);
  c8->SetGrid();


  TH1F* h8 = (TH1F*)h1->Clone("");
  h8->GetYaxis()->SetTitle("#frac{N_{match}}{N_{tot}}");
  ZoomFromTree(h8,tree,n,"NMatchClustersP");
  h8->GetXaxis()->SetTitle("RUN Index");
  h8->GetXaxis()->SetTitleOffset(1.86);
  h8->GetXaxis()->SetTitleSize(0.03);
  h8->Draw();

  tree->Draw("NextInt():NMatchClustersP:xe:NMatchClustersPRMS","","goff");
  TGraphErrors* NMatchCl = new TGraphErrors(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2(), tree->GetV3(), tree->GetV4());
  NMatchCl->SetMarkerStyle(20);
  NMatchCl->SetMarkerColor(1);
  NMatchCl->Draw("same p") ;

  c8->Update();
  if(SavePlots) c8->SaveAs(ClusterChargedvsTot);


  TCanvas* c5 = new TCanvas("Pi0Position", "Mean Pi0 Mass", 1000, 500);
  c5->SetFillColor(0);
  c5->SetBorderSize(0);
  c5->SetFrameBorderMode(0);
  c5->SetGrid();

  gStyle->SetOptStat(0);

  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.02);
  gPad->SetGrid();



  TLine* lUp = new TLine(0,148,nRun+0.5,148);
  lUp->SetLineColor(kRed);

  TLine* lDown = new TLine(0,122,nRun+0.5,122);
 lDown->SetLineColor(kRed);

  TH1F * h5 = (TH1F*)h1->Clone("");
  ZoomFromTree(h5,tree,n,"MeanPosSM");
  h5->GetXaxis()->SetTitle("RUN Index");
  h5->GetXaxis()->SetTitleOffset(1.86);
  h5->GetXaxis()->SetTitleSize(0.03);
  h5->GetYaxis()->SetTitle("Mean_{#pi^{0}}");

  h5->Draw();


  tree->Draw("NextInt():MeanPosEMCAL:xe:MeanPosEMCALErr","","goff");
  TGraphErrors * AverMeanEMCAL = new TGraphErrors(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2(),tree->GetV3(),tree->GetV4());
  AverMeanEMCAL->SetMarkerStyle(20);
  AverMeanEMCAL->SetMarkerColor(1);
  AverMeanEMCAL->Draw("same P");

  if(n>12)
    {
  tree->Draw("NextInt():MeanPosDCAL:xe:MeanPosDCALErr","","goff");
  TGraphErrors * AverMeanDCAL = new TGraphErrors(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2(),tree->GetV3(),tree->GetV4());
  AverMeanDCAL->SetMarkerStyle(20);
  AverMeanDCAL->SetMarkerColor(2);
  AverMeanDCAL->Draw("same P");
    }
  // // TLine lineLow(V1[0], 120., V1[GetSelectedRows()], 145.);
  // //  TLine lineUp(V1(1), y_0, x_1, y_1);
  // lineLow->SetLineColor(2);
  // //  lineUp->SetLineColor(2);
  // lineLow->Draw("same");
  // // lineUp->Draw("same");

  for(Int_t ism = 0 ; ism < n ; ism++){

    tree->Draw(Form("NextInt():MeanPosSM[%i]:xe:MeanPosErrSM[%i]",ism,ism),"","goff");
    AverMeanSM[ism] = new TGraphErrors(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2(),tree->GetV3(),tree->GetV4());
    if (ism !=8)AverMeanSM[ism]->SetMarkerColor(ism<10?ism+2:ism+1);else AverMeanSM[ism]->SetMarkerColor(7);
    AverMeanSM[ism]->SetMarkerStyle(21+(ism<10 ? ism: ism-10));
    AverMeanSM[ism]->Draw("same P");
  }


  TLegend* l5 = new TLegend(0.123, 0.744, 0.933, 0.894);
  l5->SetNColumns((n+1)/2.);
  l5->SetFillColor(0);
  l5->SetBorderSize(0);
  l5->SetTextSize(0.04);
  l5->SetHeader(Form("<M_{#pi^{0}}> (MeV) in %s (period %s trigger %s)",fCalorimeter->Data(),period->Data(),((*fTrigger)(r)).Data()));
  l5->AddEntry(AverMeanEMCAL,"average EMCAL", "p");
  if(n>12) l5->AddEntry(AverMeanDCAL,"average DCAL", "p");
  for(Int_t ism = 0 ; ism < n ; ism++){
    TString projname = Form("SM %d",ism);
    l5->AddEntry(AverMeanSM[ism],projname.Data(), "p");
  }
  l5->Draw("same");
  lUp->Draw("same");
  lDown->Draw("same");

  c5->Update();
  if(SavePlots) c5->SaveAs(Pi0Mass);


  TCanvas* c6 = new TCanvas("Pi0Width", "Mean Pi0 Width", 1000, 500);
  c6->SetFillColor(0);
  c6->SetBorderSize(0);
  c6->SetFrameBorderMode(0);
  c6->SetGrid();

  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.02);
  gPad->SetGrid();

  TH1F* h6 = (TH1F*)h1->Clone("");
  ZoomFromTree(h6,tree,n,"WidthSM");
  h6->GetXaxis()->SetTitle("RUN Index");
  h6->GetXaxis()->SetTitleOffset(1.86);
  h6->GetXaxis()->SetTitleSize(0.03);
  h6->GetYaxis()->SetTitle("#sigma_{#pi^{0}}");
  h6->Draw();

  tree->Draw("NextInt():WidthEMCAL:xe:WidthEMCALErr","","goff");
  TGraphErrors * AverWidthEMCAL = new TGraphErrors(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2(),tree->GetV3(),tree->GetV4());
  AverWidthEMCAL->SetMarkerStyle(20);
  AverWidthEMCAL->SetMarkerColor(1);
  AverWidthEMCAL->Draw("same P");

  if(n>12)
    {
  tree->Draw("NextInt():WidthDCAL:xe:WidthDCALErr","","goff");
  TGraphErrors * AverWidthDCAL = new TGraphErrors(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2(),tree->GetV3(),tree->GetV4());
  AverWidthDCAL->SetMarkerStyle(20);
  AverWidthDCAL->SetMarkerColor(2);
  AverWidthDCAL->Draw("same P");
    }

  for(Int_t ism = 0 ; ism < n ; ism++){
    tree->Draw(Form("NextInt():WidthSM[%i]:xe:WidthErrSM[%i]",ism,ism),"","goff");
    AverWidthSM[ism] = new TGraphErrors(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2(),tree->GetV3(),tree->GetV4());
    if (ism !=8)AverWidthSM[ism]->SetMarkerColor(ism<10?ism+2:ism+1);else AverWidthSM[ism]->SetMarkerColor(7);
    AverWidthSM[ism]->SetMarkerStyle(21+(ism<10 ? ism: ism-10));
    AverWidthSM[ism]->Draw("same P");
  }


  TLegend* l6 = new TLegend(0.123, 0.744, 0.933, 0.894);
  l6->SetNColumns((n+1)/2.);
  l6->SetFillColor(0);
  l6->SetBorderSize(0);
  l6->SetTextSize(0.04);
  l6->SetHeader(Form("#sigma_{#pi^{0}} in %s (period %s trigger %s)",fCalorimeter->Data(),period->Data(),((*fTrigger)(r)).Data()));
  l6->AddEntry(AverWidthEMCAL,"total EMCAL", "p");
  if(n>12) l6->AddEntry(AverWidthDCAL,"total DCAL", "p");
  for(Int_t ism = 0 ; ism < n ; ism++){
    TString projname = Form("SM %d",ism);
    l6->AddEntry(AverWidthSM[ism],projname.Data(), "p");
  }
  l6->Draw("same");
  c6->Update();
  if(SavePlots) c6->SaveAs(Pi0Width);

  TCanvas* c7 = new TCanvas("Npi0", "Mean Nb of Pi0", 1000, 500);
  c7->SetFillColor(0);
  c7->SetBorderSize(0);
  c7->SetFrameBorderMode(0);
  c7->SetGrid();

  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.02);
  gPad->SetGrid();

  TH1F* h7 = (TH1F*)h1->Clone("");
  ZoomFromTree(h7,tree,n,"Npi0SM");
  if (h7->GetMinimum() > 0.) {c7->SetLogy();}
  h7->GetXaxis()->SetTitle("RUN Index");
  h7->GetXaxis()->SetTitleOffset(1.86);
  h7->GetXaxis()->SetTitleSize(0.03);
  h7->GetYaxis()->SetTitle("<N_{#pi^{0}}>/event");
  h7->Draw();

  tree->Draw("NextInt():Npi0EMCAL:xe:Npi0EMCALErr","","goff");
  if (tree->GetMinimum("Npi0EMCAL") > 1) c4->SetLogy();
  TGraphErrors * AverNpi0EMCAL = new TGraphErrors(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2(),tree->GetV3(),tree->GetV4());
  AverNpi0EMCAL->SetMarkerStyle(20);
  AverNpi0EMCAL->SetMarkerColor(1);
  AverNpi0EMCAL->Draw("same P");

  if(n>12)
  {
  tree->Draw("NextInt():Npi0DCAL:xe:Npi0DCALErr","","goff");
  if (tree->GetMinimum("Npi0DCAL") > 1) c4->SetLogy();
  TGraphErrors * AverNpi0DCAL = new TGraphErrors(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2(),tree->GetV3(),tree->GetV4());
  AverNpi0DCAL->SetMarkerStyle(20);
  AverNpi0DCAL->SetMarkerColor(2);
  AverNpi0DCAL->Draw("same P");
  }

  for(Int_t ism = 0 ; ism < n ; ism++){
    tree->Draw(Form("NextInt():Npi0SM[%i]:xe:Npi0ErrSM[%i]",ism,ism),"","goff");
    AverNpi0SM[ism] = new  TGraphErrors(tree->GetSelectedRows(), tree->GetV1(), tree->GetV2(),tree->GetV3(),tree->GetV4());
    if (ism !=8)AverNpi0SM[ism]->SetMarkerColor(ism<10?ism+2:ism+1);else AverNpi0SM[ism]->SetMarkerColor(7);
    AverNpi0SM[ism]->SetMarkerStyle(21+(ism<10 ? ism: ism-10));
    AverNpi0SM[ism]->Draw("same P");
  }

  TLegend* l7 = new TLegend(0.123, 0.744, 0.933, 0.894);
  l7->SetNColumns((n+1)/2.);
  l7->SetFillColor(0);
  l7->SetBorderSize(0);
  l7->SetTextSize(0.04);
  l7->SetHeader(Form("<N_{#pi^{0}}>/event in %s (period %s trigger %s)",fCalorimeter->Data(),period->Data(),((*fTrigger)(r)).Data()));
  l7->AddEntry(AverNpi0EMCAL,"total EMCAL", "p");
  if(n>12) l7->AddEntry(AverNpi0DCAL,"total DCAL", "p");
  for(Int_t ism = 0 ; ism < n ; ism++){
    TString projname = Form("SM %d",ism);
    l7->AddEntry(AverNpi0SM[ism],projname.Data(), "p");
  }
  l7->Draw("same");
  c7->Update();
  if(SavePlots)  c7->SaveAs(Pi0Entries);
  if(SavePlots==2) c7->SaveAs(Pi0Entries2);

  fout->mkdir(Form("%s/%s/%s/%s",period->Data(),pass->Data(),"TrendingQA",fTrigger->Data()));
  fout->cd();
  fout->Cd(Form("%s/%s/%s/%s",period->Data(),pass->Data(),"TrendingQA",fTrigger->Data()));

  gROOT->GetListOfCanvases()->Write();
  gROOT->GetListOfCanvases()->Delete();

  if((!Expr.IsNull()) && (!Expr.EndsWith(".root"))) elist->Write();
  if(listNotZero) {listNotZero->Reset();}
  if(elist) {elist->Reset();}
  delete h1;

  return 0;

}

//---------------------------------------------------------------------------------------------
TH1F* ZoomFromTree(TH1F* h, TTree* atree, Int_t n, const char* aVar, UShort_t  aScaleFactor)
{

  atree->SetEventList(0);
  TEventList *listNotZero = (TEventList*)gDirectory->Get("listNotZero");
  atree->SetEventList(listNotZero);

  double treeMin = GetTreeMinimum(atree,n,aVar);
  double treeMax = GetTreeMaximum(atree,n,aVar);
  double offset =  30*((treeMax - treeMin)/(100.*aScaleFactor));

  if(treeMin != -treeMax){
   h->SetMinimum(TMath::Max(0.,treeMin-offset));
   h->SetMaximum(treeMax+2*offset);
  }

  atree->SetEventList(0);
  TEventList *elist = (TEventList*)gDirectory->Get("elist");
  atree->SetEventList(elist);

  return h;

}

//--------------------------------------------------------------------------------------------
Double_t GetTreeMinimum(TTree* aTree,Int_t n, const char* columname)
{

  TLeaf* leaf = aTree->GetLeaf(columname);
  if (!leaf) {
    return 0;
  }
  TBranch* branch = leaf->GetBranch();
  Double_t cmin = 3.40282e+38;
  for (Long64_t i = 0; i < aTree->GetEntries(); ++i) {
    Long64_t entryNumber = aTree->GetEntryNumber(i);
    if (entryNumber < 0) break;
    branch->GetEntry(entryNumber);
    for (Int_t j = 0;j < TMath::Min(leaf->GetLen(),n); ++j) {
      Double_t val = leaf->GetValue(j);
      if (val < cmin) {
	cmin = val;
      }
    }
  }

  return cmin;

}

//______________________________________________________________________________
Double_t GetTreeMaximum(TTree* aTree,Int_t n,const char* columname)
{

  TLeaf* leaf = aTree->GetLeaf(columname);
  if (!leaf) {
    return 0;
  }
  TBranch* branch = leaf->GetBranch();
  Double_t cmax = - 3.40282e+38;
  for (Long64_t i = 0; i < aTree->GetEntries(); ++i) {
    Long64_t entryNumber = aTree->GetEntryNumber(i);
    if (entryNumber < 0) break;
    branch->GetEntry(entryNumber);
    for (Int_t j = 0; j <  TMath::Min(leaf->GetLen(),n); ++j) {
      Double_t val = leaf->GetValue(j);
      if (val > cmax) {
	cmax = val;
      }
    }
  }

  return cmax;

}

