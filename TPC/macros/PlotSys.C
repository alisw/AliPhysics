/// \file PlotSys.C
/// \author marian.ivanov@cern.ch
/// \brief Make sys watch default plots
///
/// See `$ALICE_ROOT/STEER/AliSysInfo.cxx`.
/// Input   -  syswatch.log  - text log file created by process to be monitored
/// Output  -  syswatch.root - root files with default histograms
/// Number of top violators - only top consumer displayed
///         
/// Default histogram:  
///              
/// TOP violateors      - CPU and Virtual memory usage
/// Detector reports    - CPU and Virtual memory usage per detector
/// 
/// 
/// 
/// 
/// 
/// Usage example:
/// ~~~{.cpp}
/// .x ~/rootlogon.C
/// gROOT->LoadMacro("$ALICE_ROOT/macros/PlotSys.C+");
/// MakePlots("syswatch.log","syswatch.root",10);
/// TFile f("syswatch.root");
/// TBrowser b;
/// ~~~

#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TCut.h"
#include "TStyle.h"
#include "AliSysInfo.h"

TObject * htemp; 
TTree *tree=0;
TFile *fout=0;
TCut cutVM("cutVM","deltaVM>10");
TCut cutDT("cutDT","deltaT>2"); 
Int_t ctop=10;


Float_t TopUsage(TTree* tree, const char *exp, const char*cut, Int_t order);
void TopVM();
void TopCPU();
void TopVMDetector();
void TopCPUDetector();

void PInit(const char *log="syswatch.log", const char *out="syswatch.root"){
  /// Set Input output

  tree = AliSysInfo::MakeTree(log);
  fout = new TFile(out,"recreate");
}



void MakePlots(const char *log="syswatch.log", const char *out="syswatch.root", Int_t top=10){
  ///

  ctop=top;
  PInit(log,out);
  gStyle->SetOptStat(0);
  //
  // Top users
  //
  TopVM();
  TopCPU();
  //
  // Reports per detector
  //
  fout->mkdir("cpuDetector");
  fout->mkdir("VMDetector");
  //
  fout->cd("VMDetector");
  TopVMDetector();
  //
  fout->cd();
  fout->cd("cpuDetector");
  TopCPUDetector();

  //
  fout->Close();
  ctop=top;
  delete fout;
}

void TopVM(){
  /// select top user of virtual Memory
  /// MakeReport - ASCII and histogram

  TH1 * his=0;
  TH2 * his2=0;
  Float_t thVM = TopUsage(tree,"deltaVM","",ctop);
  cutVM = TCut("cutDT",Form("deltaVM>%f",thVM));
  //
  //
  printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n\n");
  printf("TOP Virtual memory user\n");
  tree->Scan("deltaVM:sname",cutVM,"colsize=20");
  printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n\n");
  //
  tree->Draw("deltaVM:sname>>hhh","1"+cutVM,"*");  
  his2 = (TH2F*)(tree->GetHistogram())->Clone("dvmsname");
  delete tree->GetHistogram();
  his2->SetYTitle("Delta Virtual Memory (MBy)");
  his2->SetMarkerStyle(22);
  his2->SetMarkerSize(1); 
  his2->Draw("l*");
  his2->Write("DVMvsName");
  delete his2;
  //
  tree->Draw("VM:sname>>hhh","id2<3"+cutVM,"*");
  his2 = (TH2F*)(tree->GetHistogram())->Clone("vmsname");
  delete tree->GetHistogram();
  his2->SetYTitle("Delta Virtual Memory (MBy)");
  his2->SetMarkerStyle(22);
  his2->SetMarkerSize(1); 
  his2->Draw("l*");
  his2->Write("VMvsName");
  delete his2;
  //
  //
  tree->Draw("VM:T>>hhh","deltaVM>1","line*");
  his = (TH1*)tree->GetHistogram()->Clone("vmt");
  delete tree->GetHistogram();
  his->SetXTitle("Time (sec)");
  his->SetYTitle("Virtual Memory (MBy)");
  his->GetYaxis()->SetTitleOffset(1.2); 
  his->SetMarkerStyle(22);
  his->SetMarkerSize(1); 
  his->Draw();
  his->Write("VMvsTime");
  delete his;
}

void TopCPU(){  
  /// select top user of CPU
  /// MakeReport - ASCII and histogram

  TH2 * his2=0;
  Float_t thDT = TopUsage(tree,"deltaT","id2<3",ctop);
  cutDT = TCut("cutDT",Form("deltaT>%f",thDT));
  //
  printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
  printf("/n/n/nTOP CPU user\n");
  tree->Scan("deltaT:sname",cutDT,"colsize=20");
  printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
  //
  tree->Draw("deltaT:sname>>hhh","id2<3"+cutDT,"*"); 
  his2 = (TH2F*)(tree->GetHistogram())->Clone("tsname"); 
  delete tree->GetHistogram();
  his2->SetName("VMsanme");
  his2->SetYTitle("Delta CPU time(sec)");
  his2->SetMarkerStyle(22);
  his2->SetMarkerSize(1); 
  his2->GetXaxis()->SetLabelSize(0.03);
  his2->Draw("l*");
  his2->Write("CPUvsName");
  delete his2;
}


void TopVMDetector(){
  /// Draw usage of VM

  TH2 * his2=0;
  //
  //
  //detector part
  //
  for (Int_t idet=0; idet<12; idet++){
    char cdet[100];
    char cdvm[100];
    sprintf(cdet,"id0==%d",idet);
    char expr[100];
    sprintf(expr,"deltaVM:sname>>hhh");  
    //
    Float_t thDVM = TopUsage(tree,"deltaVM",cdet,ctop);
    sprintf(cdvm,"%s&&deltaT>%f",cdet, thDVM);


    //
    tree->Draw(expr,cdvm,"*"); 
    his2 = (TH2F*)(tree->GetHistogram())->Clone("xxx");
    delete tree->GetHistogram();
    his2->SetYTitle("Delta Virtual Memory (MBy)");
    his2->SetMarkerStyle(22);
    his2->SetMarkerSize(1); 
    his2->Draw("l*");
    his2->Write(Form("DVMvsName_%d",idet));
    delete his2;
    //
    //    
    sprintf(expr,"VM:sname>>hhh");
    tree->Draw(expr,cdvm,"*"); 
    his2 = (TH2F*)(tree->GetHistogram())->Clone("yyy");
    delete tree->GetHistogram();
    his2->SetYTitle("Delta Virtual Memory (MBy)");
    his2->SetMarkerStyle(22);
    his2->SetMarkerSize(1); 
    his2->Draw("l*");
    his2->Write(Form("VMvsName_%d",idet));     
    delete his2;
  }
}



void TopCPUDetector(){
  /// Draw usage of CPU

  TH2 * his2=0;
  //
  //
  // CPU
  //
  for (Int_t idet=0; idet<12; idet++){
    char cdet[100];
    char cdtime[100];
    sprintf(cdet,"id0==%d",idet);
    char expr[100];
    sprintf(expr,"deltaT:sname>>hhh");  
    //
    Float_t thDT = TopUsage(tree,"deltaT",cdet,ctop);
    sprintf(cdtime,"%s&&deltaT>%f",cdet, thDT);
    //
    tree->Draw(expr,cdtime,"*"); 
    his2 = (TH2F*)(tree->GetHistogram())->Clone("dtsname");
    delete tree->GetHistogram();
    his2->SetYTitle("Delta CPU time(sec)");
    his2->SetMarkerStyle(22);
    his2->SetMarkerSize(1); 
    his2->GetXaxis()->SetLabelSize(0.03);
    his2->Draw("l*");
    his2->Write(Form("CPUvsName_%d",idet));
    delete his2;
  }
}

 




Float_t TopUsage(TTree* tree, const char *exp, const char*cut, Int_t order){
  /// Find value for given order
  /// Used to select top violator

  Int_t entries = tree->Draw(Form("%s>>hhh1",exp),cut,"goff");
  if (entries<=1) {
    if (tree->GetHistogram()) delete tree->GetHistogram(); 
    printf("%s\t No entries\n",cut);
    return -10000;
  }
  if (!tree->GetV1()) {
    printf("%s\t No entries\n",cut);
    return -10000; 
  }
  Int_t *index = new Int_t[entries];
  TMath::Sort(entries, tree->GetV1(), index);
  Int_t oindex = TMath::Min(order, entries);
  Float_t val = tree->GetV1()[index[oindex-1]];
  if (tree->GetHistogram()) delete tree->GetHistogram();
  delete [] index;
  return val;
}
