/*
  Origin:  marian.ivanov@cern.ch
  Make sys watch default plots (see $ALICE_ROOT/STEER/AliSysInfo.cxx):
  Input   -  syswatch.log  - text log file created by process to be monitored
  Output  -  syswatch.root - root files with default histograms
  Number of top violators - only top consumer displayed
  Default histogram:  
  TOP violators      - CPU and Virtual memory usage
  Detector reports    - CPU and Virtual memory usage per detector
  //
  
  Example usage:
  //  1. Initialize
  gROOT->LoadMacro("$ALICE_ROOT/macros/PlotSys.C+");
  PInit("syswatch.log","syswatch.root","syswatch.sum");
  //  2. Make ascii report.
  SumDetector()
  //  3. Make histos of top violators
   MakePlots(20);
  //  4. Browse the results
  TFile f("syswatch.root");
  TBrowser b;
 */

#include <stdio.h>
#include "AliReconstruction.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TCut.h"
#include "TStyle.h"
#include "AliSysInfo.h"

const Int_t kNDetectors=AliReconstruction::kNDetectors;
const char* fgkDetectorName[kNDetectors] = {"ITS", "TPC", "TRD", "TOF", "PHOS", "HMPID", "EMCAL", "MUON", "FMD", "ZDC", "PMD", "T0", "VZERO", "ACORDE", "HLT"};



TObject * htemp; 
TTree *tree=0;
TFile *fout=0;
TString sumFile;
TCut cutVM("cutVM","deltaVM>10");
TCut cutDT("cutDT","deltaT>2"); 
Int_t ctop=20;


Float_t TopUsage(TTree* tree, const char *exp, const char*cut, Int_t order);
Double_t SumUsage(TTree* tree, const char *exp, const char*cut);
void TopVM();
void TopCPU();
void TopVMDetector();
void TopCPUDetector();
void SumDetector();

void PInit(const char *log="syswatch.log", const char *out="syswatch.root", const char * sumName="syswatch.sum"){
  //
  // Set Input output
  //
  tree = AliSysInfo::MakeTree(log);
  fout = new TFile(out,"recreate");
  sumFile=sumName;
}



void MakePlots(Int_t top=20){
  //
  //
  //
  ctop=top;
  gStyle->SetOptStat(0);
  //
  // Top users
  //
  TopVM();
  TopCPU();
  //
  // Reports per detector
  //
  fout->cd();
  for (Int_t idet=0; idet<kNDetectors;idet++){
    fout->cd();
    fout->mkdir(fgkDetectorName[idet]);
  }
  TopVMDetector();
  TopCPUDetector();
  //
  fout->Close();
  ctop=top;
  delete fout;
}

void TopVM(){
  //
  // select top user of virtual Memory 
  // MakeReport - ASCII and histogram
  // 
  TH1 * his=0;
  TH2 * his2=0;
  Float_t thVM = TopUsage(tree,"deltaVM","",ctop);
  cutVM = TCut("cutDT",Form("deltaVM>%f",thVM));
  //
  //
  printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n\n");
  printf("TOP Virtual memory user\n");
  //tree->Scan("deltaVM:sname",cutVM,"colsize=20");
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
  //
  // select top user of CPU 
  // MakeReport - ASCII and histogram
  // 
  TH2 * his2=0;
  Float_t thDT = TopUsage(tree,"deltaT","id2<3",ctop);
  cutDT = TCut("cutDT",Form("deltaT>%f",thDT));
  //
  printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
  printf("/n/n/nTOP CPU user\n");
  //tree->Scan("deltaT:sname",cutDT,"colsize=20");
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
  //
  // Draw usage of VM
  //
  TH2 * his2=0;
  //
  //
  //detector part
  //
  for (Int_t idet=0; idet<kNDetectors; idet++){
    fout->cd();
    fout->cd(fgkDetectorName[idet]);
    char cdet[100];
    char cdvm[100];
    sprintf(cdet,"id0==%d",idet);
    char expr[100];
    sprintf(expr,"deltaVM:sname>>hhh");  
    //
    Float_t thDVM = TopUsage(tree,"deltaVM",cdet,ctop);
    sprintf(cdvm,"%s&&deltaVM>%f",cdet, thDVM);
    //
    tree->Draw(expr,cdvm,"GOFF"); 
    his2 = (TH2F*)(tree->GetHistogram())->Clone("xxx");
    delete tree->GetHistogram();
    his2->SetYTitle("Delta Virtual Memory (MBy)");
    his2->SetMarkerStyle(22);
    his2->SetMarkerSize(1); 
    //his2->Draw("l*");
    his2->Write(Form("DVMvsName_%d",idet));
    delete his2;
    //
    //    
    sprintf(expr,"VM:sname>>hhh");
    tree->Draw(expr,cdvm,"goff"); 
    his2 = (TH2F*)(tree->GetHistogram())->Clone("yyy");
    delete tree->GetHistogram();
    his2->SetYTitle("Delta Virtual Memory (MBy)");
    his2->SetMarkerStyle(22);
    his2->SetMarkerSize(1); 
    //his2->Draw("l*");
    his2->Write(Form("VMvsName_%d",idet));     
    delete his2;
  }
  fout->cd();
}



void TopCPUDetector(){
  //
  // Draw usage of CPU
  //
  TH2 * his2=0;
  //
  //
  // CPU
  //
  for (Int_t idet=0; idet<kNDetectors; idet++){
    fout->cd();
    fout->cd(fgkDetectorName[idet]);
    char cdet[100];
    char cdtime[100];
    sprintf(cdet,"id0==%d",idet);
    char expr[100];
    sprintf(expr,"deltaT:sname>>hhh");  
    //
    Float_t thDT = TopUsage(tree,"deltaT",cdet,ctop);
    sprintf(cdtime,"%s&&deltaT>%f",cdet, thDT);
    //
    tree->Draw(expr,cdtime,"goff"); 
    his2 = (TH2F*)(tree->GetHistogram())->Clone("dtsname");
    delete tree->GetHistogram();
    his2->SetYTitle("Delta CPU time(sec)");
    his2->SetMarkerStyle(22);
    his2->SetMarkerSize(1); 
    his2->GetXaxis()->SetLabelSize(0.03);
    //his2->Draw("l*");
    his2->Write(Form("CPUvsName_%d",idet));
    delete his2;
  }
  fout->cd();
}

void SumDetector(){
  //
  // Sum - detector information
  //
  FILE * pFile;
  pFile = fopen (sumFile,"w");
  char cdet[100];
  char expr[100];
  sprintf(cdet,"id0>=0&&id2>=0");
  Double_t sumdTAll  = SumUsage(tree,"deltaT",cdet);
  Double_t sumdVMAll = SumUsage(tree,"deltaVM",cdet);
  printf("%s%s%s%s%s\n","Det/C:","sumDt/F:","sumDvm/F:","fracDt/F:","fracDvm/F");
  printf("%s\t%f\t%f\t%f\t%f\t\n","all", sumdTAll,sumdVMAll,100.,100.);
  fprintf(pFile,"%s%s%s%s%s\n","Det/C:","sumDt/F:","sumDvm/F:","fracDt/F:","fracDvm/F");
  fprintf(pFile,"%s\t%f\t%f\t%f\t%f\t\n","all", sumdTAll,sumdVMAll,100.,100.);
  for (Int_t idet=0; idet<kNDetectors; idet++){
    sprintf(cdet,"id0==%d&&id2>=0",idet);
    sprintf(expr,"deltaT:sname>>hhh");  
    Double_t sumdT  = SumUsage(tree,"deltaT",cdet);
    Double_t sumdVM = SumUsage(tree,"deltaVM",cdet);
    printf("%s\t%f\t%f\t%f\t%f\t\n",fgkDetectorName[idet], sumdT,sumdVM,100.*sumdT/sumdTAll, 100.*sumdVM/sumdVMAll); 
    fprintf(pFile,"%s\t%f\t%f\t%f\t%f\t\n",fgkDetectorName[idet], sumdT,sumdVM,100.*sumdT/sumdTAll, 100.*sumdVM/sumdVMAll); 
  }
  fclose (pFile);
}






Float_t TopUsage(TTree* tree, const char *exp, const char*cut, Int_t order){
  //
  // 
  // Find value for given order
  // Used to select top violator
  //
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


Double_t SumUsage(TTree* tree, const char *exp, const char*cut){
  //
  // return sum of usage
  //
  Int_t  entries = tree->Draw(Form("%s",exp),cut,"goff");
  if (entries==0) return 0;
  Double_t mean = TMath::Mean(entries, tree->GetV1());
  return entries*mean;
}
