#include "TList.h"
#include "TTree.h"
#include "AliSysInfo.h"
#include "TObjString.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TFile.h"

const Int_t NMODS = 14;
TString folder;
const char *module_name[NMODS] = {"SPD_Performance",
                                  "CaloQA",
                                  "Vertex_Performance",
                                  "PWGPP_QAsymHists",
                                  "VZERO_Performance",
                                  "TPC_PerformanceQA",
                                  "SDD_Performance",
                                  "PWGPPdEdxSSDQA",
                                  "ITS_Performance",
                                  "TRD_Performance",
                                  "MTR_ChamberEffMap",
                                  "PWG2forwardDnDeta",
                                  "ImpParRes_Performance",
                                  "MUON_QA"
};                                  

void SetFolder(const char *task_name) 
{
   if (!strcmp(task_name, "SPD_Performance")) folder = "coutput1";
   if (!strcmp(task_name, "CaloQA")) folder = "CaloQA";
   if (!strcmp(task_name, "Vertex_Performance")) folder = "cOutputVtxESD";
   if (!strcmp(task_name, "PWGPP_QAsymHists")) folder = "QAsymHists_Global QAsymHists_ITS QAsymHists_ITS_SA QAsymHists_TPC";
   if (!strcmp(task_name, "VZERO_Performance")) folder = "QAVZEROHists";
   if (!strcmp(task_name, "TPC_PerformanceQA")) folder = "TPCQA";
   if (!strcmp(task_name, "SDD_Performance")) folder = "coutputRP";
   if (!strcmp(task_name, "PWGPPdEdxSSDQA")) folder = "SSDdEdxQA";
   if (!strcmp(task_name, "ITS_Performance")) folder = "cOutputITS";
   if (!strcmp(task_name, "TRD_Performance")) folder = "checkESD infoGen checkDET TRDefficiency TRDresolution checkPID";
   if (!strcmp(task_name, "MTR_ChamberEffMap")) folder = "testMTRChamberEff triggerChamberEff";
   if (!strcmp(task_name, "PWG2forwardDnDeta")) folder = "BackgroundCorrected";
   if (!strcmp(task_name, "ImpParRes_Performance")) folder = "coutputd0ITSpureSARec coutputd0ITSpureSASkip coutputd0allPointRec coutputd0allPointSkip coutputd0partPointRec coutputd0partPointSkip coutputd0onepointSPDRec coutputd0onepointSPDSkip coutputd0postvTracRec coutputd0postvTracSkip coutputd0negtvTracRec coutputd0negtvTracSkip coutputd0pullAllpointRec coutputd0pullAllpointSkip coutputd0onlyRefitRec coutputd0onlyRefitSkip coutputd0sinThetaRec coutputd0sinThetaSkip coutputd0allPointTrue coutputd0postvTracTrue coutputd0negtvTracTrue coutputd0pullAllpointTrue coutputd0phiAllpointSkip coutputd0phiPostvtracSkip coutputd0phiNegtvtracSkip coutputd0clusterTypeSPD01Skip coutputd0clusterTypeSPD02Skip coutputd0clusterTypeSPD03Skip coutputd0clusterTypeSPD11Skip coutputd0clusterTypeSPD12Skip coutputd0clusterTypeSPD13Skip coutputd0PID coutputd0Pt coutputNentries coutputEstimVtx";
   if (!strcmp(task_name, "MUON_QA")) folder = "general1 expert general2";
}


void testOutputSize(const char *filename="QA/QAresults.root")
{
   gSystem->Load("libANALYSIS");
   gSystem->Load("libANALYSISalice");
   gSystem->Load("libCORRFW");
// Add aditional AliRoot libraries
   gSystem->Load("libTender");
   gSystem->Load("libPWG0base");
   gSystem->Load("libPWG0dep");
   gSystem->Load("libPWG0selectors");
   gSystem->Load("libPWGPP");
   gSystem->Load("libPWG2");
   gSystem->Load("libPWG2forward");
   gSystem->Load("libEMCALUtils");
   gSystem->Load("libPWG4PartCorrBase");
   gSystem->Load("libPWG4PartCorrDep");
   gSystem->Load("libPWGHFbase");
   gSystem->Load("libPWGmuon");
   gSystem->Load("libPWGmuondep");
   gSystem->Unlink("syswatch.log");
   AliSysInfo::AddStamp("Start", 0, 0);
 
   TList temp;
   TObjString *os;
   if (!TFile::Open(filename)) return;
   TDirectory *cdir = gDirectory;

   for (Int_t imod=0; imod<NMODS; imod++) {
      SetFolder(module_name[imod]);
      TObjArray *arr = folder.Tokenize(" ");
      TIter next(arr);
      cdir->cd();
      if (!gDirectory->cd(module_name[imod])) continue;
      printf("Module: %s folder: %s\n", module_name[imod], gDirectory->GetName());
      gDirectory->ls();
      while ((os=(TObjString*)next())) {
         TSeqCollection *list = (TSeqCollection*)gDirectory->Get(os->GetString());
         if (list) list->SetOwner();
         temp.Add(list);
      }
      AliSysInfo::AddStamp(module_name[imod], imod+1,1);
      delete arr;
//      temp.Clear();
   }
   TTree *tree = AliSysInfo::MakeTree("syswatch.log");
   tree->SetName("syswatch");
   tree->SetAlias("event", "id0");
   tree->SetAlias("task",  "id1");
   tree->SetAlias("stage", "id2");
   // Already defined aliases
   // tree->SetAlias("deltaT","stampSec-stampOldSec");
   // tree->SetAlias("T","stampSec-first");
   // tree->SetAlias("deltaVM","(pI.fMemVirtual-pIOld.fMemVirtual)");
   // tree->SetAlias("VM","pI.fMemVirtual");
   TCanvas *c = new TCanvas("SysInfo","PWGPP QA train, run #127719 (2.84 M events)" ,10,10,1200,800);
   tree->SetMarkerStyle(kFullSquare);
   tree->SetMarkerColor(kRed);
   tree->SetMarkerSize(1.5);
   tree->Draw("deltaVM:sname","id1==1","", 1234567890, 0);
   TH1* hist = (TH1*)gPad->GetListOfPrimitives()->FindObject("htemp");
   if (hist) {
      hist->SetTitle("dVM[MB] output list");
      hist->GetXaxis()->SetTitle("Module");
      hist->GetYaxis()->SetTitle("deltaVM [MB]");
   }   
   c->SetGridx();
   c->SetGridy();
//   tree->SetMarkerStyle(kOpenSquare);
//   tree->SetMarkerColor(kRed);
//   tree->SetMarkerSize(1);
//   tree->Draw("VM:event","id1==2","SAME", 1234567890, 0);
   temp.Clear();
   delete tree;
/*
   TFileMerger m;
   m.AddFile("QAresults1.root");
   m.AddFile("QAresults2.root");
   m.OutputFile("QAresults.root");
   m.Merge();
*/   
}
