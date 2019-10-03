void memoryCheck()
{
// include path
   gSystem->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD ");
   gSystem->AddIncludePath("-I$ALICE_ROOT/include");

// Load analysis framework libraries
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
   TFile *f = new TFile("QA/syswatch.root");
   if (!f) return;
   TTree *t = (TTree*)f->Get("syswatch");
// read the analysis manager from file
   TFile *file = TFile::Open("QA/QA.root");
   if (!file) return;
   TIter nextkey(file->GetListOfKeys());
   AliAnalysisManager *mgr = 0;
   TKey *key;
   while ((key=(TKey*)nextkey())) {
      if (!strcmp(key->GetClassName(), "AliAnalysisManager"))
         mgr = (AliAnalysisManager*)file->Get(key->GetName());
   };
   if (!mgr) {
      ::Error("mrmoryCheck", "No analysis manager found in file QA.root");
      return;
   }
   
   TString task_name;
   TIter next(mgr->GetTasks());
   TObject *task;
   while ((task=next())) task_name += Form("%s ", task->GetName());
   t->SetAlias("event", "id0");
   t->SetAlias("RM","pI.fMemResident");
   TCanvas *canvas = new TCanvas("SysInfo",Form("sysinfo QA for %s", task_name.Data()),10,10,1200,1000);
   canvas->Divide(1,2);
   canvas->cd(1)->SetBorderMode(0);
   t->SetMarkerStyle(kCircle);
   t->SetMarkerColor(kRed);
   t->Draw("VM:event","id1==-1 && id2==-1","", 1234567890, 0);
   TH1* hist = (TH1*)gPad->GetListOfPrimitives()->FindObject("htemp");
   if (hist) {
      hist->SetTitle(Form("VM[MB] %s",task_name.Data()));
      hist->GetYaxis()->SetTitle("VM [MB]");
   }   
   canvas->cd(2)->SetBorderMode(0);
   t->SetMarkerStyle(kOpenSquare);
   t->SetMarkerColor(kBlue);
   t->Draw("RM:event","id1==-1 && id2==-1","", 1234567890, 0);
   TH1* hist = (TH1*)gPad->GetListOfPrimitives()->FindObject("htemp");
   if (hist) {
      hist->SetTitle(Form("RM[MB] %s", task_name.Data()));
      hist->GetYaxis()->SetTitle("RM [MB]");
   } 
   canvas->SaveAs("syswatch.gif");
   delete t;
}
