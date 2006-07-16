#include "TNtuple.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

//______________________________________________________________________________
class TaskPt : public AliAnalysisTask {

public:
   TaskPt(const char *name);
   virtual ~TaskPt() {;}
   
   virtual void              Exec(Option_t *option);
   
   ClassDef(TaskPt, 0); // example of analysis
};

//______________________________________________________________________________
class TaskCombine : public AliAnalysisTask {

public:
   TaskCombine(const char *name);
   virtual ~TaskCombine() {;}
   
   virtual void              Exec(Option_t *option);
   
   ClassDef(TaskCombine, 0); // example of combination
};

ClassImp(TaskPt)

//______________________________________________________________________________
TaskPt::TaskPt(const char *name)
       :AliAnalysisTask(name,"")
{
// Constructor.
   // Input slot #0 works with an Ntuple
   DefineInput(0, TNtuple::Class());
   // Output slot #0 writes into a TH1 container
   DefineOutput(0, TH1F::Class());
}

//______________________________________________________________________________
void TaskPt::Exec(Option_t *)
{
// Task making a pt distribution.
   static Int_t icalls=0;
   icalls++;
   TCanvas *c1 = new TCanvas(Form("c%i",icalls),"Dynamic Filling Example",200,10,700,500);
   TH1F *h1 = new TH1F(Form("hpt%i",icalls),"This is the pt distribution",100,0,5);
   // Get input data
   TNtuple *ntuple = (TNtuple*)GetInputData(0);
   if (!ntuple) {
      printf("WOOPS ! Where is input 0 for %s ?\n", GetName());
      return;
   }
   Float_t px, py, pt;
   const Int_t kUPDATE = 1000;
   ntuple->SetBranchAddress("px", &px);
   ntuple->SetBranchAddress("py", &py);
   Long64_t nentries = ntuple->GetEntries();
   for (Long64_t i=0; i<nentries; i++) {
      ntuple->GetEntry(i);
      pt = TMath::Sqrt(px*px+py*py);
      h1->Fill(pt);
      if (i && (i%kUPDATE) == 0) {
         if (i == kUPDATE) h1->Draw();
         c1->Modified();
         c1->Update();
         if (gSystem->ProcessEvents())
            break;
      }
   }
   // Post final data. It will be written to a file with option "RECREATE"
//   h1->Draw();
   PostData(0, h1, "RECREATE");
}      

ClassImp(TaskCombine)

//______________________________________________________________________________
TaskCombine::TaskCombine(const char *name)
       :AliAnalysisTask(name,"")
{
// Constructor.
   // Input slot #0 works with a TH1
   DefineInput(0, TH1F::Class());
   // Input slot #1 works with a TH1
   DefineInput(1, TH1F::Class());
   // Output slot #0 writes into a TH1 container
   DefineOutput(0, TH1F::Class());
}
   
//______________________________________________________________________________
void TaskCombine::Exec(Option_t *)
{
// Task combining 2 histograms.
   new TCanvas("c2","h1+h2",200,10,700,500);
   TH1F *h1 = (TH1F*)GetInputData(0);
   TH1F *h2 = (TH1F*)GetInputData(1);
   TH1F *hsum = new TH1F("hsum","h1+h2",100,0,5);
   hsum->Add(h1);
   hsum->Add(h2);
   hsum->Draw();
   PostData(0,hsum);
}   


//______________________________________________________________________________
void testAna()
{
   // Make the analysis manager
   TFile *f = 0;
   AliAnalysisManager *mgr = new AliAnalysisManager();
   // Make a task
   AliAnalysisTask *task1 = new TaskPt("PtTask1");
   mgr->AddTask(task1);
   AliAnalysisTask *task2 = new TaskPt("PtTask2");
   mgr->AddTask(task2);
   AliAnalysisTask *task = new TaskCombine("TaskCombine");
   mgr->AddTask(task);
   // Create containers for input/output
   AliAnalysisDataContainer *cinput = mgr->CreateContainer("cntuple", 
                      TNtuple::Class(),AliAnalysisManager::kInputContainer);
   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist1", TH1::Class());
   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("chist2", TH1::Class());
   AliAnalysisDataContainer *coutput = mgr->CreateContainer("chistsum", 
                      TH1::Class(),AliAnalysisManager::kOutputContainer);
   coutput1->SetFileName("output1.root");                   
   coutput2->SetFileName("output2.root");                   
   mgr->ConnectInput(task1,0,cinput);
   mgr->ConnectInput(task2,0,cinput);
   mgr->ConnectOutput(task1,0,coutput1);
   mgr->ConnectOutput(task2,0,coutput2);
   mgr->ConnectInput(task,0,coutput1);
   mgr->ConnectInput(task,1,coutput2);
   mgr->ConnectOutput(task,0,coutput);
   // Open data
   if (!gSystem->AccessPathName("hsimple.root")) {
      f = new TFile("hsimple.root");
      TNtuple *ntpl = (TNtuple*)f->Get("ntuple");
      cinput->SetData(ntpl);
   } else {
      printf("FIRST run $ROOTSYS/tutorials/hsimple.C\n");
      return;
   }   
   
   if (mgr->InitAnalysis()) {
      mgr->PrintStatus();
      mgr->ExecAnalysis();
   }   
}                         
                      
