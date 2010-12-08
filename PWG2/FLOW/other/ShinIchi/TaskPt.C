const char *anatype = "ESD";

void TaskPt()
{
// Analysis using ESD data
// Automatically generated analysis steering macro executed in grid subjobs

   TStopwatch timer;
   timer.Start();

// Set temporary merging directory to current one
   gSystem->Setenv("TMPDIR", gSystem->pwd());

// include path
   gSystem->AddIncludePath("-I$ALICE_ROOT/include");

// Load analysis framework libraries
   gSystem->Load("libANALYSIS");
   gSystem->Load("libANALYSISalice");
   gSystem->Load("libCORRFW");

// Add aditional AliRoot libraries

// analysis source to be compiled at runtime (if any)
   gROOT->ProcessLine(".L AliAnalysisTaskPt.cxx+g");

// connect to AliEn and make the chain
   if (!TGrid::Connect("alien://")) return;
// read the analysis manager from file
   TFile *file = TFile::Open("TaskPt.root");
   if (!file) return;
   TIter nextkey(file->GetListOfKeys());
   AliAnalysisManager *mgr = 0;
   TKey *key;
   while ((key=(TKey*)nextkey())) {
      if (!strcmp(key->GetClassName(), "AliAnalysisManager"))
         mgr = (AliAnalysisManager*)file->Get(key->GetName());
   };
   if (!mgr) {
      ::Error("TaskPt", "No analysis manager found in file TaskPt.root");
      return;
   }

   mgr->PrintStatus();
   AliLog::SetGlobalLogLevel(AliLog::kError);
   TChain *chain = CreateChain("wn.xml", anatype);

   mgr->StartAnalysis("localfile", chain);
   timer.Stop();
   timer.Print();
}

//________________________________________________________________________________
TChain* CreateChain(const char *xmlfile, const char *type="ESD")
{
// Create a chain using url's from xml file
   TString filename;
   Int_t run = 0;
   TString treename = type;
   treename.ToLower();
   treename += "Tree";
   printf("***************************************\n");
   printf("    Getting chain of trees %s\n", treename.Data());
   printf("***************************************\n");
   TAlienCollection *coll = TAlienCollection::Open(xmlfile);
   if (!coll) {
      ::Error("CreateChain", "Cannot create an AliEn collection from %s", xmlfile);
      return NULL;
   }
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   TChain *chain = new TChain(treename);
   coll->Reset();
   while (coll->Next()) {
      filename = coll->GetTURL();
      if (mgr) {
         Int_t nrun = AliAnalysisManager::GetRunFromAlienPath(filename);
         if (nrun && nrun != run) {
            printf("### Run number detected from chain: %d\n", nrun);
            mgr->SetRunFromPath(nrun);
            run = nrun;
         }
      }
      chain->Add(filename);
   }
   if (!chain->GetNtrees()) {
      ::Error("CreateChain", "No tree found from collection %s", xmlfile);
      return NULL;
   }
   return chain;
}

