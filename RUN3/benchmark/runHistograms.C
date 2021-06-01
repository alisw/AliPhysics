R__ADD_INCLUDE_PATH($ALICE_ROOT)
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <ANALYSIS/macros/train/AddESDHandler.C>
#include <RUN3/benchmark/AddTaskHistogram.C>

TChain* CreateChain(const char *xmlfile, const char *type="ESD");
TChain *CreateLocalChain(const char *txtfile, const char *type, int nfiles);

void runHistograms(Bool_t mc = kFALSE)
{
   const char *anatype = "ESD";

   TChain *chain = CreateLocalChain("wnlocal.txt", anatype, 10);
   if (!chain) return;
   chain->SetNotify(0x0);
   ULong64_t nentries = chain->GetEntries();
   cout << nentries << " entries in the chain." << endl;

   AliAnalysisManager *mgr = new AliAnalysisManager("Histogram example");
   AliESDInputHandler *handler = AddESDHandler();
      
   AddTaskHistogram("");
   
   if (!mgr->InitAnalysis()) return;
   mgr->PrintStatus();

   mgr->SetDebugLevel(1);
   //   mgr->StartAnalysis("localfile", chain, 123456789, 0);
   mgr->StartAnalysis("localfile", chain, nentries, 0);
}

TChain *CreateLocalChain(const char *txtfile, const char *type, int nfiles)
{
   TString treename = type;
   treename.ToLower();
   treename += "Tree";
   printf("***************************************\n");
   printf("    Getting chain of trees %s\n", treename.Data());
   printf("***************************************\n");
   // Open the file
   ifstream in;
   in.open(txtfile);
   Int_t count = 0;
    // Read the input list of files and add them to the chain
   TString line;
   TChain *chain = new TChain(treename);
   while (in.good())
   {
      in >> line;
      if (line.IsNull() || line.BeginsWith("#")) continue;
      if (count++ == nfiles) break;
      TString esdFile(line);
      TFile *file = TFile::Open(esdFile);
      if (file && !file->IsZombie()) {
         chain->Add(esdFile);
         file->Close();
      } else {
         Error("GetChainforTestMode", "Skipping un-openable file: %s", esdFile.Data());
      }   
   }
   in.close();
   if (!chain->GetListOfFiles()->GetEntries()) {
       Error("CreateLocalChain", "No file from %s could be opened", txtfile);
       delete chain;
       return nullptr;
   }
   return chain;
}

//________________________________________________________________________________
TChain* CreateChain(const char *xmlfile, const char *type)
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
   TGridCollection *coll = gGrid->OpenCollection(xmlfile);
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
