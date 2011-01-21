const char *anatype = "ESD";

void AnalysisStrange()
{
// Analysis using ESD data
// Automatically generated analysis steering macro executed in grid subjobs

   TStopwatch timer;
   timer.Start();

// Set temporary merging directory to current one
   gSystem->Setenv("TMPDIR", gSystem->pwd());

// Reset existing include path and add current directory first in the search
   gSystem->SetIncludePath("-I.");
// load base root libraries
   gSystem->Load("libTree");
   gSystem->Load("libGeom");
   gSystem->Load("libVMC");
   gSystem->Load("libPhysics");

   gSystem->Load("libMinuit");

// Load analysis framework libraries
   if (!SetupPar("STEERBase")) return;
   if (!SetupPar("ESD")) return;
   if (!SetupPar("AOD")) return;
   if (!SetupPar("ANALYSIS")) return;
   if (!SetupPar("ANALYSISalice")) return;
   if (!SetupPar("CORRFW")) return;

// Compile other par packages
   if (!SetupPar("OADB.par")) return;
// include path
   TString intPath = gInterpreter->GetIncludePath();
   TObjArray *listpaths = intPath.Tokenize(" ");
   TIter nextpath(listpaths);
   TObjString *pname;
   while ((pname=(TObjString*)nextpath())) {
      TString current = pname->GetName();
      if (current.Contains("AliRoot") || current.Contains("ALICE_ROOT")) continue;
      gSystem->AddIncludePath(current);
   }
   if (listpaths) delete listpaths;
   gROOT->ProcessLine(".include $ALICE_ROOT/include");
   printf("Include path: %s\n", gSystem->GetIncludePath());

// Add aditional AliRoot libraries

// analysis source to be compiled at runtime (if any)
   gROOT->ProcessLine(".L AliAnalysisTaskPerformanceStrange.cxx+g");

// connect to AliEn and make the chain
   if (!TGrid::Connect("alien://")) return;
// read the analysis manager from file
   TFile *file = TFile::Open("analysis.root");
   if (!file) return;
   TIter nextkey(file->GetListOfKeys());
   AliAnalysisManager *mgr = 0;
   TKey *key;
   while ((key=(TKey*)nextkey())) {
      if (!strcmp(key->GetClassName(), "AliAnalysisManager"))
         mgr = (AliAnalysisManager*)file->Get(key->GetName());
   };
   if (!mgr) {
      ::Error("AnalysisStrange", "No analysis manager found in file analysis.root");
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

//________________________________________________________________________________
Bool_t SetupPar(const char *package) {
// Compile the package and set it up.
   TString pkgdir = package;
   pkgdir.ReplaceAll(".par","");
   gSystem->Exec(TString::Format("tar xvzf %s.par", pkgdir.Data()));
   TString cdir = gSystem->WorkingDirectory();
   gSystem->ChangeDirectory(pkgdir);
   // Check for BUILD.sh and execute
   if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
      printf("*******************************\n");
      printf("*** Building PAR archive    ***\n");
      printf("*******************************\n");
      if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
         ::Error("SetupPar", "Cannot build par archive %s", pkgdir.Data());
         gSystem->ChangeDirectory(cdir);
         return kFALSE;
      }
   } else {
      ::Error("SetupPar","Cannot access PROOF-INF/BUILD.sh for package %s", pkgdir.Data());
      gSystem->ChangeDirectory(cdir);
      return kFALSE;
   }
   // Check for SETUP.C and execute
   if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
      printf("*******************************\n");
      printf("***    Setup PAR archive    ***\n");
      printf("*******************************\n");
      gROOT->Macro("PROOF-INF/SETUP.C");
   } else {
      ::Error("SetupPar","Cannot access PROOF-INF/SETUP.C for package %s", pkgdir.Data());
      gSystem->ChangeDirectory(cdir);
      return kFALSE;
   }
   // Restore original workdir
   gSystem->ChangeDirectory(cdir);
   return kTRUE;
}
