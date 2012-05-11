void HighPtDeDx_lhc10b_Data_ESDs_merge(const char *dir, Int_t stage=0)
{
// Automatically generated merging macro executed in grid subjobs

   TStopwatch timer;
   timer.Start();

// Reset existing include path and add current directory first in the search
   gSystem->SetIncludePath("-I.");
// Load analysis framework libraries
   gSystem->Load("libANALYSIS");
   gSystem->Load("libOADB");
   gSystem->Load("libANALYSISalice");
   gSystem->Load("libCORRFW");

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
   gSystem->AddIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include ");
   gROOT->ProcessLine(".include $ALICE_ROOT/include");
   printf("Include path: %s\n", gSystem->GetIncludePath());

// Add aditional AliRoot libraries

// Analysis source to be compiled at runtime (if any)
   gROOT->ProcessLine(".L DebugClasses.C+g");
   gROOT->ProcessLine(".L AliAnalysisTaskHighPtDeDx.cxx+g");
   gROOT->ProcessLine(".L AliAnalysisTaskHighPtDeDxV0.cxx+g");

// fast xrootd reading enabled
   printf("!!! You requested FastRead option. Using xrootd flags to reduce timeouts. Note that this may skip some files that could be accessed !!!");
   gEnv->SetValue("XNet.ConnectTimeout",50);
   gEnv->SetValue("XNet.RequestTimeout",50);
   gEnv->SetValue("XNet.MaxRedirectCount",2);
   gEnv->SetValue("XNet.ReconnectTimeout",50);
   gEnv->SetValue("XNet.FirstConnectMaxCnt",1);

// Set temporary merging directory to current one
   gSystem->Setenv("TMPDIR", gSystem->pwd());

// Set temporary compilation directory to current one
   gSystem->SetBuildDir(gSystem->pwd(), kTRUE);

// Connect to AliEn
   if (!TGrid::Connect("alien://")) return;
   TString outputDir = dir;
   TString outputFiles = "EventStat_temp.root,HighPtDeDx_Tree.root,HighPtDeDxV0_Tree.root";
   TString mergeExcludes = " ";
   TObjArray *list = outputFiles.Tokenize(",");
   TIter *iter = new TIter(list);
   TObjString *str;
   TString outputFile;
   Bool_t merged = kTRUE;
   while((str=(TObjString*)iter->Next())) {
      outputFile = str->GetString();
      if (outputFile.Contains("*")) continue;
      Int_t index = outputFile.Index("@");
      if (index > 0) outputFile.Remove(index);
      // Skip already merged outputs
      if (!gSystem->AccessPathName(outputFile)) {
         printf("Output file <%s> found. Not merging again.",outputFile.Data());
         continue;
      }
      if (mergeExcludes.Contains(outputFile.Data())) continue;
      merged = AliAnalysisAlien::MergeOutput(outputFile, outputDir, 100, stage);
      if (!merged) {
         printf("ERROR: Cannot merge %s\n", outputFile.Data());
         return;
      }
   }
   // all outputs merged, validate
   ofstream out;
   out.open("outputs_valid", ios::out);
   out.close();
   // read the analysis manager from file
   if (!outputDir.Contains("Stage")) return;
   AliAnalysisManager *mgr = AliAnalysisAlien::LoadAnalysisManager("HighPtDeDx_lhc10b_Data_ESDs.root");
   if (!mgr) return;
   mgr->SetRunFromPath(mgr->GetRunFromAlienPath(dir));
   mgr->SetSkipTerminate(kFALSE);
   mgr->PrintStatus();
   AliLog::SetGlobalLogLevel(AliLog::kError);
   TTree *tree = NULL;
   mgr->StartAnalysis("gridterminate", tree);
}

