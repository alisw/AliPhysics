#ifndef __CINT__
#include <Riostream.h>
#include <TGrid.h>
#include <TGridResult.h>
#endif
void RsnGridPlugin(TString analysisMode) {

   Bool_t valid = kTRUE;
   TString dsConfig = AliAnalysisManager::GetGlobalStr("rsnTrainDSConfig",valid);
   Int_t globalTrainID = AliAnalysisManager::GetGlobalInt("rsnGlobalTrainID",valid);

   Int_t numRuns = AliAnalysisManager::GetGlobalInt("rsnGridNumRuns",valid);
   Int_t numRunsSkip = AliAnalysisManager::GetGlobalInt("rsnGridNumRunsSkip",valid);

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Printf("Error[RsnGridPlugin] mgr is null !!!");
      return;
   }

   AliAnalysisAlien *plugin = (AliAnalysisAlien *) mgr->GetGridHandler();
   if (!plugin) {
      Printf("Error[RsnGridPlugin] : plugin is null !!!");
      return;
   }

   // getting latest train id
   TString rsnTrainName = gSystem->BaseName(dsConfig.Data());
   rsnTrainName.ReplaceAll(".txt","");
   rsnTrainName.Append(TString::Format("/%03d/%d_%d",globalTrainID,numRunsSkip,numRuns).Data());

   if (!gGrid) TGrid::Connect("alien://");
   if (!gGrid) return;
   TGridResult *r = gGrid->Query(TString::Format("%s/RsnTrain/%s",gGrid->GetHomeDirectory(),rsnTrainName.Data()).Data(),"*/analysis.root");
   Int_t idRsnTrain = 0;
   if (r) {
      TString s = r->GetKey(r->GetSize()-1,"lfn");
      s.ReplaceAll("/analysis.root","");
      s = gSystem->BaseName(s);
      if (!s.IsNull()) idRsnTrain = s.Atoi();
      if (!analysisMode.CompareTo("full")) idRsnTrain++;
   }
   rsnTrainName.Append(Form("/%03d",idRsnTrain));

   TString rsnTrainWkDir = TString::Format("RsnTrain/%s",rsnTrainName.Data()).Data();
   Info("RsnGridPlugin()",TString::Format("RSN Train directory : %s%s",gGrid->GetHomeDirectory(),rsnTrainWkDir.Data()).Data());

   plugin->SetGridWorkingDir(rsnTrainWkDir.Data());
   plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output

   plugin->SetAPIVersion("V1.1x");

   TString rootver = AliAnalysisManager::GetGlobalStr("rsnLegoTrainROOTversion",valid);
   plugin->SetROOTVersion(rootver.Data());

   TString alirootVersion = AliAnalysisManager::GetGlobalStr("rsnLegoTrainAliROOTversion",valid);
   if (alirootVersion.IsNull()) {
      if (gSystem->Getenv("ALICE_ROOT")) alirootVersion = gSystem->GetFromPipe("aliroot --version | awk '{print $3}'");
   }
   plugin->SetAliROOTVersion(alirootVersion.Data());

   plugin->SetExecutableCommand("aliroot -b -q");

   plugin->SetAnalysisMacro("RsnTrain.C");
   plugin->SetMasterResubmitThreshold(90);
   plugin->SetTTL(84600);
   plugin->SetInputFormat("xml-single");
   plugin->SetJDLName("RsnTrain.jdl");
   plugin->SetPrice(1);
   plugin->SetSplitMode("se");
   plugin->SetNtestFiles(2);
   plugin->SetMergeViaJDL();
   plugin->SetOverwriteMode(kFALSE);
   //    plugin->SetKeepLogs(kTRUE);

   RsnSetData(plugin,dsConfig,numRuns,numRunsSkip,1000);

   plugin->SetSplitMaxInputFileNumber(25);

   //   Fatal("RsnDataSet","No dataset found !!!");
}

void RsnSetData(AliAnalysisAlien *plugin,TString dsConf,Int_t numRuns = 1000,Int_t numRunsSkip=0,Int_t maxRunsPerMaster = 1000) {

   Bool_t dsFound = kTRUE;
   Int_t nRunsPerMaster = 0;

   Bool_t valid = kTRUE;
   TString legoTrainPath = AliAnalysisManager::GetGlobalStr("rsnLegoTrainPath",valid);

   if (gSystem->AccessPathName(dsConf.Data())) dsConf.Prepend(Form("%s/",legoTrainPath.Data()));
   dsConf = gSystem->ExpandPathName(dsConf.Data());

   if (dsConf.Contains(".txt")) {
      ifstream in;
      in.open(dsConf.Data());
      if (!in.is_open()) Fatal("RsnSetData",Form("File %s was not found !!!",dsConf.Data()));
      Printf("DS config file : %s",dsConf.Data());
      TString line;
      Bool_t isRun = kFALSE;
      while (in.good())
      {
         in >> line;
         if (line.IsNull()) continue;
         if (line.Contains("BASE")) {
            GetParameterFromConfig(line);
            plugin->SetGridDataDir(line.Data());
            Printf("BASE -> %s",line.Data());
            continue;
         }
         if (line.Contains("PREFIX")) {
            GetParameterFromConfig(line);
            plugin->SetRunPrefix(line.Data());
            Printf("PREFIX -> %s",line.Data());
            continue;
         }
         if (line.Contains("DATA_PATTERN")) {
            GetParameterFromConfig(line);
            plugin->SetDataPattern(line.Data());
            Printf("DATA_PATTERN -> %s",line.Data());
            continue;
         }
         if (!line.CompareTo("RUNS")) {
            isRun = kTRUE;
            in >> line;
         }
         if (isRun) {
            if (numRunsSkip>0) {
               numRunsSkip--;
               continue;
            } else {
               if (nRunsPerMaster < numRuns ) {
                  Printf("Adding RUN : %s",line.Data());
                  plugin->AddRunNumber(line.Data());
                  nRunsPerMaster++;
               } else {
                  break;
               }

            }
         }
      }
   } else {
      plugin->SetGridDataDir("");
      plugin->SetDataPattern("");
      Fatal("RsnDataSet","No dataset found !!!");
   }

   if (nRunsPerMaster > maxRunsPerMaster) nRunsPerMaster = maxRunsPerMaster;
   plugin->SetNrunsPerMaster(nRunsPerMaster);

}

void GetParameterFromConfig(TString &str,TString token="=") {
   TObjArray *array = str.Tokenize(token.Data());
   TObjString *strObj = (TObjString *)array->At(1);
   if (strObj) str = strObj->GetString();
   else str = "";
}

