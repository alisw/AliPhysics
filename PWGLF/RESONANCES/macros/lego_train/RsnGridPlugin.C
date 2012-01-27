#ifndef __CINT__
#include <Riostream.h>
#endif
void RsnGridPlugin() {

   Bool_t valid = kTRUE;
   Int_t idRsnTrain = AliAnalysisManager::GetGlobalInt("rsnTrainID",valid);
   TString dsConfig = AliAnalysisManager::GetGlobalStr("rsnTrainDSConfig",valid);

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) { Printf("Error[RsnGridPlugin] mgr is null !!!"); return; }

   AliAnalysisAlien *plugin = (AliAnalysisAlien *) mgr->GetGridHandler();
   if (!plugin) { Printf("Error[RsnGridPlugin] : plugin is null !!!"); return; }

   TString rsnTrainName = gSystem->BaseName(dsConfig.Data());
   rsnTrainName.ReplaceAll(".txt",Form("/%03d",idRsnTrain));

   plugin->SetGridWorkingDir(Form("RsnTrain/%s",rsnTrainName.Data()));
   plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output

   plugin->SetAPIVersion("V1.1x");
   plugin->SetROOTVersion("v5-30-03-1");

   TString alirootVersion="v5-02-16-AN";
   alirootVersion = gSystem->GetFromPipe("aliroot --version | awk '{print $3}'");

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
   //    plugin->SetKeepLogs(kTRUE);

   RsnSetData(plugin,dsConfig,1000);

   plugin->SetSplitMaxInputFileNumber(50);

   //   Fatal("RsnDataSet","No dataset found !!!");
}

void RsnSetData(AliAnalysisAlien *plugin,TString dsConf,Int_t maxRunsPerMaster = 1000) {

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
            //            Printf("Adding RUN : %s",line.Data());
            plugin->AddRunNumber(line.Data());
            nRunsPerMaster++;
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

void GetParameterFromConfig(TString &str,TString token="="){
   TObjArray *array = str.Tokenize(token.Data());
   TObjString *strObj = (TObjString *)array->At(1);
   if (strObj) str = strObj->GetString();
   else str = "";
}
