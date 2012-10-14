AliAnalysisTask *AddRsnTaskTrain(const char *commonStr,const char *rsnStr,const char *rsnCutStr) {
	// rsnStr -> <Name>
	// rsnCutStr -> <CutName>
   // This will use AddRsnPairs<Name>.C
   // and for cuts AddRsnDaughterCuts<CutName>.C
   // and <opt> string is passed to AddRsnDaughterCuts<CutName>.C
   // so you can control different cut settings
   // string "<Name>:mon" means that it will add monitoring histograms from cuts
   // Note : for now you have to set gRsnUseMiniPackage = 0 to have mon histograms
	//    return AddRsnTask("<Name>:mon","<CutName>:<opt>","");
	// or like we are using it now
	//    return AddRsnTask(rsnStr,rsnCutStr,"");


	// Creating Rsn Train Manager
   AliRsnTrainManager *rsnMgr = new AliRsnTrainManager();


   gROOT->LoadMacro("RsnTrainCommonSettings.C");
   RsnTrainCommonSettings(commonStr);

   rsnMgr->Print();

   gROOT->LoadMacro("AddRsnTask.C");
   return AddRsnTask(rsnStr,rsnCutStr,"");
}
