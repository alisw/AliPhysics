#ifndef __CINT__
#endif
AliAnalysisManager *AliRsnTrain()
{

  // debug for ANALYSIS manager
  Int_t debugAnalysisMgr = 0;

  // debug level for RESONANCE package
  AliLog::EType_t debugRsnType = AliLog::kInfo;
  Int_t debugRsn = 0;
//   debugRsnType = AliLog::kDebug+debugRsn;

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    mgr = new AliAnalysisManager ( "RSN Train" );

  // setting AnalysisManager debug level
  mgr->SetDebugLevel ( debugAnalysisMgr );

  // if runESDMCFilter.C is commented then runRsnAnalysis*.C has ESD and MC input
  // if it is included then runRsnAnalysis*.C has aod (without MC) which is produced by runESDMCFilter.C
  // do the ESD and MC filter (you can comment it if you want)
//   AddOneTask("runESDMCFilter.C",debugRsnType);

  // do Resonance analysis
  AddOneTask("runRsnAnalysisSE.C",debugRsnType);
  
  // ME is not supported yet 
  //   AddOneTask("runRsnAnalysisME.C",debugRsnType);

  return mgr;
}

void AddOneTask(TString macro,AliLog::EType_t debugRsn)
{
  gROOT->LoadMacro(Form("%s",macro.Data()));
  macro.ReplaceAll(".C","");
  gROOT->ProcessLine(Form("%s((AliLog::EType_t)%d)",macro.Data(),debugRsn));
}
