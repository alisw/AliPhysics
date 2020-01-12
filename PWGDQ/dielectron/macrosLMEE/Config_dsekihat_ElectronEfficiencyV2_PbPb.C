TString names[] = {
"ResolutionTrackCuts"
//,"cutTPC3sigma",
//,"cutITS3sigma"
};
const Int_t nCut = sizeof(names)/sizeof(names[0]);
Int_t GetN(){return nCut;}

AliAnalysisFilter *Config_dsekihat_ElectronEfficiencyV2_PbPb(
    const Int_t cutID,
    const Bool_t isAOD,
    const Float_t PtMin ,
    const Float_t PtMax ,
    const Float_t EtaMin,
    const Float_t EtaMax
    )
{

  AliAnalysisFilter *anaFilter = new AliAnalysisFilter(Form("anaFilter_%s",names[cutID].Data()),Form("anaFilter_%d",cutID)); // named constructor seems mandatory!
  LMEECutLib *lib = new LMEECutLib(names[cutID],isAOD); 
  anaFilter->AddCuts(lib->SetupTrackCuts(PtMin,PtMax,EtaMin,EtaMax));

  if(names[cutID].Contains("Resolution",TString::kIgnoreCase) || names[cutID].Contains("noPID",TString::kIgnoreCase)){
    printf("Do not add PID cut for ResolutionTrackCuts/noPID\n");
  }
  else anaFilter->AddCuts(lib->SetupPIDCuts());
  if(!isAOD) anaFilter->AddCuts(lib->SetupESDtrackCuts());

  anaFilter->Print();
  return anaFilter;

}
