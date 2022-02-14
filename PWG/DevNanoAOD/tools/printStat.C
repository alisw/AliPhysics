void printStat(const char* path = "alien:///alice/cern.ch/user/a/alitrain/PWGZZ/NanoAOD_Analysis/35_20190606-0914_child_%d/merge/AnalysisResults.root", int nChild = 10)
{
  TGrid::Connect("alien:");
  
  TString messages;
  
  for (int i=1; i<=nChild; i++)
  {
    TString folder;
    folder.Form(path, i);
    
    TFile* file = TFile::Open(folder);
    if (!file)
      continue;
    
    auto list = (TList*) file->Get("NanoAODNormalisation/Normalisation");
    
    const char* hists[] = { "fCandidateEvents_Skimming", "fSelectedEvents_Skimming", "fCandidateEvents_Filter", "fSelectedEvents_Filter" };
    
    for (int j=0; j<4; j++) {
      auto hist = (TH2*) list->FindObject(hists[j]);
      if (!hist) {
        messages += Form("Not found: %d %s\n", i, hists[j]);
        continue;
      }
      
      auto proj = hist->ProjectionY();
      int entries = (int) proj->GetBinContent((j % 2 == 0) ? 1 : 5);
      
      printf("%12d   ", (int) entries);
    }
    printf("\n");
  }
  
  printf("%s", messages.Data());
}
