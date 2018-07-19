void appendHistos(const char* inputFilename, const char* cutSetting, const char* outputFilename)
{
  TFile inputFile(inputFilename ,"read");

  TFile outputFile(outputFilename ,"update");
  outputFile.cd();
  const int nPlots = 8;

  string postFix = cutSetting;
  string histNames[nPlots] = {"momentUnfolded1", "momentUnfoldedRMS", "momentUnfolded2", "momentUnfolded3", "momentReweighted1", "momentReweightedRMS", "momentReweighted2", "momentReweighted3"};


  for( int i = 0; i < nPlots; i++){
  
    TH1D* currHist = (TH1D*)inputFile.FindObjectAny((histNames[i]).c_str());
    string currName = histNames[i] + "_" + cutSetting;
    currHist->SetName((currName).c_str());
    
    currHist->Write();
  
  }
  

  outputFile.Close();
  inputFile.Close();
}