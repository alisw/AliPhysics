void produceSys(const char* inputFilename, const char* outputFilename)
{

  const int nPlots = 8;
  string histNames[nPlots] = {"momentUnfolded1", "momentUnfoldedRMS", "momentUnfolded2", "momentUnfolded3", "momentReweighted1", "momentReweightedRMS", "momentReweighted2", "momentReweighted3"};

  TFile inputFile(inputFilename ,"READ");

  TObjArray plots(1);

  char currHistName[50];

  for( int i = 0; i < nPlots; i++){

    TList* cutVariations = new TList();
    cutVariations->SetOwner();

    // first put all histograms in array
    for(int cutSetting = 100; cutSetting < 120; cutSetting++){     
      sprintf(currHistName, "%s_%d", histNames[i].c_str(), cutSetting);
      TH1D* currHist = (TH1D*)inputFile.FindObjectAny(currHistName);
      if(!currHist) cout << histNames[i] + " not found." << endl;
      cutVariations->Add(currHist);
    }

    TList* systContrib = new TList();
    systContrib->SetName((histNames[i]).c_str());

    systContrib->Add(makeSystematics(cutVariations, 1, 2, getCutSettingName(101)));
    systContrib->Add(makeSystematics(cutVariations, 3, 4, getCutSettingName(103)));
    systContrib->Add(makeSystematics(cutVariations, 5, 6, getCutSettingName(105)));
    systContrib->Add(makeSystematics(cutVariations, 7, 8, getCutSettingName(107)));
    systContrib->Add(makeSystematics(cutVariations, 9, 10, getCutSettingName(109)));
    systContrib->Add(makeSystematics(cutVariations, 11, 12, getCutSettingName(111)));
    systContrib->Add(makeSystematics(cutVariations, 13, 14, getCutSettingName(113)));
    systContrib->Add(makeSystematics(cutVariations, 15, getCutSettingName(115)));
    systContrib->Add(makeSystematics(cutVariations, 16, 17, getCutSettingName(116)));
    systContrib->Add(makeSystematics(cutVariations, 18, 19, getCutSettingName(118)));

    systContrib->Add(totalSystematics(systContrib));
    systContrib->Add(getHistWithSysErrors(cutVariations, systContrib));
    plots.Add(systContrib);

    delete cutVariations;
  }


  TFile* outputFile =  new TFile(outputFilename, "RECREATE");
  outputFile->cd();

  // Write Histos in file
  for(Int_t i = 0; i < plots.GetEntries(); i++){

    TList* currList = plots.At(i);
    string currListName = currList->GetName();
    currList->Write(currListName.c_str(), TObject::kSingleKey);
  }

  outputFile->Close();
  inputFile.Close();

}


TH1D* totalSystematics(TList* systContrib){

  TH1D* syst = systContrib->At(0)->Clone("total");
  syst->Reset();

  for(Int_t i = 0; i < systContrib->GetEntries(); i++){

    TH1D* currContrib = systContrib->At(i);

    for(Int_t j = 1; j < currContrib->GetXaxis()->GetNbins(); j++){

      Double_t currentValue = syst->GetBinContent(j);
      Double_t newValue = currContrib->GetBinContent(j);

      syst->SetBinContent(j, currentValue + newValue*newValue);
    }

  }

  for(Int_t j = 1; j < syst->GetXaxis()->GetNbins(); j++){

    Double_t currentValue = syst->GetBinContent(j);
    syst->SetBinContent(j, TMath::Sqrt(currentValue));
  }

  return syst;

}


TH1D* makeSystematics(TList* cutVariations, Int_t upID, Int_t downID, string cutSettingName){

  TH1D* nominal = (TH1D*) cutVariations->At(0);
  TH1D* syst = nominal->Clone(cutSettingName.c_str());
  syst->Reset();

  TH1D* up = (TH1D*) (cutVariations->At(upID)->Clone());
  TH1D* down = (TH1D*) (cutVariations->At(downID)->Clone());

  up->Divide(nominal);
  for(Int_t i = 1; i <= up->GetXaxis()->GetNbins(); i++){
    Double_t value = up->GetBinContent(i);
    if(value) value = TMath::Abs(value - 1);
    up->SetBinContent(i, value);
  }
  down->Divide(nominal);
  for(Int_t i = 1; i <= down->GetXaxis()->GetNbins(); i++){
    Double_t value = down->GetBinContent(i);
    if(value) value = TMath::Abs(value - 1);
    down->SetBinContent(i, value);
  }

  for(Int_t i = 1; i <= syst->GetXaxis()->GetNbins(); i++){

    Double_t upValue = up->GetBinContent(i);
    Double_t downValue = down->GetBinContent(i);
    if(upValue > downValue) syst->SetBinContent(i, upValue);
    else syst->SetBinContent(i, downValue);
  }

  syst->GetYaxis()->SetTitle("rel. sys. unc.");
  delete up;
  delete down;
  return syst;
}

TH1D* makeSystematics(TList* cutVariations, Int_t ID, string cutSettingName){

  TH1D* nominal = (TH1D*) cutVariations->At(0);
  TH1D* syst = (TH1D*) (cutVariations->At(ID)->Clone(cutSettingName.c_str()));

  syst->Divide(nominal);
  for(Int_t i = 1; i < syst->GetXaxis()->GetNbins(); i++){
    Double_t value = syst->GetBinContent(i);
    if(value) value = TMath::Abs(value - 1);
    syst->SetBinContent(i, value);
    syst->SetBinError(i, 0);
  }

  syst->GetYaxis()->SetTitle("rel. sys. unc.");
  return syst;
}

TH1D* getHistWithSysErrors(TList* cutVariations, TList* systContrib){

  TH1D* hist = (TH1D*) (cutVariations->At(0))->Clone();
  TH1D* sysErr = (TH1D*) systContrib->Last();

  for(Int_t i = 1; i < hist->GetXaxis()->GetNbins(); i++){
    Double_t err = sysErr->GetBinContent(i);
    hist->SetBinError(i, err);
  }

  return hist;
}


string getCutSettingName(Int_t cutSetting){

  string name = "";
  if(cutSetting == 101 || cutSetting == 102) name = "Maxchi2perITSclu";
  if(cutSetting == 103 || cutSetting == 104) name = "Maxchi2perTPCclu";
  if(cutSetting == 105 || cutSetting == 106) name = "RatioCrossedRowsOverFindableClustersTPC";
  if(cutSetting == 107 || cutSetting == 108) name = "FractionSharedClustersTPC";
  if(cutSetting == 109 || cutSetting == 110) name = "MaxChi2TPCConstrained";
  if(cutSetting == 111 || cutSetting == 112) name = "DCAtoVertexXYPtDep";
  if(cutSetting == 113 || cutSetting == 114) name = "DCAtoVertexZ";
  if(cutSetting == 115) name = "ClusterReqITS";
  if(cutSetting == 116 || cutSetting == 117) name = "Ncrnclgeomlength";
  if(cutSetting == 118 || cutSetting == 119) name = "DeadzoneWidth";

  return name;
}
