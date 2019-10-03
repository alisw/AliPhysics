
void makeCorrectionFile()
{
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include -I$ALICE_ROOT/include");
  string outputFileName = "CorrectionFile.root";

  string parCompCorrFolder            = "/Users/mkrueger/Documents/EarlyRAAPbPb5TeV/data/ParticleCompositionCorrection/";
  string parCompCorrFileNamePP        = parCompCorrFolder + "pp_Monach13.root";
  string parCompCorrFileNamePP7       = parCompCorrFolder + "pp_Perugia11.root";
  string parCompCorrFileNamePPb       = parCompCorrFolder + "pPb_5TeVDPMJET.root";
  string parCompCorrFileNamePbPb_0005 = parCompCorrFolder + "PbPb_c0005_5TeVHIJING.root";
  string parCompCorrFileNamePbPb_0510 = parCompCorrFolder + "PbPb_c0510_5TeVHIJING.root";
  string parCompCorrFileNamePbPb_1020 = parCompCorrFolder + "PbPb_c1020_5TeVHIJING.root";
  string parCompCorrFileNamePbPb_2040 = parCompCorrFolder + "PbPb_c2040_5TeVHIJING.root";
  string parCompCorrFileNamePbPb_4060 = parCompCorrFolder + "PbPb_c4060_5TeVHIJING.root";
  string parCompCorrFileNamePbPb_6080 = parCompCorrFolder + "PbPb_c6080_5TeVHIJING.root";

  parCompCorrFolder                   = "/Users/mkrueger/Documents/XeXeRAA/data/ParticleCompositionCorrection/";
  string parCompCorrFileNameXeXe_0005 = parCompCorrFolder + "XeXe_c0005_5TeVXeXeHIJING.root";
  string parCompCorrFileNameXeXe_0510 = parCompCorrFolder + "XeXe_c0510_5TeVXeXeHIJING.root";
  string parCompCorrFileNameXeXe_1020 = parCompCorrFolder + "XeXe_c1020_5TeVXeXeHIJING.root";
  string parCompCorrFileNameXeXe_2040 = parCompCorrFolder + "XeXe_c2040_5TeVXeXeHIJING.root";
  string parCompCorrFileNameXeXe_4060 = parCompCorrFolder + "XeXe_c4060_5TeVXeXeHIJING.root";
  string parCompCorrFileNameXeXe_6080 = parCompCorrFolder + "XeXe_c6080_5TeVXeXeHIJING.root";

  string parCompCorrHistName          = "ParComCorrFactor";

  string accCorrFilePPb               = "pPb/pPb_accCorr.root";
  string accCorrHistName              = "acc_corr_pt_etaFull";



  TObjArray* outputHistos = new TObjArray(1);

  outputHistos->Add(getHistogram(parCompCorrFileNamePP, parCompCorrHistName, "partCompCorr_pp"));
  outputHistos->Add(getHistogram(parCompCorrFileNamePP7, parCompCorrHistName, "partCompCorr_pp_Perugia11"));

  outputHistos->Add(getHistogram(parCompCorrFileNamePPb, parCompCorrHistName, "partCompCorr_pPb"));

  outputHistos->Add(getHistogram(parCompCorrFileNamePbPb_0005, parCompCorrHistName, "partCompCorr_PbPb_0005"));
  outputHistos->Add(getHistogram(parCompCorrFileNamePbPb_0510, parCompCorrHistName, "partCompCorr_PbPb_0510"));
  outputHistos->Add(getHistogram(parCompCorrFileNamePbPb_1020, parCompCorrHistName, "partCompCorr_PbPb_1020"));
  outputHistos->Add(getHistogram(parCompCorrFileNamePbPb_2040, parCompCorrHistName, "partCompCorr_PbPb_2040"));
  outputHistos->Add(getHistogram(parCompCorrFileNamePbPb_4060, parCompCorrHistName, "partCompCorr_PbPb_4060"));
  outputHistos->Add(getHistogram(parCompCorrFileNamePbPb_6080, parCompCorrHistName, "partCompCorr_PbPb_6080"));

  outputHistos->Add(getHistogram(parCompCorrFileNameXeXe_0005, parCompCorrHistName, "partCompCorr_XeXe_0005"));
  outputHistos->Add(getHistogram(parCompCorrFileNameXeXe_0510, parCompCorrHistName, "partCompCorr_XeXe_0510"));
  outputHistos->Add(getHistogram(parCompCorrFileNameXeXe_1020, parCompCorrHistName, "partCompCorr_XeXe_1020"));
  outputHistos->Add(getHistogram(parCompCorrFileNameXeXe_2040, parCompCorrHistName, "partCompCorr_XeXe_2040"));
  outputHistos->Add(getHistogram(parCompCorrFileNameXeXe_4060, parCompCorrHistName, "partCompCorr_XeXe_4060"));
  outputHistos->Add(getHistogram(parCompCorrFileNameXeXe_6080, parCompCorrHistName, "partCompCorr_XeXe_6080"));

  outputHistos->Add(getHistogram(accCorrFilePPb, accCorrHistName, "accCorr_pPb"));

  outputHistos->Add(getSecondaryScalingFactor("pp"));
  outputHistos->Add(getSecondaryScalingFactor("pPb"));

  outputHistos->Add(getSecondaryScalingFactor("PbPb_0005"));
  outputHistos->Add(getSecondaryScalingFactor("PbPb_0510"));
  outputHistos->Add(getSecondaryScalingFactor("PbPb_1020"));
  outputHistos->Add(getSecondaryScalingFactor("PbPb_2040"));
  outputHistos->Add(getSecondaryScalingFactor("PbPb_4060"));
  outputHistos->Add(getSecondaryScalingFactor("PbPb_6080"));

  outputHistos->Add(getSecondaryScalingFactor("XeXe_0005"));
  outputHistos->Add(getSecondaryScalingFactor("XeXe_0510"));
  outputHistos->Add(getSecondaryScalingFactor("XeXe_1020"));
  outputHistos->Add(getSecondaryScalingFactor("XeXe_2040"));
  outputHistos->Add(getSecondaryScalingFactor("XeXe_4060"));
  outputHistos->Add(getSecondaryScalingFactor("XeXe_6080"));

  // Write output to file
  TFile* outputFile =  new TFile(outputFileName.c_str(),"RECREATE");
  outputFile->cd();
  outputHistos->Write();
}




///------------------------------------------------------------------------------
///----------------------- Open Files -------------------------------------------
///------------------------------------------------------------------------------
TH1D* getHistogram(string fileName, string histName, string newName){

  TFile* inputFile = TFile::Open(fileName.c_str(),"READ");
  if(!inputFile) cout << "Error: " << fileName << " not found." << endl;

  TH1D* corrHist = (TH1D*)inputFile->FindObjectAny(histName.c_str())->Clone(newName.c_str());
  if(!corrHist) cout << "Error: " << histName << " not found in " << fileName << "." << endl;
//  inputFile->Close();
  return corrHist;
}



TH1D* getSecondaryScalingFactor(string collisionSystem)
{
  string name = "secScalingFactors_" + collisionSystem;
  TH1D* secScalingFactors = NULL;
  Double_t binsPt[6] = {0, 0.15, 0.5, 1.0 , 1.5, 10};
//  Double_t binsPtXe[7] = {0, 0.1, 0.2, 0.5, 1.0, 2.0, 10};
//  if(collisionSystem.find("XeXe") != std::string::npos) secScalingFactors = new TH1D(name.c_str(),name.c_str(), 6, binsPtXe);
//  else 
  secScalingFactors = new TH1D(name.c_str(),name.c_str(), 5, binsPt);

  if(collisionSystem == "pp"){ //TODO 7TeV, 13TeV?
    secScalingFactors->SetBinContent(1,1.28);
    secScalingFactors->SetBinContent(2,1.28);
    secScalingFactors->SetBinContent(3,1.52);
    secScalingFactors->SetBinContent(4,1.56);
    secScalingFactors->SetBinContent(5,1.56);
  }
  if (collisionSystem == "pPb"){
    secScalingFactors->SetBinContent(1,1.37);
    secScalingFactors->SetBinContent(2,1.37);
    secScalingFactors->SetBinContent(3,1.70);
    secScalingFactors->SetBinContent(4,1.71);
    secScalingFactors->SetBinContent(5,1.71);
  }
  if(collisionSystem == "PbPb_0005"){
    secScalingFactors->SetBinContent(1,1.47);
    secScalingFactors->SetBinContent(2,1.47);
    secScalingFactors->SetBinContent(3,1.50);
    secScalingFactors->SetBinContent(4,1.58);
    secScalingFactors->SetBinContent(5,1.58);
  }
  if(collisionSystem == "PbPb_0510"){
    secScalingFactors->SetBinContent(1,1.44);
    secScalingFactors->SetBinContent(2,1.44);
    secScalingFactors->SetBinContent(3,1.49);
    secScalingFactors->SetBinContent(4,1.60);
    secScalingFactors->SetBinContent(5,1.60);
  }
  if(collisionSystem == "PbPb_1020"){
    secScalingFactors->SetBinContent(1,1.42);
    secScalingFactors->SetBinContent(2,1.42);
    secScalingFactors->SetBinContent(3,1.50);
    secScalingFactors->SetBinContent(4,1.61);
    secScalingFactors->SetBinContent(5,1.61);
  }
  if(collisionSystem == "PbPb_2040"){
    secScalingFactors->SetBinContent(1,1.41);
    secScalingFactors->SetBinContent(2,1.41);
    secScalingFactors->SetBinContent(3,1.50);
    secScalingFactors->SetBinContent(4,1.61);
    secScalingFactors->SetBinContent(5,1.61);
  }
  if(collisionSystem == "PbPb_4060"){
    secScalingFactors->SetBinContent(1,1.38);
    secScalingFactors->SetBinContent(2,1.38);
    secScalingFactors->SetBinContent(3,1.51);
    secScalingFactors->SetBinContent(4,1.61);
    secScalingFactors->SetBinContent(5,1.61);
  }
  if(collisionSystem == "PbPb_6080"){
    secScalingFactors->SetBinContent(1,1.33);
    secScalingFactors->SetBinContent(2,1.33);
    secScalingFactors->SetBinContent(3,1.47);
    secScalingFactors->SetBinContent(4,1.55);
    secScalingFactors->SetBinContent(5,1.55);
  }
  if(collisionSystem == "XeXe_0005"){
    secScalingFactors->SetBinContent(1,0.86);
    secScalingFactors->SetBinContent(2,0.86);
    secScalingFactors->SetBinContent(3,1.02);
    secScalingFactors->SetBinContent(5,1.54);
    secScalingFactors->SetBinContent(6,1.54);
  }
  if(collisionSystem == "XeXe_0510"){
    secScalingFactors->SetBinContent(1,0.86);
    secScalingFactors->SetBinContent(2,0.86);
    secScalingFactors->SetBinContent(3,1.02);
    secScalingFactors->SetBinContent(5,1.54);
    secScalingFactors->SetBinContent(6,1.54);
  }
  if(collisionSystem == "XeXe_1020"){
    secScalingFactors->SetBinContent(1,0.86);
    secScalingFactors->SetBinContent(2,0.86);
    secScalingFactors->SetBinContent(3,1.02);
    secScalingFactors->SetBinContent(5,1.54);
    secScalingFactors->SetBinContent(6,1.54);
  }
  if(collisionSystem == "XeXe_2040"){
    secScalingFactors->SetBinContent(1,0.86);
    secScalingFactors->SetBinContent(2,0.86);
    secScalingFactors->SetBinContent(3,1.02);
    secScalingFactors->SetBinContent(5,1.54);
    secScalingFactors->SetBinContent(6,1.54);
  }
  if(collisionSystem == "XeXe_4060"){
    secScalingFactors->SetBinContent(1,0.86);
    secScalingFactors->SetBinContent(2,0.86);
    secScalingFactors->SetBinContent(3,1.02);
    secScalingFactors->SetBinContent(5,1.54);
    secScalingFactors->SetBinContent(6,1.54);
  }
  if(collisionSystem == "XeXe_6080"){
    secScalingFactors->SetBinContent(1,0.86);
    secScalingFactors->SetBinContent(2,0.86);
    secScalingFactors->SetBinContent(3,1.02);
    secScalingFactors->SetBinContent(5,1.54);
    secScalingFactors->SetBinContent(6,1.54);
  }

  Double_t binsPtDefault[49] = {0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0,11.0};
  TH1D* scalingFactorsLinearized = new TH1D((name + "_interpol").c_str(), "scalingFactorsLinearized", 48, binsPtDefault);
  linInterpolation(secScalingFactors, scalingFactorsLinearized);
  delete secScalingFactors;
  return scalingFactorsLinearized;
}


void linInterpolation(TH1D* hToInterpol, TH1D* hInterpol){
  Int_t oldBins = hToInterpol->GetNbinsX();
  Int_t newBins = hInterpol->GetNbinsX();

  //First and last bins extra
  Int_t lowBin=hInterpol->FindBin(hToInterpol->GetBinCenter(1));
  Int_t upBin=hInterpol->FindBin(hToInterpol->GetBinCenter(oldBins));

  Double_t lowContent = hToInterpol->GetBinContent(1);
  Double_t upperContent = hToInterpol->GetBinContent(oldBins);

  for(int i=1; i<lowBin; i++){hInterpol->SetBinContent(i,lowContent);}
  for(int i=upBin; i<=newBins; i++){hInterpol->SetBinContent(i,upperContent);}

  //Linear interpolation between bin centers
  for (int k=1; k<oldBins; k++ ){
    Int_t lowerEge = hInterpol->FindBin(hToInterpol->GetBinCenter(k));
    Int_t upperEge = hInterpol->FindBin(hToInterpol->GetBinCenter(k+1));

    for (int i=lowerEge; i<upperEge; i++ ){

      Double_t pt = hInterpol->GetBinCenter(i);
      Double_t x1 = hToInterpol->GetBinCenter(k);
      Double_t x2 = hToInterpol->GetBinCenter(k+1);
      Double_t y1 = hToInterpol->GetBinContent(k);
      Double_t y2 = hToInterpol->GetBinContent(k+1);

      Double_t interVal = ((y2-y1)/(x2-x1)) * pt + ( y2- ( ((y2-y1)/(x2-x1)) * x2 ) );
      hInterpol->SetBinContent(i,interVal);
    }
  }
}
