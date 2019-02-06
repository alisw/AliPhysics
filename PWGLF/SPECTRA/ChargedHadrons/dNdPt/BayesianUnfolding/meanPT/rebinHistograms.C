
Int_t gNMultBins = 1170 +1; //300
Int_t gBinWidth = 3;    //15

Double_t* getMultBinning(Int_t nBinsMult = 100, Double_t binWidth = 2){

  Double_t* multBinEdges = new Double_t [nBinsMult+1];
  multBinEdges[0] = -0.5;
  multBinEdges[1] = 0.5;
  for(Int_t multBin = 1; multBin < nBinsMult; multBin++){
    multBinEdges[multBin+1] = multBinEdges[multBin] + binWidth;
  }
  return multBinEdges;
}

Double_t* gMultBinEdges = getMultBinning(gNMultBins, gBinWidth);



void rebinHistograms(string fileName, string inputPath = "Input/Input_", string outputPath = "Input/"){

  /// Define input and output files
  string inputFileName = inputPath + fileName + ".root";
  string outputFileName = outputPath + "Input_" + fileName + "_rebin.root";
  TObjArray* outputHistos = new TObjArray(1);

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");

  TH1D* multDistMeasured = (TH1D*) inputFile->FindObjectAny("multDistMeasured");
  TH1D* multDistMeasuredRebin = getRebinnedMultDist(multDistMeasured);
  outputHistos->Add(multDistMeasuredRebin);

  TH2D* multPtMeasured = (TH2D*) inputFile->FindObjectAny("multPtMeasured");
  TH2D* multPtMeasuredRebin = getRebinnedMultPt(multPtMeasured);
  outputHistos->Add(multPtMeasuredRebin);

  TH2D* responseMatrixOrig = (TH2D*) inputFile->FindObjectAny("responseMatrixOrig");
  TH2D* responseMatrixOrigRebin = getRebinnedMultMult(responseMatrixOrig);
  outputHistos->Add(responseMatrixOrigRebin);

  TH2D* responseMatrixRebin = getResponseMatrix(responseMatrixOrigRebin);
  responseMatrixRebin->SetName("responseMatrix");
  outputHistos->Add(responseMatrixRebin);

  TH2D* responseMatrixTracksOrig = (TH2D*) inputFile->FindObjectAny("responseMatrixTracksOrig");
  TH2D* responseMatrixTracksOrigRebin = getRebinnedMultMult(responseMatrixTracksOrig);
  outputHistos->Add(responseMatrixTracksOrigRebin);

  TH2D* responseMatrixTracksRebin = getResponseMatrix(responseMatrixTracksOrigRebin);
  responseMatrixTracksRebin->SetName("responseMatrixTracks");
  outputHistos->Add(responseMatrixTracksRebin);

  TH2D* multPtUncorrMC = (TH2D*) inputFile->FindObjectAny("multPtUncorrectedMC");
  TH2D* multPtUncorrMCRebin = getRebinnedMultPt(multPtUncorrMC);
  outputHistos->Add(multPtUncorrMCRebin);

  TH2D* multPtMeasuredMC = (TH2D*) inputFile->FindObjectAny("multPtMeasuredMC");
  TH2D* multPtMeasuredMCRebin = getRebinnedMultPt(multPtMeasuredMC);
  outputHistos->Add(multPtMeasuredMCRebin);

  TH2D* multPtGeneratedMC = (TH2D*) inputFile->FindObjectAny("multPtGeneratedMC");
  TH2D* multPtGeneratedMCRebin = getRebinnedMultPt(multPtGeneratedMC);
  outputHistos->Add(multPtGeneratedMCRebin);


  TH1D* ratio = getRebinRatio(multDistMeasured, multDistMeasuredRebin);
  outputHistos->Add(ratio);

  TFile* outputFile =  TFile::Open(outputFileName.c_str(),"RECREATE");
  outputFile->cd();

  outputHistos->Write();
//  outputFile->Close();
//  delete outputFile;
}


TH1D* getRebinnedMultDist(TH1D* multDist) {

  string name = multDist->GetName();
  multDist->SetName((name + "_old").c_str());
  TH1D* multDistRebin = new TH1D(name.c_str(),multDist->GetTitle(), gNMultBins, gMultBinEdges);
  multDistRebin->Sumw2();

  for(Int_t multBin = 1; multBin <= multDist->GetNbinsX(); multBin++){

    Double_t content = multDist->GetBinContent(multBin);
    if(!content) continue;
    Double_t error = multDist->GetBinError(multBin);

    Double_t multValue = multDist->GetXaxis()->GetBinCenter(multBin);
    Int_t bin = multDistRebin->FindBin(multValue);
    Double_t currentContent = multDistRebin->GetBinContent(bin);
    Double_t currentError = multDistRebin->GetBinError(bin);
    multDistRebin->SetBinContent(bin, content + currentContent);
    multDistRebin->SetBinError(bin, TMath::Sqrt(currentError*currentError + error*error));

  }
  return multDistRebin;
}


TH2D* getRebinnedMultPt(TH2D* multPt) {

  string name = multPt->GetName();
  multPt->SetName((name + "_old").c_str());

  Double_t* ptBinEdges = multPt->GetYaxis()->GetXbins()->GetArray();
  Int_t nBinsPt = multPt->GetNbinsY();
  TH2D* multPtRebin = new TH2D(name.c_str(),multPt->GetTitle(), gNMultBins, gMultBinEdges, nBinsPt, ptBinEdges);
  multPtRebin->Sumw2();

  for(Int_t multBin = 1; multBin <= multPt->GetNbinsX(); multBin++){
    Double_t multValue = multPt->GetXaxis()->GetBinCenter(multBin);
    for(Int_t ptBin = 1; ptBin <= multPt->GetNbinsY(); ptBin++){

      Double_t content = multPt->GetBinContent(multBin, ptBin);
      if(!content) continue;
      Double_t error = multPt->GetBinError(multBin, ptBin);

      Double_t ptValue = multPt->GetYaxis()->GetBinCenter(ptBin);
      Int_t bin = multPtRebin->FindBin(multValue, ptValue);

      Double_t currentContent = multPtRebin->GetBinContent(bin);
      Double_t currentError = multPtRebin->GetBinError(bin);
      multPtRebin->SetBinContent(bin, content + currentContent);
      multPtRebin->SetBinError(bin, TMath::Sqrt(currentError*currentError + error*error));

    }
  }
  return multPtRebin;
}

TH2D* getRebinnedMultMult(TH2D* multMult) {

  string name = multMult->GetName();
  multMult->SetName((name + "_old").c_str());

  TH2D* multMultRebin = new TH2D(name.c_str(), multMult->GetTitle(), gNMultBins, gMultBinEdges, gNMultBins, gMultBinEdges);
  multMultRebin->Sumw2();

  for(Int_t multBin = 1; multBin <= multMult->GetNbinsX(); multBin++){
    Double_t multValue = multMult->GetXaxis()->GetBinCenter(multBin);
    for(Int_t multBin2 = 1; multBin2 <= multMult->GetNbinsY(); multBin2++){

      Double_t content = multMult->GetBinContent(multBin, multBin2);
      if(!content) continue;
      Double_t error = multMult->GetBinError(multBin, multBin2);

      Double_t multValue2 = multMult->GetYaxis()->GetBinCenter(multBin2);
      Int_t bin = multMultRebin->FindBin(multValue, multValue2);

      Double_t currentContent = multMultRebin->GetBinContent(bin);
      Double_t currentError = multMultRebin->GetBinError(bin);
      multMultRebin->SetBinContent(bin, content + currentContent);
      multMultRebin->SetBinError(bin, TMath::Sqrt(currentError*currentError + error*error));

    }
  }
  return multMultRebin;
}

TH3D* getRebinnedMultPtMult(TH3D* multPtMult) {

  string name = multPtMult->GetName();
  multPtMult->SetName((name + "_old").c_str());

  Double_t* ptBinEdges = multPtMult->GetYaxis()->GetXbins()->GetArray();
  Int_t nBinsPt = multPtMult->GetNbinsY();
  TH3D* multPtMultRebin = new TH3D(name.c_str(),multPtMult->GetTitle(), gNMultBins, gMultBinEdges, nBinsPt, ptBinEdges, gNMultBins, gMultBinEdges);
  multPtMultRebin->Sumw2();

  for(Int_t multBin = 1; multBin <= multPtMult->GetNbinsX(); multBin++){
    Double_t multValue = multPtMult->GetXaxis()->GetBinCenter(multBin);
    for(Int_t ptBin = 1; ptBin <= multPtMult->GetNbinsY(); ptBin++){
      Double_t ptValue = multPtMult->GetYaxis()->GetBinCenter(ptBin);
      for(Int_t multBin2 = 1; multBin2 <= multPtMult->GetNbinsZ(); multBin2++){

        Double_t content = multPtMult->GetBinContent(multBin, ptBin, multBin2);
        if(!content) continue;
        Double_t error = multPtMult->GetBinError(multBin, ptBin, multBin2);

        Double_t multValue2 = multPtMult->GetZaxis()->GetBinCenter(multBin2);

        Int_t bin = multPtMultRebin->FindBin(multValue, ptValue, multValue2);

        Double_t currentContent = multPtMultRebin->GetBinContent(bin);
        Double_t currentError = multPtMultRebin->GetBinError(bin);
        multPtMultRebin->SetBinContent(bin, content + currentContent);
        multPtMultRebin->SetBinError(bin, TMath::Sqrt(currentError*currentError + error*error));

      }
    }
  }
  return multPtMultRebin;
}


TH1D* getRebinRatio(TH1D* fineBinnedHist, TH1D* coarseBinnedHist){

  TH1D* original = coarseBinnedHist->Clone("original");
  original->Reset();

  for(Int_t i = 1; i <= original->GetNbinsX(); i++){

    Double_t xValue = original->GetXaxis()->GetBinCenter(i);
    Int_t bin = fineBinnedHist->FindBin(xValue);

    Double_t content = fineBinnedHist->GetBinContent(bin);
    Double_t error = fineBinnedHist->GetBinError(bin);
    original->SetBinContent(i, content);
    original->SetBinError(i, error);
  }

  TH1D* ratio = coarseBinnedHist->Clone("ratio");
  ratio->Divide(original);
  delete original;
  return ratio;
}



/// ---------------------------------------------------------------------------
/// Function to normalize response Matrix
/// ---------------------------------------------------------------------------
TH2D* getResponseMatrix(TH2D* responseMatrixOrig){

  TH2D* responseMatrix = responseMatrixOrig->Clone("responseMatrix");
  responseMatrix->Reset();

  for (Int_t NchBin = 1 ; NchBin <= responseMatrix->GetNbinsY() ; NchBin++){
    Double_t integral = 0;
    for (Int_t NaccBin = 1 ; NaccBin <= responseMatrix->GetNbinsX() ; NaccBin++)
      integral += responseMatrixOrig->GetBinContent(NaccBin, NchBin);

    for (Int_t NaccBin = 1 ; NaccBin <= responseMatrix->GetNbinsX() ; NaccBin++){
      Double_t value = responseMatrixOrig->GetBinContent(NaccBin, NchBin);
      if (integral) responseMatrix->SetBinContent(NaccBin, NchBin, value/integral);
    }
  }
  return responseMatrix;
}
