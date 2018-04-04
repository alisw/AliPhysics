void applyCorrections(string dataSet = "pp_5TeV_08", string listName = "", string outputFileName = "")
{
//  dataSet   = "pp_5TeV_08";
//  dataSet   = "pp_7TeV_08";
  dataSet   = "pp_8TeV_08";
//  dataSet   = "pPb_5TeV_08";
//  dataSet   = "PbPb_5TeV_08";
//  dataSet   = "XeXe_5TeV_08";
//  dataSet   = "pp_13TeV_08";

  string colSys = dataSet.substr(0, dataSet.find("_")); //"pp";//
  cout << colSys << endl;
  if(listName.empty()) listName = "mkrueger_" + colSys + "_eta_0.80_cutMode_100";
cout << listName << endl;
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include -I$ALICE_ROOT/include");

  Bool_t includeCrosscheckHistos = kFALSE;
  Bool_t includeSimulation = kFALSE;
  Double_t etaRange = 2*0.8; // for pPb -0.76..0.3

  string mcFileName = dataSet + "/MC.root";
  string dataFileName = dataSet + "/Data.root";

  if(outputFileName.empty()) outputFileName = dataSet + "/Input_" + dataSet + ".root";

  string correctionFileName = "CorrectionFiles/CorrectionFile.root";


  string parCompCorrHistName = "ParCompCorrFactor";
  string parCompCorrFileNames[6];

  string generator = "";
  if(dataSet.find("pp_7TeV") != string::npos) generator = "_Perugia11";
cout << generator << endl;
  /// Container to store histograms
  TObjArray* outputHistos = new TObjArray(1);

  ///------------------------------------------------------------------------------
  ///----------------------- Open Files -------------------------------------------
  ///------------------------------------------------------------------------------
  cout << endl;

  TFile* dataFile = TFile::Open(dataFileName.c_str(),"READ");

  string dataListName = listName;
  if(listName.empty()) dataListName = dataFile->GetListOfKeys()->At(0)->GetName();
  TList* dataList = (TList*)dataFile->Get(dataListName.c_str());

  TFile* mcFile = TFile::Open(mcFileName.c_str(),"READ");
  string mcListName = listName;
  if(listName.empty()) mcListName = mcFile->GetListOfKeys()->At(0)->GetName();
  TList* mcList = (TList*)mcFile->Get(mcListName.c_str());

  TFile* correctionFile = TFile::Open(correctionFileName.c_str(),"READ");

  cout << "---------------------------------------------------------------------------" << endl;
  if(dataList) cout << "Input Data: " << dataFileName << " : " << dataListName  << " : " << dataList->GetEntries() << endl;
  if(mcList) cout << "Input MC  : " << mcFileName << "   : " << mcListName << " : " << mcList->GetEntries() << endl;
  cout << "Correction File    : " << correctionFileName << "   : " << correctionFile->GetNkeys() << endl;
  cout << "Output    : " << outputFileName << endl;
  cout << "---------------------------------------------------------------------------" << endl;

  if(includeSimulation){

    TFile* simulationFileNoCR = TFile::Open("pp_5TeV_08/pythiaSimWoCR.root","READ");
    TFile* simulationFileMPI = TFile::Open("pp_5TeV_08/pythiaMPI.root","READ");
    THnD* pythiaNoCR = (THnD*)simulationFileNoCR->FindObjectAny("fHistMCMultPtGenerated");
    pythiaNoCR->SetName("pythiaNoCR");
    THnD* pythiaMPI = (THnD*)simulationFileMPI->FindObjectAny("fHistMCMultPtGenerated");
    pythiaMPI->SetName("pythiaMPI");

    TH1D* simEvents = (TH1D*)simulationFileNoCR->FindObjectAny("h1Event");
    TH1D* simEventsMPI = (TH1D*)simulationFileMPI->FindObjectAny("h1Event");
    Double_t nSimEvents = simEvents->GetBinContent(1);
    Double_t nSimEventsMPI = simEventsMPI->GetBinContent(1);
    Double_t normFactorSim = 1 / (2.0 * TMath::Pi() * etaRange * nSimEvents);
    Double_t normFactorSimMPI = 1 / (2.0 * TMath::Pi() * etaRange * nSimEventsMPI);

    TH2D* multPtSimulatedMC = pythiaNoCR->Projection(1,0);
    multPtSimulatedMC->SetName("multPtSimulatedMC");
    normalizeGeneratedMultPt(multPtSimulatedMC, normFactorSim);

    pythiaMPI->GetAxis(2)->SetRange(2,2);
    TH2D* multPtSimulatedMPI1MC = pythiaMPI->Projection(1,0);
    multPtSimulatedMPI1MC->SetName("multPtSimulatedMPI1MC");
    normalizeGeneratedMultPt(multPtSimulatedMPI1MC, normFactorSimMPI);

    pythiaMPI->GetAxis(2)->SetRange(3,3);
    TH2D* multPtSimulatedMPI2MC = pythiaMPI->Projection(1,0);
    multPtSimulatedMPI2MC->SetName("multPtSimulatedMPI2MC");
    normalizeGeneratedMultPt(multPtSimulatedMPI2MC, normFactorSimMPI);

    pythiaMPI->GetAxis(2)->SetRange(4,4);
    TH2D* multPtSimulatedMPI3MC = pythiaMPI->Projection(1,0);
    multPtSimulatedMPI3MC->SetName("multPtSimulatedMPI3MC");
    normalizeGeneratedMultPt(multPtSimulatedMPI3MC, normFactorSimMPI);

    delete pythiaNoCR;
    delete pythiaMPI;
    delete simEvents;
    delete simEventsMPI;

  }

  ///------------------------------------------------------------------------------
  ///----------------------- Control Histograms -----------------------------------
  ///------------------------------------------------------------------------------

  TH1D* eventCountData = dataList->FindObject("fEventCount");
  eventCountData->SetName("eventCountData");
  TH1D* eventCountMC = mcList->FindObject("fEventCount");
  eventCountMC->SetName("eventCountMC");
  TH1D* trackVsParticle = mcList->FindObject("fHistMCTrackParticle");

  TH2D* ptResolution = (TH2D*)((THnD*) mcList->FindObject("fHistMCPtRes"))->Projection(1,0);
  ptResolution->SetName("ptResolution");
  TH2D* etaResolution = (TH2D*)((THnD*) mcList->FindObject("fHistMCEtaRes"))->Projection(1,0);
  etaResolution->SetName("etaResolution");

  if(includeCrosscheckHistos){
    TH2D* multResolution = (TH2D*)((THnD*) mcList->FindObject("fHistMCMultRes"))->Projection(1,0);
    multResolution->SetName("multResolution");
  }

  ///------------------------------------------------------------------------------
  ///----------------------- MC Histograms ----------------------------------------
  ///------------------------------------------------------------------------------

  // Response Matrix Orig //TODO hier ist keine centrality selection gemacht!!!??
  // TODO nacc einschraenken
  THnD* histResponseMatrixOrig = (THnD*)mcList->FindObject("fHistMCResponseMat");
  TH2D* responseMatrixOrig = (TH2D*) histResponseMatrixOrig->Projection(1,0);
  responseMatrixOrig->SetName("responseMatrixOrig");
  responseMatrixOrig->SetTitle("");
  cleanResponseMatrix(responseMatrixOrig);
  TH2D* responseMatrix = getResponseMatrix(responseMatrixOrig);

  // Response Matrix Tracks Orig
  THnD* histResponseMatrixTracksOrig = (THnD*)mcList->FindObject("fHistMCResponseMatTracks");
  TH2D* responseMatrixTracksOrig = (TH2D*) histResponseMatrixTracksOrig->Projection(1,0);
  responseMatrixTracksOrig->SetName("responseMatrixTracksOrig");
  responseMatrixTracksOrig->SetTitle("");
//  cleanResponseMatrix(responseMatrixTracksOrig);
  TH2D* responseMatrixTracks = getResponseMatrix(responseMatrixTracksOrig);
  responseMatrixTracks->SetName("responseMatrixTracks");






  ///------------------------------------------------------------------------------
  ///----------------------- Apply Corrections ------------------------------------
  ///------------------------------------------------------------------------------
  Int_t minCentBin = 1;
  Int_t maxCentBin = 1;
  if(colSys == "PbPb" || colSys == "XeXe") maxCentBin = 6;

  // restrict centrality
  THnD* histEventOriginal = (THnD*)dataList->FindObject("fHistEvent");
  histEventOriginal->GetAxis(1)->SetRange(minCentBin, maxCentBin);

  THnD* histEventOriginalMC = (THnD*)mcList->FindObject("fHistEvent");
  histEventOriginalMC->GetAxis(1)->SetRange(minCentBin, maxCentBin);

  THnD* histTrack = (THnD*)dataList->FindObject("fHistTrack");
  histTrack->GetAxis(3)->SetRange(minCentBin, maxCentBin);

  THnD* histTrackGeneratedOriginalMC = (THnD*)mcList->FindObject("fHistMCMultPtGenerated");
  histTrackGeneratedOriginalMC->GetAxis(2)->SetRange(minCentBin, maxCentBin);

  THnD* histTrackMC = (THnD*)mcList->FindObject("fHistTrack");
  histTrackMC->GetAxis(3)->SetRange(minCentBin, maxCentBin);

  // do projections
  TH1D* multDistMeasured = (TH1D*) histEventOriginal->Projection(0);
  multDistMeasured->SetName("multDistMeasured");
  multDistMeasured->SetTitle("");
  multDistMeasured->GetXaxis()->SetTitle("#it{N}_{acc}");
  multDistMeasured->GetYaxis()->SetTitle("#it{n}(#it{N}_{acc})");

  TH1D* multDistMeasuredMC = (TH1D*) histEventOriginalMC->Projection(0);
  multDistMeasuredMC->SetName("multDistMeasuredMC");
  multDistMeasuredMC->SetTitle("");
  multDistMeasuredMC->GetXaxis()->SetTitle("#it{N}_{acc}");
  multDistMeasuredMC->GetYaxis()->SetTitle("#it{n}_{MC}(#it{N}_{acc})");

  Double_t nEventsData = multDistMeasured->Integral() - multDistMeasured->GetBinContent(1);
  Double_t normFactorData = 1 / (2.0 * TMath::Pi() * etaRange * nEventsData);
  // TODO event inefficiency: after corrections the pt spectra contain info of genEvents/measEvents times more events
  Double_t nEventsMC = multDistMeasuredMC->Integral() - multDistMeasuredMC->GetBinContent(1);
  Double_t normFactorMC = 1 / (2.0 * TMath::Pi() * etaRange * nEventsMC);

  cout << "---------------------------------------------------------------------------" << endl;
  cout << "Events data:                   " << nEventsData << endl;
  cout << "Events MC:                     " << nEventsMC   << endl;
  cout << "---------------------------------------------------------------------------" << endl;
  cout << endl;

  // Measured MultPt Histogram before correction
  TH2D* multPtUncorrected = (TH2D*) histTrack->Projection(0, 2);
  multPtUncorrected->SetName("multPtUncorrected");
  multPtUncorrected->SetTitle("");

  TH2D* multPtUncorrectedMC = (TH2D*) histTrackMC->Projection(0, 2);
  multPtUncorrectedMC->SetName("multPtUncorrectedMC");
  multPtUncorrectedMC->SetTitle("");

  TH2D* multPtGeneratedMC = (TH2D*) histTrackGeneratedOriginalMC->Projection(1,0);
  multPtGeneratedMC->SetName("multPtGeneratedMC");
  multPtGeneratedMC->SetTitle("");
  multPtGeneratedMC->GetXaxis()->SetTitle("#it{N}_{ch}");
  multPtGeneratedMC->GetYaxis()->SetTitle("#it{p}_{T} (GeV/c)");


  // get histos for corrections
  TH2D* multPtMeasured = multPtUncorrected->Clone();
  multPtMeasured->SetName("multPtMeasured");
  multPtMeasured->Reset();

  TH2D* multPtMeasuredMC = multPtUncorrectedMC->Clone();
  multPtMeasuredMC->SetName("multPtMeasuredMC");
  multPtMeasuredMC->Reset();

  THnD* genPrimTracksOrig = (THnD*) mcList->FindObject("fHistMCGenPrimTrack");
  THnD* recPrimTracksOrig = (THnD*) mcList->FindObject("fHistMCRecPrimTrack");
  THnD* recSecTracksOrig  = (THnD*) mcList->FindObject("fHistMCRecSecTrack");
  THnD* recTracksOrig     = (THnD*) mcList->FindObject("fHistMCRecTrack");


  for(Int_t centBin = minCentBin; centBin <= maxCentBin; centBin++){

    string centSuffix = "";
    if(maxCentBin > 1) {
      centSuffix = getCentSuffix(centBin);
      cout << "Applying corrections for centrality " << histTrackMC->GetAxis(3)->GetBinLowEdge(centBin) << "% - " << histTrackMC->GetAxis(3)->GetBinUpEdge(centBin) << "%" << endl;
    }

    // get mult pt for current centrality
    histTrack->GetAxis(3)->SetRange(centBin, centBin);
    histTrackMC->GetAxis(3)->SetRange(centBin, centBin);
    TH2D* multPtUncorrectedCent = (TH2D*) histTrack->Projection(0, 2);
    multPtUncorrectedCent->SetName("dummy1");
    TH2D* multPtUncorrectedCentMC = (TH2D*) histTrackMC->Projection(0, 2);
    multPtUncorrectedCentMC->SetName("dummy2");

    genPrimTracksOrig->GetAxis(2)->SetRange(centBin, centBin);
    recPrimTracksOrig->GetAxis(2)->SetRange(centBin, centBin);
    recSecTracksOrig->GetAxis(2)->SetRange(centBin, centBin);
    recTracksOrig->GetAxis(2)->SetRange(centBin, centBin);

    // Get histograms for efficiency correction
    TH2D* genPrimTracks = (TH2D*) genPrimTracksOrig->Projection(1,0);
    genPrimTracks->SetName((string("genPrimTracks") + centSuffix).c_str());
    TH2D* recPrimTracks = (TH2D*)recPrimTracksOrig->Projection(1,0);
    recPrimTracks->SetName((string("recPrimTracks") + centSuffix).c_str());
    TH2D* recSecTracks  = (TH2D*)recSecTracksOrig->Projection(1,0);
    recSecTracks->SetName((string("recSecTracks") + centSuffix).c_str());
    TH2D* recTracks     = (TH2D*)recTracksOrig->Projection(1,0);
    recTracks->SetName((string("recTracks") + centSuffix).c_str());

    // Calculate 1d efficiencies (pt)
    TH1D* tempHist = genPrimTracks->ProjectionX();
    TH1D* efficiencyPt =  recPrimTracks->ProjectionX((string("efficiencyPt") + centSuffix).c_str());
    efficiencyPt->Divide(tempHist);
    delete tempHist;

    tempHist = genPrimTracks->ProjectionY();
    TH1D* efficiencyEta =  recPrimTracks->ProjectionY((string("efficiencyEta") + centSuffix).c_str());
    efficiencyEta->Divide(tempHist);
    delete tempHist;

    // Calculate 1d correction factor and secondary contamination (pt)
    tempHist = recPrimTracks->ProjectionX();
    TH1D* primCorrFactorPt =  genPrimTracks->ProjectionX((string("primCorrFactorPt") + centSuffix).c_str());
    primCorrFactorPt->Divide(tempHist);
    delete tempHist;

    tempHist = recTracks->ProjectionX();
    TH1D* secContPt =  recSecTracks->ProjectionX((string("secContPt") + centSuffix).c_str());
    secContPt->Divide(tempHist);
    delete tempHist;

    // Calculate correctionFactorMC = effCorr*(1-secCont)
    tempHist = (TH1D*)primCorrFactorPt->Clone("temp");
    tempHist->Multiply(secContPt);
    TH1D* corrFactorMC = primCorrFactorPt->Clone((string("corrFactorMC") + centSuffix).c_str());
    corrFactorMC->Add(tempHist, -1);
    delete tempHist;
    // ... and compare if it yields the same as genPrimTracksPt/recTracksPt
    tempHist = recTracks->ProjectionX();
    TH1D* corrFactorMCOrig =  genPrimTracks->ProjectionX((string("corrFactorMCOrig") + centSuffix).c_str());
    corrFactorMCOrig->Divide(tempHist);
    delete tempHist;


    // Calculate 2d correction factor and secondary contamination (pt,eta)
    TH2D* primCorrFactor =  genPrimTracks->Clone((string("primCorrFactor") + centSuffix).c_str());
    primCorrFactor->Divide(recPrimTracks);

    TH2D* secCont =  recSecTracks->Clone((string("secCont") + centSuffix).c_str());
    secCont->Divide(recTracks);

    // Calculate correctionFactorMC = effCorr*(1-secCont)
    TH2D* tempHist2d = primCorrFactor->Clone("temp2d");
    tempHist2d->Multiply(secCont);
    TH2D* corrFactorMC2d = primCorrFactor->Clone((string("corrFactorMC2d") + centSuffix).c_str());

    // to crosscheck if it worked;
    TH2D* corrFactorOrig2d =  genPrimTracks->Clone((string("corrFactorOrig2d") + centSuffix).c_str());
    corrFactorOrig2d->Divide(recTracks);


    // Get Correction histograms (pt)
    TH1D* partComCorr = correctionFile->FindObjectAny((string("partCompCorr_") + colSys + generator + centSuffix).c_str());
    cout << "  -> " << partComCorr->GetName() << endl;
    TH1D* secScalingFactors = correctionFile->FindObjectAny((string("secScalingFactors_") + colSys + centSuffix + "_interpol").c_str());
    cout << "  -> " << secScalingFactors->GetName() << endl;
    TH1D* accCorr = NULL;
    if(colSys == "pPb") {
      accCorr = (TH1D*) correctionFile->FindObjectAny("accCorr_pPb");
      cout << "  -> " << accCorr->GetName() << endl;
    }
    cout << endl;
    // apply corrections
    TH2D* multPtMeasuredCent = getCorrectedHist(multPtUncorrectedCent, normFactorData, primCorrFactorPt, secContPt, secScalingFactors, partComCorr);
    multPtMeasured->Add(multPtMeasuredCent);

    TH2D* multPtMeasuredCentMC= getCorrectedHist(multPtUncorrectedCentMC, normFactorMC, primCorrFactorPt, secContPt);
    multPtMeasuredMC->Add(multPtMeasuredCentMC);

    delete multPtMeasuredCent;
    delete multPtMeasuredCentMC;
    delete multPtUncorrectedCentMC;
    delete multPtUncorrectedCent;

    outputHistos->Add(genPrimTracks);
    outputHistos->Add(recPrimTracks);
    outputHistos->Add(recSecTracks);
    outputHistos->Add(recTracks);

    outputHistos->Add(efficiencyEta);
    outputHistos->Add(efficiencyPt);

    outputHistos->Add(primCorrFactorPt);
    outputHistos->Add(secContPt);

    outputHistos->Add(corrFactorMCOrig);
    outputHistos->Add(corrFactorMC);

    outputHistos->Add(partComCorr);
    outputHistos->Add(secScalingFactors);
  }


  normalizeGeneratedMultPt(multPtGeneratedMC, normFactorMC);


  ///------------------------------------------------------------------------------
  ///----------------------- Add Output Histograms --------------------------------
  ///------------------------------------------------------------------------------


  if(includeSimulation){
    outputHistos->Add(multPtSimulatedMC);
    outputHistos->Add(multPtSimulatedMPI1MC);
    outputHistos->Add(multPtSimulatedMPI2MC);
    outputHistos->Add(multPtSimulatedMPI3MC);
  }

  // Control Histos
  outputHistos->Add(eventCountData);
  outputHistos->Add(eventCountMC);
  outputHistos->Add(ptResolution);
  outputHistos->Add(etaResolution);

  // Data
  outputHistos->Add(multDistMeasured);

  // MC
  outputHistos->Add(responseMatrixOrig);
  outputHistos->Add(responseMatrixTracksOrig);
  outputHistos->Add(responseMatrix);
  outputHistos->Add(responseMatrixTracks);

  outputHistos->Add(multDistMeasuredMC);

  outputHistos->Add(multPtUncorrected);
  outputHistos->Add(multPtMeasured);
  outputHistos->Add(multPtUncorrectedMC);
  outputHistos->Add(multPtMeasuredMC);
  outputHistos->Add(multPtGeneratedMC);

//  outputHistos->Add(multPtMultGen);

  // Write output to file
  TFile* outputFile =  new TFile(outputFileName.c_str(),"RECREATE");
  outputFile->cd();
  outputHistos->Write();
  cout << endl << endl;
}

void printHistogramStats(TH2D* hist){

  Int_t allBins = 0;
  Int_t emptyBins = 0;

  for(Int_t x = 1; x <= hist->GetNbinsX(); x++){
    for(Int_t y = 1; y <= hist->GetNbinsY(); y++){
      allBins++;
      if(!hist->GetBinContent(x,y)) emptyBins += 1;
    }
  }
  cout << hist->GetName() << " has " << allBins << " bins, where "<< emptyBins << " are empty (" << 100*emptyBins/allBins << " %)" << endl << endl;
}


TH2D* getCorrectedHist(TH2D* multPtUncorrected, Double_t normFactor, TH1D* primCorrFactor, TH1D* secContam, TH1D* secScalingFactors = NULL, TH1D* partComCorr = NULL, TH1D* accCorrFactors = NULL){

  string originalName = multPtUncorrected->GetName();
  TH2D* multPt = multPtUncorrected->Clone("dummyName");
  multPt->Reset();
  multPt->SetName((originalName + "_Corrected").c_str());

  for(Int_t multBin = 1; multBin <= multPt->GetNbinsX(); multBin++){
    for(Int_t ptBin = 1; ptBin <= multPt->GetNbinsY(); ptBin++){

      Double_t currentContent = multPtUncorrected->GetBinContent(multBin, ptBin);
      Double_t currentError = multPtUncorrected->GetBinError(multBin, ptBin);

      Double_t pt = multPtUncorrected->GetYaxis()->GetBinCenter(ptBin);
      Double_t width_pt = multPtUncorrected->GetYaxis()->GetBinWidth(ptBin);

      Double_t primCorr = primCorrFactor->GetBinContent(ptBin);
      Double_t secConta = secContam->GetBinContent(ptBin);

      Double_t pt = multPt->GetYaxis()->GetBinCenter(ptBin);

      Double_t secScaling = 1;
      Double_t partComp = 1;
      Double_t accCorr = 1;

      if(pt < 0.15 || pt > 10.) continue;
      if(secScalingFactors){
        if((pt > secScalingFactors->GetXaxis()->GetXmin()) && (pt < secScalingFactors->GetXaxis()->GetXmax()))
          secScaling = secScalingFactors->GetBinContent(secScalingFactors->FindBin(pt));
      }

      if(partComCorr){
        if((pt > partComCorr->GetXaxis()->GetXmin()) && (pt < partComCorr->GetXaxis()->GetXmax()))
          partComp = partComCorr->GetBinContent(partComCorr->FindBin(pt));
      }

      if(accCorrFactors){
        if((pt > accCorrFactors->GetXaxis()->GetXmin()) && (pt < accCorrFactors->GetXaxis()->GetXmax()))
          accCorr = accCorrFactors->GetBinContent(accCorrFactors->FindBin(pt));
      }

      if(!partComp || !secScaling || !accCorr) {cout << "SOMETHING IS TERRIBLY WRONG" << endl; partComp = 1;}
      Double_t corrFactor = accCorr * normFactor * primCorr * (1 - secScaling*secConta) / (width_pt * pt * partComp);
//      Double_t corrFactor = primCorr * (1 - secConta);

      if(currentContent){
        // also update entries?
        multPt->SetBinContent(multBin, ptBin, currentContent * corrFactor);
        multPt->SetBinError(multBin, ptBin, currentError * corrFactor);
      }
    }
  }
  return multPt;
}

// 3D DPG pt, nacc, nch
void normalizeGeneratedMultPt(THnD* multPtGeneratedMC, Double_t normFactor){

    for(Int_t ptBin = 1; ptBin <= multPtGeneratedMC->GetAxis(0)->GetNbins(); ptBin++){
      for(Int_t multBin = 1; multBin <= multPtGeneratedMC->GetAxis(1)->GetNbins(); multBin++){
        for(Int_t multBin2 = 1; multBin2 <= multPtGeneratedMC->GetAxis(2)->GetNbins(); multBin2++){

          const Int_t index[3] = {ptBin, multBin, multBin2};
          Double_t currentContent = multPtGeneratedMC->GetBinContent(index);
          Double_t currentError = multPtGeneratedMC->GetBinError(index);

          Double_t pt = multPtGeneratedMC->GetAxis(0)->GetBinCenter(ptBin);
          Double_t width_pt = multPtGeneratedMC->GetAxis(0)->GetBinWidth(ptBin);

          Double_t corrFactor = normFactor / (width_pt * pt);
    //      Double_t corrFactor = 1;

          if(currentContent){
            // also update entries?
            multPtGeneratedMC->SetBinContent(index, currentContent * corrFactor);
            multPtGeneratedMC->SetBinError(index, currentError * corrFactor);
          }
        }
      }
    }
}




void normalizeGeneratedMultPt(TH2D* multPtGeneratedMC, Double_t normFactor){

  for(Int_t multBin = 1; multBin <= multPtGeneratedMC->GetNbinsX(); multBin++){
    for(Int_t ptBin = 1; ptBin <= multPtGeneratedMC->GetNbinsY(); ptBin++){

      Double_t currentContent = multPtGeneratedMC->GetBinContent(multBin, ptBin);
      Double_t currentError = multPtGeneratedMC->GetBinError(multBin, ptBin);

      Double_t pt = multPtGeneratedMC->GetYaxis()->GetBinCenter(ptBin);
      Double_t width_pt = multPtGeneratedMC->GetYaxis()->GetBinWidth(ptBin);

      Double_t corrFactor = normFactor / (width_pt * pt);
//      Double_t corrFactor = 1;

      if(currentContent){
        // also update entries?
        multPtGeneratedMC->SetBinContent(multBin, ptBin, currentContent * corrFactor);
        multPtGeneratedMC->SetBinError(multBin, ptBin, currentError * corrFactor);
      }
    }
  }
}




TH3D* getCorrectedHist(TH3D* multPtEtaUncorrected, Double_t normFactor, TH2D* primCorrFactor, TH2D* secContam){

  string originalName = multPtUncorrected->GetName();
  TH3D* multPtEta = multPtEtaUncorrected->Clone("dummyName");
  multPtEta->Reset();
  multPtEta->SetName((originalName + "_Corrected").c_str());

  for(Int_t multBin = 1; multBin <= multPtEta->GetNbinsX(); multBin++){
    for(Int_t ptBin = 1; ptBin <= multPtEta->GetNbinsY(); ptBin++){

      Double_t pt = multPtEtaUncorrected->GetYaxis()->GetBinCenter(ptBin);
      Double_t width_pt = multPtEtaUncorrected->GetYaxis()->GetBinWidth(ptBin);

      for(Int_t etaBin = 1; etaBin <= multPtEta->GetNbinsZ(); etaBin++){

        Double_t currentContent = multPtEtaUncorrected->GetBinContent(multBin, ptBin, etaBin);
        Double_t currentError = multPtEtaUncorrected->GetBinError(multBin, ptBin, etaBin);


        Double_t primCorr = primCorrFactor->GetBinContent(ptBin, etaBin);
        Double_t secConta = secContam->GetBinContent(ptBin, etaBin);

//        Double_t corrFactor = normFactor * primCorr * (1 - secConta) / (width_pt * pt);
        Double_t corrFactor = primCorr * (1 - secConta);

        if(currentContent){
          // also update entries?
          multPtEta->SetBinContent(multBin, ptBin, etaBin, currentContent * corrFactor);
          multPtEta->SetBinError(multBin, ptBin, etaBin, currentError * corrFactor);
        }
      }
    }
  }
  //projection
  return multPtEta;
}








void printHistInfo(THnD* hist){
  Int_t nDimensions = hist->GetNdimensions();
  cout << "   -> " << hist->GetName() << " ( ";
  for(Int_t i = 0; i < nDimensions; i++){
    cout << hist->GetAxis(i)->GetTitle();
    if(i != nDimensions-1) cout << ", ";
  }
  cout << " )" << endl;
}

TH1D* getEfficiencyCorrection(TH2D* multPtUncorrMC, TH2D* multPtGeneratedMC){

  cout << "WARNING: delete Nacc=0 tracks before dividing!!" << endl;
  TH1D* correctionFactors = multPtGeneratedMC->ProjectionY("correctionFactors");
  TH1D* measured = multPtUncorrMC->ProjectionY("dummyName");
  corr->Divide(corr2);
  delete measured;
  return correctionFactors;
}



// function to apply 1/dpt 1/pt and corrhist
void applyCorrections(TH2D* multPt, TH1D* corrHist=NULL){
  return;
  string originalName = multPt->GetName();

  TH2D* yieldSpectra = multPt->Clone("dummyName");
  yieldSpectra->Reset();

  for(Int_t multBin = 1; multBin <= yieldSpectra->GetNbinsX(); multBin++){
    for(Int_t ptBin = 1; ptBin <= yieldSpectra->GetNbinsY(); ptBin++){
      Double_t currentContent = multPt->GetBinContent(multBin, ptBin);
      Double_t currentError = multPt->GetBinError(multBin, ptBin);
      Double_t pt = yieldSpectra->GetYaxis()->GetBinCenter(ptBin);
      Double_t width_pt = yieldSpectra->GetYaxis()->GetBinWidth(ptBin);

      yieldSpectra->SetBinContent(multBin, ptBin, currentContent / (width_pt * pt));
      yieldSpectra->SetBinError(multBin, ptBin, currentError / (width_pt * pt));
    }
  }

  yieldSpectra->SetName(originalName.c_str());
  return yieldSpectra;
}


void applyCorrections(TH3D* multPt){
  return;
  string originalName = multPt->GetName();

  TH3D* yieldSpectra = multPt->Clone("dummyName");
  yieldSpectra->Reset();

  for(Int_t multBin = 1; multBin <= yieldSpectra->GetNbinsX(); multBin++){
    for(Int_t Nch = 1; Nch <= yieldSpectra->GetNbinsZ(); Nch++){
      for(Int_t ptBin = 1; ptBin <= yieldSpectra->GetNbinsY(); ptBin++){
	Double_t currentContent = multPt->GetBinContent(multBin, ptBin, Nch);
	Double_t currentError = multPt->GetBinError(multBin, ptBin, Nch);
	Double_t pt = yieldSpectra->GetYaxis()->GetBinCenter(ptBin);
	Double_t width_pt = yieldSpectra->GetYaxis()->GetBinWidth(ptBin);

	yieldSpectra->SetBinContent(multBin, ptBin, Nch, currentContent / (width_pt * pt));
	yieldSpectra->SetBinError(multBin, ptBin, Nch, currentError / (width_pt * pt));
      }
    }
  }

  yieldSpectra->SetName(originalName.c_str());
  return yieldSpectra;
}

/*
const Int_t Nbins = 3;
Double_t binsPt[Nbins+1] = {0.15,0.5,1,1.5};
TH1D* hSecScaling = new TH1D("hSecScaling","hSecScaling",Nbins,binsPt);
hSecScaling->SetBinContent(1,1.37);
hSecScaling->SetBinContent(2,1.70);
hSecScaling->SetBinContent(3,1.71);
correction->SetSecRatio(hSecScaling);

//  correction->SetAccaptanceCorrFile("$TOOLS/frameWork/pPb/corrFiles/acceptance_corrections_nominal/acc_corr_pt_etaNew.root");
  correction->SetPtResCorrFile("$TOOLS/frameWork/pPb/corrFiles/ptResolCorr.root");
  //TODO what is difference to Patricks File??
  correction->SetParCompCorrFile("$TOOLS/frameWork/pPb/corrFiles/parCompCorr_etaF.root");

*/




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

void cleanResponseMatrix(TH2D* responseMatrixOrig){
  for(Int_t i = 1; i <=  responseMatrixOrig->GetNbinsX(); i++){
    responseMatrixOrig->SetBinContent(i, 1, 0);
    responseMatrixOrig->SetBinContent(1, i, 0);
    responseMatrixOrig->SetBinError(i, 1, 0);
    responseMatrixOrig->SetBinError(1, i, 0);
  }
}

string getCentSuffix(Int_t centBin){
  string centSuffix = "";
  switch (centBin) {
    case 1: centSuffix = "_0005"; break;
    case 2: centSuffix = "_0510"; break;
    case 3: centSuffix = "_1020"; break;
    case 4: centSuffix = "_2040"; break;
    case 5: centSuffix = "_4060"; break;
    case 6: centSuffix = "_6080"; break;
    default: centSuffix = "??";
  }
  return centSuffix;
}
