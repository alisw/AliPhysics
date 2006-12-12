changeProductionCrossSections() {
  //
  // calculates dN/dEta from the standard correction and from
  // corrections that has been evaluated with different relative cross
  // sections (of single diff, double diff and non diff). The ratios
  // (standard to changed x-section) of the different dN/deta
  // distributions are saved to a file.
  //

  gSystem->Load("libPWG0base");

  // folder names in systematics root file
  const Char_t* folderNames[] = {"triggerBiasDD","triggerBiasSD","triggerBiasND","vertexRecoDD","vertexRecoSD","vertexRecoND"};

  const Char_t* changes[]  = {"pythia","ddmore","ddless","sdmore","sdless", "dmore", "dless"};

  Float_t scalesDD[] = {1.0, 1.5, 0.5, 1.0, 1.0, 1.5, 0.5};
  Float_t scalesSD[] = {1.0, 1.0, 1.0, 1.5, 0.5, 1.5, 0.5};

  // cross section from Pythia
  //Float_t sigmaND = 55.2;
  //Float_t sigmaDD = 9.78;
  //Float_t sigmaSD = 14.30;

  // getting data
  TFile* finCorr = TFile::Open("analysis_esd.root");
  dNdEtaAnalysis* analysis = new dNdEtaAnalysis("dndeta", "dndeta");
  analysis->LoadHistograms();

  // getting corrections
  AlidNdEtaCorrection* correctionStandard = new AlidNdEtaCorrection("dndeta_correction","dndeta_correction");

  TFile* finCorr = TFile::Open("correction_map.root");
  correctionStandard->LoadHistograms();

  // getting corrections for pure SD, DD and ND
  TFile* finSyst = TFile::Open("systematics.root");
  AlidNdEtaCorrection* correctionsSyst[7];
  for (Int_t i=0; i<6; i++) {
    correctionsSyst[i] = new AlidNdEtaCorrection(folderNames[i],folderNames[i]);
    correctionsSyst[i]->LoadHistograms();
  }

  TH1F* hRatios[21];
  for (Int_t j=0; j<3; j++) { // j = 0 (change vtx), j = 1 (change trg), j = 2 (change both)    

    for (Int_t i=0; i<7; i++) { // changes in cross section
      
      // temporary correction with changed cross sections
      AlidNdEtaCorrection* correctionTmp = (AlidNdEtaCorrection*)correctionStandard->Clone();;

      correctionTmp->Reset();

      if (j==0 || j==2) {
	//correctionTmp->GetVertexRecoCorrection()->Reset();
	//correctionTmp->GetVertexRecoCorrection()->Add(correctionsSyst[3]->GetVertexRecoCorrection(),scalesDD[i]);
	//correctionTmp->GetVertexRecoCorrection()->Add(correctionsSyst[4]->GetVertexRecoCorrection(),scalesSD[i]);
	//correctionTmp->GetVertexRecoCorrection()->Add(correctionsSyst[5]->GetVertexRecoCorrection(),1);

	correctionTmp->Add(correctionsSyst[3], scalesDD[i]);
	correctionTmp->Add(correctionsSyst[4], scalesSD[i]);
	correctionTmp->Add(correctionsSyst[5], 1);	
      }
      if (j==1 || j==2) {
	//correctionTmp->GetTriggerBiasCorrectionINEL()->Reset();
	//correctionTmp->GetTriggerBiasCorrectionINEL()->Add(correctionsSyst[0]->GetTriggerBiasCorrectionINEL(),scalesDD[i]);
	//correctionTmp->GetTriggerBiasCorrectionINEL()->Add(correctionsSyst[1]->GetTriggerBiasCorrectionINEL(),scalesSD[i]);
	//correctionTmp->GetTriggerBiasCorrectionINEL()->Add(correctionsSyst[2]->GetTriggerBiasCorrectionINEL(),1);

	correctionTmp->Add(correctionsSyst[0], scalesDD[i]);
	correctionTmp->Add(correctionsSyst[1], scalesSD[i]);
	correctionTmp->Add(correctionsSyst[2], 1);
      }
      correctionTmp->Finish();

      dNdEtaAnalysis* analysisTmp = (dNdEtaAnalysis*)analysis->Clone();
      analysisTmp->Finish(correctionTmp, 0.3, AlidNdEtaCorrection::kINEL);

      Int_t counter = i + j*7;
      
      hRatios[counter] = (TH1F*)analysisTmp->GetdNdEtaHistogram(2)->Clone();
      
      TString name("ratio");
      if (j==0) name.Append("_vetexReco_");
      if (j==1) name.Append("_triggerBias_");
      if (j==2) name.Append("_vertexReco_triggerBias_");
      name.Append(changes[i]);
      hRatios[counter]->SetName(name.Data());
      name.Append(Form(" (DD #times %0.1f, SD #times %0.1f)",scalesDD[i],scalesSD[i]));
      hRatios[counter]->SetTitle(name.Data());
      hRatios[counter]->SetYTitle("ratio (Pythia x-sections)/(changed x-sections)");
      
      delete analysisTmp;
      delete correctionTmp;
    }
  }

  for (Int_t i=1; i<21; i++) 
    hRatios[i]->Divide(hRatios[0],hRatios[i],1,1);
  hRatios[0]->Divide(hRatios[0],hRatios[0],1,1);

  TFile* fout = new TFile("systematics_xsections.root","RECREATE");
  
  for (Int_t i=0; i<21; i++) 
    hRatios[i]->Write();
  
  fout->Write();
  fout->Close();
}

changeParticleComposition() {
  //
  // calculates dN/dEta from the standard correction and from
  // corrections that has been evaluated with different relative
  // abundancies of kaons and protons and save the ratios in a file

  gSystem->Load("libPWG0base");

  const Char_t* folderNames[] = {"correction_0","correction_1","correction_2","correction_3"};
  Float_t  scalesPi[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  Float_t  scalesK[]  = {0.5, 1.5, 1.0, 1.0, 0.5, 1.5};
  Float_t  scalesP[]  = {1.0, 1.0, 0.5, 1.5, 0.5, 1.5};

  // getting data
  TFile* finCorr = TFile::Open("analysis_esd.root");
  dNdEtaAnalysis* analysis = new dNdEtaAnalysis("dndeta", "dndeta");
  analysis->LoadHistograms();

  // getting corrections for pi, K, p and other particles
  TFile* finSyst = TFile::Open("systematics.root");
  AlidNdEtaCorrection* correctionsSyst[4];
  for (Int_t i=0; i<4; i++) {
    correctionsSyst[i] = new AlidNdEtaCorrection(folderNames[i],folderNames[i]);
    correctionsSyst[i]->LoadHistograms();
  }
  AlidNdEtaCorrection* correctionAll = (AlidNdEtaCorrection*) correctionsSyst[0]->Clone();
  correctionAll->Reset();
  correctionAll->Add(correctionsSyst[0]);
  correctionAll->Add(correctionsSyst[1]);
  correctionAll->Add(correctionsSyst[2]);
  correctionAll->Add(correctionsSyst[3]);
  correctionAll->Finish();
  analysis->Finish(correctionAll, 0.3, AlidNdEtaCorrection::kINEL);

  TH1F* hRatios[6];
  for (Int_t i=0; i<6; i++) { 
    // temporary correction with changed particle composistion
    AlidNdEtaCorrection* correctionTmp = (AlidNdEtaCorrection*)correctionAll->Clone();
    
    correctionTmp->Reset();
    
    correctionTmp->Add(correctionsSyst[0],scalesPi[i]);    
    correctionTmp->Add(correctionsSyst[1],scalesK[i]);    
    correctionTmp->Add(correctionsSyst[2],scalesP[i]);    
    correctionTmp->Add(correctionsSyst[3],1);    
    
    correctionTmp->Finish();

    dNdEtaAnalysis* analysisTmp = (dNdEtaAnalysis*)analysis->Clone();
    analysisTmp->Finish(correctionTmp, 0.3, AlidNdEtaCorrection::kINEL);
    
    hRatios[i] = (TH1F*)analysisTmp->GetdNdEtaHistogram(2)->Clone();
    hRatios[i]->Divide((TH1F*)analysis->GetdNdEtaHistogram(2)->Clone(),hRatios[i],1,1,"B"); 

    TString name(Form("ratio_%d",i));
    hRatios[i]->SetName(name.Data());
    name = TString(Form("#pi #times %0.1f, K #times %0.1f, p #times %0.1f",scalesPi[i],scalesK[i],scalesP[i]));
    hRatios[i]->SetTitle(name.Data());
    hRatios[i]->SetYTitle("ratio (standard/changed compositions)");    
  }

  TFile* fout = new TFile("systematics_composition.root","RECREATE");
  
  for (Int_t i=0; i<6; i++) 
    hRatios[i]->Write();
  
  fout->Write();
  fout->Close();


}
