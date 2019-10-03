// AliFemtoSplittingMergingQA.C - macro to create the splitting/merging
// test with the pion correlation function.
// As a default the anti-splitting and anti-merging cuts are open
// and the two correlation functions:
// AliFemtoShareQualityCorrFctn and AliFemtoTPCInnerCorrFctn 
// can be used to study the splitting (former) and merging (latter) effect
// If ones needs to produce a "clean" sample with both effects removed, 
// one needs to change the cut values to the "reasonable" ones, or perform
// the full systematic analysis with the above-mentioned functions 

// Author: Adam Kisiel. Adam.Kisiel@cern.ch

// parameters:
// listFileName - a text file containing a list of ESDs (with full paths) 
// to analyze
void AliFemtoSplittingMergingQA(const char *listFileName)
{
  double PionMass = 0.13956995;
  double chargePi = 1.0;
  
  // Load the neccessary libraries
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libEG");
  gSystem->Load("libCint");
  gSystem->Load("libESD");
  gSystem->Load("libPWG2femtoscopy");
  gSystem->Load("libPWG2femtoscopyUser");

  // Set-up the reader for ALICE ESD
  AliFemtoEventReaderESD* Reader=new AliFemtoEventReaderESD();
  // Read only constrained momenta - primordial particles
  Reader->SetConstrained(true);
  Reader->SetReadTPCInner(true);
  // Read the file list from the filename supplied by the user
  Reader->SetInputFile(listFileName);
  
  // Setup the manager 
  AliFemtoManager* Manager=new AliFemtoManager();
  // Point to the data source - the reader
  Manager->SetEventReader(Reader);
  
  // Setup the analysis
  AliFemtoSimpleAnalysis* an =new AliFemtoSimpleAnalysis();
  // Number of events to construct the background
  an->SetNumEventsToMix(3);

  // The event selector
  AliFemtoBasicEventCut* mec = new AliFemtoBasicEventCut();
  // Accept events with the given multiplicity
  mec->SetEventMult(0,100000);
  // and z-vertex distance to the center of the TPC
  mec->SetVertZPos(-1000,1000);
	
  // The track selector
  AliFemtoESDTrackCut* dtc = new AliFemtoESDTrackCut();
  // We want positive pions
  dtc->SetPidProbPion(0.2,1.001);
  dtc->SetPidProbMuon(0.0,0.8);
  dtc->SetPidProbKaon(0.0,0.1);
  dtc->SetPidProbProton(0.0,0.1);
  dtc->SetMostProbablePion();
  dtc->SetCharge(1.0);
  // so we set the correct mass
  dtc->SetMass(PionMass);
  // we select low pt
  dtc->SetPt(0.1,0.7);
  dtc->SetStatus(AliESDtrack::kTPCrefit|AliESDtrack::kITSrefit);
  dtc->SetminTPCncls(95);
  dtc->SetRemoveKinks(kTRUE);
  dtc->SetLabel(kFALSE);
  dtc->SetMaxITSChiNdof(3.0);
  dtc->SetMaxTPCChiNdof(2.0);
  dtc->SetMaxSigmaToVertex(3.0);

  AliFemtoCutMonitorParticleYPt *cutPass = new AliFemtoCutMonitorParticleYPt("cutPass", 0.13957);
  AliFemtoCutMonitorParticleYPt *cutFail = new AliFemtoCutMonitorParticleYPt("cutFail", 0.13957);
  dtc->AddCutMonitor(cutPass, cutFail);

  // Pair selector
  AliFemtoShareQualityTPCEntranceSepPairCut *sqpc = new AliFemtoShareQualityTPCEntranceSepPairCut();
  // remove split track pairs and pairs that share hits
  
  // Set maximim allowed "quality" for the pair
  //  1.0 - accept all pairs
  // -0.5 - reject all pairs
  // a reasonable value should lie between 0.0 and 0.5
  sqpc->SetShareQualityMax(1.0);

  // Set maximum allowed shared hits fraction per pair
  //  1.0 - accept all pairs
  //  0.0 - reject all pairs
  // a reasonable value is small but nno-zero (0.05)
  sqpc->SetShareFractionMax(1.0);

  // Set minimum allowed separation between nominal TPC entrance points
  // of the two tracks in the pair
  // 0.0 - accept all pairs
  // a reasonable value is 3.0 [cm]
  sqpc->SetTPCEntranceSepMinimum(0.0);
  sqpc->SetRemoveSameLabel(kFALSE);

  // Add the cuts to the analysis
  an->SetEventCut(mec);
  an->SetFirstParticleCut(dtc);
  an->SetSecondParticleCut(dtc);
  an->SetPairCut(sqpc);
  
  // Setup correlation functions
  // A simple qinv correlation function
  AliFemtoQinvCorrFctn *cqinv= new AliFemtoQinvCorrFctn("qinvcf",75,0.0,0.4);
  
  // A correlation function to monitor the splitting and cluster sharing in TPC
  AliFemtoShareQualityCorrFctn *csqqinv= new AliFemtoShareQualityCorrFctn("sqqinvcf",75,0.0,0.4);
  
  // A correlation function to monitor the distance at the entrance to the TPC
  AliFemtoTPCInnerCorrFctn *tpcin = new AliFemtoTPCInnerCorrFctn("tpcin",80, 0.0, 0.4);
  
  // add the correlation functions to the analysis
  an->AddCorrFctn(cqinv);
  an->AddCorrFctn(csqqinv);
  an->AddCorrFctn(tpcin);

  // Add the analysis to the manager
  Manager->AddAnalysis(an);	

  // Run the event loop
  long nE= 100000;
  if (Manager->Init())
    cout<<" Problem"<<endl;
	
  int Status=0;
  long int nEP=0;
  while ((!Status)&&(nEP<nE))
    {
      nEP++;
      cout<<" next event "<<nEP<<endl;
      Status=Manager->ProcessEvent();
    } 
	
  // Save the results
  char ofname[200];
  sprintf(ofname, "QinvCF.In.root");

  TFile f1 (ofname,"RECREATE","Data");
  cqinv->Numerator()->Write();
  cqinv->Denominator()->Write();
  csqqinv->WriteHistos();
  tpcin->WriteHistos();
  cutPass->Write();
  cutFail->Write();

  // Save the cut settings in the output file 
  f1.mkdir("Settings");
  f1.cd("Settings");
  
  TList *tListSettings = an->ListSettings();
  tListSettings->Write();

  f1.Close();

  TFile f2("Listout.root","RECREATE");
  TList *tOutList = an->GetOutputList();
  tOutList->Write();
}
