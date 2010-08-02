//
// This function configures the entire task for all resonances the user is interested in.
// This is done by creating all configuration objects which are defined in the package.
//
// Generally speaking, one has to define the following objects for each resonance:
//
//  1 - an AliRsnPairDef to define the resonance decay channel to be studied
//  2 - an AliRsnPair{Ntuple|Functions} where the output is stored
//  3 - one or more AliRsnCut objects to define track selections
//      which will have then to be organized into AliRsnCutSet objects
//  4 - an AliRsnCutManager to include all cuts to be applied (see point 3)
//  5 - definitions to build the TNtuple or histograms which are returned
//
// The return value is used to know if the configuration was successful
//
Bool_t RsnConfigTask(AliRsnAnalysisSE* &task, const char *dataLabel)
{
  // for safety, return if no task is passed
  if (!task)
  {
    Error("ConfigTaskRsn", "Task not found");
    return kFALSE;
  }
  
  // interpret the useful information from second argument
  TString strDataLabel(dataLabel);
  Bool_t isESD   = strDataLabel.Contains("ESD");
  Bool_t isAOD   = strDataLabel.Contains("AOD");
  Bool_t isSim   = strDataLabel.Contains("sim");
  Bool_t isData  = strDataLabel.Contains("data");
  Bool_t isPass1 = strDataLabel.Contains("pass1");
  Bool_t isPass2 = strDataLabel.Contains("pass2");
  Bool_t isPID   = strDataLabel.Contains("PID");

  //
  // -- Set cuts for events (applied to all analyses) -----------------------------------------------
  //
  
  // primary vertex range
  AliRsnCutPrimaryVertex *cutVertex   = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 0, kFALSE);
  AliRsnCutSet           *cutSetEvent = new AliRsnCutSet("eventCuts", AliRsnCut::kEvent);
  cutSetEvent->AddCut(cutVertex);
  cutSetEvent->SetCutScheme("cutVertex");
  task->SetEventCuts(cutSetEvent);

  //
  // -- Setup pairs ---------------------------------------------------------------------------------
  //

  // decay channels
  AliRsnPairDef         *pairDefPM = new AliRsnPairDef(AliPID::kKaon, '+', AliPID::kKaon, '-', 333, 1.019455);
  AliRsnPairDef         *pairDefPP = new AliRsnPairDef(AliPID::kKaon, '+', AliPID::kKaon, '+', 333, 1.019455);
  AliRsnPairDef         *pairDefMM = new AliRsnPairDef(AliPID::kKaon, '-', AliPID::kKaon, '-', 333, 1.019455);

  // computation objects
  //AliRsnPairNtuple      *pairPM = new AliRsnPairNtuple("PairPM" , pairDefPM);
  //AliRsnPairNtuple      *truePM = new AliRsnPairNtuple("TruePM" , pairDefPM);
  //AliRsnPairNtuple      *pairPP = new AliRsnPairNtuple("PairPP" , pairDefPP);
  //AliRsnPairNtuple      *pairMM = new AliRsnPairNtuple("PairMM" , pairDefMM);
  AliRsnPairFunctions    *pairPM = new AliRsnPairFunctions(Form("PairPM%s", (isPID ? "pid" : "nopid")), pairDefPM);
  AliRsnPairFunctions    *truePM = new AliRsnPairFunctions(Form("TruePM%s", (isPID ? "pid" : "nopid")), pairDefPM);
  AliRsnPairFunctions    *pairPP = new AliRsnPairFunctions(Form("PairPP%s", (isPID ? "pid" : "nopid")), pairDefPP);
  AliRsnPairFunctions    *pairMM = new AliRsnPairFunctions(Form("PairMM%s", (isPID ? "pid" : "nopid")), pairDefMM);

  //
  // -- Setup cuts ----------------------------------------------------------------------------------
  //
  
  // track cut -----------------------
  // --> global cuts for 2010 analysis
  AliRsnCutESD2010 *cuts2010 = new AliRsnCutESD2010(Form("cutESD2010%s", (isPID ? "pid" : "nopid")));
  // ----> set the flag for sim/data management
  cuts2010->SetMC(isSim);
  // ----> require to check PID or not, depending on the label
  cuts2010->SetCheckITS(isPID);
  cuts2010->SetCheckTPC(isPID);
  cuts2010->SetCheckTOF(isPID);
  // ----> set TPC rangeisPID and calibration
  cuts2010->SetTPCrange(5.0, 3.0);
  cuts2010->SetTPCpLimit(0.35);
  cuts2010->SetITSband(4.0);
  if (isSim) cuts2010->SetTPCpar(2.15898 / 50.0, 1.75295E1, 3.40030E-9, 1.96178, 3.91720);
  else       cuts2010->SetTPCpar(1.41543 / 50.0, 2.63394E1, 5.0411E-11, 2.12543, 4.88663);
  // ----> set standard quality cuts for TPC global tracks
  cuts2010->GetCutsTPC()->SetRequireTPCStandAlone(kTRUE); // require to have the projection at inner TPC wall
  cuts2010->GetCutsTPC()->SetMinNClustersTPC(70);
  cuts2010->GetCutsTPC()->SetMaxChi2PerClusterTPC(4);
  cuts2010->GetCutsTPC()->SetAcceptKinkDaughters(kFALSE);
  cuts2010->GetCutsTPC()->SetRequireTPCRefit(kTRUE);
  cuts2010->GetCutsTPC()->SetRequireITSRefit(kTRUE);
  cuts2010->GetCutsTPC()->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  cuts2010->GetCutsTPC()->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9"); // DCA pt dependent: 7*(0.0050+0.0060/pt0.9)
  cuts2010->GetCutsTPC()->SetMaxDCAToVertexZ(1e6); // disabled
  cuts2010->GetCutsTPC()->SetDCAToVertex2D(kFALSE); // each DCA is checked separately
  cuts2010->GetCutsTPC()->SetRequireSigmaToVertex(kFALSE);
  // ----> set standard quality cuts for ITS standalone tracks
  cuts2010->GetCutsITS()->SetRequireITSStandAlone(kTRUE, kTRUE);
  cuts2010->GetCutsITS()->SetRequireITSRefit(kTRUE);
  cuts2010->GetCutsITS()->SetMinNClustersITS(4);
  cuts2010->GetCutsITS()->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  cuts2010->GetCutsITS()->SetMaxChi2PerClusterITS(1.);
  cuts2010->GetCutsITS()->SetMaxDCAToVertexXYPtDep("0.0595+0.0182/pt^1.55"); // DCA pt dependent
  cuts2010->GetCutsITS()->SetMaxDCAToVertexZ(1e6); // disabled
  cuts2010->GetCutsITS()->SetDCAToVertex2D(kFALSE); // each DCA is checked separately
  // ----> set the configuration for TOF PID checks
  if (isData && (isPass1 || isPass2))
  {
    cuts2010->SetTOFcalibrateESD(kTRUE);
    //if (isPass2) cuts2010->SetTOFcalibrateESD(kFALSE); // potrebbe anche essere kFALSE
    cuts2010->SetTOFcorrectTExp(kTRUE);
    cuts2010->SetTOFuseT0(kTRUE);
    cuts2010->SetTOFtuneMC(kFALSE);
    cuts2010->SetTOFresolution(100.0);
  }
  else if (isSim)
  {
    cuts2010->SetTOFcalibrateESD(kFALSE);
    cuts2010->SetTOFcorrectTExp(kTRUE);
    cuts2010->SetTOFuseT0(kTRUE);
    cuts2010->SetTOFtuneMC(kTRUE);
    cuts2010->SetTOFresolution(100.0);
  }
  // ----> initialize cut (creates the AliESDpid, AliTOFcalib and AliTOFT0maker inside)
  cuts2010->Initialize();
  
  // cut sets ---------------------------------
  // --> only common cuts for tracks are needed
  AliRsnCutSet *cutSetDaughterCommon = new AliRsnCutSet("commonDaughterCuts", AliRsnCut::kDaughter);
  cutSetDaughterCommon->AddCut(cuts2010);
  cutSetDaughterCommon->SetCutScheme(cuts2010->GetName());
   
  // configure cut managers -------------------
  pairPM->GetCutManager()->SetCommonDaughterCuts(cutSetDaughterCommon);
  truePM->GetCutManager()->SetCommonDaughterCuts(cutSetDaughterCommon);
  pairPP->GetCutManager()->SetCommonDaughterCuts(cutSetDaughterCommon);
  pairMM->GetCutManager()->SetCommonDaughterCuts(cutSetDaughterCommon);
  
  // set additional option for true pairs when needed
  truePM->SetOnlyTrue();

  //
  // -- Setup functions -----------------------------------------------------------------------------
  //
  
  // function axes
  AliRsnValue *axisIM = new AliRsnValue("IM", AliRsnValue::kPairInvMass,   4000,  0.9,  4.9);
  AliRsnValue *axisPt = new AliRsnValue("PT", AliRsnValue::kPairPt,         100,  0.0, 10.0);
  AliRsnValue *axisY[4];
  axisY[0] = new AliRsnValue("Y1", AliRsnValue::kPairY, 1, -0.5,  0.5);
  axisY[1] = new AliRsnValue("Y2", AliRsnValue::kPairY, 1, -0.6,  0.6);
  axisY[2] = new AliRsnValue("Y3", AliRsnValue::kPairY, 1, -0.7,  0.7);
  axisY[3] = new AliRsnValue("Y4", AliRsnValue::kPairY, 1, -0.8,  0.8);
  
  // the ntuple output requires to get directly the values
  //pairPM->AddValue(axisIM);
  //pairPM->AddValue(axisPt);
  //pairPM->AddValue(axisEta);
  //truePM->AddValue(axisIM);
  //truePM->AddValue(axisPt);
  //truePM->AddValue(axisEta);
  //pairPP->AddValue(axisIM);
  //pairPP->AddValue(axisPt);
  //pairPP->AddValue(axisEta);
  //pairMM->AddValue(axisIM);
  //pairMM->AddValue(axisPt);
  //pairMM->AddValue(axisEta);
  
  // functions
  AliRsnFunction *fcnImPtY[4];
  for (Int_t i = 0; i < 4; i++)
  {
    fcnImPtY[i] = new AliRsnFunction;
    //fcnImPtY->UseSparse();
    // --> add axes
    fcnImPtY[i]->AddAxis(axisIM);
    fcnImPtY[i]->AddAxis(axisPt);
    fcnImPtY[i]->AddAxis(axisY[i]);
  
    // add functions to pairs
    pairPM->AddFunction(fcnImPtY[i]);
    truePM->AddFunction(fcnImPtY[i]);
    pairPP->AddFunction(fcnImPtY[i]);
    pairMM->AddFunction(fcnImPtY[i]);
  }
  
  //
  // -- Conclusion ----------------------------------------------------------------------------------
  //
  
  // add all created AliRsnPair objects to the AliRsnAnalysisManager in the task
  task->GetAnalysisManager()->Add(pairPM);
  task->GetAnalysisManager()->Add(truePM);
  task->GetAnalysisManager()->Add(pairPP);
  task->GetAnalysisManager()->Add(pairMM);

  return kTRUE;
}
