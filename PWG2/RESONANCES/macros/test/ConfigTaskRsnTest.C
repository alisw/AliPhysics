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
  Bool_t isESD   = strDataLabel->Contains("ESD");
  Bool_t isAOD   = strDataLabel->Contains("AOD");
  Bool_t isSim   = strDataLabel->Contains("sim");
  Bool_t isData  = strDataLabel->Contains("data");
  Bool_t isPass1 = strDataLabel->Contains("pass1");
  Bool_t isPass2 = strDataLabel->Contains("pass2");

  //
  // -- Set cuts for events (applied to all analyses) -----------------------------------------------
  //
  
  // primary vertex range
  AliRsnCutPrimaryVertex *cutVertex   = new AliRsnCutPrimaryVertex("cutVertex", 3);
  AliRsnCutSet           *cutSetEvent = new AliRsnCutSet("eventCuts", AliRsnCut::kEvent);
  cutSetEvent->AddCut(cutVertex);
  cutSetEvent->SetCutScheme("cutVertex");
  task->SetEventCuts(cutSetEvent);

  //
  // -- Setup pairs ---------------------------------------------------------------------------------
  //

  // decay channels
  AliRsnPairDef         *pairDef = new AliRsnPairDef(AliPID::kKaon, '+', AliPID::kKaon, '-', 333, 1.019455);

  // computation objects
  AliRsnPairFunctions   *pairPMhist = new AliRsnPairFunctions("pairHist", pairDef);
  AliRsnPairNtuple      *pairPMntp  = new AliRsnPairNtuple   ("pairNtp" , pairDef);

  //
  // -- Setup cuts ----------------------------------------------------------------------------------
  //
  
  // -- track cut --
  // --> global cuts for 2010 analysis
  AliRsnCutESD2010 *cuts2010 = new AliRsnCutESD2010("cuts2010");
  // ----> set the flag for sim/data management
  cuts2010->SetMC(isSim);
  // ----> require to check PID
  cuts2010->SetCheckITS(1);
  cuts2010->SetCheckTPC(1);
  cuts2010->SetCheckTOF(1);
  // ----> set TPC ranges and calibration
  cuts2010->SetTPCrange(5.0, 3.0);
  cuts2010->SetTPCpLimit(0.35);
  cuts2010->SetITSband(4.0);
  if (isMC) cuts2010->SetTPCpar(2.15898 / 50.0, 1.75295E1, 3.40030E-9, 1.96178, 3.91720);
  else      cuts2010->SetTPCpar(1.41543 / 50.0, 2.63394E1, 5.0411E-11, 2.12543, 4.88663);
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
  cuts2010->GetCutsITS()->SetRequireITSStandAlone(kTRUE);
  cuts2010->GetCutsITS()->SetRequireITSPureStandAlone(kFALSE);
  cuts2010->GetCutsITS()->SetRequireITSRefit(kTRUE);
  cuts2010->GetCutsITS()->SetMinNClustersITS(4);
  cuts2010->GetCutsITS()->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  cuts2010->GetCutsITS()->SetMaxChi2PerClusterITS(1.);
  cuts2010->GetCutsITS()->SetMaxDCAToVertexXYPtDep("0.0595+0.0182/pt^1.55"); // DCA pt dependent
  cuts2010->GetCutsITS()->SetMaxDCAToVertexZ(1e6); // disabled
  cuts2010->GetCutsITS()->SetDCAToVertex2D(kFALSE); // each DCA is checked separately
  cuts2010->GetCutsITS()->SetRequireITSPid(kTRUE);
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
  cuts2010->Initialize();
  
  // -- tracks --> TPC global tracks selection
  AliRsnCutDaughterType *cutTypeTPC = new AliRsnCutDaughterType("cutTypeTPC", AliRsnCutDaughterType::kTrackTPC);
  
  // -- tracks --> ITS standalone tracks selection
  AliRsnCutDaughterType *cutTypeITS = new AliRsnCutDaughterType("cutTypeITS", AliRsnCutDaughterType::kTrackITSSA);
  
  // -- tracks --> primarity/quality (global tracks)
  AliRsnCutESDPrimary *cutQualityTPC = new AliRsnCutESDPrimary("cutQualityTPC");
  cutQualityTPC->GetCuts()->SetRequireTPCStandAlone(kTRUE);
  cutQualityTPC->GetCuts()->SetMinNClustersTPC(70);
  cutQualityTPC->GetCuts()->SetAcceptKinkDaughters(kFALSE);
  cutQualityTPC->GetCuts()->SetRequireTPCRefit(kTRUE);
  cutQualityTPC->GetCuts()->SetRequireITSRefit(kTRUE);
  cutQualityTPC->GetCuts()->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  cutQualityTPC->GetCuts()->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9"); // DCA pt dependent: 7*(0.0050+0.0060/pt0.9)
  cutQualityTPC->GetCuts()->SetMaxDCAToVertexZ(1e6);
  cutQualityTPC->GetCuts()->SetDCAToVertex2D(kFALSE);
  cutQualityTPC->GetCuts()->SetRequireSigmaToVertex(kFALSE);
  
  // -- tracks --> primary/quality (ITS standalone)
  AliRsnCutESDPrimary *cutQualityITS = new AliRsnCutESDPrimary("cutQualityITS");
  cutQualityITS->GetCuts()->SetRequireITSStandAlone(kTRUE);
  //cutQualityITS->GetCuts()->SetRequireITSPureStandAlone(kFALSE);
  cutQualityITS->GetCuts()->SetRequireITSRefit(kTRUE); 
  cutQualityITS->GetCuts()->SetMinNClustersITS(4);
  cutQualityITS->GetCuts()->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  cutQualityITS->GetCuts()->SetMaxChi2PerClusterITS(1.);
  cutQualityITS->GetCuts()->SetMaxDCAToVertexXYPtDep("0.0595+0.0182/pt^1.55");
  //cutQualityITS->GetCuts()->SetRequireITSPid(kTRUE);
  
  // -- tracks --> realistic PID
  AliRsnCutPID *cutPID = new AliRsnCutPID("cutPIDKaon", AliPID::kKaon, 0.9);
  cutPID->SetPrior(AliPID::kElectron, 0.02);
  cutPID->SetPrior(AliPID::kMuon,     0.02);
  cutPID->SetPrior(AliPID::kPion,     0.83);
  cutPID->SetPrior(AliPID::kKaon,     0.07);
  cutPID->SetPrior(AliPID::kProton,   0.06);
  
  // -- pairs --> window in eta
  AliRsnCutStd *cutEtaPair = new AliRsnCutStd("cutEtaPair", AliRsnCut::kMother, AliRsnCutStd::kEta, -0.9,  0.9);
   
  // cut sets
  AliRsnCutSet *cutSetDaughterCommon = new AliRsnCutSet("commonDaughterCuts", AliRsnCut::kDaughter);
  AliRsnCutSet *cutSetMother         = new AliRsnCutSet("motherCuts"        , AliRsnCut::kMother);

  // --> add related cuts
  cutSetDaughterCommon->AddCut(cuts2010);
  cutSetDaughterCommon->AddCut(cutTypeITS);
  cutSetDaughterCommon->AddCut(cutTypeTPC);
  cutSetDaughterCommon->AddCut(cutQualityITS);
  cutSetDaughterCommon->AddCut(cutQualityTPC);
  cutSetDaughterCommon->AddCut(cutPID);
  cutSetMother        ->AddCut(cutEtaPair);
  // --> define schemes
  cutSetDaughterCommon->SetCutScheme("cuts2010");
  cutSetMother        ->SetCutScheme("");
   
  // cut managers
  // define a proper name for each mult bin, to avoid omonyme output histos
  pairPMhist->GetCutManager()->SetCommonDaughterCuts(cutSetDaughterCommon);
  //pairPMhist->GetCutManager()->SetMotherCuts        (cutSetMother);
  pairPMntp ->GetCutManager()->SetCommonDaughterCuts(cutSetDaughterCommon);
  //pairPMntp ->GetCutManager()->SetMotherCuts        (cutSetMother);

  // function axes
  AliRsnValue *axisIM   = new AliRsnValue(AliRsnValue::kPairInvMass,   50,  0.9,  1.4);
  AliRsnValue *axisPt   = new AliRsnValue(AliRsnValue::kPairPt,        50,  0.0, 10.0);
  AliRsnValue *axisEta  = new AliRsnValue(AliRsnValue::kPairEta,       30, -1.5,  1.5);
  
  pairPMntp->AddValue(axisIM);
  pairPMntp->AddValue(axisPt);
  pairPMntp->AddValue(axisEta);

  // functions
  AliRsnFunction *fcn      = new AliRsnFunction;
  AliRsnFunction *fcnPt    = new AliRsnFunction;
  AliRsnFunction *fcnPtEta = new AliRsnFunction;
  // --> add axes
  fcn     ->AddAxis(axisIM);
  fcnPt   ->AddAxis(axisIM);
  fcnPt   ->AddAxis(axisPt);
  fcnPtEta->AddAxis(axisIM);
  fcnPtEta->AddAxis(axisPt);
  fcnPtEta->AddAxis(axisEta);
  
  pairPMhist->AddFunction(fcn);
  pairPMhist->AddFunction(fcnPt);
  pairPMhist->AddFunction(fcnPtEta);
  
  // add everything to pair manager
  task->GetAnalysisManager()->Add(pairPMhist);
  task->GetAnalysisManager()->Add(pairPMntp);

  return kTRUE;
}
