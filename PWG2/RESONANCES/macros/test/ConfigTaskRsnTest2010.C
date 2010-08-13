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
Bool_t RsnConfigTask(AliRsnAnalysisSE* &task, const char *dataLabel, Bool_t perfectPID, Bool_t useITS, Bool_t usePID)
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
  AliRsnPairDef         *pairDefpm = new AliRsnPairDef(AliPID::kKaon, '+', AliPID::kKaon, '-', 333, 1.019455);

  // computation objects
  AliRsnPairFunctions   *pairPMhist = new AliRsnPairFunctions("pairPMHist", pairDefpm);
  AliRsnPairNtuple      *pairPMntp  = new AliRsnPairNtuple   ("pairPMNtp" , pairDefpm);

  //
  // -- Setup cuts ----------------------------------------------------------------------------------
  //
  
  // -- track cut --
  // --> global cuts for 2010 analysis
  AliRsnCutESD2010 *cuts2010 = new AliRsnCutESD2010("cuts2010");
  // ----> set the flag for sim/data management
  cuts2010->SetMC(isSim);
  // ----> require to check PID
  cuts2010->SetCheckITS(kFALSE);
  cuts2010->SetCheckTPC(kFALSE);
  cuts2010->SetCheckTOF(kFALSE);
  // ----> set TPC ranges and calibration
  cuts2010->SetTPCrange(5.0, 3.0);
  cuts2010->SetTPCpLimit(0.35);
  cuts2010->SetITSband(4.0);
  if (isSim) cuts2010->SetTPCpar(2.15898 / 50.0, 1.75295E1, 3.40030E-9, 1.96178, 3.91720);
  else       cuts2010->SetTPCpar(1.41543 / 50.0, 2.63394E1, 5.0411E-11, 2.12543, 4.88663);
  // ----> set standard quality cuts for TPC global tracks
  //cuts2010->GetCutsTPC()->SetRequireTPCStandAlone(kTRUE); // require to have the projection at inner TPC wall
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
  cuts2010->Initialize();
  
  // -- tracks --> perfect PID
  AliRsnCutPID *cutPID = new AliRsnCutPID("cutPID", AliPID::kKaon, 0.0, kTRUE);
  
  // cut sets
  AliRsnCutSet *cutSetDaughterCommon = new AliRsnCutSet("commonDaughterCuts", AliRsnCut::kDaughter);

  // --> add related cuts
  cutSetDaughterCommon->AddCut(cuts2010);
  cutSetDaughterCommon->AddCut(cutPID);

  // --> define schemes
  if (perfectPID)
  {
    if (useITS)
    {
      pairPMhist->SetName("pmPerfectITS");
      cuts2010->SetUseITSSA(kTRUE);
      cuts2010->SetCheckITS(kFALSE);
      cuts2010->SetCheckTPC(kFALSE);
      cuts2010->SetCheckTOF(kFALSE);
      cutSetDaughterCommon->SetCutScheme("cuts2010&cutPID");
    }
    else
    {
      pairPMhist->SetName("pmPerfectNoITS");
      cuts2010->SetUseITSSA(kFALSE);
      cuts2010->SetCheckITS(kFALSE);
      cuts2010->SetCheckTPC(kFALSE);
      cuts2010->SetCheckTOF(kFALSE);
      cutSetDaughterCommon->SetCutScheme("cuts2010&cutPID");
    }
  }
  else
  {
    if (useITS)
    {
      if (usePID)
      {
        pairPMhist->SetName("pmExpITSpid");
        cuts2010->SetUseITSSA(kTRUE);
        cuts2010->SetCheckITS(kTRUE);
        cuts2010->SetCheckTPC(kTRUE);
        cuts2010->SetCheckTOF(kTRUE);
        cutSetDaughterCommon->SetCutScheme("cuts2010");
      }
      else
      {
        pairPMhist->SetName("pmExpITSnopid");
        cuts2010->SetUseITSSA(kTRUE);
        cuts2010->SetCheckITS(kFALSE);
        cuts2010->SetCheckTPC(kFALSE);
        cuts2010->SetCheckTOF(kFALSE);
        cutSetDaughterCommon->SetCutScheme("cuts2010");
      }
    }
    else
    {
      if (usePID)
      {
        pairPMhist->SetName("pmExpNoITSpid");
        cuts2010->SetUseITSSA(kFALSE);
        cuts2010->SetCheckITS(kTRUE);
        cuts2010->SetCheckTPC(kTRUE);
        cuts2010->SetCheckTOF(kTRUE);
        cutSetDaughterCommon->SetCutScheme("cuts2010");
      }
      else
      {
        pairPMhist->SetName("pmExpNoITSnopid");
        cuts2010->SetUseITSSA(kFALSE);
        cuts2010->SetCheckITS(kFALSE);
        cuts2010->SetCheckTPC(kFALSE);
        cuts2010->SetCheckTOF(kFALSE);
        cutSetDaughterCommon->SetCutScheme("cuts2010");
      }
    }
  }
  
  // cut managers
  // define a proper name for each mult bin, to avoid omonyme output histos
  pairPMhist->GetCutManager()->SetCommonDaughterCuts(cutSetDaughterCommon);
  pairPMntp ->GetCutManager()->SetCommonDaughterCuts(cutSetDaughterCommon);

  // function axes
  AliRsnValue *axisIM = new AliRsnValue("IM", AliRsnValue::kPairInvMass,        50,  0.9,  1.4);
  AliRsnValue *axisPt = new AliRsnValue("PT", AliRsnValue::kPairPt,             50,  0.0, 20.0);

  // functions for TH1-like output
  AliRsnFunction *fcnPt    = new AliRsnFunction;
  // --> add axes
  fcnPt   ->AddAxis(axisIM);
  fcnPt   ->AddAxis(axisPt);
  
  // add functions to TH1-like output
  pairPMhist->AddFunction(fcnPt);
  
  // add values to TNtuple-like output
  pairPMntp->AddValue(axisIM);
  pairPMntp->AddValue(axisPt);
  
  // add everything to analysis manager
  task->GetAnalysisManager()->Add(pairPMhist);
  //task->GetAnalysisManager()->Add(pairPMntp);

  return kTRUE;
}
