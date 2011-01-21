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
Bool_t RsnConfigKStar
(
  const char *taskName, 
  const char *options,
  const char *path
)
{
  // retrieve analysis manager & task and exit in case of failure
  AliAnalysisManager *mgr  = AliAnalysisManager::GetAnalysisManager();
  AliRsnAnalysisSE   *task = (AliRsnAnalysisSE*)mgr->GetTask(taskName);
  if (!task)
  {
    Error("RsnConfigPhi", "Task not found");
    return kFALSE;
  }
  
  // load useful macros
  gROOT->LoadMacro(Form("%s/QualityCutsITS.C", path));
  gROOT->LoadMacro(Form("%s/QualityCutsTPC.C", path));
  
  // interpret the useful information from second argument
  TString opt(options);
  TString suffix("_pid");
  Bool_t isSim    = opt.Contains("sim");
  Bool_t isData   = opt.Contains("data");
  Bool_t addITSSA = opt.Contains("its");
  if (!isSim && !isData)
  {
    Error("RsnConfigPhi", "Required to know if working on data or MC");
    return kFALSE;
  }
  if (addITSSA) suffix += "_its";
  
  // define some names using a standard part plus the above suffix
  TString pairPMname(suffix);
  TString pairMPname(suffix);
  TString truePMname(suffix);
  TString trueMPname(suffix);
  TString pairPPname(suffix);
  TString pairMMname(suffix);
  TString esdCutName(suffix);
  pairPMname.Prepend("PairPM");
  pairMPname.Prepend("PairMP");
  truePMname.Prepend("TruePM");
  trueMPname.Prepend("TrueMP");
  pairPPname.Prepend("PairPP");
  pairMMname.Prepend("PairMM");
  esdCutName.Prepend("cutESD2010");
  
  //
  // -- Setup pair definition (phi decay trees)  and pairs ------------------------------------------
  //

  // decay channels
  AliRsnPairDef *pairDefPM = new AliRsnPairDef(AliPID::kKaon, '+', AliPID::kPion, '-', 313, 0.896);
  AliRsnPairDef *pairDefMP = new AliRsnPairDef(AliPID::kKaon, '+', AliPID::kPion, '-', 313, 0.896);
  AliRsnPairDef *pairDefPP = new AliRsnPairDef(AliPID::kKaon, '+', AliPID::kPion, '+', 313, 0.896);
  AliRsnPairDef *pairDefMM = new AliRsnPairDef(AliPID::kKaon, '-', AliPID::kPion, '-', 313, 0.896);

  // computation objects
  AliRsnPairFunctions *pairPM = new AliRsnPairFunctions(pairPMname.Data(), pairDefPM);
  AliRsnPairFunctions *pairMP = new AliRsnPairFunctions(pairMPname.Data(), pairDefMP);
  AliRsnPairFunctions *truePM = new AliRsnPairFunctions(truePMname.Data(), pairDefPM);
  AliRsnPairFunctions *trueMP = new AliRsnPairFunctions(trueMPname.Data(), pairDefMP);
  AliRsnPairFunctions *pairPP = new AliRsnPairFunctions(pairPPname.Data(), pairDefPP);
  AliRsnPairFunctions *pairMM = new AliRsnPairFunctions(pairMMname.Data(), pairDefMM);
  
  // set additional option for true pairs
  truePM->SetOnlyTrue  (kTRUE);
  truePM->SetCheckDecay(kTRUE);
  trueMP->SetOnlyTrue  (kTRUE);
  trueMP->SetCheckDecay(kTRUE);
  
  //
  // -- Setup event cuts ----------------------------------------------------------------------------
  //
  
  // cut on primary vertex:
  // - 2nd argument = 10.0     --> require |Vz| <= 10 cm
  // - 3rd argument = 0        --> disable check on number of contributors
  // - 4th argument = 'kFALSE' --> reject TPC stand-alone primary vertex
  AliRsnCutPrimaryVertex *cutVertex  = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 0, kFALSE);
  
  // check pile-up with SPD
  cutVertex->SetCheckPileUp(kTRUE);

  //
  // -- Setup single track cuts ---------------------------------------------------------------------
  //
  
  // set Bethe-Bloch parameterizations
  // that depend on the sample
  Double_t bbPar[5];
  if (isSim)
  {
    bbPar[0] = 2.15898 / 50.0;
    bbPar[1] = 1.75295E1;
    bbPar[2] = 3.40030E-9;
    bbPar[3] = 1.96178;
    bbPar[4] = 3.91720;
  }
  else
  {
    bbPar[0] = 1.41543 / 50.0;
    bbPar[1] = 2.63394E1;
    bbPar[2] = 5.0411E-11;
    bbPar[3] = 2.12543;
    bbPar[4] = 4.88663;
  }

  // cut on track quality/PID:
  // 2nd argument (variable) --> 'kTRUE' for MonteCarlo, 'kFALSE' for data
  AliRsnCutESD2010 *cuts2010[2];
  cuts2010[0] = new AliRsnCutESD2010(esdCutName.Data(), isSim);
  cuts2010[1] = new AliRsnCutESD2010(esdCutName.Data(), isSim);
  
  // specify that PID is for kaons or pions
  cuts2010[0]->SetPID(AliPID::kKaon);
  cuts2010[1]->SetPID(AliPID::kPion);
  
  // other specifications are common
  for (Int_t i = 0; i < 2; i++)
  {
    // TPC+ITS tracks are always added, ITS-SA depends on options
    cuts2010[i]->SetUseITSTPC(kTRUE);
    cuts2010[i]->SetUseITSSA (addITSSA);
    
    // set quality cuts according to standard macro (defined outside)
    cuts2010[i]->CopyCutsTPC(QualityCutsTPC());
    cuts2010[i]->CopyCutsITS(QualityCutsITS());
    
    // add always the PID check
    cuts2010[i]->SetCheckITS(kTRUE);
    cuts2010[i]->SetCheckTPC(kTRUE);
    cuts2010[i]->SetCheckTOF(kTRUE);
    
    // set the ITS PID-related variables
    cuts2010[i]->SetITSband(3.0);
    
    // set the TPC PID-related variables
    cuts2010[i]->SetTPCrange(3.0, 5.0);
    cuts2010[i]->SetTPCpLimit(0.35);
    cuts2010[i]->GetESDpid()->GetTPCResponse().SetBetheBlochParameters(bbPar[0], bbPar[1], bbPar[2], bbPar[3], bbPar[4]);
    
    // set the TOF PID-related variables
    cuts2010[i]->SetTOFrange(-3.0, 3.0);
  }
  
  //
  // -- Setup track pair cuts -----------------------------------------------------------------------
  //
  
  // rapidity window
  AliRsnCutValue *cutY = new AliRsnCutValue("cutY", AliRsnValue::kPairY, -0.5, 0.5);
  
  // this cut requires a support object to retrieve default mass
  cutY->GetValueObj()->SetSupportObject(pairDefPM);
  
  //
  // -- Add cuts to cut collections -----------------------------------------------------------------
  //
  
  // event cuts are added directly to task, and only the first time
  if (task->GetEventCuts()->GetCuts()->FindObject("cutVertex") == 0x0)
  {
    task->GetEventCuts()->AddCut(cutVertex);
    task->GetEventCuts()->SetCutScheme("cutVertex");
  }
  
  // single track cuts are added to separate sets, in this case
  
  // first candidate is always the kaon
  pairPM->GetCutManager()->GetDaughter1Cuts()->AddCut(cuts2010[0]);
  truePM->GetCutManager()->GetDaughter1Cuts()->AddCut(cuts2010[0]);
  pairMP->GetCutManager()->GetDaughter1Cuts()->AddCut(cuts2010[0]);
  trueMP->GetCutManager()->GetDaughter1Cuts()->AddCut(cuts2010[0]);
  pairPP->GetCutManager()->GetDaughter1Cuts()->AddCut(cuts2010[0]);
  pairMM->GetCutManager()->GetDaughter1Cuts()->AddCut(cuts2010[0]);
  
  pairPM->GetCutManager()->GetDaughter1Cuts()->SetCutScheme(cuts2010[0]->GetName());
  truePM->GetCutManager()->GetDaughter1Cuts()->SetCutScheme(cuts2010[0]->GetName());
  pairMP->GetCutManager()->GetDaughter1Cuts()->SetCutScheme(cuts2010[0]->GetName());
  trueMP->GetCutManager()->GetDaughter1Cuts()->SetCutScheme(cuts2010[0]->GetName());
  pairPP->GetCutManager()->GetDaughter1Cuts()->SetCutScheme(cuts2010[0]->GetName());
  pairMM->GetCutManager()->GetDaughter1Cuts()->SetCutScheme(cuts2010[0]->GetName());
                                       
  // second candidate is always the pion
  pairPM->GetCutManager()->GetDaughter2Cuts()->AddCut(cuts2010[1]);
  truePM->GetCutManager()->GetDaughter2Cuts()->AddCut(cuts2010[1]);
  pairMP->GetCutManager()->GetDaughter2Cuts()->AddCut(cuts2010[1]);
  trueMP->GetCutManager()->GetDaughter2Cuts()->AddCut(cuts2010[1]);
  pairPP->GetCutManager()->GetDaughter2Cuts()->AddCut(cuts2010[1]);
  pairMM->GetCutManager()->GetDaughter2Cuts()->AddCut(cuts2010[1]);
  
  pairPM->GetCutManager()->GetDaughter2Cuts()->SetCutScheme(cuts2010[1]->GetName());
  truePM->GetCutManager()->GetDaughter2Cuts()->SetCutScheme(cuts2010[1]->GetName());
  pairMP->GetCutManager()->GetDaughter2Cuts()->SetCutScheme(cuts2010[1]->GetName());
  trueMP->GetCutManager()->GetDaughter2Cuts()->SetCutScheme(cuts2010[1]->GetName());
  pairPP->GetCutManager()->GetDaughter2Cuts()->SetCutScheme(cuts2010[1]->GetName());
  pairMM->GetCutManager()->GetDaughter2Cuts()->SetCutScheme(cuts2010[1]->GetName());
  
  // pair cuts are added to each pair
  TString scheme(cutY->GetName());
  pairPM->GetCutManager()->GetMotherCuts()->AddCut(cutY);
  truePM->GetCutManager()->GetMotherCuts()->AddCut(cutY);
  pairMP->GetCutManager()->GetMotherCuts()->AddCut(cutY);
  trueMP->GetCutManager()->GetMotherCuts()->AddCut(cutY);
  pairPP->GetCutManager()->GetMotherCuts()->AddCut(cutY);
  pairMM->GetCutManager()->GetMotherCuts()->AddCut(cutY);
  
  pairPM->GetCutManager()->GetMotherCuts()->SetCutScheme(cutY->GetName());
  truePM->GetCutManager()->GetMotherCuts()->SetCutScheme(cutY->GetName());
  pairMP->GetCutManager()->GetMotherCuts()->SetCutScheme(cutY->GetName());
  trueMP->GetCutManager()->GetMotherCuts()->SetCutScheme(cutY->GetName());
  pairPP->GetCutManager()->GetMotherCuts()->SetCutScheme(cutY->GetName());
  pairMM->GetCutManager()->GetMotherCuts()->SetCutScheme(cutY->GetName());

  //
  // -- Setup functions -----------------------------------------------------------------------------
  //

  // axis definition
  // 0) invariant mass
  // 1) transverse momentum
  // 2) multiplicity (variable binning)
  Double_t     mult[]   = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 30., 35., 40., 50., 60., 70., 80., 90., 100., 120., 140., 160., 180., 200., 1E+8};
  Int_t        nmult    = sizeof(mult) / sizeof(mult[0]);
  AliRsnValue *axisIM   = new AliRsnValue("IM", AliRsnValue::kPairInvMass     , 0.9,   1.4, 0.001);
  AliRsnValue *axisPt   = new AliRsnValue("PT", AliRsnValue::kPairPt          , 0.0,   5.0, 0.100);
  AliRsnValue *axisMult = new AliRsnValue("M" , AliRsnValue::kEventMultESDCuts, nmult, mult);

  // multiplicity axis needs a support object
  // of type AliESDtrackCuts, correctly configured
  AliESDtrackCuts *cuts = new AliESDtrackCuts(QualityCutsTPC());
  axisMult->SetSupportObject(cuts);

  // create function and add axes
  AliRsnFunction *fcn = new AliRsnFunction;
  if ( !fcn->AddAxis(axisIM  ) ) return kFALSE;
  if ( !fcn->AddAxis(axisPt  ) ) return kFALSE;
  if ( !fcn->AddAxis(axisMult) ) return kFALSE;

  // add functions to pairs
  pairPM->AddFunction(fcn);
  truePM->AddFunction(fcn);
  pairMP->AddFunction(fcn);
  trueMP->AddFunction(fcn);
  pairPP->AddFunction(fcn);
  pairMM->AddFunction(fcn);

  //
  // -- Conclusion ----------------------------------------------------------------------------------
  //

  // add all created AliRsnPair objects to the AliRsnAnalysisManager in the task
  task->GetAnalysisManager()->Add(pairPM);
  task->GetAnalysisManager()->Add(pairMP);
  task->GetAnalysisManager()->Add(pairPP);
  task->GetAnalysisManager()->Add(pairMM);
  if (isSim) 
  {
    task->GetAnalysisManager()->Add(truePM);
    task->GetAnalysisManager()->Add(trueMP);
  }

  return kTRUE;
}
