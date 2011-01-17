/*
#include <TROOT.h>
#include <TString.h>
#include <AliAnalysisManager.h>
#include <AliRsnAnalysisSE.h>
#include <AliRsnCutESD2010.h>
#include <AliRsnCutValue.h>
#include <AliRsnPairFunctions.h>
#include <AliRsnFunction.h>
#include <AliRsnCutPrimaryVertex.h>

#include "config/QualityCutsITS.C"
#include "config/QualityCutsTPC.C"
*/

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
Bool_t RsnConfig
(
  const char *taskName, 
  const char *options,
  const char *config,
  const char *path,
  Int_t       multMin = 0,
  Int_t       multMax = 0
)
{
  // load useful macros
  gROOT->LoadMacro(Form("%s/QualityCutsITS.C", path));
  gROOT->LoadMacro(Form("%s/QualityCutsTPC.C", path));
  
  // interpret the useful information from second argument
  TString opt(options);
  Bool_t isSim  = opt.Contains("sim");
  Bool_t isData = opt.Contains("data");
  if (!isSim && !isData)
  {
    Error("RsnConfig", "Required to know if working on data or MC");
    return kFALSE;
  }
  
  // interpret the specific info from third argument
  // which should be fixed in the various calls to this function
  TString conf(config);
  Bool_t addPID    = conf.Contains("pid");
  Bool_t addITSSA  = conf.Contains("its");
  Bool_t addDipCut = conf.Contains("dip");
      
  // generate a common suffix depending on chosen options
  TString suffix;
  if (addPID)    suffix += "_pid";
  if (addITSSA)  suffix += "_its";
  if (addDipCut) suffix += "_dip";
  Info("RsnConfig", "=== Specific configuration: %s ====================================================", config);
  Info("RsnConfig", "=== suffix used           : %s ====================================================", suffix.Data());

  // retrieve analysis manager & task
  AliAnalysisManager *mgr  = AliAnalysisManager::GetAnalysisManager();
  AliRsnAnalysisSE   *task = (AliRsnAnalysisSE*)mgr->GetTask(taskName);

  // for safety, return if no task is passed
  if (!task)
  {
    Error("RsnConfig2010PhiFcn", "Task not found");
    return kFALSE;
  }
  
  //
  // -- Setup event cuts (added directly to task) ---------------------------------------------------
  //
  
  // define a common cut on primary vertex, which also checks pile-up
  AliRsnCutPrimaryVertex *cutVertex  = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 0, kFALSE);
  cutVertex->SetCheckPileUp(kTRUE);
  task->GetEventCuts()->AddCut(cutVertex);
  task->GetEventCuts()->SetCutScheme("cutVertex");

  //
  // -- Setup pairs ---------------------------------------------------------------------------------
  //

  // decay channels
  AliRsnPairDef *pairDefPM = new AliRsnPairDef(AliPID::kKaon, '+', AliPID::kKaon, '-', 333, 1.019455);
  AliRsnPairDef *pairDefPP = new AliRsnPairDef(AliPID::kKaon, '+', AliPID::kKaon, '+', 333, 1.019455);
  AliRsnPairDef *pairDefMM = new AliRsnPairDef(AliPID::kKaon, '-', AliPID::kKaon, '-', 333, 1.019455);

  // computation objects
  AliRsnPairFunctions *pairPM = new AliRsnPairFunctions(Form("PairPM%s", suffix.Data()), pairDefPM);
  AliRsnPairFunctions *truePM = new AliRsnPairFunctions(Form("TruePM%s", suffix.Data()), pairDefPM);
  AliRsnPairFunctions *pairPP = new AliRsnPairFunctions(Form("PairPP%s", suffix.Data()), pairDefPP);
  AliRsnPairFunctions *pairMM = new AliRsnPairFunctions(Form("PairMM%s", suffix.Data()), pairDefMM);

  //
  // -- Setup cuts ----------------------------------------------------------------------------------
  //

  // track cut -----------------------
  // --> global cuts for 2010 analysis
  // --> most options are set to right values by default
  // --> second argument in constructor tells if we are working in simulation or not
  AliRsnCutESD2010 *cuts2010 = new AliRsnCutESD2010(Form("cuts2010%s", suffix.Data()), isSim);
  // --> set the reference particle for PID
  cuts2010->SetPID(AliPID::kKaon);
  // --> include or not the ITS standalone (TPC is always in)
  cuts2010->SetUseITSTPC(kTRUE);
  cuts2010->SetUseITSSA (addITSSA);
  // --> set the quality cuts using the general macro and using the 'Copy()' method in AliESDtrackCuts
  cuts2010->CopyCutsTPC(QualityCutsTPC());
  cuts2010->CopyCutsITS(QualityCutsITS());
  // --> set values for PID flags, depending on the choice expressed in the options
  cuts2010->SetCheckITS (addPID);
  cuts2010->SetCheckTPC (addPID);
  cuts2010->SetCheckTOF (addPID);
  // --> set the ITS PID-related variables
  cuts2010->SetITSband(3.0);
  // --> set the TPC PID-related variables
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
  cuts2010->SetTPCrange(3.0, 5.0);
  cuts2010->SetTPCpLimit(0.35);
  cuts2010->GetESDpid()->GetTPCResponse().SetBetheBlochParameters(bbPar[0], bbPar[1], bbPar[2], bbPar[3], bbPar[4]);
  // --> set the TOF PID-related variables
  cuts2010->SetTOFrange(-3.0, 3.0);
  
  // pair cut ----------------------------------------
  // --> dip angle between daughters: (it is a cosine)
  AliRsnCutValue *cutDip = new AliRsnCutValue("cutDip", AliRsnValue::kPairDipAngle,  0.02, 1.01);
  AliRsnCutValue *cutY   = new AliRsnCutValue("cutY"  , AliRsnValue::kPairY       , -0.5 , 0.5 );
  cutY->GetValueObj()->SetSupportObject(pairDefPM);

  // setup cut set for tracks------------------------------------------------------------
  // --> in this case, only common cuts are applied, depending if working with ESD or AOD
  // --> these cuts are added always
  pairPM->GetCutManager()->GetCommonDaughterCuts()->AddCut(cuts2010);
  truePM->GetCutManager()->GetCommonDaughterCuts()->AddCut(cuts2010);
  pairPP->GetCutManager()->GetCommonDaughterCuts()->AddCut(cuts2010);
  pairMM->GetCutManager()->GetCommonDaughterCuts()->AddCut(cuts2010);
  
  pairPM->GetCutManager()->GetCommonDaughterCuts()->SetCutScheme(cuts2010->GetName());
  truePM->GetCutManager()->GetCommonDaughterCuts()->SetCutScheme(cuts2010->GetName());
  pairPP->GetCutManager()->GetCommonDaughterCuts()->SetCutScheme(cuts2010->GetName());
  pairMM->GetCutManager()->GetCommonDaughterCuts()->SetCutScheme(cuts2010->GetName());
  
  // setup cut set for pairs---------------
  // --> add rapidity range cut
  TString scheme(cutY->GetName());
  pairPM->GetCutManager()->GetMotherCuts()->AddCut(cutY);
  truePM->GetCutManager()->GetMotherCuts()->AddCut(cutY);
  pairPP->GetCutManager()->GetMotherCuts()->AddCut(cutY);
  pairMM->GetCutManager()->GetMotherCuts()->AddCut(cutY);
  // --> add dip angle cut only if required
  if (addDipCut)
  {
    pairPM->GetCutManager()->GetMotherCuts()->AddCut(cutDip);
    truePM->GetCutManager()->GetMotherCuts()->AddCut(cutDip);
    pairPP->GetCutManager()->GetMotherCuts()->AddCut(cutDip);
    pairMM->GetCutManager()->GetMotherCuts()->AddCut(cutDip);
    scheme.Append(Form("&%s", cutDip->GetName()));
  }
  pairPM->GetCutManager()->GetMotherCuts()->SetCutScheme(scheme.Data());
  truePM->GetCutManager()->GetMotherCuts()->SetCutScheme(scheme.Data());
  pairPP->GetCutManager()->GetMotherCuts()->SetCutScheme(scheme.Data());
  pairMM->GetCutManager()->GetMotherCuts()->SetCutScheme(scheme.Data());
  
  // set additional option for true pairs
  truePM->SetOnlyTrue  (kTRUE);
  truePM->SetCheckDecay(kTRUE);

  //
  // -- Setup functions -----------------------------------------------------------------------------
  //

  // axis definition
  // 0) invariant mass
  // 1) transverse momentum
  // 2) multiplicity
  Double_t mult[] = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 30., 35., 40., 50., 60., 70., 80., 90., 100., 120., 140., 160., 180., 200., 1E+8};
  Int_t    nmult  = sizeof(mult) / sizeof(mult[0]);
  AliRsnValue *axisIM   = new AliRsnValue("IM", AliRsnValue::kPairInvMass     , 0.9,   1.4, 0.001);
  AliRsnValue *axisPt   = new AliRsnValue("PT", AliRsnValue::kPairPt          , 0.0,   5.0, 0.100);
  AliRsnValue *axisMult = new AliRsnValue("M" , AliRsnValue::kEventMultESDCuts, nmult, mult);

  // initialize the support object: AliESDtrackCuts
  // configured using the standard values
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
  pairPP->AddFunction(fcn);
  pairMM->AddFunction(fcn);

  //
  // -- Conclusion ----------------------------------------------------------------------------------
  //

  // add all created AliRsnPair objects to the AliRsnAnalysisManager in the task
  task->GetAnalysisManager()->Add(pairPM);
  task->GetAnalysisManager()->Add(pairPP);
  task->GetAnalysisManager()->Add(pairMM);
  if (isSim) task->GetAnalysisManager()->Add(truePM);

  return kTRUE;
}
