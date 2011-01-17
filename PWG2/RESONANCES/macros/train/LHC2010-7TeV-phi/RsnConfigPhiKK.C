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
// This function configures a phi analysis task object.
// This is simpler than configuring an AliRsnAnalysiSE task, since it is aimed
// at a single type of resonance, but using all facilities of the framework.
//
// One has to define only the cuts and the function templates, which are 
// then replicated always for all possible types (+-, ++ and -- and trues)
//
// The return value is used to know if the configuration was successful
//
Bool_t RsnConfigPhiKK
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
  AliAnalysisManager  *mgr  = AliAnalysisManager::GetAnalysisManager();
  AliRsnAnalysisPhiKK *task = (AliRsnAnalysisPhiKK*)mgr->GetTask(taskName);

  // for safety, return if no task is passed
  if (!task)
  {
    Error("RsnConfigPhiKK", "Task not found");
    return kFALSE;
  }
  
  //
  // -- Setup event cuts ----------------------------------------------------------------------------
  //
  
  // define a common cut on primary vertex, which also checks pile-up
  AliRsnCutPrimaryVertex *cutVertex  = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 0, kFALSE);
  cutVertex->SetCheckPileUp(kTRUE);
  task->GetEventCuts()->AddCut(cutVertex);
  
  // if at least one of last two arguments is not zero, 
  // add a multiplicity cut using those arguments as range limits
  if (multMin > 0 || multMax > 0)
  {
    ::Info("RsnConfig.C", "Adding multiplicity cut: %d --> %d", multMin, multMax);
    AliRsnCutValue *cutMult = new AliRsnCutValue(Form("cutMult_%d-%d", multMin, multMax), AliRsnValue::kEventMultESDCuts, (Double_t)multMin, (Double_t)multMax);
    
    // initialize the support object: AliESDtrackCuts
    // configured using the standard values
    AliESDtrackCuts *cuts = new AliESDtrackCuts(QualityCutsTPC());
    cutMult->GetValueObj()->SetSupportObject(cuts);
    
    // add the cut and set the cut string
    task->GetEventCuts()->AddCut(cutMult);
    task->GetEventCuts()->SetCutScheme(Form("cutVertex&%s", cutMult->GetName()));
  }
  else
  {
    // if no mult cut is added, only primary vertex cut is used
    task->GetEventCuts()->SetCutScheme("cutVertex");
  }

  //
  // -- Setup other cuts (for tracks and pairs) -----------------------------------------------------
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
  cuts2010->SetUseITSSA(addITSSA);
  // --> set the quality cuts using the general macro and using the 'Copy()' method in AliESDtrackCuts
  cuts2010->CopyCutsTPC(QualityCutsTPC());
  cuts2010->CopyCutsITS(QualityCutsITS());
  // --> set values for PID flags, depending on the choice expressed in the options
  cuts2010->SetCheckITS(addPID);
  cuts2010->SetCheckTPC(addPID);
  cuts2010->SetCheckTOF(addPID);
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
  AliRsnCutValue *cutDip = new AliRsnCutValue("cutDip", AliRsnValue::kPairDipAngle, 0.02, 1.01);

  // setup cut set for tracks------------------------------------------------------------
  // --> in this case, only common cuts are applied, depending if working with ESD or AOD
  // --> these cuts are added always
  task->GetCommonDaughterCuts()->AddCut(cuts2010);
  task->GetCommonDaughterCuts()->SetCutScheme(cuts2010->GetName());
  
  // cut set for pairs---------------------
  // --> add dip angle cut only if required
  if (addDipCut)
  {
    task->GetMotherCuts()->AddCut(cutDip);
    task->GetMotherCuts()->SetCutScheme(cutDip->GetName());
  }

  //
  // -- Setup functions -----------------------------------------------------------------------------
  //

  // axis definition
  // 0) invariant mass
  // 1) transverse momentum
  // 2) rapidity
  AliRsnValue *axisIM = new AliRsnValue("IM", AliRsnValue::kPairInvMass,  0.9, 1.4, 0.001);
  AliRsnValue *axisPt = new AliRsnValue("PT", AliRsnValue::kPairPt     ,  0.0, 5.0, 0.100);
  AliRsnValue *axisY  = new AliRsnValue("Y" , AliRsnValue::kPairY      , -1.0, 1.0, 0.100);

  // create function and add axes
  AliRsnFunction *fcn = new AliRsnFunction;
  if ( !fcn->AddAxis(axisIM  ) ) return kFALSE;
  if ( !fcn->AddAxis(axisPt  ) ) return kFALSE;
  if ( !fcn->AddAxis(axisY   ) ) return kFALSE;

  // add functions to pairs
  task->AddFunction(fcn);

  //
  // -- Conclusion ----------------------------------------------------------------------------------
  //

  return kTRUE;
}
