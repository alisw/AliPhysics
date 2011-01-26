//
// This function configures the entire task for all resonances the user is interested in.
// This is done by creating all configuration objects which are defined in the package.
//
#include <TString.h>

#include "AliPID.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisManager.h"

#include "AliRsnValue.h"
#include "AliRsnPairDef.h"
#include "AliRsnFunction.h"
#include "AliRsnPairFunctions.h"
#include "AliRsnAnalysisSE.h"

#include "AliRsnCutSet.h"
#include "AliRsnCutValue.h"
#include "AliRsnCutESD2010.h"
#include "AliRsnCutPrimaryVertex.h"

#include "QualityCutsITS.C"
#include "QualityCutsTPC.C"

AliRsnCutESD2010 *ConfigCutESD2010(Bool_t isSim, Bool_t addITSTPC, Bool_t addITSSA, AliPID::EParticleType pidType);
void              ConfigPairCuts(AliRsnCutSet* cutSet, AliRsnPairDef *pairDef);
Bool_t            ConfigFunctionIM(AliRsnPairFunctions* pair);
Bool_t            ConfigFunctionRes(AliRsnPairFunctions* pair);

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
Bool_t RsnConfigPhi
(
  const char *taskName, 
  const char *options,
  Bool_t      addEventCuts = kFALSE
)
{
  // retrieve analysis manager & task and exit in case of failure
  AliAnalysisManager *mgr  = AliAnalysisManager::GetAnalysisManager();
  AliRsnAnalysisSE   *task = (AliRsnAnalysisSE*)mgr->GetTask(taskName);
  if (!task)
  {
    Error("RsnConfigPhiNew", "Task not found");
    return kFALSE;
  }
  
  // interpret the useful information from second argument
  TString opt(options);
  TString suffix("_pid");
  Bool_t isSim    = opt.Contains("sim");
  Bool_t addITSSA = opt.Contains("its");
  if (isSim) suffix += "_sim";
  if (addITSSA) suffix += "_its";
  
  // define some names using a standard part plus the above suffix
  TString pairPMname(suffix);
  TString truePMname(suffix);
  TString pairPPname(suffix);
  TString pairMMname(suffix);
  pairPMname.Prepend("PairPM");
  truePMname.Prepend("TruePM");
  pairPPname.Prepend("PairPP");
  pairMMname.Prepend("PairMM");
  
  //
  // -- Setup pair definition (phi decay trees)  and pairs ------------------------------------------
  //

  // decay channels
  AliRsnPairDef *pairDefPM = new AliRsnPairDef(AliPID::kKaon, '+', AliPID::kKaon, '-', 333, 1.019455);
  AliRsnPairDef *pairDefPP = new AliRsnPairDef(AliPID::kKaon, '+', AliPID::kKaon, '+', 333, 1.019455);
  AliRsnPairDef *pairDefMM = new AliRsnPairDef(AliPID::kKaon, '-', AliPID::kKaon, '-', 333, 1.019455);

  // computation objects
  AliRsnPairFunctions *pairPM = new AliRsnPairFunctions(pairPMname.Data(), pairDefPM);
  AliRsnPairFunctions *truePM = new AliRsnPairFunctions(truePMname.Data(), pairDefPM);
  AliRsnPairFunctions *pairPP = new AliRsnPairFunctions(pairPPname.Data(), pairDefPP);
  AliRsnPairFunctions *pairMM = new AliRsnPairFunctions(pairMMname.Data(), pairDefMM);
  
  // set additional option for true pairs
  truePM->SetOnlyTrue  (kTRUE);
  truePM->SetCheckDecay(kTRUE);
  
  //
  // -- Setup event cuts ----------------------------------------------------------------------------
  //
  
  if (addEventCuts)
  {
    // cut on primary vertex:
    // - 2nd argument = 10.0     --> require |Vz| <= 10 cm
    // - 3rd argument = 0        --> disable check on number of contributors
    // - 4th argument = 'kFALSE' --> reject TPC stand-alone primary vertex
    AliRsnCutPrimaryVertex *cutVertex  = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 0, kFALSE);
    
    // check pile-up with SPD
    cutVertex->SetCheckPileUp(kTRUE);
    
    // event cuts are added to their cut set
    task->GetEventCuts()->AddCut(cutVertex);
    task->GetEventCuts()->SetCutScheme("cutVertex");
  }

  //
  // -- Setup single track cuts ---------------------------------------------------------------------
  //

  // cut on track quality/PID:
  // 2nd argument (variable) --> 'kTRUE' for MonteCarlo, 'kFALSE' for data
  AliRsnCutESD2010 *cuts2010 = ConfigCutESD2010(isSim, kTRUE, addITSSA, AliPID::kKaon);
  
  // single track cuts are added to each pair as common ones for both daughters
  pairPM->GetCutManager()->GetCommonDaughterCuts()->AddCut(cuts2010);
  truePM->GetCutManager()->GetCommonDaughterCuts()->AddCut(cuts2010);
  pairPP->GetCutManager()->GetCommonDaughterCuts()->AddCut(cuts2010);
  pairMM->GetCutManager()->GetCommonDaughterCuts()->AddCut(cuts2010);
  
  pairPM->GetCutManager()->GetCommonDaughterCuts()->SetCutScheme(cuts2010->GetName());
  truePM->GetCutManager()->GetCommonDaughterCuts()->SetCutScheme(cuts2010->GetName());
  pairPP->GetCutManager()->GetCommonDaughterCuts()->SetCutScheme(cuts2010->GetName());
  pairMM->GetCutManager()->GetCommonDaughterCuts()->SetCutScheme(cuts2010->GetName());
  
  //
  // -- Setup track pair cuts -----------------------------------------------------------------------
  //
  
  ConfigPairCuts(pairPM->GetCutManager()->GetMotherCuts(), pairDefPM);
  ConfigPairCuts(truePM->GetCutManager()->GetMotherCuts(), pairDefPM);
  ConfigPairCuts(pairPP->GetCutManager()->GetMotherCuts(), pairDefPP);
  ConfigPairCuts(pairMM->GetCutManager()->GetMotherCuts(), pairDefMM);
  
  //
  // -- Setup functions -----------------------------------------------------------------------------
  //

  ConfigFunctionIM(pairPM);
  ConfigFunctionIM(truePM);
  ConfigFunctionIM(pairPP);
  ConfigFunctionIM(pairMM);

  ConfigFunctionRes(truePM);

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

//
// Macro for configuring the global set of cuts for 2010/2011 standard analysis
// for resonances, compliant with what proposed in the SPECTRA topical group
//
AliRsnCutESD2010 *ConfigCutESD2010(Bool_t isSim, Bool_t addITSTPC, Bool_t addITSSA, AliPID::EParticleType pidType)
{
  // interpret the options
  TString cutName("cutESD2010_pid");
  if (isSim) cutName += "_sim"; else cutName += "_data";
  if (addITSSA) cutName += "_its";
  cutName += AliPID::ParticleName(pidType);
  
  // instantiate cut object
  // 2nd argument (variable) --> 'kTRUE' for MonteCarlo, 'kFALSE' for data
  AliRsnCutESD2010 *cuts2010 = new AliRsnCutESD2010(cutName.Data(), isSim);
  
  // specify that PID is for kaons
  cuts2010->SetPID(pidType);
  
  // TPC+ITS tracks are always added, ITS-SA depends on options
  cuts2010->SetUseITSTPC(addITSTPC);
  cuts2010->SetUseITSSA (addITSSA);
  
  // set quality cuts according to standard macro (defined outside)
  cuts2010->CopyCutsTPC(QualityCutsTPC());
  cuts2010->CopyCutsITS(QualityCutsITS());
  
  // add always the PID check
  cuts2010->SetCheckITS(kTRUE);
  cuts2010->SetCheckTPC(kTRUE);
  cuts2010->SetCheckTOF(kTRUE);
  
  // set the ITS PID-related variables
  cuts2010->SetITSband(3.0);
  
  // set the TPC PID-related variables
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
  
  // set the TOF PID-related variables
  cuts2010->SetTOFrange(-3.0, 3.0);
  
  // now, cut is initialized properly
  return cuts2010;
}

//
// Utility macro to set the pair cuts
// (currently, only rapidity window)
//
void ConfigPairCuts(AliRsnCutSet* cutSet, AliRsnPairDef *pairDef)
{
  // rapidity window
  AliRsnCutValue *cutY = new AliRsnCutValue("cutY", AliRsnValue::kPairY, -0.5, 0.5);
  
  // this cut requires a support object to retrieve default mass
  cutY->GetValueObj()->SetSupportObject(pairDef);
  
  // add to the cut set and set the expression
  cutSet->AddCut(cutY);
  cutSet->SetCutScheme(cutY->GetName());
}

//
// Utility macro to setup the functions to be added
// ad add to the required pairs
//
Bool_t ConfigFunctionIM(AliRsnPairFunctions* pair)
{
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
  pair->AddFunction(fcn);
  
  // initialization OK
  return kTRUE;
}

//
// Utility macro to setup the functions to be added
// ad add to the required pairs
//
Bool_t ConfigFunctionRes(AliRsnPairFunctions* pair)
{
  // axis definition
  // 0) invariant mass
  // 1) transverse momentum
  // 2) multiplicity (variable binning)
  Double_t     mult[]   = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 30., 35., 40., 50., 60., 70., 80., 90., 100., 120., 140., 160., 180., 200., 1E+8};
  Int_t        nmult    = sizeof(mult) / sizeof(mult[0]);
  AliRsnValue *axisRes  = new AliRsnValue("RES", AliRsnValue::kPairInvMassRes  , -5.0,   5.0, 0.001);
  AliRsnValue *axisPt   = new AliRsnValue("PT" , AliRsnValue::kPairPt          ,  0.0,   5.0, 0.100);
  AliRsnValue *axisMult = new AliRsnValue("M"  , AliRsnValue::kEventMultESDCuts,  nmult, mult);

  // multiplicity axis needs a support object
  // of type AliESDtrackCuts, correctly configured
  AliESDtrackCuts *cuts = new AliESDtrackCuts(QualityCutsTPC());
  axisMult->SetSupportObject(cuts);

  // create function and add axes
  AliRsnFunction *fcn = new AliRsnFunction;
  if ( !fcn->AddAxis(axisRes ) ) return kFALSE;
  if ( !fcn->AddAxis(axisPt  ) ) return kFALSE;
  if ( !fcn->AddAxis(axisMult) ) return kFALSE;

  // add functions to pairs
  pair->AddFunction(fcn);
  
  // initialization OK
  return kTRUE;
}
