//
// CONFIGURATION SCRIPT FOR RESONANCES ANALYSIS TRAIN
//
// The configuration script made of two functions:
//   1) RsnConfigPairManager
//   2) RsnConfigTask
//
// Function 'RsnConfigPairManager' defines a common scheme used
// for axis definition, cuts, and whatever is in common among
// all the output histograms.
// Function 'RsnConfigTask' calls the above function for all
// specific resonance pairs one wants to create.
// 

//_________________________________________________________________________________________________
AliRsnPairManager* RsnConfigPairManager
(
  const char            *pairMgrName,
  Int_t                  resonancePDG,
  Double_t               resonanceMass,
  AliPID::EParticleType  type1,
  AliPID::EParticleType  type2,
  Bool_t                 useY,
  Bool_t                 addTruePairs
)
//
// -- RsnConfigPairManager --
//
// This function defines:
//  - what pairs one wants
//  - what cuts have to be applied
//  - what axes the output histogram will have
//
// Output:
//  - an AliRsnPairManager containing all pairs related to one analysis
//   (usually, all pairs in the same manager refer to the same resonance
//    and this is assumed in the philosophy of the function itself)
//
// Arguments:
//  - pairMgrName   = name assigned to output object
//  - resonancePDG  = PDG code of resonance in study
//  - resonanceMass = nominal (PDG) mass of resonance in study (required for mass-based axes, like Y and Mt)
//  - type1         = first daughter type
//  - type2         = second daughter type
//  - addTruePairs  = flag to decide if true pairs histogram must be computed
//                    (requires MC info, so it is not feasible with AOD only)
//
// == NOTE 1 ==
// in this version of the configuration macro, the name given to PairManager must contain
// some keywords which define a choice of cuts which will be added to the analysis:
//  - "NOPID": completely no PID analysis (only primary track cuts are applied)
//  - "BB"   : all tracks are used, but the TPC Bethe-Bloch cut is applied (cut value = 0.2)
//  - "PID"  : realistic PID is used
//
// == NOTE 2 ==
// in this version of the configuration macro, the output histograms have the following axes:
//  0 = inv. mass
//  1 = transverse mass
//  2 = rapidity
//
{
  // === NAME DEFINITIONS =========================================================================
   
  AliRsnPairManager  *pairMgr  = new AliRsnPairManager(pairMgrName);
   
  // examines the given name to define details about track selection and cuts
  TString               str(pairMgrName);
  AliRsnPair::EPairType pidType;
  Bool_t                useBBCut;
  if (str.Contains("NOPID"))
  {
    pidType  = AliRsnPair::kNoPID;
    useBBCut = kFALSE;
    Info("RsnConfig", "PID TYPE = No PID        -- BB CUT: not used");
  }
  else if (str.Contains("BB"))
  {
    pidType  = AliRsnPair::kNoPID;
    useBBCut = kTRUE;
    Info("RsnConfig", "PID TYPE = No PID        -- BB CUT: used");
  }
  else if (str.Contains("PID"))
  {
    pidType  = AliRsnPair::kRealisticPID;
    useBBCut = kFALSE;
    Info("RsnConfig", "PID TYPE = Realistic PID -- BB CUT: not used");
  }
  else
  {
    Error("RsnConfig", "Unrecognized keywords in the name. Can't continue");
    return 0x0;
  }
   
  // === PAIR DEFINITIONS =========================================================================
   
  // if particle #1 and #2 are different, two histograms must be built
  // for each scheme (signal, true, mixed, like-sign) exchanging both particles and signs
  Int_t i, j, nArray = 1;
  if (type1 != type2) nArray = 2;
   
  AliRsnPairDef *defUnlike[2] = {0, 0};
  AliRsnPairDef *defLikePP[2] = {0, 0};
  AliRsnPairDef *defLikeMM[2] = {0, 0};
   
  defUnlike[0] = new AliRsnPairDef(type1, '+', type2, '-', resonancePDG, resonanceMass);
  defLikePP[0] = new AliRsnPairDef(type1, '+', type2, '+', resonancePDG, resonanceMass);
  defLikeMM[0] = new AliRsnPairDef(type1, '-', type2, '-', resonancePDG, resonanceMass);
   
  defUnlike[1] = new AliRsnPairDef(type2, '+', type1, '-', resonancePDG, resonanceMass);
  defLikePP[1] = new AliRsnPairDef(type2, '+', type1, '+', resonancePDG, resonanceMass);
  defLikeMM[1] = new AliRsnPairDef(type2, '-', type1, '-', resonancePDG, resonanceMass);
   
  // === PAIR ANALYSIS ENGINES ====================================================================
   
  // define null (dummy) objects and initialize only the ones which are needed,
  // depending again on particle types;
  // array is organized as follows:
  // [0] - true pairs
  // [1] - signal
  // [2] - like PP
  // [3] - like MM
  AliRsnPair *pair[2][4];
  for (i = 0; i < 2; i++) for (j = 0; j < 4; j++) pair[i][j] = 0x0;
   
  for (i = 0; i < nArray; i++) {
    pair[i][0] = new AliRsnPair(pidType, defUnlike[i]);
    pair[i][1] = new AliRsnPair(pidType, defUnlike[i]);
    pair[i][2] = new AliRsnPair(pidType, defLikePP[i]);
    pair[i][3] = new AliRsnPair(pidType, defLikeMM[i]);
  }
   
  // === CUTS =====================================================================================
   
  // cuts for tracks:
  // -- primary track quality
  AliRsnCutESDPrimary *cutESDPrimary = new AliRsnCutESDPrimary("cutESDPrimary");
  cutESDPrimary->GetCuts()->SetMaxCovDiagonalElements(2.0, 2.0, 0.5, 0.5, 2.0);
  cutESDPrimary->GetCuts()->SetRequireSigmaToVertex(kTRUE);
  cutESDPrimary->GetCuts()->SetMaxNsigmaToVertex(3.0);
  cutESDPrimary->GetCuts()->SetRequireTPCRefit(kTRUE);
  cutESDPrimary->GetCuts()->SetAcceptKinkDaughters(kFALSE);
  cutESDPrimary->GetCuts()->SetMinNClustersTPC(50);
  cutESDPrimary->GetCuts()->SetMaxChi2PerClusterTPC(3.5);
  // -- Bethe-Bloch with kaon mass hypothesis
  Double_t sigmaTPC = 0.065;
  AliRsnCutBetheBloch *cutKaonBB = new AliRsnCutBetheBloch("cutKaonBB", 3.0 * sigmaTPC, AliPID::kKaon);
  cutKaonBB->SetCalibConstant(0, 0.76176e-1);
  cutKaonBB->SetCalibConstant(1, 10.632);
  cutKaonBB->SetCalibConstant(2, 0.13279e-4);
  cutKaonBB->SetCalibConstant(3, 1.8631);
  cutKaonBB->SetCalibConstant(4, 1.9479);
   
  // cuts on pairs:
  // -- true daughters of a phi resonance (only for true pairs histogram)
  AliRsnCutStd *cutTruePair = new AliRsnCutStd("cutTrue", AliRsnCutStd::kTruePair, resonancePDG);
  // -- whole interval in Pt and Eta
  //AliRsnCutStd *cutEtaPair = new AliRsnCutStd("cutEtaPair", AliRsnCutStd::kEta, -0.9,  0.9);
  //AliRsnCutStd *cutPtPair  = new AliRsnCutStd("cutPtPair" , AliRsnCutStd::kPt ,  0.0, 10.0);
  //AliRsnCutStd *cutYPair   = new AliRsnCutStd("cutYPair"  , AliRsnCutStd::kEta, -0.9,  0.9);
  //cutYPair->SetMass(resonanceMass);
   
  // cuts on event (specific for this analysis):
  // -- whole interval in multiplicity
  AliRsnCutStd *cutMultiplicity = new AliRsnCutStd("cutMult", AliRsnCutStd::kMult, 0, 200);
   
  // cut set definition for all pairs
  AliRsnCutSet *cutSetParticle = new AliRsnCutSet("trackCuts");
  cutSetParticle->AddCut(cutESDPrimary);
  if (useBBCut)
  {
    cutSetParticle->AddCut(cutKaonBB);
    cutSetParticle->SetCutScheme("cutKaonBB&cutESDPrimary");
  }
  else 
  {
    cutSetParticle->SetCutScheme("cutESDPrimary");
  }
   
  // cut set definition for true pairs
  AliRsnCutSet *cutSetPairTrue = new AliRsnCutSet("truePairs");
  //cutSetPairTrue->AddCut(cutPtPair);
  //cutSetPairTrue->AddCut(cutEtaPair);
  //cutSetPairTrue->AddCut(cutTruePair);
  //cutSetPairTrue->SetCutScheme("cutPtPair&cutEtaPair&cutTrue");
  cutSetPairTrue->SetCutScheme("cutTrue");
   
  // cut set definition for all pairs
  //AliRsnCutSet *cutSetPairAll = new AliRsnCutSet("allPairs");
  //cutSetPairAll->AddCut(cutPtPair);
  //if (useY) 
  //{
  //  cutSetPairAll->AddCut(cutYPair); 
  //  cutSetPairAll->SetCutScheme("cutPtPair&cutYPair");
  //}
  //else 
  //{
  //  cutSetPairAll->AddCut(cutEtaPair);
  //  cutSetPairAll->SetCutScheme("cutPtPair&cutEtaPair");
  //}
   
  // cut set definition for events
  //AliRsnCutSet *cutSetEvent = new AliRsnCutSet("cutSetMult");
  //cutSetEvent->AddCut(cutMultiplicity);
  //cutSetEvent->SetCutScheme("cutMult");
   
  // cut manager for all pairs
  // define a proper name for each mult bin, to avoid omonyme output histos
  AliRsnCutMgr *cutMgrAll = new AliRsnCutMgr("std", "All");
  cutMgrAll->SetCutSet(AliRsnCut::kParticle, cutSetParticle);
  //cutMgrAll->SetCutSet(AliRsnCut::kPair, cutSetPairAll);
  //cutMgrAll->SetCutSet(AliRsnCut::kEvent, cutSetEvent);
  
  // cut manager for true pairs
  AliRsnCutMgr *cutMgrTrue = new AliRsnCutMgr("true", "True");
  cutMgrTrue->SetCutSet(AliRsnCut::kParticle, cutSetParticle);
  cutMgrTrue->SetCutSet(AliRsnCut::kPair, cutSetPairTrue);
  //cutMgrAll->SetCutSet(AliRsnCut::kEvent, cutSetEvent);
   
  for (i = 0; i < nArray; i++) {
    pair[i][0]->SetCutMgr(cutMgrTrue);
    pair[i][0]->SetOnlyTrue();
    for (j = 1; j < 4; j++) {
      pair[i][j]->SetCutMgr(cutMgrAll);
    }
  }
   
  // === FUNCTIONS ================================================================================
   
  // define all binnings
  AliRsnFunctionAxis *axisIM   = new AliRsnFunctionAxis(AliRsnFunctionAxis::kPairInvMass,    4000,  0.0,   2.0);
  AliRsnFunctionAxis *axisPt   = new AliRsnFunctionAxis(AliRsnFunctionAxis::kPairPt,           50,  0.0,  10.0);
  AliRsnFunctionAxis *axisY    = new AliRsnFunctionAxis(AliRsnFunctionAxis::kPairY,            20, -1.0,   1.0);
  AliRsnFunctionAxis *axisEta  = new AliRsnFunctionAxis(AliRsnFunctionAxis::kPairEta,          20, -1.0,   1.0);
  AliRsnFunctionAxis *axisMult = new AliRsnFunctionAxis(AliRsnFunctionAxis::kEventMult,         8,  0.0, 200.0);
   
  AliRsnFunction *fcnPt   = new AliRsnFunction();
  AliRsnFunction *fcnEta  = new AliRsnFunction();
  AliRsnFunction *fcnY    = new AliRsnFunction();
  AliRsnFunction *fcnMult = new AliRsnFunction();
  
  fcnPt  ->AddAxis(axisIM);
  fcnEta ->AddAxis(axisIM);
  fcnY   ->AddAxis(axisIM);
  fcnMult->AddAxis(axisIM);
   
  fcnPt  ->AddAxis(axisPt);
  fcnY   ->AddAxis(axisY); 
  fcnEta ->AddAxis(axisEta);
  fcnMult->AddAxis(axisMult);
   
  // add functions to pairs
  for (i = 0; i < nArray; i++) {
    for (j = 0; j < 4; j++) {
      pair[i][j]->AddFunction(fcnPt);
      pair[i][j]->AddFunction(fcnEta);
      pair[i][j]->AddFunction(fcnY);
      pair[i][j]->AddFunction(fcnMult);
    }
  }
   
  // === ADD TO PAIR MANAGER ======================================================================
  
  //
  // true pairs are array [0] element:
  // if they must be excluded, this element is skipped everywhere
  //
   
  for (i = 0; i < nArray; i++) {
    Int_t start = (addTruePairs ? 0 : 1);
    for (j = start; j < 4; j++) {
      pairMgr->AddPair(pair[i][j]);
    }
  }
   
  return pairMgr;
}

//_________________________________________________________________________________________________
void RsnConfigTask(AliRsnAnalysisSE* &task, Bool_t useY = kTRUE)
//
// -- RsnConfigTask --
//
// This function configures the entire task for all resonances the user is interested in.
// It recalls the above function many times to configure all required pair managers.
//
// By default, this function has the used task as argument.
//
{
  if (!task)
  {
    Error("ConfigTaskRsn", "Task not found");
    return;
  }

  // set prior probabilities for PID
  task->SetPriorProbability(AliPID::kElectron, 0.02);
  task->SetPriorProbability(AliPID::kMuon,     0.02);
  task->SetPriorProbability(AliPID::kPion,     0.83);
  task->SetPriorProbability(AliPID::kKaon,     0.07);
  task->SetPriorProbability(AliPID::kProton,   0.06);
  task->DumpPriors();
  
  // initialize analysis manager with pairs from config
  AliRsnAnalysisManager *anaMgr = task->GetAnalysisManager(0);
    
  // create pair managers for phi
  anaMgr->Add(RsnConfigPairManager("PHI_NOPID", 333, 1.0193, AliPID::kKaon, AliPID::kKaon, useY, kTRUE));
  anaMgr->Add(RsnConfigPairManager("PHI_BB"   , 333, 1.0193, AliPID::kKaon, AliPID::kKaon, useY, kTRUE));
  anaMgr->Add(RsnConfigPairManager("PHI_PID"  , 333, 1.0193, AliPID::kKaon, AliPID::kKaon, useY, kTRUE));
  // create pair managers for kstar
  anaMgr->Add(RsnConfigPairManager("KSTAR_NOPID", 313, 0.896, AliPID::kPion, AliPID::kKaon, useY, kTRUE));
  anaMgr->Add(RsnConfigPairManager("KSTAR_BB"   , 313, 0.896, AliPID::kPion, AliPID::kKaon, useY, kTRUE));
  anaMgr->Add(RsnConfigPairManager("KSTAR_PID"  , 313, 0.896, AliPID::kPion, AliPID::kKaon, useY, kTRUE));
  
  // setup cuts for events (good primary vertex)
  AliRsnCutPrimaryVertex *cutVertex   = new AliRsnCutPrimaryVertex("cutVertex", 3);
  AliRsnCutSet           *cutSetEvent = new AliRsnCutSet("eventCuts");
  cutSetEvent->AddCut(cutVertex);
  cutSetEvent->SetCutScheme("cutVertex");
  task->SetEventCuts(cutSetEvent);
}
