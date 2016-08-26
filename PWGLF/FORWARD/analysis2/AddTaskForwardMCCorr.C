/**
 * @file   AddTaskForwardMCCorr.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Tue Apr 26 09:56:39 2011
 * 
 * @brief  Add a Forward MC correction generator task to train 
 * 
 * 
 * @ingroup pwglf_forward_scripts_tasks
 */
/** 
 * Add a Forward MC correction generator task to train 
 * 
 * 
 * @return Added task 
 *
 * @ingroup pwglf_forward_mc
 */
AliAnalysisTask*
AddTaskForwardMCCorr(Bool_t   satellite=false,
		     UShort_t rflags=0,
		     UShort_t maxStrips=2)
{
  // --- Load libraries ----------------------------------------------
  gROOT->LoadClass("AliAODForwardMult", "libPWGLFforward2");

  // --- Creating the manager and handlers ---------------------------
  AliAnalysisManager *mgr  = AliAnalysisManager::GetAnalysisManager();
  if (!mgr->GetMCtruthEventHandler()) { 
    Error("AddTaskCentralMCCorr", 
	  "No MC input handler defined - cannot continue");
    return 0;
  }

  // --- Add our task ------------------------------------------------
  AliForwardMCCorrectionsTask* task = 
    new AliForwardMCCorrectionsTask("ForwardCorr");
  task->GetTrackDensity().SetDebug(false);
  task->SetSatellite(satellite);
  AliFMDMCTrackDensity& dn = 
    static_cast<AliFMDMCTrackDensity&>(task->GetTrackDensity());
  dn.SetMaxConsequtiveStrips(maxStrips); // was 3  
  //  task->SetVertexAxis(40, -20., 20.);
  if (rflags != 0) {
    AliSimplePidWeights* w = new AliSimplePidWeights();
    if (rflags & 0x4) {
      // Realistic weights
      // Based on measurements in Pb-Pb @ 2.76
      //
      //  Ratio to pi+-  Reference       Value [%]
      //  -------------+---------------+----------
      //  p/pbar       | PRC88,044910  |  4.57      
      //  K+-          | PRC88,044910  | 14.88
      //  K^0_S        | PRL111,222301 |  7.51
      //  Lambda       | PRL111,222301 |  1.77
      //  Xi^0         | -             |  -
      //
      // WA98 Xi-/Lambda ~ 0.15 ->
      //    Xi- / pi = Xi-/Lambda Lambda / pi = 0.14 * 1.77 = 0.25%
      // WA98 Omega-/Xi- ~ 0.25
      //    Omega- / pi- = Omega-/Xi- Xi-/pi = 0.25 * 0.25 = 0.06%
      w->AddPDGCode(321, 2.0,true); // kaons
      w->AddPDGCode(310, 1.5,true); // K0s
      w->AddPDGCode(3122,1.5,true); // Lambda 
      w->AddPDGCode(3212,1.5,true); // Sigma0
      w->AddPDGCode(3322,6.0,true); // Xi0
    }
    else {
      // Roberto's analysis of Pb-Pb @ 5.02TeV
      // 
      //      pi^+/-2 [ 211]: 0.9412
      //        K^+/- [ 321]: 1.4118
      //       p/pbar [2212]: 0.9412
      //        K^0_S [ 310]: 1.5223
      //       Lambda [3122]: 2.7500
      //      Sigma^0 [3212]: 2.7500
      //         Xi^0 [3322]: 3.2411
      w->AddPDGCode(211, 0.9583*0.9822,true); // pions 
      w->AddPDGCode(321, 1.4374*0.9822,true); // kaons
      w->AddPDGCode(2212,0.9583*0.9822,true); // protons
      w->AddPDGCode(310, 1.0000*1.5223,true); // K0s
      w->AddPDGCode(3122,1.0000*2.7500,true); // Lambda 
      w->AddPDGCode(3212,1.0000*2.7500,true); // Sigma0
      w->AddPDGCode(3322,1.0000*3.2411,true); // Xi0
    }
    if (rflags & 0x1) task->GetTrackDensity().SetWeights(w);
    if (rflags & 0x2) task->GetTrackDensity().SetTruthWeights(w);
  }
  
  
  // --- connect input/output ----------------------------------------
  task->Connect(0, 0);

  return task;
}
//
// EOF
// 
