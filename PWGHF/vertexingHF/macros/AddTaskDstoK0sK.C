#include <algorithm>


AliAnalysisTaskSEDstoK0sK *AddTaskDstoK0sK(TString cutFileName,
                                           Bool_t  fReadMC,
                                           Bool_t  fFillNtuple,
                                           Bool_t  fUseSelectionBit,
                                           Int_t   AODProtection = 0,
                                           TString suffixName = "")
{
   //
   // AddTask for the AliAnalysisTaskSE for the D+ full resolution:
   // loop on V0s and bachelors, then on-the-fly D+ reconstruction and analysis
   // Authors:  J.Hamon, julien.hamon@cern.ch
   //


   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskResolutionK0sReco", "No analysis manager to connect to.");
      return NULL;
   }



   // Create the analysis cuts.
   //==============================================================================
   TFile *filecuts;
   filecuts = TFile::Open(cutFileName.Data());
   if (!filecuts || (filecuts && !(filecuts->IsOpen())))
      AliFatal("Cut object not found: analysis will not start!\n");

   AliRDHFCutsDstoK0sK* analysiscuts = new AliRDHFCutsDstoK0sK();
   analysiscuts = (AliRDHFCutsDstoK0sK*) filecuts->Get("AnalysisCuts");



   // Create the analysis task.
   //==============================================================================
   const Int_t nCutsTuple = 28;
   Float_t minCutsTuple[nCutsTuple]={0.}; std::fill_n(minCutsTuple, nCutsTuple, -1000.);
   Float_t maxCutsTuple[nCutsTuple]={0.}; std::fill_n(maxCutsTuple, nCutsTuple,  1000.);
   if (fFillNtuple) {
      Float_t massK0 = TDatabasePDG::Instance()->GetParticle(310)->Mass();
      Float_t massDs = TDatabasePDG::Instance()->GetParticle(431)->Mass();
      Float_t massDp = TDatabasePDG::Instance()->GetParticle(411)->Mass();

      minCutsTuple[ 0] = massK0-0.05;  maxCutsTuple[ 0] = massK0+0.05;  // "K0sInvMass"
      minCutsTuple[ 1] =  0.;          maxCutsTuple[ 1] =   40.;        // "K0sPt"
      minCutsTuple[ 2] =  0.;          maxCutsTuple[ 2] =  140.;        // "K0sRxy"
      minCutsTuple[ 3] =  0.;          maxCutsTuple[ 3] =   20.;        // "K0sCTau"
      minCutsTuple[ 4] =  0.99;        maxCutsTuple[ 4] =    1.;        // "K0sCosPA"
      minCutsTuple[ 5] =  0.;          maxCutsTuple[ 5] = 1000.;        // "K0sd0"
      minCutsTuple[ 6] = -40.;         maxCutsTuple[ 6] =   40.;        // "K0sd0DauPos"
      minCutsTuple[ 7] = -40.;         maxCutsTuple[ 7] =   40.;        // "K0sd0DauNeg"
      minCutsTuple[ 8] =  0.;          maxCutsTuple[ 8] =    3.;        // "K0sDCAdau"
      minCutsTuple[ 9] =  0.1;         maxCutsTuple[ 9] =   40.;        // "K0sPtDauPos"
      minCutsTuple[10] =  0.1;         maxCutsTuple[10] =   40.;        // "K0sPtDauNeg"
      minCutsTuple[11] =  0.3;         maxCutsTuple[11] =   40.;        // "BachPt"
      minCutsTuple[12] = -3.;          maxCutsTuple[12] =    3.;        // "Bachd0"
      minCutsTuple[13] = massDs-0.2;   maxCutsTuple[13] = massDs+0.2;   // "CanInvMassDs"
      minCutsTuple[14] = massDp-0.2;   maxCutsTuple[14] = massDp+0.2;   // "CanInvMassDplus"
      minCutsTuple[15] =  1.;          maxCutsTuple[15] =  24.;         // "CanPt"
      minCutsTuple[16] =  0.;          maxCutsTuple[16] =   8.;         // "CanDCAProngToProng"
      minCutsTuple[17] = -1.;          maxCutsTuple[17] =   1.;         // "CanCosThetaStarK0s"
      minCutsTuple[18] = -1.;          maxCutsTuple[18] =   1.;         // "CanCosThetaStarBach"
      minCutsTuple[19] =  0.;          maxCutsTuple[19] =   1.;         // "CanCosPA"
      minCutsTuple[20] = -1.;          maxCutsTuple[20] =   1.;         // "CanCosPAxy"
      minCutsTuple[21] =  0.;          maxCutsTuple[21] =   3.;         // "CanDLengthXY"
      minCutsTuple[22] =  0.;          maxCutsTuple[22] = 1000.;        // "CanNormDLengthXY"
      minCutsTuple[23] =  0.;          maxCutsTuple[23] =   3.;         // "CanDLength3D"
      minCutsTuple[24] =  0.;          maxCutsTuple[24] = 1000.;        // "CanNormDLength3D"
      minCutsTuple[25] =  0.;          maxCutsTuple[25] = 1000.;        // "CanSigmaVtx"
      minCutsTuple[26] = -1000.;       maxCutsTuple[26] = 1000.;        // "CanNormTopoBach"
      minCutsTuple[27] = -1000.;       maxCutsTuple[27] = 1000.;        // "CanNormTopoK0s"
   }


   AliAnalysisTaskSEDstoK0sK *task = new AliAnalysisTaskSEDstoK0sK("AnalysisDs", analysiscuts, fReadMC, fFillNtuple, nCutsTuple, minCutsTuple, maxCutsTuple);
   if (!task)
      AliFatal("Analysis task AliAnalysisTaskSEDstoK0sK not found (NULL pointer)");
   task->SetAODMismatchProtection(AODProtection);
   task->SetUseSelectionBit(fUseSelectionBit);
   mgr->AddTask(task);



   // Create the input/output containers.
   //==============================================================================
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());


   // Output slots #0 already used by the base class
   TString commonOutputName = AliAnalysisManager::GetCommonFileName();
   commonOutputName += ":PWGHF_D2H_AnalysisDstoK0sK";
   commonOutputName += suffixName;


   AliAnalysisDataContainer *outputCuts = mgr->CreateContainer("AnalysisCuts", TList::Class(),
                                          AliAnalysisManager::kOutputContainer, commonOutputName.Data());
   mgr->ConnectOutput(task, 1, outputCuts);

   AliAnalysisDataContainer *outputNorm = mgr->CreateContainer("NormalisationCounter", AliNormalizationCounter::Class(),
                                          AliAnalysisManager::kOutputContainer, commonOutputName.Data());
   mgr->ConnectOutput(task, 2, outputNorm);

   AliAnalysisDataContainer *outputSele = mgr->CreateContainer("DstoK0sSelections", TList::Class(),
                                          AliAnalysisManager::kOutputContainer, commonOutputName.Data());
   mgr->ConnectOutput(task, 3, outputSele);


   if (!fFillNtuple) {
      AliAnalysisDataContainer *outputCand = mgr->CreateContainer("DstoK0sVariablesCandidates", TList::Class(),
                                             AliAnalysisManager::kOutputContainer, commonOutputName.Data());
      mgr->ConnectOutput(task, 4, outputCand);

      AliAnalysisDataContainer *outputPIDs = mgr->CreateContainer("DstoK0sVariablesPID", TList::Class(),
                                             AliAnalysisManager::kOutputContainer, commonOutputName.Data());
      mgr->ConnectOutput(task, 5, outputPIDs);

      if (fReadMC) {
         AliAnalysisDataContainer *outputMC   = mgr->CreateContainer("DstoK0sMonteCarlo", TList::Class(),
                                                AliAnalysisManager::kOutputContainer, commonOutputName.Data());
         mgr->ConnectOutput(task, 6, outputMC);
      }
   } else {
      AliAnalysisDataContainer *outputTupl = mgr->CreateContainer("DstoK0sKTuple", TNtuple::Class(),
                                             AliAnalysisManager::kOutputContainer, commonOutputName.Data());
      mgr->ConnectOutput(task, 4, outputTupl);
   }


   return task;
}
