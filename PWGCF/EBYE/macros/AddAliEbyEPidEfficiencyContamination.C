//=========================================================================//
//                                                                         //
//           Analysis AddTask for Net-particle higher moments study        //
//              Author: Deepika Rathee  || Satyajit Jena                   //
//                      drathee@cern.ch || sjena@cern.ch                   //
//                        Thu Jun 19 11:44:51 CEST 2015                    //
//                              Nirbhay  K. Behera                         //
//                         (Last Modified: 2017/04/07)                     //
//=========================================================================//
TString fileNameBase="AnalysisResults.root";

AliAnalysisTask *AddAliEbyEPidEfficiencyContamination(Bool_t isModeAOD = 0,
						      Int_t aodFilterBit = 768, 
						      Bool_t IsMC  = 0.,
						      Bool_t IsQA = 0,
						      Int_t pidtype = 1,//pidtype:charge-0, pion-1, kaon-2, proton-3
						      Double_t DCAxy = 2.,
						      Double_t DCAz = 2.,
						      Double_t Vx = 0.3,
						      Double_t Vy = 0.3,
						      Double_t Vz = 10.,
						      Double_t Eta = 0.8,
						      Int_t TPCCrossRow = 80,
						      Double_t Chi2NDF = 4.,
						      const char* CentEstimator = "V0M",
						      Float_t nSigmaITS = 2.,
						      Float_t nSigmaTPC = 2.,
						      Float_t nSigmaTOF = 3.,
						      Float_t nSigmaTPClow = 2.,
						      Float_t minPtTOF = 0.69,
						      Float_t minPtTPClow = 0.69,
						      TString taskname = "TestHM") {

  Double_t ptl, pth;
 
  //Set the track cuts here ---
  if( pidtype == 0){
    ptl = 0.18;
    pth = 2.1;
  }
  else{
    ptl = 0.35;
    pth = 1.6;
  }


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddAliEbyEPidEfficiencyContamination", "No analysis manager to connect to.");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddAliEbyEPidEfficiencyContamination", "This task requires an input event handler");
    return NULL;
  }
 
  const Char_t *pidname[4] = {"Nch","Npi","Nka","Npr"};
  
  const Char_t *ctsk = Form("%sNET%s",pidname[pidtype], taskname);
  
  AliEbyEPidEfficiencyContamination *task = new AliEbyEPidEfficiencyContamination(ctsk);
  
  if(isModeAOD) {
    task->SetIsAOD(isModeAOD);                       
    task->SetAODtrackCutBit(aodFilterBit);
  }

   task->SetIsMC(IsMC);
   task->RunQA(IsQA);

   task->SetDca( DCAxy,DCAz );
   task->SetVertexDiamond( Vx, Vy, Vz );
   task->SetKinematicsCuts( ptl, pth, Eta );
   if( pidtype == 0){
     task->SetNumberOfPtBins( 16 );
   }
   else{
     task->SetNumberOfPtBins( 19 );
   }
   task->SetTPCTrackQualityCuts( TPCCrossRow, Chi2NDF );
   task->SetCentralityEstimator( CentEstimator );
   
   //--------- PID -----
   task->SetPidType(pidtype);
   task->SetNSigmaMaxITS(nSigmaITS);
   task->SetNSigmaMaxTPC(nSigmaTPC);
   task->SetNSigmaMaxTOF(nSigmaTOF);
   task->SetNSigmaMaxTPClow(nSigmaTPClow);
   task->SetMinPtForTOFRequired(minPtTOF);
   task->SetMaxPtForTPClow(minPtTPClow);
   
   task->SelectCollisionCandidates(AliVEvent::kMB);
   mgr->AddTask(task);
 
   
   //TString basefilename = AliAnalysisManager::GetCommonFileName();
   //TString basefilename = Form("AnalysisResult.%s.root",taskname.Data());
   
   AliAnalysisDataContainer *coutt   = mgr->CreateContainer(Form("Net%s%s",pidname[pidtype], taskname.Data()),
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    fileNameBase.Data() );
   mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, coutt);
   
  return task;
}

