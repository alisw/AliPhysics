/**************************************************************************
 * Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//=========================================================================//
//             AliEbyE Analysis for Net-Particle Higher Moment study       //
//                           Nirbhay K. Behera                             //
//                           nbehera@cern.ch                               //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//                                                                         //
//                        (Last Modified 2018/08/27)                       //
//                 Dealing with Wide pT Window Modified to ESDs            //
//Some parts of the code are taken from J. Thaeder/ M. Weber NetParticle   //
//analysis task.                                                           //
//=========================================================================//
TString fileNameBase="AnalysisResults.root";

//Caution-> runName: LHC10h, LHC11h, LHC15o and LHC10hAMPT (only supported)
//v7: support only for proton. 

AliAnalysisTask *AddAliEbyEPidEfficiencyContamination(
						      TString runName = "LHC10h",
						      Bool_t isModeAOD = 0,
						      Int_t aodFilterBit = 768, 
						      Bool_t IsMC  = 0.,
						      Int_t pidtype = 1,//pidtype:charge-0, pion-1, kaon-2, proton-3
						      Double_t DCAxy = 2.,
						      Double_t DCAz = 2.,
						      Double_t Vx = 0.3,
						      Double_t Vy = 0.3,
						      Double_t Vz = 10.,
						      Bool_t IsRapCut = 0,
						      Bool_t IsTotalMom = 0,
						      Double_t Eta = 0.8,
						      Int_t TPCCrossRow = 80,
						      Double_t Chi2NDF = 4.,
						      const char* CentEstimator = "V0M",
						      Int_t pidstrategy = 0,
						      Float_t nSigmaITS = 2.5,
						      Float_t nSigmaTPC = 2.5,
						      Float_t nSigmaTOF = 2.5,
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


  task->SetRunPeriod(runName);
  
  if(isModeAOD) {
    task->SetIsAOD(isModeAOD);                       
    task->SetAODtrackCutBit(aodFilterBit);
  }

   task->SetIsMC(IsMC);
   task->SetDca( DCAxy,DCAz );
   task->SetVertexDiamond( Vx, Vy, Vz );
   task->SetIsRapidityCut( IsRapCut );
   task->SetUseTotalMomentumCut( IsTotalMom );
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
   task->SetPidStrategy(pidstrategy);
   task->SetNSigmaMaxITS(nSigmaITS);
   task->SetNSigmaMaxTPC(nSigmaTPC);
   task->SetNSigmaMaxTOF(nSigmaTOF);

   if( runName == "LHC10h") task->SelectCollisionCandidates(AliVEvent::kMB);
   else if ( runName == "LHC15o") task->SelectCollisionCandidates(AliVEvent::kINT7);
   
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

