#if !defined (__CINT__) || defined (__CLING__)
#include "AliATDJetCorrDev.h"
#include "AliAnalysisManager.h"
#include "TFile.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#endif

AliATDJetCorrDev* AddTaskDJetCorrDev(
  AliATDJetCorrDev::ECandidateType cand = AliATDJetCorrDev::kDstartoKpipi,
  TString filename = "DStartoKpipiCuts.root",
  Bool_t theMCon = kFALSE,
  Bool_t reco = kTRUE /*must be true if theMCon is false*/,
  TString suffix = "",
  TString jetArrname = "",
  TString trackArrname = "tracks",
  TString candArrName = "",
  TString rhoname="",
  TString jetArrBkgname = "",
  TString trackArrBkgname = "",
  TString sbArrName = "",
  TString rhonameBkg="",
  Float_t R = 0.4,
  Float_t jptcut = 10.,
  const char *cutType = "TPC",
  Double_t percjetareacut = -1.,
  AliATDJetCorrDev::ECorrelationMethod CorrMethod = AliATDJetCorrDev::kConstituent,
  Bool_t isPPData = kTRUE
)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskSEDmesonsFilterCJ", "No analysis manager to connect to.");
    return NULL;
  } 

  Bool_t useStdC = kFALSE;
  TFile* filecuts=TFile::Open(filename);
  if(!filecuts || (filecuts && !filecuts->IsOpen())) {
    std::cout<<"Input file not found: use std cuts"<<std::endl;
    useStdC = kTRUE;
  }

  AliRDHFCuts *analysiscuts=0x0;
  switch (cand) {
  case 0 :
    if(useStdC) {
      analysiscuts = new AliRDHFCutsD0toKpi();
      analysiscuts->SetStandardCutsPP2010();
    } else
      analysiscuts = (AliRDHFCutsD0toKpi*)filecuts->Get("D0toKpiCuts");
    break;
  case 1 :
    if(useStdC) {
      analysiscuts = new AliRDHFCutsDStartoKpipi();
      analysiscuts->SetStandardCutsPP2010();
    } else
      analysiscuts = (AliRDHFCutsDStartoKpipi*)filecuts->Get("DStartoKpipiCuts");
    analysiscuts->SetName("DStartoKpipiCuts");
    break;
  }
  
  if (!analysiscuts) { // mm let's see if everything is ok
    printf("Specific AliRDHFCuts not found");
    return NULL;
  } 

  printf("CREATE TASK\n"); //CREATE THE TASK

  TString candname="DStar"; 
  if(cand==0)  candname="D0";
  TString sR = Form("R%.0f",R*10);

  // create the task
  TString taskCorrName="TaskFlavourJetCorrelations";
  taskCorrName+=candname;
  taskCorrName+=suffix;
  if(theMCon) taskCorrName+="MC";
  if(!reco)   taskCorrName+="gen";
  taskCorrName+=cutType;
  taskCorrName+=Form("PTj%.0f",jptcut);
  taskCorrName+=sR;

  AliATDJetCorrDev* taskCorr =  new AliATDJetCorrDev(taskCorrName.Data(), analysiscuts, cand);

  taskCorr->SetCorrelationMethod(CorrMethod);
  taskCorr->SetMC(theMCon);
  taskCorr->SetUseReco(reco);
  taskCorr->SetIsPPData(isPPData);
  
  AliParticleContainer *trackCont  = taskCorr->AddParticleContainer(trackArrname);
  
  AliJetContainer *jetCont = taskCorr->AddJetContainer(jetArrname,cutType,R);
  if(jetCont) {
     jetCont->ConnectParticleContainer(trackCont);
     jetCont->SetJetPtCut(jptcut);
     jetCont->SetPercAreaCut(percjetareacut);
     jetCont->SetRhoName(rhoname); 
  }
    
  if(!jetArrBkgname.IsNull())
  {
      AliParticleContainer *trackContBkg  = taskCorr->AddParticleContainer(trackArrBkgname);
     
      AliJetContainer *jetContBkg = taskCorr->AddJetContainer(jetArrBkgname,cutType,R);
      if(jetContBkg) {
          jetContBkg->ConnectParticleContainer(trackContBkg);
          jetContBkg->SetJetPtCut(jptcut);
          jetContBkg->SetPercAreaCut(percjetareacut);
          jetContBkg->SetRhoName(rhonameBkg);
      }
  }

  mgr->AddTask(taskCorr);

  if(theMCon)
  {
     suffix+="MC";
     if(reco) suffix+="rec";  
  }
  
  // Create and connect containers for input/output
  TString nameContainerF0="histos";
  TString nameContainerF1="cutfile";
  
  TString nameContainerFC2="Dcandidates";
  TString nameContainerFC3="DSBcandidates";

  nameContainerF0  += candname;
  nameContainerF1  += candname;
  nameContainerFC2 += candname;
  nameContainerFC3 += candname;

  nameContainerF0  += suffix;
  nameContainerF1  += suffix;
  nameContainerFC2 += suffix;
  nameContainerFC3 += suffix;
    
    AliAnalysisDataContainer* coutput1 = mgr->CreateContainer(nameContainerF0, TList::Class(), AliAnalysisManager::kOutputContainer,Form("%s:DmesonsForJetCorrelations",AliAnalysisManager::GetCommonFileName()));


    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(nameContainerF1,AliRDHFCuts::Class(), AliAnalysisManager::kOutputContainer,Form("%s:DmesonsForJetCorrelations", AliAnalysisManager::GetCommonFileName()));

    mgr->ConnectInput(taskCorr,0,mgr->GetCommonInputContainer());

    TObjArray * cnt = mgr->GetContainers();

    if(!candArrName.IsNull())
    {
        AliAnalysisDataContainer* coutputFC2 = reinterpret_cast<AliAnalysisDataContainer*>(cnt->FindObject(nameContainerFC2));
        if (!coutputFC2) {
            ::Error("AliATDJetCorrDev", "Could not find input container '%s'!", nameContainerFC2.Data());
        }
        taskCorr->SetUseCandArray(kTRUE);
        mgr->ConnectInput(taskCorr, 1, coutputFC2);
    }

    if(!sbArrName.IsNull())
    {
        AliAnalysisDataContainer* coutputFC3 = reinterpret_cast<AliAnalysisDataContainer*>(cnt->FindObject(nameContainerFC3));
        if (!coutputFC3) {
            ::Error("AliATDJetCorrDev", "Could not find input container '%s'!", nameContainerFC3.Data());
        }
        taskCorr->SetUseSBArray(kTRUE);
        mgr->ConnectInput(taskCorr, 2, coutputFC3);
    }

    mgr->ConnectOutput(taskCorr, 1, coutput1);
    mgr->ConnectOutput(taskCorr, 2, coutput2);

    Printf("Input and Output connected to the manager");
    return taskCorr;
}

