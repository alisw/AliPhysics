#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisTaskSEXicTopKpi.h"
#include "AliAnalysisManager.h"
#include "TFile.h"
#include "AliAnalysisDataContainer.h"
#include "TChain.h"
#endif

AliAnalysisTaskSEXicTopKpi *AddTaskSEXicTopKpi(Bool_t readMC=kFALSE,
					     Int_t system=0/*0=pp,1=PbPb*/,
					     Float_t minC=0, Float_t maxC=0,
					     TString finDirname="",
					       TString finname="",TString finObjname="D0toKpiCuts",TString strLcCutFile="",TString strLcCutObjName="",TString suffix="")
{
  //
  // AddTask for the AliAnalysisTaskSE for D0 candidates
  // invariant mass histogram and association with MC truth 
  // (using MC info in AOD) and cut variables distributions
  // C.Bianchin  chiara.bianchin@pd.infn.it
  //


  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskD0Distr", "No analysis manager to connect to.");
    return NULL;
  }   

  TString filename="",out1name="", out2name="", out3name="", out4name="" ,inname="";
  filename = AliAnalysisManager::GetCommonFileName();
  filename += ":PWG3_D2H_";

    filename+="XicpKpi";
    out1name="entriesHisto";
    //AliNormalizationCounter
    out2name="normalizationCounter";
    out3name="outputList";
    out4name="varTree";
    filename+=suffix.Data();
    
    inname="cinputmassD0_1";

    inname += finDirname.Data();
    inname += suffix.Data();
    out1name += suffix.Data();
    out2name += suffix.Data();
    out3name += suffix.Data();
    out4name += suffix.Data();

   //setting my cut values

    //cuts order
    //       printf("    |M-MD0| [GeV]    < %f\n",fD0toKpiCuts[0]);
    //     printf("    dca    [cm]  < %f\n",fD0toKpiCuts[1]);
    //     printf("    cosThetaStar     < %f\n",fD0toKpiCuts[2]);
    //     printf("    pTK     [GeV/c]    > %f\n",fD0toKpiCuts[3]);
    //     printf("    pTpi    [GeV/c]    > %f\n",fD0toKpiCuts[4]);
    //     printf("    |d0K|  [cm]  < %f\n",fD0toKpiCuts[5]);
    //     printf("    |d0pi| [cm]  < %f\n",fD0toKpiCuts[6]);
    //     printf("    d0d0  [cm^2] < %f\n",fD0toKpiCuts[7]);
    //     printf("    cosThetaPoint    > %f\n",fD0toKpiCuts[8]);

  TFile* filecuts;
  Printf("Opening file %s",finname.Data());
  filecuts=TFile::Open(finname.Data());
  if(!filecuts ||(filecuts&& !filecuts->IsOpen())){
    ::Fatal("AddTaskXictopKpi", "Input file not found : check your cut object");
  }

  
  AliRDHFCutsD0toKpi* RDHFD0toKpi=new AliRDHFCutsD0toKpi();

  RDHFD0toKpi = (AliRDHFCutsD0toKpi*)filecuts->Get(finObjname.Data());
  if(!RDHFD0toKpi){
    ::Fatal("AddTaskD0Mass", "Specific AliRDHFCuts not found");
    return NULL;
  }
  if(minC!=0 && maxC!=0) { //if centrality 0 and 0 leave the values in the cut object
      RDHFD0toKpi->SetMinCentrality(minC);
      RDHFD0toKpi->SetMaxCentrality(maxC);
    } 

  //  RDHFD0toKpi->SetName(Form("D0toKpiCuts%d",flag));

  TString centr="";
  if(minC!=0 && maxC!=0) centr = Form("%.0f%.0f",minC,maxC);
  else centr = Form("%.0f%.0f",RDHFD0toKpi->GetMinCentrality(),RDHFD0toKpi->GetMaxCentrality());
  out1name+=centr;
  inname+=centr;

  // Aanalysis task    
  TString taskname="XicPkPiMassAndDistrAnalysis";
  AliAnalysisTaskSEXicTopKpi *massXicTask = new AliAnalysisTaskSEXicTopKpi(taskname.Data(),RDHFD0toKpi);
  massXicTask->SetDebugLevel(-1);
  massXicTask->SetReadMC(readMC);
  massXicTask->SetSystem(system); //0=pp, 1=pPb, 2=PbPb

  // moved in RunAnalysis (i.e.: macro customization)
  //massXicTask->SetUseLcTrackFilteringCut(kTRUE);
  //massXicTask->SetMaxPtSPDkFirst(kTRUE,1.);
  //massXicTask->SetFillTree(2);  // moved in RunAnalysis (i.e.: macro customization)
  //massXicTask->SetMaxChi2Cut(1.5);


  if(!strLcCutFile.IsNull()){
    TFile *flc=TFile::Open(strLcCutFile.Data(),"READ");
    //printf("\n===== Lc cut file open: %s =====\n\n",strLcCutFile.Data());
    printf("\n===== Xic cut file open: %s =====\n\n",strLcCutFile.Data());
    if(!flc){
      Printf("Wrong file for Xic cuts");
      return 0x0;
    }

    //AliRDHFCutsLctopKpi* cutsLc=(AliRDHFCutsLctopKpi*)flc->Get(strLcCutObjName.Data());
    //if(!cutsLc){    
    AliRDHFCutsXictopKpi* cutsXic=(AliRDHFCutsXictopKpi*)flc->Get(strLcCutObjName.Data());
    if(!cutsXic){
      Printf("Wrong object name for Lc cuts");
      return 0x0;
    }
    //massXicTask->SetLcCuts(cutsLc);
    massXicTask->SetXicCuts(cutsXic);
  }
  if(system==0){
    massXicTask->SetRecalcOwnPrimVtx(kTRUE);
  }
  //  massXicTask->SetAODMismatchProtection(AODProtection);
  //  massD0Task->SetRejectSDDClusters(kTRUE);

  //   massD0Task->SetWriteVariableTree(kTRUE);

  mgr->AddTask(massXicTask);
  
  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinputmass = mgr->CreateContainer(inname,TChain::Class(), 
							  AliAnalysisManager::kInputContainer);

  AliAnalysisDataContainer *coutputmass1 = mgr->CreateContainer(out1name,TH1F::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //mass
  //  AliAnalysisDataContainer *coutputmassD04 = mgr->CreateContainer(out4name,AliRDHFCutsD0toKpi::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //cuts
  AliAnalysisDataContainer *coutputmass2 = mgr->CreateContainer(out2name,AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //counter
  AliAnalysisDataContainer *coutputmass3 = mgr->CreateContainer(out3name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //counter
  AliAnalysisDataContainer *coutputmass4 = mgr->CreateContainer(out4name,TTree::Class(),AliAnalysisManager::kOutputContainer, filename.Data());

  mgr->ConnectInput(massXicTask,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(massXicTask,1,coutputmass1);
  mgr->ConnectOutput(massXicTask,2,coutputmass2);
  mgr->ConnectOutput(massXicTask,3,coutputmass3);
  mgr->ConnectOutput(massXicTask,4,coutputmass4);

  return massXicTask;
}
