
// $Id: AddTaskDiJetCorrelations.C

AliAnalysisTaskDiJetCorrelations *AddTaskDiJetCorrelations(TString suffixName="",
							   Bool_t SEorME = kTRUE,
							   Bool_t ppOrPbPb = kTRUE,
							   Double_t pTrg1min = 12.0,
							   Double_t pTrg1max = 16.0,
							   Double_t pTrg2min = 5.0,
							   Double_t pTrg2max = 8.0,
                                                         //  Double_t alpha = TMath::Pi(),
                               Double_t alphaB2B = TMath::Pi(),
							   TString effLoc = "")
{
  
  Bool_t UseFbits = kTRUE;
  
  //____________________________________| Correlation class setting..
  AliAnalysisTaskDiJetCorrelations *dijetcorrelations = new AliAnalysisTaskDiJetCorrelations("");
  dijetcorrelations->SelectCollisionCandidates(AliVEvent::kMB);
  dijetcorrelations->SetSystem(ppOrPbPb); //PbPb = kTRUE
  dijetcorrelations->SetSEorME(SEorME); //kTRUE for mixed events
  if(SEorME)dijetcorrelations->SetMESettings(500, 25000, 8); //evt,track,minMixEvents
  dijetcorrelations->SetTrigger1PTValue(pTrg1min, pTrg1max); //GeV/c
  dijetcorrelations->SetTrigger2PTValue(pTrg2min, pTrg2max); //GeV/c
  dijetcorrelations->SetFilterBit(UseFbits);
  if(UseFbits)dijetcorrelations->SetFilterType(272);
  if(ppOrPbPb)dijetcorrelations->SetCentralityRange(0., 100); // 0-100%
  dijetcorrelations->SetDataType(kTRUE); //track Data/MC tracks=1 or MC Part=0?
  dijetcorrelations->SetVarCentBin(kTRUE); // kTRUE have some trouble ! FIX ME
  dijetcorrelations->SetVarPtBin(kTRUE); // kTRUE have some trouble ! FIX ME
 // dijetcorrelations->SetAlphaAngle(alphaB2B);
  if(effLoc!="")dijetcorrelations->SetEffCorrection(GetEfficiencyCorr(effLoc));
  dijetcorrelations->SetAlphaAngle(alphaB2B);
 
    

  // Create containers for input/output        
  TString finDirname         = "_DiJetMayCERN";
  TString inname             = "cinputDiJetCorrelations";
  TString outBasicname       = "coutputDiJetBasicPlots";
  TString outCorrname        = "coutputDiJetCorrHistos";

  finDirname += suffixName.Data();
  inname            +=   finDirname.Data();
  outBasicname      +=   finDirname.Data();
  outCorrname       +=   finDirname.Data();
  
  
  // Get the pointer to the existing analysis manager via the static access method.                                                     
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    cout<<"AddTaskDiJetCorrelations", "No analysis manager to connect to."<<endl;
  }
  
  mgr->AddTask(dijetcorrelations);
  
  //Input and Output Slots:
  AliAnalysisDataContainer *cinputDiJetCorrelations = mgr->CreateContainer(inname,TChain::Class(), AliAnalysisManager::kInputContainer);
  //TString outputfile = AliAnalysisManager::GetCommonFileName();
  //outputfile += ":PWGCF_Di_Jet_Corr";
  TString outputfile = AliAnalysisManager::GetCommonFileName();//"resultsDiJetCorrelationsT112to16T25to8Dec8.root";

  
  AliAnalysisDataContainer *coutputDiJetCorrelations1 = mgr->CreateContainer(outBasicname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutputDiJetCorrelations2 = mgr->CreateContainer(outCorrname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  
  
  mgr->ConnectInput(dijetcorrelations,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(dijetcorrelations,1,coutputDiJetCorrelations1);
  mgr->ConnectOutput(dijetcorrelations,2,coutputDiJetCorrelations2);
  
  return dijetcorrelations;
  
}


/*
//loads the efficiency correction
TH3F *GetEfficiencyCorr(TString effLoc){

  TFile* f = 0x0;
  f = TFile::Open(effLoc.Data());
  if(!f){
    Printf("%s%d no input data",(char*)__FILE__,__LINE__);
    return;
  }
  
  TH3F *tmp1 = (TH3F*)f->Get("hpTetaCentRec");
  if(!tmp1){
    Printf("%s%d Couldn't find hpTetaCentRec",(char*)__FILE__,__LINE__);
    return;
  }
  
  return (TH3F*)tmp1->Clone("fEffHist3D");
}
 */
