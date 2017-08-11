void Pi0SpectrumAOD12dGroup2(Bool_t recompile = true, Bool_t kTender=true, TString tag="./")
{
  /* $Id: Pi0SpectrumLHC11h.C 59580 2012-11-14 14:45:19Z hqvigsta $ */
  // Based on $ALICE_ROOT/PWGGA/PHOSTasks/PHOS_PbPb/macros/GRID/Pi0SpectrumLHC11h.C 
  
  TStopwatch timer;
  timer.Start();
  
  TStringToken libs("Core,Tree,Geom,VMC,Physics,Minuit,Gui,XMLParser,Minuit2,Proof,STEERBase,ESD,AOD,OADB,ANALYSIS,ANALYSISalice,CDB,RAWDatabase,STEER,CORRFW,PHOSUtils,PHOSbase,PHOSpi0Calib,PHOSrec,PHOSshuttle,PHOSsim",",");
  while( libs.NextToken() )
    gSystem->Load( Form("lib%s.so", libs.Data()) );
  
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  
  //load analysis framework
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libPWGGAPHOSTasks.so");
  
  //load tender libraries
  gSystem->Load("libTender.so");
  gSystem->Load("libTenderSupplies.so");

  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include "
			  "-I$ALICE_PHYSICS -I$ALICE_PHYSICS/include "
			  "-I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER "
			  "-I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD "
			  "-I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWG "
			  "-I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/ANALYSIS -I$ALICE_PHYSICS/PHOS -I$ALICE_PHYSICS "
			  "-I$ALICE_PHYSICS/PWGGA/PHOSTasks/PHOS_pp_pi0 "
			  "-I$ALICE_PHYSICS/PWGGA/PHOSTasks/PHOS_PbPb -g"); 
  
  // A task can be compiled dynamically with AClic
  gROOT->LoadMacro("./AliCaloTriggerSimulator.cxx+g");
  gROOT->LoadMacro("./AliAnalysisTaskPHOSTrigPi0.cxx+g");
  
  // Connect to alien
  TString token = gSystem->Getenv("GRID_TOKEN") ;
  TGrid::Connect("alien://");
  AliAnalysisGrid *alienHandler=0x0;
  gROOT->LoadMacro("CreateAlienHandler.C");
  AliAnalysisGrid *alienHandler = CreateAlienHandler();
  if (!alienHandler) return;
  
  // Make the analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("Manager", "Manager");
  mgr->SetGridHandler(alienHandler);
  
  AliAODInputHandler* aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);
  
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");  
  AliAnalysisTaskPIDResponse* taskPID = AddTaskPIDResponse(false,false,true,2);

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/PHOSTasks/PHOS_PbPb/AddAODPHOSTender.C");  
  AliPHOSTenderTask* tender = AddAODPHOSTender("","","",1,false);
  AliPHOSTenderSupply* supply = tender->GetPHOSTenderSupply();
  supply->ForceUsingBadMap("alien:///alice/cern.ch/user/s/syano/OADB/syano_LHC12d_PHOS_BadMap.root");

  // Debug level
  mgr->SetDebugLevel(2);
  // Add my task

  gROOT->LoadMacro("./AddTaskPHOSTrigPi0.C");
  
  AliAnalysisTaskPHOSTrigPi0* task1 = AddTaskPHOSTrigPi0("kINT7","","LHC12d",AliVEvent::kINT7);
  task1->SetMC(false);
  task1->SetAnalysisTriggerClass("CINT7");
  task1->SetModuleCondition(0,true);
  task1->SetModuleCondition(1,false);
  task1->SetModuleCondition(2,true);
  
  task1->SetTRUCondition(0,0,true);
  task1->SetTRUCondition(0,1,false);
  task1->SetTRUCondition(0,2,true);
  task1->SetTRUCondition(0,3,false);
  task1->SetTRUCondition(0,4,false);
  task1->SetTRUCondition(0,5,true);
  task1->SetTRUCondition(0,6,true);
  task1->SetTRUCondition(0,7,true);

  task1->SetTRUCondition(2,0,true);
  task1->SetTRUCondition(2,1,true);
  task1->SetTRUCondition(2,2,true);
  task1->SetTRUCondition(2,3,false);
  task1->SetTRUCondition(2,4,false);
  task1->SetTRUCondition(2,5,true);
  task1->SetTRUCondition(2,6,false);
  task1->SetTRUCondition(2,7,true);

  AliAnalysisTaskPHOSTrigPi0* task2 = AddTaskPHOSTrigPi0("kPHI7","","LHC12d",AliVEvent::kPHI7);
  task2->SetMC(false);
  task2->SetAnalysisTriggerClass("CPHI7");
  task2->SetModuleCondition(0,true);
  task2->SetModuleCondition(1,false);
  task2->SetModuleCondition(2,true);

  task2->SetTRUCondition(0,0,true);
  task2->SetTRUCondition(0,1,false);
  task2->SetTRUCondition(0,2,true);
  task2->SetTRUCondition(0,3,false);
  task2->SetTRUCondition(0,4,false);
  task2->SetTRUCondition(0,5,true);
  task2->SetTRUCondition(0,6,true);
  task2->SetTRUCondition(0,7,true);

  task2->SetTRUCondition(2,0,true);
  task2->SetTRUCondition(2,1,true);
  task2->SetTRUCondition(2,2,true);
  task2->SetTRUCondition(2,3,false);
  task2->SetTRUCondition(2,4,false);
  task2->SetTRUCondition(2,5,true);
  task2->SetTRUCondition(2,6,false);
  task2->SetTRUCondition(2,7,true);
  
  TFile* efTiming = TFile::Open("alien:///alice/cern.ch/user/s/syano/OADB/timing_calib_map.root");
  TH2F* fHistMapLGT[3];
  TH2F* fHistMapHGT[3];
  fHistMapLGT[0] = (TH2F*)efTiming->Get("fHistMapLGTM1")->Clone("fHistMapLGTM1");
  fHistMapLGT[1] = (TH2F*)efTiming->Get("fHistMapLGTM2")->Clone("fHistMapLGTM2");
  fHistMapLGT[2] = (TH2F*)efTiming->Get("fHistMapLGTM3")->Clone("fHistMapLGTM3");
  fHistMapHGT[0] = (TH2F*)efTiming->Get("fHistMapHGTM1")->Clone("fHistMapHGTM1");
  fHistMapHGT[1] = (TH2F*)efTiming->Get("fHistMapHGTM2")->Clone("fHistMapHGTM2");
  fHistMapHGT[2] = (TH2F*)efTiming->Get("fHistMapHGTM3")->Clone("fHistMapHGTM3");
  task1->SetTOFCalibration(0,fHistMapLGT[0],fHistMapHGT[0]);
  task1->SetTOFCalibration(1,fHistMapLGT[1],fHistMapHGT[1]);
  task1->SetTOFCalibration(2,fHistMapLGT[2],fHistMapHGT[2]);
  task2->SetTOFCalibration(0,fHistMapLGT[0],fHistMapHGT[0]);
  task2->SetTOFCalibration(1,fHistMapLGT[1],fHistMapHGT[1]);
  task2->SetTOFCalibration(2,fHistMapLGT[2],fHistMapHGT[2]);

  TFile *fBadMap = TFile::Open("alien:///alice/cern.ch/user/s/syano/OADB/syano_LHC12d_PHOS_BadMap.root");
  TH2I * fMap[3];
  fMap[0] = (TH2I*) fBadMap->Get("PHOS_BadMap_mod1") ;
  fMap[1] = (TH2I*) fBadMap->Get("PHOS_BadMap_mod2") ;
  fMap[2] = (TH2I*) fBadMap->Get("PHOS_BadMap_mod3") ;
  task1->SetBadCell(0,fMap[0]);
  task1->SetBadCell(1,fMap[1]);
  task1->SetBadCell(2,fMap[2]);
  task2->SetBadCell(0,fMap[0]);
  task2->SetBadCell(1,fMap[1]);
  task2->SetBadCell(2,fMap[2]);

  TFile *fBadTileMap = TFile::Open("alien:///alice/cern.ch/user/s/syano/OADB/syano_LHC12d_PHOS_BadTileMaps.v2.root");
  TH2I * fTile[3];
  fTile[0] = (TH2I*) fBadTileMap->Get("hBadTileMapM1") ;
  fTile[1] = (TH2I*) fBadTileMap->Get("hBadTileMapM2") ;
  fTile[2] = (TH2I*) fBadTileMap->Get("hBadTileMapM3") ;

  TFile* multi_mean_corr_ef = TFile::Open("alien:///alice/cern.ch/user/s/syano/OADB/MultiCorrection.root");
  TH1F* pf_fHistSPDTrackletsUnitRap_mean  = (TH1F*)multi_mean_corr_ef->Get("prof_fHistSPDTrackletsUnitRap");
  TH1F* pf_fHistSPDTrackletsRapGap1_mean  = (TH1F*)multi_mean_corr_ef->Get("prof_fHistSPDTrackletsGapRap1");
  TH1F* pf_fHistSPDTrackletsRapGap2_mean  = (TH1F*)multi_mean_corr_ef->Get("prof_fHistSPDTrackletsGapRap2");
  TFile* multi_corr_ef = TFile::Open("alien:///alice/cern.ch/user/s/syano/OADB/MultiCorrection_LHC12dG2.root");
  TH1F* pf_fHistSPDTrackletsUnitRap  = (TH1F*)multi_corr_ef->Get("prof_fHistSPDTrackletsUnitRap");
  TH1F* pf_fHistSPDTrackletsRapGap1  = (TH1F*)multi_corr_ef->Get("prof_fHistSPDTrackletsGapRap1");
  TH1F* pf_fHistSPDTrackletsRapGap2  = (TH1F*)multi_corr_ef->Get("prof_fHistSPDTrackletsGapRap2");
  TF1*  fit_fHistSPDTrackletsUnitRap = (TF1*)multi_corr_ef->Get("fit_fHistSPDTrackletsUnitRap");
  TF1*  fit_fHistSPDTrackletsRapGap1 = (TF1*)multi_corr_ef->Get("fit_fHistSPDTrackletsGapRap1");
  TF1*  fit_fHistSPDTrackletsRapGap2 = (TF1*)multi_corr_ef->Get("fit_fHistSPDTrackletsGapRap2");
  Double_t mean_unit = pf_fHistSPDTrackletsUnitRap_mean->GetBinContent(1);
  Double_t mean_gap1 = pf_fHistSPDTrackletsRapGap1_mean->GetBinContent(1);
  Double_t mean_gap2 = pf_fHistSPDTrackletsRapGap2_mean->GetBinContent(1);

  TH1F* pf_fHistV0ACMulti_mean  = (TH1F*)multi_mean_corr_ef->Get("prof_fHistV0ACMulti");
  TH1F* pf_fHistV0ACMulti  = (TH1F*)multi_corr_ef->Get("prof_fHistV0ACMulti");
  TF1*  fit_fHistV0ACMulti = (TF1*)multi_corr_ef->Get("fit_fHistV0ACMulti");
  Double_t meanV0AC = pf_fHistV0ACMulti_mean->GetBinContent(1);

  task1->SetSPDTrackletMultiCorrections(pf_fHistSPDTrackletsUnitRap,pf_fHistSPDTrackletsRapGap1,pf_fHistSPDTrackletsRapGap2,mean_unit,mean_gap1,mean_gap2);
  task1->SetV0MultiCorrections(pf_fHistV0ACMulti,meanV0AC);
  task2->SetSPDTrackletMultiCorrections(pf_fHistSPDTrackletsUnitRap,pf_fHistSPDTrackletsRapGap1,pf_fHistSPDTrackletsRapGap2,mean_unit,mean_gap1,mean_gap2);
  task2->SetV0MultiCorrections(pf_fHistV0ACMulti,meanV0AC);

  task1->SetMeanMultiForMultiBin(7.867,3.235,6.626,77.9);
  task2->SetMeanMultiForMultiBin(7.867,3.235,6.626,77.9);

  mgr->InitAnalysis();
  mgr->SetDebugLevel(2);
  mgr->PrintStatus();
  mgr->StartAnalysis("grid");
  
  timer.Stop();
  timer.Print();
}
