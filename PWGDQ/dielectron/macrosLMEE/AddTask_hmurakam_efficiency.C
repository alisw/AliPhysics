// ROOT6 modifications
#ifdef __CLING__
#include <AliAnalysisManager.h>
#include <AliAODInputHandler.h>
#include <AliDielectronVarCuts.h>

// Tell ROOT where to find AliPhysics headers
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGDQ/dielectron/macrosLMEE/Config_hmurakam_efficiency.C>

#endif

AliAnalysisTaskElectronEfficiencyV2* AddTask_hmurakam_efficiency(TString name = "name",
                                                                 Bool_t isAOD = kTRUE,
                                                                 Bool_t getFromAlien = kFALSE,
                                                                 TString configFile="Config_hmurakam_efficiency.C",
                                                                 Bool_t tofcor = kTRUE,
                                                                 TString year  = "",
                                                                 TString generator ="",
                                                                 TString outputFileName="LMEE.root",
                                                                 TString suffix="")

{

  std::cout << "########################################\nADDTASK of ANALYSIS started\n########################################" << std::endl;

  // Configuring Analysis Manager
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName = outputFileName; // create a subfolder in the file

  // Loading individual config file either local or from Alien

 // ROOT6 modifications
#if defined(__CLING__)
  // ROOT6-specific code here....
  printf("ROOT6\n");
  
#elif defined(__CINT__)
  // ROOT5-specific code here ...
  
  printf("ROOT5 !!!\n");

  TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
  //TString configBasePath= Form("%s/",gSystem->pwd());
  //Load updated macros from private ALIEN path
  if (getFromAlien //&&
      && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/h/hmurakam/PWGDQ/dielectron/macrosLMEE/%s .",configFile.Data()))))
    {
      configBasePath=Form("%s/",gSystem->pwd());
    }
  TString configFilePath(configBasePath+configFile);

  // Loading config and cutlib
  Bool_t err=kFALSE;
  err |= gROOT->LoadMacro(configFilePath.Data());
  if (err) { Error("AddTask_hmurakam_efficiency","Config(s) could not be loaded!"); return 0x0; }

  std::cout << "Configpath:  " << configFilePath << std::endl;

#endif

  // Creating an instance of the task
  AliAnalysisTaskElectronEfficiencyV2* task = new AliAnalysisTaskElectronEfficiencyV2(name.Data());

  // AddGenerators
  task->SetGeneratorMCSignalName(generator);
  task->SetGeneratorULSSignalName(generator);
  // Set TOF correction
  if(tofcor){
    SetTOFSigmaEleCorrection(task, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, year.Data());
  }

  // Event selection. Is the same for all the different cutsettings
  task->SetEnablePhysicsSelection(kTRUE);
  task->SetTriggerMask(triggerNames);
  task->SetEventFilter(SetupEventCuts(isAOD)); //returns eventCuts from Config.
  task->SetCentrality(centMin, centMax);

  // Set minimum and maximum values of generated tracks. Only used to save computing power.
  // Do not set here your analysis pt-cuts
  task->SetMinPtGen(minGenPt);
  task->SetMaxPtGen(maxGenPt);
  task->SetMinEtaGen(minGenEta);
  task->SetMaxEtaGen(maxGenEta);

  // Set minimum and maximum values of generated tracks. Only used to save computing power.
  task->SetKinematicCuts(minPtCut, maxPtCut, minEtaCut, maxEtaCut);

  // Set Binning
  if (usePtVector == true) {
    std::vector<double> ptBinsVec;
    for (unsigned int i = 0; i < nBinsPt+1; ++i){
      ptBinsVec.push_back(ptBins[i]);
    }
    task->SetPtBins(ptBinsVec);
  }
  else task->SetPtBinsLinear   (minPtBin,  maxPtBin, stepsPtBin);
  task->SetEtaBinsLinear  (minEtaBin, maxEtaBin, stepsEtaBin);
  task->SetPhiBinsLinear  (minPhiBin, maxPhiBin, stepsPhiBin);
  task->SetThetaBinsLinear(minThetaBin, maxThetaBin, stepsThetaBin);

  // mee
  const Int_t Nmee = 1391;
  Double_t mee[Nmee];  //  Double_t mee[Nmee] = {};
  for(Int_t j=0;j<1100 ;j++) mee[j] = 0.001 * (j-  0)  +  0.0;//from 0 to 1.1 GeV/c2, every 1 MeV/c2
  for(Int_t k=1100;k<Nmee;k++) mee[k] = 0.01  * (k-1100) +  1.1;//from 1.1 to 4 GeV/c2, evety 10 MeV/c2
  //  TVectorD *v_mee = new TVectorD(Nmee);
  //  for(Int_t k=0;k<Nmee;k++) (*v_mee)[k] = mee[k];
  //  task->SetMassBins(v_mee);
  std::vector<double> massBinsVec;
  for (Int_t l = 0; l < Nmee; ++l) massBinsVec.push_back(mee[l]);
  task->SetMassBins(massBinsVec);
  
  // ptee
  const Int_t Nptee = 47;
  Double_t ptee[Nptee] = {0.0,
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
    2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0,
    5.0, 6.0, 7.0, 8.0, 9.0, 10.0};

  //  TVectorD *v_ptee = new TVectorD(Nptee);
  //  for(Int_t m=0;m<Nptee;m++) (*v_ptee)[m] = ptee[m];
  //  task->SetPairPtBins(v_pTee);
  std::vector<double> pteeBinsVec;
  for (Int_t l = 0; l < Nptee; ++l) pteeBinsVec.push_back(ptee[l]);
  task->SetPairPtBins(pteeBinsVec);

  // Resolution File, If resoFilename = "" no correction is applied
  //  SetResolutionFile(year);
  std::string resoFilename = Form("%s.root",year.Data());
  std::string resoFilenameFromAlien = Form("/alice/cern.ch/user/h/hmurakam/PWGDQ/dielectron/resolution/%s.root",year.Data());
  //  printf(Form("%s.root\n",year.Data()));

  cout << resoFilename << endl;
  //  task->SetResolutionFile(resoFilename);
  //  task->SetResolutionFileFromAlien(resoFilenameFromAlien);
  task->SetResolutionFile(resoFilename,"/alice/cern.ch/user/h/hmurakam/PWGDQ/dielectron/resolution/" + resoFilename);
  //  task->SetResolutionFile(resoFilename,resoFilenameFromAlien);
  task->SetResolutionDeltaPtBinsLinear   (DeltaMomMin, DeltaMomMax, NbinsDeltaMom);
  task->SetResolutionRelPtBinsLinear   (RelMomMin, RelMomMax, NbinsRelMom);
  task->SetResolutionEtaBinsLinear  (DeltaEtaMin, DeltaEtaMax, NbinsDeltaEta);
  task->SetResolutionPhiBinsLinear  (DeltaPhiMin, DeltaPhiMax, NbinsDeltaPhi);
  task->SetResolutionThetaBinsLinear(DeltaThetaMin, DeltaThetaMax, NbinsDeltaTheta);

  // Set centrality correction. If resoFilename = "" no correction is applied
  //task->SetCentralityFile(centralityFilename);

  // Pairing related config
  task->SetDoPairing(DoPairing);
  task->SetULSandLS(DoULSLS);

  // Add MCSignals. Can be set to see differences of:
  // e.g. secondaries and primaries. or primaries from charm and resonances
  AddSingleLegMCSignal(task);
  AddPairMCSignal(task);

  // Adding cutsettings
  TObjArray*  arrNames=names.Tokenize(";");
  const Int_t nDie=arrNames->GetEntriesFast();

  printf("Add %d cuts\n",nDie);
  for (int iCut = 0; iCut < nDie; ++iCut){
    //TString cutDefinition(arrNames->At(iCut)->GetName());
    AliAnalysisFilter* filter = SetupTrackCutsAndSettings(iCut, isAOD);
    task->AddTrackCuts(filter);
  }

  TString outlistname = Form("efficiency%s",suffix.Data());
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(outlistname, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  return task;

}
