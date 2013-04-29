// $Id$

const Int_t  mode   = 0; // 0-->Local, 1-->Grid interactive
const Bool_t bIsAOD = kTRUE;
const Bool_t bIsMC  = kFALSE;

const Bool_t bIsPhysSel = kFALSE;
const Bool_t bIsCentSel = kFALSE;
const Bool_t bIsEvnPSel = kFALSE;
const Bool_t bIsRespPID = kFALSE;

const TString setData = "dataset.txt";
//=============================================================================

const Bool_t bIsEMCalAna = kFALSE;

const Bool_t   bFastOnly   = kTRUE;
const Bool_t   bHistoW     = kTRUE;
const Double_t dMinEclus   =  5.;
const Double_t dMinPtClus  =  5.;
const Double_t dVz         = 10.;
const Bool_t   bVzDiff     = kTRUE;
const Double_t dCentMin    = -1.;
const Double_t dCentMax    = -1.;
const Double_t dMinScaleCT = -1.;
const Double_t dMaxScaleCT = -1.;
//=============================================================================

void AnalysisTrainCorrJetsLocal()
{
  if (mode==1 && !TGrid::Connect("alien://")) {
    ::Error("AnalysisTrainCorrJetsLocal.C::AnalysisTrainCorrJetsLocal", "Can not connect to the Grid!");
    return;
  }
  if (LoadLibraries()) {
    ::Error("AnalysisTrainCorrJetsLocal.C::AnalysisTrainCorrJetsLocal", "Load libraries failed!");
    return;
  }
//=============================================================================

  TChain *chain = CreateAODFriendChain(setData);
  if (!chain) {
    ::Error("AnalysisTrainCorrJetsLocal.C::AnalysisTrainCorrJetsLocal", "Creating input chain failed!");
    return;
  }
//=============================================================================

  const UInt_t triggerMask = AliVEvent::kMB;
  AliAnalysisManager *mgr  = new AliAnalysisManager("AnalysisTrainCorrJetsLocal", "Analysis Train Jet Correlation Local");

  if (bIsAOD) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C");
    AliAODInputHandler *aodIH = AddAODHandler();
  } else {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C");
    AliESDInputHandler *esdIH = AddESDHandler();
//  esdIH->SetReadFriends(kFALSE);
  }
  if (bIsMC && !bIsAOD) {
    AliMCEventHandler *mctEH = new AliMCEventHandler();
    mcH->SetPreReadMode(AliMCEventHandler::kLmPreRead);
    mcH->SetReadTR(kTRUE);
    mgr->SetMCtruthEventHandler(mcH);
  }
//=============================================================================

  if (bIsPhysSel && !bIsAOD && !bIsEMCalAna) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask *taskPhysSel = AddTaskPhysicsSelection(bIsMC);
  }

  if (bIsCentSel && !bIsAOD) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
    if (bIsMC) taskCentrality->SetMCInput();
  }

  if (bIsEvnPSel) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskEventplane.C");
    AliEPSelectionTask *taskEventPlane = AddTaskEventplane();
    if (bIsMC) taskEventPlane->SetUseMCRP();
  }

  if (bIsRespPID) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTaskSE *taskRespPID = AddTaskPIDResponse(bIsMC);
  }
//=============================================================================

  if (bIsEMCalAna) {
    if (bIsEvnPSel) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPhysicsSelection.C");
      AliPhysicsSelectionTask *taskPhysSel = AddTaskEmcalPhysicsSelection(bFastOnly,   bHistoW, triggerMask,
                                                                          dMinEclus,   dMinPtClus,
                                                                          dVz,         bVzDiff,
                                                                          bIsCentSel ? dCentMin : -1.,
                                                                          bIsCentSel ? dCentMax : -1.,
                                                                          dMinScaleCT, dMaxScaleCT);
    }

    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalSetup.C");
    AliEmcalSetupTask *taskSetupEMCal = AddTaskEmcalSetup();
    taskSetupEMCal->SetGeoPath("$ALICE_ROOT/OADB/EMCAL");

    gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskJetPreparation.C");
    AddTaskJetPreparation(bIsAOD ? "AOD" : "ESD"); //TODO
  }
//=============================================================================

  gROOT->LoadMacro("AddTasksFlavorJet.C"); AddTasksFlavorJet();
  if (mgr->InitAnalysis()) { mgr->PrintStatus(); mgr->StartAnalysis("local",chain); }
  return;
}

//=============================================================================
TChain *CreateAODFriendChain(TString setData)
{ 
  if (!setData.EndsWith(".txt")) return 0;

  TChain *chain = new TChain("aodTree");
  TChain *cFrid = new TChain("aodTree");

  TString dataFile;
  ifstream dataList(setData.Data(), ios::in); 
  while (!dataList.eof()) {
    dataFile.ReadLine(dataList,kFALSE);
    if (!dataFile.EndsWith("AliAOD.root")) continue;
    if (!gSystem->AccessPathName(dataFile.Data())) chain->Add(dataFile.Data());

    dataFile.ReplaceAll("AliAOD.root","AliAOD.VertexingHF.root");
    if (!gSystem->AccessPathName(dataFile.Data())) cFrid->Add(dataFile.Data());
  } dataList.close();

  chain->AddFriend(cFrid);
  return chain;
}

//=============================================================================
Bool_t LoadLibraries()
{
  if (gSystem->Load("libTree")       <0) return kTRUE;
  if (gSystem->Load("libGeom")       <0) return kTRUE;
  if (gSystem->Load("libPhysics")    <0) return kTRUE;
  if (gSystem->Load("libVMC")        <0) return kTRUE;
  if (gSystem->Load("libMinuit")     <0) return kTRUE;
  if (gSystem->Load("libMinuit2")    <0) return kTRUE;

  if (gSystem->Load("libCore")       <0) return kTRUE;
  if (gSystem->Load("libXMLIO")      <0) return kTRUE;
  if (gSystem->Load("libXMLParser")  <0) return kTRUE;
  if (gSystem->Load("libProof")      <0) return kTRUE;
  if (gSystem->Load("libProofPlayer")<0) return kTRUE;
  if (gSystem->Load("libGui")        <0) return kTRUE;
//=============================================================================

  gSystem->AddIncludePath("-Wno-deprecated");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_ROOT/PWGHF/vertexingHF");
  gSystem->AddIncludePath("-I$ALICE_ROOT/EMCAL");
  gSystem->AddIncludePath("-I$ALICE_ROOT/JETAN");
  gSystem->AddIncludePath("-I$ALICE_ROOT/JETAN/fastjet");

  if (gSystem->Load("libSTEERBase")         <0) return kTRUE;
  if (gSystem->Load("libESD")               <0) return kTRUE;
  if (gSystem->Load("libAOD")               <0) return kTRUE;
  if (gSystem->Load("libANALYSIS")          <0) return kTRUE;
  if (gSystem->Load("libOADB")              <0) return kTRUE;
  if (gSystem->Load("libANALYSISalice")     <0) return kTRUE;
  if (gSystem->Load("libCORRFW")            <0) return kTRUE;

  if (gSystem->Load("libPWGTools")          <0) return kTRUE;
  if (gSystem->Load("libPWGflowBase")       <0) return kTRUE;
  if (gSystem->Load("libPWGflowTasks")      <0) return kTRUE;
  if (gSystem->Load("libPWGHFbase")         <0) return kTRUE;
  if (gSystem->Load("libPWGHFvertexingHF")  <0) return kTRUE;

  if (gSystem->Load("libSTAT")              <0) return kTRUE;
  if (gSystem->Load("libEMCALUtils")        <0) return kTRUE;
//if (gSystem->Load("libPHOSUtils")         <0) return kTRUE;

  if (gSystem->Load("libCDB")               <0) return kTRUE;
  if (gSystem->Load("libRAWDatabase")       <0) return kTRUE;
  if (gSystem->Load("libRAWDatarec")        <0) return kTRUE;
  if (gSystem->Load("libSTEER")             <0) return kTRUE;
  if (gSystem->Load("libITSbase")           <0) return kTRUE;
  if (gSystem->Load("libITSrec")            <0) return kTRUE;
  if (gSystem->Load("libTPCbase")           <0) return kTRUE;
  if (gSystem->Load("libTPCrec")            <0) return kTRUE;
  if (gSystem->Load("libTRDbase")           <0) return kTRUE;
  if (gSystem->Load("libTRDrec")            <0) return kTRUE;
  if (gSystem->Load("libTOFbase")           <0) return kTRUE;
//if (gSystem->Load("libTOFrec")            <0) return kTRUE;
  if (gSystem->Load("libHMPIDbase")         <0) return kTRUE;
  if (gSystem->Load("libEMCALraw")          <0) return kTRUE;
  if (gSystem->Load("libEMCALbase")         <0) return kTRUE;
  if (gSystem->Load("libEMCALrec")          <0) return kTRUE;
  if (gSystem->Load("libVZERObase")         <0) return kTRUE;
  if (gSystem->Load("libVZEROrec")          <0) return kTRUE;
  if (gSystem->Load("libTENDER")            <0) return kTRUE;
  if (gSystem->Load("libTENDERSupplies")    <0) return kTRUE;

  if (gSystem->Load("libCGAL")              <0) return kTRUE;
  if (gSystem->Load("libfastjet")           <0) return kTRUE;
  if (gSystem->Load("libsiscone")           <0) return kTRUE;
  if (gSystem->Load("libSISConePlugin")     <0) return kTRUE;

  if (gSystem->Load("libJETAN")             <0) return kTRUE;
  if (gSystem->Load("libFASTJETAN")         <0) return kTRUE;
  if (gSystem->Load("libPWGEMCAL")          <0) return kTRUE;
  if (gSystem->Load("libPWGGAEMCALTasks")   <0) return kTRUE;
  if (gSystem->Load("libPWGJEEMCALJetTasks")<0) return kTRUE;

  if (gROOT->LoadMacro("AliAnalysisTaskSEDmesonsFilterCJ.cxx+")  <0) return kTRUE;
  if (gROOT->LoadMacro("AliAnalysisTaskFlavorJetCorrelations.cxx+")<0) return kTRUE;
  return kFALSE;
}
