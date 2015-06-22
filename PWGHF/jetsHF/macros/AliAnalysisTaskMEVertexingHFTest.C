#ifndef __CINT__
#include <AliVEvent.h>
#endif

void AliAnalysisTaskMEVertexingHFTest()
{
Int_t num=0;

   if (gSystem->Load("libTree.so") < 0) {num++; return num;}
   if (gSystem->Load("libGeom.so") < 0) {num++; return num;}
   if (gSystem->Load("libVMC.so") < 0) {num++; return num;}
   if (gSystem->Load("libMinuit.so") < 0) {num++; return num;}
   if (gSystem->Load("libPhysics.so") < 0) {num++; return num;}
   if (gSystem->Load("libSTEERBase.so") < 0) {num++; return num;}
   if (gSystem->Load("libESD.so") < 0) {num++; return num;}
   if (gSystem->Load("libAOD.so") < 0) {num++; return num;}
   if (gSystem->Load("libANALYSIS.so") < 0) {num++; return num;}
   if (gSystem->Load("libOADB.so") < 0) {num++; return num;}
   if (gSystem->Load("libANALYSISalice.so") < 0) {num++; return num;}
   if (gSystem->Load("libEventMixing.so") < 0) {num++; return num;}


  //
  // Test macro for the AliAnalysisTaskME for heavy-flavour event mixing
  // r.romita@gsi.de
  //

  Bool_t useParFiles=kFALSE;
  
  //gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/LoadLibraries.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/vertexingHF/macros/LoadLibraries.C");
  LoadLibraries(useParFiles);

  // Local files 
  
  TChain* chain = new TChain("aodTree");
  Char_t fileName[100];
  //sprintf(fileName,"AliAODs.root");
  sprintf(fileName,"/Users/chitrasen/alice/PbPb2760_AODs/alice/data/2010/LHC10h/000139510/ESDs/pass2/AOD086/0001/AliAOD.root");//CJ
  chain->Add(fileName);
  //chain->Draw("fRefMult");
 
  // Create the analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager","My Manager");
  mgr->SetDebugLevel(10);
 
  // Input Handler
  //AliMultiEventInputHandler *inputHandler = new AliMultiEventInputHandler(4,1);
  AliMultiInputEventHandler *multiInputHandler = new AliMultiInputEventHandler();

  multiInputHandler->AddInputEventHandler(new AliAODInputHandler());
  mgr->SetInputEventHandler(multiInputHandler);

  const Int_t bufferSize = 1;
  const Int_t mixNum = 5;
  AliMixInputEventHandler *mixHandler = new AliMixInputEventHandler ( bufferSize, mixNum );
  //mixHandler->SetInputHandlerForMixing ( dynamic_cast<AliMultiInputEventHandler *> ( mgr->GetInputEventHandler() ) );
  mixHandler->SetInputHandlerForMixing(multiInputHandler);


  //AliEventPoolOTF* pool = new AliEventPoolOTF("event pool", "AOD");
  // apply selections
  //pool->SetMultiplicityBin(0, 100, 2);
  //pool->SetZVertexBinning(-20., 20., 2);
  //pool->Init();

  AliMixEventPool *evPool = new AliMixEventPool();
  
  //AliMixEventCutObj *centrality = new AliMixEventCutObj ( AliMixEventCutObj::kCentrality, 0, 100, 10, "V0M" );
 AliMixEventCutObj *multi = new AliMixEventCutObj ( AliMixEventCutObj::kMultiplicity, 2, 102, 10 );
 AliMixEventCutObj *zvertex = new AliMixEventCutObj ( AliMixEventCutObj::kZVertex, -10, 10, 1 );

 //evPool->AddCut(centrality);
 evPool->AddCut(multi);
 evPool->AddCut(zvertex);
 
  //set tag directory
  //Char_t tagDir[100];
  //sprintf(tagDir,".");
  //evPool->SetTagDirectory(tagDir);
  //mgr->SetInputEventHandler(inputHandler);
  //mgr->SetEventPool(evPool);
  //inputHandler->SetEventPool(evPool);

 // adds event pool (comment it and u will have default mixing)
 mixHandler->SetEventPool(evPool);

 mixHandler->SelectCollisionCandidates(AliVEvent::kAny);
 //mixHandler->SelectCollisionCandidates(AliVEvent::kAnyINT);
 //mixHandler->DoMixIfNotEnoughEvents(kFALSE);

 multiInputHandler->AddInputEventHandler(mixHandler);


 // Output 
  AliAODHandler *aodHandler = new AliAODHandler();
  aodHandler->SetOutputFileName("AliAOD.VertexingHF.root");
  aodHandler->SetCreateNonStandardAOD();
  mgr->SetOutputEventHandler(aodHandler);

  
 
  //gROOT->LoadMacro("AddTaskHFMixing.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/vertexingHF/macros/AddTaskHFMixing.C");//CJ
  //AliAnalysisTaskMEVertexingHF *hfTask = AddTaskHFMixing();

  AliAnalysisTaskMEVertexingHF *hfTask = new AliAnalysisTaskMEVertexingHF("mixing vertexing HF");
  mgr->AddTask(hfTask);
   
  //
  // Run the analysis
  //    
  printf("CHAIN HAS %d ENTRIES\n",(Int_t)chain->GetEntries());
  if(!mgr->InitAnalysis()) return;

  TString anSrc = "local"; TString anMode = "test"; Long64_t nEvents = 1e10; Long64_t nSkip = 0;
  mgr->PrintStatus();
  
  TStopwatch watch;
  watch.Start();
  //mgr->StartAnalysis("mix",chain, 1000);
  mgr->StartAnalysis("mix",chain);
  //mgr->StartAnalysis(anSrc.Data(), nEvents, nSkip);
  watch.Stop();
  watch.Print();
  delete mgr;

  return;
}
