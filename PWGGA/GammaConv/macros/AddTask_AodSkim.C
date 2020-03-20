// $Id: AddTaskSkim.C 4586 2013-01-16 15:32:16Z loizides $
AliConversionAodSkimTask *AddTask_AodSkim( const Double_t mine             = -1,
                                const Double_t mintrackpt       = -1,
                                const Double_t minconvpt        = -1,
                                const Double_t minconveta       = -999,
                                const Double_t maxconveta       =  999,
                                const Double_t minconvphi       = -999,
                                const Double_t maxconvphi       =  999,
                                const Bool_t doBothMinEandPt    = kFALSE,
                                const Bool_t doBothConvPtandAcc = kFALSE,
                                const Double_t ycutmc           = 0.9,
                                const char *gammabr             = 0,
                                const UInt_t trigsel            = AliVEvent::kAny,
                                const Bool_t doCopyTOF          = kTRUE,
                                const Bool_t doCopyTracklets    = kTRUE,
                                const Bool_t doCopyTracks       = kTRUE,
                                const Bool_t doCopyTrigger      = kTRUE,
                                const Bool_t doCopyPTrigger     = kTRUE,
                                const Bool_t doCopyCells        = kTRUE,
                                const Bool_t doCopyPCells       = kTRUE,
                                const Bool_t doCopyClusters     = kTRUE,
                                const Bool_t doCopyDiMuons      = kFALSE,
                                const Bool_t doCopyTrdTracks    = kFALSE,
                                const Bool_t doCopyCascades     = kFALSE,
                                const Bool_t doCopyV0s          = kFALSE,
                                const Bool_t doCopyMC           = kTRUE,
                                const Bool_t doQA               = kFALSE, // activate additional QA histograms
                                const char *name=0)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskAodSkim", "No analysis manager found.");
    return 0;
  }

  TString tname(name);
  if (name==0)
    tname = Form("aodskim_mine%.1f_trigsel%u",mine,trigsel);

  TString fname(tname);
  fname += ".root";

  cout << "Calling AddTaskAodSkim with " << endl;
  cout << "-> mine=" << mine << endl;
  cout << "-> trigsel=" << trigsel << endl;
  if (name)
    cout << "-> name=" << name << endl;
  else
    cout << "-> name=0" << endl;
  cout << "Writing skimmed aod tree to " << fname.Data() << endl;

  AliAODHandler *output = (AliAODHandler*)mgr->GetOutputEventHandler();
   if (output==0) {
    AliAODHandler *output = new AliAODHandler;
    output->SetOutputFileName(fname);
    output->SetFillAOD(1);
    output->SetFillAODforRun(1);
    output->SetFillExtension(0);
    output->SetTreeBuffSize(30*1024*1024);
    mgr->SetOutputEventHandler(output);
  }

  AliAODInputHandler *input = (AliAODInputHandler*)mgr->GetInputEventHandler();
  AliConversionAodSkimTask *task = new AliConversionAodSkimTask(tname);
  task->SetClusMinE(mine);
  task->SetTrackMinPt(mintrackpt);
  task->SetDoBothMinTrackAndClus(doBothMinEandPt);
  task->SetDoQA(doQA);

  task->SetConvMinPt(minconvpt);
  task->SetConvMinEta(minconveta);
  task->SetConvMaxEta(maxconveta);
  task->SetConvMinPhi(minconvphi);
  task->SetConvMaxPhi(maxconvphi);
  task->SetDoBothConvPtAndAcc(doBothConvPtandAcc);
  if (gammabr) {
    task->SetGammaBrName(gammabr);
    task->SetCopyConv(kTRUE);
    input->AddFriend((char*)"AliAODGammaConversion.root");
  } else {
    task->SetCopyConv(kFALSE);
  }
  task->SelectCollisionCandidates(trigsel);
  task->SetCopyHeader(kTRUE);
  task->SetCopyVZERO(kTRUE);
  task->SetCopyTZERO(kTRUE);
  task->SetCopyVertices(kTRUE);
  task->SetCopyZDC(kTRUE);
  if (ycutmc>0) {
    task->SetRemoveMcParts(kTRUE);
    task->SetCutMcY(ycutmc);
  } else {
    task->SetRemoveMcParts(kFALSE);
  }
  task->SetCopyTOF(doCopyTOF);
  task->SetCopyTracklets(doCopyTracklets);
  task->SetCopyTracks(doCopyTracks);
  task->SetCopyTrigger(doCopyTrigger);
  task->SetCopyPTrigger(doCopyPTrigger);
  task->SetCopyCells(doCopyCells);
  task->SetCopyPCells(doCopyPCells);
  task->SetCopyClusters(doCopyClusters);
  task->SetCopyDiMuons(doCopyDiMuons);
  task->SetCopyTrdTracks(doCopyTrdTracks);
  task->SetCopyCascades(doCopyCascades);
  task->SetCopyV0s(doCopyV0s);
  task->SetCopyMC(doCopyMC);
  task->SetCopyMCHeader(doCopyMC);

  cout << "Task configured with: " << task->Str() << endl;
  mgr->AddTask(task);

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer()  ;
  TString contName(tname);
  contName += "_histos";
  AliAnalysisDataContainer *coutput = mgr->CreateContainer( contName.Data(),
                                                            TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            Form("%s_histos.root", tname.Data()));
  mgr->ConnectInput  (task, 0,  cinput );
  mgr->ConnectOutput (task, 1, coutput );

  return task;
}
