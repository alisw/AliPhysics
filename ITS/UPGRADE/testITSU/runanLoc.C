TChain* loadChain(const char* inpData, const char* chName="esdTree");
const TObjArray* GetTrackingConditions(const char* recopar="ITS/Calib/RecoParam/Run0_999999999_v0_s0.root");

void runanLoc(const char* inpData="AliESDs.root", const char* outFile="tstPerfAn.root")
{
  // runs analysis task on local data
  //
  gROOT->ProcessLine(".x LoadLibs.C(1)");
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) mgr = new AliAnalysisManager("Test reco analysis");
  AliAnalysisDataContainer *cinp = (AliAnalysisDataContainer*)mgr->GetCommonInputContainer();
  //
  AliESDInputHandler *esdInputHandler = dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!esdInputHandler) esdInputHandler = new AliESDInputHandlerRP();
  mgr->SetInputEventHandler(esdInputHandler);
  //
  AliMCEventHandler* mcInputHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!mcInputHandler) mcInputHandler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcInputHandler);
  mcInputHandler->SetPreReadMode(AliMCEventHandler::kLmPreRead);
  //
  gROOT->ProcessLine(".L AliTaskITSUPerf.cxx++");
  AliTaskITSUPerf* task = new AliTaskITSUPerf("tstAnalysis");
  TString outFString = outFile;
  if (outFString.IsNull()) outFString = "recAnOut.root";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clist", TList::Class(),AliAnalysisManager::kOutputContainer,outFString.Data());
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0,  mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  task->SetUseMC(1);  // at the moment not irrelevant
  //
  // configure task
  task->SetEtaMin(-3);
  task->SetEtaMax(3);
  task->SetNPtBins(10);
  task->SetPtMax(2);
  //
  // optionally set reconstructability condition
  task->SetTrackingConditions(GetTrackingConditions());
  //
  TChain *chain = loadChain(inpData);
  if (!chain) return;
  mgr->InitAnalysis();
  mgr->StartAnalysis("local",chain);
  //
}

const TObjArray* GetTrackingConditions(const char* recopar)
{
  // load tracking conditions used in reco
  TString rpfstr = recopar;
  Bool_t fail = kFALSE;
  const TObjArray* trCondArr = 0;
  while(1) {
    TFile* rpf = 0;
    if (rpfstr.IsNull() || !(rpf=TFile::Open(rpfstr.Data()))) {fail = kTRUE; break;}
    AliCDBEntry* cdbe = (AliCDBEntry*)rpf->Get("AliCDBEntry");
    if (!cdbe) {fail = kTRUE; break;}
    TObjArray* recoParArr = (TObjArray*)cdbe->GetObject();
    if (!recoParArr) {fail = kTRUE; break;}
    AliITSURecoParam* recoPar = (AliITSURecoParam*) recoParArr->At(1);
    if (!recoPar) {fail = kTRUE; break;}
    trCondArr = recoPar->GetTrackingConditions();
    if (!trCondArr) {fail = kTRUE; break;}
    break;
  }
  if (fail) {printf("Tracking conditions must be extracted from valid recoparams, path is %s\n",recopar); exit(1);}
  //
  return trCondArr;
}

//________________________________________________________________________________
TChain* loadChain(const char* inpData, const char* chName)
{
  TChain* chain = new TChain(chName);
  //
  TString inpDtStr = inpData;
  if (inpDtStr.EndsWith(".root")) {
    if (!inpDtStr.Contains("/")) inpDtStr = Form("%s/%s",gSystem->pwd(),inpDtStr.Data());
    printf("Setting %s as an input\n",inpDtStr.Data());
    chain->AddFile(inpDtStr.Data());
  }
  else {
    //
    ifstream inpf(inpData);
    if (!inpf.good()) {
      printf("Failed on input filename %s\n",inpData);
      return kFALSE;
    }
    //
    TString flName;
    flName.ReadLine(inpf);
    while ( !flName.IsNull() ) {
      flName = flName.Strip(TString::kBoth,' ');
      if (flName.BeginsWith("//") || flName.BeginsWith("#")) {flName.ReadLine(inpf); continue;}
      flName = flName.Strip(TString::kBoth,',');
      flName = flName.Strip(TString::kBoth,'"');
      if (!flName.Contains("/")) flName = Form("%s/%s",gSystem->pwd(),flName);
      printf("Adding %s\n",flName.Data());
      chain->AddFile(flName.Data());
      flName.ReadLine(inpf);
    }
  }
  //
  int n = chain->GetEntries();
  if (n<1) {
    printf("Obtained chain is empty\n");
    return 0;
  }
  else printf("Opened %s chain with %d entries\n",chName,n);
  return chain;
}
