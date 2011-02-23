/* $Id:  $ */

#ifndef __runCaloAODC__
#define __runCaloAODC__
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TFile.h>
#include <TROOT.h>
#include <TEnv.h>
#include <TProof.h>
#include <TChain.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TTimeStamp.h>
#include <TSystem.h>
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEMCALClusterize.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliBackgroundSelection.h"
#include "AliCentralitySelectionTask.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliMCEventHandler.h"
#include "AliPhysicsSelection.h"
#include "AliPhysicsSelectionTask.h"
#include "EMCAL/AliEMCALRecParam.h"
#include "aodtools.C"
#endif

class MyPhysicsSelectionTask : public AliPhysicsSelectionTask {
 public:
  MyPhysicsSelectionTask() : AliPhysicsSelectionTask() {} 
  MyPhysicsSelectionTask(const char *opt) : AliPhysicsSelectionTask(opt) {}
  void Terminate(Option_t*);
  ClassDef(MyPhysicsSelectionTask, 1); // My physics selection task
};

class MyCentralitySelectionTask : public AliCentralitySelectionTask {
 public:
  MyCentralitySelectionTask() : AliCentralitySelectionTask() {} 
  MyCentralitySelectionTask(const char *opt) : AliCentralitySelectionTask(opt) {}
  void UserExec(Option_t *opt) {LoadBranches(); return AliCentralitySelectionTask::UserExec(opt); }
  ClassDef(MyCentralitySelectionTask, 1); // My centrality selection task
};

class AODFillTask : public AliAnalysisTaskSE {
 public:
  AODFillTask(const char *opt=0);
  void UserExec(Option_t *opt);
  void UserCreateOutputObjects();
 protected:
  AliAODEvent *fEvent;    // aod event
  TTree       *fTree;     // tree
/*
  Bool_t       fDoZDC;    // do zdc
  Bool_t       fDoV0;     // do vzero
  Bool_t       fDoT0;     // do tzero
  Bool_t       fDoTPCv;   // do tpc vertex
  Bool_t       fDoSPDv;   // do spd vertex
  Bool_t       fDoPriv;   // do primary vertex
  Bool_t       fDoEmCs;   // do emcal cells
  Bool_t       fDoPCs;    // do phos cells
  Bool_t       fDoEmT;    // do emcal trigger
  Bool_t       fDoPT;     // do phos trigger
  Bool_t       fDoTracks; // do tracks
*/
  ClassDef(AODFillTask, 1); // AOD fill task
};

void addAOD(const AODProdParameters &pars);
void addCSel(const AODProdParameters &pars);
void addESel(const AODProdParameters &pars);
void addPSel(const AODProdParameters &pars);
void runCaloAOD(const AODProdParameters &params);

//_________________________________________________________________________________________________
void addAOD(const AODProdParameters &pars)
{
  // Add AOD task.

//  Bool_t doPsel = pars.doPS;
  // if (!doPsel)
  //  return;

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("addAOD", "No analysis manager to connect to");
    return;
  }    
  if (!mgr->GetInputEventHandler()) {
    ::Error("addAOD", "This task requires an input event handler");
    return;
  }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType();
  if (inputDataType != "ESD") {
    ::Error("addAOD", "This task works only on ESD analysis");
    return;
  }

#if 0
  MyEsdTrimTask *ttask = new MyEsdTrimTask("ETrim");
  ttask->SetBranches(pars.branches);
  mgr->AddTask(ttask);

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(ttask,0,cinput);
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("ESDTree",
                                                           TTree::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           "AliESDs.root");
  mgr->ConnectInput(ttask,0,cinput);
  mgr->ConnectOutput(ttask,1,coutput);
#endif
}   

//_________________________________________________________________________________________________
void addCSel(const AODProdParameters &pars)
{
  // Add centrality task.

  Bool_t doCsel = pars.doCS;
  if (!doCsel)
    return;

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("addCSel", "No analysis manager to connect to");
    return;
  }    
  if (!mgr->GetInputEventHandler()) {
    ::Error("addCSel", "This task requires an input event handler");
    return;
  }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType();
  if (inputDataType != "ESD") {
    ::Error("addCSel", "This task works only on ESD analysis");
    return;
  }

  AliCentralitySelectionTask *centralityTask = new MyCentralitySelectionTask("CentralitySelection");
  if (pars.doPS)
    centralityTask->SelectCollisionCandidates(AliVEvent::kMB);
  centralityTask->SetBranches(pars.branches);
  mgr->AddTask(centralityTask);
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("CentralityStat",
                                                           TList::Class(), 
                                                           AliAnalysisManager::kOutputContainer,
                                                           "EventStat.root");
  mgr->ConnectInput(centralityTask,0,cinput);
  mgr->ConnectOutput(centralityTask,1,coutput);
}   

//_________________________________________________________________________________________________
void addESel(const AODProdParameters &pars)
{
  //Add ESD filter task.

  Bool_t doAOD = pars.doAOD;
  if (!doAOD)
    return;

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("addESel", "No analysis manager to connect to");
    return;
  }    
  if (!mgr->GetInputEventHandler()) {
    ::Error("addESel", "This task requires an input event handler");
    return;
  }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType();
  if (inputDataType != "ESD") {
    ::Error("addESel", "This task works only on ESD analysis");
    return;
  }
  AliAODHandler *outh = dynamic_cast<AliAODHandler*>(mgr->GetOutputEventHandler());
  if (!outh) {
    ::Error("addESel", "Output handler does not handle AOD");
    return;
  }

  AliAODEvent *event = outh->GetAOD();
  if (event) {
    AliAODHeader *header = new AliAODHeader();
    header->SetName("header");
    event->AddObject(header);
    TClonesArray *vertices = new TClonesArray("AliAODVertex", 0);
    vertices->SetName("vertices");
    event->AddObject(vertices);
  }

  AliAnalysisTaskESDfilter *etask = new AliAnalysisTaskESDfilter("ESDfilter");
  if (pars.doPS)
    etask->SelectCollisionCandidates(AliVEvent::kMB);
  etask->SetBranches(pars.branches);
  etask->DisableVZERO();
  etask->DisableCascades();
  etask->DisableV0s();
  etask->DisableKinks();
  etask->DisableTracks();
  etask->DisablePmdClusters();
  etask->DisableCaloClusters();
  //etask->DisableCells();
  etask->DisableTracklets();
  mgr->AddTask(etask);
  mgr->ConnectInput(etask,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(etask,0,mgr->GetCommonOutputContainer());
}

//_________________________________________________________________________________________________
void addPSel(const AODProdParameters &pars)
{
  // Add physics selection task.

  Bool_t doPsel = pars.doPS;
  if (!doPsel)
    return;

  Bool_t isMC      = pars.doMC;
  Bool_t rejectBG  = 1;
  Bool_t computeBG = 0;

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("addPSel", "No analysis manager to connect to.");
    return;
  }    
  if (!mgr->GetInputEventHandler()) {
    ::Error("addPSel", "This task requires an input event handler");
    return;
  }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType();
  if (inputDataType != "ESD") {
    ::Error("addPSel", "This task works only on ESD analysis");
    return;
  }

  AliPhysicsSelectionTask *pseltask = new MyPhysicsSelectionTask("PS");
  mgr->AddTask(pseltask);
  
  AliPhysicsSelection *physSel = pseltask->GetPhysicsSelection();
  if (rejectBG) 
    physSel->AddBackgroundIdentification(new AliBackgroundSelection());
  if (computeBG)
    physSel->SetComputeBG(computeBG);
  if (isMC)      
    physSel->SetAnalyzeMC();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(pseltask,0,cinput);
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("cstatsout",
                                                           TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           "EventStat.root");
  mgr->ConnectInput(pseltask,0,cinput);
  mgr->ConnectOutput(pseltask,1,coutput);
}   

//_________________________________________________________________________________________________
void runCaloAOD(const AODProdParameters &params)
{
  // Run calo AOD.

  AliLog::SetGlobalLogLevel(AliLog::kError);

  TChain *chain = CreateChain(params);
  if (!chain) {
    ::Error("runCaloAOD", "Could not create chain");
    return;
  }

  AliAnalysisManager *mgr  = new AliAnalysisManager("Manager", "Manager");

  if (params.doAOD) {
    AliAODHandler* aodoutHandler   = new AliAODHandler();
    aodoutHandler->SetOutputFileName("AliAOD.root");
    //aodoutHandler->SetCreateNonStandardAOD();
    aodoutHandler->SetFillAOD(1);
    //aodoutHandler->Init("");
    mgr->SetOutputEventHandler(aodoutHandler);
  }

  if (params.runmode == "esd") {
    AliESDInputHandler *esdHandler = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdHandler);
    esdHandler->SetReadFriends(kFALSE);
  } else if (params.runmode == "aod") {
    ::Error("runCaloAOD", "AOD mode (at least for now) not supported");
    return;
    AliAODInputHandler *aodHandler = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodHandler);
  }

  if (params.doMC) {
    AliMCEventHandler *mcHandler = new AliMCEventHandler();
    mcHandler->SetPreReadMode(AliMCEventHandler::kLmPreRead);
    mgr->SetMCtruthEventHandler(mcHandler);
  }

  if (params.runtype == "local") {
    mgr->SetUseProgressBar(params.dlevel==0);
    mgr->SetDebugLevel(params.dlevel);
  } else {
    mgr->SetDebugLevel(0);
  }
  if (params.runtype == "proof") {
    mgr->SetAutoBranchLoading(1);
  } else {
    mgr->SetAutoBranchLoading(0);
  }
  if (params.runtype == "grid") {
    mgr->AddStatisticsTask();
  }

  addPSel(params);
  addCSel(params);
  addESel(params);
  //addTrim(params);

  mgr->InitAnalysis();
  mgr->PrintStatus();
  mgr->StartAnalysis(params.runtype,chain);
}

//_________________________________________________________________________________________________
void MyPhysicsSelectionTask::Terminate(Option_t* /*opt*/)
{
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {
    AliError("fOutput not available");
    return;
  }

  AliAnalysisDataSlot *oslot = GetOutputSlot(1);
  if (!oslot)
    return;

  AliAnalysisDataContainer *ocont = oslot->GetContainer();
  if (!ocont)
    return;

  TFile *file = OpenFile(1);
  if (!file)
    return;

  TDirectory::TContext context(file); 
  fPhysicsSelection->Print();
  fPhysicsSelection->SaveHistograms(Form("%sHists",ocont->GetName()));

  AliInfo(Form("Writing result to %s",file->GetName()));
}

//_________________________________________________________________________________________________
AODFillTask::AODFillTask(const char *opt) :
  AliAnalysisTaskSE(opt), fEvent(0), fTree(0)
{
  // Constructor.
  if (!opt)
    return;

  DefineOutput(1, TTree::Class());
}

//_________________________________________________________________________________________________
void AODFillTask::UserExec(Option_t *opt) 
{
  AliESDEvent *esdin = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esdin)
    return;

  fEvent->Reset();

#if 0
  TList* objs = fEvent->GetList();
  AliESDHeader *header = dynamic_cast<AliESDHeader*>(objs->FindObject("AliESDHeader"));
  if (header) {
    *header = *esdin->GetHeader();
  }
  AliESDRun *run = dynamic_cast<AliESDRun*>(objs->FindObject("AliESDRun"));
  if (run) {
    *run = *esdin->GetESDRun();
  }
  AliESDZDC *zdc = dynamic_cast<AliESDZDC*>(objs->FindObject("AliESDZDC"));
  if (zdc) {
    *zdc = *esdin->GetESDZDC();
  }
  AliESDVZERO *v0 = dynamic_cast<AliESDVZERO*>(objs->FindObject("AliESDVZERO"));
  if (v0) {
    *v0 = *esdin->GetVZEROData();
  }
  AliESDTZERO *t0 = dynamic_cast<AliESDTZERO*>(objs->FindObject("AliESDTZERO"));
  if (t0) {
    *t0 = *esdin->GetESDTZERO();
  }
  AliESDVertex *tpcv = dynamic_cast<AliESDVertex*>(objs->FindObject("TPCVertex"));
  if (tpcv) {
    *tpcv = *esdin->GetPrimaryVertexTPC();
  }
  AliESDVertex *spdv = dynamic_cast<AliESDVertex*>(objs->FindObject("SPDVertex"));
  if (spdv) {
    *spdv = *esdin->GetPrimaryVertexSPD();
  }
  AliESDVertex *priv = dynamic_cast<AliESDVertex*>(objs->FindObject("PrimaryVertex"));
  if (priv) {
    *priv = *esdin->GetPrimaryVertexTracks();
  }
  AliESDCaloCells *ecells = dynamic_cast<AliESDCaloCells*>(objs->FindObject("EMCALCells"));
  if (ecells) {
    *ecells = *esdin->GetEMCALCells();
  }
  AliESDCaloCells *pcells = dynamic_cast<AliESDCaloCells*>(objs->FindObject("PHOSCells"));
  if (pcells) {
    *pcells = *esdin->GetPHOSCells();
  }
  AliESDCaloTrigger *etrig = dynamic_cast<AliESDCaloTrigger*>(objs->FindObject("EMCALTrigger"));
  if (etrig) {
    *etrig = *esdin->GetCaloTrigger("EMCAL");
  }
  AliESDCaloTrigger *ptrig = dynamic_cast<AliESDCaloTrigger*>(objs->FindObject("PHOSTrigger"));
  if (ptrig) {
    *ptrig = *esdin->GetCaloTrigger("PHOS");
  }

  if (0) {
    AliESDtrackCuts *cuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    cuts->SetMinNClustersTPC(50);
    cuts->SetPtRange(0.2);
    cuts->SetEtaRange(-1.5,1.5);

    const Int_t Ntracks = esdin->GetNumberOfTracks();
    Int_t nacc = 0;
    for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
      AliESDtrack *track = esdin->GetTrack(iTracks);
      if (!track)
        continue;
      if (!cuts->IsSelected(track))
        continue;
      fEvent->AddTrack(track);
      ++nacc;
    }
    delete cuts;
    printf("selected %d out of %d \n", nacc, Ntracks);
  }
#endif
  fTree->Fill();
}

//_________________________________________________________________________________________________
void AODFillTask::UserCreateOutputObjects() 
{
  // Create output objects.

  fTree = new TTree("esdTree", "Tree with ESD objects");
  fEvent = new AliAODEvent;
#if 0
  fEvent->AddObject(new AliESDHeader());
  fEvent->AddObject(new AliESDRun());
  if (fDoZDC) 
    fEvent->AddObject(new AliESDZDC());
  if (fDoV0)
    fEvent->AddObject(new AliESDVZERO());
  if (fDoT0)
    fEvent->AddObject(new AliESDTZERO());
  if (fDoTPCv) {
    AliESDVertex *tpcv = new AliESDVertex();
    tpcv->SetName("TPCVertex");
    fEvent->AddObject(tpcv);
  }
  if (fDoSPDv) {
    AliESDVertex *spdv = new AliESDVertex();
    spdv->SetName("SPDVertex");
    fEvent->AddObject(spdv);
  }
  if (fDoPriv) {
    AliESDVertex *priv = new AliESDVertex();
    priv->SetName("PrimaryVertex");
    fEvent->AddObject(priv);
  }
  if (fDoEmCs) {
    fEvent->AddObject(new AliESDCaloCells("EMCALCells","EMCALCells from MyEsdTrimTask"));
  }
  if (fDoPCs) {
    fEvent->AddObject(new AliESDCaloCells("PHOSCells","PHOSCells from MyEsdTrimTask"));
  }
  if (fDoEmT) {
    AliESDCaloTrigger *etrig = new AliESDCaloTrigger;
    etrig->SetName("EMCALTrigger");
    fEvent->AddObject(etrig);
  }
  if (fDoPT) {
    AliESDCaloTrigger *ptrig = new AliESDCaloTrigger;
    ptrig->SetName("PHOSTrigger");
    fEvent->AddObject(ptrig);
  }
  if (fDoTracks) {
    TClonesArray *arr = new TClonesArray("AliESDtrack",0);
    arr->SetName("Tracks");
    fEvent->AddObject(arr);
  }
  //fEvent->AddObject(new AliTOFHeader());
  fEvent->GetStdContent();
  fEvent->WriteToTree(fTree);
  fTree->GetUserInfo()->Add(fEvent);
  TFile *file = OpenFile(1);
  fTree->SetDirectory(file);
  fTree->SetAutoFlush(-1024*1024*1024);
  fTree->SetAutoSave(-1024*1024*1024);
#endif
  PostData(1,fTree);
}
#endif

#if 0
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer();
    
    AliAnalysisTaskCaloFilter * filter = new AliAnalysisTaskCaloFilter("EMCALFilter");
    if(kUsePhysSel)filter->SelectCollisionCandidates(); 
    filter->SetCaloFilter(AliAnalysisTaskCaloFilter::kBoth); //kPHOS or kBoth
    filter->SwitchOnClusterCorrection();
    //filter->SetDebugLevel(10);
    AliEMCALRecoUtils * reco = filter->GetEMCALRecoUtils();
    reco->SetParticleType(AliEMCALRecoUtils::kPhoton);
    reco->SetW0(4.5);
        
    reco->SetPositionAlgorithm(AliEMCALRecoUtils::kPosTowerGlobal);

    TGeoHMatrix *matrix[4];
    
    double rotationMatrix[4][9] = {-0.014587, -0.999892, -0.002031, 0.999892, -0.014591,  0.001979, -0.002009, -0.002002,  0.999996,
				 -0.014587,  0.999892,  0.002031, 0.999892,  0.014591, -0.001979, -0.002009,  0.002002, -0.999996,
				 -0.345864, -0.938278, -0.003412, 0.938276, -0.345874,  0.003010, -0.004004, -0.002161,  0.999990,
				 -0.345861,  0.938280,  0.003412, 0.938276,  0.345874, -0.003010, -0.004004,  0.002161, -0.999990};
    
    double translationMatrix[4][3] = {0.351659,    447.576446,  176.269742,
				      1.062577,    446.893974, -173.728870,
				      -154.213287, 419.306156,  176.753692,
				      -153.018950, 418.623681, -173.243605};
    for(int j=0; j<4; j++)
      {
	matrix[j] = new TGeoHMatrix();
	matrix[j]->SetRotation(rotationMatrix[j]);
	matrix[j]->SetTranslation(translationMatrix[j]);
	matrix[j]->Print();
	filter->SetEMCALGeometryMatrixInSM(matrix[j],j);
      }
    
    
    filter->SwitchOnLoadOwnEMCALGeometryMatrices();
    
    reco->SetNonLinearityFunction(AliEMCALRecoUtils::kNoCorrection);

    //Time dependent corrections    
    //Recover file from alien  /alice/cern.ch/user/g/gconesab/TimeDepCorrectionDB
    reco->SwitchOnTimeDepCorrection();
    char cmd[200] ;
    sprintf(cmd, ".!tar xvfz CorrectionFiles.tgz") ;
    gROOT->ProcessLine(cmd) ;

    //Recalibration factors
    //Recover the file from alien  /alice/cern.ch/user/g/gconesab/RecalDB
    reco->SwitchOnRecalibration();
    TFile * f = new TFile("RecalibrationFactors.root","read");
    TH2F * h0 = (TH2F*)f->Get("EMCALRecalFactors_SM0");
    TH2F * h1 = (TH2F*)f->Get("EMCALRecalFactors_SM1");
    TH2F * h2 = (TH2F*)f->Get("EMCALRecalFactors_SM2");
    TH2F * h3 = (TH2F*)f->Get("EMCALRecalFactors_SM3");
    
    reco->SetEMCALChannelRecalibrationFactors(0,h0);
    reco->SetEMCALChannelRecalibrationFactors(1,h1);
    reco->SetEMCALChannelRecalibrationFactors(2,h2);
    reco->SetEMCALChannelRecalibrationFactors(3,h3);
    
    //Bad channels
    //Recover the file from alien  /alice/cern.ch/user/g/gconesab/BadChannelsDB
    reco->SwitchOnBadChannelsRemoval();
    reco->SwitchOnDistToBadChannelRecalculation();
    TFile * fbad = new TFile("BadChannels.root","read");
    TH2I * hbad0 = (TH2I*)fbad->Get("EMCALBadChannelMap_Mod0");
    TH2I * hbad1 = (TH2I*)fbad->Get("EMCALBadChannelMap_Mod1");
    TH2I * hbad2 = (TH2I*)fbad->Get("EMCALBadChannelMap_Mod2");
    TH2I * hbad3 = (TH2I*)fbad->Get("EMCALBadChannelMap_Mod3");
    reco->SetEMCALChannelStatusMap(0,hbad0);
    reco->SetEMCALChannelStatusMap(1,hbad1);
    reco->SetEMCALChannelStatusMap(2,hbad2);
    reco->SetEMCALChannelStatusMap(3,hbad3);

    //reco->Print("");
    filter->PrintInfo(); 
    mgr->AddTask(filter);
        
    //AliAnalysisDataContainer *cout_cuts2 = mgr->CreateContainer("Cuts", TList::Class(), 
    //					       AliAnalysisManager::kOutputContainer, "pi0calib.root");
    
    mgr->ConnectInput  (filter,  0, cinput1);
    mgr->ConnectOutput (filter, 0, coutput1 );
    //mgr->ConnectOutput (filter, 2, cout_cuts2);
    TString outputFile = AliAnalysisManager::GetCommonFileName(); 

    if(kDoConversionAnalysis && kUsePhysSel){
      TString arguments = "-run-on-train -use-own-xyz  -force-aod -mc-off ";
      gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/ConfigGammaConversion.C");
      AliAnalysisTaskGammaConversion * taskGammaConversion = 
	ConfigGammaConversion(arguments,mgr->GetCommonInputContainer());
      taskGammaConversion->SelectCollisionCandidates();
    }

    //Collision task
    AliCollisionNormalizationTask * taskNorm = new AliCollisionNormalizationTask("TaskNormalization");
    taskNorm->SetMC(kFALSE);
    taskNorm->SelectCollisionCandidates();
    mgr->AddTask(taskNorm);
    
    mgr->ConnectInput(taskNorm,0,cinput1);
    AliAnalysisDataContainer *    cOutputNorm = mgr->CreateContainer("Normalization", TList::Class(), AliAnalysisManager::kOutputContainer,outputFile.Data());
    mgr->ConnectOutput(taskNorm, 1, cOutputNorm);
    

    //Counter task
    AliAnalysisTaskCounter * counter = new AliAnalysisTaskCounter("Counter");
    counter->SwitchOnCaloFilterPatch();
    counter->SelectCollisionCandidates();  
    AliAnalysisDataContainer *coutputCount = 
      mgr->CreateContainer("counter", TList::Class(), AliAnalysisManager::kOutputContainer, outputFile.Data());
    //AliAnalysisDataContainer *coutput2  = mgr->CreateContainer("Count", TList::Class(), 
    //                                                 AliAnalysisManager::kOutputContainer, 
    //                                                 Form("%s:Count",outputFile.Data()));
    mgr->AddTask(counter);
    mgr->ConnectInput  (counter,  0, cinput1);
    mgr->ConnectOutput (counter, 1, coutputCount);

    //gROOT->LoadMacro("AddTaskCalorimeterQA.C");
    //AliAnalysisTaskParticleCorrelation * qa = AddTaskCalorimeterQA(kInputData,kFALSE,kFALSE);
    
    //-----------------------
    // Run the analysis
    //-----------------------    
    TString smode = "";
    if (mode==mLocal || mode == mLocalCAF) 
      smode = "local";
    else if (mode==mPROOF) 
      smode = "proof";
    else if (mode==mGRID) 
      smode = "local";
    
    mgr->InitAnalysis();
    mgr->PrintStatus();
    mgr->StartAnalysis(smode.Data(),chain);

    cout <<" Analysis ended sucessfully "<< endl ;
    
  }
  else cout << "Chain was not produced ! "<<endl;
  
  //sprintf(cmd, ".! rm -rf CorrectionFiles") ;
  
}
#endif


