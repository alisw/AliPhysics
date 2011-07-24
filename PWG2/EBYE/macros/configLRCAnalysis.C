AliESDtrackCuts* createAliLRCcuts(char* mode)
{
if(mode=="Global2")
{
  AliESDtrackCuts *cuts = new AliESDtrackCuts("lrcCuts","lrcCuts");
  cuts->SetMinNClustersTPC(70);
  cuts->SetMinNClustersITS(2);
  cuts->SetMaxChi2PerClusterTPC(4.0);
  cuts->SetRequireTPCRefit();
  cuts->SetRequireITSRefit();
  cuts->SetAcceptKinkDaughters(kFALSE);
  cuts->SetMaxDCAToVertexXY(0.5);
  cuts->SetMaxDCAToVertexZ(0.5);
  cuts->SetPtRange(0.15,100);
  cuts->SetEtaRange(-1.8,1.8);
  cuts->SaveHistograms("trackCuts");

  return cuts;
}
if(mode=="Global1")
{
  AliESDtrackCuts *cuts = new AliESDtrackCuts("lrcCuts","lrcCuts");
  cuts->SetMinNClustersTPC(80);
  cuts->SetMinNClustersITS(2);
  cuts->SetMaxChi2PerClusterTPC(4.0);
  cuts->SetRequireTPCRefit();
  cuts->SetRequireITSRefit();
  cuts->SetAcceptKinkDaughters(kFALSE);
  cuts->SetMaxDCAToVertexXY(0.3);
  cuts->SetMaxDCAToVertexZ(0.3);
  cuts->SetPtRange(0.15,100);
  cuts->SetEtaRange(-1.8,1.8);
  cuts->SaveHistograms("trackCuts");

  return cuts;
}
if(mode=="GlobalNoITS1")
{
  AliESDtrackCuts *cuts = new AliESDtrackCuts("lrcCuts","lrcCuts");
  cuts->SetMinNClustersTPC(80);
  cuts->SetMaxChi2PerClusterTPC(4.0);
  cuts->SetRequireTPCRefit();
  cuts->SetAcceptKinkDaughters(kFALSE);
  cuts->SetMaxDCAToVertexXY(0.3);
  cuts->SetMaxDCAToVertexZ(0.3);
  cuts->SetPtRange(0.15,100);
  cuts->SetEtaRange(-1.8,1.8);
  cuts->SaveHistograms("trackCuts");
  return cuts;
}
}
//===========================================================================
AliESDtrackCuts* createAliLRCcutsOld()
{
    AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("AliESDtrackCutsLRC");
    esdTrackCuts->DefineHistograms(1);

    Double_t cov1, cov2, cov3, cov4, cov5;
    Double_t nSigma;
    Double_t maxDCAtoVertex, maxDCAtoVertexXY, maxDCAtoVertexZ;
    Double_t minNClustersTPC;
    Double_t maxChi2PerClusterTPC;
    Double_t minPt, maxPt;
    cov1 = 2;
    cov2 = 2;
    cov3 = 0.5;
    cov4 = 0.5;
    cov5 = 2;
    nSigma = 3;
    minNClustersTPC = 50;
    maxChi2PerClusterTPC = 3.5;

    esdTrackCuts->SetMaxCovDiagonalElements(cov1, cov2, cov3, cov4, cov5);
    esdTrackCuts->SetMaxNsigmaToVertex(nSigma);
    esdTrackCuts->SetRequireSigmaToVertex(kTRUE);
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(minNClustersTPC);
    esdTrackCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);

    return esdTrackCuts;
}

//===========================================================================
AliAnalysisTaskLRC* createLRCtaskSkeleton(char* name="Task_LRC", Bool_t RunKine=kFALSE)
{
	AliAnalysisTaskLRC *taskLRC = new AliAnalysisTaskLRC(name,RunKine);
	taskLRC->SetMinPtLimit(0.2);
        taskLRC->SetMaxPtLimit(3.5);
	taskLRC->SetCheckForVtxPosition(kTRUE);
	taskLRC->SetVtxDiamond(0.4,0.4,5.0);
	return taskLRC;
}
//===========================================================================
void addAliLRCProcessors(AliAnalysisTaskLRC* taskLRC,Bool_t AddPhiWindows=kFALSE)
{
Int_t NwindowSets=1;
if(AddPhiWindows)Int_t NwindowSets=4;

for(int i=0;i<NwindowSets;i++)
{
	//FB
	taskLRC->AddLRCProcess(new AliLRCProcess(-0.2,-0.0,0.0,0.2));
	taskLRC->AddLRCProcess(new AliLRCProcess(-0.4,-0.0,0.0,0.4));
	taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.0,0.0,0.6));
	taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.0,0.0,0.8));
	//0.2 gap
	taskLRC->AddLRCProcess(new AliLRCProcess(-0.4,-0.2,0.2,0.4));
	taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.4,0.4,0.6));
	taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.6,0.6,0.8));
	//0.4 gap
	taskLRC->AddLRCProcess(new AliLRCProcess(-0.6,-0.2,0.2,0.6));
	taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.4,0.4,0.8));
	//0.6 gap
	taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,-0.2,0.2,0.8));
	
	//FULL
	taskLRC->AddLRCProcess(new AliLRCProcess(-0.8,0.8,-0.8,0.8));
}	

if(AddPhiWindows)
{
	for(Int_t i=11; i < 22; i++)
	{
		(dynamic_cast<AliLRCProcess*> (taskLRC->Proc(i)))->SetForwardWindowPhi(0,TMath::Pi()/2.0);
		(dynamic_cast<AliLRCProcess*> (taskLRC->Proc(i)))->SetBackwardWindowPhi(TMath::Pi()/2.0,TMath::Pi());
	}
	
	for(Int_t i=22; i < 33; i++)
	{
		(dynamic_cast<AliLRCProcess*> (taskLRC->Proc(i)))->SetForwardWindowPhi(0,TMath::Pi()/2.0);
		(dynamic_cast<AliLRCProcess*> (taskLRC->Proc(i)))->SetBackwardWindowPhi(0,TMath::Pi()/2.0);
	}
	
	for(Int_t i=33; i < 44; i++)
	{
		(dynamic_cast<AliLRCProcess*> (taskLRC->Proc(i)))->SetForwardWindowPhi(0,TMath::Pi()/2.0);
		(dynamic_cast<AliLRCProcess*> (taskLRC->Proc(i)))->SetBackwardWindowPhi(TMath::Pi(),TMath::Pi()*3.0/2.0);
	}
} // endif AddPhiWindows

}

//===========================================================================
void configureLRCtaskOutput(AliAnalysisTaskLRC* taskLRC,TString OutputRootFolder=":PWG2LRC")
{
if(!taskLRC){return ;}
AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
if (!mgr) 
{
	Error("configureLRCtaskOutput", "No analysis manager to connect to.");
	return ;
}

if (!mgr->GetInputEventHandler()) 
{
      Error("AddTaskLRC", "This task requires an input event handler");
      return ;
}

TString type = mgr->GetInputEventHandler()->GetDataType(); 

TString outputFileName= mgr->GetCommonFileName();
if(outputFileName=="")outputFileName="LRC.NEWIO.root";
outputFileName += OutputRootFolder ;

TString OutName;
OutName=taskLRC->GetName();

AliAnalysisDataContainer *cout_LRC = mgr->CreateContainer(OutName, TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
mgr->ConnectInput(taskLRC, 0, mgr->GetCommonInputContainer());
mgr->ConnectOutput(taskLRC, 1, cout_LRC);


cout << "N of LRC Processors ="<< taskLRC->GetListOfProcessors()->GetEntries() <<"\n";
for(Int_t i=0; i < taskLRC->GetListOfProcessors()->GetEntries(); i++)
{
	mgr->ConnectOutput(taskLRC,taskLRC->Proc(i)->GetOutputSlotNumber(),mgr->CreateContainer(((taskLRC->Proc(i)->GetShortDef()+taskLRC->GetName())+=i),TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName));
}
}
