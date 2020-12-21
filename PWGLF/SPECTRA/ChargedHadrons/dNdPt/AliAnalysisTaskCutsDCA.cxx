#include "THnSparse.h"
#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliHeader.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliMultSelection.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisUtils.h"
#include "AliESDUtils.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliAnalysisTaskCutsDCA.h"

class AliAnalysisTaskCutsDCA;    

using namespace std;

ClassImp(AliAnalysisTaskCutsDCA) 

THnSparseD* AliAnalysisTaskCutsDCA::fSparseTmp = 0;

//_____________________________________________________________________________

AliAnalysisTaskCutsDCA::AliAnalysisTaskCutsDCA() 
    : AliAnalysisTaskSE()
    , fEvent(0)
    , fESD(0)
    , fAOD(0)
    , fMCEvent(0)
    , fMCStack(0)
    , fMCHeader(0)
    , fMCGenHeader(0)
    , fEventMCzv(0)
    , fEventCent(0)
    , fEventNtracks(0)
    , fEventMultMB(0)
    , fEventMultV0M(0)
    , fEventMCb(0)
    , fEventMCnPrim(0)
    , fEventMCnPrimV0M(0)
    , fEventIsTrigger(0)
    , fEventHasVertex(0)
    , fESDtrackCuts0(0)
    , fCentEstimator("V0M")
    , fOutputList(0)
    , fLogHist(0)
    , fLogEvent(0)
    , fTrigHist(0)
    , fTrigHistU(0)
    , fRunHist(0)
    , fHistTrack1(0)
    , fHistEvent1(0)
    , fHistVertex(0)
    , fHistVZERO(0)    
    , fHistMCTrack1(0)
    , fHistMCEvent1(0)    
    , fHistMCPrim1(0)
    , fHistMultCent(0)
    , fHistMCMultCent(0)
{
    //constructor
}

//_____________________________________________________________________________

AliAnalysisTaskCutsDCA::AliAnalysisTaskCutsDCA(const char* name)
    : AliAnalysisTaskSE(name)
    , fEvent(0)
    , fESD(0)
    , fAOD(0)
    , fMCEvent(0)
    , fMCStack(0)
    , fMCHeader(0)
    , fMCGenHeader(0)
    , fEventMCzv(0)
    , fEventCent(0)
    , fEventNtracks(0)
    , fEventMultMB(0)
    , fEventMultV0M(0)
    , fEventMCb(0)
    , fEventMCnPrim(0)
    , fEventMCnPrimV0M(0)
    , fEventIsTrigger(0)
    , fEventHasVertex(0)
    , fESDtrackCuts0(0)  
    , fCentEstimator("V0M")
    , fOutputList(0)
    , fLogHist(0)
    , fRunHist(0)
    , fLogEvent(0)
    , fTrigHist(0)
    , fTrigHistU(0)
    , fHistTrack1(0)
    , fHistEvent1(0)
    , fHistVertex(0)
    , fHistVZERO(0)
    , fHistMCTrack1(0)
    , fHistMCEvent1(0)    
    , fHistMCPrim1(0)
    , fHistMultCent(0)
    , fHistMCMultCent(0)
{
    // constructor
     DefineInput(0, TChain::Class());    // define input
     DefineOutput(1, TList::Class());    // define ouptut 
}

//_____________________________________________________________________________

AliAnalysisTaskCutsDCA::~AliAnalysisTaskCutsDCA()
{
    // destructor
    if(fOutputList) { delete fOutputList; }
    if(fSparseTmp)  { delete fSparseTmp; }
}

//_____________________________________________________________________________

void AliAnalysisTaskCutsDCA::UserCreateOutputObjects()
{
    // create output list
    fOutputList = new TList();         
    fOutputList->SetOwner(kTRUE);

    // create histograms    
    fLogHist = CreateLogHist("fLogHist");
    fOutputList->Add(fLogHist);
    
    fLogEvent = CreateLogHist("fLogEvent");
    fOutputList->Add(fLogEvent);    
    
    fTrigHist = CreateLogHist("fTrigHist");
    fOutputList->Add(fTrigHist);
    
    fTrigHistU = CreateLogHist("fTrigHistU");
    fOutputList->Add(fTrigHistU);
    
    fRunHist = CreateLogHist("fRunHist");
    fOutputList->Add(fRunHist);
    
    AddAxis("pt");
    AddAxis("eta","#eta",2,-0.8,+0.8);
    AddAxis("zv",2,-10.0,10.0);
//     AddAxis("centV0Mold",22,-5,105);
    AddAxis("centV0Mnew",22,-5,105);
    AddAxis("nTracks",122,-1.5,122.0-1.5);
    AddAxis("multMB",202,-1.5,202.0-1.5);
    AddAxis("multV0",302,-2.0,302.0-2.0);
//     AddAxis("multV0A","multV0A","varsig35");
//     AddAxis("multV0C","multV0C","varsig35");
//     AddAxis("maxpt","p_{t,max} (GeV/c)","ptfew");
    AddAxis("Q",2,-2.0,+2.0);
    fHistTrack1 = CreateHist("fHistTrack1");
    fOutputList->Add(fHistTrack1);
        
//     AddAxis("EventSpecie",100,-0.5,100.0-0.5);
    AddAxis("zv",20,-20.0,20.0);
//     AddAxis("centV0Mold",22,-5,105);
    AddAxis("centV0Mnew",22,-5,105);    
    AddAxis("nTracks",2000,-0.5,2000.0-0.5);
    AddAxis("multMB",4000,-0.5,4000.0-0.5);
    AddAxis("multV0",60000,-2,60000.0-2.);
    AddAxis("multV0A",30000,-2,30000.0-2.);
    AddAxis("multV0C",30000,-2,30000.0-2.);
//     AddAxis("maxpt","p_{t,max} (GeV/c)","ptfew");
    fHistEvent1 = CreateHist("fHistEvent1");
    fOutputList->Add(fHistEvent1);
    
    AddAxis("zvTRK",60,-30.0,30.0);
    AddAxis("statusTRK",2,-0.5,1.5);
    AddAxis("nContribTRK",250,-1.5,250-1.5);
    AddAxis("zvSPD",60,-30.0,30.0);
    AddAxis("statusSPD",2,-0.5,1.5);
    AddAxis("nContribSPD",250,-1.5,250-1.5);
    AddAxis("zvTPC",60,-30.0,30.0);
    AddAxis("statusTPC",2,-0.5,1.5);
    AddAxis("nContribTPC",250,-1.5,250-1.5);    
    fHistVertex = CreateHist("fHistVertex");
    fOutputList->Add(fHistVertex);
    
    // create histograms    
    AddAxis("pt","pt");
    fHistMCTrack1 = CreateHist("fHistMCTrack1");    
    fOutputList->Add(fHistMCTrack1);
    
    AddAxis("MCzv",20,-20.0,20.0);
    AddAxis("centV0Mnew",22,-5,105);
    AddAxis("nTracks",2000,-0.5,2000.0-0.5);
    AddAxis("multMB",4000,-0.5,4000.0-0.5);
    AddAxis("multV0",60000,-2,60000.0-2.);
    AddAxis("MCb",200,0.0,20.0);
    AddAxis("MCnPrim",4000,-0.5,4000.0-0.5);
    AddAxis("MCnPrimV0M",8000,-0.5,8000.0-0.5);    
    AddAxis("isTriggered",2,-0.5,1.5);
    AddAxis("isRecVertex",2,-0.5,1.5);
    fHistMCEvent1 = CreateHist("fHistMCEvent1");
    fOutputList->Add(fHistMCEvent1);
    
    AddAxis("nPrimMC",20000,-0.5,20000.0-0.5);
    AddAxis("nPrimMCcount",20000,-0.5,20000.0-0.5);
    AddAxis("nPhysPrimMC",20000,-0.5,20000.0-0.5);
    fHistMCPrim1 = CreateHist("fHistMCPrim1");
    fOutputList->Add(fHistMCPrim1);          

    
    PostData(1, fOutputList);           // postdata 
   
}

//_____________________________________________________________________________
void AliAnalysisTaskCutsDCA::Initialize() 
{
    fEventMCzv = 0;
    fEventCent = -1;
    fEventNtracks = 0;
    fEventMultMB = 0;
    fEventMultV0M = 0;
    fEventMCb = -1;
    fEventMCnPrim = 0;
    fEventMCnPrimV0M = 0;
    fEventIsTrigger = kFALSE;
    fEventHasVertex = kFALSE;
}

//_____________________________________________________________________________
void AliAnalysisTaskCutsDCA::UserExec(Option_t *)
{           
    Initialize();
    
    fEvent = InputEvent();
    if (!fEvent) { Log("noEvent"); return; }
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());   
       
    if (fAOD) {
        TString firedTriggerClasses=((AliAODEvent*)fAOD)->GetFiredTriggerClasses();
        if (firedTriggerClasses.Length() > 0) {
            Log(fTrigHist,firedTriggerClasses.Data());
        } else {
            Log(fTrigHist,"noTriggerString");
        }
    }
    if (fESD) {
        TString firedTriggerClasses=((AliESDEvent*)fESD)->GetFiredTriggerClasses();
        if (firedTriggerClasses.Length() > 0)  {
            Log(fTrigHist,firedTriggerClasses.Data());
        } else {
            Log(fTrigHist,"noTriggerString");
        }
    }
    
    if (fESD) {    
        Log("ESD"); 
        if (!AnalyzeESD()) { Log("noAnalyzeESD"); }
    } else if (fAOD) { 
        if (!AnalyzeAOD()) { Log("noAnalyzeAOD"); }
        Log("AOD");
    } else {
        Log("noESDorAOD");
    }

    fMCEvent = MCEvent();
    if (fMCEvent) {
      Log("MC");
      fMCHeader = fMCEvent->Header();
      if (!fMCHeader) { Log("noMCHeader"); }
      fMCGenHeader = fMCHeader->GenEventHeader();
      if (!fMCGenHeader) { 
          Log("noMCGenHeader");           
      } else { 
          TString s = "mcHeader=";
          s += fMCGenHeader->GetName();
          Log(s.Data());              
      }
      fMCStack = fMCEvent->Stack();    
      if (!fMCStack) { Log("noMCStack"); }
      if (!AnalyzeMC()) { Log("noAnalyzeMC"); }
    } 
    
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kAny) { Log(fLogEvent,"kAny"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kAnyINT) { Log(fLogEvent,"kAnyINT"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7) { Log(fLogEvent,"kINT7"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMUON) { Log(fLogEvent,"kMUON"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kHighMult) { Log(fLogEvent,"kHighMult"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kEMC1) { Log(fLogEvent,"kEMC1"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kCINT5) { Log(fLogEvent,"kCINT5"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kCMUS5) { Log(fLogEvent,"kCMUS5"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMUSH7) { Log(fLogEvent,"kMUSH7"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMUL7) { Log(fLogEvent,"kMUL7"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMUU7) { Log(fLogEvent,"kMUU7"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kEMC7) { Log(fLogEvent,"kEMC7"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kEMC8) { Log(fLogEvent,"kEMC8"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMUS7) { Log(fLogEvent,"kMUS7"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kPHI1) { Log(fLogEvent,"kPHI1"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kPHI7) { Log(fLogEvent,"kPHI7"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kEMCEJE) { Log(fLogEvent,"kEMCEJE"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kEMCEGA) { Log(fLogEvent,"kEMCEGA"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kCentral) { Log(fLogEvent,"kCentral"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kSemiCentral) { Log(fLogEvent,"kSemiCentral"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kDG5) { Log(fLogEvent,"kDG5"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kZED) { Log(fLogEvent,"kZED"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kSPI7) { Log(fLogEvent,"kSPI7"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kSPI) { Log(fLogEvent,"kSPI"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT8) { Log(fLogEvent,"kINT8"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMuonSingleLowPt8) { Log(fLogEvent,"kMuonSingleLowPt8"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMuonSingleHighPt8) { Log(fLogEvent,"kMuonSingleHighPt8"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMuonLikeLowPt8) { Log(fLogEvent,"kMuonLikeLowPt8"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMuonUnlikeLowPt8) { Log(fLogEvent,"kMuonUnlikeLowPt8"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMuonUnlikeLowPt0) { Log(fLogEvent,"kMuonUnlikeLowPt0"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kTRD) { Log(fLogEvent,"kTRD"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kFastOnly) { Log(fLogEvent,"kFastOnly"); }
    if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kUserDefined) { Log(fLogEvent,"kUserDefined"); }   
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskCutsDCA::AnalyzeAOD()
{
    // put AOD analysis here
    return kFALSE;  
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskCutsDCA::AnalyzeESD()
{    
  
    UInt_t eventSpecie = fESD->GetEventSpecie();
    Double_t centralityV0M = -1;    
    
    //Physics Selection with proper trigger
    Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7); 
    if (!isSelected) {
        Log("eventNotSelected");
        return kFALSE;
    }
    // skipincomplete daq events
    if (fESD->IsIncompleteDAQ()) {
        Log("event.IsIncompleteDAQ");
        return kFALSE;
    }      
    
    TString firedTriggerClasses = fESD->GetFiredTriggerClasses();
    if (firedTriggerClasses.Length() > 0) {
        Log(fTrigHistU,firedTriggerClasses.Data());
    } else {
        Log(fTrigHistU,"noTriggerString");
    }
    Log(fCentEstimator.Data());
    
    AliCentrality *centrality = fESD->GetCentrality();
    if (!centrality) { 
        Log("noOldCentrality");         
    } else {
        centralityV0M = centrality->GetCentralityPercentile(fCentEstimator);        
        if (centralityV0M < 0.) { Log("CentOldLess0"); }
        if (centralityV0M > 100.) { Log("CentOldGreater100"); }
        if (centralityV0M == 0.) { Log("CentOldEqual0"); }
        if (centralityV0M == 100.) { Log("CentOldEqual100"); }
    }
    
    Double_t centralityNewV0M = -1;
    AliMultSelection *MultSelection = 0; 
    MultSelection = (AliMultSelection*) fEvent->FindListObject("MultSelection");
    if(!MultSelection) {
        Log("noNewCentrality");
    } else {
        centralityNewV0M = MultSelection->GetMultiplicityPercentile(fCentEstimator);
        if (centralityNewV0M < 0.) { Log("CentNewLess0"); }
        if (centralityNewV0M > 100.) { Log("CentNewGreater100"); }
        if (centralityNewV0M == 0.) { Log("CentNewEqual0"); }
        if (centralityNewV0M == 100.) { Log("CentNewEqual100"); }        
    }
   
    //vzero
    AliESDVZERO* vzero = fESD->GetVZEROData();
    if (!vzero) { 
        Log("noVZEROData");
        return kFALSE;
    }
    
    // spd cluster vs. tracklet rejectkion (OFF)        
    AliAnalysisUtils utils;
    if (utils.IsSPDClusterVsTrackletBG(fESD)){
      Log("utils.IsSPDClusterVsTrackletBG");
    }   
    // first event in chunk 
    if (utils.IsFirstEventInChunk(fESD)){
      Log("utils.IsFirstEventInChunk");
    }
    if (utils.IsPileUpMV(fESD)){
      Log("utils.IsPileUpMV");
    }
    if (utils.IsOutOfBunchPileUp(fESD)){
      Log("utils.IsOutOfBunchPileUp");
    }
    if (utils.IsPileUpEvent(fESD)){
      Log("utils.IsPileUpEvent");
    }
    if (utils.IsPileUpSPD(fESD)){
      Log("utils.IsPileUpSPD");
    }
    if (!utils.IsVertexSelected2013pA(fESD)){
      Log("utils.NOT.IsVertexSelected2013pA");
    }
    if (fESD->IsPileupFromSPD(5,0.8)){
      Log("event.IsPileupFromSPD(5,0.8)");
    }
       
    if (fESD->IsPileupFromSPD(5,0.8)){
      Log("event.IsPileupFromSPD(5,0.8)");
    }

    

    Double_t rawv0a = vzero->GetMTotV0A();    
    Double_t rawv0c = vzero->GetMTotV0C();
    
    Int_t multMB = -1;
    Double_t zv = 0;
    Double_t zvRes = 1e10;
    Double_t zvResTRK = 1e10;
    Double_t zvResTPC = 1e10;
    Double_t zvResSPD = 1e10;
    Int_t nContribTPC = -1;
    Int_t nContribTRK = -1;
    Int_t nContribSPD = -1;
    Int_t runnumber = fESD->GetRunNumber();
    Double_t zvTPC = 0;
    Double_t zvTRK = 0;
    Double_t zvSPD = 0;
    Bool_t status = kTRUE;
    Bool_t statusTPC = kFALSE;
    Bool_t statusTRK = kFALSE;
    Bool_t statusSPD = kFALSE;    
    const AliESDVertex* vtxTRK = fESD->GetPrimaryVertexTracks();
    if (vtxTRK) { 
        nContribTRK = vtxTRK->GetNContributors();
        zvTRK = vtxTRK->GetZ();
    zvResTRK = vtxTRK->GetZRes();
    statusTRK = vtxTRK->GetStatus();
    if (!statusTRK) { Log("noVertexTRK"); }
    } else {
        Log("PrimVtxTRK==0");
    }
    const AliESDVertex* vtxSPD = fESD->GetPrimaryVertexSPD();
    if (vtxSPD) { 
        nContribSPD = vtxSPD->GetNContributors();
        zvSPD = vtxSPD->GetZ();
    zvResSPD = vtxSPD->GetZRes();
    statusSPD = vtxSPD->GetStatus();
    if (!statusSPD) { Log("noVertexSPD"); }
    } else {
        Log("PrimVtxSPD==0");
    }    
    const AliESDVertex* vtxTPC = fESD->GetPrimaryVertexTPC();
    if (vtxTPC) { 
        nContribTPC = vtxTPC->GetNContributors();
        zvTPC = vtxTPC->GetZ();
    zvResTPC = vtxTPC->GetZRes();
    statusTPC = vtxTPC->GetStatus();
    if (!statusTPC) { Log("noVertexTPC"); }
    } else {
        Log("PrimVtxTPC==0");
    }       
    if (statusTRK) {
        multMB = nContribTRK;
        zv = zvTRK;
    zvRes = zvResTRK;   
    Log("usedVertexTracks");
    } else if (statusSPD) {
        multMB = nContribSPD;
        zv = zvSPD;     
    zvRes = zvResSPD;
    Log("usedVertexSPD");
    } else if (statusTPC) {
        multMB = nContribTPC;
        zv = zvTPC; 
    zvRes = zvResTPC;
    Log("usedVertexTPC");
    } else {
        status = kFALSE;
        Log("noPrimVtx");
    }
    // fill vertex cuts only for evets with vertex
    if (status) {
      if (zvRes > 0.25) {
    Log("PrimVtx.ResolutionGreater0.25");
      }    
    }
    if (statusSPD) {
      if (vtxSPD->IsFromVertexerZ() && vtxSPD->GetDispersion() > 0.04) {
    Log("PrimVtxSPD.fromVertexerZ&DispersionGreater0.04");
      }
      if (vtxSPD->GetDispersion() > 0.04) {
    Log("PrimVtxSPD.DispersionGreater0.04");
      }
      if (vtxSPD->GetDispersion() > 0.02) {
    Log("PrimVtxSPD.DispersionGreater0.02");
      }
    }
    if (statusTRK && statusSPD) {
        if (TMath::Abs(zvTRK-zvSPD) > 0.5) {
      Log("abs(zvTRK-zvSPD)Greater0.5");
    }
    }
    
    Double_t corrv0a = AliESDUtils::GetCorrV0A(rawv0a,zv);
    Double_t corrv0c = AliESDUtils::GetCorrV0C(rawv0c,zv);
    Double_t rawv0 = rawv0a+rawv0c;
    Double_t corrv0 = corrv0a+corrv0c;
    Double_t v0a = rawv0a; //rawv0a/(1.-0.00353*zv);
    Double_t v0c = rawv0c; //rawv0c/(1.+0.00471*zv);
    Double_t v0m = ((v0a)+(v0c))/(((zv)<-14.5)*((-1.5616500000)+(-0.3058090000)*((zv)-(1.6873000000))+(-0.0094082000)*TMath::Power((zv)-(1.6873000000),2)) + (TMath::Abs((zv))<14.5)*((0.9996770000)+(0.0014770600)*((zv)-(0.0003504880))+(-0.0000044498)*TMath::Power((zv)-(0.0003504880),2)+(0.0000058547)*TMath::Power((zv)-(0.0003504880),3)+(-0.0000006434)*TMath::Power((zv)-(0.0003504880),4))+((zv)>14.5)*((0.4462690000)+(0.1055920000)*((zv)-(4.0898700000))+(-0.0049480200)*TMath::Power((zv)-(4.0898700000),2)));    

    // first loop for nacc and maxpt
    Int_t nTracks = fESD->GetNumberOfTracks();
    AliESDtrack* track = 0;    
    Double_t maxpt = 0; 
    Int_t nacc = 0;

    for(Int_t iTrack = 0; iTrack < nTracks; iTrack++) {     // track loop
        track = static_cast<AliESDtrack*>(fESD->GetTrack(iTrack));
        if (!track) {
            Log("noESDtrackLoop0");
            continue;
        }    
        if (! fESDtrackCuts0 ) { continue; } 
        if (! fESDtrackCuts0->AcceptTrack(track) ) { continue; } 
        nacc++;
        if (track->Pt() > maxpt) { maxpt = track->Pt(); }
    }

    // second loop for track hist filling
    for(Int_t iTrack = 0; iTrack < nTracks; iTrack++) {     // track loop
        track = static_cast<AliESDtrack*>(fESD->GetTrack(iTrack));
        if (!track) {
            Log("noESDtrack");
            continue;
        }
        Double_t pt = track->Pt();
        Double_t eta = track->Eta();
        Double_t q   = track->GetSign();

        if (! fESDtrackCuts0 ) { continue; } 
        if (! fESDtrackCuts0->AcceptTrack(track) ) { continue; } 
        FillHist(fHistTrack1, pt,eta,zv,centralityNewV0M,nacc,multMB,v0m,q);
        }
        
    //FillHist(fHistEvent1, eventSpecie,zv,centralityV0M,centralityNewV0M,nacc,multMB,v0m,v0a,v0c,maxpt);
    //FillHist(fHistEvent1 ,zv,centralityV0M,centralityNewV0M,nacc,multMB,v0m,v0a,v0c);
    FillHist(fHistEvent1 ,zv,centralityNewV0M,nacc,multMB,v0m,v0a,v0c);    
    FillHist(fHistVertex, zvTRK, statusTRK, nContribTRK,zvSPD,statusSPD,nContribSPD,zvTPC,statusTPC,nContribTPC);
    TString srun;
    srun += runnumber;
    Log(fRunHist,srun.Data());
    return kTRUE;
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskCutsDCA::AnalyzeMC()
{
    UInt_t eventSpecie = fESD->GetEventSpecie();
    Double_t centralityV0M = -1;    
    
    //Physics Selection with proper trigger
    fEventIsTrigger = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);
    
    AliMultSelection *MultSelection = 0; 
    MultSelection = (AliMultSelection*) fEvent->FindListObject("MultSelection");
    if(!MultSelection) {
        Log("noNewCentralityMC");
    } else {
        fEventCent = MultSelection->GetMultiplicityPercentile(fCentEstimator);
    }
   
    
    Int_t multMB = -1;
    Double_t zv = 0;
    Double_t zvRes = 1e10;
    Double_t zvResTRK = 1e10;
    Double_t zvResTPC = 1e10;
    Double_t zvResSPD = 1e10;
    Int_t nContribTPC = -1;
    Int_t nContribTRK = -1;
    Int_t nContribSPD = -1;    
    Double_t zvTPC = 0;
    Double_t zvTRK = 0;
    Double_t zvSPD = 0;
    Bool_t status = kTRUE;
    Bool_t statusTPC = kFALSE;
    Bool_t statusTRK = kFALSE;
    Bool_t statusSPD = kFALSE;    
    const AliESDVertex* vtxTRK = fESD->GetPrimaryVertexTracks();
    if (vtxTRK) { 
        nContribTRK = vtxTRK->GetNContributors();
        zvTRK = vtxTRK->GetZ();
    zvResTRK = vtxTRK->GetZRes();
    statusTRK = vtxTRK->GetStatus();
    }
    const AliESDVertex* vtxSPD = fESD->GetPrimaryVertexSPD();
    if (vtxSPD) { 
        nContribSPD = vtxSPD->GetNContributors();
        zvSPD = vtxSPD->GetZ();
    zvResSPD = vtxSPD->GetZRes();
    statusSPD = vtxSPD->GetStatus();
    }
    const AliESDVertex* vtxTPC = fESD->GetPrimaryVertexTPC();
    if (vtxTPC) { 
        nContribTPC = vtxTPC->GetNContributors();
        zvTPC = vtxTPC->GetZ();
    zvResTPC = vtxTPC->GetZRes();
    statusTPC = vtxTPC->GetStatus();
    }
    if (statusTRK) {
        multMB = nContribTRK;
        zv = zvTRK;
    zvRes = zvResTRK;       
    } else if (statusSPD) {
        multMB = nContribSPD;
        zv = zvSPD;     
    zvRes = zvResSPD;    ;
    } else if (statusTPC) {
        multMB = nContribTPC;
        zv = zvTPC; 
    zvRes = zvResTPC;
    } else {
        status = kFALSE;
    }
    
    fEventHasVertex = status;
    fEventMultMB = multMB;

    // first loop for nacc and maxpt
    Int_t nTracks = fESD->GetNumberOfTracks();
    AliESDtrack* track = 0;    
    Double_t maxpt = 0; 
    Int_t nacc = 0;

    for(Int_t iTrack = 0; iTrack < nTracks; iTrack++) {     // track loop
        track = static_cast<AliESDtrack*>(fESD->GetTrack(iTrack));
        if (!track) {
            Log("noESDtrackLoop0");
            continue;
        }    
        if (! fESDtrackCuts0 ) { continue; } 
        if (! fESDtrackCuts0->AcceptTrack(track) ) { continue; } 
        nacc++;
        if (track->Pt() > maxpt) { maxpt = track->Pt(); }
    }
    fEventNtracks = nacc;
 
        //vzero
    AliESDVZERO* vzero = fESD->GetVZEROData();
    if (vzero) { 
      Double_t rawv0a = vzero->GetMTotV0A();    
      Double_t rawv0c = vzero->GetMTotV0C();
    
    Double_t corrv0a = AliESDUtils::GetCorrV0A(rawv0a,zv);
    Double_t corrv0c = AliESDUtils::GetCorrV0C(rawv0c,zv);
    Double_t rawv0 = rawv0a+rawv0c;
    Double_t corrv0 = corrv0a+corrv0c;
    Double_t v0a = rawv0a; //rawv0a/(1.-0.00353*zv);
    Double_t v0c = rawv0c; //rawv0c/(1.+0.00471*zv);
    Double_t v0m = ((v0a)+(v0c))/(((zv)<-14.5)*((-1.5616500000)+(-0.3058090000)*((zv)-(1.6873000000))+(-0.0094082000)*TMath::Power((zv)-(1.6873000000),2)) + (TMath::Abs((zv))<14.5)*((0.9996770000)+(0.0014770600)*((zv)-(0.0003504880))+(-0.0000044498)*TMath::Power((zv)-(0.0003504880),2)+(0.0000058547)*TMath::Power((zv)-(0.0003504880),3)+(-0.0000006434)*TMath::Power((zv)-(0.0003504880),4))+((zv)>14.5)*((0.4462690000)+(0.1055920000)*((zv)-(4.0898700000))+(-0.0049480200)*TMath::Power((zv)-(4.0898700000),2)));    
        
    fEventMultV0M = v0m;
        
        
    }
    
    
    
    UInt_t eventType = fMCEvent->GetEventType();
  
    Int_t nTracksMC = fMCEvent->GetNumberOfTracks();
    Int_t nPrimariesMC = fMCEvent->GetNumberOfPrimaries();
    Int_t nPrimCount = 0;
    Int_t nPhysPrim = 0;
    
    
    /*
    TString genname;
    Bool_t hasheaderlist = fMCEvent->GetCocktailGenerator(imc,genname);
    if (hasheaderlist) { 
        genname = "mcCocktail=" + genname;
        Log(genname.Data()); 
    } else {
        Log("noCocktailHeader");
    }      
    */    
    TArrayF vtxMC(3);
    fMCGenHeader->PrimaryVertex(vtxMC);
    fEventMCzv = vtxMC[2];
            
    AliGenHijingEventHeader* hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(fMCGenHeader);    
    if(!hijingGenHeader){
      AliGenCocktailEventHeader* genCocktailHeader = dynamic_cast<AliGenCocktailEventHeader*>(fMCGenHeader);
      if(genCocktailHeader) {
          TList* headerList = genCocktailHeader->GetHeaders();
          for (Int_t iH=0; iH < headerList->GetEntries(); iH++) {
            hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(headerList->At(iH));
            if (hijingGenHeader) break;
          }
      }
    }

    if (!hijingGenHeader) { 
        Log("noHIJINGGenHeader"); 
    } else {
        fEventMCb = hijingGenHeader->ImpactParameter();
    }
    
    AliMCParticle* mcpart = 0;
    Double_t pt = 0;    
    Bool_t isPrim;
    Bool_t isPhysPrim;
    Bool_t isPrimStack;
    Bool_t isPhysPrimStack;   
    Bool_t isChargedPrim;
    Int_t label;
    Int_t abslabel;    
    // Loop over particles
    for(int imc = 1; imc < (nTracksMC); imc++) {
        mcpart  = static_cast<AliMCParticle*>(fMCEvent->GetTrack(imc));
    if (!mcpart) { 
            Log("noMCParticle");
            continue;
        }        
         pt = mcpart->Pt();
        label           = mcpart->GetLabel();
        abslabel        = TMath::Abs(label);
        isPrim          = mcpart->IsPrimary();
        isPhysPrim      = fMCEvent->IsPhysicalPrimary(imc);        
        isPhysPrimStack = fMCStack->IsPhysicalPrimary(imc);
        Double_t charge = mcpart->Charge()/3.;
        isChargedPrim = kFALSE;
        if (isPhysPrim && ( TMath::Abs(charge) > 0.001 ) ) { isChargedPrim = kTRUE; }
        if (isPhysPrim != isPhysPrimStack) { Log("PhyPrimMismatch"); }        
        if (label != imc) { Log("LabelnotIndex"); }
        if (label != imc) { Log("NegativeLabel"); }
        if (isPrim) { nPrimCount++; }
        if (isPhysPrim) { nPhysPrim++; }        
        Double_t eta = mcpart->Eta();
        if (!isChargedPrim)  { continue; }
        if ( (TMath::Abs(eta) < 0.8) && (pt > 0.15) ) { fEventMCnPrim++; }
        if ( ( 2.8 < eta) && (eta <  5.1) ) { fEventMCnPrimV0M++; } //V0A
        if ( (-3.7 < eta) && (eta < -1.7) ) { fEventMCnPrimV0M++; } //V0C 
        
        // fill track histogram
        FillHist(fHistMCTrack1, pt);
    }        
    FillHist(fHistMCEvent1, fEventMCzv, fEventCent, fEventNtracks, fEventMultMB, fEventMultV0M, fEventMCb, fEventMCnPrim, fEventMCnPrimV0M, fEventIsTrigger, fEventHasVertex);    
    return kTRUE;
}

//_____________________________________________________________________________

void AliAnalysisTaskCutsDCA::FinishTaskOutput()
{
    // finish task output
    PostData(1, fOutputList); 
}

//_____________________________________________________________________________

void AliAnalysisTaskCutsDCA::Terminate(Option_t *)
{
    // terminate
}

//____________________________________________________________________________


Long64_t AliAnalysisTaskCutsDCA::FillHist(THnSparseD* s, Double_t x1, Double_t x2, Double_t x3, Double_t x4, Double_t x5, Double_t x6, 
                     Double_t x7 , Double_t x8 , Double_t x9 , Double_t x10 , Double_t x11 , Double_t x12 )
{
    if (s->GetNdimensions() > 12) { return 0; }
    Double_t vals[12]={x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12};
    return s->Fill(vals);
}

//____________________________________________________________________________

Int_t AliAnalysisTaskCutsDCA::AddAxis(const char* label, const char* title, Int_t nbins, Double_t xmin, Double_t xmax, const char* option)
{
    Int_t n=1;
    if (fSparseTmp) {
        n += fSparseTmp->GetNdimensions();
    }
    TString s;
    TArrayI bin(n);
    TArrayD min(n);
    TArrayD max(n);
    for (int i=0; i<n-1; i++) {
        bin[i] = fSparseTmp->GetAxis(i)->GetNbins();
        min[i] = fSparseTmp->GetAxis(i)->GetXmin();
        max[i] = fSparseTmp->GetAxis(i)->GetXmax();
        s += fSparseTmp->GetAxis(i)->GetName();
        s += ":";
    }
    bin[n-1] = nbins; 
    min[n-1] = xmin;
    max[n-1] = xmax;
    s += label;
    THnSparseD* h = new THnSparseD("fSparseTmp",s.Data(),n,bin.GetArray(),min.GetArray(),max.GetArray());
    for (int i=0; i<n-1; i++) {
        if (fSparseTmp->GetAxis(i)->GetXbins() && fSparseTmp->GetAxis(i)->GetXbins()->GetSize()) { h->SetBinEdges(i,fSparseTmp->GetAxis(i)->GetXbins()->GetArray()); }
        h->GetAxis(i)->SetTitle(fSparseTmp->GetAxis(i)->GetTitle());
        h->GetAxis(i)->SetName(fSparseTmp->GetAxis(i)->GetName());
    }
    h->GetAxis(n-1)->SetTitle(title);
    h->GetAxis(n-1)->SetName(label);
    if (fSparseTmp) { delete fSparseTmp; }
    fSparseTmp = h;
    return fSparseTmp->GetNdimensions();
}

//____________________________________________________________________________

Int_t AliAnalysisTaskCutsDCA::AddAxis(const char* label, Int_t nbins, Double_t xmin, Double_t xmax, const char* option)
{
    return AddAxis(label, label, nbins, xmin, xmax, option);
}

//____________________________________________________________________________

Int_t AliAnalysisTaskCutsDCA::AddAxis(const char* label, const char* title, Int_t nbins, Double_t* xbins, const char* option)
{
    Int_t n = AddAxis(label, title, nbins, xbins[0], xbins[nbins], option);
    fSparseTmp->SetBinEdges(n-1,xbins);         
    return n;
}

//____________________________________________________________________________

Int_t AliAnalysisTaskCutsDCA::AddAxis(const char* label, Int_t nbins, Double_t* xbins, const char* option)
{
    return AddAxis(label, label, nbins, xbins, option);
}

//____________________________________________________________________________

Int_t AliAnalysisTaskCutsDCA::AddAxis(const char* label, const char* title, const char* option)
{
    TString o(option);
    o.ToLower();
    if (o.Contains("ptfew")) {
        const Int_t nbins = 21;
        Double_t xbins[22] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0,  5.0, 10.0,  20.0, 50.0,  200.0};
        return AddAxis(label, title, nbins, xbins);    
    }
    if (o.Contains("ptveryfew")) {
        const Int_t nbins = 5;
        Double_t xbins[6] = {0.0, 0.5, 1.0, 1.5, 2.0, 200.0};
        return AddAxis(label, title, nbins, xbins);    
    }    
    if (o.Contains("pt")) {
        const Int_t nbins = 81;
        Double_t xbins[82] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 180.0, 200.0};
        return AddAxis(label, title, nbins, xbins);    
    }   
    if (o.Contains("cent")) {
        const Int_t nbins = 11;
        Double_t xbins[12] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.};
        return AddAxis(label, title, nbins, xbins);
    }    
    if (o.Contains("varsig35")) {
        const Int_t nbins = 35;
        Double_t xbins[36] = {-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,30,40,50,60,70,80,90,100,200,300,400,500,1000,2000};
        return AddAxis(label, title, nbins, xbins);
    }    
    return 0;
}

//____________________________________________________________________________

Int_t AliAnalysisTaskCutsDCA::AddAxis(const char* label, const char* option)
{
    return AddAxis(label, label, option);
}

//____________________________________________________________________________

Int_t AliAnalysisTaskCutsDCA::AddAxis(const char* option)
{
    TString o(option);
    o.ToLower();
    if (o.Contains("pt"))   return AddAxis("pt","p_{T} (GeV/c)",option);    
    if (o.Contains("cent"))   return AddAxis("cent","centrality",option);
    return 0;
}

//____________________________________________________________________________

THnSparseD* AliAnalysisTaskCutsDCA::CreateHist(const char* name)
{
    if (!fSparseTmp) return 0;
    THnSparseD* h = fSparseTmp;
    h->SetName(name);
    fSparseTmp = 0;    
    return h;
}

//____________________________________________________________________________

TH1D* AliAnalysisTaskCutsDCA::CreateLogHist(const char* name, const char* title)
{   
   TH1D *h = 0;
   if (title) { 
       h = new TH1D(name,title,200,0,200);        
    } else {
       h = new TH1D(name,name,200,0,200);
    }                   
   return h;
}   

//____________________________________________________________________________
  
TH1D* AliAnalysisTaskCutsDCA::CreateLogHist(const char* name) 
{ 
    return CreateLogHist(name,name);    
}



//____________________________________________________________________________
