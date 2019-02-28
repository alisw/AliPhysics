#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"


#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

#include "AliESDtrackCuts.h"
#include "AliESDpid.h"
#include "AliESDVertex.h"
#include "AliAnalysisFilter.h"

#include "AliAnalysisTaskTGReducedTree.h"
#include "AliESDpidCuts.h"
#include "AliVEvent.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"

#include "AliMCEvent.h"
#include "AliMCEventHandler.h"

#include "AliDielectronPID.h"
#include "AliDielectronVarManager.h"

#include "AliConversionPhotonCuts.h"
#include "AliV0ReaderV1.h"

#include "AliESDv0KineCuts.h"

// single electron tree to be used for dark photon analysis 
// Authors: Taku Gunji (CNS-Tokyo, Taku.Gunji@cern.ch)


ClassImp(AliDielectronTGReducedTrack)
ClassImp(AliDielectronTGReducedPair)
ClassImp(AliDielectronTGReducedInfo)
ClassImp(AliAnalysisTaskTGReducedTree)


//_______________________________________________________________
AliDielectronTGReducedTrack::AliDielectronTGReducedTrack() :
  fPx(-999),
  fPy(-999),
  fPz(-999),
  fPt(-999),
  fXv(-999),
  fYv(-999),
  fZv(-999),
  fPhi(-999),
  fTheta(-999),
  fCharge(-999),
  fImpactParXY(-999),
  fImpactParZ(-999),
  fITSsignal(-999),
  fITSnSigmaEleRaw(-999),
  fITSnSigmaEle(-999),
  fNclsITS(-999),
  fNclsSITS(-999),
  fNclsSFracITS(-999),
  fITSchi2Cl(-999),
  fPIn(-999),
  fPOut(-999),
  fTPCsignal(-999),
  fTOFsignal(-999),
  fTPCnSigmaEleRaw(-999),
  fTPCnSigmaEle(-999),
  fTOFnSigmaEle(-999),
  fGoldenChi2(-999),
  fTPCsignalN(-999),
  fTPCsignalNfrac(-999),
  fTPCchi2Cl(-999),
  fNclsTPC(-999),
  fNclsSTPC(-999),
  fNFclsTPC(-999),
  fNFclsTPCrFrac(-999),
  fNFclsTPCfCross(-999),
  fTrackStatus(-999),
  fSelectInfo(-999),
  fContVtx(-999),
  fPdgCode(-999),
  fMCPx(-999),
  fMCPy(-999),
  fMCPz(-999),
  fMCXv(-999),
  fMCYv(-999),
  fMCZv(-999),
  fPdgMother(-999),
  fMotherID(-999),
  fMotherLabel(-999),
  fGeneratorIndex(-999),
  fStatusAnotherLeg(-999),
  fLabelAnotherLeg(-999),
  fIDAnotherLeg(-999),
  fESDV0Conv(-999),
  fESDV0ConvID2(-999),
  fESDV0Mass(-999),
  fESDV0R(-999),
  fESDV0Phi(-999),
  fESDV0PsiPair(-999),
  fESDV0PhivPair(-999),
  fESDV0CosP(-999),
  fESDV0ArmAlpha(-999),
  fESDV0ArmPt(-999),
  fESDV0OA(-999),
  fESDV0cXv(-999),
  fESDV0cYv(-999),
  fESDV0cZv(-999),
  fESDV0Xv(-999),
  fESDV0Yv(-999),
  fESDV0Zv(-999),
  fV0Conv(-999),
  fV0ConvID2(-999),
  fV0Mass(-999),
  fV0Pt(-999),
  fV0Phi(-999),
  fV0Eta(-999),
  fV0DCAr(-999),
  fV0DCAz(-999),
  fV0X(-999),
  fV0Y(-999),
  fV0Z(-999),
  fV0R(-999),
  fV0Alpha(-999),
  fV0Qt(-999),
  fV0Psi(-999),
  fID(-999),
  fLabel(-999),
  fQualityFlags(-999)
{
}

//_______________________________________________________________
AliDielectronTGReducedTrack::AliDielectronTGReducedTrack(const AliDielectronTGReducedTrack &c):
  TObject(c),
  fPx(-999),
  fPy(-999),
  fPz(-999),
  fPt(-999),
  fXv(-999),
  fYv(-999),
  fZv(-999),
  fPhi(-999),
  fTheta(-999),
  fCharge(-999),
  fImpactParXY(-999),
  fImpactParZ(-999),
  fITSsignal(-999),
  fITSnSigmaEleRaw(-999),
  fITSnSigmaEle(-999),
  fNclsITS(-999),
  fNclsSITS(-999),
  fNclsSFracITS(-999),
  fITSchi2Cl(-999),
  fPIn(-999),
  fPOut(-999),
  fTPCsignal(-999),
  fTOFsignal(-999),
  fTPCnSigmaEleRaw(-999),
  fTPCnSigmaEle(-999),
  fTOFnSigmaEle(-999),
  fGoldenChi2(-999),
  fTPCsignalN(-999),
  fTPCsignalNfrac(-999),
  fTPCchi2Cl(-999),
  fNclsTPC(-999),
  fNclsSTPC(-999),
  fNFclsTPC(-999),
  fNFclsTPCrFrac(-999),
  fNFclsTPCfCross(-999),
  fTrackStatus(-999),
  fSelectInfo(-999),
  fContVtx(-999),
  fPdgCode(-999),
  fMCPx(-999),
  fMCPy(-999),
  fMCPz(-999),
  fMCXv(-999),
  fMCYv(-999),
  fMCZv(-999),
  fPdgMother(-999),
  fMotherID(-999),
  fMotherLabel(-999),
  fGeneratorIndex(-999),
  fStatusAnotherLeg(-999),
  fLabelAnotherLeg(-999),
  fIDAnotherLeg(-999),
  fESDV0Conv(-999),
  fESDV0ConvID2(-999),
  fESDV0Mass(-999),
  fESDV0R(-999),
  fESDV0Phi(-999),
  fESDV0PsiPair(-999),
  fESDV0PhivPair(-999),
  fESDV0CosP(-999),
  fESDV0ArmAlpha(-999),
  fESDV0ArmPt(-999),
  fESDV0OA(-999),
  fESDV0cXv(-999),
  fESDV0cYv(-999),
  fESDV0cZv(-999),
  fESDV0Xv(-999),
  fESDV0Yv(-999),
  fESDV0Zv(-999),
  fV0Conv(-999),
  fV0ConvID2(-999),
  fV0Mass(-999),
  fV0Pt(-999),
  fV0Phi(-999),
  fV0Eta(-999),
  fV0DCAr(-999),
  fV0DCAz(-999),
  fV0X(-999),
  fV0Y(-999),
  fV0Z(-999),
  fV0R(-999),
  fV0Alpha(-999),
  fV0Qt(-999),
  fV0Psi(-999),
  fID(-999),
  fLabel(-999),
  fQualityFlags(-999)
{
}

//_______________________________________________________________                  
AliDielectronTGReducedTrack::~AliDielectronTGReducedTrack()
{
  // Do nothing...
}



//_______________________________________________________________
AliDielectronTGReducedPair::AliDielectronTGReducedPair() : 
  fM(-999),
  fPx(-999),
  fPy(-999),
  fPz(-999),
  fPt(-999),
  fXv(-999),
  fYv(-999),
  fZv(-999),
  fPhiv(-999),
  fOpeningAngle(-999),
  fLegDistXY(-999),
  fLabel1(-999),
  fLabel2(-999),
  fIndex1(-999),
  fIndex2(-999),
  fC1(-999),
  fC2(-999)
{

}


//_______________________________________________________________
AliDielectronTGReducedPair::AliDielectronTGReducedPair(const AliDielectronTGReducedPair &c) :
  TObject(c),
  fM(-999),
  fPx(-999),
  fPy(-999),
  fPz(-999),
  fPt(-999),
  fXv(-999),
  fYv(-999),
  fZv(-999),
  fPhiv(-999),
  fOpeningAngle(-999),
  fLegDistXY(-999),
  fLabel1(-999),
  fLabel2(-999),
  fIndex1(-999),
  fIndex2(-999),
  fC1(-999),
  fC2(-999)
{
}


//_______________________________________________________________                  
AliDielectronTGReducedPair::~AliDielectronTGReducedPair()
{
  // Do nothing...
}



//_______________________________________________________________
AliDielectronTGReducedInfo::AliDielectronTGReducedInfo() :
  TObject(), fEventTag(-999), fRun(-1), fEvt(-1), fVtx(), fVtxCont(-1), fNe(-1), fV0Cand(-1), fTracks(0x0), fPairs(0x0)
{

  fVtx[0] = -1;
  fVtx[1] = -1; 
  fVtx[2] = -1;
  fTracks = new TClonesArray("AliDielectronTGReducedTrack", 100000);
  fPairs = new TClonesArray("AliDielectronTGReducedPair", 100000);
  
}

//_______________________________________________________________
AliDielectronTGReducedInfo::AliDielectronTGReducedInfo(const AliDielectronTGReducedInfo &c) :
  TObject(c),
  fEventTag(c.fEventTag),
  fRun(c.fRun), 
  fEvt(c.fEvt), 
  fVtxCont(c.fVtxCont), 
  fNe(c.fNe), 
  fV0Cand(c.fV0Cand),
  fTracks(c.fTracks),
  fPairs(c.fPairs)
{

  fVtx[0] = c.fVtx[0];
  fVtx[1] = c.fVtx[1];
  fVtx[2] = c.fVtx[2];
  fTracks = new TClonesArray("AliDielectronTGReducedTrack", 100000);
  fPairs = new TClonesArray("AliDielectronTGReducedPair", 100000);
}



//_______________________________________________________________  
AliDielectronTGReducedInfo::~AliDielectronTGReducedInfo()
{
  // Do nothing...

}


//_______________________________________________________________  
void AliDielectronTGReducedInfo::ClearEvent()
{
  if(fTracks) fTracks->Clear("C");
  if(fPairs) fPairs->Clear("C");
  fRun = 0;
  fEvt = 0;
  fVtxCont = 0;
  fNe = 0;
  fV0Cand = 0;
}


//_______________________________________________________________  
AliAnalysisTaskTGReducedTree::AliAnalysisTaskTGReducedTree(const char *name) 
  : AliAnalysisTaskSE(name), fESD(0), fOutputList(0), fTree(0), 
    fEventStat(0),
    fESDtrackCuts(0), 
    fTrackFilter(0), fMCEvent(0),
    fPIDResponse(0), 
    fPIDCuts(0), 
    fReducedInfo(0x0), fV0OpenCuts(0),
    fTriggerMask(AliVEvent::kMB),  
    fSelectPhysics(kFALSE), fFiredTrigger(""), fFiredExclude(kFALSE),
    fConvCut(117), 
    fCutArray(0), fMesonCutArray(0), fGammaCandidates(0),
    fV0Reader(NULL),fInputEvent(NULL),fReaderGammas(NULL),
    hasMC(kFALSE), fEvalEfficiencyFlag(kFALSE), 
    fEvalEfficiencyMCset(0),
    fEvalEfficiencyIndex(0),
    fEvalEfficiencyParticle(0),
  hGen(NULL), hMul(NULL)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
  //DefineOutput(1, AliDielectronTGReducedInfo::Class());
  //DefineOutput(1, TTree::Class());

  hPtRap[0] = NULL;  hPtRap[1] = NULL;  hPtRap[2] = NULL;
  hPtRapConv[0] = NULL;  hPtRapConv[1] = NULL;  hPtRapConv[2] = NULL;
  hRPtConv[0] = NULL;  hRPtConv[1] = NULL;  hRPtConv[2] = NULL;
}


//________________________________________________________________________
AliAnalysisTaskTGReducedTree::~AliAnalysisTaskTGReducedTree()
{
  //
  // Destructor
  //
  if(fEventStat)       { delete fEventStat;       fEventStat=0; }
  if(fTree)            { delete fTree;  fTree=0;}
} 

//________________________________________________________________________
void AliAnalysisTaskTGReducedTree::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if (man->GetMCtruthEventHandler()!=0x0) hasMC=kTRUE;


  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
  
  fReducedInfo = new AliDielectronTGReducedInfo();
  fTree = new TTree("DstTree","Reduced ESD information");

  fTree->Branch("Event",&fReducedInfo,16000,99);
  
  fEventStat=new TH1D("hEventStat", "Event statistics", 5, 0, 5);
  fEventStat->GetXaxis()->SetBinLabel(1,"Before Phys. Sel.");
  fEventStat->GetXaxis()->SetBinLabel(2,"After Phys. Sel.");
  fEventStat->GetXaxis()->SetBinLabel(3,"Bin3 not used");
  fEventStat->GetXaxis()->SetBinLabel(4,"Bin4 not used");
  fEventStat->GetXaxis()->SetBinLabel(5,"Bin5 not used");


  /// get V0 reader task
  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1");
  if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader  
  
  fOutputList->Add(fTree);
  fOutputList->Add(fEventStat);

  if(fV0Reader && (AliConversionPhotonCuts*)fV0Reader->GetConversionCuts()){
    if(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms()){
      fOutputList->Add(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
    }
  }

  //PostData(1, fOutputList);

  /// TrackFilter
  if(!fESDtrackCuts){
    fESDtrackCuts = new AliESDtrackCuts;
    fESDtrackCuts->SetPtRange( 0.2 , 100. );
    fESDtrackCuts->SetEtaRange( -0.8 , 0.8 );
    fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
    fESDtrackCuts->SetRequireSigmaToVertex(kFALSE);
    fESDtrackCuts->SetDCAToVertex2D(kFALSE);
    fESDtrackCuts->SetMaxDCAToVertexZ(3.);
    fESDtrackCuts->SetMaxDCAToVertexXY(1.);
    fESDtrackCuts->SetRequireITSRefit(kTRUE);
    fESDtrackCuts->SetRequireTPCRefit(kTRUE);
    
    fESDtrackCuts->SetMinNClustersITS(3); //5
    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    
    fESDtrackCuts->SetMaxChi2PerClusterITS(4.5);
    fESDtrackCuts->SetMinNClustersTPC(80); //100
    fESDtrackCuts->SetMinNCrossedRowsTPC(100); //100
    fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fESDtrackCuts->SetMaxChi2PerClusterTPC(4);
  }

  if(!fPIDCuts){
    fPIDCuts = new AliDielectronPID();
    fPIDCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3.5,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    fPIDCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    fPIDCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    fPIDCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5,4. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }

  fTrackFilter = new AliAnalysisFilter("trackFilter");
  fTrackFilter->AddCuts(fESDtrackCuts);


  if(hasMC && fEvalEfficiencyFlag){
    hPtRap[0] = new TH2F("hPtRapAll", "Pt and Rap for generated electrons (non-conv.)", 
			 100, 0, 10, 100, -2, 2);
    hPtRap[0]->SetXTitle("p_{T} (GeV/c)");  hPtRap[0]->SetYTitle("y");
    hPtRap[1] = new TH2F("hPtRapReco", "Pt and Rap for reconstructed electrons (non conv.)", 
			  100, 0, 10, 100, -2, 2);
    hPtRap[1]->SetXTitle("p_{T} (GeV/c)");  hPtRap[1]->SetYTitle("y");
    hPtRap[2] = new TH2F("hPtRapAcc", "Pt and Rap for accepted electrons (non conv.)", 
			 100, 0, 10, 100, -2, 2);
    hPtRap[2]->SetXTitle("p_{T} (GeV/c)");  hPtRap[2]->SetYTitle("y");


    hPtRapConv[0] = new TH2F("hPtRapConvAll", "Pt and Rap for generated electrons (conv.)", 
			 100, 0, 10, 100, -2, 2);
    hPtRapConv[0]->SetXTitle("p_{T} (GeV/c)");  hPtRapConv[0]->SetYTitle("y");
    hPtRapConv[1] = new TH2F("hPtRapConvReco", "Pt and Rap for reconstructed electrons (conv.)", 
			  100, 0, 10, 100, -2, 2);
    hPtRapConv[1]->SetXTitle("p_{T} (GeV/c)");  hPtRapConv[1]->SetYTitle("y");
    hPtRapConv[2] = new TH2F("hPtRapConvAcc", "Pt and Rap for accepted electrons (conv.)", 
			 100, 0, 10, 100, -2, 2);
    hPtRapConv[2]->SetXTitle("p_{T} (GeV/c)");  hPtRapConv[2]->SetYTitle("y");


    hRPtConv[0] = new TH2F("hRPtConvAll", "R and Pt for generated electrons from conversions", 
			   100, 0, 100, 100, -100, 100);
    hRPtConv[0]->SetXTitle("R [cm]"); hRPtConv[0]->SetYTitle("p_{T} (GeV/c)"); 

    hRPtConv[1] = new TH2F("hRPtConvReco", "R and Pt for reconstructed electrons from conversions", 
			   100, 0, 100, 100, -100, 100);
    hRPtConv[1]->SetXTitle("R [cm]"); hRPtConv[1]->SetYTitle("p_{T} (GeV/c)"); 

    hRPtConv[2] = new TH2F("hRPtConvAcc", "R and Pt for accepted electrons from conversions", 
			   100, 0, 100, 100, -100, 100);
    hRPtConv[2]->SetXTitle("R [cm]"); hRPtConv[2]->SetYTitle("p_{T} (GeV/c)"); 






    fOutputList->Add(hPtRap[0]); fOutputList->Add(hPtRap[1]); fOutputList->Add(hPtRap[2]);
    fOutputList->Add(hPtRapConv[0]); fOutputList->Add(hPtRapConv[1]); fOutputList->Add(hPtRapConv[2]);
    fOutputList->Add(hRPtConv[0]); fOutputList->Add(hRPtConv[1]); fOutputList->Add(hRPtConv[2]);
  }
  if(hasMC){
    hGen = new TH1F("hGen","hGen",10,0,10);
    hMul = new TH2F("hMul","hMul",10,0,10, 200, 0, 1000);
    fOutputList->Add(hGen);
    fOutputList->Add(hMul);
  }


  PostData(1, fOutputList);


}


//________________________________________________________________________
void AliAnalysisTaskTGReducedTree::UserExec(Option_t *) 
{

  // Main loop
  // Called for each event

  fReducedInfo->ClearEvent();

  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  if (!inputHandler) return;

  fInputEvent = InputEvent();

  ///check MC
  if(hasMC){
    fMCEvent = 0x0;
    AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!mcHandler){ 
      AliError("Could not retrive MC event handler!");
      return ; 
    }
    if (!mcHandler->InitOk() ){ return ;}
    if (!mcHandler->TreeK() ) { return ;}
    
    AliMCEvent* mcEvent = mcHandler->MCEvent();
    if (!mcEvent){ 
      AliError("Could not retrieve MC event!"); 
      return;
    }
    fMCEvent = mcEvent;
    fMCEvent->InitEvent();
    fMCEvent->PreReadAll();
  }

  if ( !(inputHandler->GetPIDResponse() )){
    AliFatal("This task needs the PID response attached to the input event handler!");
  }else{
    fPIDResponse = inputHandler->GetPIDResponse();
    //if(!hasMC){
    AliDielectronVarManager::SetPIDResponse(inputHandler->GetPIDResponse());
    //}
  }

  //Was event selected ?
  ULong64_t isSelected = AliVEvent::kAny;
  Bool_t isRejected = kFALSE;

  if( inputHandler){
    if(fSelectPhysics){
      if(inputHandler->GetEventSelection()){
	isSelected = inputHandler->IsEventSelected();
	isSelected &= fTriggerMask;
      }
    }

    if( isSelected && !fFiredTrigger.IsNull() ){
      TString firedTriggerClasses=fInputEvent->GetFiredTriggerClasses();
      isSelected=(firedTriggerClasses.Contains(fFiredTrigger))^fFiredExclude;
    }
  }

  fEventStat->Fill(0);
  /// physics selection   
  if (isSelected==0){
    return;
  }
  fEventStat->Fill(1);

  // Post output data.
  fESD = dynamic_cast<AliESDEvent*>(fInputEvent);
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }
  
  printf("There are %d tracks in this event... Let's go\n", fESD->GetNumberOfTracks());

  const AliESDVertex* vtx = fESD->GetPrimaryVertex();

  const AliESDVertex* vertex = 0;
  vertex = fESD->GetPrimaryVertexTracks();
  if(!vertex)
    vertex = fESD->GetPrimaryVertexSPD();
  if(!vertex)
    vertex = fESD->GetPrimaryVertexTPC();


  Int_t fNe = 0;

  TClonesArray& tracks = *(fReducedInfo->fTracks);
  AliDielectronTGReducedTrack *rtrack = NULL;

  // Track loop
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
    }

    Double_t values[AliDielectronVarManager::kNMaxValues];
    AliDielectronVarManager::Fill(track, values);

    Int_t selectInfo = fTrackFilter->IsSelected(track);
    if(!selectInfo) continue; 

    Bool_t selectPID = kFALSE;
    selectPID = fPIDCuts->IsSelected(track);

    Bool_t isElectron = kTRUE;
    if(hasMC){
      if (!fMCEvent) {
        AliError(Form("Could not receive MC -> hasMC set to kFALSE!!"));
        continue;
      }else{
        AliMCParticle* mcTrack = dynamic_cast<AliMCParticle *>
          (fMCEvent->GetTrack(TMath::Abs(track->GetLabel())));
        if(TMath::Abs(mcTrack->PdgCode())!=11){
          isElectron = kFALSE;
        }
      }
    }else{
      isElectron = selectPID;
    }

    selectInfo = selectInfo | (selectPID <<1) | (isElectron << 2); 
    
    ///final cut to fill into tree
    if(selectInfo<3)continue;


    rtrack = new(tracks[fNe]) AliDielectronTGReducedTrack();
    //rtrack = new AliDielectronTGReducedTrack();


    rtrack->ID(iTracks); //track->GetID());
    rtrack->Label(track->GetLabel());
    rtrack->Px(values[AliDielectronVarManager::kPx]);
    rtrack->Py(values[AliDielectronVarManager::kPy]);
    rtrack->Pz(values[AliDielectronVarManager::kPz]);
    rtrack->Pt(values[AliDielectronVarManager::kPt]);
    rtrack->Xv(values[AliDielectronVarManager::kXv]);
    rtrack->Yv(values[AliDielectronVarManager::kYv]);
    rtrack->Zv(values[AliDielectronVarManager::kZv]);
    rtrack->Phi(values[AliDielectronVarManager::kPhi]);
    rtrack->Theta(values[AliDielectronVarManager::kTheta]);
    rtrack->Charge(values[AliDielectronVarManager::kCharge]);
    rtrack->ImpactParXY(values[AliDielectronVarManager::kImpactParXY]);
    rtrack->ImpactParZ(values[AliDielectronVarManager::kImpactParZ]);
    rtrack->ITSsignal(values[AliDielectronVarManager::kITSsignal]);
    rtrack->ITSnSigmaEleRaw(values[AliDielectronVarManager::kITSnSigmaEleRaw]);
    rtrack->ITSnSigmaEle(values[AliDielectronVarManager::kITSnSigmaEle]);
    rtrack->NclsITS(values[AliDielectronVarManager::kNclsITS]);
    rtrack->NclsSITS(values[AliDielectronVarManager::kNclsSITS]);
    rtrack->NclsSFracITS(values[AliDielectronVarManager::kNclsSFracITS]);
    rtrack->ITSchi2Cl(values[AliDielectronVarManager::kITSchi2Cl]);
    rtrack->PIn(values[AliDielectronVarManager::kPIn]);
    rtrack->POut(values[AliDielectronVarManager::kPOut]);
    rtrack->TPCsignal(values[AliDielectronVarManager::kTPCsignal]);
    rtrack->TOFsignal(values[AliDielectronVarManager::kTOFsignal]);
    rtrack->TPCnSigmaEleRaw(values[AliDielectronVarManager::kTPCnSigmaEleRaw]);
    rtrack->TPCnSigmaEle(values[AliDielectronVarManager::kTPCnSigmaEle]);
    rtrack->TOFnSigmaEle(values[AliDielectronVarManager::kTOFnSigmaEle]);
    rtrack->GoldenChi2(track->GetChi2TPCConstrainedVsGlobal(vertex));
    rtrack->TPCsignalN(values[AliDielectronVarManager::kTPCsignalN]);
    rtrack->TPCsignalNfrac(values[AliDielectronVarManager::kTPCsignalNfrac]);
    rtrack->TPCchi2Cl(values[AliDielectronVarManager::kTPCchi2Cl]);
    rtrack->NclsTPC(values[AliDielectronVarManager::kNclsTPC]);
    rtrack->NclsSTPC(values[AliDielectronVarManager::kNclsSTPC]);
    rtrack->NFclsTPC(values[AliDielectronVarManager::kNFclsTPC]);
    rtrack->NFclsTPCrFrac(values[AliDielectronVarManager::kNFclsTPCrFrac]);
    rtrack->NFclsTPCfCross(values[AliDielectronVarManager::kNFclsTPCfCross]);
    rtrack->TrackStatus(values[AliDielectronVarManager::kTrackStatus]);
    rtrack->SelectInfo(selectInfo);
    rtrack->ContVtx(vtx->UsesTrack(track->GetID()));

    if(hasMC){
       AliMCParticle* mcTrack = dynamic_cast<AliMCParticle *>
	 (fMCEvent->GetTrack(TMath::Abs(track->GetLabel())));

       rtrack->PdgCode(mcTrack->PdgCode());
       rtrack->MCPx(mcTrack->Px());
       rtrack->MCPy(mcTrack->Py());
       rtrack->MCPz(mcTrack->Pz());
       rtrack->MCXv(mcTrack->Xv());
       rtrack->MCYv(mcTrack->Yv());
       rtrack->MCZv(mcTrack->Zv());
       rtrack->GeneratorIndex(mcTrack->GetGeneratorIndex());
       

       if(!(mcTrack->GetMother() < 0)) {
	 AliMCParticle* mcmother = dynamic_cast<AliMCParticle *>
	   (fMCEvent->GetTrack(TMath::Abs(mcTrack->GetMother())));
	 rtrack->MotherID(mcTrack->GetMother());
	 rtrack->PdgMother(mcmother->PdgCode());
	 rtrack->MotherLabel(mcmother->GetLabel());

	 Int_t child1 = mcmother->GetDaughterFirst();
	 Int_t child2 = mcmother->GetDaughterLast();

	 AliMCParticle *mcTrackD[2] = {NULL, NULL};
	 Int_t index=0;
	 for(int ichild=child1; ichild<=child2; ichild++){
	   AliMCParticle *mcTrackTmp =  dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(ichild));
	   if(TMath::Abs(mcTrackTmp->PdgCode())==11){
	     mcTrackD[index] =  mcTrackTmp;
	     index++;
	   }
	 }

	 if(index==2){ //electron and positron pairs from this mother
	   //// get reconstructed information
	   AliMCParticle *mcTrackDtmp =  NULL;
	   if(TMath::Abs(track->GetLabel())==TMath::Abs(mcTrackD[0]->GetLabel())){
	     mcTrackDtmp = mcTrackD[1];
	   }else if(TMath::Abs(track->GetLabel())==TMath::Abs(mcTrackD[1]->GetLabel())){
	     mcTrackDtmp = mcTrackD[0];
	   }

	   ///loop over ESD track again...
	   for (Int_t itr = 0; itr < fESD->GetNumberOfTracks(); itr++) {
	     if(itr==iTracks) continue;
	     AliESDtrack* trackD = fESD->GetTrack(itr);
	     Int_t pairFlag = -99; 
	     if(TMath::Abs(trackD->GetLabel())==TMath::Abs(mcTrackDtmp->GetLabel())){
	       pairFlag = fTrackFilter->IsSelected(trackD);	     
	       rtrack->StatusAnotherLeg(pairFlag);
	       rtrack->LabelAnotherLeg(TMath::Abs(trackD->GetLabel()));
	       rtrack->IDAnotherLeg(TMath::Abs(trackD->GetID()));
	       break;
	     }
	   }//end of itr loop
	 }
       }//end of looking at mothers
    }//has MC
    fNe++;
  }// end of track loop
  

  
  ///process for pairs ///
  TClonesArray& pairs = *(fReducedInfo->fPairs);
  AliDielectronTGReducedPair *rpair = NULL;  
  Int_t fNpair = 0;

  for(Int_t ine=0; ine<fNe; ine++){
    for(Int_t jne=ine+1; jne<fNe; jne++){
      AliDielectronTGReducedTrack *trk1 = 
	(AliDielectronTGReducedTrack*)fReducedInfo->GetTrack(ine);
      AliDielectronTGReducedTrack *trk2 = 
	(AliDielectronTGReducedTrack*)fReducedInfo->GetTrack(jne);


      AliESDtrack *tr1 = fESD->GetTrack(trk1->ID());
      AliESDtrack *tr2 = fESD->GetTrack(trk2->ID());

      AliKFParticle fPair ;
      fPair.Initialize();
      AliKFParticle kf1(*tr1, 11*trk1->Charge());
      AliKFParticle kf2(*tr2, 11*trk2->Charge());
      fPair.AddDaughter(kf1);
      fPair.AddDaughter(kf2);
      
      rpair = new(pairs[fNpair]) AliDielectronTGReducedPair();

      rpair->M(fPair.GetMass());
      rpair->Px(fPair.GetPx());
      rpair->Py(fPair.GetPy());
      rpair->Pz(fPair.GetPz());
      rpair->Xv(fPair.GetX());
      rpair->Yv(fPair.GetY());
      rpair->Zv(fPair.GetZ());
      rpair->Label1(trk1->ID());
      rpair->Label2(trk2->ID());
      rpair->Index1(ine);
      rpair->Index2(jne);
      rpair->C1(trk1->Charge());
      rpair->C2(trk2->Charge());
      rpair->OpeningAngle(kf1.GetAngleXY(kf2));
      rpair->Phiv(kf1.GetAngle(kf2)); //temporary
      rpair->LegDistXY(kf1.GetAngleRZ(kf2)); //temporary

      fNpair++;
      
    }
  }

  ////////////////////////
  /// getting conversion information from on the fly ALICE V0 finder          
  /// by reading ESD            
  ///do debug for conversion    

  if(fV0OpenCuts){
    fV0OpenCuts->SetEvent(fESD);
    AliKFVertex primaryVertexKF(*vtx);
    fV0OpenCuts->SetPrimaryVertex(&primaryVertexKF);
  }    


  for(Int_t currentV0Index=0; currentV0Index<fESD->GetNumberOfV0s(); currentV0Index++){
    AliESDv0 *fCurrentV0=(AliESDv0*)(fESD->GetV0(currentV0Index));
    if(!fCurrentV0){
      printf("Requested V0 does not exist");
      continue;
    }
    // look at on the fly 
    if(!fCurrentV0->GetOnFlyStatus()) continue;

    AliESDtrack *trNeg = fESD->GetTrack(fCurrentV0->GetIndex(0));
    AliESDtrack *trPos = fESD->GetTrack(fCurrentV0->GetIndex(1));
    // protection against LS v0s
    if(trNeg->Charge() == trPos->Charge()) continue;
    
    Bool_t v0ChargesAreCorrect = (trPos->GetSign()==+1 ? kTRUE : kFALSE);
    trPos = (!v0ChargesAreCorrect ? fESD->GetTrack(fCurrentV0->GetNindex()) : trPos);
    trNeg = (!v0ChargesAreCorrect ? fESD->GetTrack(fCurrentV0->GetPindex()) : trNeg); 

    
   
    AliKFParticle fPair ;
    fPair.Initialize();
    AliKFParticle kf1(*trNeg, -11);
    AliKFParticle kf2(*trPos,  11);
    //fPair.ConstructGamma(kf1,kf2); 
    fPair.AddDaughter(kf1);
    fPair.AddDaughter(kf2);

    AliDielectronPair candidate;
    candidate.SetPdgCode(22); //assuming from photon
    candidate.SetGammaTracks(trNeg, 11, trPos, 11);

    // TrackLabels
    Int_t currentTrackLabels[2]={-1,-1};
    // Get Daughter KF Particles 
    const AliExternalTrackParam *fCurrentExternalTrackParamPositive = fCurrentV0->GetParamP();
    currentTrackLabels[0] = fCurrentV0->GetPindex(); //! this can associate to the ESDtrack (track->GetID()) 
    const AliExternalTrackParam *fCurrentExternalTrackParamNegative = fCurrentV0->GetParamN();
    currentTrackLabels[1] = fCurrentV0->GetNindex(); //! this can associate to the ESDtrack (track->GetID())
    
    currentTrackLabels[0] = trPos->GetID();
    currentTrackLabels[1] = trNeg->GetID();


    ////judge whether this is K0, Lambda, or Gamma
    Int_t pdgV0=0; Int_t pdgP=0; Int_t pdgN=0;
    Bool_t processV0;
    if(fV0OpenCuts){
      processV0 = fV0OpenCuts->ProcessV0(fCurrentV0, pdgV0, pdgP, pdgN);
    }


    ////associate with track       
    /// proper tagging is lost if one electron join 2 photon pairs   
    for (Int_t iTracks = 0; iTracks < fNe; iTracks++) {
      rtrack = fReducedInfo->GetTrack(iTracks);

      if(rtrack->ID() == currentTrackLabels[0] ||
	 rtrack->ID() == currentTrackLabels[1] ){

	if(rtrack->ID() == currentTrackLabels[0]){
	  rtrack->ESDV0Conv(0);
	  rtrack->ESDV0ConvID2(currentTrackLabels[1]);
	}else if(rtrack->ID() == currentTrackLabels[1]){
	  rtrack->ESDV0Conv(1);
	  rtrack->ESDV0ConvID2(currentTrackLabels[0]);
	}

	rtrack->ESDV0Mass(fPair.GetMass());
	rtrack->ESDV0R(fPair.GetR());
	rtrack->ESDV0Phi(fPair.GetPhi());
	rtrack->fQualityFlags = (processV0 ? pdgV0 : -999);
	
	rtrack->ESDV0PsiPair(
			     fInputEvent ? 
			     candidate.PsiPair(fInputEvent->GetMagneticField()) : -5);
	rtrack->ESDV0PhivPair(
			      fInputEvent ? 
			      candidate.PhivPair(fInputEvent->GetMagneticField()) : -5);
	rtrack->ESDV0CosP(candidate.GetCosPointingAngle(vtx));
	rtrack->ESDV0ArmAlpha(candidate.GetArmAlpha());
	rtrack->ESDV0ArmPt(candidate.GetArmPt());
	rtrack->ESDV0OA(candidate.OpeningAngle());
	rtrack->ESDV0cXv(candidate.Xv());
	rtrack->ESDV0cYv(candidate.Yv());
	rtrack->ESDV0cZv(candidate.Zv());

	rtrack->ESDV0Xv(fCurrentV0->Xv());
	rtrack->ESDV0Yv(fCurrentV0->Yv());
	rtrack->ESDV0Zv(fCurrentV0->Zv());


      }

    }// end of iTracks loop        

    


  } //end of V0 loop
	 

  
  /// by reading using AliV0ReaderV1. this has to be defined in main macro.
  fReaderGammas = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut  
  // Loop over Photon Candidates allocated by ReaderV1       
  for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
    if(!PhotonCandidate) continue;
    //if(!((AliConversionPhotonCuts*)fCutArray->At(0))->PhotonIsSelected(PhotonCandidate,fInputEvent)) continue;
    //fGammaCandidates->Add(PhotonCandidate); 


    for (Int_t iTracks = 0; iTracks < fNe; iTracks++) {
      rtrack = fReducedInfo->GetTrack(iTracks);
      if(rtrack->ID() == PhotonCandidate->GetLabel1() ||
	 rtrack->ID() ==  PhotonCandidate->GetLabel2()){
	
	if(rtrack->ID() == PhotonCandidate->GetLabel1()){
	  rtrack->V0Conv(0);
	  rtrack->V0ConvID2(PhotonCandidate->GetLabel2());
	}else if(rtrack->ID() == PhotonCandidate->GetLabel2()){
	  rtrack->V0Conv(1);
	  rtrack->V0ConvID2(PhotonCandidate->GetLabel1());
	}

	/// proper tagging is lost if one electron join 2 photon pairs    
	rtrack->V0Mass(PhotonCandidate->GetInvMassPair());
        rtrack->V0Pt(PhotonCandidate->GetPhotonPt());
        rtrack->V0Phi(PhotonCandidate->GetPhotonPhi());
        rtrack->V0Eta(PhotonCandidate->GetPhotonEta());
        rtrack->V0DCAr(PhotonCandidate->GetDCArToPrimVtx());
        rtrack->V0DCAz(PhotonCandidate->GetDCAzToPrimVtx());
        rtrack->V0X(PhotonCandidate->GetConversionX());
        rtrack->V0Y(PhotonCandidate->GetConversionY());
        rtrack->V0Z(PhotonCandidate->GetConversionZ());
        rtrack->V0R(PhotonCandidate->GetConversionRadius());
        rtrack->V0Alpha(PhotonCandidate->GetArmenterosAlpha());
        rtrack->V0Qt(PhotonCandidate->GetArmenterosQt());
        rtrack->V0Psi(PhotonCandidate->GetPsiPair());
      }

    }
  }



  Int_t isThisMC = 0;

  if(hasMC){
    ////get generator header
    TList *l = (TList*)fMCEvent->GetCocktailList();
    //cout<<" GetCocktailList entries = "<<l->GetEntries()<<endl;
    for (Int_t i = l->GetEntries()-1; i >= 0; i--){
      AliGenEventHeader* gh=(AliGenEventHeader*)l->At(i);
      TString genname=gh->GetName();

      printf("Generator Index %d - %s\n", i, genname.Data());

      if(genname.Contains("Pythia MB")){
	hGen->Fill(0., 1.0); 
	hMul->Fill(0., 1.0*gh->NProduced());
      }else if(genname.Contains("Pythia CC")){
	hGen->Fill(1., 1.0); 
	hMul->Fill(1,  1.0*gh->NProduced());
      }else if(genname.Contains("Pythia BB")){
	hGen->Fill(2., 1.0); 
	hMul->Fill(2., 1.0*gh->NProduced());
      }else if(genname.Contains("Pythia B")){
	hGen->Fill(3., 1.0); 
	hMul->Fill(3., 1.0*gh->NProduced());
      }else if(genname.Contains("Jpsi2ee")){
	hGen->Fill(4., 1.0); 
	hMul->Fill(4., 1.0*gh->NProduced());
      }else if(genname.Contains("B2Jpsi2ee")){
	hGen->Fill(5., 1.0);
	hMul->Fill(5., 1.0*gh->NProduced());
      }

      if(genname.Contains(fEvalEfficiencyMCset)){
	isThisMC = 1;
      }
    }
    
    //AliStack *stack = (AliStack*)fMCEvent->Stack();
    //cout<<"Primary summary "<<stack->GetNtrack()<<" "<<stack->GetNprimary()<<" "<<fMCEvent->GetNumberOfTracks()<<" "<<fMCEvent->GetNumberOfPrimaries()<<" "<<tot<<endl;

    if(fEvalEfficiencyFlag && isThisMC){ 
      if(fMCEvent){
	for(int iTrack=0; iTrack<fMCEvent->GetNumberOfTracks(); iTrack++){
	  AliMCParticle *part = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(iTrack));	
	  AliMCParticle *mother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(part->GetMother())));	
	  
	  Int_t genOK = -1;
	  Int_t recoHF=-1;
	  Int_t recoConv=-1;

	  if(mother && TMath::Abs(part->PdgCode())==11){

	    //cout<<" Event "<<fESD->GetEventNumberInFile()<<" "<<iTrack<<" "<<part->PdgCode()<<" "<<part->GetGeneratorIndex()<<" mother PDG="<<mother->PdgCode()<<" "<<gh->GetName()<<" "<<part->Pt()<<" "<<part->Y()<<endl;

	    //select specific event generator for the efficiency evaluation
	    if(fEvalEfficiencyIndex && part->GetGeneratorIndex()!=fEvalEfficiencyIndex) continue;
	  
	    ///select mother particles
	    if(fEvalEfficiencyParticle == 0){ 
	      ///all electrons except conversions
	      genOK = -1;
	      if(TMath::Abs(mother->PdgCode()) != 22) genOK = 1;
	    }else if(fEvalEfficiencyParticle == 4){ 
	      /// all electrons from charm 
	      genOK = -1;
	      if(TMath::Abs(mother->PdgCode())/100 == 4 ||
		 TMath::Abs(mother->PdgCode())/1000 == 4) genOK = 1;
	    }else if(fEvalEfficiencyParticle == 5){
	      /// all electrons from bottom 
	      genOK = -1;
	      if(TMath::Abs(mother->PdgCode())/100 == 5||
                 TMath::Abs(mother->PdgCode())/1000 == 5) genOK = 1;
	    }else if(fEvalEfficiencyParticle == 45){
	      /// all electrons from heavy flavors
	      genOK = -1;
	      if(TMath::Abs(mother->PdgCode())/100 == 4 ||
		 TMath::Abs(mother->PdgCode())/100 == 5 ||
                 TMath::Abs(mother->PdgCode())/1000 == 4 || 
                 TMath::Abs(mother->PdgCode())/1000 == 5 
		 ) genOK = 1;
	    }else{
	      genOK = -1;
	    }

	    if(genOK==1){
	      for (Int_t itr = 0; itr < fESD->GetNumberOfTracks(); itr++) {
		AliESDtrack* reco = fESD->GetTrack(itr);
		if(TMath::Abs(reco->GetLabel())==TMath::Abs(part->GetLabel())){
		  recoHF = fTrackFilter->IsSelected(reco);	     
		  break;
		}
	      }
	    
	      hPtRap[0]->Fill(part->Pt(), part->Y());
	      if(recoHF>=0) hPtRap[1]->Fill(part->Pt(), part->Y());
	      if(recoHF>0) hPtRap[2]->Fill(part->Pt(), part->Y());
	    }

	  ///select Conversions for efficiency calculations 
	  ///Some rapidity cuts must be applied in order to get the conversions from BP, ITS, TPC...
	  ///Reco conversions would come outside |Zv|<60 (TPC bessel). 
	  ///But don't care since what I'd like to know is the 
	  /// tagging efficiency for ITS.
	    if(TMath::Abs(mother->PdgCode()) == 22){
	      AliVParticle *grandmother = (AliVParticle*)fMCEvent->GetTrack(TMath::Abs(mother->GetMother()));	
	      
	      if(grandmother && TMath::Abs(grandmother->PdgCode()) == 111 &&
		 TMath::Abs(grandmother->Y())<0.8){
		
		///look for reconstructed pairs
		for (Int_t itr = 0; itr < fESD->GetNumberOfTracks(); itr++) {
		  AliESDtrack* reco = fESD->GetTrack(itr);
		if(TMath::Abs(reco->GetLabel())==TMath::Abs(part->GetLabel())){
		  recoConv = fTrackFilter->IsSelected(reco);	     
		  break;
		}
		}
		
		hPtRapConv[0]->Fill(part->Pt(), part->Y());
		if(recoConv>=0) hPtRapConv[1]->Fill(part->Pt(), part->Y());
		if(recoConv>0) hPtRapConv[2]->Fill(part->Pt(), part->Y());
		
		hRPtConv[0]->Fill(sqrt(pow(part->Xv(),2)+pow(part->Yv(),2)), part->Zv());
		if(recoConv>=0) hRPtConv[1]->Fill(sqrt(pow(part->Xv(),2)+pow(part->Yv(),2)), part->Zv());
		if(recoConv>0) hRPtConv[2]->Fill(sqrt(pow(part->Xv(),2)+pow(part->Yv(),2)), part->Zv());
	      } //require that grandmother is pi0 and good acceptance
	    } //conversion electrons
	  }
	} //MC track look
      }// fMCEvent
    }//hasMC
  }

  fReducedInfo->fRun = fESD->GetRunNumber();
  fReducedInfo->fEvt = fESD->GetEventNumberInFile();
  fReducedInfo->fNe = fNe;
  fReducedInfo->fVtx[0] = vtx->GetX();
  fReducedInfo->fVtx[1] = vtx->GetY();
  fReducedInfo->fVtx[2] = vtx->GetZ();
  fReducedInfo->fVtxCont = vtx->GetNContributors();
  fReducedInfo->fV0Cand = fReaderGammas->GetEntriesFast();

  fTree->Fill();
  PostData(1, fOutputList);

}


//________________________________________________________________________
void AliAnalysisTaskTGReducedTree::Terminate(Option_t *) 
{
}


//________________________________________________________________________
void AliAnalysisTaskTGReducedTree::SetGammaConvCut(TString photonCut)
{

  /*
    cuts.AddCut("00010113", "00200009227302008250400000", "0152101500000000");
    cuts.AddCut("00010113", "00200009227302008250404000", "0152101500000000");
    cuts.AddCut("00010113", "00200069227302008250400000", "0152101500000000");
    cuts.AddCut("00010113", "00200069227302008250404000", "0152101500000000");
    cuts.AddCut("00010113", "00200019227302008250400000", "0152101500000000");
    cuts.AddCut("00010113", "00200019227302008250404000", "0152101500000000");
    cuts.AddCut("00010113", "00200079227302008250400000", "0152101500000000");
    cuts.AddCut("00010113", "00200079227302008250404000", "0152101500000000");
  */

  fCutArray = new TList();
  fCutArray->SetOwner(kTRUE);
  
  AliConversionPhotonCuts *analysisCuts = new AliConversionPhotonCuts();
  //difference from above 
  // minR (2->0) : R>5cm --> R>0cm
  // single pT (0->6) : pT>50 MeV --> pT>40MeV
  // TOF PID (2->0) : -5<sigma<5 -> no TOF electron PID
  //TString photonCut = "00000069227300008250400000";
  analysisCuts->InitializeCutsFromCutString(photonCut.Data());
  analysisCuts->PrintCuts();
  fCutArray->Add(analysisCuts);
    
}






