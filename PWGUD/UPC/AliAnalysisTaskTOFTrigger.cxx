/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// c++ headers
#include <iostream>
#include <string.h>

// root headers
#include "TROOT.h"
#include "TFile.h"
#include "TColor.h"
#include "TEfficiency.h"
#include <TGeoMatrix.h>

// aliroot headers
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDEvent.h"
#include "AliESDVZERO.h"
#include "AliESDZDC.h"
#include "AliPIDResponse.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliTOFTriggerMask.h"
#include "AliAnalysisUtils.h"
#include "AliGeomManager.h"
#include "AliTrackerBase.h"
#include "AliExternalTrackParam.h"
#include "AliTOFGeometry.h"
#include "AliESDtrackCuts.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"


// my headers
#include "AliAnalysisTaskTOFTrigger.h"

ClassImp(AliAnalysisTaskTOFTrigger);

using std::cout;
using std::endl;

//analysis skeleton of UPC nano AODs,

//_____________________________________________________________________________
AliAnalysisTaskTOFTrigger::AliAnalysisTaskTOFTrigger()
  : AliAnalysisTaskSE(),fOutputList(0),fPIDResponse(0),fTrackCuts(0),
	fTOFmask(0),
	eff_MaxiPadLTM_All(0),
	eff_MaxiPadLTM_Mu(0),
	eff_MaxiPadLTM_El(0),
	eff_MaxiPadLTM_1Trk_All(0),
	eff_MaxiPadLTM_1Trk_Mu(0),
	eff_MaxiPadLTM_1Trk_El(0),
	eff_AverageTracklets(0),
	eff_AverageTrackPt(0),
	eff_MaxiPadLTM_Around(0),
	eff_MaxiPadLTM_OnlyAround(0),
	eff_MaxiPadLTM_Clusters(0),
	hTrackDistributionLTM(0),
	hTrackDistribution_Mu(0),
	hTrackDistribution_El(0),
	hTrackDistribution(0),
	hFiredMaxiPad(0),
	hFiredMaxiPadOnlyAround(0),
	hNotFiredMaxiPadCls(0),
	hExtraFiredMaxiPadCls(0),
	hNotFiredMaxiPadTrk(0),
	hExtraFiredMaxiPadTrk(0),
	hTrackPadCorrPhi(0),
	hTrackPadCorrEta(0),
	hNoiseMaxiPad(0),
	hTriggerCounter(0),
	hTriggerCounterIR1(0),
	hTriggerCounterIR2(0),
	hNFiredMaxiPads(0),
	hNFiredMaxiPadsOnlyAround(0),
	hNTracklets(0),
	hDetIn0(0),
	hDetIn1(0),
	hDetIn2(0),
	hDetIn3(0),
	hDetIn4(0),
	hPadDistance(0),
	hTrackPt(0),
	hNMaxiPadIn(0),
	hNCrossTracks(0),
	hBadMaxiPadMask(0),
	hTOFHitTime(0),
	fGeomLoaded(kFALSE),
	fMaxPt(0),
	fMinPt(0),
	fMaxMulti(0),
	fTriggerClass(0),
	fMaxBCs(0),
	fUseEventSelection(0),
	fTrackCutSet(0),
	fMaxTrackError(0),
	fMinTOF(0),
	fMaxTOF(0),
	fEventCuts(0)
{

//Dummy constructor

}//AliAnalysisTaskTOFTrigger


//_____________________________________________________________________________
AliAnalysisTaskTOFTrigger::AliAnalysisTaskTOFTrigger(const char *name,Float_t lowpt,Float_t highpt,Int_t highmult,TString trgcls,Int_t nBCs,Bool_t useEVS,Int_t cutSet,Float_t maxErr,Float_t mintof,Float_t maxtof)
  : AliAnalysisTaskSE(name),fOutputList(0),fPIDResponse(0),fTrackCuts(0),
	fTOFmask(0),
	eff_MaxiPadLTM_All(0),
	eff_MaxiPadLTM_Mu(0),
	eff_MaxiPadLTM_El(0),
	eff_MaxiPadLTM_1Trk_All(0),
	eff_MaxiPadLTM_1Trk_Mu(0),
	eff_MaxiPadLTM_1Trk_El(0),
	eff_AverageTracklets(0),
	eff_AverageTrackPt(0),
	eff_MaxiPadLTM_Around(0),
	eff_MaxiPadLTM_OnlyAround(0),
	eff_MaxiPadLTM_Clusters(0),
	hTrackDistributionLTM(0),
	hTrackDistribution_Mu(0),
	hTrackDistribution_El(0),
	hTrackDistribution(0),
	hFiredMaxiPad(0),
	hFiredMaxiPadOnlyAround(0),
	hNotFiredMaxiPadCls(0),
	hExtraFiredMaxiPadCls(0),
	hNotFiredMaxiPadTrk(0),
	hExtraFiredMaxiPadTrk(0),
	hTrackPadCorrPhi(0),
	hTrackPadCorrEta(0),
	hNoiseMaxiPad(0),
	hTriggerCounter(0),
	hTriggerCounterIR1(0),
	hTriggerCounterIR2(0),
	hNFiredMaxiPads(0),
	hNFiredMaxiPadsOnlyAround(0),
	hNTracklets(0),
	hDetIn0(0),
	hDetIn1(0),
	hDetIn2(0),
	hDetIn3(0),
	hDetIn4(0),
	hPadDistance(0),
	hTrackPt(0),
	hNMaxiPadIn(0),
	hNCrossTracks(0),
	hBadMaxiPadMask(0),
	hTOFHitTime(0),
	fGeomLoaded(kFALSE),
	fMaxPt(highpt),
	fMinPt(lowpt),
	fMaxMulti(highmult),
	fTriggerClass(trgcls),
	fMaxBCs(nBCs),
	fUseEventSelection(useEVS),
	fTrackCutSet(cutSet),
	fMaxTrackError(maxErr),
	fMinTOF(mintof),
	fMaxTOF(maxtof),
	fEventCuts(0)
{

  DefineOutput(1, TList::Class());

}//AliAnalysisTaskTOFTrigger

//_____________________________________________________________________________
AliAnalysisTaskTOFTrigger::~AliAnalysisTaskTOFTrigger()
{
  // Destructor

  // Destructor
  if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis){
     delete fOutputList;
     fOutputList = 0x0;
  }

}//~AliAnalysisTaskTOFTrigger


//_____________________________________________________________________________
void AliAnalysisTaskTOFTrigger::UserCreateOutputObjects()
{

  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  
  if(fTrackCutSet == 1)fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
  if(fTrackCutSet == 2){
  	fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
	fTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);
	}
  if(fTrackCutSet == 3)fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb();
  if(fTrackCutSet == 4){
  	fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb();
	fTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);
	}

  fOutputList = new TList();
  fOutputList ->SetOwner();

  eff_MaxiPadLTM_All = new TEfficiency("eff_MaxiPadLTM_All"," ",72,0,72,23,0,23);
  fOutputList->Add(eff_MaxiPadLTM_All);
  eff_MaxiPadLTM_Mu = new TEfficiency("eff_MaxiPadLTM_Mu"," ",72,0,72,23,0,23);
  fOutputList->Add(eff_MaxiPadLTM_Mu);
  eff_MaxiPadLTM_El = new TEfficiency("eff_MaxiPadLTM_El"," ",72,0,72,23,0,23);
  fOutputList->Add(eff_MaxiPadLTM_El);

  eff_MaxiPadLTM_1Trk_All = new TEfficiency("eff_MaxiPadLTM_1Trk_All"," ",72,0,72,23,0,23);
  fOutputList->Add(eff_MaxiPadLTM_1Trk_All);
  eff_MaxiPadLTM_1Trk_Mu = new TEfficiency("eff_MaxiPadLTM_1Trk_Mu"," ",72,0,72,23,0,23);
  fOutputList->Add(eff_MaxiPadLTM_1Trk_Mu);
  eff_MaxiPadLTM_1Trk_El = new TEfficiency("eff_MaxiPadLTM_1Trk_El"," ",72,0,72,23,0,23);
  fOutputList->Add(eff_MaxiPadLTM_1Trk_El);
  
  eff_AverageTracklets = new TEfficiency("eff_AverageTracklets"," ",1000,0,10000);
  fOutputList->Add(eff_AverageTracklets);
  
  eff_AverageTrackPt = new TEfficiency("eff_AverageTrackPt"," ",50,0,5);
  fOutputList->Add(eff_AverageTrackPt);
  
  eff_MaxiPadLTM_Around = new TEfficiency("eff_MaxiPadLTM_Around"," ",72,0,72,23,0,23);
  fOutputList->Add(eff_MaxiPadLTM_Around);
  eff_MaxiPadLTM_OnlyAround = new TEfficiency("eff_MaxiPadLTM_OnlyAround"," ",72,0,72,23,0,23);
  fOutputList->Add(eff_MaxiPadLTM_OnlyAround);
  
  eff_MaxiPadLTM_Clusters = new TEfficiency("eff_MaxiPadLTM_Clusters"," ",72,0,72,23,0,23);
  fOutputList->Add(eff_MaxiPadLTM_Clusters);

  hTrackDistributionLTM = new TH2F("hTrackDistributionLTM","hTrackDistributionLTM",72,0,72,23,0,23);
  fOutputList->Add(hTrackDistributionLTM);
  hTrackDistribution_Mu = new TH2F("hTrackDistribution_Mu","hTrackDistribution_Mu",360,0,360,100,-1,1);
  fOutputList->Add(hTrackDistribution_Mu);
  hTrackDistribution_El = new TH2F("hTrackDistribution_El","hTrackDistribution_El",360,0,360,100,-1,1);
  fOutputList->Add(hTrackDistribution_El);
  hTrackDistribution = new TH2F("hTrackDistribution","hTrackDistribution",360,0,360,100,-1,1);
  fOutputList->Add(hTrackDistribution);
  hFiredMaxiPad = new TH2F("hFiredMaxiPad","hFiredMaxiPad",72,0,72,23,0,23);
  fOutputList->Add(hFiredMaxiPad);
  hFiredMaxiPadOnlyAround = new TH2F("hFiredMaxiPadOnlyAround","hFiredMaxiPadOnlyAround",72,0,72,23,0,23);
  fOutputList->Add(hFiredMaxiPadOnlyAround);
  
  hNotFiredMaxiPadCls = new TH2F("hNotFiredMaxiPadCls","hNotFiredMaxiPadCls",72,0,72,23,0,23);
  fOutputList->Add(hNotFiredMaxiPadCls);
  hExtraFiredMaxiPadCls = new TH2F("hExtraFiredMaxiPadCls","hExtraFiredMaxiPadCls",72,0,72,23,0,23);
  fOutputList->Add(hExtraFiredMaxiPadCls);
  hNotFiredMaxiPadTrk = new TH2F("hNotFiredMaxiPadTrk","hNotFiredMaxiPadTrk",72,0,72,23,0,23);
  fOutputList->Add(hNotFiredMaxiPadTrk);
  hExtraFiredMaxiPadTrk = new TH2F("hExtraFiredMaxiPadTrk","hExtraFiredMaxiPadTrk",72,0,72,23,0,23);
  fOutputList->Add(hExtraFiredMaxiPadTrk);
  
  hTrackPadCorrPhi = new TH2F("hTrackPadCorrPhi","hTrackPadCorrPhi",1440,0,360,72,0,72);
  fOutputList->Add(hTrackPadCorrPhi);
  hTrackPadCorrEta = new TH2F("hTrackPadCorrEta","hTrackPadCorrEta",1000,-1,1,23,0,23);
  fOutputList->Add(hTrackPadCorrEta);
  hNoiseMaxiPad = new TH2F("hNoiseMaxiPad","hNoiseMaxiPad",72,0,72,23,0,23);
  fOutputList->Add(hNoiseMaxiPad);
  hTriggerCounter = new TH1I("hTriggerCounter","hTriggerCounter",2,-0.5,1.5);
  hTriggerCounter->GetXaxis()->SetBinLabel(1,"CTRUE-B");
  hTriggerCounter->GetXaxis()->SetBinLabel(2,fTriggerClass.Data());
  fOutputList->Add(hTriggerCounter);
  hTriggerCounterIR1 = new TH1I("hTriggerCounterIR1","hTriggerCounterIR1",91,-0.5,90.5);
  fOutputList->Add(hTriggerCounterIR1);
  hTriggerCounterIR2 = new TH1I("hTriggerCounterIR2","hTriggerCounterIR2",91,-0.5,90.5);
  fOutputList->Add(hTriggerCounterIR2);
  hNFiredMaxiPads = new TH1F("hNFiredMaxiPads","hNFiredMaxiPads",1657,-0.5,1656.5);
  fOutputList->Add(hNFiredMaxiPads);
  hNFiredMaxiPadsOnlyAround = new TH1F("hNFiredMaxiPadsOnlyAround","hNFiredMaxiPadsOnlyAround",1657,-0.5,1656.5);
  fOutputList->Add(hNFiredMaxiPadsOnlyAround);
  hNTracklets = new TH1F("hNTracklets","hNTracklets",1001,-0.5,1000.5);
  fOutputList->Add(hNTracklets);
  hDetIn0 = new TH1I("hDetIn0","hDetIn0",18,-0.5,17.5),
  fOutputList->Add(hDetIn0);
  hDetIn1 = new TH1I("hDetIn1","hDetIn1",5,-0.5,4.5),
  fOutputList->Add(hDetIn1);
  hDetIn2 = new TH1I("hDetIn2","hDetIn2",19,-0.5,18.5),
  fOutputList->Add(hDetIn2);
  hDetIn3 = new TH1I("hDetIn3","hDetIn3",2,-0.5,1.5),
  fOutputList->Add(hDetIn3);
  hDetIn4 = new TH1I("hDetIn4","hDetIn4",48,-0.5,47.5),
  fOutputList->Add(hDetIn4);
  hPadDistance = new TH1F("hPadDistance","hPadDistance",1000,0,10);
  fOutputList->Add(hPadDistance);
  hTrackPt = new TH1F("hTrackPt","hTrackPt",500,0,50);
  fOutputList->Add(hTrackPt);
  hNMaxiPadIn = new TH1I("hNMaxiPadIn","hNMaxiPadIn",14,-3.5,10.5);
  fOutputList->Add(hNMaxiPadIn);
  hNCrossTracks = new TH1I("hNCrossTracks","hNCrossTracks",100,0.5,100.5);
  fOutputList->Add(hNCrossTracks);
  hBadMaxiPadMask = new TH2I("hBadMaxiPadMask","hBadMaxiPadMask",72,0,72,23,0,23);
  fOutputList->Add(hBadMaxiPadMask);
  hTOFHitTime = new TH1F("hTOFHitTime","hTOFHitTime",800000,-200000,600000);
  fOutputList->Add(hTOFHitTime);
  
  if(fUseEventSelection){
  	fEventCuts.AddQAplotsToList(fOutputList);
  	fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kAny);
	}

  PostData(1, fOutputList);

}//UserCreateOutputObjects


//_____________________________________________________________________________
void AliAnalysisTaskTOFTrigger::UserExec(Option_t *)
{


  TString fileName = ((TTree*) GetInputData(0))->GetCurrentFile()->GetName();
  
  Bool_t fBadMaxiPadMask[23][72];
  if(!fGeomLoaded){

	AliCDBManager *cdb = AliCDBManager::Instance();
	AliCDBEntry *cdbe = cdb->Get("TRIGGER/TOF/TriggerMask");
        AliTOFTriggerMask *fOCDBmask = (AliTOFTriggerMask *)cdbe->GetObject();
	
	//UInt_t BadLTMs[11] = {10,12,14,15,47,64,65,66,67,68,69};
	//UInt_t BadMaxiPads[5][2] = {{19,1}, {26,4}, {33,4}, {34,1}, {34,2}};
	UInt_t fgFromTriggertoDCS[72] = {0,1,4,5, 8, 9,12,13,16,17,20,21,24,25,28,29,32,33,36,37,40,41,44,45,48,49,52,53,56,57,60,61,64,65,68,69,
                                3,2,7,6,11,10,15,14,19,18,23,22,27,26,31,30,35,34,39,38,43,42,47,46,51,50,55,54,59,58,63,62,67,66,71,70};

	for(Int_t indexLTM=0; indexLTM<72; ++indexLTM) {
    		for(Int_t channelCTTM=0; channelCTTM<23; ++channelCTTM) {
			fBadMaxiPadMask[channelCTTM][indexLTM] = !fOCDBmask->IsON(fgFromTriggertoDCS[indexLTM],channelCTTM);
			//for(Int_t j = 0; j<11; j++)if(indexLTM == BadLTMs[j])fBadMaxiPadMask[channelCTTM][indexLTM] = 1;
			//for(Int_t j = 0; j<5; j++)if(indexLTM == BadMaxiPads[j][0] && channelCTTM == BadMaxiPads[j][1])fBadMaxiPadMask[channelCTTM][indexLTM] = 1;
			if(fBadMaxiPadMask[channelCTTM][indexLTM])hBadMaxiPadMask->SetBinContent(indexLTM+1,channelCTTM+1,1);
			else hBadMaxiPadMask->SetBinContent(indexLTM+1,channelCTTM+1,0);
			}
		}
	Int_t nAliveChannels = 0;
	for(Int_t channelCTTM=0; channelCTTM<23; ++channelCTTM){
		for(Int_t indexLTM=0; indexLTM<72; ++indexLTM){
			cout<<fBadMaxiPadMask[channelCTTM][indexLTM]<<",";
			nAliveChannels += !fBadMaxiPadMask[channelCTTM][indexLTM];
			}
		cout<<endl;
		}
	cout<<"N Active channels = "<<nAliveChannels<<endl;
	
  	for (int i=0;i<18;i++) {
   		AliGeomManager::GetOrigGlobalMatrix( Form("TOF/sm%02d",i) ,matOrig[i]);
   		matCurr[i] = *AliGeomManager::GetMatrix( Form("TOF/sm%02d",i) );
   		}
  	fGeomLoaded = kTRUE;
	}

  AliESDEvent *esd = (AliESDEvent*) InputEvent();
  if(!esd) return;

  const AliTOFHeader *tofH = esd->GetTOFHeader();
  fTOFmask = tofH->GetTriggerMask();

  TString trigger = esd->GetFiredTriggerClasses();
  
  Bool_t triggeredEvent = kFALSE;
  TString token;
  Ssiz_t from = 0;
  while (fTriggerClass.Tokenize(token, from, "[/]"))if(trigger.Contains(token.Data()))triggeredEvent=kTRUE;
  
  if(!triggeredEvent && !trigger.Contains("CTRUE-B")) return;

  TBits fIR1Map = esd->GetHeader()->GetIRInt1InteractionMap();
  TBits fIR2Map = esd->GetHeader()->GetIRInt2InteractionMap();
  Int_t fClosestIR1 = 100;
  Int_t fClosestIR2 = 100;
  for(Int_t item=-1; item>=-90; item--) {
    Int_t bin = 90+item;
    Bool_t isFired = fIR1Map.TestBitNumber(bin);
    if(isFired) {
      fClosestIR1 = TMath::Abs(item);
      break;
    }
  if(fClosestIR1 == 100)fClosestIR1 = 0;
  }
  for(Int_t item=-1; item>=-90; item--) {
    Int_t bin = 90+item;
    Bool_t isFired = fIR2Map.TestBitNumber(bin);
    if(isFired) {
      fClosestIR2 = TMath::Abs(item);
      break;
    }
  }
  if(fClosestIR2 == 100)fClosestIR2 = 0;

  if(triggeredEvent){
  	hTriggerCounterIR1->Fill(fClosestIR1);
	hTriggerCounterIR2->Fill(fClosestIR2);
	}
  if(fClosestIR1 < fMaxBCs && fClosestIR1 != 0)return;
  if(fClosestIR2 < fMaxBCs && fClosestIR2 != 0)return;

  Bool_t isGoodCTRUE = trigger.Contains("CTRUE-B") && !esd->GetHeader()->IsTriggerInputFired("VBA") && !esd->GetHeader()->IsTriggerInputFired("VBC") && !esd->GetHeader()->IsTriggerInputFired("SMB");

  if(isGoodCTRUE)hTriggerCounter->Fill(0);
  if(triggeredEvent)hTriggerCounter->Fill(1);

  for (Int_t ltm=0;ltm<72;ltm++){
   	for (Int_t channel=0;channel<23;channel++){
		if(fTOFmask->IsON(ltm,channel) && isGoodCTRUE)hNoiseMaxiPad->Fill(ltm,channel);
		if(fTOFmask->IsON(ltm,channel) && triggeredEvent)hFiredMaxiPad->Fill(ltm,channel);
		}
	}

  if(!triggeredEvent) return;

  hNFiredMaxiPads->Fill(fTOFmask->GetNumberMaxiPadOn());
  Int_t fNtracklets = esd->GetMultiplicity()->GetNumberOfTracklets();
  hNTracklets->Fill(fNtracklets);
  if(fNtracklets>fMaxMulti) return;
  
  if(fUseEventSelection){
  	if(!fEventCuts.AcceptEvent(esd)){
  		PostData(1, fOutputList);
  		return;
		}
	}
	
  Int_t numTracksPerMaxiPad[72][23];
  Int_t numMuonTracksPerMaxiPad[72][23];
  Int_t numElectronTracksPerMaxiPad[72][23];
  for (Int_t indexLTM=0; indexLTM<72; ++indexLTM) {
    for (Int_t channelCTTM=0; channelCTTM<23; ++channelCTTM) {
      numTracksPerMaxiPad[indexLTM][channelCTTM] = 0;
      numMuonTracksPerMaxiPad[indexLTM][channelCTTM] = 0;
      numElectronTracksPerMaxiPad[indexLTM][channelCTTM] = 0;
    }
  }

  TClonesArray* tofClusters = ((AliESDEvent*)fInputEvent)->GetESDTOFClusters();
  
  Int_t nTOFhits = 0;
  for (Int_t icl=0;icl<tofClusters->GetEntriesFast();icl++){
     AliESDTOFCluster* cl = (AliESDTOFCluster*) tofClusters->At(icl);
     nTOFhits+=cl->GetNTOFhits();
   }
   
   TArrayI fTOFhits;
   TArrayI fTrackIndices;
   
   fTOFhits.Reset();
   fTrackIndices.Reset();
   fTOFhits.Set(nTOFhits);
   fTrackIndices.Set(nTOFhits);
  
  Int_t hitCounts = 0;
   for (Int_t icl=0;icl<tofClusters->GetEntriesFast();icl++){
     AliESDTOFCluster* cl = (AliESDTOFCluster*) tofClusters->At(icl);
     for (Int_t ihit=0;ihit<cl->GetNTOFhits();ihit++){
       AliESDTOFHit* hit = (AliESDTOFHit*) cl->GetTOFHit(ihit);
       hTOFHitTime->Fill(hit->GetTime());
       if(hit->GetTime() < fMinTOF  || hit->GetTime() > fMaxTOF)continue;
       Int_t channel = hit->GetTOFchannel();
       Int_t trackIndex = (cl->GetNMatchableTracks()==1) ? cl->GetTrackIndex(0) : -1;
       fTOFhits.AddAt(channel,hitCounts);
       fTrackIndices.AddAt(trackIndex,hitCounts);
       hitCounts++;
     }
   }  
  
  Double_t cv[21];
   
  //Track loop
  for(Int_t iTrack=0; iTrack<esd->GetNumberOfTracks(); iTrack++) {
    AliESDtrack *esdTrackOrig = dynamic_cast<AliESDtrack*>(esd->GetTrack(iTrack));
    if( !esdTrackOrig ) continue;

    if(!fTrackCuts->AcceptTrack(esdTrackOrig))continue;
    hTrackPt->Fill(esdTrackOrig->Pt());
    if(esdTrackOrig->Pt()>fMaxPt || esdTrackOrig->Pt()<fMinPt)continue;
    
    Int_t detId[5] = {-1,-1,-1,-1,-1};
    Int_t indexLTM[2] = {-1,-1};
    for(Int_t i =0; i<fTrackIndices.GetSize(); i++){
    	if(iTrack == fTrackIndices.GetAt(i)){
		AliTOFGeometry::GetVolumeIndices(fTOFhits.GetAt(i),detId);
		GetLTMIndex(detId,indexLTM);
		UInt_t channelCTTM = indexLTM[1]/2;
    		eff_MaxiPadLTM_Clusters->Fill(fTOFmask->IsON(indexLTM[0],channelCTTM),indexLTM[0],channelCTTM);
		//cout<<"Track "<<iTrack<<" Trigger pad in LTM "<<indexLTM[0]<<" CTTM "<<channelCTTM<<" Fired = "<<fTOFmask->IsON(indexLTM[0],channelCTTM)<<endl;
		if(trigger.Contains("CCUP30-B") || trigger.Contains("CCUP31-B")){
			if(!fTOFmask->IsON(indexLTM[0],channelCTTM) && (fTOFmask->GetNumberMaxiPadOn()< 2))hNotFiredMaxiPadCls->Fill(indexLTM[0],channelCTTM);
			}
		if(trigger.Contains("CCUP31-B")){
			if(fTOFmask->IsON(indexLTM[0],channelCTTM) && (fTOFmask->GetNumberMaxiPadOn()> 6))hExtraFiredMaxiPadCls->Fill(indexLTM[0],channelCTTM);
			}
		}
	}
    
    AliESDtrack *esdTrack = (AliESDtrack*)esdTrackOrig->Clone();

    AliExternalTrackParam* trc = (AliExternalTrackParam*)esdTrack->GetOuterParam();
    if (!trc){
    	cout<<"No outer param !!!!"<<endl;
    	continue; // no outer param
    }
	
    if (!AliTrackerBase::PropagateTrackToBxByBz(trc, AliTOFGeometry::RinTOF(), esdTrack->GetMassForTracking(), 1.0, kTRUE)) continue; //propagation failed
    
    // go to specific frame sector and reach target X in this frame:
    Double_t phi = trc->PhiPos()*TMath::RadToDeg();
    if (phi<0) phi += 360;
    Int_t sectOld = -1, sect = int(phi/20.);
    Bool_t failed = kFALSE;
    while(sectOld!=sect) {
    	sectOld = sect;
   	Double_t alpha = (sect*20.+10)*TMath::DegToRad();
	if (!trc->Rotate(alpha) || !AliTrackerBase::PropagateTrackToBxByBz(trc, AliTOFGeometry::RinTOF(), esdTrack->GetMassForTracking(), 1.0, kFALSE)){ 
		// don't rotate at every step anymore
        	failed = kTRUE; 
		break;
   		}
     	// make sure the propagation did not change the sector
     	phi = trc->PhiPos()*TMath::RadToDeg();
     	if (phi<0) phi += 360;
     	sect = int(phi/20.0);
    	}
    if (failed) continue;
    
    //cout<<"Track"<<endl;
    //cout<<"At RinTOF uncertainty = "<<cv[0]<<" ; "<<cv[2]<<" ; "<<cv[5]<<endl;
    
    trc->GetCovarianceXYZPxPyPz(cv);
    if (cv[0]<0 || TMath::Sqrt(cv[0])>fMaxTrackError){hNMaxiPadIn->Fill(-2); continue;}
    if (cv[2]<0 || TMath::Sqrt(cv[2])>fMaxTrackError){hNMaxiPadIn->Fill(-2); continue;}
    if (cv[5]<0 || TMath::Sqrt(cv[5])>fMaxTrackError){hNMaxiPadIn->Fill(-2); continue;}
    
    if(esdTrackOrig->GetTOFsignal()>99998){hNMaxiPadIn->Fill(-3); continue;}

    //Fine propagation from TOF radius
    Bool_t isin = kFALSE;
    Float_t dist3d[3]={-1.,-1.,-1.}; // residual to TOF channel
    Float_t rTOFused = AliTOFGeometry::RinTOF();
    Double_t pos[3]={0.0,0.0,0.0};
    Float_t posF_In[3]={0.0,0.0,0.0};
    Float_t posF_Out[3]={0.0,0.0,0.0};
    Double_t locTmp[3]={0.0,0.0,0.0};
    UInt_t nFiredPads = 0;
    UInt_t instep = 0;
    
    while (rTOFused<AliTOFGeometry::Rmax()){
    	rTOFused += 0.1; // 1 mm step
	instep++;
	
    	if(!AliTrackerBase::PropagateTrackParamOnlyToBxByBz(trc,rTOFused,1,kFALSE)){hNMaxiPadIn->Fill(-1); break;}
	phi = trc->PhiPos()*TMath::RadToDeg();
        if (phi<0) phi += 360;
        sect = Int_t(phi/20.0);
        if (sect!=sectOld){
        	sectOld = sect;
        	Double_t alpha = (sect*20.0+10)*TMath::DegToRad();
        	if (!trc->Rotate(alpha)){hNMaxiPadIn->Fill(-1); break;}
		}
	
    	trc->GetXYZ(pos);
	matCurr[sect].MasterToLocal(pos, locTmp); // go to sector local frame, accounting for misaligment
	matOrig[sect].LocalToMaster(locTmp, pos); // go back to ideal lab frame
	
	posF_In[0] = pos[0];
	posF_In[1] = pos[1];
	posF_In[2] = pos[2];
	dist3d[0] = AliTOFGeometry::GetPadDx(posF_In);
     	dist3d[1] = AliTOFGeometry::GetPadDy(posF_In);
     	dist3d[2] = AliTOFGeometry::GetPadDz(posF_In);
	
	if(TMath::Abs(dist3d[0]) < 1.25 && TMath::Abs(dist3d[1]) < 0.1 && TMath::Abs(dist3d[2]) < 1.75){
		isin= kTRUE;
		nFiredPads++;
		hPadDistance->Fill(TMath::Sqrt(dist3d[0]*dist3d[0]+dist3d[1]*dist3d[1]+dist3d[2]*dist3d[2]));

		//cout<<"Is in, radius = "<<rTOFused<<" Step = "<<instep<<endl;
		UInt_t outstep = 0;
		while (rTOFused<AliTOFGeometry::Rmax()){
    			rTOFused += 0.1; // 1 mm step
			outstep++;
			
    			if(!AliTrackerBase::PropagateTrackParamOnlyToBxByBz(trc,rTOFused,1,kFALSE)){hNMaxiPadIn->Fill(-1); break;}
			phi = trc->PhiPos()*TMath::RadToDeg();
        		if (phi<0) phi += 360;
       			sect = Int_t(phi/20.0);
        		if (sect!=sectOld){
        			sectOld = sect;
        			Double_t alpha = (sect*20.0+10)*TMath::DegToRad();
        			if (!trc->Rotate(alpha)){hNMaxiPadIn->Fill(-1); break;}
				}
	
    			trc->GetXYZ(pos);
			matCurr[sect].MasterToLocal(pos, locTmp); // go to sector local frame, accounting for misaligment
			matOrig[sect].LocalToMaster(locTmp, pos); // go back to ideal lab frame
			
			posF_Out[0] = pos[0];
			posF_Out[1] = pos[1];
			posF_Out[2] = pos[2];
			dist3d[0] = AliTOFGeometry::GetPadDx(posF_Out);
     			dist3d[1] = AliTOFGeometry::GetPadDy(posF_Out);
     			dist3d[2] = AliTOFGeometry::GetPadDz(posF_Out);
			if(dist3d[0] == -2 || dist3d[1] == -2 || dist3d[2] == -2){
				//cout<<"Is out, radius = "<<rTOFused<<" Step = "<<outstep<<endl;
				break;
				}
			}
		}
	if(isin){
		isin = kFALSE;
		detId[0]=AliTOFGeometry::GetSector(posF_In);
		detId[1]=AliTOFGeometry::GetPlate(posF_In);
		detId[2]=AliTOFGeometry::GetStrip(posF_In);
		detId[3]=AliTOFGeometry::GetPadZ(posF_In);
		detId[4]=AliTOFGeometry::GetPadX(posF_In);
		//cout<<"Crossed pad "<<nFiredPads<<" in ID "<<detId[0]<<" "<<detId[1]<<" "<<detId[2]<<" "<<detId[3]<<" "<<detId[4]<<endl;
		hDetIn0->Fill(detId[0]);
		hDetIn1->Fill(detId[1]);
		hDetIn2->Fill(detId[2]);
		hDetIn3->Fill(detId[3]);
		hDetIn4->Fill(detId[4]);

		GetLTMIndex(detId,indexLTM);
		UInt_t channelCTTM = indexLTM[1]/2;
		  
		//cout<<"Trigger pad "<<nFiredPads<<" in LTM "<<indexLTM[0]<<" CTTM "<<channelCTTM<<" Fired = "<<fTOFmask->IsON(indexLTM[0],channelCTTM)<<endl;

		hTrackDistribution->Fill(trc->Phi()*TMath::RadToDeg(),trc->Eta());
		hTrackDistributionLTM->Fill(indexLTM[0],channelCTTM);

		if(fTOFmask->IsON(indexLTM[0],channelCTTM)){
		    hTrackPadCorrPhi->Fill(trc->Phi()*TMath::RadToDeg(),indexLTM[0]);
		    hTrackPadCorrEta->Fill(trc->Eta(),channelCTTM);
		    }

		if(!fBadMaxiPadMask[channelCTTM][indexLTM[0]])eff_AverageTrackPt->Fill(fTOFmask->IsON(indexLTM[0],channelCTTM),esdTrack->Pt());
                numTracksPerMaxiPad[indexLTM[0]][channelCTTM] += 1;

		Float_t fPIDTPCMuon = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kMuon);
		Float_t fPIDTPCElectron = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kElectron);
		if(TMath::Abs(fPIDTPCMuon) < TMath::Abs(fPIDTPCElectron)){
		    hTrackDistribution_Mu->Fill(trc->Phi()*TMath::RadToDeg(),trc->Eta());
                    numMuonTracksPerMaxiPad[indexLTM[0]][channelCTTM] += 1;
		    }
		if(TMath::Abs(fPIDTPCMuon) > TMath::Abs(fPIDTPCElectron)){
		    hTrackDistribution_El->Fill(trc->Phi()*TMath::RadToDeg(),trc->Eta());
                    numElectronTracksPerMaxiPad[indexLTM[0]][channelCTTM] += 1;
		    }
		    
		if(trigger.Contains("CCUP30-B") || trigger.Contains("CCUP31-B")){
			if(!fTOFmask->IsON(indexLTM[0],channelCTTM) && (fTOFmask->GetNumberMaxiPadOn()< 2))hNotFiredMaxiPadTrk->Fill(indexLTM[0],channelCTTM);
			}
		if(trigger.Contains("CCUP31-B")){
			if(fTOFmask->IsON(indexLTM[0],channelCTTM) && (fTOFmask->GetNumberMaxiPadOn()> 6))hExtraFiredMaxiPadTrk->Fill(indexLTM[0],channelCTTM);
			}
		}
     	}
	hNMaxiPadIn->Fill(nFiredPads);
  }

  // filling TEfficiency object
  for (Int_t indexLTM=0; indexLTM<72; ++indexLTM) {
    for (Int_t channelCTTM=0; channelCTTM<23; ++channelCTTM) {
      hNCrossTracks->Fill(numTracksPerMaxiPad[indexLTM][channelCTTM]);
      const Bool_t isON = fTOFmask->IsON(indexLTM, channelCTTM);
      
      Int_t phiPlus = indexLTM + 1;
      if(phiPlus == 36) phiPlus = 0;
      if(phiPlus == 72) phiPlus = 36;
      
      Int_t phiMinus = indexLTM - 1;
      if(phiMinus == -1) phiMinus = 35;
      if(phiMinus == 35) phiMinus = 71;
      
      Int_t etaPlus = channelCTTM + 1;
      if(etaPlus > 22){
      	etaPlus = channelCTTM -1;
	if(phiPlus<36)phiPlus += 35;
	if(phiPlus>35)phiPlus -= 35;
	if(phiMinus<36)phiMinus += 35;
	if(phiMinus>35)phiMinus -= 35;
	
	}
      Int_t etaMinus = channelCTTM - 1;
      if(etaMinus < 0)etaMinus = channelCTTM + 1;
      
      const Bool_t isONaround = (fTOFmask->IsON(indexLTM, etaPlus) || fTOFmask->IsON(indexLTM, etaMinus)
      			 ||fTOFmask->IsON(phiPlus, channelCTTM) || fTOFmask->IsON(phiMinus, channelCTTM)
			 ||fTOFmask->IsON(phiPlus, etaPlus) || fTOFmask->IsON(phiPlus, etaMinus)
			 ||fTOFmask->IsON(phiMinus, etaPlus) || fTOFmask->IsON(phiMinus, etaMinus));
      
      for (Int_t l=0; l<numTracksPerMaxiPad[indexLTM][channelCTTM]; ++l){
        eff_MaxiPadLTM_All->Fill(isON, indexLTM, channelCTTM);
	eff_MaxiPadLTM_Around->Fill(isON||isONaround, indexLTM, channelCTTM);
	eff_MaxiPadLTM_OnlyAround->Fill(!isON||isONaround, indexLTM, channelCTTM);
	
	if(isONaround && !isON)hFiredMaxiPadOnlyAround->Fill(indexLTM, channelCTTM);
	if(isONaround && !isON)hNFiredMaxiPadsOnlyAround->Fill(fTOFmask->GetNumberMaxiPadOn());
	if(!fBadMaxiPadMask[channelCTTM][indexLTM])eff_AverageTracklets->Fill(isON,fNtracklets);
	}

      for (Int_t l=0; l<numMuonTracksPerMaxiPad[indexLTM][channelCTTM]; ++l)
        eff_MaxiPadLTM_Mu->Fill(isON, indexLTM, channelCTTM);

      for (Int_t l=0; l<numElectronTracksPerMaxiPad[indexLTM][channelCTTM]; ++l)
        eff_MaxiPadLTM_El->Fill(isON, indexLTM, channelCTTM);

      if (numTracksPerMaxiPad[indexLTM][channelCTTM] == 1) {
        eff_MaxiPadLTM_1Trk_All->Fill(isON, indexLTM, channelCTTM);

        if (numMuonTracksPerMaxiPad[indexLTM][channelCTTM])
          eff_MaxiPadLTM_1Trk_Mu->Fill(isON, indexLTM, channelCTTM);

        if (numElectronTracksPerMaxiPad[indexLTM][channelCTTM])
          eff_MaxiPadLTM_1Trk_El->Fill(isON, indexLTM, channelCTTM);
      }
    }
  }

  PostData(1, fOutputList);

}//UserExec


//_____________________________________________________________________________
void AliAnalysisTaskTOFTrigger::Terminate(Option_t *)
{
  cout<<"Analysis complete."<<endl;

}//Terminate


//_____________________________________________________________________________
void AliAnalysisTaskTOFTrigger::GetLTMIndex(const Int_t * const detind, Int_t *indexLTM) {
  //
  // getting LTMmatrix indexes for current digit
  //

  Int_t iStrip = 19*detind[1]+detind[2] + (detind[1]>2 ? -4: 0);
  indexLTM[0] = detind[0]*2;
  if (detind[4]<24) indexLTM[0] += iStrip<45 ? 0 : 36;
  else              indexLTM[0] += iStrip<46 ? 1 : 37;
  if (indexLTM[0]<36) indexLTM[1] = iStrip;
  else                indexLTM[1] = 90-iStrip;
  
}
