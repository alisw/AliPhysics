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
	hTrackDistributionLTM(0),
	hTrackDistribution_Mu(0),
	hTrackDistribution_El(0),
	hTrackDistribution(0),
	hFiredMaxiPad(0),
	hNotFiredMaxiPad(0),
	hTrackPadCorrPhi(0),  
	hTrackPadCorrEta(0),
	hNoiseMaxiPad(0),
	hTriggerCounter(0),
	hTriggerCounterIR1(0),
	hTriggerCounterIR2(0),
	hNFiredMaxiPads(0),
	hDetIn0(0),
	hDetIn1(0),
	hDetIn2(0),
	hDetIn3(0),
	hDetIn4(0),
	hPadDistance(0),
	hTrackPt(0),
	hNMaxiPadIn(0),
	fGeomLoaded(kFALSE),
	fMaxPt(0),
	fMinPt(0),
	fMaxMulti(0),
	fTriggerClass(0),
	fMaxBCs(0)
	   

{

//Dummy constructor

}//AliAnalysisTaskTOFTrigger


//_____________________________________________________________________________
AliAnalysisTaskTOFTrigger::AliAnalysisTaskTOFTrigger(const char *name,Float_t lowpt,Float_t highpt,Int_t highmult,TString trgcls,Int_t nBCs) 
  : AliAnalysisTaskSE(name),fOutputList(0),fPIDResponse(0),fTrackCuts(0),
	fTOFmask(0),
	eff_MaxiPadLTM_All(0),
	eff_MaxiPadLTM_Mu(0),
	eff_MaxiPadLTM_El(0),
	hTrackDistributionLTM(0),
	hTrackDistribution_Mu(0),
	hTrackDistribution_El(0),
	hTrackDistribution(0),
	hFiredMaxiPad(0),
	hNotFiredMaxiPad(0),
	hTrackPadCorrPhi(0),  
	hTrackPadCorrEta(0),
	hNoiseMaxiPad(0),
	hTriggerCounter(0),
	hTriggerCounterIR1(0),
	hTriggerCounterIR2(0),
	hNFiredMaxiPads(0),
	hDetIn0(0),
	hDetIn1(0),
	hDetIn2(0),
	hDetIn3(0),
	hDetIn4(0),
	hPadDistance(0),
	hTrackPt(0),
	hNMaxiPadIn(0),
	fGeomLoaded(kFALSE),
	fMaxPt(highpt),
	fMinPt(lowpt),
	fMaxMulti(highmult),
	fTriggerClass(trgcls),
	fMaxBCs(nBCs)    

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
  
  fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();

  fOutputList = new TList();
  fOutputList ->SetOwner();
  
  eff_MaxiPadLTM_All = new TEfficiency("eff_MaxiPadLTM_All"," ",72,0,72,23,0,23);
  fOutputList->Add(eff_MaxiPadLTM_All);
  eff_MaxiPadLTM_Mu = new TEfficiency("eff_MaxiPadLTM_Mu"," ",72,0,72,23,0,23);
  fOutputList->Add(eff_MaxiPadLTM_Mu);
  eff_MaxiPadLTM_El = new TEfficiency("eff_MaxiPadLTM_El"," ",72,0,72,23,0,23);
  fOutputList->Add(eff_MaxiPadLTM_El);
  
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
  hNotFiredMaxiPad = new TH2F("hNotFiredMaxiPad","hNotFiredMaxiPad",72,0,72,23,0,23);
  fOutputList->Add(hNotFiredMaxiPad);
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
  hNMaxiPadIn = new TH1I("hNMaxiPadIn","hNMaxiPadIn",13,-2.5,10.5);
  fOutputList->Add(hNMaxiPadIn);
  
  PostData(1, fOutputList);

}//UserCreateOutputObjects


//_____________________________________________________________________________
void AliAnalysisTaskTOFTrigger::UserExec(Option_t *) 
{

  if(!fGeomLoaded){
  	AliGeomManager::LoadGeometry();
  	AliGeomManager::ApplyAlignObjsFromCDB("ITS TRD TOF");
  	fGeomLoaded = kTRUE;
	}

  AliESDEvent *esd = (AliESDEvent*) InputEvent();
  if(!esd) return;
  
  const AliTOFHeader *tofH = esd->GetTOFHeader();
  fTOFmask = tofH->GetTriggerMask();
  
  TString trigger = esd->GetFiredTriggerClasses();
  if(!trigger.Contains(fTriggerClass.Data()) && !trigger.Contains("CTRUE-B")) return;
  
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
  
  if(trigger.Contains(fTriggerClass.Data())){
  	hTriggerCounterIR1->Fill(fClosestIR1);
	hTriggerCounterIR2->Fill(fClosestIR2);
	}
  if(fClosestIR1 < fMaxBCs && fClosestIR1 != 0)return;
  if(fClosestIR2 < fMaxBCs && fClosestIR2 != 0)return;

  Bool_t isGoodCTRUE = trigger.Contains("CTRUE-B") && !esd->GetHeader()->IsTriggerInputFired("VBA") && !esd->GetHeader()->IsTriggerInputFired("VBC") && !esd->GetHeader()->IsTriggerInputFired("SH2");

  if(isGoodCTRUE)hTriggerCounter->Fill(0); 
  if(trigger.Contains(fTriggerClass.Data()))hTriggerCounter->Fill(1);
  
  for (Int_t ltm=0;ltm<72;ltm++){
   	for (Int_t channel=0;channel<23;channel++){
		if(fTOFmask->IsON(ltm,channel) && isGoodCTRUE)hNoiseMaxiPad->Fill(ltm,channel);
		if(fTOFmask->IsON(ltm,channel) && trigger.Contains(fTriggerClass.Data()))hFiredMaxiPad->Fill(ltm,channel);
		}
	}

  if(!trigger.Contains(fTriggerClass.Data())) return;

  hNFiredMaxiPads->Fill(fTOFmask->GetNumberMaxiPadOn());
  if(fTOFmask->GetNumberMaxiPadOn()>fMaxMulti) return;

  //Track loop
  for(Int_t iTrack=0; iTrack<esd->GetNumberOfTracks(); iTrack++) {
    AliESDtrack *esdTrack = dynamic_cast<AliESDtrack*>(esd->GetTrack(iTrack));
    if( !esdTrack ) continue;
    
    if(!fTrackCuts->AcceptTrack(esdTrack))continue;
    hTrackPt->Fill(esdTrack->Pt());
    if(esdTrack->Pt()>fMaxPt || esdTrack->Pt()<fMinPt)continue;
    
    AliExternalTrackParam* trc = (AliExternalTrackParam*)esdTrack->GetOuterParam();
    if (!trc){
    	cout<<"No outer param !!!!"<<endl; 
    	continue; // no outer param
    }

    Double_t phi = trc->Phi()*TMath::RadToDeg();
    if (phi<0) phi += 360;
    Int_t sect = Int_t(phi/20.);
    Double_t alpha = (sect*20.+10)*TMath::DegToRad();
    if (!trc->Rotate(alpha)) continue; // failed

    if (!AliTrackerBase::PropagateTrackToBxByBz(trc, AliTOFGeometry::RinTOF(), esdTrack->GetMassForTracking(), 1.0, kFALSE)) continue; // propagation failed
    
    phi = trc->Phi()*TMath::RadToDeg();
    if (phi<0) phi += 360;
    Int_t sect1 = int(phi/20.);
    if (sect!=sect1) {
    alpha = (sect1*20.+10)*TMath::DegToRad();
    if (!trc->Rotate(alpha)) continue; // failed
    if (!AliTrackerBase::PropagateTrackToBxByBz(trc, AliTOFGeometry::RinTOF(),esdTrack->GetMassForTracking(), 1.0, kFALSE)) continue; // propagation failed
    sect = sect1;
    }
    
    //Fine propagation from TOF radius
    Bool_t isin = kFALSE;
    Float_t dist3d[3]={-1.,-1.,-1.}; // residual to TOF channel
    Float_t rTOFused = AliTOFGeometry::RinTOF(); 
    UInt_t nmaxstep = 500; // to be tuned
    Double_t pos[3]={0.0,0.0,0.0};
    Float_t posF_In[3]={0.0,0.0,0.0};
    Float_t posF_Out[3]={0.0,0.0,0.0};
    UInt_t nFiredPads = 0;
    //cout<<"Track"<<endl;
    for(UInt_t instep = 0; instep < nmaxstep; instep++){
    	rTOFused += 0.1; // 1 mm step     
    	if(!trc->PropagateTo(rTOFused,esd->GetMagneticField())){hNMaxiPadIn->Fill(-1); break;}
    	trc->GetXYZ(pos);
	posF_In[0] = pos[0];
	posF_In[1] = pos[1];
	posF_In[2] = pos[2];
	dist3d[0] = AliTOFGeometry::GetPadDx(posF_In);
     	dist3d[1] = AliTOFGeometry::GetPadDy(posF_In);
     	dist3d[2] = AliTOFGeometry::GetPadDz(posF_In);
	if(dist3d[0] != -2 && dist3d[1] != -2 && dist3d[2] != -2){
		isin= kTRUE;
		nFiredPads++;
		hPadDistance->Fill(TMath::Sqrt(dist3d[0]*dist3d[0]+dist3d[1]*dist3d[1]+dist3d[2]*dist3d[2]));
		
		//cout<<"Is in, radius = "<<rTOFused<<" Step = "<<instep<<endl;
		
		for(UInt_t outstep = 0; outstep < nmaxstep; outstep++){
    			rTOFused += 0.1; // 1 mm step     
    			if(!trc->PropagateTo(rTOFused,esd->GetMagneticField())){hNMaxiPadIn->Fill(-1); break;}
    			trc->GetXYZ(pos);
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
		Int_t detId[5] = {-1,-1,-1,-1,-1};
		Int_t indexLTM[2] = {-1,-1};
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
		    	 
		eff_MaxiPadLTM_All->Fill(fTOFmask->IsON(indexLTM[0],channelCTTM),indexLTM[0],channelCTTM);
		
		Float_t fPIDTPCMuon = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kMuon);
		Float_t fPIDTPCElectron = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kElectron);
		if(TMath::Abs(fPIDTPCMuon) < TMath::Abs(fPIDTPCElectron)){
		    hTrackDistribution_Mu->Fill(trc->Phi()*TMath::RadToDeg(),trc->Eta());
		    eff_MaxiPadLTM_Mu->Fill(fTOFmask->IsON(indexLTM[0],channelCTTM),indexLTM[0],channelCTTM);
		    }
		if(TMath::Abs(fPIDTPCMuon) > TMath::Abs(fPIDTPCElectron)){
		    hTrackDistribution_El->Fill(trc->Phi()*TMath::RadToDeg(),trc->Eta());
		    eff_MaxiPadLTM_El->Fill(fTOFmask->IsON(indexLTM[0],channelCTTM),indexLTM[0],channelCTTM);
		    }
		}
     	}
	hNMaxiPadIn->Fill(nFiredPads);
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
