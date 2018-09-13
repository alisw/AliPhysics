#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliITSRecPoint.h"
#include "AliESDEvent.h"
#include "AliTrackPointArray.h"
#include "AliITSgeomTGeo.h"
#include "AliESDfriend.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSresponseSDD.h"
#include "AliTriggerConfiguration.h"
#include "AliGeomManager.h"
#include <TSystem.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TChain.h>
#include <TGeoGlobalMagField.h>
#include "AliESDInputHandlerRP.h"
/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
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

//*************************************************************************
// Implementation of class AliAnalysiTaskSDDRP
// AliAnalysisTaskSE to extract from ESD + ESDfreinds + ITS rec points
// performance plots for SDD detector
//
// Author: F. Prino, prino@to.infn.it
//*************************************************************************


#include "AliAnalysisTaskSDDRP.h"

ClassImp(AliAnalysisTaskSDDRP)
//______________________________________________________________________________
AliAnalysisTaskSDDRP::AliAnalysisTaskSDDRP() : AliAnalysisTaskSE("SDD RecPoints"), 
  fOutput(0),
  fHistNEvents(0),
  fHistCluInLay(0),
  fHistAllPMod(0),
  fHistGoodPMod(0),
  fHistBadRegMod(0),
  fHistMissPMod(0),
  fHistSkippedMod(0),
  fHistOutAccMod(0),
  fHistNoRefitMod(0),
  fHistAllPXloc(0),
  fHistGoodPXloc(0),
  fHistBadRegXloc(0),
  fHistMissPXloc(0),
  fHistAllPZloc(0),
  fHistGoodPZloc(0),
  fHistBadRegZloc(0),
  fHistMissPZloc(0),
  fHistdEdxL3VsP(0),
  fHistdEdxL4VsP(0),
  fHistdEdxVsMod(0),
  fRecPMod(0),
  fTrackPMod(0),
  fGoodAnMod(0),
  fRecPLadLay3(0),
  fRecPLadLay4(0),
  fTrackPLadLay3(0),
  fTrackPLadLay4(0),
  fGoodAnLadLay3(0),
  fGoodAnLadLay4(0),
  fEtaPhiTracks(0),
  fEtaPhiTracksLay3(0),
  fEtaPhiTracksLay4(0),
  fDriftTimeRP(0),
  fDriftTimeTPAll(0),
  fDriftTimeTPAllMod(0),
  fDriftTimeTPNoExtra(0),
  fDriftTimeTPExtra(0),
  fCluSizAnVsTime(0),
  fCluSizTbVsTime(0),
  fProfRecPtsLay3VsTime(0),
  fProfRecPtsLay4VsTime(0),
  fProfTrPtsLay3VsTime(0),
  fProfTrPtsLay4VsTime(0),
  fProfFracTrRecLay3VsTime(0),
  fProfFracTrRecLay4VsTime(0),
  fProfFracTrkWithPntLay3VsTime(0),
  fProfFracTrkWithPntLay4VsTime(0),
  fResp(0),
  fTrigConfig(0),
  fUseITSsaTracks(kFALSE),
  fMinITSpts(3),
  fMinTPCpts(70),
  fMinPfordEdx(0.5),
  fTriggerClass(""),
  fOnlyEventsWithSDD(kTRUE),
  fExcludeBadMod(kFALSE),
  fReadCDB(kTRUE),
  fInitCalib(kFALSE)
{
  //
  DefineOutput(1, TList::Class());
}


//___________________________________________________________________________
AliAnalysisTaskSDDRP::~AliAnalysisTaskSDDRP(){
  //
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput;
    fOutput = 0;
  }
  
}


//___________________________________________________________________________

void AliAnalysisTaskSDDRP::UserCreateOutputObjects() {

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  fHistNEvents = new TH1F("hNEvents", "Number of processed events",5,-1.5,3.5);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"All events");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"Selected triggers"); 
  fHistNEvents->GetXaxis()->SetBinLabel(3,"Without SDD"); 
  fHistNEvents->GetXaxis()->SetBinLabel(4,"With SDD"); 
  fHistNEvents->GetXaxis()->SetBinLabel(5,"Analyzed events"); 
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);

  fHistCluInLay = new TH1F("hCluInLay","hCluInLay",7,-1.5,5.5);
  fHistCluInLay->Sumw2();
  fHistCluInLay->SetMinimum(0);
  fOutput->Add(fHistCluInLay);

  // -- Module histos

  fHistAllPMod  = new TH1F("hAllPmod","Crossing Tracks vs. Module",260,239.5,499.5);
  fHistAllPMod->Sumw2();
  fHistAllPMod->SetMinimum(0);
  fOutput->Add(fHistAllPMod);

  fHistGoodPMod  = new TH1F("hGoodPmod","PointsAssocToTrack per Module",260,239.5,499.5);
  fHistGoodPMod->Sumw2();
  fHistGoodPMod->SetMinimum(0);
  fOutput->Add(fHistGoodPMod);

  fHistBadRegMod  = new TH1F("hBadRegmod","Tracks in BadRegion per Module",260,239.5,499.5);
  fHistBadRegMod->Sumw2();
  fHistBadRegMod->SetMinimum(0);
  fOutput->Add(fHistBadRegMod);

  fHistMissPMod  = new TH1F("hMissPmod","Missing Points per Module",260,239.5,499.5);
  fHistMissPMod->Sumw2();
  fHistMissPMod->SetMinimum(0);
  fOutput->Add(fHistMissPMod);

  fHistSkippedMod  = new TH1F("hSkippedmod","Tracks in Skipped Module",260,239.5,499.5);
  fHistSkippedMod->Sumw2();
  fHistSkippedMod->SetMinimum(0);
  fOutput->Add(fHistSkippedMod);

  fHistOutAccMod  = new TH1F("hOutAccmod","Tracks outside zAcc per Module",260,239.5,499.5);
  fHistOutAccMod->Sumw2();
  fHistOutAccMod->SetMinimum(0);
  fOutput->Add(fHistOutAccMod);

  fHistNoRefitMod  = new TH1F("hNoRefitmod","Points rejected in refit per Module",260,239.5,499.5);
  fHistNoRefitMod->Sumw2();
  fHistNoRefitMod->SetMinimum(0);
  fOutput->Add(fHistNoRefitMod);



  fRecPMod = new TH1F("hRPMod","Rec Points per Module",260,239.5,499.5);
  fRecPMod->Sumw2();
  fRecPMod->SetMinimum(0);
  fOutput->Add(fRecPMod);

  fTrackPMod = new TH1F("hTPMod","Track Points per Module",260,239.5,499.5);
  fTrackPMod->Sumw2();
  fTrackPMod->SetMinimum(0);
  fOutput->Add(fTrackPMod);

  fGoodAnMod = new TH1F("hGAMod","Good Anodes per Module",260,239.5,499.5);
  fOutput->Add(fGoodAnMod);

  // -- Local coordinates

  fHistAllPXloc  = new TH1F("hAllPxloc","Crossing Tracks vs. Xloc",75, -3.75, 3.75);
  fHistAllPXloc->Sumw2();
  fHistAllPXloc->SetMinimum(0);
  fOutput->Add(fHistAllPXloc);

  fHistGoodPXloc  = new TH1F("hGoodPxloc","PointsAssocToTrack vs. Xloc",75, -3.75, 3.75);
  fHistGoodPXloc->Sumw2();
  fHistGoodPXloc->SetMinimum(0);
  fOutput->Add(fHistGoodPXloc);

  fHistBadRegXloc  = new TH1F("hBadRegxloc","Tracks in BadRegion vs. Xloc",75, -3.75, 3.75);
  fHistBadRegXloc->Sumw2();
  fHistBadRegXloc->SetMinimum(0);
  fOutput->Add(fHistBadRegXloc);

  fHistMissPXloc  = new TH1F("hMissPxloc","Missing Points vs. Xloc",75, -3.75, 3.75);
  fHistMissPXloc->Sumw2();
  fHistMissPXloc->SetMinimum(0);
  fOutput->Add(fHistMissPXloc);

  fHistAllPZloc  = new TH1F("hAllPzloc","Crossing Tracks vs. Zloc",77, -3.85, 3.85);
  fHistAllPZloc->Sumw2();
  fHistAllPZloc->SetMinimum(0);
  fOutput->Add(fHistAllPZloc);

  fHistGoodPZloc  = new TH1F("hGoodPzloc","PointsAssocToTrack vs. Zloc",77, -3.85, 3.85);
  fHistGoodPZloc->Sumw2();
  fHistGoodPZloc->SetMinimum(0);
  fOutput->Add(fHistGoodPZloc);

  fHistBadRegZloc  = new TH1F("hBadRegzloc","Tracks in BadRegion vs. Zloc",77, -3.85, 3.85);
  fHistBadRegZloc->Sumw2();
  fHistBadRegZloc->SetMinimum(0);
  fOutput->Add(fHistBadRegZloc);

  fHistMissPZloc  = new TH1F("hMissPzloc","Missing Points vs. Zloc",77, -3.85, 3.85);
  fHistMissPZloc->Sumw2();
  fHistMissPZloc->SetMinimum(0);
  fOutput->Add(fHistMissPZloc);

  // -- Ladder histos

  fRecPLadLay3 = new TH1F("hRPLad3","Rec Points per Ladder Layer 3",14,-0.5,13.5);
  fRecPLadLay3->Sumw2();
  fRecPLadLay3->SetMinimum(0);
  fOutput->Add(fRecPLadLay3);

  fRecPLadLay4 = new TH1F("hRPLad4","Rec Points per Ladder Layer 4",22,-0.5,21.5);
  fRecPLadLay4->Sumw2();
  fRecPLadLay4->SetMinimum(0);
  fOutput->Add(fRecPLadLay4);

  fTrackPLadLay3 = new TH1F("hTPLad3","Track Points per Ladder Layer 3",14,-0.5,13.5);
  fTrackPLadLay3->Sumw2();
  fTrackPLadLay3->SetMinimum(0);
  fOutput->Add(fTrackPLadLay3);

  fTrackPLadLay4 = new TH1F("hTPLad4","Track Points per Ladder Layer 4",22,-0.5,21.5);
  fTrackPLadLay4->Sumw2();
  fTrackPLadLay4->SetMinimum(0);
  fOutput->Add(fTrackPLadLay4);

  fGoodAnLadLay3 = new TH1F("hGALad3","Good Anodes per Ladder Layer 3",14,-0.5,13.5);
  fOutput->Add(fGoodAnLadLay3);

  fGoodAnLadLay4 = new TH1F("hGALad4","Good Anodes per Ladder Layer 4",22,-0.5,21.5);
  fOutput->Add(fGoodAnLadLay4);

  fDriftTimeRP=new TH1F("hDrTimRP","Drift Time from Rec Points (ns)",640,0.,6400.);
  fDriftTimeRP->Sumw2();
  fDriftTimeRP->SetMinimum(0.);
  fOutput->Add(fDriftTimeRP);

  fEtaPhiTracks=new TH2F("hEtaPhiTracks","",50,-1.,1.,200,0.,2.*TMath::Pi());
  fEtaPhiTracks->SetMinimum(0.);
  fOutput->Add(fEtaPhiTracks);

  fEtaPhiTracksLay3=new TH2F("hEtaPhiTracksLay3","",50,-1.,1.,200,0.,2.*TMath::Pi());
  fEtaPhiTracksLay3->SetMinimum(0.);
  fOutput->Add(fEtaPhiTracksLay3);

  fEtaPhiTracksLay4=new TH2F("hEtaPhiTracksLay4","",50,-1.,1.,200,0.,2.*TMath::Pi());
  fEtaPhiTracksLay4->SetMinimum(0.);
  fOutput->Add(fEtaPhiTracksLay4);

  fDriftTimeTPAll=new TH1F("hDrTimTPAll","Drift Time from Track Points (ns)",640,0.,6400.);
  fDriftTimeTPAll->Sumw2();
  fDriftTimeTPAll->SetMinimum(0.);
  fOutput->Add(fDriftTimeTPAll);

  fDriftTimeTPAllMod=new TH2F("hDrTimTPAllMod","Drift Time from Track Points (ns) vs. mod number lay ",260,239.5,499.5,640,0.,6400.);
  fOutput->Add(fDriftTimeTPAllMod);

  fDriftTimeTPNoExtra=new TH1F("hDrTimTPNoExtra","Drift Time from Track Points (ns)",640,0.,6400.);
  fDriftTimeTPNoExtra->Sumw2();
  fDriftTimeTPNoExtra->SetMinimum(0.);
  fOutput->Add(fDriftTimeTPNoExtra);

  fDriftTimeTPExtra=new TH1F("hDrTimTPExtra","Drift Time from Track Points (ns)",640,0.,6400.);
  fDriftTimeTPExtra->Sumw2();
  fDriftTimeTPExtra->SetMinimum(0.);
  fOutput->Add(fDriftTimeTPExtra);

  // dE/dx histos

  fHistdEdxL3VsP=new TH2F("hdEdxL3VsP","dE/dx vs. p lay3",40,0.,2.,100,0.,500.);
  fHistdEdxL3VsP->Sumw2();
  fHistdEdxL3VsP->SetMinimum(0);
  fOutput->Add(fHistdEdxL3VsP);

  fHistdEdxL4VsP=new TH2F("hdEdxL4VsP","dE/dx vs. p lay4",40,0.,2.,100,0.,500);
  fHistdEdxL4VsP->Sumw2();
  fHistdEdxL4VsP->SetMinimum(0);
  fOutput->Add(fHistdEdxL4VsP);

  fHistdEdxVsMod=new TH2F("hdEdxVsMod","dE/dx vs. mod",260,239.5,499.5,100,0.,500.);
  fHistdEdxVsMod->Sumw2();
  fHistdEdxVsMod->SetMinimum(0);
  fOutput->Add(fHistdEdxVsMod);

  for(Int_t it=0; it<8; it++){
    fSignalTime[it]=new TH1F(Form("hSigTimeInt%d",it),Form("hSigTimeInt%d",it),100,0.,300.);
    fSignalTime[it]->Sumw2();
    fSignalTime[it]->SetMinimum(0);
    fOutput->Add(fSignalTime[it]);
  }

  // cluster size histos
  
  fCluSizAnVsTime = new TH2F("hCluSizAn","hCluSizAn",40,0.,6400.,15,-0.5,14.5);
  fCluSizAnVsTime->Sumw2();
  fCluSizAnVsTime->SetMinimum(0);
  fOutput->Add(fCluSizAnVsTime);

  fCluSizTbVsTime = new TH2F("hCluSizTb","hCluSizTb",40,0.,6400.,15,-0.5,14.5);
  fCluSizTbVsTime->Sumw2();
  fCluSizTbVsTime->SetMinimum(0);
  fOutput->Add(fCluSizTbVsTime);


  // profiles of rec and track points vs. time in run
  fProfRecPtsLay3VsTime = new TProfile("profRecPtsLay3VsTime"," ; time (sec) ; nRecPts Layer3",500,0.,1000.*60);
  fProfRecPtsLay4VsTime = new TProfile("profRecPtsLay4VsTime"," ; time (sec) ; nRecPts Layer4",500,0.,1000.*60);
  fOutput->Add(fProfRecPtsLay3VsTime);
  fOutput->Add(fProfRecPtsLay4VsTime);
  fProfFracTrkWithPntLay3VsTime = new TProfile("profFracTrkWithPntLay3VsTime"," ; time (sec) ; frac of tracks with point in Layer3",500,0.,1000.*60);
  fProfFracTrkWithPntLay4VsTime = new TProfile("profFracTrkWithPntLay4VsTime"," ; time (sec) ; frac of tracks with point in Layer4",500,0.,1000.*60);
  fOutput->Add(fProfFracTrkWithPntLay3VsTime);
  fOutput->Add(fProfFracTrkWithPntLay4VsTime);
  fProfTrPtsLay3VsTime = new TProfile("profTrPtsLay3VsTime"," ; time (sec) ; nTrPts Layer3",500,0.,1000.*60);
  fProfTrPtsLay4VsTime = new TProfile("profTrPtsLay4VsTime"," ; time (sec) ; nTrPts Layer4",500,0.,1000.*60);
  fOutput->Add(fProfTrPtsLay3VsTime);
  fOutput->Add(fProfTrPtsLay4VsTime);
  fProfFracTrRecLay3VsTime = new TProfile("profFracTrRecLay3VsTime"," ; time (sec) ; nTrPts/nRecPts Layer3",500,0.,1000.*60);
  fProfFracTrRecLay4VsTime = new TProfile("profFracTrRecLay4VsTime"," ; time (sec) ; nTrPts/nRecPts Layer4",500,0.,1000.*60);
  fOutput->Add(fProfFracTrRecLay3VsTime);
  fOutput->Add(fProfFracTrRecLay4VsTime);

    // Read dead channels from OCDB

// == Not allowed in QA train. This is set by AliTaskCDBconnect. (A.G. 14/10/2011)
// If needed for independent use of the task, protect via a flag that is OFF by default.
/*
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) AliFatal("No analysis manager");

  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  Int_t nrun=mgr->GetRunFromPath();
  man->SetRun(nrun);
*/  
    
  PostData(1, fOutput);

}
//______________________________________________________________________________
void AliAnalysisTaskSDDRP::UserExec(Option_t *)
{
  //
  AliESDEvent *esd = (AliESDEvent*) (InputEvent());
  if(!esd) {
    printf("AliAnalysisTaskSDDRP::UserExec(): bad ESD\n");
    return;
  } 


  if(!ESDfriend()) {
    printf("AliAnalysisTaskSDDRP::UserExec(): bad ESDfriend\n");
    return;
  }
  
// Code below moved from UserCreateOutputObjects where the run for OCDB may not 
// be yet properly set. Make sure this is called only once. (A.G. 14/10/2011)

/************/

  if (fReadCDB && !fInitCalib) {
    AliCDBManager* man = AliCDBManager::Instance();
    if (!man) {
       AliFatal("CDB not set but needed by AliAnalysisTaskSDDRP");
       return;
    }   
    AliCDBEntry* eT=(AliCDBEntry*)man->Get("GRP/CTP/Config");
    if(eT){
      eT->PrintId();
      eT->PrintMetaData();
      fTrigConfig=(AliTriggerConfiguration*)eT->GetObject();
    }
    if(!eT || !fTrigConfig){
      AliError("Cannot retrieve CDB entry for GRP/CTP/Config");
      return;      
    }
    AliCDBEntry* eR=(AliCDBEntry*)man->Get("ITS/Calib/RespSDD");
    if (eR) {
      eR->PrintId();
      eR->PrintMetaData();
      fResp=(AliITSresponseSDD*)eR->GetObject();
    }else{
      AliError("Cannot retrieve CDB entry for ITS/Calib/RespSDD");
      return;
    }
    AliCDBEntry* eC=(AliCDBEntry*)man->Get("ITS/Calib/CalibSDD");
    if (!eC) {
       AliError("Cannot retrieve CDB entry for ITS/Calib/CalibSDD");
       return;
    }   
    Int_t countGood3[14];
    Int_t countGood4[22];
    Int_t countGoodMod[260];
    eC->PrintId();
    eC->PrintMetaData();
    TObjArray* calsdd=(TObjArray*)eC->GetObject();
    for(Int_t ilad=0;ilad<14;ilad++) countGood3[ilad]=0;
    for(Int_t ilad=0;ilad<22;ilad++) countGood4[ilad]=0;
    for(Int_t imod=0;imod<260;imod++) countGoodMod[imod]=0;
    for(Int_t imod=0;imod<260;imod++){
      AliITSCalibrationSDD* cal=(AliITSCalibrationSDD*)calsdd->At(imod);
      if(cal->IsBad()) continue;
	   Int_t modId=imod+AliITSgeomTGeo::GetModuleIndex(3,1,1);
      Int_t lay,lad,det;
      AliITSgeomTGeo::GetModuleId(modId,lay,lad,det);
      if(fExcludeBadMod && !CheckModule(lay,lad,det)) continue;
      for(Int_t ian=0; ian<512; ian++){
        if(cal->IsBadChannel(ian)) continue;
        countGoodMod[imod]++;
        if(lay==3) countGood3[lad-1]++;
        else if(lay==4) countGood4[lad-1]++;
      }
    }
  
    for(Int_t imod=0;imod<260;imod++) fGoodAnMod->SetBinContent(imod+1,countGoodMod[imod]);
    fGoodAnMod->SetMinimum(0);
    for(Int_t ilad=0;ilad<14;ilad++) fGoodAnLadLay3->SetBinContent(ilad+1,countGood3[ilad]);
    fGoodAnLadLay3->SetMinimum(0);    
    for(Int_t ilad=0;ilad<22;ilad++) fGoodAnLadLay4->SetBinContent(ilad+1,countGood4[ilad]);
    fGoodAnLadLay4->SetMinimum(0);
    fInitCalib = kTRUE;
  }  
/************/
  
  fProfRecPtsLay3VsTime->SetName(Form("profRecPtsLay3VsTimeRun%d",esd->GetRunNumber()));
  fProfRecPtsLay4VsTime->SetName(Form("profRecPtsLay4VsTimeRun%d",esd->GetRunNumber()));
  fProfTrPtsLay3VsTime->SetName(Form("profTrPtsLay3VsTimeRun%d",esd->GetRunNumber()));
  fProfTrPtsLay4VsTime->SetName(Form("profTrPtsLay4VsTimeRun%d",esd->GetRunNumber()));
  fProfFracTrRecLay3VsTime->SetName(Form("profFracTrRecLay3VsTimeRun%d",esd->GetRunNumber()));
  fProfFracTrRecLay4VsTime->SetName(Form("profFracTrRecLay4VsTimeRun%d",esd->GetRunNumber()));
  fProfFracTrkWithPntLay3VsTime->SetName(Form("profFracTrkWithPntLay3VsTime%d",esd->GetRunNumber()));
  fProfFracTrkWithPntLay4VsTime->SetName(Form("profFracTrkWithPntLay4VsTime%d",esd->GetRunNumber()));
  PostData(1, fOutput);

  fHistNEvents->Fill(-1);

  TString firedTriggerClasses=esd->GetFiredTriggerClasses();
  if(!firedTriggerClasses.Contains(fTriggerClass.Data())) return;
  fHistNEvents->Fill(0);

  if(fTrigConfig){
    Bool_t sddOK=esd->IsDetectorInTriggerCluster("ITSSDD",fTrigConfig);
    if(!sddOK) fHistNEvents->Fill(1.);
    else fHistNEvents->Fill(2.);
    if(fOnlyEventsWithSDD && !sddOK) return;
  }
  fHistNEvents->Fill(3.);

  const AliTimeStamp* t0=esd->GetCTPStart();
  UInt_t t0sec=t0->GetSeconds();
  AliTimeStamp tev=esd->GetAliTimeStamp();
  UInt_t tevsec=tev.GetSeconds();
  tevsec-=t0sec;

  const AliMultiplicity *alimult = esd->GetMultiplicity();
  Int_t nRecPtsLay3=0;
  Int_t nRecPtsLay4=0;
  if(alimult){
    nRecPtsLay3 = alimult->GetNumberOfITSClusters(2);
    nRecPtsLay4 = alimult->GetNumberOfITSClusters(3);
  }
  fProfRecPtsLay3VsTime->Fill(tevsec,nRecPtsLay3);
  fProfRecPtsLay4VsTime->Fill(tevsec,nRecPtsLay4);

  const AliTrackPointArray *array = 0;
  Int_t ntracks = esd->GetNumberOfTracks();
  Int_t nTrPtsLay3=0;
  Int_t nTrPtsLay4=0;
  Int_t nTrWithPtIn3=0;
  Int_t nTrWithPtIn4=0;
  Int_t nAccTr=0;
  for (Int_t itrack=0; itrack < ntracks; itrack++) {
    AliESDtrack * track = esd->GetTrack(itrack);
    if (!track) continue;

    Bool_t accept=kTRUE;
    if(fUseITSsaTracks){ 
      if(track->GetNcls(1)>0) accept=kFALSE;
    }else{
      if(track->GetNcls(1)<fMinTPCpts) accept=kFALSE;
    }
    if(track->GetNcls(0) < fMinITSpts) accept=kFALSE;    
    Int_t trstatus=track->GetStatus();
    if(!(trstatus&AliESDtrack::kITSrefit)) accept=kFALSE;
    if(!accept) continue;
    Double_t eta=track->Eta();
    Double_t phi=track->Phi();
    
    fHistCluInLay->Fill(-1.); // bin -1 counts accepted tracks
    fEtaPhiTracks->Fill(eta,phi);
    UChar_t clumap=track->GetITSClusterMap();
    for(Int_t iBit=0; iBit<6; iBit++){
      if(clumap&(1<<iBit)) fHistCluInLay->Fill(iBit);
    }
    nAccTr++;
    if(clumap&(1<<2)){
      fEtaPhiTracksLay3->Fill(eta,phi);
      nTrWithPtIn3++;
    }
    if(clumap&(1<<3)){
      fEtaPhiTracksLay4->Fill(eta,phi);
      nTrWithPtIn4++;
    }

    Double_t dedx[4];
    track->GetITSdEdxSamples(dedx);
    Float_t mom=track->P();
    Int_t iMod,status;
    Float_t xloc,zloc;
    for(Int_t iLay=2; iLay<=3; iLay++){
      Bool_t ok=track->GetITSModuleIndexInfo(iLay,iMod,status,xloc,zloc);
      if(ok){
	iMod+=240;
	fHistAllPMod->Fill(iMod);
	fHistAllPXloc->Fill(xloc);
	fHistAllPZloc->Fill(zloc);
	if(status==1){
	  fHistGoodPMod->Fill(iMod);
	  fHistGoodPXloc->Fill(xloc);
	  fHistGoodPZloc->Fill(zloc);
	  if(mom>fMinPfordEdx) fHistdEdxVsMod->Fill(iMod,dedx[iLay-2]);
	  if(iLay==2) fHistdEdxL3VsP->Fill(mom,dedx[0]);
	  else fHistdEdxL4VsP->Fill(mom,dedx[1]);
	}
	else if(status==2){ 
	  fHistBadRegMod->Fill(iMod);
	  fHistBadRegXloc->Fill(xloc);
	  fHistBadRegZloc->Fill(zloc);
	}
	else if(status==3) fHistSkippedMod->Fill(iMod);
	else if(status==4) fHistOutAccMod->Fill(iMod);
	else if(status==5){
	  fHistMissPMod->Fill(iMod);
	  fHistMissPXloc->Fill(xloc);
	  fHistMissPZloc->Fill(zloc);
	}
	else if(status==6) fHistNoRefitMod->Fill(iMod);
      }
    }


    array = track->GetTrackPointArray();
    if(!array) continue;
    for(Int_t ipt=0; ipt<array->GetNPoints(); ipt++) {
      AliTrackPoint point;
      Int_t modId;
      array->GetPoint(point,ipt);
      Int_t volId = point.GetVolumeID();
      Int_t layerId = AliGeomManager::VolUIDToLayer(volId,modId);
      if(layerId<3 || layerId>4) continue;
      if(layerId==3) nTrPtsLay3++;
      else if(layerId==4) nTrPtsLay4++;
      modId+=AliITSgeomTGeo::GetModuleIndex(layerId,1,1);
      Int_t lay,lad,det;
      AliITSgeomTGeo::GetModuleId(modId,lay,lad,det);
      if(fExcludeBadMod && !CheckModule(lay,lad,det)) continue;
      fTrackPMod->Fill(modId);
      fDriftTimeTPAll->Fill(point.GetDriftTime());
      fDriftTimeTPAllMod->Fill(modId,point.GetDriftTime());
      if(point.IsExtra()) fDriftTimeTPExtra->Fill(point.GetDriftTime());
      else fDriftTimeTPNoExtra->Fill(point.GetDriftTime());
      Float_t dtime=point.GetDriftTime();
      if(fResp) dtime=point.GetDriftTime()-fResp->GetTimeZero(modId);
      Int_t cluTyp=point.GetClusterType();
      Int_t clSizAn=(cluTyp>>8)&0xFF;
      Int_t clSizTb=cluTyp&0xFF;
      fCluSizAnVsTime->Fill(dtime,clSizAn);
      fCluSizTbVsTime->Fill(dtime,clSizTb);
      Int_t theBin=int(dtime/6500.*8.);
      if(layerId==3){
	fTrackPLadLay3->Fill(lad-1);	  
	if(dedx[0]>0. && track->P()>fMinPfordEdx) fSignalTime[theBin]->Fill(dedx[0]);
      }
      if(layerId==4){
	fTrackPLadLay4->Fill(lad-1);
	if(dedx[1]>0.&& track->P()>fMinPfordEdx) fSignalTime[theBin]->Fill(dedx[1]);
      }
    }
  }

  fProfTrPtsLay3VsTime->Fill(tevsec,nTrPtsLay3);
  fProfTrPtsLay4VsTime->Fill(tevsec,nTrPtsLay4);
  if(nRecPtsLay3>0){
    Double_t frac3=(Double_t)nTrPtsLay3/(Double_t)nRecPtsLay3;
    fProfFracTrRecLay3VsTime->Fill(tevsec,frac3);
  }
  if(nRecPtsLay4>0){
    Double_t frac4=(Double_t)nTrPtsLay4/(Double_t)nRecPtsLay4;
    fProfFracTrRecLay4VsTime->Fill(tevsec,frac4);
  }
  if(nAccTr>0){
    Double_t fracTrPt3=(Double_t)nTrWithPtIn3/(Double_t)nAccTr;
    Double_t fracTrPt4=(Double_t)nTrWithPtIn4/(Double_t)nAccTr;
    fProfFracTrkWithPntLay3VsTime->Fill(tevsec,fracTrPt3);
    fProfFracTrkWithPntLay4VsTime->Fill(tevsec,fracTrPt4);
  }

  AliESDInputHandlerRP *hand = dynamic_cast<AliESDInputHandlerRP*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  TTree* tR = 0;
  if (hand) tR = hand->GetTreeR("ITS");
  if (tR){
    TClonesArray *ITSrec= new TClonesArray("AliITSRecPoint");
    TBranch *branch =tR->GetBranch("ITSRecPoints");
    branch->SetAddress(&ITSrec);
    for (Int_t modId=240; modId<500; modId++){
      Int_t lay,lad,det;
      AliITSgeomTGeo::GetModuleId(modId,lay,lad,det);
      if(fExcludeBadMod && !CheckModule(lay,lad,det)) continue;
      branch->GetEvent(modId);
      Int_t nrecp = ITSrec->GetEntries();	
      fRecPMod->Fill(modId,nrecp);	  
      if(lay==3) fRecPLadLay3->Fill(lad-1,nrecp);
      if(lay==4) fRecPLadLay4->Fill(lad-1,nrecp);
      for(Int_t irec=0;irec<nrecp;irec++) {
	AliITSRecPoint *recp = (AliITSRecPoint*)ITSrec->At(irec);
	fDriftTimeRP->Fill(recp->GetDriftTime());
      }
    }
    ITSrec->Delete();
    delete ITSrec;
  }
  PostData(1,fOutput);
  
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskSDDRP::CheckModule(Int_t lay, Int_t lad, Int_t det) const{
  //
  if(lay==4){
    if(lad==3 && det==5) return kFALSE; // 1500 V
    if(lad==3 && det==6) return kFALSE; // 0 V
    if(lad==3 && det==7) return kFALSE; // 1500 V
    if(lad==4 && det==1) return kFALSE; // 0 V
    if(lad==4 && det==2) return kFALSE; // 1500 V
    if(lad==7 && det==5) return kFALSE; // 0 MV
    if(lad==9 && det==3) return kFALSE; // 1500 V
    if(lad==9 && det==4) return kFALSE; // 0 V
    if(lad==9 && det==5) return kFALSE; // 1500 V
    if(lad==11 && det==6) return kFALSE; // 1500 V
    if(lad==11 && det==7) return kFALSE; // 0 V
    if(lad==11 && det==8) return kFALSE; // 1500 V
    if(lad==18 && det==5) return kFALSE; // 1500 V
    if(lad==18 && det==6) return kFALSE; // 0 V
    if(lad==18 && det==7) return kFALSE; // 1500 V
    if(lad==22 && det==1) return kFALSE; // 0 V
    if(lad==22 && det==2) return kFALSE; // 1500 V
  }
  if(lay==3){
    if(lad==4 && det==4) return kFALSE; // 1500 V 
    if(lad==3) return kFALSE;  // swapped in geometry
  }
  return kTRUE;
}

//______________________________________________________________________________
void AliAnalysisTaskSDDRP::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
  fHistNEvents= dynamic_cast<TH1F*>(fOutput->FindObject("hNEvents"));
  printf("AliAnalysisTaskSDDRP::Terminate --- Number of analyzed events = %.0f\n",fHistNEvents->GetBinContent(1));
  return;
}





