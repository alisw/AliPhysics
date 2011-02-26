
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////////////////////////
//                                                                             
// AliTRDCalibTask                                                               
//                                                                             
// Offline TRD calibration task                                
//                        
// Author:
//   R. Bailhache (R.Bailhache@gsi.de)
//                            
//////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
using namespace std;
#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TH2I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
#include "TMath.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TIterator.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliExternalTrackParam.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliCentrality.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDfriendTrack.h"
#include "AliTRDtrackV1.h"
#include "AliTRDseedV1.h"
#include "AliTRDcluster.h"
#include "AliTRDgeometry.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliTRDCalDet.h"

#include "AliTRDCalibraVector.h"
#include "AliTRDCalibraFillHisto.h"
#include "AliTRDCalibraVdriftLinearFit.h" 

#include "AliTRDcalibDB.h"
#include "AliCDBId.h"
#include "AliLog.h"


#include "AliTRDCalibTask.h"


ClassImp(AliTRDCalibTask)

//________________________________________________________________________
  AliTRDCalibTask::AliTRDCalibTask(const char *name) 
    : AliAnalysisTaskSE(name), fESD(0),
      fESDfriend(0),
      fkEsdTrack(0),
      fFriendTrack(0),
      fCalibObject(0),
      fTrdTrack(0),
      fCl(0),
      fListHist(0),
      fTRDCalibraFillHisto(0),
      fNEvents(0),
      fNEventsInput(0),
      fNbTRDTrack(0),
      fNbTRDTrackOffline(0),
      fNbTRDTrackStandalone(0),
      fNbTPCTRDtrack(0),
      fNbGoodTracks(0),
      fNbTimeBin(0),
      fNbTimeBinOffline(0),
      fNbTimeBinStandalone(0),
      fNbClusters(0),
      fNbClustersOffline(0),
      fNbClustersStandalone(0),
      fNbTracklets(0),
      fNbTrackletsOffline(0),
      fNbTrackletsStandalone(0),
      fAbsoluteGain(0),
      fCH2dSum(0),
      fPH2dSum(0),
      fCH2dSM(0),
      fPH2dSM(0),
      fHisto2d(kTRUE),
      fVector2d(kFALSE),
      fVdriftLinear(kTRUE),
      fNbTimeBins(0),
      fSelectedTrigger(new TObjArray()),
      fRejected(kTRUE),
      fEsdTrackCuts(0),
      fRequirePrimaryVertex(kFALSE),
      fVtxTPC(kFALSE),
      fVtxSPD(kFALSE),
      fMinNbContributors(0),
      fRangePrimaryVertexZ(9999999.0),
      fMinNbTracks(9),
      fMaxNbTracks(500),
      fLow(0),
      fHigh(30),
      fFillZero(kFALSE),
      fNormalizeNbOfCluster(kFALSE),
      fRelativeScale(0.0),
      fMaxCluster(100.0),
      fNbMaxCluster(2),
      fOfflineTracks(kFALSE),
      fStandaloneTracks(kFALSE),
      fFirstRunGain(-1),
      fVersionGainUsed(-1),
      fSubVersionGainUsed(-1),
      fFirstRunGainLocal(-1),
      fVersionGainLocalUsed(-1),
      fSubVersionGainLocalUsed(-1),
      fFirstRunVdrift(-1),
      fVersionVdriftUsed(-1), 
      fSubVersionVdriftUsed(-1),
      fCalDetGain(0x0),
      fMaxEvent(0),
      fCounter(0),
      fDebug(0)
{
  //
  // Default constructor
  //

  fNz[0] = 0;
  fNz[1] = 0;
  fNz[2] = 0;
  
  fNrphi[0] = 0;
  fNrphi[1] = 0;
  fNrphi[2] = 0;

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
        
  // Output slot #0 writes into a TList container
  DefineOutput(1, TList::Class());
  
   
}
//____________________________________________________________________________________
AliTRDCalibTask::~AliTRDCalibTask()
{
  //
  // AliTRDCalibTask destructor
  //

  // Pointeur
  if(fNEvents) delete fNEvents;
  if(fNEventsInput) delete fNEventsInput;
  if(fNbTRDTrack) delete fNbTRDTrack;
  if(fNbTRDTrackOffline) delete fNbTRDTrackOffline;
  if(fNbTRDTrackStandalone) delete fNbTRDTrackStandalone;
  if(fNbTPCTRDtrack) delete fNbTPCTRDtrack;
  if(fNbGoodTracks) delete fNbGoodTracks;
  if(fNbTimeBin) delete fNbTimeBin;
  if(fNbTimeBinOffline) delete fNbTimeBinOffline;
  if(fNbTimeBinStandalone) delete fNbTimeBinStandalone;
  if(fNbClusters) delete fNbClusters;
  if(fNbClustersOffline) delete fNbClustersOffline;
  if(fNbClustersStandalone) delete fNbClustersStandalone;
  if(fNbTracklets) delete fNbTracklets;
  if(fNbTrackletsOffline) delete fNbTrackletsOffline;
  if(fNbTrackletsStandalone) delete fNbTrackletsStandalone;
  if(fAbsoluteGain) delete fAbsoluteGain;
  if(fCH2dSum) delete fCH2dSum;
  if(fPH2dSum) delete fPH2dSum;
  if(fCH2dSM) delete fCH2dSM;
  if(fPH2dSM) delete fPH2dSM;
  if(fCalDetGain) delete fCalDetGain;
  
  if(fSelectedTrigger) {
    fSelectedTrigger->Delete();
    delete fSelectedTrigger;
  }
  if(fEsdTrackCuts) {
    delete fEsdTrackCuts;
  }
  
}

//________________________________________________________________________
void AliTRDCalibTask::UserCreateOutputObjects() 
{
  //
  // Create the histos
  //
  //cout << "AliTRDCalibTask::CreateOutputObjects() IN" << endl;

  // Number of time bins
  if(fNbTimeBins==0) {
    AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
    fNbTimeBins = cal->GetNumberOfTimeBinsDCS();
    if(fNbTimeBins <= 0){ 
      AliWarning(Form("No of TimeBins from DB [%d] use default [30]", fNbTimeBins));
      fNbTimeBins = 30;
    }
  }

  // instance calibration 
  fTRDCalibraFillHisto = AliTRDCalibraFillHisto::Instance();
  fTRDCalibraFillHisto->SetHisto2d(fHisto2d); // choose to use histograms
  fTRDCalibraFillHisto->SetVector2d(fVector2d); // choose to use vectors
  fTRDCalibraFillHisto->SetCH2dOn();  // choose to calibrate the gain
  fTRDCalibraFillHisto->SetPH2dOn();  // choose to calibrate the drift velocity
  fTRDCalibraFillHisto->SetPRF2dOn(); // choose to look at the PRF
  fTRDCalibraFillHisto->SetLinearFitterOn(fVdriftLinear); // Other possibility vdrift VDRIFT
  fTRDCalibraFillHisto->SetLinearFitterDebugOn(fVdriftLinear); // Other possibility vdrift
  for(Int_t k = 0; k < 3; k++){
    if(((fNz[k] != 10) && (fNrphi[k] != 10)) && ((fNz[k] != 100) && (fNrphi[k] != 100))) {
      fTRDCalibraFillHisto->SetNz(k,fNz[k]);                                    // Mode calibration
      fTRDCalibraFillHisto->SetNrphi(k,fNrphi[k]);                             // Mode calibration
    }
    else {
      if((fNz[k] == 100) && (fNrphi[k] == 100))  {
	if(fVector2d) AliInfo("The mode all together is not supported by the vector method");
	fTRDCalibraFillHisto->SetAllTogether(k);
      }
      if((fNz[k] == 10) && (fNrphi[k] == 10))  {
	if(fVector2d) AliInfo("The mode per supermodule is not supported by the vector method");
	fTRDCalibraFillHisto->SetPerSuperModule(k);
      }
    }
  }
  // Variables for how to fill
  fTRDCalibraFillHisto->SetFillWithZero(fFillZero);
  fTRDCalibraFillHisto->SetNormalizeNbOfCluster(fNormalizeNbOfCluster); 
  fTRDCalibraFillHisto->SetMaxCluster(fMaxCluster);
  fTRDCalibraFillHisto->SetNbMaxCluster(fNbMaxCluster);
  
  // Init with 30 timebins
  fTRDCalibraFillHisto->Init2Dhistos(fNbTimeBins); // initialise the histos
  fTRDCalibraFillHisto->SetNumberClusters(fLow); // At least 11 clusters
  fTRDCalibraFillHisto->SetNumberClustersf(fHigh); // At least 11 clusters
  fRelativeScale = fTRDCalibraFillHisto->GetRelativeScale(); // Get the relative scale for the gain
  
  // For testing only
  if(fDebug > 2) fTRDCalibraFillHisto->SetDebugLevel(1); //debug stuff
  
  // output list
  fListHist = new TList();
  fListHist->SetOwner();
  if(fHisto2d) {  
    fListHist->Add(fTRDCalibraFillHisto->GetCH2d());
    fListHist->Add(fTRDCalibraFillHisto->GetPH2d()); 
    fListHist->Add(fTRDCalibraFillHisto->GetPRF2d());
  } 
  if(fVdriftLinear) fListHist->Add((TObject *) fTRDCalibraFillHisto->GetVdriftLinearFit());
  if(fVector2d) fListHist->Add((TObject *) fTRDCalibraFillHisto->GetCalibraVector()); //calibra vector  
  fNEvents = new TH1I("NEvents","NEvents", 2, 0, 2);
  fListHist->Add(fNEvents);
  fNEventsInput = new TH1I("NEventsInput","NEventsInput", 2, 0, 2);
  fListHist->Add(fNEventsInput);
  
  // absolute gain calibration even without AliESDfriend
  Int_t nBinsPt = 25;
  Double_t minPt = 0.001;
  Double_t maxPt = 10.0;
  
  Double_t *binLimLogPt = new Double_t[nBinsPt+1];
  Double_t *binLimPt    = new Double_t[nBinsPt+1];
  for(Int_t i=0; i<=nBinsPt; i++) binLimLogPt[i]=(Double_t)TMath::Log10(minPt) + (TMath::Log10(maxPt)-TMath::Log10(minPt))/nBinsPt*(Double_t)i ;
  for(Int_t i=0; i<=nBinsPt; i++) binLimPt[i]=(Double_t)TMath::Power(10,binLimLogPt[i]);
  
  fAbsoluteGain = new TH2F("AbsoluteGain","AbsoluteGain", 200, 0.0, 700.0, nBinsPt, binLimPt);
  fAbsoluteGain->SetYTitle("Momentum at TRD");
  fAbsoluteGain->SetXTitle("charge deposit [a.u]");
  fAbsoluteGain->SetZTitle("counts");
  fAbsoluteGain->SetStats(0);
  fAbsoluteGain->Sumw2();
  fListHist->Add(fAbsoluteGain);
  

  
  /////////////////////////////////////////
  // First debug level
  ///////////////////////////////////////
  if(fDebug > 0) {
    
    // Standart with AliESDfriend
    fPH2dSM = new TProfile2D("PH2dSM","Nz10Nrphi10"
			    ,fNbTimeBins,-0.05,(Double_t)((fNbTimeBins-0.5)/10.0)
			   ,18,0,18);
    fPH2dSM->SetYTitle("Det/pad groups");
    fPH2dSM->SetXTitle("time [#mus]");
    fPH2dSM->SetZTitle("<PH> [a.u.]");
    fPH2dSM->SetStats(0);
    //
    fCH2dSM = new TH2I("CH2dSM","Nz10Nrphi10",50,0,300,18,0,18);
    fCH2dSM->SetYTitle("Det/pad groups");
    fCH2dSM->SetXTitle("charge deposit [a.u]");
    fCH2dSM->SetZTitle("counts");
    fCH2dSM->SetStats(0);
    fCH2dSM->Sumw2();
    //
    fPH2dSum = new TProfile2D("PH2dSum","Nz100Nrphi100"
			    ,fNbTimeBins,-0.05,(Double_t)((fNbTimeBins-0.5)/10.0)
			    ,1,0,1);
    fPH2dSum->SetYTitle("Det/pad groups");
    fPH2dSum->SetXTitle("time [#mus]");
    fPH2dSum->SetZTitle("<PH> [a.u.]");
    fPH2dSum->SetStats(0);
    //
    fCH2dSum = new TH2I("CH2dSum","Nz100Nrphi100",50,0,300,1,0,1);
    fCH2dSum->SetYTitle("Det/pad groups");
    fCH2dSum->SetXTitle("charge deposit [a.u]");
    fCH2dSum->SetZTitle("counts");
    fCH2dSum->SetStats(0);
    fCH2dSum->Sumw2();
    //
    fNbGoodTracks = new TH2F("NbGoodTracks","NbGoodTracks",500,0.0,2500.0,200,0.0,100.0);
    fNbGoodTracks->SetXTitle("Nb of good tracks");
    fNbGoodTracks->SetYTitle("Centrality");
    fNbGoodTracks->SetStats(0);

    
    // Add them
    fListHist->Add(fPH2dSM);
    fListHist->Add(fCH2dSM);
    fListHist->Add(fPH2dSum);
    fListHist->Add(fCH2dSum);
    fListHist->Add(fNbGoodTracks);
  }

  /////////////////////////////////////////
  // Second debug level
  ///////////////////////////////////////
  if(fDebug > 1) {

    fNbTRDTrack = new TH1F("TRDTrack","TRDTrack",50,0,50);
    fNbTRDTrack->Sumw2();
    fNbTRDTrackOffline = new TH1F("TRDTrackOffline","TRDTrackOffline",50,0,50);
    fNbTRDTrackOffline->Sumw2();
    fNbTRDTrackStandalone = new TH1F("TRDTrackStandalone","TRDTrackStandalone",50,0,50);
    fNbTRDTrackStandalone->Sumw2();
    fNbTPCTRDtrack = new TH2F("NbTPCTRDtrack","NbTPCTRDtrack",100,0,100,100,0,100);
    fNbTPCTRDtrack->Sumw2();
    //
    fNbTimeBin = new TH1F("NbTimeBin","NbTimeBin",35,0,35);
    fNbTimeBin->Sumw2();
    fNbTimeBinOffline = new TH1F("NbTimeBinOffline","NbTimeBinOffline",35,0,35);
    fNbTimeBinOffline->Sumw2();
    fNbTimeBinStandalone = new TH1F("NbTimeBinStandalone","NbTimeBinStandalone",35,0,35);
    fNbTimeBinStandalone->Sumw2();
    //
    fNbClusters = new TH1F("NbClusters","",35,0,35);
    fNbClusters->Sumw2();
    fNbClustersOffline = new TH1F("NbClustersOffline","",35,0,35);
    fNbClustersOffline->Sumw2();
    fNbClustersStandalone = new TH1F("NbClustersStandalone","",35,0,35);
    fNbClustersStandalone->Sumw2();
    //
    fNbTracklets = new TH1F("NbTracklets","NbTracklets",540,0.,540.);
    fNbTracklets->Sumw2();
    fNbTrackletsOffline = new TH1F("NbTrackletsOffline","NbTrackletsOffline",540,0.,540.);
    fNbTrackletsOffline->Sumw2();
    fNbTrackletsStandalone = new TH1F("NbTrackletsStandalone","NbTrackletsStandalone",540,0.,540.);
    fNbTrackletsStandalone->Sumw2();
   
    fListHist->Add(fNbTRDTrack);
    fListHist->Add(fNbTRDTrackOffline);
    fListHist->Add(fNbTRDTrackStandalone);
    fListHist->Add(fNbTPCTRDtrack);
    
    fListHist->Add(fNbTimeBin);
    fListHist->Add(fNbTimeBinOffline);
    fListHist->Add(fNbTimeBinStandalone);
    fListHist->Add(fNbClusters);
    fListHist->Add(fNbClustersOffline);
    fListHist->Add(fNbClustersStandalone);
    fListHist->Add(fNbTracklets);
    fListHist->Add(fNbTrackletsOffline);
    fListHist->Add(fNbTrackletsStandalone);
    
  }

  delete [] binLimLogPt;
  delete [] binLimPt;

  PostData(1,fListHist);

  //cout << "AliTRDCalibTask::UserCreateOutputObjects() OUT" << endl;

}

//________________________________________________________________________
void AliTRDCalibTask::UserExec(Option_t *) 
{
  //
  // Filling of the histos
  //
  //cout << "AliTRDCalibTask::Exec() IN" << endl;
  
  // Init Versions and subversions used
  if((fFirstRunGain==-1) || (fVersionGainUsed==-1) || (fSubVersionGainUsed==-1) || (fFirstRunGainLocal==-1) || (fVersionGainLocalUsed==-1) || (fSubVersionGainLocalUsed==-1) || (fFirstRunVdrift==-1) || (fVersionVdriftUsed==-1) || (fSubVersionVdriftUsed==-1)) {
    if(!SetVersionSubversion()) {
      PostData(1, fListHist);
      return;
    }
  }
  if(fCounter==0) {
    fTRDCalibraFillHisto->SetFirstRunGain(fFirstRunGain); // Gain Used
    fTRDCalibraFillHisto->SetVersionGainUsed(fVersionGainUsed); // Gain Used
    fTRDCalibraFillHisto->SetSubVersionGainUsed(fSubVersionGainUsed); // Gain Used
    fTRDCalibraFillHisto->SetFirstRunGainLocal(fFirstRunGainLocal); // Gain Used
    fTRDCalibraFillHisto->SetVersionGainLocalUsed(fVersionGainLocalUsed); // Gain Used
    fTRDCalibraFillHisto->SetSubVersionGainLocalUsed(fSubVersionGainLocalUsed); // Gain Used
    fTRDCalibraFillHisto->SetFirstRunVdrift(fFirstRunVdrift); // Vdrift Used
    fTRDCalibraFillHisto->SetVersionVdriftUsed(fVersionVdriftUsed); // Vdrift Used
    fTRDCalibraFillHisto->SetSubVersionVdriftUsed(fSubVersionVdriftUsed); // Vdrift Used
    fTRDCalibraFillHisto->InitCalDet();
  }
  
  //  AliLog::SetGlobalLogLevel(AliLog::kError);
  //  cout << "AliTRDCalibTask::Exec() 1" << endl;
  fESD = dynamic_cast<AliESDEvent*>(fInputEvent);
  if(!fESD){
    AliError("ESD Event missing");
    PostData(1, fListHist);
    return;
  }

  const char* type = fESD->GetBeamType();
  
  
  //printf("Counter %d\n",fCounter);
  
  fCounter++;
  fNEventsInput->Fill(1);

  //cout << "maxEvent = " << fMaxEvent << endl;
  //if(fCounter%100==0) cout << "fCounter = " << fCounter << endl;
  if((fMaxEvent != 0) && (fMaxEvent < fCounter)) {
    PostData(1, fListHist);
    return;
  }
  //if(fCounter%100==0) cout << "fCounter1 = " << fCounter << endl;
  //cout << "event = " << fCounter << endl;
  
  //printf("Counter %d\n",fCounter);
  
  ///////////////////
  // Check trigger
  ///////////////////
  Bool_t pass = kTRUE;

  if (strstr(type,"p-p")) {
   
    //printf("Will check the triggers\n");

    Int_t numberOfTriggerSelected = fSelectedTrigger->GetEntriesFast();
    //printf("numberofTriggerSelected %d\n",numberOfTriggerSelected);
    if(fRejected) {
      pass = kTRUE;
      for(Int_t k = 0; k < numberOfTriggerSelected; k++){
	const TObjString *const obString=(TObjString*)fSelectedTrigger->At(k);
	const TString tString=obString->GetString();
	if(fESD->IsTriggerClassFired((const char*)tString)) {
	  pass = kFALSE;
	}
      }
    }
    else {
      pass = kFALSE;
      for(Int_t k = 0; k < numberOfTriggerSelected; k++){
	const TObjString *const obString=(TObjString*)fSelectedTrigger->At(k);
	const TString tString=obString->GetString();
	if(fESD->IsTriggerClassFired((const char*)tString)) {
	  pass = kTRUE;
	}
      }
    }
    if(!pass) {
      PostData(1, fListHist);
      return;
    }   

  }
    
  //printf("Class Fired %s\n",(const char*)fESD->GetFiredTriggerClasses());
  //printf("Trigger passed\n");
  
  ///////////////////////////////
  // Require a primary vertex
  //////////////////////////////
  if(fRequirePrimaryVertex) {
    const AliESDVertex* vtxESD = 0x0;
    if      (fVtxTPC) vtxESD = fESD->GetPrimaryVertexTPC() ;
    else if (fVtxSPD) vtxESD = fESD->GetPrimaryVertexSPD() ;
    else              vtxESD = fESD->GetPrimaryVertexTracks() ;
    if(!vtxESD){
      PostData(1, fListHist);
      return;
    }
    Int_t nCtrb = vtxESD->GetNContributors();
    if(nCtrb < fMinNbContributors) {
      PostData(1, fListHist);     
      return;
    }
    Double_t zPosition = vtxESD->GetZ();
    if(TMath::Abs(zPosition) > fRangePrimaryVertexZ) {
      PostData(1, fListHist);
      return;
    }     
    
  }
  
  //printf("Primary vertex passed\n");
  
  //////////////////////////////////////
  // Requirement on number of good tracks
  //////////////////////////////////////
  Int_t nGoodParticles = 0;
  Double_t nbTracks = fESD->GetNumberOfTracks();
  for(Int_t itrack = 0; itrack < nbTracks; itrack++) {
    if(ParticleGood(itrack)) nGoodParticles++;  
  }
  if(fDebug > 0)  {
    // Centrality
    AliCentrality *esdCentrality = fESD->GetCentrality();
    Float_t centrality = esdCentrality->GetCentralityPercentile("V0M");
    //Float_t centralityb = esdCentrality->GetCentralityPercentile("CL1");
    fNbGoodTracks->Fill(nGoodParticles,centrality);
    //printf("centrality %f, centralityb %f\n",centrality,centralityb);
  }
  
  if (strstr(type,"Pb-Pb")) {
    //printf("Will check the number of good tracks\n");
    if((nGoodParticles < fMinNbTracks) || (nGoodParticles > fMaxNbTracks)) {
      PostData(1, fListHist);
      return;
    }
  }
  
  fNEvents->Fill(1);
  
  // In total
  Int_t nbTrdTracks = 0;
  // standalone
  Int_t nbTrdTracksStandalone = 0;
  // offline
  Int_t nbTrdTracksOffline = 0;
  // TPC
  Int_t nbtrackTPC = 0;
  

  
  if (nbTracks <= 0.0) {
    
    if(fDebug > 1) {
      fNbTRDTrack->Fill(nbTrdTracks);
      fNbTRDTrackStandalone->Fill(nbTrdTracksStandalone);
      fNbTRDTrackOffline->Fill(nbTrdTracksOffline);
    }
    PostData(1, fListHist);
    return;
  }
  
  
  fESDfriend = dynamic_cast<AliESDfriend*> (fESD->FindListObject("AliESDfriend"));
  if(!fESDfriend){
    AliError("fESDfriend not available");
    PostData(1, fListHist);
    return;
  }

  if(fESDfriend->TestSkipBit()) {
    PostData(1, fListHist);
    return;
  }
  
  //printf("has friends\n");

  /*
  ////////////////////////////////////
   // Check the number of TPC tracks
   ///////////////////////////////////
   //printf("Nb of tracks %f\n",nbTracks);
   for(Int_t itrk = 0; itrk < nbTracks; itrk++){
     // Get ESD track
     fkEsdTrack = fESD->GetTrack(itrk);
     ULong_t status = fkEsdTrack->GetStatus(); 
     if(status&(AliESDtrack::kTPCout)) nbtrackTPC++;
     if((status&(AliESDtrack::kTRDout)) && (!(status&(AliESDtrack::kTRDin)))) {
       nbTrdTracks++;    
       nbTrdTracksStandalone++;
     }
     if((status&(AliESDtrack::kTRDin))) {
       nbTrdTracks++;    
       nbTrdTracksOffline++;
     }
   }
   
   if((nbtrackTPC>0) && (nbTrdTracks > (3.0*nbtrackTPC))) pass = kFALSE;
   
   if(fDebug > 1) {
     
     fNbTRDTrack->Fill(nbTrdTracks);
     fNbTRDTrackStandalone->Fill(nbTrdTracksStandalone);
     fNbTRDTrackOffline->Fill(nbTrdTracksOffline);
     fNbTPCTRDtrack->Fill(nbTrdTracks,nbtrackTPC);
   
   }

   if(!pass) {
     PostData(1, fListHist);
     return;
   }
  */
  
  /////////////////////////////////////
  // Loop on AliESDtrack
  ////////////////////////////////////
  //printf("Nb of tracks %f\n",nbTracks);      
  for(int itrk=0; itrk < nbTracks; ++itrk){
    
    // Get ESD track
    fkEsdTrack = fESD->GetTrack(itrk);
    if(!fkEsdTrack) continue;
    ULong_t status = fkEsdTrack->GetStatus(); 
    if(status&(AliESDtrack::kTPCout)) ++nbtrackTPC;
    
    // Quality cuts on the AliESDtrack
    if((fEsdTrackCuts) && (!fEsdTrackCuts->IsSelected((AliVParticle *)fkEsdTrack))) {
      //printf("Not a good track\n");
      continue;
    }
    
    // First Absolute gain calibration
    Int_t trdNTracklets = (Int_t) fkEsdTrack->GetTRDntracklets();
    Int_t trdNTrackletsPID = (Int_t) fkEsdTrack->GetTRDntrackletsPID(); 
    if((trdNTracklets > 0) && (trdNTrackletsPID > 0)) {
      for(Int_t iPlane = 0; iPlane < 6; ++iPlane){
	//Double_t slide = fkEsdTrack->GetTRDslice(iPlane);
	//printf("Number of slide %d\n",fkEsdTrack->GetNumberOfTRDslices());
	//Double_t momentum = fkEsdTrack->GetTRDmomentum(iPlane);
	//printf("momentum %f, slide %f\n",momentum,slide);
	if(fkEsdTrack->GetTRDslice(iPlane) > 0.0) 
	  fAbsoluteGain->Fill(fkEsdTrack->GetTRDslice(iPlane)*8.0/100.0,
			      fkEsdTrack->GetTRDmomentum(iPlane)); 
      }
    }     
    
    // Other cuts
    Bool_t good = kTRUE;
    Bool_t standalonetrack = kFALSE;
    Bool_t offlinetrack = kFALSE;
    //ULong_t status = fkEsdTrack->GetStatus();
    
    fFriendTrack = fESDfriend->GetTrack(itrk);
    if(!fFriendTrack)  {
      //printf("No friend track %d\n",itrk);
      continue;
    }
    //////////////////////////////////////
    // Loop on calibration objects
    //////////////////////////////////////
    Int_t icalib=0;
    Int_t nTRDtrackV1=0;
    while((fCalibObject = (TObject *)(fFriendTrack->GetCalibObject(icalib++)))){
      //printf("Name %s\n",fCalibObject->IsA()->GetName());
      if(strcmp(fCalibObject->IsA()->GetName(), "AliTRDtrackV1") != 0) continue;
      //printf("Find the calibration object\n");
      ++nTRDtrackV1;
      
      if((status&(AliESDtrack::kTRDout)) && (!(status&(AliESDtrack::kTRDin)))) {
	standalonetrack = kTRUE;
      }
      if((status&(AliESDtrack::kTRDin))) {
	offlinetrack = kTRUE;
      }
      if(fOfflineTracks){
	if(!offlinetrack){
	  good = kFALSE;
	}
      }
      else if(fStandaloneTracks){
	if(!standalonetrack){
	  good = kFALSE;
	}
      }
      
      fTrdTrack = (AliTRDtrackV1 *)fCalibObject;
      if(good) {
	//cout << "good" << endl;
	fTRDCalibraFillHisto->UpdateHistogramsV1(fTrdTrack);
	//printf("Fill fTRDCalibraFillHisto\n");
      }
      
      //////////////////////////////////
      // Debug 
      ////////////////////////////////
      
      if(fDebug > 0) {
	
	//printf("Enter debug\n");
	
	Int_t nbtracklets = 0;
	
	// Check some stuff
	Bool_t standalonetracklet = kFALSE;  
	const AliTRDseedV1 *tracklet = 0x0;
	//////////////////////////////////////
	// Loop tracklets
	///////////////////////////////////// 
	Int_t nbclusters=0;
	Double_t phtb[AliTRDseedV1::kNtb];
        memset(phtb, 0, AliTRDseedV1::kNtb*sizeof(Double_t));
	Double_t sum = 0.0;
	Float_t normalisation = 6.67;
	Int_t detector = 0;
	Int_t sector = 0;
	for(Int_t itr = 0; itr < 6; ++itr){
	  
	  if(!(tracklet = fTrdTrack->GetTracklet(itr))) continue;
	  if(!tracklet->IsOK()) continue;
	  ++nbtracklets;
	  standalonetracklet = kFALSE; 
	  if(tracklet->IsStandAlone()) standalonetracklet = kTRUE;
	  
	  nbclusters = 0;
	  memset(phtb, 0, AliTRDseedV1::kNtb*sizeof(Double_t));
	  sum = 0.0;
	  detector = 0;
	  sector = 0;
	  //Int_t crossrow = 0;
	  
	  // Check no shared clusters
	  //for(int icc=AliTRDseedV1::kNtb; icc<AliTRDseedV1::kNclusters; icc++){
	  //  if((fcl = tracklet->GetClusters(icc)))  crossrow = 1;
	  // }
	  
	  // Loop on clusters
	  Int_t time = 0;
	  Float_t ch = 0;
	  Float_t qcl = 0;
	  for(int ic=0; ic<AliTRDseedV1::kNtb; ++ic){
	    
	    if(!(fCl = tracklet->GetClusters(ic))) continue;
	    ++nbclusters;
	    time = fCl->GetPadTime();
	    ch =  tracklet->GetdQdl(ic);
	    qcl = TMath::Abs(fCl->GetQ());
	    detector = fCl->GetDetector();	  
	    // Add the charge if shared cluster
	    if((ic+AliTRDseedV1::kNtb) < AliTRDseedV1::kNclusters) {
	      if((fCl = tracklet->GetClusters(ic+AliTRDseedV1::kNtb))) {
		qcl += TMath::Abs(fCl->GetQ());
		//printf("Add the cluster charge\n");
	      }
	    }
	    if((time>-1) && (time<fNbTimeBins)) phtb[time]=qcl;
	    if((fCalDetGain) && (fCalDetGain->GetValue(detector) > 0.0)) sum += ch*fCalDetGain->GetValue(detector)/normalisation;	
	    else sum += ch/normalisation;
	    
	    if(fDebug > 1) {
	      fNbTimeBin->Fill(time);
	      if(tracklet->IsStandAlone()) fNbTimeBinStandalone->Fill(time);
	      else fNbTimeBinOffline->Fill(time);
	    }
	  }
	  sector = AliTRDgeometry::GetSector(detector);
	  
	  if(fDebug > 1) {
	    fNbTracklets->Fill(detector);
	    if(tracklet->IsStandAlone()) fNbTrackletsStandalone->Fill(detector);
	    else fNbTrackletsOffline->Fill(detector);
	    
	    fNbClusters->Fill(nbclusters);
	    if(tracklet->IsStandAlone())  fNbClustersStandalone->Fill(nbclusters);
	    else  fNbClustersOffline->Fill(nbclusters);
	  }	   
	  
	  if(fDebug > 0) {
	    if((nbclusters > fLow) && (nbclusters < fHigh)){
	      if(fRelativeScale > 0.0) sum = sum/fRelativeScale;	       
	      fCH2dSM->Fill(sum,sector+0.5);
	      fCH2dSum->Fill(sum,0.5);
	      Bool_t checknoise = kTRUE;
	      if(fMaxCluster > 0) {
		if(phtb[0] > fMaxCluster) checknoise = kFALSE;
		if(fNbTimeBins > fNbMaxCluster) {
		  for(Int_t k = (fNbTimeBins-fNbMaxCluster); k < fNbTimeBins; k++){
		    if(phtb[k] > fMaxCluster) checknoise = kFALSE;
		  }
		}
	      }
	      if(checknoise) {	       
		for(int ic=0; ic<fNbTimeBins; ic++){
		  if(fFillZero) {
		    fPH2dSum->Fill((Double_t)(ic/10.0),0.5,(Double_t)phtb[ic]);
		    fPH2dSM->Fill((Double_t)(ic/10.0),sector+0.5,(Double_t)phtb[ic]);
		  }
		  else {
		    if(phtb[ic] > 0.0) {
		      fPH2dSum->Fill((Double_t)(ic/10.0),0.0,(Double_t)phtb[ic]);
		      fPH2dSM->Fill((Double_t)(ic/10.0),sector+0.5,(Double_t)phtb[ic]);
		    }
		  }
		}
	      }
	    }
	  }
	} // loop on tracklets
	
      } // debug
      
    }// while calibration objects
    if(nTRDtrackV1 > 0) {
      ++nbTrdTracks;      
      if((status&(AliESDtrack::kTRDout)) && (!(status&(AliESDtrack::kTRDin)))) {
	++nbTrdTracksStandalone;
      }
      if((status&(AliESDtrack::kTRDin))) {
	++nbTrdTracksOffline;
      }
    }
    //delete fFriendTrack;
  } // loop ESD track
  
  if(fDebug > 1) {
    fNbTRDTrack->Fill(nbTrdTracks);
    fNbTRDTrackStandalone->Fill(nbTrdTracksStandalone);
    fNbTRDTrackOffline->Fill(nbTrdTracksOffline);
    fNbTPCTRDtrack->Fill(nbTrdTracks,nbtrackTPC);
  }
  
  // Post output data
  PostData(1, fListHist);
  //cout << "AliTRDCalibTask::Exec() OUT" << endl;
}
     
//________________________________________________________________________
void AliTRDCalibTask::Terminate(Option_t *) 
{
  //
  // Terminate
  //
  
  if(fTRDCalibraFillHisto) fTRDCalibraFillHisto->DestroyDebugStreamer();

 
}
//_______________________________________________________
Bool_t AliTRDCalibTask::Load(const Char_t *filename)
{
  //
  // Generic container loader
  //

  if(!TFile::Open(filename)){
    //AliWarning(Form("Couldn't open file %s.", filename));
    return kFALSE;
  }
  TList *o = 0x0;
  if(!(o = (TList*)gFile->Get(GetName()))){
    //AliWarning("Missing histogram container.");
    return kFALSE;
  }
  fListHist = (TList*)o->Clone(GetName());
  gFile->Close();
  return kTRUE;
}
//_______________________________________________________
Bool_t AliTRDCalibTask::Load(TList *lister)
{
  //
  // Generic container loader
  //

  fListHist = (TList*)lister->Clone(GetName());
  return kTRUE;
}
//________________________________________________________________________
void AliTRDCalibTask::Plot() 
{
  //
  // Plot the histos stored in the TList
  //
 
  if(!fListHist) return;

  /////////////////////////////////////
  // Take the debug stuff
  /////////////////////////////////////

  TH1I *nEvents  = (TH1I *) fListHist->FindObject("NEvents");

  TH2F *absoluteGain  = (TH2F *) fListHist->FindObject("AbsoluteGain");

  TH1F *trdTrack = (TH1F *) fListHist->FindObject("TRDTrack");
  TH1F *trdTrackOffline = (TH1F *) fListHist->FindObject("TRDTrackOffline");
  TH1F *trdTrackStandalone = (TH1F *) fListHist->FindObject("TRDTrackStandalone");

  TH2F *tpctrdTrack = (TH2F *) fListHist->FindObject("NbTPCTRDtrack");

  TH1F *nbTimeBin = (TH1F *) fListHist->FindObject("NbTimeBin");
  TH1F *nbTimeBinOffline = (TH1F *) fListHist->FindObject("NbTimeBinOffline");
  TH1F *nbTimeBinStandalone = (TH1F *) fListHist->FindObject("NbTimeBinStandalone");

  TH1F *nbClusters = (TH1F *) fListHist->FindObject("NbClusters");
  TH1F *nbClustersOffline = (TH1F *) fListHist->FindObject("NbClustersOffline");
  TH1F *nbClustersStandalone = (TH1F *) fListHist->FindObject("NbClustersStandalone");

  TH1F *nbTracklets = (TH1F *) fListHist->FindObject("NbTracklets");
  TH1F *nbTrackletsOffline = (TH1F *) fListHist->FindObject("NbTrackletsOffline");
  TH1F *nbTrackletsStandalone = (TH1F *) fListHist->FindObject("NbTrackletsStandalone");
  
  /////////////////////////////////////
  // Take the calibration objects
  /////////////////////////////////////

  TH2I *ch2d = (TH2I *) fListHist->FindObject("CH2d");
  TProfile2D *ph2d = (TProfile2D *) fListHist->FindObject("PH2d");

  TH2I *ch2dSum = (TH2I *) fListHist->FindObject("CH2dSum");
  TProfile2D *ph2dSum = (TProfile2D *) fListHist->FindObject("PH2dSum");

  TH2I *ch2dSM = (TH2I *) fListHist->FindObject("CH2dSM");
  TProfile2D *ph2dSM = (TProfile2D *) fListHist->FindObject("PH2dSM");
  
  AliTRDCalibraVdriftLinearFit *linearfit = (AliTRDCalibraVdriftLinearFit *) fListHist->FindObject("AliTRDCalibraVdriftLinearFit");
  
  ////////////////////////////////////////////////
  // Add the AliTRDCalibraVdriftLinearFit
  ///////////////////////////////////////////////
  
  Int_t first = 0;
  TH2S *histolinearfitsum = 0x0;
  
  if(linearfit) {
    for(Int_t det = 0; det < 540; det++) {
      if(linearfit->GetLinearFitterHisto(det)){
	if(TMath::Abs(first)<0.0001){
	  histolinearfitsum = linearfit->GetLinearFitterHisto(det);
	  first += 1;
	}
	else {
          if (histolinearfitsum) {
	    histolinearfitsum->Add(linearfit->GetLinearFitterHisto(det));
	  }
	}
      }
    }
  }

  ///////////////////////////
  // Style
  //////////////////////////

  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.13);

  /////////////////////////
  // Plot
  ////////////////////////

 if(nEvents) {
   
    TCanvas *debugEvents = new TCanvas("cNEvents","cNEvents",10,10,510,510);
    debugEvents->cd(1);
    if(nEvents) nEvents->Draw();
      
  }

 if(absoluteGain) {
   
    TCanvas *debugAbsoluteGain = new TCanvas("cAbsoluteGain","cAbsoluteGain",10,10,510,510);
    debugAbsoluteGain->cd(1);
    if(absoluteGain) absoluteGain->Draw();
      
  }

  if(trdTrack || tpctrdTrack) {
    
    TCanvas *debugtrdtpcTrack = new TCanvas("TRDtracktpctrdtrack","TRDtracktpctrdtrack",10,10,510,510);
    debugtrdtpcTrack->Divide(2,1);
    debugtrdtpcTrack->cd(1);
    if(trdTrack) trdTrack->Draw();
    if(trdTrackOffline) trdTrackOffline->Draw("same");
    if(trdTrackStandalone) trdTrackStandalone->Draw("same");
    TLegend *leg = new TLegend(0.4,0.6,0.89,0.89);
    if(trdTrack) leg->AddEntry(trdTrack,"All","p");
    if(trdTrackOffline) leg->AddEntry(trdTrackOffline,"Offline","p");
    if(trdTrackStandalone) leg->AddEntry(trdTrackStandalone,"Standalone","p");
    leg->Draw("same");
    debugtrdtpcTrack->cd(2);
    if(tpctrdTrack) tpctrdTrack->Draw();
    TLine *line = new TLine(0.0,0.0,100.0,100.0);
    line->Draw("same");
    
  }
 
  if(nbTimeBin || nbTracklets || nbClusters) {
    
    TCanvas *debugTracklets = new TCanvas("TRDtimebintrackletcluster","TRDtimebintrackletcluster",10,10,510,510);
    debugTracklets->Divide(3,1);
    debugTracklets->cd(1);
    if(nbTimeBin) nbTimeBin->Draw();
    if(nbTimeBinOffline) nbTimeBinOffline->Draw("same");
    if(nbTimeBinStandalone) nbTimeBinStandalone->Draw("same");
    TLegend *lega = new TLegend(0.4,0.6,0.89,0.89);
    if(nbTimeBin) lega->AddEntry(nbTimeBin,"All","p");
    if(nbTimeBinOffline) lega->AddEntry(nbTimeBinOffline,"Offline","p");
    if(nbTimeBinStandalone) lega->AddEntry(nbTimeBinStandalone,"Standalone","p");
    lega->Draw("same");
    debugTracklets->cd(2);
    if(nbTracklets) nbTracklets->Draw();
    if(nbTrackletsOffline) nbTrackletsOffline->Draw("same");
    if(nbTrackletsStandalone) nbTrackletsStandalone->Draw("same");
    TLegend *legb = new TLegend(0.4,0.6,0.89,0.89);
    if(nbTracklets) legb->AddEntry(nbTracklets,"All","p");
    if(nbTrackletsOffline) legb->AddEntry(nbTrackletsOffline,"Offline","p");
    if(nbTrackletsStandalone) legb->AddEntry(nbTrackletsStandalone,"Standalone","p");
    legb->Draw("same");
    debugTracklets->cd(3);
    if(nbClusters) nbClusters->Draw();
    if(nbClustersOffline) nbClustersOffline->Draw("same");
    if(nbClustersStandalone) nbClustersStandalone->Draw("same");
    TLegend *legc = new TLegend(0.4,0.6,0.89,0.89);
    if(nbClusters) legc->AddEntry(nbClusters,"All","p");
    if(nbClustersOffline) legc->AddEntry(nbClustersOffline,"Offline","p");
    if(nbClustersStandalone) legc->AddEntry(nbClustersStandalone,"Standalone","p");
    legc->Draw("same");
  
  }

  if(ch2dSum || ph2dSum || histolinearfitsum) {
    
    TCanvas *debugSum = new TCanvas("SumCalibrationObjects","SumCalibrationObjects",10,10,510,510);
    debugSum->Divide(3,1);
    debugSum->cd(1);
    if(ch2dSum) ch2dSum->Draw("lego");
    debugSum->cd(2);
    if(ph2dSum) ph2dSum->Draw("lego");
    debugSum->cd(3);
    if(histolinearfitsum) histolinearfitsum->Draw();
  
  }

  if(ch2dSM || ph2dSM) {
    
    TCanvas *debugSM = new TCanvas("SMCalibrationObjects","SMCalibrationObjects",10,10,510,510);
    debugSM->Divide(2,1);
    debugSM->cd(1);
    if(ch2dSM) ch2dSM->Draw("lego");
    debugSM->cd(2);
    if(ph2dSM) ph2dSM->Draw("lego");
     
  }

  if(ch2d || ph2d) {
    
    TCanvas *debug = new TCanvas("CalibrationObjects","CalibrationObjects",10,10,510,510);
    debug->Divide(2,1);
    debug->cd(1);
    if(ch2d) ch2d->Draw("lego");
    debug->cd(2);
    if(ph2d) ph2d->Draw("lego");
     
  }
 
}
//_______________________________________________________________________________________
void  AliTRDCalibTask::AddTask(const AliTRDCalibTask * calibTask) {

  //
  // Add stats
  //

  TList *listcalibTask = calibTask->GetList();
  if(!listcalibTask) return;

  TH1I *nEvents  = (TH1I *) listcalibTask->FindObject("NEvents");
  TH1I *nEventsInput  = (TH1I *) listcalibTask->FindObject("NEventsInput");
  TH2F *absoluteGain  = (TH2F *) listcalibTask->FindObject("AbsoluteGain");

  TH1F *trdTrack = (TH1F *) listcalibTask->FindObject("TRDTrack");
  TH1F *trdTrackOffline = (TH1F *) listcalibTask->FindObject("TRDTrackOffline");
  TH1F *trdTrackStandalone = (TH1F *) listcalibTask->FindObject("TRDTrackStandalone");

  TH2F *tpctrdTrack = (TH2F *) listcalibTask->FindObject("NbTPCTRDtrack");

  TH1F *nbTimeBin = (TH1F *) listcalibTask->FindObject("NbTimeBin");
  TH1F *nbTimeBinOffline = (TH1F *) listcalibTask->FindObject("NbTimeBinOffline");
  TH1F *nbTimeBinStandalone = (TH1F *) listcalibTask->FindObject("NbTimeBinStandalone");

  TH1F *nbClusters = (TH1F *) listcalibTask->FindObject("NbClusters");
  TH1F *nbClustersOffline = (TH1F *) listcalibTask->FindObject("NbClustersOffline");
  TH1F *nbClustersStandalone = (TH1F *) listcalibTask->FindObject("NbClustersStandalone");

  TH1F *nbTracklets = (TH1F *) listcalibTask->FindObject("NbTracklets");
  TH1F *nbTrackletsOffline = (TH1F *) listcalibTask->FindObject("NbTrackletsOffline");
  TH1F *nbTrackletsStandalone = (TH1F *) listcalibTask->FindObject("NbTrackletsStandalone");
  
  TH2I *ch2d = (TH2I *) listcalibTask->FindObject("CH2d");
  TProfile2D *ph2d = (TProfile2D *) listcalibTask->FindObject("PH2d");
  TProfile2D *prf2d = (TProfile2D *) listcalibTask->FindObject("PRF2d");

  TH2I *ch2dSum = (TH2I *) listcalibTask->FindObject("CH2dSum");
  TProfile2D *ph2dSum = (TProfile2D *) listcalibTask->FindObject("PH2dSum");

  TH2I *ch2dSM = (TH2I *) listcalibTask->FindObject("CH2dSM");
  TProfile2D *ph2dSM = (TProfile2D *) listcalibTask->FindObject("PH2dSM");
  
  AliTRDCalibraVdriftLinearFit *linearfit = (AliTRDCalibraVdriftLinearFit *) listcalibTask->FindObject("AliTRDCalibraVdriftLinearFit");  
  AliTRDCalibraVector *calibraVector = (AliTRDCalibraVector *) listcalibTask->FindObject("AliTRDCalibraVector");  

  //

  TH1I *inEventsInput  = (TH1I *) fListHist->FindObject("NEventsInput");
  TH1I *inEvents  = (TH1I *) fListHist->FindObject("NEvents");
  TH2F *iabsoluteGain  = (TH2F *) fListHist->FindObject("AbsoluteGain");

  TH1F *itrdTrack = (TH1F *) fListHist->FindObject("TRDTrack");
  TH1F *itrdTrackOffline = (TH1F *) fListHist->FindObject("TRDTrackOffline");
  TH1F *itrdTrackStandalone = (TH1F *) fListHist->FindObject("TRDTrackStandalone");

  TH2F *itpctrdTrack = (TH2F *) fListHist->FindObject("NbTPCTRDtrack");

  TH1F *inbTimeBin = (TH1F *) fListHist->FindObject("NbTimeBin");
  TH1F *inbTimeBinOffline = (TH1F *) fListHist->FindObject("NbTimeBinOffline");
  TH1F *inbTimeBinStandalone = (TH1F *) fListHist->FindObject("NbTimeBinStandalone");

  TH1F *inbClusters = (TH1F *) fListHist->FindObject("NbClusters");
  TH1F *inbClustersOffline = (TH1F *) fListHist->FindObject("NbClustersOffline");
  TH1F *inbClustersStandalone = (TH1F *) fListHist->FindObject("NbClustersStandalone");

  TH1F *inbTracklets = (TH1F *) fListHist->FindObject("NbTracklets");
  TH1F *inbTrackletsOffline = (TH1F *) fListHist->FindObject("NbTrackletsOffline");
  TH1F *inbTrackletsStandalone = (TH1F *) fListHist->FindObject("NbTrackletsStandalone");
  
  TH2I *ich2d = (TH2I *) fListHist->FindObject("CH2d");
  TProfile2D *iph2d = (TProfile2D *) fListHist->FindObject("PH2d");
  TProfile2D *iprf2d = (TProfile2D *) fListHist->FindObject("PRF2d");

  TH2I *ich2dSum = (TH2I *) fListHist->FindObject("CH2dSum");
  TProfile2D *iph2dSum = (TProfile2D *) fListHist->FindObject("PH2dSum");

  TH2I *ich2dSM = (TH2I *) fListHist->FindObject("CH2dSM");
  TProfile2D *iph2dSM = (TProfile2D *) fListHist->FindObject("PH2dSM");
  
  AliTRDCalibraVdriftLinearFit *ilinearfit = (AliTRDCalibraVdriftLinearFit *) fListHist->FindObject("AliTRDCalibraVdriftLinearFit");  
  AliTRDCalibraVector *icalibraVector = (AliTRDCalibraVector *) fListHist->FindObject("AliTRDCalibraVector");  


  // Add

  if(nEventsInput) {
    if(inEventsInput) {
      inEventsInput->Add(nEventsInput);
      //printf("Add Events\n");
    }
    else {
      //printf("Create new Events\n");
      inEventsInput = new TH1I(*nEventsInput);
      fListHist->Add(inEventsInput);
    }
  }
  
  if(nEvents) {
    if(inEvents) {
      inEvents->Add(nEvents);
      //printf("Add Events\n");
    }
    else {
      //printf("Create new Events\n");
      inEvents = new TH1I(*nEvents);
      fListHist->Add(inEvents);
    }
  }
  
  if(absoluteGain) {
    if(iabsoluteGain) iabsoluteGain->Add(absoluteGain);
    else {
      iabsoluteGain = new TH2F(*absoluteGain);
      fListHist->Add(iabsoluteGain);
    }
  }
  
  if(trdTrack) {
    if(itrdTrack) itrdTrack->Add(trdTrack);
    else {
     itrdTrack = new TH1F(*trdTrack);
     fListHist->Add(itrdTrack);
    }
  }

  if(trdTrackOffline) {
    if(itrdTrackOffline) itrdTrackOffline->Add(trdTrackOffline);
    else {
      itrdTrackOffline = new TH1F(*trdTrackOffline);
      fListHist->Add(itrdTrackOffline);
    }
  }

  if(trdTrackStandalone) {
    if(itrdTrackStandalone) itrdTrackStandalone->Add(trdTrackStandalone);
    else {
      itrdTrackStandalone = new TH1F(*trdTrackStandalone);
      fListHist->Add(itrdTrackStandalone);
    }
  }

  if(tpctrdTrack) {
    if(itpctrdTrack) itpctrdTrack->Add(tpctrdTrack);
    else {
      itpctrdTrack = new TH2F(*tpctrdTrack);
      fListHist->Add(itpctrdTrack);
    }
  }

  if(nbTimeBin) {
    if(inbTimeBin) inbTimeBin->Add(nbTimeBin);
    else {
      inbTimeBin = new TH1F(*inbTimeBin);
      fListHist->Add(inbTimeBin);
    }
  }

  if(nbTimeBinOffline) {
    if(inbTimeBinOffline) inbTimeBinOffline->Add(nbTimeBinOffline);
    else {
      inbTimeBinOffline = new TH1F(*nbTimeBinOffline);
      fListHist->Add(inbTimeBinOffline);
    }
  }
  
  if(nbTimeBinStandalone) {
    if(inbTimeBinStandalone) inbTimeBinStandalone->Add(nbTimeBinStandalone);
    else {
      inbTimeBinStandalone = new TH1F(*nbTimeBinStandalone);
      fListHist->Add(inbTimeBinStandalone);
    }
  }

  if(nbClusters) {
    if(inbClusters) inbClusters->Add(nbClusters);
    else {
      inbClusters = new TH1F(*nbClusters);
      fListHist->Add(inbClusters);
    }
  }
  
  if(nbClustersOffline) {
    if(inbClustersOffline) inbClustersOffline->Add(nbClustersOffline);
    else {
      inbClustersOffline = new TH1F(*nbClustersOffline);
      fListHist->Add(inbClustersOffline);
    }
  }
  
  if(nbClustersStandalone) {
    if(inbClustersStandalone) inbClustersStandalone->Add(nbClustersStandalone);
    else {
      inbClustersStandalone = new TH1F(*nbClustersStandalone);
      fListHist->Add(inbClustersStandalone);
    }
  }

  if(nbTracklets) {
    if(inbTracklets) inbTracklets->Add(nbTracklets);
    else {
      inbTracklets = new TH1F(*nbTracklets);
      fListHist->Add(inbTracklets);
    }
  }

  if(nbTrackletsOffline) {
    if(inbTrackletsOffline) inbTrackletsOffline->Add(nbTrackletsOffline);
    else {
      inbTrackletsOffline = new TH1F(*nbTrackletsOffline);
      fListHist->Add(inbTrackletsOffline);
    }
  }
  
  if(nbTrackletsStandalone) {
    if(inbTrackletsStandalone) inbTrackletsStandalone->Add(nbTrackletsStandalone);
    else {
      inbTrackletsStandalone = new TH1F(*nbTrackletsStandalone);
      fListHist->Add(inbTrackletsStandalone);
    }
  }
  
  if(ch2d) {
    if(ich2d) ich2d->Add(ch2d);
    else {
      ich2d = new TH2I(*ch2d);
      fListHist->Add(ich2d);
    }
  }

  if(ph2d) {
    if(iph2d) iph2d->Add(ph2d);
    else {
      iph2d = new TProfile2D(*ph2d);
      fListHist->Add(iph2d);
    }
  }

  if(prf2d) {
    if(iprf2d) iprf2d->Add(prf2d);
    else {
      iprf2d = new TProfile2D(*prf2d);
      fListHist->Add(iprf2d);
    }
  }

  if(ch2dSum) {
    if(ich2dSum) ich2dSum->Add(ch2dSum);
    else {
      ich2dSum = new TH2I(*ch2dSum);
      fListHist->Add(ich2dSum);
    }
  }

  if(ph2dSum) {
    if(iph2dSum) iph2dSum->Add(ph2dSum);
    else {
      iph2dSum = new TProfile2D(*ph2dSum);
      fListHist->Add(iph2dSum);
    }
  }

  if(ch2dSM) {
    if(ich2dSM) ich2dSM->Add(ch2dSM);
    else {
      ich2dSM = new TH2I(*ch2dSM);
      fListHist->Add(ich2dSM);
    }
  }

  if(ph2dSM) {
    if(iph2dSM) iph2dSM->Add(ph2dSM);
    else {
      iph2dSM = new TProfile2D(*ph2dSM);
      fListHist->Add(iph2dSM);
    }
  }
  
  if(linearfit) {
    if(ilinearfit) ilinearfit->Add(linearfit);
    else {
      ilinearfit = new AliTRDCalibraVdriftLinearFit(*linearfit);
      fListHist->Add(ilinearfit);
    }
  }

  if(calibraVector) {
    if(icalibraVector) icalibraVector->Add(calibraVector);
    else {
      icalibraVector = new AliTRDCalibraVector(*calibraVector);
      fListHist->Add(icalibraVector);
    }
  }
  
}
//________________________________________________________________________________
Long64_t AliTRDCalibTask::Merge(TCollection *li) {
  
  //
  // merge component
  //
  
  TIterator* iter = li->MakeIterator();
  AliTRDCalibTask* cal = 0;

  while ((cal = (AliTRDCalibTask*)iter->Next())) {
    if (!cal->InheritsFrom(AliTRDCalibTask::Class())) {
      Error("Merge","Attempt to add object of class %s to a %s", cal->ClassName(), this->ClassName());
      return -1;
    }

    // add histograms here...
    this->AddTask(cal);
    
  }
  
  return 0;
  
}
//_____________________________________________________
Bool_t AliTRDCalibTask::SetVersionSubversion(){
  //
  // Load Chamber Gain factors into the Tender supply
  //
  
  printf("SetVersionSubversion\n");

  //find previous entry from the UserInfo
  TTree *tree=((TChain*)GetInputData(0))->GetTree();
  if (!tree) {
    AliError("Tree not found in ESDhandler");
    return kFALSE;
  }
 	 
  TList *userInfo=(TList*)tree->GetUserInfo();
  if (!userInfo) {
    AliError("No UserInfo found in tree");
    return kFALSE;
  }

  TList *cdbList=(TList*)userInfo->FindObject("cdbList");
  if (!cdbList) {
    AliError("No cdbList found in UserInfo");
    if (AliLog::GetGlobalLogLevel()>=AliLog::kError) userInfo->Print();
    return kFALSE;
  }
 	
  TIter nextCDB(cdbList);
  TObjString *os=0x0;
  while ( (os=(TObjString*)nextCDB()) ){
    if(os->GetString().Contains("TRD/Calib/ChamberGainFactor")){
      // Get Old gain calibration
      AliCDBId *id=AliCDBId::MakeFromString(os->GetString());
      fFirstRunGain = id->GetFirstRun();
      fVersionGainUsed = id->GetVersion();
      fSubVersionGainUsed = id->GetSubVersion();
    } else if(os->GetString().Contains("TRD/Calib/ChamberVdrift")){
      // Get Old drift velocity calibration
      AliCDBId *id=AliCDBId::MakeFromString(os->GetString());
      fFirstRunVdrift = id->GetFirstRun();
      fVersionVdriftUsed = id->GetVersion();
      fSubVersionVdriftUsed = id->GetSubVersion();
    } else if(os->GetString().Contains("TRD/Calib/LocalGainFactor")){
      // Get Old drift velocity calibration
      AliCDBId *id=AliCDBId::MakeFromString(os->GetString());
      fFirstRunGainLocal = id->GetFirstRun();
      fVersionGainLocalUsed = id->GetVersion();
      fSubVersionGainLocalUsed = id->GetSubVersion();
    }
  }

  //printf("VersionGain %d, SubversionGain %d, VersionLocalGain %d, Subversionlocalgain %d, Versionvdrift %d, Subversionvdrift %d\n",fVersionGainUsed,fSubVersionGainUsed,fVersionGainLocalUsed,fSubVersionGainLocalUsed,fVersionVdriftUsed,fSubVersionVdriftUsed);

  // Check
  if((fFirstRunGain < 0)            || 
     (fFirstRunGainLocal < 0)       || 
     (fFirstRunVdrift < 0)          || 
     (fVersionGainUsed < 0)         || 
     (fVersionGainLocalUsed < 0)    || 
     (fSubVersionGainUsed < 0)      || 
     (fSubVersionGainLocalUsed < 0) || 
     (fVersionVdriftUsed < 0)       || 
     (fSubVersionVdriftUsed < 0)) {
    AliError("No recent calibration found");
    return kFALSE;
  }
  else return kTRUE;

}
//_________________________________________________________________________________________________________________________
Bool_t AliTRDCalibTask::ParticleGood(int i) const {

  //
  // Definition of good tracks
  //

  
  AliESDtrack *track = fESD->GetTrack(i);
  if (!track->IsOn(AliESDtrack::kTPCrefit)) return 0;        // TPC refit
  if (track->GetTPCNcls() < 90) return 0;                    // number of TPC clusters
  if (fabs(track->Eta())>0.8) return 0;                         // fiducial pseudorapidity
  Float_t r,z;
  track->GetImpactParametersTPC(r,z);
  if (fabs(z)>2.0) return 0;                          // impact parameter in z
  if (fabs(r)>2.0) return 0;                          // impact parameter in xy
  if (r==0) return 0;
  return 1;   


}



