
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
#include "TH2S.h"
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
#include "TLinearFitter.h"
#include "TVectorD.h"

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
#include "AliTRDCalibraExbAltFit.h" 

#include "AliTRDcalibDB.h"
#include "AliCDBId.h"
#include "AliLog.h"

#include "AliTRDCalibChamberStatus.h"

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
      fTRDChamberStatus(0),
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
      fCH2dTest(0),
      fPH2dTest(0),
      fLinearVdriftTest(0),
      fOnInstance(kTRUE),
      fHisto2d(kTRUE),
      fVector2d(kFALSE),
      fVdriftLinear(kTRUE),
      fExbAlt(kFALSE),
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
      fFirstRunExB(-1),
      fVersionExBUsed(-1), 
      fSubVersionExBUsed(-1),
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
  if(fCH2dTest) delete fCH2dTest;
  if(fPH2dTest) delete fPH2dTest;
  if(fLinearVdriftTest) delete fLinearVdriftTest;
  if(fCalDetGain) delete fCalDetGain;
  
  if(fSelectedTrigger) {
    fSelectedTrigger->Delete();
    delete fSelectedTrigger;
  }
  if(fEsdTrackCuts) {
    delete fEsdTrackCuts;
  }
  
  if(fTRDChamberStatus) delete fTRDChamberStatus;
  
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

  // output list
  fListHist = new TList();
  fListHist->SetOwner();
  
  // init chamber status
  fTRDChamberStatus = new AliTRDCalibChamberStatus();
  fTRDChamberStatus->Init();
  fListHist->Add(fTRDChamberStatus->GetSparseI());
  
  // instance calibration 
  fTRDCalibraFillHisto = AliTRDCalibraFillHisto::Instance();
  if(fOnInstance) {
    fTRDCalibraFillHisto->SetHisto2d(fHisto2d); // choose to use histograms
    fTRDCalibraFillHisto->SetVector2d(fVector2d); // choose to use vectors
    fTRDCalibraFillHisto->SetCH2dOn();  // choose to calibrate the gain
    fTRDCalibraFillHisto->SetPH2dOn();  // choose to calibrate the drift velocity
    fTRDCalibraFillHisto->SetPRF2dOn(); // choose to look at the PRF
    fTRDCalibraFillHisto->SetLinearFitterOn(fVdriftLinear); // Other possibility vdrift VDRIFT
    fTRDCalibraFillHisto->SetLinearFitterDebugOn(fVdriftLinear); // Other possibility vdrift
    fTRDCalibraFillHisto->SetExbAltFitOn(fExbAlt); // Alternative method for exb
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
    
    // For testing only
    if(fDebug > 2) fTRDCalibraFillHisto->SetDebugLevel(1); //debug stuff
    
    if(fHisto2d) {  
      fListHist->Add(fTRDCalibraFillHisto->GetCH2d());
      fListHist->Add(fTRDCalibraFillHisto->GetPH2d()); 
      fListHist->Add(fTRDCalibraFillHisto->GetPRF2d());
    } 
    if(fVdriftLinear) fListHist->Add((TObject *)fTRDCalibraFillHisto->GetVdriftLinearFit());
    if(fVector2d) fListHist->Add((TObject *) fTRDCalibraFillHisto->GetCalibraVector()); //calibra vector
    if(fExbAlt) fListHist->Add((TObject *)fTRDCalibraFillHisto->GetExbAltFit());
  }
  fRelativeScale = fTRDCalibraFillHisto->GetRelativeScale(); // Get the relative scale for the gain
  
  fNEvents = new TH1I(Form("NEvents_%s",(const char*)fName),"NEvents", 2, 0, 2);
  fListHist->Add(fNEvents);
  fNEventsInput = new TH1I(Form("NEventsInput_%s",(const char*)fName),"NEventsInput", 2, 0, 2);
  fListHist->Add(fNEventsInput);
  
  // absolute gain calibration even without AliESDfriend
  Int_t nBinsPt = 25;
  Double_t minPt = 0.001;
  Double_t maxPt = 10.0;
  
  Double_t *binLimLogPt = new Double_t[nBinsPt+1];
  Double_t *binLimPt    = new Double_t[nBinsPt+1];
  for(Int_t i=0; i<=nBinsPt; i++) binLimLogPt[i]=(Double_t)TMath::Log10(minPt) + (TMath::Log10(maxPt)-TMath::Log10(minPt))/nBinsPt*(Double_t)i ;
  for(Int_t i=0; i<=nBinsPt; i++) binLimPt[i]=(Double_t)TMath::Power(10,binLimLogPt[i]);
  
  fAbsoluteGain = new TH2F(Form("AbsoluteGain_%s",(const char*)fName),"AbsoluteGain", 200, 0.0, 700.0, nBinsPt, binLimPt);
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

    fLinearVdriftTest = new TH2S(Form("LFDV0testversion_%s",(const char*)fName),"LFDV0testversion",36,-0.9,0.9,48,-1.2,1.2);
    fLinearVdriftTest->SetXTitle("tan(phi_{track})");
    fLinearVdriftTest->SetYTitle("dy/dt");
    fLinearVdriftTest->SetZTitle("Number of tracklets");
    fLinearVdriftTest->SetStats(0);
    fLinearVdriftTest->SetDirectory(0);
    
    // Standart with AliESDfriend
    fPH2dTest = new TProfile2D(Form("PH2dTest_%s",(const char*)fName),"Nz0Nrphi0"
			    ,fNbTimeBins,-0.05,(Double_t)((fNbTimeBins-0.5)/10.0)
			   ,540,0,540);
    fPH2dTest->SetYTitle("Det/pad groups");
    fPH2dTest->SetXTitle("time [#mus]");
    fPH2dTest->SetZTitle("<PH> [a.u.]");
    fPH2dTest->SetStats(0);
    //
    fCH2dTest = new TH2I(Form("CH2dTest_%s",(const char*)fName),"Nz0Nrphi0",50,0,300,540,0,540);
    fCH2dTest->SetYTitle("Det/pad groups");
    fCH2dTest->SetXTitle("charge deposit [a.u]");
    fCH2dTest->SetZTitle("counts");
    fCH2dTest->SetStats(0);
    fCH2dTest->Sumw2();

    //
    fPH2dSM = new TProfile2D(Form("PH2dSM_%s",(const char*)fName),"Nz10Nrphi10"
			    ,fNbTimeBins,-0.05,(Double_t)((fNbTimeBins-0.5)/10.0)
			   ,18,0,18);
    fPH2dSM->SetYTitle("Det/pad groups");
    fPH2dSM->SetXTitle("time [#mus]");
    fPH2dSM->SetZTitle("<PH> [a.u.]");
    fPH2dSM->SetStats(0);
    //
    fCH2dSM = new TH2I(Form("CH2dSM_%s",(const char*)fName),"Nz10Nrphi10",50,0,300,18,0,18);
    fCH2dSM->SetYTitle("Det/pad groups");
    fCH2dSM->SetXTitle("charge deposit [a.u]");
    fCH2dSM->SetZTitle("counts");
    fCH2dSM->SetStats(0);
    fCH2dSM->Sumw2();
    //
    fPH2dSum = new TProfile2D(Form("PH2dSum_%s",(const char*)fName),"Nz100Nrphi100"
			    ,fNbTimeBins,-0.05,(Double_t)((fNbTimeBins-0.5)/10.0)
			    ,1,0,1);
    fPH2dSum->SetYTitle("Det/pad groups");
    fPH2dSum->SetXTitle("time [#mus]");
    fPH2dSum->SetZTitle("<PH> [a.u.]");
    fPH2dSum->SetStats(0);
    //
    fCH2dSum = new TH2I(Form("CH2dSum_%s",(const char*)fName),"Nz100Nrphi100",50,0,300,1,0,1);
    fCH2dSum->SetYTitle("Det/pad groups");
    fCH2dSum->SetXTitle("charge deposit [a.u]");
    fCH2dSum->SetZTitle("counts");
    fCH2dSum->SetStats(0);
    fCH2dSum->Sumw2();
    
    
    // Add them
    fListHist->Add(fLinearVdriftTest);
    fListHist->Add(fPH2dTest);
    fListHist->Add(fCH2dTest);
    fListHist->Add(fPH2dSM);
    fListHist->Add(fCH2dSM);
    fListHist->Add(fPH2dSum);
    fListHist->Add(fCH2dSum);

  }

  /////////////////////////////////////////
  // Second debug level
  ///////////////////////////////////////
  if(fDebug > 1) {

    fNbGoodTracks = new TH2F(Form("NbGoodTracks_%s",(const char*)fName),"NbGoodTracks",500,0.0,2500.0,200,0.0,100.0);
    fNbGoodTracks->SetXTitle("Nb of good tracks");
    fNbGoodTracks->SetYTitle("Centrality");
    fNbGoodTracks->SetStats(0);

    fNbTRDTrack = new TH1F(Form("TRDTrack_%s",(const char*)fName),"TRDTrack",50,0,50);
    fNbTRDTrack->Sumw2();
    fNbTRDTrackOffline = new TH1F(Form("TRDTrackOffline_%s",(const char*)fName),"TRDTrackOffline",50,0,50);
    fNbTRDTrackOffline->Sumw2();
    fNbTRDTrackStandalone = new TH1F(Form("TRDTrackStandalone_%s",(const char*)fName),"TRDTrackStandalone",50,0,50);
    fNbTRDTrackStandalone->Sumw2();
    fNbTPCTRDtrack = new TH2F(Form("NbTPCTRDtrack_%s",(const char*)fName),"NbTPCTRDtrack",100,0,100,100,0,100);
    fNbTPCTRDtrack->Sumw2();
    //
    fNbTimeBin = new TH1F(Form("NbTimeBin_%s",(const char*)fName),"NbTimeBin",35,0,35);
    fNbTimeBin->Sumw2();
    fNbTimeBinOffline = new TH1F(Form("NbTimeBinOffline_%s",(const char*)fName),"NbTimeBinOffline",35,0,35);
    fNbTimeBinOffline->Sumw2();
    fNbTimeBinStandalone = new TH1F(Form("NbTimeBinStandalone_%s",(const char*)fName),"NbTimeBinStandalone",35,0,35);
    fNbTimeBinStandalone->Sumw2();
    //
    fNbClusters = new TH1F(Form("NbClusters_%s",(const char*)fName),"",35,0,35);
    fNbClusters->Sumw2();
    fNbClustersOffline = new TH1F(Form("NbClustersOffline_%s",(const char*)fName),"",35,0,35);
    fNbClustersOffline->Sumw2();
    fNbClustersStandalone = new TH1F(Form("NbClustersStandalone_%s",(const char*)fName),"",35,0,35);
    fNbClustersStandalone->Sumw2();
    //
    fNbTracklets = new TH1F(Form("NbTracklets_%s",(const char*)fName),"NbTracklets",540,0.,540.);
    fNbTracklets->Sumw2();
    fNbTrackletsOffline = new TH1F(Form("NbTrackletsOffline_%s",(const char*)fName),"NbTrackletsOffline",540,0.,540.);
    fNbTrackletsOffline->Sumw2();
    fNbTrackletsStandalone = new TH1F(Form("NbTrackletsStandalone_%s",(const char*)fName),"NbTrackletsStandalone",540,0.,540.);
    fNbTrackletsStandalone->Sumw2();
   
    fListHist->Add(fNbGoodTracks);
   
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
    if(fOnInstance) {
      fTRDCalibraFillHisto->SetFirstRunGain(fFirstRunGain); // Gain Used
      fTRDCalibraFillHisto->SetVersionGainUsed(fVersionGainUsed); // Gain Used
      fTRDCalibraFillHisto->SetSubVersionGainUsed(fSubVersionGainUsed); // Gain Used
      fTRDCalibraFillHisto->SetFirstRunGainLocal(fFirstRunGainLocal); // Gain Used
      fTRDCalibraFillHisto->SetVersionGainLocalUsed(fVersionGainLocalUsed); // Gain Used
      fTRDCalibraFillHisto->SetSubVersionGainLocalUsed(fSubVersionGainLocalUsed); // Gain Used
      fTRDCalibraFillHisto->SetFirstRunVdrift(fFirstRunVdrift); // Vdrift Used
      fTRDCalibraFillHisto->SetVersionVdriftUsed(fVersionVdriftUsed); // Vdrift Used
      fTRDCalibraFillHisto->SetSubVersionVdriftUsed(fSubVersionVdriftUsed); // Vdrift Used
      if((fFirstRunExB != -1) && (fVersionExBUsed != -1) && (fSubVersionExBUsed != -1)){
	fTRDCalibraFillHisto->SetFirstRunExB(fFirstRunExB); // ExB Used
	fTRDCalibraFillHisto->SetVersionExBUsed(fVersionExBUsed); // ExB Used
	fTRDCalibraFillHisto->SetSubVersionExBUsed(fSubVersionExBUsed); // ExB Used
      }
      fTRDCalibraFillHisto->InitCalDet();
    }
    if(fDebug > 1){
      // title CH2dTest
      TString name("Ver");
      name += fVersionGainUsed;
      name += "Subver";
      name += fSubVersionGainUsed;
      name += "FirstRun";
      name += fFirstRunGain;
      name += "Nz0Nrphi0";
      fCH2dTest->SetTitle(name);  
      // title PH2dTest
      TString namee("Ver");
      namee += fVersionVdriftUsed;
      namee += "Subver";
      namee += fSubVersionVdriftUsed;
      namee += "FirstRun";
      namee += fFirstRunVdrift;
      namee += "Nz0Nrphi0";
      fPH2dTest->SetTitle(namee); 
    }
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
  if(fDebug > 1)  {
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
    
    fFriendTrack = fESDfriend->GetTrack(itrk);
    if(!fFriendTrack)  {
      //printf("No friend track %d\n",itrk);
      continue;
    }

    // Other cuts
    fTrdTrack = 0x0;
    Bool_t good = kTRUE;
    Bool_t standalonetrack = kFALSE;
    Bool_t offlinetrack = kFALSE;
    //ULong_t status = fkEsdTrack->GetStatus();
    
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
      // process chamberstatus
      fTRDChamberStatus->ProcessTrack(fTrdTrack);
    }

    // Quality cuts on the AliESDtrack
    if((fEsdTrackCuts) && (!fEsdTrackCuts->IsSelected((AliVParticle *)fkEsdTrack))) {
      //printf("Not a good track\n");
      continue;
    }
    
    // First Absolute gain calibration
    Int_t trdNTracklets = (Int_t) fkEsdTrack->GetTRDntracklets();
    Int_t trdNTrackletsPID = (Int_t) fkEsdTrack->GetTRDntrackletsPID(); 
    //printf("Number of trd tracklets %d and PID trd tracklets %d\n",trdNTracklets,trdNTrackletsPID);
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
    
    
    if(!fTrdTrack) continue; 

    if(good && fOnInstance) {
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
	
	if((nbclusters > fLow) && (nbclusters < fHigh)){
	  if(fRelativeScale > 0.0) sum = sum/fRelativeScale;
	  fCH2dTest->Fill(sum,detector+0.5);	       
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
		fPH2dTest->Fill((Double_t)(ic/10.0),detector+0.5,(Double_t)phtb[ic]);
		fPH2dSum->Fill((Double_t)(ic/10.0),0.5,(Double_t)phtb[ic]);
		fPH2dSM->Fill((Double_t)(ic/10.0),sector+0.5,(Double_t)phtb[ic]);
	      }
	      else {
		if(phtb[ic] > 0.0) {
		  fPH2dTest->Fill((Double_t)(ic/10.0),detector+0.5,(Double_t)phtb[ic]);
		  fPH2dSum->Fill((Double_t)(ic/10.0),0.0,(Double_t)phtb[ic]);
		  fPH2dSM->Fill((Double_t)(ic/10.0),sector+0.5,(Double_t)phtb[ic]);
		}
	      }
	    }
	  }
	}
	if(detector == 0) FindP1TrackPHtrackletV1Test(tracklet,nbclusters);
	
      } // loop on tracklets
      
      
    } // debug
    
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
//_______________________________________________________________________________________
void  AliTRDCalibTask::AddTask(const AliTRDCalibTask * calibTask) {

  //
  // Add stats
  //

  TList *listcalibTask = calibTask->GetList();
  if(!listcalibTask) return;

  THnSparseI *histoEntries = (THnSparseI *) listcalibTask->FindObject("NumberOfEntries");

  TH1I *nEvents  = (TH1I *) listcalibTask->FindObject(Form("NEvents_%s",(const char*)calibTask->GetName()));
  TH1I *nEventsInput  = (TH1I *) listcalibTask->FindObject(Form("NEventsInput_%s",(const char*)calibTask->GetName()));
  TH2F *absoluteGain  = (TH2F *) listcalibTask->FindObject(Form("AbsoluteGain_%s",(const char*)calibTask->GetName()));

  TH1F *trdTrack = (TH1F *) listcalibTask->FindObject(Form("TRDTrack_%s",(const char*)calibTask->GetName()));
  TH1F *trdTrackOffline = (TH1F *) listcalibTask->FindObject(Form("TRDTrackOffline_%s",(const char*)calibTask->GetName()));
  TH1F *trdTrackStandalone = (TH1F *) listcalibTask->FindObject(Form("TRDTrackStandalone_%s",(const char*)calibTask->GetName()));

  TH2F *tpctrdTrack = (TH2F *) listcalibTask->FindObject(Form("NbTPCTRDtrack_%s",(const char*)calibTask->GetName()));

  TH1F *nbTimeBin = (TH1F *) listcalibTask->FindObject(Form("NbTimeBin_%s",(const char*)calibTask->GetName()));
  TH1F *nbTimeBinOffline = (TH1F *) listcalibTask->FindObject(Form("NbTimeBinOffline_%s",(const char*)calibTask->GetName()));
  TH1F *nbTimeBinStandalone = (TH1F *) listcalibTask->FindObject(Form("NbTimeBinStandalone_%s",(const char*)calibTask->GetName()));

  TH1F *nbClusters = (TH1F *) listcalibTask->FindObject(Form("NbClusters_%s",(const char*)calibTask->GetName()));
  TH1F *nbClustersOffline = (TH1F *) listcalibTask->FindObject(Form("NbClustersOffline_%s",(const char*)calibTask->GetName()));
  TH1F *nbClustersStandalone = (TH1F *) listcalibTask->FindObject(Form("NbClustersStandalone_%s",(const char*)calibTask->GetName()));

  TH1F *nbTracklets = (TH1F *) listcalibTask->FindObject(Form("NbTracklets_%s",(const char*)calibTask->GetName()));
  TH1F *nbTrackletsOffline = (TH1F *) listcalibTask->FindObject(Form("NbTrackletsOffline_%s",(const char*)calibTask->GetName()));
  TH1F *nbTrackletsStandalone = (TH1F *) listcalibTask->FindObject(Form("NbTrackletsStandalone_%s",(const char*)calibTask->GetName()));
  
  TH2I *ch2d = (TH2I *) listcalibTask->FindObject("CH2d");
  TProfile2D *ph2d = (TProfile2D *) listcalibTask->FindObject("PH2d");
  TProfile2D *prf2d = (TProfile2D *) listcalibTask->FindObject("PRF2d");

  TH2I *ch2dSum = (TH2I *) listcalibTask->FindObject(Form("CH2dSum_%s",(const char*)calibTask->GetName()));
  TProfile2D *ph2dSum = (TProfile2D *) listcalibTask->FindObject(Form("PH2dSum_%s",(const char*)calibTask->GetName()));

  TH2I *ch2dSM = (TH2I *) listcalibTask->FindObject(Form("CH2dSM_%s",(const char*)calibTask->GetName()));
  TProfile2D *ph2dSM = (TProfile2D *) listcalibTask->FindObject(Form("PH2dSM_%s",(const char*)calibTask->GetName()));
  
  AliTRDCalibraVdriftLinearFit *linearfit = (AliTRDCalibraVdriftLinearFit *) listcalibTask->FindObject("AliTRDCalibraVdriftLinearFit");  
  AliTRDCalibraExbAltFit *exbaltfit = (AliTRDCalibraExbAltFit *) listcalibTask->FindObject("AliTRDCalibraExbAltFit");  
  AliTRDCalibraVector *calibraVector = (AliTRDCalibraVector *) listcalibTask->FindObject("AliTRDCalibraVector");  

  //

  THnSparseI *inhistoEntries = (THnSparseI *) fListHist->FindObject("NumberOfEntries");

  TH1I *inEventsInput  = (TH1I *) fListHist->FindObject(Form("NEventsInput_%s",(const char*)fName));
  TH1I *inEvents  = (TH1I *) fListHist->FindObject(Form("NEvents_%s",(const char*)fName));
  TH2F *iabsoluteGain  = (TH2F *) fListHist->FindObject(Form("AbsoluteGain_%s",(const char*)fName));

  TH1F *itrdTrack = (TH1F *) fListHist->FindObject(Form("TRDTrack_%s",(const char*)fName));
  TH1F *itrdTrackOffline = (TH1F *) fListHist->FindObject(Form("TRDTrackOffline_%s",(const char*)fName));
  TH1F *itrdTrackStandalone = (TH1F *) fListHist->FindObject(Form("TRDTrackStandalone_%s",(const char*)fName));

  TH2F *itpctrdTrack = (TH2F *) fListHist->FindObject(Form("NbTPCTRDtrack_%s",(const char*)fName));

  TH1F *inbTimeBin = (TH1F *) fListHist->FindObject(Form("NbTimeBin_%s",(const char*)fName));
  TH1F *inbTimeBinOffline = (TH1F *) fListHist->FindObject(Form("NbTimeBinOffline_%s",(const char*)fName));
  TH1F *inbTimeBinStandalone = (TH1F *) fListHist->FindObject(Form("NbTimeBinStandalone_%s",(const char*)fName));

  TH1F *inbClusters = (TH1F *) fListHist->FindObject(Form("NbClusters_%s",(const char*)fName));
  TH1F *inbClustersOffline = (TH1F *) fListHist->FindObject(Form("NbClustersOffline_%s",(const char*)fName));
  TH1F *inbClustersStandalone = (TH1F *) fListHist->FindObject(Form("NbClustersStandalone_%s",(const char*)fName));

  TH1F *inbTracklets = (TH1F *) fListHist->FindObject(Form("NbTracklets_%s",(const char*)fName));
  TH1F *inbTrackletsOffline = (TH1F *) fListHist->FindObject(Form("NbTrackletsOffline_%s",(const char*)fName));
  TH1F *inbTrackletsStandalone = (TH1F *) fListHist->FindObject(Form("NbTrackletsStandalone_%s",(const char*)fName));
  
  TH2I *ich2d = (TH2I *) fListHist->FindObject("CH2d");
  TProfile2D *iph2d = (TProfile2D *) fListHist->FindObject("PH2d");
  TProfile2D *iprf2d = (TProfile2D *) fListHist->FindObject("PRF2d");

  TH2I *ich2dSum = (TH2I *) fListHist->FindObject(Form("CH2dSum_%s",(const char*)fName));
  TProfile2D *iph2dSum = (TProfile2D *) fListHist->FindObject(Form("PH2dSum_%s",(const char*)fName));

  TH2I *ich2dSM = (TH2I *) fListHist->FindObject(Form("CH2dSM_%s",(const char*)fName));
  TProfile2D *iph2dSM = (TProfile2D *) fListHist->FindObject(Form("PH2dSM_%s",(const char*)fName));
  
  AliTRDCalibraVdriftLinearFit *ilinearfit = (AliTRDCalibraVdriftLinearFit *) fListHist->FindObject("AliTRDCalibraVdriftLinearFit");  
  AliTRDCalibraExbAltFit *iexbaltfit = (AliTRDCalibraExbAltFit *) fListHist->FindObject("AliTRDCalibraExbAltFit");
  AliTRDCalibraVector *icalibraVector = (AliTRDCalibraVector *) fListHist->FindObject("AliTRDCalibraVector");  


  // Add

  if(histoEntries) {
    if(inhistoEntries) {
      inhistoEntries->Add(histoEntries);
      //printf("Add Events\n");
    }
    else {
      //printf("Create new Events\n");
      inhistoEntries = (THnSparseI *) histoEntries->Clone();
      fListHist->Add(inhistoEntries);
    }
  }

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

  if(exbaltfit) {
    if(iexbaltfit) iexbaltfit->Add(exbaltfit);
    else {
      iexbaltfit = new AliTRDCalibraExbAltFit(*exbaltfit);
      fListHist->Add(iexbaltfit);
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
    } else if(os->GetString().Contains("TRD/Calib/ChamberExB")){
      // Get Old drift velocity calibration
      AliCDBId *id=AliCDBId::MakeFromString(os->GetString());
      fFirstRunExB = id->GetFirstRun();
      fVersionExBUsed = id->GetVersion();
      fSubVersionExBUsed = id->GetSubVersion();
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
//______________________________________________________________________________________________________________________
Bool_t AliTRDCalibTask::FindP1TrackPHtrackletV1Test(const AliTRDseedV1 *tracklet, Int_t nbclusters)
{
  //
  // Drift velocity calibration:
  // Fit the clusters with a straight line
  // From the slope find the drift velocity
  //

  ////////////////////////////////////////////////
  //Number of points: if less than 3 return kFALSE
  /////////////////////////////////////////////////
  if(nbclusters <= 2) return kFALSE;

  ////////////
  //Variables
  ////////////
  // results of the linear fit
  Double_t dydt                       = 0.0;                                // dydt tracklet after straight line fit
  Double_t errorpar                   = 0.0;                                // error after straight line fit on dy/dt
  Double_t pointError                 = 0.0;                                // error after straight line fit 
  // pad row problemes: avoid tracklet that cross pad rows, tilting angle in the constant
  Int_t    crossrow                   = 0;                                  // if it crosses a pad row
  Int_t    rowp                       = -1;                                 // if it crosses a pad row
  Float_t  tnt                        = tracklet->GetTilt();                // tan tiltingangle
  TLinearFitter linearFitterTracklet(2,"pol1");
  linearFitterTracklet.StoreData(kTRUE);  
 
  
  ///////////////////////////////////////////
  // Take the parameters of the track
  //////////////////////////////////////////
  // take now the snp, tnp and tgl from the track
  Double_t snp = tracklet->GetSnp();             // sin dy/dx at the end of the chamber
  Double_t tnp = 0.0;                            // dy/dx at the end of the chamber 
  if( TMath::Abs(snp) <  1.){
    tnp = snp / TMath::Sqrt((1.-snp)*(1.+snp));
  } 
  Double_t tgl  = tracklet->GetTgl();           // dz/dl
  Double_t dzdx = tgl*TMath::Sqrt(1+tnp*tnp);   // dz/dx calculated from dz/dl
  // at the entrance
  //Double_t tnp = tracklet->GetYref(1);      // dy/dx at the entrance of the chamber
  //Double_t tgl = tracklet->GetZref(1);      // dz/dl at the entrance of the chamber
  //Double_t dzdx = tgl;                      //*TMath::Sqrt(1+tnp*tnp); // dz/dx from dz/dl
  // at the end with correction due to linear fit
  //Double_t tnp = tracklet->GetYfit(1);      // dy/dx at the end of the chamber after fit correction
  //Double_t tgl = tracklet->GetZfit(1);      // dz/dl at the end of the chamber after fit correction 


  ////////////////////////////
  // loop over the clusters
  ////////////////////////////
  Int_t  nbli = 0;
  AliTRDcluster *cl                   = 0x0;
  //////////////////////////////
  // Check no shared clusters
  //////////////////////////////
  for(int icc=AliTRDseedV1::kNtb; icc<AliTRDseedV1::kNclusters; icc++){
    cl = tracklet->GetClusters(icc);
    if(cl)  crossrow = 1;
  }
  //////////////////////////////////
  // Loop clusters
  //////////////////////////////////
  for(int ic=0; ic<AliTRDseedV1::kNtb; ic++){
    if(!(cl = tracklet->GetClusters(ic))) continue;
    //if((fLimitChargeIntegration) && (!cl->IsInChamber())) continue;
    
    Double_t ycluster                 = cl->GetY();
    Int_t time                        = cl->GetPadTime();
    Double_t timeis                   = time/10.0;
    //See if cross two pad rows
    Int_t    row                      = cl->GetPadRow();
    if(rowp==-1) rowp                 = row;
    if(row != rowp) crossrow          = 1;

    linearFitterTracklet.AddPoint(&timeis,ycluster,1);
    nbli++;  

    
  }
  
  ////////////////////////////////////
  // Do the straight line fit now
  ///////////////////////////////////
  if(nbli <= 2){ 
    linearFitterTracklet.ClearPoints();  
    return kFALSE; 
  }
  TVectorD pars;
  linearFitterTracklet.Eval();
  linearFitterTracklet.GetParameters(pars);
  pointError  =  TMath::Sqrt(linearFitterTracklet.GetChisquare()/(nbli-2));
  errorpar    =  linearFitterTracklet.GetParError(1)*pointError;
  dydt        = pars[1]; 
  //printf("chis %f, nbli %d, pointError %f, parError %f, errorpar %f\n",linearFitterTracklet->GetChisquare(),nbli,pointError,linearFitterTracklet->GetParError(1),errorpar);
  linearFitterTracklet.ClearPoints();  
 
  /////////////////////////
  // Cuts quality
  ////////////////////////
  
  if(nbclusters < fLow) return kFALSE;
  if(nbclusters > fHigh) return kFALSE;
  if(pointError >= 0.3) return kFALSE;
  if(crossrow == 1) return kTRUE;
  
  ///////////////////////
  // Fill
  //////////////////////

  if(fDebug > 0){
    //Add to the linear fitter of the detector
    if( TMath::Abs(snp) <  1.){
      Double_t x = tnp-dzdx*tnt; 
      //if(!fLinearVdriftTest) printf("Not there\n");
      Double_t nbentries = fLinearVdriftTest->GetEntries();
      if(nbentries < (5.0*32767)) fLinearVdriftTest->Fill(x,dydt);
    }
  }
  
  return kTRUE;
}

