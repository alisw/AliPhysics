
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

// macro for very simple analysis



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

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
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

#include "AliTRDCalibraFillHisto.h"
#include "AliTRDCalibraVdriftLinearFit.h" 

#include "AliCDBManager.h"
#include "AliTRDcalibDB.h"


#include "AliTRDCalibTask.h"


ClassImp(AliTRDCalibTask)

//________________________________________________________________________
  AliTRDCalibTask::AliTRDCalibTask(const char *name) 
    : AliAnalysisTask(name, ""), fESD(0), fESDfriend(0),
      fEsdTrack(0),
      fFriendTrack(0),
      fCalibObject(0),
      fTrdTrack(0),
      fCl(0),
      fListHist(0),
      fTRDCalibraFillHisto(0),
      fNEvents(0),
      fNbTRDTrack(0),
      fNbTRDTrackOffline(0),
      fNbTRDTrackStandalone(0),
      fNbTPCTRDtrack(0),
      fNbTimeBin(0),
      fNbTimeBinOffline(0),
      fNbTimeBinStandalone(0),
      fNbClusters(0),
      fNbClustersOffline(0),
      fNbClustersStandalone(0),
      fNbTracklets(0),
      fNbTrackletsOffline(0),
      fNbTrackletsStandalone(0),
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
      fMinNbContributors(0),
      fLow(0),
      fHigh(30),
      fFillZero(kFALSE),
      fNormalizeNbOfCluster(kFALSE),
      fRelativeScale(0.0),
      fMaxCluster(100.0),
      fNbMaxCluster(2),
      fOfflineTracks(kFALSE),
      fStandaloneTracks(kFALSE),
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
  DefineOutput(0, TList::Class());
  
   
}
//____________________________________________________________________________________
AliTRDCalibTask::~AliTRDCalibTask()
{
  //
  // AliTRDCalibTask destructor
  //

  // Pointeur
  if(fNEvents) delete fNEvents;
  if(fNbTRDTrack) delete fNbTRDTrack;
  if(fNbTRDTrackOffline) delete fNbTRDTrackOffline;
  if(fNbTRDTrackStandalone) delete fNbTRDTrackStandalone;
  if(fNbTPCTRDtrack) delete fNbTPCTRDtrack;
  if(fNbTimeBin) delete fNbTimeBin;
  if(fNbTimeBinOffline) delete fNbTimeBinOffline;
  if(fNbTimeBinStandalone) delete fNbTimeBinStandalone;
  if(fNbClusters) delete fNbClusters;
  if(fNbClustersOffline) delete fNbClustersOffline;
  if(fNbClustersStandalone) delete fNbClustersStandalone;
  if(fNbTracklets) delete fNbTracklets;
  if(fNbTrackletsOffline) delete fNbTrackletsOffline;
  if(fNbTrackletsStandalone) delete fNbTrackletsStandalone;
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
void AliTRDCalibTask::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once
  
  TTree* tree = dynamic_cast<TTree*> (GetInputData(0)); //pointer wird "umgecastet" auf anderen Variablentyp
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    
    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else {
      fESD = esdH->GetEvent();
      esdH->SetActiveBranches("ESDfriend*");
      Printf("*** CONNECTED NEW EVENT ****");
    }
    
  }
}
//________________________________________________________________________
void AliTRDCalibTask::CreateOutputObjects() 
{
  

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
  if(fHisto2d) {  
    fListHist->Add(fTRDCalibraFillHisto->GetCH2d());
    fListHist->Add(fTRDCalibraFillHisto->GetPH2d()); 
    fListHist->Add(fTRDCalibraFillHisto->GetPRF2d());
  } 
  if(fVdriftLinear) fListHist->Add((TObject *) fTRDCalibraFillHisto->GetVdriftLinearFit());
  if(fVector2d) fListHist->Add((TObject *) fTRDCalibraFillHisto->GetCalibraVector()); //calibra vector  
  fNEvents = new TH1I("NEvents","NEvents", 2, 0, 2);
  fListHist->Add(fNEvents);
  
  /////////////////////////////////////////
  // First debug level
  ///////////////////////////////////////
  if(fDebug > 0) {
    
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
    
    // Add them
    fListHist->Add(fPH2dSM);
    fListHist->Add(fCH2dSM);
    fListHist->Add(fPH2dSum);
    fListHist->Add(fCH2dSum);
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
  
 }
 //________________________________________________________________________
 void AliTRDCalibTask::Exec(Option_t *) 
 {
   AliLog::SetGlobalLogLevel(AliLog::kError);

   if (!fESD) {
     Printf("ERROR: fESD not available");
     return;
   }

   fCounter++;
   if((fMaxEvent != 0) && (fMaxEvent < fCounter)) return;

   ///////////////////////////////
   // Require a primary vertex
   //////////////////////////////
   if(fRequirePrimaryVertex) {
     //printf("test\n");
     const AliESDVertex* vtxESD = 0x0;
     vtxESD = fESD->GetPrimaryVertexTPC() ;
     if(!vtxESD){
       //printf("No primary at all vertex\n");
       return;
     }
     Int_t nCtrb = vtxESD->GetNContributors();
     //printf("Number of contributors %d\n",nCtrb);
     if(nCtrb < fMinNbContributors) return;
   }

   //printf("Primary vertex passed\n");

   ///////////////////
   // Check trigger
   ///////////////////
   //printf("Class Fired %s\n",(const char*)fESD->GetFiredTriggerClasses());
   Bool_t pass = kTRUE;
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
   if(!pass) return;   

   //printf("Trigger passed\n");
  
   fNEvents->Fill(1);

   // In total
   Int_t nbTrdTracks = 0;
   // standalone
   Int_t nbTrdTracksStandalone = 0;
   // offline
   Int_t nbTrdTracksOffline = 0;
   // TPC
   Int_t nbtrackTPC = 0;

   Double_t nbTracks = fESD->GetNumberOfTracks();
   //printf("Number of tracks %f\n",nbTracks);  

   if (nbTracks <= 0.0) {
     
     if(fDebug > 1) {
       fNbTRDTrack->Fill(nbTrdTracks);
       fNbTRDTrackStandalone->Fill(nbTrdTracksStandalone);
       fNbTRDTrackOffline->Fill(nbTrdTracksOffline);
     }
     PostData(0, fListHist);
     return;
   }

   fESDfriend = (AliESDfriend *)fESD->FindListObject("AliESDfriend");
   fESD->SetESDfriend(fESDfriend);
   if(!fESDfriend) return;   
   
   //printf("has friends\n");

   ////////////////////////////////////
   // Check the number of TPC tracks
   ///////////////////////////////////
   for(Int_t itrk = 0; itrk < nbTracks; itrk++){
     // Get ESD track
     fEsdTrack = fESD->GetTrack(itrk);
     ULong_t status = fEsdTrack->GetStatus(); 
     if(status&(AliESDtrack::kTPCout)) nbtrackTPC++;
     // Check that the calibration object is here
     if(fESDfriend && (fESDfriend->GetTrack(itrk))) {
       fFriendTrack = new AliESDfriendTrack(*(fESDfriend->GetTrack(itrk)));
       Int_t counteer = 0;
       Int_t icalib=0;
       while((fCalibObject = (TObject *)(fFriendTrack->GetCalibObject(icalib++)))){
	 if(strcmp(fCalibObject->IsA()->GetName(), "AliTRDtrackV1") != 0) continue;
	 counteer++;
       }
       if(counteer > 0) {
	 nbTrdTracks++;
	 if((status&(AliESDtrack::kTRDout)) && (!(status&(AliESDtrack::kTRDin)))) {
	   nbTrdTracksStandalone++;
	 }
	 if((status&(AliESDtrack::kTRDin))) {
	   nbTrdTracksOffline++;
	 }
       }
       if(fFriendTrack) delete fFriendTrack;
     }
   }
   //printf("Number of TPC tracks %d, TRD %d\n",nbtrackTPC,nbTrdTracks);
   if((nbtrackTPC>0) && (nbTrdTracks > (3.0*nbtrackTPC))) pass = kFALSE;
   
   if(fDebug > 1) {
     
     fNbTRDTrack->Fill(nbTrdTracks);
     fNbTRDTrackStandalone->Fill(nbTrdTracksStandalone);
     fNbTRDTrackOffline->Fill(nbTrdTracksOffline);
     fNbTPCTRDtrack->Fill(nbTrdTracks,nbtrackTPC);
   
   }

   if(!pass) {
     PostData(0, fListHist);
     return;
   }
  
   /////////////////////////////////////
   // Loop on AliESDtrack
   ////////////////////////////////////
      
   for(int itrk=0; itrk < nbTracks; itrk++){

     // Get ESD track
     fEsdTrack = fESD->GetTrack(itrk);
     if(!fEsdTrack) continue;

     // Quality cuts on the AliESDtrack
     if((fEsdTrackCuts) && (!fEsdTrackCuts->IsSelected((AliVParticle *)fEsdTrack))) {
       //printf("Not a good track\n");
       continue;
     }
     
     // Other cuts
     Bool_t good = kTRUE;
     Bool_t standalonetrack = kFALSE;
     Bool_t offlinetrack = kFALSE;
     ULong_t status = fEsdTrack->GetStatus();
     
     if(!(fESDfriend->GetTrack(itrk)))  continue;   
     
     fFriendTrack = new AliESDfriendTrack(*(fESDfriend->GetTrack(itrk)));
     
     //////////////////////////////////////
     // Loop on calibration objects
     //////////////////////////////////////
     Int_t icalib=0;
     while((fCalibObject = (TObject *)(fFriendTrack->GetCalibObject(icalib++)))){
       if(strcmp(fCalibObject->IsA()->GetName(), "AliTRDtrackV1") != 0) continue;

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
       if(good) fTRDCalibraFillHisto->UpdateHistogramsV1(fTrdTrack);
       
       //////////////////////////////////
       // Debug 
       ////////////////////////////////

       if(fDebug > 0) {
	 
	 
	 Int_t nbtracklets = 0;
	 
	 // Check some stuff
	 Bool_t standalonetracklet = kFALSE;  
	 const AliTRDseedV1 *tracklet = 0x0;
	 //////////////////////////////////////
	 // Loop tracklets
	 ///////////////////////////////////// 
	 for(Int_t itr = 0; itr < 6; itr++){
	   
	   if(!(tracklet = fTrdTrack->GetTracklet(itr))) continue;
	   if(!tracklet->IsOK()) continue;
	   nbtracklets++;
	   standalonetracklet = kFALSE; 
	   if(tracklet->IsStandAlone()) standalonetracklet = kTRUE;

	   Int_t nbclusters = 0;
	   Double_t phtb[AliTRDseedV1::kNtb];
	   memset(phtb, 0, AliTRDseedV1::kNtb*sizeof(Double_t));
	   Double_t sum = 0.0;
	   Float_t normalisation = 6.67;
	   Int_t detector = 0;
	   Int_t sector = 0;
	   //Int_t crossrow = 0;
	   
	   // Check no shared clusters
	   //for(int icc=AliTRDseedV1::kNtb; icc<AliTRDseedV1::kNclusters; icc++){
	   //  if((fcl = tracklet->GetClusters(icc)))  crossrow = 1;
	   // }
	   
	   // Loop on clusters
	   for(int ic=0; ic<AliTRDseedV1::kNtb; ic++){
	     
	     if(!(fCl = tracklet->GetClusters(ic))) continue;
	     nbclusters++;
	     Int_t time = fCl->GetPadTime();
	     Float_t ch =  tracklet->GetdQdl(ic);
	     Float_t qcl = TMath::Abs(fCl->GetQ());
	     detector = fCl->GetDetector();	  
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
     
     delete fFriendTrack;
     
   } // loop ESD track
   
   // Post output data
   PostData(0, fListHist);
 }     
//________________________________________________________________________
void AliTRDCalibTask::Terminate(Option_t *) 
{
  
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
	  histolinearfitsum ->Add(linearfit->GetLinearFitterHisto(det));
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
