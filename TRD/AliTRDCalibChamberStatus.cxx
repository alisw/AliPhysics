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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliTRDCalibChamberStatus: to determine which half chambers are off     //
// Produce a AliTRDCalChamberStatus calibration object                    //
// Check with the AliTRDCalDCSFEEv2 info                                  //
//                                                                        //
//                                                                        //
// Authors:                                                               //
//   J. Book (jbook@ikf.uni-frankfurt.de)                                 //
//   R. Bailhache (rbailhache@ikf.uni-frankfurt.de)                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


//Root includes
#include <THnSparse.h>

#include <TDirectory.h>
#include <TFile.h>
#include <TAxis.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TCanvas.h>

//AliRoot includes
#include "AliRawReader.h"

//header file
#include "AliLog.h"
#include "AliTRDCalibChamberStatus.h"
#include "AliTRDgeometry.h"
#include "AliTRDfeeParam.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDSignalIndex.h"
#include "AliTRDpadPlane.h"
#include "./Cal/AliTRDCalChamberStatus.h"
#include "./Cal/AliTRDCalDCSv2.h"
#include "./Cal/AliTRDCalDCSFEEv2.h"

#include "AliTRDrawStream.h"
#include "AliTRDseedV1.h"
#include "AliTRDcluster.h"

#ifdef ALI_DATE
#include "event.h"
#endif

ClassImp(AliTRDCalibChamberStatus) /*FOLD00*/

//_____________________________________________________________________
AliTRDCalibChamberStatus::AliTRDCalibChamberStatus() : /*FOLD00*/
  TObject(),
  fDetector(-1),
  fNumberOfTimeBins(0),
  fCounterEventNotEmpty(0),
  fCalChamberStatus(0x0),
  fHnSparseI(0x0),
  fHnSparseHCM(0x0),
  fHnSparseEvtDet(0x0),
  fHnSparseDebug(0x0),
  fHnSparseMCM(0x0),
  fC1(0x0),
  fDebugLevel(0)
{
    //
    // default constructor
    //

}
//_____________________________________________________________________
AliTRDCalibChamberStatus::AliTRDCalibChamberStatus(const AliTRDCalibChamberStatus &ped) : /*FOLD00*/
  TObject(ped),
  fDetector(ped.fDetector),
  fNumberOfTimeBins(ped.fNumberOfTimeBins),
  fCounterEventNotEmpty(ped.fCounterEventNotEmpty),
  fCalChamberStatus(ped.fCalChamberStatus),
  fHnSparseI(ped.fHnSparseI),
  fHnSparseHCM(ped.fHnSparseHCM),
  fHnSparseEvtDet(ped.fHnSparseEvtDet),
  fHnSparseDebug(ped.fHnSparseDebug),  
  fHnSparseMCM(ped.fHnSparseMCM),
  fC1(ped.fC1),
  fDebugLevel(ped.fDebugLevel)
{
    //
    // copy constructor
    //
  
}
//_____________________________________________________________________
AliTRDCalibChamberStatus& AliTRDCalibChamberStatus::operator = (const  AliTRDCalibChamberStatus &source)
{
  //
  // assignment operator
  //
  if (&source == this) return *this;
  new (this) AliTRDCalibChamberStatus(source);

  return *this;
}
//_____________________________________________________________________
AliTRDCalibChamberStatus::~AliTRDCalibChamberStatus() /*FOLD00*/
{
  //
  // destructor
  //
  if(fCalChamberStatus){
    delete fCalChamberStatus;
  }
  if(fHnSparseI) {
    delete fHnSparseI;
  }
  if(fHnSparseHCM) {
    delete fHnSparseHCM;
  }
  if(fHnSparseEvtDet) {
    delete fHnSparseEvtDet;
  }
  if(fHnSparseDebug) {
    delete fHnSparseDebug;
  }
  if(fHnSparseMCM) {
   delete fHnSparseMCM;
  }
  if(fC1) {
   delete fC1;
  }
}

//_____________________________________________________________________
void AliTRDCalibChamberStatus::Init() 
{
  //
  // Init the different THnSparse
  //
  //

  //
  // Init the fHnSparseI
  //

  //create the map
  Int_t thnDimEvt[4]; // sm, layer, stack, halfchamber
  thnDimEvt[0] = 18;
  thnDimEvt[1] = 6;
  thnDimEvt[2] = 5;
  thnDimEvt[3] = 2;
  //arrays for lower bounds :
  Double_t* binEdgesEvt[4];
  for(Int_t ivar = 0; ivar < 4; ivar++)
    binEdgesEvt[ivar] = new Double_t[thnDimEvt[ivar] + 1];
  //values for bin lower bounds
  for(Int_t i=0; i<=thnDimEvt[0]; i++) binEdgesEvt[0][i]= 0.0  + (18.0)/thnDimEvt[0]*(Double_t)i;
  for(Int_t i=0; i<=thnDimEvt[1]; i++) binEdgesEvt[1][i]= 0.0  + (6.0)/thnDimEvt[1]*(Double_t)i;
  for(Int_t i=0; i<=thnDimEvt[2]; i++) binEdgesEvt[2][i]= 0.0  + (5.0)/thnDimEvt[2]*(Double_t)i;
  for(Int_t i=0; i<=thnDimEvt[3]; i++) binEdgesEvt[3][i]= 0.0  + (2.0)/thnDimEvt[3]*(Double_t)i;
  
  //create the THnSparse
  fHnSparseI = new THnSparseI("NumberOfEntries","NumberOfEntries",4,thnDimEvt);
  for (int k=0; k<4; k++) {
    fHnSparseI->SetBinEdges(k,binEdgesEvt[k]);
  }
  fHnSparseI->Sumw2();

  //
  // Init the fHnSparseHCM (THnSparseI)
  //

  //create the THnSparse
  fHnSparseHCM = new THnSparseI("HCMerrors","HCMerrors",4,thnDimEvt);
  for (int k=0; k<4; k++) {
    fHnSparseHCM->SetBinEdges(k,binEdgesEvt[k]);
  }
  fHnSparseHCM->Sumw2();


  //---------//
  //  Debug  //
  if(fDebugLevel > 0) {
  
    //
    // Init the fHnSparseEvtDet (THnSparseI)
    //
    
    //create the map
    Int_t thnDimEvts[3]; // event, detector, halfchamber
    thnDimEvts[0] = 10000;
    thnDimEvts[1] = 540;
    thnDimEvts[2] = 2;
    //arrays for lower bounds :
    Double_t* binEdgesEvts[3];
    for(Int_t ivar = 0; ivar < 3; ivar++)
      binEdgesEvts[ivar] = new Double_t[thnDimEvts[ivar] + 1];
    //values for bin lower bounds
    for(Int_t i=0; i<=thnDimEvts[0]; i++) binEdgesEvts[0][i]= 0.0  + (10000.0)/thnDimEvts[0]*(Double_t)i;
    for(Int_t i=0; i<=thnDimEvts[1]; i++) binEdgesEvts[1][i]= 0.0  + (540.0)/thnDimEvts[1]*(Double_t)i;
    for(Int_t i=0; i<=thnDimEvts[2]; i++) binEdgesEvts[2][i]= 0.0  + (2.0)/thnDimEvts[2]*(Double_t)i;
    
    //create the THnSparse
    fHnSparseEvtDet = new THnSparseI("NumberOfEntriesPerEvent","NumberOfEntriesPerEvent",3,thnDimEvts);
    for (int k=0; k<3; k++) {
      fHnSparseEvtDet->SetBinEdges(k,binEdgesEvts[k]);
    }
    fHnSparseEvtDet->Sumw2();

    //
    // Init the fHnSparseDebug (THnSparseI)
    //
    
    //create the THnSparse
    fHnSparseDebug = new THnSparseI("NumberOfDifferentDecisions","NumberOfDifferentDecisions",4,thnDimEvt);
    for (int k=0; k<4; k++) {
      fHnSparseDebug->SetBinEdges(k,binEdgesEvt[k]);
    }
    fHnSparseDebug->Sumw2();
        
    //
    // Init the fHnSparseMCM (THnSparseI)
    //
    
    //create the map
    Int_t thnDimEvtt[6]; // sm, layer, stack, ROB, MCM
    thnDimEvtt[0] = 18;
    thnDimEvtt[1] = 6;
    thnDimEvtt[2] = 5;
    thnDimEvtt[3] = 8;
    thnDimEvtt[4] = 18;
    thnDimEvtt[5] = 18;
    //arrays for lower bounds :
    Double_t* binEdgesEvtt[6];
    for(Int_t ivar = 0; ivar < 6; ivar++)
      binEdgesEvtt[ivar] = new Double_t[thnDimEvtt[ivar] + 1];
    //values for bin lower bounds
    for(Int_t i=0; i<=thnDimEvtt[0]; i++) binEdgesEvtt[0][i]= 0.0  + (18.0)/thnDimEvtt[0]*(Double_t)i;
    for(Int_t i=0; i<=thnDimEvtt[1]; i++) binEdgesEvtt[1][i]= 0.0  + (6.0)/thnDimEvtt[1]*(Double_t)i;
    for(Int_t i=0; i<=thnDimEvtt[2]; i++) binEdgesEvtt[2][i]= 0.0  + (5.0)/thnDimEvtt[2]*(Double_t)i;
    for(Int_t i=0; i<=thnDimEvtt[3]; i++) binEdgesEvtt[3][i]= 0.0  + (8.0)/thnDimEvtt[3]*(Double_t)i;
    for(Int_t i=0; i<=thnDimEvtt[4]; i++) binEdgesEvtt[4][i]= 0.0  + (18.0)/thnDimEvtt[4]*(Double_t)i;
    for(Int_t i=0; i<=thnDimEvtt[5]; i++) binEdgesEvtt[5][i]= 0.0  + (18.0)/thnDimEvtt[5]*(Double_t)i;
    
    //create the THnSparse
    fHnSparseMCM = new THnSparseI("MCMerrorDCS","MCMerrorDCS",6,thnDimEvtt);
    for (int k=0; k<6; k++) {
      fHnSparseMCM->SetBinEdges(k,binEdgesEvtt[k]);
    }
    fHnSparseMCM->Sumw2();
  
  }
  //  Debug  //
  //---------//

}
//_____________________________________________________________________
void AliTRDCalibChamberStatus::ProcessTrack(AliTRDtrackV1 * trdTrack)
{
  //
  // Track Processing to get half chamber status
  //
  //
  
  
  const AliTRDseedV1 *tracklet = 0x0;
  AliTRDcluster *cluster;
  //////////////////////////////////////
  // Loop tracklets
  ///////////////////////////////////// 
  for(Int_t itr = 0; itr < 6; ++itr){
	
    if(!(tracklet = trdTrack->GetTracklet(itr))) continue;
    if(!tracklet->IsOK()) continue;
  
    // Loop on clusters
    for(int ic=0; ic<AliTRDseedV1::kNtb; ++ic){
      if((cluster = tracklet->GetClusters(ic))) {
	//printf("ic %d\n",ic);
	break;
      }
    }
    if(!cluster) continue;
    
    Int_t det     = cluster->GetDetector();	  
    Int_t layer   = AliTRDgeometry::GetLayer(det);
    Int_t sm      = AliTRDgeometry::GetSector(det);
    Int_t stac    = AliTRDgeometry::GetStack(det);
    
    Int_t col     = cluster->GetPadCol();
    Int_t iMcm    = (Int_t)(col/18);   // current group of 18 col pads
    Double_t rphi = 0.5;
    if(iMcm > 3) rphi = 1.5;
	  
    Double_t val[4] = {sm,layer,stac,rphi}; 
    if(fHnSparseI->GetBinContent((const Int_t*)val)<2147483646) fHnSparseI->Fill(&val[0]); 
  }
  
}
//_____________________________________________________________________
void AliTRDCalibChamberStatus::ProcessEvent(AliRawReader * rawReader, Int_t nevents_physics)
{
  //
  // Event Processing loop with AliTRDrawStream
  //
  //
  
  Bool_t notEmpty = kFALSE;
    
  AliTRDdigitsManager *digitsManager = new AliTRDdigitsManager(kTRUE);
  digitsManager->CreateArrays();
  
  AliTRDrawStream *rawStream = new AliTRDrawStream(rawReader);
  rawStream->SetDigitsManager(digitsManager);
  //  rawStream->SetNoErrorWarning();
  //  rawStream->SetSharedPadReadout(kFALSE);

  
  Int_t det    = 0;
  while ((det = rawStream->NextChamber(digitsManager, NULL, NULL)) >= 0) { 

    //nextchamber loop
    
    // do the QA analysis
    if (digitsManager->GetIndexes(det)->HasEntry()) {//QA
      // printf("there is ADC data on this chamber!\n");
      
      AliTRDSignalIndex *indexes = digitsManager->GetIndexes(det);
      if (indexes->IsAllocated() == kFALSE) {
	// AliError("Indexes do not exist!");
	break;
      }
      
      Int_t iRow  = 0;
      Int_t iCol  = 0;
      indexes->ResetCounters();
      
      while (indexes->NextRCIndex(iRow, iCol)){
	Int_t iMcm        = (Int_t)(iCol/18);   // current group of 18 col pads
	
	Int_t layer = AliTRDgeometry::GetLayer(det);
	Int_t sm    = AliTRDgeometry::GetSector(det);
	Int_t stac  = AliTRDgeometry::GetStack(det);
	Double_t rphi = 0.5;
	if(iMcm > 3) rphi = 1.5;

	Double_t val[4] = {sm,layer,stac,rphi}; 
	fHnSparseI->Fill(&val[0]); 
	notEmpty = kTRUE;
	
	//---------//
	//  Debug  //
	if(fDebugLevel > 0) {
	  Int_t detector = AliTRDgeometry::GetDetector(layer,stac,sm);
	  Double_t valu[3] = {nevents_physics,detector,rphi};
	  fHnSparseEvtDet->Fill(&valu[0]); 
	}
	//  Debug  //
	//---------//
      }
      
    }
    digitsManager->ClearArrays(det);
  }

  if(notEmpty) fCounterEventNotEmpty++;

  delete digitsManager;
  delete rawStream;
   
}
//_____________________________________________________________________
Bool_t AliTRDCalibChamberStatus::TestEventHisto(Int_t nevent) /*FOLD00*/
{
  //
  //  Test event loop
  // fill the fHnSparseI with entries
  //
  
  AliTRDgeometry geo;


  for(Int_t ievent=0; ievent<nevent; ievent++){
    for (Int_t ism=0; ism<18; ism++){
      for (Int_t istack=0; istack<5; istack++){
	for (Int_t ipl=0; ipl<6; ipl++){
	  for (Int_t icol=0; icol<geo.GetColMax(ipl); icol++){
	    Int_t side = 0;
	    if(icol > 72) side = 1;
	    Double_t val[4] = {ism,ipl,istack,side}; 
	    fHnSparseI->Fill(&val[0]); 
	  }
	}
      }
    }
  }
  
  return kTRUE;

}
//_____________________________________________________________________
void AliTRDCalibChamberStatus::AnalyseHisto(Int_t limit) /*FOLD00*/
{
  //
  //  Create the AliTRDCalChamberStatus according to the fHnSparseI
  //
  
  if(fCalChamberStatus) delete fCalChamberStatus;
  fCalChamberStatus = new AliTRDCalChamberStatus();

  // Check if enough events/tracklets per halfchamber to say something
  Double_t mean=0.0; //number of tracklets per HCS
  Int_t coord2[4];
  for(Int_t bin = 0; bin < fHnSparseI->GetNbins(); bin++) {
    //if(fHnSparseI->GetBinContent(bin,coord2)==0.0) printf(" bin shouldnt be empty!!\n");
    mean+=fHnSparseI->GetBinContent(bin,coord2);
  }
  mean/=fHnSparseI->GetNbins();
  //printf(" mean tracklets per halfchamber %f \n",mean);
  if((fCounterEventNotEmpty < limit) && (mean < limit)) {
    // Say all good
    for (Int_t idet=0; idet<540; idet++) {
      fCalChamberStatus->SetStatus(idet,AliTRDCalChamberStatus::kGood);
    }
    return;
  }

  // set all chambers to NoData
  for (Int_t idet=0; idet<540; idet++) {
    fCalChamberStatus->SetStatus(idet,AliTRDCalChamberStatus::kNoData);
  }

  // set status according to fHnSparseI 
  Int_t coord[4];
  for(Int_t bin = 0; bin < fHnSparseI->GetNbins(); bin++) {
    
    //Double_t content = fHnSparseI->GetBinContent(bin,coord);
    fHnSparseI->GetBinContent(bin,coord);
    // layer, stack, sector
    Int_t detector = AliTRDgeometry::GetDetector(coord[1]-1,coord[2]-1,coord[0]-1);
    //
    //printf("Number of entries for detector %d: %f\n",detector,content);
    //
    // Check which halfchamber side corresponds to the bin number (0=A, 1=B)
    // Change the status accordingly
    //
    if(coord[3]-1==0) { // HCS-A
      fCalChamberStatus->SetStatus(detector,AliTRDCalChamberStatus::kGood);
      fCalChamberStatus->UnsetStatusBit(detector,AliTRDCalChamberStatus::kNoDataHalfChamberSideA); // A has data
      //fCalChamberStatus->UnsetStatusBit(detector,AliTRDCalChamberStatus::kNoData); 
    }
    else { //HCS-B
      fCalChamberStatus->SetStatus(detector,AliTRDCalChamberStatus::kGood);
      fCalChamberStatus->UnsetStatusBit(detector,AliTRDCalChamberStatus::kNoDataHalfChamberSideB); // B has data
      //fCalChamberStatus->UnsetStatusBit(detector,AliTRDCalChamberStatus::kNoData); 
    }
  }

  // printf
  //for (Int_t idet=0; idet<540; idet++) {
  //  if(fCalChamberStatus->IsNoData(idet)) printf("No Data: chamber %d\n",idet);
  //}


}
//_____________________________________________________________________
void AliTRDCalibChamberStatus::CheckEORStatus(AliTRDCalDCSv2 *calDCS) /*FOLD00*/
{
  //
  //  Correct the AliTRDCalChamberStatus according to the AliTRDCalDCSv2
  //  Using globale state of the HalfChamberMerger (HCM)
  //
  for(Int_t det = 0; det < 540; det++) {
    AliTRDCalDCSFEEv2* calDCSFEEEOR = calDCS->GetCalDCSFEEObj(det);

    if(!calDCSFEEEOR) continue;
    
    // MCM Global State Machine State Definitions
    //  low_power =  0,
    //  test      =  1,
    //  wait_pre  =  3,
    //  preproc   =  7,
    //  zero_sp   =  8,
    //  full_rd   =  9,
    //  clear_st  = 11,
    //  wait_L1   = 12,
    //  tr_send   = 14,
    //  tr_proc   = 15 

    Int_t sm   = AliTRDgeometry::GetSector(det);
    Int_t lay  = AliTRDgeometry::GetLayer(det);
    Int_t stac = AliTRDgeometry::GetStack(det);
    
    Int_t stateA = 0;   // 0=bad, 1=good state
    Int_t stateB = 0;

    // loop over all mcm to define DCS-HCS
    for(Int_t ii = 0; ii < 8; ii++) { //ROB loop
      for(Int_t i = 0; i < 18; i++) { //MCM loop
	
	Int_t side = ii%2;  // 0=sideA, 1=sideB
	Int_t cstate = calDCSFEEEOR->GetMCMGlobalState(ii,i); //current mcm state
	
	if(cstate==3) {
	  switch(side) {
	  case 0: stateA=1; break;
	  case 1: stateB=1; break; 
	  }
	}
      }
    }

    //---------//
    //  Debug  //
    if(fDebugLevel > 0) {
      if( ((fCalChamberStatus->GetStatus(det) <= 1) && (stateA==0 || stateB==0)) || 
	  ((fCalChamberStatus->GetStatus(det) == 2) && (stateA==1 || stateB==1)) || 
	  ((fCalChamberStatus->GetStatus(det) == 3) && (stateA==1 || stateB==0)) ||
	  ((fCalChamberStatus->GetStatus(det) == 4) && (stateB==0 || stateB==1))  )
	{
	  //printf(" Different half chamber status in DCS and DATA!!\n");
	  Double_t val[4] = {sm,lay,stac,1};
	  fHnSparseDebug->Fill(&val[0]); 
	  
	  // Fill MCM status map
	  for(Int_t ii = 0; ii < 8; ii++) { //ROB loop
	    for(Int_t i = 0; i < 18; i++) { //MCM loop
	      Double_t valss[6] = {sm,lay,stac,ii,i
				   ,calDCSFEEEOR->GetMCMGlobalState(ii,i)};
	      fHnSparseMCM->Fill(&valss[0]);
	      
	    } 
	  }
	}
    }
    //---------//
    //  Debug  //

    //---------------------------------------
    // Change the status according to DCS
    //---------------------------------------
    Int_t StatusData = fCalChamberStatus->GetStatus(det);
    switch(StatusData) 
      {
      case 1: 
	if(stateA==0 && stateB==0) fCalChamberStatus->SetStatus(det,2); // completely masked from DCS
	if(stateA==1 && stateB==0) fCalChamberStatus->SetStatus(det,4); // Only B side masked from DCS
	if(stateA==0 && stateB==1) fCalChamberStatus->SetStatus(det,3); // Only A side masked from DCS
	if(stateA==1 && stateB==1) fCalChamberStatus->SetStatus(det,1);
	break;
      case 2: // completely masked from DATA
	if(stateA==0 && stateB==0) fCalChamberStatus->SetStatus(det,2); // completely masked from DCS
	break;
      case 3: // Only A side masked from DATA
	if(stateA==0 && stateB==0) fCalChamberStatus->SetStatus(det,2); // completely masked from DCS
	if(stateA==1 && stateB==0) fCalChamberStatus->SetStatus(det,2); // Only B side masked from DCS
	if(stateA==0 && stateB==1) fCalChamberStatus->SetStatus(det,3); // Only A side masked from DCS
	if(stateA==1 && stateB==1) fCalChamberStatus->SetStatus(det,3);
	break;
      case 4: // Only B side masked from DATA
	if(stateA==0 && stateB==0) fCalChamberStatus->SetStatus(det,2); // completely masked from DCS
	if(stateA==1 && stateB==0) fCalChamberStatus->SetStatus(det,4); // Only B side masked from DCS
	if(stateA==0 && stateB==1) fCalChamberStatus->SetStatus(det,2); // Only A side masked from DCS
	if(stateA==1 && stateB==1) fCalChamberStatus->SetStatus(det,4);
	break;
      }
    
  }

}

//_____________________________________________________________________________________
void AliTRDCalibChamberStatus::Add(AliTRDCalibChamberStatus *calibChamberStatus) /*FOLD00*/
{
    //
    //  Add the THnSparseI of this calibChamberStatus
    //

  fCounterEventNotEmpty += calibChamberStatus->GetNumberEventNotEmpty();

  THnSparseI *hnSparseI = calibChamberStatus->GetSparseI();
  if(!hnSparseI) return;

  if(!fHnSparseI) {
    fHnSparseI = (THnSparseI *) hnSparseI->Clone();
  }
  else {
    fHnSparseI->Add(hnSparseI);
  }
  

}
//_____________________________________________________________________
void AliTRDCalibChamberStatus::DumpToFile(const Char_t *filename, const Char_t *dir, Bool_t append) /*FOLD00*/
{
    //
    //  Write class to file
    //

    TString sDir(dir);
    TString option;

    if ( append )
	option = "update";
    else
        option = "recreate";

    TDirectory *backup = gDirectory;
    TFile f(filename,option.Data());
    f.cd();
    if ( !sDir.IsNull() ){
	f.mkdir(sDir.Data());
	f.cd(sDir);
    }
    this->Write();
    f.Close();

    if ( backup ) backup->cd();
}
//_____________________________________________________________________________
TH2D* AliTRDCalibChamberStatus::PlotSparseI(Int_t sm,Int_t side) 
{
  //
  // Plot number of entries for supermodule sm 
  // as a function of layer and stack
  //

  if(!fHnSparseI) return 0x0;
  
  fHnSparseI->GetAxis(0)->SetRange(sm+1,sm+1); 
  fHnSparseI->GetAxis(3)->SetRange(side+1,side+1);  
  TH2D *h2 = fHnSparseI->Projection(1,2);
 

  return h2;

}
//_____________________________________________________________________
TH2F *AliTRDCalibChamberStatus::MakeHisto2DSmPlEORStatus(AliTRDCalDCSv2 *calDCS, Int_t sm, Int_t pl) /*FOLD00*/
{
  //
  //  Plot globale state of the HalfChamberMerger (HCM)
  //
  AliTRDfeeParam *paramfee = AliTRDfeeParam::Instance();

  AliTRDgeometry *trdGeo = new AliTRDgeometry();
  AliTRDpadPlane *padPlane0 = trdGeo->GetPadPlane(pl,0);        // layer,stack
  Double_t row0    = padPlane0->GetRow0();
  Double_t col0    = padPlane0->GetCol0();

  char  name[1000];
  snprintf(name,1000,"%s DCS status sm %d pl %d",GetTitle(),sm,pl);
  TH2F * his = new TH2F( name, name, 88,-TMath::Abs(row0),TMath::Abs(row0)
                                   ,148,-TMath::Abs(col0),TMath::Abs(col0));


  // Where we begin
  Int_t offsetsmpl = 30*sm+pl;
  Int_t nstack = 5;
  Int_t ncols = 144;

  for (Int_t k = 0; k < nstack; k++){
    Int_t det = offsetsmpl+k*6;
    Int_t stac = AliTRDgeometry::GetStack(det);
    AliTRDCalDCSFEEv2* calDCSFEEEOR = calDCS->GetCalDCSFEEObj(det);
    if(!calDCSFEEEOR) { continue;}
    for (Int_t icol=0; icol<ncols; icol++){
      Int_t nrows = 16;
      if(stac==2) nrows = 12;
      for (Int_t irow=0; irow<nrows; irow++){
	Int_t binz     = 0;
	Int_t kb       = 5-1-k;
	Int_t krow     = nrows-1-irow;
	Int_t kcol     = ncols-1-icol;
	if(kb > 2) binz = 16*(kb-1)+12+krow+1+2*(kb+1);
	else binz = 16*kb+krow+1+2*(kb+1); 
	Int_t biny = kcol+1+2;
	// Take the value
	Int_t mcm = paramfee->GetMCMfromPad(irow,icol);
	Int_t rob = paramfee->GetROBfromPad(irow,icol);
	if(mcm < 0) AliWarning("Problem with mcm number");
	Int_t state = calDCSFEEEOR->GetMCMGlobalState(rob,TMath::Abs(mcm)); 
       	his->SetBinContent(binz,biny,state);
      }
    }
    for(Int_t icol = 1; icol < 147; icol++){
      for(Int_t l = 0; l < 2; l++){
	Int_t binz     = 0;
	Int_t kb       = 5-1-k;
	if(kb > 2) binz = 16*(kb-1)+12+1+2*(kb+1)-(l+1);
	else binz = 16*kb+1+2*(kb+1)-(l+1); 
	his->SetBinContent(binz,icol,16.0);
      }
    }
  }
  
  for(Int_t icol = 1; icol < 147; icol++){
    his->SetBinContent(88,icol,16.0);
    his->SetBinContent(87,icol,16.0);
  }
  for(Int_t irow = 1; irow < 89; irow++){
    his->SetBinContent(irow,1,16.0);
    his->SetBinContent(irow,2,16.0);
    his->SetBinContent(irow,147,16.0);
    his->SetBinContent(irow,148,16.0);
  }

  his->SetXTitle("z (cm)");
  his->SetYTitle("y (cm)");
  his->SetMaximum(12);
  his->SetMinimum(0.0);
  his->SetStats(0);

  return his;

}
//_____________________________________________________________________________
TCanvas* AliTRDCalibChamberStatus::PlotHistos2DSmEORStatus(AliTRDCalDCSv2 *calDCS, Int_t sm, const Char_t *name)
{
  //
  // Make 2D graph
  //

  gStyle->SetPalette(1);
  fC1 = new TCanvas(name,name,50,50,600,800);
  fC1->Divide(3,2);
  fC1->cd(1);
  MakeHisto2DSmPlEORStatus(calDCS,sm,0)->Draw("colz");
  fC1->cd(2);
  MakeHisto2DSmPlEORStatus(calDCS,sm,1)->Draw("colz");
  fC1->cd(3);
  MakeHisto2DSmPlEORStatus(calDCS,sm,2)->Draw("colz");
  fC1->cd(4);
  MakeHisto2DSmPlEORStatus(calDCS,sm,3)->Draw("colz");
  fC1->cd(5);
  MakeHisto2DSmPlEORStatus(calDCS,sm,4)->Draw("colz");
  fC1->cd(6);
  MakeHisto2DSmPlEORStatus(calDCS,sm,5)->Draw("colz");

  return fC1;

}



