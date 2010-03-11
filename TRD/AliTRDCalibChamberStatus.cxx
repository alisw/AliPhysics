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
// Check with the AliTRDCalDCSFEE info                                    //
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

//AliRoot includes
#include "AliRawReader.h"

//header file
#include "AliLog.h"
#include "AliTRDCalibChamberStatus.h"
#include "AliTRDgeometry.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDSignalIndex.h"
#include "AliTRDrawFastStream.h"
#include "./Cal/AliTRDCalChamberStatus.h"
#include "./Cal/AliTRDCalDCS.h"
#include "./Cal/AliTRDCalDCSFEE.h"

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
    thnDimEvtt[4] = 16;
    thnDimEvtt[5] = 16;
    //arrays for lower bounds :
    Double_t* binEdgesEvtt[6];
    for(Int_t ivar = 0; ivar < 6; ivar++)
      binEdgesEvtt[ivar] = new Double_t[thnDimEvtt[ivar] + 1];
    //values for bin lower bounds
    for(Int_t i=0; i<=thnDimEvtt[0]; i++) binEdgesEvtt[0][i]= 0.0  + (18.0)/thnDimEvtt[0]*(Double_t)i;
    for(Int_t i=0; i<=thnDimEvtt[1]; i++) binEdgesEvtt[1][i]= 0.0  + (6.0)/thnDimEvtt[1]*(Double_t)i;
    for(Int_t i=0; i<=thnDimEvtt[2]; i++) binEdgesEvtt[2][i]= 0.0  + (5.0)/thnDimEvtt[2]*(Double_t)i;
    for(Int_t i=0; i<=thnDimEvtt[3]; i++) binEdgesEvtt[3][i]= 0.0  + (8.0)/thnDimEvtt[3]*(Double_t)i;
    for(Int_t i=0; i<=thnDimEvtt[4]; i++) binEdgesEvtt[4][i]= 0.0  + (16.0)/thnDimEvtt[4]*(Double_t)i;
    for(Int_t i=0; i<=thnDimEvtt[5]; i++) binEdgesEvtt[5][i]= 0.0  + (16.0)/thnDimEvtt[5]*(Double_t)i;
    
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
void AliTRDCalibChamberStatus::ProcessEvent(AliRawReader * rawReader, Int_t nevents_physics)
{
  //
  // Event Processing loop 
  //
  //
  
  Bool_t notEmpty = kFALSE;
    
  AliTRDrawFastStream *rawStream = new AliTRDrawFastStream(rawReader);
  rawStream->SetSharedPadReadout(kFALSE);

  AliTRDdigitsManager *digitsManager = new AliTRDdigitsManager(kTRUE);
  digitsManager->CreateArrays();
  
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

  if(digitsManager) delete digitsManager;
  if(rawStream) delete rawStream;
   
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
void AliTRDCalibChamberStatus::AnalyseHisto() /*FOLD00*/
{
    //
    //  Create the AliTRDCalChamberStatus according to the fHnSparseI
    //

  if(fCalChamberStatus) delete fCalChamberStatus;
  fCalChamberStatus = new AliTRDCalChamberStatus();

  // Check if enough events to say something
  if(fCounterEventNotEmpty < 30) {
    // Say all installed
    for (Int_t ism=0; ism<18; ism++) {
      for (Int_t ipl=0; ipl<6; ipl++) {
	for (Int_t istack=0; istack<5; istack++) {
	  // layer, stack, sector
	  Int_t det = AliTRDgeometry::GetDetector(ipl,istack,ism);
	  fCalChamberStatus->SetStatus(det,1);
	}
      }
    }
    return;
  }

  // Mask out all chambers
  for (Int_t ism=0; ism<18; ism++) {
    for (Int_t ipl=0; ipl<6; ipl++) {
      for (Int_t istack=0; istack<5; istack++) {
	// layer, stack, sector
	Int_t det = AliTRDgeometry::GetDetector(ipl,istack,ism);
	fCalChamberStatus->SetStatus(det,2);
      }
    }
  }

  // Unmask good chambers 
  Int_t coord[4];
  for(Int_t bin = 0; bin < fHnSparseI->GetNbins(); bin++) {
    
    fHnSparseI->GetBinContent(bin,coord);
    // layer, stack, sector
    Int_t detector = AliTRDgeometry::GetDetector(coord[1]-1,coord[2]-1,coord[0]-1);

    //
    // Check which halfchamber side corresponds to the bin number (0=A, 1=B)
    // Change the status accordingly
    //

    switch(fCalChamberStatus->GetStatus(detector)) 
      {    
      case 1: break;  // no changes
      case 2: 
	if(coord[3]-1==0) {
	  fCalChamberStatus->SetStatus(detector,4); break;      // only SideB is masked
	}
	else {
	  fCalChamberStatus->SetStatus(detector,3); break;      // only SideA is masked
	}
      case 3:  fCalChamberStatus->SetStatus(detector,1); break;  // unmask SideA
      case 4:  fCalChamberStatus->SetStatus(detector,1); break;  // unmask SideB
      }
  }


}
//_____________________________________________________________________
void AliTRDCalibChamberStatus::CheckEORStatus(AliTRDCalDCS *calDCS) /*FOLD00*/
{
  //
  //  Correct the AliTRDCalChamberStatus according to the AliTRDCalDCS
  //  Using globale state of the HalfChamberMerger (HCM)
  //
  
  for(Int_t det = 0; det < 540; det++) {
    AliTRDCalDCSFEE* calDCSFEEEOR = calDCS->GetCalDCSFEEObj(det);
    if(!calDCSFEEEOR) { continue;}
    
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
    
    Int_t stateA = calDCSFEEEOR->GetMCMGlobalState(4,17); // HCM Side A
    Int_t stateB = calDCSFEEEOR->GetMCMGlobalState(5,17); // HCM Side B
    Int_t rphi = -1;

    //printf("DCS: stateA %d \t stateB %d \n",stateA,stateB);
    if(stateA!=3 && stateA!=9) rphi = 1;
    Double_t vals[4] = {sm,lay,stac,rphi};
    if(rphi!=-1) fHnSparseHCM->Fill(&vals[0]);
    
    if(stateB!=3 && stateB!=9) rphi = 2;
    vals[3] = rphi;
    if(rphi!=-1) fHnSparseHCM->Fill(&vals[0]);
    
    //---------//
    //  Debug  //
    if(fDebugLevel > 0) {
      if( (fCalChamberStatus->GetStatus(det) <= 1) && (stateA!=3 && stateA!=9) || 
	  (fCalChamberStatus->GetStatus(det) <= 1) && (stateB!=3 && stateB!=9) || 
	  (fCalChamberStatus->GetStatus(det) >= 2) && (stateA==3 || stateA==9) || 
	  (fCalChamberStatus->GetStatus(det) >= 2) && (stateB==3 || stateB==9)  )
	{
	  //printf(" Different half chamber status in DCS and DATA!!\n");
	  Double_t val[4] = {sm,lay,stac,1};
	  fHnSparseDebug->Fill(&val[0]); 
	  
	  if(rphi!=-1) {  // error in DCS information
	    // Fill MCM status map
	    for(Int_t ii = 0; ii < 8; ii++) { //ROB loop
	      for(Int_t i = 0; i < 16; i++) { //MCM loop
		Double_t valss[6] = {sm,lay,stac,ii,i
				     ,calDCSFEEEOR->GetMCMGlobalState(ii,i)};
		fHnSparseMCM->Fill(&valss[0]);
	      }
	    } 
	  }
	}
    }
    //---------//
    //  Debug  //
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


