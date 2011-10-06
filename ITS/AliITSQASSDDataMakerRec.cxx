/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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
/* $Id$    */
//  *************************************************************
//  Checks the quality assurance 
//  by comparing with reference data
//  contained in a DB
//  -------------------------------------------------------------
//  W. Ferrarese + P. Cerello Feb 2008
//  INFN Torino

// --- ROOT system ---
#include <TH2F.h>
#include <TTree.h>
#include <TMath.h>
#include <TString.h>
#include <TSystem.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliITSQASSDDataMakerRec.h"
#include "AliQADataMakerRec.h"
#include "AliLog.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliRawReader.h"
#include "AliRawReaderRoot.h"
#include "AliITSRawStreamSSD.h"
#include "AliITSgeomTGeo.h"
#include "AliRawEventHeaderBase.h"
#include "AliITSRecPoint.h"
#include "AliITSRecPointContainer.h"
#include "AliITSdigitSSD.h"
#include "AliITSBadChannelsSSDv2.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"

ClassImp(AliITSQASSDDataMakerRec)

AliITSQASSDDataMakerRec::AliITSQASSDDataMakerRec(AliITSQADataMakerRec *aliITSQADataMakerRec, Bool_t kMode, Int_t ldc) :
TObject(),
  fAliITSQADataMakerRec(aliITSQADataMakerRec),
  fkOnline(kMode),
  fLDC(ldc),
  fSSDRawsOffset(0), fSSDRawsDAOffset(0),
  fSSDRawsCommonLevelOffset(0),
  fSSDhRawsTask(0),
  fSSDhDigitsTask(0),
  fSSDhRecPointsTask(0),
  fGenRawsOffset(0),
  fGenDigitsOffset(0),
  fGenRecPointsOffset(0),
  fCDBManager(0) {
  // Default constructor   
  //initilize the raw signal vs strip number histograms
  fGenRawsOffset = new Int_t[AliRecoParam::kNSpecies];
  fGenRecPointsOffset = new Int_t[AliRecoParam::kNSpecies];
  fGenDigitsOffset = new Int_t[AliRecoParam::kNSpecies];
  for(Int_t i=0; i<AliRecoParam::kNSpecies; i++) {
    fGenRawsOffset[i] = 0;
    fGenRecPointsOffset[i] = 0;
    fGenDigitsOffset[i]=0;
  }
  if(fkOnline) {
    fCDBManager = AliCDBManager::Instance();
    //fCDBManager->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    fCDBManager->SetDefaultStorage(gSystem->Getenv("AMORE_CDB_URI"));

    Int_t runNumber = 1;
    if(!gSystem->Getenv("DATE_RUN_NUMBER")) {
      AliWarning("DATE_RUN_NUMBER not defined!!!\n"); }
    else {
      runNumber = atoi(gSystem->Getenv("DATE_RUN_NUMBER"));}
    fCDBManager->SetRun(runNumber);
    AliCDBEntry *geomGRP = fCDBManager->Get("GRP/Geometry/Data");
    if(!geomGRP) AliWarning("GRP geometry not found!!!\n");
    
    Int_t gLayer = 0,gLadder = 0, gModule = 0;
    Int_t gHistCounterRawSignal = 0, gHistCounterCM = 0;
    TString gTitle; 

    for(Int_t iModule = 500; iModule < fgkSSDMODULES + 500; iModule++) {
      AliITSgeomTGeo::GetModuleId(iModule,gLayer,gLadder,gModule);
      gTitle = "SSD_RawSignal_Layer"; gTitle += gLayer;
      gTitle += "_Ladder"; gTitle += gLadder;
      gTitle += "_Module"; gTitle += gModule;
      fHistSSDRawSignalModule[gHistCounterRawSignal] = new TH1F(gTitle.Data(),gTitle.Data(),
								2*fgkNumberOfPSideStrips,0,2*fgkNumberOfPSideStrips);
      gHistCounterRawSignal += 1;
      
      for(Int_t iStrip = 0; iStrip < 2*fgkNumberOfPSideStrips; iStrip++)
	fOccupancyMatrix[iModule-500][iStrip] = 0; 

      //CM histograms
      gTitle = "SSD_CM_PSide_Layer"; gTitle += gLayer;
      gTitle += "_Ladder"; gTitle += gLadder;
      gTitle += "_Module"; gTitle += gModule;
      fHistSSDCMModule[gHistCounterCM] = new TH1F(gTitle.Data(),gTitle.Data(),
						  100,-50.,50.);
      fHistSSDCMModule[gHistCounterCM]->GetXaxis()->SetTitle("CM");
      gHistCounterCM += 1;
    }//module loop

    for(Int_t iModule = 500; iModule < fgkSSDMODULES + 500; iModule++) {
      AliITSgeomTGeo::GetModuleId(iModule,gLayer,gLadder,gModule);
      gTitle = "SSD_CM_NSide_Layer"; gTitle += gLayer;
      gTitle += "_Ladder"; gTitle += gLadder;
      gTitle += "_Module"; gTitle += gModule;
      fHistSSDCMModule[gHistCounterCM] = new TH1F(gTitle.Data(),gTitle.Data(),
						  100,-50.,50.);
      fHistSSDCMModule[gHistCounterCM]->GetXaxis()->SetTitle("CM");
      gHistCounterCM += 1;

    }
  }//online flag
  else {
    fCDBManager = NULL;
    for(Int_t iModule = 0; iModule < fgkSSDMODULES; iModule++) {
      fHistSSDCMModule[iModule] = NULL;
      fHistSSDCMModule[fgkSSDMODULES+iModule] = NULL;
      fHistSSDRawSignalModule[iModule] = NULL;
    }
  }
}
/*
//____________________________________________________________________________ 
AliITSQASSDDataMakerRec::AliITSQASSDDataMakerRec(const AliITSQASSDDataMakerRec& qadm) :
TObject(),
fAliITSQADataMakerRec(qadm.fAliITSQADataMakerRec),
fkOnline(qadm.fkOnline),
fLDC(qadm.fLDC),
fSSDRawsOffset(qadm.fSSDRawsOffset), fSSDRawsDAOffset(qadm.fSSDRawsDAOffset),
fSSDRawsCommonLevelOffset(qadm.fSSDRawsCommonLevelOffset),
fSSDhRawsTask(qadm.fSSDhRawsTask),
fSSDhDigitsTask(qadm.fSSDhDigitsTask),
fSSDhRecPointsTask(qadm.fSSDhRecPointsTask),
fGenRawsOffset(qadm.fGenRawsOffset),
fGenDigitsOffset(qadm.fGenDigitsOffset),
fGenRecPointsOffset(qadm.fGenRecPointsOffset),
fCDBManager(qadm.fCDBManager) {
//copy ctor 
fAliITSQADataMakerRec->SetName((const char*)qadm.fAliITSQADataMakerRec->GetName()) ; 
fAliITSQADataMakerRec->SetTitle((const char*)qadm.fAliITSQADataMakerRec->GetTitle());
}

//__________________________________________________________________
AliITSQASSDDataMakerRec& AliITSQASSDDataMakerRec::operator = (const AliITSQASSDDataMakerRec& qac )
{
// Equal operator.
this->~AliITSQASSDDataMakerRec();
new(this) AliITSQASSDDataMakerRec(qac);
return *this;
}
*/
//__________________________________________________________________
AliITSQASSDDataMakerRec::~AliITSQASSDDataMakerRec() {
  // destructor
  if(fkOnline) {
    for(Int_t iModule = 0; iModule < fgkSSDMODULES; iModule++) {
      if(fHistSSDRawSignalModule[iModule]) delete fHistSSDRawSignalModule[iModule];
      if(fHistSSDCMModule[iModule]) delete fHistSSDCMModule[iModule];
      if(fHistSSDCMModule[fgkSSDMODULES+iModule]) delete fHistSSDCMModule[fgkSSDMODULES+iModule];
    }
    
    if(fCDBManager) delete fCDBManager;
  }
}

//____________________________________________________________________________ 
void AliITSQASSDDataMakerRec::StartOfDetectorCycle()
{

  if(fAliITSQADataMakerRec->ListExists(AliQAv1::kRAWS)==kFALSE)return;

  //if ( fAliITSQADataMakerRec->GetRawsData(0) == NULL ) // Raws not defined
  //return ;
 
  //for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
  //if (!AliQAv1::Instance()->IsEventSpecieSet(specie)) continue;
  //cout << "StartOfDetectorCycle: Event specie " << specie << " is set" << endl;

  //Detector specific actions at start of cycle
  AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQADM::Start of SSD Cycle\n");
  //Int_t specie = fAliITSQADataMakerRec->GetEventSpecie();

  //}//event specie loop

}

//____________________________________________________________________________ 
void AliITSQASSDDataMakerRec::ResetRawsMonitoredObjects() 
{
  //Resetting the raw data monitored objects
  //Data size per DDL
  //for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
  //if (!AliQAv1::Instance()->IsEventSpecieSet(specie)) continue;
  /*((TH1F *)(fAliITSQADataMakerRec->ResetRawsData(fGenRawsOffset[specie]+4)));
  //Data size per LDC
  ((TH1F *)(fAliITSQADataMakerRec->ResetRawsData(fGenRawsOffset[specie]+22)));*/
  Int_t specie = fAliITSQADataMakerRec->GetEventSpecie();
  //cout << "(AliITSQASSDDataMakerRec::ResetRawsMonitoredObjects): Event specie " << specie << " is set" << endl;
  //online part
  if(fkOnline) {
    for(Int_t iModule = 500; iModule < fgkSSDMODULES + 500; iModule++) {
      for(Int_t iStrip = 0; iStrip < 2*fgkNumberOfPSideStrips; iStrip++)
	fOccupancyMatrix[iModule-500][iStrip] = 0;
    }//module loop
    
    //for(Int_t iSSDOffset = 0; iSSDOffset < fSSDRawsCommonLevelOffset; iSSDOffset++)
    //(fAliITSQADataMakerRec->ResetRawsData(fGenRawsOffset[specie]+iSSDOffset));
    //for(Int_t iSSDOffset = fSSDRawsCommonLevelOffset; iSSDOffset < fSSDRawsOffset; iSSDOffset++)
    //(fAliITSQADataMakerRec->ResetRawsData(fGenRawsOffset[specie]+iSSDOffset));

    Int_t gHistPositionOccupancyPerLadder = 0;
    Int_t gLayer = 0, gLadder = 0, gModule = 0;
    for(Int_t iModule = 0; iModule < fgkSSDMODULES; iModule++) {
      AliITSgeomTGeo::GetModuleId(iModule+500,gLayer,gLadder,gModule);
      
      gHistPositionOccupancyPerLadder = (gLayer == 5) ? 2*(gLadder - 1) : 2*(gLadder - 1 + fgkSSDLADDERSLAYER5);
      
      //P-SIDE OCCUPANCY
      int offs = fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+fgkSSDMODULES;
      fAliITSQADataMakerRec->ResetRawsData(offs+gHistPositionOccupancyPerLadder);
      //N-SIDE OCCUPANCY
      fAliITSQADataMakerRec->ResetRawsData(offs+gHistPositionOccupancyPerLadder+1);
      //
      fAliITSQADataMakerRec->ResetRawsData(offs+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6);
      fAliITSQADataMakerRec->ResetRawsData(offs+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+1);
    }//module loop
  }//online flag
  //
}

//____________________________________________________________________________ 
void AliITSQASSDDataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray** /*list*/)
{
  // finalize ssd cycle
  //
  Int_t specie = fAliITSQADataMakerRec->GetEventSpecie();
  //cout << "(AliITSQASSDDataMakerRec::EndOfDetectorCycle): Event specie " << specie << " is set" << endl;
  
  // launch the QA checking
  AliDebug(AliQAv1::GetQADebugLevel(),"AliITSDM instantiates checker with Run(AliQAv1::kITS, task, list)\n"); 
  AliDebug(AliQAv1::GetQADebugLevel(), Form("Offset: %d\n",fGenRawsOffset[specie]));
  //Printf("Offset: %d\n",fGenRawsOffset[specie]);

  if(task == AliQAv1::kRAWS) {
    //
    for (int trCl=-1;trCl<fAliITSQADataMakerRec->GetNTrigClasses();trCl++) { // RS Loop over all trigger classes (+ global counter, -1)
      //
      TObjArray &harr = *fAliITSQADataMakerRec->GetRawsDataOfTrigClass(trCl);
      int offs = fGenRawsOffset[specie];
      int nSSDEventPerCycle = fAliITSQADataMakerRec->GetEvCountCycleRaws(trCl);
      //    
      //Data size per DDL
      for(Int_t i = 0; i < fgkNumOfDDLs; i++) {
	if (!(harr[offs+5+i] && harr[offs+4+i])) continue;
	Double_t gSizePerDDL = ((TH1*)harr[offs+5+i])->GetMean();
	//cout<<"DDL: "<<i+512<<" - Size: "<<gSizePerDDL<<" - Mean: "<<gSizePerDDL<<endl;
	//cout<<"Entries: "<<((TH1*)harr[offs+5+i])->GetEntries()<<endl;
	((TH1*)harr[offs+4])->SetBinContent(i+1,gSizePerDDL);
	//cout<<"After filling DDL: "<<i+512<<" - Size: "<< ((TH1F*)harr[offs+4+i])->GetBinContent(i+1)<<endl;
      }
      //
      //Data size per LDC
      for(Int_t i = 0; i < fgkNumOfLDCs; i++) {
	if ( !(harr[offs+23+i]&&harr[offs+22+i])) continue;
	Double_t gSizePerLDC = ((TH1*)harr[offs+23+i])->GetMean();
	((TH1*)harr[offs+22])->SetBinContent(i+1,gSizePerLDC);
	//cout<<"LDC: "<<i+170<<" - Size: "<<gSizePerLDC<<" - Mean: "<<" - Size: "<<((TH1*)harr[offs+23+i])->GetMean()<<endl;
      }
      //
      //if (harr[offs+4])  cout<<"Data size/ DDL entries: "<<((TH1*)harr[offs+4 ])->GetEntries()<< " mean: " << ((TH1*)harr[offs+4])->GetMean()<<endl;   
      //if (harr[offs+22]) cout<<"Data size/ LDC entries: "<<((TH1*)harr[offs+22])->GetEntries()<< " mean: " << ((TH1*)harr[offs+22])->GetMean()<<endl;
      //
      //online part
      if(fkOnline) {
	//Output of the DA
	MonitorOCDBObjects(trCl);
	//Monitor common mode values
	MonitorCMValues(trCl);
	//
	Int_t gHistPositionOccupancyPerModule = 0;
	Int_t gLayer = 0, gLadder = 0, gModule = 0;
	//occupancy per module
	for(Int_t iModule = 0; iModule < fgkSSDMODULES; iModule++) {
	  AliITSgeomTGeo::GetModuleId(iModule+500,gLayer,gLadder,gModule);	  
	  gHistPositionOccupancyPerModule = (gLayer == 5) ? ((gLadder - 1)*fgkSSDMODULESPERLADDERLAYER5 + gModule - 1) : ((gLadder - 1)*fgkSSDMODULESPERLADDERLAYER6 + gModule + fgkSSDMODULESLAYER5 - 1);
	  TH1* htmp = (TH1*)harr[offs+fSSDRawsCommonLevelOffset+gHistPositionOccupancyPerModule];
	  if (htmp) {
	    for(Int_t iBins = 1; iBins < fHistSSDRawSignalModule[iModule]->GetXaxis()->GetNbins(); iBins++) htmp->SetBinContent(iBins,fOccupancyMatrix[iModule][iBins-1]);
	    if(nSSDEventPerCycle != 0) htmp->Scale(100./nSSDEventPerCycle);
	  }//module loop
	}
	//
	//occupancy per ladder
	Int_t gHistPositionOccupancyPerLadder = 0;
	Int_t lLadderLocationY = 0;
	Double_t occupancy = 0.0, occupancyThreshold = 0.0, occupancyAverage = 0.0;
	for(Int_t iModule = 0; iModule < fgkSSDMODULES; iModule++) {
	  AliITSgeomTGeo::GetModuleId(iModule+500,gLayer,gLadder,gModule);
	  //
	  gHistPositionOccupancyPerModule = (gLayer == 5) ? ((gLadder - 1)*fgkSSDMODULESPERLADDERLAYER5 + gModule - 1) : ((gLadder - 1)*fgkSSDMODULESPERLADDERLAYER6 + gModule + fgkSSDMODULESLAYER5 - 1);
	  gHistPositionOccupancyPerLadder = (gLayer == 5) ? 2*(gLadder - 1) : 2*(gLadder - 1 + fgkSSDLADDERSLAYER5);
	  //
	  TH1* htmpo = (TH1*)harr[offs+fSSDRawsCommonLevelOffset+gHistPositionOccupancyPerModule];
	  TH1* h1t = 0;
	  TH2* h2t = 0;
	  if (htmpo) {
	    //P-SIDE OCCUPANCY
	    occupancy          = GetOccupancyModule(htmpo,0,0,0);
	    occupancyThreshold = GetOccupancyModule(htmpo,0,1,3);
	    occupancyAverage   = GetOccupancyModule(htmpo,0,2,0);
	    if ( (h1t=(TH1*)harr[offs+fSSDRawsCommonLevelOffset+fgkSSDMODULES+gHistPositionOccupancyPerLadder]) ) h1t->SetBinContent(gModule,occupancy);
	    lLadderLocationY = 3*gLadder; // sideP=1 sideN=0 
	    if(gLayer == 5) {
	      //occupancy per module - no threshold
	      if ( (h2t=(TH2*)harr[offs+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6])   ) h2t->SetBinContent(gModule,lLadderLocationY,occupancy);
	      //occupancy per module - threshold @ 3%
	      if ( (h2t=(TH2*)harr[offs+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+2]) ) h2t->SetBinContent(gModule,lLadderLocationY,occupancyThreshold);
	      //average occupancy per module
	      if ( (h2t=(TH2*)harr[offs+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+4]) ) h2t->SetBinContent(gModule,lLadderLocationY,occupancyAverage);
	    }
	    else if(gLayer == 6) {
	      //occupancy per module - no threshold
	      if ( (h2t=(TH2*)harr[offs+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+1]) ) h2t->SetBinContent(gModule,lLadderLocationY,occupancy);
	      //occupancy per module - threshold @ 3%
	      if ( (h2t=(TH2*)harr[offs+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+3]) ) h2t->SetBinContent(gModule,lLadderLocationY,occupancyThreshold);
	      //average occupancy per module
	      if ( (h2t=(TH2*)harr[offs+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+5]) ) h2t->SetBinContent(gModule,lLadderLocationY,occupancyAverage);
	    }
	    //
	    //N-SIDE OCCUPANCY
	    //
	    occupancy          = GetOccupancyModule(htmpo,1,0,0);   
	    occupancyThreshold = GetOccupancyModule(htmpo,1,1,3);   
	    occupancyAverage   = GetOccupancyModule(htmpo,1,2,0);   
	    if ( (h1t=(TH1*)harr[offs+fSSDRawsCommonLevelOffset+fgkSSDMODULES+gHistPositionOccupancyPerLadder+1]) ) h1t->SetBinContent(gModule,occupancy);
	    if(gLayer == 5) {
	      //occupancy per module - no threshold
	      if ( (h2t=(TH2*)harr[offs+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6])   ) h2t->SetBinContent(gModule,lLadderLocationY-1,occupancy);
	      //occupancy per module - threshold @ 3%
	      if ( (h2t=(TH2*)harr[offs+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+2]) ) h2t->SetBinContent(gModule,lLadderLocationY-1,occupancyThreshold);
	      //average occupancy per module
	      if ( (h2t=(TH2*)harr[offs+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+4]) ) h2t->SetBinContent(gModule,lLadderLocationY-1,occupancyAverage);
	    }
	    else if(gLayer == 6) {
	      //occupancy per module - no threshold
	      if ( (h2t=(TH2*)harr[offs+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+1]) ) h2t->SetBinContent(gModule,lLadderLocationY-1,occupancy);
	      //occupancy per module - threshold @ 3%
	      if ( (h2t=(TH2*)harr[offs+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+3]) ) h2t->SetBinContent(gModule,lLadderLocationY-1,occupancyThreshold);
	      //average occupancy per module
	      if ( (h2t=(TH2*)harr[offs+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+5]) ) h2t->SetBinContent(gModule,lLadderLocationY-1,occupancyAverage);
	    }
	  } // htmpo
	}//module loop
      }//online flag for SSD
      //
      //TH2* h2tmp = (TH2*)harr[offs+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6];
      //if (h2tmp) AliInfo(Form("Entries 2d occupancy no thres- lay 5: %d"),h2tmp->GetEntries());
      //h2tmp = (TH2*)harr[offs+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+3];
      //cout<<"entries 2d occupancy thres- lay 6: "<<h2tmp->GetEntries()<< " mean: " << h2tmp->GetMean() << endl; //Somehow the other occupancy maps do give nonzero values for GetMean() here
      //
      // TH1* h1tmp = (TH1*)harr[fGenRawsOffset[1]+4];
      // if (h1tmp) cout<<"Data size/ DDL entries: "<<h1tmp->GetEntries()<< " mean: " << h1tmp->GetMean()<<endl;   
      //Reset of the raws
    } //  RS Loop over all trigger classes (+ global counter, -1)
    //AliQAChecker::Instance()->Run( AliQAv1::kITS , task, list);
    //
    ResetRawsMonitoredObjects();
  }//raw data end of cycle
  //
      
  // 
}

//____________________________________________________________________________ 
Int_t AliITSQASSDDataMakerRec::InitRaws() {  
  // Initialization for RAW data - SSD -
  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  Int_t rv = 0 ; 
  fSSDRawsOffset = 0;

  //for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
  //if (!AliQAv1::Instance()->IsEventSpecieSet(specie)) continue;
  Int_t specie = fAliITSQADataMakerRec->GetEventSpecie();
  //cout << "(AliITSQASSDDataMakerRec::InitRaws): Event specie " << specie << " is set" << endl;
  //cout << "(AliITSQASSDDataMakerRec::InitRaws): Offset " << offsRw << endl;
  int offsRw = fGenRawsOffset[specie];
    
    
  if(fkOnline) {
    AliDebug(AliQAv1::GetQADebugLevel(), "Book Online Histograms for SSD\n");
  }
  else {
    AliDebug(AliQAv1::GetQADebugLevel(), "Book Offline Histograms for SSD\n ");
  }
  AliDebug(AliQAv1::GetQADebugLevel(), Form("Number of histograms (SPD+SDD): %d\n",offsRw));
  
  TString gTitle;
  TString gName;
  //book online-offline QA histos
  TH1F *fHistSSDEventType = new TH1F("fHistSSDEventType",
				     "SSD Event Type;Event type;Events",
				     31,-1,30);
  fHistSSDEventType->SetStats(kFALSE);
  rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDEventType, 
					   offsRw+fSSDRawsOffset, expert, !image, !saveCorr);
  fSSDRawsOffset += 1;
  //cout<<"(AliITSQASSDDataMakerRec::InitRaws): "<<offsRw+fSSDRawsOffset-1<<" - Name: "<<(fAliITSQADataMakerRec->GetRawsData(offsRw))->GetName()<<endl;
  TH1F *fHistSSDDataSize = new TH1F("fHistSSDDataSize",
				    "SSD Data Size;(SSD data size) [KB];Events",
				    1000,0,500);
  rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDDataSize, 
					   offsRw+fSSDRawsOffset, !expert, image, !saveCorr);
  fSSDRawsOffset += 1;
  TH1F *fHistSSDDataSizePercentage = new TH1F("fHistSSDDataSizePercentage",
					      "SSD Data Size Percentage;SSD data size [%];Events",
					      1000,0,100);
  rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDDataSizePercentage, 
					   offsRw+fSSDRawsOffset, expert, !image, !saveCorr);
  fSSDRawsOffset += 1;
  TH1F *fHistSSDDDLId = new TH1F("fHistSSDDDLId",
				 "SSD DDL Id;DDL id;Events",16,511.5,527.5);
  fHistSSDDDLId->SetStats(kFALSE);
  rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDDDLId, 
					   offsRw+fSSDRawsOffset, expert, !image, !saveCorr);
  fSSDRawsOffset += 1;
  TH1F *fHistSSDDataSizePerDDL = new TH1F("fHistSSDDataSizePerDDL",
					  "SSD Data Size Per DDL;DDL id;<SSD data size> [KB]",
					  16,511.5,527.5);
  fHistSSDDataSizePerDDL->SetStats(kFALSE);
  rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDDataSizePerDDL, 
					   offsRw+fSSDRawsOffset, !expert, image, !saveCorr);

  fSSDRawsOffset += 1;
  TH1F *fHistSSDDataSizeDDL[fgkNumOfDDLs];
  for(Int_t i = 1; i < fgkNumOfDDLs+1; i++) {
    gName = Form("fHistSSDDataSizeDDL%d", i+511) ;
    gTitle = Form("SSD Data Size DDL %d", i+511) ;
    fHistSSDDataSizeDDL[i-1] = new TH1F(gName.Data(),
					Form("%s;(SSD data size) [KB];Events", gTitle.Data()),
					100,0,50);
    rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDDataSizeDDL[i-1], 
					     offsRw+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1;
  }
  
  TH1F *fHistSSDLDCId = new TH1F("fHistSSDLDCId","SSD LDC Id;LDC id;Events",8,169.5,177.5);
  fHistSSDLDCId->SetStats(kFALSE);
  rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDLDCId, 
					   offsRw+fSSDRawsOffset, expert, !image, !saveCorr);
  fSSDRawsOffset += 1;
  TH1F *fHistSSDDataSizePerLDC = new TH1F("fHistSSDDataSizePerLDC",
					  "SSD Data Size Per LDC;LDC id;<SSD data size> [KB]",
					  8,169.5,177.5);
  fHistSSDDataSizePerLDC->SetStats(kFALSE);
  rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDDataSizePerLDC, 
					   offsRw+fSSDRawsOffset, expert, !image, !saveCorr);  fSSDRawsOffset += 1;
  TH1F *fHistSSDDataSizeLDC[fgkNumOfLDCs];
  for(Int_t i = 1; i < fgkNumOfLDCs+1; i++) {
    gName = "fHistSSDDataSizeLDC"; 
    if(i == 1) gName += "170";
    if(i == 2) gName += "171";
    if(i == 3) gName += "172";
    if(i == 4) gName += "173";
    if(i == 5) gName += "174";
    if(i == 6) gName += "175";
    if(i == 7) gName += "176";
    if(i == 8) gName += "177";
    
    gTitle = "SSD Data Size LDC "; gTitle += gName.Data();
    fHistSSDDataSizeLDC[i-1] = new TH1F(gName.Data(),
					Form("%s;SSD data size [KB];Events", gTitle.Data()),
					1000,0,100);
    rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDDataSizeLDC[i-1], 
					     offsRw+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1;
  }
  fSSDRawsCommonLevelOffset = fSSDRawsOffset;
  
  if(fkOnline) {
    Int_t gLayer = 0, gLadder = 0, gModule = 0;
    //occupancy per SSD module
    TH1F *fHistSSDOccupancyModule[fgkSSDMODULES]; 
    for(Int_t i = 500; i < fgkSSDMODULES + 500; i++) {
      AliITSgeomTGeo::GetModuleId(i,gLayer,gLadder,gModule);
      gName = "fHistSSD_Occupancy_Layer";
      gTitle = "SSD Occupancy Layer";
      if(gLayer == 5) {
	gName += gLayer; gName += "_Ladder"; 
	gName += 499+gLadder;
	gTitle += gLayer; gTitle += "_Ladder"; 
	gTitle += 499+gLadder;
      }
      if(gLayer == 6) {
	gName += gLayer; gName += "_Ladder"; 
	gName += 599+gLadder;
	gTitle += gLayer; gTitle += "_Ladder"; 
	gTitle += 599+gLadder;
      }
      gName += "_Module"; gName += gModule; 
      gTitle += "_Module"; gTitle += gModule; 
      
      fHistSSDOccupancyModule[i-500] = new TH1F(gName.Data(),Form("%s;N_{strip};Occupancy [%%]", gTitle.Data()),
						2*fgkNumberOfPSideStrips,0,2*fgkNumberOfPSideStrips);
      fHistSSDOccupancyModule[i-500]->GetXaxis()->SetTitleColor(1);
      rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDOccupancyModule[i-500], 
					       offsRw+fSSDRawsOffset, expert, !image, !saveCorr);
      fSSDRawsOffset += 1;
    }
    
    //Occupancy per SSD ladder
    Int_t occupancyCounter = 0;
    TH1F *fHistSSDOccupancyLadder[2*(fgkSSDLADDERSLAYER5 + fgkSSDLADDERSLAYER6)];
    for(Int_t iLayer = 5; iLayer < 7; iLayer++) {
      for(Int_t iLadder = 1; iLadder < AliITSgeomTGeo::GetNLadders(iLayer) + 1; iLadder++) {
        //P-side occupancy plots
        gName  = "fHistSSD_Occupancy_Layer"; 
        gTitle = "SSD Occupancy Layer"; 
        if(iLayer == 5) {
          gName += iLayer; gName += "_Ladder"; 
          gName += 499+iLadder;
          gTitle += iLayer; gTitle += "_Ladder"; 
          gTitle += 499+iLadder;
        }
        if(iLayer == 6) {
          gName += iLayer; gName += "_Ladder"; 
          gName += 599+iLadder;
          gTitle += iLayer; gTitle += "_Ladder"; 
          gTitle += 599+iLadder;
        }
        gName += "_PSide";
        gTitle += "_PSide";
        fHistSSDOccupancyLadder[occupancyCounter] = new TH1F(gName.Data(),
                                                             Form("%s;Module number;Occupancy [%%]", gTitle.Data()),
                                                             AliITSgeomTGeo::GetNDetectors(iLayer),
                                                             0.5,AliITSgeomTGeo::GetNDetectors(iLayer)+0.5);
        fHistSSDOccupancyLadder[occupancyCounter]->GetXaxis()->SetTitleColor(1);
        rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDOccupancyLadder[occupancyCounter], 
						 offsRw+fSSDRawsOffset, expert, !image, !saveCorr);
        occupancyCounter += 1; fSSDRawsOffset += 1;
        //N-side occupancy plots
        gName = "fHistSSD_Occupancy_Layer"; 
        gTitle = "SSD Occupancy Layer"; 
        if(iLayer == 5) {
          gName += iLayer; gName += "_Ladder"; 
          gName += 499+iLadder;
          gTitle += iLayer; gTitle += "_Ladder"; 
          gTitle += 499+iLadder;
        }
        if(iLayer == 6) {
          gName += iLayer; gName += "_Ladder"; 
          gName += 599+iLadder;
          gTitle += iLayer; gTitle += "_Ladder"; 
          gTitle += 599+iLadder;
        }
        gName += "_NSide";
        gTitle += "_NSide";
        fHistSSDOccupancyLadder[occupancyCounter] = new TH1F(gName.Data(),
                                                             Form("%s;Module number;Occupancy [%%]", gTitle.Data()),
                                                             AliITSgeomTGeo::GetNDetectors(iLayer),
                                                             0.5,AliITSgeomTGeo::GetNDetectors(iLayer)+0.5);
        fHistSSDOccupancyLadder[occupancyCounter]->GetXaxis()->SetTitleColor(1);
        rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDOccupancyLadder[occupancyCounter], 
						 offsRw+fSSDRawsOffset, expert, !image, !saveCorr);
        occupancyCounter += 1; fSSDRawsOffset += 1;
      }//ladder loop
    }//layer loop

    //top level occupancy plots
    //occupancy per module - no threshold
    TH2F *fHistSSDOccupancyLayer5 = new TH2F("fHistSSDOccupancyLayer5",
					     "SSD Occupancy (Layer 5) - No threshold;N_{modules};N_{Ladders}",
					     fgkSSDMODULESPERLADDERLAYER5,
					     0,fgkSSDMODULESPERLADDERLAYER5,
					     3*fgkSSDLADDERSLAYER5,
					     500,500+fgkSSDLADDERSLAYER5);  
    fHistSSDOccupancyLayer5->GetZaxis()->SetRangeUser(0.0,100.0);
    Char_t fLabel[3];
    for(Int_t iBin = 1; iBin < fgkSSDMODULESPERLADDERLAYER5 + 1; iBin++){
      snprintf(fLabel,2,"%d",iBin);
      fHistSSDOccupancyLayer5->GetXaxis()->SetBinLabel(iBin,fLabel);
    }
    fHistSSDOccupancyLayer5->SetStats(kFALSE);
    rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDOccupancyLayer5, 
					     offsRw+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1;
    TH2F *fHistSSDOccupancyLayer6 = new TH2F("fHistSSDOccupancyLayer6",
					     "Occupancy per module (Layer 6) - No threshold;N_{modules};N_{Ladders}",
					     fgkSSDMODULESPERLADDERLAYER6,
					     0,fgkSSDMODULESPERLADDERLAYER6,
					     3*fgkSSDLADDERSLAYER6,
					     600,600+fgkSSDLADDERSLAYER6);
    fHistSSDOccupancyLayer6->GetZaxis()->SetRangeUser(0.0,100.0);
    for(Int_t iBin = 1; iBin < fgkSSDMODULESPERLADDERLAYER6 + 1; iBin++){
      snprintf(fLabel,2,"%d",iBin);
      fHistSSDOccupancyLayer6->GetXaxis()->SetBinLabel(iBin,fLabel);
    }
    fHistSSDOccupancyLayer6->SetStats(kFALSE);
    rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDOccupancyLayer6, 
					     offsRw+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1;

    //occupancy per module - threshold @ 3%
    TH2F *fHistSSDOccupancyThresholdLayer5 = new TH2F("fHistSSDOccupancyThresholdLayer5",
						      "Occupancy per module (Layer 5) - Threshold 3%;N_{modules};N_{Ladders};Entries",
						      fgkSSDMODULESPERLADDERLAYER5,
						      0,fgkSSDMODULESPERLADDERLAYER5,
						      3*fgkSSDLADDERSLAYER5,
						      500,500+fgkSSDLADDERSLAYER5);  
    fHistSSDOccupancyThresholdLayer5->GetZaxis()->SetRangeUser(3.0,10.0);
    for(Int_t iBin = 1; iBin < fgkSSDMODULESPERLADDERLAYER5 + 1; iBin++){
      snprintf(fLabel,2,"%d",iBin);
      fHistSSDOccupancyThresholdLayer5->GetXaxis()->SetBinLabel(iBin,fLabel);
    }
    fHistSSDOccupancyThresholdLayer5->SetStats(kFALSE);
    rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDOccupancyThresholdLayer5, 
					     offsRw+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1;
    TH2F *fHistSSDOccupancyThresholdLayer6 = new TH2F("fHistSSDOccupancyThresholdLayer6",
						      "Occupancy per module (Layer 6) - Threshold 3%;N_{modules};N_{Ladders}",
						      fgkSSDMODULESPERLADDERLAYER6,
						      0,fgkSSDMODULESPERLADDERLAYER6,
						      3*fgkSSDLADDERSLAYER6,
						      600,600+fgkSSDLADDERSLAYER6);
    fHistSSDOccupancyThresholdLayer6->GetZaxis()->SetRangeUser(3.0,10.0);
    for(Int_t iBin = 1; iBin < fgkSSDMODULESPERLADDERLAYER6 + 1; iBin++){
      snprintf(fLabel,2,"%d",iBin);
      fHistSSDOccupancyThresholdLayer6->GetXaxis()->SetBinLabel(iBin,fLabel);
    }
    fHistSSDOccupancyThresholdLayer6->SetStats(kFALSE);
    rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDOccupancyThresholdLayer6, 
					     offsRw+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1;

    //Average occupancy per module
    TH2F *fHistSSDAverageOccupancyLayer5 = new TH2F("fHistSSDAverageOccupancyLayer5",
						    "Average occupancy per module (Layer 5);N_{modules};N_{Ladders}",
						    fgkSSDMODULESPERLADDERLAYER5,
						    0,fgkSSDMODULESPERLADDERLAYER5,
						    3*fgkSSDLADDERSLAYER5,
						    500,500+fgkSSDLADDERSLAYER5);  
    fHistSSDAverageOccupancyLayer5->GetZaxis()->SetRangeUser(0.0,5.0);
    for(Int_t iBin = 1; iBin < fgkSSDMODULESPERLADDERLAYER5 + 1; iBin++){
      snprintf(fLabel,2,"%d",iBin);
      fHistSSDAverageOccupancyLayer5->GetXaxis()->SetBinLabel(iBin,fLabel);
    }
    fHistSSDAverageOccupancyLayer5->SetStats(kFALSE);
    rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDAverageOccupancyLayer5, 
					     offsRw+fSSDRawsOffset, !expert, image, !saveCorr);
    fSSDRawsOffset += 1;
    TH2F *fHistSSDAverageOccupancyLayer6 = new TH2F("fHistSSDAverageOccupancyLayer6",
						    "Average occupancy per module (Layer 6);N_{modules};N_{Ladders}",
						    fgkSSDMODULESPERLADDERLAYER6,
						    0,fgkSSDMODULESPERLADDERLAYER6,
						    3*fgkSSDLADDERSLAYER6,
						    600,600+fgkSSDLADDERSLAYER6);
    fHistSSDAverageOccupancyLayer6->GetZaxis()->SetRangeUser(0.0,5.0);
    for(Int_t iBin = 1; iBin < fgkSSDMODULESPERLADDERLAYER6 + 1; iBin++){
      snprintf(fLabel,2,"%d",iBin);
      fHistSSDAverageOccupancyLayer6->GetXaxis()->SetBinLabel(iBin,fLabel);
    }
    fHistSSDAverageOccupancyLayer6->SetStats(kFALSE);
    rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDAverageOccupancyLayer6, 
					     offsRw+fSSDRawsOffset, !expert, image, !saveCorr);
    fSSDRawsOffset += 1;

    //Output of the DA
    TH2F *fHistSSDPSideBadChannelMapLayer5 = new TH2F("fHistSSDPSideBadChannelMapLayer5",
						      "Layer 5;N_{module};N_{ladder}",
						      22,1,23,
						      34,500,534);
    fHistSSDPSideBadChannelMapLayer5->GetXaxis()->SetTitleColor(1);
    fHistSSDPSideBadChannelMapLayer5->SetStats(kFALSE);
    fHistSSDPSideBadChannelMapLayer5->GetYaxis()->SetTitleOffset(1.8);
    fHistSSDPSideBadChannelMapLayer5->GetXaxis()->SetNdivisions(22);
    fHistSSDPSideBadChannelMapLayer5->GetYaxis()->SetNdivisions(34);
    fHistSSDPSideBadChannelMapLayer5->GetXaxis()->SetLabelSize(0.03);
    fHistSSDPSideBadChannelMapLayer5->GetYaxis()->SetLabelSize(0.03);
    fHistSSDPSideBadChannelMapLayer5->GetZaxis()->SetTitleOffset(1.6);
    fHistSSDPSideBadChannelMapLayer5->GetZaxis()->SetTitle("Bad channels (p-side)[%]");
    rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDPSideBadChannelMapLayer5, 
					     offsRw+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1; fSSDRawsDAOffset += 1;
    
    TH2F *fHistSSDNSideBadChannelMapLayer5 = new TH2F("fHistSSDNSideBadChannelMapLayer5",
						      "Layer 5;N_{module};N_{ladder}",
						      22,1,23,
						      34,500,534);
    fHistSSDNSideBadChannelMapLayer5->GetXaxis()->SetTitleColor(1);
    fHistSSDNSideBadChannelMapLayer5->SetStats(kFALSE);
    fHistSSDNSideBadChannelMapLayer5->GetYaxis()->SetTitleOffset(1.8);
    fHistSSDNSideBadChannelMapLayer5->GetXaxis()->SetNdivisions(22);
    fHistSSDNSideBadChannelMapLayer5->GetYaxis()->SetNdivisions(34);
    fHistSSDNSideBadChannelMapLayer5->GetXaxis()->SetLabelSize(0.03);
    fHistSSDNSideBadChannelMapLayer5->GetYaxis()->SetLabelSize(0.03);
    fHistSSDNSideBadChannelMapLayer5->GetZaxis()->SetTitleOffset(1.6);
    fHistSSDNSideBadChannelMapLayer5->GetZaxis()->SetTitle("Bad channels (n-side)[%]");
    rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDNSideBadChannelMapLayer5, 
					     offsRw+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1; fSSDRawsDAOffset += 1;
    
    TH2F *fHistSSDPSideBadChannelMapLayer6 = new TH2F("fHistSSDPSideBadChannelMapLayer6",
						      "Layer 6;N_{module};N_{ladder}",
						      25,1,26,
						      38,600,638);
    fHistSSDPSideBadChannelMapLayer6->GetXaxis()->SetTitleColor(1);
    fHistSSDPSideBadChannelMapLayer6->SetStats(kFALSE);
    fHistSSDPSideBadChannelMapLayer6->GetYaxis()->SetTitleOffset(1.8);
    fHistSSDPSideBadChannelMapLayer6->GetXaxis()->SetNdivisions(25);
    fHistSSDPSideBadChannelMapLayer6->GetYaxis()->SetNdivisions(38);
    fHistSSDPSideBadChannelMapLayer6->GetXaxis()->SetLabelSize(0.03);
    fHistSSDPSideBadChannelMapLayer6->GetYaxis()->SetLabelSize(0.03);
    fHistSSDPSideBadChannelMapLayer6->GetZaxis()->SetTitleOffset(1.6);
    fHistSSDPSideBadChannelMapLayer6->GetZaxis()->SetTitle("Bad channels (p-side)[%]");
    rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDPSideBadChannelMapLayer6, 
					     offsRw+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1; fSSDRawsDAOffset += 1;
    
    TH2F *fHistSSDNSideBadChannelMapLayer6 = new TH2F("fHistSSDNSideBadChannelMapLayer6",
						      "Layer 6;N_{module};N_{ladder}",
						      25,1,26,
						      38,600,638);
    fHistSSDNSideBadChannelMapLayer6->GetXaxis()->SetTitleColor(1);
    fHistSSDNSideBadChannelMapLayer6->SetStats(kFALSE);
    fHistSSDNSideBadChannelMapLayer6->GetYaxis()->SetTitleOffset(1.8);
    fHistSSDNSideBadChannelMapLayer6->GetXaxis()->SetNdivisions(25);
    fHistSSDNSideBadChannelMapLayer6->GetYaxis()->SetNdivisions(38);
    fHistSSDNSideBadChannelMapLayer6->GetXaxis()->SetLabelSize(0.03);
    fHistSSDNSideBadChannelMapLayer6->GetYaxis()->SetLabelSize(0.03);
    fHistSSDNSideBadChannelMapLayer6->GetZaxis()->SetTitleOffset(1.6);
    fHistSSDNSideBadChannelMapLayer6->GetZaxis()->SetTitle("Bad channels (n-side)[%]");
    rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDNSideBadChannelMapLayer6, 
					     offsRw+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1; fSSDRawsDAOffset += 1;

    //Common mode values
    TH2F *fHistSSDPSideCommonModeMapLayer5 = new TH2F("fHistSSDPSideCommonModeMapLayer5",
						      "Layer 5 - CM (P-side);N_{module};N_{ladder}",
						      22,1,23,
						      34,500,534);
    fHistSSDPSideCommonModeMapLayer5->GetXaxis()->SetTitleColor(1);
    fHistSSDPSideCommonModeMapLayer5->SetStats(kFALSE);
    fHistSSDPSideCommonModeMapLayer5->GetYaxis()->SetTitleOffset(1.8);
    fHistSSDPSideCommonModeMapLayer5->GetXaxis()->SetNdivisions(22);
    fHistSSDPSideCommonModeMapLayer5->GetYaxis()->SetNdivisions(34);
    fHistSSDPSideCommonModeMapLayer5->GetXaxis()->SetLabelSize(0.03);
    fHistSSDPSideCommonModeMapLayer5->GetYaxis()->SetLabelSize(0.03);
    fHistSSDPSideCommonModeMapLayer5->GetZaxis()->SetTitleOffset(1.6);
    fHistSSDPSideCommonModeMapLayer5->GetZaxis()->SetTitle("RMS(CM) (P-side)");
    rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDPSideCommonModeMapLayer5, 
					     offsRw+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1; fSSDRawsDAOffset += 1;
    
    TH2F *fHistSSDNSideCommonModeMapLayer5 = new TH2F("fHistSSDNSideCommonModeMapLayer5",
						      "Layer 5 - CM (N-side);N_{module};N_{ladder}",
						      22,1,23,
						      34,500,534);
    fHistSSDNSideCommonModeMapLayer5->GetXaxis()->SetTitleColor(1);
    fHistSSDNSideCommonModeMapLayer5->SetStats(kFALSE);
    fHistSSDNSideCommonModeMapLayer5->GetYaxis()->SetTitleOffset(1.8);
    fHistSSDNSideCommonModeMapLayer5->GetXaxis()->SetNdivisions(22);
    fHistSSDNSideCommonModeMapLayer5->GetYaxis()->SetNdivisions(34);
    fHistSSDNSideCommonModeMapLayer5->GetXaxis()->SetLabelSize(0.03);
    fHistSSDNSideCommonModeMapLayer5->GetYaxis()->SetLabelSize(0.03);
    fHistSSDNSideCommonModeMapLayer5->GetZaxis()->SetTitleOffset(1.6);
    fHistSSDNSideCommonModeMapLayer5->GetZaxis()->SetTitle("RMS(CM) (N-side)");
    rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDNSideCommonModeMapLayer5, 
					     offsRw+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1; fSSDRawsDAOffset += 1;
    
    TH2F *fHistSSDPSideCommonModeMapLayer6 = new TH2F("fHistSSDPSideCommonModeMapLayer6",
						      "Layer 6 - CM (P-side);N_{module};N_{ladder}",
						      25,1,26,
						      38,600,638);
    fHistSSDPSideCommonModeMapLayer6->GetXaxis()->SetTitleColor(1);
    fHistSSDPSideCommonModeMapLayer6->SetStats(kFALSE);
    fHistSSDPSideCommonModeMapLayer6->GetYaxis()->SetTitleOffset(1.8);
    fHistSSDPSideCommonModeMapLayer6->GetXaxis()->SetNdivisions(25);
    fHistSSDPSideCommonModeMapLayer6->GetYaxis()->SetNdivisions(38);
    fHistSSDPSideCommonModeMapLayer6->GetXaxis()->SetLabelSize(0.03);
    fHistSSDPSideCommonModeMapLayer6->GetYaxis()->SetLabelSize(0.03);
    fHistSSDPSideCommonModeMapLayer6->GetZaxis()->SetTitleOffset(1.6);
    fHistSSDPSideCommonModeMapLayer6->GetZaxis()->SetTitle("RMS(CM) (P-side)");
    rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDPSideCommonModeMapLayer6, 
					     offsRw+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1; fSSDRawsDAOffset += 1;
    
    TH2F *fHistSSDNSideCommonModeMapLayer6 = new TH2F("fHistSSDNSideCommonModeMapLayer6",
						      "Layer 6 - CM (N-side);N_{module};N_{ladder}",
						      25,1,26,
						      38,600,638);
    fHistSSDNSideCommonModeMapLayer6->GetXaxis()->SetTitleColor(1);
    fHistSSDNSideCommonModeMapLayer6->SetStats(kFALSE);
    fHistSSDNSideCommonModeMapLayer6->GetYaxis()->SetTitleOffset(1.8);
    fHistSSDNSideCommonModeMapLayer6->GetXaxis()->SetNdivisions(25);
    fHistSSDNSideCommonModeMapLayer6->GetYaxis()->SetNdivisions(38);
    fHistSSDNSideCommonModeMapLayer6->GetXaxis()->SetLabelSize(0.03);
    fHistSSDNSideCommonModeMapLayer6->GetYaxis()->SetLabelSize(0.03);
    fHistSSDNSideCommonModeMapLayer6->GetZaxis()->SetTitleOffset(1.6);
    fHistSSDNSideCommonModeMapLayer6->GetZaxis()->SetTitle("RMS(CM) (N-side)");
    rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDNSideCommonModeMapLayer6, 
					     offsRw+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1; fSSDRawsDAOffset += 1;
  }//online flag
  
  fSSDhRawsTask = fSSDRawsOffset;
  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SSD Raws histograms booked\n",fSSDhRawsTask));
  AliDebug(AliQAv1::GetQADebugLevel(), Form("Number of histograms (SPD+SDD+SSD): %d\n",offsRw+fSSDhRawsTask));  
  AliDebug(AliQAv1::GetQADebugLevel(),Form("Number of histograms (SPD+SDD+SSD): %d\n",offsRw+fSSDRawsOffset));
  
  /*
    fSSDhTask = fSSDRawsOffset;
    AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SSD Raws histograms booked\n",fSSDhTask));
    AliDebug(AliQAv1::GetQADebugLevel(), Form("Number of histograms (SPD+SDD+SSD): %d\n",offsRw+fSSDhTask));  
    AliDebug(AliQAv1::GetQADebugLevel(),Form("Number of histograms (SPD+SDD+SSD): %d\n",offsRw+fSSDRawsOffset));
  */


  //} //event species loop

  return rv ; 


}

//____________________________________________________________________________
Int_t AliITSQASSDDataMakerRec::MakeRaws(AliRawReader* rawReader) 
{ 
  // Fill QA for RAW - SSD -
  Int_t rv = 0 ; 
  Int_t gStripNumber;
  Int_t gHistPosition;
  Int_t gLayer = 0,gLadder = 0, gModule = 0;

  Double_t gSizePerDDL[fgkNumOfDDLs] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t gSizePerLDC[fgkNumOfLDCs] = {0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t sumSSDDataSize = 0.0;
  Double_t eventSize = -1.0;

  // RS: I've commented this loop over species: makes no sense in the event loop
  //  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
  //    if (!AliQAv1::Instance()->IsEventSpecieSet(specie)) continue;
  int specie = fAliITSQADataMakerRec->GetEventSpecie(); 
  //cout << "MakeRaws: Event specie " << specie << " is set" << endl;
  //
  int offsRw = fGenRawsOffset[specie];
  //
  //AliInfo(Form("offsRw %d\n",offsRw));
  if(fkOnline) {
    //reset the signal vs strip number histograms
    for(Int_t iModule = 0; iModule < fgkSSDMODULES; iModule++)
      fHistSSDRawSignalModule[iModule]->Reset();
  }//online flag
  
  rawReader->Select("ITSSSD",-1,-1);  
  rawReader->Reset(); //rawReader->NextEvent();   
  //AliInfo(Form("offsRw %d\n",offsRw));
  fAliITSQADataMakerRec->FillRawsData(offsRw,rawReader->GetType());
  
  AliITSRawStreamSSD gSSDStream(rawReader);    
  AliRawReaderRoot *rootReader = (AliRawReaderRoot *)rawReader;
  if(rootReader) {
    const AliRawEventHeaderBase *header = rootReader->GetEventHeader();
    if(header)
      eventSize = header->GetEventSize();
  }
  while (gSSDStream.Next()) {
    if(gSSDStream.GetModuleID() < 0) continue;
    /*cout<<"DDL: "<<rawReader->GetDDLID()<<
      " - LDC: "<<rawReader->GetLDCId()<<
      " - Size: "<<rawReader->GetDataSize()<<
      " - Equipment size: "<<rawReader->GetEquipmentSize()<<endl;*/
    Int_t ddlid = rawReader->GetDDLID();
    if(ddlid<0){
      AliError("GetDDLID returns a negative value");
      continue;
    }
    gSizePerDDL[ddlid] = rawReader->GetDataSize();
    //gSizePerLDC[rawReader->GetLDCId()-8] = rawReader->GetDataSize();
    AliITSgeomTGeo::GetModuleId(gSSDStream.GetModuleID(),gLayer,gLadder,gModule);
    gHistPosition = (gLayer == 5) ? ((gLadder - 1)*fgkSSDMODULESPERLADDERLAYER5 + gModule - 1) : ((gLadder - 1)*fgkSSDMODULESPERLADDERLAYER6 + gModule + fgkSSDMODULESLAYER5 - 1);
    if(fkOnline) {
      if(gSSDStream.GetStrip() < 0) {
	//Printf("Layer: %d - Ladder: %d - Module: %d - Strip: %d - Signal: %lf",gLayer,gLadder,gModule,gSSDStream.GetStrip(),gSSDStream.GetSignal());
	if(TMath::Abs(gSSDStream.GetStrip()) < 7)
	  fHistSSDCMModule[gHistPosition]->Fill(gSSDStream.GetSignal());
	if(TMath::Abs(gSSDStream.GetStrip()) > 6)
	  fHistSSDCMModule[fgkSSDMODULES+gHistPosition]->Fill(gSSDStream.GetSignal());
      }//CM values
      else {
	gStripNumber = (gSSDStream.GetSideFlag() == 0) ? gSSDStream.GetStrip() : -gSSDStream.GetStrip() + 2*fgkNumberOfPSideStrips;
	//AliDebug(AliQAv1::GetQADebugLevel(), Form("ModulePosition: %d - Layer: %d - Ladder: %d - Module: %d\n",gHistPosition,gLayer,gLadder,gModule));
	fHistSSDRawSignalModule[gHistPosition]->Fill(gStripNumber,gSSDStream.GetSignal());
	//fAliITSQADataMakerRec->FillRawsData(offsRw+gHistPosition+fSSDRawsCommonLevelOffset,gStripNumber,gSSDStream.GetSignal());
      }//normal strip signal
    }//online flag
  }//streamer loop   
  
  //event size calculation and filling info
  for(Int_t i = 0; i < fgkNumOfDDLs; i++) {
    sumSSDDataSize += gSizePerDDL[i];
    if(gSizePerDDL[i] > 0) {
      fAliITSQADataMakerRec->FillRawsData(offsRw+3,i+512);
      fAliITSQADataMakerRec->FillRawsData(offsRw+5+i,gSizePerDDL[i]/1e+03);
      //if(i == 5)
      //cout<<gSizePerDDL[i]/1e+03<<endl;
      //cout<<"Event: "<<nSSDEventPerCycle<<" - DDL: "<<i+512<<
      //" - Data size: "<<gSizePerDDL[i]/1e+03<<endl;
    }
    //(fAliITSQADataMakerRec->FillRawsData(offsRw+4),i+512,gSizePerDDL[i]/1e+06);
  }
  for(Int_t i = 0; i < fgkNumOfLDCs; i++) {
    //LDC 170
    if(i == 0)
      gSizePerLDC[i] = gSizePerDDL[8] + gSizePerDDL[9];
    //LDC 171
    if(i == 1)
      gSizePerLDC[i] = gSizePerDDL[10] + gSizePerDDL[11];
    //LDC 172
    if(i == 2)
      gSizePerLDC[i] = gSizePerDDL[12] + gSizePerDDL[13];
    //LDC 173
    if(i == 3)
      gSizePerLDC[i] = gSizePerDDL[14] + gSizePerDDL[15];
    //LDC 174
    if(i == 4)
      gSizePerLDC[i] = gSizePerDDL[0] + gSizePerDDL[1];
    //LDC 175
    if(i == 5)
      gSizePerLDC[i] = gSizePerDDL[2] + gSizePerDDL[3];
      //LDC 176
    if(i == 6)
      gSizePerLDC[i] = gSizePerDDL[4] + gSizePerDDL[5];
    //LDC 177
    if(i == 7)
      gSizePerLDC[i] = gSizePerDDL[6] + gSizePerDDL[7];
    
    if(gSizePerLDC[i] > 0) fAliITSQADataMakerRec->FillRawsData(offsRw+21,i+170);
    fAliITSQADataMakerRec->FillRawsData(offsRw+23+i,gSizePerLDC[i]/1e+03);
    //cout<<"Event: "<<nSSDEventPerCycle<<" - LDC: "<<i+170<<
    //" - Data size: "<<gSizePerLDC[i]<<endl;
    
    //(fAliITSQADataMakerRec->FillRawsData(offsRw+22),i+6,gSizePerLDC[i]/1e+06);
  }
  
  // RS Invalid
  //    for(Int_t i = 0; i < fgkNumOfDDLs; i++) {
  //      Double_t gSizePerDDL = ((TH1F *)fAliITSQADataMakerRec->GetRawsData(offsRw+5+i))->GetMean();
  //if(i == 5)
  //cout<<"DDL: "<<i+512<<" - Size: "<<gSizePerDDL<<
  //" - Mean: "<<((TH1F *)fAliITSQADataMakerRec->GetRawsData(offsRw+5+i))->GetMean()<<endl;
  //    }
  
  if(sumSSDDataSize) fAliITSQADataMakerRec->FillRawsData(offsRw+1,sumSSDDataSize/1e+03);
  if(eventSize) fAliITSQADataMakerRec->FillRawsData(offsRw+2,100.*sumSSDDataSize/eventSize);
  
  //Occupancy calculation
  if(fkOnline) {
    for(Int_t iModule = 0; iModule < fgkSSDMODULES; iModule++) {
      GetOccupancyStrip(fHistSSDRawSignalModule[iModule],fOccupancyMatrix[iModule]);
      //if(iModule == 156) {    
      //cout<<"========================================"<<endl;
      //AliITSgeomTGeo::GetModuleId(656,gLayer,gLadder,gModule);
      /*for(Int_t iBin = 1; iBin < fHistSSDRawSignalModule[iModule]->GetXaxis()->GetNbins(); iBin++) {
	if((iBin >= 750)&&(iBin <= 780))
	cout<<"Event: "<<nSSDEventPerCycle<<" - Ladder: "<<gLadder+499<<
	" - Module: "<<gModule<<" - Strip: "<<iBin<<
	" - Signal: "<<fHistSSDRawSignalModule[iModule]->GetBinContent(iBin)<<endl;
	}*///strip loop --> to be removed
      //}//module cut --> to be removed
    }//module loop
  }//online flag for SSD
  
  // RS Invalid
  //cout<<"Event type entries: "<<((TH1*)fAliITSQADataMakerRec->GetRawsData(offsRw))->GetEntries()<<endl;
  //cout<<"DDL id entries at MR: "<<((TH1F *)fAliITSQADataMakerRec->GetRawsData(offsRw+3))->GetEntries()<< endl;    
  //cout<<"LDC id entries at MR: "<<((TH1F *)fAliITSQADataMakerRec->GetRawsData(offsRw+21))->GetEntries()<< endl; 
  
  //  } //event species loop // commented by RS
  //
  return rv ; 
}

//____________________________________________________________________________ //
void AliITSQASSDDataMakerRec::GetOccupancyStrip(TH1 *lHisto, Int_t *occupancyMatrix) { 
  //Increments the entries in the occupancy matrix based 
  //on whether the signal for each strip is larger than the cutoff
  //Currently the cutoff is at 0 which means that if ZS
  //works, whatever comes from the FEROM is considered as "signal"
  //Double_t cutoff = 0.0;
  TString histname = lHisto->GetName();
  //cout<<histname.Data()<<endl;
  //if(histname.Contains("Layer5_Ladder8_Module3")) {
  for(Int_t iBin = 1; iBin < lHisto->GetXaxis()->GetNbins(); iBin++) {
    Double_t y = lHisto->GetBinContent(iBin);
    if(y) {
      occupancyMatrix[iBin-1] += 1;
    }
    //if((iBin >= 750)&&(iBin <= 780))
    //cout<<"Strip: "<<iBin<<
    //" - Signal: "<<y<<
    //" - Occupancy: "<<occupancyMatrix[iBin-1]<<endl;
    
  }
  //}
}

//____________________________________________________________________________ 
Double_t AliITSQASSDDataMakerRec::GetOccupancyModule(TH1 *lHisto, 
						     Int_t stripside,
						     Int_t mode,
						     Double_t threshold) {
  //Mode 0: calculates the occupancy of a module 
  //        without a threshold in the strip occupancy
  //Mode 1: calculates the occupancy of a module 
  //        with the set threshold in the strip occupancy
  //Mode 2: calculates the average occupancy of a module 

  // bo: TDC >0 or # of sigmas wrt noise ?
  //stripside == 0 --> P-side
  //stripside == 1 --> N-side
  Int_t lNumFiredBins = 0;
  Double_t sumOccupancy = 0.0;
  TString histname = lHisto->GetName();
  for(Int_t iBin = 1 + stripside*fgkNumberOfPSideStrips; iBin < fgkNumberOfPSideStrips*(1 + stripside); iBin++){
    /*if(histname.Contains("Layer5_Ladder507_Module3")) {
      cout<<lHisto->GetName()<<
      " - Strip: "<<iBin<<
      " - Bin content: "<<lHisto->GetBinContent(iBin)<<endl;
      }*/
    if (lHisto->GetBinContent(iBin) > threshold) {
      lNumFiredBins++; 
      sumOccupancy = lHisto->GetBinContent(iBin);
    }
  }
  
  Double_t lOccupancy = 0.0;
  if((mode == 0)||(mode == 1))
    lOccupancy = (100.*lNumFiredBins)/fgkNumberOfPSideStrips; // percentage  
  else if(mode == 2)
    lOccupancy = 100.*sumOccupancy/fgkNumberOfPSideStrips;

  /* if(histname.Contains("Layer5_Ladder507_Module3"))
     cout<<"Fired strips: "<<lNumFiredBins<<
     " - Occupancy: "<<lOccupancy<<endl;*/
  ;
  //AliDebug(AliQAv1::GetQADebugLevel(), Form("Fired strips: %d - Total strips: %d - Occupancy :%lf\n",lNumFiredBins,lHisto->GetNbinsX(),lOccupancy));
  
  return lOccupancy;
}

//____________________________________________________________________________
void AliITSQASSDDataMakerRec::MonitorCMValues(Int_t trCl) 
{
  //Monitor in AMORE the CM values
  //for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
  //if (!AliQAv1::Instance()->IsEventSpecieSet(specie)) continue;
  //cout << "MonitorCMValues: Event specie " << specie << " is set" << endl;
  Int_t specie = fAliITSQADataMakerRec->GetEventSpecie();
  int offsRw = fGenRawsOffset[specie];
   
  //compute the rms of the CM values
  Int_t gLayer = 0, gLadder = 0, gModule = 0;
  Double_t rmsPsideCM = 0.0, rmsNsideCM = 0.0;
  for(Int_t i = 0; i < fgkSSDMODULES; i++) {
    rmsPsideCM = 0.0; rmsNsideCM = 0.0;
    AliITSgeomTGeo::GetModuleId(i+500,gLayer,gLadder,gModule);
    //Printf("%s - %s",((TH1*)list->At(i))->GetName(),
    //((TH1*)list->At(1698+i))->GetName());
    rmsPsideCM = fHistSSDCMModule[i]->GetRMS();
    rmsNsideCM = fHistSSDCMModule[fgkSSDMODULES+i]->GetRMS();
    //Printf("rmsPside: %lf - rmsNside: %lf",rmsPsideCM,rmsNsideCM);
    if(gLayer == 5) {
      TH2* h4 = (TH2*)fAliITSQADataMakerRec->GetRawsData(offsRw+fSSDRawsOffset-fSSDRawsDAOffset+4,trCl);
      if (h4) h4->SetBinContent(gModule,gLadder,rmsPsideCM);
      TH2* h5 = (TH2*)fAliITSQADataMakerRec->GetRawsData(offsRw+fSSDRawsOffset-fSSDRawsDAOffset+5,trCl);
      if (h5) h5->SetBinContent(gModule,gLadder,rmsNsideCM);
    }
    if(gLayer == 6) {
      TH2* h6 = (TH2*)fAliITSQADataMakerRec->GetRawsData(offsRw+fSSDRawsOffset-fSSDRawsDAOffset+6,trCl);
      if (h6) h6->SetBinContent(gModule,gLadder,rmsPsideCM);
      TH2* h7 = (TH2*)fAliITSQADataMakerRec->GetRawsData(offsRw+fSSDRawsOffset-fSSDRawsDAOffset+7,trCl);
      if (h7) h7->SetBinContent(gModule,gLadder,rmsNsideCM);
    }
  }//module loopcdcd
   //}//event species loop
}

//____________________________________________________________________________
void AliITSQASSDDataMakerRec::MonitorOCDBObjects(Int_t trCl) 
{ 
  //Monitor in AMORE the output of the DA
  //Currently only the bad channel list is monitored
  //Todo: Noise - Pedestal


  //for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
  //if (!AliQAv1::Instance()->IsEventSpecieSet(specie)) continue;
  //cout << "MonitorOCDBObjects: Event specie " << specie << " is set" << endl;
  Int_t specie = fAliITSQADataMakerRec->GetEventSpecie();
  int offsRw = fGenRawsOffset[specie];
  fAliITSQADataMakerRec->SetEventSpecie(AliRecoParam::ConvertIndex(specie));
  
  //((TH2F *)fAliITSQADataMakerRec->ResetRawsData(offsRw+fSSDRawsOffset-fSSDRawsDAOffset));
  //((TH2F *)fAliITSQADataMakerRec->ResetRawsData(offsRw+fSSDRawsOffset-fSSDRawsDAOffset+1));
  //((TH2F *)fAliITSQADataMakerRec->ResetRawsData(offsRw+fSSDRawsOffset-fSSDRawsDAOffset+2));
  //((TH2F *)fAliITSQADataMakerRec->ResetRawsData(offsRw+fSSDRawsOffset-fSSDRawsDAOffset+3));
  
  AliCDBEntry *entryBadChannelsSSD = fCDBManager->Get("ITS/Calib/BadChannelsSSD");
  if(!entryBadChannelsSSD)AliError("OCDB entry for the bad channel list is not valid!"); 
  AliITSBadChannelsSSDv2 *badchannelsSSD = NULL;
  if(entryBadChannelsSSD)badchannelsSSD = (AliITSBadChannelsSSDv2 *)entryBadChannelsSSD->GetObject();
  if(!badchannelsSSD)
    AliError("Bad channel list object is not a valid AliITSBadChannelsSSD object!");

  //_____________________________________________________________________________//                       
  Int_t nBadPSideChannels = 0, nBadNSideChannels = 0;
  Int_t layer = 0, ladder = 0, module = 0;
  Int_t nPSideChannelsLayer5 = 0, nNSideChannelsLayer5 = 0;
  Int_t nPSideChannelsLayer6 = 0, nNSideChannelsLayer6 = 0;
  //_____________________________________________________________________________//                      

  for(Int_t i = 0; i < fgkSSDMODULES; i++) {
    AliITSgeomTGeo::GetModuleId(i+500,layer,ladder,module);
    nBadPSideChannels = 0, nBadNSideChannels = 0;
    nPSideChannelsLayer5 = 0, nNSideChannelsLayer5 = 0;
    nPSideChannelsLayer6 = 0, nNSideChannelsLayer6 = 0;

    Int_t badChannel = 0;
    if(badchannelsSSD){
      for(Int_t j = 0; j < fgkNumberOfPSideStrips; j++) {
	badChannel = (Int_t)(badchannelsSSD->GetBadChannelP(i,j));
	if(badChannel != 0) {
	  if(layer == 5)
	    nPSideChannelsLayer5 += 1;
	  if(layer == 6)
	    nPSideChannelsLayer6 += 1;
	  nBadPSideChannels += 1;
	}//badchannel flag != 0
	badChannel = (Int_t)(badchannelsSSD->GetBadChannelN(i,j));
	if(badChannel != 0) {
	  if(layer == 5)
	    nNSideChannelsLayer5 += 1;
	  if(layer == 6)
	    nNSideChannelsLayer6 += 1;
	  nBadNSideChannels += 1;
	}//badchannel flag != 0
      }//loop over strips
    }

    //cout << "Bad channels P side module " << module << ": " << nBadPSideChannels << endl;
    //cout << "Bad channels N side module " << module << ": " << nBadNSideChannels << endl;
    
    if(layer == 5) {
      /*if((module == 10)&&(ladder == 10)) {
	cout<<"Npside bad: "<<nPSideChannelsLayer5<<" - Total: "<<fgkNumberOfPSideStrips<<" - Percentage: "<<(100.*nPSideChannelsLayer5/fgkNumberOfPSideStrips)<<endl;
	}*/
      TH2* h0 = (TH2*)fAliITSQADataMakerRec->GetRawsData(offsRw+fSSDRawsOffset-fSSDRawsDAOffset,trCl);
      if (h0) {
	if(nPSideChannelsLayer5 > 0) h0->SetBinContent(module,ladder,100.*nPSideChannelsLayer5/fgkNumberOfPSideStrips);
	else                         h0->SetBinContent(module,ladder,0.0001);
      }
      TH2* h1 = (TH2*)fAliITSQADataMakerRec->GetRawsData(offsRw+fSSDRawsOffset-fSSDRawsDAOffset+1,trCl);
      if (h1) {
	if(nNSideChannelsLayer5 > 0) h1->SetBinContent(module,ladder,100.*nNSideChannelsLayer5/fgkNumberOfPSideStrips);
	else                         h1->SetBinContent(module,ladder,0.0001);
      }
    }//layer 5                                                                    
    if(layer == 6) {
      TH2* h2 = (TH2*)fAliITSQADataMakerRec->GetRawsData(offsRw+fSSDRawsOffset-fSSDRawsDAOffset+2,trCl);
      if (h2) {
	if(nPSideChannelsLayer6 > 0) h2->SetBinContent(module,ladder,100.*nPSideChannelsLayer6/fgkNumberOfPSideStrips);
	else                         h2->SetBinContent(module,ladder,0.0001);
      }
      TH2* h3 = (TH2*)fAliITSQADataMakerRec->GetRawsData(offsRw+fSSDRawsOffset-fSSDRawsDAOffset+3,trCl);      
      if (h3) {
	if(nNSideChannelsLayer6 > 0) h3->SetBinContent(module,ladder,100.*nNSideChannelsLayer6/fgkNumberOfPSideStrips);
	else                         h3->SetBinContent(module,ladder,0.0001);
      }
    }//layer 6                                                              
  }//module loop
  
  //cout << "entries bad channel layer 5 n side " << ((TH2F *)fAliITSQADataMakerRec->GetRawsData(offsRw+fSSDRawsOffset-fSSDRawsDAOffset+1))->GetEntries() << " - Bad channels P side layer 5 module 10 ladder 10: " << ((TH2F *)fAliITSQADataMakerRec->GetRawsData(offsRw+fSSDRawsOffset-fSSDRawsDAOffset))->GetBinContent(10,10)<<endl;


  //} //event species loop


}

//____________________________________________________________________________ 
Int_t AliITSQASSDDataMakerRec::InitDigits() { 
  // Initialization for DIGIT data - SSD -
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  Int_t rv = 0 ; 
  //  fGenDigitsOffset = (fAliITSQADataMakerRec->fDigitsQAList[AliRecoParam::kDefault])->GetEntries();
  
  // custom code here
  TH1F *fHistSSDModule = new TH1F("fHistSSDDigitsModule",
                                  ";SSD Module Number;N_{DIGITS}",
                                  1698,499.5,2197.5);  
  rv =  fAliITSQADataMakerRec->Add2DigitsList(fHistSSDModule,
					      fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 0, !expert, image);
  fSSDhDigitsTask += 1;
  TH2F *fHistSSDModuleStrip = new TH2F("fHistSSDDigitsModuleStrip",
                                       ";N_{Strip};N_{Module}",
                                       1540,0,1540,1698,499.5,2197.5);  
  rv = fAliITSQADataMakerRec->Add2DigitsList(fHistSSDModuleStrip,
					     fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 1, !expert, image);
  fSSDhDigitsTask += 1;
  
  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SSD Digits histograms booked\n",fSSDhDigitsTask));
  return rv ; 
}

//____________________________________________________________________________
Int_t AliITSQASSDDataMakerRec::MakeDigits(TTree *digits) 
{ 
  // Fill QA for DIGIT - SSD -
  //  AliITS *fITS  = (AliITS*)gAlice->GetModule("ITS");
  //  fITS->SetTreeAddress();
  //  TClonesArray *iSSDdigits  = fITS->DigitsAddress(2);
  //  
  Int_t rv = 0 ; 
  TBranch *branchD = digits->GetBranch("ITSDigitsSSD");
  if (!branchD) { 
    AliError("can't get the branch with the SSD ITS digits !");
    return rv;
  }
  
  static TClonesArray statDigits("AliITSdigitSSD");
  TClonesArray *iSSDdigits = &statDigits;
  branchD->SetAddress(&iSSDdigits);  
  for(Int_t iModule = 500; iModule < 2198; iModule++) {
    iSSDdigits->Clear();
    digits->GetEvent(iModule);    
    Int_t ndigits = iSSDdigits->GetEntries();
    //printf("Module = %i \t ndigits = %i\t offset = %i \n",iModule,ndigits,fAliITSQADataMakerRec->GetEventSpecie() );
    fAliITSQADataMakerRec->FillDigitsData(fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 0,iModule,ndigits);
    if(ndigits != 0)
      AliDebug(AliQAv1::GetQADebugLevel(),Form("Module: %d - Digits: %d",iModule,ndigits));
    
    for (Int_t iDigit = 0; iDigit < ndigits; iDigit++) {
      AliITSdigit *dig = (AliITSdigit*)iSSDdigits->UncheckedAt(iDigit);
      Int_t fStripNumber = (dig->GetCoord1() == 0) ? dig->GetCoord2() : dig->GetCoord2() + fgkNumberOfPSideStrips;
      fAliITSQADataMakerRec->FillDigitsData(fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 1,fStripNumber,iModule,dig->GetSignal());
    }//digit loop
  }//module loop
  //
  return rv ; 
}

//____________________________________________________________________________ 
Int_t AliITSQASSDDataMakerRec::InitRecPoints()
{
  // Initialization for RECPOINTS - SSD -
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  Int_t rv = 0 ; 
  //AliInfo(Form("fAliITSQADataMakerRec->GetEventSpecie() %d\n",fAliITSQADataMakerRec->GetEventSpecie()));
  //AliInfo(Form("fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] %d\n",fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]));
  //  fGenRecPointsOffset = (fAliITSQADataMakerRec->fRecPointsQAList[AliRecoParam::kDefault])->GetEntries();
  //AliDebug(AliQAv1::GetQADebugLevel(), Form("**-------*-*-*-*-*-*-***************AliITSQASSDataMakerRec::MakeRecpoints offset %d \t %d \n",fGenOffset,fGenRecPointsOffset));
  Int_t nModuleOffset = 500;
  Int_t nITSTotalModules = AliITSgeomTGeo::GetNModules();

  TH1F *fHistSSDModuleIdLayer5 = new TH1F("fHistSSDModuleIdLayer5",
					  "Module Id - Layer 5;Module Id;Entries",
					  fgkSSDMODULESLAYER5,
					  nModuleOffset - 0.5,
					  nITSTotalModules-fgkSSDMODULESLAYER6+0.5);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDModuleIdLayer5, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 0, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDModuleIdLayer6 = new TH1F("fHistSSDModuleIdLayer6",
					  "Module Id - Layer 6;Module Id;Entries",
					  fgkSSDMODULESLAYER6,
					  nModuleOffset+fgkSSDMODULESLAYER5 - 0.5,
					  nITSTotalModules + 0.5);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDModuleIdLayer6, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 1, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDClusterPerEventLayer5 = new TH1F("fHistSSDClusterPerEventLayer5",
						 "N_{clusters} - Layer 5;N_{clusters};Entries;",
						 100,0.1,5000);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDClusterPerEventLayer5,
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 2, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDClusterPerEventLayer6 = new TH1F("fHistSSDClusterPerEventLayer6",
						 "N_{clusters} - Layer 6;N_{clusters};Entries;",
						 100,0.1,5000);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDClusterPerEventLayer6,
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 3, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDLocalXLayer5 = new TH1F("fHistSSDLocalXLayer5",
					"Local x coord.- Layer 5;x [cm];Entries;",
					100,-4.,4.);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDLocalXLayer5,
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 4, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDLocalXLayer6 = new TH1F("fHistSSDLocalXLayer6",
					"Local x coord.- Layer 6;x [cm];Entries;",
					100,-4.,4.);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDLocalXLayer6, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 5, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDLocalZLayer5 = new TH1F("fHistSSDLocalZLayer5",
					"Local z coord.- Layer 5;z [cm];Entries;",
					100,-4.,4.);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDLocalZLayer5, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 6, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDLocalZLayer6 = new TH1F("fHistSSDLocalZLayer6",
					"Local z coord.- Layer 6;z [cm];Entries;",
					100,-4.,4.);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDLocalZLayer6, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 7, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDGlobalXLayer5 = new TH1F("fHistSSDGlobalXLayer5",
					 "Global x - Layer 5;x [cm];Entries;",
					 100,-40.,40.);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDGlobalXLayer5, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 8, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDGlobalXLayer6 = new TH1F("fHistSSDGlobalXLayer6",
					 "Global x - Layer 6;x [cm];Entries;",
					 100,-45.,45.);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDGlobalXLayer6, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 9, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDGlobalYLayer5 = new TH1F("fHistSSDGlobalYLayer5",
					 "Global y - Layer 5;y [cm];Entries;",
					 100,-40.,40);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDGlobalYLayer5, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 10, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDGlobalYLayer6 = new TH1F("fHistSSDGlobalYLayer6",
					 "Global y - Layer 6;y [cm];Entries;",
					 100,-45.,45.);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDGlobalYLayer6, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 11, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDGlobalZLayer5 = new TH1F("fHistSSDGlobalZLayer5",
					 "Global z - Layer 5;z [cm];Entries;",
					 100,-45.,45);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDGlobalZLayer5, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 12, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDGlobalZLayer6 = new TH1F("fHistSSDGlobalZLayer6",
					 "Global z - Layer 6;z [cm];Entries;",
					 100,-55.,55.);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDGlobalZLayer6, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 13, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDPhiLayer5 = new TH1F("fHistSSDPhiLayer5",
				     "#phi - Layer 5;#phi [rad];Entries;",
				     100,-TMath::Pi(),TMath::Pi());
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDPhiLayer5, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 14, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDPhiLayer6 = new TH1F("fHistSSDPhiLayer6",
				     "#phi - Layer 6;#phi [rad];Entries;",
				     100,-TMath::Pi(),TMath::Pi());
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDPhiLayer6, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 15, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDThetaLayer5 = new TH1F("fHistSSDThetaLayer5",
				       "#theta - Layer 5;#theta [rad];Entries;",
				       100,-TMath::Pi(),TMath::Pi());
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDThetaLayer5, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 16, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDThetaLayer6 = new TH1F("fHistSSDThetaLayer6",
				       "#theta - Layer 6;#theta [rad];Entries;",
				       100,-TMath::Pi(),TMath::Pi());
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDThetaLayer6, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 17, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDRadiusLayer5 = new TH1F("fHistSSDRadiusLayer5",
					"r - Layer 5;r [cm];Entries;",
					100,35.,50.);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDRadiusLayer5, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 18, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDRadiusLayer6 = new TH1F("fHistSSDRadiusLayer6",
					"r - Layer 6;r [cm];Entries;",
					100,35.,50.);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDRadiusLayer6, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 19, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDClusterTypeLayer5 = new TH1F("fHistSSDClusterTypeLayer5",
					     "CL type - Layer 5;Cluster type;Entries;",
					     150,0,150);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDClusterTypeLayer5, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 20, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDClusterTypeLayer6 = new TH1F("fHistSSDClusterTypeLayer6",
					     "CL type - Layer 6;Cluster type;Entries;",
					     150,0,150);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDClusterTypeLayer6, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 21, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDChargeRatioLayer5 = new TH1F("fHistSSDChargeRatioLayer5",
					     "Charge ratio - Layer 5;q_{ratio};Entries;",
					     100,-2.0,2.0);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDChargeRatioLayer5, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 22, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDChargeRatioLayer6 = new TH1F("fHistSSDChargeRatioLayer6",
					     "Charge ratio - Layer 6;q_{ratio};Entries;",
					     100,-2.0,2.0);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDChargeRatioLayer6, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 23, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDChargekeVLayer5 = new TH1F("fHistSSDChargekeVLayer5",
					   "Charge - Layer 5;q [keV];Entries;",
					   100,0.,300.);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDChargekeVLayer5, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 24, !expert, image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDChargekeVLayer6 = new TH1F("fHistSSDChargekeVLayer6",
					   "Charge - Layer 6;q [keV];Entries;",
					   100,0.,300.);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDChargekeVLayer6, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 25, !expert, image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDChargePSideLayer5 = new TH1F("fHistSSDChargePSideLayer5",
					     "Charge P- Layer 5;q_{P} [keV];Entries;",
					     100,0.,300.);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDChargePSideLayer5,
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 26, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDChargePSideLayer6 = new TH1F("fHistSSDChargePSideLayer6",
					     "Charge P- Layer 6;q_{P} [keV];Entries;",
					     100,0.,300.);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDChargePSideLayer6,
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 27, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDChargeNSideLayer5 = new TH1F("fHistSSDChargeNSideLayer5",
					     "Charge N- Layer 5;q_{N} [keV];Entries;",
					     100,0.,300.);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDChargeNSideLayer5,
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 28, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDChargeNSideLayer6 = new TH1F("fHistSSDChargeNSideLayer6",
					     "Charge N- Layer 6;q_{N} [keV];Entries;",
					     100,0.,300.);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDChargeNSideLayer6,
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 29, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDChargeRatio2Layer5 = new TH1F("fHistSSDChargeRatio2Layer5",
					      "Charge Ratio qN/qP - Layer 5;q_{N}/q_{P};Entries;",
					      100,0,2);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDChargeRatio2Layer5,
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 30, expert, !image);
  fSSDhRecPointsTask += 1;
  TH1F *fHistSSDChargeRatio2Layer6 = new TH1F("fHistSSDChargeRatio2Layer6",
					      "Charge Ratio qN/qP - Layer 6;q_{N}/q_{P};Entries;",
					      100,0,2);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDChargeRatio2Layer6,
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 31, expert, !image);
  fSSDhRecPointsTask += 1;
  TH2F *fHistSSDChargePNSideLayer5 = new TH2F("fHistSSDChargePNSideLayer5",
					      "Charge correlation - Layer 5;q_{P} [keV];q_{N} [keV]",
					      100,0.,300.,
					      100,0.,300.);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDChargePNSideLayer5,
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 32, expert, !image);
  fSSDhRecPointsTask += 1;
  TH2F *fHistSSDChargePNSideLayer6 = new TH2F("fHistSSDChargePNSideLayer6",
					      "Charge correlation - Layer 6;q_{P} [keV];q_{N} [keV]",
					      100,0.,300.,
					      100,0.,300.);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDChargePNSideLayer6,
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 33, expert, !image);
  fSSDhRecPointsTask += 1;
  TH2F *fHistSSDChargeMapLayer5 = new TH2F("fHistSSDChargeMapLayer5",
					   "Charge map;N_{modules};N_{Ladders}",
					   fgkSSDMODULESPERLADDERLAYER5,
					   -0.5,fgkSSDMODULESPERLADDERLAYER5+0.5,
					   3*fgkSSDLADDERSLAYER5,
					   -0.5,fgkSSDLADDERSLAYER5+0.5);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDChargeMapLayer5, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 34, expert, !image);
  fSSDhRecPointsTask += 1;
  TH2F *fHistSSDChargeMapLayer6 = new TH2F("fHistSSDChargeMapLayer6",
					   "Charge map;N_{modules};N_{Ladders}",
					   fgkSSDMODULESPERLADDERLAYER6,
					   -0.5,fgkSSDMODULESPERLADDERLAYER6+0.5,
					   3*fgkSSDLADDERSLAYER6,
					   -0.5,fgkSSDLADDERSLAYER6+0.5);
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDChargeMapLayer6, 
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 35, expert, !image);
  fSSDhRecPointsTask += 1;
  TH2F *fHistSSDClusterMapLayer5 = new TH2F("fHistSSDClusterMapLayer5",
					    "Layer 5;N_{module};N_{ladder}",
					    22,1,23,
					    34,500,534);
  fHistSSDClusterMapLayer5->GetXaxis()->SetTitleColor(1);
  fHistSSDClusterMapLayer5->SetStats(kFALSE);
  fHistSSDClusterMapLayer5->GetYaxis()->SetTitleOffset(1.8);
  fHistSSDClusterMapLayer5->GetXaxis()->SetNdivisions(22);
  fHistSSDClusterMapLayer5->GetYaxis()->SetNdivisions(34);
  fHistSSDClusterMapLayer5->GetXaxis()->SetLabelSize(0.03);
  fHistSSDClusterMapLayer5->GetYaxis()->SetLabelSize(0.03);
  fHistSSDClusterMapLayer5->GetZaxis()->SetTitleOffset(1.4);
  fHistSSDClusterMapLayer5->GetZaxis()->SetTitle("N_{clusters}");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDClusterMapLayer5,
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 36, !expert, image);
  fSSDhRecPointsTask += 1;
  TH2F *fHistSSDClusterMapLayer6 = new TH2F("fHistSSDClusterMapLayer6",
					    "Layer 6;N_{module};N_{ladder}",
					    25,1,26,
					    38,600,638);
  fHistSSDClusterMapLayer6->GetXaxis()->SetTitleColor(1);
  fHistSSDClusterMapLayer6->SetStats(kFALSE);
  fHistSSDClusterMapLayer6->GetYaxis()->SetTitleOffset(1.8);
  fHistSSDClusterMapLayer6->GetXaxis()->SetNdivisions(25);
  fHistSSDClusterMapLayer6->GetYaxis()->SetNdivisions(38);
  fHistSSDClusterMapLayer6->GetXaxis()->SetLabelSize(0.03);
  fHistSSDClusterMapLayer6->GetYaxis()->SetLabelSize(0.03);
  fHistSSDClusterMapLayer6->GetZaxis()->SetTitleOffset(1.4);
  fHistSSDClusterMapLayer6->GetZaxis()->SetTitle("N_{clusters}");
  rv = fAliITSQADataMakerRec->Add2RecPointsList(fHistSSDClusterMapLayer6,
						fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 37, !expert, image);
  fSSDhRecPointsTask += 1;
  //printf ("%d SSD Recs histograms booked\n",fSSDhRecPointsTask);
  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SSD Recs histograms booked\n",fSSDhRecPointsTask));
  return rv ; 
}

//____________________________________________________________________________ 
Int_t AliITSQASSDDataMakerRec::MakeRecPoints(TTree *clustersTree)
{
  // Fill QA for recpoints - SSD -
  //printf("*-*-*-*-*-*-*---*-*-*-------*-*-*-*-*-*-***************AliITSQASSDataMakerRec::MakeRecpoints called \n");
  //
  Int_t rv = 0 ; 
  Int_t gLayer = 0, gLadder = 0, gModule = 0;
  Int_t lLadderLocationY = 0;
  TClonesArray *recpoints = NULL;
  AliITSRecPointContainer* rpcont=AliITSRecPointContainer::Instance();
  rpcont->FetchClusters(0,clustersTree); 
  if(!rpcont->GetStatusOK() || !rpcont->IsSSDActive()){
    AliError("can't get SSD clusters !");
    return rv;
  }
 
  //AliInfo(Form("fAliITSQADataMakerRec->GetEventSpecie() %d\n",fAliITSQADataMakerRec->GetEventSpecie()));
  //AliInfo(Form("fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] %d\n",fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]));
  Int_t nClustersLayer5 = 0, nClustersLayer6 = 0;
  Int_t npoints = 0;      
  Float_t cluglo[3]={0.,0.,0.}; 
  //printf("*-*-*-*-*-*-*---*-*-*-------*-*-*-*-*-*-***************AliITSQASSDataMakerRec::MakeRecpoints STEP1 \n");
  // AliITSgeomTGeo::GetModuleIndex() issues an error in case the arguments
  // are illegal and returns -1
  Int_t firMod = TMath::Max(0,AliITSgeomTGeo::GetModuleIndex(5,1,1));
  Int_t lasMod =  AliITSgeomTGeo::GetNModules();
  for(Int_t module = firMod; module < lasMod; module++){
    recpoints = rpcont->UncheckedGetClusters(module);
    npoints += recpoints->GetEntries();
    AliITSgeomTGeo::GetModuleId(module,gLayer,gLadder,gModule);
    //printf("SSDDataMAkerRec:::::::::::::::::::::::gLayer ========== %d \n\n",gLayer);
    lLadderLocationY = 3*gLadder;
    ////printf("*-*-*-*-*-*-*---*-*-*-------*-*-*-*-*-*-***************AliITSQASSDataMakerRec::MakeRecpoints inside loop \n");
    for(Int_t j = 0;j < recpoints->GetEntries(); j++){
      ////printf("*-*-*-*-*-*-*---*-*-*-------*-*-*-*-*-*-***************AliITSQASSDataMakerRec::MakeRecpoints inside loop 2\n");
      AliITSRecPoint *recp = (AliITSRecPoint*)recpoints->At(j);    
      Int_t layer = recp->GetLayer();
      //printf("SSDDataMAkerRec:::::::::::::::::::::::layer ========== %d \n\n",layer);
      recp->GetGlobalXYZ(cluglo);
      Float_t radius = TMath::Sqrt(cluglo[0]*cluglo[0]+cluglo[1]*cluglo[1]); 
      Float_t phi = TMath::ATan2(cluglo[1],cluglo[0]);
      Float_t theta = TMath::ATan2(radius,cluglo[2]);
      Double_t chargeRatio = recp->GetChargeRatio();
      Double_t clusterCharge = recp->GetQ();
      Double_t chargePSide = clusterCharge*(1. + chargeRatio);
      Double_t chargeNSide = clusterCharge*(1. - chargeRatio);
      if(layer == 4) {
	//printf("-*-*-*-*-*-***************AliITSQASSDataMakerRec::MakeRecpoints Filling 4 called \n");
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 20,recp->GetType());

	if(recp->GetType() != 1) continue;
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 0,module);
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 4,recp->GetDetLocalX());
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 6,recp->GetDetLocalZ());
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 8,cluglo[0]);
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 10,cluglo[1]);
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 12,cluglo[2]);
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 14,phi);
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 16,theta);
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 18,radius);
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 22,recp->GetChargeRatio());
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 24,recp->GetQ());
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 26,chargePSide);
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 28,chargeNSide);
	if(chargePSide != 0.) fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 30,chargeNSide/chargePSide);
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 32,chargePSide,chargeNSide);
	fAliITSQADataMakerRec->SetRecPointsDataBinContent(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 34,gModule,lLadderLocationY,recp->GetQ());
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 36,gModule,499+gLadder,1);
	nClustersLayer5 += 1;
      }//layer 5 histograms
      if(layer == 5) {
	//printf("-*-*-*-*-*-***************AliITSQASSDataMakerRec::MakeRecpoints Filling 5 called \n");
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 21,recp->GetType());

	if(recp->GetType() != 1) continue;
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 1,module);
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 5,recp->GetDetLocalX());
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 7,recp->GetDetLocalZ());
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 9,cluglo[0]);
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 11,cluglo[1]);
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 13,cluglo[2]);
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 15,phi);
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 17,theta);
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 19,radius);
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 23,recp->GetChargeRatio());
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 25,recp->GetQ());
        fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 27,chargePSide);
        fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 29,chargeNSide);
        if(chargePSide != 0.) fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 31,chargeNSide/chargePSide);
        fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 33,chargePSide,chargeNSide);
	fAliITSQADataMakerRec->SetRecPointsDataBinContent(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 35,gModule,lLadderLocationY,recp->GetQ());
	fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 37,gModule,599+gLadder,1);
	nClustersLayer6 += 1;
      }//layer 6 histograms
    }//rec. points loop
  }//module loop
  
  //printf("-*-*-*-*-*-***************AliITSQASSDataMakerRec::MakeRecpoints Filling called \n");
  fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 2,nClustersLayer5);
  fAliITSQADataMakerRec->FillRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 3,nClustersLayer6);
  //
  return rv ; 
}

//____________________________________________________________________________ 
Int_t AliITSQASSDDataMakerRec::GetOffset(AliQAv1::TASKINDEX_t task,Int_t specie) {
  // Returns offset number according to the specified task 
  Int_t offset=0;
  if( task == AliQAv1::kRAWS ) {
    offset=fGenRawsOffset[specie];  
  }
  else if( task == AliQAv1::kDIGITSR ) {
    offset=fGenDigitsOffset[specie];   
  }
  else if( task == AliQAv1::kRECPOINTS ) {
    offset=fGenRecPointsOffset[specie];   
  }

  return offset;
}

//_______________________________________________________________

void AliITSQASSDDataMakerRec::SetOffset(AliQAv1::TASKINDEX_t task, Int_t offset, Int_t specie) {
  // Returns offset number according to the specified task
  if( task == AliQAv1::kRAWS ) {
    fGenRawsOffset[specie]=offset;
  }
  else if( task == AliQAv1::kDIGITSR ) {
    fGenDigitsOffset[specie]=offset;
  }
  else if( task == AliQAv1::kRECPOINTS ) {
    fGenRecPointsOffset[specie]=offset;
  }
}

//____________________________________________________________________________ 
Int_t AliITSQASSDDataMakerRec::GetTaskHisto(AliQAv1::TASKINDEX_t task) {
  // Returns the number of histograms associated to the specified task
  Int_t histotot=0;

  if( task == AliQAv1::kRAWS ) {
    histotot=fSSDhRawsTask;  
  }
  else if( task == AliQAv1::kDIGITSR ) {
    histotot=fSSDhDigitsTask;
  }
  else if( task == AliQAv1::kRECPOINTS ){
    histotot=fSSDhRecPointsTask;   
  }
  else { 
    AliWarning("No task has been selected. TaskHisto set to zero.\n");
  }

  return histotot;
}

//____________________________________________________________________________ 
void AliITSQASSDDataMakerRec::ResetDetector(AliQAv1::TASKINDEX_t task)
{
  if(task==AliQAv1::kRAWS) {
    ResetRawsMonitoredObjects();
  }

}
