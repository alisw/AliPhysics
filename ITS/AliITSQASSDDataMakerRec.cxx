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
#include <TH2D.h>
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
#include "AliITSdigitSSD.h"
#include "AliITSBadChannelsSSDv2.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"

ClassImp(AliITSQASSDDataMakerRec)

AliITSQASSDDataMakerRec::AliITSQASSDDataMakerRec(AliITSQADataMakerRec *aliITSQADataMakerRec, Bool_t kMode, Int_t ldc) :
TObject(),
fAliITSQADataMakerRec(aliITSQADataMakerRec),
fSSDEvent(0),
fSSDEventPerCycle(0),
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
    fCDBManager->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    //fCDBManager->SetDefaultStorage(gSystem->Getenv("AMORE_CDB_URI"));
    Int_t runNumber = atoi(gSystem->Getenv("DATE_RUN_NUMBER"));
    if(!runNumber) 
      AliWarning("DATE_RUN_NUMBER not defined!!!\n");
    
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
      fHistSSDRawSignalModule[gHistCounterRawSignal] = new TH1D(gTitle.Data(),gTitle.Data(),
								2*fgkNumberOfPSideStrips,0,2*fgkNumberOfPSideStrips);
      gHistCounterRawSignal += 1;
      
      for(Int_t iStrip = 0; iStrip < 2*fgkNumberOfPSideStrips; iStrip++)
	fOccupancyMatrix[iModule-500][iStrip] = 0; 

      //CM histograms
      gTitle = "SSD_CM_PSide_Layer"; gTitle += gLayer;
      gTitle += "_Ladder"; gTitle += gLadder;
      gTitle += "_Module"; gTitle += gModule;
      fHistSSDCMModule[gHistCounterCM] = new TH1D(gTitle.Data(),gTitle.Data(),
						  100,-50.,50.);
      fHistSSDCMModule[gHistCounterCM]->GetXaxis()->SetTitle("CM");
      gHistCounterCM += 1;
    }//module loop

    for(Int_t iModule = 500; iModule < fgkSSDMODULES + 500; iModule++) {
      AliITSgeomTGeo::GetModuleId(iModule,gLayer,gLadder,gModule);
      gTitle = "SSD_CM_NSide_Layer"; gTitle += gLayer;
      gTitle += "_Ladder"; gTitle += gLadder;
      gTitle += "_Module"; gTitle += gModule;
      fHistSSDCMModule[gHistCounterCM] = new TH1D(gTitle.Data(),gTitle.Data(),
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

//____________________________________________________________________________ 
AliITSQASSDDataMakerRec::AliITSQASSDDataMakerRec(const AliITSQASSDDataMakerRec& qadm) :
TObject(),
fAliITSQADataMakerRec(qadm.fAliITSQADataMakerRec),
fSSDEvent(qadm.fSSDEvent),
fSSDEventPerCycle(qadm.fSSDEventPerCycle),
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

 if ( fAliITSQADataMakerRec->GetRawsData(0) == NULL ) // Raws not defined
 	return ;
 
 //for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
 //if (!AliQAv1::Instance()->IsEventSpecieSet(specie)) continue;
//cout << "StartOfDetectorCycle: Event specie " << specie << " is set" << endl;

 //Detector specific actions at start of cycle
 AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQADM::Start of SSD Cycle\n");
 //Int_t specie = fAliITSQADataMakerRec->GetEventSpecie();

 //}//event specie loop

}

//____________________________________________________________________________ 
void AliITSQASSDDataMakerRec::ResetRawsMonitoredObjects() {
  //Resetting the raw data monitored objects
  //Data size per DDL
  //for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
  //if (!AliQAv1::Instance()->IsEventSpecieSet(specie)) continue;
    /*((TH1D *)(fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+4)))->Reset();
    //Data size per LDC
    ((TH1D *)(fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+22)))->Reset();*/
  Int_t specie = fAliITSQADataMakerRec->GetEventSpecie();
  //cout << "(AliITSQASSDDataMakerRec::ResetRawsMonitoredObjects): Event specie " << specie << " is set" << endl;
  //online part
  if(fkOnline) {
    for(Int_t iModule = 500; iModule < fgkSSDMODULES + 500; iModule++) {
      for(Int_t iStrip = 0; iStrip < 2*fgkNumberOfPSideStrips; iStrip++)
	fOccupancyMatrix[iModule-500][iStrip] = 0;
    }//module loop
    
    //for(Int_t iSSDOffset = 0; iSSDOffset < fSSDRawsCommonLevelOffset; iSSDOffset++)
    //(fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+iSSDOffset))->Reset();
    //for(Int_t iSSDOffset = fSSDRawsCommonLevelOffset; iSSDOffset < fSSDRawsOffset; iSSDOffset++)
    //(fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+iSSDOffset))->Reset();

    Int_t gHistPositionOccupancyPerLadder = 0;
    Int_t gLayer = 0, gLadder = 0, gModule = 0;
    for(Int_t iModule = 0; iModule < fgkSSDMODULES; iModule++) {
      AliITSgeomTGeo::GetModuleId(iModule+500,gLayer,gLadder,gModule);
      
      gHistPositionOccupancyPerLadder = (gLayer == 5) ? 2*(gLadder - 1) : 2*(gLadder - 1 + fgkSSDLADDERSLAYER5);
      
      //P-SIDE OCCUPANCY
      fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+fgkSSDMODULES+gHistPositionOccupancyPerLadder)->Reset();
      //N-SIDE OCCUPANCY
      fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+fgkSSDMODULES+gHistPositionOccupancyPerLadder+1)->Reset();
      
      ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6))->Reset();
      ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+1))->Reset();
    }//module loop
  }//online flag
  //}//event species loop
}

//____________________________________________________________________________ 
void AliITSQASSDDataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray* /*list*/)
{

  //for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
  //if (!AliQAv1::Instance()->IsEventSpecieSet(specie)) continue;
    //if (!AliQAv1::Instance(AliQAv1::GetDetIndex(GetName()))->IsEventSpecieSet(AliRecoParam::ConvertIndex(specie))) continue;

  Int_t specie = fAliITSQADataMakerRec->GetEventSpecie();
  fAliITSQADataMakerRec->SetEventSpecie(AliRecoParam::ConvertIndex(specie));
  //cout << "(AliITSQASSDDataMakerRec::EndOfDetectorCycle): Event specie " << specie << " is set" << endl;
  
  // launch the QA checking
  AliDebug(AliQAv1::GetQADebugLevel(),"AliITSDM instantiates checker with Run(AliQAv1::kITS, task, list)\n"); 
  AliDebug(AliQAv1::GetQADebugLevel(), Form("Offset: %d\n",fGenRawsOffset[specie]));
  //Printf("Offset: %d\n",fGenRawsOffset[specie]);
  
  if(task == AliQAv1::kRAWS) {
    
    //  if (fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]) != NULL ) {
    //AliInfo(Form("Event type entries: %d - Physics events: %d"),((TH1D*)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]))->GetEntries(),((TH1D*)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]))->GetBinContent(((TH1D*)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]))->FindBin(7)));
    //cout<<"(EndOfDetectorCycle) Event type entries: "<<((TH1D*)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]))->GetEntries()<< " type 7: " << ((TH1D*)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]))->GetBinContent(((TH1D*)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]))->FindBin(7)) << endl;
    
    //AliInfo(Form("SSD Data Size entries: %d"),((TH1D*)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+1))->GetEntries());
    //cout<<"SSD Data Size entries: "<<((TH1D*)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+1))->GetEntries()<<endl;
    
    //Data size per DDL
    for(Int_t i = 0; i < fgkNumOfDDLs; i++) {
      Double_t gSizePerDDL = ((TH1D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+5+i))->GetMean();
      //cout<<"DDL: "<<i+512<<" - Size: "<<gSizePerDDL<<" - Mean: "<<
      //(fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+5+i))->GetMean()<<endl;
      //cout<<"Entries: "<<((TH1D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+5+i))->GetEntries()<<endl;
      ((TH1D *)(fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+4)))->SetBinContent(i+1,gSizePerDDL);
      //cout<<"After filling DDL: "<<i+512<<" - Size: "<<((TH1D *)(fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+4)))->GetBinContent(i+1)<<endl;
    }
    
    
    //Data size per LDC
    for(Int_t i = 0; i < fgkNumOfLDCs; i++) {
      Double_t gSizePerLDC = ((TH1D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+23+i))->GetMean();
      ((TH1D *)(fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+22)))->SetBinContent(i+1,gSizePerLDC);
      //cout<<"LDC: "<<i+170<<" - Size: "<<gSizePerLDC<<" - Mean: "<<
      //" - Size: "<<((TH1D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+23+i))->GetMean()<<endl;
    }
    
    //cout<<"Data size/ DDL entries: "<<((TH1D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+4))->GetEntries()<< " mean: " << ((TH1D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+4))->GetMean()<<endl;   
    //    cout<<"Data size/ LDC entries: "<<((TH1D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+22))->GetEntries()<< " mean: " << ((TH1D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+22))->GetMean()<<endl;
    
    //online part
    if(fkOnline) {
      //Output of the DA
      MonitorOCDBObjects();
      //Monitor common mode values
      MonitorCMValues();
      
      Int_t gHistPositionOccupancyPerModule = 0;
      Int_t gLayer = 0, gLadder = 0, gModule = 0;
      //occupancy per module
      for(Int_t iModule = 0; iModule < fgkSSDMODULES; iModule++) {
	AliITSgeomTGeo::GetModuleId(iModule+500,gLayer,gLadder,gModule);
	
	gHistPositionOccupancyPerModule = (gLayer == 5) ? ((gLadder - 1)*fgkSSDMODULESPERLADDERLAYER5 + gModule - 1) : ((gLadder - 1)*fgkSSDMODULESPERLADDERLAYER6 + gModule + fgkSSDMODULESLAYER5 - 1);
	for(Int_t iBins = 1; iBins < fHistSSDRawSignalModule[iModule]->GetXaxis()->GetNbins(); iBins++)
	  fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+gHistPositionOccupancyPerModule)->SetBinContent(iBins,fOccupancyMatrix[iModule][iBins-1]);
	
	if(fSSDEventPerCycle != 0)
	  ((TH1F *)(fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+gHistPositionOccupancyPerModule)))->Scale(100./fSSDEventPerCycle);
      }//module loop
      
      //AliInfo(Form("Entries occupancy 511/3: %d"),((TH1D*)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+222))->GetEntries());
      //  cout<<"Entries occupancy 511/3 : "<<((TH1D*)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+222))->GetEntries()<<endl;
      
      //occupancy per ladder
      Int_t gHistPositionOccupancyPerLadder = 0;
      Int_t lLadderLocationY = 0;
      Double_t occupancy = 0.0, occupancyThreshold = 0.0, occupancyAverage = 0.0;
      for(Int_t iModule = 0; iModule < fgkSSDMODULES; iModule++) {
	AliITSgeomTGeo::GetModuleId(iModule+500,gLayer,gLadder,gModule);
	
	gHistPositionOccupancyPerModule = (gLayer == 5) ? ((gLadder - 1)*fgkSSDMODULESPERLADDERLAYER5 + gModule - 1) : ((gLadder - 1)*fgkSSDMODULESPERLADDERLAYER6 + gModule + fgkSSDMODULESLAYER5 - 1);
	gHistPositionOccupancyPerLadder = (gLayer == 5) ? 2*(gLadder - 1) : 2*(gLadder - 1 + fgkSSDLADDERSLAYER5);
	
	//P-SIDE OCCUPANCY
	occupancy = GetOccupancyModule((TH1 *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+gHistPositionOccupancyPerModule),0,0,0);
	occupancyThreshold = GetOccupancyModule((TH1 *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+gHistPositionOccupancyPerModule),0,1,3);
	occupancyAverage = GetOccupancyModule((TH1 *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+gHistPositionOccupancyPerModule),0,2,0);
	
	fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+fgkSSDMODULES+gHistPositionOccupancyPerLadder)->SetBinContent(gModule,occupancy);
	lLadderLocationY = 3*gLadder; // sideP=1 sideN=0 
	if(gLayer == 5) {
	  //occupancy per module - no threshold
	  ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6))->SetBinContent(gModule,lLadderLocationY,occupancy);
	  //occupancy per module - threshold @ 3%
	  ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+2))->SetBinContent(gModule,lLadderLocationY,occupancyThreshold);
	  //average occupancy per module
	  ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+4))->SetBinContent(gModule,lLadderLocationY,occupancyAverage);
	}
	else if(gLayer == 6) {
	  //occupancy per module - no threshold
	  ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+1))->SetBinContent(gModule,lLadderLocationY,occupancy);
	  //occupancy per module - threshold @ 3%
	  ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+3))->SetBinContent(gModule,lLadderLocationY,occupancyThreshold);
	  //average occupancy per module
	  ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+5))->SetBinContent(gModule,lLadderLocationY,occupancyAverage);
	}
	
	//N-SIDE OCCUPANCY
	occupancy = GetOccupancyModule((TH1 *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+gHistPositionOccupancyPerModule),1,0,0);   
	occupancyThreshold = GetOccupancyModule((TH1 *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+gHistPositionOccupancyPerModule),1,1,3);   
	occupancyAverage = GetOccupancyModule((TH1 *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+gHistPositionOccupancyPerModule),1,2,0);   
	
	fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+fgkSSDMODULES+gHistPositionOccupancyPerLadder+1)->SetBinContent(gModule,occupancy);
	if(gLayer == 5) {
	  //occupancy per module - no threshold
	  ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6))->SetBinContent(gModule,lLadderLocationY-1,occupancy);
	  //occupancy per module - threshold @ 3%
	  ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+2))->SetBinContent(gModule,lLadderLocationY-1,occupancyThreshold);
	  //average occupancy per module
	  ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+4))->SetBinContent(gModule,lLadderLocationY-1,occupancyAverage);
	}
	else if(gLayer == 6) {
	  //occupancy per module - no threshold
	  ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+1))->SetBinContent(gModule,lLadderLocationY-1,occupancy);
	  //occupancy per module - threshold @ 3%
	  ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+3))->SetBinContent(gModule,lLadderLocationY-1,occupancyThreshold);
	  //average occupancy per module
	  ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+5))->SetBinContent(gModule,lLadderLocationY-1,occupancyAverage);
	}
      }//module loop
    }//online flag for SSD
    
    //AliInfo(Form("Entries 2d occupancy no thres- lay 5: %d"),((TH2D*)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6))->GetEntries());
    //cout<<"entries 2d occupancy thres- lay 6: "<<((TH2D*)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+3))->GetEntries()<< " mean: " << ((TH2D*)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+3))->GetMean() << endl; //Somehow the other occupancy maps do give nonzero values for GetMean() here
    
    fSSDEventPerCycle = 0;
    
    //cout<<"Data size/ DDL entries: "<<((TH1D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[1]+4))->GetEntries()<< " mean: " << ((TH1D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[1]+4))->GetMean()<<endl;   
    
    //Reset of the raws
    ResetRawsMonitoredObjects();
    
    //AliQAChecker::Instance()->Run( AliQAv1::kITS , task, list);
  }//raw data end of cycle
  
  //} //event specie loop
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
  //cout << "(AliITSQASSDDataMakerRec::InitRaws): Offset " << fGenRawsOffset[specie] << endl;
    
    
  if(fkOnline) {
    AliDebug(AliQAv1::GetQADebugLevel(), "Book Online Histograms for SSD\n");
  }
  else {
    AliDebug(AliQAv1::GetQADebugLevel(), "Book Offline Histograms for SSD\n ");
  }
  AliDebug(AliQAv1::GetQADebugLevel(), Form("Number of histograms (SPD+SDD): %d\n",fGenRawsOffset[specie]));
  
  TString gTitle;
  TString gName;
  //book online-offline QA histos
  TH1D *fHistSSDEventType = new TH1D("fHistSSDEventType",
				     "SSD Event Type;Event type;Events",
				     31,-1,30);
  rv = fAliITSQADataMakerRec->Add2RawsList((new TH1D(*fHistSSDEventType)), 
					   fGenRawsOffset[specie]+fSSDRawsOffset, expert, !image, !saveCorr);
  fSSDRawsOffset += 1;
  delete fHistSSDEventType;
  //cout<<"(AliITSQASSDDataMakerRec::InitRaws): "<<fGenRawsOffset[specie]+fSSDRawsOffset-1<<" - Name: "<<(fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]))->GetName()<<endl;
  TH1D *fHistSSDDataSize = new TH1D("fHistSSDDataSize",
				    "SSD Data Size;(SSD data size) [KB];Events",
				    1000,0,500);
  rv = fAliITSQADataMakerRec->Add2RawsList((new TH1D(*fHistSSDDataSize)), 
					   fGenRawsOffset[specie]+fSSDRawsOffset, expert, !image, !saveCorr);
  delete fHistSSDDataSize;
  fSSDRawsOffset += 1;
  TH1D *fHistSSDDataSizePercentage = new TH1D("fHistSSDDataSizePercentage",
					      "SSD Data Size Percentage;SSD data size [%];Events",
					      1000,0,100);
  rv = fAliITSQADataMakerRec->Add2RawsList((new TH1D(*fHistSSDDataSizePercentage)), 
					   fGenRawsOffset[specie]+fSSDRawsOffset, expert, !image, !saveCorr);
  delete fHistSSDDataSizePercentage;
  fSSDRawsOffset += 1;
  TH1D *fHistSSDDDLId = new TH1D("fHistSSDDDLId",
				 "SSD DDL Id;DDL id;Events",16,511.5,527.5);
  rv = fAliITSQADataMakerRec->Add2RawsList((new TH1D(*fHistSSDDDLId)), 
					   fGenRawsOffset[specie]+fSSDRawsOffset, expert, !image, !saveCorr);
  delete fHistSSDDDLId;
  fSSDRawsOffset += 1;
  TH1D *fHistSSDDataSizePerDDL = new TH1D("fHistSSDDataSizePerDDL",
					  "SSD Data Size Per DDL;DDL id;<SSD data size> [KB]",
					  16,511.5,527.5);
  rv = fAliITSQADataMakerRec->Add2RawsList((new TH1D(*fHistSSDDataSizePerDDL)), 
					   fGenRawsOffset[specie]+fSSDRawsOffset, !expert, image, !saveCorr);
  delete fHistSSDDataSizePerDDL;
  fSSDRawsOffset += 1;
  TH1D *fHistSSDDataSizeDDL[fgkNumOfDDLs];
  for(Int_t i = 1; i < fgkNumOfDDLs+1; i++) {
    gName = Form("fHistSSDDataSizeDDL%d", i+511) ;
    gTitle = Form("SSD Data Size DDL %d", i+511) ;
    fHistSSDDataSizeDDL[i-1] = new TH1D(gName.Data(),
					Form("%s;(SSD data size) [KB];Events", gTitle.Data()),
					100,0,50);
    rv = fAliITSQADataMakerRec->Add2RawsList((new TH1D(*fHistSSDDataSizeDDL[i-1])), 
					     fGenRawsOffset[specie]+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1;
  }
  for(Int_t i = 1; i < fgkNumOfDDLs+1; i++) delete fHistSSDDataSizeDDL[i-1];
  
  TH1D *fHistSSDLDCId = new TH1D("fHistSSDLDCId","SSD LDC Id;LDC id;Events",8,169.5,177.5);
  rv = fAliITSQADataMakerRec->Add2RawsList((new TH1D(*fHistSSDLDCId)), 
					   fGenRawsOffset[specie]+fSSDRawsOffset, expert, !image, !saveCorr);
  delete fHistSSDLDCId;
  fSSDRawsOffset += 1;
  TH1D *fHistSSDDataSizePerLDC = new TH1D("fHistSSDDataSizePerLDC",
					  "SSD Data Size Per LDC;LDC id;<SSD data size> [KB]",
					  8,169.5,177.5);
  rv = fAliITSQADataMakerRec->Add2RawsList((new TH1D(*fHistSSDDataSizePerLDC)), 
					   fGenRawsOffset[specie]+fSSDRawsOffset, !expert, image, !saveCorr);
  delete fHistSSDDataSizePerLDC;
  fSSDRawsOffset += 1;
  TH1D *fHistSSDDataSizeLDC[fgkNumOfLDCs];
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
    fHistSSDDataSizeLDC[i-1] = new TH1D(gName.Data(),
					Form("%s;SSD data size [KB];Events", gTitle.Data()),
					1000,0,100);
    rv = fAliITSQADataMakerRec->Add2RawsList((new TH1D(*fHistSSDDataSizeLDC[i-1])), 
					     fGenRawsOffset[specie]+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1;
  }
  for(Int_t i = 1; i < fgkNumOfLDCs+1; i++) delete fHistSSDDataSizeLDC[i-1];
  fSSDRawsCommonLevelOffset = fSSDRawsOffset;
  
  if(fkOnline) {
    Int_t gLayer = 0, gLadder = 0, gModule = 0;
    //occupancy per SSD module
    TH1D *fHistSSDOccupancyModule[fgkSSDMODULES]; 
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
      
      fHistSSDOccupancyModule[i-500] = new TH1D(gName.Data(),Form("%s;N_{strip};Occupancy [%]", gTitle.Data()),
						2*fgkNumberOfPSideStrips,0,2*fgkNumberOfPSideStrips);
      fHistSSDOccupancyModule[i-500]->GetXaxis()->SetTitleColor(1);
      rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDOccupancyModule[i-500], 
					       fGenRawsOffset[specie]+fSSDRawsOffset, expert, !image, !saveCorr);
      fSSDRawsOffset += 1;
    }
    
    //Occupancy per SSD ladder
    Int_t occupancyCounter = 0;
    TH1D *fHistSSDOccupancyLadder[2*(fgkSSDLADDERSLAYER5 + fgkSSDLADDERSLAYER6)];
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
        fHistSSDOccupancyLadder[occupancyCounter] = new TH1D(gName.Data(),
                                                             Form("%s;Module number;Occupancy [%]", gTitle.Data()),
                                                             AliITSgeomTGeo::GetNDetectors(iLayer),
                                                             0.5,AliITSgeomTGeo::GetNDetectors(iLayer)+0.5);
        fHistSSDOccupancyLadder[occupancyCounter]->GetXaxis()->SetTitleColor(1);
        rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDOccupancyLadder[occupancyCounter], 
                                            fGenRawsOffset[specie]+fSSDRawsOffset, expert, !image, !saveCorr);
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
        fHistSSDOccupancyLadder[occupancyCounter] = new TH1D(gName.Data(),
                                                             Form("%s;Module number;Occupancy [%]", gTitle.Data()),
                                                             AliITSgeomTGeo::GetNDetectors(iLayer),
                                                             0.5,AliITSgeomTGeo::GetNDetectors(iLayer)+0.5);
        fHistSSDOccupancyLadder[occupancyCounter]->GetXaxis()->SetTitleColor(1);
        rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDOccupancyLadder[occupancyCounter], 
                                            fGenRawsOffset[specie]+fSSDRawsOffset, expert, !image, !saveCorr);
        occupancyCounter += 1; fSSDRawsOffset += 1;
      }//ladder loop
    }//layer loop

    //top level occupancy plots
    //occupancy per module - no threshold
    TH2D *fHistSSDOccupancyLayer5 = new TH2D("fHistSSDOccupancyLayer5",
					     "SSD Occupancy (Layer 5) - No threshold;N_{modules};N_{Ladders}",
					     fgkSSDMODULESPERLADDERLAYER5,
					     0,fgkSSDMODULESPERLADDERLAYER5,
					     3*fgkSSDLADDERSLAYER5,
					     500,500+fgkSSDLADDERSLAYER5);  
    fHistSSDOccupancyLayer5->GetZaxis()->SetRangeUser(0.0,100.0);
    Char_t fLabel[3];
    for(Int_t iBin = 1; iBin < fgkSSDMODULESPERLADDERLAYER5 + 1; iBin++){
      sprintf(fLabel,"%d",iBin);
      fHistSSDOccupancyLayer5->GetXaxis()->SetBinLabel(iBin,fLabel);
    }
    rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDOccupancyLayer5, 
					fGenRawsOffset[specie]+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1;
    TH2D *fHistSSDOccupancyLayer6 = new TH2D("fHistSSDOccupancyLayer6",
					     "Occupancy per module (Layer 6) - No threshold;N_{modules};N_{Ladders}",
					     fgkSSDMODULESPERLADDERLAYER6,
					     0,fgkSSDMODULESPERLADDERLAYER6,
					     3*fgkSSDLADDERSLAYER6,
					     600,600+fgkSSDLADDERSLAYER6);
    fHistSSDOccupancyLayer6->GetZaxis()->SetRangeUser(0.0,100.0);
    for(Int_t iBin = 1; iBin < fgkSSDMODULESPERLADDERLAYER6 + 1; iBin++){
      sprintf(fLabel,"%d",iBin);
      fHistSSDOccupancyLayer6->GetXaxis()->SetBinLabel(iBin,fLabel);
    }
    rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDOccupancyLayer6, 
                                        fGenRawsOffset[specie]+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1;

    //occupancy per module - threshold @ 3%
    TH2D *fHistSSDOccupancyThresholdLayer5 = new TH2D("fHistSSDOccupancyThresholdLayer5",
						      "Occupancy per module (Layer 5) - Threshold 3%;N_{modules};N_{Ladders};Entries",
						      fgkSSDMODULESPERLADDERLAYER5,
						      0,fgkSSDMODULESPERLADDERLAYER5,
						      3*fgkSSDLADDERSLAYER5,
						      500,500+fgkSSDLADDERSLAYER5);  
    fHistSSDOccupancyThresholdLayer5->GetZaxis()->SetRangeUser(3.0,10.0);
    for(Int_t iBin = 1; iBin < fgkSSDMODULESPERLADDERLAYER5 + 1; iBin++){
      sprintf(fLabel,"%d",iBin);
      fHistSSDOccupancyThresholdLayer5->GetXaxis()->SetBinLabel(iBin,fLabel);
    }
    rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDOccupancyThresholdLayer5, 
                                        fGenRawsOffset[specie]+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1;
    TH2D *fHistSSDOccupancyThresholdLayer6 = new TH2D("fHistSSDOccupancyThresholdLayer6",
						      "Occupancy per module (Layer 6) - Threshold 3%;N_{modules};N_{Ladders}",
						      fgkSSDMODULESPERLADDERLAYER6,
						      0,fgkSSDMODULESPERLADDERLAYER6,
						      3*fgkSSDLADDERSLAYER6,
						      600,600+fgkSSDLADDERSLAYER6);
    fHistSSDOccupancyThresholdLayer6->GetZaxis()->SetRangeUser(3.0,10.0);
    for(Int_t iBin = 1; iBin < fgkSSDMODULESPERLADDERLAYER6 + 1; iBin++){
      sprintf(fLabel,"%d",iBin);
      fHistSSDOccupancyThresholdLayer6->GetXaxis()->SetBinLabel(iBin,fLabel);
    }
    rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDOccupancyThresholdLayer6, 
					fGenRawsOffset[specie]+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1;

    //Average occupancy per module
    TH2D *fHistSSDAverageOccupancyLayer5 = new TH2D("fHistSSDAverageOccupancyLayer5",
						    "Average occupancy per module (Layer 5);N_{modules};N_{Ladders}",
						    fgkSSDMODULESPERLADDERLAYER5,
						    0,fgkSSDMODULESPERLADDERLAYER5,
						    3*fgkSSDLADDERSLAYER5,
						    500,500+fgkSSDLADDERSLAYER5);  
    fHistSSDAverageOccupancyLayer5->GetZaxis()->SetRangeUser(0.0,5.0);
    for(Int_t iBin = 1; iBin < fgkSSDMODULESPERLADDERLAYER5 + 1; iBin++){
      sprintf(fLabel,"%d",iBin);
      fHistSSDAverageOccupancyLayer5->GetXaxis()->SetBinLabel(iBin,fLabel);
    }
    rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDAverageOccupancyLayer5, 
					     fGenRawsOffset[specie]+fSSDRawsOffset, !expert, image, !saveCorr);
    fSSDRawsOffset += 1;
    TH2D *fHistSSDAverageOccupancyLayer6 = new TH2D("fHistSSDAverageOccupancyLayer6",
						    "Average occupancy per module (Layer 6);N_{modules};N_{Ladders}",
						    fgkSSDMODULESPERLADDERLAYER6,
						    0,fgkSSDMODULESPERLADDERLAYER6,
						    3*fgkSSDLADDERSLAYER6,
						    600,600+fgkSSDLADDERSLAYER6);
    fHistSSDAverageOccupancyLayer6->GetZaxis()->SetRangeUser(0.0,5.0);
    for(Int_t iBin = 1; iBin < fgkSSDMODULESPERLADDERLAYER6 + 1; iBin++){
      sprintf(fLabel,"%d",iBin);
      fHistSSDAverageOccupancyLayer6->GetXaxis()->SetBinLabel(iBin,fLabel);
    }
    rv = fAliITSQADataMakerRec->Add2RawsList(fHistSSDAverageOccupancyLayer6, 
					     fGenRawsOffset[specie]+fSSDRawsOffset, !expert, image, !saveCorr);
    fSSDRawsOffset += 1;

    //Output of the DA
    TH2D *fHistSSDPSideBadChannelMapLayer5 = new TH2D("fHistSSDPSideBadChannelMapLayer5",
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
					     fGenRawsOffset[specie]+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1; fSSDRawsDAOffset += 1;
    
    TH2D *fHistSSDNSideBadChannelMapLayer5 = new TH2D("fHistSSDNSideBadChannelMapLayer5",
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
					     fGenRawsOffset[specie]+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1; fSSDRawsDAOffset += 1;
    
    TH2D *fHistSSDPSideBadChannelMapLayer6 = new TH2D("fHistSSDPSideBadChannelMapLayer6",
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
					     fGenRawsOffset[specie]+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1; fSSDRawsDAOffset += 1;
    
    TH2D *fHistSSDNSideBadChannelMapLayer6 = new TH2D("fHistSSDNSideBadChannelMapLayer6",
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
                                        fGenRawsOffset[specie]+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1; fSSDRawsDAOffset += 1;

   //Common mode values
    TH2D *fHistSSDPSideCommonModeMapLayer5 = new TH2D("fHistSSDPSideCommonModeMapLayer5",
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
					     fGenRawsOffset[specie]+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1; fSSDRawsDAOffset += 1;
    
    TH2D *fHistSSDNSideCommonModeMapLayer5 = new TH2D("fHistSSDNSideCommonModeMapLayer5",
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
					     fGenRawsOffset[specie]+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1; fSSDRawsDAOffset += 1;
    
    TH2D *fHistSSDPSideCommonModeMapLayer6 = new TH2D("fHistSSDPSideCommonModeMapLayer6",
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
					     fGenRawsOffset[specie]+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1; fSSDRawsDAOffset += 1;
    
    TH2D *fHistSSDNSideCommonModeMapLayer6 = new TH2D("fHistSSDNSideCommonModeMapLayer6",
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
					     fGenRawsOffset[specie]+fSSDRawsOffset, expert, !image, !saveCorr);
    fSSDRawsOffset += 1; fSSDRawsDAOffset += 1;
  }//online flag
  
  fSSDhRawsTask = fSSDRawsOffset;
  AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SSD Raws histograms booked\n",fSSDhRawsTask));
  AliDebug(AliQAv1::GetQADebugLevel(), Form("Number of histograms (SPD+SDD+SSD): %d\n",fGenRawsOffset[specie]+fSSDhRawsTask));  
  AliDebug(AliQAv1::GetQADebugLevel(),Form("Number of histograms (SPD+SDD+SSD): %d\n",fGenRawsOffset[specie]+fSSDRawsOffset));
  
  /*
   fSSDhTask = fSSDRawsOffset;
   AliDebug(AliQAv1::GetQADebugLevel(),Form("%d SSD Raws histograms booked\n",fSSDhTask));
   AliDebug(AliQAv1::GetQADebugLevel(), Form("Number of histograms (SPD+SDD+SSD): %d\n",fGenRawsOffset[specie]+fSSDhTask));  
   AliDebug(AliQAv1::GetQADebugLevel(),Form("Number of histograms (SPD+SDD+SSD): %d\n",fGenRawsOffset[specie]+fSSDRawsOffset));
   */


  //} //event species loop

  return rv ; 


}

//____________________________________________________________________________
Int_t AliITSQASSDDataMakerRec::MakeRaws(AliRawReader* rawReader) { 
  // Fill QA for RAW - SSD -
  Int_t rv = 0 ; 

  Int_t gStripNumber;
  Int_t gHistPosition;
  Int_t gLayer = 0,gLadder = 0, gModule = 0;

  Double_t gSizePerDDL[fgkNumOfDDLs] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t gSizePerLDC[fgkNumOfLDCs] = {0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t sumSSDDataSize = 0.0;
  Double_t eventSize = -1.0;

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    if (!AliQAv1::Instance()->IsEventSpecieSet(specie)) continue;
    //cout << "MakeRaws: Event specie " << specie << " is set" << endl;

    //AliInfo(Form("fGenRawsOffset[specie] %d\n",fGenRawsOffset[specie]));
    if(fkOnline) {
      //reset the signal vs strip number histograms
      for(Int_t iModule = 0; iModule < fgkSSDMODULES; iModule++)
	fHistSSDRawSignalModule[iModule]->Reset();
    }//online flag
    
    rawReader->Select("ITSSSD",-1,-1);  
    rawReader->Reset(); //rawReader->NextEvent();   
    //AliInfo(Form("fGenRawsOffset[specie] %d\n",fGenRawsOffset[specie]));
    (fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]))->Fill(rawReader->GetType());
    
    if(rawReader->GetType() == 7) {
      fSSDEvent += 1;
      fSSDEventPerCycle += 1;
    }
    
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
      gSizePerDDL[rawReader->GetDDLID()] = rawReader->GetDataSize();
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
	  //fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+gHistPosition+fSSDRawsCommonLevelOffset)->Fill(gStripNumber,gSSDStream.GetSignal());
	}//normal strip signal
      }//online flag
    }//streamer loop   
    
    //event size calculation and filling info
    for(Int_t i = 0; i < fgkNumOfDDLs; i++) {
      sumSSDDataSize += gSizePerDDL[i];
      if(gSizePerDDL[i] > 0) {
	(fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+3))->Fill(i+512);
	(fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+5+i))->Fill(gSizePerDDL[i]/1e+03);
	
	
	//if(i == 5)
	//cout<<gSizePerDDL[i]/1e+03<<endl;
	//cout<<"Event: "<<fSSDEventPerCycle<<" - DDL: "<<i+512<<
	//" - Data size: "<<gSizePerDDL[i]/1e+03<<endl;
      }
      //(fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+4))->Fill(i+512,gSizePerDDL[i]/1e+06);
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
      
      if(gSizePerLDC[i] > 0)
	(fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+21))->Fill(i+170);
      (fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+23+i))->Fill(gSizePerLDC[i]/1e+03);
      //cout<<"Event: "<<fSSDEventPerCycle<<" - LDC: "<<i+170<<
      //" - Data size: "<<gSizePerLDC[i]<<endl;
      
      //(fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+22))->Fill(i+6,gSizePerLDC[i]/1e+06);
    }
    
    //    for(Int_t i = 0; i < fgkNumOfDDLs; i++) {
    //      Double_t gSizePerDDL = ((TH1D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+5+i))->GetMean();
    //if(i == 5)
      //cout<<"DDL: "<<i+512<<" - Size: "<<gSizePerDDL<<
      //" - Mean: "<<((TH1D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+5+i))->GetMean()<<endl;
    //    }
    
    if(sumSSDDataSize) 
      (fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+1))->Fill(sumSSDDataSize/1e+03);
    if(eventSize)
      (fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+2))->Fill(100.*sumSSDDataSize/eventSize);
    
    //Occupancy calculation
    if(fkOnline) {
      for(Int_t iModule = 0; iModule < fgkSSDMODULES; iModule++) {
	GetOccupancyStrip(fHistSSDRawSignalModule[iModule],fOccupancyMatrix[iModule]);
	//if(iModule == 156) {    
	//cout<<"========================================"<<endl;
	//AliITSgeomTGeo::GetModuleId(656,gLayer,gLadder,gModule);
	/*for(Int_t iBin = 1; iBin < fHistSSDRawSignalModule[iModule]->GetXaxis()->GetNbins(); iBin++) {
	  if((iBin >= 750)&&(iBin <= 780))
	  cout<<"Event: "<<fSSDEventPerCycle<<" - Ladder: "<<gLadder+499<<
	  " - Module: "<<gModule<<" - Strip: "<<iBin<<
	  " - Signal: "<<fHistSSDRawSignalModule[iModule]->GetBinContent(iBin)<<endl;
	  }*///strip loop --> to be removed
	//}//module cut --> to be removed
      }//module loop
    }//online flag for SSD
    
    //cout<<"Event type entries: "<<((TH1*)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]))->GetEntries()<<endl;
    //cout<<"DDL id entries at MR: "<<((TH1D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+3))->GetEntries()<< endl;    
    //cout<<"LDC id entries at MR: "<<((TH1D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+21))->GetEntries()<< endl; 
    
  } //event species loop
  
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
 void AliITSQASSDDataMakerRec::MonitorCMValues() {
  //Monitor in AMORE the CM values
   //for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
   //if (!AliQAv1::Instance()->IsEventSpecieSet(specie)) continue;
     //cout << "MonitorCMValues: Event specie " << specie << " is set" << endl;
   Int_t specie = fAliITSQADataMakerRec->GetEventSpecie();
   fAliITSQADataMakerRec->SetEventSpecie(AliRecoParam::ConvertIndex(specie));
   
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
       ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsOffset-fSSDRawsDAOffset+4))->SetBinContent(gModule,gLadder,rmsPsideCM);
       ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsOffset-fSSDRawsDAOffset+5))->SetBinContent(gModule,gLadder,rmsNsideCM);
     }
     if(gLayer == 6) {
       ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsOffset-fSSDRawsDAOffset+6))->SetBinContent(gModule,gLadder,rmsPsideCM);
       ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsOffset-fSSDRawsDAOffset+7))->SetBinContent(gModule,gLadder,rmsNsideCM);
     }
   }//module loop
   //}//event species loop
 }

//____________________________________________________________________________
void AliITSQASSDDataMakerRec::MonitorOCDBObjects() { 
  //Monitor in AMORE the output of the DA
  //Currently only the bad channel list is monitored
  //Todo: Noise - Pedestal


  //for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
  //if (!AliQAv1::Instance()->IsEventSpecieSet(specie)) continue;
  //cout << "MonitorOCDBObjects: Event specie " << specie << " is set" << endl;
  Int_t specie = fAliITSQADataMakerRec->GetEventSpecie();
  fAliITSQADataMakerRec->SetEventSpecie(AliRecoParam::ConvertIndex(specie));
  
  //((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsOffset-fSSDRawsDAOffset))->Reset();
  //((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsOffset-fSSDRawsDAOffset+1))->Reset();
  //((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsOffset-fSSDRawsDAOffset+2))->Reset();
  //((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsOffset-fSSDRawsDAOffset+3))->Reset();
  
  AliCDBEntry *entryBadChannelsSSD = fCDBManager->Get("ITS/Calib/BadChannelsSSD");
  if(!entryBadChannelsSSD) 
    AliError("OCDB entry for the bad channel list is not valid!"); 
  AliITSBadChannelsSSDv2 *badchannelsSSD = (AliITSBadChannelsSSDv2 *)entryBadChannelsSSD->GetObject();
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


    //cout << "Bad channels P side module " << module << ": " << nBadPSideChannels << endl;
    //cout << "Bad channels N side module " << module << ": " << nBadNSideChannels << endl;
    
    if(layer == 5) {
      /*if((module == 10)&&(ladder == 10)) {
	cout<<"Npside bad: "<<nPSideChannelsLayer5<<" - Total: "<<fgkNumberOfPSideStrips<<" - Percentage: "<<(100.*nPSideChannelsLayer5/fgkNumberOfPSideStrips)<<endl;
	}*/
      if(nPSideChannelsLayer5 > 0)
	((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsOffset-fSSDRawsDAOffset))->SetBinContent(module,ladder,100.*nPSideChannelsLayer5/fgkNumberOfPSideStrips);
      else ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsOffset-fSSDRawsDAOffset))->SetBinContent(module,ladder,0.0001);
      if(nNSideChannelsLayer5 > 0)
	((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsOffset-fSSDRawsDAOffset+1))->SetBinContent(module,ladder,100.*nNSideChannelsLayer5/fgkNumberOfPSideStrips);
      else ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsOffset-fSSDRawsDAOffset+1))->SetBinContent(module,ladder,0.0001);
    }//layer 5                                                                    
    if(layer == 6) {
      if(nPSideChannelsLayer6 > 0)
        ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsOffset-fSSDRawsDAOffset+2))->SetBinContent(module,ladder,100.*nPSideChannelsLayer6/fgkNumberOfPSideStrips);
      else ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsOffset-fSSDRawsDAOffset+2))->SetBinContent(module,ladder,0.0001);
      if(nNSideChannelsLayer6 > 0)
        ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsOffset-fSSDRawsDAOffset+3))->SetBinContent(module,ladder,100.*nNSideChannelsLayer6/fgkNumberOfPSideStrips);
      else ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsOffset-fSSDRawsDAOffset+3))->SetBinContent(module,ladder,0.0001);
    }//layer 6                                                              
  }//module loop

  //cout << "entries bad channel layer 5 n side " << ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsOffset-fSSDRawsDAOffset+1))->GetEntries() << " - Bad channels P side layer 5 module 10 ladder 10: " << ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenRawsOffset[specie]+fSSDRawsOffset-fSSDRawsDAOffset))->GetBinContent(10,10)<<endl;


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
Int_t AliITSQASSDDataMakerRec::MakeDigits(TTree *digits) { 
  // Fill QA for DIGIT - SSD -
//  AliITS *fITS  = (AliITS*)gAlice->GetModule("ITS");
//  fITS->SetTreeAddress();
//  TClonesArray *iSSDdigits  = fITS->DigitsAddress(2);
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
    fAliITSQADataMakerRec->GetDigitsData(fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 0)->Fill(iModule,ndigits);
    if(ndigits != 0)
      AliDebug(AliQAv1::GetQADebugLevel(),Form("Module: %d - Digits: %d",iModule,ndigits));
    
    for (Int_t iDigit = 0; iDigit < ndigits; iDigit++) {
      AliITSdigit *dig = (AliITSdigit*)iSSDdigits->UncheckedAt(iDigit);
      Int_t fStripNumber = (dig->GetCoord1() == 0) ? dig->GetCoord2() : dig->GetCoord2() + fgkNumberOfPSideStrips;
      ((TH2F *)fAliITSQADataMakerRec->GetDigitsData(fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 1))->Fill(fStripNumber,iModule,dig->GetSignal());
    }//digit loop
  }//module loop
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
  Int_t rv = 0 ; 
  Int_t gLayer = 0, gLadder = 0, gModule = 0;
  Int_t lLadderLocationY = 0;
  TBranch *branchRecP = clustersTree->GetBranch("ITSRecPoints");
  if (!branchRecP) { 
    AliError("can't get the branch with the ITS clusters !");
    return rv;
  }
 
  //AliInfo(Form("fAliITSQADataMakerRec->GetEventSpecie() %d\n",fAliITSQADataMakerRec->GetEventSpecie()));
  //AliInfo(Form("fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] %d\n",fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()]));
  static TClonesArray statRecpoints("AliITSRecPoint");
  TClonesArray *recpoints = &statRecpoints;
  branchRecP->SetAddress(&recpoints);
  Int_t nClustersLayer5 = 0, nClustersLayer6 = 0;
  Int_t npoints = 0;      
  Float_t cluglo[3]={0.,0.,0.}; 
  //printf("*-*-*-*-*-*-*---*-*-*-------*-*-*-*-*-*-***************AliITSQASSDataMakerRec::MakeRecpoints STEP1 \n");
  for(Int_t module = 0; module < clustersTree->GetEntries(); module++){
    branchRecP->GetEvent(module);
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
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 20)->Fill(recp->GetType());

	if(recp->GetType() != 1) continue;
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 0)->Fill(module);
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 4)->Fill(recp->GetDetLocalX());
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 6)->Fill(recp->GetDetLocalZ());
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 8)->Fill(cluglo[0]);
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 10)->Fill(cluglo[1]);
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 12)->Fill(cluglo[2]);
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 14)->Fill(phi);
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 16)->Fill(theta);
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 18)->Fill(radius);
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 22)->Fill(recp->GetChargeRatio());
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 24)->Fill(recp->GetQ());
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 26)->Fill(chargePSide);
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 28)->Fill(chargeNSide);
	if(chargePSide != 0.) fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 30)->Fill(chargeNSide/chargePSide);
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 32)->Fill(chargePSide,chargeNSide);
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 34)->SetBinContent(gModule,lLadderLocationY,recp->GetQ());
	((TH2F *)fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 36))->Fill(gModule,499+gLadder,1);
	nClustersLayer5 += 1;
      }//layer 5 histograms
      if(layer == 5) {
	//printf("-*-*-*-*-*-***************AliITSQASSDataMakerRec::MakeRecpoints Filling 5 called \n");
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 21)->Fill(recp->GetType());

	if(recp->GetType() != 1) continue;
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 1)->Fill(module);
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 5)->Fill(recp->GetDetLocalX());
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 7)->Fill(recp->GetDetLocalZ());
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 9)->Fill(cluglo[0]);
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 11)->Fill(cluglo[1]);
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 13)->Fill(cluglo[2]);
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 15)->Fill(phi);
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 17)->Fill(theta);
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 19)->Fill(radius);
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 23)->Fill(recp->GetChargeRatio());
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 25)->Fill(recp->GetQ());
        fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 27)->Fill(chargePSide);
        fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 29)->Fill(chargeNSide);
        if(chargePSide != 0.) fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 31)->Fill(chargeNSide/chargePSide);
        fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 33)->Fill(chargePSide,chargeNSide);
	fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 35)->SetBinContent(gModule,lLadderLocationY,recp->GetQ());
	((TH2F *)fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 37))->Fill(gModule,599+gLadder,1);
	nClustersLayer6 += 1;
      }//layer 6 histograms
    }//rec. points loop
  }//module loop
  
  //printf("-*-*-*-*-*-***************AliITSQASSDataMakerRec::MakeRecpoints Filling called \n");
  fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 2)->Fill(nClustersLayer5);
  fAliITSQADataMakerRec->GetRecPointsData(fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()] + 3)->Fill(nClustersLayer6);

  statRecpoints.Clear();
  return rv ; 
}

//____________________________________________________________________________ 
Int_t AliITSQASSDDataMakerRec::GetOffset(AliQAv1::TASKINDEX_t task) {
  // Returns offset number according to the specified task 
  Int_t offset=0;
  if( task == AliQAv1::kRAWS ) {
    offset=fGenRawsOffset[fAliITSQADataMakerRec->GetEventSpecie()];  
  }
  else if( task == AliQAv1::kDIGITSR ) {
    offset=fGenDigitsOffset[fAliITSQADataMakerRec->GetEventSpecie()];   
  }
  else if( task == AliQAv1::kRECPOINTS ) {
    offset=fGenRecPointsOffset[fAliITSQADataMakerRec->GetEventSpecie()];   
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
