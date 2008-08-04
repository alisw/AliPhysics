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
#include "AliQA.h"
#include "AliQAChecker.h"
#include "AliRawReader.h"
#include "AliRawReaderRoot.h"
#include "AliITSRawStreamSSD.h"
#include "AliITSgeomTGeo.h"
#include "AliRawEventHeaderBase.h"
#include "AliITSRecPoint.h"
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
fSSDhTask(0),
fGenOffset(0),
fCDBManager(0) {
  // Default constructor   
  //initilize the raw signal vs strip number histograms
  if(fkOnline) {
    fCDBManager = AliCDBManager::Instance();
    //fCDBManager->SetDefaultStorage("local://$ALICE_ROOT");
    fCDBManager->SetDefaultStorage(gSystem->Getenv("AMORE_CDB_URI"));
    Int_t runNumber = atoi(gSystem->Getenv("DATE_RUN_NUMBER"));
    if(!runNumber) 
      AliInfo("DATE_RUN_NUMBER not defined!!!\n");
    
    fCDBManager->SetRun(runNumber);
    AliCDBEntry *geomGRP = fCDBManager->Get("GRP/Geometry/Data");
    if(!geomGRP) AliInfo("GRP geometry not found!!!\n");
    
    Int_t gLayer = 0,gLadder = 0, gModule = 0;
    Int_t gHistCounter = 0;
    TString gTitle; 

    for(Int_t iModule = 500; iModule < fgkSSDMODULES + 500; iModule++) {
      AliITSgeomTGeo::GetModuleId(iModule,gLayer,gLadder,gModule);
      gTitle = "SSD_RawSignal_Layer"; gTitle += gLayer;
      gTitle += "_Ladder"; gTitle += gLadder;
      gTitle += "_Module"; gTitle += gModule;
      fHistSSDRawSignalModule[gHistCounter] = new TH1D(gTitle.Data(),gTitle.Data(),
						       2*fgkNumberOfPSideStrips,0,2*fgkNumberOfPSideStrips);
      gHistCounter += 1;
      
      for(Int_t iStrip = 0; iStrip < 2*fgkNumberOfPSideStrips; iStrip++)
	fOccupancyMatrix[iModule][iStrip] = 0; 
    }//module loop
  }//online flag
  else {
    for(Int_t iModule = 0; iModule < fgkSSDMODULES; iModule++) 
      fHistSSDRawSignalModule[iModule]=NULL;
    fCDBManager = NULL;
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
fSSDhTask(qadm.fSSDhTask),
fGenOffset(qadm.fGenOffset),
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
  for(Int_t iModule = 0; iModule < fgkSSDMODULES; iModule++)
    if(fHistSSDRawSignalModule[iModule]) delete fHistSSDRawSignalModule[iModule];
  if(fCDBManager) delete fCDBManager;
}

//____________________________________________________________________________ 
void AliITSQASSDDataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  AliDebug(1,"AliITSQADM::Start of SSD Cycle\n");
}

//____________________________________________________________________________ 
void AliITSQASSDDataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray* list)
{
  // launch the QA checking
  AliDebug(1,"AliITSDM instantiates checker with Run(AliQA::kITS, task, list)\n"); 
  
  //online part
  if(fkOnline) {
    //Output of the DA
    MonitorOCDBObjects();

    Int_t gHistPositionOccupancyPerModule = 0;
    Int_t gLayer = 0, gLadder = 0, gModule = 0;
    //occupancy per module
    for(Int_t iModule = 0; iModule < fgkSSDMODULES; iModule++) {
      AliITSgeomTGeo::GetModuleId(iModule+500,gLayer,gLadder,gModule);

      gHistPositionOccupancyPerModule = (gLayer == 5) ? ((gLadder - 1)*fgkSSDMODULESPERLADDERLAYER5 + gModule - 1) : ((gLadder - 1)*fgkSSDMODULESPERLADDERLAYER6 + gModule + fgkSSDMODULESLAYER5 - 1);
      for(Int_t iBins = 1; iBins < fHistSSDRawSignalModule[iModule]->GetXaxis()->GetNbins(); iBins++)
	fAliITSQADataMakerRec->GetRawsData(fGenOffset+fSSDRawsCommonLevelOffset+gHistPositionOccupancyPerModule)->SetBinContent(iBins,fOccupancyMatrix[iModule][iBins-1]);
      
      if(fSSDEventPerCycle != 0)
	((TH1D *)(fAliITSQADataMakerRec->GetRawsData(fGenOffset+fSSDRawsCommonLevelOffset+gHistPositionOccupancyPerModule)))->Scale(100./fSSDEventPerCycle);
    }//module loop

    //occupancy per ladder
    Int_t gHistPositionOccupancyPerLadder = 0;
    Int_t lLadderLocationY = 0;
    Double_t occupancy = 0.0;
    for(Int_t iModule = 0; iModule < fgkSSDMODULES; iModule++) {
      AliITSgeomTGeo::GetModuleId(iModule+500,gLayer,gLadder,gModule);
      
      gHistPositionOccupancyPerModule = (gLayer == 5) ? ((gLadder - 1)*fgkSSDMODULESPERLADDERLAYER5 + gModule - 1) : ((gLadder - 1)*fgkSSDMODULESPERLADDERLAYER6 + gModule + fgkSSDMODULESLAYER5 - 1);
      gHistPositionOccupancyPerLadder = (gLayer == 5) ? 2*(gLadder - 1) : 2*(gLadder - 1 + fgkSSDLADDERSLAYER5);
      
      //P-SIDE OCCUPANCY                                                                     
      occupancy = GetOccupancyModule((TH1 *)fAliITSQADataMakerRec->GetRawsData(fGenOffset+fSSDRawsCommonLevelOffset+fgkSSDMODULES+gHistPositionOccupancyPerModule),0);
      fAliITSQADataMakerRec->GetRawsData(fGenOffset+fSSDRawsCommonLevelOffset+fgkSSDMODULES+fgkSSDMODULES+gHistPositionOccupancyPerLadder)->Fill(gModule,occupancy);
      lLadderLocationY = 3*gLadder; // sideP=1 sideN=0 
      if(gLayer == 5)                                                               
        ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenOffset+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6))->SetBinContent(gModule,lLadderLocationY,occupancy);
      else if(gLayer == 6)
        ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenOffset+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+1))->SetBinContent(gModule,lLadderLocationY,occupancy);

      //N-SIDE OCCUPANCY                                                                           
      occupancy = GetOccupancyModule((TH1 *)fAliITSQADataMakerRec->GetRawsData(fGenOffset+fSSDRawsCommonLevelOffset+fgkSSDMODULES+gHistPositionOccupancyPerModule),1);   
      fAliITSQADataMakerRec->GetRawsData(fGenOffset+fSSDRawsCommonLevelOffset+fgkSSDMODULES+fgkSSDMODULES+gHistPositionOccupancyPerLadder+1)->Fill(gModule,occupancy);
      if(gLayer == 5)                                                               
        ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenOffset+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6))->SetBinContent(gModule,lLadderLocationY-1,occupancy);
      else if(gLayer == 6)
        ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenOffset+fSSDRawsCommonLevelOffset+fgkSSDMODULES+2*fgkSSDLADDERSLAYER5+2*fgkSSDLADDERSLAYER6+1))->SetBinContent(gModule,lLadderLocationY-1,occupancy);
    }//module loop
  }//online flag for SSD

  fSSDEventPerCycle = 0;
  
  AliQAChecker::Instance()->Run( AliQA::kITS , task, list);
}

//____________________________________________________________________________ 
void AliITSQASSDDataMakerRec::InitRaws() {  
  // Initialization for RAW data - SSD -
    fGenOffset = (fAliITSQADataMakerRec->fRawsQAList)->GetEntries();

  if(fkOnline) {
    AliInfo("Book Online Histograms for SSD\n");
  }
  else {
    AliInfo("Book Offline Histograms for SSD\n ");
  }
  AliInfo(Form("Number of histograms (SPD+SDD): %d\n",fGenOffset));
  TString gTitle = 0;
  //book online-offline QA histos
  TH1F *fHistSSDEventType = new TH1F("fHistSSDEventType",
				     ";Event type;Events",
				     31,-1,30);
  fAliITSQADataMakerRec->Add2RawsList(fHistSSDEventType, 
				      fGenOffset+fSSDRawsOffset);
  fSSDRawsOffset += 1;
  TH1F *fHistSSDDataSize = new TH1F("fHistSSDDataSize",
				    ";log(SSD data size) [Bytes];Events",
				    100,3,8);
  fAliITSQADataMakerRec->Add2RawsList(fHistSSDDataSize, 
				      fGenOffset+fSSDRawsOffset);
  fSSDRawsOffset += 1;
  TH1F *fHistSSDDataSizePercentage = new TH1F("fHistSSDDataSizePercentage",
					      ";SSD data size [%];Events",
					      100,0,100);
  fAliITSQADataMakerRec->Add2RawsList(fHistSSDDataSizePercentage, 
				      fGenOffset+fSSDRawsOffset);
  fSSDRawsOffset += 1;
  TH1F *fHistSSDDDLId = new TH1F("fHistSSDDDLId",
				 ";DDL id;Events",20,510.5,530.5);
  fAliITSQADataMakerRec->Add2RawsList(fHistSSDDDLId, 
				      fGenOffset+fSSDRawsOffset);
  fSSDRawsOffset += 1;
  TH1F *fHistSSDDataSizePerDDL = new TH1F("fHistSSDDataSizePerDDL",
					  ";DDL id;<SSD data size> [MB]",
					  20,510.5,530.5);
  fAliITSQADataMakerRec->Add2RawsList(fHistSSDDataSizePerDDL, 
				      fGenOffset+fSSDRawsOffset);
  fSSDRawsOffset += 1;
  TH1F *fHistSSDDataSizeDDL[fgkNumOfDDLs];
  for(Int_t i = 1; i < fgkNumOfDDLs+1; i++) {
    gTitle = "fHistSSDDataSizeDDL"; gTitle += i+511;
    fHistSSDDataSizeDDL[i-1] = new TH1F(gTitle.Data(),
					";log(SSD data size) [Bytes];Events",
					100,1,8);
    fAliITSQADataMakerRec->Add2RawsList(fHistSSDDataSizeDDL[i-1], 
					fGenOffset+fSSDRawsOffset);
    fSSDRawsOffset += 1;
  }
  
  TH1F *fHistSSDLDCId = new TH1F("fHistSSDLDCId",";LDC id;Events",10,0.5,10.5);
  fAliITSQADataMakerRec->Add2RawsList(fHistSSDLDCId, 
				      fGenOffset+fSSDRawsOffset);
  fSSDRawsOffset += 1;
  TH1F *fHistSSDDataSizePerLDC = new TH1F("fHistSSDDataSizePerLDC",
					  ";LDC id;<SSD data size> [MB]",
					  100,0,20);
  fAliITSQADataMakerRec->Add2RawsList(fHistSSDDataSizePerLDC, 
				      fGenOffset+fSSDRawsOffset);
  fSSDRawsOffset += 1;
  TH1F *fHistSSDDataSizeLDC[fgkNumOfLDCs];
  for(Int_t i = 1; i < fgkNumOfLDCs+1; i++) {
    gTitle = "fHistSSDDataSizeLDC"; gTitle += i;
    fHistSSDDataSizeLDC[i-1] = new TH1F(gTitle.Data(),
					";log(SSD data size) [Bytes];Events",
					100,1,8);
    fAliITSQADataMakerRec->Add2RawsList(fHistSSDDataSizeLDC[i-1], 
					fGenOffset+fSSDRawsOffset);
    fSSDRawsOffset += 1;
  }
  fSSDRawsCommonLevelOffset = fSSDRawsOffset;

  if(fkOnline) {
    Int_t gLayer = 0, gLadder = 0, gModule = 0;
    //occupancy per SSD module
    TH1D *fHistSSDOccupancyModule[fgkSSDMODULES]; 
    for(Int_t i = 500; i < fgkSSDMODULES + 500; i++) {
      AliITSgeomTGeo::GetModuleId(i,gLayer,gLadder,gModule);
      gTitle = "fHistSSD_Occupancy_Layer"; gTitle += gLayer;
      gTitle += "_Ladder"; 
      if(gLayer == 5)
	gTitle += 499+gLadder;
      if(gLayer == 6)
	gTitle += 599+gLadder;
      gTitle += "_Module"; gTitle += gModule; 
      fHistSSDOccupancyModule[i-500] = new TH1D(gTitle.Data(),gTitle.Data(),
						2*fgkNumberOfPSideStrips,0,2*fgkNumberOfPSideStrips);
      fHistSSDOccupancyModule[i-500]->GetXaxis()->SetTitleColor(1);
      fHistSSDOccupancyModule[i-500]->GetXaxis()->SetTitle("N_{strip}");
      fHistSSDOccupancyModule[i-500]->GetYaxis()->SetTitle("Occupancy [%]");
      fAliITSQADataMakerRec->Add2RawsList(fHistSSDOccupancyModule[i-500], fGenOffset+fSSDRawsOffset);
      fSSDRawsOffset += 1;
    }

    //Occupancy per SSD ladder
    Int_t occupancyCounter = 0;
    TH1D *fHistSSDOccupancyLadder[2*(fgkSSDLADDERSLAYER5 + fgkSSDLADDERSLAYER6)];
    for(Int_t iLayer = 5; iLayer < 7; iLayer++) {
      for(Int_t iLadder = 1; iLadder < AliITSgeomTGeo::GetNLadders(iLayer) + 1; iLadder++) {
	//P-side occupancy plots
	gTitle = "fHistSSD_Occupancy_Layer"; gTitle += iLayer;
	gTitle += "_Ladder"; 
	if(iLayer == 5)
	  gTitle += 499+iLadder;
	if(iLayer == 6)
	  gTitle += 599+iLadder;
	gTitle += "_PSide";
	fHistSSDOccupancyLadder[occupancyCounter] = new TH1D(gTitle.Data(),
							     gTitle.Data(),
							     AliITSgeomTGeo::GetNDetectors(iLayer),
							     0.5,AliITSgeomTGeo::GetNDetectors(iLayer)+0.5);
	fHistSSDOccupancyLadder[occupancyCounter]->GetXaxis()->SetTitleColor(1);
	fHistSSDOccupancyLadder[occupancyCounter]->GetXaxis()->SetTitle("Module number");
	fHistSSDOccupancyLadder[occupancyCounter]->GetYaxis()->SetTitle("Occupancy [%]");
	fAliITSQADataMakerRec->Add2RawsList(fHistSSDOccupancyLadder[occupancyCounter], 
					    fGenOffset+fSSDRawsOffset);
	occupancyCounter += 1; fSSDRawsOffset += 1;
	//N-side occupancy plots
	gTitle = "fHistSSD_Occupancy_Layer"; gTitle += iLayer;
	gTitle += "_Ladder"; 
	if(iLayer == 5)
	  gTitle += 499+iLadder;
	if(iLayer == 6)
	  gTitle += 599+iLadder;
	gTitle += "_NSide";
	fHistSSDOccupancyLadder[occupancyCounter] = new TH1D(gTitle.Data(),
							     gTitle.Data(),
							     AliITSgeomTGeo::GetNDetectors(iLayer),
							     0.5,AliITSgeomTGeo::GetNDetectors(iLayer)+0.5);
	fHistSSDOccupancyLadder[occupancyCounter]->GetXaxis()->SetTitleColor(1);
	fHistSSDOccupancyLadder[occupancyCounter]->GetXaxis()->SetTitle("Module number");
	fHistSSDOccupancyLadder[occupancyCounter]->GetYaxis()->SetTitle("Occupancy [%]");
	fAliITSQADataMakerRec->Add2RawsList(fHistSSDOccupancyLadder[occupancyCounter], 
					    fGenOffset+fSSDRawsOffset);
	occupancyCounter += 1; fSSDRawsOffset += 1;
      }//ladder loop
    }//layer loop

    //top level occupancy plots
    TH2D *fHistSSDOccupancyLayer5 = new TH2D("fHistSSDOccupancyLayer5",
					     ";N_{modules};N_{Ladders}",
					     fgkSSDMODULESPERLADDERLAYER5,
					     0,fgkSSDMODULESPERLADDERLAYER5,
					     3*fgkSSDLADDERSLAYER5,
					     0,fgkSSDLADDERSLAYER5);  
    Char_t fLabel[3];
    for(Int_t iBin = 1; iBin < fgkSSDMODULESPERLADDERLAYER5 + 1; iBin++){
      sprintf(fLabel,"%d",iBin);
      fHistSSDOccupancyLayer5->GetXaxis()->SetBinLabel(iBin,fLabel);
    }
    fAliITSQADataMakerRec->Add2RawsList(fHistSSDOccupancyLayer5, fGenOffset+fSSDRawsOffset);
    fSSDRawsOffset += 1;
    TH2D *fHistSSDOccupancyLayer6 = new TH2D("fHistSSDOccupancyLayer6",
					     ";N_{modules};N_{Ladders}",
					     fgkSSDMODULESPERLADDERLAYER6,
					     0,fgkSSDMODULESPERLADDERLAYER6,
					     3*fgkSSDLADDERSLAYER6,
					     0,fgkSSDLADDERSLAYER6);
    for(Int_t iBin = 1; iBin < fgkSSDMODULESPERLADDERLAYER6 + 1; iBin++){
      sprintf(fLabel,"%d",iBin);
      fHistSSDOccupancyLayer6->GetXaxis()->SetBinLabel(iBin,fLabel);
    }
    fAliITSQADataMakerRec->Add2RawsList(fHistSSDOccupancyLayer6, fGenOffset+fSSDRawsOffset);
    fSSDRawsOffset += 1;

    //Output of the DA
    TH2D *fHistPSideBadChannelMapLayer5 = new TH2D("fHistPSideBadChannelMapLayer5",
						   "Layer 5;N_{module};N_{ladder}",
						   22,1,23,
						   34,500,534);
    fHistPSideBadChannelMapLayer5->GetXaxis()->SetTitleColor(1);
    fHistPSideBadChannelMapLayer5->SetStats(kFALSE);
    fHistPSideBadChannelMapLayer5->GetYaxis()->SetTitleOffset(1.8);
    fHistPSideBadChannelMapLayer5->GetXaxis()->SetNdivisions(22);
    fHistPSideBadChannelMapLayer5->GetYaxis()->SetNdivisions(34);
    fHistPSideBadChannelMapLayer5->GetXaxis()->SetLabelSize(0.03);
    fHistPSideBadChannelMapLayer5->GetYaxis()->SetLabelSize(0.03);
    fHistPSideBadChannelMapLayer5->GetZaxis()->SetTitleOffset(1.6);
    fHistPSideBadChannelMapLayer5->GetZaxis()->SetTitle("Bad channels (p-side)[%]");
    fAliITSQADataMakerRec->Add2RawsList(fHistPSideBadChannelMapLayer5, fGenOffset+fSSDRawsOffset);
    fSSDRawsOffset += 1; fSSDRawsDAOffset += 1;

    TH2D *fHistNSideBadChannelMapLayer5 = new TH2D("fHistNSideBadChannelMapLayer5",
						   "Layer 5;N_{module};N_{ladder}",
						   22,1,23,
						   34,500,534);
    fHistNSideBadChannelMapLayer5->GetXaxis()->SetTitleColor(1);
    fHistNSideBadChannelMapLayer5->SetStats(kFALSE);
    fHistNSideBadChannelMapLayer5->GetYaxis()->SetTitleOffset(1.8);
    fHistNSideBadChannelMapLayer5->GetXaxis()->SetNdivisions(22);
    fHistNSideBadChannelMapLayer5->GetYaxis()->SetNdivisions(34);
    fHistNSideBadChannelMapLayer5->GetXaxis()->SetLabelSize(0.03);
    fHistNSideBadChannelMapLayer5->GetYaxis()->SetLabelSize(0.03);
    fHistNSideBadChannelMapLayer5->GetZaxis()->SetTitleOffset(1.6);
    fHistNSideBadChannelMapLayer5->GetZaxis()->SetTitle("Bad channels (n-side)[%]");
    fAliITSQADataMakerRec->Add2RawsList(fHistNSideBadChannelMapLayer5, fGenOffset+fSSDRawsOffset);
    fSSDRawsOffset += 1; fSSDRawsDAOffset += 1;

    TH2D *fHistPSideBadChannelMapLayer6 = new TH2D("fHistPSideBadChannelMapLayer6",
						   "Layer 6;N_{module};N_{ladder}",
						   25,1,26,
						   38,600,638);
    fHistPSideBadChannelMapLayer6->GetXaxis()->SetTitleColor(1);
    fHistPSideBadChannelMapLayer6->SetStats(kFALSE);
    fHistPSideBadChannelMapLayer6->GetYaxis()->SetTitleOffset(1.8);
    fHistPSideBadChannelMapLayer6->GetXaxis()->SetNdivisions(25);
    fHistPSideBadChannelMapLayer6->GetYaxis()->SetNdivisions(38);
    fHistPSideBadChannelMapLayer6->GetXaxis()->SetLabelSize(0.03);
    fHistPSideBadChannelMapLayer6->GetYaxis()->SetLabelSize(0.03);
    fHistPSideBadChannelMapLayer6->GetZaxis()->SetTitleOffset(1.6);
    fHistPSideBadChannelMapLayer6->GetZaxis()->SetTitle("Bad channels (p-side)[%]");
    fAliITSQADataMakerRec->Add2RawsList(fHistPSideBadChannelMapLayer6, fGenOffset+fSSDRawsOffset);
    fSSDRawsOffset += 1; fSSDRawsDAOffset += 1;

    TH2D *fHistNSideBadChannelMapLayer6 = new TH2D("fHistNSideBadChannelMapLayer6",
						   "Layer 6;N_{module};N_{ladder}",
						   25,1,26,
						   38,600,638);
    fHistNSideBadChannelMapLayer6->GetXaxis()->SetTitleColor(1);
    fHistNSideBadChannelMapLayer6->SetStats(kFALSE);
    fHistNSideBadChannelMapLayer6->GetYaxis()->SetTitleOffset(1.8);
    fHistNSideBadChannelMapLayer6->GetXaxis()->SetNdivisions(25);
    fHistNSideBadChannelMapLayer6->GetYaxis()->SetNdivisions(38);
    fHistNSideBadChannelMapLayer6->GetXaxis()->SetLabelSize(0.03);
    fHistNSideBadChannelMapLayer6->GetYaxis()->SetLabelSize(0.03);
    fHistNSideBadChannelMapLayer6->GetZaxis()->SetTitleOffset(1.6);
    fHistNSideBadChannelMapLayer6->GetZaxis()->SetTitle("Bad channels (n-side)[%]");
    fAliITSQADataMakerRec->Add2RawsList(fHistNSideBadChannelMapLayer6, fGenOffset+fSSDRawsOffset);
    fSSDRawsOffset += 1; fSSDRawsDAOffset += 1;
  }//online flag
  fSSDhTask = fSSDRawsOffset;
  AliDebug(1,Form("%d SSD Raws histograms booked\n",fSSDhTask));
  AliInfo(Form("Number of histograms (SPD+SDD+SSD): %d\n",fGenOffset+fSSDhTask));  
  AliDebug(1,Form("Number of histograms (SPD+SDD+SSD): %d\n",fGenOffset+fSSDRawsOffset));
}

//____________________________________________________________________________
void AliITSQASSDDataMakerRec::MakeRaws(AliRawReader* rawReader) { 
  // Fill QA for RAW - SSD -
  Int_t gStripNumber;
  Int_t gHistPosition;
  Int_t gLayer = 0,gLadder = 0, gModule = 0;
  
  Double_t gSizePerDDL[fgkNumOfDDLs] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t gSizePerLDC[fgkNumOfLDCs] = {0.,0.,0.};
  Double_t sumSSDDataSize = 0.0;
  Double_t eventSize = -1.0;

  if(fkOnline) {
    //reset the signal vs strip number histograms
    for(Int_t iModule = 0; iModule < fgkSSDMODULES; iModule++)
      fHistSSDRawSignalModule[iModule]->Reset();
  }//online flag

  rawReader->Select("ITSSSD",-1,-1);  
  rawReader->Reset(); //rawReader->NextEvent();   
  (fAliITSQADataMakerRec->GetRawsData(fGenOffset))->Fill(rawReader->GetType());
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
    gSizePerDDL[rawReader->GetDDLID()] = rawReader->GetDataSize();
    gSizePerLDC[rawReader->GetLDCId()-6] = rawReader->GetDataSize();
    AliITSgeomTGeo::GetModuleId(gSSDStream.GetModuleID(),gLayer,gLadder,gModule);
    gStripNumber = (gSSDStream.GetSideFlag() == 0) ? gSSDStream.GetStrip() : gSSDStream.GetStrip() + fgkNumberOfPSideStrips;
    gHistPosition = (gLayer == 5) ? ((gLadder - 1)*fgkSSDMODULESPERLADDERLAYER5 + gModule - 1) : ((gLadder - 1)*fgkSSDMODULESPERLADDERLAYER6 + gModule + fgkSSDMODULESLAYER5 - 1);
    //AliInfo(Form("ModulePosition: %d - Layer: %d - Ladder: %d - Module: %d\n",gHistPosition,gLayer,gLadder,gModule));
    if(fkOnline)
      fHistSSDRawSignalModule[gHistPosition]->Fill(gStripNumber,gSSDStream.GetSignal());
    //fAliITSQADataMakerRec->GetRawsData(fGenOffset+gHistPosition+fSSDRawsCommonLevelOffset)->Fill(gStripNumber,gSSDStream.GetSignal());
  }//streamer loop   
  
  //event size calculation and filling info
  for(Int_t i = 0; i < fgkNumOfDDLs; i++) {
    sumSSDDataSize += gSizePerDDL[i];
    if(gSizePerDDL[i] > 0) {
      (fAliITSQADataMakerRec->GetRawsData(fGenOffset+3))->Fill(i+512);
      (fAliITSQADataMakerRec->GetRawsData(fGenOffset+5+i))->Fill(TMath::Log10(gSizePerDDL[i]));
    }
    (fAliITSQADataMakerRec->GetRawsData(fGenOffset+4))->Fill(i+512,gSizePerDDL[i]/1e+06);
  }
  for(Int_t i = 0; i < fgkNumOfLDCs; i++) {
    if(gSizePerLDC[i] > 0) {
      (fAliITSQADataMakerRec->GetRawsData(fGenOffset+21))->Fill(i+6);
      (fAliITSQADataMakerRec->GetRawsData(fGenOffset+23+i))->Fill(TMath::Log10(gSizePerLDC[i]));
    }
    (fAliITSQADataMakerRec->GetRawsData(fGenOffset+22))->Fill(i+6,gSizePerLDC[i]/1e+06);
  }
  if(sumSSDDataSize) 
    (fAliITSQADataMakerRec->GetRawsData(fGenOffset+1))->Fill(TMath::Log10(sumSSDDataSize));
  if(eventSize)
    (fAliITSQADataMakerRec->GetRawsData(fGenOffset+2))->Fill(100.*sumSSDDataSize/eventSize);

  //Occupancy calculation
  if(fkOnline) {
    for(Int_t iModule = 0; iModule < fgkSSDMODULES; iModule++)
      GetOccupancyStrip(fHistSSDRawSignalModule[iModule],fOccupancyMatrix[iModule]);
  }//online flag for SSD
}

//____________________________________________________________________________ //
void AliITSQASSDDataMakerRec::GetOccupancyStrip(TH1 *lHisto, Int_t *occupancyMatrix) { 
  //Increments the entries in the occupancy matrix based 
  //on whether the signal for each strip is larger than the cutoff
  //Currently the cutoff is at 0 which means that if ZS
  //works, whatever comes from the FEROM is considered as "signal"
  Double_t cutoff = 0.0;
  for(Int_t iBin = 1; iBin < lHisto->GetXaxis()->GetNbins(); iBin++) {
    Double_t y = lHisto->GetBinContent(iBin);
    if(y > cutoff) {
      occupancyMatrix[iBin-1] += 1;
      //cout<<"Event: "<<iEvent<<" - Strip: "<<i<<" - Signal: "<<fRawSignal->GetBinContent(i)<<                              
      //" - Occupancy: "<<fOccupancy[i-1]<<endl;                                                                  
    }
  }
}

//____________________________________________________________________________ 
Double_t AliITSQASSDDataMakerRec::GetOccupancyModule(TH1 *lHisto, Int_t stripside) { 
  // bo: TDC >0 or # of sigmas wrt noise ?
  //stripside == 0 --> P-side
  //stripside == 1 --> N-side
  Int_t lNumFiredBins = 0;
  for(Int_t iBin = 1 + stripside*fgkNumberOfPSideStrips; iBin < fgkNumberOfPSideStrips*(1 + stripside); iBin++){
    if (lHisto->GetBinContent(iBin) > 0)
      lNumFiredBins++; 
  }
  
  Double_t lOccupancy = (100.*lNumFiredBins)/fgkNumberOfPSideStrips; // percentage
  //AliInfo(Form("Fired strips: %d - Total strips: %d - Occupancy :%lf\n",lNumFiredBins,lHisto->GetNbinsX(),lOccupancy));
  
  return lOccupancy;
}

//____________________________________________________________________________
void AliITSQASSDDataMakerRec::MonitorOCDBObjects() { 
  //Monitor in AMORE the output of the DA
  //Currently only the bad channel list is monitored
  //Todo: Noise - Pedestal
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
    if(layer == 5) {
      if(nPSideChannelsLayer5 > 0)
	((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenOffset+fSSDRawsOffset-fSSDRawsDAOffset))->Fill(module,499+ladder,
									    100.*nPSideChannelsLayer5/fgkNumberOfPSideStrips);
      else ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenOffset+fSSDRawsOffset-fSSDRawsDAOffset))->Fill(module,499+ladder,0.0001);
      if(nNSideChannelsLayer5 > 0)
	((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenOffset+fSSDRawsOffset-fSSDRawsDAOffset+1))->Fill(module,499+ladder,
									    100.*nNSideChannelsLayer5/fgkNumberOfPSideStrips);
      else ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenOffset+fSSDRawsOffset-fSSDRawsDAOffset+1))->Fill(module,499+ladder,0.0001);
    }//layer 5                                                                                                                      
    if(layer == 6) {
      if(nPSideChannelsLayer6 > 0)
        ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenOffset+fSSDRawsOffset-fSSDRawsDAOffset+2))->Fill(module,599+ladder,
									    100.*nPSideChannelsLayer6/fgkNumberOfPSideStrips);
      else ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenOffset+fSSDRawsOffset-fSSDRawsDAOffset+2))->Fill(module,599+ladder,0.0001);
      if(nNSideChannelsLayer6 > 0)
        ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenOffset+fSSDRawsOffset-fSSDRawsDAOffset+3))->Fill(module,599+ladder,
									    100.*nNSideChannelsLayer6/fgkNumberOfPSideStrips);
      else ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenOffset+fSSDRawsOffset-fSSDRawsDAOffset+3))->Fill(module,599+ladder,0.0001);
    }//layer 6                                                                                                                      
  }//module loop
}

//____________________________________________________________________________ 
void AliITSQASSDDataMakerRec::InitRecPoints()
{
  // Initialization for RECPOINTS - SSD -
 
  fGenOffset = (fAliITSQADataMakerRec->fRecPointsQAList)->GetEntries();
  Int_t nModuleOffset = 500;
  Int_t nITSTotalModules = AliITSgeomTGeo::GetNModules();

  TH1F *fHistModuleIdLayer5 = new TH1F("fHistModuleIdLayer5",
				       "Module Id - Layer 5;Module Id;Entries",
				       fgkSSDMODULESLAYER5,
				       nModuleOffset - 0.5,
				       nITSTotalModules-fgkSSDMODULESLAYER6+0.5);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistModuleIdLayer5, 
					   fGenOffset + 0);
  fSSDhTask += 1;
  TH1F *fHistModuleIdLayer6 = new TH1F("fHistModuleIdLayer6",
				       "Module Id - Layer 6;Module Id;Entries",
				       fgkSSDMODULESLAYER6,
				       nModuleOffset+fgkSSDMODULESLAYER5 - 0.5,
				       nITSTotalModules + 0.5);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistModuleIdLayer6, 
					   fGenOffset + 1);
  fSSDhTask += 1;
  TH1F *fHistClusterPerEventLayer5 = new TH1F("fHistClusterPerEventLayer5",
					      "N_{clusters} - Layer 5;N_{clusters};Entries;",
					      100,0.1,5000);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistClusterPerEventLayer5,
					   fGenOffset + 2);
  fSSDhTask += 1;
  TH1F *fHistClusterPerEventLayer6 = new TH1F("fHistClusterPerEventLayer6",
					      "N_{clusters} - Layer 6;N_{clusters};Entries;",
					      100,0.1,5000);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistClusterPerEventLayer6,
					   fGenOffset + 3);
  fSSDhTask += 1;
  TH1F *fHistLocalXLayer5 = new TH1F("fHistLocalXLayer5",
				     "Local x coord.- Layer 5;x [cm];Entries;",
				     100,-4.,4.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistLocalXLayer5,
					   fGenOffset + 4);
  fSSDhTask += 1;
  TH1F *fHistLocalXLayer6 = new TH1F("fHistLocalXLayer6",
				     "Local x coord.- Layer 6;x [cm];Entries;",
				     100,-4.,4.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistLocalXLayer6, 
					   fGenOffset + 5);
  fSSDhTask += 1;
  TH1F *fHistLocalZLayer5 = new TH1F("fHistLocalZLayer5",
				     "Local z coord.- Layer 5;z [cm];Entries;",
				     100,-4.,4.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistLocalZLayer5, 
					   fGenOffset + 6);
  fSSDhTask += 1;
  TH1F *fHistLocalZLayer6 = new TH1F("fHistLocalZLayer6",
				     "Local z coord.- Layer 6;z [cm];Entries;",
				     100,-4.,4.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistLocalZLayer6, 
					   fGenOffset + 7);
  fSSDhTask += 1;
  TH1F *fHistGlobalXLayer5 = new TH1F("fHistGlobalXLayer5",
				      "Global x - Layer 5;x [cm];Entries;",
				      100,-40.,40.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistGlobalXLayer5, 
					   fGenOffset + 8);
  fSSDhTask += 1;
  TH1F *fHistGlobalXLayer6 = new TH1F("fHistGlobalXLayer6",
				      "Global x - Layer 6;x [cm];Entries;",
				      100,-45.,45.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistGlobalXLayer6, 
					   fGenOffset + 9);
  fSSDhTask += 1;
  TH1F *fHistGlobalYLayer5 = new TH1F("fHistGlobalYLayer5",
				      "Global y - Layer 5;y [cm];Entries;",
				      100,-40.,40);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistGlobalYLayer5, 
					   fGenOffset + 10);
  fSSDhTask += 1;
  TH1F *fHistGlobalYLayer6 = new TH1F("fHistGlobalYLayer6",
				      "Global y - Layer 6;y [cm];Entries;",
				      100,-45.,45.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistGlobalYLayer6, 
					   fGenOffset + 11);
  fSSDhTask += 1;
  TH1F *fHistGlobalZLayer5 = new TH1F("fHistGlobalZLayer5",
				      "Global z - Layer 5;z [cm];Entries;",
				      100,-45.,45);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistGlobalZLayer5, 
					   fGenOffset + 12);
  fSSDhTask += 1;
  TH1F *fHistGlobalZLayer6 = new TH1F("fHistGlobalZLayer6",
				      "Global z - Layer 6;z [cm];Entries;",
				      100,-55.,55.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistGlobalZLayer6, 
					   fGenOffset + 13);
  fSSDhTask += 1;
  TH1F *fHistPhiLayer5 = new TH1F("fHistPhiLayer5",
				  "#phi - Layer 5;#phi [rad];Entries;",
				  100,-TMath::Pi(),TMath::Pi());
  fAliITSQADataMakerRec->Add2RecPointsList(fHistPhiLayer5, 
					   fGenOffset + 14);
  fSSDhTask += 1;
  TH1F *fHistPhiLayer6 = new TH1F("fHistPhiLayer6",
				  "#phi - Layer 6;#phi [rad];Entries;",
				  100,-TMath::Pi(),TMath::Pi());
  fAliITSQADataMakerRec->Add2RecPointsList(fHistPhiLayer6, 
					   fGenOffset + 15);
  fSSDhTask += 1;
  TH1F *fHistThetaLayer5 = new TH1F("fHistThetaLayer5",
				    "#theta - Layer 5;#theta [rad];Entries;",
				    100,-TMath::Pi(),TMath::Pi());
  fAliITSQADataMakerRec->Add2RecPointsList(fHistThetaLayer5, 
					   fGenOffset + 16);
  fSSDhTask += 1;
  TH1F *fHistThetaLayer6 = new TH1F("fHistThetaLayer6",
				    "#theta - Layer 6;#theta [rad];Entries;",
				    100,-TMath::Pi(),TMath::Pi());
  fAliITSQADataMakerRec->Add2RecPointsList(fHistThetaLayer6, 
					   fGenOffset + 17);
  fSSDhTask += 1;
  TH1F *fHistRadiusLayer5 = new TH1F("fHistRadiusLayer5",
				     "r - Layer 5;r [cm];Entries;",
				     100,35.,50.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistRadiusLayer5, 
					   fGenOffset + 18);
  fSSDhTask += 1;
  TH1F *fHistRadiusLayer6 = new TH1F("fHistRadiusLayer6",
				     "r - Layer 6;r [cm];Entries;",
				     100,35.,50.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistRadiusLayer6, 
					   fGenOffset + 19);
  fSSDhTask += 1;
  TH1F *fHistClusterTypeLayer5 = new TH1F("fHistClusterTypeLayer5",
					  "CL type - Layer 5;Cluster type;Entries;",
					  150,0,150);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistClusterTypeLayer5, 
					   fGenOffset + 20);
  fSSDhTask += 1;
  TH1F *fHistClusterTypeLayer6 = new TH1F("fHistClusterTypeLayer6",
					  "CL type - Layer 6;Cluster type;Entries;",
					  150,0,150);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistClusterTypeLayer6, 
					   fGenOffset + 21);
  fSSDhTask += 1;
  TH1F *fHistChargeRatioLayer5 = new TH1F("fHistChargeRatioLayer5",
					  "Charge ratio - Layer 5;q_{ratio};Entries;",
					  100,-2.0,2.0);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargeRatioLayer5, 
					   fGenOffset + 22);
  fSSDhTask += 1;
  TH1F *fHistChargeRatioLayer6 = new TH1F("fHistChargeRatioLayer6",
					  "Charge ratio - Layer 6;q_{ratio};Entries;",
					  100,-2.0,2.0);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargeRatioLayer6, 
					   fGenOffset + 23);
  fSSDhTask += 1;
  TH1F *fHistChargekeVLayer5 = new TH1F("fHistChargekeVLayer5",
					"Charge - Layer 5;q [keV];Entries;",
					100,0.,300.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargekeVLayer5, 
					   fGenOffset + 24);
  fSSDhTask += 1;
  TH1F *fHistChargekeVLayer6 = new TH1F("fHistChargekeVLayer6",
					"Charge - Layer 6;q [keV];Entries;",
					100,0.,300.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargekeVLayer6, 
					   fGenOffset + 25);
  fSSDhTask += 1;
  TH1F *fHistChargePSideLayer5 = new TH1F("fHistChargePSideLayer5",
                                          "Charge P- Layer 5;q_{P} [keV];Entries;",
                                          100,0.,300.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargePSideLayer5,
					   fGenOffset + 26);
  fSSDhTask += 1;
  TH1F *fHistChargePSideLayer6 = new TH1F("fHistChargePSideLayer6",
                                          "Charge P- Layer 6;q_{P} [keV];Entries;",
                                          100,0.,300.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargePSideLayer6,
					   fGenOffset + 27);
  fSSDhTask += 1;
  TH1F *fHistChargeNSideLayer5 = new TH1F("fHistChargeNSideLayer5",
                                          "Charge N- Layer 5;q_{N} [keV];Entries;",
                                          100,0.,300.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargeNSideLayer5,
					   fGenOffset + 28);
  fSSDhTask += 1;
  TH1F *fHistChargeNSideLayer6 = new TH1F("fHistChargeNSideLayer6",
                                          "Charge N- Layer 6;q_{N} [keV];Entries;",
                                          100,0.,300.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargeNSideLayer6,
					   fGenOffset + 29);
  fSSDhTask += 1;
  TH1F *fHistChargeRatio2Layer5 = new TH1F("fHistChargeRatio2Layer5",
                                          "Charge Ratio qN/qP - Layer 5;q_{N}/q_{P};Entries;",
                                          100,0,2);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargeRatio2Layer5,
					   fGenOffset + 30);
  fSSDhTask += 1;
  TH1F *fHistChargeRatio2Layer6 = new TH1F("fHistChargeRatio2Layer6",
					   "Charge Ratio qN/qP - Layer 6;q_{N}/q_{P};Entries;",
					   100,0,2);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargeRatio2Layer6,
					   fGenOffset + 31);
  fSSDhTask += 1;
  TH2F *fHistChargePNSideLayer5 = new TH2F("fHistChargePNSideLayer5",
                                           "Charge correlation - Layer 5;q_{P} [keV];q_{N} [keV]",
                                           100,0.,300.,
                                           100,0.,300.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargePNSideLayer5,
					   fGenOffset + 32);
  fSSDhTask += 1;
  TH2F *fHistChargePNSideLayer6 = new TH2F("fHistChargePNSideLayer6",
                                           "Charge correlation - Layer 6;q_{P} [keV];q_{N} [keV]",
                                           100,0.,300.,
                                           100,0.,300.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargePNSideLayer6,
					   fGenOffset + 33);
  fSSDhTask += 1;
  TH2F *fHistChargeMapLayer5 = new TH2F("fHistChargeMapLayer5",
					"Charge map;N_{modules};N_{Ladders}",
					fgkSSDMODULESPERLADDERLAYER5,
					-0.5,fgkSSDMODULESPERLADDERLAYER5+0.5,
					3*fgkSSDLADDERSLAYER5,
					-0.5,fgkSSDLADDERSLAYER5+0.5);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargeMapLayer5, 
					   fGenOffset + 34);
  fSSDhTask += 1;
  TH2F *fHistChargeMapLayer6 = new TH2F("fHistChargeMapLayer6",
					"Charge map;N_{modules};N_{Ladders}",
					fgkSSDMODULESPERLADDERLAYER6,
					-0.5,fgkSSDMODULESPERLADDERLAYER6+0.5,
					3*fgkSSDLADDERSLAYER6,
					-0.5,fgkSSDLADDERSLAYER6+0.5);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargeMapLayer6, 
					   fGenOffset + 35);
  fSSDhTask += 1;

  AliDebug(1,Form("%d SSD Recs histograms booked\n",fSSDhTask));
}

//____________________________________________________________________________ 
void AliITSQASSDDataMakerRec::MakeRecPoints(TTree *clustersTree)
{
  // Fill QA for recpoints - SSD -

  Int_t gLayer = 0, gLadder = 0, gModule = 0;
  Int_t lLadderLocationY = 0;
  TBranch *branchRecP = clustersTree->GetBranch("ITSRecPoints");
  if (!branchRecP) { 
    AliError("can't get the branch with the ITS clusters !");
    return;
  }
  static TClonesArray statRecpoints("AliITSRecPoint");
  TClonesArray *recpoints = &statRecpoints;
  branchRecP->SetAddress(&recpoints);
  Int_t nClustersLayer5 = 0, nClustersLayer6 = 0;
  Int_t npoints = 0;      
  Float_t cluglo[3]={0.,0.,0.}; 
  for(Int_t module = 0; module < clustersTree->GetEntries(); module++){
    branchRecP->GetEvent(module);
    npoints += recpoints->GetEntries();
    AliITSgeomTGeo::GetModuleId(module,gLayer,gLadder,gModule);
    lLadderLocationY = 3*gLadder;

    for(Int_t j = 0;j < recpoints->GetEntries(); j++){
      AliITSRecPoint *recp = (AliITSRecPoint*)recpoints->At(j);    
      Int_t layer = recp->GetLayer();
      recp->GetGlobalXYZ(cluglo);
      Float_t radius = TMath::Sqrt(cluglo[0]*cluglo[0]+cluglo[1]*cluglo[1]); 
      Float_t phi = TMath::ATan2(cluglo[1],cluglo[0]);
      Float_t theta = TMath::ATan2(radius,cluglo[2]);
      Double_t chargeRatio = recp->GetChargeRatio();
      Double_t clusterCharge = recp->GetQ();
      Double_t chargePSide = clusterCharge*(1. + chargeRatio);
      Double_t chargeNSide = clusterCharge*(1. - chargeRatio);
      if(layer == 4) {
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 0)->Fill(module);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 4)->Fill(recp->GetDetLocalX());
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 6)->Fill(recp->GetDetLocalZ());
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 8)->Fill(cluglo[0]);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 10)->Fill(cluglo[1]);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 12)->Fill(cluglo[2]);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 14)->Fill(phi);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 16)->Fill(theta);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 18)->Fill(radius);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 20)->Fill(recp->GetType());
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 22)->Fill(recp->GetChargeRatio());
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 24)->Fill(recp->GetQ());
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 26)->Fill(chargePSide);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 28)->Fill(chargeNSide);
	if(chargePSide != 0.) fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 30)->Fill(chargeNSide/chargePSide);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 32)->Fill(chargePSide,chargeNSide);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 34)->SetBinContent(gModule,lLadderLocationY,recp->GetQ());
	nClustersLayer5 += 1;
      }//layer 5 histograms
      if(layer == 5) {
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 1)->Fill(module);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 5)->Fill(recp->GetDetLocalX());
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 7)->Fill(recp->GetDetLocalZ());
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 9)->Fill(cluglo[0]);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 11)->Fill(cluglo[1]);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 13)->Fill(cluglo[2]);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 15)->Fill(phi);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 17)->Fill(theta);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 19)->Fill(radius);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 21)->Fill(recp->GetType());
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 23)->Fill(recp->GetChargeRatio());
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 25)->Fill(recp->GetQ());
        fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 27)->Fill(chargePSide);
        fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 29)->Fill(chargeNSide);
        if(chargePSide != 0.) fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 31)->Fill(chargeNSide/chargePSide);
        fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 33)->Fill(chargePSide,chargeNSide);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 35)->SetBinContent(gModule,lLadderLocationY,recp->GetQ());
	nClustersLayer6 += 1;
      }//layer 6 histograms
    }//rec. points loop
  }//module loop

  fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 2)->Fill(nClustersLayer5);
  fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 3)->Fill(nClustersLayer6);

  statRecpoints.Clear();
}

