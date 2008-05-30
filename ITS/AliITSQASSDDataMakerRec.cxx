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

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"

ClassImp(AliITSQASSDDataMakerRec)

AliITSQASSDDataMakerRec::AliITSQASSDDataMakerRec(AliITSQADataMakerRec *aliITSQADataMakerRec, Bool_t kMode, Int_t ldc) :
TObject(),
fAliITSQADataMakerRec(aliITSQADataMakerRec),
fSSDEvent(0),
fkOnline(kMode),
fLDC(ldc),
fSSDRawsOffset(0),
fSSDhTask(0),
fGenOffset(0) {
  // Default constructor   
  //initilize the raw signal vs strip number histograms
  Int_t fLayer = 0,fLadder = 0, fModule = 0;
  TString fTitle; Int_t fHistCounter = 0;
  if(fkOnline) {
    AliCDBManager * man = AliCDBManager::Instance();
    man->SetDefaultStorage("local://$ALICE_ROOT");
    Int_t runNumber = atoi(gSystem->Getenv("DATE_RUN_NUMBER"));
    if(!runNumber) 
      AliInfo("DATE_RUN_NUMBER not defined!!!\n");

    man->SetRun(runNumber);
    AliCDBEntry *geomGRP = man->Get("GRP/Geometry/Data");
    if(!geomGRP) AliInfo("GRP geometry not found!!!\n");

    for(Int_t i = 500; i < fgkSSDMODULES + 500; i++) {
      AliITSgeomTGeo::GetModuleId(i,fLayer,fLadder,fModule);
      fTitle = "SSD_RawSignal_Layer"; fTitle += fLayer;
      fTitle += "_Ladder"; fTitle += fLadder;
      fTitle += "_Module"; fTitle += fModule; 
      fHistSSDRawSignalModule[fHistCounter] = new TH1D(fTitle.Data(),fTitle.Data(),1540,0,1540);
      fHistCounter += 1;
    }
  }//online flag
  else {
    for(Int_t i = 0; i < fgkSSDMODULES; i++) {
      fHistSSDRawSignalModule[i]=NULL;
    }
  }
}

//____________________________________________________________________________ 
AliITSQASSDDataMakerRec::AliITSQASSDDataMakerRec(const AliITSQASSDDataMakerRec& qadm) :
TObject(),
fAliITSQADataMakerRec(qadm.fAliITSQADataMakerRec),
fSSDEvent(qadm.fSSDEvent),
fkOnline(qadm.fkOnline),
fLDC(qadm.fLDC),
fSSDRawsOffset(qadm.fSSDRawsOffset),
fSSDhTask(qadm.fSSDhTask),
fGenOffset(qadm.fGenOffset) {
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
  for(Int_t i = 0; i < fgkSSDMODULES; i++) 
    if(fHistSSDRawSignalModule[i]) delete fHistSSDRawSignalModule[i];
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
  
  //scaling SSD occupancy plots
  if(fkOnline) {
    if(fSSDEvent != 0) {
      ((TH2F *)fAliITSQADataMakerRec->GetRawsData(fGenOffset+28))->Scale(1./fSSDEvent);
      ((TH2F *)fAliITSQADataMakerRec->GetRawsData(fGenOffset+29))->Scale(1./fSSDEvent);
    }
    Int_t fLayer = 0, fLadder = 0;
    Int_t fOccPosition = 0; 
    for(Int_t iLayer = 5; iLayer < 7; iLayer++) {
      for(Int_t iLadder = 1; iLadder < AliITSgeomTGeo::GetNLadders(iLayer) + 1; iLadder++) {
	fOccPosition  = (fLayer == 5) ? fLadder : fLadder + fgkSSDLADDERSLAYER5;
	if(fSSDEvent != 0) {
	  //P-SIDE OCCUPANCY - scaling
	  fAliITSQADataMakerRec->GetRawsData(fGenOffset + fSSDRawsOffset - 1 + 2*fOccPosition - 1)->Scale(1./fSSDEvent);
	  //N-SIDE OCCUPANCY - scaling
	  fAliITSQADataMakerRec->GetRawsData(fGenOffset + fSSDRawsOffset - 1 + 2*fOccPosition)->Scale(1./fSSDEvent);
	}//ssd events != 0
      }//ladder loop
    }//layer loop
  }//online flag for SSD
  
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
  TString fTitle = 0;
  //book online QA histos
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
      fTitle = "fHistSSDDataSizeDDL"; fTitle += i+511;
      fHistSSDDataSizeDDL[i-1] = new TH1F(fTitle.Data(),
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
      fTitle = "fHistSSDDataSizeLDC"; fTitle += i;
      fHistSSDDataSizeLDC[i-1] = new TH1F(fTitle.Data(),
					  ";log(SSD data size) [Bytes];Events",
					  100,1,8);
      fAliITSQADataMakerRec->Add2RawsList(fHistSSDDataSizeLDC[i-1], 
					  fGenOffset+fSSDRawsOffset);
      fSSDRawsOffset += 1;
    }
    fSSDhTask = fSSDRawsOffset;
  if(fkOnline) {
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

    //Occupancy per ladder
    fSSDhTask = fSSDRawsOffset;
    TH1D *fHistOccupancyLadder[2*(fgkSSDLADDERSLAYER5 + fgkSSDLADDERSLAYER6)];
    for(Int_t iLayer = 5; iLayer < 7; iLayer++) {
      for(Int_t iLadder = 1; iLadder < AliITSgeomTGeo::GetNLadders(iLayer) + 1; iLadder++) {
	//P-side occupancy plots
	fTitle = "fHistSSD_Occupancy_Layer"; fTitle += iLayer;
	fTitle += "_Ladder"; fTitle += iLadder; fTitle += "_PSide";
	fHistOccupancyLadder[fSSDhTask] = 
	  new TH1D(fTitle.Data(),
		   fTitle.Data(),
		   AliITSgeomTGeo::GetNDetectors(iLayer),
		   0.5,AliITSgeomTGeo::GetNDetectors(iLayer)+0.5);
	fHistOccupancyLadder[fSSDhTask]->GetXaxis()->SetTitleColor(1);
	fHistOccupancyLadder[fSSDhTask]->GetXaxis()->SetTitle("Module number");
	fHistOccupancyLadder[fSSDhTask]->GetYaxis()->SetTitle("Occupancy [%]");
	fAliITSQADataMakerRec->Add2RawsList(fHistOccupancyLadder[fSSDhTask], 
					    fGenOffset+fSSDhTask);	
	fSSDhTask++;
	//N-side occupancy plots
	fTitle = "fHistSSD_Occupancy_Layer"; fTitle += iLayer;
	fTitle += "_Ladder"; fTitle += iLadder; fTitle += "_NSide";
	fHistOccupancyLadder[fSSDhTask] = 
	  new TH1D(fTitle.Data(),
		   fTitle.Data(),
		   AliITSgeomTGeo::GetNDetectors(iLayer),
		   0.5,AliITSgeomTGeo::GetNDetectors(iLayer)+0.5);
	fHistOccupancyLadder[fSSDhTask]->GetXaxis()->SetTitleColor(1);
	fHistOccupancyLadder[fSSDhTask]->GetXaxis()->SetTitle("Module number");
	fHistOccupancyLadder[fSSDhTask]->GetYaxis()->SetTitle("Occupancy [%]");
	fAliITSQADataMakerRec->Add2RawsList(fHistOccupancyLadder[fSSDhTask], 
					    fGenOffset+fSSDhTask);	
	fSSDhTask++;
      }//ladder loop
    }//layer loop
  }//online flag
  AliDebug(1,Form("%d SSD Raws histograms booked\n",fSSDhTask));
  AliInfo(Form("Number of histograms (SPD+SDD+SSD): %d\n",fGenOffset+fSSDhTask));  
  AliInfo(Form("Number of histograms (SPD+SDD+SSD): %d\n",fGenOffset+fSSDRawsOffset));
}


//____________________________________________________________________________
void AliITSQASSDDataMakerRec::MakeRaws(AliRawReader* rawReader) { 
  // Fill QA for RAW - SSD -
  Int_t fStripNumber;
  Int_t fHistPosition;
  Int_t fLayer = 0,fLadder = 0, fModule = 0;
  
  Double_t fSizePerDDL[fgkNumOfDDLs] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t fSizePerLDC[fgkNumOfLDCs] = {0.,0.,0.};
  Double_t sumSSDDataSize = 0.0;
  Double_t eventSize = -1.0;

  if(fkOnline) {
    //reset the signal vs strip number histograms
    for(Int_t i = 0; i < fgkSSDMODULES; i++) 
      fHistSSDRawSignalModule[i]->Reset();
  }
  rawReader->Select("ITSSSD",-1,-1);  
  rawReader->Reset(); //rawReader->NextEvent();   
  (fAliITSQADataMakerRec->GetRawsData(fGenOffset))->Fill(rawReader->GetType());
  if(rawReader->GetType() == 7)
    fSSDEvent += 1;
  
  AliITSRawStreamSSD fSSDStream(rawReader);    
  AliRawReaderRoot *rootReader = (AliRawReaderRoot *)rawReader;
  if(rootReader) {
    const AliRawEventHeaderBase *header = rootReader->GetEventHeader();
    if(header)
      eventSize = header->GetEventSize();
  }
  while (fSSDStream.Next()) {
    if(fSSDStream.GetModuleID() < 0) continue;
    fSizePerDDL[rawReader->GetDDLID()] = rawReader->GetDataSize();
    fSizePerLDC[rawReader->GetLDCId()-6] = rawReader->GetDataSize();
    AliITSgeomTGeo::GetModuleId(fSSDStream.GetModuleID(),fLayer,fLadder,fModule);
    fStripNumber = (fSSDStream.GetSideFlag() == 0) ? fSSDStream.GetStrip() : fSSDStream.GetStrip() + fgkNumberOfPSideStrips;
    fHistPosition = (fLayer == 5) ? ((fLadder - 1)*fgkSSDMODULESPERLADDERLAYER5 + fModule - 1) : ((fLadder - 1)*fgkSSDMODULESPERLADDERLAYER6 + fModule + fgkSSDMODULESLAYER5 - 1);
    //AliInfo(Form("ModulePosition: %d - Layer: %d - Ladder: %d - Module: %d\n",fHistPosition,fLayer,fLadder,fModule));
    if(fkOnline)
      fHistSSDRawSignalModule[fHistPosition]->Fill(fStripNumber,fSSDStream.GetSignal());
  }//streamer loop   
  
  
  for(Int_t i = 0; i < fgkNumOfDDLs; i++) {
    sumSSDDataSize += fSizePerDDL[i];
    if(fSizePerDDL[i] > 0) {
      (fAliITSQADataMakerRec->GetRawsData(fGenOffset+3))->Fill(i+512);
      (fAliITSQADataMakerRec->GetRawsData(fGenOffset+5+i))->Fill(TMath::Log10(fSizePerDDL[i]));
    }
    (fAliITSQADataMakerRec->GetRawsData(fGenOffset+4))->Fill(i+512,fSizePerDDL[i]/1e+06);
  }
  for(Int_t i = 0; i < fgkNumOfLDCs; i++) {
    if(fSizePerLDC[i] > 0) {
      (fAliITSQADataMakerRec->GetRawsData(fGenOffset+21))->Fill(i+6);
      (fAliITSQADataMakerRec->GetRawsData(fGenOffset+23+i))->Fill(TMath::Log10(fSizePerLDC[0]));
    }
    (fAliITSQADataMakerRec->GetRawsData(fGenOffset+22))->Fill(i+6,fSizePerLDC[i]/1e+06);
  }
  if(sumSSDDataSize) 
    (fAliITSQADataMakerRec->GetRawsData(fGenOffset+1))->Fill(TMath::Log10(sumSSDDataSize));
  if(eventSize)
    (fAliITSQADataMakerRec->GetRawsData(fGenOffset+2))->Fill(100.*sumSSDDataSize/eventSize);
  if(fkOnline) {
    //AliInfo(Form("EVENT: %d\n",fSSDEvent));
    Double_t occupancy = 0.0;
    Int_t lLadderLocationY = 0;    
    Int_t fOccPosition = 0;
    for(Int_t i = 0; i < fgkSSDMODULES; i++) {
      AliITSgeomTGeo::GetModuleId(i+500,fLayer,fLadder,fModule); 
      fOccPosition  = (fLayer == 5) ? fLadder : fLadder + fgkSSDLADDERSLAYER5;
      //P-SIDE OCCUPANCY
      occupancy = GetSSDOccupancyRaws(fHistSSDRawSignalModule[i],0);
      fAliITSQADataMakerRec->GetRawsData(fGenOffset + fSSDRawsOffset - 1 + 2*fOccPosition - 1)->Fill(fModule,occupancy);
      lLadderLocationY = 3*fLadder; // sideP=1 sideN=0
      if(fLayer == 5)
	((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenOffset+28))->SetBinContent(fModule,lLadderLocationY,occupancy);
      else if(fLayer == 6)
	((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenOffset+29))->SetBinContent(fModule,lLadderLocationY,occupancy);
      //N-SIDE OCCUPANCY
      occupancy = GetSSDOccupancyRaws(fHistSSDRawSignalModule[i],1);
      fAliITSQADataMakerRec->GetRawsData(fGenOffset + fSSDRawsOffset - 1 + 2*fOccPosition)->Fill(fModule,occupancy);
      if(fLayer == 5)
	((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenOffset+28))->SetBinContent(fModule,lLadderLocationY-1,occupancy);
      else if(fLayer == 6)
	((TH2D *)fAliITSQADataMakerRec->GetRawsData(fGenOffset+29))->SetBinContent(fModule,lLadderLocationY-1,occupancy);
    }
  }//online flag for SSD
}




//____________________________________________________________________________ 
Double_t AliITSQASSDDataMakerRec::GetSSDOccupancyRaws(TH1 *lHisto, Int_t stripside) { 
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
  TH1F *fHistLocalXLayer5 = new TH1F("fHistLocalXLayer5",
				     "Local x coord.- Layer 5;x [cm];Entries;",
				     100,-4.,4.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistLocalXLayer5,
					   fGenOffset + 2);
  fSSDhTask += 1;
  TH1F *fHistLocalXLayer6 = new TH1F("fHistLocalXLayer6",
				     "Local x coord.- Layer 6;x [cm];Entries;",
				     100,-4.,4.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistLocalXLayer6, 
					   fGenOffset + 3);
  fSSDhTask += 1;
  TH1F *fHistLocalZLayer5 = new TH1F("fHistLocalZLayer5",
				     "Local z coord.- Layer 5;z [cm];Entries;",
				     100,-4.,4.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistLocalZLayer5, 
					   fGenOffset + 4);
  fSSDhTask += 1;
  TH1F *fHistLocalZLayer6 = new TH1F("fHistLocalZLayer6",
				     "Local z coord.- Layer 6;z [cm];Entries;",
				     100,-4.,4.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistLocalZLayer6, 
					   fGenOffset + 5);
  fSSDhTask += 1;
  TH1F *fHistGlobalXLayer5 = new TH1F("fHistGlobalXLayer5",
				      "Global x - Layer 5;x [cm];Entries;",
				      100,-40.,40.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistGlobalXLayer5, 
					   fGenOffset + 6);
  fSSDhTask += 1;
  TH1F *fHistGlobalXLayer6 = new TH1F("fHistGlobalXLayer6",
				      "Global x - Layer 6;x [cm];Entries;",
				      100,-45.,45.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistGlobalXLayer6, 
					   fGenOffset + 7);
  fSSDhTask += 1;
  TH1F *fHistGlobalYLayer5 = new TH1F("fHistGlobalYLayer5",
				      "Global y - Layer 5;y [cm];Entries;",
				      100,-40.,40);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistGlobalYLayer5, 
					   fGenOffset + 8);
  fSSDhTask += 1;
  TH1F *fHistGlobalYLayer6 = new TH1F("fHistGlobalYLayer6",
				      "Global y - Layer 6;y [cm];Entries;",
				      100,-45.,45.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistGlobalYLayer6, 
					   fGenOffset + 9);
  fSSDhTask += 1;
  TH1F *fHistGlobalZLayer5 = new TH1F("fHistGlobalZLayer5",
				      "Global z - Layer 5;z [cm];Entries;",
				      100,-45.,45);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistGlobalZLayer5, 
					   fGenOffset + 10);
  fSSDhTask += 1;
  TH1F *fHistGlobalZLayer6 = new TH1F("fHistGlobalZLayer6",
				      "Global z - Layer 6;z [cm];Entries;",
				      100,-55.,55.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistGlobalZLayer6, 
					   fGenOffset + 11);
  fSSDhTask += 1;
  TH1F *fHistPhiLayer5 = new TH1F("fHistPhiLayer5",
				  "#phi - Layer 5;#phi [rad];Entries;",
				  100,-TMath::Pi(),TMath::Pi());
  fAliITSQADataMakerRec->Add2RecPointsList(fHistPhiLayer5, 
					   fGenOffset + 12);
  fSSDhTask += 1;
  TH1F *fHistPhiLayer6 = new TH1F("fHistPhiLayer6",
				  "#phi - Layer 6;#phi [rad];Entries;",
				  100,-TMath::Pi(),TMath::Pi());
  fAliITSQADataMakerRec->Add2RecPointsList(fHistPhiLayer6, 
					   fGenOffset + 13);
  fSSDhTask += 1;
  TH1F *fHistThetaLayer5 = new TH1F("fHistThetaLayer5",
				    "#theta - Layer 5;#theta [rad];Entries;",
				    100,-TMath::Pi(),TMath::Pi());
  fAliITSQADataMakerRec->Add2RecPointsList(fHistThetaLayer5, 
					   fGenOffset + 14);
  fSSDhTask += 1;
  TH1F *fHistThetaLayer6 = new TH1F("fHistThetaLayer6",
				    "#theta - Layer 6;#theta [rad];Entries;",
				    100,-TMath::Pi(),TMath::Pi());
  fAliITSQADataMakerRec->Add2RecPointsList(fHistThetaLayer6, 
					   fGenOffset + 15);
  fSSDhTask += 1;
  TH1F *fHistRadiusLayer5 = new TH1F("fHistRadiusLayer5",
				     "r - Layer 5;r [cm];Entries;",
				     100,35.,50.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistRadiusLayer5, 
					   fGenOffset + 16);
  fSSDhTask += 1;
  TH1F *fHistRadiusLayer6 = new TH1F("fHistRadiusLayer6",
				     "r - Layer 6;r [cm];Entries;",
				     100,35.,50.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistRadiusLayer6, 
					   fGenOffset + 17);
  fSSDhTask += 1;
  TH1F *fHistClusterTypeLayer5 = new TH1F("fHistClusterTypeLayer5",
					  "CL type - Layer 5;Cluster type;Entries;",
					  150,0,150);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistClusterTypeLayer5, 
					   fGenOffset + 18);
  fSSDhTask += 1;
  TH1F *fHistClusterTypeLayer6 = new TH1F("fHistClusterTypeLayer6",
					  "CL type - Layer 6;Cluster type;Entries;",
					  150,0,150);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistClusterTypeLayer6, 
					   fGenOffset + 19);
  fSSDhTask += 1;
  TH1F *fHistChargeRatioLayer5 = new TH1F("fHistChargeRatioLayer5",
					  "Charge ratio - Layer 5;q_{ratio};Entries;",
					  100,-2.0,2.0);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargeRatioLayer5, 
					   fGenOffset + 20);
  fSSDhTask += 1;
  TH1F *fHistChargeRatioLayer6 = new TH1F("fHistChargeRatioLayer6",
					  "Charge ratio - Layer 6;q_{ratio};Entries;",
					  100,-2.0,2.0);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargeRatioLayer6, 
					   fGenOffset + 21);
  fSSDhTask += 1;
  TH1F *fHistChargekeVLayer5 = new TH1F("fHistChargekeVLayer5",
					"Charge - Layer 5;q [keV];Entries;",
					100,0.,300.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargekeVLayer5, 
					   fGenOffset + 22);
  fSSDhTask += 1;
  TH1F *fHistChargekeVLayer6 = new TH1F("fHistChargekeVLayer6",
					"Charge - Layer 6;q [keV];Entries;",
					100,0.,300.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargekeVLayer6, 
					   fGenOffset + 23);
  fSSDhTask += 1;
  TH1F *fHistChargeADCLayer5 = new TH1F("fHistChargeADCLayer5",
					"Charge - Layer 5;q [ADC];Entries;",
					100,0.,300.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargeADCLayer5, 
					   fGenOffset + 24);
  fSSDhTask += 1;
  TH1F *fHistChargeADCLayer6 = new TH1F("fHistChargeADCLayer6",
					"Charge - Layer 6;q [ADC];Entries;",
					100,0.,300.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargeADCLayer6, 
					   fGenOffset + 25);
  fSSDhTask += 1;
  TH2F *fHistChargeMapLayer5 = new TH2F("fHistChargeMapLayer5",
					"Charge map;N_{modules};N_{Ladders}",
					fgkSSDMODULESPERLADDERLAYER5,
					-0.5,fgkSSDMODULESPERLADDERLAYER5+0.5,
					3*fgkSSDLADDERSLAYER5,
					-0.5,fgkSSDLADDERSLAYER5+0.5);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargeMapLayer5, 
					   fGenOffset + 26);
  fSSDhTask += 1;
  TH2F *fHistChargeMapLayer6 = new TH2F("fHistChargeMapLayer6",
					"Charge map;N_{modules};N_{Ladders}",
					fgkSSDMODULESPERLADDERLAYER6,
					-0.5,fgkSSDMODULESPERLADDERLAYER6+0.5,
					3*fgkSSDLADDERSLAYER6,
					-0.5,fgkSSDLADDERSLAYER6+0.5);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargeMapLayer6, 
					   fGenOffset + 27);
  fSSDhTask += 1;

// custom code here
//fSSDhTask must be incremented by one unit every time a histogram is ADDED to the QA List

  AliDebug(1,Form("%d SSD Recs histograms booked\n",fSSDhTask));
}

//____________________________________________________________________________ 
void AliITSQASSDDataMakerRec::MakeRecPoints(TTree *clustersTree)
{
  // Fill QA for recpoints - SSD -
  Int_t fLayer = 0, fLadder = 0, fModule = 0;
  Int_t lLadderLocationY = 0;
  TBranch *branchRecP = clustersTree->GetBranch("ITSRecPoints");
  if (!branchRecP) { 
    AliError("can't get the branch with the ITS clusters !");
    return;
  }
  static TClonesArray statRecpoints("AliITSRecPoint");
  TClonesArray *recpoints = &statRecpoints;
  branchRecP->SetAddress(&recpoints);
  Int_t npoints = 0;      
  Float_t cluglo[3]={0.,0.,0.}; 
  for(Int_t module = 0; module < clustersTree->GetEntries(); module++){
    branchRecP->GetEvent(module);
    npoints += recpoints->GetEntries();
    AliITSgeomTGeo::GetModuleId(module,fLayer,fLadder,fModule);
    lLadderLocationY = 3*fLadder;

    for(Int_t j = 0;j < recpoints->GetEntries(); j++){
      AliITSRecPoint *recp = (AliITSRecPoint*)recpoints->At(j);    
      Int_t layer = recp->GetLayer();
      recp->GetGlobalXYZ(cluglo);
      Float_t radius = TMath::Sqrt(cluglo[0]*cluglo[0]+cluglo[1]*cluglo[1]); 
      Float_t phi = TMath::ATan2(cluglo[1],cluglo[0]);
      Float_t theta = TMath::ATan2(radius,cluglo[2]);
      if(layer == 4) {
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 0)->Fill(module);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 2)->Fill(recp->GetDetLocalX());
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 4)->Fill(recp->GetDetLocalZ());
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 6)->Fill(cluglo[0]);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 8)->Fill(cluglo[1]);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 10)->Fill(cluglo[2]);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 12)->Fill(phi);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 14)->Fill(theta);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 16)->Fill(radius);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 18)->Fill(recp->GetType());
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 20)->Fill(recp->GetChargeRatio());
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 22)->Fill(recp->GetQ());
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 26)->SetBinContent(fModule,lLadderLocationY,recp->GetQ());
      }//layer 5 histograms
      if(layer == 5) {
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 1)->Fill(module);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 3)->Fill(recp->GetDetLocalX());
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 5)->Fill(recp->GetDetLocalZ());
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 7)->Fill(cluglo[0]);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 9)->Fill(cluglo[1]);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 11)->Fill(cluglo[2]);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 13)->Fill(phi);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 15)->Fill(theta);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 17)->Fill(radius);
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 19)->Fill(recp->GetType());
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 21)->Fill(recp->GetChargeRatio());
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 23)->Fill(recp->GetQ());
	fAliITSQADataMakerRec->GetRecPointsData(fGenOffset + 27)->SetBinContent(fModule,lLadderLocationY,recp->GetQ());
      }//layer 6 histograms
    }//rec. points loop
  }//module loop

  statRecpoints.Clear();
}

//____________________________________________________________________________ 
/*void AliITSQASSDDataMakerRec::InitESDs() {  
  // Initialization for ESD - SSD -
  fESDsOffset = (fAliITSQADataMakerRec->fESDsQAList)->GetEntries();
  AliDebug(1,Form("Number of ESD histograms (SPD+SDD): %d\n",fRawsOffset));
  
  TH1F *fHistSSDTrackPt = new TH1F("fHistSSDTrackPt",
				   ";P_{T} [GeV/c];dN/dP_{T}",
				   100,0.1,10.1);
  fAliITSQADataMakerRec->Add2ESDPointsList(fHistSSDTrackPt, 
					   fESDsOffset + 0);
  fSSDhTask += 1;
  TH1F *fHistSSDTrackEta = new TH1F("fHistSSDTrackEta",
				   ";#eta;dN/d#eta",
				   40,-2.,2.);
  fAliITSQADataMakerRec->Add2ESDPointsList(fHistSSDTrackEta,
					   fESDsOffset + 1);
  fSSDhTask += 1;
  TH1F *fHistSSDTrackPhi = new TH1F("fHistSSDTrackPhi",
				   ";#phi;dN/d#phi",
				   100,0,2.*TMath::Pi());
  fAliITSQADataMakerRec->Add2ESDPointsList(fHistSSDTrackPhi,
					   fESDsOffset + 2);
  fSSDhTask += 1;
  
  AliDebug(1,Form("%d SSD ESDs histograms booked\n",fSSDhTask));
  AliInfo(Form("Number of histograms (ITS): %d\n",fESDsOffset+fSSDhTask));  
}

//____________________________________________________________________________
void AliITSQASSDDataMakerRec::MakeESDs(AliESDEvent * esd) {
  // make QA data from ESDs
  Int_t clusters[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  for(Int_t j = 0; j < esd->GetNumberOfTracks(); j++) {
    AliESDtrack *track = esd->GetTrack(j);
    if ((track->GetStatus() & AliESDtrack::kITSrefit)==0) continue;
    Int_t nClustersITS = track->GetITSclusters(clusters);

   (fAliITSQADataMakerRec->GetESDsData(fESDsOffset+0))->Fill(track->Pt());
   (fAliITSQADataMakerRec->GetESDsData(fESDsOffset+1))->Fill(track->Eta());
   (fAliITSQADataMakerRec->GetESDsData(fESDsOffset+2))->Fill(track->Phi());
  }//track loop
}
*/
