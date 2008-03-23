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
/* $Id:$    */
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

ClassImp(AliITSQASSDDataMakerRec)

AliITSQASSDDataMakerRec::AliITSQASSDDataMakerRec(AliITSQADataMakerRec *aliITSQADataMakerRec, Bool_t kMode, Int_t ldc) :
TObject(),
fAliITSQADataMakerRec(aliITSQADataMakerRec),
fSSDEvent(0),
fkOnline(kMode),
fLDC(ldc),
fSSDRawsOffset(0),
fSSDhRaws(0),
fSSDhRecs(0),
fRawsOffset(0),
fRecsOffset(0) {
  // Default constructor   
  //initilize the raw signal vs strip number histograms
  Int_t fLayer = 0,fLadder = 0, fModule = 0;
  TString fTitle; Int_t fHistCounter = 0;
  if(fkOnline) {
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
fSSDhRaws(qadm.fSSDhRaws),
fSSDhRecs(qadm.fSSDhRecs),
fRawsOffset(qadm.fRawsOffset),
fRecsOffset(qadm.fRecsOffset) {
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
void AliITSQASSDDataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX_t /*task*/, TObjArray* /*list*/)
{
  // launch the QA checking
  AliDebug(1,"AliITSDM instantiates checker with Run(AliQA::kITS, task, list)\n"); 
  
  //scaling SSD occupancy plots
  if(fkOnline) {
    if(fSSDEvent != 0) {
      ((TH2F *)fAliITSQADataMakerRec->GetRawsData(fRawsOffset+28))->Scale(1./fSSDEvent);
      ((TH2F *)fAliITSQADataMakerRec->GetRawsData(fRawsOffset+29))->Scale(1./fSSDEvent);
    }
    Int_t fLayer = 0, fLadder = 0;
    Int_t fOccPosition = 0; 
    for(Int_t iLayer = 5; iLayer < 7; iLayer++) {
      for(Int_t iLadder = 1; iLadder < AliITSgeomTGeo::GetNLadders(iLayer) + 1; iLadder++) {
	fOccPosition  = (fLayer == 5) ? fLadder : fLadder + fgkSSDLADDERSLAYER5;
	if(fSSDEvent != 0) {
	  //P-SIDE OCCUPANCY - scaling
	  fAliITSQADataMakerRec->GetRawsData(fRawsOffset + fSSDRawsOffset - 1 + 2*fOccPosition - 1)->Scale(1./fSSDEvent);
	  //N-SIDE OCCUPANCY - scaling
	  fAliITSQADataMakerRec->GetRawsData(fRawsOffset + fSSDRawsOffset - 1 + 2*fOccPosition)->Scale(1./fSSDEvent);
	}//ssd events != 0
      }//ladder loop
    }//layer loop
  }//online flag for SSD
  //AliQAChecker::Instance()->Run( AliQA::kITS , task, list);
}

//____________________________________________________________________________ 
void AliITSQASSDDataMakerRec::InitRaws() {  

  // Initialization for RAW data - SSD -
  fRawsOffset = (fAliITSQADataMakerRec->fRawsQAList)->GetEntries();

  if(fkOnline) {
    AliInfo("Book Online Histograms for SSD\n");
  }
  else {
    AliInfo("Book Offline Histograms for SSD\n ");
  }
  AliInfo(Form("Number of histograms (SPD+SDD): %d\n",fRawsOffset));
  TString fTitle = 0;
  //book online QA histos
  if(fkOnline) {
    TH1F *fHistSSDEventType = new TH1F("fHistSSDEventType",
				       ";Event type;Events",
				       31,-1,30);
    fAliITSQADataMakerRec->Add2RawsList(fHistSSDEventType, 
					fRawsOffset+fSSDRawsOffset);
    fSSDRawsOffset += 1;
    TH1F *fHistSSDDataSize = new TH1F("fHistSSDDataSize",
				      ";log(SSD data size) [Bytes];Events",
				      100,3,8);
    fAliITSQADataMakerRec->Add2RawsList(fHistSSDDataSize, 
					fRawsOffset+fSSDRawsOffset);
    fSSDRawsOffset += 1;
    TH1F *fHistSSDDataSizePercentage = new TH1F("fHistSSDDataSizePercentage",
						";SSD data size [%];Events",
						100,0,100);
    fAliITSQADataMakerRec->Add2RawsList(fHistSSDDataSizePercentage, 
					fRawsOffset+fSSDRawsOffset);
    fSSDRawsOffset += 1;
    TH1F *fHistSSDDDLId = new TH1F("fHistSSDDDLId",
				   ";DDL id;Events",20,510.5,530.5);
    fAliITSQADataMakerRec->Add2RawsList(fHistSSDDDLId, 
					fRawsOffset+fSSDRawsOffset);
    fSSDRawsOffset += 1;
    TH1F *fHistSSDDataSizePerDDL = new TH1F("fHistSSDDataSizePerDDL",
					    ";DDL id;<SSD data size> [MB]",
					    20,510.5,530.5);
    fAliITSQADataMakerRec->Add2RawsList(fHistSSDDataSizePerDDL, 
					fRawsOffset+fSSDRawsOffset);
    fSSDRawsOffset += 1;
    TH1F *fHistSSDDataSizeDDL[fgkNumOfDDLs];
    for(Int_t i = 1; i < fgkNumOfDDLs+1; i++) {
      fTitle = "fHistSSDDataSizeDDL"; fTitle += i+511;
      fHistSSDDataSizeDDL[i-1] = new TH1F(fTitle.Data(),
					  ";log(SSD data size) [Bytes];Events",
					  100,1,8);
      fAliITSQADataMakerRec->Add2RawsList(fHistSSDDataSizeDDL[i-1], 
					  fRawsOffset+fSSDRawsOffset);
      fSSDRawsOffset += 1;
    }
    
    TH1F *fHistSSDLDCId = new TH1F("fHistSSDLDCId",";LDC id;Events",10,0.5,10.5);
    fAliITSQADataMakerRec->Add2RawsList(fHistSSDLDCId, 
					fRawsOffset+fSSDRawsOffset);
    fSSDRawsOffset += 1;
    TH1F *fHistSSDDataSizePerLDC = new TH1F("fHistSSDDataSizePerLDC",
					    ";LDC id;<SSD data size> [MB]",
					    100,0,20);
    fAliITSQADataMakerRec->Add2RawsList(fHistSSDDataSizePerLDC, 
					fRawsOffset+fSSDRawsOffset);
    fSSDRawsOffset += 1;
    TH1F *fHistSSDDataSizeLDC[fgkNumOfLDCs];
    for(Int_t i = 1; i < fgkNumOfLDCs+1; i++) {
      fTitle = "fHistSSDDataSizeLDC"; fTitle += i;
      fHistSSDDataSizeLDC[i-1] = new TH1F(fTitle.Data(),
					  ";log(SSD data size) [Bytes];Events",
					  100,1,8);
      fAliITSQADataMakerRec->Add2RawsList(fHistSSDDataSizeLDC[i-1], 
					  fRawsOffset+fSSDRawsOffset);
      fSSDRawsOffset += 1;
    }

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
    fAliITSQADataMakerRec->Add2RawsList(fHistSSDOccupancyLayer5, fRawsOffset+fSSDRawsOffset);
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
    fAliITSQADataMakerRec->Add2RawsList(fHistSSDOccupancyLayer6, fRawsOffset+fSSDRawsOffset);
    fSSDRawsOffset += 1;

    //Occupancy per ladder
    fSSDhRaws = fSSDRawsOffset;
    TH1D *fHistOccupancyLadder[2*(fgkSSDLADDERSLAYER5 + fgkSSDLADDERSLAYER6)];
    for(Int_t iLayer = 5; iLayer < 7; iLayer++) {
      for(Int_t iLadder = 1; iLadder < AliITSgeomTGeo::GetNLadders(iLayer) + 1; iLadder++) {
	//P-side occupancy plots
	fTitle = "fHistSSD_Occupancy_Layer"; fTitle += iLayer;
	fTitle += "_Ladder"; fTitle += iLadder; fTitle += "_PSide";
	fHistOccupancyLadder[fSSDhRaws] = 
	  new TH1D(fTitle.Data(),
		   fTitle.Data(),
		   AliITSgeomTGeo::GetNDetectors(iLayer),
		   0.5,AliITSgeomTGeo::GetNDetectors(iLayer)+0.5);
	fHistOccupancyLadder[fSSDhRaws]->GetXaxis()->SetTitleColor(1);
	fHistOccupancyLadder[fSSDhRaws]->GetXaxis()->SetTitle("Module number");
	fHistOccupancyLadder[fSSDhRaws]->GetYaxis()->SetTitle("Occupancy [%]");
	fAliITSQADataMakerRec->Add2RawsList(fHistOccupancyLadder[fSSDhRaws], 
					    fRawsOffset+fSSDhRaws);	
	fSSDhRaws++;
	//N-side occupancy plots
	fTitle = "fHistSSD_Occupancy_Layer"; fTitle += iLayer;
	fTitle += "_Ladder"; fTitle += iLadder; fTitle += "_NSide";
	fHistOccupancyLadder[fSSDhRaws] = 
	  new TH1D(fTitle.Data(),
		   fTitle.Data(),
		   AliITSgeomTGeo::GetNDetectors(iLayer),
		   0.5,AliITSgeomTGeo::GetNDetectors(iLayer)+0.5);
	fHistOccupancyLadder[fSSDhRaws]->GetXaxis()->SetTitleColor(1);
	fHistOccupancyLadder[fSSDhRaws]->GetXaxis()->SetTitle("Module number");
	fHistOccupancyLadder[fSSDhRaws]->GetYaxis()->SetTitle("Occupancy [%]");
	fAliITSQADataMakerRec->Add2RawsList(fHistOccupancyLadder[fSSDhRaws], 
					    fRawsOffset+fSSDhRaws);	
	fSSDhRaws++;
      }//ladder loop
    }//layer loop
  }//online flag
  AliDebug(1,Form("%d SSD Raws histograms booked\n",fSSDhRaws));
  AliInfo(Form("Number of histograms (SPD+SDD+SSD): %d\n",fRawsOffset+fSSDhRaws));  
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
  
    rawReader->Select("ITSSSD",-1,-1);  
    rawReader->Reset(); //rawReader->NextEvent();   
    (fAliITSQADataMakerRec->GetRawsData(fRawsOffset))->Fill(rawReader->GetType());
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
      fHistSSDRawSignalModule[fHistPosition]->Fill(fStripNumber,fSSDStream.GetSignal());
    }//streamer loop   

    
    for(Int_t i = 0; i < fgkNumOfDDLs; i++) {
      sumSSDDataSize += fSizePerDDL[i];
      if(fSizePerDDL[i] > 0) {
	(fAliITSQADataMakerRec->GetRawsData(fRawsOffset+3))->Fill(i+512);
	(fAliITSQADataMakerRec->GetRawsData(fRawsOffset+5+i))->Fill(TMath::Log10(fSizePerDDL[i]));
      }
      (fAliITSQADataMakerRec->GetRawsData(fRawsOffset+4))->Fill(i+512,fSizePerDDL[i]/1e+06);
    }
    for(Int_t i = 0; i < fgkNumOfLDCs; i++) {
      if(fSizePerLDC[i] > 0) {
	(fAliITSQADataMakerRec->GetRawsData(fRawsOffset+21))->Fill(i+6);
	(fAliITSQADataMakerRec->GetRawsData(fRawsOffset+23+i))->Fill(TMath::Log10(fSizePerLDC[0]));
      }
      (fAliITSQADataMakerRec->GetRawsData(fRawsOffset+22))->Fill(i+6,fSizePerLDC[i]/1e+06);
    }
    if(sumSSDDataSize) 
      (fAliITSQADataMakerRec->GetRawsData(fRawsOffset+1))->Fill(TMath::Log10(sumSSDDataSize));
    if(eventSize)
      (fAliITSQADataMakerRec->GetRawsData(fRawsOffset+2))->Fill(100.*sumSSDDataSize/eventSize);
 
    //AliInfo(Form("EVENT: %d\n",fSSDEvent));
    Double_t occupancy = 0.0;
    Int_t lLadderLocationY = 0;    
    Int_t fOccPosition = 0;
    for(Int_t i = 0; i < fgkSSDMODULES; i++) {
      AliITSgeomTGeo::GetModuleId(i+500,fLayer,fLadder,fModule); 
      fOccPosition  = (fLayer == 5) ? fLadder : fLadder + fgkSSDLADDERSLAYER5;
      //P-SIDE OCCUPANCY
      occupancy = GetSSDOccupancyRaws(fHistSSDRawSignalModule[i],0);
      fAliITSQADataMakerRec->GetRawsData(fRawsOffset + fSSDRawsOffset - 1 + 2*fOccPosition - 1)->Fill(fModule,occupancy);
      lLadderLocationY = 3*fLadder; // sideP=1 sideN=0
      if(fLayer == 5)
	((TH2D *)fAliITSQADataMakerRec->GetRawsData(fRawsOffset+28))->SetBinContent(fModule,lLadderLocationY,occupancy);
      else if(fLayer == 6)
	((TH2D *)fAliITSQADataMakerRec->GetRawsData(fRawsOffset+29))->SetBinContent(fModule,lLadderLocationY,occupancy);
      //N-SIDE OCCUPANCY
      occupancy = GetSSDOccupancyRaws(fHistSSDRawSignalModule[i],1);
      fAliITSQADataMakerRec->GetRawsData(fRawsOffset + fSSDRawsOffset - 1 + 2*fOccPosition)->Fill(fModule,occupancy);
      if(fLayer == 5)
	((TH2D *)fAliITSQADataMakerRec->GetRawsData(fRawsOffset+28))->SetBinContent(fModule,lLadderLocationY-1,occupancy);
      else if(fLayer == 6)
	((TH2D *)fAliITSQADataMakerRec->GetRawsData(fRawsOffset+29))->SetBinContent(fModule,lLadderLocationY-1,occupancy);
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
  fRecsOffset = (fAliITSQADataMakerRec->fRecPointsQAList)->GetEntries();
  Int_t nModuleOffset = 500;
  Int_t nITSTotalModules = AliITSgeomTGeo::GetNModules();

  TH1F *fHistModuleIdLayer5 = new TH1F("fHistModuleIdLayer5",
				       "Module Id - Layer 5;Module Id;Entries",
				       fgkSSDMODULESLAYER5,
				       nModuleOffset - 0.5,
				       nITSTotalModules-fgkSSDMODULESLAYER6+0.5);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistModuleIdLayer5, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistModuleIdLayer6 = new TH1F("fHistModuleIdLayer6",
				       "Module Id - Layer 6;Module Id;Entries",
				       fgkSSDMODULESLAYER6,
				       nModuleOffset+fgkSSDMODULESLAYER5 - 0.5,
				       nITSTotalModules + 0.5);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistModuleIdLayer6, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistLocalXLayer5 = new TH1F("fHistLocalXLayer5",
				     "Local x coord.- Layer 5;x [cm];Entries;",
				     100,-4.,4.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistLocalXLayer5,
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistLocalXLayer6 = new TH1F("fHistLocalXLayer6",
				     "Local x coord.- Layer 6;x [cm];Entries;",
				     100,-4.,4.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistLocalXLayer6, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistLocalZLayer5 = new TH1F("fHistLocalZLayer5",
				     "Local z coord.- Layer 5;z [cm];Entries;",
				     100,-4.,4.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistLocalZLayer5, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistLocalZLayer6 = new TH1F("fHistLocalZLayer6",
				     "Local z coord.- Layer 6;z [cm];Entries;",
				     100,-4.,4.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistLocalZLayer6, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistGlobalXLayer5 = new TH1F("fHistGlobalXLayer5",
				      "Global x - Layer 5;x [cm];Entries;",
				      100,-40.,40.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistGlobalXLayer5, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistGlobalXLayer6 = new TH1F("fHistGlobalXLayer6",
				      "Global x - Layer 6;x [cm];Entries;",
				      100,-45.,45.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistGlobalXLayer6, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistGlobalYLayer5 = new TH1F("fHistGlobalYLayer5",
				      "Global y - Layer 5;y [cm];Entries;",
				      100,-40.,40);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistGlobalYLayer5, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistGlobalYLayer6 = new TH1F("fHistGlobalYLayer6",
				      "Global y - Layer 6;y [cm];Entries;",
				      100,-45.,45.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistGlobalYLayer6, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistGlobalZLayer5 = new TH1F("fHistGlobalZLayer5",
				      "Global z - Layer 5;z [cm];Entries;",
				      100,-45.,45);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistGlobalZLayer5, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistGlobalZLayer6 = new TH1F("fHistGlobalZLayer6",
				      "Global z - Layer 6;z [cm];Entries;",
				      100,-55.,55.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistGlobalZLayer6, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistPhiLayer5 = new TH1F("fHistPhiLayer5",
				  "#phi - Layer 5;#phi [rad];Entries;",
				  100,-TMath::Pi(),TMath::Pi());
  fAliITSQADataMakerRec->Add2RecPointsList(fHistPhiLayer5, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistPhiLayer6 = new TH1F("fHistPhiLayer6",
				  "#phi - Layer 6;#phi [rad];Entries;",
				  100,-TMath::Pi(),TMath::Pi());
  fAliITSQADataMakerRec->Add2RecPointsList(fHistPhiLayer6, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistThetaLayer5 = new TH1F("fHistThetaLayer5",
				    "#theta - Layer 5;#theta [rad];Entries;",
				    100,-TMath::Pi(),TMath::Pi());
  fAliITSQADataMakerRec->Add2RecPointsList(fHistThetaLayer5, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistThetaLayer6 = new TH1F("fHistThetaLayer6",
				    "#theta - Layer 6;#theta [rad];Entries;",
				    100,-TMath::Pi(),TMath::Pi());
  fAliITSQADataMakerRec->Add2RecPointsList(fHistThetaLayer6, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistRadiusLayer5 = new TH1F("fHistRadiusLayer5",
				     "r - Layer 5;r [cm];Entries;",
				     100,35.,50.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistRadiusLayer5, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistRadiusLayer6 = new TH1F("fHistRadiusLayer6",
				     "r - Layer 6;r [cm];Entries;",
				     100,35.,50.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistRadiusLayer6, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistClusterTypeLayer5 = new TH1F("fHistClusterTypeLayer5",
					  "CL type - Layer 5;Cluster type;Entries;",
					  150,0,150);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistClusterTypeLayer5, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistClusterTypeLayer6 = new TH1F("fHistClusterTypeLayer6",
					  "CL type - Layer 6;Cluster type;Entries;",
					  150,0,150);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistClusterTypeLayer6, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistChargeRatioLayer5 = new TH1F("fHistChargeRatioLayer5",
					  "Charge ratio - Layer 5;q_{ratio};Entries;",
					  100,-2.0,2.0);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargeRatioLayer5, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistChargeRatioLayer6 = new TH1F("fHistChargeRatioLayer6",
					  "Charge ratio - Layer 6;q_{ratio};Entries;",
					  100,-2.0,2.0);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargeRatioLayer6, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistChargekeVLayer5 = new TH1F("fHistChargekeVLayer5",
					"Charge - Layer 5;q [keV];Entries;",
					100,0.,300.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargekeVLayer5, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistChargekeVLayer6 = new TH1F("fHistChargekeVLayer6",
					"Charge - Layer 6;q [keV];Entries;",
					100,0.,300.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargekeVLayer6, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistChargeADCLayer5 = new TH1F("fHistChargeADCLayer5",
					"Charge - Layer 5;q [ADC];Entries;",
					100,0.,300.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargeADCLayer5, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH1F *fHistChargeADCLayer6 = new TH1F("fHistChargeADCLayer6",
					"Charge - Layer 6;q [ADC];Entries;",
					100,0.,300.);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargeADCLayer6, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH2F *fHistChargeMapLayer5 = new TH2F("fHistChargeMapLayer5",
					"Charge map;N_{modules};N_{Ladders}",
					fgkSSDMODULESPERLADDERLAYER5,
					-0.5,fgkSSDMODULESPERLADDERLAYER5+0.5,
					3*fgkSSDLADDERSLAYER5,
					-0.5,fgkSSDLADDERSLAYER5+0.5);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargeMapLayer5, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;
  TH2F *fHistChargeMapLayer6 = new TH2F("fHistChargeMapLayer6",
					"Charge map;N_{modules};N_{Ladders}",
					fgkSSDMODULESPERLADDERLAYER6,
					-0.5,fgkSSDMODULESPERLADDERLAYER6+0.5,
					3*fgkSSDLADDERSLAYER6,
					-0.5,fgkSSDLADDERSLAYER6+0.5);
  fAliITSQADataMakerRec->Add2RecPointsList(fHistChargeMapLayer6, 
					   fRecsOffset + fSSDhRecs);
  fSSDhRecs += 1;

// custom code here
//fSSDhRecs must be incremented by one unit every time a histogram is ADDED to the QA List

  AliDebug(1,Form("%d SSD Recs histograms booked\n",fSSDhRecs));
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
  TClonesArray * recpoints = new TClonesArray("AliITSRecPoint") ;
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
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 0)->Fill(module);
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 2)->Fill(recp->GetDetLocalX());
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 4)->Fill(recp->GetDetLocalZ());
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 6)->Fill(cluglo[0]);
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 8)->Fill(cluglo[1]);
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 10)->Fill(cluglo[2]);
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 12)->Fill(phi);
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 14)->Fill(theta);
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 16)->Fill(radius);
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 18)->Fill(recp->GetType());
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 20)->Fill(recp->GetChargeRatio());
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 22)->Fill(recp->GetQ());
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 26)->SetBinContent(fModule,lLadderLocationY,recp->GetQ());
      }//layer 5 histograms
      if(layer == 5) {
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 1)->Fill(module);
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 3)->Fill(recp->GetDetLocalX());
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 5)->Fill(recp->GetDetLocalZ());
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 7)->Fill(cluglo[0]);
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 9)->Fill(cluglo[1]);
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 11)->Fill(cluglo[2]);
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 13)->Fill(phi);
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 15)->Fill(theta);
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 17)->Fill(radius);
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 19)->Fill(recp->GetType());
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 21)->Fill(recp->GetChargeRatio());
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 23)->Fill(recp->GetQ());
	fAliITSQADataMakerRec->GetRecPointsData(fRecsOffset + 27)->SetBinContent(fModule,lLadderLocationY,recp->GetQ());
      }//layer 6 histograms
    }//rec. points loop
  }//module loop

  recpoints->Delete();
  delete recpoints;
}
