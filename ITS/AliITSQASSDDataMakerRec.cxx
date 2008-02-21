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

ClassImp(AliITSQASSDDataMakerRec)

AliITSQASSDDataMakerRec::AliITSQASSDDataMakerRec(AliITSQADataMakerRec *aliITSQADataMakerRec, Bool_t kMode, Int_t ldc) :
TObject(),
fAliITSQADataMakerRec(aliITSQADataMakerRec),
fSSDEvent(0),
fkOnline(kMode),
fLDC(ldc),
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
void AliITSQASSDDataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX /*task*/, TObjArray* /*list*/)
{
  // launch the QA checking
  AliDebug(1,"AliITSDM instantiates checker with Run(AliQA::kITS, task, list)\n"); 
  
  //SSD occupancy plots
  if(fkOnline) {
    Int_t fLayer = 0, fLadder = 0;
    Int_t fOccPosition = 0; 
    for(Int_t iLayer = 5; iLayer < 7; iLayer++) {
      if(fSSDEvent != 0)
	((TH2F *)fAliITSQADataMakerRec->GetRawsData(fLayer))->Scale(1./fSSDEvent);
      for(Int_t iLadder = 1; iLadder < AliITSgeomTGeo::GetNLadders(iLayer) + 1; iLadder++) {
	fOccPosition  = (fLayer == 5) ? fLadder : fLadder + fgkSSDLADDERSLAYER5;
	if(fSSDEvent != 0) {
	  //P-SIDE OCCUPANCY - scaling
	  fAliITSQADataMakerRec->GetRawsData(fRawsOffset -1 + 2*fLadder - 1)->Scale(1./fSSDEvent);

	  //N-SIDE OCCUPANCY - scaling
	  fAliITSQADataMakerRec->GetRawsData(fRawsOffset -1 + 2*fLadder)->Scale(1./fSSDEvent);
	}
      }      
    }
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
    //data size plots
    TH1F *hSSDDataSize = new TH1F("hSSDDataSize",";log(SSD data size) [Bytes];Events",100,3,8);
    fAliITSQADataMakerRec->Add2RawsList(hSSDDataSize, 0);
    TH1F *hSSDDataSizePercentage = new TH1F("hSSDDataSizePercentage",";SSD data size [%];Events",100,0,100);
    fAliITSQADataMakerRec->Add2RawsList(hSSDDataSizePercentage, 1);
    TH1F *hDDLId = new TH1F("hDDLId",";DDL id;Events",100,510,530);
    fAliITSQADataMakerRec->Add2RawsList(hDDLId, 2);
    TH1F *hSSDDataSizePerDDL = new TH1F("hSSDDataSizePerDDL",";DDL id;<SSD data size> [MB]",100,510,530);
    fAliITSQADataMakerRec->Add2RawsList(hSSDDataSizePerDDL, 3);
    TH1F *hSSDDataSizePerLDC = new TH1F("hSSDDataSizePerLDC",";LDC id;<SSD data size> [MB]",100,0,20);
    fAliITSQADataMakerRec->Add2RawsList(hSSDDataSizePerLDC, 4);

    //top level occupancy plots
    TH2D *hSSDOccupancyLayer5 = new TH2D("hSSDOccupancyLayer5",";N_{modules};N_{Ladders}",fgkSSDMODULESPERLADDERLAYER5,0,fgkSSDMODULESPERLADDERLAYER5,3*fgkSSDLADDERSLAYER5,0,fgkSSDLADDERSLAYER5);  
    Char_t fLabel[3];
    for(Int_t iBin = 1; iBin < fgkSSDMODULESPERLADDERLAYER5 + 1; iBin++){
      sprintf(fLabel,"%d",iBin);
      hSSDOccupancyLayer5->GetXaxis()->SetBinLabel(iBin,fLabel);
    }
    fAliITSQADataMakerRec->Add2RawsList(hSSDOccupancyLayer5, 5);
    TH2D *hSSDOccupancyLayer6 = new TH2D("hSSDOccupancyLayer6",";N_{modules};N_{Ladders}",fgkSSDMODULESPERLADDERLAYER6,0,fgkSSDMODULESPERLADDERLAYER6,3*fgkSSDLADDERSLAYER6,0,fgkSSDLADDERSLAYER6);
    for(Int_t iBin = 1; iBin < fgkSSDMODULESPERLADDERLAYER6 + 1; iBin++){
      sprintf(fLabel,"%d",iBin);
      hSSDOccupancyLayer6->GetXaxis()->SetBinLabel(iBin,fLabel);
    }
    fAliITSQADataMakerRec->Add2RawsList(hSSDOccupancyLayer6, 6);
    fRawsOffset = 7;
    
    TH1D *hOccupancyLadder[2*(fgkSSDLADDERSLAYER5 + fgkSSDLADDERSLAYER6)];
    for(Int_t iLayer = 5; iLayer < 7; iLayer++) {
      for(Int_t iLadder = 1; iLadder < AliITSgeomTGeo::GetNLadders(iLayer) + 1; iLadder++) {
	//P-side occupancy plots
	fTitle = "SSD_Occupancy_Layer"; fTitle += iLayer;
	fTitle += "_Ladder"; fTitle += iLadder; fTitle += "_PSide";
	hOccupancyLadder[fSSDhRaws] = new TH1D(fTitle.Data(),fTitle.Data(),AliITSgeomTGeo::GetNDetectors(iLayer),0,AliITSgeomTGeo::GetNDetectors(iLayer));
	hOccupancyLadder[fSSDhRaws]->GetXaxis()->SetTitleColor(1);
	hOccupancyLadder[fSSDhRaws]->GetXaxis()->SetTitle("Module number");
	hOccupancyLadder[fSSDhRaws]->GetYaxis()->SetTitle("Occupancy [%]");
	fAliITSQADataMakerRec->Add2RawsList(hOccupancyLadder[fSSDhRaws], fRawsOffset + fSSDhRaws);	
	fSSDhRaws++;
	//N-side occupancy plots
	fTitle = "SSD_Occupancy_Layer"; fTitle += iLayer;
	fTitle += "_Ladder"; fTitle += iLadder; fTitle += "_NSide";
	hOccupancyLadder[fSSDhRaws] = new TH1D(fTitle.Data(),fTitle.Data(),AliITSgeomTGeo::GetNDetectors(iLayer),0,AliITSgeomTGeo::GetNDetectors(iLayer));
	hOccupancyLadder[fSSDhRaws]->GetXaxis()->SetTitleColor(1);
	hOccupancyLadder[fSSDhRaws]->GetXaxis()->SetTitle("Module number");
	hOccupancyLadder[fSSDhRaws]->GetYaxis()->SetTitle("Occupancy [%]");
	fAliITSQADataMakerRec->Add2RawsList(hOccupancyLadder[fSSDhRaws], fRawsOffset + fSSDhRaws);	
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
  
  Double_t fSizePerDDL[fgkNumOfDDLs];
  Double_t fSizePerLDC[20];
  Double_t sumSSDDataSize = 0.0;
  Double_t eventSize = -1.0;

  if(fkOnline) {
    //reset the signal vs strip number histograms
    for(Int_t i = 0; i < fgkSSDMODULES; i++) 
      fHistSSDRawSignalModule[i]->Reset();
  
    rawReader->Select("ITSSSD",-1,-1);  
    rawReader->Reset(); //rawReader->NextEvent();   
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
      fSizePerLDC[rawReader->GetLDCId()] = rawReader->GetDataSize();
      AliITSgeomTGeo::GetModuleId(fSSDStream.GetModuleID(),fLayer,fLadder,fModule);
      fStripNumber = (fSSDStream.GetSideFlag() == 0) ? fSSDStream.GetStrip() : fSSDStream.GetStrip() + fgkNumberOfPSideStrips;
      fHistPosition = (fLayer == 5) ? ((fLadder - 1)*fgkSSDMODULESPERLADDERLAYER5 + fModule - 1) : ((fLadder - 1)*fgkSSDMODULESPERLADDERLAYER6 + fModule + fgkSSDMODULESLAYER5 - 1);
      //AliInfo(Form("ModulePosition: %d - Layer: %d - Ladder: %d - Module: %d\n",fHistPosition,fLayer,fLadder,fModule));
      fHistSSDRawSignalModule[fHistPosition]->Fill(fStripNumber,fSSDStream.GetSignal());
    }//streamer loop   

    
    for(Int_t i = 0; i < fgkNumOfDDLs; i++) {
      sumSSDDataSize += fSizePerDDL[i];
      if(fSizePerDDL[i] > 0) 
	(fAliITSQADataMakerRec->GetRawsData(2))->Fill(i+512);
      (fAliITSQADataMakerRec->GetRawsData(3))->Fill(i+512,fSizePerDDL[i]/1e+06);
    }
    for(Int_t i = 0; i < 20; i++) 
      (fAliITSQADataMakerRec->GetRawsData(4))->Fill(i,fSizePerLDC[i]/1e+06);
    if(sumSSDDataSize) 
      (fAliITSQADataMakerRec->GetRawsData(0))->Fill(TMath::Log10(sumSSDDataSize));
    if(eventSize)
      (fAliITSQADataMakerRec->GetRawsData(1))->Fill(100.*sumSSDDataSize/eventSize);
 
    //AliInfo(Form("EVENT: %d\n",fSSDEvent));
    Double_t occupancy = 0.0;
    Int_t lLadderLocationY = 0;
    
    Int_t fOccPosition = 0;
    for(Int_t i = 0; i < fgkSSDMODULES; i++) {
      AliITSgeomTGeo::GetModuleId(i+500,fLayer,fLadder,fModule); 
      fOccPosition  = (fLayer == 5) ? fLadder : fLadder + fgkSSDLADDERSLAYER5;
      //P-SIDE OCCUPANCY
      occupancy = GetSSDOccupancyRaws(fHistSSDRawSignalModule[i],0);
      fAliITSQADataMakerRec->GetRawsData(fRawsOffset -1 + 2*fLadder - 1)->Fill(fModule,occupancy);
      lLadderLocationY = 3*fLadder; // sideP=1 sideN=0
      ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fLayer))->SetBinContent(fModule,lLadderLocationY,occupancy);
      //N-SIDE OCCUPANCY
      occupancy = GetSSDOccupancyRaws(fHistSSDRawSignalModule[i],1);
      fAliITSQADataMakerRec->GetRawsData(fRawsOffset -1 + 2*fLadder)->Fill(fModule,occupancy);
      ((TH2D *)fAliITSQADataMakerRec->GetRawsData(fLayer))->SetBinContent(fModule,lLadderLocationY-1,occupancy);
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
    if (lHisto->GetBinContent(iBin))
      lNumFiredBins++; 
  }
  
  Double_t lOccupancy = (100.*lNumFiredBins)/lHisto->GetNbinsX(); // percentage
  //AliInfo(Form("Fired strips: %d - Total strips: %d - Occupancy :%lf\n",lNumFiredBins,lHisto->GetNbinsX(),lOccupancy));
  
  return lOccupancy;
}


//____________________________________________________________________________ 
void AliITSQASSDDataMakerRec::InitRecPoints()
{
  // Initialization for RECPOINTS - SSD -
  fRecsOffset = (fAliITSQADataMakerRec->fRecPointsQAList)->GetEntries();

// custom code here
//fSSDhRecs must be incremented by one unit every time a histogram is ADDED to the QA List

  AliDebug(1,Form("%d SSD Recs histograms booked\n",fSSDhRecs));
}

//____________________________________________________________________________ 
void AliITSQASSDDataMakerRec::MakeRecPoints(TTree * /*clustersTree*/)
{
  // Fill QA for recpoints - SSD -
}

