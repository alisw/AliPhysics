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
#include <TString.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliITSQASSDDataMakerRec.h"
#include "AliLog.h"
#include "AliQA.h"
#include "AliQAChecker.h"
#include "AliRawReader.h"
#include "AliITSRawStreamSSD.h"
#include "AliITSgeomTGeo.h"


ClassImp(AliITSQASSDDataMakerRec)

//____________________________________________________________________________ 
AliITSQASSDDataMakerRec::AliITSQASSDDataMakerRec(AliITSQADataMakerRec *aliITSQADataMakerRec, Bool_t kMode, Int_t ldc) :
TObject(),
fAliITSQADataMakerRec(aliITSQADataMakerRec),
fkOnline(kMode),
fLDC(ldc),
fSSDhRaws(0),
fSSDhRecs(0),
fRawsOffset(0),
fRecsOffset(0)
{
  //ctor used to discriminate OnLine-Offline analysis
  //initilize the raw signal vs strip number histograms
  Int_t fLayer = 0,fLadder = 0, fModule = 0;
  TString fTitle; Int_t fHistCounter = 0;
  for(Int_t i = 500; i < fgkSSDMODULES + 500; i++) {
    AliITSgeomTGeo::GetModuleId(i,fLayer,fLadder,fModule);
    fTitle = "SSD_RawSignal_Layer"; fTitle += fLayer;
    fTitle += "_Ladder"; fTitle += fLadder;
    fTitle += "_Module"; fTitle += fModule;
    fHistSSDRawSignalModule[fHistCounter] = new TH1D(fTitle.Data(),fTitle.Data(),1540,0,1540);
    fHistCounter++;
  }
}

//____________________________________________________________________________ 
AliITSQASSDDataMakerRec::AliITSQASSDDataMakerRec(const AliITSQASSDDataMakerRec& qadm) :
TObject(),
fAliITSQADataMakerRec(qadm.fAliITSQADataMakerRec),
fkOnline(qadm.fkOnline),
fLDC(qadm.fLDC),
fSSDhRaws(qadm.fSSDhRaws),
fSSDhRecs(qadm.fSSDhRecs),
fRawsOffset(qadm.fRawsOffset),
fRecsOffset(qadm.fRecsOffset)
{
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
    Int_t fLayer = 0, fLadder = 0, fModule = 0;
    Int_t fOccPosition = 0;
    for(Int_t i = 0; i < fgkSSDMODULES; i++) {
      AliITSgeomTGeo::GetModuleId(i+500,fLayer,fLadder,fModule); 
      fOccPosition  = (fLayer == 5) ? fLadder : fLadder + 34;
      Double_t occupancy = GetSSDOccupancyRaws(fHistSSDRawSignalModule[i]);
      //if(occupancy > 0) {
      //AliInfo(Form("========================================="));
      //AliInfo(Form("%s Position: %d - Initial histos: %d - Occupancy: %lf\n",fHistSSDRawSignalModule[i]->GetName(),fOccPosition,fInitialHistograms,occupancy));
      //}
      fAliITSQADataMakerRec->GetRawsData(fRawsOffset + fOccPosition - 1)->Fill(fModule,occupancy);
    }
  }//online flag for SSD
  //AliQAChecker::Instance()->Run( AliQA::kITS , task, list);
}

//____________________________________________________________________________ 
void AliITSQASSDDataMakerRec::InitRaws() {  

  // Initialization for RAW data - SSD -
  fRawsOffset = (fAliITSQADataMakerRec->fRawsQAList)->GetEntries();

  const Int_t fSSDLADDERS = 72;
  if(fkOnline) {
    AliInfo("Book Online Histograms for SSD\n");
  }
  else {
    AliInfo("Book Offline Histograms for SSD\n ");
  }
  AliInfo(Form("Number of histograms (SPD+SDD): %d\n",fRawsOffset));
  TString fTitle = 0;
  //Int_t fLayer = 0,fLadder = 0, fModule = 0;
  //book online QA histos
  if(fkOnline) {
    TH1D *hOccupancyLadder[fSSDLADDERS];
    for(Int_t iLayer = 5; iLayer < 7; iLayer++) {
      for(Int_t iLadder = 1; iLadder < AliITSgeomTGeo::GetNLadders(iLayer) + 1; iLadder++) {
	fTitle = "SSD_Occupancy_Layer"; fTitle += iLayer;
	fTitle += "_Ladder"; fTitle += iLadder;
	hOccupancyLadder[fSSDhRaws] = new TH1D(fTitle.Data(),fTitle.Data(),AliITSgeomTGeo::GetNDetectors(iLayer),0,AliITSgeomTGeo::GetNDetectors(iLayer));
	hOccupancyLadder[fSSDhRaws]->GetXaxis()->SetTitleColor(1);
	hOccupancyLadder[fSSDhRaws]->GetXaxis()->SetTitle("Module number");
	hOccupancyLadder[fSSDhRaws]->GetYaxis()->SetTitle("Occupancy [%]");
	hOccupancyLadder[fSSDhRaws]->SetMarkerStyle(kFullCircle);
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

  rawReader->Select("ITSSSD",-1,-1);  
  rawReader->Reset();                         
  AliITSRawStreamSSD fSSDStream(rawReader);    
  while (fSSDStream.Next()) {
    if(fSSDStream.GetModuleID() < 0) continue;
    AliITSgeomTGeo::GetModuleId(fSSDStream.GetModuleID(),fLayer,fLadder,fModule);
    fStripNumber  = (fSSDStream.GetSideFlag() == 0) ? fSSDStream.GetStrip() : fSSDStream.GetStrip() + 768;
    fHistPosition = (fLayer == 5) ? ((fLadder - 1)*22 + fModule - 1) : ((fLadder - 1)*25 + fModule + 748 - 1);
    //AliInfo(Form("ModulePosition: %d - Layer: %d - Ladder: %d - Module: %d\n",fHistPosition,fLayer,fLadder,fModule));
    fHistSSDRawSignalModule[fHistPosition]->Fill(fStripNumber,fSSDStream.GetSignal());
  }//streamer loop
}

//____________________________________________________________________________ 
Double_t AliITSQASSDDataMakerRec::GetSSDOccupancyRaws(TH1 *lHisto) { 
  // bo: TDC >0 or # of sigmas wrt noise ?
  Int_t lNumFiredBins = 0;
  for(Int_t iBin = 1; iBin < lHisto->GetNbinsX(); iBin++){
    if (lHisto->GetBinContent(iBin))
      lNumFiredBins++; 
  }
  Double_t lOccupancy = (100.*lNumFiredBins)/lHisto->GetNbinsX(); // percentage
  //if(lOccupancy > 0) 
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


