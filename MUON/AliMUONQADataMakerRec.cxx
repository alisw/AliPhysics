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

// $Id$

// --- MUON header files ---
#include "AliMUONQADataMakerRec.h"

#include "AliMUONCluster.h"  
#include "AliMUONRawStreamTracker.h"
#include "AliMUONRawStreamTrigger.h"
#include "AliMUONConstants.h"  
#include "AliMUONVClusterStore.h"
#include "AliMUONVCluster.h"

#include "AliMUONDigitMaker.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVTriggerStore.h"
#include "AliMUONVDigit.h"
#include "AliMUONLocalTrigger.h"

#include "AliMUONDDLTrigger.h"
#include "AliMUONRegHeader.h"
#include "AliMUONDarcHeader.h"
#include "AliMUONLocalStruct.h"

#include "AliMUONGeometryTransformer.h"

#include "AliMpDDLStore.h"
#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"
#include "AliMpConstants.h"
#include "AliMpPad.h"
#include "AliMpBusPatch.h"
#include "AliMpTriggerCrate.h"
#include "AliMpLocalBoard.h"
#include "AliMpCDB.h"

// --- AliRoot header files ---
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDMuonCluster.h"
#include "AliLog.h"
#include "AliRawReader.h"
#include "AliQAChecker.h"

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH1I.h> 
#include <TH2F.h>
#include <TH3F.h> 
#include <TLorentzVector.h>
#include <Riostream.h>

//-----------------------------------------------------------------------------
/// \class AliMUONQADataMakerRec
///
/// MUON base class for quality assurance data (histo) maker
///
/// \author C. Finck

/// \cond CLASSIMP
ClassImp(AliMUONQADataMakerRec)
/// \endcond
           
//____________________________________________________________________________ 
AliMUONQADataMakerRec::AliMUONQADataMakerRec() : 
    AliQADataMakerRec(AliQA::GetDetName(AliQA::kMUON), "MUON Quality Assurance Data Maker"),
    fIsInitRaws(kFALSE),
    fIsInitRecPoints(kFALSE),
    fIsInitESDs(kFALSE),
    fDigitStore(0x0),
    fTriggerStore(0x0),
    fDigitMaker(0x0)
{
    /// ctor
  fDigitStore = AliMUONVDigitStore::Create("AliMUONDigitStoreV1");
  fDigitMaker = new AliMUONDigitMaker(kTRUE);
}

//____________________________________________________________________________ 
AliMUONQADataMakerRec::AliMUONQADataMakerRec(const AliMUONQADataMakerRec& qadm) :
    AliQADataMakerRec(qadm),
    fIsInitRaws(kFALSE),
    fIsInitRecPoints(kFALSE),
    fIsInitESDs(kFALSE),
    fDigitStore(0x0),
    fTriggerStore(0x0),
    fDigitMaker(0x0)
{
    ///copy ctor 
    SetName((const char*)qadm.GetName()) ; 
    SetTitle((const char*)qadm.GetTitle()); 

    // Do not copy the digit store and digit maker, but create its own ones
    fDigitStore = AliMUONVDigitStore::Create("AliMUONDigitStoreV1");
    fDigitMaker = new AliMUONDigitMaker(kTRUE);
}

//__________________________________________________________________
AliMUONQADataMakerRec& AliMUONQADataMakerRec::operator = (const AliMUONQADataMakerRec& qadm )
{
  /// Assignment operator

  // check assignment to self
  if (this == &qadm) return *this;

  this->~AliMUONQADataMakerRec();
  new(this) AliMUONQADataMakerRec(qadm);
  return *this;
}

//__________________________________________________________________
AliMUONQADataMakerRec::~AliMUONQADataMakerRec()
{
    /// dtor
    delete fDigitStore;
    delete fTriggerStore;
    delete fDigitMaker;
}

//____________________________________________________________________________ 
void AliMUONQADataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray* list)
{
    ///Detector specific actions at end of cycle

    // Display trigger histos in a more user friendly way
    DisplayTriggerInfo(task);

    // do the QA checking
    AliQAChecker::Instance()->Run(AliQA::kMUON, task, list) ;
}

//____________________________________________________________________________ 
void AliMUONQADataMakerRec::InitRaws()
{
    /// create Raws histograms in Raws subdir
    TH1I* h0 = new TH1I("hRawBusPatch", "buspatch distribution",  1932, 1, 1932); 
    Add2RawsList(h0, kRawBusPatch);

    TH1I* h1 = new TH1I("hRawCharge", "Charge distribution in rawdata", 4096, 0, 4095); 
    Add2RawsList(h1, kRawCharge);
		
    for (Int_t iDDL = 0; iDDL < 20; ++iDDL) 
    {
      TH1F* h2 = new TH1F(Form("%s%d", "hRawBusPatchDDL", iDDL), Form("%s %d","RAW Buspatch distribution for DDL", iDDL), 50, 0, 49); 
      Add2RawsList(h2, kRawBuspatchDDL+iDDL);
    }

    TH3F* h3 = new TH3F("hTriggerScalersBendPlane", "Trigger scalers in bending plane",
			4, 10.5, 14.5,
			234, 0.5, 234.5,
			16, -0.5, 15.5);
    h3->GetXaxis()->SetTitle("Chamber");
    h3->GetYaxis()->SetTitle("Board");
    h3->GetZaxis()->SetTitle("Strip");
    Add2RawsList(h3, kTriggerScalersBP);

    TH3F* h4 = new TH3F("hTriggerScalersNonBendPlane", "Trigger scalers in non-bending plane",
			4, 10.5, 14.5,
			234, 0.5, 234.5,
			16, -0.5, 15.5);
    h4->GetXaxis()->SetTitle("Chamber");
    h4->GetYaxis()->SetTitle("Board");
    h4->GetZaxis()->SetTitle("Strip");
    Add2RawsList(h4, kTriggerScalersNBP);

    fIsInitRaws = InitDisplayHistos(AliQA::kRAWS);
}

//____________________________________________________________________________ 
void AliMUONQADataMakerRec::InitRecPoints()
{
    /// create Reconstructed Points histograms in RecPoints subdir
    TH3F* h0 = new TH3F("hTriggerDigitsBendPlane", "Trigger digits in bending plane",
			4, 10.5, 14.5,
			234, 0.5, 234.5,
			16, -0.5, 15.5);
    h0->GetXaxis()->SetTitle("Chamber");
    h0->GetYaxis()->SetTitle("Board");
    h0->GetZaxis()->SetTitle("Strip");
    Add2RecPointsList(h0, kTriggerDigitsBendPlane);

    TH3F* h1 = new TH3F("hTriggerDigitsNonBendPlane", "Trigger digits in non-bending plane",
			4, 10.5, 14.5,
			234, 0.5, 234.5,
			16, -0.5, 15.5);
    h1->GetXaxis()->SetTitle("Chamber");
    h1->GetYaxis()->SetTitle("Board");
    h1->GetZaxis()->SetTitle("Strip");
    Add2RecPointsList(h1, kTriggerDigitsNonBendPlane);

    TH1F* h2 = new TH1F("hTriggeredBoards", "Triggered boards", 234, 0.5, 234.5);
    Add2RecPointsList(h2, kTriggeredBoards);

    fIsInitRecPoints = InitDisplayHistos(AliQA::kRECPOINTS);
}


//____________________________________________________________________________ 
void AliMUONQADataMakerRec::InitESDs()
{
    ///create ESDs histograms in ESDs subdir
  TH1F* h0 = new TH1F("hESDnTracks", "ESDs track number distribution", 30, 0., 30.);  
  Add2ESDsList(h0, kESDnTracks);

  TH1F* h1 = new TH1F("hESDMomentum", "ESDs P distribution", 300, 0., 300) ; 
  Add2ESDsList(h1, kESDMomentum);

  TH1F* h2 = new TH1F("hESDPt", "ESDs Pt distribution", 200, 0., 50) ; 
  Add2ESDsList(h2, kESDPt);

  TH1F* h3 = new TH1F("hESDRapidity", "ESDs rapidity distribution", 200, -4.5,-2.) ; 
  Add2ESDsList(h3, kESDRapidity);

  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); ++i) 
  {
    TH2F* h4 = new TH2F(Form("%s%d", "hESDClusterHitMap", i), 
		     Form("%s %d", "ESD Clusters hit distribution for chamber", i),
		     100, -1*AliMUONConstants::Rmax(i/2), AliMUONConstants::Rmax(i/2),
		     100, -1*AliMUONConstants::Rmax(i/2), AliMUONConstants::Rmax(i/2)); 
    Add2ESDsList(h4, kESDClusterHitMap+i);
  }
  
  fIsInitESDs =  kTRUE;
}

//____________________________________________________________________________
void AliMUONQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
    /// make QA for rawdata

    if ( ! fIsInitRaws ) {
      AliWarningStream() 
        << "Skipping function due to a failure in Init" << endl;
      return;
    }    

    Int_t busPatchId;
    UShort_t manuId;  
    UChar_t channelId;
    UShort_t charge;

    rawReader->Reset();
    AliMUONRawStreamTracker rawStreamTrack(rawReader);
    rawStreamTrack.First();
    while( rawStreamTrack.Next(busPatchId, manuId, channelId, charge) ) {
  
      GetRawsData(kRawBusPatch)->Fill(busPatchId);
      GetRawsData(kRawCharge)->Fill(charge);
      Int_t iDDL = rawStreamTrack.GetCurentDDL();
      GetRawsData(kRawBuspatchDDL + iDDL)->Fill(AliMpBusPatch::GetLocalBusID(busPatchId, iDDL));
		
		  
    } // Next digit


    // Get trigger scalers

    Int_t loCircuit=0;
    AliMpCDB::LoadDDLStore();

    AliMUONRawStreamTrigger rawStreamTrig(rawReader);
    while (rawStreamTrig.NextDDL()) 
    {
      // If not a scaler event, do nothing
      Bool_t scalerEvent =  rawReader->GetDataHeader()->GetL1TriggerMessage() & 0x1;
      if(!scalerEvent) break;

      AliMUONDDLTrigger* ddlTrigger = rawStreamTrig.GetDDLTrigger();
      AliMUONDarcHeader* darcHeader = ddlTrigger->GetDarcHeader();

      Int_t nReg = darcHeader->GetRegHeaderEntries();
    
      for(Int_t iReg = 0; iReg < nReg ;iReg++)
      {   //reg loop

	// crate info  
	AliMpTriggerCrate* crate = AliMpDDLStore::Instance()->
	  GetTriggerCrate(rawStreamTrig.GetDDL(), iReg);

	AliMUONRegHeader* regHeader =  darcHeader->GetRegHeaderEntry(iReg);

	// loop over local structures
	Int_t nLocal = regHeader->GetLocalEntries();
	for(Int_t iLocal = 0; iLocal < nLocal; iLocal++) 
	{
	  AliMUONLocalStruct* localStruct = regHeader->GetLocalEntry(iLocal);
        
	  // if card exist
	  if (!localStruct) continue;
          
	  loCircuit = crate->GetLocalBoardId(localStruct->GetId());

	  if ( !loCircuit ) continue; // empty slot

	  AliMpLocalBoard* localBoard = AliMpDDLStore::Instance()->GetLocalBoard(loCircuit, false);

	  // skip copy cards
	  if( !localBoard->IsNotified()) 
	    continue;

	  Int_t cathode = localStruct->GetComptXY()%2;

	  ERaw hindex = (cathode==0) ? kTriggerScalersBP : kTriggerScalersNBP;

	  // loop over strips
	  for (Int_t ibitxy = 0; ibitxy < 16; ++ibitxy) {
	    if(localStruct->GetXY1(ibitxy) > 0)
	      ((TH3F*)GetRawsData(hindex))->Fill(11+0, loCircuit, ibitxy, 2*localStruct->GetXY1(ibitxy));
	    if(localStruct->GetXY2(ibitxy) > 0)
	      ((TH3F*)GetRawsData(hindex))->Fill(11+1, loCircuit, ibitxy, 2*localStruct->GetXY2(ibitxy));
	    if(localStruct->GetXY3(ibitxy) > 0)
	      ((TH3F*)GetRawsData(hindex))->Fill(11+2, loCircuit, ibitxy, 2*localStruct->GetXY3(ibitxy));
	    if(localStruct->GetXY4(ibitxy) > 0)
	      ((TH3F*)GetRawsData(hindex))->Fill(11+3, loCircuit, ibitxy, 2*localStruct->GetXY4(ibitxy));
	  } // loop on strips
	} // iLocal
      } // iReg
    } // NextDDL
}

//____________________________________________________________________________
void AliMUONQADataMakerRec::MakeRecPoints(TTree* clustersTree)
{
  
    /// makes data from trigger response
      
    if ( ! fIsInitRecPoints ) {
      AliWarningStream() 
        << "Skipping function due to a failure in Init" << endl;
      return;
    }    

    // Fired pads info
    fDigitStore->Clear();

    if (!fTriggerStore) fTriggerStore = AliMUONVTriggerStore::Create(*clustersTree);
    fTriggerStore->Clear();
    fTriggerStore->Connect(*clustersTree, false);
    clustersTree->GetEvent(0);

    AliMUONLocalTrigger* locTrg;
    TIter nextLoc(fTriggerStore->CreateLocalIterator());

    while ( ( locTrg = static_cast<AliMUONLocalTrigger*>(nextLoc()) ) ) 
    {
      if (locTrg->IsNull()) continue;
   
      TArrayS xyPattern[2];
      locTrg->GetXPattern(xyPattern[0]);
      locTrg->GetYPattern(xyPattern[1]);

      Int_t nBoard = locTrg->LoCircuit();

      Bool_t xTrig=locTrg->IsTrigX();
      Bool_t yTrig=locTrg->IsTrigY();
    
      if (xTrig && yTrig)
	((TH1F*)GetRecPointsData(kTriggeredBoards))->Fill(nBoard);

      fDigitMaker->TriggerDigits(nBoard, xyPattern, *fDigitStore);
    }

    TIter nextDigit(fDigitStore->CreateIterator());
    AliMUONVDigit* mDigit;
    while ( ( mDigit = static_cast<AliMUONVDigit*>(nextDigit()) ) )
    {
      Int_t detElemId = mDigit->DetElemId();
      Int_t ch = detElemId/100;
      Int_t localBoard = mDigit->ManuId();
      Int_t channel = mDigit->ManuChannel();
      Int_t cathode = mDigit->Cathode();
      ERecPoints hindex 
        = ( cathode == 0 ) ? kTriggerDigitsBendPlane : kTriggerDigitsNonBendPlane;
      
      ((TH3F*)GetRecPointsData(hindex))->Fill(ch, localBoard, channel);
    }
}

//____________________________________________________________________________
void AliMUONQADataMakerRec::MakeESDs(AliESDEvent* esd)
{
    /// make QA data from ESDs

    if ( ! fIsInitESDs ) {
      AliWarningStream() 
        << "Skipping function due to a failure in Init" << endl;
      return;
    }    

    TLorentzVector v1;

    Int_t nTracks = (Int_t)esd->GetNumberOfMuonTracks() ; 
    GetESDsData(0)->Fill(nTracks);

    for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {

      AliESDMuonTrack* muonTrack = esd->GetMuonTrack(iTrack);
      muonTrack->LorentzP(v1);

      GetESDsData(1)->Fill(v1.P());
      GetESDsData(2)->Fill(v1.Pt());
      GetESDsData(3)->Fill(v1.Rapidity());

      TClonesArray clusters =  muonTrack->GetClusters();

      for (Int_t iCluster = 0; iCluster <clusters.GetEntriesFast(); ++iCluster) {
	AliESDMuonCluster* cluster = (AliESDMuonCluster*)clusters.At(iCluster);
	GetESDsData(kESDClusterHitMap+cluster->GetChamberId())
          ->Fill(cluster->GetX(), cluster->GetY());
      }
    }
}

//____________________________________________________________________________ 
void AliMUONQADataMakerRec::StartOfDetectorCycle()
{
    /// Detector specific actions at start of cycle
  
}

//____________________________________________________________________________ 
void AliMUONQADataMakerRec::DisplayTriggerInfo(AliQA::TASKINDEX_t task)
{
  //
  /// Display trigger information in a user-friendly way:
  /// from local board and strip numbers to their position on chambers
  //
  if(task!=AliQA::kRECPOINTS && task!=AliQA::kRAWS) return;

  AliMpCDB::LoadDDLStore();

  Int_t detElemId;
  Float_t xg1, yg1, zg1, xg2, yg2, zg2, binContent;
  Float_t xLocal1=0., yLocal1=0., xLocal2=0., yLocal2=0., xWidth=0., yWidth=0.;
  Float_t x1Label, y1Label, x2Label, y2Label;
  Int_t x1Int, x2Int, y1Int, y2Int;
  Int_t binCh, binBoard, binStrip;

  AliMUONGeometryTransformer transform;
  transform.LoadGeometryData();

  TH3F* histoStrips=0x0;
  TH2F* histoDisplayStrips=0x0;
  TH2F* histoDisplayBoards=0x0;
  TH1F* histoBoards=0x0;

  const Float_t kShift = 0.15;

  for (Int_t iCath = 0; iCath < 2; ++iCath)
  {

    if(task==AliQA::kRECPOINTS){
      ERecPoints hindex 
	= ( iCath == 0 ) ? kTriggerDigitsBendPlane : kTriggerDigitsNonBendPlane;
      histoStrips = (TH3F*)GetRecPointsData(hindex);
      histoBoards = (TH1F*)GetRecPointsData(kTriggeredBoards);
    }
    else if(task==AliQA::kRAWS){
      ERaw hindex 
	= ( iCath == 0 ) ? kTriggerScalersBP : kTriggerScalersNBP;
      histoStrips = (TH3F*)GetRawsData(hindex);
      if(histoStrips->GetEntries()==0) return; // No scalers found
    }
    
    // Load mapping
    if ( ! AliMpSegmentation::Instance(kFALSE) ) {
      /// Load mapping
      if ( ! AliMpCDB::LoadDDLStore() ) {
        AliError("Could not access mapping from OCDB !");
        return;
      }
    }  

    for (Int_t iChamber = 0; iChamber < 4; ++iChamber)
    {
      Int_t iCh = iChamber + AliMpConstants::NofTrackingChambers();

      if(task==AliQA::kRECPOINTS){
	histoDisplayStrips = (TH2F*)GetRecPointsData(kTriggerDigitsDisplay + 4*iCath + iChamber);
	histoDisplayBoards = (TH2F*)GetRecPointsData(kTriggerBoardsDisplay);
      }
      else if(task==AliQA::kRAWS){
	histoDisplayStrips = (TH2F*)GetRawsData(kTriggerScalersDisplay + 4*iCath + iChamber);
      }

      for(Int_t iBoard=1; iBoard<=234; iBoard++){
	// get detElemId
	detElemId = AliMpDDLStore::Instance()->GetDEfromLocalBoard(iBoard, iCh);

	if (!detElemId) continue;

	AliMpLocalBoard* localBoard = AliMpDDLStore::Instance()->GetLocalBoard(iBoard, false);

	// skip copy cards
	if( !localBoard->IsNotified()) 
	  continue;

	const AliMpVSegmentation* seg 
	  = AliMpSegmentation::Instance()
	  ->GetMpSegmentation(detElemId, AliMp::GetCathodType(iCath));  
        
	// loop over strips
	for (Int_t ibitxy = 0; ibitxy < 16; ++ibitxy) 
	{
	  AliMpPad pad = seg->PadByLocation(AliMpIntPair(iBoard,ibitxy),kFALSE);
                        
	  if (!pad.IsValid())
	    continue;

	  if(iCath==0){ // Geometry info from bending plane only
	    if(ibitxy==0) {
	      xLocal1 = pad.Position().X();
	      yLocal1 = pad.Position().Y();
	      xWidth = pad.Dimensions().X();
	      yWidth = pad.Dimensions().Y();
	    }
	    xLocal2 = pad.Position().X();
	    yLocal2 = pad.Position().Y();
	  }

	  Float_t dimX = pad.Dimensions().X();
	  Float_t dimY = pad.Dimensions().Y();

	  Float_t stripX = pad.Position().X();
	  Float_t stripY = pad.Position().Y();
		   
	  transform.Local2Global(detElemId, stripX, stripY, 0, xg1, yg1, zg1);

	  x1Int = histoDisplayStrips->GetXaxis()->FindBin(xg1 - dimX + kShift);
	  y1Int = histoDisplayStrips->GetYaxis()->FindBin(yg1 - dimY + kShift);
	  x2Int = histoDisplayStrips->GetXaxis()->FindBin(xg1 + dimX - kShift);
	  y2Int = histoDisplayStrips->GetYaxis()->FindBin(yg1 + dimY - kShift);

	  binCh = histoStrips->GetXaxis()->FindBin(iCh+1);
	  binBoard = histoStrips->GetYaxis()->FindBin(iBoard);
	  binStrip = histoStrips->GetZaxis()->FindBin(ibitxy);
	  binContent = histoStrips->GetBinContent(binCh, binBoard, binStrip);

	  if(binContent==0) continue;

	  for(Int_t binX=x1Int; binX<=x2Int; binX++){
	    for(Int_t binY=y1Int; binY<=y2Int; binY++){	      
	      histoDisplayStrips->SetBinContent(binX, binY, binContent);
	    }
	  }
	}// ibitxy

	if(task != AliQA::kRECPOINTS) continue;
	if(iCath>0 || iChamber>0) continue; // Geometry info from bending plane only
	transform.Local2Global(detElemId, xLocal1, yLocal1, 0, xg1, yg1, zg1);
	transform.Local2Global(detElemId, xLocal2, yLocal2, 0, xg2, yg2, zg2);

	// Per board
	x1Label = TMath::Min(xg1,xg2) - xWidth + kShift;
	y1Label = TMath::Min(yg1,yg2) - yWidth + kShift;
	x2Label = TMath::Max(xg1,xg2) + xWidth - kShift;
	y2Label = TMath::Max(yg1,yg2) + yWidth - kShift;

	x1Int = histoDisplayBoards->GetXaxis()->FindBin(x1Label);
	y1Int = histoDisplayBoards->GetYaxis()->FindBin(y1Label);
	x2Int = histoDisplayBoards->GetXaxis()->FindBin(x2Label);
	y2Int = histoDisplayBoards->GetYaxis()->FindBin(y2Label);

	binBoard = histoBoards->GetXaxis()->FindBin(iBoard);
	binContent = histoBoards->GetBinContent(binBoard);

	if(binContent==0) continue;

	for(Int_t binX=x1Int; binX<=x2Int; binX++){
	  for(Int_t binY=y1Int; binY<=y2Int; binY++){
	    histoDisplayBoards->SetBinContent(binX, binY, binContent);
	  }
	}
      } // iBoard
    } // iChamber
  } // iCath
}


//____________________________________________________________________________ 
Bool_t AliMUONQADataMakerRec::InitDisplayHistos(AliQA::TASKINDEX_t task)
{
  //
  /// Initialize trigger information display histograms,
  /// using mapping and geometry
  //
  AliMpCDB::LoadDDLStore();

  const Int_t kNumOfBoards = AliMpConstants::NofLocalBoards();

  AliMUONGeometryTransformer transform;
  transform.LoadGeometryData();

  TString cathCode[2] = {"BendPlane", "NonBendPlane"};

  TArrayF xLocal1[4], yLocal1[4], xLocal2[4], yLocal2[4], xWidth[4], yWidth[4];

  TArrayF xAxisStrip[8];
  TArrayF yAxisStrip[8];
  TArrayF xAxisBoard[8];
  TArrayF yAxisBoard[8];

  TH2F* histoStrips=0x0;
  TH2F* histoBoards=0x0;

  const Float_t kResetValue=1234567.;
   
  for(Int_t ch=0; ch<4; ch++){
    xLocal1[ch].Set(kNumOfBoards);
    yLocal1[ch].Set(kNumOfBoards);
    xLocal2[ch].Set(kNumOfBoards);
    yLocal2[ch].Set(kNumOfBoards);
    xWidth[ch].Set(kNumOfBoards);
    yWidth[ch].Set(kNumOfBoards);
  }

  for(Int_t cath=0; cath<2; cath++){
    for(Int_t ch=0; ch<4; ch++){
      Int_t chCath = 4*cath + ch;
      xAxisBoard[chCath].Set(60);
      xAxisBoard[chCath].Reset(kResetValue);
      yAxisBoard[chCath].Set(60);
      yAxisBoard[chCath].Reset(kResetValue);

      xAxisStrip[chCath].Set(700);
      xAxisStrip[chCath].Reset(kResetValue);
      yAxisStrip[chCath].Set(700);
      yAxisStrip[chCath].Reset(kResetValue);
    }
  }
   
  Float_t xg1, yg1, zg1, xg2, yg2, zg2;

  TString histoName, histoTitle;

  const Float_t kShift = 0.;

  // Load mapping
  if ( ! AliMpSegmentation::Instance(kFALSE) ) {
    /// Load mapping
    if ( ! AliMpCDB::LoadDDLStore() ) {
      AliError("Could not access mapping from OCDB !");
      return kFALSE;
    }
  }  

  for(Int_t iCath=0; iCath<2; iCath++){
    for (Int_t iChamber = 0; iChamber < 4; ++iChamber) {
      Int_t iCh = iChamber + AliMpConstants::NofTrackingChambers();
      for(Int_t iLoc = 0; iLoc < 234; iLoc++) {  
	Int_t iBoard = iLoc+1;
	Int_t detElemId = AliMpDDLStore::Instance()->GetDEfromLocalBoard(iBoard, iCh);

	if (!detElemId) continue;

	AliMpLocalBoard* localBoard = AliMpDDLStore::Instance()->GetLocalBoard(iBoard, kFALSE);

	// skip copy cards
	if( !localBoard->IsNotified()) 
	  continue;

	// get segmentation
	const AliMpVSegmentation* seg = AliMpSegmentation::Instance()
	  ->GetMpSegmentation(detElemId,
			      AliMp::GetCathodType(iCath));

	Int_t chCath = 4*iCath + iChamber;
	// loop over strips
	for (Int_t ibitxy = 0; ibitxy < 16; ++ibitxy) {
	  // get pad from electronics
	  AliMpPad pad = seg->PadByLocation(AliMpIntPair(iBoard,ibitxy),kFALSE);

	  if (!pad.IsValid()) continue;

	  if(iCath==0){
	    if(ibitxy==0) {
	      xLocal1[iChamber][iLoc] = pad.Position().X();
	      yLocal1[iChamber][iLoc] = pad.Position().Y();
	      xWidth[iChamber][iLoc] = pad.Dimensions().X();
	      yWidth[iChamber][iLoc] = pad.Dimensions().Y();
	    }
	    xLocal2[iChamber][iLoc] = pad.Position().X();
	    yLocal2[iChamber][iLoc] = pad.Position().Y();
	  }

	  // Check fired pads
	  Float_t dimX = pad.Dimensions().X();
	  Float_t dimY = pad.Dimensions().Y();
     
	  Float_t stripX = pad.Position().X();
	  Float_t stripY = pad.Position().Y();

	  transform.Local2Global(detElemId, stripX, stripY, 0, xg1, yg1, zg1);

	  AddSortedPoint(xg1 - dimX + kShift, xAxisStrip[chCath], kResetValue);
	  AddSortedPoint(xg1 + dimX - kShift, xAxisStrip[chCath], kResetValue);

	  AddSortedPoint(yg1 - dimY + kShift, yAxisStrip[chCath], kResetValue);
	  AddSortedPoint(yg1 + dimY - kShift, yAxisStrip[chCath], kResetValue);
	} // loop on strips  

	transform.Local2Global(detElemId, xLocal1[iChamber][iLoc], yLocal1[iChamber][iLoc], 0, xg1, yg1, zg1);
	transform.Local2Global(detElemId, xLocal2[iChamber][iLoc], yLocal2[iChamber][iLoc], 0, xg2, yg2, zg2);

	// Per board
	AddSortedPoint(TMath::Min(xg1,xg2) - xWidth[iChamber][iLoc] + kShift, xAxisBoard[chCath], kResetValue);
	AddSortedPoint(TMath::Max(xg1,xg2) + xWidth[iChamber][iLoc] - kShift, xAxisBoard[chCath], kResetValue);

	AddSortedPoint(TMath::Min(yg1,yg2) - yWidth[iChamber][iLoc] + kShift, yAxisBoard[chCath], kResetValue);
	AddSortedPoint(TMath::Max(yg1,yg2) + yWidth[iChamber][iLoc] - kShift, yAxisBoard[chCath], kResetValue);
      } // loop on local boards 
    } // loop on chambers
  } // loop on cathodes

  const Float_t kMinDiff = 0.1;

  // Book histos
  for(Int_t iCath=0; iCath<2; iCath++){
    for(Int_t iChamber=0; iChamber<4; iChamber++){
      Int_t chCath = 4*iCath + iChamber;
      Int_t ipoint=0;
      while(TMath::Abs(xAxisStrip[chCath][ipoint]-kResetValue)>kMinDiff) { ipoint++; }
      xAxisStrip[chCath].Set(ipoint);

      ipoint = 0;
      while(TMath::Abs(yAxisStrip[chCath][ipoint]-kResetValue)>kMinDiff) { ipoint++; }
      yAxisStrip[chCath].Set(ipoint);

      ipoint = 0;
      while(TMath::Abs(xAxisBoard[chCath][ipoint]-kResetValue)>kMinDiff) { ipoint++; }
      xAxisBoard[chCath].Set(ipoint);

      ipoint = 0;
      while(TMath::Abs(yAxisBoard[chCath][ipoint]-kResetValue)>kMinDiff) { ipoint++; }
      yAxisBoard[chCath].Set(ipoint);

      if(task==AliQA::kRECPOINTS){
	histoName = Form("hTriggerDigits%sChamber%i", cathCode[iCath].Data(), 11+iChamber);
	histoTitle = Form("Chamber %i: Fired pads %s", 11+iChamber, cathCode[iCath].Data());
      }
      else if(task==AliQA::kRAWS){
	histoName = Form("hScalers%sChamber%i", cathCode[iCath].Data(), 11+iChamber);
	histoTitle = Form("Chamber %i: Scalers %s", 11+iChamber, cathCode[iCath].Data());
      }
      histoStrips = new TH2F(histoName.Data(), histoTitle.Data(), 
			     xAxisStrip[chCath].GetSize()-1, xAxisStrip[chCath].GetArray(),
			     yAxisStrip[chCath].GetSize()-1, yAxisStrip[chCath].GetArray());
      histoStrips->SetXTitle("X (cm)");
      histoStrips->SetYTitle("Y (cm)");

      if(task==AliQA::kRECPOINTS){
	Add2RecPointsList(histoStrips, kTriggerDigitsDisplay + 4*iCath + iChamber);
      }
      else if(task==AliQA::kRAWS){
	Add2RawsList(histoStrips, kTriggerScalersDisplay + 4*iCath + iChamber);
      }

      if(task != AliQA::kRECPOINTS) continue;
      if(iCath>0 || iChamber>0) continue;

      histoName = "hFiredBoardsDisplay";
      histoTitle = "Fired boards";
      histoBoards = new TH2F(histoName.Data(), histoTitle.Data(), 
			     xAxisBoard[chCath].GetSize()-1, xAxisBoard[chCath].GetArray(),
			     yAxisBoard[chCath].GetSize()-1, yAxisBoard[chCath].GetArray());
      histoBoards->SetXTitle("X (cm)");
      histoBoards->SetYTitle("Y (cm)");

      Add2RecPointsList(histoBoards, kTriggerBoardsDisplay + 4*iCath + iChamber);
    } // loop on chamber
  } // loop on cathode
  
  return kTRUE;
}


//____________________________________________________________________________ 
Bool_t AliMUONQADataMakerRec::AddSortedPoint(Float_t currVal, TArrayF& position, const Float_t kResetValue)
{
  //
  /// Add sorted point in array according to an increasing order.
  /// Used to build display histograms axis.
  //
  Int_t nEntries = position.GetSize();
  Float_t tmp1, tmp2;
  const Float_t kMinDiff = 0.1;
  for(Int_t i=0; i<nEntries; i++){
    if(TMath::Abs(position[i]-currVal)<kMinDiff) return kFALSE;
    if(TMath::Abs(position[i]-kResetValue)<kMinDiff) {
      position[i] = currVal;
      return kTRUE;
    }
    if(currVal>position[i]) continue;
    tmp1 = position[i];
    position[i] = currVal;
    for(Int_t j=i+1; j<nEntries; j++){
      tmp2 = position[j];
      position[j] = tmp1;
      tmp1 = tmp2;
      if(tmp1==kResetValue) break;
    }
    return kTRUE;
  }
  return kFALSE;
}

