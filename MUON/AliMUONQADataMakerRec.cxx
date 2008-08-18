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

#include "AliMUON2DMap.h"
#include "AliMUONCluster.h"  
#include "AliMUONConstants.h"  
#include "AliMUONDDLTrigger.h"
#include "AliMUONDarcHeader.h"
#include "AliMUONDigitMaker.h"
#include "AliMUONLocalStruct.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONRawStreamTracker.h"
#include "AliMUONRawStreamTrigger.h"
#include "AliMUONRegHeader.h"
#include "AliMUONTrackerCalibratedDataMaker.h"
#include "AliMUONTriggerDisplay.h"
#include "AliMUONVCluster.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVTrackerData.h"
#include "AliMUONVTriggerStore.h"
#include "AliMpCDB.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpLocalBoard.h"
#include "AliMpStationType.h"
#include "AliMpTriggerCrate.h"
#include "AliRawEventHeaderBase.h"

// --- AliRoot header files ---
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDMuonCluster.h"
#include "AliLog.h"
#include "AliRawReader.h"
#include "AliQAChecker.h"
#include "AliCodeTimer.h"

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
/// \author C. Finck, D. Stocco, L. Aphecetche

/// \cond CLASSIMP
ClassImp(AliMUONQADataMakerRec)
/// \endcond
           
//____________________________________________________________________________ 
AliMUONQADataMakerRec::AliMUONQADataMakerRec() : 
AliQADataMakerRec(AliQA::GetDetName(AliQA::kMUON), "MUON Quality Assurance Data Maker"),
fIsInitRaws(kFALSE),
fIsInitRecPointsTracker(kFALSE),
fIsInitRecPointsTrigger(kFALSE),
fIsInitESDs(kFALSE),
fDigitStore(0x0),
fTriggerStore(0x0),
fDigitMaker(0x0),
fClusterStore(0x0),
fTrackerDataMaker(0x0)
{
    /// ctor
	
	Ctor();
}

//____________________________________________________________________________ 
void
AliMUONQADataMakerRec::Ctor()
{
	/// Init some members
	fDigitStore = AliMUONVDigitStore::Create("AliMUONDigitStoreV1");
	fDigitMaker = new AliMUONDigitMaker(kTRUE);
}

//____________________________________________________________________________ 
AliMUONQADataMakerRec::AliMUONQADataMakerRec(const AliMUONQADataMakerRec& qadm) :
AliQADataMakerRec(qadm),
fIsInitRaws(kFALSE),
fIsInitRecPointsTracker(kFALSE),
fIsInitRecPointsTrigger(kFALSE),
fIsInitESDs(kFALSE),
fDigitStore(0x0),
fTriggerStore(0x0),
fDigitMaker(0x0),
fClusterStore(0x0),
fTrackerDataMaker(0x0)
{
    ///copy ctor 
    SetName((const char*)qadm.GetName()) ; 
    SetTitle((const char*)qadm.GetTitle()); 

	// Do not copy the digit store and digit maker, but create its own ones
	
	Ctor();
	
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
  
  AliCodeTimerAuto("");
  
  delete fDigitStore;
  delete fTriggerStore;
  delete fDigitMaker;
	delete fClusterStore;
	delete fTrackerDataMaker;
}

//____________________________________________________________________________ 
void AliMUONQADataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray* list)
{
	///Detector specific actions at end of cycle
	
	// Display trigger histos in a more user friendly way
	DisplayTriggerInfo(task);
	
	if ( task == AliQA::kRAWS && fTrackerDataMaker ) 
	{
		TIter next(list);
		TObject* o;
		Bool_t alreadyThere(kFALSE);
		while ( ( o = next() ) && !alreadyThere )
		{
			TString classname(o->ClassName());
			if ( classname.Contains("TrackerData") ) alreadyThere = kTRUE;
		}
		if (!alreadyThere && fTrackerDataMaker) list->AddAt(fTrackerDataMaker->Data(),(Int_t)kTrackerData);
	}
	
	// do the QA checking
	AliQAChecker::Instance()->Run(AliQA::kMUON, task, list) ;
}

//____________________________________________________________________________ 
void AliMUONQADataMakerRec::InitRaws()
{
    /// create Raws histograms in Raws subdir
	
	if ( ! AliCDBManager::Instance()->GetDefaultStorage() )
	{
		AliError("CDB default storage not set. Cannot work.");
		fIsInitRaws=kFALSE;
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
	
	AliMUONTriggerDisplay triggerDisplay;
	TString histoName, histoTitle;
	for(Int_t iCath=0; iCath<AliMpConstants::NofCathodes(); iCath++){
		TString cathName = ( iCath==0 ) ? "BendPlane" : "NonBendPlane";
		for(Int_t iChamber=0; iChamber<AliMpConstants::NofTriggerChambers(); iChamber++){
			histoName = Form("hScalers%sChamber%i", cathName.Data(), 11+iChamber);
			histoTitle = Form("Chamber %i: Scalers %s", 11+iChamber, cathName.Data());
			TH2F* h5 = (TH2F*)triggerDisplay.GetEmptyDisplayHisto(histoName, AliMUONTriggerDisplay::kDisplayStrips, 
																														iCath, iChamber, histoTitle);
			Add2RawsList(h5, kTriggerScalersDisplay + AliMpConstants::NofTriggerChambers()*iCath + iChamber);
		}
	}
	
	fIsInitRaws = kTRUE;
}

//____________________________________________________________________________ 
void AliMUONQADataMakerRec::InitRecPoints()
{
	/// create Reconstructed Points histograms in RecPoints subdir
	InitRecPointsTrigger();
	InitRecPointsTracker();
}

//____________________________________________________________________________ 
void AliMUONQADataMakerRec::InitRecPointsTracker()
{
	/// create Reconstructed Points histograms in RecPoints subdir for the
	/// MUON tracker subsystem.

	AliMpDEIterator it;
	
	it.First();
	
	Int_t ndes(0);
	
	while ( !it.IsDone())
	{
		Int_t detElemId = it.CurrentDEId();
		
		it.Next();

		if ( AliMpDEManager::GetStationType(detElemId) != AliMp::kStationTrigger )
		{
			ndes = TMath::Max(ndes,detElemId);

			TH1* h = new TH1I(Form("hTrackerClusterMultiplicityForDE%04d",detElemId),
												Form("Multiplicity of the clusters in detection element %d",detElemId),
												100,0,100);
			
			h->GetXaxis()->SetTitle("Detection Element Id");
			
			Add2RecPointsList(h,kTrackerClusterMultiplicityPerDE+detElemId);
			
			h =  new TH1I(Form("hTrackerClusterChargeForDE%04d",detElemId),
										Form("Charge of the clusters in detection element %d",detElemId),
										100,0,1000);

			h->GetXaxis()->SetTitle("Detection Element Id");

			Add2RecPointsList(h,kTrackerClusterChargePerDE+detElemId);

		}

	}

	TH1* h = new TH1I("hTrackerNumberOfClustersPerDE","Number of clusters per detection element",
										ndes, -0.5, ndes - 0.5);

	h->GetXaxis()->SetTitle("Detection Element Id");

	Add2RecPointsList(h, kTrackerNumberOfClustersPerDE);

	for ( Int_t i = 0; i < AliMpConstants::NofTrackingChambers(); ++i ) 
	{
		TH1* h1 = new TH1I("hTrackerNumberOfClustersPerChamber","Number of clusters per chamber",AliMpConstants::NofTrackingChambers(),-0.5,AliMpConstants::NofTrackingChambers()-0.5);
		Add2RecPointsList(h1,kTrackerNumberOfClustersPerChamber);
		h1 = new TH1I(Form("hTrackerClusterMultiplicityForChamber%d",i),
								 Form("Cluster multiplicity for chamber %d",i),
								 100,0,100);
		Add2RecPointsList(h1,kTrackerClusterMultiplicityPerChamber+i);
		h1 = new TH1I(Form("hTrackerClusterChargeForChamber%d",i),
								 Form("Cluster charge for chamber %d",i),
								 100,0,1000);
		Add2RecPointsList(h1,kTrackerClusterChargePerChamber+i);
	}
	
	fIsInitRecPointsTracker=kTRUE;
}

//____________________________________________________________________________ 
void AliMUONQADataMakerRec::InitRecPointsTrigger()
{
	/// create Reconstructed Points histograms in RecPoints subdir for the
	/// MUON Trigger subsystem.
	
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

    AliMUONTriggerDisplay triggerDisplay;
    TString histoName, histoTitle;
    for(Int_t iCath=0; iCath<AliMpConstants::NofCathodes(); iCath++){
      TString cathName = ( iCath==0 ) ? "BendPlane" : "NonBendPlane";
      for(Int_t iChamber=0; iChamber<AliMpConstants::NofTriggerChambers(); iChamber++){
	histoName = Form("hTriggerDigits%sChamber%i", cathName.Data(), 11+iChamber);
	histoTitle = Form("Chamber %i: Fired pads %s", 11+iChamber, cathName.Data());
	TH2F* h3 = (TH2F*)triggerDisplay.GetEmptyDisplayHisto(histoName, AliMUONTriggerDisplay::kDisplayStrips, 
							      iCath, iChamber, histoTitle);
	Add2RecPointsList(h3, kTriggerDigitsDisplay + AliMpConstants::NofTriggerChambers()*iCath + iChamber);
      }
    }

    TH2F* h4 = (TH2F*)triggerDisplay.GetEmptyDisplayHisto("hFiredBoardsDisplay", AliMUONTriggerDisplay::kDisplayBoards,
							  0, 0, "Fired boards");
    Add2RecPointsList(h4, kTriggerBoardsDisplay);
	
	fIsInitRecPointsTrigger = kTRUE;
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

  if ( rawReader->GetType() == AliRawEventHeaderBase::kPhysicsEvent ) 
  {
    rawReader->Reset();
    MakeRawsTracker(rawReader);
  }
  
  rawReader->Reset();    
  MakeRawsTrigger(rawReader);
}

//____________________________________________________________________________
void AliMUONQADataMakerRec::MakeRawsTracker(AliRawReader* rawReader)
{
	/// make QA for rawdata tracker
  
	if (!fTrackerDataMaker) 
	{
		const Bool_t histogram(kFALSE);
		const Bool_t fastDecoder(kTRUE);
    
//    fTrackerDataMaker = new AliMUONTrackerRawDataMaker(rawReader,histogram,fastDecoder,takeRawReaderOwnership);

		fTrackerDataMaker = new AliMUONTrackerCalibratedDataMaker(GetRecoParam(),
                                                              AliCDBManager::Instance()->GetRun(),
                                                              rawReader,
																															AliCDBManager::Instance()->GetDefaultStorage()->GetURI(),
																															"NOGAIN",
																															histogram,
																															0.0,0.0,
                                                              fastDecoder);
		
		fTrackerDataMaker->Data()->DisableChannelLevel(); // to save up disk space, we only store starting at the manu level
		
		fTrackerDataMaker->SetRunning(kTRUE);
	}
	
	((AliMUONTrackerCalibratedDataMaker*)fTrackerDataMaker)->SetRawReader(rawReader);
	
	fTrackerDataMaker->ProcessEvent();
}

//____________________________________________________________________________
void AliMUONQADataMakerRec::MakeRawsTrigger(AliRawReader* rawReader)
{
	/// make QA for rawdata trigger
	
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
	/// Fill histograms from treeR
	
	if (fIsInitRecPointsTracker) MakeRecPointsTracker(clustersTree);
	if (fIsInitRecPointsTrigger) MakeRecPointsTrigger(clustersTree);
}

//____________________________________________________________________________
void AliMUONQADataMakerRec::MakeRecPointsTracker(TTree* clustersTree)
{
	/// Fill histograms related to tracker clusters 
	
	// First things first : do we have clusters in the TreeR ?
	// In "normal" production mode, it should be perfectly normal
	// *not* to have them.
	// But if for some reason we de-activated the combined tracking,
	// then we have clusters in TreeR, so let's take that opportunity
	// to QA them...
	
	if (!fClusterStore)
	{
		AliCodeTimerAuto("ClusterStore creation");
		fClusterStore = AliMUONVClusterStore::Create(*clustersTree);
		if (!fClusterStore) 
		{
			fIsInitRecPointsTracker = kFALSE;
			return;
		}
	}
	
	AliCodeTimerAuto("");
	
	fClusterStore->Connect(*clustersTree,kFALSE);
	clustersTree->GetEvent(0);

	TIter next(fClusterStore->CreateIterator());
	AliMUONVCluster* cluster;
	
	while ( ( cluster = static_cast<AliMUONVCluster*>(next()) ) )
	{
		Int_t detElemId = cluster->GetDetElemId();
		Int_t chamberId = AliMpDEManager::GetChamberId(detElemId);
		
		GetRecPointsData(kTrackerNumberOfClustersPerDE)->Fill(detElemId);
		GetRecPointsData(kTrackerClusterChargePerDE+detElemId)->Fill(cluster->GetCharge());
		GetRecPointsData(kTrackerClusterMultiplicityPerDE+detElemId)->Fill(cluster->GetNDigits());

		GetRecPointsData(kTrackerNumberOfClustersPerChamber)->Fill(chamberId);
		GetRecPointsData(kTrackerClusterChargePerChamber+chamberId)->Fill(cluster->GetCharge());
		GetRecPointsData(kTrackerClusterMultiplicityPerChamber+chamberId)->Fill(cluster->GetNDigits());
		
	}
	
	fClusterStore->Clear();
}

//____________________________________________________________________________
void AliMUONQADataMakerRec::MakeRecPointsTrigger(TTree* clustersTree)
{
	/// makes data from trigger response
      
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
      
      // skip "ghosts"
      if (!muonTrack->ContainTrackerData()) continue;
      
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

  AliMUONTriggerDisplay triggerDisplay;
  
  TH3F* histoStrips=0x0;
  TH2F* histoDisplayStrips=0x0;

  for (Int_t iCath = 0; iCath < AliMpConstants::NofCathodes(); iCath++)
  {
    if(task==AliQA::kRECPOINTS){
      ERecPoints hindex 
	= ( iCath == 0 ) ? kTriggerDigitsBendPlane : kTriggerDigitsNonBendPlane;
      histoStrips = (TH3F*)GetRecPointsData(hindex);
    }
    else if(task==AliQA::kRAWS){
      ERaw hindex 
	= ( iCath == 0 ) ? kTriggerScalersBP : kTriggerScalersNBP;
      histoStrips = (TH3F*)GetRawsData(hindex);
      if(histoStrips->GetEntries()==0) return; // No scalers found
    }
    
    for (Int_t iChamber = 0; iChamber < AliMpConstants::NofTriggerChambers(); iChamber++)
    {
      if(task==AliQA::kRECPOINTS){
	histoDisplayStrips = (TH2F*)GetRecPointsData(kTriggerDigitsDisplay + AliMpConstants::NofTriggerChambers()*iCath + iChamber);
      }
      else if(task==AliQA::kRAWS){
	histoDisplayStrips = (TH2F*)GetRawsData(kTriggerScalersDisplay + AliMpConstants::NofTriggerChambers()*iCath + iChamber);
      }
      Int_t bin = histoStrips->GetXaxis()->FindBin(11+iChamber);
      histoStrips->GetXaxis()->SetRange(bin,bin);
      TH2F* inputHisto = (TH2F*)histoStrips->Project3D("zy");
      triggerDisplay.FillDisplayHistogram(inputHisto, histoDisplayStrips, AliMUONTriggerDisplay::kDisplayStrips, iCath, iChamber);
    } // iChamber
  } // iCath

  if(task!=AliQA::kRECPOINTS) return;
  TH1F* histoBoards = (TH1F*)GetRecPointsData(kTriggeredBoards);
  TH2F* histoDisplayBoards = (TH2F*)GetRecPointsData(kTriggerBoardsDisplay);
  triggerDisplay.FillDisplayHistogram(histoBoards, histoDisplayBoards, AliMUONTriggerDisplay::kDisplayBoards, 0, 0);
}
