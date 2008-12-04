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
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONESDInterface.h"
#include "AliMpBusPatch.h"
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
	
  AliDebug(1,"");

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

    AliDebug(1,"");


    SetName((const char*)qadm.GetName()) ; 
    SetTitle((const char*)qadm.GetTitle()); 

	// Do not copy the digit store and digit maker, but create its own ones
	
	Ctor();
	
}

//__________________________________________________________________
AliMUONQADataMakerRec& AliMUONQADataMakerRec::operator = (const AliMUONQADataMakerRec& qadm )
{
  /// Assignment operator

  AliDebug(1,"");

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
  
  AliDebug(1,"");

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
  
  AliCodeTimerAuto("");
  
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
    if (!alreadyThere && fTrackerDataMaker) 
    {
      AliInfo("Adding fTrackerDataMaker to the list of qa objects");
      list->AddAt(fTrackerDataMaker->Data(),(Int_t)kTrackerData);
    }
    if ( fTrackerDataMaker ) 
    {
      TH1* hbp = GetRawsData(kTrackerBusPatchOccupancy);
      hbp->Reset();
      TIter nextBP(AliMpDDLStore::Instance()->CreateBusPatchIterator());
      AliMpBusPatch* bp(0x0);
      AliMUONVTrackerData* data = fTrackerDataMaker->Data();
      Int_t occDim = 2;
      
      while ( ( bp = static_cast<AliMpBusPatch*>(nextBP())) )
      {
        Int_t busPatchId = bp->GetId();
        Int_t bin = hbp->FindBin(busPatchId);
        hbp->SetBinContent(bin,data->BusPatch(busPatchId,occDim));
      }
    }
  }
  
  if ( task == AliQA::kESDS ) {
  // Normalize ESD histos
    TH1* h;
    Int_t bin;
    AliMpDEIterator it;
    it.First();
    while ( !it.IsDone()) {
    
      Int_t detElemId = it.CurrentDEId();
    
      if ( AliMpDEManager::GetStationType(detElemId) != AliMp::kStationTrigger ) {
      
        h = GetESDsData(kESDnClustersPerDE);
        Double_t nClusters = h->GetBinContent(h->GetXaxis()->FindFixBin((Double_t)detElemId));
      
        if (nClusters > 0) {
	
	  h = GetESDsData(kESDClusterChargePerDE);
	  bin = h->GetXaxis()->FindFixBin((Double_t)detElemId);
          h->SetBinContent(bin, h->GetBinContent(bin)/nClusters);
	
	  h = GetESDsData(kESDClusterMultPerDE);
	  bin = h->GetXaxis()->FindFixBin((Double_t)detElemId);
          h->SetBinContent(bin, h->GetBinContent(bin)/nClusters);
	
	  h = GetESDsData(kESDResidualXPerDEMean);
	  bin = h->GetXaxis()->FindFixBin((Double_t)detElemId);
	  Double_t meanResX = h->GetBinContent(bin)/nClusters;
          h->SetBinContent(bin, meanResX);
	
	  h = GetESDsData(kESDResidualYPerDEMean);
	  bin = h->GetXaxis()->FindFixBin((Double_t)detElemId);
	  Double_t meanResY = h->GetBinContent(bin)/nClusters;
          h->SetBinContent(bin, meanResY);
	
	  h = GetESDsData(kESDResidualXPerDESigma);
	  bin = h->GetXaxis()->FindFixBin((Double_t)detElemId);
	  if (nClusters > 1) h->SetBinContent(bin, TMath::Sqrt(h->GetBinContent(bin)/nClusters - meanResX*meanResX));
	  else h->SetBinContent(bin, 0.);
	
	  h = GetESDsData(kESDResidualYPerDESigma);
	  bin = h->GetXaxis()->FindFixBin((Double_t)detElemId);
	  if (nClusters > 1) h->SetBinContent(bin, TMath::Sqrt(h->GetBinContent(bin)/nClusters - meanResY*meanResY));
	  else h->SetBinContent(bin, 0.);
	
        }
      
      }
    
      it.Next();
    }
  
    Double_t nTracks = GetESDsData(kESDnClustersPerTrack)->GetEntries();
    if (nTracks > 0) {
      GetESDsData(kESDnClustersPerCh)->Scale(1./nTracks);
      GetESDsData(kESDnClustersPerDE)->Scale(1./nTracks);
    }
  }
  
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQA::kMUON, task, list) ;
}

//____________________________________________________________________________ 
void AliMUONQADataMakerRec::InitRaws()
{
    /// create Raws histograms in Raws subdir
	
  AliCodeTimerAuto("");
  
  Bool_t forExpert(kTRUE);
  
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
	Add2RawsList(h3, kTriggerScalersBP,forExpert);
	
	TH3F* h4 = new TH3F("hTriggerScalersNonBendPlane", "Trigger scalers in non-bending plane",
											4, 10.5, 14.5,
											234, 0.5, 234.5,
											16, -0.5, 15.5);
	h4->GetXaxis()->SetTitle("Chamber");
	h4->GetYaxis()->SetTitle("Board");
	h4->GetZaxis()->SetTitle("Strip");
	Add2RawsList(h4, kTriggerScalersNBP,forExpert);
	
	AliMUONTriggerDisplay triggerDisplay;
	TString histoName, histoTitle;
	for(Int_t iCath=0; iCath<AliMpConstants::NofCathodes(); iCath++){
		TString cathName = ( iCath==0 ) ? "BendPlane" : "NonBendPlane";
		for(Int_t iChamber=0; iChamber<AliMpConstants::NofTriggerChambers(); iChamber++){
			histoName = Form("hScalers%sChamber%i", cathName.Data(), 11+iChamber);
			histoTitle = Form("Chamber %i: Scalers %s", 11+iChamber, cathName.Data());
			TH2F* h5 = (TH2F*)triggerDisplay.GetEmptyDisplayHisto(histoName, AliMUONTriggerDisplay::kDisplayStrips, 
									      iCath, iChamber, histoTitle);
			Add2RawsList(h5, kTriggerScalersDisplay + AliMpConstants::NofTriggerChambers()*iCath + iChamber,forExpert);
		}
	}
	
  Int_t nbp(0);
  TIter next(AliMpDDLStore::Instance()->CreateBusPatchIterator());
  while (next())
  {
    ++nbp;
  }
  
  TH1* hbp = new TH1F("hTrackerBusPatchOccupancy","Occupancy of bus patches",
                      nbp,-0.5,nbp-0.5);
  
  Add2RawsList(hbp,kTrackerBusPatchOccupancy,!forExpert);
  
	fIsInitRaws = kTRUE;
}

//____________________________________________________________________________ 
void AliMUONQADataMakerRec::InitRecPoints()
{
	/// create Reconstructed Points histograms in RecPoints subdir
  
  AliCodeTimerAuto("");
  
	InitRecPointsTrigger();
	InitRecPointsTracker();
}

//____________________________________________________________________________ 
void AliMUONQADataMakerRec::InitRecPointsTracker()
{
	/// create Reconstructed Points histograms in RecPoints subdir for the
	/// MUON tracker subsystem.

  AliCodeTimerAuto("");
  
  Bool_t forExpert(kTRUE);
  
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
			
			Add2RecPointsList(h,kTrackerClusterMultiplicityPerDE+detElemId,forExpert);
			
			h =  new TH1I(Form("hTrackerClusterChargeForDE%04d",detElemId),
				      Form("Charge of the clusters in detection element %d",detElemId),
				      100,0,1000);

			h->GetXaxis()->SetTitle("Detection Element Id");

			Add2RecPointsList(h,kTrackerClusterChargePerDE+detElemId,forExpert);

		}

	}

	TH1* h = new TH1I("hTrackerNumberOfClustersPerDE","Number of clusters per detection element",
										ndes, 0.5, ndes + 0.5);

	h->GetXaxis()->SetTitle("Detection Element Id");

	Add2RecPointsList(h, kTrackerNumberOfClustersPerDE,!forExpert);

	for ( Int_t i = 0; i < AliMpConstants::NofTrackingChambers(); ++i ) 
	{
		TH1* h1 = new TH1I("hTrackerNumberOfClustersPerChamber","Number of clusters per chamber",
				   AliMpConstants::NofTrackingChambers(),-0.5,AliMpConstants::NofTrackingChambers()-0.5);
		Add2RecPointsList(h1,kTrackerNumberOfClustersPerChamber,forExpert);
		h1 = new TH1I(Form("hTrackerClusterMultiplicityForChamber%d",i),
			      Form("Cluster multiplicity for chamber %d",i),
			      100,0,100);
		Add2RecPointsList(h1,kTrackerClusterMultiplicityPerChamber+i,forExpert);
		h1 = new TH1I(Form("hTrackerClusterChargeForChamber%d",i),
			      Form("Cluster charge for chamber %d",i),
			      100,0,1000);
		Add2RecPointsList(h1,kTrackerClusterChargePerChamber+i,forExpert);
	}
	
	fIsInitRecPointsTracker=kTRUE;
}

//____________________________________________________________________________ 
void AliMUONQADataMakerRec::InitRecPointsTrigger()
{
	/// create Reconstructed Points histograms in RecPoints subdir for the
	/// MUON Trigger subsystem.
	
  Bool_t forExpert(kTRUE);
  
    TH3F* h0 = new TH3F("hTriggerDigitsBendPlane", "Trigger digits in bending plane",
			4, 10.5, 14.5,
			234, 0.5, 234.5,
			16, -0.5, 15.5);
    h0->GetXaxis()->SetTitle("Chamber");
    h0->GetYaxis()->SetTitle("Board");
    h0->GetZaxis()->SetTitle("Strip");
    Add2RecPointsList(h0, kTriggerDigitsBendPlane,forExpert);

    TH3F* h1 = new TH3F("hTriggerDigitsNonBendPlane", "Trigger digits in non-bending plane",
			4, 10.5, 14.5,
			234, 0.5, 234.5,
			16, -0.5, 15.5);
    h1->GetXaxis()->SetTitle("Chamber");
    h1->GetYaxis()->SetTitle("Board");
    h1->GetZaxis()->SetTitle("Strip");
    Add2RecPointsList(h1, kTriggerDigitsNonBendPlane,forExpert);

    TH1F* h2 = new TH1F("hTriggeredBoards", "Triggered boards", 234, 0.5, 234.5);
    Add2RecPointsList(h2, kTriggeredBoards,forExpert);

    AliMUONTriggerDisplay triggerDisplay;
    TString histoName, histoTitle;
    for(Int_t iCath=0; iCath<AliMpConstants::NofCathodes(); iCath++){
      TString cathName = ( iCath==0 ) ? "BendPlane" : "NonBendPlane";
      for(Int_t iChamber=0; iChamber<AliMpConstants::NofTriggerChambers(); iChamber++){
	histoName = Form("hTriggerDigits%sChamber%i", cathName.Data(), 11+iChamber);
	histoTitle = Form("Chamber %i: Fired pads %s", 11+iChamber, cathName.Data());
	TH2F* h3 = (TH2F*)triggerDisplay.GetEmptyDisplayHisto(histoName, AliMUONTriggerDisplay::kDisplayStrips, 
							      iCath, iChamber, histoTitle);
	Add2RecPointsList(h3, kTriggerDigitsDisplay + AliMpConstants::NofTriggerChambers()*iCath + iChamber,forExpert);
      }
    }

    TH2F* h4 = (TH2F*)triggerDisplay.GetEmptyDisplayHisto("hFiredBoardsDisplay", AliMUONTriggerDisplay::kDisplayBoards,
							  0, 0, "Fired boards");
    Add2RecPointsList(h4, kTriggerBoardsDisplay,forExpert);
	
	fIsInitRecPointsTrigger = kTRUE;
}


//____________________________________________________________________________ 
void AliMUONQADataMakerRec::InitESDs()
{
  ///create ESDs histograms in ESDs subdir
  
  Bool_t forExpert(kTRUE);
  
  Int_t nCh = AliMUONConstants::NTrackingCh();
  Int_t nDE = 1100;
  
  // track info
  TH1F* hESDnTracks = new TH1F("hESDnTracks", "number of tracks", 20, 0., 20.);
  Add2ESDsList(hESDnTracks, kESDnTracks,!forExpert);

  TH1F* hESDMatchTrig = new TH1F("hESDMatchTrig", "number of tracks matched with trigger", 20, 0., 20.);
  Add2ESDsList(hESDMatchTrig, kESDMatchTrig,!forExpert);
  
  TH1F* hESDMomentum = new TH1F("hESDMomentum", "P distribution", 300, 0., 300);
  Add2ESDsList(hESDMomentum, kESDMomentum,forExpert);

  TH1F* hESDPt = new TH1F("hESDPt", "Pt distribution", 200, 0., 50);
  Add2ESDsList(hESDPt, kESDPt,forExpert);

  TH1F* hESDRapidity = new TH1F("hESDRapidity", "rapidity distribution", 200, -4.5, -2.);
  Add2ESDsList(hESDRapidity, kESDRapidity,forExpert);

  TH1F* hESDChi2 = new TH1F("hESDChi2", "normalized chi2 distribution", 500, 0., 50.);
  Add2ESDsList(hESDChi2, kESDChi2,forExpert);
  
  // cluster info
  for (Int_t i = 0; i < nCh; i++) {
    Float_t rMax = AliMUONConstants::Rmax(i/2);
    TH2F* hESDClusterHitMap = new TH2F(Form("hESDClusterHitMap%d",i+1), Form("cluster position distribution in chamber %d",i+1),
				       100, -rMax, rMax, 100, -rMax, rMax);
    Add2ESDsList(hESDClusterHitMap, kESDClusterHitMap+i,forExpert);
  }
  
  TH1F* hESDnClustersPerTrack = new TH1F("hESDnClustersPerTrack", "number of clusters per track", 20, 0., 20.);
  Add2ESDsList(hESDnClustersPerTrack, kESDnClustersPerTrack,!forExpert);
  
  TH1F* hESDnClustersPerCh = new TH1F("hESDnClustersPerCh", "number of clusters per chamber per track;chamber ID", nCh, 0, nCh);
  hESDnClustersPerCh->SetFillColor(kRed);
  Add2ESDsList(hESDnClustersPerCh, kESDnClustersPerCh,forExpert);
  
  TH1F* hESDnClustersPerDE = new TH1F("hESDnClustersPerDE", "number of clusters per DE per track;DetElem ID", nDE, 0, nDE);
  hESDnClustersPerDE->SetFillColor(kRed);
  Add2ESDsList(hESDnClustersPerDE, kESDnClustersPerDE,forExpert);
  
  TH1F* hESDClusterCharge = new TH1F("hESDClusterCharge", "cluster charge distribution", 500, 0., 5000.);
  Add2ESDsList(hESDClusterCharge, kESDClusterCharge,forExpert);
  
  for (Int_t i = 0; i < nCh; i++) {
    TH1F* hESDClusterChargeInCh = new TH1F(Form("hESDClusterChargeInCh%d",i+1), Form("cluster charge distribution in chamber %d",i+1), 500, 0., 5000.);
    Add2ESDsList(hESDClusterChargeInCh, kESDClusterChargeInCh+i,forExpert);
  }
  
  TH1F* hESDClusterChargePerDE = new TH1F("hESDClusterChargePerDE", "cluster mean charge per DE;DetElem ID", nDE, 0, nDE);
  hESDClusterChargePerDE->SetOption("P");
  hESDClusterChargePerDE->SetMarkerStyle(kFullDotMedium);
  hESDClusterChargePerDE->SetMarkerColor(kRed);
  Add2ESDsList(hESDClusterChargePerDE, kESDClusterChargePerDE,forExpert);
  
  TH1F* hESDClusterMult = new TH1F("hESDClusterMult", "cluster multiplicity distribution", 200, 0., 200.);
  Add2ESDsList(hESDClusterMult, kESDClusterMult,forExpert);
  
  for (Int_t i = 0; i < nCh; i++) {
    TH1F* hESDClusterMultInCh = new TH1F(Form("hESDClusterMultInCh%d",i+1), Form("cluster multiplicity distribution in chamber %d",i+1), 200, 0., 200.);
    Add2ESDsList(hESDClusterMultInCh, kESDClusterMultInCh+i,forExpert);
  }
  
  TH1F* hESDClusterMultPerDE = new TH1F("hESDClusterMultPerDE", "cluster mean multiplicity per DE;DetElem ID", nDE, 0, nDE);
  hESDClusterMultPerDE->SetOption("P");
  hESDClusterMultPerDE->SetMarkerStyle(kFullDotMedium);
  hESDClusterMultPerDE->SetMarkerColor(kRed);
  Add2ESDsList(hESDClusterMultPerDE, kESDClusterMultPerDE,forExpert);
  
  // cluster - track info
  TH1F* hESDResidualX = new TH1F("hESDResidualX", "cluster-track residual-X distribution", 1000, -5., 5.);
  Add2ESDsList(hESDResidualX, kESDResidualX,forExpert);
  
  TH1F* hESDResidualY = new TH1F("hESDResidualY", "cluster-track residual-Y distribution", 1000, -1., 1.);
  Add2ESDsList(hESDResidualY, kESDResidualY,forExpert);
  
  for (Int_t i = 0; i < nCh; i++) {
    TH1F* hESDResidualXInCh = new TH1F(Form("hESDResidualXInCh%d",i+1), Form("cluster-track residual-X distribution in chamber %d",i+1), 1000, -5., 5.);
    Add2ESDsList(hESDResidualXInCh, kESDResidualXInCh+i,forExpert);
    
    TH1F* hESDResidualYInCh = new TH1F(Form("hESDResidualYInCh%d",i+1), Form("cluster-track residual-Y distribution in chamber %d",i+1), 1000, -1., 1.);
    Add2ESDsList(hESDResidualYInCh, kESDResidualYInCh+i,forExpert);
  }
  
  TH1F* hESDResidualXPerDEMean = new TH1F("hESDResidualXPerDEMean", "cluster-track residual-X per DE: mean;DetElem ID", nDE, 0, nDE);
  hESDResidualXPerDEMean->SetOption("P");
  hESDResidualXPerDEMean->SetMarkerStyle(kFullDotMedium);
  hESDResidualXPerDEMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDResidualXPerDEMean, kESDResidualXPerDEMean,forExpert);
  
  TH1F* hESDResidualYPerDEMean = new TH1F("hESDResidualYPerDEMean", "cluster-track residual-Y per DE: mean;DetElem ID", nDE, 0, nDE);
  hESDResidualYPerDEMean->SetOption("P");
  hESDResidualYPerDEMean->SetMarkerStyle(kFullDotMedium);
  hESDResidualYPerDEMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDResidualYPerDEMean, kESDResidualYPerDEMean,forExpert);
  
  TH1F* hESDResidualXPerDESigma = new TH1F("hESDResidualXPerDESigma", "cluster-track residual-X per DE: sigma;DetElem ID", nDE, 0, nDE);
  hESDResidualXPerDESigma->SetOption("P");
  hESDResidualXPerDESigma->SetMarkerStyle(kFullDotMedium);
  hESDResidualXPerDESigma->SetMarkerColor(kRed);
  Add2ESDsList(hESDResidualXPerDESigma, kESDResidualXPerDESigma,forExpert);
  
  TH1F* hESDResidualYPerDESigma = new TH1F("hESDResidualYPerDESigma", "cluster-track residual-Y per DE: sigma;DetElem ID", nDE, 0, nDE);
  hESDResidualYPerDESigma->SetOption("P");
  hESDResidualYPerDESigma->SetMarkerStyle(kFullDotMedium);
  hESDResidualYPerDESigma->SetMarkerColor(kRed);
  Add2ESDsList(hESDResidualYPerDESigma, kESDResidualYPerDESigma,forExpert);
  
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
  
  // load ESD event in the interface
  AliMUONESDInterface esdInterface;
  if (GetRecoParam()) AliMUONESDInterface::ResetTracker(GetRecoParam());
  else AliError("Unable to get recoParam: use default ones for residual calculation");
  esdInterface.LoadEvent(*esd);
  
  GetESDsData(kESDnTracks)->Fill(esdInterface.GetNTracks());
  
  Int_t nTrackMatchTrig = 0;
  
  // loop over tracks
  Int_t nTracks = (Int_t) esd->GetNumberOfMuonTracks(); 
  for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {
    
    // get the ESD track and skip "ghosts"
    AliESDMuonTrack* esdTrack = esd->GetMuonTrack(iTrack);
    if (!esdTrack->ContainTrackerData()) continue;
    
    // get corresponding MUON track
    AliMUONTrack* track = esdInterface.FindTrack(esdTrack->GetUniqueID());
    
    if (esdTrack->ContainTriggerData()) nTrackMatchTrig++;
    
    GetESDsData(kESDMomentum)->Fill(esdTrack->P());
    GetESDsData(kESDPt)->Fill(esdTrack->Pt());
    GetESDsData(kESDRapidity)->Fill(esdTrack->Y());
    GetESDsData(kESDChi2)->Fill(track->GetNormalizedChi2());
    GetESDsData(kESDnClustersPerTrack)->Fill(track->GetNClusters());
    
    // loop over clusters
    AliMUONTrackParam* trackParam = static_cast<AliMUONTrackParam*>(track->GetTrackParamAtCluster()->First());
    while (trackParam) {
      
      AliMUONVCluster* cluster = trackParam->GetClusterPtr();
      Int_t chId = cluster->GetChamberId();
      Int_t deID = cluster->GetDetElemId();
      Double_t residualX = cluster->GetX() - trackParam->GetNonBendingCoor();
      Double_t residualY = cluster->GetY() - trackParam->GetBendingCoor();
      
      GetESDsData(kESDClusterHitMap+chId)->Fill(cluster->GetX(), cluster->GetY());
      
      GetESDsData(kESDnClustersPerCh)->Fill(chId);
      GetESDsData(kESDnClustersPerDE)->Fill(deID);
      
      GetESDsData(kESDClusterCharge)->Fill(cluster->GetCharge());
      GetESDsData(kESDClusterChargeInCh+chId)->Fill(cluster->GetCharge());
      GetESDsData(kESDClusterChargePerDE)->Fill(deID, cluster->GetCharge());
      
      if (cluster->GetNDigits() > 0) { // discard clusters with pad not stored in ESD
        GetESDsData(kESDClusterMult)->Fill(cluster->GetNDigits());
        GetESDsData(kESDClusterMultInCh+chId)->Fill(cluster->GetNDigits());
	GetESDsData(kESDClusterMultPerDE)->Fill(deID, cluster->GetNDigits());
      }
      
      GetESDsData(kESDResidualX)->Fill(residualX);
      GetESDsData(kESDResidualY)->Fill(residualY);
      GetESDsData(kESDResidualXInCh+chId)->Fill(residualX);
      GetESDsData(kESDResidualYInCh+chId)->Fill(residualY);
      GetESDsData(kESDResidualXPerDEMean)->Fill(deID, residualX);
      GetESDsData(kESDResidualYPerDEMean)->Fill(deID, residualY);
      GetESDsData(kESDResidualXPerDESigma)->Fill(deID, residualX*residualX);
      GetESDsData(kESDResidualYPerDESigma)->Fill(deID, residualY*residualY);
      
      trackParam = static_cast<AliMUONTrackParam*>(track->GetTrackParamAtCluster()->After(trackParam));
    }
    
  }
  
  GetESDsData(kESDMatchTrig)->Fill(nTrackMatchTrig);
  
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
