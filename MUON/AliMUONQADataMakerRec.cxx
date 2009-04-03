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
#include "AliMUONTrackerDataMaker.h"
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
#include "AliMUONCalibrationData.h"
#include "AliMpBusPatch.h"
#include "AliMpCDB.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpLocalBoard.h"
#include "AliMpStationType.h"
#include "AliMpTriggerCrate.h"
#include "AliMpDCSNamer.h"
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
#include "AliDCSValue.h"

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
void AliMUONQADataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray** list)
{
  ///Detector specific actions at end of cycle
  
  AliCodeTimerAuto("");
  
  // Display trigger histos in a more user friendly way
  DisplayTriggerInfo(task);
  
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    SetEventSpecie(specie) ; 
    if ( task == AliQA::kRAWS && fTrackerDataMaker ) 
      {
        TIter next(list[specie]);
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
            list[specie]->AddAt(fTrackerDataMaker->Data(),(Int_t)kTrackerData);
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
    
    // Normalize RecPoints histos
    if ( task == AliQA::kRECPOINTS ) {
      
      TH1* hTrackerClusterChargePerChMean = GetRecPointsData(kTrackerClusterChargePerChMean);
      TH1* hTrackerClusterChargePerChSigma = GetRecPointsData(kTrackerClusterChargePerChSigma);
      TH1* hTrackerClusterMultiplicityPerChMean = GetRecPointsData(kTrackerClusterMultiplicityPerChMean);
      TH1* hTrackerClusterMultiplicityPerChSigma = GetRecPointsData(kTrackerClusterMultiplicityPerChSigma);
      TH1* hTrackerClusterChargePerDEMean = GetRecPointsData(kTrackerClusterChargePerDEMean);
      TH1* hTrackerClusterMultiplicityPerDEMean = GetRecPointsData(kTrackerClusterMultiplicityPerDEMean);
      
      // loop over chambers
      for (Int_t iCh = 0; iCh < AliMUONConstants::NTrackingCh(); iCh++) {
	
	TH1* hTrackerClusterChargePerChamber = GetRecPointsData(kTrackerClusterChargePerChamber+iCh);
	Double_t sigmaCharge = hTrackerClusterChargePerChamber->GetRMS();
	hTrackerClusterChargePerChMean->SetBinContent(iCh+1, hTrackerClusterChargePerChamber->GetMean());
	hTrackerClusterChargePerChMean->SetBinError(iCh+1, hTrackerClusterChargePerChamber->GetMeanError());
	hTrackerClusterChargePerChSigma->SetBinContent(iCh+1, sigmaCharge);
	hTrackerClusterChargePerChSigma->SetBinError(iCh+1, hTrackerClusterChargePerChamber->GetRMSError());
	
	TH1* hTrackerClusterMultiplicityPerChamber = GetRecPointsData(kTrackerClusterMultiplicityPerChamber+iCh);
	Double_t sigmaSize = hTrackerClusterMultiplicityPerChamber->GetRMS();
	hTrackerClusterMultiplicityPerChMean->SetBinContent(iCh+1, hTrackerClusterMultiplicityPerChamber->GetMean());
	hTrackerClusterMultiplicityPerChMean->SetBinError(iCh+1, hTrackerClusterMultiplicityPerChamber->GetMeanError());
	hTrackerClusterMultiplicityPerChSigma->SetBinContent(iCh+1, sigmaSize);
	hTrackerClusterMultiplicityPerChSigma->SetBinError(iCh+1, hTrackerClusterMultiplicityPerChamber->GetRMSError());
	
	// loop over DE into chamber iCh
	AliMpDEIterator it;
	it.First(iCh);
	while ( !it.IsDone()) {
	  
	  Int_t iDE = it.CurrentDEId();
	  
	  TH1* hTrackerClusterChargePerDE = GetRecPointsData(kTrackerClusterChargePerDE+iDE);
	  hTrackerClusterChargePerDEMean->SetBinContent(iDE+1, hTrackerClusterChargePerDE->GetMean());
	  Double_t nClusters = hTrackerClusterChargePerDE->GetEntries();
	  if (nClusters > 1) hTrackerClusterChargePerDEMean->SetBinError(iDE+1, sigmaCharge/TMath::Sqrt(nClusters));
	  else hTrackerClusterChargePerDEMean->SetBinError(iDE+1, hTrackerClusterChargePerChamber->GetXaxis()->GetXmax());
	  
	  TH1* hTrackerClusterMultiplicityPerDE = GetRecPointsData(kTrackerClusterMultiplicityPerDE+iDE);
	  hTrackerClusterMultiplicityPerDEMean->SetBinContent(iDE+1, hTrackerClusterMultiplicityPerDE->GetMean());
	  nClusters = hTrackerClusterMultiplicityPerDE->GetEntries();
	  if (nClusters > 1) hTrackerClusterMultiplicityPerDEMean->SetBinError(iDE+1, sigmaSize/TMath::Sqrt(nClusters));
	  else hTrackerClusterMultiplicityPerDEMean->SetBinError(iDE+1, hTrackerClusterMultiplicityPerChamber->GetXaxis()->GetXmax());
	  
	  it.Next();
	}
	  
      }
      
    }
    
    // Normalize ESD histos
    if ( task == AliQA::kESDS ) {
      
      Double_t nTracks = GetESDsData(kESDnClustersPerTrack)->GetEntries();
      if (nTracks <= 0) continue;
      
      TH1* hESDnClustersPerCh = GetESDsData(kESDnClustersPerCh);
      TH1* hESDnClustersPerDE = GetESDsData(kESDnClustersPerDE);
      TH1* hESDClusterChargePerChMean = GetESDsData(kESDClusterChargePerChMean);
      TH1* hESDClusterChargePerChSigma = GetESDsData(kESDClusterChargePerChSigma);
      TH1* hESDClusterSizePerChMean = GetESDsData(kESDClusterSizePerChMean);
      TH1* hESDClusterSizePerChSigma = GetESDsData(kESDClusterSizePerChSigma);
      TH1* hESDResidualXPerChMean = GetESDsData(kESDResidualXPerChMean);
      TH1* hESDResidualXPerChSigma = GetESDsData(kESDResidualXPerChSigma);
      TH1* hESDResidualYPerChMean = GetESDsData(kESDResidualYPerChMean);
      TH1* hESDResidualYPerChSigma = GetESDsData(kESDResidualYPerChSigma);
      TH1* hESDLocalChi2XPerChMean = GetESDsData(kESDLocalChi2XPerChMean);
      TH1* hESDLocalChi2YPerChMean = GetESDsData(kESDLocalChi2YPerChMean);
      TH1* hESDClusterChargePerDE = GetESDsData(kESDClusterChargePerDE);
      TH1* hESDClusterSizePerDE = GetESDsData(kESDClusterSizePerDE);
      TH1* hESDResidualXPerDEMean = GetESDsData(kESDResidualXPerDEMean);
      TH1* hESDResidualXPerDESigma = GetESDsData(kESDResidualXPerDESigma);
      TH1* hESDResidualYPerDEMean = GetESDsData(kESDResidualYPerDEMean);
      TH1* hESDResidualYPerDESigma = GetESDsData(kESDResidualYPerDESigma);
      TH1* hESDLocalChi2XPerDEMean = GetESDsData(kESDLocalChi2XPerDEMean);
      TH1* hESDLocalChi2YPerDEMean = GetESDsData(kESDLocalChi2YPerDEMean);
      TH1* hESDnTotClustersPerCh = GetESDsData(kESDnTotClustersPerCh);
      TH1* hESDnTotClustersPerDE = GetESDsData(kESDnTotClustersPerDE);
      TH1* hESDnTotFullClustersPerDE = GetESDsData(kESDnTotFullClustersPerDE);
      TH1* hESDSumClusterChargePerDE = GetESDsData(kESDSumClusterChargePerDE);
      TH1* hESDSumClusterSizePerDE = GetESDsData(kESDSumClusterSizePerDE);
      TH1* hESDSumResidualXPerDE = GetESDsData(kESDSumResidualXPerDE);
      TH1* hESDSumResidualYPerDE = GetESDsData(kESDSumResidualYPerDE);
      TH1* hESDSumResidualX2PerDE = GetESDsData(kESDSumResidualX2PerDE);
      TH1* hESDSumResidualY2PerDE = GetESDsData(kESDSumResidualY2PerDE);
      TH1* hESDSumLocalChi2XPerDE = GetESDsData(kESDSumLocalChi2XPerDE);
      TH1* hESDSumLocalChi2YPerDE = GetESDsData(kESDSumLocalChi2YPerDE);
			  
      hESDnClustersPerCh->Reset();
      hESDnClustersPerDE->Reset();
      hESDnClustersPerCh->Add(hESDnTotClustersPerCh, 1./nTracks);
      hESDnClustersPerDE->Add(hESDnTotClustersPerDE, 1./nTracks);
      
      // loop over chambers
      for (Int_t iCh = 0; iCh < AliMUONConstants::NTrackingCh(); iCh++) {
	
	TH1* hESDClusterChargeInCh = GetESDsData(kESDClusterChargeInCh+iCh);
	Double_t sigmaCharge = hESDClusterChargeInCh->GetRMS();
	hESDClusterChargePerChMean->SetBinContent(iCh+1, hESDClusterChargeInCh->GetMean());
	hESDClusterChargePerChMean->SetBinError(iCh+1, hESDClusterChargeInCh->GetMeanError());
	hESDClusterChargePerChSigma->SetBinContent(iCh+1, sigmaCharge);
	hESDClusterChargePerChSigma->SetBinError(iCh+1, hESDClusterChargeInCh->GetRMSError());
	
	TH1* hESDClusterSizeInCh = GetESDsData(kESDClusterSizeInCh+iCh);
	Double_t sigmaSize = hESDClusterSizeInCh->GetRMS();
	hESDClusterSizePerChMean->SetBinContent(iCh+1, hESDClusterSizeInCh->GetMean());
	hESDClusterSizePerChMean->SetBinError(iCh+1, hESDClusterSizeInCh->GetMeanError());
	hESDClusterSizePerChSigma->SetBinContent(iCh+1, sigmaSize);
	hESDClusterSizePerChSigma->SetBinError(iCh+1, hESDClusterSizeInCh->GetRMSError());
	
	TH1* hESDResidualXInCh = GetESDsData(kESDResidualXInCh+iCh);
	Double_t sigmaResidualX = hESDResidualXInCh->GetRMS();
	hESDResidualXPerChMean->SetBinContent(iCh+1, hESDResidualXInCh->GetMean());
	hESDResidualXPerChMean->SetBinError(iCh+1, hESDResidualXInCh->GetMeanError());
	hESDResidualXPerChSigma->SetBinContent(iCh+1, sigmaResidualX);
	hESDResidualXPerChSigma->SetBinError(iCh+1, hESDResidualXInCh->GetRMSError());
	
	TH1* hESDResidualYInCh = GetESDsData(kESDResidualYInCh+iCh);
	Double_t sigmaResidualY = hESDResidualYInCh->GetRMS();
	hESDResidualYPerChMean->SetBinContent(iCh+1, hESDResidualYInCh->GetMean());
	hESDResidualYPerChMean->SetBinError(iCh+1, hESDResidualYInCh->GetMeanError());
	hESDResidualYPerChSigma->SetBinContent(iCh+1, sigmaResidualY);
	hESDResidualYPerChSigma->SetBinError(iCh+1, hESDResidualYInCh->GetRMSError());
	
	TH1* hESDLocalChi2XInCh = GetESDsData(kESDLocalChi2XInCh+iCh);
	Double_t sigmaLocalChi2X = hESDLocalChi2XInCh->GetRMS();
	hESDLocalChi2XPerChMean->SetBinContent(iCh+1, hESDLocalChi2XInCh->GetMean());
	hESDLocalChi2XPerChMean->SetBinError(iCh+1, hESDLocalChi2XInCh->GetMeanError());
	
	TH1* hESDLocalChi2YInCh = GetESDsData(kESDLocalChi2YInCh+iCh);
	Double_t sigmaLocalChi2Y = hESDLocalChi2YInCh->GetRMS();
	hESDLocalChi2YPerChMean->SetBinContent(iCh+1, hESDLocalChi2YInCh->GetMean());
	hESDLocalChi2YPerChMean->SetBinError(iCh+1, hESDLocalChi2YInCh->GetMeanError());
	
	// loop over DE into chamber iCh
	AliMpDEIterator it;
	it.First(iCh);
	while ( !it.IsDone()) {
	  
	  Int_t iDE = it.CurrentDEId();
	  
	  Double_t nClusters = hESDnTotClustersPerDE->GetBinContent(iDE+1);
	  if (nClusters > 1) {
	    
	    hESDClusterChargePerDE->SetBinContent(iDE+1, hESDSumClusterChargePerDE->GetBinContent(iDE+1)/nClusters);
	    hESDClusterChargePerDE->SetBinError(iDE+1, sigmaCharge/TMath::Sqrt(nClusters));
	    
	    Double_t meanResX = hESDSumResidualXPerDE->GetBinContent(iDE+1)/nClusters;
	    hESDResidualXPerDEMean->SetBinContent(iDE+1, meanResX);
	    hESDResidualXPerDEMean->SetBinError(iDE+1, sigmaResidualX/TMath::Sqrt(nClusters));
	    hESDResidualXPerDESigma->SetBinContent(iDE+1, TMath::Sqrt(hESDSumResidualX2PerDE->GetBinContent(iDE+1)/nClusters - meanResX*meanResX));
	    hESDResidualXPerDESigma->SetBinError(iDE+1, sigmaResidualX/TMath::Sqrt(2.*nClusters));
	    
	    Double_t meanResY = hESDSumResidualYPerDE->GetBinContent(iDE+1)/nClusters;
	    hESDResidualYPerDEMean->SetBinContent(iDE+1, meanResY);
	    hESDResidualYPerDEMean->SetBinError(iDE+1, sigmaResidualY/TMath::Sqrt(nClusters));
	    hESDResidualYPerDESigma->SetBinContent(iDE+1, TMath::Sqrt(hESDSumResidualY2PerDE->GetBinContent(iDE+1)/nClusters - meanResY*meanResY));
	    hESDResidualYPerDESigma->SetBinError(iDE+1, sigmaResidualY/TMath::Sqrt(2.*nClusters));
	    
	    hESDLocalChi2XPerDEMean->SetBinContent(iDE+1, hESDSumLocalChi2XPerDE->GetBinContent(iDE+1)/nClusters);
	    hESDLocalChi2XPerDEMean->SetBinError(iDE+1, sigmaLocalChi2X/TMath::Sqrt(nClusters));
	    
	    hESDLocalChi2YPerDEMean->SetBinContent(iDE+1, hESDSumLocalChi2YPerDE->GetBinContent(iDE+1)/nClusters);
	    hESDLocalChi2YPerDEMean->SetBinError(iDE+1, sigmaLocalChi2Y/TMath::Sqrt(nClusters));
	    
	  } else {
	    
	    hESDClusterChargePerDE->SetBinContent(iDE+1, hESDSumClusterChargePerDE->GetBinContent(iDE+1));
	    hESDClusterChargePerDE->SetBinError(iDE+1, hESDClusterChargeInCh->GetXaxis()->GetXmax());
	    
	    hESDResidualXPerDEMean->SetBinContent(iDE+1, hESDSumResidualXPerDE->GetBinContent(iDE+1));
	    hESDResidualXPerDEMean->SetBinError(iDE+1, hESDResidualXInCh->GetXaxis()->GetXmax());
	    hESDResidualXPerDESigma->SetBinContent(iDE+1, 0.);
	    hESDResidualXPerDESigma->SetBinError(iDE+1, hESDResidualXInCh->GetXaxis()->GetXmax());
	    
	    hESDResidualYPerDEMean->SetBinContent(iDE+1, hESDSumResidualYPerDE->GetBinContent(iDE+1));
	    hESDResidualYPerDEMean->SetBinError(iDE+1, hESDResidualYInCh->GetXaxis()->GetXmax());
	    hESDResidualYPerDESigma->SetBinContent(iDE+1, 0.);
	    hESDResidualYPerDESigma->SetBinError(iDE+1, hESDResidualYInCh->GetXaxis()->GetXmax());
	    
	    hESDLocalChi2XPerDEMean->SetBinContent(iDE+1, hESDSumLocalChi2XPerDE->GetBinContent(iDE+1));
	    hESDLocalChi2XPerDEMean->SetBinError(iDE+1, hESDLocalChi2XInCh->GetXaxis()->GetXmax());
	    
	    hESDLocalChi2YPerDEMean->SetBinContent(iDE+1, hESDSumLocalChi2YPerDE->GetBinContent(iDE+1));
	    hESDLocalChi2YPerDEMean->SetBinError(iDE+1, hESDLocalChi2YInCh->GetXaxis()->GetXmax());
	    
	  }
	  
	  Double_t nFullClusters = hESDnTotFullClustersPerDE->GetBinContent(iDE+1);
	  if (nFullClusters > 1) {
	    
	    hESDClusterSizePerDE->SetBinContent(iDE+1, hESDSumClusterSizePerDE->GetBinContent(iDE+1)/nFullClusters);
	    hESDClusterSizePerDE->SetBinError(iDE+1, sigmaSize/TMath::Sqrt(nFullClusters));
	    
	  } else {
	    
	    hESDClusterSizePerDE->SetBinContent(iDE+1, hESDSumClusterSizePerDE->GetBinContent(iDE+1));
	    hESDClusterSizePerDE->SetBinError(iDE+1, hESDClusterSizeInCh->GetXaxis()->GetXmax());
	    
	  }
	  
	  it.Next();
	}
	
      }
      
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
			Add2RawsList(h5, kTriggerScalersDisplay + AliMpConstants::NofTriggerChambers()*iCath + iChamber,!forExpert);
		}
	}

	TH2F* h6 = new TH2F("hTriggerRPCI", "Trigger RPC currents",
			    4, 10.5, 14.5,
			    18, -0.5, 17.5);
	h6->GetXaxis()->SetTitle("Chamber");
	h6->GetYaxis()->SetTitle("RPC");
	Add2RawsList(h6, kTriggerRPCi, forExpert);

	TH2F* h7 = new TH2F("hTriggerRPCHV", "Trigger RPC HV",
			    4, 10.5, 14.5,
			    18, -0.5, 17.5);
	h7->GetXaxis()->SetTitle("Chamber");
	h7->GetYaxis()->SetTitle("RPC");
	Add2RawsList(h7, kTriggerRPChv, forExpert);
	
	for(Int_t iChamber=0; iChamber<AliMpConstants::NofTriggerChambers(); iChamber++){
	  histoName = Form("hRPCIChamber%i", 11+iChamber);
	  histoTitle = Form("Chamber %i: RPC Currents (#muA)", 11+iChamber);
	  TH2F* h8 = (TH2F*)triggerDisplay.GetEmptyDisplayHisto(histoName, AliMUONTriggerDisplay::kDisplaySlats, 
								0, iChamber, histoTitle);
	  Add2RawsList(h8, kTriggerIDisplay + iChamber, !forExpert);
	}

	for(Int_t iChamber=0; iChamber<AliMpConstants::NofTriggerChambers(); iChamber++){
	  histoName = Form("hRPCHVChamber%i", 11+iChamber);
	  histoTitle = Form("Chamber %i: RPC HV (V)", 11+iChamber);
	  TH2F* h9 = (TH2F*)triggerDisplay.GetEmptyDisplayHisto(histoName, AliMUONTriggerDisplay::kDisplaySlats, 
								0, iChamber, histoTitle);
	  Add2RawsList(h9, kTriggerHVDisplay + iChamber, !forExpert);
	}

	TH1F* h10 = new TH1F("hTriggerScalersTime", "Trigger scalers acquisition time", 1, 0.5, 1.5);
	h10->GetXaxis()->SetBinLabel(1, "One-bin histogram: bin is filled at each scaler event.");
	h10->GetYaxis()->SetTitle("Cumulated scaler time (s)");
	Add2RawsList(h10, kTriggerScalersTime, !forExpert);
	
  Int_t nbp(0);
  TIter next(AliMpDDLStore::Instance()->CreateBusPatchIterator());
  while (next())
  {
    ++nbp;
  }
  
  TH1* hbp = new TH1F("hTrackerBusPatchOccupancy","Occupancy of bus patches",
                      nbp,-0.5,nbp-0.5);
  
  Add2RawsList(hbp,kTrackerBusPatchOccupancy,!forExpert);

  const Bool_t histogram(kFALSE);

  fTrackerDataMaker = new AliMUONTrackerDataMaker(GetRecoParam(),
                                                  AliCDBManager::Instance()->GetRun(),
                                                  0x0,
                                                  AliCDBManager::Instance()->GetDefaultStorage()->GetURI(),
                                                  "NOGAIN",
                                                  histogram,
                                                  0.0,0.0);
		
  fTrackerDataMaker->Data()->DisableChannelLevel(); // to save up disk space, we only store starting at the manu level
	
  fTrackerDataMaker->SetRunning(kTRUE);
  
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
  TH1I *h1I;
  TH1F *h1F;
  TH2F *h2F;
  
  // histograms per chamber
  Int_t nCh = AliMpConstants::NofTrackingChambers();
  for ( Int_t i = 0; i < nCh; ++i ) 
  {
    h1I = new TH1I(Form("hTrackerClusterMultiplicityForChamber%d",i+1), Form("cluster size distribution in chamber %d;size (n_{pads})",i+1), 100,0,100);
    Add2RecPointsList(h1I,kTrackerClusterMultiplicityPerChamber+i,forExpert);
    
    h1I = new TH1I(Form("hTrackerClusterChargeForChamber%d",i+1), Form("cluster charge distribution in chamber %d;charge (ADC counts)",i+1), 500,0,5000);
    Add2RecPointsList(h1I,kTrackerClusterChargePerChamber+i,forExpert);
    
    Float_t rMax = AliMUONConstants::Rmax(i/2);
    h2F = new TH2F(Form("hTrackerClusterHitMapForChamber%d",i+1), Form("cluster position distribution in chamber %d;X (cm);Y (cm)",i+1), 100, -rMax, rMax, 100, -rMax, rMax);
    Add2RecPointsList(h2F, kTrackerClusterHitMapPerChamber+i,forExpert);
  }
  
  // summary histograms per chamber
  h1I = new TH1I("hTrackerNumberOfClustersPerChamber", "Number of clusters per chamber;chamber ID;n_{clusters}", nCh,-0.5,nCh-0.5);
  Add2RecPointsList(h1I,kTrackerNumberOfClustersPerChamber,forExpert);
  
  h1F = new TH1F("hTrackerClusterMultiplicityPerChMean", "cluster mean size per chamber;chamber ID;<size> (n_{pads})", nCh,-0.5,nCh-0.5);
  h1F->SetOption("P");
  h1F->SetMarkerStyle(kFullDotMedium);
  h1F->SetMarkerColor(kRed);
  Add2RecPointsList(h1F, kTrackerClusterMultiplicityPerChMean,forExpert);
  
  h1F = new TH1F("hTrackerClusterMultiplicityPerChSigma", "cluster size dispersion per chamber;chamber ID;#sigma_{size} (n_{pads})", nCh,-0.5,nCh-0.5);
  h1F->SetOption("P");
  h1F->SetMarkerStyle(kFullDotMedium);
  h1F->SetMarkerColor(kRed);
  Add2RecPointsList(h1F, kTrackerClusterMultiplicityPerChSigma,forExpert);
  
  h1F = new TH1F("hTrackerClusterChargePerChMean", "cluster mean charge per chamber;chamber ID;<charge> (ADC counts)", nCh,-0.5,nCh-0.5);
  h1F->SetOption("P");
  h1F->SetMarkerStyle(kFullDotMedium);
  h1F->SetMarkerColor(kRed);
  Add2RecPointsList(h1F, kTrackerClusterChargePerChMean,forExpert);
  
  h1F = new TH1F("hTrackerClusterChargePerChSigma", "cluster charge dispersion per chamber;chamber ID;#sigma_{charge} (ADC counts)", nCh,-0.5,nCh-0.5);
  h1F->SetOption("P");
  h1F->SetMarkerStyle(kFullDotMedium);
  h1F->SetMarkerColor(kRed);
  Add2RecPointsList(h1F, kTrackerClusterChargePerChSigma,forExpert);
  
  // histograms per DE
  Int_t ndes(0);
  AliMpDEIterator it;
  it.First();
  while ( !it.IsDone())
  {
    Int_t detElemId = it.CurrentDEId();
    
    if ( AliMpDEManager::GetStationType(detElemId) != AliMp::kStationTrigger )
    {
      ndes = TMath::Max(ndes,detElemId);
      
      h1I = new TH1I(Form("hTrackerClusterMultiplicityForDE%04d",detElemId), Form("cluster size distribution in detection element %d;size (n_{pads})",detElemId), 100,0,100);
      Add2RecPointsList(h1I,kTrackerClusterMultiplicityPerDE+detElemId,forExpert);
      
      h1I = new TH1I(Form("hTrackerClusterChargeForDE%04d",detElemId), Form("cluster charge distribution in detection element %d;charge (ADC counts)",detElemId), 500,0,5000);
      Add2RecPointsList(h1I,kTrackerClusterChargePerDE+detElemId,forExpert);
    }
    
    it.Next();
  }
  
  // summary histograms per DE
  h1I = new TH1I("hTrackerNumberOfClustersPerDE", "Number of clusters per detection element;DetElem ID;n_{clusters}", ndes+1,-0.5,ndes+0.5);
  Add2RecPointsList(h1I, kTrackerNumberOfClustersPerDE,!forExpert);
  
  h1F = new TH1F("hTrackerClusterMultiplicityPerDEMean", "cluster mean size per DE;DetElem ID;<size> (n_{pads})", ndes+1,-0.5,ndes+0.5);
  h1F->SetOption("P");
  h1F->SetMarkerStyle(kFullDotMedium);
  h1F->SetMarkerColor(kRed);
  Add2RecPointsList(h1F, kTrackerClusterMultiplicityPerDEMean,forExpert);
  
  h1F = new TH1F("hTrackerClusterChargePerDEMean", "cluster mean charge per DE;DetElem ID;<charge> (ADC counts)", ndes+1,-0.5,ndes+0.5);
  h1F->SetOption("P");
  h1F->SetMarkerStyle(kFullDotMedium);
  h1F->SetMarkerColor(kRed);
  Add2RecPointsList(h1F, kTrackerClusterChargePerDEMean,forExpert);
  
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
	Add2RecPointsList(h3, kTriggerDigitsDisplay + AliMpConstants::NofTriggerChambers()*iCath + iChamber,!forExpert);
      }
    }

    TH2F* h4 = (TH2F*)triggerDisplay.GetEmptyDisplayHisto("hFiredBoardsDisplay", AliMUONTriggerDisplay::kDisplayBoards,
							  0, 0, "Fired boards");
    Add2RecPointsList(h4, kTriggerBoardsDisplay,!forExpert);
	
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
  TH1F* hESDnTracks = new TH1F("hESDnTracks", "number of tracks;n_{tracks}", 20, 0., 20.);
  Add2ESDsList(hESDnTracks, kESDnTracks,!forExpert);

  TH1F* hESDMatchTrig = new TH1F("hESDMatchTrig", "number of tracks matched with trigger;n_{tracks}", 20, 0., 20.);
  Add2ESDsList(hESDMatchTrig, kESDMatchTrig,!forExpert);
  
  TH1F* hESDMomentum = new TH1F("hESDMomentum", "momentum distribution;p (GeV/c)", 300, 0., 300);
  Add2ESDsList(hESDMomentum, kESDMomentum,forExpert);

  TH1F* hESDPt = new TH1F("hESDPt", "transverse momentum distribution;p_{t} (GeV/c)", 200, 0., 50);
  Add2ESDsList(hESDPt, kESDPt,forExpert);

  TH1F* hESDRapidity = new TH1F("hESDRapidity", "rapidity distribution;rapidity", 200, -4.5, -2.);
  Add2ESDsList(hESDRapidity, kESDRapidity,forExpert);

  TH1F* hESDChi2 = new TH1F("hESDChi2", "normalized #chi^{2} distribution;#chi^{2} / ndf", 500, 0., 50.);
  Add2ESDsList(hESDChi2, kESDChi2,forExpert);
  
  TH1F* hESDProbChi2 = new TH1F("hESDProbChi2", "distribution of probability of #chi^{2};prob(#chi^{2})", 100, 0., 1.);
  Add2ESDsList(hESDProbChi2, kESDProbChi2,forExpert);
  
  // cluster info
  for (Int_t i = 0; i < nCh; i++) {
    Float_t rMax = AliMUONConstants::Rmax(i/2);
    TH2F* hESDClusterHitMap = new TH2F(Form("hESDClusterHitMap%d",i+1), Form("cluster position distribution in chamber %d;X (cm);Y (cm)",i+1),
				       100, -rMax, rMax, 100, -rMax, rMax);
    Add2ESDsList(hESDClusterHitMap, kESDClusterHitMap+i,forExpert);
  }
  
  TH1F* hESDnClustersPerTrack = new TH1F("hESDnClustersPerTrack", "number of associated clusters per track;n_{clusters}", 20, 0., 20.);
  Add2ESDsList(hESDnClustersPerTrack, kESDnClustersPerTrack,!forExpert);
  
  TH1F* hESDnClustersPerCh = new TH1F("hESDnClustersPerCh", "averaged number of clusters per chamber per track;chamber ID;<n_{clusters}>", nCh, -0.5, nCh-0.5);
  hESDnClustersPerCh->SetFillColor(kRed);
  Add2ESDsList(hESDnClustersPerCh, kESDnClustersPerCh,forExpert);
  
  TH1F* hESDnClustersPerDE = new TH1F("hESDnClustersPerDE", "averaged number of clusters per DE per track;DetElem ID;<n_{clusters}>", nDE+1, -0.5, nDE+0.5);
  hESDnClustersPerDE->SetFillColor(kRed);
  Add2ESDsList(hESDnClustersPerDE, kESDnClustersPerDE,forExpert);
  
  for (Int_t i = 0; i < nCh; i++) {
    TH1F* hESDClusterChargeInCh = new TH1F(Form("hESDClusterChargeInCh%d",i+1), Form("cluster charge distribution in chamber %d;charge (ADC counts)",i+1), 500, 0., 5000.);
    Add2ESDsList(hESDClusterChargeInCh, kESDClusterChargeInCh+i,forExpert);
  }
  
  TH1F* hESDClusterChargePerChMean = new TH1F("hESDClusterChargePerChMean", "cluster mean charge per chamber;chamber ID;<charge> (ADC counts)", nCh, -0.5, nCh-0.5);
  hESDClusterChargePerChMean->SetOption("P");
  hESDClusterChargePerChMean->SetMarkerStyle(kFullDotMedium);
  hESDClusterChargePerChMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDClusterChargePerChMean, kESDClusterChargePerChMean,forExpert);
  
  TH1F* hESDClusterChargePerChSigma = new TH1F("hESDClusterChargePerChSigma", "cluster charge dispersion per chamber;chamber ID;#sigma_{charge} (ADC counts)", nCh, -0.5, nCh-0.5);
  hESDClusterChargePerChSigma->SetOption("P");
  hESDClusterChargePerChSigma->SetMarkerStyle(kFullDotMedium);
  hESDClusterChargePerChSigma->SetMarkerColor(kRed);
  Add2ESDsList(hESDClusterChargePerChSigma, kESDClusterChargePerChSigma,forExpert);
  
  TH1F* hESDClusterChargePerDE = new TH1F("hESDClusterChargePerDE", "cluster mean charge per DE;DetElem ID;<charge> (ADC counts)", nDE+1, -0.5, nDE+0.5);
  hESDClusterChargePerDE->SetOption("P");
  hESDClusterChargePerDE->SetMarkerStyle(kFullDotMedium);
  hESDClusterChargePerDE->SetMarkerColor(kRed);
  Add2ESDsList(hESDClusterChargePerDE, kESDClusterChargePerDE,forExpert);
  
  for (Int_t i = 0; i < nCh; i++) {
    TH1F* hESDClusterSizeInCh = new TH1F(Form("hESDClusterSizeInCh%d",i+1), Form("cluster size distribution in chamber %d;size (n_{pads})",i+1), 200, 0., 200.);
    Add2ESDsList(hESDClusterSizeInCh, kESDClusterSizeInCh+i,forExpert);
  }
  
  TH1F* hESDClusterSizePerChMean = new TH1F("hESDClusterSizePerChMean", "cluster mean size per chamber;chamber ID;<size> (n_{pads})", nCh, -0.5, nCh-0.5);
  hESDClusterSizePerChMean->SetOption("P");
  hESDClusterSizePerChMean->SetMarkerStyle(kFullDotMedium);
  hESDClusterSizePerChMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDClusterSizePerChMean, kESDClusterSizePerChMean,forExpert);
  
  TH1F* hESDClusterSizePerChSigma = new TH1F("hESDClusterSizePerChSigma", "cluster size dispersion per chamber;chamber ID;#sigma_{size} (n_{pads})", nCh, -0.5, nCh-0.5);
  hESDClusterSizePerChSigma->SetOption("P");
  hESDClusterSizePerChSigma->SetMarkerStyle(kFullDotMedium);
  hESDClusterSizePerChSigma->SetMarkerColor(kRed);
  Add2ESDsList(hESDClusterSizePerChSigma, kESDClusterSizePerChSigma,forExpert);
  
  TH1F* hESDClusterSizePerDE = new TH1F("hESDClusterSizePerDE", "cluster mean size per DE;DetElem ID;<size> (n_{pads})", nDE+1, -0.5, nDE+0.5);
  hESDClusterSizePerDE->SetOption("P");
  hESDClusterSizePerDE->SetMarkerStyle(kFullDotMedium);
  hESDClusterSizePerDE->SetMarkerColor(kRed);
  Add2ESDsList(hESDClusterSizePerDE, kESDClusterSizePerDE,forExpert);
  
  // cluster - track info
  for (Int_t i = 0; i < nCh; i++) {
    TH1F* hESDResidualXInCh = new TH1F(Form("hESDResidualXInCh%d",i+1), Form("cluster-track residual-X distribution in chamber %d;#Delta_{X} (cm)",i+1), 1000, -5., 5.);
    Add2ESDsList(hESDResidualXInCh, kESDResidualXInCh+i,forExpert);
    
    TH1F* hESDResidualYInCh = new TH1F(Form("hESDResidualYInCh%d",i+1), Form("cluster-track residual-Y distribution in chamber %d;#Delta_{Y} (cm)",i+1), 1000, -1., 1.);
    Add2ESDsList(hESDResidualYInCh, kESDResidualYInCh+i,forExpert);
    
    TH1F* hESDLocalChi2XInCh = new TH1F(Form("hESDLocalChi2XInCh%d",i+1), Form("local chi2-X distribution in chamber %d;local #chi^{2}_{X}",i+1), 1000, 0., 25);
    Add2ESDsList(hESDLocalChi2XInCh, kESDLocalChi2XInCh+i,forExpert);
    
    TH1F* hESDLocalChi2YInCh = new TH1F(Form("hESDLocalChi2YInCh%d",i+1), Form("local chi2-Y distribution in chamber %d;local #chi^{2}_{Y}",i+1), 1000, 0., 25);
    Add2ESDsList(hESDLocalChi2YInCh, kESDLocalChi2YInCh+i,forExpert);
  }
  
  TH1F* hESDResidualXPerChMean = new TH1F("hESDResidualXPerChMean", "cluster-track residual-X per Ch: mean;chamber ID;<#Delta_{X}> (cm)", nCh, -0.5, nCh-0.5);
  hESDResidualXPerChMean->SetOption("P");
  hESDResidualXPerChMean->SetMarkerStyle(kFullDotMedium);
  hESDResidualXPerChMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDResidualXPerChMean, kESDResidualXPerChMean,forExpert);
  
  TH1F* hESDResidualYPerChMean = new TH1F("hESDResidualYPerChMean", "cluster-track residual-Y per Ch: mean;chamber ID;<#Delta_{Y}> (cm)", nCh, -0.5, nCh-0.5);
  hESDResidualYPerChMean->SetOption("P");
  hESDResidualYPerChMean->SetMarkerStyle(kFullDotMedium);
  hESDResidualYPerChMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDResidualYPerChMean, kESDResidualYPerChMean,forExpert);
  
  TH1F* hESDResidualXPerChSigma = new TH1F("hESDResidualXPerChSigma", "cluster-track residual-X per Ch: sigma;chamber ID;#sigma_{X} (cm)", nCh, -0.5, nCh-0.5);
  hESDResidualXPerChSigma->SetOption("P");
  hESDResidualXPerChSigma->SetMarkerStyle(kFullDotMedium);
  hESDResidualXPerChSigma->SetMarkerColor(kRed);
  Add2ESDsList(hESDResidualXPerChSigma, kESDResidualXPerChSigma,forExpert);
  
  TH1F* hESDResidualYPerChSigma = new TH1F("hESDResidualYPerChSigma", "cluster-track residual-Y per Ch: sigma;chamber ID;#sigma_{Y} (cm)", nCh, -0.5, nCh-0.5);
  hESDResidualYPerChSigma->SetOption("P");
  hESDResidualYPerChSigma->SetMarkerStyle(kFullDotMedium);
  hESDResidualYPerChSigma->SetMarkerColor(kRed);
  Add2ESDsList(hESDResidualYPerChSigma, kESDResidualYPerChSigma,forExpert);
  
  TH1F* hESDLocalChi2XPerChMean = new TH1F("hESDLocalChi2XPerChMean", "local chi2-X per Ch: mean;chamber ID;<local #chi^{2}_{X}>", nCh, -0.5, nCh-0.5);
  hESDLocalChi2XPerChMean->SetOption("P");
  hESDLocalChi2XPerChMean->SetMarkerStyle(kFullDotMedium);
  hESDLocalChi2XPerChMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDLocalChi2XPerChMean, kESDLocalChi2XPerChMean,forExpert);
  
  TH1F* hESDLocalChi2YPerChMean = new TH1F("hESDLocalChi2YPerChMean", "local chi2-Y per Ch: mean;chamber ID;<local #chi^{2}_{Y}>", nCh, -0.5, nCh-0.5);
  hESDLocalChi2YPerChMean->SetOption("P");
  hESDLocalChi2YPerChMean->SetMarkerStyle(kFullDotMedium);
  hESDLocalChi2YPerChMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDLocalChi2YPerChMean, kESDLocalChi2YPerChMean,forExpert);
  
  TH1F* hESDResidualXPerDEMean = new TH1F("hESDResidualXPerDEMean", "cluster-track residual-X per DE: mean;DetElem ID;<#Delta_{X}> (cm)", nDE+1, -0.5, nDE+0.5);
  hESDResidualXPerDEMean->SetOption("P");
  hESDResidualXPerDEMean->SetMarkerStyle(kFullDotMedium);
  hESDResidualXPerDEMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDResidualXPerDEMean, kESDResidualXPerDEMean,forExpert);
  
  TH1F* hESDResidualYPerDEMean = new TH1F("hESDResidualYPerDEMean", "cluster-track residual-Y per DE: mean;DetElem ID;<#Delta_{Y}> (cm)", nDE+1, -0.5, nDE+0.5);
  hESDResidualYPerDEMean->SetOption("P");
  hESDResidualYPerDEMean->SetMarkerStyle(kFullDotMedium);
  hESDResidualYPerDEMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDResidualYPerDEMean, kESDResidualYPerDEMean,forExpert);
  
  TH1F* hESDResidualXPerDESigma = new TH1F("hESDResidualXPerDESigma", "cluster-track residual-X per DE: sigma;DetElem ID;#sigma_{X} (cm)", nDE+1, -0.5, nDE+0.5);
  hESDResidualXPerDESigma->SetOption("P");
  hESDResidualXPerDESigma->SetMarkerStyle(kFullDotMedium);
  hESDResidualXPerDESigma->SetMarkerColor(kRed);
  Add2ESDsList(hESDResidualXPerDESigma, kESDResidualXPerDESigma,forExpert);
  
  TH1F* hESDResidualYPerDESigma = new TH1F("hESDResidualYPerDESigma", "cluster-track residual-Y per DE: sigma;DetElem ID;#sigma_{Y} (cm)", nDE+1, -0.5, nDE+0.5);
  hESDResidualYPerDESigma->SetOption("P");
  hESDResidualYPerDESigma->SetMarkerStyle(kFullDotMedium);
  hESDResidualYPerDESigma->SetMarkerColor(kRed);
  Add2ESDsList(hESDResidualYPerDESigma, kESDResidualYPerDESigma,forExpert);
  
  TH1F* hESDLocalChi2XPerDEMean = new TH1F("hESDLocalChi2XPerDEMean", "local chi2-X per DE: mean;DetElem ID;<local #chi^{2}_{X}>", nDE+1, -0.5, nDE+0.5);
  hESDLocalChi2XPerDEMean->SetOption("P");
  hESDLocalChi2XPerDEMean->SetMarkerStyle(kFullDotMedium);
  hESDLocalChi2XPerDEMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDLocalChi2XPerDEMean, kESDLocalChi2XPerDEMean,forExpert);
  
  TH1F* hESDLocalChi2YPerDEMean = new TH1F("hESDLocalChi2YPerDEMean", "local chi2-Y per DE: mean;DetElem ID;<local #chi^{2}_{Y}>", nDE+1, -0.5, nDE+0.5);
  hESDLocalChi2YPerDEMean->SetOption("P");
  hESDLocalChi2YPerDEMean->SetMarkerStyle(kFullDotMedium);
  hESDLocalChi2YPerDEMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDLocalChi2YPerDEMean, kESDLocalChi2YPerDEMean,forExpert);
  
  // intermediate histograms
  TH1F* hESDnTotClustersPerCh = new TH1F("hESDnTotClustersPerCh", "total number of associated clusters per chamber;chamber ID;#Sigma(n_{clusters})", nCh, -0.5, nCh-0.5);
  Add2ESDsList(hESDnTotClustersPerCh, kESDnTotClustersPerCh, forExpert);
  TH1F* hESDnTotClustersPerDE = new TH1F("hESDnTotClustersPerDE", "total number of associated clusters per DE;DetElem ID;#Sigma(n_{clusters})", nDE+1, -0.5, nDE+0.5);
  Add2ESDsList(hESDnTotClustersPerDE, kESDnTotClustersPerDE, forExpert);
  TH1F* hESDnTotFullClustersPerDE = new TH1F("hESDnTotFullClustersPerDE", "total number of associated clusters containing pad info per DE;DetElem ID;#Sigma(n_{full clusters})", nDE+1, -0.5, nDE+0.5);
  Add2ESDsList(hESDnTotFullClustersPerDE, kESDnTotFullClustersPerDE, forExpert);
  TH1F* hESDSumClusterChargePerDE = new TH1F("hESDSumClusterChargePerDE", "sum of cluster charge per DE;DetElem ID;#Sigma(charge) (ADC counts)", nDE+1, -0.5, nDE+0.5);
  Add2ESDsList(hESDSumClusterChargePerDE, kESDSumClusterChargePerDE, forExpert);
  TH1F* hESDSumClusterSizePerDE = new TH1F("hESDSumClusterSizePerDE", "sum of cluster size per DE;DetElem ID;#Sigma(size) (n_{pads})", nDE+1, -0.5, nDE+0.5);
  Add2ESDsList(hESDSumClusterSizePerDE, kESDSumClusterSizePerDE, forExpert);
  TH1F* hESDSumResidualXPerDE = new TH1F("hESDSumResidualXPerDE", "sum of cluster-track residual-X per DE;DetElem ID;#Sigma(#Delta_{X}) (cm)", nDE+1, -0.5, nDE+0.5);
  Add2ESDsList(hESDSumResidualXPerDE, kESDSumResidualXPerDE, forExpert);
  TH1F* hESDSumResidualYPerDE = new TH1F("hESDSumResidualYPerDE", "sum of cluster-track residual-Y per DE;DetElem ID;#Sigma(#Delta_{Y}) (cm)", nDE+1, -0.5, nDE+0.5);
  Add2ESDsList(hESDSumResidualYPerDE, kESDSumResidualYPerDE, forExpert);
  TH1F* hESDSumResidualX2PerDE = new TH1F("hESDSumResidualX2PerDE", "sum of cluster-track residual-X**2 per DE;DetElem ID;#Sigma(#Delta_{X}^{2}) (cm^{2})", nDE+1, -0.5, nDE+0.5);
  Add2ESDsList(hESDSumResidualX2PerDE, kESDSumResidualX2PerDE, forExpert);
  TH1F* hESDSumResidualY2PerDE = new TH1F("hESDSumResidualY2PerDE", "sum of cluster-track residual-Y**2 per DE;DetElem ID;#Sigma(#Delta_{Y}^{2}) (cm^{2})", nDE+1, -0.5, nDE+0.5);
  Add2ESDsList(hESDSumResidualY2PerDE, kESDSumResidualY2PerDE, forExpert);
  TH1F* hESDSumLocalChi2XPerDE = new TH1F("hESDSumLocalChi2XPerDE", "sum of local chi2-X per DE;DetElem ID;#Sigma(local #chi^{2}_{X})", nDE+1, -0.5, nDE+0.5);
  Add2ESDsList(hESDSumLocalChi2XPerDE, kESDSumLocalChi2XPerDE, forExpert);
  TH1F* hESDSumLocalChi2YPerDE = new TH1F("hESDSumLocalChi2YPerDE", "sum of local chi2-Y per DE;DetElem ID;#Sigma(local #chi^{2}_{Y})", nDE+1, -0.5, nDE+0.5);
  Add2ESDsList(hESDSumLocalChi2YPerDE, kESDSumLocalChi2YPerDE, forExpert);
  
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
  	
	((AliMUONTrackerDataMaker*)fTrackerDataMaker)->SetRawReader(rawReader);
	
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

      if (darcHeader->GetGlobalFlag()){
        UInt_t nOfClocks = darcHeader->GetGlobalClock();
        Double_t nOfSeconds = ((Double_t) nOfClocks) / 40e6; // 1 clock each 25 ns
        ((TH1F*)GetRawsData(kTriggerScalersTime))->Fill(1., nOfSeconds);
      }

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
	        GetRecPointsData(kTrackerClusterHitMapPerChamber+chamberId)->Fill(cluster->GetX(),cluster->GetY());
		
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
    GetESDsData(kESDProbChi2)->Fill(TMath::Prob(track->GetGlobalChi2(),track->GetNDF()));
    GetESDsData(kESDnClustersPerTrack)->Fill(track->GetNClusters());
    
    // loop over clusters
    AliMUONTrackParam* trackParam = static_cast<AliMUONTrackParam*>(track->GetTrackParamAtCluster()->First());
    while (trackParam) {
      
      AliMUONVCluster* cluster = trackParam->GetClusterPtr();
      Int_t chId = cluster->GetChamberId();
      Int_t deID = cluster->GetDetElemId();
      Double_t residualX = cluster->GetX() - trackParam->GetNonBendingCoor();
      Double_t residualY = cluster->GetY() - trackParam->GetBendingCoor();
      Double_t localChi2X = (cluster->GetErrX2() > 0.) ? residualX*residualX/cluster->GetErrX2() : 0.;
      Double_t localChi2Y = (cluster->GetErrY2() > 0.) ? residualY*residualY/cluster->GetErrY2() : 0.;
      
      GetESDsData(kESDClusterHitMap+chId)->Fill(cluster->GetX(), cluster->GetY());
      
      GetESDsData(kESDnTotClustersPerCh)->Fill(chId);
      GetESDsData(kESDnTotClustersPerDE)->Fill(deID);
      
      GetESDsData(kESDClusterChargeInCh+chId)->Fill(cluster->GetCharge());
      GetESDsData(kESDSumClusterChargePerDE)->Fill(deID, cluster->GetCharge());
      
      if (cluster->GetNDigits() > 0) { // discard clusters with pad not stored in ESD
	GetESDsData(kESDnTotFullClustersPerDE)->Fill(deID);
        GetESDsData(kESDClusterSizeInCh+chId)->Fill(cluster->GetNDigits());
	GetESDsData(kESDSumClusterSizePerDE)->Fill(deID, cluster->GetNDigits());
      }
      
      GetESDsData(kESDResidualXInCh+chId)->Fill(residualX);
      GetESDsData(kESDResidualYInCh+chId)->Fill(residualY);
      GetESDsData(kESDSumResidualXPerDE)->Fill(deID, residualX);
      GetESDsData(kESDSumResidualYPerDE)->Fill(deID, residualY);
      GetESDsData(kESDSumResidualX2PerDE)->Fill(deID, residualX*residualX);
      GetESDsData(kESDSumResidualY2PerDE)->Fill(deID, residualY*residualY);
      
      GetESDsData(kESDLocalChi2XInCh+chId)->Fill(localChi2X);
      GetESDsData(kESDLocalChi2YInCh+chId)->Fill(localChi2Y);
      GetESDsData(kESDSumLocalChi2XPerDE)->Fill(deID, localChi2X);
      GetESDsData(kESDSumLocalChi2YPerDE)->Fill(deID, localChi2Y);
      
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
      if(histoStrips->GetEntries()==0) continue; // No scalers found
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
      if(task==AliQA::kRAWS) {
        Float_t acquisitionTime = ((TH1F*)GetRawsData(kTriggerScalersTime))->GetBinContent(1);
        histoDisplayStrips->Scale(1./acquisitionTime);
      }
    } // iChamber
  } // iCath

  if(task==AliQA::kRAWS){    
    TH2F* histoI  = (TH2F*) GetRawsData(kTriggerRPCi);
    TH2F* histoHV = (TH2F*) GetRawsData(kTriggerRPChv);
    FillTriggerDCSHistos();
    for (Int_t iChamber = 0; iChamber < AliMpConstants::NofTriggerChambers(); iChamber++) {
      Int_t bin = histoI->GetXaxis()->FindBin(11+iChamber);
      TH2F* histoDisplayI = (TH2F*)GetRawsData(kTriggerIDisplay + iChamber);
      triggerDisplay.FillDisplayHistogram(histoI->ProjectionY("_px", bin, bin), histoDisplayI, AliMUONTriggerDisplay::kDisplaySlats, 0, iChamber);
      TH2F* histoDisplayHV = (TH2F*)GetRawsData(kTriggerHVDisplay + iChamber);
      bin = histoHV->GetXaxis()->FindBin(11+iChamber);
      triggerDisplay.FillDisplayHistogram(histoHV->ProjectionY("_px", bin, bin), histoDisplayHV, AliMUONTriggerDisplay::kDisplaySlats, 0, iChamber);
    }
  }

  if(task==AliQA::kRECPOINTS){
    TH1F* histoBoards = (TH1F*)GetRecPointsData(kTriggeredBoards);
    TH2F* histoDisplayBoards = (TH2F*)GetRecPointsData(kTriggerBoardsDisplay);
    triggerDisplay.FillDisplayHistogram(histoBoards, histoDisplayBoards, AliMUONTriggerDisplay::kDisplayBoards, 0, 0);
  }
}


//_____________________________________________________________________________
Bool_t 
AliMUONQADataMakerRec::FillTriggerDCSHistos()
{
  /// Get HV and currents values for one trigger chamber
  
  AliCodeTimerAuto("");

  AliMUONCalibrationData calibrationData(AliCDBManager::Instance()->GetRun());

  TMap* triggerDcsMap = calibrationData.TriggerDCS();

  if ( !triggerDcsMap ) 
  {
    AliError("Cannot fill DCS histos, as triggerDcsMap is NULL");
    return kFALSE;
  }

  AliMpDEIterator deIt;
  
  deIt.First();
    
  AliMpDCSNamer triggerDcsNamer("TRIGGER");

  TH2* currHisto = 0x0;

  Bool_t error = kFALSE;
  
  while ( !deIt.IsDone() )
  {
    Int_t detElemId = deIt.CurrentDEId();
    
    if ( AliMpDEManager::GetStationType(detElemId) == AliMp::kStationTrigger) {

      Int_t iChamber = AliMpDEManager::GetChamberId(detElemId);
      Int_t slat = detElemId%100;

      for(Int_t iMeas=0; iMeas<AliMpDCSNamer::kNDCSMeas; iMeas++){
	TString currAlias = triggerDcsNamer.DCSChannelName(detElemId, 0, iMeas);

	AliDebug(2, Form("\nDetElemId %i   dcsAlias %s", detElemId, currAlias.Data()));

	TPair* triggerDcsPair = static_cast<TPair*>(triggerDcsMap->FindObject(currAlias.Data()));

	if (!triggerDcsPair)
	{
	  AliError(Form("Did not find expected alias (%s) for DE %d",
			currAlias.Data(),detElemId));  
	  error = kTRUE;
	}
	else
	{
	  TObjArray* values = static_cast<TObjArray*>(triggerDcsPair->Value());
	  if (!values)
	  {
	    AliError(Form("Could not get values for alias %s",currAlias.Data()));
	    error = kTRUE;
	  }
	  else
	  {
	    TIter next(values);
	    AliDCSValue* val;

	    while ( ( val = static_cast<AliDCSValue*>(next()) ) )
	    {
	      Float_t hvi = val->GetFloat();

	      AliDebug(2, Form("Value %f", hvi));

	      switch(iMeas){
	      case AliMpDCSNamer::kDCSI:
		currHisto = (TH2F*) GetRawsData(kTriggerRPCi);
		break;
	      case AliMpDCSNamer::kDCSHV:
		currHisto = (TH2F*) GetRawsData(kTriggerRPChv);
		break;
	      } 
	      Int_t binX = currHisto->GetXaxis()->FindBin(iChamber+1);
	      Int_t binY = currHisto->GetYaxis()->FindBin(slat);
	      currHisto->SetBinContent(binX, binY, hvi);
	    } // loop on values
	  } // if (!values)
	} // if (!triggerDcsPair)
      } // loop on measured types (HV and currents)
    } // if (stationType == kStationTrigger)

    deIt.Next();
  }
  return error;
}

//____________________________________________________________________________ 
AliMUONVTrackerData* AliMUONQADataMakerRec::GetTrackerData() const
{ 
  
  return fTrackerDataMaker->Data(); 
  
}
