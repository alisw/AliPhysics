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

// $Id: AliMUONTrackerQADataMakerRec.cxx 35760 2009-10-21 21:45:42Z ivana $

// --- MUON header files ---
#include "AliMUONTrackerQADataMakerRec.h"

#include "AliDAQ.h"
#include "AliQAv1.h"
#include "AliMpDDL.h"
#include "AliMUONConstants.h"  
#include "AliMUONDigitMaker.h"
#include "AliMUONQAIndices.h"
#include "AliMUONQAMappingCheck.h"
#include "AliMUONLogger.h"
#include "AliMUONTrackerDataMaker.h"
#include "AliMUONVCluster.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVTrackerData.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONESDInterface.h"
#include "AliMUONCalibrationData.h"
#include "AliMpBusPatch.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"

// --- AliRoot header files ---
#include "AliCDBManager.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliLog.h"
#include "AliRawReader.h"
#include "AliCodeTimer.h"
#include "AliMUONVDigit.h"

// --- ROOT system ---
#include <TH1F.h> 
#include <TH1I.h> 
#include <TH2F.h>
#include <Riostream.h>
#include <TMath.h>
#include <TPaveText.h>

//-----------------------------------------------------------------------------
/// \class AliMUONTrackerQADataMakerRec
///
/// MUON base class for quality assurance data (histo) maker
///
/// \author C. Finck, D. Stocco, L. Aphecetche

/// \cond CLASSIMP
ClassImp(AliMUONTrackerQADataMakerRec)
/// \endcond
           
namespace
{
  Double_t ProtectedSqrt(Double_t x)
  {
    return ( x > 0.0 ? TMath::Sqrt(x) : 0.0 );
  }
}
//____________________________________________________________________________ 
AliMUONTrackerQADataMakerRec::AliMUONTrackerQADataMakerRec(AliQADataMakerRec* master) : 
AliMUONVQADataMakerRec(master),
fDigitStore(AliMUONVDigitStore::Create("AliMUONDigitStoreV1")),
fDigitMaker(new AliMUONDigitMaker(kTRUE)),
fClusterStore(0x0),
fTrackerDataMaker(0x0),
fMappingCheckRecPoints(0x0),
fCalibrationData(new AliMUONCalibrationData(AliCDBManager::Instance()->GetRun())),
fLogger(0x0)
{
  /// ctor
}

//__________________________________________________________________
AliMUONTrackerQADataMakerRec::~AliMUONTrackerQADataMakerRec()
{
  /// dtor
  delete fDigitStore;
  delete fDigitMaker;
  delete fClusterStore;
  delete fCalibrationData;
  delete fMappingCheckRecPoints;
  if (fLogger)
  {
    if ( fLogger->NumberOfEntries() != 0 ) 
    {
      AliError("Some unprocessed logged errors are still there... Please check");
    }
    delete fLogger;    
  }
}

//____________________________________________________________________________ 
void AliMUONTrackerQADataMakerRec::InsertTrackerData(Int_t specie, TObjArray** list,
                                                     TObject* object, Int_t indexNumber,
                                                     Bool_t replace)
{
  /// Insert an object to a given list
  
  TIter next(list[specie]);
  TObject* o;
  TObject* old(0x0);
  Bool_t alreadyThere(kFALSE);
  while ( ( o = next() ) && !alreadyThere )
  {
    TString classname(o->ClassName());
    if ( classname.Contains("TrackerData") ) 
    {
      alreadyThere = kTRUE;
      old = o;
    }
  }
  if ( (!alreadyThere && object) || (alreadyThere && replace) )
  {
    delete old;
    AliDebug(AliQAv1::GetQADebugLevel(), Form("Adding %s to the list of qa objects",object->GetName()));
    TNamed* named = static_cast<TNamed*>(object);
    named->SetName(Form("%s_%s",AliRecoParam::GetEventSpecieName(specie),object->GetName()));
    object->SetBit(AliQAv1::GetExpertBit());
    list[specie]->AddAt(object,indexNumber);
  }
}

//____________________________________________________________________________ 
void AliMUONTrackerQADataMakerRec::EndOfDetectorCycleESDs(Int_t, TObjArray**)
{
  /// Normalize ESD histograms
  
  if (!GetESDsData(AliMUONQAIndices::kESDnClustersPerTrack)) return;
  
  Double_t nTracks = GetESDsData(AliMUONQAIndices::kESDnClustersPerTrack)->GetEntries();
  if (nTracks <= 0) return;
  
  TH1* hESDnClustersPerCh = GetESDsData(AliMUONQAIndices::kESDnClustersPerCh);
  TH1* hESDnClustersPerDE = GetESDsData(AliMUONQAIndices::kESDnClustersPerDE);
  TH1* hESDClusterChargePerChMean = GetESDsData(AliMUONQAIndices::kESDClusterChargePerChMean);
  TH1* hESDClusterChargePerChSigma = GetESDsData(AliMUONQAIndices::kESDClusterChargePerChSigma);
  TH1* hESDClusterSizePerChMean = GetESDsData(AliMUONQAIndices::kESDClusterSizePerChMean);
  TH1* hESDClusterSizePerChSigma = GetESDsData(AliMUONQAIndices::kESDClusterSizePerChSigma);
  TH1* hESDResidualXPerChMean = GetESDsData(AliMUONQAIndices::kESDResidualXPerChMean);
  TH1* hESDResidualXPerChSigma = GetESDsData(AliMUONQAIndices::kESDResidualXPerChSigma);
  TH1* hESDResidualYPerChMean = GetESDsData(AliMUONQAIndices::kESDResidualYPerChMean);
  TH1* hESDResidualYPerChSigma = GetESDsData(AliMUONQAIndices::kESDResidualYPerChSigma);
  TH1* hESDLocalChi2XPerChMean = GetESDsData(AliMUONQAIndices::kESDLocalChi2XPerChMean);
  TH1* hESDLocalChi2YPerChMean = GetESDsData(AliMUONQAIndices::kESDLocalChi2YPerChMean);
  TH1* hESDLocalChi2PerChMean = GetESDsData(AliMUONQAIndices::kESDLocalChi2PerChMean);
  TH1* hESDClusterChargePerDE = GetESDsData(AliMUONQAIndices::kESDClusterChargePerDE);
  TH1* hESDClusterSizePerDE = GetESDsData(AliMUONQAIndices::kESDClusterSizePerDE);
  TH1* hESDResidualXPerDEMean = GetESDsData(AliMUONQAIndices::kESDResidualXPerDEMean);
  TH1* hESDResidualXPerDESigma = GetESDsData(AliMUONQAIndices::kESDResidualXPerDESigma);
  TH1* hESDResidualYPerDEMean = GetESDsData(AliMUONQAIndices::kESDResidualYPerDEMean);
  TH1* hESDResidualYPerDESigma = GetESDsData(AliMUONQAIndices::kESDResidualYPerDESigma);
  TH1* hESDLocalChi2XPerDEMean = GetESDsData(AliMUONQAIndices::kESDLocalChi2XPerDEMean);
  TH1* hESDLocalChi2YPerDEMean = GetESDsData(AliMUONQAIndices::kESDLocalChi2YPerDEMean);
  TH1* hESDLocalChi2PerDEMean = GetESDsData(AliMUONQAIndices::kESDLocalChi2PerDEMean);
  TH1* hESDnTotClustersPerCh = GetESDsData(AliMUONQAIndices::kESDnTotClustersPerCh);
  TH1* hESDnTotClustersPerDE = GetESDsData(AliMUONQAIndices::kESDnTotClustersPerDE);
  TH1* hESDnTotFullClustersPerDE = GetESDsData(AliMUONQAIndices::kESDnTotFullClustersPerDE);
  TH1* hESDSumClusterChargePerDE = GetESDsData(AliMUONQAIndices::kESDSumClusterChargePerDE);
  TH1* hESDSumClusterSizePerDE = GetESDsData(AliMUONQAIndices::kESDSumClusterSizePerDE);
  TH1* hESDSumResidualXPerDE = GetESDsData(AliMUONQAIndices::kESDSumResidualXPerDE);
  TH1* hESDSumResidualYPerDE = GetESDsData(AliMUONQAIndices::kESDSumResidualYPerDE);
  TH1* hESDSumResidualX2PerDE = GetESDsData(AliMUONQAIndices::kESDSumResidualX2PerDE);
  TH1* hESDSumResidualY2PerDE = GetESDsData(AliMUONQAIndices::kESDSumResidualY2PerDE);
  TH1* hESDSumLocalChi2XPerDE = GetESDsData(AliMUONQAIndices::kESDSumLocalChi2XPerDE);
  TH1* hESDSumLocalChi2YPerDE = GetESDsData(AliMUONQAIndices::kESDSumLocalChi2YPerDE);
  TH1* hESDSumLocalChi2PerDE = GetESDsData(AliMUONQAIndices::kESDSumLocalChi2PerDE);
  
  hESDnClustersPerCh->Reset();
  hESDnClustersPerDE->Reset();
  hESDnClustersPerCh->Add(hESDnTotClustersPerCh, 1./nTracks);
  hESDnClustersPerDE->Add(hESDnTotClustersPerDE, 1./nTracks);
  
  // loop over chambers
  for (Int_t iCh = 0; iCh < AliMUONConstants::NTrackingCh(); iCh++) {
    
    TH1* hESDClusterChargeInCh = GetESDsData(AliMUONQAIndices::kESDClusterChargeInCh+iCh);
    Double_t sigmaCharge = hESDClusterChargeInCh->GetRMS();
    hESDClusterChargePerChMean->SetBinContent(iCh+1, hESDClusterChargeInCh->GetMean());
    hESDClusterChargePerChMean->SetBinError(iCh+1, hESDClusterChargeInCh->GetMeanError());
    hESDClusterChargePerChSigma->SetBinContent(iCh+1, sigmaCharge);
    hESDClusterChargePerChSigma->SetBinError(iCh+1, hESDClusterChargeInCh->GetRMSError());
    
    TH1* hESDClusterSizeInCh = GetESDsData(AliMUONQAIndices::kESDClusterSizeInCh+iCh);
    Double_t sigmaSize = hESDClusterSizeInCh->GetRMS();
    hESDClusterSizePerChMean->SetBinContent(iCh+1, hESDClusterSizeInCh->GetMean());
    hESDClusterSizePerChMean->SetBinError(iCh+1, hESDClusterSizeInCh->GetMeanError());
    hESDClusterSizePerChSigma->SetBinContent(iCh+1, sigmaSize);
    hESDClusterSizePerChSigma->SetBinError(iCh+1, hESDClusterSizeInCh->GetRMSError());
    
    TH1* hESDResidualXInCh = GetESDsData(AliMUONQAIndices::kESDResidualXInCh+iCh);
    Double_t sigmaResidualX = hESDResidualXInCh->GetRMS();
    hESDResidualXPerChMean->SetBinContent(iCh+1, hESDResidualXInCh->GetMean());
    hESDResidualXPerChMean->SetBinError(iCh+1, hESDResidualXInCh->GetMeanError());
    hESDResidualXPerChSigma->SetBinContent(iCh+1, sigmaResidualX);
    hESDResidualXPerChSigma->SetBinError(iCh+1, hESDResidualXInCh->GetRMSError());
    
    TH1* hESDResidualYInCh = GetESDsData(AliMUONQAIndices::kESDResidualYInCh+iCh);
    Double_t sigmaResidualY = hESDResidualYInCh->GetRMS();
    hESDResidualYPerChMean->SetBinContent(iCh+1, hESDResidualYInCh->GetMean());
    hESDResidualYPerChMean->SetBinError(iCh+1, hESDResidualYInCh->GetMeanError());
    hESDResidualYPerChSigma->SetBinContent(iCh+1, sigmaResidualY);
    hESDResidualYPerChSigma->SetBinError(iCh+1, hESDResidualYInCh->GetRMSError());
    
    TH1* hESDLocalChi2XInCh = GetESDsData(AliMUONQAIndices::kESDLocalChi2XInCh+iCh);
    Double_t sigmaLocalChi2X = hESDLocalChi2XInCh->GetRMS();
    hESDLocalChi2XPerChMean->SetBinContent(iCh+1, hESDLocalChi2XInCh->GetMean());
    hESDLocalChi2XPerChMean->SetBinError(iCh+1, hESDLocalChi2XInCh->GetMeanError());
    
    TH1* hESDLocalChi2YInCh = GetESDsData(AliMUONQAIndices::kESDLocalChi2YInCh+iCh);
    Double_t sigmaLocalChi2Y = hESDLocalChi2YInCh->GetRMS();
    hESDLocalChi2YPerChMean->SetBinContent(iCh+1, hESDLocalChi2YInCh->GetMean());
    hESDLocalChi2YPerChMean->SetBinError(iCh+1, hESDLocalChi2YInCh->GetMeanError());
    
    TH1* hESDLocalChi2InCh = GetESDsData(AliMUONQAIndices::kESDLocalChi2InCh+iCh);
    Double_t sigmaLocalChi2 = hESDLocalChi2InCh->GetRMS();
    hESDLocalChi2PerChMean->SetBinContent(iCh+1, hESDLocalChi2InCh->GetMean());
    hESDLocalChi2PerChMean->SetBinError(iCh+1, hESDLocalChi2InCh->GetMeanError());
    
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
        hESDResidualXPerDESigma->SetBinContent(iDE+1, ProtectedSqrt(hESDSumResidualX2PerDE->GetBinContent(iDE+1)/nClusters - meanResX*meanResX));
        hESDResidualXPerDESigma->SetBinError(iDE+1, sigmaResidualX/TMath::Sqrt(2.*nClusters));
        
        Double_t meanResY = hESDSumResidualYPerDE->GetBinContent(iDE+1)/nClusters;
        hESDResidualYPerDEMean->SetBinContent(iDE+1, meanResY);
        hESDResidualYPerDEMean->SetBinError(iDE+1, sigmaResidualY/TMath::Sqrt(nClusters));
        hESDResidualYPerDESigma->SetBinContent(iDE+1, ProtectedSqrt(hESDSumResidualY2PerDE->GetBinContent(iDE+1)/nClusters - meanResY*meanResY));
        hESDResidualYPerDESigma->SetBinError(iDE+1, sigmaResidualY/TMath::Sqrt(2.*nClusters));
        
        hESDLocalChi2XPerDEMean->SetBinContent(iDE+1, hESDSumLocalChi2XPerDE->GetBinContent(iDE+1)/nClusters);
        hESDLocalChi2XPerDEMean->SetBinError(iDE+1, sigmaLocalChi2X/TMath::Sqrt(nClusters));
        
        hESDLocalChi2YPerDEMean->SetBinContent(iDE+1, hESDSumLocalChi2YPerDE->GetBinContent(iDE+1)/nClusters);
        hESDLocalChi2YPerDEMean->SetBinError(iDE+1, sigmaLocalChi2Y/TMath::Sqrt(nClusters));
        
        hESDLocalChi2PerDEMean->SetBinContent(iDE+1, hESDSumLocalChi2PerDE->GetBinContent(iDE+1)/nClusters);
        hESDLocalChi2PerDEMean->SetBinError(iDE+1, sigmaLocalChi2/TMath::Sqrt(nClusters));
        
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
        
        hESDLocalChi2PerDEMean->SetBinContent(iDE+1, hESDSumLocalChi2PerDE->GetBinContent(iDE+1));
        hESDLocalChi2PerDEMean->SetBinError(iDE+1, hESDLocalChi2InCh->GetXaxis()->GetXmax());
        
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
  
//____________________________________________________________________________ 
void AliMUONTrackerQADataMakerRec::EndOfDetectorCycleRecPoints(Int_t specie, TObjArray** list)
{
  /// Normalize RecPoints histograms
  
  if (!GetRecPointsData(AliMUONQAIndices::kTrackerClusterChargePerChMean)) return;
  
  TH1* hTrackerClusterChargePerChMean = GetRecPointsData(AliMUONQAIndices::kTrackerClusterChargePerChMean);
  TH1* hTrackerClusterChargePerChSigma = GetRecPointsData(AliMUONQAIndices::kTrackerClusterChargePerChSigma);
  TH1* hTrackerClusterMultiplicityPerChMean = GetRecPointsData(AliMUONQAIndices::kTrackerClusterMultiplicityPerChMean);
  TH1* hTrackerClusterMultiplicityPerChSigma = GetRecPointsData(AliMUONQAIndices::kTrackerClusterMultiplicityPerChSigma);
  TH1* hTrackerClusterChargePerDEMean = GetRecPointsData(AliMUONQAIndices::kTrackerClusterChargePerDEMean);
  TH1* hTrackerClusterMultiplicityPerDEMean = GetRecPointsData(AliMUONQAIndices::kTrackerClusterMultiplicityPerDEMean);
  
  // loop over chambers
  for (Int_t iCh = 0; iCh < AliMUONConstants::NTrackingCh(); iCh++) {
    
    TH1* hTrackerClusterChargePerChamber = GetRecPointsData(AliMUONQAIndices::kTrackerClusterChargePerChamber+iCh);
    Double_t sigmaCharge = hTrackerClusterChargePerChamber->GetRMS();
    hTrackerClusterChargePerChMean->SetBinContent(iCh+1, hTrackerClusterChargePerChamber->GetMean());
    hTrackerClusterChargePerChMean->SetBinError(iCh+1, hTrackerClusterChargePerChamber->GetMeanError());
    hTrackerClusterChargePerChSigma->SetBinContent(iCh+1, sigmaCharge);
    hTrackerClusterChargePerChSigma->SetBinError(iCh+1, hTrackerClusterChargePerChamber->GetRMSError());
    
    TH1* hTrackerClusterMultiplicityPerChamber = GetRecPointsData(AliMUONQAIndices::kTrackerClusterMultiplicityPerChamber+iCh);
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
      
      TH1* hTrackerClusterChargePerDE = GetRecPointsData(AliMUONQAIndices::kTrackerClusterChargePerDE+iDE);
      hTrackerClusterChargePerDEMean->SetBinContent(iDE+1, hTrackerClusterChargePerDE->GetMean());
      Double_t nClusters = hTrackerClusterChargePerDE->GetEntries();
      if (nClusters > 1) hTrackerClusterChargePerDEMean->SetBinError(iDE+1, sigmaCharge/TMath::Sqrt(nClusters));
      else hTrackerClusterChargePerDEMean->SetBinError(iDE+1, hTrackerClusterChargePerChamber->GetXaxis()->GetXmax());
      
      TH1* hTrackerClusterMultiplicityPerDE = GetRecPointsData(AliMUONQAIndices::kTrackerClusterMultiplicityPerDE+iDE);
      hTrackerClusterMultiplicityPerDEMean->SetBinContent(iDE+1, hTrackerClusterMultiplicityPerDE->GetMean());
      nClusters = hTrackerClusterMultiplicityPerDE->GetEntries();
      if (nClusters > 1) hTrackerClusterMultiplicityPerDEMean->SetBinError(iDE+1, sigmaSize/TMath::Sqrt(nClusters));
      else hTrackerClusterMultiplicityPerDEMean->SetBinError(iDE+1, hTrackerClusterMultiplicityPerChamber->GetXaxis()->GetXmax());
      
      it.Next();
    }
  }
  
  if ( fMappingCheckRecPoints ) InsertTrackerData(specie,list,fMappingCheckRecPoints->CreateData("RecPoints"),AliMUONQAIndices::kTrackerRecPoints,kTRUE);
}


//____________________________________________________________________________ 
void AliMUONTrackerQADataMakerRec::EndOfDetectorCycleRaws(Int_t specie, TObjArray** list)
{
  /// create Raws histograms in Raws subdir
  
  if ( !GetRawsData(AliMUONQAIndices::kTrackerBusPatchOccupancy) ) return;

  if ( fTrackerDataMaker ) 
  {
    /// put the trackerdata in the pipeline
    InsertTrackerData(specie,list,fTrackerDataMaker->Data(),AliMUONQAIndices::kTrackerData);

    /// project the tracerdata into buspatch occupancies (for the experts)
    TH1* hbp = GetRawsData(AliMUONQAIndices::kTrackerBusPatchOccupancy);
    hbp->Reset();
    TIter nextBP(AliMpDDLStore::Instance()->CreateBusPatchIterator());
    AliMpBusPatch* bp(0x0);
    AliMUONVTrackerData* data = fTrackerDataMaker->Data();
    Int_t occDim = 2;
    
    while ( ( bp = static_cast<AliMpBusPatch*>(nextBP())) )
    {
      Int_t busPatchId = bp->GetId();
      Int_t bin = hbp->FindBin(busPatchId);
      hbp->SetBinContent(bin,data->BusPatch(busPatchId,occDim)*100.0); // occupancy, in percent
    }
    
    /// log the readout errors (for the shifter)
    TH1* hnevents = GetRawsData(AliMUONQAIndices::kTrackerNofRawEventSeen);
    hnevents->Reset();
    hnevents->Fill(0.0,fTrackerDataMaker->Data()->NumberOfEvents(-1));
    
    if ( fLogger->NumberOfEntries() > 0 )
    {
      // readout errors      
      FillErrors(*fLogger);      
      fLogger->Clear();
    }

    /// project tracker data into DDL occupancies (for the shifter)
    TH1* hddl = GetRawsData(AliMUONQAIndices::kTrackerDDLOccupancy);
    hddl->Reset();
    TH1* hddlevents = GetRawsData(AliMUONQAIndices::kTrackerDDLNofEvents);
    hddlevents->Reset();
    
    const Int_t nddls = AliDAQ::NumberOfDdls("MUONTRK");
    const Int_t offset = AliDAQ::DdlIDOffset("MUONTRK");
    
    for ( Int_t iddl = 0; iddl < nddls; ++iddl )
    {
      AliMpDDL* ddl = AliMpDDLStore::Instance()->GetDDL(iddl);
      
      Int_t ddlId = offset + ddl->GetId();
      Int_t npads = 0;
      
      Int_t nevents = data->NumberOfEvents(iddl);
      
      hddlevents->Fill(ddlId,nevents);
      
      Double_t occ(0.0);
      
      for ( Int_t ide = 0; ide < ddl->GetNofDEs(); ++ide )
      {
        Int_t de = ddl->GetDEId(ide);
        
        npads += TMath::Nint(data->DetectionElement(de,3));
        
        occ +=  data->DetectionElement(de,4);
      }
      
      if ( nevents > 0 && npads > 0 )
      {
        occ = occ/npads/nevents;
      }

      hddl->Fill(ddlId,100.0*occ); // occ in percent
    }
  }
}

//____________________________________________________________________________ 
void AliMUONTrackerQADataMakerRec::FillErrors(AliMUONLogger& log)
{
  log.ResetItr();

  TString msg;
  Int_t occurence;
  
  TH1* hparity = GetRawsData(AliMUONQAIndices::kTrackerBusPatchParityErrors);

  TH1* htoken = GetRawsData(AliMUONQAIndices::kTrackerBusPatchTokenLostErrors);

  TH1* hpadding = GetRawsData(AliMUONQAIndices::kTrackerBusPatchPaddingErrors);
  
  TH1* hroe = GetRawsData(AliMUONQAIndices::kTrackerReadoutErrors);
    
  TH1* hnevents = GetRawsData(AliMUONQAIndices::kTrackerNofRawEventSeen);
  
  Int_t nevents = TMath::Nint(hnevents->GetBinContent(1));
  
  if ( !nevents ) 
  {
    TPaveText* text = new TPaveText(0,0,0.99,0.99,"NDC");
    text->AddText("FATAL : 0 event seen ? That's NOT normal...");
    text->SetFillColor(2); // red = FATAL
    hroe->GetListOfFunctions()->Add(text);
    return;
  }
  
  while ( log.Next(msg,occurence) )
  {
    AliDebug(1,Form("msg=%s occurence=%d",msg.Data(),occurence));
             
    if ( msg.Contains("token") )
    {
      Int_t dsp(-1),iddl(-1),ecode(-1);
      
      sscanf(msg.Data(),"Lost token error detected with address 0x%X of DDL %d and code %d.",
             &dsp,&iddl,&ecode);
      Int_t localBP = ((dsp >> 16)- 4)*5 + 1;
      Int_t buspatch = localBP + iddl*100;
      
      // Let's try to get all the suspected bus patches (i.e. one full FRT, as currently
      // we don't have more precise information to locate the faulty bus patch(es)).
      
      AliMpBusPatch* bp = AliMpDDLStore::Instance()->GetBusPatch(buspatch);
      Int_t frt = bp->GetFrtId();
      AliMpDDL* ddl = AliMpDDLStore::Instance()->GetDDL(bp->GetDdlId());
      Int_t* b = new Int_t[ddl->GetMaxDsp()];
      ddl->GetBusPerDsp(b);
      Int_t nbus(0);
      for ( Int_t i = 0; i < ddl->GetNofFrts() && !nbus; ++i ) 
      {
        if ( ddl->GetFrtId(i) == frt ) 
        {
          nbus = b[i];
        }
      }
      if (nbus<=0) 
      {
        AliError("GOT NBUS<=0 ! THAT IS BAD ! CHECK !");
        nbus=1;
      }
      
      delete[] b;
      
      while (nbus) {
        htoken->Fill(buspatch+nbus-1,occurence);
        --nbus;
      }
      
      hroe->Fill(1.0*AliMUONQAIndices::kTrackerRawNofTokenLostErrors,occurence);
    }
    
    if ( msg.Contains("Parity") )
    {
      Int_t buspatch;
      sscanf(msg.Data(),"Parity error in buspatch %d (0x%X).",&buspatch,&buspatch);
      hparity->Fill(buspatch,occurence);      
      hroe->Fill(1.0*AliMUONQAIndices::kTrackerRawNofParityErrors,occurence);
    }
    
    if ( msg.Contains("Glitch") ) 
    {
      hroe->Fill(1.0*AliMUONQAIndices::kTrackerRawNofGlitchErrors,occurence);      
    }
    
    if ( msg.Contains("Padding") )
    {
      Int_t block, dsp, buspatch;      
      sscanf(msg.Data(),"Padding word error for iBlock %d, iDsp %d, iBus %d.",&block,&dsp,&buspatch);
      hpadding->Fill(buspatch,occurence);
      hroe->Fill(1.0*AliMUONQAIndices::kTrackerRawNofPaddingErrors,occurence);      
    }
  }
  
  /// Finally we normalize to the number of events, for the shifter plot only
  hroe->Scale(100.0/nevents);
}


//____________________________________________________________________________ 
void AliMUONTrackerQADataMakerRec::InitRaws()
{
    /// create Raws histograms in Raws subdir
	
  AliCodeTimerAuto("",0);
  
  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 
 
  Int_t bpmin(999999);
  Int_t bpmax(0);
  
  TIter next(AliMpDDLStore::Instance()->CreateBusPatchIterator());
  AliMpBusPatch* bp(0x0);
  while ( ( bp = static_cast<AliMpBusPatch*>(next())) )
  {
    bpmin = TMath::Min(bpmin,bp->GetId());
    bpmax = TMath::Max(bpmax,bp->GetId());
  }
  
  Double_t xmin = bpmin-0.5;
  Double_t xmax = bpmax+0.5;
  Int_t nbins = bpmax-bpmin+1;
  
  TH1* hbp = new TH1F("hTrackerBusPatchOccupancy","Occupancy of bus patches",nbins,xmin,xmax);

  TH1* hbpnpads = new TH1F("hTrackerBusPatchNofPads","Number of pads per bus patch",nbins,xmin,xmax);

  TH1* hbpnmanus = new TH1F("hTrackerBusPatchNofManus","Number of manus per bus patch",nbins,xmin,xmax);

  Add2RawsList(hbp,AliMUONQAIndices::kTrackerBusPatchOccupancy, expert, image, !saveCorr);
  Add2RawsList(hbpnpads,AliMUONQAIndices::kTrackerBusPatchNofPads, expert, !image, !saveCorr);
  Add2RawsList(hbpnmanus,AliMUONQAIndices::kTrackerBusPatchNofManus, expert, !image, !saveCorr);

  Add2RawsList(new TH1F("hTrackerBusPatchParityErrors","Number of parity errors per bus patch",nbins,xmin,xmax),
               AliMUONQAIndices::kTrackerBusPatchParityErrors,expert,!image,!saveCorr);

  Add2RawsList(new TH1F("hTrackerBusPatchTokenLostErrors","Number of token lost errors per bus patch",nbins,xmin,xmax),
               AliMUONQAIndices::kTrackerBusPatchTokenLostErrors,expert,!image,!saveCorr);

  Add2RawsList(new TH1F("hTrackerBusPatchPaddingErrors","Number of padding errors per bus patch",nbins,xmin,xmax),
               AliMUONQAIndices::kTrackerBusPatchPaddingErrors,expert,!image,!saveCorr);

  TH1* h = new TH1F("hTrackerReadoutErrors","Number of readout errors per event;;Error rate in %",4,-0.5,3.5);
  h->SetStats(kFALSE);
  
  // The QA shifter will only see the summary plot below
  TAxis* a = h->GetXaxis();
  
  a->SetBinLabel(h->FindBin(1.0*AliMUONQAIndices::kTrackerRawNofGlitchErrors),"Glitch errors");
  a->SetBinLabel(h->FindBin(1.0*AliMUONQAIndices::kTrackerRawNofTokenLostErrors),"Token lost errors");
  a->SetBinLabel(h->FindBin(1.0*AliMUONQAIndices::kTrackerRawNofParityErrors),"Parity errors");
  a->SetBinLabel(h->FindBin(1.0*AliMUONQAIndices::kTrackerRawNofPaddingErrors),"Padding errors");

  Add2RawsList(h,AliMUONQAIndices::kTrackerReadoutErrors,!expert,image,!saveCorr);

  TH1* hnevents = new TH1F("hTrackerNofRawEventSeen","Number of events seen",1,-0.5,0.5);
  a = hnevents->GetXaxis();
  a->SetBinLabel(1,"Nevents");
  hnevents->SetStats(kFALSE);
  
  Add2RawsList(hnevents,AliMUONQAIndices::kTrackerNofRawEventSeen,expert,!image,!saveCorr);

  const Bool_t histogram(kFALSE);

  if(!fTrackerDataMaker) 
  {
    
    AliMUONTrackerDataMaker* dm = new AliMUONTrackerDataMaker(GetRecoParam(),
                                                    AliCDBManager::Instance()->GetRun(),
                                                    0x0,
                                                    "",
                                                    "NOGAIN",
                                                    histogram,
                                                    0.0,0.0);

    fLogger = new AliMUONLogger(-1);
    dm->EnableErrorLogger(fLogger);
    fTrackerDataMaker = dm;
  }
  
  fTrackerDataMaker->Data()->DisableChannelLevel(); // to save up disk space, we only store starting at the manu level

  fTrackerDataMaker->SetRunning(kTRUE);
  
  next.Reset();

  AliMUONVStore* config = fCalibrationData->Config();

  TH1* hbpconfig(0x0);
  
  if (config)
  {
    hbpconfig = new TH1F("hTrackerBusPatchConfig","Configuration of bus patches",nbins,xmin,xmax);
    Add2RawsList(hbpconfig,AliMUONQAIndices::kTrackerBusPatchConfig, expert, !image, !saveCorr);
  }
  else
  {
    AliWarning("Tracker configuration not found. Will not be able to cut on low occupancies");
  }
  
  while ( ( bp = static_cast<AliMpBusPatch*>(next())) )
  {
    Int_t n(0);
    Bool_t inConfig(kTRUE);
    
    for ( Int_t imanu = 0; imanu < bp->GetNofManus(); ++imanu )
    {
      Int_t manuId = bp->GetManuId(imanu);
      AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(bp->GetDEId());      
      n += de->NofChannelsInManu(manuId);
      if ( config && !config->FindObject(de->GetId(),manuId)) inConfig=kFALSE;
    }
    hbpnpads->Fill(bp->GetId(),n*1.0);
    hbpnmanus->Fill(bp->GetId(),bp->GetNofManus()*1.0);
    if ( hbpconfig && inConfig ) 
    {
      hbpconfig->Fill(bp->GetId());
    }
  }
  
  nbins = AliDAQ::NumberOfDdls("MUONTRK");
  const Int_t offset = AliDAQ::DdlIDOffset("MUONTRK");
  
  xmin = offset - 0.5;
  xmax  = offset + nbins - 0.5;
  
  Add2RawsList(new TH1F("hTrackerDDLOccupancy",";DDLId;DDL Occupancy in %",nbins,xmin,xmax),
               AliMUONQAIndices::kTrackerDDLOccupancy,!expert,image,!saveCorr);
  Add2RawsList(new TH1F("hTrackerDDLNofEvents","Number of events seen by DDL;DDLId",nbins,xmin,xmax),
               AliMUONQAIndices::kTrackerDDLNofEvents,expert,!image,!saveCorr);
  
}

//__________________________________________________________________
void AliMUONTrackerQADataMakerRec::InitDigits() 
{
  /// Initialized Digits spectra 
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1I* h0 = new TH1I("hDigitsDetElem", "Detection element distribution in Digits;Detection element Id;Counts",  1400, 100, 1500); 
  Add2DigitsList(h0, 0, !expert, image);
  
  TH1I* h1 = new TH1I("hDigitsADC", "ADC distribution in Digits;ACD value;Counts", 4096, 0, 4095); 
  Add2DigitsList(h1, 1, !expert, image);    
} 

//____________________________________________________________________________ 
void AliMUONTrackerQADataMakerRec::InitRecPoints()
{
  /// create Reconstructed Points histograms in RecPoints subdir for the
  /// MUON tracker subsystem.
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  AliCodeTimerAuto("",0);
  
  TH1I *h1I;
  TH1F *h1F;
  TH2F *h2F;
  
  // histograms per chamber
  Int_t nCh = AliMpConstants::NofTrackingChambers();
  for ( Int_t i = 0; i < nCh; ++i ) 
  {
    h1I = new TH1I(Form("hTrackerClusterMultiplicityForChamber%d",i+1), Form("cluster size distribution in chamber %d;size (n_{pads};Counts)",i+1), 100,0,100);
    Add2RecPointsList(h1I,AliMUONQAIndices::kTrackerClusterMultiplicityPerChamber+i, expert, !image);
    
    h1I = new TH1I(Form("hTrackerClusterChargeForChamber%d",i+1), Form("cluster charge distribution in chamber %d;charge (fC);Counts",i+1), 100,0,1000);
    Add2RecPointsList(h1I,AliMUONQAIndices::kTrackerClusterChargePerChamber+i, expert, !image);
    
    Float_t rMax = AliMUONConstants::Rmax(i/2);
    h2F = new TH2F(Form("hTrackerClusterHitMapForChamber%d",i+1), Form("cluster position distribution in chamber %d;X (cm);Y (cm)",i+1), 100, -rMax, rMax, 100, -rMax, rMax);
    Add2RecPointsList(h2F, AliMUONQAIndices::kTrackerClusterHitMapPerChamber+i, expert, !image);
  }
  
  // summary histograms per chamber
  h1I = new TH1I("hTrackerNumberOfClustersPerChamber", "Number of clusters per chamber;chamber ID;n_{clusters}", nCh,-0.5,nCh-0.5);
  Add2RecPointsList(h1I,AliMUONQAIndices::kTrackerNumberOfClustersPerChamber, !expert, image);
  
  h1F = new TH1F("hTrackerClusterMultiplicityPerChMean", "cluster mean size per chamber;chamber ID;<size> (n_{pads})", nCh,-0.5,nCh-0.5);
  h1F->SetOption("P");
  h1F->SetMarkerStyle(kFullDotMedium);
  h1F->SetMarkerColor(kRed);
  Add2RecPointsList(h1F, AliMUONQAIndices::kTrackerClusterMultiplicityPerChMean, !expert, image);
  
  h1F = new TH1F("hTrackerClusterMultiplicityPerChSigma", "cluster size dispersion per chamber;chamber ID;#sigma_{size} (n_{pads})", nCh,-0.5,nCh-0.5);
  h1F->SetOption("P");
  h1F->SetMarkerStyle(kFullDotMedium);
  h1F->SetMarkerColor(kRed);
  Add2RecPointsList(h1F, AliMUONQAIndices::kTrackerClusterMultiplicityPerChSigma, !expert, image);
  
  h1F = new TH1F("hTrackerClusterChargePerChMean", "cluster mean charge per chamber;chamber ID;<charge> (fC)", nCh,-0.5,nCh-0.5);
  h1F->SetOption("P");
  h1F->SetMarkerStyle(kFullDotMedium);
  h1F->SetMarkerColor(kRed);
  Add2RecPointsList(h1F, AliMUONQAIndices::kTrackerClusterChargePerChMean, !expert, image);
  
  h1F = new TH1F("hTrackerClusterChargePerChSigma", "cluster charge dispersion per chamber;chamber ID;#sigma_{charge} (fC)", nCh,-0.5,nCh-0.5);
  h1F->SetOption("P");
  h1F->SetMarkerStyle(kFullDotMedium);
  h1F->SetMarkerColor(kRed);
  Add2RecPointsList(h1F, AliMUONQAIndices::kTrackerClusterChargePerChSigma, !expert, image);
  
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
      Add2RecPointsList(h1I,AliMUONQAIndices::kTrackerClusterMultiplicityPerDE+detElemId, expert, !image);
      
      h1I = new TH1I(Form("hTrackerClusterChargeForDE%04d",detElemId), Form("cluster charge distribution in detection element %d;charge (fC)",detElemId), 100,0,1000);
      Add2RecPointsList(h1I,AliMUONQAIndices::kTrackerClusterChargePerDE+detElemId, expert, !image);
    }
    
    it.Next();
  }
  
  // summary histograms per DE
  h1I = new TH1I("hTrackerNumberOfClustersPerDE", "Number of clusters per detection element;DetElem ID;n_{clusters}", ndes+1,-0.5,ndes+0.5);
  Add2RecPointsList(h1I, AliMUONQAIndices::kTrackerNumberOfClustersPerDE, !expert, image);
  
  h1F = new TH1F("hTrackerClusterMultiplicityPerDEMean", "cluster mean size per DE;DetElem ID;<size> (n_{pads})", ndes+1,-0.5,ndes+0.5);
  h1F->SetOption("P");
  h1F->SetMarkerStyle(kFullDotMedium);
  h1F->SetMarkerColor(kRed);
  Add2RecPointsList(h1F, AliMUONQAIndices::kTrackerClusterMultiplicityPerDEMean, !expert, image);
  
  h1F = new TH1F("hTrackerClusterChargePerDEMean", "cluster mean charge per DE;DetElem ID;<charge> (fC)", ndes+1,-0.5,ndes+0.5);
  h1F->SetOption("P");
  h1F->SetMarkerStyle(kFullDotMedium);
  h1F->SetMarkerColor(kRed);
  Add2RecPointsList(h1F, AliMUONQAIndices::kTrackerClusterChargePerDEMean, !expert, image);
  
  if (!fMappingCheckRecPoints) fMappingCheckRecPoints = new AliMUONQAMappingCheck(RunNumber());  
}

//____________________________________________________________________________ 
void AliMUONTrackerQADataMakerRec::InitESDs()
{
  ///create ESDs histograms in ESDs subdir
  
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  Int_t nCh = AliMUONConstants::NTrackingCh();
  Int_t nDE = 1100;
  
  // track info
  TH1F* hESDnTracks = new TH1F("hESDnTracks", "number of tracks;n_{tracks}", 20, 0., 20.);
  Add2ESDsList(hESDnTracks, AliMUONQAIndices::kESDnTracks, !expert, image);

  TH1F* hESDMatchTrig = new TH1F("hESDMatchTrig", "number of tracks matched with trigger;n_{tracks}", 20, 0., 20.);
  Add2ESDsList(hESDMatchTrig, AliMUONQAIndices::kESDMatchTrig, !expert, image);
  
  TH1F* hESDMomentum = new TH1F("hESDMomentum", "momentum distribution;p (GeV/c)", 300, 0., 300);
  Add2ESDsList(hESDMomentum, AliMUONQAIndices::kESDMomentum, !expert, image);

  TH1F* hESDPt = new TH1F("hESDPt", "transverse momentum distribution;p_{t} (GeV/c)", 200, 0., 50);
  Add2ESDsList(hESDPt, AliMUONQAIndices::kESDPt, !expert, image);

  TH1F* hESDRapidity = new TH1F("hESDRapidity", "rapidity distribution;rapidity", 200, -4.5, -2.);
  Add2ESDsList(hESDRapidity, AliMUONQAIndices::kESDRapidity, !expert, image);

  TH1F* hESDChi2 = new TH1F("hESDChi2", "normalized #chi^{2} distribution;#chi^{2} / ndf", 500, 0., 50.);
  Add2ESDsList(hESDChi2, AliMUONQAIndices::kESDChi2, !expert, image);
  
  TH1F* hESDProbChi2 = new TH1F("hESDProbChi2", "distribution of probability of #chi^{2};prob(#chi^{2})", 100, 0., 1.);
  Add2ESDsList(hESDProbChi2, AliMUONQAIndices::kESDProbChi2, !expert, image);
  
  TH1F* hESDThetaX = new TH1F("hESDThetaX", "#theta_{X} distribution;#theta_{X} (degree)", 360, -180., 180);
  Add2ESDsList(hESDThetaX, AliMUONQAIndices::kESDThetaX, !expert, image);
  
  TH1F* hESDThetaY = new TH1F("hESDThetaY", "#theta_{Y} distribution;#theta_{Y} (degree)", 360, -180., 180);
  Add2ESDsList(hESDThetaY, AliMUONQAIndices::kESDThetaY, !expert, image);
  
  // cluster info
  for (Int_t i = 0; i < nCh; i++) {
    Float_t rMax = AliMUONConstants::Rmax(i/2);
    TH2F* hESDClusterHitMap = new TH2F(Form("hESDClusterHitMap%d",i+1), Form("cluster position distribution in chamber %d;X (cm);Y (cm)",i+1),
				       100, -rMax, rMax, 100, -rMax, rMax);
    Add2ESDsList(hESDClusterHitMap, AliMUONQAIndices::kESDClusterHitMap+i, expert, !image);
  }
  
  TH1F* hESDnClustersPerTrack = new TH1F("hESDnClustersPerTrack", "number of associated clusters per track;n_{clusters}", 20, 0., 20.);
  Add2ESDsList(hESDnClustersPerTrack, AliMUONQAIndices::kESDnClustersPerTrack, !expert, image);
  
  TH1F* hESDnClustersPerCh = new TH1F("hESDnClustersPerCh", "averaged number of clusters per chamber per track;chamber ID;<n_{clusters}>", nCh, -0.5, nCh-0.5);
  hESDnClustersPerCh->SetFillColor(kRed);
  Add2ESDsList(hESDnClustersPerCh, AliMUONQAIndices::kESDnClustersPerCh, !expert, image);
  
  TH1F* hESDnClustersPerDE = new TH1F("hESDnClustersPerDE", "averaged number of clusters per DE per track;DetElem ID;<n_{clusters}>", nDE+1, -0.5, nDE+0.5);
  hESDnClustersPerDE->SetFillColor(kRed);
  Add2ESDsList(hESDnClustersPerDE, AliMUONQAIndices::kESDnClustersPerDE, !expert, image);
  
  for (Int_t i = 0; i < nCh; i++) {
    TH1F* hESDClusterChargeInCh = new TH1F(Form("hESDClusterChargeInCh%d",i+1), Form("cluster charge distribution in chamber %d;charge (fC)",i+1), 100, 0., 1000.);
    Add2ESDsList(hESDClusterChargeInCh, AliMUONQAIndices::kESDClusterChargeInCh+i, expert, !image);
  }
  
  TH1F* hESDClusterChargePerChMean = new TH1F("hESDClusterChargePerChMean", "cluster mean charge per chamber;chamber ID;<charge> (fC)", nCh, -0.5, nCh-0.5);
  hESDClusterChargePerChMean->SetOption("P");
  hESDClusterChargePerChMean->SetMarkerStyle(kFullDotMedium);
  hESDClusterChargePerChMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDClusterChargePerChMean, AliMUONQAIndices::kESDClusterChargePerChMean, !expert, image);
  
  TH1F* hESDClusterChargePerChSigma = new TH1F("hESDClusterChargePerChSigma", "cluster charge dispersion per chamber;chamber ID;#sigma_{charge} (fC)", nCh, -0.5, nCh-0.5);
  hESDClusterChargePerChSigma->SetOption("P");
  hESDClusterChargePerChSigma->SetMarkerStyle(kFullDotMedium);
  hESDClusterChargePerChSigma->SetMarkerColor(kRed);
  Add2ESDsList(hESDClusterChargePerChSigma, AliMUONQAIndices::kESDClusterChargePerChSigma, !expert, image);
  
  TH1F* hESDClusterChargePerDE = new TH1F("hESDClusterChargePerDE", "cluster mean charge per DE;DetElem ID;<charge> (fC)", nDE+1, -0.5, nDE+0.5);
  hESDClusterChargePerDE->SetOption("P");
  hESDClusterChargePerDE->SetMarkerStyle(kFullDotMedium);
  hESDClusterChargePerDE->SetMarkerColor(kRed);
  Add2ESDsList(hESDClusterChargePerDE, AliMUONQAIndices::kESDClusterChargePerDE, !expert, image);
  
  for (Int_t i = 0; i < nCh; i++) {
    TH1F* hESDClusterSizeInCh = new TH1F(Form("hESDClusterSizeInCh%d",i+1), Form("cluster size distribution in chamber %d;size (n_{pads})",i+1), 200, 0., 200.);
    Add2ESDsList(hESDClusterSizeInCh, AliMUONQAIndices::kESDClusterSizeInCh+i, expert, !image);
  }
  
  TH1F* hESDClusterSizePerChMean = new TH1F("hESDClusterSizePerChMean", "cluster mean size per chamber;chamber ID;<size> (n_{pads})", nCh, -0.5, nCh-0.5);
  hESDClusterSizePerChMean->SetOption("P");
  hESDClusterSizePerChMean->SetMarkerStyle(kFullDotMedium);
  hESDClusterSizePerChMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDClusterSizePerChMean, AliMUONQAIndices::kESDClusterSizePerChMean, !expert, image);
  
  TH1F* hESDClusterSizePerChSigma = new TH1F("hESDClusterSizePerChSigma", "cluster size dispersion per chamber;chamber ID;#sigma_{size} (n_{pads})", nCh, -0.5, nCh-0.5);
  hESDClusterSizePerChSigma->SetOption("P");
  hESDClusterSizePerChSigma->SetMarkerStyle(kFullDotMedium);
  hESDClusterSizePerChSigma->SetMarkerColor(kRed);
  Add2ESDsList(hESDClusterSizePerChSigma, AliMUONQAIndices::kESDClusterSizePerChSigma, !expert, image);
  
  TH1F* hESDClusterSizePerDE = new TH1F("hESDClusterSizePerDE", "cluster mean size per DE;DetElem ID;<size> (n_{pads})", nDE+1, -0.5, nDE+0.5);
  hESDClusterSizePerDE->SetOption("P");
  hESDClusterSizePerDE->SetMarkerStyle(kFullDotMedium);
  hESDClusterSizePerDE->SetMarkerColor(kRed);
  Add2ESDsList(hESDClusterSizePerDE, AliMUONQAIndices::kESDClusterSizePerDE, !expert, image);
  
  // cluster - track info
  for (Int_t i = 0; i < nCh; i++) {
    TH1F* hESDResidualXInCh = new TH1F(Form("hESDResidualXInCh%d",i+1), Form("cluster-track residual-X distribution in chamber %d;#Delta_{X} (cm)",i+1), 1000, -5., 5.);
    Add2ESDsList(hESDResidualXInCh, AliMUONQAIndices::kESDResidualXInCh+i, expert, !image);
    
    TH1F* hESDResidualYInCh = new TH1F(Form("hESDResidualYInCh%d",i+1), Form("cluster-track residual-Y distribution in chamber %d;#Delta_{Y} (cm)",i+1), 1000, -1., 1.);
    Add2ESDsList(hESDResidualYInCh, AliMUONQAIndices::kESDResidualYInCh+i, expert, !image);
    
    TH1F* hESDLocalChi2XInCh = new TH1F(Form("hESDLocalChi2XInCh%d",i+1), Form("local chi2-X distribution in chamber %d;local #chi^{2}_{X}",i+1), 1000, 0., 25);
    Add2ESDsList(hESDLocalChi2XInCh, AliMUONQAIndices::kESDLocalChi2XInCh+i, expert, !image);
    
    TH1F* hESDLocalChi2YInCh = new TH1F(Form("hESDLocalChi2YInCh%d",i+1), Form("local chi2-Y distribution in chamber %d;local #chi^{2}_{Y}",i+1), 1000, 0., 25);
    Add2ESDsList(hESDLocalChi2YInCh, AliMUONQAIndices::kESDLocalChi2YInCh+i, expert, !image);
    
    TH1F* hESDLocalChi2InCh = new TH1F(Form("hESDLocalChi2InCh%d",i+1), Form("local chi2 (~0.5*(#chi^{2}_{X}+#chi^{2}_{Y})) distribution in chamber %d;local #chi^{2}",i+1), 1000, 0., 25);
    Add2ESDsList(hESDLocalChi2InCh, AliMUONQAIndices::kESDLocalChi2InCh+i, expert, !image);
  }
  
  TH1F* hESDResidualXPerChMean = new TH1F("hESDResidualXPerChMean", "cluster-track residual-X per Ch: mean;chamber ID;<#Delta_{X}> (cm)", nCh, -0.5, nCh-0.5);
  hESDResidualXPerChMean->SetOption("P");
  hESDResidualXPerChMean->SetMarkerStyle(kFullDotMedium);
  hESDResidualXPerChMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDResidualXPerChMean, AliMUONQAIndices::kESDResidualXPerChMean, !expert, image);
  
  TH1F* hESDResidualYPerChMean = new TH1F("hESDResidualYPerChMean", "cluster-track residual-Y per Ch: mean;chamber ID;<#Delta_{Y}> (cm)", nCh, -0.5, nCh-0.5);
  hESDResidualYPerChMean->SetOption("P");
  hESDResidualYPerChMean->SetMarkerStyle(kFullDotMedium);
  hESDResidualYPerChMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDResidualYPerChMean, AliMUONQAIndices::kESDResidualYPerChMean, !expert, image);
  
  TH1F* hESDResidualXPerChSigma = new TH1F("hESDResidualXPerChSigma", "cluster-track residual-X per Ch: sigma;chamber ID;#sigma_{X} (cm)", nCh, -0.5, nCh-0.5);
  hESDResidualXPerChSigma->SetOption("P");
  hESDResidualXPerChSigma->SetMarkerStyle(kFullDotMedium);
  hESDResidualXPerChSigma->SetMarkerColor(kRed);
  Add2ESDsList(hESDResidualXPerChSigma, AliMUONQAIndices::kESDResidualXPerChSigma, !expert, image);
  
  TH1F* hESDResidualYPerChSigma = new TH1F("hESDResidualYPerChSigma", "cluster-track residual-Y per Ch: sigma;chamber ID;#sigma_{Y} (cm)", nCh, -0.5, nCh-0.5);
  hESDResidualYPerChSigma->SetOption("P");
  hESDResidualYPerChSigma->SetMarkerStyle(kFullDotMedium);
  hESDResidualYPerChSigma->SetMarkerColor(kRed);
  Add2ESDsList(hESDResidualYPerChSigma, AliMUONQAIndices::kESDResidualYPerChSigma, !expert, image);
  
  TH1F* hESDLocalChi2XPerChMean = new TH1F("hESDLocalChi2XPerCh", "local chi2-X per Ch: mean;chamber ID;<local #chi^{2}_{X}>", nCh, -0.5, nCh-0.5);
  hESDLocalChi2XPerChMean->SetOption("P");
  hESDLocalChi2XPerChMean->SetMarkerStyle(kFullDotMedium);
  hESDLocalChi2XPerChMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDLocalChi2XPerChMean, AliMUONQAIndices::kESDLocalChi2XPerChMean, !expert, image);
  
  TH1F* hESDLocalChi2YPerChMean = new TH1F("hESDLocalChi2YPerCh", "local chi2-Y per Ch: mean;chamber ID;<local #chi^{2}_{Y}>", nCh, -0.5, nCh-0.5);
  hESDLocalChi2YPerChMean->SetOption("P");
  hESDLocalChi2YPerChMean->SetMarkerStyle(kFullDotMedium);
  hESDLocalChi2YPerChMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDLocalChi2YPerChMean, AliMUONQAIndices::kESDLocalChi2YPerChMean, !expert, image);
  
  TH1F* hESDLocalChi2PerChMean = new TH1F("hESDLocalChi2PerCh", "local chi2 (~0.5*(#chi^{2}_{X}+#chi^{2}_{Y})) per Ch: mean;chamber ID;<local #chi^{2}>", nCh, -0.5, nCh-0.5);
  hESDLocalChi2PerChMean->SetOption("P");
  hESDLocalChi2PerChMean->SetMarkerStyle(kFullDotMedium);
  hESDLocalChi2PerChMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDLocalChi2PerChMean, AliMUONQAIndices::kESDLocalChi2PerChMean, !expert, image);
  
  TH1F* hESDResidualXPerDEMean = new TH1F("hESDResidualXPerDEMean", "cluster-track residual-X per DE: mean;DetElem ID;<#Delta_{X}> (cm)", nDE+1, -0.5, nDE+0.5);
  hESDResidualXPerDEMean->SetOption("P");
  hESDResidualXPerDEMean->SetMarkerStyle(kFullDotMedium);
  hESDResidualXPerDEMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDResidualXPerDEMean, AliMUONQAIndices::kESDResidualXPerDEMean, !expert, image);
  
  TH1F* hESDResidualYPerDEMean = new TH1F("hESDResidualYPerDEMean", "cluster-track residual-Y per DE: mean;DetElem ID;<#Delta_{Y}> (cm)", nDE+1, -0.5, nDE+0.5);
  hESDResidualYPerDEMean->SetOption("P");
  hESDResidualYPerDEMean->SetMarkerStyle(kFullDotMedium);
  hESDResidualYPerDEMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDResidualYPerDEMean, AliMUONQAIndices::kESDResidualYPerDEMean, !expert, image);
  
  TH1F* hESDResidualXPerDESigma = new TH1F("hESDResidualXPerDESigma", "cluster-track residual-X per DE: sigma;DetElem ID;#sigma_{X} (cm)", nDE+1, -0.5, nDE+0.5);
  hESDResidualXPerDESigma->SetOption("P");
  hESDResidualXPerDESigma->SetMarkerStyle(kFullDotMedium);
  hESDResidualXPerDESigma->SetMarkerColor(kRed);
  Add2ESDsList(hESDResidualXPerDESigma, AliMUONQAIndices::kESDResidualXPerDESigma, !expert, image);
  
  TH1F* hESDResidualYPerDESigma = new TH1F("hESDResidualYPerDESigma", "cluster-track residual-Y per DE: sigma;DetElem ID;#sigma_{Y} (cm)", nDE+1, -0.5, nDE+0.5);
  hESDResidualYPerDESigma->SetOption("P");
  hESDResidualYPerDESigma->SetMarkerStyle(kFullDotMedium);
  hESDResidualYPerDESigma->SetMarkerColor(kRed);
  Add2ESDsList(hESDResidualYPerDESigma, AliMUONQAIndices::kESDResidualYPerDESigma, !expert, image);
  
  TH1F* hESDLocalChi2XPerDEMean = new TH1F("hESDLocalChi2XPerDE", "local chi2-X per DE: mean;DetElem ID;<local #chi^{2}_{X}>", nDE+1, -0.5, nDE+0.5);
  hESDLocalChi2XPerDEMean->SetOption("P");
  hESDLocalChi2XPerDEMean->SetMarkerStyle(kFullDotMedium);
  hESDLocalChi2XPerDEMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDLocalChi2XPerDEMean, AliMUONQAIndices::kESDLocalChi2XPerDEMean, !expert, image);
  
  TH1F* hESDLocalChi2YPerDEMean = new TH1F("hESDLocalChi2YPerDE", "local chi2-Y per DE: mean;DetElem ID;<local #chi^{2}_{Y}>", nDE+1, -0.5, nDE+0.5);
  hESDLocalChi2YPerDEMean->SetOption("P");
  hESDLocalChi2YPerDEMean->SetMarkerStyle(kFullDotMedium);
  hESDLocalChi2YPerDEMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDLocalChi2YPerDEMean, AliMUONQAIndices::kESDLocalChi2YPerDEMean, !expert, image);
  
  TH1F* hESDLocalChi2PerDEMean = new TH1F("hESDLocalChi2PerDE", "local chi2 (~0.5*(#chi^{2}_{X}+#chi^{2}_{Y})) per DE: mean;DetElem ID;<local #chi^{2}>", nDE+1, -0.5, nDE+0.5);
  hESDLocalChi2PerDEMean->SetOption("P");
  hESDLocalChi2PerDEMean->SetMarkerStyle(kFullDotMedium);
  hESDLocalChi2PerDEMean->SetMarkerColor(kRed);
  Add2ESDsList(hESDLocalChi2PerDEMean, AliMUONQAIndices::kESDLocalChi2PerDEMean, !expert, image);
  
  // intermediate histograms
  TH1F* hESDnTotClustersPerCh = new TH1F("hESDnTotClustersPerCh", "total number of associated clusters per chamber;chamber ID;#Sigma(n_{clusters})", nCh, -0.5, nCh-0.5);
  Add2ESDsList(hESDnTotClustersPerCh, AliMUONQAIndices::kESDnTotClustersPerCh, expert, !image);
  TH1F* hESDnTotClustersPerDE = new TH1F("hESDnTotClustersPerDE", "total number of associated clusters per DE;DetElem ID;#Sigma(n_{clusters})", nDE+1, -0.5, nDE+0.5);
  Add2ESDsList(hESDnTotClustersPerDE, AliMUONQAIndices::kESDnTotClustersPerDE, expert, !image);
  TH1F* hESDnTotFullClustersPerDE = new TH1F("hESDnTotFullClustersPerDE", "total number of associated clusters containing pad info per DE;DetElem ID;#Sigma(n_{full clusters})", nDE+1, -0.5, nDE+0.5);
  Add2ESDsList(hESDnTotFullClustersPerDE, AliMUONQAIndices::kESDnTotFullClustersPerDE, expert, !image);
  TH1F* hESDSumClusterChargePerDE = new TH1F("hESDSumClusterChargePerDE", "sum of cluster charge per DE;DetElem ID;#Sigma(charge) (fC)", nDE+1, -0.5, nDE+0.5);
  Add2ESDsList(hESDSumClusterChargePerDE, AliMUONQAIndices::kESDSumClusterChargePerDE, expert, !image);
  TH1F* hESDSumClusterSizePerDE = new TH1F("hESDSumClusterSizePerDE", "sum of cluster size per DE;DetElem ID;#Sigma(size) (n_{pads})", nDE+1, -0.5, nDE+0.5);
  Add2ESDsList(hESDSumClusterSizePerDE, AliMUONQAIndices::kESDSumClusterSizePerDE, expert, !image);
  TH1F* hESDSumResidualXPerDE = new TH1F("hESDSumResidualXPerDE", "sum of cluster-track residual-X per DE;DetElem ID;#Sigma(#Delta_{X}) (cm)", nDE+1, -0.5, nDE+0.5);
  Add2ESDsList(hESDSumResidualXPerDE, AliMUONQAIndices::kESDSumResidualXPerDE, expert, !image);
  TH1F* hESDSumResidualYPerDE = new TH1F("hESDSumResidualYPerDE", "sum of cluster-track residual-Y per DE;DetElem ID;#Sigma(#Delta_{Y}) (cm)", nDE+1, -0.5, nDE+0.5);
  Add2ESDsList(hESDSumResidualYPerDE, AliMUONQAIndices::kESDSumResidualYPerDE, expert, !image);
  TH1F* hESDSumResidualX2PerDE = new TH1F("hESDSumResidualX2PerDE", "sum of cluster-track residual-X**2 per DE;DetElem ID;#Sigma(#Delta_{X}^{2}) (cm^{2})", nDE+1, -0.5, nDE+0.5);
  Add2ESDsList(hESDSumResidualX2PerDE, AliMUONQAIndices::kESDSumResidualX2PerDE, expert, !image);
  TH1F* hESDSumResidualY2PerDE = new TH1F("hESDSumResidualY2PerDE", "sum of cluster-track residual-Y**2 per DE;DetElem ID;#Sigma(#Delta_{Y}^{2}) (cm^{2})", nDE+1, -0.5, nDE+0.5);
  Add2ESDsList(hESDSumResidualY2PerDE, AliMUONQAIndices::kESDSumResidualY2PerDE, expert, !image);
  TH1F* hESDSumLocalChi2XPerDE = new TH1F("hESDSumLocalChi2XPerDE", "sum of local chi2-X per DE;DetElem ID;#Sigma(local #chi^{2}_{X})", nDE+1, -0.5, nDE+0.5);
  Add2ESDsList(hESDSumLocalChi2XPerDE, AliMUONQAIndices::kESDSumLocalChi2XPerDE, expert, !image);
  TH1F* hESDSumLocalChi2YPerDE = new TH1F("hESDSumLocalChi2YPerDE", "sum of local chi2-Y per DE;DetElem ID;#Sigma(local #chi^{2}_{Y})", nDE+1, -0.5, nDE+0.5);
  Add2ESDsList(hESDSumLocalChi2YPerDE, AliMUONQAIndices::kESDSumLocalChi2YPerDE, expert, !image);
  TH1F* hESDSumLocalChi2PerDE = new TH1F("hESDSumLocalChi2PerDE", "sum of local chi2 (~0.5*(#chi^{2}_{X}+#chi^{2}_{Y})) per DE;DetElem ID;#Sigma(local #chi^{2})", nDE+1, -0.5, nDE+0.5);
  Add2ESDsList(hESDSumLocalChi2PerDE, AliMUONQAIndices::kESDSumLocalChi2PerDE, expert, !image);
}

//____________________________________________________________________________
void AliMUONTrackerQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
	/// make QA for rawdata tracker
  	
  AliCodeTimerAuto("",0);

  /// forces init
  GetRawsData(AliMUONQAIndices::kTrackerBusPatchOccupancy);
  
	((AliMUONTrackerDataMaker*)fTrackerDataMaker)->SetRawReader(rawReader);
	
	fTrackerDataMaker->ProcessEvent();
  
  
}

//__________________________________________________________________
void AliMUONTrackerQADataMakerRec::MakeDigits(TTree* digitsTree)         
{
  /// makes data from Digits

  AliCodeTimerAuto("",0);
  
  // Do nothing in case of calibration event
  if ( GetRecoParam()->GetEventSpecie() == AliRecoParam::kCalib ) return;

  if (!fDigitStore)
    fDigitStore = AliMUONVDigitStore::Create(*digitsTree);
  fDigitStore->Connect(*digitsTree, false);
  digitsTree->GetEvent(0);
  
  TIter next(fDigitStore->CreateIterator());
  
  AliMUONVDigit* dig = 0x0;
  
  while ( ( dig = static_cast<AliMUONVDigit*>(next()) ) )
    {
    GetDigitsData(0)->Fill(dig->DetElemId());
    GetDigitsData(1)->Fill(dig->ADC());
    }
}


//____________________________________________________________________________
void AliMUONTrackerQADataMakerRec::MakeRecPoints(TTree* clustersTree)
{
	/// Fill histograms related to tracker clusters 
	
	// First things first : do we have clusters in the TreeR ?
	// In "normal" production mode, it should be perfectly normal
	// *not* to have them.
	// But if for some reason we de-activated the combined tracking,
	// then we have clusters in TreeR, so let's take that opportunity
	// to QA them...
	
  AliCodeTimerAuto("",0);

  // Do nothing in case of calibration event
  if ( GetRecoParam()->GetEventSpecie() == AliRecoParam::kCalib ) return;

	if (!fClusterStore)
	{
	  AliCodeTimerAuto("ClusterStore creation",1);
		fClusterStore = AliMUONVClusterStore::Create(*clustersTree);
		if (!fClusterStore) 
		{
			return;
		}
	}
	
	fClusterStore->Connect(*clustersTree,kFALSE);
	clustersTree->GetEvent(0);

	TIter next(fClusterStore->CreateIterator());
	AliMUONVCluster* cluster;
	
  if ( fMappingCheckRecPoints ) fMappingCheckRecPoints->NewEvent();
  
	while ( ( cluster = static_cast<AliMUONVCluster*>(next()) ) )
	{
		Int_t detElemId = cluster->GetDetElemId();
		Int_t chamberId = AliMpDEManager::GetChamberId(detElemId);
		
		GetRecPointsData(AliMUONQAIndices::kTrackerNumberOfClustersPerDE)->Fill(detElemId);
		GetRecPointsData(AliMUONQAIndices::kTrackerClusterChargePerDE+detElemId)->Fill(cluster->GetCharge());
		GetRecPointsData(AliMUONQAIndices::kTrackerClusterMultiplicityPerDE+detElemId)->Fill(cluster->GetNDigits());

		GetRecPointsData(AliMUONQAIndices::kTrackerNumberOfClustersPerChamber)->Fill(chamberId);
		GetRecPointsData(AliMUONQAIndices::kTrackerClusterChargePerChamber+chamberId)->Fill(cluster->GetCharge());
		GetRecPointsData(AliMUONQAIndices::kTrackerClusterMultiplicityPerChamber+chamberId)->Fill(cluster->GetNDigits());
	        GetRecPointsData(AliMUONQAIndices::kTrackerClusterHitMapPerChamber+chamberId)->Fill(cluster->GetX(),cluster->GetY());
		
    if ( fMappingCheckRecPoints ) fMappingCheckRecPoints->Store(*cluster);
    
	}
	
	fClusterStore->Clear();
}

//____________________________________________________________________________
void AliMUONTrackerQADataMakerRec::MakeESDs(AliESDEvent* esd)
{
  /// make QA data from ESDs

  AliCodeTimerAuto("",0);
  
  // Do nothing in case of calibration event
  if ( GetRecoParam()->GetEventSpecie() == AliRecoParam::kCalib ) return;
  
  // load ESD event in the interface
  AliMUONESDInterface esdInterface;
  if (GetRecoParam()) AliMUONESDInterface::ResetTracker(GetRecoParam(), kFALSE);
  else AliError("Unable to get recoParam: use default ones for residual calculation");
  esdInterface.LoadEvent(*esd);
  
  GetESDsData(AliMUONQAIndices::kESDnTracks)->Fill(esdInterface.GetNTracks());
  
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
    
    GetESDsData(AliMUONQAIndices::kESDMomentum)->Fill(esdTrack->P());
    GetESDsData(AliMUONQAIndices::kESDPt)->Fill(esdTrack->Pt());
    GetESDsData(AliMUONQAIndices::kESDRapidity)->Fill(esdTrack->Y());
    GetESDsData(AliMUONQAIndices::kESDChi2)->Fill(track->GetNormalizedChi2());
    GetESDsData(AliMUONQAIndices::kESDProbChi2)->Fill(TMath::Prob(track->GetGlobalChi2(),track->GetNDF()));
    GetESDsData(AliMUONQAIndices::kESDThetaX)->Fill(esdTrack->GetThetaXUncorrected() / TMath::Pi() * 180.);
    GetESDsData(AliMUONQAIndices::kESDThetaY)->Fill(esdTrack->GetThetaYUncorrected() / TMath::Pi() * 180.);
    GetESDsData(AliMUONQAIndices::kESDnClustersPerTrack)->Fill(track->GetNClusters());
    
    // loop over clusters
    AliMUONTrackParam* trackParam = static_cast<AliMUONTrackParam*>(track->GetTrackParamAtCluster()->First());
    while (trackParam) {
      
      AliMUONVCluster* cluster = trackParam->GetClusterPtr();
      Int_t chId = cluster->GetChamberId();
      Int_t deID = cluster->GetDetElemId();
      Double_t residualX = cluster->GetX() - trackParam->GetNonBendingCoor();
      Double_t residualY = cluster->GetY() - trackParam->GetBendingCoor();
      Double_t sigmaResidualX2 = cluster->GetErrX2() - trackParam->GetCovariances()(0,0);
      Double_t sigmaResidualY2 = cluster->GetErrY2() - trackParam->GetCovariances()(2,2);
      Double_t localChi2X = (sigmaResidualX2 > 0.) ? residualX*residualX/sigmaResidualX2 : 0.;
      Double_t localChi2Y = (sigmaResidualY2 > 0.) ? residualY*residualY/sigmaResidualY2 : 0.;
      Double_t localChi2 = 0.5 * trackParam->GetLocalChi2();
      
      GetESDsData(AliMUONQAIndices::kESDClusterHitMap+chId)->Fill(cluster->GetX(), cluster->GetY());
      
      GetESDsData(AliMUONQAIndices::kESDnTotClustersPerCh)->Fill(chId);
      GetESDsData(AliMUONQAIndices::kESDnTotClustersPerDE)->Fill(deID);
      
      GetESDsData(AliMUONQAIndices::kESDClusterChargeInCh+chId)->Fill(cluster->GetCharge());
      GetESDsData(AliMUONQAIndices::kESDSumClusterChargePerDE)->Fill(deID, cluster->GetCharge());
      
      if (cluster->GetNDigits() > 0) { // discard clusters with pad not stored in ESD
	GetESDsData(AliMUONQAIndices::kESDnTotFullClustersPerDE)->Fill(deID);
        GetESDsData(AliMUONQAIndices::kESDClusterSizeInCh+chId)->Fill(cluster->GetNDigits());
	GetESDsData(AliMUONQAIndices::kESDSumClusterSizePerDE)->Fill(deID, cluster->GetNDigits());
      }
      
      GetESDsData(AliMUONQAIndices::kESDResidualXInCh+chId)->Fill(residualX);
      GetESDsData(AliMUONQAIndices::kESDResidualYInCh+chId)->Fill(residualY);
      GetESDsData(AliMUONQAIndices::kESDSumResidualXPerDE)->Fill(deID, residualX);
      GetESDsData(AliMUONQAIndices::kESDSumResidualYPerDE)->Fill(deID, residualY);
      GetESDsData(AliMUONQAIndices::kESDSumResidualX2PerDE)->Fill(deID, residualX*residualX);
      GetESDsData(AliMUONQAIndices::kESDSumResidualY2PerDE)->Fill(deID, residualY*residualY);
      
      GetESDsData(AliMUONQAIndices::kESDLocalChi2XInCh+chId)->Fill(localChi2X);
      GetESDsData(AliMUONQAIndices::kESDLocalChi2YInCh+chId)->Fill(localChi2Y);
      GetESDsData(AliMUONQAIndices::kESDLocalChi2InCh+chId)->Fill(localChi2);
      GetESDsData(AliMUONQAIndices::kESDSumLocalChi2XPerDE)->Fill(deID, localChi2X);
      GetESDsData(AliMUONQAIndices::kESDSumLocalChi2YPerDE)->Fill(deID, localChi2Y);
      GetESDsData(AliMUONQAIndices::kESDSumLocalChi2PerDE)->Fill(deID, localChi2);
      
      trackParam = static_cast<AliMUONTrackParam*>(track->GetTrackParamAtCluster()->After(trackParam));
      
    }
    
  }

  GetESDsData(AliMUONQAIndices::kESDMatchTrig)->Fill(nTrackMatchTrig);
  
}

//____________________________________________________________________________ 
AliMUONVTrackerData* AliMUONTrackerQADataMakerRec::GetTrackerData() const
{ 
/// Return tracker data
  
  return fTrackerDataMaker->Data(); 
  
}

//____________________________________________________________________________ 
void
AliMUONTrackerQADataMakerRec::ResetDetectorRaws(TObjArray* list)
{
  /// Reset those histograms that must be reset (and only those), plus
  /// the trackerdata itself.
  
  TIter next(list);
  TObject* o;
  while ( ( o = next() ) )
  {
    TH1* h = dynamic_cast<TH1*>(o);
    if (h)
    {
      TString hn(h->GetName());
      
      if ( hn.Contains("Tracker") )
      {
        if ( !hn.Contains("hTrackerBusPatchNofPads") && 
            !hn.Contains("hTrackerBusPatchNofManus") &&
            !hn.Contains("hTrackerBusPatchConfig" ) )
        {
          AliDebug(1,Form("Resetting %s",hn.Data()));
          h->Reset();                  
        }
      }
      else
      {
        AliDebug(1,Form("Will not reset histogram %s",hn.Data()));
      }
    }
    else
    {
      AliMUONVTrackerData* d = dynamic_cast<AliMUONVTrackerData*>(o);
      if (d)
      {
        AliDebug(1,Form("Resetting %s",d->GetName()));
        d->Clear();
      }
      else
      {
        AliError(Form("Found an object of class %s. Do not know how to reset.",o->ClassName()));
      }
    }
  }
}
