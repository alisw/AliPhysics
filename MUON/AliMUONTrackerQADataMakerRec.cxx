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

#include "AliCDBManager.h"
#include "AliCodeTimer.h"
#include "AliDAQ.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliLog.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONConstants.h"  
#include "AliMUONDigitMaker.h"
#include "AliMUONESDInterface.h"
#include "AliMUONLogger.h"
#include "AliMUONQADataMakerRec.h"
#include "AliMUONQAIndices.h"
#include "AliMUONQAMappingCheck.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackerData.h"
#include "AliMUONTrackerDataMaker.h"
#include "AliMUONVCluster.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"
#include "AliMpBusPatch.h"
#include "AliMpConstants.h"
#include "AliMpDDL.h"
#include "AliMpDDLStore.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "AliMpManuIterator.h"
#include "AliQAv1.h"
#include "AliRawReader.h"
#include "AliRawEventHeaderBase.h"
#include <Riostream.h>
#include <TH1F.h> 
#include <TH1I.h> 
#include <TH2F.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPaveText.h>

//-----------------------------------------------------------------------------
/// \class AliMUONTrackerQADataMakerRec
///
/// Quality assurance data (histo) maker for MUON tracker
///
///
/// Note that all the methods of this class shoud not be called when eventSpecie is AliRecoParam::kCalib !
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
fDigitStore(0x0),
fDigitMaker(new AliMUONDigitMaker(kTRUE)),
fClusterStore(0x0),
fCalibrationData(new AliMUONCalibrationData(AliCDBManager::Instance()->GetRun())),
fLogger(0x0),
fBusPatchConfig(0x0),
fBPxmin(0),
fBPxmax(0),
fBPnbins(0),
fTrackerDataMakerArray(0x0),
fTrackerCalDataArray(0x0),
fTrackerRecDataArray(0x0),
fMappingCheckRecPointsArray(0x0)
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
  delete fMappingCheckRecPointsArray;
  if (fLogger)
  {
    if ( fLogger->NumberOfEntries() != 0 ) 
    {
      AliError("Some unprocessed logged errors are still there... Please check");
    }
    delete fLogger;    
  }
  delete fTrackerDataMakerArray;
  delete fTrackerCalDataArray;
  delete fTrackerRecDataArray;
  delete fBusPatchConfig;  
}

//____________________________________________________________________________ 
void AliMUONTrackerQADataMakerRec::InsertTrackerData(Int_t specie, 
                                                     TObjArray** list,
                                                     TObject* object, 
                                                     Int_t indexNumber,
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
      if ( !strcmp(object->GetName(),o->GetName()) ) 
      {
        alreadyThere = kTRUE;
        old = o;
      }
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
  //
  Bool_t firstFill = kTRUE;
  //
  for (int itc=-1;itc<Master()->GetNTrigClasses();itc++) { // RS: loop over eventual clones per trigger class
    //
    TH1* hESDnClustersPerTr = Master()->GetESDsData(AliMUONQAIndices::kESDnClustersPerTrack, itc);
    if (! hESDnClustersPerTr) continue;
    Double_t nTracks = hESDnClustersPerTr->GetEntries();
    if (nTracks <= 0) continue;
    //
    if (firstFill) { AliCodeTimerAuto("",0); firstFill = kFALSE;}
  
    TH1* hESDnClustersPerCh = Master()->GetESDsData(AliMUONQAIndices::kESDnClustersPerCh,itc);
    TH1* hESDnClustersPerDE = Master()->GetESDsData(AliMUONQAIndices::kESDnClustersPerDE,itc);
    TH1* hESDClusterChargePerChMean = Master()->GetESDsData(AliMUONQAIndices::kESDClusterChargePerChMean,itc);
    TH1* hESDClusterChargePerChSigma = Master()->GetESDsData(AliMUONQAIndices::kESDClusterChargePerChSigma,itc);
    TH1* hESDClusterSizePerChMean = Master()->GetESDsData(AliMUONQAIndices::kESDClusterSizePerChMean,itc);
    TH1* hESDClusterSizePerChSigma = Master()->GetESDsData(AliMUONQAIndices::kESDClusterSizePerChSigma,itc);
    TH1* hESDResidualXPerChMean = Master()->GetESDsData(AliMUONQAIndices::kESDResidualXPerChMean,itc);
    TH1* hESDResidualXPerChSigma = Master()->GetESDsData(AliMUONQAIndices::kESDResidualXPerChSigma,itc);
    TH1* hESDResidualYPerChMean = Master()->GetESDsData(AliMUONQAIndices::kESDResidualYPerChMean,itc);
    TH1* hESDResidualYPerChSigma = Master()->GetESDsData(AliMUONQAIndices::kESDResidualYPerChSigma,itc);
    TH1* hESDLocalChi2XPerChMean = Master()->GetESDsData(AliMUONQAIndices::kESDLocalChi2XPerChMean,itc);
    TH1* hESDLocalChi2YPerChMean = Master()->GetESDsData(AliMUONQAIndices::kESDLocalChi2YPerChMean,itc);
    TH1* hESDLocalChi2PerChMean = Master()->GetESDsData(AliMUONQAIndices::kESDLocalChi2PerChMean,itc);
    TH1* hESDClusterChargePerDE = Master()->GetESDsData(AliMUONQAIndices::kESDClusterChargePerDE,itc);
    TH1* hESDClusterSizePerDE = Master()->GetESDsData(AliMUONQAIndices::kESDClusterSizePerDE,itc);
    TH1* hESDResidualXPerDEMean = Master()->GetESDsData(AliMUONQAIndices::kESDResidualXPerDEMean,itc);
    TH1* hESDResidualXPerDESigma = Master()->GetESDsData(AliMUONQAIndices::kESDResidualXPerDESigma,itc);
    TH1* hESDResidualYPerDEMean = Master()->GetESDsData(AliMUONQAIndices::kESDResidualYPerDEMean,itc);
    TH1* hESDResidualYPerDESigma = Master()->GetESDsData(AliMUONQAIndices::kESDResidualYPerDESigma,itc);
    TH1* hESDLocalChi2XPerDEMean = Master()->GetESDsData(AliMUONQAIndices::kESDLocalChi2XPerDEMean,itc);
    TH1* hESDLocalChi2YPerDEMean = Master()->GetESDsData(AliMUONQAIndices::kESDLocalChi2YPerDEMean,itc);
    TH1* hESDLocalChi2PerDEMean = Master()->GetESDsData(AliMUONQAIndices::kESDLocalChi2PerDEMean,itc);
    TH1* hESDnTotClustersPerCh = Master()->GetESDsData(AliMUONQAIndices::kESDnTotClustersPerCh,itc);
    TH1* hESDnTotClustersPerDE = Master()->GetESDsData(AliMUONQAIndices::kESDnTotClustersPerDE,itc);
    TH1* hESDnTotFullClustersPerDE = Master()->GetESDsData(AliMUONQAIndices::kESDnTotFullClustersPerDE,itc);
    TH1* hESDSumClusterChargePerDE = Master()->GetESDsData(AliMUONQAIndices::kESDSumClusterChargePerDE,itc);
    TH1* hESDSumClusterSizePerDE = Master()->GetESDsData(AliMUONQAIndices::kESDSumClusterSizePerDE,itc);
    TH1* hESDSumResidualXPerDE = Master()->GetESDsData(AliMUONQAIndices::kESDSumResidualXPerDE,itc);
    TH1* hESDSumResidualYPerDE = Master()->GetESDsData(AliMUONQAIndices::kESDSumResidualYPerDE,itc);
    TH1* hESDSumResidualX2PerDE = Master()->GetESDsData(AliMUONQAIndices::kESDSumResidualX2PerDE,itc);
    TH1* hESDSumResidualY2PerDE = Master()->GetESDsData(AliMUONQAIndices::kESDSumResidualY2PerDE,itc);
    TH1* hESDSumLocalChi2XPerDE = Master()->GetESDsData(AliMUONQAIndices::kESDSumLocalChi2XPerDE,itc);
    TH1* hESDSumLocalChi2YPerDE = Master()->GetESDsData(AliMUONQAIndices::kESDSumLocalChi2YPerDE,itc);
    TH1* hESDSumLocalChi2PerDE = Master()->GetESDsData(AliMUONQAIndices::kESDSumLocalChi2PerDE,itc);
    //  
    if (hESDnClustersPerCh && hESDnTotClustersPerCh) {
      hESDnClustersPerCh->Reset(); 
      hESDnClustersPerCh->Add(hESDnTotClustersPerCh, 1./nTracks);
    }
    if (hESDnClustersPerDE && hESDnTotClustersPerDE) {
      hESDnClustersPerDE->Reset();      
      hESDnClustersPerDE->Add(hESDnTotClustersPerDE, 1./nTracks);
    }
    //
    // loop over chambers
    for (Int_t iCh = 0; iCh < AliMUONConstants::NTrackingCh(); iCh++) {
      //
      double sigmaCharge=0,sigmaSize=0,sigmaResidualX=0,sigmaResidualY=0,sigmaLocalChi2X=0,sigmaLocalChi2Y=0,sigmaLocalChi2=0;
      //
      TH1* hESDClusterChargeInCh = Master()->GetESDsData(AliMUONQAIndices::kESDClusterChargeInCh+iCh,itc);
      if (hESDClusterChargeInCh) {
	sigmaCharge = hESDClusterChargeInCh->GetRMS();
	if (hESDClusterChargePerChMean) { 
	  hESDClusterChargePerChMean->SetBinContent(iCh+1, hESDClusterChargeInCh->GetMean());
	  hESDClusterChargePerChMean->SetBinError(iCh+1, hESDClusterChargeInCh->GetMeanError());
	}
	if (hESDClusterChargePerChSigma) {
	  hESDClusterChargePerChSigma->SetBinContent(iCh+1, sigmaCharge);
	  hESDClusterChargePerChSigma->SetBinError(iCh+1, hESDClusterChargeInCh->GetRMSError());
	}    
      }
      //
      TH1* hESDClusterSizeInCh = Master()->GetESDsData(AliMUONQAIndices::kESDClusterSizeInCh+iCh,itc);
      if (hESDClusterSizeInCh) {
	sigmaSize = hESDClusterSizeInCh->GetRMS();
	if (hESDClusterSizePerChMean) {
	  hESDClusterSizePerChMean->SetBinContent(iCh+1, hESDClusterSizeInCh->GetMean());
	  hESDClusterSizePerChMean->SetBinError(iCh+1, hESDClusterSizeInCh->GetMeanError());
	}
	if (hESDClusterSizePerChSigma) {
	  hESDClusterSizePerChSigma->SetBinContent(iCh+1, sigmaSize);
	  hESDClusterSizePerChSigma->SetBinError(iCh+1, hESDClusterSizeInCh->GetRMSError());
	}
      }    
      //
      TH1* hESDResidualXInCh = Master()->GetESDsData(AliMUONQAIndices::kESDResidualXInCh+iCh,itc);
      if (hESDResidualXInCh) {
	sigmaResidualX = hESDResidualXInCh->GetRMS();
	if (hESDResidualXPerChMean) {
	  hESDResidualXPerChMean->SetBinContent(iCh+1, hESDResidualXInCh->GetMean());
	  hESDResidualXPerChMean->SetBinError(iCh+1, hESDResidualXInCh->GetMeanError());
	}
	if (hESDResidualXPerChSigma) {
	  hESDResidualXPerChSigma->SetBinContent(iCh+1, sigmaResidualX);
	  hESDResidualXPerChSigma->SetBinError(iCh+1, hESDResidualXInCh->GetRMSError());
	}
      }
      //
      TH1* hESDResidualYInCh = Master()->GetESDsData(AliMUONQAIndices::kESDResidualYInCh+iCh,itc);
      if (hESDResidualYInCh) {
	sigmaResidualY = hESDResidualYInCh->GetRMS();
	if (hESDResidualYPerChMean) {
	  hESDResidualYPerChMean->SetBinContent(iCh+1, hESDResidualYInCh->GetMean());
	  hESDResidualYPerChMean->SetBinError(iCh+1, hESDResidualYInCh->GetMeanError());
	}
	if (hESDResidualYPerChSigma) {
	  hESDResidualYPerChSigma->SetBinContent(iCh+1, sigmaResidualY);
	  hESDResidualYPerChSigma->SetBinError(iCh+1, hESDResidualYInCh->GetRMSError());
	}
      }
      //      
      TH1* hESDLocalChi2XInCh = Master()->GetESDsData(AliMUONQAIndices::kESDLocalChi2XInCh+iCh,itc);
      if (hESDLocalChi2XInCh) {
	sigmaLocalChi2X = hESDLocalChi2XInCh->GetRMS();
	if (hESDLocalChi2XPerChMean) {
	  hESDLocalChi2XPerChMean->SetBinContent(iCh+1, hESDLocalChi2XInCh->GetMean());
	  hESDLocalChi2XPerChMean->SetBinError(iCh+1, hESDLocalChi2XInCh->GetMeanError());
	}
      }
      //
      TH1* hESDLocalChi2YInCh = Master()->GetESDsData(AliMUONQAIndices::kESDLocalChi2YInCh+iCh,itc);
      if (hESDLocalChi2YInCh) {
	sigmaLocalChi2Y = hESDLocalChi2YInCh->GetRMS();
	if (hESDLocalChi2YPerChMean) {
	  hESDLocalChi2YPerChMean->SetBinContent(iCh+1, hESDLocalChi2YInCh->GetMean());
	  hESDLocalChi2YPerChMean->SetBinError(iCh+1, hESDLocalChi2YInCh->GetMeanError());
	}    
      }
      //
      TH1* hESDLocalChi2InCh = Master()->GetESDsData(AliMUONQAIndices::kESDLocalChi2InCh+iCh,itc);
      if (hESDLocalChi2InCh) {
	sigmaLocalChi2 = hESDLocalChi2InCh->GetRMS();
	if (hESDLocalChi2PerChMean) {
	  hESDLocalChi2PerChMean->SetBinContent(iCh+1, hESDLocalChi2InCh->GetMean());
	  hESDLocalChi2PerChMean->SetBinError(iCh+1, hESDLocalChi2InCh->GetMeanError());
	}
      }
      //
      // loop over DE into chamber iCh
      AliMpDEIterator it;
      it.First(iCh);
      while ( !it.IsDone()) {
	
	Int_t iDE = it.CurrentDEId();
	
	Double_t nClusters = hESDnTotClustersPerDE ? hESDnTotClustersPerDE->GetBinContent(iDE+1) : 0;
	if (nClusters > 1) {
	  
	  if (hESDClusterChargePerDE && hESDSumClusterChargePerDE) {
	    hESDClusterChargePerDE->SetBinContent(iDE+1, hESDSumClusterChargePerDE->GetBinContent(iDE+1)/nClusters);
	    hESDClusterChargePerDE->SetBinError(iDE+1, sigmaCharge/TMath::Sqrt(nClusters));
	  }
        
	  if (hESDResidualXPerDEMean && hESDResidualXPerDESigma && hESDSumResidualXPerDE) {
	    Double_t meanResX = hESDSumResidualXPerDE->GetBinContent(iDE+1)/nClusters;
	    hESDResidualXPerDEMean->SetBinContent(iDE+1, meanResX);
	    hESDResidualXPerDEMean->SetBinError(iDE+1, sigmaResidualX/TMath::Sqrt(nClusters));
	    //
	    hESDResidualXPerDESigma->SetBinContent(iDE+1, ProtectedSqrt(hESDSumResidualX2PerDE->GetBinContent(iDE+1)/nClusters - meanResX*meanResX));
	    hESDResidualXPerDESigma->SetBinError(iDE+1, sigmaResidualX/TMath::Sqrt(2.*nClusters));
	  }
	  //        
	  if (hESDResidualYPerDEMean && hESDResidualYPerDESigma && hESDSumResidualYPerDE) {
	    Double_t meanResY = hESDSumResidualYPerDE->GetBinContent(iDE+1)/nClusters;
	    hESDResidualYPerDEMean->SetBinContent(iDE+1, meanResY);
	    hESDResidualYPerDEMean->SetBinError(iDE+1, sigmaResidualY/TMath::Sqrt(nClusters));
	    hESDResidualYPerDESigma->SetBinContent(iDE+1, ProtectedSqrt(hESDSumResidualY2PerDE->GetBinContent(iDE+1)/nClusters - meanResY*meanResY));
	    hESDResidualYPerDESigma->SetBinError(iDE+1, sigmaResidualY/TMath::Sqrt(2.*nClusters));
	  }
	  //
	  if (hESDLocalChi2XPerDEMean && hESDSumLocalChi2XPerDE) {
	    hESDLocalChi2XPerDEMean->SetBinContent(iDE+1, hESDSumLocalChi2XPerDE->GetBinContent(iDE+1)/nClusters);
	    hESDLocalChi2XPerDEMean->SetBinError(iDE+1, sigmaLocalChi2X/TMath::Sqrt(nClusters));
	  }
	  //
	  if (hESDLocalChi2YPerDEMean && hESDSumLocalChi2YPerDE) {
	    hESDLocalChi2YPerDEMean->SetBinContent(iDE+1, hESDSumLocalChi2YPerDE->GetBinContent(iDE+1)/nClusters);
	    hESDLocalChi2YPerDEMean->SetBinError(iDE+1, sigmaLocalChi2Y/TMath::Sqrt(nClusters));
	  }
	  //
	  if (hESDLocalChi2PerDEMean && hESDSumLocalChi2PerDE) {
	    hESDLocalChi2PerDEMean->SetBinContent(iDE+1, hESDSumLocalChi2PerDE->GetBinContent(iDE+1)/nClusters);
	    hESDLocalChi2PerDEMean->SetBinError(iDE+1, sigmaLocalChi2/TMath::Sqrt(nClusters));
	  }
	} 
	else {
	  if (hESDClusterChargePerDE && hESDClusterChargeInCh && hESDSumClusterChargePerDE) {
	    hESDClusterChargePerDE->SetBinContent(iDE+1, hESDSumClusterChargePerDE->GetBinContent(iDE+1));
	    hESDClusterChargePerDE->SetBinError(iDE+1, hESDClusterChargeInCh->GetXaxis()->GetXmax());
	  }
	  //
	  if (hESDResidualXPerDEMean && hESDResidualXPerDESigma && hESDSumResidualXPerDE && hESDResidualXInCh) {
	    hESDResidualXPerDEMean->SetBinContent(iDE+1, hESDSumResidualXPerDE->GetBinContent(iDE+1));
	    hESDResidualXPerDEMean->SetBinError(iDE+1, hESDResidualXInCh->GetXaxis()->GetXmax());
	    hESDResidualXPerDESigma->SetBinContent(iDE+1, 0.);
	    hESDResidualXPerDESigma->SetBinError(iDE+1, hESDResidualXInCh->GetXaxis()->GetXmax());
	  }
	  //
	  if (hESDResidualYPerDEMean && hESDResidualYPerDESigma && hESDSumResidualYPerDE && hESDResidualYInCh) {
	    hESDResidualYPerDEMean->SetBinContent(iDE+1, hESDSumResidualYPerDE->GetBinContent(iDE+1));
	    hESDResidualYPerDEMean->SetBinError(iDE+1, hESDResidualYInCh->GetXaxis()->GetXmax());
	    hESDResidualYPerDESigma->SetBinContent(iDE+1, 0.);
	    hESDResidualYPerDESigma->SetBinError(iDE+1, hESDResidualYInCh->GetXaxis()->GetXmax());
	  }
	  //
	  if (hESDLocalChi2XPerDEMean && hESDSumLocalChi2XPerDE && hESDLocalChi2XInCh) {        
	    hESDLocalChi2XPerDEMean->SetBinContent(iDE+1, hESDSumLocalChi2XPerDE->GetBinContent(iDE+1));
	    hESDLocalChi2XPerDEMean->SetBinError(iDE+1, hESDLocalChi2XInCh->GetXaxis()->GetXmax());
	  }
	  //
	  if (hESDLocalChi2YPerDEMean && hESDSumLocalChi2YPerDE && hESDLocalChi2YInCh) {
	    hESDLocalChi2YPerDEMean->SetBinContent(iDE+1, hESDSumLocalChi2YPerDE->GetBinContent(iDE+1));
	    hESDLocalChi2YPerDEMean->SetBinError(iDE+1, hESDLocalChi2YInCh->GetXaxis()->GetXmax());
	  }
	  //
	  if (hESDLocalChi2PerDEMean && hESDSumLocalChi2PerDE && hESDLocalChi2InCh) {
	    hESDLocalChi2PerDEMean->SetBinContent(iDE+1, hESDSumLocalChi2PerDE->GetBinContent(iDE+1));
	    hESDLocalChi2PerDEMean->SetBinError(iDE+1, hESDLocalChi2InCh->GetXaxis()->GetXmax());
	  }
	}
      
	Double_t nFullClusters = hESDnTotFullClustersPerDE ? hESDnTotFullClustersPerDE->GetBinContent(iDE+1) : 0;
	if (nFullClusters > 1) {
	  if (hESDClusterSizePerDE && hESDSumClusterSizePerDE) {        
	    hESDClusterSizePerDE->SetBinContent(iDE+1, hESDSumClusterSizePerDE->GetBinContent(iDE+1)/nFullClusters);
	    hESDClusterSizePerDE->SetBinError(iDE+1, sigmaSize/TMath::Sqrt(nFullClusters));
	  }
	}
	else {
	  if (hESDClusterSizePerDE && hESDSumClusterSizePerDE) {
	    hESDClusterSizePerDE->SetBinContent(iDE+1, hESDSumClusterSizePerDE->GetBinContent(iDE+1));
	    hESDClusterSizePerDE->SetBinError(iDE+1, hESDClusterSizeInCh->GetXaxis()->GetXmax());
	  }
	}
	
	it.Next();
      } // while
      
    } // iCh loop

  } // trigger classes loop

}
  
//____________________________________________________________________________ 
void AliMUONTrackerQADataMakerRec::EndOfDetectorCycleRecPoints(Int_t specie, TObjArray** list)
{
  /// Normalize RecPoints histograms
  //
  Bool_t firstFill=kTRUE, needIns0=kTRUE, needIns1=kTRUE;
  //
  for (int itc=-1;itc<Master()->GetNTrigClasses();itc++) { // RS: loop over eventual clones per trigger class 

    TH1* hTrackerClusterChargePerChMean = GetRecPointsData(AliMUONQAIndices::kTrackerClusterChargePerChMean,itc);
    if (!hTrackerClusterChargePerChMean) continue;
    //
    if (firstFill) { AliCodeTimerAuto("",0); firstFill = kFALSE;}

    TH1* hTrackerClusterChargePerChSigma = GetRecPointsData(AliMUONQAIndices::kTrackerClusterChargePerChSigma,itc);
    TH1* hTrackerClusterMultiplicityPerChMean = GetRecPointsData(AliMUONQAIndices::kTrackerClusterMultiplicityPerChMean,itc);
    TH1* hTrackerClusterMultiplicityPerChSigma = GetRecPointsData(AliMUONQAIndices::kTrackerClusterMultiplicityPerChSigma,itc);
    TH1* hTrackerClusterChargePerDEMean = GetRecPointsData(AliMUONQAIndices::kTrackerClusterChargePerDEMean,itc);
    TH1* hTrackerClusterMultiplicityPerDEMean = GetRecPointsData(AliMUONQAIndices::kTrackerClusterMultiplicityPerDEMean,itc);
  
    // loop over chambers
    for (Int_t iCh = 0; iCh < AliMUONConstants::NTrackingCh(); iCh++) {
    
      TH1* hTrackerClusterChargePerChamber = GetRecPointsData(AliMUONQAIndices::kTrackerClusterChargePerChamber+iCh,itc);
      Double_t sigmaCharge = 0;
      if (hTrackerClusterChargePerChamber) {
	sigmaCharge = hTrackerClusterChargePerChamber->GetRMS();
	hTrackerClusterChargePerChMean->SetBinContent(iCh+1, hTrackerClusterChargePerChamber->GetMean());
	hTrackerClusterChargePerChMean->SetBinError(iCh+1,   hTrackerClusterChargePerChamber->GetMeanError());
	//
	if (hTrackerClusterChargePerChSigma) {
	  hTrackerClusterChargePerChSigma->SetBinContent(iCh+1, sigmaCharge);
	  hTrackerClusterChargePerChSigma->SetBinError(iCh+1, hTrackerClusterChargePerChamber->GetRMSError());
	}
      }
      //
      TH1* hTrackerClusterMultiplicityPerChamber = GetRecPointsData(AliMUONQAIndices::kTrackerClusterMultiplicityPerChamber+iCh,itc);
      Double_t sigmaSize = 0;
      if (hTrackerClusterMultiplicityPerChamber) {
	sigmaSize = hTrackerClusterMultiplicityPerChamber->GetRMS();
	if (hTrackerClusterMultiplicityPerChMean) {
	  hTrackerClusterMultiplicityPerChMean->SetBinContent(iCh+1, hTrackerClusterMultiplicityPerChamber->GetMean());
	  hTrackerClusterMultiplicityPerChMean->SetBinError(iCh+1, hTrackerClusterMultiplicityPerChamber->GetMeanError());
	}
	if (hTrackerClusterMultiplicityPerChSigma) {
	  hTrackerClusterMultiplicityPerChSigma->SetBinContent(iCh+1, sigmaSize);
	  hTrackerClusterMultiplicityPerChSigma->SetBinError(iCh+1, hTrackerClusterMultiplicityPerChamber->GetRMSError());
	}
      }
      //
      // loop over DE into chamber iCh
      AliMpDEIterator it;
      it.First(iCh);
      while ( !it.IsDone()) {
      
	Int_t iDE = it.CurrentDEId();
      
	TH1* hTrackerClusterChargePerDE = GetRecPointsData(AliMUONQAIndices::kTrackerClusterChargePerDE+iDE,itc);
	if (hTrackerClusterChargePerDEMean && hTrackerClusterChargePerDE) {
	  hTrackerClusterChargePerDEMean->SetBinContent(iDE+1, hTrackerClusterChargePerDE->GetMean());
	  Double_t nClusters = hTrackerClusterChargePerDE->GetEntries();
	  if (nClusters > 1) hTrackerClusterChargePerDEMean->SetBinError(iDE+1, sigmaCharge/TMath::Sqrt(nClusters));
	  else hTrackerClusterChargePerDEMean->SetBinError(iDE+1, hTrackerClusterChargePerChamber ? 
							   hTrackerClusterChargePerChamber->GetXaxis()->GetXmax() : 0);
	}
	TH1* hTrackerClusterMultiplicityPerDE = GetRecPointsData(AliMUONQAIndices::kTrackerClusterMultiplicityPerDE+iDE,itc);
	if (hTrackerClusterMultiplicityPerDEMean && hTrackerClusterMultiplicityPerDE) {
	  hTrackerClusterMultiplicityPerDEMean->SetBinContent(iDE+1, hTrackerClusterMultiplicityPerDE->GetMean());
	  Double_t nClusters = hTrackerClusterMultiplicityPerDE->GetEntries();
	  if (nClusters > 1) hTrackerClusterMultiplicityPerDEMean->SetBinError(iDE+1, sigmaSize/TMath::Sqrt(nClusters));
	  else hTrackerClusterMultiplicityPerDEMean->SetBinError(iDE+1, hTrackerClusterMultiplicityPerChamber ? 
								 hTrackerClusterMultiplicityPerChamber->GetXaxis()->GetXmax() : 0);
	}
	it.Next();
      }
    }
  
    if (needIns0 && MappingCheckRecPoints(specie) ) { // RS: I guess this should not be in the itc loop, do this only once
      InsertTrackerData(specie,list,MappingCheckRecPoints(specie)->CreateData("RecPoints"),AliMUONQAIndices::kTrackerRecPoints,kTRUE);    
      needIns0 = kFALSE;
    }
    
    if ( TrackerRecData(specie) ) {
      /// put the trackerdata in the pipeline
      if (needIns1) { // RS: I guess this should not be in the itc loop, do this only once
	InsertTrackerData(specie,list,TrackerRecData(specie),AliMUONQAIndices::kTrackerData);
	needIns1 = kFALSE;
      }
      TH1* hbp = GetRecPointsData(AliMUONQAIndices::kTrackerBusPatchOccupancy,itc);
      TH1* hnevents = GetRecPointsData(AliMUONQAIndices::kTrackerNofGoodPhysicsEventsUsed,itc);
      TH1* hddl = GetRecPointsData(AliMUONQAIndices::kTrackerDDLOccupancy,itc);
      TH1* hddlevents = GetRecPointsData(AliMUONQAIndices::kTrackerDDLNofEventsUsed,itc);
    
      if (itc==-1 && (!hbp || !hnevents || !hddl || !hddlevents)) { //RS: produce error only for trigger-blind class
	AliError(Form("Missing some histograms : cannot work : hbp=%p hnevents=%p hddl=%p hddlevents=%p",hbp,hnevents,hddl,hddlevents));
	continue; // return; // RS: the histos might be booked for some trigger class only
      }
      //
      ProjectTrackerData(TrackerRecData(specie),*hbp,*hnevents,*hddl,*hddlevents);    
    }
    else {
      AliError(Form("TrackerRecData is null for specie %s",AliRecoParam::GetEventSpecieName(specie)));
    }
  } //  RS: loop over eventual clones per trigger class
}

//____________________________________________________________________________ 
void AliMUONTrackerQADataMakerRec::EndOfDetectorCycleDigits(Int_t specie, TObjArray** list)
{
  /// create digits histograms in digits subdir
  
  if ( TrackerCalData(specie) )
  {
    AliCodeTimerAuto("",0);
    
    /// put the trackerdata in the pipeline
    InsertTrackerData(specie,list,TrackerCalData(specie),AliMUONQAIndices::kTrackerData);
    //
    for (int itc=-1;itc<Master()->GetNTrigClasses();itc++) { // RS: loop over eventual clones per trigger class 
      //    
      TH1* hbp = GetDigitsData(AliMUONQAIndices::kTrackerBusPatchOccupancy, itc);
      TH1* hnevents = GetDigitsData(AliMUONQAIndices::kTrackerNofGoodPhysicsEventsUsed, itc);
      TH1* hddl = GetDigitsData(AliMUONQAIndices::kTrackerDDLOccupancy, itc);
      TH1* hddlevents = GetDigitsData(AliMUONQAIndices::kTrackerDDLNofEventsUsed, itc);
    
      if ( (!hbp || !hnevents || !hddl || !hddlevents) ) 
      { 
        if (itc==-1)
        {
          // report error only for trigger-blind class
          AliError(Form("Missing some histograms : cannot work : hbp=%p hnevents=%p hddl=%p hddlevents=%p",hbp,hnevents,hddl,hddlevents));
          continue; //return; // RS
        }
        else
        {
          continue;
        }
      }
      //    
      ProjectTrackerData(TrackerCalData(specie), *hbp,*hnevents,*hddl,*hddlevents);    
    } //  RS: loop over eventual clones per trigger class 
  }
}

//____________________________________________________________________________ 
AliMUONQADataMakerRec* AliMUONTrackerQADataMakerRec::Master() const
{
  /// Get our master
  return static_cast<AliMUONQADataMakerRec*>(fMaster); 
}

//____________________________________________________________________________ 
void AliMUONTrackerQADataMakerRec::ProjectTrackerData(AliMUONVTrackerData* data, 
                                                      TH1& hbp,
                                                      TH1& hnevents,
                                                      TH1& hddl,
                                                      TH1& hddlevents)
{
  /// Project tracker data into shifter-friendly histograms
  
  AliCodeTimerAuto(Form("%s",data->Name()),0);
  
  /// project the tracerdata into buspatch occupancies (for the experts)
  hbp.Reset();
  hbp.SetStats(kFALSE);
  
  TIter nextBP(AliMpDDLStore::Instance()->CreateBusPatchIterator());
  AliMpBusPatch* bp(0x0);
  Int_t occDim = 2;
  
  while ( ( bp = static_cast<AliMpBusPatch*>(nextBP())) )
  {
    Int_t busPatchId = bp->GetId();
    Int_t bin = hbp.FindBin(busPatchId);
    if ( data->HasBusPatch(busPatchId) )
    {
      hbp.SetBinContent(bin,data->BusPatch(busPatchId,occDim)*100.0); // occupancy, in percent
    }
  }
  
  /// log the readout errors (for the shifter)
  hnevents.Reset();
  hnevents.SetStats(kFALSE);
  hnevents.SetBinContent(1,data->NumberOfEvents(-1));
  
  /// project tracker data into DDL occupancies (for the shifter)
  hddl.Reset();
  hddl.SetStats(kFALSE);
  hddlevents.Reset();
  hddlevents.SetStats(kFALSE);
  
  const Int_t nddls = AliDAQ::NumberOfDdls("MUONTRK");
  const Int_t offset = AliDAQ::DdlIDOffset("MUONTRK");
  
  for ( Int_t iddl = 0; iddl < nddls; ++iddl )
  {
    AliMpDDL* ddl = AliMpDDLStore::Instance()->GetDDL(iddl);
    
    Int_t ddlId = offset + ddl->GetId();
    Int_t npads = 0;
    
    Int_t nevents = data->NumberOfEvents(iddl);
    
    hddlevents.Fill(ddlId,nevents);
    
    Double_t occ(0.0);
    
    if ( nevents > 0 )
    {
      for ( Int_t ide = 0; ide < ddl->GetNofDEs(); ++ide )
      {
        Int_t de = ddl->GetDEId(ide);
        
        npads += TMath::Nint(data->DetectionElement(de,3));
        
        occ +=  data->DetectionElement(de,4);
      }
      
      if ( npads > 0 )
      {
        occ = occ/npads/nevents;
      }
      else 
      {
        occ = 0.0;
      }
    }
    
    hddl.Fill(ddlId,100.0*occ); // occ in percent
  }
  
  TPaveText* text = new TPaveText(0.50,0.8,0.9,0.9,"NDC");
  
  text->AddText(Form("MCH RUN %d - %d events",AliCDBManager::Instance()->GetRun(),data->NumberOfEvents(-1)));

  hddl.GetListOfFunctions()->Add(text);
}

//____________________________________________________________________________ 
void AliMUONTrackerQADataMakerRec::EndOfDetectorCycleRaws(Int_t specie, TObjArray** list)
{
  /// create Raws histograms in Raws subdir
  
  if ( TrackerDataMaker(specie) && TrackerDataMaker(specie)->Data() ) 
  {
    AliCodeTimerAuto("",0);
    
    /// put the trackerdata in the pipeline
    InsertTrackerData(specie,list,TrackerDataMaker(specie)->Data(),AliMUONQAIndices::kTrackerData);
    TrackerDataMaker(specie)->SetOwnerOfData(kFALSE); // now that it's attached to list, list will take care of the deletion
    //    
    for (int itc=-1;itc<Master()->GetNTrigClasses();itc++) { // RS: loop over eventual clones per trigger class

      TH1* hbp = GetRawsData(AliMUONQAIndices::kTrackerBusPatchOccupancy, itc);
      TH1* hnevents = GetRawsData(AliMUONQAIndices::kTrackerNofGoodPhysicsEventsUsed, itc);
      TH1* hddl = GetRawsData(AliMUONQAIndices::kTrackerDDLOccupancy, itc);
      TH1* hddlevents = GetRawsData(AliMUONQAIndices::kTrackerDDLNofEventsUsed, itc);
      //
      if (!hbp || !hnevents || !hddl || !hddlevents) {
	if (itc==-1) AliError(Form("Missing some histograms : cannot work : hbp=%p hnevents=%p hddl=%p hddlevents=%p",hbp,hnevents,hddl,hddlevents));
	continue; // return; // RS
      }
      //
      ProjectTrackerData(TrackerDataMaker(specie)->Data(), *hbp,*hnevents,*hddl,*hddlevents);
      //
      FillReadoutStatus(*fLogger,TrackerDataMaker(specie)->Data(), itc);      
    } // RS: loop over eventual clones per trigger class
    //
  }    
}

//____________________________________________________________________________ 
void AliMUONTrackerQADataMakerRec::FillReadoutStatus(AliMUONLogger& log, AliMUONVTrackerData* data, Int_t trigCl)
{
  // RS: I am not sure if this part is valid for looping over trigger classes (loggers are not cloned)
  if (trigCl!=-1) return; // For the moment only trigCl==-1 will be processed (i.e. trigger-blind histos)

  log.ResetItr();

  TString msg;
  Int_t occurence;
  
  TH1* hparity = GetRawsData(AliMUONQAIndices::kTrackerBusPatchParityErrors, trigCl);

  TH1* htoken = GetRawsData(AliMUONQAIndices::kTrackerBusPatchTokenLostErrors, trigCl);

  TH1* hpadding = GetRawsData(AliMUONQAIndices::kTrackerBusPatchPaddingErrors, trigCl);
  
  TH1* hrostatus = GetRawsData(AliMUONQAIndices::kTrackerReadoutStatus, trigCl);
    
  TH1* hnevents = GetRawsData(AliMUONQAIndices::kTrackerNofPhysicsEventsSeen, trigCl);

  //
  if (!hparity || !htoken || !hpadding || !hrostatus || !hnevents) return;
  
  Int_t nevents = TMath::Nint(hnevents->GetBinContent(1));
  
  if ( !nevents ) 
  {
    TPaveText* text = new TPaveText(0,0,0.99,0.99,"NDC");
    text->AddText("FATAL : 0 event seen ? That's NOT normal...");
    text->SetFillColor(2); // red = FATAL
    hrostatus->GetListOfFunctions()->Delete();
    hrostatus->GetListOfFunctions()->Add(text);
    return;
  }
  
  /////////////////////////////////////////////////////////////////
  /// Start by counting the number of errors
  /////////////////////////////////////////////////////////////////
  
  while ( log.Next(msg,occurence) )
  {
    AliDebug(1,Form("msg=%s occurence=%d",msg.Data(),occurence));
             
    if ( msg.Contains("token") )
    {
      Int_t dsp(-1),iddl(-1),ecode(-1);
      
      sscanf(msg.Data(),"Lost token error detected with address 0x%10X of DDL %10d and code %10d.",
             &dsp,&iddl,&ecode);
      Int_t localBP = ((dsp >> 16)- 4)*5 + 1;
      Int_t buspatch = localBP + iddl*100;
      
      // Let's try to get all the suspected bus patches (i.e. one full FRT, as currently
      // we don't have more precise information to locate the faulty bus patch(es)).
      
      AliMpBusPatch* bp = AliMpDDLStore::Instance()->GetBusPatch(buspatch,kFALSE);
      if (bp)
      {
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
      }
      hrostatus->Fill(1.0*AliMUONQAIndices::kTrackerRawNofTokenLostErrors,occurence);
    }
    
    if ( msg.Contains("Parity") )
    {
      Int_t buspatch;
      sscanf(msg.Data(),"Parity error in buspatch %d (0x%X).",&buspatch,&buspatch);
      hparity->Fill(buspatch,occurence);      
      hrostatus->Fill(1.0*AliMUONQAIndices::kTrackerRawNofParityErrors,occurence);
    }
    
    if ( msg.Contains("Glitch") ) 
    {
      hrostatus->Fill(1.0*AliMUONQAIndices::kTrackerRawNofGlitchErrors,occurence);      
    }
    
    if ( msg.Contains("Padding") )
    {
      Int_t block, dsp, buspatch;      
      sscanf(msg.Data(),"Padding word error for iBlock %d, iDsp %d, iBus %d.",&block,&dsp,&buspatch);
      hpadding->Fill(buspatch,occurence);
      hrostatus->Fill(1.0*AliMUONQAIndices::kTrackerRawNofPaddingErrors,occurence);      
    }
  }
  
  /////////////////////////////////////////////////////////////////
  ///
  /// Then make the status about number of missing bus patches
  ///
  /////////////////////////////////////////////////////////////////

  Int_t nofBusPatchesNotInConfig(0);
    
  for ( int i = 1; i <= fBusPatchConfig->GetNbinsX(); ++i )
  {
    Double_t buspatchId = fBusPatchConfig->GetBinCenter(i);
    if ( TMath::Nint(buspatchId) != i ) 
    {
      AliError(Form("buspathId=%e i=%d",buspatchId,i));
    }
    Double_t content = fBusPatchConfig->GetBinContent(i);
    
    if ( content <= 0. /* no content */
        && 
        AliMpDDLStore::Instance()->GetBusPatch(i,kFALSE) /* but a valid bus patch */ )
    {
      ++nofBusPatchesNotInConfig;
    }    
  }
  
  Double_t nbuspatches = fBusPatchConfig->GetEntries();
  
  Int_t bin = hrostatus->FindBin(1.0*AliMUONQAIndices::kTrackerRawNofMissingBusPatchesFromConfig);
  hrostatus->SetBinContent(bin,nofBusPatchesNotInConfig*nevents/nbuspatches);
  
  Double_t nofBusPatchesNotInData(0);
  
  TIter next(AliMpDDLStore::Instance()->CreateBusPatchIterator());
  AliMpBusPatch* bp;
  
  while ( ( bp = static_cast<AliMpBusPatch*>(next()) ) )
  {
    if ( !data->HasBusPatch(bp->GetId()) ) ++nofBusPatchesNotInData;
  }
  
  bin = hrostatus->FindBin(1.0*AliMUONQAIndices::kTrackerRawNofMissingBusPatchesFromDataStream);
  hrostatus->SetBinContent(bin,nofBusPatchesNotInData*nevents/nbuspatches);
}

//____________________________________________________________________________ 
void AliMUONTrackerQADataMakerRec::FillEventSize(AliRawReader* rawReader)
{
  /// Fill event size histogram(s)
  // RS: substituted explicit histo filling by QA framework filling accounting for cloned histos
  
  FillRawsData(AliMUONQAIndices::kTrackerNofPhysicsEventsSeen,0.0);
  
  Int_t offset = AliDAQ::DdlIDOffset("MUONTRK");
  
  for ( int i = 0; i < AliDAQ::NumberOfDdls("MUONTRK"); ++i )
  {
    rawReader->Reset();
    rawReader->Select("MUONTRK",i,i);
    if (rawReader->ReadHeader() )
    {
      UInt_t ddlsize = rawReader->GetEquipmentSize();
      FillRawsData(AliMUONQAIndices::kTrackerDDLEventSize,i+offset,ddlsize);
      FillRawsData(AliMUONQAIndices::kTrackerDDLNofEventsSeen,i+offset);
    }      
  }
  rawReader->Reset();
}

//____________________________________________________________________________ 
void AliMUONTrackerQADataMakerRec::InitCommon()
{  
  if (!fBusPatchConfig)
  {
    Int_t bpmin(999999);
    Int_t bpmax(0);
    
    TIter next(AliMpDDLStore::Instance()->CreateBusPatchIterator());
    AliMpBusPatch* bp(0x0);
    while ( ( bp = static_cast<AliMpBusPatch*>(next())) )
    {
      bpmin = TMath::Min(bpmin,bp->GetId());
      bpmax = TMath::Max(bpmax,bp->GetId());
    }
    
    fBPxmin = bpmin-0.5;
    fBPxmax = bpmax+0.5;
    fBPnbins = TMath::Nint(fBPxmax-fBPxmin);
        
    AliMUONVStore* config = fCalibrationData->Config();
    
    if (config)
    {
      fBusPatchConfig = new TH1F("hTrackerBusPatchConfig","Configuration of bus patches",fBPnbins,fBPxmin,fBPxmax);
      fBusPatchConfig->SetDirectory(0);
    }
    else
    {
      AliWarning("Tracker configuration not found. Will not be able to cut on low occupancies");
    }
    
    next.Reset();
    while ( ( bp = static_cast<AliMpBusPatch*>(next())) )
    {
      if ( config ) 
      {
        Int_t nofManusInConfig(0);
      
        for ( Int_t imanu = 0; imanu < bp->GetNofManus(); ++imanu )
        {
          Int_t manuId = bp->GetManuId(imanu);
          if ( config->FindObject(bp->GetDEId(),manuId)) ++nofManusInConfig;
        }
        
        if ( nofManusInConfig > 0 )
        {
          fBusPatchConfig->Fill(bp->GetId(),1.0);
        }
        else
        {
          fBusPatchConfig->Fill(bp->GetId(),0.0);          
        }
      }      
      else // no config, we assume all is there...
      {
        fBusPatchConfig->Fill(bp->GetId());
      }
    }
  }
}

//____________________________________________________________________________ 
void AliMUONTrackerQADataMakerRec::BookHistograms(AliQAv1::TASKINDEX_t task)
{
  AliCodeTimerAuto("",0);

  InitCommon();
  
  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1* hbp = new TH1F("hTrackerBusPatchOccupancy","Occupancy of bus patches",fBPnbins,fBPxmin,fBPxmax);
  
  Master()->Add2List(hbp,AliMUONQAIndices::kTrackerBusPatchOccupancy, task, expert, !image, !saveCorr);
  
  TH1* h = new TH1F("hTrackerBusPatchParityErrors","Number of parity errors per bus patch",fBPnbins,fBPxmin,fBPxmax);
  
  Master()->Add2List(h,AliMUONQAIndices::kTrackerBusPatchParityErrors,task,expert,!image,!saveCorr);
  
  h = new TH1F("hTrackerBusPatchTokenLostErrors","Number of token lost errors per bus patch",fBPnbins,fBPxmin,fBPxmax);
  Master()->Add2List(h,AliMUONQAIndices::kTrackerBusPatchTokenLostErrors,task,expert,!image,!saveCorr);
  
  h = new TH1F("hTrackerBusPatchPaddingErrors","Number of padding errors per bus patch",fBPnbins,fBPxmin,fBPxmax);
  
  Master()->Add2List(h,AliMUONQAIndices::kTrackerBusPatchPaddingErrors,task,expert,!image,!saveCorr);
  
  
  TH1* hnevents(0x0);
  
  if ( task == AliQAv1::kRAWS )
  {
    // for raw data, we differentiate events seen from events used to be able to detect
    // severe decoder errors that lead to no event decoded (i.e. zero event used) even if
    // events are there (i.e non-zero event seen).
    hnevents = new TH1F("hTrackerNofPhysicsEventsSeen","Number of physics events seen",1,-0.5,0.5);
    // this one will count the number of physics event the rawdatamaker is *seeing*
    TAxis* a = hnevents->GetXaxis();
    a->SetBinLabel(1,"NPhysicsEvents");
    hnevents->SetStats(kFALSE);
    Master()->Add2List(hnevents,AliMUONQAIndices::kTrackerNofPhysicsEventsSeen,task,expert,!image,!saveCorr);
  }
  
  hnevents = new TH1F("hTrackerNofGoodPhysicsEventsUsed","Number of good physics events used",1,-0.5,0.5);
  // this one will get its content from the TrackerData, i.e. it will count the number of *good* physics events *used*
  // (i.e. not empty and with no fatal readout error)
  TAxis* a = hnevents->GetXaxis();
  a->SetBinLabel(1,"NGoodPhysicsEvents");
  hnevents->SetStats(kFALSE);  

  Master()->Add2List(hnevents,AliMUONQAIndices::kTrackerNofGoodPhysicsEventsUsed,task,expert,!image,!saveCorr);

  Master()->Add2List(static_cast<TH1*>(fBusPatchConfig->Clone()),AliMUONQAIndices::kTrackerBusPatchConfig, task,expert, !image, !saveCorr);

  Int_t nbins = AliDAQ::NumberOfDdls("MUONTRK");
  const Int_t offset = AliDAQ::DdlIDOffset("MUONTRK");
  
  Double_t xmin = offset - 0.5;
  Double_t xmax  = offset + nbins - 0.5;
  
  TString what(AliQAv1::GetTaskName(task));
  
  h = new TH1F(Form("hTrackerDDL%sOccupancy",what.Data()),Form(";DDLId;DDL Occupancy in %% (from %s)",what.Data()),nbins,xmin,xmax);
  
  Master()->Add2List(h,AliMUONQAIndices::kTrackerDDLOccupancy,task,expert,!image,!saveCorr);

  if ( task == AliQAv1::kRAWS )
  {
    // see above the comment about why we have event seen vs used for raw data.
    h = new TH1F("hTrackerDDLNofEventsSeen","Number of events seen by DDL;DDLId",nbins,xmin,xmax);
    Master()->Add2List(h,AliMUONQAIndices::kTrackerDDLNofEventsSeen,task,expert,!image,!saveCorr);
  }
  
  h = new TH1F("hTrackerDDLNofEventsUsed","Number of events used by DDL;DDLId",nbins,xmin,xmax);
  Master()->Add2List(h,AliMUONQAIndices::kTrackerDDLNofEventsUsed,task,expert,!image,!saveCorr);
  
}

//____________________________________________________________________________ 
void AliMUONTrackerQADataMakerRec::InitRaws()
{
  /// create monitor objects for RAWS
	
  TrackerDataMaker(AliRecoParam::AConvert(Master()->GetEventSpecie()),kTRUE);
  
  /// Book histograms that are common to Raws and Digits
  BookHistograms(AliQAv1::kRAWS);
  
  /// Now the Raws specific parts
  TH1* h = new TH1F("hTrackerReadoutStatus","Readout status (x events)",7,-0.5,6.5);
  h->SetStats(kFALSE);
  
  TAxis* a = h->GetXaxis();
  
  a->SetBinLabel(h->FindBin(1.0*AliMUONQAIndices::kTrackerRawNofGlitchErrors),"Glitch errors");
  a->SetBinLabel(h->FindBin(1.0*AliMUONQAIndices::kTrackerRawNofTokenLostErrors),"Token lost errors");
  a->SetBinLabel(h->FindBin(1.0*AliMUONQAIndices::kTrackerRawNofParityErrors),"Parity errors");
  a->SetBinLabel(h->FindBin(1.0*AliMUONQAIndices::kTrackerRawNofPaddingErrors),"Padding errors");
  
  a->SetBinLabel(h->FindBin(1.0*AliMUONQAIndices::kTrackerRawNofEmptyEvents),"Empty events");
  
  a->SetBinLabel(h->FindBin(1.0*AliMUONQAIndices::kTrackerRawNofMissingBusPatchesFromConfig),"Not readout bus patches");
  a->SetBinLabel(h->FindBin(1.0*AliMUONQAIndices::kTrackerRawNofMissingBusPatchesFromDataStream),"Missing bus patches");

  TH1* h1 = static_cast<TH1*>(h->Clone("hTrackerReadoutStatusPerEvent"));
  h1->SetTitle("Readout status per event");
  h1->GetYaxis()->SetTitle("Percentage");
  
  // The QA shifter will only see the summary plot below

  Add2RawsList(h,AliMUONQAIndices::kTrackerReadoutStatus,kTRUE,kFALSE,kFALSE);    
  Add2RawsList(h1,AliMUONQAIndices::kTrackerReadoutStatusPerEvent,kFALSE,kTRUE,kFALSE);  

  // Lastly the event size histograms
  
  Int_t nbins = AliDAQ::NumberOfDdls("MUONTRK");
  const Int_t offset = AliDAQ::DdlIDOffset("MUONTRK");
  
  Double_t xmin = offset - 0.5;
  Double_t xmax  = offset + nbins - 0.5;
  
  h = new TH1F("hTrackerDDLEventSize","DDL event size (bytes);DDL Id;Data size (bytes)",nbins,xmin,xmax);  
  h->SetStats(kFALSE);
  Add2RawsList(h,AliMUONQAIndices::kTrackerDDLEventSize,kTRUE,kFALSE,kFALSE);

  h = new TH1F("hTrackerDDLMeanEventSize","DDL mean event size (KB) per event;DDL Id;Mean Event size (KB)",nbins,xmin,xmax);  
  h->SetStats(kFALSE);
  Add2RawsList(h,AliMUONQAIndices::kTrackerDDLEventSizePerEvent,kFALSE,kTRUE,kFALSE);
    
  Add2RawsList(new TH1F("hTrackerIsThere","tracker is there",1,0,1),AliMUONQAIndices::kTrackerIsThere,kTRUE,kFALSE,kFALSE);
  //
  //ClonePerTrigClass(AliQAv1::kRAWS); // this should be the last line  DONE at parent level
}

//__________________________________________________________________
void AliMUONTrackerQADataMakerRec::InitDigits() 
{
  /// create monitor objects for DIGITS
  
  AliCodeTimerAuto("",0);

  if ( GetRecoParam()->TryRecover() )
  {
    fDigitMaker->SetTryRecover(kTRUE);
  }
  else
  {
    fDigitMaker->SetTryRecover(kFALSE);    
  }
  
  TrackerCalData(AliRecoParam::AConvert(Master()->GetEventSpecie()),kTRUE);
  
  /// Book histograms that are common to Raws and Digits
  BookHistograms(AliQAv1::kDIGITSR);
  //
  //ClonePerTrigClass(AliQAv1::kDIGITS); // this should be the last line  DONE at parent level
} 

//____________________________________________________________________________ 
void AliMUONTrackerQADataMakerRec::InitRecPoints()
{
  /// create Reconstructed Points histograms in RecPoints subdir for the
  /// MUON tracker subsystem.
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  AliCodeTimerAuto("",0);

  TrackerRecData(AliRecoParam::AConvert(Master()->GetEventSpecie()),kTRUE);
  
  BookHistograms(AliQAv1::kRECPOINTS);
  
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
  
  MappingCheckRecPoints(AliRecoParam::AConvert(Master()->GetEventSpecie()),kTRUE);
  //
  //ClonePerTrigClass(AliQAv1::kRECPOINTS); // this should be the last line  DONE at parent level
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
  //
  //ClonePerTrigClass(AliQAv1::kESDS); // this should be the last line  DONE at parent level
}

//____________________________________________________________________________
void AliMUONTrackerQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
	/// make QA for rawdata tracker
  	
  AliCodeTimerAuto(Form("%s",AliRecoParam::GetEventSpecieName(AliRecoParam::AConvert(Master()->GetEventSpecie()))),0);

  /// forces init
  GetRawsData(AliMUONQAIndices::kTrackerBusPatchOccupancy);
  
  AliMUONTrackerDataMaker* dm = static_cast<AliMUONTrackerDataMaker*>(TrackerDataMaker(AliRecoParam::AConvert(Master()->GetEventSpecie())));

  dm->SetRawReader(rawReader);
	
  Int_t eventType = rawReader->GetType();
  
  if (eventType == AliRawEventHeaderBase::kPhysicsEvent ) 
  {
    dm->ProcessEvent();
    
    FillEventSize(rawReader);
        
    if ( dm->LastEventWasEmpty() )
    {
      FillRawsData(AliMUONQAIndices::kTrackerReadoutStatus,1.0*AliMUONQAIndices::kTrackerRawNofEmptyEvents);
    }    
  }
  //
}

//__________________________________________________________________
void AliMUONTrackerQADataMakerRec::MakeDigits(TTree* digitsTree)         
{
  /// makes data from Digits
  AliCodeTimerAuto(Form("%s",AliRecoParam::GetEventSpecieName(AliRecoParam::AConvert(Master()->GetEventSpecie()))),0);
  
  /// forces init

  GetDigitsData(AliMUONQAIndices::kTrackerBusPatchOccupancy);

  if (!fDigitStore)
    fDigitStore = AliMUONVDigitStore::Create(*digitsTree);
  fDigitStore->Connect(*digitsTree, false);
  digitsTree->GetEvent(0);
  
  TIter next(fDigitStore->CreateIterator());
  
  AliMUONVDigit* dig = 0x0;
  
  AliMUON2DMap oneEventData(true);
  
  while ( ( dig = static_cast<AliMUONVDigit*>(next()) ) )
  {
    if ( dig->IsTracker() )
    {
      if ( dig->Charge() > 0.0 )
      {
        
        Int_t detElemId = dig->DetElemId();
        Int_t manuId = dig->ManuId();
        
        AliMUONVCalibParam* param = static_cast<AliMUONVCalibParam*>(oneEventData.FindObject(detElemId,manuId));
        if (!param)
        {
          param = new AliMUONCalibParamND(1,AliMpConstants::ManuNofChannels(),detElemId,manuId,
                                          AliMUONVCalibParam::InvalidFloatValue());
          oneEventData.Add(param);
        }
        param->SetValueAsDouble(dig->ManuChannel(),0,dig->Charge());
      }
    }
  }
  TrackerCalData(AliRecoParam::AConvert(Master()->GetEventSpecie()))->Add(oneEventData);
  //
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
	
  AliCodeTimerAuto(Form("%s",AliRecoParam::GetEventSpecieName(AliRecoParam::AConvert(Master()->GetEventSpecie()))),0);
  // Forces init by requesting an histogram
  GetRecPointsData(AliMUONQAIndices::kTrackerBusPatchOccupancy); 

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
	
  AliMUONQAMappingCheck* mcr = MappingCheckRecPoints(AliRecoParam::AConvert(Master()->GetEventSpecie()));
  
  if ( mcr ) mcr->NewEvent();
  
  AliMUON2DMap oneEventData(true);
  
	while ( ( cluster = static_cast<AliMUONVCluster*>(next()) ) )
	{
		Int_t detElemId = cluster->GetDetElemId();
		Int_t chamberId = AliMpDEManager::GetChamberId(detElemId);
		
		FillRecPointsData(AliMUONQAIndices::kTrackerNumberOfClustersPerDE,detElemId);
		FillRecPointsData(AliMUONQAIndices::kTrackerClusterChargePerDE+detElemId,cluster->GetCharge());
		FillRecPointsData(AliMUONQAIndices::kTrackerClusterMultiplicityPerDE+detElemId,cluster->GetNDigits());

		FillRecPointsData(AliMUONQAIndices::kTrackerNumberOfClustersPerChamber,chamberId);
		FillRecPointsData(AliMUONQAIndices::kTrackerClusterChargePerChamber+chamberId,cluster->GetCharge());
		FillRecPointsData(AliMUONQAIndices::kTrackerClusterMultiplicityPerChamber+chamberId,cluster->GetNDigits());
		FillRecPointsData(AliMUONQAIndices::kTrackerClusterHitMapPerChamber+chamberId,cluster->GetX(),cluster->GetY());
		
    if ( mcr ) mcr->Store(*cluster);
    
    for ( int i = 0; i < cluster->GetNDigits(); ++i ) 
    {
      UInt_t digitId = cluster->GetDigitId(i);
      
      Int_t manuId = AliMUONVDigit::ManuId(digitId);
      Int_t manuChannel = AliMUONVDigit::ManuChannel(digitId);
      
      AliMUONVCalibParam* param = static_cast<AliMUONVCalibParam*>(oneEventData.FindObject(detElemId,manuId));
      if (!param)
      {
        param = new AliMUONCalibParamND(1,AliMpConstants::ManuNofChannels(),detElemId,manuId,AliMUONVCalibParam::InvalidFloatValue());
        oneEventData.Add(param);
      }
      param->SetValueAsDouble(manuChannel,0,1.0);
    }
	}
	
  TrackerRecData(AliRecoParam::AConvert(Master()->GetEventSpecie()))->Add(oneEventData);    

  fClusterStore->Clear();
  //
}

//____________________________________________________________________________
void AliMUONTrackerQADataMakerRec::MakeESDs(AliESDEvent* esd)
{
  /// make QA data from ESDs

  AliCodeTimerAuto(Form("%s",AliRecoParam::GetEventSpecieName(AliRecoParam::AConvert(Master()->GetEventSpecie()))),0);
   
  // load ESD event in the interface
  AliMUONESDInterface esdInterface;
  if (GetRecoParam()) AliMUONESDInterface::ResetTracker(GetRecoParam(), kFALSE);
  else AliError("Unable to get recoParam: use default ones for residual calculation");
  esdInterface.LoadEvent(*esd);
  
  FillESDsData(AliMUONQAIndices::kESDnTracks,esdInterface.GetNTracks());
  
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
    
    FillESDsData(AliMUONQAIndices::kESDMomentum,esdTrack->P());
    FillESDsData(AliMUONQAIndices::kESDPt,esdTrack->Pt());
    FillESDsData(AliMUONQAIndices::kESDRapidity,esdTrack->Y());
    FillESDsData(AliMUONQAIndices::kESDChi2,track->GetNormalizedChi2());
    FillESDsData(AliMUONQAIndices::kESDProbChi2,TMath::Prob(track->GetGlobalChi2(),track->GetNDF()));
    FillESDsData(AliMUONQAIndices::kESDThetaX,esdTrack->GetThetaXUncorrected() / TMath::Pi() * 180.);
    FillESDsData(AliMUONQAIndices::kESDThetaY,esdTrack->GetThetaYUncorrected() / TMath::Pi() * 180.);
    FillESDsData(AliMUONQAIndices::kESDnClustersPerTrack,track->GetNClusters());
    
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
      
      FillESDsData(AliMUONQAIndices::kESDClusterHitMap+chId,cluster->GetX(), cluster->GetY());
      
      FillESDsData(AliMUONQAIndices::kESDnTotClustersPerCh,chId);
      FillESDsData(AliMUONQAIndices::kESDnTotClustersPerDE,deID);
      
      FillESDsData(AliMUONQAIndices::kESDClusterChargeInCh+chId,cluster->GetCharge());
      FillESDsData(AliMUONQAIndices::kESDSumClusterChargePerDE,deID, cluster->GetCharge());
      
      if (cluster->GetNDigits() > 0) { // discard clusters with pad not stored in ESD
	FillESDsData(AliMUONQAIndices::kESDnTotFullClustersPerDE,deID);
        FillESDsData(AliMUONQAIndices::kESDClusterSizeInCh+chId,cluster->GetNDigits());
	FillESDsData(AliMUONQAIndices::kESDSumClusterSizePerDE,deID, cluster->GetNDigits());
      }
      
      FillESDsData(AliMUONQAIndices::kESDResidualXInCh+chId,residualX);
      FillESDsData(AliMUONQAIndices::kESDResidualYInCh+chId,residualY);
      FillESDsData(AliMUONQAIndices::kESDSumResidualXPerDE,deID, residualX);
      FillESDsData(AliMUONQAIndices::kESDSumResidualYPerDE,deID, residualY);
      FillESDsData(AliMUONQAIndices::kESDSumResidualX2PerDE,deID, residualX*residualX);
      FillESDsData(AliMUONQAIndices::kESDSumResidualY2PerDE,deID, residualY*residualY);
      
      FillESDsData(AliMUONQAIndices::kESDLocalChi2XInCh+chId,localChi2X);
      FillESDsData(AliMUONQAIndices::kESDLocalChi2YInCh+chId,localChi2Y);
      FillESDsData(AliMUONQAIndices::kESDLocalChi2InCh+chId,localChi2);
      FillESDsData(AliMUONQAIndices::kESDSumLocalChi2XPerDE,deID, localChi2X);
      FillESDsData(AliMUONQAIndices::kESDSumLocalChi2YPerDE,deID, localChi2Y);
      FillESDsData(AliMUONQAIndices::kESDSumLocalChi2PerDE,deID, localChi2);
      
      trackParam = static_cast<AliMUONTrackParam*>(track->GetTrackParamAtCluster()->After(trackParam));
      
    }
    
  }

  FillESDsData(AliMUONQAIndices::kESDMatchTrig,nTrackMatchTrig);
  //
}

//____________________________________________________________________________ 
AliMUONVTrackerData* AliMUONTrackerQADataMakerRec::GetTrackerData() const
{ 
  /// Return tracker data
  
  return TrackerDataMaker(AliRecoParam::AConvert(Master()->GetEventSpecie()))->Data(); 
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
      
      if ( !hn.Contains("TrackerBusPatchConfig") )
      {
        AliDebug(1,Form("Resetting %s",hn.Data()));

        h->Reset();
      }
      else
      {
        AliDebug(1,Form("Will not reset histogram %s",hn.Data()));          
      }
      continue;
    }
    // RS account for the case when this histos were cloned per trigger class
    TObjArray* harr = dynamic_cast<TObjArray*>(o);
    if (harr) {
      TString hn(harr->GetName());
      if ( !hn.Contains("TrackerBusPatchConfig") ) {

	AliDebug(1,Form("Resetting clones of %s",hn.Data()));

	TIter nextCl(harr);
	TH1* hc = 0;
	while ( (hc=dynamic_cast<TH1*>(nextCl())) ) hc->Reset();
      }
      else {
	AliDebug(1,Form("Will not reset clones of histogram %s",hn.Data()));  
      }
      continue;
    }
    //
    AliMUONVTrackerData* d = dynamic_cast<AliMUONVTrackerData*>(o);
    if (d) {
      AliDebug(1,Form("Resetting %s",d->GetName()));
      d->Clear();
    }
    else {
      AliError(Form("Found an object of class %s. Do not know how to reset.",o->ClassName()));
    }
  }

  fLogger->Clear();
}

//____________________________________________________________________________ 
TObjArray* AliMUONTrackerQADataMakerRec::GetArray(TObjArray*& array, Bool_t create)
{
  /// Get (or create) the array

  if ( ! array ) 
  {
    if ( create ) 
    {
      array = new TObjArray(AliRecoParam::kNSpecies);
    }
  }
  
  return array;
}

//____________________________________________________________________________ 
AliMUONVTrackerDataMaker* 
AliMUONTrackerQADataMakerRec::TrackerDataMaker(Int_t specieIndex)  const
{
  /// const version of the getter
  if ( fTrackerDataMakerArray )
  {
    return static_cast<AliMUONVTrackerDataMaker*>(fTrackerDataMakerArray->At(specieIndex));
  }
  return 0x0;
}

//____________________________________________________________________________ 
AliMUONVTrackerDataMaker* 
AliMUONTrackerQADataMakerRec::TrackerDataMaker(Int_t specieIndex, Bool_t create)
{
  /// Get (or create) TrackerDataMaker object for a given specie
  
  TObjArray* array = GetArray(fTrackerDataMakerArray,create);
  TObject* o(0x0);
  
  if ( array ) 
  {
    array->SetOwner(kTRUE);
    o = array->At(specieIndex);
    if (!o && create)
    {
      
      AliMUONTrackerDataMaker* dm = new AliMUONTrackerDataMaker(0x0,
                                                                AliCDBManager::Instance()->GetRun(),
                                                                0x0,
                                                                "",
                                                                "",
                                                                kFALSE,
                                                                0.0,0.0);
      
      if (!fLogger) fLogger = new AliMUONLogger(-1); // note that we share the logger between species... should not be a big deal though
      dm->EnableErrorLogger(fLogger);
      dm->Data()->DisableChannelLevel(); // to save up disk space, we only store starting at the manu level    
      dm->Data()->SetName("RawCharges");
      dm->SetRunning(kTRUE);

      o = dm;
      array->AddAt(o,specieIndex);      
    }
  }
  return static_cast<AliMUONVTrackerDataMaker*>(o);
}

//____________________________________________________________________________ 
AliMUONVTrackerData* 
AliMUONTrackerQADataMakerRec::TrackerCalData(Int_t specieIndex, Bool_t create)
{
  TObjArray* array = GetArray(fTrackerCalDataArray,create);
  TObject* o(0x0);
  
  if (array)
  {
    o = array->At(specieIndex);
    if (!o && create)
    {
      AliMUONTrackerData* data = new AliMUONTrackerData("CalCharges",Form("Calibrated charges (fC) %s",GetRecoParam()->GetCalibrationMode()),1);
      data->SetDimensionName(0,"charge");
      data->DisableChannelLevel(); // to save up disk space, we only store starting at the manu level
      o=data;
      array->AddAt(o,specieIndex);
    }
  }
  return static_cast<AliMUONVTrackerData*>(o);
}

//____________________________________________________________________________ 
AliMUONVTrackerData* 
AliMUONTrackerQADataMakerRec::TrackerRecData(Int_t specieIndex, Bool_t create)
{
  TObjArray* array = GetArray(fTrackerRecDataArray,create);
  TObject* o(0x0);
  
  if (array)
  {
    o = array->At(specieIndex);
    if (!o && create)
    {
      AliMUONTrackerData* data = new AliMUONTrackerData("RecCharges",Form("Calibrated charges (fC) %s for digits belonging to a reconstructed cluster",GetRecoParam()->GetCalibrationMode()),1);
      data->SetDimensionName(0,"one");
      data->DisableChannelLevel(); // to save up disk space, we only store starting at the manu level
      o=data;
      array->AddAt(o,specieIndex);
    }
  }
  return static_cast<AliMUONVTrackerData*>(o);
}

//____________________________________________________________________________ 
AliMUONQAMappingCheck* 
AliMUONTrackerQADataMakerRec::MappingCheckRecPoints(Int_t specieIndex, Bool_t create)
{
  TObjArray* array = GetArray(fMappingCheckRecPointsArray,create);
  TObject* o(0x0);
  
  if (array)
  {
    o = array->At(specieIndex);
    array->SetOwner(kTRUE);
    if (!o && create)
    {
      AliMUONQAMappingCheck* mcheck = new AliMUONQAMappingCheck(RunNumber()); 
      o=mcheck;
      array->AddAt(o,specieIndex);
    }
  }
  return static_cast<AliMUONQAMappingCheck*>(o);
}
