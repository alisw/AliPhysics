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

// $Id: AliMUONTriggerQADataMakerRec.cxx 35760 2009-10-21 21:45:42Z ivana $

// --- MUON header files ---
#include "AliMUONTriggerQADataMakerRec.h"

//-----------------------------------------------------------------------------
/// \class AliMUONTriggerQADataMakerRec
///
/// MUON class for quality assurance data (histo) maker
///
/// \author C. Finck, D. Stocco, L. Aphecetche, A. Blanc

/// \cond CLASSIMP
ClassImp(AliMUONTriggerQADataMakerRec)
/// \endcond
           
#include "AliCodeTimer.h"
#include "AliMUONConstants.h"
#include "AliMpConstants.h"
#include "AliMUONTriggerDisplay.h"
#include "TH2.h"
#include "TH1F.h"
#include "TString.h"
#include "AliRecoParam.h"
#include "AliMUONDigitStoreV2R.h"
#include "AliMUONTriggerStoreV1.h"
#include "AliMpCDB.h"
#include "AliMUONRawStreamTriggerHP.h"
#include "AliMpDDLStore.h"
#include "AliMpTriggerCrate.h"
#include "AliMpLocalBoard.h"
#include "AliQAv1.h"
#include "AliRawReader.h"
#include "AliMUONDigitMaker.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONRecoParam.h"
#include "AliMUONTriggerElectronics.h"
#include "AliMUONCalibrationData.h"
#include "AliDCSValue.h"
#include "AliMpDCSNamer.h"
#include "AliMpDEManager.h"
#include "AliMpDEIterator.h"
#include "AliCDBManager.h"
#include "TTree.h"
#include "AliMUONGlobalTriggerBoard.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONGlobalCrateConfig.h"
#include "AliMUONQAIndices.h"
#include "AliMpPad.h"
#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"

namespace
{
  Double_t ProtectedSqrt(Double_t x)
  {
    return ( x > 0.0 ? TMath::Sqrt(x) : 0.0 );
  }
}
//____________________________________________________________________________ 
AliMUONTriggerQADataMakerRec::AliMUONTriggerQADataMakerRec(AliQADataMakerRec* master) : 
AliMUONVQADataMakerRec(master),
fDigitMaker(new AliMUONDigitMaker(kFALSE)),
fCalibrationData(0x0),
fTriggerProcessor(0x0),
fDigitStore(0x0)
{
    /// ctor
}


//__________________________________________________________________
AliMUONTriggerQADataMakerRec::~AliMUONTriggerQADataMakerRec()
{
    /// dtor
  delete fDigitMaker;
  delete fDigitStore;
  delete fTriggerProcessor;
  delete fCalibrationData;
}

//____________________________________________________________________________ 
void AliMUONTriggerQADataMakerRec::EndOfDetectorCycleESDs(Int_t /*specie*/, TObjArray** /*list*/)
{
  /// Normalize ESD histograms
}
  
//____________________________________________________________________________ 
void AliMUONTriggerQADataMakerRec::EndOfDetectorCycleRecPoints(Int_t /*specie*/, TObjArray** /*list*/)
{
  /// Normalize RecPoints histograms
}


//____________________________________________________________________________ 
void AliMUONTriggerQADataMakerRec::EndOfDetectorCycleRaws(Int_t /*specie*/, TObjArray** /*list*/)
{
  /// create Raws histograms in Raws subdir
  
  Int_t histoRawsIndex[] = {
    AliMUONQAIndices::kTriggerErrorSummary,
    AliMUONQAIndices::kTriggerCalibSummary,
    AliMUONQAIndices::kTriggerReadOutErrors,
    AliMUONQAIndices::kTriggerGlobalOutput
  };
  Int_t histoRawsScaledIndex[] = {
    AliMUONQAIndices::kTriggerErrorSummaryNorm,
    AliMUONQAIndices::kTriggerCalibSummaryNorm,
    AliMUONQAIndices::kTriggerReadOutErrorsNorm,
    AliMUONQAIndices::kTriggerGlobalOutputNorm
  };
  
  const Int_t kNrawsHistos = sizeof(histoRawsIndex)/sizeof(histoRawsIndex[0]);
  Float_t scaleFactor[kNrawsHistos] = {100., 100., 100., 1.};

  for ( Int_t itc=-1; itc<AliQADataMakerRec::GetNTrigClasses(); itc++) { 
  
    DisplayTriggerInfo(itc);

    // Normalize RawData histos
    TH1* histo1D = GetRawsData(AliMUONQAIndices::kTriggerRawNAnalyzedEvents,itc);
    // This histogram is there for all relevant triggers
    // if it is not there, it means that the trigger is not taken into account
    // so we can skip the trigger class for all other histos
    if ( ! histo1D ) continue;
    Float_t nbevent = histo1D->GetBinContent(1);
    for(Int_t ihisto=0; ihisto<kNrawsHistos; ihisto++){
      TH1* inputHisto = GetRawsData(histoRawsIndex[ihisto],itc);
      TH1* scaledHisto = GetRawsData(histoRawsScaledIndex[ihisto],itc);
      // Check here for both since we do not clone Calib-only histograms
      if ( scaledHisto && inputHisto &&  nbevent > 0 ) {
        scaledHisto->Reset();
        scaledHisto->Add(inputHisto);
        scaledHisto->Scale(scaleFactor[ihisto]/nbevent);
      }
    } // loop on histos

    
    // The following histograms are surely there
    // if the histogram with analyzed events is there:
    // test on the existence of each histogram is not necessary
    TH1* hYCopy = GetRawsData(AliMUONQAIndices::kTriggerErrorLocalYCopy,itc); //number of YCopy error per board
    TH1* hYCopyTests = GetRawsData(AliMUONQAIndices::kTriggerErrorLocalYCopyTest,itc); //contains the number of YCopy test per board
    TH1* hYCopyNorm = GetRawsData(AliMUONQAIndices::kTriggerErrorLocalYCopyNorm,itc); 
    hYCopyNorm->Reset();
    hYCopyNorm->Divide(hYCopy, hYCopyTests, 100., 1.);
     
    Float_t mean = hYCopyNorm->Integral();
      
    TH1* hSummary = GetRawsData(AliMUONQAIndices::kTriggerErrorSummaryNorm,itc);
    hSummary->SetBinContent(AliMUONQAIndices::kAlgoLocalYCopy+1,mean/192.); //put the mean of the % of YCopy error in the kTriggerError's corresponding bin

    TH1F* hTriggerRatio = (TH1F*)GetRawsData(AliMUONQAIndices::kTriggerLocalRatio4434,itc);
    if ( hTriggerRatio ){
      hTriggerRatio->Divide(((TH1F*)GetRawsData(AliMUONQAIndices::kTriggerNumberOf44Dec,itc)),((TH1F*)GetRawsData(AliMUONQAIndices::kTriggerNumberOf34Dec,itc)));

      FillRatio4434Histos(1,itc,kTRUE);

      //reset bins temporary used to store informations
      ((TH1F*)GetRawsData(AliMUONQAIndices::kTriggerRatio4434AllEvents,itc))->SetBinContent(0,0); 
      Int_t nbins =  ((TH1F*)GetRawsData(AliMUONQAIndices::kTriggerRatio4434AllEvents,itc))->GetNbinsX();
      ((TH1F*)GetRawsData(AliMUONQAIndices::kTriggerRatio4434AllEvents,itc))->SetBinContent(nbins+1,0);

      ((TH1F*)GetRawsData(AliMUONQAIndices::kTriggerLocalRatio4434,itc))->SetMaximum(1.1);
      ((TH1F*)GetRawsData(AliMUONQAIndices::kTriggerRatio4434AllEvents,itc))->SetMaximum(1.1);
      ((TH1F*)GetRawsData(AliMUONQAIndices::kTriggerRatio4434SinceLastUpdate,itc))->SetMaximum(1.1);
    }
  
    if ( GetRawsData(AliMUONQAIndices::kTriggerGlobalScalersNorm,itc) ) {
      TH1* inputHisto = GetRawsData(AliMUONQAIndices::kTriggerGlobalScalers,itc);
      TH1* scaledHisto = GetRawsData(AliMUONQAIndices::kTriggerGlobalScalersNorm,itc);
      scaledHisto->Reset();
      scaledHisto->Add(inputHisto);
      Float_t scaleValue = ((TH1F*)GetRawsData(AliMUONQAIndices::kTriggerScalersTime,itc))->GetBinContent(1);
      if ( scaleValue > 0. ) scaledHisto->Scale(1./scaleValue);
    }
  } // loop on trigger classes
}

//____________________________________________________________________________ 
void AliMUONTriggerQADataMakerRec::InitRaws()
{
    /// create Raws histograms in Raws subdir
  
  // RS: Since there is no sense in cloning trigger scalers per trigger, I am (for the moment) forbidding their cloning

  AliCodeTimerAuto("",0);
  
  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 
 
  TString boardName = "Local board Id";

  Int_t nbLocalBoard = AliMUONConstants::NTriggerCircuit();

  TH1F* histo1D = 0x0;
  TH2F* histo2D = 0x0;

  AliMUONTriggerDisplay triggerDisplay;

  TString histoName, histoTitle;
  if ( GetRecoParam()->GetEventSpecie() == AliRecoParam::kCalib ) {
    histo1D = new TH1F("hTriggerScalersTime", "Acquisition time from trigger scalers", 1, 0.5, 1.5);
    histo1D->GetXaxis()->SetBinLabel(1, "One-bin histogram: bin is filled at each scaler event.");
    histo1D->GetYaxis()->SetTitle("Cumulated scaler time (s)");
    Add2RawsList(histo1D, AliMUONQAIndices::kTriggerScalersTime, expert, !image, !saveCorr);
    ForbidCloning(histo1D); // RS

    for(Int_t iCath=0; iCath<AliMpConstants::NofCathodes(); iCath++){
      TString cathName = ( iCath==0 ) ? "BendPlane" : "NonBendPlane";
      for(Int_t iChamber=0; iChamber<AliMpConstants::NofTriggerChambers(); iChamber++){
	histoName = Form("hTriggerScalers%sChamber%i", cathName.Data(), 11+iChamber);
	histoTitle = Form("Chamber %i - %s: trigger scaler counts", 11+iChamber, cathName.Data());
	histo2D = new TH2F(histoName.Data(), histoTitle.Data(),
			   nbLocalBoard, 0.5, (Float_t)nbLocalBoard + 0.5,
			   16, -0.5, 15.5);
	histo2D->GetXaxis()->SetTitle(boardName.Data());
	histo2D->GetYaxis()->SetTitle("Strip");	
	histo2D->SetOption("COLZ");	
	Add2RawsList(histo2D, AliMUONQAIndices::kTriggerScalers + AliMpConstants::NofTriggerChambers()*iCath + iChamber, expert, !image, !saveCorr);
	ForbidCloning(histo2D); // RS
      } // loop on chambers
    } // loop on cathodes
	
    for(Int_t iCath=0; iCath<AliMpConstants::NofCathodes(); iCath++){
      TString cathName = ( iCath==0 ) ? "BendPlane" : "NonBendPlane";
      for(Int_t iChamber=0; iChamber<AliMpConstants::NofTriggerChambers(); iChamber++){
	histoName = Form("hTriggerScalersDisplay%sChamber%i", cathName.Data(), 11+iChamber);
	histoTitle = Form("Chamber %i - %s: Hit rate from scalers (Hz/cm^{2})", 11+iChamber, cathName.Data());
	histo2D = (TH2F*)triggerDisplay.GetEmptyDisplayHisto(histoName, AliMUONTriggerDisplay::kDisplayStrips, 
							     iCath, iChamber, histoTitle);
	histo2D->SetOption("COLZ");
	Add2RawsList(histo2D, AliMUONQAIndices::kTriggerScalersDisplay + AliMpConstants::NofTriggerChambers()*iCath + iChamber, expert, !image, !saveCorr);
	ForbidCloning(histo2D); // RS
      } // loop on chambers
    } // loop on cathodes    

    TString axisLabel[AliMUONQAIndices::kNtrigCalibSummaryBins] = {"#splitline{Dead}{Channels}", "#splitline{Dead}{Local Boards}", "#splitline{Dead}{Regional Boards}", "#splitline{Dead}{Global Board}", "#splitline{Noisy}{Strips}"};

    TH1F* histoCalib = new TH1F("hTriggerCalibSummaryAll", "MTR calibration summary counts", AliMUONQAIndices::kNtrigCalibSummaryBins, -0.5, (Float_t)AliMUONQAIndices::kNtrigCalibSummaryBins - 0.5);
    for (Int_t ibin=1; ibin<=AliMUONQAIndices::kNtrigCalibSummaryBins; ibin++){
      histoCalib->GetXaxis()->SetBinLabel(ibin, axisLabel[ibin-1].Data());
    }
    histoCalib->SetFillColor(kBlue);
    histoCalib->GetYaxis()->SetTitle("Counts");
    // Copy of previous histo for scaling purposes
    TH1F* histoCalibNorm = (TH1F*)histoCalib->Clone("hTriggerCalibSummary");
    histoCalibNorm->SetTitle("MTR calibration summary");
    histoCalibNorm->SetOption("bartext0");
    histoCalibNorm->GetYaxis()->SetTitle("Percentage per event (%)");
    // Adding both histos after cloning to avoid problems with the expert bit
    Add2RawsList(histoCalib,     AliMUONQAIndices::kTriggerCalibSummary,      expert, !image, !saveCorr);
    ForbidCloning(histoCalib); // RS

    Add2RawsList(histoCalibNorm, AliMUONQAIndices::kTriggerCalibSummaryNorm, !expert,  image, !saveCorr);
    ForbidCloning(histoCalibNorm); // RS

  } // Calibration reco param
	
  const char *globalXaxisName[6] = {"US HPt", "US LPt", "LS HPt", "LS LPt", "SGL HPt", "SGL LPt"};
  const char *allLevelXaxisName[AliMUONQAIndices::kNtrigAlgoErrorBins] = {"Local algo X", "Local algo Y", "Local LUT","Local Y Copy" , "Local2Regional", "Regional", "Regional2Global", "GlobalFromInGlobal", "GlobalFromInLocal", "GlobalFromOutLocal"};
  const char *readoutErrNames[AliMUONQAIndices::kNtrigStructErrorBins]={"Local","Regional","Global","DARC"};

  TString errorAxisTitle = "Number of errors";

  histo1D = new TH1F("hTriggerErrorLocalXPos", "ErrorLocalXPos",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, AliMUONQAIndices::kTriggerErrorLocalXPos, expert, !image, !saveCorr);

  histo1D = new TH1F("hTriggerErrorLocalYPos", "ErrorLocalYPos",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, AliMUONQAIndices::kTriggerErrorLocalYPos, expert, !image, !saveCorr);

  histo1D = new TH1F("hTriggerErrorLocalDev", "ErrorLocalDev",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, AliMUONQAIndices::kTriggerErrorLocalDev, expert, !image, !saveCorr);

  histo1D = new TH1F("hTriggerErrorLocalTriggerDec", "ErrorLocalTriggerDec",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, AliMUONQAIndices::kTriggerErrorLocalTriggerDec, expert, !image, !saveCorr);

  histo1D = new TH1F("hTriggerErrorLocalLPtLSB", "ErrorLocalLPtLSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, AliMUONQAIndices::kTriggerErrorLocalLPtLSB, expert, !image, !saveCorr);

  histo1D = new TH1F("hTriggerErrorLocalLPtMSB", "ErrorLocalLPtMSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, AliMUONQAIndices::kTriggerErrorLocalLPtMSB, expert, !image, !saveCorr);

  histo1D = new TH1F("hTriggerErrorLocalHPtLSB", "ErrorLocalHPtLSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, AliMUONQAIndices::kTriggerErrorLocalHPtLSB, expert, !image, !saveCorr);

  histo1D = new TH1F("hTriggerErrorLocalHPtMSB", "ErrorLocalHPtMSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, AliMUONQAIndices::kTriggerErrorLocalHPtMSB, expert, !image, !saveCorr);

  histo1D = new TH1F("hTriggerErrorLocalTrigY", "ErrorLocalTrigY",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, AliMUONQAIndices::kTriggerErrorLocalTrigY, expert, !image, !saveCorr);

  if ( GetRecoParam()->GetEventSpecie() != AliRecoParam::kCalib ) {
    histo1D = new TH1F("hTriggerRatio4434Local", "Ratio4434Local",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
    histo1D->GetXaxis()->SetTitle(boardName.Data());
    histo1D->GetYaxis()->SetTitle("ratio 44/34");
    Add2RawsList(histo1D, AliMUONQAIndices::kTriggerLocalRatio4434, expert, !image, !saveCorr);                                               
    histo1D = new TH1F("hTriggerRatio4434AllEvents", "Ratio4434AllEvents",1,0,1);
    histo1D->GetXaxis()->SetTitle("Event number");
    histo1D->GetYaxis()->SetTitle("ratio 44/34");
    histo1D->SetLineColor(4);                           
    Add2RawsList(histo1D, AliMUONQAIndices::kTriggerRatio4434AllEvents, expert, !image, !saveCorr);                                               
    histo1D = new TH1F("hTriggerRatio4434SinceLastUpdate", "Ratio4434SinceLastUpdate",1,0,1);
    histo1D->GetXaxis()->SetTitle("Event number");
    histo1D->GetYaxis()->SetTitle("ratio 44/34");                           
    Add2RawsList(histo1D, AliMUONQAIndices::kTriggerRatio4434SinceLastUpdate, expert, !image, !saveCorr);
  }

  histo1D = new TH1F("hTriggerErrorLocal2RegionalLPtLSB", "ErrorLocal2RegionalLPtLSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, AliMUONQAIndices::kTriggerErrorLocal2RegionalLPtLSB, expert, !image, !saveCorr);

  histo1D = new TH1F("hTriggerErrorLocal2RegionalLPtMSB", "ErrorLocal2RegionalLPtMSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, AliMUONQAIndices::kTriggerErrorLocal2RegionalLPtMSB, expert, !image, !saveCorr);

  histo1D = new TH1F("hTriggerErrorLocal2RegionalHPtLSB", "ErrorLocal2RegionalHPtLSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, AliMUONQAIndices::kTriggerErrorLocal2RegionalHPtLSB, expert, !image, !saveCorr);

  histo1D = new TH1F("hTriggerErrorLocal2RegionalHPtMSB", "ErrorLocal2RegionalHPtMSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, AliMUONQAIndices::kTriggerErrorLocal2RegionalHPtMSB, expert, !image, !saveCorr);

  histo1D = new TH1F("hTriggerErrorOutGlobalFromInGlobal", "ErrorOutGlobalFromInGlobal",6,-0.5,6-0.5);
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  for (int ibin=0;ibin<6;ibin++){
    histo1D->GetXaxis()->SetBinLabel(ibin+1,globalXaxisName[ibin]);
  }
  Add2RawsList(histo1D, AliMUONQAIndices::kTriggerErrorOutGlobalFromInGlobal, expert, !image, !saveCorr);

  histo1D = new TH1F("hTriggerErrorOutGlobalFromInLocal", "ErrorOutGlobalFromInLocal",6,-0.5,6-0.5);
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  for (int ibin=0;ibin<6;ibin++){
    histo1D->GetXaxis()->SetBinLabel(ibin+1,globalXaxisName[ibin]);
  }
  Add2RawsList(histo1D, AliMUONQAIndices::kTriggerErrorOutGlobalFromInLocal, expert, !image, !saveCorr);

  TH1F* histoAlgoErr = new TH1F("hTriggerAlgoNumOfErrors", "Trigger Algorithm total errors",AliMUONQAIndices::kNtrigAlgoErrorBins,-0.5,(Float_t)AliMUONQAIndices::kNtrigAlgoErrorBins-0.5);
  histoAlgoErr->GetYaxis()->SetTitle("Number of events with errors");
  for (int ibin=0;ibin<AliMUONQAIndices::kNtrigAlgoErrorBins;ibin++){
    histoAlgoErr->GetXaxis()->SetBinLabel(ibin+1,allLevelXaxisName[ibin]);
  }
  histoAlgoErr->SetFillColor(kBlue);
  // Copy of previous histo for scaling purposes
  TH1F* histoAlgoErrNorm = (TH1F*)histoAlgoErr->Clone("hTriggerAlgoErrors");
  histoAlgoErrNorm->SetOption("bartext0");
  histoAlgoErrNorm->SetTitle("Trigger algorithm errors");
  histoAlgoErrNorm->GetYaxis()->SetTitle("% of events with errors");
  // Adding both histos after cloning to avoid problems with the expert bit
  Add2RawsList(histoAlgoErr,     AliMUONQAIndices::kTriggerErrorSummary,      expert, !image, !saveCorr);
  Add2RawsList(histoAlgoErrNorm, AliMUONQAIndices::kTriggerErrorSummaryNorm, !expert,  image, !saveCorr);  

  histo1D = new TH1F("hTriggerTriggeredBoards", "Triggered boards", nbLocalBoard, 0.5, (Float_t)nbLocalBoard + 0.5);
  Add2RawsList(histo1D, AliMUONQAIndices::kTriggeredBoards, expert, !image, !saveCorr);

  histo2D = (TH2F*)triggerDisplay.GetEmptyDisplayHisto("hTriggerFiredBoardsDisplay", AliMUONTriggerDisplay::kDisplayBoards,
						       0, 0, "Local board triggers / event");
  histo2D->SetOption("COLZ");
  Add2RawsList(histo2D, AliMUONQAIndices::kTriggerBoardsDisplay, expert, !image, !saveCorr);

  TH1F* histoYCopyErr = new TH1F("hTriggerErrorLocalYCopy", "Number of YCopy errors",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histoYCopyErr->GetXaxis()->SetTitle(boardName.Data());
  histoYCopyErr->GetYaxis()->SetTitle(errorAxisTitle.Data());
  // Copy of previous histo for scaling purposes
  TH1F* histoYCopyErrTest = (TH1F*)histoYCopyErr->Clone("hTriggerErrorLocalYCopyTest");
  histoYCopyErrTest->SetTitle("Number of YCopy tested");
  // Copy of previous histo for scaling purposes
  TH1F* histoYCopyErrNorm = (TH1F*)histoYCopyErr->Clone("hTriggerErrorLocalYCopyNorm");
  histoYCopyErrNorm->SetTitle("% of YCopy errors");
  // Adding both histos after cloning to avoid problems with the expert bit
  Add2RawsList(histoYCopyErr,     AliMUONQAIndices::kTriggerErrorLocalYCopy,     expert, !image, !saveCorr);
  Add2RawsList(histoYCopyErrTest, AliMUONQAIndices::kTriggerErrorLocalYCopyTest, expert, !image, !saveCorr);
  Add2RawsList(histoYCopyErrNorm, AliMUONQAIndices::kTriggerErrorLocalYCopyNorm, expert, !image, !saveCorr);

  TH1F* histoROerr = new TH1F("hTriggerReadoutNumOfErrors","Trigger Read-Out total errors", AliMUONQAIndices::kNtrigStructErrorBins, -0.5, (Float_t)AliMUONQAIndices::kNtrigStructErrorBins-0.5);
  histoROerr->GetYaxis()->SetTitle("Fraction of errors");
  histoROerr->SetFillColor(kBlue);
  for (int ibin=0;ibin<AliMUONQAIndices::kNtrigStructErrorBins;ibin++){
    histoROerr->GetXaxis()->SetBinLabel(ibin+1,readoutErrNames[ibin]);
  }
  // Copy of previous histo for scaling purposes
  TH1F* histoROerrNorm = (TH1F*)histoROerr->Clone("hTriggerReadoutErrors");
  histoROerrNorm->SetTitle("Trigger Read-Out errors");
  histoROerrNorm->SetOption("bartext0");
  histoROerrNorm->GetYaxis()->SetTitle("% of errors per event");
  // Adding both histos after cloning to avoid problems with the expert bit
  Add2RawsList(histoROerr,     AliMUONQAIndices::kTriggerReadOutErrors,      expert, !image, !saveCorr);
  Add2RawsList(histoROerrNorm, AliMUONQAIndices::kTriggerReadOutErrorsNorm, !expert,  image, !saveCorr);

  TH1F* histoGlobalMult = new TH1F("hTriggerGlobalOutMultiplicity","Trigger global outputs multiplicity", 6, -0.5, 6.-0.5);
  histoGlobalMult->GetYaxis()->SetTitle("Number of triggers"); 
  histoGlobalMult->GetXaxis()->SetTitle("Global output");
  for (int ibin=0;ibin<6;ibin++){
    histoGlobalMult->GetXaxis()->SetBinLabel(ibin+1,globalXaxisName[ibin]);
  }        
  histoGlobalMult->SetFillColor(kBlue);
  // Copy of previous histo for scaling purposes
  TH1F* histoGlobalMultNorm = (TH1F*)histoGlobalMult->Clone("hTriggerGlobalOutMultiplicityPerEvt");
  histoGlobalMultNorm->SetTitle("Trigger global outputs multiplicity per event");
  histoGlobalMultNorm->SetOption("bartext0");
  //histoGlobalMultNorm->SetBarWidth(0.5);
  //histoGlobalMultNorm->SetBarOffset(0.25);
  histoGlobalMultNorm->GetYaxis()->SetTitle("Triggers per event");
  // Adding both histos after cloning to avoid problems with the expert bit
  Add2RawsList(histoGlobalMult,     AliMUONQAIndices::kTriggerGlobalOutput,     expert, !image, !saveCorr);
  Add2RawsList(histoGlobalMultNorm, AliMUONQAIndices::kTriggerGlobalOutputNorm, expert, !image, !saveCorr);

  histo1D = new TH1F("hTriggerRawNAnalyzedEvents", "Number of analyzed events per specie", 1, 0.5, 1.5);
  Int_t esindex = AliRecoParam::AConvert(CurrentEventSpecie());
  histo1D->GetXaxis()->SetBinLabel(1, AliRecoParam::GetEventSpecieName(esindex));
  histo1D->GetYaxis()->SetTitle("Number of analyzed events");
  Add2RawsList(histo1D, AliMUONQAIndices::kTriggerRawNAnalyzedEvents, expert, !image, !saveCorr);

  if ( GetRecoParam()->GetEventSpecie() != AliRecoParam::kCalib ) {
    histo1D = new TH1F("hTriggerNumberOf34Dec", "Number of 3/4",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
    histo1D->GetXaxis()->SetTitle(boardName.Data());
    histo1D->GetYaxis()->SetTitle("Number of 3/4");
    Add2RawsList(histo1D, AliMUONQAIndices::kTriggerNumberOf34Dec, expert, !image, !saveCorr);

    histo1D = new TH1F("hTriggerNumberOf44Dec", "Number of 4/4",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
    histo1D->GetXaxis()->SetTitle(boardName.Data());
    histo1D->GetYaxis()->SetTitle("Number of 4/4");
    Add2RawsList(histo1D, AliMUONQAIndices::kTriggerNumberOf44Dec, expert, !image, !saveCorr);
  }
  
  histo1D = new TH1F("hTriggerIsThere","trigger is there",1,0,1);
  Add2RawsList(histo1D,AliMUONQAIndices::kTriggerIsThere,kTRUE,kFALSE,kFALSE);

  if ( GetRecoParam()->GetEventSpecie() == AliRecoParam::kCalib ) {
    TH1F* histoGlobalScalers = new TH1F("hTriggerGlobalScalers","Trigger global scalers", 6, -0.5, 6.-0.5);
    histoGlobalScalers->GetYaxis()->SetTitle("L0 counts");
    histoGlobalScalers->GetXaxis()->SetTitle("Global output");
    for (int ibin=0;ibin<6;ibin++){
      histoGlobalScalers->GetXaxis()->SetBinLabel(ibin+1,globalXaxisName[ibin]);
    }        
    // Copy of previous histo for scaling purposes
    TH1F* histoGlobalScalersNorm = (TH1F*)histoGlobalScalers->Clone("hTriggerGlobalScalersRate");
    histoGlobalScalersNorm->SetTitle("Trigger global L0 scalers rate");
    histoGlobalScalersNorm->SetOption("etext0");
    histoGlobalScalersNorm->GetYaxis()->SetTitle("L0 scalers rate (Hz)");
    // Adding both histos after cloning to avoid problems with the expert bit
    Add2RawsList(histoGlobalScalers,     AliMUONQAIndices::kTriggerGlobalScalers,     expert, !image, !saveCorr);
    ForbidCloning(histoGlobalScalers); // RS
    Add2RawsList(histoGlobalScalersNorm, AliMUONQAIndices::kTriggerGlobalScalersNorm, expert, !image, !saveCorr);
    ForbidCloning(histoGlobalScalersNorm); // RS
  }
  //
  //ClonePerTrigClass(AliQAv1::kRAWS); // RS: this should be the last line  DONE at parent level
  //
}

//__________________________________________________________________
void AliMUONTriggerQADataMakerRec::InitDigits() 
{
  /// Initialized Digits spectra 
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1I* h0 = new TH1I("hTriggerDigitsDetElem", "Detection element distribution in Digits;Detection element Id;Counts",  400, 1100, 1500); 
  Add2DigitsList(h0, 0, !expert, image);
  ForbidCloning(h0);
  //
  //ClonePerTrigClass(AliQAv1::kDIGITS); // this should be the last line  DONE at parent level
  //
} 

//____________________________________________________________________________ 
void AliMUONTriggerQADataMakerRec::InitRecPoints()
{
	/// create Reconstructed Points histograms in RecPoints subdir for the
	/// MUON Trigger subsystem.

  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  TH1F* histo1D = 0x0;

  histo1D = new TH1F("hTriggerNAnalyzedEvents", "Number of analyzed events per specie", 1, 0.5, 1.5);
  Int_t esindex = AliRecoParam::AConvert(CurrentEventSpecie());
  histo1D->GetXaxis()->SetBinLabel(1, AliRecoParam::GetEventSpecieName(esindex));
  histo1D->GetYaxis()->SetTitle("Number of analyzed events");
  Add2RecPointsList(histo1D, AliMUONQAIndices::kTriggerNAnalyzedEvents, expert, !image);
  ForbidCloning(histo1D);

  histo1D = new TH1F("hTriggerTrippedChambers", "Trigger RPCs in trip", 418, 1100-0.5, 1417+0.5);
  histo1D->GetXaxis()->SetTitle("DetElemId");
  histo1D->GetYaxis()->SetTitle("# of trips");
  histo1D->SetFillColor(kRed);
  histo1D->SetLineColor(kRed);
  Add2RecPointsList(histo1D, AliMUONQAIndices::kTriggerRPCtrips, !expert, image);
  ForbidCloning(histo1D);   // RS this histo is not cloned
  //
  FillTriggerDCSHistos();
  //
  //ClonePerTrigClass(AliQAv1::kRECPOINTS); DONE at parent level
  //
}


//____________________________________________________________________________ 
void AliMUONTriggerQADataMakerRec::InitESDs()
{
  /// Empty implementation
}

//____________________________________________________________________________
void AliMUONTriggerQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
	/// make QA for rawdata trigger

    AliCodeTimerAuto("",0);
	
    // Init Local/Regional/Global decision with fake values
    //

    UInt_t globalInput[4];
    for (Int_t bit=0; bit<4; bit++){
	globalInput[bit]=0;
    }

    //for (Int_t reg=0;reg<16;reg++){
    //fTriggerOutputRegionalData[reg]=0;
    //for (Int_t bit=0;bit<4;bit++){
    //fTriggerInputGlobalDataLPt[reg][bit]=0;
    //fTriggerInputGlobalDataHPt[reg][bit]=0;
    //}
    //}

    AliMUONDigitStoreV2R digitStore;
    
    AliMUONTriggerStoreV1 recoTriggerStore;

    AliMUONTriggerStoreV1 inputTriggerStore;

    AliMUONGlobalTrigger inputGlobalTrigger;

    UShort_t maxNcounts = 0xFFFF;
    
    // Get trigger Local, Regional, Global in/outputs and scalers

    Int_t loCircuit=0;
    AliMpCDB::LoadDDLStore();

    const AliMUONRawStreamTriggerHP::AliHeader*          darcHeader  = 0x0;
    const AliMUONRawStreamTriggerHP::AliRegionalHeader*  regHeader   = 0x0;
    const AliMUONRawStreamTriggerHP::AliLocalStruct*     localStruct = 0x0;

    Int_t nDeadLocal = 0, nDeadRegional = 0, nDeadGlobal = 0, nNoisyStrips = 0;
    Int_t nFiredStrips = 0, nStripsTot = 0;

    // When a crate is not present, the loop on boards is not performed
    // This should allow to correctly count the local boards
    Int_t countNotifiedBoards = 0, countAllBoards = 0;

    Bool_t containTriggerData = kFALSE;
    AliMUONRawStreamTriggerHP rawStreamTrig(rawReader);
    while (rawStreamTrig.NextDDL()) 
      {
       containTriggerData = kTRUE;

      Bool_t scalerEvent =  rawReader->GetDataHeader()->GetL1TriggerMessage() & 0x1;

      Bool_t fillScalerHistos = ( scalerEvent && 
				  ( GetRecoParam()->GetEventSpecie() == AliRecoParam::kCalib ) );

      if ( scalerEvent != fillScalerHistos ) {
	Int_t esindex = AliRecoParam::AConvert(CurrentEventSpecie());
	AliWarning(Form("Scaler event found but event specie is %s. Scaler histos will not be filled", AliRecoParam::GetEventSpecieName(esindex)));
      }

      darcHeader = rawStreamTrig.GetHeaders();

      if (darcHeader->GetGlobalFlag()){
        if ( fillScalerHistos ) {
          UInt_t nOfClocks = darcHeader->GetGlobalClock();
          Double_t nOfSeconds = ((Double_t) nOfClocks) / 40e6; // 1 clock each 25 ns
          FillRawsData(AliMUONQAIndices::kTriggerScalersTime, 1., nOfSeconds);
          const UInt_t* globScaler = darcHeader->GetGlobalScaler();
          Int_t bitCorr[6] = {2,0,3,1,4,5};
          for (Int_t bit=0; bit<6; bit++){
            FillRawsData(AliMUONQAIndices::kTriggerGlobalScalers, bitCorr[bit],(double)(*(globScaler+bit)));
          }
        }

        //Get Global datas
        inputGlobalTrigger.SetFromGlobalResponse(darcHeader->GetGlobalOutput());
        Int_t resp[6] = {inputGlobalTrigger.PairUnlikeHpt(), inputGlobalTrigger.PairUnlikeLpt(),
          inputGlobalTrigger.PairLikeHpt(), inputGlobalTrigger.PairLikeLpt(),
          inputGlobalTrigger.SingleHpt(), inputGlobalTrigger.SingleLpt()}; 
        for (Int_t bit=0; bit<6; bit++){
          if ( resp[bit] == 0 ){
            if ( fillScalerHistos )
              nDeadGlobal++;
          }
          else
            FillRawsData(AliMUONQAIndices::kTriggerGlobalOutput, bit, resp[bit]);
        } // loop on bits

        //for (Int_t Bit=0; Bit<32; Bit++){
        //fTriggerInputGlobalDataLPt[Bit/4][Bit%4]=((darcHeader->GetGlobalInput(0)>>Bit)&1);
        //fTriggerInputGlobalDataLPt[Bit/4+8][Bit%4]=((darcHeader->GetGlobalInput(1)>>Bit)&1);
        //fTriggerInputGlobalDataHPt[Bit/4][Bit%4]=((darcHeader->GetGlobalInput(2)>>Bit)&1);
        //fTriggerInputGlobalDataHPt[Bit/4+8][Bit%4]=((darcHeader->GetGlobalInput(3)>>Bit)&1);
        //}

        for (Int_t i=0; i<4; i++){
          globalInput[i]=darcHeader->GetGlobalInput(i);
        }
      }

      Int_t nReg = rawStreamTrig.GetRegionalHeaderCount();

      for(Int_t iReg = 0; iReg < nReg ;iReg++)
      {   //reg loop

	//Int_t regId=rawStreamTrig.GetDDL()*8+iReg;

	// crate info  
	  AliMpTriggerCrate* crate = AliMpDDLStore::Instance()->GetTriggerCrate(rawStreamTrig.GetDDL(), iReg);

	  regHeader =  rawStreamTrig.GetRegionalHeader(iReg);

	  //Get regional outputs -> not checked, hardware read-out doesn't work
	  //fTriggerOutputRegionalData[regId]=Int_t(regHeader->GetOutput());
	  // if ( ! fTriggerOutputRegionalData[regId] )
	  // nDeadRegional++;
	Int_t nBoardsInReg = 0; // Not necessary when regional output will work

	// loop over local structures
	Int_t nLocal = regHeader->GetLocalStructCount();

	for(Int_t iLocal = 0; iLocal < nLocal; iLocal++) 
	{
	    
	    localStruct = regHeader->GetLocalStruct(iLocal);

	  // if card exist
	  if (!localStruct) continue;
          
	  loCircuit = crate->GetLocalBoardId(localStruct->GetId());

	  if ( !loCircuit ) continue; // empty slot

	  AliMpLocalBoard* localBoard = AliMpDDLStore::Instance()->GetLocalBoard(loCircuit, false);

	  nBoardsInReg++; // Not necessary when regional output will work
	  countAllBoards++;

	  // skip copy cards
	  if( !localBoard->IsNotified()) 
	    continue;

	  AliMUONLocalTrigger inputLocalTrigger;
	  inputLocalTrigger.SetLocalStruct(loCircuit, *localStruct);
	  inputTriggerStore.Add(inputLocalTrigger);

	  countNotifiedBoards++;  

	  TArrayS xyPattern[2];	  
	  localStruct->GetXPattern(xyPattern[0]);
	  localStruct->GetYPattern(xyPattern[1]);
	  fDigitMaker->TriggerDigits(loCircuit, xyPattern, digitStore);

	  //Get electronic Decisions from data

	  //Get regional inputs -> not checked, hardware read-out doesn't work
	  //fTriggerInputRegionalDataLPt[0][loCircuit]=Int_t(((regHeader->GetInput(0))>>(2*iLocal))&1);
	  //fTriggerInputRegionalDataLPt[1][loCircuit]=Int_t(((regHeader->GetInput(1))>>((2*iLocal)+1))&1);

	  //Get local in/outputs
	  if (Int_t(localStruct->GetDec())!=0){
	    FillRawsData(AliMUONQAIndices::kTriggeredBoards,loCircuit);
	  }
	  else if ( fillScalerHistos ){
	    nDeadLocal++;
	  }

	  // loop over strips
	  if ( fillScalerHistos ) {
	    Int_t cathode = localStruct->GetComptXY()%2;
      
      Int_t offset = 0;
      if (cathode && localBoard->GetSwitch(AliMpLocalBoard::kZeroAllYLSB)) offset = -8;

	    for (Int_t ibitxy = 0; ibitxy < 16; ++ibitxy) {
	      if (ibitxy==0){
		AliDebug(AliQAv1::GetQADebugLevel(),"Filling trigger scalers");
	      }

	      UShort_t scalerVal[4] = {
		localStruct->GetXY1(ibitxy),
		localStruct->GetXY2(ibitxy),
		localStruct->GetXY3(ibitxy),
		localStruct->GetXY4(ibitxy)
	      };
        
        

        for(Int_t ich=0; ich<AliMpConstants::NofTriggerChambers(); ich++){
          // getDetElemId
          Int_t detElemId = AliMpDDLStore::Instance()->GetDEfromLocalBoard(loCircuit, ich);
					
          const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, AliMp::GetCathodType(cathode));
					
					
          Int_t istrip = ibitxy + offset;
					
          AliMpPad pad = seg->PadByLocation(loCircuit,istrip,kFALSE);
          if (!pad.IsValid()) continue;
          nStripsTot++;
          
          // UShort_t pattern = (UShort_t)xyPattern[cathode].At(ich); 
          // if ((pattern >> ibitxy) & 0x1) nFiredStrips++;
          
          if ( scalerVal[ich] > 0 ) {
            FillRawsData(AliMUONQAIndices::kTriggerScalers + AliMpConstants::NofTriggerChambers()*cathode + ich,
			 loCircuit, istrip, 2*(Float_t)scalerVal[ich]);
            nFiredStrips++;
          }

          if ( scalerVal[ich] >= maxNcounts )
            nNoisyStrips++;
        } // loop on chamber
	    } // loop on strips
	  } // scaler event
	} // iLocal
	if ( nBoardsInReg == 0 )
	  nDeadRegional++; // Not necessary when regional output will work
      } // iReg

      Float_t readoutErrors[AliMUONQAIndices::kNtrigStructErrorBins] = {
	((Float_t)rawStreamTrig.GetLocalEoWErrors())/((Float_t)countAllBoards),
	((Float_t)rawStreamTrig.GetRegEoWErrors())/16.,
	((Float_t)rawStreamTrig.GetGlobalEoWErrors())/6.,
	((Float_t)rawStreamTrig.GetDarcEoWErrors())/2.
      };
    
      for (Int_t ibin=0; ibin<AliMUONQAIndices::kNtrigStructErrorBins; ibin++){
	if ( readoutErrors[ibin] > 0 )
	  FillRawsData(AliMUONQAIndices::kTriggerReadOutErrors, ibin, readoutErrors[ibin]);
      }
    } // NextDDL

    if ( ! containTriggerData ) return;

    FillRawsData(AliMUONQAIndices::kTriggerRawNAnalyzedEvents,1.);

    nDeadLocal += AliMUONConstants::NTriggerCircuit() - countNotifiedBoards;
    if ( nStripsTot > 0 ) { // The value is != 0 only for scaler events
      AliDebug(AliQAv1::GetQADebugLevel(), Form("nStripsFired %i  nStripsTot %i", nFiredStrips, nStripsTot));
      Float_t fraction[AliMUONQAIndices::kNtrigCalibSummaryBins] = {
	((Float_t)(nStripsTot - nFiredStrips)) / ((Float_t)nStripsTot),
	//(Float_t)nDeadLocal / ((Float_t)countNotifiedBoards),
	(Float_t)nDeadLocal / ((Float_t)AliMUONConstants::NTriggerCircuit()),
	(Float_t)nDeadRegional / 16.,
	(Float_t)nDeadGlobal / 6., // Number of bits of global response
	(Float_t)nNoisyStrips / ((Float_t)nStripsTot),
      };

      for(Int_t ibin = 0; ibin < AliMUONQAIndices::kNtrigCalibSummaryBins; ibin++){
	if ( fraction[ibin] > 0. )
	  FillRawsData(AliMUONQAIndices::kTriggerCalibSummary,ibin, fraction[ibin]);
      }
    }

  TriggerElectronics()->Digits2Trigger(digitStore,recoTriggerStore);

  AliMUONGlobalTrigger* recoGlobalTriggerFromLocal;
  recoGlobalTriggerFromLocal = recoTriggerStore.Global();

  //Reconstruct Global decision from Global inputs
  UChar_t recoResp = RawTriggerInGlobal2OutGlobal(globalInput);
  AliMUONGlobalTrigger recoGlobalTriggerFromGlobal;
  recoGlobalTriggerFromGlobal.SetFromGlobalResponse(recoResp);

  // Compare data and reconstructed decisions and fill histos
  RawTriggerMatchOutLocal(inputTriggerStore, recoTriggerStore);
  //Fill ratio 44/34 histos
  for ( Int_t itc=-1; itc<AliQADataMakerRec::GetNEventTrigClasses(); ++itc ) FillRatio4434Histos(fgkUpdateRatio4434, itc, kFALSE);
  //RawTriggerMatchOutLocalInRegional(); // Not tested, hardware read-out doesn't work
  RawTriggerMatchOutGlobal(inputGlobalTrigger, recoGlobalTriggerFromGlobal, 'G');
  // Global, reconstruction from Local inputs: compare data and reconstructed decisions and fill histos
  RawTriggerMatchOutGlobal(inputGlobalTrigger, *recoGlobalTriggerFromLocal, 'L');
  // Global, reconstruction from Global inputs: compare data and reconstructed decisions and fill histos
  //
}

//__________________________________________________________________
void AliMUONTriggerQADataMakerRec::MakeDigits(TTree* digitsTree)         
{
  /// makes data from Digits

  // Do nothing in case of calibration event
  if ( GetRecoParam()->GetEventSpecie() == AliRecoParam::kCalib ) return;
  
  if (!fDigitStore)
    fDigitStore = AliMUONVDigitStore::Create(*digitsTree);
  
  fDigitStore->Clear();
  fDigitStore->Connect(*digitsTree, false);
  digitsTree->GetEvent(0);
  
  TIter next(fDigitStore->CreateIterator());
  
  AliMUONVDigit* dig = 0x0;
  
  while ( ( dig = static_cast<AliMUONVDigit*>(next()) ) )
    {
      FillDigitsData(0,dig->DetElemId());
    }
}

//____________________________________________________________________________
void AliMUONTriggerQADataMakerRec::MakeRecPoints(TTree* /*clustersTree*/)
{
  /// Fill histogram with total number of analyzed events for normalization purposes

  // Do nothing in case of calibration event
  if ( GetRecoParam()->GetEventSpecie() == AliRecoParam::kCalib ) return;
	
  FillRecPointsData(AliMUONQAIndices::kTriggerNAnalyzedEvents,1.);
}

//____________________________________________________________________________
void AliMUONTriggerQADataMakerRec::MakeESDs(AliESDEvent* /*esd*/)
{  
  /// Empty implementation
}


//____________________________________________________________________________ 
void AliMUONTriggerQADataMakerRec::DisplayTriggerInfo(Int_t itc)
{
  //
  /// Display trigger information in a user-friendly way:
  /// from local board and strip numbers to their position on chambers
  //

  AliMUONTriggerDisplay triggerDisplay;
  
  TH2* histoStrips=0x0;
  TH2* histoDisplayStrips=0x0;
  if ( GetRawsData(AliMUONQAIndices::kTriggerScalers, itc) ) {
    AliMUONTriggerDisplay::EDisplayOption displayOption = AliMUONTriggerDisplay::kNormalizeToArea;
    for (Int_t iCath = 0; iCath < AliMpConstants::NofCathodes(); iCath++)
      {    
	for (Int_t iChamber = 0; iChamber < AliMpConstants::NofTriggerChambers(); iChamber++)
	  {
	    histoStrips = (TH2*)GetRawsData(AliMUONQAIndices::kTriggerScalers + AliMpConstants::NofTriggerChambers()*iCath + iChamber, itc);

	    if(histoStrips->GetEntries()==0) continue; // No events found => No need to display

	    histoDisplayStrips = (TH2*)GetRawsData(AliMUONQAIndices::kTriggerScalersDisplay + AliMpConstants::NofTriggerChambers()*iCath + iChamber, itc);

	    triggerDisplay.FillDisplayHistogram(histoStrips, histoDisplayStrips,
						AliMUONTriggerDisplay::kDisplayStrips, iCath, iChamber, displayOption);

	    Float_t scaleValue = ((TH1*)GetRawsData(AliMUONQAIndices::kTriggerScalersTime, itc))->GetBinContent(1);
	    if(scaleValue>0.) histoDisplayStrips->Scale(1./scaleValue);
	  } // iChamber
      } // iCath
  }
  
  if ( GetRawsData(AliMUONQAIndices::kTriggeredBoards, itc) ){
    TH1* histoBoards = (TH1*)GetRawsData(AliMUONQAIndices::kTriggeredBoards, itc);
    TH2* histoDisplayBoards = (TH2*)GetRawsData(AliMUONQAIndices::kTriggerBoardsDisplay, itc);
    triggerDisplay.FillDisplayHistogram(histoBoards, histoDisplayBoards, AliMUONTriggerDisplay::kDisplayBoards, 0, 0);
    Float_t scaleValue = GetRawsData(AliMUONQAIndices::kTriggerRawNAnalyzedEvents, itc)->GetBinContent(1);
    if(scaleValue>0.) histoDisplayBoards->Scale(1./scaleValue);
  }
}


//_____________________________________________________________________________
Bool_t 
AliMUONTriggerQADataMakerRec::FillTriggerDCSHistos()
{
  /// Get HV and currents values for one trigger chamber
  // RS: Note: the histos involved in this routin are forbidden to be cloned, -1 in GetRawsData returns the default histos
  int itc = -1;
  //
  AliCodeTimerAuto("",0);
  
  TMap* triggerDcsMap = CalibrationData()->TriggerDCS();

  if ( !triggerDcsMap ) 
  {
    AliError("Cannot fill DCS histos, as triggerDcsMap is NULL");
    return kFALSE;
  }

  const Double_t kMaxDelay = 3.;
  const Double_t kMaxVariation = 25.; // Volts
  const Int_t kDefaultNpoints = 200;
  Double_t scaleFactor = 1./1000.;

  Bool_t error = kFALSE;
  Bool_t expert   = kTRUE;
  Bool_t image    = kTRUE;

  AliMpDEIterator deIt;
  
  AliMpDCSNamer triggerDcsNamer("TRIGGER");

  TH2F* currHisto = 0x0;
  Int_t histoIndex = 0;
  TString histoName, histoTitle;

  TArrayD axisSlat(18+1);
  for(Int_t islat=0; islat<=18; islat++){
    axisSlat[islat] = -0.5 + (Float_t)islat;
  }

  TArrayD axisTimeAux[4], axisTime[4], axisTimeDE(kDefaultNpoints);
  TArrayI index[4], npoints(4);

  // Build axis of times
  npoints.Reset();
  for(Int_t ich=0; ich<4; ich++){
    axisTimeAux[ich].Set(kDefaultNpoints);
    axisTimeAux[ich].Reset(-1.);
  }

  deIt.First();
  while ( !deIt.IsDone() )
  {
    Int_t detElemId = deIt.CurrentDEId();
    TObjArray* values = GetDCSValues(AliMpDCSNamer::kDCSHV, detElemId, triggerDcsMap, triggerDcsNamer);

    if ( values ) {

      AliDebug(AliQAv1::GetQADebugLevel(), Form("DetElemId %i", detElemId));

      axisTimeDE.Reset(-1.);

      TIter next(values);
      AliDCSValue* val = 0x0;
      Double_t previousVal = -999.;
      Int_t npointsde = 0;
      while ( ( val = static_cast<AliDCSValue*>(next()) ) )
      {
	if ( npointsde + 1 > kDefaultNpoints ) {
	  axisTimeDE.Set(npointsde + 1);
	}

	Double_t currVal = val->GetFloat();
	Double_t currTime = (Double_t)val->GetTimeStamp();
	if (npointsde > 0 ){
	  if ( TMath::Abs( currVal - previousVal ) < kMaxVariation && 
	       TMath::Abs( currTime - axisTimeDE[npointsde-1] ) < 40 ) continue;
	}

	axisTimeDE[npointsde] = currTime;
	previousVal = currVal;
	npointsde++;
      } // loop on values

      //      AliDebug(AliQAv1::GetQADebugLevel(), Form("Adding DE point %2i  (%2i)  %.2f  (%i)\n", previousBin, npointsde, axisTimeDE[previousBin], nTimesPerBin));

      Int_t iChamber = AliMpDEManager::GetChamberId(detElemId);
      Int_t ich = iChamber - AliMpConstants::NofTrackingChambers();

      for(Int_t ipde=0; ipde<npointsde; ipde++){

	if ( npoints[ich] + 1 > kDefaultNpoints ) {
	  axisTimeAux[ich].Set(npoints[ich] + 1);
	}

	for(Int_t ipoint = 0; ipoint < axisTimeAux[ich].GetSize(); ipoint++){
	  if (axisTimeAux[ich][ipoint] < 0.) {
	    axisTimeAux[ich][ipoint] = axisTimeDE[ipde];
	    npoints[ich]++;
	    AliDebug(AliQAv1::GetQADebugLevel(), Form("Adding point %2i  %.0f\n", ipoint, axisTimeAux[ich][ipoint]));
	    break;
	  }
	  if ( TMath::Abs( axisTimeDE[ipde] - axisTimeAux[ich][ipoint]) < kMaxDelay ) {
	    axisTimeAux[ich][ipoint] = TMath::Min(axisTimeAux[ich][ipoint], axisTimeDE[ipde]);
	    break;
	  }
	} // loop on points
      } // loop on reorganized values

    } // if ( values ) 
    deIt.Next();
  } // loop on DetElemId

  for(Int_t ich=0; ich<4; ich++){
    axisTimeAux[ich].Set(npoints[ich]);
    index[ich].Set(npoints[ich]);
    TMath::Sort(npoints[ich], axisTimeAux[ich].GetArray(), index[ich].GetArray(), kFALSE);

    axisTime[ich].Set(npoints[ich]+1);
    for(Int_t ipoint = 0; ipoint < axisTimeAux[ich].GetSize(); ipoint++){
      axisTime[ich][ipoint] = axisTimeAux[ich][index[ich][ipoint]];
    }
    Double_t minStartEndWidth = 0.1 * (axisTime[ich][npoints[ich]-1] - axisTime[ich][0]);
    axisTime[ich][npoints[ich]] = axisTime[ich][npoints[ich]-1] + minStartEndWidth;
    if ( npoints[ich] >= 1)
      if ( axisTime[ich][1] - axisTime[ich][0] < minStartEndWidth )
	axisTime[ich][0] = axisTime[ich][1] - minStartEndWidth;
  }


  // Loop again on detection elements: create and fill histos
  deIt.First();
  while ( !deIt.IsDone() )
  {
    Int_t detElemId = deIt.CurrentDEId();
    TObjArray* values = GetDCSValues(AliMpDCSNamer::kDCSHV, detElemId, triggerDcsMap, triggerDcsNamer);

    if ( values ) {
      Int_t iChamber = AliMpDEManager::GetChamberId(detElemId);
      Int_t ich = iChamber - AliMpConstants::NofTrackingChambers();

      histoIndex = AliMUONQAIndices::kTriggerRPChv + ich;
      histoName = Form("hTriggerRPCHVChamber%i", 11+ich);
      histoTitle = Form("Chamber %i: RPC HV (kV)", 11+ich);

      currHisto = (TH2F*)GetRecPointsData(histoIndex,itc); // RS this histo is not cloned

      if(!currHisto){
	currHisto  = new TH2F(histoName.Data(), histoTitle.Data(),
			      npoints[ich], axisTime[ich].GetArray(),
			      18, axisSlat.GetArray());
	currHisto->GetXaxis()->SetTitle("Time");
	currHisto->GetXaxis()->SetTimeDisplay(1);
	//currHisto->GetXaxis()->SetTimeFormat("%d%b%y %H:%M:%S");
	currHisto->GetXaxis()->SetLabelSize(0.03);
	currHisto->GetYaxis()->SetTitle("RPC");
	currHisto->SetOption("TEXT45COLZ");
	Add2RecPointsList(currHisto, histoIndex, expert, !image);
	ForbidCloning(currHisto); // RS
      }

      Int_t slat = detElemId%100;
      Int_t slatBin = currHisto->GetYaxis()->FindBin(slat);

      TIter next(values);
      AliDCSValue* val = 0x0;
      Double_t sumValuesPerBin = 0.;
      Int_t nValuesPerBin = 0;
      Int_t previousBin = -1;
      Double_t previousTime = -1., previousVal = -999., sumVal = 0., sumTime = 0.;
      Bool_t isTrip = kFALSE;
      Int_t nPointsForSlope = 0;
      while ( ( val = static_cast<AliDCSValue*>(next()) ) )
      {
	Double_t currTime = (Double_t)val->GetTimeStamp();
	Int_t currentBin = currHisto->GetXaxis()->FindBin(currTime+0.5);
	Double_t currVal = val->GetFloat();
	Double_t deltaVal = currVal - previousVal;
	Bool_t isRepeated = kFALSE;
	if ( previousTime > 0 ){
	  isRepeated = ( TMath::Abs( currVal - previousVal ) < kMaxVariation && 
			 TMath::Abs( currTime - previousTime ) < 40 );

	  // Check for trips
	  sumTime += currTime - previousTime;
	  sumVal += deltaVal;
	  nPointsForSlope++;

	  if ( sumTime > 0. && nPointsForSlope >= 3 ){
	    Double_t slope = sumVal / sumTime;
	    if ( slope < -10. ) // going down of more than 10V/s
	      isTrip = kTRUE;
	  }

	  if ( deltaVal * sumVal < 0. ) {
	    sumTime = 0.;
	    sumVal = 0.;
	    nPointsForSlope = 0;
	  }
	}

	if ( ! isRepeated ) {
	  if ( currentBin != previousBin ) {
	    if ( previousBin >= 0 ) {
	      currHisto->SetBinContent(previousBin, slatBin, scaleFactor*sumValuesPerBin/((Double_t)nValuesPerBin));
	      sumValuesPerBin = 0.;
	      nValuesPerBin = 0;
	    }
	    previousBin = currentBin;
	  }
	}
	  
	sumValuesPerBin += currVal;
	nValuesPerBin++;
	previousTime = currTime;
	previousVal = currVal;
      } // loop on values
      currHisto->SetBinContent(previousBin, slatBin, scaleFactor*sumValuesPerBin/((Double_t)nValuesPerBin)); // Fill last value
      if ( isTrip ) ((TH1*)GetRecPointsData(AliMUONQAIndices::kTriggerRPCtrips,itc))->Fill(detElemId);
    } // if ( values ) 
    deIt.Next();
  } // loop on detElem
  return error;
}


//____________________________________________________________________________ 
TObjArray* 
AliMUONTriggerQADataMakerRec::GetDCSValues(Int_t iMeas, Int_t detElemId,
					   TMap* triggerDcsMap, AliMpDCSNamer& triggerDcsNamer)
{
  //
  /// Get values of DCS data points from the map
  //

  if ( AliMpDEManager::GetStationType(detElemId) != AliMp::kStationTrigger) return 0x0;

  TString currAlias = triggerDcsNamer.DCSChannelName(detElemId, 0, iMeas);

  TPair* triggerDcsPair = static_cast<TPair*>(triggerDcsMap->FindObject(currAlias.Data()));

  if (!triggerDcsPair)
  {
    AliError(Form("Did not find expected alias (%s) for DE %d\n",
                  currAlias.Data(),detElemId));
    return 0x0;
  }

  TObjArray* values = static_cast<TObjArray*>(triggerDcsPair->Value());
  if (!values)
  {
    AliError(Form("Could not get values for alias %s\n",currAlias.Data()));
    return 0x0;
  }

  return values;
}


//____________________________________________________________________________ 
UChar_t AliMUONTriggerQADataMakerRec::RawTriggerInGlobal2OutGlobal(UInt_t globalInput[4])
{
  //
  /// Reconstruct Global Trigger decision using Global Inputs
  //

    AliCodeTimerAuto("",0);

    AliMUONGlobalCrateConfig* globalConfig = CalibrationData()->GlobalTriggerCrateConfig();

    AliMUONGlobalTriggerBoard globalTriggerBoard;
    globalTriggerBoard.Reset();
    for (Int_t i = 0; i < 4; i++) {
	globalTriggerBoard.Mask(i,globalConfig->GetGlobalMask(i));
    }

    globalTriggerBoard.RecomputeRegional(globalInput);
    globalTriggerBoard.Response();
    return globalTriggerBoard.GetResponse();

}

//____________________________________________________________________________ 
void AliMUONTriggerQADataMakerRec::RawTriggerMatchOutLocal(const AliMUONVTriggerStore& inputTriggerStore,
							   const AliMUONVTriggerStore& recoTriggerStore)
{
  //
  /// Match data and reconstructed Local Trigger decision

  AliCodeTimerAuto("",0);

  Bool_t skipBoard[234];
  memset(skipBoard,0,AliMUONConstants::NTriggerCircuit()*sizeof(Bool_t));

  Bool_t errorInYCopy = kFALSE;

  // First search for YCopy errors.
  Int_t loCircuit = -1;
  TIter next(recoTriggerStore.CreateLocalIterator());
  AliMUONLocalTrigger* recoLocalTrigger, *inputLocalTrigger;
  while ( ( recoLocalTrigger = static_cast<AliMUONLocalTrigger*>(next()) ) )
  {  
    loCircuit = recoLocalTrigger->LoCircuit();
    Int_t iboard = loCircuit - 1;

    FillRawsData(AliMUONQAIndices::kTriggerErrorLocalYCopyTest,loCircuit);
  
    inputLocalTrigger = inputTriggerStore.FindLocal(loCircuit);

    Int_t recoTrigPattern[4]  = {recoLocalTrigger->GetY1Pattern(), recoLocalTrigger->GetY2Pattern(), recoLocalTrigger->GetY3Pattern(), recoLocalTrigger->GetY4Pattern()};
    Int_t inputTrigPattern[4] = {inputLocalTrigger->GetY1Pattern(), inputLocalTrigger->GetY2Pattern(), inputLocalTrigger->GetY3Pattern(), inputLocalTrigger->GetY4Pattern()};

    AliMpLocalBoard* localBoardMp = AliMpDDLStore::Instance()->GetLocalBoard(loCircuit); // get local board object for switch value

    Bool_t errorInCopyBoard = kFALSE;
    for(Int_t ich=0; ich<4; ich++){
      if ( recoTrigPattern[ich] != inputTrigPattern[ich] ){
	skipBoard[iboard] = kTRUE;
	if ( ich >=2 ){
	  if ( localBoardMp->GetSwitch(AliMpLocalBoard::kOR0) )
	    skipBoard[iboard+1] = kTRUE;
	  if ( localBoardMp->GetSwitch(AliMpLocalBoard::kOR1) )
	    skipBoard[iboard-1] = kTRUE;
	}
	errorInCopyBoard = kTRUE;
	errorInYCopy = kTRUE;
      }
    } // loop on chambers
    if ( errorInCopyBoard )
      FillRawsData(AliMUONQAIndices::kTriggerErrorLocalYCopy,loCircuit);    
  } // loop on local boards

  if (errorInYCopy)
    FillRawsData(AliMUONQAIndices::kTriggerErrorSummary,AliMUONQAIndices::kAlgoLocalYCopy);
  
  Bool_t errorInXPosDev = kFALSE;
  Bool_t errorInYPosTrigY = kFALSE;
  Bool_t errorInLUT = kFALSE;

  next.Reset();
  Bool_t respBendPlane, respNonBendPlane;
  while ( ( recoLocalTrigger = static_cast<AliMUONLocalTrigger*>(next()) ) )
  {  
    loCircuit = recoLocalTrigger->LoCircuit();
    Int_t iboard = loCircuit - 1;
  
    // Fill ratio 44/34 histos (if not scaler event)
    if ( GetRecoParam()->GetEventSpecie() != AliRecoParam::kCalib ) {
      Bool_t is34 = ( recoLocalTrigger->GetLoDecision() != 0 );
      Bool_t is44 = TriggerElectronics()->ModifiedLocalResponse(loCircuit, respBendPlane, respNonBendPlane, kTRUE);
      if ( is34 ) FillRawsData(AliMUONQAIndices::kTriggerNumberOf34Dec,loCircuit);
      if ( is44 ) FillRawsData(AliMUONQAIndices::kTriggerNumberOf44Dec,loCircuit);

      if ( is44 && ! is34 )
	AliWarning("Event satisfies the 4/4 conditions but not the 3/4");
    }
    
    inputLocalTrigger = inputTriggerStore.FindLocal(loCircuit);

    if ( recoLocalTrigger->LoStripX() != inputLocalTrigger->LoStripX() ) {
      FillRawsData(AliMUONQAIndices::kTriggerErrorLocalXPos,loCircuit);
      errorInXPosDev = kTRUE;
    }
    
    if ( recoLocalTrigger->GetDeviation() != inputLocalTrigger->GetDeviation() ) {
      FillRawsData(AliMUONQAIndices::kTriggerErrorLocalDev,loCircuit);
      errorInXPosDev = kTRUE;
    }

    // Skip following checks in case we previously found YCopy error and YPos or trigY errors
    if ( (!skipBoard[iboard]) || ( (recoLocalTrigger->LoStripY() == inputLocalTrigger->LoStripY()) && (recoLocalTrigger->LoTrigY() == inputLocalTrigger->LoTrigY())) ) {
	
	if ( recoLocalTrigger->GetLoDecision() != inputLocalTrigger->GetLoDecision() ) {
	  FillRawsData(AliMUONQAIndices::kTriggerErrorLocalTriggerDec,loCircuit);
	}
	
	// Test Hpt and LPT
	Int_t recoLut[2]  = { recoLocalTrigger->LoLpt(),  recoLocalTrigger->LoHpt() };
	Int_t inputLut[2] = {inputLocalTrigger->LoLpt(), inputLocalTrigger->LoHpt() };
	Int_t currIndex[2][2] = {{AliMUONQAIndices::kTriggerErrorLocalLPtLSB, AliMUONQAIndices::kTriggerErrorLocalLPtMSB},
				 {AliMUONQAIndices::kTriggerErrorLocalHPtMSB, AliMUONQAIndices::kTriggerErrorLocalHPtMSB}};
	for (Int_t ilut=0; ilut<2; ilut++){
	    Int_t bitDiff = recoLut[ilut]^inputLut[ilut];
	    if ( bitDiff == 0 ) continue;
	    for (Int_t ibit=0; ibit<2; ibit++){
		Bool_t isBitDifferent = (bitDiff>>ibit)&1;
		if ( isBitDifferent ){
		  FillRawsData(currIndex[ilut][ibit],loCircuit);
		  errorInLUT = kTRUE;
		}
	    }
	}
    }
    
 
    // Skip following checks in case we previously found YCopy errors
    if ( skipBoard[iboard] ) continue;

    if ( recoLocalTrigger->LoStripY() != inputLocalTrigger->LoStripY() ) {
      FillRawsData(AliMUONQAIndices::kTriggerErrorLocalYPos,loCircuit);
      errorInYPosTrigY = kTRUE;
    }

    if ( recoLocalTrigger->LoTrigY() != inputLocalTrigger->LoTrigY()  ) {
      FillRawsData(AliMUONQAIndices::kTriggerErrorLocalTrigY,loCircuit);	
      errorInYPosTrigY = kTRUE;
    }
  } // loop on local boards
  
  if (errorInXPosDev)
    FillRawsData(AliMUONQAIndices::kTriggerErrorSummary,AliMUONQAIndices::kAlgoLocalX);

  if (errorInLUT)
    FillRawsData(AliMUONQAIndices::kTriggerErrorSummary,AliMUONQAIndices::kAlgoLocalLUT);

  if (errorInYPosTrigY)
    FillRawsData(AliMUONQAIndices::kTriggerErrorSummary,AliMUONQAIndices::kAlgoLocalY);

}
/*
//____________________________________________________________________________ 
void AliMUONTriggerQADataMakerRec::RawTriggerMatchOutLocalInRegional()
{
  //
  /// Match Local outputs and Regional inputs
  /// Not tested, hardware read-out doesn't work
  //

    for (int localId=1;localId<235;localId++){
	if(fTriggerOutputLocalDataLPtDec[0][localId]!=fTriggerInputRegionalDataLPt[0][localId]){
	    ((TH1F*)GetRawsData(kTriggerErrorLocal2RegionalLPtLSB))->Fill(localId);
	}
	if(fTriggerOutputLocalDataLPtDec[1][localId]!=fTriggerInputRegionalDataLPt[1][localId]){
	    ((TH1F*)GetRawsData(kTriggerErrorLocal2RegionalLPtMSB))->Fill(localId);
	}
	if(fTriggerOutputLocalDataHPtDec[0][localId]!=fTriggerInputRegionalDataHPt[0][localId]){
	    ((TH1F*)GetRawsData(kTriggerErrorLocal2RegionalHPtLSB))->Fill(localId);
	}
	if(fTriggerOutputLocalDataHPtDec[1][localId]!=fTriggerInputRegionalDataHPt[1][localId]){
	    ((TH1F*)GetRawsData(kTriggerErrorLocal2RegionalHPtMSB))->Fill(localId);
	}
    }
}
*/


//____________________________________________________________________________ 
void AliMUONTriggerQADataMakerRec::RawTriggerMatchOutGlobal(AliMUONGlobalTrigger& inputGlobalTrigger, 
									AliMUONGlobalTrigger& recoGlobalTrigger, 
									Char_t histo)
{
  //
  /// Match data and reconstructed Global Trigger decision for a reconstruction from Global inputs.
  /// histo='G': fill FromGlobalInput histo='L': fill from Local input;
  //

  if ( recoGlobalTrigger.GetGlobalResponse() == inputGlobalTrigger.GetGlobalResponse() )
    return;
  Int_t histoToFill;
  Int_t binToFill;
  
  if (histo=='G'){
      histoToFill=AliMUONQAIndices::kTriggerErrorOutGlobalFromInGlobal;
      binToFill=AliMUONQAIndices::kAlgoGlobalFromGlobal;
  }else{
      if (histo=='L'){
	  histoToFill=AliMUONQAIndices::kTriggerErrorOutGlobalFromInLocal;
	  binToFill=AliMUONQAIndices::kAlgoGlobalFromLocal;
      }else{
	  AliWarning(Form("Global histos not filled, 3rd argument must be 'G' or 'L'"));
	  return;
      } 
  }

  FillRawsData(AliMUONQAIndices::kTriggerErrorSummary,binToFill);
  
  Int_t inputResp[6] = {inputGlobalTrigger.PairUnlikeHpt(), inputGlobalTrigger.PairUnlikeLpt(),
			inputGlobalTrigger.PairLikeHpt(), inputGlobalTrigger.PairLikeLpt(),
			inputGlobalTrigger.SingleHpt(), inputGlobalTrigger.SingleLpt()};
  
  Int_t recoResp[6] = {recoGlobalTrigger.PairUnlikeHpt(), recoGlobalTrigger.PairUnlikeLpt(),
		       recoGlobalTrigger.PairLikeHpt(), recoGlobalTrigger.PairLikeLpt(),
		       recoGlobalTrigger.SingleHpt(), recoGlobalTrigger.SingleLpt()};
  
  for (int bit=0;bit<6;bit++){
    if ( recoResp[bit] != inputResp[bit] )
      FillRawsData(histoToFill,bit);
  }
}

//____________________________________________________________________________ 
void AliMUONTriggerQADataMakerRec::FillRatio4434Histos(Int_t evtInterval, Int_t itc, Bool_t isEndOfCycle)
{
  /// Fill ratio 44/34 histos
  TH1* histoEvents = ( isEndOfCycle ) ? GetRawsData(AliMUONQAIndices::kTriggerRawNAnalyzedEvents,itc) : GetMatchingRawsHisto(AliMUONQAIndices::kTriggerRawNAnalyzedEvents,itc);
  if ( ! histoEvents ) return;
  Int_t numEvent = Int_t(histoEvents->GetBinContent(1));

  // Fill every fgkUpdateRatio4434 events
  if (numEvent % evtInterval != 0)
    return;
  
  TH1* histo44dec = ( isEndOfCycle ) ? GetRawsData(AliMUONQAIndices::kTriggerNumberOf44Dec,itc) : GetMatchingRawsHisto(AliMUONQAIndices::kTriggerNumberOf44Dec,itc);
  TH1* histo34dec = ( isEndOfCycle ) ? GetRawsData(AliMUONQAIndices::kTriggerNumberOf34Dec,itc) : GetMatchingRawsHisto(AliMUONQAIndices::kTriggerNumberOf34Dec,itc);
  
  Float_t totalNumberOf44 = histo44dec->GetSumOfWeights();
  Float_t totalNumberOf34 = histo34dec->GetSumOfWeights();

  if ( totalNumberOf34 == 0 )
    return;

  TH1* histoAllEvents = ( isEndOfCycle ) ? GetRawsData(AliMUONQAIndices::kTriggerRatio4434AllEvents,itc) : GetMatchingRawsHisto(AliMUONQAIndices::kTriggerRatio4434AllEvents,itc);
  
  if ( ! histoAllEvents ) return;
  Int_t nbins =  histoAllEvents->GetNbinsX();
  Float_t maxBin = histoAllEvents->GetXaxis()->GetBinLowEdge(nbins+1);

  if ( numEvent - maxBin < 1) return;

  // Use the underflow and overflow to store the number of 34 and 44
  // in previous event
  Float_t previousNumOf34 = histoAllEvents->GetBinContent(0);
  Float_t previousNumOf44 = histoAllEvents->GetBinContent(nbins+1);

  Float_t numOf34Update = totalNumberOf34 - previousNumOf34;
  Float_t numOf44Update = totalNumberOf44 - previousNumOf44;

  // Not enough new tracks since last update
  //if ( numOf34Update == 0 && numOf44Update == 0 )
  if ( numOf34Update < evtInterval - 1 )
    return;

  Int_t newNbins = ( (Int_t)maxBin % fgkUpdateRatio4434 ) ? nbins : nbins+1;
  TString cloneName;
  
  TH1* histoRatioSinceLastUpdate = ( isEndOfCycle ) ? GetRawsData(AliMUONQAIndices::kTriggerRatio4434SinceLastUpdate,itc) : GetMatchingRawsHisto(AliMUONQAIndices::kTriggerRatio4434SinceLastUpdate,itc);

  TH1* histos[2] = {histoAllEvents, histoRatioSinceLastUpdate};
  
  for (Int_t ihisto=0; ihisto<2; ihisto++){
    TH1* currHisto = histos[ihisto];
    cloneName = Form("%sClone", currHisto->GetName());
    TArrayD newAxis(newNbins+1);
    for (Int_t ibin=0; ibin<newNbins; ibin++){
      newAxis[ibin] = currHisto->GetXaxis()->GetBinLowEdge(ibin+1);
    }
    newAxis[newNbins] = numEvent;
    TH1F* copyHisto = (TH1F*)currHisto->Clone(cloneName.Data());
    //currHisto->SetBins(newNbins, 0., fgkUpdateRatio4434*newNbins);
    currHisto->SetBins(newNbins, newAxis.GetArray());
    for (Int_t ibin=1; ibin<newNbins; ibin++){
      currHisto->SetBinContent(ibin, copyHisto->GetBinContent(ibin));
      currHisto->SetBinError(ibin, copyHisto->GetBinError(ibin));
    }
    delete copyHisto;
  }

  Float_t ratio4434 = totalNumberOf44/totalNumberOf34;
  Float_t errorRatio4434 = ProtectedSqrt(totalNumberOf44*(1-ratio4434))/totalNumberOf34;
    
  histoAllEvents->SetBinContent(newNbins,ratio4434);
  histoAllEvents->SetBinError(newNbins,errorRatio4434);

  Float_t ratio4434Update = 0.;
  Float_t errorRatio4434Update = 0.;

  if(numOf34Update!=0){
    ratio4434Update = numOf44Update/numOf34Update;
    if ( numOf44Update > numOf34Update ){
      AliWarning(Form("Number of 4/4 (%f) is higher than number of 3/4 (%f)", numOf44Update, numOf34Update));
    }
    errorRatio4434Update = ProtectedSqrt(numOf44Update*(1-ratio4434Update))/numOf34Update;
  }

  histoRatioSinceLastUpdate->SetBinContent(newNbins,ratio4434Update);
  histoRatioSinceLastUpdate->SetBinError(newNbins,errorRatio4434Update);

  histoAllEvents->SetBinContent(0,totalNumberOf34);
  histoAllEvents->SetBinContent(newNbins+1,totalNumberOf44);

}


//____________________________________________________________________________ 
AliMUONTriggerElectronics* AliMUONTriggerQADataMakerRec::TriggerElectronics()
{
  /// Return trigger electronics
  /// (create it if necessary)
  if ( ! fTriggerProcessor ) 
    fTriggerProcessor = new AliMUONTriggerElectronics(CalibrationData());
  return fTriggerProcessor;
}


//____________________________________________________________________________ 
AliMUONCalibrationData* AliMUONTriggerQADataMakerRec::CalibrationData()
{
  /// Return calibration data
  /// (create it if necessary)
  if ( ! fCalibrationData ) fCalibrationData = new AliMUONCalibrationData(AliCDBManager::Instance()->GetRun());
  return fCalibrationData;
}

//____________________________________________________________________________ 
void AliMUONTriggerQADataMakerRec::ResetDetectorRaws(TObjArray* list)
{
  /// Reset the calibration data
  ResetDetector(list);
  delete fTriggerProcessor;
  fTriggerProcessor = 0x0;
  delete fCalibrationData;
  fCalibrationData = 0x0;
}
