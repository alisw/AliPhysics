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

//____________________________________________________________________________ 
AliMUONTriggerQADataMakerRec::AliMUONTriggerQADataMakerRec(AliQADataMakerRec* master) : 
AliMUONVQADataMakerRec(master),
fDigitMaker(new AliMUONDigitMaker(kFALSE)),
fCalibrationData(new AliMUONCalibrationData(AliCDBManager::Instance()->GetRun())),
fTriggerProcessor(new AliMUONTriggerElectronics(fCalibrationData)),
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

  DisplayTriggerInfo();

  // Normalize RawData histos
  Float_t nbevent = GetRawsData(kRawNAnalyzedEvents)->GetBinContent(1);
  Int_t histoRawsIndex[] = {
    kTriggerErrorSummary,
    kTriggerCalibSummary,
    kTriggerReadOutErrors,
    kTriggerGlobalOutput
  };
  Int_t histoRawsScaledIndex[] = {
    kTriggerErrorSummaryNorm,
    kTriggerCalibSummaryNorm,
    kTriggerReadOutErrorsNorm,
    kTriggerGlobalOutputNorm
  };
  const Int_t kNrawsHistos = sizeof(histoRawsIndex)/sizeof(histoRawsIndex[0]);
  Float_t scaleFactor[kNrawsHistos] = {100., 100., 100., 1.};
  for(Int_t ihisto=0; ihisto<kNrawsHistos; ihisto++){
    TH1* inputHisto = GetRawsData(histoRawsIndex[ihisto]);
    TH1* scaledHisto = GetRawsData(histoRawsScaledIndex[ihisto]);
    if ( scaledHisto && inputHisto &&  nbevent > 0 ) {
      scaledHisto->Reset();
      scaledHisto->Add(inputHisto);
      scaledHisto->Scale(scaleFactor[ihisto]/nbevent);
    }
  } // loop on histos

  TH1* hYCopy = GetRawsData(kTriggerErrorLocalYCopy); //number of YCopy error per board
  TH1* hYCopyTests = GetRawsData(kTriggerErrorLocalYCopyTest); //contains the number of YCopy test per board
  TH1* hYCopyNorm = GetRawsData(kTriggerErrorLocalYCopyNorm); 
  hYCopyNorm->Reset();
  hYCopyNorm->Divide(hYCopy, hYCopyTests, 100., 1.);
     
  Float_t mean = hYCopyNorm->Integral();
      
  TH1* hSummary = GetRawsData(kTriggerErrorSummary);
  hSummary->SetBinContent(kAlgoLocalYCopy+1,mean/192.); //put the mean of the % of YCopy error in the kTriggerError's corresponding bin
}

//____________________________________________________________________________ 
void AliMUONTriggerQADataMakerRec::InitRaws()
{
    /// create Raws histograms in Raws subdir
	
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
    Add2RawsList(histo1D, kTriggerScalersTime, expert, !image, !saveCorr);

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
	Add2RawsList(histo2D, kTriggerScalers + AliMpConstants::NofTriggerChambers()*iCath + iChamber, expert, !image, !saveCorr);
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
	Add2RawsList(histo2D, kTriggerScalersDisplay + AliMpConstants::NofTriggerChambers()*iCath + iChamber, expert, !image, !saveCorr);
      } // loop on chambers
    } // loop on cathodes    

    TString axisLabel[kNtrigCalibSummaryBins] = {"#splitline{Dead}{Channels}", "#splitline{Dead}{Local Boards}", "#splitline{Dead}{Regional Boards}", "#splitline{Dead}{Global Board}", "#splitline{Noisy}{Strips}"};

    TH1F* histoCalib = new TH1F("hTriggerCalibSummaryAll", "MTR calibration summary counts", kNtrigCalibSummaryBins, -0.5, (Float_t)kNtrigCalibSummaryBins - 0.5);
    for (Int_t ibin=1; ibin<=kNtrigCalibSummaryBins; ibin++){
      histoCalib->GetXaxis()->SetBinLabel(ibin, axisLabel[ibin-1].Data());
    }
    histoCalib->SetFillColor(kBlue);
    histoCalib->GetYaxis()->SetTitle("Counts");
    // Copy of previous histo for scaling purposes
    TH1F* histoCalibNorm = (TH1F*)histoCalib->Clone("hTriggerCalibSummary");
    histoCalibNorm->SetTitle("MTR calibration summary");
    histoCalibNorm->SetOption("bar2");
    histoCalibNorm->GetYaxis()->SetTitle("Percentage per event (%)");
    // Adding both histos after cloning to avoid problems with the expert bit
    Add2RawsList(histoCalib,     kTriggerCalibSummary,      expert, !image, !saveCorr);
    Add2RawsList(histoCalibNorm, kTriggerCalibSummaryNorm, !expert,  image, !saveCorr);
  } // Calibration reco param
	
  Char_t *globalXaxisName[6] = {"US HPt", "US LPt", "LS HPt", "LS LPt", "SGL HPt", "SGL LPt"};
  Char_t *allLevelXaxisName[kNtrigAlgoErrorBins] = {"Local algo X", "Local algo Y", "Local LUT","Local Y Copy" , "Local2Regional", "Regional", "Regional2Global", "GlobalFromInGlobal", "GlobalFromInLocal", "GlobalFromOutLocal"};
  Char_t *readoutErrNames[kNtrigStructErrorBins]={"Local","Regional","Global","DARC"};

  TString errorAxisTitle = "Number of errors";

  histo1D = new TH1F("ErrorLocalXPos", "ErrorLocalXPos",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, kTriggerErrorLocalXPos, expert, !image, !saveCorr);

  histo1D = new TH1F("ErrorLocalYPos", "ErrorLocalYPos",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, kTriggerErrorLocalYPos, expert, !image, !saveCorr);

  histo1D = new TH1F("ErrorLocalDev", "ErrorLocalDev",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, kTriggerErrorLocalDev, expert, !image, !saveCorr);

  histo1D = new TH1F("ErrorLocalTriggerDec", "ErrorLocalTriggerDec",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, kTriggerErrorLocalTriggerDec, expert, !image, !saveCorr);

  histo1D = new TH1F("ErrorLocalLPtLSB", "ErrorLocalLPtLSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, kTriggerErrorLocalLPtLSB, expert, !image, !saveCorr);

  histo1D = new TH1F("ErrorLocalLPtMSB", "ErrorLocalLPtMSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, kTriggerErrorLocalLPtMSB, expert, !image, !saveCorr);

  histo1D = new TH1F("ErrorLocalHPtLSB", "ErrorLocalHPtLSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, kTriggerErrorLocalHPtLSB, expert, !image, !saveCorr);

  histo1D = new TH1F("ErrorLocalHPtMSB", "ErrorLocalHPtMSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, kTriggerErrorLocalHPtMSB, expert, !image, !saveCorr);

  histo1D = new TH1F("ErrorLocalTrigY", "ErrorLocalTrigY",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, kTriggerErrorLocalTrigY, expert, !image, !saveCorr);

  histo1D = new TH1F("ErrorLocal2RegionalLPtLSB", "ErrorLocal2RegionalLPtLSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, kTriggerErrorLocal2RegionalLPtLSB, expert, !image, !saveCorr);

  histo1D = new TH1F("ErrorLocal2RegionalLPtMSB", "ErrorLocal2RegionalLPtMSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, kTriggerErrorLocal2RegionalLPtMSB, expert, !image, !saveCorr);

  histo1D = new TH1F("ErrorLocal2RegionalHPtLSB", "ErrorLocal2RegionalHPtLSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, kTriggerErrorLocal2RegionalHPtLSB, expert, !image, !saveCorr);

  histo1D = new TH1F("ErrorLocal2RegionalHPtMSB", "ErrorLocal2RegionalHPtMSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histo1D->GetXaxis()->SetTitle(boardName.Data());
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(histo1D, kTriggerErrorLocal2RegionalHPtMSB, expert, !image, !saveCorr);

  histo1D = new TH1F("ErrorOutGlobalFromInGlobal", "ErrorOutGlobalFromInGlobal",6,-0.5,6-0.5);
  histo1D->GetYaxis()->SetTitle(errorAxisTitle.Data());
  for (int ibin=0;ibin<6;ibin++){
    histo1D->GetXaxis()->SetBinLabel(ibin+1,globalXaxisName[ibin]);
  }
  Add2RawsList(histo1D, kTriggerErrorOutGlobalFromInGlobal, expert, !image, !saveCorr);

  TH1F* histoAlgoErr = new TH1F("hTriggerAlgoNumOfErrors", "Trigger Algorithm total errors",kNtrigAlgoErrorBins,-0.5,(Float_t)kNtrigAlgoErrorBins-0.5);
  histoAlgoErr->GetYaxis()->SetTitle("Number of events with errors");
  for (int ibin=0;ibin<kNtrigAlgoErrorBins;ibin++){
    histoAlgoErr->GetXaxis()->SetBinLabel(ibin+1,allLevelXaxisName[ibin]);
  }
  histoAlgoErr->SetFillColor(kBlue);
  // Copy of previous histo for scaling purposes
  TH1F* histoAlgoErrNorm = (TH1F*)histoAlgoErr->Clone("hTriggerAlgoErrors");
  histoAlgoErrNorm->SetOption("bar2");
  histoAlgoErrNorm->SetTitle("Trigger algorithm errors");
  histoAlgoErrNorm->GetYaxis()->SetTitle("% of events with errors");
  // Adding both histos after cloning to avoid problems with the expert bit
  Add2RawsList(histoAlgoErr,     kTriggerErrorSummary,      expert, !image, !saveCorr);
  Add2RawsList(histoAlgoErrNorm, kTriggerErrorSummaryNorm, !expert,  image, !saveCorr);  

  histo1D = new TH1F("hTriggeredBoards", "Triggered boards", nbLocalBoard, 0.5, (Float_t)nbLocalBoard + 0.5);
  Add2RawsList(histo1D, kTriggeredBoards, expert, !image, !saveCorr);

  histo2D = (TH2F*)triggerDisplay.GetEmptyDisplayHisto("hFiredBoardsDisplay", AliMUONTriggerDisplay::kDisplayBoards,
						       0, 0, "Local board triggers / event");
  histo2D->SetOption("COLZ");
  Add2RawsList(histo2D, kTriggerBoardsDisplay, expert, !image, !saveCorr);

  TH1F* histoYCopyErr = new TH1F("ErrorLocalYCopy", "Number of YCopy errors",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  histoYCopyErr->GetXaxis()->SetTitle(boardName.Data());
  histoYCopyErr->GetYaxis()->SetTitle(errorAxisTitle.Data());
  // Copy of previous histo for scaling purposes
  TH1F* histoYCopyErrTest = (TH1F*)histoYCopyErr->Clone("ErrorLocalYCopyTest");
  histoYCopyErrTest->SetTitle("Number of YCopy tested");
  // Copy of previous histo for scaling purposes
  TH1F* histoYCopyErrNorm = (TH1F*)histoYCopyErr->Clone("ErrorLocalYCopyNorm");
  histoYCopyErrNorm->SetTitle("% of YCopy errors");
  // Adding both histos after cloning to avoid problems with the expert bit
  Add2RawsList(histoYCopyErr,     kTriggerErrorLocalYCopy,     expert, !image, !saveCorr);
  Add2RawsList(histoYCopyErrTest, kTriggerErrorLocalYCopyTest, expert, !image, !saveCorr);
  Add2RawsList(histoYCopyErrNorm, kTriggerErrorLocalYCopyNorm, expert, !image, !saveCorr);

  TH1F* histoROerr = new TH1F("hTriggerReadoutNumOfErrors","Trigger Read-Out total errors", kNtrigStructErrorBins, -0.5, (Float_t)kNtrigStructErrorBins-0.5);
  histoROerr->GetYaxis()->SetTitle("Fraction of errors");
  histoROerr->SetFillColor(kBlue);
  for (int ibin=0;ibin<kNtrigStructErrorBins;ibin++){
    histoROerr->GetXaxis()->SetBinLabel(ibin+1,readoutErrNames[ibin]);
  }
  // Copy of previous histo for scaling purposes
  TH1F* histoROerrNorm = (TH1F*)histoROerr->Clone("hTriggerReadoutErrors");
  histoROerrNorm->SetTitle("Trigger Read-Out errors");
  histoROerrNorm->SetOption("bar2");
  histoROerrNorm->GetYaxis()->SetTitle("% of errors per event");
  // Adding both histos after cloning to avoid problems with the expert bit
  Add2RawsList(histoROerr,     kTriggerReadOutErrors,      expert, !image, !saveCorr);
  Add2RawsList(histoROerrNorm, kTriggerReadOutErrorsNorm, !expert,  image, !saveCorr);

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
  histoGlobalMultNorm->SetOption("bar2");
  histoGlobalMultNorm->SetBarWidth(0.5);
  histoGlobalMultNorm->SetBarOffset(0.25);
  histoGlobalMultNorm->GetYaxis()->SetTitle("Triggers per event");
  // Adding both histos after cloning to avoid problems with the expert bit
  Add2RawsList(histoGlobalMult,     kTriggerGlobalOutput,     expert, !image, !saveCorr);
  Add2RawsList(histoGlobalMultNorm, kTriggerGlobalOutputNorm, expert, !image, !saveCorr);

  histo1D = new TH1F("hRawNAnalyzedEvents", "Number of analyzed events per specie", 1, 0.5, 1.5);
  Int_t esindex = AliRecoParam::AConvert(CurrentEventSpecie());
  histo1D->GetXaxis()->SetBinLabel(1, AliRecoParam::GetEventSpecieName(esindex));
  histo1D->GetYaxis()->SetTitle("Number of analyzed events");
  Add2RawsList(histo1D, kRawNAnalyzedEvents, expert, !image, !saveCorr);
}

//__________________________________________________________________
void AliMUONTriggerQADataMakerRec::InitDigits() 
{
  /// Initialized Digits spectra 
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1I* h0 = new TH1I("hDigitsDetElem", "Detection element distribution in Digits;Detection element Id;Counts",  400, 1100, 1500); 
  Add2DigitsList(h0, 0, !expert, image);
} 

//____________________________________________________________________________ 
void AliMUONTriggerQADataMakerRec::InitRecPoints()
{
	/// create Reconstructed Points histograms in RecPoints subdir for the
	/// MUON Trigger subsystem.

  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  TH1F* histo1D = 0x0;

  histo1D = new TH1F("hNAnalyzedEvents", "Number of analyzed events per specie", 1, 0.5, 1.5);
  Int_t esindex = AliRecoParam::AConvert(CurrentEventSpecie());
  histo1D->GetXaxis()->SetBinLabel(1, AliRecoParam::GetEventSpecieName(esindex));
  histo1D->GetYaxis()->SetTitle("Number of analyzed events");
  Add2RecPointsList(histo1D, kNAnalyzedEvents, expert, !image);

  histo1D = new TH1F("hTriggerTrippedChambers", "Trigger RPCs in trip", 418, 1100-0.5, 1417+0.5);
  histo1D->GetXaxis()->SetTitle("DetElemId");
  histo1D->GetYaxis()->SetTitle("# of trips");
  histo1D->SetFillColor(kRed);
  histo1D->SetLineColor(kRed);
  Add2RecPointsList(histo1D, kTriggerRPCtrips, !expert, image);

  FillTriggerDCSHistos();	
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
	
    GetRawsData(kRawNAnalyzedEvents)->Fill(1.);

    // Init Local/Regional/Global decision with fake values

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
    digitStore.Create();
    digitStore.Clear();

    AliMUONDigitStoreV2R digitStoreAll;
    digitStoreAll.Create();
    digitStoreAll.Clear();
    TArrayS xyPatternAll[2];
    for(Int_t icath=0; icath<AliMpConstants::NofCathodes(); icath++){
      xyPatternAll[icath].Set(AliMpConstants::NofTriggerChambers());
      xyPatternAll[icath].Reset(1);
    }
    
    AliMUONTriggerStoreV1 recoTriggerStore;
    recoTriggerStore.Create();
    recoTriggerStore.Clear();

    AliMUONTriggerStoreV1 inputTriggerStore;
    inputTriggerStore.Create();
    inputTriggerStore.Clear();

    AliMUONGlobalTrigger inputGlobalTrigger;

    UShort_t maxNcounts = 0xFFFF;
    
    // Get trigger Local, Regional, Global in/outputs and scalers

    Int_t loCircuit=0;
    AliMpCDB::LoadDDLStore();

    const AliMUONRawStreamTriggerHP::AliHeader*          darcHeader  = 0x0;
    const AliMUONRawStreamTriggerHP::AliRegionalHeader*  regHeader   = 0x0;
    const AliMUONRawStreamTriggerHP::AliLocalStruct*     localStruct = 0x0;

    Int_t nDeadLocal = 0, nDeadRegional = 0, nDeadGlobal = 0, nNoisyStrips = 0;

    // When a crate is not present, the loop on boards is not performed
    // This should allow to correctly count the local boards
    Int_t countNotifiedBoards = 0, countAllBoards = 0;

    AliMUONRawStreamTriggerHP rawStreamTrig(rawReader);
    while (rawStreamTrig.NextDDL()) 
    {
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
	  ((TH1F*)GetRawsData(kTriggerScalersTime))->Fill(1., nOfSeconds);
	}

	//Get Global datas
	inputGlobalTrigger.SetFromGlobalResponse(darcHeader->GetGlobalOutput());
	Bool_t resp[6] = {inputGlobalTrigger.PairUnlikeHpt(), inputGlobalTrigger.PairUnlikeLpt(),
			  inputGlobalTrigger.PairLikeHpt(), inputGlobalTrigger.PairLikeLpt(),
			  inputGlobalTrigger.SingleHpt(), inputGlobalTrigger.SingleHpt()}; 
	for (Int_t bit=0; bit<6; bit++){
	  if ( ! resp[bit] ){
	    if ( fillScalerHistos )
	      nDeadGlobal++;
	  }
	  else
	    ((TH1F*)GetRawsData(kTriggerGlobalOutput))->Fill(bit);
	}

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
	  if ( fillScalerHistos ) // Compute total number of strips
	    fDigitMaker->TriggerDigits(loCircuit, xyPatternAll, digitStoreAll);

	  //Get electronic Decisions from data

	  //Get regional inputs -> not checked, hardware read-out doesn't work
	  //fTriggerInputRegionalDataLPt[0][loCircuit]=Int_t(((regHeader->GetInput(0))>>(2*iLocal))&1);
	  //fTriggerInputRegionalDataLPt[1][loCircuit]=Int_t(((regHeader->GetInput(1))>>((2*iLocal)+1))&1);

	  //Get local in/outputs
	  if (Int_t(localStruct->GetDec())!=0){
	      ((TH1F*)GetRawsData(kTriggeredBoards))->Fill(loCircuit);
	  }
	  else if ( fillScalerHistos ){
	    nDeadLocal++;
	  }

	  // loop over strips
	  if ( fillScalerHistos ) {
	    Int_t cathode = localStruct->GetComptXY()%2;
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
		if ( scalerVal[ich] > 0 )
		  ((TH2F*)GetRawsData(kTriggerScalers + AliMpConstants::NofTriggerChambers()*cathode + ich))
		    ->Fill(loCircuit, ibitxy, 2*(Float_t)scalerVal[ich]);

		if ( scalerVal[ich] >= maxNcounts )
		  nNoisyStrips++;
	      } // loop on chamber
	    } // loop on strips
	  } // scaler event
	} // iLocal
	if ( nBoardsInReg == 0 )
	  nDeadRegional++; // Not necessary when regional output will work
      } // iReg

      Float_t readoutErrors[kNtrigStructErrorBins] = {
	((Float_t)rawStreamTrig.GetLocalEoWErrors())/((Float_t)countAllBoards),
	((Float_t)rawStreamTrig.GetRegEoWErrors())/16.,
	((Float_t)rawStreamTrig.GetGlobalEoWErrors())/6.,
	((Float_t)rawStreamTrig.GetDarcEoWErrors())/2.
      };
    
      for (Int_t ibin=0; ibin<kNtrigStructErrorBins; ibin++){
	if ( readoutErrors[ibin] > 0 )
	  ((TH1F*)GetRawsData(kTriggerReadOutErrors))->Fill(ibin, readoutErrors[ibin]);
      }
    } // NextDDL

    nDeadLocal += AliMUONConstants::NTriggerCircuit() - countNotifiedBoards;
    Int_t nStripsTot = digitStoreAll.GetSize();
    if ( nStripsTot > 0 ) { // The value is != 0 only for scaler events
      Float_t fraction[kNtrigCalibSummaryBins] = {
	((Float_t)(nStripsTot - digitStore.GetSize())) / ((Float_t)nStripsTot),
	//(Float_t)nDeadLocal / ((Float_t)countNotifiedBoards),
	(Float_t)nDeadLocal / ((Float_t)AliMUONConstants::NTriggerCircuit()),
	(Float_t)nDeadRegional / 16.,
	(Float_t)nDeadGlobal / 6., // Number of bits of global response
	(Float_t)nNoisyStrips / ((Float_t)nStripsTot),
      };

      for(Int_t ibin = 0; ibin < kNtrigCalibSummaryBins; ibin++){
	if ( fraction[ibin] > 0. )
	  ((TH1F*)GetRawsData(kTriggerCalibSummary))->Fill(ibin, fraction[ibin]);
      }
    }

  fTriggerProcessor->Digits2Trigger(digitStore,recoTriggerStore);

  //Reconstruct Global decision from Global inputs
  UChar_t recoResp = RawTriggerInGlobal2OutGlobal(globalInput);
  AliMUONGlobalTrigger recoGlobalTrigger;
  recoGlobalTrigger.SetFromGlobalResponse(recoResp);

  // Compare data and reconstructed decisions and fill histos
  RawTriggerMatchOutLocal(inputTriggerStore, recoTriggerStore);
  //RawTriggerMatchOutLocalInRegional(); // Not tested, hardware read-out doesn't work
  RawTriggerMatchOutGlobalFromInGlobal(inputGlobalTrigger, recoGlobalTrigger);
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
    GetDigitsData(0)->Fill(dig->DetElemId());
    GetDigitsData(1)->Fill(dig->ADC());
    }
}

//____________________________________________________________________________
void AliMUONTriggerQADataMakerRec::MakeRecPoints(TTree* /*clustersTree*/)
{
  /// Fill histogram with total number of analyzed events for normalization purposes

  // Do nothing in case of calibration event
  if ( GetRecoParam()->GetEventSpecie() == AliRecoParam::kCalib ) return;
	
  GetRecPointsData(kNAnalyzedEvents)->Fill(1.);
}

//____________________________________________________________________________
void AliMUONTriggerQADataMakerRec::MakeESDs(AliESDEvent* /*esd*/)
{  
  /// Empty implementation
}


//____________________________________________________________________________ 
void AliMUONTriggerQADataMakerRec::DisplayTriggerInfo()
{
  //
  /// Display trigger information in a user-friendly way:
  /// from local board and strip numbers to their position on chambers
  //

  AliMUONTriggerDisplay triggerDisplay;
  
  TH2F* histoStrips=0x0;
  TH2F* histoDisplayStrips=0x0;
  if ( GetRawsData(kTriggerScalers) ) {
    AliMUONTriggerDisplay::EDisplayOption displayOption = AliMUONTriggerDisplay::kNormalizeToArea;
    for (Int_t iCath = 0; iCath < AliMpConstants::NofCathodes(); iCath++)
      {    
	for (Int_t iChamber = 0; iChamber < AliMpConstants::NofTriggerChambers(); iChamber++)
	  {
	    histoStrips = (TH2F*)GetRawsData(kTriggerScalers + AliMpConstants::NofTriggerChambers()*iCath + iChamber);

	    if(histoStrips->GetEntries()==0) continue; // No events found => No need to display

	    histoDisplayStrips = (TH2F*)GetRawsData(kTriggerScalersDisplay + AliMpConstants::NofTriggerChambers()*iCath + iChamber);

	    triggerDisplay.FillDisplayHistogram(histoStrips, histoDisplayStrips,
						AliMUONTriggerDisplay::kDisplayStrips, iCath, iChamber, displayOption);

	    Float_t scaleValue = ((TH1F*)GetRawsData(kTriggerScalersTime))->GetBinContent(1);
	    if(scaleValue>0.) histoDisplayStrips->Scale(1./scaleValue);
	  } // iChamber
      } // iCath
  }

  if ( GetRawsData(kTriggeredBoards) ){
    TH1F* histoBoards = (TH1F*)GetRawsData(kTriggeredBoards);
    TH2F* histoDisplayBoards = (TH2F*)GetRawsData(kTriggerBoardsDisplay);
    triggerDisplay.FillDisplayHistogram(histoBoards, histoDisplayBoards, AliMUONTriggerDisplay::kDisplayBoards, 0, 0);
    Float_t scaleValue = GetRawsData(kRawNAnalyzedEvents)->GetBinContent(1);
    if(scaleValue>0.) histoDisplayBoards->Scale(1./scaleValue);
  }
}


//_____________________________________________________________________________
Bool_t 
AliMUONTriggerQADataMakerRec::FillTriggerDCSHistos()
{
  /// Get HV and currents values for one trigger chamber
  
  AliCodeTimerAuto("",0);

  TMap* triggerDcsMap = fCalibrationData->TriggerDCS();

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

      histoIndex = kTriggerRPChv + ich;
      histoName = Form("hRPCHVChamber%i", 11+ich);
      histoTitle = Form("Chamber %i: RPC HV (kV)", 11+ich);

      currHisto = (TH2F*)GetRecPointsData(histoIndex);

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
      if ( isTrip ) ((TH1F*)GetRecPointsData(kTriggerRPCtrips))->Fill(detElemId);
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
    printf(Form("Did not find expected alias (%s) for DE %d\n",
		currAlias.Data(),detElemId));
    return 0x0;
  }

  TObjArray* values = static_cast<TObjArray*>(triggerDcsPair->Value());
  if (!values)
  {
    printf(Form("Could not get values for alias %s\n",currAlias.Data()));
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

    AliMUONGlobalCrateConfig* globalConfig = fCalibrationData->GlobalTriggerCrateConfig();

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
void AliMUONTriggerQADataMakerRec::RawTriggerMatchOutLocal(AliMUONVTriggerStore& inputTriggerStore,
							   AliMUONVTriggerStore& recoTriggerStore)
{
  //
  /// Match data and reconstructed Local Trigger decision
  //

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

    ((TH1F*)GetRawsData(kTriggerErrorLocalYCopyTest))->Fill(loCircuit);
  
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
      ((TH1F*)GetRawsData(kTriggerErrorLocalYCopy))->Fill(loCircuit);    
  } // loop on local boards

  if (errorInYCopy)
    ((TH1F*)GetRawsData(kTriggerErrorSummary))->Fill(kAlgoLocalYCopy);
  
  Bool_t errorInXPosDev = kFALSE;
  Bool_t errorInYPosTrigY = kFALSE;
  Bool_t errorInLUT = kFALSE;

  next.Reset();
  while ( ( recoLocalTrigger = static_cast<AliMUONLocalTrigger*>(next()) ) )
  {  
    loCircuit = recoLocalTrigger->LoCircuit();
    Int_t iboard = loCircuit - 1;
  
    inputLocalTrigger = inputTriggerStore.FindLocal(loCircuit);

    if ( recoLocalTrigger->LoStripX() != inputLocalTrigger->LoStripX() ) {
      ((TH1F*)GetRawsData(kTriggerErrorLocalXPos))->Fill(loCircuit);
      errorInXPosDev = kTRUE;
    }

    if ( recoLocalTrigger->GetDeviation() != inputLocalTrigger->GetDeviation() ) {
      ((TH1F*)GetRawsData(kTriggerErrorLocalDev))->Fill(loCircuit);
      errorInXPosDev = kTRUE;
    }

    // Skip following checks in case we previously found YCopy error and YPos or trigY errors
    if ( (!skipBoard[iboard]) || ( (recoLocalTrigger->LoStripY() == inputLocalTrigger->LoStripY()) && (recoLocalTrigger->LoTrigY() == inputLocalTrigger->LoTrigY())) ) {
	
	if ( recoLocalTrigger->GetLoDecision() != inputLocalTrigger->GetLoDecision() ) {
	    ((TH1F*)GetRawsData(kTriggerErrorLocalTriggerDec))->Fill(loCircuit);
	}
	
	// Test Hpt and LPT
	Int_t recoLut[2]  = { recoLocalTrigger->LoLpt(),  recoLocalTrigger->LoHpt() };
	Int_t inputLut[2] = {inputLocalTrigger->LoLpt(), inputLocalTrigger->LoHpt() };
	Int_t currIndex[2][2] = {{kTriggerErrorLocalLPtLSB, kTriggerErrorLocalLPtMSB},
				 {kTriggerErrorLocalHPtMSB, kTriggerErrorLocalHPtMSB}};
	for (Int_t ilut=0; ilut<2; ilut++){
	    Int_t bitDiff = recoLut[ilut]^inputLut[ilut];
	    if ( bitDiff == 0 ) continue;
	    for (Int_t ibit=0; ibit<2; ibit++){
		Bool_t isBitDifferent = (bitDiff>>ibit)&1;
		if ( isBitDifferent ){
		    ((TH1F*)GetRawsData(currIndex[ilut][ibit]))->Fill(loCircuit);
		    errorInLUT = kTRUE;
		}
	    }
	}
    }
    
 
    // Skip following checks in case we previously found YCopy errors
    if ( skipBoard[iboard] ) continue;

    if ( recoLocalTrigger->LoStripY() != inputLocalTrigger->LoStripY() ) {
      ((TH1F*)GetRawsData(kTriggerErrorLocalYPos))->Fill(loCircuit);
      errorInYPosTrigY = kTRUE;
    }

    if ( recoLocalTrigger->LoTrigY() != inputLocalTrigger->LoTrigY()  ) {
      ((TH1F*)GetRawsData(kTriggerErrorLocalTrigY))->Fill(loCircuit);	
      errorInYPosTrigY = kTRUE;
    }
  } // loop on local boards

  if (errorInXPosDev)
    ((TH1F*)GetRawsData(kTriggerErrorSummary))->Fill(kAlgoLocalX);

  if (errorInLUT)
    ((TH1F*)GetRawsData(kTriggerErrorSummary))->Fill(kAlgoLocalLUT);

  if (errorInYPosTrigY)
    ((TH1F*)GetRawsData(kTriggerErrorSummary))->Fill(kAlgoLocalY);

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
void AliMUONTriggerQADataMakerRec::RawTriggerMatchOutGlobalFromInGlobal(AliMUONGlobalTrigger& inputGlobalTrigger, 
									AliMUONGlobalTrigger& recoGlobalTrigger)
{
  //
  /// Match data and reconstructed Global Trigger decision for a reconstruction from Global inputs
  //

  if ( recoGlobalTrigger.GetGlobalResponse() == inputGlobalTrigger.GetGlobalResponse() )
    return;

  ((TH1F*)GetRawsData(kTriggerErrorSummary))->Fill(kAlgoGlobalFromGlobal);

  Bool_t inputResp[6] = {inputGlobalTrigger.PairUnlikeHpt(), inputGlobalTrigger.PairUnlikeLpt(),
			 inputGlobalTrigger.PairLikeHpt(), inputGlobalTrigger.PairLikeLpt(),
			 inputGlobalTrigger.SingleHpt(), inputGlobalTrigger.SingleHpt()};

  Bool_t recoResp[6] = {recoGlobalTrigger.PairUnlikeHpt(), recoGlobalTrigger.PairUnlikeLpt(),
			recoGlobalTrigger.PairLikeHpt(), recoGlobalTrigger.PairLikeLpt(),
			recoGlobalTrigger.SingleHpt(), recoGlobalTrigger.SingleHpt()};

  for (int bit=0;bit<6;bit++){
    if ( recoResp[bit] != inputResp[bit] )
      ((TH1F*)GetRawsData(kTriggerErrorOutGlobalFromInGlobal))->Fill(bit);
  }
}
