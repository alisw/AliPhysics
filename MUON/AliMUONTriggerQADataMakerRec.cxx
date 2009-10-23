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
/// \author C. Finck, D. Stocco, L. Aphecetche

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

  // Normalize RawData histos
  Float_t nbevent = GetRawsData(kRawNAnalyzedEvents)->GetBinContent(1);
  Int_t histoRawsIndex[] = {
    kTriggerError,
    kTriggerCalibSummary,
    kTriggerReadOutErrors,
    kTriggerGlobalOutput
  };
  const Int_t kNrawsHistos = sizeof(histoRawsIndex)/sizeof(histoRawsIndex[0]);
  Float_t scaleFactor[kNrawsHistos] = {100., 100., 100., 1.};
  for(Int_t ihisto=0; ihisto<kNrawsHistos; ihisto++){
    TH1* currHisto = GetRawsData(histoRawsIndex[ihisto]);
    if ( currHisto && nbevent > 0 ){
      currHisto->Scale(scaleFactor[ihisto]/nbevent);
    }
  }
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
    
    histo1D = new TH1F("hTriggerScalersTime", "Acquisition time from trigger scalers", 1, 0.5, 1.5);
    histo1D->GetXaxis()->SetBinLabel(1, "One-bin histogram: bin is filled at each scaler event.");
    histo1D->GetYaxis()->SetTitle("Cumulated scaler time (s)");
    Add2RawsList(histo1D, kTriggerScalersTime, expert, !image, !saveCorr);

    TString axisLabel[kNtrigCalibSummaryBins] = {"#splitline{Dead}{Channels}", "#splitline{Dead}{Local Boards}", "#splitline{Dead}{Regional Boards}", "#splitline{Dead}{Global Board}", "#splitline{Noisy}{Strips}"};

    histo1D = new TH1F("hTriggerCalibSummary", "MTR calibration sumamry", kNtrigCalibSummaryBins, -0.5, (Float_t)kNtrigCalibSummaryBins - 0.5);
    for (Int_t ibin=1; ibin<=kNtrigCalibSummaryBins; ibin++){
      histo1D->GetXaxis()->SetBinLabel(ibin, axisLabel[ibin-1].Data());
    }
    histo1D->GetYaxis()->SetTitle("Percentage per event (%)");
    histo1D->SetOption("bar2");
    histo1D->SetStats(kFALSE);
    histo1D->SetFillColor(kRed);
    Add2RawsList(histo1D, kTriggerCalibSummary, !expert, image, !saveCorr);
  } // Calibration reco param

  histo1D = new TH1F("hTriggeredBoards", "Triggered boards", nbLocalBoard, 0.5, (Float_t)nbLocalBoard + 0.5);
  Add2RawsList(histo1D, kTriggeredBoards, expert, !image, !saveCorr);

  histo2D = (TH2F*)triggerDisplay.GetEmptyDisplayHisto("hFiredBoardsDisplay", AliMUONTriggerDisplay::kDisplayBoards,
						       0, 0, "Local board triggers / event");
  histo2D->SetOption("COLZ");
  Add2RawsList(histo2D, kTriggerBoardsDisplay, expert, !image, !saveCorr);
	
  Char_t *globalXaxisName[6] = {"US HPt", "US LPt", "LS HPt", "LS LPt", "SGL HPt", "SGL LPt"};
  Char_t *allLevelXaxisName[kNtrigAlgoErrorBins] = {"Local algo X", "Local algo Y", "Local LUT","Local Y Copy" , "Local2Regional", "Regional", "Regional2Global", "GlobalFromInGlobal", "GlobalFromInLocal", "GlobalFromOutLocal"};
  Char_t *readoutErrNames[kNtrigStructErrorBins]={"Local","Regional","Global","DARC"};

  TString errorAxisTitle = "Number of errors";

  TH1F* h11 = new TH1F("ErrorLocalXPos", "ErrorLocalXPos",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  h11->GetXaxis()->SetTitle(boardName.Data());
  h11->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(h11, kTriggerErrorLocalXPos, expert, !image, !saveCorr);

  TH1F* h12 = new TH1F("ErrorLocalYPos", "ErrorLocalYPos",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  h12->GetXaxis()->SetTitle(boardName.Data());
  h12->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(h12, kTriggerErrorLocalYPos, expert, !image, !saveCorr);

  TH1F* h13 = new TH1F("ErrorLocalDev", "ErrorLocalDev",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  h13->GetXaxis()->SetTitle(boardName.Data());
  h13->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(h13, kTriggerErrorLocalDev, expert, !image, !saveCorr);

  TH1F* h14 = new TH1F("ErrorLocalTriggerDec", "ErrorLocalTriggerDec",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  h14->GetXaxis()->SetTitle(boardName.Data());
  h14->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(h14, kTriggerErrorLocalTriggerDec, expert, !image, !saveCorr);

  TH1F* h15 = new TH1F("ErrorLocalLPtLSB", "ErrorLocalLPtLSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  h15->GetXaxis()->SetTitle(boardName.Data());
  h15->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(h15, kTriggerErrorLocalLPtLSB, expert, !image, !saveCorr);

  TH1F* h16 = new TH1F("ErrorLocalLPtMSB", "ErrorLocalLPtMSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  h16->GetXaxis()->SetTitle(boardName.Data());
  h16->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(h16, kTriggerErrorLocalLPtMSB, expert, !image, !saveCorr);

  TH1F* h17 = new TH1F("ErrorLocalHPtLSB", "ErrorLocalHPtLSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  h17->GetXaxis()->SetTitle(boardName.Data());
  h17->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(h17, kTriggerErrorLocalHPtLSB, expert, !image, !saveCorr);

  TH1F* h18 = new TH1F("ErrorLocalHPtMSB", "ErrorLocalHPtMSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  h18->GetXaxis()->SetTitle(boardName.Data());
  h18->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(h18, kTriggerErrorLocalHPtMSB, expert, !image, !saveCorr);

  TH1F* h19 = new TH1F("ErrorLocal2RegionalLPtLSB", "ErrorLocal2RegionalLPtLSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  h19->GetXaxis()->SetTitle(boardName.Data());
  h19->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(h19, kTriggerErrorLocal2RegionalLPtLSB, expert, !image, !saveCorr);

  TH1F* h20 = new TH1F("ErrorLocal2RegionalLPtMSB", "ErrorLocal2RegionalLPtMSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  h20->GetXaxis()->SetTitle(boardName.Data());
  h20->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(h20, kTriggerErrorLocal2RegionalLPtMSB, expert, !image, !saveCorr);

  TH1F* h21 = new TH1F("ErrorLocal2RegionalHPtLSB", "ErrorLocal2RegionalHPtLSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  h21->GetXaxis()->SetTitle(boardName.Data());
  h21->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(h21, kTriggerErrorLocal2RegionalHPtLSB, expert, !image, !saveCorr);

  TH1F* h22 = new TH1F("ErrorLocal2RegionalHPtMSB", "ErrorLocal2RegionalHPtMSB",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  h22->GetXaxis()->SetTitle(boardName.Data());
  h22->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(h22, kTriggerErrorLocal2RegionalHPtMSB, expert, !image, !saveCorr);

  TH1F* h23 = new TH1F("ErrorOutGlobalFromInGlobal", "ErrorOutGlobalFromInGlobal",6,-0.5,6-0.5);
  h23->GetYaxis()->SetTitle(errorAxisTitle.Data());
  for (int ibin=0;ibin<6;ibin++){
    h23->GetXaxis()->SetBinLabel(ibin+1,globalXaxisName[ibin]);
  }
  Add2RawsList(h23, kTriggerErrorOutGlobalFromInGlobal, expert, !image, !saveCorr);

  TH1F* h24 = new TH1F("hTriggerAlgoErrors", "Trigger Algorithm errors",kNtrigAlgoErrorBins,-0.5,(Float_t)kNtrigAlgoErrorBins-0.5);
  h24->GetYaxis()->SetTitle("% of error");
  h24->SetOption("bar2");
  for (int ibin=0;ibin<kNtrigAlgoErrorBins;ibin++){
    h24->GetXaxis()->SetBinLabel(ibin+1,allLevelXaxisName[ibin]);
  }
  h24->SetFillColor(2);
  Add2RawsList(h24, kTriggerError, !expert, image, !saveCorr);

  TH1F* h25 = new TH1F("ErrorLocalTrigY", "ErrorLocalTrigY",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  h25->GetXaxis()->SetTitle(boardName.Data());
  h25->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(h25, kTriggerErrorLocalTrigY, expert, !image, !saveCorr);

  TH1F* h26 = new TH1F("ErrorLocalYCopy", "ErrorLocalYCopy",nbLocalBoard,0.5,(Float_t)nbLocalBoard+0.5);
  h26->GetXaxis()->SetTitle(boardName.Data());
  h26->GetYaxis()->SetTitle(errorAxisTitle.Data());
  Add2RawsList(h26, kTriggerErrorLocalYCopy, expert, !image, !saveCorr);
	
  TH1F* h27 = new TH1F("hRawNAnalyzedEvents", "Number of analyzed events per specie", 1, 0.5, 1.5);
  Int_t esindex = AliRecoParam::AConvert(CurrentEventSpecie());
  h27->GetXaxis()->SetBinLabel(1, AliRecoParam::GetEventSpecieName(esindex));
  h27->GetYaxis()->SetTitle("Number of analyzed events");
  Add2RawsList(h27, kRawNAnalyzedEvents, expert, !image, !saveCorr);

  TH1F* h28 = new TH1F("hTriggerReadoutErrors","Trigger Read-Out errors", kNtrigStructErrorBins, -0.5, (Float_t)kNtrigStructErrorBins-0.5);
  h28->SetOption("bar2");
  h28->GetYaxis()->SetTitle("% of errors");
  for (int ibin=0;ibin<kNtrigStructErrorBins;ibin++){
    h28->GetXaxis()->SetBinLabel(ibin+1,readoutErrNames[ibin]);
  }
  h28->SetFillColor(2);
  Add2RawsList(h28, kTriggerReadOutErrors, !expert, image, !saveCorr);

  TH1F* h29 = new TH1F("hTriggerGlobalOutMultiplicity","Trigger global outputs multiplicity", 6, -0.5, 6.-0.5);
  h29->SetOption("bar2");
  h29->GetYaxis()->SetTitle("Number of triggers per event"); 
  h29->GetXaxis()->SetTitle("Global output");
  for (int ibin=0;ibin<6;ibin++){
    h29->GetXaxis()->SetBinLabel(ibin+1,globalXaxisName[ibin]);
  }        
  h29->SetFillColor(3);
  Add2RawsList(h29, kTriggerGlobalOutput, expert, !image, !saveCorr);
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
}

//____________________________________________________________________________
void AliMUONTriggerQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
	/// make QA for rawdata trigger
	
    GetRawsData(kRawNAnalyzedEvents)->Fill(1.);

    // Init Local/Regional/Global decision with fake values
 
    Int_t globaltemp[4];
    for (Int_t bit=0; bit<4; bit++){
	globaltemp[bit]=0;
	fgitmp[bit]=0;
    }

    for (Int_t loc=0;loc<235;loc++){
	fTriggerErrorLocalYCopy[loc]=kFALSE;

	fTriggerOutputLocalRecTriggerDec[loc]=0;
	fTriggerOutputLocalRecLPtDec[0][loc]=0;
	fTriggerOutputLocalRecLPtDec[1][loc]=0;
	fTriggerOutputLocalRecHPtDec[0][loc]=0;
	fTriggerOutputLocalRecHPtDec[1][loc]=0;
	fTriggerOutputLocalRecXPos[loc]=0;
	fTriggerOutputLocalRecYPos[loc]=15;
	fTriggerOutputLocalRecDev[loc]=0;
	fTriggerOutputLocalRecTrigY[loc]=1;

	fTriggerOutputLocalDataTriggerDec[loc]=0;
	fTriggerOutputLocalDataLPtDec[0][loc]=0;
	fTriggerOutputLocalDataLPtDec[1][loc]=0;
	fTriggerOutputLocalDataHPtDec[0][loc]=0;
	fTriggerOutputLocalDataHPtDec[1][loc]=0;
	fTriggerOutputLocalDataXPos[loc]=0;
	fTriggerOutputLocalDataYPos[loc]=15;
	fTriggerOutputLocalDataDev[loc]=0;
	fTriggerOutputLocalDataTrigY[loc]=1;
	fTriggerInputRegionalDataLPt[0][loc]=0;
	fTriggerInputRegionalDataLPt[1][loc]=0;
	fTriggerInputRegionalDataHPt[0][loc]=0;
	fTriggerInputRegionalDataHPt[1][loc]=0;	
    }

    for (Int_t reg=0;reg<16;reg++){
	fTriggerOutputRegionalData[reg]=0;
	for (Int_t bit=0;bit<4;bit++){
	    fTriggerInputGlobalDataLPt[reg][bit]=0;
	    fTriggerInputGlobalDataHPt[reg][bit]=0;
	}
    }

    for (Int_t bit=0;bit<6;bit++){
	fgotmp[bit]=0;
	fTriggerOutputGlobalData[bit]=0;
	fTriggerOutputGlobalRecFromGlobalInput[bit]=0;
    }

    for (Int_t loc=0;loc<243;loc++){
	for (Int_t bit=0;bit<16;bit++){
	    fTriggerPatternX1[loc][bit]=0;
	    fTriggerPatternX2[loc][bit]=0;
	    fTriggerPatternX3[loc][bit]=0;
	    fTriggerPatternX4[loc][bit]=0;

	    fTriggerPatternY1[loc][bit]=0;
	    fTriggerPatternY2[loc][bit]=0;
	    fTriggerPatternY3[loc][bit]=0;
	    fTriggerPatternY4[loc][bit]=0;
	}
    }

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
    
    AliMUONTriggerStoreV1 triggerStore;
    triggerStore.Create();
    triggerStore.Clear();


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
	for (Int_t bit=1; bit<7; bit++){
	  fTriggerOutputGlobalData[bit-1]=Int_t(((darcHeader->GetGlobalOutput())>>bit)&1);
	  if ( fillScalerHistos && !fTriggerOutputGlobalData[bit-1] )
	    nDeadGlobal++;
	}
	for (Int_t Bit=0; Bit<32; Bit++){
	  fTriggerInputGlobalDataLPt[Bit/4][Bit%4]=((darcHeader->GetGlobalInput(0)>>Bit)&1);
	  fTriggerInputGlobalDataLPt[Bit/4+8][Bit%4]=((darcHeader->GetGlobalInput(1)>>Bit)&1);
	  fTriggerInputGlobalDataHPt[Bit/4][Bit%4]=((darcHeader->GetGlobalInput(2)>>Bit)&1);
	  fTriggerInputGlobalDataHPt[Bit/4+8][Bit%4]=((darcHeader->GetGlobalInput(3)>>Bit)&1);
	}

	globaltemp[0]=darcHeader->GetGlobalInput(0);
	globaltemp[1]=darcHeader->GetGlobalInput(1);
	globaltemp[2]=darcHeader->GetGlobalInput(2);
	globaltemp[3]=darcHeader->GetGlobalInput(3);
      }

      Int_t nReg = rawStreamTrig.GetRegionalHeaderCount();

      for(Int_t iReg = 0; iReg < nReg ;iReg++)
      {   //reg loop

	  Int_t regId=rawStreamTrig.GetDDL()*8+iReg;

	// crate info  
	  AliMpTriggerCrate* crate = AliMpDDLStore::Instance()->GetTriggerCrate(rawStreamTrig.GetDDL(), iReg);

	  regHeader =  rawStreamTrig.GetRegionalHeader(iReg);

	//Get regional outputs -> not checked, hardware read-out doesn't work
	fTriggerOutputRegionalData[regId]=Int_t(regHeader->GetOutput());
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

	  countNotifiedBoards++;  

	  TArrayS xyPattern[2];	  
	  localStruct->GetXPattern(xyPattern[0]);
	  localStruct->GetYPattern(xyPattern[1]);
	  fDigitMaker->TriggerDigits(loCircuit, xyPattern, digitStore);
	  if ( fillScalerHistos ) // Compute total number of strips
	    fDigitMaker->TriggerDigits(loCircuit, xyPatternAll, digitStoreAll);

	  Int_t cathode = localStruct->GetComptXY()%2;

	  //Get electronic Decisions from data

	  //Get regional inputs -> not checked, hardware read-out doesn't work
	  fTriggerInputRegionalDataLPt[0][loCircuit]=Int_t(((regHeader->GetInput(0))>>(2*iLocal))&1);
	  fTriggerInputRegionalDataLPt[1][loCircuit]=Int_t(((regHeader->GetInput(1))>>((2*iLocal)+1))&1);

	  //Get local in/outputs
	  if (Int_t(localStruct->GetDec())!=0){
	      fTriggerOutputLocalDataTriggerDec[loCircuit]++;
	      ((TH1F*)GetRawsData(kTriggeredBoards))->Fill(loCircuit);
	  }
	  else if ( fillScalerHistos ){
	    nDeadLocal++;
	  }
	  
	  fTriggerOutputLocalDataLPtDec[0][loCircuit]=((localStruct->GetLpt())&1);
	  fTriggerOutputLocalDataLPtDec[1][loCircuit]=((localStruct->GetLpt()>>1)&1);
	  fTriggerOutputLocalDataHPtDec[0][loCircuit]=((localStruct->GetHpt())&1);
	  fTriggerOutputLocalDataHPtDec[1][loCircuit]=((localStruct->GetHpt()>>1)&1);
	  fTriggerOutputLocalDataXPos[loCircuit]=Int_t(localStruct->GetXPos());
	  fTriggerOutputLocalDataYPos[loCircuit]=Int_t(localStruct->GetYPos());
	  fTriggerOutputLocalDataDev[loCircuit]=Int_t((localStruct->GetXDev())*(pow(-1.0,(localStruct->GetSXDev()))));
	  fTriggerOutputLocalDataTrigY[loCircuit]=Int_t(localStruct->GetTrigY());
	  
	  UShort_t x1  = (Int_t)localStruct->GetX1();
	  UShort_t x2  = (Int_t)localStruct->GetX2();
	  UShort_t x3  = (Int_t)localStruct->GetX3();
	  UShort_t x4  = (Int_t)localStruct->GetX4();

	  UShort_t y1  = (Int_t)localStruct->GetY1();
	  UShort_t y2  = (Int_t)localStruct->GetY2();
	  UShort_t y3  = (Int_t)localStruct->GetY3();
	  UShort_t y4  = (Int_t)localStruct->GetY4();

	  // loop over strips
	  for (Int_t ibitxy = 0; ibitxy < 16; ++ibitxy) {

	      fTriggerPatternX1[loCircuit][ibitxy]=Int_t((x1>>ibitxy)&1);
	      fTriggerPatternX2[loCircuit][ibitxy]=Int_t((x2>>ibitxy)&1);
	      fTriggerPatternX3[loCircuit][ibitxy]=Int_t((x3>>ibitxy)&1);
	      fTriggerPatternX4[loCircuit][ibitxy]=Int_t((x4>>ibitxy)&1);
	      
	      fTriggerPatternY1[loCircuit][ibitxy]=Int_t((y1>>ibitxy)&1);
	      fTriggerPatternY2[loCircuit][ibitxy]=Int_t((y2>>ibitxy)&1);
	      fTriggerPatternY3[loCircuit][ibitxy]=Int_t((y3>>ibitxy)&1);
	      fTriggerPatternY4[loCircuit][ibitxy]=Int_t((y4>>ibitxy)&1);
	      
	      if ( fillScalerHistos ) {
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
	      } // scaler event
	  } // loop on strips
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

  fTriggerProcessor->Digits2Trigger(digitStore,triggerStore);

  TIter next(triggerStore.CreateLocalIterator());
  AliMUONLocalTrigger *localTrigger;

  while ( ( localTrigger = static_cast<AliMUONLocalTrigger*>(next()) ) )
  {
    
      //... extract information
      loCircuit = localTrigger->LoCircuit();

      AliMpLocalBoard* localBoardMp = AliMpDDLStore::Instance()->GetLocalBoard(loCircuit);  // get local board objectfor switch value
      if (localTrigger->GetLoDecision() != 0){
	  fTriggerOutputLocalRecTriggerDec[loCircuit]++;
      }
      
      fTriggerOutputLocalRecLPtDec[0][loCircuit]=Int_t(localTrigger->LoLpt() & 1);
      fTriggerOutputLocalRecLPtDec[1][loCircuit]=Int_t((localTrigger->LoLpt()>>1) & 1);
      fTriggerOutputLocalRecHPtDec[0][loCircuit]=Int_t(localTrigger->LoHpt() & 1);
      fTriggerOutputLocalRecHPtDec[1][loCircuit]=Int_t((localTrigger->LoHpt()>>1) & 1);
      fTriggerOutputLocalRecXPos[loCircuit]=localTrigger->LoStripX();
      fTriggerOutputLocalRecYPos[loCircuit]=localTrigger->LoStripY();
      fTriggerOutputLocalRecTrigY[loCircuit]=localTrigger->LoTrigY();
      fTriggerOutputLocalRecDev[loCircuit]=Int_t(localTrigger->LoDev()*(pow(-1.,localTrigger->LoSdev())));

      Bool_t firstFillYCopy=kTRUE;

      for (int bit=0; bit<16; bit++){
	  if (fTriggerPatternY1[loCircuit][bit]!=((localTrigger->GetY1Pattern()>>bit) & 1))
	  {
	      fTriggerErrorLocalYCopy[loCircuit]=kTRUE;
	      if (firstFillYCopy){
		  ((TH1F*)GetRawsData(kTriggerErrorLocalYCopy))->Fill(loCircuit);
		  ((TH1F*)GetRawsData(kTriggerError))->Fill(kAlgoLocalYCopy, 1./192.);
		  firstFillYCopy=kFALSE;
	      }
	  }
	  if (fTriggerPatternY2[loCircuit][bit]!=((localTrigger->GetY2Pattern()>>bit) & 1))
	  {
	      fTriggerErrorLocalYCopy[loCircuit]=kTRUE;
	      if (firstFillYCopy){
		  ((TH1F*)GetRawsData(kTriggerErrorLocalYCopy))->Fill(loCircuit);
		  ((TH1F*)GetRawsData(kTriggerError))->Fill(kAlgoLocalYCopy, 1./192.);
		  firstFillYCopy=kFALSE;
	      }
	  }
	  if (fTriggerPatternY3[loCircuit][bit]!=((localTrigger->GetY3Pattern()>>bit) & 1))
	  {
	      fTriggerErrorLocalYCopy[loCircuit]=kTRUE;
	      if (localBoardMp->GetSwitch(4)) fTriggerErrorLocalYCopy[loCircuit-1]=kTRUE;
	      if (localBoardMp->GetSwitch(3)) fTriggerErrorLocalYCopy[loCircuit+1]=kTRUE;
	      if (firstFillYCopy){
		  ((TH1F*)GetRawsData(kTriggerErrorLocalYCopy))->Fill(loCircuit);
		  ((TH1F*)GetRawsData(kTriggerError))->Fill(kAlgoLocalYCopy, 1./192.);
		  firstFillYCopy=kFALSE;
	      }
	  }
	  if (fTriggerPatternY4[loCircuit][bit]!=((localTrigger->GetY4Pattern()>>bit) & 1))
	  {
	      fTriggerErrorLocalYCopy[loCircuit]=kTRUE;
	      if (localBoardMp->GetSwitch(4)) fTriggerErrorLocalYCopy[loCircuit-1]=kTRUE;
	      if (localBoardMp->GetSwitch(3)) fTriggerErrorLocalYCopy[loCircuit+1]=kTRUE;
	      if (firstFillYCopy){
		  ((TH1F*)GetRawsData(kTriggerErrorLocalYCopy))->Fill(loCircuit);
		  ((TH1F*)GetRawsData(kTriggerError))->Fill(kAlgoLocalYCopy, 1./192.);
		  firstFillYCopy=kFALSE;
	      }
	  }
      }
  }

    //Reconstruct Global decision from Global inputs
    for (Int_t bit=0; bit<4; bit++){
	for (Int_t i=0; i<32; i=i+4){
	    fgitmp[bit]+=UInt_t(((globaltemp[bit]>>i)&1)*pow(2.0,i+1));
	    fgitmp[bit]+=UInt_t(((globaltemp[bit]>>(i+1))&1)*pow(2.0,i));
	    fgitmp[bit]+=UInt_t(((globaltemp[bit]>>(i+2))&1)*pow(2.0,i+2));
	    fgitmp[bit]+=UInt_t(((globaltemp[bit]>>(i+3))&1)*pow(2.0,i+3));
	    }
    }
    RawTriggerInGlobal2OutGlobal();
    for (Int_t bit=0; bit<6; bit++){
	fTriggerOutputGlobalRecFromGlobalInput[bit]=fgotmp[bit];
    }

    // Compare data and reconstructed decisions and fill histos
    RawTriggerMatchOutLocal();
    RawTriggerMatchOutLocalInRegional(); // Not tested, hardware read-out doesn't work
    RawTriggerMatchOutGlobalFromInGlobal();
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
  // Do nothing in case of calibration event
  if ( GetRecoParam()->GetEventSpecie() == AliRecoParam::kCalib ) return;
	
  GetRecPointsData(kNAnalyzedEvents)->Fill(1.);
}

//____________________________________________________________________________
void AliMUONTriggerQADataMakerRec::MakeESDs(AliESDEvent* /*esd*/)
{  
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
void AliMUONTriggerQADataMakerRec::RawTriggerInGlobal2OutGlobal()
{
  //
  /// Reconstruct Global Trigger decision using Global Inputs
  //

    AliMUONGlobalCrateConfig* globalConfig = fCalibrationData->GlobalTriggerCrateConfig();

    AliMUONGlobalTriggerBoard globalTriggerBoard;
    globalTriggerBoard.Reset();
    for (Int_t i = 0; i < 4; i++) {
	globalTriggerBoard.Mask(i,globalConfig->GetGlobalMask(i));
    }


    UShort_t regional[16];

    for (Int_t iReg = 0; iReg < 16; iReg++) {
      regional[iReg] = 0;
      if (iReg < 8) {    // right
	// Lpt
	regional[iReg] |=  (fgitmp[0] >> (4*iReg))     & 0xF;
	// Hpt
	regional[iReg] |= ((fgitmp[2] >> (4*iReg))     & 0xF) << 4;
      } else {           // left
	// Lpt
	regional[iReg] |=  (fgitmp[1] >> (4*(iReg-8))) & 0xF;
	// Hpt
	regional[iReg] |= ((fgitmp[3] >> (4*(iReg-8))) & 0xF) << 4;
      }
    }
    globalTriggerBoard.SetRegionalResponse(regional);
    globalTriggerBoard.Response();

    for (Int_t bit=1; bit<7; bit++){
	fgotmp[bit-1]=Int_t((globalTriggerBoard.GetResponse())>>bit&1);
    }
}

//____________________________________________________________________________ 
void AliMUONTriggerQADataMakerRec::RawTriggerMatchOutLocal()
{
  //
  /// Match data and reconstructed Local Trigger decision
  //

    Bool_t firstFillXPosDev=kTRUE;
    Bool_t firstFillYPosTrigY=kTRUE;
    Bool_t firstFillLUT=kTRUE;

    for (int localId=1;localId<235;localId++){
	if(fTriggerOutputLocalDataTriggerDec[localId]!=fTriggerOutputLocalRecTriggerDec[localId]){
	    ((TH1F*)GetRawsData(kTriggerErrorLocalTriggerDec))->Fill(localId);
	}
	if(fTriggerOutputLocalDataTrigY[localId]!=fTriggerOutputLocalRecTrigY[localId]){
	    if(fTriggerErrorLocalYCopy[localId]) continue;
	    ((TH1F*)GetRawsData(kTriggerErrorLocalTrigY))->Fill(localId);
	     if (firstFillYPosTrigY){
		 ((TH1F*)GetRawsData(kTriggerError))->Fill(kAlgoLocalY);
		 firstFillYPosTrigY=kFALSE;
	     }
	}

	if(fTriggerOutputLocalDataYPos[localId]!=fTriggerOutputLocalRecYPos[localId]){
	    if(fTriggerErrorLocalYCopy[localId]) continue;
	    ((TH1F*)GetRawsData(kTriggerErrorLocalYPos))->Fill(localId);
	    if (firstFillYPosTrigY){
		 ((TH1F*)GetRawsData(kTriggerError))->Fill(kAlgoLocalY);
		 firstFillYPosTrigY=kFALSE;
	     }
	}
	if(fTriggerOutputLocalDataXPos[localId]!=fTriggerOutputLocalRecXPos[localId]){
	    ((TH1F*)GetRawsData(kTriggerErrorLocalXPos))->Fill(localId);
	     if (firstFillXPosDev){
		 ((TH1F*)GetRawsData(kTriggerError))->Fill(kAlgoLocalX);
		 firstFillXPosDev=kFALSE;
	     }
	}
	if(fTriggerOutputLocalDataDev[localId]!=fTriggerOutputLocalRecDev[localId]){
	    ((TH1F*)GetRawsData(kTriggerErrorLocalDev))->Fill(localId);
	     if (firstFillXPosDev){
		 ((TH1F*)GetRawsData(kTriggerError))->Fill(kAlgoLocalX);
		 firstFillXPosDev=kFALSE;
	     }
	}
	if(fTriggerOutputLocalDataLPtDec[0][localId]!=fTriggerOutputLocalRecLPtDec[0][localId]){
	    ((TH1F*)GetRawsData(kTriggerErrorLocalLPtLSB))->Fill(localId);
	     if (firstFillLUT){
		 ((TH1F*)GetRawsData(kTriggerError))->Fill(kAlgoLocalLUT);
		 firstFillLUT=kFALSE;
	     }
	}
	if(fTriggerOutputLocalDataLPtDec[1][localId]!=fTriggerOutputLocalRecLPtDec[1][localId]){
	    ((TH1F*)GetRawsData(kTriggerErrorLocalLPtMSB))->Fill(localId);
	     if (firstFillLUT){
		 ((TH1F*)GetRawsData(kTriggerError))->Fill(kAlgoLocalLUT);
		 firstFillLUT=kFALSE;
	     }
	}
	if(fTriggerOutputLocalDataHPtDec[0][localId]!=fTriggerOutputLocalRecHPtDec[0][localId]){
	    ((TH1F*)GetRawsData(kTriggerErrorLocalHPtLSB))->Fill(localId);
	     if (firstFillLUT){
		 ((TH1F*)GetRawsData(kTriggerError))->Fill(kAlgoLocalLUT);
		 firstFillLUT=kFALSE;
	     }
	}
	if(fTriggerOutputLocalDataHPtDec[1][localId]!=fTriggerOutputLocalRecHPtDec[1][localId]){
	    ((TH1F*)GetRawsData(kTriggerErrorLocalHPtMSB))->Fill(localId);
	     if (firstFillLUT){
		 ((TH1F*)GetRawsData(kTriggerError))->Fill(kAlgoLocalLUT);
		 firstFillLUT=kFALSE;
	     }
	}
    } // loop over Local Boards
}

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

//____________________________________________________________________________ 
void AliMUONTriggerQADataMakerRec::RawTriggerMatchOutGlobalFromInGlobal()
{
  //
  /// Match data and reconstructed Global Trigger decision for a reconstruction from Global inputs
  //

  Bool_t firstFill=kTRUE;

  for (int bit=0;bit<6;bit++){
    if(fTriggerOutputGlobalData[bit]!=0){
      ((TH1F*)GetRawsData(kTriggerGlobalOutput))->Fill(5-bit);
    }
    if(fTriggerOutputGlobalData[bit]!=fTriggerOutputGlobalRecFromGlobalInput[bit]){
      ((TH1F*)GetRawsData(kTriggerErrorOutGlobalFromInGlobal))->Fill(5-bit);
      if (firstFill){
	((TH1F*)GetRawsData(kTriggerError))->Fill(kAlgoGlobalFromGlobal);
	firstFill=kFALSE;
      }
    }
  }
}
