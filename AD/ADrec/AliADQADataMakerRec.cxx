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


//  Produces the data needed to calculate the quality assurance 
//  All data must be mergeable objects
//  Handles ESDs and Raws
//  Histos defined will be used for Raw Data control and monitoring

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TF1.h> 
#include <TH1F.h> 
#include <TH1I.h> 
#include <TH2I.h> 
#include <TH2F.h> 
#include <TGraph.h> 
#include <TParameter.h>
#include <TTimeStamp.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliADQADataMakerRec.h"
#include "AliQAChecker.h"
#include "AliRawReader.h"
#include "AliADRawStream.h"
#include "AliADdigit.h"
#include "AliADConst.h"
#include "AliADReconstructor.h"
#include "AliADTrending.h"
#include "AliADCalibData.h"
#include "AliADRecoParam.h"
#include "AliCTPTimeParams.h"
#include "event.h"

ClassImp(AliADQADataMakerRec)
           
//____________________________________________________________________________ 
AliADQADataMakerRec::AliADQADataMakerRec() : 
AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kAD), "AD Quality Assurance Data Maker"),
  fCalibData(0x0),
  fRecoParam(0x0),
  fTrendingUpdateTime(0), 
  fCycleStartTime(0), 
  fCycleStopTime(0),
  fTimeSlewing(0)
    
{
  // Constructor
   
  AliDebug(AliQAv1::GetQADebugLevel(), "Construct AD QA Object");

  for(Int_t i=0; i<16; i++){  
    fEven[i] = 0;   
    fOdd[i]  = 0;
  }
  
  for(Int_t i=0; i<32; i++){  
    fADCmean[i] = 0.0;   }	
}

//____________________________________________________________________________ 
AliADQADataMakerRec::AliADQADataMakerRec(const AliADQADataMakerRec& qadm) :
  AliQADataMakerRec(),
  fCalibData(0x0),
  fRecoParam(0x0),
  fTrendingUpdateTime(0), 
  fCycleStartTime(0), 
  fCycleStopTime(0),
  fTimeSlewing(0)
  
{
  // Copy constructor 
  
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliADQADataMakerRec& AliADQADataMakerRec::operator = (const AliADQADataMakerRec& qadm )
{
  // Equal operator
  
  this->~AliADQADataMakerRec();
  new(this) AliADQADataMakerRec(qadm);
  return *this;
}

//____________________________________________________________________________
AliADCalibData* AliADQADataMakerRec::GetCalibData() const

{
  AliCDBManager *man = AliCDBManager::Instance();

  AliCDBEntry *entry=0;

  entry = man->Get("AD/Calib/Data",fRun);
  if(!entry){
    AliWarning("Load of calibration data from default storage failed!");
    AliWarning("Calibration data will be loaded from local storage ($ALICE_ROOT)");
	
    man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    entry = man->Get("AD/Calib/Data",fRun);
  }
  // Retrieval of data in directory AD/Calib/Data:

  AliADCalibData *calibdata = 0;

  if (entry) calibdata = (AliADCalibData*) entry->GetObject();
  if (!calibdata)  AliFatal("No calibration data from calibration database !");

  return calibdata;
}
//____________________________________________________________________________ 
void AliADQADataMakerRec::StartOfDetectorCycle()
{
  // Detector specific actions at start of cycle
  
  // Reset of the histogram used - to have the trend versus time -
 
  fCalibData = GetCalibData();
  fRecoParam = (AliADRecoParam*)GetRecoParam();
 
  AliCDBEntry *entry = AliCDBManager::Instance()->Get("GRP/CTP/CTPtiming");
  if (!entry) AliFatal("CTP timing parameters are not found in OCDB !");
  AliCTPTimeParams *ctpParams = (AliCTPTimeParams*)entry->GetObject();
  Float_t l1Delay = (Float_t)ctpParams->GetDelayL1L0()*25.0;

  AliCDBEntry *entry1 = AliCDBManager::Instance()->Get("GRP/CTP/TimeAlign");
  if (!entry1) AliFatal("CTP time-alignment is not found in OCDB !");
  AliCTPTimeParams *ctpTimeAlign = (AliCTPTimeParams*)entry1->GetObject();
  l1Delay += ((Float_t)ctpTimeAlign->GetDelayL1L0()*25.0);
  /*/
  AliCDBEntry *entry2 = AliCDBManager::Instance()->Get("AD/Calib/TimeDelays");
  if (!entry2) AliFatal("AD time delays are not found in OCDB !");
  TH1F *delays = (TH1F*)entry2->GetObject();
  /*/
  AliCDBEntry *entry3 = AliCDBManager::Instance()->Get("AD/Calib/TimeSlewing");
  if (!entry3) AliFatal("AD time slewing function is not found in OCDB !");
  fTimeSlewing = (TF1*)entry3->GetObject();
 

  for(Int_t i = 0 ; i < 16; ++i) {
    //Int_t board = AliADCalibData::GetBoardNumber(i);
    fTimeOffset[i] = (
		      //	((Float_t)fCalibData->GetTriggerCountOffset(board) -
		      //	(Float_t)fCalibData->GetRollOver(board))*25.0 +
		      //     fCalibData->GetTimeOffset(i) -
		      //     l1Delay+
		      //delays->GetBinContent(i+1)//+
		      //      kADOffset
		      0
		      );
    //		      AliInfo(Form(" fTimeOffset[%d] = %f  kADoffset %f",i,fTimeOffset[i],kADOffset));
  }

 
 
  
 	
  TTimeStamp currentTime;
  fCycleStartTime = currentTime.GetSec();
 
}
//____________________________________________________________________________ 
void AliADQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  // Detector specific actions at end of cycle
  // Does the QA checking
  ResetEventTrigClasses();
  
  if(task == AliQAv1::kRAWS){
    TTimeStamp currentTime;
    fCycleStopTime = currentTime.GetSec();
    
    Int_t nCorrelation = 0;
    Int_t nPair = 1;
    for(Int_t i=0;i<8;i++){
  	for(Int_t j=7;j>i;j--){
  		if( (j - i) == 4){ 
			Float_t Mean = ((TH1F*)GetRawsData(kNTimeDiffADC + nCorrelation))->GetMean();
			Float_t RMS = ((TH1F*)GetRawsData(kNTimeDiffADC + nCorrelation))->GetRMS();
			SetRawsDataBinContent(kPairTimeDiffMean,nPair,Mean);
			SetRawsDataBinContent(kPairTimeDiffRMS,nPair,RMS);
			nPair++;
			}
		 nCorrelation++;
		}
	}
    nCorrelation = 0;
    for(Int_t i=0;i<8;i++){
  	for(Int_t j=7;j>i;j--){
  		if( (j - i) == 4){ 
			Float_t Mean = ((TH1F*)GetRawsData(kNTimeDiffADA + nCorrelation))->GetMean();
			Float_t RMS = ((TH1F*)GetRawsData(kNTimeDiffADA + nCorrelation))->GetRMS();
			SetRawsDataBinContent(kPairTimeDiffMean,nPair,Mean);
			SetRawsDataBinContent(kPairTimeDiffRMS,nPair,RMS);
			nPair++;
			}
		 nCorrelation++;
		}
	}		
    
  }

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    if (! IsValidEventSpecie(specie, list)) continue ;
    SetEventSpecie(AliRecoParam::ConvertIndex(specie));
    if(task == AliQAv1::kRAWS) {
    } else if (task == AliQAv1::kESDS) {
    }
  }
  AliQAChecker::Instance()->Run(AliQAv1::kAD, task, list) ;
}

//____________________________________________________________________________ 
void AliADQADataMakerRec::InitESDs()
{
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  TH1I * h0 = new TH1I("H1I_Cell_Multiplicity_ADA", "Cell Multiplicity in ADA;Multiplicity (Nb of Cell);Counts", 35, 0, 35) ;  
  Add2ESDsList(h0, kCellMultiADA, !expert, image)  ;  
                                                                                                        
  TH1I * h1 = new TH1I("H1I_Cell_Multiplicity_ADC", "Cell Multiplicity in AD;Multiplicity (Nb of Cell);Counts", 35, 0, 35) ;  
  Add2ESDsList(h1, kCellMultiADC, !expert, image)  ;  
  
  TH1F * h2 = new TH1F("H1D_BBFlag_Counters", "BB Flag Counters;Channel;Counts",16, 0, 16) ;  
  Add2ESDsList(h2, kBBFlag, !expert, image)  ;  
  
  TH1F * h3 = new TH1F("H1D_BGFlag_Counters", "BG Flag Counters;Channel;Counts",16, 0, 16) ;  
  Add2ESDsList(h3, kBGFlag, !expert, image)  ;  
  
  TH2F * h4 = new TH2F("H2D_Charge_Channel", "ADC Charge per channel;Channel;Charge (ADC counts)",16, 0, 16, 1024, 0, 1024) ;  
  Add2ESDsList(h4, kChargeChannel, !expert, image)  ;  
  
  TH2F * h5 = new TH2F("H2D_Time_Channel", "Time per channel;Channel;Time (ns)",16, 0, 16, 400, -100, 100) ;  
  Add2ESDsList(h5, kTimeChannel, !expert, image)  ;  
  
  TH1F * h6 = new TH1F("H1D_ADA_Time", "Mean ADA Time;Time (ns);Counts",1000, -100., 100.);
  Add2ESDsList(h6,kESDADATime, !expert, image); 
  
  TH1F * h7 = new TH1F("H1D_ADC_Time", "Mean ADC Time;Time (ns);Counts",1000, -100., 100.);
  Add2ESDsList(h7,kESDADCTime, !expert, image); 
  
  TH1F * h8 = new TH1F("H1D_Diff_Time", "Diff Time ADA - ADC;Diff Time ADA - ADC (ns);Counts",1000, -200., 200.);
  Add2ESDsList(h8,kESDDiffTime, !expert, image); 
  //
  ClonePerTrigClass(AliQAv1::kESDS); // this should be the last line	
}

//____________________________________________________________________________ 
void AliADQADataMakerRec::InitDigits()
{
// create Digits histograms in Digits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  // create Digits histograms in Digits subdir
  TH1I * h0 = new TH1I("hDigitMultiplicity", "Digits multiplicity distribution in AD;# of Digits;Entries", 17,-0.5,16.5) ; 
  h0->Sumw2() ;
  Add2DigitsList(h0, 0, !expert, image) ;
     
  TH2D * h1 = new TH2D("hDigitLeadingTimePerPM", "Leading time distribution per PM in AD;PM number;Leading Time [ns]",16,0,16, 1000, 200, 300); 
  h1->Sumw2() ;
  Add2DigitsList(h1, 1, !expert, image) ; 
  
  TH2D * h2 = new TH2D("hDigitTimeWidthPerPM", "Time width distribution per PM in AD;PM number;Time width [ns]",16,0,16, 1000, 0, 100); 
  h2->Sumw2() ;
  Add2DigitsList(h2, 2, !expert, image) ;
  
  TH2I * h3 = new TH2I("hDigitChargePerClockPerPM", "Charge array per PM in AD;PM number; Clock",16,0,16,21, -10.5, 10.5);
  h3->Sumw2();
  Add2DigitsList(h3, 3, !expert, image) ;
  
  TH1I * h4 = new TH1I("hDigitBBflagsAD","Number of BB flags in AD; # of BB flags; Entries",17,-0.5,16.5);
  h4->Sumw2();
  Add2DigitsList(h4, 4, !expert, image) ;
  
  TH1I * h5 = new TH1I("hDigitBBflagsADA","Number of BB flags in ADA; # of BB flags; Entries",9,-0.5,8.5);
  h5->Sumw2();
  Add2DigitsList(h5, 5, !expert, image) ;
  
  TH1I * h6 = new TH1I("hDigitBBflagsADC","Number of BB flags in ADC; # of BB flags; Entries",9,-0.5,8.5);
  h6->Sumw2();
  Add2DigitsList(h6, 6, !expert, image) ;
  
  TH2D * h7 = new TH2D("hDigitTotalChargePerPM", "Total Charge per PM in AD;PM number; Charge [ADC counts]",16,0,16,10000,0,10000);
  h7->Sumw2();
  Add2DigitsList(h7, 7, !expert, image) ;
  
  TH2I * h8 = new TH2I("hDigitMaxChargeClockPerPM", "Clock with maximum charge per PM in AD;PM number; Clock ",16,0,16,21, -10.5, 10.5);
  h8->Sumw2();
  Add2DigitsList(h8, 8, !expert, image) ;
   
  //
  ClonePerTrigClass(AliQAv1::kDIGITS); // this should be the last line
}

//____________________________________________________________________________
void AliADQADataMakerRec::MakeDigits()
{
 // makes data from Digits

  FillDigitsData(0,fDigitsArray->GetEntriesFast()) ; 
  TIter next(fDigitsArray) ; 
    AliADdigit *ADDigit ; 
    Int_t nBBflagsADA = 0;
    Int_t nBBflagsADC = 0;
    
    while ( (ADDigit = dynamic_cast<AliADdigit *>(next())) ) {
         Int_t totCharge = 0;
         Int_t   PMNumber  = ADDigit->PMNumber();
	 if(PMNumber<8 && ADDigit->GetBBflag()) nBBflagsADC++;
	 if(PMNumber>7 && ADDigit->GetBBflag()) nBBflagsADA++;
	 
	 Short_t adc[21];
	 for(Int_t iClock=0; iClock<21; iClock++) { 
	 adc[iClock]= ADDigit->ChargeADC(iClock);
	 FillDigitsData(3, PMNumber,(float)iClock-10,(float)adc[iClock]);
	 totCharge += adc[iClock];
	 }
	    
         FillDigitsData(1,PMNumber,ADDigit->Time()); 
	 FillDigitsData(2,PMNumber,ADDigit->Width());
	 FillDigitsData(7,PMNumber,totCharge);
	 FillDigitsData(8,PMNumber,TMath::LocMax(21,adc)-10); 
	 
    }
    FillDigitsData(4,nBBflagsADA+nBBflagsADC);
    FillDigitsData(5,nBBflagsADA);
    FillDigitsData(6,nBBflagsADC);  
}

//____________________________________________________________________________
void AliADQADataMakerRec::MakeDigits(TTree* digitTree)
{
  // makes data from Digit Tree
	
  if (fDigitsArray)
    fDigitsArray->Clear() ; 
  else 
    fDigitsArray = new TClonesArray("AliADdigit", 1000) ; 

    TBranch * branch = digitTree->GetBranch("ADDigit") ;
    if ( ! branch ) {
         AliWarning("AD branch in Digit Tree not found") ; 
    } else {
         branch->SetAddress(&fDigitsArray) ;
         branch->GetEntry(0) ; 
         MakeDigits() ; 
    }  
    //
    IncEvCountCycleDigits();
    IncEvCountTotalDigits();
    //    
}


//____________________________________________________________________________
void AliADQADataMakerRec::MakeESDs(AliESDEvent* esd)
{
// Creates QA data from ESDs
  
  UInt_t eventType = esd->GetEventType();

  switch (eventType){
  case PHYSICS_EVENT:
    AliESDAD *esdAD=esd->GetADData();
   
    if (!esdAD) break;
		  
    FillESDsData(kCellMultiADA,esdAD->GetNbPMADA());
    FillESDsData(kCellMultiADC,esdAD->GetNbPMADC());   
    	
    for(Int_t i=0;i<16;i++) {
      FillESDsData(kChargeChannel,(Float_t) i,(Float_t) esdAD->GetAdc(i));
      if (i < 8) {
	if(esdAD->BBTriggerADC(i)) FillESDsData(kBBFlag,(Float_t) i);
	if(esdAD->BGTriggerADC(i)) FillESDsData(kBGFlag,(Float_t) i);
      }
      else {
	if(esdAD->BBTriggerADA(i-8)) FillESDsData(kBBFlag,(Float_t) i);  
	if(esdAD->BGTriggerADA(i-8)) FillESDsData(kBGFlag,(Float_t) i);
      }		  	
      Float_t time = (Float_t) esdAD->GetTime(i);
      FillESDsData(kTimeChannel,(Float_t) i,time);
    }
				
    Float_t timeADA = esdAD->GetADATime();
    Float_t timeADC = esdAD->GetADCTime();
    Float_t diffTime;

    if(timeADA<-1024.+1.e-6 || timeADC<-1024.+1.e-6) diffTime = -1024.;
    else diffTime = timeADA - timeADC;

    FillESDsData(kESDADATime,timeADA);
    FillESDsData(kESDADCTime,timeADC);
    FillESDsData(kESDDiffTime,diffTime);
		
    break;
  }  
  //
  IncEvCountCycleESDs();
  IncEvCountTotalESDs();  
  // 
}

//____________________________________________________________________________ 
void AliADQADataMakerRec::InitRaws()
{
  // Creates RAW histograms in Raws subdir
  if(!fRecoParam) fRecoParam = (AliADRecoParam*)GetRecoParam();
 
  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  const Int_t kNintegrator  =    2;
 
  const Int_t kNTdcTimeBins  = fRecoParam->GetNTdcTimeBins();
  const Float_t kTdcTimeMin    =  fRecoParam->GetTdcTimeMin();
  const Float_t kTdcTimeMax    = fRecoParam->GetTdcTimeMax();
  const Int_t kNTdcWidthBins =  fRecoParam->GetNTdcWidthBins();
  const Float_t kTdcWidthMin   =    fRecoParam->GetTdcWidthMin();
  const Float_t kTdcWidthMax   =  fRecoParam->GetTdcWidthMax();
  const Int_t kNChargeChannelBins   =  fRecoParam->GetNChargeChannelBins();
  const Int_t kNChargeSideBins   = fRecoParam->GetNChargeSideBins();
  const Int_t kNChargeCorrBins   = fRecoParam->GetNChargeCorrBins();
   
  const Float_t kChargeChannelMin     =    1;
  const Float_t kChargeChannelMax     = 1+kNChargeChannelBins;
  const Float_t kChargeSideMin     =    1;
  const Float_t kChargeSideMax     = 1+kNChargeSideBins;
  const Float_t kChargeCorrMin     =    0;
  const Float_t kChargeCorrMax     = kNChargeCorrBins;
  
  const Int_t kNTimeCorrBins = 614; 
  const Float_t kTimeCorrMin = 70.019531;
  const Float_t kTimeCorrMax =  129.980469; 
   
  const Int_t kNTimeDiffBins = 154; 
  const Float_t kTimeDiffMin = -15.039062;
  const Float_t kTimeDiffMax =  15.039062;
  
  const Int_t kNChannelBins  =   16;
  const Float_t kChannelMin    =    -0.5;
  const Float_t kChannelMax    =   15.5;
  
  const Int_t kNPedestalBins =  40;
  const Float_t kPedestalMin   =    0;
  const Float_t kPedestalMax   =  40; 
  
  const Int_t kNPairBins  =   8;
  const Float_t kPairMin    =    -0.5;
  const Float_t kPairMax    =   7.5;
  
  TH2I * h2i;
  TH2F * h2d;
  TH1I * h1i;
  TH1F * h1d;

  int iHisto =0;

  // Creation of Cell Multiplicity Histograms
  h1i = new TH1I("H1I_Multiplicity_ADA", "Number of channels with charge signal and time ADA;# of Channels;Entries", 9, -0.5, 8.5) ;  
  Add2RawsList(h1i,kMultiADA, expert, !image, !saveCorr);   iHisto++;
  h1i = new TH1I("H1I_Multiplicity_ADC", "Number of channels with charge signal and time ADC;# of Channels;Entries", 9, -0.5, 8.5) ;  
  Add2RawsList(h1i,kMultiADC, expert, !image, !saveCorr);   iHisto++;
 
  // Creation of Total Charge Histograms
  h1d = new TH1F("H1D_Charge_ADA", "Total Charge in ADA;Charge [ADC counts];Counts", kNChargeSideBins, kChargeSideMin, kChargeSideMax) ;  
  Add2RawsList(h1d,kChargeADA, !expert, image, saveCorr);   iHisto++;
  h1d = new TH1F("H1D_Charge_ADC", "Total Charge in ADC;Charge [ADC counts];Counts", kNChargeSideBins, kChargeSideMin, kChargeSideMax) ;  
  Add2RawsList(h1d,kChargeADC, !expert, image, saveCorr);   iHisto++;
  h1d = new TH1F("H1D_Charge_AD", "Total Charge in AD;Charge [ADC counts];Counts", 2*kNChargeSideBins, kChargeSideMin, 1+2*kNChargeSideBins) ;  
  Add2RawsList(h1d,kChargeAD, !expert,  !image, !saveCorr);   iHisto++;
   

  // Creation of Charge EoI histogram 
  h2d = new TH2F("H2D_ChargeEoI", "Signal charge per channel(pedestal substracted);Channel Number;Charge [ADC counts]"
		 ,kNChannelBins, kChannelMin, kChannelMax, kNChargeChannelBins, kChargeChannelMin, kChargeChannelMax);
  Add2RawsList(h2d,kChargeEoI, !expert, image, saveCorr); iHisto++;

  for(Int_t iInt=0;iInt<kNintegrator;iInt++){
    // Creation of Pedestal histograms 
    h2i = new TH2I(Form("H2I_Pedestal_Int%d",iInt), Form("Pedestal (Int%d);Channel;Pedestal [ADC counts]",iInt)
		   ,kNChannelBins, kChannelMin, kChannelMax,kNPedestalBins,kPedestalMin ,kPedestalMax );
    Add2RawsList(h2i,(iInt == 0 ? kPedestalInt0 : kPedestalInt1), expert, image, !saveCorr); iHisto++;
    
    h2d = new TH2F(Form("H2D_PedestalDiff_Int%d",iInt), Form("Pedestal difference Online - OCDB (Int%d);Channel; Pedestal Online - OCDB",iInt)
		   ,kNChannelBins, kChannelMin, kChannelMax,81,-10.5,70.5);
    Add2RawsList(h2d,(iInt == 0 ? kPedestalDiffInt0 : kPedestalDiffInt1), !expert, image, !saveCorr); iHisto++;
	

    // Creation of Charge EoI histograms 
    h2i = new TH2I(Form("H2I_ChargeEoI_Int%d",iInt), Form("Maximum charge per clock (Int%d);Channel;Charge [ADC counts]",iInt)
		   ,kNChannelBins, kChannelMin, kChannelMax, 1025, 0, 1025);
    Add2RawsList(h2i,(iInt == 0 ? kChargeEoIInt0 : kChargeEoIInt1), !expert, image, saveCorr); iHisto++;
    
  }	
  
  // Creation of Time histograms 
  h2i = new TH2I("H2I_Width", "HPTDC Width;Channel;Width [ns]",kNChannelBins, kChannelMin, kChannelMax, kNTdcWidthBins, kTdcWidthMin, kTdcWidthMax);
  Add2RawsList(h2i,kWidth, expert, image, saveCorr); iHisto++;
  
  h2i = new TH2I("H2I_Width_BB", "HPTDC Width w/ BB Flag condition;Channel;Width [ns]",kNChannelBins, kChannelMin, kChannelMax, kNTdcWidthBins, kTdcWidthMin, kTdcWidthMax);
  Add2RawsList(h2i,kWidthBB, expert, !image, !saveCorr); iHisto++;

  h2i = new TH2I("H2I_Width_BG", "HPTDC Width w/ BG Flag condition;Channel;Width [ns]",kNChannelBins, kChannelMin, kChannelMax, kNTdcWidthBins, kTdcWidthMin, kTdcWidthMax);
  h2i->SetDrawOption("colz");
  Add2RawsList(h2i,kWidthBG, expert, !image, !saveCorr); iHisto++;

  h2i = new TH2I("H2I_HPTDCTime", "HPTDC Time;Channel;Leading Time [ns]",kNChannelBins, kChannelMin, kChannelMax, kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax);
  Add2RawsList(h2i,kHPTDCTime, !expert, image, saveCorr); iHisto++;
  
  h2i = new TH2I("H2I_HPTDCTime_BB", "HPTDC Time w/ BB Flag condition;Channel;Leading Time [ns]",kNChannelBins, kChannelMin, kChannelMax, kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax);
  Add2RawsList(h2i,kHPTDCTimeBB, !expert, image, !saveCorr); iHisto++;

  h2i = new TH2I("H2I_HPTDCTime_BG", "HPTDC Time w/ BG Flag condition;Channel;Leading Time [ns]",kNChannelBins, kChannelMin, kChannelMax, kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax);
  Add2RawsList(h2i,kHPTDCTimeBG, !expert, image, !saveCorr); iHisto++;
	
  h1d = new TH1F("H1D_ADA_Time", "ADA Time;Time [ns];Counts",kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax);
  Add2RawsList(h1d,kADATime, expert, !image, !saveCorr); iHisto++;
	
  h1d = new TH1F("H1D_ADC_Time", "ADC Time;Time [ns];Counts",kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax);
  Add2RawsList(h1d,kADCTime, expert, !image, !saveCorr); iHisto++;
	
  h1d = new TH1F("H1D_Diff_Time","Diff ADA-ADC Time;Time [ns];Counts",kNTimeDiffBins,kTimeDiffMin,kTimeDiffMax);
  Add2RawsList(h1d,kDiffTime, expert, !image, !saveCorr); iHisto++;

  h2d = new TH2F("H2D_TimeADA_ADC", "Mean Time in ADC versus ADA;Time ADA [ns];Time ADC [ns]", kNTdcTimeBins/8, kTdcTimeMin,kTdcTimeMax,kNTdcTimeBins/8, kTdcTimeMin,kTdcTimeMax) ;  
  Add2RawsList(h2d,kTimeADAADC, expert, !image, !saveCorr);   iHisto++;
  
  h2d = new TH2F("H2D_TimeSlewingOff", "Time Vs Charge (no slewing correction);Leading Time[ns];Charge [ADC counts]", kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax, kNChargeCorrBins, kChargeCorrMin, kChargeCorrMax) ;  
  Add2RawsList(h2d,kTimeSlewingOff, expert, !image, !saveCorr);   iHisto++;
  
  h2d = new TH2F("H2D_TimeSlewingOn", "Time Vs Charge (after slewing correction);Leading Time[ns];Charge [ADC counts]", kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax, kNChargeCorrBins, kChargeCorrMin, kChargeCorrMax) ;  
  Add2RawsList(h2d,kTimeSlewingOn, expert, !image, !saveCorr);   iHisto++;
  
  h2d = new TH2F("H2D_WidthSlewing", "Width Vs Charge ;Time Width [ns];Charge [ADC counts]", kNTdcWidthBins, kTdcWidthMin, kTdcWidthMax, kNChargeCorrBins, kChargeCorrMin, kChargeCorrMax) ;  
  Add2RawsList(h2d,kWidthSlewing, expert, !image, !saveCorr);   iHisto++;
  
  //Creation of pair coincidence histograms
  h1i = new TH1I("H1I_MultiBBCoincidence_ADA", "Number of BB flag coincidences in ADA;# of BB Coincidences;Entries", 5, -0.5, 4.5) ;  
  Add2RawsList(h1i,kNBBCoincADA, !expert, image, saveCorr);   iHisto++;
  h1i = new TH1I("H1I_MultiBBCoincidence_ADC", "Number of BB flag coincidences in ADC;# of BB Coincidences;Entries", 5, -0.5, 4.5) ;  
  Add2RawsList(h1i,kNBBCoincADC, !expert, image, saveCorr);   iHisto++;
  
  h1i = new TH1I("H1I_MultiBGCoincidence_ADA", "Number of BG flag coincidences in ADA;# of BG Coincidences;Entries", 5, -0.5, 4.5) ;  
  Add2RawsList(h1i,kNBGCoincADA, !expert, image, saveCorr);   iHisto++;
  h1i = new TH1I("H1I_MultiBGCoincidence_ADC", "Number of BG flag coincidences in ADC;# of BG Coincidences;Entries", 5, -0.5, 4.5) ;  
  Add2RawsList(h1i,kNBGCoincADC, !expert, image, saveCorr);   iHisto++;
  
  h1d = new TH1F("H1D_Pair_TimeDiffMean","Time difference mean for coincidence pair [ns];Pair number;Time mean [ns]",kNPairBins, kPairMin, kPairMax);
  Add2RawsList(h1d,kPairTimeDiffMean, expert, !image, !saveCorr); iHisto++;
  
  h1d = new TH1F("H1D_Pair_TimeDiffRMS","Time difference RMS for coincidence pair [ns];Pair number;Time RMS [ns]",kNPairBins, kPairMin, kPairMax);
  Add2RawsList(h1d,kPairTimeDiffRMS, expert, !image, !saveCorr); iHisto++;

  //Creation of Clock histograms
  h2d = new TH2F("H2D_BBFlagVsClock", "BB-Flags Versus LHC-Clock;Channel;LHC Clocks",kNChannelBins, kChannelMin, kChannelMax,21, -10.5, 10.5 );
  Add2RawsList(h2d,kBBFlagVsClock, !expert, image, saveCorr); iHisto++;
	
  h2d = new TH2F("H2D_BGFlagVsClock", "BG-Flags Versus LHC-Clock;Channel;LHC Clocks",kNChannelBins, kChannelMin, kChannelMax,21, -10.5, 10.5 );
  Add2RawsList(h2d,kBGFlagVsClock, !expert, image, saveCorr); iHisto++;

  for(Int_t iInt=0;iInt<kNintegrator;iInt++){
  	h2d = new TH2F(Form("H2D_ChargeVsClock_Int%d",iInt), Form("Charge Versus LHC-Clock (Int%d);Channel;LHCClock;Charge [ADC counts]",iInt),kNChannelBins, kChannelMin, kChannelMax,21, -10.5, 10.5 );
  	Add2RawsList(h2d,(iInt == 0 ? kChargeVsClockInt0 : kChargeVsClockInt1 ), expert, image, saveCorr); iHisto++;
	}
  
  h2d = new TH2F("H2D_BBFlagPerChannel", "BB-Flags Versus Channel;Channel;BB Flags Count",kNChannelBins, kChannelMin, kChannelMax,22,-0.5,21.5);
  Add2RawsList(h2d,kBBFlagsPerChannel, expert, !image, saveCorr); iHisto++;
  
  h2d = new TH2F("H2D_BGFlagPerChannel", "BG-Flags Versus Channel;Channel;BG Flags Count",kNChannelBins, kChannelMin, kChannelMax,22,-0.5,21.5);
  Add2RawsList(h2d,kBGFlagsPerChannel, expert, !image, saveCorr); iHisto++;
  
  h1d = new TH1F("H1D_FlagNoTime", "Number of events with BB/BG flag but no time measurement;Channel;Entries",kNChannelBins, kChannelMin, kChannelMax);
  Add2RawsList(h1d,kFlagNoTime, !expert, image, saveCorr); iHisto++;
  
  h1d = new TH1F("H1D_TimeNoFlag", "Number of events with time measurement but no BB/BG flag;Channel;Entries",kNChannelBins, kChannelMin, kChannelMax);
  Add2RawsList(h1d,kTimeNoFlag, !expert, image, saveCorr); iHisto++;
  
  //Correlation histograms
  Int_t nCorrelation = 0;
  for(Int_t i=0;i<8;i++){
  	for(Int_t j=7;j>i;j--){
  		h2d = new TH2F(Form("ChargeCorr/H2D_kNChargeCorrADA_%d_%d",i,j),Form("Charge Correlation ADA module%d - module%d",i,j),kNChargeCorrBins, kChargeCorrMin, kChargeCorrMax,kNChargeCorrBins, kChargeCorrMin, kChargeCorrMax);
		Add2RawsList(h2d,kNChargeCorrADA+nCorrelation, expert, !image, !saveCorr); iHisto++; nCorrelation++;
		}
	}
  nCorrelation = 0;
  for(Int_t i=0;i<8;i++){
  	for(Int_t j=7;j>i;j--){
  		h2d = new TH2F(Form("ChargeCorr/H2D_kNChargeCorrADC_%d_%d",i,j),Form("Charge Correlation ADC module%d - module%d",i,j),kNChargeCorrBins, kChargeCorrMin, kChargeCorrMax,kNChargeCorrBins, kChargeCorrMin, kChargeCorrMax);
		Add2RawsList(h2d,kNChargeCorrADC+nCorrelation, expert, !image, !saveCorr); iHisto++; nCorrelation++;
		}
	}
  nCorrelation = 0;
  for(Int_t i=0;i<8;i++){
  	for(Int_t j=7;j>i;j--){
  		h2d = new TH2F(Form("TimeCorr/H2D_kNTimeCorrADA_%d_%d",i,j),Form("Time Correlation ADA module%d - module%d",i,j),kNTimeCorrBins,kTimeCorrMin,kTimeCorrMax,kNTimeCorrBins,kTimeCorrMin,kTimeCorrMax);
		Add2RawsList(h2d,kNTimeCorrADA+nCorrelation, expert, !image, !saveCorr); iHisto++; nCorrelation++;
		}
	}
  nCorrelation = 0;
  for(Int_t i=0;i<8;i++){
  	for(Int_t j=7;j>i;j--){
  		h2d = new TH2F(Form("TimeCorr/H2D_kNTimeCorrADC_%d_%d",i,j),Form("Time Correlation ADC module%d - module%d",i,j),kNTimeCorrBins,kTimeCorrMin,kTimeCorrMax,kNTimeCorrBins,kTimeCorrMin,kTimeCorrMax);
		Add2RawsList(h2d,kNTimeCorrADC+nCorrelation, expert, !image, !saveCorr); iHisto++; nCorrelation++;
		}
	}
  nCorrelation = 0;
  for(Int_t i=0;i<8;i++){
  	for(Int_t j=7;j>i;j--){
  		h1d = new TH1F(Form("TimeDiff/H1D_kNTimeDiffADA_%d_%d",i,j),Form("Time Difference ADA module%d - module%d",i,j),kNTimeDiffBins,kTimeDiffMin,kTimeDiffMax);
		Add2RawsList(h1d,kNTimeDiffADA+nCorrelation, expert, !image, !saveCorr); iHisto++; nCorrelation++;
		}
	}
  nCorrelation = 0;
  for(Int_t i=0;i<8;i++){
  	for(Int_t j=7;j>i;j--){
  		h1d = new TH1F(Form("TimeDiff/H1D_kNTimeDiffADC_%d_%d",i,j),Form("Time Difference ADC module%d - module%d",i,j),kNTimeDiffBins,kTimeDiffMin,kTimeDiffMax);
		Add2RawsList(h1d,kNTimeDiffADC+nCorrelation, expert, !image, !saveCorr); iHisto++; nCorrelation++;
		}
	}
  
  AliDebug(AliQAv1::GetQADebugLevel(), Form("%d Histograms has been added to the Raws List",iHisto));
  //
  ClonePerTrigClass(AliQAv1::kRAWS); // this should be the last line
}

//____________________________________________________________________________
void AliADQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
  // Fills histograms with Raws, computes average ADC values dynamically (pedestal subtracted)
                  
					  
  // Check id histograms already created for this Event Specie
  if ( ! GetRawsData(kPedestalInt0) )
    InitRaws() ;

  rawReader->Reset() ; 
  AliADRawStream* rawStream  = new AliADRawStream(rawReader); 
  if(!(rawStream->Next())) return;  
 
  eventTypeType eventType = rawReader->GetType();

  Int_t    mulADA = 0 ; 
  Int_t    mulADC = 0 ; 
  Double_t timeADA =0., timeADC = 0.;
  Double_t weightADA =0., weightADC = 0.;
  UInt_t   itimeADA=0, itimeADC=0;
  Double_t chargeADA=0., chargeADC=0.;

  Double_t diffTime=-100000.;
  
  Int_t	   pBBmulADA = 0;
  Int_t	   pBBmulADC = 0;
  Int_t	   pBGmulADA = 0;
  Int_t	   pBGmulADC = 0;
  Double_t pDiffTime =-100000.;

  
  switch (eventType){
  case PHYSICS_EVENT:
  
    Int_t  iFlag=0;
    Int_t  pedestal;
    Float_t OCDBdiff;
    Int_t  integrator[16];
    Bool_t flagBB[16];	 
    Bool_t flagBG[16];	 
    Float_t charge;
    Int_t  offlineCh;
    Float_t adc[16], time[16], width[16], timeCorr[16]; 
    Int_t  iPair=0;

    for(Int_t iChannel=0; iChannel<16; iChannel++) { // BEGIN : Loop over channels
		   
      offlineCh = kOfflineChannel[iChannel];
		   
      // Fill Pedestal histograms
	   
      for(Int_t j=15; j<21; j++) {
	if((rawStream->GetBGFlag(iChannel,j) || rawStream->GetBBFlag(iChannel,j))) iFlag++;
      }

      if(iFlag == 0){ //No Flag found
	for(Int_t j=15; j<21; j++){
	  pedestal= (Int_t) rawStream->GetPedestal(iChannel, j);
	  integrator[offlineCh] = rawStream->GetIntegratorFlag(iChannel, j);
	  OCDBdiff = pedestal - fCalibData->GetPedestal(offlineCh+16*integrator[offlineCh]);

	  FillRawsData((integrator[offlineCh] == 0 ? kPedestalInt0 : kPedestalInt1),offlineCh,pedestal);
	  FillRawsData((integrator[offlineCh] == 0 ? kPedestalDiffInt0 : kPedestalDiffInt1),offlineCh,OCDBdiff);
	}
      }
      // Fill Charge EoI histograms
	   
      adc[offlineCh]    = 0.0;
      // Search for the maximum charge in the train of 21 LHC clocks 
      // regardless of the integrator which has been operated:
      Float_t maxadc = 0;
      Int_t imax = -1;
      Float_t adcPedSub[21];
      for(Int_t iClock=0; iClock<21; iClock++){
	Bool_t iIntegrator = rawStream->GetIntegratorFlag(iChannel,iClock);
	Int_t k = offlineCh+16*iIntegrator;

	adcPedSub[iClock] = rawStream->GetPedestal(iChannel,iClock) - fCalibData->GetPedestal(k);
	if(adcPedSub[iClock] <= fRecoParam->GetNSigmaPed()*fCalibData->GetSigma(k)) {
	  adcPedSub[iClock] = 0;
	  continue;
	}
	
	if(iClock < fRecoParam->GetStartClock() || iClock > fRecoParam->GetEndClock()) continue;
	if(adcPedSub[iClock] > maxadc) {
	  maxadc = adcPedSub[iClock];
	  imax   = iClock;
	}
      }
      //printf(Form("Channel %d (online), %d (offline)\n",iChannel,j)); 
      if (imax != -1) {
	Int_t start = imax - fRecoParam->GetNPreClocks();
	if (start < 0) start = 0;
	Int_t end = imax + fRecoParam->GetNPostClocks();
	if (end > 20) end = 20;
	for(Int_t iClock = start; iClock <= end; iClock++) {
	  adc[offlineCh] += adcPedSub[iClock];
	}
      }
	
		
      Int_t iClock  = imax;
      charge = rawStream->GetPedestal(iChannel,iClock); // Charge at the maximum 

      integrator[offlineCh]    = rawStream->GetIntegratorFlag(iChannel,iClock);
      flagBB[offlineCh]	 = rawStream->GetBBFlag(iChannel, iClock);
      flagBG[offlineCh]	 = rawStream->GetBGFlag(iChannel,iClock );
      Int_t board = AliADCalibData::GetBoardNumber(offlineCh);
      time[offlineCh] = rawStream->GetTime(iChannel)*fCalibData->GetTimeResolution(board);
      width[offlineCh] = rawStream->GetWidth(iChannel)*fCalibData->GetWidthResolution(board);

      if (time[offlineCh] >= 1e-6) FillRawsData(kChargeEoI,offlineCh,adc[offlineCh]);

      FillRawsData((integrator[offlineCh] == 0 ? kChargeEoIInt0 : kChargeEoIInt1),offlineCh,charge);

      Float_t sigma = fCalibData->GetSigma(offlineCh+16*integrator[offlineCh]);
		  
      if((adc[offlineCh] > fRecoParam->GetNSigmaPed()*sigma) && !(time[offlineCh] <1.e-6)){ 
	if(offlineCh<8) {
	  mulADC++;
	  chargeADC += adc[offlineCh];
	  
	} else {
	  mulADA++;
	  chargeADA += adc[offlineCh];
	}
      }
		   
  
      // Fill HPTDC Time Histograms
      timeCorr[offlineCh] = CorrectLeadingTime(offlineCh,time[offlineCh],adc[offlineCh]);
      //timeCorr[offlineCh] = time[offlineCh];

      //const Float_t p1 = 2.50; // photostatistics term in the time resolution
      //const Float_t p2 = 3.00; // sleewing related term in the time resolution
      if(timeCorr[offlineCh]>-1024 + 1.e-6){
	//Float_t nphe = adc[offlineCh]*kChargePerADC/(fCalibData->GetGain(offlineCh)*TMath::Qe());
	Float_t timeErr = 1;
	/*/
	if (nphe>1.e-6) timeErr = TMath::Sqrt(kIntTimeRes*kIntTimeRes+
					      p1*p1/nphe+
					      p2*p2*(fTimeSlewing->GetParameter(0)*fTimeSlewing->GetParameter(1))*(fTimeSlewing->GetParameter(0)*fTimeSlewing->GetParameter(1))*
					      TMath::Power(adc[offlineCh]/fCalibData->GetCalibDiscriThr(offlineCh,kTRUE),2.*(fTimeSlewing->GetParameter(1)-1.))/
					      (fCalibData->GetCalibDiscriThr(offlineCh,kTRUE)*fCalibData->GetCalibDiscriThr(offlineCh,kTRUE)));/*/

	if (timeErr>1.e-6) {
	  if (offlineCh<8) {
	    itimeADC++;
	    timeADC += timeCorr[offlineCh]/(timeErr*timeErr);
	    weightADC += 1./(timeErr*timeErr);
	  }else{
	    itimeADA++;
	    timeADA += timeCorr[offlineCh]/(timeErr*timeErr);
	    weightADA += 1./(timeErr*timeErr);
	  }
	}
      }
      
      // Fill Flag and Charge Versus LHC-Clock histograms
      Int_t nbbFlag = 0;
      Int_t nbgFlag = 0;
      
      for(Int_t iEvent=0; iEvent<21; iEvent++){
	charge = rawStream->GetPedestal(iChannel,iEvent);
	Int_t intgr = rawStream->GetIntegratorFlag(iChannel,iEvent);
	Bool_t bbFlag	  = rawStream->GetBBFlag(iChannel, iEvent);
	Bool_t bgFlag	  = rawStream->GetBGFlag(iChannel,iEvent );
	if(bbFlag) nbbFlag++;
	if(bgFlag) nbgFlag++;
	
	FillRawsData((intgr == 0 ? kChargeVsClockInt0 : kChargeVsClockInt1 ), offlineCh,(float)iEvent-10,(float)charge);
	FillRawsData(kBBFlagVsClock, offlineCh,(float)iEvent-10,(float)bbFlag);
	FillRawsData(kBGFlagVsClock, offlineCh,(float)iEvent-10,(float)bgFlag);
	
      }
      FillRawsData(kBBFlagsPerChannel, offlineCh,nbbFlag);
      FillRawsData(kBGFlagsPerChannel, offlineCh,nbgFlag);
      if((nbbFlag+nbgFlag)>0 && time[offlineCh]<1e-6)FillRawsData(kFlagNoTime,offlineCh);
      if((nbbFlag+nbgFlag)==0 && time[offlineCh]>1e-6)FillRawsData(kTimeNoFlag,offlineCh);
      
      FillRawsData(kTimeSlewingOff,time[offlineCh],adc[offlineCh]);
      FillRawsData(kTimeSlewingOn,timeCorr[offlineCh],adc[offlineCh]);
      FillRawsData(kWidthSlewing,width[offlineCh],adc[offlineCh]);
      
      FillRawsData(kHPTDCTime,offlineCh,time[offlineCh]);
      FillRawsData(kWidth,offlineCh,width[offlineCh]);
      //if(flagBB[offlineCh]) {
      if(nbbFlag > 0){
	FillRawsData(kHPTDCTimeBB,offlineCh,time[offlineCh]);
	FillRawsData(kWidthBB,offlineCh,width[offlineCh]);
      }
      //if(flagBG[offlineCh]) {
      if(nbgFlag > 0){
	FillRawsData(kHPTDCTimeBG,offlineCh,time[offlineCh]);
	FillRawsData(kWidthBG,offlineCh,width[offlineCh]);
      }
      

    }// END of Loop over channels
    
    //Correlation Cside
    Int_t nCorrelation = 0;
    for(Int_t iChannel=0; iChannel<8; iChannel++) {
    	for(Int_t jChannel=7; jChannel>iChannel; jChannel--) {
		FillRawsData(kNChargeCorrADC+nCorrelation,adc[iChannel],adc[jChannel]);
		FillRawsData(kNTimeCorrADC+nCorrelation,time[iChannel],time[jChannel]);
		if(time[iChannel]>1e-6 && time[jChannel]>1e-6) FillRawsData(kNTimeDiffADC+nCorrelation,time[iChannel]-time[jChannel]);
		nCorrelation++;
		}
	}
    //Correlation Aside
    nCorrelation = 0;
    for(Int_t iChannel=8; iChannel<16; iChannel++) {
    	for(Int_t jChannel=15; jChannel>iChannel; jChannel--) {
		FillRawsData(kNChargeCorrADA+nCorrelation,adc[iChannel],adc[jChannel]);
		FillRawsData(kNTimeCorrADA+nCorrelation,time[iChannel],time[jChannel]);
		if(time[iChannel]>1e-6 && time[jChannel]>1e-6) FillRawsData(kNTimeDiffADA+nCorrelation,time[iChannel]-time[jChannel]);
		nCorrelation++;
		}
    	}
	
    for(Int_t iChannel=0; iChannel<4; iChannel++) {//Loop over pairs
    		if(flagBB[iChannel] && flagBB[iChannel+4]) pBBmulADC++;
		if(flagBB[iChannel+8] && flagBB[iChannel+12]) pBBmulADA++;
		if(flagBG[iChannel] && flagBG[iChannel+4]) pBGmulADC++;
		if(flagBG[iChannel+8] && flagBG[iChannel+12]) pBGmulADA++;
		}	
    FillRawsData(kNBBCoincADA,pBBmulADA);
    FillRawsData(kNBBCoincADC,pBBmulADC);
    FillRawsData(kNBGCoincADA,pBGmulADA);
    FillRawsData(kNBGCoincADC,pBGmulADC);
	
	
    if(weightADA>1.e-6) timeADA /= weightADA; 
    else timeADA = -1024.;
    if(weightADC>1.e-6) timeADC /= weightADC;
    else timeADC = -1024.;
    if(timeADA<-1024.+1.e-6 || timeADC<-1024.+1.e-6) diffTime = -1024.;
    else diffTime = timeADA - timeADC;
    

		
    FillRawsData(kADATime,timeADA);
    FillRawsData(kADCTime,timeADC);
    FillRawsData(kDiffTime,diffTime);
    FillRawsData(kTimeADAADC,timeADA,timeADC);

    FillRawsData(kMultiADA,mulADA);
    FillRawsData(kMultiADC,mulADC);

    FillRawsData(kChargeADA,chargeADA);
    FillRawsData(kChargeADC,chargeADC);
    FillRawsData(kChargeAD,chargeADA + chargeADC);
	    
    break;
  } // END of SWITCH : EVENT TYPE 
	
  TParameter<double> * p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), GetRawsData(kMultiADA)->GetName()))) ; 
  if (p) p->SetVal((double)mulADA) ; 

  p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), GetRawsData(kMultiADC)->GetName()))) ; 
  if (p) p->SetVal((double)mulADC) ;                     

  p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), GetRawsData(kChargeADA)->GetName()))) ; 
  if (p) p->SetVal((double)chargeADA) ; 

  p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), GetRawsData(kChargeADC)->GetName()))) ; 
  if (p) p->SetVal((double)chargeADC) ;                     

  p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), GetRawsData(kChargeAD)->GetName()))) ; 
  if (p) p->SetVal((double)(chargeADA + chargeADC)) ;                     
	                   	
  p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), GetRawsData(kADATime)->GetName()))) ; 
  if (p) p->SetVal((double)timeADA) ; 
	
  p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), GetRawsData(kADCTime)->GetName()))) ; 
  if (p) p->SetVal((double)timeADC) ;                     
	
  p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), GetRawsData(kDiffTime)->GetName()))) ; 
  if (p) p->SetVal((double)diffTime) ;                     
	
  delete rawStream; rawStream = 0x0;      
  //
  IncEvCountCycleRaws();
  IncEvCountTotalRaws();
  //
}

//____________________________________________________________________________ 
Float_t AliADQADataMakerRec::CorrectLeadingTime(Int_t /*i*/, Float_t time, Float_t adc) const
{
  // Correct the leading time
  // for slewing effect and
  // misalignment of the channels
  if (time < 1e-6) return -1024;

  // Channel alignment and general offset subtraction
  //  time -= fTimeOffset[i];
  //AliInfo(Form("time-offset %f", time));

  // In case of pathological signals
  //if (adc < 1e-6) return time;

  // Slewing correction
  //Float_t thr = fCalibData->GetCalibDiscriThr(i,kTRUE);
  //AliInfo(Form("adc %f thr %f dtime %f ", adc,thr,fTimeSlewing->Eval(adc/thr)));
  time -= fTimeSlewing->Eval(adc);

  return time;
}

