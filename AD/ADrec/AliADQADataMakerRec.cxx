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
#include <TPaveText.h>

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
#include "AliADQAParam.h"
#include "AliCTPTimeParams.h"
#include "AliLHCClockPhase.h"
#include "event.h"

ClassImp(AliADQADataMakerRec)
           
//____________________________________________________________________________ 
AliADQADataMakerRec::AliADQADataMakerRec() : 
AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kAD), "AD Quality Assurance Data Maker"),
  fCalibData(0x0),
  fRecoParam(0x0),
  fQAParam(0x0),
  fTrendingUpdateTime(0), 
  fCycleStartTime(0), 
  fCycleStopTime(0),
  fADADist(56.7),
  fADCDist(65.19),
  fOldRun(0)
    
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
  fQAParam(0x0),
  fTrendingUpdateTime(0), 
  fCycleStartTime(0), 
  fCycleStopTime(0),
  fADADist(56.7),
  fADCDist(65.19),
  fOldRun(0)
  
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
AliADQAParam* AliADQADataMakerRec::GetQAParam() const

{
  AliCDBManager *man = AliCDBManager::Instance();

  AliCDBEntry *entry=0;

  entry = man->Get("AD/Calib/QAParam",fRun);
  if(!entry){
    AliWarning("Load of QA param from default storage failed!");
    AliWarning("QA parameters will be loaded from local storage ($ALICE_ROOT)");
	
    man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    entry = man->Get("AD/Calib/QAParam",fRun);
  }
  // Retrieval of data in directory AD/Calib/QA:

  AliADQAParam *QAParam = 0;

  if (entry) QAParam = (AliADQAParam*) entry->GetObject();
  if (!QAParam)  AliFatal("No QA param from calibration database !");

  return QAParam;
}
//____________________________________________________________________________ 
void AliADQADataMakerRec::StartOfDetectorCycle()
{
  // Detector specific actions at start of cycle
  
  // Reset of the histogram used - to have the trend versus time -
 
  fCalibData = GetCalibData();
  fQAParam = GetQAParam();
  if(!fRecoParam)fRecoParam = (AliADRecoParam*)GetRecoParam();
	
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
    TDatime currentTime;
    fCycleStopTime = currentTime.GetSecond();
    
    if (fRun!=fOldRun){
    ((TH1F*)GetRawsData(kTrend_TriggerChargeQuantileADA))->SetBins(1,0,1);
    ((TH1F*)GetRawsData(kTrend_TriggerChargeQuantileADC))->SetBins(1,0,1);
    fOldRun=fRun;
    }
    
    Double_t xq[1] = {0.9};
    Double_t yq[1];
    UInt_t currentBins = ((TH1F*)GetRawsData(kTrend_TriggerChargeQuantileADC))->GetNbinsX();
    
    ((TH1F*)GetRawsData(kTrend_TriggerChargeQuantileADC))->SetBins(currentBins+1,0,currentBins+1);
    ((TH1F*)GetRawsData(kChargeADC_PC))->GetQuantiles(1,yq,xq);
    
    ((TH1F*)GetRawsData(kTrend_TriggerChargeQuantileADC))->SetBinContent(currentBins,yq[0]);
    ((TH1F*)GetRawsData(kTrend_TriggerChargeQuantileADC))->GetXaxis()->LabelsOption("v");
    if (currentBins%10 == 1)((TH1F*)GetRawsData(kTrend_TriggerChargeQuantileADC))->GetXaxis()->SetBinLabel(currentBins,Form("%d:%02d:%02d",currentTime.GetHour(),currentTime.GetMinute(),currentTime.GetSecond()));
    ((TH1F*)GetRawsData(kChargeADC_PC))->Reset("ICES");
    
    ((TH1F*)GetRawsData(kTrend_TriggerChargeQuantileADA))->SetBins(currentBins+1,0,currentBins+1);
    ((TH1F*)GetRawsData(kChargeADA_PC))->GetQuantiles(1,yq,xq);
    
    ((TH1F*)GetRawsData(kTrend_TriggerChargeQuantileADA))->SetBinContent(currentBins,yq[0]);
    ((TH1F*)GetRawsData(kTrend_TriggerChargeQuantileADA))->GetXaxis()->LabelsOption("v");
    if (currentBins%10 == 1)((TH1F*)GetRawsData(kTrend_TriggerChargeQuantileADA))->GetXaxis()->SetBinLabel(currentBins,Form("%d:%02d:%02d",currentTime.GetHour(),currentTime.GetMinute(),currentTime.GetSecond()));
    ((TH1F*)GetRawsData(kChargeADA_PC))->Reset("ICES");
    
    
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
   

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    if (! IsValidEventSpecie(specie, list)) continue ;
    SetEventSpecie(AliRecoParam::ConvertIndex(specie));
    if(task == AliQAv1::kRAWS) {
    	AliQAChecker::Instance()->Run(AliQAv1::kAD, task, list) ;
    } else if (task == AliQAv1::kESDS) {
    }
  }
  
    }
  
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
  
  TH2F * h5 = new TH2F("H2D_Time_Channel", "Time per channel;Channel;Time (ns)",16,0,16, 1638, -79.980469, 79.980469);  
  Add2ESDsList(h5, kTimeChannel, !expert, image)  ;  
  
  TH1F * h6 = new TH1F("H1D_ADA_Time", "Mean ADA Time;Time (ns);Counts",1638, -79.980469, 79.980469);
  Add2ESDsList(h6,kESDADATime, !expert, image); 
  
  TH1F * h7 = new TH1F("H1D_ADC_Time", "Mean ADC Time;Time (ns);Counts",1638, -79.980469, 79.980469);
  Add2ESDsList(h7,kESDADCTime, !expert, image); 
  
  TH1F * h8 = new TH1F("H1D_Diff_Time", "Diff Time ADA - ADC;Diff Time ADA - ADC (ns);Counts",1000, -200., 200.);
  Add2ESDsList(h8,kESDDiffTime, !expert, image); 
  
  TH2F * h9 = new TH2F("H2D_ADA_TimeVsCharge", "TimeVsCharge ADA;Time (ns); Charge(ADC counts);Counts",1638, -79.980469, 79.980469,5000,0,5000);
  Add2ESDsList(h9,kESDADATimeVsCharge, !expert, image);
  
  TH2F * h10 = new TH2F("H2D_ADC_TimeVsCharge", "TimeVsCharge ADC;Time (ns); Charge(ADC counts);Counts",1638, -79.980469, 79.980469,5000,0,5000);
  Add2ESDsList(h10,kESDADCTimeVsCharge, !expert, image);
  
  TH2F * h11 = new TH2F("H2D_ADA_PairTimeSumDiff", "Pair Time Sum Vs Diff ADA; t1+t2 (ns); t1-t2 (ns);Counts",82, 79.980469, 160.058594,410, 0.000000, 40.039062);
  Add2ESDsList(h11,kESDADAPairTimeSumDiff, !expert, image);
  
  TH2F * h12 = new TH2F("H2D_ADC_PairTimeSumDiff", "Pair Time Sum Vs Diff ADC; t1+t2 (ns); t1-t2 (ns);Counts",82, 79.980469, 160.058594,410, 0.000000, 40.039062);
  Add2ESDsList(h12,kESDADCPairTimeSumDiff, !expert, image);
  
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
     
  TH2D * h1 = new TH2D("hDigitLeadingTimePerPM", "Leading time distribution per PM in AD;PM number;Leading Time [ns]",16,0,16, 3062, 0.976562, 300); 
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
    digitTree->ResetBranchAddress(branch);  
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
	if(esdAD->BBTriggerADC(i)) {
		FillESDsData(kBBFlag,(Float_t) i);
		FillESDsData(kESDADCTimeVsCharge,esdAD->GetTime(i),esdAD->GetAdc(i));
		}
	if(esdAD->BGTriggerADC(i)) FillESDsData(kBGFlag,(Float_t) i);
	
      }
      else {
	if(esdAD->BBTriggerADA(i-8)){ 
		FillESDsData(kBBFlag,(Float_t) i); 
		FillESDsData(kESDADATimeVsCharge,esdAD->GetTime(i),esdAD->GetAdc(i));
		} 
	if(esdAD->BGTriggerADA(i-8)) FillESDsData(kBGFlag,(Float_t) i);
	
      }		  	
      Float_t time = (Float_t) esdAD->GetTime(i);
      FillESDsData(kTimeChannel,(Float_t) i,time);
    }
    
    for (Int_t i = 0; i < 4; ++i) {
      Float_t time1 = esdAD->GetTime(i);
      Float_t time2 = esdAD->GetTime(i+4);
      if(time1<-1024.+1.e-6 || time2<-1024.+1.e-6) continue;
      Float_t timeDiff = TMath::Abs(time1-time2);
      Float_t timeSum = time1+time2;
      FillESDsData(kESDADCPairTimeSumDiff,timeSum,timeDiff);
	}
		
     for (Int_t i = 8; i < 12; ++i) {
      Float_t time1 = esdAD->GetTime(i);
      Float_t time2 = esdAD->GetTime(i+4);
      if(time1<-1024.+1.e-6 || time2<-1024.+1.e-6) continue;
      Float_t timeDiff = TMath::Abs(time1-time2);
      Float_t timeSum = time1+time2;
      FillESDsData(kESDADAPairTimeSumDiff,timeSum,timeDiff);
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
  if(!fRecoParam)fRecoParam = (AliADRecoParam*)GetRecoParam();
  if(!fQAParam) fQAParam = (AliADQAParam*)GetQAParam();
 
  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  const Int_t kNintegrator  =    2;
 
  const Int_t kNTdcTimeBins  = fQAParam->GetNTdcTimeBins();
  const Float_t kTdcTimeMin    =  fQAParam->GetTdcTimeMin();
  const Float_t kTdcTimeMax    = fQAParam->GetTdcTimeMax();
  const Int_t kNTdcTimeBinsFlag  = fQAParam->GetNTdcTimeBinsFlag();
  const Float_t kTdcTimeMinBBFlag    =  fQAParam->GetTdcTimeMinBBFlag();
  const Float_t kTdcTimeMaxBBFlag    =  fQAParam->GetTdcTimeMaxBBFlag();
  const Float_t kTdcTimeMinBGFlag    =  fQAParam->GetTdcTimeMinBGFlag();
  const Float_t kTdcTimeMaxBGFlag    =  fQAParam->GetTdcTimeMaxBGFlag();
  const Int_t kNTdcTimeRatioBins  = fQAParam->GetNTdcTimeRatioBins();
  const Float_t kTdcTimeRatioMin    =  fQAParam->GetTdcTimeRatioMin();
  const Float_t kTdcTimeRatioMax    = fQAParam->GetTdcTimeRatioMax();
  
  const Int_t kNTdcWidthBins =  fQAParam->GetNTdcWidthBins();
  const Float_t kTdcWidthMin   =    fQAParam->GetTdcWidthMin();
  const Float_t kTdcWidthMax   =  fQAParam->GetTdcWidthMax();
  
  const Int_t kNChargeChannelBins   =  fQAParam->GetNChargeChannelBins();
  const Int_t kChargeChannelMin   =  fQAParam->GetChargeChannelMin();
  const Int_t kChargeChannelMax   =  fQAParam->GetChargeChannelMax();
  
  const Int_t kNChargeSideBins   = fQAParam->GetNChargeSideBins();
  const Int_t kChargeSideMin   = fQAParam->GetChargeSideMin();
  const Int_t kChargeSideMax   = fQAParam->GetChargeSideMax();
  
  const Int_t kNChargeCorrBins   = fQAParam->GetNChargeCorrBins();
  const Int_t kChargeCorrMin   = fQAParam->GetChargeCorrMin();
  const Int_t kChargeCorrMax   = fQAParam->GetChargeCorrMax();
     
  const Int_t   kNPairTimeCorrBins = fQAParam->GetNPairTimeCorrBins(); 
  const Float_t kPairTimeCorrMin =  fQAParam->GetPairTimeCorrMin();
  const Float_t kPairTimeCorrMax =  fQAParam->GetPairTimeCorrMax(); 
   
  const Int_t kNPairTimeDiffBins = fQAParam->GetNPairTimeDiffBins();
  const Float_t kPairTimeDiffMin = fQAParam->GetPairTimeDiffMin();
  const Float_t kPairTimeDiffMax = fQAParam->GetPairTimeDiffMax(); 
  
  const Int_t   kNMeanTimeCorrBins = fQAParam->GetNMeanTimeCorrBins();
  const Float_t kMeanTimeCorrMin =   fQAParam->GetMeanTimeCorrMin();
  const Float_t kMeanTimeCorrMax =   fQAParam->GetMeanTimeCorrMax();    
  
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
  
  h1d = new TH1F("H1D_Trend_TriggerChargeQuantileADA","Trigger charge quantile",1, 0, 1) ;  
  Add2RawsList(h1d,kTrend_TriggerChargeQuantileADA, !expert, image, saveCorr);   iHisto++;
  h1d->SetLineWidth(2);
  h1d->SetLineColor(kBlue);
  
  h1d = new TH1F("H1D_Charge_ADA_PC",Form("Total integrated [-%d,+%d] charge;Charge [ADC counts]",fRecoParam->GetNPreClocks(),fRecoParam->GetNPostClocks()), kNChargeSideBins, kChargeSideMin, kChargeSideMax) ;  
  Add2RawsList(h1d,kChargeADA_PC, !expert, image, saveCorr);   iHisto++;
  
  h1d = new TH1F("H1D_Trend_TriggerChargeQuantileADC","Trigger charge quantile",1, 0, 1) ;  
  Add2RawsList(h1d,kTrend_TriggerChargeQuantileADC, !expert, image, saveCorr);   iHisto++;
  h1d->SetLineWidth(2);
  h1d->SetLineColor(kRed);
  
  h1d = new TH1F("H1D_Charge_ADC_PC",Form("Total integrated [-%d,+%d] charge;Charge [ADC counts]",fRecoParam->GetNPreClocks(),fRecoParam->GetNPostClocks()), kNChargeSideBins, kChargeSideMin, kChargeSideMax) ;  
  Add2RawsList(h1d,kChargeADC_PC, !expert, image, saveCorr);   iHisto++;
  
 
  // Creation of Cell Multiplicity Histograms
  h1i = new TH1I("H1I_Multiplicity_ADA", "Number of channels with charge signal and time ADA;# of Channels;Entries", 9, -0.5, 8.5) ;  
  Add2RawsList(h1i,kMultiADA, expert, !image, !saveCorr);   iHisto++;
  h1i = new TH1I("H1I_Multiplicity_ADC", "Number of channels with charge signal and time ADC;# of Channels;Entries", 9, -0.5, 8.5) ;  
  Add2RawsList(h1i,kMultiADC, expert, !image, !saveCorr);   iHisto++;
 
  // Creation of Total Charge Histograms
  h1d = new TH1F("H1D_Charge_ADA",Form("Total integrated [-%d,+%d] charge;Charge [ADC counts]",fRecoParam->GetNPreClocks(),fRecoParam->GetNPostClocks()), kNChargeSideBins, kChargeSideMin, kChargeSideMax) ;  
  Add2RawsList(h1d,kChargeADA, !expert, image, saveCorr);   iHisto++;
  h1d->SetLineWidth(2);
  h1d->SetLineColor(kBlue);
  h1d = new TH1F("H1D_Charge_ADC",Form("Total integrated [-%d,+%d] charge;Charge [ADC counts]",fRecoParam->GetNPreClocks(),fRecoParam->GetNPostClocks()), kNChargeSideBins, kChargeSideMin, kChargeSideMax) ;  
  Add2RawsList(h1d,kChargeADC, !expert, image, saveCorr);   iHisto++;
  h1d->SetLineWidth(2);
  h1d->SetLineColor(kRed);
  h1d = new TH1F("H1D_Charge_AD",Form("Total integrated [-%d,+%d] charge;Charge [ADC counts]",fRecoParam->GetNPreClocks(),fRecoParam->GetNPostClocks()), 2*kNChargeSideBins, kChargeSideMin, 1+2*kNChargeSideBins) ;  
  Add2RawsList(h1d,kChargeAD, !expert,  !image, !saveCorr);   iHisto++;
   

  // Creation of Charge EoI histogram 
  h2d = new TH2F("H2D_ChargeEoI", Form("Integrated [-%d,+%d] charge;Channel Number;Charge [ADC counts]",fRecoParam->GetNPreClocks(),fRecoParam->GetNPostClocks())
		 ,kNChannelBins, kChannelMin, kChannelMax, kNChargeChannelBins, kChargeChannelMin, kChargeChannelMax);
  Add2RawsList(h2d,kChargeEoI, !expert, image, saveCorr); iHisto++;
  
  h2d = new TH2F("H2D_ChargeEoIBB", Form("Integrated [-%d,+%d] charge w/ BB Flag condition;Channel Number;Charge [ADC counts]",fRecoParam->GetNPreClocks(),fRecoParam->GetNPostClocks())
		 ,kNChannelBins, kChannelMin, kChannelMax, kNChargeChannelBins, kChargeChannelMin, kChargeChannelMax);
  Add2RawsList(h2d,kChargeEoIBB, !expert, image, saveCorr); iHisto++;
  
  h2d = new TH2F("H2D_ChargeEoIBG", Form("Integrated [-%d,+%d] charge w/ BG Flag condition;Channel Number;Charge [ADC counts]",fRecoParam->GetNPreClocks(),fRecoParam->GetNPostClocks())
		 ,kNChannelBins, kChannelMin, kChannelMax, kNChargeChannelBins, kChargeChannelMin, kChargeChannelMax);
  Add2RawsList(h2d,kChargeEoIBG, !expert, image, saveCorr); iHisto++;


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
  
  h2i = new TH2I("H2I_ChargeSaturation", "Maximum charge per clock, both Ints;Channel;Charge [ADC counts]",kNChannelBins, kChannelMin, kChannelMax, 1024, 1, 1025);
  Add2RawsList(h2i,kChargeSaturation, !expert, image, saveCorr); iHisto++;	
  
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
  
  h2i = new TH2I("H2I_HPTDCTime_BB", "HPTDC Time w/ BB Flag condition;Channel;Leading Time [ns]",kNChannelBins, kChannelMin, kChannelMax, kNTdcTimeBinsFlag, kTdcTimeMinBBFlag, kTdcTimeMaxBBFlag);
  Add2RawsList(h2i,kHPTDCTimeBB, !expert, image, !saveCorr); iHisto++;

  h2i = new TH2I("H2I_HPTDCTime_BG", "HPTDC Time w/ BG Flag condition;Channel;Leading Time [ns]",kNChannelBins, kChannelMin, kChannelMax, kNTdcTimeBinsFlag, kTdcTimeMinBGFlag, kTdcTimeMaxBGFlag);
  Add2RawsList(h2i,kHPTDCTimeBG, !expert, image, !saveCorr); iHisto++;
  
  //With wide binning for ratio
  h2d = new TH2F("H2D_HPTDCTimeRebin", "HPTDC Time;Channel;Leading Time [ns]",kNChannelBins, kChannelMin, kChannelMax, kNTdcTimeRatioBins, kTdcTimeRatioMin, kTdcTimeRatioMax);
  Add2RawsList(h2d,kHPTDCTimeRebin, !expert, image, saveCorr); iHisto++;
  
  h2d = new TH2F("H2D_HPTDCTimeRebin_BB", "Ratio Time w_BB_Flag/All ;Channel;Leading Time [ns]",kNChannelBins, kChannelMin, kChannelMax, kNTdcTimeRatioBins, kTdcTimeRatioMin, kTdcTimeRatioMax);
  Add2RawsList(h2d,kHPTDCTimeRebinBB, !expert, image, !saveCorr); iHisto++;

  h2d = new TH2F("H2D_HPTDCTimeRebin_BG", "Ratio Time w_BG_Flag/All;Channel;Leading Time [ns]",kNChannelBins, kChannelMin, kChannelMax, kNTdcTimeRatioBins, kTdcTimeRatioMin, kTdcTimeRatioMax);
  Add2RawsList(h2d,kHPTDCTimeRebinBG, !expert, image, !saveCorr); iHisto++;

  //Mean time histograms	
  h1d = new TH1F("H1D_MeanTimeADA", "Mean Time;Mean time [ns];Counts",kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax);
  Add2RawsList(h1d,kMeanTimeADA, expert, !image, !saveCorr); iHisto++;
  h1d->SetLineWidth(2);
  h1d->SetLineColor(kBlue);
	
  h1d = new TH1F("H1D_MeanTimeADC", "Mean Time;Mean time [ns];Counts",kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax);
  Add2RawsList(h1d,kMeanTimeADC, expert, !image, !saveCorr); iHisto++;
  h1d->SetLineWidth(2);
  h1d->SetLineColor(kRed);
	
  h1d = new TH1F("H1D_MeanTimeDifference","Mean Time Difference ADA-ADC ;AD Mean time t_{A} - t_{C} [ns];Counts",1024,-150,150);
  Add2RawsList(h1d,kMeanTimeDiff, expert, !image, !saveCorr); iHisto++;

  h2d = new TH2F("H2D_MeanTimeCorr", "AD Mean time t_{A} vs t_{C};Mean time ADA [ns];Mean time ADC [ns]", kNMeanTimeCorrBins,kMeanTimeCorrMin,kMeanTimeCorrMax,kNMeanTimeCorrBins, kMeanTimeCorrMin,kMeanTimeCorrMax) ;  
  Add2RawsList(h2d,kMeanTimeCorr, expert, !image, !saveCorr);   iHisto++;
 
  h2d = new TH2F("H2D_MeanTimeSumDiff", "AD Mean time t_{A} - t_{C} vs t_{A} + t_{C}; AD Mean time t_{A} - t_{C} [ns];AD Mean time t_{A} + t_{C} [ns]", 307, -150.000000, 149.804688, 410, 200.0, 600.390625);  
  Add2RawsList(h2d,kMeanTimeSumDiff, expert, !image, !saveCorr);   iHisto++;
 
  //Slewing histograms
  h2d = new TH2F("H2D_TimeSlewingADA", "Time Vs Charge ADA;Log10(1/Charge) [ADC counts]; Leading Time[ns]", 200,-4,0, kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax);  
  Add2RawsList(h2d,kTimeSlewingADA, expert, !image, !saveCorr);   iHisto++;
  
  h2d = new TH2F("H2D_TimeSlewingADC", "Time Vs Charge ADC;Log10(1/Charge) [ADC counts]; Leading Time[ns]", 200,-4,0, kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax) ;  
  Add2RawsList(h2d,kTimeSlewingADC, expert, !image, !saveCorr);   iHisto++;
  
  h2d = new TH2F("H2D_WidthSlewing", "Width Vs Charge ;Time Width [ns];Charge [ADC counts]", kNTdcWidthBins, kTdcWidthMin, kTdcWidthMax, kNChargeCorrBins, kChargeCorrMin, kChargeCorrMax) ;  
  Add2RawsList(h2d,kWidthSlewing, expert, !image, !saveCorr);   iHisto++;
  
  //Creation of pair coincidence histograms
  h1i = new TH1I("H1I_MultiBBCoincidence_ADA", "Number of BB flag coincidences;# of BB Coincidences;Entries", 5, -0.5, 4.5) ;  
  Add2RawsList(h1i,kNBBCoincADA, !expert, image, saveCorr);   iHisto++;
  h1i->SetLineWidth(2);
  h1i->SetLineColor(kBlue);
  h1i->GetXaxis()->SetNdivisions(505);
  h1i = new TH1I("H1I_MultiBBCoincidence_ADC", "Number of BB flag coincidences;# of BB Coincidences;Entries", 5, -0.5, 4.5) ;  
  Add2RawsList(h1i,kNBBCoincADC, !expert, image, saveCorr);   iHisto++;
  h1i->SetLineWidth(2);
  h1i->SetLineColor(kRed);
  h1i->GetXaxis()->SetNdivisions(505);
  
  h1i = new TH1I("H1I_MultiBGCoincidence_ADA", "Number of BG flag coincidences;# of BG Coincidences;Entries", 5, -0.5, 4.5) ;  
  Add2RawsList(h1i,kNBGCoincADA, !expert, image, saveCorr);   iHisto++;
  h1i->SetLineWidth(2);
  h1i->SetLineColor(kBlue);
  h1i->GetXaxis()->SetNdivisions(505);
  h1i = new TH1I("H1I_MultiBGCoincidence_ADC", "Number of BG flag coincidences;# of BG Coincidences;Entries", 5, -0.5, 4.5) ;  
  Add2RawsList(h1i,kNBGCoincADC, !expert, image, saveCorr);   iHisto++;
  h1i->SetLineWidth(2);
  h1i->SetLineColor(kRed);
  h1i->GetXaxis()->SetNdivisions(505);
  
  h2i = new TH2I("H2I_BBCoincCorr", "Number of BB flag coincidences;# of BB Coincidences ADA;# of BB Coincidences ADC",5, -0.5, 4.5, 5, -0.5, 4.5);
  Add2RawsList(h2i,kNBBCoincCorr, !expert, image, !saveCorr); iHisto++;
  h2i->GetXaxis()->SetNdivisions(505);
  h2i->GetYaxis()->SetNdivisions(505);
  
  h2i = new TH2I("H2I_BGCoincCorr", "Number of BG flag coincidences;# of BG Coincidences ADA;# of BG Coincidences ADC",5, -0.5, 4.5, 5, -0.5, 4.5);
  Add2RawsList(h2i,kNBGCoincCorr, !expert, image, !saveCorr); iHisto++;
  h2i->GetXaxis()->SetNdivisions(505);
  h2i->GetYaxis()->SetNdivisions(505);
  
  h1d = new TH1F("H1D_NEventsBBFlag", "Number of events with a BB/BG flag per channel;# of Channels;Event rate", kNChannelBins, kChannelMin, kChannelMax);  
  Add2RawsList(h1d,kNEventsBBFlag, expert, image, saveCorr);   iHisto++;
  h1d->SetLineWidth(3);
  h1d->SetLineColor(kBlue);
  
  h1d = new TH1F("H1D_NEventsBGFlag", "Number of events with a BB/BG flag per channel;# of Channels;Event rate", kNChannelBins, kChannelMin, kChannelMax);  
  Add2RawsList(h1d,kNEventsBGFlag, expert, image, saveCorr);   iHisto++;
  h1d->SetLineWidth(3);
  h1d->SetLineColor(kRed);

  //Creation of trigger histogram
  h1d = new TH1F("H1D_Trigger_Type", "AD0 Trigger Type;;Counts", 11,0 ,11) ;  
  Add2RawsList(h1d,kTriggers, !expert, image, saveCorr);   iHisto++;
  h1d->SetFillColor(kAzure-8);
  h1d->SetLineWidth(2);
  h1d->GetXaxis()->SetLabelSize(0.04);
  h1d->GetXaxis()->SetNdivisions(808,kFALSE);
  h1d->GetXaxis()->SetBinLabel(1, "UBA");
  h1d->GetXaxis()->SetBinLabel(2, "UBC");
  h1d->GetXaxis()->SetBinLabel(3, "UGA");
  h1d->GetXaxis()->SetBinLabel(4, "UGC");
  h1d->GetXaxis()->SetBinLabel(5, "UBA & UBC");
  h1d->GetXaxis()->SetBinLabel(6, "UBA || UBC");
  h1d->GetXaxis()->SetBinLabel(7, "(UBA || UBC) & !(UGA || UGC)");
  h1d->GetXaxis()->SetBinLabel(8, "UGA & UBC");
  h1d->GetXaxis()->SetBinLabel(9, "UGC & UBA");
  h1d->GetXaxis()->SetBinLabel(10, "UGA || UGC");
  h1d->GetXaxis()->SetBinLabel(11, "(UGA & UBC) || (UGC & UBA)");
  
  h2d = new TH2F("H2D_Decision", "AD Decision; ADA; ADC", 4,0 ,4,4,0,4) ;  
  Add2RawsList(h2d,kDecisions, !expert, image, saveCorr);   iHisto++;
  h2d->SetOption("coltext");
  h2d->GetXaxis()->SetLabelSize(0.06);
  h2d->GetYaxis()->SetLabelSize(0.06);
  h2d->GetXaxis()->SetNdivisions(808,kFALSE);
  h2d->GetYaxis()->SetNdivisions(808,kFALSE);
  h2d->GetXaxis()->SetBinLabel(1, "Empty");
  h2d->GetXaxis()->SetBinLabel(2, "BB");
  h2d->GetXaxis()->SetBinLabel(3, "BG");
  h2d->GetXaxis()->SetBinLabel(4, "Fake");
  h2d->GetYaxis()->SetBinLabel(1, "Empty");
  h2d->GetYaxis()->SetBinLabel(2, "BB");
  h2d->GetYaxis()->SetBinLabel(3, "BG");
  h2d->GetYaxis()->SetBinLabel(4, "Fake");

  //Creation of debug histograms
  h1d = new TH1F("H1D_Pair_TimeDiffMean","Time difference mean for coincidence pair [ns];Pair number;Time mean [ns]",kNPairBins, kPairMin, kPairMax);
  Add2RawsList(h1d,kPairTimeDiffMean, expert, !image, !saveCorr); iHisto++;
  
  h1d = new TH1F("H1D_Pair_TimeDiffRMS","Time difference RMS for coincidence pair [ns];Pair number;Time RMS [ns]",kNPairBins, kPairMin, kPairMax);
  Add2RawsList(h1d,kPairTimeDiffRMS, expert, !image, !saveCorr); iHisto++;

  //Creation of Clock histograms
  h2d = new TH2F("H2D_BBFlagVsClock", "BB-Flags vs LHC-Clock(All events);Channel;LHC Clocks",kNChannelBins, kChannelMin, kChannelMax,21, -10.5, 10.5 );
  Add2RawsList(h2d,kBBFlagVsClock, !expert, image, saveCorr); iHisto++;
  
  h2d = new TH2F("H2D_BBFlagVsClock_ADOR", "BB-Flags vs LHC-Clock(AD-OR);Channel;LHC Clocks",kNChannelBins, kChannelMin, kChannelMax,21, -10.5, 10.5 );
  Add2RawsList(h2d,kBBFlagVsClock_ADOR, !expert, image, saveCorr); iHisto++;
	
  h2d = new TH2F("H2D_BGFlagVsClock", "BG-Flags vs LHC-Clock;Channel;LHC Clocks",kNChannelBins, kChannelMin, kChannelMax,21, -10.5, 10.5 );
  Add2RawsList(h2d,kBGFlagVsClock, !expert, image, saveCorr); iHisto++;
 
  for(Int_t iInt=0;iInt<kNintegrator;iInt++){
  	h2d = new TH2F(Form("H2D_ChargeVsClock_Int%d",iInt), Form("Charge Versus LHC-Clock (Int%d);Channel;LHCClock;Charge [ADC counts]",iInt),kNChannelBins, kChannelMin, kChannelMax,21, -10.5, 10.5 );
  	Add2RawsList(h2d,(iInt == 0 ? kChargeVsClockInt0 : kChargeVsClockInt1 ), expert, image, saveCorr); iHisto++;
	}
	
  h2d = new TH2F("H2D_MaxChargeVsClock", "Maximum Charge Versus LHC-Clock;Channel;LHC Clocks",kNChannelBins, kChannelMin, kChannelMax,21, -10.5, 10.5 );
  Add2RawsList(h2d,kMaxChargeClock, !expert, image, saveCorr); iHisto++;
  
  h2d = new TH2F("H2D_BBFlagPerChannel", "BB-Flags Versus Channel;Channel;BB Flags Count",kNChannelBins, kChannelMin, kChannelMax,22,-0.5,21.5);
  Add2RawsList(h2d,kBBFlagsPerChannel, expert, !image, saveCorr); iHisto++;
  
  h2d = new TH2F("H2D_BGFlagPerChannel", "BG-Flags Versus Channel;Channel;BG Flags Count",kNChannelBins, kChannelMin, kChannelMax,22,-0.5,21.5);
  Add2RawsList(h2d,kBGFlagsPerChannel, expert, !image, saveCorr); iHisto++;
  
  h1d = new TH1F("H1D_FlagNoTime", "No Flag-No Time events;Channel;Event rate",kNChannelBins, kChannelMin, kChannelMax);
  Add2RawsList(h1d,kFlagNoTime, !expert, image, saveCorr); iHisto++;
  h1d->SetLineWidth(3);
  h1d->SetLineColor(kBlue);
  
  h1d = new TH1F("H1D_TimeNoFlag", "No Flag-No Time events;Channel;Event rate",kNChannelBins, kChannelMin, kChannelMax);
  Add2RawsList(h1d,kTimeNoFlag, !expert, image, saveCorr); iHisto++;
  h1d->SetLineWidth(3);
  h1d->SetLineColor(kRed);
  
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
  		h2d = new TH2F(Form("TimeCorr/H2D_kNTimeCorrADA_%d_%d",i,j),Form("Time Correlation ADA module%d - module%d",i,j),kNPairTimeCorrBins,kPairTimeCorrMin,kPairTimeCorrMax,kNPairTimeCorrBins,kPairTimeCorrMin,kPairTimeCorrMax);
		Add2RawsList(h2d,kNTimeCorrADA+nCorrelation, expert, !image, !saveCorr); iHisto++; nCorrelation++;
		}
	}
  nCorrelation = 0;
  for(Int_t i=0;i<8;i++){
  	for(Int_t j=7;j>i;j--){
  		h2d = new TH2F(Form("TimeCorr/H2D_kNTimeCorrADC_%d_%d",i,j),Form("Time Correlation ADC module%d - module%d",i,j),kNPairTimeCorrBins,kPairTimeCorrMin,kPairTimeCorrMax,kNPairTimeCorrBins,kPairTimeCorrMin,kPairTimeCorrMax);
		Add2RawsList(h2d,kNTimeCorrADC+nCorrelation, expert, !image, !saveCorr); iHisto++; nCorrelation++;
		}
	}
  nCorrelation = 0;
  for(Int_t i=0;i<8;i++){
  	for(Int_t j=7;j>i;j--){
  		h1d = new TH1F(Form("TimeDiff/H1D_kNTimeDiffADA_%d_%d",i,j),Form("Time Difference ADA module%d - module%d",i,j),kNPairTimeDiffBins,kPairTimeDiffMin,kPairTimeDiffMax);
		Add2RawsList(h1d,kNTimeDiffADA+nCorrelation, expert, !image, !saveCorr); iHisto++; nCorrelation++;
		}
	}
  nCorrelation = 0;
  for(Int_t i=0;i<8;i++){
  	for(Int_t j=7;j>i;j--){
  		h1d = new TH1F(Form("TimeDiff/H1D_kNTimeDiffADC_%d_%d",i,j),Form("Time Difference ADC module%d - module%d",i,j),kNPairTimeDiffBins,kPairTimeDiffMin,kPairTimeDiffMax);
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
  Double_t chargeTrigADA=0., chargeTrigADC=0.;

  Double_t diffTime=-100000, sumTime = -10000;
  
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
    Float_t adc[16], time[16], width[16], timeCorr[16], adcTrig[16]; 
    Int_t  iPair=0;

    for(Int_t iChannel=0; iChannel<16; iChannel++) { // BEGIN : Loop over channels
		   
      offlineCh = kOfflineChannel[iChannel];
		   
      // Fill Pedestal histograms
      iFlag = 0;   
      for(Int_t j=0; j<=20; j++) {
	if((rawStream->GetBGFlag(iChannel,j) || rawStream->GetBBFlag(iChannel,j))) iFlag++;
      }

      if(iFlag == 0){ //No Flag found
	for(Int_t j=11; j<=17; j++){
	  pedestal= (Int_t) rawStream->GetPedestal(iChannel, j);
	  integrator[offlineCh] = rawStream->GetIntegratorFlag(iChannel, j);
	  OCDBdiff = pedestal - fCalibData->GetPedestal(offlineCh+16*integrator[offlineCh]);

	  FillRawsData((integrator[offlineCh] == 0 ? kPedestalInt0 : kPedestalInt1),offlineCh,pedestal);
	  FillRawsData((integrator[offlineCh] == 0 ? kPedestalDiffInt0 : kPedestalDiffInt1),offlineCh,OCDBdiff);
	}
      }
      // Fill Charge EoI histograms
	   
      adc[offlineCh]    = 0.0;
      adcTrig[offlineCh] = 0.0;
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
      adcTrig[offlineCh] = adcPedSub[10];
      if (imax != -1) {
	Int_t start = imax - fRecoParam->GetNPreClocks();
	if (start < 0) start = 0;
	Int_t end = imax + fRecoParam->GetNPostClocks();
	if (end > 20) end = 20;
	for(Int_t iClock = start; iClock <= end; iClock++) {
	  adc[offlineCh] += adcPedSub[iClock];
	}
	charge = rawStream->GetPedestal(iChannel,imax); // Charge at the maximum
      }
      else charge = -1024;		
      
      if(charge != -1024) FillRawsData(kMaxChargeClock,offlineCh,imax-10);
       
      integrator[offlineCh] = rawStream->GetIntegratorFlag(iChannel,imax);
      
      Int_t board = AliADCalibData::GetBoardNumber(offlineCh);
      time[offlineCh] = rawStream->GetTime(iChannel)*fCalibData->GetTimeResolution(board);
      width[offlineCh] = rawStream->GetWidth(iChannel)*fCalibData->GetWidthResolution(board);

      FillRawsData(kChargeEoI,offlineCh,adc[offlineCh]);

      FillRawsData((integrator[offlineCh] == 0 ? kChargeEoIInt0 : kChargeEoIInt1),offlineCh,charge);
      FillRawsData(kChargeSaturation,offlineCh,charge);

      Float_t sigma = fCalibData->GetSigma(offlineCh+16*integrator[offlineCh]);
		  
      if((adc[offlineCh] > fRecoParam->GetNSigmaPed()*sigma) && !(time[offlineCh] <1.e-6)){ 
	if(offlineCh<8) {
	  mulADC++;
	  chargeADC += adc[offlineCh];
	  chargeTrigADC += adcTrig[offlineCh];
	  
	} else {
	  mulADA++;
	  chargeADA += adc[offlineCh];
	  chargeTrigADA += adcTrig[offlineCh];
	}
      }
      // Fill HPTDC Time Histograms
      timeCorr[offlineCh] = CorrectLeadingTime(offlineCh,time[offlineCh],adc[offlineCh]);
   
      if(time[offlineCh] > 1.e-6){
        FillRawsData(kWidthSlewing,width[offlineCh],adc[offlineCh]);
	Float_t timeErr = 1;
	if (adc[offlineCh]>1) timeErr = 1/adc[offlineCh];

	if (offlineCh<8) {
	    itimeADC++;
	    timeADC += time[offlineCh]/(timeErr*timeErr);
	    weightADC += 1./(timeErr*timeErr); 
      	    if (adc[offlineCh]>1) FillRawsData(kTimeSlewingADC,TMath::Log10(1.0/adc[offlineCh]),time[offlineCh]);
	  }
	else{
	    itimeADA++;
	    timeADA += time[offlineCh]/(timeErr*timeErr);
	    weightADA += 1./(timeErr*timeErr);
	    if (adc[offlineCh]>1) FillRawsData(kTimeSlewingADA,TMath::Log10(1.0/adc[offlineCh]),time[offlineCh]);
	  }
	
      }
      
      // Fill Flag and Charge Versus LHC-Clock histograms
      Int_t nbbFlag = 0;
      Int_t nbgFlag = 0;
      flagBB[offlineCh] = rawStream->GetBBFlag(iChannel,10);
      flagBG[offlineCh] = rawStream->GetBGFlag(iChannel,10);
      
      for(Int_t iEvent=0; iEvent<21; iEvent++){
	charge = rawStream->GetPedestal(iChannel,iEvent);
	Int_t intgr = rawStream->GetIntegratorFlag(iChannel,iEvent);
	Bool_t bbFlag	  = rawStream->GetBBFlag(iChannel,iEvent);
	Bool_t bgFlag	  = rawStream->GetBGFlag(iChannel,iEvent);
	if(bbFlag) nbbFlag++;
	if(bgFlag) nbgFlag++;
	
	FillRawsData((intgr == 0 ? kChargeVsClockInt0 : kChargeVsClockInt1 ), offlineCh,(float)iEvent-10,(float)charge);
	FillRawsData(kBBFlagVsClock, offlineCh,(float)iEvent-10,(float)bbFlag);
	FillRawsData(kBGFlagVsClock, offlineCh,(float)iEvent-10,(float)bgFlag);
	
      }
      FillRawsData(kBBFlagsPerChannel, offlineCh,nbbFlag);
      FillRawsData(kBGFlagsPerChannel, offlineCh,nbgFlag);
      //if((nbbFlag+nbgFlag)>0 && time[offlineCh]<1e-6)FillRawsData(kFlagNoTime,offlineCh);
      //if((nbbFlag+nbgFlag)==0 && time[offlineCh]>1e-6)FillRawsData(kTimeNoFlag,offlineCh);
      
      if((flagBB[offlineCh] || flagBG[offlineCh]) && time[offlineCh]<1e-6)FillRawsData(kFlagNoTime,offlineCh);
      if((!flagBB[offlineCh] && !flagBG[offlineCh]) && time[offlineCh]>1e-6)FillRawsData(kTimeNoFlag,offlineCh);
      
      FillRawsData(kHPTDCTime,offlineCh,time[offlineCh]);
      FillRawsData(kHPTDCTimeRebin,offlineCh,time[offlineCh]);
      FillRawsData(kWidth,offlineCh,width[offlineCh]);
      if(flagBB[offlineCh]) {
      //if(nbbFlag > 0){
	FillRawsData(kHPTDCTimeBB,offlineCh,time[offlineCh]);
	FillRawsData(kHPTDCTimeRebinBB,offlineCh,time[offlineCh]);
	FillRawsData(kWidthBB,offlineCh,width[offlineCh]);
	FillRawsData(kChargeEoIBB,offlineCh,adc[offlineCh]);
	FillRawsData(kNEventsBBFlag,offlineCh);
      }
      if(flagBG[offlineCh]) {
      //if(nbgFlag > 0){
	FillRawsData(kHPTDCTimeBG,offlineCh,time[offlineCh]);
	FillRawsData(kHPTDCTimeRebinBG,offlineCh,time[offlineCh]);
	FillRawsData(kWidthBG,offlineCh,width[offlineCh]);
	FillRawsData(kChargeEoIBG,offlineCh,adc[offlineCh]);
	FillRawsData(kNEventsBGFlag,offlineCh);
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
		
    for(Int_t iChannel=0; iChannel<4; iChannel++) {//Loop over pairs of pads
    	//Enable time is used to turn off the coincidence 
    	if((!fCalibData->GetEnableTiming(iChannel) || flagBB[iChannel]) && (!fCalibData->GetEnableTiming(iChannel+4) || flagBB[iChannel+4])) pBBmulADC++;
	if((!fCalibData->GetEnableTiming(iChannel) || flagBG[iChannel]) && (!fCalibData->GetEnableTiming(iChannel+4) || flagBG[iChannel+4])) pBGmulADC++;

	if((!fCalibData->GetEnableTiming(iChannel+8) || flagBB[iChannel+8]) && (!fCalibData->GetEnableTiming(iChannel+12) || flagBB[iChannel+12])) pBBmulADA++;
	if((!fCalibData->GetEnableTiming(iChannel+8) || flagBG[iChannel+8]) && (!fCalibData->GetEnableTiming(iChannel+12) || flagBG[iChannel+12])) pBGmulADA++;
	}
					
    FillRawsData(kNBBCoincADA,pBBmulADA);
    FillRawsData(kNBBCoincADC,pBBmulADC);
    FillRawsData(kNBGCoincADA,pBGmulADA);
    FillRawsData(kNBGCoincADC,pBGmulADC);
    FillRawsData(kNBBCoincCorr,pBBmulADA,pBBmulADC);
    FillRawsData(kNBGCoincCorr,pBGmulADA,pBGmulADC);
    
    for(Int_t iChannel=0; iChannel<16; iChannel++) {
    	for(Int_t iEvent=0; iEvent<21; iEvent++){
		offlineCh = kOfflineChannel[iChannel];
		Bool_t bbFlag = rawStream->GetBBFlag(iChannel,iEvent);
		if(pBBmulADA>0 || pBBmulADC>0)FillRawsData(kBBFlagVsClock_ADOR, offlineCh,(float)iEvent-10,(float)bbFlag);
		}
	}
    
    //Triggers
    Bool_t UBA = kFALSE;
    Bool_t UBC = kFALSE;
    Bool_t UGA = kFALSE;
    Bool_t UGC = kFALSE;
    
    if(pBBmulADA>=fCalibData->GetBBAThreshold()) UBA = kTRUE;
    if(pBBmulADC>=fCalibData->GetBBCThreshold()) UBC = kTRUE;
    if(pBGmulADA>=fCalibData->GetBGAThreshold()) UGA = kTRUE;
    if(pBGmulADC>=fCalibData->GetBGCThreshold()) UGC = kTRUE;

    if(UBA) FillRawsData(kTriggers,0);
    if(UBC) FillRawsData(kTriggers,1);
    if(UGA) FillRawsData(kTriggers,2);
    if(UGC) FillRawsData(kTriggers,3);
    if(UBA && UBC) FillRawsData(kTriggers,4);
    if(UBA || UBC) FillRawsData(kTriggers,5);
    if((UBA || UBC) && !(UGA || UGC)) FillRawsData(kTriggers,6);
    if(UGA && UBC) FillRawsData(kTriggers,7);
    if(UGC && UBA) FillRawsData(kTriggers,8);
    if(UGA || UGC) FillRawsData(kTriggers,9);
    if((UGA && UBC) || (UGC && UBA)) FillRawsData(kTriggers,10);
    	
    //Average times
    if(weightADA>1) timeADA /= weightADA; 
    else timeADA = -1024.;
    if(weightADC>1) timeADC /= weightADC;
    else timeADC = -1024.;
    if(timeADA<1.e-6 || timeADC<1.e-6) {
    	diffTime = -1024.;
	sumTime = -1024;
	}
    else {
    	diffTime = timeADA - timeADC;
	sumTime = timeADA + timeADC;
	}
    	
    FillRawsData(kMeanTimeADA,timeADA);
    FillRawsData(kMeanTimeADC,timeADC);
    FillRawsData(kMeanTimeDiff,diffTime);
    FillRawsData(kMeanTimeCorr,timeADA,timeADC);
    FillRawsData(kMeanTimeSumDiff,diffTime,sumTime);

    FillRawsData(kMultiADA,mulADA);
    FillRawsData(kMultiADC,mulADC);

    FillRawsData(kChargeADA,chargeADA);
    FillRawsData(kChargeADC,chargeADC);
    FillRawsData(kChargeAD,chargeADA + chargeADC);
    if(mulADA!=0)FillRawsData(kChargeADA_PC,chargeTrigADA/mulADA);
    if(mulADC!=0)FillRawsData(kChargeADC_PC,chargeTrigADA/mulADC);
    
    //Decisions
    Int_t windowOffset = (fCalibData->GetTriggerCountOffset(0) - 3242)*25;
    Int_t ADADecision=0;
    Int_t ADCDecision=0; 

    if(timeADA > (fADADist - windowOffset + 5*fRecoParam->GetTimeWindowBBALow()) && timeADA < (fADADist - windowOffset  + 5*fRecoParam->GetTimeWindowBBAUp())) ADADecision=1;
    else if(timeADA > (-fADADist - windowOffset + 5*fRecoParam->GetTimeWindowBGALow()) && timeADA < (-fADADist - windowOffset + 5*fRecoParam->GetTimeWindowBGAUp())) ADADecision=2;
    else if(timeADA>-1024.+1.e-6) ADADecision=3;
    else ADADecision=0;
    
    if(timeADC > (fADCDist - windowOffset + 5*fRecoParam->GetTimeWindowBBCLow()) && timeADC < (fADCDist - windowOffset + 5*fRecoParam->GetTimeWindowBBCUp())) ADCDecision=1;
    else if(timeADC > (-fADCDist - windowOffset + 5*fRecoParam->GetTimeWindowBGCLow()) && timeADC < (-fADCDist - windowOffset + 5*fRecoParam->GetTimeWindowBGCUp())) ADCDecision=2;
    else if(timeADC>-1024.+1.e-6) ADCDecision=3;
    else ADCDecision=0;

    FillRawsData(kDecisions,ADADecision,ADCDecision);
    
    break;
  } // END of SWITCH : EVENT TYPE 
	 
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
  // In case of pathological signals
  if (adc < 1)return time;

  // Slewing correction
  //time -= fTimeSlewing->Eval(adc);
  
   // Channel alignment and general offset subtraction
  //  time -= fHptdcOffset[i];
  //AliInfo(Form("time-offset %f", time));
  

  return time;
}

