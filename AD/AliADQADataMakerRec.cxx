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
//#include "AliADTrending.h"
#include "AliADCalibData.h"
#include "AliCTPTimeParams.h"
#include "event.h"

ClassImp(AliADQADataMakerRec)
           
//____________________________________________________________________________ 
AliADQADataMakerRec::AliADQADataMakerRec() : 
AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kAD), "AD Quality Assurance Data Maker"),
  fCalibData(0x0),
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
 
  AliCDBEntry *entry = AliCDBManager::Instance()->Get("GRP/CTP/CTPtiming");
  if (!entry) AliFatal("CTP timing parameters are not found in OCDB !");
  AliCTPTimeParams *ctpParams = (AliCTPTimeParams*)entry->GetObject();
  Float_t l1Delay = (Float_t)ctpParams->GetDelayL1L0()*25.0;

  AliCDBEntry *entry1 = AliCDBManager::Instance()->Get("GRP/CTP/TimeAlign");
  if (!entry1) AliFatal("CTP time-alignment is not found in OCDB !");
  AliCTPTimeParams *ctpTimeAlign = (AliCTPTimeParams*)entry1->GetObject();
  l1Delay += ((Float_t)ctpTimeAlign->GetDelayL1L0()*25.0);

  AliCDBEntry *entry2 = AliCDBManager::Instance()->Get("AD/Calib/TimeDelays");
  if (!entry2) AliFatal("AD time delays are not found in OCDB !");
  TH1F *delays = (TH1F*)entry2->GetObject();

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
		      delays->GetBinContent(i+1)//+
		      //      kV0Offset
		      );
    //		      AliInfo(Form(" fTimeOffset[%d] = %f  kV0offset %f",i,fTimeOffset[i],kV0Offset));
  }

 
 
  
 	
  TTimeStamp currentTime;
  fCycleStartTime = currentTime.GetSec();
 
  //  fNTotEvents = 0;
}
//____________________________________________________________________________ 
void AliADQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  // Detector specific actions at end of cycle
  // Does the QA checking
  ResetEventTrigClasses();
  //
  AliQAChecker::Instance()->Run(AliQAv1::kAD, task, list) ;
  
  if(task == AliQAv1::kRAWS){
    TTimeStamp currentTime;
    fCycleStopTime = currentTime.GetSec();
  }

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    if (! IsValidEventSpecie(specie, list)) continue ;
    SetEventSpecie(AliRecoParam::ConvertIndex(specie));
    if(task == AliQAv1::kRAWS) {
    } else if (task == AliQAv1::kESDS) {
    }
  }
}

//____________________________________________________________________________ 
void AliADQADataMakerRec::InitESDs()
{
  // Creates histograms to control ESDs
}

//____________________________________________________________________________ 
void AliADQADataMakerRec::InitDigits()
{

}

//____________________________________________________________________________
void AliADQADataMakerRec::MakeDigits()
{

}

//____________________________________________________________________________
void AliADQADataMakerRec::MakeDigits(TTree *digitTree)
{

}


//____________________________________________________________________________
void AliADQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
 
}

//____________________________________________________________________________ 
void AliADQADataMakerRec::InitRaws()
{
  // Creates RAW histograms in Raws subdir

  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  const Int_t kNintegrator  =    2;
 
  const Int_t kNTdcTimeBins  = 1280;
  const Float_t kTdcTimeMin    =    0.;
  const Float_t kTdcTimeMax    = 75.;
    const Int_t kNTdcWidthBins =  256;
  const Float_t kTdcWidthMin   =    0;
  const Float_t kTdcWidthMax   =  200.;
  const Int_t kNChargeBins   = 1024;
  const Float_t kChargeMin     =    0;
  const Float_t kChargeMax     = 1024;
  const Int_t kNChannelBins  =   16;
  const Float_t kChannelMin    =    0;
  const Float_t kChannelMax    =   16;
  const Int_t kNPedestalBins =  200;
  const Float_t kPedestalMin   =    0;
  const Float_t kPedestalMax   =  200;

  TH2I * h2i;
  TH2F * h2d;
  TH1I * h1i;
  TH1F * h1d;

  int iHisto =0;

  // Creation of Cell Multiplicity Histograms
  h1i = new TH1I("H1I_Multiplicity_ADA", "Cell Multiplicity in ADA;# of Cells;Entries", 35, 0, 35) ;  
  Add2RawsList(h1i,kMultiADA, !expert, image, saveCorr);   iHisto++;
  h1i = new TH1I("H1I_Multiplicity_ADC", "Cell Multiplicity in ADC;# of Cells;Entries", 35, 0, 35) ;  
  Add2RawsList(h1i,kMultiADC, !expert, image, saveCorr);   iHisto++;
 
  // Creation of Total Charge Histograms
  h1d = new TH1F("H1D_Charge_ADA", "Total Charge in ADA;Charge [ADC counts];Counts", 4000, 0, 30000) ;  
  Add2RawsList(h1d,kChargeADA, !expert, image, saveCorr);   iHisto++;
  h1d = new TH1F("H1D_Charge_ADC", "Total Charge in ADC;Charge [ADC counts];Counts", 4000, 0, 50000) ;  
  Add2RawsList(h1d,kChargeADC, !expert, image, saveCorr);   iHisto++;
  h1d = new TH1F("H1D_Charge_AD", "Total Charge in AD;Charge [ADC counts];Counts", 4000, 0, 80000) ;  
  Add2RawsList(h1d,kChargeAD, !expert,  image, saveCorr);   iHisto++;
   

  // Creation of Charge EoI histogram 
  h2d = new TH2F("H2D_ChargeEoI", "Charge Event of Interest;Channel Number;Charge [ADC counts]"
		 ,kNChannelBins, kChannelMin, kChannelMax, kNChargeBins, kChargeMin, 2.*kChargeMax);
  Add2RawsList(h2d,kChargeEoI, !expert, image, saveCorr); iHisto++;

  for(Int_t iInt=0;iInt<kNintegrator;iInt++){
    // Creation of Pedestal histograms 
    h2i = new TH2I(Form("H2I_Pedestal_Int%d",iInt), Form("Pedestal (Int%d);Channel;Pedestal [ADC counts]",iInt)
		   ,kNChannelBins, kChannelMin, kChannelMax,kNPedestalBins,kPedestalMin ,kPedestalMax );
    Add2RawsList(h2i,(iInt == 0 ? kPedestalInt0 : kPedestalInt1), !expert, image, saveCorr); iHisto++;
	

    // Creation of Charge EoI histograms 
    h2i = new TH2I(Form("H2I_ChargeEoI_Int%d",iInt), Form("Charge EoI (Int%d);Channel;Charge [ADC counts]",iInt)
		   ,kNChannelBins, kChannelMin, kChannelMax, kNChargeBins, kChargeMin, kChargeMax);
    Add2RawsList(h2i,(iInt == 0 ? kChargeEoIInt0 : kChargeEoIInt1), !expert, image, saveCorr); iHisto++;
    
  }	
  
  // Creation of Time histograms 
  h2i = new TH2I("H2I_Width", "HPTDC Width;Channel;Width [ns]",kNChannelBins, kChannelMin, kChannelMax, kNTdcWidthBins, kTdcWidthMin, kTdcWidthMax);
  Add2RawsList(h2i,kWidth, !expert, image, saveCorr); iHisto++;

  h2i = new TH2I("H2I_HPTDCTime", "HPTDC Time;Channel;Leading Time [ns]",kNChannelBins, kChannelMin, kChannelMax, kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax);
  Add2RawsList(h2i,kHPTDCTime, !expert, image, saveCorr); iHisto++;
	
  h1d = new TH1F("H1D_ADA_Time", "ADA Time;Time [ns];Counts",kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax);
  Add2RawsList(h1d,kADATime, !expert, image, saveCorr); iHisto++;
	
  h1d = new TH1F("H1D_ADC_Time", "ADC Time;Time [ns];Counts",kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax);
  Add2RawsList(h1d,kADCTime, !expert, image, saveCorr); iHisto++;
	
  h1d = new TH1F("H1D_Diff_Time","Diff ADA-ADC Time;Time [ns];Counts",kNTdcTimeBins, -50., 50.);
  Add2RawsList(h1d,kDiffTime, !expert, image, saveCorr); iHisto++;

  h2d = new TH2F("H2D_TimeADA_ADC", "Mean Time in ADC versus ADA;Time ADA [ns];Time ADC [ns]", kNTdcTimeBins/8, kTdcTimeMin,kTdcTimeMax,kNTdcTimeBins/8, kTdcTimeMin,kTdcTimeMax) ;  
  Add2RawsList(h2d,kTimeADAADC, !expert, image, saveCorr);   iHisto++;
  

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

  
  switch (eventType){
  case PHYSICS_EVENT:
  
    //    fNTotEvents++; // Use framework counters instead

    Int_t  iFlag=0;
    Int_t  pedestal;
    Int_t  integrator;
    Bool_t flagBB[16];	 
    Bool_t flagBG[16];	 
    Float_t charge;
    Int_t  offlineCh;
    Float_t adc[16], time[16], width[16], timeCorr[16]; 

    for(Int_t iChannel=0; iChannel<16; iChannel++) { // BEGIN : Loop over channels
		   
      offlineCh = rawStream->GetOfflineChannel(iChannel);
		   
      // Fill Pedestal histograms
	   
      for(Int_t j=15; j<21; j++) {
	if((rawStream->GetBGFlag(iChannel,j) || rawStream->GetBBFlag(iChannel,j))) iFlag++;
      }

      if(iFlag == 0){ //No Flag found
	for(Int_t j=15; j<21; j++){
	  pedestal= (Int_t) rawStream->GetPedestal(iChannel, j);
	  integrator = rawStream->GetIntegratorFlag(iChannel, j);

	  FillRawsData((integrator == 0 ? kPedestalInt0 : kPedestalInt1),offlineCh,pedestal);
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
	Int_t k = offlineCh+64*iIntegrator;

	//printf(Form("clock = %d adc = %f ped %f\n",iClock,rawStream->GetPedestal(iChannel,iClock),fPedestal[k]));

	adcPedSub[iClock] = rawStream->GetPedestal(iChannel,iClock) - fCalibData->GetPedestal(k);
	//				if(adcPedSub[iClock] <= GetRecoParam()->GetNSigmaPed()*fCalibData->GetSigma(k)) {
	if(adcPedSub[iClock] <= 2.*fCalibData->GetSigma(k)) {
	  adcPedSub[iClock] = 0;
	  continue;
	}
	//		   		if(iClock < GetRecoParam()->GetStartClock() || iClock > GetRecoParam()->GetEndClock()) continue;
	if(iClock < 8 || iClock > 12) continue;
	if(adcPedSub[iClock] > maxadc) {
	  maxadc = adcPedSub[iClock];
	  imax   = iClock;
	}
      }
      //printf(Form("Channel %d (online), %d (offline)\n",iChannel,j)); 
      if (imax != -1) {
	//		   		Int_t start = imax - GetRecoParam()->GetNPreClocks();
	Int_t start = imax - 2;
	if (start < 0) start = 0;
	//		   		Int_t end = imax + GetRecoParam()->GetNPostClocks();
	Int_t end = imax + 1;
	if (end > 20) end = 20;
	for(Int_t iClock = start; iClock <= end; iClock++) {
	  adc[offlineCh] += adcPedSub[iClock];
	}
      }
	
		
      Int_t iClock  = imax;
      charge = rawStream->GetPedestal(iChannel,iClock); // Charge at the maximum 

      integrator    = rawStream->GetIntegratorFlag(iChannel,iClock);
      flagBB[offlineCh]	 = rawStream->GetBBFlag(iChannel, iClock);
      flagBG[offlineCh]	 = rawStream->GetBGFlag(iChannel,iClock );
      Int_t board = AliADCalibData::GetBoardNumber(offlineCh);
      time[offlineCh] = rawStream->GetTime(iChannel)*fCalibData->GetTimeResolution(board);
      width[offlineCh] = rawStream->GetWidth(iChannel)*fCalibData->GetWidthResolution(board);

      if (time[offlineCh] >= 1e-6) FillRawsData(kChargeEoI,offlineCh,adc[offlineCh]);

      FillRawsData((integrator == 0 ? kChargeEoIInt0 : kChargeEoIInt1),offlineCh,charge);

      Float_t sigma = fCalibData->GetSigma(offlineCh+16*integrator);
		  
      if((adc[offlineCh] > 2.*sigma) && !(time[offlineCh] <1.e-6)){ 
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

      const Float_t p1 = 2.50; // photostatistics term in the time resolution
      const Float_t p2 = 3.00; // sleewing related term in the time resolution
      if(timeCorr[offlineCh]>-1024 + 1.e-6){
	Float_t nphe = adc[offlineCh]*kChargePerADC/(fCalibData->GetGain(offlineCh)*TMath::Qe());
	Float_t timeErr = 0;
	if (nphe>1.e-6) timeErr = TMath::Sqrt(kIntTimeRes*kIntTimeRes+
					      p1*p1/nphe+
					      p2*p2*(fTimeSlewing->GetParameter(0)*fTimeSlewing->GetParameter(1))*(fTimeSlewing->GetParameter(0)*fTimeSlewing->GetParameter(1))*
					      TMath::Power(adc[offlineCh]/fCalibData->GetCalibDiscriThr(offlineCh,kTRUE),2.*(fTimeSlewing->GetParameter(1)-1.))/
					      (fCalibData->GetCalibDiscriThr(offlineCh,kTRUE)*fCalibData->GetCalibDiscriThr(offlineCh,kTRUE)));

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
      FillRawsData(kHPTDCTime,offlineCh,timeCorr[offlineCh]);
      FillRawsData(kWidth,offlineCh,width[offlineCh]);

    }// END of Loop over channels

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
Float_t AliADQADataMakerRec::CorrectLeadingTime(Int_t i, Float_t time, Float_t adc) const
{
  // Correct the leading time
  // for slewing effect and
  // misalignment of the channels
  if (time < 1e-6) return -1024;

  // Channel alignment and general offset subtraction
  //  if (i < 32) time -= kV0CDelayCables;
  //  time -= fTimeOffset[i];
  //AliInfo(Form("time-offset %f", time));

  // In case of pathological signals
  if (adc < 1e-6) return time;

  // Slewing correction
  Float_t thr = fCalibData->GetCalibDiscriThr(i,kTRUE);
  //AliInfo(Form("adc %f thr %f dtime %f ", adc,thr,fTimeSlewing->Eval(adc/thr)));
  time -= fTimeSlewing->Eval(adc/thr);

  return time;
}

