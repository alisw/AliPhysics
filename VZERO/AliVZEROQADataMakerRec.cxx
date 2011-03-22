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
#include <TH2D.h> 
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
#include "AliVZEROQADataMakerRec.h"
#include "AliQAChecker.h"
#include "AliRawReader.h"
#include "AliVZERORawStream.h"
#include "AliVZEROdigit.h"
#include "AliVZEROConst.h"
#include "AliVZEROReconstructor.h"
#include "AliVZEROTrending.h"
#include "AliVZEROCalibData.h"
#include "AliCTPTimeParams.h"
#include "event.h"

 const Float_t kMinBBA = 68. ;
 const Float_t kMaxBBA = 100. ;
 const Float_t kMinBBC = 75.5 ;
 const Float_t kMaxBBC = 100. ;
 const Float_t kMinBGA = 54. ;
 const Float_t kMaxBGA = 58. ;
 const Float_t kMinBGC = 69.5 ;
 const Float_t kMaxBGC = 74. ;

 
 
 

ClassImp(AliVZEROQADataMakerRec)
           
//____________________________________________________________________________ 
  AliVZEROQADataMakerRec::AliVZEROQADataMakerRec() : 
	AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kVZERO), "VZERO Quality Assurance Data Maker"),
	fCalibData(0x0),
        fEvent(0), 
        fNTotEvents(0), 
        fNSubEvents(0), 
        fTrendingUpdateEvent(0), 
        fNTrendingUpdates(0), 
        fTrendingUpdateTime(0), 
        fCycleStartTime(0), 
        fCycleStopTime(0),
		fTimeSlewing(0)
    
{
   // Constructor
   
      AliDebug(AliQAv1::GetQADebugLevel(), "Construct VZERO QA Object");

   for(Int_t i=0; i<64; i++){  
       fEven[i] = 0;   
       fOdd[i]  = 0;
  }
  
   for(Int_t i=0; i<128; i++){  
       fADCmean[i] = 0.0;   }	
}

//____________________________________________________________________________ 
  AliVZEROQADataMakerRec::AliVZEROQADataMakerRec(const AliVZEROQADataMakerRec& qadm) :
  AliQADataMakerRec(),
	fCalibData(0x0),
        fEvent(0), 
        fNTotEvents(0), 
        fNSubEvents(0), 
        fTrendingUpdateEvent(0), 
        fNTrendingUpdates(0), 
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
AliVZEROQADataMakerRec& AliVZEROQADataMakerRec::operator = (const AliVZEROQADataMakerRec& qadm )
{
  // Equal operator
  
  this->~AliVZEROQADataMakerRec();
  new(this) AliVZEROQADataMakerRec(qadm);
  return *this;
}

//____________________________________________________________________________
AliVZEROCalibData* AliVZEROQADataMakerRec::GetCalibData() const

{
   AliCDBManager *man = AliCDBManager::Instance();

   AliCDBEntry *entry=0;

   entry = man->Get("VZERO/Calib/Data",fRun);
   if(!entry){
	AliWarning("Load of calibration data from default storage failed!");
	AliWarning("Calibration data will be loaded from local storage ($ALICE_ROOT)");
	
	man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
	entry = man->Get("VZERO/Calib/Data",fRun);
   }
   // Retrieval of data in directory VZERO/Calib/Data:

   AliVZEROCalibData *calibdata = 0;

   if (entry) calibdata = (AliVZEROCalibData*) entry->GetObject();
   if (!calibdata)  AliFatal("No calibration data from calibration database !");

   return calibdata;
}


 
//____________________________________________________________________________ 
void AliVZEROQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  // Detector specific actions at end of cycle
  // Does the QA checking
  
  AliQAChecker::Instance()->Run(AliQAv1::kVZERO, task, list) ;
  
  if(task == AliQAv1::kRAWS){
	TTimeStamp currentTime;
	fCycleStopTime = currentTime.GetSec();
  }

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    if (! IsValidEventSpecie(specie, list)) 
      continue ;
    SetEventSpecie(AliRecoParam::ConvertIndex(specie)) ; 
    if(task == AliQAv1::kRAWS){
    } else if (task == AliQAv1::kESDS) {
    }
  }
}

//____________________________________________________________________________ 
void AliVZEROQADataMakerRec::InitESDs()
{
  // Creates histograms to control ESDs
  
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
	
  TH2D * h2d;
  TH1I * h1i;
  TH1D * h1d;
		
  h1i = new TH1I("H1I_Cell_Multiplicity_V0A", "Cell Multiplicity in V0A;Multiplicity (Nb of Cell);Counts", 35, 0, 35) ;  
  Add2ESDsList(h1i, kCellMultiV0A, !expert, image)  ;  
                                                                                                        
  h1i = new TH1I("H1I_Cell_Multiplicity_V0C", "Cell Multiplicity in V0;Multiplicity (Nb of Cell);Counts", 35, 0, 35) ;  
  Add2ESDsList(h1i, kCellMultiV0C, !expert, image)  ;  
   
  h1d = new TH1D("H1D_MIP_Multiplicity_V0A", "MIP Multiplicity in V0A;Multiplicity (Nb of MIP);Counts", 1000, 0, 1000) ;  
  Add2ESDsList(h1d, kMIPMultiV0A, !expert, image)  ;  
  
  h1d = new TH1D("H1D_MIP_Multiplicity_V0C", "MIP Multiplicity in V0C;Multiplicity (Nb of MIP);Counts", 1000, 0, 1000) ;  
  Add2ESDsList(h1d, kMIPMultiV0C, !expert, image)  ;  

  h2d = new TH2D("H2D_MIP_Multiplicity_Channel", "MIP Multiplicity per Channel;Channel;Multiplicity (Nb of MIP)",64, 0, 64, 100, 0, 100) ;  
  Add2ESDsList(h2d, kMIPMultiChannel, !expert, image)  ;  
  
  h1d = new TH1D("H1D_BBFlag_Counters", "BB Flag Counters;Channel;Counts",64, 0, 64) ;  
  Add2ESDsList(h1d, kBBFlag, !expert, image)  ;  
  
  h1d = new TH1D("H1D_BGFlag_Counters", "BG Flag Counters;Channel;Counts",64, 0, 64) ;  
  Add2ESDsList(h1d, kBGFlag, !expert, image)  ;  
  
  h2d = new TH2D("H2D_Charge_Channel", "ADC Charge per channel;Channel;Charge (ADC counts)",64, 0, 64, 1024, 0, 1024) ;  
  Add2ESDsList(h2d, kChargeChannel, !expert, image)  ;  
  
  h2d = new TH2D("H2D_Time_Channel", "Time per channel;Channel;Time (ns)",64, 0, 64, 400, -100, 100) ;  
  Add2ESDsList(h2d, kTimeChannel, !expert, image)  ;  
  
  h1d = new TH1D("H1D_V0A_Time", "Mean V0A Time;Time (ns);Counts",1000, -100., 100.);
  Add2ESDsList(h1d,kESDV0ATime, !expert, image); 
  
  h1d = new TH1D("H1D_V0C_Time", "Mean V0C Time;Time (ns);Counts",1000, -100., 100.);
  Add2ESDsList(h1d,kESDV0CTime, !expert, image); 
  
  h1d = new TH1D("H1D_Diff_Time", "Diff Time V0A - V0C;Diff Time V0A - V0C (ns);Counts",1000, -200., 200.);
  Add2ESDsList(h1d,kESDDiffTime, !expert, image); 
	
}

//____________________________________________________________________________ 
 void AliVZEROQADataMakerRec::InitRaws()
 {
   // Creates RAW histograms in Raws subdir

   const Bool_t expert   = kTRUE ; 
   const Bool_t saveCorr = kTRUE ; 
   const Bool_t image    = kTRUE ; 

  const Int_t kNintegrator  =    2;
 
  const Int_t kNTdcTimeBins  = 1280;
  const Float_t kTdcTimeMin    =    0.;
  const Float_t kTdcTimeMax    = 125.;
  const Int_t kNTdcWidthBins =  128;
  const Float_t kTdcWidthMin   =    0;
  const Float_t kTdcWidthMax   =  50.;
  const Int_t kNChargeBins   = 1024;
  const Float_t kChargeMin     =    0;
  const Float_t kChargeMax     = 1024;
  const Int_t kNChannelBins  =   64;
  const Float_t kChannelMin    =    0;
  const Float_t kChannelMax    =   64;
  const Int_t kNPedestalBins =  200;
  const Float_t kPedestalMin   =    0;
  const Float_t kPedestalMax   =  200;
  const Float_t kTimeMin       =   0;
  const Float_t kTimeMax       = 100;
  const Int_t kNMIPBins      = 512;
  const Float_t kMIPMin        =   0;
  const Float_t kMIPMax        = 16;

  TH2I * h2i;
  TH2D * h2d;
  TH1I * h1i;
  TH1D * h1d;

  int iHisto =0;
    // Creation of Trigger Histogram
  h1d = new TH1D("H1D_Trigger_Type", "V0 Trigger Type;;Counts", 4,0 ,4) ;  
  Add2RawsList(h1d,kTriggers, !expert, image, saveCorr);   iHisto++;
 	h1d->SetFillColor(29);
	h1d->SetLineWidth(2);
	h1d->GetXaxis()->SetLabelSize(0.06);
    h1d->GetXaxis()->SetNdivisions(808,kFALSE);
	h1d->GetXaxis()->SetBinLabel(1, "V0-AND");
	h1d->GetXaxis()->SetBinLabel(2, "V0-OR");
   	h1d->GetXaxis()->SetBinLabel(3, "V0-BGA");
	h1d->GetXaxis()->SetBinLabel(4, "V0-BGC");

   // Creation of Cell Multiplicity Histograms
  h1i = new TH1I("H1I_Multiplicity_V0A", "Cell Multiplicity in V0A;# of Cells;Entries", 35, 0, 35) ;  
  Add2RawsList(h1i,kMultiV0A, expert, image, saveCorr);   iHisto++;
  h1i = new TH1I("H1I_Multiplicity_V0C", "Cell Multiplicity in V0C;# of Cells;Entries", 35, 0, 35) ;  
  Add2RawsList(h1i,kMultiV0C, expert, image, saveCorr);   iHisto++;
 
  // Creation of Total Charge Histograms
  h1d = new TH1D("H1D_Charge_V0A", "Total Charge in V0A;Charge [ADC counts];Counts", 2000, 0, 10000) ;  
  Add2RawsList(h1d,kChargeV0A, expert, !image, saveCorr);   iHisto++;
  h1d = new TH1D("H1D_Charge_V0C", "Total Charge in V0C;Charge [ADC counts];Counts", 2000, 0, 10000) ;  
  Add2RawsList(h1d,kChargeV0C, expert, !image, saveCorr);   iHisto++;
  h1d = new TH1D("H1D_Charge_V0", "Total Charge in V0;Charge [ADC counts];Counts", 2000, 0, 20000) ;  
  Add2RawsList(h1d,kChargeV0, expert, !image, saveCorr);   iHisto++;
  
  // Creation of MIP Histograms
  h1d = new TH1D("H1D_MIP_V0A", "Total MIP in V0A;Multiplicity [MIP];Counts", kNMIPBins,kMIPMin ,32*kMIPMax) ;  
  Add2RawsList(h1d,kRawMIPV0A, expert, !image, saveCorr);   iHisto++;
  h1d = new TH1D("H1D_MIP_V0C", "Total MIP in V0C;Multiplicity [MIP];Counts", kNMIPBins,kMIPMin ,32*kMIPMax) ;  
  Add2RawsList(h1d,kRawMIPV0C, expert, !image, saveCorr);   iHisto++;
  h1d = new TH1D("H1D_MIP_V0", "Total MIP in V0;Multiplicity [MIP];Counts", 2*kNMIPBins,kMIPMin ,64*kMIPMax) ;  
  Add2RawsList(h1d,kRawMIPV0, expert, !image, saveCorr);   iHisto++;
  h2d = new TH2D("H2D_MIP_Channel", "Nb of MIP per channel;Channel;# of Mips", kNChannelBins, kChannelMin, kChannelMax,kNMIPBins,kMIPMin ,kMIPMax) ;  
  Add2RawsList(h2d,kRawMIPChannel, expert, !image, !saveCorr);   iHisto++;
  
 

  // Creation of Charge EoI histogram 
   h2d = new TH2D("H2D_ChargeEoI", "Charge Event of Interest;Channel Number;Charge [ADC counts]"
   		,kNChannelBins, kChannelMin, kChannelMax, kNChargeBins, kChargeMin, kChargeMax);
   Add2RawsList(h2d,kChargeEoI, !expert, image, !saveCorr); iHisto++;

 for(Int_t iInt=0;iInt<kNintegrator;iInt++){
    // Creation of Pedestal histograms 
    h2i = new TH2I(Form("H2I_Pedestal_Int%d",iInt), Form("Pedestal (Int%d);Channel;Pedestal [ADC counts]",iInt)
		,kNChannelBins, kChannelMin, kChannelMax,kNPedestalBins,kPedestalMin ,kPedestalMax );
    Add2RawsList(h2i,(iInt == 0 ? kPedestalInt0 : kPedestalInt1), expert, !image, !saveCorr); iHisto++;
	

   // Creation of Charge EoI histograms 
    h2i = new TH2I(Form("H2I_ChargeEoI_Int%d",iInt), Form("Charge EoI (Int%d);Channel;Charge [ADC counts]",iInt)
		,kNChannelBins, kChannelMin, kChannelMax, kNChargeBins, kChargeMin, kChargeMax);
    Add2RawsList(h2i,(iInt == 0 ? kChargeEoIInt0 : kChargeEoIInt1), expert, image, !saveCorr); iHisto++;
    
    h2i = new TH2I(Form("H2I_ChargeEoI_BB_Int%d",iInt), Form("Charge EoI w/ BB Flag (Int%d);Channel;Charge [ADC counts]",iInt)
		,kNChannelBins, kChannelMin, kChannelMax, kNChargeBins, kChargeMin, kChargeMax);
    Add2RawsList(h2i,(iInt == 0 ? kChargeEoIBBInt0 : kChargeEoIBBInt1), expert, !image, !saveCorr); iHisto++;
    
    h2i = new TH2I(Form("H2I_ChargeEoI_BG_Int%d",iInt), Form("Charge EoI w/ BG Flag (Int%d);Channel;Charge [ADC counts]",iInt)
		,kNChannelBins, kChannelMin, kChannelMax, kNChargeBins, kChargeMin, kChargeMax);
    Add2RawsList(h2i,(iInt == 0 ?  kChargeEoIBGInt0: kChargeEoIBGInt1), expert, !image, !saveCorr); iHisto++;

    // Creation of Charge versus LHC Clock histograms 
    h2d = new TH2D(Form("H2D_ChargeVsClock_Int%d",iInt), Form("Charge Versus LHC-Clock (Int%d);Channel;LHCClock;Charge [ADC counts]",iInt)
		,kNChannelBins, kChannelMin, kChannelMax,21, -10.5, 10.5 );
    Add2RawsList(h2d,(iInt == 0 ? kChargeVsClockInt0 : kChargeVsClockInt1 ), expert, !image, !saveCorr); iHisto++;
	
    // Creation of Minimum Bias Charge histograms 
    for(Int_t iBB=0;iBB<2;iBB++){
		for(Int_t iBG=0;iBG<2;iBG++){
			h2i = new TH2I(Form("H2I_ChargeMB_BB%d_BG%d_Int%d",iBB,iBG,iInt), Form("MB Charge (BB=%d, BG=%d, Int=%d);Channel;Charge [ADC counts]",iBB,iBG,iInt)
				,kNChannelBins, kChannelMin, kChannelMax,kNChargeBins, kChargeMin, kChargeMax);
			int idx;
			if(iInt==0){
				if(iBB==0){
					if(iBG==0) idx = kChargeMBBB0BG0Int0;
					else idx = kChargeMBBB0BG1Int0;
				} else {
					if(iBG==0) idx = kChargeMBBB1BG0Int0;
					else idx = kChargeMBBB1BG1Int0;
				}
			} else {
				if(iBB==0){
					if(iBG==0) idx = kChargeMBBB0BG0Int1;
					else idx = kChargeMBBB0BG1Int1;
				} else {
					if(iBG==0) idx = kChargeMBBB1BG0Int1;
					else idx = kChargeMBBB1BG1Int1;
				}
			}
			Add2RawsList(h2i,idx, expert, !image, !saveCorr); iHisto++;
		}
    }
	
 }
 
     // Creation of Time histograms 
	h2i = new TH2I("H2I_Width", "HPTDC Width;Channel;Width [ns]",kNChannelBins, kChannelMin, kChannelMax, kNTdcWidthBins, kTdcWidthMin, kTdcWidthMax);
 	Add2RawsList(h2i,kWidth, expert, !image, !saveCorr); iHisto++;

 	h2i = new TH2I("H2I_Width_BB", "HPTDC Width w/ BB Flag condition;Channel;Width [ns]",kNChannelBins, kChannelMin, kChannelMax, kNTdcWidthBins, kTdcWidthMin, kTdcWidthMax);
 	Add2RawsList(h2i,kWidthBB, expert, !image, !saveCorr); iHisto++;

 	h2i = new TH2I("H2I_Width_BG", "HPTDC Width w/ BG Flag condition;Channel;Width [ns]",kNChannelBins, kChannelMin, kChannelMax, kNTdcWidthBins, kTdcWidthMin, kTdcWidthMax);
 	Add2RawsList(h2i,kWidthBG, expert, !image, !saveCorr); iHisto++;

 	h2i = new TH2I("H2I_HPTDCTime", "HPTDC Time;Channel;Leading Time [ns]",kNChannelBins, kChannelMin, kChannelMax, kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax);
 	Add2RawsList(h2i,kHPTDCTime, expert, image, !saveCorr); iHisto++;

 	h2i = new TH2I("H2I_HPTDCTime_BB", "HPTDC Time w/ BB Flag condition;Channel;Leading Time [ns]",kNChannelBins, kChannelMin, kChannelMax, kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax);
 	Add2RawsList(h2i,kHPTDCTimeBB, !expert, image, !saveCorr); iHisto++;

 	h2i = new TH2I("H2I_HPTDCTime_BG", "HPTDC Time w/ BG Flag condition;Channel;Leading Time [ns]",kNChannelBins, kChannelMin, kChannelMax, kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax);
 	Add2RawsList(h2i,kHPTDCTimeBG, !expert, image, !saveCorr); iHisto++;
	
 	h1d = new TH1D("H1D_V0A_Time", "V0A Time;Time [ns];Counts",kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax);
 	Add2RawsList(h1d,kV0ATime, expert, !image, saveCorr); iHisto++;
	
 	h1d = new TH1D("H1D_V0C_Time", "V0C Time;Time [ns];Counts",kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax);
 	Add2RawsList(h1d,kV0CTime, expert, !image, saveCorr); iHisto++;
	
 	h1d = new TH1D("H1D_Diff_Time","Diff V0A-V0C Time;Time [ns];Counts",2*kNTdcTimeBins, -50., 50.);
 	Add2RawsList(h1d,kDiffTime, expert, !image, saveCorr); iHisto++;

    h2d = new TH2D("H2D_TimeV0A_V0C", "Mean Time in V0C versus V0A;Time V0A [ns];Time V0C [ns]", 
  		150, kTimeMin,kTimeMax,150,kTimeMin,kTimeMax) ;  
    Add2RawsList(h2d,kTimeV0AV0C, !expert, image, !saveCorr);   iHisto++;
	
 	// Creation of Flag versus LHC Clock histograms 

 	h1d = new TH1D("H1D_BBFlagPerChannel", "BB-Flags Versus Channel;Channel;BB Flags Count",kNChannelBins, kChannelMin, kChannelMax );
	h1d->SetMinimum(0);
 	Add2RawsList(h1d,kBBFlagsPerChannel, !expert, image, !saveCorr); iHisto++;

 	h2d = new TH2D("H2D_BBFlagVsClock", "BB-Flags Versus LHC-Clock;Channel;LHC Clocks",kNChannelBins, kChannelMin, kChannelMax,21, -10.5, 10.5 );
 	Add2RawsList(h2d,kBBFlagVsClock, expert, !image, !saveCorr); iHisto++;
	
 	h2d = new TH2D("H2D_BGFlagVsClock", "BG-Flags Versus LHC-Clock;Channel;LHC Clocks",kNChannelBins, kChannelMin, kChannelMax,21, -10.5, 10.5 );
 	Add2RawsList(h2d,kBGFlagVsClock, expert, !image, !saveCorr); iHisto++;
	 
	 
 	AliDebug(AliQAv1::GetQADebugLevel(), Form("%d Histograms has been added to the Raws List",iHisto));
 }

//____________________________________________________________________________ 
void AliVZEROQADataMakerRec::InitDigits()
{
  // create Digits histograms in Digits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1I *fhDigTDC[64]; 
  TH1I *fhDigADC[64]; 
  
  // create Digits histograms in Digits subdir
  TH1I * h0 = new TH1I("hDigitMultiplicity", "Digits multiplicity distribution in VZERO;# of Digits;Entries", 100, 0, 99) ; 
  h0->Sumw2() ;
  Add2DigitsList(h0, 0, !expert, image) ;
  
  for (Int_t i=0; i<64; i++)
    {
    fhDigTDC[i] = new TH1I(Form("hDigitTDC%d", i),Form("Digit TDC in cell %d; TDC value;Entries",i),300,0.,149.);
    
    fhDigADC[i]= new TH1I(Form("hDigitADC%d",i),Form("Digit ADC in cell %d;ADC value;Entries",i),1024,0.,1023.);
    
    Add2DigitsList(fhDigTDC[i],i+1, !expert, image);
    Add2DigitsList(fhDigADC[i],i+1+64, !expert, image);  
    }  
}

//____________________________________________________________________________
void AliVZEROQADataMakerRec::MakeDigits()
{
  // makes data from Digits

  GetDigitsData(0)->Fill(fDigitsArray->GetEntriesFast()) ; 
  TIter next(fDigitsArray) ; 
  AliVZEROdigit *aVZERODigit ; 
  while ( (aVZERODigit = dynamic_cast<AliVZEROdigit *>(next())) ) {
    Int_t   aPMNumber  = aVZERODigit->PMNumber();         
    GetDigitsData(aPMNumber +1)->Fill( aVZERODigit->Time()) ;    // in 100 of picoseconds
    GetDigitsData(aPMNumber +1+64)->Fill( aVZERODigit->ADC()) ;
  }  
}


//____________________________________________________________________________
void AliVZEROQADataMakerRec::MakeDigits(TTree *digitTree)
{
  // makes data from Digit Tree
	
  if ( fDigitsArray ) 
    fDigitsArray->Clear() ; 
  else 
    fDigitsArray = new TClonesArray("AliVZEROdigit", 1000) ; 
  
  TBranch * branch = digitTree->GetBranch("VZERODigit") ;
  if ( ! branch ) {
    AliWarning("VZERO branch in Digit Tree not found") ; 
  } else {
    branch->SetAddress(&fDigitsArray) ;
    branch->GetEntry(0) ; 
    MakeDigits() ; 
  }  
}


//____________________________________________________________________________
void AliVZEROQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  // Creates QA data from ESDs
  
  UInt_t eventType = esd->GetEventType();

  switch (eventType){
  	case PHYSICS_EVENT:
  	AliESDVZERO *esdVZERO=esd->GetVZEROData();
   
	if (!esdVZERO) break;
		  
    	GetESDsData(kCellMultiV0A)->Fill(esdVZERO->GetNbPMV0A());
    	GetESDsData(kCellMultiV0C)->Fill(esdVZERO->GetNbPMV0C());  
    	GetESDsData(kMIPMultiV0A)->Fill(esdVZERO->GetMTotV0A());
    	GetESDsData(kMIPMultiV0C)->Fill(esdVZERO->GetMTotV0C());  
    	
	for(Int_t i=0;i<64;i++) {
			GetESDsData(kMIPMultiChannel)->Fill((Float_t) i,(Float_t) esdVZERO->GetMultiplicity(i));
		   	GetESDsData(kChargeChannel)->Fill((Float_t) i,(Float_t) esdVZERO->GetAdc(i));
			if (i < 32) {
			  if(esdVZERO->BBTriggerV0C(i)) GetESDsData(kBBFlag)->Fill((Float_t) i);
			  if(esdVZERO->BGTriggerV0C(i)) GetESDsData(kBGFlag)->Fill((Float_t) i);
			}
			else {
			  if(esdVZERO->BBTriggerV0A(i-32)) GetESDsData(kBBFlag)->Fill((Float_t) i);  
			  if(esdVZERO->BGTriggerV0A(i-32)) GetESDsData(kBGFlag)->Fill((Float_t) i);
			}		  	
			Float_t time = (Float_t) esdVZERO->GetTime(i); //Convert in ns:  1 TDC channel = 100ps 
			GetESDsData(kTimeChannel)->Fill((Float_t) i,time);
    	}
				
	Float_t timeV0A = esdVZERO->GetV0ATime();
	Float_t timeV0C = esdVZERO->GetV0CTime();
	Float_t diffTime;

	if(timeV0A<-1024.+1.e-6 || timeV0C<-1024.+1.e-6) diffTime = -1024.;
	else diffTime = timeV0A - timeV0C;

	GetESDsData(kESDV0ATime)->Fill(timeV0A);
	GetESDsData(kESDV0CTime)->Fill(timeV0C);
	GetESDsData(kESDDiffTime)->Fill(diffTime);
		
	break;
	}  
  
}

//____________________________________________________________________________
 void AliVZEROQADataMakerRec::MakeRaws(AliRawReader* rawReader)
 {
  // Fills histograms with Raws, computes average ADC values dynamically (pedestal subtracted)
                  
					  
   // Check id histograms already created for this Event Specie
   if ( ! GetRawsData(kPedestalInt0) )
     InitRaws() ;

   rawReader->Reset() ; 
  AliVZERORawStream* rawStream  = new AliVZERORawStream(rawReader); 
 if(!(rawStream->Next())) return;  
 
  eventTypeType eventType = rawReader->GetType();

  Int_t    mulV0A = 0 ; 
  Int_t    mulV0C = 0 ; 
  Double_t timeV0A =0., timeV0C = 0.;
  Double_t weightV0A =0., weightV0C = 0.;
  UInt_t   itimeV0A=0, itimeV0C=0;
  Double_t chargeV0A=0., chargeV0C=0.;
  Double_t mipV0A=0., mipV0C=0.;

  Double_t diffTime=-100000.;

  
  switch (eventType){
       case PHYSICS_EVENT:
  
  		fNTotEvents++;

       Int_t  iFlag=0;
       Int_t  pedestal;
       Int_t  integrator;
       Bool_t flagBB[64];	 
       Bool_t flagBG[64];	 
       Int_t  mbCharge;
	   Float_t charge;
       Int_t  offlineCh;
       Float_t adc[64], time[64], width[64], timeCorr[64]; 

       for(Int_t iChannel=0; iChannel<64; iChannel++) { // BEGIN : Loop over channels
		   
	   offlineCh = rawStream->GetOfflineChannel(iChannel);
		   
	   // Fill Pedestal histograms
	   
           for(Int_t j=15; j<21; j++) {
		       if((rawStream->GetBGFlag(iChannel,j) || rawStream->GetBBFlag(iChannel,j))) iFlag++;
           }

           if(iFlag == 0){ //No Flag found
		       for(Int_t j=15; j<21; j++){
	   		       pedestal= (Int_t) rawStream->GetPedestal(iChannel, j);
	   		       integrator = rawStream->GetIntegratorFlag(iChannel, j);

	   		       GetRawsData((integrator == 0 ? kPedestalInt0 : kPedestalInt1))->Fill(offlineCh,pedestal);
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
	       Int_t board = AliVZEROCalibData::GetBoardNumber(offlineCh);
	       time[offlineCh] = rawStream->GetTime(iChannel)*fCalibData->GetTimeResolution(board);
	       width[offlineCh] = rawStream->GetWidth(iChannel)*fCalibData->GetWidthResolution(board);

	       if (time[offlineCh] >= 1e-6) GetRawsData(kChargeEoI)->Fill(offlineCh,adc[offlineCh]);

           GetRawsData((integrator == 0 ? kChargeEoIInt0 : kChargeEoIInt1))->Fill(offlineCh,charge);
	   	   if(flagBB[offlineCh]) GetRawsData((integrator == 0 ? kChargeEoIBBInt0 : kChargeEoIBBInt1))->Fill(offlineCh,charge);
           if(flagBG[offlineCh]) GetRawsData((integrator == 0 ? kChargeEoIBGInt0 : kChargeEoIBGInt1))->Fill(offlineCh,charge);

			Float_t sigma = fCalibData->GetSigma(offlineCh+64*integrator);

		   
	   // Calculation of the number of MIP
	   Double_t mipEoI = adc[offlineCh] * fCalibData->GetMIPperADC(offlineCh);

	   if((adc[offlineCh] > 2.*sigma) && !(time[offlineCh] <1.e-6)){ 
		   ((TH2D*)GetRawsData(kRawMIPChannel))->Fill(offlineCh,mipEoI);
        	   if(offlineCh<32) {
				   mulV0C++;
				   chargeV0C += adc[offlineCh];
				   mipV0C += mipEoI;
		   } else {
				   mulV0A++;
				   chargeV0A += adc[offlineCh];
				   mipV0A += mipEoI;
		   }
	   }

	   // Fill Charge Minimum Bias Histograms
		   
	   int idx;
	   for(Int_t iBunch=0; iBunch<10; iBunch++){
			   integrator = rawStream->GetIntMBFlag(iChannel, iBunch);
			   bool bbFlag     = rawStream->GetBBMBFlag(iChannel, iBunch);
			   bool bgFlag     = rawStream->GetBGMBFlag(iChannel, iBunch);
			   mbCharge   = rawStream->GetChargeMB(iChannel, iBunch);

			   if(integrator==0){
				   if(bbFlag==0){
					   if(bgFlag==0) idx = kChargeMBBB0BG0Int0;
					   else idx = kChargeMBBB0BG1Int0;
				   } else {
					   if(bgFlag==0) idx = kChargeMBBB1BG0Int0;
					   else idx = kChargeMBBB1BG1Int0;
				   }
			   } else {
				   if(bbFlag==0){
					   if(bgFlag==0) idx = kChargeMBBB0BG0Int1;
					   else idx = kChargeMBBB0BG1Int1;
				   } else {
					   if(bgFlag==0) idx = kChargeMBBB1BG0Int1;
					   else idx = kChargeMBBB1BG1Int1;
				   }
			   }
			   GetRawsData(idx)->Fill(offlineCh,mbCharge);
       }   

	  // Fill HPTDC Time Histograms
	   timeCorr[offlineCh] = CorrectLeadingTime(offlineCh,time[offlineCh],adc[offlineCh]);

  	   const Float_t p1 = 2.50; // photostatistics term in the time resolution
 	   const Float_t p2 = 3.00; // slewing related term in the time resolution
	   if(timeCorr[offlineCh]>-1024 + 1.e-6){
			Float_t nphe = adc[offlineCh]*kChargePerADC/(fCalibData->GetGain(offlineCh)*TMath::Qe());
			Float_t timeErr = 0;
			if (nphe>1.e-6) timeErr = TMath::Sqrt(kIntTimeRes*kIntTimeRes+
				      p1*p1/nphe+
				      p2*p2*(fTimeSlewing->GetParameter(0)*fTimeSlewing->GetParameter(1))*(fTimeSlewing->GetParameter(0)*fTimeSlewing->GetParameter(1))*
				      TMath::Power(adc[offlineCh]/fCalibData->GetCalibDiscriThr(offlineCh,kTRUE),2.*(fTimeSlewing->GetParameter(1)-1.))/
				      (fCalibData->GetCalibDiscriThr(offlineCh,kTRUE)*fCalibData->GetCalibDiscriThr(offlineCh,kTRUE)));

			if (timeErr>1.e-6) {
			  if (offlineCh<32) {
			    itimeV0C++;
			    timeV0C += timeCorr[offlineCh]/(timeErr*timeErr);
			    weightV0C += 1./(timeErr*timeErr);
			  }else{
			    itimeV0A++;
			    timeV0A += timeCorr[offlineCh]/(timeErr*timeErr);
			    weightV0A += 1./(timeErr*timeErr);
			  }
			}
	   }
		GetRawsData(kHPTDCTime)->Fill(offlineCh,timeCorr[offlineCh]);
		GetRawsData(kWidth)->Fill(offlineCh,width[offlineCh]);
        if(flagBB[offlineCh]) {
			GetRawsData(kHPTDCTimeBB)->Fill(offlineCh,timeCorr[offlineCh]);
			GetRawsData(kWidthBB)->Fill(offlineCh,width[offlineCh]);
		}
		if(flagBG[offlineCh]) {
			GetRawsData(kHPTDCTimeBG)->Fill(offlineCh,timeCorr[offlineCh]);
			GetRawsData(kWidthBG)->Fill(offlineCh,width[offlineCh]);
		}

	   // Fill Flag and Charge Versus LHC-Clock histograms
	   
	   for(Int_t iEvent=0; iEvent<21; iEvent++){
               charge = rawStream->GetPedestal(iChannel,iEvent);
               integrator = rawStream->GetIntegratorFlag(iChannel,iEvent);
               bool bbFlag	  = rawStream->GetBBFlag(iChannel, iEvent);
               bool bgFlag	  = rawStream->GetBGFlag(iChannel,iEvent );

               ((TH2*) GetRawsData((integrator == 0 ? kChargeVsClockInt0 : kChargeVsClockInt1 )))->Fill(offlineCh,(float)iEvent-10,(float)charge);
               ((TH2*) GetRawsData(kBBFlagVsClock))->Fill(offlineCh,(float)iEvent-10,(float)bbFlag);
               ((TH2*) GetRawsData(kBGFlagVsClock))->Fill(offlineCh,(float)iEvent-10,(float)bgFlag);
               if(iEvent==10) ((TH1*) GetRawsData(kBBFlagsPerChannel))->Fill(offlineCh,(float)bbFlag);
           }

       }// END of Loop over channels

		if(weightV0A>1.e-6) timeV0A /= weightV0A; 
		else timeV0A = -1024.;
		if(weightV0C>1.e-6) timeV0C /= weightV0C;
		else timeV0C = -1024.;
		if(timeV0A<-1024.+1.e-6 || timeV0C<-1024.+1.e-6) diffTime = -1024.;
		else diffTime = timeV0A - timeV0C;

		Bool_t v0ABB = kFALSE;
		Bool_t v0CBB = kFALSE;
		Bool_t v0ABG = kFALSE;
		Bool_t v0CBG = kFALSE;
		
		if(timeV0A>kMinBBA && timeV0A<kMaxBBA) {
			v0ABB = kTRUE;
		} else if(timeV0A>kMinBGA && timeV0A<kMaxBGA) {
			v0ABG = kTRUE;
		}
		if(timeV0C>kMinBBC && timeV0C<kMaxBBC) {
			v0CBB = kTRUE;
		} else if(timeV0C>kMinBGC && timeV0C<kMaxBGC) {
			v0CBG = kTRUE;
		}

// Fill Trigger output histogram
		if(v0ABB && v0CBB) GetRawsData(kTriggers)->Fill(0);
		if((v0ABB || v0CBB) && !(v0ABG || v0CBG)) GetRawsData(kTriggers)->Fill(1);
		if(v0ABG && v0CBB) GetRawsData(kTriggers)->Fill(2);
		if(v0ABB && v0CBG) GetRawsData(kTriggers)->Fill(3);
		

		GetRawsData(kV0ATime)->Fill(timeV0A);
		GetRawsData(kV0CTime)->Fill(timeV0C);
		GetRawsData(kDiffTime)->Fill(diffTime);
		GetRawsData(kTimeV0AV0C)->Fill(timeV0A,timeV0C);

		GetRawsData(kMultiV0A)->Fill(mulV0A);
		GetRawsData(kMultiV0C)->Fill(mulV0C);

		GetRawsData(kChargeV0A)->Fill(chargeV0A);
		GetRawsData(kChargeV0C)->Fill(chargeV0C);
		GetRawsData(kChargeV0)->Fill(chargeV0A + chargeV0C);

		GetRawsData(kRawMIPV0A)->Fill(mipV0A);
		GetRawsData(kRawMIPV0C)->Fill(mipV0C);
		GetRawsData(kRawMIPV0)->Fill(mipV0A + mipV0C);
		break;
	    
	} // END of SWITCH : EVENT TYPE 
	
	fEvent++; 
	TParameter<double> * p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), GetRawsData(kMultiV0A)->GetName()))) ; 
	if (p) p->SetVal((double)mulV0A) ; 

	p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), GetRawsData(kMultiV0C)->GetName()))) ; 
	if (p) p->SetVal((double)mulV0C) ;                     

	p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), GetRawsData(kChargeV0A)->GetName()))) ; 
	if (p) p->SetVal((double)chargeV0A) ; 

	p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), GetRawsData(kChargeV0C)->GetName()))) ; 
	if (p) p->SetVal((double)chargeV0C) ;                     

	p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), GetRawsData(kChargeV0)->GetName()))) ; 
	if (p) p->SetVal((double)(chargeV0A + chargeV0C)) ;                     
	
	p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), GetRawsData(kRawMIPV0A)->GetName()))) ; 
	if (p) p->SetVal((double)mipV0A) ; 
	
	p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), GetRawsData(kRawMIPV0C)->GetName()))) ; 
	if (p) p->SetVal((double)mipV0C) ;                     
	
	p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), GetRawsData(kRawMIPV0)->GetName()))) ; 
	if (p) p->SetVal((double)(mipV0A + mipV0C)) ;                     
	
	p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), GetRawsData(kV0ATime)->GetName()))) ; 
	if (p) p->SetVal((double)timeV0A) ; 
	
	p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), GetRawsData(kV0CTime)->GetName()))) ; 
	if (p) p->SetVal((double)timeV0C) ;                     
	
	p = dynamic_cast<TParameter<double>*>(GetParameterList()->FindObject(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(AliQAv1::kRAWS).Data(), GetRawsData(kDiffTime)->GetName()))) ; 
	if (p) p->SetVal((double)diffTime) ;                     
	
  	delete rawStream; rawStream = 0x0;      


 }

//____________________________________________________________________________ 
void AliVZEROQADataMakerRec::StartOfDetectorCycle()
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

  AliCDBEntry *entry2 = AliCDBManager::Instance()->Get("VZERO/Calib/TimeDelays");
  if (!entry2) AliFatal("VZERO time delays are not found in OCDB !");
  TH1F *delays = (TH1F*)entry2->GetObject();

  AliCDBEntry *entry3 = AliCDBManager::Instance()->Get("VZERO/Calib/TimeSlewing");
  if (!entry3) AliFatal("VZERO time slewing function is not found in OCDB !");
  fTimeSlewing = (TF1*)entry3->GetObject();

  for(Int_t i = 0 ; i < 64; ++i) {
    //Int_t board = AliVZEROCalibData::GetBoardNumber(i);
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
 
  fNTotEvents = 0;
}


//-------------------------------------------------------------------------------------------------
Float_t AliVZEROQADataMakerRec::CorrectLeadingTime(Int_t i, Float_t time, Float_t adc) const
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

