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

///////////////////////////////////////////////////////////////////////
//                                                                   //
//  Produces the data needed to calculate the ZDC quality assurance. //
//  QA objects are 1 & 2 Dimensional histograms.                     //
//  author: C. Oppedisano                                            //
//                                                                   //
///////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH2F.h>
#include <TLine.h>
#include <TProfile.h>
#include <Riostream.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQAChecker.h"
#include "AliZDCReco.h"
#include "AliRawReader.h"
#include "AliZDCQADataMakerRec.h"
#include "AliZDCPedestals.h"
#include "AliZDCRawStream.h"
#include "AliZDCDigit.h"
#include "AliESDZDC.h"
#include "AliESDEvent.h"

ClassImp(AliZDCQADataMakerRec)
           
//____________________________________________________________________________ 
  AliZDCQADataMakerRec::AliZDCQADataMakerRec() : 
  AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kZDC), "ZDC Quality Assurance Data Maker"), 
  fPedCalibData(0x0)
{
  // ctor
}

//____________________________________________________________________________ 
AliZDCQADataMakerRec::AliZDCQADataMakerRec(const AliZDCQADataMakerRec& qadm) :
  AliQADataMakerRec(),      
  fPedCalibData(qadm.fPedCalibData)

{
  //copy ctor 
  SetName((const char*)qadm.GetName()); 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliZDCQADataMakerRec& AliZDCQADataMakerRec::operator = (const AliZDCQADataMakerRec& qadm )
{
  // Equal operator.
  this->~AliZDCQADataMakerRec();
  new(this) AliZDCQADataMakerRec(qadm);
  return *this;
}

//____________________________________________________________________________ 
AliZDCQADataMakerRec::~AliZDCQADataMakerRec()
{
  if(fPedCalibData && !(AliCDBManager::Instance()->GetCacheFlag())){
    delete fPedCalibData;
    fPedCalibData=0;
  } 
}

//____________________________________________________________________________ 
AliZDCPedestals* AliZDCQADataMakerRec::GetPedCalibData() const
{

  // Retrieving pedestal calibration object from OCDB
  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("ZDC/Calib/Pedestals");
  if(!entry) AliWarning("No calibration data loaded!");  

  AliZDCPedestals *calibdata = (AliZDCPedestals*)  (entry->GetObject());
  if(!calibdata) AliFatal("Wrong calibration object in calibration  file!");

  return calibdata;

}

//____________________________________________________________________________ 
void AliZDCQADataMakerRec::InitDigits()
{
  // create Digits histograms in Digits subdir
  //
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  // ------------------- HIGH GAIN CHAIN ---------------------------
  TH1F * hDigZNCTot = new TH1F("hDigZNCTot", "Signal in ZNC;Amplitude [ADC counts];Counts", 100, 0., 6000.);
  TH1F * hDigZNATot = new TH1F("hDigZNATot", "Signal in ZNA;Amplitude [ADC counts];Counts", 100, 0., 6000.);
  TH1F * hDigZPCTot = new TH1F("hDigZPCTot", "Signal in ZPC;Amplitude [ADC counts];Counts", 100, 0., 6000.);
  TH1F * hDigZPATot = new TH1F("hDigZPATot", "Signal in ZPA;Amplitude [ADC counts];Counts", 100, 0., 6000.);
  Add2DigitsList(hDigZNCTot, 0, !expert, image);
  Add2DigitsList(hDigZNATot, 1, !expert, image);
  Add2DigitsList(hDigZPCTot, 2, !expert, image);
  Add2DigitsList(hDigZPATot, 3, !expert, image);
  //
  TH1F * hDigSumQZNC = new TH1F("hDigSumQZNC", "Signal in 4 ZNC PMQ;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hDigSumQZNA = new TH1F("hDigSumQZNA", "Signal in 4 ZNA PMQ;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hDigSumQZPC = new TH1F("hDigSumQZPC", "Signal in 4 ZPC PMQ;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hDigSumQZPA = new TH1F("hDigSumQZPA", "Signal in 4 ZPA PMQ;Amplitude [ADC counts];Counts",100, 0., 4000.);
  Add2DigitsList(hDigSumQZNC, 4, expert, !image);
  Add2DigitsList(hDigSumQZNA, 5, expert, !image);
  Add2DigitsList(hDigSumQZPC, 6, expert, !image);
  Add2DigitsList(hDigSumQZPA, 7, expert, !image);
  //
  TH1F * hDigPMCZNC = new TH1F("hDigPMCZNC", "Signal in ZNC PMC;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hDigPMCZNA = new TH1F("hDigPMCZNA", "Signal in ZNA PMC;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hDigPMCZPC = new TH1F("hDigPMCZPC", "Signal in ZPC PMC;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hDigPMCZPA = new TH1F("hDigPMCZPA", "Signal in ZPA PMC;Amplitude [ADC counts];Counts",100, 0., 4000.);
  Add2DigitsList(hDigPMCZNC, 8, expert, !image);
  Add2DigitsList(hDigPMCZNA, 9, expert, !image);
  Add2DigitsList(hDigPMCZPC, 10, expert, !image);
  Add2DigitsList(hDigPMCZPA, 11, expert, !image);
  // 
  ClonePerTrigClass(AliQAv1::kDIGITS); // this should be the last line
}

//____________________________________________________________________________
void AliZDCQADataMakerRec::InitRaws()
{
  // create Digits histograms in Digits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  TH1F * hZNCSpectrum = new TH1F("hZNCSpectrum","ZNC spectrum;Amplitude [ADC counts];Counts",100,0.,1200.);
  TH1F * hZPCSpectrum = new TH1F("hZPCSpectrum","ZPC spectrum;Amplitude [ADC counts];Counts",100,0.,1200.);
  TH1F * hZNASpectrum = new TH1F("hZNASpectrum","ZNA spectrum;Amplitude [ADC counts];Counts",100,0.,1200.);
  TH1F * hZPASpectrum = new TH1F("hZPASpectrum","ZPA spectrum;Amplitude [ADC counts];Counts",100,0.,1200.);
  // Booking from ch. 8 for checked signals to avoid running QA on pedestals!
  TH1F * hZEM1Spectrum = new TH1F("hZEM1Spectrum","ZEM1 spectrum;Amplitude [ADC counts];Counts",100,8., 1208.);
  TH1F * hZEM2Spectrum = new TH1F("hZEM2Spectrum","ZEM2 spectrum;Amplitude [ADC counts];Counts",100,8., 1208.);
  Add2RawsList(hZNCSpectrum, 0, expert, !image);
  Add2RawsList(hZNASpectrum, 1, expert, !image);
  Add2RawsList(hZPCSpectrum, 2, expert, !image);
  Add2RawsList(hZPASpectrum, 3, expert, !image);
  Add2RawsList(hZEM1Spectrum, 4, expert, !image);
  Add2RawsList(hZEM2Spectrum, 5, expert, !image);
    
  TH1F * hRawPMCZNC = new TH1F("hRawPMCZNC", "Raw ZNC PMC;Amplitude [ADC counts];Counts",100, 8., 1208.);
  TH1F * hRawPMCZNA = new TH1F("hRawPMCZNA", "Raw ZNA PMC;Amplitude [ADC counts];Counts",100, 8., 1208.);
  TH1F * hRawPMCZPC = new TH1F("hRawPMCZPC", "Raw ZPC PMC;Amplitude [ADC counts];Counts",100, 8., 1208.);
  TH1F * hRawPMCZPA = new TH1F("hRawPMCZPA", "Raw ZPA PMC;Amplitude [ADC counts];Counts",100, 8., 1208.);
  Add2RawsList(hRawPMCZNC, 6, expert, !image);
  Add2RawsList(hRawPMCZNA, 7, expert, !image);
  Add2RawsList(hRawPMCZPC, 8, expert, !image);
  Add2RawsList(hRawPMCZPA, 9, expert, !image);
  TH1F * hRawSumQZNC = new TH1F("hRawSumQZNC", "Raw sumQ ZNC;Amplitude [ADC counts];Counts",100, 0., 1200.);
  TH1F * hRawSumQZNA = new TH1F("hRawSumQZNA", "Raw sumQ ZNA;Amplitude [ADC counts];Counts",100, 0., 1200.);
  TH1F * hRawSumQZPC = new TH1F("hRawSumQZPC", "Raw sumQ ZPC;Amplitude [ADC counts];Counts",100, 0., 1200.);
  TH1F * hRawSumQZPA = new TH1F("hRawSumQZPA", "Raw sumQ ZPA;Amplitude [ADC counts];Counts",100, 0., 1200.);
  Add2RawsList(hRawSumQZNC, 10, expert, !image);
  Add2RawsList(hRawSumQZNA, 11, expert, !image);
  Add2RawsList(hRawSumQZPC, 12, expert, !image);
  Add2RawsList(hRawSumQZPA, 13, expert, !image);
  
  TH1F * hRawTDCZEM1 = new TH1F("hRawTDCZEM1", "TDC ZEM1;TDC [ns]",160, -350., -310.);
  Add2RawsList(hRawTDCZEM1, 14, expert, !image);
  TH1F * hRawTDCZPC = new TH1F("hRawTDCZPC", "TDC ZPC;TDC [ns]",160, -350., -310.);
  Add2RawsList(hRawTDCZPC, 15, expert, !image);
  
  TProfile * hRawADCProfs = new TProfile("hRawADCProfs", "ADC profiles;ADC id;Mean ADC values",22,-0.5,21.5,10.,1210.,"");
  Add2RawsList(hRawADCProfs, 16, expert, !image);
  TProfile * hRawTDCProfs = new TProfile("hRawTDCProfs", "TDC profiles;TDC id;Mean TDC values",6,0.5,6.5,-340.,-300.,"S");
  Add2RawsList(hRawTDCProfs, 17, expert, !image);
  
  TH1F * hRawADCs = new TH1F("hRawADCs", "ADCs;ADC id;Mean ADC values",22,-0.5,21.5);
  Add2RawsList(hRawADCs, 18, !expert, image);
 
  TH1F * hRawTDCs = new TH1F("hRawTDCs", "TDCs;TDC id;Mean TDC values",6,0.5,6.5);
//  hRawTDCs->SetMaximum(-300); hRawTDCs->SetMinimum(-340);
  Add2RawsList(hRawTDCs, 19, !expert, image);
  
  TH2F *hZNCrawCentr  = new TH2F("hZNCrawCentr", "ZNC centroid;X (cm);Y(cm)", 100,-3.5,3.5,100,-3.5,3.5);
  Add2RawsList(hZNCrawCentr, 20, expert, image);
  TH2F *hZNArawCentr  = new TH2F("hZNArawCentr", "ZNA centroid;X (cm);Y(cm)", 100,-3.5,3.5,100,-3.5,3.5);
  Add2RawsList(hZNArawCentr, 21, expert, image);
  
  TH2F *hTimeZDC = new TH2F("hTimeZDC", "ZDC timing;Z_{vertex}/c (ns);time (ns)", 60,-30.,30.,120,-60,-60);
  Add2RawsList(hTimeZDC, 22, !expert, image);
  //
  ClonePerTrigClass(AliQAv1::kRAWS); // this should be the last line
}

//____________________________________________________________________________
void AliZDCQADataMakerRec::InitRecPoints()
{
  // create Digits histograms in Digits subdir
  const Bool_t expert = kTRUE ; 
  const Bool_t image  = kTRUE ; 
  //
  // ------------------- HIGH GAIN CHAIN ---------------------------
  TH1F * hRecZNCTot = new TH1F("hRecZNCTot", "Rec signal in ZNC;Amplitude [ADC counts];Counts", 100, 0., 2000.);
  TH1F * hRecZNATot = new TH1F("hRecZNATot", "Rec signal in ZNA;Amplitude [ADC counts];Counts", 100, 0., 2000.);
  TH1F * hRecZPCTot = new TH1F("hRecZPCTot", "Rec signal in ZPC;Amplitude [ADC counts];Counts", 100, 0., 8000.);
  TH1F * hRecZPATot = new TH1F("hRecZPATot", "Rec signal in ZPA;Amplitude [ADC counts];Counts", 100, 0., 8000.);
  Add2RecPointsList(hRecZNCTot, 0, !expert, image);
  Add2RecPointsList(hRecZNATot, 1, !expert, image);
  Add2RecPointsList(hRecZPCTot, 2, !expert, image);
  Add2RecPointsList(hRecZPATot, 3, !expert, image);
  //
  TH1F * hRecSumQZNC = new TH1F("hRecSumQZNC", "Rec summed 4 ZNC quadrants;Amplitude [ADC counts];Counts",100, 0., 2000.);
  TH1F * hRecSumQZNA = new TH1F("hRecSumQZNA", "Rec summed 4 ZNA quadrants;Amplitude [ADC counts];Counts",100, 0., 2000.);
  TH1F * hRecSumQZPC = new TH1F("hRecSumQZPC", "Rec summed 4 ZPC quadrants;Amplitude [ADC counts];Counts",100, 0., 2000.);
  TH1F * hRecSumQZPA = new TH1F("hRecSumQZPA", "Rec summed 4 ZPA quadrants;Amplitude [ADC counts];Counts",100, 0., 2000.);
  Add2RecPointsList(hRecSumQZNC, 4, expert, !image);
  Add2RecPointsList(hRecSumQZNA, 5, expert, !image);
  Add2RecPointsList(hRecSumQZPC, 6, expert, !image);
  Add2RecPointsList(hRecSumQZPA, 7, expert, !image);
  //
  TH1F * hRecPMCZNC = new TH1F("hRecPMCZNC", "Rec common ZNC PMT;Amplitude [ADC counts];Counts",100, 0., 2000.);
  TH1F * hRecPMCZNA = new TH1F("hRecPMCZNA", "Rec common ZNA PMT;Amplitude [ADC counts];Counts",100, 0., 2000.);
  TH1F * hRecPMCZPC = new TH1F("hRecPMCZPC", "Rec common ZPC PMT;Amplitude [ADC counts];Counts",100, 0., 2000.);
  TH1F * hRecPMCZPA = new TH1F("hRecPMCZPA", "Rec common ZPA PMT;Amplitude [ADC counts];Counts",100, 0., 2000.);
  Add2RecPointsList(hRecPMCZNC, 8 , expert, !image);
  Add2RecPointsList(hRecPMCZNA, 9 , expert, !image);
  Add2RecPointsList(hRecPMCZPC, 10, expert, !image);
  Add2RecPointsList(hRecPMCZPA, 11, expert, !image); 
  //
  ClonePerTrigClass(AliQAv1::kRECPOINTS); // this should be the last line
}

//____________________________________________________________________________
void AliZDCQADataMakerRec::InitESDs()
{
  //Booking ESDs histograms
  //
  const Bool_t expert = kTRUE ; 
  const Bool_t image  = kTRUE ; 
  
  // ------------------- HIGH GAIN CHAIN ---------------------------
  TH1F * hESDZNCTot = new TH1F("hESDZNCTot", "Energy in ZNC", 100, 0., 4000.);
  TH1F * hESDZNATot = new TH1F("hESDZNATot", "Energy in ZNA", 100, 0., 4000.);
  TH1F * hESDZPCTot = new TH1F("hESDZPCTot", "Energy in ZPC", 100, 0., 4000.);
  TH1F * hESDZPATot = new TH1F("hESDZPATot", "Energy in ZPA", 100, 0., 4000.);
  Add2ESDsList(hESDZNCTot, 0, !expert, image);
  Add2ESDsList(hESDZNATot, 1, !expert, image);
  Add2ESDsList(hESDZPCTot, 2, !expert, image);
  Add2ESDsList(hESDZPATot, 3, !expert, image);
  //
  TH1F * hESDZEM1 = new TH1F("hESDZEM1", "Energy in ZEM1", 100, 0., 2000.);
  TH1F * hESDZEM2 = new TH1F("hESDZEM2", "Energy in ZEM2", 100, 0., 2000.);
  Add2ESDsList(hESDZEM1,4, !expert, image);
  Add2ESDsList(hESDZEM2,5, !expert, image);
  //
  TH1F * hESDSumQZNC = new TH1F("hESDSumQZNC", "Sum of 4 ZNC energy",100, 0., 2000.);
  TH1F * hESDSumQZNA = new TH1F("hESDSumQZNA", "Sum of 4 ZNA energy",100, 0., 2000.);
  TH1F * hESDSumQZPC = new TH1F("hESDSumQZPC", "Sum of 4 ZPC energy",100, 0., 2000.);
  TH1F * hESDSumQZPA = new TH1F("hESDSumQZPA", "Sum of 4 ZPA energy",100, 0., 2000.);
  Add2ESDsList(hESDSumQZNC, 6, expert, !image);
  Add2ESDsList(hESDSumQZNA, 7, expert, !image);
  Add2ESDsList(hESDSumQZPC, 8, expert, !image);
  Add2ESDsList(hESDSumQZPA, 9, expert, !image);
  //
  TH1F * hESDPMCZNC = new TH1F("hESDPMCZNC", "Energy in ZNC PMC",100, 0., 2000.);
  TH1F * hESDPMCZNA = new TH1F("hESDPMCZNA", "Energy in ZNA PMC",100, 0., 2000.);
  TH1F * hESDPMCZPC = new TH1F("hESDPMCZPC", "Energy in ZPC PMC",100, 0., 2000.);
  TH1F * hESDPMCZPA = new TH1F("hESDPMCZPA", "Energy in ZPA PMC",100, 0., 2000.);
  Add2ESDsList(hESDPMCZNC, 10, expert, !image);
  Add2ESDsList(hESDPMCZNA, 11, expert, !image);
  Add2ESDsList(hESDPMCZPC, 12, expert, !image);
  Add2ESDsList(hESDPMCZPA, 13, expert, !image);
  // 
  ClonePerTrigClass(AliQAv1::kESDS); // this should be the last line
}

//___________________________________________________________________________
void AliZDCQADataMakerRec::MakeDigits(TTree *digitTree)
{
  // makes data from Digit Tree
  if(!GetDigitsData(0)) InitDigits();

  if(!digitTree){
    AliError("Can't get ZDC digit tree!!");
    return; 
  }	
   
  TBranch * branch = digitTree->GetBranch("ZDC");
  if(!branch){
    AliError("ZDC branch in digit tree not found"); 
    return;
  } 
    
  AliZDCDigit *digit = 0x0;
  branch->SetAddress(&digit);
     
  Float_t adcSum_ZNC=0., adcSum_ZNA=0., adcSum_ZPC=0., adcSum_ZPA=0.;
  Float_t adcSumQ_ZNC=0., adcSumQ_ZNA=0., adcSumQ_ZPC=0., adcSumQ_ZPA=0.;
  //Float_t adcSum_ZNC_lg=0., adcSum_ZNA_lg=0., adcSum_ZPC_lg=0., adcSum_ZPA_lg=0.;
  //Float_t adcSumQ_ZNC_lg=0., adcSumQ_ZNA_lg=0., adcSumQ_ZPC_lg=0., adcSumQ_ZPA_lg=0.;
  
  Int_t ndig = digitTree->GetEntries();
  for(Int_t i=0; i<ndig; i++){
      branch->GetEntry(i);
      
      if(digit->GetSector(0)==1 && digit->GetSector(1)!=5){
	  adcSum_ZNC += digit->GetADCValue(0);
	  //adcSum_ZNC_lg += digit->GetADCValue(1);
	  //
	  if(digit->GetSector(1)!=0){
	      adcSumQ_ZNC += digit->GetADCValue(0);
	      //adcSumQ_ZNC_lg+= digit->GetADCValue(1);
	  }
	  else{
	      FillDigitsData(8,digit->GetADCValue(0));
	      //FillDigitsData(20,digit->GetADCValue(1));
	  }
      }
      else if(digit->GetSector(0)==2){
	  adcSum_ZPC += digit->GetADCValue(0);
	  //adcSum_ZPC_lg += digit->GetADCValue(1);
	  //
	  if(digit->GetSector(1)!=0){
	      adcSumQ_ZPC += digit->GetADCValue(0);
	      //adcSumQ_ZPC_lg+= digit->GetADCValue(1);
	  }
	  else{
	      FillDigitsData(10,digit->GetADCValue(0));
	      //FillDigitsData(22,digit->GetADCValue(1));
	  }
      }
      else if(digit->GetSector(0)==4 && digit->GetSector(1)!=5){
	  adcSum_ZNA += digit->GetADCValue(0);
	  //adcSum_ZNA_lg += digit->GetADCValue(1);
	  //
	  if(digit->GetSector(1)!=0){
	      adcSumQ_ZNA += digit->GetADCValue(0);
	      //adcSumQ_ZNA_lg+= digit->GetADCValue(1);
	  }
	  else{
	      FillDigitsData(9,digit->GetADCValue(0));
	      //FillDigitsData(21,digit->GetADCValue(1));
	  }
      }
      else if(digit->GetSector(0)==5){
	  adcSum_ZPA += digit->GetADCValue(0);
	  //adcSum_ZPA_lg += digit->GetADCValue(1);
	  //
	  if(digit->GetSector(1)!=0){
	      adcSumQ_ZPA += digit->GetADCValue(0);
	      //adcSumQ_ZPA_lg+= digit->GetADCValue(1);
	  }
	  else{
	      FillDigitsData(11,digit->GetADCValue(0));
	      //FillDigitsData(23,digit->GetADCValue(1));
	  }
      }
  }
  //
  FillDigitsData(0,adcSum_ZNC);
  FillDigitsData(1,adcSum_ZNA);
  FillDigitsData(2,adcSum_ZPC);
  FillDigitsData(3,adcSum_ZPA);
  //
  FillDigitsData(4,adcSumQ_ZNC);
  FillDigitsData(5,adcSumQ_ZNA);
  FillDigitsData(6,adcSumQ_ZPC);
  FillDigitsData(7,adcSumQ_ZPA);
  
  delete digit;
  digit=0;
  //
  IncEvCountCycleDigits();
  IncEvCountTotalDigits();
  //
}


//____________________________________________________________________________
void AliZDCQADataMakerRec::MakeRaws(AliRawReader *rawReader)
{
  // Filling Raws QA histos
  //
  // Checking the event type 
//  if (rawReader->GetType()!=7){
  
    // Check if histograms already created for this Event Specie
    if(!GetRawsData(0)) InitRaws();
    // Parameters for mean value pedestal subtraction
    int const kNch = 24;
    Float_t meanPed[2*kNch];    
    for(Int_t jj=0; jj<2*kNch; jj++) meanPed[jj] = fPedCalibData->GetMeanPed(jj);
    
    Float_t zncSignal=0., znaSignal=0., zpcSignal=0., zpaSignal=0.;
    Float_t zncSumQ=0., znaSumQ=0., zpcSumQ=0., zpaSumQ=0.;
    Float_t zncpmC=0., znapmC=0., zpcpmC=0., zpapmC=0.;
    Bool_t isZNCFired=kFALSE, isZPCFired=kFALSE, isZNAFired=kFALSE, isZPAFired=kFALSE;
    Int_t  indZNC=0, indZNA=0, indZPC=0, indZPA=0;
    Float_t zncTDC[10], zpcTDC[10], zem1TDC[10], zem2TDC[10], znaTDC[10], zpaTDC[10];
    Float_t zncSumTDC[10], znaSumTDC[10];
    for(Int_t i=0; i<10; i++){
       zncTDC[i]=zpcTDC[i]=zem1TDC[i]=zem2TDC[i]=znaTDC[i]=zpaTDC[i]=zncSumTDC[i]=znaSumTDC[i]=-999.;
    }
    Float_t tdcGate=-999., l0=-999.;
    Int_t iMultZNCTDC=0, iMultZPCTDC=0, iMultZEM1TDC=0, iMultZEM2TDC=0, iMultZNATDC=0, iMultZPATDC=0;
    Int_t iMultTDCC=0, iMultTDCA=0;
    
    const Float_t x[4] = {-1.75, 1.75, -1.75, 1.75};
    const Float_t y[4] = {-1.75, -1.75, 1.75, 1.75};
    const Float_t alpha=0.5;
    Float_t numXZNC=0., numYZNC=0., denZNC=0., wZNC=0.; 
    Float_t numXZNA=0., numYZNA=0., denZNA=0., wZNA=0.; 
    
    rawReader->Reset();
    AliZDCRawStream stream(rawReader);
    while(stream.Next()){

      if(stream.IsADCDataWord() && 
    	 (stream.GetADCModule()==0 || stream.GetADCModule()==1)){
       
    	 Int_t det = stream.GetSector(0);
    	 Int_t quad = stream.GetSector(1);
    	 Int_t gain = stream.GetADCGain();
    	 Int_t pedindex=0;
    	 
    	 // Stuff for pedestal subtraction
    	 if(quad != 5){ // ZDCs (not reference PTMs)
	  Float_t rawVal=-99., pedSubVal=-99.;
    	  if(det == 1){    
    	    pedindex = quad;
    	    if(gain == 0){
	      rawVal = (Float_t) (stream.GetADCValue());
    	      pedSubVal = (Float_t) (rawVal-meanPed[pedindex]); 
    	      zncSignal  += pedSubVal; 
    	      isZNCFired = kTRUE;
    	      if(quad!=0){
	        zncSumQ += pedSubVal;
	        if(pedSubVal>0.&& zncpmC>7.){
	          wZNC = TMath::Power(pedSubVal, alpha);
	          numXZNC += x[quad-1]*wZNC;
	          numYZNC += y[quad-1]*wZNC;
	          denZNC += wZNC;
	        }
	      }
    	      else{
  		zncpmC = pedSubVal;
  		FillRawsData(6,zncpmC);
  	      }
	      indZNC++;
	      
	      FillRawsData(16, pedindex, pedSubVal);
  	    }
  	  }
  	  else if(det == 2){ 
  	    pedindex = quad+5;
  	    if(gain == 0){
	      rawVal = (Float_t) (stream.GetADCValue());
    	      pedSubVal = (Float_t) (rawVal-meanPed[pedindex]); 
  	      zpcSignal += pedSubVal; 
  	      isZPCFired = kTRUE;
  	      if(quad!=0) zpcSumQ += pedSubVal;
  	      else{
  		zpcpmC = pedSubVal;
  		FillRawsData(8,zpcpmC);
  	      }
	      indZPC++;
	      
	      FillRawsData(16, pedindex, pedSubVal);
  	    }
  	  }
  	  else if(det == 3){ 
  	    pedindex = quad+9;
  	    if(quad==1){     
  	      if(gain == 0){
	        rawVal = (Float_t) (stream.GetADCValue());
    	        pedSubVal = (Float_t) (rawVal-meanPed[pedindex]); 
  		FillRawsData(4,pedSubVal);
	        FillRawsData(16,pedindex, pedSubVal);
  	      }
  	    }
  	    else if(quad==2){ 
  	      if(gain == 0){
	        rawVal = (Float_t) (stream.GetADCValue());
    	        pedSubVal = (Float_t) (rawVal-meanPed[pedindex]); 
  		FillRawsData(5,pedSubVal); 
	        FillRawsData(16,pedindex, pedSubVal);
  	      }
  	    }
  	  }
  	  else if(det == 4){	   
  	    pedindex = quad+12;
  	    if(gain == 0){
	      rawVal = (Float_t) (stream.GetADCValue());
    	      pedSubVal = (Float_t) (rawVal-meanPed[pedindex]); 
  	      znaSignal  += pedSubVal; 
  	      isZNAFired = kTRUE;
  	      if(quad!=0){
	        znaSumQ += pedSubVal;
	        if(pedSubVal>0.&& znapmC>7.) {
	          wZNA = TMath::Power(pedSubVal, alpha);
	          numXZNA += x[quad-1]*wZNA;
	          numYZNA += y[quad-1]*wZNA;
	          denZNA += wZNA;
	        }
	      }
  	      else{
  		znapmC = pedSubVal;
  		FillRawsData(7,znapmC);
  	      }
	      indZNA++;
	      
	      FillRawsData(16,pedindex, pedSubVal);
  	    }
  	  }
  	  else if(det == 5){
  	    pedindex = quad+17;
  	    if(gain == 0){
	      rawVal = (Float_t) (stream.GetADCValue());
    	      pedSubVal = (Float_t) (rawVal-meanPed[pedindex]); 
  	      zpaSignal  += pedSubVal; 
  	      isZPAFired = kTRUE;
  	      if(quad!=0) zpaSumQ += pedSubVal;
  	      else{
  		zpapmC = pedSubVal;
  		FillRawsData(9,zpapmC);
  	      }
	      indZPA++;
	      
	      FillRawsData(16,pedindex, pedSubVal);
  	    }
  	  }
                	 
	 }

  	 if(isZNCFired && indZNC==5){
  	   FillRawsData(0,zncSignal);
  	   FillRawsData(10,zncSumQ); 
           //
	   Float_t xZNC, yZNC; 	         
	   if(denZNC!=0){
	     xZNC = numXZNC/denZNC;
	     yZNC = numYZNC/denZNC;
	   } 
	   else xZNC = yZNC = 999.;
	   FillRawsData(20,xZNC, yZNC);
  	 }
  	 if(isZPCFired && indZPC==5){
  	   FillRawsData(2,zpcSignal);
           FillRawsData(12,zpcSumQ); 
         }
         if(isZNAFired && indZNA==5){ 
           FillRawsData(1,znaSignal);
           FillRawsData(11,znaSumQ); 
	   //
	   Float_t xZNA, yZNA;
	   if(denZNA!=0){
	     xZNA = numXZNA/denZNA;
	     yZNA = numYZNA/denZNA;
	   } 
	   else xZNA = yZNA = 999.;
	   FillRawsData(21,xZNA, yZNA);
         }
         if(isZPAFired && indZPA==5){ 
           FillRawsData(3,zpaSignal);
           FillRawsData(13,zpaSumQ); 
         }
	 
	 if(indZNC==5){
	   zncSignal = zncSumQ = zncpmC = 0;
	   isZNCFired=kFALSE; indZNC=0;
	 }
	 if(indZPC==5){
	   zpcSignal = zpcSumQ = zpcpmC = 0;
	   isZPCFired=kFALSE; indZPC=0;
	 }
	 if(indZNA==5){
	   znaSignal = znaSumQ = znapmC = 0;
	   isZNAFired=kFALSE; indZNA=0;
	 }
	 if(indZPA==5){
	   zpaSignal = zpaSumQ = zpapmC = 0;
	   isZPAFired=kFALSE; indZPA=0;
	 }
	 
      } //IsADCDataWord && signal ADCs
      else if(stream.IsZDCTDCDatum()){
	 if(stream.GetChannel()==1){
	    zncTDC[iMultZNCTDC] = (0.025*stream.GetZDCTDCDatum());
	    iMultZNCTDC++;
	 }
	 else if(stream.GetChannel()==3){
	    zpcTDC[iMultZPCTDC] = (0.025*stream.GetZDCTDCDatum());
	    iMultZPCTDC++;
	 }
	 else if(stream.GetChannel()==5){
	    znaTDC[iMultZNATDC] = (0.025*stream.GetZDCTDCDatum());
	    iMultZNATDC++;
	 }
	 else if(stream.GetChannel()==7){
	    zpaTDC[iMultZPATDC] = (0.025*stream.GetZDCTDCDatum());
	    iMultZPATDC++;
	 }
	 else if(stream.GetChannel()==8){
	    zem1TDC[iMultZEM1TDC] = (0.025*stream.GetZDCTDCDatum());
	    iMultZEM1TDC++;
	 }
	 else if(stream.GetChannel()==9){
	    zem2TDC[iMultZEM2TDC] = (0.025*stream.GetZDCTDCDatum());
	    iMultZEM2TDC++;
	 }
	 else if(stream.GetChannel()==10){
	    zncSumTDC[iMultZEM2TDC] = (0.025*stream.GetZDCTDCDatum());
	    iMultTDCC++;
	 }
	 else if(stream.GetChannel()==12){
	    znaSumTDC[iMultZEM2TDC] = (0.025*stream.GetZDCTDCDatum());
	    iMultTDCA++;
	 }
	 else if(stream.GetChannel()==14) tdcGate = (0.025*stream.GetZDCTDCDatum());
	 else if(stream.GetChannel()==15) l0 = (0.025*stream.GetZDCTDCDatum());
	 
	 if(stream.GetChannel()==16 && tdcGate!=-999.){
	   for(Int_t iHit=0; iHit<10; iHit++){
	      if(zncTDC[iHit]!=-999.){
	        if(zncTDC[iHit]-tdcGate>-340. && zncTDC[iHit]-tdcGate<-300.) 
		   FillRawsData(17,1, zncTDC[iHit]-tdcGate);
	      }
	      if(zpcTDC[iHit]!=-999.){
	        Float_t diffZPC = zpcTDC[iHit]-tdcGate;
	        FillRawsData(15,diffZPC);
	        if(diffZPC>-340. && diffZPC<-300.) FillRawsData(17,2, diffZPC);
	      }
	      if(znaTDC[iHit]!=-999.){
	        if(znaTDC[iHit]-tdcGate>-340. && znaTDC[iHit]-tdcGate<-300.) 
	          FillRawsData(17,3, znaTDC[iHit]-tdcGate);
	      }
	      if(zpaTDC[iHit]!=-999.){
	        if(zpaTDC[iHit]-tdcGate>-340. && zpaTDC[iHit]-tdcGate<-300.) 
	          FillRawsData(17,4, zpaTDC[iHit]-tdcGate);
	      }
	      if(zem1TDC[iHit]!=-999.){
	        Float_t diffZEM1 = zem1TDC[iHit]-tdcGate;
	        FillRawsData(14,diffZEM1);
		if(diffZEM1>-340. && diffZEM1<-300.) FillRawsData(17,5, diffZEM1);
	      }
	      if(zem2TDC[iHit]!=-999.){
	        if(zem2TDC[iHit]-tdcGate>-340. && zem2TDC[iHit]-tdcGate<-300.) 
	          FillRawsData(17,6, zem2TDC[iHit]-tdcGate);
              }
	      if(zncSumTDC[iHit]!=-999.){
	         Float_t tdcC = zncSumTDC[iHit]-l0;
		 if(znaSumTDC[iHit]!=-999.){
		    Float_t tdcA = znaSumTDC[iHit]-l0;
	            //if (((tdcC-tdcA-refDelta)*(tdcC-tdcA-refDelta)/(sigmaDelta*sigmaDelta) +
	                //(tdcC+tdcA-refSum)*(tdcC+tdcA-refSum)/(sigmaSum*sigmaSum))< 1.0)
			FillRawsData(22,tdcC-tdcA,tdcC+tdcA);
		    
		 }
              }
	   }
	   //
	   tdcGate = -999.;
           for(Int_t i=0; i<10; i++){
              zncTDC[i] = zpcTDC[i] = zem1TDC[i] = zem2TDC[i] = znaTDC[i] = zpaTDC[i] = -999.;
	      zncSumTDC[i] = znaSumTDC[i] = -999.;
           } 
	 }
      }
    
    } //stream.Next()
//  } // check on event type
//  else{
//    AliDebug(1,Form("Skipping non-physics event for QA -> event type %d \n", rawReader->GetType())); 
//  }
//
    IncEvCountCycleRaws();
    IncEvCountTotalRaws();
    //
}

//____________________________________________________________________________
void AliZDCQADataMakerRec::MakeRecPoints(TTree * clustersTree)
{
  // Filling QA histos from RecPoints

  TBranch *branch = clustersTree->GetBranch("ZDC");
  if(!branch){ 
    AliError("Can't get the ZDC branch for rec points!");
    return;
  }
  
  if(!GetRecPointsData(0)) InitRecPoints() ;

  Float_t sum_ZNC=0., sum_ZNA=0., sum_ZPC=0., sum_ZPA=0.;
  Float_t sumQ_ZNC=0., sumQ_ZNA=0., sumQ_ZPC=0., sumQ_ZPA=0.;

  AliZDCReco reco;
  AliZDCReco* preco = &reco;
  clustersTree->SetBranchAddress("ZDC", &preco);

  clustersTree->GetEntry(0);
  for(Int_t i=0; i<5; i++){
    sum_ZNC += reco.GetZN1HREnTow(i);
    sum_ZPC += reco.GetZN2HREnTow(i);
    sum_ZNA += reco.GetZP1HREnTow(i);
    sum_ZPA += reco.GetZP2HREnTow(i);
    if(i==0){
      FillRecPointsData(8,reco.GetZN1HREnTow(i));
      FillRecPointsData(9,reco.GetZN2HREnTow(i));
      FillRecPointsData(10,reco.GetZP1HREnTow(i));
      FillRecPointsData(11,reco.GetZP2HREnTow(i));
    }
    else{
      sumQ_ZNC += reco.GetZN1HREnTow(i);
      sumQ_ZPC += reco.GetZN2HREnTow(i);
      sumQ_ZNA += reco.GetZP1HREnTow(i);
      sumQ_ZPA += reco.GetZP2HREnTow(i);
    }
  }
  
  FillRecPointsData(0,sum_ZNC);
  FillRecPointsData(1,sum_ZNA);
  FillRecPointsData(2,sum_ZPC);
  FillRecPointsData(3,sum_ZPA);
  //
  FillRecPointsData(4,sumQ_ZNC);
  FillRecPointsData(5,sumQ_ZNA);
  FillRecPointsData(6,sumQ_ZPC);
  FillRecPointsData(7,sumQ_ZPA);
  //
  IncEvCountCycleRecPoints();
  IncEvCountTotalRecPoints();
  //  
}  

//____________________________________________________________________________
void AliZDCQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  // make QA data from ESDs
  //
  
  // Check id histograms already created for this Event Specie
  if(!GetESDsData(0)) InitESDs() ;

  AliESDZDC * zdcESD =  esd->GetESDZDC();
  //
  /*TString beamType = esd->GetBeamType();
  Double_t centr_ZNC[2]={999.,999}, centr_ZNA[2]={999.,999};
  if(((beamType.CompareTo("pp"))==0) || ((beamType.CompareTo("p-p"))==0)
     ||((beamType.CompareTo("PP"))==0) || ((beamType.CompareTo("P-P"))==0)){
    zdcESD->GetZNCentroidInpp(centr_ZNC, centr_ZNA);
  }
  else if((beamType.CompareTo("A-A")) == 0 || (beamType.CompareTo("Pb-Pb")) == 0){
    Float_t beamEne = esd->GetBeamEnergy();
    zdcESD->GetZNCentroidInPbPb(beamEne, centr_ZNC, centr_ZNA);
  }
  else printf("\n WARNING!!! AliZDCQADataMakerRec::MakeESDs: can't calculate centroids for beam type: %s\n\n",beamType.Data());
  FillESDsData(0,centr_ZNC[0], centr_ZNC[1]);
  FillESDsData(1,centr_ZNA[0], centr_ZNA[1]);*/

  FillESDsData(0,esd->GetZDCN1Energy());
  FillESDsData(1,esd->GetZDCN2Energy());
  FillESDsData(2,esd->GetZDCP1Energy());
  FillESDsData(3,esd->GetZDCP2Energy());
  FillESDsData(4,esd->GetZDCEMEnergy(0));
  FillESDsData(5,esd->GetZDCEMEnergy(1));
  //
  Double_t sumQZNC=0., sumQZPC=0., sumQZNA=0., sumQZPA=0.;
  //Double_t sumQZNC_lg=0., sumQZPC_lg=0., sumQZNA_lg=0., sumQZPA_lg=0.;
  //
  const Double_t *towZNC, *towZPC, *towZNA, *towZPA;
  //const Double_t *towZNC_lg, *towZPC_lg, *towZNA_lg, *towZPA_lg;
  //
  towZNC = zdcESD->GetZN1TowerEnergy();
  towZPC = zdcESD->GetZP1TowerEnergy();
  towZNA = zdcESD->GetZN2TowerEnergy();
  towZPA = zdcESD->GetZP2TowerEnergy();
  //
  /*towZNC_lg = zdcESD->GetZN1TowerEnergyLR();
  towZPC_lg = zdcESD->GetZP1TowerEnergyLR();
  towZNA_lg = zdcESD->GetZN2TowerEnergyLR();
  towZPA_lg = zdcESD->GetZP2TowerEnergyLR();*/
  //
  for(Int_t i=0; i<5; i++){
     if(i==0){
       FillESDsData(10,towZNC[i]);
       FillESDsData(11,towZNA[i]);
       FillESDsData(12,towZPC[i]);
       FillESDsData(13,towZPA[i]);
     }
     else{
       sumQZNC += towZNC[i];
       sumQZPC += towZPC[i];
       sumQZNA += towZNA[i];
       sumQZPA += towZPA[i];
     }
  }
  FillESDsData(6,sumQZNC);
  FillESDsData(7,sumQZNA);
  FillESDsData(8,sumQZPC);
  FillESDsData(9,sumQZPA);
  //
  IncEvCountCycleESDs();
  IncEvCountTotalESDs();
  //
}

//____________________________________________________________________________
void AliZDCQADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle

  fPedCalibData = GetPedCalibData();
  
}

//____________________________________________________________________________ 
void AliZDCQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  //
  ResetEventTrigClasses();
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {  // RS: loop over event types
    //
    if (!IsValidEventSpecie(specie, list)) continue;
    SetEventSpecie(AliRecoParam::ConvertIndex(specie));
    //
    for (int itc=-1;itc<GetNTrigClasses();itc++) { // RS: loop over eventual clones per trigger class
      //
      if( task == AliQAv1::kRAWS) {
	TProfile* h16 = dynamic_cast<TProfile*> (GetRawsData(16, itc));
	TProfile* h17 =  dynamic_cast<TProfile*> (GetRawsData(17, itc));
	TH1F* h18 =  dynamic_cast<TH1F*> (GetRawsData(18, itc));
	TH1F* h19 =  dynamic_cast<TH1F*> (GetRawsData(19, itc));
	TH2F* h20 =  dynamic_cast<TH2F*> (GetRawsData(20, itc));
	TH2F* h21 =  dynamic_cast<TH2F*> (GetRawsData(21, itc));
	TH2F* h22 =  dynamic_cast<TH2F*> (GetRawsData(22, itc));
	if (!h16 || !h17 || !h18 || !h19){
	 AliWarning("AliZDCQADataMakerRec -> RAW histos 16||17||18||19 not found!"); 
	 AliWarning(Form("for specie %s and trigger class %s",
			 AliRecoParam::GetEventSpecieName(specie), AliQADataMaker::GetTrigClassName(itc)));
	}
	else{
	  //h16->Draw("");
	  for(Int_t ibin=1; ibin<=h16->GetNbinsX(); ibin++){
	    h18->SetBinContent(ibin, h16->GetBinContent(ibin)); 
	    h18->SetBinError(ibin, h16->GetBinError(ibin));
	  }
	  for(Int_t ibin=1; ibin<=h17->GetNbinsX(); ibin++){
	    h19->SetBinContent(ibin, h17->GetBinContent(ibin)); 
	    h19->SetBinError(ibin, h17->GetBinError(ibin));
	  }
	  h18->SetLineColor(kBlue); h18->SetLineWidth(2);
	  h19->SetLineColor(kAzure-3); h19->SetLineWidth(2);
        }
	if(!h20){
	 AliWarning("AliZDCQADataMakerRec -> RAW histos 20||21 not found!"); 
	 AliWarning(Form("for specie %s and trigger class %s",
			 AliRecoParam::GetEventSpecieName(specie), AliQADataMaker::GetTrigClassName(itc)));
	}
	else{
	 h20->SetMarkerColor(kPink+7); 
	 h21->SetMarkerColor(kBlue+2);
	}
	if(!h22) {
	 AliWarning("AliZDCQADataMakerRec -> RAW histo 22 not found!"); 
	 AliWarning(Form("for specie %s and trigger class %s",
			 AliRecoParam::GetEventSpecieName(specie), AliQADataMaker::GetTrigClassName(itc)));
	}
	else h22->SetMarkerColor(kAzure+7);
      }
    } // loop over t
  } //  loop over species
  	
  AliQAChecker::Instance()->Run(AliQAv1::kZDC, task, list) ;  
}
