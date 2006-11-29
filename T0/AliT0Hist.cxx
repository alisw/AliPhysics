#include "AliT0RawReader.h"
#include "AliT0RawData.h"
#include "/home/alla/AliRoot/verynew/RAW/AliRawReaderFile.h"


#include <Riostream.h>
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TArrayI.h"
#include "AliT0Hist.h"  
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "/home/alla/AliRoot/verynew/RAW/AliRawReader.h"
#include "/home/alla/AliRoot/verynew/RAW/AliRawReaderFile.h"
#include "AliT0RawReader.h"
ClassImp(AliT0Hist)

  AliT0Hist::AliT0Hist() : TObject()
{
  fhADCR = new TH1F("hADCR","ADC right",100,1.,500);
  fhADCL = new TH1F("hADCL","ADC right",100,1.,500);
  fhTimeR = new TH1F("hTimeR","Time Right",100,400,1000);
  fhTimeL = new TH1F("hTimeL","Time Left",100,10,500);
  fhBestTimeR = new TH1F("hBestTimeR","Time Right",100,400,1000);
  fhBestTimeL = new TH1F("hBestTimeL","Time Left",100,0,500);
  fhADCDet = new TH1F("hADCDet","ADC vs Ndet",30,0.,30);
  fhTimeDet = new TH1F("hTimeDet","Time vs Ndet",30,0,30);
  fhTimeDiff = new TH1F("hTimeDiff"," Time diff",100,350,450);
  fhMeanTime = new TH1F("hMeanTime"," Mean Time ",100,250,350);
  fhT0detR = new TH1F("hT0detR"," T0 vs Ndet right ",15,0,15);
  fhT0detL = new TH1F("hT0detL"," T0 vs Ndet left ",15,0,15);
  fhTimevsADCL = new TH2F("hTimevsADCL","time vs ADC left",100,0,500,100,1,500);
  fhTimevsADCR = new TH2F("hTimevsADCR","time vs ADC right",100,0,500,100,400,1000);
}


  AliT0Hist::~AliT0Hist()
{
  /*
  delete fhADCR;
  delete fhADCL;
  delete fhTimeR;
  delete fhTimeL;
  delete fhBestTimeR;
  delete fhBestTimeL;
  delete fhADCDet;
  delete fhTimeDet;
  delete fhTimeDiff;
  delete fhMeanTime;
  delete fhT0detR;
  delete fhT0detL;
  delete fhTimevsADCL;
  */
}

void AliT0Hist::FillHist(AliRunLoader* runLoader, 
			    AliRawReader* rawReader) const 
{

   TArrayI *fTime= new TArrayI(24);
   TArrayI *fADC = new TArrayI(24);
   AliT0RawReader myrawreader;
   

   
  Int_t iEvent = 0;
  while (rawReader->NextEvent()) {
      runLoader->GetEvent(iEvent++);
      
      myrawreader.NextThing( rawReader);
      myrawreader.GetTime(*fTime);
      myrawreader.GetADC(*fADC);
      
      Int_t besttimeleft=99999;
      Int_t besttimeright=99999;
      Int_t ibestR=999; Int_t ibestL=999;
      for (Int_t i=0; i<12; i++ )
	{
	  Int_t timel=fTime->At(i);
	  if(timel<besttimeleft && timel!=0) { besttimeleft=timel; ibestL=i;}
	  fhTimeL->Fill(timel);
	  fhTimeDet->Fill(i,timel);
	  Int_t adcl=fADC->At(i);
	  fhADCL->Fill(adcl);
	  fhADCDet->Fill(i,adcl);
	  Int_t timer=fTime->At(i+12);
	  if(timer<besttimeright && timer!=0){ besttimeright=timer; ibestR=i;}
	  fhTimeDet->Fill(i+12,timer);
	  fhTimeR->Fill(timer);
	  Int_t adcr=fADC->At(i+12);
	  fhADCR->Fill(adcr);
	  fhADCDet->Fill(i+12,adcr);
	  fhTimevsADCL->Fill(adcl,timel);
	  fhTimevsADCR->Fill(adcr,timer);

	}
      
      fhBestTimeR->Fill(besttimeright);
      fhBestTimeL->Fill(besttimeleft);
      Int_t timeDiff1=besttimeright-besttimeleft;
      Float_t meanTime=(besttimeright+besttimeleft)/2.;
      fhTimeDiff->Fill(timeDiff1);
      fhMeanTime->Fill(meanTime);
      //     Float_t  timeDiff2=((besttimeright-besttimeleft)-(350.-69.7))/2;
      fhT0detR->Fill(ibestR); fhT0detL->Fill(ibestL);
  }
}



