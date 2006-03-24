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

/* $Id$ */

#include <Riostream.h>

#include <TDirectory.h>

#include "AliRunLoader.h"
#include "AliRun.h"
#include <AliESD.h>
#include "AliLog.h"
#include <TClonesArray.h>
#include "AliSTARTRecPoint.h"
#include "AliRawReader.h"
#include "AliSTARTRawReader.h"
#include "AliSTARTLoader.h"
#include "AliSTARTdigit.h"
#include "AliSTARTReconstructor.h"
#include "AliSTARTParameters.h"
#include "AliSTARTAlignData.h"
#include "AliSTARTCalibData.h"
#include "AliCDBLocal.h"
#include "AliCDBStorage.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"

#include <TArrayI.h>

ClassImp(AliSTARTReconstructor)
AliSTARTAlignData* AliSTARTReconstructor::fgAlignData = 0;
AliSTARTCalibData* AliSTARTReconstructor::fgCalibData = 0;

  void  AliSTARTReconstructor::ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const
{
  //START raw data-> digits conversion
 // reconstruct time information from raw data
  AliSTARTRawReader myrawreader(rawReader,digitsTree);
  myrawreader.NextThing();
}
  void AliSTARTReconstructor::Reconstruct(TTree*digitsTree, TTree*clustersTree) const
{
// START digits reconstruction
// STARTRecPoint writing 

    //Q->T-> coefficients !!!! should be asked!!!
  //  Float_t ph2MIP=500;
  Float_t gain[24], timeDelayCFD[24], timeDelayLED[24];
  Float_t zdetA,zdetC;
  TObjArray slewingLED;
    
  TArrayI * fADC = new TArrayI(24); 
  TArrayI * fTimeCFD = new TArrayI(24); 
  TArrayI * fADCLED = new TArrayI(24); 
  TArrayI * fTimeLED = new TArrayI(24); 
  cout<<" fTimeCFD "<<fTimeCFD<<endl;

  AliSTARTParameters* param = AliSTARTParameters::Instance();
  param->Init();
  Int_t ph2MIP = param->GetPh2Mip();     
  Int_t channelWidth = param->GetChannelWidth() ;  
  
    for (Int_t i=0; i<24; i++){
      timeDelayCFD[i] = param->GetTimeDelayCFD(i);
      timeDelayLED[i] = param->GetTimeDelayLED(i);
      gain[i] = param->GetGain(i);
      slewingLED.AddAtAndExpand(param->GetSlew(i),i);
     }
    zdetC = param->GetZposition(0);
    zdetA  = param->GetZposition(1);
  
  AliDebug(1,Form("Start DIGITS reconstruction "));
 
  TBranch *brDigits=digitsTree->GetBranch("START");
  AliSTARTdigit *fDigits = new AliSTARTdigit();
  if (brDigits) {
    brDigits->SetAddress(&fDigits);
  }else{
    cerr<<"EXEC Branch START digits not found"<<endl;
    return;
  }
  brDigits->GetEntry(0);
  fDigits->GetTime(*fTimeCFD);
  fDigits->GetADC(*fADC);
  fDigits->GetTimeAmp(*fTimeLED);
  fDigits->GetADCAmp(*fADCLED);

  Float_t time[24], adc[24];
  for (Int_t ipmt=0; ipmt<24; ipmt++)
    {
      
         if(fTimeCFD->At(ipmt)>0 ){
	 time[ipmt] = channelWidth *( fTimeCFD->At(ipmt)) - 1000*timeDelayCFD[ipmt];
	 cout<<ipmt<<" "<<time[ipmt];
	 Float_t adc_digPs = channelWidth * Float_t (fADC->At(ipmt)) ;
	 //	 cout<<"  adc_digmV "<< adc_digPs<<endl;
	 adc[ipmt] = TMath::Exp(adc_digPs/1000) /gain[ipmt];
	 //	 cout<<" adc"<<adc[ipmt]<<" inMIP "<<adc[ipmt]/50<< endl;
         }
    }

  Int_t besttimeright=channelWidth * (fDigits->BestTimeRight());
  Int_t besttimeleft=channelWidth * (fDigits->BestTimeLeft());
  //folding with experimental time distribution
  //  Float_t c = 29.9792; // cm/ns
  Float_t c = 0.0299792; // cm/ps
  Float_t lenr=TMath::Sqrt(zdetA*zdetA + 6.5*6.5);
  Float_t lenl=TMath::Sqrt(zdetC*zdetC + 6.5*6.5);
  Float_t timeDiff=channelWidth * (fDigits->TimeDiff());
  Int_t meanTime=channelWidth * (fDigits->MeanTime());
  Float_t ds=(c*(timeDiff)-(lenr-lenl))/2;
  AliDebug(2,Form(" timediff in ns %f  real point%f",timeDiff,ds));
  
  /*
  fDigits->GetSumMult(*fSumMult);
  Int_t multipl[4]; 
 
 for (Int_t i=0; i<4; i++)
    {
      Float_t  mult=Float_t (fSumMult->At(i));
      Float_t   realMultmV=TMath::Exp(mult/mV2channel);
      multipl[i]=Int_t ((realMultmV/ph2mV)/500+0.5);
    }
  */

  //  AliDebug(2,Form(" multiplicity Abs side %i  multiplicity non-Abs side %i",multipl[1],multipl[2]));

  AliSTARTRecPoint* frecpoints= new AliSTARTRecPoint ();
  clustersTree->Branch( "START", "AliSTARTRecPoint" ,&frecpoints, 405,1);
  frecpoints->SetTimeBestRight(besttimeright);
  frecpoints->SetTimeBestLeft(besttimeleft);
  frecpoints->SetVertex(ds);
  frecpoints->SetMeanTime(meanTime);
  /*
  frecpoints->SetMult(multipl[0]);
  frecpoints->SetMultA(multipl[2]);
  frecpoints->SetMultC(multipl[1]);
  */
  clustersTree->Fill();
}


void AliSTARTReconstructor::FillESD(AliRunLoader* runLoader, AliESD *pESD) const
{

  /***************************************************
  Resonstruct digits to vertex position
  ****************************************************/
  
  //  Float_t c = 0.3;  //speed of light mm/ps
  // Float_t Zposition=0;
  
  if (!runLoader) {
    AliError("Reconstruct >> No run loader");
    return;
  }
  
  AliDebug(1,Form("Start FillESD START"));

  AliSTARTLoader* pStartLoader = (AliSTARTLoader*) runLoader->GetLoader("STARTLoader");
 
  pStartLoader->LoadRecPoints("READ");
  
    TTree *treeR = pStartLoader->TreeR();
  
   AliSTARTRecPoint* frecpoints= new AliSTARTRecPoint ();
    if (!frecpoints) {
    AliError("Reconstruct Fill ESD >> no recpoints found");
    return;
  }
  
  AliDebug(1,Form("Start FillESD START"));
   TBranch *brRec = treeR->GetBranch("START");
    if (brRec) {
      brRec->SetAddress(&frecpoints);
    }else{
      cerr<<"EXEC Branch START rec not found"<<endl;
      exit(111);
    } 
 
    brRec->GetEntry(0);
    Float_t Zposition=frecpoints->GetVertex();
    
    pESD->SetT0zVertex(Zposition);
    pStartLoader->UnloadRecPoints();
    
} // vertex in 3 sigma






