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
#include <AliESD.h>
#include "AliLog.h"
#include "AliT0Loader.h"
#include <TClonesArray.h>
#include "AliT0RecPoint.h"
#include "AliRawReader.h"
#include "AliT0RawReader.h"
#include "AliT0digit.h"
#include "AliT0Reconstructor.h"
#include "AliT0Parameters.h"
#include "AliT0Calibrator.h"
#include "AliCDBLocal.h"
#include "AliCDBStorage.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"

#include <TArrayI.h>
#include <TGraph.h>

ClassImp(AliT0Reconstructor)

  AliT0Reconstructor:: AliT0Reconstructor(): AliReconstructor(),
			fDigits(NULL),
			fTree(0x0),
			fZposition(0)

{
 AliDebug(1,"Start reconstructor ");
}
//____________________________________________________________________

AliT0Reconstructor::AliT0Reconstructor(const AliT0Reconstructor &r):
  fDigits(NULL),
  fTree(0x0),
  fZposition(0)
  
{
  //
  // AliT0Reconstructor copy constructor
  //

  ((AliT0Reconstructor &) r).Copy(*this);

}

//_____________________________________________________________________________
AliT0Reconstructor &AliT0Reconstructor::operator=(const AliT0Reconstructor &r)
{
  //
  // Assignment operator
  //

  if (this != &r) ((AliT0Reconstructor &) r).Copy(*this);
  return *this;

}

//_____________________________________________________________________________

void AliT0Reconstructor::Reconstruct(TTree*digitsTree, TTree*clustersTree) const

{
// T0 digits reconstruction
// T0RecPoint writing 

  
  Float_t  timeDelayLED[24];
  Float_t zdetA,zdetC;
  TObjArray slewingLEDrec;
  TObjArray walk;
    
  TArrayI * timeCFD = new TArrayI(24); 
  TArrayI * timeLED = new TArrayI(24); 
  TArrayI * chargeQT0 = new TArrayI(24); 
  TArrayI * chargeQT1 = new TArrayI(24); 
  
  AliT0Parameters* param = AliT0Parameters::Instance();
  param->Init();
  AliT0Calibrator *calib=new AliT0Calibrator(); 

  Int_t mV2Mip = param->GetmV2Mip();     
  //mV2Mip = param->GetmV2Mip();     
  Int_t channelWidth = param->GetChannelWidth() ;  
  
  for (Int_t i=0; i<24; i++){
    TGraph* gr = param ->GetSlewRec(i);
    slewingLEDrec.AddAtAndExpand(gr,i) ;  
  }
  zdetC = param->GetZposition(0);
  zdetA  = param->GetZposition(1);
    
  AliDebug(1,Form("Start DIGITS reconstruction "));

  TBranch *brDigits=digitsTree->GetBranch("T0");
  AliT0digit *fDigits = new AliT0digit() ;
  if (brDigits) {
    brDigits->SetAddress(&fDigits);
  }else{
    cerr<<"EXEC Branch T0 digits not found"<<endl;
    return;
  }
  
  digitsTree->GetEvent(0);
  digitsTree->GetEntry(0);
  brDigits->GetEntry(0);
  fDigits->GetTimeCFD(*timeCFD);
  fDigits->GetTimeLED(*timeLED);
  fDigits->GetQT0(*chargeQT0);
  fDigits->GetQT1(*chargeQT1);
  
  Float_t besttimeA=999999;
  Float_t besttimeC=999999;
  Int_t pmtBestA=99999;
  Int_t pmtBestC=99999;
  Float_t timeDiff=999999, meanTime=0;
  


   AliT0RecPoint* frecpoints= new AliT0RecPoint ();
   clustersTree->Branch( "T0", "AliT0RecPoint" ,&frecpoints, 405,1);

  Float_t time[24], adc[24];
  for (Int_t ipmt=0; ipmt<24; ipmt++) {
    if(timeCFD->At(ipmt)>0 ){
      Int_t qt0= chargeQT0->At(ipmt);
      Int_t qt1= chargeQT1->At(ipmt);
      if((qt1-qt0)>0)  adc[ipmt] = TMath::Exp( Double_t (channelWidth*(qt1-qt0)/1000));
      time[ipmt] = channelWidth * (calib-> WalkCorrection( ipmt,qt1 , timeCFD->At(ipmt) ) ) ;
      
      //LED
      Double_t sl = (timeLED->At(ipmt) - timeCFD->At(ipmt)- (1000.*timeDelayLED[ipmt]/channelWidth))*channelWidth;
      Double_t qt=((TGraph*)slewingLEDrec.At(ipmt))->Eval(sl/1000.);
      frecpoints->SetTime(ipmt,time[ipmt]);
      frecpoints->SetAmp(ipmt,adc[ipmt]);
      frecpoints->SetAmpLED(ipmt,qt);
    }
    else {
      time[ipmt] = 0;
      adc[ipmt] = 0;
    }
  }

  for (Int_t ipmt=0; ipmt<12; ipmt++){
    if(time[ipmt] > 1 ) {
      if(time[ipmt]<besttimeC){
	besttimeC=time[ipmt]; //timeC
	pmtBestC=ipmt;
      }
    }
  }
  for ( Int_t ipmt=12; ipmt<24; ipmt++){
    if(time[ipmt] > 1) {
      if(time[ipmt]<besttimeA) {
	besttimeA=time[ipmt]; //timeA
        pmtBestA=ipmt;}
    }
  }
  if(besttimeA !=999999)  frecpoints->SetTimeBestA(Int_t(besttimeA));
  if( besttimeC != 999999 ) frecpoints->SetTimeBestC(Int_t(besttimeC));
  AliDebug(1,Form(" besttimeA %f ps,  besttimeC %f ps",besttimeA, besttimeC));
  Float_t c = 0.0299792; // cm/ps
  Float_t vertex = 0;
  if(besttimeA !=999999 && besttimeC != 999999 ){
    timeDiff = besttimeC - besttimeA;
    meanTime = (besttimeA + besttimeC)/2.;
    vertex = c*(timeDiff)/2.; //-(lenr-lenl))/2;
    AliDebug(1,Form("  timeDiff %f ps,  meanTime %f ps, vertex %f cm",timeDiff, meanTime,vertex ));
    frecpoints->SetVertex(vertex);
    frecpoints->SetMeanTime(Int_t(meanTime));
    
  }
  clustersTree->Fill();
}


//_______________________________________________________________________

void AliT0Reconstructor::Reconstruct(AliRawReader* rawReader, TTree*recTree) const
{
// T0 raw ->
// T0RecPoint writing 

    //Q->T-> coefficients !!!! should be asked!!!
  Float_t  timeDelayLED[24];
  Float_t zdetA,zdetC;
  TObjArray slewingLEDrec;
  TObjArray walk;
    
  TArrayI * timeCFD = new TArrayI(24); 
  TArrayI * timeLED = new TArrayI(24); 
  TArrayI * chargeQT0 = new TArrayI(24); 
  TArrayI * chargeQT1 = new TArrayI(24); 

   AliT0RawReader myrawreader(rawReader);
   if (!myrawreader.Next())
     AliDebug(1,Form(" no raw data found!! %i", myrawreader.Next()));
   Int_t allData[110][5];
   for (Int_t i=0; i<110; i++) {
     allData[i][0]=myrawreader.GetData(i,0);
   }

  AliT0Parameters* param = AliT0Parameters::Instance();
  param->Init();
  AliT0Calibrator *calib=new AliT0Calibrator(); 

  Int_t mV2Mip = param->GetmV2Mip();     
  //mV2Mip = param->GetmV2Mip();     
  Int_t channelWidth = param->GetChannelWidth() ;  
    
  for (Int_t i=0; i<24; i++){
    TGraph* gr = param ->GetSlewRec(i);
    slewingLEDrec.AddAtAndExpand(gr,i) ;  
  }
  
  zdetC = param->GetZposition(0);
  zdetA  = param->GetZposition(1);
    
   for (Int_t in=0; in<24; in++)
     {
       timeLED->AddAt(allData[in+1][0],in);
       timeCFD->AddAt(allData[in+25][0],in);
       chargeQT1->AddAt(allData[in+55][0],in);
       chargeQT0->AddAt(allData[in+79][0],in);
       AliDebug(10, Form(" readed Raw %i %i %i %i %i", in, timeLED->At(in),timeCFD->At(in),chargeQT0->At(in),chargeQT1->At(in)));
     }

  Float_t besttimeA=999999;
  Float_t besttimeC=999999;
  Int_t pmtBestA=99999;
  Int_t pmtBestC=99999;
  Float_t timeDiff=999999, meanTime=0;

  
  AliT0RecPoint* frecpoints= new AliT0RecPoint ();
 
  recTree->Branch( "T0", "AliT0RecPoint" ,&frecpoints, 405,1);
  

  Float_t time[24], adc[24];
  for (Int_t ipmt=0; ipmt<24; ipmt++) {
    if(timeCFD->At(ipmt)>0 ){
      Int_t qt0= chargeQT0->At(ipmt);
      Int_t qt1= chargeQT1->At(ipmt);
      if((qt1-qt0)>0)  adc[ipmt] = TMath::Exp( Double_t (channelWidth*(qt1-qt0)/1000));
      time[ipmt] = channelWidth * (calib-> WalkCorrection( ipmt,qt1 , timeCFD->At(ipmt) ) ) ;
      Double_t sl = (timeLED->At(ipmt) - timeCFD->At(ipmt)- (1000.*timeDelayLED[ipmt]/channelWidth))*channelWidth;
      Double_t qt=((TGraph*)slewingLEDrec.At(ipmt))->Eval(sl/1000.);
      frecpoints->SetTime(ipmt,time[ipmt]);
      frecpoints->SetAmp(ipmt,adc[ipmt]);
      frecpoints->SetAmpLED(ipmt,qt);
    }
    else {
      time[ipmt] = 0;
      adc[ipmt] = 0;
    }
  }

  for (Int_t ipmt=0; ipmt<12; ipmt++){
    if(time[ipmt] > 1 ) {
      if(time[ipmt]<besttimeC){
	besttimeC=time[ipmt]; //timeC
	pmtBestC=ipmt;
      }
    }
  }
  for ( Int_t ipmt=12; ipmt<24; ipmt++){
    if(time[ipmt] > 1) {
      if(time[ipmt]<besttimeA) {
	besttimeA=time[ipmt]; //timeA
        pmtBestA=ipmt;}
    }
  }
  if(besttimeA !=999999)  frecpoints->SetTimeBestA(Int_t(besttimeA));
  if( besttimeC != 999999 ) frecpoints->SetTimeBestC(Int_t(besttimeC));
  AliDebug(1,Form(" besttimeA %f ps,  besttimeC %f ps",besttimeA, besttimeC));
  Float_t c = 0.0299792; // cm/ps
  Float_t vertex = 0;
  if(besttimeA !=999999 && besttimeC != 999999 ){
    timeDiff = besttimeC - besttimeA;
    meanTime = (besttimeA + besttimeC)/2.;
    vertex = c*(timeDiff)/2.; //-(lenr-lenl))/2;
    AliDebug(1,Form("  timeDiff %f ps,  meanTime %f ps, vertex %f cm",timeDiff, meanTime,vertex ));
    frecpoints->SetVertex(vertex);
    frecpoints->SetMeanTime(Int_t(meanTime));
    
  }
  recTree->Fill();
}
//____________________________________________________________

void AliT0Reconstructor::FillESD(AliRunLoader* runLoader, AliESD *pESD) const
{

  /***************************************************
  Resonstruct digits to vertex position
  ****************************************************/
  
  
  if (!runLoader) {
    AliError("Reconstruct >> No run loader");
    return;
  }
  
  AliDebug(1,Form("Start FillESD T0"));

  AliT0Loader* pStartLoader = (AliT0Loader*) runLoader->GetLoader("T0Loader");
 
  pStartLoader->LoadRecPoints("READ");

  TTree *treeR = pStartLoader->TreeR();
  
   AliT0RecPoint* frecpoints= new AliT0RecPoint ();
    if (!frecpoints) {
    AliError("Reconstruct Fill ESD >> no recpoints found");
    return;
  }
  
  AliDebug(1,Form("Start FillESD T0"));
  TBranch *brRec = treeR->GetBranch("T0");
  if (brRec) {
    brRec->SetAddress(&frecpoints);
  }else{
    cerr<<"EXEC Branch T0 rec not found"<<endl;
    // exit(111);
    return;
  } 
    
    brRec->GetEntry(0);
    Float_t timeStart, Zposition, amp[24], time[24];
    Int_t mean0 = 12450;
    Int_t i;
    Zposition = frecpoints -> GetVertex();
    timeStart = frecpoints -> GetMeanTime() - mean0;
    for ( i=0; i<24; i++) {
       time[i] = Float_t (frecpoints -> GetTime(i)) / 1000.; // ps to ns
      amp[i] = frecpoints -> GetAmp(i);
    }
    pESD->SetT0zVertex(Zposition); //vertex Z position 
    pESD->SetT0(timeStart);        // interaction time 
    pESD->SetT0time(time);         // best TOF on each PMT 
    pESD->SetT0amplitude(amp);     // number of particles(MIPs) on each PMT
    cout<<" ESD >> "<<Zposition<<" "<<timeStart<<endl;

    pStartLoader->UnloadRecPoints();
   
} // vertex in 3 sigma






