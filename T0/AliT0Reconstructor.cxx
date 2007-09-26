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
#include <AliESDEvent.h>
#include "AliLog.h"
#include "AliT0Loader.h"
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
#include <TMath.h>
//#include <TGeoManager.h>
//#include <TClonesArray.h>
//#include <TString.h>
//#include <TGeoPNEntry.h>
//#include <TGeoPhysicalNode.h>
//#include  <TGeoHMatrix.h>

ClassImp(AliT0Reconstructor)

  AliT0Reconstructor:: AliT0Reconstructor(): AliReconstructor(),
  fZposition(0),
  fParam(NULL),
  fAmpLEDrec(),
  fdZ_A(0),
  fdZ_C(0)
{
 AliDebug(1,"Start reconstructor ");
  
  fParam = AliT0Parameters::Instance();
  fParam->Init();
  for (Int_t i=0; i<24; i++){
    TGraph* gr = fParam ->GetAmpLEDRec(i);
    fAmpLEDrec.AddAtAndExpand(gr,i) ;  
    fTime0vertex[i]= fParam->GetTimeDelayDA(i);
  }
  fdZ_C = TMath::Abs(fParam->GetZPositionShift("T0/C/PMT1"));
  fdZ_A = TMath::Abs(fParam->GetZPositionShift("T0/A/PMT15"));

}
//____________________________________________________________________

AliT0Reconstructor::AliT0Reconstructor(const AliT0Reconstructor &r):
  AliReconstructor(r),
  fZposition(0),
  fParam(NULL),
  fAmpLEDrec(),
  fdZ_A(0),
  fdZ_C(0)
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
  
   
  TArrayI * timeCFD = new TArrayI(24); 
  TArrayI * timeLED = new TArrayI(24); 
  TArrayI * chargeQT0 = new TArrayI(24); 
  TArrayI * chargeQT1 = new TArrayI(24); 

  AliT0Calibrator *calib=new AliT0Calibrator(); 

  //  Int_t mV2Mip = param->GetmV2Mip();     
  //mV2Mip = param->GetmV2Mip();     
  Int_t channelWidth = fParam->GetChannelWidth() ;  
  Int_t meanT0 = fParam->GetMeanT0();
  
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
      Double_t qt0 = Double_t(chargeQT0->At(ipmt));
      Double_t qt1 = Double_t(chargeQT1->At(ipmt));
      if((qt1-qt0)>0)  adc[ipmt] = TMath::Exp( Double_t (channelWidth*(qt1-qt0)/1000));
      time[ipmt] = calib-> WalkCorrection( ipmt,Int_t(qt1) , timeCFD->At(ipmt) ) ;
      
      //LED
      Double_t sl = (timeLED->At(ipmt) - time[ipmt])*channelWidth;
      Double_t qt=((TGraph*)fAmpLEDrec.At(ipmt))->Eval(sl/1000.);
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
    timeDiff =(besttimeC - besttimeA)*channelWidth;
    meanTime = (meanT0 - (besttimeA + besttimeC)/2) * channelWidth;
    vertex = c*(timeDiff)/2. + (fdZ_A - fdZ_C)/2; //-(lenr-lenl))/2;
    AliDebug(1,Form("  timeDiff %f ps,  meanTime %f ps, vertex %f cm",timeDiff, meanTime,vertex ));
    frecpoints->SetVertex(vertex);
    frecpoints->SetMeanTime(Int_t(meanTime));
    
  }
  //time in each channel as time[ipmt]-MeanTimeinThisChannel(with vertex=0)
  for (Int_t ipmt=0; ipmt<24; ipmt++) {
    if(time[ipmt]>1) {
      time[ipmt] = (time[ipmt] - fTime0vertex[ipmt])*channelWidth;
      frecpoints->SetTime(ipmt,time[ipmt]);
    }
  }
  clustersTree->Fill();

  delete timeCFD;
  delete timeLED;
  delete chargeQT0; 
  delete chargeQT1; 
}


//_______________________________________________________________________

void AliT0Reconstructor::Reconstruct(AliRawReader* rawReader, TTree*recTree) const
{
  // T0 raw ->
  // T0RecPoint writing 
  
  //Q->T-> coefficients !!!! should be measured!!!
  Int_t allData[110][5];
  
  TArrayI * timeCFD = new TArrayI(24); 
  TArrayI * timeLED = new TArrayI(24); 
  TArrayI * chargeQT0 = new TArrayI(24); 
  TArrayI * chargeQT1 = new TArrayI(24); 

   for (Int_t i=0; i<110; i++) {
     allData[i][0]=0;
   }

   AliT0RawReader myrawreader(rawReader);
   if (!myrawreader.Next())
     AliDebug(1,Form(" no raw data found!! %i", myrawreader.Next()));
   for (Int_t i=0; i<110; i++) {
     allData[i][0]=myrawreader.GetData(i,0);
   }
  AliT0Calibrator *calib = new AliT0Calibrator(); 

  //  Int_t mV2Mip = param->GetmV2Mip();     
  //mV2Mip = param->GetmV2Mip();     
  Int_t channelWidth = fParam->GetChannelWidth() ;  
  Int_t meanT0 = fParam->GetMeanT0();
 
  for (Int_t in=0; in<24; in++)
     {
       timeLED->AddAt(allData[in+1][0],in);
       timeCFD->AddAt(allData[in+25][0],in);
       chargeQT1->AddAt(allData[in+57][0],in);
       chargeQT0->AddAt(allData[in+80][0],in);
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
      Double_t qt0 = Double_t(chargeQT0->At(ipmt));
      Double_t qt1 = Double_t(chargeQT1->At(ipmt));
      if((qt1-qt0)>0)  adc[ipmt] = TMath::Exp( Double_t (channelWidth*(qt1-qt0)/1000));
      //      time[ipmt] = channelWidth * (calib-> WalkCorrection( ipmt,qt1 , timeCFD->At(ipmt) ) ) ;
      time[ipmt] = calib-> WalkCorrection( ipmt,Int_t(qt1) , timeCFD->At(ipmt) ) ;
      Double_t sl = (timeLED->At(ipmt) - time[ipmt])*channelWidth;
      Double_t qt=((TGraph*)fAmpLEDrec.At(ipmt))->Eval(sl/1000.);
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
    timeDiff = (besttimeC - besttimeA)*channelWidth;
    //    meanTime = (besttimeA + besttimeC)/2.;
    meanTime = (meanT0 - (besttimeA + besttimeC)/2) * channelWidth;
    vertex = c*(timeDiff)/2. + (fdZ_A - fdZ_C)/2; //-(lenr-lenl))/2;
    AliDebug(1,Form("  timeDiff %f ps,  meanTime %f ps, vertex %f cm",timeDiff, meanTime,vertex ));
    frecpoints->SetVertex(vertex);
    frecpoints->SetMeanTime(Int_t(meanTime));
    
  }
  //time in each channel as time[ipmt]-MeanTimeinThisChannel(with vertex=0)
  for (Int_t ipmt=0; ipmt<24; ipmt++) {
    if(time[ipmt]>1) {
      time[ipmt] = (time[ipmt] - fTime0vertex[ipmt])*channelWidth;
      frecpoints->SetTime(ipmt,time[ipmt]);
    }
  }
  recTree->Fill();
 

  delete timeCFD;
  delete timeLED;
  delete chargeQT0; 
  delete chargeQT1; 
  
}
//____________________________________________________________

void AliT0Reconstructor::FillESD(TTree */*digitsTree*/, TTree *clustersTree, AliESDEvent *pESD) const
{

  /***************************************************
  Resonstruct digits to vertex position
  ****************************************************/
  
  AliDebug(1,Form("Start FillESD T0"));

  TTree *treeR = clustersTree;
  
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
    Int_t i;
    Zposition = frecpoints -> GetVertex();
    timeStart = frecpoints -> GetMeanTime() ;
    for ( i=0; i<24; i++) {
      time[i] = Float_t (frecpoints -> GetTime(i)); // ps to ns
      amp[i] = frecpoints -> GetAmp(i);
    }
    pESD->SetT0zVertex(Zposition); //vertex Z position 
    pESD->SetT0(timeStart);        // interaction time 
    pESD->SetT0time(time);         // best TOF on each PMT 
    pESD->SetT0amplitude(amp);     // number of particles(MIPs) on each PMT
 
    AliDebug(1,Form(" Z position %f cm,  T0  %f ps",Zposition , timeStart));

} // vertex in 3 sigma






