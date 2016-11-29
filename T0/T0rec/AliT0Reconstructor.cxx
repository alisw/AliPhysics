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
/*********************************************************************
 *  T0 reconstruction and filling ESD
 *  - reconstruct mean time (interation time) 
 *  - vertex position
 *  -  multiplicity
 ********************************************************************/

#include <AliESDEvent.h>
#include "AliLog.h"
#include "AliT0RecPoint.h"
#include "AliRawReader.h"
#include "AliT0RawReader.h"
#include "AliT0digit.h"
#include "AliT0Reconstructor.h"
#include "AliT0Parameters.h"
#include "AliT0Calibrator.h"
#include "AliESDfriend.h"
#include "AliESDTZERO.h"
//#include "AliESDTZEROfriend.h"
#include "AliLog.h"
#include "AliCDBEntry.h" 
#include "AliCDBManager.h"
#include "AliCTPTimeParams.h"
#include "AliLHCClockPhase.h"
#include "AliT0CalibSeasonTimeShift.h"
#include "AliESDRun.h"
#include "AliGRPObject.h"

#include <TArrayI.h>
#include <TGraph.h>
#include <TMath.h>
#include <Riostream.h>
#include <TBits.h>

ClassImp(AliT0Reconstructor)

  AliT0Reconstructor:: AliT0Reconstructor(): AliReconstructor(),
					     fdZonA(0),
					     fdZonC(0),
					     fZposition(0),
					     fParam(NULL),
					     fAmpLEDrec(),
					     fQTC(0),
					     fAmpLED(0),
                                             fCalib(),
                                             fLatencyHPTDC(9000),
                                             fLatencyL1(0),
                                             fLatencyL1A(0),
                                             fLatencyL1C(0),
					     fGRPdelays(0),
					     fTimeMeanShift(0x0),
					     fTimeSigmaShift(0x0),
                                             fESDTZERO(NULL),
                                             fESDTZEROfriend(NULL),
                                             fIsCDFfromGRP(kFALSE), 
                                             fMeanOrA(0),
                                             fMeanOrC(0),
                                             fMeanTVDC(0),
                                             fLHCperiod(kFALSE)
{
  for (Int_t i=0; i<24; i++)  { fTime0vertex[i] =0; fQT1mean[i]=0;}

  //constructor
  AliCDBEntry *entry = AliCDBManager::Instance()->Get("GRP/CTP/CTPtiming");
  if (!entry) AliFatal("CTP timing parameters are not found in OCDB !");
  AliCTPTimeParams *ctpParams = (AliCTPTimeParams*)entry->GetObject();
  Float_t l1Delay = (Float_t)ctpParams->GetDelayL1L0()*25.0;

  AliCDBEntry *entry1 = AliCDBManager::Instance()->Get("GRP/CTP/TimeAlign");
  if (!entry1) AliFatal("CTP time-alignment is not found in OCDB !");
  AliCTPTimeParams *ctpTimeAlign = (AliCTPTimeParams*)entry1->GetObject();
  l1Delay += ((Float_t)ctpTimeAlign->GetDelayL1L0()*25.0);
 
  AliCDBEntry *entry4 = AliCDBManager::Instance()->Get("GRP/Calib/LHCClockPhase");
  if (!entry4) AliFatal("LHC clock-phase shift is not found in OCDB !");
  AliLHCClockPhase *phase = (AliLHCClockPhase*)entry4->GetObject();
  
  fGRPdelays = l1Delay - phase->GetMeanPhase();
  
  AliCDBEntry *entry5 = AliCDBManager::Instance()->Get("T0/Calib/TimeAdjust");
  if (entry5) {
    AliT0CalibSeasonTimeShift *timeshift = (AliT0CalibSeasonTimeShift*)entry5->GetObject();
    fTimeMeanShift = timeshift->GetT0Means();
    fTimeSigmaShift  = timeshift->GetT0Sigmas();
  }
  else
    AliWarning("Time Adjust is not found in OCDB !");
  
  fParam = AliT0Parameters::Instance();
  fParam->Init();
  
  for (Int_t i=0; i<24; i++){
    TGraph* gr = fParam ->GetAmpLEDRec(i);
    if (gr) fAmpLEDrec.AddAtAndExpand(gr,i) ; 
    TGraph* gr1 = fParam ->GetAmpLED(i);
    if (gr1) fAmpLED.AddAtAndExpand(gr1,i) ; 
    TGraph* gr2 = fParam ->GetQTC(i);
    if (gr2) fQTC.AddAtAndExpand(gr2,i) ; 	
    fTime0vertex[i] = fParam->GetCFD(i);
    fQT1mean[i] = fParam->GetQT1(i);
    fPedestal[i] = fParam->GetPedestalOld(i);
    printf(" OCDB fTime0vertex %f fQT1mean %f pedestal %f \n",fTime0vertex[i], fQT1mean[i],fPedestal[i] );
  }
  fMeanOrA = fParam->GetMeanOrA();
  fMeanOrC = fParam->GetMeanOrC();
  fMeanTVDC = fParam->GetMeanVertex();


  fLatencyL1 = fParam->GetLatencyL1();
  fLatencyL1A = fParam->GetLatencyL1A(); 
  fLatencyL1C = fParam->GetLatencyL1C();
  fLatencyHPTDC = fParam->GetLatencyHPTDC();
  AliDebug(2,Form(" LatencyL1 %f latencyL1A %f latencyL1C %f latencyHPTDC %f \n",fLatencyL1, fLatencyL1A, fLatencyL1C, fLatencyHPTDC));
 
  for (Int_t i=0; i<24; i++) {
     if( fTime0vertex[i] < 500 || fTime0vertex[i] > 60000)
       { fTime0vertex[i] =( 1000.*fLatencyHPTDC - 1000.*fLatencyL1 + 1000.*fGRPdelays)/24.4;
	 fIsCDFfromGRP=kTRUE;
       }
     AliDebug(2,Form("OCDB mean CFD time %i %f \n",i, fTime0vertex[i]));
  }
  //here real Z position
  fdZonC = TMath::Abs(fParam->GetZPosition("T0/C/PMT1"));
  fdZonA = TMath::Abs(fParam->GetZPosition("T0/A/PMT15"));
  printf("!!!!fdZonC %f fdZonA %f \n",fdZonC, fdZonA);
  
  fCalib = new AliT0Calibrator();
  fESDTZERO  = new AliESDTZERO();
  //LHC period
   AliCDBEntry* entry6 = AliCDBManager::Instance()->Get("GRP/GRP/Data");
  AliGRPObject* grpData = dynamic_cast<AliGRPObject*>(entry6->GetObject());
  if (!grpData) {printf("Failed to get GRP data for run"); return;}
  TString LHCperiod = grpData->GetLHCPeriod();
  if(LHCperiod.Contains("LHC15")|| LHCperiod.Contains("LHC16")) fLHCperiod=kTRUE;
  printf(" LHCperiod %i\n",fLHCperiod);
}

//_____________________________________________________________________________
void AliT0Reconstructor::Reconstruct(TTree*digitsTree, TTree*clustersTree) const
{
  // T0 digits reconstruction
  Int_t refAmp = 0 ; /*Int_t (GetRecoParam()->GetRefAmp());*/
  
  TArrayI * timeCFD = new TArrayI(24); 
  TArrayI * timeLED = new TArrayI(24); 
  TArrayI * chargeQT0 = new TArrayI(24); 
  TArrayI * chargeQT1 = new TArrayI(24); 

 
  Float_t c = 29.9792458; // cm/ns
  Float_t channelWidth = fParam->GetChannelWidth() ;  
  Double32_t vertex = 9999999, meanVertex = 0 ;
  Double32_t timeDiff=999999, meanTime=999999, timeclock=999999;
  
  
  AliDebug(1,Form("Start DIGITS reconstruction "));
  
  Float_t lowAmpThreshold =  GetRecoParam()->GetAmpLowThreshold();  
  Float_t highAmpThreshold =  GetRecoParam()->GetAmpHighThreshold(); 
  printf("Reconstruct(TTree*digitsTree highAmpThreshold %f  lowAmpThreshold %f \n",lowAmpThreshold, highAmpThreshold);

  //shift T0A, T0C , T0AC
  Float_t shiftA = GetRecoParam() -> GetLow(310);  
  Float_t shiftC = GetRecoParam() -> GetLow(311);  
  Float_t shiftAC = GetRecoParam() -> GetLow(312);
  AliDebug(2, Form("Reconstruct(TTree*digitsTree shiftA %f shiftC %f shiftAC %f \n",shiftA, shiftC, shiftAC));

  Double32_t besttimeA=9999999;  Double32_t besttimeA_best=0;
  Double32_t besttimeC=9999999;  Double32_t besttimeC_best=0;
  Int_t timeDelayCFD[24]; 
  Int_t badpmt[24];
  //Bad channel
  for (Int_t i=0; i<24; i++) {
    badpmt[i] = GetRecoParam() -> GetBadChannels(i);
    timeDelayCFD[i] =  Int_t (fParam->GetTimeDelayCFD(i));
  }
  fCalib->SetEq(0);
  TBranch *brDigits=digitsTree->GetBranch("T0");
  AliT0digit *fDigits = new AliT0digit() ;
  if (brDigits) {
    brDigits->SetAddress(&fDigits);
  }else{
    AliError(Form("EXEC Branch T0 digits not found"));
    return;
  }
  
  digitsTree->GetEvent(0);
  digitsTree->GetEntry(0);
  brDigits->GetEntry(0);
  fDigits->GetTimeCFD(*timeCFD);
  fDigits->GetTimeLED(*timeLED);
  fDigits->GetQT0(*chargeQT0);
  fDigits->GetQT1(*chargeQT1);
  
  Int_t onlineMean =  fDigits->MeanTime();
  Int_t corridor = GetRecoParam() -> GetCorridor();  
 
  Bool_t tr[5];
  for (Int_t i=0; i<5; i++) tr[i]=false; 
  
  
  AliT0RecPoint frecpoints;
  AliT0RecPoint * pfrecpoints = &frecpoints;
  clustersTree->Branch( "T0", "AliT0RecPoint" ,&pfrecpoints);
  Int_t timecenterA = 511;
  Int_t timecenterC = 511;
  //   if(fLHCperiod) { timecenterC=512;timecenterA=516;}
  
  Float_t time[24], adc[24], adcmip[24];
  for (Int_t ipmt=0; ipmt<24; ipmt++) {
  if(timeCFD->At(ipmt)>0) printf(" pmt %i time %i low %i up %i\n",
  ipmt, timeCFD->At(ipmt)- timeDelayCFD[ipmt],511-corridor, 511+corridor);
      if( (timeCFD->At(ipmt) - timeDelayCFD[ipmt])>511-corridor &&
	  (timeCFD->At(ipmt) - timeDelayCFD[ipmt])<511+corridor  ) {
	 
      Float_t timefull = 0.001*( timeCFD->At(ipmt) - 511 - timeDelayCFD[ipmt])  * channelWidth;
      frecpoints.SetTimeFull(ipmt, 0 ,timefull) ;
      if(( chargeQT1->At(ipmt) - chargeQT0->At(ipmt))>0)  
	adc[ipmt] = chargeQT1->At(ipmt) - chargeQT0->At(ipmt);
      else
	adc[ipmt] = 0; 
      // no walk correction for 2015 data
      time[ipmt] = timeCFD->At(ipmt) -  timeDelayCFD[ipmt];
      if(fLHCperiod ) {
	if (ipmt<12) time[ipmt] =   time[ipmt] - timecenterC;
	if (ipmt>=12) time[ipmt] =   time[ipmt] -timecenterA ;
      }
      else
	{
	  time[ipmt] = fCalib-> WalkCorrection(refAmp, ipmt, Int_t(adc[ipmt]),  timeCFD->At(ipmt)) ;
	  time[ipmt] =   time[ipmt] - timecenterC;
	}
      Double_t sl = Double_t(timeLED->At (ipmt) - timeCFD->At(ipmt));
      //    time[ipmt] = fCalib-> WalkCorrection( refAmp,ipmt, Int_t(sl),  timeCFD->At(ipmt) ) ;
      //   AliDebug(5,Form(" ipmt %i QTC  %i , time in chann %i (led-cfd) %i ",
      //	    ipmt, Int_t(adc[ipmt]) ,Int_t(time[ipmt]),Int_t( sl)));
      //  printf(" ipmt %i QTC  %i , time in chann %i \n ",
      //		    ipmt, Int_t(adc[ipmt]) ,Int_t(time[ipmt]));
      
    Double_t ampMip = 0;
      TGraph* ampGraph = (TGraph*)fAmpLED.At(ipmt);
      if (ampGraph) ampMip = ampGraph->Eval(sl);
      Double_t qtMip = 0;
      TGraph* qtGraph = (TGraph*)fQTC.At(ipmt);
      if (qtGraph) qtMip = qtGraph->Eval(adc[ipmt]);
      AliDebug(5,Form("  Amlitude in MIPS LED %f ,  QTC %f in channels %f\n ",ampMip,qtMip, adc[ipmt]));
       frecpoints.SetTime(ipmt, Float_t(time[ipmt]) );
      frecpoints.SetAmpLED(ipmt, Float_t( ampMip)); 
      frecpoints.SetAmp(ipmt, Float_t(qtMip));
      adcmip[ipmt]=qtMip;
      
    }
    else {
      time[ipmt] = -99999;
      adc[ipmt] = 0;
      adcmip[ipmt] = 0;
      
    }
  }
  Int_t npmtsC=0;
  for (Int_t ipmt=0; ipmt<12; ipmt++){
    if(time[ipmt] !=0  && time[ipmt] != -99999
       &&  adcmip[ipmt]>lowAmpThreshold && adcmip[ipmt]<highAmpThreshold )
      {
	if(time[ipmt]<besttimeC) besttimeC=time[ipmt]; //timeC
	//	if(TMath::Abs(time[ipmt])<TMath::Abs(besttimeC_best)) 
	//	  besttimeC_best=time[ipmt]; //timeC
	besttimeC_best += time[ipmt];  //sum of timeC 
	npmtsC++;
      }
  }
  Int_t npmtsA=0;
  for ( Int_t ipmt=12; ipmt<24; ipmt++)
    {
      if(time[ipmt] != 0 && time[ipmt] != -99999
	 && adcmip[ipmt]>lowAmpThreshold && adcmip[ipmt]<highAmpThreshold)
	{
	  if(time[ipmt]<besttimeA) besttimeA=time[ipmt]; 
	  //	  if(TMath::Abs(time[ipmt] ) < TMath::Abs(besttimeA_best)) 
	  //    besttimeA_best=time[ipmt]; //timeA
	  besttimeA_best += time[ipmt];  //sum of timeA 
	  npmtsA++;
	}
    }
  if (npmtsC>0) besttimeC_best = besttimeC_best/npmtsC;
  if (npmtsA>0) besttimeA_best = besttimeA_best/npmtsA;
  
  if( besttimeA < 999999 && besttimeA!=0) {
    frecpoints.SetTimeBestA((besttimeA_best * channelWidth  - fdZonA/c)  );
    frecpoints.SetTime1stA((besttimeA * channelWidth  - fdZonA/c - shiftA) );
    //    if(fLHCperiod ) 
    //   frecpoints.SetTime1stA((besttimeA * channelWidth  - fdZonA/c - fTimeMeanShift[1]) );
    tr[1]=true;
  }
 printf(" 1stimeA %f besttimeA %f fdZonCA%f  shiftA %f \n",
	besttimeA * channelWidth,besttimeA_best * channelWidth, fdZonC/c, fTimeMeanShift[1]);
  
  if( besttimeC < 999999 && besttimeC!=0) {
    frecpoints.SetTimeBestC((besttimeC_best * channelWidth  - fdZonC/c) );
    frecpoints.SetTime1stC((besttimeC * channelWidth  - fdZonC/c - shiftC) );
    //   if(fLHCperiod ) 
    //    frecpoints.SetTime1stC((besttimeC * channelWidth  - fdZonC/c - fTimeMeanShift[2]) );
   tr[2]=true;
  }
  //  printf(" 1stimeC %f besttimeC %f fdZonC %f  shiftC %f \n",
  //	 besttimeC * channelWidth,besttimeC_best * channelWidth, fdZonC/c, fTimeMeanShift[2]);

  AliDebug(5,Form("1stimeA %f , besttimeA %f 1sttimeC %f besttimeC %f ",
		  besttimeA, besttimeA_best,
		  besttimeC, besttimeC_best) );

  if(besttimeA <999999 && besttimeC < 999999 ){
    //    timeDiff = (besttimeC - besttimeA)*channelWidth;
    timeDiff = (besttimeA - besttimeC)*channelWidth;
    meanTime = channelWidth * (besttimeA_best + besttimeC_best)/2. ; 
    timeclock = channelWidth * (besttimeA + besttimeC)/2. - shiftAC ;
    //  if(fLHCperiod) 
    //  timeclock = channelWidth * (besttimeA + besttimeC)/2. - fTimeMeanShift[0] ;
   vertex = meanVertex - 0.001* c*(timeDiff)/2.;// + (fdZonA - fdZonC)/2;
    tr[0]=true; 
  }
  frecpoints.SetVertex(vertex);
  frecpoints.SetMeanTime(meanTime );
  frecpoints.SetT0clock(timeclock );
  frecpoints.SetT0Trig(tr);
  
 AliDebug(5,Form("fRecPoints:::  1stimeA %f , besttimeA %f 1sttimeC %f besttimeC %f vertex %f",
		  frecpoints.Get1stTimeA(),  frecpoints.GetBestTimeA(),
		  frecpoints.Get1stTimeC(),  frecpoints.GetBestTimeC(), 
		  vertex ) );
  
  AliDebug(5,Form("T0 triggers %d %d %d %d %d",tr[0],tr[1],tr[2],tr[3],tr[4]));
  
  //online mean
  frecpoints.SetOnlineMean(Int_t(onlineMean));
  AliDebug(10,Form("  timeDiff %f #channel,  meanTime %f #channel, vertex %f cm online mean %i timeclock %f ps",timeDiff, meanTime,vertex, Int_t(onlineMean), timeclock));
  
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
  //
  Float_t meanOrA=0, meanOrC=0, meanTVDC=0, meanQT1[24]={0};
  if (fMeanOrA==0)  meanOrA = fTime0vertex[0] + 587;
  else 
    meanOrA=fMeanOrA;
  if (fMeanOrC==0) meanOrC = fTime0vertex[0] + 678;
  else 
    meanOrC=fMeanOrC;
  if(  fMeanTVDC==0)  meanTVDC = fTime0vertex[0] + 2564;
  else 
    meanTVDC=fMeanTVDC;
  for (int i=0; i<24; i++) {
    if (fQT1mean[i]==0)  meanQT1[i]= fTime0vertex[0] + 2564;
    else 
      meanQT1[i]=fQT1mean[i];
  } 
 // old or new QTC: 0 old ; >0 new
  Int_t oldORnew = GetRecoParam()->GetLow(205);
  Int_t timeDelayCFD[24]; 
  Int_t corridor = GetRecoParam() -> GetCorridor();  
  if(fIsCDFfromGRP) corridor *=5;
  Int_t badpmt[24];
  //Bad channel
  for (Int_t i=0; i<24; i++) {
    badpmt[i] = GetRecoParam() -> GetBadChannels(i);
    timeDelayCFD[i] =  Int_t (fParam->GetTimeDelayCFD(i));
  }
  Int_t equalize = GetRecoParam() -> GetEq();
  fCalib->SetEq(equalize);
  Int_t low[500], high[500];
  Float_t timefull=-99999;;
  Float_t tvdc  = -99999; Float_t ora = -99999; Float_t orc = -99999;
  
  Int_t timeCFD[24], timeLED[24], chargeQT0[24], chargeQT1[24];
  Float_t time2zero[24]; 
  Double32_t timeDiff, meanTime, timeclock;
  timeDiff =  meanTime = timeclock = 9999999;
  Float_t c = 29.9792458; // cm/ns
  Double32_t vertex = 9999999;
  Int_t amplitude[26], amplitudeNew[26];
  Int_t onlineMean=0;
  Float_t meanVertex = 0;
   for (Int_t i0=0; i0<24; i0++) {
    low[i0] = Int_t(fTime0vertex[i0]) - corridor;
    high[i0] = Int_t(fTime0vertex[i0]) + corridor;
    time2zero[i0] = 99999;
   }
  
  Int_t alldata[250][5];   // container for readed raw 
  for (Int_t i0=0; i0<250; i0++)
    for (Int_t j0=0; j0<5; j0++)  alldata[i0][j0]=0; 
  
  Float_t lowAmpThreshold =  GetRecoParam()->GetAmpLowThreshold();  
  Float_t highAmpThreshold =  GetRecoParam()->GetAmpHighThreshold();
  
  Double32_t besttimeA=9999999;  Double32_t besttimeA_best=0;
  Double32_t besttimeC=9999999;  Double32_t besttimeC_best=0;

  Float_t channelWidth = fParam->GetChannelWidth() ;  
  
  AliT0RecPoint frecpoints;
  AliT0RecPoint * pfrecpoints = &frecpoints;
  
  recTree->Branch( "T0", "AliT0RecPoint" ,&pfrecpoints);
  
  AliDebug(10," before read data ");
  AliT0RawReader myrawreader(rawReader);
  
  UInt_t type =rawReader->GetType();
  
  if (!myrawreader.Next())
    AliDebug(1,Form(" no raw data found!!"));
  else
    {  
      for (Int_t i=0; i<24; i++)
	{
	  timeCFD[i]=0; timeLED[i]=0;
	}
      
      if(type == 7  ) {  //only physics 
	for (Int_t i=0; i<226; i++) {
	  for (Int_t iHit=0; iHit<5; iHit++) {
	    alldata[i][iHit] = myrawreader.GetData(i,iHit);
	  }
	}
	
	Int_t fBCID=Int_t (rawReader->GetBCID());
	Int_t trmbunch= myrawreader.GetTRMBunchID();
	AliDebug(10,Form(" CDH BC ID %i, TRM BC ID %i \n", fBCID, trmbunch ));
	if( (trmbunch-fBCID)!=37  ) {
	  AliDebug(0,Form("wrong :::: CDH BC ID %i, TRM BC ID %i \n", fBCID, trmbunch ));
	  //	  type = -1;
	}
	for (Int_t in=0; in<12; in++)  
	  {
	    for (Int_t iHit=0; iHit<5; iHit++) {
		if(alldata[in+1][iHit] > low[in] && 
		   alldata[in+1][iHit] < high[in])
		  {
		    //		    printf(" ::Reconstruct :: readed i %i hit %i cfd %i \n",
		    //		       in+1,iHit, alldata[in+1][iHit] ); 
		    timeCFD[in] = alldata[in+1][iHit] ; 
 		    break;
		  }
		
	      }
	    for (Int_t iHit=0; iHit<5; iHit++) {
		if(alldata[in+1+56][iHit] > low[in+12] && 
		   alldata[in+1+56][iHit] < high[in+12])
		  {
		    //		    printf(" ::Reconstruct :: readed i %i hit %i cfd %i \n",
		    //		   in+12,iHit, alldata[in+1+56][iHit] ); 
		    timeCFD[in+12] = alldata[in+56+1][iHit] ;
		    break;
		  }
	    }	    
	  }
	ReadNewQTC(alldata, amplitudeNew);
	ReadOldQTC(alldata, amplitude);
	onlineMean = alldata[49][0];
	
	Double32_t time[24],  adcmip[24], ampnewmip[24];
	Int_t adc[24];
	for (Int_t ipmt=0; ipmt<24; ipmt++) {
	  if(timeCFD[ipmt] >  0 && amplitude[ipmt]>0 ){
	   //for simulated data
	     //for physics  data
	    //	   adc[ipmt] = fAmplitude[ipmt];
	   Int_t refAmp = Int_t (fTime0vertex[ipmt]);
	   adc[ipmt]=amplitude[ipmt];
	   time[ipmt] = fCalib-> WalkCorrection( refAmp, ipmt, adc[ipmt], timeCFD[ipmt] ) ;
	   Double32_t qtMip = 0;
	   TGraph * qtGraph = (TGraph*)fQTC.At(ipmt);
	   if (qtGraph) {
	     // if(oldORnew)
	     // qtMip = qtGraph->Eval(Float_t (adc[ipmt]) ) ;
	     // else 
	     qtMip = qtGraph->Eval(Float_t (adc[ipmt]) - fPedestal[ipmt] );
	   }
	   if( equalize  ==0 ) 
	     frecpoints.SetTime(ipmt, Float_t(time[ipmt]) );
	   else 
	     frecpoints.SetTime(ipmt, Float_t(time[ipmt] + fTime0vertex[ipmt]) );
	   // frecpoints.SetTime(ipmt, Float_t(time[ipmt] ) );
	   if(qtMip<0) qtMip=0;
	   frecpoints.SetAmp(ipmt,  qtMip); 
	   adcmip[ipmt]=qtMip;
	   ampnewmip[ipmt]=Double32_t(amplitudeNew[ipmt]);
	   frecpoints.SetAmpLED(ipmt,ampnewmip[ipmt] ); //new amplitude 
	  }
	  else {
	    time[ipmt] = -9999;
	    adc[ipmt] = 0;
	    adcmip[ipmt] = 0;
	  }
	}
	Int_t npmtsC=0;
	for (Int_t ipmt=0; ipmt<12; ipmt++){
	  if(time[ipmt] !=0 &&  time[ipmt] > -9000 
	     /*&& badpmt[ipmt]==0 */
	    &&  adcmip[ipmt]>lowAmpThreshold )
	    {
	      if(time[ipmt]<besttimeC) besttimeC=time[ipmt]; //timeC
	      //	     if(TMath::Abs(time[ipmt])<TMath::Abs(besttimeC_best)) 
	      //	       besttimeC_best=time[ipmt]; //timeC	     
	      besttimeC_best += time[ipmt];  //sum of timeA 
	      npmtsC++;
	    }
	}
	Int_t npmtsA=0;
	for ( Int_t ipmt=12; ipmt<24; ipmt++)
	  {
	    if(time[ipmt] != 0 &&  time[ipmt] > -9000 
	       /* && badpmt[ipmt]==0*/ 
	       && adcmip[ipmt]>lowAmpThreshold )
	      {
		if(time[ipmt]<besttimeA) besttimeA=time[ipmt]; 
		// if(TMath::Abs(time[ipmt] ) < TMath::Abs(besttimeA_best)) 
		//	 besttimeA_best=time[ipmt]; //timeA
		besttimeA_best += time[ipmt];  //sum of timeA 
		npmtsA++;
	      }
	 }
	if (npmtsC>0) besttimeC_best = besttimeC_best/npmtsC;
	if (npmtsA>0) besttimeA_best = besttimeA_best/npmtsA;
	
	if(besttimeA < 999999 && besttimeA!=0 ) {
	  if( equalize  ==0 ) 
	    frecpoints.SetTime1stA((besttimeA * channelWidth)- 1000.*fLatencyHPTDC + 1000.*fLatencyL1A - 1000.*fGRPdelays - fTimeMeanShift[1] ); 
	  else
	    {
	      frecpoints.SetTimeBestA((besttimeA_best * channelWidth )); 
	      frecpoints.SetTime1stA((besttimeA * channelWidth - fTimeMeanShift[1])); 
	    }
       }
	if( besttimeC < 999999 && besttimeC!=0) {
	  if( equalize  ==0 ) 
	    frecpoints.SetTime1stC((besttimeC * channelWidth)- 1000.*fLatencyHPTDC +1000.*fLatencyL1C - 1000.*fGRPdelays - fTimeMeanShift[2]);
	  else
	    {
	      frecpoints.SetTimeBestC((besttimeC_best * channelWidth ));
	      frecpoints.SetTime1stC((besttimeC * channelWidth - fTimeMeanShift[2]));
	    }
	}
	AliDebug(5,Form("1stimeA %f , besttimeA %f 1sttimeC %f besttimeC %f ",
			besttimeA, besttimeA_best,
		       besttimeC, besttimeC_best) );
	AliDebug(5,Form("fRecPoints:::  1stimeA %f , besttimeA %f 1sttimeC %f besttimeC %f shiftA %f shiftC %f ",
			frecpoints.Get1stTimeA(),  frecpoints.GetBestTimeA(),
			frecpoints.Get1stTimeC(),  frecpoints.GetBestTimeC(), 
			fTimeMeanShift[1], fTimeMeanShift[2] ) );
	if( besttimeC < 999999 &&  besttimeA < 999999) { 
	  if(equalize  ==0 )
	    timeclock = (channelWidth*(besttimeC + besttimeA)/2.- 1000.*fLatencyHPTDC +1000.*fLatencyL1 - 1000.*fGRPdelays - fTimeMeanShift[0]);
	 else
	   {
	     timeclock = channelWidth * Float_t( besttimeA+besttimeC)/2. - fTimeMeanShift[0];
	     meanTime = channelWidth * Float_t(besttimeA_best + besttimeC_best )/2.;
	   }
	 timeDiff = ( besttimeA - besttimeC)* 0.001* channelWidth ;
	 vertex =  meanVertex - c*(timeDiff)/2. ; //+ (fdZonA - fdZonC)/2; 
       }
       
      }  //if phys event       
       AliDebug(1,Form("  timeDiff %f #channel,  meanTime %f #ps, TOFmean%f  vertex %f cm meanVertex %f  \n",timeDiff, meanTime,timeclock, vertex,meanVertex));
       frecpoints.SetT0clock(timeclock);
       frecpoints.SetVertex(vertex);
       frecpoints.SetMeanTime(meanTime);
       frecpoints.SetOnlineMean(Int_t(onlineMean));
      
      // Set triggers
      Bool_t tr[5];
       Int_t trchan[5] = {50,51,52,55,56};
       for (Int_t i=0; i<5; i++) tr[i] = false; 
       Int_t triggername[3];//20ns around 0
       triggername[0] = (TMath::Abs(meanTVDC)<2147483647)?(Int_t)meanTVDC:0;
       triggername[1] = (TMath::Abs(meanOrA) <2147483647)?(Int_t)meanOrA:0;
       triggername[2] = (TMath::Abs(meanOrC) <2147483647)?(Int_t)meanOrC:0;
       
       for (Int_t itr=0; itr<5; itr++) {
         for (Int_t iHit=0; iHit<5; iHit++) 
           {
             Int_t trr=trchan[itr];
             if(itr<3 ) { 
              if( (alldata[trr][iHit] - triggername[itr]) > -800 &&
                   (alldata[trr][iHit] - triggername[itr]) < 800)  tr[itr]=true;
            }
             else 
            if( alldata[trr][iHit] > 0)  tr[itr]=true;
             AliDebug(5,Form("Reconstruct :::  T0 triggers iHit %i tvdc %d orA %d orC %d centr %d semicentral %d",iHit, tr[0],tr[1],tr[2],tr[3],tr[4]));
           }       
      }
       frecpoints.SetT0Trig(tr);  
      
      // all times without amplitude correction 
      Float_t timecent;
      for (Int_t iHit=0; iHit<5; iHit++) 
	{
	  timefull = timecent = -9999; 
	  tvdc = ora = orc = -9999;
	  if(alldata[50][iHit]>0) 
	    tvdc = (Float_t(alldata[50][iHit]) - meanTVDC) * channelWidth* 0.001; 
	  if(alldata[51][iHit]>0)
	    ora = (Float_t(alldata[51][iHit]) - meanOrA) * channelWidth* 0.001;
	  if(alldata[52][iHit]>0) 
	    orc = (Float_t(alldata[52][iHit]) - meanOrC) * channelWidth* 0.001;
	  
	  frecpoints.SetOrC( iHit, orc);
	  frecpoints.SetOrA( iHit, ora);
	  frecpoints.SetTVDC( iHit, tvdc);
	  for (Int_t i0=0; i0<24; i0++) {
	    if (equalize  ==0 )  
	      timecent = fTime0vertex[i0] + timeDelayCFD[i0];
	    else
	      timecent = fTime0vertex[i0];
	    timefull = -9999; 
	    if (i0<12) 
            {
	      if(alldata[i0+1][iHit]>1) {
		timefull = (Float_t(alldata[i0+1][iHit]) - timecent)* channelWidth* 0.001;
		//	if(alldata[i0+1][iHit]>1) printf ("Reconstruct ::: RAW pmt %i hit %i time %f readed %i\n",i0,iHit, timefull,alldata[i0+1][iHit] );
	      }
	      else 
		if(alldata[i0+45][iHit]>1) {
		  timefull = (Float_t(alldata[i0+45][iHit]) - timecent)* channelWidth* 0.001;
		  //	  if(alldata[i0+45][iHit]>1) printf ("Reconstruct ::: RAW pmt %i hit %i time %f readed %i\n",i0,iHit, timefull,alldata[i0+45][iHit] );
		}
            }
	    frecpoints.SetTimeFull(i0, iHit,timefull) ;
	  }
	}
      //set MPD
      frecpoints.SetMultA(Float_t(amplitudeNew[24]) );
      frecpoints.SetMultC(Float_t(amplitudeNew[25]) );
      //      printf ("Reconstruct :::  T0 MPDA %i MPDC %i", amplitude[24], amplitude[25]);

      //FIT CFD
      Double32_t timeFIT=0;
      for (int i=0; i<4; i++) {
       	for (Int_t iHit=0; iHit<5; iHit++) {	
	  if(alldata[i+211][iHit]>9000 && alldata[i+211][iHit]<12000) {
	    timeFIT=  Double32_t(alldata[i+211][iHit]);
	    frecpoints.SetFITTime(i, timeFIT);
	    break;
	  }
	}
      }
      
    } // if (else )raw data
  recTree->Fill();
}
  
  
//____________________________________________________________
  
  void AliT0Reconstructor::FillESD(TTree */*digitsTree*/, TTree *clustersTree, AliESDEvent *pESD) const
  {

  /***************************************************
  Resonstruct digits to vertex position
  ****************************************************/
  
  AliDebug(1,Form("Start FillESD T0"));
  if(!pESD) {
    AliError("No ESD Event");
    return;
  }
  pESD ->SetT0spread(fTimeSigmaShift);

  Float_t channelWidth = fParam->GetChannelWidth() ;  
  Float_t c = 0.0299792458; // cm/ps
  Float_t currentVertex=0, shift=0;
  Int_t ncont=-1;
  const AliESDVertex* vertex = pESD->GetPrimaryVertex();
  if (!vertex)        vertex = pESD->GetPrimaryVertexSPD();
  if (!vertex)        vertex = pESD->GetPrimaryVertexTPC();
  if (!vertex)        vertex = pESD->GetVertex();

  if (vertex) {
    AliDebug(2, Form("Got %s (%s) from ESD: %f", 
		    vertex->GetName(), vertex->GetTitle(), vertex->GetZ()));
    currentVertex = vertex->GetZ();
    
    ncont = vertex->GetNContributors();
    if(ncont>0 ) {
      shift = currentVertex/c;
    }
  }
  TTree *treeR = clustersTree;
  
  AliT0RecPoint frecpoints;
  AliT0RecPoint * pfrecpoints = &frecpoints;
  
  AliDebug(1,Form("Start FillESD T0"));
  TBranch *brRec = treeR->GetBranch("T0");
  if (brRec) {
    brRec->SetAddress(&pfrecpoints);
  }else{
    AliError(Form("EXEC Branch T0 rec not found"));
    return;
  } 
  
  brRec->GetEntry(0);
  Double32_t ampnew[24], time[24], ampQTC[24], timecorr[24];  
  Double32_t* tcorr;
  for(Int_t i=0; i<24; i++) 
    ampnew[i]=time[i]=ampQTC[i]=timecorr[i]=0;

  //1st time 
  Double32_t timeClock[3];
  Double32_t zPosition = frecpoints.GetVertex();

  timeClock[0] = frecpoints.GetT0clock() ;
  timeClock[1] = frecpoints.Get1stTimeA() + shift;
  timeClock[2] = frecpoints.Get1stTimeC() - shift;
  //best time
  Double32_t timemean[3];
  timemean[0] = frecpoints.GetMeanTime();
  timemean[1] = frecpoints.GetBestTimeA() + shift;
  timemean[2] = frecpoints.GetBestTimeC() - shift;

  for(Int_t i=0; i<3; i++) {
    fESDTZERO->SetT0TOF(i,timeClock[i]);   // interaction time (ns) 
    fESDTZERO->SetT0TOFbest(i,timemean[i]);   // interaction time (ns) 
  }
  for ( Int_t i=0; i<24; i++) {
    time[i] =  frecpoints.GetTime(i); // ps to ns
    if ( time[i] != 0 && time[i]>-9999) {
      ampQTC[i] = frecpoints.GetAmp(i);
      ampnew[i] = frecpoints.AmpLED(i);
      AliDebug(1,Form("T0: %i  time %f  ampold %f ampnew %f \n", i, time[i], ampQTC[i], ampnew[i]));
      }
  }
  //   for ( Int_t i=0; i<24; i++)   
  // printf("T0: %i  time %f  ampQTC %f ampNewQTC %f \n", i, time[i], ampQTC[i], ampnew[i]);
  fESDTZERO->SetT0time(time);         // best TOF on each PMT 
  fESDTZERO->SetT0amplitude(ampQTC);     // amplitude old QTC
  fESDTZERO->SetT0NewAmplitude(ampnew);     // amplitude new QTC
  
  Int_t trig= frecpoints.GetT0Trig();
  frecpoints.PrintTriggerSignals( trig);
  //  printf(" !!!!! FillESD trigger %i \n",trig);
  fESDTZERO->SetT0Trig(trig);
  fESDTZERO->SetT0zVertex(zPosition); //vertex Z position 

  Double32_t multA=frecpoints.GetMultA();
  Double32_t multC=frecpoints.GetMultC();
  fESDTZERO->SetMultA(multA); // for backward compatubility
  fESDTZERO->SetMultC(multC); // for backward compatubility


  for (Int_t iHit =0; iHit<5; iHit++ ) {
       AliDebug(10,Form("FillESD ::: iHit %i tvdc %f orA %f orC %f\n", iHit,
	   frecpoints.GetTVDC(iHit),
	   frecpoints.GetOrA(iHit),
		       frecpoints.GetOrC(iHit) ));
    fESDTZERO->SetTVDC(iHit,frecpoints.GetTVDC(iHit));
    fESDTZERO->SetOrA(iHit,frecpoints.GetOrA(iHit));
    fESDTZERO->SetOrC(iHit,frecpoints.GetOrC(iHit));
    
    for (Int_t i0=0; i0<24; i0++) 
	fESDTZERO->SetTimeFull(i0, iHit,frecpoints.GetTimeFull(i0,iHit));	
    //FIT CFD
    for (Int_t i0=0; i0<4; i0++) 
     fESDTZERO->SetPileupTime(i0, frecpoints.GetFITTime(i0)); //// 19.05.2016
   
  AliDebug(1,Form("T0: SPDshift %f Vertex %f (T0A+T0C)/2 best %f #ps T0signal %f ps OrA %f ps OrC %f ps T0trig %i\n",shift, zPosition, timemean[0], timeClock[0], timeClock[1], timeClock[2], trig));

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // background flags
  Bool_t background = BackgroundFlag();
  fESDTZERO->SetBackgroundFlag(background);
  Bool_t pileup =  PileupFlag();
  fESDTZERO->SetPileupFlag(pileup);
  TBits pileupbits = SetPileupBits();
  fESDTZERO->SetPileupBits(pileupbits);
  TBits pileout =fESDTZERO-> GetT0PileupBits();
  pileout.Print();

  //  for (Int_t i=0; i<5; i++) {
  // fESDTZERO->SetPileupTime(i, frecpoints.GetTVDC(i) ) ;
    //   printf("!!!!!! FillESD :: pileup %i %f %f \n", i,fESDTZERO->GetPileupTime(i), frecpoints.GetTVDC(i));
  }
  Bool_t sat  = SatelliteFlag();
  fESDTZERO->SetSatelliteFlag(sat);
  
  
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (pESD) 
   pESD->SetTZEROData(fESDTZERO);
 
  

} // vertex in 3 sigma

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 //____________________________________________________________
  
Bool_t AliT0Reconstructor::PileupFlag() const
{
  //
  Bool_t pileup = false;
  Float_t tvdc[5];
  for (Int_t ih=0; ih<5; ih++) 
    {
      tvdc[ih] =  fESDTZERO->GetTVDC(ih);
      
     if(  tvdc[0] !=0 &&  tvdc[0]> -10 && tvdc[0]< 10 )
       if(ih>0 && tvdc[ih]>20 ) 	pileup = true; 
     if( tvdc[0] >20 || (tvdc[0] < -20  &&  tvdc[0] > -9000) ) pileup =true;
     //    if (pileup) printf(" !!!!! pile up %i  tvdc %f \n",ih,  tvdc[ih]); 
    }


  return pileup;

}

 //____________________________________________________________
  
TBits AliT0Reconstructor::SetPileupBits() const
{
  TBits pileup ;
  Float_t tvdc[5];
  Int_t pos, bc[21];
  UInt_t ibc;
  pileup.ResetAllBits();
  for ( Int_t nbc=0; nbc<21; nbc++) bc[nbc]=0;
  for (Int_t ih=0; ih<5; ih++) 
    {
      tvdc[ih] =  fESDTZERO->GetTVDC(ih);
      if(tvdc[ih]!=0 && tvdc[ih]>-290 &&tvdc[ih]<290 ) {
	if( tvdc[ih]>0) pos = Int_t (tvdc[ih]+6)/25;
	if(tvdc[ih]<0&&tvdc[ih]>-290)  pos = Int_t (tvdc[ih]-6)/25;	
	//	printf("AliT0Reconstructor::PileupFlag():: hit %i tvdc %f pos %i bc %i\n",ih,tvdc[ih],pos, bc[pos+10]);

	bc[pos+10] = 1;
      }
    }
  for ( Int_t nbc=0; nbc<21; nbc++) {
    if(bc[10]>0) {
      ibc=UInt_t(nbc);
      if (bc[nbc]>0)  pileup.SetBitNumber(ibc,kTRUE);
    }
  }
  
  //  pileup.Print();
  return pileup;
  }
 //____________________________________________________________
  
Bool_t AliT0Reconstructor::BackgroundFlag() const
{
 
  Bool_t background = false;
  /*  
      Float_t orA = fESDTZERO->GetOrA(0);
      Float_t orC = fESDTZERO->GetOrC(0);
      Float_t tvdc =  fESDTZERO->GetTVDC(ih);
      
      if ( (orA > -5 && orA <5) && (orC > -5 && orC <5) && (tvdc < -5 || tvdc > 5)) {
      background = true;
      //   printf(" orA %f orC %f tvdc %f\n", orA, orC, tvdc);
      } 
  */ 
  return background;


}


 //____________________________________________________________
  
Bool_t  AliT0Reconstructor::SatelliteFlag() const
{

 Float_t satelliteLow = GetRecoParam() -> GetLowSatelliteThreshold();
  Float_t satelliteHigh = GetRecoParam() -> GetHighSatelliteThreshold();
  Bool_t satellite = false;
  for (Int_t i0=0; i0<24; i0++) {
    Float_t timefull =  fESDTZERO -> GetTimeFull(i0,0);
    if( timefull > satelliteLow && timefull < satelliteHigh)  satellite=true;
  }
	
  return satellite;

}
 //____________________________________________________________
void  AliT0Reconstructor::ReadNewQTC(Int_t alldata[250][5], Int_t amplitude[26]) const
{
  // QT00 -> QT11
  printf("@@ readNewQTC");
  Float_t a[26], b[26];
  Int_t qt01mean[26], qt11mean[26];
  for(int i=0; i<26; i++) {
     a[i] = GetRecoParam() -> GetLow(i+130);
     b[i] = GetRecoParam() -> GetLow(i+156);
     if(i<24) 
       qt11mean[i] =qt01mean[i] =fTime0vertex[i] + 15500;
     else
      qt11mean[i] =qt01mean[i] =fTime0vertex[0] + 15500;
 
     amplitude[i]=0;
  }
  Int_t diff[4];
  Int_t pmt;
   //new QTC C side
  for (Int_t ik=0; ik<106; ik+=4)
    {
      for(int id=0; id<2; id++) diff[id] = 0;
      if (ik<48)          pmt=ik/4;
      if (ik>47 && ik<52) pmt= 24;   
      if (ik>51 && ik<56) pmt= 25; 
      if(ik>55)           pmt=(ik-8)/4;

     for(Int_t iHt = 0; iHt<5; iHt++) {
	if(alldata[107+ik+1][iHt] > (qt01mean[pmt]-1000) &&
	   alldata[107+ik+1][iHt] < (qt01mean[pmt]+1000) ) {
	  diff[0]=alldata[107+ik][iHt] - alldata[107+ik+1][iHt];
	  //	  printf(" newQTC 00 ik %i iHt %i pmt %i  QT00 %i QT01 %i \n", ik, iHt, pmt, alldata[107+ik][iHt],  alldata[107+ik+1][iHt]);
	  break;
	}
      }
      for(Int_t iHt = 0; iHt<5; iHt++) {
	if( alldata[107+ik+3][iHt] > (qt11mean[pmt]-1000) &&
	    alldata[107+ik+3][iHt] < (qt11mean[pmt]+1000) ) {
	  diff[1]=alldata[107+ik+2][iHt] - alldata[107+ik+3][iHt];
	  //	  printf(" newQTC 11 ik %i iHt %i pmt %i QT10 %i QT11 %i \n", ik, iHt, pmt, alldata[107+ik+2][iHt],  alldata[107+ik+3][iHt]);
	  break;
	}
      }
      if(diff[0] != 0)  amplitude[pmt]=diff[0];
      if(diff[1] != 0)  {
	amplitude[pmt] = a[pmt]*diff[1] + b[pmt];  
	//	if (pmt==24 || pmt==25) printf(" @@@ new MPD pmt %i amp %f a %f b %f \n",pmt,  amplitude[pmt],a[pmt], b[pmt]);
      }
      //    if(diff[0] == 0 &&diff[1]==0) amplitude[pmt]=0;
    }
}
 //____________________________________________________________
void  AliT0Reconstructor::ReadOldQTC(Int_t alldata[250][5], Int_t amplitude[26] ) const
{
  Int_t  chargeQT0[26], chargeQT1[26], pedestal[26];
  Float_t meanQT1[26];
  for (int i=0; i<24; i++) {
    if (fQT1mean[i]==0)  meanQT1[i]= fTime0vertex[0] + 2564;
    else 
      meanQT1[i]=fQT1mean[i];
  } 
  for (int i=24; i<26; i++) meanQT1[i]= fQT1mean[23];

  for (Int_t i0=0; i0<26; i0++) {
    amplitude[i0]=0;
    chargeQT1[i0]=0;
    chargeQT0[i0]=0;
  }
  Int_t ind[26];
  for (int iii=0; iii<12; iii++) ind[iii]=25;
  for (int iii=12; iii<24; iii++) ind[iii]=57;
  ind[24]=5;
  ind[25]=55;
  for (Int_t in=0; in<24;  in++)
    {
      /*      if(in==24|| in==25)
	printf(" MPD %i %i %i data %i %i \n",
	in, 2*in+ind[in],  ind[in], alldata[2*in+ind[in]+1][0], alldata[2*in+ind[in]][0]);*/
      for (Int_t iHit=0; iHit<5; iHit++) 
	{
	  if (alldata[2*in+ind[in]+1][iHit] > meanQT1[in]-800 &&  
	      alldata[2*in+ind[in]+1][iHit] < meanQT1[in]+800 ) {
	    chargeQT1[in] = alldata[2*in+ind[in]+1][iHit];
	    break;
	  }
	}
      for (Int_t iHit=0; iHit<5; iHit++) 
	{
	  if( (alldata[2*in+ind[in]][iHit] - chargeQT1[in])>fPedestal[in] &&
	      chargeQT1[in]>0)
	    {
	      chargeQT0[in]=alldata[2*in+ind[in]][iHit];
	      //      printf(" readedOld QTC Raw %i %i %i\n",  in, chargeQT0[in],chargeQT1[in]);
 	      break;
	    }
	}
      
      if( (chargeQT0[in]-chargeQT1[in])>fPedestal[in]) {
	amplitude[in]=chargeQT0[in]-chargeQT1[in];
	//	printf(" OLD amplitude PMT %i %i \n",in, amplitude[in]);
      }
    }
}    


