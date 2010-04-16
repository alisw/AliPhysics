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

#include <TArrayI.h>
#include <TGraph.h>
#include <TMath.h>
#include <Riostream.h>

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
                                             fLatencyL1C(0)

{
  //constructor

  fParam = AliT0Parameters::Instance();
  fParam->Init();
 
  for (Int_t i=0; i<24; i++){
        TGraph* gr = fParam ->GetAmpLEDRec(i);
	if (gr) fAmpLEDrec.AddAtAndExpand(gr,i) ; 
	  TGraph* gr1 = fParam ->GetAmpLED(i);
	  if (gr1) fAmpLED.AddAtAndExpand(gr1,i) ; 
	  TGraph* gr2 = fParam ->GetQTC(i);
	  if (gr2) fQTC.AddAtAndExpand(gr2,i) ; 	
  }

  fLatencyL1 = fParam->GetLatencyL1();
  fLatencyL1A = fParam->GetLatencyL1A();
  fLatencyL1C = fParam->GetLatencyL1C();
  fLatencyHPTDC = fParam->GetLatencyHPTDC();
  AliDebug(10,Form(" LatencyL1 %f latencyL1A %f latencyL1C %f latencyHPTDC %f \n",fLatencyL1, fLatencyL1A, fLatencyL1C, fLatencyHPTDC));
  
  // fdZonC = TMath::Abs(fParam->GetZPositionShift("T0/C/PMT1"));
  //fdZonA = TMath::Abs(fParam->GetZPositionShift("T0/A/PMT15"));


  fCalib = new AliT0Calibrator();

}

//_____________________________________________________________________________
void AliT0Reconstructor::Reconstruct(TTree*digitsTree, TTree*clustersTree) const
  
{
  // T0 digits reconstruction
  Int_t refAmp = GetRecoParam()->GetRefAmp();
  
  TArrayI * timeCFD = new TArrayI(24); 
  TArrayI * timeLED = new TArrayI(24); 
  TArrayI * chargeQT0 = new TArrayI(24); 
  TArrayI * chargeQT1 = new TArrayI(24); 

 
  Float_t channelWidth = fParam->GetChannelWidth() ;  
  Float_t meanVertex = fParam->GetMeanVertex();
  Float_t c = 0.0299792; // cm/ps
  Double32_t vertex = 9999999;
  Double32_t timeDiff=999999, meanTime=999999, timeclock=999999;

  
  AliDebug(1,Form("Start DIGITS reconstruction "));
  

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

  Bool_t tr[5];
  for (Int_t i=0; i<5; i++) tr[i]=false; 
  
  Double32_t besttimeA=999999;
  Double32_t besttimeC=999999;
  Int_t pmtBestA=99999;
  Int_t pmtBestC=99999;
  
  AliT0RecPoint* frecpoints= new AliT0RecPoint ();
  clustersTree->Branch( "T0", "AliT0RecPoint" ,&frecpoints, 405,1);
  
  Float_t time[24], adc[24];
  for (Int_t ipmt=0; ipmt<24; ipmt++) {
    if(timeCFD->At(ipmt)>0 ){
     if(( chargeQT1->At(ipmt) - chargeQT0->At(ipmt))>0)  
	adc[ipmt] = chargeQT1->At(ipmt) - chargeQT0->At(ipmt);
      else
	adc[ipmt] = 0;
      
     // time[ipmt] = fCalib-> WalkCorrection(refAmp, ipmt, adc[ipmt],  timeCFD->At(ipmt)) ;
	     
      Double_t sl = Double_t(timeLED->At(ipmt) - timeCFD->At(ipmt));
      time[ipmt] = fCalib-> WalkCorrection( refAmp,ipmt, Int_t(sl),  timeCFD->At(ipmt) ) ;
      AliDebug(10,Form(" ipmt %i QTC %i , time in chann %i (led-cfd) %i ",
		       ipmt, Int_t(adc[ipmt]) ,Int_t(time[ipmt]),Int_t( sl)));

      Double_t ampMip =((TGraph*)fAmpLED.At(ipmt))->Eval(sl);
      Double_t qtMip = ((TGraph*)fQTC.At(ipmt))->Eval(adc[ipmt]);
      AliDebug(10,Form("  Amlitude in MIPS LED %f ,  QTC %f in channels %i\n ",ampMip,qtMip, adc[ipmt]));
      
      frecpoints->SetTime(ipmt, Float_t(time[ipmt]) );
      frecpoints->SetAmp(ipmt, Float_t( ampMip)); //for cosmic &pp beam 
      frecpoints->SetAmpLED(ipmt, Float_t(qtMip));
      
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
  if(besttimeA < 999999) {
    frecpoints->SetTimeBestA(Int_t(besttimeA *channelWidth));
    tr[1]=true;
  }
  if( besttimeC < 999999 ) {
    frecpoints->SetTimeBestC(Int_t(besttimeC *channelWidth));
    tr[2]=true;
  }
  AliDebug(10,Form(" besttimeA %f ch,  besttimeC %f ch",besttimeA, besttimeC));
  if(besttimeA <999999 && besttimeC < 999999 ){
    //    timeDiff = (besttimeC - besttimeA)*channelWidth;
    timeDiff = (besttimeA - besttimeC)*channelWidth;
    meanTime = (besttimeA + besttimeC)/2;// * channelWidth); 
    timeclock = meanTime *channelWidth ;
    vertex = meanVertex - c*(timeDiff)/2.;// + (fdZonA - fdZonC)/2;
    tr[0]=true; 
  }
  frecpoints->SetVertex(vertex);
  frecpoints->SetMeanTime(meanTime);
  frecpoints->SetT0clock(timeclock);
  frecpoints->SetT0Trig(tr);

  for (Int_t i=0; i<5; i++) {
    printf(" T0 trigers %i ",tr[i]);
  }
    printf(" \n ");

  //online mean
  frecpoints->SetOnlineMean(Int_t(onlineMean));
  AliDebug(10,Form("  timeDiff %i #channel,  meanTime %i #channel, vertex %f cm online mean %i timeclock %i ps",timeDiff, meanTime,vertex, Int_t(onlineMean), timeclock));
  
  

   
  
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
  // reference amplitude and time ref. point from reco param

  Int_t refAmp = GetRecoParam()->GetRefAmp();
  Int_t refPoint = GetRecoParam()->GetRefPoint();

  Int_t allData[110][5];
  
  Int_t timeCFD[24], timeLED[24], chargeQT0[24], chargeQT1[24];
  Double32_t timeDiff=999999, meanTime=999999, timeclock=999999;
  Float_t c = 29.9792458; // cm/ns
  Double32_t vertex = 9999999;
  Int_t onlineMean=0;
  // Float_t meanVertex = fParam->GetMeanVertex();
  Float_t meanVertex = 0;

  AliDebug(1,Form(" @@@@ Latency ",fLatencyL1));
  for (Int_t i0=0; i0<105; i0++)
    {
      for (Int_t j0=0; j0<5; j0++) allData[i0][j0]=0; 	
    }
   
  Double32_t besttimeA=9999999;
  Double32_t besttimeC=9999999;
  Int_t pmtBestA=99999;
  Int_t pmtBestC=99999;
   
  AliT0RecPoint* frecpoints= new AliT0RecPoint ();
  
  recTree->Branch( "T0", "AliT0RecPoint" ,&frecpoints, 405,1);
   
  AliDebug(10," before read data ");
  AliT0RawReader myrawreader(rawReader);

  UInt_t type =rawReader->GetType();

  if (!myrawreader.Next())
    AliDebug(1,Form(" no raw data found!!"));
  else
    {  
      if(type == 7) {  //only physics 
      for (Int_t i=0; i<105; i++) {
	for (Int_t iHit=0; iHit<5; iHit++) 
	  {
	    allData[i][iHit] = myrawreader.GetData(i,iHit);
	  }
      }
      Int_t ref=0;
      if (refPoint>0) 
      ref = allData[refPoint][0]-5000;

      Float_t channelWidth = fParam->GetChannelWidth() ;  
      
      //       Int_t meanT0 = fParam->GetMeanT0();
       
	  
	  for (Int_t in=0; in<12; in++)  
	    {
	      timeCFD[in] = allData[in+1][0] ;
	      timeCFD[in+12] = allData[in+56+1][0] ;
	      timeLED[in] = allData[in+12+1][0] ;
	      timeLED[in+12] = allData[in+68+1][0] ;
	      AliDebug(10, Form(" readed i %i cfdC %i cfdA %i ledC %i ledA%i ",
				in, timeCFD[in],timeCFD[in+12],timeLED[in], 
				timeLED[in+12]));   
	    }
	  
	  for (Int_t in=0; in<12;  in++)
	    {
	      chargeQT0[in]=allData[2*in+25][0];
	      chargeQT1[in]=allData[2*in+26][0];
	    }
	  
	   for (Int_t in=12; in<24;  in++)
	     {
	       chargeQT0[in]=allData[2*in+57][0];
	       chargeQT1[in]=allData[2*in+58][0];
	     }
	   
	   //	 } //cosmic with physics event
       for (Int_t in=0; in<24; in++)  
	 AliDebug(10, Form(" readed Raw %i %i %i %i %i",
			   in, timeLED[in],timeCFD[in],chargeQT0[in],chargeQT1[in]));
        onlineMean = allData[49][0];       
       
       Double32_t time[24], adc[24];
       for (Int_t ipmt=0; ipmt<24; ipmt++) {
	 if(timeCFD[ipmt]>0 && timeLED[ipmt]>0){
	   //for simulated data
	     //for physics  data
	   if(( chargeQT1[ipmt] - chargeQT0[ipmt])>0)  
	     adc[ipmt] = chargeQT0[ipmt] - chargeQT1[ipmt];
	   else
	     adc[ipmt] = 0;
	   

	   //	   time[ipmt] = fCalib-> WalkCorrection(refAmp, ipmt, adc[ipmt], timeCFD[ipmt] ) ;
	   
      	   Double_t sl = timeLED[ipmt] - timeCFD[ipmt];
	   time[ipmt] = fCalib-> WalkCorrection( refAmp,ipmt, Int_t(sl), timeCFD[ipmt] ) ;
	   AliDebug(10,Form(" ipmt %i QTC %i , time in chann %i (led-cfd) %i ",
			    ipmt, Int_t(adc[ipmt]) ,Int_t(time[ipmt]),Int_t( sl)));
	   Double_t ampMip =( (TGraph*)fAmpLED.At(ipmt))->Eval(sl);
	   Double_t qtMip = ((TGraph*)fQTC.At(ipmt))->Eval(adc[ipmt]);
	   AliDebug(10,Form("  Amlitude in MIPS LED %f ; QTC %f;  in channels %i\n ",ampMip,qtMip, adc[ipmt]));
	     
	   frecpoints->SetTime(ipmt, Float_t(time[ipmt]) );
	   // frecpoints->SetTime(ipmt,Double32_t(timeCFD[ipmt]));
	   frecpoints->SetAmpLED(ipmt, Double32_t( qtMip)); //for cosmic &pp beam 
	   frecpoints->SetAmp(ipmt, Double32_t(ampMip));
	     
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
       if(besttimeA < 999999) 
	 frecpoints->SetTimeBestA(0.001* besttimeA * channelWidth - fLatencyHPTDC + fLatencyL1A);
       if( besttimeC < 999999 ) 
	 frecpoints->SetTimeBestC( 0.001 *besttimeC * channelWidth - fLatencyHPTDC + fLatencyL1C);
       AliDebug(10,Form(" pmtA %i besttimeA %f ps, pmtC %i besttimeC %f ps",
		       pmtBestA,besttimeA, pmtBestC,  besttimeC));
        if(besttimeA <999999 && besttimeC < 999999 ){
	 timeDiff = ( besttimeA - besttimeC) *0.001 * channelWidth + fLatencyL1A - fLatencyL1C;
	 timeclock = 0.001*channelWidth * Float_t( besttimeA+besttimeC)/2. - fLatencyHPTDC + fLatencyL1;  
	 meanTime = (besttimeA+besttimeC-2.*Float_t(ref))/2.;
	 vertex =  meanVertex - c*(timeDiff)/2. ; //+ (fdZonA - fdZonC)/2; 
	}
      }  //if phys event       
      AliDebug(5,Form("  timeDiff %f #channel,  meanTime %f #channel, TOFmean%f  vertex %f cm meanVertex %f online mean %i \n",timeDiff, meanTime,timeclock, vertex,meanVertex, onlineMean));
      frecpoints->SetT0clock(timeclock);
      frecpoints->SetVertex(vertex);
      frecpoints->SetMeanTime(meanTime);
      frecpoints->SetOnlineMean(Int_t(onlineMean));
	// Set triggers
      
      Bool_t tr[5];
      Int_t trchan[5]= {50,51,52,55,56};
      for (Int_t i=0; i<5; i++) tr[i]=false; 
      for (Int_t itr=0; itr<5; itr++) {
	if(allData[trchan[itr]][0]>0) tr[itr]=true;
	frecpoints->SetT0Trig(tr);
      }
      
    } // if (else )raw data
  recTree->Fill();
  if(frecpoints) delete frecpoints;
}
  
  
  //____________________________________________________________
  
  void AliT0Reconstructor::FillESD(TTree */*digitsTree*/, TTree *clustersTree, AliESDEvent *pESD) const
  {

  /***************************************************
  Resonstruct digits to vertex position
  ****************************************************/
  
  AliDebug(1,Form("Start FillESD T0"));
  Float_t channelWidth = fParam->GetChannelWidth() ;  
  Float_t c = 29.9792458; // cm/ns
  Float_t currentVertex=0, shift=0;
  Int_t ncont=0;
  const AliESDVertex* vertex = pESD->GetPrimaryVertex();
  if (!vertex)        vertex = pESD->GetPrimaryVertexSPD();
  if (!vertex)        vertex = pESD->GetPrimaryVertexTPC();
  if (!vertex)        vertex = pESD->GetVertex();

  if (vertex) {
    AliDebug(2, Form("Got %s (%s) from ESD: %f", 
		    vertex->GetName(), vertex->GetTitle(), vertex->GetZ()));
    currentVertex = vertex->GetZ();
    
    ncont = vertex->GetNContributors();
    // cout<<" spdver "<<spdver<<" ncont "<<ncont<<endl;
    if(ncont>2 ) {
      shift = currentVertex/c;
      //	  cout<<" vertex shif "<<shift<<" vertex "<<spdver<<" IsFromVertexer3D  "<<fverSPD->IsFromVertexer3D()<<endl;
    }
  }
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
    AliError(Form("EXEC Branch T0 rec not found"));
    return;
  } 
    
    brRec->GetEntry(0);
    Double32_t amp[24], time[24], ampQTC[24], timecorr[24];  
    Double32_t timeClock[3];
    Double32_t zPosition = frecpoints -> GetVertex();
    Double32_t timeStart = frecpoints -> GetMeanTime();
    timeClock[0] = frecpoints -> GetT0clock() ;
    timeClock[1] = frecpoints -> GetBestTimeA() + shift;
    timeClock[2] = frecpoints -> GetBestTimeC() - shift;
    for ( Int_t i=0; i<24; i++) {
      time[i] =  frecpoints -> GetTime(i); // ps to ns
      amp[i] = frecpoints -> GetAmp(i);
      ampQTC[i] = frecpoints -> AmpLED(i);
    }
    Int_t trig= frecpoints ->GetT0Trig();
    pESD->SetT0Trig(trig);
    
    pESD->SetT0zVertex(zPosition); //vertex Z position 
    pESD->SetT0(timeStart);        // interaction time 
    for(Int_t i=0; i<3; i++) 
      pESD->SetT0TOF(i,timeClock[i]);   // interaction time (ns) 
    pESD->SetT0time(time);         // best TOF on each PMT 
    pESD->SetT0amplitude(amp);     // number of particles(MIPs) on each PMT
    
    AliDebug(1,Form("T0: Vertex %f (T0A+T0C)/2 %f #channels T0signal %f ns OrA %f ns OrC %f T0trig %i\n",zPosition, timeStart, timeClock[0], timeClock[1], timeClock[2], trig));

    /* if (pESD) {
     AliESDfriend *fr = (AliESDfriend*)pESD->FindListObject("AliESDfriend");
     if (fr) {
        AliDebug(1, Form("Writing TZERO friend data to ESD tree"));

	if (ncont>2) {
	  for ( Int_t i=0; i<24; i++) {
	    if(i<12) timecorr[i]=time[i]-shift*channelWidth;
	    if(i>11) timecorr[i]=time[i]-shift*channelWidth;
	    fr->SetT0timeCorr(timecorr) ;
	  }
            fr->SetT0ampLEDminCFD(amp);
            fr->SetT0ampQTC(ampQTC);
	}
    }
  }
    */
    
    
} // vertex in 3 sigma






