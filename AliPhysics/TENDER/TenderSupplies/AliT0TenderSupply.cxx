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

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//  T0 Tender supply    //
//  //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <AliESDEvent.h>
#include <AliESDtrack.h>
#include <AliTender.h>
#include <AliT0TenderSupply.h>
#include <AliCDBManager.h>
#include <AliCDBEntry.h>
#include <AliT0CalibSeasonTimeShift.h>
#include <AliESDInputHandler.h>

ClassImp(AliT0TenderSupply)

//________________________________________________________________________
AliT0TenderSupply::AliT0TenderSupply():
  AliTenderSupply(),
  fCorrectMeanTime(kFALSE),
  fCorrectStartTimeOnAmplSatur(kFALSE),
  fAmplitudeThreshold(100), 
  fPass4LHC11aCorrection(kFALSE),
  fLHC16j(kFALSE),
  fEvent(0),
  fnEvent(0)
{
  //
  // default constructor
  //
  for(int i=0; i<4; i++) fTimeOffset[i]=0;
  for(int i=0; i<24; i++) fFixMeanCFD[i]=0;
  
}

//________________________________________________________________________
AliT0TenderSupply::AliT0TenderSupply(const char *name, const AliTender *tender):
  AliTenderSupply(name,tender),
  fCorrectMeanTime(kFALSE),
  fCorrectStartTimeOnAmplSatur(kFALSE),
  fAmplitudeThreshold(100),
  fPass4LHC11aCorrection(kFALSE),
  fLHC16j(kFALSE),
  fEvent(0)  ,
  fnEvent(0)
{
  //
  // constructor
  //
  for(int i=0; i<4; i++) fTimeOffset[i]=0;
  for(int i=0; i<24; i++) fFixMeanCFD[i]=0;

}

//________________________________________________________________________
AliT0TenderSupply::~AliT0TenderSupply(){
  //
  // destructor
  //
  
}

//________________________________________________________________________
void AliT0TenderSupply::Init(){
  // Init
  //
  fCorrectStartTimeOnAmplSatur = kFALSE;
  fAmplitudeThreshold = 100; //in mips
   Int_t run = fTender->GetRun();
  if (run == 0) return;    // to skip first init, when we don't have yet a run number
  Printf("----------- TZERO Tender ----------------");

  fCorrectMeanTime = kFALSE; //reset
  for(int i=0; i<4; i++) fTimeOffset[i]=0;

  // align T0s for LHC10def periods 
 // LHC11h

  //  if (fTender->GetRun()>=167693 &&  fTender->GetRun()<=170593){
    //saturation
    //    fCorrectStartTimeOnAmplSatur = kTRUE;
    // fAmplitudeThreshold = 100; //in mips
    Printf("Loading TZERO OCBD entries for run %i\n", fTender->GetRun() );
    AliESDInputHandler *esdIH = dynamic_cast<AliESDInputHandler*>  (fTender->GetESDhandler());
    TTree *tree= (TTree*)esdIH->GetTree();
    TFile *file= (TFile*)tree->GetCurrentFile();
    TString fileName(file->GetName());
    if (fTender->GetRun()>=167693 &&  fTender->GetRun()<=170593 && 
	fileName.Contains("pass2") ) {
	fCorrectMeanTime=kTRUE;
	AliCDBManager* ocdbMan = AliCDBManager::Instance();
        ocdbMan->SetRun(fTender->GetRun());    
        AliCDBEntry *entry = ocdbMan->Get("T0/Calib/TimeAdjust/");
 //   AliCDBEntry *entry = ocdbMan->Get("T0/Calib/TimeOffsetAOD");
        if(entry) {
            AliT0CalibSeasonTimeShift *clb = (AliT0CalibSeasonTimeShift*) entry->GetObject();
            Float_t *t0means = clb->GetT0Means();
            for (Int_t i=0;i<4;i++) fTimeOffset[i] = t0means[i];
       } else {
            for (Int_t i=0;i<4;i++) fTimeOffset[i] = 0;
             AliWarning("T0Tender no T0 entry found T0shift set to 0");
              }
    }

    if ( fileName.Contains("LHC16d3") || fileName.Contains("LHC16a2d2") ) {
      fCorrectMeanTime=kTRUE;
      fTimeOffset[0]=35; fTimeOffset[1]=25; fTimeOffset[2]=40; 
	 }
    
    fnEvent=0;
    if (fileName.Contains("LHC16j2a") ) {
      fLHC16j=kTRUE;
      Float_t corr[24] = {1.5, 3.4, 8.4, -2.5, -5, 1.8, 0.6, 2.0, -19, -11, -7, 7,
			  -4, 2., -5, 0, 2.5, 3.2, 3.2, 2.2, -2, 3.8, 3.1, 6};
      for(int i=0; i<24; i++) 	 fFixMeanCFD[i]=corr[i];
      fTimeOffset[0]=-57; fTimeOffset[1]=-34; fTimeOffset[2]=-72;
    }
    
    if (fileName.Contains("LHC16j2b") ) {
      fLHC16j=kTRUE;
      Float_t corr[24] = {2.6, 5.6, 8.2, -0.3, -6., 3.7,
			  0., 5., -17., -11, -7, 8.5,
			  -0.8, 8., -0., 2, 6., 7.,
			  6., 4.4, 2., 6.38, 6., 8.98};
      for(int i=0; i<24; i++) fFixMeanCFD[i]=corr[i];     
      fTimeOffset[0]=-76; fTimeOffset[1]=-80; fTimeOffset[2]=-67;
    }
    
    if (fileName.Contains("LHC16j4") ) {
      fLHC16j=kTRUE;
      Float_t corr[24] = {1.5, 1.9, 7, -1.8, -6, 2, -1, 1.3, -17, -8.6, -7.8,5.5,
			  -4, 3, -4.6, 0, 3.6, 3.4, 3, 2.4, 0, 4.2, 3, 6};
      for(int i=0; i<24; i++) fFixMeanCFD[i]=corr[i];     
      fTimeOffset[0]=-36; fTimeOffset[1]=-31; fTimeOffset[2]=-42;
    }
}

//________________________________________________________________________
void AliT0TenderSupply::ProcessEvent(){

    //
    // loop over all online T0 candidates and flag
    // selected daughter tracks using the status bis of the TObject
    //

    fEvent=fTender->GetEvent();
    if (!fEvent) return;
     //...........................................
   //Do something when the run number changed, like loading OCDB entries etc.
     if(fTender->RunChanged()) Init();
    
   if(fTender->RunChanged()){
      if (fTender->GetRun()>=139699&&  fTender->GetRun()<=146860){
     Init();
         AliESDInputHandler *esdIH = dynamic_cast<AliESDInputHandler*>  (fTender->GetESDhandler());
        if (esdIH) {
          TTree *tree= (TTree*)esdIH->GetTree();
          TFile *file= (TFile*)tree->GetCurrentFile();
         if (file){
            TString fileName(file->GetName());
	    if (fileName.Contains("pass4") ) fPass4LHC11aCorrection=kTRUE;
	  }
	}
      }
    }
    
    if(fPass4LHC11aCorrection) {
      const Double32_t* mean = fEvent->GetT0TOF();
      fEvent->SetT0TOF(0, (mean[1]+mean[2])/2.);
 
    }
    //...........................................
    if(fCorrectStartTimeOnAmplSatur){
        //correct A side ORA on amplitude saturation
 	printf(" correct A side ORA on amplitude saturation\n");
       const Double32_t* time = fEvent->GetT0time();
        const Double32_t* amplitude = fEvent->GetT0amplitude();

        Int_t idxOfFirstPmtA = -1;
        Double32_t timeOrA   = 99999;
        for(int ipmt=12; ipmt<24; ipmt++){ //loop over A side
            if( amplitude[ipmt] < fAmplitudeThreshold){
                if( time[ipmt] > -200 && time[ipmt]!=0 && time[ipmt] < timeOrA ){ 
                    timeOrA        = time[ipmt];
                    idxOfFirstPmtA = ipmt;
                }
            }
        }

        if(idxOfFirstPmtA>-1){ //a hit in aside with less than 40 mips
            const Double32_t* mean = fEvent->GetT0TOF();
            Double32_t timeOrC = mean[2];
            Double32_t timeOrAplusOrC = (timeOrA+timeOrC)/2;
            fEvent->SetT0TOF(0, timeOrAplusOrC);
         }
    }

    //...........................................
    if(fCorrectMeanTime) {
      // correct mean time offsets  
      const Double32_t* mean = fEvent->GetT0TOF();
      for(int it0=0; it0<3; it0++){
	if( mean[it0] > -2000  && mean[it0]<2000 ) {
	  fEvent->SetT0TOF(it0, mean[it0] - fTimeOffset[it0]);
	  }
      }
    }
    //...........................................
    if( fLHC16j) {
      const Double32_t* time=fEvent->GetT0time();
      Float_t c = 0.0299792458; // cm/ps
      Float_t currentVertex=fEvent->GetPrimaryVertex()->GetZ();
      Float_t shift=currentVertex/c;
      SetLHC16jMC(time, shift);
      fnEvent++;
     }
}

//-------------------------------------------------------
void  AliT0TenderSupply::SetLHC16jMC(const Double32_t* time, Float_t shift)
{

  Float_t besttimeA=999999; Float_t besttimeC=999999;
  Float_t t0A, t0C, timecorr[24];
  for (int i=0; i<12; i++) {
    if ( time[i]!=0) {
      timecorr[i] = time[i]- fFixMeanCFD[i];
      if(timecorr[i]<besttimeC)  besttimeC=timecorr[i];
    }
  }
  for (int i=12; i<24; i++) {
    if( time[i]!=0) {
       timecorr[i] = time[i]- fFixMeanCFD[i];
      if(timecorr[i]<besttimeA)  besttimeA=timecorr[i];
    }
  }
  if(besttimeC<1000 && besttimeC>-1000) {
    t0C=24.4*besttimeC - shift - fTimeOffset[2];
    fEvent->SetT0TOF(2, t0C);
  }	    
  if(besttimeA<1000 && besttimeA>-1000) {
    t0A=24.4*besttimeA + shift - fTimeOffset[1];
    fEvent->SetT0TOF(1, t0A);
  }	    
  if(besttimeC<1000 && besttimeC>-1000
     && besttimeA<1000 && besttimeA>-1000) 
    fEvent->SetT0TOF(0, 24.4*(besttimeC+besttimeA)/2. - fTimeOffset[0]);

}


