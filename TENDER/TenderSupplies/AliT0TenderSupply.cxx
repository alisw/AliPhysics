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
  fPass4LHC11aCorrection(kFALSE)
{
  //
  // default constructor
  //
  for(int i=0; i<4; i++) fTimeOffset[i]=0;
}

//________________________________________________________________________
AliT0TenderSupply::AliT0TenderSupply(const char *name, const AliTender *tender):
  AliTenderSupply(name,tender),
  fCorrectMeanTime(kFALSE),
  fCorrectStartTimeOnAmplSatur(kFALSE),
  fAmplitudeThreshold(100),
  fPass4LHC11aCorrection(kFALSE)
{
  //
  // constructor
  //
  for(int i=0; i<4; i++) fTimeOffset[i]=0;

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

 //   TString filename=Form("/tzero/alla/alice/ESDtree/2013/195483/alice_cern.ch_user_a_alla_treeLHC16c_output195483TOFcuts_000_%.3i_AnalysisResults.root",ifile);
 	
    if ( fileName.Contains("LHC16d3") || fileName.Contains("LHC16a2d2") ) {
   	 fCorrectMeanTime=kTRUE;
         fTimeOffset[0]=35; fTimeOffset[1]=25; fTimeOffset[2]=40; 
	 }
    Printf("fCorrectMeanTime %i \n", fCorrectMeanTime);
    // }  

}

//________________________________________________________________________
void AliT0TenderSupply::ProcessEvent(){

    //
    // loop over all online T0 candidates and flag
    // selected daughter tracks using the status bis of the TObject
    //

    AliESDEvent *event=fTender->GetEvent();
    if (!event) return;
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
      const Double32_t* mean = event->GetT0TOF();
      event->SetT0TOF(0, (mean[1]+mean[2])/2.);
 
    }
    //...........................................
    if(fCorrectStartTimeOnAmplSatur){
        //correct A side ORA on amplitude saturation
 	printf(" correct A side ORA on amplitude saturation\n");
       const Double32_t* time = event->GetT0time();
        const Double32_t* amplitude = event->GetT0amplitude();

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
            const Double32_t* mean = event->GetT0TOF();
            Double32_t timeOrC = mean[2];
            Double32_t timeOrAplusOrC = (timeOrA+timeOrC)/2;
            event->SetT0TOF(0, timeOrAplusOrC);
         }
    }

    //...........................................
    if(fCorrectMeanTime) {
      // correct mean time offsets  
      const Double32_t* mean = event->GetT0TOF();
      for(int it0=0; it0<3; it0++){
	if( mean[it0] < 2000  && mean[it0]>-2000 ) {
	  event->SetT0TOF(it0, mean[it0] - fTimeOffset[it0]);
	  printf("T0 %i mean[it0] %f fTimeOffset[it0] %f\n",it0, mean[it0],
	  fTimeOffset[it0]);
	  }
     }
    }
    //...........................................


}


