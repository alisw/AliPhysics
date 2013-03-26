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
#include <AliAnalysisManager.h>
#include <AliProdInfo.h>

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
  AliLog::SetClassDebugLevel("AliT0TenderSupply",10); 
  Int_t run = fTender->GetRun();
  if (run == 0) return;    // to skip first init, when we don't have yet a run number

  // reset to no-action
  fCorrectMeanTime = kFALSE; //reset
  for(int i=0; i<4; i++) fTimeOffset[i]=0;
  fPass4LHC11aCorrection=kFALSE;
  fCorrectStartTimeOnAmplSatur = kFALSE;

  // check if MC
  AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  if (!inputHandler) return;
  TList *uiList = inputHandler->GetUserInfo();
  AliProdInfo prodInfo(uiList);
  prodInfo.List();
  Bool_t fIsMc = kFALSE;
  if ( mgr->GetMCtruthEventHandler() ) fIsMc=kTRUE;
  if (prodInfo.IsMC() == kTRUE) fIsMc=kTRUE;                // protection

  if (!fIsMc) {  // we consider actions only for data, not for MC

    // align T0s for LHC10d/LHC10e
    if ((run>=122195 &&  run<=126437) || (run>=127712 && run<=130850)) {
      AliInfo("Loading TZERO OCBD entries");
      fCorrectMeanTime=kTRUE;
      AliCDBManager* ocdbMan = AliCDBManager::Instance();
      ocdbMan->SetRun(fTender->GetRun());    
      AliCDBEntry *entry = ocdbMan->Get("T0/Calib/TimeAdjust/");
      if(entry) {
      AliT0CalibSeasonTimeShift *clb = (AliT0CalibSeasonTimeShift*) entry->GetObject();
      Float_t *t0means = clb->GetT0Means();
      for (Int_t i=0;i<4;i++) fTimeOffset[i] = t0means[i];
      } else {
	for (Int_t i=0;i<4;i++) fTimeOffset[i] = 0;
	AliWarning("T0Tender no T0 entry found T0shift set to 0");
      }
    }  
    
    // LHC11a pass4
    if (run>=139699 && run<=146860) {
      Int_t recoPass = prodInfo.GetRecoPass();
      if (recoPass < 0) {
	AliESDInputHandler *esdIH = dynamic_cast<AliESDInputHandler*>(fTender->GetESDhandler());
	if (esdIH) {
	  TTree *tree= (TTree*)esdIH->GetTree();
	  TFile *file= (TFile*)tree->GetCurrentFile();
	  if (file) {
	    TString fileName(file->GetName());
	    if (fileName.Contains("pass4") ) recoPass=4;
	  }
	}
      }
      if (recoPass == 4) fPass4LHC11aCorrection=kTRUE;
    }
	
    // LHC11h
    fAmplitudeThreshold = 100; //in mips
    if(167693<= run && run<=170593){  
      fCorrectStartTimeOnAmplSatur = kTRUE;
      fAmplitudeThreshold = 50; //in mips
    }

  }


  AliInfo("|******************************************************|");
  AliInfo(Form("|    Alice T0 Tender Initialisation (Run %d)       |",run));
  AliInfo("|    Settings:                                         |");
  AliInfo(Form("|    Monte Carlo flag               :  %d               |",fIsMc));
  AliInfo(Form("|    Adjust Offsets (LHC10d/LHC10e) :  %d               |",fCorrectMeanTime));
  AliInfo(Form("|    LHC11a pass4 patch   (LHC11a)  :  %d               |",fPass4LHC11aCorrection));
  AliInfo(Form("|    Amplitude Correction (LHC11h)  :  %d               |",fCorrectStartTimeOnAmplSatur));
  AliInfo("|******************************************************|");

}

//________________________________________________________________________
void AliT0TenderSupply::ProcessEvent(){

    //
    // loop over all online T0 candidates and flag
    // selected daughter tracks using the status bis of the TObject
    //

    AliESDEvent *event=fTender->GetEvent();
    if (!event) return;

    //Do something when the run number changed, like loading OCDB entries etc.
    if(fTender->RunChanged()) Init();

    
    if(fPass4LHC11aCorrection) {
      const Double32_t* mean = event->GetT0TOF();
      event->SetT0TOF(0, (mean[1]+mean[2])/2.);
     }
 
    if(fCorrectStartTimeOnAmplSatur){
        //correct A side ORA on amplitude saturation
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
            event->SetT0TOF(1, timeOrA);
        }
    }

    if(fCorrectMeanTime) {
      // correct mean time offsets  
      const Double32_t* mean = event->GetT0TOF();
      for(int it0=0; it0<3; it0++){
	if(-2000 < mean[it0]){
	  event->SetT0TOF(it0, mean[it0] - fTimeOffset[it0]); 
	}
      }
    }


}


