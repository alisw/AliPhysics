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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TOF tender: - load updated calibrations (for TOF and T0)
//             - set tofHeader if missing in ESDs (2010 data)
//             - re-apply PID 
//             
// Contacts: Pietro.Antonioli@bo.infn.it                                     //
//           Francesco.Noferini@bo.infn.it                                   //
///////////////////////////////////////////////////////////////////////////////
#include <TFile.h>
#include <TChain.h>

#include <TMath.h>
#include <TRandom.h>
#include <AliLog.h>
#include <AliESDEvent.h>
#include <AliESDtrack.h>
#include <AliESDInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliESDpid.h>
#include <AliTender.h>

#include <AliTOFcalib.h>
#include <AliTOFT0maker.h>

#include <AliGeomManager.h>
#include <AliCDBManager.h>
#include <AliCDBEntry.h>

#include <AliOADBContainer.h>
#include <AliTOFPIDParams.h>

#include <AliT0CalibSeasonTimeShift.h>

#include "AliTOFTenderSupply.h"

ClassImp(AliTOFTenderSupply)

Float_t AliTOFTenderSupply::fgT0Aresolution = 75.;
Float_t AliTOFTenderSupply::fgT0Cresolution = 65.;

AliTOFTenderSupply::AliTOFTenderSupply() :
  AliTenderSupply(),
  fESDpid(0x0),
  fTenderNoAction(kTRUE),
  fIsMC(kFALSE),
  fCorrectExpTimes(kTRUE),
  fCorrectTRDBug(kFALSE),
  fLHC10dPatch(kFALSE),
  fT0DetectorAdjust(kFALSE),
  fDebugLevel(0),
  fAutomaticSettings(kTRUE),
  fRecoPass(0),
  fUserRecoPass(0),
  fForceCorrectTRDBug(kFALSE),
  fTOFPIDParams(0x0),
  fTOFCalib(0x0),
  fTOFT0maker(0x0),
  fT0IntercalibrationShift(0),
  fGeomSet(kFALSE),
  fIsEnteringInTRD(kFALSE),
  fInTRD(kFALSE),
  fIsComingOutTRD(kFALSE),
  fOutTRD(kFALSE),
  fRhoTRDin(288.38), // cm
  fRhoTRDout(366.38), // cm
  fStep(0.5),
  fMagField(0.),
  fCDBkey(0)



{
  //
  // default ctor
  //
  fT0shift[0] = 0;
  fT0shift[1] = 0;
  fT0shift[2] = 0;
  fT0shift[3] = 0;
}

//_____________________________________________________
AliTOFTenderSupply::AliTOFTenderSupply(const char *name, const AliTender *tender) :
  AliTenderSupply(name,tender),
  fESDpid(0x0),
  fTenderNoAction(kTRUE),
  fIsMC(kFALSE),
  fCorrectExpTimes(kTRUE),
  fCorrectTRDBug(kFALSE),
  fLHC10dPatch(kFALSE),
  fT0DetectorAdjust(kFALSE),
  fDebugLevel(0),
  fAutomaticSettings(kTRUE),
  fRecoPass(0),
  fUserRecoPass(0),
  fForceCorrectTRDBug(kFALSE),
  fTOFPIDParams(0x0),
  fTOFCalib(0x0),
  fTOFT0maker(0x0),
  fT0IntercalibrationShift(0),
  fGeomSet(kFALSE),
  fIsEnteringInTRD(kFALSE),
  fInTRD(kFALSE),
  fIsComingOutTRD(kFALSE),
  fOutTRD(kFALSE),
  fRhoTRDin(288.38), // cm
  fRhoTRDout(366.38), // cm
  fStep(0.5),
  fMagField(0.),
  fCDBkey(0)
 
{
  //
  // named ctor
  //

  fT0shift[0] = 0;
  fT0shift[1] = 0;
  fT0shift[2] = 0;
  fT0shift[3] = 0;
}

//_____________________________________________________
void AliTOFTenderSupply::Init()
{

  fTenderNoAction = kFALSE;
  // Initialise TOF tender (this is called at each detected run change)
  AliLog::SetClassDebugLevel("AliTOFTenderSupply",10); 

  // Setup PID object, check for MC, set AliTOFcalib and TOFT0 maker conf
  Int_t run = fTender->GetRun();
  if (run == 0) return;                // to skip first init, when we don't have yet a run number

  fGeomSet=kFALSE;                         
  // Even if the user didn't set fIsMC, we force it on if we find the MC handler  
  AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();
  if (mgr->GetMCtruthEventHandler() && !(fIsMC) ) {
    AliWarning("This ESD is MC, fIsMC found OFF: fIsMC turned ON");
    fIsMC=kTRUE;
  }

  if (fAutomaticSettings) {
    if (!fIsMC) {
      if (fUserRecoPass == 0) DetectRecoPass();
      else fRecoPass = fUserRecoPass;
    }
    if (run<114737) {
      fTenderNoAction = kTRUE;
    }
    else if (run>=114737&&run<=117223) {      //period="LHC10B";
      if (fRecoPass == 2) {fCorrectExpTimes=kTRUE; fCorrectTRDBug=kFALSE;}
      else if (fRecoPass == 3) {fCorrectExpTimes=kFALSE; fCorrectTRDBug=kTRUE;}
      fLHC10dPatch=kFALSE;
      fT0IntercalibrationShift = 0;
      fT0DetectorAdjust=kFALSE;   // it was kTRUE
    }
    else if (run>=118503&&run<=121040) { //period="LHC10C";
      if (fRecoPass == 2) {fCorrectExpTimes=kTRUE; fCorrectTRDBug=kFALSE;}
      else if (fRecoPass == 3) {fCorrectExpTimes=kFALSE; fCorrectTRDBug=kTRUE;}
      fLHC10dPatch=kFALSE;
      fT0IntercalibrationShift = 0;
      fT0DetectorAdjust=kFALSE;
    }
    else if (run>=122195&&run<=126437) { //period="LHC10D";
      fCorrectExpTimes=kFALSE;
      fLHC10dPatch=kTRUE;
      fT0IntercalibrationShift = 0;
      fT0DetectorAdjust=kFALSE;     // it was kTRUE
    }
    else if (run>=127719&&run<=130850) { //period="LHC10E";
      fCorrectExpTimes=kFALSE;
      fLHC10dPatch=kFALSE;
      fT0IntercalibrationShift = 30.;
      fT0DetectorAdjust=kFALSE;      // it was kTRUE
    }
    else if (run>=133004&&run<=135029) { //period="LHC10F";
      fTenderNoAction=kTRUE;
    }
    else if (run>=135654&&run<=136377) { //period="LHC10G";
      fTenderNoAction=kTRUE;
    }
    else if (run>=136851&&run<=139517) { //period="LHC10H" - pass2;
      fCorrectExpTimes=kFALSE;
      fLHC10dPatch=kFALSE;                
      fT0IntercalibrationShift = 0.;
      fT0DetectorAdjust=kFALSE;      // it was kTRUE
    }
    else if (run>=139699) {              //period="LHC11A";
      fTenderNoAction=kTRUE;
    }
  }

  if (fTenderNoAction) {
    AliInfo(" |---------------------------------------------------------------------------|");
    AliInfo(" |                                                                           |");
    AliInfo(Form(" |  TOF tender is not supported for run %d                           |",run));
    AliInfo(" |  TOF tender will do nothing.                                              |");
    AliInfo(" |  Check TOF tender usage for run/periods at:                               |");
    AliInfo(" |  https://twiki.cern.ch/twiki/bin/view/ALICE/TOF.                          |");
    AliInfo(" |---------------------------------------------------------------------------|");
    AliInfo(" ");
    return;
  }


  // Load from OADB TOF resolution
  LoadTOFPIDParams(run);

  // Check if another tender wagon already created the esd pid object
  // if not we create it and set it to the ESD input handler
  fESDpid=fTender->GetESDhandler()->GetESDpid();
  if (!fESDpid) {
    fESDpid=new AliESDpid;
    fTender->GetESDhandler()->SetESDpid(fESDpid);
  }


  // Configure TOF calibration class
  if (!fTOFCalib)fTOFCalib=new AliTOFcalib();  // create if needed
  fTOFCalib->SetRemoveMeanT0(!(fIsMC));        // must be kFALSE on MC (default is kTRUE)
  fTOFCalib->SetCalibrateTOFsignal(!(fIsMC));  // must be kFALSE on MC (no new calibration) (default is kTRUE)
  fTOFCalib->SetCorrectTExp(fCorrectExpTimes); // apply a fine tuning on the expected times at low momenta
                                               // (this is done for LHC10b/LHC10c pass2)

  // Configure TOFT0 maker class
  //  if (!fTOFT0maker) fTOFT0maker = new AliTOFT0maker(fESDpid,fTOFCalib); // create if needed
  if (!fTOFT0maker) fTOFT0maker = new AliTOFT0maker(fESDpid); // without passing AliTOFCalib it uses the diamond
  fTOFT0maker->SetTimeResolution(fTOFPIDParams->GetTOFresolution());     // set TOF resolution for the PID
  fTOFT0maker->SetTOFT0algorithm(2);

  AliInfo("|******************************************************|");
  AliInfo(Form("|    Alice TOF Tender Initialisation (Run %d)  |",fTender->GetRun()));
  AliInfo("|    Settings:                                         |");
  AliInfo(Form("|    Correct Exp Times              :  %d               |",fCorrectExpTimes));
  AliInfo(Form("|    Correct TRD Bug                :  %d               |",fCorrectTRDBug));
  AliInfo(Form("|    LHC10d patch                   :  %d               |",fLHC10dPatch));
  AliInfo(Form("|    TOF resolution for TOFT0 maker :  %5.2f (ps)     |",fTOFPIDParams->GetTOFresolution()));
  AliInfo(Form("|    MC flag                        :  %d               |",fIsMC));
  AliInfo(Form("|    T0 detector offsets applied    :  %d               |",fT0DetectorAdjust));
  AliInfo(Form("|    TOF/T0 intecalibration shift   :  %5.2f (ps)     |",fT0IntercalibrationShift));
  AliInfo("|******************************************************|");


}

//_____________________________________________________
void AliTOFTenderSupply::ProcessEvent()
{
  //
  // Use updated calibrations for TOF and T0, reapply PID information
  // For MC: timeZero sampling and additional smearing for T0

  if (fDebugLevel > 1) AliInfo("process event");

  AliESDEvent *event=fTender->GetEvent();
  if (!event) return;
  if (fDebugLevel > 1) AliInfo("event read");


    
  if (fTender->RunChanged()){ 

    Init();
    if (fTenderNoAction) return;            
    Int_t versionNumber = GetOCDBVersion(fTender->GetRun());
    fTOFCalib->SetRunParamsSpecificVersion(versionNumber);
    fTOFCalib->Init(fTender->GetRun());
    
    if(event->GetT0TOF()){ // read T0 detector correction from OCDB
	// OCDB instance
	if (fT0DetectorAdjust) {
	  AliCDBManager* ocdbMan = AliCDBManager::Instance();
	  ocdbMan->SetRun(fTender->GetRun());    
	  AliCDBEntry *entry = ocdbMan->Get("T0/Calib/TimeAdjust/");
	  if(entry) {
	    AliT0CalibSeasonTimeShift *clb = (AliT0CalibSeasonTimeShift*) entry->GetObject();
	    Float_t *t0means= clb->GetT0Means();
	    //      Float_t *t0sigmas = clb->GetT0Sigmas();
	    fT0shift[0] = t0means[0] + fT0IntercalibrationShift;
	    fT0shift[1] = t0means[1] + fT0IntercalibrationShift;
	    fT0shift[2] = t0means[2] + fT0IntercalibrationShift;
	    fT0shift[3] = t0means[3] + fT0IntercalibrationShift;
	  } else {
	    for (Int_t i=0;i<4;i++) fT0shift[i]=0;
	    AliWarning("TofTender no T0 entry found T0shift set to 0");
	  }
	} else {
	  for (Int_t i=0;i<4;i++) fT0shift[i]=0;
	}
    }
  }

  if (fTenderNoAction) return;

  fTOFCalib->CalibrateESD(event);   //recalculate TOF signal (no harm for MC, see settings inside init)


  // patches for various reconstruction bugs
  if (fLHC10dPatch && !(fIsMC)) RecomputeTExp(event);   // LHC10d pass2: fake full TRD geometry
  if ( (fCorrectTRDBug && !(fIsMC)) || (fForceCorrectTRDBug)) FixTRDBug(event);     // LHC10b,c pass3: wrong TRD dE/dx 

  Double_t startTime = 0.;
  if (fIsMC) startTime = fTOFCalib->TuneForMC(event,fTOFPIDParams->GetTOFresolution());   // this is for old MC when we didn't jitter startTime in MC

  if (fDebugLevel > 1) Printf(" TofTender: startTime %f",startTime);
  if (fDebugLevel > 1) Printf(" TofTender: T0 time (orig) %f %f %f",event->GetT0TOF(0),event->GetT0TOF(1),event->GetT0TOF(2));

  // event by event TO detector treatment  
  if(event->GetT0TOF()){   // protection: we adjust T0 only if it is there....

    if (event->GetT0TOF(0) == 0) event->SetT0TOF(0, 9999999.); // in case no information we set to unknown
    if (event->GetT0TOF(1) == 0) event->SetT0TOF(1, 99999.);
    if (event->GetT0TOF(2) == 0) event->SetT0TOF(2, 99999.);

    if (fT0DetectorAdjust) {
      if(!fIsMC){   // data: apply shifts to align around Zero
	event->SetT0TOF(0,event->GetT0TOF(0) - fT0shift[0]);
	event->SetT0TOF(1,event->GetT0TOF(1) - fT0shift[1]);
	event->SetT0TOF(2,event->GetT0TOF(2) - fT0shift[2]);
      } else {
      // MC: add smearing for realistic T0A and T0C resolution
	Double_t defResolutionT0A = 33.;   // in future we will get this from ESDrun data structure or via OCDB
	Double_t defResolutionT0C = 30.;   // for the moment we don't trust them
	if ( (fgT0Aresolution > defResolutionT0A) && (event->GetT0TOF(1)<90000.) ) { // add smearing only if signal is there
	  Double_t addedSmearingT0A = TMath::Sqrt(fgT0Aresolution*fgT0Aresolution - defResolutionT0A*defResolutionT0A);
	  Double_t smearingT0A = gRandom->Gaus(0.,addedSmearingT0A);
	  event->SetT0TOF(1,event->GetT0TOF(1) + smearingT0A);
	}
	if ( (fgT0Cresolution > defResolutionT0C) && (event->GetT0TOF(2)<90000.) ) { // add smearing only if signal is there
	Double_t addedSmearingT0C = TMath::Sqrt(fgT0Cresolution*fgT0Cresolution - defResolutionT0C*defResolutionT0C);
	Double_t smearingT0C = gRandom->Gaus(0.,addedSmearingT0C);
        event->SetT0TOF(2,event->GetT0TOF(2) + smearingT0C);
	}
	if (event->GetT0TOF(0)<90000.) { // we recompute the AND only if it is already there...
	  Double_t smearedT0AC = (event->GetT0TOF(1)+event->GetT0TOF(2))/2.;
	  event->SetT0TOF(0,smearedT0AC); 
	}
	if (fDebugLevel > 1) Printf(" TofTender: T0 time (postSmear) %f %f %f",event->GetT0TOF(0),event->GetT0TOF(1),event->GetT0TOF(2));
	// add finally the timeZero offset also to the T0 detector information
	event->SetT0TOF(0,event->GetT0TOF(0) + startTime);
	event->SetT0TOF(1,event->GetT0TOF(1) + startTime);
	event->SetT0TOF(2,event->GetT0TOF(2) + startTime);  
	if (fDebugLevel > 1) Printf(" TofTender: T0 time (postStart) %f %f %f",event->GetT0TOF(0),event->GetT0TOF(1),event->GetT0TOF(2));
      }
    }
    // after shifts adjust (data) or smearing+offset (MC) we 'clean' to default if signals not there 
    if(event->GetT0TOF(0) > 900000) event->SetT0TOF(0, 999999.);
    if(event->GetT0TOF(1) > 90000)  event->SetT0TOF(1, 99999.);
    if(event->GetT0TOF(2) > 90000)  event->SetT0TOF(2, 99999.);
  }
  if (fDebugLevel > 1) Printf(" TofTender: T0 time (FINAL) %f %f %f",event->GetT0TOF(0),event->GetT0TOF(1),event->GetT0TOF(2));
  
  //compute timeZero of the event via TOF-TO
  fTOFT0maker->ComputeT0TOF(event);
  fTOFT0maker->WriteInESD(event);

  //  set preferred startTime: this is now done via AliPIDResponseTask
  fESDpid->SetTOFResponse(event, (AliESDpid::EStartTimeType_t)fTOFPIDParams->GetStartTimeMethod());

  // recalculate PID probabilities
  // this is for safety, especially if the user doesn't attach a PID tender after TOF tender  
  Int_t ntracks=event->GetNumberOfTracks();
  //  AliESDtrack *track = NULL;
  //  Float_t tzeroTrack = 0;
  for(Int_t itrack = 0; itrack < ntracks; itrack++){
    //    track = event->GetTrack(itrack);
    //    tzeroTrack = fESDpid->GetTOFResponse().GetStartTime(track->P());
    //    Printf("================> Track # %d mom: %f tzeroTrack %f",itrack,track->P(),tzeroTrack);
    fESDpid->MakeTOFPID(event->GetTrack(itrack),0);   
  }
  
  
}


//_____________________________________________________
void AliTOFTenderSupply::RecomputeTExp(AliESDEvent *event) const
{
  /*
   * calibrate TExp
   */

  
  /* loop over tracks */
  AliESDtrack *track = NULL;
  for (Int_t itrk = 0; itrk < event->GetNumberOfTracks(); itrk++) {
    /* get track and calibrate */
    track = event->GetTrack(itrk);
    RecomputeTExp(track);
  }
  
}

//_____________________________________________________
void AliTOFTenderSupply::RecomputeTExp(AliESDtrack *track) const
{
  /*** 
       THIS METHOD IS BASED ON THEORETICAL EXPECTED TIME COMPUTED
       USING AVERAGE MOMENTUM BETWEEN INNER/OUTER TRACK PARAMS 
       IT IS A ROUGH APPROXIMATION APPLIED TO FIX LHC10d-pass2 DATA
       WHERE A WRONG GEOMETRY (FULL TRD) WAS INSERTED
  ***/

  Double_t texp[AliPID::kSPECIES];
  if (!track || !(track->GetStatus() & AliESDtrack::kTOFout)) return;


  /* get track params */
  Float_t l = track->GetIntegratedLength();
  Float_t p = track->P();
  if (track->GetInnerParam() && track->GetOuterParam()) {
    Float_t pin = track->GetInnerParam()->P();
    Float_t pout = track->GetOuterParam()->P();
    p = 0.5 * (pin + pout);
  }
  /* loop over particle types and compute expected time */
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++)
    texp[ipart] = GetExpTimeTh(AliPID::ParticleMass(ipart), p, l) - 37.; 
  // 37 is a final semiempirical offset to further adjust (calibrations were
  // done with "standard" integratedTimes)
  /* set integrated times */
  track->SetIntegratedTimes(texp);

}


//______________________________________________________________________________
void AliTOFTenderSupply::DetectRecoPass()
{
  //
  // Detect reconstruction information
  //
  
  //reset information
  fRecoPass=0;
  
  //Get the current file to check the reconstruction pass (UGLY, but not stored in ESD... )
  AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  if (!inputHandler) return;
  
  TTree *tree= (TTree*)inputHandler->GetTree();
  TFile *file= (TFile*)tree->GetCurrentFile();
  
  if (!file) {
    AliFatal("Current file not found");
    return; // coverity
  }
  
  //find pass from file name (UGLY, but not stored in ESD... )
  TString fileName(file->GetName());
  if (fileName.Contains("pass1") ) {
    fRecoPass=1;
  } else if (fileName.Contains("pass2") ) {
    fRecoPass=2;
  } else if (fileName.Contains("pass3") ) {
    fRecoPass=3;
  } else if (fileName.Contains("pass4") ) {
    fRecoPass=4;
  } else if (fileName.Contains("pass5") ) {
    fRecoPass=5;
  } else if (fileName.Contains("pass6") ) {
    fRecoPass=6;
  }
  if (fRecoPass == 0) {
    AliInfo(Form("From file name %s reco pass cannot be detected",fileName.Data()));
    AliInfo("Change file name or use SetUserRecoPass method");
    AliFatal("------------- TOF tender cannot run with reco pass unspecified, issuing FATAL error ---------- ");
  }
}


//______________________________________________________________________________
void AliTOFTenderSupply::InitGeom()
{

  if (fGeomSet == kTRUE) return;

  //  Printf("\n \n ----- calling InitGeom to fix TRD Bug ----- \n \n");
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  fCDBkey = man->SetLock(kFALSE, fCDBkey);
  Int_t run = fTender->GetRun();
  //  Printf(" ---------> run is %d",run);
  man->SetRun(run);
  AliCDBEntry *entry = (AliCDBEntry*)man->Get("GRP/Geometry/Data");
  if (entry) {
    AliGeomManager::LoadGeometry();
    AliGeomManager::ApplyAlignObjsFromCDB("ITS TPC TRD TOF");
    //    fCDBkey = man->SetLock(kTRUE, fCDBkey);
    //    Printf("\n \n ----- Geometry loaded ------ \n \n");
  }
  fGeomSet=kTRUE;

}


//______________________________________________________________________________
void AliTOFTenderSupply::FixTRDBug(AliESDEvent* event)
//
// recompute texp fixing wrong dE/dx from TRD (LHC10b,c pass3)
//
{

  if (fGeomSet == kFALSE) InitGeom();

  //  Printf("Running FixTRD bug ");
  /* loop over tracks */
  AliESDtrack *track = NULL;
  for (Int_t itrk = 0; itrk < event->GetNumberOfTracks(); itrk++) {
    track = event->GetTrack(itrk);
    FixTRDBug(track);
  }
}


//_____________________________________________________
void AliTOFTenderSupply::FixTRDBug(AliESDtrack *track)
{
  // 
  //
  //


    ULong_t status=track->GetStatus();
    if (!( ( (status & AliVTrack::kITSrefit)==AliVTrack::kITSrefit ) &&
	   ( (status & AliVTrack::kTPCrefit)==AliVTrack::kTPCrefit ) &&
	   ( (status & AliVTrack::kTPCout)==AliVTrack::kTPCout ) &&
	   ( (status & AliVTrack::kTOFout)==AliVTrack::kTOFout ) &&
	   ( (status & AliVTrack::kTIME)==AliVTrack::kTIME ) ) ) return;

    fIsEnteringInTRD=kFALSE;
    fInTRD=kFALSE;
    fIsComingOutTRD=kFALSE;
    fOutTRD=kFALSE;

    //    Printf("Track reached TOF %f",track->P());
    Double_t correctionTimes[5] = {0.,0.,0.,0.,0.}; // to be added to the expected times
    FindTRDFix(track, correctionTimes);
    Double_t expectedTimes[5] = {0.,0.,0.,0.,0.}; track->GetIntegratedTimes(expectedTimes);
    //    Printf("Exp. times: %f %f %f %f %f",
    //	   expectedTimes[0],expectedTimes[1],expectedTimes[2],expectedTimes[3],expectedTimes[4]);
    //    Printf("Corr. times: %f %f %f %f %f",
    //	   correctionTimes[0],correctionTimes[1],correctionTimes[2],correctionTimes[3],correctionTimes[4]);

    for (Int_t jj=0; jj<5; jj++) expectedTimes[jj]+=correctionTimes[jj];
    track->SetIntegratedTimes(expectedTimes);
  
}


//________________________________________________________________________
void AliTOFTenderSupply::FindTRDFix(AliESDtrack *track,Double_t *corrections)
{

  Double_t pT = track->Pt();
  ULong_t status=track->GetStatus();
  Bool_t isTRDout = (status & AliVTrack::kTRDout)==AliVTrack::kTRDout;

  Double_t length = 0.;

  Double_t xyzIN[3]={0.,0.,0.};
  fIsEnteringInTRD = track->GetXYZAt(fRhoTRDin,fMagField,xyzIN);

  Double_t xyzOUT[3]={0.,0.,0.};
  fIsComingOutTRD = track->GetXYZAt(fRhoTRDout,fMagField,xyzOUT);

  if (fIsEnteringInTRD && fIsComingOutTRD) {


    Double_t phiIN = TMath::Pi()+TMath::ATan2(-xyzIN[1],-xyzIN[0]);
    phiIN *= TMath::RadToDeg();
    fInTRD = ( (phiIN>=  0. && phiIN<= 40.) ||
	       (phiIN>=140. && phiIN<=220.) ||
	       (phiIN>=340. && phiIN<=360.) ); // TRD SMs installed @ 2010

    Double_t phiOUT = TMath::Pi()+TMath::ATan2(-xyzOUT[1],-xyzOUT[0]);
    phiOUT *= TMath::RadToDeg();
    fOutTRD = ( (phiOUT>=  0. && phiOUT<= 40.) ||
		(phiOUT>=140. && phiOUT<=220.) ||
		(phiOUT>=340. && phiOUT<=360.) ); // TRD SMs installed @ 2010

    length = 0.;

    if (fInTRD || fOutTRD) {

      if ( ( fInTRD && fOutTRD ) || ( fInTRD && !fOutTRD ) ) {
	length = EstimateLengthInTRD1(track);
      } else if ( !fInTRD && fOutTRD ) {
	length = EstimateLengthInTRD2(track);
      }

    } else { // ( !fInTRD && !fOutTRD )

      length = EstimateLengthOutTRD(track);

    }

  }
  //  Printf("estimated length in TRD %f [isTRDout %d]",length,isTRDout);
  CorrectDeltaTimes(pT,length,isTRDout,corrections);

}

//________________________________________________________________________
void AliTOFTenderSupply::CorrectDeltaTimes(Double_t pT,
						Double_t length,
						Bool_t flagTRDout,
						Double_t *corrections)
{

  corrections[2] = CorrectExpectedPionTime(pT,length,flagTRDout);
  corrections[0] = corrections[2]; // x electrons used pion corrections
  corrections[1] = corrections[2]; // x muons used pion corrections
  corrections[3] = CorrectExpectedKaonTime(pT,length,flagTRDout);
  corrections[4] = CorrectExpectedProtonTime(pT,length,flagTRDout);

}

//________________________________________________________________________
Double_t AliTOFTenderSupply::CorrectExpectedPionTime(Double_t pT,
							  Double_t length,
							  Bool_t isTRDout)
{
  // correction for expected time for pions

  Double_t delta=0.;

  //  Printf("Flags Ent In ComingOut Out %d %d %d %d",fIsEnteringInTRD,fInTRD,fIsComingOutTRD,fOutTRD);
  if (!fIsEnteringInTRD || !fIsComingOutTRD) { // zone 5

    Float_t p[2]={0.,0.};

    if (isTRDout) {

      if (pT<0.30) {
	p[0] = 180.; p[1] = 0.;
      } else if (pT>=0.30 && pT<0.35) {
	p[0] = 740.; p[1] = -1800.;
      } else if (pT>=0.35 && pT<0.40) {
	p[0] = 488.; p[1] =-1080.;
      } else if (pT>=0.40 && pT<0.45) {
	p[0] = 179.; p[1] = -307.;
      } else if (pT>=0.45 && pT<0.50) {
	p[0] =  97.; p[1] = -123.;
      } else { //if (pT>=0.50)
	p[0] = 120.; p[1] = -172.;
      }

    } else {

      if (pT<0.30) {
	p[0] =  70.; p[1] =    0.;
      } else if (pT>=0.30 && pT<0.35) {
	p[0] = 339.; p[1] = -927.;
      } else if (pT>=0.35 && pT<0.40) {
	p[0] =  59.; p[1] = -121.;
      } else if (pT>=0.40 && pT<0.50) {
	p[0] =  21.; p[1] =  -24.;
      } else { //if (pT>=0.50)
	p[0] =  42.; p[1] =  -67.;
      }

    }

    delta = p[0]+p[1]*pT;
    //    Printf("Pion time: %f %f %f %f",p[0],p[1],length,delta);

  } else {

    Float_t p[2] = {0.,0.};

    if ( fInTRD &&  fOutTRD) { // zone 1

      if (isTRDout) {

	if (length<130.) {
	  p[0] = 0.; p[1] = 0.;
	} else if (length>=130. && length<170.) {
	  p[0] = -20.5; p[1] = 0.25;
	} else {//if (length>=170.)
	  p[0] = 22.; p[1] = 0.;
	}

      } else { // !isTRDout

	p[0] = 20.; p[1] = 0.;

      }

    } else if (!fInTRD && !fOutTRD) { // zone 2

      p[0] = 0.; p[1] = 0.;

    } else if ( fInTRD &&  !fOutTRD) { // zone 3

      if (isTRDout) {

	if (length< 75.) {
	  p[0] = 17.; p[1] =  0.;
	} else if (length>= 75. && length< 95.) {
	  p[0] = 81.; p[1] = -0.85;
	} else if (length>= 95. && length<155.) {
	  p[0] =  0.; p[1] =  0.;
	} else {//if (length>=155.)
	  p[0] = 10.; p[1] =  0.;
	}

      } else { // !isTRDout

	p[0] = 0.; p[1] = 0.;

      }

    } else if (!fInTRD &&  fOutTRD) { // zone 4

      if (isTRDout) {

	if (length<80.) {
	  p[0] =  0.; p[1] =  0.;
	} else {//if (length>=80.)
	  p[0] = 10.; p[1] =  0.;
	}

      } else { // !isTRDout

	if (length<30.) {
	  p[0] =  0.; p[1] =  0.;
	} else {//if (length>=30.)
	  p[0] =  6.; p[1] =  0.;
	}

      }

    }

    delta = p[0]+p[1]*length;
    //    Printf("Pion time: %f %f %f %f",p[0],p[1],length,delta);

  }

  return delta;

}

//________________________________________________________________________
Double_t AliTOFTenderSupply::CorrectExpectedKaonTime(Double_t pT,
							  Double_t length,
							  Bool_t isTRDout)
{
  // correction for expected time for kaons

  Double_t delta=0.;
  //  Printf("Flags Ent In ComingOut Out %d %d %d %d",fIsEnteringInTRD,fInTRD,fIsComingOutTRD,fOutTRD);

  if (!fIsEnteringInTRD || !fIsComingOutTRD) { // zone 5

    Float_t p[2]={0.,0.};

    if (isTRDout) {

      if (pT<0.4) {
	p[0] =  900.; p[1] = 0.;
      } else if (pT>=0.40 && pT<0.45) {
	p[0] = 3100.; p[1] = -6000.;
      } else if (pT>=0.45 && pT<0.50) {
	p[0] = 1660.; p[1] = -2800.;
      } else if (pT>=0.50 && pT<0.55) {
	p[0] =  860.; p[1] = -1200.;
      } else { //if (pT>=0.55)
	p[0] =  200.; p[1] = 0.;
      }

    } else {

      if (pT<0.30) {
	p[0] =   0.; p[1] =    0.;
      } else if (pT>=0.30 && pT<0.32) {
	p[0] = 570.; p[1] =    0.;
      } else if (pT>=0.32 && pT<0.35) {
	p[0] = 3171.; p[1] = -8133.;
      } else if (pT>=0.35 && pT<0.40) {
	p[0] = 1815.; p[1] = -4260.;
      } else if (pT>=0.40 && pT<0.45) {
	p[0] =  715.; p[1] =  -1471.;
      } else if (pT>=0.45 && pT<0.50) {
	p[0] =  233.; p[1] =  -407.;
      } else if (pT>=0.50 && pT<0.55) {
	p[0] =  408.; p[1] =  -752.;
      } else { //if (pT>=0.55)
	p[0] =  408.-752.*0.55; p[1] = 0.;
      }

    }

    delta = p[0]+p[1]*pT;
    //    Printf("Kaon time: %f %f %f %f",p[0],p[1],length,delta);

  } else {

    Float_t p[2] = {0.,0.};

    if ( fInTRD &&  fOutTRD) { // zone 1

      if (isTRDout) {

	if (length<95.) {
	  p[0] = 20.; p[1] = 0.;
	} else if (length>=95. && length<195.) {
	  p[0] = -24.0; p[1] = 0.10+0.0041*length;
	} else {//if (length>=195.)
	  p[0] =  150.; p[1] = 0.;
	}

      } else { // !isTRDout

	p[0] = 40.; p[1] = 0.;

      }

    } else if (!fInTRD && !fOutTRD) { // zone 2

      p[0] = 0.; p[1] = 0.;

    } else if ( fInTRD &&  !fOutTRD) { // zone 3

      if (isTRDout) {

	if (length< 15.) {
	  p[0] = 180.; p[1] =  0.;
	} else if (length>= 15. && length< 55.) {
	  p[0] = 215.; p[1] = -2.5;
	} else {//if (length>=55.)
	  p[0] = 78.; p[1] =  0.;
	}

      } else { // !isTRDout

	p[0] = 0.; p[1] = 0.;

      }

    } else if (!fInTRD &&  fOutTRD) { // zone 4

      if (isTRDout) {

	if (length< 55.) {
	  p[0] =  0.; p[1] =  0.;
	} else if (length>= 55. && length<115.) {
	  p[0] = -85.; p[1] = 1.9;
	} else {//if (length>=115.)
	  p[0] = 100.; p[1] =  0.;
	}

      } else { // !isTRDout

	p[0] =  0.; p[1] =  0.;

      }

    }

    delta = p[0]+p[1]*length;
    //    Printf("Kaon time: %f %f %f %f",p[0],p[1],length,delta);

  }

  return delta;

}

//________________________________________________________________________
Double_t AliTOFTenderSupply::CorrectExpectedProtonTime(Double_t pT,
							    Double_t length,
							    Bool_t isTRDout)
{
  // correction for expected time for protons

  Double_t delta=0.;
  //  Printf("Flags Ent In ComingOut Out %d %d %d %d",fIsEnteringInTRD,fInTRD,fIsComingOutTRD,fOutTRD);

  if (!fIsEnteringInTRD || !fIsComingOutTRD) { // zone 5
    Float_t p[2]={0.,0.};


    if (isTRDout) {

      if (pT<0.375) {
	p[0] = 1000.; p[1] = 0.;
      }	else if (pT>=0.375 && pT<0.45) {
	p[0] = 1500.; p[1] = 0.;
      } else if (pT>=0.45 && pT<0.50) {
	p[0] = 4650.; p[1] = -7000.;
      } else if (pT>=0.50 && pT<0.55) {
	p[0] = 3150.; p[1] = -4000.;
      } else { //if (pT>=0.55)
	p[0] = 3150. -4000.*0.55; p[1] = 0.;
      }

    } else {

      if (pT<0.32) {
	p[0] = 2963.-5670.*0.032; p[1] = 0.;
      } else if (pT>=0.32 && pT<0.35) {
	p[0] = 2963.; p[1] =  -5670.;
      } else if (pT>=0.35 && pT<0.40) {
	p[0] = 4270.; p[1] =  -9400.;
      } else if (pT>=0.40 && pT<0.45) {
	p[0] = 1550.; p[1] =  -2600.;
      } else if (pT>=0.45 && pT<0.50) {
	p[0] = 1946.; p[1] =  -3480.;
      } else if (pT>=0.50 && pT<0.55) {
	p[0] = 1193.; p[1] =  -1974.;
      } else { //if (pT>=0.55)
	p[0] = 1193.-1974.*0.55; p[1] = 0.;
      }

    }

    delta = p[0]+p[1]*pT;
    //    Printf("Proton time: %f %f %f %f",p[0],p[1],length,delta);

  } else {

    Float_t p[2] = {0.,0.};

    if ( fInTRD &&  fOutTRD) { // zone 1

      if (isTRDout) {

	if (length<90.) {
	  p[0] = 0.; p[1] = 0.;
	} else if (length>=90. && length<200.) {
	  p[0] = 1063.; p[1] = -32.+0.30*length-0.00072*length*length;
	} else {//if (length>=200.)
	  p[0] =  900.; p[1] = 0.;
	}

      } else { // !isTRDout

	p[0] = 80.; p[1] = 0.;

      }

    } else if (!fInTRD && !fOutTRD) { // zone 2

      if (isTRDout) {
	p[0] = 0.; p[1] = 0.;
      } else {
	if (length<125.) {
	  p[0] = 0.; p[1] = 0.;
	} else if (length>=125. && length<180.) {
	  p[0] = -132.; p[1] = 1.3;
	} else {
	  p[0] = 100.; p[1] = 0.;
	}

      }

    } else if ( fInTRD &&  !fOutTRD) { // zone 3

      if (isTRDout) {

	if (length< 30.) {
	  p[0] = 670.; p[1] =  0.;
	} else if (length>= 30. && length<155.) {
	  p[0] = 944.; p[1] = -11.+0.064*length;
	} else {//if (length>=155.)
	  p[0] = 780.; p[1] =  0.;
	}

      } else { // !isTRDout

	if (length< 30.) {
	  p[0] = 140.; p[1] = -4.5;
	} else {
	  p[0] = 0.; p[1] = 0.;
	}
      
      }

    } else if (!fInTRD &&  fOutTRD) { // zone 4

      if (isTRDout) {

	if (length< 45.) {
	  p[0] = 130.; p[1] =  0.;
	} else if (length>= 45. && length<120.) {
	  p[0] = -190.; p[1] = 6.5;
	} else {//if (length>=120.)
	  p[0] = 750.; p[1] =  0.;
	}

      } else { // !isTRDout

	if (length<75.5) {
	  p[0] =    0.; p[1] =  0.;
	} else if (length>= 75.5 && length<90.) {
	  p[0] = -830.; p[1] = 11.;
	} else {
	  p[0] =  160.; p[1] =  0.;
	}

      }

    }

    delta = p[0]+p[1]*length;
    //    Printf("Proton time: %f %f %f %f",p[0],p[1],length,delta);

  }

  return delta;

}

//________________________________________________________________________
Double_t AliTOFTenderSupply::EstimateLengthInTRD1(AliESDtrack *track)
{

  Double_t xyz0[3]={0.,0.,0.};
  Bool_t stayInTRD = track->GetXYZAt(fRhoTRDin,fMagField,xyz0);

  Double_t phi0 = TMath::Pi()+TMath::ATan2(-xyz0[1],-xyz0[0]);
  phi0 *= TMath::RadToDeg();
  stayInTRD = stayInTRD && ( (phi0>=  0. && phi0<= 40.) ||
			     (phi0>=140. && phi0<=220.) ||
			     (phi0>=340. && phi0<=360.) );

  Double_t trackLengthInTRD = 0.;
  Int_t iStep=0;

  Double_t b[3];track->GetBxByBz(b);

  Double_t xyz1[3]={0.,0.,0.};
  Double_t rho = fRhoTRDin;
  while (stayInTRD && rho<=fRhoTRDout) {
    iStep++;
    rho += fStep;

    for (Int_t ii=0; ii<3; ii++) xyz1[ii]=0.;
    stayInTRD = track->GetXYZAt(rho,fMagField,xyz1);
    Double_t phi1 = TMath::Pi()+TMath::ATan2(-xyz1[1],-xyz1[0]);
    phi1 *= TMath::RadToDeg();
    stayInTRD = stayInTRD && ( (phi1>=  0. && phi1<= 40.) ||
			       (phi1>=140. && phi1<=220.) ||
			       (phi1>=340. && phi1<=360.) );

    Double_t l2  = TMath::Sqrt((xyz1[0]-xyz0[0])*(xyz1[0]-xyz0[0]) +
			       (xyz1[1]-xyz0[1])*(xyz1[1]-xyz0[1]) +
			       (xyz1[2]-xyz0[2])*(xyz1[2]-xyz0[2]));
    trackLengthInTRD += l2;

    for (Int_t ii=0; ii<3; ii++) xyz0[ii]=xyz1[ii];
  }

  return trackLengthInTRD;

}

//________________________________________________________________________
Double_t AliTOFTenderSupply::EstimateLengthInTRD2(AliESDtrack *track)
{

  Double_t xyz0[3]={0.,0.,0.};
  Bool_t stayInTRD = track->GetXYZAt(fRhoTRDout,fMagField,xyz0);

  Double_t phi0 = TMath::Pi()+TMath::ATan2(-xyz0[1],-xyz0[0]);
  phi0 *= TMath::RadToDeg();
  stayInTRD = stayInTRD && ( (phi0>=  0. && phi0<= 40.) ||
			     (phi0>=140. && phi0<=220.) ||
			     (phi0>=340. && phi0<=360.) );

  Double_t trackLengthInTRD = 0.;
  Int_t iStep=0;

  Double_t b[3];track->GetBxByBz(b);

  Double_t xyz1[3]={0.,0.,0.};
  Double_t rho = fRhoTRDout;
  while (stayInTRD && rho>=fRhoTRDin) {
    iStep++;
    rho -= fStep;

    for (Int_t ii=0; ii<3; ii++) xyz1[ii]=0.;
    stayInTRD = track->GetXYZAt(rho,fMagField,xyz1);
    Double_t phi1 = TMath::Pi()+TMath::ATan2(-xyz1[1],-xyz1[0]);
    phi1 *= TMath::RadToDeg();
    stayInTRD = stayInTRD && ( (phi1>=  0. && phi1<= 40.) ||
			       (phi1>=140. && phi1<=220.) ||
			       (phi1>=340. && phi1<=360.) );

    Double_t l2  = TMath::Sqrt((xyz0[0]-xyz1[0])*(xyz0[0]-xyz1[0]) +
			       (xyz0[1]-xyz1[1])*(xyz0[1]-xyz1[1]) +
			       (xyz0[2]-xyz1[2])*(xyz0[2]-xyz1[2]));
    trackLengthInTRD += l2;

    for (Int_t ii=0; ii<3; ii++) xyz0[ii]=xyz1[ii];
  }

  return trackLengthInTRD;

}

//________________________________________________________________________
Double_t AliTOFTenderSupply::EstimateLengthOutTRD(AliESDtrack *track)
{

  Double_t xyz0[3]={0.,0.,0.};
  Bool_t stayInTRD = track->GetXYZAt(fRhoTRDin,fMagField,xyz0);

  Double_t phi0 = TMath::Pi()+TMath::ATan2(-xyz0[1],-xyz0[0]);
  phi0 *= TMath::RadToDeg();
  stayInTRD = stayInTRD && !( (phi0>=  0. && phi0<= 40.) ||
			      (phi0>=140. && phi0<=220.) ||
			      (phi0>=340. && phi0<=360.) );

  Double_t trackLengthInTRD = 0.;
  Int_t iStep=0;

  Double_t b[3];track->GetBxByBz(b);

  Double_t xyz1[3]={0.,0.,0.};
  Double_t rho = fRhoTRDin;
  while (stayInTRD && rho<=fRhoTRDout) {
    iStep++;
    rho += fStep;

    for (Int_t ii=0; ii<3; ii++) xyz1[ii]=0.;
    stayInTRD = track->GetXYZAt(rho,fMagField,xyz1);
    Double_t phi1 = TMath::Pi()+TMath::ATan2(-xyz1[1],-xyz1[0]);
    phi1 *= TMath::RadToDeg();
    stayInTRD = stayInTRD && !( (phi1>=  0. && phi1<= 40.) ||
				(phi1>=140. && phi1<=220.) ||
				(phi1>=340. && phi1<=360.) );

    Double_t l2  = TMath::Sqrt((xyz1[0]-xyz0[0])*(xyz1[0]-xyz0[0]) +
			       (xyz1[1]-xyz0[1])*(xyz1[1]-xyz0[1]) +
			       (xyz1[2]-xyz0[2])*(xyz1[2]-xyz0[2]));
    trackLengthInTRD += l2;

    for (Int_t ii=0; ii<3; ii++) xyz0[ii]=xyz1[ii];
  }

  return trackLengthInTRD;

}

//________________________________________________________________________
Int_t AliTOFTenderSupply::GetOCDBVersion(Int_t runNo)
{
  Int_t verNo = -1;
  if ( (runNo>=118503 && runNo <=121040) )  { // LHC10C
    if (fRecoPass == 2) {                     // on pass2
      if (runNo >= 119159 && runNo <= 119163) verNo=3;
      else if (runNo >= 119837 && runNo <= 119934) verNo=4;
      else if (runNo >= 120067 && runNo <= 120244) verNo=4;
      else if (runNo >= 120503 && runNo <= 120505) verNo=4;
      else if (runNo >= 120616 && runNo <= 120671) verNo=4;
      else if (runNo >= 120741 && runNo <= 120829) verNo=4;
    }
  }
  return verNo;
}


//__________________________________________________________________________
void AliTOFTenderSupply::LoadTOFPIDParams(Int_t runNumber)
{
  //
  // Load the TOF pid params from the OADB
  //

  if (fTOFPIDParams) delete fTOFPIDParams;
  fTOFPIDParams=0x0;
  
  TFile *oadbf = new TFile("$ALICE_ROOT/OADB/COMMON/PID/data/TOFPIDParams.root");
  if (oadbf && oadbf->IsOpen()) {
    AliInfo("Loading TOF Params from $ALICE_ROOT/OADB/COMMON/PID/data/TOFPIDParams.root");
    AliOADBContainer *oadbc = (AliOADBContainer *)oadbf->Get("TOFoadb");
    if (oadbc) fTOFPIDParams = dynamic_cast<AliTOFPIDParams *>(oadbc->GetObject(runNumber,"TOFparams"));
    oadbf->Close();
    delete oadbc;
  }
  delete oadbf;

  if (!fTOFPIDParams) {
    AliError("TOFPIDParams.root not found in $ALICE_ROOT/OADB/COMMON/PID/data !!");
    fTOFPIDParams = new AliTOFPIDParams;
    fTOFPIDParams->SetTOFresolution(90.);
    fTOFPIDParams->SetStartTimeMethod(AliESDpid::kTOF_T0);
  }  
}
