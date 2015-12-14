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



/*
  Responsible: Raphaelle Bailhache (rbailhache@ikf.uni-frankfurt.de) 
  Code to analyze the TRD calibration and to produce OCDB entries  


   .x ~/rootlogon.C
   gSystem->Load("libANALYSIS");
   gSystem->Load("libTRDcalib");

   AliTRDPreprocessorOffline proces;
   TString ocdbPath="local:////"
   ocdbPath+=gSystem->GetFromPipe("pwd");

   proces.CalibTimeGain("CalibObjects.root",run0,run1,ocdbPath);
   proces.CalibTimeVdrift("CalibObjects.root",run0,run1,ocdbPath);
  // take the raw calibration data from the file CalibObjects.root 
  // and make a OCDB entry with run  validity run0-run1
  // results are stored at the ocdbPath - local or alien ...
  // default storage ""- data stored at current working directory 
 
*/
#include "AliLog.h"
#include "Riostream.h"
#include <fstream>
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH2I.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TMath.h"
#include "THnSparse.h"
#include "TProfile2D.h"
#include "AliTRDCalDet.h"
#include "AliTRDCalPad.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliTRDCalibraMode.h"
#include "AliTRDCalibraFit.h"
#include "AliTRDCalibraVdriftLinearFit.h"
#include "AliTRDCalibraExbAltFit.h"
#include "AliTRDPreprocessorOffline.h"
#include "AliTRDCalChamberStatus.h"
#include "AliTRDCalibChamberStatus.h"
#include "AliTRDCommonParam.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliTRDdEdxBaseUtils.h"
#include "AliTRDdEdxCalibHistArray.h"
#include "AliTRDdEdxCalibUtils.h"

ClassImp(AliTRDPreprocessorOffline)

  AliTRDPreprocessorOffline::AliTRDPreprocessorOffline():
  TNamed("TPCPreprocessorOffline","TPCPreprocessorOffline"),
  fMethodSecond(kTRUE),
  fNameList("TRDCalib"),
  fCalDetGainUsed(0x0),
  fCalDetVdriftUsed(0x0),
  fCalDetExBUsed(0x0),
  fCH2d(0x0),
  fPH2d(0x0),
  fPRF2d(0x0),
  fSparse(0x0),
  fAliTRDCalibraVdriftLinearFit(0x0),
  fAliTRDCalibraExbAltFit(0x0),
  fNEvents(0x0),
  fAbsoluteGain(0x0),
  fPlots(new TObjArray(kNumCalibObjs)),
  fCalibObjects(new TObjArray(kNumCalibObjs)),
  fFirstRunGainUsed(0),
  fVersionGainUsed(0),
  fSubVersionGainUsed(0),
  fFirstRunVdriftUsed(0),
  fVersionVdriftUsed(0), 
  fSubVersionVdriftUsed(0),
  fFirstRunExBUsed(0),
  fVersionExBUsed(0), 
  fSubVersionExBUsed(0),
  fNoExBUsedInReco(kFALSE),
  fSwitchOnValidation(kTRUE),
  fSwitchOnChamberStatus(kTRUE),
  fVdriftValidated(kFALSE),
  fExBValidated(kFALSE),
  fT0Validated(kFALSE),
  fMinStatsVdriftT0PH(800*20),
  fMinStatsVdriftLinear(800),
  fMinStatsGain(800),
  fMinStatsPRF(600),
  fMinStatsChamberStatus(20),
  fMinSingleStatsChamberStatus(0.05),
  fBackCorrectGain(kFALSE),  
  fBackCorrectVdrift(kTRUE),
  fNotEnoughStatisticsForTheGain(kFALSE),
  fNotEnoughStatisticsForTheVdriftLinear(kFALSE),
  fStatusNeg(0),
  fStatusPos(0),
  fBadCalibValidate(40),
  fNoDataValidate(40),
  fRMSBadCalibratedGain(15.0),
  fRMSBadCalibratedVdrift(20.0),
  fRMSBadCalibratedExB(20.0),
  fMinTimeOffsetValidate(-1.6),
  fRobustFitDriftVelocity(kTRUE),
  fRobustFitExbAlt(kFALSE),
  fAlternativeVdrfitFit(kFALSE),
  fAlternativeExbAltFit(kFALSE),
  fMinNbOfPointVdriftFit(11),
  fMethodeGain(0),
  fOutliersFitChargeLow(0.03),
  fOutliersFitChargeHigh(0.7),
  fBeginFitCharge(3.5),
  fT0Shift0(0.124797),
  fT0Shift1(0.267451),
  fMaxValueT0(5.),
  fPHQon(kTRUE),
  fDebugPHQon(kFALSE)
{
  //
  // default constructor
  //
  
  memset(fNotCalib, 0, sizeof(Int_t) * 18);
  memset(fNotGood, 0, sizeof(Int_t) * 18);
  memset(fBadCalib, 0, sizeof(Int_t) * 18);
  memset(fNoData, 0, sizeof(Int_t) * 18);
  memset(fNoDataA, 0, sizeof(Int_t) * 18);
  memset(fNoDataB, 0, sizeof(Int_t) * 18);
}
//_________________________________________________________________________________________________________________
AliTRDPreprocessorOffline::~AliTRDPreprocessorOffline() {
  //
  // Destructor
  //

  if(fCalDetGainUsed) delete fCalDetGainUsed;
  if(fCalDetVdriftUsed) delete fCalDetVdriftUsed;
  if(fCalDetExBUsed) delete fCalDetExBUsed;
  if(fCH2d) delete fCH2d;
  if(fPH2d) delete fPH2d;
  if(fPRF2d) delete fPRF2d;
  if(fSparse) delete fSparse;

  if(IsPHQon()){
    AliTRDdEdxCalibUtils::DeleteHistArray();
    AliTRDdEdxCalibUtils::DeleteObjArray();
  }

  if(fAliTRDCalibraVdriftLinearFit) delete fAliTRDCalibraVdriftLinearFit;
  if(fAliTRDCalibraExbAltFit) delete fAliTRDCalibraExbAltFit;
  if(fNEvents) delete fNEvents;
  if(fAbsoluteGain) delete fAbsoluteGain;
  if(fPlots) delete fPlots;
  if(fCalibObjects) delete fCalibObjects;
  
}
//___________________________________________________________________________________
void AliTRDPreprocessorOffline::Process(const Char_t* file, Int_t startRunNumber, Int_t endRunNumber, AliCDBStorage* ocdbStorage) 
{
  //
  // Process to the gain, vdrift, timeoffset, exb and chamber status calibration
  //
  
  if(SetCalDetGain(startRunNumber,fVersionGainUsed,fSubVersionGainUsed) && SetCalDetVdriftExB(startRunNumber,fVersionVdriftUsed,fSubVersionVdriftUsed,fVersionExBUsed,fSubVersionExBUsed)) {
    
    CalibVdriftT0(file,startRunNumber,endRunNumber,ocdbStorage);
    CalibGain(file,startRunNumber,endRunNumber,ocdbStorage);
    if(fSwitchOnChamberStatus) CalibChamberStatus(file,startRunNumber,endRunNumber,ocdbStorage);
    CalibExbAlt(file,startRunNumber,endRunNumber,ocdbStorage);

  }

  if(IsPHQon()){
    printf("\n                  AliTRDPreprocessorOffline PHQ on!!\n\n");
    AliTRDdEdxBaseUtils::PrintControl();
    CalibPHQ(file, startRunNumber, endRunNumber, ocdbStorage);
  }
  else{
    printf("\n                  AliTRDPreprocessorOffline PHQ off!!\n\n");
  }

  PrintStatus();
  
}
//___________________________________________________________________________________________________________________

void AliTRDPreprocessorOffline::CalibVdriftT0(const Char_t* file, Int_t startRunNumber, Int_t endRunNumber, AliCDBStorage* ocdbStorage){
  //
  // make calibration of the drift velocity
  // Input parameters:
  //      file                             - the location of input file
  //      startRunNumber, endRunNumber     - run validity period 
  //      ocdbStorage                      - OCDB storage
  //                                       - if empty - local storage 'pwd' uesed
  if (ocdbStorage==0x0) {
    TString localStorage = "local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
    ocdbStorage=AliCDBManager::Instance()->GetStorage(localStorage.Data());
  }
  //
  // 1. Initialization 
  //
  fVdriftValidated = kTRUE;
  fT0Validated = kTRUE;
  fExBValidated = kTRUE;
  fNotEnoughStatisticsForTheVdriftLinear = kFALSE;
  //
  // 2. extraction of the information
  //
  if(ReadVdriftLinearFitGlobal(file) && fCalDetVdriftUsed && fCalDetExBUsed) {
    if(!AnalyzeVdriftLinearFit()) return;
  }
  if(ReadVdriftT0Global(file)) {
    if(!AnalyzeVdriftT0()) return;
  }
  //
  // 3. Append QA plots
  //
  //MakeDefaultPlots(fVdriftArray,fVdriftArray);
  //
  //
  // 4. validate OCDB entries
  //
  if(fSwitchOnValidation==kTRUE && ValidateVdrift()==kFALSE) { 
    //AliError("TRD vdrift OCDB parameters out of range!");
    fVdriftValidated = kFALSE;
  }
  if(fSwitchOnValidation==kTRUE && ValidateT0()==kFALSE) { 
    //AliError("TRD t0 OCDB parameters out of range!");
    fT0Validated = kFALSE;
  }
  if(fSwitchOnValidation==kTRUE && ValidateExB()==kFALSE) { 
    //AliError("TRD t0 OCDB parameters out of range!");
    fExBValidated = kFALSE;
  }
  //
  // 5. update of OCDB
  //
  //
  if(fVdriftValidated) UpdateOCDBVdrift(startRunNumber,endRunNumber,ocdbStorage);
  if(fT0Validated) UpdateOCDBT0(startRunNumber,endRunNumber,ocdbStorage);
  if(fExBValidated) UpdateOCDBExB(startRunNumber,endRunNumber,ocdbStorage);
  
}
//___________________________________________________________________________________________________________________

void AliTRDPreprocessorOffline::CalibExbAlt(const Char_t* file, Int_t startRunNumber, Int_t endRunNumber, AliCDBStorage* ocdbStorage){
  //
  // make calibration of the drift velocity
  // Input parameters:
  //      file                             - the location of input file
  //      startRunNumber, endRunNumber     - run validity period 
  //      ocdbStorage                      - OCDB storage
  //                                       - if empty - local storage 'pwd' uesed
  if (ocdbStorage==0x0) {
    TString localStorage = "local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
    ocdbStorage=AliCDBManager::Instance()->GetStorage(localStorage.Data());
  }
  //
  // 1. Initialization 
  //

  //
  // 2. extraction of the information
  //
  if(ReadExbAltFitGlobal(file)) {
    if(!AnalyzeExbAltFit()) return;
  }
  //
  // 3. Append QA plots
  //
  //MakeDefaultPlots(fVdriftArray,fVdriftArray);
  //
  //
  // 4. validate OCDB entries
  //
  //
  // 5. update of OCDB
  //
  //
  UpdateOCDBExBAlt(startRunNumber,endRunNumber,ocdbStorage);
  
}

//_________________________________________________________________________________________________________________

void AliTRDPreprocessorOffline::CalibGain(const Char_t* file, Int_t startRunNumber, Int_t endRunNumber, AliCDBStorage* ocdbStorage){
  //
  // make calibration of the drift velocity
  // Input parameters:
  //      file                             - the location of input file
  //      startRunNumber, endRunNumber     - run validity period 
  //      ocdbStorage                      - OCDB storage
  //                                       - if empty - local storage 'pwd' uesed
  if (ocdbStorage==0x0) {
    TString localStorage = "local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
    ocdbStorage=AliCDBManager::Instance()->GetStorage(localStorage.Data());
  }
  //
  fNotEnoughStatisticsForTheGain = kFALSE;
  //
  // 1. Initialization 
  if(!ReadGainGlobal(file)) return;
  //
  //
  // 2. extraction of the information
  //
  if(!AnalyzeGain()) return;
  if(fBackCorrectGain) CorrectFromDetGainUsed();
  //if(fBackCorrectVdrift) CorrectFromDetVdriftUsed();
  //
  // 3. Append QA plots
  //
  //MakeDefaultPlots(fVdriftArray,fVdriftArray);
  //
  //
  // 4. validate OCDB entries
  //
  if(fSwitchOnValidation==kTRUE && ValidateGain()==kFALSE) { 
    //AliError("TRD gain OCDB parameters out of range!");
    return;
  }
  //
  // 5. update of OCDB
  //
  //
  if((!fCalDetVdriftUsed) || (fCalDetVdriftUsed && fVdriftValidated)) UpdateOCDBGain(startRunNumber,endRunNumber,ocdbStorage);
  
  
}
//________________________________________________________________________________________________________________

void AliTRDPreprocessorOffline::CalibPRF(const Char_t* file, Int_t startRunNumber, Int_t endRunNumber, AliCDBStorage* ocdbStorage){
  //
  // make calibration of the drift velocity
  // Input parameters:
  //      file                             - the location of input file
  //      startRunNumber, endRunNumber     - run validity period 
  //      ocdbStorage                      - OCDB storage
  //                                       - if empty - local storage 'pwd' uesed
  if (ocdbStorage==0x0) {
    TString localStorage = "local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
    ocdbStorage=AliCDBManager::Instance()->GetStorage(localStorage.Data());
  }
  //
  // 1. Initialization 
  if(!ReadPRFGlobal(file)) return;
  //
  //
  // 2. extraction of the information
  //
  if(!AnalyzePRF()) return;
  //
  // 3. Append QA plots
  //
  //MakeDefaultPlots(fVdriftArray,fVdriftArray);
  //
  //
  //
  // 4. validate OCDB entries
  //
  if(fSwitchOnValidation==kTRUE && ValidatePRF()==kFALSE) { 
    //AliError("TRD prf OCDB parameters out of range!");
    return;
  }
  //
  // 5. update of OCDB
  //
  //
  UpdateOCDBPRF(startRunNumber,endRunNumber,ocdbStorage);
  
}
//________________________________________________________________________________________________________________
void AliTRDPreprocessorOffline::CalibPHQ(const Char_t* file, Int_t startRunNumber, Int_t endRunNumber, AliCDBStorage* ocdbStorage)
{
  //
  // make calibration of puls height Q
  // Input parameters:
  //      startRunNumber, endRunNumber     - run validity period 
  //      ocdbStorage                      - OCDB storage
  //                                       - if empty - local storage 'pwd' uesed
  //

  if (ocdbStorage==0x0) {
    TString localStorage = "local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
    ocdbStorage=AliCDBManager::Instance()->GetStorage(localStorage.Data());
  }
  //printf("test %s\n", ocdbStorage.Data());

  if(!ReadPHQGlobal(file)) return;

  if(!AnalyzePHQ(startRunNumber)) return;

  UpdateOCDBPHQ(startRunNumber,endRunNumber,ocdbStorage);
}

//________________________________________________________________________________________________________________

void AliTRDPreprocessorOffline::CalibChamberStatus(const Char_t* file, Int_t startRunNumber, Int_t endRunNumber, AliCDBStorage* ocdbStorage){
  //
  // make calibration of the chamber status
  // Input parameters:
  //      startRunNumber, endRunNumber     - run validity period 
  //      ocdbStorage                      - OCDB storage
  //                                       - if empty - local storage 'pwd' uesed
  if (ocdbStorage==0x0) {
    TString localStorage = "local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
    ocdbStorage=AliCDBManager::Instance()->GetStorage(localStorage.Data());
  }
  //
  //
  // 1. Initialization  
  if(!ReadStatusGlobal(file)) return;
  //
  //
  //
  // 2. extraction of the information
  //
  if(!AnalyzeChamberStatus()) return;
  //
  // 3. Append QA plots
  //
  //MakeDefaultPlots(fVdriftArray,fVdriftArray);
  //
  //
  //
  // 4. validate OCDB entries
  //
  //printf("Enough stats for vdrift? %d\n",(Int_t)fNotEnoughStatisticsForTheVdriftLinear);
  //printf("Enough stats for gain? %d\n",(Int_t)fNotEnoughStatisticsForTheGain); 
  if((!fNotEnoughStatisticsForTheVdriftLinear) && (!fNotEnoughStatisticsForTheGain)) {
    if(fSwitchOnValidation==kTRUE && ValidateChamberStatus()==kFALSE) { 
      //AliError("TRD Chamber status OCDB parameters not ok!");
      return;
    }
    //
    // 5. update of OCDB
    //
    //
    UpdateOCDBChamberStatus(startRunNumber,endRunNumber,ocdbStorage);
  }
  
}
//______________________________________________________________________________________________________
Bool_t AliTRDPreprocessorOffline::Init(const Char_t* fileName){
  //
  // read the calibration used during the reconstruction
  // 

  if(ReadVdriftT0Global(fileName)) {
    
    TString nameph = fPH2d->GetTitle();
    fFirstRunVdriftUsed = GetFirstRun(nameph); 
    fVersionVdriftUsed = GetVersion(nameph);  
    fSubVersionVdriftUsed = GetSubVersion(nameph);    

    //printf("Found Version %d, Subversion %d for vdrift\n",fVersionVdriftUsed,fSubVersionVdriftUsed);
  
  }

  if(ReadGainGlobal(fileName)) {

    TString namech = fCH2d->GetTitle();
    fFirstRunGainUsed = GetFirstRun(namech); 
    fVersionGainUsed = GetVersion(namech);  
    fSubVersionGainUsed = GetSubVersion(namech);    

    //printf("Found Version %d, Subversion %d for gain\n",fVersionGainUsed,fSubVersionGainUsed);

  }
  
  if(ReadVdriftLinearFitGlobal(fileName)) {

    TString namelinear = fAliTRDCalibraVdriftLinearFit->GetNameCalibUsed();
    fFirstRunExBUsed = GetFirstRun(namelinear); 
    fVersionExBUsed = GetVersion(namelinear);  
    fSubVersionExBUsed = GetSubVersion(namelinear);   

    //printf("Found Version %d, Subversion %d, run %d for ExB\n",fVersionExBUsed,fSubVersionExBUsed,fFirstRunExBUsed);
    
  }
   
  if(fVersionVdriftUsed == 0) fStatusPos = fStatusPos |kVdriftErrorOld;
  if(fVersionGainUsed == 0) fStatusPos = fStatusPos | kGainErrorOld;

  return kTRUE;
  
}
//___________________________________________________________________________________________________________________

Bool_t AliTRDPreprocessorOffline::ReadStatusGlobal(const Char_t* fileName){
  //
  // read calibration entries from file
  // 
  if(fSparse) return kTRUE;
  TFile fcalib(fileName);
  TList * array = (TList*)fcalib.Get(fNameList);
  if (array){
    fSparse = (THnSparseI *) array->FindObject("NumberOfEntries");
    if(!fSparse) return kFALSE;
  }
  else 
    return kFALSE;
  
  return kTRUE;
  
}
//___________________________________________________________________________________________________________________

Bool_t AliTRDPreprocessorOffline::ReadPHQGlobal(const Char_t* fileName)
{
  //
  // read calibration entries from file
  //

  return AliTRDdEdxCalibUtils::ReadHistArray(fileName, fNameList);
}

//___________________________________________________________________________________________________________________

Bool_t AliTRDPreprocessorOffline::ReadGainGlobal(const Char_t* fileName){
  //
  // read calibration entries from file
  // 
  if(fCH2d) return kTRUE;
  TFile fcalib(fileName);
  TList * array = (TList*)fcalib.Get(fNameList);
  if (array){
    TH2I *ch2d = (TH2I *) array->FindObject("CH2d");
    if(!ch2d) {
      delete array;
      return kFALSE;
    }
    fCH2d = (TH2I*)ch2d->Clone();
    delete array;
    //fNEvents = (TH1I *) array->FindObject("NEvents");
    //fAbsoluteGain = (TH2F *) array->FindObject("AbsoluteGain");
  }else{
    TH2I *ch2d = (TH2I *) fcalib.Get("CH2d");
    if(!ch2d) return kFALSE;
    fCH2d = (TH2I*)ch2d->Clone();
    //fNEvents = (TH1I *) fcalib.Get("NEvents");
    //fAbsoluteGain = (TH2F *) fcalib.Get("AbsoluteGain");
  }
  fCH2d->SetDirectory(0);
  //printf("title of CH2d %s\n",fCH2d->GetTitle());

  return kTRUE;
  
}
//_________________________________________________________________________________________________________________

Bool_t AliTRDPreprocessorOffline::ReadVdriftT0Global(const Char_t* fileName){
  //
  // read calibration entries from file
  // 
  if(fPH2d) return kTRUE;
  TFile fcalib(fileName);
  TList * array = (TList*)fcalib.Get(fNameList);
  if (array){
    TProfile2D *ph2d = (TProfile2D *) array->FindObject("PH2d");
    if(!ph2d) {
      delete array;
      return kFALSE;
    }
    fPH2d = (TProfile2D*)ph2d->Clone();
    //fNEvents = (TH1I *) array->FindObject("NEvents");
    delete array;
  }else{
    TProfile2D *ph2d = (TProfile2D *) fcalib.Get("PH2d");
    if(!ph2d) return kFALSE;
    fPH2d = (TProfile2D*)ph2d->Clone();
    //fNEvents = (TH1I *) fcalib.Get("NEvents");
  }
  fPH2d->SetDirectory(0);
  //printf("title of PH2d %s\n",fPH2d->GetTitle());
  
  return kTRUE;
  
}
//___________________________________________________________________________________________________________________

Bool_t AliTRDPreprocessorOffline::ReadVdriftLinearFitGlobal(const Char_t* fileName){
  //
  // read calibration entries from file
  // 
  if(fAliTRDCalibraVdriftLinearFit) return kTRUE;
  TFile fcalib(fileName);
  TList * array = (TList*)fcalib.Get(fNameList);
  if (array){
    AliTRDCalibraVdriftLinearFit * dummy = (AliTRDCalibraVdriftLinearFit *) array->FindObject("AliTRDCalibraVdriftLinearFit");
    fAliTRDCalibraVdriftLinearFit = dummy ? (AliTRDCalibraVdriftLinearFit *) dummy->Clone() : 0x0;
    //fNEvents = (TH1I *) array->FindObject("NEvents");
    delete array;
  }else{
    fAliTRDCalibraVdriftLinearFit = (AliTRDCalibraVdriftLinearFit *) fcalib.Get("AliTRDCalibraVdriftLinearFit");
    //fNEvents = (TH1I *) fcalib.Get("NEvents");
  }
  if(!fAliTRDCalibraVdriftLinearFit) {
    //printf("No AliTRDCalibraVdriftLinearFit\n");
    return kFALSE;
  }
  return kTRUE;
  
}
//_____________________________________________________________________________________________________________
Bool_t AliTRDPreprocessorOffline::ReadExbAltFitGlobal(const Char_t* fileName){
  //
  // read calibration entries from file
  // 
  if(fAliTRDCalibraExbAltFit) return kTRUE;
  TFile fcalib(fileName);
  TList * array = (TList*)fcalib.Get(fNameList);
  if (array){
     AliTRDCalibraExbAltFit * dummy = (AliTRDCalibraExbAltFit *) array->FindObject("AliTRDCalibraExbAltFit");
     fAliTRDCalibraExbAltFit = dummy ? (AliTRDCalibraExbAltFit *)dummy->Clone() : 0x0;
    //fNEvents = (TH1I *) array->FindObject("NEvents");
    delete array;
  }else{
    fAliTRDCalibraExbAltFit = (AliTRDCalibraExbAltFit *) fcalib.Get("AliTRDCalibraExbAltFit");
    //fNEvents = (TH1I *) fcalib.Get("NEvents");
  }
  if(!fAliTRDCalibraExbAltFit) {
    //printf("No AliTRDCalibraExbAltFit\n");
    return kFALSE;
  }
  return kTRUE;
  
}
//_____________________________________________________________________________________________________________

Bool_t AliTRDPreprocessorOffline::ReadPRFGlobal(const Char_t* fileName){
  //
  // read calibration entries from file
  // 
  if(fPRF2d) return kTRUE;
  TFile fcalib(fileName);
  TList * array = (TList*)fcalib.Get(fNameList);
  if (array){
    TProfile2D *prf2d = (TProfile2D *) array->FindObject("PRF2d");
    if(!prf2d) {
      delete array;
      return kFALSE;
    }
    fPRF2d = (TProfile2D*)prf2d->Clone();
    delete array;
    //fNEvents = (TH1I *) array->FindObject("NEvents");
  }else{
    TProfile2D *prf2d = (TProfile2D *) fcalib.Get("PRF2d");
    if(!prf2d) return kFALSE;
    fPRF2d = (TProfile2D*)prf2d->Clone();
    //fNEvents = (TH1I *) fcalib.Get("NEvents");
  }
  fPRF2d->SetDirectory(0);
  //printf("title of PRF2d %s\n",fPRF2d->GetTitle());
  
  return kTRUE;

}
//__________________________________________________________________________________________________________

Bool_t AliTRDPreprocessorOffline::AnalyzeGain(){
  //
  // Analyze gain - produce the calibration objects
  //

  AliTRDCalibraFit *calibra = AliTRDCalibraFit::Instance();
  calibra->ChooseMethod(fMethodeGain);
  calibra->SetBeginFitCharge(fBeginFitCharge);
  calibra->SetFitOutliersChargeLow(fOutliersFitChargeLow);
  calibra->SetFitOutliersChargeHigh(fOutliersFitChargeHigh);
  calibra->SetMinEntries(fMinStatsGain); // If there is less than 1000 entries in the histo: no fit
  if(!calibra->AnalyseCH(fCH2d)) {
    return kFALSE;
  }

  Int_t nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(0))
    + 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(0));
  Int_t nbfit       = calibra->GetNumberFit();
  Int_t nbE         = calibra->GetNumberEnt();


  //Bool_t ok = kFALSE;
  Bool_t meanother = kFALSE;
  // enough statistics
  if ((nbtg >                  0) && 
      (nbfit        >= 0.5*nbE) && (nbE > 30)) {
    // create the cal objects
    if(!fBackCorrectGain) {
      calibra->PutMeanValueOtherVectorFit(1,kTRUE);
      meanother = kTRUE;
    }
    TObjArray object           = calibra->GetVectorFit();
    AliTRDCalDet *calDetGain   = calibra->CreateDetObjectGain(&object,meanother);
    TH1F *coefGain  = calDetGain->MakeHisto1DAsFunctionOfDet();
    // Put them in the array
    fCalibObjects->AddAt(calDetGain,kGain);
    fPlots->AddAt(coefGain,kGain);
    //
    //ok = kTRUE;
  }
  else {
    fNotEnoughStatisticsForTheGain = kTRUE;
    Int_t minStatsGain = fMinStatsGain*30;
    calibra->SetMinEntries(minStatsGain); // Because we do it for all, we increase this
    Double_t gainoverallnotnormalized =  calibra->AnalyseCHAllTogether(fCH2d);
    if(fCalDetGainUsed && (gainoverallnotnormalized > 0.0)) {
      AliTRDCalDet *calDetGain = new AliTRDCalDet(*fCalDetGainUsed);
      Int_t ndetu = 0;
      Double_t oldmean = fCalDetGainUsed->CalcMean(kFALSE,ndetu);
      //printf("oldmean %f and ndetu %f\n",oldmean,ndetu);
      if((oldmean > 0.0) && (ndetu>0))  {
      	Double_t scalefactor = calibra->GetScaleFactorGain();
	//printf("Correction factor %f\n",gainoverallnotnormalized*scalefactor);
	calDetGain->Multiply(gainoverallnotnormalized*scalefactor/oldmean);
	//printf("newmean %f\n",calDetGain->CalcMean(kFALSE));
	TH1F *coefGain  = calDetGain->MakeHisto1DAsFunctionOfDet();
	fCalibObjects->AddAt(calDetGain,kGain);
	fPlots->AddAt(coefGain,kGain);
	// 
	//ok = kTRUE;
	fStatusNeg = fStatusNeg | kGainNotEnoughStatsButFill;
      }
      else {
	fStatusPos = fStatusPos | kGainErrorOld;
      }      
    }
    else {
      if(gainoverallnotnormalized <= 0.0) fStatusNeg = fStatusNeg | kGainNotEnoughStatsNotFill;
      if(!fCalDetGainUsed) fStatusPos = fStatusPos | kGainErrorOld;
    }
  }
  
  calibra->ResetVectorFit();
  
  //return ok;
  return kTRUE;
  
}
//_____________________________________________________________________________________________________
Bool_t AliTRDPreprocessorOffline::AnalyzeVdriftT0(){
  //
  // Analyze VdriftT0 - produce the calibration objects
  //

  AliTRDCalibraFit *calibra = AliTRDCalibraFit::Instance();
  calibra->SetT0Shift0(fT0Shift0);
  calibra->SetT0Shift1(fT0Shift1);
  calibra->SetMaxValueT0(fMaxValueT0);
  calibra->SetMinEntries(fMinStatsVdriftT0PH); // If there is less than 1000 entries in the histo: no fit
  if(!calibra->AnalysePH(fPH2d)) {
    return kFALSE;
  }
  //calibra->SetDebugLevel(2);

  Int_t nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(1))
    + 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(1));
  Int_t nbfit       = calibra->GetNumberFit();
  Int_t nbfitSuccess = calibra->GetNumberFitSuccess();
  Int_t nbE         = calibra->GetNumberEnt();

  //printf("nbtg %d, nbfit %d, nbE %d, nbfitSuccess %d\n",nbtg,nbfit,nbE,nbfitSuccess);

  //Bool_t ok = kFALSE;
  if ((nbtg >                  0) && 
      (nbfit        >= 0.5*nbE) && (nbE > 30) && (nbfitSuccess > 30)) {
    //printf("nbtg %d, nbfit %d, nbE %d, nbfitSucess %d\n",nbtg,nbfit,nbE,nbfitSuccess);
    //printf("Pass the cut for VdriftT0\n");
    // create the cal objects
    calibra->RemoveOutliers(1,kFALSE);
    calibra->PutMeanValueOtherVectorFit(1,kFALSE);
    calibra->RemoveOutliers2(kFALSE);
    calibra->PutMeanValueOtherVectorFit2(1,kFALSE);
    //
    TObjArray object = calibra->GetVectorFit();
    AliTRDCalDet *calDetVdrift = calibra->CreateDetObjectVdrift(&object);
    TH1F *coefVdriftPH  = calDetVdrift->MakeHisto1DAsFunctionOfDet();
    AliTRDCalPad *calPadVdrift = (AliTRDCalPad *)calibra->CreatePadObjectVdrift(&object,calDetVdrift);
    TH1F *coefPadVdrift   = calPadVdrift->MakeHisto1D();
    object       = calibra->GetVectorFit2();
    AliTRDCalDet *calDetT0  = calibra->CreateDetObjectT0(&object);
    TH1F *coefT0  = calDetT0->MakeHisto1DAsFunctionOfDet();
    AliTRDCalPad *calPadT0 = (AliTRDCalPad *)calibra->CreatePadObjectT0(&object,calDetT0);
    TH1F *coefPadT0  = calPadT0->MakeHisto1D();
    // Put them in the array
    fCalibObjects->AddAt(calDetT0,kT0PHDet);
    fCalibObjects->AddAt(calDetVdrift,kVdriftPHDet);
    fCalibObjects->AddAt(calPadT0,kT0PHPad);
    fCalibObjects->AddAt(calPadVdrift,kVdriftPHPad);
    fPlots->AddAt(coefVdriftPH,kVdriftPHDet);
    fPlots->AddAt(coefT0,kT0PHDet);
    fPlots->AddAt(coefPadVdrift,kVdriftPHPad);
    fPlots->AddAt(coefPadT0,kT0PHPad);
    //
    //ok = kTRUE;
  }
  else {
    //printf("Not enough stats timeoffset\n");
    fStatusNeg = fStatusNeg | kTimeOffsetNotEnoughStatsNotFill;
  }
  calibra->ResetVectorFit();
 
  //return ok;

  return kTRUE;
  
}
//____________________________________________________________________________________________________________________
Bool_t AliTRDPreprocessorOffline::AnalyzeVdriftLinearFit(){
  //
  // Analyze vdrift linear fit - produce the calibration objects
  //

  //printf("Analyse linear fit\n");

  
  AliTRDCalibraFit *calibra = AliTRDCalibraFit::Instance();
  calibra->SetCalDetVdriftExB(fCalDetVdriftUsed,fCalDetExBUsed);
  calibra->SetMinEntries(fMinStatsVdriftLinear); // If there is less than 1000 entries in the histo: no fit
  //printf("The mean stat is by %d for VdriftLinear\n",fMinStatsVdriftLinear);
  //fAliTRDCalibraVdriftLinearFit->SetSeeDetector(0);
  //fAliTRDCalibraVdriftLinearFit->SetDebugLevel(1);
  //printf("Fill PE Array\n");
  fAliTRDCalibraVdriftLinearFit->SetRobustFit(fRobustFitDriftVelocity);
  fAliTRDCalibraVdriftLinearFit->SetMinNumberOfPointsForFit(fMinNbOfPointVdriftFit);
  if(!fAlternativeVdrfitFit)
    fAliTRDCalibraVdriftLinearFit->FillPEArray();
  else
    fAliTRDCalibraVdriftLinearFit->FillPEArray2();
  //printf("AliTRDCalibraFit\n");
  if(!calibra->AnalyseLinearFitters(fAliTRDCalibraVdriftLinearFit)) {
    return kFALSE;
  }
  //printf("After\n");

  //Int_t nbtg        = 540;
  Int_t nbfit       = calibra->GetNumberFit();
  Int_t nbE         = calibra->GetNumberEnt();

  
  //Bool_t ok = kFALSE;
  // enough statistics
  if ((nbfit        >= 0.5*nbE) && (nbE > 30)) {
    // create the cal objects
    //calibra->RemoveOutliers(1,kTRUE);
    calibra->PutMeanValueOtherVectorFit(1,kTRUE);
    //calibra->RemoveOutliers2(kTRUE);
    calibra->PutMeanValueOtherVectorFit2(1,kTRUE);
    //
    TObjArray object  = calibra->GetVectorFit();
    AliTRDCalDet *calDetVdrift = calibra->CreateDetObjectVdrift(&object,kTRUE);
    TH1F *coefDriftLinear      = calDetVdrift->MakeHisto1DAsFunctionOfDet();
    object                     = calibra->GetVectorFit2();
    AliTRDCalDet *calDetLorentz = calibra->CreateDetObjectLorentzAngle(&object);
    TH1F *coefLorentzAngle = calDetLorentz->MakeHisto1DAsFunctionOfDet();
    //if(!calDetLorentz) printf("No lorentz created\n");
    // Put them in the array
    fCalibObjects->AddAt(calDetVdrift,kVdriftLinear);
    fCalibObjects->AddAt(calDetLorentz,kLorentzLinear);
    fPlots->AddAt(coefDriftLinear,kVdriftLinear);
    fPlots->AddAt(coefLorentzAngle,kLorentzLinear);
    //
    //ok = kTRUE;
  }
  else {
    fNotEnoughStatisticsForTheVdriftLinear = kTRUE;
    Int_t minNumberOfEntriesForAll = fMinStatsVdriftLinear*30;
    calibra->SetMinEntries(minNumberOfEntriesForAll); // Because we do it for all, we increase this
    Double_t vdriftoverall = -100.0;
    Double_t exboverall = 100.0;
    calibra->AnalyseLinearFittersAllTogether(fAliTRDCalibraVdriftLinearFit,vdriftoverall,exboverall);
    //printf("Found mean vdrift %f and exb %f\n",vdriftoverall,exboverall);
    if(fCalDetVdriftUsed && (vdriftoverall > 0.0) && (exboverall < 70.0)) {
      AliTRDCalDet *calDetVdrift = new AliTRDCalDet(*fCalDetVdriftUsed);
      AliTRDCalDet *calDetLorentz = new AliTRDCalDet(*fCalDetExBUsed);
      Int_t ndetuv = 0;
      Int_t ndetue = 0;
      Double_t oldmeanvdrift = fCalDetVdriftUsed->CalcMean(kFALSE,ndetuv);
      Double_t oldmeanexb = fCalDetExBUsed->CalcMean(kFALSE,ndetue);
      Double_t oldmeanexbk = fCalDetExBUsed->GetMean();
      //printf("oldmeanvdrift %f and ndetuv %d\n",oldmeanvdrift,ndetuv);
      //printf("oldmeanexb %f and ndetue %d\n",oldmeanexb,ndetue);
      if(((oldmeanvdrift > 0.0) && (oldmeanexb < 70.0) && (ndetuv > 0) && (ndetue>0)) ||
	 ((oldmeanvdrift > 0.0) && (TMath::Abs(oldmeanexbk) < 0.0001) && (ndetuv > 0) && (ndetue==0)))  {
	//printf("Correction factor %f\n",vdriftoverall);
	calDetVdrift->Multiply(vdriftoverall/oldmeanvdrift);
	if(TMath::Abs(oldmeanexb) > 0.0001) calDetLorentz->Multiply(exboverall/oldmeanexb);
	//printf("newmean %f\n",calDetVdrift->CalcMean(kFALSE));
	TH1F *coefDriftLinear  = calDetVdrift->MakeHisto1DAsFunctionOfDet();
	TH1F *coefLorentzAngle = calDetLorentz->MakeHisto1DAsFunctionOfDet();
	// Put them in the array
	fCalibObjects->AddAt(calDetVdrift,kVdriftLinear);
	fCalibObjects->AddAt(calDetLorentz,kLorentzLinear);
	fPlots->AddAt(coefDriftLinear,kVdriftLinear);
	fPlots->AddAt(coefLorentzAngle,kLorentzLinear);
	// 
	//ok = kTRUE;
	fStatusNeg = fStatusNeg | kVdriftNotEnoughStatsButFill;
      }
      else {
	if(oldmeanvdrift) fStatusPos = fStatusPos | kVdriftErrorOld;
	if(oldmeanexb) fStatusPos = fStatusPos | kExBErrorOld;
      }      
    }
    else {
      if((vdriftoverall <= 0.0) && (exboverall > 70.0)) fStatusNeg = fStatusNeg | kVdriftNotEnoughStatsNotFill;
      if(!fCalDetVdriftUsed) fStatusPos = fStatusPos | kVdriftErrorOld;
      if(!fCalDetExBUsed) fStatusPos = fStatusPos | kExBErrorOld;
    }
  }
  
  calibra->ResetVectorFit();
  
  //return ok;
  return kTRUE;
  
}
//________________________________________________________________________________________________________________

Bool_t AliTRDPreprocessorOffline::AnalyzeExbAltFit(){
  //
  // Analyze vdrift linear fit - produce the calibration objects
  //

  //printf("Analyse linear fit\n");

  
  AliTRDCalibraFit *calibra = AliTRDCalibraFit::Instance();
  calibra->SetMinEntries(fMinStatsVdriftLinear); // If there is less than 1000 entries in the histo: no fit
  //printf("Fill PE Array\n");
  fAliTRDCalibraExbAltFit->SetRobustFit(fRobustFitExbAlt);
  if(!fAlternativeExbAltFit)
    fAliTRDCalibraExbAltFit->FillPEArray();
  else
    fAliTRDCalibraExbAltFit->FillPEArray2();
  //printf("AliTRDCalibraFit\n");
  if(!calibra->AnalyseExbAltFit(fAliTRDCalibraExbAltFit)) return kFALSE;
  //printf("After\n");

  //Int_t nbtg        = 540;
  Int_t nbfit       = calibra->GetNumberFit();
  Int_t nbE         = calibra->GetNumberEnt();

  
  //Bool_t ok = kFALSE;
  // enough statistics
  if ((nbfit        >= 0.5*nbE) && (nbE > 30)) {
    // create the cal objects
    //calibra->RemoveOutliers(1,kTRUE);
    calibra->PutMeanValueOtherVectorFit2(1,kTRUE);
    //
    TObjArray object  = calibra->GetVectorFit2();
    AliTRDCalDet *calDetLorentz = calibra->CreateDetObjectExbAlt(&object);
    TH1F *coefLorentzAngle = calDetLorentz->MakeHisto1DAsFunctionOfDet();
    //if(!calDetLorentz) printf("No lorentz created\n");
    // Put them in the array
    fCalibObjects->AddAt(calDetLorentz,kExbAlt);
    fPlots->AddAt(coefLorentzAngle,kExbAlt);
    //
    //ok = kTRUE;
  }
  
  calibra->ResetVectorFit();
  
  //return ok;

  return kTRUE;
  
}
//________________________________________________________________________________________________________________

Bool_t AliTRDPreprocessorOffline::AnalyzePRF(){
  //
  // Analyze PRF - produce the calibration objects
  //

  AliTRDCalibraFit *calibra = AliTRDCalibraFit::Instance();
  calibra->SetMinEntries(fMinStatsPRF); // If there is less than 1000 entries in the histo: no fit
  if(!calibra->AnalysePRFMarianFit(fPRF2d)) return kFALSE;

  Int_t nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(2))
    + 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(2));
  Int_t nbfit       = calibra->GetNumberFit();
  Int_t nbE         = calibra->GetNumberEnt();

  
  //Bool_t ok = kFALSE;
  // enough statistics
  if ((nbtg >                  0) && 
      (nbfit        >= 0.95*nbE) && (nbE > 30)) {
    // create the cal objects
    TObjArray object  = calibra->GetVectorFit();
    AliTRDCalPad *calPadPRF = (AliTRDCalPad*) calibra->CreatePadObjectPRF(&object);
    TH1F *coefPRF           = calPadPRF->MakeHisto1D();
    // Put them in the array
    fCalibObjects->AddAt(calPadPRF,kPRF);
    fPlots->AddAt(coefPRF,kPRF);
    //
    //ok = kTRUE;
  }
  
  calibra->ResetVectorFit();
  
  //return ok;
  return kFALSE;

  
}

//_____________________________________________________________________________
Bool_t AliTRDPreprocessorOffline::AnalyzePHQ(Int_t startRunNumber)
{
  //
  //Produce PHQ calibration results
  //
  TList *lout = 0x0;
  TTreeSRedirector *calibStream = 0x0;
  if(IsDebugPHQon()){
    lout = new TList;
    lout->SetOwner();

    calibStream = new TTreeSRedirector(Form("TRDCalibStream_%010d.root", startRunNumber));
  }

  for(Int_t iter=0; iter<AliTRDdEdxCalibUtils::GetHistArray()->GetSize(); iter++){
    THnBase *hi = (THnBase*) AliTRDdEdxCalibUtils::GetHistAt(iter);
    TObjArray *obji = AliTRDdEdxCalibUtils::HistToObj(hi, startRunNumber, lout, calibStream);
    //printf("test analyze %s\n", obji->GetName());
    AliTRDdEdxCalibUtils::GetObjArray()->AddAt(obji, iter);
  }

  fCalibObjects->AddAt(AliTRDdEdxCalibUtils::GetObjArray(), kPHQ);

  if(lout){
    TFile *fout=new TFile(Form("TRDCalibList_%010d.root", startRunNumber),"recreate");
    fout->cd();
    lout->Write();
    fout->Save();
    fout->Close();
    delete fout;
  }
  delete calibStream;
  delete lout;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliTRDPreprocessorOffline::AnalyzeChamberStatus()
{
  //
  // Produce AliTRDCalChamberStatus out of calibration results
  //
  
  // set up AliTRDCalibChamberStatus
  AliTRDCalibChamberStatus *chamberStatus = new AliTRDCalibChamberStatus();
  chamberStatus->SetSparseI(fSparse);
  chamberStatus->AnalyseHisto(fMinStatsChamberStatus, fMinSingleStatsChamberStatus);
  // get AliTRDCalChamberStatus
  AliTRDCalChamberStatus *calChamberStatus = chamberStatus->GetCalChamberStatus();

  // get calibration objects
  AliTRDCalDet *calDetGain   = (AliTRDCalDet *) fCalibObjects->At(kGain);
  AliTRDCalDet *calDetVDrift = (AliTRDCalDet *) fCalibObjects->At(kVdriftLinear);
  AliTRDCalDet *calDetExB    = (AliTRDCalDet *) fCalibObjects->At(kLorentzLinear);

  // Check
  if((!calDetGain) || (!calDetVDrift) || (!fCH2d) || (!calDetExB) || (!calChamberStatus)) return kFALSE;

  // Gain
  Double_t gainmean   = calDetGain->GetMean();
  Double_t vdriftmean = calDetVDrift->GetMean();
  Double_t exbmean    = calDetExB->GetMean();

  Double_t gainrms    = calDetGain->GetRMSRobust();
  Double_t vdriftrms  = calDetVDrift->GetRMSRobust();
  Double_t exbrms     = calDetExB->GetRMSRobust();

  //printf("Gain mean: %f, rms: %f\n",gainmean,gainrms);
  //printf("Vdrift mean: %f, rms: %f\n",vdriftmean,vdriftrms);
  //printf("ExB mean: %f, rms: %f\n",exbmean,exbrms);

  // Check
  if((TMath::Abs(gainrms) < 0.001) || (TMath::Abs(vdriftrms) < 0.001) || (TMath::Abs(exbrms) < 0.0000001)) return kFALSE;

  // Take mean each SM
  Double_t *gainmeanSM = new Double_t[18];
  Double_t *vdriftmeanSM = new Double_t[18];
  Double_t *exbmeanSM = new Double_t[18];
  //Double_t *t0meanSM = new Double_t[18];
  for(Int_t sm=0; sm< 18; sm++) {
    gainmeanSM[sm] = calDetGain->GetMeanSM(kFALSE,sm);
    vdriftmeanSM[sm] = calDetVDrift->GetMeanSM(kFALSE,sm);
    exbmeanSM[sm] = calDetExB->GetMeanSM(kFALSE,sm);
    //t0meanSM[sm] = calDetGain->GetMeanSM(kFALSE);
  }


  // mask chambers with empty gain entries
  //Int_t counter = 0;
  for (Int_t idet = 0; idet < 540; idet++) {

    // ch2d
    TH1I *projch =  (TH1I *) fCH2d->ProjectionX("projch",idet+1,idet+1,(Option_t *)"e");
    Double_t entries = projch->GetEntries();
    //printf("Number of entries %f for det %d\n",entries,idet);

    TVectorD error(3);
    Bool_t heree    = fAliTRDCalibraVdriftLinearFit->GetError(idet,&error);
    Double_t entriesvd = 0.;
    if(heree) {
      entriesvd = error[2];
    }
    //printf("Number of entries for detector %d for gain %f and for vd %f\n",idet,entries,entriesvd);
    
    // sm number
    Int_t smnumber = (Int_t) idet/30;

    // gain
    Double_t gain = calDetGain->GetValue(idet);

    // vdrift
    Double_t vdrift = calDetVDrift->GetValue(idet);

    // exb
    Double_t exb = calDetExB->GetValue(idet);

    
    if( (entries<50 && !calChamberStatus->IsNoData(idet))  ||
        TMath::Abs(gainmean-gain) > (fRMSBadCalibratedGain*gainrms)          ||
        TMath::Abs(vdriftmean-vdrift) > (fRMSBadCalibratedVdrift*vdriftrms)    ||
        TMath::Abs(exbmean-exb) > (fRMSBadCalibratedExB*exbrms) ) {
     
      //printf(" chamber det %03d masked \n",idet);
      //printf(" gainmean %f and gain %f, gainrms %f \n",gainmean,gain,gainrms);
      //printf(" vdriftmean %f and vdrift %f, vdriftrms %f \n",vdriftmean,vdrift,vdriftrms);
      //printf(" exbmean %f and exb %f, exbrms %f \n",exbmean,exb,exbrms);
      
      calChamberStatus->SetStatus(idet,AliTRDCalChamberStatus::kBadCalibrated);
      //counter++;
    }

    if(TMath::Abs(gainmeanSM[smnumber]-gain) < 0.000001  ||
       TMath::Abs(vdriftmeanSM[smnumber]-vdrift) < 0.000001 ||
       entries < fMinStatsGain                              ||
       entriesvd < fMinStatsVdriftLinear                    ||
       TMath::Abs(exbmeanSM[smnumber]-exb) < 0.000001) {
      
      //printf(" chamber det %03d notcalibrated sm %d \n",idet,smnumber);
      //printf(" gainmeanSM %f and gain %f\n",gainmeanSM[smnumber],gain);
      //printf(" vdriftmeanSM %f and vdrift %f \n",vdriftmeanSM[smnumber],vdrift);
      //printf(" exbmeanSM %f and exb %f \n",exbmeanSM[smnumber],exb);

      calChamberStatus->SetStatus(idet,AliTRDCalChamberStatus::kNotCalibrated);
    }


    delete projch;
    
   }

  // Security
  for(Int_t sm=0; sm < 18; sm++) {
    Int_t smnodata   = 0;
    Int_t smbadcalib = 0;
    Int_t smnodataA   = 0;
    Int_t smnodataB   = 0;
    Int_t smnotgood   = 0;
    Int_t smnotcalib  = 0;
    for(Int_t det = 0; det < 30; det++){
      Int_t detector = sm*30+det;
      if(calChamberStatus->IsNoData(detector)) smnodata++;
      else {
	if(calChamberStatus->IsBadCalibrated(detector)) smbadcalib++;
      }

      if(calChamberStatus->IsNoDataSideA(detector)) smnodataA++;
      if(calChamberStatus->IsNoDataSideB(detector)) smnodataB++;
      if(calChamberStatus->IsNotCalibrated(detector)) smnotcalib++;
      if(!calChamberStatus->IsGood(detector)) smnotgood++;
    }
    fNotGood[sm]  = smnotgood;
    fNotCalib[sm]  = smnotcalib;
    fNoDataA[sm]  = smnodataA;
    fNoDataB[sm]  = smnodataA;
    fNoData[sm]  = smnodata;
    fBadCalib[sm]= smbadcalib;
    //printf("No Data %d, No Data A %d, No Data B %d, bad calibrated %d, not calibrated %d and not good %d for %d\n",fNoData[sm],fNoDataA[sm],fNoDataB[sm],fBadCalib[sm],fNotCalib[sm],fNotGood[sm],sm);
  }

  // delete
  delete []gainmeanSM;
  delete []vdriftmeanSM;
  delete []exbmeanSM;
  
  // Security
  //   for(Int_t sm=0; sm < 18; sm++) {
  //     Int_t counter = 0;
  //     for(Int_t det = 0; det < 30; det++){
  //       Int_t detector = sm*30+det;
  //       if(calChamberStatus->IsBadCalibrated(detector)) counter++;
  //     }
  //     if(counter >= 20) {
  //       for(Int_t det = 0; det < 30; det++){
  // 	Int_t detector = sm*30+det;
  // 	calChamberStatus->SetStatus(detector,AliTRDCalChamberStatus::kGood);
  //       }
  //     }
  //  }

   fCalibObjects->AddAt(calChamberStatus,kChamberStatus);
   return kTRUE;

 }


 //________________________________________________________________________________________________
 void AliTRDPreprocessorOffline::CorrectFromDetGainUsed() {
   //
   // Correct from the gas gain used afterwards
   //
   AliTRDCalDet *calDetGain = (AliTRDCalDet *) fCalibObjects->At(kGain);
   if(!calDetGain) return;

   // Calculate mean
   Double_t mean = 0.0;
   Int_t nbdet = 0;

   for(Int_t det = 0; det < 540; det++) {

     Float_t gaininit = fCalDetGainUsed->GetValue(det);
     Float_t gainout = calDetGain->GetValue(det);


     if(TMath::Abs(gainout-1.0) > 0.000001) {
       mean += (gaininit*gainout);
       nbdet++;
     }  
   }
   if(nbdet > 0) mean = mean/nbdet;

   for(Int_t det = 0; det < 540; det++) {

     Float_t gaininit = fCalDetGainUsed->GetValue(det);
     Float_t gainout = calDetGain->GetValue(det);

     if(TMath::Abs(gainout-1.0) > 0.000001) {
       Double_t newgain = gaininit*gainout;
       if(newgain < 0.1) newgain = 0.1;
       if(newgain > 1.9) newgain = 1.9;
       calDetGain->SetValue(det,newgain);
     }
     else {
       Double_t newgain = mean;
       if(newgain < 0.1) newgain = 0.1;
       if(newgain > 1.9) newgain = 1.9;
       calDetGain->SetValue(det,newgain);
     }
   }


 }
 //________________________________________________________________________________________________
 void AliTRDPreprocessorOffline::CorrectFromDetVdriftUsed() {
   //
   // Correct from the drift velocity
   //

   //printf("Correct for vdrift\n");

   AliTRDCalDet *calDetGain = (AliTRDCalDet *) fCalibObjects->At(kGain);
   if(!calDetGain) return;

   Int_t detVdrift = kVdriftPHDet;
   if(fMethodSecond) detVdrift = kVdriftLinear;
   AliTRDCalDet *calDetVdrift = (AliTRDCalDet *) fCalibObjects->At(detVdrift);
   if(!calDetVdrift) return;

   // Calculate mean
   if(!fNotEnoughStatisticsForTheVdriftLinear) {
     for(Int_t det = 0; det < 540; det++) {
       
       Float_t vdriftinit = fCalDetVdriftUsed->GetValue(det);
       Float_t vdriftout = calDetVdrift->GetValue(det);
       
       Float_t gain = calDetGain->GetValue(det);
       if(vdriftout > 0.0) gain = gain*vdriftinit/vdriftout;
       if(gain < 0.1) gain = 0.1;
       if(gain > 1.9) gain = 1.9;
       calDetGain->SetValue(det,gain);
     }
   }
   else {
     
     Float_t vdriftinit = fCalDetVdriftUsed->CalcMean(kFALSE);
     Float_t vdriftout = calDetVdrift->CalcMean(kFALSE);
     Float_t factorcorrectif = 1.0;
     if(vdriftout > 0.0) factorcorrectif = vdriftinit/vdriftout;
     for(Int_t det = 0; det < 540; det++) {
       Float_t gain = calDetGain->GetValue(det);
       gain = gain*factorcorrectif;
       if(gain < 0.1) gain = 0.1;
       if(gain > 1.9) gain = 1.9;
       calDetGain->SetValue(det,gain);
     }
     
   }
   
 }
//_________________________________________________________________________________________________________________
 void AliTRDPreprocessorOffline::UpdateOCDBGain(Int_t startRunNumber, Int_t endRunNumber, AliCDBStorage *storage){
   //
   // Update OCDB entry
   //

   Bool_t status = kTRUE;

   AliCDBMetaData *metaData= new AliCDBMetaData();
   metaData->SetObjectClassName("AliTRDCalDet");
   metaData->SetResponsible("Raphaelle Bailhache");
   metaData->AddDateToComment();
   metaData->SetBeamPeriod(1);

   AliCDBId id1("TRD/Calib/ChamberGainFactor", startRunNumber, endRunNumber);
   AliTRDCalDet *calDet = (AliTRDCalDet *) fCalibObjects->At(kGain);
   if(calDet) status = storage->Put(calDet, id1, metaData);
   if (status==kFALSE) fStatusPos = fStatusPos | kCalibFailedExport;

 }
 //___________________________________________________________________________________________________________________
 void AliTRDPreprocessorOffline::UpdateOCDBExB(Int_t startRunNumber, Int_t endRunNumber, AliCDBStorage *storage){
   //
   // Update OCDB entry
   //

   Bool_t status = kTRUE;

   Int_t detExB = kLorentzLinear;
   if(!fMethodSecond) return;

   //printf("Pass\n");

   AliCDBMetaData *metaData= new AliCDBMetaData();
   metaData->SetObjectClassName("AliTRDCalDet");
   metaData->SetResponsible("Raphaelle Bailhache");
   metaData->AddDateToComment();
   metaData->SetBeamPeriod(1);

   AliCDBId id1("TRD/Calib/ChamberExB", startRunNumber, endRunNumber);
   AliTRDCalDet *calDet = (AliTRDCalDet *) fCalibObjects->At(detExB);
   if(calDet) status = storage->Put(calDet, id1, metaData);
   if (status==kFALSE) fStatusPos = fStatusPos | kCalibFailedExport;
   //if(!calDet) printf("No caldet\n");

 }
//___________________________________________________________________________________________________________________
void AliTRDPreprocessorOffline::UpdateOCDBExBAlt(Int_t startRunNumber, Int_t endRunNumber, AliCDBStorage *storage){
  //
  // Update OCDB entry
  //

  Bool_t status = kTRUE;

  Int_t detExB = kExbAlt;
  if(!fMethodSecond) return;

  //printf("Pass\n");

  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("AliTRDCalDet");
  metaData->SetResponsible("Theo Rascanu");
  metaData->AddDateToComment();
  metaData->SetBeamPeriod(1);

  AliCDBId id1("TRD/Calib/ChamberExBAlt", startRunNumber, endRunNumber);
  AliTRDCalDet *calDet = (AliTRDCalDet *) fCalibObjects->At(detExB);
  if(calDet) status = storage->Put(calDet, id1, metaData);
  if (status==kFALSE) fStatusPos = fStatusPos | kCalibFailedExport;
  //if(!calDet) printf("No caldet\n");

}
 //___________________________________________________________________________________________________________________
 void AliTRDPreprocessorOffline::UpdateOCDBVdrift(Int_t startRunNumber, Int_t endRunNumber, AliCDBStorage* storage){
   //
   // Update OCDB entry
   //

   Bool_t status = kTRUE;

   Int_t detVdrift = kVdriftPHDet;

   if(fMethodSecond) detVdrift = kVdriftLinear;

   AliCDBMetaData *metaData= new AliCDBMetaData();
   metaData->SetObjectClassName("AliTRDCalDet");
   metaData->SetResponsible("Raphaelle Bailhache");
   metaData->AddDateToComment();
   metaData->SetBeamPeriod(1);

   AliCDBId id1("TRD/Calib/ChamberVdrift", startRunNumber, endRunNumber);
   AliTRDCalDet *calDet = (AliTRDCalDet *) fCalibObjects->At(detVdrift);
   if(calDet) status = storage->Put(calDet, id1, metaData);
   if (status==kFALSE) fStatusPos = fStatusPos | kCalibFailedExport;

   //

   if(!fMethodSecond) {

     AliCDBMetaData *metaDataPad= new AliCDBMetaData();
     metaDataPad->SetObjectClassName("AliTRDCalPad");
     metaDataPad->SetResponsible("Raphaelle Bailhache");
     metaDataPad->AddDateToComment();
     metaDataPad->SetBeamPeriod(1);

     AliCDBId id1Pad("TRD/Calib/LocalVdrift", startRunNumber, endRunNumber);
     AliTRDCalPad *calPad = (AliTRDCalPad *) fCalibObjects->At(kVdriftPHPad);
     if(calPad) status = storage->Put(calPad, id1Pad, metaDataPad);
     if (status==kFALSE) fStatusPos = fStatusPos | kCalibFailedExport;

   }

 }
 //________________________________________________________________________________________________________________________
 void AliTRDPreprocessorOffline::UpdateOCDBT0(Int_t startRunNumber, Int_t endRunNumber, AliCDBStorage* storage){
   //
   // Update OCDB entry
   //

   Bool_t status = kTRUE;

   AliCDBMetaData *metaData= new AliCDBMetaData();
   metaData->SetObjectClassName("AliTRDCalDet");
   metaData->SetResponsible("Raphaelle Bailhache");
   metaData->AddDateToComment();
   metaData->SetBeamPeriod(1);

   AliCDBId id1("TRD/Calib/ChamberT0", startRunNumber, endRunNumber);
   AliTRDCalDet *calDet = (AliTRDCalDet *) fCalibObjects->At(kT0PHDet);
   if(calDet) storage->Put(calDet, id1, metaData);

   //

   AliCDBMetaData *metaDataPad= new AliCDBMetaData();
   metaDataPad->SetObjectClassName("AliTRDCalPad");
   metaDataPad->SetResponsible("Raphaelle Bailhache");
   metaDataPad->AddDateToComment();
   metaDataPad->SetBeamPeriod(1);

   AliCDBId id1Pad("TRD/Calib/LocalT0", startRunNumber, endRunNumber);
   AliTRDCalPad *calPad = (AliTRDCalPad *) fCalibObjects->At(kT0PHPad);
   if(calPad) status = storage->Put(calPad, id1Pad, metaDataPad);
   if (status==kFALSE) fStatusPos = fStatusPos | kCalibFailedExport;


 }
 //_________________________________________________________________________________________________________________
 void AliTRDPreprocessorOffline::UpdateOCDBPRF(Int_t startRunNumber, Int_t endRunNumber, AliCDBStorage* storage){
   //
   // Update OCDB entry
   //

   Bool_t status = kTRUE;

   AliCDBMetaData *metaData= new AliCDBMetaData();
   metaData->SetObjectClassName("AliTRDCalPad");
   metaData->SetResponsible("Raphaelle Bailhache");
   metaData->AddDateToComment();
   metaData->SetBeamPeriod(1);


   AliCDBId id1("TRD/Calib/PRFWidth", startRunNumber, endRunNumber);
   AliTRDCalPad *calPad = (AliTRDCalPad *) fCalibObjects->At(kPRF);
   if(calPad) status = storage->Put(calPad, id1, metaData);
   if (status==kFALSE) fStatusPos = fStatusPos | kCalibFailedExport;


 }
//_________________________________________________________________________________________________________________
void AliTRDPreprocessorOffline::UpdateOCDBPHQ(Int_t startRunNumber, Int_t endRunNumber, AliCDBStorage* storage)
{
  //
  // Update OCDB entry
  //
  Bool_t status = kTRUE;
  
  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Raphaelle Bailhache and Xianguo Lu");
  metaData->AddDateToComment();
  metaData->SetBeamPeriod(1);

  AliCDBId id1("TRD/Calib/PHQ", startRunNumber, endRunNumber);
  TObjArray *cobj = (TObjArray *) fCalibObjects->At(kPHQ);
  if(cobj){
    //cobj->Print();
     status = storage->Put(cobj, id1, metaData);
    if (status==kFALSE) fStatusPos = fStatusPos | kCalibFailedExport;
  }
}

 //_________________________________________________________________________________________________________________
 void AliTRDPreprocessorOffline::UpdateOCDBChamberStatus(Int_t startRunNumber, Int_t endRunNumber, AliCDBStorage* storage){
   //
   // Update OCDB entry
   //

   Bool_t status = kTRUE;

   AliCDBMetaData *metaData= new AliCDBMetaData();
   metaData->SetObjectClassName("AliTRDCalChamberStatus");
   metaData->SetResponsible("Raphaelle Bailhache and Julian Book");
   metaData->AddDateToComment();
   metaData->SetBeamPeriod(1);

   AliCDBId id1("TRD/Calib/ChamberStatus", startRunNumber, endRunNumber);
   AliTRDCalChamberStatus *calChamberStatus = (AliTRDCalChamberStatus *) fCalibObjects->At(kChamberStatus);
   if(calChamberStatus) status = storage->Put(calChamberStatus, id1, metaData);
   if (status==kFALSE) fStatusPos = fStatusPos | kCalibFailedExport;


 }
 //__________________________________________________________________________________________________________________________
 Bool_t AliTRDPreprocessorOffline::ValidateGain() {
   //
   // Validate OCDB entry
   //

   AliTRDCalDet *calDet = (AliTRDCalDet *) fCalibObjects->At(kGain);
   if(calDet) {
     Double_t mean = calDet->GetMean();
     Double_t rms = calDet->GetRMSRobust();
     if((mean > 0.2) && (mean < 1.4) && (rms < 0.5)) return kTRUE;
     //if((mean > 0.2) && (mean < 1.4)) return kTRUE;
     else {
       fStatusPos = fStatusPos | kGainErrorRange;
       return kFALSE;
     }
   }
   else return kFALSE;
   


 }
 //__________________________________________________________________________________________________________________________
 Bool_t AliTRDPreprocessorOffline::ValidateVdrift(){
   //
   // Update OCDB entry
   //

   Int_t detVdrift = kVdriftPHDet;
   Bool_t ok = kTRUE;

   if(fMethodSecond) detVdrift = kVdriftLinear;

   AliTRDCalDet *calDet = (AliTRDCalDet *) fCalibObjects->At(detVdrift);
   if(calDet) {
     Double_t mean = calDet->GetMean();
     Double_t rms = calDet->GetRMSRobust();
     //printf("Vdrift::mean %f, rms %f\n",mean,rms);
     if(!((mean > 1.0) && (mean < 2.0) && (rms < 0.5))) {
       fStatusPos = fStatusPos | kVdriftErrorRange;
       ok = kFALSE;
     }
   }
   else return kFALSE; 

   if(!fMethodSecond) {
     AliTRDCalPad *calPad = (AliTRDCalPad *) fCalibObjects->At(kVdriftPHPad);
     if(calPad) {
       Double_t mean = calPad->GetMean();
       Double_t rms = calPad->GetRMS();
       //printf("Vdrift::meanpad %f, rmspad %f\n",mean,rms);
       if(!((mean > 0.9) && (mean < 1.1) && (rms < 0.6))) {
	 fStatusPos = fStatusPos | kVdriftErrorRange;
	 ok = kFALSE;
       }
     }
     else return kFALSE;
   }

   return ok;

 }
 //__________________________________________________________________________________________________________________________
 Bool_t AliTRDPreprocessorOffline::ValidateExB(){
   //
   // Update OCDB entry
   //

   AliTRDCalDet *calDet = (AliTRDCalDet *) fCalibObjects->At(kLorentzLinear);
   if(calDet) {
     Double_t mean = calDet->GetMean();
     Double_t rms = calDet->GetRMSRobust();
     //printf("Vdrift::mean %f, rms %f\n",mean,rms);
     if(!((mean > -1.0) && (mean < 1.0) && (rms < 0.5))) {
       fStatusNeg = fStatusNeg | kExBErrorRange;
       return kFALSE;
     }
     else return kTRUE;
   }
   else return kFALSE; 
   
 }
 //__________________________________________________________________________________________________________________________
 Bool_t AliTRDPreprocessorOffline::ValidateT0(){
   //
   // Update OCDB entry
   //

   AliTRDCalDet *calDet = (AliTRDCalDet *) fCalibObjects->At(kT0PHDet);
   AliTRDCalPad *calPad = (AliTRDCalPad *) fCalibObjects->At(kT0PHPad);
   if(calDet && calPad) {
     Double_t meandet = calDet->GetMean();
     Double_t rmsdet = calDet->GetRMSRobust();
     Double_t meanpad = calPad->GetMean();
     //Double_t rmspad = calPad->GetRMS();
     //printf("T0::meandet %f, rmsdet %f,meanpad %f\n",meandet,rmsdet,meanpad);
     if((meandet >   fMinTimeOffsetValidate) && (meandet < 5.0) && (rmsdet < 4.0) && (meanpad < 5.0) && (meanpad > -0.5)) return kTRUE;
     else {
       fStatusPos = fStatusPos | kTimeOffsetErrorRange;
       return kFALSE;
     }
   }
   else return kFALSE;

 }
 //__________________________________________________________________________________________________________________________
 Bool_t AliTRDPreprocessorOffline::ValidatePRF() const{
   //
   // Update OCDB entry
   //

   AliTRDCalPad *calPad = (AliTRDCalPad *) fCalibObjects->At(kPRF);
   if(calPad) {
     Double_t meanpad = calPad->GetMean();
     Double_t rmspad = calPad->GetRMS();
     //printf("PRF::meanpad %f, rmspad %f\n",meanpad,rmspad);
     if((meanpad < 1.0) && (rmspad < 0.8)) return kTRUE;
     else return kFALSE;
   }
   else return kFALSE;


 }
 //__________________________________________________________________________________________________________________________
Bool_t AliTRDPreprocessorOffline::ValidateChamberStatus(){
  //
  // Update OCDB entry
  //
  
  AliTRDCalChamberStatus *calChamberStatus = (AliTRDCalChamberStatus *) fCalibObjects->At(kChamberStatus);
  if(calChamberStatus) {

    Int_t detectornodataA    = 0;
    Int_t detectornodataB    = 0;
    Int_t detectornodata     = 0;
    Int_t detectorbadcalib   = 0;
    Int_t detectornotgood    = 0;
    Int_t detectornotcalib   = 0;
    
    for(Int_t sm=0; sm < 18; sm++) {
      //printf("%d chambers w/o data in sm %d\n",fNoData[sm],sm);
      //printf("%d bad calibrated chambers in sm %d\n",fBadCalib[sm],sm);
      if(fNoData[sm] != 30) detectornodata += fNoData[sm];
      if(fNoDataA[sm] != 30) detectornodataA += fNoDataA[sm];
      if(fNoDataB[sm] != 30) detectornodataB += fNoDataB[sm];
      if(fNotCalib[sm] != 30) detectornotcalib += fNotCalib[sm];
      if(fNotGood[sm] != 30) detectornotgood += fNotGood[sm];
      detectorbadcalib+=fBadCalib[sm];
    }
    //printf("Number of chambers w/o data %d\n",detectornodata);
    //printf("Number of chambers w/o data A %d\n",detectornodataA);
    //printf("Number of chambers w/o data B %d\n",detectornodataB);
    //printf("Number of chambers bad calibrated %d\n",detectorbadcalib);
    //printf("Number of chambers not good %d\n",detectornotgood);
    //printf("Number of chambers not calibrated %d\n",detectornotcalib);

    if((detectornodata > fNoDataValidate) ||
       (detectornodataA > fNoDataValidate) ||
       (detectornodataB > fNoDataValidate) ||
       (detectorbadcalib > fBadCalibValidate)){
      fStatusPos = fStatusPos | kChamberStatusErrorRange;
      return kFALSE;
    }
    if((detectornotcalib > fNoDataValidate) ||
       (detectornotgood > fNoDataValidate)){
      fStatusNeg = fStatusNeg | kChamberStatusTooFewGood;
      return kFALSE;
    }
    return kTRUE;
  }
  else return kFALSE;
  
}
//_____________________________________________________________________________
Int_t AliTRDPreprocessorOffline::GetVersion(TString name) const
{
  //
  // Get version from the title
  //
  
  // Some patterns
  const Char_t *version = "Ver";
  if(!strstr(name.Data(),version)) return -1;
  const Char_t *after = "Subver";  
  if(!strstr(name.Data(),after)) return -1;

  for(Int_t ver = 0; ver < 999999999; ver++) {

    TString vertry(version);
    vertry += ver;
    vertry += after;

    //printf("vertry %s and name %s\n",vertry.Data(),name.Data());

    if(strstr(name.Data(),vertry.Data())) return ver;
    
  }
  
  return -1;

}

//_____________________________________________________________________________
Int_t AliTRDPreprocessorOffline::GetSubVersion(TString name) const
{
  //
  // Get subversion from the title
  //
  
  // Some patterns
  const Char_t *subversion = "Subver";
  if(!strstr(name.Data(),subversion)) return -1;
  const Char_t *after = "FirstRun";
  if(!strstr(name.Data(),after)) {
    after = "Nz";
  }
  if(!strstr(name.Data(),after)) return -1;

  
  for(Int_t ver = 0; ver < 999999999; ver++) {
    
    TString vertry(subversion);
    vertry += ver;
    vertry += after;

    //printf("vertry %s and name %s\n",vertry.Data(),name.Data());

    if(strstr(name.Data(),vertry.Data())) return ver;
    
  }
  
  return -1;

}

//_____________________________________________________________________________
Int_t AliTRDPreprocessorOffline::GetFirstRun(TString name) const
{
  //
  // Get first run from the title
  //
  
  // Some patterns
  const Char_t *firstrun = "FirstRun";
  if(!strstr(name.Data(),firstrun)) return -1;
  const Char_t *after = "Nz";  
  if(!strstr(name.Data(),after)) return -1;
  
  
  for(Int_t ver = 0; ver < 999999999; ver++) {

    TString vertry(firstrun);
    vertry += ver;
    vertry += after;

    //printf("vertry %s and name %s\n",vertry.Data(),name.Data());

    if(strstr(name.Data(),vertry.Data())) return ver;
    
  }
  
  return -1;

}
//_____________________________________________________________________________
Bool_t AliTRDPreprocessorOffline::CheckStatus(Int_t status, Int_t bitMask) const
{
  //
  // Checks the status
  //

  return (status & bitMask) ? kTRUE : kFALSE;
  
}
//_____________________________________________________________________________
Int_t AliTRDPreprocessorOffline::GetStatus() const
{
  //
  // Checks the status
  // fStatusPos: errors
  // fStatusNeg: only info
  //

  if(fStatusPos > 0) return fStatusPos;
  else return (-TMath::Abs(fStatusNeg));
  
}
//_____________________________________________________________________________
void AliTRDPreprocessorOffline::PrintStatus() const
{
  //
  // Do Summary
  //

  AliInfo(Form("The error status is %d",fStatusPos));
  AliInfo(Form("IsGainErrorOld? %d",(Int_t)IsGainErrorOld()));
  AliInfo(Form("IsVdriftErrorOld? %d",(Int_t)IsVdriftErrorOld()));
  AliInfo(Form("IsGainErrorRange? %d",(Int_t)IsGainErrorRange()));
  AliInfo(Form("IsVdriftErrorRange? %d",(Int_t)IsVdriftErrorRange()));
  AliInfo(Form("IsTimeOffsetErrorRange? %d",(Int_t)IsTimeOffsetErrorRange()));
  AliInfo(Form("IsChamberStatusErrorRange? %d",(Int_t)IsChamberStatusErrorRange()));
  AliInfo(Form("IsCalibFailedExport? %d",(Int_t)IsCalibFailedExport()));

 
  AliInfo(Form("The info status is %d",fStatusNeg));
  AliInfo(Form("IsGainNotEnoughStatsButFill? %d",(Int_t)IsGainNotEnoughStatsButFill()));
  AliInfo(Form("IsVdriftNotEnoughStatsButFill? %d",(Int_t)IsVdriftNotEnoughStatsButFill()));
  AliInfo(Form("IsGainNotEnoughStatsNotFill? %d",(Int_t)IsGainNotEnoughStatsNotFill()));
  AliInfo(Form("IsVdriftNotEnoughStatsNotFill? %d",(Int_t)IsVdriftNotEnoughStatsNotFill()));
  AliInfo(Form("IsTimeOffsetNotEnoughStatsNotFill? %d",(Int_t)IsTimeOffsetNotEnoughStatsNotFill()));
  AliInfo(Form("IsExBErrorRange? %d",(Int_t)IsExBErrorRange()));
  AliInfo(Form("IsExBErrorOld? %d",(Int_t)IsExBErrorOld()));
  AliInfo(Form("IsChamberStatusTooFewGood? %d",(Int_t)IsChamberStatusTooFewGood()));
 
  
}
//___________________________________________________________________________________
void AliTRDPreprocessorOffline::SetCalDetVdrift(AliTRDCalDet *calDetVdriftUsed) 
{

  fCalDetVdriftUsed = calDetVdriftUsed;

  fCalDetExBUsed = new AliTRDCalDet("lorentz angle tan","lorentz angle tan (detector value)");
  for(Int_t k = 0; k < 540; k++){
    fCalDetExBUsed->SetValue(k,AliTRDCommonParam::Instance()->GetOmegaTau(fCalDetVdriftUsed->GetValue(k)));
    //printf("Set the exb object for detector %d, vdrift %f and exb %f\n",k,fCalDetVdriftUsed->GetValue(k),fCalDetExBUsed->GetValue(k));
  }
  
};
//___________________________________________________________________________________
Bool_t AliTRDPreprocessorOffline::SetCalDetGain(Int_t runNumber, Int_t version, Int_t subversion) 
{
  //
  // Set the fCalDetGainUsed
  //

  if(version < 1) {
    fStatusPos = fStatusPos | kGainErrorOld;
    return kFALSE;
  }

  AliCDBEntry *entry = AliCDBManager::Instance()->Get("TRD/Calib/ChamberGainFactor",runNumber, version, subversion);
  if(!entry) {
    AliError("Found no entry\n");
    fStatusPos = fStatusPos | kGainErrorOld;
    return kFALSE;
  }
  //const AliCDBId id = entry->GetId();
  //version = id.GetVersion();
  //subversion = id.GetSubVersion();
  //printf("Found version %d and subversion %d for vdrift\n",version,subversion);
  AliTRDCalDet* calDet = (AliTRDCalDet *)entry->GetObject();
  if(calDet) fCalDetGainUsed = calDet;
  else {
    fStatusPos = fStatusPos | kGainErrorOld;
    return kFALSE;
  }
  
  return kTRUE;

}
//___________________________________________________________________________________
Bool_t AliTRDPreprocessorOffline::SetCalDetVdriftExB(Int_t runNumber, Int_t versionv, Int_t subversionv, Int_t versionexb, Int_t subversionexb) 
{
  //
  // Set the fCalDetVdriftUsed and fCalDetExBUsed
  //

  if(versionv < 1) {
    fStatusPos = fStatusPos | kVdriftErrorOld;
    fStatusPos = fStatusPos | kExBErrorOld;  
    return kFALSE;
  }

  AliCDBEntry *entry = AliCDBManager::Instance()->Get("TRD/Calib/ChamberVdrift",runNumber, versionv, subversionv);
  if(!entry) {
    AliError("Found no entry\n");
    fStatusPos = fStatusPos | kVdriftErrorOld;
    return kFALSE;
  }
  AliTRDCalDet* calDet = (AliTRDCalDet *)entry->GetObject();
  if(calDet) fCalDetVdriftUsed = calDet;
  else {
    fStatusPos = fStatusPos | kVdriftErrorOld;
    return kFALSE;
  }

  // ExB object

  if((versionexb == 0) && (subversionexb == 0)) {
    
    fCalDetExBUsed = new AliTRDCalDet("lorentz angle tan","lorentz angle tan (detector value)");
    for(Int_t k = 0; k < 540; k++){
      fCalDetExBUsed->SetValue(k,AliTRDCommonParam::Instance()->GetOmegaTau(fCalDetVdriftUsed->GetValue(k)));
      //printf("Nothing found: set the exb object for detector %d, vdrift %f and exb %f\n",k,fCalDetVdriftUsed->GetValue(k),fCalDetExBUsed->GetValue(k));
    }
  }
  else {

    entry = 0x0;
    entry = AliCDBManager::Instance()->Get("TRD/Calib/ChamberExB",runNumber, versionexb, subversionexb);
    if(!entry) {
      //printf("Found no entry\n");
      fStatusPos = fStatusPos | kExBErrorOld;   
      return kFALSE;
    }
    AliTRDCalDet* calDetexb = (AliTRDCalDet *)entry->GetObject();
    if(!calDetexb) {
      fStatusPos = fStatusPos | kExBErrorOld;   
      return kFALSE;
    }
    
    Double_t meanexb = calDetexb->GetMean();
    //printf("Mean value %f\n",meanexb);
    if((meanexb > 70) || (fNoExBUsedInReco)) {
      fCalDetExBUsed = new AliTRDCalDet("lorentz angle tan","lorentz angle tan (detector value)");
      for(Int_t k = 0; k < 540; k++){
	fCalDetExBUsed->SetValue(k,AliTRDCommonParam::Instance()->GetOmegaTau(fCalDetVdriftUsed->GetValue(k)));
	//printf("Found but: set the exb object for detector %d, vdrift %f and exb %f\n",k,fCalDetVdriftUsed->GetValue(k),fCalDetExBUsed->GetValue(k));
      }
    }
    else {
      fCalDetExBUsed = calDetexb;
    }
    
  }

  
  return kTRUE;

}

