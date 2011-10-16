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
#include "AliTRDPreprocessorOffline.h"
#include "AliTRDCalChamberStatus.h"
#include "AliTRDCommonParam.h"


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
  fAliTRDCalibraVdriftLinearFit(0x0),
  fNEvents(0x0),
  fAbsoluteGain(0x0),
  fPlots(new TObjArray(9)),
  fCalibObjects(new TObjArray(9)),
  fVersionGainUsed(0),
  fSubVersionGainUsed(0),
  fFirstRunVdriftUsed(0),
  fVersionVdriftUsed(0), 
  fSubVersionVdriftUsed(0),
  fSwitchOnValidation(kTRUE),
  fVdriftValidated(kFALSE),
  fExBValidated(kFALSE),
  fT0Validated(kFALSE),
  fMinStatsVdriftT0PH(800*20),
  fMinStatsVdriftLinear(800),
  fMinStatsGain(800),
  fMinStatsPRF(600),
  fBackCorrectGain(kFALSE),  
  fBackCorrectVdrift(kTRUE),
  fNotEnoughStatisticsForTheGain(kFALSE),
  fNotEnoughStatisticsForTheVdriftLinear(kFALSE),
  fStatusNeg(0),
  fStatusPos(0)
{
  //
  // default constructor
  //
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
  if(fAliTRDCalibraVdriftLinearFit) delete fAliTRDCalibraVdriftLinearFit;
  if(fNEvents) delete fNEvents;
  if(fAbsoluteGain) delete fAbsoluteGain;
  if(fPlots) delete fPlots;
  if(fCalibObjects) delete fCalibObjects;
  
}
//___________________________________________________________________________________________________________________

void AliTRDPreprocessorOffline::CalibVdriftT0(const Char_t* file, Int_t startRunNumber, Int_t endRunNumber, TString ocdbStorage){
  //
  // make calibration of the drift velocity
  // Input parameters:
  //      file                             - the location of input file
  //      startRunNumber, endRunNumber     - run validity period 
  //      ocdbStorage                      - path to the OCDB storage
  //                                       - if empty - local storage 'pwd' uesed
  if (ocdbStorage.Length()<=0) ocdbStorage="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
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
  if(ReadVdriftLinearFitGlobal(file) && fCalDetVdriftUsed && fCalDetExBUsed) AnalyzeVdriftLinearFit();
  if(ReadVdriftT0Global(file)) AnalyzeVdriftT0();
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
//_________________________________________________________________________________________________________________

void AliTRDPreprocessorOffline::CalibGain(const Char_t* file, Int_t startRunNumber, Int_t endRunNumber, TString ocdbStorage){
  //
  // make calibration of the drift velocity
  // Input parameters:
  //      file                             - the location of input file
  //      startRunNumber, endRunNumber     - run validity period 
  //      ocdbStorage                      - path to the OCDB storage
  //                                       - if empty - local storage 'pwd' uesed
  if (ocdbStorage.Length()<=0) ocdbStorage="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
  //
  fNotEnoughStatisticsForTheGain = kFALSE;
  //
  // 1. Initialization 
  if(!ReadGainGlobal(file)) return;
  //
  //
  // 2. extraction of the information
  //
  AnalyzeGain();
  if(fBackCorrectGain) CorrectFromDetGainUsed();
  if(fBackCorrectVdrift) CorrectFromDetVdriftUsed();
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

void AliTRDPreprocessorOffline::CalibPRF(const Char_t* file, Int_t startRunNumber, Int_t endRunNumber, TString ocdbStorage){
  //
  // make calibration of the drift velocity
  // Input parameters:
  //      file                             - the location of input file
  //      startRunNumber, endRunNumber     - run validity period 
  //      ocdbStorage                      - path to the OCDB storage
  //                                       - if empty - local storage 'pwd' uesed
  if (ocdbStorage.Length()<=0) ocdbStorage="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
  //
  // 1. Initialization 
  if(!ReadPRFGlobal(file)) return;
  //
  //
  // 2. extraction of the information
  //
  AnalyzePRF();
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

void AliTRDPreprocessorOffline::CalibChamberStatus(Int_t startRunNumber, Int_t endRunNumber, TString ocdbStorage){
  //
  // make calibration of the chamber status
  // Input parameters:
  //      startRunNumber, endRunNumber     - run validity period 
  //      ocdbStorage                      - path to the OCDB storage
  //                                       - if empty - local storage 'pwd' uesed
  if (ocdbStorage.Length()<=0) ocdbStorage="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
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
  if(fSwitchOnValidation==kTRUE && ValidateChamberStatus()==kFALSE) { 
    //AliError("TRD Chamber status OCDB parameters not ok!");
    return;
  }
  //
  // 5. update of OCDB
  //
  //
  if((!fNotEnoughStatisticsForTheGain) && (!fNotEnoughStatisticsForTheGain)) UpdateOCDBChamberStatus(startRunNumber,endRunNumber,ocdbStorage);
  
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
    fVersionGainUsed = GetVersion(namech);  
    fSubVersionGainUsed = GetSubVersion(namech);    

    //printf("Found Version %d, Subversion %d for gain\n",fVersionGainUsed,fSubVersionGainUsed);

  }
   
  if(fVersionVdriftUsed == 0) fStatusPos = fStatusPos |kVdriftErrorOld;
  if(fVersionGainUsed == 0) fStatusPos = fStatusPos | kGainErrorOld;
 
  return kTRUE;
  
}
//___________________________________________________________________________________________________________________

Bool_t AliTRDPreprocessorOffline::ReadGainGlobal(const Char_t* fileName){
  //
  // read calibration entries from file
  // 
  if(fCH2d) return kTRUE;
  TFile fcalib(fileName);
  TObjArray * array = (TObjArray*)fcalib.Get(fNameList);
  if (array){
    TH2I *ch2d = (TH2I *) array->FindObject("CH2d");
    if(!ch2d) return kFALSE;
    fCH2d = (TH2I*)ch2d->Clone();
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
  TObjArray * array = (TObjArray*)fcalib.Get(fNameList);
  if (array){
    TProfile2D *ph2d = (TProfile2D *) array->FindObject("PH2d");
    if(!ph2d) return kFALSE;
    fPH2d = (TProfile2D*)ph2d->Clone();
    //fNEvents = (TH1I *) array->FindObject("NEvents");
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
  TObjArray * array = (TObjArray*)fcalib.Get(fNameList);
  if (array){
    fAliTRDCalibraVdriftLinearFit = (AliTRDCalibraVdriftLinearFit *) array->FindObject("AliTRDCalibraVdriftLinearFit");
    //fNEvents = (TH1I *) array->FindObject("NEvents");
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

Bool_t AliTRDPreprocessorOffline::ReadPRFGlobal(const Char_t* fileName){
  //
  // read calibration entries from file
  // 
  if(fPRF2d) return kTRUE;
  TFile fcalib(fileName);
  TObjArray * array = (TObjArray*)fcalib.Get(fNameList);
  if (array){
    TProfile2D *prf2d = (TProfile2D *) array->FindObject("PRF2d");
    if(!prf2d) return kFALSE;
    fPRF2d = (TProfile2D*)prf2d->Clone();
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
  calibra->SetMinEntries(fMinStatsGain); // If there is less than 1000 entries in the histo: no fit
  calibra->AnalyseCH(fCH2d);

  Int_t nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(0))
    + 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(0));
  Int_t nbfit       = calibra->GetNumberFit();
  Int_t nbE         = calibra->GetNumberEnt();


  Bool_t ok = kFALSE;
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
    ok = kTRUE;
  }
  else {
    fNotEnoughStatisticsForTheGain = kTRUE;
    Int_t minStatsGain = fMinStatsGain*30;
    calibra->SetMinEntries(minStatsGain); // Because we do it for all, we increase this
    Double_t gainoverallnotnormalized =  calibra->AnalyseCHAllTogether(fCH2d);
    if(fCalDetGainUsed && (gainoverallnotnormalized > 0.0)) {
      AliTRDCalDet *calDetGain = new AliTRDCalDet(*fCalDetGainUsed);
      Double_t oldmean = fCalDetGainUsed->CalcMean(kFALSE);
      //printf("oldmean %f\n",oldmean);
      if(oldmean > 0.0)  {
      	Double_t scalefactor = calibra->GetScaleFactorGain();
	//printf("Correction factor %f\n",gainoverallnotnormalized*scalefactor);
	calDetGain->Multiply(gainoverallnotnormalized*scalefactor/oldmean);
	//printf("newmean %f\n",calDetGain->CalcMean(kFALSE));
	TH1F *coefGain  = calDetGain->MakeHisto1DAsFunctionOfDet();
	fCalibObjects->AddAt(calDetGain,kGain);
	fPlots->AddAt(coefGain,kGain);
	// 
	ok = kTRUE;
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
  
  return ok;
  
}
//_____________________________________________________________________________________________________
Bool_t AliTRDPreprocessorOffline::AnalyzeVdriftT0(){
  //
  // Analyze VdriftT0 - produce the calibration objects
  //

  AliTRDCalibraFit *calibra = AliTRDCalibraFit::Instance();
  calibra->SetMinEntries(fMinStatsVdriftT0PH); // If there is less than 1000 entries in the histo: no fit
  calibra->AnalysePH(fPH2d);

  Int_t nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(1))
    + 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(1));
  Int_t nbfit       = calibra->GetNumberFit();
  Int_t nbfitSuccess = calibra->GetNumberFitSuccess();
  Int_t nbE         = calibra->GetNumberEnt();

  //printf("nbtg %d, nbfit %d, nbE %d, nbfitSuccess %d\n",nbtg,nbfit,nbE,nbfitSuccess);

  Bool_t ok = kFALSE;
  if ((nbtg >                  0) && 
      (nbfit        >= 0.5*nbE) && (nbE > 30) && (nbfitSuccess > 30)) {
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
    ok = kTRUE;
  }
  else {
    //printf("Not enough stats timeoffset\n");
    fStatusNeg = fStatusNeg | kTimeOffsetNotEnoughStatsNotFill;
  }
  calibra->ResetVectorFit();
 
  return ok;
  
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
  //printf("Fill PE Array\n");
  fAliTRDCalibraVdriftLinearFit->FillPEArray();
  //printf("AliTRDCalibraFit\n");
  calibra->AnalyseLinearFitters(fAliTRDCalibraVdriftLinearFit);
  //printf("After\n");

  //Int_t nbtg        = 540;
  Int_t nbfit       = calibra->GetNumberFit();
  Int_t nbE         = calibra->GetNumberEnt();

  
  Bool_t ok = kFALSE;
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
    ok = kTRUE;
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
      Double_t oldmeanvdrift = fCalDetVdriftUsed->CalcMean(kFALSE);
      Double_t oldmeanexb = fCalDetExBUsed->CalcMean(kFALSE);
      //printf("oldmean %f\n",oldmean);
      if((oldmeanvdrift > 0.0) && (oldmeanexb < 70.0))  {
	//printf("Correction factor %f\n",vdriftoverall);
	calDetVdrift->Multiply(vdriftoverall/oldmeanvdrift);
	calDetLorentz->Multiply(exboverall/oldmeanexb);
	//printf("newmean %f\n",calDetVdrift->CalcMean(kFALSE));
	TH1F *coefDriftLinear  = calDetVdrift->MakeHisto1DAsFunctionOfDet();
	TH1F *coefLorentzAngle = calDetLorentz->MakeHisto1DAsFunctionOfDet();
	// Put them in the array
	fCalibObjects->AddAt(calDetVdrift,kVdriftLinear);
	fCalibObjects->AddAt(calDetLorentz,kLorentzLinear);
	fPlots->AddAt(coefDriftLinear,kVdriftLinear);
	fPlots->AddAt(coefLorentzAngle,kLorentzLinear);
	// 
	ok = kTRUE;
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
  
  return ok;
  
}
//________________________________________________________________________________________________________________

Bool_t AliTRDPreprocessorOffline::AnalyzePRF(){
  //
  // Analyze PRF - produce the calibration objects
  //

  AliTRDCalibraFit *calibra = AliTRDCalibraFit::Instance();
  calibra->SetMinEntries(fMinStatsPRF); // If there is less than 1000 entries in the histo: no fit
  calibra->AnalysePRFMarianFit(fPRF2d);

  Int_t nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(2))
    + 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(2));
  Int_t nbfit       = calibra->GetNumberFit();
  Int_t nbE         = calibra->GetNumberEnt();

  
  Bool_t ok = kFALSE;
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
    ok = kTRUE;
  }
  
  calibra->ResetVectorFit();
  
  return ok;
  
}

//_____________________________________________________________________________
Bool_t AliTRDPreprocessorOffline::AnalyzeChamberStatus()
{
  //
  // Produce AliTRDCalChamberStatus out of calibration results
  //
  
  // set up AliTRDCalChamberStatus
  AliTRDCalChamberStatus *CalChamberStatus = new AliTRDCalChamberStatus();
  for(Int_t det = 0; det < 540; det++) CalChamberStatus->SetStatus(det,1);

  // get calibration objects
  AliTRDCalDet *calDetGain   = (AliTRDCalDet *) fCalibObjects->At(kGain);
  AliTRDCalDet *calDetVDrift = (AliTRDCalDet *) fCalibObjects->At(kVdriftLinear);
  AliTRDCalDet *calDetExB    = (AliTRDCalDet *) fCalibObjects->At(kLorentzLinear);

  // Check
  if((!calDetGain) || (!calDetVDrift) || (!fCH2d) || (!calDetExB)) return kFALSE;

  // Gain
  Double_t gainmean = calDetGain->GetMean();
  Double_t vdriftmean = calDetVDrift->GetMean();
  Double_t exbmean = calDetExB->GetMean();

  Double_t gainrms = calDetGain->GetRMSRobust();
  Double_t vdriftrms = calDetVDrift->GetRMSRobust();
  Double_t exbrms = calDetExB->GetRMSRobust();

  //printf("Gain mean: %f, rms: %f\n",gainmean,gainrms);
  //printf("Vdrift mean: %f, rms: %f\n",vdriftmean,vdriftrms);
  //printf("ExB mean: %f, rms: %f\n",exbmean,exbrms);

  // Check
  if((TMath::Abs(gainrms) < 0.001) || (TMath::Abs(vdriftrms) < 0.001) || (TMath::Abs(exbrms) < 0.0000001)) return kFALSE;

  // mask chambers with empty gain entries
  //Int_t counter = 0;
  for (Int_t idet = 0; idet < 540; idet++) {

    // ch2d
    TH1I *projch =  (TH1I *) fCH2d->ProjectionX("projch",idet+1,idet+1,(Option_t *)"e");
    Double_t entries = projch->GetEntries();

    // gain
    Double_t gain = calDetGain->GetValue(idet);

    // vdrift
    Double_t vdrift = calDetVDrift->GetValue(idet);

    // exb
    Double_t exb = calDetExB->GetValue(idet);


    if(entries<=0.5 ||
       TMath::Abs(gainmean-gain) > (15.0*gainrms) ||
       TMath::Abs(vdriftmean-vdrift) > (15.0*vdriftrms) ||
       TMath::Abs(exbmean-exb) > (50.0*exbrms)) {
     
      //printf(" chamber det %03d masked \n",idet);
      //printf(" gainmean %f and gain %f, gainrms %f \n",gainmean,gain,gainrms);
      //printf(" vdriftmean %f and vdrift %f, vdriftrms %f \n",vdriftmean,vdrift,vdriftrms);
      //printf(" exbmean %f and exb %f, exbrms %f \n",exbmean,exb,exbrms);
      
      CalChamberStatus->SetStatus(idet,AliTRDCalChamberStatus::kMasked);
      //counter++;
    }

     /*
     // installed supermodules+1 -> abort
     if(counter > (7+1)*30) {
       printf("ERROR: more than one SM to be masked!! \n Abort...\n");
       if(projch) delete projch;
       return 0x0;
     }
     */

    delete projch;
    
   }

   // Security
   for(Int_t sm=0; sm < 18; sm++) {
     Int_t counter = 0;
     for(Int_t det = 0; det < 30; det++){
       Int_t detector = sm*30+det;
       if(CalChamberStatus->IsMasked(detector)) counter++;
     }
     if(counter >= 10) {
       for(Int_t det = 0; det < 30; det++){
	 Int_t detector = sm*30+det;
	 CalChamberStatus->SetStatus(detector,AliTRDCalChamberStatus::kInstalled);
       }
     }
   }

   fCalibObjects->AddAt(CalChamberStatus,kChamberStatus);
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
 void AliTRDPreprocessorOffline::UpdateOCDBGain(Int_t startRunNumber, Int_t endRunNumber, const Char_t *storagePath){
   //
   // Update OCDB entry
   //

   AliCDBMetaData *metaData= new AliCDBMetaData();
   metaData->SetObjectClassName("AliTRDCalDet");
   metaData->SetResponsible("Raphaelle Bailhache");
   metaData->SetBeamPeriod(1);

   AliCDBId id1("TRD/Calib/ChamberGainFactor", startRunNumber, endRunNumber);
   AliCDBStorage * gStorage = AliCDBManager::Instance()->GetStorage(storagePath);
   AliTRDCalDet *calDet = (AliTRDCalDet *) fCalibObjects->At(kGain);
   if(calDet) gStorage->Put(calDet, id1, metaData);


 }
 //___________________________________________________________________________________________________________________
 void AliTRDPreprocessorOffline::UpdateOCDBExB(Int_t startRunNumber, Int_t endRunNumber, const Char_t *storagePath){
   //
   // Update OCDB entry
   //

   Int_t detExB = kLorentzLinear;
   if(!fMethodSecond) return;

   //printf("Pass\n");

   AliCDBMetaData *metaData= new AliCDBMetaData();
   metaData->SetObjectClassName("AliTRDCalDet");
   metaData->SetResponsible("Raphaelle Bailhache");
   metaData->SetBeamPeriod(1);

   AliCDBId id1("TRD/Calib/ChamberExB", startRunNumber, endRunNumber);
   AliCDBStorage * gStorage = AliCDBManager::Instance()->GetStorage(storagePath);
   AliTRDCalDet *calDet = (AliTRDCalDet *) fCalibObjects->At(detExB);
   if(calDet) gStorage->Put(calDet, id1, metaData);
   //if(!calDet) printf("No caldet\n");

 }
 //___________________________________________________________________________________________________________________
 void AliTRDPreprocessorOffline::UpdateOCDBVdrift(Int_t startRunNumber, Int_t endRunNumber, const Char_t *storagePath){
   //
   // Update OCDB entry
   //

   Int_t detVdrift = kVdriftPHDet;

   if(fMethodSecond) detVdrift = kVdriftLinear;

   AliCDBMetaData *metaData= new AliCDBMetaData();
   metaData->SetObjectClassName("AliTRDCalDet");
   metaData->SetResponsible("Raphaelle Bailhache");
   metaData->SetBeamPeriod(1);

   AliCDBId id1("TRD/Calib/ChamberVdrift", startRunNumber, endRunNumber);
   AliCDBStorage * gStorage = AliCDBManager::Instance()->GetStorage(storagePath);
   AliTRDCalDet *calDet = (AliTRDCalDet *) fCalibObjects->At(detVdrift);
   if(calDet) gStorage->Put(calDet, id1, metaData);

   //

   if(!fMethodSecond) {

     AliCDBMetaData *metaDataPad= new AliCDBMetaData();
     metaDataPad->SetObjectClassName("AliTRDCalPad");
     metaDataPad->SetResponsible("Raphaelle Bailhache");
     metaDataPad->SetBeamPeriod(1);

     AliCDBId id1Pad("TRD/Calib/LocalVdrift", startRunNumber, endRunNumber);
     AliTRDCalPad *calPad = (AliTRDCalPad *) fCalibObjects->At(kVdriftPHPad);
     if(calPad) gStorage->Put(calPad, id1Pad, metaDataPad);

   }

 }
 //________________________________________________________________________________________________________________________
 void AliTRDPreprocessorOffline::UpdateOCDBT0(Int_t startRunNumber, Int_t endRunNumber, const Char_t *storagePath){
   //
   // Update OCDB entry
   //

   AliCDBMetaData *metaData= new AliCDBMetaData();
   metaData->SetObjectClassName("AliTRDCalDet");
   metaData->SetResponsible("Raphaelle Bailhache");
   metaData->SetBeamPeriod(1);

   AliCDBId id1("TRD/Calib/ChamberT0", startRunNumber, endRunNumber);
   AliCDBStorage * gStorage = AliCDBManager::Instance()->GetStorage(storagePath);
   AliTRDCalDet *calDet = (AliTRDCalDet *) fCalibObjects->At(kT0PHDet);
   if(calDet) gStorage->Put(calDet, id1, metaData);

   //

   AliCDBMetaData *metaDataPad= new AliCDBMetaData();
   metaDataPad->SetObjectClassName("AliTRDCalPad");
   metaDataPad->SetResponsible("Raphaelle Bailhache");
   metaDataPad->SetBeamPeriod(1);

   AliCDBId id1Pad("TRD/Calib/LocalT0", startRunNumber, endRunNumber);
   AliTRDCalPad *calPad = (AliTRDCalPad *) fCalibObjects->At(kT0PHPad);
   if(calPad) gStorage->Put(calPad, id1Pad, metaDataPad);



 }
 //_________________________________________________________________________________________________________________
 void AliTRDPreprocessorOffline::UpdateOCDBPRF(Int_t startRunNumber, Int_t endRunNumber, const Char_t *storagePath){
   //
   // Update OCDB entry
   //

   AliCDBMetaData *metaData= new AliCDBMetaData();
   metaData->SetObjectClassName("AliTRDCalPad");
   metaData->SetResponsible("Raphaelle Bailhache");
   metaData->SetBeamPeriod(1);

   AliCDBId id1("TRD/Calib/PRFWidth", startRunNumber, endRunNumber);
   AliCDBStorage * gStorage = AliCDBManager::Instance()->GetStorage(storagePath);
   AliTRDCalPad *calPad = (AliTRDCalPad *) fCalibObjects->At(kPRF);
   if(calPad) gStorage->Put(calPad, id1, metaData);


 }
 //_________________________________________________________________________________________________________________
 void AliTRDPreprocessorOffline::UpdateOCDBChamberStatus(Int_t startRunNumber, Int_t endRunNumber, const Char_t *storagePath){
   //
   // Update OCDB entry
   //

   AliCDBMetaData *metaData= new AliCDBMetaData();
   metaData->SetObjectClassName("AliTRDCalChamberStatus");
   metaData->SetResponsible("Raphaelle Bailhache");
   metaData->SetBeamPeriod(1);

   AliCDBId id1("TRD/Calib/ChamberStatus", startRunNumber, endRunNumber);
   AliCDBStorage * gStorage = AliCDBManager::Instance()->GetStorage(storagePath);
   AliTRDCalChamberStatus *calChamberStatus = (AliTRDCalChamberStatus *) fCalibObjects->At(kChamberStatus);
   if(calChamberStatus) gStorage->Put(calChamberStatus, id1, metaData);


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
     //printf("T0::minimum %f, rmsdet %f,meanpad %f, rmspad %f\n",meandet,rmsdet,meanpad,rmspad);
     if((meandet > -1.5) && (meandet < 5.0) && (rmsdet < 4.0) && (meanpad < 5.0) && (meanpad > -0.5)) return kTRUE;
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
    Int_t detectormasked = 0;
    for(Int_t det = 0; det < 540; det++) {
      if(calChamberStatus->IsMasked(det)) detectormasked++;
    }
    //printf("Number of chambers masked %d\n",detectormasked);
    if(detectormasked > 40) {
      fStatusPos = fStatusPos | kChamberStatusErrorRange;
      return kFALSE;
    }
    else return kTRUE;
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

 
  AliInfo(Form("The info status is %d",fStatusNeg));
  AliInfo(Form("IsGainNotEnoughStatsButFill? %d",(Int_t)IsGainNotEnoughStatsButFill()));
  AliInfo(Form("IsVdriftNotEnoughStatsButFill? %d",(Int_t)IsVdriftNotEnoughStatsButFill()));
  AliInfo(Form("IsGainNotEnoughStatsNotFill? %d",(Int_t)IsGainNotEnoughStatsNotFill()));
  AliInfo(Form("IsVdriftNotEnoughStatsNotFill? %d",(Int_t)IsVdriftNotEnoughStatsNotFill()));
  AliInfo(Form("IsTimeOffsetNotEnoughStatsNotFill? %d",(Int_t)IsTimeOffsetNotEnoughStatsNotFill()));

  AliInfo(Form("IsExBErrorRange? %d",(Int_t)IsExBErrorRange()));
  AliInfo(Form("IsExBErrorOld? %d",(Int_t)IsExBErrorOld()));
  
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


