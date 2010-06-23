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
#include "Riostream.h"
#include <fstream>
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH2I.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TH1F.h"
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


ClassImp(AliTRDPreprocessorOffline)

AliTRDPreprocessorOffline::AliTRDPreprocessorOffline():
  TNamed("TPCPreprocessorOffline","TPCPreprocessorOffline"),
  fMethodSecond(kTRUE),
  fNameList("TRDCalib"),
  fCalDetGainUsed(0x0),
  fCalDetVdriftUsed(0x0),
  fCH2d(0x0),
  fPH2d(0x0),
  fPRF2d(0x0),
  fAliTRDCalibraVdriftLinearFit(0x0),
  fNEvents(0x0),
  fAbsoluteGain(0x0),
  fPlots(new TObjArray(8)),
  fCalibObjects(new TObjArray(8)),
  fVersionGainUsed(0),
  fSubVersionGainUsed(0),
  fVersionVdriftUsed(0), 
  fSubVersionVdriftUsed(0)
{
  //
  // default constructor
  //
}

AliTRDPreprocessorOffline::~AliTRDPreprocessorOffline() {
  //
  // Destructor
  //

  if(fCalDetGainUsed) delete fCalDetGainUsed;
  if(fCalDetVdriftUsed) delete fCalDetVdriftUsed;
  if(fCH2d) delete fCH2d;
  if(fPH2d) delete fPH2d;
  if(fPRF2d) delete fPRF2d;
  if(fAliTRDCalibraVdriftLinearFit) delete fAliTRDCalibraVdriftLinearFit;
  if(fNEvents) delete fNEvents;
  if(fAbsoluteGain) delete fAbsoluteGain;
  if(fPlots) delete fPlots;
  if(fCalibObjects) delete fCalibObjects;
  
}

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
  //
  // 2. extraction of the information
  //
  if(ReadVdriftT0Global(file)) AnalyzeVdriftT0();
  if(ReadVdriftLinearFitGlobal(file)) AnalyzeVdriftLinearFit();
  //
  // 3. Append QA plots
  //
  //MakeDefaultPlots(fVdriftArray,fVdriftArray);
  //
  //
  // 4. update of OCDB
  //
  //
  UpdateOCDBVdrift(startRunNumber,endRunNumber,ocdbStorage);
  UpdateOCDBT0(startRunNumber,endRunNumber,ocdbStorage);

}


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
  // 1. Initialization 
  if(!ReadGainGlobal(file)) return;
  //
  //
  // 2. extraction of the information
  //
  AnalyzeGain();
  if(fCalDetGainUsed) CorrectFromDetGainUsed();
  if(fCalDetVdriftUsed) CorrectFromDetVdriftUsed();
  //
  // 3. Append QA plots
  //
  //MakeDefaultPlots(fVdriftArray,fVdriftArray);
  //
  //
  // 4. update of OCDB
  //
  //
  UpdateOCDBGain(startRunNumber,endRunNumber,ocdbStorage);
  
}


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
  // 4. update of OCDB
  //
  //
  UpdateOCDBPRF(startRunNumber,endRunNumber,ocdbStorage);
  
}

Bool_t AliTRDPreprocessorOffline::Init(const Char_t* fileName){
  //
  // read the calibration used during the reconstruction
  // 

  if(ReadVdriftT0Global(fileName)) {
    
    TString nameph = fPH2d->GetTitle();
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
   
  return kTRUE;
  
}


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



Bool_t AliTRDPreprocessorOffline::AnalyzeGain(){
  //
  // Analyze gain - produce the calibration objects
  //

  AliTRDCalibraFit *calibra = AliTRDCalibraFit::Instance();
  calibra->SetMinEntries(800); // If there is less than 1000 entries in the histo: no fit
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
    if(!fCalDetGainUsed) {
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
  
  calibra->ResetVectorFit();
  
  return ok;
  
}

Bool_t AliTRDPreprocessorOffline::AnalyzeVdriftT0(){
  //
  // Analyze VdriftT0 - produce the calibration objects
  //

  AliTRDCalibraFit *calibra = AliTRDCalibraFit::Instance();
  calibra->SetMinEntries(800*20); // If there is less than 1000 entries in the histo: no fit
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
  calibra->ResetVectorFit();
 
  return ok;
  
}

Bool_t AliTRDPreprocessorOffline::AnalyzeVdriftLinearFit(){
  //
  // Analyze vdrift linear fit - produce the calibration objects
  //

  AliTRDCalibraFit *calibra = AliTRDCalibraFit::Instance();
  calibra->SetMinEntries(800); // If there is less than 1000 entries in the histo: no fit
  fAliTRDCalibraVdriftLinearFit->FillPEArray();
  calibra->AnalyseLinearFitters(fAliTRDCalibraVdriftLinearFit);

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
    // Put them in the array
    fCalibObjects->AddAt(calDetVdrift,kVdriftLinear);
    fCalibObjects->AddAt(calDetLorentz,kLorentzLinear);
    fPlots->AddAt(coefDriftLinear,kVdriftLinear);
    fPlots->AddAt(coefLorentzAngle,kLorentzLinear);
    //
    ok = kTRUE;
  }
  
  calibra->ResetVectorFit();
  
  return ok;
  
}


Bool_t AliTRDPreprocessorOffline::AnalyzePRF(){
  //
  // Analyze PRF - produce the calibration objects
  //

  AliTRDCalibraFit *calibra = AliTRDCalibraFit::Instance();
  calibra->SetMinEntries(600); // If there is less than 1000 entries in the histo: no fit
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
  
    if(TMath::Abs(gainout-1.0) > 0.000001) calDetGain->SetValue(det,gaininit*gainout);
    else calDetGain->SetValue(det,mean);
  }


}

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
  
  for(Int_t det = 0; det < 540; det++) {
    
    Float_t vdriftinit = fCalDetVdriftUsed->GetValue(det);
    Float_t vdriftout = calDetVdrift->GetValue(det);

    Float_t gain = calDetGain->GetValue(det);
    if(vdriftout > 0.0) gain = gain*vdriftinit/vdriftout;
    calDetGain->SetValue(det,gain);

    
  }
 
}

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
//_____________________________________________________________________________
Int_t AliTRDPreprocessorOffline::GetVersion(TString name) const
{
  //
  // Get version from the title
  //
  
  // Some patterns
  const Char_t *version = "Ver";
  if(!strstr(name.Data(),version)) return -1;
  
  for(Int_t ver = 0; ver < 999999999; ver++) {

    TString vertry(version);
    vertry += ver;
    vertry += "Subver";

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
  
  for(Int_t ver = 0; ver < 999999999; ver++) {
    
    TString vertry(subversion);
    vertry += ver;
    vertry += "Nz";

    //printf("vertry %s and name %s\n",vertry.Data(),name.Data());

    if(strstr(name.Data(),vertry.Data())) return ver;
    
  }
  
  return -1;

}

