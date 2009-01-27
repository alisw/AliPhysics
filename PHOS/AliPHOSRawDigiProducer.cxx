/**************************************************************************
 * Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved.      *
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

//This class produces PHOS digits of one event
//using AliPHOSRawDecoder. 
//
//   For example:
//   TClonesArray *digits = new TClonesArray("AliPHOSDigit",100);
//   AliRawReader* rawReader = new AliRawReaderDate("2006run2211.raw");
//   AliPHOSRawDecoder dc(rawReader);
//   while (rawReader->NextEvent()) {
//     AliPHOSRawDigiProducer producer;
//     producer.MakeDigits(digits,&dc);
//   }

// Author: Boris Polichtchouk

// --- ROOT system ---
#include "TClonesArray.h"

// --- AliRoot header files ---
#include "AliPHOSRawDigiProducer.h"
#include "AliPHOSRawDecoder.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSDigit.h"
#include "AliPHOSRecoParam.h"
#include "AliPHOSCalibData.h"
#include "AliPHOSPulseGenerator.h"
#include "AliLog.h"

ClassImp(AliPHOSRawDigiProducer)

AliPHOSCalibData * AliPHOSRawDigiProducer::fgCalibData  = 0 ; 

//--------------------------------------------------------------------------------------
AliPHOSRawDigiProducer::AliPHOSRawDigiProducer():
  TObject(),
  fEmcMinE(0.),
  fCpvMinE(0.),
  fSampleQualityCut(1.),
  fGlobalAltroOffset(0),
  fEmcCrystals(0),
  fGeom(0),
  fPulseGenerator(0)
{
  // Default constructor
}
//--------------------------------------------------------------------------------------
AliPHOSRawDigiProducer::AliPHOSRawDigiProducer(const AliPHOSRecoParam* recoParam):
  TObject(),
  fEmcMinE(0.),
  fCpvMinE(0.),
  fSampleQualityCut(1.),
  fGlobalAltroOffset(0),
  fEmcCrystals(0),
  fGeom(0),
  fPulseGenerator(0)
{
  // Constructor takes paramerters from recoParam

  if(!recoParam) AliFatal("Reconstruction parameters are not set!");

  fEmcMinE = recoParam->GetEMCRawDigitThreshold();
  fCpvMinE = recoParam->GetCPVMinE();
  fSampleQualityCut = recoParam->GetEMCSampleQualityCut() ;
  fGlobalAltroOffset = recoParam->GetGlobalAltroOffset() ;

  fGeom=AliPHOSGeometry::GetInstance() ;
  if(!fGeom) fGeom = AliPHOSGeometry::GetInstance("IHEP");

  fEmcCrystals=fGeom->GetNCristalsInModule()*fGeom->GetNModules() ;

  fPulseGenerator = new AliPHOSPulseGenerator();

  GetCalibrationParameters() ; 
}
//--------------------------------------------------------------------------------------                       
AliPHOSRawDigiProducer::AliPHOSRawDigiProducer(const AliPHOSRawDigiProducer &dp):
  TObject(),
  fEmcMinE(0.),
  fCpvMinE(0.),
  fSampleQualityCut(1.),
  fGlobalAltroOffset(0),
  fEmcCrystals(0),
  fGeom(0),
  fPulseGenerator(0)
{                                                          
  // Copy constructor

  fEmcMinE = dp.fEmcMinE ;
  fCpvMinE = dp.fCpvMinE ;
  fEmcCrystals = dp.fEmcCrystals ;
  fPulseGenerator = new AliPHOSPulseGenerator();
  fGeom = dp.fGeom ;
}
//--------------------------------------------------------------------------------------
AliPHOSRawDigiProducer& AliPHOSRawDigiProducer::operator= (const AliPHOSRawDigiProducer &dp)
{
  // Assign operator

  if(&dp == this) return *this;

  fEmcMinE = dp.fEmcMinE ;
  fCpvMinE = dp.fCpvMinE ;
  fSampleQualityCut = dp.fSampleQualityCut ;
  fGlobalAltroOffset = dp.fGlobalAltroOffset ;
  fEmcCrystals = dp.fEmcCrystals ;
  fGeom = dp.fGeom ;
  if(fPulseGenerator) delete fPulseGenerator ;
  fPulseGenerator = new AliPHOSPulseGenerator();
  return  *this;
} 
//--------------------------------------------------------------------------------------                                                   
AliPHOSRawDigiProducer::~AliPHOSRawDigiProducer()
{
  // Desctructor
  if(fPulseGenerator) delete fPulseGenerator ;
  fPulseGenerator=0 ;
}
//--------------------------------------------------------------------------------------
void AliPHOSRawDigiProducer::MakeDigits(TClonesArray *digits, AliPHOSRawDecoder* decoder) 
{
  //Makes the job.
  //TClonesArray *digits and raw data decoder should be provided by calling function.

  digits->Clear();
 
  Int_t    iDigit   = 0 ;
  Int_t relId[4], absId =0;

  const Double_t baseLine=1. ; //Minimal energy of digit in ADC ch. 
  const Double_t highLowDiff=2.; //Maximal difference between High and Low channels in LG adc channels 

  //Temporary array for LowGain digits
  TClonesArray tmpLG("AliPHOSDigit",10000) ;
  Int_t ilgDigit=0 ;


  //Read current altro offcet from RCU
  decoder->SetAmpOffset(fGlobalAltroOffset) ;

  //Let decoder subtract pedestals in case of ZS
  decoder->SetCalibData(fgCalibData) ;
  
  while (decoder->NextDigit()) {

    Double_t energy=decoder->GetEnergy() ; 
    if(energy<=baseLine) //in ADC channels
      continue ;

    //remove digits with bad shape. Decoder should calculate quality so that 
    //in default case quality [0,1], while larger values of quality mean somehow 
    //corrupted samples, 999 means obviously corrupted sample.
    //It is difficult to fit samples with overflow (even setting cut on overflow values)
    //because too few points are left to fit. So we do not evaluate samples with overflow
    if(decoder->GetSampleQuality() > fSampleQualityCut && !(decoder->IsOverflow()))
       continue ;

    Bool_t lowGainFlag = decoder->IsLowGain();

    relId[0] = decoder->GetModule();
    relId[1] = 0;
    relId[2] = decoder->GetRow();
    relId[3] = decoder->GetColumn();
    fGeom->RelToAbsNumbering(relId, absId);

    Double_t time = decoder->GetTime() ;
    time = CalibrateT(time,relId,lowGainFlag) ;

    energy = CalibrateE(energy,relId,lowGainFlag) ;

    if(energy <= 0.) 
       continue;

    if(lowGainFlag){
      new(tmpLG[ilgDigit]) AliPHOSDigit(-1,absId,(Float_t)energy,(Float_t)time);
      ilgDigit++ ; 
    }
    else{ 
      if(decoder->IsOverflow()) //Keep this digit to replace it by Low Gain later.
                                //If there is no LogGain it wil be removed by cut on Min E
        new((*digits)[iDigit]) AliPHOSDigit(-1,absId,-1.f,(Float_t)time);
      else
        new((*digits)[iDigit]) AliPHOSDigit(-1,absId,(Float_t)energy,(Float_t)time);
      iDigit++;
    }
  }

  //Now scan created LG and HG digits and keep only those which appeared in both lists 
  //replace energy of HighGain digits only if there is overflow
  //negative energy (overflow)
  digits->Sort() ;
  tmpLG.Sort() ;
  Int_t iLG = 0;
  Int_t nLG1 = tmpLG.GetEntriesFast()-1 ;

  for(Int_t iDig=0 ; iDig < digits->GetEntriesFast() ; iDig++) { 
    AliPHOSDigit * digHG = dynamic_cast<AliPHOSDigit*>(digits->At(iDig)) ;
    if (!digHG) continue;
    AliPHOSDigit * digLG = dynamic_cast<AliPHOSDigit*>(tmpLG.At(iLG)) ;
    while(digLG && iLG<nLG1 && digHG->GetId()> digLG->GetId()){
      iLG++ ;
      digLG = dynamic_cast<AliPHOSDigit*>(tmpLG.At(iLG)) ;
    }
    absId=digHG->GetId() ;                                                                                                         
    fGeom->AbsToRelNumbering(absId,relId) ;                                                                                              
 
    if(digLG && digHG->GetId() == digLG->GetId()){ //we found pair
      if(digHG->GetEnergy()<0.){ //This is overflow in HG
        digHG->SetTime(digLG->GetTime()) ;
        digHG->SetEnergy(digLG->GetEnergy()) ;
      } 
      else{ //Make approximate comparison of HG and LG energies
        Double_t de = (digHG->GetEnergy()-digLG->GetEnergy()) ; 
        if(TMath::Abs(de)>CalibrateE(double(highLowDiff),relId,1)){ //too strong difference, remove digit
          digits->RemoveAt(iDig) ;
        }
      }
    }
    else{ //no pair - remove
      // temporary fix for dead LG channels
      if(relId[2]%2==1 && relId[3]%16==4) 
        continue ;
      if(digHG->GetEnergy()>CalibrateE(double(5.),relId,1)) //One can not always find LG with Amp<5 ADC ch.
        digits->RemoveAt(iDig) ;                                                                                                            
    }
  }

  CleanDigits(digits) ;

}
//____________________________________________________________________________
Double_t AliPHOSRawDigiProducer::CalibrateE(Double_t amp, Int_t* relId, Bool_t isLowGain)
{
  // Convert EMC LG amplitude to HG (multipli by ~16)
  // Calibration parameters are taken from calibration data base 
  if(fgCalibData){ 
    Int_t   module = relId[0];  
    Int_t   column = relId[3];
    Int_t   row    = relId[2];
    if(relId[1]==0) { // this is EMC 
      if(isLowGain){
        amp*= fgCalibData->GetHighLowRatioEmc(module,column,row);
      }
      return amp ;         
    }         
  }          
  return 0;        
}
//____________________________________________________________________________
Double_t AliPHOSRawDigiProducer::CalibrateT(Double_t time, Int_t * relId, Bool_t /* isLowGain */)
{
  //Calibrate time
  time*=fPulseGenerator->GetRawFormatTimeTrigger() ;
  if(fgCalibData){
    Int_t   module = relId[0];
    Int_t   column = relId[3];
    Int_t   row    = relId[2];
    if(relId[1]==0) { // this is EMC
      time += fgCalibData->GetTimeShiftEmc(module,column,row);                   
      return time ;             
    }
  }
 
  return -999.;
}
//____________________________________________________________________________
void AliPHOSRawDigiProducer::CleanDigits(TClonesArray * digits)
{
  // remove digits with amplitudes below threshold.
  // remove digits in bad channels
  // sort digits with icreasing AbsId
  
  //remove digits in bad map and below threshold
  Bool_t isBadMap = 0 ;
  if(fgCalibData->GetNumOfEmcBadChannels()){
    isBadMap=1 ;
  }
  
  for(Int_t i=0; i<digits->GetEntriesFast(); i++){
    AliPHOSDigit * digit = static_cast<AliPHOSDigit*>(digits->At(i)) ;
    if(!digit)
      continue  ;
    if ( (IsInEMC(digit) && digit->GetEnergy() < fEmcMinE) ||
	 (IsInCPV(digit) && digit->GetEnergy() < fCpvMinE) ){
      digits->RemoveAt(i) ;
      continue ;
    }
    if(isBadMap){ //check bad map now
      Int_t relid[4] ;
      fGeom->AbsToRelNumbering(digit->GetId(), relid) ; 
      if(fgCalibData->IsBadChannelEmc(relid[0],relid[3],relid[2])){
	digits->RemoveAt(i) ;
      }
    }
  }

  //Compress, sort and set indexes
  digits->Compress() ;
//  digits->Sort(); already sorted earlier
  for (Int_t i = 0 ; i < digits->GetEntriesFast() ; i++) { 
    AliPHOSDigit *digit = static_cast<AliPHOSDigit*>( digits->At(i) ) ; 
    digit->SetIndexInList(i) ;     
  }
}
//____________________________________________________________________________
Bool_t AliPHOSRawDigiProducer::IsInEMC(AliPHOSDigit * digit) const
{
  // Tells if (true) or not (false) the digit is in a PHOS-EMC module
  return digit->GetId() <= fEmcCrystals ;

}

//____________________________________________________________________________
Bool_t AliPHOSRawDigiProducer::IsInCPV(AliPHOSDigit * digit) const
{
  // Tells if (true) or not (false) the digit is in a PHOS-CPV module
  return digit->GetId() > fEmcCrystals ;
}
//____________________________________________________________________________
void AliPHOSRawDigiProducer::GetCalibrationParameters() 
{
  // Set calibration parameters:
  // if calibration database exists, they are read from database,
  // otherwise, reconstruction stops in the constructor of AliPHOSCalibData
  //
  // It is a user responsilibity to open CDB before reconstruction, for example: 
  // AliCDBStorage* storage = AliCDBManager::Instance()->GetStorage("local://CalibDB");

  if (!fgCalibData){
    fgCalibData = new AliPHOSCalibData(-1); //use AliCDBManager's run number
  }
  if (fgCalibData->GetCalibDataEmc() == 0)
    AliFatal("Calibration parameters for PHOS EMC not found. Stop reconstruction.\n");
  if (fgCalibData->GetCalibDataCpv() == 0)
    AliFatal("Calibration parameters for PHOS CPV not found. Stop reconstruction.\n");
}
