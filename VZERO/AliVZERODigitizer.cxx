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
 
///_________________________________________________________________________
///
/// This class constructs Digits out of Hits
///
///

// --- Standard library ---

// --- ROOT system ---
#include <TTree.h>
#include <TRandom.h>

// --- AliRoot header files ---
#include "AliVZEROConst.h"
#include "AliRun.h"
#include "AliVZERO.h"
#include "AliVZEROhit.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliRunDigitizer.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliVZEROCalibData.h"

#include "AliVZEROdigit.h"
#include "AliVZERODigitizer.h"

ClassImp(AliVZERODigitizer)

 AliVZERODigitizer::AliVZERODigitizer()
{
  // default constructor

   fNdigits = 0;
   fDigits  = 0;
  
   fPhotoCathodeEfficiency =   0.18;
   fPMVoltage              =  768.0;
   fPMGain = TMath::Power((fPMVoltage / 112.5) ,7.04277); 
   
   fCalibData = GetCalibData();
}

//____________________________________________________________________________ 
  AliVZERODigitizer::AliVZERODigitizer(AliRunDigitizer* manager)
                    :AliDigitizer(manager)
                    
{
  // constructor
  
  fNdigits = 0;
  fDigits  = 0;
  
  fPhotoCathodeEfficiency =   0.18;
  fPMVoltage              =  768.0;
  fPMGain = TMath::Power( (fPMVoltage / 112.5) ,7.04277 );
  
  fCalibData = GetCalibData();
  
}
           
//____________________________________________________________________________ 
  AliVZERODigitizer::~AliVZERODigitizer()
{
  // destructor
  
  if (fDigits) {
    fDigits->Delete();
    delete fDigits;
    fDigits=0; 
  }
}

//_____________________________________________________________________________
Bool_t AliVZERODigitizer::Init()
{
  // Initialises the digitizer

  // Initialises the Digit array
  fDigits = new TClonesArray ("AliVZEROdigit", 1000);

  return kTRUE;
}

//____________________________________________________________________________
void AliVZERODigitizer::Exec(Option_t* /*option*/) 
{   
  // Creates digits from hits
     
  Int_t       map[80];    // 48 values on V0C + 32 on V0A
  Int_t       adc[64];    // 32 PMs on V0C + 32 PMs on V0A
  Float_t    time[80], time2[64], adc_gain[80];    
  fNdigits     =    0;  
  Float_t cPM  = fPhotoCathodeEfficiency * fPMGain;

  // Retrieval of ADC gain values from local CDB 
  // I use only the first 64th values of the calibration array in CDB 
  // as I have no beam burst structure - odd or even beam burst number
  
  // Reminder : We have 16 scintillating cells mounted on 8 PMs 
  // on Ring 3 and Ring 4 in V0C -  added to produce  ADC outputs 
  // on these rings... 
   
  for(Int_t i=0; i<16; i++) { adc_gain[i]  = fCalibData->GetGain(i) ; };
  
  for(Int_t j=16; j<48; j=j+2) { 
            Int_t i=(j+17)/2;
            adc_gain[j]   = fCalibData->GetGain(i) ; 
	    adc_gain[j+1] = fCalibData->GetGain(i) ; }
  for(Int_t i=48; i<80; i++) { adc_gain[i]  = fCalibData->GetGain(i-16) ; }
  
//  for(Int_t i=0; i<80; i++) { printf(" i = %d gain = %f\n\n", i, adc_gain[i] );} 
            
  AliRunLoader* outRunLoader = 
    AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());    
  if (!outRunLoader) {
    Error("Exec", "Can not get output Run Loader");
    return;}
    
  AliLoader* outLoader = outRunLoader->GetLoader("VZEROLoader");
  if (!outLoader) {
    Error("Exec", "Can not get output VZERO Loader");
    return;}

  outLoader->LoadDigits("update");
  if (!outLoader->TreeD()) outLoader->MakeTree("D");
  outLoader->MakeDigitsContainer();
  TTree* treeD  = outLoader->TreeD();
  Int_t bufsize = 16000;
  treeD->Branch("VZERODigit", &fDigits, bufsize); 

  for (Int_t iInput = 0; iInput < fManager->GetNinputs(); iInput++) {
    AliRunLoader* runLoader = 
      AliRunLoader::GetRunLoader(fManager->GetInputFolderName(iInput));
    AliLoader* loader = runLoader->GetLoader("VZEROLoader");
    if (!loader) {
      Error("Exec", "Can not get VZERO Loader for input %d", iInput);
      continue;}
      
    if (!runLoader->GetAliRun()) runLoader->LoadgAlice();

    AliVZERO* vzero = (AliVZERO*) runLoader->GetAliRun()->GetDetector("VZERO");
    if (!vzero) {
      Error("Exec", "No VZERO detector for input %d", iInput);
      continue;}
      
    loader->LoadHits();
    TTree* treeH = loader->TreeH();
    if (!treeH) {
      Error("Exec", "Cannot get TreeH for input %d", iInput);
      continue; }
    
    Float_t timeV0 = 1e12;      
    for(Int_t i=0; i<80; i++) { map[i]  = 0; time[i] = 0.0; }
	      
    TClonesArray* hits = vzero->Hits();
             
//  Now makes Digits from hits
         
    Int_t nTracks = (Int_t) treeH->GetEntries();
    for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
      vzero->ResetHits();
      treeH->GetEvent(iTrack);
      Int_t nHits = hits->GetEntriesFast();
      for (Int_t iHit = 0; iHit < nHits; iHit++) {
	AliVZEROhit* hit = (AliVZEROhit *)hits->UncheckedAt(iHit);
	Int_t nPhot = hit->Nphot();
	Int_t cell  = hit->Cell();                                    
	map[cell] += nPhot;
	Float_t dt_scintillator = gRandom->Gaus(0,0.7);
	time[cell] = dt_scintillator + 1e9*hit->Tof();
	if(time[cell] < timeV0) timeV0 = time[cell];
      }           // hit   loop
    }             // track loop

    loader->UnloadHits();

  }               // input loop

// Now builds the scintillator cell response (80 cells i.e. 80 responses)
         
   for (Int_t i=0; i<80; i++) {    
       Float_t q1 = Float_t ( map[i] )* cPM * kQe;
       Float_t noise = gRandom->Gaus(10.5,3.22);
       Float_t pmResponse  =  q1/kC*TMath::Power(ktheta/kthau,1/(1-ktheta/kthau)) 
        + noise*1e-3;
       map[i] = Int_t( pmResponse * adc_gain[i]);
     }
     
 
// Now transforms 80 cell responses into 64 photomultiplier responses 
	
   for (Int_t j=0; j<32; j++){
	adc[j]  = map [j];
	time2[j]= time[j];}
	
   for (Int_t j=48; j<80; j++){
	adc[j-16]  = map [j];
	time2[j-16]= time[j];}
	
   for (Int_t j=0; j<16; j++){
	adc[16+j] = map [16 + 2*j]+ map [16 + 2*j + 1];
	time2[16+j] = TMath::Min(time [16 + 2*j], time [16 + 2*j + 1]);}
	
// Now add digits to the digit Tree 
        
   for (Int_t i=0; i<64; i++) {      
      if(adc[i] > 0) {
//    printf(" Event, cell, adc, tof = %d %d %d %f\n", 
//                  outRunLoader->GetEventNumber(),i, map[i], time[i]*100.0);
//   multiply by 10 to have 100 ps per channel :
      AddDigit(i, adc[i], Int_t(time2[i]*10.0) );
     } 
   }

  
  treeD->Fill();
  outLoader->WriteDigits("OVERWRITE");  
  outLoader->UnloadDigits();     
  ResetDigit();
}

//____________________________________________________________________________
void AliVZERODigitizer::AddDigit(Int_t PMnumber, Int_t adc, Int_t time) 
 { 
 
// Adds Digit 
 
  TClonesArray &ldigits = *fDigits;  
  new(ldigits[fNdigits++]) AliVZEROdigit(PMnumber,adc,time);
}
//____________________________________________________________________________
void AliVZERODigitizer::ResetDigit()
{
//
// Clears Digits
//
  fNdigits = 0;
  if (fDigits) fDigits->Delete();
}

//____________________________________________________________________________
AliVZEROCalibData* AliVZERODigitizer::GetCalibData() const

{
  AliCDBManager *man = AliCDBManager::Instance();

  AliCDBEntry *entry=0;

  entry = man->Get("VZERO/Calib/Data");

  if(!entry){
    AliWarning("Load of calibration data from default storage failed!");
    AliWarning("Calibration data will be loaded from local storage ($ALICE_ROOT)");
    Int_t runNumber = man->GetRun();
    entry = man->GetStorage("local://$ALICE_ROOT")
      ->Get("VZERO/Calib/Data",runNumber);
	
  }

  // Retrieval of data in directory VZERO/Calib/Data:


  AliVZEROCalibData *calibdata = 0;

  if (entry) calibdata = (AliVZEROCalibData*) entry->GetObject();
  if (!calibdata)  AliError("No calibration data from calibration database !");

  return calibdata;

}

