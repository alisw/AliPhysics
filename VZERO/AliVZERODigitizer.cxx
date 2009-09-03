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

///_________________________________________________________________________
///
/// This class constructs Digits out of Hits
///
///

// --- Standard library ---

// --- ROOT system ---
#include <TMath.h>
#include <TTree.h>
#include <TMap.h>
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>
#include <AliGeomManager.h>
#include <TRandom.h>

// --- AliRoot header files ---
#include "AliVZEROConst.h"
#include "AliRun.h"
#include "AliVZERO.h"
#include "AliVZEROhit.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliGRPObject.h"
#include "AliRunDigitizer.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliVZEROCalibData.h"

#include "AliVZEROdigit.h"
#include "AliVZERODigitizer.h"

ClassImp(AliVZERODigitizer)

 AliVZERODigitizer::AliVZERODigitizer()
                   :AliDigitizer(),
                    fCalibData(GetCalibData()),
                    fPhotoCathodeEfficiency(0.18),
                    fPMVoltage(768.0),
                    fPMGain(TMath::Power((fPMVoltage / 112.5) ,7.04277)),
                    fNdigits(0),
                    fDigits(0),
                    fCollisionMode(0),
                    fBeamEnergy(0.)
   
{
  // default constructor

//    fNdigits = 0;
//    fDigits  = 0;
//   
//    fPhotoCathodeEfficiency =   0.18;
//    fPMVoltage              =  768.0;
//    fPMGain = TMath::Power((fPMVoltage / 112.5) ,7.04277); 
   
//   fCalibData = GetCalibData();
}

//____________________________________________________________________________ 
  AliVZERODigitizer::AliVZERODigitizer(AliRunDigitizer* manager)
                    :AliDigitizer(manager),
		     fCalibData(GetCalibData()),
                     fPhotoCathodeEfficiency(0.18),
                     fPMVoltage(768.0),
                     fPMGain(TMath::Power((fPMVoltage / 112.5) ,7.04277)),
		     fNdigits(0),
                     fDigits(0),
		     fCollisionMode(0),
                     fBeamEnergy(0.)
		                        
{
  // constructor
  
//   fNdigits = 0;
//   fDigits  = 0;
//   
//   fPhotoCathodeEfficiency =   0.18;
//   fPMVoltage              =  768.0;
//   fPMGain = TMath::Power( (fPMVoltage / 112.5) ,7.04277 );
  
//  fCalibData = GetCalibData();
  
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
  
  //  TGeoHMatrix *im = AliGeomManager::GetMatrix("VZERO/V0C");
  //  im->Print();

  GetCollisionMode();
  return kTRUE;
}

//____________________________________________________________________________
void AliVZERODigitizer::Exec(Option_t* /*option*/) 
{   
  // Creates digits from hits
     
  Float_t     map[80];    // 48 values on V0C + 32 on V0A
//  Int_t       pmNumber[80];
  Float_t     adc[64];    // 32 PMs on V0C + 32 PMs on V0A
  Float_t     time[80], time_ref[80], time2[64];
  Float_t     adc_gain[80]; 
  Float_t     adc_pedestal[64],adc_sigma[64];    
  fNdigits     =    0;  
  Float_t     pmGain_smeared[64];  
  Float_t     cPM[80];
  
  // Smearing of the PM tubes intrinsic characteristics
  
  for(Int_t ii=0; ii<64; ii++) {
      pmGain_smeared[ii] = gRandom->Gaus(fPMGain, fPMGain/5); }
              
  // Retrieval of ADC gain values and pedestal information from local CDB 
  // I use only the first 64th values of the calibration array in CDB 
  // as I have no beam burst structure - odd or even beam burst number
  
  // Reminder : We have 16 scintillating cells mounted on 8 PMs 
  // on Ring 3 and Ring 4 in V0C -  added to produce  ADC outputs 
  // on these rings... 
   
  for(Int_t i=0; i<16; i++) { 
	adc_gain[i] = fCalibData->GetGain(i); 
	cPM[i]      = fPhotoCathodeEfficiency * pmGain_smeared[i];
  }
  
  for(Int_t j=16; j<48; j=j+2) { 
	Int_t i=(j+17)/2;
	adc_gain[j]   = fCalibData->GetGain(i);	    
	adc_gain[j+1] = fCalibData->GetGain(i); 
	cPM[j]        = fPhotoCathodeEfficiency * pmGain_smeared[i];   
	cPM[j+1]      = fPhotoCathodeEfficiency * pmGain_smeared[i]; 
  }
	    
  for(Int_t i=48; i<80; i++){ 
	adc_gain[i] = fCalibData->GetGain(i-16); 
	cPM[i]      = fPhotoCathodeEfficiency * pmGain_smeared[i-16];
  };
  
  for(Int_t  i=0; i<64; i++){ 
	  adc_pedestal[i] = fCalibData->GetPedestal(i);
	  adc_sigma[i]    = fCalibData->GetSigma(i); 
  }; 
                                
//  for(Int_t i=0; i<64; i++) { printf(" i = %d pedestal = %f sigma = %f \n\n", 
//                                       i, adc_pedestal[i], adc_sigma[i] );} 
            
  AliRunLoader* outRunLoader =  AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());    
  if (!outRunLoader) {
    Error("Exec", "Can not get output Run Loader");
    return;
  }
    
  AliLoader* outLoader = outRunLoader->GetLoader("VZEROLoader");

  if (!outLoader) {
    Error("Exec", "Can not get output VZERO Loader");
    return;
  }

  const char* mode = "update";
  if(outRunLoader->GetEventNumber() == 0) mode = "recreate";
  outLoader->LoadDigits(mode);

  if (!outLoader->TreeD()) outLoader->MakeTree("D");
  outLoader->MakeDigitsContainer();
  TTree* treeD  = outLoader->TreeD();
  Int_t bufsize = 16000;
  treeD->Branch("VZERODigit", &fDigits, bufsize); 
  
  for (Int_t iInput = 0; iInput < fManager->GetNinputs(); iInput++) {
     AliRunLoader* runLoader = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(iInput));
     AliLoader* loader = runLoader->GetLoader("VZEROLoader");
     if (!loader) {
       Error("Exec", "Can not get VZERO Loader for input %d", iInput);
       continue;
	 }
      
     if (!runLoader->GetAliRun()) runLoader->LoadgAlice();

     AliVZERO* vzero = (AliVZERO*) runLoader->GetAliRun()->GetDetector("VZERO");
     if (!vzero) {
       Error("Exec", "No VZERO detector for input %d", iInput);
       continue;
	 }
      
     loader->LoadHits();
     TTree* treeH = loader->TreeH();
     if (!treeH) {
       Error("Exec", "Cannot get TreeH for input %d", iInput);
       continue; 
	 }
       
     for(Int_t i=0; i<80; i++) {map[i] = 0; time[i] = 0.0;}
     
     TClonesArray* hits = vzero->Hits();
             
//  Now makes Digits from hits
     
     Int_t nTracks = (Int_t) treeH->GetEntries();
     for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
         for (Int_t i=0; i<80; i++) {time_ref[i] = 999999.0;}   
         vzero->ResetHits();
         treeH->GetEvent(iTrack);
         Int_t nHits = hits->GetEntriesFast();
         for (Int_t iHit = 0; iHit < nHits; iHit++) {
			 AliVZEROhit* hit = (AliVZEROhit *)hits->UncheckedAt(iHit);
			 Int_t nPhot = hit->Nphot();
			 Int_t cell  = hit->Cell();                          
			 map[cell] += Float_t(nPhot);
			 Float_t dt_scintillator = gRandom->Gaus(0,0.7);
			 Float_t t = dt_scintillator + 1e9*hit->Tof();
			 if (t > 0.0) {
				 if(t < time_ref[cell]) time_ref[cell] = t;
				 time[cell] = TMath::Min(t,time_ref[cell]); 
			 }
         }           // hit   loop      
     }             // track loop

     loader->UnloadHits();

  }               // input loop

// Now builds the scintillator cell response (80 cells i.e. 80 responses)
         
   for (Int_t i=0; i<80; i++) {    
        Float_t q1 = Float_t ( map[i] )* cPM[i] * kQe;
        Float_t noise = gRandom->Gaus(10.5,3.22);
        Float_t pmResponse  =  q1/kC*TMath::Power(ktheta/kthau,1/(1-ktheta/kthau)) 
        + noise*1e-3; 	
	if(fCollisionMode >0) adc_gain[i] = adc_gain[i]/70.0; // reduce dynamics in Ion Collision Mode
        map[i] =  pmResponse * adc_gain[i];
        Float_t MIP = 1.0/fCalibData->GetMIPperADC(GetPMNumber(i));
	if(fCollisionMode >0) MIP=2.0;
//	printf("cell = %d,  ADC = %d, TDC = %f \n",i,map[i], time[i]*10.0 );
        if(map[i] > (MIP/2.) )
	          {map[i] = gRandom->Gaus(map[i], (MIP/6.) );}
   }
      
// Now transforms 80 cell responses into 64 photomultiplier responses
// Also adds the ADC pedestals taken out of the calibration data base
	
   for (Int_t j=0; j<16; j++){
        adc[j]  = map [j] + gRandom->Gaus(adc_pedestal[j], adc_sigma[j]);
	time2[j]= time[j];}
	
   for (Int_t j=48; j<80; j++){
        adc[j-16]  = map [j] + gRandom->Gaus(adc_pedestal[j-16],adc_sigma[j-16]);
	time2[j-16]= time[j]; }
	
   for (Int_t j=0; j<16; j++){
        adc[16+j] = map [16+2*j]+ map [16+2*j+1] + gRandom->Gaus(adc_pedestal[16+j], adc_sigma[16+j]);
	Float_t min_time = TMath::Min(time [16+2*j],time [16+2*j+1]);
	time2[16+j] = min_time;
	if(min_time==0.0){time2[16+j]=TMath::Max(time[16+2*j],time[16+2*j+1]);}
   }
   	

// Now add digits to the digit Tree 
        
   for (Int_t i=0; i<64; i++) {      
      if(adc[i] > 0) {
//           printf(" Event, cell, adc, tof = %d %d %d %d\n", 
//                    outRunLoader->GetEventNumber(),i, adc[i], Int_t((time2[i]*10.0) +0.5));
//           multiply by 10 to have 100 ps per channel :
		  
		  AddDigit(i, adc[i], (time2[i]*10.0) ) ;
	  }      
   }
  treeD->Fill();
  outLoader->WriteDigits("OVERWRITE");  
  outLoader->UnloadDigits();     
  ResetDigit();
}

//____________________________________________________________________________
void AliVZERODigitizer::AddDigit(Int_t PMnumber, Float_t adc, Float_t time) 
 { 
 
// Adds Digit 
 
  TClonesArray &ldigits = *fDigits;  
  Bool_t integrator;
  if (((Int_t) gRandom->Uniform(2))<1) integrator = kFALSE;
  else integrator = kTRUE;
	 
  new(ldigits[fNdigits++]) AliVZEROdigit(PMnumber,adc,time,0,kFALSE,kFALSE,integrator);
	 
}
//____________________________________________________________________________
void AliVZERODigitizer::ResetDigit()
{

// Clears Digits

  fNdigits = 0;
  if (fDigits) fDigits->Delete();
}

//____________________________________________________________________________
void AliVZERODigitizer::GetCollisionMode()
{
// Retrieves the collision mode from GRP data

// Initialization of the GRP entry 

   Int_t run = AliCDBManager::Instance()->GetRun();
  
//   printf("\n ++++++ Run Number retrieved as %d \n",run); 
 
  AliCDBEntry*  entry = AliCDBManager::Instance()->Get("GRP/GRP/Data",run);
  AliGRPObject* grpData = 0x0;
   
  if(entry){
    TMap* m = dynamic_cast<TMap*>(entry->GetObject());  // old GRP entry
    if(m){
       m->Print();
       grpData = new AliGRPObject();
       grpData->ReadValuesFromMap(m);
    }
    else{
       grpData = dynamic_cast<AliGRPObject*>(entry->GetObject());  // new GRP entry
       entry->SetOwner(0);
    }
    AliCDBManager::Instance()->UnloadFromCache("GRP/GRP/Data");
  }

  if(!grpData) AliError("No GRP entry found in OCDB!");

// Retrieval of collision mode

  TString beamType = grpData->GetBeamType();
  if(beamType==AliGRPObject::GetInvalidString()){
     AliError("GRP/GRP/Data entry:  missing value for the beam type !");
     AliError("\t VZERO cannot retrieve beam type\n");
     return;
  }

   if( (beamType.CompareTo("P-P") ==0)  || (beamType.CompareTo("p-p") ==0) ){
       fCollisionMode=0;
  }
   else if( (beamType.CompareTo("Pb-Pb") ==0)  || (beamType.CompareTo("A-A") ==0) ){
       fCollisionMode=1;
   }
    
  fBeamEnergy = grpData->GetBeamEnergy();
  if(fBeamEnergy==AliGRPObject::GetInvalidFloat()) {
     AliError("GRP/GRP/Data entry:  missing value for the beam energy ! Using 0");
     fBeamEnergy = 0.;
  }
  
//     printf("\n ++++++ Beam type and collision mode retrieved as %s %d @ %1.3f GeV ++++++\n\n",beamType.Data(), fCollisionMode, fBeamEnergy);

}

//____________________________________________________________________________
AliVZEROCalibData* AliVZERODigitizer::GetCalibData() const

{
  AliCDBManager *man = AliCDBManager::Instance();

  AliCDBEntry *entry=0;

  entry = man->Get("VZERO/Calib/Data");

//   if(!entry){
//     AliWarning("Load of calibration data from default storage failed!");
//     AliWarning("Calibration data will be loaded from local storage ($ALICE_ROOT)");
//     Int_t runNumber = man->GetRun();
//     entry = man->GetStorage("local://$ALICE_ROOT/OCDB")
//       ->Get("VZERO/Calib/Data",runNumber);
// 	
//   }

  // Retrieval of data in directory VZERO/Calib/Data:


  AliVZEROCalibData *calibdata = 0;

  if (entry) calibdata = (AliVZEROCalibData*) entry->GetObject();
  if (!calibdata)  AliFatal("No calibration data from calibration database !");

  return calibdata;

}

//____________________________________________________________________________
Int_t AliVZERODigitizer::GetPMNumber(Int_t cell) const

{
   
  Int_t pmNumber[80] = { 0,  1,  2,  3,  4,  5,  6,  7,
                          8,  9, 10, 11, 12, 13, 14, 15, 
			 16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 
			 24, 24, 25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31,
			 32, 33, 34, 35, 36, 37, 38, 39,
		         40, 41, 42, 43, 44, 45, 46, 47, 
			 48, 49, 50, 51, 52, 53, 54, 55,
			 56, 57, 58, 59, 60, 61, 62, 63 };
			      
  return pmNumber[cell];	
}


