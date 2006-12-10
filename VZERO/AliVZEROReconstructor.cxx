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

///////////////////////////////////////////////////////////////////////////////
///                                                                          //
/// class for VZERO reconstruction                                           //
///                                                                          //
///////////////////////////////////////////////////////////////////////////////

#include "AliRunLoader.h"
#include "AliVZEROReconstructor.h"
#include "AliESD.h"

ClassImp(AliVZEROReconstructor)

//_____________________________________________________________________________
AliVZEROReconstructor:: AliVZEROReconstructor(): AliReconstructor(),
   fESDVZERO(0x0),
   fESD(0x0),
   fRunLoader(0x0),
   fCalibData(GetCalibData())
{
  // Default constructor  
  // Get calibration data
  
  // fCalibData = GetCalibData(); 
}


//_____________________________________________________________________________
AliVZEROReconstructor& AliVZEROReconstructor::operator = 
  (const AliVZEROReconstructor& /*reconstructor*/)
{
// assignment operator

  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliVZEROReconstructor::~AliVZEROReconstructor()
{
// destructor
   delete fESDVZERO; 
   
}

//_____________________________________________________________________________
void AliVZEROReconstructor::Init(AliRunLoader* runLoader)
{
/// initializer

  fRunLoader = runLoader;
  fESDVZERO  = new AliESDVZERO;
}

//______________________________________________________________________
void AliVZEROReconstructor::Reconstruct(AliRunLoader* runLoader) const
{

  AliVZEROLoader* loader = (AliVZEROLoader* )runLoader->GetLoader( "VZEROLoader" );
//  AliVZEROLoader* loader = (AliVZEROLoader* )fRunLoader->GetLoader( "VZEROLoader" );
 
  loader->LoadDigits("READ");
  TTree* vzeroDigitsTree = loader->TreeD();
  if (!vzeroDigitsTree) return;

  TClonesArray* vzeroDigits = new TClonesArray("AliVZEROdigit",1000);
  TBranch* digitBranch = vzeroDigitsTree->GetBranch("VZERODigit");
  digitBranch->SetAddress(&vzeroDigits);

  Int_t   NbPMV0A = 0;
  Int_t   NbPMV0C = 0;
  Int_t   MTotV0A = 0;
  Int_t   MTotV0C = 0;
  Float_t ADCV0A  = 0.0;
  Float_t ADCV0C  = 0.0;
  Float_t MultV0A[4];
  Float_t MultV0C[4];
  Int_t   MRingV0A[4];
  Int_t   MRingV0C[4];
  
  Int_t   ADC[64]; 
  Float_t MIP[64];
  for (Int_t i=0; i<64; i++){
       ADC[i] = 0;
       MIP[i] = 110.0;}
  for (Int_t j=0; j<4; j++){
       MultV0A[j]  = 0.0;
       MultV0C[j]  = 0.0;
       MRingV0A[j] = 0;
       MRingV0C[j] = 0;}
     
  // loop over VZERO entries
  Int_t nEntries = (Int_t)vzeroDigitsTree->GetEntries();
  for (Int_t e=0; e<nEntries; e++) {
    vzeroDigitsTree->GetEvent(e);

    Int_t nDigits = vzeroDigits->GetEntriesFast();
    
    for (Int_t d=0; d<nDigits; d++) {    
      AliVZEROdigit* digit = (AliVZEROdigit*)vzeroDigits->At(d);      
      Int_t  PMNumber      = digit->PMNumber();  
      ADC[PMNumber] = digit->ADC();  
      if (PMNumber<=31) {
        if (PMNumber<=7) MultV0C[0]=MultV0C[0]+ float(ADC[PMNumber])/MIP[PMNumber];
	if (PMNumber>=8  && PMNumber<=15) MultV0C[1]=MultV0C[1]+ float(ADC[PMNumber])/MIP[PMNumber];
	if (PMNumber>=16 && PMNumber<=23) MultV0C[2]=MultV0C[2]+ float(ADC[PMNumber])/MIP[PMNumber];
	if (PMNumber>=24 && PMNumber<=31) MultV0C[3]=MultV0C[3]+ float(ADC[PMNumber])/MIP[PMNumber];
        ADCV0C = ADCV0C + float(ADC[PMNumber])/MIP[PMNumber];
	if(ADC[PMNumber] > 4) NbPMV0C++;
      }	
      if (PMNumber>=32) {
        if (PMNumber>=32 && PMNumber<=39) MultV0A[0]=MultV0A[0]+ float(ADC[PMNumber])/MIP[PMNumber];
	if (PMNumber>=40 && PMNumber<=47) MultV0A[1]=MultV0A[1]+ float(ADC[PMNumber])/MIP[PMNumber];
	if (PMNumber>=48 && PMNumber<=55) MultV0A[2]=MultV0A[2]+ float(ADC[PMNumber])/MIP[PMNumber];
	if (PMNumber>=56 && PMNumber<=63) MultV0A[3]=MultV0A[3]+ float(ADC[PMNumber])/MIP[PMNumber];
        ADCV0A = ADCV0A + float(ADC[PMNumber])/MIP[PMNumber];
	if(ADC[PMNumber] > 4) NbPMV0A++;
      }
    } // end of loop over digits
    
  } // end of loop over events in digits tree
  
  MTotV0A = int(ADCV0A + 0.5);
  MTotV0C = int(ADCV0C + 0.5);
  for (Int_t j=0; j<4; j++){       
       MRingV0A[j] = int(MultV0A[j] + 0.5);
       MRingV0C[j] = int(MultV0C[j] + 0.5);}
     
  AliDebug(1,Form("VZERO multiplicities : %d (V0A) %d (V0C)", MTotV0A, MTotV0C));
  AliDebug(1,Form("Number of PMs fired  : %d (V0A) %d (V0C)", NbPMV0A, NbPMV0C));

  fESDVZERO->SetNbPMV0A(NbPMV0A);
  fESDVZERO->SetNbPMV0C(NbPMV0C);
  fESDVZERO->SetMTotV0A(MTotV0A);
  fESDVZERO->SetMTotV0C(MTotV0C);
  fESDVZERO->SetMRingV0A(MRingV0A);
  fESDVZERO->SetMRingV0C(MRingV0C);
  
}


//_____________________________________________________________________________
void AliVZEROReconstructor::FillESD(AliRunLoader* /*runLoader*/, 
				    AliESD* esd) const
{
// fill ESD 

   if (esd) { 
      AliDebug(1, Form("Writing VZERO data to ESD tree"));
      esd->SetVZEROData(fESDVZERO);
   }
  
}

//_____________________________________________________________________________
AliCDBStorage* AliVZEROReconstructor::SetStorage(const char *uri) 
{
  
  Bool_t deleteManager = kFALSE;
  
  AliCDBManager *manager = AliCDBManager::Instance();
  AliCDBStorage *defstorage = manager->GetDefaultStorage();
  
  if(!defstorage || !(defstorage->Contains("VZERO"))){ 
     AliWarning("No default storage set or default storage doesn't contain VZERO!");
     manager->SetDefaultStorage(uri);
     deleteManager = kTRUE;
  }
 
  AliCDBStorage *storage = manager->GetDefaultStorage();

  if(deleteManager){
    AliCDBManager::Instance()->UnsetDefaultStorage();
    defstorage = 0;   // the storage is killed by AliCDBManager::Instance()->Destroy()
  }

  return storage; 
}

//_____________________________________________________________________________
AliVZEROCalibData* AliVZEROReconstructor::GetCalibData() const
{

  // Getting calibration object for VZERO set

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
