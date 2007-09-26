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
#include "AliRawReader.h"
#include "AliVZEROReconstructor.h"
#include "AliVZERORawStream.h"
#include "AliESDEvent.h"

ClassImp(AliVZEROReconstructor)

//_____________________________________________________________________________
AliVZEROReconstructor:: AliVZEROReconstructor(): AliReconstructor(),
   fESDVZERO(0x0),
   fESD(0x0),
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
void AliVZEROReconstructor::Init()
{
// initializer

  fESDVZERO  = new AliESDVZERO;
}

//______________________________________________________________________
void AliVZEROReconstructor::ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const
{
// converts to digits

  if (!digitsTree) {
    AliError("No digits tree!");
    return;
  }

  TClonesArray* digitsArray = new TClonesArray("AliVZEROdigit");
  digitsTree->Branch("VZERODigit", &digitsArray);

  rawReader->Reset();
  AliVZERORawStream rawStream(rawReader);
  if (rawStream.Next()) {
    for(Int_t iChannel = 0; iChannel < 64; iChannel++) {
    Int_t adc = rawStream.GetADC(iChannel);  
    Int_t time = rawStream.GetTime(iChannel);
    new ((*digitsArray)[digitsArray->GetEntriesFast()])
      AliVZEROdigit(iChannel,adc,time);
    }
  }

  digitsTree->Fill();
}

//______________________________________________________________________
void AliVZEROReconstructor::FillESD(TTree* digitsTree, TTree* /*clustersTree*/,
				    AliESDEvent* esd) const
{
// fills multiplicities to the ESD

  if (!digitsTree) {
    AliError("No digits tree!");
    return;
  }

  TClonesArray* digitsArray = NULL;
  TBranch* digitBranch = digitsTree->GetBranch("VZERODigit");
  digitBranch->SetAddress(&digitsArray);

  Int_t   nbPMV0A = 0;
  Int_t   nbPMV0C = 0;
  Int_t   mTotV0A = 0;
  Int_t   mTotV0C = 0;
  Float_t adcV0A  = 0.0;
  Float_t adcV0C  = 0.0;
  Float_t multV0A[4];
  Float_t multV0C[4];
  Int_t   mRingV0A[4];
  Int_t   mRingV0C[4];
  
  Int_t   adc[64]; 
  Float_t mip[64];
  for (Int_t i=0; i<64; i++){
       adc[i] = 0;
       mip[i] = 110.0;}
  for (Int_t j=0; j<4; j++){
       multV0A[j]  = 0.0;
       multV0C[j]  = 0.0;
       mRingV0A[j] = 0;
       mRingV0C[j] = 0;}
     
  // loop over VZERO entries
  Int_t nEntries = (Int_t)digitsTree->GetEntries();
  for (Int_t e=0; e<nEntries; e++) {
    digitsTree->GetEvent(e);

    Int_t nDigits = digitsArray->GetEntriesFast();
    
    for (Int_t d=0; d<nDigits; d++) {    
      AliVZEROdigit* digit = (AliVZEROdigit*)digitsArray->At(d);      
      Int_t  pmNumber      = digit->PMNumber();  
      adc[pmNumber] = digit->ADC(); 
      // cut of ADC at MIP/2
      if  (adc[pmNumber] > (mip[pmNumber]/2)) { 
        if (pmNumber<=31) {
          if (pmNumber<=7) multV0C[0]=multV0C[0]+ float(adc[pmNumber])/mip[pmNumber];
	  if (pmNumber>=8  && pmNumber<=15) multV0C[1]=multV0C[1]+ float(adc[pmNumber])/mip[pmNumber];
	  if (pmNumber>=16 && pmNumber<=23) multV0C[2]=multV0C[2]+ float(adc[pmNumber])/mip[pmNumber];
	  if (pmNumber>=24 && pmNumber<=31) multV0C[3]=multV0C[3]+ float(adc[pmNumber])/mip[pmNumber];
          adcV0C = adcV0C + float(adc[pmNumber])/mip[pmNumber];
	  nbPMV0C++;
        }	
        if (pmNumber>=32 ) {
          if (pmNumber>=32 && pmNumber<=39) multV0A[0]=multV0A[0]+ float(adc[pmNumber])/mip[pmNumber];
	  if (pmNumber>=40 && pmNumber<=47) multV0A[1]=multV0A[1]+ float(adc[pmNumber])/mip[pmNumber];
	  if (pmNumber>=48 && pmNumber<=55) multV0A[2]=multV0A[2]+ float(adc[pmNumber])/mip[pmNumber];
	  if (pmNumber>=56 && pmNumber<=63) multV0A[3]=multV0A[3]+ float(adc[pmNumber])/mip[pmNumber];
          adcV0A = adcV0A + float(adc[pmNumber])/mip[pmNumber];
	  nbPMV0A++;
        }
      }
    } // end of loop over digits
    
  } // end of loop over events in digits tree
  
  mTotV0A = int(adcV0A + 0.5);
  mTotV0C = int(adcV0C + 0.5);
  for (Int_t j=0; j<4; j++){       
       mRingV0A[j] = int(multV0A[j] + 0.5);
       mRingV0C[j] = int(multV0C[j] + 0.5);}
     
  AliDebug(1,Form("VZERO multiplicities : %d (V0A) %d (V0C)", mTotV0A, mTotV0C));
  AliDebug(1,Form("Number of PMs fired  : %d (V0A) %d (V0C)", nbPMV0A, nbPMV0C));

  fESDVZERO->SetNbPMV0A(nbPMV0A);
  fESDVZERO->SetNbPMV0C(nbPMV0C);
  fESDVZERO->SetMTotV0A(mTotV0A);
  fESDVZERO->SetMTotV0C(mTotV0C);
  fESDVZERO->SetMRingV0A(mRingV0A);
  fESDVZERO->SetMRingV0C(mRingV0C);
  
  if (esd) { 
    AliDebug(1, Form("Writing VZERO data to ESD tree"));
    esd->SetVZEROData(fESDVZERO);
  }
}

//_____________________________________________________________________________
AliCDBStorage* AliVZEROReconstructor::SetStorage(const char *uri) 
{
// Sets the storage  

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
  // Gets calibration object for VZERO set

  AliCDBManager *man = AliCDBManager::Instance();

  AliCDBEntry *entry=0;

  entry = man->Get("VZERO/Calib/Data");

//   if(!entry){
//     AliWarning("Load of calibration data from default storage failed!");
//     AliWarning("Calibration data will be loaded from local storage ($ALICE_ROOT)");
//     Int_t runNumber = man->GetRun();
//     entry = man->GetStorage("local://$ALICE_ROOT")
//       ->Get("VZERO/Calib/Data",runNumber);
// 	
//   }

  // Retrieval of data in directory VZERO/Calib/Data:

  AliVZEROCalibData *calibdata = 0;

  if (entry) calibdata = (AliVZEROCalibData*) entry->GetObject();
  if (!calibdata)  AliFatal("No calibration data from calibration database !");

  return calibdata;
}
