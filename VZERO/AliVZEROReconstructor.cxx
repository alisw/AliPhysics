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
#include "AliGRPObject.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliVZEROReconstructor.h"
#include "AliVZERORawStream.h"
#include "AliESDEvent.h"
#include "AliVZEROTriggerMask.h"
#include "AliESDfriend.h"
#include "AliESDVZEROfriend.h"
#include "AliVZEROdigit.h"

ClassImp(AliVZEROReconstructor)

//_____________________________________________________________________________
AliVZEROReconstructor:: AliVZEROReconstructor(): AliReconstructor(),
                        fESDVZERO(0x0),
                        fESD(0x0),
                        fESDVZEROfriend(0x0),
                        fCalibData(GetCalibData()),
                        fCollisionMode(0),
                        fBeamEnergy(0.)
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
   delete fESDVZEROfriend;
}

//_____________________________________________________________________________
void AliVZEROReconstructor::Init()
{
// initializer

  fESDVZERO  = new AliESDVZERO;
  fESDVZEROfriend = new AliESDVZEROfriend;
  
  GetCollisionMode();  // fCollisionMode =1 for Pb-Pb simulated data
}

//______________________________________________________________________
void AliVZEROReconstructor::ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const
{
// converts RAW to digits - pedestal is subtracted 

  if (!digitsTree) {
    AliError("No digits tree!");
    return;
  }

  TClonesArray* digitsArray = new TClonesArray("AliVZEROdigit");
  digitsTree->Branch("VZERODigit", &digitsArray);

  fESDVZEROfriend->Reset();

  rawReader->Reset();
  AliVZERORawStream rawStream(rawReader);
  if (rawStream.Next()) {  
     Float_t ADC_max[64], adc[64]; 
     Bool_t BBFlag[64], BGFlag[64], integrator[64]; 
     Float_t time[64], width[64];  
     for(Int_t i=0; i<64; i++) {
         // Search for the maximum charge in the train of 21 LHC clocks 
         // regardless of the integrator which has been operated:
         ADC_max[i] = 0.0;
	 Int_t imax = 0;
         for(Int_t iClock=0; iClock<21; iClock++){
             if(rawStream.GetPedestal(i,iClock) > ADC_max[i])  
	        {ADC_max[i]= rawStream.GetPedestal(i,iClock);
		 imax      = iClock;}
         }
	 // Convert i (FEE channel numbering) to j (aliroot channel numbering)
	 Int_t j       =  rawStream.GetOfflineChannel(i);
	 adc[j]        =  ADC_max[i];
	 time[j]       =  rawStream.GetTime(i);
	 width[j]      =  rawStream.GetWidth(i);
	 BBFlag[j]     =  rawStream.GetBBFlag(i,imax);
	 BGFlag[j]     =  rawStream.GetBGFlag(i,imax); 
	 integrator[j] =  rawStream.GetIntegratorFlag(i,imax); 

	 // Filling the esd friend object
	 fESDVZEROfriend->SetBBScalers(j,rawStream.GetBBScalers(i));
	 fESDVZEROfriend->SetBGScalers(j,rawStream.GetBGScalers(i));
	 for (Int_t iBunch = 0; iBunch < AliESDVZEROfriend::kNBunches; iBunch++) {
	     fESDVZEROfriend->SetChargeMB(j,iBunch,rawStream.GetChargeMB(i,iBunch));
	     fESDVZEROfriend->SetIntMBFlag(j,iBunch,rawStream.GetIntMBFlag(i,iBunch));
	     fESDVZEROfriend->SetBBMBFlag(j,iBunch,rawStream.GetBBMBFlag(i,iBunch));
	     fESDVZEROfriend->SetBGMBFlag(j,iBunch,rawStream.GetBGMBFlag(i,iBunch));
	 }
	 for (Int_t iEv = 0; iEv < AliESDVZEROfriend::kNEvOfInt; iEv++) {
	     fESDVZEROfriend->SetPedestal(j,iEv,rawStream.GetPedestal(i,iEv));
	     fESDVZEROfriend->SetIntegratorFlag(j,iEv,rawStream.GetIntegratorFlag(i,iEv));
	     fESDVZEROfriend->SetBBFlag(j,iEv,rawStream.GetBBFlag(i,iEv));
	     fESDVZEROfriend->SetBGFlag(j,iEv,rawStream.GetBGFlag(i,iEv));
	 }
	 fESDVZEROfriend->SetTime(j,rawStream.GetTime(i));
	 fESDVZEROfriend->SetWidth(j,rawStream.GetWidth(i));
     }  

     // Filling the esd friend object
     fESDVZEROfriend->SetTriggerInputs(rawStream.GetTriggerInputs());
     fESDVZEROfriend->SetTriggerInputsMask(rawStream.GetTriggerInputsMask());

     for(Int_t iScaler = 0; iScaler < AliESDVZEROfriend::kNScalers; iScaler++)
         fESDVZEROfriend->SetTriggerScalers(iScaler,rawStream.GetTriggerScalers(iScaler));

     for(Int_t iBunch = 0; iBunch < AliESDVZEROfriend::kNBunches; iBunch++)
         fESDVZEROfriend->SetBunchNumbersMB(iBunch,rawStream.GetBunchNumbersMB(iBunch));
     

     // Channels(aliroot numbering) will be ordered in the tree
     for(Int_t iChannel = 0; iChannel < 64; iChannel++) {
         new ((*digitsArray)[digitsArray->GetEntriesFast()])
             AliVZEROdigit(iChannel, adc[iChannel], time[iChannel],
	                   width[iChannel], BBFlag[iChannel], BGFlag[iChannel],integrator[iChannel]);
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

  Float_t   mult[64];  
  Float_t    adc[64]; 
  Float_t   time[64]; 
  Float_t  width[64];
  Bool_t  BBFlag[64];
  Bool_t  BGFlag[64];
   
  for (Int_t i=0; i<64; i++){
       adc[i]    = 0.0;
       mult[i]   = 0.0;
       time[i]   = 0.0;
       width[i]  = 0.0;
       BBFlag[i] = kFALSE;
       BGFlag[i] = kFALSE;
  }
     
  // loop over VZERO entries to get multiplicity
  Int_t nEntries = (Int_t)digitsTree->GetEntries();
  for (Int_t e=0; e<nEntries; e++) {
    digitsTree->GetEvent(e);

    Int_t nDigits = digitsArray->GetEntriesFast();
    
    for (Int_t d=0; d<nDigits; d++) {    
        AliVZEROdigit* digit = (AliVZEROdigit*)digitsArray->At(d);      
        Int_t  pmNumber      = digit->PMNumber(); 
        // Pedestal retrieval and suppression: 
        Float_t  pedestal      = fCalibData->GetPedestal(d);
        adc[pmNumber]   =  digit->ADC() - pedestal; 
        time[pmNumber]  =  digit->Time();
	width[pmNumber] =  digit->Width();
	BBFlag[pmNumber]=  digit->BBFlag();
	BGFlag[pmNumber]=  digit->BGFlag();
        // printf("PM = %d,  MIP per ADC channel = %f \n",pmNumber, fCalibData->GetMIPperADC(pmNumber));
		//AliInfo(Form("PM = %d,  ADC = %f TDC %f",pmNumber, digit->ADC(),digit->Time()));
        // cut of ADC at 1MIP/2 
	if(fCollisionMode >0) { 
	   Float_t MIP = 2.0;  
           if (adc[pmNumber] > (MIP/2.) ) mult[pmNumber] += adc[pmNumber]*(1.0/MIP) ;
	   }
           else{    
           if (adc[pmNumber] > (1.0/fCalibData->GetMIPperADC(pmNumber) /2.) ) 
	       mult[pmNumber] += adc[pmNumber]*fCalibData->GetMIPperADC(pmNumber);
        } 	    
    } // end of loop over digits
  } // end of loop over events in digits tree
         
  fESDVZERO->SetMultiplicity(mult);
  fESDVZERO->SetADC(adc);
  fESDVZERO->SetTime(time);
  fESDVZERO->SetWidth(width);
  fESDVZERO->SetBBFlag(BBFlag);
  fESDVZERO->SetBGFlag(BGFlag);

  // now get the trigger mask

  AliVZEROTriggerMask *TriggerMask = new AliVZEROTriggerMask();
  TriggerMask->SetAdcThreshold(10.0/2.0);
  TriggerMask->SetTimeWindowWidthBBA(50);
  TriggerMask->SetTimeWindowWidthBGA(20);
  TriggerMask->SetTimeWindowWidthBBC(50);
  TriggerMask->SetTimeWindowWidthBGC(20);
  TriggerMask->FillMasks(digitsTree,digitsArray);

  fESDVZERO->SetBBtriggerV0A(TriggerMask->GetBBtriggerV0A());
  fESDVZERO->SetBGtriggerV0A(TriggerMask->GetBGtriggerV0A());
  fESDVZERO->SetBBtriggerV0C(TriggerMask->GetBBtriggerV0C());
  fESDVZERO->SetBGtriggerV0C(TriggerMask->GetBGtriggerV0C());
  
  if (esd) { 
     AliDebug(1, Form("Writing VZERO data to ESD tree"));
     esd->SetVZEROData(fESDVZERO);
  }

  if (esd) {
     AliESDfriend *fr = (AliESDfriend*)esd->FindListObject("AliESDfriend");
     if (fr) {
        AliDebug(1, Form("Writing VZERO friend data to ESD tree"));
        fr->SetVZEROfriend(fESDVZEROfriend);
    }
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

//____________________________________________________________________________
void AliVZEROReconstructor::GetCollisionMode()
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

   if(!grpData) { AliError("No GRP entry found in OCDB!");
                  return; }

// Retrieval of simulated collision mode 

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
  
//    printf("\n ++++++ Beam type and collision mode retrieved as %s %d @ %1.3f GeV ++++++\n\n",beamType.Data(), fCollisionMode, fBeamEnergy);

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

