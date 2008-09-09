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
#include "AliVZEROTriggerMask.h"
#include "AliESDfriend.h"
#include "AliESDVZEROfriend.h"

ClassImp(AliVZEROReconstructor)

//_____________________________________________________________________________
AliVZEROReconstructor:: AliVZEROReconstructor(): AliReconstructor(),
   fESDVZERO(0x0),
   fESD(0x0),
   fESDVZEROfriend(0x0),
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
   delete fESDVZEROfriend;
}

//_____________________________________________________________________________
void AliVZEROReconstructor::Init()
{
// initializer

  fESDVZERO  = new AliESDVZERO;
  fESDVZEROfriend = new AliESDVZEROfriend;
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
     Int_t ADC_max[64], adc[64], time[64], width[64], BBFlag[64], BGFlag[64];   
     for(Int_t i=0; i<64; i++) {
         // Search for the maximum charge in the train of 21 LHC clocks 
         // regardless of the integrator which has been operated:
         ADC_max[i] = 0;
	 Int_t imax = 0;
         for(Int_t iClock=0; iClock<21; iClock++){
             if((Int_t)rawStream.GetPedestal(i,iClock) > ADC_max[i])  
	        {ADC_max[i]=(Int_t)rawStream.GetPedestal(i,iClock);
		 imax      = iClock;}
         }
	 // Convert i (FEE channel numbering) to j (aliroot channel numbering)
	 Int_t j   =  rawStream.GetOfflineChannel(i);
	 adc[j]    =  ADC_max[i];
	 time[j]   =  rawStream.GetTime(i);
	 width[j]  =  rawStream.GetWidth(i);
	 BBFlag[j] =  rawStream.GetBBFlag(i,imax);
	 BGFlag[j] =  rawStream.GetBGFlag(i,imax); 

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

     for (Int_t iBunch = 0; iBunch < AliESDVZEROfriend::kNBunches; iBunch++)
       fESDVZEROfriend->SetBunchNumbersMB(iBunch,rawStream.GetBunchNumbersMB(iBunch));
     

     // Channels(aliroot numbering) will be ordered in the tree
     for(Int_t iChannel = 0; iChannel < 64; iChannel++) {
         new ((*digitsArray)[digitsArray->GetEntriesFast()])
             AliVZEROdigit(iChannel, adc[iChannel], time[iChannel],
	                   width[iChannel], BBFlag[iChannel], BGFlag[iChannel]);
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

  Short_t Multiplicity[64];
  Float_t   mult[64];  
  Short_t    adc[64]; 
  Short_t   time[64]; 
  Short_t  width[64];
  Bool_t  BBFlag[64];
  Bool_t  BGFlag[64];
   
  for (Int_t i=0; i<64; i++){
       adc[i]    = 0;
       mult[i]   = 0.0;
       time[i]   = 0;
       width[i]  = 0;
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
        Int_t  pedestal      = int(fCalibData->GetPedestal(d));
        adc[pmNumber]   = (Short_t) digit->ADC() - pedestal; 
        time[pmNumber]  = (Short_t) digit->Time();
	width[pmNumber] = (Short_t) digit->Width();
	BBFlag[pmNumber]= digit->BBFlag();
	BGFlag[pmNumber]= digit->BGFlag();
        // printf("PM = %d,  MIP per ADC channel = %f \n",pmNumber, GetMIP(pmNumber));
        // cut of ADC at 1MIP/2 
        if (adc[pmNumber] > (int(1.0/GetMIP(pmNumber)) /2) ) 
	    mult[pmNumber] += float(adc[pmNumber])*GetMIP(pmNumber);
    } // end of loop over digits
  } // end of loop over events in digits tree
  
  for (Int_t j=0; j<64; j++) Multiplicity[j] = short(mult[j]+0.5); 
        
  fESDVZERO->SetMultiplicity(Multiplicity);
  fESDVZERO->SetADC(adc);
  fESDVZERO->SetTime(time);
  fESDVZERO->SetWidth(width);
  fESDVZERO->SetBBFlag(BBFlag);
  fESDVZERO->SetBGFlag(BGFlag);

  // now get the trigger mask

  AliVZEROTriggerMask *TriggerMask = new AliVZEROTriggerMask();
  TriggerMask->SetAdcThreshold(20.0/2.0);
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

//_____________________________________________________________________________
Float_t AliVZEROReconstructor::GetMIP(Int_t channel) const {

// Computes the MIP conversion factor - MIP per ADC channel - 
// Argument passed is the PM number (aliroot numbering)

  Float_t P0[64] = {
  7.094891, 7.124938, 7.089708, 7.098169, 7.094482, 7.147250, 7.170978, 7.183392, 
  7.145760, 7.148096, 7.153840, 7.143544, 7.186069, 7.194580, 7.203516, 7.195176, 
  7.188333, 7.198607, 7.209412, 7.226565, 7.221695, 7.205132, 7.191238, 7.227724, 
  7.232810, 7.252655, 7.230309, 7.273518, 7.273518, 7.242969, 7.252859, 7.252655, 
  7.026802, 7.079913, 7.134147, 7.092387, 7.079561, 7.072848, 7.123192, 7.003141, 
  7.024667, 7.124784, 7.123442, 7.129744, 7.110671, 7.143031, 7.139439, 7.178109, 
  7.247803, 7.139396, 7.293809, 7.094454, 6.992198, 7.206448, 7.244765, 7.056197, 
  7.263595, 7.138569, 7.089582, 7.215683, 7.266183, 7.165123, 7.243276, 7.235135 };
  Float_t P1[64] = {
  0.135569, 0.146405, 0.142425, 0.144278, 0.142307, 0.141648, 0.128477, 0.138239, 
  0.144173, 0.143419, 0.143572, 0.144482, 0.138024, 0.136542, 0.135955, 0.138537, 
  0.148521, 0.141999, 0.139627, 0.130014, 0.134970, 0.135635, 0.139094, 0.140634, 
  0.137971, 0.142080, 0.142793, 0.142778, 0.142778, 0.146045, 0.139133, 0.142080, 
  0.144121, 0.142311, 0.136564, 0.142686, 0.138792, 0.166285, 0.136387, 0.155391, 
  0.176082, 0.140408, 0.164738, 0.144270, 0.142766, 0.147486, 0.141951, 0.138012, 
  0.132394, 0.142849, 0.140477, 0.144592, 0.141558, 0.157646, 0.143758, 0.173385, 
  0.146489, 0.143279, 0.145230, 0.147203, 0.147333, 0.144979, 0.148597, 0.138985 };

// High Voltage retrieval from Calibration Data Base:  
  Float_t  HV = fCalibData->GetMeanHV(channel);  
  Float_t MIP = 0.5/TMath::Exp((TMath::Log(HV) - P0[channel] )/P1[channel]);
  return MIP; 

}
