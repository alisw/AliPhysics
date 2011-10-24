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

#include <TH1F.h>
#include <TF1.h>
#include <TParameter.h>

#include "AliRunLoader.h"
#include "AliRawReader.h"
#include "AliGRPObject.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliVZEROReconstructor.h"
#include "AliVZERORawStream.h"
#include "AliVZEROConst.h"
#include "AliESDEvent.h"
#include "AliVZEROTriggerMask.h"
#include "AliESDfriend.h"
#include "AliESDVZEROfriend.h"
#include "AliVZEROdigit.h"
#include "AliVZEROCalibData.h"
#include "AliRunInfo.h"
#include "AliCTPTimeParams.h"
#include "AliLHCClockPhase.h"

ClassImp(AliVZEROReconstructor)

//_____________________________________________________________________________
AliVZEROReconstructor:: AliVZEROReconstructor(): AliReconstructor(),
                        fESDVZERO(0x0),
                        fESD(0x0),
                        fESDVZEROfriend(0x0),
                        fCalibData(NULL),
                        fTriggerData(NULL),
                        fTimeSlewing(NULL),
                        fSaturationCorr(NULL),
                        fEqFactors(NULL),
                        fCollisionMode(0),
                        fBeamEnergy(0.),
                        fDigitsArray(0)
{
  // Default constructor  
  // Get calibration data
  
  fCalibData = GetCalibData();

  AliCDBEntry *entry = AliCDBManager::Instance()->Get("GRP/CTP/CTPtiming");
  if (!entry) AliFatal("CTP timing parameters are not found in OCDB !");
  AliCTPTimeParams *ctpParams = (AliCTPTimeParams*)entry->GetObject();
  Float_t l1Delay = (Float_t)ctpParams->GetDelayL1L0()*25.0;

  AliCDBEntry *entry1 = AliCDBManager::Instance()->Get("GRP/CTP/TimeAlign");
  if (!entry1) AliFatal("CTP time-alignment is not found in OCDB !");
  AliCTPTimeParams *ctpTimeAlign = (AliCTPTimeParams*)entry1->GetObject();
  l1Delay += ((Float_t)ctpTimeAlign->GetDelayL1L0()*25.0);

  AliCDBEntry *entry2 = AliCDBManager::Instance()->Get("VZERO/Calib/TimeDelays");
  if (!entry2) AliFatal("VZERO time delays are not found in OCDB !");
  TH1F *delays = (TH1F*)entry2->GetObject();

  AliCDBEntry *entry3 = AliCDBManager::Instance()->Get("VZERO/Calib/TimeSlewing");
  if (!entry3) AliFatal("VZERO time slewing function is not found in OCDB !");
  fTimeSlewing = (TF1*)entry3->GetObject();

  AliCDBEntry *entry4 = AliCDBManager::Instance()->Get("GRP/Calib/LHCClockPhase");
  if (!entry4) AliFatal("LHC clock-phase shift is not found in OCDB !");
  AliLHCClockPhase *phase = (AliLHCClockPhase*)entry4->GetObject();

  for(Int_t i = 0 ; i < 64; ++i) {
    Int_t board = AliVZEROCalibData::GetBoardNumber(i);
    fTimeOffset[i] = (((Float_t)fCalibData->GetRollOver(board)-
		       (Float_t)fCalibData->GetTriggerCountOffset(board))*25.0+
		       fCalibData->GetTimeOffset(i)-
		       l1Delay-
		       phase->GetMeanPhase()+
		       delays->GetBinContent(i+1)+
		       kV0Offset);
  }

  AliCDBEntry *entry5 =  AliCDBManager::Instance()->Get("VZERO/Calib/Saturation");
  if (!entry5) AliFatal("Saturation entry is not found in OCDB !");
  fSaturationCorr = (TObjArray*)entry5->GetObject();

  AliCDBEntry *entry6 = AliCDBManager::Instance()->Get("VZERO/Trigger/Data");
  if (!entry6) AliFatal("VZERO trigger config data is not found in OCDB !");
  fTriggerData = (AliVZEROTriggerData*)entry6->GetObject();

  AliCDBEntry *entry7 = AliCDBManager::Instance()->Get("VZERO/Calib/EqualizationFactors");
  if (!entry7) AliFatal("VZERO equalization factors are not found in OCDB !");
  fEqFactors = (TH1F*)entry7->GetObject();
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

  if(fESDVZERO)
    delete fESDVZERO;
  if(fESDVZEROfriend)
   delete fESDVZEROfriend;
  if(fDigitsArray)
   delete fDigitsArray;
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
// converts RAW to digits 

  if (!digitsTree) {
    AliError("No digits tree!");
    return;
  }

  if (!fDigitsArray)
    fDigitsArray = new TClonesArray("AliVZEROdigit", 64);
  digitsTree->Branch("VZERODigit", &fDigitsArray);

  fESDVZEROfriend->Reset();

  rawReader->Reset();
  AliVZERORawStream rawStream(rawReader);
  if (rawStream.Next()) { 

    Int_t aBBflagsV0A = 0;
    Int_t aBBflagsV0C = 0;
    Int_t aBGflagsV0A = 0;
    Int_t aBGflagsV0C = 0;

    for(Int_t iChannel=0; iChannel < 64; ++iChannel) {
      Int_t offlineCh = rawStream.GetOfflineChannel(iChannel);
      // ADC charge samples
      Short_t chargeADC[AliVZEROdigit::kNClocks];
      for(Int_t iClock=0; iClock < AliVZEROdigit::kNClocks; ++iClock) {
	chargeADC[iClock] = rawStream.GetPedestal(iChannel,iClock);
      }
      // Integrator flag
      Bool_t integrator = rawStream.GetIntegratorFlag(iChannel,AliVZEROdigit::kNClocks/2);
      // Beam-beam and beam-gas flags
      if(offlineCh<32) {
	if (rawStream.GetBBFlag(iChannel,AliVZEROdigit::kNClocks/2)) aBBflagsV0C |= (1 << offlineCh);
	if (rawStream.GetBGFlag(iChannel,AliVZEROdigit::kNClocks/2)) aBGflagsV0C |= (1 << offlineCh);
      } else {
	if (rawStream.GetBBFlag(iChannel,AliVZEROdigit::kNClocks/2)) aBBflagsV0A |= (1 << (offlineCh-32));
	if (rawStream.GetBGFlag(iChannel,AliVZEROdigit::kNClocks/2)) aBGflagsV0A |= (1 << (offlineCh-32));
      }
      // HPTDC data (leading time and width)
      Int_t board = AliVZEROCalibData::GetBoardNumber(offlineCh);
      Float_t time = rawStream.GetTime(iChannel)*fCalibData->GetTimeResolution(board);
      Float_t width = rawStream.GetWidth(iChannel)*fCalibData->GetWidthResolution(board);
      // Add a digit
      if(!fCalibData->IsChannelDead(iChannel)){
	  new ((*fDigitsArray)[fDigitsArray->GetEntriesFast()])
	    AliVZEROdigit(offlineCh, time,
			  width,integrator,
			  chargeADC);
      }

      // Filling the part of esd friend object that is available only for raw data
      fESDVZEROfriend->SetBBScalers(offlineCh,rawStream.GetBBScalers(iChannel));
      fESDVZEROfriend->SetBGScalers(offlineCh,rawStream.GetBGScalers(iChannel));
      for (Int_t iBunch = 0; iBunch < AliESDVZEROfriend::kNBunches; iBunch++) {
	fESDVZEROfriend->SetChargeMB(offlineCh,iBunch,rawStream.GetChargeMB(iChannel,iBunch));
	fESDVZEROfriend->SetIntMBFlag(offlineCh,iBunch,rawStream.GetIntMBFlag(iChannel,iBunch));
	fESDVZEROfriend->SetBBMBFlag(offlineCh,iBunch,rawStream.GetBBMBFlag(iChannel,iBunch));
	fESDVZEROfriend->SetBGMBFlag(offlineCh,iBunch,rawStream.GetBGMBFlag(iChannel,iBunch));
      }
      for (Int_t iEv = 0; iEv < AliESDVZEROfriend::kNEvOfInt; iEv++) {
	  fESDVZEROfriend->SetBBFlag(offlineCh,iEv,rawStream.GetBBFlag(iChannel,iEv));
	  fESDVZEROfriend->SetBGFlag(offlineCh,iEv,rawStream.GetBGFlag(iChannel,iEv));
      }
    }  

    // Filling the global part of esd friend object that is available only for raw data
    rawStream.FillTriggerBits(fTriggerData);
    fESDVZEROfriend->SetTriggerInputs(rawStream.GetTriggerInputs());
    fESDVZEROfriend->SetTriggerInputsMask(rawStream.GetTriggerInputsMask());

    for(Int_t iScaler = 0; iScaler < AliESDVZEROfriend::kNScalers; iScaler++)
      fESDVZEROfriend->SetTriggerScalers(iScaler,rawStream.GetTriggerScalers(iScaler));

    for(Int_t iBunch = 0; iBunch < AliESDVZEROfriend::kNBunches; iBunch++)
      fESDVZEROfriend->SetBunchNumbersMB(iBunch,rawStream.GetBunchNumbersMB(iBunch));
     
    // Store the BB and BG flags in the digits tree (user info)
    digitsTree->GetUserInfo()->Add(new TParameter<int>("BBflagsV0A",aBBflagsV0A));
    digitsTree->GetUserInfo()->Add(new TParameter<int>("BBflagsV0C",aBBflagsV0C));
    digitsTree->GetUserInfo()->Add(new TParameter<int>("BGflagsV0A",aBGflagsV0A));
    digitsTree->GetUserInfo()->Add(new TParameter<int>("BGflagsV0C",aBGflagsV0C));

    UShort_t chargeA,chargeC;
    rawStream.CalculateChargeForCentrTriggers(fTriggerData,chargeA,chargeC);
    digitsTree->GetUserInfo()->Add(new TParameter<int>("ChargeA",(Int_t)chargeA));
    digitsTree->GetUserInfo()->Add(new TParameter<int>("ChargeC",(Int_t)chargeC));
    digitsTree->GetUserInfo()->Add(new TParameter<int>("TriggerInputs",(Int_t)rawStream.GetTriggerInputs()));

    digitsTree->Fill();
  }

  fDigitsArray->Clear();
}      

//______________________________________________________________________
void AliVZEROReconstructor::FillESD(TTree* digitsTree, TTree* /*clustersTree*/,
				    AliESDEvent* esd) const
{
// fills multiplicities to the ESD - pedestal is now subtracted
    
  if (!digitsTree) {
      AliError("No digits tree!");
      return;
  }

  TBranch* digitBranch = digitsTree->GetBranch("VZERODigit");
  digitBranch->SetAddress(&fDigitsArray);

  Float_t   mult[64];  
  Float_t    adc[64]; 
  Float_t   time[64]; 
  Float_t  width[64];
  Bool_t aBBflag[64];
  Bool_t aBGflag[64];
   
  for (Int_t i=0; i<64; i++){
       adc[i]    = 0.0;
       mult[i]   = 0.0;
       time[i]   = kInvalidTime;
       width[i]  = 0.0;
       aBBflag[i] = kFALSE;
       aBGflag[i] = kFALSE;
  }
     
  Int_t aBBflagsV0A = 0;
  Int_t aBBflagsV0C = 0;
  Int_t aBGflagsV0A = 0;
  Int_t aBGflagsV0C = 0;

  if (digitsTree->GetUserInfo()->FindObject("BBflagsV0A")) {
    aBBflagsV0A = ((TParameter<int>*)digitsTree->GetUserInfo()->FindObject("BBflagsV0A"))->GetVal();
  }
  else {
    if (esd && (esd->GetEventType() == 7))
      AliWarning("V0A beam-beam flags not found in digits tree UserInfo!");
  }

  if (digitsTree->GetUserInfo()->FindObject("BBflagsV0C")) {
    aBBflagsV0C = ((TParameter<int>*)digitsTree->GetUserInfo()->FindObject("BBflagsV0C"))->GetVal();
  }
  else {
    if (esd && (esd->GetEventType() == 7))
      AliWarning("V0C beam-beam flags not found in digits tree UserInfo!");
  }

  if (digitsTree->GetUserInfo()->FindObject("BGflagsV0A")) {
    aBGflagsV0A = ((TParameter<int>*)digitsTree->GetUserInfo()->FindObject("BGflagsV0A"))->GetVal();
  }
  else {
    if (esd && (esd->GetEventType() == 7))
      AliWarning("V0A beam-gas flags not found in digits tree UserInfo!");
  }

  if (digitsTree->GetUserInfo()->FindObject("BGflagsV0C")) {
    aBGflagsV0C = ((TParameter<int>*)digitsTree->GetUserInfo()->FindObject("BGflagsV0C"))->GetVal();
  }
  else {
    if (esd && (esd->GetEventType() == 7))
      AliWarning("V0C beam-gas flags not found in digits tree UserInfo!");
  }

  // Beam-beam and beam-gas flags (hardware)
  for (Int_t iChannel = 0; iChannel < 64; ++iChannel) {
    if(iChannel < 32) {
      aBBflag[iChannel] = (aBBflagsV0C >> iChannel) & 0x1;
      aBGflag[iChannel] = (aBGflagsV0C >> iChannel) & 0x1;
    }
    else {
      aBBflag[iChannel] = (aBBflagsV0A >> (iChannel-32)) & 0x1;
      aBGflag[iChannel] = (aBGflagsV0A >> (iChannel-32)) & 0x1;
    }
  }

  // Fill the trigger charges and bits
  UShort_t chargeA = 0;
  UShort_t chargeC = 0;
  UShort_t triggerInputs = 0;
  if (digitsTree->GetUserInfo()->FindObject("ChargeA")) {
    chargeA = (UShort_t)(((TParameter<int>*)digitsTree->GetUserInfo()->FindObject("ChargeA"))->GetVal());
  }
  else {
    if (esd && (esd->GetEventType() == 7))
      AliWarning("V0A trigger charge not found in digits tree UserInfo!");
  }

  if (digitsTree->GetUserInfo()->FindObject("ChargeC")) {
    chargeC = (UShort_t)(((TParameter<int>*)digitsTree->GetUserInfo()->FindObject("ChargeC"))->GetVal());
  }
  else {
    if (esd && (esd->GetEventType() == 7))
      AliWarning("V0C trigger charge not found in digits tree UserInfo!");
  }

  if (digitsTree->GetUserInfo()->FindObject("TriggerInputs")) {
    triggerInputs = (UShort_t)(((TParameter<int>*)digitsTree->GetUserInfo()->FindObject("TriggerInputs"))->GetVal());
  }
  else {
    if (esd && (esd->GetEventType() == 7))
      AliWarning("V0C trigger charge not found in digits tree UserInfo!");
  }

  fESDVZERO->SetTriggerChargeA(chargeA);
  fESDVZERO->SetTriggerChargeC(chargeC);
  fESDVZERO->SetTriggerBits(triggerInputs);
  fESDVZERO->SetBit(AliESDVZERO::kTriggerChargeBitsFilled,kTRUE);

  Int_t nEntries = (Int_t)digitsTree->GetEntries();
  for (Int_t e=0; e<nEntries; e++) {
    digitsTree->GetEvent(e);

    Int_t nDigits = fDigitsArray->GetEntriesFast();
    
    for (Int_t d=0; d<nDigits; d++) {    
        AliVZEROdigit* digit = (AliVZEROdigit*) fDigitsArray->At(d);      
        Int_t  pmNumber = digit->PMNumber();

        // Pedestal retrieval and suppression
	Bool_t integrator = digit->Integrator();
        Float_t maxadc = 0;
        Int_t imax = -1;
        Float_t adcPedSub[AliVZEROdigit::kNClocks];
        for(Int_t iClock=0; iClock < AliVZEROdigit::kNClocks; ++iClock) {
	  Short_t charge = digit->ChargeADC(iClock);
	  Bool_t iIntegrator = (iClock%2 == 0) ? integrator : !integrator;
	  Int_t k = pmNumber + 64*iIntegrator;
	  adcPedSub[iClock] = (Float_t)charge - fCalibData->GetPedestal(k);
	  if(adcPedSub[iClock] <= GetRecoParam()->GetNSigmaPed()*fCalibData->GetSigma(k)) {
	    adcPedSub[iClock] = 0;
	    continue;
	  }
	  if(iClock < GetRecoParam()->GetStartClock() || iClock > GetRecoParam()->GetEndClock()) continue;
	  if(adcPedSub[iClock] > maxadc) {
	    maxadc = adcPedSub[iClock];
	    imax   = iClock;
	  }
	}

	if (imax != -1) {
	  Int_t start = imax - GetRecoParam()->GetNPreClocks();
	  if (start < 0) start = 0;
	  Int_t end = imax + GetRecoParam()->GetNPostClocks();
	  if (end > 20) end = 20;
	  for(Int_t iClock = start; iClock <= end; iClock++) {
	    adc[pmNumber] += adcPedSub[iClock];
	  }
	}

	// HPTDC leading time and width
	// Correction for slewing and various time delays
        time[pmNumber]  =  CorrectLeadingTime(pmNumber,digit->Time(),adc[pmNumber]);
	width[pmNumber] =  digit->Width();

	if (adc[pmNumber] > 0) {
	  AliDebug(1,Form("PM = %d ADC = %f TDC %f (%f)   Int %d (%d %d %d %d %d)    %f %f   %f %f    %d %d",pmNumber, adc[pmNumber],
			  digit->Time(),time[pmNumber],
			  integrator,
			  digit->ChargeADC(8),digit->ChargeADC(9),digit->ChargeADC(10),
			  digit->ChargeADC(11),digit->ChargeADC(12),
			  fCalibData->GetPedestal(pmNumber),fCalibData->GetSigma(pmNumber),
			  fCalibData->GetPedestal(pmNumber+64),fCalibData->GetSigma(pmNumber+64),
			  aBBflag[pmNumber],aBGflag[pmNumber]));
	    };

	TF1 *saturationFunc = (TF1*)fSaturationCorr->UncheckedAt(pmNumber);
	if (!saturationFunc) AliFatal(Form("Saturation correction for channel %d is not found!",pmNumber));
	AliDebug(1,Form("Saturation PM=%d   %f %f",pmNumber,adc[pmNumber],saturationFunc->Eval(adc[pmNumber])));
	mult[pmNumber] = saturationFunc->Eval(adc[pmNumber])*fCalibData->GetMIPperADC(pmNumber);

	// Fill ESD friend object
	for (Int_t iEv = 0; iEv < AliESDVZEROfriend::kNEvOfInt; iEv++) {
	  fESDVZEROfriend->SetPedestal(pmNumber,iEv,(Float_t)digit->ChargeADC(iEv));
	  fESDVZEROfriend->SetIntegratorFlag(pmNumber,iEv,(iEv%2 == 0) ? integrator : !integrator);
	}
	fESDVZEROfriend->SetTime(pmNumber,digit->Time());
	fESDVZEROfriend->SetWidth(pmNumber,digit->Width());

    } // end of loop over digits
  } // end of loop over events in digits tree
         
  fESDVZERO->SetBit(AliESDVZERO::kCorrectedLeadingTime,kTRUE);
  fESDVZERO->SetMultiplicity(mult);
  fESDVZERO->SetADC(adc);
  fESDVZERO->SetTime(time);
  fESDVZERO->SetWidth(width);
  fESDVZERO->SetBit(AliESDVZERO::kOnlineBitsFilled,kTRUE);
  fESDVZERO->SetBBFlag(aBBflag);
  fESDVZERO->SetBGFlag(aBGflag);
  fESDVZERO->SetBit(AliESDVZERO::kCorrectedForSaturation,kTRUE);

  // now fill the V0 decision and channel flags
  {
    AliVZEROTriggerMask triggerMask;
    triggerMask.SetRecoParam(GetRecoParam());
    triggerMask.FillMasks(fESDVZERO, fCalibData, fTimeSlewing);
  }

  if (esd) { 
     AliDebug(1, Form("Writing VZERO data to ESD tree"));
     esd->SetVZEROData(fESDVZERO);
     const AliESDRun *esdRun = esd->GetESDRun();
     if (esdRun) {
       Float_t factors[64];
       Float_t factorSum = 0;
       for(Int_t i = 0; i < 64; ++i) {
	 factors[i] = fEqFactors->GetBinContent(i+1)*fCalibData->GetMIPperADC(i);
	 factorSum += factors[i];
       }
       for(Int_t i = 0; i < 64; ++i) factors[i] *= (64./factorSum);
       
       esd->SetVZEROEqFactors(factors);
     }
     else
       AliError("AliESDRun object is not available! Cannot write the equalization factors!");
  }

  if (esd) {
     AliESDfriend *fr = (AliESDfriend*)esd->FindListObject("AliESDfriend");
     if (fr) {
        AliDebug(1, Form("Writing VZERO friend data to ESD tree"));
        fr->SetVZEROfriend(fESDVZEROfriend);
    }
  }

  fDigitsArray->Clear();
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
  // Retrieval of collision mode 

  TString beamType = GetRunInfo()->GetBeamType();
  if(beamType==AliGRPObject::GetInvalidString()){
     AliError("VZERO cannot retrieve beam type");
     return;
  }

  if( (beamType.CompareTo("P-P") ==0)  || (beamType.CompareTo("p-p") ==0) ){
    fCollisionMode=0;
  }
  else if( (beamType.CompareTo("Pb-Pb") ==0)  || (beamType.CompareTo("A-A") ==0) ){
    fCollisionMode=1;
  }
    
  fBeamEnergy = GetRunInfo()->GetBeamEnergy();
  if(fBeamEnergy==AliGRPObject::GetInvalidFloat()) {
     AliError("Missing value for the beam energy ! Using 0");
     fBeamEnergy = 0.;
  }
  
  AliDebug(1,Form("\n ++++++ Beam type and collision mode retrieved as %s %d @ %1.3f GeV ++++++\n\n",beamType.Data(), fCollisionMode, fBeamEnergy));

}

//_____________________________________________________________________________
AliVZEROCalibData* AliVZEROReconstructor::GetCalibData() const
{
  // Gets calibration object for VZERO set

  AliCDBManager *man = AliCDBManager::Instance();

  AliCDBEntry *entry=0;

  entry = man->Get("VZERO/Calib/Data");

  AliVZEROCalibData *calibdata = 0;

  if (entry) calibdata = (AliVZEROCalibData*) entry->GetObject();
  if (!calibdata)  AliFatal("No calibration data from calibration database !");

  return calibdata;
}

Float_t AliVZEROReconstructor::CorrectLeadingTime(Int_t i, Float_t time, Float_t adc) const
{
  // Correct the leading time
  // for slewing effect and
  // misalignment of the channels
  if (time < 1e-6) return kInvalidTime;

  // Channel alignment and general offset subtraction
  if (i < 32) time -= kV0CDelayCables;
  time -= fTimeOffset[i];

  // In case of pathological signals
  if (adc < 1e-6) return time;

  // Slewing correction
  Float_t thr = fCalibData->GetCalibDiscriThr(i,kTRUE);
  time -= fTimeSlewing->Eval(adc/thr);

  return time;
}
