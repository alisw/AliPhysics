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

// $Id$

#include "AliMUONDigitizerV3.h"

#include "AliLog.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONChamber.h"
#include "AliMUONConstants.h"
#include "AliMUONData.h"
#include "AliMUONDigit.h"
#include "AliMUONTriggerDecisionV1.h"
#include "AliMUONTriggerElectronics.h"
#include "AliMUONVCalibParam.h"
#include "AliMpDEManager.h"
#include "AliMpStationType.h"
#include "AliRun.h"
#include "AliRunDigitizer.h"
#include "AliRunLoader.h"
#include "Riostream.h"
#include "TRandom.h"
#include "TString.h"

ClassImp(AliMUONDigitizerV3)

//_____________________________________________________________________________
AliMUONDigitizerV3::AliMUONDigitizerV3(AliRunDigitizer* manager, 
                                       ETriggerCodeVersion triggerCodeVersion)
: AliDigitizer(manager),
fIsInitialized(kFALSE),
fOutputData(0x0),
fCalibrationData(0x0),
fTriggerProcessor(0x0),
fTriggerCodeVersion(triggerCodeVersion)
{
  //
  // Ctor.
  //
  AliDebug(1,Form("AliRunDigitizer=%p",fManager));
}

//______________________________________________________________________________
AliMUONDigitizerV3::AliMUONDigitizerV3(const AliMUONDigitizerV3& right) 
  : AliDigitizer(right) 
{  
/// Protected copy constructor (not implemented)

  AliFatal("Copy constructor not provided.");
}

//_____________________________________________________________________________
AliMUONDigitizerV3::~AliMUONDigitizerV3()
{
  //
  // Dtor. Note we're the owner of some pointers.
  // 
  AliDebug(1,"dtor");
  delete fOutputData;
  delete fCalibrationData;
  delete fTriggerProcessor;
}

//______________________________________________________________________________
AliMUONDigitizerV3& 
AliMUONDigitizerV3::operator=(const AliMUONDigitizerV3& right)
{
/// Protected assignement operator (not implemented)

  // check assignement to self
  if (this == &right) return *this;

  AliFatal("Assignement operator not provided.");
    
  return *this;  
}    

//_____________________________________________________________________________
void 
AliMUONDigitizerV3::ApplyResponseToDigit(AliMUONDigit& digit)
{
  // For trigger digits, simply does nothing.
  //
  // For tracking digits, starting from an ideal digit's charge, we :
  //
  // - add some noise (thus leading to a realistic charge)
  // - divide by a gain (thus decalibrating the digit)
  // - add a pedestal (thus decalibrating the digit)
  // - sets the signal to zero if below 3*sigma of the noise
  //
  
  static const Int_t MAXADC = (1<<12)-1; // We code the charge on a 12 bits ADC.
  
  Int_t detElemId = digit.DetElemId();
  AliMpStationType stationType = AliMpDEManager::GetStationType(detElemId);
  if ( stationType == kStationTrigger )
  {
    return;    
  }
  
  // The following is for tracking digits only.
  
  Int_t manuId = digit.ManuId();
  Int_t manuChannel = digit.ManuChannel();
  
  AliMUONVCalibParam* pedestal = fCalibrationData->Pedestal(detElemId,manuId);
  if (!pedestal)
  {
    AliFatal(Form("Could not get pedestal for DE=%d manuId=%d",
                  detElemId,manuId));    
  }
  Float_t pedestalMean = pedestal->ValueAsFloat(manuChannel,0);
  Float_t pedestalSigma = pedestal->ValueAsFloat(manuChannel,1);
  
  Float_t adc_noise = gRandom->Gaus(0.0,pedestalSigma);
  
  AliMUONVCalibParam* gain = fCalibrationData->Gain(detElemId,manuId);
  if (!gain)
  {
    AliFatal(Form("Could not get gain for DE=%d manuId=%d",
                  detElemId,manuId));    
  }
  
  Float_t gainMean = gain->ValueAsFloat(manuChannel,0);
  Float_t signal_noise = adc_noise*gainMean;
  
  Float_t signal = digit.Signal() + signal_noise;
  Int_t adc;
  
  if ( gainMean < 1E-6 )
  {
    AliError(Form("Got a too small gain %e for DE=%d manuId=%d manuChannel=%d. "
                  "Setting signal to 0.",
                  gainMean,detElemId,manuId,manuChannel));
    adc = 0;
  }
  else
  {
    adc = TMath::Nint( signal / gainMean + pedestalMean );
    
    if ( adc <= pedestalMean + 3.0*pedestalSigma ) 
    {
      adc = 0;
    }
  }
  // be sure we stick to 12 bits.
  if ( adc > MAXADC )
  {
    adc = MAXADC;
  }
  digit.SetPhysicsSignal(TMath::Nint(signal));
  digit.SetSignal(adc);
  digit.SetADC(adc);
}

//_____________________________________________________________________________
void
AliMUONDigitizerV3::ApplyResponse()
{
  //
  // Loop over all chamber digits, and apply the response to them
  // Note that this method may remove digits.
  //
  for ( Int_t ich = 0; ich < AliMUONConstants::NCh(); ++ich )
	{
    TClonesArray* digits = fOutputData->Digits(ich);
    Int_t n = digits->GetEntriesFast();
    for ( Int_t i = 0; i < n; ++i )
    {
      AliMUONDigit* d = static_cast<AliMUONDigit*>(digits->UncheckedAt(i));
      ApplyResponseToDigit(*d);
      if ( d->Signal() <= 0 )
      {
        digits->RemoveAt(i);
      }
    }
    digits->Compress();
  }    
}

//_____________________________________________________________________________
void
AliMUONDigitizerV3::AddOrUpdateDigit(TClonesArray& array, 
                                     const AliMUONDigit& digit)
{
  //
  // Add or update a digit, depending on whether there's already a digit
  // for the corresponding channel.
  //
  Int_t ix = FindDigitIndex(array,digit);
  
  if (ix>=0)
  {
    AliMUONDigit* d = static_cast<AliMUONDigit*>(array.UncheckedAt(ix));
    Bool_t ok = MergeDigits(digit,*d);
    if (!ok)
    {
      AliError("Digits are not mergeable !");
    }
  }
  else
  {
    ix = array.GetLast() + 1;
    new(array[ix]) AliMUONDigit(digit);
  }
  
}

//_____________________________________________________________________________
void
AliMUONDigitizerV3::Exec(Option_t*)
{
  //
  // Main method.
  // We first loop over input files, and merge the sdigits we found there.
  // We then digitize all the resulting sdigits
  // And we finally generate the trigger outputs.
  //
  AliDebug(1, "Running digitizer.");
  
  if ( fManager->GetNinputs() == 0 )
  {
    AliWarning("No input set. Nothing to do.");
    return;
  }
  
  if ( !fIsInitialized )
  {
    AliError("Not initialized. Cannot perform the work. Sorry");
    return;
  }
  
  Int_t nInputFiles = fManager->GetNinputs();
  
  if ( fOutputData->TreeD() == 0x0 )
  {
    AliDebug(1,"Calling MakeDigitsContainer");
    fOutputData->GetLoader()->MakeDigitsContainer();
  }
  fOutputData->MakeBranch("D,GLT");
  fOutputData->SetTreeAddress("D,GLT");
  
  // Loop over all the input files, and merge the sdigits found in those
  // files.
  for ( Int_t iFile = 0; iFile < nInputFiles; ++iFile )
  {    
    AliMUONData* inputData = GetDataAccess(fManager->GetInputFolderName(iFile));
    if (!inputData)
    {
      AliFatal(Form("Could not get access to input file #%d",iFile));
    }

    inputData->GetLoader()->LoadSDigits("READ");
    inputData->SetTreeAddress("S");
    inputData->GetSDigits();

    MergeWithSDigits(*fOutputData,*inputData,fManager->GetMask(iFile));
    
    inputData->ResetSDigits();
    inputData->GetLoader()->UnloadSDigits();
    delete inputData;
  }
  
  // At this point, we do have digit arrays (one per chamber) which contains 
  // the merging of all the sdigits of the input file(s).
  // We now massage them to apply the detector response, i.e. this
  // is here that we do the "digitization" work.
  
  ApplyResponse();

  // We generate the global and local trigger decisions.
  fTriggerProcessor->ExecuteTask();
  
  // Fill the output treeD
  fOutputData->Fill("D,GLT");
  
  // Write to the output tree(D).
  // Please note that as GlobalTrigger, LocalTrigger and Digits are in the same
  // tree (=TreeD) in different branches, this WriteDigits in fact writes all of 
  // the 3 branches.
  fOutputData->GetLoader()->WriteDigits("OVERWRITE");
  
  // Finally, we clean up after ourselves.
  fOutputData->ResetDigits();
  fOutputData->ResetTrigger();
  fOutputData->GetLoader()->UnloadDigits();
}

//_____________________________________________________________________________
Int_t
AliMUONDigitizerV3::FindDigitIndex(TClonesArray& array, const AliMUONDigit& digit)
{
  // 
  // Return the index of digit within array, if that digit is there, 
  // otherwise returns -1
  //
  // FIXME: this is of course not the best implementation you can think of.
  // Reconsider the use of hit/digit map... ? (but be sure it's needed!)
  //
  Int_t n = array.GetEntriesFast();
  for ( Int_t i = 0; i < n; ++i )
  {
    AliMUONDigit* d = static_cast<AliMUONDigit*>(array.UncheckedAt(i));
    if ( d->DetElemId() == digit.DetElemId() &&
         d->PadX() == digit.PadX() &&
         d->PadY() == digit.PadY() && 
         d->Cathode() == digit.Cathode() )
    {
      return i;
    }
  }
  return -1;
}

//_____________________________________________________________________________
AliMUONData* 
AliMUONDigitizerV3::GetDataAccess(const TString& folderName)
{
  //
  // Create an AliMUONData to deal with data found in folderName.
  //
  AliDebug(1,Form("Getting access to folder %s",folderName.Data()));
  AliRunLoader* runLoader = AliRunLoader::GetRunLoader(folderName);
  if (!runLoader)
  {
    AliError(Form("Could not get RunLoader from folder %s",folderName.Data()));
    return 0x0;
  }
  AliLoader* loader = static_cast<AliLoader*>(runLoader->GetLoader("MUONLoader"));
  if (!loader)
  {
    AliError(Form("Could not get MuonLoader from folder %s",folderName.Data()));
    return 0x0;
  }
  AliMUONData* data = new AliMUONData(loader,"MUON","MUONDataForDigitOutput");
  AliDebug(1,Form("AliMUONData=%p loader=%p",data,loader));
  return data;
}

//_____________________________________________________________________________
Bool_t
AliMUONDigitizerV3::Init()
{
  //
  // Initialization of the TTask :
  // a) set the outputData pointer
  // b) create the calibrationData, according to run number
  // c) create the trigger processing task
  //
  AliDebug(1,"");
  
  if ( fIsInitialized )
  {
    AliError("Object already initialized.");
    return kFALSE;
  }
  
  if (!fManager)
  {
    AliError("fManager is null !");
    return kFALSE;
  }
  
  fOutputData = GetDataAccess(fManager->GetOutputFolderName());
  if (!fOutputData)
  {
    AliError("Can not perform digitization. I'm sorry");
    return kFALSE;
  }
  AliDebug(1,Form("fOutputData=%p",fOutputData));
  
  AliRunLoader* runLoader = fOutputData->GetLoader()->GetRunLoader();
  AliRun* galice = runLoader->GetAliRun();  
  Int_t runnumber = galice->GetRunNumber();
  
  fCalibrationData = new AliMUONCalibrationData(runnumber);
  
  switch (fTriggerCodeVersion)
  {
    case kTriggerDecision:
      fTriggerProcessor = new AliMUONTriggerDecisionV1(fOutputData);
      break;
    case kTriggerElectronics:
      fTriggerProcessor = new AliMUONTriggerElectronics(fOutputData);
      break;
    default:
      AliFatal("Unknown trigger processor type");
      break;
  }
  AliDebug(1,Form("Using the following trigger code %s - %s",
                  fTriggerProcessor->GetName(),fTriggerProcessor->GetTitle()));
  fIsInitialized = kTRUE;
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t
AliMUONDigitizerV3::MergeDigits(const AliMUONDigit& src, 
                                AliMUONDigit& srcAndDest)
{
  //
  // Merge 2 digits (src and srcAndDest) into srcAndDest.
  //
  AliDebug(1,"Merging the following digits:");
  StdoutToAliDebug(1,src.Print("tracks"););
  StdoutToAliDebug(1,srcAndDest.Print("tracks"););
  
  Bool_t check = ( src.DetElemId() == srcAndDest.DetElemId() &&
                   src.PadX() == srcAndDest.PadX() &&
                   src.PadY() == srcAndDest.PadY() &&
                   src.Cathode() == srcAndDest.Cathode() );
  if (!check)
  {
    return kFALSE;
  }
  
  srcAndDest.AddSignal(src.Signal());
  srcAndDest.AddPhysicsSignal(src.Physics());
  for ( Int_t i = 0; i < src.Ntracks(); ++i )
  {
    srcAndDest.AddTrack(src.Track(i),src.TrackCharge(i));
  }
  StdoutToAliDebug(1,cout << "result:"; srcAndDest.Print("tracks"););
  return kTRUE;
}

//_____________________________________________________________________________
void 
AliMUONDigitizerV3::MergeWithSDigits(AliMUONData& outputData, 
                                     const AliMUONData& inputData, Int_t mask)
{
  //
  // Merge the sdigits in inputData with the digits already present in outputData
  //
  AliDebug(1,"");
  
	for ( Int_t ich = 0; ich < AliMUONConstants::NCh(); ++ich )
	{
    TClonesArray* iDigits = inputData.SDigits(ich); 
    TClonesArray* oDigits = outputData.Digits(ich);
    if (!iDigits)
    {
      AliError(Form("Could not get sdigits for ich=%d",ich));
      return;
    }
    Int_t nSDigits = iDigits->GetEntriesFast();
    for ( Int_t k = 0; k < nSDigits; ++k )
		{
			AliMUONDigit* sdigit = static_cast<AliMUONDigit*>(iDigits->UncheckedAt(k));
      if (!sdigit)
      {
        AliError(Form("Could not get sdigit for ich=%d and k=%d",ich,k));
      }
      else
      {
        // Update the track references using the mask.
        // FIXME: this is dirty, for backward compatibility only.
        // Should re-design all this way of keeping track of MC information...
        if ( mask ) sdigit->PatchTracks(mask);
        // Then add or update the digit to the output.
        AddOrUpdateDigit(*oDigits,*sdigit);
      }
		}   
  }
}
