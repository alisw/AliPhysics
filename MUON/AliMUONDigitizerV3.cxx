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
#include "AliMUON.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONConstants.h"
#include "AliMUONData.h"
#include "AliMUONDataIterator.h"
#include "AliMUONDigit.h"
#include "AliMUONSegmentation.h"
#include "AliMUONTriggerDecisionV1.h"
#include "AliMUONTriggerEfficiencyCells.h"
#include "AliMUONTriggerElectronics.h"
#include "AliMUONVCalibParam.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpIntPair.h"
#include "AliMpPad.h"
#include "AliMpStationType.h"
#include "AliMpVSegmentation.h"
#include "AliRun.h"
#include "AliRunDigitizer.h"
#include "AliRunLoader.h"
#include "Riostream.h"
#include "TF1.h"
#include "TRandom.h"
#include "TString.h"

///
/// The digitizer is performing the transformation to go from SDigits (digits
/// w/o any electronic noise) to Digits (w/ electronic noise, and decalibration)
/// 
/// The decalibration is performed by doing the reverse operation of the
/// calibration, that is we do (Signal+pedestal)/gain -> ADC
///
/// Note also that the digitizer takes care of merging sdigits that belongs
/// to the same pad, either because we're merging several input sdigit files
/// or with a single file because the sdigitizer does not merge sdigits itself
/// (for performance reason mainly, and because anyway we know we have to do it
/// here, at the digitization level).
///

namespace
{
  AliMUON* muon()
  {
    return static_cast<AliMUON*>(gAlice->GetModule("MUON"));
  }

  AliMUONSegmentation* Segmentation()
  {
    static AliMUONSegmentation* segmentation = muon()->GetSegmentation();
    return segmentation;
  }
}

const Double_t AliMUONDigitizerV3::fgkNSigmas=3;

ClassImp(AliMUONDigitizerV3)

//_____________________________________________________________________________
AliMUONDigitizerV3::AliMUONDigitizerV3(AliRunDigitizer* manager, 
                                       ETriggerCodeVersion triggerCodeVersion,
                                       Bool_t generateNoisyDigits)
: AliDigitizer(manager),
fIsInitialized(kFALSE),
fOutputData(0x0),
fCalibrationData(0x0),
fTriggerProcessor(0x0),
fTriggerCodeVersion(triggerCodeVersion),
fTriggerEfficiency(0x0),
fFindDigitIndexTimer(),
fGenerateNoisyDigitsTimer(),
fExecTimer(),
fNoiseFunction(0x0),
fGenerateNoisyDigits(generateNoisyDigits)
{
  //
  // Ctor.
  //
  AliDebug(1,Form("AliRunDigitizer=%p",fManager));
  fGenerateNoisyDigitsTimer.Start(kTRUE); fGenerateNoisyDigitsTimer.Stop();
  fExecTimer.Start(kTRUE); fExecTimer.Stop();
  fFindDigitIndexTimer.Start(kTRUE); fFindDigitIndexTimer.Stop();
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
  delete fNoiseFunction;
  
  AliInfo(Form("Execution time for FindDigitIndex() : R:%.2fs C:%.2fs",
               fFindDigitIndexTimer.RealTime(),fFindDigitIndexTimer.CpuTime()));
  if ( fGenerateNoisyDigits )
  {
    AliInfo(Form("Execution time for GenerateNoisyDigits() : R:%.2fs C:%.2fs",
                 fGenerateNoisyDigitsTimer.RealTime(),
                 fGenerateNoisyDigitsTimer.CpuTime()));
  }
  AliInfo(Form("Execution time for Exec() : R:%.2fs C:%.2fs",
               fExecTimer.RealTime(),fExecTimer.CpuTime()));
  
}

//_____________________________________________________________________________
void 
AliMUONDigitizerV3::ApplyResponseToTrackerDigit(AliMUONDigit& digit, Bool_t addNoise)
{
  // For tracking digits, starting from an ideal digit's charge, we :
  //
  // - add some noise (thus leading to a realistic charge), if requested to do so
  // - divide by a gain (thus decalibrating the digit)
  // - add a pedestal (thus decalibrating the digit)
  // - sets the signal to zero if below 3*sigma of the noise
  //

  static const Int_t kMaxADC = (1<<12)-1; // We code the charge on a 12 bits ADC.

  Float_t signal = digit.Signal();

  if ( !addNoise )
    {
      digit.SetADC(TMath::Nint(signal));
      return;
    }
  
  Int_t detElemId = digit.DetElemId();
  
  Int_t manuId = digit.ManuId();
  Int_t manuChannel = digit.ManuChannel();
  
  AliMUONVCalibParam* pedestal = fCalibrationData->Pedestals(detElemId,manuId);
  if (!pedestal)
    {
      AliFatal(Form("Could not get pedestal for DE=%d manuId=%d",
		    detElemId,manuId));    
    }
  Float_t pedestalMean = pedestal->ValueAsFloat(manuChannel,0);
  Float_t pedestalSigma = pedestal->ValueAsFloat(manuChannel,1);
  
  AliMUONVCalibParam* gain = fCalibrationData->Gains(detElemId,manuId);
  if (!gain)
    {
      AliFatal(Form("Could not get gain for DE=%d manuId=%d",
		    detElemId,manuId));    
    }    
  Float_t gainMean = gain->ValueAsFloat(manuChannel,0);

  Float_t adcNoise = gRandom->Gaus(0.0,pedestalSigma);
     
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
      adc = TMath::Nint( signal / gainMean + pedestalMean + adcNoise);///
      
      if ( adc <= pedestalMean + fgkNSigmas*pedestalSigma ) 
	{
	  adc = 0;
	}
    }
  
  // be sure we stick to 12 bits.
  if ( adc > kMaxADC )
    {
      adc = kMaxADC;
    }
  
  digit.SetPhysicsSignal(TMath::Nint(signal));
  digit.SetSignal(adc);
  digit.SetADC(adc);
}

//_____________________________________________________________________________
void 
AliMUONDigitizerV3::ApplyResponseToTriggerDigit(AliMUONDigit& digit,
                                                AliMUONData* data)
{
  if ( !fTriggerEfficiency ) return;

  AliMUONDigit* correspondingDigit = FindCorrespondingDigit(digit,data);

  if(!correspondingDigit)return;//reject bad correspondences

  Int_t detElemId = digit.DetElemId();
  
  const AliMpVSegmentation* segment[2] = 
  {
    Segmentation()->GetMpSegmentation(detElemId,digit.Cathode()), 
    Segmentation()->GetMpSegmentation(detElemId,correspondingDigit->Cathode())
  };

  AliMpPad pad[2] = 
  {
    segment[0]->PadByIndices(AliMpIntPair(digit.PadX(),digit.PadY()),kTRUE), 
    segment[1]->PadByIndices(AliMpIntPair(correspondingDigit->PadX(),correspondingDigit->PadY()),kTRUE)
  };

  Int_t ix(0);
  Int_t iy(1);

  if (digit.Cathode()==0)
  {
	  ix=1;
	  iy=0;
  }
  
  Float_t x = pad[ix].Position().X();
  Float_t y = pad[iy].Position().Y();
  if ( x==-1 && y==-1 )
  {
	  x=-9999.;
	  y=-9999.;
    AliError(Form("Got an unknown position for a digit in DE %d at (ix,iy)=(%d,%d)",
             detElemId,pad[ix].GetIndices().GetFirst(),pad[iy].GetIndices().GetSecond()));
  }
  Float_t x0 = segment[0]->Dimensions().X();
  Float_t y0 = segment[1]->Dimensions().Y();
  TVector2 newCoord = fTriggerEfficiency->ChangeReferenceFrame(x, y, x0, y0);
  
  Bool_t isTrig[2];
  fTriggerEfficiency->IsTriggered(detElemId, newCoord.Px(), newCoord.Py(), 
                                  isTrig[0], isTrig[1]);
  
  if (!isTrig[digit.Cathode()])
  {
	  digit.SetSignal(0);
  }
  
  if ( &digit != correspondingDigit )
  {
	  if (!isTrig[correspondingDigit->Cathode()])
    {
      correspondingDigit->SetSignal(0);
	  }
  }
}

//_____________________________________________________________________________
void
AliMUONDigitizerV3::ApplyResponse()
{
  //
  // Loop over all chamber digits, and apply the response to them
  // Note that this method may remove digits.
  //
  const Bool_t kAddNoise = kTRUE;
  
  for ( Int_t ich = 0; ich < AliMUONConstants::NCh(); ++ich )
  {
    TClonesArray* digits = fOutputData->Digits(ich);
    Int_t n = digits->GetEntriesFast();
    Bool_t trackingChamber = ( ich < AliMUONConstants::NTrackingCh() );
    for ( Int_t i = 0; i < n; ++i )
    {
      AliMUONDigit* d = static_cast<AliMUONDigit*>(digits->UncheckedAt(i));
      if ( trackingChamber )
      {
        ApplyResponseToTrackerDigit(*d,kAddNoise);
      }
      else
      {
        ApplyResponseToTriggerDigit(*d,fOutputData);
      }
      if ( d->Signal() <= 0 )
      {
        digits->RemoveAt(i);
      }
    }
    digits->Compress();
  }    
  
// The version below, using iterator, does not yet work (as the iterator
// assumes it is reading digits from the tree, while in this case it's
// writing...)
//
//  AliMUONDigit* digit(0x0);
//
//  // First loop on tracker digits
//  AliMUONDataIterator tracker(fOutputData,"D",AliMUONDataIterator::kTrackingChambers);
//  
//  while ( ( digit = static_cast<AliMUONDigit*>(tracker.Next()) ) )
//  {
//    ApplyResponseToTrackerDigit(*digit);
//    if ( digit->Signal() <= 0 )
//    {
//      tracker.Remove();
//    }    
//    
//  }
//
//  // Then loop on trigger digits
//  AliMUONDataIterator trigger(fOutputData,"D",AliMUONDataIterator::kTriggerChambers);
//  
//  while ( ( digit = static_cast<AliMUONDigit*>(trigger.Next()) ) )
//  {
//    ApplyResponseToTriggerDigit(*digit,fOutputData);
//    if ( digit->Signal() <= 0 )
//    {
//      trigger.Remove();
//    }    
//  }
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
  // Second, we digitize all the resulting sdigits
  // Then we generate noise-only digits (for tracker only)
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
  
  fExecTimer.Start(kFALSE);

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
  
  if ( fGenerateNoisyDigits )
  {
    // Generate noise-only digits for tracker.
    GenerateNoisyDigits();
  }
  
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
  
  fExecTimer.Stop();
}

//_____________________________________________________________________________
AliMUONDigit* 
AliMUONDigitizerV3::FindCorrespondingDigit(AliMUONDigit& digit,
                                           AliMUONData* data) const
{                                                
  AliMUONDataIterator it(data,"D",AliMUONDataIterator::kTriggerChambers);
  AliMUONDigit* cd;

  for(;;){
      cd = static_cast<AliMUONDigit*>(it.Next());
      if(!cd)continue;
      if (cd->DetElemId() == digit.DetElemId() &&
	  cd->PadX() == digit.PadX() &&
	  cd->PadY() == digit.PadY() && 
	  cd->Cathode() == digit.Cathode()){
	  break;
      }
  }
  //The corresponding digit is searched only forward in the AliMUONData.
  //In this way when the first digit of the couple is given, the second digit is found and the efficiency is applied to both. 
  //Afterwards, when the second digit of the couple is given the first one is not found and hence efficiency is not applied again
  
  while ( ( cd = static_cast<AliMUONDigit*>(it.Next()) ) )
  {
    if(!cd)continue;//avoid problems when 1 digit is removed
    if ( cd->DetElemId() == digit.DetElemId() &&
         cd->Hit() == digit.Hit() &&
         cd->Cathode() != digit.Cathode() )
    {
      return cd;
    }
  }
  return 0x0;
}


//_____________________________________________________________________________
Int_t
AliMUONDigitizerV3::FindDigitIndex(TClonesArray& array, 
                                   const AliMUONDigit& digit) const
{
  // 
  // Return the index of digit within array, if that digit is there, 
  // otherwise returns -1
  //
  // FIXME: this is of course not the best implementation you can think of.
  // Reconsider the use of hit/digit map... ? (but be sure it's needed!)
  //
  
  fFindDigitIndexTimer.Start(kFALSE);
  
  Int_t n = array.GetEntriesFast();
  for ( Int_t i = 0; i < n; ++i )
  {
    AliMUONDigit* d = static_cast<AliMUONDigit*>(array.UncheckedAt(i));
    if ( d->DetElemId() == digit.DetElemId() &&
         d->PadX() == digit.PadX() &&
         d->PadY() == digit.PadY() && 
         d->Cathode() == digit.Cathode() )
    {
      fFindDigitIndexTimer.Stop();
      return i;
    }
  }
  fFindDigitIndexTimer.Stop();
  return -1;
}

//_____________________________________________________________________________
void
AliMUONDigitizerV3::GenerateNoisyDigits()
{
  //
  // According to a given probability, generate digits that
  // have a signal above the noise cut (ped+n*sigma_ped), i.e. digits
  // that are "only noise".
  //
  
  if ( !fNoiseFunction )
  {
    fNoiseFunction = new TF1("AliMUONDigitizerV3::fNoiseFunction","gaus",
                             fgkNSigmas,fgkNSigmas*10);
    
    fNoiseFunction->SetParameters(1,0,1);
  }
  
  fGenerateNoisyDigitsTimer.Start(kFALSE);
  
  for ( Int_t i = 0; i < AliMUONConstants::NTrackingCh(); ++i )
  {
    AliMpDEIterator it;
  
    it.First(i);
  
    while ( !it.IsDone() )
    {
      for ( Int_t cathode = 0; cathode < 2; ++cathode )
      {
        GenerateNoisyDigitsForOneCathode(it.CurrentDE(),cathode);
      }
      it.Next();
    }
  }
  
  fGenerateNoisyDigitsTimer.Stop();
}
 
//_____________________________________________________________________________
void
AliMUONDigitizerV3::GenerateNoisyDigitsForOneCathode(Int_t detElemId, Int_t cathode)
{
  //
  // Generate noise-only digits for one cathode of one detection element.
  // Called by GenerateNoisyDigits()
  //
  
  TClonesArray* digits = fOutputData->Digits(detElemId/100-1);
  
  const AliMpVSegmentation* seg = Segmentation()->GetMpSegmentation(detElemId,cathode);
  Int_t nofPads = seg->NofPads();
  
  Int_t maxIx = seg->MaxPadIndexX();
  Int_t maxIy = seg->MaxPadIndexY();
  
  static const Double_t kProbToBeOutsideNsigmas = TMath::Erfc(fgkNSigmas/TMath::Sqrt(2.0)) / 2. ;
  
  Int_t nofNoisyPads = TMath::Nint(kProbToBeOutsideNsigmas*nofPads);
  if ( !nofNoisyPads ) return;
  
  nofNoisyPads = 
    TMath::Nint(gRandom->Gaus(nofNoisyPads,
                              nofNoisyPads/TMath::Sqrt(nofNoisyPads)));
  
  AliDebug(3,Form("DE %d cath %d nofNoisyPads %d",detElemId,cathode,nofNoisyPads));
  
  for ( Int_t i = 0; i < nofNoisyPads; ++i )
  {
    Int_t ix(-1);
    Int_t iy(-1);
    do {
      ix = gRandom->Integer(maxIx+1);
      iy = gRandom->Integer(maxIy+1);
    } while ( !seg->HasPad(AliMpIntPair(ix,iy)) );
    AliMUONDigit d;
    d.SetDetElemId(detElemId);
    d.SetCathode(cathode);
    d.SetPadX(ix);
    d.SetPadY(iy);
    if ( FindDigitIndex(*digits,d) >= 0 )
    {
      // this digit is already there, and not noise-only, we simply skip it
      continue;
    }
    AliMpPad pad = seg->PadByIndices(AliMpIntPair(ix,iy));
    Int_t manuId = pad.GetLocation().GetFirst();
    Int_t manuChannel = pad.GetLocation().GetSecond();
    
    d.SetElectronics(manuId,manuChannel);
    
    AliMUONVCalibParam* pedestals = fCalibrationData->Pedestals(detElemId,manuId);
    
    Float_t pedestalMean = pedestals->ValueAsFloat(manuChannel,0);
    Float_t pedestalSigma = pedestals->ValueAsFloat(manuChannel,1);
    
    Double_t ped = fNoiseFunction->GetRandom()*pedestalSigma;

    d.SetSignal(TMath::Nint(ped+pedestalMean+0.5));
    d.SetPhysicsSignal(0);
    d.NoiseOnly(kTRUE);
    AliDebug(3,Form("Adding a pure noise digit :"));
    StdoutToAliDebug(3,cout << "Before Response: " << endl; 
                     d.Print(););
    ApplyResponseToTrackerDigit(d,kFALSE);
    if ( d.Signal() > 0 )
    {
      AddOrUpdateDigit(*digits,d);
    }
    else
    {
      AliError("Pure noise below threshold. This should not happen. Not adding "
               "this digit.");
    }
    StdoutToAliDebug(3,cout << "After Response: " << endl; 
                     d.Print(););
  }
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
      fTriggerProcessor = new AliMUONTriggerElectronics(fOutputData,fCalibrationData);
      break;
    default:
      AliFatal("Unknown trigger processor type");
      break;
  }
  AliDebug(1,Form("Using the following trigger code %s - %s",
                  fTriggerProcessor->GetName(),fTriggerProcessor->GetTitle()));
  
  if ( muon()->GetTriggerEffCells() )
  {
    fTriggerEfficiency = fCalibrationData->TriggerEfficiency();
    if ( fTriggerEfficiency )
    {
      AliInfo("Will apply trigger efficiency");
    }
    else
    {
      AliError("I was requested to apply trigger efficiency, but I could "
               "not get it !");
    }
  }
  
  if ( fGenerateNoisyDigits )
  {
    AliInfo("Will generate noise-only digits for tracker");
  }
  
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
