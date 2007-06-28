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

#include "AliCDBManager.h"
#include "AliLog.h"
#include "AliMUON.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONConstants.h"
#include "AliMUONDigit.h"
#include "AliMUONLogger.h"
#include "AliMUONTriggerEfficiencyCells.h"
#include "AliMUONTriggerElectronics.h"
#include "AliMUONTriggerStoreV1.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVDigitStore.h"
#include "AliMpCathodType.h"
#include "AliMpConstants.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpDEManager.h"
#include "AliMpIntPair.h"
#include "AliMpPad.h"
#include "AliMpSegmentation.h"
#include "AliMpStationType.h"
#include "AliMpVSegmentation.h"
#include "AliRun.h"
#include "AliRunDigitizer.h"
#include "AliRunLoader.h"
#include <Riostream.h>
#include <TF1.h>
#include <TFile.h>
#include <TMath.h>
#include <TRandom.h>
#include <TString.h>
#include <TSystem.h>

#include "AliMUONGeometryTransformer.h" //ADDED for trigger noise
///
/// \class AliMUONDigitizerV3
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
/// \author Laurent Aphecetche

namespace
{
  AliMUON* muon()
  {
    return static_cast<AliMUON*>(gAlice->GetModule("MUON"));
  }

  //ADDED for trigger noise
  const AliMUONGeometryTransformer* GetTransformer()
  {
      return muon()->GetGeometryTransformer();
  }
}

const Double_t AliMUONDigitizerV3::fgkNSigmas=3;

/// \cond CLASSIMP
ClassImp(AliMUONDigitizerV3)
/// \endcond

//_____________________________________________________________________________
AliMUONDigitizerV3::AliMUONDigitizerV3(AliRunDigitizer* manager, 
                                       Int_t generateNoisyDigits)
: AliDigitizer(manager),
fIsInitialized(kFALSE),
fCalibrationData(0x0),
fTriggerProcessor(0x0),
fTriggerEfficiency(0x0),
fGenerateNoisyDigitsTimer(),
fExecTimer(),
fNoiseFunction(0x0),
fNoiseFunctionTrig(0x0),
  fGenerateNoisyDigits(generateNoisyDigits),
  fLogger(new AliMUONLogger(1000)),
fTriggerStore(new AliMUONTriggerStoreV1),
fDigitStore(0x0),
fOutputDigitStore(0x0)
{
  /// Ctor.

  AliDebug(1,Form("AliRunDigitizer=%p",fManager));
  fGenerateNoisyDigitsTimer.Start(kTRUE); fGenerateNoisyDigitsTimer.Stop();
  fExecTimer.Start(kTRUE); fExecTimer.Stop();
}

//_____________________________________________________________________________
AliMUONDigitizerV3::~AliMUONDigitizerV3()
{
  /// Dtor. Note we're the owner of some pointers.

  AliDebug(1,"dtor");

  delete fCalibrationData;
  delete fTriggerProcessor;
  delete fNoiseFunction;
  delete fNoiseFunctionTrig;
  delete fTriggerStore;
  delete fDigitStore;
  delete fOutputDigitStore;
  
  if ( fGenerateNoisyDigits )
  {
    AliDebug(1, Form("Execution time for GenerateNoisyDigits() : R:%.2fs C:%.2fs",
                 fGenerateNoisyDigitsTimer.RealTime(),
                 fGenerateNoisyDigitsTimer.CpuTime()));
  }
  AliDebug(1, Form("Execution time for Exec() : R:%.2fs C:%.2fs",
               fExecTimer.RealTime(),fExecTimer.CpuTime()));
 
  AliInfo("Summary of messages");
  fLogger->Print();
  
  delete fLogger;
}

//_____________________________________________________________________________
void 
AliMUONDigitizerV3::ApplyResponseToTrackerDigit(AliMUONVDigit& digit, Bool_t addNoise)
{
  /// For tracking digits, starting from an ideal digit's charge, we :
  ///
  /// - add some noise (thus leading to a realistic charge), if requested to do so
  /// - "divide" by a gain (thus decalibrating the digit)
  /// - add a pedestal (thus decalibrating the digit)
  /// - sets the signal to zero if below 3*sigma of the noise

  static const Int_t kMaxADC = (1<<12)-1; // We code the charge on a 12 bits ADC.
  
  Float_t charge = digit.Charge();
  
  // We set the charge to 0, as the only relevant piece of information
  // after Digitization is the ADC value.  
  digit.SetCharge(0);
  
  if ( !addNoise )
  {
    digit.SetADC(TMath::Min(kMaxADC,TMath::Nint(charge)));
    return;
  }
  
  Int_t detElemId = digit.DetElemId();
  
  Int_t manuId = digit.ManuId();
  Int_t manuChannel = digit.ManuChannel();
  
  AliMUONVCalibParam* pedestal = fCalibrationData->Pedestals(detElemId,manuId);
  if (!pedestal)
  {
    fLogger->Log(Form("%s:%d:Could not get pedestal for DE=%4d manuId=%4d. Disabling.",
                      __FILE__,__LINE__,
                      detElemId,manuId));
    digit.SetADC(0);
    return;    
  }
  Float_t pedestalMean = pedestal->ValueAsFloat(manuChannel,0);
  Float_t pedestalSigma = pedestal->ValueAsFloat(manuChannel,1);
  
  AliMUONVCalibParam* gain = fCalibrationData->Gains(detElemId,manuId);
  if (!gain)
  {
    fLogger->Log(Form("%s:%d:Could not get gain for DE=%4d manuId=%4d. Disabling.",
                      __FILE__,__LINE__,
                      detElemId,manuId));
    digit.SetADC(0);
    return;        
  }    

  Float_t a0 = gain->ValueAsFloat(manuChannel,0);
  Float_t a1 = gain->ValueAsFloat(manuChannel,1);

  Int_t thres = gain->ValueAsInt(manuChannel,2);
  Float_t chargeThres = a0*thres;
  
  Float_t padc(0); // (adc - ped) value
  
  if ( charge <= chargeThres || TMath::Abs(a1) < 1E-12 ) 
  {
    // linear part only
    
    if ( TMath::Abs(a0) > 1E-12 ) 
    {
      padc = charge/a0;    
    }
  }
  else // charge > chargeThres && a1 not zero
  {
    // parabolic part
    padc = TMath::Sqrt((chargeThres-charge)/a1) + thres;
  }

  Float_t adcNoise = gRandom->Gaus(0.0,pedestalSigma);
  Int_t adc(0);
  
  if ( padc > 0 ) 
  {
    adc = TMath::Nint(padc + pedestalMean + adcNoise);
  }
  
  // be sure we stick to 12 bits.
  if ( adc > kMaxADC )
  {
    adc = kMaxADC;
  }
  
  AliDebug(3,Form("DE %4d Manu %4d Ch %2d Charge %e A0 %e A1 %e Thres %d padc %e ADC %4d",
                  detElemId,manuId,manuChannel,charge,
                  gain->ValueAsFloat(manuChannel,0),
                  gain->ValueAsFloat(manuChannel,1),
                  gain->ValueAsInt(manuChannel,2),
                  padc,
                  adc));
                  
  digit.SetADC(adc);
}

//_____________________________________________________________________________
void 
AliMUONDigitizerV3::ApplyResponseToTriggerDigit(const AliMUONVDigitStore& digitStore,
                                                AliMUONVDigit& digit)
{
  /// \todo add comment

  if ( !fTriggerEfficiency ) return;

  if (digit.IsEfficiencyApplied()) return;

  AliMUONVDigit* correspondingDigit = FindCorrespondingDigit(digitStore,digit);

  if (!correspondingDigit) return; //reject bad correspondences

  Int_t detElemId = digit.DetElemId();

  AliMpSegmentation* segmentation = AliMpSegmentation::Instance();
  const AliMpVSegmentation* segment[2] = 
  {
    segmentation->GetMpSegmentation(detElemId,AliMp::GetCathodType(digit.Cathode())), 
    segmentation->GetMpSegmentation(detElemId,AliMp::GetCathodType(correspondingDigit->Cathode()))
  };

  AliMpPad pad[2] = 
  {
    segment[0]->PadByIndices(AliMpIntPair(digit.PadX(),digit.PadY()),kTRUE), 
    segment[1]->PadByIndices(AliMpIntPair(correspondingDigit->PadX(),correspondingDigit->PadY()),kTRUE)
  };

  Int_t p0(1);
  if (digit.Cathode()==0) p0=0;

  AliMpIntPair location = pad[p0].GetLocation(0);
  Int_t nboard = location.GetFirst();

  Bool_t isTrig[2];

  fTriggerEfficiency->IsTriggered(detElemId, nboard, 
                                  isTrig[0], isTrig[1]);
  digit.EfficiencyApplied(kTRUE);
  correspondingDigit->EfficiencyApplied(kTRUE);

  if (!isTrig[digit.Cathode()])
  {
	  digit.SetCharge(0);
  }
  
  if ( &digit != correspondingDigit )
  {
	  if (!isTrig[correspondingDigit->Cathode()])
    {
      correspondingDigit->SetCharge(0);
	  }
  }
}

//_____________________________________________________________________________
void
AliMUONDigitizerV3::ApplyResponse(const AliMUONVDigitStore& store,
                                  AliMUONVDigitStore& filteredStore)
{
  /// Loop over all chamber digits, and apply the response to them
  /// Note that this method may remove digits.

  filteredStore.Clear();
  
  const Bool_t kAddNoise = kTRUE;
  
  TIter next(store.CreateIterator());
  AliMUONVDigit* digit;
  
  while ( ( digit = static_cast<AliMUONVDigit*>(next()) ) )
  {
    AliMp::StationType stationType = AliMpDEManager::GetStationType(digit->DetElemId());
    
    if ( stationType != AliMp::kStationTrigger )
    {
      ApplyResponseToTrackerDigit(*digit,kAddNoise);
    }
    else
    {
      ApplyResponseToTriggerDigit(store,*digit);
    }
    if ( digit->ADC() > 0  || digit->Charge() > 0 )
    {
      filteredStore.Add(*digit,AliMUONVDigitStore::kIgnore);
    }
  }
}    

//_____________________________________________________________________________
void
AliMUONDigitizerV3::Exec(Option_t*)
{
  /// Main method.
  /// We first loop over input files, and merge the sdigits we found there.
  /// Second, we digitize all the resulting sdigits
  /// Then we generate noise-only digits (for tracker only)
  /// And we finally generate the trigger outputs.
    
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
  
  AliLoader* outputLoader = GetLoader(fManager->GetOutputFolderName());
  
  outputLoader->MakeDigitsContainer();
  
  TTree* oTreeD = outputLoader->TreeD();
  
  if (!oTreeD) 
  {
    AliFatal("Cannot create output TreeD");
  }

  // Loop over all the input files, and merge the sdigits found in those
  // files.
  
  for ( Int_t iFile = 0; iFile < nInputFiles; ++iFile )
  {    
    AliLoader* inputLoader = GetLoader(fManager->GetInputFolderName(iFile));

    inputLoader->LoadSDigits("READ");

    TTree* iTreeS = inputLoader->TreeS();
    if (!iTreeS)
    {
      AliFatal(Form("Could not get access to input file #%d",iFile));
    }
    
    AliMUONVDigitStore* inputStore = AliMUONVDigitStore::Create(*iTreeS);
    inputStore->Connect(*iTreeS);
    
    iTreeS->GetEvent(0);
    
    MergeWithSDigits(fDigitStore,*inputStore,fManager->GetMask(iFile));

    inputLoader->UnloadSDigits();
    
    delete inputStore;
  }
  
  // At this point, we do have digit arrays (one per chamber) which contains 
  // the merging of all the sdigits of the input file(s).
  // We now massage them to apply the detector response, i.e. this
  // is here that we do the "digitization" work.
  
  if (!fOutputDigitStore)
  {
    fOutputDigitStore = fDigitStore->Create();
  }
  
  if ( fGenerateNoisyDigits>=2 )
  {
    // Generate noise-only digits for trigger.
    GenerateNoisyDigitsForTrigger(*fDigitStore);
  }

  ApplyResponse(*fDigitStore,*fOutputDigitStore);
  
  if ( fGenerateNoisyDigits )
  {
    // Generate noise-only digits for tracker.
    GenerateNoisyDigits(*fOutputDigitStore);
  }
  
  // We generate the global and local trigger decisions.
  fTriggerProcessor->Digits2Trigger(*fOutputDigitStore,*fTriggerStore);

  // Prepare output tree
  Bool_t okD = fOutputDigitStore->Connect(*oTreeD,kFALSE);
  Bool_t okT = fTriggerStore->Connect(*oTreeD,kFALSE);
  if (!okD || !okT)
  {
    AliError(Form("Could not make branch : Digit %d Trigger %d",okD,okT));
    return;
  }
  
  // Fill the output treeD
  oTreeD->Fill();
  
  // Write to the output tree(D).
  // Please note that as GlobalTrigger, LocalTrigger and Digits are in the same
  // tree (=TreeD) in different branches, this WriteDigits in fact writes all of 
  // the 3 branches.
  outputLoader->WriteDigits("OVERWRITE");  
  
  outputLoader->UnloadDigits();
  
  // Finally, we clean up after ourselves.
  fTriggerStore->Clear();
  fDigitStore->Clear();
  fOutputDigitStore->Clear();
  fExecTimer.Stop();
}

//_____________________________________________________________________________
AliMUONVDigit* 
AliMUONDigitizerV3::FindCorrespondingDigit(const AliMUONVDigitStore& digitStore,
                                           AliMUONVDigit& digit) const
{                                                
  /// Find, if it exists, the digit corresponding to digit.Hit(), in the 
  /// other cathode

  TIter next(digitStore.CreateIterator());
  AliMUONVDigit* d;
  
  while ( ( d = static_cast<AliMUONVDigit*>(next()) ) )
  {
    if ( d->DetElemId() == digit.DetElemId() &&
         d->Hit() == digit.Hit() &&
         d->Cathode() != digit.Cathode() )
    {
      return d;
    }      
  }    
  return 0x0;
}


//_____________________________________________________________________________
void
AliMUONDigitizerV3::GenerateNoisyDigits(AliMUONVDigitStore& digitStore)
{
  /// According to a given probability, generate digits that
  /// have a signal above the noise cut (ped+n*sigma_ped), i.e. digits
  /// that are "only noise".
  
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
        GenerateNoisyDigitsForOneCathode(digitStore,it.CurrentDEId(),cathode);
      }
      it.Next();
    }
  }
  
  fGenerateNoisyDigitsTimer.Stop();
}
 
//_____________________________________________________________________________
void
AliMUONDigitizerV3::GenerateNoisyDigitsForOneCathode(AliMUONVDigitStore& digitStore,
                                                     Int_t detElemId, Int_t cathode)
{
  /// Generate noise-only digits for one cathode of one detection element.
  /// Called by GenerateNoisyDigits()
  
  const AliMpVSegmentation* seg 
    = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,AliMp::GetCathodType(cathode));
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
    AliMpPad pad;
    
    do {
      ix = gRandom->Integer(maxIx+1);
      iy = gRandom->Integer(maxIy+1);
      pad = seg->PadByIndices(AliMpIntPair(ix,iy),kFALSE);
    } while ( !pad.IsValid() );

    Int_t manuId = pad.GetLocation().GetFirst();
    Int_t manuChannel = pad.GetLocation().GetSecond();    

    AliMUONVCalibParam* pedestals = fCalibrationData->Pedestals(detElemId,manuId);
    
    if (!pedestals) 
    {
      // no pedestal available for this channel, simply give up
      continue;
    }
    
    AliMUONVDigit* d = digitStore.CreateDigit(detElemId,manuId,manuChannel,cathode);
    
    d->SetPadXY(ix,iy);
    
    Float_t pedestalMean = pedestals->ValueAsFloat(manuChannel,0);
    Float_t pedestalSigma = pedestals->ValueAsFloat(manuChannel,1);
    
    Double_t ped = fNoiseFunction->GetRandom()*pedestalSigma;

    d->SetCharge(TMath::Nint(ped+pedestalMean+0.5));
    d->NoiseOnly(kTRUE);
    ApplyResponseToTrackerDigit(*d,kFALSE);
    if ( d->ADC() > 0 )
    {
      Bool_t ok = digitStore.Add(*d,AliMUONVDigitStore::kDeny);
      // this can happen (that we randomly chose a digit that is
      // already there). We simply ignore this, but log the occurence
      // to cross-check that it's not too frequent.
      if (!ok)
      {
        fLogger->Log("Collision while adding noiseOnly digit");
      }
      else
      {
        fLogger->Log("Added noiseOnly digit");
      }
    }
    else
    {
      AliError("Pure noise below threshold. This should not happen. Not adding "
               "this digit.");
    }
    delete d;
  }
}


//_____________________________________________________________________________
void
AliMUONDigitizerV3::GenerateNoisyDigitsForTrigger(AliMUONVDigitStore& digitStore)
{
  /// Generate noise-only digits for one cathode of one detection element.
  /// Called by GenerateNoisyDigits()

  if ( !fNoiseFunctionTrig )
  {
    fNoiseFunctionTrig = new TF1("AliMUONDigitizerV3::fNoiseFunctionTrig","landau",
				 50.,270.);
    
    fNoiseFunctionTrig->SetParameters(3.91070e+02, 9.85026, 9.35881e-02);
  }

  AliMpPad pad[2];
  AliMUONVDigit *d[2]={0x0};

  for ( Int_t chamberId = AliMUONConstants::NTrackingCh(); chamberId < AliMUONConstants::NCh(); ++chamberId )
  {
  
    Int_t nofNoisyPads = 50;

    Float_t r=-1, fi = 0., gx, gy, x, y, z, xg01, yg01, zg, xg02, yg02;
    AliMpDEIterator it;
  
    AliDebug(3,Form("Chamber %d nofNoisyPads %d",chamberId,nofNoisyPads));

    for ( Int_t i = 0; i < nofNoisyPads; ++i )
    {
      //printf("Generating noise %i\n",i);
	Int_t ix(-1);
	Int_t iy(-1);
	Bool_t isOk = kFALSE;
	Int_t detElemId = -1;
	do {
	  //r = gRandom->Landau(9.85026, 9.35881e-02);
	    r = fNoiseFunctionTrig->GetRandom();
	    fi = 2. * TMath::Pi() * gRandom->Rndm();
	    //printf("r = %f\tfi = %f\n", r, fi);
	    gx = r * TMath::Cos(fi);
	    gy = r * TMath::Sin(fi);

	    for ( it.First(chamberId); ! it.IsDone(); it.Next() ){
		Int_t currDetElemId = it.CurrentDEId();
		const AliMpVSegmentation* seg
		    = AliMpSegmentation::Instance()->GetMpSegmentation(currDetElemId,AliMp::GetCathodType(0));
		if (!seg) continue;
		Float_t deltax = seg->Dimensions().X();
		Float_t deltay = seg->Dimensions().Y();
		GetTransformer()->Local2Global(currDetElemId, -deltax, -deltay, 0, xg01, yg01, zg);
		GetTransformer()->Local2Global(currDetElemId,  deltax,  deltay, 0, xg02, yg02, zg);
		Float_t xg1 = xg01, xg2 = xg02, yg1 = yg01, yg2 = yg02;
		if(xg01>xg02){
		    xg1 = xg02;
		    xg2 = xg01;
		}
		if(yg01>yg02){
		    yg1 = yg02;
		    yg2 = yg01;
		}
		if(gx>=xg1 && gx<=xg2 && gy>=yg1 && gy<=yg2){
		    detElemId = currDetElemId;
		    GetTransformer()->Global2Local(detElemId, gx, gy, 0, x, y, z);
		    pad[0] = seg->PadByPosition(TVector2(x,y),kFALSE);
		    if(!pad[0].IsValid()) continue;
		    isOk = kTRUE;
		    break;
		}
	    } // loop on slats
	} while ( !isOk );

	const AliMpVSegmentation* seg1
	    = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,AliMp::GetCathodType(1));
	pad[1] = seg1->PadByPosition(TVector2(x,y),kFALSE);

	for ( Int_t cathode = 0; cathode < 2; ++cathode ){
	  Int_t manuId = pad[cathode].GetLocation(0).GetFirst();
	  Int_t manuChannel = pad[cathode].GetLocation(0).GetSecond();    
	  d[cathode] = digitStore.CreateDigit(detElemId,manuId,manuChannel,cathode);
	  ix = pad[cathode].GetIndices().GetFirst();
	  iy = pad[cathode].GetIndices().GetSecond();
	  d[cathode]->SetPadXY(ix,iy);
	  //d[cathode].SetSignal(1);
	  //d[cathode].SetPhysicsSignal(0);
	  d[cathode]->SetCharge(1);
	  d[cathode]->NoiseOnly(kTRUE);
	  AliDebug(3,Form("Adding a pure noise digit :"));

	  Bool_t ok = digitStore.Add(*d[cathode],AliMUONVDigitStore::kDeny);
	  if (!ok)
	  {
	      fLogger->Log("Collision while adding TriggerNoise digit");
	  }
	  else
	  {
	      fLogger->Log("Added triggerNoise digit");
	  }
	} //loop on cathodes
    } // loop on noisy pads
  } // loop on chambers
}


//_____________________________________________________________________________
AliLoader*
AliMUONDigitizerV3::GetLoader(const TString& folderName)
{
  /// Get a MUON loader

  AliDebug(2,Form("Getting access to folder %s",folderName.Data()));
  AliLoader* loader = AliRunLoader::GetDetectorLoader("MUON",folderName.Data());
  if (!loader)
  {
    AliError(Form("Could not get MuonLoader from folder %s",folderName.Data()));
    return 0x0;
  }
  return loader;
}

//_____________________________________________________________________________
Bool_t
AliMUONDigitizerV3::Init()
{
  /// Initialization of the TTask :
  /// a) create the calibrationData, according to run number
  /// b) create the trigger processing task

  AliDebug(2,"");
  
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
  
  Int_t runnumber = AliCDBManager::Instance()->GetRun();
  
  fCalibrationData = new AliMUONCalibrationData(runnumber);
  if ( !fCalibrationData->Pedestals() )
  {
    AliFatal("Could not access pedestals from OCDB !");
  }
  if ( !fCalibrationData->Gains() )
  {
    AliFatal("Could not access gains from OCDB !");
  }
  fTriggerProcessor = new AliMUONTriggerElectronics(fCalibrationData);
  
  if ( muon()->GetTriggerEffCells() )
  {
    fTriggerEfficiency = fCalibrationData->TriggerEfficiency();
    if ( fTriggerEfficiency )
    {
      AliDebug(1, "Will apply trigger efficiency");
    }
    else
    {
      AliFatal("I was requested to apply trigger efficiency, but I could "
               "not get it !");
    }
  }
  
  AliDebug(1, Form("Will %s generate noise-only digits for tracker",
                     (fGenerateNoisyDigits ? "":"NOT")));

  fIsInitialized = kTRUE;
  return kTRUE;
}

//_____________________________________________________________________________
void 
AliMUONDigitizerV3::MergeWithSDigits(AliMUONVDigitStore*& outputStore,
                                     const AliMUONVDigitStore& input,
                                     Int_t mask)
{
  /// Merge the sdigits in inputData with the digits already present in outputData
  
  if ( !outputStore ) outputStore = input.Create();
  
  TIter next(input.CreateIterator());
  AliMUONVDigit* sdigit;
  
  while ( ( sdigit = static_cast<AliMUONVDigit*>(next()) ) )
  {
    // Update the track references using the mask.
    // FIXME: this is dirty, for backward compatibility only.
    // Should re-design all this way of keeping track of MC information...
    if ( mask ) sdigit->PatchTracks(mask);
    // Then add or update the digit to the output.
    AliMUONVDigit* added = outputStore->Add(*sdigit,AliMUONVDigitStore::kMerge);
    if (!added)
    {
      AliError("Could not add digit in merge mode");
    }
  }
}
