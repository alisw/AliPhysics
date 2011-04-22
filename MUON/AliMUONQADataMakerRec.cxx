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

// --- MUON header files ---
#include "AliMUONQADataMakerRec.h"

//-----------------------------------------------------------------------------
/// \class AliMUONQADataMakerRec
///
/// MUON base class for quality assurance data (histo) maker
///
/// It is now only a steering class for the two subclasses AliMUONTrackerQADataMakerRec
/// and AliMUONTriggerQADataMakerRec
///
/// \author C. Finck, D. Stocco, L. Aphecetche, A. Blanc

#include "AliDAQ.h"
#include "AliMUONTrackerQADataMakerRec.h"
#include "AliMUONTriggerQADataMakerRec.h"
#include "AliQAChecker.h"
#include "AliRawReader.h"
#include "AliRawEventHeaderBase.h"

/// \cond CLASSIMP
ClassImp(AliMUONQADataMakerRec)
/// \endcond
           
//____________________________________________________________________________ 
AliMUONQADataMakerRec::AliMUONQADataMakerRec(Bool_t tracker, Bool_t trigger) : 
AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kMUON), "MUON Quality Assurance Data Maker"),
fTracker(tracker ? new AliMUONTrackerQADataMakerRec(this) : 0x0),
fTrigger(trigger ? new AliMUONTriggerQADataMakerRec(this) : 0x0)
{
  /// ctor
}

//__________________________________________________________________
AliMUONQADataMakerRec::~AliMUONQADataMakerRec()
{
    /// dtor
  delete fTracker;
  delete fTrigger;
}

//____________________________________________________________________________ 
Int_t AliMUONQADataMakerRec::Add2List(TH1 * hist, const Int_t index, AliQAv1::TASKINDEX_t task, const Bool_t expert, const Bool_t image, const Bool_t saveForCorr)
{
  TObjArray** list = GetList(task);
  if (list)
  {
    return Add2List(hist,index,list,expert,image,saveForCorr);
  }
  return -1;
}

//____________________________________________________________________________ 
void AliMUONQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray** list)
{
  /// Detector specific actions at end of cycle
  
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
  {
    if (! IsValidEventSpecie(specie, list)  ) continue;
    
    SetEventSpecie(AliRecoParam::ConvertIndex(specie)); // needed by the GetXXXData methods
    
    if ( task == AliQAv1::kRAWS ) 
    {
      if ( fTracker ) fTracker->EndOfDetectorCycleRaws(specie,list);
      if ( fTrigger ) fTrigger->EndOfDetectorCycleRaws(specie,list);
    }  
    else if ( task == AliQAv1::kRECPOINTS )
    {
      // normalize recpoints histograms
      if ( fTracker ) fTracker->EndOfDetectorCycleRecPoints(specie,list);
      if ( fTrigger ) fTrigger->EndOfDetectorCycleRecPoints(specie,list);
    }
    else if ( task == AliQAv1::kESDS ) 
    {
      // normalize esds histograms
      if ( fTracker ) fTracker->EndOfDetectorCycleESDs(specie,list);
      if ( fTrigger ) fTrigger->EndOfDetectorCycleESDs(specie,list);
    }
    else if ( task == AliQAv1::kDIGITSR ) 
    {
      if ( fTracker ) fTracker->EndOfDetectorCycleDigits(specie,list);        
      if ( fTrigger ) fTrigger->EndOfDetectorCycleDigits(specie,list);
    }
    else
    {
      AliFatal(Form("Not implemented for task %s",AliQAv1::GetTaskName(task).Data()));
    }
  } // loop on specie
    
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQAv1::kMUON,task,list,const_cast<AliDetectorRecoParam*>(GetRecoParam()));
}

//____________________________________________________________________________ 
TObject* AliMUONQADataMakerRec::GetData(AliQAv1::TASKINDEX_t task, const Int_t index)
{
  TObjArray** list = GetList(task);
  if (list) return GetData(list,index);
  return 0x0;
}

//____________________________________________________________________________ 
TObjArray** AliMUONQADataMakerRec::GetList(AliQAv1::TASKINDEX_t task)
{
  //  enum TASKINDEX_t {
  //    kNULLTASKINDEX=-1, kRAWS, kHITS, kSDIGITS, kDIGITS, kDIGITSR, kRECPOINTS, kTRACKSEGMENTS, kRECPARTICLES, kESDS, kNTASKINDEX };
  if ( task == AliQAv1::kRAWS ) 
  {
      return fRawsQAList;
  }
  else if ( task == AliQAv1::kDIGITS || task == AliQAv1::kDIGITSR )
  {
    return fDigitsQAList;
  }
  else if ( task == AliQAv1::kRECPOINTS ) 
  {
    return fRecPointsQAList;
  }
  else
  {
      AliFatal(Form("task %s not supported here yet",AliQAv1::GetTaskName(task).Data()));
  }
  return 0x0;
}

//____________________________________________________________________________ 
void AliMUONQADataMakerRec::InitRaws()
{
  /// create Raws histograms in Raws subdir
	
  if ( fTracker ) fTracker->InitRaws();
  if ( fTrigger ) fTrigger->InitRaws();
}

//__________________________________________________________________
void AliMUONQADataMakerRec::InitDigits() 
{
  /// Initialized Digits spectra 
  if ( fTracker ) fTracker->InitDigits();
  if ( fTrigger ) fTrigger->InitDigits();
} 

//____________________________________________________________________________ 
void AliMUONQADataMakerRec::InitRecPoints()
{
	/// create Reconstructed Points histograms in RecPoints subdir
  if ( fTracker ) fTracker->InitRecPoints();
  if ( fTrigger ) fTrigger->InitRecPoints();
}


//____________________________________________________________________________ 
void AliMUONQADataMakerRec::InitESDs()
{
  ///create ESDs histograms in ESDs subdir
  if ( fTracker ) fTracker->InitESDs();
  if ( fTrigger ) fTrigger->InitESDs();
}

//____________________________________________________________________________
void AliMUONQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
  /// make QA for rawdata
  /// Note that we do not call the sub-datamaker MakeRaws method
  /// for events where the MCH or MTR is not part of the readout...
  
  if ( !rawReader || !rawReader->GetDetectorPattern() ) return;
  
  UInt_t clmask = rawReader->GetDetectorPattern()[0];
    
  if ( fTracker && rawReader->GetType() == AliRawEventHeaderBase::kPhysicsEvent ) 
  {
    UInt_t mchMask = AliDAQ::DetectorPattern(" MUONTRK ");
    if ( clmask & mchMask ) 
    {
      rawReader->Reset();
      fTracker->MakeRaws(rawReader);
    }
  }
  
  if ( fTrigger && (rawReader->GetType() == AliRawEventHeaderBase::kPhysicsEvent ||
                    rawReader->GetType() == AliRawEventHeaderBase::kCalibrationEvent ) )
  {
    UInt_t mtrMask = AliDAQ::DetectorPattern(" MUONTRG ");
    if ( clmask & mtrMask )
    {
      rawReader->Reset();    
      fTrigger->MakeRaws(rawReader);
    }
  }
}

//__________________________________________________________________
void AliMUONQADataMakerRec::MakeDigits()         
{
  /// makes data from Digits
  
  AliFatal("Not implemented");
}

//__________________________________________________________________
void AliMUONQADataMakerRec::MakeDigits(TTree* digitsTree)         
{
  /// makes data from Digits

  // Do nothing in case of calibration event
  if ( GetEventSpecie() == AliRecoParam::kCalib ) return;

  if ( fTracker ) fTracker->MakeDigits(digitsTree);
  if ( fTrigger ) fTrigger->MakeDigits(digitsTree);  
}

//____________________________________________________________________________
void AliMUONQADataMakerRec::MakeRecPoints(TTree* clustersTree)
{
	/// Fill histograms from treeR

  // Do nothing in case of calibration event
  if ( GetEventSpecie() == AliRecoParam::kCalib ) return;
	
  if ( fTracker ) fTracker->MakeRecPoints(clustersTree);
  if ( fTrigger ) fTrigger->MakeRecPoints(clustersTree);  
}

//____________________________________________________________________________
void AliMUONQADataMakerRec::MakeESDs(AliESDEvent* esd)
{
  /// make QA data from ESDs

  // Do nothing in case of calibration event
  if ( GetEventSpecie() == AliRecoParam::kCalib ) return;
  
  if ( fTracker ) fTracker->MakeESDs(esd);
  if ( fTrigger ) fTrigger->MakeESDs(esd);  

 }

//____________________________________________________________________________ 
void AliMUONQADataMakerRec::ResetDetector(AliQAv1::TASKINDEX_t task)
{
  /// Reset internals
  
  for (int spec = 0; spec < AliRecoParam::kNSpecies; ++spec) 
  {
    if (!AliQAv1::Instance()->IsEventSpecieSet(AliRecoParam::ConvertIndex(spec)))
      continue;
    
    if ( task == AliQAv1::kRAWS ) 
    {
      if (fTracker) fTracker->ResetDetectorRaws(fRawsQAList[spec]);
      if (fTrigger) fTrigger->ResetDetectorRaws(fRawsQAList[spec]);
    }
    else if ( task == AliQAv1::kRECPOINTS )
    {
      if (fTracker) fTracker->ResetDetectorRecPoints(fRecPointsQAList[spec]);
      if (fTrigger) fTrigger->ResetDetectorRecPoints(fRecPointsQAList[spec]);
    }
    else if ( task == AliQAv1::kESDS ) 
    {
      if (fTracker) fTracker->ResetDetectorESDs(fESDsQAList[spec]);
      if (fTrigger) fTrigger->ResetDetectorESDs(fESDsQAList[spec]);
    }
    else if ( task == AliQAv1::kDIGITS ) 
    {
      if (fTracker) fTracker->ResetDetectorDigits(fDigitsQAList[spec]);
      if (fTrigger) fTrigger->ResetDetectorDigits(fDigitsQAList[spec]);
    }
    else
    {
      AliFatal(Form("Not implemented for task %s",AliQAv1::GetTaskName(task).Data()));
    }
  }
}

//____________________________________________________________________________ 
void AliMUONQADataMakerRec::StartOfDetectorCycle()
{
    /// Detector specific actions at start of cycle  
}
