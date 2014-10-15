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

#include "AliMUONQAChecker.h"

/// \class AliMUONQAChecker
///
/// Implementation of AliQACheckerBase for MCH and MTR
///
/// \author Laurent Aphecetche, Subatech

#include "AliMUONRecoParam.h"
#include "AliMUONTrackerQAChecker.h"
#include "AliMUONTriggerQAChecker.h"
#include "AliCodeTimer.h"
#include "AliMUONQAIndices.h"

/// \cond CLASSIMP
ClassImp(AliMUONQAChecker)
/// \endcond

namespace
{
  const Int_t TRACKER=0;
  const Int_t TRIGGER=1;
}

//__________________________________________________________________
AliMUONQAChecker::AliMUONQAChecker() : 
    AliQACheckerBase("MUON","MUON Quality Assurance Data Maker"),
fCheckers(new TObjArray)
{
	/// ctor
  fCheckers->SetOwner(kTRUE);
  fCheckers->AddAt(new AliMUONTrackerQAChecker(),TRACKER);
  fCheckers->AddAt(new AliMUONTriggerQAChecker(),TRIGGER);
}          

//__________________________________________________________________
AliMUONQAChecker::~AliMUONQAChecker() 
{
	/// dtor
  delete fCheckers;
}

//______________________________________________________________________________
void
AliMUONQAChecker::Check(Double_t* rv, AliQAv1::ALITASK_t index, 
                        TObjArray** list, 
                        const AliDetectorRecoParam * recoParam)
{
  /// Check objects in list
  
  AliCodeTimerAuto(AliQAv1::GetTaskName(index),0);
  
  const AliMUONRecoParam* muonRecoParam = static_cast<const AliMUONRecoParam*>(recoParam);
  AliMUONVQAChecker::ECheckCode* ecc(0x0);

  for ( Int_t i = 0; i < AliRecoParam::kNSpecies; ++i ) 
  {
    rv[i] = -1.0;
  }
  
  for ( Int_t ic = 0; ic <= fCheckers->GetLast(); ++ic )
  {
    if ( ic != TRACKER && ic != TRIGGER ) continue;

    Bool_t trackerRequested(kFALSE);
    Bool_t triggerRequested(kFALSE);
    
    for ( Int_t i = 0; i < AliRecoParam::kNSpecies; ++i ) 
    {
      // no need to take into account detector that was not requested
      if ( ic == TRACKER && AliQAv1::GetData(list,AliMUONQAIndices::kTrackerIsThere,AliRecoParam::ConvertIndex(i)) ) trackerRequested=kTRUE;
      if ( ic == TRIGGER && AliQAv1::GetData(list,AliMUONQAIndices::kTriggerIsThere,AliRecoParam::ConvertIndex(i)) ) triggerRequested=kTRUE;
    }
    
    if ( ic == TRACKER && !trackerRequested ) 
    {
      AliInfo("Skipping tracker check as tracker not requested");
      continue;      
    }
    
    if ( ic == TRIGGER && !triggerRequested ) 
    {
      AliInfo("Skipping trigger check as trigger not requested");
      continue;      
    }
    
    AliMUONVQAChecker* qac = static_cast<AliMUONVQAChecker*>(fCheckers->At(ic));
    
    if ( index == AliQAv1::kRAW ) 
    {    
      ecc = qac->CheckRaws(list,muonRecoParam);
    }

    if ( index == AliQAv1::kREC)
    {
      ecc = qac->CheckRecPoints(list,muonRecoParam);
    }
    
    if ( index == AliQAv1::kESD )
    {
      ecc = qac->CheckESD(list,muonRecoParam);
    }
    
    if ( ecc ) 
    {
      for ( Int_t i = 0; i < AliRecoParam::kNSpecies; ++i ) 
      {
        // no need to take into account detector that was not requested
        if ( ic == TRACKER && AliQAv1::GetData(list,AliMUONQAIndices::kTrackerIsThere,AliRecoParam::ConvertIndex(i))==0x0 ) continue;
        if ( ic == TRIGGER && AliQAv1::GetData(list,AliMUONQAIndices::kTriggerIsThere,AliRecoParam::ConvertIndex(i))==0x0 ) continue;
                
        switch ( ecc[i] ) 
        {
          case AliMUONVQAChecker::kInfo:
            rv[i] = 1.0;
            break;
          case AliMUONVQAChecker::kWarning:
            rv[i] = 0.75;
            break;
          case AliMUONVQAChecker::kError:
            rv[i] = 0.25;
            break;
          case AliMUONVQAChecker::kFatal:
            rv[i] = -1.0;
            break;
          default:
            AliError("Invalid ecc value. FIXME !");
            rv[i] = -1.0;
            break;
        }
      }
    }

    delete[] ecc;
  }
}

//______________________________________________________________________________
void AliMUONQAChecker::Init(const AliQAv1::DETECTORINDEX_t det) 
{
  /// intialises QA and QA checker settings
  AliQAv1::Instance(det) ; 
  Float_t hiValue[AliQAv1::kNBIT] ; 
  Float_t lowValue[AliQAv1::kNBIT] ;
  lowValue[AliQAv1::kINFO]      = 0.999   ; 
  hiValue[AliQAv1::kINFO]       = 1.0 ; 
  hiValue[AliQAv1::kWARNING]    = 0.99 ; 
  lowValue[AliQAv1::kWARNING]   = 0.5 ; 
  lowValue[AliQAv1::kERROR]     = 0.0   ; 
  hiValue[AliQAv1::kERROR]      = 0.5 ; 
  lowValue[AliQAv1::kFATAL]     = -1.0   ; 
  hiValue[AliQAv1::kFATAL]      = 0.0 ; 
  SetHiLo(&hiValue[0], &lowValue[0]) ; 
}
