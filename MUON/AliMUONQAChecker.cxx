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

/// \cond CLASSIMP
ClassImp(AliMUONQAChecker)
/// \endcond

//__________________________________________________________________
AliMUONQAChecker::AliMUONQAChecker() : 
    AliQACheckerBase("MUON","MUON Quality Assurance Data Maker"),
fCheckers(new TObjArray)
{
	/// ctor
  fCheckers->SetOwner(kTRUE);
  fCheckers->Add(new AliMUONTrackerQAChecker());
  fCheckers->Add(new AliMUONTriggerQAChecker());
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
  
  TIter next(fCheckers);
  AliMUONVQAChecker* qac;
  const AliMUONRecoParam* muonRecoParam = static_cast<const AliMUONRecoParam*>(recoParam);
  AliMUONVQAChecker::ECheckCode* ecc(0x0);

  for ( Int_t i = 0; i < AliRecoParam::kNSpecies; ++i ) 
  {
    rv[i] = -1.0;
  }
  
  while ( ( qac = static_cast<AliMUONVQAChecker*>(next())) )
  {
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

//______________________________________________________________________________
void 
AliMUONQAChecker::SetQA(AliQAv1::ALITASK_t index, Double_t * value) const
{
  /// sets the QA according the return value of the Check

  AliQAv1 * qa = AliQAv1::Instance(index);
  
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    qa->UnSet(AliQAv1::kFATAL, specie);
    qa->UnSet(AliQAv1::kWARNING, specie);
    qa->UnSet(AliQAv1::kERROR, specie);
    qa->UnSet(AliQAv1::kINFO, specie);

    if ( ! value ) { // No checker is implemented, set all QA to Fatal
      qa->Set(AliQAv1::kFATAL, specie) ; 
    } else {
      if ( value[specie] >= fLowTestValue[AliQAv1::kFATAL] && value[specie] < fUpTestValue[AliQAv1::kFATAL] ) 
        qa->Set(AliQAv1::kFATAL, specie) ; 
      else if ( value[specie] > fLowTestValue[AliQAv1::kERROR] && value[specie] <= fUpTestValue[AliQAv1::kERROR]  )
        qa->Set(AliQAv1::kERROR, specie) ; 
      else if ( value[specie] > fLowTestValue[AliQAv1::kWARNING] && value[specie] <= fUpTestValue[AliQAv1::kWARNING]  )
        qa->Set(AliQAv1::kWARNING, specie) ;
      else if ( value[specie] > fLowTestValue[AliQAv1::kINFO] && value[specie] <= fUpTestValue[AliQAv1::kINFO] ) 
        qa->Set(AliQAv1::kINFO, specie) ; 	
    }
  }
}
