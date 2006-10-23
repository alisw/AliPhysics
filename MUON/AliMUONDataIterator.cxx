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

#include "AliMUONDataIterator.h"

#include "AliMUONConstants.h"
#include "AliMUONData.h"
#include "AliMUONDataDigitIterator.h"
#include "TString.h"

/// \class AliMUONDataIterator
/// A wrapper to various iterators used to loop over
/// objects handled by AliMUONData, like sdigits, digits, rawclusters, 
/// and so on.
/// Currently only implemented for digits, as a proof-of-principle.
///
/// \author Laurent Aphecetche

/// \cond CLASSIMP
ClassImp(AliMUONDataIterator)
/// \endcond

namespace
{
  void GetChamberNumbers(AliMUONDataIterator::EIterationStyle type, 
                         Int_t& firstChamber, Int_t& lastChamber)
{
    switch ( type )
    {
      case AliMUONDataIterator::kAllChambers:
        firstChamber=0;
        lastChamber=AliMUONConstants::NCh()-1;
        break;
      case AliMUONDataIterator::kTrackingChambers:
        firstChamber=0;
        lastChamber=AliMUONConstants::NCh()-1;
        break;
      case AliMUONDataIterator::kTriggerChambers:
        firstChamber=AliMUONConstants::NTrackingCh();
        lastChamber=AliMUONConstants::NCh()-1;
        break;
      default:
        firstChamber=lastChamber=-1;
        break;
    }
}
}

//_____________________________________________________________________________
AliMUONDataIterator::AliMUONDataIterator() 
: 
TObject(),
fIterator(0x0)
{
}

//_____________________________________________________________________________
AliMUONDataIterator::AliMUONDataIterator(AliMUONData* data,
                                         const char* onWhatToIterate,
                                         EIterationStyle howToIterate) 
: 
TObject(),
fIterator(0x0)
{
  TString opt(onWhatToIterate);
  opt.ToLower();
  if ( opt.Contains("digit") || opt.Contains("d") )
  {
      Int_t firstChamber;
      Int_t lastChamber;
      GetChamberNumbers(howToIterate,firstChamber,lastChamber); 
      if ( firstChamber >= 0 && lastChamber >= 0 )
      {
        data->GetLoader()->LoadDigits("READ");
        data->SetTreeAddress("D,GLT");
        fIterator = new AliMUONDataDigitIterator(data,firstChamber,lastChamber);
      }
  }
}

//_____________________________________________________________________________
AliMUONDataIterator::~AliMUONDataIterator()
{
  delete fIterator;
}

//_____________________________________________________________________________
TObject* 
AliMUONDataIterator::Next() 
{ 
  if (fIterator) return fIterator->Next(); 
  return 0x0;
}

//_____________________________________________________________________________
Bool_t 
AliMUONDataIterator::Remove() 
{ 
  if (fIterator) return fIterator->Remove(); 
  return kFALSE;
}

//_____________________________________________________________________________
void 
AliMUONDataIterator::Reset() 
{ 
  if (fIterator) fIterator->Reset(); 
}
