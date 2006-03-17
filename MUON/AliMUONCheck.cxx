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

#include "AliMUONCheck.h"

#include "AliLoader.h"
#include "AliLog.h"
#include "AliRunLoader.h"
#include "AliMUONData.h"
#include "AliMUONDataIterator.h"
#include "AliMUONDigit.h"

/// A helper class to dump data from AliRoot-generated root files.
/// 
/// Only implemented for digits so far, it is meant as a replacement
/// of the MUONCheck.C macro, or to be used by this macro to simplify it.
///

ClassImp(AliMUONCheck)

//_____________________________________________________________________________
AliMUONCheck::AliMUONCheck(const char* galiceFile, 
                           Int_t firstEvent, Int_t lastEvent) 
: TObject(),
fFileName(galiceFile),
fFirstEvent(firstEvent),
fLastEvent(lastEvent),
fRunLoader(0x0),
fData(0x0)
{
  // ctor
  
  fRunLoader = AliRunLoader::Open(fFileName.Data(),"MUONFolder","READ");
  if (!fRunLoader) 
  {
    AliError(Form("Error opening %s file \n",fFileName.Data()));
  }  
  else
  {
    AliLoader* loader = fRunLoader->GetLoader("MUONLoader");
    if ( loader )
    {
      fData = new AliMUONData(loader,"MUON","MUON");
    }
    else
    {
      AliError(Form("Could get MUONLoader"));
    }
  }
}

//_____________________________________________________________________________
AliMUONCheck::AliMUONCheck(const AliMUONCheck& rhs) : TObject(rhs)
{
  // copy ctor
  AliFatal("Implement me if needed");
}

//_____________________________________________________________________________
AliMUONCheck& 
AliMUONCheck::operator=(const AliMUONCheck&)
{
  // assignement operator
  AliFatal("Implement me if needed")
  return *this;
}

//_____________________________________________________________________________
AliMUONCheck::~AliMUONCheck()
{
  // dtor
  delete fData;
}

//_____________________________________________________________________________
void
AliMUONCheck::DumpDigits(Option_t* opt) const
{
  // Dump the digits to screen
  if ( !IsValid() ) return;
  
  Int_t nevents = fRunLoader->GetNumberOfEvents();
  Int_t endOfLoop = fLastEvent+1;
  
  if ( fLastEvent == -1 ) endOfLoop = nevents;
  
  for ( Int_t ievent = fFirstEvent; ievent < endOfLoop; ++ievent ) 
  {
    fRunLoader->GetEvent(ievent);

    AliMUONDataIterator it(fData,"digit",AliMUONDataIterator::kTrackingChambers);
    AliMUONDigit* digit;
 
     while ( ( digit = (AliMUONDigit*)it.Next() ) )
     {
       digit->Print(opt);
     }
  } 
}
