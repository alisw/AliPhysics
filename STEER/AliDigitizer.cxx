/**************************************************************************
 * Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
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

//----------------------------------------------------------------------
//  Base Class for Detector specific Merging/Digitization   
//  Collaborates with AliRunDigitizer class                
//  Author: Jiri Chudoba (CERN)
//----------------------------------------------------------------------

// system includes
#include <Riostream.h>

// ROOT includes

// AliROOT includes
#include "AliLog.h"
#include "AliDigitizer.h"
#include "AliRunDigitizer.h"

ClassImp(AliDigitizer)

//_______________________________________________________________________
AliDigitizer::AliDigitizer(const Text_t* name, const Text_t* title):
  TTask(name,title),
  fManager(0),
  fRegionOfInterest(kTRUE)
{
  //
  // Default ctor with name and title
  //
}

//_______________________________________________________________________
AliDigitizer::AliDigitizer(const AliDigitizer &dig):
  TTask(dig.GetName(),dig.GetTitle()),
  fManager(0),
  fRegionOfInterest(kTRUE)
{
  //
  // Copy ctor with
  //
  dig.Copy(*this);
}

//_______________________________________________________________________
void AliDigitizer::Copy(TObject &) const
{
  AliFatal("Not yet implemented");
}

//_______________________________________________________________________
AliDigitizer::AliDigitizer(AliRunDigitizer *manager, 
                           const Text_t* name, const Text_t* title):
  TTask(name,title),
  fManager(manager),
  fRegionOfInterest(kFALSE)
{
  //
  // ctor with name and title
  //
  fManager->AddDigitizer(this);
}

//_______________________________________________________________________
AliDigitizer::~AliDigitizer() 
{
}

//_______________________________________________________________________
Int_t AliDigitizer::GetNInputStreams() const
{
  //
  // return number of input streams
  //
  Int_t nInputStreams = 0 ;
  if (fManager)
    nInputStreams = fManager->GetNinputs() ;
  return nInputStreams ; 
}
