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

////////////////////////////////////////////////////////////////////////
//
//  Pure Virtual Base Class for Detector specific Merging/Digitization   
//                  
//  Author: Jiri Chudoba (CERN)
//
////////////////////////////////////////////////////////////////////////

/*
$Log$
*/

// system includes
#include <iostream.h>

// ROOT includes

// AliROOT includes
#include "AliDigitizer.h"
#include "AliRunDigitizer.h"

ClassImp(AliDigitizer)

AliDigitizer::AliDigitizer() :TNamed("AliDigitizer","")
{
// dummy default ctor
  ;
}

AliDigitizer::AliDigitizer(AliRunDigitizer* manager) 
  :TNamed("AliDigitizer","")
{
  manager->AddDigitizer(this);
  fDebug=1;
}

AliDigitizer::~AliDigitizer() {;}
