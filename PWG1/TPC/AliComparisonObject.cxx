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

//------------------------------------------------------------------------------
// Implementation of abstract AliComparisonObject class. It keeps information from 
// comparison of reconstructed and MC particle tracks. 
//
// Author: J.Otwinowski 14/04/2008 
//------------------------------------------------------------------------------

#include <iostream>

#include "AliMCInfo.h" 
#include "AliESDRecInfo.h" 
#include "AliComparisonObject.h" 

using namespace std;

ClassImp(AliComparisonObject)

//_____________________________________________________________________________
AliComparisonObject::AliComparisonObject():
  TNamed("AliComparisonObject","AliComparisonObject"),
  fAnalysisMode(-1),
  fHptGenerator(kFALSE)
{
  // constructor
}

//_____________________________________________________________________________
AliComparisonObject::AliComparisonObject(const char* name, const char* title):
  TNamed(name,title),
  fAnalysisMode(-1),
  fHptGenerator(kFALSE)
{
  // constructor
}

//_____________________________________________________________________________
AliComparisonObject::~AliComparisonObject(){
  // destructor 
}

