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
// Base class for selection classes for the correction framework. 
// Inherits from AliAnalysisCuts. It includes additional methods 
// to export QA histograms & study the cut statistics & correlations 
// through the bitmap of the cuts embedded in each class, if needed
// Author S.Arcelli
// silvia.Arcelli@cern.ch

#include "AliCFCutBase.h"


ClassImp(AliCFCutBase)


//___________________________________________________________________________
AliCFCutBase::AliCFCutBase():
  AliAnalysisCuts(),
  fIsQAOn(kFALSE)
{
  //
  // Default constructor
  //
}

//___________________________________________________________________________
AliCFCutBase::AliCFCutBase(const char* name, const char* title):
  AliAnalysisCuts(name, title),
  fIsQAOn(kFALSE)
{
  //
  // Constructor
  //
}

//___________________________________________________________________________
AliCFCutBase::AliCFCutBase(const AliCFCutBase& obj):
  AliAnalysisCuts(obj),
  fIsQAOn(obj.fIsQAOn)
{
  //
  // Copy Constructor
  //
}
