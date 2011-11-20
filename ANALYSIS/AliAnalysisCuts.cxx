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

/* $Id$ */
 
// Base class for analysis cuts
// Author Andreas Morsch
// andreas.morsch@cern.ch

#include <TObject.h>
#include "AliAnalysisCuts.h"


ClassImp(AliAnalysisCuts)


////////////////////////////////////////////////////////////////////////

AliAnalysisCuts::AliAnalysisCuts():
    AliVCuts(), fFilterMask(0), fSelected(kFALSE)
{
  // Default constructor
}

AliAnalysisCuts::AliAnalysisCuts(const char* name, const char* title):
    AliVCuts(name, title), fFilterMask(0), fSelected(kFALSE)
{
  // Constructor
}

AliAnalysisCuts::AliAnalysisCuts(const AliAnalysisCuts& obj):
    AliVCuts(obj), fFilterMask(obj.fFilterMask), fSelected(obj.fSelected)
{
}
