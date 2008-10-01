/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Author: Per Thomas Hille  <perthi@fys.uio.no>                  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliAltroBunch.h"

ClassImp(AliAltroBunch)

AliAltroBunch::AliAltroBunch() : 
  fData(NULL),
  fBunchSize(-1),
  fEndTimeBin(0),
  fStartTimeBin(999),
  fPrewBunchSize(0),
  fPrevEndTimeBin(0)
//  fIsFirstBunch(true)
{
  // Default constructor
}


AliAltroBunch::~AliAltroBunch()
{
  // Default destructor
}

