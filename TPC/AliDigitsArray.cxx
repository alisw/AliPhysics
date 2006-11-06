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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  AliDigitsArray  object                                //
//
//  Origin: Marian Ivanov , GSI Darmstadt
//                                                                           //
//                                                                          //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "AliSegmentID.h"
#include "TObjArray.h"

#include "TError.h"
#include "AliDigits.h"
#include "AliDetectorParam.h"
#include "AliDigitsArray.h"



ClassImp(AliDigitsArray)
//

AliDigitsArray::AliDigitsArray()
               :AliSegmentArray(),
		fParam(0)
{
  //
  // default constructor
  //
}
AliDigitsArray::AliDigitsArray(const AliDigitsArray &param)
             :AliSegmentArray(),
	      fParam(0)
{
  //
  // dummy
  //
  fParam = param.fParam;
}
//
AliDigitsArray & AliDigitsArray::operator =(const AliDigitsArray & param)
{
  //
  // dummy
  //
 fParam = param.fParam;
 return (*this);
}
//
AliDigitsArray::~AliDigitsArray()
{
  // if (fParam != 0) delete fParam;
}  

Bool_t AliDigitsArray::Setup(AliDetectorParam *param)
{
  //
  //setup array according parameters
  SetParam(param);
  return kTRUE;
}


Bool_t AliDigitsArray::SetParam(AliDetectorParam * param)
{
  fParam = param;
  return kTRUE;
}
