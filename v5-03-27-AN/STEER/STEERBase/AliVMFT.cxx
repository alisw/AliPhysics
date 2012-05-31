/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

//-------------------------------------------------------------------------
//     Base class for ESD and AOD MFT data
//     Author: Cvetan Cheshkov
//     cvetan.cheshkov@cern.ch 2/02/2011
//-------------------------------------------------------------------------

#include "AliVMFT.h"
#include "AliLog.h"

ClassImp(AliVMFT)

//__________________________________________________________________________
AliVMFT::AliVMFT(const AliVMFT& source) :
  TObject(source) { } // Copy constructor

//__________________________________________________________________________
AliVMFT& AliVMFT::operator=(const AliVMFT& source)
{
  // Assignment operator
  //
  if (this!=&source) { 
    TObject::operator=(source); 
  }
  
  return *this; 
}

