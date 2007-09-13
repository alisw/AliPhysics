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

//...
//  Checks the quality assurance. 
//  By comparing with reference data
//  Skeleton for HMPID
//...

// --- ROOT system ---
#include <TClass.h>
#include <TH1F.h> 
#include <TH1I.h> 
#include <TIterator.h> 
#include <TKey.h> 
#include <TFile.h> 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQualAss.h"
#include "AliQualAssChecker.h"
#include "AliHMPIDQualAssChecker.h"

ClassImp(AliHMPIDQualAssChecker)

//__________________________________________________________________
AliHMPIDQualAssChecker& AliHMPIDQualAssChecker::operator = (const AliHMPIDQualAssChecker& qac )
{
  // Equal operator.
  this->~AliHMPIDQualAssChecker();
  new(this) AliHMPIDQualAssChecker(qac);
  return *this;
}

