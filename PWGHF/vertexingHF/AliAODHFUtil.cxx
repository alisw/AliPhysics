/**************************************************************************
 * Copyright(c) 1998-2011, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 * *************************************************************************/

/* $Id$ */

//***********************************************************
// Class AliAODHFUtil
// class for enabling access to data not available for the moment in AODs
// Author: Carlos Perez
//***********************************************************
#include "AliAODHFUtil.h"


ClassImp(AliAODHFUtil)

//------------------------------
AliAODHFUtil::AliAODHFUtil():
  TNamed()
{
 // default ctor
 for(int i=0; i!=64; ++i)
   fVZERO[i] = 0.0;
}

//------------------------------
AliAODHFUtil::AliAODHFUtil(const char* pName):
  TNamed(pName, "")
{
  // standard ctor
 for(int i=0; i!=64; ++i)
   fVZERO[i] = 0.0;
} 
//----------------------
AliAODHFUtil::~AliAODHFUtil()
{
 // default dtor
}
//------------------------
AliAODHFUtil::AliAODHFUtil(const AliAODHFUtil& pCopy) :
  TNamed(pCopy)
{
 // copy ctor
   for(int i=0; i!=64; ++i)
     fVZERO[i] = pCopy.fVZERO[i];
}
//----------------------
void AliAODHFUtil::SetVZERO(Float_t *pVzero) {
  // set VZERO channel
   for(int i=0; i!=64; ++i)
     fVZERO[i] = pVzero[i];
}
//---------------------------
Float_t AliAODHFUtil::GetVZEROChannel(Int_t pCh) const {
  // get VZERO channel
  return fVZERO[pCh];
}

