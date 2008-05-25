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
 
// Base class for event pool.
// This class is needed by the AnalysisManager to steer a mixing analysis.
// Author Andreas Morsch
// andreas.morsch@cern.ch

#include "AliVEventPool.h"

ClassImp(AliVEventPool)


////////////////////////////////////////////////////////////////////////

AliVEventPool::AliVEventPool():
    TNamed(),
    fChain(0)
{
  // Default constructor
}

AliVEventPool::AliVEventPool(const char* name, const char* title):
    TNamed(name, title), 
    fChain()
{
  // Constructor
}


AliVEventPool::AliVEventPool(const AliVEventPool& obj):
    TNamed(obj),
    fChain(0)
{
    // Copy constructor
    fChain = obj.fChain;
}

AliVEventPool& AliVEventPool::operator=(const AliVEventPool& other)
{
// Assignment operator
    TNamed::operator=(other);
    fChain = other.fChain;
    return *this;
}
