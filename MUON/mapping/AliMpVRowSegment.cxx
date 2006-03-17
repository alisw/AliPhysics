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

// $Id$
// $MpId: AliMpVRowSegment.cxx,v 1.5 2006/03/17 11:38:43 ivana Exp $
// Category: sector
//
// Class AliMpVRowSegment
// ----------------------
// Class describing an interface for a row segment.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpVRowSegment.h"

ClassImp(AliMpVRowSegment)

//_____________________________________________________________________________
AliMpVRowSegment::AliMpVRowSegment()
  : AliMpVIndexed()
{
/// Default constructor
}

//_____________________________________________________________________________
AliMpVRowSegment::~AliMpVRowSegment() 
{
/// Destructor  
}

//_____________________________________________________________________________
AliMpVPadIterator* AliMpVRowSegment::CreateIterator() const
{
/// Give Fatal if iterator is not implemented in the derived class 

  Fatal("CreateIterator", "Iterator is not yet implemented.");
  
  return 0;
}  



