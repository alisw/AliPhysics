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
// $MpId: AliMpVRowSegment.cxx,v 1.6 2006/05/24 13:58:46 ivana Exp $
// Category: sector

//-----------------------------------------------------------------------------
// Class AliMpVRowSegment
// ----------------------
// Class describing an interface for a row segment.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMpVRowSegment.h"

#include <iostream>
#include "TString.h"

/// \cond CLASSIMP
ClassImp(AliMpVRowSegment)
/// \endcond

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

//_____________________________________________________________________________
void AliMpVRowSegment::Print(Option_t* opt) const
{
    std::cout << opt << Form("%30s has %2d motifs",ClassName(),
            GetNofMotifs());
    std::cout << Form(" left,right=%7.2f,%7.2f bottom,top=%7.2f,%7.2f",
        GetPositionX()-GetDimensionX(),
        GetPositionX()+GetDimensionX(),
        GetPositionY()-GetDimensionY(),
        GetPositionY()+GetDimensionY());
    std::cout << Form(" ix=%3d,%3d iy=%3d,%3d",
            GetLowLimitIx(),
            GetHighLimitIx(),
            GetLowLimitIy(),
            GetHighLimitIy())
    << std::endl;
}

