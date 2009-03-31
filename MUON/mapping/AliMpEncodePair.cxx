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

#include "AliMpEncodePair.h"

#include <Riostream.h>

//_______________________________________________________________________
MpPair_t AliMp::Pair(Int_t first, Int_t second)
{
/// See also AliMp::PairFirst(), AliMp::PairSecond()
/// \author L. Aphecetche, SUBATECH

  return (( first << 16 ) | second);
}

//_______________________________________________________________________
Int_t AliMp::PairFirst(MpPair_t pair )
{
/// See also AliMp::Pair(), AliMp::PairSecond()
/// \author L. Aphecetche, SUBATECH

  return ( pair & 0xFFFF0000 ) >> 16;
}

//_______________________________________________________________________
Int_t AliMp::PairSecond(MpPair_t pair)
{
/// See also AliMp::Pair(), AliMp::PairFirst()
/// \author L. Aphecetche, SUBATECH

  return pair & 0xFFFF;
}

//_______________________________________________________________________
ostream& AliMp::PairPut(ostream& stream, MpPair_t pair)
{
  if ( pair >= 0 ) {
    stream << '(' << AliMp::PairFirst(pair) 
           << ',' << AliMp::PairSecond(pair) << ')';
    return stream;
  }  
  else { 
    stream << "AliMpIntPair::Invalid";
    return stream;
  }
}    
