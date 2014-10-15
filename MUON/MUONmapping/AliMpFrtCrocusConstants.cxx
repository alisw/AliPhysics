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
// $MpId: AliMpFrtCrocusConstants.cxx,v 1.4 2006/05/24 13:58:34 ivana Exp $

//-----------------------------------------------------------------------------
// Class AliMpFrtCrocusConstants
// --------------------
// The class defines the constants for FRT Crocus
// Author: Ch. Finck, Subatech Nantes
//-----------------------------------------------------------------------------

#include "AliMpFrtCrocusConstants.h"

/// \cond CLASSIMP
ClassImp(AliMpFrtCrocusConstants)
/// \endcond

const Int_t   AliMpFrtCrocusConstants::fgkLinkPorts[10] = {0, 1, 2, 3, 5, 0, 1, 2, 3, 5};
const Int_t   AliMpFrtCrocusConstants::fgkOffset  = 5;
const Int_t   AliMpFrtCrocusConstants::fgkNofDsps = 2;
const Int_t   AliMpFrtCrocusConstants::fgkNofBusPatches  = 10;
const UInt_t  AliMpFrtCrocusConstants::fgkBaseAddress   = 0x00040000;
const UInt_t  AliMpFrtCrocusConstants::fgkAddressOffset = 0x00010000;

//____________________________________________________________________
Int_t AliMpFrtCrocusConstants::GetGlobalFrtID(Int_t localID, Int_t ddlID)
{
  /// return global bus id from local frt and ddl id

  return ddlID*fgkOffset + localID;

}
//____________________________________________________________________
Int_t AliMpFrtCrocusConstants::GetLocalFrtID(Int_t globalID, Int_t ddlID)
{
  /// return local bus id from local frt id

  return globalID - ddlID*fgkOffset;

}

//______________________________________________________________________________
AliMpFrtCrocusConstants::AliMpFrtCrocusConstants()
  : TObject()
{
/// Standard constructor
}


//______________________________________________________________________________
AliMpFrtCrocusConstants::~AliMpFrtCrocusConstants()
{
/// Destructor
}

//
// public methods
//

//______________________________________________________________________________
UInt_t  AliMpFrtCrocusConstants::GetTopAddress(Int_t id) 
{
/// return WME top address

  Int_t localFrtId = id % fgkOffset;
  
  return  fgkBaseAddress + 2*localFrtId * fgkAddressOffset;

}

//______________________________________________________________________________
Int_t  AliMpFrtCrocusConstants::GetIdFromTopAddress(UInt_t add) 
{
/// return id from WME top address
  
  return  (add - fgkBaseAddress)/(2*fgkAddressOffset);

}


//______________________________________________________________________________
UInt_t  AliMpFrtCrocusConstants::GetBotAddress(Int_t id)
{
/// return WME bottom address

  Int_t localFrtId = id % fgkOffset; 
  
  return  fgkBaseAddress + (2*localFrtId+1) * fgkAddressOffset;
 
}

//______________________________________________________________________________
Int_t  AliMpFrtCrocusConstants::GetIdFromBotAddress(UInt_t add)
{
/// return id from WME bottom address
  
  return  (add - fgkBaseAddress - fgkAddressOffset)/(2*fgkAddressOffset);
 
}


//______________________________________________________________________________
MpPair_t  AliMpFrtCrocusConstants::GetLinkPortId(Int_t index) 
{  
/// Return the linkPort/dspId by index

  if ( index >= fgkNofBusPatches ) return -1;
  
  Int_t dspId;
  if (index < fgkOffset)
    dspId = 0;
  else
    dspId = 1;  
  
  return AliMp::Pair(dspId, fgkLinkPorts[index]);
 
}

//______________________________________________________________________________
Int_t AliMpFrtCrocusConstants::GetNofDsps()
{  
/// Return the number of DSPs connected to this FRT

  return fgkNofDsps; 
}

//______________________________________________________________________________
Int_t AliMpFrtCrocusConstants::GetNofBusPatches()
{  
/// Return the number of BusPatches connected to this FRT

  return fgkNofBusPatches; 
}


