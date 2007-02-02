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
// $MpId: AliMpBusPatch.cxx,v 1.4 2006/05/24 13:58:34 ivana Exp $
//
// --------------------
// Class AliMpBusPatch
// --------------------
// The class defines the properties of BusPatch
// Author: Ivana Hrivnacova, IPN Orsay

#include "AliMpBusPatch.h"

#include "AliLog.h"

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMpBusPatch)
/// \endcond

const Int_t  AliMpBusPatch::fgkOffset = 100;
//
// static methods
//

//____________________________________________________________________
Int_t AliMpBusPatch::GetGlobalBusID(Int_t localID, Int_t ddlID)
{
  /// return global bus id from local bus and ddl id

  return ddlID*fgkOffset + localID;

}
//____________________________________________________________________
Int_t AliMpBusPatch::GetLocalBusID(Int_t globalID, Int_t ddlID)
{
  /// return local bus id from local bus id

  return globalID - ddlID*fgkOffset;

}

//______________________________________________________________________________
AliMpBusPatch::AliMpBusPatch(Int_t id, Int_t detElemId, Int_t ddlId)
  : TObject(),
    fId(id),
    fDEId(detElemId),
    fDdlId(ddlId),
    fManus(false)
{
/// Standard constructor
}

//______________________________________________________________________________
AliMpBusPatch::AliMpBusPatch(TRootIOCtor* /*ioCtor*/)
  : TObject(),
    fId(),
    fDEId(),
    fDdlId(),
    fManus()
{
/// Root IO constructor
}

//______________________________________________________________________________
AliMpBusPatch::~AliMpBusPatch()
{
/// Destructor
}

//
// public methods
//

//______________________________________________________________________________
Bool_t AliMpBusPatch::AddManu(Int_t manuId)
{
/// Add detection element with given detElemId.
/// Return true if the detection element was added

  if ( HasManu(manuId) ) {
    AliWarningStream() 
      << "Manu with manuId=" << manuId << " already present."
      << endl;
    return false;
  }    

  fManus.Add(manuId);
  return true;
}   

//______________________________________________________________________________
Int_t AliMpBusPatch::GetNofManus() const
{  
/// Return the number of detection elements connected to this DDL

  return fManus.GetSize(); 
}

//______________________________________________________________________________
Int_t  AliMpBusPatch::GetManuId(Int_t index) const
{  
/// Return the detection element by index (in loop)

  return fManus.GetValue(index); 
}

//______________________________________________________________________________
Bool_t  AliMpBusPatch::HasManu(Int_t manuId) const
{  
/// Return true if bus patch has manu with given manuId

  return fManus.HasValue(manuId); 
}


