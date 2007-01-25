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
// $MpId: AliMpDEManager.cxx,v 1.4 2006/05/24 13:58:34 ivana Exp $
// Category: management
//
// Class AliMpDEManager
// --------------------
// The manager class for definition of detection element types
// Authors: Ivana Hrivnacova, IPN Orsay
//          Laurent Aphecetche, SUBATECH Nantes

#include "AliMpDEManager.h"
#include "AliMpDEStore.h"
#include "AliMpDetElement.h"
#include "AliMpConstants.h"
#include "AliMpCathodType.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TClass.h>
//#include <TSystem.h>
//#include <TObjString.h>
//#include <TMap.h>

/// \cond CLASSIMP
ClassImp(AliMpDEManager)
/// \endcond

const Int_t AliMpDEManager::fgkCoefficient = 100;
TArrayI     AliMpDEManager::fgNofDEPerChamber;

//______________________________________________________________________________

AliMpDEManager::~AliMpDEManager()
{
/// Destructor
}

//
// static private methods
//

//______________________________________________________________________________
AliMpDetElement* AliMpDEManager::GetDetElement(Int_t detElemId, Bool_t warn)
{
/// Return det element for given detElemId

  return AliMpDEStore::Instance()->GetDetElement(detElemId, warn);
}    

//
// static public methods
//

//______________________________________________________________________________
Bool_t AliMpDEManager::IsValidDetElemId(Int_t detElemId, Bool_t warn)
{
/// Return true if detElemId is valid
/// (is present in the DE map)

  if ( GetDetElement(detElemId, warn) ) return true;

  return false;
}    

//______________________________________________________________________________
Bool_t AliMpDEManager::IsValidChamberId(Int_t chamberId, Bool_t warn)
{
/// Return true if chamberId is valid

  if ( chamberId >= 0 && chamberId < AliMpConstants::NofChambers() ) 
    return true;
 
  if (warn) 
    AliErrorClassStream() << "Wrong chamber Id " << chamberId << endl;
  
  return false;
}    

//______________________________________________________________________________
Bool_t AliMpDEManager::IsValidGeomModuleId(Int_t moduleId, Bool_t warn)
{
/// Return true if moduleId is valid

  if ( moduleId >= 0 && moduleId < AliMpConstants::NofGeomModules() ) 
    return true;
 
  if (warn) 
    AliErrorClassStream() << "Wrong module Id " << moduleId << endl;
  
  return false;
}    

//______________________________________________________________________________
Int_t  AliMpDEManager::GetChamberId(Int_t detElemId, Bool_t warn)
{
/// Return chamber Id for given detElemId

  if ( ! IsValidDetElemId(detElemId, warn) ) return -1;
  
  return detElemId/fgkCoefficient - 1;
}  

//______________________________________________________________________________
Int_t AliMpDEManager::GetGeomModuleId(Int_t detElemId, Bool_t warn)
{
/// <pre>
/// Get module Id from detection element Id                 
/// !!! moduleId != chamberId
/// Station 1:   Chamber:  1   Module:  0   Det elements:  100-103
///              Chamber:  2   Module:  1   Det elements:  200-203
/// Station 2:   Chamber:  3   Module:  2   Det elements:  300-303
///              Chamber:  4   Module:  3   Det elements:  400-403
/// Station 3:   Chamber:  5   Module:  4   Det elements:  500-504, 514-517
///                            Module:  5   Det elements:  505-513,
///              Chamber:  6   Module:  6   Det elements:  600-604, 614-617
///                            Module:  7   Det elements:  605-613
/// Station 4:   Chamber:  7   Module:  8   Det elements:  700-706, 720-725
///                            Module:  9   Det elements:  707-719
///              Chamber:  8   Module: 10   Det elements:  800-806, 820-825
///                            Module: 11   Det elements:  807-819
/// Station 5:   Chamber:  9   Module: 12   Det elements:  900-906, 920-925
///                            Module: 13   Det elements:  907-919        
///              Chamber: 10   Module: 14   Det elements: 1000-1006,1020-1025
///                            Module: 15   Det elements: 1007-1019
/// Station 6:   Chamber: 11   Module: 16   Det elements: 1100-1117
///              Chamber: 12   Module: 17   Det elements: 1200-1217
/// Station 7:   Chamber: 13   Module: 18   Det elements: 1300-1317
///              Chamber: 14   Module: 19   Det elements: 1400-1417
/// </pre>

  if ( ! IsValidDetElemId(detElemId, warn) ) return -1;
  
  return detElemId/fgkCoefficient 
           + (detElemId >=  505 && detElemId <=  513 || detElemId >= 600 )
	   + (detElemId >=  605 && detElemId <=  613 || detElemId >= 700 )
	   + (detElemId >=  707 && detElemId <=  719 || detElemId >= 800 )
	   + (detElemId >=  807 && detElemId <=  819 || detElemId >= 900 )
	   + (detElemId >=  907 && detElemId <=  919 || detElemId >= 1000 )
	   + (detElemId >= 1007 && detElemId <= 1019 || detElemId >= 1100 ) - 1;
}  

//______________________________________________________________________________
AliMp::PlaneType  AliMpDEManager::GetPlaneType(Int_t detElemId, AliMp::CathodType cath)
{
/// Return plane type                                                      \n
/// Failure causes Fatal error - as AliMp::PlaneType has no possibility
/// to return undefined value

  if ( ! IsValidDetElemId(detElemId, true) ) {
    AliFatalClass("Cannot return AliMp::PlaneType value.");
    return AliMp::kBendingPlane;
  }  

  return GetDetElement(detElemId)->GetPlaneType(cath);
}    

//______________________________________________________________________________
AliMp::StationType AliMpDEManager::GetStationType(Int_t detElemId)
{
/// Return station type                                                      \n
/// Failure causes Fatal error - as AliMp::StationType has no possibility
/// to return undefined value

  if ( ! IsValidDetElemId(detElemId, true) ) {
    AliFatalClass("Cannot return AliMp::StationType value.");
    return AliMp::kStation1;
  }  
  
  return GetDetElement(detElemId)->GetStationType();
}

//______________________________________________________________________________
AliMp::CathodType 
AliMpDEManager::GetCathod(Int_t detElemId, AliMp::PlaneType planeType)
{
/// Return cathod number for given detElemId and planeType

  if ( ! IsValidDetElemId(detElemId, true) ) {
    AliFatalClass("Cannot return AliMp::CathodType value.");
    return AliMp::kCath0;
  }  
  
  return GetDetElement(detElemId)->GetCathodType(planeType);
}

//______________________________________________________________________________
Int_t AliMpDEManager::GetNofDEInChamber(Int_t chamberId, Bool_t warn)
{
/// Return the number of detection elements in the chamber with the given 
/// chamberId

  if ( ! IsValidChamberId(chamberId,warn) ) return 0;

  // Fill array if it is empty
  if ( ! fgNofDEPerChamber.GetSize() ) {
    fgNofDEPerChamber.Set(AliMpConstants::NofChambers());
    AliMpDEIterator it;
    for ( Int_t i=0; i<AliMpConstants::NofChambers(); i++ ) {
      Int_t counter = 0;
      for ( it.First(i); ! it.IsDone(); it.Next() ) ++counter;
      fgNofDEPerChamber[i] = counter;
    }  
  }
  
  return fgNofDEPerChamber[chamberId];    

}

