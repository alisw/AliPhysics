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
// $MpId: AliMpDetElement.cxx,v 1.4 2006/05/24 13:58:34 ivana Exp $
// Category: management
//
// Class AliMpDetElement
// --------------------
// The class defines the electronics properties of detection element
// Authors: Ivana Hrivnacova, IPN Orsay
//          Laurent Aphecetche, Christian Finck, SUBATECH Nantes

#include "AliMpDetElement.h"
#include "AliMpDEManager.h"

#include "AliLog.h"

#include <TObjString.h>
#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMpDetElement)
/// \endcond

const char  AliMpDetElement::fgkNameSeparator = '_'; 

//______________________________________________________________________________
AliMpDetElement::AliMpDetElement(Int_t id, const TString& name,
                    const TString& segType, AliMp::PlaneType planeType)
  : TObject(),
    fId(id),
    fDdlId(-1),
    fName(name),
    fSegType(segType),
    fPlaneType(planeType),
    fBusPatchIds(),
    fManuToSerialNbs(1700),
    fSerialNbToManus(1700)
{
/// Standard constructor

}

//______________________________________________________________________________
AliMpDetElement::AliMpDetElement(TRootIOCtor* /*ioCtor*/)
  : TObject(),
    fId(0),
    fDdlId(-1),
    fName(),
    fSegType(),
    fPlaneType(),
    fBusPatchIds(),
    fManuToSerialNbs(),
    fSerialNbToManus()
{
/// Root IO constructor
}

//______________________________________________________________________________
AliMpDetElement::~AliMpDetElement()
{
/// Destructor
}


//
// public methods
//

//______________________________________________________________________________
Bool_t AliMpDetElement::AddBusPatch(Int_t busPatchId)
{
/// Add bus patch Id if a bus patch with the same Id is not yet present;
/// return false if bus patch was not added

  if ( HasBusPatchId(busPatchId) ) {
    AliWarningStream() 
      << "Bus patch Id = " << busPatchId << " already present."
      << endl;
    return false;
  } 

  fBusPatchIds.Add(busPatchId); 
  return true;
}  
 
//______________________________________________________________________________
void AliMpDetElement::AddManuSerial(Int_t manuId, Int_t serialNb)
{
/// Map the serial manu number 
/// (Eventually add check if the given pair already present)

  fManuToSerialNbs.Add(Long_t(manuId), Long_t(serialNb)); 
  fSerialNbToManus.Add(Long_t(serialNb), Long_t(manuId));
}      

//______________________________________________________________________________
TString AliMpDetElement::GetSegName(AliMp::CathodType cathType) const
{
/// Return the segmentation name for the given catod type

  return fSegType + fgkNameSeparator + PlaneTypeName(GetPlaneType(cathType));
}	  

//______________________________________________________________________________
AliMp::PlaneType  AliMpDetElement::GetPlaneType(AliMp::CathodType cath) const 
{
/// Return plane type                                                      \n

  if ( cath == AliMp::kCath0 ) return fPlaneType;
  else                         return AliMp::OtherPlaneType(fPlaneType); 
}    

//______________________________________________________________________________
AliMp::CathodType AliMpDetElement::GetCathodType(AliMp::PlaneType planeType) const
{
/// Return cathod type for given planeType

  if ( fPlaneType == planeType ) return AliMp::kCath0;
  else                           return AliMp::kCath1;
}

//______________________________________________________________________________
AliMp::StationType AliMpDetElement::GetStationType() const
{
/// Return station type                                                      \n
/// Failure causes Fatal error - as AliMp::StationType has no possibility
/// to return undefined value

  Int_t chamberId = AliMpDEManager::GetChamberId(fId, false);
  if ( ! AliMpDEManager::IsValidChamberId(chamberId, true) ) {
    AliFatal("Cannot return AliMp::StationType value.");
    return AliMp::kStation1;
  }  
  
  if ( chamberId ==  0 || chamberId ==  1 )  return AliMp::kStation1;
  if ( chamberId ==  2 || chamberId ==  3 )  return AliMp::kStation2;
  if ( chamberId >=  4 && chamberId <=  9 )  return AliMp::kStation345;
  if ( chamberId >= 10 && chamberId <= 13 )  return AliMp::kStationTrigger;

  // Should never get to this line
  AliFatal("Cannot return AliMp::StationType value.");
  return AliMp::kStation1;
}

//______________________________________________________________________________
Int_t AliMpDetElement::GetNofBusPatches() const
{
/// Return the number of bus patches in this detection element

  return fBusPatchIds.GetSize();
}  

//______________________________________________________________________________
Int_t AliMpDetElement::GetBusPatchId(Int_t index) const
{
/// Return the index-th bus patch

  if ( index < 0 || index > GetNofBusPatches() ) {
    AliErrorStream()
      << "In DE = " << fId << ": Index " << index << " outside limits." << endl;
    return 0;
  }     

  return  fBusPatchIds.GetValue(index);
}   


//______________________________________________________________________________
Bool_t  AliMpDetElement::HasBusPatchId(Int_t busPatchId) const
{  
/// Return true if the bus patch Id is present

  return fBusPatchIds.HasValue(busPatchId);; 
}

//______________________________________________________________________________
Int_t  AliMpDetElement::GetNofManus() const
{
/// Return the number of manus in this detection element  

  return fManuToSerialNbs.GetSize();
}   

//______________________________________________________________________________
Int_t  AliMpDetElement::GetManuSerialFromId(Int_t manuId) const
{
/// Return manu serial number from manuId

  return (Int_t)fManuToSerialNbs.GetValue(Long_t(manuId));
}

//______________________________________________________________________________
Int_t  AliMpDetElement::GetManuIdFromSerial(Int_t serialNb) const
{
/// Return manuId from manu serial number
  
  return (Int_t)fSerialNbToManus.GetValue(Long_t(serialNb));
}

