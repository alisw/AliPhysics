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
// $MpId: AliMpDEStore.cxx,v 1.4 2006/05/24 13:58:34 ivana Exp $
// Category: management
//
// Class AliMpDEStore
// --------------------
// The container class for detection element objects
// Authors: Ivana Hrivnacova, IPN Orsay
//          Laurent Aphecetche, Christian Finck, SUBATECH Nantes

#include "AliMpDEStore.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "AliMpConstants.h"
#include "AliMpFiles.h"
#include "AliMpHelper.h"
#include "AliMpIntPair.h"
#include "AliMpConstants.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TClass.h>
#include <TSystem.h>
#include <TObjString.h>
#include <TMap.h>

/// \cond CLASSIMP
ClassImp(AliMpDEStore)
/// \endcond

AliMpDEStore* AliMpDEStore::fgInstance = 0;
const char    AliMpDEStore::fgkCommentPrefix = '#'; 

//
// static methods
//

//______________________________________________________________________________
AliMpDEStore* AliMpDEStore::Instance()
{
/// Create the DE store if it does not yet exist
/// and return its instance

  if ( ! fgInstance )
    fgInstance = new AliMpDEStore();
    
  return fgInstance;
}    

//
// ctors, dtor
//

//______________________________________________________________________________
AliMpDEStore::AliMpDEStore()
: TObject(),
  fDetElements(true)
{  
/// Standard constructor

  AliDebug(1,"");
  fDetElements.SetOwner(true);

  // Create all detection elements
  FillDEs();
}

//______________________________________________________________________________
AliMpDEStore::AliMpDEStore(TRootIOCtor* /*ioCtor*/)
: TObject(),
  fDetElements()
{  
/// Constructor for IO

  AliDebug(1,"");

  fgInstance = this;
}


//______________________________________________________________________________
AliMpDEStore::~AliMpDEStore()
{
/// Destructor

  AliDebug(1,"");

  // Segmentations are deleted with fMpSegmentations 
  // El cards arrays are deleted with fElCardsMap
  
  fgInstance = 0;
}

//
// private methods
//

//______________________________________________________________________________
Bool_t AliMpDEStore::IsPlaneType(const TString& planeTypeName)
{
/// Return true if the planeTypeName corresponds to a valid plane type

  if ( planeTypeName == PlaneTypeName(AliMp::kBendingPlane) ||
       planeTypeName == PlaneTypeName(AliMp::kNonBendingPlane) ) 
    return true;   

  return false;
}  

//______________________________________________________________________________
AliMp::PlaneType AliMpDEStore::PlaneType(const TString& planeTypeName)
{
/// Return plane type for the given planeTypeName                            \n
/// Fatal error if planeTypeName is wrong 

  if ( planeTypeName == PlaneTypeName(AliMp::kBendingPlane) ) 
    return AliMp::kBendingPlane;

  if ( planeTypeName == PlaneTypeName(AliMp::kNonBendingPlane) ) 
    return AliMp::kNonBendingPlane;

  // Should never reach this line
  AliFatalClass(Form("No plane type defined for %s", planeTypeName.Data()));
  return AliMp::kBendingPlane;
}       

//______________________________________________________________________________
AliMp::StationType AliMpDEStore::StationType(const TString& stationTypeName)
{
/// Return station type for the given stationTypeName                        \n
/// Fatal error if stationTypeName is wrong 

  if ( stationTypeName == StationTypeName(AliMp::kStation1) )
    return AliMp::kStation1;

  if ( stationTypeName == StationTypeName(AliMp::kStation2) )
    return AliMp::kStation2;

  if ( stationTypeName == StationTypeName(AliMp::kStation345) )
    return AliMp::kStation345;

  if ( stationTypeName == StationTypeName(AliMp::kStationTrigger) ) 
    return AliMp::kStationTrigger;

  // Should never reach this line
  AliFatalClass(Form("No station type defined for ", stationTypeName.Data()));
  return AliMp::kStation1;
}

//______________________________________________________________________________
Bool_t AliMpDEStore::ReadManuToSerialNbs(AliMpDetElement* detElement, 
                                       AliMp::StationType stationType)
{
/// Read manu serial numbers for the given detection element
  
  TString deName = detElement->GetDEName();

  TString infile = AliMpFiles::ManuToSerialPath(deName, stationType);
  ifstream in(infile, ios::in);
  
  // Change to Error when all files available
  //if ( !in.is_open() && stationType == AliMp::kStation345 ) {
  //   AliWarningStream() << "File " << infile << " not found." << endl;
  //  return false;
  //}   
       
  char line[80];

  while ( in.getline(line,80) ) {

    if ( line[0] == '#' ) continue;

    TString tmp(AliMpHelper::Normalize(line));

    Int_t blankPos  = tmp.First(' ');

    TString sManuId(tmp(0, blankPos));

    Int_t manuId = atoi(sManuId.Data());

    TString sManuSerial(tmp(blankPos + 1, tmp.Length()-blankPos));

    Int_t manuSerial = atoi(sManuSerial.Data());
      
    // filling manuId <> manuSerial
    detElement->AddManuSerial(manuId, manuSerial); 
  }
   
  in.close();
  return true;
}

//______________________________________________________________________________
Bool_t
AliMpDEStore::ReadDENames(AliMp::StationType station)
{ 
/// Read det element names for cath = 0 from the file specified by name
/// and fill the map 

  // Open file
  TString filePath = AliMpFiles::DENamesFilePath(station);
  std::ifstream in(filePath);
  if (!in.good()) {
    AliErrorClassStream() << "Cannot open file " << filePath << endl;;
    return false;
  }
  
  // Read plane types per cathods
  //
  char line[80];
  TString word;
  TString cathName1, cathName2;
  in >> word;
  while ( ! in.eof() && cathName1.Length() == 0 ) {
    if ( word[0] == '#' ) 
      in.getline(line, 80);
    else { 
      cathName1 = word;
      in >> cathName2;
    }
    in >> word;
  }
  
  Bool_t isCathNameDefined = false;
  if ( IsPlaneType(cathName1) &&  IsPlaneType(cathName2) )
    isCathNameDefined = true;
    
  // Read DE names
  //
  Int_t detElemId;
  TString name, name0, name1, name2;
  AliMp::PlaneType planeForCathode[2];
  
  while ( ! in.eof() ) 
  {
    if ( word[0] == '#' ) 
    {
      in.getline(line, 80);
    }
    else 
    {  
      detElemId = word.Atoi();
      in >> name;
      in >> name0;
      // warning : important to check non bending first (=nbp),
      // as bp is contained within nbp...
      if ( name0.Contains(PlaneTypeName(AliMp::kNonBendingPlane)) )
      {
        planeForCathode[0] = AliMp::kNonBendingPlane;
      }
      else
      {
        planeForCathode[0] = AliMp::kBendingPlane;
      }
 
      if ( !isCathNameDefined ) 
      { 
        in >> name2;
	name1 = name0; 
        Ssiz_t pos = name1.First(AliMpDetElement::GetNameSeparator());
        name0 = name1(0,pos);

        // Other cathode is other plane...
        planeForCathode[1] = OtherPlaneType(planeForCathode[0]);
      }
      else 
      {
        name1 = name0 + AliMpDetElement::GetNameSeparator() + cathName1;
        name2 = name0 + AliMpDetElement::GetNameSeparator() + cathName2;
        if ( name2.Contains(PlaneTypeName(AliMp::kNonBendingPlane)) )
        {
          planeForCathode[1] = AliMp::kNonBendingPlane;
        }
        else
        {
          planeForCathode[1] = AliMp::kBendingPlane;
        }        
      }   

      if ( planeForCathode[0]==planeForCathode[1] )
      {
        AliFatalClass(Form("Got the same cathode type for both planes"
                      " of DetElemId %d",detElemId));
      }
      
      AliMpDetElement* detElement 
        = new AliMpDetElement(detElemId, name, name0, planeForCathode[0]);
      if ( ! fDetElements.GetValue(detElemId) ) 
      {
        AliDebugClassStream(3)  
          << "Adding DE name "  << detElemId << "  " << name << endl;
        fDetElements.Add(detElemId, detElement); 
        
        // Read manu serial numbers for this det element
        ReadManuToSerialNbs(detElement, station);
      } 
      else 
      {
        AliWarningClassStream()
          << "Det element "  << detElemId << "  " << name << " already defined." << endl;
      }	
    } 
    in >> word;
  }

  // Close file
  in.close();
  
  return true;
}

//______________________________________________________________________________
void AliMpDEStore::FillDEs()
{
/// Fill DE names from files
  AliDebugClass(2,"");
  Bool_t result1 = ReadDENames(AliMp::kStation1);
  Bool_t result2 = ReadDENames(AliMp::kStation2);
  Bool_t result3 = ReadDENames(AliMp::kStation345);
  Bool_t result4 = ReadDENames(AliMp::kStationTrigger);
  
  Bool_t result = result1 && result2 && result3 && result4;
  if ( ! result ) {
    AliErrorClassStream() << "Error in reading DE names files" << endl;
  }  
}

//
// public methods
//


//______________________________________________________________________________
AliMpDetElement* AliMpDEStore::GetDetElement(Int_t detElemId, Bool_t warn) const
{
/// Return det element for given detElemId

  AliMpDetElement* detElement
    = (AliMpDetElement*)fDetElements.GetValue(detElemId);
    
  if ( ! detElement && warn ) {  
    AliErrorClassStream() 
        << "Detection element " << detElemId << " not defined." << endl;
  }	

  return detElement;
}    

//______________________________________________________________________________
AliMpIntPair  AliMpDEStore::GetDetElemIdManu(Int_t manuSerial) const
{
/// Return the detElemId and manuId for given serial manu number

  for ( Int_t i = 0; i < fDetElements.GetSize(); i++ ) {
    
    AliMpDetElement* detElement 
      = (AliMpDetElement*)fDetElements.GetObject(i);
      
    Int_t manuId = detElement->GetManuIdFromSerial(manuSerial);
    if ( manuId ) return AliMpIntPair(detElement->GetId(), manuId);
  }    

  // manu with this serial number does not exist
  return AliMpIntPair::Invalid();
}  
