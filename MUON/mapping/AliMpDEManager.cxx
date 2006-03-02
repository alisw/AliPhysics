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
// $MpId: AliMpDEManager.cxx,v 1.2 2006/03/02 16:30:09 ivana Exp $
// Category: management
//
// Class AliMpDEManager
// --------------------
// The manager class for definition of detection element types
// Authors: Ivana Hrivnacova, IPN Orsay
//          Laurent Aphecetche, SUBATECH Nantes

#include "AliMpDEManager.h"
#include "AliMpConstants.h"
#include "AliMpFiles.h"
#include "AliMpIntPair.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TSystem.h>
#include <TObjString.h>
#include <TMap.h>

const char  AliMpDEManager::fgkNameSeparator = '_'; 
const char  AliMpDEManager::fgkCommentPrefix = '#'; 
const Int_t AliMpDEManager::fgkCoefficient = 100;
AliMpExMap  AliMpDEManager::fgDENamesMap(true);
AliMpExMap  AliMpDEManager::fgDECathBNBMap(true);

ClassImp(AliMpDEManager)

//______________________________________________________________________________
AliMpDEManager::AliMpDEManager()
    : TObject()
{  
/// Protected default/standard constructor
}

//______________________________________________________________________________
AliMpDEManager::AliMpDEManager(const AliMpDEManager& rhs)
 : TObject(rhs)
{
/// Protected copy constructor

  AliFatal("Not implemented.");
}

//______________________________________________________________________________

AliMpDEManager::~AliMpDEManager()
{
/// Destructor
}

//______________________________________________________________________________
AliMpDEManager&  AliMpDEManager::operator=(const AliMpDEManager& rhs)
{
/// Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal("Not implemented.");
    
  return *this;  
}    
          
//
// static private methods
//

//______________________________________________________________________________
Bool_t AliMpDEManager::IsPlaneType(const TString& planeTypeName)
{
/// Return true if the planeTypeName corresponds to a valid plane type

  if ( planeTypeName == PlaneTypeName(kBendingPlane) ||
       planeTypeName == PlaneTypeName(kNonBendingPlane) ) 
    return true;   

  return false;
}  

//______________________________________________________________________________
AliMpPlaneType AliMpDEManager::PlaneType(const TString& planeTypeName)
{
/// Return plane type for the given planeTypeName                            \n
/// Fatal error if planeTypeName is wrong 

  if ( planeTypeName == PlaneTypeName(kBendingPlane) ) 
    return kBendingPlane;

  if ( planeTypeName == PlaneTypeName(kNonBendingPlane) ) 
    return kNonBendingPlane;

  // Should never reach this line
  AliFatalClass(Form("No plane type defined for %s", planeTypeName.Data()));
  return kBendingPlane;
}       

//______________________________________________________________________________
AliMpStationType AliMpDEManager::StationType(const TString& stationTypeName)
{
/// Return station type for the given stationTypeName                        \n
/// Fatal error if stationTypeName is wrong 

  if ( stationTypeName == StationTypeName(kStation1) )
    return kStation1;

  if ( stationTypeName == StationTypeName(kStation2) )
    return kStation2;

  if ( stationTypeName == StationTypeName(kStation345) )
    return kStation345;

  if ( stationTypeName == StationTypeName(kStationTrigger) ) 
    return kStationTrigger;

  // Should never reach this line
  AliFatalClass(Form("No station type defined for ", stationTypeName.Data()));
  return kStation1;
}

//______________________________________________________________________________
Bool_t
AliMpDEManager::ReadDENames(AliMpStationType station)
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
  TString name1, name2;
  AliMpPlaneType planeForCathode[2];
  
  while ( ! in.eof() ) 
  {
    if ( word[0] == '#' ) 
    {
      in.getline(line, 80);
    }
    else 
    {  
      detElemId = word.Atoi();
      in >> name1;
      // warning : important to check non bending first (=nbp),
      // as bp is contained within nbp...
      if ( name1.Contains(PlaneTypeName(kNonBendingPlane)) )
      {
        planeForCathode[0] = kNonBendingPlane;
      }
      else
      {
        planeForCathode[0] = kBendingPlane;
      }
      if ( !isCathNameDefined ) 
      {       
        in >> name2;
        // Other cathode is other plane...
        if ( planeForCathode[0] == kBendingPlane ) 
        {
          planeForCathode[1]=kNonBendingPlane;
        }
        else
        {
          planeForCathode[1]=kBendingPlane;
        }
      }
      else 
      {
        name1 += fgkNameSeparator;
        name2 = name1;
        name1 += cathName1;
        name2 += cathName2;
        if ( name2.Contains(PlaneTypeName(kNonBendingPlane)) )
        {
          planeForCathode[1] = kNonBendingPlane;
        }
        else
        {
          planeForCathode[1] = kBendingPlane;
        }        
      }   

      if ( planeForCathode[0]==planeForCathode[1] )
      {
        AliFatalClass(Form("Got the same cathode type for both planes"
                      " of DetElemId %d",detElemId));
      }
      
      if ( ! fgDENamesMap.GetValue(detElemId) ) 
      {
        AliDebugClassStream(1)  
        << "Adding  "  << detElemId << "  " << name1 << "  " << name2 << endl;
        fgDENamesMap.Add(detElemId, 
                         new TPair(new TObjString(name1), new TObjString(name2)));
        fgDECathBNBMap.Add(detElemId,
                           new AliMpIntPair(planeForCathode[0],planeForCathode[1]));
      } 
    } 
    in >> word;
  }

  // Close file
  in.close();
  
  return true;
}

//______________________________________________________________________________
void AliMpDEManager::FillDENames()
{
/// Fill DE names from files

  Bool_t result1 = ReadDENames(kStation1);
  Bool_t result2 = ReadDENames(kStation2);
  Bool_t result3 = ReadDENames(kStation345);
  Bool_t result4 = ReadDENames(kStationTrigger);
  
  Bool_t result = result1 && result2 && result3 && result4;
  if ( ! result ) {
    AliErrorClassStream() << "Error in reading DE names files" << endl;
  }  
}

//
// static public methods
//

//______________________________________________________________________________
Bool_t AliMpDEManager::IsValidDetElemId(Int_t detElemId, Bool_t warn)
{
/// Return true if detElemId is valid
/// (is present in the DE names files)

  if ( fgDENamesMap.GetSize() == 0 ) FillDENames();

  if ( fgDENamesMap.GetValue(detElemId) ) return true;

  if (warn) {
    AliErrorClassStream() 
        << "Detection element " << detElemId << " not defined." << endl;
  }	
  return false;
}    

//______________________________________________________________________________
Bool_t AliMpDEManager::IsValidCathod(Int_t cath, Bool_t warn)
{
/// Return true if cath is 0 or 1 
/// (Better solution would be to use systematically enum)

  if (cath == 0 || cath == 1 ) return true;  

  if (warn)
    AliErrorClassStream() << "Wrong cathod number " << cath << endl;
     
  return false;
}    


//______________________________________________________________________________
Bool_t AliMpDEManager::IsValid(Int_t detElemId, Int_t cath, Bool_t warn)
{
/// Return true if both detElemId and cathod number are valid

  return ( IsValidDetElemId(detElemId, warn) && IsValidCathod(cath, warn) );
}    

//______________________________________________________________________________
Bool_t AliMpDEManager::IsValidModuleId(Int_t moduleId, Bool_t warn)
{
/// Return true if moduleId is valid

  if ( moduleId >= 0 && moduleId < AliMpConstants::NCh() ) 
    return true;
 
  if (warn) 
    AliErrorClassStream() << "Wrong module Id " << moduleId << endl;
  
  return false;
}    

//______________________________________________________________________________
Int_t 
AliMpDEManager::GetCathod(Int_t detElemId, AliMpPlaneType planeType)
{
  if ( !IsValidDetElemId(detElemId) ) return -1;
  AliMpIntPair* pair = 
    static_cast<AliMpIntPair*>(fgDECathBNBMap.GetValue(detElemId));
  if ( planeType == pair->GetFirst() )
  {
    return 0;
  }
  return 1;
}

//______________________________________________________________________________
TString AliMpDEManager::GetDEName(Int_t detElemId, Int_t cath, Bool_t warn)
{
/// Return det element type name

  if ( ! IsValid(detElemId, cath, warn) ) return "undefined";

  TPair* namePair = (TPair*)fgDENamesMap.GetValue(detElemId);

  if (cath == 0) return ((TObjString*)namePair->Key())->GetString();
  if (cath == 1) return ((TObjString*)namePair->Value())->GetString();
  
  return "undefined";
}    

//______________________________________________________________________________
TString AliMpDEManager::GetDETypeName(Int_t detElemId, Int_t cath, Bool_t warn) 
{
/// Return det element type name

  TString fullName = GetDEName(detElemId, cath, warn);  
  
  // cut plane type extension
  Ssiz_t pos = fullName.First(fgkNameSeparator);
  return fullName(0,pos);
}    

//______________________________________________________________________________
Int_t  AliMpDEManager::GetModuleId(Int_t detElemId, Bool_t warn)
{
/// Return module Id for given detElemId

  if ( ! IsValidDetElemId(detElemId, warn) ) return -1;
  
  return detElemId/fgkCoefficient - 1;
}  

//______________________________________________________________________________
AliMpPlaneType  AliMpDEManager::GetPlaneType(Int_t detElemId, Int_t cath)
{
/// Return plane type                                                      \n
/// Failure causes Fatal error - as AliMpPlaneType has no possibility
/// to return undefined value

  if ( ! IsValid(detElemId, cath, true) ) {
    AliFatalClass("Cannot return AliMpPlaneType value.");
    return kBendingPlane;
  }  

  TPair* namePair = (TPair*)fgDENamesMap.GetValue(detElemId);

  TString fullName;  
  if (cath == 0) fullName = ((TObjString*)namePair->Key())->GetString();
  if (cath == 1) fullName = ((TObjString*)namePair->Value())->GetString();
  
  // Get plane type name
  Ssiz_t pos = fullName.First(fgkNameSeparator);
  TString planeTypeName = fullName(pos+1,fullName.Length()-pos);
  
  return PlaneType(planeTypeName);
}    

//______________________________________________________________________________
AliMpStationType AliMpDEManager::GetStationType(Int_t detElemId)
{
/// Return station type                                                      \n
/// Failure causes Fatal error - as AliMpStationType has no possibility
/// to return undefined value

  if ( ! IsValidDetElemId(detElemId, true) ) {
    AliFatalClass("Cannot return AliMpStationType value.");
    return kStation1;
  }  
  
  Int_t moduleId = GetModuleId(detElemId, false);
  if ( ! IsValidModuleId(moduleId, true) ) {
    AliFatalClass("Cannot return AliMpStationType value.");
    return kStation1;
  }  
  
  if ( moduleId ==  0 || moduleId ==  1 )  return kStation1;
  if ( moduleId ==  2 || moduleId ==  3 )  return kStation2;
  if ( moduleId >=  4 && moduleId <=  9 )  return kStation345;
  if ( moduleId >= 10 && moduleId <= 13 )  return kStationTrigger;

  // Should never get to this line
  AliFatalClass("Cannot return AliMpStationType value.");
  return kStation1;
}

