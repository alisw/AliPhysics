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
#include "AliMpConstants.h"
#include "AliMpFiles.h"
#include "AliMpIntPair.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TSystem.h>
#include <TObjString.h>
#include <TMap.h>

/// \cond CLASSIMP
ClassImp(AliMpDEManager)
/// \endcond

const char  AliMpDEManager::fgkNameSeparator = '_'; 
const char  AliMpDEManager::fgkCommentPrefix = '#'; 
const Int_t AliMpDEManager::fgkCoefficient = 100;
AliMpExMap  AliMpDEManager::fgDENamesMap(true);
AliMpExMap  AliMpDEManager::fgDECathBNBMap(true);

//______________________________________________________________________________

AliMpDEManager::~AliMpDEManager()
{
/// Destructor
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
Int_t 
AliMpDEManager::GetCathod(Int_t detElemId, AliMpPlaneType planeType)
{
/// Return cathod number for given detElemId and planeType

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
  
  Int_t chamberId = GetChamberId(detElemId, false);
  if ( ! IsValidChamberId(chamberId, true) ) {
    AliFatalClass("Cannot return AliMpStationType value.");
    return kStation1;
  }  
  
  if ( chamberId ==  0 || chamberId ==  1 )  return kStation1;
  if ( chamberId ==  2 || chamberId ==  3 )  return kStation2;
  if ( chamberId >=  4 && chamberId <=  9 )  return kStation345;
  if ( chamberId >= 10 && chamberId <= 13 )  return kStationTrigger;

  // Should never get to this line
  AliFatalClass("Cannot return AliMpStationType value.");
  return kStation1;
}

