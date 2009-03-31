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

//-----------------------------------------------------------------------------
// Class AliMpDEStore
// --------------------
// The container class for detection element objects
// Authors: Ivana Hrivnacova, IPN Orsay
//          Laurent Aphecetche, Christian Finck, SUBATECH Nantes
//-----------------------------------------------------------------------------

#include <cstdlib>
#include "AliMpDEStore.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "AliMpConstants.h"
#include "AliMpFiles.h"
#include "AliMpDataStreams.h"
#include "AliMpHelper.h"
#include "AliMpConstants.h"
#include "AliMpExMapIterator.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TClass.h>
#include <TSystem.h>
#include <TObjString.h>
#include <TObjArray.h>
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
AliMpDEStore* AliMpDEStore::Instance(Bool_t warn)
{
/// Create the DE store if it does not yet exist
/// and return its instance

  if ( ! fgInstance && warn  ) {
    AliWarningClass("DE Store has not been loaded");
  }  
     
  return fgInstance;
}    

//______________________________________________________________________________
AliMpDEStore* AliMpDEStore::ReadData(const AliMpDataStreams& dataStreams, 
                                     Bool_t warn)
{
/// Load the DE store data from ASCII data files
/// and return its instance

  if ( fgInstance ) {
    if ( warn )
      AliWarningClass("DE Store has been already loaded");
    return fgInstance;
  }  
  
  if ( dataStreams.GetReadFromFiles() )
    AliInfoClass("Reading DE Store from ASCII files.");

  fgInstance = new AliMpDEStore(dataStreams);
  return fgInstance;
}    

//
// ctors, dtor
//

//______________________________________________________________________________
AliMpDEStore::AliMpDEStore(const AliMpDataStreams& dataStreams)
: TObject(),
  fkDataStreams(dataStreams),
  fDetElements()
{  
/// Standard constructor

  AliDebug(1,"");
  fDetElements.SetOwner(true);

  // Create all detection elements
  FillDEs();
}

//______________________________________________________________________________
AliMpDEStore::AliMpDEStore(TRootIOCtor* ioCtor)
: TObject(),
  fkDataStreams(ioCtor),
  fDetElements(ioCtor)
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
Bool_t
AliMpDEStore::ReadDENames(AliMp::StationType station, 
                          AliMq::Station12Type station12)
{ 
/// Read det element names for cath = 0 from the file specified by name
/// and fill the map 

  // Open stream
  istream& in 
    = fkDataStreams.
        CreateDataStream(AliMpFiles::DENamesFilePath(station, station12));
  
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
      
      AliMpDetElement* detElement = new AliMpDetElement(detElemId, name, name0, planeForCathode[0]);
      
      if ( ! fDetElements.GetValue(detElemId) ) 
      {
        AliDebugClassStream(3)  
          << "Adding DE name "  << detElemId << "  " << name << endl;
        fDetElements.Add(detElemId, detElement); 
      } 
      else 
      {
        AliWarningClassStream()
          << "Det element "  << detElemId << "  " << name << " already defined." << endl;
      }	
    } 
    in >> word;
  }
  
  delete &in;

  return true;
}

//______________________________________________________________________________
void AliMpDEStore::FillDEs()
{
/// Fill DE names from files
  AliDebugClass(2,"");
  Bool_t result1 = ReadDENames(AliMp::kStation12, AliMq::kStation1);
  Bool_t result2 = ReadDENames(AliMp::kStation12, AliMq::kStation2);
  Bool_t result3 = ReadDENames(AliMp::kStation345);
  Bool_t result4 = ReadDENames(AliMp::kStationTrigger);
  
  Bool_t result = result1 && result2 && result3 && result4;
  if ( ! result ) {
    AliErrorClassStream() << "Error in reading DE names files" << endl;
  }  
  AliDebug(1,Form("%d detection elements were read in",fDetElements.GetSize()));
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
AliMpDetElement* AliMpDEStore::GetDetElement(const TString& deName, Bool_t warn) const
{
/// Return det element for given deName

  TIter next(fDetElements.CreateIterator());
  AliMpDetElement* detElement;
  
  while ( ( detElement = static_cast<AliMpDetElement*>(next()) ) )
  {
             
    if (deName.CompareTo(detElement->GetDEName()) == 0) 

      return detElement;
  }

  if (warn) {  
    AliErrorClassStream() 
	<< "Detection element with name" << deName.Data() << " not defined." << endl;
  }	

  return 0x0;   

}
