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
// $MpId: AliMpFiles.cxx,v 1.4 2005/08/26 15:43:36 ivana Exp $
// Category: basic
//
// Class AliMpFiles
// ----------------
// Class for generating file names and paths.
// The input files:
// zones.dat, zones_special.dat - sector description
// motif*.dat   - motif description (generated from Exceed)
// padPos*.dat  - pad positions in motif
//
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <stdlib.h>

#include "AliMpFiles.h"
#include "AliLog.h"
#include "TClass.h"

ClassImp(AliMpFiles)

//
// static

// static data members

const TString AliMpFiles::fgkDefaultTop = GetDefaultTop();
const TString AliMpFiles::fgkDataDir = "/data";
const TString AliMpFiles::fgkStationDir = "/station";
const TString AliMpFiles::fgkBendingDir = "/bending_plane/";
const TString AliMpFiles::fgkNonBendingDir = "/non-bending_plane/";
const TString AliMpFiles::fgkSector  = "zones"; 
const TString AliMpFiles::fgkSectorSpecial = "zones_special";
const TString AliMpFiles::fgkSectorSpecial2 = "zones_special_outer";
const TString AliMpFiles::fgkMotifPrefix   = "motif";  
const TString AliMpFiles::fgkMotifSpecialPrefix ="motifSpecial";
const TString AliMpFiles::fgkPadPosPrefix  = "padPos"; 
const TString AliMpFiles::fgkDataExt = ".dat";      
const TString AliMpFiles::fgkBergToGCFileName = "/bergToGC"; 
const TString AliMpFiles::fgkTriggerLocalBoards = "MUONLocalTriggerBoard";

TString AliMpFiles::fgTop = AliMpFiles::fgkDefaultTop;

//______________________________________________________________________________
AliMpFiles::AliMpFiles()
  : TObject()
{
/// Default constructor
}
  
//______________________________________________________________________________
AliMpFiles::AliMpFiles(const AliMpFiles& right)
  : TObject(right) 
{
/// Protected copy constructor 

  AliFatalClass("Attempt to copy AliMpFiles singleton.");
}


//______________________________________________________________________________
AliMpFiles::~AliMpFiles() 
{
/// Destructor
}

// operators

//______________________________________________________________________________
AliMpFiles& AliMpFiles::operator=(const AliMpFiles& right)
{
/// Assignment operator

  // check assignment to self
  if (this == &right) return *this;

  AliFatalClass("Attempt to assign AliMpFiles singleton.");
    
  return *this;  
}    
          
//
// private methods
//

//______________________________________________________________________________
const char* AliMpFiles::GetDefaultTop()
{
  const char* top = getenv("MINSTALL");    
  if (!top)
  {
    const char* ntop = getenv("ALICE_ROOT");
    if (!ntop) return 0;
    TString dirPath(ntop);
    dirPath += "/MUON/mapping"; 
    return dirPath.Data();
  }
  return top;
}

//______________________________________________________________________________
TString AliMpFiles::PlaneDataDir(AliMpStationType station, 
                                 AliMpPlaneType plane)
{
/// Returns path to data files with sector description
/// for a specified plane.

  switch (station) {
  case kStation1:
  case kStation2:
    switch (plane) {
    case kBendingPlane:
      return fgTop + fgkDataDir + StationDataDir(station) + fgkBendingDir;
      ;;
    case kNonBendingPlane:   
      return fgTop + fgkDataDir + StationDataDir(station) + fgkNonBendingDir;
      ;;
    }   
    break;
  case kStation345:
  case kStationTrigger:  
    return fgTop + fgkDataDir + StationDataDir(station) + "/";
    break;
  default:  
    AliFatalClass("Incomplete switch on AliMpPlaneType");
    break;
  }
  return TString();
}

//______________________________________________________________________________
TString AliMpFiles::StationDataDir(AliMpStationType station)
{
/// Returns the station directory name for the specified station number.

  TString stationDataDir(fgkStationDir);
  switch (station) {
  case kStation1: 
    stationDataDir += 1;
    break;
    ;;
  case kStation2: 
    stationDataDir += 2;
    break;
    ;;
  case kStation345: 
    stationDataDir += "345/";
    break;
    ;;      
  case kStationTrigger:
    stationDataDir += "Trigger/";
    break;
    ;;
  default:
    stationDataDir += "Invalid/";
    break;
  }   
  return stationDataDir;
}

//
// public methods
//

//_____________________________________________________________________________
TString AliMpFiles::SlatFilePath(AliMpStationType stationType,
                                 const char* slatType,
                                 AliMpPlaneType plane)
{
/// \todo add ..

  return TString(PlaneDataDir(stationType,plane) + slatType + "." +
		 ( plane == kNonBendingPlane ? "NonBending":"Bending" ) + ".slat");
}

//_____________________________________________________________________________
TString AliMpFiles::SlatPCBFilePath(AliMpStationType stationType,
                                    const char* pcbType)
{
/// Get the full path for a given PCB (only relevant to stations 3,
/// 4, 5 and trigger). The bending parameter below is of no use in this case, but
/// we use it to re-use the PlaneDataDir() method untouched.

  return TString(PlaneDataDir(stationType,kNonBendingPlane) + pcbType +
                 ".pcb");
}

//______________________________________________________________________________
TString
AliMpFiles::LocalTriggerBoardMapping()
{
  return TString(PlaneDataDir(kStationTrigger,kNonBendingPlane) + fgkTriggerLocalBoards
                 + fgkDataExt);
}

//______________________________________________________________________________
TString AliMpFiles::SectorFilePath(AliMpStationType station, 
                                   AliMpPlaneType plane)
{
/// Return path to data file with sector description.
 
  return TString(PlaneDataDir(station, plane) + fgkSector + fgkDataExt);
}
    
//______________________________________________________________________________
TString AliMpFiles::SectorSpecialFilePath(AliMpStationType station, 
                                          AliMpPlaneType plane)
{
/// Return path to data file with sector special description (irregular motifs).

  return TString(PlaneDataDir(station, plane) + fgkSectorSpecial + fgkDataExt);
}
    
//______________________________________________________________________________
TString AliMpFiles::SectorSpecialFilePath2(AliMpStationType station, 
                                           AliMpPlaneType plane)
{
/// Returns path to data file with sector special description (irregular motifs).

  return TString(PlaneDataDir(station, plane) + fgkSectorSpecial2 + fgkDataExt);
}
    
//______________________________________________________________________________
TString AliMpFiles::MotifFilePath(AliMpStationType station, 
                                  AliMpPlaneType plane, 
                                  const TString& motifTypeID)
{
/// Returns path to data file for a given motif type.

  return TString(PlaneDataDir(station, plane) 
                 + fgkMotifPrefix +  motifTypeID + fgkDataExt);
}
    
//______________________________________________________________________________
TString AliMpFiles::PadPosFilePath(AliMpStationType station, 
                                   AliMpPlaneType plane, 
                                   const TString& motifTypeID)
{
/// Returns path to data file with pad positions for a given motif type.

  return TString(PlaneDataDir(station, plane) 
                 + fgkPadPosPrefix +  motifTypeID + fgkDataExt);
}

//______________________________________________________________________________ 
TString AliMpFiles::MotifSpecialFilePath(AliMpStationType station, 
                                         AliMpPlaneType plane,
                                         const TString& motifID)
{
/// Returns path to data file with pad dimensions for a given motif ID.

  return TString(PlaneDataDir(station, plane) 
                 + fgkMotifSpecialPrefix + motifID + fgkDataExt);

}

//______________________________________________________________________________ 
TString AliMpFiles::BergToGCFilePath(AliMpStationType station)
{
/// Returns the path of the file which describes the correspondance between
/// the berg number and the gassiplex channel.

  return fgTop + fgkDataDir + StationDataDir(station)
              + fgkBergToGCFileName + fgkDataExt;
}

//______________________________________________________________________________ 
void 
AliMpFiles::SetTopPath(const TString& topPath)
{ 
  fgTop = topPath; 
}

