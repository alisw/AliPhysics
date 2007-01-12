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
// $MpId: AliMpFiles.cxx,v 1.12 2006/05/23 13:09:54 ivana Exp $
// Category: basic
// ----------------
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

#include "AliMpFiles.h"

#include "AliLog.h"

#include <TClass.h>
#include <Riostream.h>

#include <stdlib.h>

/// \cond CLASSIMP
ClassImp(AliMpFiles)
/// \endcond

//
// static

// static data members

const TString AliMpFiles::fgkDataDir = "/data";
const TString AliMpFiles::fgkStationDir = "/station";
const TString AliMpFiles::fgkBendingDir = "/bending_plane/";
const TString AliMpFiles::fgkNonBendingDir = "/non-bending_plane/";
const TString AliMpFiles::fgkDENames = "denames"; 
const TString AliMpFiles::fgkSector  = "zones"; 
const TString AliMpFiles::fgkSectorSpecial = "zones_special";
const TString AliMpFiles::fgkSectorSpecial2 = "zones_special_outer";
const TString AliMpFiles::fgkMotifPrefix   = "motif";  
const TString AliMpFiles::fgkMotifSpecialPrefix ="motifSpecial";
const TString AliMpFiles::fgkManuToSerialDir ="manu_serial/";
const TString AliMpFiles::fgkManuToSerial ="_manu";
const TString AliMpFiles::fgkPadPosPrefix  = "padPos"; 
const TString AliMpFiles::fgkDataExt = ".dat";      
const TString AliMpFiles::fgkBergToGCFileName = "/bergToGC"; 
const TString AliMpFiles::fgkTriggerLocalBoards = "crate";
const TString AliMpFiles::fgkBusPatchFileName = "DetElemIdToBusPatch";

//______________________________________________________________________________
AliMpFiles::~AliMpFiles() 
{
/// Destructor
}

//
// private methods
//

//______________________________________________________________________________
TString AliMpFiles::GetTop()
{
/// Return top path to mapping data defined either via MINSTALL
/// or ALICE_ROOT environment variable.                                      \n
/// If both variables are defined, MINSTALL is used.

  TString top = getenv("MINSTALL");    
  if ( ! top.IsNull() ) return top;

  TString ntop = getenv("ALICE_ROOT");
  if ( ntop.IsNull() ) {
    AliErrorClassStream() << "Cannot find path to mapping data." << endl;
    return ntop;
  }  
  ntop += "/MUON/mapping";
  return ntop;
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
      return GetTop() + fgkDataDir + StationDataDir(station) + fgkBendingDir;
      ;;
    case kNonBendingPlane:   
      return GetTop() + fgkDataDir + StationDataDir(station) + fgkNonBendingDir;
      ;;
    }   
    break;
  case kStation345:
  case kStationTrigger:  
    return GetTop() + fgkDataDir + StationDataDir(station) + "/";
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
    stationDataDir += "1/";
    break;
    ;;
  case kStation2: 
    stationDataDir += "2/";
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

//______________________________________________________________________________
TString AliMpFiles::BusPatchFilePath()
{
/// Return path to data file with bus patch mapping.

  return GetTop() + fgkDataDir + "/" + fgkBusPatchFileName + fgkDataExt;
}  

//______________________________________________________________________________
TString AliMpFiles::DENamesFilePath(AliMpStationType station)
{
/// Return path to data file with DE names for given station.
 
  return GetTop() + fgkDataDir + StationDataDir(station) + fgkDENames + fgkDataExt;
}

//______________________________________________________________________________
TString AliMpFiles::LocalTriggerBoardMapping()
{
/// Return path to data file with local trigger board mapping.

  return GetTop() + fgkDataDir + StationDataDir(kStationTrigger) 
          + fgkTriggerLocalBoards + fgkDataExt;;
}

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
TString AliMpFiles::MotifFileName(const TString& motifTypeID)
{
  /// Returns name of data file for a given motif type.
  
  return TString(fgkMotifPrefix +  motifTypeID + fgkDataExt);
}

//______________________________________________________________________________
TString AliMpFiles::MotifFilePath(AliMpStationType station, 
                                  AliMpPlaneType plane, 
                                  const TString& motifTypeID)
{
/// Returns path to data file for a given motif type.

  return TString(PlaneDataDir(station, plane) + MotifFileName(motifTypeID));
}

//______________________________________________________________________________
TString AliMpFiles::PadPosFileName(const TString& motifTypeID)
{
  /// Returns name of data file with pad positions for a given motif type.
  
  return TString(fgkPadPosPrefix +  motifTypeID + fgkDataExt);
}

//______________________________________________________________________________
TString AliMpFiles::PadPosFilePath(AliMpStationType station, 
                                   AliMpPlaneType plane, 
                                   const TString& motifTypeID)
{
/// Returns path to data file with pad positions for a given motif type.

  return TString(PlaneDataDir(station, plane) + PadPosFileName(motifTypeID));
}

//______________________________________________________________________________ 
TString AliMpFiles::MotifSpecialFileName(const TString& motifID)
{
  /// Returns name of data file with pad dimensions for a given motif ID.
  
  return TString(fgkMotifSpecialPrefix + motifID + fgkDataExt);
  
}

//______________________________________________________________________________ 
TString AliMpFiles::MotifSpecialFilePath(AliMpStationType station, 
                                         AliMpPlaneType plane,
                                         const TString& motifID)
{
/// Returns path to data file with pad dimensions for a given motif ID.

  return TString(PlaneDataDir(station, plane) + MotifSpecialFileName(motifID));
}

//______________________________________________________________________________ 
TString AliMpFiles::BergToGCFilePath(AliMpStationType station)
{
/// Returns the path of the file which describes the correspondance between
/// the berg number and the gassiplex channel.

  return GetTop() + fgkDataDir + StationDataDir(station)
              + fgkBergToGCFileName + fgkDataExt;
}

//______________________________________________________________________________ 
TString AliMpFiles::ManuToSerialPath(const TString& deName, AliMpStationType station)
{
/// Returns the path of the file for the manu id to their serial number

  return  GetTop() + fgkDataDir + StationDataDir(station)
              + fgkManuToSerialDir + deName + fgkManuToSerial + fgkDataExt; 
}


//______________________________________________________________________________ 
void 
AliMpFiles::SetTopPath(const TString& topPath)
{ 
  GetTop() = topPath; 
}

