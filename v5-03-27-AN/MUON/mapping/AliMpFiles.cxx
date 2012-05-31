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

//-----------------------------------------------------------------------------
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
//-----------------------------------------------------------------------------

#include "AliMpFiles.h"

#include "AliLog.h"

#include <TClass.h>
#include <Riostream.h>

#include <stdlib.h>

/// \cond CLASSIMP
ClassImp(AliMpFiles)
/// \endcond

//
// static private methods
//

//______________________________________________________________________________
const TString& AliMpFiles::GetDataDir() 
{
  /// data directory
  static const TString kDataDir = "/data";
  return kDataDir;
}  
  
//______________________________________________________________________________
const TString& AliMpFiles::GetDataRunDir() 
{
  /// directory for run dependent data
  static const TString kDataRunDir = "/data_run";
  return kDataRunDir;
}

//______________________________________________________________________________
const TString& AliMpFiles::GetStationDir() 
{
  /// station directory
  static const TString kStationDir = "/station";
  return kStationDir;
}  
  
//______________________________________________________________________________
const TString& AliMpFiles::GetBendingDir() 
{
  /// bending plane directory
  static const TString kBendingDir = "bending_plane/";
  return kBendingDir;
}  
  
//______________________________________________________________________________
const TString& AliMpFiles::GetNonBendingDir() 
{
  /// non-bending plane directory
  static const TString kNonBendingDir = "non-bending_plane/";
  return kNonBendingDir;
}  

//______________________________________________________________________________
const TString& AliMpFiles::GetDENames() 
{
  /// DE names data file name
  static const TString kDENames = "denames"; 
  return kDENames;
}  
  
//______________________________________________________________________________
const TString& AliMpFiles::GetSector() 
{
  /// sector data file name
  static const TString kSector  = "zones"; 
  return kSector;
}  
  
//______________________________________________________________________________
const TString& AliMpFiles::GetSectorSpecial() 
{
  /// sector special data file name
  static const TString kSectorSpecial = "zones_special";
  return kSectorSpecial;
}  
  
//______________________________________________________________________________
const TString& AliMpFiles::GetSectorSpecial2() 
{
  /// sector special data file name
  static const TString kSectorSpecial2 = "zones_special_outer";
  return kSectorSpecial2;
}  

//______________________________________________________________________________
const TString& AliMpFiles::GetMotifPrefix() 
{
  /// motif data file name
  static const TString kMotifPrefix   = "motif"; 
  return kMotifPrefix;
}   
    

//______________________________________________________________________________
const TString& AliMpFiles::GetMotifSpecialPrefix() 
{
  /// special motif data file name 
  static const TString kMotifSpecialPrefix ="motifSpecial";
  return kMotifSpecialPrefix;
}  

//______________________________________________________________________________
const TString& AliMpFiles::GetManuToSerial() 
{
  /// manu to serial file name suffix
  static const TString kManuToSerial ="_manu";
  return kManuToSerial;
}  

//______________________________________________________________________________
const TString& AliMpFiles::GetPadPosPrefix() 
{  
  /// pad position data file name
  static const TString kPadPosPrefix  = "padPos"; 
  return kPadPosPrefix;
}  

//______________________________________________________________________________
const TString& AliMpFiles::GetDataExt() 
{
  /// file extension
  static const TString kDataExt = ".dat"; 
  return kDataExt;
}
       
//______________________________________________________________________________
const TString& AliMpFiles::GetBergToGCFileName() 
{
  /// BergToGC mapping file name
  static const TString kBergToGCFileName = "bergToGC"; 
  return kBergToGCFileName;
}  

//______________________________________________________________________________
const TString& AliMpFiles::GetTriggerLocalBoards() 
{
  ///  local board name to id mapping
  static const TString kTriggerLocalBoards = "RegionalCrate";
  return kTriggerLocalBoards;
}  

//______________________________________________________________________________
const TString& AliMpFiles::GetTriggerGlobalBoards() 
{
  ///  global board name to id mapping
  static const TString kTriggerGlobalBoards = "GlobalCrate";
  return kTriggerGlobalBoards;
}  
  
//______________________________________________________________________________
const TString& AliMpFiles::GetBusPatchFileName() 
{
  /// DetElemIdToBusPatch file name
  static const TString kBusPatchFileName = "DetElemIdToBusPatch";
  return kBusPatchFileName;
}  
  
//______________________________________________________________________________
const TString& AliMpFiles::GetBusPatchInfoFileName() 
{
  /// BusPatch length file name
  static const TString kBusPatchInfoFileName = "BusPatchInfo";
  return kBusPatchInfoFileName;
}  
  
//______________________________________________________________________________
const TString& AliMpFiles::GetBusPatchSpecialFileName() 
{
  /// BusPatch special file name
  static const TString kBusPatchSpecialFileName = "BusPatchSpecial";
  return kBusPatchSpecialFileName;
}  
  
//______________________________________________________________________________
const TString& AliMpFiles::GetSerialToBinFileName() 
{
  /// serial to bin  number file name
  static const TString kSerialToBinFileName = "ManuSerialToBin";
  return kSerialToBinFileName;
}  
  
//
// static public methods
//

//______________________________________________________________________________
TString AliMpFiles::PlaneDataDir(AliMp::StationType station, 
                                 AliMq::Station12Type station12Type,
                                 AliMp::PlaneType plane)
{
/// Returns path to data files with sector description
/// for a specified plane.

  switch (station) {
  case AliMp::kStation12:
    switch (plane) {
    case AliMp::kBendingPlane:
      return GetTop() + GetDataDir() + StationDataDir(station, station12Type) + GetBendingDir();
      ;;
    case AliMp::kNonBendingPlane:   
      return GetTop() + GetDataDir() + StationDataDir(station, station12Type) + GetNonBendingDir();
      ;;
    }   
    break;
  case AliMp::kStation345:
  case AliMp::kStationTrigger:  
    return GetTop() + GetDataDir() + StationDataDir(station, AliMq::kNotSt12);
    break;
  default:  
    AliFatalClass("Incomplete switch on AliMp::PlaneType");
    break;
  }
  return TString();
}

//______________________________________________________________________________
TString AliMpFiles::StationDataDir(AliMp::StationType station,
                                   AliMq::Station12Type station12Type)
{
/// Returns the station directory name for the specified station number.

  TString stationDataDir(GetStationDir());
  switch (station) {
  case AliMp::kStation12: 
    switch (station12Type) {
    case AliMq::kStation1:
      stationDataDir += "1/";
      break;
      ;;
    case AliMq::kStation2:   
      stationDataDir += "2/";
      break;
      ;;
    case AliMq::kNotSt12:   
      AliFatalClass("Incorrect switch on AliMq::Station12Type");
      break;
    }   
    break;
    ;;
  case AliMp::kStation345: 
    stationDataDir += "345/";
    break;
    ;;      
  case AliMp::kStationTrigger:
    stationDataDir += "Trigger/";
    break;
    ;;
  default:
    stationDataDir += "Invalid/";
    break;
  }   
  return stationDataDir;
}

//______________________________________________________________________________
TString AliMpFiles::BusPatchFilePath()
{
/// Return path to data file with bus patch mapping.

  return GetTop() + GetDataDir() + "/" + GetBusPatchFileName() + GetDataExt();
}  

//______________________________________________________________________________
TString AliMpFiles::BusPatchInfoFilePath()
{
/// Return path to data file with bus patch mapping.

  return GetTop() + GetDataDir() + "/" + GetBusPatchInfoFileName() + GetDataExt();
}  

//______________________________________________________________________________
TString AliMpFiles::BusPatchSpecialFilePath()
{
/// Return path to data file with special bus patch mapping.

  return GetTop() + GetDataDir() + "/" + GetBusPatchSpecialFileName() + GetDataExt();
}  

//______________________________________________________________________________
TString AliMpFiles::SerialToBinFilePath()
{
/// Return path to data file containing manu serial numbers with their bin.

  return GetTop() + GetDataDir() + "/" + GetSerialToBinFileName() + GetDataExt();
}  


//______________________________________________________________________________
TString AliMpFiles::DENamesFilePath(AliMp::StationType station,
                                    AliMq::Station12Type station12Type)
{
/// Return path to data file with DE names for given station.
 
  return GetTop() + GetDataDir() + StationDataDir(station, station12Type) 
           + GetDENames() + GetDataExt();
}

//______________________________________________________________________________
TString AliMpFiles::LocalTriggerBoardMapping()
{
/// Return path to data file with local trigger board mapping.

  return GetTop() + GetDataDir() 
          + StationDataDir(AliMp::kStationTrigger, AliMq::kNotSt12) 
          + GetTriggerLocalBoards() + GetDataExt();;
}

//______________________________________________________________________________
TString AliMpFiles::GlobalTriggerBoardMapping()
{
/// Return path to data file with local trigger board mapping.

  return GetTop() + GetDataDir() 
      + StationDataDir(AliMp::kStationTrigger, AliMq::kNotSt12) 
      + GetTriggerGlobalBoards() + GetDataExt();;
}

//_____________________________________________________________________________
TString AliMpFiles::SlatFilePath(AliMp::StationType stationType,
                                 const char* slatType,
                                 AliMp::PlaneType plane)
{
/// \todo add ..

  return TString(PlaneDataDir(stationType, AliMq::kNotSt12, plane) 
                 + slatType + "." +
		 ( plane == AliMp::kNonBendingPlane ? "NonBending":"Bending" ) + ".slat");
}

//_____________________________________________________________________________
TString AliMpFiles::SlatPCBFilePath(AliMp::StationType stationType,
                                    const char* pcbType)
{
/// Get the full path for a given PCB (only relevant to stations 3,
/// 4, 5 and trigger). The bending parameter below is of no use in this case, but
/// we use it to re-use the PlaneDataDir() method untouched.

  return TString(PlaneDataDir(stationType, AliMq::kNotSt12, AliMp::kNonBendingPlane) 
                 + pcbType + ".pcb");
}

//______________________________________________________________________________
TString AliMpFiles::SectorFilePath(AliMq::Station12Type station12Type, 
                                   AliMp::PlaneType plane)
{
/// Return path to data file with sector description.
 
  return TString(PlaneDataDir(AliMp::kStation12, station12Type, plane) 
                 + GetSector() + GetDataExt());
}
    
//______________________________________________________________________________
TString AliMpFiles::SectorSpecialFilePath(AliMq::Station12Type station12Type,
                                          AliMp::PlaneType plane)
{
/// Return path to data file with sector special description (irregular motifs).

  return TString(PlaneDataDir(AliMp::kStation12, station12Type, plane) 
                 + GetSectorSpecial() + GetDataExt());
}
    
//______________________________________________________________________________
TString AliMpFiles::SectorSpecialFilePath2(AliMq::Station12Type station12Type, 
                                           AliMp::PlaneType plane)
{
/// Returns path to data file with sector special description (irregular motifs).

  return TString(PlaneDataDir(AliMp::kStation12, station12Type, plane) 
                 + GetSectorSpecial2() + GetDataExt());
}

//______________________________________________________________________________
TString AliMpFiles::MotifFileName(const TString& motifTypeID)
{
  /// Returns name of data file for a given motif type.
  
  return TString(GetMotifPrefix() +  motifTypeID + GetDataExt());
}

//______________________________________________________________________________
TString AliMpFiles::MotifFilePath(AliMp::StationType station, 
                                  AliMq::Station12Type station12Type,
                                  AliMp::PlaneType plane, 
                                  const TString& motifTypeID)
{
/// Returns path to data file for a given motif type.

  return TString(PlaneDataDir(station, station12Type, plane) 
                 + MotifFileName(motifTypeID));
}

//______________________________________________________________________________
TString AliMpFiles::PadPosFileName(const TString& motifTypeID)
{
  /// Returns name of data file with pad positions for a given motif type.
  
  return TString(GetPadPosPrefix() +  motifTypeID + GetDataExt());
}

//______________________________________________________________________________
TString AliMpFiles::PadPosFilePath(AliMp::StationType station, 
                                   AliMq::Station12Type station12Type,
                                   AliMp::PlaneType plane, 
                                   const TString& motifTypeID)
{
/// Returns path to data file with pad positions for a given motif type.

  return TString(PlaneDataDir(station, station12Type, plane) 
                 + PadPosFileName(motifTypeID));
}

//______________________________________________________________________________ 
TString AliMpFiles::MotifSpecialFileName(const TString& motifID)
{
  /// Returns name of data file with pad dimensions for a given motif ID.
  
  return TString(GetMotifSpecialPrefix() + motifID + GetDataExt());
  
}

//______________________________________________________________________________ 
TString AliMpFiles::MotifSpecialFilePath(AliMp::StationType station, 
                                         AliMq::Station12Type station12Type,
                                         AliMp::PlaneType plane,
                                         const TString& motifID)
{
/// Returns path to data file with pad dimensions for a given motif ID.

  return TString(PlaneDataDir(station, station12Type, plane) 
                 + MotifSpecialFileName(motifID));
}

//______________________________________________________________________________ 
TString AliMpFiles::BergToGCFilePath(AliMp::StationType station,
                                     AliMq::Station12Type station12Type)
{
/// Returns the path of the file which describes the correspondance between
/// the berg number and the gassiplex channel.

  return GetTop() + GetDataDir() + StationDataDir(station, station12Type)
              + GetBergToGCFileName() + GetDataExt();
}

//______________________________________________________________________________ 
TString AliMpFiles::ManuToSerialPath(const TString& deName, 
                                     AliMp::StationType station,
                                     AliMq::Station12Type station12Type)
{
/// Returns the path of the file for the manu id to their serial number

  return  GetTop() + GetDataRunDir() + StationDataDir(station, station12Type)
              + deName + GetManuToSerial() + GetDataExt(); 
}


//______________________________________________________________________________ 
void 
AliMpFiles::SetTopPath(const TString& topPath)
{ 
/// Set top file path

  GetTop() = topPath; 
}

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

//
// ctors, dtor
//

//______________________________________________________________________________
AliMpFiles::~AliMpFiles() 
{
/// Destructor
}

