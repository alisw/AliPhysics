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

ClassImp(AliMpFiles)

// static data members

AliMpFiles* AliMpFiles::fgInstance = 0;
const TString AliMpFiles::fgkDefaultTop = getenv("MINSTALL");    
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

//______________________________________________________________________________
AliMpFiles::AliMpFiles()
  : TObject(),
    fTop(fgkDefaultTop)
{
/// Default constructor
    
  if (fgInstance) {
    Fatal("AliMpFiles", 
          "AliMpFiles: attempt to create two instances of singleton.");
  }
      
  fgInstance = this;      
}
  
//______________________________________________________________________________
AliMpFiles::AliMpFiles(const AliMpFiles& right)
  : TObject(right) 
{
/// Protected copy constructor 

  Fatal("AliMpFiles", "Attempt to copy AliMpFiles singleton.");
}


//______________________________________________________________________________
AliMpFiles::~AliMpFiles() 
{
/// Destructor

  fgInstance = 0;      
}

// operators

//______________________________________________________________________________
AliMpFiles& AliMpFiles::operator=(const AliMpFiles& right)
{
/// Assignment operator

  // check assignment to self
  if (this == &right) return *this;

  Fatal("operator=", "Attempt to assign AliMpFiles singleton.");
    
  return *this;  
}    
          
//
// private methods
//

//______________________________________________________________________________
TString AliMpFiles::PlaneDataDir(AliMpStationType station, 
                                 AliMpPlaneType plane) const
{
/// Returns path to data files with sector description
/// for a specified plane.

  switch (station) {
  case kStation1:
  case kStation2:
    switch (plane) {
    case kBendingPlane:
      return fTop + fgkDataDir + StationDataDir(station) + fgkBendingDir;
      ;;
    case kNonBendingPlane:   
      return fTop + fgkDataDir + StationDataDir(station) + fgkNonBendingDir;
      ;;
    }   
    break;
  case kStation345:
    return fTop + fgkDataDir + StationDataDir(station) + "/";
    break;
  }

  Fatal("PlaneDataDir", "Incomplete switch on AliMpPlaneType");
  return TString();
}

//______________________________________________________________________________
TString AliMpFiles::StationDataDir(AliMpStationType station) const
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
  }   
  return stationDataDir;
}

//
// public methods
//

//______________________________________________________________________________
AliMpFiles* AliMpFiles::Instance() 
{ 
/// Return the singleton instance;
/// Creates it if it does not yet exist,

  if (!fgInstance)  fgInstance = new AliMpFiles(); 
  
  return fgInstance; 
}


//_____________________________________________________________________________
TString AliMpFiles::SlatFilePath(const char* slatType,
				 AliMpPlaneType plane) const
{
/// \todo add ..

  return TString(PlaneDataDir(kStation345,plane) + slatType + "." +
		 ( plane == kNonBendingPlane ? "NonBending":"Bending" ) + ".slat");
}

//_____________________________________________________________________________
TString AliMpFiles::SlatPCBFilePath(const char* pcbType) const
{
/// Get the full path for a given PCB (only relevant to stations 3,
/// 4 and 5). The bending parameter below is of no use in this case, but
/// we use it to re-use the PlaneDataDir() method untouched.

  return TString(PlaneDataDir(kStation345,kNonBendingPlane) + pcbType +
		 ".pcb");
}

//_____________________________________________________________________________
TString 
AliMpFiles::DetElemIdToSlatTypeFilePath() const
{
/// Get the full path of the file containing the mapping detElemId <->
/// SlatType.
/// The bending parameter below is of no use in this case, but
/// we use it to re-use the PlaneDataDir() method untouched.

  return TString(PlaneDataDir(kStation345,kNonBendingPlane) + 
		 "DetElemIdToSlatType.dat");
}
//______________________________________________________________________________
TString AliMpFiles::SectorFilePath(AliMpStationType station, 
                                   AliMpPlaneType plane) const
{
/// Return path to data file with sector description.
 
  return TString(PlaneDataDir(station, plane) + fgkSector + fgkDataExt);
}
    
//______________________________________________________________________________
TString AliMpFiles::SectorSpecialFilePath(AliMpStationType station, 
                                          AliMpPlaneType plane) const
{
/// Return path to data file with sector special description (irregular motifs).

  return TString(PlaneDataDir(station, plane) + fgkSectorSpecial + fgkDataExt);
}
    
//______________________________________________________________________________
TString AliMpFiles::SectorSpecialFilePath2(AliMpStationType station, 
                                           AliMpPlaneType plane) const
{
/// Returns path to data file with sector special description (irregular motifs).

  return TString(PlaneDataDir(station, plane) + fgkSectorSpecial2 + fgkDataExt);
}
    
//______________________________________________________________________________
TString AliMpFiles::MotifFilePath(AliMpStationType station, 
                                  AliMpPlaneType plane, 
                                  const TString& motifTypeID) const
{
/// Returns path to data file for a given motif type.

  return TString(PlaneDataDir(station, plane) 
                 + fgkMotifPrefix +  motifTypeID + fgkDataExt);
}
    
//______________________________________________________________________________
TString AliMpFiles::PadPosFilePath(AliMpStationType station, 
                                   AliMpPlaneType plane, 
                                   const TString& motifTypeID) const
{
/// Returns path to data file with pad positions for a given motif type.

  return TString(PlaneDataDir(station, plane) 
                 + fgkPadPosPrefix +  motifTypeID + fgkDataExt);
}

//______________________________________________________________________________ 
TString AliMpFiles::MotifSpecialFilePath(AliMpStationType station, 
                                         AliMpPlaneType plane,
                                         const TString& motifID) const
{
/// Returns path to data file with pad dimensions for a given motif ID.

  return TString(PlaneDataDir(station, plane) 
                 + fgkMotifSpecialPrefix + motifID + fgkDataExt);

}

//______________________________________________________________________________ 
TString AliMpFiles::BergToGCFilePath(AliMpStationType station) const
{
/// Returns the path of the file which describes the correspondance between
/// the berg number and the gassiplex channel.

  return fTop + fgkDataDir + StationDataDir(station)
              + fgkBergToGCFileName + fgkDataExt;
}
