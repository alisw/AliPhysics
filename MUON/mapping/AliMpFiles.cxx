// $Id$
// Category: sector
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
//    
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
// 
  Fatal("AliMpFiles", "Attempt to copy AliMpFiles singleton.");
}


//______________________________________________________________________________
AliMpFiles::~AliMpFiles() {
//

  fgInstance = 0;      
}

// operators

//______________________________________________________________________________
AliMpFiles& AliMpFiles::operator=(const AliMpFiles& right)
{
  // check assignement to self
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
// Returns path to data files with sector description
// for a specified plane.
// ---

  switch (plane) {
    case kBendingPlane:
       return fTop + fgkDataDir + StationDataDir(station) + fgkBendingDir;
       ;;
    case kNonBendingPlane:   
       return fTop + fgkDataDir + StationDataDir(station) + fgkNonBendingDir;
       ;;
  }   
  
  Fatal("PlaneDataDir", "Incomplete switch on AliMpPlaneType");
  return TString();
}

//______________________________________________________________________________
TString AliMpFiles::StationDataDir(AliMpStationType station) const
{
// Returns the station directory name for the specified station number.
// ---

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
  }   
  return stationDataDir;
}

//
// public methods
//

//______________________________________________________________________________
AliMpFiles* AliMpFiles::Instance() 
{ 
// Return the singleton instance;
// Creates it if it does not yet exist,
//
// ---

  if (!fgInstance)  fgInstance = new AliMpFiles(); 
  
  return fgInstance; 
}

//______________________________________________________________________________
TString AliMpFiles::SectorFilePath(AliMpStationType station, 
                                   AliMpPlaneType plane) const
{
// Returns path to data file with sector description.
// ---
 
  return TString(PlaneDataDir(station, plane) + fgkSector + fgkDataExt);
}
    
//______________________________________________________________________________
TString AliMpFiles::SectorSpecialFilePath(AliMpStationType station, 
                                          AliMpPlaneType plane) const
{
// Returns path to data file with sector special description (irregular motifs).
// ---

  return TString(PlaneDataDir(station, plane) + fgkSectorSpecial + fgkDataExt);
}
    
//______________________________________________________________________________
TString AliMpFiles::SectorSpecialFilePath2(AliMpStationType station, 
                                           AliMpPlaneType plane) const
{
// Returns path to data file with sector special description (irregular motifs).
// ---

  return TString(PlaneDataDir(station, plane) + fgkSectorSpecial2 + fgkDataExt);
}
    
//______________________________________________________________________________
TString AliMpFiles::MotifFilePath(AliMpStationType station, 
                                  AliMpPlaneType plane, 
                                  const TString& motifTypeID) const
{
// Returns path to data file for a given motif type.
// ---

  return TString(PlaneDataDir(station, plane) 
                 + fgkMotifPrefix +  motifTypeID + fgkDataExt);
}
    
//______________________________________________________________________________
TString AliMpFiles::PadPosFilePath(AliMpStationType station, 
                                   AliMpPlaneType plane, 
                                   const TString& motifTypeID) const
{
// Returns path to data file with pad positions for a given motif type.
// ---

  return TString(PlaneDataDir(station, plane) 
                 + fgkPadPosPrefix +  motifTypeID + fgkDataExt);
}

//______________________________________________________________________________ 
TString AliMpFiles::MotifSpecialFilePath(AliMpStationType station, 
                                         AliMpPlaneType plane,
                                         const TString& motifID) const
{
// Returns path to data file with pad dimensions for a given motif ID.
// ---

  return TString(PlaneDataDir(station, plane) 
                 + fgkMotifSpecialPrefix + motifID + fgkDataExt);

}

//______________________________________________________________________________ 
TString AliMpFiles::BergToGCFilePath(AliMpStationType station) const
{
// Returns the path of the file which describes the correspondance between
// the berg number and the gassiplex channel.
// ---

  return fTop + fgkDataDir + StationDataDir(station)
              + fgkBergToGCFileName + fgkDataExt;
}
