/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpFiles.h,v 1.4 2005/08/26 15:43:36 ivana Exp $

/// \ingroup basic
/// \class AliMpFiles
/// \brief Class for generating file names and paths.
///
/// The input files:
/// - zones.dat, zones_special.dat - sector description
/// - motif*.dat   - motif description (generated from Exceed)
/// - padPos*.dat  - pad positions in motif
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_FILES_H
#define ALI_MP_FILES_H

#include <TObject.h>
#include <TString.h>

#include "AliMpStationType.h"
#include "AliMpPlaneType.h"

class AliMpFiles : public TObject
{
  public:
    AliMpFiles();
    // --> protected
    //AliMpFiles(const AliMpFiles& right);
    virtual ~AliMpFiles();
    
    // static access method
    static AliMpFiles* Instance();

    // methods
    TString SlatFilePath(const char* slatType,
			 AliMpPlaneType plane) const;
    TString SlatPCBFilePath(const char* pcbType) const;

    TString DetElemIdToSlatTypeFilePath() const;

    TString SectorFilePath(AliMpStationType station, 
                           AliMpPlaneType plane) const;
    TString SectorSpecialFilePath(AliMpStationType station, 
                           AliMpPlaneType plane) const;
    TString SectorSpecialFilePath2(AliMpStationType station, 
                           AliMpPlaneType plane) const;
    TString MotifFilePath(AliMpStationType station, AliMpPlaneType plane, 
                          const TString& motifTypeID) const;
    TString MotifSpecialFilePath(AliMpStationType station, 
                          AliMpPlaneType plane, const TString& motifID) const;
    TString PadPosFilePath(AliMpStationType station, 
                          AliMpPlaneType plane, const TString& motifTypeID) const;
    TString BergToGCFilePath(AliMpStationType station) const;

    // set methods
    void SetTopPath(const TString& topPath);
			      
  protected:
    AliMpFiles(const AliMpFiles& right);
    
    // operators
    AliMpFiles& operator=(const AliMpFiles& right);    
    
  private: 
    // methods
    TString PlaneDataDir(AliMpStationType station, AliMpPlaneType plane) const; 
    TString StationDataDir(AliMpStationType station) const; 
  
    // static data members  
    static AliMpFiles*   fgInstance;       //this instance
    static const TString fgkDefaultTop;    //top directory path (default)
    static const TString fgkDataDir;       //data directory
    static const TString fgkStationDir;    //station directory
    static const TString fgkBendingDir;    //bending plane directory
    static const TString fgkNonBendingDir; //non-bending plane directory
    static const TString fgkSector;        //sector data file name
    static const TString fgkSectorSpecial; //sector special data file name
    static const TString fgkSectorSpecial2;//sector special data file name
    static const TString fgkMotifPrefix;   //motif data file name
    static const TString fgkMotifSpecialPrefix; //special motif data file name 
    static const TString fgkPadPosPrefix;  //pad position data file name
    static const TString fgkDataExt;       //file extension
    static const TString fgkBergToGCFileName; //name of the file with BergToGC mapping
    
    // data members
    TString  fTop; // top directory path
    

  ClassDef(AliMpFiles, 1) //File names and paths 
};  

// inline functions

inline void AliMpFiles::SetTopPath(const TString& topPath)
{ fTop = topPath; }

#endif //ALI_MP_FILES_H
