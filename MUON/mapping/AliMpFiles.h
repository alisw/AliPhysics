/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpFiles.h,v 1.10 2006/05/24 13:58:07 ivana Exp $

/// \ingroup basic
/// \class AliMpFiles
/// \brief Class for generating file names and paths.
///
/// The input files:
/// - zones.dat, zones_special.dat - sector description
/// - motif*.dat   - motif description (generated from Exceed)
/// - padPos*.dat  - pad positions in motif
///
/// \author David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_FILES_H
#define ALI_MP_FILES_H

#include <TObject.h>

#include "AliMpStationType.h"
#include "AliMpPlaneType.h"

#include <TString.h>

class AliMpFiles : public TObject
{
  public:
    // --> protected
    //AliMpFiles();
    //AliMpFiles(const AliMpFiles& right);
    virtual ~AliMpFiles();
  
    //
    // methods
    //
    
    // bus patch
    //
    static TString BusPatchFilePath(); 

    // de names
    //
    static TString DENamesFilePath(AliMp::StationType stationType);

    // trigger
    //
    static TString LocalTriggerBoardMapping();
  
    // slats
    //
    static TString SlatFilePath(AliMp::StationType stationType, 
                                const char* slatType, AliMp::PlaneType plane);
    static TString SlatPCBFilePath(AliMp::StationType stationType, 
                                const char* pcbType);
    // sectors
    //
    static TString SectorFilePath(AliMp::StationType station, 
                                  AliMp::PlaneType plane);
    static TString SectorSpecialFilePath(AliMp::StationType station, 
                                  AliMp::PlaneType plane);
    static TString SectorSpecialFilePath2(AliMp::StationType station, 
                                  AliMp::PlaneType plane);
    // motifs
    //
    static TString MotifFilePath(AliMp::StationType station, 
                                 AliMp::PlaneType plane, 
                                 const TString& motifTypeID);
    static TString MotifFileName(const TString& motifTypeID);
    static TString MotifSpecialFilePath(AliMp::StationType station,
                                 AliMp::PlaneType plane, const TString& motifID);
    static TString MotifSpecialFileName(const TString& motifID);
    static TString PadPosFilePath(AliMp::StationType station, 
                                 AliMp::PlaneType plane, const TString& motifTypeID);
    static TString PadPosFileName(const TString& motifTypeID);

    static TString BergToGCFilePath(AliMp::StationType station);

    static TString ManuToSerialPath(const TString& deName, AliMp::StationType station);

  
    // set methods
    static void SetTopPath(const TString& topPath);
  
  private: 
    AliMpFiles();
    AliMpFiles(const AliMpFiles& right);
  
    // operators
    AliMpFiles& operator=(const AliMpFiles& right);    
    // methods
    static TString GetTop();
    static TString PlaneDataDir(AliMp::StationType station, AliMp::PlaneType plane); 
    static TString StationDataDir(AliMp::StationType station); 
  
    // static data members  
    static const TString fgkDataDir;       ///< data directory
    static const TString fgkStationDir;    ///< station directory
    static const TString fgkBendingDir;    ///< bending plane directory
    static const TString fgkNonBendingDir; ///< non-bending plane directory
    static const TString fgkDENames;       ///< DE names data file name
    static const TString fgkSector;        ///< sector data file name
    static const TString fgkSectorSpecial; ///< sector special data file name
    static const TString fgkSectorSpecial2;///< sector special data file name
    static const TString fgkMotifPrefix;   ///< motif data file name
    static const TString fgkMotifSpecialPrefix; ///< special motif data file name 
    static const TString fgkManuToSerialDir;///< manu to serial file directory
    static const TString fgkManuToSerial;  ///< manu to serial file name suffix
    static const TString fgkPadPosPrefix;  ///< pad position data file name
    static const TString fgkDataExt;       ///< file extension
    static const TString fgkBergToGCFileName;  ///< BergToGC mapping file name
    static const TString fgkTriggerLocalBoards;///<  local board name to id mapping
    static const TString fgkBusPatchFileName;  ///< DetElemIdToBusPatch file name

  ClassDef(AliMpFiles, 0) //File names and paths 
};  

#endif //ALI_MP_FILES_H
