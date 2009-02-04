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
#include "AliMpStation12Type.h"
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
    
    static TString PlaneDataDir(AliMp::StationType station, 
                                AliMq::Station12Type station12Type,
                                AliMp::PlaneType plane);
    static TString StationDataDir(AliMp::StationType station,
                                AliMq::Station12Type station12Type);
  
    // bus patch
    //
    static TString BusPatchFilePath(); 
    static TString BusPatchInfoFilePath(); 
    static TString BusPatchSpecialFilePath(); 

    // de names
    //
    static TString DENamesFilePath(AliMp::StationType stationType,
                                   AliMq::Station12Type station12Type);

    // trigger
    //
    static TString LocalTriggerBoardMapping();
    static TString GlobalTriggerBoardMapping();
    
    // slats
    //
    static TString SlatFilePath(AliMp::StationType stationType, 
                                const char* slatType, AliMp::PlaneType plane);
    static TString SlatPCBFilePath(AliMp::StationType stationType, 
                                const char* pcbType);
    // sectors
    //
    static TString SectorFilePath(AliMq::Station12Type station, 
                                  AliMp::PlaneType plane);
    static TString SectorSpecialFilePath(AliMq::Station12Type station, 
                                  AliMp::PlaneType plane);
    static TString SectorSpecialFilePath2(AliMq::Station12Type station, 
                                  AliMp::PlaneType plane);
    // motifs
    //
    static TString MotifFilePath(AliMp::StationType station, 
                                 AliMq::Station12Type station12Type,
                                 AliMp::PlaneType plane, 
                                 const TString& motifTypeID);
    static TString MotifFileName(const TString& motifTypeID);
    static TString MotifSpecialFilePath(AliMp::StationType station,
                                 AliMq::Station12Type station12Type,
                                 AliMp::PlaneType plane, const TString& motifID);
    static TString MotifSpecialFileName(const TString& motifID);
    static TString PadPosFilePath(AliMp::StationType station, 
                                 AliMq::Station12Type station12Type,
                                 AliMp::PlaneType plane, const TString& motifTypeID);
    static TString PadPosFileName(const TString& motifTypeID);

    static TString BergToGCFilePath(AliMp::StationType station,
                                 AliMq::Station12Type station12Type);

    static TString ManuToSerialPath(const TString& deName, 
                                 AliMp::StationType station,
                                 AliMq::Station12Type station12Type);

    static TString SerialToBinFilePath();

    // set methods
    static void SetTopPath(const TString& topPath);
    static TString GetTop();
  
  private: 
    /// Not implemented
    AliMpFiles();
    /// Not implemented
    AliMpFiles(const AliMpFiles& right);
    /// Not implemented
    AliMpFiles& operator=(const AliMpFiles& right);    

    // static data members  
    static const TString fgkDataDir;       ///< data directory
    static const TString fgkDataRunDir;    ///< directory for run dependent data
    static const TString fgkStationDir;    ///< station directory
    static const TString fgkBendingDir;    ///< bending plane directory
    static const TString fgkNonBendingDir; ///< non-bending plane directory
    static const TString fgkDENames;       ///< DE names data file name
    static const TString fgkSector;        ///< sector data file name
    static const TString fgkSectorSpecial; ///< sector special data file name
    static const TString fgkSectorSpecial2;///< sector special data file name
    static const TString fgkMotifPrefix;   ///< motif data file name
    static const TString fgkMotifSpecialPrefix; ///< special motif data file name 
    static const TString fgkManuToSerial;  ///< manu to serial file name suffix
    static const TString fgkPadPosPrefix;  ///< pad position data file name
    static const TString fgkDataExt;       ///< file extension
    static const TString fgkBergToGCFileName;  ///< BergToGC mapping file name
    static const TString fgkTriggerLocalBoards;///<  local board name to id mapping
    static const TString fgkTriggerGlobalBoards;///<  global board name to id mapping
    static const TString fgkBusPatchFileName;  ///< DetElemIdToBusPatch file name
    static const TString fgkBusPatchInfoFileName;///< BusPatch length file name
    static const TString fgkBusPatchSpecialFileName;///< BusPatch special file name
    static const TString fgkSerialToBinFileName; ///< serial to bin  number file name
    
  ClassDef(AliMpFiles, 0) //File names and paths 
};  

#endif //ALI_MP_FILES_H
