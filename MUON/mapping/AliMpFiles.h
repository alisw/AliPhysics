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
    // --> protected
    //AliMpFiles();
    //AliMpFiles(const AliMpFiles& right);
    virtual ~AliMpFiles();
  
    //
    // methods
    //
    
    // trigger
    //
    static TString LocalTriggerBoardMapping();
  
    // slats
    //
    static TString SlatFilePath(AliMpStationType stationType, 
                                const char* slatType, AliMpPlaneType plane);
    static TString SlatPCBFilePath(AliMpStationType stationType, 
                                const char* pcbType);
    // sectors
    //
    static TString SectorFilePath(AliMpStationType station, 
                                  AliMpPlaneType plane);
    static TString SectorSpecialFilePath(AliMpStationType station, 
                                  AliMpPlaneType plane);
    static TString SectorSpecialFilePath2(AliMpStationType station, 
                                  AliMpPlaneType plane);
    // motifs
    //
    static TString MotifFilePath(AliMpStationType station, 
                                 AliMpPlaneType plane, 
                                 const TString& motifTypeID);
    static TString MotifSpecialFilePath(AliMpStationType station,
                                 AliMpPlaneType plane, const TString& motifID);
    static TString PadPosFilePath(AliMpStationType station, 
                                 AliMpPlaneType plane, const TString& motifTypeID);
    static TString BergToGCFilePath(AliMpStationType station);
  
    // set methods
    static void SetTopPath(const TString& topPath);
  
  protected:
    AliMpFiles();
    AliMpFiles(const AliMpFiles& right);
  
    // operators
    AliMpFiles& operator=(const AliMpFiles& right);    
  
  private: 
    // methods
    static const char* GetDefaultTop();
    static TString PlaneDataDir(AliMpStationType station, AliMpPlaneType plane); 
    static TString StationDataDir(AliMpStationType station); 
  
    // static data members  
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
    static const TString fgkTriggerLocalBoards; // local board name to id mapping
  
    static TString  fgTop; // top directory path
    

  ClassDef(AliMpFiles, 0) //File names and paths 
};  

#endif //ALI_MP_FILES_H
