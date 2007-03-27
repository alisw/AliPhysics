#ifndef ALIITSINITGEOMETRY_H
#define ALIITSINITGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
$Id$
*/

/////////////////////////////////////////////////////////////////////
// Class to inilize AliITSgeom and the like for both simulation
//  and reconstruction.
/////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TString.h>

typedef enum {
  kvPPRasymmFMD=10,kv11=11,kv11Hybrid=110,kvDefault=1
} AliITSVersion_t;


class AliITSgeom;
class TArrayD;
class TGeoHMatrix;

class AliITSInitGeometry : public TObject{
 public:
    AliITSInitGeometry(AliITSVersion_t 
       version=kvPPRasymmFMD,Int_t minorversion=2);//Standard Constructor
    //virtual ~AliITSInitGeometry(); // Destructor
    //
    AliITSgeom* CreateAliITSgeom(); // Create and intilize geometry from TGeom
    Bool_t InitAliITSgeom(AliITSgeom *geom);//Initilize goemetry from gGeoManager
    // Getters and Setters
    TString GetGeometryName()const {return fName;}// Return geometry name
    void    SetGeometryName(const Char_t *name){fName = name;}// Set Geometry name
    Int_t   GetMajorVersion()const {return (Int_t)fMajorVersion;} // Return geometry major version
    Int_t   GetMinorVersion()const{return fMinorVersion;}// Return geometry minor version
    Bool_t  GetTiming()const{return fTiming;} // return routine timing flag
    void    SetTiming(Bool_t time=kTRUE){fTiming=time;}// Set routine timing (on)
    Bool_t  GetSegGeom()const{return fSegGeom;} // return flag indecating the use of AliITSsegmentation or AliITSgeomS?D class in fShape.
    void    SetSegGeom(Bool_t seg=kTRUE){fSegGeom = seg;}// Set the use of AliITSsegmentation class' instead of AliITSgeomS?D class in fShape
    Bool_t  GetDecoding()const{return fDecode;}// Return flag indecating wether to use new/old decoding
    void    SetDecoding(Bool_t newdec=kFALSE){fDecode = newdec;}// Set flag to use new/old decoding

    static const Bool_t SPDIsTGeoNative() {return !fgkOldSPDbarrel;}
    static const Bool_t SDDIsTGeoNative() {return !fgkOldSDDbarrel;}
    static const Bool_t SSDIsTGeoNative() {return !fgkOldSSDbarrel;}

    static const Bool_t SDDconeIsTGeoNative()   {return ! fgkOldSDDcone;} 
    static const Bool_t SSDconeIsTGeoNative()   {return ! fgkOldSSDcone;}
    static const Bool_t SPDshieldIsTGeoNative() {return ! fgkOldSPDshield;}
    static const Bool_t SDDshieldIsTGeoNative() {return ! fgkOldSDDshield; }
    static const Bool_t SSDshieldIsTGeoNative() {return ! fgkOldSSDshield;}
    static const Bool_t ServicesAreTGeoNative() {return ! fgkOldServices;}
    static const Bool_t SupportIsTGeoNative()   {return ! fgkOldSupports;}

 private:		   
    // Virtual MC code reproduction
    Bool_t InitAliITSgeomPPRasymmFMD(AliITSgeom *geom);
    Bool_t InitAliITSgeomV11Hybrid(AliITSgeom *geom);
    Bool_t InitAliITSgeomV11(AliITSgeom *geom);
    Bool_t InitGeomShapePPRasymmFMD(AliITSDetector idet,Bool_t *initSeg,
				       TArrayD &shapePar,AliITSgeom *geom);
    Bool_t InitSegmentationPPRasymmFMD(AliITSDetector idet,Bool_t *initSeg,
				       TArrayD &shapePar,AliITSgeom *geom);
    Bool_t GetTransformation(const TString &volumePath,TGeoHMatrix &mat);
    Bool_t GetShape(const TString &volumePath,TString &shapeType,TArrayD &par);
    void DecodeDetectorLayers(Int_t mod,Int_t &lay,Int_t &lad,Int_t &det);
    void DecodeDetector(Int_t &mod,Int_t lay,Int_t cpn0,
                        Int_t cpn1,Int_t cpn2) const;
    void RecodeDetector(Int_t mod,Int_t &cpn0,Int_t &cpn1,Int_t &cpn2);

    TString   fName;          // Geometry name
    Int_t           fMinorVersion;  // Geometry minor version
    AliITSVersion_t fMajorVersion;  // Geometry swich value
    Bool_t          fTiming;        // Flag to start inilization timing
    Bool_t          fSegGeom;       // Flag to switch between the old use of
                              // AliITSgeomS?D class, or AliITSsegmentation
                              // class in fShape of AliITSgeom class.
    Bool_t          fDecode;        // Flag for new/old decoding

    static const Bool_t fgkOldSPDbarrel;   // use old geo for SPD ?
    static const Bool_t fgkOldSDDbarrel;   // use old geo for SDD ?
    static const Bool_t fgkOldSSDbarrel;   // use old geo for SSD ?
    static const Bool_t fgkOldSDDcone;    // use old geo for SDD cone ?
    static const Bool_t fgkOldSSDcone;    // use old geo for SSD cone?
    static const Bool_t fgkOldSPDshield;  // use old geo for SPD shield ?
    static const Bool_t fgkOldSDDshield;  // use old geo for SDD shield ?
    static const Bool_t fgkOldSSDshield;  // use old geo for SDD shield ?
    static const Bool_t fgkOldServices;  // use old geo for services ?
    static const Bool_t fgkOldSupports;  // use old geo for supports ?

    ClassDef(AliITSInitGeometry,0) // create/Init AliITSgeom
    // 0 in ClassDef indicates that this class will not be "saved" in a file.
};

#endif

