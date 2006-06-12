#ifndef ALIITSINITGEOMETRY_H
#define ALIITSINITGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
$Id$
*/

/*
  Class to inilize AliITSgeom and the like for both simulation
  and reconstriction.
 */
#include <TObject.h>
#include <TArrayD.h>
#include <TString.h>
#include <TGeoMatrix.h>

#include "AliITSgeom.h"

class AliITSInitGeometry : public TObject{
 public:
    AliITSInitGeometry(); // Default Creator
    AliITSInitGeometry(const Char_t *name,Int_t minorversion);//Standard Creator
    //virtual ~AliITSInitGeometry(); // Destructor
    //
    AliITSgeom* CreateAliITSgeom(); // Create and intilize geometry from TGeom
    Bool_t InitAliITSgeom(AliITSgeom *geom);//Initilize goemetry from gGeoManager
    // Getters and Setters
    TString GetGeometryName()const {return fName;}// Return geometry name
    void    SetGeometryName(const Char_t *name){fName = name;}// Set Geometry name
    Int_t   GetMajorVersion()const {return fMajorVersion;} // Return geometry major version
    void    SetMajorVersion(Int_t majorVersion){fMajorVersion = majorVersion;} // Set geometry major version
    Int_t   GetMinorVersion()const{return fMinorVersion;}// Return geometry minor version
    void    SetMinorVersion(Int_t minorVersion){fMinorVersion = minorVersion;}
    Bool_t  GetTiming()const{return fTiming;} // return routine timing flag
    void    SetTiming(Bool_t time=kTRUE){fTiming=time;}// Set routine timing (on)
    Bool_t  GetSegGeom()const{return fSegGeom;} // return flag indecating the use of AliITSsegmentation or AliITSgeomS?D class in fShape.
    void    SetSegGeom(Bool_t seg=kTRUE){fSegGeom = seg;}// Set the use of AliITSsegmentation class' instead of AliITSgeomS?D class in fShape
    Bool_t  GetDecoding()const{return fDecode;}// Return flag indecating wether to use new/old decoding
    void    SetDecoding(Bool_t newdec=kFALSE){fDecode = newdec;}// Set flag to use new/old decoding

 private:
    // Virtual MC code reproduction
    Bool_t InitAliITSgeomPPRasymmFMD(AliITSgeom *geom);
    Bool_t InitGeomShapePPRasymmFMD(AliITSDetector idet,Bool_t *initSeg,
				       TArrayD &shapePar,AliITSgeom *geom);
    Bool_t InitSegmentationPPRasymmFMD(AliITSDetector idet,Bool_t *initSeg,
				       TArrayD &shapePar,AliITSgeom *geom);
    Bool_t GetTransformation(const TString &volumePath,TGeoHMatrix &mat);
    Bool_t GetShape(const TString &volumePath,TString &shapeType,TArrayD &par);
    void DecodeDetectorLayers(Int_t mod,Int_t &lay,Int_t &lad,Int_t &det);
    void DecodeDetector(Int_t &mod,Int_t lay,Int_t cpn0,Int_t cpn1,Int_t cpn2);
    void RecodeDetector(Int_t mod,Int_t &cpn0,Int_t &cpn1,Int_t &cpn2);

    TString   fName;          // Geometry name
    Int_t     fMinorVersion;  // Geometry minor version
    Int_t     fMajorVersion;  // Geometry swich value
    Bool_t    fTiming;        // Flag to start inilization timing
    Bool_t    fSegGeom;       // Flag to switch between the old use of
                              // AliITSgeomS?D class, or AliITSsegmentation
                              // class in fShape of AliITSgeom class.
    Bool_t    fDecode;        // Flag for new/old decoding

    ClassDef(AliITSInitGeometry,0) // create/Init AliITSgeom
    // 0 in ClassDef indicates that this class will not be "saved" in a file.
};

#endif

