#ifndef ALIITSINITGEOMETRYUPG_H
#define ALIITSINITGEOMETRYUPG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
$Id: AliITSInitGeometryUpg.h $
*/

/////////////////////////////////////////////////////////////////////
// Class to inilize AliITSgeom and the like for both simulation
//  and reconstruction.
/////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TString.h>
#include "AliITSgeom.h"

typedef enum {
  kvDefault=0,kv11=11,kvUpgrade=20
} AliITSVersion_t;

class TArrayD;
class TGeoHMatrix;
class TDatime;

class AliITSInitGeometryUpg : public TObject{
 public:

    AliITSInitGeometryUpg();//Default Constructor
    AliITSInitGeometryUpg(AliITSVersion_t version,
		       Int_t minorversion=2);//Standard Constructor
    //virtual ~AliITSInitGeometryUpg(); // Destructor
    //
    // Create and initialize geometry from TGeo
    AliITSgeom* CreateAliITSgeom();
    AliITSgeom* CreateAliITSgeom(Int_t major,Int_t minor); 
    Bool_t InitAliITSgeom(AliITSgeom *geom);//Initilize geometry from gGeoManager
    // Getters and Setters
    // Getters and Setters
    void    SetVersion(AliITSVersion_t maj,Int_t min) {// Set Major and Minor versions
        fMajorVersion=maj;fMinorVersion=min;}
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
     // Set debug level. debug=0 no debug info outputted.
    void    SetDebug(Int_t debug=0){fDebug=debug;};
    // Retrun debug value
    Int_t   GetDebug()const{return fDebug;};
    // Decode module number into old layer, ladder, and detector numbers
    void DecodeDetectorLayers(Int_t mod,Int_t &lay,Int_t &lad,Int_t &det);
    // find module number by layer, and copy numbers
    void DecodeDetector(Int_t &mod,Int_t lay,Int_t cpn0,
                        Int_t cpn1,Int_t cpn2) const;
    // Given module number, find copy numbers.
    void RecodeDetector(Int_t mod,Int_t &cpn0,Int_t &cpn1,Int_t &cpn2);
   // fills the string str with the major and minor version number
    Bool_t WriteVersionString(Char_t *str,Int_t length,
                              AliITSVersion_t maj,Int_t min,
                              const Char_t *cvsDate,const Char_t *cvsRev)const;
    // decodes the string str with the major and minor version number
    Bool_t ReadVersionString(const Char_t *str,Int_t length,
                             AliITSVersion_t &maj,Int_t &min,TDatime &dt)const;

    Int_t GetNumberOfLayers();
    Int_t GetNumberOfLadders(const Int_t lay)const;
    Int_t GetNumberOfModules(const Int_t lay)const;

 private:
    // Decode module number into old layer, ladder, and detector numbers
    void DecodeDetectorLayersv11(Int_t mod,Int_t &lay,
				 Int_t &lad,Int_t &det);
    // find module number by layer, and copy numbers
    void DecodeDetectorv11(Int_t &mod,Int_t lay,Int_t cpn0,Int_t cpn1,
			   Int_t cpn2)const;
    // Given module number, find copy numbers.
    void RecodeDetectorv11(Int_t mod,Int_t &cpn0,Int_t &cpn1,
			   Int_t &cpn2);
    // Decode module number into old layer, ladder, and detector numbers
    void DecodeDetectorvUpgrade(Int_t &mod,Int_t lay,Int_t cpn0,Int_t cpn1,
                                    Int_t cpn2)const;
    // find module number by layer, and copy numbers
    void RecodeDetectorvUpgrade(Int_t mod,Int_t &cpn0,Int_t &cpn1,
                                    Int_t &cpn2);
    // Given module number, find copy numbers.
    void DecodeDetectorLayersvUpgrade(Int_t mod,Int_t &lay,Int_t &lad,
                                    Int_t &det);
		   
    // Virtual MC code 
    Bool_t InitAliITSgeomV11(AliITSgeom *geom);
    Bool_t InitAliITSgeomVUpgrade(AliITSgeom *geom);
    Bool_t GetTransformation(const TString &volumePath,TGeoHMatrix &mat);
    Bool_t GetShape(const TString &volumePath,TString &shapeType,TArrayD &par);
    void TransposeTGeoHMatrix(TGeoHMatrix *m) const;

    TString         fName;         // Geometry name
    Int_t           fMinorVersion; // Geometry minor version
    AliITSVersion_t fMajorVersion; // Geometry swich value
    Bool_t          fTiming;       // Flag to start inilization timing
    Bool_t          fSegGeom;      // Flag to switch between the old use of
                                   // AliITSgeomS?D class, or AliITSsegmentation
                                   // class in fShape of AliITSgeom class.
    Bool_t          fDecode;       // Flag for new/old decoding
    Int_t           fDebug;        // Debug flag

    ClassDef(AliITSInitGeometryUpg,0) // create/Init AliITSgeom
    // 0 in ClassDef indicates that this class will not be "saved" in a file.
};

#endif

