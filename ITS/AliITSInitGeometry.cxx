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

/*
$Id$
*/
////////////////////////////////////////////////////////////////
//  This class initializes the class AliITSgeom
//  The initialization is done starting from 
//  a geometry coded by means of the ROOT geometrical modeler
//  This initialization can be used both for simulation and reconstruction
///////////////////////////////////////////////////////////////

#include <TArrayD.h>
#include <TArrayF.h>
#include <TStopwatch.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoShape.h>
#include <TGeoBBox.h>
#include <TGeoTrd1.h>
#include <TGeoTrd2.h>
#include <TGeoArb8.h>
#include <TGeoTube.h>
#include <TGeoCone.h>
#include <TGeoSphere.h>
#include <TGeoPara.h>
#include <TGeoPgon.h>
#include <TGeoPcon.h>
#include <TGeoEltu.h>
#include <TGeoHype.h>
#include <TMath.h>

#include "AliLog.h"
#include "AliITSgeomSPD.h"
#include "AliITSgeomSDD.h"
#include "AliITSgeomSSD.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSsegmentationSSD.h"
#include "AliITSgeom.h"
#include "AliITSInitGeometry.h"
#include <TDatime.h>

ClassImp(AliITSInitGeometry)

const Bool_t AliITSInitGeometry::fgkOldSPDbarrel = kFALSE;
const Bool_t AliITSInitGeometry::fgkOldSDDbarrel = kFALSE;
const Bool_t AliITSInitGeometry::fgkOldSSDbarrel = kFALSE;
const Bool_t AliITSInitGeometry::fgkOldSDDcone   = kTRUE;
const Bool_t AliITSInitGeometry::fgkOldSSDcone   = kTRUE;
const Bool_t AliITSInitGeometry::fgkOldSPDshield = kTRUE;
const Bool_t AliITSInitGeometry::fgkOldSDDshield = kTRUE;
const Bool_t AliITSInitGeometry::fgkOldSSDshield = kTRUE;
const Bool_t AliITSInitGeometry::fgkOldServices  = kTRUE;
const Bool_t AliITSInitGeometry::fgkOldSupports  = kTRUE;
//______________________________________________________________________
AliITSInitGeometry::AliITSInitGeometry():
TObject(),                   // Base Class
fName(0),                    // Geometry name
fMinorVersion(-1),           // Minor version number/type
fMajorVersion(kvDefault),    // Major versin number
fTiming(kFALSE),             // Flag to start inilization timing
fSegGeom(kFALSE),            // Flag to switch between the old use of
                             // AliITSgeomS?D class, or AliITSsegmentation
                             // class in fShape of AliITSgeom class.
fDecode(kFALSE),             // Flag for new/old decoding
fDebug(0){                   // Debug flag
    // Default Creator
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   A default inilized AliITSInitGeometry object

    fName = "Undefined";
}
//______________________________________________________________________
AliITSInitGeometry::AliITSInitGeometry(AliITSVersion_t version,
                                       Int_t minorversion):
TObject(),                   // Base Class
fName(0),                    // Geometry name
fMinorVersion(minorversion), // Minor version number/type
fMajorVersion(version),      // Major versin number
fTiming(kFALSE),             // Flag to start inilization timing
fSegGeom(kFALSE),            // Flag to switch between the old use of
                             // AliITSgeomS?D class, or AliITSsegmentation
                             // class in fShape of AliITSgeom class.
fDecode(kFALSE),             // Flag for new/old decoding
fDebug(0){                   // Debug flag
    // Default Creator
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   A default inilized AliITSInitGeometry object

    if(version == kvPPRasymmFMD && (fMinorVersion==1|| fMinorVersion==2)){
        fName="AliITSvPPRasymmFMD";
    }else if(version == kv11Hybrid){
        fName="AliITSv11Hybrid";
    }else {
        AliFatal(Form("Undefined geometry: fMajorVersion=%d, "
                      "fMinorVersion= %d",(Int_t)fMajorVersion,fMinorVersion));
        fName = "Undefined";
    } // end if
    return;
}
//______________________________________________________________________
AliITSgeom* AliITSInitGeometry::CreateAliITSgeom(){
    // Creates and Initilizes the geometry transformation class AliITSgeom
    // to values appropreate to this specific geometry. Now that
    // the segmentation is part of AliITSgeom, the detector
    // segmentations are also defined here.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   A pointer to a new properly inilized AliITSgeom class. If
    //   pointer = 0 then failed to init.


  AliITSVersion_t version = kvDefault;
  Int_t minor = 0;
  TDatime datetime;
  TGeoVolume *itsV = gGeoManager->GetVolume("ITSV");
  if(!itsV){
    Error("CreateAliITSgeom","Can't find ITS volume ITSV, aborting");
    return 0;
  }// end if
  const Char_t *title = itsV->GetTitle();
  if(!ReadVersionString(title,(Int_t)strlen(title),version,minor,
			datetime))
    Warning("UpdateInternalGeometry","Can't read title=%s\n",title);
  SetTiming(kFALSE);
  SetSegGeom(kFALSE);
  SetDecoding(kFALSE);
  AliITSgeom *geom = CreateAliITSgeom(version,minor);
  AliDebug(1,"AliITSgeom object has been initialized from TGeo\n");
  return geom;
}
//______________________________________________________________________
AliITSgeom* AliITSInitGeometry::CreateAliITSgeom(Int_t major,Int_t minor){
    // Creates and Initilizes the geometry transformation class AliITSgeom
    // to values appropreate to this specific geometry. Now that
    // the segmentation is part of AliITSgeom, the detector
    // segmentations are also defined here.
    // Inputs:
    //   Int_t major   major version, see AliITSVersion_t
    //   Int_t minor   minor version
    // Outputs:
    //   none.
    // Return:
    //   A pointer to a new properly inilized AliITSgeom class. If
    //   pointer = 0 then failed to init.

    switch(major){
    case kvtest:
        SetGeometryName("AliITSvtest");
        SetVersion(kvtest,minor);
        break;
    case kvSPD02:
        SetGeometryName("AliITSvSPD02");
        SetVersion(kvSPD02,minor);
        break;
    case kvSDD03:
        SetGeometryName("AliITSvSDD03");
        SetVersion(kvSDD03,minor);
        break;
    case kvSSD03:
        SetGeometryName("AliITSvSSD03");
        SetVersion(kvSSD03,minor);
        break;
    case kvITS04:
        SetGeometryName("AliITSvBeamTest03");
        SetVersion(kvITS04,minor);
        break;
    case kvPPRcourseasymm:
        SetGeometryName("AliITSvPPRcourseasymm");
        SetVersion(kvPPRcourseasymm,minor);
        break;
    case kvPPRasymmFMD:
        SetGeometryName("AliITSvPPRasymmFMD");
        SetVersion(kvPPRasymmFMD,minor);
        break;
    case kv11:
        SetGeometryName("AliITSv11");
        SetVersion(kv11,minor);
        break;
    case kv11Hybrid:
        SetGeometryName("AliITSv11Hybrid");
        SetVersion(kv11Hybrid,minor);
        break;
    case kvDefault:
    default:
        SetGeometryName("Undefined");
        SetVersion(kvDefault,minor);
        break;
    } // end switch
    AliITSgeom *geom = new AliITSgeom();
    if(!InitAliITSgeom(geom)){ // Error initilization failed
	delete geom;
	geom = 0;
    } // end if
    return geom;
}
//______________________________________________________________________
Bool_t AliITSInitGeometry::InitAliITSgeom(AliITSgeom *geom){
  // Initilizes the geometry transformation class AliITSgeom
  // to values appropreate to this specific geometry. Now that
  // the segmentation is part of AliITSgeom, the detector
  // segmentations are also defined here.
  // Inputs:
  //   AliITSgeom *geom  A pointer to the AliITSgeom class
  // Outputs:
  //   AliITSgeom *geom  This pointer recreated and properly inilized.
  // Return:
  //   none.

    if(!gGeoManager){
        AliFatal("The geometry manager has not been initialized (e.g. "
                 "TGeoManager::Import(\"geometry.root\")should be "
                 "called in advance) - exit forced");
        return kFALSE;
    } // end if
    switch(fMajorVersion) {
    case kvtest: {
        if(GetMinorVersion()==1) return InitAliITSgeomPPRasymmFMD(geom);
        else if(GetMinorVersion()==2) return InitAliITSgeomtest2(geom);
    } break; // end case
    case kvSPD02: { 
        return InitAliITSgeomSPD02(geom);
    } break; // end case
    case kvSDD03: { 
        return InitAliITSgeomSDD03(geom);
    } break; // end case
    case kvSSD03: { 
        return InitAliITSgeomSSD03(geom);
    } break; // end case
    case kvITS04: { 
        return InitAliITSgeomITS04(geom);
    } break; // end case
    case kvPPRasymmFMD: { 
        return InitAliITSgeomPPRasymmFMD(geom);
    } break; // end case
    case kvPPRcourseasymm: { 
        return kTRUE; // No sensitive detectors in course geometry
    } break; // end case
    case kv11Hybrid: { 
        return InitAliITSgeomV11Hybrid(geom);
    } break; // end case
    case kv11: {
        return InitAliITSgeomV11(geom);
    } break; // end case
    case kvDefault: default: {
        AliFatal("Undefined geometry");
        return kFALSE;
    } break; // end case
    } // end switch
    return kFALSE;
}
//______________________________________________________________________
void AliITSInitGeometry::TransposeTGeoHMatrix(TGeoHMatrix *m)const{
    // Transpose the rotation matrix part of a TGeoHMatrix. This
    // is needed because TGeo stores the transpose of the rotation
    // matrix as compared to what AliITSgeomMatrix uses (and Geant3).
    // Inputs:
    //    TGeoHMatrix *m  The matrix to be transposed
    // Outputs:
    //    TGEoHMatrix *m  The transposed matrix
    // Return:
    //    none.
    Int_t i;
    Double_t r[9];

    if(m==0) return; // no matrix to transpose.
    for(i=0;i<9;i += 4) r[i] = m->GetRotationMatrix()[i]; // diagonals
    r[1] = m->GetRotationMatrix()[3];
    r[2] = m->GetRotationMatrix()[6];
    r[3] = m->GetRotationMatrix()[1];
    r[5] = m->GetRotationMatrix()[7];
    r[6] = m->GetRotationMatrix()[2];
    r[7] = m->GetRotationMatrix()[5];
    m->SetRotation(r);
    return;
}
//______________________________________________________________________
Bool_t AliITSInitGeometry::InitAliITSgeomtest2(AliITSgeom *geom){
    // Initilizes the geometry transformation class AliITSgeom
    // to values appropreate to this specific geometry. Now that
    // the segmentation is part of AliITSgeom, the detector
    // segmentations are also defined here.
    // Inputs:
    //   AliITSgeom *geom  A pointer to the AliITSgeom class
    // Outputs:
    //   AliITSgeom *geom  This pointer recreated and properly inilized.
    // Return:
    //   none.
    //  const Double_t kcm2micron = 1.0E4;
    const Int_t kItype=0; // Type of transormation defined 0=> Geant
    const Int_t klayers = 6; // number of layers in the ITS
    const Int_t kladders[klayers]   = {1,1,1,1,1,1}; // Number of ladders
    const Int_t kdetectors[klayers] = {1,1,1,1,1,1};// number of detector/lad
    const AliITSDetector kIdet[6]   = {kSPD,kSPD,kSDD,kSDD,kSSD,kSSD};
    const TString kNames[klayers] = {
	 "/ALIC_1/ITSV_1/ITSspd1_1/ITS1_1", // lay=1
	 "/ALIC_1/ITSV_1/ITSspd2_1/ITS2_1", // lay=2
	 "/ALIC_1/ITSV_1/ITSsdd1_1/ITS3_1", // lay=3
	 "/ALIC_1/ITSV_1/ITSsdd2_1/ITS4_1", // lay=4
	 "/ALIC_1/ITSV_1/ITSssd1_1/ITS5_1", // lay=5
	 "/ALIC_1/ITSV_1/ITSssd2_1/ITS6_1"};// Lay=6
    Int_t mod,nmods=0,lay,lad,det,cpn0,cpn1,cpn2;
    Double_t tran[3]={0.0,0.0,0.0},rot[10]={9*0.0,1.0};
    TArrayD shapePar;
    TString shapeName;
    TGeoHMatrix matrix;
    Bool_t initSeg[3]={kFALSE,kFALSE,kFALSE};
    TStopwatch *time = 0x0;if(fTiming) time=new TStopwatch();

    if(fTiming) time->Start();
    for(mod=0;mod<klayers;mod++) nmods += kladders[mod]*kdetectors[mod];
    geom->Init(kItype,klayers,kladders,kdetectors,nmods);
    for(mod=0;mod<nmods;mod++){
        DecodeDetectorLayers(mod,lay,lad,det); // Write
        geom->CreateMatrix(mod,lay,lad,det,kIdet[lay-1],tran,rot);
        RecodeDetector(mod,cpn0,cpn1,cpn2); // Write reusing lay,lad,det.
        geom->GetGeomMatrix(mod)->SetPath(kNames[lay-1]);
        GetTransformation(kNames[lay-1].Data(),matrix);
        geom->SetTrans(mod,matrix.GetTranslation());
        TransposeTGeoHMatrix(&matrix); // Transpose TGeo's rotation matrixes
        geom->SetRotMatrix(mod,matrix.GetRotationMatrix());
        if(initSeg[kIdet[lay-1]]) continue;
        GetShape(kNames[lay-1],shapeName,shapePar);
        if(shapeName.CompareTo("BOX")){
            Error("InitITSgeom2","Geometry changed without proper code update"
                  "or error in reading geometry. Shape is not BOX shape is %s",
                  shapeName.Data());
            return kFALSE;
        } // end if
	InitGeomShapePPRasymmFMD(kIdet[lay-1],initSeg,shapePar,geom);
    } // end for module
    if(fTiming){
        time->Stop();
        time->Print();
        delete time;
    } // end if
    return kTRUE;
}
//______________________________________________________________________
Bool_t AliITSInitGeometry::InitAliITSgeomSPD02(AliITSgeom *geom){
    // Initilizes the geometry transformation class AliITSgeom
    // to values appropreate to this specific geometry. Now that
    // the segmentation is part of AliITSgeom, the detector
    // segmentations are also defined here.
    // Inputs:
    //   AliITSgeom *geom  A pointer to the AliITSgeom class
    // Outputs:
    //   AliITSgeom *geom  This pointer recreated and properly inilized.
    // Return:
    //   none.
    const Int_t kltypess=2;
    const Int_t knlayers=5;
    const TString knames[kltypess]=
        {"ALIC_1/ITSV_1/ITEL_%d/IMB0_1/IMBS_1",//lay=1,2,4,5
         "ALIC_1/ITSV_1/IDET_%d/ITS0_1/ITST_1"};// lay=3
    const Int_t kitsGeomTreeCopys[2]={4,1};
    const Int_t knlad[knlayers]={knlayers*1},kndet[knlayers]={knlayers*1};
    TString path,shapeName;
    TGeoHMatrix matrix;
    TArrayD shapePar;
    TArrayF shapeParF;
    Double_t trans[3]={3*0.0},rot[10]={10*0.0};
    Int_t npar=3,mod,i,j,lay,lad,det,cpy;
    Float_t par[20];
    TStopwatch *time = 0x0;if(fTiming) time=new TStopwatch();

    par[0]=0.64;par[1]=0.5*300.0E-4;par[2]=3.48;
    mod=5;;
    geom->Init(0,knlayers,knlad,kndet,mod);

    if(fTiming) time->Start();
    for(i=0;i<kltypess;i++)for(cpy=1;cpy<=kitsGeomTreeCopys[i];cpy++){
        path.Form(knames[i].Data(),cpy);
        GetTransformation(path.Data(),matrix);
        GetShape(path.Data(),shapeName,shapePar);
        shapeParF.Set(shapePar.GetSize());
        for(j=0;j<shapePar.GetSize();j++) shapeParF[j]=shapePar[j];
        lay = cpy;
        if(i==0&&cpy>2) lay=cpy+1;
        if(i==1) lay=3;
        DecodeDetector(mod,kitsGeomTreeCopys[i],1,cpy,0);
        DecodeDetectorLayers(mod,lay,lad,det);
        geom->CreateMatrix(mod,lay,lad,det,kSPD,trans,rot);
        geom->SetTrans(mod,matrix.GetTranslation());
        geom->SetRotMatrix(mod,matrix.GetRotationMatrix());
        geom->GetGeomMatrix(mod)->SetPath(path.Data());
        if(!(geom->IsShapeDefined((Int_t)kSPD)))
            geom->ReSetShape(kSPD,new AliITSgeomSPD425Short(npar,par));
    } // end for i,cpy/
    if(fTiming){
        time->Stop();
        time->Print();
        delete time;
    } // end if
    return kTRUE;
}
//______________________________________________________________________
Bool_t AliITSInitGeometry::InitAliITSgeomSDD03(AliITSgeom *geom){
    // Initilizes the geometry transformation class AliITSgeom
    // to values appropreate to this specific geometry. Now that
    // the segmentation is part of AliITSgeom, the detector
    // segmentations are also defined here.
    // Inputs:
    //   AliITSgeom *geom  A pointer to the AliITSgeom class
    // Outputs:
    //   AliITSgeom *geom  This pointer recreated and properly inilized.
    // Return:
    //   none
    const Int_t knlayers=12;
    //   const Int_t kndeep=6;
    const Int_t kltypess=2;
    const AliITSDetector kidet[knlayers]={kSSD,kSDD};
    const TString knames[kltypess]={
        "/ALIC_1/ITSV_1/ITEL_%d/ITAI_1/IMB0_1/IMBS_1",
        "/ALIC_1/ITSV_1/IDET_%d/IDAI_1/ITS0_1/ITST_1"};
    const Int_t kitsGeomTreeCopys[kltypess]={10,2};
    const Int_t knp=384;
    const Float_t kpitch=50.E-4;/*cm*/
    Float_t box[3]={0.5*kpitch*(Float_t)knp,150.E-4,1.0},p[knp+1],n[knp+1];
    Int_t nlad[knlayers]={knlayers*1};
    Int_t ndet[knlayers]={knlayers*1};
    Int_t mod=knlayers,lay=0,lad=0,det=0,i,j,cp0;
    TString path,shapeName;
    TGeoHMatrix matrix;
    Double_t trans[3]={3*0.0},rot[10]={10*0.0};
    TArrayD shapePar;
    TArrayF shapeParF;
    Bool_t isShapeDefined[kltypess]={kltypess*kFALSE};

    geom->Init(0,knlayers,nlad,ndet,mod);
    p[0]=-box[0];
    n[0]=box[0];
    // Fill in anode and cathode strip locations (lower edge)
    for(i=1;i<knp;i++){
        p[i] =p[i-1]+kpitch;
        n[i] =n[i-1]-kpitch;
    } // end for i
    p[knp]=box[0];
    n[knp]=-box[0];
    for(i=0;i<kltypess;i++)for(cp0=1;cp0<=kitsGeomTreeCopys[i];cp0++){
        DecodeDetector(mod,kitsGeomTreeCopys[i],cp0,1,2);
        DecodeDetectorLayers(mod,lay,lad,det);
        path.Form(knames[i].Data(),cp0);
        GetTransformation(path.Data(),matrix);
        GetShape(path.Data(),shapeName,shapePar);
        shapeParF.Set(shapePar.GetSize());
        for(j=0;j<shapePar.GetSize();j++)shapeParF[j]=shapePar[j];
        geom->CreateMatrix(mod,lay,lad,det,kidet[i],trans,rot);
        geom->SetTrans(mod,matrix.GetTranslation());
        geom->SetRotMatrix(mod,matrix.GetRotationMatrix());
        geom->GetGeomMatrix(mod)->SetPath(path.Data());
        switch (kidet[i]){
        case kSDD: if(!(geom->IsShapeDefined((Int_t)kSDD))){
            geom->ReSetShape(kSDD,new AliITSgeomSDD256(shapeParF.GetSize(),
                                                       shapeParF.GetArray()));
            isShapeDefined[i]=kTRUE;
        } break;
        case kSSD:if(!(geom->IsShapeDefined((Int_t)kSSD))){
            geom->ReSetShape(kSSD,new AliITSgeomSSD(box,0.0,0.0,
                                                    knp+1,p,knp+1,n));
            isShapeDefined[i]=kTRUE;
        } break;
        default:{} break;
        } // end switch
    } // end for i,cp0

    return kTRUE;
}
//______________________________________________________________________
Bool_t AliITSInitGeometry::InitAliITSgeomSSD03(AliITSgeom *geom){
    // Initilizes the geometry transformation class AliITSgeom
    // to values appropreate to this specific geometry. Now that
    // the segmentation is part of AliITSgeom, the detector
    // segmentations are also defined here.
    // Inputs:
    //   AliITSgeom *geom  A pointer to the AliITSgeom class
    // Outputs:
    //   AliITSgeom *geom  This pointer recreated and properly inilized.
    // Return:
    //   none.
    const Int_t knlayers=5;
    //   const Int_t kndeep=6;
    const Int_t kltypess=3;
    const AliITSDetector kIdet[knlayers]={kND,kSSD,kND};
    const TString knames[kltypess]={
        "/ALIC_1/ITSV_1/ITSA_%d/ITSS_1",
        "/ALIC_1/ITSV_1/IGAR_%d/IAIR_1/ITST_1",
        "/ALIC_1/ITSV_1/IFRA_%d/IFRS_1"};
    const Int_t kitsGeomTreeCopys[kltypess]={3,1,1};
    const Int_t kitsGeomDetTypes[kltypess]={1,2,3};
    const Int_t knp=384;
    const Float_t kpitch=50.E-4;//cm
    Bool_t initSeg[3]={kFALSE, kFALSE, kFALSE};
    Float_t box[3]={0.5*kpitch*(Float_t)knp,150.E-4,1.0},p[knp+1],n[knp+1];
    Int_t nlad[knlayers]={knlayers*1};
    Int_t ndet[knlayers]={knlayers*1};
    Int_t mod=knlayers,lay=0,lad=0,det=0,i,j,cp0;
    TString path,shapeName;
    TGeoHMatrix matrix;
    Double_t trans[3]={3*0.0},rot[10]={10*0.0};
    TArrayD shapePar;
    TArrayF shapeParF;
    Bool_t isShapeDefined[kltypess]={kltypess*kFALSE};

    geom->Init(0,knlayers,nlad,ndet,mod);
    p[0]=-box[0];
    n[0]=box[0];
    // Fill in anode and cathode strip locations (lower edge)
    for(i=1;i<knp;i++){
        p[i] =p[i-1]+kpitch;
        n[i] =n[i-1]-kpitch;
    } // end for i
    p[knp]=box[0];
    n[knp]=-box[0];
    for(i=0;i<kltypess;i++)for(cp0=1;cp0<=kitsGeomTreeCopys[i];cp0++){
        DecodeDetector(mod,kitsGeomDetTypes[i],cp0,1,0);
        DecodeDetectorLayers(mod,lay,lad,det);
        path.Form(knames[i].Data(),cp0);
        GetTransformation(path.Data(),matrix);
        GetShape(path.Data(),shapeName,shapePar);
        shapeParF.Set(shapePar.GetSize());
        for(j=0;j<shapePar.GetSize();j++)shapeParF[j]=shapePar[j];
        geom->CreateMatrix(mod,lay,lad,det,kIdet[i],trans,rot);
        geom->SetTrans(mod,matrix.GetTranslation());
        geom->SetRotMatrix(mod,matrix.GetRotationMatrix());
        geom->GetGeomMatrix(mod)->SetPath(path.Data());
        switch (kIdet[i]){
        case kSSD:if(!(geom->IsShapeDefined((Int_t)kSSD))){
            InitGeomShapePPRasymmFMD(kIdet[lay-1],initSeg,shapePar,geom);
            isShapeDefined[i]=kTRUE;
        } break;
        default:{} break;
        } // end switch
    } // end for i,cp0

    return kTRUE;
}
//______________________________________________________________________
Bool_t AliITSInitGeometry::InitAliITSgeomITS04(AliITSgeom *geom){
    // Initilizes the geometry transformation class AliITSgeom
    // to values appropreate to this specific geometry. Now that
    // the segmentation is part of AliITSgeom, the detector
    // segmentations are also defined here.
    // Inputs:
    //   AliITSgeom *geom  A pointer to the AliITSgeom class
    // Outputs:
    //   AliITSgeom *geom  This pointer recreated and properly inilized.
    // Return:
    //   none.

  // We can not use AliITSvBeamTestITS04::fgk... data members because
  // AliITSInitGeometry is part of the base library while AliITSvBeamTestITS04
  // is part of the simulation library. This would introduce a dependance
  // between the 2 libraries


    const Int_t knlayers = 6;
    Int_t nlad[knlayers], ndet[knlayers];
    
    nlad[0] = 1; ndet[0] = 2;
    nlad[1] = 1; ndet[1] = 2;
    nlad[2] = 1; ndet[2] = 1;
    nlad[3] = 1; ndet[3] = 1;
    nlad[4] = 1; ndet[4] = 2;
    nlad[5] = 1; ndet[5] = 2;

    Int_t nModTot = 10;
    geom->Init(0,knlayers,nlad,ndet,nModTot);

    /*
    //=== Set default shapes 
    const Float_t kDxyzSPD[] = {AliITSvBeamTestITS04::fgkSPDwidthSens/2,
                                AliITSvBeamTestITS04::fgkSPDthickSens/2,
                                AliITSvBeamTestITS04::fgkSPDlengthSens/2};  
    if(!(geom->IsShapeDefined(kSPD)))
       geom->ReSetShape(kSPD,new AliITSgeomSPD425Short(3,(Float_t *)kDxyzSPD));
    
    const Float_t kDxyzSDD[] = {AliITSvBeamTestITS04::fgkSDDwidthSens/2.,
                                AliITSvBeamTestITS04::fgkSDDthickSens/2.,
                                AliITSvBeamTestITS04::fgkSDDlengthSens/2.};
    if(!(geom->IsShapeDefined(kSDD)))
	geom->ReSetShape(kSDD, new AliITSgeomSDD256(3,(Float_t *)kDxyzSDD));
    
    const Float_t kDxyzSSD[] = {AliITSvBeamTestITS04::fgkSSDlengthSens/2,
                                AliITSvBeamTestITS04::fgkSSDthickSens/2,
                                AliITSvBeamTestITS04::fgkSSDwidthSens/2};
    if(!(geom->IsShapeDefined(kSSD)))
       geom->ReSetShape(kSSD,new AliITSgeomSSD75and275(3,(Float_t *)kDxyzSSD));
    
    // Creating the matrices in AliITSgeom for each sensitive volume
    // (like in AliITSv11GeometrySDD) mln
    // Here, each layer is one detector
    
    char layerName[30];
    Int_t startMod = 0,mod;
    TGeoVolume *itsmotherVolume = gGeoManager->GetVolume("ITSV");
    // SPD
    for (Int_t i=0; i<4;i++) {
	sprintf(layerName, "ITSspdWafer_%i",i+1);
	TGeoNode *layNode = itsmotherVolume->GetNode(layerName);
	if (layNode) {
	    TGeoHMatrix layMatrix(*layNode->GetMatrix());	    
	    Double_t *trans  = layMatrix.GetTranslation();
	    Double_t *r      = layMatrix.GetRotationMatrix();
	    Double_t rot[10] = {r[0],r[1],r[2],
				r[3],r[4],r[5],
				r[6],r[7],r[8], 1.0};
	    Int_t iDet = 1;
	    Int_t iLad = 1;
	    Int_t iLay = 1;
            DecodeDetector(mod,layNode->GetNumber(),i+1,0,0);
            DecodeDetectorLayers(mod,iLay,iLad,iDet);
	    geom->CreateMatrix(startMod,iLay,iLad,iDet,kSPD,trans,rot);
	    startMod++;
	};
    };
    
    // SDD
    for (Int_t i=0; i<2;i++) {
	sprintf(layerName, "ITSsddWafer_%i",i+4+1);
	TGeoNode *layNode = itsmotherVolume->GetNode(layerName);
	if (layNode) {
	    TGeoHMatrix layMatrix(*layNode->GetMatrix());
	    Double_t *trans  = layMatrix.GetTranslation();
	    Double_t *r      = layMatrix.GetRotationMatrix();
	    Double_t rot[10] = {r[0],r[1],r[2],
				r[3],r[4],r[5],
				r[6],r[7],r[8], 1.0};
	    Int_t iDet = 1;
	    Int_t iLad = 1;
	    Int_t iLay = 1;
            DecodeDetector(mod,layNode->GetNumber(),i+1,0,0);
            DecodeDetectorLayers(mod,iLay,iLad,iDet);
	    geom->CreateMatrix(startMod,iLay,iLad,iDet,kSDD,trans,rot);
	    startMod++;
	};
    };
    
    // SSD
    for (Int_t i=0; i<4;i++) {
	sprintf(layerName, "ITSssdWafer_%i",i+4+2+1);
	TGeoNode *layNode = itsmotherVolume->GetNode(layerName);
	if (layNode) {
	    TGeoHMatrix layMatrix(*layNode->GetMatrix());	    
	    Double_t *trans  = layMatrix.GetTranslation();
	    Double_t *r      = layMatrix.GetRotationMatrix();
	    Double_t rot[10] = {r[0],r[1],r[2],
				r[3],r[4],r[5],
				r[6],r[7],r[8], 1.0};
	    Int_t iDet = 1;
	    Int_t iLad = 1;
	    Int_t iLay = 5;
            DecodeDetector(mod,layNode->GetNumber(),i+1,0,0);
            DecodeDetectorLayers(mod,iLay,iLad,iDet);
	    geom->CreateMatrix(startMod,iLay,iLad,iDet,kSSD,trans,rot);
	    startMod++;
	};
    };

    return kTRUE;
  */
    return kFALSE;
}
//______________________________________________________________________
Bool_t AliITSInitGeometry::InitAliITSgeomPPRasymmFMD(AliITSgeom *geom){
    // Initilizes the geometry transformation class AliITSgeom
    // to values appropreate to this specific geometry. Now that
    // the segmentation is part of AliITSgeom, the detector
    // segmentations are also defined here.
    // Inputs:
    //   AliITSgeom *geom  A pointer to the AliITSgeom class
    // Outputs:
    //   AliITSgeom *geom  This pointer recreated and properly inilized.
    // Return:
    //   none.
  //    const Double_t kcm2micron = 1.0E4;
    const Int_t kItype=0; // Type of transormation defined 0=> Geant
    const Int_t klayers = 6; // number of layers in the ITS
    const Int_t kladders[klayers]   = {20,40,14,22,34,38}; // Number of ladders
    const Int_t kdetectors[klayers] = {4,4,6,8,22,25};// number of detector/lad
    const AliITSDetector kIdet[6]   = {kSPD,kSPD,kSDD,kSDD,kSSD,kSSD};
    const TString kPathbase = "/ALIC_1/ITSV_1/ITSD_1/";
    const TString kNames[2][klayers] = {
	{"%sIT12_1/I12A_%d/I10A_%d/I103_%d/I101_1/ITS1_1", // lay=1
	 "%sIT12_1/I12A_%d/I20A_%d/I1D3_%d/I1D1_1/ITS2_1", // lay=2
	 "%sIT34_1/I004_%d/I302_%d/ITS3_%d/", // lay=3
	 "%sIT34_1/I005_%d/I402_%d/ITS4_%d/", // lay=4
	 "%sIT56_1/I565_%d/I562_%d/ITS5_%d/", // lay=5
	 "%sIT56_1/I569_%d/I566_%d/ITS6_%d/"},// lay=6
// 	{"%sIT12_1/I12B_%d/I10B_%d/I107_%d/I101_1/ITS1_1", // lay=1
// 	 "%sIT12_1/I12B_%d/I20B_%d/I1D7_%d/I1D1_1/ITS2_1", // lay=2
	{"%sIT12_1/I12B_%d/I10B_%d/L1H-STAVE%d_1/I107_%d/I101_1/ITS1_1", // lay=1
	 "%sIT12_1/I12B_%d/I20B_%d/L2H-STAVE%d_1/I1D7_%d/I1D1_1/ITS2_1", // lay=2
	 "%sIT34_1/I004_%d/I302_%d/ITS3_%d", // lay=3
	 "%sIT34_1/I005_%d/I402_%d/ITS4_%d", // lay=4
	 "%sIT56_1/I565_%d/I562_%d/ITS5_%d", // lay=5
	 "%sIT56_1/I569_%d/I566_%d/ITS6_%d"}};// Lay=6
    /*
      Int_t itsGeomTreeCopys[knlayers][3]= {{10, 2, 4},// lay=1
      {10, 4, 4},// lay=2
      {14, 6, 1},// lay=3
      {22, 8, 1},// lay=4
      {34,22, 1},// lay=5
      {38,25, 1}};//lay=6
    */
    Int_t mod,nmods=0,lay,lad,det,cpn0,cpn1,cpn2, cpnHS;
    Double_t tran[3]={0.0,0.0,0.0},rot[10]={9*0.0,1.0};
    TArrayD shapePar;
    TString path,shapeName;
    TGeoHMatrix matrix;
    Bool_t initSeg[3]={kFALSE,kFALSE,kFALSE};
    TStopwatch *time = 0x0;if(fTiming) time=new TStopwatch();

    if(fTiming) time->Start();
    for(mod=0;mod<klayers;mod++) nmods += kladders[mod]*kdetectors[mod];
    geom->Init(kItype,klayers,kladders,kdetectors,nmods);
    for(mod=0;mod<nmods;mod++){
        DecodeDetectorLayers(mod,lay,lad,det); // Write
        geom->CreateMatrix(mod,lay,lad,det,kIdet[lay-1],tran,rot);
        RecodeDetector(mod,cpn0,cpn1,cpn2); // Write reusing lay,lad,det.

	if (kIdet[lay-1]==kSPD) { // we need 1 more copy number because of the half-stave
	  if (det<3) cpnHS = 0; else cpnHS = 1;
	  path.Form(kNames[fMinorVersion-1][lay-1].Data(),kPathbase.Data(),
		    cpn0,cpn1,cpnHS,cpn2);
	} else {
	  path.Form(kNames[fMinorVersion-1][lay-1].Data(),kPathbase.Data(),
		    cpn0,cpn1,cpn2);
	};
//         path.Form(kNames[fMinorVersion-1][lay-1].Data(),
//                   kPathbase.Data(),cpn0,cpn1,cpn2);

        geom->GetGeomMatrix(mod)->SetPath(path);
        GetTransformation(path.Data(),matrix);
        geom->SetTrans(mod,matrix.GetTranslation());
	TransposeTGeoHMatrix(&matrix); //Transpose TGeo's rotation matrixes
        geom->SetRotMatrix(mod,matrix.GetRotationMatrix());
        if(initSeg[kIdet[lay-1]]) continue;
        GetShape(path,shapeName,shapePar);
        if(shapeName.CompareTo("BOX")){
	  Error("InitITSgeomPPRasymmFMD",
		"Geometry changed without proper code update or error "
		"in reading geometry. Shape is not BOX. Shape is %s",
		shapeName.Data());
	  return kFALSE;
        } // end if
	InitGeomShapePPRasymmFMD(kIdet[lay-1],initSeg,shapePar,geom);
    } // end for module
    if(fTiming){
        time->Stop();
        time->Print();
        delete time;
    } // end if
    return kTRUE;
}
//______________________________________________________________________
Bool_t AliITSInitGeometry::InitAliITSgeomV11Hybrid(AliITSgeom *geom){
    // Initilizes the geometry transformation class AliITSgeom
    // to values appropreate to this specific geometry. Now that
    // the segmentation is part of AliITSgeom, the detector
    // segmentations are also defined here.
    // Inputs:
    //   AliITSgeom *geom  A pointer to the AliITSgeom class
    // Outputs:
    //   AliITSgeom *geom  This pointer recreated and properly inilized.
    // Return:
    //   none.

  const Int_t kItype  = 0; // Type of transformation defined 0=> Geant
  const Int_t klayers = 6; // number of layers in the ITS
  const Int_t kladders[klayers]   = {20,40,14,22,34,38}; // Number of ladders
  const Int_t kdetectors[klayers] = {4,4,6,8,22,25};// number of detector/lad
  const AliITSDetector kIdet[6]   = {kSPD,kSPD,kSDD,kSDD,kSSD,kSSD};
  const TString kPathbase = "/ALIC_1/ITSV_1/";

  char *pathSPDsens1, *pathSPDsens2;
  if (SPDIsTGeoNative()) {
    pathSPDsens1="%sITSSPDCarbonFiberSectorV_%d/ITSSPDSensitiveVirtualvolumeM0_1/LAY1_STAVE_%d/HALF-STAVE%d_1/LAY1_LADDER_%d/LAY1_SENSOR_1";
    pathSPDsens2="%sITSSPDCarbonFiberSectorV_%d/ITSSPDSensitiveVirtualvolumeM0_1/LAY2_STAVE_%d/HALF-STAVE%d_1/LAY2_LADDER_%d/LAY2_SENSOR_1";
  } else{
    pathSPDsens1 = "%sITSD_1/IT12_1/I12B_%d/I10B_%d/L1H-STAVE%d_1/I107_%d/I101_1/ITS1_1";
    pathSPDsens2 = "%sITSD_1/IT12_1/I12B_%d/I20B_%d/L2H-STAVE%d_1/I1D7_%d/I1D1_1/ITS2_1";
  }

  char *pathSDDsens1, *pathSDDsens2;
  if (SDDIsTGeoNative()) {
    pathSDDsens1 = "%sITSsddLayer3_1/ITSsddLadd_%d/ITSsddSensor3_%d/ITSsddWafer3_%d/ITSsddSensitivL3_1";
    pathSDDsens2 = "%sITSsddLayer4_1/ITSsddLadd_%d/ITSsddSensor4_%d/ITSsddWafer4_%d/ITSsddSensitivL4_1";
  } else{
    pathSDDsens1 = "%sITSD_1/IT34_1/I004_%d/I302_%d/ITS3_%d";
    pathSDDsens2 = "%sITSD_1/IT34_1/I005_%d/I402_%d/ITS4_%d";
  }

  char *pathSSDsens1, *pathSSDsens2;
  if (SSDIsTGeoNative()) {
    pathSSDsens1 = "%sITSssdLayer5_1/ITSssdLay5Ladd_%d/ITSsddSensor5_%d/ITSsddSensitivL5_1";
    pathSSDsens2 = "%sITSssdLayer6_1/ITSssdLay6Ladd_%d/ITSsddSensor6_%d/ITSsddSensitivL6_1";
  } else{
    pathSSDsens1 = "%sITSD_1/IT56_1/I565_%d/I562_%d/ITS5_%d";
    pathSSDsens2 = "%sITSD_1/IT56_1/I569_%d/I566_%d/ITS6_%d";
  }

  const TString kNames[klayers] = {
    pathSPDsens1, // lay=1
    pathSPDsens2, // lay=2
    pathSDDsens1, // lay=3
    pathSDDsens2, // lay=4
    pathSSDsens1, // lay=5
    pathSSDsens2};// Lay=6
  
  Int_t mod,nmods=0, lay, lad, det, cpn0, cpn1, cpn2, cpnHS=1;
  Double_t tran[3]={0.,0.,0.}, rot[10]={9*0.0,1.0};
  TArrayD shapePar;
  TString path, shapeName;
  TGeoHMatrix matrix;
  Bool_t initSeg[3]={kFALSE, kFALSE, kFALSE};
  TStopwatch *time = 0x0;
  if(fTiming) time = new TStopwatch();

  if(fTiming) time->Start();
  for(mod=0;mod<klayers;mod++) nmods += kladders[mod]*kdetectors[mod];
  geom->Init(kItype,klayers,kladders,kdetectors,nmods);

  for(mod=0; mod<nmods; mod++) {

    DecodeDetectorLayers(mod,lay,lad,det);
    geom->CreateMatrix(mod,lay,lad,det,kIdet[lay-1],tran,rot);
    RecodeDetectorv11Hybrid(mod,cpn0,cpn1,cpn2);

//     if (SPDIsTGeoNative())
//       if (kIdet[lay-1]==kSPD) {
// 	cpn0 = lad-1;
// 	cpn1 = det-1;
// 	cpn2 = 1;
//       }
//     if (SDDIsTGeoNative())
//       if (kIdet[lay-1]==kSDD) {
// 	cpn0 = lad-1;
// 	cpn1 = det-1;
// 	cpn2 = 1;
//       }
//     if (SSDIsTGeoNative())
//       if (kIdet[lay-1]==kSSD) {
// 	cpn0 = lad-1;
// 	cpn1 = det-1;
// 	cpn2 = 1;
//       }

    if (kIdet[lay-1]==kSPD) { // we need 1 more copy number because of the half-stave
      if (det<3) cpnHS = 0; else cpnHS = 1;
      path.Form(kNames[lay-1].Data(),kPathbase.Data(),cpn0,cpn1,cpnHS,cpn2);
    } else {
      path.Form(kNames[lay-1].Data(),kPathbase.Data(),cpn0,cpn1,cpn2);
    };

    geom->GetGeomMatrix(mod)->SetPath(path);
    GetTransformation(path.Data(),matrix);
    geom->SetTrans(mod,matrix.GetTranslation());
    TransposeTGeoHMatrix(&matrix); //Transpose TGeo's rotation matrixes
    geom->SetRotMatrix(mod,matrix.GetRotationMatrix());
    if(initSeg[kIdet[lay-1]]) continue;
    GetShape(path,shapeName,shapePar);
    if(shapeName.CompareTo("BOX")){
      Error("InitITSgeom","Geometry changed without proper code update"
	    "or error in reading geometry. Shape is not BOX.");
      return kFALSE;
    } // end if
    InitGeomShapePPRasymmFMD(kIdet[lay-1],initSeg,shapePar,geom);
  } // end for module

  if(fTiming){
    time->Stop();
    time->Print();
    delete time;
  } // end if
  return kTRUE;
}
//______________________________________________________________________
Bool_t AliITSInitGeometry::InitAliITSgeomV11(AliITSgeom *geom){
  // Initilizes the geometry transformation class AliITSgeom
  // Now that the segmentation is part of AliITSgeom, the detector
  // segmentations are also defined here.
  //
  // Inputs:
  //   AliITSgeom *geom  A pointer to the AliITSgeom class
  // Outputs:
  //   AliITSgeom *geom  This pointer recreated and properly inilized.
  // LG


  const Int_t kItype=0; // Type of transormation defined 0=> Geant
  const Int_t klayers = 6; // number of layers in the ITS
  const Int_t kladders[klayers]   = {20,40,14,22,34,38}; // Number of ladders
  const Int_t kdetectors[klayers] = {4,4,6,8,22,25};// number of detector/lad
  const AliITSDetector kIdet[6]   = {kSPD,kSPD,kSDD,kSDD,kSSD,kSSD};

  const TString kPathbase = "/ALIC_1/ITSV_1/";
  const TString kNames[klayers] =
    {"AliITSInitGeometry:spd missing", // lay=1
     "AliITSInitGeometry:spd missing", // lay=2
     "%sITSsddLayer3_1/ITSsddLadd_%d/ITSsddSensor_%d/ITSsddWafer_1/ITSsddSensitiv_1", // lay=3
     "%sITSsddLayer4_1/ITSsddLadd_%d/ITSsddSensor_%d/ITSsddWafer_1/ITSsddSensitiv_1", // lay=4
     "AliITSInitGeometry:ssd missing", // lay=5
     "AliITSInitGeometry:ssd missing"};// lay=6
 
  Int_t mod,nmods=0,lay,lad,det,cpn0,cpn1,cpn2;
  Double_t tran[3]={0.0,0.0,0.0},rot[10]={9*0.0,1.0};
  TArrayD shapePar;
  TString path,shapeName;
  TGeoHMatrix matrix;
  Bool_t initSeg[3]={kFALSE,kFALSE,kFALSE};
  TStopwatch *time = 0x0;if(fTiming) time=new TStopwatch();
  
  if(fTiming) time->Start();
  for(mod=0;mod<klayers;mod++) nmods += kladders[mod]*kdetectors[mod];
  
  geom->Init(kItype,klayers,kladders,kdetectors,nmods);
  for(mod=0;mod<nmods;mod++) {
    
    DecodeDetectorLayers(mod,lay,lad,det); // Write
    geom->CreateMatrix(mod,lay,lad,det,kIdet[lay-1],tran,rot);
    RecodeDetector(mod,cpn0,cpn1,cpn2); // Write reusing lay,lad,det.
    path.Form(kNames[lay-1].Data(),
	      kPathbase.Data(),cpn0,cpn1,cpn2);
    geom->GetGeomMatrix(mod)->SetPath(path);
    if (GetTransformation(path.Data(),matrix)) {
      geom->SetTrans(mod,matrix.GetTranslation());
      TransposeTGeoHMatrix(&matrix); //Transpose TGeo's rotation matrixes
      geom->SetRotMatrix(mod,matrix.GetRotationMatrix());
    }
    
    if(initSeg[kIdet[lay-1]]) continue;
    GetShape(path,shapeName,shapePar);
    if(shapeName.CompareTo("BOX")){
      Error("InitAliITSgeomV11","Geometry changed without proper code update"
	    "or error in reading geometry. Shape is not BOX.");
      return kFALSE;
    } // end if
    InitGeomShapePPRasymmFMD(kIdet[lay-1],initSeg,shapePar,geom);
    
  } // end for module
  
  if(fTiming){
    time->Stop();
    time->Print();
    delete time;
  } // end if
  return kTRUE;
}

//______________________________________________________________________
Bool_t AliITSInitGeometry::InitGeomShapePPRasymmFMD(AliITSDetector idet,
						       Bool_t *initSeg,
						       TArrayD &shapePar,
						       AliITSgeom *geom){
    // Initilizes the geometry segmentation class AliITSgeomS?D, or
    // AliITSsegmentationS?D depending on the vaule of fSegGeom,
    // to values appropreate to this specific geometry. Now that
    // the segmentation is part of AliITSgeom, the detector
    // segmentations are also defined here.
    // Inputs:
    //   Int_t      lay    The layer number/name.
    //   AliITSgeom *geom  A pointer to the AliITSgeom class
    // Outputs:
    //   AliITSgeom *geom  This pointer recreated and properly inilized.
    // Return:
    //   none.
  //   const Double_t kcm2micron = 1.0E4;
    const Double_t kmicron2cm = 1.0E-4;
    Int_t i;
    TArrayF shapeParF;

    shapeParF.Set(shapePar.GetSize());
    for(i=0;i<shapePar.GetSize();i++) shapeParF[i]=shapePar[i];
    switch (idet){
    case kSPD:{
	initSeg[idet] = kTRUE;
	AliITSgeomSPD *geomSPD = new AliITSgeomSPD425Short();
	Float_t bx[256],bz[280];
	for(i=000;i<256;i++) bx[i] =  50.0*kmicron2cm; // in x all are 50 microns.
	for(i=000;i<160;i++) bz[i] = 425.0*kmicron2cm; // most are 425 microns
	// except below
	for(i=160;i<280;i++) bz[i] =   0.0*kmicron2cm; // Outside of detector.
	bz[ 31] = bz[ 32] = 625.0*kmicron2cm; // first chip boundry
	bz[ 63] = bz[ 64] = 625.0*kmicron2cm; // first chip boundry
	bz[ 95] = bz[ 96] = 625.0*kmicron2cm; // first chip boundry
	bz[127] = bz[128] = 625.0*kmicron2cm; // first chip boundry
	bz[160] = 425.0*kmicron2cm;// Set so that there is no zero pixel size for fNz.
	geomSPD->ReSetBins(shapeParF[1],256,bx,160,bz);
	geom->ReSetShape(idet,geomSPD);
    }break;
    case kSDD:{
	initSeg[idet] = kTRUE;
	AliITSgeomSDD *geomSDD = new AliITSgeomSDD256(shapeParF.GetSize(),
						      shapeParF.GetArray());
	geom->ReSetShape(idet,geomSDD);
    }break;
    case kSSD:{
	initSeg[idet] = kTRUE;
	AliITSgeomSSD *geomSSD = new AliITSgeomSSD275and75(
	    shapeParF.GetSize(),shapeParF.GetArray());
	geom->ReSetShape(idet,geomSSD);
    }break;
    default:{// Others, Note no kSDDp or kSSDp in this geometry.
	geom->ReSetShape(idet,0);
	Info("InitGeomShapePPRasymmFMD",
	     "default Dx=%f Dy=%f Dz=%f default=%d",
	     shapePar[0],shapePar[1],shapePar[2],idet);
    }break;
    } // end switch
    return kTRUE;
}
//______________________________________________________________________
Bool_t AliITSInitGeometry::InitSegmentationPPRasymmFMD(AliITSDetector idet,
						       Bool_t *initSeg,
						       TArrayD &shapePar,
						       AliITSgeom *geom){
    // Initilizes the geometry segmentation class AliITSgeomS?D, or
    // AliITSsegmentationS?D depending on the vaule of fSegGeom,
    // to values appropreate to this specific geometry. Now that
    // the segmentation is part of AliITSgeom, the detector
    // segmentations are also defined here.
    // Inputs:
    //   Int_t      lay    The layer number/name.
    //   AliITSgeom *geom  A pointer to the AliITSgeom class
    // Outputs:
    //   AliITSgeom *geom  This pointer recreated and properly inilized.
    // Return:
    //   none.
    const Double_t kcm2micron = 1.0E4;
    Int_t i;

    switch (idet){
    case kSPD:{
	initSeg[idet] = kTRUE;
	AliITSsegmentationSPD *segSPD = new AliITSsegmentationSPD();
	segSPD->SetDetSize(2.*shapePar[0]*kcm2micron, // X
			   2.*shapePar[2]*kcm2micron, // Z
			   2.*shapePar[1]*kcm2micron);// Y  Microns
	segSPD->SetNPads(256,160);// Number of Bins in x and z
	Float_t bx[256],bz[280];
	for(i=000;i<256;i++) bx[i] =  50.0; // in x all are 50 microns.
	for(i=000;i<160;i++) bz[i] = 425.0; // most are 425 microns
	// except below
	for(i=160;i<280;i++) bz[i] =   0.0; // Outside of detector.
	bz[ 31] = bz[ 32] = 625.0; // first chip boundry
	bz[ 63] = bz[ 64] = 625.0; // first chip boundry
	bz[ 95] = bz[ 96] = 625.0; // first chip boundry
	bz[127] = bz[128] = 625.0; // first chip boundry
	bz[160] = 425.0;// Set so that there is no zero pixel size for fNz.
	segSPD->SetBinSize(bx,bz); // Based on AliITSgeomSPD for now.
	geom->ReSetShape(idet,segSPD);
    }break;
    case kSDD:{
	initSeg[idet] = kTRUE;
	AliITSsegmentationSDD *segSDD = new AliITSsegmentationSDD();
	segSDD->SetDetSize(shapePar[0]*kcm2micron, // X
			   2.*shapePar[2]*kcm2micron, // Z
			   2.*shapePar[1]*kcm2micron);// Y  Microns
	segSDD->SetNPads(256,256);// Anodes, Samples
	geom->ReSetShape(idet,segSDD);
    }break;
    case kSSD:{
	initSeg[idet] = kTRUE;
	AliITSsegmentationSSD *segSSD = new AliITSsegmentationSSD();
	segSSD->SetDetSize(2.*shapePar[0]*kcm2micron, // X
			   2.*shapePar[2]*kcm2micron, // Z
			   2.*shapePar[1]*kcm2micron);// Y  Microns.
	segSSD->SetPadSize(95.,0.); // strip x pitch in microns
	segSSD->SetNPads(768,2); // number of strips on each side, sides.
	segSSD->SetAngles(0.0075,0.0275); // strip angels rad P and N side.
	segSSD->SetAnglesLay5(0.0075,0.0275);//strip angels rad P and N
	segSSD->SetAnglesLay6(0.0275,0.0075);//strip angels rad P and N
	geom->ReSetShape(idet,segSSD);
    }break;
    default:{// Others, Note no kSDDp or kSSDp in this geometry.
	geom->ReSetShape(idet,0);
	Info("InitSegmentationPPRasymmFMD",
	     "default segmentation Dx=%f Dy=%f Dz=%f default=%d",
	     shapePar[0],shapePar[1],shapePar[2],idet);
    }break;
    } // end switch
    return kTRUE;
}
//______________________________________________________________________
Bool_t AliITSInitGeometry::GetTransformation(const TString &volumePath,
					     TGeoHMatrix &mat){
    // Returns the Transformation matrix between the volume specified
    // by the path volumePath and the Top or mater volume. The format
    // of the path volumePath is as follows (assuming ALIC is the Top volume)
    // "/ALIC_1/DDIP_1/S05I_2/S05H_1/S05G_3". Here ALIC is the top most
    // or master volume which has only 1 instance of. Of all of the daughter
    // volumes of ALICE, DDIP volume copy #1 is indicated. Similarly for
    // the daughter volume of DDIP is S05I copy #2 and so on.
    // Inputs:
    //   TString& volumePath  The volume path to the specific volume
    //                        for which you want the matrix. Volume name
    //                        hierarchy is separated by "/" while the
    //                        copy number is appended using a "_".
    // Outputs:
    //  TGeoHMatrix &mat      A matrix with its values set to those
    //                        appropriate to the Local to Master transformation
    // Return:
    //   A logical value if kFALSE then an error occurred and no change to
    //   mat was made.

    // We have to preserve the modeler state

    // Preserve the modeler state.
    gGeoManager->PushPath();
    if (!gGeoManager->cd(volumePath.Data())) {
      gGeoManager->PopPath();
      Error("GetTransformation","Error in cd-ing to ",volumePath.Data());
      return kFALSE;
    } // end if !gGeoManager
    mat = *gGeoManager->GetCurrentMatrix();
    // Retstore the modeler state.
    gGeoManager->PopPath();
    return kTRUE;
}
//______________________________________________________________________
Bool_t AliITSInitGeometry::GetShape(const TString &volumePath,
				    TString &shapeType,TArrayD &par){
    // Returns the shape and its parameters for the volume specified
    // by volumeName.
    // Inputs:
    //   TString& volumeName  The volume name
    // Outputs:
    //   TString &shapeType   Shape type
    //   TArrayD &par         A TArrayD of parameters with all of the
    //                        parameters of the specified shape.
    // Return:
    //   A logical indicating whether there was an error in getting this
    //   information
    Int_t npar;
    gGeoManager->PushPath();
    if (!gGeoManager->cd(volumePath.Data())) {
	gGeoManager->PopPath();
	return kFALSE;
    }
    TGeoVolume * vol = gGeoManager->GetCurrentVolume();
    gGeoManager->PopPath();
    if (!vol) return kFALSE;
    TGeoShape *shape = vol->GetShape();
    TClass *classType = shape->IsA();
    if (classType==TGeoBBox::Class()) {
	shapeType = "BOX";
	npar = 3;
	par.Set(npar);
	TGeoBBox *box = (TGeoBBox*)shape;
	par.AddAt(box->GetDX(),0);
	par.AddAt(box->GetDY(),1);
	par.AddAt(box->GetDZ(),2);
	return kTRUE;
    } // end if
    if (classType==TGeoTrd1::Class()) {
	shapeType = "TRD1";
	npar = 4;
	par.Set(npar);
	TGeoTrd1 *trd1 = (TGeoTrd1*)shape;
	par.AddAt(trd1->GetDx1(),0);
	par.AddAt(trd1->GetDx2(),1);
	par.AddAt(trd1->GetDy(), 2);
	par.AddAt(trd1->GetDz(), 3);
	return kTRUE;
    } // end if
    if (classType==TGeoTrd2::Class()) {
	shapeType = "TRD2";
	npar = 5;
	par.Set(npar);
	TGeoTrd2 *trd2 = (TGeoTrd2*)shape;
	par.AddAt(trd2->GetDx1(),0);
	par.AddAt(trd2->GetDx2(),1);
	par.AddAt(trd2->GetDy1(),2);
	par.AddAt(trd2->GetDy2(),3);
	par.AddAt(trd2->GetDz(), 4);
	return kTRUE;
    } // end if
    if (classType==TGeoTrap::Class()) {
	shapeType = "TRAP";
	npar = 11;
	par.Set(npar);
	TGeoTrap *trap = (TGeoTrap*)shape;
	Double_t tth = TMath::Tan(trap->GetTheta()*TMath::DegToRad());
	par.AddAt(trap->GetDz(),0);
	par.AddAt(tth*TMath::Cos(trap->GetPhi()*TMath::DegToRad()),1);
	par.AddAt(tth*TMath::Sin(trap->GetPhi()*TMath::DegToRad()),2);
	par.AddAt(trap->GetH1(),3);
	par.AddAt(trap->GetBl1(),4);
	par.AddAt(trap->GetTl1(),5);
	par.AddAt(TMath::Tan(trap->GetAlpha1()*TMath::DegToRad()),6);
	par.AddAt(trap->GetH2(),7);
	par.AddAt(trap->GetBl2(),8);
	par.AddAt(trap->GetTl2(),9);
	par.AddAt(TMath::Tan(trap->GetAlpha2()*TMath::DegToRad()),10);
	return kTRUE;
    } // end if
    if (classType==TGeoTube::Class()) {
	shapeType = "TUBE";
	npar = 3;
	par.Set(npar);
	TGeoTube *tube = (TGeoTube*)shape;
	par.AddAt(tube->GetRmin(),0);
	par.AddAt(tube->GetRmax(),1);
	par.AddAt(tube->GetDz(),2);
	return kTRUE;
    } // end if
    if (classType==TGeoTubeSeg::Class()) {
	shapeType = "TUBS";
	npar = 5;
	par.Set(npar);
	TGeoTubeSeg *tubs = (TGeoTubeSeg*)shape;
	par.AddAt(tubs->GetRmin(),0);
	par.AddAt(tubs->GetRmax(),1);
	par.AddAt(tubs->GetDz(),2);
	par.AddAt(tubs->GetPhi1(),3);
	par.AddAt(tubs->GetPhi2(),4);
	return kTRUE;
    } // end if
    if (classType==TGeoCone::Class()) {
	shapeType = "CONE";
	npar = 5;
	par.Set(npar);
	TGeoCone *cone = (TGeoCone*)shape;
	par.AddAt(cone->GetDz(),0);
	par.AddAt(cone->GetRmin1(),1);
	par.AddAt(cone->GetRmax1(),2);
	par.AddAt(cone->GetRmin2(),3);
	par.AddAt(cone->GetRmax2(),4);
	return kTRUE;
    } // end if
    if (classType==TGeoConeSeg::Class()) {
	shapeType = "CONS";
	npar = 7;
	par.Set(npar);
	TGeoConeSeg *cons = (TGeoConeSeg*)shape;
	par.AddAt(cons->GetDz(),0);
	par.AddAt(cons->GetRmin1(),1);
	par.AddAt(cons->GetRmax1(),2);
	par.AddAt(cons->GetRmin2(),3);
	par.AddAt(cons->GetRmax2(),4);
	par.AddAt(cons->GetPhi1(),5);
	par.AddAt(cons->GetPhi2(),6);
	return kTRUE;
    } // end if
    if (classType==TGeoSphere::Class()) {
	shapeType = "SPHE";
	npar = 6;
	par.Set(npar);
	
	TGeoSphere *sphe = (TGeoSphere*)shape;
	par.AddAt(sphe->GetRmin(),0);
	par.AddAt(sphe->GetRmax(),1);
	par.AddAt(sphe->GetTheta1(),2);
	par.AddAt(sphe->GetTheta2(),3);
	par.AddAt(sphe->GetPhi1(),4);
	par.AddAt(sphe->GetPhi2(),5);
	return kTRUE;
    } // end if
    if (classType==TGeoPara::Class()) {
	shapeType = "PARA";
	npar = 6;
	par.Set(npar);
	TGeoPara *para = (TGeoPara*)shape;
	par.AddAt(para->GetX(),0);
	par.AddAt(para->GetY(),1);
	par.AddAt(para->GetZ(),2);
	par.AddAt(para->GetTxy(),3);
	par.AddAt(para->GetTxz(),4);
	par.AddAt(para->GetTyz(),5);
	return kTRUE;
    } // end if
    if (classType==TGeoPgon::Class()) {
	shapeType = "PGON";
	TGeoPgon *pgon = (TGeoPgon*)shape;
	Int_t nz = pgon->GetNz();
	const Double_t *rmin = pgon->GetRmin();
	const Double_t *rmax = pgon->GetRmax();
	const Double_t *z = pgon->GetZ();
	npar = 4 + 3*nz;
	par.Set(npar);
	par.AddAt(pgon->GetPhi1(),0);
	par.AddAt(pgon->GetDphi(),1);
	par.AddAt(pgon->GetNedges(),2);
	par.AddAt(pgon->GetNz(),3);
	for (Int_t i=0; i<nz; i++) {
	    par.AddAt(z[i], 4+3*i);
	    par.AddAt(rmin[i], 4+3*i+1);
	    par.AddAt(rmax[i], 4+3*i+2);
	}
	return kTRUE;
    } // end if
    if (classType==TGeoPcon::Class()) {
	shapeType = "PCON";
	TGeoPcon *pcon = (TGeoPcon*)shape;
	Int_t nz = pcon->GetNz();
	const Double_t *rmin = pcon->GetRmin();
	const Double_t *rmax = pcon->GetRmax();
	const Double_t *z = pcon->GetZ();
	npar = 3 + 3*nz;
	par.Set(npar);
	par.AddAt(pcon->GetPhi1(),0);
	par.AddAt(pcon->GetDphi(),1);
	par.AddAt(pcon->GetNz(),2);
	for (Int_t i=0; i<nz; i++) {
	    par.AddAt(z[i], 3+3*i);
	    
	    par.AddAt(rmin[i], 3+3*i+1);
	    par.AddAt(rmax[i], 3+3*i+2);
	}
	return kTRUE;
    } // end if
    if (classType==TGeoEltu::Class()) {
	shapeType = "ELTU";
	npar = 3;
	par.Set(npar);
	TGeoEltu *eltu = (TGeoEltu*)shape;
	par.AddAt(eltu->GetA(),0);
	par.AddAt(eltu->GetB(),1);
	par.AddAt(eltu->GetDz(),2);
	return kTRUE;
    } // end if
    if (classType==TGeoHype::Class()) {
	shapeType = "HYPE";
	npar = 5;
	par.Set(npar);
	TGeoHype *hype = (TGeoHype*)shape;
	par.AddAt(TMath::Sqrt(hype->RadiusHypeSq(0.,kTRUE)),0);
	par.AddAt(TMath::Sqrt(hype->RadiusHypeSq(0.,kFALSE)),1);
	par.AddAt(hype->GetDZ(),2);
	par.AddAt(hype->GetStIn(),3);
	par.AddAt(hype->GetStOut(),4);
	return kTRUE;
    } // end if
    if (classType==TGeoGtra::Class()) {
	shapeType = "GTRA";
	npar = 12;
	par.Set(npar);
	TGeoGtra *trap = (TGeoGtra*)shape;
	Double_t tth = TMath::Tan(trap->GetTheta()*TMath::DegToRad());
	par.AddAt(trap->GetDz(),0);
	par.AddAt(tth*TMath::Cos(trap->GetPhi()*TMath::DegToRad()),1);
	par.AddAt(tth*TMath::Sin(trap->GetPhi()*TMath::DegToRad()),2);
	par.AddAt(trap->GetH1(),3);
	par.AddAt(trap->GetBl1(),4);
	par.AddAt(trap->GetTl1(),5);
	par.AddAt(TMath::Tan(trap->GetAlpha1()*TMath::DegToRad()),6);
	par.AddAt(trap->GetH2(),7);
	par.AddAt(trap->GetBl2(),8);
	par.AddAt(trap->GetTl2(),9);
	par.AddAt(TMath::Tan(trap->GetAlpha2()*TMath::DegToRad()),10);
	par.AddAt(trap->GetTwistAngle(),11);
	return kTRUE;
    } // end if
    if (classType==TGeoCtub::Class()) {
	shapeType = "CTUB";
	npar = 11;
	par.Set(npar);
	TGeoCtub *ctub = (TGeoCtub*)shape;
	const Double_t *lx = ctub->GetNlow();
	const Double_t *tx = ctub->GetNhigh();
	par.AddAt(ctub->GetRmin(),0);
	par.AddAt(ctub->GetRmax(),1);
	par.AddAt(ctub->GetDz(),2);
	par.AddAt(ctub->GetPhi1(),3);
	par.AddAt(ctub->GetPhi2(),4);
	par.AddAt(lx[0],5);
	par.AddAt(lx[1],6);
	par.AddAt(lx[2],7);
	par.AddAt(tx[0],8);
	par.AddAt(tx[1],9);
	par.AddAt(tx[2],10);
	return kTRUE;
    } // end if
    Error("GetShape","Getting shape parameters for shape %s not implemented",
	  shape->ClassName());
    shapeType = "Unknown";
    return kFALSE;
}
//______________________________________________________________________
void AliITSInitGeometry::DecodeDetector(
    Int_t &mod,Int_t layer,Int_t cpn0,Int_t cpn1,Int_t cpn2) const {
    // decode geometry into detector module number. There are two decoding
    // Scheams. Old which does not follow the ALICE coordinate system
    // requirements, and New which dose.
    // Inputs:
    //    Int_t layer    The ITS layer
    //    Int_t cpn0     The lowest copy number
    //    Int_t cpn1     The middle copy number
    //    Int_t cpn2     the highest copy number
    // Output:
    //    Int_t &mod     The module number assoicated with this set
    //                   of copy numbers.
    // Return:
    //    none.

    // This is a FIXED switch yard function. I (Bjorn Nilsen) Don't 
    // like them but I see not better way for the moment.
    switch (fMajorVersion){
    case kvtest:{
        if(GetMinorVersion()==1)
            return DecodeDetectorvPPRasymmFMD(mod,layer,cpn0,cpn1,cpn2);
        else if(GetMinorVersion()==2)
            return DecodeDetectorvtest2(mod,layer,cpn0,cpn1,cpn2);
        Warning("DecodeDetector",
                "Geometry is kvtest minor version=%d is not defined",
                GetMinorVersion());
    }break;
    case kvDefault:{
        Error("DecodeDetector","Major version = kvDefault, not supported");
    }break;
    case kvSPD02:{
        return DecodeDetectorvSPD02(mod,layer,cpn0,cpn1,cpn2);
    }break;
    case kvSDD03:{
        return DecodeDetectorvSDD03(mod,layer,cpn0,cpn1,cpn2);
    }break;
    case kvSSD03:{
        return DecodeDetectorvSSD03(mod,layer,cpn0,cpn1,cpn2);
    }break;
    case kvITS04:{
        return DecodeDetectorvITS04(mod,layer,cpn0,cpn1,cpn2);
    }break;
    case kvPPRcourseasymm:{
        return DecodeDetectorvPPRcourseasymm(mod,layer,cpn0,cpn1,cpn2);
    }break;
    case kvPPRasymmFMD:{
        return DecodeDetectorvPPRasymmFMD(mod,layer,cpn0,cpn1,cpn2);
    }break;
    case kv11:{
        return DecodeDetectorv11(mod,layer,cpn0,cpn1,cpn2);
    }break;
    case kv11Hybrid:{
        return DecodeDetectorv11Hybrid(mod,layer,cpn0,cpn1,cpn2);
    }break;
    default:{
        Error("DecodeDetector","Major version = %d, not supported",
              (Int_t)fMajorVersion);
        return;
    }break;
    } // end switch
    return;
}
//______________________________________________________________________
void AliITSInitGeometry::RecodeDetector(Int_t mod,Int_t &cpn0,
                                        Int_t &cpn1,Int_t &cpn2){
    // decode geometry into detector module number. There are two decoding
    // Scheams. Old which does not follow the ALICE coordinate system
    // requirements, and New which dose.
    // Inputs:
    //    Int_t mod      The module number assoicated with this set
    //                   of copy numbers.
    // Output:
    //    Int_t cpn0     The lowest copy number
    //    Int_t cpn1     The middle copy number
    //    Int_t cpn2     the highest copy number
    // Return:
    //    none.

    // This is a FIXED switch yard function. I (Bjorn Nilsen) Don't 
    // like them but I see not better way for the moment.
    switch (fMajorVersion){
    case kvtest:{
        if(GetMinorVersion()==1) 
            return RecodeDetectorvPPRasymmFMD(mod,cpn0,cpn1,cpn2);
        else if(GetMinorVersion()==2)
            return RecodeDetectorvtest2(mod,cpn0,cpn1,cpn2);
        Warning("RecodeDetector",
                "Geometry is kvtest minor version=%d is not defined",
                GetMinorVersion());
        return;
    }break;
    case kvDefault:{
        Error("RecodeDetector","Major version = kvDefault, not supported");
        return;
    }break;
    case kvSPD02:{
        return RecodeDetectorvSPD02(mod,cpn0,cpn1,cpn2);
    }break;
    case kvSDD03:{
        return RecodeDetectorvSDD03(mod,cpn0,cpn1,cpn2);
    }break;
    case kvSSD03:{
        return RecodeDetectorvSSD03(mod,cpn0,cpn1,cpn2);
    }break;
    case kvITS04:{
        return RecodeDetectorvITS04(mod,cpn0,cpn1,cpn2);
    }break;
    case kvPPRcourseasymm:{
        return RecodeDetectorvPPRcourseasymm(mod,cpn0,cpn1,cpn2);
    }break;
    case kvPPRasymmFMD:{
        return RecodeDetectorvPPRasymmFMD(mod,cpn0,cpn1,cpn2);
    }break;
    case kv11:{
        return RecodeDetectorv11(mod,cpn0,cpn1,cpn2);
    }break;
    case kv11Hybrid:{
        return RecodeDetectorv11Hybrid(mod,cpn0,cpn1,cpn2);
    }break;
    default:{
        Error("RecodeDetector","Major version = %d, not supported",
              (Int_t)fMajorVersion);
        return;
    }break;
    } // end switch
    return;
}
//______________________________________________________________________
void AliITSInitGeometry::DecodeDetectorLayers(Int_t mod,Int_t &layer,
                                              Int_t &lad,Int_t &det){
    // decode geometry into detector module number. There are two decoding
    // Scheams. Old which does not follow the ALICE coordinate system
    // requirements, and New which dose. Note, this use of layer ladder
    // and detector numbers are strictly for internal use of this
    // specific code. They do not represent the "standard" layer ladder
    // or detector numbering except in a very old and obsoleate sence.
    // Inputs:
    //    Int_t mod      The module number assoicated with this set
    //                   of copy numbers.
    // Output:
    //    Int_t lay     The layer number
    //    Int_t lad     The ladder number
    //    Int_t det     the dettector number
    // Return:
    //    none.

    // This is a FIXED switch yard function. I (Bjorn Nilsen) Don't 
    // like them but I see not better way for the moment.
    switch (fMajorVersion) {
    case kvtest:{
        if(GetMinorVersion()==1) 
            return DecodeDetectorLayersvPPRasymmFMD(mod,layer,lad,det);
        else if(GetMinorVersion()==2)
            return DecodeDetectorLayersvtest2(mod,layer,lad,det);
        Warning("DecodeDetectorLayers",
                "Geometry is kvtest minor version=%d is not defined",
                GetMinorVersion());
        return;
    } break;
    case kvDefault:{
        Error("DecodeDetectorLayers",
              "Major version = kvDefault, not supported");
        return;
    }break;
    case kvSPD02:{
        return DecodeDetectorLayersvSPD02(mod,layer,lad,det);
    }break;
    case kvSDD03:{
        return DecodeDetectorLayersvSDD03(mod,layer,lad,det);
    }break;
    case kvSSD03:{
        return DecodeDetectorLayersvSSD03(mod,layer,lad,det);
    }break;
    case kvITS04:{
        return DecodeDetectorLayersvITS04(mod,layer,lad,det);
    }break;
    case kvPPRcourseasymm:{
        return DecodeDetectorLayersvPPRcourseasymm(mod,layer,lad,det);
    }break;
    case kvPPRasymmFMD:{
        return DecodeDetectorLayersvPPRasymmFMD(mod,layer,lad,det);
    }break;
    case kv11:{
        return DecodeDetectorLayersv11(mod,layer,lad,det);
    }break;
    case kv11Hybrid:{
        return DecodeDetectorLayersv11Hybrid(mod,layer,lad,det);
    }break;
    default:{
        Error("DecodeDetectorLayers","Major version = %d, not supported",
              (Int_t)fMajorVersion);
        return;
    }break;
    } // end switch
    return;
}
//______________________________________________________________________
void AliITSInitGeometry::DecodeDetectorvSPD02(
    Int_t &mod,Int_t ncpn,Int_t cpy0,Int_t cpy1,Int_t cpy2) const {
    // decode geometry into detector module number
    // Inputs:
    //    Int_t ncpn     The Number of copies of this volume
    //    Int_t cpy0     The lowest copy number
    //    Int_t cpy1     The middle copy number
    //    Int_t cpy2     the highest copy number
    // Output:
    //    Int_t &mod     The module number assoicated with this set
    //                   of copy numbers.
    // Return:
    //    none.

    // detector = ladder = 1
    if(ncpn==4 && cpy1>2) mod = cpy1; // layer = 1,2
    else mod = cpy1-1; // layer = 4,5
    if(ncpn==1) mod = 2; // layer=3
    cpy0 = cpy2;
    return;
}
//______________________________________________________________________
void AliITSInitGeometry::RecodeDetectorvSPD02(Int_t mod,Int_t &cpn0,
                                        Int_t &cpn1,Int_t &cpn2){
    // decode geometry into detector module number. There are two decoding
    // Scheams. Old which does not follow the ALICE coordinate system
    // requirements, and New which dose.
    // Inputs:
    //    Int_t mod      The module number assoicated with this set
    //                   of copy numbers.
    // Output:
    //    Int_t cpn0     The lowest copy number
    //    Int_t cpn1     The middle copy number
    //    Int_t cpn2     the highest copy number
    // Return:
    //    none.

    cpn2 = 0;
    if(mod==2){
        cpn0 = 1;
        cpn1 = 1;
        return;
    } else if(mod<2){
        cpn0 = 1;
        cpn1 = mod+1;
    }else{
        cpn0 = 1;
        cpn1 = mod;
    } // end if
    return;
}
//______________________________________________________________________
void AliITSInitGeometry::DecodeDetectorLayersvSPD02(Int_t mod,Int_t &lay,
                                                    Int_t &lad,Int_t &det){
    // decode geometry into detector module number. There are two decoding
    // Scheams. Old which does not follow the ALICE coordinate system
    // requirements, and New which dose. Note, this use of layer ladder
    // and detector numbers are strictly for internal use of this
    // specific code. They do not represent the "standard" layer ladder
    // or detector numbering except in a very old and obsoleate sence.
    // Inputs:
    //    Int_t mod      The module number assoicated with this set
    //                   of copy numbers.
    // Output:
    //    Int_t lay     The layer number
    //    Int_t lad     The ladder number
    //    Int_t det     the dettector number
    // Return:
    //    none.

    lay = mod+1;
    lad = det = 1;
    return;
}
//______________________________________________________________________
void AliITSInitGeometry::DecodeDetectorvSDD03(
    Int_t &mod,Int_t ncpys,Int_t cpy0,Int_t cpy1,Int_t cpy2) const {
    // decode geometry into detector module number. There are two decoding
    // Scheams. Old which does not follow the ALICE coordinate system
    // requirements, and New which dose.
    // Inputs:
    //    Int_t ncpys    The number of posible copies cpn1
    //    Int_t cpy0     The lowest copy number
    //    Int_t cpy1     The middle copy number
    //    Int_t cpy2     the highest copy number
    // Output:
    //    Int_t &mod     The module number assoicated with this set
    //                   of copy numbers.
    // Return:
    //    none.

    if(ncpys==10){ // ITEL detectors
        if(cpy1>4) mod = cpy1+1;
        else mod = cpy1-1;
    }else{ // IDET detectors
        if(cpy1==1) mod = 4;
        else mod = 5;
    } // end if
    cpy0=cpy2;
    return;
}
//______________________________________________________________________
void AliITSInitGeometry::RecodeDetectorvSDD03(Int_t mod,Int_t &cpn0,
                                        Int_t &cpn1,Int_t &cpn2){
    // decode geometry into detector module number. There are two decoding
    // Scheams. Old which does not follow the ALICE coordinate system
    // requirements, and New which dose.
    // Inputs:
    //    Int_t mod      The module number assoicated with this set
    //                   of copy numbers.
    // Output:
    //    Int_t cpn0     The lowest copy number
    //    Int_t cpn1     The middle copy number
    //    Int_t cpn2     the highest copy number
    // Return:
    //    none.

    cpn0 = 1;
    cpn2 = 0;
    if(mod<4) cpn1 = mod+1;
    else if(mod==4||mod==5) cpn1 = mod-3;
    else cpn1 = mod-1;
    return;
}
//______________________________________________________________________
void AliITSInitGeometry::DecodeDetectorLayersvSDD03(Int_t mod,Int_t &lay,
                                                    Int_t &lad,Int_t &det){
    // decode geometry into detector module number. There are two decoding
    // Scheams. Old which does not follow the ALICE coordinate system
    // requirements, and New which dose. Note, this use of layer ladder
    // and detector numbers are strictly for internal use of this
    // specific code. They do not represent the "standard" layer ladder
    // or detector numbering except in a very old and obsoleate sence.
    // Inputs:
    //    Int_t mod      The module number assoicated with this set
    //                   of copy numbers.
    // Output:
    //    Int_t lay     The layer number
    //    Int_t lad     The ladder number
    //    Int_t det     the dettector number
    // Return:
    //    none.

    lad = det = 1;
    lay = mod+1;
    return;
}
//______________________________________________________________________
void AliITSInitGeometry::DecodeDetectorvSSD03(
    Int_t &mod,Int_t dtype,Int_t cpn0,Int_t cpn1,Int_t cpn2) const {
    // decode geometry into detector module number. There are two decoding
    // Scheams. Old which does not follow the ALICE coordinate system
    // requirements, and New which dose.
    // Inputs:
    //    Int_t dtype    The detector type 1=ITSA 2=IGAR 3=IFRA
    //    Int_t cpn0     The lowest copy number
    //    Int_t cpn1     The middle copy number
    //    Int_t cpn2     the highest copy number
    // Output:
    //    Int_t &mod     The module number assoicated with this set
    //                   of copy numbers.
    // Return:
    //    none.

    if(dtype==2){mod=2; return;}
    if(dtype==3){mod=3; return;}
    mod = cpn0-1;
    if(cpn0==3) mod = 4;
    cpn1=cpn2;
    return;
}
//______________________________________________________________________
void AliITSInitGeometry::RecodeDetectorvSSD03(Int_t mod,Int_t &cpn0,
                                        Int_t &cpn1,Int_t &cpn2){
    // decode geometry into detector module number. There are two decoding
    // Scheams. Old which does not follow the ALICE coordinate system
    // requirements, and New which dose.
    // Inputs:
    //    Int_t mod      The module number assoicated with this set
    //                   of copy numbers.
    // Output:
    //    Int_t cpn0     The lowest copy number
    //    Int_t cpn1     The middle copy number
    //    Int_t cpn2     the highest copy number
    // Return:
    //    none.

    cpn1=1;
    cpn2=0;
    if(mod<2) cpn0=mod+1;
    else if (mod==2||mod==3) cpn0=1;
    else cpn0 = 3;
    return;
}
//______________________________________________________________________
void AliITSInitGeometry::DecodeDetectorLayersvSSD03(Int_t mod,Int_t &lay,
                                                    Int_t &lad,Int_t &det){
    // decode geometry into detector module number. There are two decoding
    // Scheams. Old which does not follow the ALICE coordinate system
    // requirements, and New which dose. Note, this use of layer ladder
    // and detector numbers are strictly for internal use of this
    // specific code. They do not represent the "standard" layer ladder
    // or detector numbering except in a very old and obsoleate sence.
    // Inputs:
    //    Int_t mod      The module number assoicated with this set
    //                   of copy numbers.
    // Output:
    //    Int_t lay     The layer number
    //    Int_t lad     The ladder number
    //    Int_t det     the dettector number
    // Return:
    //    none.

    lad = det = 1;
    lay = mod+1;
    return;
}
//______________________________________________________________________
void AliITSInitGeometry::DecodeDetectorvITS04(
    Int_t &mod,Int_t dtype,Int_t cpn0,Int_t cpn1,Int_t cpn2) const {
    // decode geometry into detector module number. There are two decoding
    // Scheams. Old which does not follow the ALICE coordinate system
    // requirements, and New which dose.
    // Inputs:
    //    Int_t dtype    The detector type 1=ITSA 2=IGAR 3=IFRA
    //    Int_t cpn0     The lowest copy number
    //    Int_t cpn1     The middle copy number
    //    Int_t cpn2     the highest copy number
    // Output:
    //    Int_t &mod     The module number assoicated with this set
    //                   of copy numbers.
    // Return:
    //    none.

    mod = dtype-1;
    cpn0 = cpn1 = cpn2;
    return;
}
//______________________________________________________________________
void AliITSInitGeometry::RecodeDetectorvITS04(Int_t mod,Int_t &cpn0,
                                        Int_t &cpn1,Int_t &cpn2){
    // decode geometry into detector module number. There are two decoding
    // Scheams. Old which does not follow the ALICE coordinate system
    // requirements, and New which dose.
    // Inputs:
    //    Int_t mod      The module number assoicated with this set
    //                   of copy numbers.
    // Output:
    //    Int_t cpn0     The lowest copy number
    //    Int_t cpn1     The middle copy number
    //    Int_t cpn2     the highest copy number
    // Return:
    //    none.

    cpn1 = cpn2 = 0;
    switch(mod){
    case 0:case 1:case 2:case 3:{
        cpn0 = mod+1;
    }break;
    case 4: case 5:{
        cpn0 = mod-3;
    }break;
    case 6:case 7:case 8:case 9:{
        cpn0 = mod-5;
    } break;
    default:
        cpn0 = 0;
        break;
    }// end switch
    return;
}
//______________________________________________________________________
void AliITSInitGeometry::DecodeDetectorLayersvITS04(Int_t mod,Int_t &lay,
                                                    Int_t &lad,Int_t &det){
    // decode geometry into detector module number. There are two decoding
    // Scheams. Old which does not follow the ALICE coordinate system
    // requirements, and New which dose. Note, this use of layer ladder
    // and detector numbers are strictly for internal use of this
    // specific code. They do not represent the "standard" layer ladder
    // or detector numbering except in a very old and obsoleate sence.
    // Inputs:
    //    Int_t mod      The module number assoicated with this set
    //                   of copy numbers.
    // Output:
    //    Int_t lay     The layer number
    //    Int_t lad     The ladder number
    //    Int_t det     the dettector number
    // Return:
    //    none.

    lad = 1;
    switch(mod){
    case 0:case 1:case 2:case 3:{
        lay = mod/2 +1;
        det = mod%2 +1;
    }break;
    case 4: case 5:{
        lay = mod -1;
    }break;
    case 6:case 7:case 8:case 9:{
        lay = mod/2 +2;
        det = mod%2 +1;
    }break;
    default:
        lay = 0;
        det = 0;
        break;
    } // end switch
    return;
}
//______________________________________________________________________
void AliITSInitGeometry::DecodeDetectorvPPRasymmFMD(Int_t &mod,Int_t layer,Int_t cpn0,
                                        Int_t cpn1,Int_t cpn2) const {
    // decode geometry into detector module number. There are two decoding
    // Scheams. Old which does not follow the ALICE coordinate system
    // requirements, and New which dose.
    // Inputs:
    //    Int_t layer    The ITS layer
    //    Int_t cpn0     The lowest copy number
    //    Int_t cpn1     The middle copy number
    //    Int_t cpn2     the highest copy number
    // Output:
    //    Int_t &mod     The module number assoicated with this set
    //                   of copy numbers.
    // Return:
    //    none.
    const Int_t kDetPerLadderSPD[2]={2,4};
    const Int_t kDetPerLadder[6]={4,4,6,8,22,25};
    const Int_t kLadPerLayer[6]={20,40,14,22,34,38};
    Int_t lay=-1,lad=-1,det=-1,i;

    if(fDecode){ // New decoding scheam
        switch (layer){
        case 1:{
            lay = layer;
            det = 5-cpn2;
            if(cpn0==4&&cpn1==1) lad=1;
            else if(cpn0==4&&cpn1==2) lad=20;
            else if(cpn0<4){
                lad = 8-cpn1-kDetPerLadderSPD[layer-1]*(cpn0-1);
            }else{ // cpn0>4
                lad = 28-cpn1-kDetPerLadderSPD[layer-1]*(cpn0-1);
            } // end if
        } break;
        case 2:{
            lay = layer;
            det = 5-cpn2;
            if(cpn0==4&&cpn1==1) lad=1;
            else if(cpn0<4){
                lad = 14-cpn1-kDetPerLadderSPD[layer-1]*(cpn0-1);
            }else{ // cpn0>4
                lad = 54-cpn1-kDetPerLadderSPD[layer-1]*(cpn0-1);
            } // end if
        } break;
        case 3:{
            lay = layer;
            if(cpn0<5) lad = 5-cpn0;
            else lad = 19-cpn0;
            det = 7-cpn1;
        } break;
        case 4:{
            lay = layer;
            if(cpn0<7) lad = 7-cpn0;
            else lad = 29-cpn0;
            det = 9-cpn1;
        } break;
        case 5:{
            lay = layer;
            if(cpn0<10) lad = 10-cpn0;
            else lad = 44-cpn0;
            det = 23-cpn1;
        } break;
        case 6:{
            lay = layer;
            if(cpn0<9) lad = 9-cpn0;
            else lad = 47-cpn0;
            det = 26-cpn1;
        } break;
        } // end switch
        mod = 0;
        for(i=0;i<layer-1;i++) mod += kLadPerLayer[i]*kDetPerLadder[i];
        mod += kDetPerLadder[layer-1]*(lad-1)+det-1;// module start at zero.
        return;
    } // end if
    // Old decoding scheam
    switch(layer){
    case 1: case 2:{
        lay = layer;
        lad = cpn1+kDetPerLadderSPD[layer-1]*(cpn0-1);
        det = cpn2;
        }break;
    case 3: case 4:{
        lay = layer;
        lad = cpn0;
        det = cpn1;
        }break;
    case 5: case 6:{
        lay = layer;
        lad = cpn0;
        det = cpn1;
        }break;
    default:{
        }break;
    } // end switch
    mod = 0;
    for(i=0;i<layer-1;i++) mod += kLadPerLayer[i]*kDetPerLadder[i];
    mod += kDetPerLadder[layer-1]*(lad-1)+det-1;// module start at zero.
    return;
}
//______________________________________________________________________
void AliITSInitGeometry::RecodeDetectorvPPRasymmFMD(Int_t mod,Int_t &cpn0,
                                        Int_t &cpn1,Int_t &cpn2){
    // decode geometry into detector module number. There are two decoding
    // Scheams. Old which does not follow the ALICE coordinate system
    // requirements, and New which dose.
    // Inputs:
    //    Int_t mod      The module number assoicated with this set
    //                   of copy numbers.
    // Output:
    //    Int_t cpn0     The lowest copy number
    //    Int_t cpn1     The middle copy number
    //    Int_t cpn2     the highest copy number
    // Return:
    //    none.
    const Int_t kITSgeoTreeCopys[6][3]= {{10, 2, 4},// lay=1
                                         {10, 4, 4},// lay=2
                                         {14, 6, 1},// lay=3
                                         {22, 8, 1},// lay=4
                                         {34,22, 1},// lay=5
                                         {38,25, 1}};//lay=6
    const Int_t kDetPerLadderSPD[2]={2,4};
    //    const Int_t kDetPerLadder[6]={4,4,6,8,22,25};
    //    const Int_t kLadPerLayer[6]={20,40,14,22,34,38};
    Int_t lay,lad,det;

    cpn0 = cpn1 = cpn2 = 0;
    DecodeDetectorLayers(mod,lay,lad,det);
    if(fDecode){ // New decoding scheam
        switch (lay){
        case 1:{
            cpn2 = 5-det;     // Detector 1-4
            cpn1 = 1+(lad-1)%kDetPerLadderSPD[lay-1];
            cpn0 = 5-(lad+kDetPerLadderSPD[lay-1])/kDetPerLadderSPD[lay-1];
            if(mod>27) cpn0 = 15-(lad+kDetPerLadderSPD[lay-1])/
			   kDetPerLadderSPD[lay-1];
        } break;
        case 2:{
            cpn2 = 5-det;     // Detector 1-4
            cpn1 = 4-(lad+2)%kDetPerLadderSPD[lay-1];
            cpn0 = 1+(14-cpn1-lad)/kDetPerLadderSPD[lay-1];
            if(mod>131) cpn0 = 1+(54-lad-cpn1)/kDetPerLadderSPD[lay-1];
        } break;
        case 3:{
            cpn2 = 1;
            if(lad<5) cpn0 = 5-lad;
            else cpn0 = 19-lad;
            cpn1 = 7-det;
        } break;
        case 4:{
            cpn2 = 1;
            if(lad<7) cpn0 = 7-lad;
            else cpn0 = 29-lad;
            cpn1 = 9-det;
        } break;
        case 5:{
            cpn2 = 1;
            if(lad<10) cpn0 = 10-lad;
            else cpn0 = 44-lad;
            cpn1 = 23-det;
        } break;
        case 6:{
            cpn2 = 1;
            if(lad<9) cpn0 = 9-lad;
            else cpn0 = 47-lad;
            cpn1 = 26-det;
        } break;
        default:{
            Error("RecodeDetector","New: mod=%d lay=%d not 1-6.");
            return;
        } break;
        } // end switch
        if(cpn0<1||cpn1<1||cpn2<1||
           cpn0>kITSgeoTreeCopys[lay-1][0]||
           cpn1>kITSgeoTreeCopys[lay-1][1]||
           cpn2>kITSgeoTreeCopys[lay-1][2])
            Error("RecodeDetector",
                  "cpn0=%d cpn1=%d cpn2=%d mod=%d lay=%d lad=%d det=%d",
                  cpn0,cpn1,cpn2,mod,lay,lad,det);
        return;
    } // end if
    // Old encoding
    switch (lay){
    case 1: case 2:{
        cpn2 = det;     // Detector 1-4
        cpn0 = (lad+kDetPerLadderSPD[lay-1]-1)/kDetPerLadderSPD[lay-1];
        cpn1 = (lad+kDetPerLadderSPD[lay-1]-1)%kDetPerLadderSPD[lay-1] + 1;
    } break;
    case 3: case 4: case 5 : case 6:{
        cpn2 = 1;
        cpn1 = det;
        cpn0 = lad;
    } break;
    default:{
        Error("RecodeDetector","Old: mod=%d lay=%d not 1-6.");
        return;
    } break;
    } // end switch
    if(cpn0<1||cpn1<1||cpn2<1||
       cpn0>kITSgeoTreeCopys[lay-1][0]||
       cpn1>kITSgeoTreeCopys[lay-1][1]||
       cpn2>kITSgeoTreeCopys[lay-1][2])
        Error("RecodeDetector",
              "cpn0=%d cpn1=%d cpn2=%d mod=%d lay=%d lad=%d det=%d",
              cpn0,cpn1,cpn2,mod,lay,lad,det);
    return;
}
//______________________________________________________________________
void AliITSInitGeometry::DecodeDetectorLayersvPPRasymmFMD(Int_t mod,Int_t &lay,
                                              Int_t &lad,Int_t &det){
    // decode geometry into detector module number. There are two decoding
    // Scheams. Old which does not follow the ALICE coordinate system
    // requirements, and New which dose. Note, this use of layer ladder
    // and detector numbers are strictly for internal use of this
    // specific code. They do not represent the "standard" layer ladder
    // or detector numbering except in a very old and obsoleate sence.
    // Inputs:
    //    Int_t mod      The module number assoicated with this set
    //                   of copy numbers.
    // Output:
    //    Int_t lay     The layer number
    //    Int_t lad     The ladder number
    //    Int_t det     the dettector number
    // Return:
    //    none.
  //    const Int_t kDetPerLadderSPD[2]={2,4};
    const Int_t kDetPerLadder[6]={4,4,6,8,22,25};
    const Int_t kLadPerLayer[6]={20,40,14,22,34,38};
    Int_t mod2;

    det  = 0;
    lad  = 0;
    lay  = 0;
    mod2 = 0;
    do{
        mod2 += kLadPerLayer[lay]*kDetPerLadder[lay];
        lay++;
    }while(mod2<=mod); // end while
    if(lay>6||lay<1) Error("DecodeDetectorLayers","0<lay=%d>6",lay);
    mod2 -= kLadPerLayer[lay-1]*kDetPerLadder[lay-1];
    do{
        lad++;
        mod2 += kDetPerLadder[lay-1];
    }while(mod2<=mod); // end while
    if(lad>kLadPerLayer[lay-1]||lad<1) Error("DecodeDetectorLayers",
            "lad=%d>kLadPerLayer[lay-1=%d]=%d mod=%d mod2=%d",lad,lay-1,
                                            kLadPerLayer[lay-1],mod,mod2);
    mod2 -= kDetPerLadder[lay-1];
    det = mod-mod2+1;
    if(det>kDetPerLadder[lay-1]||det<1) Error("DecodeDetectorLayers",
           "det=%d>detPerLayer[lay-1=%d]=%d mod=%d mod2=%d lad=%d",det,
                                  lay-1,kDetPerLadder[lay-1],mod,mod2,lad);
    return;
}
//______________________________________________________________________
void AliITSInitGeometry::DecodeDetectorv11Hybrid(Int_t &mod,Int_t layer,Int_t cpn0,
                                        Int_t cpn1,Int_t cpn2) const {
    // decode geometry into detector module number
    // Inputs:
    //    Int_t layer    The ITS layer
    //    Int_t cpn0     The lowest copy number
    //    Int_t cpn1     The middle copy number
    //    Int_t cpn2     the highest copy number
    // Output:
    //    Int_t &mod     The module number assoicated with this set
    //                   of copy numbers.
    // Return:
    //    none.
  const Int_t kDetPerLadderSPD[2]={2,4};
  const Int_t kDetPerLadder[6]={4,4,6,8,22,25};
  const Int_t kLadPerLayer[6]={20,40,14,22,34,38};
  Int_t lad=-1,det=-1;
  
  switch(layer) {
  case 1: case 2:{
    if (SPDIsTGeoNative()) {
      lad = cpn1+kDetPerLadderSPD[layer-1]*(cpn0-1)+1;
      det = cpn2 + 1;
    } else {
      lad = cpn1+kDetPerLadderSPD[layer-1]*(cpn0-1);
      det = cpn2;
    }
  } break;
  case 3: case 4:{
    if (SDDIsTGeoNative()) {
      lad = cpn0+1;
      det = cpn1+1;
    } else {
      lad = cpn0;
      det = cpn1;
    }
  } break;
  case 5: case 6:{
    if (SSDIsTGeoNative()) {
      lad = cpn0+1;
      det = cpn1+1;
    } else {
      lad = cpn0;
      det = cpn1;
    }
  } break;
  default:{
  } break;
  } // end switch
  mod = 0;
  for(Int_t i=0;i<layer-1;i++) mod += kLadPerLayer[i]*kDetPerLadder[i];
  mod += kDetPerLadder[layer-1]*(lad-1)+det-1;// module start at zero.
  return;
}

/*
//______________________________________________________________________
void AliITSInitGeometry::RecodeDetectorv11Hybrid(Int_t mod,Int_t &cpn0,
                                        Int_t &cpn1,Int_t &cpn2) {
    // decode geometry into detector module number. There are two decoding
    // Scheams. Old which does not follow the ALICE coordinate system
    // requirements, and New which dose.
    // Inputs:
    //    Int_t mod      The module number assoicated with this set
    //                   of copy numbers.
    // Output:
    //    Int_t cpn0     The lowest copy number
    //    Int_t cpn1     The middle copy number
    //    Int_t cpn2     the highest copy number
    // Return:
    //    none.
    const Int_t kITSgeoTreeCopys[6][3]= {{10, 2, 4},// lay=1
                                         {10, 4, 4},// lay=2
                                         {14, 6, 1},// lay=3
                                         {22, 8, 1},// lay=4
                                         {34,22, 1},// lay=5
                                         {38,25, 1}};//lay=6
    const Int_t kDetPerLadderSPD[2]={2,4};
    Int_t lay,lad,det;

    cpn0 = cpn1 = cpn2 = 0;
    DecodeDetectorLayersv11Hybrid(mod,lay,lad,det);
    // Old encoding
    switch (lay){
    case 1: case 2:{
        cpn2 = det;     // Detector 1-4
        cpn0 = (lad+kDetPerLadderSPD[lay-1]-1)/kDetPerLadderSPD[lay-1];
        cpn1 = (lad+kDetPerLadderSPD[lay-1]-1)%kDetPerLadderSPD[lay-1] + 1;
    } break;
    case 3: case 4: case 5 : case 6:{
        cpn2 = 1;
        cpn1 = det;
        cpn0 = lad;
    } break;
    default:{
        Error("RecodeDetector","Old: mod=%d lay=%d not 1-6.");
        return;
    } break;
    } // end switch
    if(cpn0<1||cpn1<1||cpn2<1||
       cpn0>kITSgeoTreeCopys[lay-1][0]||
       cpn1>kITSgeoTreeCopys[lay-1][1]||
       cpn2>kITSgeoTreeCopys[lay-1][2])
        Error("RecodeDetector",
              "cpn0=%d cpn1=%d cpn2=%d mod=%d lay=%d lad=%d det=%d",
              cpn0,cpn1,cpn2,mod,lay,lad,det);
    return;
}
*/


//______________________________________________________________________
void AliITSInitGeometry::RecodeDetectorv11Hybrid(Int_t mod,Int_t &cpn0,
                                        Int_t &cpn1,Int_t &cpn2) {
    // decode geometry into detector module number. There are two decoding
    // Scheams. Old which does not follow the ALICE coordinate system
    // requirements, and New which does.
    // Inputs:
    //    Int_t mod      The module number assoicated with this set
    //                   of copy numbers.
    // Output:
    //    Int_t cpn0     The lowest copy number  (SPD sector or SDD/SSD ladder)
    //    Int_t cpn1     The middle copy number  (SPD stave or SDD/SSD module)
    //    Int_t cpn2     the highest copy number (SPD ladder or 1 for SDD/SSD)
    // Return:
    //    none.

  const Int_t kDetPerLadderSPD[2]={2,4};
  Int_t lay,lad,det;
  DecodeDetectorLayersv11Hybrid(mod,lay,lad,det);

  if (lay<3) { // SPD
    cpn2 = det;     // Detector 1-4
    cpn0 = (lad+kDetPerLadderSPD[lay-1]-1)/kDetPerLadderSPD[lay-1];
    cpn1 = (lad+kDetPerLadderSPD[lay-1]-1)%kDetPerLadderSPD[lay-1] + 1;
    if (SPDIsTGeoNative()) {
      cpn2--;
      cpn1--;
    }
  } else { // SDD and SSD
    cpn2 = 1;
    cpn1 = det;
    cpn0 = lad;
    if (lay<5) { // SDD
      if (SDDIsTGeoNative()) {
	cpn1--;
	cpn0--;
      }
    } else { //SSD
      if (SSDIsTGeoNative()) {
	cpn1--;
	cpn0--;
      }
    }
  }
}




// //______________________________________________________________________
// void AliITSInitGeometry::DecodeDetectorLayersv11Hybrid(Int_t mod,Int_t &lay,
//                                               Int_t &lad,Int_t &det) {

//     // decode module number into detector indices for v11Hybrid
//     // Inputs:
//     //    Int_t mod      The module number associated with this set
//     //                   of copy numbers.
//     // Output:
//     //    Int_t lay     The layer number
//     //    Int_t lad     The ladder number
//     //    Int_t det     the dettector number
//     // Return:
//     //    none.

//     const Int_t kDetPerLadder[6]={4,4,6,8,22,25};
//     const Int_t kLadPerLayer[6]={20,40,14,22,34,38};
//     Int_t mod2 = 0;
//     det  = 0;
//     lad  = 0;
//     lay  = 0;

//     do{
//         mod2 += kLadPerLayer[lay]*kDetPerLadder[lay];
//         lay++;
//     } while(mod2<=mod); // end while
//     if(lay>6||lay<1) Error("DecodeDetectorLayers","0<lay=%d>6",lay);
//     mod2 -= kLadPerLayer[lay-1]*kDetPerLadder[lay-1];
//     do{
//         lad++;
//         mod2 += kDetPerLadder[lay-1];
//     } while(mod2<=mod); // end while
//     if(lad>kLadPerLayer[lay-1]||lad<1) Error("DecodeDetectorLayers",
//             "lad=%d>kLadPerLayer[lay-1=%d]=%d mod=%d mod2=%d",lad,lay-1,
//                                             kLadPerLayer[lay-1],mod,mod2);
//     mod2 -= kDetPerLadder[lay-1];
//     det = mod-mod2+1;
//     if(det>kDetPerLadder[lay-1]||det<1) Error("DecodeDetectorLayers",
//            "det=%d>detPerLayer[lay-1=%d]=%d mod=%d mod2=%d lad=%d",det,
//                                   lay-1,kDetPerLadder[lay-1],mod,mod2,lad);
//     return;
// }

//______________________________________________________________________
void AliITSInitGeometry::DecodeDetectorLayersv11Hybrid(Int_t mod,Int_t &lay,
                                              Int_t &lad,Int_t &det) {

  // decode module number into detector indices for v11Hybrid
  // mod starts from 0
  // lay, lad, det start from 1

  // Inputs:
  //    Int_t mod      The module number associated with this set
  //                   of copy numbers.
  // Output:
  //    Int_t lay     The layer number
  //    Int_t lad     The ladder number
  //    Int_t det     the dettector number

  const Int_t kDetPerLadder[6] = {4,4,6,8,22,25};
  const Int_t kLadPerLayer[6]  = {20,40,14,22,34,38};
  
  Int_t mod2 = 0;
  lay  = 0;
  
  do {
    mod2 += kLadPerLayer[lay]*kDetPerLadder[lay];
    lay++;
  } while(mod2<=mod); // end while
  if(lay>6) Error("DecodeDetectorLayers","lay=%d>6",lay);

  mod2 = kLadPerLayer[lay-1]*kDetPerLadder[lay-1] - mod2+mod;
  lad = mod2/kDetPerLadder[lay-1];

  if(lad>=kLadPerLayer[lay-1]||lad<0) Error("DecodeDetectorLayers",
					    "lad=%d not in the correct range",lad);
  det = (mod2 - lad*kDetPerLadder[lay-1])+1;
  if(det>kDetPerLadder[lay-1]||det<1) Error("DecodeDetectorLayers",
					    "det=%d not in the correct range",det);
  lad++;
}

//______________________________________________________________________
Bool_t AliITSInitGeometry::WriteVersionString(Char_t *str,Int_t length,
                        AliITSVersion_t maj,Int_t min,
                        const Char_t *cvsDate,const Char_t *cvsRevision)const{
    // fills the string str with the major and minor version number
    // Inputs:
    //   Char_t *str          The character string to hold the major 
    //                        and minor version numbers in
    //   Int_t  length        The maximum number of characters which 
    //                        can be accomidated by this string. 
    //                        str[length-1] must exist and will be set to zero
    //   AliITSVersion_t maj  The major number
    //   Int_t           min  The minor number
    //   Char_t *cvsDate      The date string from cvs
    //   Char_t *cvsRevision  The Revision string from cvs
    // Outputs:
    //   Char_t *str          The character string holding the major and minor
    //                        version numbers. str[length-1] must exist
    //                        and will be set to zero
    // Return:
    //   kTRUE if no errors
    Int_t i,n,cvsDateLength,cvsRevisionLength;

    cvsDateLength = (Int_t)strlen(cvsDate);
    cvsRevisionLength = (Int_t)strlen(cvsRevision);
    i = (Int_t)maj;
    n = 50+(Int_t)(TMath::Log10(TMath::Abs((Double_t)i)))+1+
        (Int_t)(TMath::Log10(TMath::Abs((Double_t)min)))+1
        +cvsDateLength-6+cvsRevisionLength-10;
    if(GetDebug()>1) printf("AliITSInitGeometry::WriteVersionString:"
                        "length=%d major=%d minor=%d cvsDate=%s[%d] "
                        "cvsRevision=%s[%d] n=%d\n",length,i,min,cvsDate,
                        cvsDateLength,cvsRevision,cvsRevisionLength,n);
    if(i<0) n++;
    if(min<0) n++;
    if(length<n){// not enough space to write in output string.
        Warning("WriteVersionString","Output string not long enough "
                "lenght=%d must be at least %d long\n",length,n);
        return kFALSE;
    } // end if length<n
    char *cvsrevision = new char[cvsRevisionLength-10];
    char *cvsdate = new char[cvsDateLength-6];
    for(i=0;i<cvsRevisionLength-10;i++)
        if(10+i<cvsRevisionLength-1)
            cvsrevision[i] = cvsRevision[10+i]; else cvsrevision[i] = 0;
    for(i=0;i<cvsDateLength-6;i++) if(6+i<cvsDateLength-1)
        cvsdate[i] = cvsDate[6+i]; else cvsdate[i] = 0;
    for(i=0;i<length;i++) str[i] = 0; // zero it out for now.
    i = (Int_t)maj;
    sprintf(str,"Major Version= %d Minor Version= %d Revision: %s Date: %s",
            i,min,cvsrevision,cvsdate);
    if(GetDebug()>1)printf("AliITSInitGeometry::WriteVersionString: "
                       "n=%d str=%s revision[%zu] date[%zu]\n",
                       n,str,strlen(cvsrevision),strlen(cvsdate));
    delete[] cvsrevision;
    delete[] cvsdate;
    return kTRUE;
}
//______________________________________________________________________
Bool_t AliITSInitGeometry::ReadVersionString(const Char_t *str,Int_t length,
                                             AliITSVersion_t &maj,Int_t &min,
                                             TDatime &dt)const{
    // fills the string str with the major and minor version number
    // Inputs:
    //   Char_t *str   The character string to holding the major and minor
    //                 version numbers in
    //   Int_t  length The maximum number of characters which can be
    //                 accomidated by this string. str[length-1] must exist
    // Outputs:
    //   Char_t *str   The character string holding the major and minor
    //                 version numbers unchanged. str[length-1] must exist.
    //   AliITSVersion_t maj  The major number
    //   Int_t           min  The minor number
    //   TDatime         dt   The date and time of the cvs commit
    // Return:
    //   kTRUE if no errors
    Bool_t ok;
    Char_t cvsRevision[10],cvsDate[11],cvsTime[9];
    Int_t i,m,n=strlen(str),year,month,day,hours,minuits,seconds;

    if(GetDebug()>1)printf("AliITSInitGeometry::ReadVersionString:"
                       "str=%s length=%d\n",
                       str,length);
    if(n<35) return kFALSE; // not enough space for numbers
    m = sscanf(str,"Major Version= %d  Minor Version= %d Revision: %s "
               "Date: %s %s",&i,&min,cvsRevision,cvsDate,cvsTime);
    ok = m==5;
    if(!ok) return !ok;
    m = sscanf(cvsDate,"%d/%d/%d",&year,&month,&day);
    ok = m==3;
    if(!ok) return !ok;
    m = sscanf(cvsTime,"%d:%d:%d",&hours,&minuits,&seconds);
    ok = m==3;
    if(!ok) return !ok;
    dt.Set(year,month,day,hours,minuits,seconds);
    if(GetDebug()>1)printf("AliITSInitGeometry::ReadVersionString: i=%d min=%d "
                       "cvsRevision=%s cvsDate=%s cvsTime=%s m=%d\n",
                       i,min,cvsRevision,cvsDate,cvsTime,m);
    if(GetDebug()>1)printf("AliITSInitGeometry::ReadVersionString: year=%d"
                       " month=%d day=%d hours=%d minuits=%d seconds=%d\n",
                       year,month,day,hours,minuits,seconds);
    switch (i){
    case kvITS04:{
        maj = kvITS04;
    } break;
    case kvSPD02:{
        maj = kvSPD02;
    } break;
    case kvSDD03:{
        maj = kvSDD03;
    } break;
    case kvSSD03:{
        maj = kvSSD03;
    } break;
    case kvPPRasymmFMD:{
        maj = kvPPRasymmFMD;
    } break;
    case kv11:{
        maj = kv11;
    } break;
    case kv11Hybrid:{
        maj = kv11Hybrid;
    } break;
    default:{
        maj = kvDefault;
    } break;
    } // end switch
    return ok;
}
