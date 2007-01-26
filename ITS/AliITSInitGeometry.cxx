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

#include "AliLog.h"
#include "AliITSgeomSPD.h"
#include "AliITSgeomSDD.h"
#include "AliITSgeomSSD.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSsegmentationSSD.h"
#include "AliITSgeom.h"
#include "AliITSInitGeometry.h"

ClassImp(AliITSInitGeometry)
//______________________________________________________________________
AliITSInitGeometry::AliITSInitGeometry():
TObject(),
fName(),
fMinorVersion(0),
fMajorVersion(0),
fTiming(kFALSE),
fSegGeom(kFALSE),
fDecode(kFALSE){
    // Default Creator
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   A default inilized AliITSInitGeometry object
}
//______________________________________________________________________
AliITSInitGeometry::AliITSInitGeometry(const Char_t *name,Int_t minorversion):
TObject(),
fName(name),
fMinorVersion(minorversion),
fMajorVersion(0),
fTiming(kFALSE),
fSegGeom(kFALSE),
fDecode(kFALSE){
    // Default Creator
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   A default inilized AliITSInitGeometry object

    if(fName.CompareTo("AliITSvPPRasymmFMD")==0)if(fMinorVersion==1||
						   fMinorVersion==2){
	fMajorVersion=10;
	return;
    } // end if
    // if not defined geometry error
    Error("AliITSInitGeometry(name,version)"," Name must be AliITSvPPRasymmFMD"
	" and version must be 1 or 2 for now.");
    fMinorVersion = 0;
    fName = "";
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
    AliFatal("The geometry manager has not been initialized (e.g. TGeoManager::Import(\"geometry.root\")should be called in advance) - exit forced");
    return kFALSE;
  }
  switch(fMajorVersion){
  case 10:{ // only case defined so far
    return InitAliITSgeomPPRasymmFMD(geom);
  }break; // end case
  default:{
    Error("InitAliITSgeom","Undefined geomtery");
    return kFALSE;
  } break; // end case
  } // end switch
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
	{"%sIT12_1/I12B_%d/I10B_%d/I107_%d/I101_1/ITS1_1", // lay=1
	 "%sIT12_1/I12B_%d/I20B_%d/I1D7_%d/I1D1_1/ITS2_1", // lay=2
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
    Int_t mod,nmods=0,lay,lad,det,cpn0,cpn1,cpn2;
    Double_t tran[3]={0.0,0.0,0.0},rot[10]={9*0.0,1.0};
    TArrayD shapePar;
    TString path,shapeName;
    TGeoHMatrix materix;
    Bool_t initSeg[3]={kFALSE,kFALSE,kFALSE};
    TStopwatch *time = 0x0;if(fTiming) time=new TStopwatch();

    if(fTiming) time->Start();
    for(mod=0;mod<klayers;mod++) nmods += kladders[mod]*kdetectors[mod];
    geom->Init(kItype,klayers,kladders,kdetectors,nmods);
    for(mod=0;mod<nmods;mod++){
        DecodeDetectorLayers(mod,lay,lad,det); // Write
        geom->CreateMatrix(mod,lay,lad,det,kIdet[lay-1],tran,rot);
        RecodeDetector(mod,cpn0,cpn1,cpn2); // Write reusing lay,lad,det.
        path.Form(kNames[fMinorVersion-1][lay-1].Data(),
                  kPathbase.Data(),cpn0,cpn1,cpn2);
        geom->GetGeomMatrix(mod)->SetPath(path);
        GetTransformation(path.Data(),materix);
        geom->SetTrans(mod,materix.GetTranslation());
        geom->SetRotMatrix(mod,materix.GetRotationMatrix());
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
    }
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
    }
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
    }
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
    }
    if (classType==TGeoTube::Class()) {
	shapeType = "TUBE";
	npar = 3;
	par.Set(npar);
	TGeoTube *tube = (TGeoTube*)shape;
	par.AddAt(tube->GetRmin(),0);
	par.AddAt(tube->GetRmax(),1);
	par.AddAt(tube->GetDz(),2);
	return kTRUE;
    }
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
    }
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
    }
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
    }
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
    }
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
    }
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
    }
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
    }
    if (classType==TGeoEltu::Class()) {
	shapeType = "ELTU";
	npar = 3;
	par.Set(npar);
	TGeoEltu *eltu = (TGeoEltu*)shape;
	par.AddAt(eltu->GetA(),0);
	par.AddAt(eltu->GetB(),1);
	par.AddAt(eltu->GetDz(),2);
	return kTRUE;
    }
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
    }
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
    }
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
    }
    Error("GetShape","Getting shape parameters for shape %s not implemented",
	  shape->ClassName());
    return kFALSE;
}
//______________________________________________________________________
void AliITSInitGeometry::DecodeDetector(Int_t &mod,Int_t layer,Int_t cpn0,
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
void AliITSInitGeometry::DecodeDetectorLayers(Int_t mod,Int_t &lay,
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
    if(lad>kLadPerLayer[lay-1]||lad<1) Error("DecodeDetectorLayera",
            "lad=%d>kLadPerLayer[lay-1=%d]=%d mod=%d mod2=%d",lad,lay-1,
                                            kLadPerLayer[lay-1],mod,mod2);
    mod2 -= kDetPerLadder[lay-1];
    det = mod-mod2+1;
    if(det>kDetPerLadder[lay-1]||det<1) Error("DecodeDetectorLayers",
           "det=%d>detPerLayer[lay-1=%d]=%d mod=%d mod2=%d lad=%d",det,
                                  lay-1,kDetPerLadder[lay-1],mod,mod2,lad);
    return;
}

