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
#include <TGeoMatrix.h>
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
#include "AliITSsegmentationSPD.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSsegmentationSSD.h"
#include "AliITSInitGeometry.h"
#include <TDatime.h>

ClassImp(AliITSInitGeometry)

//______________________________________________________________________
AliITSInitGeometry::AliITSInitGeometry():
TObject(),                   // Base Class
fName(0),                    // Geometry name
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
AliITSInitGeometry::AliITSInitGeometry(AliITSVersion_t version):
TObject(),                   // Base Class
fName(0),                    // Geometry name
fMajorVersion(version),      // Major version number
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

  switch (version) {
    case kv11:
        fName="AliITSv11";
	break;
    case kvDefault:
    default:
        AliFatal(Form("Undefined geometry: fMajorVersion=%d, ",(Int_t)fMajorVersion));
        fName = "Undefined";
	break;
    } // switch
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
  TDatime datetime;
  TGeoVolume *itsV = gGeoManager->GetVolume("ITSV");
  if(!itsV){
    AliError("Can't find ITS volume ITSV, exiting - nothing done!");
    return 0;
  }// end if
  const Char_t *title = itsV->GetTitle();
  if(!ReadVersionString(title,version))
    Warning("UpdateInternalGeometry","Can't read title=%s\n",title);
  SetTiming(kFALSE);
  SetSegGeom(kFALSE);
  SetDecoding(kFALSE);
  AliITSgeom *geom = CreateAliITSgeom(version);
  AliDebug(1,"AliITSgeom object has been initialized from TGeo\n");
  return geom;
}
//______________________________________________________________________
AliITSgeom* AliITSInitGeometry::CreateAliITSgeom(Int_t major){
    // Creates and Initilizes the geometry transformation class AliITSgeom
    // to values appropreate to this specific geometry. Now that
    // the segmentation is part of AliITSgeom, the detector
    // segmentations are also defined here.
    // Inputs:
    //   Int_t major   major version, see AliITSVersion_t
    //   
    // Outputs:
    //   none.
    // Return:
    //   A pointer to a new properly inilized AliITSgeom class. If
    //   pointer = 0 then failed to init.

    switch(major){
    case kv11:
        SetGeometryName("AliITSv11");
        SetVersion(kv11);
        break;
    case kvDefault:
    default:
        SetGeometryName("Undefined");
        SetVersion(kvDefault);
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

  const Int_t kItype  = 0; // Type of transformation defined 0=> Geant
  const Int_t klayers = 6; // number of layers in the ITS
  const Int_t kladders[klayers]   = {20,40,14,22,34,38}; // Number of ladders
  const Int_t kdetectors[klayers] = {4,4,6,8,22,25};// number of detector/lad
  const AliITSDetector kIdet[6]   = {kSPD,kSPD,kSDD,kSDD,kSSD,kSSD};
  const TString kPathbase = "/ALIC_1/ITSV_1/";

  const char *pathSPDsens1, *pathSPDsens2;
  pathSPDsens1="%sITSSPD_1/ITSSPDCarbonFiberSectorV_%d/ITSSPDSensitiveVirtualvolumeM0_1/ITSSPDlay1-Stave_%d/ITSSPDhalf-Stave%d_1/ITSSPDlay1-Ladder_%d/ITSSPDlay1-sensor_1";
  pathSPDsens2="%sITSSPD_1/ITSSPDCarbonFiberSectorV_%d/ITSSPDSensitiveVirtualvolumeM0_1/ITSSPDlay2-Stave_%d/ITSSPDhalf-Stave%d_1/ITSSPDlay2-Ladder_%d/ITSSPDlay2-sensor_1";

  const char *pathSDDsens1, *pathSDDsens2;
  pathSDDsens1 = "%sITSsddLayer3_1/ITSsddLadd_%d/ITSsddSensor3_%d/ITSsddWafer3_%d/ITSsddSensitivL3_1";
  pathSDDsens2 = "%sITSsddLayer4_1/ITSsddLadd_%d/ITSsddSensor4_%d/ITSsddWafer4_%d/ITSsddSensitivL4_1";

  const char *pathSSDsens1, *pathSSDsens2;
  pathSSDsens1 = "%sITSssdLayer5_1/ITSssdLay5Ladd_%d/ITSssdSensor5_%d/ITSssdSensitivL5_1";
  pathSSDsens2 = "%sITSssdLayer6_1/ITSssdLay6Ladd_%d/ITSssdSensor6_%d/ITSssdSensitivL6_1";

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
    RecodeDetector(mod,cpn0,cpn1,cpn2);

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
  } // end for module

  if(fTiming){
    time->Stop();
    time->Print();
    delete time;
  } // end if
  return kTRUE;
}

//_______________________________________________________________________
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
      Error("GetTransformation","Error in cd-ing to %s",volumePath.Data());
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
    case kvDefault:{
        Error("DecodeDetector","Major version = kvDefault, not supported");
    }break;
    case kv11:{
        return DecodeDetectorv11(mod,layer,cpn0,cpn1,cpn2);
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
    case kvDefault:{
        Error("RecodeDetector","Major version = kvDefault, not supported");
        return;
    }
    case kv11:{
        return RecodeDetectorv11(mod,cpn0,cpn1,cpn2);
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
    case kvDefault:{
        Error("DecodeDetectorLayers",
              "Major version = kvDefault, not supported");
        return;
    }break;
    case kv11:{
        return DecodeDetectorLayersv11(mod,layer,lad,det);
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
void AliITSInitGeometry::DecodeDetectorv11(Int_t &mod,Int_t layer,
                                 Int_t cpn0,Int_t cpn1,Int_t cpn2) const {
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
    lad = cpn1+kDetPerLadderSPD[layer-1]*(cpn0-1);
    det = cpn2;
  } break;
  case 3: case 4:{
    lad = cpn0+1;
    det = cpn1+1;
  } break;
  case 5: case 6:{
    lad = cpn0+1;
    det = cpn1+1;
  } break;
  default:{
  } break;
  } // end switch
  mod = 0;
  for(Int_t i=0;i<layer-1;i++) mod += kLadPerLayer[i]*kDetPerLadder[i];
  mod += kDetPerLadder[layer-1]*(lad-1)+det-1;// module start at zero.
  return;
}


//______________________________________________________________________
void AliITSInitGeometry::RecodeDetectorv11(Int_t mod,Int_t &cpn0,
					   Int_t &cpn1,Int_t &cpn2) {
    // decode geometry into detector module number using the new decoding
    // Scheme.
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

    DecodeDetectorLayersv11(mod,lay,lad,det);
    if (lay<3) { // SPD
        cpn2 = det;     // Detector 1-4
        cpn0 = (lad+kDetPerLadderSPD[lay-1]-1)/kDetPerLadderSPD[lay-1];
        cpn1 = (lad+kDetPerLadderSPD[lay-1]-1)%kDetPerLadderSPD[lay-1] + 1;
    } else { // SDD and SSD
        cpn2 = 1;
        cpn1 = det;
        cpn0 = lad;
        if (lay<5) { // SDD
	  cpn1--;
	  cpn0--;
        } else { //SSD
	  cpn1--;
	  cpn0--;
        } // end if Lay<5/else
    } // end if lay<3/else

}


//______________________________________________________________________
void AliITSInitGeometry::DecodeDetectorLayersv11(Int_t mod,Int_t &lay,
						 Int_t &lad,Int_t &det) {

  // decode module number into detector indices for v11
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
Bool_t AliITSInitGeometry::WriteVersionString(Char_t *str,Int_t length,AliITSVersion_t maj)const{
    // fills the string str with the major version number
    // Inputs:
    //   Char_t *str          The character string to hold the major version number
    //   Int_t  length        The maximum number of characters which 
    //                        can be accommodated by this string. 
    //                        str[length-1] must exist
    //   AliITSVersion_t maj  The major number


    Int_t i = (Int_t)maj;
 
    snprintf(str,length-1,"Major Version= %d",i);
    return kTRUE;
}
//______________________________________________________________________
Bool_t AliITSInitGeometry::ReadVersionString(const Char_t *str,AliITSVersion_t &maj)const{
    // fills the string str with the major and minor version number
    // Inputs:
    //   Char_t *str   The character string to holding the major version number
    //   Int_t  length The maximum number of characters which can be
    //                 accommodated by this string. str[length-1] must exist
    // Outputs:
    //   AliITSVersion_t maj  The major number

    // Return:
    //   kTRUE if no errors

  Bool_t retcode=kFALSE;
  Int_t n=strlen(str);
  if(n<15) return retcode; // not enough space for numbers
  Int_t m,i;
  m = sscanf(str,"Major Version= %2d",&i);
  maj = kvDefault;
  if(m>0){
    retcode = kTRUE;
    if(i==11){
      maj = kv11;
    }
  }
  return retcode;
}


