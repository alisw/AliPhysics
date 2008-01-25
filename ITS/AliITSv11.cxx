/**************************************************************************
 * Copyright(c) 2007-2008, ALICE Experiment at CERN, All rights reserved. *
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


//************************************************************************
//
//                 Inner Traking System geometry v11
//
//  Based on ROOT geometrical modeler
//
// B. Nilsen, L. Gaudichet
//************************************************************************


#include <TClonesArray.h>
#include <TLorentzVector.h>

#include "AliITS.h"
#include "AliITSDetTypeSim.h"
#include <TVirtualMC.h>

#include "AliITSgeom.h"
#include "AliITSgeomSDD.h"
#include "AliITSgeomSPD.h"
#include "AliITSgeomSSD.h"
#include "AliITShit.h"

#include "AliITSCalibrationSDD.h"

#include "AliITSsegmentationSDD.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSsegmentationSSD.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliTrackReference.h"
#include "AliMC.h"

#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoPcon.h>
#include "AliITSv11.h"
#include "AliITSv11GeometrySPD.h"
#include "AliITSv11GeometrySDD.h"
#include "AliITSv11GeometrySSD.h"
#include "AliITSv11GeometrySupport.h"



ClassImp(AliITSv11)
 


//______________________________________________________________________
AliITSv11::AliITSv11() : 
AliITS(),
fGeomDetOut(kFALSE),
fGeomDetIn(kFALSE),
fByThick(kTRUE),
fMajorVersion(IsVersion()),
fMinorVersion(0),
fEuclidGeomDet(),
fRead(),
fWrite(),
fSPDgeom(),
fSDDgeom(0),
fSSDgeom(),
fSupgeom(),
fIgm(kv11)
{
  //    Standard default constructor for the ITS version 11.

    fIdN          = 0;
    fIdName       = 0;
    fIdSens       = 0;
    Int_t i;
    for(i=0;i<60;i++) fRead[i] = '\0';
    for(i=0;i<60;i++) fWrite[i] = '\0';
    for(i=0;i<60;i++) fEuclidGeomDet[i] = '\0';
    strncpy(fRead,"$ALICE_ROOT/ITS/ITSgeometry_vPPRasymmFMD.det",60);
}


//______________________________________________________________________
AliITSv11::AliITSv11(const char *name, const char *title): 
AliITS("ITS", title),
fGeomDetOut(kFALSE),
fGeomDetIn(kFALSE),
fByThick(kTRUE),
fMajorVersion(IsVersion()),
fMinorVersion(0),
fEuclidGeomDet(),
fRead(),
fWrite(),
fSPDgeom(),
fSDDgeom(0),
fSSDgeom(),
fSupgeom(),
fIgm(kv11)
{
  //    Standard constructor for the ITS version 11.

  fSDDgeom = new AliITSv11GeometrySDD(0);

  Int_t i;
  fIdN = 6;
  fIdName = new TString[fIdN];
  fIdName[0] = name; // removes warning message
  fIdName[0] = "ITS1";
  fIdName[1] = "ITS2";
  fIdName[2] = fSDDgeom->GetSenstiveVolumeName3();
  fIdName[3] = fSDDgeom->GetSenstiveVolumeName4();
  fIdName[4] = "ITS5";
  fIdName[5] = "ITS6";
  fIdSens    = new Int_t[fIdN];
  for(i=0;i<fIdN;i++) fIdSens[i] = 0;
  // not needed, fByThick set to kTRUE in in the member initialization lis
  

  fEuclidGeometry="$ALICE_ROOT/ITS/ITSgeometry_vPPRasymm2.euc";
  strncpy(fEuclidGeomDet,"$ALICE_ROOT/ITS/ITSgeometry_vPPRasymm2.det",60);
  strncpy(fRead,fEuclidGeomDet,60);
  strncpy(fWrite,fEuclidGeomDet,60);
  strncpy(fRead,"$ALICE_ROOT/ITS/ITSgeometry_vPPRasymmFMD.det",60);
}
//______________________________________________________________________
AliITSv11::AliITSv11(Int_t debugITS,Int_t debugSPD,Int_t debugSDD,
		   Int_t debugSSD,Int_t debugSUP) :
AliITS("ITS","ITS geometry v11"),
fGeomDetOut(kFALSE),
fGeomDetIn(kFALSE),
fByThick(kTRUE),
fMajorVersion(IsVersion()),
fMinorVersion(0),
fEuclidGeomDet(),
fRead(),
fWrite(),
fSPDgeom(),
fSDDgeom(0),
fSSDgeom(),
fSupgeom(),
fIgm(kv11)
{
  // Standard default constructor for the ITS version 11.


  //   fSPDgeom = new AliITSv11GeometrySPD(debugSPD);
  fSDDgeom = new AliITSv11GeometrySDD(debugSDD);
  fSDDgeom->SetDebug(debugSDD);
  //   fSupgeom = new AliITSv11GeometrySupport(debugSUP);

  Int_t i;
  fIdN = 6;
  fIdName = new TString[fIdN];
  fIdName[0] = fSPDgeom->GetSenstiveVolumeName1();
  fIdName[1] = fSPDgeom->GetSenstiveVolumeName2();
  fIdName[2] = fSDDgeom->GetSenstiveVolumeName3();
  fIdName[3] = fSDDgeom->GetSenstiveVolumeName4();
  fIdName[4] = fSSDgeom->GetSenstiveVolumeName5();
  fIdName[5] = fSSDgeom->GetSenstiveVolumeName6();
  fIdSens    = new Int_t[fIdN];
  for(i=0;i<fIdN;i++) fIdSens[i] = 0;
  fEuclidOut    = kFALSE; // Don't write Euclide file
  
  fEuclidGeometry="$ALICE_ROOT/ITS/ITSgeometry_vPPRasymm2.euc";
  strncpy(fEuclidGeomDet,"$ALICE_ROOT/ITS/ITSgeometry_vPPRasymm2.det",60);
  strncpy(fRead,fEuclidGeomDet,60);
  strncpy(fWrite,fEuclidGeomDet,60);
  strncpy(fRead,"$ALICE_ROOT/ITS/ITSgeometry_vPPRasymmFMD.det",60);

  debugITS = (debugSPD && debugSSD && debugSUP && debugSDD); //remove temp. warnings
}
//______________________________________________________________________
AliITSv11::~AliITSv11() {
  delete fSDDgeom;
}
//______________________________________________________________________
void AliITSv11::BuildGeometry(){

}
//______________________________________________________________________
void AliITSv11::CreateGeometry(){
    //
    // Create ROOT geometry
    //
    // These constant character strings are set by cvs during commit
    // do not change them unless you know what you are doing!
    const Char_t *cvsDate="$Date$";
    const Char_t *cvsRevision="$Revision$";

    TGeoManager *geoManager = gGeoManager;
    TGeoVolume *vALIC = geoManager->GetTopVolume();

    TGeoPcon *sITS = new TGeoPcon("ITS Top Volume",0.0,360.0,2);

    // DefineSection(section number, Z, Rmin, Rmax).
    const Double_t kcm = 1.0;
    sITS->DefineSection(0,-300.0*kcm,0.01*kcm,50.0*kcm);
    sITS->DefineSection(1,+300.0*kcm,0.01*kcm,50.0*kcm);

    TGeoMedium *air = gGeoManager->GetMedium("ITS_AIR$");
    TGeoVolume *vITS = new TGeoVolume("ITSV",sITS,air);
    vITS->SetVisibility(kFALSE);
    const Int_t length=100;
    Char_t vstrng[length];
    if(fIgm.WriteVersionString(vstrng,length,(AliITSVersion_t)IsVersion(),
                               fMinorVersion,cvsDate,cvsRevision))
        vITS->SetTitle(vstrng);
    //printf("Title set to %s\n",vstrng);
    vALIC->AddNode(vITS,1,0);

//   fSPDgeom->CenteralSPD(vITS);

  fSDDgeom->Layer3(vITS);
  fSDDgeom->Layer4(vITS);

//     fSupgeom->SPDCone(vITS);
//     fSupgeom->SPDThermalSheald(vITS);
//     fSupgeom->SDDCone(vITS);
//     fSupgeom->SSDCone(vITS);
//     fSupgeom->ServicesCableSupport(vITS);

}
//______________________________________________________________________
void AliITSv11::CreateMaterials(){
    // Create Standard ITS Materials
    // Inputs:
    //  none.
    // Outputs:
    //  none.
    // Return:
    // none.

    
    //
    fSPDgeom->AliITSv11Geometry::CreateDefaultMaterials();
    // Detector specific material definistions
    fSPDgeom->CreateMaterials();
    fSDDgeom->CreateMaterials();
    fSSDgeom->CreateMaterials();
    fSupgeom->CreateMaterials();
}
/*
//______________________________________________________________________
void AliITSv11::InitAliITSgeom(){
  //
  // Fill fITSgeom with the 3 sub-detector geometries
  //

  if (gGeoManager) gGeoManager->Export("geometry.root");

    const Int_t knlayers = 6;
    const Int_t kndeep = 3;
    const AliITSDetector kidet[knlayers]={kSPD,kSPD,kSDD,kSDD,kSSD,kSSD};
    const TString knames[knlayers] = {
      "AliITSv11:spd missing",  // lay=1
      "AliITSv11:spd missing",  // lay=2
      "/ALIC_1/ITSV_1/ITSsddLayer3_1/ITSsddLadd_%d/ITSsddSensor_%d/ITSsddWafer_%d", // lay=3
      "/ALIC_1/ITSV_1/ITSsddLayer4_1/ITSsddLadd_%d/ITSsddSensor_%d/ITSsddWafer_%d", // lay=4
      "AliITSv11:ssd missing",  // lay=5
      "AliITSv11:ssd missing"   // lay=6
    };

    const Int_t kitsGeomTreeCopys[knlayers][kndeep]= {{10, 2, 4},// lay=1
                                                     {10, 4, 4}, // lay=2
                                                     {14, 6, 1}, // lay=3
                                                     {22, 8, 1}, // lay=4
                                                     {34,22, 1}, // lay=5
                                                     {38,25, 1}};// lay=6
    Int_t       nlad[knlayers],ndet[knlayers];
    Int_t       mod,lay,lad=0,det=0,i,j,k,cp0,cp1,cp2;
    TString path,shapeName;
    TGeoHMatrix materix;
    Double_t trans[3]={3*0.0},rot[10]={9*0.0,1.0};
    TArrayD shapePar;
    TArrayF shapeParF;
    Bool_t shapeDefined[3]={kFALSE,kFALSE,kFALSE};

    AliDebug(1,"Reading Geometry transformation directly from Modler.");
    mod = 0;
    for(i=0;i<knlayers;i++){
        k = 1;
        for(j=0;j<kndeep;j++) if(kitsGeomTreeCopys[i][j]!=0)
            k *= TMath::Abs(kitsGeomTreeCopys[i][j]);
        mod += k;
    } // end for i

    SetITSgeom(0);
    nlad[0]=20;nlad[1]=40;nlad[2]=14;nlad[3]=22;nlad[4]=34;nlad[5]=38;
    ndet[0]= 4;ndet[1]= 4;ndet[2]= 6;ndet[3]= 8;ndet[4]=22;ndet[5]=25;
    AliITSgeom* geom = new AliITSgeom(0,6,nlad,ndet,mod);
    SetITSgeom(geom);
    mod = 0;
    for(lay=1;lay<=knlayers;lay++){

        for(cp0=0; cp0<kitsGeomTreeCopys[lay-1][0]; cp0++){
            for(cp1=0; cp1<kitsGeomTreeCopys[lay-1][1]; cp1++){
                for(cp2=1; cp2<=kitsGeomTreeCopys[lay-1][2]; cp2++){

                    path.Form(knames[lay-1].Data(),
                              cp0,cp1,cp2);
                    switch (lay){
                    case 1:{
                        det = cp2;
                        lad = cp1+2*(cp0-1);
                    }break;
                    case 2:{
                        det = cp2;
                        lad = cp1+4*(cp0-1);
                    } break;
                    case 3: case 4: case 5: case 6:{
                        det = cp1;
                        lad = cp0;
                    } break;
                    } // end switch
                         //AliInfo(Form("path=%s lay=%d lad=%d det=%d",
                         //             path.Data(),lay,lad,det));
                    gMC->GetTransformation(path.Data(),materix);
                    gMC->GetShape(path.Data(),shapeName,shapePar);
                    shapeParF.Set(shapePar.GetSize());
                    for(i=0;i<shapePar.GetSize();i++) shapeParF[i]=shapePar[i];
                    geom->CreateMatrix(mod,lay,lad,det,kidet[lay-1],trans,rot);
                    geom->SetTrans(mod,materix.GetTranslation());
                    geom->SetRotMatrix(mod,materix.GetRotationMatrix());
		    geom->GetGeomMatrix(mod)->SetPath(path.Data());
                    switch (lay){
                    case 1: case 2:
			if(!shapeDefined[kSPD]){
                        geom->ReSetShape(kSPD,new AliITSgeomSPD425Short(
                                shapeParF.GetSize(),shapeParF.GetArray()));
			shapeDefined[kSPD] = kTRUE;
                    }break;
                    case 3: case 4:
			if(!shapeDefined[kSDD]){
                        geom->ReSetShape(kSDD,new AliITSgeomSDD256(
                                shapeParF.GetSize(),shapeParF.GetArray()));
			shapeDefined[kSDD] = kTRUE;
                    }break;
                    case 5: case 6:
			if(!shapeDefined[kSSD]){
                        geom->ReSetShape(kSSD,new AliITSgeomSSD75and275(
                                shapeParF.GetSize(),shapeParF.GetArray()));
			shapeDefined[kSSD] = kTRUE;
                    }break;
                    default:{
                    }break;
                    } // end switch
                    mod++;
                } /// end for cp2
            } // end for cp1
        } // end for cp0
    } // end for lay

//   fSDDgeom->ExportSensorGeometry(GetITSgeom(), +3, 0);  //SDD
}
*/
//______________________________________________________________________
void AliITSv11::Init(){
  //
  //     Initialise the ITS after it has been created.
  //

  //AliInfo(Form("Minor version %d",fMinorVersion));
    //
    UpdateInternalGeometry();
    AliITS::Init();
    if(fGeomDetOut) GetITSgeom()->WriteNewFile(fWrite);

    //
/*
    if(fRead[0]=='\0') strncpy(fRead,fEuclidGeomDet,60);
    if(fWrite[0]=='\0') strncpy(fWrite,fEuclidGeomDet,60);
    if(GetITSgeom()!=0) SetITSgeom(0x0);
    AliITSgeom* geom = new AliITSgeom();
    SetITSgeom(geom);
    if(fGeomDetIn) GetITSgeom()->ReadNewFile(fRead);
    else this->InitAliITSgeom();
    if(fGeomDetOut) GetITSgeom()->WriteNewFile(fWrite);
    AliITS::Init();
*/    //
}

// //______________________________________________________________________
// void AliITSv11::SetDefaults(){
//   //
//   // Set response ans segmentation models for SPD, SDD and SSD
//   //
//      const Float_t kconv = 1.0e+04; // convert cm to microns
//     AliInfo("Called");    

//     if(!fDetTypeSim) fDetTypeSim = new AliITSDetTypeSim();
//     fDetTypeSim->SetITSgeom(GetITSgeom());
  
//     AliITSgeomSPD  *s0;
//     AliITSgeomSDD  *s1;
//     AliITSgeomSSD  *s2;
//     Int_t i;
//     Float_t bx[256],bz[280];
   
//     fDetTypeSim->ResetCalibrationArray();
//     fDetTypeSim->ResetSegmentation();
//     fDetTypeSim->SetDefaults();
    
//     //SPD
//     s0 = (AliITSgeomSPD*) GetITSgeom()->GetShape(kSPD);// Get shape info. Do it this way for now.
//     AliITSsegmentationSPD* seg0 = (AliITSsegmentationSPD*)fDetTypeSim->GetSegmentationModel(0);
//     seg0->SetDetSize(s0->GetDx()*2.*kconv, // base this on AliITSgeomSPD
// 		     s0->GetDz()*2.*kconv, // for now.
// 		     s0->GetDy()*2.*kconv); // x,z,y full width in microns.
//     seg0->SetNPads(256,160);// Number of Bins in x and z
//     for(i=000;i<256;i++) bx[i] =  50.0; // in x all are 50 microns.
//     for(i=000;i<160;i++) bz[i] = 425.0; // most are 425 microns except below
//     for(i=160;i<280;i++) bz[i] =   0.0; // Outside of detector.
//     bz[ 31] = bz[ 32] = 625.0; // first chip boundry
//     bz[ 63] = bz[ 64] = 625.0; // first chip boundry
//     bz[ 95] = bz[ 96] = 625.0; // first chip boundry
//     bz[127] = bz[128] = 625.0; // first chip boundry
//     bz[160] = 425.0; // Set so that there is no zero pixel size for fNz.
//     seg0->SetBinSize(bx,bz); // Based on AliITSgeomSPD for now.
//     SetSegmentationModel(kSPD,seg0);
//     // set digit and raw cluster classes to be used
//     const char *kData0=(fDetTypeSim->GetCalibrationModel(GetITSgeom()->GetStartSPD()))->DataType();
//     if (strstr(kData0,"real")) fDetTypeSim->SetDigitClassName(kSPD,"AliITSdigit");
//     else fDetTypeSim->SetDigitClassName(kSPD,"AliITSdigitSPD");



//     // SDD
//     s1 = (AliITSgeomSDD*) GetITSgeom()->GetShape(kSDD);// Get shape info. Do it this way for now.
//     AliITSsegmentationSDD* seg1 = (AliITSsegmentationSDD*)fDetTypeSim->GetSegmentationModel(1);
//     seg1->SetDetSize(s1->GetDx()*kconv, // base this on AliITSgeomSDD
// 		     s1->GetDz()*2.*kconv, // for now.
// 		     s1->GetDy()*2.*kconv); // x,z,y full width in microns.
//     seg1->SetNPads(256,256);// Use AliITSgeomSDD for now
//     SetSegmentationModel(kSDD,seg1);
//     const char *kData1=(fDetTypeSim->GetCalibrationModel(GetITSgeom()->GetStartSDD()))->DataType();
//     AliITSCalibrationSDD* rsp = (AliITSCalibrationSDD*)fDetTypeSim->GetCalibrationModel(GetITSgeom()->GetStartSDD());
//     const char *kopt=rsp->GetZeroSuppOption();
//     if((!strstr(kopt,"2D")) && (!strstr(kopt,"1D")) || strstr(kData1,"real") ){
// 	fDetTypeSim->SetDigitClassName(kSDD,"AliITSdigit");
//     } else fDetTypeSim->SetDigitClassName(kSDD,"AliITSdigitSDD");


//     // SSD
//     s2 = (AliITSgeomSSD*) GetITSgeom()->GetShape(kSSD);// Get shape info. Do it this way for now.
//     AliITSsegmentationSSD* seg2 = (AliITSsegmentationSSD*)fDetTypeSim->GetSegmentationModel(2);
//     seg2->SetDetSize(s2->GetDx()*2.*kconv, // base this on AliITSgeomSSD
// 		     s2->GetDz()*2.*kconv, // for now.
// 		     s2->GetDy()*2.*kconv); // x,z,y full width in microns.
//     seg2->SetPadSize(95.,0.); // strip x pitch in microns
//     seg2->SetNPads(768,0); // number of strips on each side.
//     seg2->SetAngles(0.0075,0.0275); // strip angels rad P and N side.
//     seg2->SetAnglesLay5(0.0075,0.0275); // strip angels rad P and N side.
//     seg2->SetAnglesLay6(0.0275,0.0075); // strip angels rad P and N side.
//     SetSegmentationModel(kSSD,seg2); 
//         const char *kData2=(fDetTypeSim->GetCalibrationModel(GetITSgeom()->GetStartSSD()))->DataType();
//     if(strstr(kData2,"real") ) fDetTypeSim->SetDigitClassName(kSSD,"AliITSdigit");
//     else fDetTypeSim->SetDigitClassName(kSSD,"AliITSdigitSSD");
//     if(fgkNTYPES>3){
// 	Warning("SetDefaults",
// 		"Only the four basic detector types are initialised!");
//     }// end if

    
//     return;
// }


//______________________________________________________________________
void AliITSv11::SetDefaults(){
  //
  // Set response and segmentation models for SPD, SDD and SSD
  //
     const Float_t kconv = 1.0e+04; // convert cm to microns
    AliInfo("Called");    

//     if(!fDetTypeSim) fDetTypeSim = new AliITSDetTypeSim();
//     fDetTypeSim->SetITSgeom(GetITSgeom());
    if(!fDetTypeSim) {
      Warning("SetDefaults","Error fDetTypeSim not defined");
      return;
    }
  
    AliITSgeomSPD  *s0;
    AliITSgeomSDD  *s1;
    AliITSgeomSSD  *s2;
    Int_t i;
    Float_t bx[256],bz[280];
   
    fDetTypeSim->ResetCalibrationArray();
    fDetTypeSim->ResetSegmentation();
    fDetTypeSim->SetDefaults();
    
    //SPD
    s0 = (AliITSgeomSPD*) GetITSgeom()->GetShape(kSPD);// Get shape info. Do it this way for now.
    AliITSsegmentationSPD* seg0 = (AliITSsegmentationSPD*)fDetTypeSim->GetSegmentationModel(0);
    seg0->SetDetSize(s0->GetDx()*2.*kconv, // base this on AliITSgeomSPD
		     s0->GetDz()*2.*kconv, // for now.
		     s0->GetDy()*2.*kconv); // x,z,y full width in microns.
    seg0->SetNPads(256,160);// Number of Bins in x and z
    for(i=000;i<256;i++) bx[i] =  50.0; // in x all are 50 microns.
    for(i=000;i<160;i++) bz[i] = 425.0; // most are 425 microns except below
    for(i=160;i<280;i++) bz[i] =   0.0; // Outside of detector.
    bz[ 31] = bz[ 32] = 625.0; // first chip boundry
    bz[ 63] = bz[ 64] = 625.0; // first chip boundry
    bz[ 95] = bz[ 96] = 625.0; // first chip boundry
    bz[127] = bz[128] = 625.0; // first chip boundry
    bz[160] = 425.0; // Set so that there is no zero pixel size for fNz.
    seg0->SetBinSize(bx,bz); // Based on AliITSgeomSPD for now.
    SetSegmentationModel(kSPD,seg0);
    // set digit and raw cluster classes to be used
    const char *kData0=(fDetTypeSim->GetCalibrationModel(GetITSgeom()->GetStartSPD()))->DataType();
    if (strstr(kData0,"real")) fDetTypeSim->SetDigitClassName(kSPD,"AliITSdigit");
    else fDetTypeSim->SetDigitClassName(kSPD,"AliITSdigitSPD");



    // SDD
    s1 = (AliITSgeomSDD*) GetITSgeom()->GetShape(kSDD);// Get shape info. Do it this way for now.
    AliITSsegmentationSDD* seg1 = (AliITSsegmentationSDD*)fDetTypeSim->GetSegmentationModel(1);
    seg1->SetDetSize(s1->GetDx()*kconv, // base this on AliITSgeomSDD
		     s1->GetDz()*2.*kconv, // for now.
		     s1->GetDy()*2.*kconv); // x,z,y full width in microns.
    seg1->SetNPads(256,256);// Use AliITSgeomSDD for now
    SetSegmentationModel(kSDD,seg1);
    const char *kData1=(fDetTypeSim->GetCalibrationModel(GetITSgeom()->GetStartSDD()))->DataType();
    AliITSCalibrationSDD* rsp = (AliITSCalibrationSDD*)fDetTypeSim->GetCalibrationModel(GetITSgeom()->GetStartSDD());
    const char *kopt=rsp->GetZeroSuppOption();
    if((!strstr(kopt,"2D")) && (!strstr(kopt,"1D")) || strstr(kData1,"real") ){
	fDetTypeSim->SetDigitClassName(kSDD,"AliITSdigit");
    } else fDetTypeSim->SetDigitClassName(kSDD,"AliITSdigitSDD");


    // SSD
    s2 = (AliITSgeomSSD*) GetITSgeom()->GetShape(kSSD);// Get shape info. Do it this way for now.
    AliITSsegmentationSSD* seg2 = (AliITSsegmentationSSD*)fDetTypeSim->GetSegmentationModel(2);
    seg2->SetDetSize(s2->GetDx()*2.*kconv, // base this on AliITSgeomSSD
		     s2->GetDz()*2.*kconv, // for now.
		     s2->GetDy()*2.*kconv); // x,z,y full width in microns.
    seg2->SetPadSize(95.,0.); // strip x pitch in microns
    seg2->SetNPads(768,0); // number of strips on each side.
    seg2->SetAngles(0.0075,0.0275); // strip angels rad P and N side.
    seg2->SetAnglesLay5(0.0075,0.0275); // strip angels rad P and N side.
    seg2->SetAnglesLay6(0.0275,0.0075); // strip angels rad P and N side.
    SetSegmentationModel(kSSD,seg2); 
        const char *kData2=(fDetTypeSim->GetCalibrationModel(GetITSgeom()->GetStartSSD()))->DataType();
    if(strstr(kData2,"real") ) fDetTypeSim->SetDigitClassName(kSSD,"AliITSdigit");
    else fDetTypeSim->SetDigitClassName(kSSD,"AliITSdigitSSD");
    if(fgkNTYPES>3){
	Warning("SetDefaults",
		"Only the four basic detector types are initialised!");
    }// end if

    
    return;
}





//______________________________________________________________________
void AliITSv11::DrawModule() const{

}

// //______________________________________________________________________
// void AliITSv11::StepManager(){
//   //
//   //    Called for every step in the ITS, then calles the AliITShit class
//   // creator with the information to be recoreded about that hit.
//   //
//     Int_t         copy, id;
//     TLorentzVector position, momentum;
//     static TLorentzVector position0;
//     static Int_t stat0=0;

//     if(!(this->IsActive())){
// 	return;
//     } // end if !Active volume.

//     if(!(gMC->TrackCharge())) return;

//     id=gMC->CurrentVolID(copy);

//     Bool_t sensvol = kFALSE;
//     for(Int_t kk=0;kk<6;kk++)if(id == fIdSens[kk])sensvol=kTRUE;
//     if(sensvol && (gMC->IsTrackExiting())){
// 	copy = fTrackReferences->GetEntriesFast();
// 	TClonesArray &lTR = *fTrackReferences;
// 	// Fill TrackReference structure with this new TrackReference.
// 	new(lTR[copy]) AliTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber());
//     } // if Outer ITS mother Volume


//     Int_t   copy1,copy2;  
//     Int_t   vol[5];
//     TClonesArray &lhits = *fHits;
//     //
//     // Track status
//     vol[3] = 0;
//     vol[4] = 0;
//     if(gMC->IsTrackInside())      vol[3] +=  1;
//     if(gMC->IsTrackEntering())    vol[3] +=  2;
//     if(gMC->IsTrackExiting())     vol[3] +=  4;
//     if(gMC->IsTrackOut())         vol[3] +=  8;
//     if(gMC->IsTrackDisappeared()) vol[3] += 16;
//     if(gMC->IsTrackStop())        vol[3] += 32;
//     if(gMC->IsTrackAlive())       vol[3] += 64;
//     //
//     // Fill hit structure.
//     if(!(gMC->TrackCharge())) return;
//     //
//     // Only entering charged tracks
//     if((id = gMC->CurrentVolID(copy)) == fIdSens[0]) {
// 	vol[0] = 1;
// 	id = gMC->CurrentVolOffID(2,copy);
// 	//detector copy in the ladder = 1<->4  (ITS1 < I101 < I103 < I10A)
// 	vol[1] = copy;
// 	gMC->CurrentVolOffID(3,copy1);
// 	//ladder copy in the module   = 1<->2  (I10A < I12A)
// 	gMC->CurrentVolOffID(4,copy2);
// 	//module copy in the layer    = 1<->10 (I12A < IT12)
// 	vol[2] = copy1+(copy2-1)*2;//# of ladders in one module  = 2
//     } else if(id == fIdSens[1]){
// 	vol[0] = 2;
// 	id = gMC->CurrentVolOffID(2,copy);
// 	//detector copy in the ladder = 1<->4  (ITS2 < I1D1 < I1D3 < I20A)
// 	vol[1] = copy;
// 	gMC->CurrentVolOffID(3,copy1);
// 	//ladder copy in the module   = 1<->4  (I20A < I12A)
// 	gMC->CurrentVolOffID(4,copy2);
// 	//module copy in the layer    = 1<->10 (I12A < IT12)
// 	vol[2] = copy1+(copy2-1)*4;//# of ladders in one module  = 4
//     } else if(id == fIdSens[2]){
// 	vol[0] = 3;
// 	id = gMC->CurrentVolOffID(1,copy);
// 	//detector copy in the ladder = 1<->6  (ITS3 < I302 < I004)
// 	vol[1] = copy;
// 	id = gMC->CurrentVolOffID(2,copy);
// 	//ladder copy in the layer    = 1<->14 (I004 < IT34)
// 	vol[2] = copy;
//     } else if(id == fIdSens[3]){
// 	vol[0] = 4;
// 	id = gMC->CurrentVolOffID(1,copy);
// 	//detector copy in the ladder = 1<->8  (ITS4 < I402 < I005)
// 	vol[1] = copy;
// 	id = gMC->CurrentVolOffID(2,copy);
// 	//ladder copy in the layer    = 1<->22 (I005 < IT34))
// 	vol[2] = copy;
//     }else if(id == fIdSens[4]){
// 	vol[0] = 5;
// 	id = gMC->CurrentVolOffID(1,copy);
// 	//detector copy in the ladder = 1<->22  (ITS5 < I562 < I565)
// 	vol[1] = copy;
// 	id = gMC->CurrentVolOffID(2,copy);
// 	//ladder copy in the layer    = 1<->34 (I565 < IT56)
// 	vol[2] = copy;
//     }else if(id == fIdSens[5]){
// 	vol[0] = 6;
// 	id = gMC->CurrentVolOffID(1,copy);
// 	//detector copy in the ladder = 1<->25  (ITS6 < I566 < I569)
// 	vol[1] = copy;
// 	id = gMC->CurrentVolOffID(2,copy);
// 	//ladder copy in the layer = 1<->38 (I569 < IT56)
// 	vol[2] = copy;
//     } else {
// 	return; // not an ITS volume?
//     } // end if/else if (gMC->CurentVolID(copy) == fIdSens[i])
//     //
//     gMC->TrackPosition(position);
//     gMC->TrackMomentum(momentum);
//     vol[4] = stat0;
//     if(gMC->IsTrackEntering()){
// 	position0 = position;
// 	stat0 = vol[3];
// 	return;
//     } // end if IsEntering
//     // Fill hit structure with this new hit.
    
//     new(lhits[fNhits++]) AliITShit(fIshunt,gAlice->GetMCApp()->GetCurrentTrackNumber(),vol,
// 				   gMC->Edep(),gMC->TrackTime(),position,
// 				   position0,momentum);

//     position0 = position;
//     stat0 = vol[3];

//     return;
// }


//______________________________________________________________________
void AliITSv11::StepManager(){
  //
  //    Called for every step in the ITS, then calles the AliITShit class
  // creator with the information to be recoreded about that hit.
  //
    Int_t         copy, id;
    TLorentzVector position, momentum;
    static TLorentzVector position0;
    static Int_t stat0=0;

    if(!(this->IsActive())){
	return;
    } // end if !Active volume.

    if(!(gMC->TrackCharge())) return;

    id=gMC->CurrentVolID(copy);

    Bool_t sensvol = kFALSE;
    for(Int_t kk=0;kk<6;kk++)if(id == fIdSens[kk])sensvol=kTRUE;
    if(sensvol && (gMC->IsTrackExiting())){
	AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber(), AliTrackReference::kITS);
    } // if Outer ITS mother Volume


    Int_t   copy1,copy2;  
    Int_t   vol[5];
    TClonesArray &lhits = *fHits;
    //
    // Track status
    vol[3] = 0;
    vol[4] = 0;
    if(gMC->IsTrackInside())      vol[3] +=  1;
    if(gMC->IsTrackEntering())    vol[3] +=  2;
    if(gMC->IsTrackExiting())     vol[3] +=  4;
    if(gMC->IsTrackOut())         vol[3] +=  8;
    if(gMC->IsTrackDisappeared()) vol[3] += 16;
    if(gMC->IsTrackStop())        vol[3] += 32;
    if(gMC->IsTrackAlive())       vol[3] += 64;
    //
    // Fill hit structure.
    if(!(gMC->TrackCharge())) return;

    // Only entering charged tracks
    if((id = gMC->CurrentVolID(copy)) == fIdSens[0]) {
	vol[0] = 1;
	id = gMC->CurrentVolOffID(2,copy);
	//detector copy in the ladder = 1<->4  (ITS1 < I101 < I103 < I10A)
	vol[1] = copy;
	gMC->CurrentVolOffID(3,copy1);
	//ladder copy in the module   = 1<->2  (I10A < I12A)
	gMC->CurrentVolOffID(4,copy2);
	//module copy in the layer    = 1<->10 (I12A < IT12)
	vol[2] = copy1+(copy2-1)*2;//# of ladders in one module  = 2

    } else if(id == fIdSens[1]){
	vol[0] = 2;
	id = gMC->CurrentVolOffID(2,copy);
	//detector copy in the ladder = 1<->4  (ITS2 < I1D1 < I1D3 < I20A)
	vol[1] = copy;
	gMC->CurrentVolOffID(3,copy1);
	//ladder copy in the module   = 1<->4  (I20A < I12A)
	gMC->CurrentVolOffID(4,copy2);
	//module copy in the layer    = 1<->10 (I12A < IT12)
	vol[2] = copy1+(copy2-1)*4;//# of ladders in one module  = 4

    } else if(id == fIdSens[2]){
	vol[0] = 3;
	id = gMC->CurrentVolOffID(1,copy);
	//detector copy in the ladder = 1<->6  (ITS3 < I302 < I004)
	vol[1] = copy;
	id = gMC->CurrentVolOffID(2,copy);
	//ladder copy in the layer    = 1<->14 (I004 < IT34)
	vol[2] = copy;

    } else if(id == fIdSens[3]){
	vol[0] = 4;
	id = gMC->CurrentVolOffID(1,copy);
	//detector copy in the ladder = 1<->8  (ITS4 < I402 < I005)
	vol[1] = copy;
	id = gMC->CurrentVolOffID(2,copy);
	//ladder copy in the layer    = 1<->22 (I005 < IT34))
	vol[2] = copy;

    }else if(id == fIdSens[4]){
	vol[0] = 5;
	id = gMC->CurrentVolOffID(1,copy);
	//detector copy in the ladder = 1<->22  (ITS5 < I562 < I565)
	vol[1] = copy;
	id = gMC->CurrentVolOffID(2,copy);
	//ladder copy in the layer    = 1<->34 (I565 < IT56)
	vol[2] = copy;

    }else if(id == fIdSens[5]){
	vol[0] = 6;
	id = gMC->CurrentVolOffID(1,copy);
	//detector copy in the ladder = 1<->25  (ITS6 < I566 < I569)
	vol[1] = copy;
	id = gMC->CurrentVolOffID(2,copy);
	//ladder copy in the layer = 1<->38 (I569 < IT56)
	vol[2] = copy;
    } else {
	return; // not an ITS volume?
    } // end if/else if (gMC->CurentVolID(copy) == fIdSens[i])
    //
    gMC->TrackPosition(position);
    gMC->TrackMomentum(momentum);
    vol[4] = stat0;
    if(gMC->IsTrackEntering()){
	position0 = position;
	stat0 = vol[3];
	return;
    } // end if IsEntering
    // Fill hit structure with this new hit.
    
    new(lhits[fNhits++]) AliITShit(fIshunt,gAlice->GetMCApp()->GetCurrentTrackNumber(),vol,
				   gMC->Edep(),gMC->TrackTime(),position,
				   position0,momentum);

    position0 = position;
    stat0 = vol[3];

    return;
}

