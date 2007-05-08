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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Inner Traking System version Test                                        //
//  This class contains the base procedures for the Inner Tracking System    //
//                                                                           //
// Authors: R. Barbera, B. S. Nilsen.                                        //
// version  Test                                                             //
// Created October 16 1999.                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TSystem.h>
#include <TGeometry.h>
#include <TVirtualMC.h>
#include <TGeoMatrix.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoBBox.h>
#include <TGeoVolume.h>

#include "AliRun.h"
#include "AliMC.h"
#include "AliMagF.h"
#include "AliITSgeom.h"
#include "AliITSgeomSDD.h"
#include "AliITSgeomSPD.h"
#include "AliITSgeomSSD.h"
#include "AliITShit.h"
#include "AliITSvtest.h"

ClassImp(AliITSvtest)

const Double_t AliITSvtest::fgkmicron = 1.0E-4;
const Double_t AliITSvtest::fgkmm = 0.10;
const Double_t AliITSvtest::fgkcm = 1.00;
const Double_t AliITSvtest::fgkDegree = 1.0;
const Double_t AliITSvtest::fgkRadian = 180./3.14159265358979323846;
const Double_t AliITSvtest::fgkgcm3 = 1.0; // assume default is g/cm^3
const Double_t AliITSvtest::fgkCelsius = 1.0; // Assume default is C
const Double_t AliITSvtest::fgkPascal  = 1.0E-3; // Assume kPascal
const Double_t AliITSvtest::fgkKPascal = 1.0;    // Asume kPascal
const Double_t AliITSvtest::fgkeV      = 1.0E-9; // GeV default
const Double_t AliITSvtest::fgkKeV     = 1.0e-6; // GeV default
const Double_t AliITSvtest::fgkMeV     = 1.0e-3; // GeV default
const Double_t AliITSvtest::fgkGeV     = 1.0;    // GeV default

 
//_____________________________________________________________________________
AliITSvtest::AliITSvtest() :
AliITS(),                   // Base Class
fGeomDetOut(kFALSE),       // Flag to write .det file out
fGeomDetIn(kFALSE),         // Flag to read .det file or directly from Geat.
fMajorVersion(IsVersion()), // Major version number == IsVersion
fMinorVersion(-1),          // Minor version number
fEuclidGeomDet(),           // file where detector transormation are define.
fRead(),                    //! file name to read .det file
fWrite(),                   //! file name to write .det file
fIgm()                      //! Geometry initilization object
{
    // Default constructor for the ITS
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.
} 
//_____________________________________________________________________________
AliITSvtest::AliITSvtest(const Char_t *title,Int_t version) :
AliITS("ITS",title),        // Base Class
fGeomDetOut(kFALSE),       // Flag to write .det file out
fGeomDetIn(kFALSE),         // Flag to read .det file or directly from Geat.
fMajorVersion(IsVersion()), // Major version number == IsVersion
fMinorVersion(version),     // Minor version number
fEuclidGeomDet("$ALICE_ROOT/ITS/ITSgeometry_test.det"),// file where detector transormation are define.
fRead("$ALICE_ROOT/ITS/ITSgeometry_test.det"),//! file name to read .det file
fWrite("$ALICE_ROOT/ITS/ITSgeometry_test.det"),//! file name to write .det file
fIgm()                      //! Geometry initilization object
{
    // Default constructor for the ITS. version=1 reads Euclide file for
    // geometry. version=2 use's internal geometry
    // Inputs:
    //   Char_t *title   Geomety title
    //   Int_t   version Minor version number to use.
    //                       =-1 Not defined
    //                       = 1 read Euclid geometry
    //                       = 2 use internal geometry minor verion 2.
    // Outputs:
    //   none.
    // Return:
    //   none.
    Int_t i;

    fIdN    = 6;
    fIdName    = new TString[fIdN];
    fIdName[0] = "ITS1";
    fIdName[1] = "ITS2";
    fIdName[2] = "ITS3";
    fIdName[3] = "ITS4";
    fIdName[4] = "ITS5";
    fIdName[5] = "ITS6";
    fIdSens    = new Int_t[fIdN];
    for(i=0;i<fIdN;i++) fIdSens[i] = 0;
}
//_____________________________________________________________________________
AliITSvtest::~AliITSvtest() {
    // Standard destructor for the ITS
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.
}
//_____________________________________________________________________________
AliITSvtest::AliITSvtest(const char *fileeuc,const char *filetme,
			 const char *name, const char *title) :
AliITS(name, title),        // Base Class
fGeomDetOut(kFALSE),       // Flag to write .det file out
fGeomDetIn(kFALSE),         // Flag to read .det file or directly from Geat.
fMajorVersion(IsVersion()), // Major version number == IsVersion
fMinorVersion(1),           // Minor version number
fEuclidGeomDet("$ALICE_ROOT/ITS/ITSgeometry_test.det"),// file where detector transormation are define.
fRead("$ALICE_ROOT/ITS/ITSgeometry_test.det"),//! file name to read .det file
fWrite("$ALICE_ROOT/ITS/ITSgeometry_test.det"),//! file name to write .det file
fIgm()                      //! Geometry initilization object
{
    // Standard constructor for the ITS
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.
    Int_t i;

    fIdN    = 6;
    fIdName    = new TString[fIdN];
    fIdName[0] = "ITS1";
    fIdName[1] = "ITS2";
    fIdName[2] = "ITS3";
    fIdName[3] = "ITS4";
    fIdName[4] = "ITS5";
    fIdName[5] = "ITS6";
    fIdSens    = new Int_t[fIdN];
    for(i=0;i<fIdN;i++) fIdSens[i] = 0;

    fEuclidMaterial = filetme;
    fEuclidGeometry = fileeuc;
    //  The .det file for the geometry must have the same name as 
    // fileeuc with .euc replaced by .det.
}
//_____________________________________________________________________________
void AliITSvtest::CreateMaterials(){
    // Read materials for the ITS
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.

    switch(GetMinorVersion()){
    case 1:
        CreateMaterialsEuclid();
        break;
    case 2:
        CreateMaterials2();
        break;
    default:
        Warning("CreateMaterials","No CreateMaterials for minor version=%d",
                 GetMinorVersion());
        break;
    } // end switch
    return;
}
//_____________________________________________________________________________
void AliITSvtest::CreateMaterialsEuclid(){
    // Read materials for the ITS
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.
    char *filtmp;

    filtmp = gSystem->ExpandPathName(fEuclidMaterial.Data());
    //  FILE *file = fopen(fEuclidMaterial.Data(),"r");
    FILE *file = fopen(filtmp,"r");
    if(file) {
        fclose(file);
        //    ReadEuclidMedia(fEuclidMaterial.Data(),this);
        ReadEuclidMedia(filtmp);
    } else {
        Error("CreateMaterials"," THE MEDIA FILE %s DOES NOT EXIST !",
              //	  fEuclidMaterial.Data());
              filtmp);
        exit(1);
    } // end if(file)
}
//_____________________________________________________________________________
void AliITSvtest::CreateMaterials2(){
    // Read materials for the ITS
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.
    //TGeoManager *mgr = gGeoManager;
    TGeoMaterial *si,*n2;
    TGeoMedium *sims,*simn,*simN2;
    Double_t params[8];
    Int_t   ifield = gAlice->Field()->Integ();
    Double_t fieldm = gAlice->Field()->Max(); // [kilogauss]

    params[0] = 1.0; // sensitive volume flag
    params[1] = (Double_t)ifield; // magnetic field type
    params[2] = fieldm; // magnetic field stregth
    params[3] = 0.1*fgkDegree; // tmaxfd  Theta max deviation over step
    params[4] = 0.0075*fgkcm; // maximum step size
    params[5] = 0.1; // Maximum fractional energy loss over a step
    params[6] = 1.0E-4*fgkcm; // tracking precision
    params[7] = 0.0*fgkcm; // Minimum step (=0 compute automatically)
                           // must always be =0!

    si = new TGeoMaterial("SI",28.86,14.0,2.33*fgkgcm3,
                          TGeoMaterial::kMatStateSolid,25.0*fgkCelsius,
                          0.0*fgkPascal);
    sims = new TGeoMedium("ITSsensitiveSi",4,si,params);
    params[0] = 0.0; // non sesitive.
    simn = new TGeoMedium("ITSnonsensitiveSi",5,si,params);
    //
    n2 = new TGeoMaterial("Nitrogen Gas",14.00674,7.0,1.250E-3*fgkgcm3,
                          TGeoMaterial::kMatStateGas,25.0*fgkCelsius,
                          101325.0*fgkPascal);
    simN2 = new TGeoMedium("ITSN2",6,n2,params);
    //
    if(sims==0 || simn==0 || simN2==0)
        Error("CreateMaterial","Error getting medium ITSsensitiveSi=%p"
              " ITSnonsensitiveSi=%p ITSN2=%p",sims,simn,simN2);
}
//_____________________________________________________________________________
void AliITSvtest::CreateGeometry(){
    // Read geometry for the ITS
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.

    switch(GetMinorVersion()){
    case 1:
        CreateGeometryEuclid();
        break;
    case 2:
        CreateGeometry2();
        break;
    default:
        Warning("CreateGeometry","No CreateMaterials for minor version=%d",
                 GetMinorVersion());
        break;
    } // end switch
    return;
}
//_____________________________________________________________________________
void AliITSvtest::CreateGeometryEuclid(){
    // Read geometry for the ITS
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.
    char topvol[5];
    char *filtmp;

    filtmp = gSystem->ExpandPathName(fEuclidGeometry.Data());
    FILE *file = fopen(filtmp,"r");
    delete [] filtmp;
    if(file) {
        fclose(file);
        printf("Ready to read Euclid geometry file\n");
        ReadEuclid(fEuclidGeometry.Data(),topvol);
        printf("Read in euclid geometries\n");
    } else {
        Error("CreateGeometry"," THE GEOM FILE %s DOES NOT EXIST !",
              fEuclidGeometry.Data());
        exit(1);
    } // end if(file)
    //
    //---Place the ITS logical volume ITSV in its mother volume (ALIC) 
    //   and make it invisible
    //
    gMC->Gspos("ITSV",1,"ALIC",0,0,0,0,"ONLY");
    //
    //---Outputs the geometry tree in the EUCLID/CAD format
    
    if (fEuclidOut) {
        gMC->WriteEuclid("ITSgeometry", "ITSV", 1, 5);
    } // end if (fEuclidOut)
    cout <<"finished with euclid geometrys"<< endl;
}
//_____________________________________________________________________________
void AliITSvtest::CreateGeometry2(){
    // Test geometry verion 2.
    //  /ALIC_1/ITSV_1/ITSspd1_1/ITS1_1/   lay=1
    //  /ALIC_1/ITSV_1/ITSspd2_1/ITS2_1/   lay=2
    //  /ALIC_1/ITSV_1/ITSsdd1_1/ITS3_1/   lay=3
    //  /ALIC_1/ITSV_1/ITSsdd2_1/ITS4_1/   lay=4
    //  /ALIC_1/ITSV_1/ITSssd1_1/ITS5_1/   lay=5
    //  /ALIC_1/ITSV_1/ITSssd2_1/ITS6_1/   Lay=6
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.
    // These constant character strings are set by cvs during commit
    // do not change them unless you know what you are doing!
    const Char_t *cvsDate="$Date$";
    const Char_t *cvsRevision="$Revision$";
    const Double_t ksensitiveSPD[3]={0.6*fgkcm,0.01*fgkcm,3.5*fgkcm};
    const Double_t ksensitiveSDD[3]={3.5085*fgkcm,0.01499*fgkcm,3.7485*fgkcm};
    const Double_t ksensitiveSSD[3]={3.75*fgkcm,0.0150*fgkcm,2.1*fgkcm};
    const Double_t kwaferSPD[3]={0.7*fgkcm,0.01*fgkcm,3.6*fgkcm};
    const Double_t kwaferSDD[3]={3.61*fgkcm,0.0150*fgkcm,4.38*fgkcm};
    const Double_t kwaferSSD[3]={3.85*fgkcm,0.0150*fgkcm,2.2*fgkcm};
    TGeoManager *mgr = gGeoManager;
    TGeoVolume *vALIC=0;
    TGeoMedium *sensitiveSi=0,*bulckSi=0,*n2=0;
    TGeoBBox *sITS,*sSPD,*sSDD,*sSSD;
    TGeoBBox *sITS1,*sITS2,*sITS3,*sITS4,*sITS5,*sITS6;
    TGeoVolume *vITS,*vSPD1,*vSPD2,*vSDD1,*vSDD2,*vSSD1,*vSSD2;
    TGeoVolume *vITS1,*vITS2,*vITS3,*vITS4,*vITS5,*vITS6;

    vALIC = mgr->GetTopVolume();
    if(vALIC==0) {
        vALIC = mgr->GetVolume("ALIC");
        if(vALIC==0) {
            Error("CreateGeometry2","vALIC=0");
            return;
        }// end if
    } // end if
    //sensitiveSi = mgr->GetMedium("ITSsensitiveSi");
    sensitiveSi = mgr->GetMedium(4);
    //bulckSi     = mgr->GetMedium("ITSnonesensitiveSi");
    bulckSi     = mgr->GetMedium(5);
    //n2          = mgr->GetMedium("ITSN2");
    n2          = mgr->GetMedium(6);
    if(sensitiveSi==0 || bulckSi==0 || n2==0){
        Error("CreateGeometry2","Error getting medium sensitiveSi=%p"
              " bulckSi=%p n2=%p",sensitiveSi,bulckSi,n2);
        TList *lmed = mgr->GetListOfMedia();
        TIter next(lmed);
        while(TGeoMedium *med = (TGeoMedium*) next())
            med->Inspect();
    } //
    sITS1 = new TGeoBBox((Double_t*)ksensitiveSPD);
    vITS1 = new TGeoVolume("ITS1",sITS1,sensitiveSi);
    sITS2 = new TGeoBBox((Double_t*)ksensitiveSPD);
    vITS2 = new TGeoVolume("ITS2",sITS2,sensitiveSi);
    sSPD  = new TGeoBBox((Double_t*)kwaferSPD);
    vSPD1 = new TGeoVolume("ITSspd1",sSPD,bulckSi);
    vSPD2 = new TGeoVolume("ITSspd2",sSPD,bulckSi);
    vSPD1->AddNode(vITS1,1,0); // one copy no translation/rotation
    vSPD2->AddNode(vITS2,1,0); // one copy no translation/rotation
    //
    sITS3 = new TGeoBBox((Double_t*)ksensitiveSDD);
    vITS3 = new TGeoVolume("ITS3",sITS3,sensitiveSi);
    sITS4 = new TGeoBBox((Double_t*)ksensitiveSDD);
    vITS4 = new TGeoVolume("ITS4",sITS4,sensitiveSi);
    sSDD  = new TGeoBBox((Double_t*)kwaferSDD);
    vSDD1 = new TGeoVolume("ITSsdd1",sSDD,bulckSi);
    vSDD2 = new TGeoVolume("ITSsdd2",sSDD,bulckSi);
    vSDD1->AddNode(vITS3,1,0); // one copy no translation/rotation
    vSDD2->AddNode(vITS4,1,0); // one copy no translation/rotation
    //
    sITS5 = new TGeoBBox((Double_t*)ksensitiveSSD);
    vITS5 = new TGeoVolume("ITS5",sITS5,sensitiveSi);
    sITS6 = new TGeoBBox((Double_t*)ksensitiveSSD);
    vITS6 = new TGeoVolume("ITS6",sITS6,sensitiveSi);
    sSSD  = new TGeoBBox((Double_t*)kwaferSSD);
    vSSD1 = new TGeoVolume("ITSssd1",sSSD,bulckSi);
    vSSD2 = new TGeoVolume("ITSssd2",sSSD,bulckSi);
    vSSD1->AddNode(vITS5,1,0); // one copy no translation/rotation
    vSSD2->AddNode(vITS6,1,0); // one copy no translation/rotation
    //
    sITS = new TGeoBBox(100.0,100.,100.0);
    vITS = new TGeoVolume("ITSV",sITS,n2); // one copy of vacume ITSV
    const Int_t length=100;
    Char_t vstrng[length];
    if(fIgm.WriteVersionString(vstrng,length,(AliITSVersion_t)IsVersion(),
                               fMinorVersion,cvsDate,cvsRevision))
        vITS->SetTitle(vstrng);
    else Error("CreateGeometry2","Error writing/setting version string");

    vALIC->AddNode(vITS,1,0);// one copy no translation/rotation
    //
    TGeoTranslation *t1,*t2,*t3,*t4,*t5,*t6;
    TGeoRotation    *r1,*r2,*r3,*r4,*r5,*r6;
    TGeoCombiTrans  *tr1,*tr2,*tr3,*tr4,*tr5,*tr6;
    t1 = new TGeoTranslation( +4.0,+0.0,+0.0); // "perfect" location
    r1 = new TGeoRotation("",90.0,0.0,0.0);    // "perfect" location
    tr1= new TGeoCombiTrans(*t1,*r1);// "perfect" location
    t2 = new TGeoTranslation( +7.0,+0.2,-0.5);
    r2 = new TGeoRotation("",-91.0,10.0,-5.0);
    tr2= new TGeoCombiTrans(*t2,*r2);
    t3 = new TGeoTranslation(+14.9,-0.6,+0.1);
    r3 = new TGeoRotation("",93.0,-7.0,5.0);
    tr3= new TGeoCombiTrans(*t3,*r3);
    t4 = new TGeoTranslation(+23.8,+0.3,-0.2);
    r4 = new TGeoRotation("",91.0,10.0,-5.0);
    tr4= new TGeoCombiTrans(*t4,*r4);
    t5 = new TGeoTranslation(+39.1,+0.1,+0.4);
    r5 = new TGeoRotation("",88.0,1.0,5.0);
    tr5= new TGeoCombiTrans(*t5,*r5);
    t6 = new TGeoTranslation(+43.6,-0.5,+0.2);
    r6 = new TGeoRotation("",92.0,0.0,-5.0);
    tr6= new TGeoCombiTrans(*t6,*r6);
    //
    vITS->AddNode(vSPD1,1,tr1);
    vITS->AddNode(vSPD2,1,tr2);
    vITS->AddNode(vSDD1,1,tr3);
    vITS->AddNode(vSDD2,1,tr4);
    vITS->AddNode(vSSD1,1,tr5);
    vITS->AddNode(vSSD2,1,tr6);
    //
    return;
}
//_____________________________________________________________________________
void AliITSvtest::Init(){
    // Initialise the ITS after it has been created.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.
    Int_t i;

    cout << endl;
    for(i=0;i<29;i++) cout << "*";cout << " ITSvtest_Init ";
    for(i=0;i<28;i++) cout << "*";cout << endl;

    switch(GetMinorVersion()){
    case 1:
        InitEuclid();
        break;
    case 2:
        Init2();
        break;
    default:
        break;
    } // end switch
    UpdateInternalGeometry();
    AliITS::Init();
    if(fGeomDetOut) GetITSgeom()->WriteNewFile(fWrite.Data());

    for(i=0;i<72;i++) cout << "*";
    cout << endl;
}
//_____________________________________________________________________________
void AliITSvtest::InitEuclid(){
    // Initialise the ITS after it has been created. Euclid version
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.
}
//_____________________________________________________________________________
void AliITSvtest::Init2(){
    // Initialise the ITS after it has been created. Geometry 2 verion
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.
    Int_t i,n=3,imed[3];
    TGeoManager *mgr = gGeoManager;
    if(mgr==0) Error("Init2","mgr=0");
    TGeoMedium *sensitiveSi = mgr->GetMedium("ITSsensitiveSi");
    TGeoMedium *bulckSi = mgr->GetMedium("ITSnonsensitiveSi");
    TGeoMedium *n2 = mgr->GetMedium("ITSN2");

    if(sensitiveSi==0) Error("Init2","sensitiveSi=0");
    if(bulckSi==0) Error("Init2","bulckSi=0");
    if(n2==0) Error("Init2","n2=0");
    imed[0] = sensitiveSi->GetId();
    imed[1] = bulckSi->GetId();
    imed[2] = n2->GetId();
    for(i=0;i<n;i++){
        if(imed[i]<=0){
            Error("Init2","GetId failed for imed[i=%d]=%d",i,imed[i]);
            return;
        } // end if
        gMC->Gstpar(imed[i],"CUTGAM",30.0*fgkKeV);
        gMC->Gstpar(imed[i],"CUTELE",30.0*fgkKeV);
        gMC->Gstpar(imed[i],"CUTNEU",30.0*fgkKeV);
        gMC->Gstpar(imed[i],"CUTHAD",30.0*fgkKeV);
        gMC->Gstpar(imed[i],"CUTMUO",30.0*fgkKeV);
        gMC->Gstpar(imed[i],"BCUTE",30.0*fgkKeV);
        gMC->Gstpar(imed[i],"BCUTM",30.0*fgkKeV);
        gMC->Gstpar(imed[i],"DCUTE",30.0*fgkKeV);
        gMC->Gstpar(imed[i],"DCUTM",30.0*fgkKeV);
        //gMC->Gstpar(imed[i],"PPCUTM",);
        //gMC->Gstpar(imed[i],"PAIR",);
        //gMC->Gstpar(imed[i],"COMPT",);
        //gMC->Gstpar(imed[i],"PHOT",);
        //gMC->Gstpar(imed[i],"PFIS",);
        gMC->Gstpar(imed[i],"DRAY",1);
        //gMC->Gstpar(imed[i],"ANNI",);
        //gMC->Gstpar(imed[i],"BREM",);
        //gMC->Gstpar(imed[i],"HADR",);
        //gMC->Gstpar(imed[i],"MUNU",);
        //gMC->Gstpar(imed[i],"DCAY",);
        gMC->Gstpar(imed[i],"LOSS",1);
        //gMC->Gstpar(imed[i],"MULS",);
        //gMC->Gstpar(imed[i],"GHCOR1",);
        //gMC->Gstpar(imed[i],"BIRK1",);
        //gMC->Gstpar(imed[i],"BRIK2",);
        //gMC->Gstpar(imed[i],"BRIK3",);
        //gMC->Gstpar(imed[i],"LABS",);
        //gMC->Gstpar(imed[i],"SYNC",);
        //gMC->Gstpar(imed[i],"STRA",);
    } // end for i
    return;
}
//_____________________________________________________________________________
void AliITSvtest::StepManager(){
    // Called for every step in the ITS
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.
    //
    // Fill hit structure.
    static TLorentzVector position, momentum; // Saves on calls to construtors
    static AliITShit hit;                     // Saves on calls to construtors
    Int_t cpn0,cpn1,cpn2,status,lay,mod,id;

    if(!(gMC->TrackCharge())) return;
    if(!(this->IsActive())) return;
    TClonesArray &lhits = *(Hits());
    // Track status
    status = 0;
    if(gMC->IsTrackInside())      status +=  1;
    if(gMC->IsTrackEntering())    status +=  2;
    if(gMC->IsTrackExiting())     status +=  4;
    if(gMC->IsTrackOut())         status +=  8;
    if(gMC->IsTrackDisappeared()) status += 16;
    if(gMC->IsTrackStop())        status += 32;
    if(gMC->IsTrackAlive())       status += 64;
    // Only entering charged tracks
    id = gMC->CurrentVolID(cpn0);
    for(lay=0;lay<6;lay++) if(id == fIdSens[lay]) break;
    lay++;
    if(lay>6) return; // not in detector
    cpn0=cpn1=cpn2=1;
    fIgm.DecodeDetector(mod,lay,cpn0,cpn1,cpn2);
    // Fill hit structure.
    hit.SetModule(mod);
    hit.SetTrack(gAlice->GetMCApp()->GetCurrentTrackNumber());
    gMC->TrackPosition(position);
    gMC->TrackMomentum(momentum);
    hit.SetPosition(position);
    hit.SetTime(gMC->TrackTime());
    hit.SetMomentum(momentum);
    hit.SetStatus(status);
    hit.SetEdep(gMC->Edep());
    hit.SetShunt(GetIshunt());
    if(gMC->IsTrackEntering()){
        hit.SetStartPosition(position);
        hit.SetStartTime(gMC->TrackTime());
        hit.SetStartStatus(status);
        return; // don't save entering hit.
    } // end if IsEntering
    // Fill hit structure with this new hit.
    new(lhits[fNhits++]) AliITShit(hit); // Use Copy Construtor.
    // Save old position... for next hit.
    hit.SetStartPosition(position);
    hit.SetStartTime(gMC->TrackTime());
    hit.SetStartStatus(status);
    //
    Double_t g0[4],l0[4],g1[4];
    position.GetXYZT(g0);
    gMC->Gmtod(g0,l0,1); // flag=1 convert coordiantes
    gMC->Gdtom(l0,g1,1); // flag=1 convert coordinates
    printf("    gMC: mod=%d g=%g %g %g %g -> l= %g %g %g %g ->g=%g %g %g %g\n",
           mod,g0[0],g0[1],g0[2],g0[3],l0[0],l0[1],l0[2],l0[3],g1[0],g1[2],g1[2],g1[3]);
    GetITSgeom()->GtoL(mod,g0,l0);
    GetITSgeom()->LtoG(mod,l0,g1);
    printf("ITSgeom: mod=%d g=%g %g %g %g -> l= %g %g %g %g ->g=%g %g %g %g\n",
           mod,g0[0],g0[1],g0[2],g0[3],l0[0],l0[1],l0[2],l0[3],g1[0],g1[2],g1[2],g1[3]);
    TGeoNode *cur = gGeoManager->GetCurrentNode();
    cur->MasterToLocal(g0,l0);
    cur->LocalToMaster(l0,g1);
    printf("   TGeo: mod=%d g=%g %g %g %g -> l= %g %g %g %g ->g=%g %g %g %g\n",
           mod,g0[0],g0[1],g0[2],g0[3],l0[0],l0[1],l0[2],l0[3],g1[0],g1[2],g1[2],g1[3]);
    printf("=====================\n");
    //

  return;
}
//______________________________________________________________________
void AliITSvtest::PrintAscii(ostream *os)const{
    // Print out class data values in Ascii Form to output stream
    // Inputs:
    //   ostream *os   Output stream where Ascii data is to be writen
    // Outputs:
    //   none.
    // Return:
    //   none.
#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC || defined __xlC__
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif

    *os << fGeomDetOut << " " << fGeomDetIn  << " ";
    *os << fMajorVersion << " " << fMinorVersion << " ";
    *os << "\"" << fEuclidGeomDet.Data() << "\"" << " ";
    *os << "\"" << fRead.Data() << "\"" << " ";
    *os << "\"" << fWrite.Data() << "\"" << " ";
    fmt = os->setf(ios::scientific); // set scientific floating point output
    os->flags(fmt); // reset back to old Formating.
    return;
}
//______________________________________________________________________
void AliITSvtest::ReadAscii(istream *is){
    // Read in class data values in Ascii Form to output stream
    // Inputs:
    //   istream *is   Input stream where Ascii data is to be read in from
    // Outputs:
    //   none.
    // Return:
    //   none.
    Char_t name[120];

    *is >> fGeomDetOut >> fGeomDetIn ;
    *is >> fMajorVersion >> fMinorVersion;
    *is >> name;
    fEuclidGeomDet = name;
    *is >> name;
    fRead = name;
    *is >> name;
    fWrite = name;
    fIgm.SetVersion((AliITSVersion_t)fMajorVersion,fMinorVersion);
    fIgm.SetGeometryName("ITS test geometry");
}
//______________________________________________________________________
ostream &operator<<(ostream &os,const AliITSvtest &s){
    // Standard output streaming function
    // Inputs:
    //   ostream            &os  output steam
    //   AliITSvtest &s class to be streamed.
    // Output:
    //   none.
    // Return:
    //   ostream &os  The stream pointer

    s.PrintAscii(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &is,AliITSvtest &s){
    // Standard inputput streaming function
    // Inputs:
    //   istream            &is  input steam
    //   AliITSvtest &s class to be streamed.
    // Output:
    //   none.
    // Return:
    //   ostream &os  The stream pointer

    s.ReadAscii(&is);
    return is;
}
