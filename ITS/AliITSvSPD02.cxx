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

#include <TGeometry.h>
#include <TNode.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TBRIK.h>
#include <TGeoMatrix.h>
#include <TVirtualMC.h>

#include "AliRun.h"
#include "AliMagF.h"
#include "AliTrackReference.h"

#include "AliITShit.h"
#include "AliITS.h"
#include "AliITSvSPD02.h"
#include "AliITSgeom.h"
#include "AliITSgeomSPD.h"
#include "AliITSDetTypeSim.h"
#include "AliITSCalibrationSPD.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSsimulationSPD.h"
#include "AliMC.h"


///////////////////////////////////////////////////////////////////////
// Step manager and 
// geometry class
// for the ITS 
// SPD test beam
// geometry of summer 2002
// 
///////////////////////////////////////////////////////////////////////
ClassImp(AliITSvSPD02)

//______________________________________________________________________
AliITSvSPD02::AliITSvSPD02():
AliITS(),
fGeomDetOut(kFALSE),
fGeomDetIn(kFALSE),
fMajorVersion(1),
fMinorVersion(2),
fGeomNumber(2002),
fEuclidGeomDet(),
fRead(),
fWrite(),
fDet1(300.0),
fDet2(300.0),
fChip1(300.0),
fChip2(300.0){
    ////////////////////////////////////////////////////////////////////////
    // Standard default constructor for the ITS SPD test beam 2002 version 1.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    A default created class.
    ////////////////////////////////////////////////////////////////////////
    Int_t i;

    for(i=0;i<60;i++) fRead[i] = '\0';
    for(i=0;i<60;i++) fWrite[i] = '\0';
    for(i=0;i<60;i++) fEuclidGeomDet[i] = '\0';
}
//______________________________________________________________________
AliITSvSPD02::AliITSvSPD02(const char *title,Int_t gn) : AliITS("ITS", title),
fGeomDetOut(kFALSE),
fGeomDetIn(kFALSE),
fMajorVersion(1),
fMinorVersion(2),
fGeomNumber(2002),
fEuclidGeomDet(),
fRead(),
fWrite(),
fDet1(300.0),
fDet2(300.0),
fChip1(300.0),
fChip2(300.0){
    ////////////////////////////////////////////////////////////////////////
    //    Standard constructor for the ITS SPD testbeam 2002 version 1.
    // Inputs:
    //    const char *title    title for this ITS geometry.
    //    Int_t      gn        Geometry version number (year) default 2002.
    // Outputs:
    //    none.
    // Return:
    //    A standard created class.
    ////////////////////////////////////////////////////////////////////////
    Int_t i;

    fGeomNumber = gn;
    fIdN = 2;
    fIdName = new TString[fIdN];
    fIdName[0] = "IMBS";
    fIdName[1] = "ITST";
    fIdSens    = new Int_t[fIdN];
    for(i=0;i<fIdN;i++) fIdSens[i] = 0;
    SetThicknessDet1();
    SetThicknessDet2();
    SetThicknessChip1();
    SetThicknessChip2();	 	 	 

    fEuclidGeometry="$ALICE_ROOT/ITS/ITSgeometry_vSPD022.euc";
    strncpy(fEuclidGeomDet,"$ALICE_ROOT/ITS/ITSgeometry_vSPD022.det",60);
    strncpy(fRead,fEuclidGeomDet,60);
    strncpy(fWrite,fEuclidGeomDet,60);
}
//______________________________________________________________________
AliITSvSPD02::AliITSvSPD02(const AliITSvSPD02 &source) : AliITS(source){
    ////////////////////////////////////////////////////////////////////////
    //     Copy Constructor for ITS SPD test beam 2002 version 1.
    // This class is not to be copied. Function only dummy.
    // Inputs:
    //    const AliITSvSPD02 &source   The class to be copied
    // Outputs:
    //    none.
    // Return:
    //    A warning message.
    ////////////////////////////////////////////////////////////////////////
    if(&source == this) return;
    Warning("Copy Constructor","Not allowed to copy AliITSvSPD02");
    return;
}
//______________________________________________________________________
AliITSvSPD02& AliITSvSPD02::operator=(const AliITSvSPD02 &source){
    ////////////////////////////////////////////////////////////////////////
    //    Assignment operator for the ITS SPD test beam 2002 version 1.
    // This class is not to be copied. Function only dummy.
    // Inputs:
    //    const AliITSvSPD02 &source   The class to be copied
    // Outputs:
    //    none.
    // Return:
    //    A Warning message
    ////////////////////////////////////////////////////////////////////////
    if(&source == this) return *this;
    Warning("= operator","Not allowed to copy AliITSvSPD02");
    return *this;
}
//______________________________________________________________________
AliITSvSPD02::~AliITSvSPD02() {
    ////////////////////////////////////////////////////////////////////////
    //    Standard destructor for the ITS SPD test beam 2002 version 1.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    ////////////////////////////////////////////////////////////////////////
}
//______________________________________________________________________
void AliITSvSPD02::BuildGeometry(){
    ////////////////////////////////////////////////////////////////////////
    //    Geometry builder for the ITS SPD test beam 2002 version 1.
    //    ALIC    ALICE Mother Volume
    //     |- ITSV     ITS Mother Volume
    //         |- IDET       Detector under Test
    //         |   |- ITS0       SPD Si Chip
    //         |   |  |- ITST      SPD Sensitivve Volume
    //         |   |- IPC0 *5    Readout chip
    //         |- ITEL *4    SPD Telescope
    //             |- IMB0       SPD Si Chip
    //             |   |- IMBS     SPD Sensitive volume
    //             |- ICMB       Chip MiniBus.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    ////////////////////////////////////////////////////////////////////////
    // Get the top alice volume.

    switch (fGeomNumber){
    case 2002:
        BuildGeometry2002();
        break;
    default:
        BuildGeometry2002();
        break;
    } // end switch
}
//______________________________________________________________________
void AliITSvSPD02::BuildGeometry2002(){
    ////////////////////////////////////////////////////////////////////////
    //    Geometry builder for the ITS SPD test beam 2002 version 1.
    //    ALIC    ALICE Mother Volume
    //     |- ITSV     ITS Mother Volume
    //         |- IDET       Detector under Test
    //         |   |- ITS0       SPD Si Chip
    //         |   |  |- ITST      SPD Sensitivve Volume
    //         |   |- IPC0 *5    Readout chip
    //         |- ITEL *4    SPD Telescope
    //             |- IMB0       SPD Si Chip
    //             |   |- IMBS     SPD Sensitive volume
    //             |- ICMB       Chip MiniBus.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    ////////////////////////////////////////////////////////////////////////
    // Get the top alice volume.
    TNode *aALIC = gAlice->GetGeometry()->GetNode("alice");
    aALIC->cd();

    // Define ITS Mother Volume
    Float_t data[3];
    Float_t ddettest=200.0E-4,ddettelescope=300.0E-4;
    Float_t dchipMiniBus=750.0E-4,dchiptest=300.0E-4;
    //Float_t yposition= 0.0;
    TRotMatrix *r0 = new TRotMatrix("ITSidrotm0","ITSidrotm0",
                                    90.0,0,0.0,0,90.0,270.0);
    data[0] = 10.0;
    data[1] = 50.0;
    data[2] = 100.0;
    TBRIK *iITSVshape = new TBRIK("ITSVshape",
                                  "ITS Logical Mother Volume","Air",
                                  data[0],data[1],data[2]);
    TNode *iITSV = new TNode("ITSV","ITS Mother Volume",iITSVshape,
                             0.0,0.0,0.0,0,0);
    iITSV->cd(); // set ourselve into ITSV subvolume of aALIC

    // SPD part of telescope (MiniBuS)
    data[0] = 0.705;
    data[1] = 0.5*ddettelescope;
    data[2] = 3.536;
    TBRIK *iIMB0shape = new TBRIK("IMB0shape","SPD wafer","Si",
				 data[0],data[1],data[2]);
    Float_t detMiniBusX,detMiniBusY,detMiniBusZ;
    data[0] = detMiniBusX = 0.64;
    data[1] = detMiniBusY = 0.5*ddettelescope;
    data[2] = detMiniBusZ = 3.48;
    TBRIK *iIMBSshape = new TBRIK("IMBSshape","SPD Sensitive volume","Si",
				 data[0],data[1],data[2]);
    Float_t chipMiniBusX,chipMiniBusY,chipMiniBusZ;
    data[0] = chipMiniBusX = 0.793;
    data[1] = chipMiniBusY = 0.5*dchipMiniBus;
    data[2] = chipMiniBusZ = 0.68;
    TBRIK *iICMBshape = new TBRIK("ICMBshape","chip Minibus","Si",
				 data[0],data[1],data[2]);
    data[0] = TMath::Max(detMiniBusX,chipMiniBusX);
    data[1] = detMiniBusY+chipMiniBusY;
    data[2] = TMath::Max(detMiniBusZ,chipMiniBusZ);
    TBRIK *iITELshape = new TBRIK("ITELshape","ITELshape","Air",
				 data[0],data[1],data[2]);

    // SPD under test
    Float_t spdX,spdY,spdZ,spdchipX,spdchipY,spdchipZ;
    data[0] = 0.705;
    data[1] = ddettest;
    data[2] = 3.536;
    TBRIK *iITS0shape = new TBRIK("ITS0shape","SPD wafer","Si",
				 data[0],data[1],data[2]); // contains detector
    data[0] = spdX = 0.64;
    data[1] = spdY = ddettest;
    data[2] = spdZ = 3.48;
    TBRIK *iITSTshape = new TBRIK("ITSTshape","SPD sensitive volume","Si",
				 data[0],data[1],data[2]);
    // ITS0 with no translation and unit rotation matrix.
    data[0] = spdchipX = 0.793;
    data[1] = spdchipY = dchiptest;
    data[2] = spdchipZ = 0.68;
    TBRIK *iIPC0shape = new TBRIK("IPC0shape","Readout Chips","Si",
				 data[0],data[1],data[2]); // chip under test
    data[0] = TMath::Max(spdchipX,spdX);
    data[1] = spdY+spdchipY;
    data[2] = TMath::Max(spdchipZ,spdZ);
    TBRIK *iIDETshape = new TBRIK("IDETshape","Detector Under Test","Air",
				 data[0],data[1],data[2]);
    // Place volumes in geometry
    Int_t i,j;
    char name[20],title[50];
    Double_t px=0.0,py=0.0,pz[4]={-38.0,0.0,0.0,0.0};
    pz[1] = pz[0]+2.0;
    pz[2] = pz[1]+38.0+spdY+spdchipY+34.5;
    pz[3] = pz[2]+2.0;
    TNode *iITEL[4],*iICMB[4],*iIMB0[4],*iIMBS[4];
    TNode *iIDET = new TNode("IDET","Detector Under Test",iIDETshape,
			    0.0,0.0,pz[1]+38.0,r0,0);
    iIDET->cd();
    TNode *iITS0 = new TNode("ITS0","SPD Chip",iITS0shape,
			    0.0,iIDETshape->GetDy()-spdY,0.0,0,0);
    TNode *iIPC0[5];
    for(i=0;i<5;i++) { //place readout chips on the back of SPD chip under test
        sprintf(name,"IPC0%d",i);
        sprintf(title,"Readout chip #%d",i+1);
        j = i-2;
        iIPC0[i] = new TNode(name,title,iIPC0shape,
                             0.0,spdchipY-iIDETshape->GetDy(),
                             j*2.0*spdchipZ+j*0.25*(spdZ-5.*spdchipZ),0,0);
    } // end for i
    iITS0->cd();
    TNode *iITST = new TNode("ITST","SPD sensitive volume",iITSTshape,
                             0.0,0.0,0.0,0,0);
    for(Int_t i=0;i<4;i++){
        iITSV->cd();
        sprintf(name,"ITEL%d",i);
        sprintf(title,"Test beam telescope element #%d",i+1);
        iITEL[i] = new TNode(name,title,iITELshape,px,py,pz[i],r0,0);
        iITEL[i]->cd();
        iICMB[i] = new TNode("ICMB","Chip MiniBus",iICMBshape,
                             0.0,-iITELshape->GetDy()+detMiniBusY,0.0,0,0);
        iIMB0[i] = new TNode("IMB0","Chip MiniBus",iIMB0shape,
                             0.0, iITELshape->GetDy()-detMiniBusY,0.0,0,0);
        iIMB0[i]->cd();
        iIMBS[i] = new TNode("IMBS","IMBS",iIMBSshape,0.0,0.0,0.0,0,0);
        // place IMBS inside IMB0 with no translation and unit rotation matrix.
    } // end for i
    aALIC->cd();
    iITST->SetLineColor(kYellow);
    fNodes->Add(iITST);
    for(i=0;i<4;i++){
        iIMBS[i]->SetLineColor(kGreen);
        fNodes->Add(iIMBS[i]);
    } // end for i
}
//______________________________________________________________________
void AliITSvSPD02::CreateGeometry(){
    ////////////////////////////////////////////////////////////////////////
    //  This routine defines and Creates the geometry for version 1 of the ITS.
    //    ALIC    ALICE Mother Volume
    //     |- ITSV     ITS Mother Volume
    //         |- IDET       Detector under Test
    //         |   |- ITS0       SPD Si Chip
    //         |   |  |- ITST      SPD Sensitivve Volume
    //         |   |- IPC0 *5    Readout chip
    //         |- ITEL *4    SPD Telescope
    //             |- IMB0       SPD Si Chip
    //             |   |- IMBS     SPD Sensitive volume
    //             |- ICMB       Chip MiniBus.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    ////////////////////////////////////////////////////////////////////////

    switch (fGeomNumber){
    case 2002:
        CreateGeometry2002();
        break;
    default:
        CreateGeometry2002();
        break;
    } // end switch
}
//______________________________________________________________________
void AliITSvSPD02::CreateGeometry2002(){
    ////////////////////////////////////////////////////////////////////////
    //  This routine defines and Creates the geometry for version 1 of the ITS.
    //    ALIC    ALICE Mother Volume
    //     |- ITSV     ITS Mother Volume
    //         |- IDET       Detector under Test
    //         |   |- ITS0       SPD Si Chip
    //         |   |  |- ITST      SPD Sensitivve Volume
    //         |   |- IPC0 *5    Readout chip
    //         |- ITEL *4    SPD Telescope
    //             |- IMB0       SPD Si Chip
    //             |   |- IMBS     SPD Sensitive volume
    //             |- ICMB       Chip MiniBus.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    ////////////////////////////////////////////////////////////////////////
    Float_t data[49];
    // Define media off-set
    Int_t *idtmed = fIdtmed->GetArray()+1; // array of media indexes
    Int_t idrotm[4]; // Array of rotation matrix indexes
    Float_t ddettest=200.0E-4,ddettelescope=300.0E-4;
    Float_t dchipMiniBus=750.0E-4,dchiptest=300.0E-4;
    Float_t yposition= 0.0;

    if(gMC==0) return;
    // Define Rotation-reflextion Matrixes needed
    // 0 is the unit matrix
    AliMatrix(idrotm[0], 90.0,0.0, 0.0,0.0, 90.0,270.0);
    data[0] = 10.0;
    data[1] = 50.0;
    data[2] = 100.0;
    gMC->Gsvolu("ITSV","BOX",idtmed[0],data,3);
    gMC->Gspos("ITSV",1,"ALIC",0.0,0.0,0.0,0,"ONLY");

    //cout << "idtmed[0]=" << idtmed[0]<<endl;
    //cout << "idtmed[1]=" << idtmed[1]<<endl;
    Float_t detMiniBusX,detMiniBusY,detMiniBusZ;
    // SPD part of telescope (MiniBuS)
    data[0] = detMiniBusX = 0.705;
    data[1] = detMiniBusY = 0.5*ddettelescope;
    data[2] = detMiniBusZ = 3.536;
    gMC->Gsvolu("IMB0", "BOX ", idtmed[1], data, 3);   // contains detector
    data[0] = 0.64;
    data[1] = 0.5*ddettelescope;
    data[2] = 3.48;
    gMC->Gsvolu("IMBS","BOX ",idtmed[1],data,3); // sensitive detecor volulme
    gMC->Gspos("IMBS",1,"IMB0",0.0,0.0,0.0,0,"ONLY"); // place IMBS inside
    // IMB0 with no translation and unit rotation matrix.
    Float_t chipMiniBusX,chipMiniBusY,chipMiniBusZ;
    data[0] = chipMiniBusX = 0.793;
    data[1] = chipMiniBusY = 0.5*dchipMiniBus;
    data[2] = chipMiniBusZ = 0.68;
    gMC->Gsvolu("ICMB","BOX ",idtmed[1],data, 3);   // chip Minibus
    data[0] = TMath::Max(detMiniBusX,chipMiniBusX);
    data[1] = detMiniBusY+chipMiniBusY;
    data[2] = TMath::Max(detMiniBusZ,chipMiniBusZ);
    gMC->Gsvolu("ITEL","BOX ",idtmed[0],data,3);
    gMC->Gspos("IMB0",1,"ITEL",0.0,data[1]-detMiniBusY,0.0,0,"ONLY");
    gMC->Gspos("ICMB",1,"ITEL",0.0,-data[1]+chipMiniBusY,0.0,0,"ONLY");

    // SPD under test
    Float_t spdX,spdY,spdZ,spdchipX,spdchipY,spdchipZ;
    data[0] = spdX = 0.705;
    data[1] = spdY = 0.5*ddettest;
    data[2] = spdZ = 3.536;
    gMC->Gsvolu("ITS0", "BOX ", idtmed[1], data, 3);   // contains detector
    data[0] = 0.64;
    data[1] = 0.5*ddettest;
    data[2] = 3.48;
    gMC->Gsvolu("ITST","BOX ",idtmed[1],data,3);// sensitive detecor volume
    gMC->Gspos("ITST",1,"ITS0",0.0,0.0,0.0,0,"ONLY"); // place ITST inside
    // ITS0 with no translation and unit rotation matrix.
    data[0] = spdchipX = 0.793;
    data[1] = spdchipY = 0.5*dchiptest;
    data[2] = spdchipZ = 0.68;
    gMC->Gsvolu("IPC0", "BOX ", idtmed[1],data,3);   // chip under test
    data[0] = TMath::Max(spdchipX,spdX);
    data[1] = spdY+spdchipY;
    data[2] = TMath::Max(spdchipZ,spdZ);
    gMC->Gsvolu("IDET","BOX ",idtmed[0],data,3);
    gMC->Gspos("ITS0",1,"IDET",0.0,data[1]-spdY,0.0,0,"ONLY");
    for(Int_t i=-2;i<3;i++) gMC->Gspos("IPC0",i+3,"IDET",0.0,-data[1]+spdchipY,
              i*2.*spdchipZ+i*0.25*(spdZ-5.*spdchipZ),0,"ONLY");

    // Positions detectors, Beam Axis Z, X to the right, Y up to the sky.
    Float_t p00X,p00Y,p00Z,p01X,p01Y,p01Z,p10X,p10Y,p10Z,p11X,p11Y,p11Z;
    p00X = 0.0;
    p00Y = 0.0;
    p00Z = -38.0;
    gMC->Gspos("ITEL",1,"ITSV",p00X,p00Y,p00Z,idrotm[0],"ONLY");
    p01X = 0.0;
    p01Y = 0.0;
    p01Z = p00Z+2.0;
    gMC->Gspos("ITEL",2,"ITSV",p01X,p01Y,p01Z,idrotm[0],"ONLY");
    Float_t pdetX,pdetY,pdetZ;
    pdetX = 0.0;
    pdetY = 0.0+yposition;
    pdetZ = p01Z+38.0;
    gMC->Gspos("IDET",1,"ITSV",pdetX,pdetY,pdetZ,idrotm[0],"ONLY");
    p10X = 0.0;
    p10Y = 0.0;
    p10Z = pdetZ + 34.5;
    gMC->Gspos("ITEL",3,"ITSV",p10X,p10Y,p10Z,idrotm[0],"ONLY");
    p11X = 0.0;
    p11Y = 0.0;
    p11Z = p10Z+2.0;
    gMC->Gspos("ITEL",4,"ITSV",p11X,p11Y,p11Z,idrotm[0],"ONLY");
}
//______________________________________________________________________
void AliITSvSPD02::CreateMaterials(){
    ////////////////////////////////////////////////////////////////////////
    //
    // Create ITS SPD test beam materials
    //     This function defines the default materials used in the Geant
    // Monte Carlo simulations for the geometries AliITSv1, AliITSv3,
    // AliITSvSPD02.
    // In general it is automatically replaced by
    // the CreateMaterials routine defined in AliITSv?. Should the function
    // CreateMaterials not exist for the geometry version you are using this
    // one is used. See the definition found in AliITSv5 or the other routine
    // for a complete definition.
    //
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    /////////////////////////////////////////////////////////////////////////

    switch (fGeomNumber){
    case 2002:
        CreateMaterials2002();
        break;
    default:
        CreateMaterials2002();
        break;
    } // end switch
}
//______________________________________________________________________
void AliITSvSPD02::CreateMaterials2002(){
    ////////////////////////////////////////////////////////////////////////
    //
    // Create ITS SPD test beam materials
    //     This function defines the default materials used in the Geant
    // Monte Carlo simulations for the geometries AliITSv1, AliITSv3,
    // AliITSvSPD02.
    // In general it is automatically replaced by
    // the CreateMaterials routine defined in AliITSv?. Should the function
    // CreateMaterials not exist for the geometry version you are using this
    // one is used. See the definition found in AliITSv5 or the other routine
    // for a complete definition.
    //
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    /////////////////////////////////////////////////////////////////////////
    Float_t tmaxfdSi = 0.1; // Degree
    Float_t stemaxSi = 0.0075; // cm
    Float_t deemaxSi = 0.1; // Fraction of particle's energy 0<deemax<=1
    Float_t epsilSi  = 1.0E-4;//
    Float_t stminSi  = 0.0; // cm "Default value used"

    Float_t tmaxfdAir = 0.1; // Degree
    Float_t stemaxAir = .10000E+01; // cm
    Float_t deemaxAir = 0.1; // Fraction of particle's energy 0<deemax<=1
    Float_t epsilAir  = 1.0E-4;//
    Float_t stminAir  = 0.0; // cm "Default value used"
    Int_t   ifield = gAlice->Field()->Integ();
    Float_t fieldm = gAlice->Field()->Max();

    AliMaterial(1,"AIR$",0.14610E+02,0.73000E+01,0.12050E-02,
                0.30423E+05,0.99900E+03);
    AliMedium(1,"AIR$",1,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,
              epsilAir,stminAir);

    AliMaterial(2,"SI$",0.28086E+02,0.14000E+02,0.23300E+01,
                0.93600E+01,0.99900E+03);
    AliMedium(2,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
              epsilSi,stminSi);
}
//______________________________________________________________________
Int_t AliITSvSPD02::DecodeDetector(Int_t id,Int_t cpy,Int_t &lay,Int_t &lad,
    Int_t &det)const{
    //     Based on the geometry tree defined in Geant 3.21, this
    // routine initilizes the Class AliITSgeom from the Geant 3.21 ITS geometry
    // sturture.
    // Inputs:
    //    Int_t id   Detector volume id
    //    Int_t cpy  Detector copy number
    // Outputs:
    //    Int_t lay  layer number
    //    Int_t lad  ladder number
    //    Int_t det  detector number
    // Return:
    //    The module number
    Int_t mod;

    lad = det = 1;
    lay = cpy;
    if(cpy>2 && id==fIdSens[0]) lay = cpy + 1;
    if(id==fIdSens[1]) lay = 3;
    mod = lay - 1;
    return mod;
}
//______________________________________________________________________
void AliITSvSPD02::InitAliITSgeom(){
    //     Based on the geometry tree defined in Geant 3.21, this
    // routine initilizes the Class AliITSgeom from the Geant 3.21 ITS geometry
    // sturture.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
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

    par[0]=0.64;par[1]=0.5*300.0E-4;par[2]=3.48;
    mod=5;
    AliITSgeom* geom = new AliITSgeom(0,knlayers,knlad,kndet,mod);
    SetITSgeom(geom);
    for(i=0;i<kltypess;i++)for(cpy=1;cpy<=kitsGeomTreeCopys[i];cpy++){
        path.Form(knames[i].Data(),cpy);
        gMC->GetTransformation(path.Data(),matrix);
        gMC->GetShape(path.Data(),shapeName,shapePar);
        shapeParF.Set(shapePar.GetSize());
        for(j=0;j<shapePar.GetSize();j++) shapeParF[j]=shapePar[j];
        mod = DecodeDetector(fIdSens[i],cpy,lay,lad,det);
        geom->CreateMatrix(mod,lay,lad,det,kSPD,trans,rot);
        geom->SetTrans(mod,matrix.GetTranslation());
        geom->SetRotMatrix(mod,matrix.GetRotationMatrix());
        geom->GetGeomMatrix(mod)->SetPath(path.Data());
        if(!(geom->IsShapeDefined((Int_t)kSPD)))
            geom->ReSetShape(kSPD,new AliITSgeomSPD425Short(npar,par));
    } // end for i,cpy/
    return;
}
//______________________________________________________________________
void AliITSvSPD02::Init(){
    ////////////////////////////////////////////////////////////////////////
    //     Initialise the ITS after it has been created.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    ////////////////////////////////////////////////////////////////////////
    Int_t i;

    cout << endl;
    for(i=0;i<26;i++) cout << "*";
    cout << " ITSvSPD02" << fMinorVersion << "_Init ";
    for(i=0;i<25;i++) cout << "*";cout << endl;

    if(fRead[0]=='\0') strncpy(fRead,fEuclidGeomDet,60);
    if(fWrite[0]=='\0') strncpy(fWrite,fEuclidGeomDet,60);
    if(GetITSgeom()!=0) SetITSgeom(0x0);
    AliITSgeom* geom = new AliITSgeom();
    SetITSgeom(geom);
    if(fGeomDetIn) GetITSgeom()->ReadNewFile(fRead);
    if(!fGeomDetIn) this->InitAliITSgeom();
    if(fGeomDetOut) GetITSgeom()->WriteNewFile(fWrite);
    AliITS::Init();

    for(i=0;i<72;i++) cout << "*";
    cout << endl;
    if(gMC) fIDMother = gMC->VolId("ITSV"); // ITS Mother Volume ID.
    else fIDMother = 0;
}
//______________________________________________________________________
void AliITSvSPD02::SetDefaults(){
    // sets the default segmentation, response, digit and raw cluster classes
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    const Float_t kconv = 1.0e+04; // convert cm to microns

    Info("SetDefaults","Setting up only SPD detector");

    if(!fDetTypeSim) fDetTypeSim = new AliITSDetTypeSim();
    fDetTypeSim->SetITSgeom(GetITSgeom());
    AliITSgeomSPD  *s0;
    Int_t i;
    Float_t bx[256],bz[280];
    fDetTypeSim->ResetCalibrationArray();
    fDetTypeSim->ResetSegmentation();
 
    //SPD
    // Get shape info. Do it this way for now.
    s0 = (AliITSgeomSPD*) GetITSgeom()->GetShape(kSPD);
    AliITSCalibration *resp0=new AliITSCalibrationSPD();
    resp0->SetTemperature();
    resp0->SetDistanceOverVoltage();
    SetCalibrationModel(0,resp0); 
	
    AliITSsegmentationSPD *seg0=new AliITSsegmentationSPD();
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
    const char *kData0=(fDetTypeSim->GetCalibrationModel(0))->DataType();
    if (strstr(kData0,"real")) fDetTypeSim->SetDigitClassName(kSPD,"AliITSdigit");
    else fDetTypeSim->SetDigitClassName(kSPD,"AliITSdigitSPD");
//    SetSimulationModel(kSPD,new AliITSsimulationSPDdubna(seg0,resp0));
//    iDetType->ReconstructionModel(new AliITSClusterFinderSPD());
   
/*
    SetResponseModel(kSDD,new AliITSCalibrationSDD());
    SetSegmentationModel(kSDD,new AliITSsegmentationSDD());
    DetType(kSDD)->ClassNames("AliITSdigitSDD","AliITSRawClusterSDD");

    SetResponseModel(kSSD,new AliITSCalibrationSSD());
    SetSegmentationModel(kSSD,new AliITSsegmentationSSD());
    DetType(kSSD)->ClassNames("AliITSdigitSSD","AliITSRawClusterSSD");
*/
    if(fgkNTYPES>3){
	Warning("SetDefaults",
		"Only the four basic detector types are initialised!");
    }// end if
    return;
}
//______________________________________________________________________
void AliITSvSPD02::SetDefaultSimulation(){
    // sets the default simulation.
    // Inputs:
    //      none.
    // Outputs:
    //      none.
    // Return:
    //      none.

  if(GetITSgeom()==0){
    Warning("SetDefaultSimulation","ITS geometry is null!");
    return;
  }

  if(!fDetTypeSim) fDetTypeSim = new AliITSDetTypeSim();
  AliITSsimulation *sim;
  //AliITSsegmentation *seg;
  //AliITSCalibration *res;
  if(fDetTypeSim){
    sim = fDetTypeSim->GetSimulationModel(kSPD);
    if (!sim) {
      //seg = (AliITSsegmentation*)fDetTypeSim->GetSegmentationModel(kSPD);
      //res = (AliITSCalibration*)fDetTypeSim->GetResponseModel(GetITSgeom()->GetStartSPD());
      sim = new AliITSsimulationSPD(fDetTypeSim);
      SetSimulationModel(kSPD,sim);
    }else{ // simulation exists, make sure it is set up properly.
      sim->SetCalibrationModel(GetITSgeom()->GetStartSPD(),(AliITSCalibration*)fDetTypeSim->GetCalibrationModel(GetITSgeom()->GetStartSPD()));
      sim->SetSegmentationModel(kSPD,(AliITSsegmentation*)fDetTypeSim->GetSegmentationModel(kSPD));
      sim->Init();
    } // end if
  } // end if iDetType
  
    /*
      if(fDetTypeSim){
        sim = fDetTypeSim->GetSimulationModel(kSDD);
        if (!sim) {
            seg = (AliITSsegmentation*)fDetTypeSim->GetSegmentationModel(kSDD);
            res = (AliITSCalibration*)fDetTypeSim->GetResponseModel(GetITSgeom()->GetStartSDD());
            sim = new AliITSsimulationSDD(seg,res);
            SetSimulationModel(kSDD,sim);
        }else{ // simulation exists, make sure it is set up properly.
	  sim->SetResponseModel((AliITSCalibration*)fDetTypeSim->GetResponseModel(GetITSgeom()->GetStartSDD()));
	  sim->SetSegmentationModel((AliITSsegmentation*)fDetTypeSim->GetSegmentationModel(kSDD));
	  sim->Init();
        } //end if
    } // end if iDetType
    if(fDetTypeSim){
        sim = fDetTypeSim->GetSimulationModel(kSSD);
        if (!sim) {
            seg = (AliITSsegmentation*)fDetTypeSim->GetSegmentationModel(kSSD);
            res = (AliITSCalibration*)fDetTypeSim->GetResponseModel(GetITSgeom()->GetStartSSD());
            sim = new AliITSsimulationSSD(seg,res);
            SetSimulationModel(kSSD,sim);
        }else{ // simulation exists, make sure it is set up properly.
	  sim->SetResponseModel((AliITSCalibration*)fDetTypeSim->GetResponseModel(GetITSgeom()->GetStartSSD()));
	  sim->SetSegmentationModel((AliITSsegmentation*)fDetTypeSim->GetSegmentationModel(kSSD));
	  sim->Init();
        } // end if
    } //
    */ 
}
//______________________________________________________________________
void AliITSvSPD02::DrawModule() const {
    ////////////////////////////////////////////////////////////////////////
    //     Draw a shaded view of the ITS SPD test beam version 1.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    ////////////////////////////////////////////////////////////////////////
    // Set everything unseen
    gMC->Gsatt("*", "seen", -1);
    // Set ALIC mother visible
    gMC->Gsatt("ALIC","SEEN",0);
    // Set ALIC ITS visible
    gMC->Gsatt("ITSV","SEEN",0);
    // Set ALIC Telescopes visible
    gMC->Gsatt("ITEL","SEEN",0);
    // Set ALIC detetcor visible
    gMC->Gsatt("IDET","SEEN",0);
    // Set Detector chip mother visible and drawn
    gMC->Gsatt("IPC0","SEEN",1);
    // Set Detector mother visible and drawn
    gMC->Gsatt("ITS0","SEEN",1);
    // Set minibus chip mother visible and drawn
    gMC->Gsatt("ICMB","SEEN",1);
    // Set minibus mother visible and drawn
    gMC->Gsatt("IMB0","SEEN",1);
}
//______________________________________________________________________
void AliITSvSPD02::StepManager(){
    ////////////////////////////////////////////////////////////////////////
    //    Called for every step in the ITS SPD test beam, then calles the 
    // AliITShit class  creator with the information to be recoreded about
    //  that hit.
    //     The value of the macro ALIITSPRINTGEOM if set to 1 will allow the
    // printing of information to a file which can be used to create a .det
    // file read in by the routine CreateGeometry(). If set to 0 or any other
    // value except 1, the default behavior, then no such file is created nor
    // it the extra variables and the like used in the printing allocated.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    ////////////////////////////////////////////////////////////////////////
    Int_t         copy=0, id;
    TLorentzVector position, momentum;
    static TLorentzVector position0;
    static Int_t stat0=0;
    if((id=gMC->CurrentVolID(copy) == fIDMother)&&
       (gMC->IsTrackEntering()||gMC->IsTrackExiting())){
        copy = fTrackReferences->GetEntriesFast();
        TClonesArray &lTR = *fTrackReferences;
        // Fill TrackReference structure with this new TrackReference.
        new(lTR[copy]) AliTrackReference(gAlice->GetMCApp()->
                                         GetCurrentTrackNumber());
    } // if Outer ITS mother Volume
    if(!(this->IsActive())){
        return;
    } // end if !Active volume.
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
    id = gMC->CurrentVolID(copy);
    if(id==fIdSens[0]){  // Volume name "IMBS"
        vol[2] = vol[1] = 1; // Det, ladder
        id = gMC->CurrentVolOffID(2,copy);
        //detector copy in the ladder = 1<->4  (ITS1 < I101 < I103 < I10A)
        vol[0] = copy; // Lay
        if(copy>2) vol[0]++;
    } else if(id == fIdSens[1]){ // Volume name "ITST"
        vol[0] = 3; // layer
        vol[1] = 1; // ladder
        id = gMC->CurrentVolOffID(2,copy);
        //detector copy in the ladder = 1<->4  (ITS2 < I1D1 < I1D3 < I20A)
        vol[2] = 1;  // detector
    } else return; // end if
    //
    gMC->TrackPosition(position);
    gMC->TrackMomentum(momentum);
    vol[4] = stat0;
    if(gMC->IsTrackEntering()){
        position0 = position;
        stat0 = vol[3];
        return;
    } // end if IsEntering
    // Fill hit structure with this new hit only for non-entrerance hits.
    else new(lhits[fNhits++]) AliITShit(fIshunt,
                                  gAlice->GetMCApp()->GetCurrentTrackNumber(),
                                        vol,gMC->Edep(),gMC->TrackTime(),
                                        position,position0,momentum);
    //
    position0 = position;
    stat0 = vol[3];

    return;
}

