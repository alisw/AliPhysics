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

#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <TMath.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TTUBE.h>
#include <TTUBS.h>
#include <TPCON.h>
#include <TFile.h>    // only required for Tracking function?
#include <TCanvas.h>
#include <TObjArray.h>
#include <TLorentzVector.h>
#include <TObjString.h>
#include <TClonesArray.h>
#include <TBRIK.h>
#include <TSystem.h>
#include <TVirtualMC.h>

#include "AliMC.h"
#include "AliRun.h"
#include "AliMagF.h"
#include "AliConst.h"
#include "AliITSGeant3Geometry.h"
#include "AliTrackReference.h"
#include "AliITShit.h"
#include "AliITS.h"
#include "AliITSvSDD03.h"
#include "AliITSgeom.h"
#include "AliITSgeomSPD.h"
#include "AliITSgeomSDD.h"
#include "AliITSgeomSSD.h"
#include "AliITSDetType.h"
#include "AliITSresponseSPD.h"
#include "AliITSresponseSDD.h"
#include "AliITSresponseSSD.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSsegmentationSSD.h"
#include "AliITSsimulationSPD.h"
#include "AliITSsimulationSDD.h"
#include "AliITSsimulationSSD.h"
#include "AliITSClusterFinderSPD.h"
#include "AliITSClusterFinderSDD.h"
#include "AliITSClusterFinderSSD.h"


ClassImp(AliITSvSDD03)

//______________________________________________________________________
AliITSvSDD03::AliITSvSDD03() {
    ////////////////////////////////////////////////////////////////////////
    // Standard default constructor for the ITS SDD test beam 2002 version 1.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    A default created class.
    ////////////////////////////////////////////////////////////////////////
    Int_t i;

    fIdN          = 0;
    fIdName       = 0;
    fIdSens       = 0;
    fEuclidOut    = kFALSE; // Don't write Euclide file
    fGeomDetOut   = kFALSE; // Don't write .det file
    fGeomDetIn    = kFALSE; // Don't Read .det file
    fMajorVersion = IsVersion();
    fMinorVersion = -1;
    for(i=0;i<60;i++) fRead[i] = '\0';
    for(i=0;i<60;i++) fWrite[i] = '\0';
    for(i=0;i<60;i++) fEuclidGeomDet[i] = '\0';
}
//______________________________________________________________________
AliITSvSDD03::AliITSvSDD03(const char *title,Int_t year) : AliITS("ITS", title){
    ////////////////////////////////////////////////////////////////////////
    //    Standard constructor for the ITS SDD testbeam 2002 version 1.
    // Inputs:
    //    const char *title    title for this ITS geometry.
    // Outputs:
    //    none.
    // Return:
    //    A standard created class.
    ////////////////////////////////////////////////////////////////////////
    Int_t i;

    fIdN = 2;
    fIdName = new TString[fIdN];
    fIdName[0] = "IMBS";
    fIdName[1] = "ITST";
    fIdSens    = new Int_t[fIdN];
    for(i=0;i<fIdN;i++) fIdSens[i] = 0;
    fMajorVersion = IsVersion();
    fMinorVersion = 2;
    fEuclidOut    = kFALSE; // Don't write Euclide file
    fGeomDetOut   = kFALSE; // Don't write .det file
    fGeomDetIn    = kFALSE; // Don't Read .det file
    fYear         = year;
    SetThicknessDet1();
    SetThicknessDet2();
    SetThicknessChip1();
    SetThicknessChip2();	 	 	 

    fEuclidGeometry="$ALICE_ROOT/ITS/ITSgeometry_vSDD032.euc";
    strncpy(fEuclidGeomDet,"$ALICE_ROOT/ITS/ITSgeometry_vSDD032.det",60);
    strncpy(fRead,fEuclidGeomDet,60);
    strncpy(fWrite,fEuclidGeomDet,60);
}
//______________________________________________________________________
AliITSvSDD03::AliITSvSDD03(const AliITSvSDD03 &source) : AliITS(source){
    ////////////////////////////////////////////////////////////////////////
    //     Copy Constructor for ITS SDD test beam 2002 version 1.
    // This class is not to be copied. Function only dummy.
    // Inputs:
    //    const AliITSvSDD03 &source   The class to be copied
    // Outputs:
    //    none.
    // Return:
    //    A warning message.
    ////////////////////////////////////////////////////////////////////////
    if(&source == this) return;
    Warning("Copy Constructor","Not allowed to copy AliITSvSDD03");
    return;
}
//______________________________________________________________________
AliITSvSDD03& AliITSvSDD03::operator=(const AliITSvSDD03 &source){
    ////////////////////////////////////////////////////////////////////////
    //    Assignment operator for the ITS SDD test beam 2002 version 1.
    // This class is not to be copied. Function only dummy.
    // Inputs:
    //    const AliITSvSDD03 &source   The class to be copied
    // Outputs:
    //    none.
    // Return:
    //    A Warning message
    ////////////////////////////////////////////////////////////////////////
    if(&source == this) return *this;
    Warning("= operator","Not allowed to copy AliITSvSDD03");
    return *this;
}
//______________________________________________________________________
AliITSvSDD03::~AliITSvSDD03() {
    ////////////////////////////////////////////////////////////////////////
    //    Standard destructor for the ITS SDD test beam 2002 version 1.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    ////////////////////////////////////////////////////////////////////////
}
//______________________________________________________________________
void AliITSvSDD03::BuildGeometry(){
    ////////////////////////////////////////////////////////////////////////
    //    Geometry builder for the ITS SDD test beam 2002 version 1.
    //    ALIC    ALICE Mother Volume
    //     |- ITSV     ITS Mother Volume
    //         |- IDET       Detector under Test
    //         |   |- ITS0       SDD Si Chip
    //         |   |  |- ITST      SDD Sensitivve Volume
    //         |   |- IPC0 *5    Readout chip
    //         |- ITEL *4    SDD Telescope
    //             |- IMB0       SDD Si Chip
    //             |   |- IMBS     SDD Sensitive volume
    //             |- ICMB       Chip MiniBus.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    ////////////////////////////////////////////////////////////////////////
    // Get the top alice volume.
    TNode *ALIC = gAlice->GetGeometry()->GetNode("alice");
    ALIC->cd();

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
    TBRIK *ITSVshape = new TBRIK("ITSVshape","ITS Logical Mother Volume","Air",
				 data[0],data[1],data[2]);
    TNode *ITSV = new TNode("ITSV","ITS Mother Volume",ITSVshape,
			    0.0,0.0,0.0,0,0);
    ITSV->cd(); // set ourselve into ITSV subvolume of ALIC

    // SDD part of telescope (MiniBuS)
    data[0] = 0.705;
    data[1] = 0.5*ddettelescope;
    data[2] = 3.536;
    TBRIK *IMB0shape = new TBRIK("IMB0shape","SDD wafer","Si",
				 data[0],data[1],data[2]);
    Float_t detMiniBusX,detMiniBusY,detMiniBusZ;
    data[0] = detMiniBusX = 0.64;
    data[1] = detMiniBusY = 0.5*ddettelescope;
    data[2] = detMiniBusZ = 3.48;
    TBRIK *IMBSshape = new TBRIK("IMBSshape","SDD Sensitive volume","Si",
				 data[0],data[1],data[2]);
    Float_t chipMiniBusX,chipMiniBusY,chipMiniBusZ;
    data[0] = chipMiniBusX = 0.793;
    data[1] = chipMiniBusY = 0.5*dchipMiniBus;
    data[2] = chipMiniBusZ = 0.68;
    TBRIK *ICMBshape = new TBRIK("ICMBshape","chip Minibus","Si",
				 data[0],data[1],data[2]);
    data[0] = TMath::Max(detMiniBusX,chipMiniBusX);
    data[1] = detMiniBusY+chipMiniBusY;
    data[2] = TMath::Max(detMiniBusZ,chipMiniBusZ);
    TBRIK *ITELshape = new TBRIK("ITELshape","ITELshape","Air",
				 data[0],data[1],data[2]);

    // SDD under test
    Float_t spdX,spdY,spdZ,spdchipX,spdchipY,spdchipZ;
    data[0] = 0.705;
    data[1] = ddettest;
    data[2] = 3.536;
    TBRIK *ITS0shape = new TBRIK("ITS0shape","SDD wafer","Si",
				 data[0],data[1],data[2]); // contains detector
    data[0] = spdX = 0.64;
    data[1] = spdY = ddettest;
    data[2] = spdZ = 3.48;
    TBRIK *ITSTshape = new TBRIK("ITSTshape","SDD sensitive volume","Si",
				 data[0],data[1],data[2]);
    // ITS0 with no translation and unit rotation matrix.
    data[0] = spdchipX = 0.793;
    data[1] = spdchipY = dchiptest;
    data[2] = spdchipZ = 0.68;
    TBRIK *IPC0shape = new TBRIK("IPC0shape","Readout Chips","Si",
				 data[0],data[1],data[2]); // chip under test
    data[0] = TMath::Max(spdchipX,spdX);
    data[1] = spdY+spdchipY;
    data[2] = TMath::Max(spdchipZ,spdZ);
    TBRIK *IDETshape = new TBRIK("IDETshape","Detector Under Test","Air",
				 data[0],data[1],data[2]);
    // Place volumes in geometry
    Int_t i,j;
    char name[20],title[50];
    Double_t px=0.0,py=0.0,pz[4]={-38.0,0.0,0.0,0.0};
    pz[1] = pz[0]+2.0;
    pz[2] = pz[1]+38.0+spdY+spdchipY+34.5;
    pz[3] = pz[2]+2.0;
    TNode *ITEL[4],*ICMB[4],*IMB0[4],*IMBS[4];
    TNode *IDET = new TNode("IDET","Detector Under Test",IDETshape,
			    0.0,0.0,pz[1]+38.0,r0,0);
    IDET->cd();
    TNode *ITS0 = new TNode("ITS0","SDD Chip",ITS0shape,
			    0.0,IDETshape->GetDy()-spdY,0.0,0,0);
    TNode *IPC0[5];
    for(i=0;i<5;i++) { //place readout chips on the back of SDD chip under test
	sprintf(name,"IPC0%d",i);
	sprintf(title,"Readout chip #%d",i+1);
	j = i-2;
	IPC0[i] = new TNode(name,title,IPC0shape,
			    0.0,spdchipY-IDETshape->GetDy(),
			    j*2.0*spdchipZ+j*0.25*(spdZ-5.*spdchipZ),0,0);
    } // end for i
    ITS0->cd();
    TNode *ITST = new TNode("ITST","SDD sensitive volume",ITSTshape,
			    0.0,0.0,0.0,0,0);
    for(Int_t i=0;i<4;i++){
	ITSV->cd();
	sprintf(name,"ITEL%d",i);
	sprintf(title,"Test beam telescope element #%d",i+1);
	ITEL[i] = new TNode(name,title,ITELshape,px,py,pz[i],r0,0);
	ITEL[i]->cd();
	ICMB[i] = new TNode("ICMB","Chip MiniBus",ICMBshape,
			    0.0,-ITELshape->GetDy()+detMiniBusY,0.0,0,0);
	IMB0[i] = new TNode("IMB0","Chip MiniBus",IMB0shape,
			    0.0, ITELshape->GetDy()-detMiniBusY,0.0,0,0);
	IMB0[i]->cd();
	IMBS[i] = new TNode("IMBS","IMBS",IMBSshape,0.0,0.0,0.0,0,0);
	// place IMBS inside IMB0 with no translation and unit rotation matrix.
    } // end for i
    ALIC->cd();
    ITST->SetLineColor(kYellow);
    fNodes->Add(ITST);
    for(i=0;i<4;i++){
	IMBS[i]->SetLineColor(kGreen);
	fNodes->Add(IMBS[i]);
    } // end for i
}
//______________________________________________________________________
void AliITSvSDD03::CreateGeometry(){
    ////////////////////////////////////////////////////////////////////////
    //  This routine defines and Creates the geometry for version 1 of the ITS.
    //    ALIC    ALICE Mother Volume
    //     |- ITSV     ITS Mother Volume
    //         |- IDET       Detector under Test (box containing SDD)
    //         |   |-IDAI        Air inside box
    //         |       |- ITS0       SDD Si Chip
    //         |          |- ITST      SDD Sensitivve Volume
    //         |- ITEL *10   SSD Telescope (plastic box containting SSD's)
    //         |   |- ITAI       Air inside box
    //         |       |- IMB0       SDD Si Chip
    //         |           |- IMBS     SDD Sensitive volume
    //         |-ISNT*4    Sintilator triggers
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
    //Float_t ddettest=200.0E-4,ddettelescope=300.0E-4;
    //Float_t dchipMiniBus=750.0E-4,dchiptest=300.0E-4;
    //Float_t yposition= 0.0;
    const Float_t kmm=0.1,kcm=1.0,kmicm=0.001;

    // Define Rotation-reflextion Matrixes needed
    // 0 is the unit matrix
    AliMatrix(idrotm[0], 90.0,0.0, 0.0,0.0, 90.0,270.0); // SDD and SSD X
    AliMatrix(idrotm[1], 90.0,0.0, 0.0,0.0, 90.0,270.0); // SSD Y
    data[0] = 100.0*kmm;
    data[1] = 100.0*kmm;
    data[2] = 800.0*kmm;
    gMC->Gsvolu("ITSV","BOX ",idtmed[0],data,3);
    gMC->Gspos("ITSV",1,"ALIC",0.0,0.0,0.0,0,"ONLY");

    //cout << "idtmed[0]=" << idtmed[0]<<endl;
    //cout << "idtmed[1]=" << idtmed[1]<<endl;
    // Crossed sintilator triggers (2 in front 2 in back)
    AliMatrix(idrotm[2],90.0,0.0,90.0,90.0,90.0,0.0);//Rotate about Z 90 degree
    data[0] = 10.0*kcm;
    data[1] = 2.0*kcm;
    data[2] = 2.0*kmm;
    gMC->Gsvolu("ISNT","BOX ",idtmed[3],data,3);
    gMC->Gspos("ISNT",1,"ITSV",0.0,0.0,800.0*kmm+data[2],0,"ONLY");
    gMC->Gspos("ISNT",2,"ITSV",0.0,0.0,800.0*kmm,idrotm[2],"ONLY");
    gMC->Gspos("ISNT",3,"ITSV",0.0,0.0,-800.0*kmm,0,"ONLY");
    gMC->Gspos("ISNT",4,"ITSV",0.0,0.0,-800.0*kmm-data[2],idrotm[2],"ONLY");
    Float_t detMiniBusX,detMiniBusY,detMiniBusZ;
    // SSD part of telescope (MiniBuS)
    data[0] = detMiniBusX = 10600.0*kmicm;
    data[1] = detMiniBusY = 0.0150*kcm;
    data[2] = detMiniBusZ = 1.1*kcm;
    gMC->Gsvolu("IMB0", "BOX ", idtmed[1], data, 3);   // contains detector
    data[0] = 0.5*384*50*kmicm;
    data[1] = 0.1499*kcm;
    data[2] = 1.0*kcm;
    gMC->Gsvolu("IMBS","BOX ",idtmed[1],data,3); // sensitive detecor volulme
    gMC->Gspos("IMBS",1,"IMB0",0.0,0.0,0.0,0,"ONLY"); // place IMBS inside
    // Box containing SSD's
    data[0] = 4.0*kcm;
    data[1] = 0.500*kcm;
    data[2] = 5.0*kcm;
    gMC->Gsvolu("ITAI","BOX ",idtmed[0],data,3);
    // Plastic box size = insize + thickness.
    data[0] = data[0] + 2.0*kmm;
    data[1] = data[1] + 200.0*kmicm;
    data[2] = data[2] + 2.0*kmm;
    gMC->Gsvolu("ITEL","BOX ",idtmed[4],data,3);
    gMC->Gspos("ITAI",1,"ITEL",0.0,0.0,0.0,0,"ONLY");
    gMC->Gspos("IMB0",1,"ITAI",0.0,data[1]-detMiniBusY,0.0,0,"ONLY");

    // SDD under test
    Float_t sddX,sddY,sddZ;
    data[0] = sddX = 3.62500*kcm;
    data[1] = sddY = 0.01500*kcm;
    data[2] = sddZ = 4.37940*kcm;
    gMC->Gsvolu("ITS0", "BOX ", idtmed[1], data, 3);   // contains detector
    data[0] = 3.50860*kcm;
    data[1] = 0.01499*kcm;
    data[2] = 3.76320*kcm;
    gMC->Gsvolu("ITST","BOX ",idtmed[1],data,3);// sensitive detecor volume
    gMC->Gspos("ITST",1,"ITS0",0.0,0.0,0.0,0,"ONLY"); // place ITST inside
    // Box containing SDD under test
    data[0] = 4.0*kcm;
    data[1] = 0.5*kcm;
    data[2] = 5.0*kcm;
    gMC->Gsvolu("IDAI","BOX ",idtmed[0],data,3);
    data[0] = data[0] + 2.0*kmm;
    data[1] = data[1] + 200.0*kmicm;
    data[2] = data[2] = 2.0*kmm;
    gMC->Gsvolu("IDET","BOX ",idtmed[0],data,3);
    gMC->Gspos("IDAI",1,"IDET",0.0,0.0,0.0,0,"ONLY");
    gMC->Gspos("ITS0",1,"IDAI",0.0,0.0,0.0,0,"ONLY");

    // Positions detectors, Beam Axis Z, X to the right, Y up to the sky.
    Float_t p00X,p00Y,p00Z,p01X,p01Y,p01Z,p10X,p10Y,p10Z,p11X,p11Y,p11Z;
    p00X = 0.0*kcm;
    p00Y = 0.0*kcm;
    p00Z = -694*kmm;
    gMC->Gspos("ITEL",1,"ITSV",p00X,p00Y,p00Z,idrotm[0],"ONLY");//SSD X
    p01X = 0.0*kcm;
    p01Y = 0.0*kcm;
    p01Z = -684*kmm;
    gMC->Gspos("ITEL",2,"ITSV",p01X,p01Y,p01Z,idrotm[1],"ONLY");//SSD Y
    p01X = 0.0*kcm;
    p01Y = 0.0*kcm;
    p01Z = -612*kmm;
    gMC->Gspos("ITEL",3,"ITSV",p01X,p01Y,p01Z,idrotm[0],"ONLY");//SSD X
    p01X = 0.0*kcm;
    p01Y = 0.0*kcm;
    p01Z = -602*kmm;
    Float_t pdetX,pdetY,pdetZ;
    gMC->Gspos("ITEL",4,"ITSV",p01X,p01Y,p01Z,idrotm[1],"ONLY");//SSD Y
    pdetX = 0.0*kcm;
    pdetY = 0.0*kcm;
    pdetZ = 0.0*kcm;
    gMC->Gspos("IDET",1,"ITSV",pdetX,pdetY,pdetZ,idrotm[0],"ONLY");// Detecor
    p10X = 0.0*kcm;
    p10Y = 0.0*kcm;
    p10Z = +450.0*kmm;
    gMC->Gspos("ITEL",5,"ITSV",p10X,p10Y,p10Z,idrotm[0],"ONLY");//SSD X
    p11X = 0.0*kcm;
    p11Y = 0.0*kcm;
    p11Z = +460.0*kcm;
    gMC->Gspos("ITEL",6,"ITSV",p11X,p11Y,p11Z,idrotm[1],"ONLY");//SSD Y
    p11X = 0.0*kcm;
    p11Y = 0.0*kcm;
    p11Z = +540.0*kcm;
    gMC->Gspos("ITEL",7,"ITSV",p11X,p11Y,p11Z,idrotm[0],"ONLY");//SSD X
    p11X = 0.0*kcm;
    p11Y = 0.0*kcm;
    p11Z = +550.0*kcm;
    gMC->Gspos("ITEL",8,"ITSV",p11X,p11Y,p11Z,idrotm[1],"ONLY");//SSD Y
    p11X = 0.0*kcm;
    p11Y = 0.0*kcm;
    p11Z = +737.0*kcm;
    gMC->Gspos("ITEL",9,"ITSV",p11X,p11Y,p11Z,idrotm[0],"ONLY");//SSD X
    p11X = 0.0*kcm;
    p11Y = 0.0*kcm;
    p11Z = +747.0*kcm;
    gMC->Gspos("ITEL",10,"ITSV",p11X,p11Y,p11Z,idrotm[1],"ONLY");//SSD Y
}
//______________________________________________________________________
void AliITSvSDD03::CreateMaterials(){
    ////////////////////////////////////////////////////////////////////////
    //
    // Create ITS SDD test beam materials
    //     This function defines the default materials used in the Geant
    // Monte Carlo simulations for the geometries AliITSv1, AliITSv3,
    // AliITSvSDD03.
    // In general it is automatically replaced by
    // the CreatMaterials routine defined in AliITSv?. Should the function
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
    //
    const Float_t kgpcm3=1.0,kcm=1.0;
    //
    Float_t z[10],a[10],w[10];

    z[0] = 7.0; a[0] = 14.00674; w[0] = 0.80;
    z[1] = 8.0; a[1] = 15.99940; w[1] = 0.20;
    AliMixture(1,"AIR$",a,z,0.12050E-02*kgpcm3,2,w);
    AliMedium(1,"AIR$",1,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,
	      epsilAir,stminAir);

    AliMaterial(2,"SI$",0.28086E+02,0.14000E+02,0.23300E+01,
		0.93600E+01*kcm,0.99900E+03);
    AliMedium(2,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
    // sintilator is basicaly polystyrene
    z[0] = 7.0; a[0] = 14.00674; w[0] = 0.80;
    AliMixture(3,"Sintilator$",a,z,1.032*kgpcm3,0,w);
    AliMedium(3,"Sintilator$",3,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
    // assumed to be Lucite/Plexiglas
    z[0] = 7.0; a[0] = 14.00674; w[0] = 0.80;
    AliMixture(4,"PlasticBox$",a,z,1.18*kgpcm3,0,w);
    AliMedium(4,"PlasticBox$",4,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
}
//______________________________________________________________________
void AliITSvSDD03::InitAliITSgeom(){
    //     Based on the geometry tree defined in Geant 3.21, this
    // routine initilizes the Class AliITSgeom from the Geant 3.21 ITS geometry
    // sturture.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.

    if(strcmp(gMC->GetName(),"TGeant3")) {
        Error("InitAliITSgeom",
              "Wrong Monte Carlo. InitAliITSgeom uses TGeant3 calls");
        return;
    } // end if
    cout << "Reading Geometry transformation directly from Geant 3." << endl;
    const Int_t np=384;
    const Float_t pitch=50.E-4;/*cm*/
    Float_t box[3]={0.5*pitch*(Float_t)np,150.E-4,1.0},p[np],n[np];
    const Int_t ltypess = 2;
    const Int_t nlayers = 6;
    const Int_t ndeep = 6;
    Int_t itsGeomTreeNames[ltypess][ndeep],lnam[20],lnum[20];
    Int_t nlad[nlayers],ndet[nlayers];
    Double_t t[3],r[10];
    Float_t  par[20],att[20];
    Int_t    npar,natt,idshape,imat,imed;
    AliITSGeant3Geometry *ig = new AliITSGeant3Geometry();
    Int_t mod,typ,lay,lad,det,cpy,i,j,k;
    Char_t names[ltypess][ndeep][4];
    Int_t itsGeomTreeCopys[ltypess][ndeep];
    Char_t *namesA[ltypess][ndeep] = {
        {"ALIC","ITSV","ITEL","ITAI","IMB0","IMBS"}, // lay=5
        {"ALIC","ITSV","IDET","IDAI","ITS0","ITST"}};// Test SDD
    Int_t itsGeomTreeCopysA[ltypess][ndeep]= {{1,1,10,1,1,1},// lay=5
                                              {1,1,1,1,1,1}};//lay=3 TestSDD
    for(i=0;i<ltypess;i++)for(j=0;j<ndeep;j++){
        for(k=0;k<4;k++) names[i][j][k] = namesA[i][j][k];
        itsGeomTreeCopys[i][j] = itsGeomTreeCopysA[i][j];
    } // end for i,j
    for(i=0;i<np;i++){// Fill in anode and cathode strip locations (lower edge)
        p[i] = 0.5*pitch*(Float_t)np + pitch*(Float_t)i;
        n[i] = pitch*(Float_t)np - p[i];
    } // end for i
    // Sorry, but this is not very pritty code. It should be replaced
    // at some point with a version that can search through the geometry
    // tree its self.
    cout << "Reading Geometry informaton from Geant3 common blocks" << endl;
    for(i=0;i<20;i++) lnam[i] = lnum[i] = 0;
    for(i=0;i<ltypess;i++)for(j=0;j<ndeep;j++) 
        strncpy((char*) &itsGeomTreeNames[i][j],names[i][j],4);
    //	itsGeomTreeNames[i][j] = ig->StringToInt(names[i][j]);
    mod = 5;
    if(fITSgeom!=0) delete fITSgeom;
    nlad[0]=1;nlad[1]=1;nlad[2]=1;nlad[3]=1;nlad[4]=1;
    ndet[0]=1;ndet[1]=1;ndet[2]=1;ndet[3]=1;ndet[4]=1;
    fITSgeom = new AliITSgeom(0,nlayers,nlad,ndet,mod);
    for(typ=1;typ<=ltypess;typ++){
        for(j=0;j<ndeep;j++) lnam[j] = itsGeomTreeNames[typ-1][j];
        for(j=0;j<ndeep;j++) lnum[j] = itsGeomTreeCopys[typ-1][j];
        lad = 1;
        det = 1;
        for(cpy=1;cpy<=itsGeomTreeCopys[typ-1][2];cpy++){
            lnum[2] = cpy;
            lay = cpy;
            if(cpy>2 && typ==1) lay = cpy +1;
            if(typ==2) lay = 3;
            mod = lay-1;
            ig->GetGeometry(ndeep,lnam,lnum,t,r,idshape,npar,natt,par,att,
                            imat,imed);
            switch (typ){
            case 1:
                if(!(fITSgeom->IsShapeDefined((Int_t)kSDD))){
                    fITSgeom->CreatMatrix(mod,lay,lad,det,kSDD,t,r);
                    fITSgeom->ReSetShape(kSDD,new AliITSgeomSDD256(npar,par));
                } // end if
                break;
            case 2:
                if(!(fITSgeom->IsShapeDefined((Int_t)kSSD))){
                    fITSgeom->CreatMatrix(mod,lay,lad,det,kSSD,t,r);
                    fITSgeom->ReSetShape(kSSD,new AliITSgeomSSD(box,0.0,0.0,
                                                              np+1,p,np+1,n));
                } // end if
                break;
            } // end switch
        } // end for cpy
    } // end for typ
    return;
}
//______________________________________________________________________
void AliITSvSDD03::Init(){
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
    cout << " AliITSvSDD03" << fMinorVersion << "_Init ";
    for(i=0;i<25;i++) cout << "*";cout << endl;
//
    if(fRead[0]=='\0') strncpy(fRead,fEuclidGeomDet,60);
    if(fWrite[0]=='\0') strncpy(fWrite,fEuclidGeomDet,60);
    if(fITSgeom!=0) delete fITSgeom;
    fITSgeom = new AliITSgeom();
    if(fGeomDetIn) fITSgeom->ReadNewFile(fRead);
    if(!fGeomDetIn) this->InitAliITSgeom();
    if(fGeomDetOut) fITSgeom->WriteNewFile(fWrite);
    AliITS::Init();
//
    for(i=0;i<72;i++) cout << "*";
    cout << endl;
    fIDMother = gMC->VolId("ITSV"); // ITS Mother Volume ID.
}
//______________________________________________________________________
void AliITSvSDD03::SetDefaults(){
    // sets the default segmentation, response, digit and raw cluster classes
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    const Float_t kconv = 1.0e+04; // convert cm to microns

    Info("SetDefaults","Setting up only SDD detector");

    AliITSDetType *iDetType;
    AliITSgeomSDD *s1;
    AliITSgeomSSD *s2;
    //Int_t i;
    //Float_t bx[256],bz[280];

    // SDD
    iDetType=DetType(kSDD);
    s1 = (AliITSgeomSDD*) fITSgeom->GetShape(kSDD);// Get shape info. Do it this way for now.
    AliITSresponseSDD *resp1=new AliITSresponseSDD("simulated");
    SetResponseModel(kSDD,resp1);
    AliITSsegmentationSDD *seg1=new AliITSsegmentationSDD(fITSgeom,resp1);
    seg1->SetDetSize(s1->GetDx()*kconv, // base this on AliITSgeomSDD
		     s1->GetDz()*2.*kconv, // for now.
		     s1->GetDy()*2.*kconv); // x,z,y full width in microns.
    seg1->SetNPads(256,256);// Use AliITSgeomSDD for now
    SetSegmentationModel(kSDD,seg1);
    const char *kData1=(iDetType->GetResponseModel())->DataType();
    const char *kopt=iDetType->GetResponseModel()->ZeroSuppOption();
    if((!strstr(kopt,"2D")) && (!strstr(kopt,"1D")) || strstr(kData1,"real") ){
	iDetType->ClassNames("AliITSdigit","AliITSRawClusterSDD");
    } else iDetType->ClassNames("AliITSdigitSDD","AliITSRawClusterSDD");
//    SetSimulationModel(kSDD,new AliITSsimulationSDD(seg1,resp1));
//    iDetType->ReconstructionModel(new AliITSClusterFinderSDD());

    // SSD  Layer 5
    iDetType=DetType(kSSD);
    s2 = (AliITSgeomSSD*) fITSgeom->GetShape(kSSD);// Get shape info. Do it this way for now.
    AliITSresponse *resp2=new AliITSresponseSSD("simulated");
    SetResponseModel(kSSD,resp2);
    AliITSsegmentationSSD *seg2=new AliITSsegmentationSSD(fITSgeom);
    seg2->SetDetSize(s2->GetDx()*2.*kconv, // base this on AliITSgeomSSD
		     s2->GetDz()*2.*kconv, // for now.
		     s2->GetDy()*2.*kconv); // x,z,y full width in microns.
    seg2->SetPadSize(95.,0.); // strip x pitch in microns
    seg2->SetNPads(768,0); // number of strips on each side.
    seg2->SetAngles(0.0075,0.0275); // strip angels rad P and N side.
    seg2->SetAnglesLay5(0.0075,0.0275); // strip angels rad P and N side.
    seg2->SetAnglesLay6(0.0275,0.0075); // strip angels rad P and N side.
    SetSegmentationModel(kSSD,seg2); 
    const char *kData2=(iDetType->GetResponseModel())->DataType();
    if(strstr(kData2,"real") ) iDetType->ClassNames("AliITSdigit",
						    "AliITSRawClusterSSD");
    else iDetType->ClassNames("AliITSdigitSSD","AliITSRawClusterSSD");
//    SetSimulationModel(kSSD,new AliITSsimulationSSD(seg2,resp2));
//    iDetType->ReconstructionModel(new AliITSClusterFinderSSD());

    if(kNTYPES>3){
	Warning("SetDefaults",
		"Only the four basic detector types are initialised!");
    }// end if
    return;
}
//______________________________________________________________________
void AliITSvSDD03::SetDefaultSimulation(){
    // sets the default simulation.
    // Inputs:
    //      none.
    // Outputs:
    //      none.
    // Return:
    //      none.

    AliITSDetType *iDetType;
    AliITSsimulation *sim;
    iDetType=DetType(kSPD);
    sim = iDetType->GetSimulationModel();
    if (!sim) {
        AliITSsegmentation *seg0=
            (AliITSsegmentation*)iDetType->GetSegmentationModel();
        AliITSresponse *res0 = (AliITSresponse*)iDetType->GetResponseModel();
        AliITSsimulationSPD *sim0=new AliITSsimulationSPD(seg0,res0);
        SetSimulationModel(kSPD,sim0);
    }else{ // simulation exists, make sure it is set up properly.
        ((AliITSsimulationSPD*)sim)->Init(
            (AliITSsegmentationSPD*) iDetType->GetSegmentationModel(),
            (AliITSresponseSPD*) iDetType->GetResponseModel());
//        if(sim->GetResponseModel()==0) sim->SetResponseModel(
//            (AliITSresponse*)iDetType->GetResponseModel());
//        if(sim->GetSegmentationModel()==0) sim->SetSegmentationModel(
//            (AliITSsegmentation*)iDetType->GetSegmentationModel());
    } // end if
    iDetType=DetType(kSDD);
    sim = iDetType->GetSimulationModel();
    if (!sim) {
        AliITSsegmentation *seg1=
            (AliITSsegmentation*)iDetType->GetSegmentationModel();
        AliITSresponse *res1 = (AliITSresponse*)iDetType->GetResponseModel();
        AliITSsimulationSDD *sim1=new AliITSsimulationSDD(seg1,res1);
        SetSimulationModel(kSDD,sim1);
    }else{ // simulation exists, make sure it is set up properly.
        ((AliITSsimulationSDD*)sim)->Init(
            (AliITSsegmentationSDD*) iDetType->GetSegmentationModel(),
            (AliITSresponseSDD*) iDetType->GetResponseModel());
//        if(sim->GetResponseModel()==0) sim->SetResponseModel(
//            (AliITSresponse*)iDetType->GetResponseModel());
//        if(sim->GetSegmentationModel()==0) sim->SetSegmentationModel(
//            (AliITSsegmentation*)iDetType->GetSegmentationModel());
    } //end if
    iDetType=DetType(kSSD);
    sim = iDetType->GetSimulationModel();
    if (!sim) {
        AliITSsegmentation *seg2=
            (AliITSsegmentation*)iDetType->GetSegmentationModel();
        AliITSresponse *res2 = (AliITSresponse*)iDetType->GetResponseModel();
        AliITSsimulationSSD *sim2=new AliITSsimulationSSD(seg2,res2);
        SetSimulationModel(kSDD,sim2);
    }else{ // simulation exists, make sure it is set up properly.
        ((AliITSsimulationSSD*)sim)->Init(
            (AliITSsegmentationSSD*) iDetType->GetSegmentationModel(),
            (AliITSresponseSSD*) iDetType->GetResponseModel());
//        if(sim->GetResponseModel()==0) sim->SetResponseModel(
//            (AliITSresponse*)iDetType->GetResponseModel());
//        if(sim->GetSegmentationModel()==0) sim->SetSegmentationModel(
//            (AliITSsegmentation*)iDetType->GetSegmentationModel());
    } // end if
}
//______________________________________________________________________
void AliITSvSDD03::DrawModule(){
    ////////////////////////////////////////////////////////////////////////
    //     Draw a shaded view of the ITS SDD test beam version 1.
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
void AliITSvSDD03::StepManager(){
    ////////////////////////////////////////////////////////////////////////
    //    Called for every step in the ITS SDD test beam, then calles the 
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
    Int_t         copy, id;
    TLorentzVector position, momentum;
    static TLorentzVector position0;
    static Int_t stat0=0;
    if((id=gMC->CurrentVolID(copy) == fIDMother)&&
       (gMC->IsTrackEntering()||gMC->IsTrackExiting())){
	copy = fTrackReferences->GetEntriesFast();
	TClonesArray &lTR = *fTrackReferences;
	// Fill TrackReference structure with this new TrackReference.
	new(lTR[copy]) AliTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber());
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
	vol[2] = 1; // detector
	id = gMC->CurrentVolOffID(3,copy);
	//detector copy in the ladder = 1<->4  (ITS1 < I101 < I103 < I10A)
	vol[1] = copy; // ladder
	vol[0] = 5; // Lay
	if(copy>4){
	    vol[0] = 6;
	    vol[1] = copy - 4;
	} // end if
	if(copy>2) vol[0]++;
    } else if(id == fIdSens[1]){ // Volume name "ITST"
	vol[0] = 3; // layer
	vol[1] = 1; // ladder
	id = gMC->CurrentVolOffID(3,copy);
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
					gAlice->GetMCApp()->GetCurrentTrackNumber(),vol,
					gMC->Edep(),gMC->TrackTime(),position,
					position0,momentum);
    //
    position0 = position;
    stat0 = vol[3];

    return;
}

