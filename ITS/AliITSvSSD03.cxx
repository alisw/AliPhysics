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

// Enrico Fragiacomo - 15/03/2004
// Geometry for the June 2003 SSD beam test

#include <TGeometry.h>
#include <TNode.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TBRIK.h>
#include <TVirtualMC.h>
#include <TGeoMatrix.h>

#include "AliRun.h"
#include "AliMagF.h"
#include "AliTrackReference.h"
#include "AliITShit.h"
#include "AliITS.h"
#include "AliITSvSSD03.h"
#include "AliITSgeom.h"
#include "AliITSgeomSSD.h"
#include "AliITSDetTypeSim.h"
#include "AliITSCalibrationSSD.h"
#include "AliITSsegmentationSSD.h"
#include "AliITSsimulationSSD.h"
#include "AliMC.h"


///////////////////////////////////////////////////////////////////////
// Step manager and 
// geometry class
// for the ITS 
// SSD test beam
// geometry of June 2003
// 
///////////////////////////////////////////////////////////////////////
ClassImp(AliITSvSSD03)

//______________________________________________________________________
AliITSvSSD03::AliITSvSSD03() {
    ////////////////////////////////////////////////////////////////////////
    // Standard default constructor for the ITS SSD test beam 2003 version 1.
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
    fGeomNumber   = 2003; // default value
    for(i=0;i<60;i++) fRead[i] = '\0';
    for(i=0;i<60;i++) fWrite[i] = '\0';
    for(i=0;i<60;i++) fEuclidGeomDet[i] = '\0';
}
//______________________________________________________________________
AliITSvSSD03::AliITSvSSD03(const char *title,Int_t gn) : AliITS("ITS", title){
    ////////////////////////////////////////////////////////////////////////
    //    Standard constructor for the ITS SSD testbeam 2003 version 1.
    // Inputs:
    //    const char *title    title for this ITS geometry.
    //    Int_t      gn        Geometry version number (year) default 2003.
    // Outputs:
    //    none.
    // Return:
    //    A standard created class.
    ////////////////////////////////////////////////////////////////////////
    Int_t i;

    fGeomNumber = gn;

    fIdN = 1; 
    fIdName = new TString[fIdN];
    fIdName[0] = "ITST";
    fIdSens    = new Int_t[fIdN];
    for(i=0;i<fIdN;i++) fIdSens[i] = 0;

    fMajorVersion = IsVersion();
    fMinorVersion = 2;
    fEuclidOut    = kFALSE; // Don't write Euclide file
    fGeomDetOut   = kFALSE; // Don't write .det file
    fGeomDetIn    = kFALSE; // Don't Read .det file

    SetThicknessDet1();
    SetThicknessDet2();
    SetThicknessChip1();
    SetThicknessChip2();	 	 	 

    fEuclidGeometry="$ALICE_ROOT/ITS/ITSgeometry_vSSD03.euc";
    strncpy(fEuclidGeomDet,"$ALICE_ROOT/ITS/ITSgeometry_vSSD03.det",60);
    strncpy(fRead,fEuclidGeomDet,60);
    strncpy(fWrite,fEuclidGeomDet,60);
}
//______________________________________________________________________
AliITSvSSD03::AliITSvSSD03(const AliITSvSSD03 &source) : AliITS(source){
    ////////////////////////////////////////////////////////////////////////
    //     Copy Constructor for ITS SSD test beam 2003 version 1.
    // This class is not to be copied. Function only dummy.
    // Inputs:
    //    const AliITSvSSD03 &source   The class to be copied
    // Outputs:
    //    none.
    // Return:
    //    A warning message.
    ////////////////////////////////////////////////////////////////////////
    if(&source == this) return;
    Warning("Copy Constructor","Not allowed to copy AliITSvSSD03");
    return;
}
//______________________________________________________________________
AliITSvSSD03& AliITSvSSD03::operator=(const AliITSvSSD03 &source){
    ////////////////////////////////////////////////////////////////////////
    //    Assignment operator for the ITS SSD test beam 2003 version 1.
    // This class is not to be copied. Function only dummy.
    // Inputs:
    //    const AliITSvSSD03 &source   The class to be copied
    // Outputs:
    //    none.
    // Return:
    //    A Warning message
    ////////////////////////////////////////////////////////////////////////
    if(&source == this) return *this;
    Warning("= operator","Not allowed to copy AliITSvSSD03");
    return *this;
}
//______________________________________________________________________
AliITSvSSD03::~AliITSvSSD03() {
    ////////////////////////////////////////////////////////////////////////
    //    Standard destructor for the ITS SSD test beam 2003 version 1.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    ////////////////////////////////////////////////////////////////////////
}
//______________________________________________________________________
void AliITSvSSD03::BuildGeometry(){
    ////////////////////////////////////////////////////////////////////////
    //    Geometry builder for the ITS SSD test beam 2003 version 1.
    //    ALIC    ALICE Mother Volume
    //     |- ITSV     ITS Mother Volume
    //         |- ITST       Detector under Test
    //
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    ////////////////////////////////////////////////////////////////////////
    // Get the top alice volume.

    switch (fGeomNumber){
    case 2003:
        BuildGeometry2003();
        break;
    default:
        BuildGeometry2003();
        break;
    } // end switch
}
//______________________________________________________________________
void AliITSvSSD03::BuildGeometry2003(){
    ////////////////////////////////////////////////////////////////////////
    //    Geometry builder for the ITS SSD test beam 2003 version 1.
    //    ALIC    ALICE Mother Volume
    //     |- ITSV     ITS Mother Volume
    //         |- ITST       Detector under Test
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
    Float_t ddettest=300.0E-4;
    //Float_t yposition= 0.0;
    TRotMatrix *r0 = new TRotMatrix("ITSidrotm0","ITSidrotm0",
				    90.0,0,0.0,0,90.0,270.0);

    // Mother volume (beam along z)
    data[0] = 10.0;   // in centimeter
    data[1] = 50.0;
    data[2] = 100.0;
    TBRIK *iITSVshape = new TBRIK("ITSVshape",
                                  "ITS Logical Mother Volume","Air",
                                  data[0],data[1],data[2]);
    TNode *iITSV = new TNode("ITSV","ITS Mother Volume",iITSVshape,
			    0.0,0.0,0.0,0,0);
    iITSV->cd(); // set ourselve into ITSV subvolume of aALIC

    // SSD part of telescope  (Note. strips in local xz plan)
    data[0] = 3.5;   // half-length of the SSD module
    data[1] = 0.5*ddettest;    // half-width of the SSD module 
    data[2] = 2.0;   // half-heigth of the SSD module
    TBRIK *iITSTshape = new TBRIK("ITSTshape","SSD sensitive volume","Si",
				 data[0],data[1],data[2]);
    TNode *iITST = new TNode("ITST","SSD sensitive volume",iITSTshape,
			    0.0,0.0,0.0,r0,0);

    aALIC->cd();
    iITST->SetLineColor(kYellow);
    fNodes->Add(iITST);
}
//______________________________________________________________________
void AliITSvSSD03::CreateGeometry(){
    ////////////////////////////////////////////////////////////////////////
    //   Geometry builder for the ITS SSD test beam 2003 version 1.
    //    ALIC    ALICE Mother Volume
    //     |- ITSV     ITS Mother Volume
    //         |- IDET       Detector under Test
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    ////////////////////////////////////////////////////////////////////////

    switch (fGeomNumber){
    case 2003:
        CreateGeometry2003();
        break;
    default:
        CreateGeometry2003();
        break;
    } // end switch
}
//______________________________________________________________________
void AliITSvSSD03::CreateGeometry2003(){
  ////////////////////////////////////////////////////////////////////////
  //
  //    ALIC    ALICE Mother Volume
  //     |- ITSV     Beamtest Mother Volume
  //         |
  //         |- ITSA       Aluminum cover for scintillator
  //         |    |-ITSS    first Trieste trigger plastic scintillator 
  //         |- ITSA       Aluminum cover for scintillator
  //         |    |-ITSS    second Trieste's trigger plastic scintillator
  //         |
  //         |- IGAR       Black box around ITST       
  //         |    |-IAIR    Air inside the black box
  //         |        |-ITST    Detector under Test 
  //         |
  //         |- IFRA       Aluminum cover for scintillator
  //         |    |-IFRS    French plastic scintillator 
  //         |
  //         |- ITSA       Aluminum cover for scintillator
  //         |    |-ITSS    third Trieste's plastic scintillator
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
   //Float_t yposition= 0.0;
  
  if(gMC==0) return;
  // Define Rotation-reflextion Matrixes needed
  // 0 is the unit matrix

  // Beamtest mother volume (air) positioned in ALIC mother volume
  data[0] = 500.0;
  data[1] = 500.0;
  data[2] = 1000.0;
  gMC->Gsvolu("ITSV","BOX",idtmed[0],data,3);
  gMC->Gspos("ITSV",1,"ALIC",0.0,0.0,0.0,0,"ONLY");
  
  // Trieste's plastic scintillators for the trigger (2 at beam enter)
  // ...define them (aluminum cover + scintillator inside)
  // aluminum cover
  data[0] = 30.01; // size+2x50 microns Kapton
  data[1] = 1.01;
  data[2] = 20.01;
  //gMC->Gsvolu("ITSA","BOX ",idtmed[3],data,3);// 
  gMC->Gsvolu("ITSA","BOX ",idtmed[4],data,3);// 
  data[0] = 30.0;
  data[1] = 1.0;
  data[2] = 20.0;
  // plastic scintillator
  gMC->Gsvolu("ITSS","BOX ",idtmed[2],data,3);
  gMC->Gspos("ITSS",1,"ITSA",0.0,0.0,0.0,0,"ONLY"); 
  // ... and place them inside ITSV
  AliMatrix(idrotm[0], 90.0,0.0, 0.0,0.0, 90.0,270.0);
  // first scintillator 
  gMC->Gspos("ITSA",1,"ITSV",0.0,0.0,-282.0,idrotm[0],"ONLY"); 
  // second scintillator 
  gMC->Gspos("ITSA",2,"ITSV",0.0,0.0,-280.0,idrotm[0],"ONLY"); 

  // black kapton box with the SSD sensor inside (width 50 microns)
  data[0] = 20.0;
  data[1] = 20.0;
  data[2] = 20.0;
  gMC->Gsvolu("IGAR","BOX ",idtmed[4],data,3); //
  // air in the black kapton box 
  data[0] = 19.99;
  data[1] = 19.99;
  data[2] = 19.99;
  gMC->Gsvolu("IAIR","BOX ",idtmed[0],data,3); //
  // SSD sensor 
  Float_t ddettest=300.0E-4;
  data[0] = 3.5;
  data[1] = 0.5*ddettest;
  data[2] = 2.0;
  gMC->Gsvolu("ITST","BOX ",idtmed[1],data,3);// sensitive detector volume
  // place ITST inside IAIR (no rotation: it will be rotated with IGAR)
  gMC->Gspos("ITST",1,"IAIR",0.0,0.0,0.0,0,"ONLY"); 
  // place IAIR inside IGAR
  gMC->Gspos("IAIR",1,"IGAR",0.0,0.0,0.0,0,"ONLY"); 
  // place IGAR inside ITSV
  AliMatrix(idrotm[0], 90.0,0.0, 0.0,0.0, 90.0,270.0);
  gMC->Gspos("IGAR",1,"ITSV",0.0,0.0,0.0,idrotm[0],"ONLY"); 
  //gMC->Gspos("IGAR",1,"ITSV",0.0,0.0,0.0,0,"ONLY"); 
 
  // The so called French detector 
  // ...define it (Kapton cover + scintillator inside)
  // Kapton cover
  data[0] = 2.01; // size+2x50 microns Kapton width
  data[1] = 1.01;
  data[2] = 1.01;
  gMC->Gsvolu("IFRA","BOX ",idtmed[4],data,3);// 
  data[0] = 2.0;
  data[1] = 1.0;
  data[2] = 1.0;
  // plastic scintillator
  gMC->Gsvolu("IFRS","BOX ",idtmed[2],data,3);
  gMC->Gspos("IFRS",1,"IFRA",0.0,0.0,0.0,0,"ONLY"); 
  // ... and place it inside ITSV
  AliMatrix(idrotm[0], 90.0,0.0, 0.0,0.0, 90.0,270.0); 
  gMC->Gspos("IFRA",1,"ITSV",0.0,0.0,16.0,idrotm[0],"ONLY"); 

  // An other Trieste's plastic scintillator for the trigger 
  // ...just place an other copy inside ITSV
  AliMatrix(idrotm[0], 90.0,0.0, 0.0,0.0, 90.0,270.0);
  gMC->Gspos("ITSA",3,"ITSV",0.0,0.0,270.0,idrotm[0],"ONLY"); 
 
}

//______________________________________________________________________
void AliITSvSSD03::CreateMaterials(){
    ////////////////////////////////////////////////////////////////////////
    //
    // Create ITS SSD test beam materials
    //     This function defines the default materials used in the Geant
    // Monte Carlo simulations for the geometries AliITSv1, AliITSv3,
    // AliITSvSSD03.
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
    case 2003:
        CreateMaterials2003();
        break;
    default:
        CreateMaterials2003();
        break;
    } // end switch
}
//______________________________________________________________________
void AliITSvSSD03::CreateMaterials2003(){
    ////////////////////////////////////////////////////////////////////////
    //
    // Create ITS SSD test beam materials
    //     This function defines the default materials used in the Geant
    // Monte Carlo simulations for the geometries AliITSv1, AliITSv3,
    // AliITSvSSD03.
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

    Int_t   ifield = gAlice->Field()->Integ();
    Float_t fieldm = gAlice->Field()->Max();   

    // Scintillator CH
    Float_t ascin[2]={1.01,12.01};
    Float_t zscin[2]={1,6};
    Float_t wscin[2]={1,1};
    Float_t denscin=1.03;
    AliMixture( 3, "Scintillator$",ascin,zscin,denscin,-2,wscin);
    AliMedium(3, "Scintillator$", 3, 1, ifield, fieldm, 0.1, .01, 
	      0.1, .0001, 0.0);
    
    // Aluminum
    Float_t tmaxfdAl = 0.1; // Degree
    Float_t stemaxAl = 0.01; // cm
    Float_t deemaxAl = 0.1; // Fraction of particle's energy 0<deemax<=1
    Float_t epsilAl  = 1.0E-4;//
    Float_t stminAl  = 0.0; // cm "Default value used"
    AliMaterial(4,  "Al$", 26.98, 13., 2.7, 8.9, 37.2);
    AliMedium(4,  "Al$",  4, 0, ifield, fieldm, tmaxfdAl, stemaxAl, 
	      deemaxAl, epsilAl, stminAl);

    // Air
    Float_t tmaxfdAir = 0.1; // Degree
    Float_t stemaxAir = .10000E+01; // cm
    Float_t deemaxAir = 0.1; // Fraction of particle's energy 0<deemax<=1
    Float_t epsilAir  = 1.0E-4;//
    Float_t stminAir  = 0.0; // cm "Default value used"
    AliMaterial(1,"AIR$",0.14610E+03,0.73000E+01,0.12050E-03,
		0.30423E+05,0.99900E+03);
    AliMedium(1,"AIR$",1,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,
	      epsilAir,stminAir);
    
    // Silicon
    Float_t tmaxfdSi = 0.1; // Degree
    //Float_t stemaxSi = 0.0075; // cm
    //Float_t deemaxSi = 0.1; // Fraction of particle's energy 0<deemax<=1
    //Float_t stminSi  = 0.0; // cm "Default value used"
    //Float_t tmaxfdSi = 10; // Degree
    Float_t stemaxSi = 0.01; // cm
    Float_t deemaxSi = 0.1; // Fraction of particle's energy 0<deemax<=1
    Float_t epsilSi  = 1.0E-4;//
    //Float_t epsilSi  = 0.003;//
    Float_t stminSi  = 0.003; // cm "Default value used"
    AliMaterial(2,"SSD SI$",0.28086E+02,0.14000E+02,0.23300E+01,
		0.93600E+01,0.99900E+03);
    AliMedium(2,"SSD SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);

    // Kapton
    AliMaterial(5, "Kapton$", 12.011, 6., 1.3, 31.27, 999.);
    AliMedium(5, "Kapton$",    5, 0,ifield,fieldm, 10., .01, .1, .003, .003);
}
//______________________________________________________________________
void AliITSvSSD03::InitAliITSgeom(){
  // sturture.
  // Inputs:
  //    none.
  // Outputs:
  //    none.
  // Return:
  //    none.
    const Int_t knlayers=1;
    const TString kname="ALIC_1/ITSV_1/IGAR_1/IAIR_1/ITST_1";
    const Int_t knlad[knlayers]={knlayers*1},kndet[knlayers]={knlayers*1};
    Int_t npar;
    Float_t par[20];
    Double_t trans[3]={3*0.0},rot[10]={10*0.0};
    TGeoHMatrix materix;

    AliITSgeom* geom = new AliITSgeom(0,knlayers,knlad,kndet,1);
    SetITSgeom(geom);
    npar=3;par[0]=3.5;par[1]=0.5*300.0E-4;par[2]=2.0;
    geom->ReSetShape(kSSD,new AliITSgeomSSD275and75(npar,par));
    gMC->GetTransformation(kname.Data(),materix);
    geom->CreateMatrix(0,1,1,1,kSSD,trans,rot);
    geom->SetTrans(0,materix.GetTranslation());
    geom->SetRotMatrix(0,materix.GetRotationMatrix());
    geom->GetGeomMatrix(0)->SetPath(kname.Data());
    return;
}
//______________________________________________________________________
void AliITSvSSD03::Init(){
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
    cout << " ITSvSSD03" << fMinorVersion << "_Init ";
    for(i=0;i<25;i++) cout << "*";cout << endl;
//
    if(fRead[0]=='\0') strncpy(fRead,fEuclidGeomDet,60);
    if(fWrite[0]=='\0') strncpy(fWrite,fEuclidGeomDet,60);
    if(GetITSgeom()!=0) SetITSgeom(0x0);
    AliITSgeom* geom = new AliITSgeom();
    SetITSgeom(geom);
    if(fGeomDetIn) GetITSgeom()->ReadNewFile(fRead);
    if(!fGeomDetIn) this->InitAliITSgeom();
    if(fGeomDetOut) GetITSgeom()->WriteNewFile(fWrite);
    AliITS::Init();
//
    for(i=0;i<72;i++) cout << "*";
    cout << endl;
    if(gMC) fIDMother = gMC->VolId("ITSV"); // ITS Mother Volume ID.
    else fIDMother = 0;
}
//______________________________________________________________________
void AliITSvSSD03::SetDefaults(){
    // sets the default segmentation, rerponse, digit and raw cluster classes
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    const Float_t kconv = 1.0e+04; // convert cm to microns

    Info("SetDefaults","Setting up only SSD detector");

    if(!fDetTypeSim) fDetTypeSim = new AliITSDetTypeSim();
    fDetTypeSim->SetITSgeom(GetITSgeom());
    AliITSgeomSSD  *s0;
    fDetTypeSim->ResetCalibrationArray();
    fDetTypeSim->ResetSegmentation();
 
    //SSD

    // Get shape info. Do it this way for now.
    s0 = (AliITSgeomSSD*) GetITSgeom()->GetShape(kSSD);
    AliITSCalibration *resp0=new AliITSCalibrationSSD("simulated");
    SetCalibrationModel(GetITSgeom()->GetStartSSD(),resp0);	

    AliITSsegmentationSSD *seg0=new AliITSsegmentationSSD();
    seg0->SetDetSize(s0->GetDx()*2.*kconv, // base this on AliITSgeomSSD
		     s0->GetDz()*2.*kconv, // for now.
		     s0->GetDy()*2.*kconv); // x,z,y full width in microns.
    //seg0->SetNPads(256,160);// Number of Bins in x and z

    SetSegmentationModel(kSSD,seg0);

    // set digit and raw cluster classes to be used
    const char *kData0=(fDetTypeSim->GetCalibrationModel(GetITSgeom()->GetStartSSD()))->DataType();
    if (strstr(kData0,"real")) fDetTypeSim->SetDigitClassName(kSSD,"AliITSdigit");
    else fDetTypeSim->SetDigitClassName(kSSD,"AliITSdigitSSD");
//    SetSimulationModel(kSSD,new AliITSsimulationSSD(seg0,resp0));
//    iDetType->ReconstructionModel(new AliITSClusterFinderSSD());

/*
    SetResponseModel(kSPD,new AliITSCalibrationSPD());
    SetSegmentationModel(kSPD,new AliITSsegmentationSPD());
    fDetTypeSim->SetDigitClassName(kSPD,"AliITSdigitSPD");

    SetResponseModel(kSDD,new AliITSCalibrationSDD());
    SetSegmentationModel(kSDD,new AliITSsegmentationSDD());
    fDetTypeSim->SetDigitClassName(kSDD,"AliITSdigitSDD");
*/


    if(fgkNTYPES>3){
	Warning("SetDefaults",
		"Only the four basic detector types are initialised!");
    }// end if
    return;
}
//______________________________________________________________________
void AliITSvSSD03::SetDefaultSimulation(){
    // sets the default simulation.
    // Inputs:
    //      none.
    // Outputs:
    //      none.
    // Return:
    //      none.

  if(!fDetTypeSim) fDetTypeSim = new AliITSDetTypeSim();

  AliITSsimulation *sim;
  //  AliITSsegmentation *seg;
  // AliITSCalibration *res;
  /*
  //SPD
  if(fDetTypeSim){
    sim = fDetTypeSim->GetSimulationModel(kSPD);
    if (!sim) {
      seg = (AliITSsegmentation*)fDetTypeSim->GetSegmentationModel(kSPD);
      res = (AliITSCalibration*)fDetTypeSim->GetResponseModel(nspd);
      sim = new AliITSsimulationSPDdubna(seg,res,1);
      SetSimulationModel(kSPD,sim);
    }else{ // simulation exists, make sure it is set up properly.
      sim->SetSegmentationModel((AliITSsegmentation*)fDetTypeSim->GetSegmentationModel(kSPD));
      sim->SetResponseModel((AliITSCalibration*)fDetTypeSim->GetResponseModel(nspd));
      ((AliITSsimulation*)sim)->Init();
      //        if(sim->GetResponseModel()==0) sim->SetResponseModel(
      //            (AliITSCalibration*)iDetType->GetResponseModel());
      //        if(sim->GetSegmentationModel()==0) sim->SetSegmentationModel(
      //            (AliITSsegmentation*)iDetType->GetSegmentationModel());
    } // end if
  } // end if !fDetTypeSim

  //SDD
  if(fDetTypeSim){
    sim = fDetTypeSim->GetSimulationModel(kSDD);
    if (!sim) {
      seg = (AliITSsegmentation*)fDetTypeSim->GetSegmentationModel(kSDD);
      res = (AliITSCalibration*)fDetTypeSim->GetResponseModel(nsdd);
      sim = new AliITSsimulationSDD(seg,res);
      SetSimulationModel(kSDD,sim);
    }else{ // simulation exists, make sure it is set up properly.
      sim->SetSegmentationModel((AliITSsegmentation*)fDetTypeSim->GetSegmentationModel(kSDD));
      sim->SetResponseModel((AliITSCalibration*)fDetTypeSim->GetResponseModel(nsdd));

      ((AliITSsimulation*)sim)->Init();
      //        if(sim->GetResponseModel()==0) sim->SetResponseModel(
      //            (AliITSCalibration*)iDetType->GetResponseModel());
      //        if(sim->GetSegmentationModel()==0) sim->SetSegmentationModel(
      //            (AliITSsegmentation*)iDetType->GetSegmentationModel());
    } //end if
  } // end if !iDetType
  */
  //SSD
  if(fDetTypeSim){
    sim = fDetTypeSim->GetSimulationModel(kSSD);
    if (!sim) {
      //  seg = (AliITSsegmentation*)fDetTypeSim->GetSegmentationModel(kSSD);
      // res = (AliITSCalibration*)fDetTypeSim->GetResponseModel(GetITSgeom()->GetStartSSD());
      sim = new AliITSsimulationSSD(fDetTypeSim);
      SetSimulationModel(kSSD,sim);
    }else{ // simulation exists, make sure it is set up properly.
      sim->SetSegmentationModel(kSSD,(AliITSsegmentation*)fDetTypeSim->GetSegmentationModel(kSSD));
      sim->SetCalibrationModel(GetITSgeom()->GetStartSSD(),(AliITSCalibration*)fDetTypeSim->GetCalibrationModel(GetITSgeom()->GetStartSSD()));
      ((AliITSsimulation*)sim)->Init();
      //        if(sim->GetResponseModel()==0) sim->SetResponseModel(
      //            (AliITSCalibration*)iDetType->GetResponseModel());
      //        if(sim->GetSegmentationModel()==0) sim->SetSegmentationModel(
      //            (AliITSsegmentation*)iDetType->GetSegmentationModel());
    } // end if
  } // end if !iDetType
}

//______________________________________________________________________
void AliITSvSSD03::DrawModule() const {
    ////////////////////////////////////////////////////////////////////////
    //     Draw a shaded view of the ITS SSD test beam version 1.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    ////////////////////////////////////////////////////////////////////////

    gMC->Gsatt("*", "seen", -1);
    gMC->Gsatt("ALIC","SEEN",0);
    gMC->Gsatt("ITSV","SEEN",0);
    gMC->Gsatt("ITSA","SEEN",1);
    gMC->Gsatt("ITSS","SEEN",1);
    gMC->Gsatt("IGAR","SEEN",1);
    gMC->Gsatt("IAIR","SEEN",0);
    gMC->Gsatt("ITST","SEEN",1);
    gMC->Gsatt("IFRA","SEEN",1);
    gMC->Gsatt("IFRS","SEEN",1);
}
//______________________________________________________________________
void AliITSvSSD03::StepManager(){
    ////////////////////////////////////////////////////////////////////////
    //    Called for every step in the ITS SSD, then calles the 
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
    if(id==fIdSens[0]){  // Volume name "ITST"
	vol[2] = vol[1] = vol[0] = 1; // Det, ladder
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
    else new(lhits[fNhits++]) AliITShit(fIshunt,gAlice->GetMCApp()->GetCurrentTrackNumber(),vol,
					gMC->Edep(),gMC->TrackTime(),position,
					position0,momentum);
    //cout<<gMC->Edep()<<endl;
    //
    position0 = position;
    stat0 = vol[3];

    return;
}

