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
/////////////////////////////////////////////////////////////////
//  Class for the SDD beam test August2004                     //
//                                                             //
//                                                             //
/////////////////////////////////////////////////////////////////

#include <TGeoManager.h>
#include <TLorentzVector.h>
#include <TVirtualMC.h>
#include <TGeoMatrix.h>

#include "AliMC.h"
#include "AliRun.h"
#include "AliMagF.h"
#include "AliTrackReference.h"
#include "AliITShit.h"
#include "AliITSgeom.h"
#include "AliITSgeomSDD.h"
#include "AliITSgeomSSD.h"
#include "AliITSDetTypeSim.h"
#include "AliITSCalibrationSPD.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSCalibrationSSD.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSsegmentationSSD.h"
#include "AliITSsimulationSPD.h"
#include "AliITSsimulationSDD.h"
#include "AliITSsimulationSSD.h"

#include "AliITSvSDD03.h"

ClassImp(AliITSvSDD03)

//______________________________________________________________________
AliITSvSDD03::AliITSvSDD03() :
AliITS(),
fMajorVersion(IsVersion()),
fMinorVersion(2),
fIDMother(0),
fYear(2003),
fTarg(kNoTarg),
fTargThick(0.0),
fIgm(kvSDD03){
    ////////////////////////////////////////////////////////////////////////
    // Standard default constructor for the ITS SDD test beam 2002 version 1.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    A default created class.
    ////////////////////////////////////////////////////////////////////////

    fIdN          = 0;
    fIdName       = 0;
    fIdSens       = 0;
}
//______________________________________________________________________
AliITSvSDD03::AliITSvSDD03(const char *title,Int_t year):
AliITS("ITS", title),
fMajorVersion(IsVersion()),
fMinorVersion(2),
fIDMother(0),
fYear(year),
fTarg(kNoTarg),
fTargThick(0.0),
fIgm(kvSDD03){
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

    fIdN = 3;
    fIdName = new TString[fIdN];
    fIdName[0] = "IMBS";
    fIdName[1] = "ITST";
    fIdName[2] = "ISNT";
    fIdSens    = new Int_t[fIdN];
    for(i=0;i<fIdN;i++) fIdSens[i] = 0;

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
/*
//______________________________________________________________________
Int_t AliITSvSDD03::DecodeDetector(Int_t id,Int_t cpy,Int_t &lay,
                                   Int_t &det,Int_t &lad) const{
    // Given the Geant id and copy volume number, returns the layer, ladder,
    // and detector number, allong with the module number of the detector
    // involved. Returns -1 and lay=0, lad=0, and det=0 if not a sensitive 
    // volume.
    // Inputs:
    //    Int_t id    Geometry volume id number
    //    Int_t cpy   Geometry copy number
    // Outputs:
    //    Int_t lay   ITS layer number
    //    Int_t lad   ITS ladder number
    //    Int_t det   ITS detector number
    // Return:
    //    Int_t module number.
    Int_t mod;

    lay = 0; lad = 0; det = 0; mod = -1;
    if(id==fIdSens[0]){ // Volume name is IMBS (ITEL)
        lad = 1; det = 1;
        lay = cpy;
        if(cpy>4) lay+=2;
        mod = lay-1;
        return mod;
    }// end if
    if(id==fIdSens[1]){ // Volume name is ITST (IDet)
        lad = 1; det = 1;lay = cpy+4; 
	mod = lay-1;
        return mod;
    }// end if
    return mod;
}
*/
//______________________________________________________________________
void AliITSvSDD03::CreateGeometry(){
    ////////////////////////////////////////////////////////////////////////
    //  This routine defines and Creates the geometry for version 1 of the ITS.
    //    ALIC    ALICE Mother Volume
    //     |- ITSV     ITS Mother Volume
    //         |- IDET *2    Detector under Test (boxcontaining SDD)
    //         |   |-IDAI        Air inside box
    //         |       |- ITS0       SDD Si Chip
    //         |          |- ITST      SDD Sensitivve Volume
    //         |- ITEL *10   SSD Telescope (plastic box containting SSD's)
    //         |   |- ITAI       Air inside box
    //         |       |- IMB0       SDD Si Chip
    //         |           |- IMBS     SDD Sensitive volume
    //         |-ISNT*4    Sintilator triggers
    //
    //      ITEL ITEL ITEL ITEL IDET IDET ITEL ITEL ITEL ITEL ITEL ITEL
    // Z->  -584 -574 -504 -494  000 +052 +601 +610 +684 +694 +877 +887
    //        |    |    |    |    |    |    |    |    |    |    |    |
    // cpn1   1    2    3    4    1    2    5    6    7    8    9   10
    //
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
    const Float_t kmm=0.1,kcm=1.0,kmicm=kmm/1000.;
    // These constant character strings are set by cvs during commit
    // do not change them unless you know what you are doing!
    const Char_t *cvsDate="$Date$";
    const Char_t *cvsRevision="$Revision$";
    // Define Rotation-reflextion Matrixes needed
    // 0 is the unit matrix
    AliMatrix(idrotm[0], 90.0,0.0, 0.0,0.0, 90.0,270.0); // SDD and SSD X
    AliMatrix(idrotm[1], 90.0,90.0, 0.0,180.0, 90.0,270.0); // SSD Y
    AliMatrix(idrotm[2],90.0,90.0,90.0,180.0,0.0,0.0);  //Rotate about Z 90 degree
    /*
    data[0] = 150.0*kmm;
    data[1] = 150.0*kmm;
    data[2] = 1100.0*kmm;
    gMC->Gsvolu("ITSV","BOX ",idtmed[0],data,3);
    gMC->Gspos("ITSV",1,"ALIC",0.0,0.0,0.0,0,"ONLY");
    */
    TGeoVolumeAssembly *itsV = gGeoManager->MakeVolumeAssembly("ITSV");
    const Int_t kLength=100;
    Char_t vstrng[kLength];
    if(fIgm.WriteVersionString(vstrng,kLength,(AliITSVersion_t)IsVersion(),
                               fMinorVersion,cvsDate,cvsRevision))
        itsV->SetTitle(vstrng);
    else Error("CreateGeometry","Error writing/setting version string");
    //printf("Title set to %s\n",vstrng);
    TGeoVolume *alic = gGeoManager->GetVolume("ALIC");
    if(alic==0) {
        Error("CreateGeometry","alic=0");
        return;
    } // end if
    alic->AddNode(itsV,1,0);

    // Crossed sintilator triggers (2 in front 2 in back)
    data[0] = 10.0*kcm;
    data[1] = 2.0*kcm;
    data[2] = 2.0*kmm;
    gMC->Gsvolu("ISNT","BOX ",idtmed[2],data,3);
    gMC->Gspos("ISNT",1,"ITSV",0.0,0.0,-950.0*kmm,0,"ONLY");
    gMC->Gspos("ISNT",2,"ITSV",0.0,0.0,-950.0*kmm-data[2],idrotm[2],"ONLY");
    gMC->Gspos("ISNT",3,"ITSV",0.0,0.0,950.0*kmm+data[2],0,"ONLY");
    gMC->Gspos("ISNT",4,"ITSV",0.0,0.0,950.0*kmm,idrotm[2],"ONLY");


////Create Volumes

    // SSD part of telescope (MiniBuS)
    Float_t detMiniBusX,detMiniBusY,detMiniBusZ;
    data[0] = detMiniBusX = 10600.0*kmicm;
    data[1] = detMiniBusY = 0.150*kmm;
    data[2] = detMiniBusZ = 1.1*kcm;
    gMC->Gsvolu("IMB0", "BOX ", idtmed[1], data, 3);   // contains detector
    data[0] = 0.5*384*50*kmicm;
    data[1] = 0.1499*kmm;
    data[2] = 1.0*kcm;
    gMC->Gsvolu("IMBS","BOX ",idtmed[1],data,3); // sensitive detector volume
    gMC->Gspos("IMBS",1,"IMB0",0.0,0.0,0.0,0,"ONLY"); // place IMBS inside
    // Box containing SSD's
    data[0] = 11600.0*kmicm;
    data[1] = 0.450*kcm;
    data[2] = 1.16*kcm;
    gMC->Gsvolu("ITAI","BOX ",idtmed[0],data,3);
    // Plastic box size = insize + thickness.
    data[0] = data[0] + 2.0*kmm;
    data[1] = data[1] + 200.0*kmicm;
    data[2] = data[2] + 2.0*kmm;
    gMC->Gsvolu("ITEL","BOX ",idtmed[3],data,3);
    gMC->Gspos("ITAI",1,"ITEL",0.0,0.0,0.0,0,"ONLY");
    gMC->Gspos("IMB0",1,"ITAI",0.0,0.0,0.0,0,"ONLY");

    // SDD under test
    Float_t sddX,sddY,sddZ;
    data[0] = sddX = 3.62500*kcm;
    data[1] = sddY = 0.1500*kmm;
    data[2] = sddZ = 4.37940*kcm;
    gMC->Gsvolu("ITS0", "BOX ", idtmed[1], data, 3);   // contains detector
    data[0] = 3.50860*kcm;
    data[1] = 0.1499*kmm;
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
    data[2] = data[2] + 2.0*kmm;
    gMC->Gsvolu("IDET","BOX ",idtmed[3],data,3);
    gMC->Gspos("IDAI",1,"IDET",0.0,0.0,0.0,0,"ONLY");
    gMC->Gspos("ITS0",1,"IDAI",0.0,0.0,0.0,0,"ONLY");


//// Position detectors, Beam Axis Z, X to the right, Y up to the sky.
    // Upsteram planes of the telescope    
    Float_t p00X,p00Y,p00Z,p01X,p01Y,p01Z,p10X,p10Y,p10Z,p11X,p11Y,p11Z;
    p00X = 0.0*kcm;
    p00Y = 0.0*kcm;
    p00Z = -584*kmm;
    gMC->Gspos("ITEL",1,"ITSV",p00X,p00Y,p00Z,idrotm[0],"ONLY");//SSD X
    p01X = 0.0*kcm;
    p01Y = 0.0*kcm;
    p01Z = -574*kmm;
    gMC->Gspos("ITEL",2,"ITSV",p01X,p01Y,p01Z,idrotm[1],"ONLY");//SSD Y
    p01X = 0.0*kcm;
    p01Y = 0.0*kcm;
    p01Z = -504*kmm;
    gMC->Gspos("ITEL",3,"ITSV",p01X,p01Y,p01Z,idrotm[0],"ONLY");//SSD X
    p01X = 0.0*kcm;
    p01Y = 0.0*kcm;
    p01Z = -494*kmm;
    gMC->Gspos("ITEL",4,"ITSV",p01X,p01Y,p01Z,idrotm[1],"ONLY");//SSD Y

    // Downstream planes of the telescope
    p10X = 0.0*kcm;
    p10Y = 0.0*kcm;
    p10Z = +601.0*kmm; 
    gMC->Gspos("ITEL",5,"ITSV",p10X,p10Y,p10Z,idrotm[0],"ONLY");//SSD X
    p11X = 0.0*kcm;
    p11Y = 0.0*kcm;
    p11Z = +610.0*kmm; //611.0
    gMC->Gspos("ITEL",6,"ITSV",p11X,p11Y,p11Z,idrotm[1],"ONLY");//SSD Y
    p11X = 0.0*kcm;
    p11Y = 0.0*kcm;
    p11Z = +684.0*kmm;
    gMC->Gspos("ITEL",7,"ITSV",p11X,p11Y,p11Z,idrotm[0],"ONLY");//SSD X
    p11X = 0.0*kcm;
    p11Y = 0.0*kcm;
    p11Z = +694.0*kmm;
    gMC->Gspos("ITEL",8,"ITSV",p11X,p11Y,p11Z,idrotm[1],"ONLY");//SSD Y
    p11X = 0.0*kcm;
    p11Y = 0.0*kcm;
    p11Z = +877.0*kmm;
    gMC->Gspos("ITEL",9,"ITSV",p11X,p11Y,p11Z,idrotm[0],"ONLY");//SSD X
    p11X = 0.0*kcm;
    p11Y = 0.0*kcm;
    p11Z = +887.0*kmm;
    gMC->Gspos("ITEL",10,"ITSV",p11X,p11Y,p11Z,idrotm[1],"ONLY");//SSD Y

    // SDDs 
    Float_t pdet1X,pdet1Y,pdet1Z;
    Float_t pdet2X,pdet2Y,pdet2Z;
    pdet1X = 0.0*kcm;
    pdet1Y = 0.0*kcm;
    pdet1Z = 0.0*kcm;
    gMC->Gspos("IDET",1,"ITSV",pdet1X,pdet1Y,pdet1Z,idrotm[0],"ONLY");// Detector1
    pdet2X = 0.0*kcm;
    pdet2Y = 0.0*kcm;
    pdet2Z = 52*kmm; //52
    gMC->Gspos("IDET",2,"ITSV",pdet2X,pdet2Y,pdet2Z,idrotm[0],"ONLY");// Detector2

// Target definition and placement
    if(fTarg){
      data[0] = 30*kmm;
      data[1] = fTargThick*kmm;  // Target thickness
      data[2] = 30*kmm;
      gMC->Gsvolu("ITGT","BOX ",idtmed[fTarg],data,3);

      Float_t a,z,dens,radl,absl;
      Float_t* ubuf=0; Int_t nbuf;
      char* ssss=0;
      gMC->Gfmate(idtmed[fTarg],ssss,a,z,dens,radl,absl,ubuf,nbuf);

      Info("CreateGeometry","Target A=%f,  Z=%f,  dens=%f",a,z,dens);
      Info("Creategeometry","Target thickness=%f mm",fTargThick);

      Float_t ptgtX,ptgtY,ptgtZ;
      ptgtX = 0.0*kcm;
      ptgtY = 0.0*kcm;
      ptgtZ = -50*kmm;
      gMC->Gspos("ITGT",1,"ITSV",ptgtX,ptgtY,ptgtZ,idrotm[0],"ONLY");// Target
    }else{
      Info("CreateGeometry","No target defined");
    }
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
    Float_t stemaxSi = 0.0075; // cm //0.0075
    Float_t deemaxSi = 0.1; // Fraction of particle's energy 0<deemax<=1
    Float_t epsilSi  = 1.0E-4;//
    Float_t stminSi  = 0.0; // cm "Default value used"

    Float_t tmaxfdAir = 0.1; // Degree
    Float_t stemaxAir = .10000E+01; // 1 cm  //cm
    Float_t deemaxAir = 0.1; // Fraction of particle's energy 0<deemax<=1
    Float_t epsilAir  = 1.0E-4;//
    Float_t stminAir  = 0.0; // cm "Default value used"
    Int_t   ifield = gAlice->Field()->Integ();
    Float_t fieldm = gAlice->Field()->Max();
    //

    // AIR
    Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
    Float_t zAir[4]={6.,7.,8.,18.};
    Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
    Float_t dAir = 1.20479E-3;
    // Lucite/Plexiglass
    Float_t aLuc[3] = {1.,12.,16.};
    Float_t zLuc[3] = {1.,6.,8.};
    Float_t wLuc[3] = {8.,5.,2.};
    Float_t dLuc = 1.19;
    // stainless steel
    Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
    Float_t zsteel[4] = { 26.,24.,28.,14. };
    Float_t wsteel[4] = { .715,.18,.1,.005 };
    Float_t dsteel = 7.88;

    AliMixture(1, "AIR$",aAir,zAir,dAir,4,wAir);
    AliMaterial(2,"SI$",28.086,14.0,2.3300,9.3600,999.00);
    AliMixture(3,"Sintilator$",aLuc,zLuc,dLuc,-3,wLuc);
    AliMixture(4,"PlasticBox$",aLuc,zLuc,dLuc,-3,wLuc);
    AliMaterial(5, "IRON$", 55.85, 26., 7.87, 1.76, 999.00);
    AliMaterial(6, "LEAD$", 207.19, 82., 11.35, .56, 999.00);
    AliMixture(7, "STAINLESS STEEL$", asteel, zsteel,dsteel, 4, wsteel);
    AliMaterial(9, "C$", 12.011, 6., 2.265, 18.8, 999.00);
    AliMaterial(10, "Al$", 26.98, 13., 2.70, 8.9, 999.00);
    AliMaterial(11, "Be$", 9.012, 4., 1.848, 35.3, 999.00);
    AliMaterial(12, "Ti$", 47.88, 22., 4.54, 3.56, 999.00);
    AliMaterial(13, "Sn$", 118.69, 50., 7.31, 1.21, 999.00); 
    AliMaterial(14, "Cu$", 63.55, 29., 8.96, 1.43, 999.00);
    AliMaterial(15, "Ge$", 72.59, 32., 5.323, 2.30, 999.00);
    AliMaterial(20, "W$", 183.85, 74., 19.3, 0.35, 999.00);
 
    AliMedium(1,"AIR$",1,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,
	      epsilAir,stminAir);
    AliMedium(2,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
    AliMedium(3,"Scintillator$",3,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
    AliMedium(4,"PlasticBox$",4,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
    AliMedium(5,"IRON$",5,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
    AliMedium(6,"LEAD$",6,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
    AliMedium(7,"StainlessSteel$",7,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);

    AliMedium(9,"C$",9,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);
    AliMedium(10,"Al$",10,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);
    AliMedium(11,"Be$",11,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);
    AliMedium(12,"Ti$",12,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);
    AliMedium(13,"Sn$",13,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);
    AliMedium(14,"Cu$",14,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);
    AliMedium(15,"Ge$",15,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);
    AliMedium(20,"W$",20,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);
    //dummy materials to avoid warning during simulation (galice.cuts)

   AliMedium(21,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(25,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(26,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(27,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(51,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(52,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(53,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(54,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(55,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(56,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(61,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(62,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(63,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(64,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(65,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(68,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(69,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(70,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(71,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(72,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(73,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(74,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(75,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(76,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(77,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(78,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(79,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(80,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(81,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(82,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(83,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(84,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(85,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(90,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(91,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(92,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(93,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(94,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(95,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);
   AliMedium(96,"SI$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,
	      epsilSi,stminSi);

}/*
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
    const Int_t knlayers=12;
    //   const Int_t kndeep=6;
    const Int_t kltypess=2;
    const AliITSDetector kidet[knlayers]={kSSD,kSDD};
    const TString knames[kltypess]={
        "/ALIC_1/ITSV_1/ITEL_%d/ITAI_1/IMB0_1/IMBS_1",
        "/ALIC_1/ITSV_1/IDET_%d/IDAI_1/ITS0_1/ITST_1"};
    const Int_t kitsGeomTreeCopys[kltypess]={10,2};
    const Int_t knp=384;
    const Float_t kpitch=50.E-4;//cm
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
    AliITSgeom *geom = new AliITSgeom(0,knlayers,nlad,ndet,mod);
    if(GetITSgeom()!=0) SetITSgeom(0x0);// delet existing if there.
    SetITSgeom(geom);

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
        mod = DecodeDetector(fIdSens[i],cp0,lay,lad,det);
        path.Form(knames[i].Data(),cp0);
        gMC->GetTransformation(path.Data(),matrix);
        gMC->GetShape(path.Data(),shapeName,shapePar);
        shapeParF.Set(shapePar.GetSize());
        for(j=0;j<shapePar.GetSize();j++)shapeParF[j]=shapePar[j];
        geom->CreateMatrix(mod,lay,lad,det,kidet[i],trans,rot);
        geom->SetTrans(mod,matrix.GetTranslation());
        geom->SetRotMatrix(mod,matrix.GetRotationMatrix());
        geom->GetGeomMatrix(mod)->SetPath(path.Data());
        switch (kidet[i]){
        case kSDD: if(!(GetITSgeom()->IsShapeDefined((Int_t)kSDD))){
            geom->ReSetShape(kSDD,new AliITSgeomSDD256(shapeParF.GetSize(),
                                                       shapeParF.GetArray()));
            isShapeDefined[i]=kTRUE;
        } break;
        case kSSD:if(!(GetITSgeom()->IsShapeDefined((Int_t)kSSD))){
            geom->ReSetShape(kSSD,new AliITSgeomSSD(box,0.0,0.0,
                                                    knp+1,p,knp+1,n));
            isShapeDefined[i]=kTRUE;
        } break;
        default:{} break;
        } // end switch
    } // end for i,cp0
    return;
}*/
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


    Info("Init","**********AliITSvSDD03 %d _Init *************",fMinorVersion);

    AliDebug(1,Form("Init: Major version %d Minor version %d",fMajorVersion,
                 fMinorVersion));
    //
    UpdateInternalGeometry();
    AliITS::Init();
    //
    fIDMother = gMC->VolId("ITSV"); // ITS Mother Volume ID.

}/*
//______________________________________________________________________
void AliITSvSDD03::SetDefaults(){
    // sets the default segmentation, response, digit and raw cluster classes
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.

  //    const Float_t kconv = 1.0e+04; // convert cm to microns

    if(!fDetTypeSim) fDetTypeSim = new AliITSDetTypeSim();
    fDetTypeSim->SetITSgeom(GetITSgeom());
    fDetTypeSim->ResetCalibrationArray();
    fDetTypeSim->ResetSegmentation();
 
    AliITSgeomSDD *s1;
    AliITSgeomSSD *s2;
    SetCalibrationModel(GetITSgeom()->GetStartSPD(),new AliITSCalibrationSPD());
    SetSegmentationModel(kSPD,(AliITSsegmentationSPD*)
			 (GetITSgeom()->GetShape(kSPD)));
    fDetTypeSim->SetDigitClassName(kSPD,"AliITSdigitSPD");

    // SDD
    s1 = (AliITSgeomSDD*) GetITSgeom()->GetShape(kSDD);// Get shape info. Do it this way for now.
    AliITSCalibrationSDD *resp1=new AliITSCalibrationSDD("simulated");
    SetCalibrationModel(GetITSgeom()->GetStartSDD(),resp1);

    AliITSsegmentationSDD *seg1 = (AliITSsegmentationSDD*)
			 (GetITSgeom()->GetShape(kSDD));
    seg1->SetDriftSpeed(AliITSDriftSpeedSDD::DefaultDriftSpeed());
    seg1->SetNPads(256,256);// Use AliITSgeomSDD for now
    SetSegmentationModel(kSDD,seg1);
    const char *kData1=(fDetTypeSim->GetCalibrationModel(GetITSgeom()->GetStartSDD()))->DataType();

    // SSD  Layer 5

    s2 = (AliITSgeomSSD*) GetITSgeom()->GetShape(kSSD);// Get shape info. Do it this way for now.
   
    AliITSCalibration *resp2= new AliITSCalibrationSSD("simulated");
    SetCalibrationModel(GetITSgeom()->GetStartSSD(),resp2);

    AliITSsegmentationSSD *seg2 = (AliITSsegmentationSSD*)
			 (GetITSgeom()->GetShape(kSSD));
    seg2->SetPadSize(50.,0.); // strip x pitch in microns
    seg2->SetNPads(384,0); // number of strips on each side.
    seg2->SetLayer(5);
    seg2->SetAngles(0.,0.); // strip angles rad P and N side.
    seg2->SetAnglesLay5(0.,0.); // strip angles rad P and N side.
    seg2->SetAnglesLay6(0.,0.); // strip angles rad P and N side.

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
void AliITSvSDD03::SetDefaultSimulation(){
    // sets the default simulation.
    // Inputs:
    //      none.
    // Outputs:
    //      none.
    // Return:
    //      none.

  if(!fDetTypeSim) fDetTypeSim = new AliITSDetTypeSim();
  AliITSsimulation *sim;
  //AliITSsegmentation *seg;
  //AliITSCalibration *res;
  //SPD
  if(fDetTypeSim){
    sim = fDetTypeSim->GetSimulationModel(kSPD);
    if (!sim) {
      //seg =(AliITSsegmentation*)fDetTypeSim->GetSegmentationModel(kSPD);
      //if(seg==0) seg = new AliITSsegmentationSPD();
      //res = (AliITSCalibration*)fDetTypeSim->GetResponseModel(GetITSgeom()->GetStartSPD());
      //if(res==0) res = new AliITSCalibrationSPD();
      sim = new AliITSsimulationSPD(fDetTypeSim);
      SetSimulationModel(kSPD,sim);
    }else{ // simulation exists, make sure it is set up properly.
      sim->SetSegmentationModel(kSPD,(AliITSsegmentation*)fDetTypeSim->GetSegmentationModel(kSPD));
      sim->SetCalibrationModel(GetITSgeom()->GetStartSPD(),(AliITSCalibration*)fDetTypeSim->GetCalibrationModel(GetITSgeom()->GetStartSPD()));
      sim->Init();
    } // end if
  } // end if iDetType
  //SDD
  if(fDetTypeSim){
    sim = fDetTypeSim->GetSimulationModel(kSDD);
    if (!sim) {
      //      seg = (AliITSsegmentation*)fDetTypeSim->GetSegmentationModel(kSDD);
      //res = (AliITSCalibration*)fDetTypeSim->GetResponseModel(GetITSgeom()->GetStartSDD());
      sim = new AliITSsimulationSDD(fDetTypeSim);
      SetSimulationModel(kSDD,sim);
    }else{ // simulation exists, make sure it is set up properly.
      sim->SetSegmentationModel(kSDD,(AliITSsegmentation*)fDetTypeSim->GetSegmentationModel(kSDD));
      sim->SetCalibrationModel(GetITSgeom()->GetStartSDD(),(AliITSCalibration*)fDetTypeSim->GetCalibrationModel(GetITSgeom()->GetStartSDD()));
      
      sim->Init();
    } //end if
  } // end if iDetType
  //SSD
  if(fDetTypeSim){
    sim = fDetTypeSim->GetSimulationModel(kSSD);
    if (!sim) {
      //      seg = (AliITSsegmentation*)fDetTypeSim->GetSegmentationModel(kSSD);
      // res = (AliITSCalibration*)fDetTypeSim->GetResponseModel(GetITSgeom()->GetStartSSD());
      sim = new AliITSsimulationSSD(fDetTypeSim);
      SetSimulationModel(kSSD,sim);
    }else{ // simulation exists, make sure it is set up properly.
      sim->SetSegmentationModel(kSSD,(AliITSsegmentation*)fDetTypeSim->GetSegmentationModel(kSSD));
      sim->SetCalibrationModel(GetITSgeom()->GetStartSSD(),(AliITSCalibration*)fDetTypeSim->GetCalibrationModel(GetITSgeom()->GetStartSSD()));
      sim->Init();
    } // end if
  } // end if iDetType
  }*/
//______________________________________________________________________
void AliITSvSDD03::DrawModule() const{
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
    gMC->Gsatt("ITSV","SEEN",1);
    // Set ALIC Telescopes visible
    gMC->Gsatt("ITEL","SEEN",1);
    gMC->Gsatt("ITEL","colo",2);
    // Set ALIC detetcor visible
    gMC->Gsatt("IDET","SEEN",1);
    gMC->Gsatt("IDET","colo",4);
    // Set ALIC Scintillator visible
    gMC->Gsatt("ISNT","SEEN",1);
    gMC->Gsatt("ISNT","colo",3);
    // Set Detector mother visible and drawn
//    gMC->Gsatt("ITS0","SEEN",1);
    // Set minibus mother visible and drawn
//    gMC->Gsatt("IMB0","SEEN",1);

    // Draw
    gMC->Gdraw("alic", 60, 30, 180, 10,10, .12, .12);
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
    if(!(this->IsActive())) return;
    if(!(gMC->TrackCharge())) return;

    Int_t  cpy0,cpy1,ncpys=0,status,id,mod;
    TLorentzVector position, momentum;
    static AliITShit hit;// Saves on calls to construtors
    //TClonesArray &lhits = *(GetDetTypeSim()->GetHits());
    TClonesArray &lhits = *(Hits());
    //
    // Track status
    // Track status
    status = 0;
    if(gMC->IsTrackInside())      status +=  1;
    if(gMC->IsTrackEntering())    status +=  2;
    if(gMC->IsTrackExiting())     status +=  4;
    if(gMC->IsTrackOut())         status +=  8;
    if(gMC->IsTrackDisappeared()) status += 16;
    if(gMC->IsTrackStop())        status += 32;
    if(gMC->IsTrackAlive())       status += 64;
    //
    // Fill hit structure.
    gMC->TrackPosition(position);
    gMC->TrackMomentum(momentum);
    id   = gMC->CurrentVolID(cpy0);
    gMC->CurrentVolOffID(3,cpy1);
    if(id==fIdSens[0])ncpys=10;
    if(id==fIdSens[1])ncpys=2;
    fIgm.DecodeDetector(mod,ncpys,cpy0,cpy1,1);
    //
    // Fill hit structure.
    //
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
    //Info("StepManager","Calling Copy Constructor");
    new(lhits[fNhits++]) AliITShit(hit); // Use Copy Construtor.
    // Save old position... for next hit.
    hit.SetStartPosition(position);
    hit.SetStartTime(gMC->TrackTime());
    hit.SetStartStatus(status);
    return;
}

