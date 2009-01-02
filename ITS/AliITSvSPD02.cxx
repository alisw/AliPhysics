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

#include <TGeoManager.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
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
fMajorVersion((Int_t)kvSPD02),
fMinorVersion(2),
fGeomNumber(2002),
fDet1(300.0),
fDet2(300.0),
fChip1(300.0),
fChip2(300.0),
fIDMother(0),
fIgm(kvSPD02){
    ////////////////////////////////////////////////////////////////////////
    // Standard default constructor for the ITS SPD test beam 2002 version 1.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    A default created class.
    ////////////////////////////////////////////////////////////////////////
}
//______________________________________________________________________
AliITSvSPD02::AliITSvSPD02(const char *title,Int_t gn) : AliITS("ITS", title),
fMajorVersion(1),
fMinorVersion(2),
fGeomNumber(2002),
fDet1(300.0),
fDet2(300.0),
fChip1(300.0),
fChip2(300.0),
fIDMother(0),
fIgm(kvSPD02){
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
    //
    //      ITEL ITEL IDET ITEL ITEL
    // Z->  -38  -36   02  36.5 38.5
    //       |    |     |    |    |
    // cpn1  1    2     1    3    4
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
    // These constant character strings are set by cvs during commit
    // do not change them unless you know what you are doing!
    const Char_t *cvsDate="$Date$";
    const Char_t *cvsRevision="$Revision$";

    if(gMC==0) return;
    // Define Rotation-reflextion Matrixes needed
    // 0 is the unit matrix
    AliMatrix(idrotm[0], 90.0,0.0, 0.0,0.0, 90.0,270.0);
    /*
    data[0] = 10.0;
    data[1] = 50.0;
    data[2] = 100.0;
    gMC->Gsvolu("ITSV","BOX",idtmed[0],data,3);
    gMC->Gspos("ITSV",1,"ALIC",0.0,0.0,0.0,0,"ONLY");
    */
    TGeoVolumeAssembly *itsV = gGeoManager->MakeVolumeAssembly("ITSV");
    const Int_t length=100;
    Char_t vstrng[length];
    if(fIgm.WriteVersionString(vstrng,length,(AliITSVersion_t)IsVersion(),
                               fMinorVersion,cvsDate,cvsRevision))
        itsV->SetTitle(vstrng);
    else Error("CreateGeometry","Error writing/setting version string");
    TGeoVolume *alic = gGeoManager->GetVolume("ALIC");
    if(alic==0) {
        Error("CreateGeometry","alic=0");
        return;
    } // end if
    // See idrotm[199] for angle definitions.
    alic->AddNode(itsV,1,0);

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
void AliITSvSPD02::Init(){
    //     Initialise the ITS after it has been created.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.

    AliDebug(1,Form("Init: Major version %d Minor version %d",fMajorVersion,
                 fMinorVersion));
    //
    UpdateInternalGeometry();
    AliITS::Init();

    //
    fIDMother = gMC->VolId("ITSV"); // ITS Mother Volume ID.

}
/*
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
   

//    SetResponseModel(kSDD,new AliITSCalibrationSDD());
//    SetSegmentationModel(kSDD,new AliITSsegmentationSDD());
//    DetType(kSDD)->ClassNames("AliITSdigitSDD","AliITSRawClusterSDD");

//    SetResponseModel(kSSD,new AliITSCalibrationSSD());
//    SetSegmentationModel(kSSD,new AliITSsegmentationSSD());
//    DetType(kSSD)->ClassNames("AliITSdigitSSD","AliITSRawClusterSSD");

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

//      if(fDetTypeSim){
//        sim = fDetTypeSim->GetSimulationModel(kSDD);
//        if (!sim) {
//            seg = (AliITSsegmentation*)fDetTypeSim->GetSegmentationModel(kSDD);
//            res = (AliITSCalibration*)fDetTypeSim->GetResponseModel(GetITSgeom()->GetStartSDD());
//            sim = new AliITSsimulationSDD(seg,res);
//           SetSimulationModel(kSDD,sim);
//        }else{ // simulation exists, make sure it is set up properly.
//	  sim->SetResponseModel((AliITSCalibration*)fDetTypeSim->GetResponseModel(GetITSgeom()->GetStartSDD()));
//	  sim->SetSegmentationModel((AliITSsegmentation*)fDetTypeSim->GetSegmentationModel(kSDD));
//	  sim->Init();
//        } //end if
//    } // end if iDetType
//    if(fDetTypeSim){
//        sim = fDetTypeSim->GetSimulationModel(kSSD);
//        if (!sim) {
//            seg = (AliITSsegmentation*)fDetTypeSim->GetSegmentationModel(kSSD);
//            res = (AliITSCalibration*)fDetTypeSim->GetResponseModel(GetITSgeom()->GetStartSSD());
//            sim = new AliITSsimulationSSD(seg,res);
//            SetSimulationModel(kSSD,sim);
//        }else{ // simulation exists, make sure it is set up properly.
//	  sim->SetResponseModel((AliITSCalibration*)fDetTypeSim->GetResponseModel(GetITSgeom()->GetStartSSD()));
//	  sim->SetSegmentationModel((AliITSsegmentation*)fDetTypeSim->GetSegmentationModel(kSSD));
//	  sim->Init();
//        } // end if
//    } //
}
*/
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

    if(!(this->IsActive())){
        return;
    } // end if !Active volume.
    if(!(gMC->TrackCharge())) return;

    Int_t cpy0=0,cpy1=0,id,mod,ncpys,status;
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
    id = gMC->CurrentVolID(cpy0);
    if(id==fIdSens[0]){  // Volume name "IMBS". Det=1, ladder=1
        ncpys = 4;
    } else if(id == fIdSens[1]){ // Volume name "ITST"
        ncpys = 1;
    } else return; // end if
    id = gMC->CurrentVolOffID(2,cpy1);
    fIgm.DecodeDetector(mod,ncpys,cpy0,cpy1,0);
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

