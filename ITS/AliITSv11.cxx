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

//************************************************************************
//                 Inner Traking System geometry v11
//
//  Based on ROOT geometrical modeler
//
// B. Nilsen, L. Gaudichet
//************************************************************************


//#include <Riostream.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include "AliITS.h"
#include "AliITSDetTypeSim.h"
#include "AliITSgeom.h"
#include "AliITSgeomSDD.h"
#include "AliITSgeomSPD.h"
#include "AliITSgeomSSD.h"
#include "AliITShit.h"

#include "AliITSCalibrationSDD.h"
#include "AliITSCalibrationSPD.h"
#include "AliITSCalibrationSSD.h"

#include "AliITSsegmentationSDD.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSsegmentationSSD.h"
#include "AliITSv11.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliTrackReference.h"
#include "AliMC.h"

#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoPcon.h>
//#include "AliITSv11GeometrySPD.h"
#include "AliITSv11GeometrySDD.h"
//#include "AliITSv11GeometrySupport.h"



ClassImp(AliITSv11)
 


//______________________________________________________________________
AliITSv11::AliITSv11() : AliITS(),
  fGeomDetOut(kFALSE), // Don't write .det file
  fGeomDetIn(kTRUE),   // Read .det file
  fByThick(kTRUE),
  fMajorVersion(0),
  fMinorVersion(0),
  fSDDgeom(0)
{
  //    Standard default constructor for the ITS version 11.

    fIdN          = 0;
    fIdName       = 0;
    fIdSens       = 0;
    fEuclidOut    = kFALSE; // Don't write Euclide file
    Int_t i;
    for(i=0;i<60;i++) fRead[i] = '\0';
    for(i=0;i<60;i++) fWrite[i] = '\0';
    for(i=0;i<60;i++) fEuclidGeomDet[i] = '\0';
    strncpy(fRead,"$ALICE_ROOT/ITS/ITSgeometry_vPPRasymmFMD.det",60);
}


//______________________________________________________________________
AliITSv11::AliITSv11(const char *name, const char *title)
  : AliITS("ITS", title),
    fGeomDetOut(kFALSE), // Don't write .det file
    fGeomDetIn(kTRUE),   // Read .det file
    fByThick(kTRUE),
    fMajorVersion(0),
    fMinorVersion(0),
    fSDDgeom(0)
{
  //    Standard constructor for the ITS version 11.

  fSDDgeom = new AliITSv11GeometrySDD(0);

  Int_t i;
  fIdN = 6;
  fIdName = new TString[fIdN];
  fIdName[0] = name; // removes warning message
  fIdName[0] = "ITS1";
  fIdName[1] = "ITS2";
  fIdName[2] = fSDDgeom->GetSenstiveVolumeMame();
  fIdName[3] = fSDDgeom->GetSenstiveVolumeMame();
  fIdName[4] = "ITS5";
  fIdName[5] = "ITS6";
  fIdSens    = new Int_t[fIdN];
  for(i=0;i<fIdN;i++) fIdSens[i] = 0;
  fEuclidOut    = kFALSE; // Don't write Euclide file
  //SetDensityServicesByThickness();
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
    fGeomDetOut(kFALSE), // Don't write .det file
    fGeomDetIn(kTRUE),   // Read .det file
    fByThick(kTRUE),
    fMajorVersion(0),
    fMinorVersion(0),
    fSDDgeom(0)
 {
  // Standard default constructor for the ITS version 11.


  //   fSPDgeom = new AliITSv11GeometrySPD(debugSPD);
  fSDDgeom = new AliITSv11GeometrySDD(debugSDD);
  fSDDgeom->SetDebug(debugSDD);
  //   fSupgeom = new AliITSv11GeometrySupport(debugSUP);

  Int_t i;
  fIdN = 6;
  fIdName = new TString[fIdN];
  fIdName[0] = "ITS1";
  fIdName[1] = "ITS2";
  fIdName[2] = fSDDgeom->GetSenstiveVolumeMame();
  fIdName[3] = fSDDgeom->GetSenstiveVolumeMame();
  fIdName[4] = "ITS5";
  fIdName[5] = "ITS6";
  fIdSens    = new Int_t[fIdN];
  for(i=0;i<fIdN;i++) fIdSens[i] = 0;
  fEuclidOut    = kFALSE; // Don't write Euclide file
  //SetDensityServicesByThickness();
  
  fEuclidGeometry="$ALICE_ROOT/ITS/ITSgeometry_vPPRasymm2.euc";
  strncpy(fEuclidGeomDet,"$ALICE_ROOT/ITS/ITSgeometry_vPPRasymm2.det",60);
  strncpy(fRead,fEuclidGeomDet,60);
  strncpy(fWrite,fEuclidGeomDet,60);
  strncpy(fRead,"$ALICE_ROOT/ITS/ITSgeometry_vPPRasymmFMD.det",60);

  debugITS = (debugSPD && debugSSD && debugSUP && debugSDD); //remove temp. warnings
}


//______________________________________________________________________
AliITSv11::AliITSv11(const AliITSv11 &source) :
  AliITS(source),
  fGeomDetOut(kFALSE), // Don't write .det file
  fGeomDetIn(kTRUE),   // Read .det file
  fByThick(kTRUE),
  fMajorVersion(0),
  fMinorVersion(0),
  fSDDgeom(0)
{
    //     Copy Constructor for ITS version 11. This function is not to be
    // used. If any other instance of this function, other than "this" is
    // passed, an error message is returned.
    // Inputs:
    //   const AliITSv11 &source This class

    if(&source == this) return;
    Warning("Copy Constructor","Not allowed to copy AliITSv11");
    return;
}


//______________________________________________________________________
AliITSv11& AliITSv11::operator=(const AliITSv11 
						  &source){
    //     Assignment operator for the ITS version 11. This function is not 
    // to be used. If any other instance of this function, other than "this" 
    // is passed, an error message is returned.
    // Inputs:
    //   const AliITSv11 &source This class

    if(&source == this) return *this;
    Warning("= operator","Not allowed to copy AliITSv11");
    return *this;
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

  TGeoManager *geoManager = gGeoManager;
  TGeoVolume *vALIC = geoManager->GetTopVolume();

  TGeoPcon *sITS = new TGeoPcon("ITS Top Volume",0.0,360.0,2);

  // DefineSection(section number, Z, Rmin, Rmax).
  const Double_t kcm = 1.0;
  sITS->DefineSection(0,-300.0*kcm,0.01*kcm,50.0*kcm);
  sITS->DefineSection(1,+300.0*kcm,0.01*kcm,50.0*kcm);

  TGeoMedium *air = gGeoManager->GetMedium("ITS_ITSair");
  TGeoVolume *vITS = new TGeoVolume("vITS",sITS,air);
  vITS->SetVisibility(kFALSE);
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
  //
  // Create ITS materials
  // Defined media here should correspond to the one defined in galice.cuts
  // File which is red in (AliMC*) fMCApp::Init() { ReadTransPar(); }
  //

//     Int_t   ifield = gAlice->Field()->Integ();
//     Float_t fieldm = gAlice->Field()->Max();

//     Float_t tmaxfd = 0.1; // 1.0; // Degree
//     Float_t stemax = 1.0; // cm
//     Float_t deemax = 0.1; // 30.0; // Fraction of particle's energy 0<deemax<=1
//     Float_t epsil  = 1.0E-4; // 1.0; // cm
//     Float_t stmin  = 0.0; // cm "Default value used"

//     Float_t tmaxfdSi = 0.1; // .10000E+01; // Degree
//     Float_t stemaxSi = 0.0075; //  .10000E+01; // cm
//     Float_t deemaxSi = 0.1; // 0.30000E-02; // Fraction of particle's energy 0<deemax<=1
//     Float_t epsilSi  = 1.0E-4;// .10000E+01;
//     Float_t stminSi  = 0.0; // cm "Default value used"

//     Float_t tmaxfdAir = 0.1; // .10000E+01; // Degree
//     Float_t stemaxAir = .10000E+01; // cm
//     Float_t deemaxAir = 0.1; // 0.30000E-02; // Fraction of particle's energy 0<deemax<=1
//     Float_t epsilAir  = 1.0E-4;// .10000E+01;
//     Float_t stminAir  = 0.0; // cm "Default value used"

//     Float_t tmaxfdServ = 1.0; // 10.0; // Degree
//     Float_t stemaxServ = 1.0; // 0.01; // cm
//     Float_t deemaxServ = 0.5; // 0.1; // Fraction of particle's energy 0<deemax<=1
//     Float_t epsilServ  = 1.0E-3; // 0.003; // cm
//     Float_t stminServ  = 0.0; //0.003; // cm "Default value used"

//     // Freon PerFluorobuthane C4F10 see 
//     // http://st-support-cooling-electronics.web.cern.ch/
//     //        st-support-cooling-electronics/default.htm
//     Float_t afre[2]  = { 12.011,18.9984032 };
//     Float_t zfre[2]  = { 6., 9. };
//     Float_t wfre[2]  = { 4.,10. };
//     Float_t densfre  = 1.52;


//     //CM55J
//     Float_t aCM55J[4]={12.0107,14.0067,15.9994,1.00794};
//     Float_t zCM55J[4]={6.,7.,8.,1.};
//     Float_t wCM55J[4]={0.908508078,0.010387573,0.055957585,0.025146765};
//     Float_t dCM55J = 1.63;

//     //ALCM55J
//     Float_t aALCM55J[5]={12.0107,14.0067,15.9994,1.00794,26.981538};
//     Float_t zALCM55J[5]={6.,7.,8.,1.,13.};
//     Float_t wALCM55J[5]={0.817657902,0.0093488157,0.0503618265,0.0226320885,0.1};
//     Float_t dALCM55J = 1.9866;

//     //Si Chips
//     Float_t aSICHIP[6]={12.0107,14.0067,15.9994,1.00794,28.0855,107.8682};
//     Float_t zSICHIP[6]={6.,7.,8.,1.,14., 47.};
//     Float_t wSICHIP[6]={0.039730642,0.001396798,0.01169634,
// 			0.004367771,0.844665,0.09814344903};
//     Float_t dSICHIP = 2.36436;

//     //Inox
//     Float_t aINOX[9]={12.0107,54.9380, 28.0855,30.9738,32.066,
// 		      58.6928,55.9961,95.94,55.845};
//     Float_t zINOX[9]={6.,25.,14.,15.,16., 28.,24.,42.,26.};
//     Float_t wINOX[9]={0.0003,0.02,0.01,0.00045,0.0003,0.12,0.17,0.025,0.654};
//     Float_t dINOX = 8.03;

//     //SDD HV microcable
//     Float_t aHVm[5]={12.0107,1.00794,14.0067,15.9994,26.981538};
//     Float_t zHVm[5]={6.,1.,7.,8.,13.};
//     Float_t wHVm[5]={0.520088819984,0.01983871336,0.0551367996,0.157399667056, 0.247536};
//     Float_t dHVm = 1.6087;

//     //SDD LV+signal cable
//     Float_t aLVm[5]={12.0107,1.00794,14.0067,15.9994,26.981538};
//     Float_t zLVm[5]={6.,1.,7.,8.,13.};
//     Float_t wLVm[5]={0.21722436468,0.0082859922,0.023028867,0.06574077612, 0.68572};
//     Float_t dLVm = 2.1035;

//     //SDD hybrid microcab
//     Float_t aHLVm[5]={12.0107,1.00794,14.0067,15.9994,26.981538};
//     Float_t zHLVm[5]={6.,1.,7.,8.,13.};
//     Float_t wHLVm[5]={0.24281879711,0.00926228815,0.02574224025,0.07348667449, 0.64869};
//     Float_t dHLVm = 2.0502;

//     //SDD anode microcab
//     Float_t aALVm[5]={12.0107,1.00794,14.0067,15.9994,26.981538};
//     Float_t zALVm[5]={6.,1.,7.,8.,13.};
//     Float_t wALVm[5]={0.392653705471,0.0128595919215,
// 		      0.041626868025,0.118832707289, 0.431909};
//     Float_t dALVm = 2.0502;

//     //X7R capacitors
//     Float_t aX7R[7]={137.327,47.867,15.9994,58.6928,63.5460,118.710,207.2};
//     Float_t zX7R[7]={56.,22.,8.,28.,29.,50.,82.};
//     Float_t wX7R[7]={0.251639432,0.084755042,0.085975822,
// 		     0.038244751,0.009471271,0.321736471,0.2081768};
//     Float_t dX7R = 7.14567;

//     // AIR
//     Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
//     Float_t zAir[4]={6.,7.,8.,18.};
//     Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
//     Float_t dAir = 1.20479E-3;

//     // Water
//     Float_t aWater[2]={1.00794,15.9994};
//     Float_t zWater[2]={1.,8.};
//     Float_t wWater[2]={0.111894,0.888106};
//     Float_t dWater   = 1.0;

//     // CERAMICS
//   //     94.4% Al2O3 , 2.8% SiO2 , 2.3% MnO , 0.5% Cr2O3
//     Float_t acer[5]  = { 26.981539,15.9994,28.0855,54.93805,51.9961 };
//     Float_t zcer[5]  = {       13.,     8.,    14.,     25.,    24. };
//     Float_t wcer[5]  = {.4443408,.5213375,.0130872,.0178135,.003421};
//     Float_t denscer  = 3.6;

//     //G10FR4
//     Float_t zG10FR4[14] = {14.00,	20.00,	13.00,	12.00,	5.00,
// 			   22.00,	11.00,	19.00,	26.00,	9.00,
// 			   8.00,	6.00,	7.00,	1.00};
//     Float_t aG10FR4[14] = {28.0855000,40.0780000,26.9815380,24.3050000,
// 			   10.8110000,47.8670000,22.9897700,39.0983000,
// 			   55.8450000,18.9984000,15.9994000,12.0107000,
// 			   14.0067000,1.0079400};
//     Float_t wG10FR4[14] = {0.15144894,0.08147477,0.04128158,0.00904554,
// 			   0.01397570,0.00287685,0.00445114,0.00498089,
// 			   0.00209828,0.00420000,0.36043788,0.27529426,
// 			   0.01415852,0.03427566};
//     Float_t densG10FR4= 1.8;
    
//      //--- EPOXY  --- C18 H19 O3
//       Float_t aEpoxy[3] = {15.9994, 1.00794, 12.0107} ; 
//       Float_t zEpoxy[3] = {     8.,      1.,      6.} ; 
//       Float_t wEpoxy[3] = {     3.,     19.,     18.} ; 
//       Float_t dEpoxy = 1.8 ;

//       // rohacell: C9 H13 N1 O2
//     Float_t arohac[4] = {12.01,  1.01, 14.010, 16.};
//     Float_t zrohac[4] = { 6.,    1.,    7.,     8.};
//     Float_t wrohac[4] = { 9.,   13.,    1.,     2.};
//     Float_t drohac    = 0.05;

//     // If he/she means stainless steel (inox) + Aluminium and Zeff=15.3383 then
// //
// // %Al=81.6164 %inox=100-%Al

//     Float_t aInAl[5] = {27., 55.847,51.9961,58.6934,28.0855 };
//     Float_t zInAl[5] = {13., 26.,24.,28.,14. };
//     Float_t wInAl[5] = {.816164, .131443,.0330906,.0183836,.000919182};
//     Float_t dInAl    = 3.075;

//     // Kapton
//     Float_t aKapton[4]={1.00794,12.0107, 14.010,15.9994};
//     Float_t zKapton[4]={1.,6.,7.,8.};
//     Float_t wKapton[4]={0.026362,0.69113,0.07327,0.209235};
//     Float_t dKapton   = 1.42;

//     //SDD ruby sph.
//     Float_t aAlOxide[2]  = { 26.981539,15.9994};
//     Float_t zAlOxide[2]  = {       13.,     8.};
//     Float_t wAlOxide[2]  = {0.4707, 0.5293};
//     Float_t dAlOxide     = 3.97;

//     //---------
//     AliMaterial(1,"ITSsddSi",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
//     AliMedium(1,"ITSsddSi",1,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);
    
//     AliMixture(5,"ITSair",aAir,zAir,dAir,4,wAir);
//     AliMedium(5,"ITSair",5,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);
    
//     AliMixture(7,"ITSsddSiChip",aSICHIP,zSICHIP,dSICHIP,6,wSICHIP);
//     AliMedium(7,"ITSsddSiChip",7,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

//     AliMaterial(79,"SDD SI insensitive$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
//     AliMedium(79,"SDD SI insensitive$",79,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMaterial(11,"ITSal",0.26982E+02,0.13000E+02,0.26989E+01,0.89000E+01,0.99900E+03);
//     AliMedium(11,"ITSal",11,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMixture(9,"ITSsddCarbonM55J",aCM55J,zCM55J,dCM55J,4,wCM55J);
//     AliMedium(9,"ITSsddCarbonM55J",9,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMixture(10,"SDD AIR$",aAir,zAir,dAir,4,wAir);
//     AliMedium(10,"SDD AIR$",10,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

//     AliMixture(12, "WATER",aWater,zWater,dWater,2,wWater);
//     AliMedium(12,"WATER",12,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//      AliMixture(69,"ITSsddCAlM55J",aALCM55J,zALCM55J,dALCM55J,5,wALCM55J);
//     AliMedium(69,"ITSsddCAlM55J",69,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
  
//     AliMixture(70, "ITSsddKAPTON_POLYCH2", aKapton, zKapton, dKapton, 4, wKapton);
//     AliMedium(70,"ITSsddKAPTON_POLYCH2",70,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMixture(77,"SDDX7Rcapacitors",aX7R,zX7R,dX7R,7,wX7R);
//     AliMedium(77,"SDDX7Rcapacitors",77,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMixture(78,"SDD ruby sph. Al2O3$",aAlOxide,zAlOxide,dAlOxide,2,wAlOxide);
//     AliMedium(78,"SDD ruby sph. Al2O3$",78,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);


//     AliMaterial(64,"ALUMINUM$",0.26982E+02,0.13000E+02,0.26989E+01,0.89000E+01,0.99900E+03);
//     AliMedium(64,"ALUMINUM$",64,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMaterial(14,"COPPER",0.63546E+02,0.29000E+02,0.89600E+01,0.14300E+01,0.99900E+03);
//     AliMedium(14,"COPPER",14,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMaterial(2,"SPD SI CHIP$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
//     AliMedium(2,"SPD SI CHIP$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

//     AliMaterial(3,"SPD SI BUS$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
//     AliMedium(3,"SPD SI BUS$",3,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

//     AliMixture(4,"C (M55J)$",aCM55J,zCM55J,dCM55J,4,wCM55J);
//     AliMedium(4,"C (M55J)$",4,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

 
//     AliMixture(6,"GEN AIR$",aAir,zAir,dAir,4,wAir);
//     AliMedium(6,"GEN AIR$",6,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

//     AliMixture(13,"Freon$",afre,zfre,densfre,-2,wfre);
//     AliMedium(13,"Freon$",13,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);


//     AliMixture(15,"CERAMICS$",acer,zcer,denscer,5,wcer);
//     AliMedium(15,"CERAMICS$",15,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMixture(20,"SSD C (M55J)$",aCM55J,zCM55J,dCM55J,4,wCM55J);
//     AliMedium(20,"SSD C (M55J)$",20,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMixture(21,"SSD AIR$",aAir,zAir,dAir,4,wAir);
//     AliMedium(21,"SSD AIR$",21,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

//     AliMixture(25,"G10FR4$",aG10FR4,zG10FR4,densG10FR4,14,wG10FR4);
//     AliMedium(25,"G10FR4$",25,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//      AliMixture(26,"GEN C (M55J)$",aCM55J,zCM55J,dCM55J,4,wCM55J);
//     AliMedium(26,"GEN C (M55J)$",26,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMixture(27,"GEN Air$",aAir,zAir,dAir,4,wAir);
//     AliMedium(27,"GEN Air$",27,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

//     AliMaterial(51,"SPD SI$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
//     AliMedium(51,"SPD SI$",51,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

//     AliMaterial(52,"SPD SI CHIP$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
//     AliMedium(52,"SPD SI CHIP$",52,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

//     AliMaterial(53,"SPD SI BUS$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
//     AliMedium(53,"SPD SI BUS$",53,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

//     AliMixture(54,"SPD C (M55J)$",aCM55J,zCM55J,dCM55J,4,wCM55J);
//     AliMedium(54,"SPD C (M55J)$",54,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMixture(55,"SPD AIR$",aAir,zAir,dAir,4,wAir);
//     AliMedium(55,"SPD AIR$",55,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

//     AliMixture(56, "SPD KAPTON(POLYCH2)", aKapton, zKapton, dKapton, 4, wKapton);
//     AliMedium(56,"SPD KAPTON(POLYCH2)$",56,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMixture(61,"EPOXY$",aEpoxy,zEpoxy,dEpoxy,-3,wEpoxy);
//     AliMedium(61,"EPOXY$",61,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMaterial(62,"SILICON$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
//     AliMedium(62,"SILICON$",62,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

//     AliMixture(63, "KAPTONH(POLYCH2)", aKapton, zKapton, dKapton, 4, wKapton);
//     AliMedium(63,"KAPTONH(POLYCH2)$",63,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);


//     AliMixture(65,"INOX$",aINOX,zINOX,dINOX,9,wINOX);
//     AliMedium(65,"INOX$",65,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMixture(68,"ROHACELL$",arohac,zrohac,drohac,-4,wrohac);
//     AliMedium(68,"ROHACELL$",68,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);


//      AliMaterial(71,"ITS SANDW A$",0.12011E+02,0.60000E+01,0.2115E+00,0.17479E+03,0.99900E+03);
//     AliMedium(71,"ITS SANDW A$",71,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMaterial(72,"ITS SANDW B$",0.12011E+02,0.60000E+01,0.27000E+00,0.18956E+03,0.99900E+03);
//     AliMedium(72,"ITS SANDW B$",72,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMaterial(73,"ITS SANDW C$",0.12011E+02,0.60000E+01,0.41000E+00,0.90868E+02,0.99900E+03);
//     AliMedium(73,"ITS SANDW C$",73,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMaterial(74,"HEAT COND GLUE$",0.12011E+02,0.60000E+01,0.1930E+01,0.22100E+02,0.99900E+03);
//     AliMedium(74,"HEAT COND GLUE$",74,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMaterial(75,"ELASTO SIL$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
//     AliMedium(75,"ELASTO SIL$",75,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMaterial(76,"SPDBUS(AL+KPT+EPOX)$",0.19509E+02,0.96502E+01,0.19060E+01,0.15413E+02,0.99900E+03);
//     AliMedium(76,"SPDBUS(AL+KPT+EPOX)$",76,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
               

//     AliMixture(80,"SDD HV microcable$",aHVm,zHVm,dHVm,5,wHVm);
//     AliMedium(80,"SDD HV microcable$",80,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMixture(81,"SDD LV+signal cable$",aLVm,zLVm,dLVm,5,wLVm);
//     AliMedium(81,"SDD LV+signal cable$",81,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMixture(82,"SDD hybrid microcab$",aHLVm, zHLVm,dHLVm,5,wHLVm);
//     AliMedium(82,"SDD hybrid microcab$",82,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMixture(83,"SDD anode microcab$",aALVm,zALVm,dALVm,5,wALVm);
//     AliMedium(83,"SDD anode microcab$",83,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMaterial(84,"SDD/SSD rings$",0.123565E+02,0.64561E+01,0.18097E+01,0.229570E+02,0.99900E+03);
//     AliMedium(84,"SDD/SSD rings$",84,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

//     AliMixture(85,"inox/alum$",aInAl,zInAl,dInAl,5,wInAl);
//     AliMedium(85,"inox/alum$",85,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);


//     // special media to take into account services in the SDD and SSD 
//     // cones for the FMD

//        Float_t aA[13],zZ[13],wW[13],den;
//     // From Pierluigi Barberis calculations of 2SPD+1SDD October 2 2002.
//     zZ[0] = 1.0; aA[0] = 1.00794; // Hydrogen
//     zZ[1] = 6.0; aA[1] = 12.011; // Carbon
//     zZ[2] = 7.0; aA[2] = 14.00674; // Nitrogen
//     zZ[3] = 8.0; aA[3] = 15.9994; // Oxigen
//     zZ[4] = 14.0; aA[4] = 28.0855; // Silicon
//     zZ[5] = 24.0; aA[5] = 51.9961; //Cromium
//     zZ[6] = 25.0; aA[6] = 54.938049; // Manganese
//     zZ[7] = 26.0; aA[7] = 55.845; // Iron
//     zZ[8] = 28.0; aA[8] = 58.6934; // Nickle
//     zZ[9] = 29.0; aA[9] = 63.546; // Copper
//     zZ[10] = 13.0; aA[10] = 26.981539; // Alulminum
//     zZ[11] = 47.0; aA[11] = 107.8682; // Silver
//     zZ[12] = 27.0; aA[12] = 58.9332; // Cobolt
//     wW[0] = 0.019965;
//     wW[1] = 0.340961;
//     wW[2] = 0.041225;
//     wW[3] = 0.200352;
//     wW[4] = 0.000386;
//     wW[5] = 0.001467;
//     wW[6] = 0.000155;
//     wW[7] = 0.005113;
//     wW[8] = 0.000993;
//     wW[9] = 0.381262;
//     wW[10] = 0.008121;
//     wW[11] = 0.000000;
//     wW[12] = 0.000000;
//     if(fByThick){// New values seeITS_MatBudget_4B.xls
// 	den = 1.5253276; // g/cm^3  Cell O370
//     }else{
// 	den = 2.58423412; // g/cm^3 Cell L370
//     } // end if fByThick
//     //den = 6161.7/(3671.58978);//g/cm^3 Volume does not exclude holes
//     AliMixture(86,"AIRFMDSDD$",aA,zZ,den,+11,wW);
//     AliMedium(86,"AIRFMDSDD$",86,0,ifield,fieldm,tmaxfdAir,stemaxAir,
// 	      deemaxAir,epsilAir,stminAir);


//     wW[0] = 0.019777;
//     wW[1] = 0.325901;
//     wW[2] = 0.031848;
//     wW[3] = 0.147668;
//     wW[4] = 0.030609;
//     wW[5] = 0.013993;
//     wW[6] = 0.001479;
//     wW[7] = 0.048792;
//     wW[8] = 0.009477;
//     wW[9] = 0.350697;
//     wW[10] = 0.014546;
//     wW[11] = 0.005213;
//     wW[12] = 0.000000;
//     if(fByThick){// New values seeITS_MatBudget_4B.xls
// 	den = 1.2464275; // g/cm^3   Cell O403
//     }else{
// 	den = 1.28134409; // g/cm^3  Cell L403
//     } // end if fByThick
//     //den = 7666.3/(9753.553259); // volume does not exclude holes
//     AliMixture(87,"AIRFMDSSD$",aA,zZ,den,+12,wW); 
//     AliMedium(87,"AIRFMDSSD$",87,0,ifield,fieldm,tmaxfdAir,stemaxAir,
// 	      deemaxAir,epsilAir,stminAir);

//     wW[0] = 0.016302;
//     wW[1] = 0.461870;
//     wW[2] = 0.033662;
//     wW[3] = 0.163595;
//     wW[4] = 0.000315;
//     wW[5] = 0.001197;
//     wW[6] = 0.000127;
//     wW[7] = 0.004175;
//     wW[8] = 0.000811;
//     wW[9] = 0.311315;
//     wW[10] = 0.006631;
//     wW[11] = 0.000000;
//     wW[12] = 0.000000;
//     if(fByThick){// New values seeITS_MatBudget_4B.xls
// 	den = 1.9353276; // g/cm^3  Cell N370
//     }else{
// 	den = 3.2788626; // g/cm^3 Cell F370
//     } // end if fByThick
//     //den = 7667.1/(3671.58978); // Volume does not excludeholes
//     AliMixture(88,"ITS SANDW CFMDSDD$",aA,zZ,den,+11,wW); 
//     AliMedium(88,"ITS SANDW CFMDSDD$",88,0,ifield,fieldm,tmaxfd,stemax,
// 	      deemax,epsil,stmin);

//     wW[0] = 0.014065;
//     wW[1] = 0.520598;
//     wW[2] = 0.022650;
//     wW[3] = 0.105018;
//     wW[4] = 0.021768;
//     wW[5] = 0.009952;
//     wW[6] = 0.001051;
//     wW[7] = 0.034700;
//     wW[8] = 0.006740;
//     wW[9] = 0.249406;
//     wW[10] = 0.010345;
//     wW[11] = 0.0003707;
//     wW[12] = 0.000000;
//     if(fByThick){// New values seeITS_MatBudget_4B.xls
// 	den = 1.6564275; // g/cm^3  Cell N304
//     }else{
// 	den = 1.7028296; // g/cm^3  Cell F304
//     } // end if fByThick
//     //den = 1166.5/(3671.58978); // Volume does not exclude holes
//     AliMixture(89,"ITS SANDW CFMDSSD$",aA,zZ,den,+12,wW); 
//     AliMedium(89,"ITS SANDW CFMDSSD$",89,0,ifield,fieldm,tmaxfd,stemax,
// 	      deemax,epsil,stmin);

//     wW[0] = 0.005970;
//     wW[1] = 0.304704;
//     wW[2] = 0.042510;
//     wW[3] = 0.121715;
//     wW[4] = 0.001118;
//     wW[5] = 0.030948;
//     wW[6] = 0.003270;
//     wW[7] = 0.107910;
//     wW[8] = 0.020960;
//     wW[9] = 0.360895;
//     wW[10] = 0.000000;
//     wW[11] = 0.000000;
//     wW[12] = 0.000000;
//     if(fByThick){// New values seeITS_MatBudget_4B.xls
// 	den = 80.31136576; // g/cm^3 Cell H329
//     }else{
// 	den = 87.13062; // g/cm^3  Cell G329
//     } // end if fByThick
//     //den = 1251.3/(0.05*2.0*TMath::Pi()*(7.75*7.75 - 3.7*3.7)); // g/cm^3
//     AliMixture(97,"SPD SERVICES$",aA,zZ,den,+10,wW); 
//     AliMedium(97,"SPD SERVICES$",97,0,ifield,fieldm,tmaxfd,stemax,
// 	      deemax,epsil,stmin);

//     // Special media

//     AliMaterial(90,"SPD shield$", 12.011, 6., 1.93/10. , 22.1*10., 999);
//     AliMedium(90,"SPD shield$",90,0,ifield,fieldm,tmaxfdServ,stemaxServ,
// 	      deemaxServ,epsilServ,stminServ);

//     AliMaterial(91, "SPD End ladder$", 47.0447, 21.7963, 3.6374, 4.4711, 999); 
//     AliMedium(91,"SPD End ladder$",91,0,ifield,fieldm,tmaxfdServ,stemaxServ,
// 	      deemaxServ,epsilServ,stminServ);

//     AliMaterial(92, "SPD cone$",28.0855, 14., 2.33, 9.36, 999);    
//     AliMedium(92,"SPD cone$",92,0,ifield,fieldm,tmaxfdServ,stemaxServ,
// 	      deemaxServ,epsilServ,stminServ);

// //     Material with fractional Z not actually used
// //     AliMaterial(93, "SDD End ladder$", 69.9298, 29.8246, 0.3824, 36.5103, 999);
// //     AliMedium(93,"SDD End ladder$",93,0,ifield,fieldm,tmaxfdServ,stemaxServ,
// //               deemaxServ,epsilServ,stminServ);

//     AliMaterial(94, "SDD cone$",63.546, 29., 1.15, 1.265, 999);
//     AliMedium(94,"SDD cone$",94,0,ifield,fieldm,tmaxfdServ,stemaxServ,
// 	      deemaxServ,epsilServ,stminServ);

// //     Material with fractional Z not actually used
// //     AliMaterial(95, "SSD End ladder$", 32.0988, 15.4021, 0.68, 35.3238, 999); 
// //     AliMedium(95,"SSD End ladder$",95,0,ifield,fieldm,tmaxfdServ,stemaxServ,
// //     deemaxServ,epsilServ,stminServ);

//     AliMaterial(96, "SSD cone$",63.546, 29., 1.15, 1.265, 999);
//     AliMedium(96,"SSD cone$",96,0,ifield,fieldm,tmaxfdServ,stemaxServ,
// 	      deemaxServ,epsilServ,stminServ);


    Int_t   ifield = gAlice->Field()->Integ();
    Float_t fieldm = gAlice->Field()->Max();

    Float_t tmaxfd = 0.1; // 1.0; // Degree
    Float_t stemax = 1.0; // cm
    Float_t deemax = 0.1; // 30.0; // Fraction of particle's energy 0<deemax<=1
    Float_t epsil  = 1.0E-4; // 1.0; // cm
    Float_t stmin  = 0.0; // cm "Default value used"

    Float_t tmaxfdSi = 0.1; // .10000E+01; // Degree
    Float_t stemaxSi = 0.0075; //  .10000E+01; // cm
    Float_t deemaxSi = 0.1; // 0.30000E-02; // Fraction of particle's energy 0<deemax<=1
    Float_t epsilSi  = 1.0E-4;// .10000E+01;
    Float_t stminSi  = 0.0; // cm "Default value used"

    Float_t tmaxfdAir = 0.1; // .10000E+01; // Degree
    Float_t stemaxAir = .10000E+01; // cm
    Float_t deemaxAir = 0.1; // 0.30000E-02; // Fraction of particle's energy 0<deemax<=1
    Float_t epsilAir  = 1.0E-4;// .10000E+01;
    Float_t stminAir  = 0.0; // cm "Default value used"

    Float_t tmaxfdServ = 1.0; // 10.0; // Degree
    Float_t stemaxServ = 1.0; // 0.01; // cm
    Float_t deemaxServ = 0.5; // 0.1; // Fraction of particle's energy 0<deemax<=1
    Float_t epsilServ  = 1.0E-3; // 0.003; // cm
    Float_t stminServ  = 0.0; //0.003; // cm "Default value used"

    // Freon PerFluorobuthane C4F10 see 
    // http://st-support-cooling-electronics.web.cern.ch/
    //        st-support-cooling-electronics/default.htm
    Float_t afre[2]  = { 12.011,18.9984032 };
    Float_t zfre[2]  = { 6., 9. };
    Float_t wfre[2]  = { 4.,10. };
    Float_t densfre  = 1.52;


    //CM55J

    Float_t aCM55J[4]={12.0107,14.0067,15.9994,1.00794};
    Float_t zCM55J[4]={6.,7.,8.,1.};
    Float_t wCM55J[4]={0.908508078,0.010387573,0.055957585,0.025146765};
    Float_t dCM55J = 1.63;

    //ALCM55J

    Float_t aALCM55J[5]={12.0107,14.0067,15.9994,1.00794,26.981538};
    Float_t zALCM55J[5]={6.,7.,8.,1.,13.};
    Float_t wALCM55J[5]={0.817657902,0.0093488157,0.0503618265,0.0226320885,0.1};
    Float_t dALCM55J = 1.9866;

    //Si Chips

    Float_t aSICHIP[6]={12.0107,14.0067,15.9994,1.00794,28.0855,107.8682};
    Float_t zSICHIP[6]={6.,7.,8.,1.,14., 47.};
    Float_t wSICHIP[6]={0.039730642,0.001396798,0.01169634,0.004367771,0.844665,0.09814344903};
    Float_t dSICHIP = 2.36436;

    //Inox
    
    Float_t aINOX[9]={12.0107,54.9380, 28.0855,30.9738,32.066,58.6928,55.9961,95.94,55.845};
    Float_t zINOX[9]={6.,25.,14.,15.,16., 28.,24.,42.,26.};
    Float_t wINOX[9]={0.0003,0.02,0.01,0.00045,0.0003,0.12,0.17,0.025,0.654};
    Float_t dINOX = 8.03;

    //SDD HV microcable

    Float_t aHVm[5]={12.0107,1.00794,14.0067,15.9994,26.981538};
    Float_t zHVm[5]={6.,1.,7.,8.,13.};
    Float_t wHVm[5]={0.520088819984,0.01983871336,0.0551367996,0.157399667056, 0.247536};
    Float_t dHVm = 1.6087;

    //SDD LV+signal cable

    Float_t aLVm[5]={12.0107,1.00794,14.0067,15.9994,26.981538};
    Float_t zLVm[5]={6.,1.,7.,8.,13.};
    Float_t wLVm[5]={0.21722436468,0.0082859922,0.023028867,0.06574077612, 0.68572};
    Float_t dLVm = 2.1035;

    //SDD hybrid microcab

    Float_t aHLVm[5]={12.0107,1.00794,14.0067,15.9994,26.981538};
    Float_t zHLVm[5]={6.,1.,7.,8.,13.};
    Float_t wHLVm[5]={0.24281879711,0.00926228815,0.02574224025,0.07348667449, 0.64869};
    Float_t dHLVm = 2.0502;

    //SDD anode microcab

    Float_t aALVm[5]={12.0107,1.00794,14.0067,15.9994,26.981538};
    Float_t zALVm[5]={6.,1.,7.,8.,13.};
    Float_t wALVm[5]={0.392653705471,0.0128595919215,0.041626868025,0.118832707289, 0.431909};
    Float_t dALVm = 2.0502;

    //X7R capacitors

    Float_t aX7R[7]={137.327,47.867,15.9994,58.6928,63.5460,118.710,207.2};
    Float_t zX7R[7]={56.,22.,8.,28.,29.,50.,82.};
    Float_t wX7R[7]={0.251639432,0.084755042,0.085975822,0.038244751,0.009471271,0.321736471,0.2081768};
    Float_t dX7R = 7.14567;

    // AIR

    Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
    Float_t zAir[4]={6.,7.,8.,18.};
    Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
    Float_t dAir = 1.20479E-3;

    // Water

    Float_t aWater[2]={1.00794,15.9994};
    Float_t zWater[2]={1.,8.};
    Float_t wWater[2]={0.111894,0.888106};
    Float_t dWater   = 1.0;

    // CERAMICS
  //     94.4% Al2O3 , 2.8% SiO2 , 2.3% MnO , 0.5% Cr2O3
    Float_t acer[5]  = { 26.981539,15.9994,28.0855,54.93805,51.9961 };
    Float_t zcer[5]  = {       13.,     8.,    14.,     25.,    24. };
    Float_t wcer[5]  = {.4443408,.5213375,.0130872,.0178135,.003421};
    Float_t denscer  = 3.6;

    //G10FR4

    Float_t zG10FR4[14] = {14.00,	20.00,	13.00,	12.00,	5.00,	22.00,	11.00,	19.00,	26.00,	9.00,	8.00,	6.00,	7.00,	1.00};
    Float_t aG10FR4[14] = {28.0855000,40.0780000,26.9815380,24.3050000,10.8110000,47.8670000,22.9897700,39.0983000,55.8450000,18.9984000,15.9994000,12.0107000,14.0067000,1.0079400};
    Float_t wG10FR4[14] = {0.15144894,0.08147477,0.04128158,0.00904554,0.01397570,0.00287685,0.00445114,0.00498089,0.00209828,0.00420000,0.36043788,0.27529426,0.01415852,0.03427566};
    Float_t densG10FR4= 1.8;
    
     //--- EPOXY  --- C18 H19 O3
      Float_t aEpoxy[3] = {15.9994, 1.00794, 12.0107} ; 
      Float_t zEpoxy[3] = {     8.,      1.,      6.} ; 
      Float_t wEpoxy[3] = {     3.,     19.,     18.} ; 
      Float_t dEpoxy = 1.8 ;

      // rohacell: C9 H13 N1 O2
    Float_t arohac[4] = {12.01,  1.01, 14.010, 16.};
    Float_t zrohac[4] = { 6.,    1.,    7.,     8.};
    Float_t wrohac[4] = { 9.,   13.,    1.,     2.};
    Float_t drohac    = 0.05;

    // If he/she means stainless steel (inox) + Aluminium and Zeff=15.3383 then
//
// %Al=81.6164 %inox=100-%Al

    Float_t aInAl[5] = {27., 55.847,51.9961,58.6934,28.0855 };
    Float_t zInAl[5] = {13., 26.,24.,28.,14. };
    Float_t wInAl[5] = {.816164, .131443,.0330906,.0183836,.000919182};
    Float_t dInAl    = 3.075;

    // Kapton

    Float_t aKapton[4]={1.00794,12.0107, 14.010,15.9994};
    Float_t zKapton[4]={1.,6.,7.,8.};
    Float_t wKapton[4]={0.026362,0.69113,0.07327,0.209235};
    Float_t dKapton   = 1.42;

    //SDD ruby sph.
    Float_t aAlOxide[2]  = { 26.981539,15.9994};
    Float_t zAlOxide[2]  = {       13.,     8.};
    Float_t wAlOxide[2]  = {0.4707, 0.5293};
    Float_t dAlOxide     = 3.97;


    AliMaterial(1,"SI$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(1,"SI$",1,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMaterial(2,"SPD SI CHIP$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(2,"SPD SI CHIP$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMaterial(3,"SPD SI BUS$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(3,"SPD SI BUS$",3,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMixture(4,"C (M55J)$",aCM55J,zCM55J,dCM55J,4,wCM55J);
    AliMedium(4,"C (M55J)$",4,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(5,"AIR$",aAir,zAir,dAir,4,wAir);
    AliMedium(5,"AIR$",5,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

    AliMixture(6,"GEN AIR$",aAir,zAir,dAir,4,wAir);
    AliMedium(6,"GEN AIR$",6,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

    AliMixture(7,"SDD SI CHIP$",aSICHIP,zSICHIP,dSICHIP,6,wSICHIP);
    AliMedium(7,"SDD SI CHIP$",7,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMixture(9,"SDD C (M55J)$",aCM55J,zCM55J,dCM55J,4,wCM55J);
    AliMedium(9,"SDD C (M55J)$",9,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(10,"SDD AIR$",aAir,zAir,dAir,4,wAir);
    AliMedium(10,"SDD AIR$",10,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

    AliMaterial(11,"AL$",0.26982E+02,0.13000E+02,0.26989E+01,0.89000E+01,0.99900E+03);
    AliMedium(11,"AL$",11,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(12, "Water$",aWater,zWater,dWater,2,wWater);
    AliMedium(12,"WATER$",12,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(13,"Freon$",afre,zfre,densfre,-2,wfre);
    AliMedium(13,"Freon$",13,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(14,"COPPER$",0.63546E+02,0.29000E+02,0.89600E+01,0.14300E+01,0.99900E+03);
    AliMedium(14,"COPPER$",14,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
    AliMixture(15,"CERAMICS$",acer,zcer,denscer,5,wcer);
    AliMedium(15,"CERAMICS$",15,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(20,"SSD C (M55J)$",aCM55J,zCM55J,dCM55J,4,wCM55J);
    AliMedium(20,"SSD C (M55J)$",20,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(21,"SSD AIR$",aAir,zAir,dAir,4,wAir);
    AliMedium(21,"SSD AIR$",21,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

    AliMixture(25,"G10FR4$",aG10FR4,zG10FR4,densG10FR4,14,wG10FR4);
    AliMedium(25,"G10FR4$",25,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

     AliMixture(26,"GEN C (M55J)$",aCM55J,zCM55J,dCM55J,4,wCM55J);
    AliMedium(26,"GEN C (M55J)$",26,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(27,"GEN Air$",aAir,zAir,dAir,4,wAir);
    AliMedium(27,"GEN Air$",27,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

    AliMaterial(51,"SPD SI$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(51,"SPD SI$",51,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMaterial(52,"SPD SI CHIP$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(52,"SPD SI CHIP$",52,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMaterial(53,"SPD SI BUS$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(53,"SPD SI BUS$",53,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMixture(54,"SPD C (M55J)$",aCM55J,zCM55J,dCM55J,4,wCM55J);
    AliMedium(54,"SPD C (M55J)$",54,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(55,"SPD AIR$",aAir,zAir,dAir,4,wAir);
    AliMedium(55,"SPD AIR$",55,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

    AliMixture(56, "SPD KAPTON(POLYCH2)", aKapton, zKapton, dKapton, 4, wKapton);
    AliMedium(56,"SPD KAPTON(POLYCH2)$",56,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(61,"EPOXY$",aEpoxy,zEpoxy,dEpoxy,-3,wEpoxy);
    AliMedium(61,"EPOXY$",61,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(62,"SILICON$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(62,"SILICON$",62,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

    AliMixture(63, "KAPTONH(POLYCH2)", aKapton, zKapton, dKapton, 4, wKapton);
    AliMedium(63,"KAPTONH(POLYCH2)$",63,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(64,"ALUMINUM$",0.26982E+02,0.13000E+02,0.26989E+01,0.89000E+01,0.99900E+03);
    AliMedium(64,"ALUMINUM$",64,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(65,"INOX$",aINOX,zINOX,dINOX,9,wINOX);
    AliMedium(65,"INOX$",65,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(68,"ROHACELL$",arohac,zrohac,drohac,-4,wrohac);
    AliMedium(68,"ROHACELL$",68,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

     AliMixture(69,"SDD C AL (M55J)$",aALCM55J,zALCM55J,dALCM55J,5,wALCM55J);
    AliMedium(69,"SDD C AL (M55J)$",69,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
  
    AliMixture(70, "SDDKAPTON (POLYCH2)", aKapton, zKapton, dKapton, 4, wKapton);
    AliMedium(70,"SDDKAPTON (POLYCH2)$",70,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

     AliMaterial(71,"ITS SANDW A$",0.12011E+02,0.60000E+01,0.2115E+00,0.17479E+03,0.99900E+03);
    AliMedium(71,"ITS SANDW A$",71,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(72,"ITS SANDW B$",0.12011E+02,0.60000E+01,0.27000E+00,0.18956E+03,0.99900E+03);
    AliMedium(72,"ITS SANDW B$",72,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(73,"ITS SANDW C$",0.12011E+02,0.60000E+01,0.41000E+00,0.90868E+02,0.99900E+03);
    AliMedium(73,"ITS SANDW C$",73,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(74,"HEAT COND GLUE$",0.12011E+02,0.60000E+01,0.1930E+01,0.22100E+02,0.99900E+03);
    AliMedium(74,"HEAT COND GLUE$",74,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(75,"ELASTO SIL$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(75,"ELASTO SIL$",75,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    // SPD bus (data from Petra Riedler)
    Float_t aSPDbus[5] = {1.00794,12.0107,14.01,15.9994,26.982 };
    Float_t zSPDbus[5] = {1.,6.,7.,8.,13.};
    Float_t wSPDbus[5] = {0.023523,0.318053,0.009776,0.078057,0.570591};
    Float_t dSPDbus    = 2.128505;

    //   AliMaterial(76,"SPDBUS(AL+KPT+EPOX)$",0.19509E+02,0.96502E+01,0.19060E+01,0.15413E+02,0.99900E+03);
    AliMixture(76,"SPDBUS(AL+KPT+EPOX)$",aSPDbus,zSPDbus,dSPDbus,5,wSPDbus);
    AliMedium(76,"SPDBUS(AL+KPT+EPOX)$",76,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
               
    AliMixture(77,"SDD X7R capacitors$",aX7R,zX7R,dX7R,7,wX7R);
    AliMedium(77,"SDD X7R capacitors$",77,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(78,"SDD ruby sph. Al2O3$",aAlOxide,zAlOxide,dAlOxide,2,wAlOxide);
    AliMedium(78,"SDD ruby sph. Al2O3$",78,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMaterial(79,"SDD SI insensitive$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(79,"SDD SI insensitive$",79,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(80,"SDD HV microcable$",aHVm,zHVm,dHVm,5,wHVm);
    AliMedium(80,"SDD HV microcable$",80,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(81,"SDD LV+signal cable$",aLVm,zLVm,dLVm,5,wLVm);
    AliMedium(81,"SDD LV+signal cable$",81,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(82,"SDD hybrid microcab$",aHLVm, zHLVm,dHLVm,5,wHLVm);
    AliMedium(82,"SDD hybrid microcab$",82,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(83,"SDD anode microcab$",aALVm,zALVm,dALVm,5,wALVm);
    AliMedium(83,"SDD anode microcab$",83,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
    Float_t aDSring[4]={12.0107,      1.00794,     14.0067,      15.9994};
    Float_t zDSring[4]={ 6.,          1.,           7.,           8.};
    Float_t wDSring[4]={ 0.854323888, 0.026408778,  0.023050265,  0.096217069};
    Float_t dDSring = 0.2875;
    AliMixture(84,"SDD/SSD rings$",aDSring,zDSring,dDSring,4,wDSring);
    AliMedium(84,"SDD/SSD rings$",84,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    AliMixture(85,"inox/alum$",aInAl,zInAl,dInAl,5,wInAl);
    AliMedium(85,"inox/alum$",85,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

    // special media to take into account services in the SDD and SSD 
    // cones for the FMD
    //Begin_Html
    /*
      <A HREF="http://www.Physics.ohio-state.edu/~nilsen/ITS/ITS_MatBudget_4B.xls">
      </pre>
      <br clear=left>
      <font size=+2 color=blue>
      <p> The Exel spread sheet from which these density number come from.
      </font></A>
    */
    //End_Html

    //  AliMaterial(86,"AIRFMDSDD$",0.14610E+02,0.73000E+01,0.12050E-02,0.30423E+05,0.99900E+03);
    Float_t aA[13],zZ[13],wW[13],den;
    // From Pierluigi Barberis calculations of 2SPD+1SDD October 2 2002.
    zZ[0] = 1.0; aA[0] = 1.00794; // Hydrogen
    zZ[1] = 6.0; aA[1] = 12.011; // Carbon
    zZ[2] = 7.0; aA[2] = 14.00674; // Nitrogen
    zZ[3] = 8.0; aA[3] = 15.9994; // Oxigen
    zZ[4] = 14.0; aA[4] = 28.0855; // Silicon
    zZ[5] = 24.0; aA[5] = 51.9961; //Cromium
    zZ[6] = 25.0; aA[6] = 54.938049; // Manganese
    zZ[7] = 26.0; aA[7] = 55.845; // Iron
    zZ[8] = 28.0; aA[8] = 58.6934; // Nickle
    zZ[9] = 29.0; aA[9] = 63.546; // Copper
    zZ[10] = 13.0; aA[10] = 26.981539; // Alulminum
    zZ[11] = 47.0; aA[11] = 107.8682; // Silver
    zZ[12] = 27.0; aA[12] = 58.9332; // Cobolt
    wW[0] = 0.019965;
    wW[1] = 0.340961;
    wW[2] = 0.041225;
    wW[3] = 0.200352;
    wW[4] = 0.000386;
    wW[5] = 0.001467;
    wW[6] = 0.000155;
    wW[7] = 0.005113;
    wW[8] = 0.000993;
    wW[9] = 0.381262;
    wW[10] = 0.008121;
    wW[11] = 0.000000;
    wW[12] = 0.000000;
    if(fByThick){// New values seeITS_MatBudget_4B.xls
	den = 1.5253276; // g/cm^3  Cell O370
    }else{
	den = 2.58423412; // g/cm^3 Cell L370
    } // end if fByThick
    //den = 6161.7/(3671.58978);//g/cm^3 Volume does not exclude holes
    AliMixture(86,"AIRFMDSDD$",aA,zZ,den,+11,wW);
    AliMedium(86,"AIRFMDSDD$",86,0,ifield,fieldm,tmaxfdAir,stemaxAir,
	      deemaxAir,epsilAir,stminAir);

    //AliMaterial(87,"AIRFMDSSD$",0.14610E+02,0.73000E+01,0.12050E-02,0.30423E+05,0.99900E+03);
    // From Pierluigi Barberis calculations of SSD October 2 2002.
    wW[0] = 0.019777;
    wW[1] = 0.325901;
    wW[2] = 0.031848;
    wW[3] = 0.147668;
    wW[4] = 0.030609;
    wW[5] = 0.013993;
    wW[6] = 0.001479;
    wW[7] = 0.048792;
    wW[8] = 0.009477;
    wW[9] = 0.350697;
    wW[10] = 0.014546;
    wW[11] = 0.005213;
    wW[12] = 0.000000;
    if(fByThick){// New values seeITS_MatBudget_4B.xls
	den = 1.2464275; // g/cm^3   Cell O403
    }else{
	den = 1.28134409; // g/cm^3  Cell L403
    } // end if fByThick
    //den = 7666.3/(9753.553259); // volume does not exclude holes
    AliMixture(87,"AIRFMDSSD$",aA,zZ,den,+12,wW); 
    AliMedium(87,"AIRFMDSSD$",87,0,ifield,fieldm,tmaxfdAir,stemaxAir,
	      deemaxAir,epsilAir,stminAir);

    //AliMaterial(88,"ITS SANDW CFMDSDD$",0.12011E+02,0.60000E+01,0.41000E+00,0.90868E+02,0.99900E+03);
    // From Pierluigi Barberis calculations of 1SDD+Carbon fiber October 2 2002
    wW[0] = 0.016302;
    wW[1] = 0.461870;
    wW[2] = 0.033662;
    wW[3] = 0.163595;
    wW[4] = 0.000315;
    wW[5] = 0.001197;
    wW[6] = 0.000127;
    wW[7] = 0.004175;
    wW[8] = 0.000811;
    wW[9] = 0.311315;
    wW[10] = 0.006631;
    wW[11] = 0.000000;
    wW[12] = 0.000000;
    if(fByThick){// New values seeITS_MatBudget_4B.xls
	den = 1.9353276; // g/cm^3  Cell N370
    }else{
	den = 3.2788626; // g/cm^3 Cell F370
    } // end if fByThick
    //den = 7667.1/(3671.58978); // Volume does not excludeholes
    AliMixture(88,"ITS SANDW CFMDSDD$",aA,zZ,den,+11,wW); 
    AliMedium(88,"ITS SANDW CFMDSDD$",88,0,ifield,fieldm,tmaxfd,stemax,
	      deemax,epsil,stmin);

    //AliMaterial(89,"ITS SANDW CFMDSSD$",0.12011E+02,0.60000E+01,0.41000E+00,0.90868E+02,0.99900E+03);
    // From Pierluigi Barberis calculations of SSD+Carbon fiber October 2 2002.
    wW[0] = 0.014065;
    wW[1] = 0.520598;
    wW[2] = 0.022650;
    wW[3] = 0.105018;
    wW[4] = 0.021768;
    wW[5] = 0.009952;
    wW[6] = 0.001051;
    wW[7] = 0.034700;
    wW[8] = 0.006740;
    wW[9] = 0.249406;
    wW[10] = 0.010345;
    wW[11] = 0.0003707;
    wW[12] = 0.000000;
    if(fByThick){// New values seeITS_MatBudget_4B.xls
	den = 1.6564275; // g/cm^3  Cell N304
    }else{
	den = 1.7028296; // g/cm^3  Cell F304
    } // end if fByThick
    //den = 1166.5/(3671.58978); // Volume does not exclude holes
    AliMixture(89,"ITS SANDW CFMDSSD$",aA,zZ,den,+12,wW); 
    AliMedium(89,"ITS SANDW CFMDSSD$",89,0,ifield,fieldm,tmaxfd,stemax,
	      deemax,epsil,stmin);

    //AliMaterial(97,"SPD SERVICES$",0.12011E+02,0.60000E+01,0.41000E+00,0.90868E+02,0.99900E+03);
    // From Pierluigi Barberis calculations of 1SPD October 2 2002.
    wW[0] = 0.005970;
    wW[1] = 0.304704;
    wW[2] = 0.042510;
    wW[3] = 0.121715;
    wW[4] = 0.001118;
    wW[5] = 0.030948;
    wW[6] = 0.003270;
    wW[7] = 0.107910;
    wW[8] = 0.020960;
    wW[9] = 0.360895;
    wW[10] = 0.000000;
    wW[11] = 0.000000;
    wW[12] = 0.000000;
    if(fByThick){// New values seeITS_MatBudget_4B.xls
	den = 80.31136576; // g/cm^3 Cell H329
    }else{
	den = 87.13062; // g/cm^3  Cell G329
    } // end if fByThick
    //den = 1251.3/(0.05*2.0*TMath::Pi()*(7.75*7.75 - 3.7*3.7)); // g/cm^3
    AliMixture(97,"SPD SERVICES$",aA,zZ,den,+10,wW); 
    AliMedium(97,"SPD SERVICES$",97,0,ifield,fieldm,tmaxfd,stemax,
	      deemax,epsil,stmin);


    // Special media

    AliMaterial(90,"SPD shield$", 12.011, 6., 1.93/10. , 22.1*10., 999);
    AliMedium(90,"SPD shield$",90,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);

    // SPD End Ladder (data from Petra Riedler)
    Float_t aSPDel[5] = {1.00794,12.0107,14.01,15.9994,63.54 };
    Float_t zSPDel[5] = {1.,6.,7.,8.,29.};
    Float_t wSPDel[5] = {0.004092,0.107274,0.011438,0.032476,0.844719};
    Float_t dSPDel    = 3.903403;

    //   AliMaterial(91, "SPD End ladder$", 47.0447, 21.7963, 3.6374, 4.4711, 999); 
    AliMixture(91,"SPD End ladder$",aSPDel,zSPDel,dSPDel,5,wSPDel);
    AliMedium(91,"SPD End ladder$",91,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);

    AliMaterial(92, "SPD cone$",28.0855, 14., 2.33, 9.36, 999);    
    AliMedium(92,"SPD cone$",92,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);
    /*  Material with fractional Z not actually used
    AliMaterial(93, "SDD End ladder$", 69.9298, 29.8246, 0.3824, 36.5103, 999);
    AliMedium(93,"SDD End ladder$",93,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);
    */
    AliMaterial(94, "SDD cone$",63.546, 29., 1.15, 1.265, 999);
    AliMedium(94,"SDD cone$",94,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);
    /* Material with fractional Z not actually used
    AliMaterial(95, "SSD End ladder$", 32.0988, 15.4021, 0.68, 35.3238, 999); 
    AliMedium(95,"SSD End ladder$",95,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);
    */
    AliMaterial(96, "SSD cone$",63.546, 29., 1.15, 1.265, 999);
    AliMedium(96,"SSD cone$",96,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);

}


//______________________________________________________________________
void AliITSv11::InitAliITSgeom(){
  //
  // Fill fITSgeom with the 3 sub-detector geometries
  //

    const Int_t knlayers = 6;
    const Int_t kndeep = 3;
    const AliITSDetector idet[knlayers]={kSPD,kSPD,kSDD,kSDD,kSSD,kSSD};
    const TString names[knlayers] = {

      "",  // lay=1
      "",  // lay=2
      "/ALIC_1/ITSV_1/ITSsddLayer3_1/ITSsddLadd_%d/ITSsddSensor_%d/ITSsddWafer_%d", // lay=3
      "/ALIC_1/ITSV_1/ITSsddLayer4_1/ITSsddLadd_%d/ITSsddSensor_%d/ITSsddWafer_%d", // lay=4
      "",  // lay=5
      ""}; // lay=6

    const Int_t itsGeomTreeCopys[knlayers][kndeep]= {{10, 2, 4},// lay=1
                                                     {10, 4, 4},// lay=2
                                                     {14, 6, 1},// lay=3
                                                     {22, 8, 1},// lay=4
                                                     {34,22, 1},// lay=5
                                                     {38,25, 1}};//lay=6
    Int_t       nlad[knlayers],ndet[knlayers];
    Int_t       mod,lay,lad=0,det=0,i,j,k,cp0,cp1,cp2;
    TString path,shapeName;
    TGeoHMatrix materix;
    Double_t trans[3]={3*0.0},rot[10]={9*0.0,1.0};
    TArrayD shapePar;
    TArrayF shapeParF;

    AliDebug(1,"Reading Geometry transformation directly from Modler.");
    mod = 0;
    for(i=0;i<knlayers;i++){
        k = 1;
        for(j=0;j<kndeep;j++) if(itsGeomTreeCopys[i][j]!=0)
            k *= TMath::Abs(itsGeomTreeCopys[i][j]);
        mod += k;
    } // end for i

    if(GetITSgeom()!=0) delete GetITSgeom();
    SetITSgeom(0);
    nlad[0]=20;nlad[1]=40;nlad[2]=14;nlad[3]=22;nlad[4]=34;nlad[5]=38;
    ndet[0]= 4;ndet[1]= 4;ndet[2]= 6;ndet[3]= 8;ndet[4]=22;ndet[5]=25;
    AliITSgeom* geom = new AliITSgeom(0,6,nlad,ndet,mod);
    SetITSgeom(geom);
    mod = 0;
    for(lay=1;lay<=knlayers;lay++){
        for(cp0=1;cp0<=itsGeomTreeCopys[lay-1][0];cp0++){
            for(cp1=1;cp1<=itsGeomTreeCopys[lay-1][1];cp1++){
                for(cp2=1;cp2<=itsGeomTreeCopys[lay-1][2];cp2++){
                    path.Form(names[lay-1].Data(),
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
                    geom->CreateMatrix(mod,lay,lad,det,idet[lay-1],trans,rot);
                    geom->SetTrans(mod,materix.GetTranslation());
                    geom->SetRotMatrix(mod,materix.GetRotationMatrix());
                    switch (lay){
                    case 1: case 2:{
                        geom->ReSetShape(kSPD,new AliITSgeomSPD425Short(
                                shapeParF.GetSize(),shapeParF.GetArray()));
                    }break;
                    case 3: case 4:{
                        geom->ReSetShape(kSDD,new AliITSgeomSDD256(
                                shapeParF.GetSize(),shapeParF.GetArray()));
                    }break;
                    case 5: case 6:{
                        geom->ReSetShape(kSSD,new AliITSgeomSSD75and275(
                                shapeParF.GetSize(),shapeParF.GetArray()));
                    }break;
                    default:{
                    }break;
                    } // end switch
                    mod++;
                } /// end for cp2
            } // end for cp1
        } // end for cp0
    } // end for lay
    return;

//   fSDDgeom->ExportSensorGeometry(GetITSgeom(), +3, 0);  //SDD

}


//______________________________________________________________________
void AliITSv11::Init(){
  //
  //     Initialise the ITS after it has been created.
  //

  //AliInfo(Form("Minor version %d",fMinorVersion));
    //
    if(fRead[0]=='\0') strncpy(fRead,fEuclidGeomDet,60);
    if(fWrite[0]=='\0') strncpy(fWrite,fEuclidGeomDet,60);
    if(GetITSgeom()!=0) SetITSgeom(0x0);
    AliITSgeom* geom = new AliITSgeom();
    SetITSgeom(geom);
    if(fGeomDetIn) GetITSgeom()->ReadNewFile(fRead);
    else this->InitAliITSgeom();
    if(fGeomDetOut) GetITSgeom()->WriteNewFile(fWrite);
    AliITS::Init();
    //
}

//______________________________________________________________________
void AliITSv11::SetDefaults(){
  //
  // Set response ans segmentation models for SPD, SDD and SSD
  //
     const Float_t kconv = 1.0e+04; // convert cm to microns
    AliInfo("Called");    

    if(!fDetTypeSim) fDetTypeSim = new AliITSDetTypeSim();
    fDetTypeSim->SetITSgeom(GetITSgeom());
  
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
	copy = fTrackReferences->GetEntriesFast();
	TClonesArray &lTR = *fTrackReferences;
	// Fill TrackReference structure with this new TrackReference.
	new(lTR[copy]) AliTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber());
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
    //
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
