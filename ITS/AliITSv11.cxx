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

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  Inner Traking System version 11                                         //
//  This class contains the base procedures for the Inner Tracking System   //
//                                                                          //
// Authors: R. Barbera                                                      //
// version 6.                                                               //
// Created  2000.                                                           //
//                                                                          //
//  NOTE: THIS IS THE  SYMMETRIC PPR geometry of the ITS.                   //
// THIS WILL NOT WORK                                                       //
// with the geometry or module classes or any analysis classes. You are     //
// strongly encouraged to uses AliITSv5.                                    //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
// See AliITSv11::StepManager().
// General C/C++ includes
#include <stdio.h>
#include <stdlib.h>
// General Root includes
#include <Riostream.h>
#include <TMath.h>
#include <float.h>
#include <TFile.h>    // only required for Tracking function?
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TObjString.h>
// Root Geometry includes
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoPcon.h>
#include <TGeoCone.h>
#include <TGeoTube.h> // contaings TGeoTubeSeg
#include <TGeoArb8.h>
#include <TGeoCompositeShape.h>
#include <TGeoMatrix.h>
#include <TGeoNode.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include "AliITSv11GeometrySupport.h"
// General AliRoot includes
#include "AliRun.h"
#include "AliMagF.h"
#include "AliConst.h"
// ITS specific includes
#include "AliITShit.h"
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
#include "AliITSBaseGeometry.h"
#include "AliITSv11.h"
#include "AliITSv11GeometrySupport.h"

ClassImp(AliITSv11)

/*
  Some temparary #define's used untill ROOT has addoppted the proper
  Getter in it's classes.
  These Below are for TGeoPcon functions.
*/

//______________________________________________________________________
AliITSv11::AliITSv11() : AliITS(),
fEuclidOut(kFALSE),
fGeomDetOut(kFALSE),
fGeomDetIn(kFALSE),
fMajorVersion(11),
fMinorVersion(0),
fDet1(0.0),
fDet2(0.0),
fChip1(0.0),
fChip2(0.0),
fRails(0),
fFluid(1){
    // Standard default constructor for the ITS version 11.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return
    //   A default constructed AliITSv11 class.

    //fITSV = 0;
    //fcS = 0;
//   fcD = 0;
}
//______________________________________________________________________
AliITSv11::AliITSv11(const char *title) : AliITS("ITS", title),
fEuclidOut(kFALSE),
fGeomDetOut(kFALSE),
fGeomDetIn(kFALSE),
fMajorVersion(11),
fMinorVersion(0),
fDet1(0.0),
fDet2(0.0),
fChip1(0.0),
fChip2(0.0),
fRails(0),
fFluid(1){
    // Standard constructor for the ITS version 11.
    // Inputs:
    //   const char *title  The title of for this geometry.
    // Outputs:
    //   none.
    // Return
    //   A Standard constructed AliITSv11 class.
}
//______________________________________________________________________
AliITSv11::~AliITSv11() {
    // Standard destructor for the ITS version 11.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return
    //   none.

//    if(fITSV!=0) delete fITSV;
//    if(fcS!=0) delete fcS;
//    if(fcD!=0) delete fcD;
}
//______________________________________________________________________
AliITSv11::AliITSv11(const AliITSv11 &source) : AliITS(source){
    //     Copy Constructor for ITS version 11.
    // Inputs:
    //   AliITSv11 &source  class to be copied from.
    // Outputs:
    //   none.
    // Return
    //   none.

    if(&source == this) return;
    Error("Copy Constructor","Not allowed to copy AliITSv11");
    return;
}
//______________________________________________________________________
AliITSv11& AliITSv11::operator=(const AliITSv11 &source){
    //    Assignment operator for the ITS version 11.
    // Inputs:
    //   AliITSv11 &source  class to be copied from.
    // Outputs:
    //   none.
    // Return
    //   none.

    if(&source == this) return *this;
    Error("= operator","Not allowed to copy AliITSv11");
    return *this;
}
//______________________________________________________________________
void AliITSv11::BuildGeometry(){
    // This routine defines and Creates the geometry for version 11 of 
    // the ITS for use in the simulation display routines. This is a 
    // very simplified geometry for speed of viewing.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return
    //   none.
    //TVector3 t(0.0,0.0,0.0);

    //if(fITSV==0) fITSV = new AliITSGeometryITSV(this,"ALIC");
    //if(fcS==0) fcS = new AliITSGeometrySSDCone(this,t,"TSV",1);

    //fcS->BuildDisplayGeometry();
}
//______________________________________________________________________
void AliITSv11::CreateGeometry(){
    // This routine defines and Creates the geometry for version 11 of 
    // the ITS. The geometry is used by the particle trasport routines,
    // and therefore, is very detailed.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return
    //   none.
    TVector3 t(0.0,0.0,0.0);
    const Double_t kcm = 1.0;

    TGeoManager *mgr = gGeoManager;
    TGeoVolume *vALIC = mgr->GetTopVolume();

    TGeoPcon *sITS = new TGeoPcon("ITS Top Volume, Daughter of ALIC",
                                  0.0,360.0,2);
    // DefineSection(section number, Z, Rmin, Rmax).
    sITS->DefineSection(0,-300.0*kcm,0.01*kcm,50.0*kcm);
    sITS->DefineSection(1,+300.0*kcm,0.01*kcm,50.0*kcm);
    TGeoVolume *vITS = new TGeoVolume("ITSV",sITS,0);
    mgr->AddVolume(vITS);
    vITS->SetVisibility(kFALSE);
    vALIC->AddNode(vITS,1,0);
    //
    AliITSv11GeometrySupport *sup = new AliITSv11GeometrySupport(GetDebug());
    sup->SPDCone(vITS);
    sup->SPDThermalSheald(vITS);
    sup->SDDCone(vITS);
    sup->SSDCone(vITS);
    sup->ServicesCableSupport(vITS);
}
//______________________________________________________________________
void AliITSv11::CreateMaterials(){
    // Create ITS materials
    //     This function defines the default materials used in the Geant
    // Monte Carlo simulations for the geometries AliITSv11.
    // In general it is automatically replaced by
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return
    //   none.

    //TGeoMaterial *C  = new TGeoMaterial("ITSCarbon",12.0,6.0,2.265);
    TGeoMaterial *matAl = new TGeoMaterial("ITSAluminum",26.981539,13.0,2.07);
    TGeoMixture *matCfiber = new TGeoMixture("ITSCarbonFiber",6,1.930);
    TGeoMixture *matRohacell = new TGeoMixture("ITSRohacell",6,1.930);
    TGeoMixture *matStaselite = new TGeoMixture("ITSStaselite4411w",6,1.930);
    TGeoMixture *matAir = new TGeoMixture("ITSAir",6,1.205*1.E-3);
    TGeoMixture *matStainless = new TGeoMixture("ITSStainless",6,1.930);
    //
    Double_t medSPDcone[20];
    medSPDcone[0] = 1.0; // imat
    medSPDcone[1] = 0.0; // isvol
    medSPDcone[2] = gAlice->Field()->Integ(); // ifield
    medSPDcone[3] = gAlice->Field()->Max(); // fieldm
    medSPDcone[4] = 1.0; // tmaxfd [degrees]
    medSPDcone[5] = 1.0; // stemax [cm]
    medSPDcone[6] = 0.5; // deemax [fraction]
    medSPDcone[7] = 1.0E-3; // epsil [cm]
    medSPDcone[8] = 0.0; // stmin [cm]
    new TGeoMedium("ITSspdCarbonFiber",1,matCfiber,medSPDcone);
    medSPDcone[0] += 1.0;
    new TGeoMedium("ITSspdStaselite4411w",2,matStaselite,medSPDcone);
    medSPDcone[0] += 1.0;
    new TGeoMedium("ITSspdRohacell50A",3,matRohacell,medSPDcone);
    medSPDcone[0] += 1.0;
    new TGeoMedium("ITSspdStainlesSteal",4,matStainless,medSPDcone);
    medSPDcone[0] += 1.0;
    new TGeoMedium("ITSspdAir",5,matAir,medSPDcone);
    medSPDcone[0] += 1.0;
    new TGeoMedium("ITSspdAl",6,matAl,medSPDcone);
    //
    Double_t medSSDcone[20];
    medSSDcone[0] = 1.0; // imat
    medSSDcone[1] = 0.0; // isvol
    medSSDcone[2] = gAlice->Field()->Integ(); // ifield
    medSSDcone[3] = gAlice->Field()->Max(); // fieldm
    medSSDcone[4] = 1.0; // tmaxfd [degrees]
    medSSDcone[5] = 1.0; // stemax [cm]
    medSSDcone[6] = 0.5; // deemax [fraction]
    medSSDcone[7] = 1.0E-3; // epsil [cm]
    medSSDcone[8] = 0.0; // stmin [cm]
    new TGeoMedium("ITSssdCarbonFiber",1,matCfiber,medSSDcone);
    medSSDcone[0] += 1.0;
    new TGeoMedium("ITSssdStaselite4411w",2,matStaselite,medSSDcone);
    medSSDcone[0] += 1.0;
    new TGeoMedium("ITSssdRohacell50A",3,matRohacell,medSSDcone);
    medSSDcone[0] += 1.0;
    new TGeoMedium("ITSssdStainlesSteal",4,matStainless,medSSDcone);
    medSSDcone[0] += 1.0;
    new TGeoMedium("ITSssdAir",5,matAir,medSSDcone);
    medSSDcone[0] += 1.0;
    new TGeoMedium("ITSssdAl",6,matAl,medSSDcone);
}
//______________________________________________________________________
void AliITSv11::InitAliITSgeom(){
    // Based on the geometry tree defined in Geant 3.21, this
    // routine initilizes the Class AliITSgeom from the Geant 3.21 ITS 
    // geometry sturture.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return
    //   none.
}
//______________________________________________________________________
void AliITSv11::Init(){
    // Initialise the ITS after it has been created.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return
    //   none.
}
//______________________________________________________________________
void AliITSv11::SetDefaults(){
    // Sets the default segmentation, response, digit and raw cluster 
    // classes to be used. These defaults can be overwritten in the
    // macros that do these later steps. Defaults are give hear for the
    // general user.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return
    //   none.
}
//______________________________________________________________________
void AliITSv11::DrawModule(){
    // Draw a standard set of shaded view of the ITS version 11.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return
    //   none.
}
//______________________________________________________________________
void AliITSv11::StepManager(){
    // Called for every step in the ITS, then calles the AliITShit class
    // creator with the information to be recoreded about that hit.
    //  The value of the macro ALIITSPRINTGEOM if set to 1 will allow the
    // printing of information to a file which can be used to create a .det
    // file read in by the routine CreateGeometry(). If set to 0 or any other
    // value except 1, the default behavior, then no such file is created nor
    // is the extra variables and the like used in the printing allocated.
}
 
