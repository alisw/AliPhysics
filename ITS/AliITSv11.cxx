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
AliITSv11::AliITSv11() : AliITS() {
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
AliITSv11::AliITSv11(const char *title) : AliITS("ITS", title){
    // Standard constructor for the ITS version 11.
    // Inputs:
    //   const char *title  The title of for this geometry.
    // Outputs:
    //   none.
    // Return
    //   A Standard constructed AliITSv11 class.

    //fITSV = 0;
    //fcS = 0;
//    fcD = 0;
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
    TVector3 t(0.0,0.0,0.0);

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
    TGeoVolume *ALIC = mgr->GetTopVolume();

    TGeoPcon *itsv = new TGeoPcon("ITS Top Volume, Daughter of ALIC",
                                  0.0,360.0,2);
    // DefineSection(section number, Z, Rmin, Rmax).
    itsv->DefineSection(0,-300.0*kcm,0.01*kcm,50.0*kcm);
    itsv->DefineSection(1,+300.0*kcm,0.01*kcm,50.0*kcm);
    TGeoVolume *ITSV = new TGeoVolume("ITSV",itsv,0);
    //mgr->AddVolume(ITSV);
    ITSV->SetVisibility(kFALSE);
    ALIC->AddNode(ITSV,1,0);
    //
    AliITSv11GeometrySupport *sup = new AliITSv11GeometrySupport(GetDebug());
    //sup->SPDCone(ITSV);
    //sup->SDDCone(ITSV);
    sup->SSDCone(ITSV);
    //sup->ServicesCableSupport(ITSV);
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
    TGeoMaterial *Al = new TGeoMaterial("ITSAluminum",26.981539,13.0,2.07);
    TGeoMixture *Cfiber = new TGeoMixture("ITSCarbonFiber",6,1.930);
    TGeoMixture *Rohacell = new TGeoMixture("ITSRohacell",6,1.930);
    TGeoMixture *Staselite = new TGeoMixture("ITSStaselite4411w",6,1.930);
    TGeoMixture *Air = new TGeoMixture("ITSAir",6,1.205*1.E-3);
    TGeoMixture *Stainless = new TGeoMixture("ITSStainless",6,1.930);
    //
    Double_t SPDcone[20];
    SPDcone[0] = 1.0; // imat
    SPDcone[1] = 0.0; // isvol
    SPDcone[2] = gAlice->Field()->Integ(); // ifield
    SPDcone[3] = gAlice->Field()->Max(); // fieldm
    SPDcone[4] = 1.0; // tmaxfd [degrees]
    SPDcone[5] = 1.0; // stemax [cm]
    SPDcone[6] = 0.5; // deemax [fraction]
    SPDcone[7] = 1.0E-3; // epsil [cm]
    SPDcone[8] = 0.0; // stmin [cm]
    new TGeoMedium("ITSspdCarbonFiber",1,Cfiber,SPDcone);
    SPDcone[0] += 1.0;
    new TGeoMedium("ITSspdStaselite4411w",2,Staselite,SPDcone);
    SPDcone[0] += 1.0;
    new TGeoMedium("ITSspdRohacell50A",3,Rohacell,SPDcone);
    SPDcone[0] += 1.0;
    new TGeoMedium("ITSspdStainlesSteal",4,Stainless,SPDcone);
    SPDcone[0] += 1.0;
    new TGeoMedium("ITSspdAir",5,Air,SPDcone);
    SPDcone[0] += 1.0;
    new TGeoMedium("ITSspdAl",6,Al,SPDcone);
    //
    Double_t SSDcone[20];
    SSDcone[0] = 1.0; // imat
    SSDcone[1] = 0.0; // isvol
    SSDcone[2] = gAlice->Field()->Integ(); // ifield
    SSDcone[3] = gAlice->Field()->Max(); // fieldm
    SSDcone[4] = 1.0; // tmaxfd [degrees]
    SSDcone[5] = 1.0; // stemax [cm]
    SSDcone[6] = 0.5; // deemax [fraction]
    SSDcone[7] = 1.0E-3; // epsil [cm]
    SSDcone[8] = 0.0; // stmin [cm]
    new TGeoMedium("ITSssdCarbonFiber",1,Cfiber,SSDcone);
    SSDcone[0] += 1.0;
    new TGeoMedium("ITSssdStaselite4411w",2,Staselite,SSDcone);
    SSDcone[0] += 1.0;
    new TGeoMedium("ITSssdRohacell50A",3,Rohacell,SSDcone);
    SSDcone[0] += 1.0;
    new TGeoMedium("ITSssdStainlesSteal",4,Stainless,SSDcone);
    SSDcone[0] += 1.0;
    new TGeoMedium("ITSssdAir",5,Air,SSDcone);
    SSDcone[0] += 1.0;
    new TGeoMedium("ITSssdAl",6,Al,SSDcone);
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
 
