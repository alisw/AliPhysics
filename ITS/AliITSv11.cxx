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
#include <TFile.h>    // only required for Tracking function?
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TObjString.h>
// Root Geometry includes
#include <TGeoManager.h>
#include <TGeoPcon.h>
#include <TGeoTube.h>
#include <TGeoNode.h>
#include <TGeoMaterial.h>
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

// Units, Convert from k?? to cm,degree,GeV,seconds,
const Double_t kmm = 0.10; // Convert mm to TGeom's cm.
const Double_t kcm = 1.00; // Convert cv to TGeom's cm.

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

    TGeoManager *mgr = gGeoManager;
    TGeoVolume *ALIC = mgr->GetTopVolume();

    TGeoPcon *itsv = new TGeoPcon("ITS Top Volume, Daughter of ALIC",0.0,360.0,2);
    // DefineSection(section number, Z, Rmin, Rmax).
    itsv->DefineSection(0,-100.0*kcm,0.01*kcm,50.0*kcm);
    itsv->DefineSection(1,+100.0*kcm,0.01*kcm,50.0*kcm);
    TGeoVolume *ITSV = new TGeoVolume("ITSV",itsv,0);
    mgr->AddVolume(ITSV);
    ALIC->AddNode(ITSV,1,0);
    //
    SSDCone(ITSV);
}
//______________________________________________________________________
Double_t AliITSv11::RmaxFrom2Points(TGeoPcon *p,Int_t i1,Int_t i2,Double_t z){
    // functions Require at parts of Volume A to be already defined.
    // Retruns the value of Rmax corresponding to point z alone the line
    // defined by the two points p.Rmax(i1),p-GetZ(i1) and p->GetRmax(i2),
    // p->GetZ(i2).

    return p->GetRmax(i2)+(p->GetRmax(i1)-p->GetRmax(i2))*(z-p->GetZ(i2))/
	(p->GetZ(i1)-p->GetZ(i2));
}
//______________________________________________________________________
Double_t AliITSv11::RminFrom2Points(TGeoPcon *p,Int_t i1,Int_t i2,Double_t z){
    // Retruns the value of Rmin corresponding to point z alone the line
    // defined by the two points p->GetRmin(i1),p->GetZ(i1) and p->GetRmin(i2),
    // p->GetZ(i2).

    return p->GetRmin(i2)+(p->GetRmin(i1)-p->GetRmin(i2))*(z-p->GetZ(i2))/
	(p->GetZ(i1)-p->GetZ(i2));
}
//______________________________________________________________________
Double_t AliITSv11::RFrom2Points(Double_t *p,Double_t *Z,Int_t i1,Int_t i2,Double_t z){
    // Retruns the value of Rmin corresponding to point z alone the line
    // defined by the two points p->GetRmin(i1),p->GetZ(i1) and p->GetRmin(i2),
    // p->GetZ(i2).

    return p[i2]+(p[i1]-p[i2])*(z-Z[i2])/(Z[i1]-Z[i2]);
}
//______________________________________________________________________
Double_t AliITSv11::Zfrom2MinPoints(TGeoPcon *p,Int_t i1,Int_t i2,Double_t r){
    // Retruns the value of Z corresponding to point R alone the line
    // defined by the two points p->GetRmin(i1),p->GetZ(i1) and 
    // p->GetRmin(i2),p->GetZ(i2)

    return p->GetZ(i2)+(p->GetZ(i1)-p->GetZ(i2))*(r-p->GetRmin(i2))/
	(p->GetRmin(i1)-p->GetRmin(i2));
}
//______________________________________________________________________
Double_t AliITSv11::Zfrom2MaxPoints(TGeoPcon *p,Int_t i1,Int_t i2,Double_t r){
    // Retruns the value of Z corresponding to point R alone the line
    // defined by the two points p->GetRmax(i1),p->GetZ(i1) and 
    // p->GetRmax(i2),p->GetZ(i2)

    return p->GetZ(i2)+(p->GetZ(i1)-p->GetZ(i2))*(r-p->GetRmax(i2))/
	(p->GetRmax(i1)-p->GetRmax(i2));
}
//______________________________________________________________________
Double_t AliITSv11::Zfrom2Points(Double_t *Z,Double_t *p,Int_t i1,Int_t i2,Double_t r){
    // Retruns the value of Z corresponding to point R alone the line
    // defined by the two points p->GetRmax(i1),p->GetZ(i1) and 
    // p->GetRmax(i2),p->GetZ(i2)

    return Z[i2]+(Z[i1]-Z[i2])*(r-p[i2])/(p[i1]-p[i2]);
}
//______________________________________________________________________
Double_t AliITSv11::RmaxFromZpCone(TGeoPcon *p,Double_t tc,Double_t z,Double_t th){
    // General SSD Outer Cone surface equation Rmax.
    Double_t tantc = TMath::Tan(tc*TMath::DegToRad());
    Double_t costc = TMath::Cos(tc*TMath::DegToRad());

    return -tantc*(z-p->GetZ(4))+p->GetRmax(4)+th/costc;
}
//______________________________________________________________________
Double_t AliITSv11::RmaxFromZpCone(Double_t *GetRmax,Double_t *GetZ,Double_t tc,Double_t z,Double_t th){
    // General SSD Outer Cone surface equation Rmax.
    Double_t tantc = TMath::Tan(tc*TMath::DegToRad());
    Double_t costc = TMath::Cos(tc*TMath::DegToRad());

    return -tantc*(z-GetZ[4])+GetRmax[4]+th/costc;
}
//______________________________________________________________________
Double_t AliITSv11::RminFromZpCone(TGeoPcon *p,Double_t tc,Double_t z,Double_t th){
    // General SSD Inner Cone surface equation Rmin.
    Double_t tantc = TMath::Tan(tc*TMath::DegToRad());
    Double_t costc = TMath::Cos(tc*TMath::DegToRad());

    return -tantc*(z-p->GetZ(3))+p->GetRmin(3)+th/costc;
}
//______________________________________________________________________
Double_t AliITSv11::RminFromZpCone(Double_t *GetRmin,Double_t *GetZ,Double_t tc,Double_t z,Double_t th){
    // General SSD Inner Cone surface equation Rmin.
    Double_t tantc = TMath::Tan(tc*TMath::DegToRad());
    Double_t costc = TMath::Cos(tc*TMath::DegToRad());

    return -tantc*(z-GetZ[3])+GetRmin[3]+th/costc;
}
//______________________________________________________________________
Double_t AliITSv11::ZFromRmaxpCone(TGeoPcon *p,Double_t tc,Double_t r,Double_t th){
    // General SSD Outer cone Surface equation for z.
    Double_t tantc = TMath::Tan(tc*TMath::DegToRad());
    Double_t costc = TMath::Cos(tc*TMath::DegToRad());

    return p->GetZ(4)+(p->GetRmax(4)+th/costc-r)/tantc;
}
//______________________________________________________________________
Double_t AliITSv11::ZFromRmaxpCone(Double_t *GetRmax,Double_t *GetZ,Double_t tc,Double_t r,Double_t th){
    // General SSD Outer cone Surface equation for z.
    Double_t tantc = TMath::Tan(tc*TMath::DegToRad());
    Double_t costc = TMath::Cos(tc*TMath::DegToRad());

    return GetZ[4]+(GetRmax[4]+th/costc-r)/tantc;
}
//______________________________________________________________________
Double_t AliITSv11::ZFromRminpCone(TGeoPcon *p,Double_t tc,Double_t r,Double_t th){
    // General SSD Inner cone Surface equation for z.
    Double_t tantc = TMath::Tan(tc*TMath::DegToRad());
    Double_t costc = TMath::Cos(tc*TMath::DegToRad());

    return p->GetZ(3)+(p->GetRmin(3)+th/costc-r)/tantc;
}
//______________________________________________________________________
void AliITSv11::RadiusOfCurvature(Double_t rc,Double_t theta0,Double_t z0,
		       Double_t r0,Double_t theta1,Double_t &z1,
		       Double_t &r1){
    // Given a initial point z0,r0, the initial angle theta0, and the radius
    // of curvature, returns the point z1, r1 at the angle theta1. Theta
    // measured from the r axis in the clock wise direction [degrees].
    Double_t sin0 = TMath::Sin(theta0*TMath::DegToRad());
    Double_t cos0 = TMath::Cos(theta0*TMath::DegToRad());
    Double_t sin1 = TMath::Sin(theta1*TMath::DegToRad());
    Double_t cos1 = TMath::Cos(theta1*TMath::DegToRad());

    z1 = rc*(sin1-sin0)+z0;
    r1 = rc*(cos1-cos0)+r0;
    return;
}
//______________________________________________________________________
void AliITSv11::SPDCone(TGeoVolume *Moth){
    // Define the detail SPD support cone geometry.
    // Inputs:
    //   none.
    // Outputs:
    //  none.
    // Return:
    //  none.
}
//______________________________________________________________________
void AliITSv11::SDDCone(TGeoVolume *Moth){
    // Define the detail SDD support cone geometry.
    // Inputs:
    //   none.
    // Outputs:
    //  none.
    // Return:
    //  none.
} 
//______________________________________________________________________
void AliITSv11::SSDCone(TGeoVolume *Moth){
    // Define the detail SSD support cone geometry.
    // Inputs:
    //   none.
    // Outputs:
    //  none.
    // Return:
    //  none.
    const Double_t Thickness = 13.0*kmm; // Thickness of Rohacell+carbon fiber
    const Double_t Cthick    = 1.5*kmm; // Carbon finber thickness
    const Double_t Rcurv     = 15.0*kmm; // Radius of curvature.
    const Double_t Tc        = 51.0; // angle of SSD cone [degrees].
    const Double_t Sintc = TMath::Sin(Tc*TMath::DegToRad());
    const Double_t Costc = TMath::Cos(Tc*TMath::DegToRad());
    const Double_t Tantc = TMath::Tan(Tc*TMath::DegToRad());
    const Double_t ZouterMilled = (13.5-5.0)*kmm;
    const Double_t Zcylinder    = 170.0*kmm;
    const Double_t Z0           = Zcylinder + 100.0*kmm;
    const Int_t Nspoaks         = 12;
    const Int_t Nmounts         = 4;
    const Double_t DmountAngle  = 9.0; // degrees
    const Double_t RoutMax      = 0.5*985.0*kmm;
    const Double_t RoutHole     = 0.5*965.0*kmm;
    const Double_t RoutMin      = 0.5*945.0*kmm;
    const Double_t RholeMax     = 0.5*890.0*kmm;
    const Double_t RholeMin     = 0.5*740.0*kmm;
    const Double_t RpostMin     = 316.0*kmm;
    const Double_t ZpostMax     = 196.0*kmm;
    const Int_t Nposts          = 6;
    const Double_t Phi0Post     = 0.0; // degree
    const Double_t dRpost       = 23.0*kmm;
    const Double_t RinMax       = 0.5*590.0*kmm;
    const Double_t RinCylinder  = 0.5*597.0*kmm;
    const Double_t RinHole      = 0.5*575.0*kmm;
    const Double_t RinMin       = 0.5*562.0*kmm;
    const Double_t dZin         = 15.0*kmm;
    // SSD-SDD Thermal/Mechanical cylinder mounts
    const Int_t NinScrews          = 40;
    const Double_t Phi0Screws      = 0.5*360.0/((const Double_t)NinScrews);//d
    const Double_t RcylinderScrews = 0.5*570.0*kmm;//from older drawing????
    const Double_t DscrewHead      = 8.0*kmm;
    const Double_t DscrewShaft     = 4.6*kmm;
    const Double_t ThScrewHeadHole = 8.5*kmm;
    // SDD mounting bracket, SSD part
    const Double_t NssdSupports      = 3;// mounting of U and T
    const Double_t DssdsddBracketAngle = 9.0; // degrees
    const Double_t Phi0SDDsupports   = 0.0; // degree
    const Double_t RsddSupportPlate  = 0.5*585.0*kmm;
    const Double_t ThSDDsupportPlate = 4.0*kmm;
    const Double_t WsddSupportPlate  = 70.0*kmm;
    TGeoMedium *SSDcf  = 0; // SSD support cone Carbon Fiber materal number.
    TGeoMedium *SSDfs  = 0; // SSD support cone inserto stesalite 4411w.
    TGeoMedium *SSDfo  = 0; // SSD support cone foam, Rohacell 50A.
    TGeoMedium *SSDss  = 0; // SSD support cone screw material,Stainless steal
    TGeoMedium *SSDair = 0; // SSD support cone Air
    TGeoMedium *SSDal  = 0; // SSD support cone SDD mounting bracket Al

    // Lets start with the upper left outer carbon fiber surface.
    // Between za[2],rmaxa[2] and za[4],rmaxa[4] there is a curved section
    // given by rmaxa = rmaxa[2]-r*Sind(t) for 0<=t<=Tc and 
    // za = za[2] + r*Cosd(t) for 0<=t<=Tc. Simularly between za[1],rmina[1
    // and za[3],rmina[3] there is a curve section given by
    // rmina = rmina[1]-r*Sind(t) for 0<=t<=Tc and za = za[1]+r&Sind(t)
    // for t<=0<=Tc. These curves have been replaced by straight lines
    // between the equivelent points for simplicity.
    Double_t dza = Thickness/Sintc-(RoutMax-RoutMin)/Tantc;
    Int_t i,j;
    Double_t z[9],rn[9],rx[9],phi,dphi;
    Double_t t,t0,Z,Rmin,Rmax; // Temp variables.
    if(dza<=0){ // The number or order of the points are in error for a proper
	// call to pcons!
	Error("SSDcone","The definition of the points for a call to PCONS is"
	      " in error. abort.");
	return;
    } // end if
    // Poly-cone Volume A. Top part of SSD cone Carbon Fiber.
    phi   = 0.0;
    dphi  = 360.0;
    z[0]  = 0.0;
    rn[0] = RoutMin;
    rx[0] = RoutMax;
    z[1]  = z[0]+ZouterMilled - dza; // za[2] - dza.
    rn[1] = rn[0];
    rx[1] = rx[0];
    z[2]  = z[0]+ZouterMilled;//From Drawing ALR-0767 and ALR-0767/3
    rx[2] = rx[0];
    RadiusOfCurvature(Rcurv,0.,z[1],rn[1],Tc,z[3],rn[3]);
    rn[2] = RFrom2Points(rn,z,3,1,z[2]);
    RadiusOfCurvature(Rcurv,0.,z[2],rx[2],Tc,z[4],rx[4]);
    rn[4] = RminFromZpCone(rn,z,Tc,z[4]);
    rx[3] = RFrom2Points(rx,z,4,2,z[3]);
    rn[5] = RholeMax;
    z[5]  = Zfrom2Points(z,rn,4,3,rn[5]);
    rx[5] = RmaxFromZpCone(rx,z,Tc,z[5]);
    rn[6] = RholeMax;
    rx[6] = rn[6];
    z[6]  = ZFromRmaxpCone(z,rx,Tc,rx[6]);
    TGeoPcon *A = new TGeoPcon("ITS SSD Suport cone Carbon Fiber "
				   "Surface outer left",phi,dphi,7);
    for(i=0;i<A->GetNz();i++){
	if(fDebug) cout<<i<<"A: z="<<z[i]<<" Rmin="<<rn[i]<<" Rmax="<<rx[i]<<endl;
	A->DefineSection(i,z[i],rn[i],rx[i]);
    } // end for i
    //
    // Poly-cone Volume B. Stesalite inside volume A.
    // Now lets define the Inserto Stesalite 4411w material volume.
    phi   = 0.0;
    dphi  = 360.0;
    z[0]  = A->GetZ(0);
    rn[0] = A->GetRmin(0)+Cthick;
    rx[0] = A->GetRmax(0)-Cthick;
    z[1]  = A->GetZ(1);
    rn[1] = rn[0];
    rx[1] = rx[0];
    z[2]  = A->GetZ(2);
    rx[2] = rx[1];
    RadiusOfCurvature(Rcurv-Cthick,0.,z[2],rx[2],Tc,z[3],rx[3]);
    RadiusOfCurvature(Rcurv+Cthick,0.,z[1],rn[1],Tc,z[4],rn[4]);
    rn[2] = RFrom2Points(rn,z,4,1,z[2]);
    rn[3] = RFrom2Points(rn,z,4,1,z[3]);
    z[5]  = z[4]+(Thickness-2.0*Cthick)/Sintc;
    rn[5] = RmaxFromZpCone(A,Tc,z[5],-Cthick);
    rx[5] = rn[5];
    rx[4] = RFrom2Points(rx,z,5,3,z[4]);
    TGeoPcon *B = new TGeoPcon("ITS SSD Suport cone Inserto Stesalite "
				   "left edge",phi,dphi,6);
    for(i=0;i<B->GetNz();i++){
	if(fDebug) cout<<i<<"B: z="<<z[i]<<" Rmin="<<rn[i]<<" Rmax="<<rx[i]<<endl;
	B->DefineSection(i,z[i],rn[i],rx[i]);
    } // end for i
    //
    // Poly-cone Volume C. Foam inside volume A.
    // Now lets define the Rohacell foam material volume.
    phi   = 0.0;
    dphi  = 360.0;
    z[0]  = B->GetZ(4);
    rn[0] = B->GetRmin(4);
    rx[0] = rn[0];
    z[1]  = B->GetZ(5);
    rx[1] = B->GetRmin(5);
    rn[2] = A->GetRmin(5)+Cthick;//space for carbon fiber covering hole
    z[2]  = ZFromRminpCone(A,Tc,rn[2],+Cthick);
    rn[1] = RFrom2Points(rn,z,2,0,z[1]);
    rx[3] = A->GetRmin(6)+Cthick;
    rn[3] = rx[3];
    z[3]  = ZFromRmaxpCone(A,Tc,rx[3],-Cthick);
    rx[2] = RFrom2Points(rx,z,3,1,z[2]);
    TGeoPcon *C = new TGeoPcon("ITS SSD Suport cone Rohacell foam "
				   "left edge",phi,dphi,4);
    for(i=0;i<C->GetNz();i++){
	if(fDebug) cout<<i<<"C: z="<<z[i]<<" Rmin="<<rn[i]<<" Rmax="<<rx[i]<<endl;
	C->DefineSection(i,z[i],rn[i],rx[i]);
    } // end for i
    //
    // In volume SCB, th Inserto Stesalite 4411w material volume, there
    // are a number of Stainless steel screw and pin studs which will be
    // filled with screws/studs.
    rn[0] = 0.0,rx[0] = 6.0,z[0] = 0.5*10.0; // mm
    TGeoTube *D = new TGeoTube("ITS Screw+stud used to mount things to "
				   "the SSD support cone",rn[0],rx[0],z[0]);
    rn[0] = 0.0;rx[0] = 6.0;z[0] = 0.5*12.0; // mm
    TGeoTube *E = new TGeoTube("ITS pin used to mount things to the "
				   "SSD support cone",rn[0],rx[0],z[0]);
    //
    // Poly-cone Volume F. Foam in spoak reagion, inside volume A.
    // There is no carbon fiber between this upper left section and the
    // SSD spoaks. We remove it by replacing it with Rohacell foam.
    t = Cthick/(0.5*(RholeMax+RholeMin));// It is not posible to get
    // the carbon fiber thickness uniform in this phi direction. We can only
    // make it a fixed angular thickness.
    t *= 180.0/TMath::Pi();
    phi  = 12.5+t; // degrees see drawing ALR-0767.
    dphi  = 5.0 - 2.0*t; // degrees
    z[0]  = C->GetZ(2);
    rn[0] = C->GetRmin(3);
    rx[0] = rn[0];
    rn[1] = A->GetRmin(5);
    rx[1] = rn[0];
    z[1]  = ZFromRminpCone(A,Tc,rn[1],+Cthick);
    z[2]  = C->GetZ(3);
    rn[2] = rn[1];
    rx[2] = rx[1];
    rn[3] = A->GetRmin(6);
    rx[3] = rn[3];
    z[3]  = ZFromRmaxpCone(A,Tc,rx[3],-Cthick);
    TGeoPcon *F = new TGeoPcon("ITS SSD Top Suport cone Rohacell foam "
				   "Spoak",phi,dphi,4);
    for(i=0;i<F->GetNz();i++){
	if(fDebug) cout<<i<<"F: z="<<z[i]<<" Rmin="<<rn[i]<<" Rmax="<<rx[i]<<endl;
	F->DefineSection(i,z[i],rn[i],rx[i]);
    } // end for i
    //=================================================================
     // Poly-cone Volume G.
    // Now for the spoak part of the SSD cone.
    // It is not posible to inclue the radius of curvature between
    // the spoak part and the upper left part of the SSD cone or lowwer right
    // part. This would be discribed by the following curves.
    // R = Rmax - (5mm)*Sin(t) phi = phi0+(5mm*180/(Pi*RoutHole))*Sin(t) 
    // where 0<=t<=90 For the inner curve a simular equiation holds.
    phi   = 12.5; // degrees see drawing ALR-0767.
    dphi  = 5.0; // degrees
    z[0]  = A->GetZ(5);
    rn[0] = A->GetRmin(5);
    rx[0] = rn[0];
    z[1]  = A->GetZ(6);
    rn[1] = RminFromZpCone(A,Tc,z[1]);
    rx[1] = rx[0];
    rn[2] = RholeMin;
    z[2]  = ZFromRminpCone(A,Tc,rn[2]);
    rx[2] = RmaxFromZpCone(A,Tc,z[2]);
    rn[3] = rn[2];
    rx[3] = rn[3];
    z[3]  = ZFromRmaxpCone(A,Tc,rx[3]);
    TGeoPcon *G = new TGeoPcon("ITS SSD spoak carbon fiber surfaces",
				   phi,dphi,4);
    for(i=0;i<G->GetNz();i++){
	if(fDebug) cout<<i<<"G: z="<<z[i]<<" Rmin="<<rn[i]<<" Rmax="<<rx[i]<<endl;
	G->DefineSection(i,z[i],rn[i],rx[i]);
    } // end for i
    // For the foam core.
    // Poly-cone Volume H.
    t = Cthick/(0.5*(RholeMax+RholeMin));// It is not posible to get the
    // carbon fiber thickness uniform in this phi direction. We can only
    // make it a fixed angular thickness.
    t *= 180.0/TMath::Pi();
    phi   = 12.5+t; // degrees
    dphi  = 5.0 - 2.0*t; // degrees see drawing ALR-0767.
    z[0]  = F->GetZ(1);
    rn[0] = G->GetRmin(0);
    rx[0] = rn[0];
    z[1]  = F->GetZ(3);
    rn[1] = RminFromZpCone(A,Tc,z[1],+Cthick);
    rx[1] = rx[0];
    z[2]  = ZFromRminpCone(A,Tc,G->GetRmin(2),+Cthick);
    rn[2] = G->GetRmin(2);
    rx[2] = RmaxFromZpCone(A,Tc,z[2],-Cthick);
    z[3]  = ZFromRmaxpCone(A,Tc,G->GetRmin(3),-Cthick);
    rn[3] = G->GetRmin(3);
    rx[3] = rn[3];
    TGeoPcon *H = new TGeoPcon("ITS SSD support cone Rohacell foam Spoak",
				   phi,dphi,4); 
    for(i=0;i<H->GetNz();i++){
	if(fDebug) cout<<i<<"H: z="<<z[i]<<" Rmin="<<rn[i]<<" Rmax="<<rx[i]<<endl;
	H->DefineSection(i,z[i],rn[i],rx[i]);
    } // end for i
    //
    //==================================================================
    // Now for the Inner most part of the SSD cone.
    //Poly-cone Volume I.
    phi   = 0.0;
    dphi  = 360.0;
    z[0]  = G->GetZ(2);
    rn[0] = G->GetRmin(2);
    rx[0] = rn[0];
    z[1]  = G->GetZ(3);
    rn[1] = RminFromZpCone(A,Tc,z[1]);
    rx[1] = rx[0];
    rn[4] = RinMin;
    rn[5] = RinMin;
    RadiusOfCurvature(Rcurv,90.0,0.0,RinMax,90.0-Tc,Z,rx[5]); // z dummy
    z[5]  = ZFromRmaxpCone(A,Tc,rx[5]);
    z[6]  = Zcylinder;
    rn[6] = RinMin;
    z[7]  = z[6];
    rn[7] = RinCylinder;
    rn[8] = RinCylinder;
    rx[8] = rn[8];
    Rmin   = rn[5];
    RadiusOfCurvature(Rcurv,90.0-Tc,z[5],rx[5],90.0,Z,Rmax);
    Rmax   = RinMax;
    z[8]  = Z+(z[5]-Z)*(rx[8]-Rmax)/(rx[5]-Rmax);
    rx[6] = RFrom2Points(rx,z,8,5,z[6]);
    rx[7] = rx[6];
    z[3]  = Z-dZin;
    z[4]  = z[3];
    rx[3] = RmaxFromZpCone(A,Tc,z[3]);
    rx[4] = rx[3];
    //rmin dummy
    RadiusOfCurvature(Rcurv,90.,z[3],0.,90.-Tc,z[2],Rmin);
    rn[2] = RminFromZpCone(A,Tc,z[2]);
    rx[2] = RmaxFromZpCone(A,Tc,z[2]);
    // z dummy
    RadiusOfCurvature(Rcurv,90.-Tc,0.0,rn[2],90.0,Z,rn[3]);
    TGeoPcon *I = new TGeoPcon("ITS SSD lower/inner right part of SSD "
				   "cone",phi,dphi,9);
    for(i=0;i<I->GetNz();i++){
	if(fDebug) cout<<i<<"I: z="<<z[i]<<" Rmin="<<rn[i]<<" Rmax="<<rx[i]<<endl;
	I->DefineSection(i,z[i],rn[i],rx[i]);
    } // end for i
    // Now for Inserto volume at the inner most radius.
    // Poly-cone Volume K.
    phi   = 0.0;
    dphi  = 360.0;
    z[1]  = I->GetZ(3)+Cthick;
    rn[1] = I->GetRmin(3);
    z[2]  = z[1];
    rn[2] = I->GetRmin(4);
    rn[3] = rn[2];
    rn[4] = rn[2];
    rx[4] = I->GetRmax(5)-Cthick*Sintc;
    RadiusOfCurvature(Rcurv+Cthick,90.0,z[1],rn[1],90.0-Tc,z[0],rn[0]);
    rx[0] = rn[0];
    z[3]  = z[0]+(Thickness-2.0*Cthick)*Costc;;
    rx[3] = rx[0]+(Thickness-2.0*Cthick)*Sintc;
    rx[1] = RFrom2Points(rx,z,3,0,z[1]);
    rx[2] = rx[1];
    z[4]  = ZFromRmaxpCone(A,Tc,rx[4],-Cthick);
    rn[5] = rn[2];
    z[5]  = I->GetZ(6);
    rx[5] = (I->GetRmax(5)-I->GetRmax(8))/(I->GetZ(5)-I->GetZ(8))*(z[5]-z[4])+
	     rx[4];
    TGeoPcon *K = new TGeoPcon("ITS SSD inner most inserto material",
				   phi,dphi,6);
    for(i=0;i<K->GetNz();i++){
	if(fDebug) cout<<i<<"K: z="<<z[i]<<" Rmin="<<rn[i]<<" Rmax="<<rx[i]<<endl;
	K->DefineSection(i,z[i],rn[i],rx[i]);
    } // end for i
    // Now for foam core at the inner most radius.
    // Poly-cone Volume J.
    phi   = 0.0;
    dphi  = 360.0;
    rn[0] = I->GetRmin(0)-Cthick;
    z[0]  = ZFromRminpCone(A,Tc,rn[0],+Cthick);
    rx[0] = rn[0];
    rx[1] = rx[0];
    z[1]  = ZFromRmaxpCone(A,Tc,rx[1],-Cthick);
    rn[1] = RminFromZpCone(A,Tc,z[1],+Cthick);
    z[2]  = K->GetZ(0);
    rn[2] = K->GetRmin(0);
    rx[2] = RmaxFromZpCone(A,Tc,z[2],-Cthick);
    z[3]  = K->GetZ(3);
    rn[3] = K->GetRmax(3);
    rx[3] = rn[3];
    TGeoPcon *J = new TGeoPcon("ITS SSD inner most foam core",phi,dphi,4); 
    for(i=0;i<J->GetNz();i++){
	if(fDebug) cout<<i<<"J: z="<<z[i]<<" Rmin="<<rn[i]<<" Rmax="<<rx[i]<<endl;
	J->DefineSection(i,z[i],rn[i],rx[i]);
    } // end for i
    // Now for foam core at the top of the inner most radius where 
    // the spoaks are.
    t = Cthick/(0.5*(RholeMax+RholeMin));// It is not posible to get the
    // carbon fiber thickness uniform in this phi direction. We can only
    // make it a fixed angular thickness.
    // Poly-cone Volume L.
    t *= 180.0/TMath::Pi();
    phi   = 12.5+t; // degrees
    dphi  = 5.0 - 2.0*t; // degrees see drawing ALR-0767.
    z[0]  = H->GetZ(2);
    rn[0] = H->GetRmin(2);
    rx[0] = rn[0];
    z[1]  = J->GetZ(0);
    rn[1] = J->GetRmin(0);
    rx[1] = I->GetRmax(1);
    z[2]  = H->GetZ(3);
    rn[2] = rn[1];
    rx[2] = rx[1];
    z[3]  = J->GetZ(1);
    rn[3] = rn[2];
    rx[3] = rn[3];
    TGeoPcon *L = new TGeoPcon("ITS SSD Bottom cone Rohacell foam Spoak",
				   phi,dphi,4);
    for(i=0;i<L->GetNz();i++){
	if(fDebug) cout<<i<<"L: z="<<z[i]<<" Rmin="<<rn[i]<<" Rmax="<<rx[i]<<endl;
	L->DefineSection(i,z[i],rn[i],rx[i]);
    } // end for i
    // Now for the SSD mounting posts
    // Poly-cone Volume O.
    dphi  = 180.0*dRpost/(RpostMin+0.5*dRpost)/TMath::Pi(); //
    phi   = Phi0Post-0.5*dphi; // degrees
    rn[0] = RpostMin+dRpost;
    rx[0] = rn[0];
    z[0]  = ZFromRmaxpCone(A,Tc,rx[0]);
    rn[1] = RpostMin;
    z[1]  = ZFromRmaxpCone(A,Tc,rn[1]);
    rx[1] = rx[0];
    z[2]  = ZpostMax;
    rn[2] = RpostMin;
    rx[2] = rn[2]+dRpost;
    TGeoPcon *O = new TGeoPcon("ITS SSD mounting post, carbon fiber",
				   phi,dphi,3);
    for(i=0;i<O->GetNz();i++){
	if(fDebug) cout<<i<<"O: z="<<z[i]<<" Rmin="<<rn[i]<<" Rmax="<<rx[i]<<endl;
	O->DefineSection(i,z[i],rn[i],rx[i]);
    } // end for i
    // Now for the SSD mounting posts
    // Poly-cone Volume P.
    t = 180.0*Cthick/(RpostMin+0.5*dRpost)/TMath::Pi();
    dphi  = O->GetDphi()-2.0*t; // degrees
    phi   = O->GetPhi1()+t; //
    rn[0] = O->GetRmin(0)-Cthick;
    rx[0] = rn[0];
    z[0]  = ZFromRmaxpCone(A,Tc,rx[0]);
    rn[1] = O->GetRmin(1)+Cthick;
    rx[1] = O->GetRmin(0)-Cthick;
    z[1]  = ZFromRmaxpCone(A,Tc,rn[1]);
    rn[2] = rn[1];
    rx[2] = rx[1];
    z[2]  = ZpostMax;
    TGeoPcon *P = new TGeoPcon("ITS SSD mounting post, Inserto",
				   phi,dphi,3);
    for(i=0;i<P->GetNz();i++){
	if(fDebug) cout<<i<<"P: z="<<z[i]<<" Rmin="<<rn[i]<<" Rmax="<<rx[i]<<endl;
	P->DefineSection(i,z[i],rn[i],rx[i]);
    } // end for i
    // This insrto continues into the SSD cone displacing the foam
    // and the carbon fiber surface at those points where the posts are.
    //Poly-cone Vol. M
    phi   = P->GetPhi1();
    dphi  = P->GetDphi();
    rn[0] = RpostMin+dRpost-Cthick;
    rx[0] = rn[0];
    z[0]  = ZFromRminpCone(A,Tc,rn[0],+Cthick);
    rx[1] = rx[0];
    z[1]  = ZFromRmaxpCone(A,Tc,rx[1],-Cthick);
    rn[1] = RminFromZpCone(A,Tc,z[1],+Cthick);
    rn[2] = RpostMin+Cthick;
    z[2]  = ZFromRminpCone(A,Tc,rn[2],+Cthick);
    rx[2] = RmaxFromZpCone(A,Tc,z[2],-Cthick);
    rn[3] = rn[2];
    rx[3] = rn[3];
    z[3]  = ZFromRmaxpCone(A,Tc,rx[3],-Cthick);
    TGeoPcon *M = new TGeoPcon("ITS SSD mounting post foam substitute, "
				   "Inserto",phi,dphi,4);
    for(i=0;i<M->GetNz();i++){
	if(fDebug) cout<<i<<"M: z="<<z[i]<<" Rmin="<<rn[i]<<" Rmax="<<rx[i]<<endl;
	M->DefineSection(i,z[i],rn[i],rx[i]);
    } // end for i
    //
    //Poly-cone Vol. N
    phi   = P->GetPhi1();
    dphi  = P->GetDphi();
    z[0]  = M->GetZ(1);
    rn[0] = M->GetRmax(1);
    rx[0] = rn[0];
    rx[1] = rx[0];
    z[1]  = ZFromRmaxpCone(A,Tc,rx[1]);
    rn[1] = RmaxFromZpCone(A,Tc,z[1],-Cthick);
    z[2]  = M->GetZ(3);
    rn[2] = M->GetRmin(3);
    rx[2] = RmaxFromZpCone(A,Tc,z[2]);
    rn[3] = rn[2];
    rx[3] = rn[3];
    z[3]  = ZFromRmaxpCone(A,Tc,rx[3]);
    TGeoPcon *N = new TGeoPcon("ITS SSD mounting post CF subsititute, "
				   "Inserto",phi,dphi,4);
    for(i=0;i<N->GetNz();i++){ 
	if(fDebug) cout<<i<<"N: z="<<z[i]<<" Rmin="<<rn[i]<<" Rmax="<<rx[i]<<endl;
	N->DefineSection(i,z[i],rn[i],rx[i]);
    } // end for i
    // Bolt heads holding the SSD-SDD tube to the SSD cone.
    // Bolt -- PolyCone
    //Poly-cone Volume Q.
    phi   = 0.0;
    dphi  = 360.0;
    z[0]  = I->GetZ(4)-ThSDDsupportPlate;
    rn[0] = 0.0;
    rx[0] = 0.5*DscrewHead;
    z[1]  = I->GetZ(4)-ThScrewHeadHole;
    rn[1] = 0.0;
    rx[1] = 0.5*DscrewHead;
    z[2]  = z[1];
    rn[2] = 0.0;
    rx[2] = 0.5*DscrewShaft;
    z[3]  = z[2];
    rn[3] = 0.0;
    rx[3] = rx[2];
    TGeoPcon *Q = new TGeoPcon("ITS SSD Thermal sheal stainless steel "
				   "bolts",phi,dphi,4);
    for(i=0;i<Q->GetNz();i++){
	if(fDebug) cout<<i<<"Q: z="<<z[i]<<" Rmin="<<rn[i]<<" Rmax="<<rx[i]<<endl;
	Q->DefineSection(i,z[i],rn[i],rx[i]);
    } // end for i
    // air infront of bolt (stasolit Volume K) -- Tube
    z[0]  = 0.5*(Thickness-ThScrewHeadHole);
    rn[0] = 0.0;
    rx[0] = 0.5*DscrewHead;
    TGeoTube *R = new TGeoTube("ITS Air in front of bolt (in stasolit)",
				   rn[0],rx[0],z[0]);
    // air infront of bolt (carbon fiber volume I) -- Tube
    z[0]  = 0.5*Thickness;
    rn[0] = 0.0;
    rx[0] = R->GetRmax();
    TGeoTube *S = new TGeoTube("ITS Air in front of Stainless Steal "
				   "Screw end, N6",rn[0],rx[0],z[0]);
    // SDD support plate, SSD side.
    //Poly-cone Volume T.
    dphi  = 180.0*WsddSupportPlate/(RsddSupportPlate*TMath::Pi());
    phi   = Phi0SDDsupports-0.5*dphi;
    z[0]  = K->GetZ(2);
    rn[0] = I->GetRmin(4);
    rx[0] = RsddSupportPlate;
    z[1]  = I->GetZ(4) - ThSDDsupportPlate;
    rn[1] = rn[0];
    rx[1] = rx[0];
    TGeoPcon *T = new TGeoPcon("ITS SSD-SDD mounting bracket Inserto->Al.",
				   phi,dphi,2);
    for(i=0;i<T->GetNz();i++){
	if(fDebug) cout<<i<<"T: z="<<z[i]<<" Rmin="<<rn[i]<<" Rmax="<<rx[i]<<endl;
	T->DefineSection(i,z[i],rn[i],rx[i]);
    } // end for i
    //
    // Poly-cone Volume U.
    TGeoPcon *U;
    if(I->GetRmin(3)<T->GetRmax(0)){
	dphi  = T->GetDphi();
	phi   = T->GetPhi1();
	z[2]  = I->GetZ(4);
	rn[2] = T->GetRmin(0);
	rx[2] = T->GetRmax(0);
	z[3]  = K->GetZ(2);
	rn[3] = rn[2];
	rx[3] = rx[2];
	z[1]  = z[2];
	rn[1] = I->GetRmin(3);
	rx[1] = rx[3];
	rx[0] = T->GetRmax(0);
	rn[0] = rx[0];
	z[0]  = Zfrom2MinPoints(I,2,3,rn[0]);
	U = new TGeoPcon("ITS SSD-SDD mounting bracket CF->Al.",phi,dphi,4);
    }else{
	dphi  = T->GetDphi();
	phi   = T->GetPhi1();
	z[0]  = I->GetZ(4);
	rn[0] = T->GetRmin(0);
	rx[0] = T->GetRmax(0);
	z[1]  = K->GetZ(2);
	rn[1] = rn[0];
	rx[1] = rx[0];
	U = new TGeoPcon("ITS SSD-SDD mounting bracket CF->Al.",phi,dphi,2);
    }// end if
    for(i=0;i<U->GetNz();i++){
	if(fDebug) cout<<i<<"U: z="<<z[i]<<" Rmin="<<rn[i]<<" Rmax="<<rx[i]<<endl;
	U->DefineSection(i,z[i],rn[i],rx[i]);
    } // end for i
    //
    TGeoManager *mgr = gGeoManager;
    SSDcf = mgr->GetMedium("ITSssdCarbonFiber");
    SSDfs = mgr->GetMedium("ITSssdStaselite4411w");
    SSDfo = mgr->GetMedium("ITSssdRohacell50A");
    SSDss = mgr->GetMedium("ITSssdStainlessSteal");
    SSDair= mgr->GetMedium("ITSssdAir");
    SSDal = mgr->GetMedium("ITSssdAl");
    TGeoVolume *Av,*Bv,*Cv,*Dv,*Ev,*Fv,*Gv,*Hv,*Iv,*Jv,*Kv,*Lv,*Mv,*Nv,
	       *Ov,*Pv,*Qv,*Rv,*Sv,*Tv,*Uv;
    Av = new TGeoVolume("ITSssdConeA",A,SSDcf);
    mgr->AddVolume(Av);
    Av->SetLineColor(1);
    Av->SetLineWidth(1);
    Bv = new TGeoVolume("ITSssdConeB",B,SSDfs);
    mgr->AddVolume(Bv);
    Cv = new TGeoVolume("ITSssdConeC",C,SSDfo);
    mgr->AddVolume(Cv);
    Dv = new TGeoVolume("ITSssdConeD",D,SSDss);
    mgr->AddVolume(Dv);
    Ev = new TGeoVolume("ITSssdConeE",E,SSDss);
    mgr->AddVolume(Ev);
    Fv = new TGeoVolume("ITSssdConeF",F,SSDfo);
    mgr->AddVolume(Fv);
    Gv = new TGeoVolume("ITSssdConeG",G,SSDcf);
    mgr->AddVolume(Gv);
    Gv->SetLineColor(2);
    Gv->SetLineWidth(2);
    Hv = new TGeoVolume("ITSssdConeH",H,SSDfo);
    mgr->AddVolume(Hv);
    Iv = new TGeoVolume("ITSssdConeI",I,SSDcf);
    mgr->AddVolume(Iv);
    Iv->SetLineColor(3);
    Iv->SetLineWidth(3);
    Jv = new TGeoVolume("ITSssdConeJ",J,SSDfo);
    mgr->AddVolume(Jv);
    Kv = new TGeoVolume("ITSssdConeK",K,SSDfs);
    mgr->AddVolume(Kv);
    Lv = new TGeoVolume("ITSssdConeL",L,SSDfo);
    mgr->AddVolume(Lv);
    Mv = new TGeoVolume("ITSssdConeM",M,SSDfs);
    mgr->AddVolume(Mv);
    Nv = new TGeoVolume("ITSssdConeN",N,SSDfs);
    mgr->AddVolume(Nv);
    Ov = new TGeoVolume("ITSssdConeO",O,SSDcf);
    mgr->AddVolume(Ov);
    Iv->SetLineColor(4);
    Iv->SetLineWidth(4);
    Pv = new TGeoVolume("ITSssdConeP",P,SSDfs);
    mgr->AddVolume(Pv);
    Qv = new TGeoVolume("ITSssdConeQ",Q,SSDss);
    mgr->AddVolume(Qv);
    Rv = new TGeoVolume("ITSssdConeR",R,SSDair);
    mgr->AddVolume(Rv);
    Sv = new TGeoVolume("ITSssdConeS",S,SSDair);
    mgr->AddVolume(Sv);
    Tv = new TGeoVolume("ITSssdConeT",T,SSDal);
    mgr->AddVolume(Tv);
    Uv = new TGeoVolume("ITSssdConeU",U,SSDal);
    mgr->AddVolume(Uv);
    //
    TGeoTranslation *tran = new TGeoTranslation("ITSssdConeTrans",0.0,0.0,-Z0);
    TGeoRotation *rot180  = new TGeoRotation("ITSssdConeRot180",0.0,180.0,0.0);
    TGeoCombiTrans *flip  = new TGeoCombiTrans("ITSssdConeFlip",0.0,0.0,Z0,rot180);
    TGeoTranslation *tranR,*tranS;
    TGeoCombiTrans *fliptran,*rottran;
    TGeoRotation *rot,*zspoaks,*zspoaks180;
    Av->AddNode(Bv,1,0);
    Av->AddNode(Cv,1,0);
    Moth->AddNode(Av,1,tran); // RB24 side
    Moth->AddNode(Av,2,flip); // RB26 side (Absorber)
    Moth->AddNode(Iv,1,tran); // RB24 side
    Moth->AddNode(Iv,2,flip); // RB26 side (Absorber)
    Gv->AddNode(Hv,1,0);
    for(i=0;i<Nspoaks;i++){ // SSD Cone Spoaks
	zspoaks = new TGeoRotation("",0.0,0.0,
				   ((Double_t)i*360.)/((Double_t)Nspoaks));
	rottran = new TGeoCombiTrans("",0.0,0.0,-Z0,zspoaks);
	Moth->AddNode(Gv,i+1,rottran); // RB24 side
	Av->AddNode(Fv,i+1,zspoaks);
	Iv->AddNode(Lv,i+1,zspoaks);
	zspoaks180 =  new TGeoRotation("",0.0,180.0,
				       ((Double_t)i*360.)/((Double_t)Nspoaks));
	fliptran = new TGeoCombiTrans("",0.0,0.0,Z0,zspoaks180);
	Moth->AddNode(Gv,Nspoaks+i+1,fliptran); // RB26 side
    } // end for i
    Iv->AddNode(Jv,1,0);
    Iv->AddNode(Kv,1,0);
    Ov->AddNode(Pv,1,0);
    //Pv->AddNode(Qv,2,?); // Screw head
    //Pv->AddNode(Qv,3,?); // Screw head
    //Pv->AddNode(Vv,1,?); // Air hole in Posts
    //Pv->AddNode(Vv,2,?); // Air hole in Posts
    //Mv->AddNode(Wv,1,?); // Air hole in Posts
    //Mv->AddNode(Wv,2,?); // Air hole in Posts
    //Nv->AddNode(Xv,1,?); // Air hole in Posts
    //Nv->AddNode(Xv,2,?); // Air hole in Posts
    TGeoRotation *zposts,*zposts180;
    for(i=0;i<Nposts;i++){ // SSD Cone mounting posts
	zposts = new TGeoRotation("",0.0,0.0,
				  ((Double_t)i*360.)/((Double_t)Nposts));
	rottran = new TGeoCombiTrans("",0.0,0.0,-Z0,zposts);
	Moth->AddNode(Ov,i+1,rottran); // RB24 side
	Jv->AddNode(Mv,i+1,zposts);
	Iv->AddNode(Nv,i+1,zposts);
	//Jv->AddNode(Xv,2*i+3,?); // Air hole in Posts
	//Jv->AddNode(Xv,2*i+4,?); // Air hole in Posts
	zposts180 = new TGeoRotation("",0.0,180.0,
				  ((Double_t)i*360.)/((Double_t)Nposts));
	fliptran = new TGeoCombiTrans("",0.0,0.0,Z0,zposts180);
	Moth->AddNode(Ov,Nposts+i+1,fliptran); // RB26 side
    } // end for i
    //
    for(i=0;i<NinScrews;i++){
	t = Phi0Screws+360.*((Double_t)i)/((Double_t)NinScrews);
	t *= TMath::DegToRad();
	tran= new TGeoTranslation("",RcylinderScrews*TMath::Cos(t),
				  RcylinderScrews*TMath::Sin(t),0.0);
	Kv->AddNode(Qv,i+4,rottran);
	if(/*not where volumes U and T are*/kTRUE){
	    tranR = new TGeoTranslation("",RinHole*TMath::Cos(t),
					RinHole*TMath::Sin(t),
					K->GetZ(2)+R->GetDz());
	    tranS = new TGeoTranslation("",RinHole*TMath::Cos(t),
					RinHole*TMath::Sin(t),
					I->GetZ(4)+S->GetDz());
	    Kv->AddNode(Rv,i,tranR);
	    Iv->AddNode(Sv,i,tranS);
	} // end if
    } // end for i
    Int_t NcD=1,NcE=1,NcR=1,NcS=1;
    const Int_t Nbscrew=2,Nbpins=3,Nrailsc=4,Nrailp=2;
    Double_t da[] = {-3.5,-1.5,1.5,3.5};
    for(i=0;i<2;i++){ // Mounting for ITS-TPC bracket or ITS-Rails
	t0 = TMath::Pi()*((Double_t)i);
	for(j=-Nbscrew/2;j<=Nbscrew/2;j++)if(j!=0){//screws per ITS-TPC bracket
	    t = t0 + 5.0*((Double_t)j)*TMath::DegToRad();
	    tran = new TGeoTranslation("",RoutHole*TMath::Cos(t),
				      RoutHole*TMath::Sin(t),
				      B->GetZ(0)-D->GetDz());
	    Bv->AddNode(Dv,NcD,tran);
	    if(fDebug) cout << "D: NcD="<<NcD<<endl;
	    NcD++;
	} // end or j
	for(j=-Nbpins/2;j<=Nbpins/2;j++){ // pins per ITS-TPC bracket
	    t = t0 + 3.0*((Double_t)j)*TMath::DegToRad();
	    tran = new TGeoTranslation("",RoutHole*TMath::Cos(t),
				      RoutHole*TMath::Sin(t),
				      B->GetZ(0)-D->GetDz());
	    Bv->AddNode(Ev,NcE,tran);
	    if(fDebug) cout << "E: NcE="<<NcE<<endl;
	    NcE++;
	} // end or j
	t0 = (96.5+187.*((Double_t)i))*TMath::DegToRad();
	for(j=0;j<Nrailsc;j++){ // screws per ITS-rail bracket
	    t = t0+da[j]*TMath::DegToRad();
	    tran = new TGeoTranslation("",RoutHole*TMath::Cos(t),
				      RoutHole*TMath::Sin(t),
				      B->GetZ(0)-D->GetDz());
	    Bv->AddNode(Dv,NcD,tran);
	    if(fDebug) cout << "D2: NcD="<<NcD<<endl;
	    NcD++;
	} // end or j
	t0 = (91.5+184.*((Double_t)i))*TMath::DegToRad();
	for(j=-Nrailp/2;j<=Nrailp/2;j++)if(j!=0){ // pins per ITS-rail bracket
	    t = t0+(7.0*((Double_t)j))*TMath::DegToRad();
	    tran = new TGeoTranslation("",RoutHole*TMath::Cos(t),
				      RoutHole*TMath::Sin(t),
				      B->GetZ(0)-D->GetDz());
	    Bv->AddNode(Ev,NcE,tran);
	    if(fDebug) cout << "E2: NcE="<<NcE<<endl;
	    NcE++;
	} // end or j
    } // end for i
    for(i=0;i<Nmounts;i++){ // mounting points for SPD-cone+Beam-pipe support
	t0 = (45.0+((Double_t)i)*360./((Double_t)Nmounts))*TMath::DegToRad();
	for(j=-1;j<=1;j++)if(j!=0){ // 2 screws per bracket
	    t = t0+((Double_t)j)*0.5*DmountAngle;
	    tran = new TGeoTranslation("",RoutHole*TMath::Cos(t),
				      RoutHole*TMath::Sin(t),
				      B->GetZ(0)-D->GetDz());
	    Bv->AddNode(Dv,NcD,tran);
	    if(fDebug) cout << "D3: NcD="<<NcD<<endl;
	    NcD++;
	} // end for j
	for(j=0;j<1;j++){ // 1 pin per bracket
	    t = t0;
	    tran = new TGeoTranslation("",RoutHole*TMath::Cos(t),
				      RoutHole*TMath::Sin(t),
				      B->GetZ(0)-D->GetDz());
	    Bv->AddNode(Ev,NcE,tran);
	    if(fDebug) cout << "E3: NcE="<<NcE<<endl;
	    NcE++;
	} // end for j
    } // end for i
    tran = new TGeoTranslation("",TMath::Cos(T->GetPhi1()+0.5*T->GetDphi()),
			      TMath::Sin(T->GetPhi1()+0.5*T->GetDphi()),
			      T->GetZ(T->GetNz()-1)+R->GetDz());
    Tv->AddNode(Rv,NcR++,tran);
    tran = new TGeoTranslation("",TMath::Cos(U->GetPhi1()+0.5*U->GetDphi()),
			      TMath::Sin(U->GetPhi1()+0.5*U->GetDphi()),
			      U->GetZ(U->GetNz()-1)+S->GetDz());
    Uv->AddNode(Sv,NcS++,tran);
    for(i=0;i<NssdSupports;i++){ // mounting braclets for SSD/SDD 
	t0 = ((Double_t)i*360./((Double_t)NssdSupports));
	rot = new TGeoRotation("",0.0,0.0,t0);
	Kv->AddNode(Tv,i+1,rot);
	Iv->AddNode(Uv,i+1,rot);
	if(fDebug) cout << "T/U: copy number="<<i+1<<endl;
	//for(j=0;j<1;j++){ // 1 screws per bracket
	//    t = t0;
	//} // end for j
	for(j=0;j<2;j++)if(j!=0){ // 2 pin per bracket
	    t = t0 + ((Double_t)j)*0.5*DssdsddBracketAngle;
	    tran = new TGeoTranslation("",RinHole*TMath::Cos(t),
			      RinHole*TMath::Sin(t),
			      T->GetZ(T->GetNz()-1)-E->GetDz());
	    Kv->AddNode(Ev,NcE++,tran);
	} // end for j
    } // end for i
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
 
