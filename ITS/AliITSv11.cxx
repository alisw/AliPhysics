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
#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <TMath.h>
#include <TFile.h>    // only required for Tracking function?
#include <TObjArray.h>
#include <TLorentzVector.h>
#include <TObjString.h>
#include <TClonesArray.h>
#include <TBRIK.h>
#include <TSystem.h>


#include "AliRun.h"
#include "AliMagF.h"
#include "AliConst.h"
#include "AliITSGeant3Geometry.h"
#include "AliITShit.h"
#include "AliITSv11.h"
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
//
#include "AliITSGeometryITSV.h"
#include "AliITSGeometrySSDCone.h"
#include "AliITSGeometrySDDCone.h"

ClassImp(AliITSv11)

//______________________________________________________________________
AliITSv11::AliITSv11() : AliITS() {
    // Standard default constructor for the ITS version 11.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return
    //   A default constructed AliITSv11 class.

    fITSV = 0;
    fcS = 0;
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

    fITSV = 0;
    fcS = 0;
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

    if(fITSV!=0) delete fITSV;
    if(fcS!=0) delete fcS;
//    if(fcD!=0) delete fcD;
}
//______________________________________________________________________
AliITSv11::AliITSv11(const AliITSv11 &source){
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

    if(fITSV==0) fITSV = new AliITSGeometryITSV(this,"ALIC");
    if(fcS==0) fcS = new AliITSGeometrySSDCone(this,t,"TSV",1);

    fcS->BuildDisplayGeometry();
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

    if(fITSV==0) fITSV = new AliITSGeometryITSV(this,"ALIC");
    if(fcS==0) fcS = new AliITSGeometrySSDCone(this,t,"TSV",1);
    //
    fITSV->CreateG3Geometry();
    fcS->CreateG3Geometry("TSV",t);
    //
    fITSV->PositionGeometry("ALIC",1,t,0);
    fcS->PositionG3Geometry(fITSV->GetParams(),1,t,0);
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
/*
    TVector3 t(0.0,0.0,0.0);

    if(fITSV==0) fITSV = new AliITSGeometryITSV(this,"ALIC");
    if(fcS==0) fcS = new AliITSGeometrySSDCone(this,t,"TSV",1);

    fITSV->CreateG3Materials();
    fcS->CreateG3Materials();
*/


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

  // Freon
  Float_t afre[2]  = { 12.011,18.9984032 };
  Float_t zfre[2]  = { 6., 9. };
  Float_t wfre[2]  = { 5.,12. };
  Float_t densfre  = 1.5;

  // --- Define the various materials and media for GEANT --- 
  // AliMaterial(Int_t imat, const char* name, Float_t a, Float_t z,
  //              Float_t dens, Float_t radl, Float_t absl,
  //              Float_t *buf=0, Int_t nwbuf=0)
  //AliMedium(Int_t numed, const char *name, Int_t nmat,
  //          Int_t isvol, Int_t ifield, Float_t fieldm,
  //          Float_t tmaxfd, Float_t stemax, Float_t deemax,
  //          Float_t epsil, Float_t stmin, Float_t *ubuf=0, Int_t nbuf=0)
  AliMaterial(1,"SI$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
  AliMedium(1,"SI$",1,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

  AliMaterial(2,"SPD SI CHIP$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
  AliMedium(2,"SPD SI CHIP$",2,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

  AliMaterial(3,"SPD SI BUS$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
  AliMedium(3,"SPD SI BUS$",3,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

  AliMaterial(4,"C (M55J)$",0.12011E+02,0.60000E+01,0.1930E+01,0.22100E+02,0.99900E+03);
  AliMedium(4,"C (M55J)$",4,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(5,"AIR$",0.14610E+02,0.73000E+01,0.12050E-02,0.30423E+05,0.99900E+03);
  AliMedium(5,"AIR$",5,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

  AliMaterial(6,"GEN AIR$",0.14610E+02,0.73000E+01,0.12050E-02,0.30423E+05,0.99900E+03);
  AliMedium(6,"GEN AIR$",6,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

  AliMaterial(7,"SDD SI CHIP$",0.374952E+02,0.178184E+02,0.24485E+01,0.76931E+01,0.99900E+03);
  AliMedium(7,"SDD SI CHIP$",7,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

  AliMaterial(9,"SDD C (M55J)$",0.123565E+02,0.64561E+01,0.18097E+01,0.229570E+02,0.99900E+03);
  AliMedium(9,"SDD C (M55J)$",9,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(10,"SDD AIR$",0.14610E+02,0.73000E+01,0.12050E-02,0.30423E+05,0.99900E+03);
  AliMedium(10,"SDD AIR$",10,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

  AliMaterial(11,"AL$",0.26982E+02,0.13000E+02,0.26989E+01,0.89000E+01,0.99900E+03);
  AliMedium(11,"AL$",11,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(12,"WATER$",0.14322E+02,0.72167E+01,0.10000E+01,0.35759E+02,0.94951E+02);
  AliMedium(12,"WATER$",12,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMixture(13,"Freon$",afre,zfre,densfre,-2,wfre);
  AliMedium(13,"Freon$",13,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(14,"COPPER$",0.63546E+02,0.29000E+02,0.89600E+01,0.14300E+01,0.99900E+03);
  AliMedium(14,"COPPER$",14,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(15,"CERAMICS$",0.22314E+02,0.10856E+02,0.36000E+01,0.76200E+01,0.31901E+02);
  AliMedium(15,"CERAMICS$",15,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(20,"SSD C (M55J)$",0.12011E+02,0.60000E+01,0.1930E+01,0.22100E+02,0.99900E+03);
  AliMedium(20,"SSD C (M55J)$",20,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(21,"SSD AIR$",0.14610E+02,0.73000E+01,0.12050E-02,0.30423E+05,0.99900E+03);
  AliMedium(21,"SSD AIR$",21,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

  AliMaterial(25,"G10FR4$",0.17749E+02,0.88750E+01,0.18000E+01,0.21822E+02,0.99900E+03);
  AliMedium(25,"G10FR4$",25,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(26,"GEN C (M55J)$",0.12011E+02,0.60000E+01,0.1930E+01,0.22100E+02,0.99900E+03);
  AliMedium(26,"GEN C (M55J)$",26,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(27,"GEN Air$",0.14610E+02,0.73000E+01,0.12050E-02,0.30423E+05,0.99900E+03);
  AliMedium(27,"GEN Air$",27,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

  AliMaterial(51,"SPD SI$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
  AliMedium(51,"SPD SI$",51,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

  AliMaterial(52,"SPD SI CHIP$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
  AliMedium(52,"SPD SI CHIP$",52,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

  AliMaterial(53,"SPD SI BUS$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
  AliMedium(53,"SPD SI BUS$",53,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

  AliMaterial(54,"SPD C (M55J)$",0.12011E+02,0.60000E+01,0.1930E+01,0.22100E+02,0.99900E+03);
  AliMedium(54,"SPD C (M55J)$",54,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(55,"SPD AIR$",0.14610E+02,0.73000E+01,0.12050E-02,0.30423E+05,0.99900E+03);
  AliMedium(55,"SPD AIR$",55,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);

  AliMaterial(56,"SPD KAPTON(POLYCH2)$",0.14000E+02,0.71770E+01,0.13000E+01,0.31270E+02,0.99900E+03);
  AliMedium(56,"SPD KAPTON(POLYCH2)$",56,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(61,"EPOXY$",0.17749E+02,0.88750E+01,0.18000E+01,0.21822E+02,0.99900E+03);
  AliMedium(61,"EPOXY$",61,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(62,"SILICON$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
  AliMedium(62,"SILICON$",62,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);

  AliMaterial(63,"KAPTONH(POLYCH2)$",0.14000E+02,0.71770E+01,0.13000E+01,0.31270E+02,0.99900E+03);
  AliMedium(63,"KAPTONH(POLYCH2)$",63,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(64,"ALUMINUM$",0.26982E+02,0.13000E+02,0.26989E+01,0.89000E+01,0.99900E+03);
  AliMedium(64,"ALUMINUM$",64,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(65,"INOX$",0.55098E+02,0.2572E+02,0.7900E+01,0.17800E+01,0.99900E+03);
  AliMedium(65,"INOX$",65,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(68,"ROHACELL$",0.123974E+02,0.62363E+01,0.500E-01,0.80986E+03,0.99900E+03);
  AliMedium(68,"ROHACELL$",68,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(69,"SDD C AL (M55J)$",0.138802E+02,0.71315E+01,0.19837E+01,0.176542E+02,0.99900E+03);
  AliMedium(69,"SDD C AL (M55J)$",69,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(70,"SDDKAPTON (POLYCH2)$",0.14000E+02,0.71770E+01,0.13000E+01,0.31270E+02,0.99900E+03);
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

  AliMaterial(76,"SPDBUS(AL+KPT+EPOX)$",0.19509E+02,0.96502E+01,0.19060E+01,0.15413E+02,0.99900E+03);
  AliMedium(76,"SPDBUS(AL+KPT+EPOX)$",76,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);
               
  AliMaterial(77,"SDD X7R capacitors$",0.1157516E+03,0.477056E+02,0.67200E+01,0.14236E+01,0.99900E+03);
  AliMedium(77,"SDD X7R capacitors$",77,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(78,"SDD ruby sph. Al2O3$",0.218101E+02,0.106467E+02,0.39700E+01,0.48539E+01,0.99900E+03);
  AliMedium(78,"SDD ruby sph. Al2O3$",78,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(79,"SDD SI insensitive$",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
  AliMedium(79,"SDD SI insensitive$",79,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(80,"SDD HV microcable$",0.159379E+02,0.78598E+01,0.16087E+01,0.217906E+02,0.99900E+03);
  AliMedium(80,"SDD HV microcable$",80,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(81,"SDD LV+signal cable$",0.223689E+02,0.108531+02,0.21035E+01,0.13440E+02,0.99900E+03);
  AliMedium(81,"SDD LV+signal cable$",81,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(82,"SDD hybrid microcab$",0.218254E+02,0.106001E+02,0.20502E+01,0.137308E+02,0.99900E+03);
  AliMedium(82,"SDD hybrid microcab$",82,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(83,"SDD anode microcab$",0.186438E+02,0.91193E+01,0.17854E+01,0.176451E+02,0.99900E+03);
  AliMedium(83,"SDD anode microcab$",83,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(84,"SDD/SSD rings$",0.123565E+02,0.64561E+01,0.18097E+01,0.229570E+02,0.99900E+03);
  AliMedium(84,"SDD/SSD rings$",84,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);

  AliMaterial(85,"inox/alum$",0.321502E+02,0.153383E+02,0.30705E+01,0.69197E+01,0.99900E+03);
  AliMedium(85,"inox/alum$",85,0,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin);


  // Special media

  AliMaterial(90,"SPD shield$", 12.011, 6., 1.93/10. , 22.1*10., 999);
  AliMedium(90,"SPD shield$",90,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);

  AliMaterial(91, "SPD End ladder$", 47.0447, 21.7963, 3.6374, 4.4711, 999); 
  AliMedium(91,"SPD End ladder$",91,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);

  AliMaterial(92, "SPD cone$",28.0855, 14., 2.33, 9.36, 999);    
  AliMedium(92,"SPD cone$",92,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);

  AliMaterial(93, "SDD End ladder$", 69.9298, 29.8246, 0.3824, 36.5103, 999); 
  AliMedium(93,"SDD End ladder$",93,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);

  AliMaterial(94, "SDD cone$",63.546, 29., 1.15, 1.265, 999);
  AliMedium(94,"SDD cone$",94,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);

  AliMaterial(95, "SSD End ladder$", 32.0988, 15.4021, 0.68, 35.3238, 999); 
  AliMedium(95,"SSD End ladder$",95,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);
  
  AliMaterial(96, "SSD cone$",63.546, 29., 1.15, 1.265, 999);
  AliMedium(96,"SSD cone$",96,0,ifield,fieldm,tmaxfdServ,stemaxServ,deemaxServ,epsilServ,stminServ);
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
 
