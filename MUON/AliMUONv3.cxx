/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *      SigmaEffect_thetadegrees                                                                  *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpeateose. It is      *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 3    //
/////////////////////////////////////////////////////////

// Old MUONv1 class
// (AliMUONv1.h 1.11, AliMUONv1.cxx 1.60)
// - now replaced with a new one where geometry and materials
// are created using new geometry builders
// (See ALIMUON*GeometryBuilder classes)
// To be removed later

#include <TRandom.h>
#include <TF1.h>
#include <TClonesArray.h>
#include <TRandom.h> 
#include <TVirtualMC.h>

#include "AliMUONv3.h"
#include "AliConst.h" 
#include "AliMUONChamber.h"
#include "AliMUONConstants.h"
#include "AliMUONFactory.h"
#include "AliMUONHit.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONChamberGeometry.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliMC.h"

ClassImp(AliMUONv3)
 
//___________________________________________
AliMUONv3::AliMUONv3()
 : AliMUON(),
   fTrackMomentum(), fTrackPosition()
{
// Constructor
    fChambers   = 0;
    fStations   = 0;
    fStepManagerVersionOld  = kFALSE;
    fAngleEffect = kTRUE;
    fStepMaxInActiveGas     = 0.6;
    fStepSum    =  0x0;
    fDestepSum  =  0x0;
    fElossRatio =  0x0;
    fAngleEffect10   = 0x0;
    fAngleEffectNorma= 0x0;
} 
//___________________________________________
AliMUONv3::AliMUONv3(const char *name, const char *title)
  : AliMUON(name,title), 
    fTrackMomentum(), fTrackPosition()
{
// Constructor
    // By default include all stations
    fStations = new Int_t[5];
    for (Int_t i=0; i<5; i++) fStations[i] = 1;

    AliMUONFactory factory;
    factory.Build(this, title);

    fStepManagerVersionOld = kFALSE;
    fAngleEffect = kTRUE;
    fStepMaxInActiveGas = 0.6;

    fStepSum   = new Float_t [AliMUONConstants::NCh()];
    fDestepSum = new Float_t [AliMUONConstants::NCh()];
    for (Int_t i=0; i<AliMUONConstants::NCh(); i++) {
      fStepSum[i] =0.0;
      fDestepSum[i]=0.0;
    }
    // Ratio of particle mean eloss with respect MIP's Khalil Boudjemline, sep 2003, PhD.Thesis and Particle Data Book
    fElossRatio = new TF1("ElossRatio","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",0.5,5.); 
    fElossRatio->SetParameter(0,1.02138);
    fElossRatio->SetParameter(1,-9.54149e-02);
    fElossRatio->SetParameter(2,+7.83433e-02); 
    fElossRatio->SetParameter(3,-9.98208e-03);
    fElossRatio->SetParameter(4,+3.83279e-04);

    // Angle effect in tracking chambers at theta =10 degres as a function of ElossRatio (Khalil BOUDJEMLINE sep 2003 Ph.D Thesis) (in micrometers)
    fAngleEffect10 = new TF1("AngleEffect10","[0]+[1]*x+[2]*x*x",0.5,3.0);
    fAngleEffect10->SetParameter(0, 1.90691e+02);
    fAngleEffect10->SetParameter(1,-6.62258e+01);
    fAngleEffect10->SetParameter(2,+1.28247e+01);
    // Angle effect: Normalisation form theta=10 degres to theta between 0 and 10 (Khalil BOUDJEMLINE sep 2003 Ph.D Thesis)  
    // Angle with respect to the wires assuming that chambers are perpendicular to the z axis.
    fAngleEffectNorma = new TF1("AngleEffectNorma","[0]+[1]*x+[2]*x*x+[3]*x*x*x",0.0,10.0);
    fAngleEffectNorma->SetParameter(0,4.148);
    fAngleEffectNorma->SetParameter(1,-6.809e-01);
    fAngleEffectNorma->SetParameter(2,5.151e-02);
    fAngleEffectNorma->SetParameter(3,-1.490e-03);
}

//_____________________________________________________________________________
AliMUONv3::AliMUONv3(const AliMUONv3& right) 
  : AliMUON(right) 
{  
  // copy constructor (not implemented)

  Fatal("AliMUONv3", "Copy constructor not provided.");
}

//_____________________________________________________________________________
AliMUONv3& AliMUONv3::operator=(const AliMUONv3& right)
{
  // assignement operator (not implemented)

  // check assignement to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignement operator not provided.");
    
  return *this;  
}    

//___________________________________________
void AliMUONv3::CreateGeometry()
{
//
//   Note: all chambers have the same structure, which could be 
//   easily parameterised. This was intentionally not done in order
//   to give a starting point for the implementation of the actual 
//   design of each station. 
  Int_t *idtmed = fIdtmed->GetArray()-1099;

//   Distance between Stations
//
     Float_t bpar[3];
     Float_t tpar[3];
//      Float_t pgpar[10];
     Float_t zpos1, zpos2, zfpos;
     // Outer excess and inner recess for mother volume radius
     // with respect to ROuter and RInner
     Float_t dframep=.001; // Value for station 3 should be 6 ...
     // Width (RdPhi) of the frame crosses for stations 1 and 2 (cm)
//      Float_t dframep1=.001;
     Float_t dframep1 = 11.0;
//      Bool_t frameCrosses=kFALSE;     
     Bool_t frameCrosses=kTRUE;     
     Float_t *dum=0;
     
//      Float_t dframez=0.9;
     // Half of the total thickness of frame crosses (including DAlu)
     // for each chamber in stations 1 and 2:
     // 3% of X0 of composite material,
     // but taken as Aluminium here, with same thickness in number of X0
     Float_t dframez = 3. * 8.9 / 100;
//      Float_t dr;
     Float_t dstation;

//
//   Rotation matrices in the x-y plane  
     Int_t idrotm[1199];
//   phi=   0 deg
     AliMatrix(idrotm[1100],  90.,   0., 90.,  90., 0., 0.);
//   phi=  90 deg
     AliMatrix(idrotm[1101],  90.,  90., 90., 180., 0., 0.);
//   phi= 180 deg
     AliMatrix(idrotm[1102],  90., 180., 90., 270., 0., 0.);
//   phi= 270 deg
     AliMatrix(idrotm[1103],  90., 270., 90.,   0., 0., 0.);
//
     Float_t phi=2*TMath::Pi()/12/2;

//
//   pointer to the current chamber
//   pointer to the current chamber
     Int_t idAlu1=idtmed[1103]; // medium 4
     Int_t idAlu2=idtmed[1104]; // medium 5
//     Int_t idAlu1=idtmed[1100];
//     Int_t idAlu2=idtmed[1100];
     Int_t idAir=idtmed[1100]; // medium 1
//      Int_t idGas=idtmed[1105]; // medium 6 = Ar-isoC4H10 gas
     Int_t idGas=idtmed[1108]; // medium 9 = Ar-CO2 gas (80%+20%)
     

     AliMUONChamber *iChamber, *iChamber1, *iChamber2;

     if (fStations[0]) {
	 
//********************************************************************
//                            Station 1                             **
//********************************************************************
//  CONCENTRIC
     // indices 1 and 2 for first and second chambers in the station
     // iChamber (first chamber) kept for other quanties than Z,
     // assumed to be the same in both chambers
     iChamber1 = iChamber = (AliMUONChamber*) (*fChambers)[0];
     iChamber2 =(AliMUONChamber*) (*fChambers)[1];
     zpos1=iChamber1->Z(); 
     zpos2=iChamber2->Z();
     dstation = TMath::Abs(zpos2 - zpos1);
     // DGas decreased from standard one (0.5)
     iChamber->SetDGas(0.4); iChamber2->SetDGas(0.4);
     // DAlu increased from standard one (3% of X0),
     // because more electronics with smaller pads
     iChamber->SetDAlu(3.5 * 8.9 / 100.); iChamber2->SetDAlu(3.5 * 8.9 / 100.);
     zfpos=-(iChamber->DGas()+dframez+iChamber->DAlu())/2;
     
//
//   Mother volume
     tpar[0] = iChamber->RInner()-dframep; 
     tpar[1] = (iChamber->ROuter()+dframep)/TMath::Cos(phi);
     tpar[2] = dstation/5;

     gMC->Gsvolu("S01M", "TUBE", idAir, tpar, 3);
     gMC->Gsvolu("S02M", "TUBE", idAir, tpar, 3);
     gMC->Gspos("S01M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");
     gMC->Gspos("S02M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");     
// // Aluminium frames
// // Outer frames
//      pgpar[0] = 360/12/2;
//      pgpar[1] = 360.;
//      pgpar[2] = 12.;
//      pgpar[3] =   2;
//      pgpar[4] = -dframez/2;
//      pgpar[5] = iChamber->ROuter();
//      pgpar[6] = pgpar[5]+dframep1;
//      pgpar[7] = +dframez/2;
//      pgpar[8] = pgpar[5];
//      pgpar[9] = pgpar[6];
//      gMC->Gsvolu("S01O", "PGON", idAlu1, pgpar, 10);
//      gMC->Gsvolu("S02O", "PGON", idAlu1, pgpar, 10);
//      gMC->Gspos("S01O",1,"S01M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("S01O",2,"S01M", 0.,0.,+zfpos,  0,"ONLY");
//      gMC->Gspos("S02O",1,"S02M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("S02O",2,"S02M", 0.,0.,+zfpos,  0,"ONLY");
// //
// // Inner frame
//      tpar[0]= iChamber->RInner()-dframep1;
//      tpar[1]= iChamber->RInner();
//      tpar[2]= dframez/2;
//      gMC->Gsvolu("S01I", "TUBE", idAlu1, tpar, 3);
//      gMC->Gsvolu("S02I", "TUBE", idAlu1, tpar, 3);

//      gMC->Gspos("S01I",1,"S01M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("S01I",2,"S01M", 0.,0.,+zfpos,  0,"ONLY");
//      gMC->Gspos("S02I",1,"S02M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("S02I",2,"S02M", 0.,0.,+zfpos,  0,"ONLY");
//
// Frame Crosses
     if (frameCrosses) {
         // outside gas
         // security for inside mother volume
	 bpar[0] = (iChamber->ROuter() - iChamber->RInner())
	   * TMath::Cos(TMath::ASin(dframep1 /
				   (iChamber->ROuter() - iChamber->RInner())))
	   / 2.0;
	 bpar[1] = dframep1/2;
	 // total thickness will be (4 * bpar[2]) for each chamber,
	 // which has to be equal to (2 * dframez) - DAlu
	 bpar[2] = (2.0 * dframez - iChamber->DAlu()) / 4.0;
	 gMC->Gsvolu("S01B", "BOX", idAlu1, bpar, 3);
	 gMC->Gsvolu("S02B", "BOX", idAlu1, bpar, 3);
	 
	 gMC->Gspos("S01B",1,"S01M", -iChamber->RInner()-bpar[0] , 0, zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("S01B",2,"S01M",  iChamber->RInner()+bpar[0] , 0, zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("S01B",3,"S01M", 0, -iChamber->RInner()-bpar[0] , zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("S01B",4,"S01M", 0,  iChamber->RInner()+bpar[0] , zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("S01B",5,"S01M", -iChamber->RInner()-bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("S01B",6,"S01M", +iChamber->RInner()+bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("S01B",7,"S01M", 0, -iChamber->RInner()-bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("S01B",8,"S01M", 0, +iChamber->RInner()+bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
	 
	 gMC->Gspos("S02B",1,"S02M", -iChamber->RInner()-bpar[0] , 0, zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("S02B",2,"S02M",  iChamber->RInner()+bpar[0] , 0, zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("S02B",3,"S02M", 0, -iChamber->RInner()-bpar[0] , zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("S02B",4,"S02M", 0,  iChamber->RInner()+bpar[0] , zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("S02B",5,"S02M", -iChamber->RInner()-bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("S02B",6,"S02M", +iChamber->RInner()+bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("S02B",7,"S02M", 0, -iChamber->RInner()-bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("S02B",8,"S02M", 0, +iChamber->RInner()+bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
     }
//
//   Chamber Material represented by Alu sheet
     tpar[0]= iChamber->RInner();
     tpar[1]= iChamber->ROuter();
     tpar[2] = (iChamber->DGas()+iChamber->DAlu())/2;
     gMC->Gsvolu("S01A", "TUBE",  idAlu2, tpar, 3);
     gMC->Gsvolu("S02A", "TUBE",idAlu2, tpar, 3);
     gMC->Gspos("S01A", 1, "S01M", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("S02A", 1, "S02M", 0., 0., 0.,  0, "ONLY");
//     
//   Sensitive volumes
     // tpar[2] = iChamber->DGas();
     tpar[2] = iChamber->DGas()/2;
     gMC->Gsvolu("S01G", "TUBE", idGas, tpar, 3);
     gMC->Gsvolu("S02G", "TUBE", idGas, tpar, 3);
     gMC->Gspos("S01G", 1, "S01A", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("S02G", 1, "S02A", 0., 0., 0.,  0, "ONLY");
//
// Frame Crosses to be placed inside gas
     // NONE: chambers are sensitive everywhere
//      if (frameCrosses) {

// 	 dr = (iChamber->ROuter() - iChamber->RInner());
// 	 bpar[0] = TMath::Sqrt(dr*dr-dframep1*dframep1/4)/2;
// 	 bpar[1] = dframep1/2;
// 	 bpar[2] = iChamber->DGas()/2;
// 	 gMC->Gsvolu("S01F", "BOX", idAlu1, bpar, 3);
// 	 gMC->Gsvolu("S02F", "BOX", idAlu1, bpar, 3);
	 
// 	 gMC->Gspos("S01F",1,"S01G", +iChamber->RInner()+bpar[0] , 0, 0, 
// 		    idrotm[1100],"ONLY");
// 	 gMC->Gspos("S01F",2,"S01G", -iChamber->RInner()-bpar[0] , 0, 0, 
// 		    idrotm[1100],"ONLY");
// 	 gMC->Gspos("S01F",3,"S01G", 0, +iChamber->RInner()+bpar[0] , 0, 
// 		    idrotm[1101],"ONLY");
// 	 gMC->Gspos("S01F",4,"S01G", 0, -iChamber->RInner()-bpar[0] , 0, 
// 		    idrotm[1101],"ONLY");
	 
// 	 gMC->Gspos("S02F",1,"S02G", +iChamber->RInner()+bpar[0] , 0, 0, 
// 		    idrotm[1100],"ONLY");
// 	 gMC->Gspos("S02F",2,"S02G", -iChamber->RInner()-bpar[0] , 0, 0, 
// 		    idrotm[1100],"ONLY");
// 	 gMC->Gspos("S02F",3,"S02G", 0, +iChamber->RInner()+bpar[0] , 0, 
// 		    idrotm[1101],"ONLY");
// 	 gMC->Gspos("S02F",4,"S02G", 0, -iChamber->RInner()-bpar[0] , 0, 
// 		    idrotm[1101],"ONLY");
//      }
     }
     if (fStations[1]) {
	 
//********************************************************************
//                            Station 2                             **
//********************************************************************
     // indices 1 and 2 for first and second chambers in the station
     // iChamber (first chamber) kept for other quanties than Z,
     // assumed to be the same in both chambers
     iChamber1 = iChamber = (AliMUONChamber*) (*fChambers)[2];
     iChamber2 =(AliMUONChamber*) (*fChambers)[3];
     zpos1=iChamber1->Z(); 
     zpos2=iChamber2->Z();
     dstation = TMath::Abs(zpos2 - zpos1);
     // DGas and DAlu not changed from standard values
     zfpos=-(iChamber->DGas()+dframez+iChamber->DAlu())/2;
     
//
//   Mother volume
     tpar[0] = iChamber->RInner()-dframep; 
     tpar[1] = (iChamber->ROuter()+dframep)/TMath::Cos(phi);
     tpar[2] = dstation/5;

     gMC->Gsvolu("S03M", "TUBE", idAir, tpar, 3);
     gMC->Gsvolu("S04M", "TUBE", idAir, tpar, 3);
     gMC->Gspos("S03M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");
     gMC->Gspos("S04M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");
     gMC->Gsbool("S03M", "L3DO");
     gMC->Gsbool("S03M", "L3O1");
     gMC->Gsbool("S03M", "L3O2");
     gMC->Gsbool("S04M", "L3DO");
     gMC->Gsbool("S04M", "L3O1");
     gMC->Gsbool("S04M", "L3O2");

// // Aluminium frames
// // Outer frames
//      pgpar[0] = 360/12/2;
//      pgpar[1] = 360.;
//      pgpar[2] = 12.;
//      pgpar[3] =   2;
//      pgpar[4] = -dframez/2;
//      pgpar[5] = iChamber->ROuter();
//      pgpar[6] = pgpar[5]+dframep;
//      pgpar[7] = +dframez/2;
//      pgpar[8] = pgpar[5];
//      pgpar[9] = pgpar[6];
//      gMC->Gsvolu("S03O", "PGON", idAlu1, pgpar, 10);
//      gMC->Gsvolu("S04O", "PGON", idAlu1, pgpar, 10);
//      gMC->Gspos("S03O",1,"S03M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("S03O",2,"S03M", 0.,0.,+zfpos,  0,"ONLY");
//      gMC->Gspos("S04O",1,"S04M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("S04O",2,"S04M", 0.,0.,+zfpos,  0,"ONLY");
// //
// // Inner frame
//      tpar[0]= iChamber->RInner()-dframep;
//      tpar[1]= iChamber->RInner();
//      tpar[2]= dframez/2;
//      gMC->Gsvolu("S03I", "TUBE", idAlu1, tpar, 3);
//      gMC->Gsvolu("S04I", "TUBE", idAlu1, tpar, 3);

//      gMC->Gspos("S03I",1,"S03M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("S03I",2,"S03M", 0.,0.,+zfpos,  0,"ONLY");
//      gMC->Gspos("S04I",1,"S04M", 0.,0.,-zfpos,  0,"ONLY");
//      gMC->Gspos("S04I",2,"S04M", 0.,0.,+zfpos,  0,"ONLY");
//
// Frame Crosses
     if (frameCrosses) {
         // outside gas
         // security for inside mother volume
	 bpar[0] = (iChamber->ROuter() - iChamber->RInner())
	   * TMath::Cos(TMath::ASin(dframep1 /
				   (iChamber->ROuter() - iChamber->RInner())))
	   / 2.0;
	 bpar[1] = dframep1/2;
	 // total thickness will be (4 * bpar[2]) for each chamber,
	 // which has to be equal to (2 * dframez) - DAlu
	 bpar[2] = (2.0 * dframez - iChamber->DAlu()) / 4.0;
	 gMC->Gsvolu("S03B", "BOX", idAlu1, bpar, 3);
	 gMC->Gsvolu("S04B", "BOX", idAlu1, bpar, 3);
	 
	 gMC->Gspos("S03B",1,"S03M", -iChamber->RInner()-bpar[0] , 0, zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("S03B",2,"S03M", +iChamber->RInner()+bpar[0] , 0, zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("S03B",3,"S03M", 0, -iChamber->RInner()-bpar[0] , zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("S03B",4,"S03M", 0, +iChamber->RInner()+bpar[0] , zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("S03B",5,"S03M", -iChamber->RInner()-bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("S03B",6,"S03M", +iChamber->RInner()+bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("S03B",7,"S03M", 0, -iChamber->RInner()-bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("S03B",8,"S03M", 0, +iChamber->RInner()+bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
	 
	 gMC->Gspos("S04B",1,"S04M", -iChamber->RInner()-bpar[0] , 0, zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("S04B",2,"S04M", +iChamber->RInner()+bpar[0] , 0, zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("S04B",3,"S04M", 0, -iChamber->RInner()-bpar[0] , zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("S04B",4,"S04M", 0, +iChamber->RInner()+bpar[0] , zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("S04B",5,"S04M", -iChamber->RInner()-bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("S04B",6,"S04M", +iChamber->RInner()+bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("S04B",7,"S04M", 0, -iChamber->RInner()-bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("S04B",8,"S04M", 0, +iChamber->RInner()+bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
     }
//
//   Chamber Material represented by Alu sheet
     tpar[0]= iChamber->RInner();
     tpar[1]= iChamber->ROuter();
     tpar[2] = (iChamber->DGas()+iChamber->DAlu())/2;
     gMC->Gsvolu("S03A", "TUBE", idAlu2, tpar, 3);
     gMC->Gsvolu("S04A", "TUBE", idAlu2, tpar, 3);
     gMC->Gspos("S03A", 1, "S03M", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("S04A", 1, "S04M", 0., 0., 0.,  0, "ONLY");
//     
//   Sensitive volumes
     // tpar[2] = iChamber->DGas();
     tpar[2] = iChamber->DGas()/2;
     gMC->Gsvolu("S03G", "TUBE", idGas, tpar, 3);
     gMC->Gsvolu("S04G", "TUBE", idGas, tpar, 3);
     gMC->Gspos("S03G", 1, "S03A", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("S04G", 1, "S04A", 0., 0., 0.,  0, "ONLY");
//
// Frame Crosses to be placed inside gas 
     // NONE: chambers are sensitive everywhere
//      if (frameCrosses) {

// 	 dr = (iChamber->ROuter() - iChamber->RInner());
// 	 bpar[0] = TMath::Sqrt(dr*dr-dframep1*dframep1/4)/2;
// 	 bpar[1] = dframep1/2;
// 	 bpar[2] = iChamber->DGas()/2;
// 	 gMC->Gsvolu("S03F", "BOX", idAlu1, bpar, 3);
// 	 gMC->Gsvolu("S04F", "BOX", idAlu1, bpar, 3);
	 
// 	 gMC->Gspos("S03F",1,"S03G", +iChamber->RInner()+bpar[0] , 0, 0, 
// 		    idrotm[1100],"ONLY");
// 	 gMC->Gspos("S03F",2,"S03G", -iChamber->RInner()-bpar[0] , 0, 0, 
// 		    idrotm[1100],"ONLY");
// 	 gMC->Gspos("S03F",3,"S03G", 0, +iChamber->RInner()+bpar[0] , 0, 
// 		    idrotm[1101],"ONLY");
// 	 gMC->Gspos("S03F",4,"S03G", 0, -iChamber->RInner()-bpar[0] , 0, 
// 		    idrotm[1101],"ONLY");
	 
// 	 gMC->Gspos("S04F",1,"S04G", +iChamber->RInner()+bpar[0] , 0, 0, 
// 		    idrotm[1100],"ONLY");
// 	 gMC->Gspos("S04F",2,"S04G", -iChamber->RInner()-bpar[0] , 0, 0, 
// 		    idrotm[1100],"ONLY");
// 	 gMC->Gspos("S04F",3,"S04G", 0, +iChamber->RInner()+bpar[0] , 0, 
// 		    idrotm[1101],"ONLY");
// 	 gMC->Gspos("S04F",4,"S04G", 0, -iChamber->RInner()-bpar[0] , 0, 
// 		    idrotm[1101],"ONLY");
//      }
     }
     // define the id of tracking media:
     Int_t idCopper = idtmed[1110];
     Int_t idGlass  = idtmed[1111];
     Int_t idCarbon = idtmed[1112];
     Int_t idRoha   = idtmed[1113];

      // sensitive area: 40*40 cm**2
     const Float_t ksensLength = 40.; 
     const Float_t ksensHeight = 40.; 
     const Float_t ksensWidth  = 0.5; // according to TDR fig 2.120 
     const Int_t ksensMaterial = idGas;
     const Float_t kyOverlap   = 1.5; 

     // PCB dimensions in cm; width: 30 mum copper   
     const Float_t kpcbLength  = ksensLength; 
     const Float_t kpcbHeight  = 60.; 
     const Float_t kpcbWidth   = 0.003;   
     const Int_t   kpcbMaterial= idCopper;

     // Insulating material: 200 mum glass fiber glued to pcb  
     const Float_t kinsuLength = kpcbLength; 
     const Float_t kinsuHeight = kpcbHeight; 
     const Float_t kinsuWidth  = 0.020;   
     const Int_t kinsuMaterial = idGlass;

     // Carbon fiber panels: 200mum carbon/epoxy skin   
     const Float_t kpanelLength = ksensLength; 
     const Float_t kpanelHeight = ksensHeight; 
     const Float_t kpanelWidth  = 0.020;      
     const Int_t kpanelMaterial = idCarbon;

     // rohacell between the two carbon panels   
     const Float_t krohaLength = ksensLength; 
     const Float_t krohaHeight = ksensHeight; 
     const Float_t krohaWidth  = 0.5;
     const Int_t krohaMaterial = idRoha;

     // Frame around the slat: 2 sticks along length,2 along height  
     // H: the horizontal ones 
     const Float_t khFrameLength = kpcbLength; 
     const Float_t khFrameHeight = 1.5; 
     const Float_t khFrameWidth  = ksensWidth; 
     const Int_t khFrameMaterial = idGlass;

     // V: the vertical ones 
     const Float_t kvFrameLength = 4.0; 
     const Float_t kvFrameHeight = ksensHeight + khFrameHeight; 
     const Float_t kvFrameWidth  = ksensWidth;
     const Int_t kvFrameMaterial = idGlass;

     // B: the horizontal border filled with rohacell 
     const Float_t kbFrameLength = khFrameLength; 
     const Float_t kbFrameHeight = (kpcbHeight - ksensHeight)/2. - khFrameHeight; 
     const Float_t kbFrameWidth  = khFrameWidth;
     const Int_t kbFrameMaterial = idRoha;

     // NULOC: 30 mum copper + 200 mum vetronite (same radiation length as 14mum copper)
     const Float_t knulocLength = 2.5; 
     const Float_t knulocHeight = 7.5; 
     const Float_t knulocWidth  = 0.0030 + 0.0014; // equivalent copper width of vetronite; 
     const Int_t   knulocMaterial = idCopper;

     const Float_t kslatHeight = kpcbHeight; 
     const Float_t kslatWidth = ksensWidth + 2.*(kpcbWidth + kinsuWidth + 
					       2.* kpanelWidth + krohaWidth);
     const Int_t kslatMaterial = idAir;
     const Float_t kdSlatLength = kvFrameLength; // border on left and right 

     Float_t spar[3];  
     Int_t i, j;

     // the panel volume contains the rohacell

     Float_t twidth = 2 * kpanelWidth + krohaWidth; 
     Float_t panelpar[3] = { kpanelLength/2., kpanelHeight/2., twidth/2. }; 
     Float_t rohapar[3] = { krohaLength/2., krohaHeight/2., krohaWidth/2. }; 

     // insulating material contains PCB-> gas-> 2 borders filled with rohacell

     twidth = 2*(kinsuWidth + kpcbWidth) + ksensWidth;  
     Float_t insupar[3] = { kinsuLength/2., kinsuHeight/2., twidth/2. }; 
     twidth -= 2 * kinsuWidth; 
     Float_t pcbpar[3] = { kpcbLength/2., kpcbHeight/2., twidth/2. }; 
     Float_t senspar[3] = { ksensLength/2., ksensHeight/2., ksensWidth/2. }; 
     Float_t theight = 2*khFrameHeight + ksensHeight;
     Float_t hFramepar[3]={khFrameLength/2., theight/2., khFrameWidth/2.}; 
     Float_t bFramepar[3]={kbFrameLength/2., kbFrameHeight/2., kbFrameWidth/2.}; 
     Float_t vFramepar[3]={kvFrameLength/2., kvFrameHeight/2., kvFrameWidth/2.}; 
     Float_t nulocpar[3]={knulocLength/2., knulocHeight/2., knulocWidth/2.}; 
     Float_t xx;
     Float_t xxmax = (kbFrameLength - knulocLength)/2.; 
     Int_t index=0;
     
     if (fStations[2]) {
	 
//********************************************************************
//                            Station 3                             **
//********************************************************************
     // indices 1 and 2 for first and second chambers in the station
     // iChamber (first chamber) kept for other quanties than Z,
     // assumed to be the same in both chambers
     iChamber1 = iChamber = (AliMUONChamber*) (*fChambers)[4];
     iChamber2 =(AliMUONChamber*) (*fChambers)[5];
     zpos1=iChamber1->Z(); 
     zpos2=iChamber2->Z();
     dstation = TMath::Abs(zpos2 - zpos1);

//
//   Mother volume
     tpar[0] = iChamber->RInner()-dframep; 
     tpar[1] = (iChamber->ROuter()+dframep)/TMath::Cos(phi);
     tpar[2] = dstation/5;

     const char *slats5Mother = "S05M";
     const char *slats6Mother = "S06M";
     Float_t zoffs5 = 0;
     Float_t zoffs6 = 0;

     if (gAlice->GetModule("DIPO")) {
       slats5Mother="DDIP";
       slats6Mother="DDIP";

       zoffs5 = TMath::Abs(zpos1);
       zoffs6 = TMath::Abs(zpos2);
     }

     else {
       gMC->Gsvolu("S05M", "TUBE", idAir, tpar, 3);
       gMC->Gsvolu("S06M", "TUBE", idAir, tpar, 3);
       gMC->Gspos("S05M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");
       gMC->Gspos("S06M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");
     }

     // volumes for slat geometry (xx=5,..,10 chamber id): 
     // Sxx0 Sxx1 Sxx2 Sxx3  -->   Slat Mother volumes 
     // SxxG                          -->   Sensitive volume (gas)
     // SxxP                          -->   PCB (copper) 
     // SxxI                          -->   Insulator (vetronite) 
     // SxxC                          -->   Carbon panel 
     // SxxR                          -->   Rohacell
     // SxxH, SxxV                    -->   Horizontal and Vertical frames (vetronite)
     // SB5x                          -->   Volumes for the 35 cm long PCB
     // slat dimensions: slat is a MOTHER volume!!! made of air

     // only for chamber 5: slat 1 has a PCB shorter by 5cm!

     Float_t tlength = 35.;
     Float_t panelpar2[3]  = { tlength/2., panelpar[1],  panelpar[2]}; 
     Float_t rohapar2[3]   = { tlength/2., rohapar[1],   rohapar[2]}; 
     Float_t insupar2[3]   = { tlength/2., insupar[1],   insupar[2]}; 
     Float_t pcbpar2[3]    = { tlength/2., pcbpar[1],    pcbpar[2]}; 
     Float_t senspar2[3]   = { tlength/2., senspar[1],   senspar[2]}; 
     Float_t hFramepar2[3] = { tlength/2., hFramepar[1], hFramepar[2]}; 
     Float_t bFramepar2[3] = { tlength/2., bFramepar[1], bFramepar[2]}; 

     const Int_t knSlats3 = 5;  // number of slats per quadrant
     const Int_t knPCB3[knSlats3] = {3,3,4,3,2}; // n PCB per slat
     const Float_t kxpos3[knSlats3] = {31., 40., 0., 0., 0.};
     Float_t slatLength3[knSlats3]; 

     // create and position the slat (mother) volumes 

     char volNam5[5];
     char volNam6[5];
     Float_t xSlat3;

     Float_t spar2[3];
     for (i = 0; i<knSlats3; i++){
       slatLength3[i] = kpcbLength * knPCB3[i] + 2. * kdSlatLength; 
       xSlat3 = slatLength3[i]/2. - kvFrameLength/2. + kxpos3[i]; 
       if (i==1 || i==0) slatLength3[i] -=  2. *kdSlatLength; // frame out in PCB with circular border 
       Float_t ySlat31 =  ksensHeight * i - kyOverlap * i; 
       Float_t ySlat32 = -ksensHeight * i + kyOverlap * i; 
       spar[0] = slatLength3[i]/2.; 
       spar[1] = kslatHeight/2.;
       spar[2] = kslatWidth/2. * 1.01; 
       // take away 5 cm from the first slat in chamber 5
       Float_t xSlat32 = 0;
       if (i==1 || i==2) { // 1 pcb is shortened by 5cm
	 spar2[0] = spar[0]-5./2.;
	 xSlat32 = xSlat3 - 5/2.;
       }
       else {
	 spar2[0] = spar[0];
	 xSlat32 = xSlat3;
       }
       spar2[1] = spar[1];
       spar2[2] = spar[2]; 
       Float_t dzCh3=spar[2] * 1.01;
       // zSlat to be checked (odd downstream or upstream?)
       Float_t zSlat = (i%2 ==0)? -spar[2] : spar[2];

	if (gAlice->GetModule("DIPO")) {zSlat*=-1.;}

       sprintf(volNam5,"S05%d",i);
       gMC->Gsvolu(volNam5,"BOX",kslatMaterial,spar2,3);
       gMC->Gspos(volNam5, i*4+1,slats5Mother, -xSlat32, ySlat31, zoffs5-zSlat-2.*dzCh3, 0, "ONLY");
       gMC->Gspos(volNam5, i*4+2,slats5Mother, +xSlat32, ySlat31, zoffs5-zSlat+2.*dzCh3, 0, "ONLY");
       
	if (i>0) { 
	 gMC->Gspos(volNam5, i*4+3,slats5Mother,-xSlat32, ySlat32, zoffs5-zSlat-2.*dzCh3, 0, "ONLY");
	 gMC->Gspos(volNam5, i*4+4,slats5Mother,+xSlat32, ySlat32, zoffs5-zSlat+2.*dzCh3, 0, "ONLY");
       }
       sprintf(volNam6,"S06%d",i);
       gMC->Gsvolu(volNam6,"BOX",kslatMaterial,spar,3);
       gMC->Gspos(volNam6, i*4+1,slats6Mother,-xSlat3, ySlat31, zoffs6-zSlat-2.*dzCh3, 0, "ONLY");
       gMC->Gspos(volNam6, i*4+2,slats6Mother,+xSlat3, ySlat31, zoffs6-zSlat+2.*dzCh3, 0, "ONLY");
       if (i>0) { 
	 gMC->Gspos(volNam6, i*4+3,slats6Mother,-xSlat3, ySlat32, zoffs6-zSlat-2.*dzCh3, 0, "ONLY");
	 gMC->Gspos(volNam6, i*4+4,slats6Mother,+xSlat3, ySlat32, zoffs6-zSlat+2.*dzCh3, 0, "ONLY");
       }
     }

     // create the panel volume 
 
     gMC->Gsvolu("S05C","BOX",kpanelMaterial,panelpar,3);
     gMC->Gsvolu("SB5C","BOX",kpanelMaterial,panelpar2,3);
     gMC->Gsvolu("S06C","BOX",kpanelMaterial,panelpar,3);

     // create the rohacell volume 

     gMC->Gsvolu("S05R","BOX",krohaMaterial,rohapar,3);
     gMC->Gsvolu("SB5R","BOX",krohaMaterial,rohapar2,3);
     gMC->Gsvolu("S06R","BOX",krohaMaterial,rohapar,3);

     // create the insulating material volume 

     gMC->Gsvolu("S05I","BOX",kinsuMaterial,insupar,3);
     gMC->Gsvolu("SB5I","BOX",kinsuMaterial,insupar2,3);
     gMC->Gsvolu("S06I","BOX",kinsuMaterial,insupar,3);

     // create the PCB volume 

     gMC->Gsvolu("S05P","BOX",kpcbMaterial,pcbpar,3);
     gMC->Gsvolu("SB5P","BOX",kpcbMaterial,pcbpar2,3);
     gMC->Gsvolu("S06P","BOX",kpcbMaterial,pcbpar,3);
 
     // create the sensitive volumes,
     gMC->Gsvolu("S05G","BOX",ksensMaterial,dum,0);
     gMC->Gsvolu("S06G","BOX",ksensMaterial,dum,0);


     // create the vertical frame volume 

     gMC->Gsvolu("S05V","BOX",kvFrameMaterial,vFramepar,3);
     gMC->Gsvolu("S06V","BOX",kvFrameMaterial,vFramepar,3);

     // create the horizontal frame volume 

     gMC->Gsvolu("S05H","BOX",khFrameMaterial,hFramepar,3);
     gMC->Gsvolu("SB5H","BOX",khFrameMaterial,hFramepar2,3);
     gMC->Gsvolu("S06H","BOX",khFrameMaterial,hFramepar,3);

     // create the horizontal border volume 

     gMC->Gsvolu("S05B","BOX",kbFrameMaterial,bFramepar,3);
     gMC->Gsvolu("SB5B","BOX",kbFrameMaterial,bFramepar2,3);
     gMC->Gsvolu("S06B","BOX",kbFrameMaterial,bFramepar,3);

     index=0; 
     for (i = 0; i<knSlats3; i++){
       sprintf(volNam5,"S05%d",i);
       sprintf(volNam6,"S06%d",i);
       Float_t xvFrame  = (slatLength3[i] - kvFrameLength)/2.;
       Float_t xvFrame2  = xvFrame;
       if ( i==1 || i ==2 ) xvFrame2 -= 5./2.;
       // position the vertical frames 
       if (i!=1 && i!=0) { 
	 gMC->Gspos("S05V",2*i-1,volNam5, xvFrame2, 0., 0. , 0, "ONLY");
	 gMC->Gspos("S05V",2*i  ,volNam5,-xvFrame2, 0., 0. , 0, "ONLY");
	 gMC->Gspos("S06V",2*i-1,volNam6, xvFrame, 0., 0. , 0, "ONLY");
	 gMC->Gspos("S06V",2*i  ,volNam6,-xvFrame, 0., 0. , 0, "ONLY");
       }       
       // position the panels and the insulating material 
       for (j=0; j<knPCB3[i]; j++){
	 index++;
	 Float_t xx = ksensLength * (-knPCB3[i]/2.+j+.5); 
	 Float_t xx2 = xx + 5/2.; 
	 
	 Float_t zPanel = spar[2] - panelpar[2]; 
	 if ( (i==1 || i==2) && j == knPCB3[i]-1) { // 1 pcb is shortened by 5cm 
	   gMC->Gspos("SB5C",2*index-1,volNam5, xx, 0., zPanel , 0, "ONLY");
	   gMC->Gspos("SB5C",2*index  ,volNam5, xx, 0.,-zPanel , 0, "ONLY");
	   gMC->Gspos("SB5I",index    ,volNam5, xx, 0., 0      , 0, "ONLY");
	 }
	 else if ( (i==1 || i==2) && j < knPCB3[i]-1) {
	   gMC->Gspos("S05C",2*index-1,volNam5, xx2, 0., zPanel , 0, "ONLY");
	   gMC->Gspos("S05C",2*index  ,volNam5, xx2, 0.,-zPanel , 0, "ONLY");
	   gMC->Gspos("S05I",index    ,volNam5, xx2, 0., 0 , 0, "ONLY");
	 }
	 else {
	   gMC->Gspos("S05C",2*index-1,volNam5, xx, 0., zPanel , 0, "ONLY");
	   gMC->Gspos("S05C",2*index  ,volNam5, xx, 0.,-zPanel , 0, "ONLY");
	   gMC->Gspos("S05I",index    ,volNam5, xx, 0., 0 , 0, "ONLY");
	 }
	 gMC->Gspos("S06C",2*index-1,volNam6, xx, 0., zPanel , 0, "ONLY");
	 gMC->Gspos("S06C",2*index  ,volNam6, xx, 0.,-zPanel , 0, "ONLY");
	 gMC->Gspos("S06I",index,volNam6, xx, 0., 0 , 0, "ONLY");
       } 
     }
     
     // position the rohacell volume inside the panel volume
     gMC->Gspos("S05R",1,"S05C",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("SB5R",1,"SB5C",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("S06R",1,"S06C",0.,0.,0.,0,"ONLY"); 

     // position the PCB volume inside the insulating material volume
     gMC->Gspos("S05P",1,"S05I",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("SB5P",1,"SB5I",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("S06P",1,"S06I",0.,0.,0.,0,"ONLY"); 
     // position the horizontal frame volume inside the PCB volume
     gMC->Gspos("S05H",1,"S05P",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("SB5H",1,"SB5P",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("S06H",1,"S06P",0.,0.,0.,0,"ONLY"); 
     // position the sensitive volume inside the horizontal frame volume
     gMC->Gsposp("S05G",1,"S05H",0.,0.,0.,0,"ONLY",senspar,3); 
     gMC->Gsposp("S05G",1,"SB5H",0.,0.,0.,0,"ONLY",senspar2,3); 
     gMC->Gsposp("S06G",1,"S06H",0.,0.,0.,0,"ONLY",senspar,3); 
     // position the border volumes inside the PCB volume
     Float_t yborder = ( kpcbHeight - kbFrameHeight ) / 2.; 
     gMC->Gspos("S05B",1,"S05P",0., yborder,0.,0,"ONLY"); 
     gMC->Gspos("S05B",2,"S05P",0.,-yborder,0.,0,"ONLY"); 
     gMC->Gspos("SB5B",1,"SB5P",0., yborder,0.,0,"ONLY"); 
     gMC->Gspos("SB5B",2,"SB5P",0.,-yborder,0.,0,"ONLY"); 
     gMC->Gspos("S06B",1,"S06P",0., yborder,0.,0,"ONLY"); 
     gMC->Gspos("S06B",2,"S06P",0.,-yborder,0.,0,"ONLY"); 

     // create the NULOC volume and position it in the horizontal frame

     gMC->Gsvolu("S05N","BOX",knulocMaterial,nulocpar,3);
     gMC->Gsvolu("S06N","BOX",knulocMaterial,nulocpar,3);
     index = 0;
     Float_t xxmax2 = xxmax - 5./2.;
     for (xx = -xxmax; xx<=xxmax; xx+=2*knulocLength) { 
       index++; 
       gMC->Gspos("S05N",2*index-1,"S05B", xx, 0.,-kbFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S05N",2*index  ,"S05B", xx, 0., kbFrameWidth/4., 0, "ONLY");
       if (xx > -xxmax2 && xx< xxmax2) {
	 gMC->Gspos("S05N",2*index-1,"SB5B", xx, 0.,-kbFrameWidth/4., 0, "ONLY");
	 gMC->Gspos("S05N",2*index  ,"SB5B", xx, 0., kbFrameWidth/4., 0, "ONLY");
       }
       gMC->Gspos("S06N",2*index-1,"S06B", xx, 0.,-kbFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S06N",2*index  ,"S06B", xx, 0., kbFrameWidth/4., 0, "ONLY");
     }
     
     // position the volumes approximating the circular section of the pipe
     Float_t yoffs = ksensHeight/2. - kyOverlap; 
     Float_t epsilon = 0.001; 
     Int_t ndiv=6;
     Float_t divpar[3];
     Double_t dydiv= ksensHeight/ndiv;
     Double_t ydiv = yoffs -dydiv;
     Int_t imax=0; 
     imax = 1; 
     Float_t rmin = 33.; 
     Float_t z1 = spar[2], z2=2*spar[2]*1.01; 
     if (gAlice->GetModule("DIPO")) {z1*=-1.;}
     for (Int_t idiv=0;idiv<ndiv; idiv++){ 
       ydiv+= dydiv;
       Float_t xdiv = 0.; 
       if (ydiv<rmin) xdiv= rmin * TMath::Sin( TMath::ACos(ydiv/rmin) );
       divpar[0] = (kpcbLength-xdiv)/2.; 
       divpar[1] = dydiv/2. - epsilon;
       divpar[2] = ksensWidth/2.; 
       Float_t xvol=(kpcbLength+xdiv)/2.+1.999;
       Float_t yvol=ydiv + dydiv/2.; 
       //printf ("y ll = %f y ur = %f \n",yvol - divpar[1], yvol + divpar[1]); 
       gMC->Gsposp("S05G",imax+4*idiv+1,slats5Mother,-xvol, yvol, zoffs5-z1-z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S06G",imax+4*idiv+1,slats6Mother,-xvol, yvol, zoffs6-z1-z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S05G",imax+4*idiv+2,slats5Mother,-xvol,-yvol, zoffs5-z1-z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S06G",imax+4*idiv+2,slats6Mother,-xvol,-yvol, zoffs6-z1-z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S05G",imax+4*idiv+3,slats5Mother,+xvol, yvol, zoffs5-z1+z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S06G",imax+4*idiv+3,slats6Mother,+xvol, yvol, zoffs6-z1+z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S05G",imax+4*idiv+4,slats5Mother,+xvol,-yvol, zoffs5-z1+z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S06G",imax+4*idiv+4,slats6Mother,+xvol,-yvol, zoffs6-z1+z2, 0, "ONLY",divpar,3);
     }
     }
     
 if (fStations[3]) {

//********************************************************************
//                            Station 4                             **
//********************************************************************
     // indices 1 and 2 for first and second chambers in the station
     // iChamber (first chamber) kept for other quanties than Z,
     // assumed to be the same in both chambers
     iChamber1 = iChamber = (AliMUONChamber*) (*fChambers)[6];
     iChamber2 =(AliMUONChamber*) (*fChambers)[7];
     zpos1=iChamber1->Z(); 
     zpos2=iChamber2->Z();
     dstation = TMath::Abs(zpos2 - zpos1);
//      zfpos=-(iChamber->DGas()+dframez+iChamber->DAlu())/2; // not used any more
     
//
//   Mother volume
     tpar[0] = iChamber->RInner()-dframep; 
     tpar[1] = (iChamber->ROuter()+dframep)/TMath::Cos(phi);
     tpar[2] = dstation/4;

     gMC->Gsvolu("S07M", "TUBE", idAir, tpar, 3);
     gMC->Gsvolu("S08M", "TUBE", idAir, tpar, 3);
     gMC->Gspos("S07M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");
     gMC->Gspos("S08M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");
     

     const Int_t knSlats4 = 6;  // number of slats per quadrant
     const Int_t knPCB4[knSlats4] = {4,4,5,5,4,3}; // n PCB per slat
     const Float_t kxpos4[knSlats4] = {38.5, 40., 0., 0., 0., 0.};
     Float_t slatLength4[knSlats4];     

     // create and position the slat (mother) volumes 

     char volNam7[5];
     char volNam8[5];
     Float_t xSlat4;
     Float_t ySlat4;

     for (i = 0; i<knSlats4; i++){
       slatLength4[i] = kpcbLength * knPCB4[i] + 2. * kdSlatLength; 
       xSlat4 = slatLength4[i]/2. - kvFrameLength/2. + kxpos4[i]; 
       if (i==1) slatLength4[i] -=  2. *kdSlatLength; // frame out in PCB with circular border 
       ySlat4 =  ksensHeight * i - kyOverlap *i;
       
       spar[0] = slatLength4[i]/2.; 
       spar[1] = kslatHeight/2.;
       spar[2] = kslatWidth/2.*1.01; 
       Float_t dzCh4=spar[2]*1.01;
       // zSlat to be checked (odd downstream or upstream?)
       Float_t zSlat = (i%2 ==0)? spar[2] : -spar[2]; 
       sprintf(volNam7,"S07%d",i);
       gMC->Gsvolu(volNam7,"BOX",kslatMaterial,spar,3);
       gMC->Gspos(volNam7, i*4+1,"S07M",-xSlat4, ySlat4, -zSlat-2.*dzCh4, 0, "ONLY");
       gMC->Gspos(volNam7, i*4+2,"S07M",+xSlat4, ySlat4, -zSlat+2.*dzCh4, 0, "ONLY");
       if (i>0) { 
	 gMC->Gspos(volNam7, i*4+3,"S07M",-xSlat4,-ySlat4, -zSlat-2.*dzCh4, 0, "ONLY");
	 gMC->Gspos(volNam7, i*4+4,"S07M",+xSlat4,-ySlat4, -zSlat+2.*dzCh4, 0, "ONLY");
       }
       sprintf(volNam8,"S08%d",i);
       gMC->Gsvolu(volNam8,"BOX",kslatMaterial,spar,3);
       gMC->Gspos(volNam8, i*4+1,"S08M",-xSlat4, ySlat4, -zSlat-2.*dzCh4, 0, "ONLY");
       gMC->Gspos(volNam8, i*4+2,"S08M",+xSlat4, ySlat4, -zSlat+2.*dzCh4, 0, "ONLY");
       if (i>0) { 
	 gMC->Gspos(volNam8, i*4+3,"S08M",-xSlat4,-ySlat4, -zSlat-2.*dzCh4, 0, "ONLY");
	 gMC->Gspos(volNam8, i*4+4,"S08M",+xSlat4,-ySlat4, -zSlat+2.*dzCh4, 0, "ONLY");
       }
     }
     

     // create the panel volume 
 
     gMC->Gsvolu("S07C","BOX",kpanelMaterial,panelpar,3);
     gMC->Gsvolu("S08C","BOX",kpanelMaterial,panelpar,3);

     // create the rohacell volume 

     gMC->Gsvolu("S07R","BOX",krohaMaterial,rohapar,3);
     gMC->Gsvolu("S08R","BOX",krohaMaterial,rohapar,3);

     // create the insulating material volume 

     gMC->Gsvolu("S07I","BOX",kinsuMaterial,insupar,3);
     gMC->Gsvolu("S08I","BOX",kinsuMaterial,insupar,3);

     // create the PCB volume 

     gMC->Gsvolu("S07P","BOX",kpcbMaterial,pcbpar,3);
     gMC->Gsvolu("S08P","BOX",kpcbMaterial,pcbpar,3);
 
     // create the sensitive volumes,

     gMC->Gsvolu("S07G","BOX",ksensMaterial,dum,0);
     gMC->Gsvolu("S08G","BOX",ksensMaterial,dum,0);

     // create the vertical frame volume 

     gMC->Gsvolu("S07V","BOX",kvFrameMaterial,vFramepar,3);
     gMC->Gsvolu("S08V","BOX",kvFrameMaterial,vFramepar,3);

     // create the horizontal frame volume 

     gMC->Gsvolu("S07H","BOX",khFrameMaterial,hFramepar,3);
     gMC->Gsvolu("S08H","BOX",khFrameMaterial,hFramepar,3);

     // create the horizontal border volume 

     gMC->Gsvolu("S07B","BOX",kbFrameMaterial,bFramepar,3);
     gMC->Gsvolu("S08B","BOX",kbFrameMaterial,bFramepar,3);

     index=0; 
     for (i = 0; i<knSlats4; i++){
       sprintf(volNam7,"S07%d",i);
       sprintf(volNam8,"S08%d",i);
       Float_t xvFrame  = (slatLength4[i] - kvFrameLength)/2.;
       // position the vertical frames 
       if (i!=1 && i!=0) { 
	 gMC->Gspos("S07V",2*i-1,volNam7, xvFrame, 0., 0. , 0, "ONLY");
	 gMC->Gspos("S07V",2*i  ,volNam7,-xvFrame, 0., 0. , 0, "ONLY");
	 gMC->Gspos("S08V",2*i-1,volNam8, xvFrame, 0., 0. , 0, "ONLY");
	 gMC->Gspos("S08V",2*i  ,volNam8,-xvFrame, 0., 0. , 0, "ONLY");
       }
       // position the panels and the insulating material 
       for (j=0; j<knPCB4[i]; j++){
	 index++;
	 Float_t xx = ksensLength * (-knPCB4[i]/2.+j+.5); 

	 Float_t zPanel = spar[2] - panelpar[2]; 
	 gMC->Gspos("S07C",2*index-1,volNam7, xx, 0., zPanel , 0, "ONLY");
	 gMC->Gspos("S07C",2*index  ,volNam7, xx, 0.,-zPanel , 0, "ONLY");
	 gMC->Gspos("S08C",2*index-1,volNam8, xx, 0., zPanel , 0, "ONLY");
	 gMC->Gspos("S08C",2*index  ,volNam8, xx, 0.,-zPanel , 0, "ONLY");

	 gMC->Gspos("S07I",index,volNam7, xx, 0., 0 , 0, "ONLY");
	 gMC->Gspos("S08I",index,volNam8, xx, 0., 0 , 0, "ONLY");
       } 
     }

     // position the rohacell volume inside the panel volume
     gMC->Gspos("S07R",1,"S07C",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("S08R",1,"S08C",0.,0.,0.,0,"ONLY"); 

     // position the PCB volume inside the insulating material volume
     gMC->Gspos("S07P",1,"S07I",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("S08P",1,"S08I",0.,0.,0.,0,"ONLY"); 
     // position the horizontal frame volume inside the PCB volume
     gMC->Gspos("S07H",1,"S07P",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("S08H",1,"S08P",0.,0.,0.,0,"ONLY"); 
     // position the sensitive volume inside the horizontal frame volume
     gMC->Gsposp("S07G",1,"S07H",0.,0.,0.,0,"ONLY",senspar,3); 
     gMC->Gsposp("S08G",1,"S08H",0.,0.,0.,0,"ONLY",senspar,3); 
     // position the border volumes inside the PCB volume
     Float_t yborder = ( kpcbHeight - kbFrameHeight ) / 2.; 
     gMC->Gspos("S07B",1,"S07P",0., yborder,0.,0,"ONLY"); 
     gMC->Gspos("S07B",2,"S07P",0.,-yborder,0.,0,"ONLY"); 
     gMC->Gspos("S08B",1,"S08P",0., yborder,0.,0,"ONLY"); 
     gMC->Gspos("S08B",2,"S08P",0.,-yborder,0.,0,"ONLY"); 

     // create the NULOC volume and position it in the horizontal frame

     gMC->Gsvolu("S07N","BOX",knulocMaterial,nulocpar,3);
     gMC->Gsvolu("S08N","BOX",knulocMaterial,nulocpar,3);
     index = 0;
     for (xx = -xxmax; xx<=xxmax; xx+=2*knulocLength) { 
       index++; 
       gMC->Gspos("S07N",2*index-1,"S07B", xx, 0.,-kbFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S07N",2*index  ,"S07B", xx, 0., kbFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S08N",2*index-1,"S08B", xx, 0.,-kbFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S08N",2*index  ,"S08B", xx, 0., kbFrameWidth/4., 0, "ONLY");
     }

     // position the volumes approximating the circular section of the pipe
     Float_t yoffs = ksensHeight/2. - kyOverlap; 
     Float_t epsilon = 0.001; 
     Int_t ndiv=6;
     Float_t divpar[3];
     Double_t dydiv= ksensHeight/ndiv;
     Double_t ydiv = yoffs -dydiv;
     Int_t imax=0; 
     imax = 1; 
     Float_t rmin = 40.; 
     Float_t z1 = -spar[2], z2=2*spar[2]*1.01; 
     for (Int_t idiv=0;idiv<ndiv; idiv++){ 
       ydiv+= dydiv;
       Float_t xdiv = 0.; 
       if (ydiv<rmin) xdiv= rmin * TMath::Sin( TMath::ACos(ydiv/rmin) );
       divpar[0] = (kpcbLength-xdiv)/2.; 
       divpar[1] = dydiv/2. - epsilon;
       divpar[2] = ksensWidth/2.; 
       Float_t xvol=(kpcbLength+xdiv)/2.+1.999;
       Float_t yvol=ydiv + dydiv/2.;
       gMC->Gsposp("S07G",imax+4*idiv+1,"S07M", -xvol, yvol, -z1-z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S08G",imax+4*idiv+1,"S08M", -xvol, yvol, -z1-z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S07G",imax+4*idiv+2,"S07M", -xvol,-yvol, -z1-z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S08G",imax+4*idiv+2,"S08M", -xvol,-yvol, -z1-z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S07G",imax+4*idiv+3,"S07M", xvol, yvol, -z1+z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S08G",imax+4*idiv+3,"S08M", xvol, yvol, -z1+z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S07G",imax+4*idiv+4,"S07M", xvol,-yvol, -z1+z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S08G",imax+4*idiv+4,"S08M", xvol,-yvol, -z1+z2, 0, "ONLY",divpar,3);
     }





 }

 if (fStations[4]) {
     

//********************************************************************
//                            Station 5                             **
//********************************************************************
     // indices 1 and 2 for first and second chambers in the station
     // iChamber (first chamber) kept for other quanties than Z,
     // assumed to be the same in both chambers
     iChamber1 = iChamber = (AliMUONChamber*) (*fChambers)[8];
     iChamber2 =(AliMUONChamber*) (*fChambers)[9];
     zpos1=iChamber1->Z(); 
     zpos2=iChamber2->Z();
     dstation = TMath::Abs(zpos2 - zpos1);
//      zfpos=-(iChamber->DGas()+dframez+iChamber->DAlu())/2; // not used any more
     
//
//   Mother volume
     tpar[0] = iChamber->RInner()-dframep; 
     tpar[1] = (iChamber->ROuter()+dframep)/TMath::Cos(phi);
     tpar[2] = dstation/5.;

     gMC->Gsvolu("S09M", "TUBE", idAir, tpar, 3);
     gMC->Gsvolu("S10M", "TUBE", idAir, tpar, 3);
     gMC->Gspos("S09M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");
     gMC->Gspos("S10M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");


     const Int_t knSlats5 = 7;  // number of slats per quadrant
     const Int_t knPCB5[knSlats5] = {5,5,6,6,5,4,3}; // n PCB per slat
     const Float_t kxpos5[knSlats5] = {38.5, 40., 0., 0., 0., 0., 0.};
     Float_t slatLength5[knSlats5]; 
     char volNam9[5];
     char volNam10[5];
     Float_t xSlat5;
     Float_t ySlat5;

     for (i = 0; i<knSlats5; i++){
       slatLength5[i] = kpcbLength * knPCB5[i] + 2. * kdSlatLength; 
       xSlat5 = slatLength5[i]/2. - kvFrameLength/2. +kxpos5[i]; 
       if (i==1 || i==0) slatLength5[i] -=  2. *kdSlatLength; // frame out in PCB with circular border 
       ySlat5 = ksensHeight * i - kyOverlap * i; 
       spar[0] = slatLength5[i]/2.; 
       spar[1] = kslatHeight/2.;
       spar[2] = kslatWidth/2. * 1.01; 
       Float_t dzCh5=spar[2]*1.01;
       // zSlat to be checked (odd downstream or upstream?)
       Float_t zSlat = (i%2 ==0)? -spar[2] : spar[2]; 
       sprintf(volNam9,"S09%d",i);
       gMC->Gsvolu(volNam9,"BOX",kslatMaterial,spar,3);
       gMC->Gspos(volNam9, i*4+1,"S09M",-xSlat5, ySlat5, -zSlat-2.*dzCh5, 0, "ONLY");
       gMC->Gspos(volNam9, i*4+2,"S09M",+xSlat5, ySlat5, -zSlat+2.*dzCh5, 0, "ONLY");
       if (i>0) { 
	   gMC->Gspos(volNam9, i*4+3,"S09M",-xSlat5,-ySlat5, -zSlat-2.*dzCh5, 0, "ONLY");
	   gMC->Gspos(volNam9, i*4+4,"S09M",+xSlat5,-ySlat5, -zSlat+2.*dzCh5, 0, "ONLY");
       }
       sprintf(volNam10,"S10%d",i);
       gMC->Gsvolu(volNam10,"BOX",kslatMaterial,spar,3);
       gMC->Gspos(volNam10, i*4+1,"S10M",-xSlat5, ySlat5, -zSlat-2.*dzCh5, 0, "ONLY");
       gMC->Gspos(volNam10, i*4+2,"S10M",+xSlat5, ySlat5, -zSlat+2.*dzCh5, 0, "ONLY");
       if (i>0) { 
	   gMC->Gspos(volNam10, i*4+3,"S10M",-xSlat5,-ySlat5, -zSlat-2.*dzCh5, 0, "ONLY");
	   gMC->Gspos(volNam10, i*4+4,"S10M",+xSlat5,-ySlat5, -zSlat+2.*dzCh5, 0, "ONLY");
       }
     }

     // create the panel volume 
 
     gMC->Gsvolu("S09C","BOX",kpanelMaterial,panelpar,3);
     gMC->Gsvolu("S10C","BOX",kpanelMaterial,panelpar,3);

     // create the rohacell volume 

     gMC->Gsvolu("S09R","BOX",krohaMaterial,rohapar,3);
     gMC->Gsvolu("S10R","BOX",krohaMaterial,rohapar,3);

     // create the insulating material volume 

     gMC->Gsvolu("S09I","BOX",kinsuMaterial,insupar,3);
     gMC->Gsvolu("S10I","BOX",kinsuMaterial,insupar,3);

     // create the PCB volume 

     gMC->Gsvolu("S09P","BOX",kpcbMaterial,pcbpar,3);
     gMC->Gsvolu("S10P","BOX",kpcbMaterial,pcbpar,3);
 
     // create the sensitive volumes,

     gMC->Gsvolu("S09G","BOX",ksensMaterial,dum,0);
     gMC->Gsvolu("S10G","BOX",ksensMaterial,dum,0);

     // create the vertical frame volume 

     gMC->Gsvolu("S09V","BOX",kvFrameMaterial,vFramepar,3);
     gMC->Gsvolu("S10V","BOX",kvFrameMaterial,vFramepar,3);

     // create the horizontal frame volume 

     gMC->Gsvolu("S09H","BOX",khFrameMaterial,hFramepar,3);
     gMC->Gsvolu("S10H","BOX",khFrameMaterial,hFramepar,3);

     // create the horizontal border volume 

     gMC->Gsvolu("S09B","BOX",kbFrameMaterial,bFramepar,3);
     gMC->Gsvolu("S10B","BOX",kbFrameMaterial,bFramepar,3);

     index=0; 
     for (i = 0; i<knSlats5; i++){
       sprintf(volNam9,"S09%d",i);
       sprintf(volNam10,"S10%d",i);
       Float_t xvFrame  = (slatLength5[i] - kvFrameLength)/2.;
       // position the vertical frames 
       if (i!=1 && i!=0) { 
	 gMC->Gspos("S09V",2*i-1,volNam9, xvFrame, 0., 0. , 0, "ONLY");
	 gMC->Gspos("S09V",2*i  ,volNam9,-xvFrame, 0., 0. , 0, "ONLY");
	 gMC->Gspos("S10V",2*i-1,volNam10, xvFrame, 0., 0. , 0, "ONLY");
	 gMC->Gspos("S10V",2*i  ,volNam10,-xvFrame, 0., 0. , 0, "ONLY");
       }
       
       // position the panels and the insulating material 
       for (j=0; j<knPCB5[i]; j++){
	 index++;
	 Float_t xx = ksensLength * (-knPCB5[i]/2.+j+.5); 

	 Float_t zPanel = spar[2] - panelpar[2]; 
	 gMC->Gspos("S09C",2*index-1,volNam9, xx, 0., zPanel , 0, "ONLY");
	 gMC->Gspos("S09C",2*index  ,volNam9, xx, 0.,-zPanel , 0, "ONLY");
	 gMC->Gspos("S10C",2*index-1,volNam10, xx, 0., zPanel , 0, "ONLY");
	 gMC->Gspos("S10C",2*index  ,volNam10, xx, 0.,-zPanel , 0, "ONLY");

	 gMC->Gspos("S09I",index,volNam9, xx, 0., 0 , 0, "ONLY");
	 gMC->Gspos("S10I",index,volNam10, xx, 0., 0 , 0, "ONLY");
       } 
     }

     // position the rohacell volume inside the panel volume
     gMC->Gspos("S09R",1,"S09C",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("S10R",1,"S10C",0.,0.,0.,0,"ONLY"); 

     // position the PCB volume inside the insulating material volume
     gMC->Gspos("S09P",1,"S09I",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("S10P",1,"S10I",0.,0.,0.,0,"ONLY"); 
     // position the horizontal frame volume inside the PCB volume
     gMC->Gspos("S09H",1,"S09P",0.,0.,0.,0,"ONLY"); 
     gMC->Gspos("S10H",1,"S10P",0.,0.,0.,0,"ONLY"); 
     // position the sensitive volume inside the horizontal frame volume
     gMC->Gsposp("S09G",1,"S09H",0.,0.,0.,0,"ONLY",senspar,3); 
     gMC->Gsposp("S10G",1,"S10H",0.,0.,0.,0,"ONLY",senspar,3); 
     // position the border volumes inside the PCB volume
     Float_t yborder = ( kpcbHeight - kbFrameHeight ) / 2.; 
     gMC->Gspos("S09B",1,"S09P",0., yborder,0.,0,"ONLY"); 
     gMC->Gspos("S09B",2,"S09P",0.,-yborder,0.,0,"ONLY"); 
     gMC->Gspos("S10B",1,"S10P",0., yborder,0.,0,"ONLY"); 
     gMC->Gspos("S10B",2,"S10P",0.,-yborder,0.,0,"ONLY"); 

     // create the NULOC volume and position it in the horizontal frame

     gMC->Gsvolu("S09N","BOX",knulocMaterial,nulocpar,3);
     gMC->Gsvolu("S10N","BOX",knulocMaterial,nulocpar,3);
     index = 0;
     for (xx = -xxmax; xx<=xxmax; xx+=2*knulocLength) { 
       index++; 
       gMC->Gspos("S09N",2*index-1,"S09B", xx, 0.,-kbFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S09N",2*index  ,"S09B", xx, 0., kbFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S10N",2*index-1,"S10B", xx, 0.,-kbFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S10N",2*index  ,"S10B", xx, 0., kbFrameWidth/4., 0, "ONLY");
     }
     // position the volumes approximating the circular section of the pipe
     Float_t yoffs = ksensHeight/2. - kyOverlap; 
     Float_t epsilon = 0.001; 
     Int_t ndiv=6;
     Float_t divpar[3];
     Double_t dydiv= ksensHeight/ndiv;
     Double_t ydiv = yoffs -dydiv;
     Int_t imax=0; 
     //     for (Int_t islat=0; islat<knSlats3; islat++) imax += knPCB3[islat]; 
     imax = 1; 
     Float_t rmin = 40.; 
     Float_t z1 = spar[2], z2=2*spar[2]*1.01; 
     for (Int_t idiv=0;idiv<ndiv; idiv++){ 
       ydiv+= dydiv;
       Float_t xdiv = 0.; 
       if (ydiv<rmin) xdiv= rmin * TMath::Sin( TMath::ACos(ydiv/rmin) );
       divpar[0] = (kpcbLength-xdiv)/2.; 
       divpar[1] = dydiv/2. - epsilon;
       divpar[2] = ksensWidth/2.; 
       Float_t xvol=(kpcbLength+xdiv)/2. + 1.999;
       Float_t yvol=ydiv + dydiv/2.;
       gMC->Gsposp("S09G",imax+4*idiv+1,"S09M", -xvol, yvol, -z1-z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S10G",imax+4*idiv+1,"S10M", -xvol, yvol, -z1-z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S09G",imax+4*idiv+2,"S09M", -xvol,-yvol, -z1-z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S10G",imax+4*idiv+2,"S10M", -xvol,-yvol, -z1-z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S09G",imax+4*idiv+3,"S09M", +xvol, yvol, -z1+z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S10G",imax+4*idiv+3,"S10M", +xvol, yvol, -z1+z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S09G",imax+4*idiv+4,"S09M", +xvol,-yvol, -z1+z2, 0, "ONLY",divpar,3);
       gMC->Gsposp("S10G",imax+4*idiv+4,"S10M", +xvol,-yvol, -z1+z2, 0, "ONLY",divpar,3);
     }

 }
 
//********************************************************************
//                            Trigger                               **
//******************************************************************** 
 /* 
    zpos1 and zpos2 are the middle of the first and second
    planes of station 1 (+1m for second station):
    zpos1=(zpos1m+zpos1p)/2=(15999+16071)/2=16035 mm, thick/2=40 mm
    zpos2=(zpos2m+zpos2p)/2=(16169+16241)/2=16205 mm, thick/2=40 mm
    zposxm and zposxp= middles of gaz gaps within a detection plane
    rem: the total thickness accounts for 1 mm of al on both
    side of the RPCs (see zpos1 and zpos2)
 */

// vertical gap between right and left chambers (kDXZERO*2=4cm)
 const Float_t kDXZERO=2.; 
// main distances for chamber definition in first plane/first station
 const Float_t kXMIN=34.;       
 const Float_t kXMED=51.;                                
 const Float_t kXMAX=272.; 
// kXMAX will become 255. in real life. segmentation to be updated accordingly
// (see fig.2-4 & 2-5 of Local Trigger Board PRR)
 const Float_t kYMIN=34.;                              
 const Float_t kYMAX=51.;                              
// inner/outer radius of flange between beam shield. and chambers (1/station)
 const Float_t kRMIN[2]={50.,50.};
 const Float_t kRMAX[2]={64.,68.};
// z position of the middle of the gas gap in mother vol 
 const Float_t kZm=-3.6;
 const Float_t kZp=+3.6;     
 
 iChamber1 = (AliMUONChamber*) (*fChambers)[10];     
 zpos1 = iChamber1->Z();

// ratio of zpos1m/zpos1p and inverse for first plane
 Float_t zmp=(zpos1+3.6)/(zpos1-3.6);
 Float_t zpm=1./zmp;
 
 Int_t icount=0; // chamber counter (0 1 2 3)
 
 for (Int_t istation=0; istation<2; istation++) { // loop on stations	 
     for (Int_t iplane=0; iplane<2; iplane++) {	  // loop on detection planes
	 
	 Int_t iVolNum=1; // counter Volume Number
	 icount = Int_t(iplane*TMath::Power(2,0))+
	     Int_t(istation*TMath::Power(2,1));
	 
	 char volPlane[5]; 
	 sprintf(volPlane,"SM%d%d",istation+1,iplane+1);
	 
	 iChamber = (AliMUONChamber*) (*fChambers)[10+icount];
	 Float_t zpos = iChamber->Z();	     
	 
// mother volume 
	 tpar[0] = iChamber->RInner(); 
	 tpar[1] = iChamber->ROuter(); 
	 tpar[2] = 4.0;    
	 gMC->Gsvolu(volPlane,"TUBE",idAir,tpar,3);
	 
// Flange between beam shielding and RPC 
	 tpar[0]= kRMIN[istation];
	 tpar[1]= kRMAX[istation];
	 tpar[2]= 4.0;
	 
	 char volFlange[5];
	 sprintf(volFlange,"SF%dA",icount+1);	 
	 gMC->Gsvolu(volFlange,"TUBE",idAlu1,tpar,3);     //Al
	 gMC->Gspos(volFlange,1,volPlane,0.,0.,0.,0,"MANY");
	 
// scaling factor
	 Float_t zRatio = zpos / zpos1;
	 
// chamber prototype
	 tpar[0]= 0.;
	 tpar[1]= 0.;
	 tpar[2]= 0.;
	 
	 char volAlu[5]; // Alu
	 char volBak[5]; // Bakelite
	 char volGaz[5]; // Gas streamer
	 
	 sprintf(volAlu,"SC%dA",icount+1);
	 sprintf(volBak,"SB%dA",icount+1);
	 sprintf(volGaz,"SG%dA",icount+1);
	 
	 gMC->Gsvolu(volAlu,"BOX",idAlu1,tpar,0);           // Al
	 gMC->Gsvolu(volBak,"BOX",idtmed[1107],tpar,0);     // Bakelite
	 gMC->Gsvolu(volGaz,"BOX",idtmed[1106],tpar,0);     // Gas streamer
	 
// chamber type A
	 tpar[0] = -1.;
	 tpar[1] = -1.;
	 
	 Float_t xA=(kDXZERO+kXMED+(kXMAX-kXMED)/2.)*zRatio;
	 Float_t yAm=0.;
	 Float_t yAp=0.;
	 
	 tpar[2] = 0.1;    
	 gMC->Gsposp(volGaz,1,volBak,0.,0.,0.,0,"ONLY",tpar,3);
	 tpar[2] = 0.3;
	 gMC->Gsposp(volBak,1,volAlu,0.,0.,0.,0,"ONLY",tpar,3);
	 
	 tpar[2] = 0.4;
	 tpar[0] = ((kXMAX-kXMED)/2.)*zRatio;
	 tpar[1] = kYMIN*zRatio;
	 
	 gMC->Gsposp(volAlu,iVolNum++,volPlane, -xA,yAm,-kZm,0,"ONLY",tpar,3);
	 gMC->Gsposp(volAlu,iVolNum++,volPlane,  xA,yAp,-kZp,0,"ONLY",tpar,3);
	 gMC->Gsbool(volAlu,volFlange);
	 
// chamber type B    
	 Float_t tpar1save=tpar[1];
	 Float_t y1msave=yAm;
	 Float_t y1psave=yAp;
	 
	 tpar[0] = ((kXMAX-kXMIN)/2.) * zRatio;
	 tpar[1] = ((kYMAX-kYMIN)/2.) * zRatio;
	 
	 Float_t xB=(kDXZERO+kXMIN)*zRatio+tpar[0];
	 Float_t yBp=(y1msave+tpar1save)*zpm+tpar[1];
	 Float_t yBm=(y1psave+tpar1save)*zmp+tpar[1];	 

	 gMC->Gsposp(volAlu,iVolNum++,volPlane, -xB, yBp,-kZp,0,"ONLY",tpar,3);
	 gMC->Gsposp(volAlu,iVolNum++,volPlane,  xB, yBm,-kZm,0,"ONLY",tpar,3);
	 gMC->Gsposp(volAlu,iVolNum++,volPlane, -xB,-yBp,-kZp,0,"ONLY",tpar,3);
	 gMC->Gsposp(volAlu,iVolNum++,volPlane,  xB,-yBm,-kZm,0,"ONLY",tpar,3);
	 
// chamber type C (note : same Z than type B)
	 tpar1save=tpar[1];
	 y1msave=yBm;
	 y1psave=yBp;
	 
	 tpar[0] = (kXMAX/2)*zRatio;
	 tpar[1] = (kYMAX/2)*zRatio;
	 
	 Float_t xC=kDXZERO*zRatio+tpar[0];
	 Float_t yCp=(y1psave+tpar1save)*1.+tpar[1];
	 Float_t yCm=(y1msave+tpar1save)*1.+tpar[1];
	 
	 gMC->Gsposp(volAlu,iVolNum++,volPlane,-xC, yCp,-kZp,0,"ONLY",tpar,3);
	 gMC->Gsposp(volAlu,iVolNum++,volPlane, xC, yCm,-kZm,0,"ONLY",tpar,3);
	 gMC->Gsposp(volAlu,iVolNum++,volPlane,-xC,-yCp,-kZp,0,"ONLY",tpar,3);
	 gMC->Gsposp(volAlu,iVolNum++,volPlane, xC,-yCm,-kZm,0,"ONLY",tpar,3);
	 	 
// chamber type D, E and F (same size)        
	 tpar1save=tpar[1];
	 y1msave=yCm;
	 y1psave=yCp;
	 
	 tpar[0] = (kXMAX/2.)*zRatio;
	 tpar[1] =  kYMIN*zRatio;
	 
	 Float_t xD=kDXZERO*zRatio+tpar[0];
	 Float_t yDp=(y1msave+tpar1save)*zpm+tpar[1];
	 Float_t yDm=(y1psave+tpar1save)*zmp+tpar[1];
	 
	 gMC->Gsposp(volAlu,iVolNum++,volPlane, -xD, yDm,-kZm,0,"ONLY",tpar,3);
	 gMC->Gsposp(volAlu,iVolNum++,volPlane,  xD, yDp,-kZp,0,"ONLY",tpar,3);
	 gMC->Gsposp(volAlu,iVolNum++,volPlane, -xD,-yDm,-kZm,0,"ONLY",tpar,3);
	 gMC->Gsposp(volAlu,iVolNum++,volPlane,  xD,-yDp,-kZp,0,"ONLY",tpar,3);
	 
	 tpar1save=tpar[1];
	 y1msave=yDm;
	 y1psave=yDp;
	 Float_t yEp=(y1msave+tpar1save)*zpm+tpar[1];
	 Float_t yEm=(y1psave+tpar1save)*zmp+tpar[1];
	 
	 gMC->Gsposp(volAlu,iVolNum++,volPlane, -xD, yEp,-kZp,0,"ONLY",tpar,3);
	 gMC->Gsposp(volAlu,iVolNum++,volPlane,  xD, yEm,-kZm,0,"ONLY",tpar,3);
	 gMC->Gsposp(volAlu,iVolNum++,volPlane, -xD,-yEp,-kZp,0,"ONLY",tpar,3);
	 gMC->Gsposp(volAlu,iVolNum++,volPlane,  xD,-yEm,-kZm,0,"ONLY",tpar,3);
	 
	 tpar1save=tpar[1];
	 y1msave=yEm;
	 y1psave=yEp;
	 Float_t yFp=(y1msave+tpar1save)*zpm+tpar[1];
	 Float_t yFm=(y1psave+tpar1save)*zmp+tpar[1];
	 
	 gMC->Gsposp(volAlu,iVolNum++,volPlane, -xD, yFm,-kZm,0,"ONLY",tpar,3);
	 gMC->Gsposp(volAlu,iVolNum++,volPlane,  xD, yFp,-kZp,0,"ONLY",tpar,3);
	 gMC->Gsposp(volAlu,iVolNum++,volPlane, -xD,-yFm,-kZm,0,"ONLY",tpar,3);
	 gMC->Gsposp(volAlu,iVolNum++,volPlane,  xD,-yFp,-kZp,0,"ONLY",tpar,3);

// Positioning plane in ALICE     
	 gMC->Gspos(volPlane,1,"ALIC",0.,0.,zpos,0,"ONLY");
	 
     } // end loop on detection planes
 } // end loop on stations

}

 
//___________________________________________
void AliMUONv3::CreateMaterials()
{
  // *** DEFINITION OF AVAILABLE MUON MATERIALS *** 
  //
  //     Ar-CO2 gas (80%+20%)
    Float_t ag1[3]   = { 39.95,12.01,16. };
    Float_t zg1[3]   = { 18.,6.,8. };
    Float_t wg1[3]   = { .8,.0667,.13333 };
    Float_t dg1      = .001821;
    //
    //     Ar-buthane-freon gas -- trigger chambers 
    Float_t atr1[4]  = { 39.95,12.01,1.01,19. };
    Float_t ztr1[4]  = { 18.,6.,1.,9. };
    Float_t wtr1[4]  = { .56,.1262857,.2857143,.028 };
    Float_t dtr1     = .002599;
    //
    //     Ar-CO2 gas 
    Float_t agas[3]  = { 39.95,12.01,16. };
    Float_t zgas[3]  = { 18.,6.,8. };
    Float_t wgas[3]  = { .74,.086684,.173316 };
    Float_t dgas     = .0018327;
    //
    //     Ar-Isobutane gas (80%+20%) -- tracking 
    Float_t ag[3]    = { 39.95,12.01,1.01 };
    Float_t zg[3]    = { 18.,6.,1. };
    Float_t wg[3]    = { .8,.057,.143 };
    Float_t dg       = .0019596;
    //
    //     Ar-Isobutane-Forane-SF6 gas (49%+7%+40%+4%) -- trigger 
    Float_t atrig[5] = { 39.95,12.01,1.01,19.,32.066 };
    Float_t ztrig[5] = { 18.,6.,1.,9.,16. };
    Float_t wtrig[5] = { .49,1.08,1.5,1.84,0.04 };
    Float_t dtrig    = .0031463;
    //
    //     bakelite 

    Float_t abak[3] = {12.01 , 1.01 , 16.};
    Float_t zbak[3] = {6.     , 1.   , 8.};
    Float_t wbak[3] = {6.     , 6.   , 1.}; 
    Float_t dbak = 1.4;

    Float_t epsil, stmin, deemax, tmaxfd, stemax;

    Int_t iSXFLD   = gAlice->Field()->Integ();
    Float_t sXMGMX = gAlice->Field()->Max();
    //
    // --- Define the various materials for GEANT --- 
    AliMaterial(9, "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
    AliMaterial(10, "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
    AliMaterial(15, "AIR$      ", 14.61, 7.3, .001205, 30423.24, 67500);
    AliMixture(19, "Bakelite$", abak, zbak, dbak, -3, wbak);
    AliMixture(20, "ArC4H10 GAS$", ag, zg, dg, 3, wg);
    AliMixture(21, "TRIG GAS$", atrig, ztrig, dtrig, -5, wtrig);
    AliMixture(22, "ArCO2 80%$", ag1, zg1, dg1, 3, wg1);
    AliMixture(23, "Ar-freon $", atr1, ztr1, dtr1, 4, wtr1);
    AliMixture(24, "ArCO2 GAS$", agas, zgas, dgas, 3, wgas);
    // materials for slat: 
    //     Sensitive area: gas (already defined) 
    //     PCB: copper 
    //     insulating material and frame: vetronite
    //     walls: carbon, rohacell, carbon 
  Float_t aglass[5]={12.01, 28.09, 16.,   10.8,  23.};
  Float_t zglass[5]={ 6.,   14.,    8.,    5.,   11.};
  Float_t wglass[5]={ 0.5,  0.105, 0.355, 0.03,  0.01};
  Float_t dglass=1.74;

  // rohacell: C9 H13 N1 O2
  Float_t arohac[4] = {12.01,  1.01, 14.010, 16.};
  Float_t zrohac[4] = { 6.,    1.,    7.,     8.};
  Float_t wrohac[4] = { 9.,   13.,    1.,     2.};
  Float_t drohac    = 0.03;

  AliMaterial(31, "COPPER$",   63.54,    29.,   8.96,  1.4, 0.);
  AliMixture(32, "Vetronite$",aglass, zglass, dglass,    5, wglass);
  AliMaterial(33, "Carbon$",   12.01,     6.,  2.265, 18.8, 49.9);
  AliMixture(34, "Rohacell$", arohac, zrohac, drohac,   -4, wrohac); 


    epsil  = .001; // Tracking precision, 
    stemax = -1.;  // Maximum displacement for multiple scat 
    tmaxfd = -20.; // Maximum angle due to field deflection 
    deemax = -.3;  // Maximum fractional energy loss, DLS 
    stmin  = -.8;
    //
    //    Air 
    AliMedium(1, "AIR_CH_US         ", 15, 1, iSXFLD, sXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
    //
    //    Aluminum 

    AliMedium(4, "ALU_CH_US          ", 9, 0, iSXFLD, sXMGMX, tmaxfd, fMaxStepAlu, 
	    fMaxDestepAlu, epsil, stmin);
    AliMedium(5, "ALU_CH_US          ", 10, 0, iSXFLD, sXMGMX, tmaxfd, fMaxStepAlu, 
	    fMaxDestepAlu, epsil, stmin);
    //
    //    Ar-isoC4H10 gas 

    AliMedium(6, "AR_CH_US          ", 20, 1, iSXFLD, sXMGMX, tmaxfd, fMaxStepGas, 
	    fMaxDestepGas, epsil, stmin);
//
    //    Ar-Isobuthane-Forane-SF6 gas 

    AliMedium(7, "GAS_CH_TRIGGER    ", 21, 1, iSXFLD, sXMGMX, tmaxfd, stemax, deemax, epsil, stmin);

    AliMedium(8, "BAKE_CH_TRIGGER   ", 19, 0, iSXFLD, sXMGMX, tmaxfd, fMaxStepAlu, 
	    fMaxDestepAlu, epsil, stmin);

    AliMedium(9, "ARG_CO2   ", 22, 1, iSXFLD, sXMGMX, tmaxfd, fMaxStepGas, 
	    fMaxDestepAlu, epsil, stmin);
    // tracking media for slats: check the parameters!! 
    AliMedium(11, "PCB_COPPER        ", 31, 0, iSXFLD, sXMGMX, tmaxfd, 
	      fMaxStepAlu, fMaxDestepAlu, epsil, stmin);
    AliMedium(12, "VETRONITE         ", 32, 0, iSXFLD, sXMGMX, tmaxfd, 
	      fMaxStepAlu, fMaxDestepAlu, epsil, stmin);
    AliMedium(13, "CARBON            ", 33, 0, iSXFLD, sXMGMX, tmaxfd, 
	      fMaxStepAlu, fMaxDestepAlu, epsil, stmin);
    AliMedium(14, "Rohacell          ", 34, 0, iSXFLD, sXMGMX, tmaxfd, 
	      fMaxStepAlu, fMaxDestepAlu, epsil, stmin);
}

//___________________________________________

void AliMUONv3::Init()
{
   // 
   // Initialize Tracking Chambers
   //

   if(fDebug) printf("\n%s: Start Init for version 1 - CPC chamber type\n\n",ClassName());
   Int_t i;
   for (i=0; i<AliMUONConstants::NCh(); i++) {
       ( (AliMUONChamber*) (*fChambers)[i])->Init();
   }
   
   //
   // Set the chamber (sensitive region) GEANT identifier
   ((AliMUONChamber*)(*fChambers)[0])->GetGeometry()->SetSensitiveVolume("S01G");
   ((AliMUONChamber*)(*fChambers)[1])->GetGeometry()->SetSensitiveVolume("S02G");

   ((AliMUONChamber*)(*fChambers)[2])->GetGeometry()->SetSensitiveVolume("S03G");
   ((AliMUONChamber*)(*fChambers)[3])->GetGeometry()->SetSensitiveVolume("S04G");

   ((AliMUONChamber*)(*fChambers)[4])->GetGeometry()->SetSensitiveVolume("S05G");
   ((AliMUONChamber*)(*fChambers)[5])->GetGeometry()->SetSensitiveVolume("S06G");

   ((AliMUONChamber*)(*fChambers)[6])->GetGeometry()->SetSensitiveVolume("S07G");
   ((AliMUONChamber*)(*fChambers)[7])->GetGeometry()->SetSensitiveVolume("S08G");

   ((AliMUONChamber*)(*fChambers)[8])->GetGeometry()->SetSensitiveVolume("S09G");
   ((AliMUONChamber*)(*fChambers)[9])->GetGeometry()->SetSensitiveVolume("S10G");

   ((AliMUONChamber*)(*fChambers)[10])->GetGeometry()->SetSensitiveVolume("SG1A");
   ((AliMUONChamber*)(*fChambers)[11])->GetGeometry()->SetSensitiveVolume("SG2A");
   ((AliMUONChamber*)(*fChambers)[12])->GetGeometry()->SetSensitiveVolume("SG3A");
   ((AliMUONChamber*)(*fChambers)[13])->GetGeometry()->SetSensitiveVolume("SG4A");

   if(fDebug) printf("\n%s: Finished Init for version 1 - CPC chamber type\n",ClassName());

   //cp 
   if(fDebug) printf("\n%s: Start Init for Trigger Circuits\n",ClassName());
   for (i=0; i<AliMUONConstants::NTriggerCircuit(); i++) {
     ( (AliMUONTriggerCircuit*) (*fTriggerCircuits)[i])->Init(i);
   }
   if(fDebug) printf("%s: Finished Init for Trigger Circuits\n",ClassName());
   //cp

}

//_______________________________________________________________________________
Int_t  AliMUONv3::GetChamberId(Int_t volId) const
{
// Check if the volume with specified  volId is a sensitive volume (gas) 
// of some chamber and returns the chamber number;
// if not sensitive volume - return 0.
// ---

  for (Int_t i = 1; i <= AliMUONConstants::NCh(); i++)
    if ( ((AliMUONChamber*)(*fChambers)[i-1])->IsSensId(volId) ) return i;

  return 0;
}
//_______________________________________________________________________________
void AliMUONv3::StepManager()
{
  // Stepmanager for the chambers

 if (fStepManagerVersionOld) {
    StepManagerOld();
    return;
  }

  // Only charged tracks
  if( !(gMC->TrackCharge()) ) return; 
  // Only charged tracks
  
  // Only gas gap inside chamber
  // Tag chambers and record hits when track enters 
  Int_t   idvol=-1;
  Int_t   iChamber=0;
  Int_t   id=0;
  Int_t   copy;
  const  Float_t kBig = 1.e10;

  id=gMC->CurrentVolID(copy);
  iChamber = GetChamberId(id);
  idvol=GetChamberId(id)-1;

  if (idvol == -1) return;

   if( gMC->IsTrackEntering() ) {
     Float_t theta = fTrackMomentum.Theta();
     if ((TMath::Pi()-theta)*kRaddeg>=15.) gMC->SetMaxStep(fStepMaxInActiveGas); // We use Pi-theta because z is negative
  }

//  if (GetDebug()) {
//     Float_t z = ( (AliMUONChamber*)(*fChambers)[idvol])->Z() ;
//      Info("StepManager Step","Active volume found %d chamber %d Z chamber is %f ",idvol,iChamber, z);
//   }  
  // Particule id and mass, 
  Int_t     ipart = gMC->TrackPid();
  Float_t   mass  = gMC->TrackMass();

  fDestepSum[idvol]+=gMC->Edep();
  // Get current particle id (ipart), track position (pos)  and momentum (mom)
  if ( fStepSum[idvol]==0.0 )  gMC->TrackMomentum(fTrackMomentum);
  fStepSum[idvol]+=gMC->TrackStep();
  
//   if (GetDebug()) {
//     Info("StepManager Step","iChamber %d, Particle %d, theta %f phi %f mass %f StepSum %f eloss %g",
// 	 iChamber,ipart, fTrackMomentum.Theta()*kRaddeg, fTrackMomentum.Phi()*kRaddeg, mass, fStepSum[idvol], gMC->Edep());
//     Info("StepManager Step","Track Momentum %f %f %f", fTrackMomentum.X(), fTrackMomentum.Y(), fTrackMomentum.Z()) ;
//     gMC->TrackPosition(fTrackPosition);
//     Info("StepManager Step","Track Position %f %f %f",fTrackPosition.X(),fTrackPosition.Y(),fTrackPosition.Z()) ;
//   }

  // Track left chamber or StepSum larger than fStepMaxInActiveGas
  if ( gMC->IsTrackExiting() || 
       gMC->IsTrackStop() || 
       gMC->IsTrackDisappeared()||
       (fStepSum[idvol]>fStepMaxInActiveGas) ) {
    
    if   ( gMC->IsTrackExiting() || 
	   gMC->IsTrackStop() || 
	   gMC->IsTrackDisappeared() ) gMC->SetMaxStep(kBig);

    gMC->TrackPosition(fTrackPosition);
    Float_t theta = fTrackMomentum.Theta();
    Float_t phi   = fTrackMomentum.Phi();
    
    TLorentzVector backToWire( fStepSum[idvol]/2.*sin(theta)*cos(phi),
			       fStepSum[idvol]/2.*sin(theta)*sin(phi),
			       fStepSum[idvol]/2.*cos(theta),0.0       );
    //     if (GetDebug()) 
    //       Info("StepManager Exit","Track Position %f %f %f",fTrackPosition.X(),fTrackPosition.Y(),fTrackPosition.Z()) ;
    //     if (GetDebug()) 
    //        Info("StepManager Exit ","Track backToWire %f %f %f",backToWire.X(),backToWire.Y(),backToWire.Z()) ;
    fTrackPosition-=backToWire;
    
    //-------------- Angle effect 
    // Ratio between energy loss of particle and Mip as a function of BetaGamma of particle (Energy/Mass)
    
    Float_t betaxGamma    = fTrackMomentum.P()/mass;//  pc/mc2
    Float_t sigmaEffect10degrees;
    Float_t sigmaEffectThetadegrees;
    Float_t eLossParticleELossMip;
    Float_t yAngleEffect=0.;
    Float_t thetawires      =  TMath::Abs( TMath::ASin( TMath::Sin(TMath::Pi()-theta) * TMath::Sin(phi) ) );// We use Pi-theta because z is negative


    if (fAngleEffect){
    if ( (betaxGamma >3.2)   &&  (thetawires*kRaddeg<=15.) ) {
      betaxGamma=TMath::Log(betaxGamma);
      eLossParticleELossMip = fElossRatio->Eval(betaxGamma);
      // 10 degrees is a reference for a model (arbitrary)
      sigmaEffect10degrees=fAngleEffect10->Eval(eLossParticleELossMip);// in micrometers
      // Angle with respect to the wires assuming that chambers are perpendicular to the z axis.
      sigmaEffectThetadegrees =  sigmaEffect10degrees/fAngleEffectNorma->Eval(thetawires*kRaddeg);  // For 5mm gap  
      if ( (iChamber==1)  ||  (iChamber==2) )  
	sigmaEffectThetadegrees/=(1.09833e+00+1.70000e-02*(thetawires*kRaddeg)); // The gap is different (4mm)
      yAngleEffect=1.e-04*gRandom->Gaus(0,sigmaEffectThetadegrees); // Error due to the angle effect in cm
    }
    }
    
    // One hit per chamber
    GetMUONData()->AddHit(fIshunt, gAlice->GetMCApp()->GetCurrentTrackNumber(), iChamber, ipart, 
			  fTrackPosition.X(), fTrackPosition.Y()+yAngleEffect, fTrackPosition.Z(), 0.0, 
			  fTrackMomentum.P(),theta, phi, fStepSum[idvol], fDestepSum[idvol],
			  fTrackPosition.X(),fTrackPosition.Y(),fTrackPosition.Z());
//     if (GetDebug()){
//       Info("StepManager Exit","Particle exiting from chamber %d",iChamber);
//       Info("StepManager Exit","StepSum %f eloss geant %g ",fStepSum[idvol],fDestepSum[idvol]);
//       Info("StepManager Exit","Track Position %f %f %f",fTrackPosition.X(),fTrackPosition.Y(),fTrackPosition.Z()) ;
//     }
    fStepSum[idvol]  =0; // Reset for the next event
    fDestepSum[idvol]=0; // Reset for the next event
  }
}

//__________________________________________
void AliMUONv3::StepManagerOld()
{
  // Old Stepmanager for the chambers
  Int_t          copy, id;
  static Int_t   idvol;
  static Int_t   vol[2];
  Int_t          ipart;
  TLorentzVector pos;
  TLorentzVector mom;
  Float_t        theta,phi;
  Float_t        destep, step;
  
  static Float_t sstep;
  static Float_t eloss, eloss2, xhit, yhit, zhit, tof, tlength;
  const  Float_t kBig = 1.e10;
  static Float_t hits[15];

  TClonesArray &lhits = *fHits;

  //
  //
  // Only charged tracks
  if( !(gMC->TrackCharge()) ) return; 
  //
  // Only gas gap inside chamber
  // Tag chambers and record hits when track enters 
  id=gMC->CurrentVolID(copy);
  vol[0] = GetChamberId(id);
  idvol = vol[0] -1;

  if (idvol == -1) return;

  //
  // Get current particle id (ipart), track position (pos)  and momentum (mom) 
  gMC->TrackPosition(pos);
  gMC->TrackMomentum(mom);

  ipart  = gMC->TrackPid();

  //
  // momentum loss and steplength in last step
  destep = gMC->Edep();
  step   = gMC->TrackStep();
  // cout<<"------------"<<step<<endl;
  //
  // record hits when track enters ...
  if( gMC->IsTrackEntering()) {

      gMC->SetMaxStep(fMaxStepGas);
      Double_t tc = mom[0]*mom[0]+mom[1]*mom[1];
      Double_t rt = TMath::Sqrt(tc);
      Double_t pmom = TMath::Sqrt(tc+mom[2]*mom[2]);
      Double_t tx = mom[0]/pmom;
      Double_t ty = mom[1]/pmom;
      Double_t tz = mom[2]/pmom;
      Double_t s  = ((AliMUONChamber*)(*fChambers)[idvol])
	  ->ResponseModel()
	  ->Pitch()/tz;
      theta   = Float_t(TMath::ATan2(rt,Double_t(mom[2])))*kRaddeg;
      phi     = Float_t(TMath::ATan2(Double_t(mom[1]),Double_t(mom[0])))*kRaddeg;
      hits[0] = Float_t(ipart);         // Geant3 particle type
      hits[1] = pos[0]+s*tx;            // X-position for hit
      hits[2] = pos[1]+s*ty;            // Y-position for hit
      hits[3] = pos[2]+s*tz;            // Z-position for hit
      hits[4] = theta;                  // theta angle of incidence
      hits[5] = phi;                    // phi angle of incidence 
      hits[8] = 0;//PadHits does not exist anymore  (Float_t) fNPadHits;    // first padhit
      hits[9] = -1;                     // last pad hit
      hits[10] = mom[3];                // hit momentum P
      hits[11] = mom[0];                // Px
      hits[12] = mom[1];                // Py
      hits[13] = mom[2];                // Pz
      tof=gMC->TrackTime();
      hits[14] = tof;                   // Time of flight
      tlength  = 0;
      eloss    = 0;
      eloss2   = 0;
      sstep=0;
      xhit     = pos[0];
      yhit     = pos[1];      
      zhit     = pos[2];      
      Chamber(idvol).ChargeCorrelationInit();
      // Only if not trigger chamber

//       printf("---------------------------\n");
//       printf(">>>> Y =  %f \n",hits[2]);
//       printf("---------------------------\n");
    
      

     //  if(idvol < AliMUONConstants::NTrackingCh()) {
// 	  //
// 	  //  Initialize hit position (cursor) in the segmentation model 
// 	  ((AliMUONChamber*) (*fChambers)[idvol])
// 	      ->SigGenInit(pos[0], pos[1], pos[2]);
//       } else {
// 	  //geant3->Gpcxyz();
// 	  //printf("In the Trigger Chamber #%d\n",idvol-9);
//       }
  }
  eloss2+=destep;
  sstep+=step;

  // cout<<sstep<<endl;

  // 
  // Calculate the charge induced on a pad (disintegration) in case 
  //
  // Mip left chamber ...
  if( gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared()){
      gMC->SetMaxStep(kBig);
      eloss   += destep;
      tlength += step;
      
      Float_t x0,y0,z0;
      Float_t localPos[3];
      Float_t globalPos[3] = {pos[0], pos[1], pos[2]};
      gMC->Gmtod(globalPos,localPos,1); 

      if(idvol < AliMUONConstants::NTrackingCh()) {
// tracking chambers
	  x0 = 0.5*(xhit+pos[0]);
	  y0 = 0.5*(yhit+pos[1]);
	  z0 = 0.5*(zhit+pos[2]);
      } else {
// trigger chambers
	  x0 = xhit;
	  y0 = yhit;
	  z0 = 0.;
      }
      

      //      if (eloss >0)  MakePadHits(x0,y0,z0,eloss,tof,idvol);
      
	  
      hits[6] = tlength;   // track length
      hits[7] = eloss2;    // de/dx energy loss


      //      if (fNPadHits > (Int_t)hits[8]) {
      //	  hits[8] = hits[8]+1;
      //	  hits[9] = 0: // PadHits does not exist anymore (Float_t) fNPadHits;
      //}
//
//    new hit 
      
      new(lhits[fNhits++]) 
	  AliMUONHit(fIshunt, gAlice->GetMCApp()->GetCurrentTrackNumber(), vol,hits);
      eloss = 0; 
      //
      // Check additional signal generation conditions 
      // defined by the segmentation
      // model (boundary crossing conditions)
      // only for tracking chambers
  } else if 
      ((idvol < AliMUONConstants::NTrackingCh()) &&
       ((AliMUONChamber*) (*fChambers)[idvol])->SigGenCond(pos[0], pos[1], pos[2]))
  {
      ((AliMUONChamber*) (*fChambers)[idvol])
	  ->SigGenInit(pos[0], pos[1], pos[2]);
      
      Float_t localPos[3];
      Float_t globalPos[3] = {pos[0], pos[1], pos[2]};
      gMC->Gmtod(globalPos,localPos,1); 

      eloss    += destep;

      // if (eloss > 0 && idvol < AliMUONConstants::NTrackingCh())
      //	MakePadHits(0.5*(xhit+pos[0]),0.5*(yhit+pos[1]),pos[2],eloss,tof,idvol);
      xhit     = pos[0];
      yhit     = pos[1]; 
      zhit     = pos[2];
      eloss = 0;
      tlength += step ;
      //
      // nothing special  happened, add up energy loss
  } else {        
      eloss   += destep;
      tlength += step ;
  }
}


