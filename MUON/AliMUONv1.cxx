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
 * about the suitability of this software for any purpeateose. It is      *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log$
Revision 1.14  2000/11/06 09:16:50  morsch
Avoid overlap of slat volumes.

Revision 1.13  2000/10/26 07:33:44  morsch
Correct x-position of slats in station 5.

Revision 1.12  2000/10/25 19:55:35  morsch
Switches for each station individually for debug and lego.

Revision 1.11  2000/10/22 16:44:01  morsch
Update of slat geometry for stations 3,4,5 (A. deFalco)

Revision 1.10  2000/10/12 16:07:04  gosset
StepManager:
* SigGenCond only called for tracking chambers,
  hence no more division by 0,
  and may use last ALIROOT/dummies.C with exception handling;
* "10" replaced by "AliMUONConstants::NTrackingCh()".

Revision 1.9  2000/10/06 15:37:22  morsch
Problems with variable redefinition in for-loop solved.
Variable names starting with u-case letters changed to l-case.

Revision 1.8  2000/10/06 09:06:31  morsch
Include Slat chambers (stations 3-5) into geometry (A. de Falco)

Revision 1.7  2000/10/02 21:28:09  fca
Removal of useless dependecies via forward declarations

Revision 1.6  2000/10/02 17:20:45  egangler
Cleaning of the code (continued ) :
-> coding conventions
-> void Streamers
-> some useless includes removed or replaced by "class" statement

Revision 1.5  2000/06/28 15:16:35  morsch
(1) Client code adapted to new method signatures in AliMUONSegmentation (see comments there)
to allow development of slat-muon chamber simulation and reconstruction code in the MUON
framework. The changes should have no side effects (mostly dummy arguments).
(2) Hit disintegration uses 3-dim hit coordinates to allow simulation
of chambers with overlapping modules (MakePadHits, Disintegration).

Revision 1.4  2000/06/26 14:02:38  morsch
Add class AliMUONConstants with MUON specific constants using static memeber data and access methods.

Revision 1.3  2000/06/22 14:10:05  morsch
HP scope problems corrected (PH)

Revision 1.2  2000/06/15 07:58:49  morsch
Code from MUON-dev joined

Revision 1.1.2.14  2000/06/14 14:37:25  morsch
Initialization of TriggerCircuit added (PC)

Revision 1.1.2.13  2000/06/09 21:55:47  morsch
Most coding rule violations corrected.

Revision 1.1.2.12  2000/05/05 11:34:29  morsch
Log inside comments.

Revision 1.1.2.11  2000/05/05 10:06:48  morsch
Coding Rule violations regarding trigger section corrected (CP)
Log messages included.
*/

/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 0    //
/////////////////////////////////////////////////////////

#include <TTUBE.h>
#include <TNode.h> 
#include <TRandom.h> 
#include <TLorentzVector.h> 
#include <iostream.h>

#include "AliMUONv1.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliMagF.h"
#include "AliCallf77.h"
#include "AliConst.h" 
#include "AliMUONChamber.h"
#include "AliMUONHit.h"
#include "AliMUONPadHit.h"
#include "AliMUONConstants.h"
#include "AliMUONTriggerCircuit.h"

ClassImp(AliMUONv1)
 
//___________________________________________
AliMUONv1::AliMUONv1() : AliMUON()
{
// Constructor
    fChambers = 0;
}
 
//___________________________________________
AliMUONv1::AliMUONv1(const char *name, const char *title)
       : AliMUON(name,title)
{
// Constructor
}

//___________________________________________
void AliMUONv1::CreateGeometry()
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
     Float_t pgpar[10];
     Float_t zpos1, zpos2, zfpos;
     Float_t dframep=.001; // Value for station 3 should be 6 ...
     Float_t dframep1=.001;
//     Bool_t frames=kTRUE;
     Bool_t frames=kFALSE;     
     
     Float_t dframez=0.9;
     Float_t dr;
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
     Int_t idAlu1=idtmed[1103];
     Int_t idAlu2=idtmed[1104];
//     Int_t idAlu1=idtmed[1100];
//     Int_t idAlu2=idtmed[1100];
     Int_t idAir=idtmed[1100];
     Int_t idGas=idtmed[1105];
     

     AliMUONChamber *iChamber, *iChamber1, *iChamber2;
     Int_t stations[5] = {1, 1, 1, 1, 1};
     
     if (stations[0]) {
	 
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
     dstation = zpos2 - zpos1;
     zfpos=-(iChamber->DGas()+dframez+iChamber->DAlu())/2;
     
//
//   Mother volume
     tpar[0] = iChamber->RInner()-dframep1; 
     tpar[1] = (iChamber->ROuter()+dframep1)/TMath::Cos(phi);
     tpar[2] = dstation/5;

     gMC->Gsvolu("C01M", "TUBE", idAir, tpar, 3);
     gMC->Gsvolu("C02M", "TUBE", idAir, tpar, 3);
     gMC->Gspos("C01M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");
     gMC->Gspos("C02M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");     
// Aluminium frames
// Outer frames
     pgpar[0] = 360/12/2;
     pgpar[1] = 360.;
     pgpar[2] = 12.;
     pgpar[3] =   2;
     pgpar[4] = -dframez/2;
     pgpar[5] = iChamber->ROuter();
     pgpar[6] = pgpar[5]+dframep1;
     pgpar[7] = +dframez/2;
     pgpar[8] = pgpar[5];
     pgpar[9] = pgpar[6];
     gMC->Gsvolu("C01O", "PGON", idAlu1, pgpar, 10);
     gMC->Gsvolu("C02O", "PGON", idAlu1, pgpar, 10);
     gMC->Gspos("C01O",1,"C01M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C01O",2,"C01M", 0.,0.,+zfpos,  0,"ONLY");
     gMC->Gspos("C02O",1,"C02M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C02O",2,"C02M", 0.,0.,+zfpos,  0,"ONLY");
//
// Inner frame
     tpar[0]= iChamber->RInner()-dframep1;
     tpar[1]= iChamber->RInner();
     tpar[2]= dframez/2;
     gMC->Gsvolu("C01I", "TUBE", idAlu1, tpar, 3);
     gMC->Gsvolu("C02I", "TUBE", idAlu1, tpar, 3);

     gMC->Gspos("C01I",1,"C01M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C01I",2,"C01M", 0.,0.,+zfpos,  0,"ONLY");
     gMC->Gspos("C02I",1,"C02M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C02I",2,"C02M", 0.,0.,+zfpos,  0,"ONLY");
//
// Frame Crosses
     if (frames) {

	 bpar[0] = (iChamber->ROuter() - iChamber->RInner())/2;
	 bpar[1] = dframep1/2;
	 bpar[2] = dframez/2;
	 gMC->Gsvolu("C01B", "BOX", idAlu1, bpar, 3);
	 gMC->Gsvolu("C02B", "BOX", idAlu1, bpar, 3);
	 
	 gMC->Gspos("C01B",1,"C01M", +iChamber->RInner()+bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C01B",2,"C01M", -iChamber->RInner()-bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C01B",3,"C01M", 0, +iChamber->RInner()+bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C01B",4,"C01M", 0, -iChamber->RInner()-bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C01B",5,"C01M", +iChamber->RInner()+bpar[0] , 0,+zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C01B",6,"C01M", -iChamber->RInner()-bpar[0] , 0,+zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C01B",7,"C01M", 0, +iChamber->RInner()+bpar[0] ,+zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C01B",8,"C01M", 0, -iChamber->RInner()-bpar[0] ,+zfpos, 
		    idrotm[1101],"ONLY");
	 
	 gMC->Gspos("C02B",1,"C02M", +iChamber->RInner()+bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C02B",2,"C02M", -iChamber->RInner()-bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C02B",3,"C02M", 0, +iChamber->RInner()+bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C02B",4,"C02M", 0, -iChamber->RInner()-bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C02B",5,"C02M", +iChamber->RInner()+bpar[0] , 0,+zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C02B",6,"C02M", -iChamber->RInner()-bpar[0] , 0,+zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C02B",7,"C02M", 0, +iChamber->RInner()+bpar[0] ,+zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C02B",8,"C02M", 0, -iChamber->RInner()-bpar[0] ,+zfpos, 
		    idrotm[1101],"ONLY");
     }
//
//   Chamber Material represented by Alu sheet
     tpar[0]= iChamber->RInner();
     tpar[1]= iChamber->ROuter();
     tpar[2] = (iChamber->DGas()+iChamber->DAlu())/2;
     gMC->Gsvolu("C01A", "TUBE",  idAlu2, tpar, 3);
     gMC->Gsvolu("C02A", "TUBE",idAlu2, tpar, 3);
     gMC->Gspos("C01A", 1, "C01M", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("C02A", 1, "C02M", 0., 0., 0.,  0, "ONLY");
//     
//   Sensitive volumes
     // tpar[2] = iChamber->DGas();
     tpar[2] = iChamber->DGas()/2;
     gMC->Gsvolu("C01G", "TUBE", idtmed[1108], tpar, 3);
     gMC->Gsvolu("C02G", "TUBE", idtmed[1108], tpar, 3);
     gMC->Gspos("C01G", 1, "C01A", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("C02G", 1, "C02A", 0., 0., 0.,  0, "ONLY");
//
// Frame Crosses to be placed inside gas 
     if (frames) {

	 dr = (iChamber->ROuter() - iChamber->RInner());
	 bpar[0] = TMath::Sqrt(dr*dr-dframep1*dframep1/4)/2;
	 bpar[1] = dframep1/2;
	 bpar[2] = iChamber->DGas()/2;
	 gMC->Gsvolu("C01F", "BOX", idAlu1, bpar, 3);
	 gMC->Gsvolu("C02F", "BOX", idAlu1, bpar, 3);
	 
	 gMC->Gspos("C01F",1,"C01G", +iChamber->RInner()+bpar[0] , 0, 0, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C01F",2,"C01G", -iChamber->RInner()-bpar[0] , 0, 0, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C01F",3,"C01G", 0, +iChamber->RInner()+bpar[0] , 0, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C01F",4,"C01G", 0, -iChamber->RInner()-bpar[0] , 0, 
		    idrotm[1101],"ONLY");
	 
	 gMC->Gspos("C02F",1,"C02G", +iChamber->RInner()+bpar[0] , 0, 0, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C02F",2,"C02G", -iChamber->RInner()-bpar[0] , 0, 0, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C02F",3,"C02G", 0, +iChamber->RInner()+bpar[0] , 0, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C02F",4,"C02G", 0, -iChamber->RInner()-bpar[0] , 0, 
		    idrotm[1101],"ONLY");
     }
     }
     if (stations[1]) {
	 
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
     dstation = zpos2 - zpos1;
     zfpos=-(iChamber->DGas()+dframez+iChamber->DAlu())/2;
     
//
//   Mother volume
     tpar[0] = iChamber->RInner()-dframep; 
     tpar[1] = (iChamber->ROuter()+dframep)/TMath::Cos(phi);
     tpar[2] = dstation/5;

     gMC->Gsvolu("C03M", "TUBE", idAir, tpar, 3);
     gMC->Gsvolu("C04M", "TUBE", idAir, tpar, 3);
     gMC->Gspos("C03M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");
     gMC->Gspos("C04M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");

// Aluminium frames
// Outer frames
     pgpar[0] = 360/12/2;
     pgpar[1] = 360.;
     pgpar[2] = 12.;
     pgpar[3] =   2;
     pgpar[4] = -dframez/2;
     pgpar[5] = iChamber->ROuter();
     pgpar[6] = pgpar[5]+dframep;
     pgpar[7] = +dframez/2;
     pgpar[8] = pgpar[5];
     pgpar[9] = pgpar[6];
     gMC->Gsvolu("C03O", "PGON", idAlu1, pgpar, 10);
     gMC->Gsvolu("C04O", "PGON", idAlu1, pgpar, 10);
     gMC->Gspos("C03O",1,"C03M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C03O",2,"C03M", 0.,0.,+zfpos,  0,"ONLY");
     gMC->Gspos("C04O",1,"C04M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C04O",2,"C04M", 0.,0.,+zfpos,  0,"ONLY");
//
// Inner frame
     tpar[0]= iChamber->RInner()-dframep;
     tpar[1]= iChamber->RInner();
     tpar[2]= dframez/2;
     gMC->Gsvolu("C03I", "TUBE", idAlu1, tpar, 3);
     gMC->Gsvolu("C04I", "TUBE", idAlu1, tpar, 3);

     gMC->Gspos("C03I",1,"C03M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C03I",2,"C03M", 0.,0.,+zfpos,  0,"ONLY");
     gMC->Gspos("C04I",1,"C04M", 0.,0.,-zfpos,  0,"ONLY");
     gMC->Gspos("C04I",2,"C04M", 0.,0.,+zfpos,  0,"ONLY");
//
// Frame Crosses
     if (frames) {

	 bpar[0] = (iChamber->ROuter() - iChamber->RInner())/2;
	 bpar[1] = dframep/2;
	 bpar[2] = dframez/2;
	 gMC->Gsvolu("C03B", "BOX", idAlu1, bpar, 3);
	 gMC->Gsvolu("C04B", "BOX", idAlu1, bpar, 3);
	 
	 gMC->Gspos("C03B",1,"C03M", +iChamber->RInner()+bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C03B",2,"C03M", -iChamber->RInner()-bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C03B",3,"C03M", 0, +iChamber->RInner()+bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C03B",4,"C03M", 0, -iChamber->RInner()-bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C03B",5,"C03M", +iChamber->RInner()+bpar[0] , 0,+zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C03B",6,"C03M", -iChamber->RInner()-bpar[0] , 0,+zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C03B",7,"C03M", 0, +iChamber->RInner()+bpar[0] ,+zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C03B",8,"C03M", 0, -iChamber->RInner()-bpar[0] ,+zfpos, 
		    idrotm[1101],"ONLY");
	 
	 gMC->Gspos("C04B",1,"C04M", +iChamber->RInner()+bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C04B",2,"C04M", -iChamber->RInner()-bpar[0] , 0,-zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C04B",3,"C04M", 0, +iChamber->RInner()+bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C04B",4,"C04M", 0, -iChamber->RInner()-bpar[0] ,-zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C04B",5,"C04M", +iChamber->RInner()+bpar[0] , 0,+zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C04B",6,"C04M", -iChamber->RInner()-bpar[0] , 0,+zfpos, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C04B",7,"C04M", 0, +iChamber->RInner()+bpar[0] ,+zfpos, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C04B",8,"C04M", 0, -iChamber->RInner()-bpar[0] ,+zfpos, 
		    idrotm[1101],"ONLY");
     }
//
//   Chamber Material represented by Alu sheet
     tpar[0]= iChamber->RInner();
     tpar[1]= iChamber->ROuter();
     tpar[2] = (iChamber->DGas()+iChamber->DAlu())/2;
     gMC->Gsvolu("C03A", "TUBE", idAlu2, tpar, 3);
     gMC->Gsvolu("C04A", "TUBE", idAlu2, tpar, 3);
     gMC->Gspos("C03A", 1, "C03M", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("C04A", 1, "C04M", 0., 0., 0.,  0, "ONLY");
//     
//   Sensitive volumes
     // tpar[2] = iChamber->DGas();
     tpar[2] = iChamber->DGas()/2;
     gMC->Gsvolu("C03G", "TUBE", idGas, tpar, 3);
     gMC->Gsvolu("C04G", "TUBE", idGas, tpar, 3);
     gMC->Gspos("C03G", 1, "C03A", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("C04G", 1, "C04A", 0., 0., 0.,  0, "ONLY");

     if (frames) {
//
// Frame Crosses to be placed inside gas 
	 dr = (iChamber->ROuter() - iChamber->RInner());
	 bpar[0] = TMath::Sqrt(dr*dr-dframep*dframep/4)/2;
	 bpar[1] = dframep/2;
	 bpar[2] = iChamber->DGas()/2;
	 gMC->Gsvolu("C03F", "BOX", idAlu1, bpar, 3);
	 gMC->Gsvolu("C04F", "BOX", idAlu1, bpar, 3);
	 
	 gMC->Gspos("C03F",1,"C03G", +iChamber->RInner()+bpar[0] , 0, 0, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C03F",2,"C03G", -iChamber->RInner()-bpar[0] , 0, 0, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C03F",3,"C03G", 0, +iChamber->RInner()+bpar[0] , 0, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C03F",4,"C03G", 0, -iChamber->RInner()-bpar[0] , 0, 
		    idrotm[1101],"ONLY");
	 
	 gMC->Gspos("C04F",1,"C04G", +iChamber->RInner()+bpar[0] , 0, 0, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C04F",2,"C04G", -iChamber->RInner()-bpar[0] , 0, 0, 
		    idrotm[1100],"ONLY");
	 gMC->Gspos("C04F",3,"C04G", 0, +iChamber->RInner()+bpar[0] , 0, 
		    idrotm[1101],"ONLY");
	 gMC->Gspos("C04F",4,"C04G", 0, -iChamber->RInner()-bpar[0] , 0, 
		    idrotm[1101],"ONLY");
     }
     }
     // define the id of tracking media:
     Int_t idCopper = idtmed[1110];
     Int_t idGlass  = idtmed[1111];
     Int_t idCarbon = idtmed[1112];
     Int_t idRoha   = idtmed[1113];

      // sensitive area: 40*40 cm**2
     const Float_t sensLength = 40.; 
     const Float_t sensHeight = 40.; 
     const Float_t sensWidth  = 0.5; // according to TDR fig 2.120 
     const Int_t sensMaterial = idGas;
     const Float_t yOverlap   = 1.5; 

     // PCB dimensions in cm; width: 30 mum copper   
     const Float_t pcbLength  = sensLength; 
     const Float_t pcbHeight  = 60.; 
     const Float_t pcbWidth   = 0.003;   
     const Int_t pcbMaterial  = idCopper;

     // Insulating material: 200 mum glass fiber glued to pcb  
     const Float_t insuLength = pcbLength; 
     const Float_t insuHeight = pcbHeight; 
     const Float_t insuWidth  = 0.020;   
     const Int_t insuMaterial = idGlass;

     // Carbon fiber panels: 200mum carbon/epoxy skin   
     const Float_t panelLength = sensLength; 
     const Float_t panelHeight = sensHeight; 
     const Float_t panelWidth  = 0.020;      
     const Int_t panelMaterial = idCarbon;

     // rohacell between the two carbon panels   
     const Float_t rohaLength = sensLength; 
     const Float_t rohaHeight = sensHeight; 
     const Float_t rohaWidth  = 0.5;
     const Int_t rohaMaterial = idRoha;

     // Frame around the slat: 2 sticks along length,2 along height  
     // H: the horizontal ones 
     const Float_t hFrameLength = pcbLength; 
     const Float_t hFrameHeight = 1.5; 
     const Float_t hFrameWidth  = sensWidth; 
     const Int_t hFrameMaterial = idGlass;

     // V: the vertical ones 
     const Float_t vFrameLength = 4.0; 
     const Float_t vFrameHeight = sensHeight + hFrameHeight; 
     const Float_t vFrameWidth  = sensWidth;
     const Int_t vFrameMaterial = idGlass;

     // B: the horizontal border filled with rohacell 
     const Float_t bFrameLength = hFrameLength; 
     const Float_t bFrameHeight = (pcbHeight - sensHeight)/2. - hFrameHeight; 
     const Float_t bFrameWidth  = hFrameWidth;
     const Int_t bFrameMaterial = idRoha;

     // NULOC: 30 mum copper + 200 mum vetronite (same radiation length as 14mum copper)
     const Float_t nulocLength = 2.5; 
     const Float_t nulocHeight = 7.5; 
     const Float_t nulocWidth  = 0.0030 + 0.0014; // equivalent copper width of vetronite; 
     const Int_t   nulocMaterial = idCopper;

     // Gassiplex package 
     const Float_t gassiLength   = 1.0; 
     const Float_t gassiHeight   = 1.0; 
     const Float_t gassiWidth    = 0.15; // check it !!!
     const Int_t   gassiMaterial = idGlass; 

     const Float_t slatHeight = pcbHeight; 
     const Float_t slatWidth = sensWidth + 2.*(pcbWidth + insuWidth + 
					       2.* panelWidth + rohaWidth);
     const Int_t slatMaterial = idAir;
     const Float_t dSlatLength = vFrameLength; // border on left and right 

     Float_t spar[3];  
     Int_t i, j;

     Float_t sensPar[3] = { sensLength/2., sensHeight/2., sensWidth/2. }; 
     Float_t pcbpar[3] = { pcbLength/2., pcbHeight/2., pcbWidth/2. }; 
     Float_t insupar[3] = { insuLength/2., insuHeight/2., insuWidth/2. }; 
     Float_t panelpar[3] = { panelLength/2., panelHeight/2., panelWidth/2. }; 
     Float_t rohapar[3] = { rohaLength/2., rohaHeight/2., rohaWidth/2. }; 
     Float_t vFramepar[3]={vFrameLength/2., vFrameHeight/2., vFrameWidth/2.}; 
     Float_t hFramepar[3]={hFrameLength/2., hFrameHeight/2., hFrameWidth/2.}; 
     Float_t bFramepar[3]={bFrameLength/2., bFrameHeight/2., bFrameWidth/2.}; 
     Float_t nulocpar[3]={nulocLength/2., nulocHeight/2., nulocWidth/2.}; 
     Float_t gassipar[3]={gassiLength/2., gassiHeight/2., gassiWidth/2.}; 
     Float_t xx;
     Float_t xxmax = (bFrameLength - nulocLength)/2.; 
     Int_t index=0;
     
     if (stations[2]) {
	 
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
     dstation = zpos2 - zpos1;

     zfpos=-(iChamber->DGas()+dframez+iChamber->DAlu())/2;
//
//   Mother volume
     tpar[0] = iChamber->RInner()-dframep; 
     tpar[1] = (iChamber->ROuter()+dframep)/TMath::Cos(phi);
     tpar[2] = dstation/5;
     gMC->Gsvolu("C05M", "TUBE", idAir, tpar, 3);
     gMC->Gsvolu("C06M", "TUBE", idAir, tpar, 3);
     gMC->Gspos("C05M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");
     gMC->Gspos("C06M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");
 
     // volumes for slat geometry (xx=5,..,10 chamber id): 
     // Sxx0 Sxx1 Sxx2 Sxx3  -->   Slat Mother volumes 
     // SxxG                          -->   Sensitive volume (gas)
     // SxxP                          -->   PCB (copper) 
     // SxxI                          -->   Insulator (vetronite) 
     // SxxC                          -->   Carbon panel 
     // SxxR                          -->   Rohacell
     // SxxH, SxxV                    -->   Horizontal and Vertical frames (vetronite)

     // slat dimensions: slat is a MOTHER volume!!! made of air


     const Int_t nSlats3 = 4;  // number of slats per quadrant
     const Int_t nPCB3[nSlats3] = {3,4,3,2}; // n PCB per slat
     Float_t slatLength3[nSlats3]; 

     // create and position the slat (mother) volumes 

     char volNam5[5];
     char volDiv5[5];
     char volNam6[5];
     char volDiv6[5];
     Float_t xSlat3;

     for (i = 0; i<nSlats3; i++){
	 slatLength3[i] = pcbLength * nPCB3[i] + 2. * dSlatLength; 
	 xSlat3 = slatLength3[i]/2. - vFrameLength/2.; 
	 if (i==0) xSlat3 += 40.;
	 
	 Float_t ySlat31 =  sensHeight * (i+0.5) - yOverlap * i - yOverlap/2.; 
	 Float_t ySlat32 = -sensHeight * (i+0.5) + yOverlap * i + yOverlap/2.; 
	 spar[0] = slatLength3[i]/2.; 
	 spar[1] = slatHeight/2.;
	 spar[2] = slatWidth/2.; 
	 // zSlat to be checked (odd downstream or upstream?)
	 Float_t zSlat = (i%2 ==0)? -slatWidth/2. : slatWidth/2.; 
	 zSlat*=1.01;
	 sprintf(volNam5,"S05%d",i);
	 gMC->Gsvolu(volNam5,"BOX",slatMaterial,spar,3);
	 gMC->Gspos(volNam5, i*4+1,"C05M", xSlat3, ySlat31, zSlat, 0, "ONLY");
	 gMC->Gspos(volNam5, i*4+2,"C05M",-xSlat3, ySlat31, zSlat, 0, "ONLY");
	 gMC->Gspos(volNam5, i*4+3,"C05M", xSlat3, ySlat32,-zSlat, 0, "ONLY");
	 gMC->Gspos(volNam5, i*4+4,"C05M",-xSlat3, ySlat32,-zSlat, 0, "ONLY");
	 sprintf(volNam6,"S06%d",i);
	 gMC->Gsvolu(volNam6,"BOX",slatMaterial,spar,3);
	 gMC->Gspos(volNam6, i*4+1,"C06M", xSlat3, ySlat31, zSlat, 0, "ONLY");
	 gMC->Gspos(volNam6, i*4+2,"C06M",-xSlat3, ySlat31, zSlat, 0, "ONLY");
	 gMC->Gspos(volNam6, i*4+3,"C06M", xSlat3, ySlat32,-zSlat, 0, "ONLY");
	 gMC->Gspos(volNam6, i*4+4,"C06M",-xSlat3, ySlat32,-zSlat, 0, "ONLY");
	 // 1st pcb in 1st slat made by some rectangular divisions
/*
	 if (i==0) {                        
	     Int_t ndiv=8;
	     Double_t dydiv= sensHeight/ndiv;
	     Double_t ydiv = -dydiv; 
	     for (Int_t idiv=0;idiv<ndiv; idiv++){ 
		 ydiv+= dydiv;
		 Float_t xdiv =0; 
		 if (ydiv<30) xdiv= 30. * TMath::Sin( TMath::ACos(ydiv/30.) );
		 spar[0] = (pcbLength-xdiv)/2.; 
		 spar[1] = dydiv/2.;
		 spar[2] = slatWidth/2.; 
		 
		 sprintf(volDiv5,"D05%d",idiv);
		 sprintf(volDiv6,"D06%d",idiv);
		 
		 gMC->Gsvolu(volDiv5,"BOX",sensMaterial,spar,3);
		 Float_t xvol=(pcbLength+xdiv)/2.;
		 Float_t yvol=ydiv+dydiv/2.;
		 gMC->Gspos(volDiv5, 1,"C05M", xvol, yvol, zSlat, 0, "ONLY");
		 gMC->Gspos(volDiv5, 2,"C05M",-xvol, yvol, zSlat, 0, "ONLY");
		 gMC->Gspos(volDiv5, 3,"C05M", xvol,-yvol,-zSlat, 0, "ONLY");
		 gMC->Gspos(volDiv5, 4,"C05M",-xvol,-yvol,-zSlat, 0, "ONLY");
		 gMC->Gspos(volDiv6, 1,"C06M", xvol, yvol, zSlat, 0, "ONLY");
		 gMC->Gspos(volDiv6, 2,"C06M",-xvol, yvol, zSlat, 0, "ONLY");
		 gMC->Gspos(volDiv6, 3,"C06M", xvol,-yvol,-zSlat, 0, "ONLY");
		 gMC->Gspos(volDiv6, 4,"C06M",-xvol,-yvol,-zSlat, 0, "ONLY");
	     }
	 }
 */
     }

     // create the sensitive volumes (subdivided as the PCBs),

     gMC->Gsvolu("S05G","BOX",sensMaterial,sensPar,3);
     gMC->Gsvolu("S06G","BOX",sensMaterial,sensPar,3);

     // create the PCB volume 

     gMC->Gsvolu("S05P","BOX",pcbMaterial,pcbpar,3);
     gMC->Gsvolu("S06P","BOX",pcbMaterial,pcbpar,3);
 
     // create the insulating material volume 

     gMC->Gsvolu("S05I","BOX",insuMaterial,insupar,3);
     gMC->Gsvolu("S06I","BOX",insuMaterial,insupar,3);

     // create the panel volume 
 
     gMC->Gsvolu("S05C","BOX",panelMaterial,panelpar,3);
     gMC->Gsvolu("S06C","BOX",panelMaterial,panelpar,3);

     // create the rohacell volume 

     gMC->Gsvolu("S05R","BOX",rohaMaterial,rohapar,3);
     gMC->Gsvolu("S06R","BOX",rohaMaterial,rohapar,3);

     // create the vertical frame volume 

     gMC->Gsvolu("S05V","BOX",vFrameMaterial,vFramepar,3);
     gMC->Gsvolu("S06V","BOX",vFrameMaterial,vFramepar,3);

     // create the horizontal frame volume 

     gMC->Gsvolu("S05H","BOX",hFrameMaterial,hFramepar,3);
     gMC->Gsvolu("S06H","BOX",hFrameMaterial,hFramepar,3);

     // create the horizontal border volume 

     gMC->Gsvolu("S05B","BOX",bFrameMaterial,bFramepar,3);
     gMC->Gsvolu("S06B","BOX",bFrameMaterial,bFramepar,3);

     index=0; 
     for (i = 0; i<nSlats3; i++){
       sprintf(volNam5,"S05%d",i);
       sprintf(volNam6,"S06%d",i);
       Float_t xvFrame  = (slatLength3[i] - vFrameLength)/2.;
       gMC->Gspos("S05V",2*i-1,volNam5, xvFrame, 0., 0. , 0, "ONLY");
       gMC->Gspos("S05V",2*i  ,volNam5,-xvFrame, 0., 0. , 0, "ONLY");
       gMC->Gspos("S06V",2*i-1,volNam6, xvFrame, 0., 0. , 0, "ONLY");
       gMC->Gspos("S06V",2*i  ,volNam6,-xvFrame, 0., 0. , 0, "ONLY");
       for (j=0; j<nPCB3[i]; j++){
	 index++;
	 Float_t xx = sensLength * (-nPCB3[i]/2.+j+.5); 
	 Float_t yy = 0.;
	 Float_t zSens = 0.;  
	 gMC->Gspos("S05G",index,volNam5, xx, yy, zSens , 0, "ONLY");
	 gMC->Gspos("S06G",index,volNam6, xx, yy, zSens , 0, "ONLY");
	 Float_t zPCB = (sensWidth+pcbWidth)/2.; 
	 gMC->Gspos("S05P",2*index-1,volNam5, xx, yy, zPCB , 0, "ONLY");
	 gMC->Gspos("S05P",2*index  ,volNam5, xx, yy,-zPCB , 0, "ONLY");
	 gMC->Gspos("S06P",2*index-1,volNam6, xx, yy, zPCB , 0, "ONLY");
	 gMC->Gspos("S06P",2*index  ,volNam6, xx, yy,-zPCB , 0, "ONLY");
	 Float_t zInsu = (insuWidth+pcbWidth)/2. + zPCB; 
	 gMC->Gspos("S05I",2*index-1,volNam5, xx, yy, zInsu , 0, "ONLY");
	 gMC->Gspos("S05I",2*index  ,volNam5, xx, yy,-zInsu , 0, "ONLY");
	 gMC->Gspos("S06I",2*index-1,volNam6, xx, yy, zInsu , 0, "ONLY");
	 gMC->Gspos("S06I",2*index  ,volNam6, xx, yy,-zInsu , 0, "ONLY");
	 Float_t zPanel1 = (insuWidth+panelWidth)/2. + zInsu; 
	 gMC->Gspos("S05C",4*index-3,volNam5, xx, yy, zPanel1 , 0, "ONLY");
	 gMC->Gspos("S05C",4*index-2,volNam5, xx, yy,-zPanel1 , 0, "ONLY");
	 gMC->Gspos("S06C",4*index-3,volNam6, xx, yy, zPanel1 , 0, "ONLY");
	 gMC->Gspos("S06C",4*index-2,volNam6, xx, yy,-zPanel1 , 0, "ONLY");
	 Float_t zRoha = (rohaWidth+panelWidth)/2. + zPanel1; 
	 gMC->Gspos("S05R",2*index-1,volNam5, xx, yy, zRoha , 0, "ONLY");
	 gMC->Gspos("S05R",2*index  ,volNam5, xx, yy,-zRoha , 0, "ONLY");
	 gMC->Gspos("S06R",2*index-1,volNam6, xx, yy, zRoha , 0, "ONLY");
	 gMC->Gspos("S06R",2*index  ,volNam6, xx, yy,-zRoha , 0, "ONLY");
	 Float_t zPanel2 = (rohaWidth+panelWidth)/2. + zRoha; 
	 gMC->Gspos("S05C",4*index-1,volNam5, xx, yy, zPanel2 , 0, "ONLY");
	 gMC->Gspos("S05C",4*index  ,volNam5, xx, yy,-zPanel2 , 0, "ONLY");
	 gMC->Gspos("S06C",4*index-1,volNam6, xx, yy, zPanel2 , 0, "ONLY");
	 gMC->Gspos("S06C",4*index  ,volNam6, xx, yy,-zPanel2 , 0, "ONLY");
	 Float_t yframe = (sensHeight + hFrameHeight)/2.;
	 gMC->Gspos("S05H",2*index-1,volNam5, xx, yframe, 0. , 0, "ONLY");
	 gMC->Gspos("S05H",2*index  ,volNam5, xx,-yframe, 0. , 0, "ONLY");
	 gMC->Gspos("S06H",2*index-1,volNam6, xx, yframe, 0. , 0, "ONLY");
	 gMC->Gspos("S06H",2*index  ,volNam6, xx,-yframe, 0. , 0, "ONLY");
	 Float_t yborder = (bFrameHeight + hFrameHeight)/2. + yframe;
	 gMC->Gspos("S05B",2*index-1,volNam5, xx, yborder, 0. , 0, "ONLY");
	 gMC->Gspos("S05B",2*index  ,volNam5, xx,-yborder, 0. , 0, "ONLY");
	 gMC->Gspos("S06B",2*index-1,volNam6, xx, yborder, 0. , 0, "ONLY");
	 gMC->Gspos("S06B",2*index  ,volNam6, xx,-yborder, 0. , 0, "ONLY");
       } 
     }

     // create the NULOC volume and position it in the horizontal frame

     gMC->Gsvolu("S05N","BOX",nulocMaterial,nulocpar,3);
     gMC->Gsvolu("S06N","BOX",nulocMaterial,nulocpar,3);


     index = 0;


     for (xx = -xxmax; xx<=xxmax; xx+=3*nulocLength) { 
       index++; 
       gMC->Gspos("S05N",2*index-1,"S05B", xx, 0.,-bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S05N",2*index  ,"S05B", xx, 0., bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S06N",2*index-1,"S06B", xx, 0.,-bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S06N",2*index  ,"S06B", xx, 0., bFrameWidth/4., 0, "ONLY");
     }

     // create the gassiplex volume 

     gMC->Gsvolu("S05E","BOX",gassiMaterial,gassipar,3);
     gMC->Gsvolu("S06E","BOX",gassiMaterial,gassipar,3);


     // position 4 gassiplex in the nuloc

     gMC->Gspos("S05E",1,"S05N", 0., -3 * nulocHeight/8., 0. , 0, "ONLY");
     gMC->Gspos("S05E",2,"S05N", 0.,    - nulocHeight/8., 0. , 0, "ONLY");
     gMC->Gspos("S05E",3,"S05N", 0.,      nulocHeight/8., 0. , 0, "ONLY");
     gMC->Gspos("S05E",4,"S05N", 0.,  3 * nulocHeight/8., 0. , 0, "ONLY");
     gMC->Gspos("S06E",1,"S06N", 0., -3 * nulocHeight/8., 0. , 0, "ONLY");
     gMC->Gspos("S06E",2,"S06N", 0.,    - nulocHeight/8., 0. , 0, "ONLY");
     gMC->Gspos("S06E",3,"S06N", 0.,      nulocHeight/8., 0. , 0, "ONLY");
     gMC->Gspos("S06E",4,"S06N", 0.,  3 * nulocHeight/8., 0. , 0, "ONLY");
     }
 if (stations[3]) {
     

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
     dstation = zpos2 - zpos1;
     zfpos=-(iChamber->DGas()+dframez+iChamber->DAlu())/2;
     
//
//   Mother volume
     tpar[0] = iChamber->RInner()-dframep; 
     tpar[1] = (iChamber->ROuter()+dframep)/TMath::Cos(phi);
     tpar[2] = dstation/5;

     gMC->Gsvolu("C07M", "TUBE", idAir, tpar, 3);
     gMC->Gsvolu("C08M", "TUBE", idAir, tpar, 3);
     gMC->Gspos("C07M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");
     gMC->Gspos("C08M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");
     

     const Int_t nSlats4 = 6;  // number of slats per quadrant
     const Int_t nPCB4[nSlats4] = {4,5,5,4,3,2}; // n PCB per slat

     // slat dimensions: slat is a MOTHER volume!!! made of air
     Float_t slatLength4[nSlats4];     

     // create and position the slat (mother) volumes 

     char volNam7[5];
     char volNam8[5];
     Float_t xSlat4;
     Float_t ySlat4;


     for (i = 0; i<nSlats4; i++){
	 slatLength4[i] = pcbLength * nPCB4[i] + 2. * dSlatLength; 
	 xSlat4 = slatLength4[i]/2. - vFrameLength/2.; 
	 if (i==0) xSlat4 += 37.5;
	 ySlat4 =  sensHeight * i - yOverlap *i;
	 
	 spar[0] = slatLength4[i]/2.; 
	 spar[1] = slatHeight/2.;
	 spar[2] = slatWidth/2.; 
	 // zSlat to be checked (odd downstream or upstream?)
	 Float_t zSlat = (i%2 ==0)? slatWidth/2. : -slatWidth/2.; 
	 zSlat*=1.01;
	 sprintf(volNam7,"S07%d",i);
	 gMC->Gsvolu(volNam7,"BOX",slatMaterial,spar,3);
	 gMC->Gspos(volNam7, i*4+1,"C07M", xSlat4, ySlat4, zSlat, 0, "ONLY");
	 gMC->Gspos(volNam7, i*4+2,"C07M",-xSlat4, ySlat4, zSlat, 0, "ONLY");
	 if (i>0) { 
	   gMC->Gspos(volNam7, i*4+3,"C07M", xSlat4,-ySlat4, zSlat, 0, "ONLY");
	   gMC->Gspos(volNam7, i*4+4,"C07M",-xSlat4,-ySlat4, zSlat, 0, "ONLY");
	 }
	 sprintf(volNam8,"S08%d",i);
	 gMC->Gsvolu(volNam8,"BOX",slatMaterial,spar,3);
	 gMC->Gspos(volNam8, i*4+1,"C08M", xSlat4, ySlat4, zSlat, 0, "ONLY");
	 gMC->Gspos(volNam8, i*4+2,"C08M",-xSlat4, ySlat4, zSlat, 0, "ONLY");
	 if (i>0) { 
	   gMC->Gspos(volNam8, i*4+3,"C08M", xSlat4,-ySlat4, zSlat, 0, "ONLY");
	   gMC->Gspos(volNam8, i*4+4,"C08M",-xSlat4,-ySlat4, zSlat, 0, "ONLY");
	 }
     }

     // create the sensitive volumes (subdivided as the PCBs),
 
     gMC->Gsvolu("S07G","BOX",sensMaterial,sensPar,3);
     gMC->Gsvolu("S08G","BOX",sensMaterial,sensPar,3);

     // create the PCB volume 

     gMC->Gsvolu("S07P","BOX",pcbMaterial,pcbpar,3);
     gMC->Gsvolu("S08P","BOX",pcbMaterial,pcbpar,3);
 
     // create the insulating material volume 

     gMC->Gsvolu("S07I","BOX",insuMaterial,insupar,3);
     gMC->Gsvolu("S08I","BOX",insuMaterial,insupar,3);

     // create the panel volume 

     gMC->Gsvolu("S07C","BOX",panelMaterial,panelpar,3);
     gMC->Gsvolu("S08C","BOX",panelMaterial,panelpar,3);

     // create the rohacell volume 
 
     gMC->Gsvolu("S07R","BOX",rohaMaterial,rohapar,3);
     gMC->Gsvolu("S08R","BOX",rohaMaterial,rohapar,3);

     // create the vertical frame volume 

     gMC->Gsvolu("S07V","BOX",vFrameMaterial,vFramepar,3);
     gMC->Gsvolu("S08V","BOX",vFrameMaterial,vFramepar,3);

     // create the horizontal frame volume 

     gMC->Gsvolu("S07H","BOX",hFrameMaterial,hFramepar,3);
     gMC->Gsvolu("S08H","BOX",hFrameMaterial,hFramepar,3);

     // create the horizontal border volume 

     gMC->Gsvolu("S07B","BOX",bFrameMaterial,bFramepar,3);
     gMC->Gsvolu("S08B","BOX",bFrameMaterial,bFramepar,3);
 
     for (i = 0; i<nSlats4; i++){
       sprintf(volNam7,"S07%d",i);
       sprintf(volNam8,"S08%d",i);
       Float_t xvFrame  = (slatLength4[i] - vFrameLength)/2.;
       gMC->Gspos("S07V",2*i-1,volNam7, xvFrame, 0., 0. , 0, "ONLY");
       gMC->Gspos("S07V",2*i  ,volNam7,-xvFrame, 0., 0. , 0, "ONLY");
       gMC->Gspos("S08V",2*i-1,volNam8, xvFrame, 0., 0. , 0, "ONLY");
       gMC->Gspos("S08V",2*i  ,volNam8,-xvFrame, 0., 0. , 0, "ONLY");
       for (j=0; j<nPCB4[i]; j++){
	 index++;
	 Float_t xx = sensLength * (-nPCB4[i]/2.+j+.5); 
	 Float_t yy = 0.;
	 Float_t zSens = 0.;  
	 gMC->Gspos("S07G",index,volNam7, xx, yy, zSens , 0, "ONLY");
	 gMC->Gspos("S08G",index,volNam8, xx, yy, zSens , 0, "ONLY");
	 Float_t zPCB = (sensWidth+pcbWidth)/2.; 
	 gMC->Gspos("S07P",2*index-1,volNam7, xx, yy, zPCB , 0, "ONLY");
	 gMC->Gspos("S07P",2*index  ,volNam7, xx, yy,-zPCB , 0, "ONLY");
	 gMC->Gspos("S08P",2*index-1,volNam8, xx, yy, zPCB , 0, "ONLY");
	 gMC->Gspos("S08P",2*index  ,volNam8, xx, yy,-zPCB , 0, "ONLY");
	 Float_t zInsu = (insuWidth+pcbWidth)/2. + zPCB; 
	 gMC->Gspos("S07I",2*index-1,volNam7, xx, yy, zInsu , 0, "ONLY");
	 gMC->Gspos("S07I",2*index  ,volNam7, xx, yy,-zInsu , 0, "ONLY");
	 gMC->Gspos("S08I",2*index-1,volNam8, xx, yy, zInsu , 0, "ONLY");
	 gMC->Gspos("S08I",2*index  ,volNam8, xx, yy,-zInsu , 0, "ONLY");
	 Float_t zPanel1 = (insuWidth+panelWidth)/2. + zInsu; 
	 gMC->Gspos("S07C",4*index-3,volNam7, xx, yy, zPanel1 , 0, "ONLY");
	 gMC->Gspos("S07C",4*index-2,volNam7, xx, yy,-zPanel1 , 0, "ONLY");
	 gMC->Gspos("S08C",4*index-3,volNam8, xx, yy, zPanel1 , 0, "ONLY");
	 gMC->Gspos("S08C",4*index-2,volNam8, xx, yy,-zPanel1 , 0, "ONLY");
	 Float_t zRoha = (rohaWidth+panelWidth)/2. + zPanel1; 
	 gMC->Gspos("S07R",2*index-1,volNam7, xx, yy, zRoha , 0, "ONLY");
	 gMC->Gspos("S07R",2*index  ,volNam7, xx, yy,-zRoha , 0, "ONLY");
	 gMC->Gspos("S08R",2*index-1,volNam8, xx, yy, zRoha , 0, "ONLY");
	 gMC->Gspos("S08R",2*index  ,volNam8, xx, yy,-zRoha , 0, "ONLY");
	 Float_t zPanel2 = (rohaWidth+panelWidth)/2. + zRoha; 
	 gMC->Gspos("S07C",4*index-1,volNam7, xx, yy, zPanel2 , 0, "ONLY");
	 gMC->Gspos("S07C",4*index  ,volNam7, xx, yy,-zPanel2 , 0, "ONLY");
	 gMC->Gspos("S08C",4*index-1,volNam8, xx, yy, zPanel2 , 0, "ONLY");
	 gMC->Gspos("S08C",4*index  ,volNam8, xx, yy,-zPanel2 , 0, "ONLY");
	 Float_t yframe = (sensHeight + hFrameHeight)/2.;
	 gMC->Gspos("S07H",2*index-1,volNam7, xx, yframe, 0. , 0, "ONLY");
	 gMC->Gspos("S07H",2*index  ,volNam7, xx,-yframe, 0. , 0, "ONLY");
	 gMC->Gspos("S08H",2*index-1,volNam8, xx, yframe, 0. , 0, "ONLY");
	 gMC->Gspos("S08H",2*index  ,volNam8, xx,-yframe, 0. , 0, "ONLY");
	 Float_t yborder = (bFrameHeight + hFrameHeight)/2. + yframe;
	 gMC->Gspos("S07B",2*index-1,volNam7, xx, yborder, 0. , 0, "ONLY");
	 gMC->Gspos("S07B",2*index  ,volNam7, xx,-yborder, 0. , 0, "ONLY");
	 gMC->Gspos("S08B",2*index-1,volNam8, xx, yborder, 0. , 0, "ONLY");
	 gMC->Gspos("S08B",2*index  ,volNam8, xx,-yborder, 0. , 0, "ONLY");
       } 
     }

     // create the NULOC volume and position it in the horizontal frame
 
     gMC->Gsvolu("S07N","BOX",nulocMaterial,nulocpar,3);
     gMC->Gsvolu("S08N","BOX",nulocMaterial,nulocpar,3);


     index = 0; 
     for (xx = -xxmax; xx<=xxmax; xx+=3*nulocLength) { 
       index++; 
       gMC->Gspos("S07N",2*index-1,"S07B", xx, 0.,-bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S07N",2*index  ,"S07B", xx, 0., bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S08N",2*index-1,"S08B", xx, 0.,-bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S08N",2*index  ,"S08B", xx, 0., bFrameWidth/4., 0, "ONLY");
     }

     // create the gassiplex volume 

     gMC->Gsvolu("S07E","BOX",gassiMaterial,gassipar,3);
     gMC->Gsvolu("S08E","BOX",gassiMaterial,gassipar,3);


     // position 4 gassiplex in the nuloc

     gMC->Gspos("S07E",1,"S07N", 0., -3 * nulocHeight/8., 0. , 0, "ONLY");
     gMC->Gspos("S07E",2,"S07N", 0.,    - nulocHeight/8., 0. , 0, "ONLY");
     gMC->Gspos("S07E",3,"S07N", 0.,      nulocHeight/8., 0. , 0, "ONLY");
     gMC->Gspos("S07E",4,"S07N", 0.,  3 * nulocHeight/8., 0. , 0, "ONLY");
     gMC->Gspos("S08E",1,"S08N", 0., -3 * nulocHeight/8., 0. , 0, "ONLY");
     gMC->Gspos("S08E",2,"S08N", 0.,    - nulocHeight/8., 0. , 0, "ONLY");
     gMC->Gspos("S08E",3,"S08N", 0.,      nulocHeight/8., 0. , 0, "ONLY");
     gMC->Gspos("S08E",4,"S08N", 0.,  3 * nulocHeight/8., 0. , 0, "ONLY");
     
 }
 if (stations[4]) {
     

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
     dstation = zpos2 - zpos1;
     zfpos=-(iChamber->DGas()+dframez+iChamber->DAlu())/2;
     
//
//   Mother volume
     tpar[0] = iChamber->RInner()-dframep; 
     tpar[1] = (iChamber->ROuter()+dframep)/TMath::Cos(phi);
     tpar[2] = dstation/5;

     gMC->Gsvolu("C09M", "TUBE", idAir, tpar, 3);
     gMC->Gsvolu("C10M", "TUBE", idAir, tpar, 3);
     gMC->Gspos("C09M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");
     gMC->Gspos("C10M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");


     const Int_t nSlats5 = 7;  // number of slats per quadrant
     const Int_t nPCB5[nSlats5] = {7,7,6,6,5,4,2}; // n PCB per slat

     // slat dimensions: slat is a MOTHER volume!!! made of air
     Float_t slatLength5[nSlats5]; 
     char volNam9[5];
     char volNam10[5];
     Float_t xSlat5;
     Float_t ySlat5;

     for (i = 0; i<nSlats5; i++){
       slatLength5[i] = pcbLength * nPCB5[i] + 2. * dSlatLength; 
       xSlat5 = slatLength5[i]/2. - vFrameLength/2.; 
       if (i==0) xSlat5 += 37.5;
       ySlat5 = sensHeight * i - yOverlap * i; 
       spar[0] = slatLength5[i]/2.; 
       spar[1] = slatHeight/2.;
       spar[2] = slatWidth/2.; 
       // zSlat to be checked (odd downstream or upstream?)
       Float_t zSlat = (i%2 ==0)? -slatWidth/2. : slatWidth/2.; 
       zSlat*=1.01;
       sprintf(volNam9,"S09%d",i);
       gMC->Gsvolu(volNam9,"BOX",slatMaterial,spar,3);
       gMC->Gspos(volNam9, i*4+1,"C09M", xSlat5, ySlat5, zSlat, 0, "ONLY");
       gMC->Gspos(volNam9, i*4+2,"C09M",-xSlat5, ySlat5, zSlat, 0, "ONLY");
       if (i>0) { 
	 gMC->Gspos(volNam9, i*4+3,"C09M", xSlat5,-ySlat5, zSlat, 0, "ONLY");
	 gMC->Gspos(volNam9, i*4+4,"C09M",-xSlat5,-ySlat5, zSlat, 0, "ONLY");
       }
       sprintf(volNam10,"S10%d",i);
       gMC->Gsvolu(volNam10,"BOX",slatMaterial,spar,3);
       gMC->Gspos(volNam10, i*4+1,"C10M", xSlat5, ySlat5, zSlat, 0, "ONLY");
       gMC->Gspos(volNam10, i*4+2,"C10M",-xSlat5, ySlat5, zSlat, 0, "ONLY");
       if (i>0) { 
	 gMC->Gspos(volNam10, i*4+3,"C10M", xSlat5,-ySlat5, zSlat, 0, "ONLY");
	 gMC->Gspos(volNam10, i*4+4,"C10M",-xSlat5,-ySlat5, zSlat, 0, "ONLY");
       }
     }

     // create the sensitive volumes (subdivided as the PCBs),
      
     gMC->Gsvolu("S09G","BOX",sensMaterial,sensPar,3);
     gMC->Gsvolu("S10G","BOX",sensMaterial,sensPar,3);
     
     // create the PCB volume 
     
     gMC->Gsvolu("S09P","BOX",pcbMaterial,pcbpar,3);
     gMC->Gsvolu("S10P","BOX",pcbMaterial,pcbpar,3);
 
     // create the insulating material volume 
     
     gMC->Gsvolu("S09I","BOX",insuMaterial,insupar,3);
     gMC->Gsvolu("S10I","BOX",insuMaterial,insupar,3);

     // create the panel volume 
  
     gMC->Gsvolu("S09C","BOX",panelMaterial,panelpar,3);
     gMC->Gsvolu("S10C","BOX",panelMaterial,panelpar,3);
     
     // create the rohacell volume 
      
     gMC->Gsvolu("S09R","BOX",rohaMaterial,rohapar,3);
     gMC->Gsvolu("S10R","BOX",rohaMaterial,rohapar,3);
     
     // create the vertical frame volume 
   
     gMC->Gsvolu("S09V","BOX",vFrameMaterial,vFramepar,3);
     gMC->Gsvolu("S10V","BOX",vFrameMaterial,vFramepar,3);

     // create the horizontal frame volume 
  
     gMC->Gsvolu("S09H","BOX",hFrameMaterial,hFramepar,3);
     gMC->Gsvolu("S10H","BOX",hFrameMaterial,hFramepar,3);

     // create the horizontal border volume 

     gMC->Gsvolu("S09B","BOX",bFrameMaterial,bFramepar,3);
     gMC->Gsvolu("S10B","BOX",bFrameMaterial,bFramepar,3);

     
     for (i = 0; i<nSlats5; i++){
       sprintf(volNam9,"S09%d",i);
       sprintf(volNam10,"S10%d",i);
       Float_t xvFrame  = (slatLength5[i] - vFrameLength)/2.;
       gMC->Gspos("S09V",2*i-1,volNam9, xvFrame, 0., 0. , 0, "ONLY");
       gMC->Gspos("S09V",2*i  ,volNam9,-xvFrame, 0., 0. , 0, "ONLY");
       gMC->Gspos("S10V",2*i-1,volNam10, xvFrame, 0., 0. , 0, "ONLY");
       gMC->Gspos("S10V",2*i  ,volNam10,-xvFrame, 0., 0. , 0, "ONLY");
       for (j=0; j<nPCB5[i]; j++){
	 index++;
	 Float_t xx = sensLength/2. * (-nPCB5[i]+2*j+1); 
	 Float_t yy = 0.;
	 Float_t zSens = 0.;  
	 gMC->Gspos("S09G",index,volNam9, xx, yy, zSens , 0, "ONLY");
	 gMC->Gspos("S10G",index,volNam10, xx, yy, zSens , 0, "ONLY");
	 Float_t zPCB = (sensWidth+pcbWidth)/2.; 
	 gMC->Gspos("S09P",2*index-1,volNam9, xx, yy, zPCB , 0, "ONLY");
	 gMC->Gspos("S09P",2*index  ,volNam9, xx, yy,-zPCB , 0, "ONLY");
	 gMC->Gspos("S10P",2*index-1,volNam10, xx, yy, zPCB , 0, "ONLY");
	 gMC->Gspos("S10P",2*index  ,volNam10, xx, yy,-zPCB , 0, "ONLY");
	 Float_t zInsu = (insuWidth+pcbWidth)/2. + zPCB; 
	 gMC->Gspos("S09I",2*index-1,volNam9, xx, yy, zInsu , 0, "ONLY");
	 gMC->Gspos("S09I",2*index  ,volNam9, xx, yy,-zInsu , 0, "ONLY");
	 gMC->Gspos("S10I",2*index-1,volNam10, xx, yy, zInsu , 0, "ONLY");
	 gMC->Gspos("S10I",2*index  ,volNam10, xx, yy,-zInsu , 0, "ONLY");
	 Float_t zPanel1 = (insuWidth+panelWidth)/2. + zInsu; 
	 gMC->Gspos("S09C",4*index-3,volNam9, xx, yy, zPanel1 , 0, "ONLY");
	 gMC->Gspos("S09C",4*index-2,volNam9, xx, yy,-zPanel1 , 0, "ONLY");
	 gMC->Gspos("S10C",4*index-3,volNam10, xx, yy, zPanel1 , 0, "ONLY");
	 gMC->Gspos("S10C",4*index-2,volNam10, xx, yy,-zPanel1 , 0, "ONLY");
	 Float_t zRoha = (rohaWidth+panelWidth)/2. + zPanel1; 
	 gMC->Gspos("S09R",2*index-1,volNam9, xx, yy, zRoha , 0, "ONLY");
	 gMC->Gspos("S09R",2*index  ,volNam9, xx, yy,-zRoha , 0, "ONLY");
	 gMC->Gspos("S10R",2*index-1,volNam10, xx, yy, zRoha , 0, "ONLY");
	 gMC->Gspos("S10R",2*index  ,volNam10, xx, yy,-zRoha , 0, "ONLY");
	 Float_t zPanel2 = (rohaWidth+panelWidth)/2. + zRoha; 
	 gMC->Gspos("S09C",4*index-1,volNam9, xx, yy, zPanel2 , 0, "ONLY");
	 gMC->Gspos("S09C",4*index  ,volNam9, xx, yy,-zPanel2 , 0, "ONLY");
	 gMC->Gspos("S10C",4*index-1,volNam10, xx, yy, zPanel2 , 0, "ONLY");
	 gMC->Gspos("S10C",4*index  ,volNam10, xx, yy,-zPanel2 , 0, "ONLY");
	 Float_t yframe = (sensHeight + hFrameHeight)/2.;
	 gMC->Gspos("S09H",2*index-1,volNam9, xx, yframe, 0. , 0, "ONLY");
	 gMC->Gspos("S09H",2*index  ,volNam9, xx,-yframe, 0. , 0, "ONLY");
	 gMC->Gspos("S10H",2*index-1,volNam10, xx, yframe, 0. , 0, "ONLY");
	 gMC->Gspos("S10H",2*index  ,volNam10, xx,-yframe, 0. , 0, "ONLY");
	 Float_t yborder = (bFrameHeight + hFrameHeight)/2. + yframe;
	 gMC->Gspos("S09B",2*index-1,volNam9, xx, yborder, 0. , 0, "ONLY");
	 gMC->Gspos("S09B",2*index  ,volNam9, xx,-yborder, 0. , 0, "ONLY");
	 gMC->Gspos("S10B",2*index-1,volNam10, xx, yborder, 0. , 0, "ONLY");
	 gMC->Gspos("S10B",2*index  ,volNam10, xx,-yborder, 0. , 0, "ONLY");
       } 
     }

     // create the NULOC volume and position it in the horizontal frame
     
     gMC->Gsvolu("S09N","BOX",nulocMaterial,nulocpar,3);
     gMC->Gsvolu("S10N","BOX",nulocMaterial,nulocpar,3);
 
     index = 0; 
     for (xx = -xxmax; xx<=xxmax; xx+=3*nulocLength) { 
       index++; 
       gMC->Gspos("S09N",2*index-1,"S09B", xx, 0.,-bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S09N",2*index  ,"S09B", xx, 0., bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S10N",2*index-1,"S10B", xx, 0.,-bFrameWidth/4., 0, "ONLY");
       gMC->Gspos("S10N",2*index  ,"S10B", xx, 0., bFrameWidth/4., 0, "ONLY");
     }

     // create the gassiplex volume 
     
     gMC->Gsvolu("S09E","BOX",gassiMaterial,gassipar,3);
     gMC->Gsvolu("S10E","BOX",gassiMaterial,gassipar,3);


     // position 4 gassiplex in the nuloc
     
     gMC->Gspos("S09E",1,"S09N", 0., -3 * nulocHeight/8., 0. , 0, "ONLY");
     gMC->Gspos("S09E",2,"S09N", 0.,    - nulocHeight/8., 0. , 0, "ONLY");
     gMC->Gspos("S09E",3,"S09N", 0.,      nulocHeight/8., 0. , 0, "ONLY");
     gMC->Gspos("S09E",4,"S09N", 0.,  3 * nulocHeight/8., 0. , 0, "ONLY");
     gMC->Gspos("S10E",1,"S10N", 0., -3 * nulocHeight/8., 0. , 0, "ONLY");
     gMC->Gspos("S10E",2,"S10N", 0.,    - nulocHeight/8., 0. , 0, "ONLY");
     gMC->Gspos("S10E",3,"S10N", 0.,      nulocHeight/8., 0. , 0, "ONLY");
     gMC->Gspos("S10E",4,"S10N", 0.,  3 * nulocHeight/8., 0. , 0, "ONLY");
     
 }
 

///////////////////////////////////////
// GEOMETRY FOR THE TRIGGER CHAMBERS //
///////////////////////////////////////

// 03/00 P. Dupieux : introduce a slighly more realistic  
//                    geom. of the trigger readout planes with
//                    2 Zpos per trigger plane (alternate
//                    between left and right of the trigger)  

//  Parameters of the Trigger Chambers

    		
     const Float_t kXMC1MIN=34.;       
     const Float_t kXMC1MED=51.;                                
     const Float_t kXMC1MAX=272.;                               
     const Float_t kYMC1MIN=34.;                              
     const Float_t kYMC1MAX=51.;                              
     const Float_t kRMIN1=50.;
     const Float_t kRMAX1=62.;
     const Float_t kRMIN2=50.;
     const Float_t kRMAX2=66.;

//   zposition of the middle of the gas gap in mother vol 
     const Float_t kZMCm=-3.6;
     const Float_t kZMCp=+3.6;


// TRIGGER STATION 1 - TRIGGER STATION 1 - TRIGGER STATION 1

     // iChamber 1 and 2 for first and second chambers in the station
     // iChamber (first chamber) kept for other quanties than Z,
     // assumed to be the same in both chambers
     iChamber1 = iChamber = (AliMUONChamber*) (*fChambers)[10];
     iChamber2 =(AliMUONChamber*) (*fChambers)[11]; 

     // 03/00 
     // zpos1 and zpos2 are now the middle of the first and second
     // plane of station 1 : 
     // zpos1=(16075+15995)/2=16035 mm, thick/2=40 mm
     // zpos2=(16225+16145)/2=16185 mm, thick/2=40 mm
     //
     // zpos1m=15999 mm , zpos1p=16071 mm (middles of gas gaps)
     // zpos2m=16149 mm , zpos2p=16221 mm (middles of gas gaps)
     // rem : the total thickness accounts for 1 mm of al on both 
     // side of the RPCs (see zpos1 and zpos2), as previously

     zpos1=iChamber1->Z();
     zpos2=iChamber2->Z();


// Mother volume definition     
     tpar[0] = iChamber->RInner(); 
     tpar[1] = iChamber->ROuter();
     tpar[2] = 4.0;    
     gMC->Gsvolu("CM11", "TUBE", idAir, tpar, 3);
     gMC->Gsvolu("CM12", "TUBE", idAir, tpar, 3);
     
// Definition of the flange between the beam shielding and the RPC 
     tpar[0]= kRMIN1;
     tpar[1]= kRMAX1;
     tpar[2]= 4.0;
   
     gMC->Gsvolu("CF1A", "TUBE", idAlu1, tpar, 3);     //Al
     gMC->Gspos("CF1A", 1, "CM11", 0., 0., 0., 0, "MANY");
     gMC->Gspos("CF1A", 2, "CM12", 0., 0., 0., 0, "MANY");


// FIRST PLANE OF STATION 1

//   ratios of zpos1m/zpos1p and inverse for first plane
     Float_t zmp=(zpos1-3.6)/(zpos1+3.6);
     Float_t zpm=1./zmp;
   

// Definition of prototype for chambers in the first plane     
          
     tpar[0]= 0.;
     tpar[1]= 0.;
     tpar[2]= 0.;
          
     gMC->Gsvolu("CC1A", "BOX ", idAlu1, tpar, 0);           //Al    
     gMC->Gsvolu("CB1A", "BOX ", idtmed[1107], tpar, 0);     //Bakelite 
     gMC->Gsvolu("CG1A", "BOX ", idtmed[1106], tpar, 0);     //Gas streamer

// chamber type A
     tpar[0] = -1.;
     tpar[1] = -1.;
     
     const Float_t kXMC1A=kXMC1MED+(kXMC1MAX-kXMC1MED)/2.;
     const Float_t kYMC1Am=0.;
     const Float_t kYMC1Ap=0.;
          
     tpar[2] = 0.1;    
     gMC->Gsposp("CG1A", 1, "CB1A", 0., 0., 0., 0, "ONLY",tpar,3);
     tpar[2] = 0.3;
     gMC->Gsposp("CB1A", 1, "CC1A", 0., 0., 0., 0, "ONLY",tpar,3);

     tpar[2] = 0.4;
     tpar[0] = (kXMC1MAX-kXMC1MED)/2.;
     tpar[1] = kYMC1MIN;

     gMC->Gsposp("CC1A", 1, "CM11",kXMC1A,kYMC1Am,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 2, "CM11",-kXMC1A,kYMC1Ap,kZMCp, 0, "ONLY", tpar, 3);
     
//  chamber type B    
     Float_t tpar1save=tpar[1];
     Float_t y1msave=kYMC1Am;
     Float_t y1psave=kYMC1Ap;
 
     tpar[0] = (kXMC1MAX-kXMC1MIN)/2.;
     tpar[1] = (kYMC1MAX-kYMC1MIN)/2.;
     
     const Float_t kXMC1B=kXMC1MIN+tpar[0];
     const Float_t kYMC1Bp=(y1msave+tpar1save)*zpm+tpar[1];
     const Float_t kYMC1Bm=(y1psave+tpar1save)*zmp+tpar[1];

     gMC->Gsposp("CC1A", 3, "CM11",kXMC1B,kYMC1Bp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 4, "CM11",-kXMC1B,kYMC1Bm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 5, "CM11",kXMC1B,-kYMC1Bp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 6, "CM11",-kXMC1B,-kYMC1Bm,kZMCm, 0, "ONLY", tpar, 3);
     
//  chamber type C  (end of type B !!)      
     tpar1save=tpar[1];
     y1msave=kYMC1Bm;
     y1psave=kYMC1Bp;

     tpar[0] = kXMC1MAX/2;
     tpar[1] = kYMC1MAX/2;
     
     const Float_t kXMC1C=tpar[0];
// warning : same Z than type B
     const Float_t kYMC1Cp=(y1psave+tpar1save)*1.+tpar[1];
     const Float_t kYMC1Cm=(y1msave+tpar1save)*1.+tpar[1];
     
     gMC->Gsposp("CC1A", 7, "CM11",kXMC1C,kYMC1Cp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 8, "CM11",-kXMC1C,kYMC1Cm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 9, "CM11",kXMC1C,-kYMC1Cp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 10, "CM11",-kXMC1C,-kYMC1Cm,kZMCm, 0, "ONLY", tpar, 3);
     
//  chamber type D, E and F (same size)        
     tpar1save=tpar[1];
     y1msave=kYMC1Cm;
     y1psave=kYMC1Cp;

     tpar[0] = kXMC1MAX/2.;
     tpar[1] = kYMC1MIN;
     
     const Float_t kXMC1D=tpar[0];
     const Float_t kYMC1Dp=(y1msave+tpar1save)*zpm+tpar[1];
     const Float_t kYMC1Dm=(y1psave+tpar1save)*zmp+tpar[1];
     
     gMC->Gsposp("CC1A", 11, "CM11",kXMC1D,kYMC1Dm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 12, "CM11",-kXMC1D,kYMC1Dp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 13, "CM11",kXMC1D,-kYMC1Dm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 14, "CM11",-kXMC1D,-kYMC1Dp,kZMCp, 0, "ONLY", tpar, 3);


     tpar1save=tpar[1];
     y1msave=kYMC1Dm;
     y1psave=kYMC1Dp;
     const Float_t kYMC1Ep=(y1msave+tpar1save)*zpm+tpar[1];
     const Float_t kYMC1Em=(y1psave+tpar1save)*zmp+tpar[1];
     
     gMC->Gsposp("CC1A", 15, "CM11",kXMC1D,kYMC1Ep,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 16, "CM11",-kXMC1D,kYMC1Em,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 17, "CM11",kXMC1D,-kYMC1Ep,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 18, "CM11",-kXMC1D,-kYMC1Em,kZMCm, 0, "ONLY", tpar, 3);

     tpar1save=tpar[1];
     y1msave=kYMC1Em;
     y1psave=kYMC1Ep;
     const Float_t kYMC1Fp=(y1msave+tpar1save)*zpm+tpar[1];
     const Float_t kYMC1Fm=(y1psave+tpar1save)*zmp+tpar[1];
    
     gMC->Gsposp("CC1A", 19, "CM11",kXMC1D,kYMC1Fm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 20, "CM11",-kXMC1D,kYMC1Fp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 21, "CM11",kXMC1D,-kYMC1Fm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC1A", 22, "CM11",-kXMC1D,-kYMC1Fp,kZMCp, 0, "ONLY", tpar, 3);

// Positioning first plane in ALICE     
     gMC->Gspos("CM11", 1, "ALIC", 0., 0., zpos1, 0, "ONLY");

// End of geometry definition for the first plane of station 1



// SECOND PLANE OF STATION 1 : proj ratio = zpos2/zpos1

     const Float_t kZ12=zpos2/zpos1;
      
// Definition of prototype for chambers in the second plane of station 1    
          
     tpar[0]= 0.;
     tpar[1]= 0.;
     tpar[2]= 0.;
          
     gMC->Gsvolu("CC2A", "BOX ", idAlu1, tpar, 0);           //Al    
     gMC->Gsvolu("CB2A", "BOX ", idtmed[1107], tpar, 0);     //Bakelite 
     gMC->Gsvolu("CG2A", "BOX ", idtmed[1106], tpar, 0);     //Gas streamer

// chamber type A
     tpar[0] = -1.;
     tpar[1] = -1.;
     
     const Float_t kXMC2A=kXMC1A*kZ12;
     const Float_t kYMC2Am=0.;
     const Float_t kYMC2Ap=0.;
          
     tpar[2] = 0.1;    
     gMC->Gsposp("CG2A", 1, "CB2A", 0., 0., 0., 0, "ONLY",tpar,3);
     tpar[2] = 0.3;
     gMC->Gsposp("CB2A", 1, "CC2A", 0., 0., 0., 0, "ONLY",tpar,3);

     tpar[2] = 0.4;
     tpar[0] = ((kXMC1MAX-kXMC1MED)/2.)*kZ12;
     tpar[1] = kYMC1MIN*kZ12;

     gMC->Gsposp("CC2A", 1, "CM12",kXMC2A,kYMC2Am,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 2, "CM12",-kXMC2A,kYMC2Ap,kZMCp, 0, "ONLY", tpar, 3);
     

//  chamber type B    

     tpar[0] = ((kXMC1MAX-kXMC1MIN)/2.)*kZ12;
     tpar[1] = ((kYMC1MAX-kYMC1MIN)/2.)*kZ12;
     
     const Float_t kXMC2B=kXMC1B*kZ12;
     const Float_t kYMC2Bp=kYMC1Bp*kZ12;
     const Float_t kYMC2Bm=kYMC1Bm*kZ12;
     gMC->Gsposp("CC2A", 3, "CM12",kXMC2B,kYMC2Bp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 4, "CM12",-kXMC2B,kYMC2Bm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 5, "CM12",kXMC2B,-kYMC2Bp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 6, "CM12",-kXMC2B,-kYMC2Bm,kZMCm, 0, "ONLY", tpar, 3);

     
//  chamber type C   (end of type B !!)     

     tpar[0] = (kXMC1MAX/2)*kZ12;
     tpar[1] = (kYMC1MAX/2)*kZ12;
     
     const Float_t kXMC2C=kXMC1C*kZ12;
     const Float_t kYMC2Cp=kYMC1Cp*kZ12;
     const Float_t kYMC2Cm=kYMC1Cm*kZ12;     
     gMC->Gsposp("CC2A", 7, "CM12",kXMC2C,kYMC2Cp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 8, "CM12",-kXMC2C,kYMC2Cm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 9, "CM12",kXMC2C,-kYMC2Cp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 10, "CM12",-kXMC2C,-kYMC2Cm,kZMCm, 0, "ONLY", tpar, 3);
     
//  chamber type D, E and F (same size)        

     tpar[0] = (kXMC1MAX/2.)*kZ12;
     tpar[1] = kYMC1MIN*kZ12;
     
     const Float_t kXMC2D=kXMC1D*kZ12;
     const Float_t kYMC2Dp=kYMC1Dp*kZ12;
     const Float_t kYMC2Dm=kYMC1Dm*kZ12;     
     gMC->Gsposp("CC2A", 11, "CM12",kXMC2D,kYMC2Dm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 12, "CM12",-kXMC2D,kYMC2Dp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 13, "CM12",kXMC2D,-kYMC2Dm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 14, "CM12",-kXMC2D,-kYMC2Dp,kZMCp, 0, "ONLY", tpar, 3);

     const Float_t kYMC2Ep=kYMC1Ep*kZ12;
     const Float_t kYMC2Em=kYMC1Em*kZ12;
     gMC->Gsposp("CC2A", 15, "CM12",kXMC2D,kYMC2Ep,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 16, "CM12",-kXMC2D,kYMC2Em,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 17, "CM12",kXMC2D,-kYMC2Ep,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 18, "CM12",-kXMC2D,-kYMC2Em,kZMCm, 0, "ONLY", tpar, 3);


     const Float_t kYMC2Fp=kYMC1Fp*kZ12;
     const Float_t kYMC2Fm=kYMC1Fm*kZ12;
     gMC->Gsposp("CC2A", 19, "CM12",kXMC2D,kYMC2Fm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 20, "CM12",-kXMC2D,kYMC2Fp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 21, "CM12",kXMC2D,-kYMC2Fm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC2A", 22, "CM12",-kXMC2D,-kYMC2Fp,kZMCp, 0, "ONLY", tpar, 3);

// Positioning second plane of station 1 in ALICE     
     
     gMC->Gspos("CM12", 1, "ALIC", 0., 0., zpos2, 0, "ONLY");

// End of geometry definition for the second plane of station 1



// TRIGGER STATION 2 - TRIGGER STATION 2 - TRIGGER STATION 2    

     // 03/00 
     // zpos3 and zpos4 are now the middle of the first and second
     // plane of station 2 : 
     // zpos3=(17075+16995)/2=17035 mm, thick/2=40 mm
     // zpos4=(17225+17145)/2=17185 mm, thick/2=40 mm
     //
     // zpos3m=16999 mm , zpos3p=17071 mm (middles of gas gaps)
     // zpos4m=17149 mm , zpos4p=17221 mm (middles of gas gaps)
     // rem : the total thickness accounts for 1 mm of al on both 
     // side of the RPCs (see zpos3 and zpos4), as previously
     iChamber1 = iChamber = (AliMUONChamber*) (*fChambers)[12];
     iChamber2 =(AliMUONChamber*) (*fChambers)[13];
     Float_t zpos3=iChamber1->Z();
     Float_t zpos4=iChamber2->Z();


// Mother volume definition     
     tpar[0] = iChamber->RInner(); 
     tpar[1] = iChamber->ROuter();
     tpar[2] = 4.0;    
 
     gMC->Gsvolu("CM21", "TUBE", idAir, tpar, 3);
     gMC->Gsvolu("CM22", "TUBE", idAir, tpar, 3);
     
// Definition of the flange between the beam shielding and the RPC 
//  ???? interface shielding

     tpar[0]= kRMIN2;
     tpar[1]= kRMAX2;
     tpar[2]= 4.0;
   
     gMC->Gsvolu("CF2A", "TUBE", idAlu1, tpar, 3);            //Al
     gMC->Gspos("CF2A", 1, "CM21", 0., 0., 0., 0, "MANY");
     gMC->Gspos("CF2A", 2, "CM22", 0., 0., 0., 0, "MANY");
    


// FIRST PLANE OF STATION 2 : proj ratio = zpos3/zpos1

     const Float_t kZ13=zpos3/zpos1; 

// Definition of prototype for chambers in the first plane of station 2       
     tpar[0]= 0.;
     tpar[1]= 0.;
     tpar[2]= 0.;
          
     gMC->Gsvolu("CC3A", "BOX ", idAlu1, tpar, 0);           //Al  
     gMC->Gsvolu("CB3A", "BOX ", idtmed[1107], tpar, 0);     //Bakelite 
     gMC->Gsvolu("CG3A", "BOX ", idtmed[1106], tpar, 0);     //Gas streamer


// chamber type A
     tpar[0] = -1.;
     tpar[1] = -1.;
     
     const Float_t kXMC3A=kXMC1A*kZ13;
     const Float_t kYMC3Am=0.;
     const Float_t kYMC3Ap=0.;
          
     tpar[2] = 0.1;    
     gMC->Gsposp("CG3A", 1, "CB3A", 0., 0., 0., 0, "ONLY",tpar,3);
     tpar[2] = 0.3;
     gMC->Gsposp("CB3A", 1, "CC3A", 0., 0., 0., 0, "ONLY",tpar,3);

     tpar[2] = 0.4;
     tpar[0] = ((kXMC1MAX-kXMC1MED)/2.)*kZ13;
     tpar[1] = kYMC1MIN*kZ13;
     gMC->Gsposp("CC3A", 1, "CM21",kXMC3A,kYMC3Am,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 2, "CM21",-kXMC3A,kYMC3Ap,kZMCp, 0, "ONLY", tpar, 3);

     
//  chamber type B    
     tpar[0] = ((kXMC1MAX-kXMC1MIN)/2.)*kZ13;
     tpar[1] = ((kYMC1MAX-kYMC1MIN)/2.)*kZ13;
     
     const Float_t kXMC3B=kXMC1B*kZ13;
     const Float_t kYMC3Bp=kYMC1Bp*kZ13;
     const Float_t kYMC3Bm=kYMC1Bm*kZ13;
     gMC->Gsposp("CC3A", 3, "CM21",kXMC3B,kYMC3Bp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 4, "CM21",-kXMC3B,kYMC3Bm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 5, "CM21",kXMC3B,-kYMC3Bp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 6, "CM21",-kXMC3B,-kYMC3Bm,kZMCm, 0, "ONLY", tpar, 3);

     
//  chamber type C  (end of type B !!)      
     tpar[0] = (kXMC1MAX/2)*kZ13;
     tpar[1] = (kYMC1MAX/2)*kZ13;
     
     const Float_t kXMC3C=kXMC1C*kZ13;
     const Float_t kYMC3Cp=kYMC1Cp*kZ13;
     const Float_t kYMC3Cm=kYMC1Cm*kZ13;     
     gMC->Gsposp("CC3A", 7, "CM21",kXMC3C,kYMC3Cp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 8, "CM21",-kXMC3C,kYMC3Cm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 9, "CM21",kXMC3C,-kYMC3Cp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 10, "CM21",-kXMC3C,-kYMC3Cm,kZMCm, 0, "ONLY", tpar, 3);
     

//  chamber type D, E and F (same size)         

     tpar[0] = (kXMC1MAX/2.)*kZ13;
     tpar[1] = kYMC1MIN*kZ13;
     
     const Float_t kXMC3D=kXMC1D*kZ13;
     const Float_t kYMC3Dp=kYMC1Dp*kZ13;
     const Float_t kYMC3Dm=kYMC1Dm*kZ13;          
     gMC->Gsposp("CC3A", 11, "CM21",kXMC3D,kYMC3Dm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 12, "CM21",-kXMC3D,kYMC3Dp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 13, "CM21",kXMC3D,-kYMC3Dm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 14, "CM21",-kXMC3D,-kYMC3Dp,kZMCp, 0, "ONLY", tpar, 3);

     const Float_t kYMC3Ep=kYMC1Ep*kZ13;
     const Float_t kYMC3Em=kYMC1Em*kZ13;
     gMC->Gsposp("CC3A", 15, "CM21",kXMC3D,kYMC3Ep,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 16, "CM21",-kXMC3D,kYMC3Em,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 17, "CM21",kXMC3D,-kYMC3Ep,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 18, "CM21",-kXMC3D,-kYMC3Em,kZMCm, 0, "ONLY", tpar, 3);

     const Float_t kYMC3Fp=kYMC1Fp*kZ13;
     const Float_t kYMC3Fm=kYMC1Fm*kZ13;
     gMC->Gsposp("CC3A", 19, "CM21",kXMC3D,kYMC3Fm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 20, "CM21",-kXMC3D,kYMC3Fp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 21, "CM21",kXMC3D,-kYMC3Fm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC3A", 22, "CM21",-kXMC3D,-kYMC3Fp,kZMCp, 0, "ONLY", tpar, 3);
       

// Positioning first plane of station 2 in ALICE
     
     gMC->Gspos("CM21", 1, "ALIC", 0., 0., zpos3, 0, "ONLY");

// End of geometry definition for the first plane of station 2




// SECOND PLANE OF STATION 2 : proj ratio = zpos4/zpos1

     const Float_t kZ14=zpos4/zpos1;
     
// Definition of prototype for chambers in the second plane of station 2    
          
     tpar[0]= 0.;
     tpar[1]= 0.;
     tpar[2]= 0.;
          
     gMC->Gsvolu("CC4A", "BOX ", idAlu1, tpar, 0);           //Al      
     gMC->Gsvolu("CB4A", "BOX ", idtmed[1107], tpar, 0);     //Bakelite 
     gMC->Gsvolu("CG4A", "BOX ", idtmed[1106], tpar, 0);     //Gas streamer

// chamber type A
     tpar[0] = -1.;
     tpar[1] = -1.;
     
     const Float_t kXMC4A=kXMC1A*kZ14;
     const Float_t kYMC4Am=0.;
     const Float_t kYMC4Ap=0.;
          
     tpar[2] = 0.1;    
     gMC->Gsposp("CG4A", 1, "CB4A", 0., 0., 0., 0, "ONLY",tpar,3);
     tpar[2] = 0.3;
     gMC->Gsposp("CB4A", 1, "CC4A", 0., 0., 0., 0, "ONLY",tpar,3);

     tpar[2] = 0.4;
     tpar[0] = ((kXMC1MAX-kXMC1MED)/2.)*kZ14;
     tpar[1] = kYMC1MIN*kZ14;
     gMC->Gsposp("CC4A", 1, "CM22",kXMC4A,kYMC4Am,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 2, "CM22",-kXMC4A,kYMC4Ap,kZMCp, 0, "ONLY", tpar, 3);
     

//  chamber type B    
     tpar[0] = ((kXMC1MAX-kXMC1MIN)/2.)*kZ14;
     tpar[1] = ((kYMC1MAX-kYMC1MIN)/2.)*kZ14;
     
     const Float_t kXMC4B=kXMC1B*kZ14;
     const Float_t kYMC4Bp=kYMC1Bp*kZ14;
     const Float_t kYMC4Bm=kYMC1Bm*kZ14;
     gMC->Gsposp("CC4A", 3, "CM22",kXMC4B,kYMC4Bp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 4, "CM22",-kXMC4B,kYMC4Bm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 5, "CM22",kXMC4B,-kYMC4Bp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 6, "CM22",-kXMC4B,-kYMC4Bm,kZMCm, 0, "ONLY", tpar, 3);

     
//  chamber type C   (end of type B !!)      
     tpar[0] =(kXMC1MAX/2)*kZ14;
     tpar[1] =  (kYMC1MAX/2)*kZ14;
     
     const Float_t kXMC4C=kXMC1C*kZ14;
     const Float_t kYMC4Cp=kYMC1Cp*kZ14;
     const Float_t kYMC4Cm=kYMC1Cm*kZ14;     
     gMC->Gsposp("CC4A", 7, "CM22",kXMC4C,kYMC4Cp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 8, "CM22",-kXMC4C,kYMC4Cm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 9, "CM22",kXMC4C,-kYMC4Cp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 10, "CM22",-kXMC4C,-kYMC4Cm,kZMCm, 0, "ONLY", tpar, 3);

     
//  chamber type D, E and F (same size)      
     tpar[0] = (kXMC1MAX/2.)*kZ14;
     tpar[1] =  kYMC1MIN*kZ14;
     
     const Float_t kXMC4D=kXMC1D*kZ14;
     const Float_t kYMC4Dp=kYMC1Dp*kZ14;
     const Float_t kYMC4Dm=kYMC1Dm*kZ14;          
     gMC->Gsposp("CC4A", 11, "CM22",kXMC4D,kYMC4Dm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 12, "CM22",-kXMC4D,kYMC4Dp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 13, "CM22",kXMC4D,-kYMC4Dm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 14, "CM22",-kXMC4D,-kYMC4Dp,kZMCp, 0, "ONLY", tpar, 3);

     const Float_t kYMC4Ep=kYMC1Ep*kZ14;
     const Float_t kYMC4Em=kYMC1Em*kZ14;          
     gMC->Gsposp("CC4A", 15, "CM22",kXMC4D,kYMC4Ep,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 16, "CM22",-kXMC4D,kYMC4Em,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 17, "CM22",kXMC4D,-kYMC4Ep,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 18, "CM22",-kXMC4D,-kYMC4Em,kZMCm, 0, "ONLY", tpar, 3);

     const Float_t kYMC4Fp=kYMC1Fp*kZ14;
     const Float_t kYMC4Fm=kYMC1Fm*kZ14;          
     gMC->Gsposp("CC4A", 19, "CM22",kXMC4D,kYMC4Fm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 20, "CM22",-kXMC4D,kYMC4Fp,kZMCp, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 21, "CM22",kXMC4D,-kYMC4Fm,kZMCm, 0, "ONLY", tpar, 3);
     gMC->Gsposp("CC4A", 22, "CM22",-kXMC4D,-kYMC4Fp,kZMCp, 0, "ONLY", tpar, 3);
     

// Positioning second plane of station 2 in ALICE
     
     gMC->Gspos("CM22", 1, "ALIC", 0., 0., zpos4, 0, "ONLY");

// End of geometry definition for the second plane of station 2

// End of trigger geometry definition

}


 
//___________________________________________
void AliMUONv1::CreateMaterials()
{
  // *** DEFINITION OF AVAILABLE MUON MATERIALS *** 
  //
  //     Ar-CO2 gas 
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

void AliMUONv1::Init()
{
   // 
   // Initialize Tracking Chambers
   //

   printf("\n\n\n Start Init for version 1 - CPC chamber type\n\n\n");
   Int_t i;
   for (i=0; i<AliMUONConstants::NCh(); i++) {
       ( (AliMUONChamber*) (*fChambers)[i])->Init();
   }
   
   //
   // Set the chamber (sensitive region) GEANT identifier
   AliMC* gMC = AliMC::GetMC(); 
   ((AliMUONChamber*)(*fChambers)[0])->SetGid(gMC->VolId("C01G"));
   ((AliMUONChamber*)(*fChambers)[1])->SetGid(gMC->VolId("C02G"));

   ((AliMUONChamber*)(*fChambers)[2])->SetGid(gMC->VolId("C03G"));
   ((AliMUONChamber*)(*fChambers)[3])->SetGid(gMC->VolId("C04G"));

   ((AliMUONChamber*)(*fChambers)[4])->SetGid(gMC->VolId("S05G"));
   ((AliMUONChamber*)(*fChambers)[5])->SetGid(gMC->VolId("S06G"));

   ((AliMUONChamber*)(*fChambers)[6])->SetGid(gMC->VolId("S07G"));
   ((AliMUONChamber*)(*fChambers)[7])->SetGid(gMC->VolId("S08G"));

   ((AliMUONChamber*)(*fChambers)[8])->SetGid(gMC->VolId("S09G"));
   ((AliMUONChamber*)(*fChambers)[9])->SetGid(gMC->VolId("S10G"));

   ((AliMUONChamber*)(*fChambers)[10])->SetGid(gMC->VolId("CG1A"));
   ((AliMUONChamber*)(*fChambers)[11])->SetGid(gMC->VolId("CG2A"));
   ((AliMUONChamber*)(*fChambers)[12])->SetGid(gMC->VolId("CG3A"));
   ((AliMUONChamber*)(*fChambers)[13])->SetGid(gMC->VolId("CG4A"));

   printf("\n\n\n Finished Init for version 0 - CPC chamber type\n\n\n");

   //cp 
   printf("\n\n\n Start Init for Trigger Circuits\n\n\n");
   for (i=0; i<AliMUONConstants::NTriggerCircuit(); i++) {
     ( (AliMUONTriggerCircuit*) (*fTriggerCircuits)[i])->Init(i);
   }
   printf(" Finished Init for Trigger Circuits\n\n\n");
   //cp

}

//___________________________________________
void AliMUONv1::StepManager()
{
  Int_t          copy, id;
  static Int_t   idvol;
  static Int_t   vol[2];
  Int_t          ipart;
  TLorentzVector pos;
  TLorentzVector mom;
  Float_t        theta,phi;
  Float_t        destep, step;
  
  static Float_t eloss, eloss2, xhit, yhit, zhit, tof, tlength;
  const  Float_t kBig=1.e10;
  //  modifs perso
  static Float_t hits[15];

  TClonesArray &lhits = *fHits;

  //
  // Set maximum step size for gas
  // numed=gMC->GetMedium();
  //
  // Only charged tracks
  if( !(gMC->TrackCharge()) ) return; 
  //
  // Only gas gap inside chamber
  // Tag chambers and record hits when track enters 
  idvol=-1;
  id=gMC->CurrentVolID(copy);
  
    for (Int_t i=1; i<=AliMUONConstants::NCh(); i++) {
      if(id==((AliMUONChamber*)(*fChambers)[i-1])->GetGid()){ 
	  vol[0]=i; 
	  idvol=i-1;
      }
    }
    if (idvol == -1) return;
  //
  // Get current particle id (ipart), track position (pos)  and momentum (mom) 
  gMC->TrackPosition(pos);
  gMC->TrackMomentum(mom);

  ipart  = gMC->TrackPid();
  //Int_t ipart1 = gMC->IdFromPDG(ipart);
  //printf("ich, ipart %d %d \n",vol[0],ipart1);

  //
  // momentum loss and steplength in last step
  destep = gMC->Edep();
  step   = gMC->TrackStep();
  
  //
  // record hits when track enters ...
  if( gMC->IsTrackEntering()) {
      gMC->SetMaxStep(fMaxStepGas);
      Double_t tc = mom[0]*mom[0]+mom[1]*mom[1];
      Double_t rt = TMath::Sqrt(tc);
      Double_t pmom = TMath::Sqrt(tc+mom[2]*mom[2]);
      Double_t tx=mom[0]/pmom;
      Double_t ty=mom[1]/pmom;
      Double_t tz=mom[2]/pmom;
      Double_t s=((AliMUONChamber*)(*fChambers)[idvol])
	  ->ResponseModel()
	  ->Pitch()/tz;
      theta   = Float_t(TMath::ATan2(rt,Double_t(mom[2])))*kRaddeg;
      phi     = Float_t(TMath::ATan2(Double_t(mom[1]),Double_t(mom[0])))*kRaddeg;
      hits[0] = Float_t(ipart);         // Geant3 particle type
      hits[1] = pos[0]+s*tx;                 // X-position for hit
      hits[2] = pos[1]+s*ty;                 // Y-position for hit
      hits[3] = pos[2]+s*tz;                 // Z-position for hit
      hits[4] = theta;                  // theta angle of incidence
      hits[5] = phi;                    // phi angle of incidence 
      hits[8] = (Float_t) fNPadHits;   // first padhit
      hits[9] = -1;                     // last pad hit

      // modifs perso
      hits[10] = mom[3]; // hit momentum P
      hits[11] = mom[0]; // Px/P
      hits[12] = mom[1]; // Py/P
      hits[13] = mom[2]; // Pz/P
      // fin modifs perso
      tof=gMC->TrackTime();
      hits[14] = tof;    // Time of flight
      // phi angle of incidence
      tlength = 0;
      eloss   = 0;
      eloss2  = 0;
      xhit    = pos[0];
      yhit    = pos[1];      
      zhit    = pos[2];      
      // Only if not trigger chamber

      
      

      if(idvol<AliMUONConstants::NTrackingCh()) {
	  //
	  //  Initialize hit position (cursor) in the segmentation model 
	  ((AliMUONChamber*) (*fChambers)[idvol])
	      ->SigGenInit(pos[0], pos[1], pos[2]);
      } else {
	  //geant3->Gpcxyz();
	  //printf("In the Trigger Chamber #%d\n",idvol-9);
      }
  }
  eloss2+=destep;
  
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

      if(idvol<AliMUONConstants::NTrackingCh()) {
// tracking chambers
	  x0 = 0.5*(xhit+pos[0]);
	  y0 = 0.5*(yhit+pos[1]);
	  z0 = 0.5*(zhit+pos[2]);
	  //	  z0 = localPos[2];
      } else {
// trigger chambers
	  x0=xhit;
	  y0=yhit;
//	  z0=yhit;
	  z0=0.;
      }
      

      if (eloss >0)  MakePadHits(x0,y0,z0,eloss,tof,idvol);
      
	  
      hits[6]=tlength;
      hits[7]=eloss2;
      if (fNPadHits > (Int_t)hits[8]) {
	  hits[8]= hits[8]+1;
	  hits[9]= (Float_t) fNPadHits;
      }
    
      new(lhits[fNhits++]) 
	  AliMUONHit(fIshunt,gAlice->CurrentTrack(),vol,hits);
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


      if (eloss > 0 && idvol < AliMUONConstants::NTrackingCh())
	MakePadHits(0.5*(xhit+pos[0]),0.5*(yhit+pos[1]),pos[2],eloss,tof,idvol);
      xhit     = pos[0];
      yhit     = pos[1]; 
      zhit     = pos[2]; 
      eloss    = destep;
      tlength += step ;
      //
      // nothing special  happened, add up energy loss
  } else {        
      eloss   += destep;
      tlength += step ;
  }
}


