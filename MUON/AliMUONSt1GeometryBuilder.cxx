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

// $Id$
//
// Class AliMUONSt1GeometryBuilder
// -------------------------------
// MUON Station1 coarse geometry construction class.
// Extracted from AliMUONv1
// by Ivana Hrivnacova, IPN Orsay
// Included in AliRoot 2004/01/23

#include <TVirtualMC.h>
#include <TGeoMatrix.h>

#include "AliLog.h"

#include "AliMUONSt1GeometryBuilder.h"
#include "AliMUON.h"
#include "AliMUONConstants.h"
#include "AliMUONGeometryModule.h"
#include "AliMUONGeometryEnvelopeStore.h"

/// \cond CLASSIMP
ClassImp(AliMUONSt1GeometryBuilder)
/// \endcond

//______________________________________________________________________________
AliMUONSt1GeometryBuilder::AliMUONSt1GeometryBuilder(AliMUON* muon)
 : AliMUONVGeometryBuilder(0, 2),
   fMUON(muon)
{
/// Standard constructor

}

//______________________________________________________________________________
AliMUONSt1GeometryBuilder::AliMUONSt1GeometryBuilder()
 : AliMUONVGeometryBuilder(),
   fMUON(0)
{
/// Default constructor
}

//______________________________________________________________________________
AliMUONSt1GeometryBuilder::~AliMUONSt1GeometryBuilder() 
{
/// Destructor
}

//
// public methods
//

//______________________________________________________________________________
void AliMUONSt1GeometryBuilder::CreateGeometry()
{
/// From AliMUONv1::CreateGeometry()

//********************************************************************
//                            Station 1                             **
//********************************************************************
//  CONCENTRIC
     // indices 1 and 2 for first and second chambers in the station
     // iChamber (first chamber) kept for other quanties than Z,
     // assumed to be the same in both chambers

     // Get tracking medias Ids     
     Int_t *idtmed = fMUON->GetIdtmed()->GetArray()-1099;
     Int_t idAir= idtmed[1100]; // medium 1
     Int_t idAlu1=idtmed[1103]; // medium 4
     Int_t idAlu2=idtmed[1104]; // medium 5
     Int_t idGas=idtmed[1108];  // medium 9 = Ar-CO2 gas (80%+20%)
     Bool_t frameCrosses=kTRUE;     

     // Rotation matrices in the x-y plane  
     // phi=   0 deg
     Int_t irot1;
     fMUON->AliMatrix(irot1,  90.,   0., 90.,  90., 0., 0.);
     // phi=  90 deg
     Int_t irot2;
     fMUON->AliMatrix(irot2,  90.,  90., 90., 180., 0., 0.);

     // DGas decreased from standard one (0.5)
     const Float_t kDGas = 0.4;

     // DAlu increased from standard one (3% of X0),
     // because more electronics with smaller pads
     const Float_t kDAlu = 3.5 * 8.9 / 100.;

     // Half of the total thickness of frame crosses (including DAlu)
     // for each chamber in stations 1 and 2:
     // 3% of X0 of composite material,
     // but taken as Aluminium here, with same thickness in number of X0
     Float_t dframez = 3. * 8.9 / 100;
     Float_t zfpos=-(kDGas+dframez+kDAlu)/2;
             // The same parameters are defined in builder for station 2 

     // Mother volume
     // Outer excess and inner recess for mother volume radius
     // with respect to ROuter and RInner
     Float_t dframep=.001; // Value for station 3 should be 6 ...
     // Width (RdPhi) of the frame crosses for stations 1 and 2 (cm)
     // Float_t dframep1=.001;
     Float_t dframep1 = 11.0;
     Float_t phi=2*TMath::Pi()/12/2;
             // The same parameters are defined in builder for station 2 

     Float_t tpar[3];
     Float_t dstation = (-AliMUONConstants::DefaultChamberZ(1)) - 
                        (-AliMUONConstants::DefaultChamberZ(0));
     tpar[0] = AliMUONConstants::Rmin(0)-dframep; 
     tpar[1] = (AliMUONConstants::Rmax(0)+dframep)/TMath::Cos(phi);
     tpar[2] = dstation/5;

     gMC->Gsvolu("S01M", "TUBE", idAir, tpar, 3);
     gMC->Gsvolu("S02M", "TUBE", idAir, tpar, 3);

     // CHANGED
     //gMC->Gspos("S01M", 1, "ALIC", 0., 0., zpos1 , 0, "ONLY");
     //gMC->Gspos("S02M", 1, "ALIC", 0., 0., zpos2 , 0, "ONLY");
     
     GetEnvelopes(0)->AddEnvelope("S01M", 100, false);
     GetEnvelopes(1)->AddEnvelope("S02M", 200, false);
         

// // Aluminium frames
// // Outer frames
//      pgpar[0] = 360/12/2;
//      pgpar[1] = 360.;
//      pgpar[2] = 12.;
//      pgpar[3] =   2;
//      pgpar[4] = -dframez/2;
//      pgpar[5] = AliMUONConstants::Rmax(0);
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
//      tpar[0]= AliMUONConstants::Rmin(0)-dframep1;
//      tpar[1]= AliMUONConstants::Rmin(0);
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
         Float_t bpar[3];
	 bpar[0] = (AliMUONConstants::Rmax(0) - AliMUONConstants::Rmin(0))
	   * TMath::Cos(TMath::ASin(dframep1 /
		   (AliMUONConstants::Rmax(0) - AliMUONConstants::Rmin(0))))
	   / 2.0;
	 bpar[1] = dframep1/2;
	 // total thickness will be (4 * bpar[2]) for each chamber,
	 // which has to be equal to (2 * dframez) - DAlu
	 bpar[2] = (2.0 * dframez - kDAlu) / 4.0;
	 gMC->Gsvolu("S01B", "BOX", idAlu1, bpar, 3);
	 gMC->Gsvolu("S02B", "BOX", idAlu1, bpar, 3);
	 
	 gMC->Gspos("S01B",1,"S01M", +AliMUONConstants::Rmin(0)+bpar[0] , 0,-zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S01B",2,"S01M", -AliMUONConstants::Rmin(0)-bpar[0] , 0,-zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S01B",3,"S01M", 0, +AliMUONConstants::Rmin(0)+bpar[0] ,-zfpos, 
		    irot2,"ONLY");
	 gMC->Gspos("S01B",4,"S01M", 0, -AliMUONConstants::Rmin(0)-bpar[0] ,-zfpos, 
		    irot2,"ONLY");
	 gMC->Gspos("S01B",5,"S01M", +AliMUONConstants::Rmin(0)+bpar[0] , 0,+zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S01B",6,"S01M", -AliMUONConstants::Rmin(0)-bpar[0] , 0,+zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S01B",7,"S01M", 0, +AliMUONConstants::Rmin(0)+bpar[0] ,+zfpos, 
		    irot2,"ONLY");
	 gMC->Gspos("S01B",8,"S01M", 0, -AliMUONConstants::Rmin(0)-bpar[0] ,+zfpos, 
		    irot2,"ONLY");
	 
	 gMC->Gspos("S02B",1,"S02M", +AliMUONConstants::Rmin(0)+bpar[0] , 0,-zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S02B",2,"S02M", -AliMUONConstants::Rmin(0)-bpar[0] , 0,-zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S02B",3,"S02M", 0, +AliMUONConstants::Rmin(0)+bpar[0] ,-zfpos, 
		    irot2,"ONLY");
	 gMC->Gspos("S02B",4,"S02M", 0, -AliMUONConstants::Rmin(0)-bpar[0] ,-zfpos, 
		    irot2,"ONLY");
	 gMC->Gspos("S02B",5,"S02M", +AliMUONConstants::Rmin(0)+bpar[0] , 0,+zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S02B",6,"S02M", -AliMUONConstants::Rmin(0)-bpar[0] , 0,+zfpos, 
		    irot1,"ONLY");
	 gMC->Gspos("S02B",7,"S02M", 0, +AliMUONConstants::Rmin(0)+bpar[0] ,+zfpos, 
		    irot2,"ONLY");
	 gMC->Gspos("S02B",8,"S02M", 0, -AliMUONConstants::Rmin(0)-bpar[0] ,+zfpos, 
		    irot2,"ONLY");
     }
//
//   Chamber Material represented by Alu sheet
     tpar[0]= AliMUONConstants::Rmin(0);
     tpar[1]= AliMUONConstants::Rmax(0);
     tpar[2] = (kDGas+kDAlu)/2;
     gMC->Gsvolu("S01A", "TUBE",  idAlu2, tpar, 3);
     gMC->Gsvolu("S02A", "TUBE",idAlu2, tpar, 3);
     gMC->Gspos("S01A", 1, "S01M", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("S02A", 1, "S02M", 0., 0., 0.,  0, "ONLY");
//     
//   Sensitive volumes
     // tpar[2] = kDGas;
     tpar[2] = kDGas/2;
     gMC->Gsvolu("S01G", "TUBE", idGas, tpar, 3);
     gMC->Gsvolu("S02G", "TUBE", idGas, tpar, 3);
     gMC->Gspos("S01G", 1, "S01A", 0., 0., 0.,  0, "ONLY");
     gMC->Gspos("S02G", 1, "S02A", 0., 0., 0.,  0, "ONLY");
//
// Frame Crosses to be placed inside gas
     // NONE: chambers are sensitive everywhere
//      if (frameCrosses) {

// 	 dr = (AliMUONConstants::Rmax(0) - AliMUONConstants::Rmin(0));
// 	 bpar[0] = TMath::Sqrt(dr*dr-dframep1*dframep1/4)/2;
// 	 bpar[1] = dframep1/2;
// 	 bpar[2] = kDGas/2;
// 	 gMC->Gsvolu("S01F", "BOX", idAlu1, bpar, 3);
// 	 gMC->Gsvolu("S02F", "BOX", idAlu1, bpar, 3);
	 
// 	 gMC->Gspos("S01F",1,"S01G", +AliMUONConstants::Rmin(0)+bpar[0] , 0, 0, 
// 		    irot1,"ONLY");
// 	 gMC->Gspos("S01F",2,"S01G", -AliMUONConstants::Rmin(0)-bpar[0] , 0, 0, 
// 		    irot1,"ONLY");
// 	 gMC->Gspos("S01F",3,"S01G", 0, +AliMUONConstants::Rmin(0)+bpar[0] , 0, 
// 		    irot2,"ONLY");
// 	 gMC->Gspos("S01F",4,"S01G", 0, -AliMUONConstants::Rmin(0)-bpar[0] , 0, 
// 		    irot2,"ONLY");
	 
// 	 gMC->Gspos("S02F",1,"S02G", +AliMUONConstants::Rmin(0)+bpar[0] , 0, 0, 
// 		    irot1,"ONLY");
// 	 gMC->Gspos("S02F",2,"S02G", -AliMUONConstants::Rmin(0)-bpar[0] , 0, 0, 
// 		    irot1,"ONLY");
// 	 gMC->Gspos("S02F",3,"S02G", 0, +AliMUONConstants::Rmin(0)+bpar[0] , 0, 
// 		    irot2,"ONLY");
// 	 gMC->Gspos("S02F",4,"S02G", 0, -AliMUONConstants::Rmin(0)-bpar[0] , 0, 
// 		    irot2,"ONLY");
//      }
}

//______________________________________________________________________________
void AliMUONSt1GeometryBuilder::SetTransformations()
{
/// Define the transformations for the station2 chambers.

  Double_t zpos1= - AliMUONConstants::DefaultChamberZ(0); 
  SetTranslation(0, TGeoTranslation(0., 0., zpos1));

  Double_t zpos2 = - AliMUONConstants::DefaultChamberZ(1); 
  SetTranslation(0, TGeoTranslation(0., 0., zpos2));
}

//______________________________________________________________________________
void AliMUONSt1GeometryBuilder::SetSensitiveVolumes()
{
/// Define the sensitive volumes for station1 chambers.

  GetGeometry(0)->SetSensitiveVolume("S01G");
  GetGeometry(1)->SetSensitiveVolume("S02G");
}
