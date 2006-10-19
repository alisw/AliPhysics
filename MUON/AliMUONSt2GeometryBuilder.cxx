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
// Class AliMUONSt2GeometryBuilder
// -------------------------------
// MUON Station2 coarse geometry construction class.
// Extracted from AliMUONv1
// Dummy version of station 2 with the right DE id (Ch. Finck)

#include <TVirtualMC.h>
#include <TGeoMatrix.h>

#include "AliLog.h"

#include "AliMUONSt2GeometryBuilder.h"
#include "AliMUON.h"
#include "AliMUONConstants.h"
#include "AliMUONGeometryModule.h"
#include "AliMUONGeometryEnvelopeStore.h"

/// \cond CLASSIMP
ClassImp(AliMUONSt2GeometryBuilder)
/// \endcond

//______________________________________________________________________________
AliMUONSt2GeometryBuilder::AliMUONSt2GeometryBuilder(AliMUON* muon)
 : AliMUONVGeometryBuilder(2, 2), 
   fMUON(muon)
{
// Standard constructor

}

//______________________________________________________________________________
AliMUONSt2GeometryBuilder::AliMUONSt2GeometryBuilder()
 : AliMUONVGeometryBuilder(),
   fMUON(0)
{
// Default constructor
}


//______________________________________________________________________________
AliMUONSt2GeometryBuilder::~AliMUONSt2GeometryBuilder() {
//
}

//
// public methods
//

//______________________________________________________________________________
void AliMUONSt2GeometryBuilder::CreateGeometry() 
{
// From AliMUONv1::CreateGeometry()

//
//********************************************************************
//                            Station 2                             **
//********************************************************************
     // indices 1 and 2 for first and second chambers in the station
     // iChamber (first chamber) kept for other quanties than Z,
     // assumed to be the same in both chambers

     // Get tracking medias Ids     
     Int_t *idtmed = fMUON->GetIdtmed()->GetArray()-1099;
     Int_t idAir= idtmed[1100]; // medium 1
  //     Int_t idAlu1=idtmed[1103]; // medium 4
     //     Int_t idAlu2=idtmed[1104]; // medium 5
     Int_t idGas=idtmed[1108];  // medium 9 = Ar-CO2 gas (80%+20%)

     const Float_t kDeltaQuad = 0.01;//2.6; dummy value til we find the good value
     const Float_t kDeltaZ = 6.5/2.;

     // Rotation matrices in the x-y plane  
     // phi= 0 deg
     Int_t irot1;
     fMUON->AliMatrix(irot1,  90.,   0., 90.,  90., 0., 0.);
     // phi= 90 deg
     Int_t irot2;
     fMUON->AliMatrix(irot2,  90.,  90., 90., 180., 0., 0.);

     // Half of the total thickness of frame crosses (including DAlu)
     // for each chamber in stations 1 and 2:
     // 3% of X0 of composite material,
     // but taken as Aluminium here, with same thickness in number of X0
     //     Float_t dframez = 3. * 8.9 / 100;
     // DGas and DAlu not changed from standard values
     //     Double_t zfpos=-(iChamber->DGas()+dframez+iChamber->DAlu())/2;
             // The same parameters are defined in builder for station 1 

     //    sensitive gas gap
     const Float_t kDGas = 0.5;

     //    3% radiation length of aluminum (X0=8.9 cm)      
     // const Float_t kDAlu = 3.5 * 8.9 / 100.;
     
     // Mother volume
     // Outer excess and inner recess for mother volume radius
     // with respect to ROuter and RInner
     //     Float_t dframep=.001; // Value for station 3 should be 6 ...
     // Width (RdPhi) of the frame crosses for stations 1 and 2 (cm)
     // Float_t dframep1=.001;
     //     Float_t phi=2*TMath::Pi()/12/2;
             // The same parameters are defined in builder for station 1 (30deg)

     Float_t tpar[5];
     //     Double_t dstation = (-iChamber2->Z()) - (-iChamber1->Z());


     Float_t posx, posy, posz;
//   Chamber Material represented by Alu sheet
     tpar[0]= AliMUONConstants::Rmin(1);
     tpar[1]= AliMUONConstants::Rmax(1);
     tpar[2] = kDGas/2;
     tpar[3] = 0.;
     tpar[4] = 90.;
  
//     
//   Sensitive volumes
     gMC->Gsvolu("S03G", "TUBS", idGas, tpar, 5);
     gMC->Gsvolu("S04G", "TUBS", idGas, tpar, 5);

     Int_t detElemId;

     posx = kDeltaQuad;
     posy = kDeltaQuad;
     posz = -kDeltaZ;

     detElemId = 301;
     gMC->Gsvolu("LE01", "TUBS", idAir, tpar, 5);
     GetEnvelopes(2)->AddEnvelope("LE01", detElemId, true, TGeoTranslation(posx, posy, posz),
				   TGeoRotation("rot1",90,0,90,90,0,0) );
     detElemId = 401;
     gMC->Gsvolu("LF01", "TUBS", idAir, tpar, 5);
     GetEnvelopes(3)->AddEnvelope("LF01", detElemId, true, TGeoTranslation(posx, posy, posz),
				   TGeoRotation("rot1",90,0,90,90,0,0) );
     detElemId = 300;
     gMC->Gsvolu("LE02", "TUBS", idAir, tpar, 5);
     GetEnvelopes(2)->AddEnvelope("LE02", detElemId, true, TGeoTranslation(-posx, posy,-posz),
				  TGeoRotation("rot2",90,180,90,90,180,0) );
     detElemId = 400;
     gMC->Gsvolu("LF02", "TUBS", idAir, tpar, 5);
     GetEnvelopes(3)->AddEnvelope("LF02", detElemId, true, TGeoTranslation(-posx, posy,-posz),
				  TGeoRotation("rot2",90,180,90,90,180,0) );
     detElemId = 302;
     gMC->Gsvolu("LE03", "TUBS", idAir, tpar, 5);
     GetEnvelopes(2)->AddEnvelope("LE03", detElemId, true, TGeoTranslation(posx, -posy, -posz),
				  TGeoRotation("rot3",90,0,90,270,180,0) );
     detElemId = 402;
     gMC->Gsvolu("LF03", "TUBS", idAir, tpar, 5);
     GetEnvelopes(3)->AddEnvelope("LF03", detElemId, true, TGeoTranslation(posx, -posy, -posz),
				  TGeoRotation("rot3",90,0,90,270,180,0) );
     detElemId = 303;
     gMC->Gsvolu("LE04", "TUBS", idAir, tpar, 5);
     GetEnvelopes(2)->AddEnvelope("LE04", detElemId, true, TGeoTranslation(-posx, -posy, posz),
				  TGeoRotation("rot4",90,180,90,270,0,0) );
     detElemId = 403;
     gMC->Gsvolu("LF04", "TUBS", idAir, tpar, 5);
     GetEnvelopes(3)->AddEnvelope("LF04", detElemId, true, TGeoTranslation(-posx, -posy, posz),
				  TGeoRotation("rot4",90,180,90,270,0,0) );

     GetEnvelopes(2)->AddEnvelopeConstituent("S03G", "LE01", 1);
     GetEnvelopes(3)->AddEnvelopeConstituent("S04G", "LF01", 1);

     GetEnvelopes(2)->AddEnvelopeConstituent("S03G", "LE02", 2);
     GetEnvelopes(3)->AddEnvelopeConstituent("S04G", "LF02", 2);

     GetEnvelopes(2)->AddEnvelopeConstituent("S03G", "LE03", 3);
     GetEnvelopes(3)->AddEnvelopeConstituent("S04G", "LF03", 3);

     GetEnvelopes(2)->AddEnvelopeConstituent("S03G", "LE04", 4);
     GetEnvelopes(3)->AddEnvelopeConstituent("S04G", "LF04", 4);
}

//______________________________________________________________________________
void AliMUONSt2GeometryBuilder::SetTransformations()
{
// Defines the transformations for the station2 chambers.
// ---

  // Define chamber volumes as virtual
  SetVolume(2, "SC03", true);
  SetVolume(3, "SC04", true);

  Double_t zpos1 = - AliMUONConstants::DefaultChamberZ(2); 
  SetTranslation(2, TGeoTranslation(0., 0., zpos1));

  Double_t zpos2 = - AliMUONConstants::DefaultChamberZ(3); 
  SetTranslation(3, TGeoTranslation(0., 0., zpos2));
}

//______________________________________________________________________________
void AliMUONSt2GeometryBuilder::SetSensitiveVolumes()
{
// Defines the sensitive volumes for station2 chambers.
// ---

  GetGeometry(2)->SetSensitiveVolume("S03G");
  GetGeometry(3)->SetSensitiveVolume("S04G");
}
