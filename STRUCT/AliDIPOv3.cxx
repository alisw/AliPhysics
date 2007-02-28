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

#include <TVirtualMC.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoMedium.h>
#include <TGeoMatrix.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoCone.h>
#include <TGeoPcon.h>

#include "AliConst.h"
#include "AliDIPOv3.h"
#include "AliMagF.h"
#include "AliRun.h"
 
ClassImp(AliDIPOv3)
 
//_____________________________________________________________________________
AliDIPOv3::AliDIPOv3() 
{
  //
  // Last design of magnetic dipole version 2
  //
}
 
//_____________________________________________________________________________
AliDIPOv3::AliDIPOv3(const char *name, const char *title)
  : AliDIPOv2(name,title)
{
  //
  // Standard constructor for the magnetic dipole version 2    
}


//_____________________________________________________________________________
void AliDIPOv3::CreateSpectrometerDipole()
{
// Detailed dipole geometry as built
//
//
// The top volume
//
    TGeoVolume* top = gGeoManager->GetVolume("ALIC");
//
// Media
//
    TGeoMedium* kMedSteel      = gGeoManager->GetMedium("DIPO_ST_C3");
    TGeoMedium* kMedCoil       = gGeoManager->GetMedium("DIPO_Coil_C1");
    TGeoMedium* kMedCoilSh     = gGeoManager->GetMedium("DIPO_Coil_C3");
    TGeoMedium* kMedCable      = gGeoManager->GetMedium("DIPO_ALU_C2");
    TGeoMedium* kMedAlu        = gGeoManager->GetMedium("DIPO_ALU_C2");
    TGeoMedium* kMedAir        = gGeoManager->GetMedium("DIPO_AIR_MUON");
//
// Rotations
// 
    TGeoRotation* rotxz      = new TGeoRotation("rotxz",    90.,   0., 90.,  90.,  180.,   0.);
    TGeoRotation* rotxz108   = new TGeoRotation("rotxz108", 90., 108., 90., 198.,  180.,   0.);
    TGeoRotation* rotxz180   = new TGeoRotation("rotxz180", 90., 180., 90., 270.,  180.,   0.);
    TGeoRotation* rotxz288   = new TGeoRotation("rotxz288", 90., 288., 90.,  18.,  180.,   0.);

    TGeoRotation* rotxy180   = new TGeoRotation("rotxy180", 90., 180., 90., 270.,    0.,   0.);
    TGeoRotation* rotxy108   = new TGeoRotation("rotxy108", 90., 108., 90., 198.,    0.,   0.);
    TGeoRotation* rotxy288   = new TGeoRotation("rotxy288", 90., 288., 90.,  18.,    0.,   0.);

    TGeoRotation* rot00   = new TGeoRotation("rot00",   180.,   0., 90., 151.,   90.,  61.);
    TGeoRotation* rot01   = new TGeoRotation("rot01",   180.,   0., 90.,  29.,-  90., -61.);
    TGeoRotation* rot02   = new TGeoRotation("rot02",     0.,   0., 90., 151.,   90.,  61.);
    TGeoRotation* rot03   = new TGeoRotation("rot03",     0.,   0., 90.,  29.,-  90., -61.);
    TGeoRotation* rot04   = new TGeoRotation("rot04",    90.,  61., 90., 151.,    0.,   0.);
    TGeoRotation* rot05   = new TGeoRotation("rot05",    90., -61., 90.,-151.,    0.,   0.);
    TGeoRotation* rot06   = new TGeoRotation("rot06",    90., 119., 90., 209.,    0.,   0.);
    TGeoRotation* rot07   = new TGeoRotation("rot07",    90.,-119., 90.,-209.,    0.,   0.);

    const Float_t kZDipole = 975.; 
    Float_t dx, dy, dz;
    //
    // Mother volume
    TGeoPcon* shDDIP = new TGeoPcon(0., 360., 4);
    shDDIP->DefineSection(0, -250.55, 30.1, 570.);
    shDDIP->DefineSection(1,   37.00, 30.1, 570.);
    shDDIP->DefineSection(2,   37.00, ( 37.00 + kZDipole) * TMath::Tan(2. * kDegrad), 570.);
    shDDIP->DefineSection(3,  260.55, (250.55 + kZDipole) * TMath::Tan(2. * kDegrad), 570.);
    TGeoVolume* voDDIP = new TGeoVolume("DDIP", shDDIP, kMedAir);
//
// Yoke
// 
    Float_t yokeLength    = 309.4;
    Float_t blockLength   = yokeLength / 7.;
    Float_t gapWidthFront = 297.6;
    Float_t gapWidthRear  = 395.4;
    Float_t dGap          = (gapWidthRear - gapWidthFront) / 12.;
    Float_t gapHeight     = 609.1;
    Float_t blockHeight   = 147.1;
    

    TGeoVolumeAssembly* asYoke = new TGeoVolumeAssembly("DYoke");	
// Base
    char name[32];
    Float_t lx0 = gapWidthFront + 2. * blockHeight;
    Float_t lx  = lx0;
    
    TGeoVolumeAssembly* asYokeBase = new TGeoVolumeAssembly("DYokeBase");	
    for (Int_t i = 0; i < 7; i++) {
	sprintf(name, "DYokeBaseBlock%1d", i);
	TGeoVolume*  voBaseBlock = new TGeoVolume(name,
						  new TGeoBBox(lx/2., blockHeight/2., blockLength/2.),
						  kMedSteel);
	asYokeBase->AddNode(voBaseBlock, 1, new TGeoTranslation(0., 0., Float_t(i - 3) * blockLength));
	lx += 2. * dGap;
    }

    asYoke->AddNode(asYokeBase, 1, new TGeoTranslation(0., -(gapHeight + blockHeight)/2. , 0.));
    asYoke->AddNode(asYokeBase, 2, new TGeoTranslation(0., +(gapHeight + blockHeight)/2. , 0.));

 
// Side Wall
    TGeoVolumeAssembly* asYokeSide = new TGeoVolumeAssembly("DYokeSide");	
    TGeoVolume*  voSideBlock = new TGeoVolume("DSideBlock",
					      new TGeoBBox(blockHeight/2., gapHeight/2., blockLength/2.),
					      kMedSteel);
    
    for (Int_t i = 0; i < 7; i++) {
	asYokeSide->AddNode(voSideBlock, i, new TGeoTranslation(Float_t(i - 3) * dGap, 0., Float_t(i - 3) * blockLength));
    }


    asYoke->AddNode(asYokeSide, 1, new TGeoTranslation(+lx0/2. + 3. * dGap - blockHeight/2., 0., 0.));
    asYoke->AddNode(asYokeSide, 2, new TGeoCombiTrans( -lx0/2. - 3. * dGap + blockHeight/2., 0., 0., rotxz));

//    
// Coils
//
    Float_t coilRi   = 207.;
    Float_t coilD    =  70.;
    Float_t coilRo   = coilRi + coilD;
    Float_t coilH    =  77.;
    Float_t phiMin   = -61.;
    Float_t phiMax   =  61.;
    Float_t lengthSt = 240. + 48.;    
    Float_t phiKnee  = phiMax * kDegrad;
    Float_t  rKnee   = 25.;
    
//  Circular sections
    TGeoVolumeAssembly* asCoil = new TGeoVolumeAssembly("DCoil");	

    TGeoVolume*  voDC1 = new TGeoVolume("DC1",
					      new TGeoTubeSeg(coilRi, coilRo, coilH / 2., phiMin, phiMax),
					      kMedCoil);
    TGeoVolume*  voDC2 = new TGeoVolume("DC2",
					      new TGeoTubeSeg(coilRi + 5., coilRo - 5., coilH / 2., phiMin, phiMax),
					      kMedCoilSh);

    voDC1->AddNode(voDC2, 1, gGeoIdentity);
    voDC2->SetVisibility(0);

    dz = lengthSt / 2. + coilH / 2. + rKnee;
    dx = 0.;
    
    asCoil->AddNode(voDC1, 1, new TGeoTranslation(-dx, 0., -dz));
    asCoil->AddNode(voDC1, 2, new TGeoCombiTrans(  dx, 0., -dz, rotxy180));
    asCoil->AddNode(voDC1, 3, new TGeoTranslation(-dx, 0.,  dz));
    asCoil->AddNode(voDC1, 4, new TGeoCombiTrans(  dx, 0.,  dz, rotxz180));

    
// 90deg Knees

    
    TGeoVolume* voDC11 = new TGeoVolume("DC11", 
				       new TGeoTubeSeg(rKnee, rKnee + coilH, coilD/2., 270., 360.),
				       kMedCoil);
    
    
    dx = - TMath::Cos(phiKnee) * (coilRi + coilD/2.); 
    dy = - TMath::Sin(phiKnee) * (coilRi + coilD/2.);  
    dz = lengthSt / 2.;

    asCoil->AddNode(voDC11, 1, new TGeoCombiTrans( dx, dy,  -dz, rot00));
    asCoil->AddNode(voDC11, 2, new TGeoCombiTrans( dx, dy,   dz, rot02));
    asCoil->AddNode(voDC11, 3, new TGeoCombiTrans(-dx, dy,  -dz, rot01));
    asCoil->AddNode(voDC11, 4, new TGeoCombiTrans(-dx, dy,   dz, rot03));

    TGeoVolume* voDC12 = new TGeoVolume("DC12", 
				       new TGeoTubeSeg(rKnee, rKnee + coilH, coilD/2., 0., 90.),
				       kMedCoil);
    

    asCoil->AddNode(voDC12, 1, new TGeoCombiTrans( dx, -dy,  -dz, rot01));
    asCoil->AddNode(voDC12, 2, new TGeoCombiTrans( dx, -dy,   dz, rot03));
    asCoil->AddNode(voDC12, 3, new TGeoCombiTrans(-dx, -dy,  -dz, rot00));
    asCoil->AddNode(voDC12, 4, new TGeoCombiTrans(-dx, -dy,   dz, rot02));

// Straight sections

    
    TGeoVolume* voDL0 = new TGeoVolume("DL0", 
				       new TGeoBBox(coilD / 2. + 2., coilH / 2. + 2., lengthSt / 2.),
				       kMedCoil);

    TGeoVolume* voDL1 = new TGeoVolume("DL1", 
				       new TGeoBBox(coilD / 2., coilH / 2., lengthSt / 2.),
				       kMedCoil);
    

    TGeoVolume* voDL2 = new TGeoVolume("DL2", 
				       new TGeoBBox(coilD / 2. - 5., coilH / 2. - 5., lengthSt / 2. - 5.),
				       kMedCoilSh);
    // Sleeves
    TGeoVolume* voDL3 = new TGeoVolume("DL3", 
				       new TGeoBBox(1., coilH / 2., 120.),
				       kMedAlu);

    TGeoVolume* voDL4 = new TGeoVolume("DL4", 
				       new TGeoBBox(coilD/2., 1., 120.),
				       kMedAlu);
    
    voDL0->SetVisibility(0);
    voDL1->AddNode(voDL2, 1, gGeoIdentity);
    voDL0->AddNode(voDL1, 1, gGeoIdentity);
    voDL0->AddNode(voDL3, 1, new TGeoTranslation(-coilD/2. - 1., 0., 0.));
    voDL0->AddNode(voDL3, 2, new TGeoTranslation(+coilD/2. + 1., 0., 0.));
    voDL0->AddNode(voDL4, 1, new TGeoTranslation(0., -coilH/2. - 1., 0.));
    voDL0->AddNode(voDL4, 2, new TGeoTranslation(0., +coilH/2. + 1., 0.));
    
    
    dx += (rKnee + coilH/2.)  * TMath::Sin(phiKnee);
    dy -= (rKnee + coilH/2.)  * TMath::Cos(phiKnee);
    dz = 0.; 
    
    asCoil->AddNode(voDL0, 1, new TGeoCombiTrans( dx,  dy, dz, rot04));
    asCoil->AddNode(voDL0, 2, new TGeoCombiTrans( dx, -dy, dz, rot05));
    asCoil->AddNode(voDL0, 3, new TGeoCombiTrans(-dx,  dy, dz, rot06));
    asCoil->AddNode(voDL0, 4, new TGeoCombiTrans(-dx, -dy, dz, rot07));

// Contactor
// Outer face planes
    
    TGeoVolumeAssembly* asContactor = new TGeoVolumeAssembly("DContactor");
    dx = -5.;
    TGeoVolume* voDC10 = new TGeoVolume("DC10", 
					new TGeoTubeSeg(coilRo + 3.5, coilRo + 83.5, 1., -20., 20.),
					kMedCable);
    asContactor->AddNode(voDC10, 1, new TGeoTranslation(dx, 0, -32.325));
    asContactor->AddNode(voDC10, 2, new TGeoTranslation(dx, 0, +32.325));
 

// Coil Support
// 
    Float_t sW = 89.;
    
    TGeoVolumeAssembly* asDCoilSupport = new TGeoVolumeAssembly("DCoilSupport");

    // Steel fixed to the yoke
    TGeoVolume* voDCS01 = new TGeoVolume("DCS01", 
					 new TGeoTubeSeg(coilRo, 325.,  1., 21., 51.),
					 kMedAlu);
    
    // Steel on the coil
    TGeoVolume* voDCS02 = new TGeoVolume("DCS02", 
					 new TGeoTubeSeg(coilRo, coilRo + 2., sW/2., 21., 51.),
					 kMedAlu);
    TGeoVolume* voDCS021 = new TGeoVolume("DCS021", 
					 new TGeoConeSeg(sW/2., coilRo, 320., coilRo, coilRo + 2., 21., 21.4),
					 kMedAlu);
    

    // Sleeves
    TGeoVolume* voDCS03 = new TGeoVolume("DCS03", 
					 new TGeoTubeSeg(coilRi - 2., coilRo + 2.,  1., 21., 51.),
					 kMedAlu);

    TGeoVolume* voDCS04 = new TGeoVolume("DCS04", 
					 new TGeoTubeSeg(coilRi - 2., coilRi,  coilH/2., 21., 51.),
					 kMedAlu);
    

    TGeoVolume* voDCS05 = new TGeoVolume("DCS05", 
					 new TGeoTubeSeg(coilRi - 2., coilRo,  1., 21., 51.),
					 kMedAlu);
    // 
    
    

    asDCoilSupport->AddNode(voDCS02, 1, new TGeoTranslation(0., 0., -(sW - coilH)/2.));
    asDCoilSupport->AddNode(voDCS04, 1, gGeoIdentity);    
    for (Int_t i = 0; i < 9; i++) 
    {
	char name[16];
	sprintf(name, "rotdcs%1d", i);
	Float_t phi = Float_t(i) * 3.75;
	TGeoRotation* rot   = new TGeoRotation(name, 90., phi, 90., 90. + phi,    0.,   0.);	
	asDCoilSupport->AddNode(voDCS021, i, new TGeoCombiTrans(0., 0., -(sW - coilH)/2., rot));    
    }
    


    asDCoilSupport->AddNode(voDCS01, 1, new TGeoTranslation(0., 0., -sW/2. - (sW - coilH)/2. - 1.));    
    asDCoilSupport->AddNode(voDCS03, 1, new TGeoTranslation(0., 0., +coilH/2. + 1.));    
    asDCoilSupport->AddNode(voDCS05, 1, new TGeoTranslation(0., 0., -coilH/2. - 1.));    

    TGeoVolumeAssembly* asDipole = new TGeoVolumeAssembly("Dipole");
    
    asDipole->AddNode(asYoke, 1, gGeoIdentity);
    asDipole->AddNode(asCoil, 1, gGeoIdentity);
    dz = lengthSt / 2. + coilH / 2. + rKnee;
    asDipole->AddNode(asContactor, 1, new TGeoTranslation(0., 0., dz));
    asDipole->AddNode(asContactor, 2, new TGeoCombiTrans( 0., 0., dz, rotxy180));

    asDipole->AddNode(asDCoilSupport, 1, new TGeoTranslation(0., 0., dz));
    asDipole->AddNode(asDCoilSupport, 2, new TGeoCombiTrans( 0., 0., dz, rotxy180));
    asDipole->AddNode(asDCoilSupport, 3, new TGeoCombiTrans( 0., 0., dz, rotxy108));
    asDipole->AddNode(asDCoilSupport, 4, new TGeoCombiTrans( 0., 0., dz, rotxy288));

    asDipole->AddNode(asDCoilSupport, 5, new TGeoCombiTrans( 0., 0., -dz, rotxz));
    asDipole->AddNode(asDCoilSupport, 6, new TGeoCombiTrans( 0., 0., -dz, rotxz108));
    asDipole->AddNode(asDCoilSupport, 7, new TGeoCombiTrans( 0., 0., -dz, rotxz180));
    asDipole->AddNode(asDCoilSupport, 8, new TGeoCombiTrans( 0., 0., -dz, rotxz288));

    voDDIP->SetVisContainers();
    voDDIP->SetVisibility(0);
    voDDIP->AddNode(asDipole, 1, new TGeoTranslation(0., 0., -0.55));
    top->AddNode(voDDIP, 1, new TGeoCombiTrans(0., 0., -kZDipole, rotxz));

}







