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
// Class AliMUONSlatGeometryBuilder
// -------------------------------
// Abstract base class for geometry construction per chamber.
//
// Author: Eric Dumonteil (dumontei@cea.fr)


// This Builder is designed according to the enveloppe methode. The basic idea is to be able to allow moves 
// of the slats on the support panels. 
// Those moves can be described with a simple set of parameters. The next step should be now to describe all 
// the slats and their places by a unique 
// class, which would make the SlatBuilder far more compact since now only three parameters can define a slat 
// and its position, like:
//   * Bool_t rounded_shape_slat
//   * Float_t slat_length
//   * Float_t slat_number or Float_t slat_position

#include <TVirtualMC.h>
#include <TGeoMatrix.h>
#include <Riostream.h>

#include "AliRun.h"
#include "AliLog.h"

#include "AliMUONSlatGeometryBuilder.h"
#include "AliMUON.h"
#include "AliMUONChamber.h"
#include "AliMUONGeometryModule.h"
#include "AliMUONGeometryEnvelopeStore.h"
#include "AliMUONConstants.h"

ClassImp(AliMUONSlatGeometryBuilder)


//______________________________________________________________________________
AliMUONSlatGeometryBuilder::AliMUONSlatGeometryBuilder(AliMUON* muon)
 : AliMUONVGeometryBuilder("slat.dat",
                           muon->Chamber(4).GetGeometry(), 
			   muon->Chamber(5).GetGeometry(), 
                           muon->Chamber(6).GetGeometry(), 
			   muon->Chamber(7).GetGeometry(), 
			   muon->Chamber(8).GetGeometry(), 
			   muon->Chamber(9).GetGeometry()),
   fMUON(muon)
{
// Standard constructor

}

//______________________________________________________________________________
AliMUONSlatGeometryBuilder::AliMUONSlatGeometryBuilder() 
 : AliMUONVGeometryBuilder(),
   fMUON(0)
{
// Default constructor
}


//______________________________________________________________________________
AliMUONSlatGeometryBuilder::AliMUONSlatGeometryBuilder(const AliMUONSlatGeometryBuilder& rhs)
  : AliMUONVGeometryBuilder(rhs)
{
  AliFatal("Copy constructor is not implemented.");
}

//______________________________________________________________________________
AliMUONSlatGeometryBuilder::~AliMUONSlatGeometryBuilder() {
//
}

//______________________________________________________________________________
AliMUONSlatGeometryBuilder& 
AliMUONSlatGeometryBuilder::operator = (const AliMUONSlatGeometryBuilder& rhs) 
{
  // check assignement to self
  if (this == &rhs) return *this;

  AliFatal("Assignment operator is not implemented.");
    
  return *this;  
}

//
// public methods
//

//______________________________________________________________________________
void AliMUONSlatGeometryBuilder::CreateGeometry()
{
  // CreateGeometry is the method containing all the informations concerning Stations 345 geometry.
  // It includes description and placements of support panels and slats.
  // The code comes directly from what was written in AliMUONv1.cxx before, with modifications concerning 
  // the use of Enveloppe method to place the Geant volumes.
  // Now, few changes would allow the creation of a Slat methode where slat could be described by few parameters, 
  // and this builder would then be dedicated only to the
  // placements of the slats. Those modifications could shorten the Station 345 geometry by a non-negligeable factor...
 
  Int_t *idtmed = fMUON->GetIdtmed()->GetArray()-1099;

  Float_t angle;
  Float_t *dum=0;

  // define the id of tracking media:
  Int_t idAir    = idtmed[1100]; // medium 1
  Int_t idGas    = idtmed[1108]; // medium 9 = Ar-CO2 gas (80%+20%)
  Int_t idCopper = idtmed[1110];
  Int_t idG10    = idtmed[1111];
  Int_t idCarbon = idtmed[1112];
  Int_t idRoha   = idtmed[1113];
  Int_t idNomex  = idtmed[1114]; // honey comb
  Int_t idNoryl  = idtmed[1115]; 
  Int_t idNomexB = idtmed[1116]; // bulk material 

  // sensitive area: 40*40 cm**2
  const Float_t kSensLength = 40.; 
  const Float_t kSensHeight = 40.; 
  const Float_t kSensWidth  = 0.5; // according to TDR fig 2.120 
  const Int_t kSensMaterial = idGas;
  //     const Float_t kYoverlap   = 1.5; 

  // PCB dimensions in cm; width: 30 mum copper   
  const Float_t kPcbLength  = kSensLength; 
  const Float_t kPcbHeight  = 58.; // updated Ch. Finck 
  const Float_t kPcbWidth   = 0.003; 
  const Int_t kPcbMaterial  = idCopper;

  // Insulating material: 220 mum G10 fiber  glued to pcb  
  const Float_t kInsuLength = kPcbLength; 
  const Float_t kInsuHeight = kPcbHeight; 
  const Float_t kInsuWidth  = 0.022;  // updated Ch. Finck 
  const Int_t kInsuMaterial = idG10;

  // Carbon fiber panels: 200mum carbon/epoxy skin   
  const Float_t kCarbonWidth  = 0.020;      
  const Int_t kCarbonMaterial = idCarbon;

  // Nomex (honey comb) between the two panel carbon skins    
  const Float_t kNomexLength = kSensLength; 
  const Float_t kNomexHeight = kSensHeight; 
  const Float_t kNomexWidth  = 0.8; // updated Ch. Finck 
  const Int_t kNomexMaterial = idNomex;
 
  // Bulk Nomex under panel sandwich Ch. Finck    
  const Float_t kNomexBWidth  = 0.025; 
  const Int_t kNomexBMaterial = idNomexB;

  // Panel sandwich 0.02 carbon*2 + 0.8 nomex     
  const Float_t kPanelLength = kSensLength; 
  const Float_t kPanelHeight = kSensHeight; 
  const Float_t kPanelWidth  = 2 * kCarbonWidth + kNomexWidth;

  // spacer around the slat: 2 sticks along length,2 along height  
  // H: the horizontal ones 
  const Float_t kHframeLength = kPcbLength; 
  const Float_t kHframeHeight = 1.95; // updated Ch. Finck 
  const Float_t kHframeWidth  = kSensWidth; 
  const Int_t kHframeMaterial = idNoryl;

  // V: the vertical ones; vertical spacers 
  const Float_t kVframeLength = 2.5; 
  const Float_t kVframeHeight = kSensHeight + kHframeHeight; 
  const Float_t kVframeWidth  = kSensWidth;
  const Int_t kVframeMaterial = idNoryl;

  // B: the horizontal border filled with rohacell: ok Ch. Finck
  const Float_t kBframeLength = kHframeLength; 
  const Float_t kBframeHeight = (kPcbHeight - kSensHeight)/2. - kHframeHeight; 
  const Float_t kBframeWidth  = kHframeWidth;
  const Int_t kBframeMaterial = idRoha;

  // NULOC: 30 mum copper + 200 mum vetronite (same radiation length as 14mum copper) for electronics
  const Float_t kNulocLength   = 2.5; 
  const Float_t kNulocHeight   = kBframeHeight;
  const Float_t kNulocWidth    = 0.0030 + 0.0014; // equivalent copper width of vetronite; 
  const Int_t   kNulocMaterial = idCopper;

  // Slat parameters
  const Float_t kSlatHeight = kPcbHeight; 
  const Float_t kSlatWidth  = kSensWidth + 2.*(kPcbWidth + kInsuWidth + kPanelWidth 
					       + kNomexBWidth); //replaced rohacell with Nomex Ch. Finck 
  const Int_t   kSlatMaterial = idAir;
  const Float_t kDslatLength  = -1.25; // position of the slat respect to the beam plane (half vertical spacer) Ch. Finck
  Float_t zSlat               = AliMUONConstants::DzSlat();// implemented Ch. Finck
  Float_t dzCh                = AliMUONConstants::DzCh();

  Float_t spar[3];  
  Int_t i, j;
  Int_t detElemId;

  // the panel volume contains the nomex
  Float_t panelpar[3] = { kPanelLength/2., kPanelHeight/2., kPanelWidth/2. }; 
  Float_t nomexpar[3] = { kNomexLength/2., kNomexHeight/2., kNomexWidth/2. }; 
  Float_t twidth =  kPanelWidth +  kNomexBWidth; 
  Float_t nomexbpar[3] = {kNomexLength/2., kNomexHeight/2.,twidth/2. };// bulk nomex 

  // insulating material contains PCB-> gas   
  twidth = 2*(kInsuWidth + kPcbWidth) + kSensWidth ; 
  Float_t insupar[3] = {kInsuLength/2., kInsuHeight/2., twidth/2. }; 
  twidth -= 2 * kInsuWidth; 
  Float_t pcbpar[3]  = {kPcbLength/2., kPcbHeight/2., twidth/2. }; 
  Float_t senspar[3] = {kSensLength/2., kSensHeight/2., kSensWidth/2. }; 
  Float_t theight    = 2 * kHframeHeight + kSensHeight;
  Float_t hFramepar[3] = {kHframeLength/2., theight/2., kHframeWidth/2.}; 
  Float_t bFramepar[3] = {kBframeLength/2., kBframeHeight/2., kBframeWidth/2.}; 
  Float_t vFramepar[3] = {kVframeLength/2., kVframeHeight/2., kVframeWidth/2.};
  Float_t nulocpar[3]  = {kNulocLength/2.,  kNulocHeight/2.,  kNulocWidth/2.}; 

  Float_t xx;
  Float_t xxmax = (kBframeLength - kNulocLength)/2.; 
  Int_t index=0;
      
  AliMUONChamber *iChamber, *iChamber1, *iChamber2;

  Int_t* fStations = new Int_t[5];
  for (Int_t i=0; i<5; i++) fStations[i] = 1;
  fStations[2] = 1;
     
  if (fStations[2])
    {
      //********************************************************************
      //                            Station 3                             **
      //********************************************************************
      // indices 1 and 2 for first and second chambers in the station
      // iChamber (first chamber) kept for other quanties than Z,
      // assumed to be the same in both chambers

      iChamber = &fMUON->Chamber(4);
      iChamber1 = iChamber;
      iChamber2 = &fMUON->Chamber(5);
     
      //iChamber1->GetGeometry()->SetDebug(kTRUE);
      //iChamber2->GetGeometry()->SetDebug(kTRUE);
 
      if (gAlice->GetModule("DIPO")) {
	// if DIPO is preset, the whole station will be placed in DDIP volume
	iChamber1->GetGeometry()->SetMotherVolume("DDIP");
	iChamber2->GetGeometry()->SetMotherVolume("DDIP");
      }


      // volumes for slat geometry (xx=5,..,10 chamber id): 
      // Sxx0 Sxx1 Sxx2 Sxx3  -->   Slat Mother volumes 
      // SxxG                          -->   Sensitive volume (gas)
      // SxxP                          -->   PCB (copper) 
      // SxxI                          -->   Insulator (G10) 
      // SxxC                          -->   Carbon panel 
      // SxxN                          -->   Nomex comb
      // SxxX                          -->   Nomex bulk
      // SxxH, SxxV                    -->   Horizontal and Vertical frames (Noryl)
      // SB5x                          -->   Volumes for the 35 cm long PCB
      // slat dimensions: slat is a MOTHER volume!!! made of air

      // only for chamber 5: slat 1 has a PCB shorter by 5cm!

      Float_t tlength = 35.;
      Float_t panelpar2[3]  = { tlength/2., panelpar[1],  panelpar[2]}; 
      Float_t nomexpar2[3]  = { tlength/2., nomexpar[1],  nomexpar[2]}; 
      Float_t nomexbpar2[3] = { tlength/2., nomexbpar[1],  nomexbpar[2]}; 
      Float_t insupar2[3]   = { tlength/2., insupar[1],   insupar[2]}; 
      Float_t pcbpar2[3]    = { tlength/2., pcbpar[1],    pcbpar[2]}; 
      Float_t senspar2[3]   = { tlength/2., senspar[1],   senspar[2]}; 
      Float_t hFramepar2[3] = { tlength/2., hFramepar[1], hFramepar[2]}; 
      Float_t bFramepar2[3] = { tlength/2., bFramepar[1], bFramepar[2]}; 
      Float_t *dum=0;
      Float_t pcbDLength3   = (kPcbLength - tlength);

      const Int_t   kNslats3         = 5;  // number of slats per quadrant
      const Int_t   kNPCB3[kNslats3] = {4, 4, 4, 3, 2}; // n PCB per slat
      const Float_t kXpos3[kNslats3] = {0., 0., 0., 0., 0.};//{31., 0., 0., 0., 0.};
      const Float_t kYpos3[kNslats3] = {0, 37.8, 37.7, 37.3, 33.7};
      Float_t slatLength3[kNslats3]; 

      // create and position the slat (mother) volumes 

      char idSlatCh5[5];
      char idSlatCh6[5];
      Float_t xSlat3;
      Float_t ySlat3 = 0;
      Float_t angle = 0.;
      Float_t spar2[3];
      for (i = 0; i < kNslats3; i++){

	slatLength3[i] = kPcbLength * kNPCB3[i] + 2.* kVframeLength; 
	xSlat3 = slatLength3[i]/2. +  kDslatLength + kXpos3[i]; 
	ySlat3 += kYpos3[i];

	spar[0] = slatLength3[i]/2.; 
	spar[1] = kSlatHeight/2.;
	spar[2] = kSlatWidth/2.; 
	// take away 5 cm from the first slat in chamber 5
        if (i == 0 || i == 1 || i == 2) { // 1 pcb is shortened by 5cm
	  spar2[0] = spar[0] - pcbDLength3/2.;
	} else {
	  spar2[0] = spar[0];
	}
	spar2[1] = spar[1];
	spar2[2] = spar[2]; 
	Float_t dzCh3 = dzCh; 
	Float_t zSlat3 = (i%2 ==0)? -zSlat : zSlat; // seems not that zSlat3 = zSlat4 & 5 refering to plan PQ7EN345-6 ?

	sprintf(idSlatCh5,"LA%d",kNslats3-1+i);
	gMC->Gsvolu(idSlatCh5,"BOX",kSlatMaterial,spar2,3);
	detElemId = 500 + i + kNslats3-1;
	GetEnvelopes(4)->AddEnvelope(idSlatCh5, detElemId, true, TGeoTranslation(xSlat3, ySlat3, -zSlat3 + dzCh3),
				     TGeoRotation("rot1",90,angle,90,90+angle,0,0) );

	sprintf(idSlatCh5,"LA%d",3*kNslats3-2+i);
	gMC->Gsvolu(idSlatCh5,"BOX",kSlatMaterial,spar2,3);
	detElemId = 550 + i + kNslats3-1;
	GetEnvelopes(4)->AddEnvelope(idSlatCh5, detElemId, true, TGeoTranslation(-xSlat3, ySlat3, zSlat3 - dzCh3),
				     TGeoRotation("rot2",90,180+angle,90,90+angle,180,0) );

	if (i > 0) { 
	  sprintf(idSlatCh5,"LA%d",kNslats3-1-i);
	  gMC->Gsvolu(idSlatCh5,"BOX",kSlatMaterial,spar2,3);
	  detElemId = 500 - i + kNslats3-1;
	  GetEnvelopes(4)->AddEnvelope(idSlatCh5, detElemId, true, TGeoTranslation(xSlat3, -ySlat3, -zSlat3 + dzCh3), 
				       TGeoRotation("rot3",90,angle,90,270+angle,180,0) );

	  sprintf(idSlatCh5,"LA%d",3*kNslats3-2-i);
	  gMC->Gsvolu(idSlatCh5,"BOX",kSlatMaterial,spar2,3);
	  detElemId = 550 - i + kNslats3-1;
	  GetEnvelopes(4)->AddEnvelope(idSlatCh5, detElemId, true, TGeoTranslation(-xSlat3, -ySlat3, zSlat3 - dzCh3),
				       TGeoRotation("rot4",90,180+angle,90,270+angle,0,0) );
	}

	sprintf(idSlatCh6,"LB%d",kNslats3-1+i);  
	gMC->Gsvolu(idSlatCh6,"BOX",kSlatMaterial,spar,3);
	detElemId = 600 + i  + kNslats3-1;
	GetEnvelopes(5)->AddEnvelope(idSlatCh6, detElemId, true, TGeoTranslation(xSlat3, ySlat3, -zSlat3 + dzCh3),
				     TGeoRotation("rot5",90,angle,90,90+angle,0,0) );
	sprintf(idSlatCh6,"LB%d",3*kNslats3-2+i);
	gMC->Gsvolu(idSlatCh6,"BOX",kSlatMaterial,spar,3);
	detElemId = 650 + i + kNslats3-1;
	GetEnvelopes(5)->AddEnvelope(idSlatCh6, detElemId, true, TGeoTranslation(-xSlat3, ySlat3, zSlat3 - dzCh3),
				     TGeoRotation("rot6",90,180+angle,90,90+angle,180,0) );

	if (i > 0) { 
	  sprintf(idSlatCh6,"LB%d",kNslats3-1-i);
	  gMC->Gsvolu(idSlatCh6,"BOX",kSlatMaterial,spar,3);
	  detElemId = 600 - i + kNslats3-1;
	  GetEnvelopes(5)->AddEnvelope(idSlatCh6, detElemId, true, TGeoTranslation(xSlat3, -ySlat3, -zSlat3 + dzCh3),
				       TGeoRotation("rot7",90,angle,90,270+angle,180,0) );

	  sprintf(idSlatCh6,"LB%d",3*kNslats3-2-i);
	  gMC->Gsvolu(idSlatCh6,"BOX",kSlatMaterial,spar,3);
	  detElemId = 650 - i + kNslats3-1;
	  GetEnvelopes(5)->AddEnvelope(idSlatCh6, detElemId, true, TGeoTranslation(-xSlat3, -ySlat3, zSlat3 - dzCh3),
				       TGeoRotation("rot8",90,180+angle,90,270+angle,0,0) );
	}
      }
     
      // create the panel volume 
 
      gMC->Gsvolu("S05C","BOX",kCarbonMaterial,panelpar,3);
      gMC->Gsvolu("SB5C","BOX",kCarbonMaterial,panelpar2,3);
      gMC->Gsvolu("S06C","BOX",kCarbonMaterial,panelpar,3);
 
      // create the nomex volume (honey comb)

      gMC->Gsvolu("S05N","BOX",kNomexMaterial,nomexpar,3);
      gMC->Gsvolu("SB5N","BOX",kNomexMaterial,nomexpar2,3);
      gMC->Gsvolu("S06N","BOX",kNomexMaterial,nomexpar,3);
 
      // create the nomex volume (bulk)

      gMC->Gsvolu("S05X","BOX",kNomexBMaterial,nomexbpar,3);
      gMC->Gsvolu("SB5X","BOX",kNomexBMaterial,nomexbpar2,3);
      gMC->Gsvolu("S06X","BOX",kNomexBMaterial,nomexbpar,3);

      // create the insulating material volume 

      gMC->Gsvolu("S05I","BOX",kInsuMaterial,insupar,3);
      gMC->Gsvolu("SB5I","BOX",kInsuMaterial,insupar2,3);
      gMC->Gsvolu("S06I","BOX",kInsuMaterial,insupar,3);
 
      // create the PCB volume 

      gMC->Gsvolu("S05P","BOX",kPcbMaterial,pcbpar,3);
      gMC->Gsvolu("SB5P","BOX",kPcbMaterial,pcbpar2,3);
      gMC->Gsvolu("S06P","BOX",kPcbMaterial,pcbpar,3);
 
      // create the sensitive volumes,

      gMC->Gsvolu("S05G","BOX",kSensMaterial,dum,0);
      gMC->Gsvolu("S06G","BOX",kSensMaterial,dum,0);

      // create the vertical frame volume 

      gMC->Gsvolu("S05V","BOX",kVframeMaterial,vFramepar,3);
      gMC->Gsvolu("S06V","BOX",kVframeMaterial,vFramepar,3);

      // create the horizontal frame volume 

      gMC->Gsvolu("S05H","BOX",kHframeMaterial,hFramepar,3);
      gMC->Gsvolu("SB5H","BOX",kHframeMaterial,hFramepar2,3);
      gMC->Gsvolu("S06H","BOX",kHframeMaterial,hFramepar,3);
 
      // create the horizontal border volume 

      gMC->Gsvolu("S05B","BOX",kBframeMaterial,bFramepar,3);
      gMC->Gsvolu("SB5B","BOX",kBframeMaterial,bFramepar2,3);
      gMC->Gsvolu("S06B","BOX",kBframeMaterial,bFramepar,3);
 
      index = 0; 
      for (i = 0; i<kNslats3; i++){
	for (Int_t quadrant = 1; quadrant <= 4; quadrant++) {

	  if (i == 0 && quadrant == 2) continue;
	  if (i == 0 && quadrant == 4) continue;

	  sprintf(idSlatCh5,"LA%d",ConvertSlatNum(i,quadrant,kNslats3-1));
	  sprintf(idSlatCh6,"LB%d",ConvertSlatNum(i,quadrant,kNslats3-1));
	  Float_t xvFrame  = (slatLength3[i] - kVframeLength)/2.;
	  Float_t xvFrame2  = xvFrame;

	  if (i == 0 || i == 1 || i == 2) xvFrame2 -= pcbDLength3/2.;

	  // position the vertical frames 
	  if ( i > 2) { 
	    GetEnvelopes(4)->AddEnvelopeConstituent("S05V", idSlatCh5, 
						    (2*i-1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
	    GetEnvelopes(4)->AddEnvelopeConstituent("S05V", idSlatCh5, 
						    (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));
	    GetEnvelopes(5)->AddEnvelopeConstituent("S06V", idSlatCh6, 
						    (2*i-1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
	    GetEnvelopes(5)->AddEnvelopeConstituent("S06V", idSlatCh6, 
						    (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));	  
	  } 

	  if (i == 2) {
	    GetEnvelopes(4)->AddEnvelopeConstituent("S05V", idSlatCh5, 
						    (2*i-1)*10+quadrant,TGeoTranslation(xvFrame2,0.,0.));
	    GetEnvelopes(4)->AddEnvelopeConstituent("S05V", idSlatCh5, 
						    (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));
	    GetEnvelopes(5)->AddEnvelopeConstituent("S06V", idSlatCh6, 
						    (2*i-1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
	    GetEnvelopes(5)->AddEnvelopeConstituent("S06V", idSlatCh6, 
						    (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));
	  }

	  if (i == 0 || i == 1) { // no rounded spacer for the moment (Ch. Finck)
	    GetEnvelopes(4)->AddEnvelopeConstituent("S05V", idSlatCh5, 
						    (2*i-1)*10+quadrant,TGeoTranslation(xvFrame2,0.,0.));
	    GetEnvelopes(5)->AddEnvelopeConstituent("S06V", idSlatCh6, 
						    (2*i-1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
	  }

	  // position the panels and the insulating material 
	  for (j = 0; j < kNPCB3[i]; j++){
	    if (i == 1 && j == 0) continue;
	    if (i == 0 && j == 0) continue;
	    index++;
	    Float_t xx = kSensLength * (-kNPCB3[i]/2. + j + 0.5); 
	    Float_t xx2 = xx - pcbDLength3/2.; 
	 
	    Float_t zPanel = spar[2] - nomexbpar[2]; 

	    if ( (i == 0 || i == 1 || i == 2) && j == kNPCB3[i]-1) { // 1 pcb is shortened by 5cm 
	      GetEnvelopes(4)->AddEnvelopeConstituent("SB5X", idSlatCh5, 2*index-1,TGeoTranslation(xx2,0.,zPanel));
	      GetEnvelopes(4)->AddEnvelopeConstituent("SB5X", idSlatCh5, 2*index,TGeoTranslation(xx2,0.,-zPanel));
	      GetEnvelopes(4)->AddEnvelopeConstituent("SB5I", idSlatCh5, index,TGeoTranslation(xx2,0.,0.));
	    } else {
	      GetEnvelopes(4)->AddEnvelopeConstituent("S05X", idSlatCh5, 2*index-1,TGeoTranslation(xx,0.,zPanel));
	      GetEnvelopes(4)->AddEnvelopeConstituent("S05X", idSlatCh5, 2*index,TGeoTranslation(xx,0.,-zPanel));
	      GetEnvelopes(4)->AddEnvelopeConstituent("S05I", idSlatCh5, index,TGeoTranslation(xx,0.,0.));
	    }
	    GetEnvelopes(5)->AddEnvelopeConstituent("S06X", idSlatCh6, 2*index-1,TGeoTranslation(xx,0.,zPanel));
	    GetEnvelopes(5)->AddEnvelopeConstituent("S06X", idSlatCh6, 2*index,TGeoTranslation(xx,0.,-zPanel));
	    GetEnvelopes(5)->AddEnvelopeConstituent("S06I", idSlatCh6, index,TGeoTranslation(xx,0.,0.));
 
	  } 
	}
      }

      // position the nomex volume inside the panel volume
      gMC->Gspos("S05N",1,"S05C",0.,0.,0.,0,"ONLY"); 
      gMC->Gspos("SB5N",1,"SB5C",0.,0.,0.,0,"ONLY"); 
      gMC->Gspos("S06N",1,"S06C",0.,0.,0.,0,"ONLY"); 
  
      // position panel volume inside the bulk nomex material volume
      gMC->Gspos("S05C",1,"S05X",0.,0.,kNomexBWidth/2.,0,"ONLY"); 
      gMC->Gspos("SB5C",1,"SB5X",0.,0.,kNomexBWidth/2.,0,"ONLY"); 
      gMC->Gspos("S06C",1,"S06X",0.,0.,kNomexBWidth/2.,0,"ONLY"); 

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
      Float_t yborder = ( kPcbHeight - kBframeHeight ) / 2.; 
      gMC->Gspos("S05B",1,"S05P",0., yborder,0.,0,"ONLY"); 
      gMC->Gspos("S05B",2,"S05P",0.,-yborder,0.,0,"ONLY"); 
      gMC->Gspos("SB5B",1,"SB5P",0., yborder,0.,0,"ONLY"); 
      gMC->Gspos("SB5B",2,"SB5P",0.,-yborder,0.,0,"ONLY"); 

      gMC->Gspos("S06B",1,"S06P",0., yborder,0.,0,"ONLY"); 
      gMC->Gspos("S06B",2,"S06P",0.,-yborder,0.,0,"ONLY"); 
  
      // create the NULOC volume and position it in the horizontal frame
      gMC->Gsvolu("S05E","BOX",kNulocMaterial,nulocpar,3);
      gMC->Gsvolu("S06E","BOX",kNulocMaterial,nulocpar,3);
      index = 0;
      Float_t xxmax2 = xxmax - pcbDLength3/2.;
      for (xx = -xxmax; xx <= xxmax; xx += 2*kNulocLength) { 
	index++; 
	gMC->Gspos("S05E",2*index-1,"S05B", xx, 0.,-kBframeWidth/2. + kNulocWidth/2, 0, "ONLY");
	gMC->Gspos("S05E",2*index  ,"S05B", xx, 0., kBframeWidth/2. - kNulocWidth/2, 0, "ONLY");
	gMC->Gspos("S06E",2*index-1,"S06B", xx, 0.,-kBframeWidth/2. + kNulocWidth/2, 0, "ONLY");
	gMC->Gspos("S06E",2*index  ,"S06B", xx, 0., kBframeWidth/2.-  kNulocWidth/2, 0, "ONLY");
	if (xx > -xxmax2 && xx< xxmax2) {
	  gMC->Gspos("S05E",2*index-1,"SB5B", xx, 0.,-kBframeWidth/2.+ kNulocWidth/2, 0, "ONLY");
	  gMC->Gspos("S05E",2*index  ,"SB5B", xx, 0., kBframeWidth/2.- kNulocWidth/2, 0, "ONLY");
	}
      }

      // position the volumes approximating the circular section of the pipe
      Float_t epsilon = 0.001; 
      Int_t ndiv = 6;
      Int_t imax = 1;
      Double_t divpar[3];
      Double_t dydiv = kSensHeight/ndiv;
      Double_t ydiv  = (kSensHeight - dydiv)/2.;
      Double_t rmin  = 31.5;  // Corrected in sep04 from PQ-LAT-SR2 de CEA-DSM-DAPNIA-SIS/BE ph HARDY 19-Oct-2002 slat 
      Double_t xdiv  = 0.;
      Float_t xvol;
      Float_t yvol;

      for (Int_t idiv = 0; idiv < ndiv; idiv++){ 
	ydiv += dydiv;
	xdiv = 0.; 
	if (ydiv < rmin) xdiv = rmin * TMath::Sin( TMath::ACos(ydiv/rmin) );
	divpar[0] = (kPcbLength - xdiv)/2.; 
	divpar[1] = dydiv/2. - epsilon;
	divpar[2] = kSensWidth/2.; 
	xvol = (kPcbLength + xdiv)/2.;
	yvol = ydiv; 

	// Volumes close to the beam pipe for slat i=1 so 4 slats per chamber
	for (Int_t quadrant = 1; quadrant <= 4; quadrant++) {
	  sprintf(idSlatCh5,"LA%d",ConvertSlatNum(1,quadrant,kNslats3-1));
	  sprintf(idSlatCh6,"LB%d",ConvertSlatNum(1,quadrant,kNslats3-1));

	  GetEnvelopes(4)->AddEnvelopeConstituentParam("S05G", idSlatCh5, quadrant*100+imax+4*idiv+1,
						       TGeoTranslation(xvol-(kPcbLength * kNPCB3[1]/2.),yvol-kPcbLength,0.),3,divpar);

	  GetEnvelopes(5)->AddEnvelopeConstituentParam("S06G", idSlatCh6,  quadrant*100+imax+4*idiv+1,
						       TGeoTranslation(xvol-(kPcbLength * kNPCB3[1]/2.),yvol-kPcbLength,0.),3,divpar);
	}
      }

      // Volumes close to the beam pipe for slat i=0 so 2 slats per chamber (central slat for station 3)
      //      Gines Martinez, Subatech sep 04
      // 9 box volumes are used to define the PCB closed to the beam pipe of the slat 122000SR1 of chamber 5 and 6 of St3
      // Accordingly to plan PQ-LAT-SR1 of CEA-DSM-DAPNIA-SIS/BE ph HARDY 8-Oct-2002
      // Rmin = 31.5 cm
      rmin = 31.5; //in cm  
      ndiv  = 9; 
      dydiv = kSensHeight/ndiv;           // Vertical size of the box volume approximating the rounded PCB
      ydiv  = -kSensHeight/2 + dydiv/2.;   // Initializing vertical position of the volume from bottom
      xdiv  = 0.;                         // Initializing horizontal position of the box volumes

      for (Int_t idiv = 0; idiv < ndiv; idiv++){ 
	xdiv = TMath::Abs( rmin * TMath::Sin( TMath::ACos(ydiv/rmin) ) );
	divpar[0] = (kPcbLength - xdiv)/2.; // Dimension of the box volume
	divpar[1] = dydiv/2. - epsilon;
	divpar[2] = kSensWidth/2.; 
	xvol = (kPcbLength + xdiv)/2.; //2D traslition for positionning of box volume
	yvol =  ydiv;
	Int_t side;
	for (side = 1; side <= 2; side++) {
	  sprintf(idSlatCh5,"LA%d",4); 	   
	  sprintf(idSlatCh6,"LB%d",4);
	  if(side == 2) {
	    sprintf(idSlatCh5,"LA%d",13); 	   
	    sprintf(idSlatCh6,"LB%d",13);
	  }	   
	  GetEnvelopes(4)->AddEnvelopeConstituentParam("S05G", idSlatCh5,500+side*100+imax+4*idiv+1,
						       TGeoTranslation(xvol-(kPcbLength * kNPCB3[0]/2.),yvol,0.),3,divpar);

	  GetEnvelopes(5)->AddEnvelopeConstituentParam("S06G", idSlatCh6,500+side*100+imax+4*idiv+1,
						       TGeoTranslation(xvol-(kPcbLength * kNPCB3[0]/2.),yvol,0.),3,divpar);
	}
	ydiv += dydiv; // Going from bottom to top
      }
      // cout << "Geometry for Station 3...... done" << endl;	
    }
    
  if (fStations[3]) {


    // //********************************************************************
    // //                            Station 4                             **
    // //********************************************************************
    //      // indices 1 and 2 for first and second chambers in the station
    //      // iChamber (first chamber) kept for other quanties than Z,
    //      // assumed to be the same in both chambers
    //      corrected geometry (JP. Cussonneau, Ch. Finck)
 
    iChamber = &fMUON->Chamber(6);
    iChamber1 = iChamber;
    iChamber2 = &fMUON->Chamber(7);

    const Int_t   kNslats4          = 7;  // number of slats per quadrant
    const Int_t   kNPCB4[kNslats4]  = {5, 6, 5, 5, 4, 3, 2}; // n PCB per slat
    const Float_t kXpos4[kNslats4]  = {38.2, 0., 0., 0., 0., 0., 0.};
    const Float_t kYpos41[kNslats4] = {0., 38.2, 34.40, 36.60, 29.3, 37.0, 28.6};
    const Float_t kYpos42[kNslats4] = {0., 38.2, 37.85, 37.55, 29.4, 37.0, 28.6};

    Float_t slatLength4[kNslats4];     

    // create and position the slat (mother) volumes 

    char idSlatCh7[5];
    char idSlatCh8[5];
    Float_t xSlat4;
    Float_t ySlat41 = 0;
    Float_t ySlat42 = 0;

    angle = 0.;

    for (i = 0; i<kNslats4; i++){
      slatLength4[i] = kPcbLength * kNPCB4[i] + 2. * kVframeLength; 
      xSlat4 = slatLength4[i]/2. + kDslatLength + kXpos4[i]; 
      ySlat41 += kYpos41[i];
      ySlat42 += kYpos42[i];

      spar[0] = slatLength4[i]/2.; 
      spar[1] = kSlatHeight/2.;
      spar[2] = kSlatWidth/2.; 
      Float_t dzCh4 = dzCh;
      Float_t zSlat4 = (i%2 ==0)? -zSlat : zSlat; 

      sprintf(idSlatCh7,"LC%d",kNslats4-1+i);
      gMC->Gsvolu(idSlatCh7,"BOX",kSlatMaterial,spar,3);
      detElemId = 700 + i + kNslats4-1;
      GetEnvelopes(6)->AddEnvelope(idSlatCh7, detElemId, true, TGeoTranslation(xSlat4, ySlat41, -zSlat4 + dzCh4),
				   TGeoRotation("rot1",90,angle,90,90+angle,0,0) );

      sprintf(idSlatCh7,"LC%d",3*kNslats4-2+i);
      gMC->Gsvolu(idSlatCh7,"BOX",kSlatMaterial,spar,3);
      detElemId = 750 + i + kNslats4-1;
      GetEnvelopes(6)->AddEnvelope(idSlatCh7, detElemId, true, TGeoTranslation(-xSlat4, ySlat41, zSlat4 - dzCh4),
				   TGeoRotation("rot2",90,180+angle,90,90+angle,180,0) );
 
      if (i > 0) { 
	sprintf(idSlatCh7,"LC%d",kNslats4-1-i);
	gMC->Gsvolu(idSlatCh7,"BOX",kSlatMaterial,spar,3);
	detElemId = 700 - i + kNslats4-1;
	GetEnvelopes(6)->AddEnvelope(idSlatCh7, detElemId, true, TGeoTranslation(xSlat4, -ySlat41, -zSlat4 + dzCh4),
				     TGeoRotation("rot3",90,angle,90,270+angle,180,0) );

	sprintf(idSlatCh7,"LC%d",3*kNslats4-2-i);
	detElemId = 750 - i + kNslats4-1;
	gMC->Gsvolu(idSlatCh7,"BOX",kSlatMaterial,spar,3);
	GetEnvelopes(6)->AddEnvelope(idSlatCh7, detElemId, true, 
				     TGeoTranslation(-xSlat4, -ySlat41, zSlat4 - dzCh4),
				     TGeoRotation("rot4",90,180+angle,90,270+angle,0,0) );
      }

      sprintf(idSlatCh8,"LD%d",kNslats4-1+i);
      gMC->Gsvolu(idSlatCh8,"BOX",kSlatMaterial,spar,3);
      detElemId = 800 + i + kNslats4-1;
      GetEnvelopes(7)->AddEnvelope(idSlatCh8, detElemId, true, TGeoTranslation(xSlat4, ySlat42, -zSlat4 + dzCh4),
				   TGeoRotation("rot5",90,angle,90,90+angle,0,0) );

      sprintf(idSlatCh8,"LD%d",3*kNslats4-2+i);
      detElemId = 850 + i + kNslats4-1;
      gMC->Gsvolu(idSlatCh8,"BOX",kSlatMaterial,spar,3);
      GetEnvelopes(7)->AddEnvelope(idSlatCh8, detElemId, true, TGeoTranslation(-xSlat4, ySlat42, zSlat4 - dzCh4),
				   TGeoRotation("rot6",90,180+angle,90,90+angle,180,0) );
      if (i > 0) { 
	sprintf(idSlatCh8,"LD%d",kNslats4-1-i);
	detElemId = 800 - i + kNslats4-1;
	gMC->Gsvolu(idSlatCh8,"BOX",kSlatMaterial,spar,3);
	GetEnvelopes(7)->AddEnvelope(idSlatCh8, detElemId, true, TGeoTranslation(xSlat4, -ySlat42, -zSlat4 + dzCh4),
				     TGeoRotation("rot7",90,angle,90,270+angle,180,0) );
	sprintf(idSlatCh8,"LD%d",3*kNslats4-2-i);
	detElemId = 850 - i + kNslats4-1;
	gMC->Gsvolu(idSlatCh8,"BOX",kSlatMaterial,spar,3);
	GetEnvelopes(7)->AddEnvelope(idSlatCh8, detElemId, true, TGeoTranslation(-xSlat4, -ySlat42, zSlat4 - dzCh4),
				     TGeoRotation("rot8",90,180+angle,90,270+angle,0,0) );
      }
    }
     
    // create the panel volume 
 
    gMC->Gsvolu("S07C","BOX",kCarbonMaterial,panelpar,3);
    gMC->Gsvolu("S08C","BOX",kCarbonMaterial,panelpar,3);

    // create the nomex volume 

    gMC->Gsvolu("S07N","BOX",kNomexMaterial,nomexpar,3);
    gMC->Gsvolu("S08N","BOX",kNomexMaterial,nomexpar,3);


    // create the nomex volume (bulk)

    gMC->Gsvolu("S07X","BOX",kNomexBMaterial,nomexbpar,3);
    gMC->Gsvolu("S08X","BOX",kNomexBMaterial,nomexbpar,3);

    // create the insulating material volume 

    gMC->Gsvolu("S07I","BOX",kInsuMaterial,insupar,3);
    gMC->Gsvolu("S08I","BOX",kInsuMaterial,insupar,3);

    // create the PCB volume 

    gMC->Gsvolu("S07P","BOX",kPcbMaterial,pcbpar,3);
    gMC->Gsvolu("S08P","BOX",kPcbMaterial,pcbpar,3);
 
    // create the sensitive volumes,

    gMC->Gsvolu("S07G","BOX",kSensMaterial,dum,0);
    gMC->Gsvolu("S08G","BOX",kSensMaterial,dum,0);

    // create the vertical frame volume 

    gMC->Gsvolu("S07V","BOX",kVframeMaterial,vFramepar,3);
    gMC->Gsvolu("S08V","BOX",kVframeMaterial,vFramepar,3);

    // create the horizontal frame volume 

    gMC->Gsvolu("S07H","BOX",kHframeMaterial,hFramepar,3);
    gMC->Gsvolu("S08H","BOX",kHframeMaterial,hFramepar,3);

    // create the horizontal border volume 

    gMC->Gsvolu("S07B","BOX",kBframeMaterial,bFramepar,3);
    gMC->Gsvolu("S08B","BOX",kBframeMaterial,bFramepar,3);

    index = 0; 
    for (i = 0; i < kNslats4; i++){
      for (Int_t quadrant = 1; quadrant <= 4; quadrant++) {

	if (i == 0 && quadrant == 2) continue;
	if (i == 0 && quadrant == 4) continue;

	sprintf(idSlatCh7,"LC%d",ConvertSlatNum(i,quadrant,kNslats4-1));
	sprintf(idSlatCh8,"LD%d",ConvertSlatNum(i,quadrant,kNslats4-1));
	Float_t xvFrame  = (slatLength4[i] - kVframeLength)/2.;

	// position the vertical frames 
	if (i != 1) { 
	  GetEnvelopes(6)->AddEnvelopeConstituent("S07V", idSlatCh7, (2*i-1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
	  GetEnvelopes(6)->AddEnvelopeConstituent("S07V", idSlatCh7, (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));
	  GetEnvelopes(7)->AddEnvelopeConstituent("S08V", idSlatCh8, (2*i-1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
	  GetEnvelopes(7)->AddEnvelopeConstituent("S08V", idSlatCh8, (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));
	} else { // no rounded spacer yet
	  GetEnvelopes(6)->AddEnvelopeConstituent("S07V", idSlatCh7, (2*i-1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
	  // GetEnvelopes(6)->AddEnvelopeConstituent("S07V", idSlatCh7, (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));
	  GetEnvelopes(7)->AddEnvelopeConstituent("S08V", idSlatCh8, (2*i-1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
	  // GetEnvelopes(7)->AddEnvelopeConstituent("S08V", idSlatCh8, (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));
	}
	// position the panels and the insulating material 
	for (j = 0; j < kNPCB4[i]; j++){
	  if (i == 1 && j == 0) continue;
	  index++;
	  Float_t xx = kSensLength * (-kNPCB4[i]/2.+j+.5); 

	  Float_t zPanel = spar[2] - nomexbpar[2]; 
	  GetEnvelopes(6)->AddEnvelopeConstituent("S07X", idSlatCh7, 2*index-1,TGeoTranslation(xx,0.,zPanel));
	  GetEnvelopes(6)->AddEnvelopeConstituent("S07X", idSlatCh7, 2*index,TGeoTranslation(xx,0.,-zPanel));
	  GetEnvelopes(6)->AddEnvelopeConstituent("S07I", idSlatCh7, index,TGeoTranslation(xx,0.,0.));
	  GetEnvelopes(7)->AddEnvelopeConstituent("S08X", idSlatCh8, 2*index-1,TGeoTranslation(xx,0.,zPanel));
	  GetEnvelopes(7)->AddEnvelopeConstituent("S08X", idSlatCh8, 2*index,TGeoTranslation(xx,0.,-zPanel));
	  GetEnvelopes(7)->AddEnvelopeConstituent("S08I", idSlatCh8, index,TGeoTranslation(xx,0.,0.));
	}
      } 
    }

    // position the nomex volume inside the panel volume
    gMC->Gspos("S07N",1,"S07C",0.,0.,0.,0,"ONLY"); 
    gMC->Gspos("S08N",1,"S08C",0.,0.,0.,0,"ONLY"); 

    // position panel volume inside the bulk nomex material volume
    gMC->Gspos("S07C",1,"S07X",0.,0.,kNomexBWidth/2.,0,"ONLY"); 
    gMC->Gspos("S08C",1,"S08X",0.,0.,kNomexBWidth/2.,0,"ONLY"); 

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
    Float_t yborder = ( kPcbHeight - kBframeHeight ) / 2.; 
    gMC->Gspos("S07B",1,"S07P",0., yborder,0.,0,"ONLY"); 
    gMC->Gspos("S07B",2,"S07P",0.,-yborder,0.,0,"ONLY"); 
    gMC->Gspos("S08B",1,"S08P",0., yborder,0.,0,"ONLY"); 
    gMC->Gspos("S08B",2,"S08P",0.,-yborder,0.,0,"ONLY"); 

    // create the NULOC volume and position it in the horizontal frame

    gMC->Gsvolu("S07E","BOX",kNulocMaterial,nulocpar,3);
    gMC->Gsvolu("S08E","BOX",kNulocMaterial,nulocpar,3);
    index = 0;
    for (xx = -xxmax; xx <= xxmax; xx += 2*kNulocLength) { 
      index++; 
      gMC->Gspos("S07E",2*index-1,"S07B", xx, 0.,-kBframeWidth/2. + kNulocWidth/2, 0, "ONLY");
      gMC->Gspos("S07E",2*index  ,"S07B", xx, 0., kBframeWidth/2. - kNulocWidth/2, 0, "ONLY");
      gMC->Gspos("S08E",2*index-1,"S08B", xx, 0.,-kBframeWidth/2. + kNulocWidth/2, 0, "ONLY");
      gMC->Gspos("S08E",2*index  ,"S08B", xx, 0., kBframeWidth/2. - kNulocWidth/2, 0, "ONLY");
    }

    // position the volumes approximating the circular section of the pipe

    Float_t epsilon = 0.001; 
    Int_t ndiv = 10;
    Int_t imax = 1; 
    Double_t divpar[3];
    Double_t dydiv = kSensHeight/ndiv;
    Double_t ydiv  = (kSensHeight - dydiv)/2.;
    Float_t rmin   = 39.5;// Corrected in sep04 from PQ-LAT-NR3 de CEA-DSM-DAPNIA-SIS/BE ph HARDY 19-Oct-2002 slat 
    Float_t xdiv   = 0.; 
    Float_t xvol;
    Float_t yvol;

    for (Int_t idiv = 0; idiv < ndiv; idiv++){ 
      ydiv += dydiv;
      xdiv = 0.; 
      if (ydiv < rmin) xdiv = rmin * TMath::Sin( TMath::ACos(ydiv/rmin) );
      divpar[0] = (kPcbLength - xdiv)/2.; 
      divpar[1] = dydiv/2. - epsilon;
      divpar[2] = kSensWidth/2.; 
      xvol = (kPcbLength + xdiv)/2.;
      yvol = ydiv ;
       
      for (Int_t quadrant = 1; quadrant <= 4; quadrant++) {
	sprintf(idSlatCh7,"LC%d",ConvertSlatNum(1,quadrant,kNslats4-1));
	sprintf(idSlatCh8,"LD%d",ConvertSlatNum(1,quadrant,kNslats4-1));
	 
	GetEnvelopes(6)->AddEnvelopeConstituentParam("S07G",idSlatCh7, quadrant*100+imax+4*idiv+1,
						     TGeoTranslation(xvol-kPcbLength * kNPCB4[1]/2.,yvol-kPcbLength,0.),3,divpar);
	 
	GetEnvelopes(7)->AddEnvelopeConstituentParam("S08G", idSlatCh8, quadrant*100+imax+4*idiv+1,
						     TGeoTranslation(xvol-kPcbLength * kNPCB4[1]/2.,yvol-kPcbLength,0.),3,divpar);
      }
    }
    // cout << "Geometry for Station 4...... done" << endl;

  }
    
  if (fStations[4]) {
      

    // //********************************************************************
    // //                            Station 5                             **
    // //********************************************************************
    //      // indices 1 and 2 for first and second chambers in the station
    //      // iChamber (first chamber) kept for other quanties than Z,
    //      // assumed to be the same in both chambers
    //      corrected geometry (JP. Cussonneau, Ch. Finck)

    iChamber = &fMUON->Chamber(8);
    iChamber1 = iChamber;
    iChamber2 = &fMUON->Chamber(9);
 
    const Int_t   kNslats5         = 7;  // number of slats per quadrant
    const Int_t   kNPCB5[kNslats5] = {5, 6, 6, 6, 5, 4, 3}; // n PCB per slat
    const Float_t kXpos5[kNslats5] = {38.2, 0., 0., 0., 0., 0., 0.};
    const Float_t kYpos5[kNslats5] = {0., 38.2, 37.9, 37.6, 37.3, 37.05, 36.75};
    Float_t slatLength5[kNslats5]; 

    // create and position the slat (mother) volumes 

    char idSlatCh9[5];
    char idSlatCh10[5];
    Float_t xSlat5;
    Float_t ySlat5 = 0;
    angle = 0.;

    for (i = 0; i < kNslats5; i++){

      slatLength5[i] = kPcbLength * kNPCB5[i] + 2.* kVframeLength; 
      xSlat5 = slatLength5[i]/2. + kDslatLength + kXpos5[i]; 
      ySlat5 += kYpos5[i];

      spar[0] = slatLength5[i]/2.; 
      spar[1] = kSlatHeight/2.;
      spar[2] = kSlatWidth/2.; 

      Float_t dzCh5  = dzCh;
      Float_t zSlat5 = (i%2 ==0)? -zSlat : zSlat; 

      sprintf(idSlatCh9,"LE%d",kNslats5-1+i);
      detElemId = 900 + i + kNslats5-1;
      gMC->Gsvolu(idSlatCh9,"BOX",kSlatMaterial,spar,3);
      GetEnvelopes(8)->AddEnvelope(idSlatCh9, detElemId, true, TGeoTranslation(xSlat5, ySlat5, -zSlat5 + dzCh5),
				   TGeoRotation("rot1",90,angle,90,90+angle,0,0) );

      sprintf(idSlatCh9,"LE%d",3*kNslats5-2+i);
      detElemId = 950 + i + kNslats5-1;
      gMC->Gsvolu(idSlatCh9,"BOX",kSlatMaterial,spar,3);
      GetEnvelopes(8)->AddEnvelope(idSlatCh9, detElemId, true, TGeoTranslation(-xSlat5, ySlat5, zSlat5 - dzCh5),
				   TGeoRotation("rot2",90,180+angle,90,90+angle,180,0) );
 
      if (i > 0) { 
	sprintf(idSlatCh9,"LE%d",kNslats5-1-i);
	detElemId = 900 - i + kNslats5-1;
	gMC->Gsvolu(idSlatCh9,"BOX",kSlatMaterial,spar,3);
	GetEnvelopes(8)->AddEnvelope(idSlatCh9, detElemId, true, TGeoTranslation(xSlat5, -ySlat5, -zSlat5 + dzCh5),
				     TGeoRotation("rot3",90,angle,90,270+angle,180,0) );

	sprintf(idSlatCh9,"LE%d",3*kNslats5-2-i);
	detElemId = 950 - i + kNslats5-1;
	gMC->Gsvolu(idSlatCh9,"BOX",kSlatMaterial,spar,3);
	GetEnvelopes(8)->AddEnvelope(idSlatCh9, detElemId, true, TGeoTranslation(-xSlat5, -ySlat5, zSlat5 - dzCh5),
				     TGeoRotation("rot4",90,180+angle,90,270+angle,0,0)  );
      }

      sprintf(idSlatCh10,"LF%d",kNslats5-1+i);
      detElemId = 1000 + i + kNslats5-1;
      gMC->Gsvolu(idSlatCh10,"BOX",kSlatMaterial,spar,3);
      GetEnvelopes(9)->AddEnvelope(idSlatCh10, detElemId, true, TGeoTranslation(xSlat5, ySlat5, -zSlat5 + dzCh5),
				   TGeoRotation("rot5",90,angle,90,90+angle,0,0) );

      sprintf(idSlatCh10,"LF%d",3*kNslats5-2+i);
      detElemId = 1050 + i + kNslats5-1;
      gMC->Gsvolu(idSlatCh10,"BOX",kSlatMaterial,spar,3);
      GetEnvelopes(9)->AddEnvelope(idSlatCh10, detElemId, true, TGeoTranslation(-xSlat5, ySlat5, zSlat5 - dzCh5),
				   TGeoRotation("rot6",90,180+angle,90,90+angle,180,0) );

      if (i > 0) { 
	sprintf(idSlatCh10,"LF%d",kNslats5-1-i);
	detElemId = 1000 - i + kNslats5-1;
	gMC->Gsvolu(idSlatCh10,"BOX",kSlatMaterial,spar,3);
	GetEnvelopes(9)->AddEnvelope(idSlatCh10, detElemId, true, TGeoTranslation(xSlat5, -ySlat5, -zSlat5 + dzCh5),
				     TGeoRotation("rot7",90,angle,90,270+angle,180,0) );
	sprintf(idSlatCh10,"LF%d",3*kNslats5-2-i);
	detElemId = 1050 - i + kNslats5-1;
	gMC->Gsvolu(idSlatCh10,"BOX",kSlatMaterial,spar,3);
	GetEnvelopes(9)->AddEnvelope(idSlatCh10, detElemId, true, TGeoTranslation(-xSlat5, -ySlat5, zSlat5 - dzCh5),
				     TGeoRotation("rot8",90,180+angle,90,270+angle,0,0) );
      }
    }

    // create the panel volume 
 
    gMC->Gsvolu("S09C","BOX",kCarbonMaterial,panelpar,3);
    gMC->Gsvolu("S10C","BOX",kCarbonMaterial,panelpar,3);

    // create the nomex volume 

    gMC->Gsvolu("S09N","BOX",kNomexMaterial,nomexpar,3);
    gMC->Gsvolu("S10N","BOX",kNomexMaterial,nomexpar,3);


    // create the nomex volume (bulk)

    gMC->Gsvolu("S09X","BOX",kNomexBMaterial,nomexbpar,3);
    gMC->Gsvolu("S10X","BOX",kNomexBMaterial,nomexbpar,3);

    // create the insulating material volume 

    gMC->Gsvolu("S09I","BOX",kInsuMaterial,insupar,3);
    gMC->Gsvolu("S10I","BOX",kInsuMaterial,insupar,3);

    // create the PCB volume 

    gMC->Gsvolu("S09P","BOX",kPcbMaterial,pcbpar,3);
    gMC->Gsvolu("S10P","BOX",kPcbMaterial,pcbpar,3);
 
    // create the sensitive volumes,

    gMC->Gsvolu("S09G","BOX",kSensMaterial,dum,0);
    gMC->Gsvolu("S10G","BOX",kSensMaterial,dum,0);

    // create the vertical frame volume 

    gMC->Gsvolu("S09V","BOX",kVframeMaterial,vFramepar,3);
    gMC->Gsvolu("S10V","BOX",kVframeMaterial,vFramepar,3);

    // create the horizontal frame volume 

    gMC->Gsvolu("S09H","BOX",kHframeMaterial,hFramepar,3);
    gMC->Gsvolu("S10H","BOX",kHframeMaterial,hFramepar,3);

    // create the horizontal border volume 

    gMC->Gsvolu("S09B","BOX",kBframeMaterial,bFramepar,3);
    gMC->Gsvolu("S10B","BOX",kBframeMaterial,bFramepar,3);

    index = 0; 
    for (i = 0; i < kNslats5; i++){
      for (Int_t quadrant = 1; quadrant <= 4; quadrant++) {

	if (i == 0 && quadrant == 2) continue;
	if (i == 0 && quadrant == 4) continue;

	sprintf(idSlatCh9,"LE%d",ConvertSlatNum(i,quadrant,kNslats5-1));
	sprintf(idSlatCh10,"LF%d",ConvertSlatNum(i,quadrant,kNslats5-1));
	Float_t xvFrame  = (slatLength5[i] - kVframeLength)/2.; // ok

	// position the vertical frames (spacers)
	if (i != 1) { 
	  GetEnvelopes(8)->AddEnvelopeConstituent("S09V", idSlatCh9, (2*i-1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
	  GetEnvelopes(8)->AddEnvelopeConstituent("S09V", idSlatCh9, (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));
	  GetEnvelopes(9)->AddEnvelopeConstituent("S10V", idSlatCh10, (2*i-1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
	  GetEnvelopes(9)->AddEnvelopeConstituent("S10V", idSlatCh10, (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));
	} else {  // no rounded spacer yet
	  GetEnvelopes(8)->AddEnvelopeConstituent("S09V", idSlatCh9, (2*i-1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
	  //	   GetEnvelopes(8)->AddEnvelopeConstituent("S09V", idSlatCh9, (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));
	  GetEnvelopes(9)->AddEnvelopeConstituent("S10V", idSlatCh10, (2*i-1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
	  //	   GetEnvelopes(9)->AddEnvelopeConstituent("S10V", idSlatCh10, (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));
	}

	// position the panels and the insulating material 
	for (j = 0; j < kNPCB5[i]; j++){
	  if (i == 1 && j == 0) continue;
	  index++;
	  Float_t xx = kSensLength * (-kNPCB5[i]/2.+j+.5); 

	  Float_t zPanel = spar[2] - nomexbpar[2]; 
	  GetEnvelopes(8)->AddEnvelopeConstituent("S09X", idSlatCh9, 2*index-1,TGeoTranslation(xx,0.,zPanel));
	  GetEnvelopes(8)->AddEnvelopeConstituent("S09X", idSlatCh9, 2*index,TGeoTranslation(xx,0.,-zPanel));
	  GetEnvelopes(8)->AddEnvelopeConstituent("S09I", idSlatCh9, index,TGeoTranslation(xx,0.,0.));

	  GetEnvelopes(9)->AddEnvelopeConstituent("S10X", idSlatCh10, 2*index-1,TGeoTranslation(xx,0.,zPanel));
	  GetEnvelopes(9)->AddEnvelopeConstituent("S10X", idSlatCh10, 2*index,TGeoTranslation(xx,0.,-zPanel));
	  GetEnvelopes(9)->AddEnvelopeConstituent("S10I", idSlatCh10, index,TGeoTranslation(xx,0.,0.));
	}
      } 
    }

    // position the nomex volume inside the panel volume
    gMC->Gspos("S09N",1,"S09C",0.,0.,0.,0,"ONLY"); 
    gMC->Gspos("S10N",1,"S10C",0.,0.,0.,0,"ONLY"); 

    // position panel  volume inside the bulk nomex material volume
    gMC->Gspos("S09C",1,"S09X",0.,0.,kNomexBWidth/2.,0,"ONLY"); 
    gMC->Gspos("S10C",1,"S10X",0.,0.,kNomexBWidth/2.,0,"ONLY"); 

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
    Float_t yborder = ( kPcbHeight - kBframeHeight ) / 2.; 
    gMC->Gspos("S09B",1,"S09P",0., yborder,0.,0,"ONLY"); 
    gMC->Gspos("S09B",2,"S09P",0.,-yborder,0.,0,"ONLY"); 
    gMC->Gspos("S10B",1,"S10P",0., yborder,0.,0,"ONLY"); 
    gMC->Gspos("S10B",2,"S10P",0.,-yborder,0.,0,"ONLY"); 

    //      // create the NULOC volume and position it in the horizontal frame

    gMC->Gsvolu("S09E","BOX",kNulocMaterial,nulocpar,3);
    gMC->Gsvolu("S10E","BOX",kNulocMaterial,nulocpar,3);
    index = 0;
    for (xx = -xxmax; xx <= xxmax; xx += 2*kNulocLength) { 
      index++; 
      gMC->Gspos("S09E",2*index-1,"S09B", xx, 0.,-kBframeWidth/2. + kNulocWidth/2, 0, "ONLY");
      gMC->Gspos("S09E",2*index  ,"S09B", xx, 0., kBframeWidth/2. - kNulocWidth/2, 0, "ONLY");
      gMC->Gspos("S10E",2*index-1,"S10B", xx, 0.,-kBframeWidth/2. + kNulocWidth/2, 0, "ONLY");
      gMC->Gspos("S10E",2*index  ,"S10B", xx, 0., kBframeWidth/2. - kNulocWidth/2, 0, "ONLY");
    }


    // position the volumes approximating the circular section of the pipe
    Float_t epsilon = 0.001; 
    Int_t ndiv = 10;
    Int_t imax = 1; 
    Double_t divpar[3];
    Double_t dydiv = kSensHeight/ndiv;
    Double_t ydiv  = (kSensHeight - dydiv)/2.;
    Float_t rmin   = 39.5;
    Float_t xdiv   = 0.; 
    Float_t xvol;
    Float_t yvol; 

    for (Int_t idiv = 0; idiv < ndiv; idiv++){ 
      ydiv += dydiv;
      xdiv = 0.; 
      if (ydiv < rmin) xdiv = rmin * TMath::Sin( TMath::ACos(ydiv/rmin) );
      divpar[0] = (kPcbLength - xdiv)/2.; 
      divpar[1] = dydiv/2. - epsilon;
      divpar[2] = kSensWidth/2.; 
      xvol = (kPcbLength + xdiv)/2.;
      yvol = ydiv;

      for (Int_t quadrant = 1; quadrant <= 4; quadrant++) {
	sprintf(idSlatCh9,"LE%d",ConvertSlatNum(1,quadrant,kNslats5-1));
	sprintf(idSlatCh10,"LF%d",ConvertSlatNum(1,quadrant,kNslats5-1));

	GetEnvelopes(8)->AddEnvelopeConstituentParam("S09G", idSlatCh9, quadrant*100+imax+4*idiv+1,
						     TGeoTranslation(xvol-kPcbLength * kNPCB5[1]/2.,yvol-kPcbLength,0.),3,divpar);
	GetEnvelopes(9)->AddEnvelopeConstituentParam("S10G", idSlatCh10,  quadrant*100+imax+4*idiv+1,
						     TGeoTranslation(xvol-kPcbLength * kNPCB5[1]/2.,yvol-kPcbLength,0.),3,divpar);
      }
    }
    // cout << "Geometry for Station 5...... done" << endl;

  }
}


//______________________________________________________________________________
void AliMUONSlatGeometryBuilder::SetTransformations()
{
// Defines the transformations for the station2 chambers.
// ---

  AliMUONChamber* iChamber1 = &fMUON->Chamber(4);
  Double_t zpos1 = - iChamber1->Z(); 
  iChamber1->GetGeometry()
    ->SetTranslation(TGeoTranslation(0., 0., zpos1));

  AliMUONChamber* iChamber2 = &fMUON->Chamber(5);
  Double_t zpos2 = - iChamber2->Z(); 
  iChamber2->GetGeometry()
    ->SetTranslation(TGeoTranslation(0., 0., zpos2));

 iChamber1 = &fMUON->Chamber(6);
  zpos1 = - iChamber1->Z(); 
  iChamber1->GetGeometry()
    ->SetTranslation(TGeoTranslation(0., 0., zpos1));

  iChamber2 = &fMUON->Chamber(7);
  zpos2 = - iChamber2->Z(); 
  iChamber2->GetGeometry()
    ->SetTranslation(TGeoTranslation(0., 0., zpos2));

 iChamber1 = &fMUON->Chamber(8);
  zpos1 = - iChamber1->Z(); 
  iChamber1->GetGeometry()
    ->SetTranslation(TGeoTranslation(0., 0., zpos1));

  iChamber2 = &fMUON->Chamber(9);
  zpos2 = - iChamber2->Z(); 
  iChamber2->GetGeometry()
    ->SetTranslation(TGeoTranslation(0., 0., zpos2));

}

//______________________________________________________________________________
void AliMUONSlatGeometryBuilder::SetSensitiveVolumes()
{
// Defines the sensitive volumes for slat stations chambers.
// ---

  GetGeometry(4)->SetSensitiveVolume("S05G");
  GetGeometry(5)->SetSensitiveVolume("S06G");
  GetGeometry(6)->SetSensitiveVolume("S07G");
  GetGeometry(7)->SetSensitiveVolume("S08G");
  GetGeometry(8)->SetSensitiveVolume("S09G");
  GetGeometry(9)->SetSensitiveVolume("S10G");
}

//______________________________________________________________________________
Int_t  AliMUONSlatGeometryBuilder::ConvertSlatNum(Int_t numslat, Int_t quadnum, Int_t fspq) const
{
// On-line function establishing the correspondance between numslat (the slat number on a particular quadrant (numslat->0....4 for St3))
// and slatnum (the slat number on the whole panel (slatnum->1...18 for St3)
  numslat += 1;
  if (quadnum==2 || quadnum==3) 
    numslat += fspq;
  else
    numslat = fspq + 2-numslat;
  numslat -= 1;
	      
  if (quadnum==3 || quadnum==4) numslat += 2*fspq+1;

  return numslat;
}
