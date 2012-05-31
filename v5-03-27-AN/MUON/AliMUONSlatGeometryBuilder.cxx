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

//-----------------------------------------------------------------------------
/// \class AliMUONSlatGeometryBuilder
/// This Builder is designed according to the enveloppe methode. The basic idea is to be able to allow moves 
/// of the slats on the support panels. 
/// Those moves can be described with a simple set of parameters. The next step should be now to describe all 
/// the slats and their places by a unique 
/// class, which would make the SlatBuilder far more compact since now only three parameters can define a slat 
/// and its position, like:
///   - Bool_t rounded_shape_slat
///   - Float_t slat_length
///   - Float_t slat_number or Float_t slat_position
/// Reference system is the one described in the note ALICE-INT-2003-038  v.2  EDMS Id 406391 
///
/// \author Eric Dumonteil (dumontei@cea.fr)
//-----------------------------------------------------------------------------

#include "AliMUONSlatGeometryBuilder.h"
#include "AliMUON.h"
#include "AliMUONConstants.h"
#include "AliMUONGeometryModule.h"
#include "AliMUONGeometryEnvelopeStore.h"
#include "AliMUONConstants.h"

#include "AliMpDEManager.h"

#include "AliRun.h"
#include "AliLog.h"

#include <TVirtualMC.h>
#include <TGeoBBox.h>
#include <TGeoVolume.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoCompositeShape.h>
#include <TGeoTube.h>
#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONSlatGeometryBuilder)
/// \endcond

//______________________________________________________________________________
AliMUONSlatGeometryBuilder::AliMUONSlatGeometryBuilder(AliMUON* muon)
 : AliMUONVGeometryBuilder(4, 12),
   fMUON(muon)
{
/// Standard constructor

}

//______________________________________________________________________________
AliMUONSlatGeometryBuilder::AliMUONSlatGeometryBuilder() 
 : AliMUONVGeometryBuilder(),
   fMUON(0)
{
/// Default constructor
}

//______________________________________________________________________________
AliMUONSlatGeometryBuilder::~AliMUONSlatGeometryBuilder() 
{
/// Destructor
}

//
// public methods
//

//______________________________________________________________________________
void AliMUONSlatGeometryBuilder::CreateGeometry()
{
  /// CreateGeometry is the method containing all the informations concerning Stations 345 geometry.
  /// It includes description and placements of support panels and slats.
  /// The code comes directly from what was written in AliMUONv1.cxx before, with modifications concerning 
  /// the use of Enveloppe method to place the Geant volumes.
  /// Now, few changes would allow the creation of a Slat methode where slat could be described by few parameters, 
  /// and this builder would then be dedicated only to the
  /// placements of the slats. Those modifications could shorten the Station 345 geometry by a non-negligeable factor...
 
  Int_t *idtmed = fMUON->GetIdtmed()->GetArray()-1099;

  Float_t angle;
  Float_t *dum=0;

  // define the id of tracking media:
  //  Int_t idAir    = idtmed[1100]; // medium 1
  Int_t idGas    = idtmed[1108]; // medium 9 = Ar-CO2 gas (80%+20%)
  Int_t idCopper = idtmed[1110];
  Int_t idG10    = idtmed[1111];
  Int_t idCarbon = idtmed[1112];
  Int_t idRoha   = idtmed[1113];
  Int_t idNomex  = idtmed[1114]; // honey comb
  Int_t idNoryl  = idtmed[1115]; 
  Int_t idNomexB = idtmed[1116]; // bulk material 
  
  // Getting mediums for pannel support geometry
  TGeoMedium* kMedNomex     = gGeoManager->GetMedium("MUON_Nomex");
  TGeoMedium* kMedCarbon    = gGeoManager->GetMedium("MUON_CARBON");

  // sensitive area: 40*40 cm**2
  const Float_t kSensLength = 40.; 
  const Float_t kSensHeight = 40.; 
  const Float_t kSensWidth  = AliMUONConstants::Pitch()*2;// 0.5 cm, according to TDR fig 2.120 
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

  // Frame along the rounded (spacers) slats 
  const Float_t kRframeHeight = 2.00; 

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

  // R: rounded part of vertical spacers 
  const Float_t kRframeLength = 2.0; 
  const Float_t kRframeWidth  = kSensWidth;
  const Int_t kRframeMaterial = idNoryl;

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

  // Readout cables: Equivalent to 260 mum copper   
  const Float_t kCableHeight = 2.6; 
  const Float_t kCableWidth  = 0.026;  
  const Int_t kCableMaterial = idCopper;

  // Slat parameters
  const Float_t kSlatHeight = kPcbHeight; 
  const Float_t kSlatWidth  = kSensWidth + 2.*(kPcbWidth + kInsuWidth + kPanelWidth 
					       + kNomexBWidth); //replaced rohacell with Nomex Ch. Finck 
  // const Int_t   kSlatMaterial = idAir;
  const Float_t kDslatLength  = -1.25; // position of the slat respect to the beam plane (half vertical spacer) Ch. Finck
  Float_t zSlat               = AliMUONConstants::DzSlat();// implemented Ch. Finck
  Float_t dzCh                = AliMUONConstants::DzCh();

  Float_t spar[3];  
  Int_t i, j;
  Int_t detElemId;
  Int_t moduleId;

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
  Int_t* fStations = new Int_t[5];
  for (Int_t iSt=0; iSt<5; iSt++) fStations[iSt] = 1;
  fStations[2] = 1;
     
  if (fStations[2])
    {
      //********************************************************************
      //                            Station 3                             **
      //********************************************************************
      // Mother volume for each chamber in St3 is an envelop (or assembly)
      // There is one assembly mother per half a chamber
      // Mother volume for each chamber in St3 is an envelop (or assembly)
      // There is one assembly mother per half a chamber  called SC05I, SC05O, SC06I and SC06O
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
      // Only for chamber 5: slat 1 has a PCB shorter by 5cm!

      Float_t tlength = 35.;
      Float_t panelpar2[3]  = { tlength/2., panelpar[1],  panelpar[2]}; 
      Float_t nomexpar2[3]  = { tlength/2., nomexpar[1],  nomexpar[2]}; 
      Float_t nomexbpar2[3] = { tlength/2., nomexbpar[1],  nomexbpar[2]}; 
      Float_t insupar2[3]   = { tlength/2., insupar[1],   insupar[2]}; 
      Float_t pcbpar2[3]    = { tlength/2., pcbpar[1],    pcbpar[2]}; 
      Float_t senspar2[3]   = { tlength/2., senspar[1],   senspar[2]}; 
      Float_t hFramepar2[3] = { tlength/2., hFramepar[1], hFramepar[2]}; 
      Float_t bFramepar2[3] = { tlength/2., bFramepar[1], bFramepar[2]}; 
      Float_t pcbDLength3   = (kPcbLength - tlength);
      
      // For rounded pcb of central slat
      Float_t csvPcbLength = 59.25-40.; // PQ-LAT-SR1
      Float_t panelpar3[3]  = { csvPcbLength/2., panelpar[1],  panelpar[2]}; 
      Float_t nomexpar3[3]  = { csvPcbLength/2., nomexpar[1],  nomexpar[2]}; 
      Float_t nomexbpar3[3] = { csvPcbLength/2., nomexbpar[1],  nomexbpar[2]}; 
      Float_t insupar3[3]   = { csvPcbLength/2., insupar[1],   insupar[2]}; 
      Float_t pcbpar3[3]    = { csvPcbLength/2., pcbpar[1],    pcbpar[2]}; 
      Float_t senspar3[3]   = { csvPcbLength/2., senspar[1],   senspar[2]}; 
      Float_t hFramepar3[3] = { csvPcbLength/2., hFramepar[1], hFramepar[2]}; 
      Float_t bFramepar3[3] = { csvPcbLength/2., bFramepar[1], bFramepar[2]}; 
      Float_t cPhi = TMath::RadToDeg()*(TMath::Pi()/2.-TMath::ACos(hFramepar3[1]/(AliMUONConstants::Rmin(2)-kRframeLength)));
      Float_t cFramepar3[5] = { AliMUONConstants::Rmin(2)-kRframeLength, AliMUONConstants::Rmin(2), kRframeWidth, -cPhi, cPhi}; 

      const Int_t   kNslats3         = 5;  // number of slats per quadrant
      const Int_t   kNPCB3[kNslats3] = {4, 4, 4, 3, 2}; // n PCB per slat
      const Float_t kXpos3[kNslats3] = {0., 0., 0., 0., 0.};//{31., 0., 0., 0., 0.};
      const Float_t kYpos3[kNslats3] = {0, 37.8, 37.7, 37.3, 33.7};
      Float_t slatLength3[kNslats3]; 

      Float_t rPhi1 = TMath::RadToDeg()*(TMath::ASin((kYpos3[1]-hFramepar3[1])/(AliMUONConstants::Rmin(2))));
      Float_t rPhi2 = TMath::RadToDeg()*(TMath::ACos(-vFramepar[0]/(AliMUONConstants::Rmin(2)-kRframeLength)));
      Float_t rFramepar3[5] = { AliMUONConstants::Rmin(2)-kRframeLength, AliMUONConstants::Rmin(2), kRframeWidth, rPhi1, rPhi2}; 
      Float_t vrFrameHeight = hFramepar3[1]+kYpos3[1]-AliMUONConstants::Rmin(2)+kRframeLength;

      // create and position the slat (mother) volumes 

      char idSlatCh5[6];
      char idSlatCh6[6];
      Float_t xSlat3;
      Float_t ySlat3 = 0;
      angle = 0.;

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
	Float_t dzSlat3 = -0.25; // see drawing PQ7EN345-6 (Delta_slat=80mm instead 85mm)
	Float_t zSlat3 = (i%2 ==0)? -(zSlat+dzSlat3) : (zSlat+dzSlat3); // seems not that zSlat3 = zSlat4 & 5 refering to plan PQ7EN345-6 ? -> Indeed, fixed J.C.

	sprintf(idSlatCh5,"SLA%d",i+kNslats3-1);
	detElemId = 509 - (i + kNslats3-1-4);
	moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
	if (detElemId == 508 || detElemId == 509) // Round slat, new rotation due to mapping convention
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh5, detElemId, true, TGeoTranslation(xSlat3, ySlat3, -zSlat3 + dzCh3),
					      TGeoRotation("rot1",90,180+angle,90,90+angle,180,0) );
	else {
	  if (detElemId % 2 == 0)
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh5, detElemId, true, TGeoTranslation(xSlat3, ySlat3, -zSlat3 + dzCh3),
						TGeoRotation("rot1",90,angle,90,90+angle,0,0) );	  
	  else 
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh5, detElemId, true, TGeoTranslation(xSlat3, ySlat3, -zSlat3 + dzCh3),
						TGeoRotation("rot1",90,angle,90,270+angle,180,0) ); 
	}
     
	sprintf(idSlatCh5,"SLA%d",3*kNslats3-2+i);
	detElemId = 500 + (i + kNslats3-1-4);
	moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
	if (detElemId == 500 || detElemId == 501) // Round slat, new rotation due to mapping convention
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh5, detElemId, true, TGeoTranslation(-xSlat3, ySlat3, zSlat3 - dzCh3),
					      TGeoRotation("rot2",90,angle,90,90+angle,0,0) );
	else {
	  if (detElemId % 2 == 1) 
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh5, detElemId, true, TGeoTranslation(-xSlat3, ySlat3, zSlat3 - dzCh3),
						TGeoRotation("rot2",90,180+angle,90,90+angle,180,0) );	   
	  else
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh5, detElemId, true, TGeoTranslation(-xSlat3, ySlat3, zSlat3 - dzCh3),
						TGeoRotation("rot2",90,180+angle,90,270+angle,0,0) );
	}

	if (i > 0) { 
	  sprintf(idSlatCh5,"SLA%d",kNslats3-1-i);
	  detElemId = 509 + (i + kNslats3-1-4);
  	  moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
	  if (detElemId % 2 == 0 ) {
	    if (detElemId == 510) // Round slat, new rotation due to mapping convention
	      GetEnvelopes(moduleId)->AddEnvelope(idSlatCh5, detElemId, true, TGeoTranslation(xSlat3, -ySlat3, -zSlat3 + dzCh3), 
						  TGeoRotation("rot3",90,180+angle,90,270+angle,0,0) );
	    else
	      GetEnvelopes(moduleId)->AddEnvelope(idSlatCh5, detElemId, true, TGeoTranslation(xSlat3, -ySlat3, -zSlat3 + dzCh3), 
						  TGeoRotation("rot3",90,angle,90,90+angle,0,0) );
	  }
	  else
	      GetEnvelopes(moduleId)->AddEnvelope(idSlatCh5, detElemId, true, TGeoTranslation(xSlat3, -ySlat3, -zSlat3 + dzCh3), 
						  TGeoRotation("rot3",90,angle,90,270+angle,180,0) );

	  sprintf(idSlatCh5,"SLA%d",3*kNslats3-2-i);
	  detElemId = 518 - (i + kNslats3-1-4);
  	  moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
	  if (detElemId % 2 == 1) {
	    if (detElemId == 517) // Round slat, new rotation due to mapping convention
	      GetEnvelopes(moduleId)->AddEnvelope(idSlatCh5, detElemId, true, TGeoTranslation(-xSlat3, -ySlat3, zSlat3 - dzCh3),
						  TGeoRotation("rot4",90,angle,90,270+angle,180,0) );
	    else
	      GetEnvelopes(moduleId)->AddEnvelope(idSlatCh5, detElemId, true, TGeoTranslation(-xSlat3, -ySlat3, zSlat3 - dzCh3),
						  TGeoRotation("rot4",90,180+angle,90,90+angle,180,0) );
	  }
	  else
	      GetEnvelopes(moduleId)->AddEnvelope(idSlatCh5, detElemId, true, TGeoTranslation(-xSlat3, -ySlat3, zSlat3 - dzCh3),
						  TGeoRotation("rot4",90,180+angle,90,270+angle,0,0) );   
	}

	sprintf(idSlatCh6,"SLB%d",kNslats3-1+i);  
	detElemId = 609 - (i  + kNslats3-1-4);
  	moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
	if (detElemId == 608 || detElemId == 609) // Round slat, new rotation due to mapping convention
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh6, detElemId, true, TGeoTranslation(xSlat3, ySlat3, -zSlat3 + dzCh3),
					      TGeoRotation("rot5",90,180+angle,90,90+angle,180,0));
	else {
	  if (detElemId % 2 == 0) 
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh6, detElemId, true, TGeoTranslation(xSlat3, ySlat3, -zSlat3 + dzCh3),
						TGeoRotation("rot5",90,angle,90,90+angle,0,0));
	  else
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh6, detElemId, true, TGeoTranslation(xSlat3, ySlat3, -zSlat3 + dzCh3),
						TGeoRotation("rot5",90,angle,90,270+angle,180,0));
	}

	sprintf(idSlatCh6,"SLB%d",3*kNslats3-2+i);
	detElemId = 600 + (i + kNslats3-1-4);
  	moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
	if (detElemId == 600 || detElemId == 601) // Round slat, new rotation due to mapping convention
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh6, detElemId, true, TGeoTranslation(-xSlat3, ySlat3, zSlat3 - dzCh3),
					      TGeoRotation("rot6",90,angle,90,90+angle,0,0) );
	else {
	  if (detElemId % 2 == 1) 
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh6, detElemId, true, TGeoTranslation(-xSlat3, ySlat3, zSlat3 - dzCh3),
						TGeoRotation("rot6",90,180+angle,90,90+angle,180,0) ); 
	  else
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh6, detElemId, true, TGeoTranslation(-xSlat3, ySlat3, zSlat3 - dzCh3),
						TGeoRotation("rot6",90,180+angle,90,270+angle,0,0) );
	}

	if (i > 0) { 
	  sprintf(idSlatCh6,"SLB%d",kNslats3-1-i);
	  detElemId = 609 + (i + kNslats3-1-4);
  	  moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
	  if (detElemId % 2 == 0 ) {
	    if (detElemId == 610) // Round slat, new rotation due to mapping convention
	      GetEnvelopes(moduleId)->AddEnvelope(idSlatCh6, detElemId, true, TGeoTranslation(xSlat3, -ySlat3, -zSlat3 + dzCh3),
						  TGeoRotation("rot7",90,180+angle,90,270+angle,0,0) );
	    else
	      GetEnvelopes(moduleId)->AddEnvelope(idSlatCh6, detElemId, true, TGeoTranslation(xSlat3, -ySlat3, -zSlat3 + dzCh3),
						  TGeoRotation("rot7",90,angle,90,90+angle,0,0) );
	  }
	  else
	      GetEnvelopes(moduleId)->AddEnvelope(idSlatCh6, detElemId, true, TGeoTranslation(xSlat3, -ySlat3, -zSlat3 + dzCh3),
						  TGeoRotation("rot7",90,angle,90,270+angle,180,0) );

	  sprintf(idSlatCh6,"SLB%d",3*kNslats3-2-i);
	  detElemId = 618 - (i + kNslats3-1-4);
  	  moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
	  if (detElemId % 2 == 1) {
	    if (detElemId == 617) // Round slat, new rotation due to mapping convention
	      GetEnvelopes(moduleId)->AddEnvelope(idSlatCh6, detElemId, true, TGeoTranslation(-xSlat3, -ySlat3, zSlat3 - dzCh3),
						  TGeoRotation("rot8",90,angle,90,270+angle,180,0) );
	    else
	      GetEnvelopes(moduleId)->AddEnvelope(idSlatCh6, detElemId, true, TGeoTranslation(-xSlat3, -ySlat3, zSlat3 - dzCh3),
						  TGeoRotation("rot8",90,180+angle,90,90+angle,180,0) );
	  }
	  else
	      GetEnvelopes(moduleId)->AddEnvelope(idSlatCh6, detElemId, true, TGeoTranslation(-xSlat3, -ySlat3, zSlat3 - dzCh3),
						  TGeoRotation("rot8",90,180+angle,90,270+angle,0,0) );   
	}
      }
     
      // create the panel volume 
 
      gMC->Gsvolu("S05C","BOX",kCarbonMaterial,panelpar,3);
      gMC->Gsvolu("SB5C","BOX",kCarbonMaterial,panelpar2,3);
      gMC->Gsvolu("SC5C","BOX",kCarbonMaterial,panelpar3,3);
      gMC->Gsvolu("SD5C","BOX",kCarbonMaterial,panelpar,3);
      gMC->Gsvolu("S06C","BOX",kCarbonMaterial,panelpar,3);
      gMC->Gsvolu("SC6C","BOX",kCarbonMaterial,panelpar3,3);
      gMC->Gsvolu("SD6C","BOX",kCarbonMaterial,panelpar,3);
 
      // create the nomex volume (honey comb)

      gMC->Gsvolu("S05N","BOX",kNomexMaterial,nomexpar,3);
      gMC->Gsvolu("SB5N","BOX",kNomexMaterial,nomexpar2,3);
      gMC->Gsvolu("SC5N","BOX",kNomexMaterial,nomexpar3,3);
      gMC->Gsvolu("SD5N","BOX",kNomexMaterial,nomexpar,3);
      gMC->Gsvolu("S06N","BOX",kNomexMaterial,nomexpar,3);
      gMC->Gsvolu("SC6N","BOX",kNomexMaterial,nomexpar3,3);
      gMC->Gsvolu("SD6N","BOX",kNomexMaterial,nomexpar,3);
 
      // create the nomex volume (bulk)

      gMC->Gsvolu("S05X","BOX",kNomexBMaterial,nomexbpar,3);
      gMC->Gsvolu("SB5X","BOX",kNomexBMaterial,nomexbpar2,3);
      gMC->Gsvolu("SC5X","BOX",kNomexBMaterial,nomexbpar3,3);
      gMC->Gsvolu("SD5X","BOX",kNomexBMaterial,nomexbpar,3);
      gMC->Gsvolu("S06X","BOX",kNomexBMaterial,nomexbpar,3);
      gMC->Gsvolu("SC6X","BOX",kNomexBMaterial,nomexbpar3,3);
      gMC->Gsvolu("SD6X","BOX",kNomexBMaterial,nomexbpar,3);

      // create the insulating material volume 

      gMC->Gsvolu("S05I","BOX",kInsuMaterial,insupar,3);
      gMC->Gsvolu("SB5I","BOX",kInsuMaterial,insupar2,3);
      gMC->Gsvolu("SC5I","BOX",kInsuMaterial,insupar3,3);
      gMC->Gsvolu("SD5I","BOX",kInsuMaterial,insupar,3);
      gMC->Gsvolu("S06I","BOX",kInsuMaterial,insupar,3);
      gMC->Gsvolu("SC6I","BOX",kInsuMaterial,insupar3,3);
      gMC->Gsvolu("SD6I","BOX",kInsuMaterial,insupar,3);
 
      // create the PCB volume 

      gMC->Gsvolu("S05P","BOX",kPcbMaterial,pcbpar,3);
      gMC->Gsvolu("SB5P","BOX",kPcbMaterial,pcbpar2,3);
      gMC->Gsvolu("SC5P","BOX",kPcbMaterial,pcbpar3,3);
      gMC->Gsvolu("SD5P","BOX",kPcbMaterial,pcbpar,3);
      gMC->Gsvolu("S06P","BOX",kPcbMaterial,pcbpar,3);
      gMC->Gsvolu("SC6P","BOX",kPcbMaterial,pcbpar3,3);
      gMC->Gsvolu("SD6P","BOX",kPcbMaterial,pcbpar,3);
 
      // create the sensitive volumes,

      gMC->Gsvolu("S05G","BOX",kSensMaterial,dum,0);
      gMC->Gsvolu("SC5G","BOX",kSensMaterial,senspar3,3);
      gMC->Gsvolu("SD5G","BOX",kSensMaterial,senspar,3);
      gMC->Gsvolu("S06G","BOX",kSensMaterial,dum,0);
      gMC->Gsvolu("SC6G","BOX",kSensMaterial,senspar3,3);
      gMC->Gsvolu("SD6G","BOX",kSensMaterial,senspar,3);

      // create the vertical frame volume 

      gMC->Gsvolu("S05V","BOX",kVframeMaterial,vFramepar,3);
      gMC->Gsvolu("S06V","BOX",kVframeMaterial,vFramepar,3);

      // create the rounded vertical frame volume 

      gMC->Gsvolu("SC5D","TUBS",kRframeMaterial,cFramepar3,5);
      gMC->Gsvolu("SD5D","TUBS",kRframeMaterial,rFramepar3,5);
      gMC->Gsvolu("SC6D","TUBS",kRframeMaterial,cFramepar3,5);
      gMC->Gsvolu("SD6D","TUBS",kRframeMaterial,rFramepar3,5);

      // create the horizontal frame volume 

      gMC->Gsvolu("S05H","BOX",kHframeMaterial,hFramepar,3);
      gMC->Gsvolu("SB5H","BOX",kHframeMaterial,hFramepar2,3);
      gMC->Gsvolu("SC5H","BOX",kHframeMaterial,hFramepar3,3);
      gMC->Gsvolu("SD5H","BOX",kHframeMaterial,hFramepar,3);
      gMC->Gsvolu("S06H","BOX",kHframeMaterial,hFramepar,3);
      gMC->Gsvolu("SC6H","BOX",kHframeMaterial,hFramepar3,3);
      gMC->Gsvolu("SD6H","BOX",kHframeMaterial,hFramepar,3);
 
      // create the horizontal border volume 

      gMC->Gsvolu("S05B","BOX",kBframeMaterial,bFramepar,3);
      gMC->Gsvolu("SB5B","BOX",kBframeMaterial,bFramepar2,3);
      gMC->Gsvolu("SC5B","BOX",kBframeMaterial,bFramepar3,3);
      gMC->Gsvolu("SD5B","BOX",kBframeMaterial,bFramepar,3);
      gMC->Gsvolu("S06B","BOX",kBframeMaterial,bFramepar,3);
      gMC->Gsvolu("SC6B","BOX",kBframeMaterial,bFramepar3,3);
      gMC->Gsvolu("SD6B","BOX",kBframeMaterial,bFramepar,3);

      // Replace the volume shape with a composite shape
      // with substracted overlap with beam shield     
      if ( gMC->IsRootGeometrySupported() ) { 
	
	// Get shape
	Int_t nSlatType = 2;
	Int_t nVol = 8;
	const char* slatType = "CD"; // C: central slat; D: rounded slat
	const char* volLetter = "CNXIPHBG";
	TString volName;
	TString compName;
	TString csName;
	TGeoVolume *mVol = 0x0;
	TObjArray centerSlat(nSlatType*((nVol+1)*2));	
	TObjArray composite(nSlatType*((nVol+1)*2));

	// Beam shield recess
	new TGeoTube("tubeCut", 0., AliMUONConstants::Rmin(2), kSlatWidth/2.+0.001);
	// Displacement
	TGeoTranslation* trCTube = new TGeoTranslation("trCTube", -(kPcbLength-csvPcbLength/2.+kVframeLength/2.), 0., 0.);
	trCTube->RegisterYourself();
	TGeoTranslation* trDTube = new TGeoTranslation("trDTube", -(kPcbLength+kVframeLength)/2., -kYpos3[1], 0.);
	trDTube->RegisterYourself();
	TGeoTranslation* trCBTube = new TGeoTranslation("trCBTube", 0., ( kPcbHeight - kBframeHeight ) / 2., 0.);
	trCBTube->Add(trCTube);
	trCBTube->RegisterYourself();
	TGeoTranslation* trDBTube = new TGeoTranslation("trDBTube", 0., ( kPcbHeight - kBframeHeight ) / 2., 0.);
	trDBTube->Add(trDTube);
	trDBTube->RegisterYourself();

	Float_t cPhi2 = (TMath::Pi()/2.-TMath::ACos((kSensHeight/2.)/(AliMUONConstants::Rmin(2)-kRframeLength)));
	TGeoBBox *boxCCut = new TGeoBBox("boxCCut",(cFramepar3[1]-cFramepar3[0]*TMath::Cos(cPhi2))/2., hFramepar3[1], cFramepar3[2]+0.001);
	// Displacement
	TGeoTranslation* trCBox = new TGeoTranslation("trCBox",cFramepar3[0]*TMath::Cos(cPhi2)+boxCCut->GetDX(), 0., 0.);
	trCBox->RegisterYourself();
	new TGeoBBox("boxDCut",(kPcbLength+kVframeLength)/2., hFramepar3[1], vFramepar[2]+0.001);
	// Displacement
	TGeoTranslation* trDBox = new TGeoTranslation("trDBox",kPcbLength/2., kYpos3[1], 0.);
	trDBox->RegisterYourself();

	TGeoBBox *boxVframe = new TGeoBBox("boxVframe",vFramepar[0],vrFrameHeight/2., vFramepar[2]);
	TGeoTranslation* trVBox = new TGeoTranslation("trVBox", 0., AliMUONConstants::Rmin(2)-kRframeLength + boxVframe->GetDY(), 0.);
	trVBox->RegisterYourself();

	for(int iCh=5; iCh<=6; iCh++){
	  for (int iSlatType = 0; iSlatType<nSlatType; iSlatType++) {
	    for (int iVol = 0; iVol<nVol; iVol++){
	      Int_t lIndex = (iCh-5)*(nSlatType*(nVol+1))+iSlatType*(nVol+1)+iVol;
	      volName=Form("S%c%d%c",slatType[iSlatType],iCh,volLetter[iVol]);
	      mVol = gGeoManager->FindVolumeFast(volName);
	      if ( !mVol ) {
		AliErrorStream() 
		  << "Slat volume " << volName << " not found" << endl;	 
	      }
	      else {
		centerSlat[lIndex] = mVol->GetShape();
		csName=Form("centerSlat%c%d%c",slatType[iSlatType],iCh,volLetter[iVol]);
		((TGeoShape*)centerSlat[lIndex])->SetName(csName); 

		// Composite shape
		TString compOperation(csName);
		compOperation+="-tubeCut:tr";
		compOperation+=slatType[iSlatType];
		if (strstr(volName,"B")){
		  compOperation+="B";
		}
		compOperation+="Tube";
		compName=Form("composite%c%d%c",slatType[iSlatType],iCh,volLetter[iVol]);
		composite[lIndex] = new TGeoCompositeShape(compName, compOperation.Data()); 
		// Reset shape to volume      
		mVol->SetShape((TGeoShape*)composite[lIndex]);
	      }
	    }

	    // For rounded spacer
	    Int_t lIndex = (iCh-5)*(nSlatType*(nVol+1))+iSlatType*(nVol+1)+nVol;
	    volName=Form("S%c%dD",slatType[iSlatType],iCh);
	    mVol = gGeoManager->FindVolumeFast(volName);
	    if ( !mVol ) {
	      AliErrorStream() 
		<< "Slat volume " << volName << " not found" << endl;	 
	    }
	    else {
	      centerSlat[lIndex] = mVol->GetShape();
	      csName=Form("centerSlat%c%dD",slatType[iSlatType],iCh);
	      ((TGeoShape*)centerSlat[lIndex])->SetName(csName);	 	  
	      
	      // Composite shape
	      TString compOperation(csName);
	      if (strstr(volName,"SC")){
		compOperation+="*boxCCut:trCBox";
	      }
	      if (strstr(volName,"SD")){
		compOperation.Prepend("(");
		compOperation+="+boxVframe:trVBox)*boxDCut:trDBox";
	      }
	      compName=Form("composite%c%dD",slatType[iSlatType],iCh);
	      composite[lIndex] = new TGeoCompositeShape(compName, compOperation.Data()); 	      
	      // Reset shape to volume      
	      mVol->SetShape((TGeoShape*)composite[lIndex]);
	    }
	  }
	}
      }
      	
      index = 0; 
      for (i = 0; i<kNslats3; i++){
	for (Int_t quadrant = 1; quadrant <= 4; quadrant++) {

	  if (i == 0 && quadrant == 2) continue;
	  if (i == 0 && quadrant == 4) continue;

	  sprintf(idSlatCh5,"SLA%d",ConvertSlatNum(i,quadrant,kNslats3-1));
	  sprintf(idSlatCh6,"SLB%d",ConvertSlatNum(i,quadrant,kNslats3-1));
	  Int_t moduleSlatCh5 = GetModuleId(idSlatCh5);
	  Int_t moduleSlatCh6 = GetModuleId(idSlatCh6);
	  Float_t xvFrame  = (slatLength3[i] - kVframeLength)/2.;
	  Float_t xvFrame2  = xvFrame;	  

	  if (i == 0 || i == 1 || i == 2) xvFrame2 -= pcbDLength3; // Correct position (J.C.)

	  // position the vertical frames 
	  if ( i > 2) { 
	    GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("S05V", idSlatCh5, 
						    (2*i+1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
	    GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("S05V", idSlatCh5, 
						    (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));
	    GetEnvelopes(moduleSlatCh6)->AddEnvelopeConstituent("S06V", idSlatCh6, 
						    (2*i+1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
	    GetEnvelopes(moduleSlatCh6)->AddEnvelopeConstituent("S06V", idSlatCh6, 
						    (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));	  
	  } 

	  if (i == 2) {
	    GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("S05V", idSlatCh5, 
						    (2*i+1)*10+quadrant,TGeoTranslation(xvFrame2,0.,0.)); 
	    GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("S05V", idSlatCh5, 
						    (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));
	    GetEnvelopes(moduleSlatCh6)->AddEnvelopeConstituent("S06V", idSlatCh6, 
						    (2*i+1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
	    GetEnvelopes(moduleSlatCh6)->AddEnvelopeConstituent("S06V", idSlatCh6, 
						    (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));
	  }

	  // Different rotation due to new mapping convention
	  if (i == 0 || i == 1) { // first vertical spacers
	    GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("S05V", idSlatCh5, 
						    (2*i+1)*10+quadrant,TGeoTranslation(-xvFrame2,0.,0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0)); 
	    GetEnvelopes(moduleSlatCh6)->AddEnvelopeConstituent("S06V", idSlatCh6, 
						    (2*i+1)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));	  
	    if (i == 0) { // rounded spacer for central slat (J.C.)
	      GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("SC5D", idSlatCh5, 
								  (2*i)*10+quadrant,TGeoTranslation(xvFrame,0.,0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	      GetEnvelopes(moduleSlatCh6)->AddEnvelopeConstituent("SC6D", idSlatCh6, 
								  (2*i)*10+quadrant,TGeoTranslation(xvFrame,0.,0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));

	    }
	    if (i == 1) { // rounded + vertical spacer for rounded slat (J.C.)
	      GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("SD5D", idSlatCh5, 
								  (2*i)*10+quadrant,TGeoTranslation(xvFrame,-kYpos3[1],0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	      GetEnvelopes(moduleSlatCh6)->AddEnvelopeConstituent("SD6D", idSlatCh6, 
								  (2*i)*10+quadrant,TGeoTranslation(xvFrame,-kYpos3[1],0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	    }
	  }

	  // position the panels and the insulating material 
	  for (j = 0; j < kNPCB3[i]; j++){
 	    index++;
	    xx = kSensLength * (-kNPCB3[i]/2. + j + 0.5); 
	    Float_t xx2 = xx - pcbDLength3/2.; 
	    Float_t xx3 = xx + (kSensLength-csvPcbLength)/2.;

	    Float_t zPanel = spar[2] - nomexbpar[2]; 

	    if (i==0 || i==1) { // Different rotation due to new mapping convention
	      if (i==0 && j==0) { // Rounded pcb of central slat (SR1, NR1)		
		GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("SC5X", idSlatCh5, 2*index-1,TGeoTranslation(-xx3,0.,zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
		GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("SC5X", idSlatCh5, 2*index,TGeoTranslation(-xx3,0.,-zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
		GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("SC5I", idSlatCh5, index,TGeoTranslation(-xx3,0.,0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
		GetEnvelopes(moduleSlatCh6)->AddEnvelopeConstituent("SC6X", idSlatCh6, 2*index-1,TGeoTranslation(-xx3,0.,zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
		GetEnvelopes(moduleSlatCh6)->AddEnvelopeConstituent("SC6X", idSlatCh6, 2*index,TGeoTranslation(-xx3,0.,-zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
		GetEnvelopes(moduleSlatCh6)->AddEnvelopeConstituent("SC6I", idSlatCh6, index,TGeoTranslation(-xx3,0.,0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	      } 
	      else {
		if (i==1 && j==0){ // Rounded pcb of rounded slats (SR2. NR2)
		  GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("SD5X", idSlatCh5, 2*index-1,TGeoTranslation(-xx,0.,zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
		  GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("SD5X", idSlatCh5, 2*index,TGeoTranslation(-xx,0.,-zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
		  GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("SD5I", idSlatCh5, index,TGeoTranslation(-xx,0.,0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
		  GetEnvelopes(moduleSlatCh6)->AddEnvelopeConstituent("SD6X", idSlatCh6, 2*index-1,TGeoTranslation(-xx,0.,zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
		  GetEnvelopes(moduleSlatCh6)->AddEnvelopeConstituent("SD6X", idSlatCh6, 2*index,TGeoTranslation(-xx,0.,-zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
		  GetEnvelopes(moduleSlatCh6)->AddEnvelopeConstituent("SD6I", idSlatCh6, index,TGeoTranslation(-xx,0.,0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
		}
		else {
		  if (j == kNPCB3[i]-1) { // 1 pcb is shortened by 5cm 
		    GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("SB5X", idSlatCh5, 2*index-1,TGeoTranslation(-xx2,0.,zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
		    GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("SB5X", idSlatCh5, 2*index,TGeoTranslation(-xx2,0.,-zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
		    GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("SB5I", idSlatCh5, index,TGeoTranslation(-xx2,0.,0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
		  } 
		  else {		
		    GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("S05X", idSlatCh5, 2*index-1,TGeoTranslation(-xx,0.,zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
		    GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("S05X", idSlatCh5, 2*index,TGeoTranslation(-xx,0.,-zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
		    GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("S05I", idSlatCh5, index,TGeoTranslation(-xx,0.,0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
		  }
		  GetEnvelopes(moduleSlatCh6)->AddEnvelopeConstituent("S06X", idSlatCh6, 2*index-1,TGeoTranslation(-xx,0.,zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
		  GetEnvelopes(moduleSlatCh6)->AddEnvelopeConstituent("S06X", idSlatCh6, 2*index,TGeoTranslation(-xx,0.,-zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
		  GetEnvelopes(moduleSlatCh6)->AddEnvelopeConstituent("S06I", idSlatCh6, index,TGeoTranslation(-xx,0.,0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
		}
	      }
	    }
	    else { 	      
	      if (i==2 && j == kNPCB3[i]-1) { // 1 pcb is shortened by 5cm 
		GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("SB5X", idSlatCh5, 2*index-1,TGeoTranslation(xx2,0.,zPanel));
		GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("SB5X", idSlatCh5, 2*index,TGeoTranslation(xx2,0.,-zPanel));
		GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("SB5I", idSlatCh5, index,TGeoTranslation(xx2,0.,0.));
	      } else {		
		GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("S05X", idSlatCh5, 2*index-1,TGeoTranslation(xx,0.,zPanel));
		GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("S05X", idSlatCh5, 2*index,TGeoTranslation(xx,0.,-zPanel));
		GetEnvelopes(moduleSlatCh5)->AddEnvelopeConstituent("S05I", idSlatCh5, index,TGeoTranslation(xx,0.,0.));
	      }
	      GetEnvelopes(moduleSlatCh6)->AddEnvelopeConstituent("S06X", idSlatCh6, 2*index-1,TGeoTranslation(xx,0.,zPanel));
	      GetEnvelopes(moduleSlatCh6)->AddEnvelopeConstituent("S06X", idSlatCh6, 2*index,TGeoTranslation(xx,0.,-zPanel));
	      GetEnvelopes(moduleSlatCh6)->AddEnvelopeConstituent("S06I", idSlatCh6, index,TGeoTranslation(xx,0.,0.));	    		
	    }
	  }
	}
      }

      
      // position the nomex volume inside the panel volume
      gMC->Gspos("S05N",1,"S05C",0.,0.,0.,0,"ONLY"); 
      gMC->Gspos("SB5N",1,"SB5C",0.,0.,0.,0,"ONLY"); 
      gMC->Gspos("SC5N",1,"SC5C",0.,0.,0.,0,"ONLY"); 
      gMC->Gspos("SD5N",1,"SD5C",0.,0.,0.,0,"ONLY"); 
      gMC->Gspos("S06N",1,"S06C",0.,0.,0.,0,"ONLY"); 
      gMC->Gspos("SC6N",1,"SC6C",0.,0.,0.,0,"ONLY"); 
      gMC->Gspos("SD6N",1,"SD6C",0.,0.,0.,0,"ONLY"); 
  
      // position panel volume inside the bulk nomex material volume
      gMC->Gspos("S05C",1,"S05X",0.,0.,kNomexBWidth/2.,0,"ONLY"); 
      gMC->Gspos("SB5C",1,"SB5X",0.,0.,kNomexBWidth/2.,0,"ONLY"); 
      gMC->Gspos("SC5C",1,"SC5X",0.,0.,kNomexBWidth/2.,0,"ONLY"); 
      gMC->Gspos("SD5C",1,"SD5X",0.,0.,kNomexBWidth/2.,0,"ONLY"); 
      gMC->Gspos("S06C",1,"S06X",0.,0.,kNomexBWidth/2.,0,"ONLY"); 
      gMC->Gspos("SC6C",1,"SC6X",0.,0.,kNomexBWidth/2.,0,"ONLY"); 
      gMC->Gspos("SD6C",1,"SD6X",0.,0.,kNomexBWidth/2.,0,"ONLY"); 

      // position the PCB volume inside the insulating material volume
      gMC->Gspos("S05P",1,"S05I",0.,0.,0.,0,"ONLY"); 
      gMC->Gspos("SB5P",1,"SB5I",0.,0.,0.,0,"ONLY"); 
      gMC->Gspos("SC5P",1,"SC5I",0.,0.,0.,0,"ONLY"); 
      gMC->Gspos("SD5P",1,"SD5I",0.,0.,0.,0,"ONLY"); 
      gMC->Gspos("S06P",1,"S06I",0.,0.,0.,0,"ONLY"); 
      gMC->Gspos("SC6P",1,"SC6I",0.,0.,0.,0,"ONLY"); 
      gMC->Gspos("SD6P",1,"SD6I",0.,0.,0.,0,"ONLY"); 
  
      // position the horizontal frame volume inside the PCB volume
      gMC->Gspos("S05H",1,"S05P",0.,0.,0.,0,"ONLY"); 
      gMC->Gspos("SB5H",1,"SB5P",0.,0.,0.,0,"ONLY"); 
      gMC->Gspos("SC5H",1,"SC5P",0.,0.,0.,0,"ONLY"); 
      gMC->Gspos("SD5H",1,"SD5P",0.,0.,0.,0,"ONLY"); 
      gMC->Gspos("S06H",1,"S06P",0.,0.,0.,0,"ONLY"); 
      gMC->Gspos("SC6H",1,"SC6P",0.,0.,0.,0,"ONLY"); 
      gMC->Gspos("SD6H",1,"SD6P",0.,0.,0.,0,"ONLY"); 
  
      // position the sensitive volume inside the horizontal frame volume
      gMC->Gsposp("S05G",1,"S05H",0.,0.,0.,0,"ONLY",senspar,3); 
      gMC->Gsposp("S05G",1,"SB5H",0.,0.,0.,0,"ONLY",senspar2,3); 
      gMC->Gspos("SC5G",1,"SC5H",0.,0.,0.,0,"ONLY"); 
      gMC->Gspos("SD5G",1,"SD5H",0.,0.,0.,0,"ONLY"); 
      gMC->Gsposp("S06G",1,"S06H",0.,0.,0.,0,"ONLY",senspar,3); 
      gMC->Gspos("SC6G",1,"SC6H",0.,0.,0.,0,"ONLY"); 
      gMC->Gspos("SD6G",1,"SD6H",0.,0.,0.,0,"ONLY"); 
  
 
      // position the border volumes inside the PCB volume
      Float_t yborder = ( kPcbHeight - kBframeHeight ) / 2.; 
      Int_t rotB = 0;
      gMC->Matrix(rotB,90,0,90,270,180,0); // rotation around x for second border

      gMC->Gspos("S05B",1,"S05P",0., yborder,0.,0,"ONLY"); 
      gMC->Gspos("S05B",2,"S05P",0.,-yborder,0.,0,"ONLY"); 
      gMC->Gspos("SB5B",1,"SB5P",0., yborder,0.,0,"ONLY"); 
      gMC->Gspos("SB5B",2,"SB5P",0.,-yborder,0.,0,"ONLY"); 
      gMC->Gspos("SC5B",1,"SC5P",0., yborder,0.,rotB,"ONLY"); 
      gMC->Gspos("SC5B",2,"SC5P",0.,-yborder,0.,0,"ONLY"); 
      gMC->Gspos("S05B",1,"SD5P",0., yborder,0.,0,"ONLY"); 
      gMC->Gspos("SD5B",1,"SD5P",0.,-yborder,0.,0,"ONLY"); 

      gMC->Gspos("S06B",1,"S06P",0., yborder,0.,0,"ONLY"); 
      gMC->Gspos("S06B",2,"S06P",0.,-yborder,0.,0,"ONLY"); 
      gMC->Gspos("SC6B",1,"SC6P",0., yborder,0.,rotB,"ONLY"); 
      gMC->Gspos("SC6B",2,"SC6P",0.,-yborder,0.,0,"ONLY"); 
      gMC->Gspos("S06B",1,"SD6P",0., yborder,0.,0,"ONLY"); 
      gMC->Gspos("SD6B",1,"SD6P",0.,-yborder,0.,0,"ONLY"); 
  
      // create the NULOC volume and position it in the horizontal frame
      gMC->Gsvolu("S05E","BOX",kNulocMaterial,nulocpar,3);
      gMC->Gsvolu("S06E","BOX",kNulocMaterial,nulocpar,3);
      index = 0;
      Float_t xxmax2 = xxmax - pcbDLength3/2.;
      Float_t xxmax3 = xxmax - (kPcbLength-csvPcbLength)/2.;
      Float_t rPhi3 = TMath::ASin((kYpos3[1]-kPcbHeight/2.)/AliMUONConstants::Rmin(2));
      Float_t xxmax4 = (AliMUONConstants::Rmin(2)*TMath::Cos(rPhi3)-kVframeLength/2.) - (kBframeLength - kNulocLength)/2.;
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
	if (xx > -xxmax3 && xx< xxmax3) {
	  gMC->Gspos("S05E",2*index-1,"SC5B", xx, 0.,-kBframeWidth/2.+ kNulocWidth/2., 0, "ONLY");
 	  gMC->Gspos("S05E",2*index  ,"SC5B", xx, 0., kBframeWidth/2.- kNulocWidth/2., 0, "ONLY");
 	  gMC->Gspos("S06E",2*index-1,"SC6B", xx, 0.,-kBframeWidth/2.+ kNulocWidth/2, 0, "ONLY");
 	  gMC->Gspos("S06E",2*index  ,"SC6B", xx, 0., kBframeWidth/2.- kNulocWidth/2, 0, "ONLY");
	}
	if (xx > xxmax4 && xx< xxmax) {
	  gMC->Gspos("S05E",2*index-1,"SD5B", xx, 0.,-kBframeWidth/2.+ kNulocWidth/2, 0, "ONLY");
	  gMC->Gspos("S05E",2*index  ,"SD5B", xx, 0., kBframeWidth/2.- kNulocWidth/2, 0, "ONLY");
	  gMC->Gspos("S06E",2*index-1,"SD6B", xx, 0.,-kBframeWidth/2.+ kNulocWidth/2, 0, "ONLY");
	  gMC->Gspos("S06E",2*index  ,"SD6B", xx, 0., kBframeWidth/2.- kNulocWidth/2, 0, "ONLY");
	}
      }      
      
      //
      //Geometry of the support pannel Verticla length 3.62m, horizontal length 1.62m, internal radius  dMotherInner of SC05 and SC06  (F. Orsini, Saclay)
      //Carbon fiber of 0.3 mm thick (2 layers) and a central layer of Nomex of 15mm thick. 
      // Outer excess and inner recess for mother volume radius
      // with respect to ROuter and RInner
      Float_t dMotherInner = AliMUONConstants::Rmin(2)-kRframeHeight; 
      Float_t nomexthickness = 1.5;
      Float_t carbonthickness = 0.03;
      Float_t supporthlength =  162.;  //    chamber 5 
      Float_t supporthlengthCh6 =  167.;  // chamber 6 
      Float_t supportvlength =  362.; 

      // Generating the composite shape of the carbon and nomex pannels
      new TGeoBBox("shNomexBoxSt3",supporthlength/2., supportvlength/2. ,nomexthickness/2.+carbonthickness+3*kCableWidth);
      new TGeoBBox("shCarbonBoxSt3",supporthlength/2., supportvlength/2. ,carbonthickness/2.); 
      new TGeoBBox("shNomexBoxSt3Ch6",(supporthlengthCh6)/2., supportvlength/2. ,nomexthickness/2.+carbonthickness+3*kCableWidth);
      new TGeoBBox("shCarbonBoxSt3Ch6",(supporthlengthCh6)/2., supportvlength/2. ,carbonthickness/2.); 
      new TGeoTubeSeg("shNomexHoleSt3",0., dMotherInner, nomexthickness/2.+carbonthickness+3*kCableWidth+0.001, -90. ,90.);
      new TGeoTubeSeg("shCarbonHoleSt3",0., dMotherInner, carbonthickness/2.+0.001, -90. ,90.);
      TGeoTranslation* trHoleSt3 = new TGeoTranslation("trHoleSt3",-supporthlength/2.,0.,0.); 
      trHoleSt3->RegisterYourself();
      TGeoTranslation* trHoleSt3Ch6 = new TGeoTranslation("trHoleSt3Ch6",-(supporthlengthCh6)/2.,0.,0.); 
      trHoleSt3Ch6->RegisterYourself();
      TGeoCompositeShape* shNomexSupportSt3  = new TGeoCompositeShape("shNomexSupportSt3","shNomexBoxSt3-shNomexHoleSt3:trHoleSt3");
      TGeoCompositeShape* shCarbonSupportSt3 = new TGeoCompositeShape("shCarbonSupportSt3","shCarbonBoxSt3-shCarbonHoleSt3:trHoleSt3");
      TGeoCompositeShape* shNomexSupportSt3Ch6  = new TGeoCompositeShape("shNomexSupportSt3Ch6","shNomexBoxSt3Ch6-shNomexHoleSt3:trHoleSt3Ch6");
      TGeoCompositeShape* shCarbonSupportSt3Ch6 = new TGeoCompositeShape("shCarbonSupportSt3Ch6","shCarbonBoxSt3Ch6-shCarbonHoleSt3:trHoleSt3Ch6");

      // Generating Nomex and Carbon pannel volumes
      TGeoVolume * voNomexSupportSt3  = new TGeoVolume("S05S", shNomexSupportSt3, kMedNomex);
      TGeoVolume * voCarbonSupportSt3 = new TGeoVolume("S05K", shCarbonSupportSt3, kMedCarbon);
      TGeoVolume * voNomexSupportSt3Ch6  = new TGeoVolume("S06S", shNomexSupportSt3Ch6, kMedNomex);
      TGeoVolume * voCarbonSupportSt3Ch6 = new TGeoVolume("S06K", shCarbonSupportSt3Ch6, kMedCarbon);

      TGeoTranslation *trCarbon1St3   = new TGeoTranslation("trCarbon1St3",0.,0., -(nomexthickness+carbonthickness)/2.);
      TGeoTranslation *trCarbon2St3   = new TGeoTranslation("trCarbon2St3",0.,0.,  (nomexthickness+carbonthickness)/2.);
      voNomexSupportSt3->AddNode(voCarbonSupportSt3,1,trCarbon1St3);
      voNomexSupportSt3->AddNode(voCarbonSupportSt3,2,trCarbon2St3);
      voNomexSupportSt3Ch6->AddNode(voCarbonSupportSt3Ch6,1,trCarbon1St3);
      voNomexSupportSt3Ch6->AddNode(voCarbonSupportSt3Ch6,2,trCarbon2St3);


      // Add readout cables
      gMC->Gsvolu("S05L","BOX",kCableMaterial,dum,0);
      gMC->Gsvolu("S06L","BOX",kCableMaterial,dum,0);

      ySlat3 = 0.;
      Float_t lCableX = 0.;
      Float_t lCableX6 = 0.;
      Float_t lCableY = 0.;
      Float_t lCableZ = 0.;
      Float_t cablepar[3] = {supporthlength/2., kCableHeight/2., kCableWidth/2.};
      Float_t cablepar6[3] = {supporthlengthCh6/2., kCableHeight/2., kCableWidth/2.};
      Float_t lCableDY = 0.;
      Int_t cIndex = 0;
      Int_t cIndex6 = 0;
      for (i = 0; i<kNslats3; i++){
	Int_t iCable = 1;
	cIndex = 0;
	cIndex6 = 0;
	ySlat3 += kYpos3[i];
	lCableY = ySlat3;

	// Cables going out from the start of slat
	if(kNPCB3[i]>=4 && i<kNslats3-2){ // Only if 4 or more pcb
	  // First top cables
	  cablepar[0] = supporthlength/2.;
	  lCableX = 0.;
	  cablepar6[0] = supporthlengthCh6/2.;
	  lCableX6 = 0.;
	  lCableDY = (kYpos3[i+1]+kYpos3[i+2])/2.-cablepar[1]; // half way between 2 slats on same side
	  lCableZ = TMath::Power(-1,i)*(nomexthickness/2.+carbonthickness+(-1+iCable++)*kCableWidth+kCableWidth/2.);
	  if(i==0){ // central slat is shorter (rounded)
	    cablepar[0] -= (kPcbLength-csvPcbLength)/2.;
	    lCableX = (kPcbLength-csvPcbLength)/2.;
	    cablepar6[0] -= (kPcbLength-csvPcbLength)/2.;
	    lCableX6 = (kPcbLength-csvPcbLength)/2.;
	  }
	  gMC->Gsposp("S05L",100*i+cIndex++,"S05S",lCableX,lCableY+lCableDY,lCableZ,0,"ONLY",cablepar,3); 	
	  gMC->Gsposp("S05L",100*i+cIndex++,"S05S",lCableX,-(lCableY+lCableDY),lCableZ,0,"ONLY",cablepar,3);
	  gMC->Gsposp("S06L",100*i+cIndex6++,"S06S",lCableX6,lCableY+lCableDY,lCableZ,0,"ONLY",cablepar6,3); 	
	  gMC->Gsposp("S06L",100*i+cIndex6++,"S06S",lCableX6,-(lCableY+lCableDY),lCableZ,0,"ONLY",cablepar6,3);

	  // Then bottom cables
	  if(i>0){
	    if(i==1){ // Rounded slat. Bottom cable starts at dMotherInner (beam pipe)
	      cablepar[0] -= dMotherInner/2.;
	      lCableX += dMotherInner/2.;
	      cablepar6[0] -= dMotherInner/2.;
	      lCableX6 += dMotherInner/2.;
	      lCableDY = (kYpos3[i]+kYpos3[i])/2.-cablepar[1];
	    }
	    if(i>=2){ 
	      lCableDY = (kYpos3[i]+kYpos3[i-1])/2.-cablepar[1]; // half way between 2 slats on same side
	      if ((lCableY-lCableDY)<(dMotherInner+cablepar[1])){
		lCableDY = lCableY - dMotherInner - cablepar[1];
	      }
	    }
	    gMC->Gsposp("S05L",100*i+cIndex++,"S05S",lCableX,lCableY-lCableDY,lCableZ,0,"ONLY",cablepar,3); 
	    gMC->Gsposp("S05L",100*i+cIndex++,"S05S",lCableX,-(lCableY-lCableDY),lCableZ,0,"ONLY",cablepar,3); 
	    gMC->Gsposp("S06L",100*i+cIndex6++,"S06S",lCableX6,lCableY-lCableDY,lCableZ,0,"ONLY",cablepar6,3); 
	    gMC->Gsposp("S06L",100*i+cIndex6++,"S06S",lCableX6,-(lCableY-lCableDY),lCableZ,0,"ONLY",cablepar6,3); 
	  }
	}
	
	// Rounded slats have an extra cable starting at second pcb
	if(i==1){ 
	  // First top cables
	  cablepar[0] = (supporthlength-kPcbLength-kVframeLength)/2.;
	  lCableX = (kPcbLength+kVframeLength)/2.;
	  cablepar6[0] = (supporthlengthCh6-kPcbLength-kVframeLength)/2.;
	  lCableX6 = (kPcbLength+kVframeLength)/2.;
	  lCableDY = (kYpos3[i+1]+kYpos3[i+2])/2.-cablepar[1]; // half way between 2 slats on same side
	  lCableZ = TMath::Power(-1,i)*(nomexthickness/2.+carbonthickness+(-1+iCable++)*kCableWidth+kCableWidth/2.);
	  gMC->Gsposp("S05L",100*i+cIndex++,"S05S",lCableX,lCableY+lCableDY,lCableZ,0,"ONLY",cablepar,3); 	
	  gMC->Gsposp("S05L",100*i+cIndex++,"S05S",lCableX,-(lCableY+lCableDY),lCableZ,0,"ONLY",cablepar,3);
	  gMC->Gsposp("S06L",100*i+cIndex6++,"S06S",lCableX6,lCableY+lCableDY,lCableZ,0,"ONLY",cablepar6,3); 	
	  gMC->Gsposp("S06L",100*i+cIndex6++,"S06S",lCableX6,-(lCableY+lCableDY),lCableZ,0,"ONLY",cablepar6,3);
	  // Then bottom cables
	  lCableDY = (kYpos3[i]+kYpos3[i])/2.-cablepar[1];
	  gMC->Gsposp("S05L",100*i+cIndex++,"S05S",lCableX,lCableY-lCableDY,lCableZ,0,"ONLY",cablepar,3); 
	  gMC->Gsposp("S05L",100*i+cIndex++,"S05S",lCableX,-(lCableY-lCableDY),lCableZ,0,"ONLY",cablepar,3); 
	  gMC->Gsposp("S06L",100*i+cIndex6++,"S06S",lCableX6,lCableY-lCableDY,lCableZ,0,"ONLY",cablepar6,3); 
	  gMC->Gsposp("S06L",100*i+cIndex6++,"S06S",lCableX6,-(lCableY-lCableDY),lCableZ,0,"ONLY",cablepar6,3); 
	}
	
	// Cables going out from the end of the slats
	// First top cables
	cablepar[0] = (supporthlength-(slatLength3[i]+kDslatLength)+kVframeLength)/2.;
	lCableX = slatLength3[i]-kVframeLength+kDslatLength+cablepar[0]-supporthlength/2.;
	cablepar6[0] = (supporthlengthCh6-(slatLength3[i]+kDslatLength)+kVframeLength)/2.;
	lCableX6 = slatLength3[i]-kVframeLength+kDslatLength+cablepar6[0]-supporthlengthCh6/2.;
	if(i+1>=kNslats3 || i+2>=kNslats3){ // If no more higher slats, then use distance to lower slat
	  lCableDY = kPcbHeight/2.+cablepar[1];
	}
	else {
	  lCableDY = (kYpos3[i+1]+kYpos3[i+2])/2.-cablepar[1];
	}
	lCableZ = TMath::Power(-1,i)*(nomexthickness/2.+carbonthickness+(-1+iCable++)*kCableWidth+kCableWidth/2.);

	if (i<=2){ // shortened pcb
	  cablepar[0] += pcbDLength3/2.;
	  lCableX -=  pcbDLength3/2.;
	} 
	gMC->Gsposp("S05L",100*i+cIndex++,"S05S",lCableX,lCableY+lCableDY,lCableZ,0,"ONLY",cablepar,3); 
	gMC->Gsposp("S05L",100*i+cIndex++,"S05S",lCableX,-(lCableY+lCableDY),lCableZ,0,"ONLY",cablepar,3); 
	gMC->Gsposp("S06L",100*i+cIndex6++,"S06S",lCableX6,lCableY+lCableDY,lCableZ,0,"ONLY",cablepar6,3); 
	gMC->Gsposp("S06L",100*i+cIndex6++,"S06S",lCableX6,-(lCableY+lCableDY),lCableZ,0,"ONLY",cablepar6,3); 
	// Then bottom cables
	if(i>0){ // Loop is over top half of slats, lower half are symmetric 
	  if (i==1) {
	    lCableDY = (kYpos3[i]+kYpos3[i])/2.-cablepar[1];
	  }
	  else{
	    lCableDY = (kYpos3[i]+kYpos3[i-1])/2.-cablepar[1]; // half way between 2 slats on same side
	  }
	  gMC->Gsposp("S05L",100*i+cIndex++,"S05S",lCableX,lCableY-lCableDY,lCableZ,0,"ONLY",cablepar,3); 
	  gMC->Gsposp("S05L",100*i+cIndex++,"S05S",lCableX,-(lCableY-lCableDY),lCableZ,0,"ONLY",cablepar,3); 
	  gMC->Gsposp("S06L",100*i+cIndex6++,"S06S",lCableX6,lCableY-lCableDY,lCableZ,0,"ONLY",cablepar6,3); 
	  gMC->Gsposp("S06L",100*i+cIndex6++,"S06S",lCableX6,-(lCableY-lCableDY),lCableZ,0,"ONLY",cablepar6,3); 
	}
      }

      Float_t dzCh5  = dzCh;
      TGeoTranslation* trSupport1St3   = new TGeoTranslation("trSupport1St3", supporthlength/2., 0. , dzCh5);
      TGeoRotation*    roSupportSt3    = new TGeoRotation("roSupportSt3",90.,180.,-90.);
      TGeoCombiTrans*  coSupport2St3   = new TGeoCombiTrans(-supporthlength/2., 0., -dzCh5, roSupportSt3);       
      TGeoTranslation* trSupport1St3Ch6   = new TGeoTranslation("trSupport1St3Ch6", supporthlengthCh6/2., 0. , dzCh5);
      TGeoCombiTrans*  coSupport2St3Ch6   = new TGeoCombiTrans(-supporthlengthCh6/2., 0., -dzCh5, roSupportSt3);       
      GetEnvelopes(5)->AddEnvelope("S05S", 0, 1, *trSupport1St3);  
      GetEnvelopes(4)->AddEnvelope("S05S", 0, 2, *coSupport2St3);  
      GetEnvelopes(7)->AddEnvelope("S06S", 0, 1, *trSupport1St3Ch6);   
      GetEnvelopes(6)->AddEnvelope("S06S", 0, 2, *coSupport2St3Ch6);  
      // End of pannel support geometry          

      // cout << "Geometry for Station 3...... done" << endl;	
    }
  if (fStations[3]) {


    // //********************************************************************
    // //                            Station 4                             **
    // //********************************************************************
    // Mother volume for each chamber in St4 is an envelop (or assembly)
    // There is one assembly mother per half a chamber  called SC07I, SC07O, SC08I and SC08O
    // Same volume name definitions as in St3
    const Int_t   kNslats4          = 7;  // number of slats per quadrant
    const Int_t   kNPCB4[kNslats4]  = {5, 6, 5, 5, 4, 3, 2}; // n PCB per slat
    const Float_t kXpos4[kNslats4]  = {38.75, 0., 0., 0., 0., 0., 0.}; // J.C. Correct value
    const Float_t kYpos41[kNslats4] = {0., 38.2, 34.40, 36.60, 29.3, 37.0, 28.6};
    const Float_t kYpos42[kNslats4] = {0., 38.2, 37.85, 37.55, 29.4, 37.0, 28.6};
    Float_t slatLength4[kNslats4];     

    Float_t rPhi1 = TMath::RadToDeg()*(TMath::ASin((kYpos41[1]-hFramepar[1])/(AliMUONConstants::Rmin(3))));
    Float_t rPhi2 = TMath::RadToDeg()*(TMath::ACos(-vFramepar[0]/(AliMUONConstants::Rmin(3)-kRframeLength)));
    Float_t rFramepar4[5] = { AliMUONConstants::Rmin(3)-kRframeLength, AliMUONConstants::Rmin(3), kRframeWidth, rPhi1, rPhi2}; 
    Float_t vrFrameHeight = hFramepar[1]+kYpos41[1]-AliMUONConstants::Rmin(3)+kRframeLength;
    
    char idSlatCh7[6];
    char idSlatCh8[6];
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

      sprintf(idSlatCh7,"SLC%d",kNslats4-1+i);
      detElemId = 713 - (i + kNslats4-1-6);
      moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
      if (detElemId % 2 == 0) {
	if (detElemId == 712) // Round slat, new rotation due to mapping convention
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh7, detElemId, true, TGeoTranslation(xSlat4, ySlat41, -zSlat4 + dzCh4),
					      TGeoRotation("rot1",90,180+angle,90,90+angle,180,0) );
	else
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh7, detElemId, true, TGeoTranslation(xSlat4, ySlat41, -zSlat4 + dzCh4),
					      TGeoRotation("rot1",90,angle,90,90+angle,0,0) );
      }
      else
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh7, detElemId, true, TGeoTranslation(xSlat4, ySlat41, -zSlat4 + dzCh4),
				   TGeoRotation("rot1",90,angle,90,270+angle,180,0) );

      sprintf(idSlatCh7,"SLC%d",3*kNslats4-2+i);
      detElemId = 700 + (i + kNslats4-1-6);
      moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
      if (detElemId % 2 == 1) {
	if (detElemId == 701) // Round slat, new rotation due to mapping convention
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh7, detElemId, true, TGeoTranslation(-xSlat4, ySlat41, zSlat4 - dzCh4),
					      TGeoRotation("rot2",90,angle,90,90+angle,0,0) );
	else
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh7, detElemId, true, TGeoTranslation(-xSlat4, ySlat41, zSlat4 - dzCh4),
					      TGeoRotation("rot2",90,180+angle,90,90+angle,180,0) );
      }
      else
	GetEnvelopes(moduleId)->AddEnvelope(idSlatCh7, detElemId, true, TGeoTranslation(-xSlat4, ySlat41, zSlat4 - dzCh4),
					    TGeoRotation("rot2",90,180+angle,90,270+angle,0,0) );

      if (i > 0) { 
	sprintf(idSlatCh7,"SLC%d",kNslats4-1-i);
	detElemId = 713 + (i + kNslats4-1-6);
        moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
	if (detElemId % 2 == 0) {
	  if (detElemId == 714) // Round slat, new rotation due to mapping convention
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh7, detElemId, true, TGeoTranslation(xSlat4, -ySlat41, -zSlat4 + dzCh4),
						TGeoRotation("rot3",90,180+angle,90,270+angle,0,0) );
	  else
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh7, detElemId, true, TGeoTranslation(xSlat4, -ySlat41, -zSlat4 + dzCh4),
						TGeoRotation("rot3",90,angle,90,90+angle,0,0) );
	}
	else
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh7, detElemId, true, TGeoTranslation(xSlat4, -ySlat41, -zSlat4 + dzCh4),
						TGeoRotation("rot3",90,angle,90,270+angle,180,0) );

	sprintf(idSlatCh7,"SLC%d",3*kNslats4-2-i);
	detElemId = 726 - (i + kNslats4-1-6);
        moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
	if (detElemId % 2 == 1) {
	  if (detElemId == 725) // Round slat, new rotation due to mapping convention
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh7, detElemId, true, TGeoTranslation(-xSlat4, -ySlat41, zSlat4 - dzCh4),
						TGeoRotation("rot4",90,angle,90,270+angle,180,0) );
	  else
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh7, detElemId, true, TGeoTranslation(-xSlat4, -ySlat41, zSlat4 - dzCh4),
						TGeoRotation("rot4",90,180+angle,90,90+angle,180,0) );
	}
	else
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh7, detElemId, true, TGeoTranslation(-xSlat4, -ySlat41, zSlat4 - dzCh4),
					      TGeoRotation("rot4",90,180+angle,90,270+angle,0,0) ); 
      }

      sprintf(idSlatCh8,"SLD%d",kNslats4-1+i);
      detElemId = 813 - (i + kNslats4-1-6);
      moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
      if (detElemId % 2 == 0) {
	if (detElemId == 812) // Round slat, new rotation due to mapping convention
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh8, detElemId, true, TGeoTranslation(xSlat4, ySlat42, -zSlat4 + dzCh4),
					      TGeoRotation("rot5",90,180+angle,90,90+angle,180,0) );
	else
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh8, detElemId, true, TGeoTranslation(xSlat4, ySlat42, -zSlat4 + dzCh4),
					      TGeoRotation("rot5",90,angle,90,90+angle,0,0) );
      }
      else
	GetEnvelopes(moduleId)->AddEnvelope(idSlatCh8, detElemId, true, TGeoTranslation(xSlat4, ySlat42, -zSlat4 + dzCh4),
					    TGeoRotation("rot5",90,angle,90,270+angle,180,0) ); 

      sprintf(idSlatCh8,"SLD%d",3*kNslats4-2+i);
      detElemId = 800 + (i + kNslats4-1-6);
      moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
      if (detElemId % 2 == 1) {
	if (detElemId == 801) // Round slat, new rotation due to mapping convention
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh8, detElemId, true, TGeoTranslation(-xSlat4, ySlat42, zSlat4 - dzCh4),
					      TGeoRotation("rot6",90,angle,90,90+angle,0,0) );
	else
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh8, detElemId, true, TGeoTranslation(-xSlat4, ySlat42, zSlat4 - dzCh4),
					      TGeoRotation("rot6",90,180+angle,90,90+angle,180,0) );
      }
      else
	GetEnvelopes(moduleId)->AddEnvelope(idSlatCh8, detElemId, true, TGeoTranslation(-xSlat4, ySlat42, zSlat4 - dzCh4),
					    TGeoRotation("rot6",90,180+angle,90,270+angle,0,0) );
      if (i > 0) { 
	sprintf(idSlatCh8,"SLD%d",kNslats4-1-i);
	detElemId = 813 + (i + kNslats4-1-6);
        moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
	if (detElemId % 2 == 0) {
	  if (detElemId == 814) // Round slat, new rotation due to mapping convention
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh8, detElemId, true, TGeoTranslation(xSlat4, -ySlat42, -zSlat4 + dzCh4),
						TGeoRotation("rot7",90,180+angle,90,270+angle,0,0) );
	  else
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh8, detElemId, true, TGeoTranslation(xSlat4, -ySlat42, -zSlat4 + dzCh4),
						TGeoRotation("rot7",90,angle,90,90+angle,0,0) );
	}
	else
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh8, detElemId, true, TGeoTranslation(xSlat4, -ySlat42, -zSlat4 + dzCh4),
					      TGeoRotation("rot7",90,angle,90,270+angle,180,0) );

	sprintf(idSlatCh8,"SLD%d",3*kNslats4-2-i);
	detElemId = 826 - (i + kNslats4-1-6);
        moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
	if (detElemId % 2 == 1) {
	  if (detElemId == 825 ) // Round slat, new rotation due to mapping convention
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh8, detElemId, true, TGeoTranslation(-xSlat4, -ySlat42, zSlat4 - dzCh4),
						TGeoRotation("rot8",90,angle,90,270+angle,180,0) ); 
	  else
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh8, detElemId, true, TGeoTranslation(-xSlat4, -ySlat42, zSlat4 - dzCh4),
						TGeoRotation("rot8",90,180+angle,90,90+angle,180,0) ); 
	}
	else
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh8, detElemId, true, TGeoTranslation(-xSlat4, -ySlat42, zSlat4 - dzCh4),
					      TGeoRotation("rot8",90,180+angle,90,270+angle,0,0) );
       
      }
    }
     
    // create the panel volume 
 
    gMC->Gsvolu("S07C","BOX",kCarbonMaterial,panelpar,3);
    gMC->Gsvolu("SD7C","BOX",kCarbonMaterial,panelpar,3);
    gMC->Gsvolu("S08C","BOX",kCarbonMaterial,panelpar,3);
    gMC->Gsvolu("SD8C","BOX",kCarbonMaterial,panelpar,3);

    // create the nomex volume 

    gMC->Gsvolu("S07N","BOX",kNomexMaterial,nomexpar,3);
    gMC->Gsvolu("SD7N","BOX",kNomexMaterial,nomexpar,3);
    gMC->Gsvolu("S08N","BOX",kNomexMaterial,nomexpar,3);
    gMC->Gsvolu("SD8N","BOX",kNomexMaterial,nomexpar,3);


    // create the nomex volume (bulk)

    gMC->Gsvolu("S07X","BOX",kNomexBMaterial,nomexbpar,3);
    gMC->Gsvolu("SD7X","BOX",kNomexBMaterial,nomexbpar,3);
    gMC->Gsvolu("S08X","BOX",kNomexBMaterial,nomexbpar,3);
    gMC->Gsvolu("SD8X","BOX",kNomexBMaterial,nomexbpar,3);

    // create the insulating material volume 

    gMC->Gsvolu("S07I","BOX",kInsuMaterial,insupar,3);
    gMC->Gsvolu("SD7I","BOX",kInsuMaterial,insupar,3);
    gMC->Gsvolu("S08I","BOX",kInsuMaterial,insupar,3);
    gMC->Gsvolu("SD8I","BOX",kInsuMaterial,insupar,3);

    // create the PCB volume 

    gMC->Gsvolu("S07P","BOX",kPcbMaterial,pcbpar,3);
    gMC->Gsvolu("SD7P","BOX",kPcbMaterial,pcbpar,3);
    gMC->Gsvolu("S08P","BOX",kPcbMaterial,pcbpar,3);
    gMC->Gsvolu("SD8P","BOX",kPcbMaterial,pcbpar,3);
 
    // create the sensitive volumes,

    gMC->Gsvolu("S07G","BOX",kSensMaterial,dum,0);
    gMC->Gsvolu("SD7G","BOX",kSensMaterial,senspar,3);
    gMC->Gsvolu("S08G","BOX",kSensMaterial,dum,0);
    gMC->Gsvolu("SD8G","BOX",kSensMaterial,senspar,3);

    // create the vertical frame volume 

    gMC->Gsvolu("S07V","BOX",kVframeMaterial,vFramepar,3);
    gMC->Gsvolu("S08V","BOX",kVframeMaterial,vFramepar,3);

    // create the rounded vertical frame volume 

    gMC->Gsvolu("SD7D","TUBS",kRframeMaterial,rFramepar4,5);
    gMC->Gsvolu("SD8D","TUBS",kRframeMaterial,rFramepar4,5);
    
    // create the horizontal frame volume 

    gMC->Gsvolu("S07H","BOX",kHframeMaterial,hFramepar,3);
    gMC->Gsvolu("SD7H","BOX",kHframeMaterial,hFramepar,3);
    gMC->Gsvolu("S08H","BOX",kHframeMaterial,hFramepar,3);
    gMC->Gsvolu("SD8H","BOX",kHframeMaterial,hFramepar,3);

    // create the horizontal border volume 

    gMC->Gsvolu("S07B","BOX",kBframeMaterial,bFramepar,3);
    gMC->Gsvolu("SD7B","BOX",kBframeMaterial,bFramepar,3);
    gMC->Gsvolu("S08B","BOX",kBframeMaterial,bFramepar,3);
    gMC->Gsvolu("SD8B","BOX",kBframeMaterial,bFramepar,3);

    // Replace the volume shape with a composite shape
    // with substracted overlap with beam shield     
    if ( gMC->IsRootGeometrySupported() ) { 
	
      // Get shape
      Int_t nSlatType = 1;
      Int_t nVol = 8;
      const char* slatType = "D"; // D: Rounded slat
      const char* volLetter = "CNXIPHBG";
      TString volName;
      TString compName;
      TString csName;
      TGeoVolume *mVol = 0x0;
      // Beam shield recess
      new TGeoTube("tube4Cut", 0., AliMUONConstants::Rmin(3), kSlatWidth/2.+0.001);
      TObjArray rounded4Slat(nSlatType*((nVol+1)*2));	
      // Displacement
      TGeoTranslation* trDTube4 = new TGeoTranslation("trDTube4", -(kPcbLength+kVframeLength)/2., -kYpos41[1], 0.);
      trDTube4->RegisterYourself();
      TGeoTranslation* trDBTube4 = new TGeoTranslation("trDBTube4", 0., ( kPcbHeight - kBframeHeight ) / 2., 0.);
      trDBTube4->Add(trDTube4);
      trDBTube4->RegisterYourself();

      TObjArray composite4(nSlatType*((nVol+1)*2));
      new TGeoBBox("box4DCut",(kPcbLength+kVframeLength)/2., hFramepar[1], vFramepar[2]+0.001);
      // Displacement
      TGeoTranslation* trDBox4 = new TGeoTranslation("trDBox4",kPcbLength/2., kYpos41[1], 0.);
      trDBox4->RegisterYourself();      

      TGeoBBox *box4Vframe = new TGeoBBox("box4Vframe",vFramepar[0],vrFrameHeight/2., vFramepar[2]);
      TGeoTranslation* trVBox4 = new TGeoTranslation("trVBox4", 0., AliMUONConstants::Rmin(3)-kRframeLength + box4Vframe->GetDY(), 0.);
      trVBox4->RegisterYourself();
      
      for(int iCh=7; iCh<=8; iCh++){
	for (int iSlatType = 0; iSlatType<nSlatType; iSlatType++) {
	  for (int iVol = 0; iVol<nVol; iVol++){
	    Int_t lIndex = (iCh-7)*(nSlatType*(nVol+1))+iSlatType*(nVol+1)+iVol;
	    volName=Form("S%c%d%c",slatType[iSlatType],iCh,volLetter[iVol]);
	    mVol = gGeoManager->FindVolumeFast(volName);
	    if ( !mVol ) {
	      AliErrorStream() 
		<< "Slat volume " << volName << " not found" << endl;	 
	    }
	    else {
	      rounded4Slat[lIndex] = mVol->GetShape();
	      csName=Form("rounded4Slat%c%d%c",slatType[iSlatType],iCh,volLetter[iVol]);
	      ((TGeoShape*)rounded4Slat[lIndex])->SetName(csName);
	      
	      // Composite shape
	      TString compOperation(csName);
	      compOperation+="-tube4Cut:tr";
	      compOperation+=slatType[iSlatType];
	      if (strstr(volName,"B")){
		compOperation+="B";
	      }	      
	      compOperation+="Tube4";
	      compName=Form("composite4%c%d%c",slatType[iSlatType],iCh,volLetter[iVol]);
	      composite4[lIndex] = new TGeoCompositeShape(compName, compOperation.Data()); 
	      
	      // Reset shape to volume      
	      mVol->SetShape((TGeoShape*)composite4[lIndex]);
	    }
	  }

	  // For rounded spacer
	  Int_t lIndex = (iCh-7)*(nSlatType*(nVol+1))+iSlatType*(nVol+1)+nVol;
	  volName=Form("S%c%dD",slatType[iSlatType],iCh);
	  mVol = gGeoManager->FindVolumeFast(volName);
	  if ( !mVol ) {
	    AliErrorStream() 
	      << "Slat volume " << volName << " not found" << endl;	 
	  }
	  else {
	    rounded4Slat[lIndex] = mVol->GetShape();
	    csName=Form("rounded4Slat%c%dD",slatType[iSlatType],iCh);
	    ((TGeoShape*)rounded4Slat[lIndex])->SetName(csName);
	    
	    // Composite shape
	    TString compOperation(csName);
	    if (strstr(volName,"SD")){
	      compOperation.Prepend("(");
	      compOperation+="+box4Vframe:trVBox4)*box4DCut:trDBox4";
	    }
	    compName=Form("composite4%c%dD",slatType[iSlatType],iCh);
	    composite4[lIndex] = new TGeoCompositeShape(compName, compOperation.Data()); 	      
	    // Reset shape to volume      
	    mVol->SetShape((TGeoShape*)composite4[lIndex]);
	  }
	}
      }
    }


    index = 0; 
    for (i = 0; i < kNslats4; i++){
      for (Int_t quadrant = 1; quadrant <= 4; quadrant++) {

	if (i == 0 && quadrant == 2) continue;
	if (i == 0 && quadrant == 4) continue;

	sprintf(idSlatCh7,"SLC%d",ConvertSlatNum(i,quadrant,kNslats4-1));
	sprintf(idSlatCh8,"SLD%d",ConvertSlatNum(i,quadrant,kNslats4-1));
	Int_t moduleSlatCh7 = GetModuleId(idSlatCh7);
	Int_t moduleSlatCh8 = GetModuleId(idSlatCh8);

	Float_t xvFrame  = (slatLength4[i] - kVframeLength)/2.;

	// position the vertical frames 
	if (i != 1) { 
	  GetEnvelopes(moduleSlatCh7)->AddEnvelopeConstituent("S07V", idSlatCh7, (2*i+1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
	  GetEnvelopes(moduleSlatCh7)->AddEnvelopeConstituent("S07V", idSlatCh7, (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));
	  GetEnvelopes(moduleSlatCh8)->AddEnvelopeConstituent("S08V", idSlatCh8, (2*i+1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
	  GetEnvelopes(moduleSlatCh8)->AddEnvelopeConstituent("S08V", idSlatCh8, (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));
	} else { // Vertical and Rounded+Vertical spacer - Different rotation due to new mapping convention
	  GetEnvelopes(moduleSlatCh7)->AddEnvelopeConstituent("S07V", idSlatCh7, (2*i+1)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	  GetEnvelopes(moduleSlatCh7)->AddEnvelopeConstituent("SD7D", idSlatCh7, (2*i)*10+quadrant,TGeoTranslation(xvFrame,-kYpos41[1],0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	  GetEnvelopes(moduleSlatCh8)->AddEnvelopeConstituent("S08V", idSlatCh8, (2*i+1)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	  GetEnvelopes(moduleSlatCh8)->AddEnvelopeConstituent("SD8D", idSlatCh8, (2*i)*10+quadrant,TGeoTranslation(+xvFrame,-kYpos42[1],0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	}
	// position the panels and the insulating material 
	for (j = 0; j < kNPCB4[i]; j++){
	  index++;
	  xx = kSensLength * (-kNPCB4[i]/2.+j+.5); 
	  Float_t zPanel = spar[2] - nomexbpar[2]; 
	  if (i==1) { // Different rotation due to new mapping convention
	    if (j==0){ // Rounded pcb of rounded slat
	      GetEnvelopes(moduleSlatCh7)->AddEnvelopeConstituent("SD7X", idSlatCh7, 2*index-1,TGeoTranslation(-xx,0.,zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	      GetEnvelopes(moduleSlatCh7)->AddEnvelopeConstituent("SD7X", idSlatCh7, 2*index,TGeoTranslation(-xx,0.,-zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	      GetEnvelopes(moduleSlatCh7)->AddEnvelopeConstituent("SD7I", idSlatCh7, index,TGeoTranslation(-xx,0.,0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	      GetEnvelopes(moduleSlatCh8)->AddEnvelopeConstituent("SD8X", idSlatCh8, 2*index-1,TGeoTranslation(-xx,0.,zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	      GetEnvelopes(moduleSlatCh8)->AddEnvelopeConstituent("SD8X", idSlatCh8, 2*index,TGeoTranslation(-xx,0.,-zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	      GetEnvelopes(moduleSlatCh8)->AddEnvelopeConstituent("SD8I", idSlatCh8, index,TGeoTranslation(-xx,0.,0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	    } else { 	      
	      GetEnvelopes(moduleSlatCh7)->AddEnvelopeConstituent("S07X", idSlatCh7, 2*index-1,TGeoTranslation(-xx,0.,zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	      GetEnvelopes(moduleSlatCh7)->AddEnvelopeConstituent("S07X", idSlatCh7, 2*index,TGeoTranslation(-xx,0.,-zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	      GetEnvelopes(moduleSlatCh7)->AddEnvelopeConstituent("S07I", idSlatCh7, index,TGeoTranslation(-xx,0.,0.));
	      GetEnvelopes(moduleSlatCh8)->AddEnvelopeConstituent("S08X", idSlatCh8, 2*index-1,TGeoTranslation(-xx,0.,zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	      GetEnvelopes(moduleSlatCh8)->AddEnvelopeConstituent("S08X", idSlatCh8, 2*index,TGeoTranslation(-xx,0.,-zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	      GetEnvelopes(moduleSlatCh8)->AddEnvelopeConstituent("S08I", idSlatCh8, index,TGeoTranslation(-xx,0.,0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	    }
	  } else { 	      
	    GetEnvelopes(moduleSlatCh7)->AddEnvelopeConstituent("S07X", idSlatCh7, 2*index-1,TGeoTranslation(xx,0.,zPanel));
	    GetEnvelopes(moduleSlatCh7)->AddEnvelopeConstituent("S07X", idSlatCh7, 2*index,TGeoTranslation(xx,0.,-zPanel));
	    GetEnvelopes(moduleSlatCh7)->AddEnvelopeConstituent("S07I", idSlatCh7, index,TGeoTranslation(xx,0.,0.));
	    GetEnvelopes(moduleSlatCh8)->AddEnvelopeConstituent("S08X", idSlatCh8, 2*index-1,TGeoTranslation(xx,0.,zPanel));
	    GetEnvelopes(moduleSlatCh8)->AddEnvelopeConstituent("S08X", idSlatCh8, 2*index,TGeoTranslation(xx,0.,-zPanel));
	    GetEnvelopes(moduleSlatCh8)->AddEnvelopeConstituent("S08I", idSlatCh8, index,TGeoTranslation(xx,0.,0.));
	  }
	}
      } 
    }

    // position the nomex volume inside the panel volume
    gMC->Gspos("S07N",1,"S07C",0.,0.,0.,0,"ONLY"); 
    gMC->Gspos("SD7N",1,"SD7C",0.,0.,0.,0,"ONLY"); 
    gMC->Gspos("S08N",1,"S08C",0.,0.,0.,0,"ONLY"); 
    gMC->Gspos("SD8N",1,"SD8C",0.,0.,0.,0,"ONLY"); 

    // position panel volume inside the bulk nomex material volume
    gMC->Gspos("S07C",1,"S07X",0.,0.,kNomexBWidth/2.,0,"ONLY"); 
    gMC->Gspos("SD7C",1,"SD7X",0.,0.,kNomexBWidth/2.,0,"ONLY"); 
    gMC->Gspos("S08C",1,"S08X",0.,0.,kNomexBWidth/2.,0,"ONLY"); 
    gMC->Gspos("SD8C",1,"SD8X",0.,0.,kNomexBWidth/2.,0,"ONLY"); 

    // position the PCB volume inside the insulating material volume
    gMC->Gspos("S07P",1,"S07I",0.,0.,0.,0,"ONLY"); 
    gMC->Gspos("SD7P",1,"SD7I",0.,0.,0.,0,"ONLY"); 
    gMC->Gspos("S08P",1,"S08I",0.,0.,0.,0,"ONLY"); 
    gMC->Gspos("SD8P",1,"SD8I",0.,0.,0.,0,"ONLY"); 

    // position the horizontal frame volume inside the PCB volume
    gMC->Gspos("S07H",1,"S07P",0.,0.,0.,0,"ONLY"); 
    gMC->Gspos("SD7H",1,"SD7P",0.,0.,0.,0,"ONLY"); 
    gMC->Gspos("S08H",1,"S08P",0.,0.,0.,0,"ONLY"); 
    gMC->Gspos("SD8H",1,"SD8P",0.,0.,0.,0,"ONLY"); 

    // position the sensitive volume inside the horizontal frame volume
    gMC->Gsposp("S07G",1,"S07H",0.,0.,0.,0,"ONLY",senspar,3); 
    gMC->Gspos("SD7G",1,"SD7H",0.,0.,0.,0,"ONLY"); 
    gMC->Gsposp("S08G",1,"S08H",0.,0.,0.,0,"ONLY",senspar,3); 
    gMC->Gspos("SD8G",1,"SD8H",0.,0.,0.,0,"ONLY"); 

    // position the border volumes inside the PCB volume
    Float_t yborder = ( kPcbHeight - kBframeHeight ) / 2.; 
    gMC->Gspos("S07B",1,"S07P",0., yborder,0.,0,"ONLY"); 
    gMC->Gspos("S07B",2,"S07P",0.,-yborder,0.,0,"ONLY");
    gMC->Gspos("S07B",1,"SD7P",0., yborder,0.,0,"ONLY");
    gMC->Gspos("SD7B",1,"SD7P",0.,-yborder,0.,0,"ONLY");  
    gMC->Gspos("S08B",1,"S08P",0., yborder,0.,0,"ONLY"); 
    gMC->Gspos("S08B",2,"S08P",0.,-yborder,0.,0,"ONLY"); 
    gMC->Gspos("S08B",1,"SD8P",0., yborder,0.,0,"ONLY"); 
    gMC->Gspos("SD8B",1,"SD8P",0.,-yborder,0.,0,"ONLY"); 

    // create the NULOC volume and position it in the horizontal frame

    gMC->Gsvolu("S07E","BOX",kNulocMaterial,nulocpar,3);
    gMC->Gsvolu("S08E","BOX",kNulocMaterial,nulocpar,3);
    index = 0;
    Float_t rPhi3 = TMath::ASin((kYpos41[1]-kPcbHeight/2.)/AliMUONConstants::Rmin(3));
    Float_t xxmax4 = (AliMUONConstants::Rmin(3)*TMath::Cos(rPhi3)-kVframeLength/2.) - (kBframeLength - kNulocLength)/2.;
    for (xx = -xxmax; xx <= xxmax; xx += 2*kNulocLength) { 
      index++; 
      gMC->Gspos("S07E",2*index-1,"S07B", xx, 0.,-kBframeWidth/2. + kNulocWidth/2, 0, "ONLY");
      gMC->Gspos("S07E",2*index  ,"S07B", xx, 0., kBframeWidth/2. - kNulocWidth/2, 0, "ONLY");
      gMC->Gspos("S08E",2*index-1,"S08B", xx, 0.,-kBframeWidth/2. + kNulocWidth/2, 0, "ONLY");
      gMC->Gspos("S08E",2*index  ,"S08B", xx, 0., kBframeWidth/2. - kNulocWidth/2, 0, "ONLY");
    }
    if (xx > xxmax4 && xx< xxmax) {
      gMC->Gspos("S07E",2*index-1,"SD7B", xx, 0.,-kBframeWidth/2.+ kNulocWidth/2, 0, "ONLY");
      gMC->Gspos("S07E",2*index  ,"SD7B", xx, 0., kBframeWidth/2.- kNulocWidth/2, 0, "ONLY");
      gMC->Gspos("S08E",2*index-1,"SD8B", xx, 0.,-kBframeWidth/2.+ kNulocWidth/2, 0, "ONLY");
      gMC->Gspos("S08E",2*index  ,"SD8B", xx, 0., kBframeWidth/2.- kNulocWidth/2, 0, "ONLY");
    }

    //
    //Geometry of the support pannel Verticla length 5.3m, horizontal length 2.6m, internal radius  dMotherInner o SC07 and SC08  (F. Orsini, Saclay)
    //Carbon fiber of 0.3 mm thick (2 layers) and a central layer of Nomex of 15mm thick. 
    Float_t dMotherInner =  AliMUONConstants::Rmin(3)-kRframeHeight; 
    Float_t nomexthickness = 1.5;
    Float_t carbonthickness = 0.03;
    Float_t supporthlength =  260.;  
    Float_t supportvlength =  530.;  
    // Generating the composite shape of the carbon and nomex pannels
    new TGeoBBox("shNomexBoxSt4",supporthlength/2., supportvlength/2. ,nomexthickness/2.+carbonthickness+3*kCableWidth);
    new TGeoBBox("shCarbonBoxSt4",supporthlength/2., supportvlength/2. ,carbonthickness/2.); 
    new TGeoTubeSeg("shNomexHoleSt4",0., dMotherInner, nomexthickness/2.+carbonthickness+3*kCableWidth+0.001, -90. ,90.);
    new TGeoTubeSeg("shCarbonHoleSt4",0., dMotherInner, carbonthickness/2.+0.001, -90. ,90.);
    TGeoTranslation* trHoleSt4 = new TGeoTranslation("trHoleSt4",-supporthlength/2.,0.,0.); 
    trHoleSt4->RegisterYourself();
    TGeoCompositeShape* shNomexSupportSt4  = new TGeoCompositeShape("shNomexSupportSt4","shNomexBoxSt4-shNomexHoleSt4:trHoleSt4");
    TGeoCompositeShape* shCarbonSupportSt4 = new TGeoCompositeShape("shCarbonSupportSt4","shCarbonBoxSt4-shCarbonHoleSt4:trHoleSt4");
 
   // Generating Nomex and Carbon pannel volumes
    TGeoVolume* voNomexSupportSt4   = new TGeoVolume("S07S", shNomexSupportSt4, kMedNomex);
    TGeoVolume* voCarbonSupportSt4  = new TGeoVolume("S07K", shCarbonSupportSt4, kMedCarbon);
    TGeoVolume* voNomexSupportSt4Ch8   = new TGeoVolume("S08S", shNomexSupportSt4, kMedNomex);
    TGeoVolume* voCarbonSupportSt4Ch8  = new TGeoVolume("S08K", shCarbonSupportSt4, kMedCarbon);
    TGeoTranslation* trCarbon1St4   = new TGeoTranslation("trCarbon1St4",0.,0., -(nomexthickness+carbonthickness)/2.);
    TGeoTranslation* trCarbon2St4   = new TGeoTranslation("trCarbon2St4",0.,0.,  (nomexthickness+carbonthickness)/2.);
    voNomexSupportSt4->AddNode(voCarbonSupportSt4,1,trCarbon1St4);
    voNomexSupportSt4->AddNode(voCarbonSupportSt4,2,trCarbon2St4);
    voNomexSupportSt4Ch8->AddNode(voCarbonSupportSt4Ch8,1,trCarbon1St4);
    voNomexSupportSt4Ch8->AddNode(voCarbonSupportSt4Ch8,2,trCarbon2St4);
 
    // Add readout cables
    gMC->Gsvolu("S07L","BOX",kCableMaterial,dum,0);
    gMC->Gsvolu("S08L","BOX",kCableMaterial,dum,0);

    ySlat41 = 0.;
    ySlat42 = 0.;
    Float_t lCableX = 0.;
    Float_t lCableY = 0.;
    Float_t lCableY8 = 0.;
    Float_t lCableZ = 0.;
    Float_t cablepar[3] = {supporthlength/2., kCableHeight/2., kCableWidth/2.};
    Float_t lCableDY = 0.;
    Float_t lCableDY8 = 0.;
    for (i = 0; i<kNslats4; i++){
      Int_t iCable = 1;
      Int_t cIndex = 0;
      Int_t cIndex8 = 0;
      ySlat41 += kYpos41[i];
      ySlat42 += kYpos42[i];

      lCableY = ySlat41;
      lCableY8 = ySlat42;

      // Cables going out from the start of slat
      if(kNPCB4[i]>=4 && i<kNslats4-2){ // Only if 4 or more pcb
	// First top cables
	cablepar[0] = (supporthlength-kXpos4[i])/2.;
	lCableX = kXpos4[i]/2.;
	lCableDY = (kYpos41[i+1]+kYpos41[i+2])/2.-cablepar[1];
	lCableDY8 = (kYpos42[i+1]+kYpos42[i+2])/2.-cablepar[1];
	lCableZ = TMath::Power(-1,i)*(nomexthickness/2.+carbonthickness+(-1+iCable++)*kCableWidth+kCableWidth/2.);
	gMC->Gsposp("S07L",10*i+cIndex++,"S07S",lCableX,lCableY+lCableDY,lCableZ,0,"ONLY",cablepar,3); 	
	gMC->Gsposp("S07L",10*i+cIndex++,"S07S",lCableX,-(lCableY+lCableDY),lCableZ,0,"ONLY",cablepar,3);
	gMC->Gsposp("S08L",10*i+cIndex8++,"S08S",lCableX,lCableY8+lCableDY8,lCableZ,0,"ONLY",cablepar,3); 	
	gMC->Gsposp("S08L",10*i+cIndex8++,"S08S",lCableX,-(lCableY8+lCableDY8),lCableZ,0,"ONLY",cablepar,3);
	// Then bottom cables
	if (i>0){
	  if (i==1) { // Rounded slat. Bottom cable starts at dMotherInner (beam pipe)
	    cablepar[0] = (supporthlength-kXpos4[i]-dMotherInner)/2.;
	    lCableX = (kXpos4[i]+dMotherInner)/2.;
	    lCableDY = (kYpos41[i]+kYpos41[i])/2.-cablepar[1];
	    lCableDY8 = (kYpos42[i]+kYpos42[i])/2.-cablepar[1];
	  }
	  if (i>=2) {
	    lCableDY = (kYpos41[i]+kYpos41[i-1])/2.-cablepar[1];
	    if ((lCableY-lCableDY)<(dMotherInner+cablepar[1])){
	      lCableDY = lCableY - dMotherInner - cablepar[1];
	    }
	    lCableDY8 = (kYpos42[i]+kYpos42[i-1])/2.-cablepar[1];
	    if ((lCableY8-lCableDY8)<(dMotherInner+cablepar[1])){
	      lCableDY8 = lCableY8 - dMotherInner - cablepar[1];
	    }
	  }
	  gMC->Gsposp("S07L",10*i+cIndex++,"S07S",lCableX,lCableY-lCableDY,lCableZ,0,"ONLY",cablepar,3); 
	  gMC->Gsposp("S07L",10*i+cIndex++,"S07S",lCableX,-(lCableY-lCableDY),lCableZ,0,"ONLY",cablepar,3); 
	  gMC->Gsposp("S08L",10*i+cIndex8++,"S08S",lCableX,lCableY8-lCableDY8,lCableZ,0,"ONLY",cablepar,3); 
	  gMC->Gsposp("S08L",10*i+cIndex8++,"S08S",lCableX,-(lCableY8-lCableDY8),lCableZ,0,"ONLY",cablepar,3); 
	}
      }

      // Rounded slats have an extra cable starting at second pcb
      if(i==1){ 
	// Only on top
	cablepar[0] = (supporthlength-kPcbLength-kVframeLength)/2.;
	lCableX = (kPcbLength+kVframeLength)/2.;
	lCableDY = (kYpos41[i+1]+kYpos41[i+2])/2.-cablepar[1]; // half way between 2 slats on same side
	lCableDY8 = (kYpos42[i+1]+kYpos42[i+2])/2.-cablepar[1]; // half way between 2 slats on same side
	lCableZ = TMath::Power(-1,i)*(nomexthickness/2.+carbonthickness+(-1+iCable++)*kCableWidth+kCableWidth/2.);
	gMC->Gsposp("S07L",10*i+cIndex++,"S07S",lCableX,lCableY+lCableDY,lCableZ,0,"ONLY",cablepar,3); 
	gMC->Gsposp("S07L",10*i+cIndex++,"S07S",lCableX,-(lCableY+lCableDY),lCableZ,0,"ONLY",cablepar,3); 
	gMC->Gsposp("S08L",10*i+cIndex8++,"S08S",lCableX,lCableY8+lCableDY8,lCableZ,0,"ONLY",cablepar,3); 
	gMC->Gsposp("S08L",10*i+cIndex8++,"S08S",lCableX,-(lCableY8+lCableDY8),lCableZ,0,"ONLY",cablepar,3); 
      }	

      // Cables going out from the end of the slats
      cablepar[0] = (supporthlength-(slatLength4[i]+kXpos4[i]+kDslatLength)+kVframeLength)/2.;
      lCableX = slatLength4[i]+kXpos4[i]-kVframeLength+kDslatLength+cablepar[0]-supporthlength/2.;
      if(i+1>=kNslats4 || i+2>=kNslats4){ // If no more higher slats, then use distance to lower slat
	lCableDY = kPcbHeight/2.+cablepar[1];
	lCableDY8 = lCableDY;
      }
      else {
	lCableDY = (kYpos41[i+1]+kYpos41[i+2])/2.-cablepar[1];
	lCableDY8 = (kYpos42[i+1]+kYpos42[i+2])/2.-cablepar[1];
      }
      lCableZ = TMath::Power(-1,i)*(nomexthickness/2.+carbonthickness+(-1+iCable++)*kCableWidth+kCableWidth/2.);
      gMC->Gsposp("S07L",10*i+cIndex++,"S07S",lCableX,lCableY+lCableDY,lCableZ,0,"ONLY",cablepar,3); 
      gMC->Gsposp("S07L",10*i+cIndex++,"S07S",lCableX,-(lCableY+lCableDY),lCableZ,0,"ONLY",cablepar,3); 
      gMC->Gsposp("S08L",10*i+cIndex8++,"S08S",lCableX,lCableY8+lCableDY8,lCableZ,0,"ONLY",cablepar,3); 
      gMC->Gsposp("S08L",10*i+cIndex8++,"S08S",lCableX,-(lCableY8+lCableDY8),lCableZ,0,"ONLY",cablepar,3); 
      // Then bottom cables
      if(i>0){
	if (i==1) {
	  lCableDY = (kYpos41[i]+kYpos41[i])/2.-cablepar[1];
	  lCableDY8 = (kYpos42[i]+kYpos42[i])/2.-cablepar[1];
	}
	else{
	  lCableDY = (kYpos41[i]+kYpos41[i-1])/2.-cablepar[1]; // half way between 2 slats on same side
	  if ((lCableY-lCableDY)<(dMotherInner+cablepar[1])){
	    lCableDY = lCableY - dMotherInner - cablepar[1];
	  }
	  lCableDY8 = (kYpos42[i]+kYpos42[i-1])/2.-cablepar[1]; // half way between 2 slats on same side
	  if ((lCableY8-lCableDY8)<(dMotherInner+cablepar[1])){
	    lCableDY8 = lCableY8 - dMotherInner - cablepar[1];
	  }
	}
	gMC->Gsposp("S07L",10*i+cIndex++,"S07S",lCableX,lCableY-lCableDY,lCableZ,0,"ONLY",cablepar,3); 
	gMC->Gsposp("S07L",10*i+cIndex++,"S07S",lCableX,-(lCableY-lCableDY),lCableZ,0,"ONLY",cablepar,3); 
	gMC->Gsposp("S08L",10*i+cIndex8++,"S08S",lCableX,lCableY8-lCableDY8,lCableZ,0,"ONLY",cablepar,3); 
	gMC->Gsposp("S08L",10*i+cIndex8++,"S08S",lCableX,-(lCableY8-lCableDY8),lCableZ,0,"ONLY",cablepar,3); 
      }	
    }
    
    Float_t dzCh7  = dzCh;
    TGeoTranslation* trSupport1St4   = new TGeoTranslation("trSupport1St4", supporthlength/2., 0. , dzCh7);
    TGeoRotation*    roSupportSt4    = new TGeoRotation("roSupportSt4",90.,180.,-90.);
    TGeoCombiTrans*  coSupport2St4   = new TGeoCombiTrans(-supporthlength/2., 0., -dzCh7, roSupportSt4); 
    GetEnvelopes(9)->AddEnvelope("S07S", 0, 1, *trSupport1St4);  
    GetEnvelopes(8)->AddEnvelope("S07S", 0, 2, *coSupport2St4);  
    GetEnvelopes(11)->AddEnvelope("S08S", 0, 1, *trSupport1St4);   
    GetEnvelopes(10)->AddEnvelope("S08S", 0, 2, *coSupport2St4);

    // End of pannel support geometry    

    // cout << "Geometry for Station 4...... done" << endl;

  }
    
  if (fStations[4]) {
      

    // //********************************************************************
    // //                            Station 5                             **
    // //********************************************************************
    // Mother volume for each chamber in St4 is an envelop (or assembly)
    // There is one assembly mother per half a chamber  called SC09I, SC09O, SC10I and SC10O 
    // Same volume name definitions as in St3
    
    const Int_t   kNslats5         = 7;  // number of slats per quadrant
    const Int_t   kNPCB5[kNslats5] = {5, 6, 6, 6, 5, 4, 3}; // n PCB per slat
    const Float_t kXpos5[kNslats5] = {38.75, 0., 0., 0., 0., 0., 0.}; // J.C. Correct value  
    const Float_t kYpos5[kNslats5] = {0., 38.2, 37.9, 37.6, 37.3, 37.05, 36.75};
    Float_t slatLength5[kNslats5]; 

    Float_t rPhi1 = TMath::RadToDeg()*(TMath::ASin((kYpos5[1]-hFramepar[1])/(AliMUONConstants::Rmin(4))));
    Float_t rPhi2 = TMath::RadToDeg()*(TMath::ACos(-vFramepar[0]/(AliMUONConstants::Rmin(4)-kRframeLength)));
    Float_t rFramepar5[5] = { AliMUONConstants::Rmin(4)-kRframeLength, AliMUONConstants::Rmin(4), kRframeWidth, rPhi1, rPhi2}; 
    Float_t vrFrameHeight = hFramepar[1]+kYpos5[1]-AliMUONConstants::Rmin(4)+kRframeLength;

    char idSlatCh9[6];
    char idSlatCh10[6];
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

      sprintf(idSlatCh9,"SLE%d",kNslats5-1+i);
      detElemId = 913 - (i + kNslats5-1-6);
      moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
      if (detElemId % 2 == 0) {
	if (detElemId == 912) // Round slat, new rotation due to mapping convention
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh9, detElemId, true, TGeoTranslation(xSlat5, ySlat5, -zSlat5 + dzCh5),
					      TGeoRotation("rot1",90,180+angle,90,90+angle,180,0) );
	else
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh9, detElemId, true, TGeoTranslation(xSlat5, ySlat5, -zSlat5 + dzCh5),
					      TGeoRotation("rot1",90,angle,90,90+angle,0,0) );
      }
      else
	GetEnvelopes(moduleId)->AddEnvelope(idSlatCh9, detElemId, true, TGeoTranslation(xSlat5, ySlat5, -zSlat5 + dzCh5),
					    TGeoRotation("rot1",90,angle,90,270+angle,180,0) );
      sprintf(idSlatCh9,"SLE%d",3*kNslats5-2+i);
      detElemId = 900 + (i + kNslats5-1-6);
      moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
      if (detElemId % 2 == 1) {
	if (detElemId == 901) // Round slat, new rotation due to mapping convention
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh9, detElemId, true, TGeoTranslation(-xSlat5, ySlat5, zSlat5 - dzCh5),
					      TGeoRotation("rot2",90,angle,90,90+angle,0,0) );
	else
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh9, detElemId, true, TGeoTranslation(-xSlat5, ySlat5, zSlat5 - dzCh5),
					      TGeoRotation("rot2",90,180+angle,90,90+angle,180,0) );
      }
      else
	GetEnvelopes(moduleId)->AddEnvelope(idSlatCh9, detElemId, true, TGeoTranslation(-xSlat5, ySlat5, zSlat5 - dzCh5),
					    TGeoRotation("rot2",90,180+angle,90,270+angle,0,0) );      

      if (i > 0) { 
	sprintf(idSlatCh9,"SLE%d",kNslats5-1-i);
	detElemId = 913 + (i + kNslats5-1-6);
	moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
	if (detElemId % 2 == 0) {
	  if (detElemId == 914) // Round slat, new rotation due to mapping convention
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh9, detElemId, true, TGeoTranslation(xSlat5, -ySlat5, -zSlat5 + dzCh5),
					      TGeoRotation("rot3",90,180+angle,90,270+angle,0,0) );
	  else
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh9, detElemId, true, TGeoTranslation(xSlat5, -ySlat5, -zSlat5 + dzCh5),
						TGeoRotation("rot3",90,angle,90,90+angle,0,0) );
	}
	else
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh9, detElemId, true, TGeoTranslation(xSlat5, -ySlat5, -zSlat5 + dzCh5),
				     TGeoRotation("rot3",90,angle,90,270+angle,180,0) );

	sprintf(idSlatCh9,"SLE%d",3*kNslats5-2-i);
	detElemId = 926 - (i + kNslats5-1-6);
	moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
	if (detElemId % 2 == 1) {
	  if (detElemId == 925) // Round slat, new rotation due to mapping convention
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh9, detElemId, true, TGeoTranslation(-xSlat5, -ySlat5, zSlat5 - dzCh5),
						TGeoRotation("rot4",90,angle,90,270+angle,180,0) );
	  else
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh9, detElemId, true, TGeoTranslation(-xSlat5, -ySlat5, zSlat5 - dzCh5),
						TGeoRotation("rot4",90,180+angle,90,90+angle,180,0)  );
	}
	else
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh9, detElemId, true, TGeoTranslation(-xSlat5, -ySlat5, zSlat5 - dzCh5),
					      TGeoRotation("rot4",90,180+angle,90,270+angle,0,0)  );   
      }

      sprintf(idSlatCh10,"SLF%d",kNslats5-1+i);
      detElemId = 1013 - (i + kNslats5-1-6);
      moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
      if (detElemId % 2 == 0) {
	if (detElemId == 1012) // Round slat, new rotation due to mapping convention
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh10, detElemId, true, TGeoTranslation(xSlat5, ySlat5, -zSlat5 + dzCh5),
					      TGeoRotation("rot5",90,180+angle,90,90+angle,180,0) );
	else
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh10, detElemId, true, TGeoTranslation(xSlat5, ySlat5, -zSlat5 + dzCh5),
					      TGeoRotation("rot5",90,angle,90,90+angle,0,0) );
      }
      else
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh10, detElemId, true, TGeoTranslation(xSlat5, ySlat5, -zSlat5 + dzCh5),
					      TGeoRotation("rot5",90,angle,90,270+angle,180,0) );

      sprintf(idSlatCh10,"SLF%d",3*kNslats5-2+i);
      detElemId = 1000 + (i + kNslats5-1-6);
      moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
      if (detElemId % 2 == 1) {
	if (detElemId == 1001) // Round slat, new rotation due to mapping convention
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh10, detElemId, true, TGeoTranslation(-xSlat5, ySlat5, zSlat5 - dzCh5),
					      TGeoRotation("rot6",90,angle,90,90+angle,0,0) ); 
	else
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh10, detElemId, true, TGeoTranslation(-xSlat5, ySlat5, zSlat5 - dzCh5),
					      TGeoRotation("rot6",90,180+angle,90,90+angle,180,0) ); 
      }
      else
	GetEnvelopes(moduleId)->AddEnvelope(idSlatCh10, detElemId, true, TGeoTranslation(-xSlat5, ySlat5, zSlat5 - dzCh5),
					    TGeoRotation("rot6",90,180+angle,90,270+angle,0,0) );

      if (i > 0) { 
	sprintf(idSlatCh10,"SLF%d",kNslats5-1-i);
	detElemId = 1013 + (i + kNslats5-1-6);
        moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
	if (detElemId % 2 == 0) {  
	  if (detElemId == 1014) // Round slat, new rotation due to mapping convention
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh10, detElemId, true, TGeoTranslation(xSlat5, -ySlat5, -zSlat5 + dzCh5),
						TGeoRotation("rot7",90,180+angle,90,270+angle,0,0) );
	  else
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh10, detElemId, true, TGeoTranslation(xSlat5, -ySlat5, -zSlat5 + dzCh5),
						TGeoRotation("rot7",90,angle,90,90+angle,0,0) );
	}
	else
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh10, detElemId, true, TGeoTranslation(xSlat5, -ySlat5, -zSlat5 + dzCh5),
					      TGeoRotation("rot7",90,angle,90,270+angle,180,0) );

	sprintf(idSlatCh10,"SLF%d",3*kNslats5-2-i);
	detElemId = 1026 - (i + kNslats5-1-6);
        moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
	if (detElemId % 2 == 1) {
	  if (detElemId == 1025) // Round slat, new rotation due to mapping convention
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh10, detElemId, true, TGeoTranslation(-xSlat5, -ySlat5, zSlat5 - dzCh5),
						TGeoRotation("rot8",90,angle,90,270+angle,180,0) );
	  else
	    GetEnvelopes(moduleId)->AddEnvelope(idSlatCh10, detElemId, true, TGeoTranslation(-xSlat5, -ySlat5, zSlat5 - dzCh5),
						TGeoRotation("rot8",90,180+angle,90,90+angle,180,0) );
	}
	else
	  GetEnvelopes(moduleId)->AddEnvelope(idSlatCh10, detElemId, true, TGeoTranslation(-xSlat5, -ySlat5, zSlat5 - dzCh5),
					      TGeoRotation("rot8",90,180+angle,90,270+angle,0,0) ); 
      }
    }

    // create the panel volume 
 
    gMC->Gsvolu("S09C","BOX",kCarbonMaterial,panelpar,3);
    gMC->Gsvolu("SD9C","BOX",kCarbonMaterial,panelpar,3);
    gMC->Gsvolu("S10C","BOX",kCarbonMaterial,panelpar,3);
    gMC->Gsvolu("SD0C","BOX",kCarbonMaterial,panelpar,3);

    // create the nomex volume 

    gMC->Gsvolu("S09N","BOX",kNomexMaterial,nomexpar,3);
    gMC->Gsvolu("SD9N","BOX",kNomexMaterial,nomexpar,3);
    gMC->Gsvolu("S10N","BOX",kNomexMaterial,nomexpar,3);
    gMC->Gsvolu("SD0N","BOX",kNomexMaterial,nomexpar,3);


    // create the nomex volume (bulk)

    gMC->Gsvolu("S09X","BOX",kNomexBMaterial,nomexbpar,3);
    gMC->Gsvolu("SD9X","BOX",kNomexBMaterial,nomexbpar,3);
    gMC->Gsvolu("S10X","BOX",kNomexBMaterial,nomexbpar,3);
    gMC->Gsvolu("SD0X","BOX",kNomexBMaterial,nomexbpar,3);

    // create the insulating material volume 

    gMC->Gsvolu("S09I","BOX",kInsuMaterial,insupar,3);
    gMC->Gsvolu("SD9I","BOX",kInsuMaterial,insupar,3);
    gMC->Gsvolu("S10I","BOX",kInsuMaterial,insupar,3);
    gMC->Gsvolu("SD0I","BOX",kInsuMaterial,insupar,3);

    // create the PCB volume 

    gMC->Gsvolu("S09P","BOX",kPcbMaterial,pcbpar,3);
    gMC->Gsvolu("SD9P","BOX",kPcbMaterial,pcbpar,3);
    gMC->Gsvolu("S10P","BOX",kPcbMaterial,pcbpar,3);
    gMC->Gsvolu("SD0P","BOX",kPcbMaterial,pcbpar,3);
 
    // create the sensitive volumes,

    gMC->Gsvolu("S09G","BOX",kSensMaterial,dum,0);
    gMC->Gsvolu("SD9G","BOX",kSensMaterial,senspar,3);
    gMC->Gsvolu("S10G","BOX",kSensMaterial,dum,0);
    gMC->Gsvolu("SD0G","BOX",kSensMaterial,senspar,3);

    // create the vertical frame volume 

    gMC->Gsvolu("S09V","BOX",kVframeMaterial,vFramepar,3);
    gMC->Gsvolu("S10V","BOX",kVframeMaterial,vFramepar,3);

    // create the rounded vertical frame volume 

    gMC->Gsvolu("SD9D","TUBS",kRframeMaterial,rFramepar5,5);
    gMC->Gsvolu("SD0D","TUBS",kRframeMaterial,rFramepar5,5);

    // create the horizontal frame volume 

    gMC->Gsvolu("S09H","BOX",kHframeMaterial,hFramepar,3);
    gMC->Gsvolu("SD9H","BOX",kHframeMaterial,hFramepar,3);
    gMC->Gsvolu("S10H","BOX",kHframeMaterial,hFramepar,3);
    gMC->Gsvolu("SD0H","BOX",kHframeMaterial,hFramepar,3);

    // create the horizontal border volume 

    gMC->Gsvolu("S09B","BOX",kBframeMaterial,bFramepar,3);
    gMC->Gsvolu("SD9B","BOX",kBframeMaterial,bFramepar,3);
    gMC->Gsvolu("S10B","BOX",kBframeMaterial,bFramepar,3);
    gMC->Gsvolu("SD0B","BOX",kBframeMaterial,bFramepar,3);

    // Replace the volume shape with a composite shape
    // with substracted overlap with beam shield     
    if ( gMC->IsRootGeometrySupported() ) { 
	
      // Get shape
      Int_t nSlatType = 1;
      Int_t nVol = 8;
      const char* slatType = "D"; // D: Rounde slat
      const char* volLetter = "CNXIPHBG";
      TString volName;
      TString compName;
      TString csName;
      TGeoVolume *mVol = 0x0;
      // Beam shield recess
      new TGeoTube("tube5Cut", 0., AliMUONConstants::Rmin(4), kSlatWidth/2.+0.001);
      TObjArray rounded5Slat(nSlatType*((nVol+1)*2));	
      // Displacement
      TGeoTranslation* trDTube5 = new TGeoTranslation("trDTube5", -(kPcbLength+kVframeLength)/2., -kYpos5[1], 0.);
      trDTube5->RegisterYourself();
      TGeoTranslation* trDBTube5 = new TGeoTranslation("trDBTube5", 0., ( kPcbHeight - kBframeHeight ) / 2., 0.);
      trDBTube5->Add(trDTube5);
      trDBTube5->RegisterYourself();

      TObjArray composite5(nSlatType*((nVol+1)*2));
      new TGeoBBox("box5DCut",(kPcbLength+kVframeLength)/2., hFramepar[1], vFramepar[2]+0.001);
      // Displacement
      TGeoTranslation* trDBox5 = new TGeoTranslation("trDBox5",kPcbLength/2., kYpos5[1], 0.);
      trDBox5->RegisterYourself();
      
      TGeoBBox *box5Vframe = new TGeoBBox("box5Vframe",vFramepar[0],vrFrameHeight/2., vFramepar[2]);
      TGeoTranslation* trVBox5 = new TGeoTranslation("trVBox5", 0., AliMUONConstants::Rmin(4)-kRframeLength + box5Vframe->GetDY(), 0.);
      trVBox5->RegisterYourself();
      
      for(int iCh=9; iCh<=10; iCh++){
	for (int iSlatType = 0; iSlatType<nSlatType; iSlatType++) {
	  for (int iVol = 0; iVol<nVol; iVol++){
	    Int_t lIndex = (iCh-9)*(nSlatType*(nVol+1))+iSlatType*(nVol+1)+iVol;
	    volName=Form("S%c%d%c",slatType[iSlatType],iCh%10,volLetter[iVol]);
	    mVol = gGeoManager->FindVolumeFast(volName);
	    if ( !mVol ) {
	      AliErrorStream() 
		<< "Slat volume " << volName << " not found" << endl;	 
	    }
	    else {
	      rounded5Slat[lIndex] = mVol->GetShape();
	      csName=Form("rounded5Slat%c%d%c",slatType[iSlatType],iCh%10,volLetter[iVol]);
	      ((TGeoShape*)rounded5Slat[lIndex])->SetName(csName);  
	      
	      // Composite shape
	      TString compOperation(csName);
	      compOperation+="-tube5Cut:tr";
	      compOperation+=slatType[iSlatType];
	      if (strstr(volName,"B")){
		compOperation+="B";
	      }
	      compOperation+="Tube5";
	      compName=Form("composite5%c%d%c",slatType[iSlatType],iCh,volLetter[iVol]);
	      composite5[lIndex] = new TGeoCompositeShape(compName, compOperation.Data()); 
	      
	      // Reset shape to volume      
	      mVol->SetShape((TGeoShape*)composite5[lIndex]);
	    }
	  }

	  // For rounded spacer
	  Int_t lIndex = (iCh-9)*(nSlatType*(nVol+1))+iSlatType*(nVol+1)+nVol;
	  volName=Form("S%c%dD",slatType[iSlatType],iCh%10);
	  mVol = gGeoManager->FindVolumeFast(volName);
	  if ( !mVol ) {
	    AliErrorStream() 
	      << "Slat volume " << volName << " not found" << endl;	 
	  }
	  else {
	    rounded5Slat[lIndex] = mVol->GetShape();
	    csName=Form("rounded5Slat%c%dD",slatType[iSlatType],iCh%10);
	    ((TGeoShape*)rounded5Slat[lIndex])->SetName(csName);	 	  
	    
	    // Composite shape
	    TString compOperation(csName);
	    if (strstr(volName,"SD")){
	      compOperation.Prepend("(");
	      compOperation+="+box5Vframe:trVBox5)*box5DCut:trDBox5";
	    }
	    compName=Form("composite5%c%dD",slatType[iSlatType],iCh%10);
	    composite5[lIndex] = new TGeoCompositeShape(compName, compOperation.Data()); 	      
	    // Reset shape to volume      
	    mVol->SetShape((TGeoShape*)composite5[lIndex]);
	  }
	}
      }
    }
    
    index = 0; 
    for (i = 0; i < kNslats5; i++){
      for (Int_t quadrant = 1; quadrant <= 4; quadrant++) {

	if (i == 0 && quadrant == 2) continue;
	if (i == 0 && quadrant == 4) continue;

	sprintf(idSlatCh9,"SLE%d",ConvertSlatNum(i,quadrant,kNslats5-1));
	sprintf(idSlatCh10,"SLF%d",ConvertSlatNum(i,quadrant,kNslats5-1));
	Int_t moduleSlatCh9 = GetModuleId(idSlatCh9);
	Int_t moduleSlatCh10 = GetModuleId(idSlatCh10);
	Float_t xvFrame  = (slatLength5[i] - kVframeLength)/2.; // ok

	// position the vertical frames (spacers)
	if (i != 1) { 
	  GetEnvelopes(moduleSlatCh9)->AddEnvelopeConstituent("S09V", idSlatCh9, (2*i+1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
	  GetEnvelopes(moduleSlatCh9)->AddEnvelopeConstituent("S09V", idSlatCh9, (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));
	  GetEnvelopes(moduleSlatCh10)->AddEnvelopeConstituent("S10V", idSlatCh10, (2*i+1)*10+quadrant,TGeoTranslation(xvFrame,0.,0.));
	  GetEnvelopes(moduleSlatCh10)->AddEnvelopeConstituent("S10V", idSlatCh10, (2*i)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.));
	} else {  // Vertical and Rounded+Vertical spacer - Different rotation due to new mapping convention
	  GetEnvelopes(moduleSlatCh9)->AddEnvelopeConstituent("S09V", idSlatCh9, (2*i+1)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	  GetEnvelopes(moduleSlatCh9)->AddEnvelopeConstituent("SD9D", idSlatCh9, (2*i)*10+quadrant,TGeoTranslation(xvFrame,-kYpos5[1],0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	  GetEnvelopes(moduleSlatCh10)->AddEnvelopeConstituent("S10V", idSlatCh10, (2*i+1)*10+quadrant,TGeoTranslation(-xvFrame,0.,0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	  GetEnvelopes(moduleSlatCh10)->AddEnvelopeConstituent("SD0D", idSlatCh10, (2*i)*10+quadrant,TGeoTranslation(xvFrame,-kYpos5[1],0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	}

	// position the panels and the insulating material 
	for (j = 0; j < kNPCB5[i]; j++){
	  index++;
	  xx = kSensLength * (-kNPCB5[i]/2.+j+.5); 

	  Float_t zPanel = spar[2] - nomexbpar[2]; 
	  if (i==1) { // Different rotation due to new mapping convention
	    if (j==0) { // Rounded pcb of rounded slat 
	    GetEnvelopes(moduleSlatCh9)->AddEnvelopeConstituent("SD9X", idSlatCh9, 2*index-1,TGeoTranslation(-xx,0.,zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	    GetEnvelopes(moduleSlatCh9)->AddEnvelopeConstituent("SD9X", idSlatCh9, 2*index,TGeoTranslation(-xx,0.,-zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	    GetEnvelopes(moduleSlatCh9)->AddEnvelopeConstituent("SD9I", idSlatCh9, index,TGeoTranslation(-xx,0.,0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	    GetEnvelopes(moduleSlatCh10)->AddEnvelopeConstituent("SD0X", idSlatCh10, 2*index-1,TGeoTranslation(-xx,0.,zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	    GetEnvelopes(moduleSlatCh10)->AddEnvelopeConstituent("SD0X", idSlatCh10, 2*index,TGeoTranslation(-xx,0.,-zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	    GetEnvelopes(moduleSlatCh10)->AddEnvelopeConstituent("SD0I", idSlatCh10, index,TGeoTranslation(-xx,0.,0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	    } else { 
	      GetEnvelopes(moduleSlatCh9)->AddEnvelopeConstituent("S09X", idSlatCh9, 2*index-1,TGeoTranslation(-xx,0.,zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	      GetEnvelopes(moduleSlatCh9)->AddEnvelopeConstituent("S09X", idSlatCh9, 2*index,TGeoTranslation(-xx,0.,-zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	      GetEnvelopes(moduleSlatCh9)->AddEnvelopeConstituent("S09I", idSlatCh9, index,TGeoTranslation(-xx,0.,0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));	    
	      GetEnvelopes(moduleSlatCh10)->AddEnvelopeConstituent("S10X", idSlatCh10, 2*index-1,TGeoTranslation(-xx,0.,zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	      GetEnvelopes(moduleSlatCh10)->AddEnvelopeConstituent("S10X", idSlatCh10, 2*index,TGeoTranslation(-xx,0.,-zPanel),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	      GetEnvelopes(moduleSlatCh10)->AddEnvelopeConstituent("S10I", idSlatCh10, index,TGeoTranslation(-xx,0.,0.),TGeoRotation("rotAbX",90,180+angle,90,90+angle,180,0));
	    }
	  } else { 
	    GetEnvelopes(moduleSlatCh9)->AddEnvelopeConstituent("S09X", idSlatCh9, 2*index-1,TGeoTranslation(xx,0.,zPanel));
	    GetEnvelopes(moduleSlatCh9)->AddEnvelopeConstituent("S09X", idSlatCh9, 2*index,TGeoTranslation(xx,0.,-zPanel));
	    GetEnvelopes(moduleSlatCh9)->AddEnvelopeConstituent("S09I", idSlatCh9, index,TGeoTranslation(xx,0.,0.));	    
	    GetEnvelopes(moduleSlatCh10)->AddEnvelopeConstituent("S10X", idSlatCh10, 2*index-1,TGeoTranslation(xx,0.,zPanel));
	    GetEnvelopes(moduleSlatCh10)->AddEnvelopeConstituent("S10X", idSlatCh10, 2*index,TGeoTranslation(xx,0.,-zPanel));
	    GetEnvelopes(moduleSlatCh10)->AddEnvelopeConstituent("S10I", idSlatCh10, index,TGeoTranslation(xx,0.,0.));
	  }
	} 
      }
    }

    // position the nomex volume inside the panel volume
    gMC->Gspos("S09N",1,"S09C",0.,0.,0.,0,"ONLY"); 
    gMC->Gspos("SD9N",1,"SD9C",0.,0.,0.,0,"ONLY"); 
    gMC->Gspos("S10N",1,"S10C",0.,0.,0.,0,"ONLY"); 
    gMC->Gspos("SD0N",1,"SD0C",0.,0.,0.,0,"ONLY"); 

    // position panel  volume inside the bulk nomex material volume
    gMC->Gspos("S09C",1,"S09X",0.,0.,kNomexBWidth/2.,0,"ONLY"); 
    gMC->Gspos("SD9C",1,"SD9X",0.,0.,kNomexBWidth/2.,0,"ONLY"); 
    gMC->Gspos("S10C",1,"S10X",0.,0.,kNomexBWidth/2.,0,"ONLY"); 
    gMC->Gspos("SD0C",1,"SD0X",0.,0.,kNomexBWidth/2.,0,"ONLY"); 

    // position the PCB volume inside the insulating material volume
    gMC->Gspos("S09P",1,"S09I",0.,0.,0.,0,"ONLY"); 
    gMC->Gspos("SD9P",1,"SD9I",0.,0.,0.,0,"ONLY"); 
    gMC->Gspos("S10P",1,"S10I",0.,0.,0.,0,"ONLY"); 
    gMC->Gspos("SD0P",1,"SD0I",0.,0.,0.,0,"ONLY"); 

    // position the horizontal frame volume inside the PCB volume
    gMC->Gspos("S09H",1,"S09P",0.,0.,0.,0,"ONLY"); 
    gMC->Gspos("SD9H",1,"SD9P",0.,0.,0.,0,"ONLY"); 
    gMC->Gspos("S10H",1,"S10P",0.,0.,0.,0,"ONLY"); 
    gMC->Gspos("SD0H",1,"SD0P",0.,0.,0.,0,"ONLY"); 

    // position the sensitive volume inside the horizontal frame volume
    gMC->Gsposp("S09G",1,"S09H",0.,0.,0.,0,"ONLY",senspar,3); 
    gMC->Gspos("SD9G",1,"SD9H",0.,0.,0.,0,"ONLY"); 
    gMC->Gsposp("S10G",1,"S10H",0.,0.,0.,0,"ONLY",senspar,3); 
    gMC->Gspos("SD0G",1,"SD0H",0.,0.,0.,0,"ONLY"); 

    // position the border volumes inside the PCB volume
    Float_t yborder = ( kPcbHeight - kBframeHeight ) / 2.; 
    gMC->Gspos("S09B",1,"S09P",0., yborder,0.,0,"ONLY"); 
    gMC->Gspos("S09B",2,"S09P",0.,-yborder,0.,0,"ONLY"); 
    gMC->Gspos("S09B",1,"SD9P",0., yborder,0.,0,"ONLY"); 
    gMC->Gspos("SD9B",1,"SD9P",0.,-yborder,0.,0,"ONLY"); 
    gMC->Gspos("S10B",1,"S10P",0., yborder,0.,0,"ONLY"); 
    gMC->Gspos("S10B",2,"S10P",0.,-yborder,0.,0,"ONLY"); 
    gMC->Gspos("S10B",1,"SD0P",0., yborder,0.,0,"ONLY"); 
    gMC->Gspos("SD0B",1,"SD0P",0.,-yborder,0.,0,"ONLY"); 

    //      // create the NULOC volume and position it in the horizontal frame

    gMC->Gsvolu("S09E","BOX",kNulocMaterial,nulocpar,3);
    gMC->Gsvolu("S10E","BOX",kNulocMaterial,nulocpar,3);
    index = 0;
    Float_t rPhi3 = TMath::ASin((kYpos5[1]-kPcbHeight/2.)/AliMUONConstants::Rmin(4));
    Float_t xxmax4 = (AliMUONConstants::Rmin(4)*TMath::Cos(rPhi3)-kVframeLength/2.) - (kBframeLength - kNulocLength)/2.;
    for (xx = -xxmax; xx <= xxmax; xx += 2*kNulocLength) { 
      index++; 
      gMC->Gspos("S09E",2*index-1,"S09B", xx, 0.,-kBframeWidth/2. + kNulocWidth/2, 0, "ONLY");
      gMC->Gspos("S09E",2*index  ,"S09B", xx, 0., kBframeWidth/2. - kNulocWidth/2, 0, "ONLY");
      gMC->Gspos("S10E",2*index-1,"S10B", xx, 0.,-kBframeWidth/2. + kNulocWidth/2, 0, "ONLY");
      gMC->Gspos("S10E",2*index  ,"S10B", xx, 0., kBframeWidth/2. - kNulocWidth/2, 0, "ONLY");
    }
    if (xx > xxmax4 && xx< xxmax) {
      gMC->Gspos("S09E",2*index-1,"SD9B", xx, 0.,-kBframeWidth/2.+ kNulocWidth/2, 0, "ONLY");
      gMC->Gspos("S09E",2*index  ,"SD9B", xx, 0., kBframeWidth/2.- kNulocWidth/2, 0, "ONLY");
      gMC->Gspos("S10E",2*index-1,"SD0B", xx, 0.,-kBframeWidth/2.+ kNulocWidth/2, 0, "ONLY");
      gMC->Gspos("S10E",2*index  ,"SD0B", xx, 0., kBframeWidth/2.- kNulocWidth/2, 0, "ONLY");
    }

    //    
    //Geometry of the support pannel Verticla length 5.7m, horizontal length 2.6m, internal radius  dMotherInner o SC09 and SC10  (F. Orsini, Saclay)
    //Carbon fiber of 0.3 mm thick (2 layers) and a central layer of Nomex of 15mm thick. 
    Float_t dMotherInner =  AliMUONConstants::Rmin(4)-kRframeHeight; 
    Float_t nomexthickness = 1.5;
    Float_t carbonthickness = 0.03;
    Float_t supporthlength =  260.;  
    Float_t supportvlength =  570.;  
    // Generating the composite shape of the carbon and nomex pannels
    new TGeoBBox("shNomexBoxSt5",supporthlength/2., supportvlength/2. ,nomexthickness/2.+carbonthickness+3*kCableWidth);
    new TGeoBBox("shCarbonBoxSt5",supporthlength/2., supportvlength/2. ,carbonthickness/2.); 
    new TGeoTubeSeg("shNomexHoleSt5",0., dMotherInner, nomexthickness/2.+carbonthickness+3*kCableWidth+0.001, -90. ,90.);
    new TGeoTubeSeg("shCarbonHoleSt5",0., dMotherInner, carbonthickness/2.+0.001, -90. ,90.);
    TGeoTranslation* trHoleSt5 = new TGeoTranslation("trHoleSt5",-supporthlength/2.,0.,0.); 
    trHoleSt5->RegisterYourself();
    TGeoCompositeShape* shNomexSupportSt5  = new TGeoCompositeShape("shNomexSupportSt5","shNomexBoxSt5-shNomexHoleSt5:trHoleSt5");
    TGeoCompositeShape* shCarbonSupportSt5 = new TGeoCompositeShape("shCarbonSupportSt5","shCarbonBoxSt5-shCarbonHoleSt5:trHoleSt5");

   // Generating Nomex and Carbon pannel volumes
    TGeoVolume* voNomexSupportSt5  = new TGeoVolume("S09S", shNomexSupportSt5, kMedNomex);
    TGeoVolume* voCarbonSupportSt5 = new TGeoVolume("S09K", shCarbonSupportSt5, kMedCarbon);
    TGeoTranslation* trCarbon1St5   = new TGeoTranslation("trCarbon1St5",0.,0., -(nomexthickness+carbonthickness)/2.);
    TGeoTranslation* trCarbon2St5   = new TGeoTranslation("trCarbon2St5",0.,0.,  (nomexthickness+carbonthickness)/2.);
    voNomexSupportSt5->AddNode(voCarbonSupportSt5,1,trCarbon1St5);
    voNomexSupportSt5->AddNode(voCarbonSupportSt5,2,trCarbon2St5);

    // Add readout cables
    gMC->Gsvolu("S09L","BOX",kCableMaterial,dum,0);

    ySlat5 = 0.;
    Float_t lCableX = 0.;
    Float_t lCableY = 0.;
    Float_t lCableZ = 0.;
    Float_t cablepar[3] = {supporthlength/2., kCableHeight/2., kCableWidth/2.};
    Float_t lCableDY = 0.;
    for (i = 0; i<kNslats5; i++){
      Int_t iCable = 1;
      Int_t cIndex = 0;
      ySlat5 += kYpos5[i];

      lCableY = ySlat5;

      // Cables going out from the start of slat
      if(kNPCB5[i]>=4){ // Only if 4 or more pcb
	// First top cables
	cablepar[0] = (supporthlength-kXpos5[i])/2.;
	lCableX = kXpos5[i]/2.;
	if(i+1>=kNslats5 || i+2>=kNslats5){ // If no more higher slats, then use distance to lower slat
	  lCableDY = (kYpos5[i]+kYpos5[i-1])/2.-cablepar[1];
	}
	else {
	  lCableDY = (kYpos5[i+1]+kYpos5[i+2])/2.-cablepar[1];
	}
	lCableZ = TMath::Power(-1,i)*(nomexthickness/2.+carbonthickness+(-1+iCable++)*kCableWidth+kCableWidth/2.);
	gMC->Gsposp("S09L",10*i+cIndex++,"S09S",lCableX,lCableY+lCableDY,lCableZ,0,"ONLY",cablepar,3); 	
	gMC->Gsposp("S09L",10*i+cIndex++,"S09S",lCableX,-(lCableY+lCableDY),lCableZ,0,"ONLY",cablepar,3);
	// Then bottom cables
	if (i>0) {
	  if (i==1) { // Rounded slat. Bottom cable starts at dMotherInner (beam pipe)
	    cablepar[0] = (supporthlength-kXpos5[i]-dMotherInner)/2.;
	    lCableX = (kXpos5[i]+dMotherInner)/2.;
	    lCableDY = (kYpos5[i]+kYpos5[i])/2.-cablepar[1];
	  }
	  else {
	    lCableDY = (kYpos5[i]+kYpos5[i-1])/2.-cablepar[1];
	    if ((lCableY-lCableDY)<(dMotherInner+cablepar[1])){
	      lCableDY = lCableY - dMotherInner - cablepar[1];
	    }
	  }
	  gMC->Gsposp("S09L",10*i+cIndex++,"S09S",lCableX,lCableY-lCableDY,lCableZ,0,"ONLY",cablepar,3); 
	  gMC->Gsposp("S09L",10*i+cIndex++,"S09S",lCableX,-(lCableY-lCableDY),lCableZ,0,"ONLY",cablepar,3); 
	}
      }
      
      // Rounded slats have an extra cable starting at second pcb
      if(i==1){ 
	// Only on top
	cablepar[0] = (supporthlength-kPcbLength-kVframeLength)/2.;
	lCableX = (kPcbLength+kVframeLength)/2.;
	lCableDY = (kYpos5[i+1]+kYpos5[i+2])/2.-cablepar[1]; // half way between 2 slats on same side
	lCableZ = TMath::Power(-1,i)*(nomexthickness/2.+carbonthickness+(-1+iCable++)*kCableWidth+kCableWidth/2.);
	gMC->Gsposp("S09L",10*i+cIndex++,"S09S",lCableX,lCableY+lCableDY,lCableZ,0,"ONLY",cablepar,3); 
	gMC->Gsposp("S09L",10*i+cIndex++,"S09S",lCableX,-(lCableY+lCableDY),lCableZ,0,"ONLY",cablepar,3); 
      }	

      // Cables going out from the end of the slats
      // First top cables
      cablepar[0] = (supporthlength-(slatLength5[i]+kXpos5[i]+kDslatLength)+kVframeLength)/2.;
      lCableX = slatLength5[i]+kXpos5[i]-kVframeLength+kDslatLength+cablepar[0]-supporthlength/2.;
      if(i+1>=kNslats5 || i+2>=kNslats5){ // If no more higher slats, then use distance to lower slat
	lCableDY = (kYpos5[i]+kYpos5[i-1])/2.-cablepar[1];
      }
      else {
	lCableDY = (kYpos5[i+1]+kYpos5[i+2])/2.-cablepar[1];
      }
      lCableZ = TMath::Power(-1,i)*(nomexthickness/2.+carbonthickness+(-1+iCable++)*kCableWidth+kCableWidth/2.);	
      gMC->Gsposp("S09L",10*i+cIndex++,"S09S",lCableX,lCableY+lCableDY,lCableZ,0,"ONLY",cablepar,3); 
      gMC->Gsposp("S09L",10*i+cIndex++,"S09S",lCableX,-(lCableY+lCableDY),lCableZ,0,"ONLY",cablepar,3); 
      if(i>0){
	if (i==1) { // Rounded slat. Bottom cable starts at dMotherInner (beam pipe)
	  lCableDY = (kYpos5[i]+kYpos5[i])/2.-cablepar[1];
	}
	else {
	  lCableDY = (kYpos5[i]+kYpos5[i-1])/2.-cablepar[1];
	  if ((lCableY-lCableDY)<(dMotherInner+cablepar[1])){
	      lCableDY = lCableY - dMotherInner - cablepar[1];
	  }
	}
	gMC->Gsposp("S09L",10*i+cIndex++,"S09S",lCableX,lCableY-lCableDY,lCableZ,0,"ONLY",cablepar,3); 
	gMC->Gsposp("S09L",10*i+cIndex++,"S09S",lCableX,-(lCableY-lCableDY),lCableZ,0,"ONLY",cablepar,3); 
      }
    }

    Float_t dzCh9  = dzCh;
    TGeoTranslation* trSupport1St5   = new TGeoTranslation("trSupport1St5", supporthlength/2., 0. , dzCh9);
    TGeoRotation*    roSupportSt5    = new TGeoRotation("roSupportSt5",90.,180.,-90.);
    TGeoCombiTrans*  coSupport2St5   = new TGeoCombiTrans(-supporthlength/2., 0., -dzCh9, roSupportSt5);
    GetEnvelopes(13)->AddEnvelope("S09S", 0, 1, *trSupport1St5);  
    GetEnvelopes(12)->AddEnvelope("S09S", 0, 2, *coSupport2St5);  
    GetEnvelopes(15)->AddEnvelope("S09S", 0, 3, *trSupport1St5);   
    GetEnvelopes(14)->AddEnvelope("S09S", 0, 4, *coSupport2St5);


    // End of pannel support geometry    

    // cout << "Geometry for Station 5...... done" << endl;

  }

  delete [] fStations;

}

//______________________________________________________________________________
void AliMUONSlatGeometryBuilder::SetVolumes()
{
/// Defines the volumes for the station345 chambers.

  if (gAlice->GetModule("DIPO")) {
    // if DIPO is preset, the whole station will be placed in DDIP volume
    SetMotherVolume(4, "DDIP");
    SetMotherVolume(5, "DDIP");
    SetMotherVolume(6, "DDIP");
    SetMotherVolume(7, "DDIP");
  } 	
  SetVolume(4, "SC05I", true);
  SetVolume(5, "SC05O", true);
  SetVolume(6, "SC06I", true);
  SetVolume(7, "SC06O", true);
     
  if (gAlice->GetModule("SHIL")) {
    SetMotherVolume(8, "YOUT2");
    SetMotherVolume(9, "YOUT2");
    SetMotherVolume(10, "YOUT2");
    SetMotherVolume(11, "YOUT2");
    SetMotherVolume(12, "YOUT2");
    SetMotherVolume(13, "YOUT2");
    SetMotherVolume(14, "YOUT2");
    SetMotherVolume(15, "YOUT2");
  }  

  SetVolume( 8, "SC07I", true);
  SetVolume( 9, "SC07O", true);
  SetVolume(10, "SC08I", true);
  SetVolume(11, "SC08O", true);
  SetVolume(12, "SC09I", true);
  SetVolume(13, "SC09O", true);
  SetVolume(14, "SC10I", true);
  SetVolume(15, "SC10O", true);
}


//______________________________________________________________________________
void AliMUONSlatGeometryBuilder::SetTransformations()
{
/// Defines the transformations for the station345 chambers.

// Stations 345 are not perpendicular to the beam axis
// See AliMUONConstants class
  TGeoRotation st345inclination("rot99");
  st345inclination.RotateX(AliMUONConstants::St345Inclination());
  
// The rotation of the half-chamber is done with respect the center of the chamber.
// the distance beween the roation axis and the chamber position is 
// AliMUONConstants::DzCh()+AliMUONConstants::DzSlat()
// Therefore the position of the half-chamber has to be corrected by a traslation in Z and Y axis
  Double_t deltaY = (AliMUONConstants::DzCh()+AliMUONConstants::DzSlat())*
    TMath::Sin(AliMUONConstants::St345Inclination() * TMath::Pi()/180.);
  Double_t deltaZ = (AliMUONConstants::DzCh()+AliMUONConstants::DzSlat())*
    (1.-TMath::Cos(AliMUONConstants::St345Inclination() * TMath::Pi()/180.));


  Double_t zpos1= - AliMUONConstants::DefaultChamberZ(4); 
  SetTransformation(4, TGeoTranslation(0., -deltaY, -deltaZ+zpos1), st345inclination);
  SetTransformation(5, TGeoTranslation(0.,  deltaY,  deltaZ+zpos1), st345inclination);

  zpos1= - AliMUONConstants::DefaultChamberZ(5); 
  SetTransformation(6, TGeoTranslation(0., -deltaY, -deltaZ+zpos1), st345inclination);
  SetTransformation(7, TGeoTranslation(0.,  deltaY,  deltaZ+zpos1), st345inclination);

  zpos1 = - AliMUONConstants::DefaultChamberZ(6); 
  SetTransformation(8, TGeoTranslation(0., -deltaY, -deltaZ+zpos1), st345inclination);
  SetTransformation(9, TGeoTranslation(0.,  deltaY,  deltaZ+zpos1), st345inclination);

  zpos1 = - AliMUONConstants::DefaultChamberZ(7); 
  SetTransformation(10, TGeoTranslation(0., -deltaY, -deltaZ+zpos1), st345inclination );
  SetTransformation(11, TGeoTranslation(0.,  deltaY,  deltaZ+zpos1), st345inclination );

  zpos1 = - AliMUONConstants::DefaultChamberZ(8); 
  SetTransformation(12, TGeoTranslation(0., -deltaY, -deltaZ+zpos1), st345inclination);
  SetTransformation(13, TGeoTranslation(0.,  deltaY,  deltaZ+zpos1), st345inclination);

  zpos1 = - AliMUONConstants::DefaultChamberZ(9); 
  SetTransformation(14, TGeoTranslation(0., -deltaY, -deltaZ+zpos1), st345inclination);
  SetTransformation(15, TGeoTranslation(0.,  deltaY,  deltaZ+zpos1), st345inclination);

}

//______________________________________________________________________________
void AliMUONSlatGeometryBuilder::SetSensitiveVolumes()
{
/// Defines the sensitive volumes for slat stations chambers.

  GetGeometry( 4)->SetSensitiveVolume("S05G");
  GetGeometry( 4)->SetSensitiveVolume("SC5G");
  GetGeometry( 4)->SetSensitiveVolume("SD5G");
  GetGeometry( 5)->SetSensitiveVolume("S05G");
  GetGeometry( 5)->SetSensitiveVolume("SC5G");
  GetGeometry( 5)->SetSensitiveVolume("SD5G");
  GetGeometry( 6)->SetSensitiveVolume("S06G");
  GetGeometry( 6)->SetSensitiveVolume("SC6G");
  GetGeometry( 6)->SetSensitiveVolume("SD6G");
  GetGeometry( 7)->SetSensitiveVolume("S06G");
  GetGeometry( 7)->SetSensitiveVolume("SC6G");
  GetGeometry( 7)->SetSensitiveVolume("SD6G");
  GetGeometry( 8)->SetSensitiveVolume("S07G");
  GetGeometry( 8)->SetSensitiveVolume("SD7G");
  GetGeometry( 9)->SetSensitiveVolume("S07G");
  GetGeometry( 9)->SetSensitiveVolume("SD7G");
  GetGeometry(10)->SetSensitiveVolume("S08G");
  GetGeometry(10)->SetSensitiveVolume("SD8G");
  GetGeometry(11)->SetSensitiveVolume("S08G");
  GetGeometry(11)->SetSensitiveVolume("SD8G");
  GetGeometry(12)->SetSensitiveVolume("S09G");
  GetGeometry(12)->SetSensitiveVolume("SD9G");
  GetGeometry(13)->SetSensitiveVolume("S09G");
  GetGeometry(13)->SetSensitiveVolume("SD9G");
  GetGeometry(14)->SetSensitiveVolume("S10G");
  GetGeometry(14)->SetSensitiveVolume("SD0G");
  GetGeometry(15)->SetSensitiveVolume("S10G");
  GetGeometry(15)->SetSensitiveVolume("SD0G");
}

//______________________________________________________________________________
Int_t  AliMUONSlatGeometryBuilder::ConvertSlatNum(Int_t numslat, Int_t quadnum, Int_t fspq) const
{
/// On-line function establishing the correspondance between numslat (the slat number on a particular quadrant (numslat->0....4 for St3))
/// and slatnum (the slat number on the whole panel (slatnum->1...18 for St3)
  numslat += 1;
  if (quadnum==2 || quadnum==3) 
    numslat += fspq;
  else
    numslat = fspq + 2-numslat;
  numslat -= 1;
	      
  if (quadnum==3 || quadnum==4) numslat += 2*fspq+1;

  return numslat;
}
