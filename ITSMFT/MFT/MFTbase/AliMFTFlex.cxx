
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

//-------------------------------------------------------------------------
//  Flex class for ALICE MFT upgrade
//  This version uses TGeo
//  Authors:
//  F. Manso
//-------------------------------------------------------------------------


#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"
#include "TGeoBBox.h"

#include "AliLog.h"

#include "AliMFTConstants.h"
#include "AliMFTLadderSegmentation.h"
#include "AliMFTFlex.h"

ClassImp(AliMFTFlex)

//_____________________________________________________________________________
AliMFTFlex::AliMFTFlex():
TNamed(),
fTrackWidth(200e-4),      // 200 microns
fTrackThickness(25e-4),   // thickness of a sensor track
fTrackNb(2),              // 2 tracks per sensor
fGlueThickness(25e-4),    // 25 microns
fKaptonThickness(25e-4),  // 25 microns
fVarnishThickness(30e-4), // 30 microns
fNLayer(0),               // number of layers
fNSensors(0),
fChipWidth(0),
fChipInterspace(0),
fChipSideOffset(0),
fIsLeftType(kFALSE),
fLadderSeg(NULL)
{
  // Constructor
  fFlexDimensions = new Double_t[3];
  fFlexOrigin = new Double_t[3];
  for(int i =0;i<3;i++) {
    fFlexDimensions[i]=0.;
    fFlexOrigin[i]=0.;
  }
}
//=============================================================================================

AliMFTFlex::~AliMFTFlex() {
  delete fFlexDimensions;
  delete fFlexOrigin;
  
}
//=============================================================================================


AliMFTFlex::AliMFTFlex(AliMFTLadderSegmentation *ladder):
TNamed(),
fTrackWidth(2000e-4),      // 200 microns
fTrackThickness(250e-4),   // thickness of a sensor track
fTrackNb(10),             // track number per sensor
fGlueThickness(250e-4),    // 25 microns
fKaptonThickness(250e-4),  // 25 microns
fVarnishThickness(300e-4), // 30 microns
fNLayer(2),               // number of layers
fNSensors(ladder->GetNSensors()),
fChipWidth(AliMFTConstants::fChipWidth),
fChipInterspace(AliMFTConstants::fChipInterspace),
fChipSideOffset(AliMFTConstants::fChipSideOffset),
fLadderSeg(ladder)
{
  // Constructor
  fFlexDimensions = new Double_t[3];
  fFlexOrigin = new Double_t[3];
}

//=============================================================================================

TGeoVolumeAssembly* AliMFTFlex::MakeFlex(Int_t nsensors)
{
  
  AliInfo(Form("Start Creating the flex for  %d Sensors",nsensors));
  // Flex Container
  TGeoVolumeAssembly* Flex = new TGeoVolumeAssembly("Flex");
  
  // Building of each flex layer
  TGeoVolumeAssembly*  alutrack = MakeAluTrack(nsensors);
  TGeoVolume*         gluelayer = MakeGlueLayer(fNLayer, fGlueThickness, fFlexDimensions[1], fFlexDimensions[0]);
  TGeoVolume*       kaptonlayer = MakeKaptonLayer(fNLayer, fKaptonThickness, fFlexDimensions[1], fFlexDimensions[0]);
  TGeoVolume*      varnishlayer = MakeVarnishLayer(fVarnishThickness, fFlexDimensions[1], fFlexDimensions[0]);
  
  // Local frame: x=lenght, y=width, z=thickness, left flex
  Double_t zkapton, zglue, ztrack, zvarnish;
  // Don't change the order of the next 4 lines..
  zkapton  = fKaptonThickness*fNLayer/2;
  zglue    = 2*zkapton + fGlueThickness*fNLayer/2;
  ztrack   = zglue + fGlueThickness*fNLayer/2;
  zvarnish = ztrack + fTrackWidth*fTrackThickness*fTrackNb/fFlexDimensions[1] + fVarnishThickness/2;
  
  fFlexDimensions[2] = zvarnish + fVarnishThickness/2;
//  fLadderSeg->SetFlexLength(fFlexDimensions);
  
  // Final flex building
  Flex->AddNode(kaptonlayer,  1,  new TGeoTranslation((1.-2.*!fIsLeftType)*fFlexDimensions[0]/2, -fFlexDimensions[1]/2, zkapton));
  Flex->AddNode(gluelayer,    1,  new TGeoTranslation((1.-2.*!fIsLeftType)*fFlexDimensions[0]/2, -fFlexDimensions[1]/2, zglue));
  Flex->AddNode(alutrack,     1,  new TGeoTranslation( 0. , 0., ztrack));
  Flex->AddNode(varnishlayer, 1,  new TGeoTranslation((1.-2.*!fIsLeftType)*fFlexDimensions[0]/2, -fFlexDimensions[1]/2, zvarnish));
  return Flex;
}

//=============================================================================================

TGeoVolumeAssembly* AliMFTFlex::MakeAluTrack(Int_t nsensors)
{
  // For the al tracks --> nsensors+1 volumes (sensors + flex end) with an al thickness according to the number of tracks
  // Local frame: x=lenght, y=width, z=thickness, building of a LEFT flex
  char name[10];
  Double_t xpos,ypos,zpos,section,lenghtfe,widthalu;
  //Int_t ntracktotal;
  //ntracktotal = nsensors*fTrackNb;

  section=fTrackWidth*fTrackThickness;  // straight-edged of a track
  widthalu = section*fTrackNb/fFlexDimensions[1]; // width of alu for the first sensor
 
  TGeoMedium* kMedAlu[5];
  for (Int_t i=1; i <= nsensors; i++) {
    sprintf(name, "MFT_Alu%d", i);
    kMedAlu[i-1] = gGeoManager->GetMedium(name);
  }
  
  // Overall container for the flex tracks
  TGeoVolumeAssembly* alutrack = new TGeoVolumeAssembly("AluTrack");

  // tracks from each sensor
  for (Int_t i=1; i <= nsensors; i++) {
    AliInfo(Form("Sensor %d",i));
    xpos = fChipSideOffset + (fChipWidth+fChipInterspace)/2 + (i-1)*(fChipWidth + fChipInterspace);
    ypos = - fFlexDimensions[1]/2;
    zpos = widthalu/2;
    AliInfo(Form("X,Y,Z =  %f %f %f",xpos,ypos,zpos));
    AliInfo(Form("dX,dY,dZ =  %f %f %f",(fChipWidth+fChipInterspace)/2, fFlexDimensions[1]/2, widthalu/2));
    
    sprintf(name, "AluTrack%d", i);
    alutrack->AddNode(new TGeoVolume(name, new TGeoBBox("BOX", (fChipWidth+fChipInterspace)/2, fFlexDimensions[1]/2, widthalu/2), kMedAlu[i-1]),
                      i, new TGeoTranslation((1.-2.*!fIsLeftType)*xpos, ypos, zpos));
  }
  // tracks at the flex end = total tracks from sensors
  lenghtfe = fFlexDimensions[0] - fChipSideOffset - nsensors*(fChipWidth + fChipInterspace);
  alutrack->AddNode(new TGeoVolume("AluTrackFE", new TGeoBBox("BOX", lenghtfe/2, fFlexDimensions[1]/2, widthalu/2), kMedAlu[nsensors-1]),
                    nsensors+1, new TGeoTranslation((1.-2.*!fIsLeftType)*(xpos + (fChipWidth+fChipInterspace)/2 + lenghtfe/2), ypos, zpos));
  alutrack->SetVisibility(1);
  alutrack->SetLineColor(kBlack);
  return alutrack;
}
//=============================================================================================

TGeoVolume* AliMFTFlex::MakeGlueLayer(Int_t nlayer, Double_t thickness, Double_t widthflex, Double_t lenghtflex)
{
  TGeoMedium *kMedEpoxy = gGeoManager->GetMedium("MFT_Epoxy");
  TGeoVolume* gluelayer = new TGeoVolume("gluelayer", new TGeoBBox("BOX", lenghtflex/2, widthflex/2, thickness*nlayer/2), kMedEpoxy);
  gluelayer->SetVisibility(1);
  gluelayer->SetLineColor(kBlue+1);
  return gluelayer;
}
//=============================================================================================

TGeoVolume* AliMFTFlex::MakeKaptonLayer(Int_t nlayer, Double_t thickness, Double_t widthflex, Double_t lenghtflex)
{
  TGeoMedium *kMedKapton = gGeoManager->GetMedium("MFT_Kapton");
  AliInfo(Form("dX,dY,dZ =  %f %f %f",lenghtflex/2, widthflex/2,  thickness*nlayer/2));
  TGeoVolume* kaptonlayer = new TGeoVolume("kaptonlayer", new TGeoBBox("BOX", lenghtflex/2, widthflex/2,  thickness*nlayer/2), kMedKapton);
  kaptonlayer->SetVisibility(1);
  kaptonlayer->SetLineColor(kOrange);
  return kaptonlayer;
}
//=============================================================================================

TGeoVolume* AliMFTFlex::MakeVarnishLayer(Double_t thickness, Double_t widthflex,  Double_t lenghtflex)
{
  // one varnish layer only
  TGeoMedium *kMedVarnish = gGeoManager->GetMedium("MFT_Epoxy");  // varnish = epoxy ???
  TGeoVolume* varnishlayer = new TGeoVolume("varnish", new TGeoBBox("BOX", lenghtflex/2, widthflex/2, thickness/2), kMedVarnish);
  varnishlayer->SetVisibility(1);
  varnishlayer->SetLineColor(kYellow);
  return varnishlayer;
}


