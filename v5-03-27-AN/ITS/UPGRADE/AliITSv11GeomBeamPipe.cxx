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

//*************************************************************************
// This class Defines the Geometry of the Beam Pipe
// for the ITS Upgrade using TGeo
//
//  Mario Sitta <sitta@to.infn.it>
//*************************************************************************


/* $Id: AliITSv11GeomBeamPipe.cxx  */
// General Root includes
#include <TMath.h>
// Root Geometry includes
//#include <AliLog.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoPcon.h>
#include <TGeoCone.h>
#include <TGeoTube.h> // contaings TGeoTubeSeg
#include <TGeoMatrix.h>
#include "AliITSv11GeomBeamPipe.h"

ClassImp(AliITSv11GeomBeamPipe)

//________________________________________________________________________
AliITSv11GeomBeamPipe::AliITSv11GeomBeamPipe(): 
  AliITSv11Geometry(),
  fBeamPipeRmin(0),
  fBeamPipeRmax(0),
  fBeamPipeZlen(0)
{
  //
  // Standard constructor
  //
}

//________________________________________________________________________
AliITSv11GeomBeamPipe::AliITSv11GeomBeamPipe(Double_t rmin, Double_t rmax, Double_t zlen, Int_t debug): 
  AliITSv11Geometry(debug),
  fBeamPipeRmin(rmin),
  fBeamPipeRmax(rmax),
  fBeamPipeZlen(zlen)
{
  //
  // Constructor setting parameters and debugging level
  //
}

//________________________________________________________________________
AliITSv11GeomBeamPipe::AliITSv11GeomBeamPipe(const AliITSv11GeomBeamPipe &s):
  AliITSv11Geometry(s.GetDebug()),
  fBeamPipeRmin(s.fBeamPipeRmin),
  fBeamPipeRmax(s.fBeamPipeRmax),
  fBeamPipeZlen(s.fBeamPipeZlen)
{
  //
  // Copy constructor
  //
}

//________________________________________________________________________
AliITSv11GeomBeamPipe& AliITSv11GeomBeamPipe::operator=(const AliITSv11GeomBeamPipe &s)
{
  //
  // Assignment operator 
  //
  if(&s == this) return *this;

  fBeamPipeRmin = s.fBeamPipeRmin;
  fBeamPipeRmax = s.fBeamPipeRmax;
  fBeamPipeZlen = s.fBeamPipeZlen;

  return *this;
}

//________________________________________________________________________
AliITSv11GeomBeamPipe::~AliITSv11GeomBeamPipe() {
  //
  // Destructor
  //
}

//________________________________________________________________________
void AliITSv11GeomBeamPipe::CreateBeamPipe(TGeoVolume *moth,
					   const TGeoManager *mgr){
//
// Creates the actual Beam Pipe and places it inside the mother volume
//
// Input:
//         moth : the TGeoVolume owing the volume structure
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      03 Mar 2012  Mario Sitta
//


  // Local variables
//none

  // Check if the user set the proper parameters
  if (fBeamPipeRmin <= 0)
    AliFatal(Form("Wrong beam pipe internal radius (%f)",fBeamPipeRmin));
  if (fBeamPipeRmax <= 0)
    AliFatal(Form("Wrong beam pipe external radius (%f)",fBeamPipeRmax));
  if (fBeamPipeZlen <= 0)
    AliFatal(Form("Wrong beam pipe half length (%f)",fBeamPipeZlen));

  if (fBeamPipeRmin >= fBeamPipeRmax)
    AliFatal(Form("Internal beam pipe radius (%f) greater than external radius (%f)",fBeamPipeRmin,fBeamPipeRmax));


  // The shape of the beam pipe: a simple cylinder
  TGeoTube *bpShape = new TGeoTube(fBeamPipeRmin,fBeamPipeRmax,fBeamPipeZlen);


  // We have all shapes: now create the real volumes
  TGeoMedium *medBe = mgr->GetMedium("ITS_BERILLIUM$");

  TGeoVolume *bpVolume = new TGeoVolume("UpgradeBeamPipe", bpShape, medBe);


  // Finally put everything in the mother volume
  moth->AddNode(bpVolume, 1, 0);


  // Upgrade beam pipe geometry is served
  return;
}

