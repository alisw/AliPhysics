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


/////////////////////////////////////////////////////////////
//  Manager and hits classes for set: RICH default version //
/////////////////////////////////////////////////////////////

#include "Riostream.h"

#include <TNode.h> 
#include <TRandom.h> 
#include <TTUBE.h>
#include <TVirtualMC.h>

#include "AliConst.h" 
#include "AliPDG.h" 
#include "AliRICHGeometry.h"
#include "AliRICHResponse.h"
#include "AliRICHResponseV0.h"
#include "AliRICHSegmentationV0.h"
#include "AliRICHv0.h"
#include "AliRun.h"

ClassImp(AliRICHv0)
    
//___________________________________________
AliRICHv0::AliRICHv0(const char *name, const char *title)
          :AliRICH(name,title)
{
//
// Version 0
// Default Segmentation, no hits
    AliRICHSegmentationV0* segmentation = new AliRICHSegmentationV0;
//
//  Segmentation parameters
    segmentation->SetPadSize(0.84,0.80);
    segmentation->SetDAnod(0.84/2);
//
//  Geometry parameters
    AliRICHGeometry* geometry = new AliRICHGeometry;
    geometry->SetGapThickness(8);
    geometry->SetProximityGapThickness(.4);
    geometry->SetQuartzLength(131);
    geometry->SetQuartzWidth(126.2);
    geometry->SetQuartzThickness(.5);
    geometry->SetOuterFreonLength(131);
    geometry->SetOuterFreonWidth(40.3);
    geometry->SetInnerFreonLength(131);
    geometry->SetInnerFreonWidth(40.3);
    geometry->SetFreonThickness(1.5);
//
//  Response parameters
    AliRICHResponseV0*  response   = new AliRICHResponseV0;
    response->SetSigmaIntegration(5.);
    response->SetChargeSlope(27.);
    response->SetChargeSpread(0.18, 0.18);
    response->SetMaxAdc(4096);
    response->SetAlphaFeedback(0.036);
    response->SetEIonisation(26.e-9);
    response->SetSqrtKx3(0.77459667);
    response->SetKx2(0.962);
    response->SetKx4(0.379);
    response->SetSqrtKy3(0.77459667);
    response->SetKy2(0.962);
    response->SetKy4(0.379);
    response->SetPitch(0.25);
    response->SetWireSag(0);                     // 1->On, 0->Off
    response->SetVoltage(2150);                  // Should only be 2000, 2050, 2100 or 2150
//
//    AliRICH *RICH = (AliRICH *) gAlice->GetDetector("RICH"); 
    
    fCkovNumber=0;
    fFreonProd=0;
    Int_t i=0;
    
    fChambers = new TObjArray(kNCH);
    for (i=0; i<kNCH; i++) {
      
      //PH      (*fChambers)[i] = new AliRICHChamber();  
      fChambers->AddAt(new AliRICHChamber(), i);  
      
    }
  
    for (i=0; i<kNCH; i++) {
      SetGeometryModel(i,geometry);
      SetSegmentationModel(i, segmentation);
      SetResponseModel(i, response);
    }
}//name ctor
//______________________________________________________________________________
void AliRICHv0::StepManager()
{//
}//AliRICHv0::StepManager()
//______________________________________________________________________________
