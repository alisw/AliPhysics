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

/*
  $Log$
  Revision 1.7  2000/10/03 21:44:09  morsch
  Use AliSegmentation and AliHit abstract base classes.

  Revision 1.6  2000/07/10 15:28:39  fca
  Correction of the inheritance scheme

  Revision 1.5  2000/06/30 16:38:51  dibari
  Removed setters.

  Revision 1.4  2000/06/13 13:13:40  jbarbosa
  Correcting previous correction...

  Revision 1.3  2000/06/13 13:06:38  jbarbosa
  Fixed compiling error for HP (multiple declaration)

  Revision 1.2  2000/06/12 15:36:16  jbarbosa
  Cleaned up version.

  Revision 1.1  2000/06/09 15:00:31  jbarbosa
  New full version. All parameters configurable.

  Revision 1.9  2000/05/31 08:19:38  jbarbosa
  Fixed bug in StepManager

  Revision 1.8  2000/05/26 17:30:08  jbarbosa
  Cerenkov angle now stored within cerenkov data structure.

  Revision 1.7  2000/05/18 10:31:36  jbarbosa
  Fixed positioning of spacers inside freon.
  Fixed positioning of proximity gap
  inside methane.
  Fixed cut on neutral particles in the StepManager.

  Revision 1.6  2000/04/28 11:51:58  morsch
   Dimensions of arrays hits and Ckov_data corrected.

  Revision 1.5  2000/04/19 13:28:46  morsch
  Major changes in geometry (parametrised), materials (updated) and
  step manager (diagnostics) (JB, AM)

*/



//////////////////////////////////////////////////////////
//  Manager and hits classes for set: RICH full version //
//////////////////////////////////////////////////////////

#include <TTUBE.h>
#include <TNode.h> 
#include <TRandom.h> 
#include <TParticle.h> 

#include "AliRICHv1.h"
#include "AliRICHHit.h"
#include "AliSegmentation.h"
#include "AliRICHResponse.h"
#include "AliRICHSegmentationV0.h"
#include "AliRICHResponseV0.h"
#include "AliRICHGeometry.h"
#include "AliRun.h"
#include "AliMC.h"
#include "iostream.h"
#include "AliCallf77.h"
#include "AliConst.h" 
#include "AliPDG.h" 
#include "TGeant3.h"

ClassImp(AliRICHv1)
    
//___________________________________________
AliRICHv1::AliRICHv1()
{

// Default constructor fo AliRICHvv1 (full version)

    //fChambers = 0;
}

//___________________________________________
AliRICHv1::AliRICHv1(const char *name, const char *title)
    : AliRICH(name,title)
{

// Full version of RICH with hits and diagnostics

  // Version 0
// Default Segmentation, no hits
    AliRICHSegmentationV0* segmentationV0 = new AliRICHSegmentationV0;
//
//  Segmentation parameters
    segmentationV0->SetPadSize(0.84,0.80);
    segmentationV0->SetDAnod(0.84/2);
//
//  Geometry parameters
    AliRICHGeometry* geometry = new AliRICHGeometry;
    geometry->SetGapThickness(8);
    geometry->SetProximityGapThickness(.4);
    geometry->SetQuartzLength(133);
    geometry->SetQuartzWidth(127.9);
    geometry->SetQuartzThickness(.5);
    geometry->SetOuterFreonLength(133);
    geometry->SetOuterFreonWidth(41.3);
    geometry->SetInnerFreonLength(133);
    geometry->SetInnerFreonWidth(41.3);
    geometry->SetFreonThickness(1.5);
//
//  Response parameters
    AliRICHResponseV0*  responseV0   = new AliRICHResponseV0;
    responseV0->SetSigmaIntegration(5.);
    responseV0->SetChargeSlope(27.);
    responseV0->SetChargeSpread(0.18, 0.18);
    responseV0->SetMaxAdc(4096);
    responseV0->SetAlphaFeedback(0.036);
    responseV0->SetEIonisation(26.e-9);
    responseV0->SetSqrtKx3(0.77459667);
    responseV0->SetKx2(0.962);
    responseV0->SetKx4(0.379);
    responseV0->SetSqrtKy3(0.77459667);
    responseV0->SetKy2(0.962);
    responseV0->SetKy4(0.379);
    responseV0->SetPitch(0.25);
//
//
//    AliRICH *RICH = (AliRICH *) gAlice->GetDetector("RICH"); 
    
    fCkovNumber=0;
    fFreonProd=0;
    Int_t i=0;
    
    fChambers = new TObjArray(kNCH);
    for (i=0; i<kNCH; i++) {
      
      (*fChambers)[i] = new AliRICHChamber();  
      
    }
  
    for (i=0; i<kNCH; i++) {
      SetGeometryModel(i,geometry);
      SetSegmentationModel(i, segmentationV0);
      SetResponseModel(i, responseV0);
      SetNsec(i,1);
      SetDebugLevel(0);
    }


}

void AliRICHv1::Init()
{

  printf("*********************************** RICH_INIT ***********************************\n");
  printf("*                                                                               *\n");
  printf("*                        AliRICHv1 Full version started                         *\n");
  printf("*                                                                               *\n");

  
  AliSegmentation*  segmentation;
  AliRICHGeometry*  geometry;
  AliRICHResponse*  response;


    // 
    // Initialize Tracking Chambers
    //
    for (Int_t i=1; i<kNCH; i++) {
	//printf ("i:%d",i);
	( (AliRICHChamber*) (*fChambers)[i])->Init(i);  
    }  
    
    //
    // Set the chamber (sensitive region) GEANT identifier
    
    ((AliRICHChamber*)(*fChambers)[0])->SetGid(1);  
    ((AliRICHChamber*)(*fChambers)[1])->SetGid(2);  
    ((AliRICHChamber*)(*fChambers)[2])->SetGid(3);  
    ((AliRICHChamber*)(*fChambers)[3])->SetGid(4);  
    ((AliRICHChamber*)(*fChambers)[4])->SetGid(5);  
    ((AliRICHChamber*)(*fChambers)[5])->SetGid(6);  
    ((AliRICHChamber*)(*fChambers)[6])->SetGid(7); 

    Float_t pos1[3]={0,471.8999,165.2599};
    Chamber(0).SetChamberTransform(pos1[0],pos1[1],pos1[2],new TRotMatrix("rot993","rot993",90,0,70.69,90,19.30999,-90));

    Float_t pos2[3]={171,470,0};
    Chamber(1).SetChamberTransform(pos2[0],pos2[1],pos2[2],new TRotMatrix("rot994","rot994",90,-20,90,70,0,0));

    Float_t pos3[3]={0,500,0};
    Chamber(2).SetChamberTransform(pos3[0],pos3[1],pos3[2],new TRotMatrix("rot995","rot995",90,0,90,90,0,0));
    
    Float_t pos4[3]={-171,470,0};
    Chamber(3).SetChamberTransform(pos4[0],pos4[1],pos4[2], new TRotMatrix("rot996","rot996",90,20,90,110,0,0));  

    Float_t pos5[3]={161.3999,443.3999,-165.3};
    Chamber(4).SetChamberTransform(pos5[0],pos5[1],pos5[2],new TRotMatrix("rot997","rot997",90,340,108.1999,70,18.2,70));

    Float_t pos6[3]={0., 471.9, -165.3,};
    Chamber(5).SetChamberTransform(pos6[0],pos6[1],pos6[2],new TRotMatrix("rot998","rot998",90,0,109.3099,90,19.30999,90));

    Float_t pos7[3]={-161.399,443.3999,-165.3};
    Chamber(6).SetChamberTransform(pos7[0],pos7[1],pos7[2],new TRotMatrix("rot999","rot999",90,20,108.1999,110,18.2,110));
    
    segmentation=Chamber(0).GetSegmentationModel(0);
    geometry=Chamber(0).GetGeometryModel();
    response=Chamber(0).GetResponseModel();
    
     
    printf("*                            Pads            : %3dx%3d                          *\n",segmentation->Npx(),segmentation->Npy());
    printf("*                            Pad size        : %5.2f x%5.2f mm2                 *\n",segmentation->Dpx(),segmentation->Dpy()); 
    printf("*                            Gap Thickness   : %5.1f cm                         *\n",geometry->GetGapThickness());
    printf("*                            Radiator Width  : %5.1f cm                         *\n",geometry->GetQuartzWidth());
    printf("*                            Radiator Length : %5.1f cm                         *\n",geometry->GetQuartzLength());
    printf("*                            Freon Thickness : %5.1f cm                         *\n",geometry->GetFreonThickness());
    printf("*                            Charge Slope    : %5.1f ADC                        *\n",response->ChargeSlope());
    printf("*                            Feedback Prob.  : %5.2f %%                         *\n",response->AlphaFeedback()*100);
    printf("*                            Debug Level     : %3d                              *\n",GetDebugLevel());
    printf("*                                                                               *\n");
    printf("*                                   Success!                                    *\n");
    printf("*                                                                               *\n");
    printf("*********************************************************************************\n");

}








