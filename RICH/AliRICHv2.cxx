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


//////////////////////////////////////////////////////////
//  Manager and hits classes for set: RICH full version //
//////////////////////////////////////////////////////////

#include "Riostream.h"

#include <TNode.h> 
#include <TParticle.h> 
#include <TRandom.h> 
#include <TTUBE.h>
#include <TVirtualMC.h>

#include "AliConst.h" 
#include "AliPDG.h" 
#include "AliRICHHit.h"
#include "AliRICHv2.h"
#include "AliRun.h"

ClassImp(AliRICHv2)
    
//___________________________________________
AliRICHv2::AliRICHv2()
{

// Default constructor fo AliRICHvv2 (full version)

    //fChambers = 0;
}

//___________________________________________
AliRICHv2::AliRICHv2(const char *name, const char *title)
    : AliRICH(name,title)
{

// Full version of RICH with hits and diagnostics, CONFIURABLE

    fCkovNumber=0;
    fFreonProd=0;
  
    fChambers = new TObjArray(kNCH);
    for (Int_t i=0; i<kNCH; i++) {
    
	(*fChambers)[i] = new AliRICHChamber();  
	
    }
}

void AliRICHv2::Init()
{

  if(fDebug) {
    printf("%s: *********************************** RICH_INIT ***********************************\n",ClassName());
    printf("%s: *                                                                               *\n",ClassName());
    printf("%s: *                    AliRICHv2 Configurable version started                     *\n",ClassName());
    printf("%s: *                                                                               *\n",ClassName());
  }

  
  AliSegmentation*  segmentation;
  AliRICHGeometry*  geometry;
  AliRICHResponse*  response;


    // 
    // Initialize Tracking Chambers
    //
    for (Int_t i=0; i<kNCH; i++) {
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

    segmentation=Chamber(0).GetSegmentationModel(0);
    geometry=Chamber(0).GetGeometryModel();
    response=Chamber(0).GetResponseModel();

    Float_t offset       = 490 + 1.276 - geometry->GetGapThickness()/2;        //distance from center of mother volume to methane
    Float_t deltaphi     = 19.5;                                               //phi angle between center of chambers - z direction
    Float_t deltatheta   = 20;                                                 //theta angle between center of chambers - x direction
    Float_t cosphi       = TMath::Cos(deltaphi*TMath::Pi()/180);
    Float_t sinphi       = TMath::Sin(deltaphi*TMath::Pi()/180);
    Float_t costheta     = TMath::Cos(deltatheta*TMath::Pi()/180);
    Float_t sintheta     = TMath::Sin(deltatheta*TMath::Pi()/180);

    Float_t pos1[3]={0.                , offset*cosphi         , offset*sinphi};
    Float_t pos2[3]={offset*sintheta   , offset*costheta       , 0. };
    Float_t pos3[3]={0.                , offset                , 0.};
    Float_t pos4[3]={-offset*sintheta  , offset*costheta       , 0.};
    Float_t pos5[3]={offset*sinphi     , offset*costheta*cosphi, -offset*sinphi};
    Float_t pos6[3]={0.                , offset*cosphi         , -offset*sinphi};
    Float_t pos7[3]={ -offset*sinphi   , offset*costheta*cosphi, -offset*sinphi};

    Chamber(0).SetChamberTransform(pos1[0],pos1[1],pos1[2],new TRotMatrix("rot993","rot993",90., 0.               , 90. - deltaphi, 90.             , deltaphi, -90.           ));
    Chamber(1).SetChamberTransform(pos2[0],pos2[1],pos2[2],new TRotMatrix("rot994","rot994",90., -deltatheta      , 90.           , 90.- deltatheta , 0.      , 0.             ));
    Chamber(2).SetChamberTransform(pos3[0],pos3[1],pos3[2],new TRotMatrix("rot995","rot995",90., 0.               , 90.           , 90.             , 0.      , 0.             ));
    Chamber(3).SetChamberTransform(pos4[0],pos4[1],pos4[2],new TRotMatrix("rot996","rot996",90.,  deltatheta      , 90.           , 90 + deltatheta , 0.      , 0.             ));
    Chamber(4).SetChamberTransform(pos5[0],pos5[1],pos5[2],new TRotMatrix("rot997","rot997",90., 360. - deltatheta, 108.2         , 90.- deltatheta ,18.2     , 90 - deltatheta));
    Chamber(5).SetChamberTransform(pos6[0],pos6[1],pos6[2],new TRotMatrix("rot998","rot998",90., 0.               , 90 + deltaphi , 90.             , deltaphi, 90.            ));
    Chamber(6).SetChamberTransform(pos7[0],pos7[1],pos7[2],new TRotMatrix("rot999","rot999",90., deltatheta       , 108.2         , 90.+ deltatheta ,18.2     , 90 + deltatheta));
    
    if(fDebug) {    
      printf("%s: *                            Pads            : %3dx%3d                          *\n",
	     ClassName(),segmentation->Npx(),segmentation->Npy());
      printf("%s: *                            Pad size        : %5.2f x%5.2f mm2                 *\n",
	     ClassName(),segmentation->Dpx(),segmentation->Dpy()); 
      printf("%s: *                            Gap Thickness   : %5.1f cm                         *\n",
	     ClassName(),geometry->GetGapThickness());
      printf("%s: *                            Radiator Width  : %5.1f cm                         *\n",
	     ClassName(),geometry->GetQuartzWidth());
      printf("%s: *                            Radiator Length : %5.1f cm                         *\n",
	     ClassName(),geometry->GetQuartzLength());
      printf("%s: *                            Freon Thickness : %5.1f cm                         *\n",
	     ClassName(),geometry->GetFreonThickness());
      printf("%s: *                            Charge Slope    : %5.1f ADC                        *\n",
	     ClassName(),response->ChargeSlope());
      printf("%s: *                            Feedback Prob.  : %5.2f %%                          *\n",
	     ClassName(),response->AlphaFeedback()*100);
      printf("%s: *                            Debug Level     : %3d                              *\n",
	     ClassName(),GetDebugLevel());
      printf("%s: *                                                                               *\n",
	     ClassName());
      printf("%s: *********************************************************************************\n",
	     ClassName());
    }
}

