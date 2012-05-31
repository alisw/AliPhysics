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

#include <TVirtualMC.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoMedium.h>
#include <TGeoMatrix.h>
#include <TGeoVolume.h>
#include <TGeoTube.h>
#include <TGeoCone.h>
#include <TGeoPcon.h>
#include <TGeoCompositeShape.h>

#include "AliConst.h"
#include "AliEMCALSpaceFrame.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliLog.h"
 
ClassImp(AliEMCALSpaceFrame)
 
//_____________________________________________________________________________
AliEMCALSpaceFrame::AliEMCALSpaceFrame() 
  : TNamed("EMCALSpaceFrame","Steel Space Frame that supports EMCAL"),
    fNumCross(0),
    fNumSubSets(0),
    fTotalHalfWidth(0.),
    fBeginPhi(0.),
    fEndPhi(0.),
    fTotalPhi(0.),
    fBeginRadius(0.),
    fHalfFrameTrans(0.),
    fFlangeHeight(0.),
    fFlangeWidth(0.),
    fRibHeight(0.),
    fRibWidth(0.),
    fCrossBottomWidth(0.),
    fCrossTopWidth(0.),
    fCrossBottomHeight(0.),
    fCrossBottomRadThick(0.),
    fCrossBeamArcLength(0.),
    fCrossBottomStartRadius(0.),
    fCrossTopHeight(0.),
    fCrossTopRadThick(0.),
    fCrossTopStart(0.),
    fEndRadius(0.),
    fEndBeamRadThick(0),
    fEndBeamBeginRadius(0)
{
  // default constructor for EMCAL Space Frame
  //initialize parameters
  fNumCross = 12; 
  fNumSubSets = 3;
  fTotalHalfWidth = 152.3;   // Half Width of a Half Frame
			     // (CalFrame comes in 2 sections) 
  fBeginPhi = 76.8;
  fEndPhi = 193.03;
  fBeginRadius = 490.;
  
  fHalfFrameTrans = fTotalHalfWidth+57.2/2.;  // Half Frame Connector is 57.2cm wide,
			    // Supermodule is 340cm wide
                            // Sources: HALF-FRAME-CONNECTOR-27E226A.pdf
                            // provided by LBL

  fFlangeWidth = 15.2;
  fRibWidth = 1.5;
  fCrossBottomHeight = 15.2;
  fCrossBottomRadThick = 1.5;
  fCrossTopHeight = 1.5;
  fCrossTopRadThick = 35./2.;

  fTotalPhi = fEndPhi - fBeginPhi;
  fFlangeHeight = fBeginRadius + 3.;
  fRibHeight = fFlangeHeight + 35;
  fCrossBottomWidth = 0.5/(Double_t)fNumSubSets * (2.*fTotalHalfWidth - 8. * fFlangeWidth);
  fCrossTopWidth = fCrossBottomWidth; // fCrossBottomWidth + fFlangeWidth - fRibWidth;
                                      // for future release pending
                                      // overlap correction - new TGeoVolume creation

  fCrossBeamArcLength = (112.62597)/(fNumCross-1)-.001; // To account for shape of TGeoBBox
  fCrossBottomStartRadius = fBeginRadius + fCrossBottomRadThick;
  fCrossTopStart = fBeginRadius + 2.*fCrossBottomRadThick + fCrossTopRadThick+0.015; // 0.015 is a 
										     // bubblegum and duct tape
										     // fix for an overlap problem
										     // will be worked out in future releases 
  fEndRadius = fRibHeight+1.15;
  fEndBeamRadThick = fCrossBottomRadThick+fCrossTopRadThick;
  fEndBeamBeginRadius = fBeginRadius + fEndBeamRadThick;
}

//_____________________________________________________________________________
AliEMCALSpaceFrame::AliEMCALSpaceFrame(const AliEMCALSpaceFrame &frame) 
  : TNamed(frame.GetName(),frame.GetTitle()),
    fNumCross(frame.fNumCross),
    fNumSubSets(frame.fNumSubSets),
    fTotalHalfWidth(frame.fTotalHalfWidth),
    fBeginPhi(frame.fBeginPhi),
    fEndPhi(frame.fEndPhi),
    fTotalPhi(frame.fTotalPhi),
    fBeginRadius(frame.fBeginRadius),
    fHalfFrameTrans(frame.fHalfFrameTrans),
    fFlangeHeight(frame.fFlangeHeight),
    fFlangeWidth(frame.fFlangeWidth),
    fRibHeight(frame.fRibHeight),
    fRibWidth(frame.fRibWidth),
    fCrossBottomWidth(frame.fCrossBottomWidth),
    fCrossTopWidth(frame.fCrossTopWidth),
    fCrossBottomHeight(frame.fCrossBottomHeight),
    fCrossBottomRadThick(frame.fCrossBottomRadThick),
    fCrossBeamArcLength(frame.fCrossBeamArcLength),
    fCrossBottomStartRadius(frame.fCrossBottomStartRadius),
    fCrossTopHeight(frame.fCrossTopHeight),
    fCrossTopRadThick(frame.fCrossTopRadThick),
    fCrossTopStart(frame.fCrossTopStart),
    fEndRadius(frame.fEndRadius),
    fEndBeamRadThick(frame.fEndBeamRadThick),
    fEndBeamBeginRadius(frame.fEndBeamBeginRadius)
{
  // copy constructor for EMCAL Space Frame

}

//_____________________________________________________________________________
void AliEMCALSpaceFrame::CreateGeometry()
{
  AliDebug(1,"Create CalFrame Geometry");
  //////////////////////////////////////Setup/////////////////////////////////////////
  TGeoVolume* top = gGeoManager->GetVolume("ALIC");
  TGeoMedium *steel = gGeoManager->GetMedium("EMCAL_S steel$");
  TGeoMedium *air = gGeoManager->GetMedium("EMCAL_Air$");

	
  //////////////////////////////////// Volumes ///////////////////////////////////////  
  TGeoVolume *calFrameMO = 
    gGeoManager->MakeTubs("CalFrame", air, fBeginRadius-2.1,fEndRadius,
			  fTotalHalfWidth*3,fBeginPhi-3,fEndPhi+3);	// Mother Volume

  calFrameMO->SetVisibility(kFALSE);

  // Half Frame Mother Volume
  TGeoVolume *calHalfFrameMO = 
    gGeoManager->MakeTubs("HalfFrame", air, fBeginRadius-2,fEndRadius,
			  fTotalHalfWidth,fBeginPhi-2.9,fEndPhi+2.9);

  calHalfFrameMO->SetVisibility(kFALSE);

  TGeoVolume *endBeams = 
    gGeoManager->MakeBox("End Beams", steel, fEndBeamRadThick, fCrossTopHeight, fTotalHalfWidth); // End Beams

  TGeoVolume *skin = 
    gGeoManager->MakeTubs("skin", steel, fRibHeight+0.15, fEndRadius, 
  			  fTotalHalfWidth, fBeginPhi, fEndPhi);// back frame

  TGeoVolume *flangeVolume = 
    gGeoManager->MakeTubs("supportBottom", steel, fBeginRadius, fFlangeHeight, 
			  fFlangeWidth, fBeginPhi, fEndPhi);   		// FlangeVolume Beams

  TGeoVolume *ribVolume = 
    gGeoManager->MakeTubs("RibVolume", steel, fFlangeHeight, fRibHeight, fRibWidth, fBeginPhi, fEndPhi);

  TGeoVolume *subSetCross = 
    gGeoManager->MakeTubs("subSetCross", air, fBeginRadius-1,  fBeginRadius+2*fCrossBottomRadThick+
                          2*fCrossTopRadThick+0.15, fCrossBottomWidth, fBeginPhi, fEndPhi); 		// Cross Beam Containers
  subSetCross->SetVisibility(kFALSE);
  /*						// Obsolete for now
  TGeoVolume *subSetCrossTop = 
    gGeoManager->MakeTubs("SubSetCrossTop", air, fBeginRadius+2*fCrossBottomRadThick-1, fBeginRadius+2*fCrossBottomRadThick+
                          2*fCrossTopRadThick+1, fCrossTopWidth, fBeginPhi, fEndPhi);	// Cross 
  subSetCrossTop->SetVisibility(kFALSE);
  */                    
  TGeoVolume *crossBottomBeams = 
    gGeoManager->MakeBox("crossBottom", steel, fCrossBottomRadThick, fCrossBottomHeight, fCrossBottomWidth); // Cross Beams

  TGeoVolume *crossTopBeams = 
    gGeoManager->MakeBox("crossTop", steel, fCrossTopRadThick, fCrossTopHeight, fCrossTopWidth); // Cross Beams
  
  TGeoTranslation *trTEST = new TGeoTranslation();
  TGeoRotation *rotTEST = new TGeoRotation();
  
  Double_t conv = TMath::Pi()/180.;
  Double_t radAngle = 0;
  Double_t endBeamParam=.4;
  //cout<<"\nfCrossBottomStartRadius: "<<fCrossBottomStartRadius<<"\n";
  
  for(Int_t i = 0; i < fNumCross; i++){
    
    Double_t loopPhi = fBeginPhi + 1.8;

    // Cross Bottom Beams
    
    radAngle = (loopPhi + i*fCrossBeamArcLength)*conv; 
    
    rotTEST->SetAngles(fBeginPhi + i*fCrossBeamArcLength, 0, 0); //  SetTranslation(Double_t dx, Double_t dy, Double_t dz);
    trTEST->SetTranslation(cos(radAngle)*fCrossBottomStartRadius, sin(radAngle)*fCrossBottomStartRadius,0);

    TGeoCombiTrans *combo = new TGeoCombiTrans(*trTEST, *rotTEST); 	// TGeoTranslation &tr, const TGeoRotation &rot);
    combo->RegisterYourself();
    crossBottomBeams->SetVisibility(1);
    subSetCross->AddNode(crossBottomBeams, i+1, combo);
    if (i != 0 && i!=fNumCross-1){
    // Cross Bottom Beams
    rotTEST->SetAngles(fBeginPhi + i*fCrossBeamArcLength, 0, 0); //  SetTranslation(Double_t dx, Double_t dy, Double_t dz);
    trTEST->SetTranslation(cos(radAngle)*fCrossTopStart, sin(radAngle)*fCrossTopStart,0);
    crossTopBeams->SetVisibility(1);
    subSetCross->AddNode(crossTopBeams, i+1,  new TGeoCombiTrans(*trTEST, *rotTEST));
    }
    

    else if(i ==0){
		rotTEST->SetAngles(fBeginPhi + i*fCrossBeamArcLength, 0, 0); //  SetTranslation(Double_t dx, Double_t dy, Double_t dz);
    trTEST->SetTranslation(cos((77-endBeamParam)*conv)*(fEndBeamBeginRadius), sin((77-endBeamParam)*conv)*(fEndBeamBeginRadius),0);
    endBeams->SetVisibility(1);
    calHalfFrameMO->AddNode(endBeams, 1,  new TGeoCombiTrans(*trTEST, *rotTEST));
    }
    else{
    rotTEST->SetAngles(193.03, 0, 0); //  SetTranslation(Double_t dx, Double_t dy, Double_t dz);
    trTEST->SetTranslation(cos((193.03+endBeamParam)*conv)*(fEndBeamBeginRadius)/*more duct tape*/, sin((193.03+endBeamParam)*conv)*(fEndBeamBeginRadius),0);
    endBeams->SetVisibility(1);
    calHalfFrameMO->AddNode(endBeams, 2,  new TGeoCombiTrans(*trTEST, *rotTEST));
    }  
  }
  
  //Beam Containers

  // Translations 

  TGeoTranslation *origin1 = new TGeoTranslation(0,0,0); // Equivalent to gGeoIdentity
  TGeoTranslation *origin2 = new TGeoTranslation(0,0,2*(fCrossBottomWidth+fFlangeWidth));
  TGeoTranslation *origin3 = new TGeoTranslation(0,0,-2*(fCrossBottomWidth+fFlangeWidth));

  // FlangeVolume translations  
  TGeoTranslation *str1 = new TGeoTranslation(0,0,-3*(fCrossBottomWidth+fFlangeWidth));
  TGeoTranslation *str2 = new TGeoTranslation(0,0,-(fCrossBottomWidth+fFlangeWidth));
  TGeoTranslation *str3 = new TGeoTranslation(0,0,(fCrossBottomWidth+fFlangeWidth));
  TGeoTranslation *str4 = new TGeoTranslation(0,0,3*(fCrossBottomWidth+fFlangeWidth));

  // Half Frame Translations
  TGeoTranslation *halfTrans1 =  new TGeoTranslation(0,0,fHalfFrameTrans);
  TGeoTranslation *halfTrans2 =  new TGeoTranslation(0,0,-fHalfFrameTrans);

  // Beams Volume 
  calHalfFrameMO->AddNode(flangeVolume, 1, str1);
  calHalfFrameMO->AddNode(flangeVolume, 2, str2);
  calHalfFrameMO->AddNode(flangeVolume, 3, str3);
  calHalfFrameMO->AddNode(flangeVolume, 4, str4);
  
  calHalfFrameMO->AddNode(ribVolume, 1, str1);
  calHalfFrameMO->AddNode(ribVolume, 2, str2);
  calHalfFrameMO->AddNode(ribVolume, 3, str3);
  calHalfFrameMO->AddNode(ribVolume, 4, str4);
  
  // Cross Beams  
  calHalfFrameMO->AddNode(subSetCross, 1, origin1);
  calHalfFrameMO->AddNode(subSetCross, 2, origin2);
  calHalfFrameMO->AddNode(subSetCross, 3, origin3);
/*					// Obsolete for now
  calHalfFrameMO->AddNode(subSetCrossTop, 1, origin1);
  calHalfFrameMO->AddNode(subSetCrossTop, 2, origin2);
  calHalfFrameMO->AddNode(subSetCrossTop, 3, origin3);
*/

  calHalfFrameMO->AddNode(skin, 1, gGeoIdentity);
  
  calFrameMO->AddNode(calHalfFrameMO, 1, halfTrans1);
  calFrameMO->AddNode(calHalfFrameMO, 2, halfTrans2);
 
  top->AddNode(calFrameMO,1,gGeoIdentity);
//  cout<<"**********************************\nfEndRadius:\t"<<fEndRadius;
}

