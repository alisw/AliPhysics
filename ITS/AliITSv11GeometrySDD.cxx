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
// SDD geometry, based on ROOT geometrical modeler
//
// Ludovic Gaudichet (Ludovic.Gaudichet@to.infn.it)
//*************************************************************************


#include <stdio.h>
#include <stdlib.h>

// General Root includes
#include <Riostream.h>
#include <TMath.h>

// Root Geometry includes
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoPcon.h>
#include <TGeoCone.h>
#include <TGeoTube.h>
#include <TGeoArb8.h>
#include <TGeoCompositeShape.h>
#include <TGeoMatrix.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>

#include "AliITSv11Geometry.h"
#include "AliITSv11GeometrySDD.h"



ClassImp(AliITSv11GeometrySDD)


AliITSv11GeometrySDD::AliITSv11GeometrySDD():AliITSv11Geometry() {
  fAddOnlyLadder3 = -1;
  fAddOnlyLadder4 = -1;
  fAddOnlySegment = 0;
  SetGeomParameters();
};
//----------------------------------------------------------------------
AliITSv11GeometrySDD::AliITSv11GeometrySDD(Int_t debug):AliITSv11Geometry(debug) {
  fAddOnlyLadder3 = -1;
  fAddOnlyLadder4 = -1;
  fAddOnlySegment = 0;
  SetGeomParameters();
};
//----------------------------------------------------------------------
void AliITSv11GeometrySDD::SetGeomParameters() {

  fSegmentLength       = 37.2*2*fgkmm;
  fLadderWidth         = 50.0*fgkmm;
  fLadderHeight        = 30.0*fgkmm;
  fLadderBeamRadius    =  0.6*fgkmm;
  fLadderLa            =  3.*fgkmm;
  fLadderHa            =  0.6*fgkmm; //total pifometer
  fLadderLb            =  3.7*fgkmm;
  fLadderHb            =  0.6*fgkmm; //total pifometer
  fLadderl             =  0.25*fgkmm;

  fBottomBeamAngle     = 56.5;
  fBeamSidePhi         = 65;

  fWaferThickness      = 0.3*fgkmm;
  fWaferWidth          = 72.5*fgkmm;
  fWaferLength         = 87.6*fgkmm;

  fHybridLength        = (37.2*2-15)*fgkmm; //total pifometer
  fHybridWidth         = 35*fgkmm; //total pifometer

  fLadWaferSep         = 2*fgkmm;
  fPinSuppWidth        = 2.5*fgkmm; // pifometer
  fPinSuppHeight       = 2.*fgkmm; // pifometer
  fPinSuppRmax         = 2.5/2.*fgkmm;
  fPinR                = 1.5/2.*fgkmm;
  fPinSuppLength       = 5.*fgkmm;
  fPinSuppThickness    = 0.5*fgkmm;
  fPinSuppConeAngle    = 4;

  fLay3Rmin            = 130.*fgkmm; //not min! Rmin virtual tube
  fLay3Rmax            = 190.*fgkmm; //not min! Rmax virtual tube
  fLay3Length          = (524.+0.)*fgkmm; // ladder+supporting rings (length of the virtual tube)
  fLay3LadderLength    = 524.*fgkmm;
  fLay3DetShortRadius  = 146.0*fgkmm;
  fLay3DetLongRadius   = 152.0*fgkmm;
  fLay3LaddShortRadius = fLay3DetShortRadius-fLadderBeamRadius+(8-1)*fgkmm; // radius from the center to the CF ladder
  fLay3LaddLongRadius  = fLay3DetLongRadius-fLadderBeamRadius+(8-1)*fgkmm; // radius from the center to the CF ladder
  fLay3LaddTopCornerEnd = 15.6*fgkmm;
  fLay3Ndet            = 6;
  fLay3Nladd           = 14;
  
  fLay4Rmin            = 220.*fgkmm; //not min! Rmin virtual tube
  fLay4Rmax            = 290.*fgkmm; //not min! Rmax virtual tube
  fLay4Length          = (671.+0.)*fgkmm; // ladder+supporting rings (length of the virtual tube)
  fLay4LadderLength    = 671.*fgkmm;
  fLay4DetShortRadius  = 235.0*fgkmm;
  fLay4DetLongRadius   = 240.5*fgkmm;
  fLay4LaddShortRadius = fLay4DetShortRadius-fLadderBeamRadius+(8-1)*fgkmm; // radius from the center to the CF ladder
  fLay4LaddLongRadius  = fLay4DetLongRadius-fLadderBeamRadius+(8-1)*fgkmm;  // radius from the center to the CF ladder
  fLay4LaddTopCornerEnd = 15.6*fgkmm;
  fLay4Ndet            = 8;
  fLay4Nladd           = 22;

  for (Int_t i=0; i<fLay3Ndet; i++)
    fLay3sensorZPos[i] = -fSegmentLength*fLay3Ndet/2.+fSegmentLength/2.+i*fSegmentLength;

  for (Int_t i=0; i<fLay4Ndet; i++)
    fLay4sensorZPos[i] = -fSegmentLength*fLay4Ndet/2.+fSegmentLength/2.+i*fSegmentLength;

};
//________________________________________________________________________
TGeoCombiTrans *AliITSv11GeometrySDD::
CreateCombiTrans(const char *name, Double_t dy, Double_t dz, Double_t dphi) {
    //
    // return the TGeoCombiTrans which make a translation in y and z
    // and a rotation in phi in the global coord system
    //

    TGeoTranslation t1(dy*CosD(90.+dphi),dy*SinD(90.+dphi), dz);
    TGeoRotation r1("",0.,0.,dphi);

    TGeoCombiTrans *combiTrans1 = new TGeoCombiTrans(name);
    combiTrans1->SetTranslation(t1);
    combiTrans1->SetRotation(r1);
    return combiTrans1;
};
//________________________________________________________________________
void AliITSv11GeometrySDD::AddTranslationTotCombiTrans(TGeoCombiTrans* ct,
                                                       Double_t dx,
                                                       Double_t dy,
                                                       Double_t dz) {
    // Add a dx,dy,dz translation to the initial TGeoCombiTrans
    const Double_t *vect = ct->GetTranslation();
    Double_t newVect[3];
    newVect[0] = vect[0]+dx; 
    newVect[1] = vect[1]+dy; 
    newVect[2] = vect[2]+dz;
    ct->SetTranslation(newVect);
};
//________________________________________________________________________
void AliITSv11GeometrySDD::Layer3(TGeoVolume *Moth) {
    // Insert the layer 3 in the mother volume. This is a virtual volume
    // containing ladders of layer 3 and the supporting rings

    TGeoTube *virtualLayer3Shape = new TGeoTube("ITSsddLayer3Shape",
                                       fLay3Rmin,fLay3Rmax,fLay3Length*0.5);
    TGeoMedium *airSDD = gGeoManager->GetMedium("ITSsddAir");
    TGeoVolume *virtualLayer3 = new TGeoVolume("ITSsddLayer3",
                                               virtualLayer3Shape, airSDD);
    TGeoVolume *lay3Ladder = CreateLay3Ladder();
    TGeoVolume *lay3Detectors = CreateLay3Detectors();

    Double_t dPhi = 360./fLay3Nladd;
    Double_t detBoxThickness = fLadWaferSep + 2*fWaferThickness;
    // placing virtual ladder and detectors volumes following ladder 
    // ordering convention
    char rotName[20];
    Int_t iLaddMin = 0;
    Int_t iLaddMax = fLay3Nladd;
    if ((fAddOnlyLadder3>=0)&&(fAddOnlyLadder3<fLay3Nladd)) {
        iLaddMin = fAddOnlyLadder3;
        iLaddMax = fAddOnlyLadder3+1;
    }
    for (Int_t iLadd = iLaddMin; iLadd < iLaddMax; iLadd++) {
        sprintf(rotName, "ITSsddLay3Ladd%i",iLadd);
        Double_t minRadiusLadBox = fLay3LaddShortRadius;
        if (iLadd%2 != 0) minRadiusLadBox = fLay3LaddLongRadius;
        minRadiusLadBox += ((TGeoBBox*)lay3Ladder->GetShape())->GetDY();
        
        TGeoCombiTrans *ctLadd = CreateCombiTrans(rotName,minRadiusLadBox,0,
                                                  -90+iLadd*dPhi);
        virtualLayer3->AddNode(lay3Ladder,iLadd,ctLadd);
        sprintf(rotName, "ITSsddLay3DetBox%i",iLadd);
        Double_t minRadiusDetBox = fLay3DetShortRadius;
        if (iLadd%2 != 0) minRadiusDetBox = fLay3DetLongRadius;
        minRadiusDetBox += detBoxThickness/2;
        TGeoCombiTrans *ctDet = CreateCombiTrans(rotName, minRadiusDetBox,0,
                                                 -90+iLadd*dPhi);
        virtualLayer3->AddNode(lay3Detectors, iLadd, ctDet);
    }
    virtualLayer3->SetVisibility(kFALSE);
    Moth->AddNode(virtualLayer3,1,0);
};
//________________________________________________________________________
TGeoVolume *AliITSv11GeometrySDD::CreateLay3Ladder() {
    // return a box volume containing the CF ladder
    TGeoVolume *laddSegment = CreateLadderSegment();
    TGeoBBox *ladBox = new TGeoBBox("ITSsddLadBox",
                               ((TGeoBBox*)laddSegment->GetShape())->GetDX(),
                               ((TGeoBBox*)laddSegment->GetShape())->GetDY(),
                                    //dX,dY = dX,dY of the segment
                               fLay3LadderLength/2);  
    TGeoMedium *airSDD = gGeoManager->GetMedium("ITSsddAir");
    TGeoVolume *virtualLadder = new TGeoVolume("ITSsddLadder",ladBox, airSDD);
    Double_t segmentLength = fSegmentLength;
    char transName[20];
    // placing virtual ladder segment following detector ordering convention
    //=======================================================================
    Int_t iSegmentMin = 1;
    Int_t iSegmentMax = fLay3Ndet;
    if (fAddOnlySegment) {
        iSegmentMin = fAddOnlySegment;
        iSegmentMax = fAddOnlySegment;
    }
    for (Int_t iSegment = iSegmentMin; iSegment <= iSegmentMax; iSegment++ ) {
        sprintf(transName, "ITSsddLay3LaddSeg%i", iSegment);
        Double_t segmentPos = segmentLength*(3-iSegment) + segmentLength/2;
        TGeoTranslation *segTr = new TGeoTranslation(transName,0,0,segmentPos);
        virtualLadder->AddNode(laddSegment, iSegment, segTr);
    }
    // putting virtual volume corresponding to the end of ladder
    //=======================================================================
    Double_t endLength = (fLay3LadderLength-fLay3Ndet*fSegmentLength)/2.;
    TGeoVolume *endLadder = CreateEndLadder( endLength );
    TGeoTranslation *endTrZPos = new TGeoTranslation("ITSsddEndTrZPos",0,0,
                                   fSegmentLength*(fLay3Ndet/2)+endLength/2.);
    //Euler rotation : about Z, then new X, then new Z
    TGeoRotation *endZNegRot = new TGeoRotation("",90, 180, -90);
    TGeoCombiTrans *endTrZNeg = new TGeoCombiTrans(0,0,
                    -fSegmentLength*(fLay3Ndet/2)-endLength/2., endZNegRot);
    if ((fAddOnlySegment==0)||(fAddOnlySegment==1))
        virtualLadder->AddNode(endLadder, 1, endTrZPos);
    if ((fAddOnlySegment==0)||(fAddOnlySegment==fLay3Ndet))
        virtualLadder->AddNode(endLadder, 2, endTrZNeg);
    virtualLadder->SetVisibility(kFALSE);
    return virtualLadder;
};
//________________________________________________________________________
TGeoArb8 *AliITSv11GeometrySDD::CreateLadderSide(Double_t dz,Double_t angle,
                          Double_t xSign,Double_t L, Double_t H, Double_t l) {
    // Create one half of the V shape corner of CF ladder
  
    TGeoArb8 *cfLaddSide = new TGeoArb8(dz);
    cfLaddSide->SetVertex(0, 0,0);
    cfLaddSide->SetVertex(1, 0, -H);
    cfLaddSide->SetVertex(2,xSign*(L*TMath::Sin(angle)-l*TMath::Cos(angle)),
                          -L*TMath::Cos(angle)-l*TMath::Sin(angle));
    cfLaddSide->SetVertex(3,xSign*L*TMath::Sin(angle),-L*TMath::Cos(angle));
    cfLaddSide->SetVertex(4, 0,0);
    cfLaddSide->SetVertex(5, 0, -H);
    cfLaddSide->SetVertex(6,xSign*(L*TMath::Sin(angle)-l*TMath::Cos(angle)),
                          -L*TMath::Cos(angle)-l*TMath::Sin(angle));
    cfLaddSide->SetVertex(7,xSign*L*TMath::Sin(angle),-L*TMath::Cos(angle));
    return cfLaddSide;
};
//________________________________________________________________________
void AliITSv11GeometrySDD::AddLadderCFstruct(Double_t dy, TGeoVolume* vol) {
    // fill a volume (segment) with the CF structure of a ladder

    TGeoMedium *carbonFiberLadderStruct = gGeoManager->GetMedium(
                                                        "ITSsddCarbonFiber");
    Double_t segmentLength = fSegmentLength;
    Double_t triangleHeight = fLadderHeight - fLadderBeamRadius;
    Double_t halfTheta = TMath::ATan( 0.5*fLadderWidth/triangleHeight );
    Double_t beta = (TMath::Pi()-2.*halfTheta)/4.;
    Double_t alpha = TMath::Pi()*3./4. - halfTheta/2.;
    Int_t colorCarbonFiber = 4;

    //--- The 3 V shape corners of the Carbon Fiber Ladder
    //--- the top V
    TGeoArb8 *cfLaddTop1 = CreateLadderSide(segmentLength/2., halfTheta, -1,
                                            fLadderLa, fLadderHa, fLadderl);
    TGeoVolume *cfLaddTopVol1 = new TGeoVolume("ITSsddCFladdTopCornerVol1",
                                           cfLaddTop1,carbonFiberLadderStruct);
    cfLaddTopVol1->SetLineColor(colorCarbonFiber);
    TGeoArb8 *cfLaddTop2 = CreateLadderSide(segmentLength/2., halfTheta, 1,
                                            fLadderLa, fLadderHa, fLadderl);
    TGeoVolume *cfLaddTopVol2 = new TGeoVolume("ITSsddCFladdTopCornerV2",
                                           cfLaddTop2,carbonFiberLadderStruct);
    cfLaddTopVol2->SetLineColor(colorCarbonFiber);
    TGeoTranslation *trTop1 = new TGeoTranslation(0, fLadderHeight/2+dy, 0);
    vol->AddNode(cfLaddTopVol1, 1, trTop1);
    vol->AddNode(cfLaddTopVol2, 1, trTop1);
    //--- the 2 side V
    TGeoArb8 *cfLaddSide1 = CreateLadderSide( segmentLength/2., beta, -1,
                                              fLadderLb, fLadderHb, fLadderl);
    TGeoVolume *cfLaddSideVol1 = new TGeoVolume("ITSsddCFladdSideCornerV1",
                                         cfLaddSide1,carbonFiberLadderStruct);
    cfLaddSideVol1->SetLineColor(colorCarbonFiber);
    TGeoArb8 *cfLaddSide2 = CreateLadderSide( segmentLength/2., beta, 1,
                                              fLadderLb, fLadderHb, fLadderl);
    TGeoVolume *cfLaddSideVol2 = new TGeoVolume("ITSsddCFladdSideCornerV2",
                                         cfLaddSide2,carbonFiberLadderStruct);
    cfLaddSideVol2->SetLineColor(colorCarbonFiber);

    Double_t dYTranslation = fLadderHeight/2. - 
                      0.5*fLadderWidth*TMath::Tan(beta) - fLadderBeamRadius;
    // because center of the triangle doesn't correspond to virtual vol. center
    Double_t distCenterSideDown =  0.5*fLadderWidth/TMath::Cos(beta);
    TGeoCombiTrans *ctSideR = CreateCombiTrans("", distCenterSideDown, 0,
                                               alpha*TMath::RadToDeg());
    AddTranslationTotCombiTrans(ctSideR, 0, -dYTranslation+dy, 0);
    TGeoCombiTrans *ctSideL = CreateCombiTrans("", distCenterSideDown,0,
                                               -alpha*TMath::RadToDeg());
    AddTranslationTotCombiTrans(ctSideL, 0, -dYTranslation+dy, 0);
    vol->AddNode(cfLaddSideVol1, 1, ctSideR);
    vol->AddNode(cfLaddSideVol2, 1, ctSideR);
    vol->AddNode(cfLaddSideVol1, 2, ctSideL);
    vol->AddNode(cfLaddSideVol2, 2, ctSideL);
    //--- The beams
    // Beams on the sides
    Double_t beamPhiPrime = TMath::ASin(1./TMath::Sqrt( (1+TMath::Sin(2*beta)*
                TMath::Sin(2*beta)/(TanD(fBeamSidePhi)*TanD(fBeamSidePhi))) ));
    if(GetDebug(1)) cout<<"Phi prime = "<<beamPhiPrime*TMath::RadToDeg()<<endl;
    Double_t beamLength = TMath::Sqrt( fLadderHeight*fLadderHeight/
                           (TMath::Sin(beamPhiPrime)*TMath::Sin(beamPhiPrime))
			       +fLadderWidth*fLadderWidth/4.)-fLadderLa/2-fLadderLb/2;
    TGeoTubeSeg *sideBeam = new TGeoTubeSeg(0,fLadderBeamRadius,beamLength/2.,
                                          0, 180);
    TGeoVolume *cfSideBeamVol = new TGeoVolume("ITSsddCFSideBeamVol", sideBeam,
                                               carbonFiberLadderStruct);
    cfSideBeamVol->SetLineColor(colorCarbonFiber);
    //Euler rotation : about Z, then new X, then new Z
    TGeoRotation *beamRot1 = new TGeoRotation("",90-2.*beta*TMath::RadToDeg(),
                                         -beamPhiPrime*TMath::RadToDeg(),-90);
    TGeoCombiTrans *beamTransf1 = new TGeoCombiTrans(0.5*triangleHeight*
    TMath::Tan(halfTheta),fLadderBeamRadius/2.+dy,-3*segmentLength/8,beamRot1);
    TGeoCombiTrans *beamTransf2 = new TGeoCombiTrans(0.5*triangleHeight*
       TMath::Tan(halfTheta),fLadderBeamRadius/2.+dy,segmentLength/8,beamRot1);
    TGeoRotation *beamRot2 = new TGeoRotation("",90-2.*beta*TMath::RadToDeg(),
					    beamPhiPrime*TMath::RadToDeg(), -90);
    TGeoCombiTrans *beamTransf3 = new TGeoCombiTrans(0.5*triangleHeight*
     TMath::Tan(halfTheta),fLadderBeamRadius/2.+dy,-segmentLength/8,beamRot2);
    TGeoCombiTrans *beamTransf4 = new TGeoCombiTrans(0.5*triangleHeight*
     TMath::Tan(halfTheta),fLadderBeamRadius/2.+dy,3*segmentLength/8,beamRot2);
    TGeoRotation *beamRot3 = new TGeoRotation("",90+2.*beta*TMath::RadToDeg(),
					    beamPhiPrime*TMath::RadToDeg(), -90);
    TGeoCombiTrans *beamTransf5 = new TGeoCombiTrans(-0.5*triangleHeight*
                                TMath::Tan(halfTheta),fLadderBeamRadius/2.+dy,
                                              -3*segmentLength/8,beamRot3);
    TGeoCombiTrans *beamTransf6 = new TGeoCombiTrans(-0.5*triangleHeight*
       TMath::Tan(halfTheta),fLadderBeamRadius/2.+dy,segmentLength/8,beamRot3);
    TGeoRotation *beamRot4 = new TGeoRotation("",90+2.*beta*TMath::RadToDeg(),
					    -beamPhiPrime*TMath::RadToDeg(), -90);
    TGeoCombiTrans *beamTransf7 = new TGeoCombiTrans(-0.5*triangleHeight*
      TMath::Tan(halfTheta),fLadderBeamRadius/2.+dy,-segmentLength/8,beamRot4);
    TGeoCombiTrans *beamTransf8 = new TGeoCombiTrans(-0.5*triangleHeight*
     TMath::Tan(halfTheta),fLadderBeamRadius/2.+dy,3*segmentLength/8,beamRot4);

    vol->AddNode(cfSideBeamVol, 1, beamTransf1);
    vol->AddNode(cfSideBeamVol, 2, beamTransf2);
    vol->AddNode(cfSideBeamVol, 3, beamTransf3);
    vol->AddNode(cfSideBeamVol, 4, beamTransf4);
    vol->AddNode(cfSideBeamVol, 5, beamTransf5);
    vol->AddNode(cfSideBeamVol, 6, beamTransf6);
    vol->AddNode(cfSideBeamVol, 7, beamTransf7);
    vol->AddNode(cfSideBeamVol, 8, beamTransf8);
    //--- Beams of the bottom
    TGeoTubeSeg *bottomBeam1 = new TGeoTubeSeg(0, fLadderBeamRadius,
					     fLadderWidth/2.-fLadderLb/3, 0, 180);
    TGeoVolume *bottomBeam1Vol = new TGeoVolume("ITSsddBottomBeam1Vol",
                                     bottomBeam1, carbonFiberLadderStruct);
    bottomBeam1Vol->SetLineColor(colorCarbonFiber);

    TGeoRotation *bottomBeamRot1 = new TGeoRotation("",90, 90, 90);
    TGeoCombiTrans *bottomBeamTransf1 = new TGeoCombiTrans(0,
                    -(fLadderHeight/2-fLadderBeamRadius)+dy,0, bottomBeamRot1);
    vol->AddNode(bottomBeam1Vol, 1, bottomBeamTransf1);
    TGeoTubeSeg *bottomBeam2 = new TGeoTubeSeg(0, fLadderBeamRadius,
					               fLadderWidth/2.-fLadderLb/3, 0, 90);
    TGeoVolume *bottomBeam2Vol = new TGeoVolume("ITSsddBottomBeam2Vol",
                                bottomBeam2, carbonFiberLadderStruct);
    bottomBeam2Vol->SetLineColor(colorCarbonFiber);
    TGeoCombiTrans *bottomBeamTransf2 = new TGeoCombiTrans(0,
                                     -(fLadderHeight/2-fLadderBeamRadius)+dy,
                                     -segmentLength/2, bottomBeamRot1);
    TGeoRotation *bottomBeamRot2 = new TGeoRotation("",-90, 90, -90);
    TGeoCombiTrans *bottomBeamTransf3 = new TGeoCombiTrans(0,
                                    -(fLadderHeight/2-fLadderBeamRadius)+dy,
							 segmentLength/2, bottomBeamRot2);
    vol->AddNode(bottomBeam2Vol, 1, bottomBeamTransf2);
    vol->AddNode(bottomBeam2Vol, 2, bottomBeamTransf3);
    TGeoTubeSeg *bottomBeam3 = new TGeoTubeSeg(0, fLadderBeamRadius,
			    0.5*fLadderWidth/SinD(fBottomBeamAngle)-fLadderLb/3,0,180);
    TGeoVolume *bottomBeam3Vol = new TGeoVolume("ITSsddBottomBeam3Vol",
                  bottomBeam3, carbonFiberLadderStruct);
    bottomBeam3Vol->SetLineColor(colorCarbonFiber);
    //bottomBeam3Vol->SetLineColor(2);
    // be careful on the next 2 beams : when "reading" from -z to +z and 
    // from the bottom of the ladder,
    // it should draw a Lambda, and not a V
    TGeoRotation *bottomBeamRot4 = new TGeoRotation("", -90, fBottomBeamAngle,
                                                    -90);
    TGeoCombiTrans *bottomBeamTransf4 = new TGeoCombiTrans(0, 
      -(fLadderHeight/2-fLadderBeamRadius)+dy,-segmentLength/4,bottomBeamRot4);
    TGeoRotation *bottomBeamRot5 = new TGeoRotation("",-90,-fBottomBeamAngle,
                                                    -90);
    TGeoCombiTrans *bottomBeamTransf5 = new TGeoCombiTrans(0,
     -(fLadderHeight/2-fLadderBeamRadius)+dy,segmentLength/4, bottomBeamRot5);
    vol->AddNode(bottomBeam3Vol, 1, bottomBeamTransf4);
    vol->AddNode(bottomBeam3Vol, 2, bottomBeamTransf5);
};
//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreateHybrid() {
    // return a box containing the front-end hybrid
    TGeoBBox *hybridBox = new TGeoBBox("ITSsddHybridBox",fHybridWidth/2,0.3/2,
                                       fHybridLength/2); // <===== 0.3 tempo
    TGeoMedium *airSDD = gGeoManager->GetMedium("ITSsddAir");
    TGeoVolume *VirtualHybrid =new TGeoVolume("ITSsddHybrid",hybridBox,airSDD);

    VirtualHybrid->SetVisibility(kFALSE);
    return VirtualHybrid;
};
//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreateLadderSegment() {
    // Return a box volume containing a segment of a ladder.

    TGeoMedium *airSDD = gGeoManager->GetMedium("ITSsddAir");

    Double_t tDY = -0.5; //space left on top of the ladder
    Double_t segmentLength = fSegmentLength;

    TGeoBBox *segBox = new TGeoBBox("ITSsddSegBox",
                                    (fLadderWidth+fPinSuppWidth)/2,
                                    fLadderHeight/2+TMath::Abs(tDY),
                                    segmentLength/2);
    TGeoVolume *VirtualSeg = new TGeoVolume("ITSsddSegment",segBox, airSDD);
    // the carbon fiber structure :
    AddLadderCFstruct(tDY, VirtualSeg);
    // the 2 hybrids :
    TGeoVolume *Hybrid = CreateHybrid();
    Double_t hybDx = ((TGeoBBox*)Hybrid->GetShape())->GetDX();
    Double_t hybDy = ((TGeoBBox*)Hybrid->GetShape())->GetDY();
    Double_t halfDiag = TMath::Sqrt(hybDx*hybDx+hybDy*hybDy);
    Double_t theta = 50;  //angle in the transverse plane of hybrids 
    // laying on the CF structure  // <========== 50° tempo
    Double_t thetaPrim = theta*TMath::DegToRad()+TMath::ATan(hybDy/hybDx);
    TGeoRotation *HybridRot1 = new TGeoRotation("",   90-theta, 0, 0);
    TGeoRotation *HybridRot2 = new TGeoRotation("", -(90-theta), 0, 0);
    Double_t hybZ = 0;
    Double_t hybY = fLadderHeight/2-halfDiag*TMath::Cos(thetaPrim)+tDY;
    Double_t hybX = halfDiag*TMath::Sin(thetaPrim);
    TGeoCombiTrans *HybTransf1 = new TGeoCombiTrans(-hybX,hybY,hybZ,
                                                    HybridRot1);
    VirtualSeg->AddNode(Hybrid, 1, HybTransf1);
    TGeoCombiTrans *HybTransf2 = new TGeoCombiTrans(hybX,hybY,hybZ,HybridRot2);
    VirtualSeg->AddNode(Hybrid, 2, HybTransf2);
    //**********************************
    TGeoVolume *pinSupport = CreatePinSupport(0);
    VirtualSeg->AddNode(pinSupport, 1, 0);
    VirtualSeg->SetVisibility(kFALSE);
    return VirtualSeg;
};
//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreatePinSupport(Double_t rotY) {
    //medium = ryton ? To code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    TGeoMedium *rytonSDD = gGeoManager->GetMedium("ITSsddStaselite4411w");
    TGeoCone *cone = new TGeoCone("ITSsddPinSuppCone",fPinSuppHeight/2.,
                                  0,fPinSuppRmax,0,fPinSuppRmax-
                                  fPinSuppHeight*TanD(fPinSuppConeAngle) );
    TGeoBBox *tong = new TGeoBBox("ITSsddPinSuppTong",fPinSuppRmax,
                                  fPinSuppLength/2.,fPinSuppThickness/2.);
    TGeoTube *hole = new TGeoTube("ITSsddPinSuppHole",0,fPinR,
                                  fPinSuppHeight/2.);
    rotY = 0.0; // Remove compiler warning.
    if(GetDebug(3)){
        cone->InspectShape();
        tong->InspectShape();
        hole->InspectShape();
    }

    TGeoTranslation *tongTrans = new TGeoTranslation("ITSsddPinSuppTongTr",0,
                   fPinSuppLength/2.,-fPinSuppHeight/2.+fPinSuppThickness/2.);
    tongTrans->RegisterYourself();
    TGeoCompositeShape *pinSupportShape = new TGeoCompositeShape(
               "ITSssdPinSupportShape","(ITSsddPinSuppCone+"
               "ITSsddPinSuppTong:ITSsddPinSuppTongTr)-ITSsddPinSuppHole");

    TGeoVolume *pinSupport = new TGeoVolume("ITSssdPinSupport",pinSupportShape,
                                            rytonSDD);
    return pinSupport;
    // include the pin itself !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // le centre de rotation, translation est le centre du cone
    // (puisque celui-ci n'a pas bouger dans la composition)
};
//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreateEndLadder(Double_t length) {
    // Return a box volume containing a end of a CF ladder.
    Double_t tDY = -0.5; //space left on top of the ladder
    TGeoMedium *airSDD = gGeoManager->GetMedium("ITSsddAir");
    TGeoMedium *carbonFiberLadderStruct = gGeoManager->GetMedium(
                                                       "ITSsddCarbonFiber");
    Double_t segmentLength = fSegmentLength;
    Double_t topCornerLength = fSegmentLength/2.-fLay4LaddTopCornerEnd;

    TGeoBBox *endBox = new TGeoBBox("ITSsddEndLaddBox",
                                    //(fLadderWidth+fPinSuppWidth)/2,
                                    (fLadderWidth)/2,
                                    fLadderHeight/2+TMath::Abs(tDY),
                                    length/2);
    TGeoVolume *virtualEnd = new TGeoVolume("ITSsddEnd",endBox, airSDD);
    //**********************************
    // coding real matter :
    //**********************************
    Double_t triangleHeight = fLadderHeight - fLadderBeamRadius;
    Double_t halfTheta = TMath::ATan( 0.5*fLadderWidth/triangleHeight );
    Double_t beta = (TMath::Pi()-2.*halfTheta)/4.;
    Double_t alpha = TMath::Pi()*3./4. - halfTheta/2.;
    Int_t colorCarbonFiber = 4;
    //--- The 3 V shape corners of the Carbon Fiber Ladder
    //--- the top V
    TGeoArb8 *cfLaddTop1 = CreateLadderSide( topCornerLength/2., halfTheta, -1,
                                             fLadderLa, fLadderHa, fLadderl);
    TGeoVolume *cfLaddTopVol1 = new TGeoVolume("ITSsddCFladdTopCornerVol1",
                                          cfLaddTop1,carbonFiberLadderStruct);
    cfLaddTopVol1->SetLineColor(colorCarbonFiber);
    TGeoArb8 *cfLaddTop2 = CreateLadderSide( topCornerLength/2., halfTheta, 1,
                                             fLadderLa, fLadderHa, fLadderl);
    TGeoVolume *cfLaddTopVol2 = new TGeoVolume("ITSsddCFladdTopCornerV2",
                                        cfLaddTop2,carbonFiberLadderStruct);
    cfLaddTopVol2->SetLineColor(colorCarbonFiber);
    TGeoTranslation *trTop1 = new TGeoTranslation(0, fLadderHeight/2+tDY,
                                               -(length-topCornerLength)/2.);
    virtualEnd->AddNode(cfLaddTopVol1, 1, trTop1);
    virtualEnd->AddNode(cfLaddTopVol2, 1, trTop1);
    //--- the 2 side V
    TGeoArb8 *cfLaddSide1 = CreateLadderSide( length/2., beta, -1,
                                              fLadderLb, fLadderHb, fLadderl);
    TGeoVolume *cfLaddSideVol1 = new TGeoVolume("ITSsddCFladdSideCornerV1",
                                          cfLaddSide1,carbonFiberLadderStruct);
    cfLaddSideVol1->SetLineColor(colorCarbonFiber);
    TGeoArb8 *cfLaddSide2 = CreateLadderSide( length/2., beta, 1,
                                              fLadderLb, fLadderHb, fLadderl);
    TGeoVolume *cfLaddSideVol2 = new TGeoVolume("ITSsddCFladdSideCornerV2",
                                          cfLaddSide2,carbonFiberLadderStruct);
    cfLaddSideVol2->SetLineColor(colorCarbonFiber);
    Double_t dYTranslation = ( fLadderHeight/2. - 0.5*fLadderWidth*
                               TMath::Tan(beta)- fLadderBeamRadius );
    // because center of the triangle doesn't correspond to virtual vol. center
    Double_t distCenterSideDown =  0.5*fLadderWidth/TMath::Cos(beta);
    TGeoCombiTrans *ctSideR = CreateCombiTrans("", distCenterSideDown, 0,
                                               alpha*TMath::RadToDeg());
    AddTranslationTotCombiTrans(ctSideR, 0, -dYTranslation+tDY, 0);
    TGeoCombiTrans *ctSideL = CreateCombiTrans("", distCenterSideDown, 0, 
                                               -alpha*TMath::RadToDeg());
    AddTranslationTotCombiTrans(ctSideL, 0, -dYTranslation+tDY, 0);
    virtualEnd->AddNode(cfLaddSideVol1, 1, ctSideR);
    virtualEnd->AddNode(cfLaddSideVol2, 1, ctSideR);
    virtualEnd->AddNode(cfLaddSideVol1, 2, ctSideL);
    virtualEnd->AddNode(cfLaddSideVol2, 2, ctSideL);
    //--- The beams
    // Beams on the sides
    Double_t beamPhiPrime = TMath::ASin(1./TMath::Sqrt( (1+TMath::Sin(2*beta)*
                TMath::Sin(2*beta)/(TanD(fBeamSidePhi)*TanD(fBeamSidePhi))) ));
    if(GetDebug(3)) cout<<"Phi prime = "<<beamPhiPrime*TMath::RadToDeg()<<endl;
     Double_t beamLength = TMath::Sqrt( fLadderHeight*fLadderHeight/
				     (TMath::Sin(beamPhiPrime)*TMath::Sin(beamPhiPrime))
				   +fLadderWidth*fLadderWidth/4.)-fLadderLa/2-fLadderLb/2;
     TGeoTubeSeg *sideBeam = new TGeoTubeSeg(0, fLadderBeamRadius,
                                             beamLength/2., 0, 180);
     TGeoVolume *cfSideBeamVol = new TGeoVolume("ITSsddCFSideBeamVol",
                                      sideBeam, carbonFiberLadderStruct);
     cfSideBeamVol->SetLineColor(colorCarbonFiber);
     //Euler rotation : about Z, then new X, then new Z
     TGeoRotation *beamRot1 = new TGeoRotation("",90-2.*beta*TMath::RadToDeg(),
					                -beamPhiPrime*TMath::RadToDeg(), -90);
     TGeoCombiTrans *beamTransf1 = new TGeoCombiTrans(0.5*triangleHeight*
                   TMath::Tan(halfTheta),fLadderBeamRadius/2.+tDY,
						   -length/2+segmentLength/8,beamRot1);
     TGeoRotation *beamRot2 = new TGeoRotation("",90-2.*beta*TMath::RadToDeg(),
					    beamPhiPrime*TMath::RadToDeg(), -90);
     TGeoCombiTrans *beamTransf3 = new TGeoCombiTrans(0.5*triangleHeight*
                     TMath::Tan(halfTheta),fLadderBeamRadius/2.+tDY,
						   -length/2+3*segmentLength/8,beamRot2);
     TGeoRotation *beamRot3 = new TGeoRotation("",90+2.*beta*TMath::RadToDeg(),
					    beamPhiPrime*TMath::RadToDeg(), -90);
     TGeoCombiTrans *beamTransf5 = new TGeoCombiTrans(-0.5*triangleHeight*
                          TMath::Tan(halfTheta),fLadderBeamRadius/2.+tDY,
						   -length/2+segmentLength/8,beamRot3);
     TGeoRotation *beamRot4 = new TGeoRotation("",90+2.*beta*TMath::RadToDeg(),
					    -beamPhiPrime*TMath::RadToDeg(), -90);
     TGeoCombiTrans *beamTransf7 = new TGeoCombiTrans(-0.5*triangleHeight*
            TMath::Tan(halfTheta),fLadderBeamRadius/2.+tDY,
						   -length/2+3*segmentLength/8,beamRot4);
     virtualEnd->AddNode(cfSideBeamVol, 1, beamTransf1);
     virtualEnd->AddNode(cfSideBeamVol, 2, beamTransf3);
     virtualEnd->AddNode(cfSideBeamVol, 3, beamTransf5);
     virtualEnd->AddNode(cfSideBeamVol, 4, beamTransf7);
     //--- Beams of the bottom
     TGeoTubeSeg *bottomBeam1 = new TGeoTubeSeg(0, fLadderBeamRadius,
					     fLadderWidth/2.-fLadderLb/3, 0, 180);
     TGeoVolume *bottomBeam1Vol = new TGeoVolume("ITSsddBottomBeam1Vol",
                                      bottomBeam1, carbonFiberLadderStruct);
     bottomBeam1Vol->SetLineColor(colorCarbonFiber);
     TGeoRotation *bottomBeamRot1 = new TGeoRotation("",90, 90, 90);
     TGeoCombiTrans *bottomBeamTransf1 = new TGeoCombiTrans(0,
                            -(fLadderHeight/2-fLadderBeamRadius)+tDY,
						 -length/2+fSegmentLength/2, bottomBeamRot1);
     virtualEnd->AddNode(bottomBeam1Vol, 1, bottomBeamTransf1);
     TGeoTubeSeg *bottomBeam2 = new TGeoTubeSeg(0, fLadderBeamRadius,
					     fLadderWidth/2.-fLadderLb/3, 0, 90);
     TGeoVolume *bottomBeam2Vol = new TGeoVolume("ITSsddBottomBeam2Vol",
                                   bottomBeam2, carbonFiberLadderStruct);
     bottomBeam2Vol->SetLineColor(colorCarbonFiber);
     TGeoCombiTrans *bottomBeamTransf2 = new TGeoCombiTrans(0,
           -(fLadderHeight/2-fLadderBeamRadius)+tDY,-length/2, bottomBeamRot1);
     virtualEnd->AddNode(bottomBeam2Vol, 1, bottomBeamTransf2);
     //**********************************
     virtualEnd->SetVisibility(kFALSE);
     return virtualEnd;
};
//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreateSDDsensor() {
    // return a box containing the SDD sensor

    TGeoBBox *sensorBox = new TGeoBBox("ITSsddSensorBox",
				  fWaferWidth/2, fWaferThickness/2,fWaferLength/2);
    TGeoMedium *airSDD = gGeoManager->GetMedium("ITSsddAir");
    TGeoVolume *virtualSensor =new TGeoVolume("ITSsddSensor",sensorBox,airSDD);

    //virtualSensor->SetVisibility(kFALSE);
    return virtualSensor;
};
//________________________________________________________________________
TGeoVolume *AliITSv11GeometrySDD::CreateLay3Detectors() {
    // return a box volume containing the detectors

    TGeoBBox *detBox = new TGeoBBox("ITSssdDetBox3",
                                    fWaferWidth/2,
                                    (fLadWaferSep + 2*fWaferThickness)/2,
                             fLay3LadderLength*((fLay3Ndet-0.5)/fLay3Ndet)/2);
    TGeoMedium *airSDD = gGeoManager->GetMedium("ITSsddAir");
    TGeoVolume *virtualDet = new TGeoVolume("ITSsddDet3",detBox, airSDD);
    TGeoVolume *vSDD = CreateSDDsensor();
    char transName[30];
    for (Int_t i=0; i<fLay3Ndet; i++) {
        Double_t localZ = fLay3sensorZPos[i];
        Double_t localY = fLadWaferSep/2+fWaferThickness/2;
        if (i%2!=0)
            localY = -localY;
        sprintf(transName, "ITSsddLay3SensorPos%i", i+1);
        TGeoTranslation *sensorPos = new TGeoTranslation(transName,0,localY,
                                                         localZ);
        virtualDet->AddNode(vSDD, i+1, sensorPos);
    }
    virtualDet->SetVisibility(kFALSE);
    return virtualDet;
};
//________________________________________________________________________
void AliITSv11GeometrySDD::Layer4(TGeoVolume *Moth) {
    // Insert the layer 4 in the mother volume. This is a virtual volume
    // containing ladders of layer 4 and the supporting rings

    TGeoTube *virtualLayer4Shape =new TGeoTube("ITSsddLayer4Shape",
                                        fLay4Rmin,fLay4Rmax,fLay4Length*0.5);
    TGeoMedium *airSDD = gGeoManager->GetMedium("ITSsddAir");
    TGeoVolume *virtualLayer4 = new TGeoVolume("ITSsddLayer4",
                                          virtualLayer4Shape, airSDD);
    TGeoVolume *lay4Ladder = CreateLay4Ladder();
    TGeoVolume *lay4Detectors = CreateLay4Detectors();
    Double_t dPhi = 360./fLay4Nladd;
    Double_t detBoxThickness = fLadWaferSep + 2*fWaferThickness;
    // placing virtual ladder and detectors volumes following ladder 
    // ordering convention
    char rotName[20];
    Int_t iLaddMin = 0;
    Int_t iLaddMax = fLay4Nladd;
    if ((fAddOnlyLadder4>=0)&&(fAddOnlyLadder4<fLay4Nladd)) {
        iLaddMin = fAddOnlyLadder4;
        iLaddMax = fAddOnlyLadder4+1;
    }
    for (Int_t iLadd = iLaddMin; iLadd < iLaddMax; iLadd++) {
        sprintf(rotName, "ITSsddLay4Ladd%i",iLadd);
        Double_t minRadiusLadBox = fLay4LaddShortRadius;
        if (iLadd%2 != 0)
            minRadiusLadBox = fLay4LaddLongRadius;
        minRadiusLadBox += ((TGeoBBox*)lay4Ladder->GetShape())->GetDY();
        TGeoCombiTrans *ctLadd = CreateCombiTrans(rotName, minRadiusLadBox,
                                                  0, -90+iLadd*dPhi);
        virtualLayer4->AddNode(lay4Ladder, iLadd, ctLadd);
        sprintf(rotName, "ITSsddLay4DetBox%i",iLadd);
        Double_t minRadiusDetBox = fLay4DetShortRadius;
        if (iLadd%2 != 0)
            minRadiusDetBox = fLay4DetLongRadius;
        minRadiusDetBox += detBoxThickness/2;
        TGeoCombiTrans *ctDet = CreateCombiTrans(rotName, minRadiusDetBox,
                                                 0, -90+iLadd*dPhi);
        virtualLayer4->AddNode(lay4Detectors, iLadd, ctDet);
    }
    virtualLayer4->SetVisibility(kFALSE);
    Moth->AddNode(virtualLayer4,1,0);
};
//________________________________________________________________________
TGeoVolume *AliITSv11GeometrySDD::CreateLay4Ladder() {
    // return a box volume containing the CF ladder and all pieces on it

    TGeoVolume *laddSegment = CreateLadderSegment();
    TGeoBBox *ladBox = new TGeoBBox("ITSsddLadBox",
                               ((TGeoBBox*)laddSegment->GetShape())->GetDX(),
				  ((TGeoBBox*)laddSegment->GetShape())->GetDY(),
				  //dX,dY = dX,dY of the segment
				  fLay4LadderLength/2);  
    TGeoMedium *airSDD = gGeoManager->GetMedium("ITSsddAir");
    TGeoVolume *virtualLadder = new TGeoVolume("ITSsddLadder",ladBox, airSDD);
    Double_t segmentLength = fSegmentLength;
    char transName[20];
    // placing virtual ladder segment following detector ordering convention
    //=======================================================================
    Int_t iSegmentMin = 1;
    Int_t iSegmentMax = fLay4Ndet;
    if (fAddOnlySegment) {
        iSegmentMin = fAddOnlySegment;
        iSegmentMax = fAddOnlySegment;
    }
    for (Int_t iSegment = iSegmentMin; iSegment <= iSegmentMax; iSegment++ ) {
        sprintf(transName, "ITSsddLay4LaddSeg%i", iSegment);
        Double_t segmentPos = segmentLength*(fLay4Ndet/2-iSegment) +
            segmentLength/2;
        TGeoTranslation *segTr = new TGeoTranslation(transName,0,0,segmentPos);
        virtualLadder->AddNode(laddSegment, iSegment, segTr);
    }
    // putting virtual volume corresponding to the end of ladder
    //=======================================================================
    Double_t endLength = (fLay4LadderLength-fLay4Ndet*fSegmentLength)/2.;
    TGeoVolume *endLadder = CreateEndLadder( endLength );
    TGeoTranslation *endTrZPos =
        new TGeoTranslation("ITSsddEndTrZPos",0,0,
                            fSegmentLength*(fLay4Ndet/2)+endLength/2.);
    //Euler rotation : about Z, then new X, then new Z
    TGeoRotation *endZNegRot = new TGeoRotation("",90, 180, -90);
    TGeoCombiTrans *endTrZNeg =
        new TGeoCombiTrans(0,0,-fSegmentLength*(fLay4Ndet/2)-endLength/2.,
                           endZNegRot);
    if ((fAddOnlySegment==0)||(fAddOnlySegment==1))
        virtualLadder->AddNode(endLadder, 1, endTrZPos);
    if ((fAddOnlySegment==0)||(fAddOnlySegment==fLay4Ndet))
        virtualLadder->AddNode(endLadder, 2, endTrZNeg);
    virtualLadder->SetVisibility(kFALSE);
    return virtualLadder;
};
//________________________________________________________________________
TGeoVolume *AliITSv11GeometrySDD::CreateLay4Detectors() {
    // return a box volume containing the detectors

    TGeoBBox *detBox = new TGeoBBox("ITSssdDetBox4",
                                    fWaferWidth/2,
                                    (fLadWaferSep + 2*fWaferThickness)/2,
                             fLay4LadderLength*((fLay4Ndet-0.5)/fLay4Ndet)/2);
    TGeoMedium *airSDD = gGeoManager->GetMedium("ITSsddAir");
    TGeoVolume *virtualDet = new TGeoVolume("ITSsddDet4",detBox, airSDD);
    TGeoVolume *vSDD = CreateSDDsensor();
    char transName[30];
    for (Int_t i=0; i<fLay4Ndet; i++) {
        Double_t localZ = fLay4sensorZPos[i];
        Double_t localY = fLadWaferSep/2+fWaferThickness/2;
        if (i%2!=0)
            localY = -localY;
        sprintf(transName, "ITSsddLay4SensorPos%i", i+1);
        TGeoTranslation *sensorPos = new TGeoTranslation(transName,0,localY,
                                                         localZ);
        virtualDet->AddNode(vSDD, i+1, sensorPos);
    }
    virtualDet->SetVisibility(kFALSE);
    return virtualDet;
};
