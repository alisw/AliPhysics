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
#include <TVectorD.h>

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
#include <TGeoNode.h>

#include "AliITSv11Geometry.h"
#include "AliITSv11GeometrySDD.h"

extern TGeoManager* gGeoManager;

ClassImp(AliITSv11GeomSDDcable)

//----------------------------------------------------------------------
AliITSv11GeomSDDcable::~AliITSv11GeomSDDcable() { 
  if (fInitialNode) delete fInitialNode; };


//----------------------------------------------------------------------
void AliITSv11GeomSDDcable::AddCheckPoint( TGeoVolume *vol, Int_t iCheckPt,
					   Double_t *coord)
{
  if (iCheckPt>=fVolumeArray.GetEntriesFast()) {
    fVolumeArray.AddLast(vol);
    TVectorD *point = new TVectorD(3,coord);
    fPointArray.AddLast(point);

  } else if ((iCheckPt >= 0)&&(iCheckPt < fVolumeArray.GetEntriesFast())) {
    fVolumeArray.AddAt(vol, iCheckPt);
    TVectorD *point = new TVectorD(3,coord);
    fPointArray.AddAt(point, iCheckPt);
  };
};


//----------------------------------------------------------------------
Int_t AliITSv11GeomSDDcable::
GetCheckPoint( Int_t iCheckPt, Int_t iOccur, Int_t motherLevel,
	       Double_t *coord ) {
// Get the coordinate of the #i check point of the #iOccur occurrence of
// its containing volume in the node tree. Coord. are given in the coordinate
// system of the #motherLevel mother level of this volume
  
  if (iCheckPt >= fVolumeArray.GetEntriesFast()) return kFALSE;
  fCurrentVol = (TGeoVolume *) fVolumeArray.UncheckedAt(iCheckPt);

  TGeoNode *mainNode;
  if (fInitialNode==0) {
    TObjArray *nodes = gGeoManager->GetListOfNodes();
    if (nodes->GetEntriesFast()==0) return kFALSE;
    mainNode = (TGeoNode *) nodes->UncheckedAt(0);
  } else {
    mainNode = fInitialNode;
  };

  for (Int_t i=0; i<50; i++) fNodeInd[i]=-1;  
  Int_t currentOccur = 0;

  // loop to get the volume position in the tree of nodes
  while ( CheckDaughter(0, mainNode) && (currentOccur!=iOccur) ) {
    currentOccur++;
    Int_t maxLevel = 0;
    while (fNodeInd[maxLevel]!=-1) maxLevel++;
    fNodeInd[maxLevel-1]++;
  };

  Int_t maxLevel = 0;
  while (fNodeInd[maxLevel]!=-1) maxLevel++;
  if (maxLevel==0) return kFALSE;

  if (motherLevel>maxLevel) motherLevel = maxLevel;
  TGeoNode *pathNode[50];
  pathNode[0] = mainNode;
  for (Int_t i=0; i<motherLevel; i++)
    pathNode[i+1] = pathNode[i]->GetDaughter(fNodeInd[i]);


  TVectorD *coordVector = (TVectorD *)fPointArray.UncheckedAt(iCheckPt);
  Double_t localCoord[3], globalCoord[3];
  localCoord[0] = (coordVector->GetMatrixArray())[0];
  localCoord[1] = (coordVector->GetMatrixArray())[1];
  localCoord[2] = (coordVector->GetMatrixArray())[2];

//   cout << '(' << localCoord[0] << ',' << localCoord[1] << ','
//        << localCoord[2] << ") "<< endl;

  for (Int_t i=motherLevel; i>=0; i--) {
//     cout << pathNode[i]->GetName() << " > ";
    pathNode[i]->GetMatrix()->LocalToMaster(localCoord, globalCoord);
    localCoord[0] = globalCoord[0];
    localCoord[1] = globalCoord[1];
    localCoord[2] = globalCoord[2];
//     cout << '(' << localCoord[0] << ',' << localCoord[1] << ','
// 	 << localCoord[2] << ") "<< endl;
  };
  coord[0] = globalCoord[0];
  coord[1] = globalCoord[1];
  coord[2] = globalCoord[2];

  return kTRUE;
};


//----------------------------------------------------------------------
void AliITSv11GeomSDDcable::SetInitialNode(TGeoVolume *vol) {
  if (fInitialNode) delete fInitialNode;
  fInitialNode = new TGeoNodeMatrix(vol,0);
  fInitialNode->SetName("nodeInConstruction");
};


//----------------------------------------------------------------------
void AliITSv11GeomSDDcable:: ResetInitialNode() {
  if (fInitialNode) delete fInitialNode;
  fInitialNode = 0;
};


//----------------------------------------------------------------------
bool AliITSv11GeomSDDcable::CheckDaughter(Int_t i, TGeoNode* node) {
// Search where is the current volume in the tree of nodes
// stop each time it find the pointer of the current volume
// the path is recorded in fNodeInd[]
// !!! recursiv function !!!

  Int_t j = fNodeInd[i];
  if (node->GetVolume()==fCurrentVol) return kTRUE;
  TObjArray *array = node->GetNodes();
  if (array) {
    Int_t nDaughters = array->GetEntriesFast();
    if (j==-1) j++;
    while (j<nDaughters) {
      TGeoNode *subNode = (TGeoNode *) array->UncheckedAt(j);
      //cout << "level " <<  i << "  " <<  subNode->GetName()<< endl;
      fNodeInd[i] = j;
      if (CheckDaughter(i+1, subNode)) return kTRUE;
      j++;
    };
    fNodeInd[i] = -1;
  };
  return kFALSE;
};



ClassImp(AliITSv11GeomSDDcableNap)

//----------------------------------------------------------------------
//----------------------------------------------------------------------
AliITSv11GeomSDDcableNap::
AliITSv11GeomSDDcableNap(Double_t width, Double_t thick):AliITSv11GeomSDDcable()
{
  fWidth = width;
  fThick = thick;
};


#include <TGeoSphere.h>
//----------------------------------------------------------------------
Int_t AliITSv11GeomSDDcableNap::
CreateAndInsertCableSegment(Int_t p1, Int_t p2, TGeoVolume *motherVol) {

  Double_t coord1[3], coord2[3];
  GetCheckPoint(p1, 0, 0, coord1);
  GetCheckPoint(p2, 0, 0, coord2);
  cout << coord1[0] << ','<< coord1[1] << ','<< coord1[2] << endl;
  cout << coord2[0] << ','<< coord2[1] << ','<< coord2[2] << endl;
  coord2[0] = 5;
  coord2[1] = 0;

  Double_t cx = (coord1[0]+coord2[0])/2;
  Double_t cy = (coord1[1]+coord2[1])/2;
  Double_t cz = (coord1[2]+coord2[2])/2;
  Double_t dx = coord2[0]-coord1[0];
  Double_t dy = coord2[1]-coord1[1];
  Double_t dz = coord2[2]-coord1[2];

  Double_t length = TMath::Sqrt(dx*dx+dy*dy+dz*dz);
  TGeoBBox *cableSeg = new TGeoBBox("", fWidth/2,fThick/2,length/2); 

  TGeoMedium *airSDD = gGeoManager->GetMedium("ITSsddAir");
  TGeoVolume *vCableSeg = new TGeoVolume("vCableSeg",cableSeg, airSDD);

  TGeoTranslation *trans = new TGeoTranslation("",cx, cy, cz);
  TGeoRotation*rot = new TGeoRotation("",-TMath::ATan(dx/dy)*TMath::RadToDeg(),
				      -TMath::ATan(dy/dz)*TMath::RadToDeg(),0);
  //rot->RotateZ( -TMath::ATan(dx/dy)*TMath::RadToDeg() );
  //rot->RotateX( -TMath::ATan(dy/dz)*TMath::RadToDeg() );
  cout << TMath::ATan(dy/dz)*TMath::RadToDeg() << " ######### " 
       << TMath::ATan(dx/dy)*TMath::RadToDeg() << endl;

  TGeoCombiTrans *combi = new TGeoCombiTrans(*trans, *rot);
  motherVol->AddNode(vCableSeg, 1, combi);



  //------------
  TGeoSphere *sphere = new TGeoSphere(0, 0.1);
  TGeoVolume *vSphere = new TGeoVolume("", sphere, airSDD);
  TGeoTranslation *trC = new TGeoTranslation("", cx, cy, cz);
  TGeoTranslation *tr1 = new TGeoTranslation("",coord1[0],coord1[1],coord1[2]);
  TGeoTranslation *tr2 = new TGeoTranslation("",coord2[0],coord2[1],coord2[2]);
  motherVol->AddNode(vSphere, 1, trC);
  motherVol->AddNode(vSphere, 2, tr1);
  motherVol->AddNode(vSphere, 3, tr2);

  return kTRUE;
};



//________________________________________________________________________


const Double_t AliITSv11GeometrySDD::fSegmentLength    = 37.2*2*fgkmm;
const Double_t AliITSv11GeometrySDD::fLadderWidth      = 50.0*fgkmm;
const Double_t AliITSv11GeometrySDD::fLadderHeight     = 30.0*fgkmm;
const Double_t AliITSv11GeometrySDD::fLadderBeamRadius =  0.6*fgkmm;
const Double_t AliITSv11GeometrySDD::fLadderLa         =  3.*fgkmm;
const Double_t AliITSv11GeometrySDD::fLadderHa         =  0.6*fgkmm; //total pifometer
const Double_t AliITSv11GeometrySDD::fLadderLb         =  3.7*fgkmm;
const Double_t AliITSv11GeometrySDD::fLadderHb         =  0.6*fgkmm; //total pifometer
const Double_t AliITSv11GeometrySDD::fLadderl          =  0.25*fgkmm;

const Double_t AliITSv11GeometrySDD::fBottomBeamAngle  = 56.5;
const Double_t AliITSv11GeometrySDD::fBeamSidePhi      = 65;

const Double_t AliITSv11GeometrySDD::fWaferThickness   = 0.3*fgkmm;
const Double_t AliITSv11GeometrySDD::fWaferWidth       = 72.5*fgkmm;
const Double_t AliITSv11GeometrySDD::fWaferLength      = 87.6*fgkmm;

const Double_t AliITSv11GeometrySDD::fHybridLength     = 65*fgkmm;
const Double_t AliITSv11GeometrySDD::fHybridWidth      = 41*fgkmm;  // pifometer
const Double_t AliITSv11GeometrySDD::fHybridThBridgeThick = 0.25*fgkmm; // pifometer

const Double_t AliITSv11GeometrySDD::fLadWaferSep      = 2*fgkmm;
const Double_t AliITSv11GeometrySDD::fPinSuppWidth     = 2.5*fgkmm; // pifometer
const Double_t AliITSv11GeometrySDD::fPinSuppHeight    = 2.*fgkmm; // pifometer
const Double_t AliITSv11GeometrySDD::fPinSuppRmax      = 2.5/2.*fgkmm;
const Double_t AliITSv11GeometrySDD::fPinR             = 1.5/2.*fgkmm;
const Double_t AliITSv11GeometrySDD::fPinSuppLength    = 5.*fgkmm;
const Double_t AliITSv11GeometrySDD::fPinSuppThickness = 0.5*fgkmm;
const Double_t AliITSv11GeometrySDD::fPinSuppConeAngle = 4;
const Double_t AliITSv11GeometrySDD::fPinDXminOnSensor = (39./2.)*fgkmm; //placement of pins on sensor
const Double_t AliITSv11GeometrySDD::fPinPinDDXOnSensor = 3*fgkmm;
const Double_t AliITSv11GeometrySDD::fPinDYOnSensor    = (52.5/2.)*fgkmm;

// parameters from ALR-0752/3
const Double_t AliITSv11GeometrySDD::fCoolPipeSuppHeight    = 3.2*fgkmm;  
const Double_t AliITSv11GeometrySDD::fCoolPipeSuppMaxLength = 14*fgkmm;
const Double_t AliITSv11GeometrySDD::fCoolPipeSuppWidthExt  = 0.4*fgkmm;
const Double_t AliITSv11GeometrySDD::fCoolPipeSuppWidthIn   = 0.65*fgkmm;
const Double_t AliITSv11GeometrySDD::fCoolPipeSuppHoleDiam  = 2*fgkmm;
const Double_t AliITSv11GeometrySDD::fCoolPipeSuppFulWidth  = 5.15*fgkmm;
const Double_t AliITSv11GeometrySDD::fCoolPipeSuppTongW     = 0.8*fgkmm; 
const Double_t AliITSv11GeometrySDD::fCoolPipeSuppAngle     = 22.5;
const Double_t AliITSv11GeometrySDD::fCoolPipeSuppSlitL     = 4.9*fgkmm;
const Double_t AliITSv11GeometrySDD::fCoolPipeSuppAxeDist   = 3.05*fgkmm;
const Double_t AliITSv11GeometrySDD::fCoolPipeInnerDiam     = 1.84*fgkmm;
const Double_t AliITSv11GeometrySDD::fCoolPipeOuterDiam     = 2.*fgkmm;

const Double_t AliITSv11GeometrySDD::fBTBthick          = 0.25 *fgkmm;
const Double_t AliITSv11GeometrySDD::fBTBlength         = 55. *fgkmm;
const Double_t AliITSv11GeometrySDD::fBTBwidth          = 18*fgkmm;
const Double_t AliITSv11GeometrySDD::fBTBaxisAtoBottom  = 4*fgkmm;
const Double_t AliITSv11GeometrySDD::fBTBaxisAtoBase    = 1.2*fgkmm;
const Double_t AliITSv11GeometrySDD::fRadiusAminBTB     = 1. *fgkmm;
const Double_t AliITSv11GeometrySDD::fRadiusBminBTB     = 0.53 *fgkmm;
const Double_t AliITSv11GeometrySDD::fBTBHoleLength     = 15 *fgkmm;
const Double_t AliITSv11GeometrySDD::fBTBHolewidth      =  6 *fgkmm;
const Double_t AliITSv11GeometrySDD::fBTBHoleRefX       = 10 *fgkmm;
const Double_t AliITSv11GeometrySDD::fBTBHoleRefY       = 6.5 *fgkmm;

const Double_t AliITSv11GeometrySDD::fLay3Rmin          = 130.*fgkmm; //not min! Rmin virtual tube
const Double_t AliITSv11GeometrySDD::fLay3Rmax          = 190.*fgkmm; //not min! Rmax virtual tube
const Double_t AliITSv11GeometrySDD::fLay3Length        = (524.+0.)*fgkmm; // ladder+supporting rings (length of the virtual tube)
const Double_t AliITSv11GeometrySDD::fLay3LadderLength  = 524.*fgkmm;
const Double_t AliITSv11GeometrySDD::fLay3DetShortRadius = 146.0*fgkmm;
const Double_t AliITSv11GeometrySDD::fLay3DetLongRadius  = 152.0*fgkmm;
const Double_t AliITSv11GeometrySDD::fLay3LaddTopCornerEnd = 15.6*fgkmm;
const Int_t    AliITSv11GeometrySDD::fLay3Ndet          = 6;
const Int_t    AliITSv11GeometrySDD::fLay3Nladd         = 14;
const Double_t AliITSv11GeometrySDD::fLay3CoolPipeSuppH = 7.5*fgkmm;

const Double_t AliITSv11GeometrySDD::fLay4Rmin          = 220.*fgkmm; //not min! Rmin virtual tube
const Double_t AliITSv11GeometrySDD::fLay4Rmax          = 290.*fgkmm; //not min! Rmax virtual tube
const Double_t AliITSv11GeometrySDD::fLay4Length        = (671.+0.)*fgkmm; // ladder+supporting rings (length of the virtual tube)
const Double_t AliITSv11GeometrySDD::fLay4LadderLength   = 671.*fgkmm;
const Double_t AliITSv11GeometrySDD::fLay4DetShortRadius = 235.0*fgkmm;
const Double_t AliITSv11GeometrySDD::fLay4DetLongRadius  = 240.5*fgkmm;
const Double_t AliITSv11GeometrySDD::fLay4LaddTopCornerEnd = 15.6*fgkmm;
const Int_t    AliITSv11GeometrySDD::fLay4Ndet           = 8;
const Int_t    AliITSv11GeometrySDD::fLay4Nladd          = 22;
const Double_t AliITSv11GeometrySDD::fLay4CoolPipeSuppH  = 10*fgkmm;


ClassImp(AliITSv11GeometrySDD)


//________________________________________________________________________
AliITSv11GeometrySDD::AliITSv11GeometrySDD():AliITSv11Geometry() {
  SetGeomParameters();
  fCoolingOn = kTRUE;
  fAddOnlyLadder3min = -1;
  fAddOnlyLadder3max = -1;
  fAddOnlyLadder4min = -1;
  fAddOnlyLadder4max = -1;
  fAddOnlySegment = -1;
  fColorCarbonFiber = 4;
  fColorRyton = 5;
  fColorPhynox = 7;
  fColorSilicon = 3;
};


//________________________________________________________________________
AliITSv11GeometrySDD::AliITSv11GeometrySDD(Int_t debug)
  : AliITSv11Geometry(debug) {
  SetGeomParameters();
  fCoolingOn = kTRUE;
  fAddOnlyLadder3min = -1;
  fAddOnlyLadder3max = -1;
  fAddOnlyLadder4min = -1;
  fAddOnlyLadder4max = -1;
  fAddOnlySegment = -1;
  fColorCarbonFiber = 4;
  fColorRyton = 5;
  fColorPhynox = 7;
  fColorSilicon = 3;
};


//________________________________________________________________________
void AliITSv11GeometrySDD::SetGeomParameters() {

  fLay3LaddShortRadius  = fLay3DetShortRadius-fLadderBeamRadius+(8-1)*fgkmm; // radius from the center to the CF ladder
  fLay3LaddLongRadius   = fLay3DetLongRadius-fLadderBeamRadius+(8-1)*fgkmm; // radius from the center to the CF ladder
  fLay4LaddShortRadius  = fLay4DetShortRadius-fLadderBeamRadius+(8-1)*fgkmm; // radius from the center to the CF ladder
  fLay4LaddLongRadius   = fLay4DetLongRadius-fLadderBeamRadius+(8-1)*fgkmm;  // radius from the center to the CF ladder


  fLay3sensorZPos[0]= (  35.8+72.4+75.8 )*fgkmm;
  fLay3sensorZPos[1]= (  35.8+72.4      )*fgkmm;
  fLay3sensorZPos[2]= (  35.8           )*fgkmm;
  fLay3sensorZPos[3]= ( -37.9           )*fgkmm;
  fLay3sensorZPos[4]= ( -37.9-74.9      )*fgkmm;
  fLay3sensorZPos[5]= ( -37.9-74.9-71.1 )*fgkmm;

  fLay4sensorZPos[0] = (  38.5+73.2+75.4+71.6 )*fgkmm;
  fLay4sensorZPos[1] = (  38.5+73.2+75.4      )*fgkmm;
  fLay4sensorZPos[2] = (  38.5+73.2           )*fgkmm;
  fLay4sensorZPos[3] = (  38.5                )*fgkmm;
  fLay4sensorZPos[4] = ( -35.6                )*fgkmm;
  fLay4sensorZPos[5] = ( -35.6-74.8           )*fgkmm;
  fLay4sensorZPos[6] = ( -35.6-74.8-72.4      )*fgkmm;
  fLay4sensorZPos[7] = ( -35.6-74.8-72.4-76.  )*fgkmm;

};


//________________________________________________________________________
void AliITSv11GeometrySDD::CheckOverlaps(Double_t precision){

  TGeoVolume *segment = CreateLadderSegment(4,1);
  segment->CheckOverlaps(precision);

  TGeoVolume *ladd3 = CreateLay3Ladder();
  ladd3->CheckOverlaps(precision);

  TGeoVolume *ladd4 = CreateLay4Ladder();
  ladd4->CheckOverlaps(precision);

  TGeoVolume *endLad = CreateEndLadder(3,0);
  endLad->CheckOverlaps(precision);

  TGeoVolume *det3 = CreateLay3Detectors();
  det3->CheckOverlaps(precision);

  TGeoVolume *det4 = CreateLay4Detectors();
  det4->CheckOverlaps(precision);

  TGeoVolume *hyb = CreateHybrid(0);
  hyb->CheckOverlaps(precision);

  TGeoMedium *airSDD = gGeoManager->GetMedium("ITSsddAir");
  TGeoBBox *layBox = new  TGeoBBox("", 300, 300, 300);
  TGeoVolume *lay = new TGeoVolume("lay", layBox, airSDD);
  Layer3(lay);
  Layer4(lay);
  lay->CheckOverlaps(precision);

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
void AliITSv11GeometrySDD::AddTranslationToCombiTrans(TGeoCombiTrans* ct,
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
void AliITSv11GeometrySDD::ShowOnePiece(TGeoVolume *moth) {
// for code developpment and debugging purposes

//     TGeoVolume *vBaseThermalBridge = CreateBaseThermalBridge();
//     Moth->AddNode(vBaseThermalBridge, 1, 0);

//   TGeoVolume* seg = CreateLadderSegment( 4, 0); //lay 4
//   Moth->AddNode(seg, 1, 0);

  AliITSv11GeomSDDcableNap napCable(2.0,0.2);
  Double_t coord1[3] = {0,0,0};
  Double_t coord2[3] = {5,5,5};

  napCable.AddCheckPoint( moth, 0, coord1);
  napCable.AddCheckPoint( moth, 1, coord2);
  napCable.CreateAndInsertCableSegment( 0, 1, moth);

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
    // Placing virtual ladder and detectors volumes following
    // ladder ordering convention
    char rotName[20];
    Int_t iLaddMin = 0;
    Int_t iLaddMax = fLay3Nladd;
    if ((fAddOnlyLadder3min>=0)&&(fAddOnlyLadder3max<fLay3Nladd)) {
        iLaddMin = fAddOnlyLadder3min;
        iLaddMax = fAddOnlyLadder3max+1;
    }
    for (Int_t iLadd = iLaddMin; iLadd < iLaddMax; iLadd++) {
        sprintf(rotName, "ITSsddLay3Ladd%i",iLadd);
        Double_t minRadiusLadBox = fLay3LaddShortRadius;
        if (iLadd%2 != 0) minRadiusLadBox = fLay3LaddLongRadius;
        minRadiusLadBox += ((TGeoBBox*)lay3Ladder->GetShape())->GetDY();
        
        TGeoCombiTrans *ctLadd = CreateCombiTrans(rotName,minRadiusLadBox,
						  0,-90+iLadd*dPhi);
        virtualLayer3->AddNode(lay3Ladder,iLadd,ctLadd);
        sprintf(rotName, "ITSsddLay3DetBox%i",iLadd);
        Double_t minRadiusDetBox = fLay3DetShortRadius;
        if (iLadd%2 != 0) minRadiusDetBox = fLay3DetLongRadius;
        minRadiusDetBox += detBoxThickness/2;
        TGeoCombiTrans *ctDet = CreateCombiTrans(rotName, minRadiusDetBox,
						 0,-90+iLadd*dPhi);
        virtualLayer3->AddNode(lay3Detectors, iLadd, ctDet);
    }
    virtualLayer3->SetVisibility(kFALSE);
    Moth->AddNode(virtualLayer3,1,0);
};


//________________________________________________________________________
TGeoVolume *AliITSv11GeometrySDD::CreateLay3Ladder() {
    // return a box volume containing the CF ladder

    Double_t segmentLength = fSegmentLength;
    TGeoVolume *laddSegmentTemp = CreateLadderSegment(3,0);
    TGeoBBox *ladBox = new TGeoBBox("ITSsddLadBox",
			   ((TGeoBBox*)laddSegmentTemp->GetShape())->GetDX(),
                           ((TGeoBBox*)laddSegmentTemp->GetShape())->GetDY(),
                                    //dX,dY = dX,dY of the segment
                               fLay3LadderLength/2);
    TGeoMedium *airSDD = gGeoManager->GetMedium("ITSsddAir");
    TGeoVolume *virtualLadder = new TGeoVolume("ITSsddLadder",ladBox, airSDD);
    cable.SetInitialNode(virtualLadder);

    // placing virtual ladder segment following detector ordering convention
    //=======================================================================
    char transName[20];
    Int_t iSegmentMin = 0;
    Int_t iSegmentMax = fLay3Ndet;
    if (fAddOnlySegment>=0) {
        iSegmentMin = fAddOnlySegment;
        iSegmentMax = fAddOnlySegment+1;
    }
    for (Int_t iSegment = iSegmentMin; iSegment < iSegmentMax; iSegment++ ) {

        TGeoVolume *laddSegment = CreateLadderSegment(3, iSegment);
        sprintf(transName, "ITSsddLay3LaddSeg%i", iSegment);
        Double_t segmentPos = segmentLength*(fLay3Ndet/2-1-iSegment) 
	                      + segmentLength/2;
        TGeoTranslation *segTr = new TGeoTranslation(transName,0,0,segmentPos);
        virtualLadder->AddNode(laddSegment, iSegment, segTr);
    }

    // putting virtual volume corresponding to the end of ladder
    //=======================================================================
    TGeoVolume *endLadder = CreateEndLadder( 3,0 );
    Double_t endLength = (fLay3LadderLength-fLay3Ndet*fSegmentLength)/2.;
    TGeoTranslation *endTrZPos = new TGeoTranslation("ITSsddEndTrZPos",0,0,
                                   fSegmentLength*(fLay3Ndet/2)+endLength/2.);
    //Euler rotation : about Z, then new X, then new Z
    TGeoRotation *endZNegRot = new TGeoRotation("",90, 180, -90);
    TGeoCombiTrans *endTrZNeg = new TGeoCombiTrans(0,0,
                    -fSegmentLength*(fLay3Ndet/2)-endLength/2., endZNegRot);
    if ((fAddOnlySegment==-1)||(fAddOnlySegment==0))
        virtualLadder->AddNode(endLadder, 1, endTrZPos);
    if ((fAddOnlySegment==-1)||(fAddOnlySegment==fLay3Ndet-1))
        virtualLadder->AddNode(endLadder, 2, endTrZNeg);

    cable.ResetInitialNode();

    virtualLadder->SetVisibility(kFALSE);
    return virtualLadder;
};


//________________________________________________________________________
TGeoArb8 *AliITSv11GeometrySDD::CreateLadderSide(Double_t dz,Double_t angle,
                          Double_t xSign,Double_t L, Double_t H, Double_t l) {
    // Create one half of the V shape corner of CF ladder
  
    TGeoArb8 *cfLaddSide = new TGeoArb8(dz);
    cfLaddSide->SetVertex( 0, 0,  0);
    cfLaddSide->SetVertex( 1, 0, -H);
    cfLaddSide->SetVertex( 2, xSign*(L*TMath::Sin(angle)-l*TMath::Cos(angle)),
			   -L*TMath::Cos(angle)-l*TMath::Sin(angle));
    cfLaddSide->SetVertex( 3, xSign*L*TMath::Sin(angle), -L*TMath::Cos(angle));
    cfLaddSide->SetVertex( 4, 0,  0);
    cfLaddSide->SetVertex( 5, 0, -H);
    cfLaddSide->SetVertex( 6, xSign*(L*TMath::Sin(angle)-l*TMath::Cos(angle)),
			   -L*TMath::Cos(angle)-l*TMath::Sin(angle));
    cfLaddSide->SetVertex(7, xSign*L*TMath::Sin(angle), -L*TMath::Cos(angle));
    return cfLaddSide;
};


//________________________________________________________________________
void AliITSv11GeometrySDD::AddLadderCFstruct(Double_t dy, TGeoVolume* vol) {
    // fill a volume (segment) with the CF structure of a ladder

    TGeoMedium *carbonFiberLadderStruct = 
      gGeoManager->GetMedium("ITSsddCarbonFiber");
    Double_t segmentLength = fSegmentLength;
    Double_t triangleHeight = fLadderHeight - fLadderBeamRadius;
    Double_t halfTheta = TMath::ATan( 0.5*fLadderWidth/triangleHeight );
    Double_t beta = (TMath::Pi()-2.*halfTheta)/4.;
    Double_t alpha = TMath::Pi()*3./4. - halfTheta/2.;

    //--- The 3 V shape corners of the Carbon Fiber Ladder
    //--- the top V
    TGeoArb8 *cfLaddTop1 = CreateLadderSide( segmentLength/2., halfTheta, -1,
					     fLadderLa, fLadderHa, fLadderl);
    TGeoVolume *cfLaddTopVol1 = new TGeoVolume("ITSsddCFladdTopCornerVol1",
				    cfLaddTop1,carbonFiberLadderStruct);
    TGeoArb8 *cfLaddTop2 = CreateLadderSide( segmentLength/2., halfTheta, 1,
					     fLadderLa, fLadderHa, fLadderl);
    TGeoVolume *cfLaddTopVol2 = new TGeoVolume("ITSsddCFladdTopCornerV2",
			            cfLaddTop2,carbonFiberLadderStruct);
    TGeoTranslation *trTop1 = new TGeoTranslation(0, fLadderHeight/2+dy, 0);
    cfLaddTopVol1->SetLineColor(fColorCarbonFiber);
    cfLaddTopVol2->SetLineColor(fColorCarbonFiber);
    vol->AddNode(cfLaddTopVol1, 1, trTop1);
    vol->AddNode(cfLaddTopVol2, 1, trTop1);

    //--- the 2 side V
    TGeoArb8 *cfLaddSide1 = CreateLadderSide( segmentLength/2., beta, -1,
					      fLadderLb, fLadderHb, fLadderl);
    TGeoVolume *cfLaddSideVol1 = new TGeoVolume( "ITSsddCFladdSideCornerV1",
                                     cfLaddSide1,carbonFiberLadderStruct);
    TGeoArb8 *cfLaddSide2 = CreateLadderSide( segmentLength/2., beta, 1,
					      fLadderLb, fLadderHb, fLadderl);
    TGeoVolume *cfLaddSideVol2 = new TGeoVolume( "ITSsddCFladdSideCornerV2",
				     cfLaddSide2,carbonFiberLadderStruct);
    cfLaddSideVol1->SetLineColor(fColorCarbonFiber);
    cfLaddSideVol2->SetLineColor(fColorCarbonFiber);

    Double_t dYTranslation = fLadderHeight/2. - 
                      0.5*fLadderWidth*TMath::Tan(beta) - fLadderBeamRadius;
    // because center of the triangle doesn't correspond to virtual vol. center
    Double_t distCenterSideDown =  0.5*fLadderWidth/TMath::Cos(beta);
    TGeoCombiTrans *ctSideR = CreateCombiTrans("", distCenterSideDown, 0,
                                               alpha*TMath::RadToDeg());
    AddTranslationToCombiTrans(ctSideR, 0, -dYTranslation+dy, 0);
    TGeoCombiTrans *ctSideL = CreateCombiTrans("", distCenterSideDown,0,
                                               -alpha*TMath::RadToDeg());
    AddTranslationToCombiTrans(ctSideL, 0, -dYTranslation+dy, 0);
    vol->AddNode(cfLaddSideVol1, 1, ctSideR);
    vol->AddNode(cfLaddSideVol2, 1, ctSideR);
    vol->AddNode(cfLaddSideVol1, 2, ctSideL);
    vol->AddNode(cfLaddSideVol2, 2, ctSideL);

    //--- The beams
    // Beams on the sides
    Double_t beamPhiPrime = TMath::ASin(1./TMath::Sqrt( (1+TMath::Sin(2*beta)*
                TMath::Sin(2*beta)/(TanD(fBeamSidePhi)*TanD(fBeamSidePhi))) ));
    if(GetDebug(1))
      cout<<"Phi prime = "<<beamPhiPrime*TMath::RadToDeg()<<endl;
    Double_t beamLength = TMath::Sqrt( fLadderHeight*fLadderHeight/
                       ( TMath::Sin(beamPhiPrime)*TMath::Sin(beamPhiPrime))
		       + fLadderWidth*fLadderWidth/4.)-fLadderLa/2-fLadderLb/2;
    TGeoTubeSeg *sideBeam= new TGeoTubeSeg(0, fLadderBeamRadius, beamLength/2.,
					    0, 180);
    TGeoVolume *cfSideBeamVol = new TGeoVolume("ITSsddCFSideBeamVol", sideBeam,
                                               carbonFiberLadderStruct);
    cfSideBeamVol->SetLineColor(fColorCarbonFiber);

    //Euler rotation : about Z, then new X, then new Z
    TGeoRotation *beamRot1 = new TGeoRotation("",90-2.*beta*TMath::RadToDeg(),
				 -beamPhiPrime*TMath::RadToDeg(),-90);
    TGeoCombiTrans *beamTransf1 = new TGeoCombiTrans( 0.5*triangleHeight*
				      TMath::Tan(halfTheta),
				      fLadderBeamRadius/2. + dy,
			       	      -3*segmentLength/8, beamRot1);
    TGeoCombiTrans *beamTransf2 = new TGeoCombiTrans( 0.5*triangleHeight*
				      TMath::Tan(halfTheta),
				      fLadderBeamRadius/2. + dy,
				      segmentLength/8, beamRot1);
    TGeoRotation *beamRot2 = new TGeoRotation("", 90-2.*beta*TMath::RadToDeg(),
				 beamPhiPrime*TMath::RadToDeg(), -90);

    TGeoCombiTrans *beamTransf3 = new TGeoCombiTrans(0.5*triangleHeight*
				      TMath::Tan(halfTheta),
				      fLadderBeamRadius/2. + dy,
				      -segmentLength/8, beamRot2);
    TGeoCombiTrans *beamTransf4 = new TGeoCombiTrans(0.5*triangleHeight*
				      TMath::Tan(halfTheta),
				      fLadderBeamRadius/2. + dy,
				      3*segmentLength/8, beamRot2);
    TGeoRotation *beamRot3 = new TGeoRotation("", 90+2.*beta*TMath::RadToDeg(),
				 beamPhiPrime*TMath::RadToDeg(), -90);

    TGeoCombiTrans *beamTransf5 = new TGeoCombiTrans(-0.5*triangleHeight*
				      TMath::Tan(halfTheta),
			              fLadderBeamRadius/2. + dy,
				      -3*segmentLength/8,beamRot3);
    TGeoCombiTrans *beamTransf6 = new TGeoCombiTrans(-0.5*triangleHeight*
				      TMath::Tan(halfTheta),
			              fLadderBeamRadius/2. + dy,
				      segmentLength/8, beamRot3);
    TGeoRotation *beamRot4 = new TGeoRotation("", 90+2.*beta*TMath::RadToDeg(),
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
    bottomBeam1Vol->SetLineColor(fColorCarbonFiber);

    TGeoRotation *bottomBeamRot1 = new TGeoRotation("",90, 90, 90);
    TGeoCombiTrans *bottomBeamTransf1 = new TGeoCombiTrans(0,
                    -(fLadderHeight/2-fLadderBeamRadius)+dy,0, bottomBeamRot1);
    TGeoTubeSeg *bottomBeam2 = new TGeoTubeSeg(0, fLadderBeamRadius,
			           fLadderWidth/2.-fLadderLb/3, 0, 90);
    TGeoVolume *bottomBeam2Vol = new TGeoVolume("ITSsddBottomBeam2Vol",
                                bottomBeam2, carbonFiberLadderStruct);
    bottomBeam2Vol->SetLineColor(fColorCarbonFiber);
    TGeoCombiTrans *bottomBeamTransf2 = new TGeoCombiTrans(0,
                                     -(fLadderHeight/2-fLadderBeamRadius)+dy,
                                     -segmentLength/2, bottomBeamRot1);
    TGeoRotation *bottomBeamRot2 = new TGeoRotation("",-90, 90, -90);
    TGeoCombiTrans *bottomBeamTransf3 = new TGeoCombiTrans(0,
                                        -(fLadderHeight/2 - fLadderBeamRadius)
					+ dy, segmentLength/2, bottomBeamRot2);

    TGeoTubeSeg *bottomBeam3 = new TGeoTubeSeg(0, fLadderBeamRadius,
			       0.5*fLadderWidth/SinD(fBottomBeamAngle)
			       - fLadderLb/3, 0, 180);
    TGeoVolume *bottomBeam3Vol = new TGeoVolume("ITSsddBottomBeam3Vol",
                                 bottomBeam3, carbonFiberLadderStruct);
    bottomBeam3Vol->SetLineColor(fColorCarbonFiber);
    //bottomBeam3Vol->SetLineColor(2);

    // be careful on the next 2 beams : when "reading" from -z to +z and 
    // from the bottom of the ladder, it should draw a Lambda, and not a V
    TGeoRotation *bottomBeamRot4 = new TGeoRotation("", -90, fBottomBeamAngle,
                                                    -90);
    TGeoCombiTrans *bottomBeamTransf4 = new TGeoCombiTrans(0, 
      -(fLadderHeight/2-fLadderBeamRadius)+dy,-segmentLength/4,bottomBeamRot4);
    TGeoRotation *bottomBeamRot5 = new TGeoRotation("",-90,-fBottomBeamAngle,
                                                    -90);
    TGeoCombiTrans *bottomBeamTransf5 = new TGeoCombiTrans(0,
     -(fLadderHeight/2-fLadderBeamRadius)+dy,segmentLength/4, bottomBeamRot5);

    vol->AddNode(bottomBeam1Vol, 1, bottomBeamTransf1);
    vol->AddNode(bottomBeam2Vol, 1, bottomBeamTransf2);
    vol->AddNode(bottomBeam2Vol, 2, bottomBeamTransf3);
    vol->AddNode(bottomBeam3Vol, 1, bottomBeamTransf4);
    vol->AddNode(bottomBeam3Vol, 2, bottomBeamTransf5);
};


//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreateHybrid(Int_t iSeg) {
// return a box containing the front-end hybrid
// Be careful : electronics side is at y<0, thermal bridge at y>0

  Double_t volumeThick =  0.2;                                  // <===== 0.2 tmp
  TGeoBBox *hybridBox = new TGeoBBox("",fHybridWidth/2, volumeThick/2,
				     (fHybridLength)/2);  // include space on each side for cables ???  
  TGeoMedium *airSDD = gGeoManager->GetMedium("ITSsddAir");
  TGeoVolume *VirtualHybrid = new TGeoVolume( "ITSsddHybridVol", hybridBox,
					      airSDD);
  
  TGeoBBox *sThermalBridge = new TGeoBBox( "", fHybridWidth/2,
					   fHybridThBridgeThick/2,
					   fHybridLength/2);
  TGeoMedium *carbonFiberLadderStruct =
    gGeoManager->GetMedium("ITSsddCarbonFiber");
  
  TGeoVolume *vThermalBridge = new TGeoVolume("ITSsddHybridThBridge",
					      sThermalBridge,
					      carbonFiberLadderStruct);
  vThermalBridge->SetLineColor(fColorCarbonFiber);

  TGeoTranslation *thBridgeTr = new TGeoTranslation("", 0,
				    volumeThick/2-fHybridThBridgeThick/2, 0);
  VirtualHybrid->AddNode(vThermalBridge, 1, thBridgeTr);

  Double_t coord[3];
  coord[0] = 0;coord[1] = 0;coord[2] = 0;
  cable.AddCheckPoint(VirtualHybrid, iSeg, coord);
  VirtualHybrid->SetVisibility(kFALSE);
  cable.GetCheckPoint( cable.GetNCheckPoints()-2, 0, 100, coord);
  return VirtualHybrid;
};


//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreateLadderSegment(Int_t iLay, Int_t iSeg) {
    // Return a box volume containing a segment of a ladder.

    //***
    Double_t tDY = -0.5; //space left on top of the ladder
    Double_t segmentLength = fSegmentLength;

    Double_t coolPipeSuppH = 0;
    Double_t sensorCenterZPos = 0; // z in segment local coord syst.
    if (iLay==3) {
      coolPipeSuppH = fLay3CoolPipeSuppH;
      sensorCenterZPos = fLay3sensorZPos[iSeg]-
	                 (fSegmentLength*fLay3Ndet/2. - 
			  fSegmentLength/2-(iSeg)*fSegmentLength);
    } else {
      coolPipeSuppH = fLay4CoolPipeSuppH;
      sensorCenterZPos = fLay4sensorZPos[iSeg]-
	                 (fSegmentLength*fLay4Ndet/2. -
			  fSegmentLength/2-(iSeg)*fSegmentLength);
    };
    if(GetDebug(1)){
      cout << "Segment ("<< iLay <<',' << iSeg 
	   << ") : sensor z shift in local segment coord.=" 
	   << sensorCenterZPos << endl;
    };
    //***

    TGeoBBox *segBox = new TGeoBBox("ITSsddSegBox",
                                    fLadderWidth/2+fPinSuppWidth+0.5,  // +0.5 is for include volume of hybrids !
                                    fLadderHeight/2+TMath::Abs(tDY),
                                    segmentLength/2);
    TGeoMedium *airSDD = gGeoManager->GetMedium("ITSsddAir");
    TGeoVolume *virtualSeg = new TGeoVolume("ITSsddSegment",
					    segBox, airSDD);

    //**********************************
    // Carbon fiber structure :
    AddLadderCFstruct(tDY, virtualSeg);

    //**********************************
    // Pine support of the sensors :
    TGeoRotation *rotPS1 = new TGeoRotation("",0,-90,90);
    TGeoRotation *rotPS2 = new TGeoRotation("",0,-90,-90);
    TGeoCombiTrans *transPS1 = new TGeoCombiTrans( fPinDYOnSensor,
				   - fLadderHeight/2.-TMath::Abs(tDY)
				   + fPinSuppHeight/2.,
				   sensorCenterZPos + fPinDXminOnSensor,rotPS1);
    TGeoCombiTrans *transPS2 = new TGeoCombiTrans( fPinDYOnSensor,
				   - fLadderHeight/2.-TMath::Abs(tDY)
				   + fPinSuppHeight/2.,
				   sensorCenterZPos + fPinDXminOnSensor + 
				   fPinPinDDXOnSensor, rotPS1);
    TGeoCombiTrans *transPS3 = new TGeoCombiTrans( fPinDYOnSensor,
				   - fLadderHeight/2.-TMath::Abs(tDY)
				   + fPinSuppHeight/2.,
				   sensorCenterZPos-fPinDXminOnSensor,rotPS1);
    TGeoCombiTrans *transPS4 = new TGeoCombiTrans( fPinDYOnSensor,
				   - fLadderHeight/2.-TMath::Abs(tDY)
				   + fPinSuppHeight/2.,
				   sensorCenterZPos-fPinDXminOnSensor -
			           fPinPinDDXOnSensor, rotPS1);
    TGeoCombiTrans *transPS5 = new TGeoCombiTrans( -fPinDYOnSensor,
				   - fLadderHeight/2. - TMath::Abs(tDY)
				   + fPinSuppHeight/2.,
				   sensorCenterZPos+fPinDXminOnSensor,rotPS2);
    TGeoCombiTrans *transPS6 = new TGeoCombiTrans( -fPinDYOnSensor,
				   - fLadderHeight/2. - TMath::Abs(tDY)
				   + fPinSuppHeight/2.,
				   sensorCenterZPos+fPinDXminOnSensor +
                                   fPinPinDDXOnSensor, rotPS2);
    TGeoCombiTrans *transPS7 = new TGeoCombiTrans( -fPinDYOnSensor,
				   - fLadderHeight/2.-TMath::Abs(tDY)
				   + fPinSuppHeight/2.,
				   sensorCenterZPos-fPinDXminOnSensor,rotPS2);
    TGeoCombiTrans *transPS8 = new TGeoCombiTrans( -fPinDYOnSensor,
				   - fLadderHeight/2.-TMath::Abs(tDY)
				   + fPinSuppHeight/2.,
				   sensorCenterZPos - fPinDXminOnSensor
				   - fPinPinDDXOnSensor, rotPS2);
    TGeoVolume *pinSupport = CreatePinSupport();
    virtualSeg->AddNode(pinSupport, 1, transPS1);
    virtualSeg->AddNode(pinSupport, 2, transPS2);
    virtualSeg->AddNode(pinSupport, 3, transPS3);
    virtualSeg->AddNode(pinSupport, 4, transPS4);
    virtualSeg->AddNode(pinSupport, 5, transPS5);
    virtualSeg->AddNode(pinSupport, 6, transPS6);
    virtualSeg->AddNode(pinSupport, 7, transPS7);
    virtualSeg->AddNode(pinSupport, 8, transPS8);

    //**********************************
    // Cooling pipe supports :
    Double_t triangleHeight = fLadderHeight - fLadderBeamRadius;
    Double_t halfTheta = TMath::ATan( 0.5*fLadderWidth/triangleHeight );
    Double_t triangleCPaxeDist = fCoolPipeSuppAxeDist-fCoolPipeSuppWidthExt-
                                 fCoolPipeSuppWidthIn+fLadderBeamRadius;

    Double_t coolPipeSuppL = TMath::Tan(halfTheta)*
                             (triangleHeight+triangleCPaxeDist/
			      TMath::Sin(halfTheta)-coolPipeSuppH);

    TGeoRotation *rotCPS2 = new TGeoRotation("",-halfTheta*TMath::RadToDeg(),
					     -90, 90);
    TGeoRotation *rotCPS1 = new TGeoRotation("",halfTheta*TMath::RadToDeg(),
					     -90, -90);
    TGeoCombiTrans *transCPS1 = new TGeoCombiTrans(coolPipeSuppL,
				    -fLadderHeight/2.-TMath::Abs(tDY)
				    +coolPipeSuppH+fLadderBeamRadius,
				    -segmentLength/2., rotCPS1);
    TGeoCombiTrans *transCPS2 = new TGeoCombiTrans(-coolPipeSuppL,
				    -fLadderHeight/2.-TMath::Abs(tDY)
				    +coolPipeSuppH+fLadderBeamRadius,
				    segmentLength/2., rotCPS2);
    TGeoCombiTrans *transCPS3 = new TGeoCombiTrans(coolPipeSuppL,
				    -fLadderHeight/2.-TMath::Abs(tDY)
				    +coolPipeSuppH+fLadderBeamRadius,
				    segmentLength/2., rotCPS1);
    TGeoCombiTrans *transCPS4 = new TGeoCombiTrans(-coolPipeSuppL,
				    -fLadderHeight/2.-TMath::Abs(tDY)
				    +coolPipeSuppH+fLadderBeamRadius,
				    -segmentLength/2., rotCPS2);

    TGeoVolume *coolPipeSuppLeft = CreateCoolPipeSupportL();
    TGeoVolume *coolPipeSuppRight = CreateCoolPipeSupportR();
    virtualSeg->AddNode(coolPipeSuppLeft,  1, transCPS1);
    virtualSeg->AddNode(coolPipeSuppLeft,  2, transCPS2);
    virtualSeg->AddNode(coolPipeSuppRight, 1, transCPS3);
    virtualSeg->AddNode(coolPipeSuppRight, 2, transCPS4);

    //**********************************
    // Cooling pipes :
    TGeoTube *coolingPipeShape = new TGeoTube( fCoolPipeInnerDiam/2,
					       fCoolPipeOuterDiam/2,
					       segmentLength/2);
    TGeoTube *coolerShape = new TGeoTube( 0, fCoolPipeInnerDiam/2,
					  segmentLength/2);

    //medium = phynox ?                              To code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //medium = water ?                               To code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TGeoMedium *phynoxSDD = gGeoManager->GetMedium("ITSsddStaselite4411w");
    TGeoMedium *coolerMediumSDD = gGeoManager->GetMedium("ITSsddStaselite4411w");
    TGeoVolume *coolingPipe = new TGeoVolume("ITSsddCoolingPipe",
					     coolingPipeShape, phynoxSDD );
    coolingPipe->SetLineColor(fColorPhynox);
    TGeoVolume *cooler = new  TGeoVolume("ITSsddCoolingLiquid",coolerShape,
					 coolerMediumSDD );

    TGeoTranslation *pipeTr1 = new TGeoTranslation(coolPipeSuppL,
				  -fLadderHeight/2.-TMath::Abs(tDY)+
				   fLadderBeamRadius+coolPipeSuppH, 0);
    TGeoTranslation *pipeTr2 = new TGeoTranslation(-coolPipeSuppL,
				  -fLadderHeight/2.-TMath::Abs(tDY)+
				   fLadderBeamRadius+coolPipeSuppH, 0);

    virtualSeg->AddNode(coolingPipe, 1, pipeTr1);
    virtualSeg->AddNode(coolingPipe, 2, pipeTr2);
    if (fCoolingOn) {
      virtualSeg->AddNode(cooler, 1, pipeTr1);
      virtualSeg->AddNode(cooler, 2, pipeTr2);
    };

    //**********************************
    // Bases of hybrid thermal bridges
    Double_t hybridAngle = 46;                                           // tmp !!!
    Double_t shiftHyb = 0.9-0.2; // shift in comparison with center of thermal Br. base // tmp !!! not clear on 0752/14-A
    // sur 0752/14-A c'est environ 0.9
    // en fait comme la hauteur des cooling pipes depends de la couche, ces deux variables varient probablement
    // suivant la "layer" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    TGeoVolume *baseThermalBridge = CreateBaseThermalBridge();
    TGeoRotation rotHybrid1("",0, 0,  -90 - hybridAngle);
    TGeoRotation rotHybrid2("",0 ,180, 90 - hybridAngle);
    TGeoCombiTrans *baseTr1 = new TGeoCombiTrans(*pipeTr2, rotHybrid1);
    TGeoCombiTrans *baseTr2 = new TGeoCombiTrans(*pipeTr1, rotHybrid2);
    virtualSeg->AddNode(baseThermalBridge, 1, baseTr1);
    virtualSeg->AddNode(baseThermalBridge, 2, baseTr2);

    //**********************************
    // the 2 hybrids :
    TGeoVolume *Hybrid = CreateHybrid(iSeg);

    Double_t hybDy = ((TGeoBBox*)Hybrid->GetShape())->GetDY();
    Double_t axeToHybridCenterDist = fBTBaxisAtoBase+hybDy;

    Double_t hybrVolDX = ( axeToHybridCenterDist*CosD(hybridAngle) 
			   - shiftHyb*SinD(hybridAngle) );
    Double_t hybrVolDY = ( axeToHybridCenterDist*SinD(hybridAngle)
			   + shiftHyb*CosD(hybridAngle) );

    TGeoCombiTrans *hybTr1 = new TGeoCombiTrans(*baseTr1);
    AddTranslationToCombiTrans( hybTr1, -hybrVolDX, hybrVolDY, 0);
    TGeoCombiTrans *hybTr2 = new TGeoCombiTrans(*baseTr2);
    AddTranslationToCombiTrans( hybTr2, hybrVolDX, hybrVolDY, 0);

    virtualSeg->AddNode(Hybrid, 1, hybTr1);
    virtualSeg->AddNode(Hybrid, 2, hybTr2);


    //**********************************
    virtualSeg->SetVisibility(kFALSE);
    return virtualSeg;
};


//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreatePinSupport() {
//
// Create a pine support
// axis of rotation is the cone axis, center in its middle
//
    TGeoCone *cone = new TGeoCone("ITSsddPinSuppCone",fPinSuppHeight/2.,
                                  0,fPinSuppRmax,0,fPinSuppRmax-
                                  fPinSuppHeight*TanD(fPinSuppConeAngle) );
    TGeoBBox *tong = new TGeoBBox("ITSsddPinSuppTong",fPinSuppRmax,
                                  fPinSuppLength/2.,fPinSuppThickness/2.);
    TGeoTube *hole = new TGeoTube("ITSsddPinSuppHole",0,fPinR,
                                  fPinSuppHeight/2.);
    if(GetDebug(3)){// Remove compiler warning.
        cone->InspectShape();
        tong->InspectShape();
        hole->InspectShape();
    };

    TGeoTranslation *tongTrans = new TGeoTranslation("ITSsddPinSuppTongTr",0,
                   fPinSuppLength/2.,-fPinSuppHeight/2.+fPinSuppThickness/2.);
    tongTrans->RegisterYourself();
    TGeoCompositeShape *pinSupportShape = new TGeoCompositeShape(
               "ITSssdPinSupportShape","(ITSsddPinSuppCone+"
               "ITSsddPinSuppTong:ITSsddPinSuppTongTr)-ITSsddPinSuppHole");

    //medium = ryton ?                           To code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TGeoMedium *rytonSDD = gGeoManager->GetMedium("ITSsddStaselite4411w");
    TGeoVolume *pinSupport = new TGeoVolume("ITSssdPinSupport",pinSupportShape,
                                            rytonSDD);
    pinSupport->SetLineColor(fColorRyton);
    return pinSupport;
    // include the pin itself                           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
};

//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreateCoolPipeSupportL() {
//
// Create half of the cooling pipe support (ALR-0752/3)
//

  Double_t diffX = fCoolPipeSuppHeight*TanD(fCoolPipeSuppAngle);
  
  TGeoArb8 *side1 = new TGeoArb8(fCoolPipeSuppHeight/2.);
  side1->SetVertex( 0, 0, -fCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 1, fCoolPipeSuppMaxLength/2.-diffX,
		       -fCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 2, fCoolPipeSuppMaxLength/2.-diffX,
		       fCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 3, 0,  fCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 4, 0, -fCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 5, fCoolPipeSuppMaxLength/2.,
		       -fCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 6, fCoolPipeSuppMaxLength/2.,
		       fCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 7, 0,  fCoolPipeSuppWidthExt/2.);
  side1->SetName("ITSsddCPSside1");

  TGeoTranslation *side1Tr = new TGeoTranslation("ITSsddCPStr1",0,
				 - fCoolPipeSuppAxeDist
				 + fCoolPipeSuppWidthExt/2., 0);
  side1Tr->RegisterYourself();
  TGeoTranslation *side2Tr = new TGeoTranslation("ITSsddCPStr2",0,
				 - fCoolPipeSuppAxeDist
                                 + fCoolPipeSuppWidthExt*3/2.
		      		 + fCoolPipeSuppWidthIn,0);
  side2Tr->RegisterYourself();
  
  TGeoBBox *middle = new TGeoBBox("ITSsddCPSmiddle",
			 (fCoolPipeSuppMaxLength/2.-fCoolPipeSuppSlitL)/2.,
			 fCoolPipeSuppWidthIn/2., fCoolPipeSuppHeight/2.);
  TGeoTranslation *middleTr = 
    new TGeoTranslation("ITSsddCPStr3",
			(fCoolPipeSuppMaxLength/2.-fCoolPipeSuppSlitL)/2.,
			-fCoolPipeSuppAxeDist+fCoolPipeSuppWidthExt
			+fCoolPipeSuppWidthIn/2., 0);
  middleTr->RegisterYourself();
  
  TGeoBBox *axeBox = new TGeoBBox("ITSsddCPSaxeBox",
				  fCoolPipeSuppTongW/4.,
				  (fCoolPipeSuppFulWidth
				   - 2*fCoolPipeSuppWidthExt
				   - fCoolPipeSuppWidthIn)/2,
				  fCoolPipeSuppHeight/2.);
  
  TGeoTranslation *axeBoxTr = new TGeoTranslation("ITSsddCPSAxBoxTr",
				  fCoolPipeSuppTongW/4.,
			       	  - fCoolPipeSuppAxeDist
				  + fCoolPipeSuppFulWidth
				  - axeBox->GetDY(), 0);
  axeBoxTr->RegisterYourself();

  TGeoTube *axe = new TGeoTube("ITSsddCPSaxe",0,fCoolPipeSuppHoleDiam/2.,
			       fCoolPipeSuppTongW/4.);
  TGeoRotation *axeRot = new TGeoRotation("ITSsddCPSaxeRot",90,90,0);

  TGeoCombiTrans *axeTrans = new TGeoCombiTrans("ITSsddCPSaxeTr",
				 fCoolPipeSuppTongW/4.,0,0,axeRot);
  axeTrans->RegisterYourself();

  if(GetDebug(3)){
    middle->InspectShape();
    axe->InspectShape();
  };

  //medium = ryton ? To code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TGeoMedium *rytonSDD = gGeoManager->GetMedium("ITSsddStaselite4411w");
  
  TGeoCompositeShape *coolPipeSuppShape = new TGeoCompositeShape(
				        "ITSsddCoolPipeSuppShapeL",
				        "ITSsddCPSmiddle:ITSsddCPStr3"
				        "+ITSsddCPSside1:ITSsddCPStr1"
				        "+ITSsddCPSside1:ITSsddCPStr2"
				        "+ITSsddCPSaxeBox:ITSsddCPSAxBoxTr"
			                "-ITSsddCPSaxe:ITSsddCPSaxeTr");
  TGeoVolume *coolPipeSupp = new  TGeoVolume("ITSsddCoolPipeSupportL",
					     coolPipeSuppShape, rytonSDD);

  coolPipeSupp->SetLineColor(fColorRyton);
  return coolPipeSupp;
};

//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreateCoolPipeSupportR() {
//
//Create half of the cooling pipe support (ALR-0752/3)
//

  Double_t diffX = fCoolPipeSuppHeight*TanD(fCoolPipeSuppAngle);
  
  TGeoArb8 *side1 = new TGeoArb8(fCoolPipeSuppHeight/2.);
  side1->SetVertex( 0, 0, -fCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 1, -(fCoolPipeSuppMaxLength/2.-diffX),
		       -fCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 2, -(fCoolPipeSuppMaxLength/2.-diffX),
		       fCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 3, 0,  fCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 4, 0, -fCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 5, -fCoolPipeSuppMaxLength/2.,
		       -fCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 6, -fCoolPipeSuppMaxLength/2.,
		       fCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 7, 0,  fCoolPipeSuppWidthExt/2.);
  side1->SetName("ITSsddCPSside1R");

  TGeoTranslation *side1Tr = new TGeoTranslation("ITSsddCPStr1R",0,
				 - fCoolPipeSuppAxeDist
				 + fCoolPipeSuppWidthExt/2., 0);
  side1Tr->RegisterYourself();
  TGeoTranslation *side2Tr = new TGeoTranslation("ITSsddCPStr2R",0,
				 - fCoolPipeSuppAxeDist
				 + fCoolPipeSuppWidthExt*3/2.
		      		 + fCoolPipeSuppWidthIn, 0);
  side2Tr->RegisterYourself();
  
  TGeoBBox *middle = new TGeoBBox("ITSsddCPSmiddleR",
				  (fCoolPipeSuppMaxLength/2.
				   - fCoolPipeSuppSlitL)/2.,
				  fCoolPipeSuppWidthIn/2., 
				  fCoolPipeSuppHeight/2.);
  TGeoTranslation *middleTr = 
    new TGeoTranslation("ITSsddCPStr3R",
			-( fCoolPipeSuppMaxLength/2.
			   -fCoolPipeSuppSlitL)/2.,
			-fCoolPipeSuppAxeDist + fCoolPipeSuppWidthExt
			+ fCoolPipeSuppWidthIn/2.,0);
  middleTr->RegisterYourself();
  
  TGeoBBox *axeBox = new TGeoBBox("ITSsddCPSaxeBoxR",
				  fCoolPipeSuppTongW/4.,
				  (fCoolPipeSuppFulWidth
				   - 2*fCoolPipeSuppWidthExt
				   - fCoolPipeSuppWidthIn)/2,
				  fCoolPipeSuppHeight/2.);
  
  TGeoTranslation *axeBoxTr = new TGeoTranslation("ITSsddCPSAxBoxTrR",
				  - fCoolPipeSuppTongW/4.,
			       	  - fCoolPipeSuppAxeDist
				  + fCoolPipeSuppFulWidth
				  - axeBox->GetDY(),0);
  axeBoxTr->RegisterYourself();

  TGeoTube *axe = new TGeoTube("ITSsddCPSaxeR",0,fCoolPipeSuppHoleDiam/2.,
			       fCoolPipeSuppTongW/4.);
  TGeoRotation *axeRot = new TGeoRotation("ITSsddCPSaxeRotR",90,90,0);

  TGeoCombiTrans *axeTrans = new TGeoCombiTrans("ITSsddCPSaxeTrR",
				 -fCoolPipeSuppTongW/4.,0,0,axeRot);
  axeTrans->RegisterYourself();

  if(GetDebug(3)){
    middle->InspectShape();
    axe->InspectShape();
  };
  
  TGeoCompositeShape *coolPipeSuppShape = new TGeoCompositeShape(
			              "ITSsddCoolPipeSuppShapeR",
			              "ITSsddCPSmiddleR:ITSsddCPStr3R"
			              "+ITSsddCPSside1R:ITSsddCPStr1R"
			              "+ITSsddCPSside1R:ITSsddCPStr2R"
				      "+ITSsddCPSaxeBoxR:ITSsddCPSAxBoxTrR"
				      "-ITSsddCPSaxeR:ITSsddCPSaxeTrR");

  //medium = ryton ? To code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TGeoMedium *rytonSDD = gGeoManager->GetMedium("ITSsddStaselite4411w");
  TGeoVolume *coolPipeSupp = new TGeoVolume( "ITSsddCoolPipeSupportR",
					     coolPipeSuppShape, rytonSDD);
  coolPipeSupp->SetLineColor(fColorRyton);

  return coolPipeSupp;
};

//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreateBaseThermalBridge() {
// ALR 0752/8

  Double_t dy = fBTBaxisAtoBase - fRadiusBminBTB - fBTBthick;

  Double_t base1width = fBTBwidth - fBTBaxisAtoBottom - fRadiusBminBTB
                        - (fRadiusAminBTB+fBTBthick);
  TGeoBBox *base1 = new TGeoBBox( "ITSsddBTBbase1", base1width/2.,
				  fBTBthick/2., fBTBlength/2.);
  TGeoTranslation *base1Tr = new TGeoTranslation("ITSsddBTBtr1",
				 fBTBaxisAtoBottom-fBTBwidth+base1width/2.,
				 -(fBTBaxisAtoBase-fBTBthick/2.), 0);
  base1Tr->RegisterYourself();

  Double_t base2width = fBTBaxisAtoBottom - fRadiusAminBTB - fBTBthick
                        - fRadiusBminBTB;
  TGeoBBox *base2 = new TGeoBBox( "ITSsddBTBbase2", base2width/2.,
				  fBTBthick/2., fBTBlength/2.);
  TGeoTranslation *base2Tr = new TGeoTranslation("ITSsddBTBtr2",
				 fBTBaxisAtoBottom - base2width/2.,
				 -(fBTBaxisAtoBase-fBTBthick/2.), 0);
  base2Tr->RegisterYourself();

  TGeoBBox *side = new TGeoBBox( "ITSsddBTBside",
				 fBTBthick/2., dy/2., fBTBlength/2.);
  TGeoTranslation *sideTr1 = new TGeoTranslation("ITSsddBTBsideTr1",
				 -fRadiusAminBTB-fBTBthick/2., -dy/2., 0);
  TGeoTranslation *sideTr2 = new TGeoTranslation("ITSsddBTBsideTr2",
				 fRadiusAminBTB+fBTBthick/2., -dy/2., 0);
  sideTr1->RegisterYourself();
  sideTr2->RegisterYourself();

  TGeoBBox *hole = new TGeoBBox( "ITSsddBTBhole", fBTBHolewidth/2.,
				 fBTBthick/2., fBTBHoleLength/2.);
  TGeoTranslation *holeTr1 = new TGeoTranslation("ITSsddBTBholeTr1",
				 - fBTBHoleRefX + fBTBHolewidth/2.,
				 - (fBTBaxisAtoBase-fBTBthick/2.),
				 fBTBHoleRefY+(fBTBHoleLength-fBTBlength)/2.);
  TGeoTranslation *holeTr2 = new TGeoTranslation("ITSsddBTBholeTr2",
				 - fBTBHoleRefX + fBTBHolewidth/2.,
				 - (fBTBaxisAtoBase-fBTBthick/2.),
				 - fBTBHoleRefY-(fBTBHoleLength-fBTBlength)/2.);
  holeTr1->RegisterYourself();
  holeTr2->RegisterYourself();

  Double_t radiusAmaxBTB = fRadiusAminBTB + fBTBthick;
  TGeoTubeSeg *mainAxis = new TGeoTubeSeg( "ITSsddBTBmainAxis",
					   fRadiusAminBTB, radiusAmaxBTB,
					   fBTBlength/2., 0., 180.);
  TGeoTubeSeg *round1 = new TGeoTubeSeg( "ITSsddBTBround1",
			   fRadiusBminBTB, fRadiusBminBTB+fBTBthick,
			   fBTBlength/2., 270., 360.);
  TGeoTranslation *roundTr1 = new TGeoTranslation("ITSsddBTBround1Tr",
				  -(fRadiusAminBTB+fBTBthick+fRadiusBminBTB),
				  -dy, 0);
  roundTr1->RegisterYourself();

  TGeoTubeSeg *round2 = new TGeoTubeSeg( "ITSsddBTBround2",
			   fRadiusBminBTB, fRadiusBminBTB+fBTBthick,
			   fBTBlength/2., 180., 270.);
  TGeoTranslation *roundTr2 = new TGeoTranslation("ITSsddBTBround2Tr",
				  (fRadiusAminBTB+fBTBthick+fRadiusBminBTB),
				  -dy, 0);
  roundTr2->RegisterYourself();

  TGeoCompositeShape *sBaseThermalBridge = new TGeoCompositeShape(
				      "ITSsddBaseThermalBridgeShape",
				      "ITSsddBTBbase1:ITSsddBTBtr1"
				      "+ ITSsddBTBbase2:ITSsddBTBtr2"
				      "+ ITSsddBTBround1:ITSsddBTBround1Tr"
				      "+ ITSsddBTBround2:ITSsddBTBround2Tr"
				      "+ ITSsddBTBside:ITSsddBTBsideTr1"
				      "+ ITSsddBTBside:ITSsddBTBsideTr2"
				      "- ITSsddBTBhole:ITSsddBTBholeTr1"
				      "- ITSsddBTBhole:ITSsddBTBholeTr2"
				      "+ ITSsddBTBmainAxis");

    if(GetDebug(3)){// Remove compiler warning.
        base1->InspectShape();
        base2->InspectShape();
        side->InspectShape();
        hole->InspectShape();
        mainAxis->InspectShape();
        round1->InspectShape();
        round2->InspectShape();
    };

  TGeoMedium *carbonFiberLadderStruct = gGeoManager->GetMedium("ITSsddCarbonFiber");
  TGeoVolume *vBaseThermalBridge = new TGeoVolume( "ITSsddBaseThermalBridge",
						   sBaseThermalBridge,
						   carbonFiberLadderStruct);

  vBaseThermalBridge->SetLineColor(fColorCarbonFiber);
  return vBaseThermalBridge;
};




//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreateEndLadder(Int_t iLay, Int_t) {
    // Return a box volume containing a end of a CF ladder.

    Double_t tDY = -0.5; //space left on top of the ladder

    TGeoMedium *airSDD = gGeoManager->GetMedium("ITSsddAir");
    TGeoMedium *carbonFiberLadderStruct = gGeoManager->GetMedium(
                                          "ITSsddCarbonFiber");
    Double_t length = 0;
    Double_t coolPipeSuppH = 0;

     if (iLay==3) {
       length = (fLay3LadderLength-fLay3Ndet*fSegmentLength)/2.;
       coolPipeSuppH = fLay3CoolPipeSuppH;
     } else {
       length = (fLay4LadderLength-fLay4Ndet*fSegmentLength)/2.;
       coolPipeSuppH = fLay4CoolPipeSuppH;
     };
    
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

    //--- The 3 V shape corners of the Carbon Fiber Ladder
    //--- the top V
    TGeoArb8 *cfLaddTop1 = CreateLadderSide( topCornerLength/2., halfTheta, -1,
                                             fLadderLa, fLadderHa, fLadderl);
    TGeoVolume *cfLaddTopVol1 = new TGeoVolume("ITSsddCFladdTopCornerVol1",
                                          cfLaddTop1,carbonFiberLadderStruct);
    cfLaddTopVol1->SetLineColor(fColorCarbonFiber);
    TGeoArb8 *cfLaddTop2 = CreateLadderSide( topCornerLength/2., halfTheta, 1,
                                             fLadderLa, fLadderHa, fLadderl);
    TGeoVolume *cfLaddTopVol2 = new TGeoVolume("ITSsddCFladdTopCornerV2",
                                        cfLaddTop2,carbonFiberLadderStruct);
    cfLaddTopVol2->SetLineColor(fColorCarbonFiber);
    TGeoTranslation *trTop1 = new TGeoTranslation(0, fLadderHeight/2+tDY,
                                               -(length-topCornerLength)/2.);
    virtualEnd->AddNode(cfLaddTopVol1, 1, trTop1);
    virtualEnd->AddNode(cfLaddTopVol2, 1, trTop1);

    //--- the 2 side V
    TGeoArb8 *cfLaddSide1 = CreateLadderSide( length/2., beta, -1,
                                              fLadderLb, fLadderHb, fLadderl);
    TGeoVolume *cfLaddSideVol1 = new TGeoVolume("ITSsddCFladdSideCornerV1",
                                          cfLaddSide1,carbonFiberLadderStruct);
    cfLaddSideVol1->SetLineColor(fColorCarbonFiber);
    TGeoArb8 *cfLaddSide2 = CreateLadderSide( length/2., beta, 1,
                                              fLadderLb, fLadderHb, fLadderl);
    TGeoVolume *cfLaddSideVol2 = new TGeoVolume("ITSsddCFladdSideCornerV2",
                                          cfLaddSide2,carbonFiberLadderStruct);
    cfLaddSideVol2->SetLineColor(fColorCarbonFiber);
    Double_t dYTranslation = ( fLadderHeight/2. - 0.5*fLadderWidth*
                               TMath::Tan(beta) - fLadderBeamRadius );

    // because center of the triangle doesn't correspond to virtual vol. center
    Double_t distCenterSideDown =  0.5*fLadderWidth/TMath::Cos(beta);
    TGeoCombiTrans *ctSideR = CreateCombiTrans("", distCenterSideDown, 0,
                                               alpha*TMath::RadToDeg());
    AddTranslationToCombiTrans(ctSideR, 0, -dYTranslation+tDY, 0);
    TGeoCombiTrans *ctSideL = CreateCombiTrans("", distCenterSideDown, 0, 
                                               -alpha*TMath::RadToDeg());
    AddTranslationToCombiTrans(ctSideL, 0, -dYTranslation+tDY, 0);
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
     TGeoTubeSeg *sideBeam = new TGeoTubeSeg( 0, fLadderBeamRadius,
					      beamLength/2., 0, 180);
     TGeoVolume *cfSideBeamVol = new TGeoVolume("ITSsddCFSideBeamVol",
                                      sideBeam, carbonFiberLadderStruct);
     cfSideBeamVol->SetLineColor(fColorCarbonFiber);

     //Euler rotation : about Z, then new X, then new Z
     TGeoRotation *beamRot1 = new TGeoRotation("",90-2.*beta*TMath::RadToDeg(),
				  -beamPhiPrime*TMath::RadToDeg(), -90);
     TGeoCombiTrans *beamTransf1 = new TGeoCombiTrans(0.5*triangleHeight*
				       TMath::Tan(halfTheta),
				       fLadderBeamRadius/2. + tDY,
				       -length/2 + segmentLength/8, beamRot1);
     TGeoRotation *beamRot2 = new TGeoRotation("", 90-2.*beta*TMath::RadToDeg(),
				  beamPhiPrime*TMath::RadToDeg(), -90);
     TGeoCombiTrans *beamTransf3 = new TGeoCombiTrans( 0.5*triangleHeight*
				       TMath::Tan(halfTheta),
				       fLadderBeamRadius/2.+tDY,
				       -length/2 + 3*segmentLength/8, beamRot2);
     TGeoRotation *beamRot3 = new TGeoRotation("",90+2.*beta*TMath::RadToDeg(),
				  beamPhiPrime*TMath::RadToDeg(), -90);
     TGeoCombiTrans *beamTransf5 = new TGeoCombiTrans(-0.5*triangleHeight*
				       TMath::Tan(halfTheta),
				       fLadderBeamRadius/2.+tDY,
			       	       -length/2 + segmentLength/8, beamRot3);
     TGeoRotation *beamRot4 = new TGeoRotation("",90 + 2.*beta*TMath::RadToDeg(),
				  -beamPhiPrime*TMath::RadToDeg(), -90);
     TGeoCombiTrans *beamTransf7 = new TGeoCombiTrans(-0.5*triangleHeight*
				       TMath::Tan(halfTheta),
				       fLadderBeamRadius/2. + tDY,
				       -length/2+3*segmentLength/8, beamRot4);
     virtualEnd->AddNode(cfSideBeamVol, 1, beamTransf1);
     virtualEnd->AddNode(cfSideBeamVol, 2, beamTransf3);
     virtualEnd->AddNode(cfSideBeamVol, 3, beamTransf5);
     virtualEnd->AddNode(cfSideBeamVol, 4, beamTransf7);

     //--- Beams of the bottom
     TGeoTubeSeg *bottomBeam1 = new TGeoTubeSeg(0, fLadderBeamRadius,
				    fLadderWidth/2.-fLadderLb/3, 0, 180);
     TGeoVolume *bottomBeam1Vol = new TGeoVolume("ITSsddBottomBeam1Vol",
                                      bottomBeam1, carbonFiberLadderStruct);
     bottomBeam1Vol->SetLineColor(fColorCarbonFiber);
     TGeoRotation *bottomBeamRot1 = new TGeoRotation("",90, 90, 90);
     TGeoCombiTrans *bottomBeamTransf1 = new TGeoCombiTrans(0,
                            -(fLadderHeight/2-fLadderBeamRadius)+tDY,
			    -length/2+fSegmentLength/2, bottomBeamRot1);
     virtualEnd->AddNode(bottomBeam1Vol, 1, bottomBeamTransf1);
     TGeoTubeSeg *bottomBeam2 = new TGeoTubeSeg(0, fLadderBeamRadius,
				    fLadderWidth/2.-fLadderLb/3, 0, 90);
     TGeoVolume *bottomBeam2Vol = new TGeoVolume("ITSsddBottomBeam2Vol",
				      bottomBeam2, carbonFiberLadderStruct);
     bottomBeam2Vol->SetLineColor(fColorCarbonFiber);
     TGeoCombiTrans *bottomBeamTransf2 = new TGeoCombiTrans(0,
           -(fLadderHeight/2-fLadderBeamRadius)+tDY,-length/2, bottomBeamRot1);
     virtualEnd->AddNode(bottomBeam2Vol, 1, bottomBeamTransf2);

     //**********************************
     //the cooling pipe supports
     Double_t triangleCPaxeDist = fCoolPipeSuppAxeDist-fCoolPipeSuppWidthExt-
                                  fCoolPipeSuppWidthIn+fLadderBeamRadius;

     Double_t coolPipeSuppL = TMath::Tan(halfTheta)*
                              (triangleHeight+triangleCPaxeDist/
			       TMath::Sin(halfTheta) - coolPipeSuppH);
     TGeoRotation *rotCPS2 = new TGeoRotation("",-halfTheta*TMath::RadToDeg(),-90,90);
     TGeoRotation *rotCPS1 = new TGeoRotation("",halfTheta*TMath::RadToDeg(),-90,-90);
     TGeoCombiTrans *transCPS1 = new TGeoCombiTrans(coolPipeSuppL,
				     -fLadderHeight/2.-TMath::Abs(tDY)+
				     coolPipeSuppH+fLadderBeamRadius,
				     -length/2., rotCPS1);
     TGeoCombiTrans *transCPS4 = new TGeoCombiTrans(-coolPipeSuppL,
				     -fLadderHeight/2.-TMath::Abs(tDY)+
				     coolPipeSuppH+fLadderBeamRadius,
				     -length/2., rotCPS2);

     TGeoVolume *coolPipeSuppLeft = CreateCoolPipeSupportL();
     virtualEnd->AddNode(coolPipeSuppLeft, 1, transCPS1);
     TGeoVolume *coolPipeSuppRight = CreateCoolPipeSupportR();
     virtualEnd->AddNode(coolPipeSuppRight, 1, transCPS4);
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
    TGeoVolume *virtualSensor = new TGeoVolume("ITSsddSensor",sensorBox,airSDD);

    // the virtual volume is already following the cood convention
    // for the local wafer coord. system : no need for additional rotation matrix
    TGeoBBox *waferShape = new TGeoBBox("ITSsddWaferShape",
				  fWaferWidth/2, fWaferThickness/2,fWaferLength/2);

    //medium = silicon ?                              To code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TGeoMedium *siliconSDD = gGeoManager->GetMedium("ITSsddStaselite4411w");
    TGeoVolume *sensor = new TGeoVolume("ITSsddWafer", waferShape, siliconSDD);
    sensor->SetLineColor(fColorSilicon);
    virtualSensor->AddNode(sensor, 1, 0);

    virtualSensor->SetVisibility(kFALSE);
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
 
    Int_t iSegmentMin = 0;
    Int_t iSegmentMax = fLay3Ndet;
    if (fAddOnlySegment>=0) {
        iSegmentMin = fAddOnlySegment;
        iSegmentMax = fAddOnlySegment+1;
    };

    for (Int_t i=iSegmentMin; i<iSegmentMax; i++) {
        Double_t localZ = fLay3sensorZPos[i];
        Double_t localY = fLadWaferSep/2+fWaferThickness/2;
        if (i%2!=0)
            localY = -localY;
        sprintf(transName, "ITSsddLay3SensorPos%i", i);
        TGeoTranslation *sensorPos = new TGeoTranslation(transName,0,localY,
                                                         localZ);
        virtualDet->AddNode(vSDD, i, sensorPos);
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
    if ((fAddOnlyLadder4min >= 0)&&(fAddOnlyLadder4max < fLay4Nladd)) {
        iLaddMin = fAddOnlyLadder4min;
        iLaddMax = fAddOnlyLadder4max+1;
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

    TGeoVolume *laddSegmentTemp = CreateLadderSegment(4,0);
    TGeoBBox *ladBox = new TGeoBBox("ITSsddLadBox",
                           ((TGeoBBox*)laddSegmentTemp->GetShape())->GetDX(),
			   ((TGeoBBox*)laddSegmentTemp->GetShape())->GetDY(),
			   //dX,dY = dX,dY of the segment
			   fLay4LadderLength/2);  
    TGeoMedium *airSDD = gGeoManager->GetMedium("ITSsddAir");
    TGeoVolume *virtualLadder = new TGeoVolume("ITSsddLadder",ladBox, airSDD);
    Double_t segmentLength = fSegmentLength;
    char transName[20];
    // placing virtual ladder segment following detector ordering convention
    //=======================================================================
    Int_t iSegmentMin = 0;
    Int_t iSegmentMax = fLay4Ndet;
    if (fAddOnlySegment>=0) {
        iSegmentMin = fAddOnlySegment;
        iSegmentMax = fAddOnlySegment+1;
    }
    for (Int_t iSegment= iSegmentMin; iSegment < iSegmentMax; iSegment++ ) {

        TGeoVolume *laddSegment = CreateLadderSegment(4, iSegment);
        sprintf(transName, "ITSsddLay4LaddSeg%i", iSegment);
        Double_t segmentPos = segmentLength*(fLay4Ndet/2-1-iSegment) +
	                      segmentLength/2;
        TGeoTranslation *segTr = new TGeoTranslation(transName,0,0,segmentPos);
        virtualLadder->AddNode(laddSegment, iSegment, segTr);
    }
    // putting virtual volume corresponding to the end of ladder
    //=======================================================================
    TGeoVolume *endLadder = CreateEndLadder( 4,-1 );
    Double_t endLength = (fLay4LadderLength-fLay4Ndet*fSegmentLength)/2.;
    TGeoTranslation *endTrZPos =
        new TGeoTranslation("ITSsddEndTrZPos",0,0,
                            fSegmentLength*(fLay4Ndet/2)+endLength/2.);
    //Euler rotation : about Z, then new X, then new Z
    TGeoRotation *endZNegRot = new TGeoRotation("",90, 180, -90);
    TGeoCombiTrans *endTrZNeg =
        new TGeoCombiTrans(0,0,-fSegmentLength*(fLay4Ndet/2)-endLength/2.,
                           endZNegRot);
    if ((fAddOnlySegment==-1)||(fAddOnlySegment==0))
        virtualLadder->AddNode(endLadder, 1, endTrZPos);
    if ((fAddOnlySegment==-1)||(fAddOnlySegment==fLay4Ndet-1))
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
    Int_t iSegmentMin = 0;
    Int_t iSegmentMax = fLay4Ndet;
    if (fAddOnlySegment>=0) {
        iSegmentMin = fAddOnlySegment;
        iSegmentMax = fAddOnlySegment+1;
    }

    for (Int_t i=iSegmentMin; i<iSegmentMax; i++) {

        Double_t localZ = fLay4sensorZPos[i];
        Double_t localY = fLadWaferSep/2+fWaferThickness/2;
        if (i%2==0)
            localY = -localY;
        sprintf(transName, "ITSsddLay4SensorPos%i", i);
        TGeoTranslation *sensorPos = new TGeoTranslation(transName, 0,
							 localY, localZ);
        virtualDet->AddNode(vSDD, i, sensorPos);
    }
    virtualDet->SetVisibility(kFALSE);
    return virtualDet;
};
