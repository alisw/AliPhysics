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



// General Root includes
//#include <Riostream.h>
#include <TMath.h>
#include <TVectorD.h>

// Root Geometry includes
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoTube.h>
#include <TGeoTorus.h>
#include <TGeoMatrix.h>

#include "AliITSv11GeomCableRound.h"

//*************************************************************************
//   Class for round cables
//
// Ludovic Gaudichet                                   gaudichet@to.infn.it
//*************************************************************************

/*
// ************************************************************************
// Here is a example on how to use this class
// ************************************************************************

  // Getting some media 
  TGeoMedium *air   = gGeoManager->GetMedium("ITS_ITSair");
  TGeoMedium *water = gGeoManager->GetMedium("ITS_WATER");
  TGeoMedium *alu   = gGeoManager->GetMedium("ITS_ITSal"); 

  // Creating a small box inside a bigger one (containers)
  TGeoBBox *box1      = new TGeoBBox("box1", 6,10,10);
  TGeoBBox *bigBox    = new TGeoBBox("bigBox", 20,10,10);
  TGeoVolume *vbox1   = new TGeoVolume("vbox1", box1, air);
  TGeoVolume *vBigBox = new TGeoVolume("vBigBox", bigBox, air);
  vbox1->SetVisibility(kFALSE);
  vBigBox->SetVisibility(kFALSE);

  TGeoTranslation *tr1 = new TGeoTranslation("negTr",-14,0,0);
  vBigBox->AddNode(vbox1, 1, tr1);
  moth->AddNode(vBigBox, 1, 0);

  // **************************************************
  // Inserting a round cable (or here a water pipe...)
  // **************************************************

  Int_t waterColor = 7;
  Int_t aluColor = 5;
  AliITSv11GeomCableRound roundCable("waterPipe", 0.9); //radius of 0.9cm
  roundCable.SetNLayers(2); 
  roundCable.SetLayer(0, 0.7, water, waterColor); // radius of 0.7cm
  roundCable.SetLayer(1, 0.2, alu, aluColor);     // thickness of 0.2cm

  // ****** Set check points and their containers ******
  // The 2 first points are in the small box (vbox1)
  // The second point is at the boundary 

  Double_t coord0[3] = {0,-2,-2};
  Double_t coord1[3] = {6,2,1};
  Double_t vect0[3]  = {1,1,0};
  Double_t vect1[3]  = {1,0,0};
  // coordinates have to be given in the specified container
  // reference system (here it's going to be vbox1).
  // vect1 and vect2 are vectors perpendicular to the segment ends
  // (These vectors don't need to be normalized)
  roundCable.AddCheckPoint( vbox1, 0, coord0, vect0);
  roundCable.AddCheckPoint( vbox1, 1, coord1, vect1);

  // Then, let's cross the boundary ! You just need
  // to put the next point in the other volume, vBigBox.
  // At the moment of creating the second segment, it will
  // be inserted in this volume. That is why the point 1 had to
  // be at the boundary, because otherwise the second segment
  // between de points 1 and 2 would have been inserted in the
  // vBigBox but in the same time would have cross its
  // boundary ...
  Double_t coord2[3] = {-2,6,4}; // coord. syst. of vBigBox !
  Double_t vect2[3]= {1,1,0.5};
  roundCable.AddCheckPoint( vBigBox, 2, coord2, vect2);

  Double_t coord3[3] = {4,6,4};
  Double_t vect3[3]= {-1,0,0};
  roundCable.AddCheckPoint( vBigBox, 3, coord3, vect3);

  Double_t coord4[3] = {4,0,-4};
  Double_t vect4[3]= {1,0,0};
  roundCable.AddCheckPoint( vBigBox, 4, coord4, vect4);

  Double_t coord5[3] = {4,-6,4};
  Double_t vect5[3]= {1,0,0};
  roundCable.AddCheckPoint( vBigBox, 5, coord5, vect5);

  Double_t coord6[3] = {7,-6,4};
  Double_t vect6[3]= {1,0,0};
  roundCable.AddCheckPoint( vBigBox, 6, coord6, vect6);

  Double_t r = 7;
  Double_t angle = 70*TMath::DegToRad(); 
  Double_t coord7[3] = {coord6[0] +r*sin(angle), coord6[1],
			coord6[2] -r*(1-cos(angle)) };
  Double_t vect7[3]= {r*cos(angle),0,-r*sin(angle)};
  roundCable.AddCheckPoint( vBigBox, 7, coord7, vect7);

  Double_t coord8[3] = { coord7[0]+vect7[0], coord7[1]+vect7[1],-10};
  Double_t vect8[3]= {0,0,1};
  roundCable.AddCheckPoint( vBigBox, 8, coord8, vect8);

  // ****** Creating the corresponding volume ******
  // Since the container volumes of the check points have
  // been recorded, this can be done at any moments, providing
  // that the container volumes are found in the sub-nodes
  // of the initial node (the top volume of the TGeoManager or
  // the volume set in SetInitialNode(TGeoVolume*) function)

  roundCable.SetInitialNode(vBigBox); //Set the root node
  roundCable.CreateAndInsertCableSegment( 1);
  // This command means : create the segment between point 0
  // and point 1. The segment is automatically inserted in the
  // container volume of point 1.

  roundCable.CreateAndInsertCableSegment( 2);
  roundCable.CreateAndInsertCableSegment( 3);

  // The following segment is going to be a torus segment.
  // The radius and position of the torus is defined by the
  // orthogonal vector of point 4 (the orientation of this vector
  // and the position of the 2 check points are enough to define
  // completely the torus)
  roundCable.CreateAndInsertTorusSegment( 4, 180);
  // The second argument is an additionnal rotation of the
  // segment around the axis defined by the 2 check points.

  roundCable.CreateAndInsertTorusSegment( 5);
  roundCable.CreateAndInsertCableSegment( 6);
  roundCable.CreateAndInsertTorusSegment( 7,180);
  roundCable.CreateAndInsertCableSegment( 8);

*/



ClassImp(AliITSv11GeomCableRound)

//________________________________________________________________________
AliITSv11GeomCableRound::
AliITSv11GeomCableRound(const char* name, Double_t radius) :
  AliITSv11GeomCable(name),
  fRadius(radius),
  fNlayer(0),
  fPhiMin(0),
  fPhiMax(360)
 {
  // Constructor
  for (Int_t i=0; i<fgkCableMaxLayer ; i++) {
    fLayThickness[i] = 0;
    fLayColor[i] = 0;
    fLayMedia[i] = 0;
  };
}

//________________________________________________________________________
AliITSv11GeomCableRound::AliITSv11GeomCableRound(const AliITSv11GeomCableRound &s) :
  AliITSv11GeomCable(s),fRadius(s.fRadius),fNlayer(s.fNlayer),fPhiMin(s.fPhiMin),
  fPhiMax(s.fPhiMax)
{
  //     Copy Constructor 
  for (Int_t i=0; i<s.fNlayer; i++) {
    fLayThickness[i] = s.fLayThickness[i];
    fLayMedia[i] = s.fLayMedia[i];
    fLayColor[i] = s.fLayColor[i];
  }
}

//________________________________________________________________________
AliITSv11GeomCableRound& AliITSv11GeomCableRound::
operator=(const AliITSv11GeomCableRound &s) {
  //     Assignment operator
  // Not fully inplemented yet !!!

  if(&s == this) return *this;
  *this = s;
  fRadius = s.fRadius;
  fPhiMin = s.fPhiMin;
  fPhiMax = s.fPhiMax;
  fNlayer = s.fNlayer;
  for (Int_t i=0; i<s.fNlayer; i++) {
    fLayThickness[i] = s.fLayThickness[i];
    fLayMedia[i] = s.fLayMedia[i];
    fLayColor[i] = s.fLayColor[i];
  };
  return *this;
}

//________________________________________________________________________
Int_t AliITSv11GeomCableRound::GetPoint( Int_t iCheckPt, Double_t *coord)
  const {
  // Get check point #iCheckPt
  TVectorD *coordVector =(TVectorD *)fPointArray.UncheckedAt(2*iCheckPt);
#if ROOT_VERSION_CODE < ROOT_VERSION(4,0,0)
  CopyFrom(coord, coordVector->GetElements());
#else
  CopyFrom(coord, coordVector->GetMatrixArray());
#endif 
  return kTRUE;
}

//________________________________________________________________________
Int_t AliITSv11GeomCableRound::GetVect( Int_t iCheckPt, Double_t *coord)
  const {
  //
  // Get vector transverse to the section at point #iCheckPt
  //

  TVectorD *coordVector =(TVectorD *)fPointArray.UncheckedAt(2*iCheckPt+1);
#if ROOT_VERSION_CODE < ROOT_VERSION(4,0,0)
  CopyFrom(coord, coordVector->GetElements());
#else
  CopyFrom(coord, coordVector->GetMatrixArray());
#endif 
  return kTRUE;
}

//________________________________________________________________________
void AliITSv11GeomCableRound::AddCheckPoint( TGeoVolume *vol, Int_t iCheckPt,
					       Double_t *coord, Double_t *orthVect)
{
  //
  // Add point #iCheckPt and its transverse vector. Point is added at (i) in
  // fPointArray and the vector is added at (i+1)
  //


  if (iCheckPt>=fVolumeArray.GetEntriesFast()) {
    fVolumeArray.AddLast(vol);
    TVectorD *point = new TVectorD(3,coord);
    TVectorD *vect  = new TVectorD(3,orthVect);
    fPointArray.AddLast(point);
    fPointArray.AddLast(vect);

  } else if ((iCheckPt >= 0)&&(iCheckPt < fVolumeArray.GetEntriesFast())) {
    fVolumeArray.AddAt(vol, iCheckPt);
    TVectorD *point = new TVectorD(3,coord);
    TVectorD *vect  = new TVectorD(3,orthVect);
    fPointArray.AddAt(point, iCheckPt*2  );
    fPointArray.AddAt(vect,  iCheckPt*2+1);
  };
}

//________________________________________________________________________
void AliITSv11GeomCableRound::PrintCheckPoints() const {
  // Print all check points

  printf("  ---\n  Printing all check points of the round cable\n");
  for (Int_t i = 0; i<fVolumeArray.GetEntriesFast(); i++) {
    TVectorD *coordVector = (TVectorD *)fPointArray.UncheckedAt(i*2);
    //TVectorD *vectVector = (TVectorD *)fPointArray.UncheckedAt(i*2+1);
    Double_t coord[3];
#if ROOT_VERSION_CODE < ROOT_VERSION(4,0,0)
    CopyFrom(coord, coordVector->GetElements());
#else
    CopyFrom(coord, coordVector->GetMatrixArray());
#endif 
    printf("   ( %.2f, %.2f, %.2f )\n", coord[0], coord[1], coord[2]);
  };
}


//________________________________________________________________________
Int_t AliITSv11GeomCableRound::CreateAndInsertCableSegment(Int_t p2)
{
//    Creates a cable segment between points p1 and p2.
//
// The segment volume is created inside the volume containing point2
// Therefore this segment should be defined in this volume only.
// I mean here that, if the previous point is in another volume,
// it should be just at the border between the 2 volumes. Also the
// orientation vector of the previous point should be othogonal to
// the surface between the 2 volumes.

  TGeoNode *mainNode;
  if (fInitialNode==0) {
    TObjArray *nodes = gGeoManager->GetListOfNodes();
    if (nodes->GetEntriesFast()==0) return kFALSE;
    mainNode = (TGeoNode *) nodes->UncheckedAt(0);
  } else {
    mainNode = fInitialNode;
  };

  Int_t p1 = p2 - 1;
  TGeoVolume *p1Vol = GetVolume(p1);
  TGeoVolume *p2Vol = GetVolume(p2);

  ResetCheckDaughter();
  fCurrentVol = p1Vol;
  if (! CheckDaughter(mainNode)) {
    printf("Error::volume containing point is not visible in node tree!\n");
    return kFALSE;
  };

  Double_t coord1[3], coord2[3], vect1[3], vect2[3];
  //=================================================
  // Get p1 position in the systeme of p2
  if (p1Vol!=p2Vol) {

    Int_t p1nodeInd[fgkCableMaxNodeLevel]; 
    for (Int_t i=0; i<fgkCableMaxNodeLevel; i++) p1nodeInd[i]=fNodeInd[i];
    Int_t p1volLevel = 0;
    while (p1nodeInd[p1volLevel]!=-1) p1volLevel++;
    p1volLevel--;

    ResetCheckDaughter();
    fCurrentVol = p2Vol;
    if (! CheckDaughter(mainNode)) {
      printf("Error::volume containing point is not visible in node tree!\n");
      return kFALSE;
    };
    Int_t p2nodeInd[fgkCableMaxNodeLevel];
    for (Int_t i=0; i<fgkCableMaxNodeLevel; i++) p2nodeInd[i]=fNodeInd[i];
    Int_t commonMotherLevel = 0;
    while (p1nodeInd[commonMotherLevel]==fNodeInd[commonMotherLevel])
      commonMotherLevel++;
    commonMotherLevel--;
    Int_t p2volLevel = 0;
    while (fNodeInd[p2volLevel]!=-1) p2volLevel++;
    p2volLevel--;

    // Get coord and vect of p1 in the common mother reference system
    GetCheckPoint(p1, 0, p1volLevel-commonMotherLevel, coord1);
    GetCheckVect( p1, 0, p1volLevel-commonMotherLevel, vect1);
    // Translate them in the reference system of the volume containing p2    
    TGeoNode *pathNode[fgkCableMaxNodeLevel];
    pathNode[0] = mainNode;
    for (Int_t i=0; i<=p2volLevel; i++) {
      pathNode[i+1] = pathNode[i]->GetDaughter(p2nodeInd[i]);
    };
    Double_t globalCoord1[3] = {coord1[0], coord1[1], coord1[2]}; 
    Double_t globalVect1[3]  = {vect1[0], vect1[1], vect1[2]};

    for (Int_t i = commonMotherLevel+1; i<=p2volLevel; i++) {
      pathNode[i+1]->GetMatrix()->MasterToLocal(globalCoord1, coord1);
      pathNode[i+1]->GetMatrix()->MasterToLocalVect(globalVect1, vect1);
      CopyFrom(globalCoord1, coord1);
      CopyFrom(globalVect1, vect1);
    };
  } else {
    GetCheckPoint(p1, 0, 0, coord1);
    GetCheckVect(p1, 0, 0, vect1);
  };
  
  //=================================================
  // Get p2 position in the systeme of p2
  GetCheckPoint(p2, 0, 0, coord2);
  GetCheckVect(p2, 0, 0, vect2);

  Double_t cx = (coord1[0]+coord2[0])/2;
  Double_t cy = (coord1[1]+coord2[1])/2;
  Double_t cz = (coord1[2]+coord2[2])/2;
  Double_t dx = coord2[0]-coord1[0];
  Double_t dy = coord2[1]-coord1[1];
  Double_t dz = coord2[2]-coord1[2];

  //=================================================
  // Positionning of the segment between the 2 points
  if ((dy<1e-31)&&(dy>0)) dy = 1e-31;
  if ((dz<1e-31)&&(dz>0)) dz = 1e-31;
  if ((dy>-1e-31)&&(dy<0)) dy = -1e-31;
  if ((dz>-1e-31)&&(dz<0)) dz = -1e-31;

  Double_t angleRot1 = -TMath::ATan2(dx,dy);
  Double_t planDiagL = TMath::Sqrt(dy*dy+dx*dx);
  Double_t angleRotDiag = -TMath::ATan2(planDiagL,dz);
  TGeoRotation *rot = new TGeoRotation("",angleRot1*TMath::RadToDeg(),
				       angleRotDiag*TMath::RadToDeg(),
				       0);
  Double_t localVect1[3], localVect2[3];
  rot->MasterToLocalVect(vect1, localVect1);
  rot->MasterToLocalVect(vect2, localVect2);
  TGeoTranslation *trans = new TGeoTranslation("",cx, cy, cz);

  //=================================================
  // Create the segment and add it to the mother volume
  TGeoVolume *vCableSeg = CreateSegment(coord1, coord2,
					localVect1, localVect2, p2);

  TGeoCombiTrans  *combi = new TGeoCombiTrans(*trans, *rot);
  p2Vol->AddNode(vCableSeg, p2, combi);
  //=================================================
  delete rot;
  delete trans;

  if (fDebug) {
    printf("---\n  Cable segment points : ");
    printf("%f, %f, %f\n",coord1[0], coord1[1], coord1[2]);
    printf("%f, %f, %f\n",coord2[0], coord2[1], coord2[2]);
  };
//   #include <TGeoSphere.h>
//   TGeoMedium *airSDD = gGeoManager->GetMedium("ITS_ITSsddAir");
//   TGeoSphere *sphere = new TGeoSphere(0, 0.15);
//   TGeoVolume *vSphere = new TGeoVolume("", sphere, airSDD);
//   TGeoTranslation *trC = new TGeoTranslation("", cx, cy, cz);
//   TGeoTranslation *tr1 = new TGeoTranslation("",coord1[0],
// 					     coord1[1],coord1[2]);
//   TGeoTranslation *tr2 = new TGeoTranslation("",coord2[0],
// 					     coord2[1],coord2[2]);
//   p2Vol->AddNode(vSphere, p2*3-2, trC);
//   p2Vol->AddNode(vSphere, p2*3-1, tr1);
//   p2Vol->AddNode(vSphere, p2*3  , tr2);

  return kTRUE;
}


//________________________________________________________________________
Int_t AliITSv11GeomCableRound::CreateAndInsertTorusSegment(Int_t p2, Double_t rotation)
{
  // Create a torus cable segment between points p1 and p2.
  // The radius and position of the torus is defined by the
  // perpendicular vector of point p2 (the orientation of this vector
  // and the position of the 2 check points are enough to completely
  // define the torus)

  TGeoNode *mainNode;
  if (fInitialNode==0) {
    TObjArray *nodes = gGeoManager->GetListOfNodes();
    if (nodes->GetEntriesFast()==0) return kFALSE;
    mainNode = (TGeoNode *) nodes->UncheckedAt(0);
  } else {
    mainNode = fInitialNode;
  };

  Int_t p1 = p2 - 1;
  TGeoVolume *p1Vol = GetVolume(p1);
  TGeoVolume *p2Vol = GetVolume(p2);

  ResetCheckDaughter();
  fCurrentVol = p1Vol;
  if (! CheckDaughter(mainNode)) {
    printf("Error::volume containing point is not visible in node tree!\n");
    return kFALSE;
  };

  Double_t coord1[3], coord2[3], vect1[3], vect2[3];
  //=================================================
  // Get p1 position in the systeme of p2
  if (p1Vol!=p2Vol) {

    Int_t p1nodeInd[fgkCableMaxNodeLevel]; 
    for (Int_t i=0; i<fgkCableMaxNodeLevel; i++) p1nodeInd[i]=fNodeInd[i];
    Int_t p1volLevel = 0;
    while (p1nodeInd[p1volLevel]!=-1) p1volLevel++;
    p1volLevel--;

    ResetCheckDaughter();
    fCurrentVol = p2Vol;
    if (! CheckDaughter(mainNode)) {
      printf("Error::volume containing point is not visible in node tree!\n");
      return kFALSE;
    };
    Int_t p2nodeInd[fgkCableMaxNodeLevel];
    for (Int_t i=0; i<fgkCableMaxNodeLevel; i++) p2nodeInd[i]=fNodeInd[i];
    Int_t commonMotherLevel = 0;
    while (p1nodeInd[commonMotherLevel]==fNodeInd[commonMotherLevel])
      commonMotherLevel++;
    commonMotherLevel--;
    Int_t p2volLevel = 0;
    while (fNodeInd[p2volLevel]!=-1) p2volLevel++;
    p2volLevel--;

    // Get coord and vect of p1 in the common mother reference system
    GetCheckPoint(p1, 0, p1volLevel-commonMotherLevel, coord1);
    GetCheckVect( p1, 0, p1volLevel-commonMotherLevel, vect1);
    // Translate them in the reference system of the volume containing p2    
    TGeoNode *pathNode[fgkCableMaxNodeLevel];
    pathNode[0] = mainNode;
    for (Int_t i=0; i<=p2volLevel; i++) {
      pathNode[i+1] = pathNode[i]->GetDaughter(p2nodeInd[i]);
    };
    Double_t globalCoord1[3] = {coord1[0], coord1[1], coord1[2]}; 
    Double_t globalVect1[3]  = {vect1[0], vect1[1], vect1[2]};

    for (Int_t i = commonMotherLevel+1; i<=p2volLevel; i++) {
      pathNode[i+1]->GetMatrix()->MasterToLocal(globalCoord1, coord1);
      pathNode[i+1]->GetMatrix()->MasterToLocalVect(globalVect1, vect1);
      CopyFrom(globalCoord1, coord1);
      CopyFrom(globalVect1, vect1);
    };
  } else {
    GetCheckPoint(p1, 0, 0, coord1);
    GetCheckVect(p1, 0, 0, vect1);
  };
  
  //=================================================
  // Get p2 position in the systeme of p2
  GetCheckPoint(p2, 0, 0, coord2);
  GetCheckVect(p2, 0, 0, vect2);

  Double_t cx = (coord1[0]+coord2[0])/2;
  Double_t cy = (coord1[1]+coord2[1])/2;
  Double_t cz = (coord1[2]+coord2[2])/2;
  Double_t dx = coord2[0]-coord1[0];
  Double_t dy = coord2[1]-coord1[1];
  Double_t dz = coord2[2]-coord1[2];
  Double_t length = TMath::Sqrt(dx*dx+dy*dy+dz*dz);

  //=================================================
  // Positionning of the segment between the 2 points
  if ((dy<1e-31)&&(dy>0)) dy = 1e-31;
  if ((dz<1e-31)&&(dz>0)) dz = 1e-31;
  if ((dy>-1e-31)&&(dy<0)) dy = -1e-31;
  if ((dz>-1e-31)&&(dz<0)) dz = -1e-31;

  Double_t angleRot1 = -TMath::ATan2(dx,dy);
  Double_t planDiagL = TMath::Sqrt(dy*dy+dx*dx);
  Double_t angleRotDiag = -TMath::ATan2(planDiagL,dz);

  TGeoRotation rotTorusTemp("",angleRot1*TMath::RadToDeg(),
			    angleRotDiag*TMath::RadToDeg(),0);
  TGeoRotation rotTorusToZ("",0,90,0);
  rotTorusTemp.MultiplyBy(&rotTorusToZ, kTRUE);
  Double_t localVect2[3];
  rotTorusTemp.MasterToLocalVect(vect2, localVect2);
  if (localVect2[1]<0) {
    localVect2[0] = -localVect2[0];
    localVect2[1] = -localVect2[1];
    localVect2[2] = -localVect2[2];
  };
  Double_t normVect2 = TMath::Sqrt(localVect2[0]*localVect2[0]+
				   localVect2[1]*localVect2[1]+
				   localVect2[2]*localVect2[2]);
  Double_t axisX[3] = {1,0,0};
  Double_t cosangleTorusSeg = (localVect2[0]*axisX[0]+
			       localVect2[1]*axisX[1]+
			       localVect2[2]*axisX[2])/normVect2;
  Double_t angleTorusSeg = TMath::ACos(cosangleTorusSeg)*TMath::RadToDeg();
  TGeoRotation rotTorus("",angleRot1*TMath::RadToDeg(),
			angleRotDiag*TMath::RadToDeg(),
			180-angleTorusSeg+rotation);
  rotTorus.MultiplyBy(&rotTorusToZ, kTRUE);
  rotTorus.MasterToLocalVect(vect2, localVect2);
  if (localVect2[1]<0) {
    localVect2[0] = -localVect2[0];
    localVect2[1] = -localVect2[1];
    localVect2[2] = -localVect2[2];
  };
  normVect2 = TMath::Sqrt(localVect2[0]*localVect2[0]+
			  localVect2[1]*localVect2[1]+
			  localVect2[2]*localVect2[2]);
  Double_t axisY[3] = {0,1,0};
  Double_t cosPhi = (localVect2[0]*axisY[0]+localVect2[1]*axisY[1]+
		     localVect2[2]*axisY[2])/normVect2;
  Double_t torusPhi1 = TMath::ACos(cosPhi);
  Double_t torusR = (length/2)/TMath::Sin(torusPhi1);
  torusPhi1 = torusPhi1*TMath::RadToDeg();
  Double_t perpLength = TMath::Sqrt(torusR*torusR-length*length/4);
  Double_t localTransT[3] = {-perpLength,0,0};
  Double_t globalTransT[3];
  rotTorus.LocalToMasterVect(localTransT, globalTransT);
  TGeoTranslation transTorus("",cx+globalTransT[0],cy+globalTransT[1],
			     cz+globalTransT[2]);

  TGeoCombiTrans  *combiTorus = new TGeoCombiTrans(transTorus, rotTorus);

  //=================================================
  // Create the segment and add it to the mother volume
  TGeoVolume *vCableSegT = CreateTorus(torusPhi1, torusR, p2);
  p2Vol->AddNode(vCableSegT, p2, combiTorus);

  if (fDebug) {
    printf("---\n  Cable segment points : ");
    printf("%f, %f, %f\n",coord1[0], coord1[1], coord1[2]);
    printf("%f, %f, %f\n",coord2[0], coord2[1], coord2[2]);
  };

  return kTRUE;
}

//________________________________________________________________________
TGeoVolume *AliITSv11GeomCableRound::CreateSegment( Double_t *coord1,
						      Double_t *coord2,
						      Double_t *localVect1,
						      Double_t *localVect2, Int_t p)
{
  // Create one cylindrical segment and its layers

  //=================================================
  // Calculate segment "deformation"
  Double_t dx = coord2[0]-coord1[0];
  Double_t dy = coord2[1]-coord1[1];
  Double_t dz = coord2[2]-coord1[2];
  Double_t length = TMath::Sqrt(dx*dx+dy*dy+dz*dz);

  // normal vectors have to point outside the TGeoCtub :
  if (-localVect1[2]<0) {
    localVect1[0] = -localVect1[0];
    localVect1[1] = -localVect1[1];
    localVect1[2] = -localVect1[2];
  };
  if (localVect2[2]<0) {
    localVect2[0] = -localVect2[0];
    localVect2[1] = -localVect2[1];
    localVect2[2] = -localVect2[2];
  };
  //=================================================
  // Create the segment
  TGeoCtub *cableSeg = new TGeoCtub(0, fRadius, length/2, fPhiMin, fPhiMax,
				    localVect1[0],localVect1[1],localVect1[2],
				    localVect2[0],localVect2[1],localVect2[2]);

  TGeoMedium *airSDD = gGeoManager->GetMedium("ITS_ITSair");
  char name[100];
  sprintf(name, "%s_%i",GetName(),p);
  TGeoVolume *vCableSeg = new TGeoVolume(name, cableSeg, airSDD);

  // add all cable layers
  Double_t layThickness[100+1];                        // 100 layers max !!!
  layThickness[0] = 0;
  for (Int_t iLay=0; iLay<fNlayer; iLay++) {
    
    layThickness[iLay+1] = fLayThickness[iLay]+layThickness[iLay];
    TGeoCtub *lay = new TGeoCtub(layThickness[iLay], layThickness[iLay+1],
				 length/2, fPhiMin, fPhiMax,
				 localVect1[0],localVect1[1],localVect1[2],
				 localVect2[0],localVect2[1],localVect2[2]);

    TGeoVolume *vLay = new TGeoVolume("vCableSegLay", lay, fLayMedia[iLay]);
    vLay->SetLineColor(fLayColor[iLay]);
    vCableSeg->AddNode(vLay, iLay+1, 0);
  };

  vCableSeg->SetVisibility(kFALSE);
  return vCableSeg;
}


//________________________________________________________________________
TGeoVolume *AliITSv11GeomCableRound::CreateTorus( Double_t &phi,
						  Double_t &r, Int_t p)
{
  // Create one torus segment and its layers

  Double_t torusR = r;
//   Double_t torusPhi1 = phi;
//   Double_t torusDPhi = -2*torusPhi1;  // bug in root ...
  Double_t torusPhi1 = 360-phi;
  Double_t torusDPhi = 2*phi;

  //=================================================
  // Create the segment
  TGeoTorus *cableSeg = new TGeoTorus(torusR, 0, fRadius, torusPhi1, torusDPhi);
  TGeoMedium *airSDD = gGeoManager->GetMedium("ITS_ITSair");
  char name[100];
  sprintf(name, "%s_%i",GetName(),p);
  TGeoVolume *vCableSeg = new TGeoVolume(name, cableSeg, airSDD);

  // add all cable layers
  Double_t layThickness[100+1];                        // 100 layers max !!!
  layThickness[0] = 0;
  for (Int_t iLay=0; iLay<fNlayer; iLay++) {
    
    layThickness[iLay+1] = fLayThickness[iLay]+layThickness[iLay];
    TGeoTorus *lay = new TGeoTorus(torusR, layThickness[iLay],
				   layThickness[iLay+1],
				   torusPhi1,torusDPhi);

    TGeoVolume *vLay = new TGeoVolume("vCableSegLay",lay,fLayMedia[iLay]);
    vLay->SetLineColor(fLayColor[iLay]);
    vCableSeg->AddNode(vLay, iLay+1, 0);
  };

  vCableSeg->SetVisibility(kFALSE);
  return vCableSeg;
}


//________________________________________________________________________
void AliITSv11GeomCableRound::SetNLayers(Int_t nLayers) {
  // Set the total number of layers
  if((nLayers>0) &&(nLayers<=fgkCableMaxLayer)) {
    fNlayer = nLayers;
    for (Int_t i = 0; i<fNlayer; i++) {
      fLayThickness[i] = 0;
      fLayMedia[i] = 0;
    };
  };
}

//________________________________________________________________________
Int_t AliITSv11GeomCableRound::SetLayer(Int_t nLayer, Double_t thick,
					   TGeoMedium *medium, Int_t color) {
  // Set layer #nLayer
  if ((nLayer<0)||(nLayer>=fNlayer)) {
    printf("Set wrong layer number of the cable\n");
    return kFALSE;
  };
  if (nLayer>0)
    if (fLayThickness[nLayer-1]<=0) {
      printf("You must define cable layer %i first !",nLayer-1);
      return kFALSE;
    };

  Double_t thickTot = 0;
  for (Int_t i=0; i<nLayer; i++) thickTot += fLayThickness[i];
  thickTot += thick;
  if (thickTot-1e-10>fRadius) {
    printf("Can't add this layer, cable thickness would be higher than total\n");
    return kFALSE;
  };

  fLayThickness[nLayer] = thick;
  fLayMedia[nLayer] = medium;
  fLayColor[nLayer] = color;

  return kTRUE;
}
