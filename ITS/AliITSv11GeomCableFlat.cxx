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
//   Class for flat cables
//
// Ludovic Gaudichet                                   gaudichet@to.infn.it
//*************************************************************************



// General Root includes
//#include <Riostream.h>
#include <TMath.h>
#include <TVectorD.h>

// Root Geometry includes
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoArb8.h>
#include <TGeoTube.h>
#include <TGeoMatrix.h>
#include <TGeoNode.h>

#include "AliITSv11GeomCableFlat.h"


ClassImp(AliITSv11GeomCableFlat)

//________________________________________________________________________
AliITSv11GeomCableFlat::AliITSv11GeomCableFlat():
  AliITSv11GeomCable(),
  fWidth(0),
  fThick(0),
  fNlayer(0)
{
  // constructor
  for (Int_t i=0; i<fgkCableMaxLayer ; i++) {
    fLayThickness[i] = 0;
    fTranslation[i]  = 0;
    fLayColor[i]     = 0;
    fLayMedia[i]     = 0;  
 };
}

//________________________________________________________________________
AliITSv11GeomCableFlat::
AliITSv11GeomCableFlat(const char* name, Double_t width, Double_t thick) :
  AliITSv11GeomCable(name),
  fWidth(width),
  fThick(thick),
  fNlayer(0)
 {
  // standard constructor
  for (Int_t i=0; i<fgkCableMaxLayer ; i++) {
    fLayThickness[i] = 0;
    fTranslation[i]  = 0;
    fLayColor[i]     = 0;
    fLayMedia[i]     = 0;  
  }; 
}

//________________________________________________________________________
AliITSv11GeomCableFlat::AliITSv11GeomCableFlat(const AliITSv11GeomCableFlat &s) :
  AliITSv11GeomCable(s),fWidth(s.fWidth),fThick(s.fThick),fNlayer(s.fNlayer)
{
  //     Copy Constructor 
  for (Int_t i=0; i<s.fNlayer; i++) {
    fLayThickness[i] = s.fLayThickness[i];
    fTranslation[i] = s.fTranslation[i];
    fLayMedia[i] = s.fLayMedia[i];
    fLayColor[i] = s.fLayColor[i];
  }
}

//________________________________________________________________________
AliITSv11GeomCableFlat& AliITSv11GeomCableFlat::
operator=(const AliITSv11GeomCableFlat &s) {
  //     Assignment operator
  // Not fully inplemented yet !!!

  if(&s == this) return *this;
  *this = s;
  fWidth = s.fWidth;
  fThick = s.fThick;
  fNlayer = s.fNlayer;
  for (Int_t i=0; i<s.fNlayer; i++) {
    fLayThickness[i] = s.fLayThickness[i];
    fTranslation[i] = s.fTranslation[i];
    fLayMedia[i] = s.fLayMedia[i];
    fLayColor[i] = s.fLayColor[i];
  };
  return *this;
}

//________________________________________________________________________
Int_t AliITSv11GeomCableFlat::GetPoint( Int_t iCheckPt, Double_t *coord)
  const {
  // Get the correct point #iCheckPt
  TVectorD *coordVector =(TVectorD *)fPointArray.At(2*iCheckPt);
  if (coordVector) {
#if ROOT_VERSION_CODE < ROOT_VERSION(4,0,0)
    CopyFrom(coord, coordVector->GetElements());
#else
    CopyFrom(coord, coordVector->GetMatrixArray());
#endif 
    return kTRUE;
  } else {
    return kFALSE;
  };
}


//________________________________________________________________________
Int_t AliITSv11GeomCableFlat::GetVect( Int_t iCheckPt, Double_t *coord)
  const {
  // Get the correct vect corresponding to point #iCheckPt

  TVectorD *coordVector =(TVectorD *)fPointArray.At(2*iCheckPt+1);
  if (coordVector) {
#if ROOT_VERSION_CODE < ROOT_VERSION(4,0,0)
    CopyFrom(coord, coordVector->GetElements());
#else
    CopyFrom(coord, coordVector->GetMatrixArray());
#endif 
    return kTRUE;
  } else {
    return kFALSE;
  };
}


//________________________________________________________________________
void AliITSv11GeomCableFlat::AddCheckPoint( TGeoVolume *vol, Int_t iCheckPt,
					    Double_t *coord, Double_t *orthVect)
{
  //
  // Add a check point. In the fPointArray, the point is at i and its vector
  // is at i+1.
  //

//   if (iCheckPt>=fVolumeArray.GetEntriesFast()) {
//     fVolumeArray.AddLast(vol);
//     TVectorD *point = new TVectorD(3,coord);
//     TVectorD *vect  = new TVectorD(3,orthVect);
//     fPointArray.AddLast(point);
//     fPointArray.AddLast(vect);

//   } else if ((iCheckPt >= 0)&&(iCheckPt < fVolumeArray.GetEntriesFast())) {
//     fVolumeArray.AddAt(vol, iCheckPt);
//     TVectorD *point = new TVectorD(3,coord);
//     TVectorD *vect  = new TVectorD(3,orthVect);
//     fPointArray.AddAt(point, iCheckPt*2  );
//     fPointArray.AddAt(vect,  iCheckPt*2+1);
//   };
  fVolumeArray.AddAtAndExpand(vol, iCheckPt);
  TVectorD *point = new TVectorD(3,coord);
  TVectorD *vect  = new TVectorD(3,orthVect);
  fPointArray.AddAtAndExpand(point, iCheckPt*2  );
  fPointArray.AddAtAndExpand(vect,  iCheckPt*2+1);
}

//________________________________________________________________________
void AliITSv11GeomCableFlat::PrintCheckPoints() const {
  // print all check points of the cable
  printf("  ---\n  Printing all check points of the flat cable\n");
  for (Int_t i = 0; i<fVolumeArray.GetEntriesFast(); i++) {
     Double_t coord[3];
     if (GetPoint( i, coord))
       printf("   ( %.2f, %.2f, %.2f )\n", coord[0], coord[1], coord[2]);
  };
}

//________________________________________________________________________
TGeoVolume* AliITSv11GeomCableFlat::CreateAndInsertCableSegment(Int_t p2,
							  Double_t rotation)
{
//    Creates a cable segment between points p1 and p2.
//    Rotation is the eventual rotation of the flat cable
//    along its length axis
//
// The segment volume is created inside the volume containing point2
// Therefore this segment should be defined in this volume only.
// I mean here that, if the previous point is in another volume,
// it should be just at the border between the 2 volumes. Also the
// orientation vector of the previous point should be orthogonal to
// the surface between the 2 volumes.

  TGeoNode *mainNode;
  if (fInitialNode==0) {
    TObjArray *nodes = gGeoManager->GetListOfNodes();
    if (nodes->GetEntriesFast()==0) return 0;
    mainNode = (TGeoNode *) nodes->UncheckedAt(0);
  } else {
    mainNode = fInitialNode;
  };

  Int_t p1 = p2 - 1;
  TGeoVolume *p2Vol = GetVolume(p2);
  TGeoVolume *p1Vol = GetVolume(p1);

  ResetCheckDaughter();
  fCurrentVol = p1Vol;
  if (! CheckDaughter(mainNode)) {
    printf("Error::volume containing point is not visible in node tree!\n");
    return 0;
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
      return 0;
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
    if (! GetCheckPoint(p1, 0, p1volLevel-commonMotherLevel, coord1) )
      return 0;
    if (! GetCheckVect( p1, 0, p1volLevel-commonMotherLevel, vect1) )
      return 0;

    // Translate them in the reference system of the volume containing p2    
    TGeoNode *pathNode[fgkCableMaxNodeLevel];
    pathNode[0] = mainNode;
    for (Int_t i=0; i<=p2volLevel; i++) {
      pathNode[i+1] = pathNode[i]->GetDaughter(p2nodeInd[i]);
    };
    Double_t globalCoord1[3] = {coord1[0], coord1[1], coord1[2]}; 
    Double_t globalVect1[3]  = {vect1[0], vect1[1], vect1[2]};

    for (Int_t i = commonMotherLevel+1; i <= p2volLevel; i++) {
      pathNode[i+1]->GetMatrix()->MasterToLocal(globalCoord1, coord1);
      pathNode[i+1]->GetMatrix()->MasterToLocalVect(globalVect1, vect1);
      CopyFrom(globalCoord1, coord1);
      CopyFrom(globalVect1, vect1);
    };
  } else {
    if (! GetCheckPoint(p1, 0, 0, coord1) ) return 0;
    if (! GetCheckVect(p1, 0, 0, vect1) ) return 0;
  };
  
  //=================================================
  // Get p2 position in the systeme of p2
  if (! GetCheckPoint(p2, 0, 0, coord2) ) return 0;
  if (! GetCheckVect(p2, 0, 0, vect2) ) return 0;

  Double_t cx = (coord1[0]+coord2[0])/2;
  Double_t cy = (coord1[1]+coord2[1])/2;
  Double_t cz = (coord1[2]+coord2[2])/2;
  Double_t dx = coord2[0]-coord1[0];
  Double_t dy = coord2[1]-coord1[1];
  Double_t dz = coord2[2]-coord1[2];

  //=================================================
  // Positionning of the segment between the 2 points
  if (TMath::Abs(dy)<1e-231) dy = 1e-231;
  if (TMath::Abs(dz)<1e-231) dz = 1e-231;
  //Double_t angleRot1 = -TMath::ATan(dx/dy);
  //Double_t planDiagL = -TMath::Sqrt(dy*dy+dx*dx);
  //if (dy<0) planDiagL = -planDiagL;
  //Double_t angleRotDiag = TMath::ATan(planDiagL/dz);

  Double_t angleRot1    = -TMath::ATan2(dx,dy);
  Double_t planDiagL    =  TMath::Sqrt(dy*dy+dx*dx);
  Double_t angleRotDiag = -TMath::ATan2(planDiagL,dz);
  //--- (Calculate rotation of segment on the Z axis)
  //-- Here I'm trying to calculate the rotation to be applied in
  //-- order to match as closer as possible this segment and the 
  //-- previous one. 
  //-- It seems that some times it doesn't work ...
  Double_t angleRotZ = 0;
  TGeoRotation rotTemp("",angleRot1*TMath::RadToDeg(),
		       angleRotDiag*TMath::RadToDeg(), rotation);
  Double_t localX[3] = {0,1,0};
  Double_t globalX[3];
  rotTemp.LocalToMasterVect(localX, globalX);
  CopyFrom(localX, globalX);
  GetCheckVect(localX, p2Vol, 0, fgkCableMaxNodeLevel+1, globalX);
  Double_t orthVect[3];
  GetCheckVect(vect1, p2Vol, 0, fgkCableMaxNodeLevel+1, orthVect);
  if (p2>1) {
    Double_t orthVectNorm2 = ScalProd(orthVect,orthVect);
    Double_t alpha1 = ScalProd(fPreviousX,orthVect)/orthVectNorm2;
    Double_t alpha2 = ScalProd(globalX,orthVect)/orthVectNorm2;
    Double_t globalX1p[3], globalX2p[3];
    globalX1p[0] = fPreviousX[0] - alpha1*orthVect[0];
    globalX1p[1] = fPreviousX[1] - alpha1*orthVect[1];
    globalX1p[2] = fPreviousX[2] - alpha1*orthVect[2];
    globalX2p[0] = globalX[0] - alpha2*orthVect[0];
    globalX2p[1] = globalX[1] - alpha2*orthVect[1];
    globalX2p[2] = globalX[2] - alpha2*orthVect[2];
    //-- now I'm searching the 3th vect which makes an orthogonal base
    //-- with orthVect and globalX1p ...
    Double_t nulVect[3] = {0,0,0};
    Double_t axis3[3];
    TMath::Normal2Plane(nulVect, orthVect, globalX1p, axis3);
    Double_t globalX1pNorm2 = ScalProd(globalX1p, globalX1p);
    Double_t beta = ScalProd(globalX2p, globalX1p)/globalX1pNorm2;
    Double_t gamma = ScalProd(globalX2p, axis3);
    angleRotZ = (TMath::ATan2(1,0) - TMath::ATan2(beta, gamma))
                *TMath::RadToDeg();
  };
  //   cout << "!!!!!!!!!!!!!!!!!!!  angle = " <<angleRotZ << endl;
  CopyFrom(fPreviousX, globalX);
  //---
  Double_t localVect1[3], localVect2[3];
  TGeoRotation rot("",angleRot1*TMath::RadToDeg(),
		   angleRotDiag*TMath::RadToDeg(),
		   rotation);
// 		   rotation-angleRotZ);
// since angleRotZ doesn't always work, I won't use it ...

  rot.MasterToLocalVect(vect1, localVect1);
  rot.MasterToLocalVect(vect2, localVect2);

  //=================================================
  // Create the segment and add it to the mother volume
  TGeoVolume *vCableSegB = CreateSegment(coord1, coord2,
					 localVect1, localVect2);
//   TGeoVolume *vCableSegB = CreateBoxSegment(coord1, coord2);

  TGeoRotation rotArbSeg("", 0, 90, 0);
  rotArbSeg.MultiplyBy(&rot, kFALSE);
  TGeoTranslation trans("",cx, cy, cz);
  TGeoCombiTrans  *combiB = new TGeoCombiTrans(trans, rotArbSeg);
  p2Vol->AddNode(vCableSegB, p2, combiB);
  //=================================================;

  if (fDebug) {
    printf("---\n  Cable segment points : ");
    printf("%f, %f, %f\n",coord1[0], coord1[1], coord1[2]);
    printf("%f, %f, %f\n",coord2[0], coord2[1], coord2[2]);
  };

//   #include <TGeoSphere.h>
//   TGeoMedium *airSDD = gGeoManager->GetMedium("ITS_ITSsddAir");
//   TGeoSphere *sphere = new TGeoSphere(0, 0.05);
//   TGeoVolume *vSphere = new TGeoVolume("", sphere, airSDD);
//   TGeoTranslation *trC = new TGeoTranslation("", cx, cy, cz);
//   TGeoTranslation *tr1 = new TGeoTranslation("",coord1[0],
// 					     coord1[1],coord1[2]);
//   TGeoTranslation *tr2 = new TGeoTranslation("",coord2[0],
// 					     coord2[1],coord2[2]);
//   p2Vol->AddNode(vSphere, p2*3-2, trC);
//   p2Vol->AddNode(vSphere, p2*3-1, tr1);
//   p2Vol->AddNode(vSphere, p2*3  , tr2);

  return vCableSegB;
}

//________________________________________________________________________
TGeoVolume* AliITSv11GeomCableFlat::CreateAndInsertBoxCableSegment(Int_t p2,
							  Double_t rotation)
{
  // This function is to be use only when the segment has the shape
  // of a simple box, i.e. the normal vector to its end is perpendicular
  // to the segment own axis
//    Creates a cable segment between points p1 and p2.
//    Rotation is the eventual rotation of the flat cable
//    along its length axis
//
// The segment volume is created inside the volume containing point2
// Therefore this segment should be defined in this volume only.
// I mean here that, if the previous point is in another volume,
// it should be just at the border between the 2 volumes. Also the
// orientation vector of the previous point should be orthogonal to
// the surface between the 2 volumes.

  TGeoNode *mainNode;
  if (fInitialNode==0) {
    TObjArray *nodes = gGeoManager->GetListOfNodes();
    if (nodes->GetEntriesFast()==0) return 0;
    mainNode = (TGeoNode *) nodes->UncheckedAt(0);
  } else {
    mainNode = fInitialNode;
  };

  Int_t p1 = p2 - 1;
  TGeoVolume *p2Vol = GetVolume(p2);
  TGeoVolume *p1Vol = GetVolume(p1);

  ResetCheckDaughter();
  fCurrentVol = p1Vol;
  if (! CheckDaughter(mainNode)) {
    printf("Error::volume containing point is not visible in node tree!\n");
    return 0;
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
      return 0;
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
    if (! GetCheckPoint(p1, 0, p1volLevel-commonMotherLevel, coord1) )
      return 0;
    if (! GetCheckVect( p1, 0, p1volLevel-commonMotherLevel, vect1) )
      return 0;

    // Translate them in the reference system of the volume containing p2    
    TGeoNode *pathNode[fgkCableMaxNodeLevel];
    pathNode[0] = mainNode;
    for (Int_t i=0; i<=p2volLevel; i++) {
      pathNode[i+1] = pathNode[i]->GetDaughter(p2nodeInd[i]);
    };
    Double_t globalCoord1[3] = {coord1[0], coord1[1], coord1[2]}; 
    Double_t globalVect1[3]  = {vect1[0], vect1[1], vect1[2]};

    for (Int_t i = commonMotherLevel+1; i <= p2volLevel; i++) {
      pathNode[i+1]->GetMatrix()->MasterToLocal(globalCoord1, coord1);
      pathNode[i+1]->GetMatrix()->MasterToLocalVect(globalVect1, vect1);
      CopyFrom(globalCoord1, coord1);
      CopyFrom(globalVect1, vect1);
    };
  } else {
    if (! GetCheckPoint(p1, 0, 0, coord1) ) return 0;
    if (! GetCheckVect(p1, 0, 0, vect1) ) return 0;
  };
  
  //=================================================
  // Get p2 position in the systeme of p2
  if (! GetCheckPoint(p2, 0, 0, coord2) ) return 0;
  if (! GetCheckVect(p2, 0, 0, vect2) ) return 0;

  Double_t cx = (coord1[0]+coord2[0])/2;
  Double_t cy = (coord1[1]+coord2[1])/2;
  Double_t cz = (coord1[2]+coord2[2])/2;
  Double_t dx = coord2[0]-coord1[0];
  Double_t dy = coord2[1]-coord1[1];
  Double_t dz = coord2[2]-coord1[2];

  //=================================================
  // Positionning of the segment between the 2 points
  if (TMath::Abs(dy)<1e-231) dy = 1e-231;
  if (TMath::Abs(dz)<1e-231) dz = 1e-231;
  //Double_t angleRot1 = -TMath::ATan(dx/dy);
  //Double_t planDiagL = -TMath::Sqrt(dy*dy+dx*dx);
  //if (dy<0) planDiagL = -planDiagL;
  //Double_t angleRotDiag = TMath::ATan(planDiagL/dz);

  Double_t angleRot1    = -TMath::ATan2(dx,dy);
  Double_t planDiagL    =  TMath::Sqrt(dy*dy+dx*dx);
  Double_t angleRotDiag = -TMath::ATan2(planDiagL,dz);
  //--- (Calculate rotation of segment on the Z axis)
  //-- Here I'm trying to calculate the rotation to be applied in
  //-- order to match as closer as possible this segment and the 
  //-- previous one. 
  //-- It seems that some times it doesn't work ...
  Double_t angleRotZ = 0;
  TGeoRotation rotTemp("",angleRot1*TMath::RadToDeg(),
		       angleRotDiag*TMath::RadToDeg(), rotation);
  Double_t localX[3] = {0,1,0};
  Double_t globalX[3];
  rotTemp.LocalToMasterVect(localX, globalX);
  CopyFrom(localX, globalX);
  GetCheckVect(localX, p2Vol, 0, fgkCableMaxNodeLevel+1, globalX);
  Double_t orthVect[3];
  GetCheckVect(vect1, p2Vol, 0, fgkCableMaxNodeLevel+1, orthVect);
  if (p2>1) {
    Double_t orthVectNorm2 = ScalProd(orthVect,orthVect);
    Double_t alpha1 = ScalProd(fPreviousX,orthVect)/orthVectNorm2;
    Double_t alpha2 = ScalProd(globalX,orthVect)/orthVectNorm2;
    Double_t globalX1p[3], globalX2p[3];
    globalX1p[0] = fPreviousX[0] - alpha1*orthVect[0];
    globalX1p[1] = fPreviousX[1] - alpha1*orthVect[1];
    globalX1p[2] = fPreviousX[2] - alpha1*orthVect[2];
    globalX2p[0] = globalX[0] - alpha2*orthVect[0];
    globalX2p[1] = globalX[1] - alpha2*orthVect[1];
    globalX2p[2] = globalX[2] - alpha2*orthVect[2];
    //-- now I'm searching the 3th vect which makes an orthogonal base
    //-- with orthVect and globalX1p ...
    Double_t nulVect[3] = {0,0,0};
    Double_t axis3[3];
    TMath::Normal2Plane(nulVect, orthVect, globalX1p, axis3);
    Double_t globalX1pNorm2 = ScalProd(globalX1p, globalX1p);
    Double_t beta = ScalProd(globalX2p, globalX1p)/globalX1pNorm2;
    Double_t gamma = ScalProd(globalX2p, axis3);
    angleRotZ = (TMath::ATan2(1,0) - TMath::ATan2(beta, gamma))
                *TMath::RadToDeg();
  };
  //   cout << "!!!!!!!!!!!!!!!!!!!  angle = " <<angleRotZ << endl;
  CopyFrom(fPreviousX, globalX);
  //---
  Double_t localVect1[3], localVect2[3];
  TGeoRotation rot("",angleRot1*TMath::RadToDeg(),
		   angleRotDiag*TMath::RadToDeg(),
		   rotation);
// 		   rotation-angleRotZ);
// since angleRotZ doesn't always work, I won't use it ...

  rot.MasterToLocalVect(vect1, localVect1);
  rot.MasterToLocalVect(vect2, localVect2);

  //=================================================
  // Create the segment and add it to the mother volume
  TGeoVolume *vCableSegB = CreateBoxSegment(coord1, coord2);

  TGeoRotation rotArbSeg("", 0, 90, 0);
  rotArbSeg.MultiplyBy(&rot, kFALSE);
  TGeoTranslation trans("",cx, cy, cz);
  TGeoCombiTrans  *combiB = new TGeoCombiTrans(trans, rotArbSeg);
  p2Vol->AddNode(vCableSegB, p2, combiB);
  //=================================================;

  if (fDebug) {
    printf("---\n  Cable segment points : ");
    printf("%f, %f, %f\n",coord1[0], coord1[1], coord1[2]);
    printf("%f, %f, %f\n",coord2[0], coord2[1], coord2[2]);
  };

  return vCableSegB;
}

//________________________________________________________________________
TGeoVolume* AliITSv11GeomCableFlat::CreateAndInsertCableCylSegment(Int_t p2,
							  Double_t rotation)
{
  // Create a flat cable segment with a curvature between points p1 and p2.
  // The radius and position of the curve is defined by the
  // perpendicular vector of point p2 (the orientation of this vector
  // and the position of the 2 check points are enough to completely
  // define the curve)
  //    Rotation is the eventual rotation of the flat cable
  //    along its length axis
  //

  TGeoNode *mainNode;
  if (fInitialNode==0) {
    TObjArray *nodes = gGeoManager->GetListOfNodes();
    if (nodes->GetEntriesFast()==0) return 0;
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
    return 0;
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
      return 0;
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
			45-angleTorusSeg+rotation);
			//180-angleTorusSeg+rotation);
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
  TGeoVolume *vCableSegT = CreateCylSegment(torusPhi1, torusR);
  p2Vol->AddNode(vCableSegT, p2, combiTorus);

  if (fDebug) {
    printf("---\n  Cable segment points : ");
    printf("%f, %f, %f\n",coord1[0], coord1[1], coord1[2]);
    printf("%f, %f, %f\n",coord2[0], coord2[1], coord2[2]);
  };

  return vCableSegT;
}


//________________________________________________________________________
TGeoVolume *AliITSv11GeomCableFlat::CreateSegment( Double_t *coord1,
						   Double_t *coord2,
						   Double_t *localVect1,
						   Double_t *localVect2 )
{

  //=================================================
  // Calculate segment "deformation"
  Double_t dx = coord2[0]-coord1[0];
  Double_t dy = coord2[1]-coord1[1];
  Double_t dz = coord2[2]-coord1[2];
  Double_t length = TMath::Sqrt(dx*dx+dy*dy+dz*dz);

  Double_t cosTheta1 = -1./TMath::Sqrt( 1 + localVect1[0]*localVect1[0]
					/localVect1[2]/localVect1[2] );
  Double_t cosTheta2 = 1./TMath::Sqrt( 1 + localVect2[0]*localVect2[0]
					/localVect2[2]/localVect2[2] );
  if (localVect1[2]<0) cosTheta1 = -cosTheta1;
  if (localVect2[2]<0) cosTheta2 = -cosTheta2;

  Double_t dL1 = 0.5*fWidth*TMath::Tan(TMath::ACos(cosTheta1));
  Double_t dL2 = 0.5*fWidth*TMath::Tan(TMath::ACos(cosTheta2));
  if (localVect1[0]<0) dL1 = - dL1;
  if (localVect2[0]<0) dL2 = - dL2;
  //---
  Double_t cosPhi1 = -1./TMath::Sqrt( 1 + localVect1[1]*localVect1[1]
					/localVect1[2]/localVect1[2] );
  Double_t cosPhi2 = 1./TMath::Sqrt( 1 + localVect2[1]*localVect2[1]
					/localVect2[2]/localVect2[2] );
  if (localVect1[2]<0) cosPhi1 = -cosPhi1;
  if (localVect2[2]<0) cosPhi2 = -cosPhi2;

  Double_t tanACosCosPhi1 = TMath::Tan(TMath::ACos(cosPhi1));
  Double_t tanACosCosPhi2 = TMath::Tan(TMath::ACos(cosPhi2));
  if (localVect1[1]<0) tanACosCosPhi1 = -tanACosCosPhi1;
  if (localVect2[1]<0) tanACosCosPhi2 = -tanACosCosPhi2;

  Double_t dl1 = 0.5*fThick*tanACosCosPhi1;
  Double_t dl2 = 0.5*fThick*tanACosCosPhi2;

  //=================================================
  // Create the segment
  TGeoArb8 *cableSeg = new TGeoArb8(fThick/2);
  cableSeg->SetVertex( 0, -fWidth/2, -length/2 - dL1 + dl1);
  cableSeg->SetVertex( 1, -fWidth/2,  length/2 + dL2 - dl2);
  cableSeg->SetVertex( 2,  fWidth/2,  length/2 - dL2 - dl2);
  cableSeg->SetVertex( 3,  fWidth/2, -length/2 + dL1 + dl1);
  cableSeg->SetVertex( 4, -fWidth/2, -length/2 - dL1 - dl1);
  cableSeg->SetVertex( 5, -fWidth/2,  length/2 + dL2 + dl2);
  cableSeg->SetVertex( 6,  fWidth/2,  length/2 - dL2 + dl2);
  cableSeg->SetVertex( 7,  fWidth/2, -length/2 + dL1 - dl1);

  TGeoMedium *airSDD = gGeoManager->GetMedium("ITS_ITSair");
  TGeoVolume *vCableSeg = new TGeoVolume(GetName(), cableSeg, airSDD);

  // add all cable layers
  for (Int_t iLay=0; iLay<fNlayer; iLay++) {

    Double_t dl1Lay = 0.5*fLayThickness[iLay]*tanACosCosPhi1;
    Double_t dl2Lay = 0.5*fLayThickness[iLay]*tanACosCosPhi2;
 
    Double_t ztr = -fThick/2;
    for (Int_t i=0;i<iLay; i++) ztr+= fLayThickness[i];
    ztr+= fLayThickness[iLay]/2;

    Double_t dl1LayS = ztr*tanACosCosPhi1;
    Double_t dl2LayS = ztr*tanACosCosPhi2;

    TGeoArb8 *lay = new TGeoArb8(fLayThickness[iLay]/2);
    lay->SetVertex( 0, -fWidth/2, -length/2 - dL1 + dl1Lay - dl1LayS);
    lay->SetVertex( 1, -fWidth/2,  length/2 + dL2 - dl2Lay + dl2LayS);
    lay->SetVertex( 2,  fWidth/2,  length/2 - dL2 - dl2Lay + dl2LayS);
    lay->SetVertex( 3,  fWidth/2, -length/2 + dL1 + dl1Lay - dl1LayS);
    lay->SetVertex( 4, -fWidth/2, -length/2 - dL1 - dl1Lay - dl1LayS);
    lay->SetVertex( 5, -fWidth/2,  length/2 + dL2 + dl2Lay + dl2LayS);
    lay->SetVertex( 6,  fWidth/2,  length/2 - dL2 + dl2Lay + dl2LayS);
    lay->SetVertex( 7,  fWidth/2, -length/2 + dL1 - dl1Lay - dl1LayS);
    TGeoVolume *vLay = new TGeoVolume("vCableSegLay", lay, fLayMedia[iLay]);
    vLay->SetLineColor(fLayColor[iLay]);
    
    if (fTranslation[iLay]==0)
      fTranslation[iLay] = new TGeoTranslation(0, 0, ztr);
    vCableSeg->AddNode(vLay, iLay+1, fTranslation[iLay]);
  };

  vCableSeg->SetVisibility(kFALSE);
  return vCableSeg;
}


//________________________________________________________________________
TGeoVolume *AliITSv11GeomCableFlat::CreateCylSegment(Double_t &phi,
						     Double_t &r)
{

  Double_t phi1 = 360-phi;
  Double_t phi2 = 360+phi;

  Double_t rMin = r-fThick/2;
  Double_t rMax = r+fThick/2;
  //=================================================
  // Create the segment

  TGeoTubeSeg *cableSeg = new TGeoTubeSeg(rMin, rMax, fWidth/2,
					  phi1, phi2);
  TGeoMedium *airSDD = gGeoManager->GetMedium("ITS_ITSair");
  TGeoVolume *vCableSeg = new TGeoVolume(GetName(), cableSeg, airSDD);

  // add all cable layers
  for (Int_t iLay=0; iLay<fNlayer; iLay++) {
 
    Double_t ztr = -fThick/2;
    for (Int_t i=0;i<iLay; i++) ztr+= fLayThickness[i];

    rMin = r  + ztr;
    rMax = r  + ztr + fLayThickness[iLay];
    TGeoTubeSeg *lay = new TGeoTubeSeg(rMin, rMax, fWidth/2,
				       phi1, phi2);

    TGeoVolume *vLay = new TGeoVolume("vCableSegLay", lay, fLayMedia[iLay]);
    vLay->SetLineColor(fLayColor[iLay]);
    
    vCableSeg->AddNode(vLay, iLay+1, 0);
  };

  vCableSeg->SetVisibility(kFALSE);
  return vCableSeg;
}


//________________________________________________________________________
TGeoVolume *AliITSv11GeomCableFlat::CreateBoxSegment( Double_t *coord1,
						      Double_t *coord2)
{

  //=================================================
  // Create a segment for the case it is a simple box
  Double_t dx = coord2[0]-coord1[0];
  Double_t dy = coord2[1]-coord1[1];
  Double_t dz = coord2[2]-coord1[2];
  Double_t length = TMath::Sqrt(dx*dx+dy*dy+dz*dz);

  TGeoBBox *cableSeg = new  TGeoBBox(fWidth/2, length/2, fThick/2);

  TGeoMedium *airSDD = gGeoManager->GetMedium("ITS_ITSair");
  TGeoVolume *vCableSeg = new TGeoVolume(GetName(), cableSeg, airSDD);

  // add all cable layers
  for (Int_t iLay=0; iLay<fNlayer; iLay++) {
 
    Double_t ztr = -fThick/2;
    for (Int_t i=0;i<iLay; i++) ztr+= fLayThickness[i];
    ztr+= fLayThickness[iLay]/2;

    TGeoBBox *lay = new  TGeoBBox(fWidth/2, length/2, fLayThickness[iLay]/2);


    TGeoVolume *vLay = new TGeoVolume("vCableSegLay", lay, fLayMedia[iLay]);
    vLay->SetLineColor(fLayColor[iLay]);
    
    if (fTranslation[iLay]==0)
      fTranslation[iLay] = new TGeoTranslation(0, 0, ztr);
    vCableSeg->AddNode(vLay, iLay+1, fTranslation[iLay]);
  };

  vCableSeg->SetVisibility(kFALSE);
  return vCableSeg;
}

//________________________________________________________________________
void AliITSv11GeomCableFlat::SetNLayers(Int_t nLayers) {
  // Set the number of layers
  if((nLayers>0) &&(nLayers<=fgkCableMaxLayer)) {

    fNlayer = nLayers;
    for (Int_t i=0; i<fgkCableMaxLayer ; i++) {
      fLayThickness[i] = 0;
      fTranslation[i]  = 0;
      fLayColor[i]     = 0;
      fLayMedia[i]     = 0;  
    }; 
  };
}

//________________________________________________________________________
Int_t AliITSv11GeomCableFlat::SetLayer(Int_t nLayer, Double_t thick,
					  TGeoMedium *medium, Int_t color) {
  // Set the layer number nLayer
  if ((nLayer<0)||(nLayer>=fNlayer)) {
    printf("Set wrong layer number of the cable\n");
    return kFALSE;
  };
  if (nLayer>0)
    if (fLayThickness[nLayer-1]<=0) {
      printf("AliITSv11GeomCableFlat::SetLayer():"
	     " You must define cable layer %i first !",nLayer-1);
      return kFALSE;
    };

  Double_t thickTot = 0;
  for (Int_t i=0; i<nLayer; i++) thickTot += fLayThickness[i];
  thickTot += thick;
  if (thickTot-1e-10>fThick) {
    printf("Can't add this layer, cable thickness would be higher than total\n");
    return kFALSE;
  };

  fLayThickness[nLayer] = thick;
  fLayMedia[nLayer] = medium;
  fLayColor[nLayer] = color;
  fTranslation[nLayer]  = 0;
  return kTRUE;
}
