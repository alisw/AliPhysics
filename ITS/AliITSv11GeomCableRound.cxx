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
#include <TGeoMatrix.h>

#include "AliITSv11GeomCableRound.h"

//*************************************************************************
//   Class for round cables
//
// Ludovic Gaudichet                                   gaudichet@to.infn.it
//*************************************************************************

ClassImp(AliITSv11GeomCableRound)

//________________________________________________________________________
AliITSv11GeomCableRound::
AliITSv11GeomCableRound(const char* name, Double_t radius) :
AliITSv11GeomCable(name) {
  // Constructor
  fRadius = radius;
  fNlayer = 0;
  for (Int_t i=0; i<fgkCableMaxLayer ; i++) {
    fLayThickness[i] = 0;
    fLayColor[i] = 0;
    fLayMedia[i] = 0;
  };
  fPhiMin = 0;
  fPhiMax = 360;
};

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
};

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
};
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
};

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

};

//________________________________________________________________________
Int_t AliITSv11GeomCableRound::CreateAndInsertCableSegment(Int_t p2)
{
//    Creates a cable segment between points p1 and p2.
//    Rotation is the eventual rotation of the flat cable
//    along its length axis
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
					localVect1, localVect2);

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
//   TGeoMedium *airSDD = gGeoManager->GetMedium("ITSsddAir");
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
};

//________________________________________________________________________
TGeoVolume *AliITSv11GeomCableRound::CreateSegment( Double_t *coord1,
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
  TGeoCtub *cableSeg = new TGeoCtub(0., fRadius, length/2, fPhiMin, fPhiMax,
				    localVect1[0],localVect1[1],localVect1[2],
				    localVect2[0],localVect2[1],localVect2[2]);

  TGeoMedium *airSDD = gGeoManager->GetMedium("ITSair");
  TGeoVolume *vCableSeg = new TGeoVolume(GetName(), cableSeg, airSDD);

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
};


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
};

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
};
