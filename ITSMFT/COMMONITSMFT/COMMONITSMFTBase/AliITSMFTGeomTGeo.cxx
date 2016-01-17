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

///////////////////////////////////////////////////////////////////////////
//    AliITSMFTGeomTGeo is a simple interface class to TGeoManager         //
//    It is used in the simulation and reconstruction in order to        //
//    query the TGeo geometry                                            //
//                                                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TClass.h>
#include <TString.h>
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>
#include <TGeoShape.h>
#include <TGeoBBox.h>
#include <TDatime.h>
#include <TMath.h>
#include <TSystem.h>

#include "AliITSMFTGeomTGeo.h"
#include "AliLog.h"
#include "AliAlignObj.h"
#include "AliITSMFTSegmentationPix.h"
using namespace TMath;

ClassImp(AliITSMFTGeomTGeo)

//______________________________________________________________________
AliITSMFTGeomTGeo::AliITSMFTGeomTGeo()
  :fNChips(0)
  ,fMatSens(0)
  ,fMatT2L(0)
  ,fSegm(0)
{
  // default c-tor
}

//______________________________________________________________________
AliITSMFTGeomTGeo::AliITSMFTGeomTGeo(const AliITSMFTGeomTGeo &src)
  :TObject(src)
  ,fNChips(src.fNChips)
  ,fMatSens(0)
  ,fMatT2L(0)
  ,fSegm(0)
{
  // copy c-tor
    if (src.fMatSens) {
      fMatSens = new TObjArray(fNChips);
      fMatSens->SetOwner(kTRUE);
      for (int i=0;i<fNChips;i++) {
	const TGeoHMatrix* mat = (TGeoHMatrix*)src.fMatSens->At(i);
	fMatSens->AddAt(new TGeoHMatrix(*mat),i);
      }
    }
    if (src.fMatT2L) {
      fMatT2L = new TObjArray(fNChips);
      fMatT2L->SetOwner(kTRUE);
      for (int i=0;i<fNChips;i++) {
	const TGeoHMatrix* mat =(TGeoHMatrix*) src.fMatT2L->At(i);
	fMatT2L->AddAt(new TGeoHMatrix(*mat),i);
      }
    }
    if (src.fSegm) {
      int sz = src.fSegm->GetEntriesFast();
      fSegm = new TObjArray(sz);
      fSegm->SetOwner(kTRUE);
      for (int i=0;i<sz;i++) {
	AliITSMFTSegmentationPix* sg = (AliITSMFTSegmentationPix*)src.fSegm->UncheckedAt(i);
	if (!sg) continue;
	fSegm->AddAt(sg->Clone(),i);
      }
    }
}

//______________________________________________________________________
AliITSMFTGeomTGeo::~AliITSMFTGeomTGeo()
{
  //d-tor
  delete fMatT2L;
  delete fMatSens;
  delete fSegm;
}

//______________________________________________________________________
TGeoHMatrix* AliITSMFTGeomTGeo::GetMatrix(Int_t index)  const
{
  // Get the transformation matrix for a given chip 'index'
  // by quering the TGeoManager
  static TGeoHMatrix matTmp;
  TGeoPNEntry *pne = GetPNEntry(index);
  if (!pne) return NULL;

  TGeoPhysicalNode *pnode = pne->GetPhysicalNode();
  if (pnode) return pnode->GetMatrix();

  const char* path = pne->GetTitle();
  gGeoManager->PushPath(); // Preserve the modeler state.
  if (!gGeoManager->cd(path)) {
    gGeoManager->PopPath();
    AliError(Form("Volume path %s not valid!",path));
    return NULL;
  }
  matTmp = *gGeoManager->GetCurrentMatrix();
  gGeoManager->PopPath();
  return &matTmp;
}

//______________________________________________________________________
Bool_t AliITSMFTGeomTGeo::GetTranslation(Int_t index, Double_t t[3])  const
{
  // Get the translation vector for a given chip 'index'
  // by quering the TGeoManager
  TGeoHMatrix *m = GetMatrix(index);
  if (!m) return kFALSE;

  Double_t *trans = m->GetTranslation();
  for (Int_t i = 0; i < 3; i++) t[i] = trans[i];

  return kTRUE;
}

//______________________________________________________________________
Bool_t AliITSMFTGeomTGeo::GetRotation(Int_t index, Double_t r[9])  const
{
  // Get the rotation matrix for a given chip 'index'
  // by quering the TGeoManager
  TGeoHMatrix *m = GetMatrix(index);
  if (!m) return kFALSE;

  Double_t *rot = m->GetRotationMatrix();
  for (Int_t i = 0; i < 9; i++) r[i] = rot[i];

  return kTRUE;
}

//______________________________________________________________________
Bool_t AliITSMFTGeomTGeo::GetOrigTranslation(Int_t index, Double_t t[3])  const
{
  // Get the original translation vector (ideal geometry)
  // for a given chip 'index' by quering the TGeoManager
  TGeoHMatrix m;
  if (!GetOrigMatrix(index,m)) return kFALSE;

  Double_t *trans = m.GetTranslation();
  for (Int_t i = 0; i < 3; i++) t[i] = trans[i];

  return kTRUE;
}

//______________________________________________________________________
Bool_t AliITSMFTGeomTGeo::GetOrigRotation(Int_t index, Double_t r[9])  const
{
  // Get the original rotation matrix (ideal geometry)
  // for a given chip 'index' by quering the TGeoManager
  TGeoHMatrix m;
  if (!GetOrigMatrix(index,m)) return kFALSE;

  Double_t *rot = m.GetRotationMatrix();
  for (Int_t i = 0; i < 9; i++) r[i] = rot[i];

  return kTRUE;
}

//______________________________________________________________________
TGeoHMatrix* AliITSMFTGeomTGeo::ExtractMatrixT2L(Int_t index) const
{
  // Get the matrix which transforms from the tracking to local r.s.
  // The method queries directly the TGeoPNEntry
  TGeoPNEntry *pne = GetPNEntry(index);
  if (!pne) return NULL;

  TGeoHMatrix *m = (TGeoHMatrix*) pne->GetMatrix();
  if (!m) AliError(Form("TGeoPNEntry (%s) contains no matrix !",pne->GetName()));

  return m;
}

//______________________________________________________________________
Bool_t AliITSMFTGeomTGeo::GetTrackingMatrix(Int_t index, TGeoHMatrix &m)
{
  // Get the matrix which transforms from the tracking r.s. to
  // the global one.
  // Returns kFALSE in case of error.
  m.Clear();

  TGeoHMatrix *m1 = GetMatrix(index);
  if (!m1) return kFALSE;

  const TGeoHMatrix *m2 = GetMatrixT2L(index);
  if (!m2) return kFALSE;

  m = *m1;
  m.Multiply(m2);

  return kTRUE;
}

//______________________________________________________________________
void AliITSMFTGeomTGeo::FetchMatrices()
{
  // store pointer on often used matrices for faster access
  if (!gGeoManager) AliFatal("Geometry is not loaded");
  fMatSens = new TObjArray(fNChips);
  fMatSens->SetOwner(kTRUE);
  for (int i=0;i<fNChips;i++) fMatSens->AddAt(new TGeoHMatrix(*ExtractMatrixSens(i)),i);
  CreateT2LMatrices();
}

//______________________________________________________________________
void AliITSMFTGeomTGeo::CreateT2LMatrices()
{
  // create tracking to local (Sensor!) matrices
  fMatT2L  = new TObjArray(fNChips);  
  fMatT2L->SetOwner(kTRUE);
  TGeoHMatrix matLtoT;
  double locA[3]={-100,0,0},locB[3]={100,0,0},gloA[3],gloB[3];
  for (int isn=0;isn<fNChips;isn++) {
    const TGeoHMatrix* matSens = GetMatrixSens(isn);
    if (!matSens) {AliFatal(Form("Failed to get matrix for sensor %d",isn)); return;}
    matSens->LocalToMaster(locA,gloA);
    matSens->LocalToMaster(locB,gloB);
    double dx = gloB[0]-gloA[0];
    double dy = gloB[1]-gloA[1];
    double t = (gloB[0]*dx+gloB[1]*dy)/(dx*dx+dy*dy),x=gloB[0]-dx*t,y=gloB[1]-dy*t;
    TGeoHMatrix* t2l = new TGeoHMatrix();
    t2l->RotateZ(ATan2(y,x)*RadToDeg()); // rotate in direction of normal to the sensor plane
    t2l->SetDx(x);
    t2l->SetDy(y);
    t2l->MultiplyLeft(&matSens->Inverse());
    fMatT2L->AddAt(t2l,isn);
    /*
    const double *gtrans = matSens->GetTranslation();
    memcpy(&rotMatrix[0], matSens->GetRotationMatrix(), 9*sizeof(Double_t));
    Double_t al = -ATan2(rotMatrix[1],rotMatrix[0]);
    Double_t rSens = Sqrt(gtrans[0]*gtrans[0] + gtrans[1]*gtrans[1]);
    Double_t tanAl = ATan2(gtrans[1],gtrans[0]) - Pi()/2; //angle of tangent
    Double_t alTr = tanAl - al;
    //
    // The X axis of tracking frame must always look outward
    loc[1] = rSens/2;
    matSens->LocalToMaster(loc,glo);
    double rPos = Sqrt(glo[0]*glo[0] + glo[1]*glo[1]);
    Bool_t rotOutward = rPos>rSens ? kFALSE : kTRUE;
    //
    // Transformation matrix
    matLtoT.Clear();
    matLtoT.SetDx(-rSens*Sin(alTr)); // translation
    matLtoT.SetDy(0.);
    matLtoT.SetDz(gtrans[2]);
    // Rotation matrix
    rotMatrix[0]= 0;  rotMatrix[1]= 1;  rotMatrix[2]= 0; // + rotation
    rotMatrix[3]=-1;  rotMatrix[4]= 0;  rotMatrix[5]= 0;
    rotMatrix[6]= 0;  rotMatrix[7]= 0;  rotMatrix[8]= 1;
    //
    TGeoRotation rot;
    rot.SetMatrix(rotMatrix);
    matLtoT.MultiplyLeft(&rot);
    if (rotOutward) matLtoT.RotateZ(180.);
    // Inverse transformation Matrix
    fMatT2L->AddAt(new TGeoHMatrix(matLtoT.Inverse()),isn);
    */
  }
  //
}

