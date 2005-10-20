/**************************************************************************
 * Copyright(c) 1998-2004, ALICE Experiment at CERN, All rights reserved. *
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
// GeometryOfflineTrd1 class  for EMCAL : singleton
//*-- Author: Aleksei Pavlinov (WSU)

/* $Id$*/

#include <TBrowser.h>
#include <TMath.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TObjArray.h>

#include "AliEMCALGeometryOfflineTrd1.h"
#include "AliEMCALShishKebabTrd1Module.h"
#include "AliEMCALGeometry.h"

ClassImp(AliEMCALGeometryOfflineTrd1)

AliEMCALGeometryOfflineTrd1* AliEMCALGeometryOfflineTrd1::fgGeomOfflineTrd1=0;

AliEMCALGeometryOfflineTrd1* AliEMCALGeometryOfflineTrd1::GetInstance()
{
  if(fgGeomOfflineTrd1==0) {
    fgGeomOfflineTrd1 = new AliEMCALGeometryOfflineTrd1();
  }
  return fgGeomOfflineTrd1;
}

AliEMCALGeometryOfflineTrd1::AliEMCALGeometryOfflineTrd1() : TNamed("geomTRD1","")
{ // this private constarctor
  fGeometry = AliEMCALGeometry::GetInstance("SHISH_62_TRD1");
  Init();
}

void AliEMCALGeometryOfflineTrd1::Init()
{
  // Super module
  // ETA direction
  fMaxInEta = fGeometry->GetNZ();
  fTrd1Modules[0] =  new AliEMCALShishKebabTrd1Module();
  fSMMaxEta = 2*fMaxInEta;
  fSMPositionEta = new TVector2[fSMMaxEta];
  fSMPositionEta[0] = fTrd1Modules[0]->GetCenterOfCell(1);
  fSMPositionEta[1] = fTrd1Modules[0]->GetCenterOfCell(2);
  for(Int_t i=1; i<fMaxInEta; i++) {
    fTrd1Modules[i] = new AliEMCALShishKebabTrd1Module(*(fTrd1Modules[i-1]));
    fSMPositionEta[2*i]   = fTrd1Modules[i]->GetCenterOfCell(1);
    fSMPositionEta[2*i+1] = fTrd1Modules[i]->GetCenterOfCell(2);
  }
  // PHI direction
  fSMPositionPhi.Set(2*fGeometry->GetNPhi());
  fShiftOnPhi = -fGeometry->GetPhiModuleSize()*fGeometry->GetNPhi()/2;
  for(Int_t i=0; i<fGeometry->GetNPhi(); i++) {
    fSMPositionPhi[2*i]   = fGeometry->GetPhiModuleSize() * (double(i) + 0.25);
    fSMPositionPhi[2*i+1] = fGeometry->GetPhiTileSize()   + fSMPositionPhi[2*i];
  }
  
  // Super Module rotations
  fNPhiSuperModule = fGeometry->GetNPhiSuperModule(); // see AliEMCALv0
  double dphi = (fGeometry->GetArm1PhiMax() - fGeometry->GetArm1PhiMin()) / fNPhiSuperModule;
  double phi, phiRad;
  fSuperModuleRotationX.RotateX(TMath::Pi()); // matrix looks not so nice
  for(Int_t i=0; i<fNPhiSuperModule; i++){
    // rotations arround Z
    phi    = fGeometry->GetArm1PhiMin() + dphi*(2*i+1)/2.; // phi= 70, 90, 110, 130, 150, 170
    phiRad = phi*TMath::DegToRad();
    if(i==1) phiRad = TMath::PiOver2();
    fSuperModuleRotationZ[i].RotateZ(phiRad);
    TString ntmp("rotationZ_");
    ntmp += int(phi);
    fNameSuperModuleRotationZ[i] = ntmp;
    // Super Module rotation
    fSuperModuleRotation[2*i]   = fSuperModuleRotationZ[i]; // Z
    fSuperModuleRotation[2*i+1] = fSuperModuleRotationZ[i] * fSuperModuleRotationX; // Z*X
  }
  // Fill fXYZofCells
  fXYZofCells = new TObjArray(fGeometry->GetNCells());
  fXYZofCells->SetName("CellsInGC"); 
  Int_t nSupMod, nTower, nIphi, nIeta, iphi, ieta;
  for(Int_t absId=1; absId<=fGeometry->GetNCells(); absId++){
    if(fGeometry->GetCellIndex(absId, nSupMod,nTower,nIphi,nIeta)){
      fGeometry->GetCellPhiEtaIndexInSModule(nSupMod,nTower,nIphi,nIeta, iphi,ieta);
      TVector3 *v = new TVector3;
      v->SetX(fSMPositionEta[ieta-1].Y()); 
      v->SetZ(fSMPositionEta[ieta-1].X()); 
      v->SetY(fSMPositionPhi[iphi-1] + fShiftOnPhi);
      v->Transform(fSuperModuleRotation[nSupMod-1]);
      fXYZofCells->AddAt(v,absId-1);
    }
  }
}

TVector3& AliEMCALGeometryOfflineTrd1::PosInSuperModule(const int nSupMod, const Int_t nTower,const Int_t nIphi,const Int_t nIeta)
{ // 10-nov-04
  static Int_t iphi, ieta;
  static TVector3 v;
  fGeometry->GetCellPhiEtaIndexInSModule(nSupMod, nTower,nIphi,nIeta, iphi,ieta);

  // x-radius; y-phi; eta-z;
  v.SetXYZ(fSMPositionEta[ieta].Y(), fSMPositionPhi[iphi], fSMPositionEta[ieta].X());
  return v;
} 

void AliEMCALGeometryOfflineTrd1::PositionInSuperModule(const Int_t iphi, const Int_t ieta, 
double &lphi, double &leta)
{ 
  static Int_t ie=0;
  lphi = fSMPositionPhi[iphi-1];
  ie = ieta - 1;
  if(ie<0) ie = 0;
  if(ie>fSMMaxEta-1) ie = fSMMaxEta-1;
  leta = fSMPositionEta[ie].X();
}

void AliEMCALGeometryOfflineTrd1::PositionInSuperModule(const int nSupMod, const Int_t nTower, const Int_t nIphi, const Int_t nIeta,
double &lphi, double &leta)
{
  static Int_t iphi,ieta;
  fGeometry->GetCellPhiEtaIndexInSModule(nSupMod, nTower,nIphi,nIeta,iphi,ieta);
  PositionInSuperModule(iphi,ieta, lphi,leta);
}

TRotation* AliEMCALGeometryOfflineTrd1::Rotation(Int_t module)
{ // module chabge from 1 to 12
  if(module<1)  module=1;
  if(module>12) module=12;
  return &fSuperModuleRotation[module];
}

TVector3* AliEMCALGeometryOfflineTrd1::CellPosition(int absId)
{ // 15-nov-04
  if(absId<1 || absId>fXYZofCells->GetSize()) return 0;
  return (TVector3*)fXYZofCells->At(absId-1);
}

void AliEMCALGeometryOfflineTrd1::PrintSuperModule()
{ // 12-nov-04
  printf(" ** Super module ** fSMMaxEta %i fSMMaxPHI %i\n ETA     eta(Z)          (X)\n",
  fSMMaxEta,fSMPositionPhi.GetSize());
  for(Int_t i=0; i<fSMMaxEta; i++) {
    printf("%3i | %8.3f         %8.3f\n", i+1,fSMPositionEta[i].X(), fSMPositionEta[i].Y());
  }
  printf("\n PHI      (Y)\n");
  for(Int_t i=0; i<fSMPositionPhi.GetSize(); i++) {
    printf("%3i | %8.3f\n",i+1, fSMPositionPhi[i]);
  }
}

void AliEMCALGeometryOfflineTrd1::PrintCell(Int_t absId)
{
  Int_t nSupMod, nTower, nIphi, nIeta, iphi, ieta;
  if(fGeometry->GetCellIndex(absId, nSupMod,nTower,nIphi,nIeta)){
     fGeometry->GetCellPhiEtaIndexInSModule(nSupMod,nTower,nIphi,nIeta, iphi,ieta);
     TVector3 *v = CellPosition(absId);
     printf("(%5i) X %8.3f Y %8.3f Z %8.3f | #sup.Mod %2i #tower %3i nIphi %1i nIeta %1i | iphi %2i ieta %2i\n",
	    absId, v->X(),v->Y(),v->Z(), nSupMod,nTower,nIphi,nIeta, iphi,ieta);
  } else {
    Warning("PrintCell","Wrong abs id %i\n",absId); 
  }
}

void AliEMCALGeometryOfflineTrd1::Browse(TBrowser* b)
{
  if(fGeometry) b->Add(fGeometry);

  for(Int_t i=0; i<fMaxInEta; i++)  if(fTrd1Modules[i]>0) b->Add(fTrd1Modules[i]);

  for(Int_t i=0; i<fNPhiSuperModule; i++){ 
    b->Add(&fSuperModuleRotationZ[i], fNameSuperModuleRotationZ[i].Data());
    for(Int_t j=0; j<2; j++) {
      TString name("rotationM_"); name += (2*i+j);
      b->Add(&fSuperModuleRotation[2*i+j], name.Data());
    }
  }

  b->Add(&fSuperModuleRotationX, "rotationX_180");

  b->Add(fXYZofCells);
}

Bool_t AliEMCALGeometryOfflineTrd1::IsFolder() const
{
  return kTRUE;
}
