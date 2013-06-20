/**************************************************************************
 * Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//
// class used to extract and store reco info of calo cluster
//
// Author: H-S. Zhu, hongsheng.zhu@cern.ch
//                   hszhu@iopp.ccnu.edu.cn
/////////////////////////////////////////////////////////////

#include <iostream>

#include <TF1.h>
#include <TVector3.h>

#include "AliESDEvent.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliVCluster.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliPHOSGeoUtils.h"
#include "AliPHOSEsdCluster.h"
#include "AliPHOSCalibData.h"
#include "AliCaloClusterInfo.h"

class TObject;

ClassImp(AliCaloClusterInfo) 

Float_t AliCaloClusterInfo::fgLogWeight = 0.07;

//-----------------------------------------------------------------------------
Double_t NonLinear(Double_t *x, Double_t * /*par*/)
{ 
  // a = par[0], b = par[1].
  // 1+a*exp(-e/b)
  
// return 0.0241+1.0504*x[0]+0.000249*x[0]*x[0];
 return 1.015*(0.0241+1.0504*x[0]+0.000249*x[0]*x[0]);
  
}

//-------------------------------------------------------------------------------------------
AliCaloClusterInfo::AliCaloClusterInfo():
  TObject(), fLorentzVector(), fModule(0), fNCells(0), fTRUNumber(0), fNTracksMatched(0), fTrackCharge(0), fPIDBit(0x0),
  fDistToBad(0.), fEmcCpvDistance(0.), fM02(0.), fM20(0.), fTOF(0.), fTrackDz(0.), fTrackDx(0.), fTrackPt(0.)
{
  //
  // Default constructor
  //
} 

//-------------------------------------------------------------------------------------------
AliCaloClusterInfo::AliCaloClusterInfo(AliVCluster* const clust, AliESDEvent* const esd, AliPHOSGeoUtils* const phosGeo, Double_t vtx[3]):
  TObject(), fLorentzVector(), fModule(0), fNCells(0), fTRUNumber(0), fNTracksMatched(0), fTrackCharge(0), fPIDBit(0x0),
  fDistToBad(0.), fEmcCpvDistance(0.), fM02(0.), fM20(0.), fTOF(0.), fTrackDz(0.), fTrackDx(0.), fTrackPt(0.)
{
  //
  // constructor
  //
  this->FillCaloClusterInfo(clust, esd, phosGeo, vtx);
} 

//-------------------------------------------------------------------------------------------
AliCaloClusterInfo::AliCaloClusterInfo(const AliCaloClusterInfo &src):
  TObject(src), fLorentzVector(src.fLorentzVector), fModule(src.fModule), fNCells(src.fNCells), fTRUNumber(src.fTRUNumber),
  fNTracksMatched(src.fNTracksMatched), fTrackCharge(src.fTrackCharge), fPIDBit(src.fPIDBit), fDistToBad(src.fDistToBad),
  fEmcCpvDistance(src.fEmcCpvDistance), fM02(src.fM02), fM20(src.fM20), fTOF(src.fTOF), fTrackDz(src.fTrackDz),
  fTrackDx(src.fTrackDx), fTrackPt(src.fTrackPt)
{
  //
  // copy constructor
  //
} 

//-------------------------------------------------------------------------------------------
AliCaloClusterInfo& AliCaloClusterInfo::operator=(const AliCaloClusterInfo &src)
{
  //
  // assignment constructor
  //
  if(&src==this) return *this;

  fLorentzVector   = src.fLorentzVector;
  fModule          = src.fModule;
  fNCells          = src.fNCells;
  fTRUNumber       = src.fTRUNumber;
  fNTracksMatched  = src.fNTracksMatched;
  fTrackCharge     = src.fTrackCharge;
  fPIDBit          = src.fPIDBit;
  fDistToBad       = src.fDistToBad;
  fEmcCpvDistance  = src.fEmcCpvDistance;
  fM02             = src.fM02;
  fM20             = src.fM20;
  fTOF             = src.fTOF;
  fTrackDz         = src.fTrackDz;
  fTrackDx         = src.fTrackDx;
  fTrackPt         = src.fTrackPt;

  return *this;
} 

//-------------------------------------------------------------------------------------------
AliCaloClusterInfo::~AliCaloClusterInfo()
{
  //
  // destructor
  //
}

//-----------------------------------------------------------------------------
void AliCaloClusterInfo::FillCaloClusterInfo(AliVCluster* const clust, AliESDEvent* const esd, AliPHOSGeoUtils* const phosGeo, Double_t vtx[3])
{
  // extract information of calo clusters
  AliAODTrack *trkAOD = 0x0;
  AliESDtrack *trkESD = 0x0;

  Int_t   relId[4]    = {0,0,0,0}; // module = relId[0]; cellX = relId[2]; cellZ = relId[3];
  Float_t position[3] = {0,0,0};
  clust->GetPosition(position); TVector3 global(position);
  phosGeo->GlobalPos2RelId(global,relId);
  fModule         = relId[0];
  fTRUNumber      = GetTRUNumber(relId[2], relId[3]);

  if (esd) { // TODO recalibration for ESD
    TVector3 vtxVector(vtx);
    AliPHOSCalibData *calibData = new AliPHOSCalibData();
    TF1 *nonLinCorr = new TF1("Non-linear", NonLinear, 0., 40., 0);

    AliPHOSEsdCluster phosClust( *(AliESDCaloCluster*) (clust) );
    phosClust.Recalibrate(calibData, esd->GetPHOSCells());
    Reclusterize(&phosClust, phosGeo);               // reclusterize 
    phosClust.EvalAll(fgLogWeight, vtxVector);       // recalculate all cluster parameters
    phosClust.SetE(nonLinCorr->Eval(phosClust.E())); // non-linear correction

    TVector3 localPos;
    const Float_t shiftX[6] = { 0.,-2.3,-2.11,-1.53, 0., 0. };
    const Float_t shiftZ[6] = { 0.,-0.4, 0.52, 0.8,  0., 0. };
    phosGeo->Global2Local(localPos, global, fModule);
    phosGeo->Local2Global(fModule, localPos.X()+shiftX[fModule], localPos.Z()+shiftZ[fModule], global);
    position[0] = global.X();
    position[1] = global.Y();
    position[2] = global.Z();
    phosClust.SetPosition(position);

    phosClust.GetMomentum(fLorentzVector, vtx);

    //TODO: Check, this may be LHC10h specific:
    if (fModule==2) fLorentzVector *= 135.5/134.0;
    if (fModule==3) fLorentzVector *= 135.5/137.2;

    Int_t iESDtrack = clust->GetTrackMatchedIndex();
    if (iESDtrack>-1) { 
      trkESD = esd->GetTrack(iESDtrack);
      if (trkESD) {
        fTrackPt     = trkESD->Pt();
        fTrackCharge = trkESD->Charge();
      }
    }
    else { fTrackPt = -1.; fTrackCharge = 0; }
  }else {
    clust->GetMomentum(fLorentzVector, vtx);

    if (clust->GetNTracksMatched()>0) {
      trkAOD = dynamic_cast <AliAODTrack*> (clust->GetTrackMatched(0));
      if (trkAOD) {
        fTrackPt     = trkAOD->Pt();
        fTrackCharge = trkAOD->Charge();
      }
    }
    else { fTrackPt = -1.; fTrackCharge = 0; }
  }

  fNCells         = clust->GetNCells();
  fNTracksMatched = clust->GetNTracksMatched();
  fDistToBad      = clust->GetDistanceToBadChannel();
  fEmcCpvDistance = clust->GetEmcCpvDistance();
  fM02            = clust->GetM02();
  fM20            = clust->GetM20();
  fTOF            = clust->GetTOF();
  fTrackDz        = clust->GetTrackDz();
  fTrackDx        = clust->GetTrackDx();
  if (TestDisp())                             fPIDBit |= BIT(1); // Disp
  if (IsInFiducialRegion(relId[2], relId[3])) fPIDBit |= BIT(2); // Fiducial

  return;
}

//-----------------------------------------------------------------------------
void AliCaloClusterInfo::Reclusterize(AliVCluster *clust, AliPHOSGeoUtils* const phosGeo)
{
  // Reclusterize to have continuous cluster

  const Int_t oldMulDigit = clust->GetNCells();
  Double32_t       *elist = clust->GetCellsAmplitudeFraction();
  UShort_t         *dlist = clust->GetCellsAbsId();

  Int_t index[oldMulDigit];
  Bool_t used[oldMulDigit];
  for (Int_t i=0; i<oldMulDigit; i++) used[i] = 0;
  Int_t   inClu = 0;
  Double_t eMax = 0.;
  for (Int_t iDigit=0; iDigit<oldMulDigit; iDigit++) { // find maximum
    if (eMax<elist[iDigit]) {
      eMax     = elist[iDigit];
      index[0] = iDigit;
      inClu    = 1;
    }
  }
  if (inClu==0) return;    // empty cluster
  used[index[0]] = kTRUE ; // mark as used
  for (Int_t i=0; i<inClu; i++) {
    for (Int_t iDigit=0; iDigit<oldMulDigit; iDigit++) {
      if (used[iDigit]) continue; // already used
      if (AreNeighbors(dlist[index[i]], dlist[iDigit], phosGeo)) {
        index[inClu] = iDigit;
        inClu++;
        used[iDigit] = kTRUE;
      }
    }
  }

  if (inClu==oldMulDigit) return; // no need to modify

  // copy
  clust->SetNCells(inClu);
  UShort_t tmpD[oldMulDigit];
  Double_t tmpE[oldMulDigit];
  for (Int_t i=0; i<oldMulDigit; i++) {
    tmpD[i]=dlist[i];
    tmpE[i]=elist[i];
  }
  // change order of digits in list so that
  // first inClu cells were true ones
  for (Int_t i=0; i<inClu; i++) {
    dlist[i]=tmpD[index[i]];
    elist[i]=tmpE[index[i]];
  }

  return;
}

//-----------------------------------------------------------------------------
Bool_t AliCaloClusterInfo::TestCPV(Double_t mf)
{
  // Parameterization of LHC10h period

  Double_t meanX=0;
  Double_t meanZ=0.;
  Double_t sx = TMath::Min(5.4,2.59719e+02*TMath::Exp(-fTrackPt/1.02053e-01)+
                           6.58365e-01*5.91917e-01*5.91917e-01/((fTrackPt-9.61306e-01)*(fTrackPt-9.61306e-01)+5.91917e-01*5.91917e-01)+1.59219);
  Double_t sz = TMath::Min(2.75,4.90341e+02*1.91456e-02*1.91456e-02/(fTrackPt*fTrackPt+1.91456e-02*1.91456e-02)+1.60);

  if(mf<0.) { // Field --
    meanZ = -0.468318;
    if (fTrackCharge>0)
      meanX = TMath::Min(7.3, 3.89994*1.20679*1.20679/(fTrackPt*fTrackPt+1.20679*1.20679)+0.249029+2.49088e+07*TMath::Exp(-fTrackPt*3.33650e+01));
    else
      meanX =-TMath::Min(7.7,3.86040*0.912499*0.912499/(fTrackPt*fTrackPt+0.912499*0.912499)+1.23114+4.48277e+05*TMath::Exp(-fTrackPt*2.57070e+01));
  }
  else {      // Field ++
    meanZ = -0.468318;
    if (fTrackCharge>0)
      meanX =-TMath::Min(8.0,3.86040*1.31357*1.31357/(fTrackPt*fTrackPt+1.31357*1.31357)+0.880579+7.56199e+06*TMath::Exp(-fTrackPt*3.08451e+01));
    else
      meanX = TMath::Min(6.85, 3.89994*1.16240*1.16240/(fTrackPt*fTrackPt+1.16240*1.16240)-0.120787+2.20275e+05*TMath::Exp(-fTrackPt*2.40913e+01));
  }

  Double_t rx = (fTrackDx-meanX)/sx;
  Double_t rz = (fTrackDz-meanZ)/sz;

  return TMath::Sqrt(rx*rx+rz*rz) > 2.;
}

//-----------------------------------------------------------------------------
Bool_t AliCaloClusterInfo::TestDisp()
{
  Double_t pt = fLorentzVector.E();
  Double_t l0 = fM02;
  Double_t l1 = fM20;

  Double_t l0Mean  = 1.53126+9.50835e+06/(1.+1.08728e+07*pt+1.73420e+06*pt*pt);
  Double_t l1Mean  = 1.12365+0.123770*TMath::Exp(-pt*0.246551)+5.30000e-03*pt;
  Double_t l0Sigma = 6.48260e-02+7.60261e+10/(1.+1.53012e+11*pt+5.01265e+05*pt*pt)+9.00000e-03*pt;
  Double_t l1Sigma = 4.44719e-04+6.99839e-01/(1.+1.22497e+00*pt+6.78604e-07*pt*pt)+9.00000e-03*pt;
  Double_t c       =-0.35-0.550*TMath::Exp(-0.390730*pt);
  Double_t R2      = 0.5*(l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma +
                     0.5*(l0-l0Mean)*(l0-l0Mean)/l0Sigma/l0Sigma +
                     0.5*c*(l1-l1Mean)*(l0-l0Mean)/l1Sigma/l0Sigma;

  return R2 < 2.5*2.5;
}

//-----------------------------------------------------------------------------
Bool_t AliCaloClusterInfo::IsInFiducialRegion(Int_t cellX, Int_t cellZ)
{
  const Int_t edgeX = 2;
  const Int_t edgeZ = 2;
  if (cellX >edgeX && cellX<(65-edgeX) && cellZ>edgeZ && cellZ <(57-edgeZ)) return kTRUE;
  return kFALSE;

//Double_t eta = TMath::Abs(fLorentzVector.Eta());                                // abs eta
//Double_t phi = (fLorentzVector.Phi())*TMath::RadToDeg(); if (phi<0.) phi+=360.; // phi in degree

//return (Bool_t)(eta<0.13 /*&& ((phi>260.25 && phi<279.75) || (phi>300.25 && phi<319.75))*/);
}

//-----------------------------------------------------------------------------
Bool_t AliCaloClusterInfo::AreNeighbors(Int_t id1, Int_t id2, AliPHOSGeoUtils* const phosGeo)
{
  // return true if absId are "Neighbors" (adjacent, including diagornaly,)

  Int_t relid1[4];   phosGeo->AbsToRelNumbering(id1, relid1);
  Int_t relid2[4];   phosGeo->AbsToRelNumbering(id2, relid2);

  // if inside the same PHOS module
  if ((relid1[0] == relid2[0]) && (relid1[1]==relid2[1])) {
    const Int_t rowdiff = TMath::Abs(relid1[2]-relid2[2]);
    const Int_t coldiff = TMath::Abs(relid1[3]-relid2[3]);

    // and if diff in both direction is 1 or less
    if ((coldiff<2)  && (rowdiff<2)) return kTRUE; // are neighbors
  }

  return kFALSE;
}

//-----------------------------------------------------------------------------
Int_t AliCaloClusterInfo::GetTRUNumber(Int_t cellX, Int_t cellZ)
{
  // Return TRU region number for given cell.
  // cellX: [1-64], cellZ: [1-56]
  
  // RCU0: TRU 1,2
  if (0<cellX && cellX<17) {
    if (0<cellZ && cellZ<29) return 2;
    else                     return 1;
  }

  // RCU1: TRU 3,4
  if (16<cellX && cellX<33) {
    if (0<cellZ && cellZ<29) return 4;
    else                     return 3;
  }

  // RCU2: TRU 5,6
  if (32<cellX && cellX<49) {
    if (0<cellZ && cellZ<29) return 6;
    else                     return 5;
  }

  // RCU3: TRU 7,8
  if (48<cellX && cellX<65) {
    if (0<cellZ && cellZ<29) return 8;
    else                     return 7;
  }

  return -1;
}
