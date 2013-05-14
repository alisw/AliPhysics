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

#include "AliESDEvent.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliVCluster.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliCaloClusterInfo.h"

class TObject;

ClassImp(AliCaloClusterInfo) 

//-------------------------------------------------------------------------------------------
AliCaloClusterInfo::AliCaloClusterInfo():
  TObject(), fLorentzVector(), fModule(0), fTRUNumber(0), fNCells(0), fPIDBit(0x0),
  fDistToBad(0.), fEmcCpvDistance(0.), fM02(0.), fM20(0.), fTOF(0.)
{
  //
  // Default constructor
  //
} 

//-------------------------------------------------------------------------------------------
AliCaloClusterInfo::AliCaloClusterInfo(AliVCluster* const clust, AliESDEvent* const esd, Int_t relID[4], Double_t mf):
  TObject(), fLorentzVector(), fModule(0), fTRUNumber(0), fNCells(0), fPIDBit(0x0),
  fDistToBad(0.), fEmcCpvDistance(0.), fM02(0.), fM20(0.), fTOF(0.)
{
  //
  // constructor
  //
  this->FillCaloClusterInfo(clust, esd, relID, mf);
} 

//-------------------------------------------------------------------------------------------
AliCaloClusterInfo::AliCaloClusterInfo(const AliCaloClusterInfo &src):
  TObject(src), fLorentzVector(src.fLorentzVector), fModule(src.fModule), fTRUNumber(src.fTRUNumber),fNCells(src.fNCells), fPIDBit(src.fPIDBit),
  fDistToBad(src.fDistToBad), fEmcCpvDistance(src.fEmcCpvDistance), fM02(src.fM02), fM20(src.fM20), fTOF(src.fTOF)
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
  fTRUNumber       = src.fTRUNumber;
  fNCells          = src.fNCells;
  fPIDBit          = src.fPIDBit;
  fDistToBad       = src.fDistToBad;
  fEmcCpvDistance  = src.fEmcCpvDistance;
  fM02             = src.fM02;
  fM20             = src.fM20;
  fTOF             = src.fTOF;

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
void AliCaloClusterInfo::FillCaloClusterInfo(AliVCluster* const clust, AliESDEvent* const esd, Int_t relID[4], Double_t mf)
{
  // extract information of calo clusters

  Short_t  trkCharge = 0;
  Double_t trkDz = 0., trkDx = 0., trkPt = -1.;

  if (!esd) { // for AOD
    AliAODTrack *trkAOD = 0x0;
    if (clust->GetNTracksMatched()>0) {
      trkAOD = dynamic_cast <AliAODTrack*> (clust->GetTrackMatched(0));
      if (trkAOD) {
        trkPt     = trkAOD->Pt();
        trkCharge = trkAOD->Charge();
      }
    }
  } else { // for ESD
    Int_t iESDtrack = clust->GetTrackMatchedIndex();
    AliESDtrack *trkESD = 0x0;
    if (iESDtrack>-1) { 
      trkESD = esd->GetTrack(iESDtrack);
      if (trkESD) {
        trkPt     = trkESD->Pt();
        trkCharge = trkESD->Charge();
      }
    }
  }
  trkDz = clust->GetTrackDz();
  trkDx = clust->GetTrackDx();

  fModule         = relID[0];
  fTRUNumber      = GetTRUNumber(relID[2], relID[3]);
  fNCells         = clust->GetNCells();
  fDistToBad      = clust->GetDistanceToBadChannel();
  fEmcCpvDistance = clust->GetEmcCpvDistance();
  fM02            = clust->GetM02();
  fM20            = clust->GetM20();
  fTOF            = clust->GetTOF();

  Double_t cpv = TestCpv(trkPt, trkCharge, trkDz, trkDx, mf);
  if (trkPt == -1. || cpv>2.) this->SetPIDBit(BIT(0));
  if (trkPt == -1. || cpv>4.) this->SetPIDBit(BIT(2));

  return;
}

//-----------------------------------------------------------------------------
Double_t AliCaloClusterInfo::TestCpv(Double_t trkPt, Short_t trkCharge, Double_t trkDz, Double_t trkDx, Double_t mf)
{
  // Parameterization of LHC10h period

  Double_t meanX=0;
  Double_t meanZ=0.;
  Double_t sx = TMath::Min(5.4,2.59719e+02*TMath::Exp(-trkPt/1.02053e-01)+
                           6.58365e-01*5.91917e-01*5.91917e-01/((trkPt-9.61306e-01)*(trkPt-9.61306e-01)+5.91917e-01*5.91917e-01)+1.59219);
  Double_t sz = TMath::Min(2.75,4.90341e+02*1.91456e-02*1.91456e-02/(trkPt*trkPt+1.91456e-02*1.91456e-02)+1.60);

  if(mf<0.) { // Field --
    meanZ = -0.468318;
    if (trkCharge>0)
      meanX = TMath::Min(7.3, 3.89994*1.20679*1.20679/(trkPt*trkPt+1.20679*1.20679)+0.249029+2.49088e+07*TMath::Exp(-trkPt*3.33650e+01));
    else
      meanX =-TMath::Min(7.7,3.86040*0.912499*0.912499/(trkPt*trkPt+0.912499*0.912499)+1.23114+4.48277e+05*TMath::Exp(-trkPt*2.57070e+01));
  }
  else {      // Field ++
    meanZ = -0.468318;
    if (trkCharge>0)
      meanX =-TMath::Min(8.0,3.86040*1.31357*1.31357/(trkPt*trkPt+1.31357*1.31357)+0.880579+7.56199e+06*TMath::Exp(-trkPt*3.08451e+01));
    else
      meanX = TMath::Min(6.85, 3.89994*1.16240*1.16240/(trkPt*trkPt+1.16240*1.16240)-0.120787+2.20275e+05*TMath::Exp(-trkPt*2.40913e+01));
  }

  Double_t rx = (trkDx-meanX)/sx;
  Double_t rz = (trkDz-meanZ)/sz;

  return TMath::Sqrt(rx*rx+rz*rz);
}

//-----------------------------------------------------------------------------
Double_t AliCaloClusterInfo::TestDisp()
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

  return R2;
}

//-----------------------------------------------------------------------------
Bool_t AliCaloClusterInfo::IsInFiducialRegion(Int_t cellX, Int_t cellZ)
{
  const Int_t edgeX = 2;
  const Int_t edgeZ = 2;
  if (cellX >edgeX && cellX<(65-edgeX) && cellZ>edgeZ && cellZ <(57-edgeZ)) return kTRUE;

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
