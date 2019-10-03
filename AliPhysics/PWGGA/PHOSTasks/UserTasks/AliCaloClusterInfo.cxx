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
#include <TParticle.h>

#include "AliStack.h"
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
  fLabels(0), fDistToBad(0.), fM02(0.), fM20(0.), fTOF(0.)
{
  //
  // Default constructor
  //
} 

//-------------------------------------------------------------------------------------------
AliCaloClusterInfo::AliCaloClusterInfo(AliVCluster* const clust, Int_t relID[4]):
  TObject(), fLorentzVector(), fModule(0), fTRUNumber(0), fNCells(0), fPIDBit(0x0),
  fLabels(0), fDistToBad(0.), fM02(0.), fM20(0.), fTOF(0.)
{
  //
  // constructor
  //
  this->FillCaloClusterInfo(clust, relID);
} 

//-------------------------------------------------------------------------------------------
AliCaloClusterInfo::AliCaloClusterInfo(const AliCaloClusterInfo &src):
  TObject(src), fLorentzVector(src.fLorentzVector), fModule(src.fModule), fTRUNumber(src.fTRUNumber),fNCells(src.fNCells), fPIDBit(src.fPIDBit),
  fLabels(src.fLabels), fDistToBad(src.fDistToBad), fM02(src.fM02), fM20(src.fM20), fTOF(src.fTOF)
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
  fLabels          = src.fLabels;
  fDistToBad       = src.fDistToBad;
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
void AliCaloClusterInfo::FillCaloClusterInfo(AliVCluster* const clust, Int_t relID[4])
{
  // extract information of calo clusters

  fModule         = relID[0];
  fTRUNumber      = GetTRUNumber(relID[2], relID[3]);
  fNCells         = clust->GetNCells();
  fDistToBad      = clust->GetDistanceToBadChannel();
  fM02            = clust->GetM02();
  fM20            = clust->GetM20();
  fTOF            = clust->GetTOF();

  Double_t cpv  = clust->GetEmcCpvDistance(); // Distance in sigmas filled by tender
  Double_t disp = TestDisp();                 // Dispersion in sigmas filled by tender
  if (cpv  > 2.)      this->SetPIDBit(BIT(0));
  if (cpv  > 4.)      this->SetPIDBit(BIT(2));
  if (disp < 2.5*2.5) this->SetPIDBit(BIT(1));
  if (disp < 1.5*1.5) this->SetPIDBit(BIT(3));

  return;
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
Bool_t AliCaloClusterInfo::IsMergedClusterFromPi0(AliStack* const stack, Int_t &pi0Indx)
{
  pi0Indx = -1;
  if (GetNLabels()>1) {
    Int_t pi0Indx1 = -1, pi0Indx2 = -2;
    if (IsClusterFromPi0(stack, GetLabelAt(0), pi0Indx1) &&
        IsClusterFromPi0(stack, GetLabelAt(1), pi0Indx2) && pi0Indx1 == pi0Indx2) {
      pi0Indx = pi0Indx1;
      return kTRUE;
    } else return kFALSE;
  } else return kFALSE;
}

//-----------------------------------------------------------------------------
Bool_t AliCaloClusterInfo::IsClusterFromCvtedPi0(AliStack* const stack, Bool_t &isConverted, Int_t &pi0Indx)
{

  if (IsClusterFromPi0Pure(stack, GetLabel(), pi0Indx))
    return kTRUE;
  else { 
    isConverted = IsClusterFromPi0Converted(stack, GetLabel(), pi0Indx);
    if (isConverted) return kTRUE;
    else return kFALSE;
  }
}

//-----------------------------------------------------------------------------
Bool_t AliCaloClusterInfo::IsClusterFromPi0(AliStack* const stack, Int_t label, Int_t &pi0Indx)
{
  if (IsClusterFromPi0Pure(stack, label, pi0Indx))
    return kTRUE;
  else if (IsClusterFromPi0Converted(stack, label, pi0Indx)) 
    return kTRUE;
  else return kFALSE;
}

//-----------------------------------------------------------------------------
Bool_t AliCaloClusterInfo::IsClusterFromPi0Pure(AliStack* const stack, Int_t label, Int_t &pi0Indx)
{
  if (label>-1) {
    TParticle* track = stack->Particle(label);
    if (track->GetPdgCode() == 22) {
      pi0Indx = track->GetFirstMother();
      if (pi0Indx>-1 && ((TParticle*)stack->Particle(pi0Indx))->GetPdgCode() == 111)
        return kTRUE;
      else return kFALSE;
    }
    else return kFALSE;
  } else return kFALSE;
}

//-----------------------------------------------------------------------------
Bool_t AliCaloClusterInfo::IsClusterFromPi0Converted(AliStack* const stack, Int_t label, Int_t &pi0Indx)
{
  if (label>-1) {
    TParticle* track = stack->Particle(label);
    if (TMath::Abs(track->GetPdgCode()) == 11)
      if (track->GetFirstMother()>-1 && ((TParticle*)stack->Particle(track->GetFirstMother()))->GetPdgCode() == 22) {
        TParticle *gamma = stack->Particle(track->GetFirstMother());
        pi0Indx = gamma->GetFirstMother();
        if (pi0Indx>-1 && ((TParticle*)stack->Particle(pi0Indx))->GetPdgCode() == 111) {
          Int_t e1 = gamma->GetFirstDaughter();
          Int_t e2 = gamma->GetLastDaughter();
          if (label == (((TParticle*)stack->Particle(e1))->Pt()>((TParticle*)stack->Particle(e2))->Pt() ? e1 : e2))
            return kTRUE;
          else return kFALSE;
        } else return kFALSE;
      } else return kFALSE;
    else return kFALSE;
  } else return kFALSE;
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
