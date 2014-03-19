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

//====================================================================================================================================================
//
//      Support class for various common operation on MFT objects
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliAODTrack.h"
#include "AliAODDimuon.h"
#include "TLorentzVector.h"
#include "AliMFTConstants.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "AliLog.h"

#include "AliMFTSupport.h"

ClassImp(AliMFTSupport)

//====================================================================================================================================================

Bool_t AliMFTSupport::ExtrapAODMuonToZ(AliAODTrack *muon, Double_t z, Double_t xy[2]) {

  if (!(muon->Pz()!=0)) return kFALSE;

  AliMUONTrackParam *param = new AliMUONTrackParam();

  param -> SetNonBendingCoor(muon->XAtDCA());
  param -> SetBendingCoor(muon->YAtDCA());
  param -> SetZ(AliMFTConstants::fZEvalKinem);
  param -> SetNonBendingSlope(muon->Px()/muon->Pz());
  param -> SetBendingSlope(muon->Py()/muon->Pz());
  param -> SetInverseBendingMomentum( muon->Charge() * (1./muon->Pz()) / (TMath::Sqrt(1+TMath::Power(muon->Py()/muon->Pz(),2))) );

  AliMUONTrackExtrap::ExtrapToZ(param, z);
  xy[0] = param->GetNonBendingCoor();
  xy[1] = param->GetBendingCoor();

  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMFTSupport::RefitAODDimuonWithCommonVertex(AliAODDimuon *dimuon, Double_t *vertex, TLorentzVector &kinem) {

  Double_t fXPointOfClosestApproach=0, fYPointOfClosestApproach=0, fZPointOfClosestApproach=0;

  AliAODTrack *muon0 = dimuon->GetMu(0);
  AliAODTrack *muon1 = dimuon->GetMu(1);

  if (!(muon0->Pz()!=0 && muon1->Pz()!=0)) return kFALSE;

  AliMUONTrackParam *param0 = new AliMUONTrackParam();
  AliMUONTrackParam *param1 = new AliMUONTrackParam();

  param0 -> SetNonBendingCoor(muon0->XAtDCA());
  param1 -> SetNonBendingCoor(muon1->XAtDCA());

  param0 -> SetBendingCoor(muon0->YAtDCA());
  param1 -> SetBendingCoor(muon1->YAtDCA());

  param0 -> SetZ(AliMFTConstants::fZEvalKinem);
  param1 -> SetZ(AliMFTConstants::fZEvalKinem);
 
  param0 -> SetNonBendingSlope(muon0->Px()/muon0->Pz());
  param1 -> SetNonBendingSlope(muon1->Px()/muon1->Pz());

  param0 -> SetBendingSlope(muon0->Py()/muon0->Pz());
  param1 -> SetBendingSlope(muon1->Py()/muon1->Pz());

  param0 -> SetInverseBendingMomentum( muon0->Charge() * (1./muon0->Pz()) / (TMath::Sqrt(1+TMath::Power(muon0->Py()/muon0->Pz(),2))) );
  param1 -> SetInverseBendingMomentum( muon1->Charge() * (1./muon1->Pz()) / (TMath::Sqrt(1+TMath::Power(muon1->Py()/muon1->Pz(),2))) );

  // here we understand in which direction we have to search the minimum...

  Double_t step = 1.;  // in cm
  Double_t startPoint = 0.;

  Double_t r[3]={0}, z[3]={startPoint, startPoint+step, startPoint+2*step};
  
  for (Int_t i=0; i<3; i++) {
    AliMUONTrackExtrap::ExtrapToZ(param0, z[i]);
    AliMUONTrackExtrap::ExtrapToZ(param1, z[i]);
    Double_t dX = param0->GetNonBendingCoor() - param1->GetNonBendingCoor();
    Double_t dY = param0->GetBendingCoor()    - param1->GetBendingCoor();
    r[i] = TMath::Sqrt(dX*dX + dY*dY);
  }
  
  Double_t researchDirection=0.;
  
  if      (r[0]>r[1] && r[1]>r[2]) researchDirection = +1.;   // towards z positive
  else if (r[0]<r[1] && r[1]<r[2]) researchDirection = -1.;   // towards z negative
  else if (r[0]<r[1] && r[1]>r[2]) { 
    printf("E-AliMFTSupport::RefitAODDimuonWithCommonVertex: Point of closest approach cannot be found for dimuon (no minima)\n");
    return kFALSE;
  }
  
  while (TMath::Abs(researchDirection)>0.5) {
    if (researchDirection>0.) {
      z[0] = z[1];
      z[1] = z[2];
      z[2] = z[1]+researchDirection*step;
    }
    else {
      z[2] = z[1];
      z[1] = z[0];
      z[0] = z[1]+researchDirection*step;
    }
    if (TMath::Abs(z[0])>900.) {
      printf("E-AliMFTSupport::RefitAODDimuonWithCommonVertex: Point of closest approach cannot be found for dimuon (no minima)\n");
      return kFALSE;
    }
    for (Int_t i=0; i<3; i++) {
      AliMUONTrackExtrap::ExtrapToZ(param0, z[i]);
      AliMUONTrackExtrap::ExtrapToZ(param1, z[i]);
      Double_t dX = param0->GetNonBendingCoor() - param1->GetNonBendingCoor();
      Double_t dY = param0->GetBendingCoor()    - param1->GetBendingCoor();
      r[i] = TMath::Sqrt(dX*dX + dY*dY);
    }
    researchDirection=0.;
    if      (r[0]>r[1] && r[1]>r[2]) researchDirection = +1.;   // towards z positive
    else if (r[0]<r[1] && r[1]<r[2]) researchDirection = -1.;   // towards z negative
  }

  // now we now that the minimum is between z[0] and z[2] and we search it
    
  step *= 0.5;
  while (step>AliMFTConstants::fPrecisionPointOfClosestApproach) {
    z[0] = z[1]-step;
    z[2] = z[1]+step;
    for (Int_t i=0; i<3; i++) {
      AliMUONTrackExtrap::ExtrapToZ(param0, z[i]);
      AliMUONTrackExtrap::ExtrapToZ(param1, z[i]);
      Double_t dX = param0->GetNonBendingCoor() - param1->GetNonBendingCoor();
      Double_t dY = param0->GetBendingCoor()    - param1->GetBendingCoor();
      r[i] = TMath::Sqrt(dX*dX + dY*dY);
    }
    if      (r[0]<r[1]) z[1] = z[0];
    else if (r[2]<r[1]) z[1] = z[2];
    else step *= 0.5;
  }
  
  fZPointOfClosestApproach = z[1];
  AliMUONTrackExtrap::ExtrapToZ(param0, fZPointOfClosestApproach);
  AliMUONTrackExtrap::ExtrapToZ(param1, fZPointOfClosestApproach);  
  fXPointOfClosestApproach = 0.5*(param0->GetNonBendingCoor() + param1->GetNonBendingCoor());
  fYPointOfClosestApproach = 0.5*(param0->GetBendingCoor()    + param1->GetBendingCoor());
  
  //-------------------------------------------------------------------------------------------------------------------------------------------------

  vertex[0] = fXPointOfClosestApproach;
  vertex[1] = fYPointOfClosestApproach;
  vertex[2] = fZPointOfClosestApproach;

  Double_t pTot[3] = {0};
  pTot[0] = param0->Px() + param1->Px();
  pTot[1] = param0->Py() + param1->Py();
  pTot[2] = param0->Pz() + param1->Pz();

  Double_t massMu = TDatabasePDG::Instance()->GetParticle("mu-")->Mass();

  Double_t energy0 = TMath::Sqrt(massMu*massMu + param0->Px()*param0->Px() + param0->Py()*param0->Py() + param0->Pz()*param0->Pz());
  Double_t energy1 = TMath::Sqrt(massMu*massMu + param1->Px()*param1->Px() + param1->Py()*param1->Py() + param1->Pz()*param1->Pz());

  kinem.SetPxPyPzE(pTot[0], pTot[1], pTot[2], energy0+energy1);

  return kTRUE;

}

//====================================================================================================================================================
