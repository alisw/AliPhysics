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
#include "TObjArray.h"
#include "TDecompLU.h"
#include "TRandom.h"

#include "AliMFTAnalysisTools.h"

ClassImp(AliMFTAnalysisTools)

//====================================================================================================================================================

Bool_t AliMFTAnalysisTools::ExtrapAODMuonToZ(AliAODTrack *muon, Double_t z, Double_t xy[2]) {

  if (!(muon->PzAtDCA()!=0)) return kFALSE;

  AliMUONTrackParam *param = new AliMUONTrackParam();

  param -> SetNonBendingCoor(muon->XAtDCA());
  param -> SetBendingCoor(muon->YAtDCA());
  param -> SetZ(AliMFTConstants::fZEvalKinem);
  param -> SetNonBendingSlope(muon->PxAtDCA()/muon->PzAtDCA());
  param -> SetBendingSlope(muon->PyAtDCA()/muon->PzAtDCA());
  param -> SetInverseBendingMomentum( muon->Charge() * (1./muon->PzAtDCA()) / (TMath::Sqrt(1+TMath::Power(muon->PyAtDCA()/muon->PzAtDCA(),2))) );

  AliMUONTrackExtrap::ExtrapToZ(param, z);
  xy[0] = param->GetNonBendingCoor();
  xy[1] = param->GetBendingCoor();

  delete param;

  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMFTAnalysisTools::ExtrapAODMuonToZ(AliAODTrack *muon, Double_t z, Double_t xy[2], TLorentzVector &kinem) {

  if (!(muon->PzAtDCA()!=0)) return kFALSE;

  AliMUONTrackParam *param = new AliMUONTrackParam();

  param -> SetNonBendingCoor(muon->XAtDCA());
  param -> SetBendingCoor(muon->YAtDCA());
  param -> SetZ(AliMFTConstants::fZEvalKinem);
  param -> SetNonBendingSlope(muon->PxAtDCA()/muon->PzAtDCA());
  param -> SetBendingSlope(muon->PyAtDCA()/muon->PzAtDCA());
  param -> SetInverseBendingMomentum( muon->Charge() * (1./muon->PzAtDCA()) / (TMath::Sqrt(1+TMath::Power(muon->PyAtDCA()/muon->PzAtDCA(),2))) );

  AliMUONTrackExtrap::ExtrapToZ(param, z);
  xy[0] = param->GetNonBendingCoor();
  xy[1] = param->GetBendingCoor();

  Double_t massMu = TDatabasePDG::Instance()->GetParticle("mu-")->Mass();
  Double_t energy = TMath::Sqrt(massMu*massMu + param->Px()*param->Px() + param->Py()*param->Py() + param->Pz()*param->Pz());
  
  kinem.SetPxPyPzE(param->Px(), param->Py(), param->Pz(), energy);

  delete param;

  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMFTAnalysisTools::ExtrapAODMuonToZ(AliAODTrack *muon, Double_t z, Double_t xy[2], TLorentzVector &kinem, TMatrixD &cov) {

  // Extrapolate muon to a given z providing the corresponding (x,y) position and updating kinematics and covariance matrix

  if (!(muon->PzAtDCA()!=0)) return kFALSE;

  AliMUONTrackParam *param = new AliMUONTrackParam();

  param -> SetNonBendingCoor(muon->XAtDCA());
  param -> SetBendingCoor(muon->YAtDCA());
  param -> SetZ(AliMFTConstants::fZEvalKinem);
  param -> SetNonBendingSlope(muon->PxAtDCA()/muon->PzAtDCA());
  param -> SetBendingSlope(muon->PyAtDCA()/muon->PzAtDCA());
  param -> SetInverseBendingMomentum( muon->Charge() * (1./muon->PzAtDCA()) / (TMath::Sqrt(1+TMath::Power(muon->PyAtDCA()/muon->PzAtDCA(),2))) );

  param -> SetCovariances(ConvertCovMatrixAOD2MUON(muon));

  AliMUONTrackExtrap::ExtrapToZCov(param, z);
  xy[0] = param->GetNonBendingCoor();
  xy[1] = param->GetBendingCoor();

  Double_t massMu = TDatabasePDG::Instance()->GetParticle("mu-")->Mass();
  Double_t energy = TMath::Sqrt(massMu*massMu + param->Px()*param->Px() + param->Py()*param->Py() + param->Pz()*param->Pz());
  
  kinem.SetPxPyPzE(param->Px(), param->Py(), param->Pz(), energy);

  cov = param->GetCovariances();

  delete param;

  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMFTAnalysisTools::ExtrapAODMuonToXY(AliAODTrack *muon, Double_t xy[2], Double_t &zFinal, TLorentzVector &kinem, TMatrixD &cov) {

  // Find the point of closest approach between the muon and the direction parallel to the z-axis defined by the given (x,y)
  // Provide the z of the above point as weel as the updated kinematics and covariance matrix

  // We look for the above-defined PCA

  AliMUONTrackParam *param = new AliMUONTrackParam();
  param -> SetNonBendingCoor(muon->XAtDCA());
  param -> SetBendingCoor(muon->YAtDCA());
  param -> SetZ(AliMFTConstants::fZEvalKinem);
  param -> SetNonBendingSlope(muon->PxAtDCA()/muon->PzAtDCA());
  param -> SetBendingSlope(muon->PyAtDCA()/muon->PzAtDCA());
  param -> SetInverseBendingMomentum( muon->Charge() * (1./muon->PzAtDCA()) / (TMath::Sqrt(1+TMath::Power(muon->PyAtDCA()/muon->PzAtDCA(),2))) );
  
  // here we want to understand in which direction we have to search the minimum...
  
  Double_t step = 1.;  // initial step, in cm
  Double_t startPoint = 0.;
  
  Double_t r[3]={0}, z[3]={startPoint-step, startPoint, startPoint+step};
  
  TVector3 **points = new TVector3*[2];     // points[0] for the muon, points[1] for the direction parallel to the z-axis defined by the given (x,y)
  
  for (Int_t i=0; i<3; i++) {
    AliMUONTrackExtrap::ExtrapToZ(param, z[i]);
    points[0] = new TVector3(param->GetNonBendingCoor(),param->GetBendingCoor(),z[i]);
    points[1] = new TVector3(xy[0],xy[1],z[i]);
    r[i] = GetDistanceBetweenPoints(points,2);
    for (Int_t iMu=0; iMu<2; iMu++) delete points[iMu];
  }
  
  Int_t researchDirection = 0;
  
  if      (r[0]>r[1] && r[1]>r[2]) researchDirection = +1;   // towards z positive
  else if (r[0]<r[1] && r[1]<r[2]) researchDirection = -1;   // towards z negative
  else if (r[0]<r[1] && r[1]>r[2]) {
    printf("E-AliMFTAnalysisTools::ExtrapAODMuonToXY: Point of closest approach cannot be found (no minima)\n");
    delete param;
    delete points;
    return kFALSE;
  }
  
  while (TMath::Abs(researchDirection)>0.5) {
      
    if (researchDirection>0) {
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
      printf("E-AliMFTAnalysisTools::ExtrapAODMuonToXY: Point of closest approach cannot be found (no minima in the fiducial region)\n");
      delete param;
      delete points;
      return kFALSE;
    }
    
    for (Int_t i=0; i<3; i++) {
      AliMUONTrackExtrap::ExtrapToZ(param, z[i]);
      points[0] = new TVector3(param->GetNonBendingCoor(),param->GetBendingCoor(),z[i]);
      points[1] = new TVector3(xy[0],xy[1],z[i]);
      r[i] = GetDistanceBetweenPoints(points,2);
      for (Int_t iMu=0; iMu<2; iMu++) delete points[iMu];
    }

    researchDirection=0;
    if      (r[0]>r[1] && r[1]>r[2]) researchDirection = +1;   // towards z positive
    else if (r[0]<r[1] && r[1]<r[2]) researchDirection = -1;   // towards z negative
    
  }
  
  // now we now that the minimum is between z[0] and z[2] and we search for it
  
  step *= 0.5;
  while (step>AliMFTConstants::fPrecisionPointOfClosestApproach) {
    z[0] = z[1]-step;
    z[2] = z[1]+step;
    for (Int_t i=0; i<3; i++) {
      AliMUONTrackExtrap::ExtrapToZ(param, z[i]);
      points[0] = new TVector3(param->GetNonBendingCoor(),param->GetBendingCoor(),z[i]);
      points[1] = new TVector3(xy[0],xy[1],z[i]);
      r[i] = GetDistanceBetweenPoints(points,2);
      for (Int_t iMu=0; iMu<2; iMu++) delete points[iMu];
    }
    if      (r[0]<r[1]) z[1] = z[0];
    else if (r[2]<r[1]) z[1] = z[2];
    else step *= 0.5;
  }
  
  zFinal = z[1];

  Double_t xyMuon[2] = {0};
  ExtrapAODMuonToZ(muon, zFinal, xyMuon, kinem, cov);

  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMFTAnalysisTools::GetAODMuonOffset(AliAODTrack *muon, Double_t xv, Double_t yv, Double_t zv, Double_t &offset) {

  Double_t xy[2] = {0};
  ExtrapAODMuonToZ(muon, zv, xy);
  
  offset = TMath::Sqrt((xv-xy[0])*(xv-xy[0]) + (yv-xy[1])*(yv-xy[1]));

  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMFTAnalysisTools::GetAODMuonOffsetSmeared(AliAODTrack *muon, Double_t xv, Double_t yv, Double_t zv,
						    Double_t smearOffsetX, Double_t smearOffsetY, Double_t &offset) {

  // Evaluate transverse offset adding to it an additional smearing (independently along the x and y directions)
  
  Double_t xy[2] = {0};
  ExtrapAODMuonToZ(muon, zv, xy);

  xy[0] = gRandom->Gaus(xy[0], smearOffsetX);
  xy[1] = gRandom->Gaus(xy[1], smearOffsetY);
  
  offset = TMath::Sqrt((xv-xy[0])*(xv-xy[0]) + (yv-xy[1])*(yv-xy[1]));

  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMFTAnalysisTools::GetAODMuonOffsetZ(AliAODTrack *muon, Double_t xv, Double_t yv, Double_t zv, Double_t &offset) {

  Double_t xy[2] = {xv, yv};
  Double_t zFinal = 0;
  TLorentzVector kinem(0,0,0,0);
  TMatrixD cov(5,5);

  ExtrapAODMuonToXY(muon, xy, zFinal, kinem, cov);

  offset = TMath::Abs(zFinal - zv);

  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMFTAnalysisTools::GetAODMuonWeightedOffset(AliAODTrack *muon, Double_t xv, Double_t yv, Double_t zv, Double_t &offset) {

  Double_t xy[2] = {0};
  TLorentzVector kinem(0,0,0,0);
  TMatrixD cov(5,5);

  ExtrapAODMuonToZ(muon, zv, xy, kinem, cov);

  TMatrixD covCoordinates(2,2);
  covCoordinates(0,0) = cov(0,0);
  covCoordinates(0,1) = cov(0,2);
  covCoordinates(1,0) = cov(2,0);
  covCoordinates(1,1) = cov(2,2);
  
  if (covCoordinates.Determinant() < covCoordinates.GetTol()) return kFALSE;

  if (TDecompLU::InvertLU(covCoordinates,covCoordinates.GetTol(),0)) {

    TMatrixD covCoordinatesInverse = covCoordinates;
    Double_t dX = xy[0] - xv;
    Double_t dY = xy[1] - yv;
    
    offset = TMath::Sqrt(0.5*(dX*dX*covCoordinatesInverse(0,0) +
			      dY*dY*covCoordinatesInverse(1,1) +
			      2.*dX*dY*covCoordinatesInverse(0,1)));
    
    return kTRUE;

  }
  
  return kFALSE;

}

//====================================================================================================================================================

Double_t AliMFTAnalysisTools::GetPseudoProperDecayTimeXY(Double_t xVtx,  Double_t yVtx, 
							 Double_t xDimu, Double_t yDimu, 
							 Double_t mDimu, Double_t ptDimu) {
  
  // pseudo-proper decay time of a particle produced in the primary vertex and decaying into a dimuon (+ X)
  // evaluated using the transverse degree of freedom of the decay topology 

  if (ptDimu != 0) {
    Double_t decayLengthXY = TMath::Sqrt((xVtx-xDimu)*(xVtx-xDimu)+(yVtx-yDimu)*(yVtx-yDimu));
    return (decayLengthXY * mDimu/ptDimu)/TMath::Ccgs()*1E12;   // in ps
  }
  
  return -99999999;
  
}

//====================================================================================================================================================

Double_t AliMFTAnalysisTools::GetPseudoProperDecayTimeZ(Double_t zVtx, 
							Double_t zDimu, 
							Double_t mDimu, Double_t pzDimu) {

  // pseudo-proper decay time of a particle produced in the primary vertex and decaying into a dimuon (+ X)
  // evaluated using the longitudinal degree of freedom of the decay topology 

  if (pzDimu != 0) {
    Double_t decayLengthZ = zDimu - zVtx;
    return (decayLengthZ * mDimu/pzDimu)/TMath::Ccgs()*1E12;  // in ps
  }

  return -99999999;

}

//====================================================================================================================================================

Bool_t AliMFTAnalysisTools::CalculatePCA(AliAODDimuon *dimuon, Double_t *pca, Double_t &pcaQuality, TLorentzVector &kinem) {

  TObjArray *muons = new TObjArray();
  muons -> Add(dimuon->GetMu(0));
  muons -> Add(dimuon->GetMu(1));
  
  Bool_t result = CalculatePCA(muons, pca, pcaQuality, kinem);
  delete muons;
  return result;

}

//====================================================================================================================================================

Bool_t AliMFTAnalysisTools::CalculatePCA(TObjArray *muons, Double_t *pca, Double_t &pcaQuality, TLorentzVector &kinem) {
  
  const Int_t nMuons = muons->GetEntriesFast();
  if (nMuons<2 || nMuons>AliMFTConstants::fNMaxMuonsForPCA) {
    printf("W-AliMFTAnalysisTools::CalculatePCA: number of muons not valid\n");
    return kFALSE;
  }

  Double_t fXPointOfClosestApproach=0, fYPointOfClosestApproach=0, fZPointOfClosestApproach=0;

  AliAODTrack *muon[AliMFTConstants::fNMaxMuonsForPCA]        = {0};
  AliMUONTrackParam *param[AliMFTConstants::fNMaxMuonsForPCA] = {0};

  // Finding AliMUONTrackParam objects for each muon
  
  for (Int_t iMu=0; iMu<nMuons; iMu++) {
    muon[iMu] = (AliAODTrack*) muons->At(iMu);
    if (TMath::Abs(muon[iMu]->PzAtDCA())<1.e-6) {
      for(Int_t i=0;i<iMu;i++) delete param[i];
      return kFALSE;
    }
    param[iMu] = new AliMUONTrackParam();
    param[iMu] -> SetNonBendingCoor(muon[iMu]->XAtDCA());
    param[iMu] -> SetBendingCoor(muon[iMu]->YAtDCA());
    param[iMu] -> SetZ(AliMFTConstants::fZEvalKinem);
    param[iMu] -> SetNonBendingSlope(muon[iMu]->PxAtDCA()/muon[iMu]->PzAtDCA());
    param[iMu] -> SetBendingSlope(muon[iMu]->PyAtDCA()/muon[iMu]->PzAtDCA());
    param[iMu] -> SetInverseBendingMomentum( muon[iMu]->Charge() * (1./muon[iMu]->PzAtDCA()) / (TMath::Sqrt(1+TMath::Power(muon[iMu]->PyAtDCA()/muon[iMu]->PzAtDCA(),2))) );
  }
  
  // here we want to understand in which direction we have to search the minimum...
  
  Double_t step = 1.;  // initial step, in cm
  Double_t startPoint = 0.;
  
  Double_t r[3]={0}, z[3]={startPoint-step, startPoint, startPoint+step};
  
  TVector3 **points = new TVector3*[AliMFTConstants::fNMaxMuonsForPCA];
  
  for (Int_t i=0; i<3; i++) {
    for (Int_t iMu=0; iMu<nMuons; iMu++) {
      // if (TMath::Abs(param[iMu]->GetInverseBendingMomentum())<1.) {
      // 	printf("W-AliMFTAnalysisTools::CalculatePCA: Evoiding floating point exception in PCA finding\n");
      // 	return kFALSE;
      // }
      AliMUONTrackExtrap::ExtrapToZ(param[iMu], z[i]);
      points[iMu] = new TVector3(param[iMu]->GetNonBendingCoor(),param[iMu]->GetBendingCoor(),z[i]);
    }
    r[i] = GetDistanceBetweenPoints(points,nMuons);
    for (Int_t iMu=0; iMu<nMuons; iMu++) delete points[iMu]; 
  }
  
  Int_t researchDirection = 0;
  
  if      (r[0]>r[1] && r[1]>r[2]) researchDirection = +1;   // towards z positive
  else if (r[0]<r[1] && r[1]<r[2]) researchDirection = -1;   // towards z negative
  else if (r[0]<r[1] && r[1]>r[2]) {
    printf("E-AliMFTAnalysisTools::CalculatePCA: Point of closest approach cannot be found for dimuon (no minima)\n");
    for (Int_t iMu=0;iMu<nMuons;iMu++) delete param[iMu];
    delete points;
    return kFALSE;
  }	
  
  while (TMath::Abs(researchDirection)>0.5) {
      
    if (researchDirection>0) {
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
      printf("E-AliMFTAnalysisTools::CalculatePCA: Point of closest approach cannot be found for dimuon (no minima in the fiducial region)\n");
      for (Int_t iMu=0;iMu<nMuons;iMu++) delete param[iMu];
      delete points;
      return kFALSE;
    }
    
    for (Int_t i=0; i<3; i++) {
      for (Int_t iMu=0; iMu<nMuons; iMu++) {
	AliMUONTrackExtrap::ExtrapToZ(param[iMu], z[i]);
	points[iMu] = new TVector3(param[iMu]->GetNonBendingCoor(),param[iMu]->GetBendingCoor(),z[i]);
      }      
      r[i] = GetDistanceBetweenPoints(points,nMuons);
      for (Int_t iMu=0;iMu<nMuons;iMu++) delete points[iMu];
    }
    researchDirection=0;
    if      (r[0]>r[1] && r[1]>r[2]) researchDirection = +1;   // towards z positive
    else if (r[0]<r[1] && r[1]<r[2]) researchDirection = -1;   // towards z negative
    
  }
  
  // now we now that the minimum is between z[0] and z[2] and we search for it
  
  Int_t nSteps = 0;

  step *= 0.5;
  while (step>AliMFTConstants::fPrecisionPointOfClosestApproach) {
    z[0] = z[1]-step;
    z[2] = z[1]+step;
    for (Int_t i=0; i<3; i++) {
      for (Int_t iMu=0; iMu<nMuons; iMu++) {
	AliMUONTrackExtrap::ExtrapToZ(param[iMu], z[i]);
	points[iMu] = new TVector3(param[iMu]->GetNonBendingCoor(),param[iMu]->GetBendingCoor(),z[i]);
      }      
      r[i] = GetDistanceBetweenPoints(points,nMuons);
      for (Int_t iMu=0;iMu<nMuons;iMu++) delete points[iMu];
    }
    //    printf("Step #%d : %f  %f  %f\n",nSteps,r[0],r[1],r[2]);
    if      (r[0]<r[1]) z[1] = z[0];
    else if (r[2]<r[1]) z[1] = z[2];
    else step *= 0.5;
    nSteps++;
  }

  // if (TMath::Abs(z[1]-1.)<0.1) {
  //   printf("Minimum found in %f in %d steps. Step = %f. p1->X() = %f, p2->X() = %f, p1->Y() = %f, p2->Y() = %f\n",
  // 	   z[1], nSteps, step, points[0]->X(), points[1]->X(), points[0]->Y(), points[1]->Y());
  //   printf("m1->X() = %f, m2->X() = %f, m1->Y() = %f, m2->Y() = %f\n",
  // 	   muon[0]->XAtDCA(),muon[1]->XAtDCA(), muon[0]->YAtDCA(),muon[1]->YAtDCA());
  // }

  // Once z of minimum is found, we evaluate the x and y coordinates by averaging over the contributing tracks
  
  fZPointOfClosestApproach = z[1];
  fXPointOfClosestApproach = 0.;
  fYPointOfClosestApproach = 0.;
  for (Int_t iMu=0; iMu<nMuons; iMu++) {
    AliMUONTrackExtrap::ExtrapToZ(param[iMu], fZPointOfClosestApproach);
    fXPointOfClosestApproach += param[iMu]->GetNonBendingCoor();
    fYPointOfClosestApproach += param[iMu]->GetBendingCoor();
  }
  fXPointOfClosestApproach /= Double_t(nMuons);
  fYPointOfClosestApproach /= Double_t(nMuons);
  
  pca[0] = fXPointOfClosestApproach;
  pca[1] = fYPointOfClosestApproach;
  pca[2] = fZPointOfClosestApproach;
  
  // Evaluating the kinematics of the N-muon
  
  Double_t pTot[3] = {0};
  Double_t ene = 0.;
  Double_t massMu = TDatabasePDG::Instance()->GetParticle("mu-")->Mass();
  for (Int_t iMu=0; iMu<nMuons; iMu++) {
    pTot[0] += param[iMu]->Px();
    pTot[1] += param[iMu]->Py();
    pTot[2] += param[iMu]->Pz();
    ene += TMath::Sqrt(massMu*massMu + param[iMu]->Px()*param[iMu]->Px() + param[iMu]->Py()*param[iMu]->Py() + param[iMu]->Pz()*param[iMu]->Pz());
  }
  
  kinem.SetPxPyPzE(pTot[0], pTot[1], pTot[2], ene);
  
  // Evaluating the PCA quality of the N-muon
  
  Double_t sum=0.,squareSum=0.;
  for (Int_t iMu=0; iMu<nMuons; iMu++) {
    Double_t wOffset = 0;
    if (!GetAODMuonWeightedOffset(muon[iMu],fXPointOfClosestApproach, fYPointOfClosestApproach, fZPointOfClosestApproach, wOffset)) {
      for(Int_t jMu=0;jMu<nMuons;jMu++) delete param[jMu];
      delete points;
      return kFALSE;
    }
    Double_t f = TMath::Exp(-0.5 * wOffset);
    sum += f;
    squareSum += f*f;
  }
  if (sum > 0.) pcaQuality =  (sum-squareSum/sum) / (nMuons-1);
  else pcaQuality = 0.;
  
  for(Int_t iMu=0;iMu<nMuons;iMu++) delete param[iMu];
  delete points;
  return kTRUE;
  
}

//=========================================================================================================================

Double_t AliMFTAnalysisTools::GetDistanceBetweenPoints(TVector3 **points, Int_t nPoints) {
  
  if (nPoints>AliMFTConstants::fNMaxMuonsForPCA) {
    printf("W-AliMFTAnalysisTools::GetDistanceBetweenPoints: number of points not valid\n");
    return 1.e9;
  }
  
  if (nPoints<2) return 0.;
  if (nPoints<3) return TMath::Sqrt( (points[0]->X()-points[1]->X()) * (points[0]->X()-points[1]->X()) +
				     (points[0]->Y()-points[1]->Y()) * (points[0]->Y()-points[1]->Y()) );
  //				     (points[0]->Z()-points[1]->Z()) * (points[0]->Z()-points[1]->Z()) );

  const Int_t nEdgesMax = ((AliMFTConstants::fNMaxMuonsForPCA) * (AliMFTConstants::fNMaxMuonsForPCA - 1)) / 2;
  
  Int_t startID[nEdgesMax]       = {0};
  Int_t stopID[nEdgesMax]        = {0};
  Double_t edgeLength[nEdgesMax] = {0};

  Bool_t pointStatus[AliMFTConstants::fNMaxMuonsForPCA] = {0};
  
  Int_t nEdges=0;
  for (Int_t i=0; i<nPoints-1; i++) {
    for (Int_t j=i+1; j<nPoints; j++) {
      edgeLength[nEdges] = TMath::Sqrt( (points[i]->X()-points[j]->X()) * (points[i]->X()-points[j]->X()) +
                                        (points[i]->Y()-points[j]->Y()) * (points[i]->Y()-points[j]->Y()) +
                                        (points[i]->Z()-points[j]->Z()) * (points[i]->Z()-points[j]->Z()) );
      stopID[nEdges]  = i;
      startID[nEdges] = j;
      nEdges++;
    }
  }
 
  // Order Edges

  Double_t min = 0;
  Int_t   iMin = 0;

  for (Int_t iEdge=0; iEdge<nEdges-1; iEdge++) {
    min  = edgeLength[iEdge];
    iMin = iEdge;
    for (Int_t j=iEdge+1; j<nEdges; j++) {
      if (edgeLength[j]<min) {
        min  = edgeLength[j];
        iMin = j;
      }
    }
    
    if (iMin != iEdge) {

      Double_t edgeLengthMin = edgeLength[iMin];
      Int_t startIDmin = startID[iMin];
      Int_t stopIDmin  = stopID[iMin];
      
      edgeLength[iMin] = edgeLength[iEdge];
      startID[iMin]    = startID[iEdge];
      stopID[iMin]     = stopID[iEdge];

      edgeLength[iEdge] = edgeLengthMin;
      startID[iEdge]    = startIDmin;
      stopID[iEdge]     = stopIDmin;

    }
    
  }
  
  // Connect

  Double_t length = 0.;
  for (Int_t i=0; i<nEdges; i++) {
    if (!(pointStatus[startID[i]] && pointStatus[stopID[i]])) {
      pointStatus[startID[i]] = kTRUE;
      pointStatus[stopID[i]]  = kTRUE;
      length += edgeLength[i];
    }
  }
  
  return length;
  
}

//====================================================================================================================================================

void AliMFTAnalysisTools::ConvertCovMatrixMUON2AOD(const TMatrixD& covMUON, Double_t covAOD[21]) {

  // Converts the cov matrix from the MUON format (TMatrixD) to the AOD one (Double_t[21])
  // 
  // Cov(x,x)       ... :   cv[0]
  // Cov(x,slopeX)  ... :   cv[1]  cv[2]
  // Cov(x,y)       ... :   cv[3]  cv[4]  cv[5]
  // Cov(x,slopeY)  ... :   cv[6]  cv[7]  cv[8]  cv[9]
  // Cov(x,invP_yz) ... :   cv[10] cv[11] cv[12] cv[13] cv[14]
  // not-used       ... :   cv[15] cv[16] cv[17] cv[18] cv[19] cv[20]

  covAOD[0]  = covMUON(0,0);

  covAOD[1]  = covMUON(1,0);
  covAOD[2]  = covMUON(1,1);

  covAOD[3]  = covMUON(2,0);  
  covAOD[4]  = covMUON(2,1);  
  covAOD[5]  = covMUON(2,2);  

  covAOD[6]  = covMUON(3,0);  
  covAOD[7]  = covMUON(3,1);  
  covAOD[8]  = covMUON(3,2);  
  covAOD[9]  = covMUON(3,3);  

  covAOD[10] = covMUON(4,0);  
  covAOD[11] = covMUON(4,1);  
  covAOD[12] = covMUON(4,2);  
  covAOD[13] = covMUON(4,3);  
  covAOD[14] = covMUON(4,4);  

  covAOD[15] = 0;  
  covAOD[16] = 0;  
  covAOD[17] = 0;  
  covAOD[18] = 0;  
  covAOD[19] = 0;  
  covAOD[20] = 0;  

}

//====================================================================================================================================================

const TMatrixD AliMFTAnalysisTools::ConvertCovMatrixAOD2MUON(AliAODTrack *muon) {

  Double_t covAOD[21] = {0};
  muon -> GetCovarianceXYZPxPyPz(covAOD);

  TMatrixD covMUON(5,5);

  covMUON(0,0) = covAOD[0];
  	         		
  covMUON(1,0) = covAOD[1];
  covMUON(1,1) = covAOD[2];
  	         		
  covMUON(2,0) = covAOD[3];  
  covMUON(2,1) = covAOD[4];  
  covMUON(2,2) = covAOD[5];  
  	         		
  covMUON(3,0) = covAOD[6];  
  covMUON(3,1) = covAOD[7];  
  covMUON(3,2) = covAOD[8];  
  covMUON(3,3) = covAOD[9];  
  	         		
  covMUON(4,0) = covAOD[10];  
  covMUON(4,1) = covAOD[11];  
  covMUON(4,2) = covAOD[12];  
  covMUON(4,3) = covAOD[13];  
  covMUON(4,4) = covAOD[14]; 

  return covMUON;

}

//====================================================================================================================================================

Bool_t AliMFTAnalysisTools::IsCorrectMatch(AliAODTrack *muon) {

  for (Int_t iPlane=0; iPlane<AliMFTConstants::fNMaxPlanes; iPlane++) if (IsWrongCluster(muon, iPlane)) return kFALSE;
  return kTRUE;

}

//====================================================================================================================================================

TString AliMFTAnalysisTools::GetGenerator(Int_t label, AliAODMCHeader* header) {

  // get the name of the generator that produced a given particle

  Int_t partCounter = 0;
  TList *genHeaders = header->GetCocktailHeaders();
  Int_t nGenHeaders = genHeaders->GetEntries();

  for (Int_t i=0; i<nGenHeaders; i++){
    AliGenEventHeader *gh = (AliGenEventHeader*) genHeaders->At(i);
    TString genName = gh->GetName();
    Int_t nPart = gh->NProduced();
    if (label>=partCounter && label<(partCounter+nPart)) return genName;
    partCounter += nPart;
  }

  TString empty="";
  return empty;

}

//====================================================================================================================================================

void AliMFTAnalysisTools::GetTrackPrimaryGenerator(AliAODTrack *track, AliAODMCHeader *header, TClonesArray *arrayMC, TString &nameGen) {

  // method to check if a track comes from a given generator

  Int_t label = TMath::Abs(track->GetLabel());
  nameGen = GetGenerator(label,header);
  
  // In case the particle is not primary nameGen will contain blank spaces. In this case, we search backward for the primary which originated the chain
  
  while (nameGen.IsWhitespace()) {
    AliAODMCParticle *mcPart = (AliAODMCParticle*) arrayMC->At(label);
    if (!mcPart) {
      printf("AliMFTAnalysisTools::GetTrackPrimaryGenerator - BREAK: No valid AliAODMCParticle at label %i\n",label);
      break;
    }
    Int_t motherLabel = mcPart->GetMother();
    if (motherLabel < 0) {
      printf("AliMFTAnalysisTools::GetTrackPrimaryGenerator - BREAK: Reached primary particle without valid mother\n");
      break;
    }
    label = motherLabel;
    nameGen = GetGenerator(label,header);
  }
  
  return;

}

//====================================================================================================================================================

Bool_t AliMFTAnalysisTools::IsTrackInjected(AliAODTrack *track, AliAODMCHeader *header, TClonesArray *arrayMC) {

  // method to check if a track comes from the signal event or from the underlying Hijing event

  TString nameGen;

  GetTrackPrimaryGenerator(track,header,arrayMC,nameGen);
  
  if (nameGen.IsWhitespace() || nameGen.Contains("ijing")) return kFALSE;
  
  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMFTAnalysisTools::TranslateMuon(AliAODTrack *muon, Double_t vtxInitial[3], Double_t vtxFinal[3]) {

  if (!(muon->PzAtDCA()!=0)) return kFALSE;

  AliMUONTrackParam *param = new AliMUONTrackParam();

  Double_t deltaVtx[3] = {0};
  for (Int_t i=0; i<3; i++) deltaVtx[i] = vtxInitial[i] - vtxFinal[i];

  param -> SetNonBendingCoor(muon->XAtDCA());
  param -> SetBendingCoor(muon->YAtDCA());
  param -> SetZ(AliMFTConstants::fZEvalKinem);
  param -> SetNonBendingSlope(muon->PxAtDCA()/muon->PzAtDCA());
  param -> SetBendingSlope(muon->PyAtDCA()/muon->PzAtDCA());
  param -> SetInverseBendingMomentum( muon->Charge() * (1./muon->PzAtDCA()) / (TMath::Sqrt(1+TMath::Power(muon->PyAtDCA()/muon->PzAtDCA(),2))) );

  // This will be interpreted as if the track (produced in an event having the primary vertex in vtxInitial) 
  // were produced in an event having the primary vertex in vtxFinal

  AliMUONTrackExtrap::ExtrapToZ(param, deltaVtx[2]);
  muon->SetXYAtDCA(param->GetNonBendingCoor() - deltaVtx[0], param->GetBendingCoor() - deltaVtx[1]);
  muon->SetPxPyPzAtDCA(param->Px(), param->Py(), param->Pz());

  delete param;

  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMFTAnalysisTools::TranslateMuonToOrigin(AliAODTrack *muon, Double_t vtx[3]) {

  Double_t origin[3] = {0,0,0};

  return TranslateMuon(muon, vtx, origin);

}

//====================================================================================================================================================

Bool_t AliMFTAnalysisTools::IsPDGCharm(Int_t pdgCode) {

  pdgCode = TMath::Abs(pdgCode/100);
  if (pdgCode>9) pdgCode /= 10;
  if (pdgCode == 4 ) return kTRUE;
  else return kFALSE;
  
}

//====================================================================================================================================================

Bool_t AliMFTAnalysisTools::IsPDGBeauty(Int_t pdgCode) {

  pdgCode = TMath::Abs(pdgCode/100);
  if (pdgCode>9) pdgCode /= 10;
  if (pdgCode == 5) return kTRUE;
  else return kFALSE;

}

//====================================================================================================================================================

Bool_t AliMFTAnalysisTools::IsPDGResonance(Int_t pdgCode) {

  Int_t id = pdgCode%100000;
  return (!((id-id%10)%110));

} 

//====================================================================================================================================================
