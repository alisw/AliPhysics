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

#include "AliMFTAnalysisTools.h"

ClassImp(AliMFTAnalysisTools)

//====================================================================================================================================================

Bool_t AliMFTAnalysisTools::ExtrapAODMuonToZ(AliAODTrack *muon, Double_t z, Double_t xy[2]) {

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

  delete param;

  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMFTAnalysisTools::ExtrapAODMuonToZ(AliAODTrack *muon, Double_t z, Double_t xy[2], TLorentzVector &kinem) {

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

  Double_t massMu = TDatabasePDG::Instance()->GetParticle("mu-")->Mass();
  Double_t energy = TMath::Sqrt(massMu*massMu + param->Px()*param->Px() + param->Py()*param->Py() + param->Pz()*param->Pz());
  
  kinem.SetPxPyPzE(param->Px(), param->Py(), param->Pz(), energy);

  delete param;

  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMFTAnalysisTools::ExtrapAODMuonToZ(AliAODTrack *muon, Double_t z, Double_t xy[2], TLorentzVector &kinem, TMatrixD &cov) {

  if (!(muon->Pz()!=0)) return kFALSE;

  AliMUONTrackParam *param = new AliMUONTrackParam();

  param -> SetNonBendingCoor(muon->XAtDCA());
  param -> SetBendingCoor(muon->YAtDCA());
  param -> SetZ(AliMFTConstants::fZEvalKinem);
  param -> SetNonBendingSlope(muon->Px()/muon->Pz());
  param -> SetBendingSlope(muon->Py()/muon->Pz());
  param -> SetInverseBendingMomentum( muon->Charge() * (1./muon->Pz()) / (TMath::Sqrt(1+TMath::Power(muon->Py()/muon->Pz(),2))) );

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

Double_t AliMFTAnalysisTools::GetAODMuonOffset(AliAODTrack *muon, Double_t xv, Double_t yv, Double_t zv) {

  Double_t xy[2] = {0};
  ExtrapAODMuonToZ(muon, zv, xy);
  
  return TMath::Sqrt((xv-xy[0])*(xv-xy[0]) + (yv-xy[1])*(yv-xy[1]));

}

//====================================================================================================================================================

Double_t AliMFTAnalysisTools::GetAODMuonWeightedOffset(AliAODTrack *muon, Double_t xv, Double_t yv, Double_t zv) {

  Double_t xy[2] = {0};
  TLorentzVector kinem(0,0,0,0);
  TMatrixD cov(5,5);

  ExtrapAODMuonToZ(muon, zv, xy, kinem, cov);

  TMatrixD covCoordinates(2,2);
  covCoordinates(0,0) = cov(0,0);
  covCoordinates(0,1) = cov(0,2);
  covCoordinates(1,0) = cov(2,0);
  covCoordinates(1,1) = cov(2,2);
  
  TMatrixD covCoordinatesInverse = covCoordinates.Invert();

  Double_t dX = xy[0] - xv;
  Double_t dY = xy[1] - yv;
  
  Double_t weightedOffset = TMath::Sqrt(0.5*(dX*dX*covCoordinatesInverse(0,0) + 
					     dY*dY*covCoordinatesInverse(1,1) + 
					     2.*dX*dY*covCoordinatesInverse(0,1)));

  return weightedOffset;

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

  return CalculatePCA(muons, pca, pcaQuality, kinem);

}

//====================================================================================================================================================

Bool_t AliMFTAnalysisTools::CalculatePCA(TObjArray *muons, Double_t *pca, Double_t &pcaQuality, TLorentzVector &kinem) {
  
  const Int_t nMuons = muons->GetEntriesFast();
  if (nMuons<2 || nMuons>AliMFTConstants::fNMaxMuonsForPCA) {
    printf("W-AliMFTAnalysisTools::CalculatePCA: number of muons not valid\n");
    return kFALSE;
  }
  
  Double_t fXPointOfClosestApproach=0, fYPointOfClosestApproach=0, fZPointOfClosestApproach=0;
  
  AliAODTrack *muon[nMuons];
  AliMUONTrackParam *param[nMuons];

  // Finding AliMUONTrackParam objects for each muon
  
  for (Int_t iMu=0; iMu<nMuons; iMu++) {
    muon[iMu] = (AliAODTrack*) muons->At(iMu);
    if (TMath::Abs(muon[iMu]->Pz())<1.e-6) return kFALSE;
    param[iMu] = new AliMUONTrackParam();
    param[iMu] -> SetNonBendingCoor(muon[iMu]->XAtDCA());
    param[iMu] -> SetBendingCoor(muon[iMu]->YAtDCA());
    param[iMu] -> SetZ(0.);
    param[iMu] -> SetNonBendingSlope(muon[iMu]->Px()/muon[iMu]->Pz());
    param[iMu] -> SetBendingSlope(muon[iMu]->Py()/muon[iMu]->Pz());
    param[iMu] -> SetInverseBendingMomentum( muon[iMu]->Charge() * (1./muon[iMu]->Pz()) / (TMath::Sqrt(1+TMath::Power(muon[iMu]->Py()/muon[iMu]->Pz(),2))) );
  }

  // here we want to understand in which direction we have to search the minimum...
  
  Double_t step = 1.;  // initial step, in cm
  Double_t startPoint = 0.;

  Double_t r[3]={0}, z[3]={startPoint, startPoint+step, startPoint+2*step};
  
  TVector3 **points = new TVector3*[nMuons];
  
  for (Int_t i=0; i<3; i++) {
    for (Int_t iMu=0; iMu<nMuons; iMu++) {
      AliMUONTrackExtrap::ExtrapToZ(param[iMu], z[i]);
      points[iMu] = new TVector3(param[iMu]->GetNonBendingCoor(),param[iMu]->GetBendingCoor(),z[i]);
    }
    r[i] = GetDistanceBetweenPoints(points,nMuons);
  }

  Int_t researchDirection = 0;
  
  if      (r[0]>r[1] && r[1]>r[2]) researchDirection = +1;   // towards z positive
  else if (r[0]<r[1] && r[1]<r[2]) researchDirection = -1;   // towards z negative
  else if (r[0]<r[1] && r[1]>r[2]) {
    printf("E-AliMFTAnalysisTools::CalculatePCA: Point of closest approach cannot be found for dimuon (no minima)\n");
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
      return kFALSE;
    }
    
    for (Int_t iMu=0; iMu<nMuons; iMu++) delete points[iMu];
    
    for (Int_t i=0; i<3; i++) {
      for (Int_t iMu=0; iMu<nMuons; iMu++) {
        AliMUONTrackExtrap::ExtrapToZ(param[iMu], z[i]);
        points[iMu] = new TVector3(param[iMu]->GetNonBendingCoor(),param[iMu]->GetBendingCoor(),z[i]);
      }      
      r[i] = GetDistanceBetweenPoints(points,nMuons);
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
      for (Int_t iMu=0; iMu<nMuons; iMu++) {
        AliMUONTrackExtrap::ExtrapToZ(param[iMu], z[i]);
        points[iMu] = new TVector3(param[iMu]->GetNonBendingCoor(),param[iMu]->GetBendingCoor(),z[i]);
      }      
      r[i] = GetDistanceBetweenPoints(points,nMuons);
    }
    if      (r[0]<r[1]) z[1] = z[0];
    else if (r[2]<r[1]) z[1] = z[2];
    else step *= 0.5;
  }
  
  delete [] points;

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
    Double_t wOffset = AliMFTAnalysisTools::GetAODMuonWeightedOffset(muon[iMu],fXPointOfClosestApproach, fYPointOfClosestApproach, fZPointOfClosestApproach);
    Double_t f = TMath::Exp(-0.5 * wOffset);
    sum += f;
    squareSum += f*f;
  }
  if (sum > 0.) pcaQuality =  (sum-squareSum/sum) / (nMuons-1);
  else pcaQuality = 0.;

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
				     (points[0]->Y()-points[1]->Y()) * (points[0]->Y()-points[1]->Y()) +
				     (points[0]->Z()-points[1]->Z()) * (points[0]->Z()-points[1]->Z()) );

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

