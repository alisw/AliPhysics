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
                                                      
/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  The seed of a local TRD track                                            //  
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include "TLinearFitter.h"

#include "AliTRDseed.h"
#include "AliTRDcluster.h"
#include "AliTRDtracker.h"

ClassImp(AliTRDseed)

//_____________________________________________________________________________
AliTRDseed::AliTRDseed() 
  :TObject()
  ,fTilt(0)
  ,fPadLength(0)
  ,fX0(0)
  ,fSigmaY(0)
  ,fSigmaY2(0)
  ,fMeanz(0)
  ,fZProb(0)
  ,fN(0)
  ,fN2(0)
  ,fNUsed(0)
  ,fFreq(0)
  ,fNChange(0)
  ,fMPads(0)
  ,fC(0)
  ,fCC(0)
  ,fChi2(1.0e10)
  ,fChi2Z(1.0e10)
{
  //
  // Default constructor  
  //

  for (Int_t i = 0; i < 25; i++) {
    fX[i]        = 0;   // x position
    fY[i]        = 0;   // y position
    fZ[i]        = 0;   // z position
    fIndexes[i]  = 0;   // Indexes
    fClusters[i] = 0x0; // Clusters
    fUsable[i]   = 0;   // Indication  - usable cluster
  }

  for (Int_t i = 0; i < 2; i++) {
    fYref[i]     = 0;   // Reference y
    fZref[i]     = 0;   // Reference z
    fYfit[i]     = 0;   // Y fit position +derivation
    fYfitR[i]    = 0;   // Y fit position +derivation
    fZfit[i]     = 0;   // Z fit position
    fZfitR[i]    = 0;   // Z fit position
    fLabels[i]   = 0;   // Labels
  }

}

//_____________________________________________________________________________
AliTRDseed::AliTRDseed(const AliTRDseed &s)
  :TObject(s)
  ,fTilt(0)
  ,fPadLength(0)
  ,fX0(0)
  ,fSigmaY(0)
  ,fSigmaY2(0)
  ,fMeanz(0)
  ,fZProb(0)
  ,fN(0)
  ,fN2(0)
  ,fNUsed(0)
  ,fFreq(0)
  ,fNChange(0)
  ,fMPads(0)
  ,fC(0)
  ,fCC(0)
  ,fChi2(0)
  ,fChi2Z(0)
{
  //
  // Copy constructor  
  //

  for (Int_t i = 0; i < 25; i++) {
    fX[i]        = 0;   // x position
    fY[i]        = 0;   // y position
    fZ[i]        = 0;   // z position
    fIndexes[i]  = 0;   // Indexes
    fClusters[i] = 0x0; // Clusters
    fUsable[i]   = 0;   // Indication  - usable cluster
  }

  for (Int_t i = 0; i < 2; i++) {
    fYref[i]     = 0;   // Reference y
    fZref[i]     = 0;   // Reference z
    fYfit[i]     = 0;   // Y fit position +derivation
    fYfitR[i]    = 0;   // Y fit position +derivation
    fZfit[i]     = 0;   // Z fit position
    fZfitR[i]    = 0;   // Z fit position
    fLabels[i]   = 0;   // Labels
  }

}

//_____________________________________________________________________________
void AliTRDseed::Reset()
{
  //
  // Reset seed
  //

  for (Int_t i = 0; i < 25; i++) {
    fX[i]        = 0;  // X position
    fY[i]        = 0;  // Y position
    fZ[i]        = 0;  // Z position
    fIndexes[i]  = 0;  // Indexes
    fClusters[i] = 0;  // Clusters
    fUsable[i]   = kFALSE;    
  }

  for (Int_t i = 0; i < 2; i++) {
    fYref[i]     = 0;  // Reference y
    fZref[i]     = 0;  // Reference z
    fYfit[i]     = 0;  // Y fit position +derivation
    fYfitR[i]    = 0;  // Y fit position +derivation
    fZfit[i]     = 0;  // Z fit position
    fZfitR[i]    = 0;  // Z fit position
    fLabels[i]   = -1; // Labels
  }
  fSigmaY  = 0;        // "Robust" sigma in y
  fSigmaY2 = 0;        // "Robust" sigma in y
  fMeanz   = 0;        // Mean vaue of z
  fZProb   = 0;        // Max probbable z
  fMPads   = 0;
  fN       = 0;        // Number of associated clusters
  fN2      = 0;        // Number of not crossed
  fNUsed   = 0;        // Number of used clusters
  fNChange = 0;        // Change z counter

}

//_____________________________________________________________________________
void AliTRDseed::CookLabels()
{
  //
  // Cook 2 labels for seed
  //

  Int_t labels[200];
  Int_t out[200];
  Int_t nlab = 0;

  for (Int_t i = 0; i < 25; i++) {
    if (!fClusters[i]) continue;
    for (Int_t ilab = 0; ilab < 3; ilab++) {
      if (fClusters[i]->GetLabel(ilab) >= 0) {
	labels[nlab] = fClusters[i]->GetLabel(ilab);
	nlab++;
      }
    }
  }

  Int_t nlab2 = AliTRDtracker::Freq(nlab,labels,out,kTRUE);
  fLabels[0] = out[0];
  if ((nlab2  > 1) && 
      (out[3] > 1)) {
    fLabels[1] = out[2];
  }

}

//_____________________________________________________________________________
void AliTRDseed::UseClusters()
{
  //
  // Use clusters
  //

  for (Int_t i = 0; i < 25; i++) {
    if (!fClusters[i]) continue;
    if (!(fClusters[i]->IsUsed())) fClusters[i]->Use();
  }

}

//_____________________________________________________________________________
void AliTRDseed::Update()
{
  //
  // Update the seed.
  //

  const Float_t kRatio  = 0.8;
  const Int_t   kClmin  = 6;
  const Float_t kmaxtan = 2;

  if (TMath::Abs(fYref[1]) > kmaxtan) return;              // Track inclined too much

  Float_t  sigmaexp  = 0.05 + TMath::Abs(fYref[1] * 0.25); // Expected r.m.s in y direction
  Float_t  ycrosscor = fPadLength * fTilt * 0.5;           // Y correction for crossing 
  fNChange = 0;

  Double_t sumw;
  Double_t sumwx;
  Double_t sumwx2;
  Double_t sumwy;
  Double_t sumwxy;
  Double_t sumwz;
  Double_t sumwxz;

  Int_t    zints[25];                    // Histograming of the z coordinate 
                                         // Get 1 and second max probable coodinates in z
  Int_t    zouts[50];       
  Float_t  allowedz[25];                 // Allowed z for given time bin
  Float_t  yres[25];                     // Residuals from reference
  Float_t  anglecor = fTilt * fZref[1];  // Correction to the angle
  
  
  fN  = 0; 
  fN2 = 0;
  for (Int_t i = 0; i < 25; i++) {
    yres[i] = 10000.0;
    if (!fClusters[i]) continue;
    yres[i] = fY[i] - fYref[0] - (fYref[1] + anglecor) * fX[i];   // Residual y
    zints[fN] = Int_t(fZ[i]);
    fN++;    
  }

  if (fN < kClmin) return;
  Int_t nz = AliTRDtracker::Freq(fN,zints,zouts,kFALSE);
  fZProb   = zouts[0];
  if (nz <= 1) zouts[3] = 0;
  if (zouts[1] + zouts[3] < kClmin) return;
  
  // Z distance bigger than pad - length
  if (TMath::Abs(zouts[0]-zouts[2]) > 12.0) {
    zouts[3]=0;           
  }
  
  Int_t  breaktime = -1;
  Bool_t mbefore   = kFALSE;
  Int_t  cumul[25][2];
  Int_t  counts[2] = { 0, 0 };
  
  if (zouts[3] >= 3) {

    //
    // Find the break time allowing one chage on pad-rows
    // with maximal numebr of accepted clusters
    //
    fNChange = 1;
    for (Int_t i = 0; i < 25; i++) {
      cumul[i][0] = counts[0];
      cumul[i][1] = counts[1];
      if (TMath::Abs(fZ[i]-zouts[0]) < 2) counts[0]++;
      if (TMath::Abs(fZ[i]-zouts[2]) < 2) counts[1]++;
    }
    Int_t  maxcount = 0;
    for (Int_t i = 0; i < 24; i++) {
      Int_t after  = cumul[24][0] - cumul[i][0];
      Int_t before = cumul[i][1];
      if (after + before > maxcount) { 
	maxcount  = after + before; 
	breaktime = i;
	mbefore   = kFALSE;
      }
      after  = cumul[24][1] - cumul[i][1];
      before = cumul[i][0];
      if (after + before > maxcount) { 
	maxcount  = after + before; 
	breaktime = i;
	mbefore   = kTRUE;
      }
    }

    breaktime -= 1;

  }

  for (Int_t i = 0; i < 25; i++) {
    if (i >  breaktime) allowedz[i] =   mbefore  ? zouts[2] : zouts[0];
    if (i <= breaktime) allowedz[i] = (!mbefore) ? zouts[2] : zouts[0];
  }  

  if (((allowedz[0] > allowedz[24]) && (fZref[1] < 0)) || 
      ((allowedz[0] < allowedz[24]) && (fZref[1] > 0))) {
    //
    // Tracklet z-direction not in correspondance with track z direction 
    //
    fNChange = 0;
    for (Int_t i = 0; i < 25; i++) {
      allowedz[i] = zouts[0];  // Only longest taken
    } 
  }
  
  if (fNChange > 0) {
    //
    // Cross pad -row tracklet  - take the step change into account
    //
    for (Int_t i = 0; i < 25; i++) {
      if (!fClusters[i]) continue; 
      if (TMath::Abs(fZ[i] - allowedz[i]) > 2) continue;
      yres[i] = fY[i] - fYref[0] - (fYref[1] + anglecor) * fX[i];   // Residual y
      if (TMath::Abs(fZ[i] - fZProb) > 2) {
	if (fZ[i] > fZProb) yres[i] += fTilt * fPadLength;
	if (fZ[i] < fZProb) yres[i] -= fTilt * fPadLength;
      }
    }
  }
  
  Double_t yres2[25];
  Double_t mean;
  Double_t sigma;
  for (Int_t i = 0; i < 25; i++) {
    if (!fClusters[i]) continue;
    if (TMath::Abs(fZ[i] - allowedz[i]) > 2) continue;
    yres2[fN2] = yres[i];
    fN2++;
  }
  if (fN2 < kClmin) {
    fN2 = 0;
    return;
  }
  EvaluateUni(fN2,yres2,mean,sigma,Int_t(fN2*kRatio-2));
  if (sigma < sigmaexp * 0.8) {
    sigma = sigmaexp;
  }
  fSigmaY = sigma;

  // Reset sums
  sumw   = 0; 
  sumwx  = 0; 
  sumwx2 = 0;
  sumwy  = 0; 
  sumwxy = 0; 
  sumwz  = 0;
  sumwxz = 0;

  fN2    = 0;
  fMeanz = 0;
  fMPads = 0;

  for (Int_t i = 0; i < 25; i++) {

    fUsable[i] = kFALSE;
    if (!fClusters[i]) continue;
    if (TMath::Abs(fZ[i] - allowedz[i]) > 2)  continue;
    if (TMath::Abs(yres[i] - mean) > 4.0 * sigma) continue;
    fUsable[i] = kTRUE;
    fN2++;
    fMPads += fClusters[i]->GetNPads();
    Float_t weight = 1.0;
    if (fClusters[i]->GetNPads() > 4) weight = 0.5;
    if (fClusters[i]->GetNPads() > 5) weight = 0.2;
   
    Double_t x = fX[i];
    sumw   += weight; 
    sumwx  += x * weight; 
    sumwx2 += x*x * weight;
    sumwy  += weight * yres[i];  
    sumwxy += weight * (yres[i]) * x;
    sumwz  += weight * fZ[i];    
    sumwxz += weight * fZ[i] * x;

  }

  if (fN2 < kClmin){
    fN2 = 0;
    return;
  }
  fMeanz = sumwz / sumw;
  Float_t correction = 0;
  if (fNChange > 0) {
    // Tracklet on boundary
    if (fMeanz < fZProb) correction =  ycrosscor;
    if (fMeanz > fZProb) correction = -ycrosscor;
  }

  Double_t det = sumw * sumwx2 - sumwx * sumwx;
  fYfitR[0]    = (sumwx2 * sumwy  - sumwx * sumwxy) / det;
  fYfitR[1]    = (sumw   * sumwxy - sumwx * sumwy)  / det;
  
  fSigmaY2 = 0;
  for (Int_t i = 0; i < 25; i++) {    
    if (!fUsable[i]) continue;
    Float_t delta = yres[i] - fYfitR[0] - fYfitR[1] * fX[i];
    fSigmaY2 += delta*delta;
  }
  fSigmaY2 = TMath::Sqrt(fSigmaY2 / Float_t(fN2-2));
  
  fZfitR[0]  = (sumwx2 * sumwz  - sumwx * sumwxz) / det;
  fZfitR[1]  = (sumw   * sumwxz - sumwx * sumwz)  / det;
  fZfit[0]   = (sumwx2 * sumwz  - sumwx * sumwxz) / det;
  fZfit[1]   = (sumw   * sumwxz - sumwx * sumwz)  / det;
  fYfitR[0] += fYref[0] + correction;
  fYfitR[1] += fYref[1];
  fYfit[0]   = fYfitR[0];
  fYfit[1]   = fYfitR[1];
    
  UpdateUsed();

}

//_____________________________________________________________________________
void AliTRDseed::UpdateUsed()
{
  //
  // Update used seed
  //

  fNUsed = 0;
  for (Int_t i = 0; i < 25; i++) {
    if (!fClusters[i]) continue;
    if ((fClusters[i]->IsUsed())) fNUsed++;
  }

}

//_____________________________________________________________________________
void AliTRDseed::EvaluateUni(Int_t nvectors, Double_t *data, Double_t &mean
                           , Double_t &sigma, Int_t hh)
{
  //
  // Robust estimator in 1D case MI version
  //
  // For the univariate case
  // estimates of location and scatter are returned in mean and sigma parameters
  // the algorithm works on the same principle as in multivariate case -
  // it finds a subset of size hh with smallest sigma, and then returns mean and
  // sigma of this subset
  //

  if (hh == 0) {
    hh = (nvectors + 2) / 2;
  }

  Double_t faclts[] = { 2.6477, 2.5092, 2.3826, 2.2662, 2.1587
                      , 2.0589, 1.9660, 1.879,  1.7973, 1.7203
                      , 1.6473 };
  Int_t    *index   = new Int_t[nvectors];
  TMath::Sort(nvectors, data, index, kFALSE);
  
  Int_t    nquant = TMath::Min(Int_t(Double_t(((hh * 1.0 / nvectors) - 0.5) * 40)) + 1,11);
  Double_t factor = faclts[nquant-1];
  
  Double_t sumx      = 0.0;
  Double_t sumx2     = 0.0;
  Int_t    bestindex = -1;
  Double_t bestmean  = 0.0; 
  Double_t bestsigma = data[index[nvectors-1]] - data[index[0]];   // Maximal possible sigma
  for (Int_t i = 0; i < hh; i++) {
    sumx  += data[index[i]];
    sumx2 += data[index[i]]*data[index[i]];
  }
  
  Double_t norm  = 1.0 / Double_t(hh);
  Double_t norm2 = 1.0 / Double_t(hh - 1);
  for (Int_t i = hh; i < nvectors; i++) {

    Double_t cmean  = sumx*norm;
    Double_t csigma = (sumx2 - hh*cmean*cmean) * norm2;
    if (csigma < bestsigma) {
      bestmean  = cmean;
      bestsigma = csigma;
      bestindex = i - hh;
    }
    
    sumx  += data[index[i]] - data[index[i-hh]];
    sumx2 += data[index[i]]*data[index[i]] - data[index[i-hh]]*data[index[i-hh]];

  }
  
  Double_t bstd = factor * TMath::Sqrt(TMath::Abs(bestsigma));
  mean  = bestmean;
  sigma = bstd;

  delete [] index;

}

//_____________________________________________________________________________
Float_t AliTRDseed::FitRiemanTilt(AliTRDseed * cseed, Bool_t terror)
{
  //
  // Fit the Rieman tilt
  //

  // Fitting with tilting pads - kz not fixed
  TLinearFitter fitterT2(4,"hyp4");  
  fitterT2.StoreData(kTRUE);
  Float_t xref2 = (cseed[2].fX0 + cseed[3].fX0) * 0.5; // Reference x0 for z
  
  Int_t npointsT = 0;
  fitterT2.ClearPoints();

  for (Int_t iLayer = 0; iLayer < 6; iLayer++) {

    if (!cseed[iLayer].IsOK()) continue;
    Double_t tilt = cseed[iLayer].fTilt;

    for (Int_t itime = 0; itime < 25; itime++) {

      if (!cseed[iLayer].fUsable[itime]) continue;
      // x relative to the midle chamber
      Double_t x = cseed[iLayer].fX[itime] + cseed[iLayer].fX0 - xref2;  
      Double_t y = cseed[iLayer].fY[itime];
      Double_t z = cseed[iLayer].fZ[itime];

      //
      // Tilted rieman
      //
      Double_t uvt[6];
      Double_t x2 = cseed[iLayer].fX[itime] + cseed[iLayer].fX0;      // Global x
      Double_t t  = 1.0 / (x2*x2 + y*y);
      uvt[1]  = t;
      uvt[0]  = 2.0 * x2   * uvt[1];
      uvt[2]  = 2.0 * tilt * uvt[1];
      uvt[3]  = 2.0 * tilt *uvt[1] * x;	      
      uvt[4]  = 2.0 * (y + tilt * z) * uvt[1];
      
      Double_t error = 2.0 * uvt[1];
      if (terror) {
        error *= cseed[iLayer].fSigmaY;
      }
      else {
        error *= 0.2; //Default error
      } 
      fitterT2.AddPoint(uvt,uvt[4],error);
      npointsT++;

    }

  }

  fitterT2.Eval();
  Double_t rpolz0 = fitterT2.GetParameter(3);
  Double_t rpolz1 = fitterT2.GetParameter(4);	    

  //
  // Linear fitter  - not possible to make boundaries
  // non accept non possible z and dzdx combination
  // 	    
  Bool_t acceptablez = kTRUE;
  for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
    if (cseed[iLayer].IsOK()) {
      Double_t zT2 = rpolz0 + rpolz1 * (cseed[iLayer].fX0 - xref2);
      if (TMath::Abs(cseed[iLayer].fZProb - zT2) > cseed[iLayer].fPadLength * 0.5 + 1.0) {
	acceptablez = kFALSE;
      }
    }
  }
  if (!acceptablez) {
    Double_t zmf  = cseed[2].fZref[0] + cseed[2].fZref[1] * (xref2 - cseed[2].fX0);
    Double_t dzmf = (cseed[2].fZref[1] + cseed[3].fZref[1]) * 0.5;
    fitterT2.FixParameter(3,zmf);
    fitterT2.FixParameter(4,dzmf);
    fitterT2.Eval();
    fitterT2.ReleaseParameter(3);
    fitterT2.ReleaseParameter(4);
    rpolz0 = fitterT2.GetParameter(3);
    rpolz1 = fitterT2.GetParameter(4);
  }
  
  Double_t chi2TR = fitterT2.GetChisquare() / Float_t(npointsT);  
  Double_t params[3];
  params[0] =  fitterT2.GetParameter(0);
  params[1] =  fitterT2.GetParameter(1);
  params[2] =  fitterT2.GetParameter(2);	    
  Double_t curvature =  1.0 + params[1] * params[1] - params[2] * params[0];

  for (Int_t iLayer = 0; iLayer < 6; iLayer++) {

    Double_t  x  = cseed[iLayer].fX0;
    Double_t  y  = 0;
    Double_t  dy = 0;
    Double_t  z  = 0;
    Double_t  dz = 0;

    // y
    Double_t res2 = (x * params[0] + params[1]);
    res2 *= res2;
    res2  = 1.0 - params[2]*params[0] + params[1]*params[1] - res2;
    if (res2 >= 0) {
      res2 = TMath::Sqrt(res2);
      y    = (1.0 - res2) / params[0];
    }

    //dy
    Double_t x0 = -params[1] / params[0];
    if (-params[2]*params[0] + params[1]*params[1] + 1 > 0) {
      Double_t rm1 = params[0] / TMath::Sqrt(-params[2]*params[0] + params[1]*params[1] + 1); 
      if (1.0/(rm1*rm1) - (x-x0) * (x-x0) > 0.0) {
	Double_t res = (x - x0) / TMath::Sqrt(1.0 / (rm1*rm1) - (x-x0)*(x-x0));
	if (params[0] < 0) res *= -1.0;
	dy = res;
      }
    }
    z  = rpolz0 + rpolz1 * (x - xref2);
    dz = rpolz1;
    cseed[iLayer].fYref[0] = y;
    cseed[iLayer].fYref[1] = dy;
    cseed[iLayer].fZref[0] = z;
    cseed[iLayer].fZref[1] = dz;
    cseed[iLayer].fC       = curvature;
    
  }

  return chi2TR;

}
