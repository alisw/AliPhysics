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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  The TRD track seed                                                    //
//                                                                        //
//  Authors:                                                              //
//    Alex Bercuci <A.Bercuci@gsi.de>                                     //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include "TLinearFitter.h"

#include "AliLog.h"
#include "AliMathBase.h"

#include "AliTRDseedV1.h"
#include "AliTRDcluster.h"
#include "AliTRDcalibDB.h"
#include "AliTRDstackLayer.h"
#include "AliTRDrecoParam.h"

#define SEED_DEBUG

ClassImp(AliTRDseedV1)

//____________________________________________________________________
AliTRDseedV1::AliTRDseedV1(Int_t layer, AliTRDrecoParam *p) 
  :AliTRDseed()
  ,fLayer(layer)
  ,fTimeBins(0)
  ,fOwner(kFALSE)
  ,fRecoParam(p)
{
  //
  // Constructor
  //

	//AliInfo("");
	AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
	fTimeBins = cal->GetNumberOfTimeBins();

}

//____________________________________________________________________
AliTRDseedV1::AliTRDseedV1(const AliTRDseedV1 &ref, Bool_t owner)
  :AliTRDseed((AliTRDseed&)ref)
  ,fLayer(ref.fLayer)
  ,fTimeBins(ref.fTimeBins)
  ,fOwner(kFALSE)
  ,fRecoParam(ref.fRecoParam)
{
  //
  // Copy Constructor performing a deep copy
  //

	//AliInfo("");

	if(owner){
		for(int ic=0; ic<fTimeBins; ic++){
			if(!fClusters[ic]) continue;
			fClusters[ic] = new AliTRDcluster(*fClusters[ic]);
		}
		fOwner = kTRUE;
	}

}

//____________________________________________________________________
AliTRDseedV1& AliTRDseedV1::operator=(const AliTRDseedV1 &ref)
{
  //
  // Assignment Operator using the copy function
  //

	//AliInfo("");
	if(this != &ref){
		ref.Copy(*this);
	}
	return *this;

}

//____________________________________________________________________
AliTRDseedV1::~AliTRDseedV1()
{
  //
  // Destructor. The RecoParam object belongs to the underlying tracker.
  //

	//AliInfo(Form("fOwner[%s]", fOwner?"YES":"NO"));

	if(fOwner) delete [] fClusters;
}

//____________________________________________________________________
void AliTRDseedV1::Copy(TObject &ref) const
{
  //
  // Copy function
  //

	//AliInfo("");
	AliTRDseedV1 &target = (AliTRDseedV1 &)ref; 
	
	target.fLayer     = fLayer;
	target.fTimeBins  = fTimeBins;
	target.fRecoParam = fRecoParam;
	AliTRDseed::Copy(target);
}

//____________________________________________________________________
Float_t AliTRDseedV1::GetQuality(Bool_t kZcorr) const
{
  //
  // Returns a quality measurement of the current seed
  //

	Float_t zcorr = kZcorr ? fTilt * (fZProb - fZref[0]) : 0.;
	return .5 * (18.0 - fN2)
		+ 10.* TMath::Abs(fYfit[1] - fYref[1])
		+ 5.* TMath::Abs(fYfit[0] - fYref[0] + zcorr)
		+ 2. * TMath::Abs(fMeanz - fZref[0]) / fPadLength;
}

//____________________________________________________________________
Bool_t	AliTRDseedV1::AttachClustersIter(AliTRDstackLayer *layer
                                       , Float_t quality
                                       , Bool_t kZcorr
                                       , AliTRDcluster *c)
{
  //
  // Iterative process to register clusters to the seed.
  // In iteration 0 we try only one pad-row and if quality not
  // sufficient we try 2 pad-rows (about 5% of tracks cross 2 pad-rows)
  //
	
	if(!fRecoParam){
		AliError("Seed can not be used without a valid RecoParam.");
		return kFALSE;
	}
	
	Float_t  tquality;
	Double_t kroady = fRecoParam->GetRoad1y();
	Double_t kroadz = fPadLength * .5 + 1.;
	
	// initialize configuration parameters
	Float_t zcorr = kZcorr ? fTilt * (fZProb - fZref[0]) : 0.;
	Int_t   niter = kZcorr ? 1 : 2;
	
	Double_t yexp, zexp;
	Int_t ncl = 0;
	// start seed update
	for (Int_t iter = 0; iter < niter; iter++) {
  	//AliInfo(Form("iter = %i", iter));
		ncl = 0;
		for (Int_t iTime = 0; iTime < fTimeBins; iTime++) {
			// define searching configuration
			Double_t dxlayer = layer[iTime].GetX() - fX0;
			if(c){
				zexp = c->GetZ();
				//Try 2 pad-rows in second iteration
				if (iter > 0) {
					zexp = fZref[0] + fZref[1] * dxlayer - zcorr;
					if (zexp > c->GetZ()) zexp = c->GetZ() + fPadLength*0.5;
					if (zexp < c->GetZ()) zexp = c->GetZ() - fPadLength*0.5;
				}
			} else zexp = fZref[0];
			yexp  = fYref[0] + fYref[1] * dxlayer - zcorr;
			// get  cluster
// 			printf("xexp = %3.3f ,yexp = %3.3f, zexp = %3.3f\n",layer[iTime].GetX(),yexp,zexp);
// 			printf("layer[%i].GetNClusters() = %i\n", iTime, layer[iTime].GetNClusters());
			Int_t    index = layer[iTime].SearchNearestCluster(yexp, zexp, kroady, kroadz);
// 			for(Int_t iclk = 0; iclk < layer[iTime].GetNClusters(); iclk++){
// 				AliTRDcluster *testcl = layer[iTime].GetCluster(iclk);
// 				printf("Cluster %i: x = %3.3f, y = %3.3f, z = %3.3f\n",iclk,testcl->GetX(), testcl->GetY(), testcl->GetZ());
// 			}
// 			printf("Index = %i\n",index);
			if (index < 0) continue;
			
			// Register cluster
			AliTRDcluster *cl = (AliTRDcluster*) layer[iTime].GetCluster(index);
			
 			//printf("Cluster %i(0x%x): x = %3.3f, y = %3.3f, z = %3.3f\n", index, cl, cl->GetX(), cl->GetY(), cl->GetZ());

			Int_t globalIndex = layer[iTime].GetGlobalIndex(index);
			fIndexes[iTime]  = globalIndex;
			fClusters[iTime] = cl;
			fX[iTime]        = dxlayer;
			fY[iTime]        = cl->GetY();
			fZ[iTime]        = cl->GetZ();
	
			// Debugging
			ncl++;
		}

#ifdef SEED_DEBUG
// 		Int_t nclusters = 0;
// 		Float_t fD[iter] = 0.;
// 		for(int ic=0; ic<fTimeBins+1; ic++){
// 			AliTRDcluster *ci = fClusters[ic];
// 			if(!ci) continue;
// 			for(int jc=ic+1; jc<fTimeBins+1; jc++){
// 				AliTRDcluster *cj = fClusters[jc];
// 				if(!cj) continue;
// 				fD[iter] += TMath::Sqrt((ci->GetY()-cj->GetY())*(ci->GetY()-cj->GetY())+
// 				(ci->GetZ()-cj->GetZ())*(ci->GetZ()-cj->GetZ()));
// 				nclusters++;
// 			}
// 		}
// 		if(nclusters) fD[iter] /= float(nclusters);
#endif

		AliTRDseed::Update();

		if(IsOK()){
			tquality = GetQuality(kZcorr);
			if(tquality < quality) break;
			else quality = tquality;
		}
		kroadz *= 2.;
	} // Loop: iter
	if (!IsOK()) return kFALSE;

	CookLabels();
	UpdateUsed();
	return kTRUE;	
}

//____________________________________________________________________
Bool_t	AliTRDseedV1::AttachClustersProj(AliTRDstackLayer *layer
                                       , Float_t /*quality*/
                                       , Bool_t kZcorr
                                       , AliTRDcluster *c)
{
  //
  // Projective algorithm to attach clusters to seeding tracklets
  //
  // Parameters
  //
  // Output
  //
  // Detailed description
  // 1. Collapse x coordinate for the full detector plane
  // 2. truncated mean on y (r-phi) direction
  // 3. purge clusters
  // 4. truncated mean on z direction
  // 5. purge clusters
  // 6. fit tracklet
  //	

	if(!fRecoParam){
		AliError("Seed can not be used without a valid RecoParam.");
		return kFALSE;
	}

	const Int_t knTimeBins = 35;
	const Int_t kClusterCandidates = 2 * knTimeBins;
	
	//define roads
	Double_t kroady = fRecoParam->GetRoad1y();
	Double_t kroadz = fPadLength * 1.5 + 1.;
	// correction to y for the tilting angle
	Float_t zcorr = kZcorr ? fTilt * (fZProb - fZref[0]) : 0.;

	// working variables
	AliTRDcluster *clusters[kClusterCandidates];
	Double_t cond[4], yexp[knTimeBins], zexp[knTimeBins],
		yres[kClusterCandidates], zres[kClusterCandidates];
	Int_t ncl, *index = 0x0, tboundary[knTimeBins];
	
	// Do cluster projection
	Int_t nYclusters = 0; Bool_t kEXIT = kFALSE;
	for (Int_t iTime = 0; iTime < fTimeBins; iTime++) {
		fX[iTime] = layer[iTime].GetX() - fX0;
		zexp[iTime] = fZref[0] + fZref[1] * fX[iTime];
		yexp[iTime] = fYref[0] + fYref[1] * fX[iTime] - zcorr;
		
		// build condition and process clusters
		cond[0] = yexp[iTime] - kroady; cond[1] = yexp[iTime] + kroady;
		cond[2] = zexp[iTime] - kroadz; cond[3] = zexp[iTime] + kroadz;
		layer[iTime].GetClusters(cond, index, ncl);
		for(Int_t ic = 0; ic<ncl; ic++){
			c = layer[iTime].GetCluster(index[ic]);
			clusters[nYclusters] = c;
			yres[nYclusters++] = c->GetY() - yexp[iTime];
			if(nYclusters >= kClusterCandidates) {
				AliWarning(Form("Cluster candidates reached limit %d. Some may be lost.", kClusterCandidates));
				kEXIT = kTRUE;
				break;
			}
		}
		tboundary[iTime] = nYclusters;
		if(kEXIT) break;
	}
	
	// Evaluate truncated mean on the y direction
	Double_t mean, sigma;
	AliMathBase::EvaluateUni(nYclusters, yres, mean, sigma, Int_t(nYclusters*.8)-2);
	//purge cluster candidates
	Int_t nZclusters = 0;
	for(Int_t ic = 0; ic<nYclusters; ic++){
		if(yres[ic] - mean > 4. * sigma){
			clusters[ic] = 0x0;
			continue;
		}
		zres[nZclusters++] = clusters[ic]->GetZ() - zexp[clusters[ic]->GetLocalTimeBin()];
	}
	
	// Evaluate truncated mean on the z direction
	AliMathBase::EvaluateUni(nZclusters, zres, mean, sigma, Int_t(nZclusters*.8)-2);
	//purge cluster candidates
	for(Int_t ic = 0; ic<nZclusters; ic++){
		if(zres[ic] - mean > 4. * sigma){
			clusters[ic] = 0x0;
			continue;
		}
	}

	
	// Select only one cluster/TimeBin
	Int_t lastCluster = 0;
	fN2 = 0;
	for (Int_t iTime = 0; iTime < fTimeBins; iTime++) {
		ncl = tboundary[iTime] - lastCluster;
		if(!ncl) continue;
		if(ncl == 1){
			c = clusters[lastCluster];
		} else if(ncl > 1){
			Float_t dold = 9999.; Int_t iptr = lastCluster;
			for(int ic=lastCluster; ic<tboundary[iTime]; ic++){
				if(!clusters[ic]) continue;
				Float_t y = yexp[iTime] - clusters[ic]->GetY();
				Float_t z = zexp[iTime] - clusters[ic]->GetZ();
				Float_t d = y * y + z * z;
				if(d > dold) continue;
				dold = d;
				iptr = ic;
			}
			c = clusters[iptr];
		}
		//Int_t globalIndex = layer[iTime].GetGlobalIndex(index);
		//fIndexes[iTime]  = globalIndex;
		fClusters[iTime] = c;
		fY[iTime]        = c->GetY();
		fZ[iTime]        = c->GetZ();
		lastCluster = tboundary[iTime];
		fN2++;
	}
	
	// number of minimum numbers of clusters expected for the tracklet
	Int_t kClmin = Int_t(fRecoParam->GetFindableClusters()*fTimeBins);
  if (fN2 < kClmin){
		AliWarning(Form("Not enough clusters to fit the tracklet %d [%d].", fN2, kClmin));
    fN2 = 0;
    return kFALSE;
  }
	AliTRDseed::Update();
	
// 	// fit tracklet and update clusters
// 	if(!FitTracklet()) return kFALSE;
// 	UpdateUsed();
	return kTRUE;
}

//____________________________________________________________________
Bool_t AliTRDseedV1::FitTracklet()
{
  //
  // Linear fit of the tracklet
  //
  // Parameters :
  //
  // Output :
  //  True if successful
  //
  // Detailed description
  // 2. Check if tracklet crosses pad row boundary
  // 1. Calculate residuals in the y (r-phi) direction
  // 3. Do a Least Square Fit to the data
  //

	//Float_t  sigmaexp  = 0.05 + TMath::Abs(fYref[1] * 0.25); // Expected r.m.s in y direction
  Float_t  ycrosscor = fPadLength * fTilt * 0.5; // Y correction for crossing
  Float_t  anglecor = fTilt * fZref[1];  // Correction to the angle

	// calculate residuals
	const Int_t knTimeBins = 35;
	Float_t yres[knTimeBins]; // y (r-phi) residuals
	Int_t zint[knTimeBins],   // Histograming of the z coordinate 
	      zout[2*knTimeBins];//
	
	fN = 0;
	for (Int_t iTime = 0; iTime < fTimeBins; iTime++) {
    if (!fClusters[iTime]) continue;
    yres[iTime] = fY[iTime] - fYref[0] - (fYref[1] + anglecor) * fX[iTime];
		zint[fN++] = Int_t(fZ[iTime]);
	}

	// calculate pad row boundary crosses
	Int_t kClmin = Int_t(fRecoParam->GetFindableClusters()*fTimeBins);
	Int_t nz = AliMathBase::Freq(fN, zint, zout, kFALSE);
  fZProb   = zout[0];
  if(nz <= 1) zout[3] = 0;
  if(zout[1] + zout[3] < kClmin) {
		AliWarning(Form("Not enough clusters to fit the cross boundary tracklet %d [%d].", zout[1]+zout[3], kClmin));
		return kFALSE;
	}
  // Z distance bigger than pad - length
  if (TMath::Abs(zout[0]-zout[2]) > fPadLength) zout[3]=0;
 

  Double_t sumw   = 0., 
		sumwx  = 0.,
		sumwx2 = 0.,
		sumwy  = 0.,
		sumwxy = 0.,
		sumwz  = 0.,
		sumwxz = 0.;
	Int_t npads;
  fMPads = 0;
	fMeanz = 0.;
	for(int iTime=0; iTime<fTimeBins; iTime++){
    fUsable[iTime] = kFALSE;
    if (!fClusters[iTime]) continue;
		npads = fClusters[iTime]->GetNPads();

		fUsable[iTime] = kTRUE;
    fN2++;
    fMPads += npads;
    Float_t weight = 1.0;
    if(npads > 5) weight = 0.2;
    else if(npads > 4) weight = 0.5;
    sumw   += weight; 
    sumwx  += fX[iTime] * weight;
    sumwx2 += fX[iTime] * fX[iTime] * weight;
    sumwy  += weight * yres[iTime];
    sumwxy += weight * yres[iTime] * fX[iTime];
    sumwz  += weight * fZ[iTime];    
    sumwxz += weight * fZ[iTime] * fX[iTime];
	}
  if (fN2 < kClmin){
		AliWarning(Form("Not enough clusters to fit the tracklet %d [%d].", fN2, kClmin));
    fN2 = 0;
    return kFALSE;
  }
  fMeanz = sumwz / sumw;
	fNChange = 0;

	// Tracklet on boundary
  Float_t correction = 0;
  if (fNChange > 0) {
    if (fMeanz < fZProb) correction =  ycrosscor;
    if (fMeanz > fZProb) correction = -ycrosscor;
  }

  Double_t det = sumw * sumwx2 - sumwx * sumwx;
  fYfitR[0]    = (sumwx2 * sumwy  - sumwx * sumwxy) / det;
  fYfitR[1]    = (sumw   * sumwxy - sumwx * sumwy)  / det;
  
  fSigmaY2 = 0;
  for (Int_t i = 0; i < fTimeBins+1; i++) {
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

	return kTRUE;
}

//_____________________________________________________________________________
Float_t AliTRDseedV1::FitRiemanTilt(AliTRDseedV1 *cseed, Bool_t terror)
{
  //
  // Fit the Rieman tilt
  //

  // Fitting with tilting pads - kz not fixed
	AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
	Int_t nTimeBins = cal->GetNumberOfTimeBins();
  TLinearFitter fitterT2(4,"hyp4");  
  fitterT2.StoreData(kTRUE);
  Float_t xref2 = (cseed[2].fX0 + cseed[3].fX0) * 0.5; // Reference x0 for z
  
  Int_t npointsT = 0;
  fitterT2.ClearPoints();

  for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
// 		printf("\nLayer %d\n", iLayer);
//     cseed[iLayer].Print();
		if (!cseed[iLayer].IsOK()) continue;
    Double_t tilt = cseed[iLayer].fTilt;

    for (Int_t itime = 0; itime < nTimeBins+1; itime++) {
// 			printf("\ttime %d\n", itime);
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
// 			printf("\tadd point :\n");
// 			for(int i=0; i<5; i++) printf("%f ", uvt[i]);
// 			printf("\n");
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

//___________________________________________________________________
void AliTRDseedV1::Print()
{
  //
  // Printing the seedstatus
  //

	AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
	Int_t nTimeBins = cal->GetNumberOfTimeBins();
	
	printf("Seed status :\n");
	printf("  fTilt      = %f\n", fTilt);
	printf("  fPadLength = %f\n", fPadLength);
	printf("  fX0        = %f\n", fX0);
	for(int ic=0; ic<nTimeBins; ic++) {
          const Char_t *isUsable = fUsable[ic]?"Yes":"No";
	  printf("  %d X[%f] Y[%f] Z[%f] Indexes[%d] clusters[%#x] usable[%s]\n"
                , ic
                , fX[ic]
                , fY[ic]
                , fZ[ic]
                , fIndexes[ic]
                , ((void *) fClusters[ic])
                , isUsable);
        }

	printf("  fYref[0] =%f fYref[1] =%f\n", fYref[0], fYref[1]);
	printf("  fZref[0] =%f fZref[1] =%f\n", fZref[0], fZref[1]);
	printf("  fYfit[0] =%f fYfit[1] =%f\n", fYfit[0], fYfit[1]);
	printf("  fYfitR[0]=%f fYfitR[1]=%f\n", fYfitR[0], fYfitR[1]);
	printf("  fZfit[0] =%f fZfit[1] =%f\n", fZfit[0], fZfit[1]);
	printf("  fZfitR[0]=%f fZfitR[1]=%f\n", fZfitR[0], fZfitR[1]);
	printf("  fSigmaY =%f\n", fSigmaY);
	printf("  fSigmaY2=%f\n", fSigmaY2);            
	printf("  fMeanz  =%f\n", fMeanz);
	printf("  fZProb  =%f\n", fZProb);
	printf("  fLabels[0]=%d fLabels[1]=%d\n", fLabels[0], fLabels[1]);
	printf("  fN      =%d\n", fN);
	printf("  fN2     =%d (>8 isOK)\n",fN2);
	printf("  fNUsed  =%d\n", fNUsed);
	printf("  fFreq   =%d\n", fFreq);
	printf("  fNChange=%d\n",  fNChange);
	printf("  fMPads  =%f\n", fMPads);
	
	printf("  fC      =%f\n", fC);        
	printf("  fCC     =%f\n",fCC);      
	printf("  fChi2   =%f\n", fChi2);  
	printf("  fChi2Z  =%f\n", fChi2Z);

}
