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

//-------------------------------------------------------------------------
//   AliESDUtils - This is a namespace that temporary provides general 
//                 purpose ESD utilities to avoid unnecessary dependencies
//                 between PWG libraries. To be removed/replaced by AODB
//                 framework.
//      
//-------------------------------------------------------------------------
// Author: Andrei Gheata

#include "AliESDUtils.h"

#include "AliESDEvent.h"
#include "AliESDVZERO.h"
#include "AliVertexerTracks.h"

//______________________________________________________________________________
Float_t AliESDUtils::GetCorrV0(const AliVEvent* esd, Float_t &v0CorrResc, Float_t *v0multChCorr, Float_t *v0multChCorrResc)
{
  // Correct V0 non-linearity, prepare a version rescaled to SPD2 corr.
  // Please describe better parameters...
  Double_t *par0;
  Double_t *par1;
  Double_t *par2;
  
  Double_t par0_137161[64] = { 6.71e-02 , 6.86e-02 , 7.06e-02 , 6.32e-02 , 
			5.91e-02 , 6.07e-02 , 5.78e-02 , 5.73e-02 , 5.91e-02 , 6.22e-02 , 
			5.90e-02 , 6.11e-02 , 5.55e-02 , 5.29e-02 , 5.19e-02 , 5.56e-02 , 
			6.25e-02 , 7.03e-02 , 5.64e-02 , 5.81e-02 , 4.57e-02 , 5.30e-02 , 
			5.13e-02 , 6.43e-02 , 6.27e-02 , 6.48e-02 , 6.07e-02 , 1.01e-01 , 
			6.68e-02 , 7.16e-02 , 6.36e-02 , 5.95e-02 , 2.52e-02 , 2.82e-02 , 
			2.56e-02 , 2.86e-02 , 2.82e-02 , 2.10e-02 , 2.13e-02 , 2.32e-02 , 
			2.75e-02 , 4.34e-02 , 3.78e-02 , 4.52e-02 , 4.11e-02 , 3.89e-02 , 
			4.10e-02 , 3.73e-02 , 4.51e-02 , 5.07e-02 , 5.42e-02 , 4.74e-02 , 
			4.33e-02 , 4.44e-02 , 4.64e-02 , 3.01e-02 , 6.38e-02 , 5.26e-02 , 
			4.99e-02 , 5.26e-02 , 5.47e-02 , 3.84e-02 , 5.00e-02 , 5.20e-02 };
  Double_t par1_137161[64] = { -6.68e-05 , -7.78e-05 , -6.88e-05 , -5.92e-05 , 
			-2.43e-05 , -3.54e-05 , -2.91e-05 , -1.99e-05 , -1.40e-05 , -4.01e-05 , 
			-2.29e-05 , -3.68e-05 , -2.53e-05 , -2.44e-06 , -9.22e-06 , -1.51e-05 , 
			-2.80e-05 , -2.34e-05 , -1.72e-05 , -1.81e-05 , -1.29e-05 , -2.65e-05 , 
			-1.61e-05 , -2.86e-05 , -1.74e-05 , -4.23e-05 , -3.41e-05 , -1.05e-04 , 
			-2.76e-05 , -4.71e-05 , -3.06e-05 , -2.32e-05 , -1.55e-06 , 2.15e-05 , 
			1.40e-05 , 2.16e-05 , 1.21e-05 , 3.05e-06 , 1.67e-05 , -3.84e-06 , 
			3.09e-06 , 1.50e-05 , 3.47e-06 , 4.87e-06 , -3.71e-07 , -1.75e-06 , 
			-1.80e-06 , 9.99e-06 , -6.46e-06 , -4.91e-06 , 1.33e-05 , -2.52e-07 , 
			-3.85e-06 , 4.94e-06 , -2.48e-07 , -1.20e-05 , 2.07e-06 , 6.12e-06 , 
			-1.18e-06 , 4.54e-06 , -1.54e-05 , -1.25e-05 , 1.46e-06 , -6.67e-06 };
  Double_t par2_137161[64] = { 1.29e-08 , 1.51e-08 , 1.43e-08 , 1.11e-08 , 
			5.04e-09 , 6.99e-09 , 5.58e-09 , 4.15e-09 , 4.00e-09 , 8.22e-09 , 
			4.97e-09 , 7.66e-09 , 4.91e-09 , 1.10e-09 , 2.64e-09 , 3.64e-09 , 
			5.76e-09 , 5.46e-09 , 3.38e-09 , 3.47e-09 , 2.43e-09 , 4.13e-09 , 
			2.80e-09 , 5.80e-09 , 3.86e-09 , 7.46e-09 , 5.98e-09 , 2.58e-08 , 
			5.50e-09 , 8.72e-09 , 5.23e-09 , 4.37e-09 , 2.33e-09 , -6.01e-10 , 
			3.99e-11 , -2.02e-10 , 7.67e-10 , 2.03e-09 , 1.17e-10 , 2.56e-09 , 
			1.16e-09 , -4.75e-10 , 1.28e-09 , 1.23e-09 , 1.62e-09 , 1.61e-09 , 
			1.93e-09 , 2.97e-10 , 2.21e-09 , 2.16e-09 , 5.22e-10 , 1.03e-09 , 
			1.56e-09 , 5.00e-10 , 1.01e-09 , 2.93e-09 , 1.05e-09 , 9.96e-11 , 
			1.21e-09 , 7.45e-10 , 3.07e-09 , 2.31e-09 , 6.70e-10 , 1.89e-09 };

  Double_t par0_137366[64] = { 7.12e-02 , 7.34e-02 , 7.39e-02 , 6.54e-02 , 6.11e-02 , 6.31e-02 , 6.15e-02 , 
			       6.00e-02 , 6.10e-02 , 6.49e-02 , 6.17e-02 , 6.33e-02 , 6.00e-02 , 5.48e-02 , 
			       5.44e-02 , 5.81e-02 , 6.49e-02 , 7.07e-02 , 5.91e-02 , 6.18e-02 , 4.82e-02 , 
			       5.67e-02 , 5.36e-02 , 6.60e-02 , 6.37e-02 , 6.78e-02 , 6.31e-02 , 1.04e-01 , 
			       6.91e-02 , 7.32e-02 , 6.61e-02 , 6.16e-02 , 2.64e-02 , 2.81e-02 , 2.64e-02 , 
			       2.85e-02 , 2.87e-02 , 2.18e-02 , 2.19e-02 , 2.43e-02 , 2.81e-02 , 4.37e-02 , 
			       3.90e-02 , 4.66e-02 , 4.24e-02 , 4.09e-02 , 4.21e-02 , 3.88e-02 , 4.83e-02 , 
			       5.23e-02 , 5.44e-02 , 4.85e-02 , 4.42e-02 , 4.58e-02 , 4.74e-02 , 3.14e-02 , 
			       6.31e-02 , 5.30e-02 , 5.01e-02 , 5.33e-02 , 5.70e-02 , 3.95e-02 , 4.98e-02 , 5.31e-02 };
  Double_t par1_137366[64] = { -6.99e-05 , -6.99e-05 , -6.94e-05 , -6.55e-05 , -3.55e-05 , -4.50e-05 , 
			       -3.10e-05 , -2.81e-05 , -2.29e-05 , -3.89e-05 , -2.53e-05 , -4.25e-05 ,
			       -1.87e-05 , -2.01e-05 , -1.53e-05 , -2.14e-05 , -2.86e-05 , -4.70e-05 ,
			       -2.23e-05 , -3.30e-05 ,-9.74e-06 , -2.62e-05 , -1.76e-05 , -2.38e-05 , 
			       -2.40e-05 , -3.43e-05 , -2.75e-05 , -6.86e-05 ,-2.35e-05 , -4.45e-05 , 
			       -2.51e-05 , -2.20e-05 , -1.25e-16 , -2.04e-17 , -2.06e-17 , -3.74e-19 ,
			       -1.18e-18 , -2.02e-15 , -3.78e-06 , -1.26e-06 , -2.71e-06 , -6.23e-17 , 
			       -7.39e-08 , -1.76e-16 , -8.98e-06 , -4.10e-18 , -1.34e-05 , -1.06e-16 , 
			       -3.34e-06 , -1.04e-05 , -5.28e-06 , -7.34e-06 , -1.05e-05 , -7.68e-06 ,
			       -1.78e-05 , -1.19e-05 , -1.78e-05 , -1.34e-06 , -9.23e-06 , -3.34e-06 ,
			       -8.02e-06 , -1.39e-05 , -1.38e-05 , -1.40e-05 };
  Double_t par2_137366[64] = { 1.41e-08 , 1.47e-08 , 1.48e-08 , 1.24e-08 , 6.82e-09 , 8.73e-09 , 6.26e-09 , 
			       5.53e-09 , 5.40e-09 , 7.93e-09 , 5.49e-09 , 8.77e-09 , 4.21e-09 , 3.93e-09 , 
			       3.60e-09 , 4.67e-09 , 5.59e-09 , 8.81e-09 , 3.89e-09 , 6.19e-09 , 1.97e-09 , 
			       4.38e-09 , 3.26e-09 , 5.00e-09 , 4.58e-09 , 6.39e-09 , 5.03e-09 , 1.30e-08 , 
			       4.95e-09 , 8.26e-09 , 4.57e-09 , 4.10e-09 , 2.35e-09 , 2.30e-09 , 2.15e-09 , 
			       2.27e-09 , 2.17e-09 , 2.27e-09 , 2.97e-09 , 2.25e-09 , 1.69e-09 , 1.44e-09 , 
			       1.66e-09 , 1.75e-09 , 2.88e-09 , 1.82e-09 , 3.64e-09 , 1.80e-09 , 1.71e-09 , 
			       2.66e-09 , 3.01e-09 , 1.95e-09 , 2.64e-09 , 2.42e-09 , 3.68e-09 , 2.66e-09 , 
			       3.92e-09 , 1.18e-09 , 2.26e-09 , 1.57e-09 , 2.02e-09 , 2.71e-09 , 2.99e-09 , 3.04e-09 }; 
  
  
  if (esd->GetRunNumber() <= 137165) {
    par0=par0_137161;
    par1=par1_137161;
    par2=par2_137161;
  }  else  {
    par0=par0_137366;
    par1=par1_137366;
    par2=par2_137366;
 }
  //
  Float_t multCorr = 0;
  Float_t multCorr2 = 0;
  Float_t multChCorr[64];
  Float_t multChCorrResc[64];
  AliVVZERO* esdV0 = esd->GetVZEROData();
  for(Int_t i = 0; i < 64; ++i) {
    if (esdV0->TestBit(AliESDVZERO::kCorrectedForSaturation)) {
      multChCorr[i] = esdV0->GetMultiplicity(i);
    }
    else {
      Double_t b = (esdV0->GetMultiplicity(i)*par1[i]-par0[i]);
      Double_t s = (b*b-4.*par2[i]*esdV0->GetMultiplicity(i)*esdV0->GetMultiplicity(i));
      Double_t n = (s<0) ? -b : (-b + TMath::Sqrt(s));
      multChCorr[i] = 2.*esdV0->GetMultiplicity(i)/n*par0[i];
    }
    multCorr += multChCorr[i];
    multChCorrResc[i] = multChCorr[i]/par0[i]/64.;
    multCorr2 += multChCorrResc[i];
  }
  v0CorrResc =  multCorr2;
  if (v0multChCorr)
    for(Int_t i = 0; i < 64; ++i) v0multChCorr[i] = multChCorr[i];
  if (v0multChCorrResc)
    for(Int_t i = 0; i < 64; ++i) v0multChCorrResc[i] = multChCorrResc[i];

  return multCorr;
}

//______________________________________________________________________________
Float_t AliESDUtils::GetCorrSPD2(Float_t spd2raw,Float_t zv)
{
  // renormalize N spd2 clusters at given Zv to acceptance at Zv=0
  const Double_t pars[] = {8.10030e-01,-2.80364e-03,-7.19504e-04};
  zv -= pars[0];
  Float_t corr = 1 + zv*(pars[1] + zv*pars[2]);
  return corr>0 ? spd2raw/corr : -1;
}  

//______________________________________________________________________________
TObjArray*  AliESDUtils::RefitESDVertexTracks(AliESDEvent* esdEv, Int_t algo, const Double_t *cuts)
{
  // Refit ESD VertexTracks and redo tracks->RelateToVertex
  // Default vertexin algorithm is 6 (multivertexer). To use old vertexed, use algo=1
  //
  static AliVertexerTracks* vtFinder = 0;
  static int currRun = 0;
  static int defAlgo = -1;
  static double bkgauss = 0; 
  const Bool_t kVtxConstr = kTRUE;
  //
  if (!vtFinder) { // create vertexer
    vtFinder = new AliVertexerTracks(esdEv->GetMagneticField());
    printf("Initialized vertexer for VertexTracks refit with field %f kG\n",esdEv->GetMagneticField());
    //
    vtFinder->SetITSMode();
    vtFinder->SetConstraintOff();
  }
  //
  if ( (cuts && algo>11) || algo!=defAlgo) { 
    // if cuts array is provided, then interpret algo as the number of parameters in the cuts.
    // otherwise, interpret it as an algorithm ID for hardwired cuts below
    if (cuts) {
      vtFinder->SetCuts((double*)cuts,algo);
      defAlgo = (Int_t)(cuts[10]);
    }
    else {
      const int kNCuts=21;
      double vtCuts[kNCuts] = 
	{1.00e-01,1.00e-01,5.00e-01,3.00e+00,1.00e+00,3.00e+00,1.00e+02,
	 1.00e+03,3.00e+00,3.00e+01,6.00e+00,4.00e+00,7.00e+00,1.00e+03,
	 5.00e+00,5.00e-02,1.00e-03,2.00e+00,1.00e+01,1.00e+00,5.00e+01};
      //
      vtCuts[10] = algo;    
      defAlgo = algo;
      vtFinder->SetCuts(vtCuts,kNCuts);
      printf("Setting vertexing algorithm to %d\n",defAlgo);
    }
  }
  if (defAlgo<0 || defAlgo>AliVertexerTracks::kMultiVertexer) {
    printf("Vertexer algorithms 0:%d are supported... \n",defAlgo);
    return nullptr;
  }
  //
  if (currRun!=esdEv->GetRunNumber() && kVtxConstr) { // update diamond for this run
    double pos[3]={esdEv->GetDiamondX(),esdEv->GetDiamondY(),0};    
    Float_t diamondcovxy[3]={0};
    esdEv->GetDiamondCovXY(diamondcovxy);
    Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,7.*7.};
    AliESDVertex initVertex(pos,cov,1.,1);
    vtFinder->SetVtxStart(&initVertex);
    bkgauss = esdEv->GetMagneticField();
    vtFinder->SetFieldkG(bkgauss);
    currRun = esdEv->GetRunNumber();
    printf("Imposed Vtx constraint for run %d\n",currRun);
    initVertex.Print();
  }
  //
  // reset old vertex info
  if (esdEv->GetPileupVerticesTracks()) esdEv->GetPileupVerticesTracks()->Clear();
  ((AliESDVertex*)esdEv->GetPrimaryVertexTracks())->SetNContributors(-1);
  //
  AliESDVertex *pvtx=vtFinder->FindPrimaryVertex(esdEv);
  if (pvtx) {
    if (pvtx->GetStatus()) {
      esdEv->SetPrimaryVertexTracks(pvtx);
      for (Int_t i=esdEv->GetNumberOfTracks(); i--;) {
	AliESDtrack *t = esdEv->GetTrack(i);
	Double_t x[3]; t->GetXYZ(x);
	t->RelateToVertex(pvtx, bkgauss, kVeryBig);
      }
    }
    delete pvtx;
  }
  else return nullptr;
  //
  return vtFinder->GetVerticesArray();
}
//________________________________________________________________________
Float_t AliESDUtils::GetCorrV0A(Float_t  v0araw, Float_t zv)
{
  // renormalize v0A signal at given Zv to acceptance at Zv=0
  const Double_t pars[] = {0.998864,-0.00407311,-2.47408e-06};
  zv -= pars[0];
  Float_t corr = 1 + zv*(pars[1] + zv*pars[2]);
  return corr>0 ? v0araw/corr : -1;
}
//________________________________________________________________________
Float_t AliESDUtils::GetCorrV0C(Float_t  v0craw, Float_t zv)
{
  // renormalize v0C signal at given Zv to acceptance at Zv=0
  const Double_t pars[] = {1.00083,0.00427623,-2.69047e-05};
  zv -= pars[0];
  Float_t corr = 1 + zv*(pars[1] + zv*pars[2]);
  return corr>0 ? v0craw/corr : -1;
}
//________________________________________________________________________
Float_t AliESDUtils::GetCorrV0A0(Float_t  v0a0raw, Float_t zv)
{
  // renormalize v0A-0 signal at given Zv to acceptance at Zv=0
  const Double_t pars[] = {0.998243,-0.00209013,-6.97686e-06};
  zv -= pars[0];
  Float_t corr = 1 + zv*(pars[1] + zv*pars[2]);
  return corr>0 ? v0a0raw/corr : -1;
}
