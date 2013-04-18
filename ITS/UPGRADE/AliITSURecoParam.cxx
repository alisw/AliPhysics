/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliITSURecoParam.h"
#include "AliLog.h"
#include "AliITSUTrackCond.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with ITS reconstruction parameters                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

ClassImp(AliITSURecoParam)


const Double_t AliITSURecoParam::fgkMaxDforV0dghtrForProlongation = 30;
const Double_t AliITSURecoParam::fgkMaxDForProlongation           = 40; 
const Double_t AliITSURecoParam::fgkMaxDZForProlongation          = 60;      
const Double_t AliITSURecoParam::fgkMinPtForProlongation          = 0.0; 
const Double_t AliITSURecoParam::fgkNSigmaRoadY                   = 5.;
const Double_t AliITSURecoParam::fgkNSigmaRoadZ                   = 5.; 
const Double_t AliITSURecoParam::fgkSigmaRoadY                    = 100e-4;//1000e-4;
const Double_t AliITSURecoParam::fgkSigmaRoadZ                    = 100e-4;//1000e-4;
const Double_t AliITSURecoParam::fgkMaxTr2ClChi2                  = 15.;
const Double_t AliITSURecoParam::fgkTanLorentzAngle               = 0;
const Double_t AliITSURecoParam::fgkMissPenalty                   = 2.0;
const Bool_t   AliITSURecoParam::fgkAllowDiagonalClusterization   = kFALSE;
//
// hardwired params for TPC-ITS border layer
const Double_t AliITSURecoParam::fgkTPCITSWallRMin                = 50.;
const Double_t AliITSURecoParam::fgkTPCITSWallRMax                = 80.;
const Double_t AliITSURecoParam::fgkTPCITSWallZSpanH              = 250.;
const Double_t AliITSURecoParam::fgkTPCITSWallMaxStep             = 6.;


//
//_____________________________________________________________________________
AliITSURecoParam::AliITSURecoParam()
  :  fNLayers(0)
  ,fMaxDforV0dghtrForProlongation(fgkMaxDforV0dghtrForProlongation)
  ,fMaxDForProlongation(fgkMaxDForProlongation)
  ,fMaxDZForProlongation(fgkMaxDZForProlongation)
  ,fMinPtForProlongation(fgkMinPtForProlongation)
  ,fNSigmaRoadY(fgkNSigmaRoadY)
  ,fNSigmaRoadZ(fgkNSigmaRoadZ)
     //
  ,fTPCITSWallRMin(fgkTPCITSWallRMin)
  ,fTPCITSWallRMax(fgkTPCITSWallRMax)
  ,fTPCITSWallZSpanH(fgkTPCITSWallZSpanH)
  ,fTPCITSWallMaxStep(fgkTPCITSWallMaxStep)
     //
  ,fAllowDiagonalClusterization(0)
  ,fTanLorentzAngle(0)
  ,fSigmaY2(0)
  ,fSigmaZ2(0)
  ,fMaxTr2ClChi2(0)
  ,fMissPenalty(0)
  ,fTrackingConditions(0)
{
  // def c-tor
  SetName("ITS");
  SetTitle("ITS");
}

//_____________________________________________________________________________
AliITSURecoParam::AliITSURecoParam(Int_t nLr)
  :  fNLayers(0)
  ,fMaxDforV0dghtrForProlongation(fgkMaxDforV0dghtrForProlongation)
  ,fMaxDForProlongation(fgkMaxDForProlongation)
  ,fMaxDZForProlongation(fgkMaxDZForProlongation)
  ,fMinPtForProlongation(fgkMinPtForProlongation)
  ,fNSigmaRoadY(fgkNSigmaRoadY)
  ,fNSigmaRoadZ(fgkNSigmaRoadZ)
     //
  ,fTPCITSWallRMin(fgkTPCITSWallRMin)
  ,fTPCITSWallRMax(fgkTPCITSWallRMax)
  ,fTPCITSWallZSpanH(fgkTPCITSWallZSpanH)
  ,fTPCITSWallMaxStep(fgkTPCITSWallMaxStep)
     //
  ,fAllowDiagonalClusterization(0)
  ,fTanLorentzAngle(0)
  ,fSigmaY2(0)
  ,fSigmaZ2(0)
  ,fMaxTr2ClChi2(0)
  ,fMissPenalty(0)
  ,fTrackingConditions(0)
{
  // def c-tor
  SetName("ITS");
  SetTitle("ITS");
  SetNLayers(nLr);
}

//_____________________________________________________________________________
AliITSURecoParam::~AliITSURecoParam() 
{
  // destructor
  delete[] fTanLorentzAngle;
  delete[] fSigmaY2;
  delete[] fSigmaZ2;
  delete[] fMaxTr2ClChi2;
  delete[] fMissPenalty;
  delete[] fAllowDiagonalClusterization;
  fTrackingConditions.Delete();
}

//_____________________________________________________________________________
AliITSURecoParam *AliITSURecoParam::GetHighFluxParam() 
{
  // make default reconstruction  parameters for hig  flux env.
  AliITSURecoParam *param = new AliITSURecoParam(); 
  //
  // put here params
  return param;
}

//_____________________________________________________________________________
AliITSURecoParam *AliITSURecoParam::GetLowFluxParam() 
{
  // make default reconstruction  parameters for low  flux env.
  AliITSURecoParam *param = new AliITSURecoParam();
  // put here params
  return param;
}

//_____________________________________________________________________________
AliITSURecoParam *AliITSURecoParam::GetCosmicTestParam() 
{
  // make default reconstruction  parameters for cosmics
  AliITSURecoParam *param = new AliITSURecoParam();
  // put here params
  return param;
}

//_____________________________________________________________________________
void  AliITSURecoParam::SetNLayers(Int_t n)
{
  // set n layers and init all layer dependent arrays
  if (fNLayers>0) AliFatal(Form("Number of layers was already set to %d",fNLayers));
  if (n<1) n = 1; // in case we want to have dummy params
  fNLayers = n;
  fTanLorentzAngle = new Double_t[n];
  fSigmaY2 = new Double_t[n];
  fSigmaZ2 = new Double_t[n];
  fMaxTr2ClChi2 = new Double_t[n];
  fMissPenalty  = new Double_t[n];
  fAllowDiagonalClusterization = new Bool_t[n];
  //
  for (int i=n;i--;) {
    fAllowDiagonalClusterization[i] = fgkAllowDiagonalClusterization;
    fTanLorentzAngle[i] = fgkTanLorentzAngle;
    fSigmaY2[i] = fgkSigmaRoadY*fgkSigmaRoadY;
    fSigmaZ2[i] = fgkSigmaRoadZ*fgkSigmaRoadZ;
    fMaxTr2ClChi2[i] = fgkMaxTr2ClChi2;
    fMissPenalty[i]  = fgkMissPenalty;
  }
  //
}

//_____________________________________________________________________________
void  AliITSURecoParam::SetTanLorentzAngle(Int_t lr, Double_t v)
{
  // set Lorentz angle value
  if (lr>=fNLayers) AliFatal(Form("Number of defined layers is %d",fNLayers));
  fTanLorentzAngle[lr] = v;
}

//_____________________________________________________________________________
void  AliITSURecoParam::SetSigmaY2(Int_t lr, Double_t v)
{
  // set Lorentz angle value
  if (lr>=fNLayers) AliFatal(Form("Number of defined layers is %d",fNLayers));
  fSigmaY2[lr] = v;
}

//_____________________________________________________________________________
void  AliITSURecoParam::SetSigmaZ2(Int_t lr, Double_t v)
{
  // set Lorentz angle value
  if (lr>=fNLayers) AliFatal(Form("Number of defined layers is %d",fNLayers));
  fSigmaZ2[lr] = v;
}

//_____________________________________________________________________________
void  AliITSURecoParam::SetMaxTr2ClChi2(Int_t lr, Double_t v)
{
  // set Lorentz angle value
  if (lr>=fNLayers) AliFatal(Form("Number of defined layers is %d",fNLayers));
  fMaxTr2ClChi2[lr] = v;
}

//_____________________________________________________________________________
void  AliITSURecoParam::SetMissPenalty(Int_t lr, Double_t v)
{
  // set Lorentz angle value
  if (lr>=fNLayers) AliFatal(Form("Number of defined layers is %d",fNLayers));
  fMissPenalty[lr] = v;
}

//_____________________________________________________________________________
void  AliITSURecoParam::SetAllowDiagonalClusterization(Int_t lr, Bool_t v)
{
  // set Lorentz angle value
  if (lr>=fNLayers) AliFatal(Form("Number of defined layers is %d",fNLayers));
  fAllowDiagonalClusterization[lr] = v;
}

//========================================================================
//_____________________________________________________________________________
void AliITSURecoParam::Print(Option_t *) const
{
  // print params
  printf("%s: %s %s\n",ClassName(),GetName(),GetTitle());
  printf("%-30s\t%f\n","fMaxDforV0dghtrForProlongation",fMaxDforV0dghtrForProlongation);
  printf("%-30s\t%f\n","fMaxDForProlongation",fMaxDForProlongation); 
  printf("%-30s\t%f\n","fMaxDZForProlongation",fMaxDZForProlongation);
  printf("%-30s\t%f\n","fMinPtForProlongation",fMinPtForProlongation);
  printf("%-30s\t%f\n","fNSigmaRoadY",fNSigmaRoadY);
  printf("%-30s\t%f\n","fNSigmaRoadZ",fNSigmaRoadZ);
  //
  printf("TPC-ITS wall: %.3f<R<%.3f DZ/2=%.3f MaxStep=%.3f\n",
	 fTPCITSWallRMin,fTPCITSWallRMax,fTPCITSWallZSpanH,fTPCITSWallMaxStep);
  //
  printf("N.Layers: %d\n",fNLayers);
  if (fNLayers>0) {
    printf("Layer-wise data:\n");
    printf("%-30s\t:","fTanLorentzAngle");  for (int i=0;i<fNLayers;i++) printf(" %+.2e",fTanLorentzAngle[i]); printf("\n");
    printf("%-30s\t:","fSigmaY2");          for (int i=0;i<fNLayers;i++) printf(" %+.2e",fSigmaY2[i]); printf("\n");
    printf("%-30s\t:","fSigmaZ2");          for (int i=0;i<fNLayers;i++) printf(" %+.2e",fSigmaZ2[i]); printf("\n");
    printf("%-30s\t:","fMaxTr2ClChi2");     for (int i=0;i<fNLayers;i++) printf(" %+.2e",fMaxTr2ClChi2[i]); printf("\n");
    printf("%-30s\t:","fMissPenalty");      for (int i=0;i<fNLayers;i++) printf(" %+.2e",fMissPenalty[i]); printf("\n");
  }
  //
  int nTrCond = GetNTrackingConditions();
  printf("%d tracking conditions defined\n",nTrCond);
  for (int itc=0;itc<nTrCond;itc++) {
    printf("Tracking condition %d\n",itc);
    GetTrackingCondition(itc)->Print();
  }
}

//__________________________________________________
void  AliITSURecoParam::AddTrackingCondition(AliITSUTrackCond* cond)
{
  // Add new tracking condition
  fTrackingConditions.AddLast(cond);
}
