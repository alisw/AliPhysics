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
const Double_t AliITSURecoParam::fgkSigmaRoadY                    = 20.E-4;//1000e-4;
const Double_t AliITSURecoParam::fgkSigmaRoadZ                    = 20.E-4;//1000e-4;
const Double_t AliITSURecoParam::fgkTanLorentzAngle               = 0;
const Bool_t   AliITSURecoParam::fgkAllowDiagonalClusterization   = kFALSE;
//
// hardwired params for TPC-ITS border layer
const Double_t AliITSURecoParam::fgkTPCITSWallRMin                = 50.;
const Double_t AliITSURecoParam::fgkTPCITSWallRMax                = 80.;
const Double_t AliITSURecoParam::fgkTPCITSWallZSpanH              = 250.;
const Double_t AliITSURecoParam::fgkTPCITSWallMaxStep             = 6.;
//
const Bool_t   AliITSURecoParam::fgkUseMatLUT[kNTrackingPhases] = {kFALSE,kFALSE,kFALSE};

//
//_____________________________________________________________________________
AliITSURecoParam::AliITSURecoParam()
  :  fNLayers(0)
  ,fMaxDforV0dghtrForProlongation(fgkMaxDforV0dghtrForProlongation)
  ,fMaxDForProlongation(fgkMaxDForProlongation)
  ,fMaxDZForProlongation(fgkMaxDZForProlongation)
  ,fMinPtForProlongation(fgkMinPtForProlongation)
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
  ,fTrackingConditions(0)
  ,fTracker(0)
  ,fSAonly(kFALSE)
  ,fMaxROCycle(126) // like in AliITSUSimulation::kMaxROCycleAccept
{
  // def c-tor
  SetName("ITS");
  SetTitle("ITS");
  for (int i=kNTrackingPhases;i--;) fUseMatLUT[i] = fgkUseMatLUT[i];
}

//_____________________________________________________________________________
AliITSURecoParam::AliITSURecoParam(Int_t nLr)
  :  fNLayers(0)
  ,fMaxDforV0dghtrForProlongation(fgkMaxDforV0dghtrForProlongation)
  ,fMaxDForProlongation(fgkMaxDForProlongation)
  ,fMaxDZForProlongation(fgkMaxDZForProlongation)
  ,fMinPtForProlongation(fgkMinPtForProlongation)
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
  ,fTrackingConditions(0)
  ,fTracker(0)
  ,fSAonly(kFALSE)
  ,fMaxROCycle(126) // like in AliITSUSimulation::kMaxROCycleAccept
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
  delete[] fAllowDiagonalClusterization;
  fTrackingConditions.Delete();
}

//_____________________________________________________________________________
void AliITSURecoParam::SetDefaultSettings(AliITSURecoParam* itsRecoParam) 
{
  // make default reconstruction  parameters for hig  flux env.
  // The settings below are taken from Ruben's MakeITSRecoParam.C
  enum {
    kBit0=0x1<<0, kBit1=0x1<<1, kBit2=0x1<<2, kBit3=0x1<<3,
    kBit4=0x1<<4, kBit5=0x1<<5, kBit6=0x1<<6, kBit7=0x7<<2, kBit8=0x1<<8
  };
  const Bool_t kAllowDiagCl = kFALSE;
  const Bool_t kUseLUT[3]={kTRUE,kTRUE,kFALSE};
  //Use TGeo mat.queries only for RefitInward
  Int_t nLr=7;
  itsRecoParam->SetNLayers(nLr);
  //
  for (int i=0; i<nLr; i++) 
    itsRecoParam->SetAllowDiagonalClusterization(i,kAllowDiagCl);
  for (int i=AliITSURecoParam::kNTrackingPhases; i--;) 
    itsRecoParam->SetUseMatLUT(i,kUseLUT[i]);
  
  // Add tracking conditions >>>
  AliITSUTrackCond *trCond=0;
  {
    int c0nBranch[7] = {3,9,15,4,5,7,10}; // max branching for the seed on layer
    int c0nCands[7]  = {10,15,45,20,60,20,10};// max candidates for the TPC seed
    float c0tr2clChi2[7] ={20,25,30,40,45,45,70};//cut on cluster to track chi2 
    float c0gloChi2[7]   = {6,10,20,30,60,60,70}; //cut on seed global norm chi2
    float c0missPen[7] = {2.,2.,2.,2.,2.,2.,2.};  // missing cluster penalty
    float c0maxChi2SA[14]={0.,0.,0.,0.,2.5,5.,10.,20.,20.,20.,20.,20.,20.,20.};
    // chi2SA vs Nclus
    float c0maxChi2Match = 10.;
    
    trCond = new AliITSUTrackCond();
    trCond->SetNLayers(nLr); 
    trCond->SetMaxITSTPCMatchChi2(c0maxChi2Match);
    //
    for (int i=0; i<nLr; i++) {
      trCond->SetMaxBranches(i,c0nBranch[i]);    // each seed propagated to given layer can produce max nBranch branches
      trCond->SetMaxCandidates(i,c0nCands[i]);   // each tpc track may have at most nCands prolongations
      trCond->SetMaxTr2ClChi2(i,c0tr2clChi2[i]);   // cut on cluster to track chi2
      trCond->SetMaxChi2GloNrm(i,c0gloChi2[i]);  // cut on cluster to track global chi2
      trCond->SetMissPenalty(i,c0missPen[i]);    // missing cluster penalty
    }
    
    for (int i=1; i<=2*nLr; i++) trCond->SetMaxITSSAChi2(i,c0maxChi2SA[i-1]);
    
    trCond->AddNewCondition(5); // min hits
    trCond->AddGroupPattern( kBit0|kBit1|kBit2, 2); // at least 2 hits in 3 inner layers
    trCond->AddGroupPattern( kBit3|kBit4      , 1); // at least 1 hit in 2 middle layers
    trCond->AddGroupPattern( kBit5|kBit6      , 1); // at least 1 hit in 2 outer layers
    //
    trCond->Init();
    //
    itsRecoParam->AddTrackingCondition(trCond);
  }
  
  //-----------------------------------------------------------
  // short tracks
  {
    int c1nBranch[7] = {0,0,0,4,6,6,10}; // max branching for the seed on layer
    int c1nCands[7]  = {0,0,0,5,5,5,8}; // max candidates for the TPC seed
    float c1tr2clChi2[7]= {0,0,0,20,20,20,30}; // cut on cluster to track chi2 
    float c1gloChi2[7]  = {0,0,0,16,40,35,30}; // cut on seed global norm chi2
    float c1missPen[7]  = {0.,0.,0.,2.,2.,2.,2.};    // missing cluster penalty
    float c1maxChi2SA[14]={0.,0.,0.,5.,13.,13.,18.,10.,10.,10.,10.,10.,10.,10.};
    // chi2SA vs Nclus
    float c1maxChi2Match = 10.;
    
    trCond = new AliITSUTrackCond();
    trCond->SetNLayers(nLr); 
    //
    trCond->ExcludeLayer(0);
    trCond->ExcludeLayer(1);
    trCond->ExcludeLayer(2);
    //
    trCond->SetMaxITSTPCMatchChi2(c1maxChi2Match);
    //
    // to exclude some layer use trCon->ExcludeLayer(lrID);
    //
    for (int i=0; i<nLr; i++) {
      trCond->SetMaxBranches(i,c1nBranch[i]);    // each seed propagated to given layer can produce max nBranch branches
      trCond->SetMaxCandidates(i,c1nCands[i]);   // each tpc track may have at most nCands prolongations
      trCond->SetMaxTr2ClChi2(i,c1tr2clChi2[i]); // cut on cluster to track chi2
      trCond->SetMaxChi2GloNrm(i,c1gloChi2[i]);  // cut on cluster to track global chi2
      trCond->SetMissPenalty(i,c1missPen[i]);    // missing cluster penalty
    }
    
    for (int i=1; i<=2*nLr; i++) trCond->SetMaxITSSAChi2(i,c1maxChi2SA[i-1]);
    
    trCond->AddNewCondition(4); // min hits
    trCond->AddGroupPattern( kBit3|kBit4|kBit5|kBit6, 4);
    
    trCond->Init();
    
    itsRecoParam->AddTrackingCondition(trCond);
  }
  
  //-----------------------------------------------------------
  // very short tracks
  {
    int c2nBranch[7] = {0,0,0,0,0,6,10}; // max branching for the seed on layer
    int c2nCands[7]  = {0,0,0,0,0,5,8}; // max candidates for the TPC seed
    float c2tr2clChi2[7]= {0,0,0,0,0,15,20}; // cut on cluster to track chi2
    float c2gloChi2[7]  = {0,0,0,0,0,15,20}; // cut on seed global norm chi2
    float c2missPen[7]  = {0.,0.,0.,0.,0.,2.,2.};    // missing cluster penalty
    float c2maxChi2SA[14]={0.,5.,5.,5.,13.,13.,18.,10.,10.,10.,10.,10.,10.,10.};
    // chi2SA vs Nclus, meaningless for 2 point tracks 
    float c2maxChi2Match = 6.;
    
    trCond = new AliITSUTrackCond();
    trCond->SetNLayers(nLr); 
    //
    trCond->ExcludeLayer(0);
    trCond->ExcludeLayer(1);
    trCond->ExcludeLayer(2);
    trCond->ExcludeLayer(3);
    trCond->ExcludeLayer(4);
    
    trCond->SetMaxITSTPCMatchChi2(c2maxChi2Match);
    //
    // to exclude some layer use trCon->ExcludeLayer(lrID);
    //
    for (int i=0; i<nLr; i++) {
      trCond->SetMaxBranches(i,c2nBranch[i]);    // each seed propagated to given layer can produce max nBranch branches
      trCond->SetMaxCandidates(i,c2nCands[i]);   // each tpc track may have at most nCands prolongations
      trCond->SetMaxTr2ClChi2(i,c2tr2clChi2[i]); // cut on cluster to track chi2
      trCond->SetMaxChi2GloNrm(i,c2gloChi2[i]);  // cut on cluster to track global chi2
      trCond->SetMissPenalty(i,c2missPen[i]);    // missing cluster penalty
    }
    
    for (int i=1;i<=2*nLr;i++) trCond->SetMaxITSSAChi2(i,c2maxChi2SA[i-1]);
    //
    trCond->AddNewCondition(2); // min hits
    trCond->AddGroupPattern( kBit5|kBit6, 2);
    //
    trCond->Init();
    //
    itsRecoParam->AddTrackingCondition(trCond);
  }
  //
}

//_____________________________________________________________________________
AliITSURecoParam *AliITSURecoParam::GetHighFluxParam(Bool_t init) 
{
  // make default reconstruction  parameters for low  flux env.
  AliITSURecoParam *param = new AliITSURecoParam();
  param->SetEventSpecie(AliRecoParam::kHighMult);
  param->SetTitle("HighMult");
  if (init) SetDefaultSettings(param);
  // put here params
  return param;
}

//_____________________________________________________________________________
AliITSURecoParam *AliITSURecoParam::GetLowFluxParam(Bool_t init) 
{
  // make default reconstruction  parameters for low  flux env.
  AliITSURecoParam *param = new AliITSURecoParam();
  param->SetEventSpecie(AliRecoParam::kLowMult);
  param->SetTitle("LowMult");
  if (init) SetDefaultSettings(param);
  // put here params
  return param;
}

//_____________________________________________________________________________
AliITSURecoParam *AliITSURecoParam::GetCosmicTestParam(Bool_t init) 
{
  // make default reconstruction  parameters for cosmics
  AliITSURecoParam *param = new AliITSURecoParam();
  param->SetEventSpecie(AliRecoParam::kCosmic);
  param->SetTitle("Cosmic");
  if (init) SetDefaultSettings(param);
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
  fAllowDiagonalClusterization = new Bool_t[n];
  //
  for (int i=n;i--;) {
    fAllowDiagonalClusterization[i] = fgkAllowDiagonalClusterization;
    fTanLorentzAngle[i] = fgkTanLorentzAngle;
    fSigmaY2[i] = fgkSigmaRoadY*fgkSigmaRoadY;
    fSigmaZ2[i] = fgkSigmaRoadZ*fgkSigmaRoadZ;
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
  //
  printf("Use material LUT at steps: "); 
  for (int i=0;i<kNTrackingPhases;i++) printf("%d : %s |",i,GetUseMatLUT(i) ? "ON":"OFF"); printf("\n");
  printf("TPC-ITS wall: %.3f<R<%.3f DZ/2=%.3f MaxStep=%.3f\n",
	 fTPCITSWallRMin,fTPCITSWallRMax,fTPCITSWallZSpanH,fTPCITSWallMaxStep);
  //
  printf("N.Layers: %d\n",fNLayers);
  if (fNLayers>0) {
    printf("Layer-wise data:\n");
    printf("%-30s\t:","fTanLorentzAngle");  for (int i=0;i<fNLayers;i++) printf(" %+.2e",fTanLorentzAngle[i]); printf("\n");
    printf("%-30s\t:","fSigmaY2");          for (int i=0;i<fNLayers;i++) printf(" %+.2e",fSigmaY2[i]); printf("\n");
    printf("%-30s\t:","fSigmaZ2");          for (int i=0;i<fNLayers;i++) printf(" %+.2e",fSigmaZ2[i]); printf("\n");
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
