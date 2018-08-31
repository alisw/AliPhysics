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

#include <TVector3.h>
#include <TDatabasePDG.h>
#include <TPDGCode.h>
#include <TParticlePDG.h>

#include "AliAODv0.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliAODEvent.h"
#include "AliESDtrack.h"
#include "AliVEvent.h"
#include "AliLog.h"
#include "AliKFParticle.h"
#include "AliVTrack.h"
#include "AliKFVertex.h"

#include "AliAODv0KineCuts.h"

ClassImp(AliAODv0KineCuts)

/// \file AliAODv0KineCuts.cxx

AliAODv0KineCuts::AliAODv0KineCuts() :
  fEvent(0x0)
  , fPrimaryVertex(0x0)
  , fType(0)
  , fMode(0)
  , fTPCNcls(30)
  , fTPCrefit(kTRUE)
  , fTPCclsRatio(0.6)
  , fNoKinks(kTRUE)
  , fKinkMotherList(1000)
  , fNumberKinkMothers(0)
  , fGcutChi2NDF(10)
  , fGcutInvMass(0.05)
  , fK0cutChi2NDF(10)
  , fLcutChi2NDF(10)
  , fUseExternalVertex(kFALSE)
  , fDeleteVertex(kFALSE)
{
  //
  // Default constructor
  //

  // default single track cuts
  fTPCNcls = 30;                // minimal number of the TPC clusters
  fTPCrefit = kTRUE;           // TPC refit
  fTPCclsRatio = 0.6;          // minimal foun/findable TPC cluster ratio
  fNoKinks = kTRUE;            // kinks - no [kTRUE] or do not care [kFalse]


  // default gamma cuts values
  fGcutChi2NDF = 10;           // Chi2NF cut value for the AliKFparticle gamma
  fGcutCosPoint[0] = 0;        // cos of the pointing angle [min, max]
  fGcutCosPoint[1] = 0.02;     // cos of the pointing angle [min, max]
  fGcutDCA[0] = 0.;            // DCA between the daughter tracks [min, max]
  fGcutDCA[1] = 0.25;          // DCA between the daughter tracks [min, max]
  fGcutVertexR[0] = 3.;        // radius of the conversion point [min, max]
  fGcutVertexR[1] = 90.;       // radius of the conversion point [min, max]
  fGcutPsiPair[0] = 0.;        // value of the psi pair cut [min, max]
  fGcutPsiPair[1] = 0.05;      // value of the psi pair cut [min, max]
  fGcutInvMass = 0.05;         // upper value on the gamma invariant mass
  // default K0 cuts
  fK0cutChi2NDF = 10;          // Chi2NF cut value for the AliKFparticle K0
  fK0cutCosPoint[0] = 0.;      // cos of the pointing angle [min, max]
  fK0cutCosPoint[1] = 0.02;    // cos of the pointing angle [min, max]
  fK0cutDCA[0] = 0.;           // DCA between the daughter tracks [min, max]
  fK0cutDCA[1] = 0.2;          // DCA between the daughter tracks [min, max]
  fK0cutVertexR[0] = 2.0;      // radius of the decay point [min, max]
  fK0cutVertexR[1] = 30.0;     // radius of the decay point [min, max]
  fK0cutInvMass[0] = 0.486;    // invariant mass window
  fK0cutInvMass[1] = 0.508;    // invariant mass window
  // Lambda & anti-Lambda cut values
  fLcutChi2NDF = 10;           // Chi2NF cut value for the AliKFparticle K0
  fLcutCosPoint[0] = 0.;       // cos of the pointing angle [min, max]
  fLcutCosPoint[1] = 0.02;     // cos of the pointing angle [min, max]
  fLcutDCA[0] = 0.;            // DCA between the daughter tracks [min, max]
  fLcutDCA[1] = 0.2;           // DCA between the daughter tracks [min, max]
  fLcutVertexR[0] = 2.0;       // radius of the decay point [min, max]
  fLcutVertexR[1] = 40.0;      // radius of the decay point [min, max]
  fLcutInvMass[0] = 1.11;      // invariant mass window
  fLcutInvMass[1] = 1.12;      // invariant mass window
    
}
//____________________________________________________________________
AliAODv0KineCuts::~AliAODv0KineCuts(){
  /// Destructor


}
//____________________________________________________________________
AliAODv0KineCuts::AliAODv0KineCuts(const AliAODv0KineCuts &ref):
  TObject(ref)
  , fEvent(0x0)
  , fPrimaryVertex(0x0)
  , fType(0)
  , fMode(0)
  , fTPCNcls(30)
  , fTPCrefit(kTRUE)
  , fTPCclsRatio(0.6)
  , fNoKinks(kTRUE)
  , fKinkMotherList(ref.fKinkMotherList)
  , fNumberKinkMothers(ref.fNumberKinkMothers)
  , fGcutChi2NDF(10)
  , fGcutInvMass(0.05)
  , fK0cutChi2NDF(10)
  , fLcutChi2NDF(10)
  , fUseExternalVertex(kFALSE)
  , fDeleteVertex(kFALSE)
{
  /// Copy operator

  ref.Copy(*this);
}
//____________________________________________________________________
AliAODv0KineCuts &AliAODv0KineCuts::operator=(const AliAODv0KineCuts &ref){
  /// assignment operator

  if(this != &ref)
    ref.Copy(*this);
  return *this; 
}
//____________________________________________________________________
void AliAODv0KineCuts::Copy(TObject &ref) const {
  /// Performs the copying of the object

  TObject::Copy(ref);

  AliAODv0KineCuts &target = dynamic_cast<AliAODv0KineCuts &>(ref);

  // default single track cuts
  target.fTPCNcls = fTPCNcls;
  target.fTPCrefit = fTPCrefit;
  target.fTPCclsRatio = fTPCclsRatio;
  target.fNoKinks = fNoKinks;
  target.fKinkMotherList = fKinkMotherList;
  target.fNumberKinkMothers = fNumberKinkMothers;
  target.fUseExternalVertex = fUseExternalVertex;  //added december 2nd 2011
  target.fDeleteVertex = fDeleteVertex;  //added december 2nd 2011

  // default gamma cuts values
  target.fGcutChi2NDF = fGcutChi2NDF;
  memcpy(target.fGcutCosPoint, fGcutCosPoint, sizeof(Float_t) * 2);
  memcpy(target.fGcutDCA, fGcutDCA, sizeof(Float_t) * 2); 
  memcpy(target.fGcutVertexR, fGcutVertexR, sizeof(Float_t) * 2);
  memcpy(target.fGcutPsiPair, fGcutPsiPair, sizeof(Float_t) * 2);
  target.fGcutInvMass = fGcutInvMass;
  // default K0 cuts
  target.fK0cutChi2NDF = fK0cutChi2NDF;
  memcpy(target.fK0cutCosPoint, fK0cutCosPoint, sizeof(Float_t) * 2);
  memcpy(target.fK0cutDCA, fK0cutDCA, sizeof(Float_t) * 2);
  memcpy(target.fK0cutVertexR, fK0cutVertexR, sizeof(Float_t) * 2);
  memcpy(target.fK0cutInvMass, fK0cutInvMass, sizeof(Float_t) * 2);
  // Lambda & anti-Lambda cut values
  target.fLcutChi2NDF = fLcutChi2NDF;
  memcpy(target.fLcutCosPoint, fLcutCosPoint, sizeof(Float_t) * 2);
  memcpy(target.fLcutDCA, fLcutDCA, sizeof(Float_t) * 2);
  memcpy(target.fLcutVertexR, fLcutVertexR, sizeof(Float_t) * 2);
  memcpy(target.fLcutInvMass, fLcutInvMass, sizeof(Float_t) * 2);
  
}
//____________________________________________________________________
Bool_t  AliAODv0KineCuts::ProcessV0(AliAODv0* const v0, Int_t &pdgV0, Int_t &pdgP, Int_t &pdgN) const 
{
  /// main user function

  if(!v0) return kFALSE;
  if(!fEvent){
    AliErrorClass("No valid Event pointer available, provide it first");
    return kFALSE;
  }

  if(!V0CutsCommon(v0)) return kFALSE;

  const Int_t id = PreselectV0(v0);

  if(!SingleTrackCuts(v0)) return kFALSE;

  switch(id){
  case kUndef:
    return kFALSE;
  case kGamma:
    return CaseGamma(v0, pdgV0, pdgP, pdgN);
  case kK0:
    return CaseK0(v0, pdgV0, pdgP, pdgN);
  case kLambda:
    return CaseLambda(v0, pdgV0, pdgP, pdgN, 0);
  case kALambda:
    return CaseLambda(v0, pdgV0, pdgP, pdgN, 1);
  default:
    return kFALSE; 
  }

  return kFALSE;
}
//____________________________________________________________________
Bool_t  AliAODv0KineCuts::ProcessV0(AliAODv0* const v0, Int_t &pdgP, Int_t &pdgN) const 
{
  /// main user function, simplified if the V0 identity is not necessary

  if(!v0) return kFALSE;
  if(!fEvent){
    AliErrorClass("No valid Event pointer available, provide it first");
    return kFALSE;
  }

  Int_t idV0 = -1;
  return ProcessV0(v0, idV0, pdgP, pdgN);

}
//____________________________________________________________________
Int_t AliAODv0KineCuts::PreselectV0(AliAODv0* const v0) const 
{
  /// Make a preselection (exclusive) of the V0 cadidates based on
  /// Armenteros plot
  /// the armenteros cut values are currently fixed and user is not able to set them via
  /// set funcions. The reason is that these cuts are optimized and furneter changes should
  /// not be necessary. To prove otherwise please study in detail before changing the values

  // for clarity
  const Float_t alpha =v0->AlphaV0(); //ap[0];
  const Float_t qt = v0->PtArmV0();//ap[1];

  // selection cuts 
  // - the reagions for different candidates must not overlap 

  // Gamma cuts
  const Double_t cutAlphaG = 0.35; 
  const Double_t cutQTG = 0.05;
  const Double_t cutAlphaG2[2] = {0.6, 0.8};
  const Double_t cutQTG2 = 0.04;

  // K0 cuts
  const Float_t cutQTK0[2] = {0.1075, 0.215};
  const Float_t cutAPK0[2] = {0.199, 0.8};   // parameters for curved QT cut
  
  // Lambda & A-Lambda cuts
  const Float_t cutQTL = 0.03;
  const Float_t cutAlphaL[2] = {0.35, 0.7};
  const Float_t cutAlphaAL[2] = {-0.7,  -0.35};
  const Float_t cutAPL[3] = {0.107, -0.69, 0.5};  // parameters fir curved QT cut


  if(kPurity == fMode){
  // Check for Gamma candidates
    if(qt < cutQTG){
      if( (TMath::Abs(alpha) < cutAlphaG) ) return kGamma;
    }
    // additional region - should help high pT gammas
    if(qt < cutQTG2){
      if( (TMath::Abs(alpha) > cutAlphaG2[0]) &&  (TMath::Abs(alpha) < cutAlphaG2[1]) ) return kGamma;
    }
  }
  if(kEffGamma == fMode){
    if(qt < cutQTG) return kGamma;
  }

  
  // Check for K0 candidates
  Float_t q = cutAPK0[0] * TMath::Sqrt(TMath::Abs(1 - alpha*alpha/(cutAPK0[1]*cutAPK0[1])));
  if( (qt > cutQTK0[0]) && (qt < cutQTK0[1]) && (qt > q) ){
    return kK0;
  }

  // Check for Lambda candidates
  q = cutAPL[0] * TMath::Sqrt(TMath::Abs(1 - ( (alpha + cutAPL[1]) * (alpha + cutAPL[1]) ) / (cutAPL[2]*cutAPL[2]) ));
  if( (alpha > cutAlphaL[0]) && (alpha < cutAlphaL[1]) && (qt > cutQTL) && (qt < q)  ){
    return kLambda;
  }

  // Check for A-Lambda candidates
  q = cutAPL[0] * TMath::Sqrt(TMath::Abs(1 - ( (alpha - cutAPL[1]) * (alpha - cutAPL[1]) ) / (cutAPL[2]*cutAPL[2]) ));
  if( (alpha > cutAlphaAL[0]) && (alpha < cutAlphaAL[1]) && (qt > cutQTL) && (qt < q)  ){
    return kALambda;
  }
  
  return kUndef;
}
//____________________________________________________________________
Bool_t  AliAODv0KineCuts::SingleTrackCuts(AliAODv0 * const v0) const 
{
  /// apply single track cuts
  /// correct sign not relevat here

  if(!v0) return kFALSE;
  
  AliAODTrack* d[2] = {
    dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(0)),
    dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(1))
  };
  
  for(Int_t i=0; i<2; ++i){
    if(!d[i]) return kFALSE;
    
    // status word
    ULong64_t status = d[i]->GetStatus();

    // No. of TPC clusters leave to the users
    if(d[i]->GetTPCNcls() < 1) return kFALSE;

    // TPC refit
    if(!(status & AliAODTrack::kTPCrefit)) return kFALSE;
  
    // TPC cluster ratio
    Float_t cRatioTPC = d[i]->GetTPCNclsF() > 0. ? static_cast<Float_t>(d[i]->GetTPCNcls())/static_cast<Float_t> (d[i]->GetTPCNclsF()) : 1.;
    if(cRatioTPC < 0.6) return kFALSE;
    
    // kinks
    if(fNoKinks && (IsKinkDaughter(d[i]) || IsKinkMother(d[i]))) return kFALSE;
    
  }

  return kTRUE;
}
//____________________________________________________________________
Bool_t AliAODv0KineCuts::CaseGamma(AliAODv0* const v0, Int_t &pdgV0, Int_t &pdgP, Int_t &pdgN) const 
{
  /// process the gamma conversion candidate

  if(!v0) return kFALSE;

  AliVTrack* daughter[2];
  Int_t pIndex = 0, nIndex = 0;

  Bool_t sign = CheckSigns(v0);
  if(sign){
    pIndex = 0;
    nIndex = 1;
  }
  else{
    pIndex = 1;
    nIndex = 0;    
  }
  daughter[0] = dynamic_cast<AliVTrack *>(v0->GetSecondaryVtx()->GetDaughter(pIndex));
  daughter[1] = dynamic_cast<AliVTrack *>(v0->GetSecondaryVtx()->GetDaughter(nIndex));
  if(!daughter[0] || !daughter[1]) return kFALSE;

  AliKFParticle *kfMother = CreateMotherParticle(daughter[0], daughter[1], TMath::Abs(kElectron), TMath::Abs(kElectron));
  if(!kfMother) return kFALSE;

//  AliAODTrack* d[2];
//  d[0] = dynamic_cast<AliAODTrack*>(fEvent->GetTrack(pIndex));
//  d[1] = dynamic_cast<AliAODTrack*>(fEvent->GetTrack(nIndex));

  Float_t iMass = v0->InvMass2Prongs(0,1,11,11);

  // cos pointing angle
  Double_t cosPoint = v0->CosPointingAngle(dynamic_cast<AliAODVertex *>(fEvent->GetPrimaryVertex()));
  cosPoint = TMath::ACos(cosPoint);

  // DCA between daughters
  Double_t dca = v0->DcaV0Daughters();

  // Production vertex
  Double_t xyz[3]; 
  v0->GetXYZ(xyz);
  Double_t r = TMath::Sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);

  Double_t xy[2];
  Double_t r2 = -1.;
  if ( GetConvPosXY(static_cast<AliAODTrack *>(daughter[0]), static_cast<AliAODTrack*>(daughter[1]), xy) ){
    r2 = TMath::Sqrt(xy[0]*xy[0] + xy[1]*xy[1]);
  }

  // psi pair 
  Double_t psiPair = PsiPair(v0);
  
  // V0 chi2/ndf
  Double_t chi2ndf = kfMother->GetChi2()/kfMother->GetNDF();

  if(kfMother) delete kfMother; 
  
  // apply the cuts

  if(iMass > fGcutInvMass) return kFALSE;

  if(chi2ndf > fGcutChi2NDF) return kFALSE;

  if(cosPoint < fGcutCosPoint[0] || cosPoint > fGcutCosPoint[1]) return kFALSE;

  if(dca < fGcutDCA[0] || dca > fGcutDCA[1]) return kFALSE;

  if(r < fGcutVertexR[0] || r > fGcutVertexR[1]) return kFALSE;

  if(psiPair < fGcutPsiPair[0] || psiPair > fGcutPsiPair[1]) return kFALSE;
  
  // all cuts passed

  pdgV0 = 22;
  if(sign){
    pdgP = -11;
    pdgN = 11;
  }
  else{
    pdgP = 11;
    pdgN = -11;
  }

  return kTRUE;
}
//____________________________________________________________________
Bool_t  AliAODv0KineCuts::CaseK0(AliAODv0* const v0, Int_t &pdgV0, Int_t &pdgP, Int_t &pdgN) const {
  /// process the K0 candidate

  if(!v0) return kFALSE;
  
  AliVTrack* daughter[2];
  Int_t pIndex = 0, nIndex = 0;
  Bool_t sign = CheckSigns(v0);
  if(sign){
    pIndex = 0;
    nIndex = 1;
  }
  else{
    pIndex = 1;
    nIndex = 0;    
  }
  daughter[0] = dynamic_cast<AliVTrack *>(v0->GetSecondaryVtx()->GetDaughter(pIndex));
  daughter[1] = dynamic_cast<AliVTrack *>(v0->GetSecondaryVtx()->GetDaughter(nIndex));

  if(!daughter[0] || !daughter[1]) return kFALSE;

  AliKFParticle *kfMother = CreateMotherParticle(daughter[0], daughter[1], TMath::Abs(kPiPlus), TMath::Abs(kPiPlus));
  if(!kfMother) return kFALSE;

//  AliAODTrack* d[2];
//  d[0] = dynamic_cast<AliAODTrack*>(fEvent->GetTrack(pIndex));
//  d[1] = dynamic_cast<AliAODTrack*>(fEvent->GetTrack(nIndex));

  Float_t iMass = v0->MassK0Short();

  // cos pointing angle
  Double_t cosPoint = v0->CosPointingAngle(dynamic_cast<AliAODVertex *>(fEvent->GetPrimaryVertex()));
  cosPoint = TMath::ACos(cosPoint);

  // DCA between daughters
  Double_t dca = v0->DcaV0Daughters();

  // Production vertex
  Double_t xyz[3]; 
  v0->GetXYZ(xyz);

  Double_t r = TMath::Sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);  

  // V0 chi2/ndf
  Double_t chi2ndf = kfMother->GetChi2()/kfMother->GetNDF();
  
  if(kfMother) delete kfMother; 

  //
  // apply the cuts
  //
  if(iMass < fK0cutInvMass[0] || iMass > fK0cutInvMass[1]) return kFALSE;

  if(chi2ndf > fK0cutChi2NDF) return kFALSE;

  if(cosPoint < fK0cutCosPoint[0] || cosPoint > fK0cutCosPoint[1]) return kFALSE;

  if(dca < fK0cutDCA[0] || dca > fK0cutDCA[1]) return kFALSE;

  if(r < fK0cutVertexR[0] || r > fK0cutVertexR[1]) return kFALSE;

  // all cuts passed
  pdgV0 = 310;
  if(sign){
    pdgP = 211;
    pdgN = -211;
  }
  else{
    pdgP = -211;
    pdgN = 211;
  }

  return kTRUE;
}
//____________________________________________________________________
Bool_t  AliAODv0KineCuts::CaseLambda(AliAODv0* const v0, Int_t &pdgV0, Int_t &pdgP, Int_t &pdgN, Int_t id) const {
  /// process teh Lambda and Anti-Lambda candidate

  if(!v0) return kFALSE;

    const Double_t cL0mass=TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();  // PDG lambda mass

  AliVTrack* daughter[2];
  Int_t pIndex = 0, nIndex = 0;
  Float_t mMass[2] = {-1., -1.};
  Bool_t sign = CheckSigns(v0);
  if(sign){
    pIndex = 0;
    nIndex = 1;
    mMass[0] = v0->MassLambda();
    mMass[1] = v0->MassAntiLambda();
  }
  else{
    pIndex = 1;
    nIndex = 0;    
    mMass[0] = v0->MassAntiLambda();
    mMass[1] = v0->MassLambda();
  }
  daughter[0] = dynamic_cast<AliVTrack *>(v0->GetSecondaryVtx()->GetDaughter(pIndex));
  daughter[1] = dynamic_cast<AliVTrack *>(v0->GetSecondaryVtx()->GetDaughter(nIndex));

  if(!daughter[0] || !daughter[1]) return kFALSE;

  AliKFParticle *kfMother[2] = {0x0, 0x0};
  // Lambda
  kfMother[0] = CreateMotherParticle(daughter[0], daughter[1], TMath::Abs(kProton), TMath::Abs(kPiPlus));
  if(!kfMother[0]) return kFALSE;
  
  // Anti-Lambda
  kfMother[1] = CreateMotherParticle(daughter[0], daughter[1], TMath::Abs(kPiPlus), TMath::Abs(kProton));
  if(!kfMother[1]) return kFALSE;

  Float_t dMass[2] = {static_cast<Float_t>(TMath::Abs(mMass[0] - cL0mass)), static_cast<Float_t>(TMath::Abs(mMass[1] - cL0mass))};
  
  Float_t p[2] = {static_cast<Float_t>(daughter[0]->P()), static_cast<Float_t>(daughter[1]->P())}; 

  // check the 3 lambda - antilambda variables
  Int_t check[2] = {-1, -1};   // 0 : lambda, 1 : antilambda
  // 1) momentum of the daughter particles - proton is expected to have higher momentum than pion
  check[0] = (p[0] > p[1]) ? 0 : 1;
  // 2) mass of the mother particle
  check[1] = (dMass[0] < dMass[1]) ? 0 : 1;
 
  // require positive correlation of (1) and (2)
  if(check[0] != check[1]){
    if(kfMother[0]) delete kfMother[0]; 
    if(kfMother[1]) delete kfMother[1]; 
    return kFALSE;
  }

  // now that the check[0] == check[1]
  const Int_t type = check[0];

  // require that the input armenteros preselection agree:
  if(type != id) return kFALSE;

  Float_t iMass =0.;
  if(sign){
    iMass = (type == 0) ? v0->MassLambda() : v0->MassAntiLambda();
  } else{
    iMass = (type == 0) ? v0->MassAntiLambda() : v0->MassLambda();
  }

  // cos pointing angle
  Double_t cosPoint = v0->CosPointingAngle(dynamic_cast<AliAODVertex *>(fEvent->GetPrimaryVertex()));
  cosPoint = TMath::ACos(cosPoint);

  // DCA between daughters
  Double_t dca = v0->DcaV0Daughters();
  
  // Production vertex
  Double_t xyz[3]; 
  v0->GetXYZ(xyz);
  Double_t r = TMath::Sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);

  // proton - pion indices
  Int_t ix[2] = {0, 1};
  if(1 == type){
    ix[0] = 1;
    ix[1] = 0;
  }

  // V0 chi2/ndf
  Double_t chi2ndf = kfMother[type]->GetChi2()/kfMother[type]->GetNDF();

  if(kfMother[0]) delete kfMother[0]; 
  if(kfMother[1]) delete kfMother[1]; 

  //
  // apply the cuts
  //

  if(iMass < fLcutInvMass[0] || iMass > fLcutInvMass[1]) return kFALSE;

  if(chi2ndf > fLcutChi2NDF) return kFALSE;

  if(cosPoint < fLcutCosPoint[0] || cosPoint > fLcutCosPoint[1]) return kFALSE;

  if(dca < fLcutDCA[0] || dca > fLcutDCA[1]) return kFALSE;

  if(r < fLcutVertexR[0] || r > fLcutVertexR[1]) return kFALSE;

  // all cuts passed

  if(0 == type){
    pdgV0 = 3122;
    if(sign){
      pdgP = 2212;
      pdgN = -211;
    }
    else{
      pdgP = -211;
      pdgN = 2212;
    }
  }
  else{
    pdgV0 = -3122;
    if(sign){
      pdgP = 211;
      pdgN = -2212;
    }
    else{
      pdgP = -2212;
      pdgN = 211;
    }
  }

  return kTRUE;
}
//____________________________________________________________________
Bool_t  AliAODv0KineCuts::V0CutsCommon(const AliAODv0 * const v0) const 
{
  /// V0 cuts common to all V0s

  AliAODTrack* dN, *dP; 

  dP = dynamic_cast<AliAODTrack *>(v0->GetSecondaryVtx()->GetDaughter(0));
  dN = dynamic_cast<AliAODTrack *>(v0->GetSecondaryVtx()->GetDaughter(1));
  
  if(!dN || !dP) return kFALSE;

  Int_t qP = dP->Charge();
  Int_t qN = dN->Charge();

  if((qP*qN) != -1) return kFALSE;

  return kTRUE;
}
//____________________________________________________________________
Bool_t AliAODv0KineCuts::CheckSigns(AliAODv0* const v0) const 
{
  /// check wheter the sign was correctly applied to
  /// V0 daughter tracks

  Bool_t correct = kFALSE;

  AliAODTrack* d[2] = {
    dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(0)),
    dynamic_cast<AliAODTrack*>(v0->GetSecondaryVtx()->GetDaughter(1))
  ,};

  Int_t sign[2];
  sign[0] = d[0]->Charge() > 0. ? 1 : -1;
  sign[1] = d[1]->Charge() > 0. ? 1 : -1;
  
  if(-1 == sign[0] && 1 == sign[1]){
    correct = kFALSE;
  }
  else{
    correct = kTRUE;
  }
  
  return correct;
}
//________________________________________________________________
Double_t AliAODv0KineCuts::PsiPair(AliAODv0* const v0) const 
{
  /// Angle between daughter momentum plane and plane

  if(!fEvent) return -1.;

  Float_t magField = fEvent->GetMagneticField();

  Int_t pIndex = -1;
  Int_t nIndex = -1;
  if(CheckSigns(v0)){
    pIndex = 0;
    nIndex = 1;
  }
  else{
    pIndex = 1;
    nIndex = 0;    
  }
 

  AliESDtrack* daughter[2];

  daughter[0] = new AliESDtrack(dynamic_cast<AliAODTrack *>(v0->GetSecondaryVtx()->GetDaughter(pIndex)));
  daughter[1] = new AliESDtrack(dynamic_cast<AliAODTrack *>(v0->GetSecondaryVtx()->GetDaughter(nIndex)));

  Double_t xyz[3];
  v0->GetXYZ(xyz);//Reconstructed coordinates of V0; to be replaced by Markus Rammler's method in case of conversions!
  
  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0};

  //reconstructed cartesian momentum components of negative daughter;
  mn[0] = v0->MomNegX();
  mn[1] = v0->MomNegY();
  mn[2] = v0->MomNegZ();
  //reconstructed cartesian momentum components of positive daughter;
  mp[0] = v0->MomPosX();
  mp[1] = v0->MomPosY();
  mp[2] = v0->MomPosZ();

  Double_t deltat = 1.;
  deltat = TMath::ATan(mp[2]/(TMath::Sqrt(mp[0]*mp[0] + mp[1]*mp[1])+1.e-13)) -  TMath::ATan(mn[2]/(TMath::Sqrt(mn[0]*mn[0] + mn[1]*mn[1])+1.e-13));//difference of angles of the two daughter tracks with z-axis

  Double_t radiussum = TMath::Sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]) + 50;//radius to which tracks shall be propagated

  Double_t momPosProp[3];
  Double_t momNegProp[3];
    
  AliExternalTrackParam pt(*daughter[0]), nt(*daughter[1]);
    
  Double_t psiPair = 4.;

  if(nt.PropagateTo(radiussum,magField) == 0)//propagate tracks to the outside
    psiPair =  -5.;
  if(pt.PropagateTo(radiussum,magField) == 0)
    psiPair = -5.;
  pt.GetPxPyPz(momPosProp);//Get momentum vectors of tracks after propagation
  nt.GetPxPyPz(momNegProp);
  
  Double_t pEle =
    TMath::Sqrt(momNegProp[0]*momNegProp[0]+momNegProp[1]*momNegProp[1]+momNegProp[2]*momNegProp[2]);//absolute momentum value of negative daughter
  Double_t pPos =
    TMath::Sqrt(momPosProp[0]*momPosProp[0]+momPosProp[1]*momPosProp[1]+momPosProp[2]*momPosProp[2]);//absolute momentum value of positive daughter
    
  Double_t scalarproduct =
    momPosProp[0]*momNegProp[0]+momPosProp[1]*momNegProp[1]+momPosProp[2]*momNegProp[2];//scalar product of propagated positive and negative daughters' momenta
    
  Double_t chipair = TMath::ACos(scalarproduct/(pEle*pPos));//Angle between propagated daughter tracks

  psiPair =  TMath::Abs(TMath::ASin(deltat/chipair));  

  delete daughter[0]; delete daughter[1];
  return psiPair; 
}
//___________________________________________________________________
Bool_t  AliAODv0KineCuts::GetConvPosXY(AliAODTrack * const ptrack, AliAODTrack * const ntrack, Double_t convpos[2]) const
{
  /// recalculate the gamma conversion XY postition

  const Double_t b = fEvent->GetMagneticField();

  AliESDtrack posESD(ptrack), negESD(ntrack);

  Double_t helixcenterpos[2];
  GetHelixCenter(ptrack,b,ptrack->Charge(),helixcenterpos);

  Double_t helixcenterneg[2];
  GetHelixCenter(ntrack,b,ntrack->Charge(),helixcenterneg);

  Double_t  poshelix[6];
  posESD.GetHelixParameters(poshelix,b);
  Double_t posradius = TMath::Abs(1./poshelix[4]);

  Double_t  neghelix[6];
  negESD.GetHelixParameters(neghelix,b);
  Double_t negradius = TMath::Abs(1./neghelix[4]);

  Double_t xpos = helixcenterpos[0];
  Double_t ypos = helixcenterpos[1];
  Double_t xneg = helixcenterneg[0];
  Double_t yneg = helixcenterneg[1];

  convpos[0] = (xpos*negradius + xneg*posradius)/(negradius+posradius);
  convpos[1] = (ypos*negradius+  yneg*posradius)/(negradius+posradius);

  return 1;
}
//___________________________________________________________________
Bool_t  AliAODv0KineCuts::GetHelixCenter(AliAODTrack * const track, Double_t b,Int_t charge, Double_t center[2]) const
{
  /// computes the center of the track helix

  Double_t pi = TMath::Pi();
  
  Double_t  helix[6];
  AliESDtrack esddaughter(track);
  esddaughter.GetHelixParameters(helix,b);
  
  Double_t xpos =  helix[5];
  Double_t ypos =  helix[0];
  Double_t radius = TMath::Abs(1./helix[4]);
  Double_t phi = helix[2];

  if(phi < 0){
    phi = phi + 2*pi;
  }

  phi -= pi/2.;
  Double_t xpoint =  radius * TMath::Cos(phi);
  Double_t ypoint =  radius * TMath::Sin(phi);

  if(b<0){
    if(charge > 0){
      xpoint = - xpoint;
      ypoint = - ypoint;
    }
    /* avoid self assignment
    if(charge < 0){
      xpoint =  xpoint;
      ypoint =  ypoint;
    }
    */
  }
  if(b>0){
    /* avoid self assignment
    if(charge > 0){
      xpoint =  xpoint;
      ypoint =  ypoint;
    }
    */
    if(charge < 0){
      xpoint = - xpoint;
      ypoint = - ypoint;
    }
  }
  center[0] =  xpos + xpoint;
  center[1] =  ypos + ypoint;

  return 1;
}
//___________________________________________________________________
AliKFParticle *AliAODv0KineCuts::CreateMotherParticle(const AliVTrack* const pdaughter, const AliVTrack* const ndaughter, Int_t pspec, Int_t nspec) const
{
  /// Creates a mother particle

  AliKFParticle pkfdaughter(*pdaughter, pspec);
  AliKFParticle nkfdaughter(*ndaughter, nspec);
  
  
  // Create the mother particle 
  AliKFParticle *m = new AliKFParticle(pkfdaughter, nkfdaughter);
  m->SetField(fEvent->GetMagneticField());
  if(TMath::Abs(kElectron) == pspec && TMath::Abs(kElectron) == nspec) m->SetMassConstraint(0, 0.001);
  else if(TMath::Abs(kPiPlus) == pspec && TMath::Abs(kPiPlus) == nspec) m->SetMassConstraint(TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass(), 0.);
  else if(TMath::Abs(kProton) == pspec && TMath::Abs(kPiPlus) == nspec) m->SetMassConstraint(TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass(), 0.);
  else if(TMath::Abs(kPiPlus) == pspec && TMath::Abs(kProton) == nspec) m->SetMassConstraint(TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass(), 0.);
  else{
    AliErrorClass("Wrong daughter ID - mass constraint can not be set");
  }

  AliKFVertex improvedVertex = *fPrimaryVertex;
  improvedVertex += *m;
  m->SetProductionVertex(improvedVertex);
  
  // update 15/06/2010
  // mother particle will not be added to primary vertex but only to its copy 
  // as this confilcts with calling
  // m->SetPrimaryVertex() function and
  // subsequently removing the mother particle afterwards
  // Source: Sergey Gorbunov

  return m;
}
//____________________________________________________________________
void  AliAODv0KineCuts::SetEvent(AliAODEvent* const event){
  /// direct setter of AOD event

  fEvent = event;
  if(!fEvent){
    AliErrorClass("Invalid input event pointer");
    return;
  }

  // Set Mother vertex List
  fNumberKinkMothers = 0;
  if(fEvent->GetNumberOfVertices() > fKinkMotherList.GetSize()) fKinkMotherList.Set(fEvent->GetNumberOfVertices());
  for(Int_t ivertex=0; ivertex < fEvent->GetNumberOfVertices(); ivertex++) {
    AliAODVertex *aodvertex = fEvent->GetVertex(ivertex);
    if(!aodvertex) continue;
    if(aodvertex->GetType()==AliAODVertex::kKink) {
      AliAODTrack *mother = (AliAODTrack *) aodvertex->GetParent();
      if(!mother) continue;
      Int_t idmother = mother->GetID();
      fKinkMotherList[fNumberKinkMothers++] = idmother;
    }
  }

  if (fUseExternalVertex) return;
  else{
	if(fPrimaryVertex && fDeleteVertex){
		delete 	fPrimaryVertex;
		fPrimaryVertex=0x0;
	}
	fPrimaryVertex = new AliKFVertex(*(fEvent->GetPrimaryVertex()));
	fDeleteVertex=kTRUE;
	}



}
//____________________________________________________________________
void  AliAODv0KineCuts::SetEvent(AliVEvent* const event){
  /// direct setter of AOD event

  fEvent = dynamic_cast<AliAODEvent*>(event);
  if(!fEvent){
    AliErrorClass("Invalid input event pointer");
    return;
  }
  
  // Set Mother vertex List
  fNumberKinkMothers = 0;
  if(fEvent->GetNumberOfVertices() > fKinkMotherList.GetSize()) fKinkMotherList.Set(fEvent->GetNumberOfVertices());
  for(Int_t ivertex=0; ivertex < fEvent->GetNumberOfVertices(); ivertex++) {
    AliAODVertex *aodvertex = fEvent->GetVertex(ivertex);
    if(!aodvertex) continue;
    if(aodvertex->GetType()==AliAODVertex::kKink) {
      AliAODTrack *mother = (AliAODTrack *) aodvertex->GetParent();
      if(!mother) continue;
      Int_t idmother = mother->GetID();
      fKinkMotherList[fNumberKinkMothers++] = idmother;
    }
  }

  if (fUseExternalVertex) return;
  else{
    if(fPrimaryVertex && fDeleteVertex){
      delete 	fPrimaryVertex;
      fPrimaryVertex=0x0;
      }
    fPrimaryVertex = new AliKFVertex(*(fEvent->GetPrimaryVertex()));
    fDeleteVertex=kTRUE;
  }
}


//________________________________________________________________
void	 AliAODv0KineCuts::UseExternalVertex(Bool_t use_external){
 /// Reenable primary Vertex from AOD event

	if (use_external) fUseExternalVertex =kTRUE;
	else fUseExternalVertex =kFALSE;
}




//________________________________________________________________
void AliAODv0KineCuts::SetPrimaryVertex(AliKFVertex* const v){
  /// set the primary vertex of the event

  	if(fPrimaryVertex && fDeleteVertex){   
		delete 	fPrimaryVertex;
		fPrimaryVertex =0x0;
		fDeleteVertex = kFALSE;
		}  
  fUseExternalVertex=kTRUE; 
  fPrimaryVertex = v; // set primary Vertex
  if(!fPrimaryVertex){
    AliErrorClass("Failed to initialize the primary vertex");
    return;
  }
}
//___________________________________________________________________
void AliAODv0KineCuts::SetMode(Int_t mode, Int_t type){
  /// this function allows the user to select (prior running the 'ProcessV0' function)
  /// to select different approaches to V0 selection - the 'mode'
  /// - and -
  /// different systems (pp, PbPb) - 'type'
  ///
  /// To see the cut values for different modes please refer to the
  /// function SetCuts()
  ///
  /// Important notice: based on the parameters particular sets of cuts will
  /// be activated for teh V0 selection. If some additional changes to single
  /// cuts are needed please us the SetXXXcut function (see the header file)

  switch(mode){
  case kPurity:
    fMode = kPurity;  // used to obtain highest purity possible - the efficiency may be low
    break;
  case kEffGamma:
    fMode = kEffGamma; // used to obtain highes efficiency possible - the purity may be worse
    break;
  default:
    AliError("V0 selection mode not recognozed, setting 'kPurity'");
    fMode = kPurity;
  }

  switch(type){
  case kPP:
    fType = kPP;  // cuts optimized for low multiplicity 
    break;
  case kPbPb:
    fType = kPbPb;  // cuts optimized for high multiplicity
    break;
  }
  
  // setup the cut values for selected mode & type
  SetCuts();

}
//___________________________________________________________________
void AliAODv0KineCuts::SetMode(Int_t mode, const char* type){
  /// overloaded function - please see above

  Int_t t = -1;

  if(!strcmp("pp", type)) t = kPP;
  else if(!(strcmp("PbPb", type))) t = kPbPb;
  else{
    AliError("data type not recognized, setting 'pp'");
    t = kPP;    
  }

  SetMode(mode, t);

}
//___________________________________________________________________
void AliAODv0KineCuts::SetCuts(){
  /// this funciton sets the default cut values based on the selected
  /// fMode and fType.
  /// please note that only the cuts that have different values than the default
  /// cuts are updated here

  // last update: 14/02/2011
  // as a very preliminary  - the only change to default cuts is to apply
  // less restricting gamma conversion selection in PreselectV0() function
  

  
}

//___________________________________________________________________
Bool_t AliAODv0KineCuts::IsKinkMother(const AliAODTrack * const track) const {
  /// Check if track is a kink mother

  for(int ivtx = 0; ivtx < fNumberKinkMothers; ivtx++){
    if(track->GetID() == fKinkMotherList[ivtx]) return kTRUE;
  }
  return kFALSE;
}

//___________________________________________________________________
Bool_t AliAODv0KineCuts::IsKinkDaughter(const AliAODTrack * const track) const {
  /// Check if track is a kink daughter

  AliAODVertex *vtx = track->GetProdVertex();
  if(!vtx) return kFALSE;
  if(vtx->GetType()==AliAODVertex::kKink) return kTRUE;
  return kFALSE;
}

