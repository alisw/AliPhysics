////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelWeightGeneratorBasic -  basic femtoscopic weight generator  ///
/// only return a simple                                                          ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifdef __ROOT__
  ClassImp(AliFemtoModelWeightGeneratorBasic, 1)
#endif

#include "AliFemtoModelWeightGeneratorBasic.h"
#include "AliFemtoModelHiddenInfo.h"

//________________________
AliFemtoModelWeightGeneratorBasic::AliFemtoModelWeightGeneratorBasic():
  AliFemtoModelWeightGenerator()
{
  /* no-op */
}
//________________________
AliFemtoModelWeightGeneratorBasic::AliFemtoModelWeightGeneratorBasic(const AliFemtoModelWeightGeneratorBasic &aModel) :
  AliFemtoModelWeightGenerator(aModel)
{
  /* no-op */
}

//________________________
AliFemtoModelWeightGeneratorBasic::~AliFemtoModelWeightGeneratorBasic()
{
  /* no-op */
}

AliFemtoModelWeightGeneratorBasic& AliFemtoModelWeightGeneratorBasic::operator=(const AliFemtoModelWeightGeneratorBasic &aModel)
{
  if (this != &aModel) {
    AliFemtoModelWeightGenerator::operator=(aModel);
  }
  
  return *this;
}

//________________________
Double_t AliFemtoModelWeightGeneratorBasic::GenerateWeight(AliFemtoPair *aPair)
{
  // Generate a simple femtoscopic weight coming only from 
  // quantum statistics - symmetrization or anti-symmetrization
  // of the pair wave function

  // Get hidden information pointers
  //AliFemtoModelHiddenInfo *inf1 = (AliFemtoModelHiddenInfo *) aPair->Track1()->HiddenInfo();
  //AliFemtoModelHiddenInfo *inf2 = (AliFemtoModelHiddenInfo *) aPair->Track2()->HiddenInfo();
  AliFemtoTrack *inf1 = (AliFemtoTrack *) aPair->Track1()->Track();
  AliFemtoTrack *inf2 = (AliFemtoTrack *) aPair->Track2()->Track();

  // Calculate pair variables
  Double_t tPx = inf1->GetTrueMomentum()->x()+inf2->GetTrueMomentum()->x();
  Double_t tPy = inf1->GetTrueMomentum()->y()+inf2->GetTrueMomentum()->y();
  Double_t tPz = inf1->GetTrueMomentum()->z()+inf2->GetTrueMomentum()->z();
  //  double tE  = inf1->GetTrueMomentum()->e +inf2->GetTrueMomentum()->.e;
  Double_t tM1 = inf1->GetMass();
  Double_t tM2 = inf2->GetMass();
  Double_t tE1 = sqrt(tM1*tM1 + inf1->GetTrueMomentum()->Mag2());
  Double_t tE2 = sqrt(tM2*tM2 + inf2->GetTrueMomentum()->Mag2());
  Double_t tE  = tE1 + tE2;
  Double_t tPt = tPx*tPx + tPy*tPy;
  Double_t tMt = tE*tE - tPz*tPz;//mCVK;
  Double_t tM  = sqrt(tMt - tPt);
  tMt = sqrt(tMt);
  tPt = sqrt(tPt);
  Double_t tBetat = tPt/tMt;

  // Boost to LCMS
  Double_t tBeta = tPz/tE;
  Double_t tGamma = tE/tMt;	    
  fKStarLong = tGamma * (inf1->GetTrueMomentum()->z() - tBeta * tE1);
  Double_t tE1L = tGamma * (tE1  - tBeta * inf1->GetTrueMomentum()->z());
  
  // Transform positions to LCMS
//   Double_t tP1zl = tGamma * (inf1->GetEmissionPoint()->z() - tBeta * inf1->GetEmissionPoint()->t());
//   Double_t tP1tl = tGamma * (inf1->GetEmissionPoint()->t() - tBeta * inf1->GetEmissionPoint()->z());
  
//   Double_t tP2zl = tGamma * (inf2->GetEmissionPoint()->z() - tBeta * inf2->GetEmissionPoint()->t());
//   Double_t tP2tl = tGamma * (inf2->GetEmissionPoint()->t() - tBeta * inf2->GetEmissionPoint()->z());
  
//   Double_t tP1pzl = tGamma * (inf1->GetTrueMomentum()->z() - tBeta * tE1);
//   Double_t tP1el  = tGamma * (tE1  - tBeta * inf1->GetTrueMomentum()->z());
  
//   Double_t tP2pzl = tGamma * (inf2->GetTrueMomentum()->z() - tBeta * tE2);
//   Double_t tP2el  = tGamma * (tE2  - tBeta * inf2->GetTrueMomentum()->z());
  
  // Rotate in transverse plane
  fKStarOut  = ( inf1->GetTrueMomentum()->x()*tPx + inf1->GetTrueMomentum()->y()*tPy)/tPt;
  fKStarSide = (-inf1->GetTrueMomentum()->x()*tPy + inf1->GetTrueMomentum()->y()*tPx)/tPt;
      
//   Double_t tP1pxl = fKStarOut;
//   Double_t tP1pyl = fKStarSide;
  
//   Double_t tP2pxl = (inf2->GetTrueMomentum()->x()*tPx + inf2->GetTrueMomentum()->y()*tPy)/tPt;
//   Double_t tP2pyl = (inf2->GetTrueMomentum()->y()*tPx - inf2->GetTrueMomentum()->x()*tPy)/tPt;;
  
//   Double_t tKO = tP1pxl - tP2pxl;
//   Double_t tKS = tP1pyl - tP2pyl;
//   Double_t tKL = tP1pzl - tP2pzl;
//   Double_t tDE = tP1el  - tP2el;
  
  // save the rotated coordinates in LCMS variables
//   Double_t tP1xl = ( inf1->GetEmissionPoint()->x()*tPx + inf1->GetEmissionPoint()->y()*tPy)/tPt;
//   Double_t tP1yl = (-inf1->GetEmissionPoint()->x()*tPy + inf1->GetEmissionPoint()->y()*tPx)/tPt;

//   Double_t tP2xl = ( inf2->GetEmissionPoint()->x()*tPx + inf2->GetEmissionPoint()->y()*tPy)/tPt;
//   Double_t tP2yl = (-inf2->GetEmissionPoint()->x()*tPy + inf2->GetEmissionPoint()->y()*tPx)/tPt;
  
  // Boost to pair cms
  fKStarOut = tMt/tM * (fKStarOut - tPt/tMt * tE1L);
  
  tBetat = tPt/tMt;
//   Double_t tGammat = 1.0/sqrt(1.0-tBetat*tBetat);
  
//   Double_t tP1xp = tGammat*(tP1xl - tBetat*tP1tl);
//   Double_t tP1tp = tGammat*(tP1tl - tBetat*tP1xl);
  
//   Double_t tP2xp = tGammat*(tP2xl - tBetat*tP2tl);
//   Double_t tP2tp = tGammat*(tP2tl - tBetat*tP2xl);
  
//   Double_t tRO = (tP1xl - tP2xl)/0.197327;
//   Double_t tRS = (tP1yl - tP2yl)/0.197327;
//   Double_t tRL = (tP1zl - tP2zl)/0.197327;
//   Double_t tDT = (tP1tl - tP2tl)/0.197327;
  
  Double_t tDX = inf1->GetEmissionPoint()->x()-inf2->GetEmissionPoint()->x();
  Double_t tDY = inf1->GetEmissionPoint()->y()-inf2->GetEmissionPoint()->y();
  Double_t tRLong = inf1->GetEmissionPoint()->z()-inf2->GetEmissionPoint()->z();
  Double_t tDTime = inf1->GetEmissionPoint()->t()-inf2->GetEmissionPoint()->t();

  Double_t tROut = (tDX*tPx + tDY*tPy)/tPt;
  Double_t tRSide = (-tDX*tPy + tDY*tPx)/tPt;

  fRStarSide = tRSide;
  Double_t tRSS = fRStarSide/0.197327;

  fRStarLong = tGamma*(tRLong - tBeta* tDTime);
  Double_t tDTimePairLCMS = tGamma*(tDTime - tBeta* tRLong);

  Double_t tRLS = fRStarLong/0.197327;
  tBeta = tPt/tMt;
  tGamma = tMt/tM;

  fRStarOut = tGamma*(tROut - tBeta* tDTimePairLCMS);
  Double_t tROS = fRStarOut/0.197327;
//   Double_t tDTimePairCMS = tGamma*(tDTimePairLCMS - tBeta* tROut);
  fRStar = ::sqrt(fRStarOut*fRStarOut + fRStarSide*fRStarSide +
			   fRStarLong*fRStarLong);
  fKStar = ::sqrt(fKStarOut*fKStarOut + fKStarSide*fKStarSide + fKStarLong*fKStarLong);
//   Double_t tRSt = fRStar/0.197327;

  if (fPairType != fgkPairTypeNone) {
    if ((fPairType == PionPlusPionPlus()) || (fPairType == KaonPlusKaonPlus()))
      return 1.0 + cos (2*(fKStarOut * tROS + fKStarSide * tRSS + fKStarLong * tRLS));
    else if (fPairType == ProtonProton())
      return 1.0 - 0.5 * cos (2*(fKStarOut * tROS + fKStarSide * tRSS + fKStarLong * tRLS));
    else 
      return 1.0;
  }
  else {
    Int_t tPairType = GetPairTypeFromPair(aPair);
    if ((tPairType == PionPlusPionPlus()) || (tPairType == KaonPlusKaonPlus()))
      return 1.0 + cos (2*(fKStarOut * tROS + fKStarSide * tRSS + fKStarLong * tRLS));
    else if (tPairType == ProtonProton())
      return 1.0 - 0.5 * cos (2*(fKStarOut * tROS + fKStarSide * tRSS + fKStarLong * tRLS));
    else 
      return 1.0;
    
  }
}

//________________________
void     AliFemtoModelWeightGeneratorBasic::SetPairType(Int_t aPairType)
{
  AliFemtoModelWeightGenerator::SetPairType(aPairType);
}
//________________________
void     AliFemtoModelWeightGeneratorBasic::SetPairTypeFromPair(AliFemtoPair *aPair)
{
  AliFemtoModelWeightGenerator::SetPairTypeFromPair(aPair);
}
//________________________
Int_t    AliFemtoModelWeightGeneratorBasic::GetPairType() const
{
  return AliFemtoModelWeightGenerator::GetPairType();
}
//________________________
AliFemtoModelWeightGenerator* AliFemtoModelWeightGeneratorBasic::Clone() const
{
  return GetGenerator();
}
//________________________
AliFemtoModelWeightGenerator* AliFemtoModelWeightGeneratorBasic::GetGenerator() const
{
  AliFemtoModelWeightGeneratorBasic *tGen = new AliFemtoModelWeightGeneratorBasic(*this);
  return tGen;
}
