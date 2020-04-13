///
/// \file AliFemtoModelWeightGeneratorBasic.cxx
/// \author Adam Kisiel kisiel@mps.ohio-state.edu
///

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoModelWeightGeneratorBasic, 1);
  /// \endcond
#endif

#include "AliFemtoModelWeightGeneratorBasic.h"
#include "AliFemtoModelHiddenInfo.h"

//________________________
AliFemtoModelWeightGeneratorBasic::AliFemtoModelWeightGeneratorBasic():
  AliFemtoModelWeightGenerator(),
  fPrintEmptyParticleNotification(true)
{
  /* no-op */
}
//________________________
AliFemtoModelWeightGeneratorBasic::AliFemtoModelWeightGeneratorBasic(const AliFemtoModelWeightGeneratorBasic &aModel) :
  AliFemtoModelWeightGenerator(aModel),
  fPrintEmptyParticleNotification(aModel.fPrintEmptyParticleNotification)
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

  fPrintEmptyParticleNotification = aModel.fPrintEmptyParticleNotification;

  return *this;
}

//________________________
Double_t AliFemtoModelWeightGeneratorBasic::GenerateWeight(AliFemtoPair *aPair)
{
  // Generate a simple femtoscopic weight coming only from
  // quantum statistics - symmetrization or anti-symmetrization
  // of the pair wave function
  {
    double cached_weight = aPair->LookupFemtoWeightCache(this);
    if (!std::isnan(cached_weight)) {
      return cached_weight;
    }
  }

  // Get hidden information pointers
  //AliFemtoModelHiddenInfo *inf1 = (AliFemtoModelHiddenInfo *) aPair->Track1()->HiddenInfo();
  //AliFemtoModelHiddenInfo *inf2 = (AliFemtoModelHiddenInfo *) aPair->Track2()->HiddenInfo();
  AliFemtoTrack *inf1 = (AliFemtoTrack *) aPair->Track1()->Track();
  AliFemtoTrack *inf2 = (AliFemtoTrack *) aPair->Track2()->Track();

  AliFemtoModelHiddenInfo &info1 = *(AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo(),
                          &info2 = *(AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo();

  const auto &p1 = *info1.GetTrueMomentum(),
             &p2 = *info2.GetTrueMomentum();

  const auto P = p1 + p2;

  // Calculate pair variables

  const Double_t tPx = P.x(),
                 tPy = P.y(),
                 tPz = P.z();

  Double_t tM1 = info1.GetMass();
  Double_t tM2 = info2.GetMass();

  Double_t tE1 = sqrt(tM1*tM1 + p1.Mag2());
  Double_t tE2 = sqrt(tM2*tM2 + p2.Mag2());

  Double_t tE  = tE1 + tE2;
  Double_t tPt = tPx*tPx + tPy*tPy;
  Double_t tMt = tE*tE - tPz*tPz;//mCVK;
  Double_t tM  = (tMt - tPt > 0.0) ? sqrt(tMt - tPt) : 0.0;

  if (tMt == 0 || tE == 0 || tM == 0 || tPt == 0 ) {
    if (fPrintEmptyParticleNotification) {
      cout << " weight generator zero tPt || tMt || tM || tPt " << tM1 << " " << tM2 << endl;
    }
    return 0.0;
  }

  tMt = sqrt(tMt);
  tPt = sqrt(tPt);

  Double_t pX = p1.x();
  Double_t pY = p1.y();
  Double_t pZ = p1.z();

  // Boost to LCMS
  Double_t tBeta = tPz/tE;
  Double_t tGamma = tE/tMt;
  fKStarLong = tGamma * (pZ - tBeta * tE1);
  Double_t tE1L = tGamma * (tE1  - tBeta * pZ);

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
  fKStarOut  = ( pX*tPx + pY*tPy)/tPt;
  fKStarSide = (-pX*tPy + pY*tPx)/tPt;

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

//   Double_t tBetat = tPt/tMt;
//   Double_t tGammat = 1.0/sqrt(1.0-tBetat*tBetat);

//   Double_t tP1xp = tGammat*(tP1xl - tBetat*tP1tl);
//   Double_t tP1tp = tGammat*(tP1tl - tBetat*tP1xl);

//   Double_t tP2xp = tGammat*(tP2xl - tBetat*tP2tl);
//   Double_t tP2tp = tGammat*(tP2tl - tBetat*tP2xl);

//   Double_t tRO = (tP1xl - tP2xl)/0.197327;
//   Double_t tRS = (tP1yl - tP2yl)/0.197327;
//   Double_t tRL = (tP1zl - tP2zl)/0.197327;
//   Double_t tDT = (tP1tl - tP2tl)/0.197327;

  // separation distance
  const auto D = *info1.GetEmissionPoint() - *info2.GetEmissionPoint();

  Double_t tDX = D.x();
  Double_t tDY = D.y();
  Double_t tRLong = D.z();
  Double_t tDTime = D.t();

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
  fRStar = ::sqrt(fRStarOut*fRStarOut + fRStarSide*fRStarSide + fRStarLong*fRStarLong);
  fKStar = ::sqrt(fKStarOut*fKStarOut + fKStarSide*fKStarSide + fKStarLong*fKStarLong);
//   Double_t tRSt = fRStar/0.197327;

  // if type not set, use pair to determine type
  auto pair_type = (fPairType == fgkPairTypeNone)
                 ? GetPairTypeFromPair(aPair)
                 : fPairType;

  double weight = ((pair_type == PionPlusPionPlus()) || (pair_type == KaonPlusKaonPlus()))
                  ? 1.0 + cos (2*(fKStarOut * tROS + fKStarSide * tRSS + fKStarLong * tRLS))
                : pair_type == ProtonProton()
                  ? 1.0 - 0.5 * cos (2*(fKStarOut * tROS + fKStarSide * tRSS + fKStarLong * tRLS))
                : 1.0;

  aPair->AddWeightToCache(this, weight);
  return weight;
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
