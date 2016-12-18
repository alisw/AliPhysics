////////////////////////////////////////////////////////////////////////////////////
//                                                                                //
//  AliFemtoModelWeightGeneratorBasicLednicky.cxx                                 //
//                                                                                //
////////////////////////////////////////////////////////////////////////////////////


#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoModelWeightGeneratorBasicLednicky);
  /// \endcond
#endif

#include "AliFemtoModelWeightGeneratorBasicLednicky.h"
#include "AliFemtoModelHiddenInfo.h"

//_________________________________________
AliFemtoModelWeightGeneratorBasicLednicky::AliFemtoModelWeightGeneratorBasicLednicky() :
  AliFemtoModelWeightGenerator(),
  fIdenticalParticles(false),
  fParamLambda(0),
  fParamAlpha(0),
  fParamRadius(0),
  fParamRef0(0),
  fParamImf0(0),
  fParamd0(0),
  fParamNorm(1.)
{
  //default constructor
}

//_________________________________________
AliFemtoModelWeightGeneratorBasicLednicky::AliFemtoModelWeightGeneratorBasicLednicky(bool aAreIdentical, double aLambda, double aAlpha, double aRadius, double aRef0, double aImf0, double ad0, double aNorm) :
  AliFemtoModelWeightGenerator(),
  fIdenticalParticles(aAreIdentical),
  fParamLambda(aLambda),
  fParamAlpha(aAlpha),
  fParamRadius(aRadius),
  fParamRef0(aRef0),
  fParamImf0(aImf0),
  fParamd0(ad0),
  fParamNorm(aNorm)
{
  //normal constructor
}


//_________________________________________
AliFemtoModelWeightGeneratorBasicLednicky::~AliFemtoModelWeightGeneratorBasicLednicky()
{
  //destructor
}


//_________________________________________
AliFemtoModelWeightGeneratorBasicLednicky::AliFemtoModelWeightGeneratorBasicLednicky(const AliFemtoModelWeightGeneratorBasicLednicky& aGenerator) :
  AliFemtoModelWeightGenerator(aGenerator),
  fIdenticalParticles(aGenerator.fIdenticalParticles),
  fParamLambda(aGenerator.fParamLambda),
  fParamAlpha(aGenerator.fParamAlpha),
  fParamRadius(aGenerator.fParamRadius),
  fParamRef0(aGenerator.fParamRef0),
  fParamImf0(aGenerator.fParamImf0),
  fParamd0(aGenerator.fParamd0),
  fParamNorm(aGenerator.fParamNorm)
{
  //copy constructor
}


//_________________________________________
AliFemtoModelWeightGeneratorBasicLednicky& AliFemtoModelWeightGeneratorBasicLednicky::operator=(const AliFemtoModelWeightGeneratorBasicLednicky& aGenerator)
{
  //assignment operator
  if (this == &aGenerator) return *this;

  AliFemtoModelWeightGenerator::operator=(aGenerator);

  fIdenticalParticles = aGenerator.fIdenticalParticles;
  fParamLambda = aGenerator.fParamLambda;
  fParamAlpha = aGenerator.fParamAlpha;
  fParamRadius = aGenerator.fParamRadius;
  fParamRef0 = aGenerator.fParamRef0;
  fParamImf0 = aGenerator.fParamImf0;
  fParamd0 = aGenerator.fParamd0;
  fParamNorm = aGenerator.fParamNorm;

  return *this;
}

//_________________________________________
double AliFemtoModelWeightGeneratorBasicLednicky::GetLednickyF1(double z)
{
  double result = (1./z)*Faddeeva::Dawson(z);
  return result;
}


//_________________________________________
double AliFemtoModelWeightGeneratorBasicLednicky::GetLednickyF2(double z)
{
  double result = (1./z)*(1.-exp(-z*z));
  return result;
}

//_________________________________________
double AliFemtoModelWeightGeneratorBasicLednicky::GetLednickyWeight(double aKStar)
{
  std::complex<double> tf0 (fParamRef0,fParamImf0);

  double tAlpha;
  if(!fIdenticalParticles) tAlpha=0.;
  else tAlpha = fParamAlpha;

  double tKStar = aKStar/hbarc;

  double tZ = 2.0*tKStar*fParamRadius;  //tZ = 2k*R, to be fed to GetLednickyF1 and GetLednickyF2

  double tCQS = tAlpha*exp(-tZ*tZ);  //quantum statistics term

  std::complex<double> tScattAmp = pow( (1./tf0) + 0.5*fParamd0*tKStar*tKStar - ImI*tKStar,-1);

  double tCFSI = (1.+tAlpha) * (0.5*std::norm(tScattAmp)/(fParamRadius*fParamRadius)*(1.-1./(2*std::sqrt(TMath::Pi()))*(fParamd0/fParamRadius)) + 2.*std::real(tScattAmp)/(fParamRadius*std::sqrt(TMath::Pi()))*GetLednickyF1(tZ) - (std::imag(tScattAmp)/fParamRadius)*GetLednickyF2(tZ));

  double tC = 1. + fParamLambda*(tCQS+tCFSI);
  tC *= fParamNorm;

  return tC;
}

//_________________________________________
double AliFemtoModelWeightGeneratorBasicLednicky::GetKStarTrue(AliFemtoPair *aPair)
{
  double tVerySmall = std::numeric_limits < double >::min();

  AliFemtoParticle *tPart1 = (AliFemtoParticle*)aPair->Track1();
  AliFemtoParticle *tPart2 = (AliFemtoParticle*)aPair->Track2();

  if(tPart1 != NULL && tPart2 != NULL)
  {
    AliFemtoModelHiddenInfo *tPart1HiddenInfo = (AliFemtoModelHiddenInfo*)tPart1->GetHiddenInfo();
    AliFemtoModelHiddenInfo *tPart2HiddenInfo = (AliFemtoModelHiddenInfo*)tPart2->GetHiddenInfo();

    if(tPart1HiddenInfo != NULL && tPart2HiddenInfo != NULL)
    {
      double px1 = tPart1HiddenInfo->GetTrueMomentum()->x();
      double py1 = tPart1HiddenInfo->GetTrueMomentum()->y();
      double pz1 = tPart1HiddenInfo->GetTrueMomentum()->z();
      double mass1 = tPart1HiddenInfo->GetMass();
      double E1 = std::sqrt(mass1*mass1 + px1*px1 + py1*py1 + pz1*pz1);
      if( (E1-mass1) < tVerySmall) return -999;  //p1 has zero momentum!

      double px2 = tPart2HiddenInfo->GetTrueMomentum()->x();
      double py2 = tPart2HiddenInfo->GetTrueMomentum()->y();
      double pz2 = tPart2HiddenInfo->GetTrueMomentum()->z();
      double mass2 = tPart2HiddenInfo->GetMass();
      double E2 = std::sqrt(mass2*mass2 + px2*px2 + py2*py2 + pz2*pz2);
      if( (E2-mass2) < tVerySmall) return -999;  //p2 has zero momentum!
      //------------------------------------------------------------

      double tMinvSq = (E1+E2)*(E1+E2) - (px1+px2)*(px1+px2) - (py1+py2)*(py1+py2) - (pz1+pz2)*(pz1+pz2);

      double tQinvSq = ((mass1*mass1 - mass2*mass2)*(mass1*mass1 - mass2*mass2))/tMinvSq + tMinvSq - 2.0*(mass1*mass1 + mass2*mass2);

      double tKStarTrue = 0.5*sqrt(tQinvSq);
/*
      //Go to LCMS
      double tPx = px1+px2;
      double tPy = py1+py2;
      double tPz = pz1+pz2;
      double tPE = E1+E2;

      double tPtrans = std::sqrt(tPx*tPx + tPy*tPy);
      double tMtrans = std::sqrt(tPE*tPE - tPz*tPz); 

      double beta = tPz/tPE;
      double gamma = tPE/tMtrans;

      double pz1L = gamma*(pz1 - beta*E1);
      double pE1L = gamma*(E1 - beta*pz1);

      double tKStarLong = pz1L;

      //rotation px->tPt
      double px1R = (px1*tPx + py1*tPy)/tPtrans;
      double py1R = (-px1*tPy + py1*tPx)/tPtrans;

      double tKStarSide = py1R;

      //go from LCMS to CMS
      beta = tPtrans/tMtrans;
      gamma = tMtrans/std::sqrt(tMinvSq);

      double px1C = gamma * (px1R - beta*pE1L);
      double tKStarOut = px1C;
*/
      return tKStarTrue;
    }
  }
  return -999;
}


//_________________________________________
Double_t AliFemtoModelWeightGeneratorBasicLednicky::GenerateWeight(AliFemtoPair *aPair)
{
  double tKStarTrue = GetKStarTrue(aPair);
  if(tKStarTrue == 0 || tKStarTrue == -999) return 0;  //bad particle, so fill with weight 0

  return GetLednickyWeight(tKStarTrue);

}

