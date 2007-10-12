////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelWeightGenerator - abstract base class for femtoscopic       ///
/// weight generator                                                         ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifdef __ROOT__
  ClassImp(AliFemtoModelGausLCMSFreezeOutGenerator, 1)
#endif

#include "AliFemtoPair.h"

#include "AliFemtoModelWeightGenerator.h"
#include "AliFemtoModelHiddenInfo.h"

const Int_t AliFemtoModelWeightGenerator::fgkPionPlusPionPlus = 1;
const Int_t AliFemtoModelWeightGenerator::fgkPionPlusPionMinus = 2;
const Int_t AliFemtoModelWeightGenerator::fgkKaonPlusKaonPlus = 3;
const Int_t AliFemtoModelWeightGenerator::fgkKaonPlusKaonMinus = 4;
const Int_t AliFemtoModelWeightGenerator::fgkProtonProton = 5;
const Int_t AliFemtoModelWeightGenerator::fgkProtonAntiproton = 6;
const Int_t AliFemtoModelWeightGenerator::fgkPionPlusKaonPlus = 7;
const Int_t AliFemtoModelWeightGenerator::fgkPionPlusKaonMinus = 8;
const Int_t AliFemtoModelWeightGenerator::fgkPionPlusProton = 9;
const Int_t AliFemtoModelWeightGenerator::fgkPionPlusAntiproton = 10;
const Int_t AliFemtoModelWeightGenerator::fgkKaonPlusProton = 11;
const Int_t AliFemtoModelWeightGenerator::fgkKaonPlusAntiproton = 12;

//_____________________________________________
AliFemtoModelWeightGenerator::AliFemtoModelWeightGenerator() :
  fPairType(0), 
  fKStarOut(0), fKStarSide(0), fKStarLong(0), fKStar(0), 
  fRStarOut(0), fRStarSide(0), fRStarLong(0), fRStar(0)
{}
//_____________________________________________
AliFemtoModelWeightGenerator::AliFemtoModelWeightGenerator(const AliFemtoModelWeightGenerator &aModel) :
  fPairType(0), 
  fKStarOut(0), fKStarSide(0), fKStarLong(0), fKStar(0), 
  fRStarOut(0), fRStarSide(0), fRStarLong(0), fRStar(0)
{
  fPairType = aModel.fPairType;
}
//_____________________________________________
AliFemtoModelWeightGenerator::~AliFemtoModelWeightGenerator(){/* no-op */}
//_____________________________________________
void     AliFemtoModelWeightGenerator::SetPairType(Int_t aPairType)
{
  fPairType = aPairType;
}

//_____________________________________________
Int_t    AliFemtoModelWeightGenerator::GetPairType() const
{
  return fPairType;
}

//_____________________________________________
void     AliFemtoModelWeightGenerator::SetPairTypeFromPair(AliFemtoPair *aPair)
{
  // Get the type of pair from PID of particles in the pair
  AliFemtoModelHiddenInfo *inf1 = ( AliFemtoModelHiddenInfo *) aPair->Track1()->HiddenInfo();
  AliFemtoModelHiddenInfo *inf2 = ( AliFemtoModelHiddenInfo *) aPair->Track2()->HiddenInfo();

  const Int_t ktPid1 = inf1->GetPDGPid();
  const Int_t ktPid2 = inf2->GetPDGPid();

  if      (((ktPid1 ==   211) && (ktPid2 ==   211)) ||
           ((ktPid1 ==  -211) && (ktPid2 ==  -211)))
    fPairType = fgkPionPlusPionPlus;
  else if (((ktPid1 ==  -211) && (ktPid2 ==   211)) ||
           ((ktPid1 ==   211) && (ktPid2 ==  -211)))
    fPairType = fgkPionPlusPionMinus;
  else if (((ktPid1 ==   321) && (ktPid2 ==   321)) ||
           ((ktPid1 ==  -321) && (ktPid2 ==  -321)))
    fPairType = fgkKaonPlusKaonPlus;
  else if (((ktPid1 ==  -321) && (ktPid2 ==   321)) ||
           ((ktPid1 ==   321) && (ktPid2 ==  -321)))
    fPairType = fgkKaonPlusKaonMinus;
  else if (((ktPid1 ==  2212) && (ktPid2 ==  2212)) ||
           ((ktPid1 == -2212) && (ktPid2 == -2212)))
    fPairType = fgkProtonProton;
  else if (((ktPid1 == -2212) && (ktPid2 ==  2212)) ||
           ((ktPid1 ==  2212) && (ktPid2 == -2212)))
    fPairType = fgkProtonAntiproton;
  else if (((ktPid1 ==   211) && (ktPid2 ==   321)) ||
           ((ktPid1 ==  -211) && (ktPid2 ==  -321)))
    fPairType = fgkPionPlusKaonPlus;
  else if (((ktPid1 ==  -211) && (ktPid2 ==   321)) ||
           ((ktPid1 ==   211) && (ktPid2 ==  -321)))
    fPairType = fgkPionPlusKaonMinus;
  else if (((ktPid1 ==   211) && (ktPid2 ==  2212)) ||
           ((ktPid1 ==  -211) && (ktPid2 == -2212)))
    fPairType = fgkPionPlusProton;
  else if (((ktPid1 ==  -211) && (ktPid2 ==  2212)) ||
           ((ktPid1 ==   211) && (ktPid2 == -2212)))
    fPairType = fgkPionPlusAntiproton;
  else if (((ktPid1 ==   321) && (ktPid2 ==  2212)) ||
           ((ktPid1 ==  -321) && (ktPid2 == -2212)))
    fPairType = fgkKaonPlusProton;
  else if (((ktPid1 ==  -321) && (ktPid2 ==  2212)) ||
           ((ktPid1 ==   321) && (ktPid2 == -2212)))
    fPairType = fgkKaonPlusAntiproton;
}

//_____________________________________________
AliFemtoModelWeightGenerator* AliFemtoModelWeightGenerator::Clone() const
{
  return 0;
}
