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

#include "AliFemtoModelWeightGenerator.h"
#include "AliFemtoModelHiddenInfo.h"

const Int_t AliFemtoModelWeightGenerator::kPionPlusPionPlus = 1;
const Int_t AliFemtoModelWeightGenerator::kPionPlusPionMinus = 2;
const Int_t AliFemtoModelWeightGenerator::kKaonPlusKaonPlus = 3;
const Int_t AliFemtoModelWeightGenerator::kKaonPlusKaonMinus = 4;
const Int_t AliFemtoModelWeightGenerator::kProtonProton = 5;
const Int_t AliFemtoModelWeightGenerator::kProtonAntiproton = 6;
const Int_t AliFemtoModelWeightGenerator::kPionPlusKaonPlus = 7;
const Int_t AliFemtoModelWeightGenerator::kPionPlusKaonMinus = 8;
const Int_t AliFemtoModelWeightGenerator::kPionPlusProton = 9;
const Int_t AliFemtoModelWeightGenerator::kPionPlusAntiproton = 10;
const Int_t AliFemtoModelWeightGenerator::kKaonPlusProton = 11;
const Int_t AliFemtoModelWeightGenerator::kKaonPlusAntiproton = 12;

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
Int_t    AliFemtoModelWeightGenerator::GetPairType()
{
  return fPairType;
}

//_____________________________________________
void     AliFemtoModelWeightGenerator::SetPairTypeFromPair(AliFemtoPair *aPair)
{
  AliFemtoModelHiddenInfo *inf1 = ( AliFemtoModelHiddenInfo *) aPair->track1()->HiddenInfo();
  AliFemtoModelHiddenInfo *inf2 = ( AliFemtoModelHiddenInfo *) aPair->track2()->HiddenInfo();

  const Int_t tPid1 = inf1->GetPDGPid();
  const Int_t tPid2 = inf2->GetPDGPid();

  if      (((tPid1 ==   211) && (tPid2 ==   211)) ||
           ((tPid1 ==  -211) && (tPid2 ==  -211)))
    fPairType = kPionPlusPionPlus;
  else if (((tPid1 ==  -211) && (tPid2 ==   211)) ||
           ((tPid1 ==   211) && (tPid2 ==  -211)))
    fPairType = kPionPlusPionMinus;
  else if (((tPid1 ==   321) && (tPid2 ==   321)) ||
           ((tPid1 ==  -321) && (tPid2 ==  -321)))
    fPairType = kKaonPlusKaonPlus;
  else if (((tPid1 ==  -321) && (tPid2 ==   321)) ||
           ((tPid1 ==   321) && (tPid2 ==  -321)))
    fPairType = kKaonPlusKaonMinus;
  else if (((tPid1 ==  2212) && (tPid2 ==  2212)) ||
           ((tPid1 == -2212) && (tPid2 == -2212)))
    fPairType = kProtonProton;
  else if (((tPid1 == -2212) && (tPid2 ==  2212)) ||
           ((tPid1 ==  2212) && (tPid2 == -2212)))
    fPairType = kProtonAntiproton;
  else if (((tPid1 ==   211) && (tPid2 ==   321)) ||
           ((tPid1 ==  -211) && (tPid2 ==  -321)))
    fPairType = kPionPlusKaonPlus;
  else if (((tPid1 ==  -211) && (tPid2 ==   321)) ||
           ((tPid1 ==   211) && (tPid2 ==  -321)))
    fPairType = kPionPlusKaonMinus;
  else if (((tPid1 ==   211) && (tPid2 ==  2212)) ||
           ((tPid1 ==  -211) && (tPid2 == -2212)))
    fPairType = kPionPlusProton;
  else if (((tPid1 ==  -211) && (tPid2 ==  2212)) ||
           ((tPid1 ==   211) && (tPid2 == -2212)))
    fPairType = kPionPlusAntiproton;
  else if (((tPid1 ==   321) && (tPid2 ==  2212)) ||
           ((tPid1 ==  -321) && (tPid2 == -2212)))
    fPairType = kKaonPlusProton;
  else if (((tPid1 ==  -321) && (tPid2 ==  2212)) ||
           ((tPid1 ==   321) && (tPid2 == -2212)))
    fPairType = kKaonPlusAntiproton;
}

//_____________________________________________
AliFemtoModelWeightGenerator* AliFemtoModelWeightGenerator::Clone() const
{
  return 0;
}
