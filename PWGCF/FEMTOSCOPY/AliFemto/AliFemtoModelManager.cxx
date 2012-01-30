////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelManager - main helper class for femtoscopy calculations     ///
/// Manages weight generation, freeze-out coordinates generation             ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifdef __ROOT__
  ClassImp(AliFemtoModelManager, 1)
#endif

#include "AliFemtoModelManager.h"
#include "AliFemtoModelHiddenInfo.h"

//_____________________________________________
AliFemtoModelManager::AliFemtoModelManager():
  fFreezeOutGenerator(0),
  fWeightGenerator(0),
  fCreateCopyHiddenInfo(kFALSE)
{
}
//_____________________________________________
AliFemtoModelManager::AliFemtoModelManager(const AliFemtoModelManager& aManager):
  fFreezeOutGenerator(0),
  fWeightGenerator(0),
  fCreateCopyHiddenInfo(aManager.fCreateCopyHiddenInfo)
{
  if (aManager.fFreezeOutGenerator) {
    fFreezeOutGenerator = aManager.fFreezeOutGenerator->Clone();
  }
  if (aManager.fWeightGenerator) {
    fWeightGenerator = aManager.fWeightGenerator->Clone();
  }
}
//_____________________________________________
AliFemtoModelManager::~AliFemtoModelManager()
{
  if (fFreezeOutGenerator) delete fFreezeOutGenerator;
  if (fWeightGenerator) delete fWeightGenerator;
}
//_____________________________________________
AliFemtoModelManager& AliFemtoModelManager::operator=(const AliFemtoModelManager& aManager)
{
  if (this == &aManager)
    return *this;
  if (aManager.fFreezeOutGenerator) {
    fFreezeOutGenerator = aManager.fFreezeOutGenerator->Clone();
  }
  else fFreezeOutGenerator = 0;
  if (aManager.fWeightGenerator) {
    fWeightGenerator = aManager.fWeightGenerator->Clone();
  }
  else fWeightGenerator = 0;
  fCreateCopyHiddenInfo = aManager.fCreateCopyHiddenInfo;
  
  return *this;
}
//_____________________________________________
void AliFemtoModelManager::AcceptFreezeOutGenerator(AliFemtoModelFreezeOutGenerator *aFreeze)
{
  fFreezeOutGenerator = aFreeze;
}
//_____________________________________________
void AliFemtoModelManager::AcceptWeightGenerator(AliFemtoModelWeightGenerator *aWeight)
{
  fWeightGenerator = aWeight;
}
//_____________________________________________
Double_t AliFemtoModelManager::GetWeight(AliFemtoPair *aPair)
{
  if (!fWeightGenerator) {
    cout << "No weight generator set! Cannot calculate weight" << endl;
    exit(0);
  }
  // Return femtoscopic weight for a given pair
  if (fCreateCopyHiddenInfo) {
    // Try to guess particle masses and pid from the weight generator
    Double_t tMass1=0.0001, tMass2=0.0001;
    Int_t tPid1=0, tPid2=0;
    if (fWeightGenerator->GetPairType() == AliFemtoModelWeightGenerator::PionPlusPionPlus()) {
      tMass1 = 0.13957;
      tMass2 = 0.13957;
      tPid1 = 211;
      tPid2 = 211;
    }
    if (fWeightGenerator->GetPairType() == AliFemtoModelWeightGenerator::PionPlusPionMinus()) {
      tMass1 = 0.13957;
      tMass2 = 0.13957;
      tPid1 = 211;
      tPid2 = -211;
    }
    if (fWeightGenerator->GetPairType() == AliFemtoModelWeightGenerator::KaonPlusKaonPlus()) {
      tMass1 = 0.493677;
      tMass2 = 0.493677;
      tPid1 = 321;
      tPid2 = 321;
    }
    if (fWeightGenerator->GetPairType() == AliFemtoModelWeightGenerator::KaonPlusKaonMinus()) {
      tMass1 = 0.493677;
      tMass2 = 0.493677;
      tPid1 = 321;
      tPid2 = -321;
    }
    if (fWeightGenerator->GetPairType() == AliFemtoModelWeightGenerator::ProtonProton()) {
      tMass1 = 0.938272;
      tMass2 = 0.938272;
      tPid1 = 2212;
      tPid2 = 2212;
    }
    if (fWeightGenerator->GetPairType() == AliFemtoModelWeightGenerator::ProtonAntiproton()) {
      tMass1 = 0.938272;
      tMass2 = 0.938272;
      tPid1 = 2212;
      tPid2 = -2212;
    }
    if (fWeightGenerator->GetPairType() == AliFemtoModelWeightGenerator::PionPlusKaonPlus()) {
      tMass1 = 0.13957;
      tMass2 = 0.493677;
      tPid1 = 211;
      tPid2 = 321;
    }
    if (fWeightGenerator->GetPairType() == AliFemtoModelWeightGenerator::PionPlusKaonMinus()) {
      tMass1 = 0.13957;
      tMass2 = 0.493677;
      tPid1 = 211;
      tPid2 = -321;
    }
    if (fWeightGenerator->GetPairType() == AliFemtoModelWeightGenerator::PionPlusProton()) {
      tMass1 = 0.13957;
      tMass2 = 0.938272;
      tPid1 = 211;
      tPid2 = 2212;
    }
    if (fWeightGenerator->GetPairType() == AliFemtoModelWeightGenerator::PionPlusAntiproton()) {
      tMass1 = 0.13957;
      tMass2 = 0.938272;
      tPid1 = 211;
      tPid2 = -2212;
    }
    if (fWeightGenerator->GetPairType() == AliFemtoModelWeightGenerator::KaonPlusProton()) {
      tMass1 = 0.493677;
      tMass2 = 0.938272;
      tPid1 = 321;
      tPid2 = 2212;
    }
    if (fWeightGenerator->GetPairType() == AliFemtoModelWeightGenerator::KaonPlusAntiproton()) {
      tMass1 = 0.493677;
      tMass2 = 0.938272;
      tPid1 = 321;
      tPid2 = -2212;
    }

    if (!(aPair->Track1()->HiddenInfo())) {
      AliFemtoModelHiddenInfo *inf1 = new AliFemtoModelHiddenInfo();
      inf1->SetTrueMomentum(aPair->Track1()->Track()->P());
      inf1->SetMass(tMass1);
      inf1->SetPDGPid(tPid1);
      aPair->Track1()->SetHiddenInfo(inf1);
      delete inf1;
    }
    if (!(aPair->Track2()->HiddenInfo())) {
      AliFemtoModelHiddenInfo *inf2 = new AliFemtoModelHiddenInfo();
      inf2->SetTrueMomentum(aPair->Track2()->Track()->P());
      inf2->SetMass(tMass2);
      inf2->SetPDGPid(tPid2);
      aPair->Track2()->SetHiddenInfo(inf2);
      delete inf2;
    }
  }

  if (fFreezeOutGenerator) {
    fFreezeOutGenerator->GenerateFreezeOut(aPair);
  }
  return fWeightGenerator->GenerateWeight(aPair);
}
//_____________________________________________
void AliFemtoModelManager::CreateCopyHiddenInfo(Bool_t aCopy)
{
  fCreateCopyHiddenInfo = aCopy;
}
//_____________________________________________
AliFemtoModelFreezeOutGenerator* AliFemtoModelManager::GetFreezeOutGenerator()
{
  return fFreezeOutGenerator;
}
//_____________________________________________
AliFemtoModelWeightGenerator*    AliFemtoModelManager::GetWeightGenerator()
{
  return fWeightGenerator;
}
