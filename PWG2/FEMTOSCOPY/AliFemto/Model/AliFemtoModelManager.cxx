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
  // Return femtoscopic weight for a fiven pair
  if (fCreateCopyHiddenInfo) {
    if (!(aPair->track1()->HiddenInfo())) {
      AliFemtoModelHiddenInfo *inf1 = new AliFemtoModelHiddenInfo();
      inf1->SetTrueMomentum(aPair->track1()->Track()->P());
      inf1->SetMass(0.13957);
      aPair->track1()->SetHiddenInfo(inf1);
      delete inf1;
    }
    if (!(aPair->track2()->HiddenInfo())) {
      AliFemtoModelHiddenInfo *inf2 = new AliFemtoModelHiddenInfo();
      inf2->SetTrueMomentum(aPair->track2()->Track()->P());
      inf2->SetMass(0.13957);
      aPair->track2()->SetHiddenInfo(inf2);
      delete inf2;
    }
  }

  if (fFreezeOutGenerator) {
    fFreezeOutGenerator->GenerateFreezeOut(aPair);
  }
  return fWeightGenerator->GenerateWeight(aPair);
}
//_____________________________________________
void AliFemtoModelManager::CreateCopyHiddenInfo()
{
  fCreateCopyHiddenInfo = kTRUE;
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
