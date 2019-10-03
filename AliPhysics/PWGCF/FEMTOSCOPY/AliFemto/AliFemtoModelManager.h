////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelManager - main helper class for femtoscopy calculations     ///
/// Manages weight generation, freeze-out coordinates generation             ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoModelManager_hh
#define AliFemtoModelManager_hh

#include "AliFemtoEnumeration.h"
#include "AliFemtoModelWeightGenerator.h"
#include "AliFemtoModelFreezeOutGenerator.h"

class AliFemtoModelManager 
{
 public:
  AliFemtoModelManager();
  AliFemtoModelManager(const AliFemtoModelManager& aManager);
  virtual  ~AliFemtoModelManager();

  AliFemtoModelManager& operator=(const AliFemtoModelManager& aManager);

  void AcceptFreezeOutGenerator(AliFemtoModelFreezeOutGenerator *aFreeze);
  void AcceptWeightGenerator(AliFemtoModelWeightGenerator *aWeight);
  void CreateCopyHiddenInfo(Bool_t aCopy=kTRUE);

  AliFemtoModelFreezeOutGenerator* GetFreezeOutGenerator();
  AliFemtoModelWeightGenerator*    GetWeightGenerator();

  virtual Double_t GetWeight(AliFemtoPair *aPair);
  
 protected:
  AliFemtoModelFreezeOutGenerator *fFreezeOutGenerator;   // Freeze-out coordinates generator
  AliFemtoModelWeightGenerator    *fWeightGenerator;      // Femtoscopic weight generator
  Bool_t                           fCreateCopyHiddenInfo; // Switch to turn on hidden-info generation

 private:
		
#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoModelManager, 1);
  /// \endcond
#endif

    };
  
#endif


