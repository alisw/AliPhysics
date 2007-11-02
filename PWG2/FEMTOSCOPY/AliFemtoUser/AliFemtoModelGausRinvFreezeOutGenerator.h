////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelGausRinvFreezeOutGenerator - freeze-out                     ///
/// coordinates generator, generating a 3D gaussian spheroid in PRF          ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoModelGausRinvFreezeOutGenerator_hh
#define AliFemtoModelGausRinvFreezeOutGenerator_hh

#include "AliFemtoModelFreezeOutGenerator.h"

#include "TRandom.h"

class AliFemtoModelGausRinvFreezeOutGenerator : public AliFemtoModelFreezeOutGenerator
{
 public:
  AliFemtoModelGausRinvFreezeOutGenerator();
  AliFemtoModelGausRinvFreezeOutGenerator(const AliFemtoModelGausRinvFreezeOutGenerator &aModel);
  virtual ~AliFemtoModelGausRinvFreezeOutGenerator();
  virtual void GenerateFreezeOut(AliFemtoPair *aPair);

  void SetSizeInv(Double_t aSizeInv);
  
  Double_t GetSizeInv() const;

  virtual AliFemtoModelFreezeOutGenerator* Clone() const;

 protected:
  Double_t fSizeInv;

 private:
  AliFemtoModelFreezeOutGenerator* GetGenerator() const;
		
#ifdef __ROOT__
  ClassDef(AliFemtoModelGausRinvFreezeOutGenerator, 1)
#endif

    };
  
#endif


