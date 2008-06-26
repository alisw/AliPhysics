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

  void SetSelectPrimaryFromHidden(bool aUse);
  Bool_t GetSelectPrimaryFromHidden();

  void SetSizeInv(Double_t aSizeInv);
  
  Double_t GetSizeInv() const;

  virtual AliFemtoModelFreezeOutGenerator* Clone() const;

 protected:
  Double_t fSizeInv;        // Size of the source
  Bool_t fSelectPrimary;    // If set to true, the existing hidden info is assumed
                            // to contain the particle creation point (in cm)
                            // and the model will try to guess whether the particle
                            // is primary based on that and assign creation point
                            // only for primary particles

 private:
  AliFemtoModelFreezeOutGenerator* GetGenerator() const;
		
#ifdef __ROOT__
  ClassDef(AliFemtoModelGausRinvFreezeOutGenerator, 1)
#endif

    };
  
#endif


