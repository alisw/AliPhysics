////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelGausLCMSFreezeOutGenerator - freeze-out                     ///
/// coordinates generator, generating a 3D gaussian ellipsoid in LCMS        ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOMODELGAUSLCMSFREEZEOUTGENERATOR_H
#define ALIFEMTOMODELGAUSLCMSFREEZEOUTGENERATOR_H

#include "AliFemtoModelFreezeOutGenerator.h"

#include "TRandom.h"

class AliFemtoModelGausLCMSFreezeOutGenerator : public AliFemtoModelFreezeOutGenerator
{
 public:
  AliFemtoModelGausLCMSFreezeOutGenerator();
  AliFemtoModelGausLCMSFreezeOutGenerator(const AliFemtoModelGausLCMSFreezeOutGenerator &aModel);
  virtual ~AliFemtoModelGausLCMSFreezeOutGenerator();
  virtual void GenerateFreezeOut(AliFemtoPair *aPair);

  void SetSizeOut(Double_t aSizeOut);
  void SetSizeSide(Double_t aSizeSide);
  void SetSizeLong(Double_t aSizeLong);
  
  Double_t GetSizeOut() const;
  Double_t GetSizeSide() const;
  Double_t GetSizeLong() const;

  virtual AliFemtoModelFreezeOutGenerator* Clone() const;

 protected:
  Double_t fSizeOut;  // Size of the source in the out direction
  Double_t fSizeSide; // Size of the source in the side direction
  Double_t fSizeLong; // Size of the source in the long direction

 private:
  AliFemtoModelFreezeOutGenerator* GetGenerator() const;
		
#ifdef __ROOT__
  ClassDef(AliFemtoModelGausLCMSFreezeOutGenerator, 1)
#endif

    };
  
#endif


