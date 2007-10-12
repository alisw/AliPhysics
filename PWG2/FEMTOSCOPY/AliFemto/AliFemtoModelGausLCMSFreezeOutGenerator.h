////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelGausLCMSFreezeOutGenerator - freeze-out                     ///
/// coordinates generator, generating a 3D gaussian ellipsoid in LCMS        ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoModelGausLCMSFreezeOutGenerator_hh
#define AliFemtoModelGausLCMSFreezeOutGenerator_hh

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
  Double_t fSizeOut;
  Double_t fSizeSide;
  Double_t fSizeLong;

 private:
  AliFemtoModelFreezeOutGenerator* GetGenerator() const;
		
#ifdef __ROOT__
  ClassDef(AliFemtoModelGausLCMSFreezeOutGenerator, 1)
#endif

    };
  
#endif


