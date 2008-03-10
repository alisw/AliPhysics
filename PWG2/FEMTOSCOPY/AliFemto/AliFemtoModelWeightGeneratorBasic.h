////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelWeightGeneratorBasic -  basic femtoscopic weight generator  ///
/// only return a simple                                                          ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoModelWeightGeneratorBasic_hh
#define AliFemtoModelWeightGeneratorBasic_hh

#include "TRandom2.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelWeightGenerator.h"

class AliFemtoModelWeightGeneratorBasic : public AliFemtoModelWeightGenerator
{
 public:
  AliFemtoModelWeightGeneratorBasic();
  AliFemtoModelWeightGeneratorBasic(const AliFemtoModelWeightGeneratorBasic &aModel);
  virtual ~AliFemtoModelWeightGeneratorBasic();
  virtual Double_t GenerateWeight(AliFemtoPair *aPair);

  virtual void     SetPairType(Int_t aPairType);
  virtual void     SetPairTypeFromPair(AliFemtoPair *aPair);
  virtual Int_t    GetPairType() const; 

  virtual AliFemtoModelWeightGenerator* Clone() const;
 protected:
  
 private:
  AliFemtoModelWeightGenerator* GetGenerator() const;

#ifdef __ROOT__
  ClassDef(AliFemtoModelWeightGeneratorBasic, 1)
#endif

    };
  
#endif


