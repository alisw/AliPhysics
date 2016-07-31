///
/// \file AliFemtoModelWeightGeneratorBasic.h
///

#ifndef AliFemtoModelWeightGeneratorBasic_hh
#define AliFemtoModelWeightGeneratorBasic_hh

#include "TRandom2.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelWeightGenerator.h"

/// \class AliFemtoModelWeightGeneratorBasic
/// \brief Basic femtoscopic weight generator only return a simple
/// \author Adam Kisiel kisiel@mps.ohio-state.edu
///
class AliFemtoModelWeightGeneratorBasic : public AliFemtoModelWeightGenerator {
public:
  /// Defaut Constructor - pair type is set to 'None' by default (causes
  /// type lookup for each pair that is called with GenerateWeight)
  ///
  AliFemtoModelWeightGeneratorBasic();
  AliFemtoModelWeightGeneratorBasic(const AliFemtoModelWeightGeneratorBasic &aModel);
  virtual ~AliFemtoModelWeightGeneratorBasic();
  AliFemtoModelWeightGeneratorBasic& operator=(const AliFemtoModelWeightGeneratorBasic &aModel);
  virtual Double_t GenerateWeight(AliFemtoPair *aPair);

  virtual void     SetPairType(Int_t aPairType);
  virtual void     SetPairTypeFromPair(AliFemtoPair *aPair);
  virtual Int_t    GetPairType() const; 

  virtual AliFemtoModelWeightGenerator* Clone() const;
protected:
  
private:
  AliFemtoModelWeightGenerator* GetGenerator() const;

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoModelWeightGeneratorBasic, 1);
  /// \endcond
#endif

};
  
#endif


