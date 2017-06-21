///
/// \file AliFemtoModelWeightGeneratorBasic.h
///

#ifndef AliFemtoModelWeightGeneratorBasic_hh
#define AliFemtoModelWeightGeneratorBasic_hh

#include "TRandom2.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelWeightGenerator.h"

/// \class AliFemtoModelWeightGeneratorBasic
/// \brief Basic femtoscopic weight generator to determine a
///        femtoscopic weight using quantum statistics
///
/// Generates a simple femtoscopic weight coming only from quantum
/// statistics - symmetrization or anti-symmetrization of the pair
/// wave function.
///
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

  /// Sets whether or not to print a notification when the generator
  /// is asked to find weight of a pair with zero mass, energy, or
  /// momentum (i.e. "empty")
  ///
  void ShouldPrintEmptyParticleNotification(bool option);

  virtual AliFemtoModelWeightGenerator* Clone() const;

protected:

  /// Internal flag to determine whether or not this generator should
  /// a print notification when asked to generate weight for "empty"
  /// AliFemtoPair
  bool fPrintEmptyParticleNotification;

private:
  AliFemtoModelWeightGenerator* GetGenerator() const;

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoModelWeightGeneratorBasic, 1);
  /// \endcond
#endif
};

inline void AliFemtoModelWeightGeneratorBasic::ShouldPrintEmptyParticleNotification(bool option)
{
  fPrintEmptyParticleNotification = option;
}

#endif


