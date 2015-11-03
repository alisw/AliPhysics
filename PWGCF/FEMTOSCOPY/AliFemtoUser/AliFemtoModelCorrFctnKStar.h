///
/// \file AliFemtoModelCorrFctnKStar.h
/// \author Andrew Kubera
///

#ifndef ALIFEMTOMODELCORRFCTNKSTAR_H
#define ALIFEMTOMODELCORRFCTNKSTAR_H

class AliFemtoPair;

#include "AliFemtoModelCorrFctn.h"

/// \class AliFemtoModelCorrFctnKStar
/// \brief The correlation function which plots numerators and denominator plots from
///        the real monte caro 
///
/// \author: Andrew Kubera, andrew.kubera@cern.ch
///
///
class AliFemtoModelCorrFctnKStar : public AliFemtoModelCorrFctn {
public:

  /**
   * Default constructor
   * 
   */
  AliFemtoModelCorrFctnKStar();
  
  /**
   * Construct with histogram parameters
   */
  AliFemtoModelCorrFctnKStar(const char *title, Int_t aNbins, Double_t aQinvLo, Double_t aQinvHi);
  
  /// Copy Constructor
  AliFemtoModelCorrFctnKStar(const AliFemtoModelCorrFctnKStar&);
  
  /// Destructor
  virtual ~AliFemtoModelCorrFctnTrueQ();
  
  /// Return information about the run of the correlation function
  virtual AliFemtoString Report();

  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPir);

protected:

  TH1F *fTrueNum;   ///< Numerator in KStar
  TH1F *fTrueDen;   ///< Denominator in Kstar

private:

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoModelCorrFctnKStar, 0);
  /// \endcond
#endif

};

#endif
