///
/// \file AliFemtoUser/AliFemtoCorrFctn3DLCMSPosQuad.h
///

#ifndef ALIFEMTOCORRFCTN3DLCMS_H
#define ALIFEMTOCORRFCTN3DLCMS_H

// forward declare classes
class TH3F;
class AliFemtoPairCut;

// include headers
#include "AliFemtoCorrFctn.h"

/// \class AliFemtoCorrFctn3DLCMSPosQuad
/// \brief Calculates the 3D correlation function for pairs of identical
///        particles in Bertsh-Pratt coordinates. All counts are stored
///        in the positive quadrant - symmetry is assumed.
///
/// Use the AliFemtoCorrFctn3DLCMSSym or AliFemtoBPLCMS3DCorrFctn for more
/// general purpose correlation functions - this is designed for
///
///
class AliFemtoCorrFctn3DLCMSPosQuad : public AliFemtoCorrFctn {
public:

  /// \class AliFemtoCorrFctn3DLCMSPosQuad::Parameters
  /// \brief Simple struct used to construct correlation function objects
  struct Parameters {
    float QHi;    // 0.5 GeV
    int nbins;    // 100
    TString name; // ""

    Parameters();
  };

  AliFemtoCorrFctn3DLCMSPosQuad(const Parameters &);

  /// Build the correlation function with parameters.
  ///
  /// \param title The title with which to give the output
  /// \param nbins The number of bins in each direction of , and q
  ///
  AliFemtoCorrFctn3DLCMSPosQuad(const TString &name, const int nbins, const float QHi);

  /// Copy Constructor
  AliFemtoCorrFctn3DLCMSPosQuad(const AliFemtoCorrFctn3DLCMSPosQuad& aCorrFctn);

  /// Deletes histograms
  virtual ~AliFemtoCorrFctn3DLCMSPosQuad();

  AliFemtoCorrFctn3DLCMSPosQuad& operator=(const AliFemtoCorrFctn3DLCMSPosQuad& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* pair) { AddRealPair(*pair); }
  virtual void AddRealPair(const AliFemtoPair &);

  virtual void AddMixedPair(AliFemtoPair* pair) { AddMixedPair(*pair); }
  virtual void AddMixedPair(const AliFemtoPair &);

  virtual void Finish();

  virtual TList* GetOutputList();

  virtual AliFemtoCorrFctn* Clone() const;

protected:

  TString name;  ///< Name of the correlation function

  TH3F* fNumerator;    ///< Numerator
  TH3F* fDenominator;  ///< Denominator
  TH3F* fQinvWeight;   ///< Qinv-Weighted denominator
};

inline AliFemtoCorrFctn* AliFemtoCorrFctn3DLCMSPosQuad::Clone() const
{
  return new AliFemtoCorrFctn3DLCMSPosQuad(*this);
}

inline void AliFemtoCorrFctn3DLCMSPosQuad::Finish()
{}

#endif
