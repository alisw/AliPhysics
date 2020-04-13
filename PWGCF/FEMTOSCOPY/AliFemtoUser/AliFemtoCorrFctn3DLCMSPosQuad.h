///
/// \file AliFemtoUser/AliFemtoCorrFctn3DLCMSPosQuad.h
///



#ifndef ALIFEMTOCORRFCTN3DLCMSPOSQUAD_H
#define ALIFEMTOCORRFCTN3DLCMSPOSQUAD_H

// forward declare classes
class TH3F;
class AliFemtoPairCut;

// include headers
#include "AliFemtoCorrFctn.h"

/// \class AliFemtoCorrFctn3DLCMSPosQuad
/// \brief Calculates the 3D correlation function for pairs of identical
///        particles in Bertsh-Pratt coordinates. All counts are stored
///        in the positive q_out quadrant - symmetry is assumed.
///
/// Use the AliFemtoCorrFctn3DLCMSSym or AliFemtoBPLCMS3DCorrFctn for more
/// general purpose correlation functions - this is designed for saving space
///
///
class AliFemtoCorrFctn3DLCMSPosQuad : public AliFemtoCorrFctn {
public:

  /// \class AliFemtoCorrFctn3DLCMSPosQuad::Parameters
  /// \brief Simple struct used to construct correlation function objects
  struct Parameters {

    double QoHi,
           QsHi,
           QlHi;

    size_t nbins_out,
           nbins_side,
           nbins_long;

    TString prefix; // name prefix
    TString suffix; // name suffix

    /// Default: 5MeV bins between +- 160 MeV
    Parameters()
      : QoHi(0.16)
      , QsHi(0.16)
      , QlHi(0.16)
      , nbins_out(QoHi / 0.005)
      , nbins_side(QsHi * 2 / 0.005 + 1)
      , nbins_long(QlHi * 2 / 0.005 + 1)
      , prefix("")
      , suffix("")
      { }

    AliFemtoCorrFctn3DLCMSPosQuad into() const
      { return AliFemtoCorrFctn3DLCMSPosQuad(*this); }

    AliFemtoCorrFctn3DLCMSPosQuad* into_ptr() const
      { return new AliFemtoCorrFctn3DLCMSPosQuad(*this); }

    #define SETTER(__name, __target, __type) \
      Parameters __name(const __type &v) \
        { Parameters p(*this); p.__target = v; return p; }

    SETTER(WithSuffix, suffix, TString);
    SETTER(WithPrefix, prefix, TString);
    SETTER(WithNbinsOut, nbins_out, int);
    SETTER(WithNbinsSide, nbins_side, int);
    SETTER(WithNbinsLong, nbins_long, int);

    #undef SETTER

    Parameters WithMaxQ(double q)
      {
        Parameters p(*this);
        p.QoHi = q;
        p.QsHi = q;
        p.QlHi = q;
        return p;
      }

    Parameters WithMaxQ(double qo, double qsl)
      {
        Parameters p(*this);
        p.QoHi = qo;
        p.QsHi = qsl;
        p.QlHi = qsl;
        return p;
      }

    /// Recalcuate number of bins such that this bin size (in GeV) is
    /// approximated in each direction
    Parameters WithBinSize(double binsize)
      {
        Parameters p(*this);
        p.nbins_out = p.QoHi / binsize;
        p.nbins_side = p.QsHi * 2 / binsize;
        p.nbins_long = p.QlHi * 2 / binsize;
        return p;
      }

    Parameters WithBins(int nbins)
      {
        Parameters p(*this);
        p.nbins_out = nbins / 2;
        p.nbins_side = nbins;
        p.nbins_long = nbins;
        return p;
      }

    Parameters WithBins(int nbins_o, int nbins_sl)
      {
        Parameters p(*this);
        p.nbins_out = nbins_o;
        p.nbins_side = nbins_sl;
        p.nbins_long = nbins_sl;
        return p;
      }

  };

  static Parameters Build()
    { return Parameters(); }

  AliFemtoCorrFctn3DLCMSPosQuad(const Parameters &);

  /// Symmetric bins (Qout uses nbins / 2)
  ///
  /// \param title The title with which to give the output
  /// \param nbins The number of bins in each direction of , and q
  ///
  AliFemtoCorrFctn3DLCMSPosQuad(const TString &prefix, const int nbins, const float QHi);

  /// Big constructor - The "Parameters" should be prefered
  AliFemtoCorrFctn3DLCMSPosQuad(const TString &prefix,
                                const TString &suffix,
                                const size_t nbins_out,
                                const size_t nbins_side,
                                const size_t nbins_long,
                                const float QoHi,
                                const float QsHi,
                                const float QlHi);

  /// Variable Bin Size
  AliFemtoCorrFctn3DLCMSPosQuad(const TString &prefix,
                                const TString &suffix,
                                const std::vector<double> &out_bins,
                                const std::vector<double> &side_bins,
                                const std::vector<double> &long_bins);

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
