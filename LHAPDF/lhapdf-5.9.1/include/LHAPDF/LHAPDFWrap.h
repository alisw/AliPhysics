#ifndef LHAPDFWRAP_H
#define LHAPDFWRAP_H

#include <string>
#include <vector>

// WARNING! This C++ interface is deprecated in favour of the 
// interface declared in LHAPDF/LHAPDF.h

// This class is a wrapper around the LHAPDF package for parton
// distribution functions of the proton.
//
// Originally by Stefan Gieseke.
// Adapted for LHAPDFv4 by Mike Whalley.
// Adapted for LHAPDFv5 by Craig Group/Mike Whalley.
// Fortran portability and interface improvements by Andy Buckley.


/////////////////////////////////////////////////////////////////


/// Wrapper class used to contain all the wrapper functions.
/// @deprecated 
/// The class-based C++ wrapper on LHAPDF will be retired in a forthcoming
/// release of LHAPDF in favour of the namespaced wrapper declared in @c
/// LHAPDF.h. Please convert client code which uses the class interface to use
/// the new interface instead.  Typically, this will just involve changing the
/// header include from @c "LHAPDF/LHAPDFWrap.h" to @c "LHAPDF/LHAPDFWrap.h",
/// changing any constructors to initialisation functions, and replacing @c
/// LHAPDFWrap objects with an @c LHAPDF namespace. For example,
/// @code
///  #include "LHAPDF/LHAPDFWrap.h"
///  LHAPDFWrap pdf = LHAPDFWrap("MRST2004qed.LHgrid", 0);
///  pdf.getDescription();
/// @endcode
/// would be replaced by
/// @code
///  #include "LHAPDF/LHAPDF.h"
///  LHAPDF::initPDFByName("MRST2004qed.LHgrid", 0);
///  LHAPDF::getDescription();
/// @endcode
class LHAPDFWrap {

public:
  /// Do-nothing constructor.
  LHAPDFWrap();
  
  /// Typical constructor with PDF set 'name' 
  /// 'name' is the name of the grid or data file of the desired set.  
  LHAPDFWrap(const std::string& name);

  /// Typical constructor with PDF set 'name' and subset 'memset' 
  /// 'name' is the name of the grid or data file of the desired set.  
  LHAPDFWrap(const std::string& name, int memset);

  /// Typical constructor (when multiple PDF sets need to be initialized) 
  /// with pdfset 'name' and subset 'memset'.
  /// 'name' is the name of the grid or data file of the desired set.  
  /// int nset specifies the reference number for the set to be initialized
  LHAPDFWrap(int nset, const std::string& name);

  /// Typical constructor (when multiple PDF sets need to be initialized) 
  /// with PDF set 'name' and subset 'memset'. 
  /// 'name' is the name of the grid or data file of the desired set.  
  /// int nset specifies the reference number for the set to be initialized
  LHAPDFWrap(int nset, const std::string& name, int memset);

  /// Returns a vector xf(x, Q) with index 0 < i < 12.
  /// 0..5 = tbar, ..., ubar, dbar; 
  /// 6 = g; 
  /// 7..12 = d, u, ..., t
  std::vector<double> xfx(double x, double Q);

  /// Returns xf(x, Q) for flavour fl - this time the flavour encoding
  /// is as in the LHAPDF manual...
  /// -6..-1 = tbar,...,ubar, dbar
  /// 1..6 = duscbt
  /// 0 = g
  double xfx(double x, double Q, int fl);

  std::vector<double> xfxp(double x, double Q, double P2, int ip);
  double xfxp(double x, double Q, double P2, int ip, int fl);

  std::vector<double> xfxa(double x, double Q, double a);
  double xfxa(double x, double Q, double a, int fl);

  std::vector<double> xfxphoton(double x, double Q);
  double xfxphoton(double x, double Q, int fl);

  /// The PDF set by name, see subdir 'PDFset' of LHAPDFv2 for choices
  void initPDFSet(const std::string& name);

  /// The PDF set by name, see subdir 'PDFset' of LHAPDFv2 for choices
  void initPDFSetByName(const std::string& name);

  /// The choice of PDF subset out of one distribution
  void initPDF(int memset);

  /// Prints a brief description of the current pdf set to stdout
  void getDescription();

  /// Number of subsets available in the current distribution.
  int numberPDF();

  /// \f$ \alpha_\mathrm{s} \f$ used by the current PDF.
  double alphasPDF(double Q);

  int getOrderPDF();

  /// Perturbative order of parton evolution and \f$ \alpha_\mathrm{s} \f$ respectively.
  int getOrderAlphaS();

  /// Quark mass used for flavour f.
  double getQMass(int f);

  /// Threshold for flavour f.
  double getThreshold(int f);

  /// Number of flavours used in the current PDF set.
  int getNf();

  /// Value of QCD lambda4 for member m 
  double getLam4(int m);

  /// Value of QCD lambda5 for member m 
  double getLam5(int m);

  double getXmin(int m);
  double getXmax(int m);
  double getQ2min(int m);
  double getQ2max(int m);

  void extrapolate();

  // Additional functions for when more than 1 PDF set is being stored in memory

  // Returns a vector xf(x, Q) with index 0 < i < 12.
  // 0..5 = tbar, ..., ubar, dbar; 
  // 6 = g; 
  // 7..12 = d, u, ..., t
  std::vector<double> xfxM(int nset, double x, double Q);

  // Returns xf(x, Q) for flavour fl - this time the flavour encoding
  // is as in the LHAPDF manual...
  // -6..-1 = tbar,...,ubar, dbar
  // 1..6 = duscbt
  // 0 = g
  double xfxM(int nset, double x, double Q, int fl);

  std::vector<double> xfxpM(int nset, double x, double Q, double P2, int ip);
  double xfxpM(int nset, double x, double Q, double P2, int ip, int fl);

  std::vector<double> xfxaM(int nset, double x, double Q, double a);
  double xfxaM(int nset, double x, double Q, double a, int fl);

  std::vector<double> xfxphotonM(int nset, double x, double Q);
  double xfxphotonM(int nset, double x, double Q, int fl);

  /// The PDF set by name, see subdir 'PDFset' of LHAPDFv2 for choices
  void initPDFSetM(int nset, const std::string& name);

  /// The PDF set by name, see subdir 'PDFset' of LHAPDFv2 for choices
  void initPDFSetByNameM(int nset, const std::string& name);

  /// The choice of PDF subset out of one distribution
  void initPDFM(int nset, int memset);

  /// Prints a brief description of the current PDF set to stdout
  void getDescriptionM(int nset);

  /// Number of subsets available in the current distribution.
  int numberPDFM(int nset);

  /// \f$ \alpha_\mathrm{s} \f$ used by the current PDF.
  double alphasPDFM(int nset, double Q);

  int getOrderPDFM(int nset);

  /// Perturbative order of parton evolution and \f$ \alpha_\mathrm{s} \f$ respectively.
  int getOrderAlphaSM(int nset);

  /// Quark mass used for flavour f.
  double getQMassM(int nset, int f);

  /// Threshold for flavour f.
  double getThresholdM(int nset, int f);

  /// Number of flavours used in the current PDF set.
  int getNfM(int nset);

  /// Value of QCD lambda4 for member m 
  double getLam4M(int nset, int m);

  /// Value of QCD lambda5 for member m 
  double getLam5M(int nset, int m);

  double getXminM(int nset, int m);
  double getXmaxM(int nset, int m);
  double getQ2minM(int nset, int m);
  double getQ2maxM(int nset, int m);

};

#endif
