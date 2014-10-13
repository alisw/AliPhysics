#ifndef LHAPDF_H
#define LHAPDF_H

#include "LHAPDF/LHAPDFConfig.h"

#include <string>
#include <vector>
#include <iostream>
#include <sstream>

/**
 * @mainpage A C++ wrapper for the LHAPDF library
 *
 * @section intro Introduction
 * The LHAPDF library provides a set of C++ wrapper functions
 * for its Fortran subroutines. New users should browse this
 * documentation and take a look at the @c CCTest1.cc and
 * @c CCTest2.cc example program source files, which are good
 * examples of how the wrapper is used.
 *
 * @section changes Recent changes
 * The LHAPDF wrapper has been improved in several ways for
 * the LHAPDF v5.4 release:
 *
 * @li String passing to the Fortran functions from C++ now
 *  correctly passes the hidden length argument, fixing problems
 *  on 64 bit systems;
 * @li The @c LHAPDFWrap class has been deprecated in favour of
 *  a set of wrapper functions in the @c LHAPDF namespace. The
 *  class interface was misleading, since no persistent state
 *  was involved and two class instances would not have been
 *  independent;
 * @li Proper C++ @c std::string arguments can now be used for
 *  set names: @c char* string arguments can still be passed, due
 *  to implicit conversion via the string(char*) constructor.
 *
 * @section credits Credits
 *
 * @li Originally by Stefan Gieseke.
 * @li Adapted for LHAPDFv4 by Mike Whalley.
 * @li Adapted for LHAPDFv5 by Craig Group/Mike Whalley.
 * @li v5.4: Fortran portability, tidying, extensions
 *  and conversion to namespaced functions by Andy Buckley.
 * @li v5.4.1: Rationalised init functions and deprecated "M"
 *  functions by Andy Buckley.
 * @li v5.5.1: Added PDFSetInfo set metadata struct, and
 *  associated querying based on reading the PDFsets.index file.
 *
 * @example ../examples/CCTest1.cc
 * This is an example of a program using the recommended C++
 * interface to LHAPDF.
 *
 * @example ../examples/CCTest2.cc
 * An example of a program using the C++ interface to LHAPDF to
 * calculate PDF errors.
 */


// Compatibility preprocessing of deprecated "M" function names
#define initPDFSetM initPDFSet
#define initPDFSetByNameM initPDFSetByName
#define initPDFM initPDF
#define initPDFByNameM initPDFByName
#define getDescriptionM getDescription
#define xfxM xfx
#define xfxpM xfxp
#define xfxaM xfxa
#define xfxphotonM xfxphoton
#define numberPDFM numberPDF
#define alphasPDFM alphasPDF
#define getOrderPDFM getOrderPDF
#define getOrderAlphaSM getOrderAlphaS
#define getQMassM getQMass
#define getThresholdM getThreshold
#define getNfM getNf
#define getLam4M getLam4
#define getLam5M getLam5
#define getXminM getXmin
#define getXmaxM getXmax
#define getQ2minM getQ2min
#define getQ2maxM getQ2max


/// Namespace containing all the LHAPDF wrapper functions.
namespace LHAPDF {

  /// @brief Enum of flavours which map to LHAPDF integer codes.
  /// Useful for improving readability of client code. Note that these codes
  /// can't be used to access elements of returned @c vector<double>, which
  /// don't use the LHAPDF scheme (they use "LHAPDF code + 6").
  enum Flavour {
    TBAR= -6, BBAR = -5, CBAR = -4, SBAR = -3, UBAR = -2, DBAR = -1,
    GLUON = 0,
    DOWN = 1, UP = 2, STRANGE = 3, CHARM = 4, BOTTOM = 5, TOP= 6,
    PHOTON = 7
  };

  /// @brief Distinction between evolution or interpolation PDF sets.
  /// Enum to choose whether evolution (i.e. @c LHpdf data file) or
  /// interpolation (i.e. @c LHgrid data file) is used.
  enum SetType {
    EVOLVE = 0, LHPDF = 0,
    INTERPOLATE = 1, LHGRID = 1
  };

  /// Level of noisiness.
  enum Verbosity { SILENT=0, LOWKEY=1, DEFAULT=2 };


  /// @name Global functions
  //@{

  /// Get LHAPDF version string.
  std::string getVersion();

  /// Get max allowed number of concurrent sets.
  int getMaxNumSets();

  /// Global initialisation.
  void initLHAPDF();

  /// Choose level of noisiness.
  void setVerbosity(Verbosity noiselevel);

  /// Extrapolate beyond grid edges.
  void extrapolate(bool extrapolate=true);

  /// Set the LHAPATH variable (the location of the PDF sets directory).
  void setPDFPath(const std::string& path);

  /// Set a steering parameter (direct map to Fortran @c setlhaparm(parm) function).
  void setParameter(const std::string& parm);

  //@}


  /// @name Set metadata
  //@{

  /// Structure containing metadata about a PDF set.
  class PDFSetInfo {
  public:
    std::string file;
    std::string description;
    int id;
    int pdflibNType, pdflibNGroup, pdflibNSet;
    int memberId;
    double lowx, highx;
    double lowQ2, highQ2;

    /// Render a standard representation of a PDF set's metadata.
    std::string toString() const {
      std::ostringstream os;
      os << "PDF set #" << id
         << " {"
         << " file='" << file << "',"
         << " description='" << description << "',"
         << " x = ["  << lowx  << ", " << highx << "],"
         << " Q2 = [" << lowQ2 << ", " << highQ2 << "]"
         << " }";
      return os.str();
    }
  };


  inline std::ostream& operator<<(std::ostream& os, const PDFSetInfo& info) {
    os << info.toString();
    return os;
  }

  /// Get a PDF set info object by filename and member number.
  PDFSetInfo getPDFSetInfo(const std::string& filename, int memid);

  /// Get a PDF set info object by the LHAPDF ID number.
  PDFSetInfo getPDFSetInfo(int id);

  /// Get a vector of PDF set info objects for all known sets.
  std::vector<PDFSetInfo> getAllPDFSetInfo();
  //@}


  /// @name Path info functions
  //@{

  /// Get path to LHAPDF installation (the "prefix" path).
  std::string prefixPath();

  /// Get path to LHAPDF PDF sets directory.
  std::string pdfsetsPath();

  /// Get path to LHAPDF PDF sets index file.
  std::string pdfsetsIndexPath();

  //@}


  /// @name Initialisation functions
  /// LHAPDF functions for initialising PDF sets. If you need to use
  /// more than one set simultaneously, use the multi-set functions, which
  /// have a integer @c nset first argument.
  //@{

  /// Initialise @a member in PDF set @a setid.
  void initPDFSet(int setid, int member);
  /// Initialise @a member in PDF set @a setid (multi-set version).
  void initPDFSet(int nset, int setid, int member); // can't have a default 3rd arg

  /// Initialise @a member in PDF set @a name, of type @a type.
  void initPDFSet(const std::string& name, SetType type, int member=0);
  /// Initialise @a member in PDF set @a name, of type @a type (multi-set version).
  void initPDFSet(int nset, const std::string& name, SetType type, int member=0);

  /// @brief Initialise @a member in PDF set file @a filename.
  /// If @a filename contains a "/" character, it will be used as a path,
  /// otherwise it will be assumed to be a PDF file in the LHAPDF @c PDFsets directory.
  void initPDFSet(const std::string& filename, int member=0);
  /// @brief Initialise @a member in PDF set file @a filename (multi-set version).
  /// If @a filename contains a "/" character, it will be used as a path,
  /// otherwise it will be assumed to be a PDF file in the LHAPDF @c PDFsets directory.
  void initPDFSet(int nset, const std::string& filename, int member=0);

  /// @brief Use @a member in current PDF set.
  /// This operation is computationally cheap.
  void usePDFMember(int member);
  /// @brief Use @a member in PDF set @a nset (multi-set version).
  /// This operation is computationally cheap.
  void usePDFMember(int nset, int member);
  //@}


  /// @name PDF set information
  //@{

  /// Prints a brief description of the current PDF set to stdout.
  void getDescription();
  /// Prints a brief description of the current PDF set to stdout.
  void getDescription(int nset);

  /// Does the current set have a photon member?
  bool hasPhoton();

  /// Number of members available in the current set.
  int numberPDF();
  /// Number of members available in the current set.
  int numberPDF(int nset);

  /// \f$ \alpha_\mathrm{s} \f$ used by the current PDF.
  double alphasPDF(double Q);
  /// \f$ \alpha_\mathrm{s} \f$ used by the current PDF.
  double alphasPDF(int nset, double Q);

  /// Get order at which the PDF was fitted.
  int getOrderPDF();
  /// Get order at which the PDF was fitted.
  int getOrderPDF(int nset);

  /// Perturbative order of parton evolution and \f$ \alpha_\mathrm{s} \f$ respectively.
  int getOrderAlphaS();
  /// Perturbative order of parton evolution and \f$ \alpha_\mathrm{s} \f$ respectively.
  int getOrderAlphaS(int nset);

  /// Quark mass used for flavour @a f.
  double getQMass(int f);
  /// Quark mass used for flavour @a f.
  double getQMass(int nset, int f);

  /// Threshold for flavour @a f.
  double getThreshold(int f);
  /// Threshold for flavour @a f.
  double getThreshold(int nset, int f);

  /// Number of flavours used in the current PDF set.
  int getNf();
  /// Number of flavours used in the current PDF set.
  int getNf(int nset);

  /// Value of QCD \f$ \lambda_4 \f$ for member @a m.
  double getLam4(int m);
  /// Value of QCD \f$ \lambda_4 \f$ for member @a m.
  double getLam4(int nset, int m);

  /// Value of QCD \f$ \lambda_5 \f$ for member @a m.
  double getLam5(int m);
  /// Value of QCD \f$ \lambda_5 \f$ for member @a m.
  double getLam5(int nset, int m);

  /// Minimum \f$ x \f$ value considered valid for this set, as specified by the set authors.
  double getXmin(int m);
  /// Minimum \f$ x \f$ value considered valid for this set, as specified by the set authors.
  double getXmin(int nset, int m);

  /// Maximum \f$ x \f$ value considered valid for this set, as specified by the set authors.
  double getXmax(int m);
  /// Maximum \f$ x \f$ value considered valid for this set, as specified by the set authors.
  double getXmax(int nset, int m);

  /// Minimum \f$ Q^2 \f$ value considered valid for this set, as specified by the set authors.
  double getQ2min(int m);
  /// Minimum \f$ Q^2 \f$ value considered valid for this set, as specified by the set authors.
  double getQ2min(int nset, int m);

  /// Maximum \f$ Q^2 \f$ value considered valid for this set, as specified by the set authors.
  double getQ2max(int m);
  /// Maximum \f$ Q^2 \f$ value considered valid for this set, as specified by the set authors.
  double getQ2max(int nset, int m);
  //@}


  /// @name Nucleon PDFs
  /// These PDFs are defined for protons --- neutron PDFs are usually obtained by
  /// isospin conjugation.
  //@{

  /// Nucleon PDF: returns a vector \f$ x f_i(x, Q) \f$ with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  std::vector<double> xfx(double x, double Q);
  /// Nucleon PDF: returns a vector @c x f_i(x, Q) with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  std::vector<double> xfx(int nset, double x, double Q);

  /// Nucleon PDF: fills primitive 13 element array pointed at by @a results with
  /// \f$ x f(x, Q) \f$ with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  void xfx(double x, double Q, double* results);
  /// Nucleon PDF: fills primitive 13 element array pointed at by @a results with
  /// \f$ x f(x, Q) \f$ with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  void xfx(int nset, double x, double Q, double* results);


  /// Nucleon PDF: returns \f$ x f(x, Q) \f$ for flavour @a fl - this time the flavour encoding
  /// is as in the LHAPDF manual.
  /// @arg -6..-1 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 0 = \f$ g \f$
  /// @arg 1..6 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  double xfx(double x, double Q, int fl);
  /// Nucleon PDF: returns @c x f(x, Q) for flavour @a fl - this time the flavour encoding
  /// is as in the LHAPDF manual.
  /// @arg -6..-1 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 0 = \f$ g \f$
  /// @arg 1..6 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  double xfx(int nset, double x, double Q, int fl);
  //@}


  /// @name Photon PDFs
  //@{

  /// Photon PDF: returns a vector \f$ x f_i(x, Q) \f$ with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  ///
  /// NB. Extra @a P2 and @a ip params.
  std::vector<double> xfxp(double x, double Q, double P2, int ip);
  /// Photon PDF: returns a vector \f$ x f_i(x, Q) \f$ with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  ///
  /// NB. Extra @a P2 and @a ip params.
  std::vector<double> xfxp(int nset, double x, double Q, double P2, int ip);

  /// Photon PDF: fills primitive 13 element array pointed at by @a results with
  /// \f$ x f(x, Q) \f$ with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  ///
  /// NB. Extra @a P2 and @a ip params.
  void xfxp(double x, double Q, double P2, int ip, double* results);
  /// Photon PDF: fills primitive 13 element array pointed at by @a results with
  /// \f$ x f(x, Q) \f$ with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  ///
  /// NB. Extra @a P2 and @a ip params.
  void xfxp(int nset, double x, double Q, double P2, int ip, double* results);


  /// Photon PDF: returns \f$ x f(x, Q) \f$ for flavour @a fl - this time the flavour encoding
  /// is as in the LHAPDF manual.
  /// @arg -6..-1 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 0 = \f$ g \f$
  /// @arg 1..6 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  ///
  /// NB. Extra @a P2 and @a ip params.
  double xfxp(double x, double Q, double P2, int ip, int fl);
  /// Photon PDF: returns \f$ x f(x, Q) \f$ for flavour @a fl - this time the flavour encoding
  /// is as in the LHAPDF manual.
  /// @arg -6..-1 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 0 = \f$ g \f$
  /// @arg 1..6 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  ///
  /// NB. Extra @a P2 and @a ip params.
  double xfxp(int nset, double x, double Q, double P2, int ip, int fl);
  //@}


  /// @name Nuclear PDFs
  //@{

  /// Nuclear PDF: returns a vector \f$ x f_i(x, Q) \f$ with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  ///
  /// NB. Extra @a a param for atomic mass number.
  std::vector<double> xfxa(double x, double Q, double a);
  /// Nuclear PDF: returns a vector \f$ x f_i(x, Q) \f$ with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  ///
  /// NB. Extra @a a param for atomic mass number.
  std::vector<double> xfxa(int nset, double x, double Q, double a);

  /// Nuclear PDF: fills primitive 13 element array pointed at by @a results with
  /// \f$ x f(x, Q) \f$ with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  ///
  /// NB. Extra @a a param for atomic mass number.
  void xfxa(double x, double Q, double a, double* results);
  /// Nuclear PDF: fills primitive 13 element array pointed at by @a results with
  /// \f$ x f(x, Q) \f$ with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  ///
  /// NB. Extra @a a param for atomic mass number.
  void xfxa(int nset, double x, double Q, double a, double* results);

  /// Nuclear PDF: returns \f$ x f(x, Q) \f$ for flavour @a fl - this time the flavour encoding
  /// is as in the LHAPDF manual.
  /// @arg -6..-1 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 0 = \f$ g \f$
  /// @arg 1..6 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  ///
  /// NB. Extra @a a param for atomic mass number.
  double xfxa(double x, double Q, double a, int fl);
  /// Nuclear PDF: returns \f$ x f(x, Q) \f$ for flavour @a fl - this time the flavour encoding
  /// is as in the LHAPDF manual.
  /// @arg -6..-1 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 0 = \f$ g \f$
  /// @arg 1..6 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  ///
  /// NB. Extra @a a param for atomic mass number.
  double xfxa(int nset, double x, double Q, double a, int fl);
  //@}


  /// @name Nucleon MRST QED PDF
  /// These functions only apply to the MRST QED PDF set, since they return an extra element
  /// for the additional photon.
  //@{

  /// MRST QED PDF: returns a vector \f$ x f_i(x, Q) \f$ with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$;
  /// @arg 13 = \f$ \gamma \f$.
  ///
  /// NB. Note extra element in this set for MRST photon.
  std::vector<double> xfxphoton(double x, double Q);
  /// MRST QED PDF: returns a vector \f$ x f_i(x, Q) \f$ with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$;
  /// @arg 13 = \f$ \gamma \f$.
  std::vector<double> xfxphoton(int nset, double x, double Q);


  /// MRST QED PDF: fills primitive 14 element array pointed at by @a results with
  /// \f$ x f(x, Q) \f$ with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  /// @arg 13 = \f$ \gamma \f$.
  ///
  /// NB. Note extra element in this set for MRST photon.
  void xfxphoton(double x, double Q, double* results);
  /// MRST QED PDF: fills primitive 14 element array pointed at by @a results with
  /// \f$ x f(x, Q) \f$ with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  /// @arg 13 = \f$ \gamma \f$.
  ///
  /// NB. Note extra element in this set for MRST photon.
  void xfxphoton(int nset, double x, double Q, double* results);


  /// MRST QED PDF: returns \f$ x f(x, Q) \f$ for flavour @a fl - this time the flavour encoding
  /// is as in the LHAPDF manual.
  /// @arg -6..-1 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 0 = \f$ g \f$
  /// @arg 1..6 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$;
  /// @arg 7 = \f$ \gamma \f$.
  ///
  /// NB. Note extra element in this set for MRST photon.
  double xfxphoton(double x, double Q, int fl);
  /// MRST QED PDF: returns \f$ x f(x, Q) \f$ for flavour @a fl - this time the flavour encoding
  /// is as in the LHAPDF manual.
  /// @arg -6..-1 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 0 = \f$ g \f$
  /// @arg 1..6 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$;
  /// @arg 7 = \f$ \gamma \f$.
  double xfxphoton(int nset, double x, double Q, int fl);
  //@}


  /// @name Deprecated initialisation functions
  /// LHAPDF functions for initialising PDF sets. If you need to use
  /// more than one set simultaneously, use the multi-set functions, which
  /// have a integer @c nset first argument.
  /// @deprecated These init methods are deprecated.
  //@{

  /// The PDF set by file path, see subdir @c PDFsets of LHAPDF for choices.
  //void initPDFSet(const std::string& path);
  /// The PDF set by file path, see subdir @c PDFsets of LHAPDF for choices.
  //void initPDFSet(int nset, const std::string& path);

  /// The PDF set by name and type, see subdir @c PDFsets of LHAPDF for choices.
  void initPDFSetByName(const std::string& name, SetType type);
  /// The PDF set by name and type, see subdir @c PDFsets of LHAPDF for choices.
  void initPDFSetByName(int nset, const std::string& name, SetType type);

  /// The PDF set by filename, see subdir @c PDFsets of LHAPDF for choices.
  void initPDFSetByName(const std::string& filename);
  /// The PDF set by filename, see subdir @c PDFsets of LHAPDF for choices.
  void initPDFSetByName(int nset, const std::string& filename);

  /// The choice of PDF member out of one distribution.
  void initPDF(int memset);
  /// The choice of PDF member out of one distribution.
  void initPDF(int nset, int memset);

  /// @brief Convenient initializer with PDF set @a name, set type @a type and member @a memset.
  /// @param name The name of the desired set.
  /// @param type The type of PDF set (grid or data) by enum.
  /// @param memset PDF number within set @a name.
  /// Equivalent to @c initPDFSetByName + @c initPDF.
  void initPDFByName(const std::string& name, SetType type, int memset);

  /// @brief Typical initializer for multiple PDF sets with PDF set @a name and member @a memset.
  /// @param nset Specifies the reference number for the set to be initialized.
  /// @param name Name of the desired set.
  /// @param type The type of PDF set (grid or data) by enum.
  /// @param memset PDF number within set @a name.
  /// Equivalent to @c initPDFSetByNameM + @c initPDFM.
  void initPDFByName(int nset, const std::string& name, SetType type, int memset);

  /// @brief Convenient initializer with PDF set @a filename and member @a memset.
  /// @param filename The name of the grid or data file of the desired set.
  /// @param memset PDF number within set @a name.
  /// Equivalent to @c initPDFSetByName + @c initPDF.
  void initPDFByName(const std::string& filename, int memset);
  /// @brief Typical initializer for multiple PDF sets with PDF set @a name and member @a memset.
  /// @param nset Specifies the reference number for the set to be initialized.
  /// @param filename Name of the grid or data file of the desired set.
  /// @param memset PDF number within set @a name.
  /// Equivalent to @c initPDFSetByNameM + @c initPDFM.
  void initPDFByName(int nset, const std::string& filename, int memset);
  //@}


}

#endif
