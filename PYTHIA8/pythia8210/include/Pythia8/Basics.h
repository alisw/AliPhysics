// Basics.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for basic, often-used helper classes.
// RndmEngine: base class for external random number generators.
// Rndm: random number generator.
// Vec4: simple four-vectors.
// RotBstMatrix: matrices encoding rotations and boosts of Vec4 objects.
// Hist: simple one-dimensional histograms.

#ifndef Pythia8_Basics_H
#define Pythia8_Basics_H

#include "Pythia8/PythiaStdlib.h"

namespace Pythia8 {

//==========================================================================

// RndmEngine is the base class for external random number generators.
// There is only one pure virtual method, that should do the generation.

class RndmEngine {

public:

  // Destructor.
  virtual ~RndmEngine() {}

  // A pure virtual method, wherein the derived class method
  // generates a random number uniformly distributed between 1 and 1.
  virtual double flat() = 0;

};

//==========================================================================

// Rndm class.
// This class handles random number generation according to the
// Marsaglia-Zaman-Tsang algorithm.

class Rndm {

public:

  // Constructors.
  Rndm() : initRndm(false), seedSave(0), sequence(0),
    useExternalRndm(false), rndmEngPtr(0) { }
  Rndm(int seedIn) : initRndm(false), seedSave(0), sequence(0),
    useExternalRndm(false), rndmEngPtr(0) { init(seedIn);}

  // Possibility to pass in pointer for external random number generation.
  bool rndmEnginePtr( RndmEngine* rndmEngPtrIn);

  // Initialize, normally at construction or in first call.
  void init(int seedIn = 0) ;

  // Generate next random number uniformly between 0 and 1.
  double flat() ;

  // Generate random numbers according to exp(-x).
  double exp() { return -log(flat()) ;}

  // Generate random numbers according to x * exp(-x).
  double xexp() { return -log(flat() * flat()) ;}

  // Generate random numbers according to exp(-x^2/2).
  double gauss() {return sqrt(-2. * log(flat())) * cos(M_PI * flat());}

  // Generate two random numbers according to exp(-x^2/2-y^2/2).
  pair<double, double> gauss2() {double r = sqrt(-2. * log(flat()));
    double phi = 2. * M_PI * flat();
    return pair<double, double>(r * sin(phi), r * cos(phi));}

  // Pick one option among  vector of (positive) probabilities.
  int pick(const vector<double>& prob) ;

  // Save or read current state to or from a binary file.
  bool dumpState(string fileName);
  bool readState(string fileName);

private:

  // Default random number sequence.
  static const int DEFAULTSEED;

  // State of the random number generator.
  bool   initRndm;
  int    i97, j97, seedSave;
  long   sequence;
  double u[97], c, cd, cm;

  // Pointer for external random number generation.
  bool   useExternalRndm;
  RndmEngine* rndmEngPtr;

};

//==========================================================================

// Forward reference to RotBstMatrix class; needed in Vec4 class.
class RotBstMatrix;

//--------------------------------------------------------------------------

// Vec4 class.
// This class implements four-vectors, in energy-momentum space.
// (But can equally well be used to hold space-time four-vectors.)

class Vec4 {

public:

  // Constructors.
  Vec4(double xIn = 0., double yIn = 0., double zIn = 0., double tIn = 0.)
    : xx(xIn), yy(yIn), zz(zIn), tt(tIn) { }
  Vec4(const Vec4& v) : xx(v.xx), yy(v.yy), zz(v.zz), tt(v.tt) { }
  Vec4& operator=(const Vec4& v) { if (this != &v) { xx = v.xx; yy = v.yy;
    zz = v.zz; tt = v.tt; } return *this; }
  Vec4& operator=(double value) { xx = value; yy = value; zz = value;
    tt = value; return *this; }

  // Member functions for input.
  void reset() {xx = 0.; yy = 0.; zz = 0.; tt = 0.;}
  void p(double xIn, double yIn, double zIn, double tIn)
    {xx = xIn; yy = yIn; zz = zIn; tt = tIn;}
  void p(Vec4 pIn) {xx = pIn.xx; yy = pIn.yy; zz = pIn.zz; tt = pIn.tt;}
  void px(double xIn) {xx = xIn;}
  void py(double yIn) {yy = yIn;}
  void pz(double zIn) {zz = zIn;}
  void e(double tIn) {tt = tIn;}

  // Member functions for output.
  double px() const {return xx;}
  double py() const {return yy;}
  double pz() const {return zz;}
  double e() const {return tt;}
  double& operator[](int i) {
    if      (i == 1) return xx;
    else if (i == 2) return yy;
    else if (i == 3) return zz;
    else             return tt;
  }
  double mCalc() const {double temp = tt*tt - xx*xx - yy*yy - zz*zz;
    return (temp >= 0.) ? sqrt(temp) : -sqrt(-temp);}
  double m2Calc() const {return tt*tt - xx*xx - yy*yy - zz*zz;}
  double pT() const {return sqrt(xx*xx + yy*yy);}
  double pT2() const {return xx*xx + yy*yy;}
  double pAbs() const {return sqrt(xx*xx + yy*yy + zz*zz);}
  double pAbs2() const {return xx*xx + yy*yy + zz*zz;}
  double eT() const {double temp = xx*xx + yy*yy;
    return tt * sqrt( temp / (temp + zz*zz) );}
  double eT2() const {double temp = xx*xx + yy*yy;
    return tt*tt * temp / (temp + zz*zz);}
  double theta() const {return atan2(sqrt(xx*xx + yy*yy), zz);}
  double phi() const {return atan2(yy,xx);}
  double thetaXZ() const {return atan2(xx,zz);}
  double pPos() const {return tt + zz;}
  double pNeg() const {return tt - zz;}
  double rap() const {return 0.5 * log( (tt + zz) / (tt - zz) );}
  double eta() const {double xyz = sqrt(xx*xx + yy*yy + zz*zz);
    return 0.5 * log( (xyz + zz) / (xyz - zz) );}

  // Member functions that perform operations.
  void rescale3(double fac) {xx *= fac; yy *= fac; zz *= fac;}
  void rescale4(double fac) {xx *= fac; yy *= fac; zz *= fac; tt *= fac;}
  void flip3() {xx = -xx; yy = -yy; zz = -zz;}
  void flip4() {xx = -xx; yy = -yy; zz = -zz; tt = -tt;}
  void rot(double thetaIn, double phiIn);
  void rotaxis(double phiIn, double nx, double ny, double nz);
  void rotaxis(double phiIn, const Vec4& n);
  void bst(double betaX, double betaY, double betaZ);
  void bst(double betaX, double betaY, double betaZ, double gamma);
  void bst(const Vec4& pIn);
  void bst(const Vec4& pIn, double mIn);
  void bstback(const Vec4& pIn);
  void bstback(const Vec4& pIn, double mIn);
  void rotbst(const RotBstMatrix& M);

  // Operator overloading with member functions
  Vec4 operator-() {Vec4 tmp; tmp.xx = -xx; tmp.yy = -yy; tmp.zz = -zz;
    tmp.tt = -tt; return tmp;}
  Vec4& operator+=(const Vec4& v) {xx += v.xx; yy += v.yy; zz += v.zz;
    tt += v.tt; return *this;}
  Vec4& operator-=(const Vec4& v) {xx -= v.xx; yy -= v.yy; zz -= v.zz;
    tt -= v.tt; return *this;}
  Vec4& operator*=(double f) {xx *= f; yy *= f; zz *= f;
    tt *= f; return *this;}
  Vec4& operator/=(double f) {xx /= f; yy /= f; zz /= f;
    tt /= f; return *this;}

  // Operator overloading with friends
  friend Vec4 operator+(const Vec4& v1, const Vec4& v2);
  friend Vec4 operator-(const Vec4& v1, const Vec4& v2);
  friend Vec4 operator*(double f, const Vec4& v1);
  friend Vec4 operator*(const Vec4& v1, double f);
  friend Vec4 operator/(const Vec4& v1, double f);
  friend double operator*(const Vec4& v1, const Vec4& v2);

  // Invariant mass of a pair and its square.
  friend double m(const Vec4& v1, const Vec4& v2);
  friend double m2(const Vec4& v1, const Vec4& v2);

  // Scalar and cross product of 3-vector parts.
  friend double dot3(const Vec4& v1, const Vec4& v2);
  friend Vec4 cross3(const Vec4& v1, const Vec4& v2);

  // theta is polar angle between v1 and v2.
  friend double theta(const Vec4& v1, const Vec4& v2);
  friend double costheta(const Vec4& v1, const Vec4& v2);

  // phi is azimuthal angle between v1 and v2 around z axis.
  friend double phi(const Vec4& v1, const Vec4& v2);
  friend double cosphi(const Vec4& v1, const Vec4& v2);

  // phi is azimuthal angle between v1 and v2 around n axis.
  friend double phi(const Vec4& v1, const Vec4& v2, const Vec4& n);
  friend double cosphi(const Vec4& v1, const Vec4& v2, const Vec4& n);

  // R is distance in cylindrical (y/eta, phi) coordinates.
  friend double RRapPhi(const Vec4& v1, const Vec4& v2);
  friend double REtaPhi(const Vec4& v1, const Vec4& v2);

  // Print a four-vector.
  friend ostream& operator<<(ostream&, const Vec4& v) ;

private:

  // Constants: could only be changed in the code itself.
  static const double TINY;

  // The four-vector data members.
  double xx, yy, zz, tt;

};

//--------------------------------------------------------------------------

// Namespace function declarations; friends of Vec4 class.

// Implementation of operator overloading with friends.

inline Vec4 operator+(const Vec4& v1, const Vec4& v2)
  {Vec4 v = v1 ; return v += v2;}

inline Vec4 operator-(const Vec4& v1, const Vec4& v2)
  {Vec4 v = v1 ; return v -= v2;}

inline Vec4 operator*(double f, const Vec4& v1)
  {Vec4 v = v1; return v *= f;}

inline Vec4 operator*(const Vec4& v1, double f)
  {Vec4 v = v1; return v *= f;}

inline Vec4 operator/(const Vec4& v1, double f)
  {Vec4 v = v1; return v /= f;}

inline double operator*(const Vec4& v1, const Vec4& v2)
  {return v1.tt*v2.tt - v1.xx*v2.xx - v1.yy*v2.yy - v1.zz*v2.zz;}

// Invariant mass of a pair and its square.
double m(const Vec4& v1, const Vec4& v2);
double m2(const Vec4& v1, const Vec4& v2);

// Scalar and cross product of 3-vector parts.
double dot3(const Vec4& v1, const Vec4& v2);
Vec4 cross3(const Vec4& v1, const Vec4& v2);

// theta is polar angle between v1 and v2.
double theta(const Vec4& v1, const Vec4& v2);
double costheta(const Vec4& v1, const Vec4& v2);

// phi is azimuthal angle between v1 and v2 around z axis.
double phi(const Vec4& v1, const Vec4& v2);
double cosphi(const Vec4& v1, const Vec4& v2);

// phi is azimuthal angle between v1 and v2 around n axis.
double phi(const Vec4& v1, const Vec4& v2, const Vec4& n);
double cosphi(const Vec4& v1, const Vec4& v2, const Vec4& n);

// R is distance in cylindrical (y/eta, phi) coordinates.
double RRapPhi(const Vec4& v1, const Vec4& v2);
double REtaPhi(const Vec4& v1, const Vec4& v2);

// Print a four-vector.
ostream& operator<<(ostream&, const Vec4& v) ;

//==========================================================================

// RotBstMatrix class.
// This class implements 4 * 4 matrices that encode an arbitrary combination
// of rotations and boosts, that can be applied to Vec4 four-vectors.

class RotBstMatrix {

public:

  // Constructors.
  RotBstMatrix() {for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j)
    { M[i][j] = (i==j) ? 1. : 0.; } } }
  RotBstMatrix(const RotBstMatrix& Min) {
    for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) {
    M[i][j] = Min.M[i][j]; } } }
  RotBstMatrix& operator=(const RotBstMatrix& Min) {if (this != &Min) {
    for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) {
    M[i][j] = Min.M[i][j]; } } } return *this; }

  // Member functions.
  void rot(double = 0., double = 0.);
  void rot(const Vec4& p);
  void bst(double = 0., double = 0., double = 0.);
  void bst(const Vec4&);
  void bstback(const Vec4&);
  void bst(const Vec4&, const Vec4&);
  void toCMframe(const Vec4&, const Vec4&);
  void fromCMframe(const Vec4&, const Vec4&);
  void rotbst(const RotBstMatrix&);
  void invert();
  void reset();

  // Crude estimate deviation from unit matrix.
  double deviation() const;

  // Print a transformation matrix.
  friend ostream& operator<<(ostream&, const RotBstMatrix&) ;

  // Private members to be accessible from Vec4.
  friend class Vec4;

private:

  // Constants: could only be changed in the code itself.
  static const double TINY;

  // The rotation-and-boost matrix data members.
  double M[4][4];

};

//--------------------------------------------------------------------------

// Namespace function declaration; friend of RotBstMatrix class.

// Print a transformation matrix.
ostream& operator<<(ostream&, const RotBstMatrix&) ;

//==========================================================================

// Hist class.
// This class handles a single histogram at a time.

class Hist{

public:

  // Constructors, including copy constructors.
  Hist() {;}
  Hist(string titleIn, int nBinIn = 100, double xMinIn = 0.,
    double xMaxIn = 1.) {
    book(titleIn, nBinIn, xMinIn, xMaxIn);}
  Hist(const Hist& h)
    : title(h.title), nBin(h.nBin), nFill(h.nFill), xMin(h.xMin),
    xMax(h.xMax), dx(h.dx), under(h.under), inside(h.inside),
    over(h.over), res(h.res) { }
  Hist(string titleIn, const Hist& h)
    : title(titleIn), nBin(h.nBin), nFill(h.nFill), xMin(h.xMin),
    xMax(h.xMax), dx(h.dx), under(h.under), inside(h.inside),
    over(h.over), res(h.res) { }
  Hist& operator=(const Hist& h) { if(this != &h) {
    nBin = h.nBin; nFill = h.nFill; xMin = h.xMin; xMax = h.xMax;
    dx = h.dx;  under = h.under; inside = h.inside; over = h.over;
    res = h.res; } return *this; }

  // Book a histogram.
  void book(string titleIn = "  ", int nBinIn = 100, double xMinIn = 0.,
    double xMaxIn = 1.) ;

  // Set title of a histogram.
  void name(string titleIn = "  ") {title = titleIn; }

  // Reset bin contents.
  void null() ;

  // Fill bin with weight.
  void fill(double x, double w = 1.) ;

  // Print a histogram with overloaded << operator.
  friend ostream& operator<<(ostream& os, const Hist& h) ;

  // Print histogram contents as a table (e.g. for Gnuplot).
  void table(ostream& os = cout, bool printOverUnder = false,
    bool xMidBin = true) const ;
  void table(string fileName, bool printOverUnder = false,
    bool xMidBin = true) const { ofstream streamName(fileName.c_str());
    table(streamName, printOverUnder, xMidBin);}

  // Print a table out of two histograms with same x axis.
  friend void table(const Hist& h1, const Hist& h2, ostream& os,
    bool printOverUnder, bool xMidBin) ;
  friend void table(const Hist& h1, const Hist& h2, string fileName,
    bool printOverUnder, bool xMidBin) ;

  // Return content of specific bin: 0 gives underflow and nBin+1 overflow.
  double getBinContent(int iBin) const;

  // Return number of entries.
  int getEntries() const {return nFill; }

  // Check whether another histogram has same size and limits.
  bool sameSize(const Hist& h) const ;

  // Take logarithm (base 10 or e) of bin contents.
  void takeLog(bool tenLog = true) ;

  // Take square root of bin contents.
  void takeSqrt() ;

  // Operator overloading with member functions
  Hist& operator+=(const Hist& h) ;
  Hist& operator-=(const Hist& h) ;
  Hist& operator*=(const Hist& h) ;
  Hist& operator/=(const Hist& h) ;
  Hist& operator+=(double f) ;
  Hist& operator-=(double f) ;
  Hist& operator*=(double f) ;
  Hist& operator/=(double f) ;

  // Operator overloading with friends
  friend Hist operator+(double f, const Hist& h1);
  friend Hist operator+(const Hist& h1, double f);
  friend Hist operator+(const Hist& h1, const Hist& h2);
  friend Hist operator-(double f, const Hist& h1);
  friend Hist operator-(const Hist& h1, double f);
  friend Hist operator-(const Hist& h1, const Hist& h2);
  friend Hist operator*(double f, const Hist& h1);
  friend Hist operator*(const Hist& h1, double f);
  friend Hist operator*(const Hist& h1, const Hist& h2);
  friend Hist operator/(double f, const Hist& h1);
  friend Hist operator/(const Hist& h1, double f);
  friend Hist operator/(const Hist& h1, const Hist& h2);

private:

  // Constants: could only be changed in the code itself.
  static const int    NBINMAX, NCOLMAX, NLINES;
  static const double TOLERANCE, TINY, LARGE, SMALLFRAC, DYAC[];
  static const char   NUMBER[];

  // Properties and contents of a histogram.
  string title;
  int    nBin, nFill;
  double xMin, xMax, dx, under, inside, over;
  vector<double> res;

};

//--------------------------------------------------------------------------

// Namespace function declarations; friends of Hist class.

// Print a histogram with overloaded << operator.
ostream& operator<<(ostream& os, const Hist& h) ;

// Print a table out of two histograms with same x axis.
void table(const Hist& h1, const Hist& h2, ostream& os = cout,
  bool printOverUnder = false, bool xMidBin = true) ;
void table(const Hist& h1, const Hist& h2, string fileName,
  bool printOverUnder = false, bool xMidBin = true) ;

// Operator overloading with friends
Hist operator+(double f, const Hist& h1);
Hist operator+(const Hist& h1, double f);
Hist operator+(const Hist& h1, const Hist& h2);
Hist operator-(double f, const Hist& h1);
Hist operator-(const Hist& h1, double f);
Hist operator-(const Hist& h1, const Hist& h2);
Hist operator*(double f, const Hist& h1);
Hist operator*(const Hist& h1, double f);
Hist operator*(const Hist& h1, const Hist& h2);
Hist operator/(double f, const Hist& h1);
Hist operator/(const Hist& h1, double f);
Hist operator/(const Hist& h1, const Hist& h2);

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_Basics_H
