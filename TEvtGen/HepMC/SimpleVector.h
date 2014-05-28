//////////////////////////////////////////////////////////////////////////
// SimpleVector.h
//////////////////////////////////////////////////////////////////////////
#ifndef  HEPMC_SIMPLEVECTOR_H
#define  HEPMC_SIMPLEVECTOR_H

//////////////////////////////////////////////////////////////////////////
// garren@fnal.gov, July 2006
//
// This header provides a place to hold the doubles which are part of one of 
// three types of physics vectors:
//    momentum 4 vector 
//    position or displacement 4 vector
//    position or displacement 3 vector
//
// For compatibility with existing code, 
// the basic expected geometrical access methods are povided 
// Also, both FourVector and ThreeVector have a templated constructor that will 
// take another vector (HepLorentzVector, GenVector, ...)
//    --> this vector must have the following methods: x(), y(), z()
//    -->  FourVector also requires the t() method
//
//////////////////////////////////////////////////////////////////////////


#include "HepMC/enable_if.h"
#include "HepMC/is_arithmetic.h"


namespace HepMC {

//! FourVector is a simple representation of a physics 4 vector

///
/// \class  FourVector
/// For compatibility with existing code, 
/// the basic expected geometrical access methods are povided.
/// Also, there is a templated constructor that will 
/// take another vector (HepLorentzVector, GenVector, ...)
/// which must have the following methods: x(), y(), z(), t().
///
class FourVector {

public:

  /// constructor requiring at least x, y, and z
  FourVector( double xin, double yin, double zin, double tin=0) 
  : m_x(xin), m_y(yin), m_z(zin), m_t(tin) {}

  /// constructor requiring only t 
  FourVector(double tin)
  : m_x(0), m_y(0), m_z(0), m_t(tin) {}

  FourVector() 
  : m_x(0), m_y(0), m_z(0), m_t(0) {}

  /// templated constructor
  /// this is used ONLY if T is not arithmetic
  template <class T >
  FourVector( const T& v,
         typename detail::disable_if< detail::is_arithmetic<T>::value, void >::type * = 0 )
  : m_x(v.x()), m_y(v.y()), m_z(v.z()), m_t(v.t()) {}

  /// copy constructor
  FourVector(const FourVector & v)
  : m_x(v.x()), m_y(v.y()), m_z(v.z()), m_t(v.t()) {}

  void swap( FourVector & other );  //!< swap

  double px() const { return m_x; }  //!< return px
  double py() const { return m_y; }  //!< return py
  double pz() const { return m_z; }  //!< return pz
  double e()  const { return m_t; }  //!< return E

  double x() const { return m_x; }  //!< return x
  double y() const { return m_y; }  //!< return y
  double z() const { return m_z; }  //!< return z
  double t() const { return m_t; }  //!< return t

  double m2() const;  //!< Invariant mass squared.
  double m() const;   //!< Invariant mass. If m2() is negative then -sqrt(-m2()) is returned.

  double perp2() const;  //!< Transverse component of the spatial vector squared.
  double perp() const;   //!< Transverse component of the spatial vector (R in cylindrical system).

  // Get spatial vector components in spherical coordinate system.
  double theta() const;  //!< The polar angle.
  double phi() const;  //!< The azimuth angle.
  double rho() const;  //!< spatial vector component magnitude

  FourVector & operator = (const FourVector &); //!< make a copy

  bool operator == (const FourVector &) const; //!< equality
  bool operator != (const FourVector &) const; //!< inequality

  double pseudoRapidity() const;  //!< Returns the pseudo-rapidity, i.e. -ln(tan(theta/2))
  double eta() const;             //!< Pseudorapidity (of the space part)

  /// set x, y, z, and t
  void set        (double x, double y, double z, double  t);

  void setX(double xin) { m_x=xin; }  //!< set x
  void setY(double yin) { m_y=yin; }  //!< set y
  void setZ(double zin) { m_z=zin; }  //!< set z
  void setT(double tin) { m_t=tin; }  //!< set t

  void setPx(double xin) { m_x=xin; }  //!< set px
  void setPy(double yin) { m_y=yin; }  //!< set py
  void setPz(double zin) { m_z=zin; }  //!< set pz
  void setE(double tin)  { m_t=tin; }  //!< set E

private:

  double m_x;
  double m_y;
  double m_z;
  double m_t;
  
};

//! ThreeVector is a simple representation of a position or displacement 3 vector

///
/// \class  ThreeVector
/// For compatibility with existing code, 
/// the basic expected geometrical access methods are povided.
/// Also, there is a templated constructor that will 
/// take another vector (HepLorentzVector, GenVector, ...)
/// which must have the following methods: x(), y(), z().
///
class ThreeVector {

public:

  /// construct using x, y, and z (only x is required)
  ThreeVector( double xin, double yin =0, double zin =0 ) 
  : m_x(xin), m_y(yin), m_z(zin) {}

  ThreeVector( ) 
  : m_x(0), m_y(0), m_z(0) {}
  
  /// templated constructor
  /// this is used ONLY if T is not arithmetic
  template <class T >
  ThreeVector( const T& v,
         typename detail::disable_if< detail::is_arithmetic<T>::value, void >::type * = 0 )
  : m_x(v.x()), m_y(v.y()), m_z(v.z()) {}

  /// copy constructor
  ThreeVector(const ThreeVector & v)
  : m_x(v.x()), m_y(v.y()), m_z(v.z()) {}

  void swap( ThreeVector & other );  //!< swap

  double x() const { return m_x; }  //!< return x
  double y() const { return m_y; }  //!< return y
  double z() const { return m_z; }  //!< return z

  void setX(double xin) { m_x=xin; }  //!< set x
  void setY(double yin) { m_y=yin; }  //!< set y
  void setZ(double zin) { m_z=zin; }  //!< set z
  void set( double x, double y, double z);   //!< set x, y, and z

  double phi()   const;  //!< The azimuth angle.
  double theta() const;  //!< The polar angle.
  double r()     const;  //!< The magnitude

  void setPhi(double);  //!< Set phi keeping magnitude and theta constant (BaBar).
  void setTheta(double);  //!< Set theta keeping magnitude and phi constant (BaBar).

  double perp2() const;  //!< The transverse component squared (rho^2 in cylindrical coordinate system).
  double perp() const;  //!< The transverse component (rho in cylindrical coordinate system).

  ThreeVector & operator = (const ThreeVector &); //!< make a copy

  bool operator == (const ThreeVector &) const; //!< equality
  bool operator != (const ThreeVector &) const; //!< inequality

private:

  double m_x;
  double m_y;
  double m_z;

};  


} // HepMC

#include "HepMC/SimpleVector.icc"

#endif  // HEPMC_SIMPLEVECTOR_H

