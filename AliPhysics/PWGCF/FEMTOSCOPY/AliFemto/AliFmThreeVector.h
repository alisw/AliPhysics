///
/// \file AliFmThreeVector.h
/// \authors Brian Lasiuk, Thomas Ullrich
/// \date April 1998
///
/// \class AliFmThreeVector
/// \brief A templated class for three-vector objects, with common geometric
///        vector methods. Common usage is via the typedef'd AliFemtoThreeVector
///        which fills the template with the float type.
///


#ifndef ALIFMTHREEVECTOR_H
#define ALIFMTHREEVECTOR_H

#ifdef __ROOT__
#include "Rtypes.h"
#endif
#ifndef __CINT__
#include <iostream>
#include <fstream>
#include <math.h>
#ifdef GNU_GCC
#    include <stddef.h>
#endif
#if defined (__SUNPRO_CC) && __SUNPRO_CC < 0x500
#    include <stdcomp.h>
#endif
#ifndef ST_NO_EXCEPTIONS
#    include <stdexcept>
#    if !defined(ST_NO_NAMESPACES)
using std::out_of_range;
#    endif
#endif
#endif // __CINT__

#ifdef WIN32
#include "gcc2vs.h"
#endif

class TRootIOCtor;//nic nie rozumiem
using namespace std;


template<class T>
class AliFmThreeVector {
public:
  AliFmThreeVector(T=0, T=0, T=0);
  //                     ROOT_VERSION(5,03,01)
#if ROOT_VERSION_CODE >= 328449
  AliFmThreeVector(TRootIOCtor *) : mX1(0), mX2(0), mX3(0) {}
#endif
  virtual ~AliFmThreeVector();

#ifndef ST_NO_MEMBER_TEMPLATES
  template<class X> AliFmThreeVector(const AliFmThreeVector<X> &);
  template<class X> AliFmThreeVector(const X *);
  template<class X> AliFmThreeVector<T> &operator=(const AliFmThreeVector<X> &);
  // AliFmThreeVector(const AliFmThreeVector<T>&);                use default
  // AliFmThreeVector<T>& operator=(const AliFmThreeVector<T>&);  use default
#else
  AliFmThreeVector(const AliFmThreeVector<float> &);
  AliFmThreeVector(const AliFmThreeVector<double> &);

  AliFmThreeVector(const float *);
  AliFmThreeVector(const double *);

  AliFmThreeVector<T> &operator=(const AliFmThreeVector<float> &);
  AliFmThreeVector<T> &operator=(const AliFmThreeVector<double> &);
#endif

  void SetX(T);
  void SetY(T);
  void SetZ(T);

  void SetPhi(T);       ///< Calculates and sets the X & Y components using the current vector magitude and the provided angle (in radians)
  void SetTheta(T);     ///< Calculates and sets the caresian components using the current vector magitude and the provided angle (in radians)
  void SetMag(T);       ///< Alias of SetMagnitude
  void SetMagnitude(T); ///< Calculates and sets the caresian components, keeping the same direction but scaled to have a magnitude set to the method parameter

  T x() const;                  ///< First cartesian component
  T y() const;                  ///< Second cartesian component
  T z() const;                  ///< Third cartesian component
  T Theta() const;              ///< Polar angle calculated with respect to the Z-axis
  T CosTheta() const;           ///< Cosine of the polar angle
  T Phi() const;                ///< Azimuthal angle calculated from X & Y components
  T Perp() const;               ///< Return the magnitude of the perpendicular component of this vector (X & Y components)
  T Perp2() const;              ///< Returns the square of the perpendicular component
  T Magnitude() const;          ///< Return the magnitude of this vector
  T Mag() const;                ///< Alias of Magnitude
  T Mag2() const;               ///< Returns the square of the magnitude of this vector
  T PseudoRapidity() const;     ///< Calculates the PseudoRapidity from the vector's polar angle
  T operator()(size_t) const;
  T operator[](size_t) const;

  T &operator()(size_t);
  T &operator[](size_t);

  T MassHypothesis(T mass) const; ///< Return the equivalent energy of a relativistic particle with provided mass moving with this vector as 3-momentum

  AliFmThreeVector<T> unit() const;         ///< Return a vector pointing in same direction but scaled to unit length
  AliFmThreeVector<T> orthogonal() const;   ///< Return a vector orthongonal to this one

  void RotateX(T);  ///< Rotate this vector around X-axis (angle in radians)
  void RotateY(T);  ///< Rotate this vector around Y-axis (angle in radians)
  void RotateZ(T);  ///< Rotate this vector around Z-axis (angle in radians)

  AliFmThreeVector<T> operator-();  ///< Return a vector with each component flipped
  AliFmThreeVector<T> operator+();
  AliFmThreeVector<T> &operator*= (double); ///< Scale this vector by a factor
  AliFmThreeVector<T> &operator/= (double); ///< Scale this vector by an inverted factor
  AliFmThreeVector<T> PseudoProduct(double, double, double) const;  ///< Returns a vector with each X-Y-Z component scaled by its respective parameter

#ifndef ST_NO_MEMBER_TEMPLATES
  template<class X> T Angle(const AliFmThreeVector<X> &) const;
  template<class X> AliFmThreeVector<T> Cross(const AliFmThreeVector<X> &) const;
  template<class X> T Dot(const AliFmThreeVector<X> &) const;
  template<class X> AliFmThreeVector<T> PseudoProduct(const AliFmThreeVector<X> &) const;

  template<class X> bool operator == (const AliFmThreeVector<X> &v) const;
  template<class X> bool operator != (const AliFmThreeVector<X> &v) const;

  template<class X> AliFmThreeVector<T> &operator+= (const AliFmThreeVector<X> &);
  template<class X> AliFmThreeVector<T> &operator-= (const AliFmThreeVector<X> &);
#else
  T Angle(const AliFmThreeVector<float> &) const;
  AliFmThreeVector<T> Cross(const AliFmThreeVector<float> &) const;
  T Dot(const AliFmThreeVector<float> &) const;
  AliFmThreeVector<T> PseudoProduct(const AliFmThreeVector<float> &) const;

  T Angle(const AliFmThreeVector<double> &) const;
  T Dot(const AliFmThreeVector<double> &) const;
  AliFmThreeVector<T> Cross(const AliFmThreeVector<double> &) const;
  AliFmThreeVector<T> PseudoProduct(const AliFmThreeVector<double> &) const;

  bool operator == (const AliFmThreeVector<float> &v) const;
  bool operator != (const AliFmThreeVector<float> &v) const;
  AliFmThreeVector<T> &operator+= (const AliFmThreeVector<float> &);
  AliFmThreeVector<T> &operator-= (const AliFmThreeVector<float> &);

  bool operator == (const AliFmThreeVector<double> &v) const;
  bool operator != (const AliFmThreeVector<double> &v) const;
  AliFmThreeVector<T> &operator+= (const AliFmThreeVector<double> &);
  AliFmThreeVector<T> &operator-= (const AliFmThreeVector<double> &);
#endif
  int Valid(double world = 1.e+5) const;
  int Bad(double world = 1.e+5) const;
protected:
  T    mX1, mX2, mX3;  // Vector components
#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFmThreeVector, 3);
  /// \endcond
#endif /* __ROOT__ */
};

#ifndef __CINT__
//
//        Implementation of member functions
//
template<class T>
inline AliFmThreeVector<T>::AliFmThreeVector(T ax, T ay, T az): mX1(ax), mX2(ay), mX3(az)
{
  /* nop */
}
template<class T>
inline AliFmThreeVector<T>::~AliFmThreeVector()
{
  /* nop */
}

template<class T>
inline void AliFmThreeVector<T>::SetX(T ax)
{
  mX1 = ax;
}

template<class T>
inline void AliFmThreeVector<T>::SetY(T ay)
{
  mX2 = ay;
}

template<class T>
inline void AliFmThreeVector<T>::SetZ(T az)
{
  mX3 = az;
}

template<class T>
void AliFmThreeVector<T>::SetPhi(T aAngle)
{
  double  r = Magnitude();
  double th = Theta();

  mX1 = r * sin(th) * cos(aAngle);
  mX2 = r * sin(th) * sin(aAngle);
}

template <class T>
void AliFmThreeVector<T>::SetTheta(T aAngle)
{
  double r  = Magnitude();
  double ph = Phi();

  mX1 = r * sin(aAngle) * cos(ph);
  mX2 = r * sin(aAngle) * sin(ph);
  mX3 = r * cos(aAngle);
}

template <class T>
void AliFmThreeVector<T>::SetMagnitude(T r)
{
  double th = Theta();
  double ph = Phi();

  mX1 = r * sin(th) * cos(ph);
  mX2 = r * sin(th) * sin(ph);
  mX3 = r * cos(th);
}

template <class T>
void AliFmThreeVector<T>::SetMag(T amag)
{
  SetMagnitude(amag);
}

template<class T>
inline T AliFmThreeVector<T>::x() const
{
  return mX1;
}

template<class T>
inline T AliFmThreeVector<T>::y() const
{
  return mX2;
}

template<class T>
inline T AliFmThreeVector<T>::z() const
{
  return mX3;
}

template<class T>
inline T AliFmThreeVector<T>::Theta() const
{
  return acos(CosTheta());
}

template<class T>
inline T AliFmThreeVector<T>::CosTheta() const
{
  return mX3 / (Mag() + 1e-20);
}

template<class T>
inline T AliFmThreeVector<T>::Phi() const
{
  return atan2(mX2, mX1);
}

template<class T>
inline T AliFmThreeVector<T>::PseudoRapidity() const
{
  //
  // change code to more optimal:
  // double m = Mag();
  // return 0.5*::log( (m+z())/(m-z()) );
  double tmp = tan(Theta() / 2.0);
  if (tmp <= 0.) return 1e20;
  return -::log(tmp);
}

template<class T>
inline AliFmThreeVector<T> AliFmThreeVector<T>::unit() const
{
  double tmp = Mag();
  if (tmp <= 0.) tmp = 1e-20;
  return *this / tmp;
}

template <class T>
T AliFmThreeVector<T>::MassHypothesis(T mass) const
{
  return ::sqrt((*this) * (*this) + mass * mass);
}

template <class T>
AliFmThreeVector<T> AliFmThreeVector<T>::orthogonal() const
{
  // Direct copy from CLHEP--it is probably better to
  // use your own dot/cross product code...
  double ax = (mX1 < 0.0) ? -mX1 : mX1;
  double ay = (mX2 < 0.0) ? -mX2 : mX2;
  double az = (mX3 < 0.0) ? -mX3 : mX3;

  if (ax < ay)
    return ax < az ? AliFmThreeVector<T>(0, mX3, -mX2) :  AliFmThreeVector<T>(mX2, -mX1, 0);
  else
    return  mX2 < mX3 ? AliFmThreeVector<T>(-mX3, 0, mX1) :  AliFmThreeVector<T>(mX2, -mX1, 0);
}

template <class T>
void AliFmThreeVector<T>::RotateX(T aAngle)
{
  // may in the future make use of the AliFmRotation class!
  double yPrime = cos(aAngle) * mX2 - sin(aAngle) * mX3;
  double zPrime = sin(aAngle) * mX2 + cos(aAngle) * mX3;

  mX2 = yPrime;
  mX3 = zPrime;
}

template <class T>
void AliFmThreeVector<T>::RotateY(T aAngle)
{
  // may in the future make use of the AliFmRotation class!
  double zPrime = cos(aAngle) * mX3 - sin(aAngle) * mX1;
  double xPrime = sin(aAngle) * mX3 + cos(aAngle) * mX1;

  mX1 = xPrime;
  mX3 = zPrime;
}

template <class T>
void AliFmThreeVector<T>::RotateZ(T aAngle)
{
  // may in the future make use of the AliFmRotation class!
  double xPrime = cos(aAngle) * mX1 - sin(aAngle) * mX2;
  double yPrime = sin(aAngle) * mX1 + cos(aAngle) * mX2;

  mX1 = xPrime;
  mX2 = yPrime;
}

template<class T>
inline T AliFmThreeVector<T>::Perp() const
{
  return ::sqrt(mX1 * mX1 + mX2 * mX2);
}

template<class T>
inline T AliFmThreeVector<T>::Perp2() const
{
  return mX1 * mX1 + mX2 * mX2;
}

template<class T>
inline T AliFmThreeVector<T>::Magnitude() const
{
  return Mag();
}

template<class T>
inline T AliFmThreeVector<T>::Mag() const
{
  return ::sqrt(mX1 * mX1 + mX2 * mX2 + mX3 * mX3);
}

template<class T>
inline T AliFmThreeVector<T>::Mag2() const
{
  return mX1 * mX1 + mX2 * mX2 + mX3 * mX3;
}

template<class T>
inline T AliFmThreeVector<T>::operator()(size_t i) const
{
  if (i <= 2)  return (&mX1)[i];
#ifndef ST_NO_EXCEPTIONS
  throw out_of_range("AliFmThreeVector<T>::operator(): bad index");
#else
  cerr << "AliFmThreeVector<T>::operator(): bad index" << endl;
#endif
  return 0;
}

template<class T>
inline T &AliFmThreeVector<T>::operator()(size_t i)
{
  if (i <= 2)  return (&mX1)[i];
#ifndef ST_NO_EXCEPTIONS
  throw out_of_range("AliFmThreeVector<T>::operator(): bad index");
#else
  cerr << "AliFmThreeVector<T>::operator(): bad index" << endl;
  return mX1;
#endif
}

template<class T>
inline T AliFmThreeVector<T>::operator[](size_t i) const
{
  if (i <= 2)  return (&mX1)[i];
#ifndef ST_NO_EXCEPTIONS
  throw out_of_range("AliFmThreeVector<T>::operator[]: bad index");
#else
  cerr << "AliFmThreeVector<T>::operator[]: bad index" << endl;
#endif
  return 0;
}

template<class T>
inline T &AliFmThreeVector<T>::operator[](size_t i)
{
  if (i <= 2)  return (&mX1)[i];
#ifndef ST_NO_EXCEPTIONS
  throw out_of_range("AliFmThreeVector<T>::operator[]: bad index");
#else
  cerr << "AliFmThreeVector<T>::operator[]: bad index" << endl;
  return mX1;
#endif
}

template<class T>
inline AliFmThreeVector<T> &AliFmThreeVector<T>::operator*= (double c)
{
  mX1 *= c;
  mX2 *= c;
  mX3 *= c;
  return *this;
}

template<class T>
inline AliFmThreeVector<T> &AliFmThreeVector<T>::operator/= (double c)
{
  mX1 /= c;
  mX2 /= c;
  mX3 /= c;
  return *this;
}

template<class T>
inline AliFmThreeVector<T>
AliFmThreeVector<T>::PseudoProduct(double ax, double ay, double az) const
{
  return AliFmThreeVector<T>(mX1 * ax, mX2 * ay, mX3 * az);
}

template<class T>
AliFmThreeVector<T> AliFmThreeVector<T>::operator- ()
{
  return AliFmThreeVector<T>(-mX1, -mX2, -mX3);
}

template<class T>
AliFmThreeVector<T> AliFmThreeVector<T>::operator+ ()
{
  return *this;
}

#ifndef ST_NO_MEMBER_TEMPLATES
#ifndef WIN32

template<class T>
template<class X>
inline AliFmThreeVector<T>::AliFmThreeVector(const AliFmThreeVector<X> &v)
  : mX1(v.x()), mX2(v.y()), mX3(v.z())
{
  /* nop */
}

template<class T>
template<class X>
inline AliFmThreeVector<T>::AliFmThreeVector(const X *a)
{
  mX1 = a[0];
  mX2 = a[1];
  mX3 = a[2];
}

template<class T>
template<class X>
inline AliFmThreeVector<T> &
AliFmThreeVector<T>::operator=(const AliFmThreeVector<X> &v)
{
  mX1 = v.x();
  mX2 = v.y();
  mX3 = v.z();
  return *this;
}

template<class T>
template<class X>
inline bool AliFmThreeVector<T>::operator== (const AliFmThreeVector<X> &v) const
{
  return mX1 == v.x() && mX2 == v.y() && mX3 == v.z();
}

template<class T>
template<class X>
inline bool AliFmThreeVector<T>::operator!= (const AliFmThreeVector<X> &v) const
{
  return !(*this == v);
}

template<class T>
template<class X>
inline AliFmThreeVector<T> &
AliFmThreeVector<T>::operator+= (const AliFmThreeVector<X> &v)
{
  mX1 += v.x();
  mX2 += v.y();
  mX3 += v.z();
  return *this;
}

template<class T>
template<class X>
inline AliFmThreeVector<T> &
AliFmThreeVector<T>::operator-= (const AliFmThreeVector<X> &v)
{
  mX1 -= v.x();
  mX2 -= v.y();
  mX3 -= v.z();
  return *this;
}

template<class T>
template<class X>
inline T AliFmThreeVector<T>::Dot(const AliFmThreeVector<X> &v) const
{
  return mX1 * v.x() + mX2 * v.y() + mX3 * v.z();
}

template<class T>
template<class X>
inline AliFmThreeVector<T>
AliFmThreeVector<T>::Cross(const AliFmThreeVector<X> &v) const
{
  return AliFmThreeVector<T>(mX2 * v.z() - mX3 * v.y(),
                             mX3 * v.x() - mX1 * v.z(),
                             mX1 * v.y() - mX2 * v.x());
}

template<class T>
template<class X>
inline T AliFmThreeVector<T>::Angle(const AliFmThreeVector<X> &vec) const
{
  double norm = this->Mag2() * vec.Mag2();

  return norm > 0 ? acos(this->Dot(vec) / (::sqrt(norm))) : 0;
}

template<class T>
template<class X>
inline AliFmThreeVector<T>
AliFmThreeVector<T>::PseudoProduct(const AliFmThreeVector<X> &v) const
{
  return this->PseudoProduct(v.x(), v.y(), v.z());
}

#endif
#else

template<class T>
inline AliFmThreeVector<T>::AliFmThreeVector(const AliFmThreeVector<float> &v)
  : mX1(v.x()), mX2(v.y()), mX3(v.z())
{
  /* nop */
}

template<class T>
inline AliFmThreeVector<T>::AliFmThreeVector(const AliFmThreeVector<double> &v)
  : mX1(v.x()), mX2(v.y()), mX3(v.z())
{
  /* nop */
}

template<class T>
inline AliFmThreeVector<T>::AliFmThreeVector(const float *a)
{
  mX1 = a[0];
  mX2 = a[1];
  mX3 = a[2];
}

template<class T>
inline AliFmThreeVector<T>::AliFmThreeVector(const double *a)
{
  mX1 = a[0];
  mX2 = a[1];
  mX3 = a[2];
}

template<class T>
inline AliFmThreeVector<T> &
AliFmThreeVector<T>::operator=(const AliFmThreeVector<float> &v)
{
  mX1 = v.x();
  mX2 = v.y();
  mX3 = v.z();
  return *this;
}

template<class T>
inline AliFmThreeVector<T> &
AliFmThreeVector<T>::operator=(const AliFmThreeVector<double> &v)
{
  mX1 = v.x();
  mX2 = v.y();
  mX3 = v.z();
  return *this;
}

template<class T>
inline bool
AliFmThreeVector<T>::operator== (const AliFmThreeVector<float> &v) const
{
  return mX1 == v.x() && mX2 == v.y() && mX3 == v.z();
}

template<class T>
inline bool
AliFmThreeVector<T>::operator== (const AliFmThreeVector<double> &v) const
{
  return mX1 == v.x() && mX2 == v.y() && mX3 == v.z();
}

template<class T>
inline bool
AliFmThreeVector<T>::operator!= (const AliFmThreeVector<float> &v) const
{
  return !(*this == v);
}

template<class T>
inline bool
AliFmThreeVector<T>::operator!= (const AliFmThreeVector<double> &v) const
{
  return !(*this == v);
}

template<class T>
inline AliFmThreeVector<T> &
AliFmThreeVector<T>::operator+= (const AliFmThreeVector<float> &v)
{
  mX1 += v.x();
  mX2 += v.y();
  mX3 += v.z();
  return *this;
}

template<class T>
inline AliFmThreeVector<T> &
AliFmThreeVector<T>::operator+= (const AliFmThreeVector<double> &v)
{
  mX1 += v.x();
  mX2 += v.y();
  mX3 += v.z();
  return *this;
}

template<class T>
inline AliFmThreeVector<T> &
AliFmThreeVector<T>::operator-= (const AliFmThreeVector<float> &v)
{
  mX1 -= v.x();
  mX2 -= v.y();
  mX3 -= v.z();
  return *this;
}

template<class T>
inline AliFmThreeVector<T> &
AliFmThreeVector<T>::operator-= (const AliFmThreeVector<double> &v)
{
  mX1 -= v.x();
  mX2 -= v.y();
  mX3 -= v.z();
  return *this;
}

template<class T>
inline T AliFmThreeVector<T>::Dot(const AliFmThreeVector<float> &v) const
{
  return mX1 * v.x() + mX2 * v.y() + mX3 * v.z();
}

template<class T>
inline T AliFmThreeVector<T>::Dot(const AliFmThreeVector<double> &v) const
{
  return mX1 * v.x() + mX2 * v.y() + mX3 * v.z();
}

template<class T>
inline AliFmThreeVector<T>
AliFmThreeVector<T>::Cross(const AliFmThreeVector<float> &v) const
{
  return AliFmThreeVector<T>(mX2 * v.z() - mX3 * v.y(),
                             mX3 * v.x() - mX1 * v.z(),
                             mX1 * v.y() - mX2 * v.x());
}

template<class T>
inline AliFmThreeVector<T>
AliFmThreeVector<T>::Cross(const AliFmThreeVector<double> &v) const
{
  return AliFmThreeVector<T>(mX2 * v.z() - mX3 * v.y(),
                             mX3 * v.x() - mX1 * v.z(),
                             mX1 * v.y() - mX2 * v.x());
}

template<class T>
inline T AliFmThreeVector<T>::Angle(const AliFmThreeVector<float> &v) const
{
  double tmp = Mag() * v.Mag();
  if (tmp <= 0) tmp = 1e-20;
  return acos(this->Dot(v) / tmp);
}

template<class T>
inline T AliFmThreeVector<T>::Angle(const AliFmThreeVector<double> &v) const
{
  double tmp = Mag() * v.Mag();
  if (tmp <= 0) tmp = 1e-20;
  return acos(this->Dot(v) / tmp);
}

template<class T>
inline AliFmThreeVector<T>
AliFmThreeVector<T>::PseudoProduct(const AliFmThreeVector<float> &v) const
{
  return this->PseudoProduct(v.x(), v.y(), v.z());
}

template<class T>
inline AliFmThreeVector<T>
AliFmThreeVector<T>::PseudoProduct(const AliFmThreeVector<double> &v) const
{
  return this->PseudoProduct(v.x(), v.y(), v.z());
}
#endif  // ST_NO_MEMBER_TEMPLATES
template<class T>
inline int
AliFmThreeVector<T>::Valid(double world) const
{
  return !Bad(world);
}

template<class T>
inline int
AliFmThreeVector<T>::Bad(double world) const
{
  for (int i = 0; i < 3; i++) {
    if (!isfinite((&mX1)[i])) return 10 + i;
    if (fabs((&mX1)[i]) > world) return 20 + i;
  }
  return 0;
}
#endif /*! __CINT__ */
#ifdef __CINT__
template<> float abs(const AliFmThreeVector<float> &v);
template<> double abs(const AliFmThreeVector<double> &v);
template<> AliFmThreeVector<double> cross_product(const AliFmThreeVector<double> &v1, const AliFmThreeVector<double> &v2);
template<> AliFmThreeVector<float>  cross_product(const AliFmThreeVector<float>  &v1, const AliFmThreeVector<float> &v2);
template<> AliFmThreeVector<double> cross_product(const AliFmThreeVector<float>  &v1, const AliFmThreeVector<double> &v2);
template<> AliFmThreeVector<double> cross_product(const AliFmThreeVector<double> &v1, const AliFmThreeVector<float> &v2);
template<> AliFmThreeVector<double> operator+ (const AliFmThreeVector<double> &v1, const AliFmThreeVector<double> &v2);
template<> AliFmThreeVector<float>  operator+ (const AliFmThreeVector<float>  &v1, const AliFmThreeVector<float> &v2);
template<> AliFmThreeVector<double> operator+ (const AliFmThreeVector<double> &v1, const AliFmThreeVector<float> &v2);
template<> AliFmThreeVector<double> operator+ (const AliFmThreeVector<float>  &v1, const AliFmThreeVector<double> &v2);
template<> AliFmThreeVector<double> operator- (const AliFmThreeVector<double> &v1, const AliFmThreeVector<double> &v2);
template<> AliFmThreeVector<float>  operator- (const AliFmThreeVector<float>  &v1, const AliFmThreeVector<float> &v2);
template<> AliFmThreeVector<double> operator- (const AliFmThreeVector<double> &v1, const AliFmThreeVector<float> &v2);
template<> AliFmThreeVector<double> operator- (const AliFmThreeVector<float>  &v1, const AliFmThreeVector<double> &v2);
template<> AliFmThreeVector<double> operator* (const AliFmThreeVector<double> &v1, const AliFmThreeVector<double> &v2);
template<> AliFmThreeVector<float>  operator* (const AliFmThreeVector<float>  &v1, const AliFmThreeVector<float> &v2);
template<> AliFmThreeVector<double> operator* (const AliFmThreeVector<double> &v1, const AliFmThreeVector<float> &v2);
template<> AliFmThreeVector<double> operator* (const AliFmThreeVector<float>  &v1, const AliFmThreeVector<double> &v2);
template<> AliFmThreeVector<double> operator* (const double                 v1, const AliFmThreeVector<float> &v2);
template<> AliFmThreeVector<double> operator* (const AliFmThreeVector<float>  &v1, const double v2);
template<> AliFmThreeVector<double> operator* (const double                 v1, const AliFmThreeVector<double> &v2);
template<> AliFmThreeVector<double> operator* (const AliFmThreeVector<double> &v1, const double v2);
template<> AliFmThreeVector<double> operator/ (const AliFmThreeVector<double> &v1, const AliFmThreeVector<double> &v2);
template<> AliFmThreeVector<float>  operator/ (const AliFmThreeVector<float>  &v1, const AliFmThreeVector<float> &v2);
template<> AliFmThreeVector<double> operator/ (const AliFmThreeVector<double> &v1, const AliFmThreeVector<float> &v2);
template<> AliFmThreeVector<double> operator/ (const AliFmThreeVector<float>  &v1, const AliFmThreeVector<double> &v2);
template<> AliFmThreeVector<double> operator/ (const                double  v1, const AliFmThreeVector<double> &v2);
template<> AliFmThreeVector<float>  operator/ (const                double  v1, const AliFmThreeVector<float> &v2);
template<> AliFmThreeVector<double> operator/ (const AliFmThreeVector<double> &v1, const double v2);
template<> AliFmThreeVector<double> operator/ (const AliFmThreeVector<float>  &v1, const double v2);
template<> istream  &operator>>(istream &is, const AliFmThreeVector<double> &v);
template<> istream  &operator>>(istream &is, const AliFmThreeVector<float> &v);
template<> ostream  &operator<<(ostream &os, const AliFmThreeVector<double> &v);
template<> ostream  &operator<<(ostream &os, const AliFmThreeVector<float> &v);
#else

//
//        Non-member functions
//

template<class T>
inline T abs(const AliFmThreeVector<T> &v)
{
  return v.Mag();
}

template<class T, class X>
inline AliFmThreeVector<T>
cross_product(const AliFmThreeVector<T> &v1, const AliFmThreeVector<X> &v2)
{
  return v1.Cross(v2);
}


//
//        Non-member operators
//
template<class T, class X>
inline AliFmThreeVector<T>
operator+ (const AliFmThreeVector<T> &v1, const AliFmThreeVector<X> &v2)
{
  return AliFmThreeVector<T>(v1) += v2;
}

template<class T, class X>
inline AliFmThreeVector<T>
operator- (const AliFmThreeVector<T> &v1, const AliFmThreeVector<X> &v2)
{
  return AliFmThreeVector<T>(v1) -= v2;
}

template<class T, class X>
inline T operator* (const AliFmThreeVector<T> &v1, const AliFmThreeVector<X> &v2)
{
  return AliFmThreeVector<T>(v1).Dot(v2);
}

template<class T>
inline AliFmThreeVector<T> operator* (const AliFmThreeVector<T> &v, double c)
{
  return AliFmThreeVector<T>(v) *= c;
}

template<class T>
inline AliFmThreeVector<T> operator* (double c, const AliFmThreeVector<T> &v)
{
  return AliFmThreeVector<T>(v) *= c;
}

template<class T, class X>
inline AliFmThreeVector<T> operator/ (const AliFmThreeVector<T> &v, X c)
{
  return AliFmThreeVector<T>(v) /= c;
}

template<class T>
ostream  &operator<<(ostream &os, const AliFmThreeVector<T> &v)
{
  return os << v.x() << '\t' << v.y() << '\t' << v.z();
}

template<class T>
istream  &operator>>(istream &is, AliFmThreeVector<T> &v)
{
  T  x, y, z;
  is >> x >> y >> z;
  v.SetX(x);
  v.SetY(y);
  v.SetZ(z);
  return is;
}


#endif /* ! __CINT__ */

#endif /* ALIFMTHREEVECTOR_H */
