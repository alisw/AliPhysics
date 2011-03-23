/***************************************************************************
 *
 * $Id$
 *
 * Author: Brian Lasiuk, Thomas Ullrich, April 1998
 ***************************************************************************
 *
 * Description:
 *
 * Remarks:   Since not all compilers support member templates
 *            we have to specialize the templated member on these
 *            platforms. If member templates are not supported the
 *            ST_NO_MEMBER_TEMPLATES flag has to be set. tu.
 *
 *            In the near future when all compilers can handle member
 *            templates this class should be cleaned up. A lot of
 *            redundant code can be removed as soon as the compilers
 *            are up-to-date. tu
 *
 ***************************************************************************
 *
 * $Log$
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.11  2005/09/22 20:09:20  fisyak
 * Make AliFemtoLorentzVector persistent
 *
 * Revision 1.10  2005/07/06 18:49:56  fisyak
 * Replace AliFemtoHelixD, AliFemtoLorentzVectorD,AliFemtoLorentzVectorF,AliFemtoMatrixD,AliFemtoMatrixF,AliFemtoPhysicalHelixD,AliFemtoThreeVectorD,AliFemtoThreeVectorF by templated version
 *
 * Revision 1.9  2005/03/28 06:02:45  perev
 * Defence FPE added
 *
 * Revision 1.8  2003/09/02 17:59:35  perev
 * gcc 3.2 updates + WarnOff
 *
 * Revision 1.7  2003/05/01 19:24:31  ullrich
 * Corrected problem in boost().
 *
 * Revision 1.6  1999/10/15 15:56:36  ullrich
 * Changed output format in operator<<, added operator>>
 *
 * Revision 1.5  1999/06/04 18:01:36  ullrich
 * New operators operator() and operator[] which can be used
 * as lvalues.
 *
 * Revision 1.4  1999/04/14 23:12:07  fisyak
 * Add __CINT__ to handle references
 *
 * Revision 1.3  1999/02/17 11:38:36  ullrich
 * Removed specialization for 'long double'.
 *
 * Revision 1.2  1999/02/14 23:11:42  fisyak
 * Fixes for Rootcint
 *
 * Revision 1.1  1999/01/30 03:59:02  fisyak
 * Root Version of AliFemtoarClassLibrary
 *
 * Revision 1.1  1999/01/23 00:27:52  ullrich
 * Initial Revision
 *
 **************************************************************************/
/*//
//// General class for a Lorentz four-vector
///*/
#ifndef ST_LORENTZ_VECTOR_HH
#define ST_LORENTZ_VECTOR_HH

#include "AliFemtoThreeVector.h"
template<class T> class AliFemtoLorentzVector {
public:
    AliFemtoLorentzVector(T = 0, T = 0, T = 0, T = 0);
    virtual ~AliFemtoLorentzVector();
    
#ifndef ST_NO_MEMBER_TEMPLATES
    template<class X> AliFemtoLorentzVector(const AliFemtoThreeVector<X>&, T);
    template<class X> AliFemtoLorentzVector(T, const AliFemtoThreeVector<X>&);   

    template<class X> AliFemtoLorentzVector(const AliFemtoLorentzVector<X>&);
    template<class X> AliFemtoLorentzVector<T>& operator=(const AliFemtoLorentzVector<X>&);
    // AliFemtoLorentzVector(const AliFemtoLorentzVector<T>&);                use default
    // AliFemtoLorentzVector<T>& operator=(const AliFemtoLorentzVector<T>&);  use default
#else
    AliFemtoLorentzVector(const AliFemtoThreeVector<float>&, T);
    AliFemtoLorentzVector(T, const AliFemtoThreeVector<float>&);   
    AliFemtoLorentzVector(const AliFemtoLorentzVector<float>&);
    
    AliFemtoLorentzVector(const AliFemtoThreeVector<double>&, T);
    AliFemtoLorentzVector(T, const AliFemtoThreeVector<double>&);   
    AliFemtoLorentzVector(const AliFemtoLorentzVector<double>&);
        
    AliFemtoLorentzVector<T>& operator=(const AliFemtoLorentzVector<float>&);
    AliFemtoLorentzVector<T>& operator=(const AliFemtoLorentzVector<double>&);
#endif
    
    T x()                     const;
    T y()                     const;
    T z()                     const;
    T t()                     const;
    T px()                    const;
    T py()                    const;
    T pz()                    const;
    T e()                     const;
    T operator()  (size_t)    const;
    T operator[]  (size_t)    const;
    
    T& operator()  (size_t);
    T& operator[]  (size_t);

    const AliFemtoThreeVector<T>& vect() const;    
    
    void SetX(T);
    void SetY(T);
    void SetZ(T);
    void SetPx(T);
    void SetPy(T);
    void SetPz(T);
    void SetE(T);
    void SetT(T);
    
#ifndef ST_NO_MEMBER_TEMPLATES
    template <class X> void SetVect(const AliFemtoThreeVector<X>&);
#else
    void SetVect(const AliFemtoThreeVector<float>&);
    void SetVect(const AliFemtoThreeVector<double>&);
#endif   

    T Perp()               const;
    T Perp2()              const;
    T PseudoRapidity()     const;
    T Phi()                const;
    T Theta()              const;
    T CosTheta()           const;
    
    T Plus()               const;
    T Minus()              const;
    
    T m()                  const; 
    T m2()                 const; 
    T mt()                 const;
    T mt2()                const;
    T Rapidity()           const;
    
#ifndef ST_NO_MEMBER_TEMPLATES
    template<class X> AliFemtoLorentzVector<T> boost(const AliFemtoLorentzVector<X>&) const;
#else
    AliFemtoLorentzVector<T> boost(const AliFemtoLorentzVector<float>&) const;
    AliFemtoLorentzVector<T> boost(const AliFemtoLorentzVector<double>&) const;
#endif   
    
    AliFemtoLorentzVector<T>  operator- ();
    AliFemtoLorentzVector<T>  operator+ ();
    AliFemtoLorentzVector<T>& operator*= (double);
    AliFemtoLorentzVector<T>& operator/= (double);

#ifndef ST_NO_MEMBER_TEMPLATES
    template<class X> bool operator == (const AliFemtoLorentzVector<X>&) const;
    template<class X> bool operator != (const AliFemtoLorentzVector<X>&) const;
    template<class X> AliFemtoLorentzVector<T>& operator+= (const AliFemtoLorentzVector<X>&);
    template<class X> AliFemtoLorentzVector<T>& operator-= (const AliFemtoLorentzVector<X>&);
#else    
    bool operator == (const AliFemtoLorentzVector<float>&) const;
    bool operator != (const AliFemtoLorentzVector<float>&) const;
    bool operator == (const AliFemtoLorentzVector<double>&) const;
    bool operator != (const AliFemtoLorentzVector<double>&) const;

    AliFemtoLorentzVector<T>& operator+= (const AliFemtoLorentzVector<float>&);
    AliFemtoLorentzVector<T>& operator-= (const AliFemtoLorentzVector<float>&);
    AliFemtoLorentzVector<T>& operator+= (const AliFemtoLorentzVector<double>&);
    AliFemtoLorentzVector<T>& operator-= (const AliFemtoLorentzVector<double>&);
#endif

protected:
    AliFemtoThreeVector<T> fThreeVector; // The three-vector component
    T	             fX4;                // The fourth component
#ifdef __ROOT__
  ClassDef(AliFemtoLorentzVector,3)
#endif
};
#ifndef __CINT__
//
//        Implementation of member functions
//
template<class T>
AliFemtoLorentzVector<T>::AliFemtoLorentzVector(T x, T y, T z, T t)
    : fThreeVector(x, y, z), fX4(t) { /* nop */ }

template<class T>
AliFemtoLorentzVector<T>::~AliFemtoLorentzVector() { /* nopt */ }    

template<class T>
const AliFemtoThreeVector<T>& AliFemtoLorentzVector<T>::vect() const 
{
    return fThreeVector;
}

template<class T>
T AliFemtoLorentzVector<T>::m2() const
{
    return (fX4*fX4 - fThreeVector*fThreeVector);    
}

template<class T>
T AliFemtoLorentzVector<T>::Plus() const { return (e() + pz()); }

template<class T>
T AliFemtoLorentzVector<T>::Minus() const { return (e() - pz()); }

template<class T>
T AliFemtoLorentzVector<T>::m() const
{
    T mass2 = m2();
    if (mass2 < 0)
	return -::sqrt(-mass2);
    else
	return ::sqrt(mass2);
}

template<class T>
T AliFemtoLorentzVector<T>::mt2() const
{
    return this->Perp2() + m2();
}

template<class T>
T AliFemtoLorentzVector<T>::mt() const
{
    //
    // change to more optimal code ?
    // return e()*e() - pz()*pz();
    T massPerp2 = mt2();
    if (massPerp2 < 0)
	return -::sqrt(-massPerp2);
    else
	return ::sqrt(massPerp2);
}

template<class T>
void AliFemtoLorentzVector<T>::SetPx(T x) {fThreeVector.SetX(x);}

template<class T>
void AliFemtoLorentzVector<T>::SetPy(T y) {fThreeVector.SetY(y);}

template<class T>
void AliFemtoLorentzVector<T>::SetPz(T z) {fThreeVector.SetZ(z);}

template<class T>
void AliFemtoLorentzVector<T>::SetX(T x) {fThreeVector.SetX(x);}

template<class T>
void AliFemtoLorentzVector<T>::SetY(T y) {fThreeVector.SetY(y);}

template<class T>
void AliFemtoLorentzVector<T>::SetZ(T z) {fThreeVector.SetZ(z);}

template<class T>
void AliFemtoLorentzVector<T>::SetT(T t) {fX4 = t;}

template<class T>
void AliFemtoLorentzVector<T>::SetE(T e) {fX4 = e;}

template<class T>
T AliFemtoLorentzVector<T>::x() const {return fThreeVector.x();}

template<class T>
T AliFemtoLorentzVector<T>::y() const {return fThreeVector.y();}

template<class T>
T AliFemtoLorentzVector<T>::z() const {return fThreeVector.z();}

template<class T>
T AliFemtoLorentzVector<T>::px() const {return fThreeVector.x();}

template<class T>
T AliFemtoLorentzVector<T>::py() const {return fThreeVector.y();}

template<class T>
T AliFemtoLorentzVector<T>::pz() const {return fThreeVector.z();}

template<class T>
T AliFemtoLorentzVector<T>::e() const {return fX4;}

template<class T>
T AliFemtoLorentzVector<T>::t() const {return fX4;}

template<class T>
T AliFemtoLorentzVector<T>::Perp() const {return fThreeVector.Perp();}

template<class T>
T AliFemtoLorentzVector<T>::Perp2() const {return fThreeVector.Perp2();}

template<class T>
T AliFemtoLorentzVector<T>::PseudoRapidity() const {return fThreeVector.PseudoRapidity();}

template<class T>
T AliFemtoLorentzVector<T>::Phi() const {return fThreeVector.Phi();}

template<class T>
T AliFemtoLorentzVector<T>::Theta() const {return fThreeVector.Theta();}

template<class T>
T AliFemtoLorentzVector<T>::CosTheta() const {return fThreeVector.CosTheta();}

template<class T>
T AliFemtoLorentzVector<T>::operator() (size_t i) const
{
    if (i < 3)
        return fThreeVector(i);
    else if (i == 3)
        return fX4;
    else {
#ifndef ST_NO_EXCEPTIONS
      throw out_of_range("AliFemtoLorentzVector<T>::operator(): bad index");  
#else
      cerr << "AliFemtoLorentzVector<T>::operator(): bad index." << endl;
#endif
      return 0;
    }
}

template<class T>
T& AliFemtoLorentzVector<T>::operator() (size_t i)
{
    if (i < 3)
        return fThreeVector(i);
    else if (i == 3)
        return fX4;
    else {
#ifndef ST_NO_EXCEPTIONS
      throw out_of_range("AliFemtoLorentzVector<T>::operator(): bad index");  
#else
      cerr << "AliFemtoLorentzVector<T>::operator(): bad index." << endl;
#endif
      return fX4;
    }
}

template<class T>
T AliFemtoLorentzVector<T>::operator[] (size_t i) const
{
    if (i < 3)
        return fThreeVector[i];
    else if (i == 3)
        return fX4;
    else {
#ifndef ST_NO_EXCEPTIONS
      throw out_of_range("AliFemtoLorentzVector<T>::operator[]: bad index"); 
#else
      cerr << "AliFemtoLorentzVector<T>::operator[]: bad index." << endl;
#endif
      return 0;
    }
}

template<class T>
T& AliFemtoLorentzVector<T>::operator[] (size_t i)
{
    if (i < 3)
        return fThreeVector[i];
    else if (i == 3)
        return fX4;
    else {
#ifndef ST_NO_EXCEPTIONS
      throw out_of_range("AliFemtoLorentzVector<T>::operator[]: bad index"); 
#else
      cerr << "AliFemtoLorentzVector<T>::operator[]: bad index." << endl;
#endif
      return fX4;
    }
}

template<class T>
T AliFemtoLorentzVector<T>::Rapidity() const
{
    return 0.5*::log((fX4+fThreeVector.z())/(fX4-fThreeVector.z())+1e-20);
}

template<class T>
AliFemtoLorentzVector<T> AliFemtoLorentzVector<T>::operator- ()
{
    return AliFemtoLorentzVector<T>(-fX4,-fThreeVector);
}

template<class T>
AliFemtoLorentzVector<T> AliFemtoLorentzVector<T>::operator+ ()
{
    return *this;
}

template<class T>
AliFemtoLorentzVector<T>& AliFemtoLorentzVector<T>::operator*= (double c)
{
    fThreeVector *= c;
    fX4 *= c;
    return *this;
}

template<class T>
AliFemtoLorentzVector<T>& AliFemtoLorentzVector<T>::operator/= (double c)
{
    fThreeVector /= c;
    fX4 /= c;
    return *this;
}

#ifndef ST_NO_MEMBER_TEMPLATES
#ifndef WIN32

template<class T>
template<class X>
AliFemtoLorentzVector<T>::AliFemtoLorentzVector(const AliFemtoThreeVector<X> &vec, T t)
	: fThreeVector(vec), fX4(t) { /* nop */ }

template<class T>
template<class X>
AliFemtoLorentzVector<T>::AliFemtoLorentzVector(T t, const AliFemtoThreeVector<X> &vec)
	: fThreeVector(vec), fX4(t) { /* nop */ }

template<class T>
template<class X>
AliFemtoLorentzVector<T>::AliFemtoLorentzVector(const AliFemtoLorentzVector<X> &vec)
	: fThreeVector(vec.vect()), fX4(vec.t()) { /* nop */ }

template<class T>
template<class X>
AliFemtoLorentzVector<T>
AliFemtoLorentzVector<T>::boost(const AliFemtoLorentzVector<X>& pframe) const
{
    T mass               = abs(pframe);
    AliFemtoThreeVector<T> eta = (-1./mass)*pframe.vect();            // gamma*beta
    T gamma              = fabs(pframe.e())/mass;
    AliFemtoThreeVector<T> pl  = ((this->vect()*eta)/(eta*eta))*eta;  // longitudinal momentum
    return AliFemtoLorentzVector<T>(gamma*this->e() - this->vect()*eta,
                              this->vect() + (gamma-1.)*pl - this->e()*eta);
}

template<class T>
template<class X>
void AliFemtoLorentzVector<T>::SetVect(const AliFemtoThreeVector<X>& v)
{
    fThreeVector = v;
}

template<class T>
template<class X>
AliFemtoLorentzVector<T>&
AliFemtoLorentzVector<T>::operator=(const AliFemtoLorentzVector<X>& vec)
{
    fThreeVector = vec.vect();
    fX4 = vec.t();
    return *this;
}

template<class T>
template<class X>
bool
AliFemtoLorentzVector<T>::operator== (const AliFemtoLorentzVector<X>& v) const
{
    return (fThreeVector == v.vect()) && (fX4 == v.t());
}

template<class T>
template<class X>
bool
AliFemtoLorentzVector<T>::operator!= (const AliFemtoLorentzVector<X>& v) const
{
    return !(*this == v);
}

template<class T>
template<class X>
AliFemtoLorentzVector<T>&
AliFemtoLorentzVector<T>::operator+= (const AliFemtoLorentzVector<X>& v)
{
    fThreeVector += v.vect();
    fX4 += v.t();
    return *this;
}

template<class T>
template<class X>
AliFemtoLorentzVector<T>&
AliFemtoLorentzVector<T>::operator-= (const AliFemtoLorentzVector<X>& v)
{
    fThreeVector -= v.vect();
    fX4 -= v.t();
    return *this;
}

#endif 
#else

template<class T>
AliFemtoLorentzVector<T>::AliFemtoLorentzVector(const AliFemtoThreeVector<float> &vec, T t)
	: fThreeVector(vec), fX4(t) { /* nop */ }

template<class T>
AliFemtoLorentzVector<T>::AliFemtoLorentzVector(const AliFemtoThreeVector<double> &vec, T t)
	: fThreeVector(vec), fX4(t) { /* nop */ }

template<class T>
AliFemtoLorentzVector<T>::AliFemtoLorentzVector(T t, const AliFemtoThreeVector<float> &vec)
	: fThreeVector(vec), fX4(t) { /* nop */ }

template<class T>
AliFemtoLorentzVector<T>::AliFemtoLorentzVector(T t, const AliFemtoThreeVector<double> &vec)
	: fThreeVector(vec), fX4(t) { /* nop */ }

template<class T>
AliFemtoLorentzVector<T>::AliFemtoLorentzVector(const AliFemtoLorentzVector<float> &vec)
	: fThreeVector(vec.vect()), fX4(vec.t()) { /* nop */ }
    
template<class T>
AliFemtoLorentzVector<T>::AliFemtoLorentzVector(const AliFemtoLorentzVector<double> &vec)
	: fThreeVector(vec.vect()), fX4(vec.t()) { /* nop */ }
    
template<class T>
AliFemtoLorentzVector<T>
AliFemtoLorentzVector<T>::boost(const AliFemtoLorentzVector<float>& pframe) const
{
    T mass               = abs(pframe);
    AliFemtoThreeVector<T> eta = (-1./mass)*pframe.vect();            // gamma*beta
    T gamma              = fabs(pframe.e())/mass;
    AliFemtoThreeVector<T> pl  = ((this->vect()*eta)/(eta*eta))*eta;  // longitudinal momentum
    return AliFemtoLorentzVector<T>(gamma*this->e() - this->vect()*eta,
                              this->vect() + (gamma-1.)*pl - this->e()*eta);
}

template<class T>
AliFemtoLorentzVector<T>
AliFemtoLorentzVector<T>::boost(const AliFemtoLorentzVector<double>& pframe) const
{
    T mass               = abs(pframe);
    AliFemtoThreeVector<T> eta = (-1./mass)*pframe.vect();            // gamma*beta
    T gamma              = fabs(pframe.e())/mass;
    AliFemtoThreeVector<T> pl  = ((this->vect()*eta)/(eta*eta))*eta;  // longitudinal momentum
    return AliFemtoLorentzVector<T>(gamma*this->e() - this->vect()*eta,
                              this->vect() + (gamma-1.)*pl - this->e()*eta);
}

template<class T>
void AliFemtoLorentzVector<T>::SetVect(const AliFemtoThreeVector<float>& v)
{
    fThreeVector = v;
}

template<class T>
void AliFemtoLorentzVector<T>::SetVect(const AliFemtoThreeVector<double>& v)
{
    fThreeVector = v;
}

template<class T>
AliFemtoLorentzVector<T>&
AliFemtoLorentzVector<T>::operator=(const AliFemtoLorentzVector<float>& vec)
{
    fThreeVector = vec.vect();
    fX4 = vec.t();
    return *this;
}

template<class T>
AliFemtoLorentzVector<T>&
AliFemtoLorentzVector<T>::operator=(const AliFemtoLorentzVector<double>& vec)
{
    fThreeVector = vec.vect();
    fX4 = vec.t();
    return *this;
}

template<class T>
bool
AliFemtoLorentzVector<T>::operator== (const AliFemtoLorentzVector<float>& v) const
{
    return (this->vect() == v.vect()) && (fX4 == v.t());
}

template<class T>
bool
AliFemtoLorentzVector<T>::operator== (const AliFemtoLorentzVector<double>& v) const
{
    return (fThreeVector == v.vect()) && (fX4 == v.t());
}

template<class T>
bool
AliFemtoLorentzVector<T>::operator!= (const AliFemtoLorentzVector<float>& v) const
{
    return !(*this == v);
}

template<class T>
bool
AliFemtoLorentzVector<T>::operator!= (const AliFemtoLorentzVector<double>& v) const
{
    return !(*this == v);
}

template<class T>
AliFemtoLorentzVector<T>&
AliFemtoLorentzVector<T>::operator+= (const AliFemtoLorentzVector<float>& v)
{
    fThreeVector += v.vect();
    fX4 += v.t();
    return *this;
}

template<class T>
AliFemtoLorentzVector<T>&
AliFemtoLorentzVector<T>::operator+= (const AliFemtoLorentzVector<double>& v)
{
    fThreeVector += v.vect();
    fX4 += v.t();
    return *this;
}

template<class T>
AliFemtoLorentzVector<T>&
AliFemtoLorentzVector<T>::operator-= (const AliFemtoLorentzVector<float>& v)
{
    fThreeVector -= v.vect();
    fX4 -= v.t();
    return *this;
}

template<class T>
AliFemtoLorentzVector<T>&
AliFemtoLorentzVector<T>::operator-= (const AliFemtoLorentzVector<double>& v)
{
    fThreeVector -= v.vect();
    fX4 -= v.t();
    return *this;
}

#endif // ST_NO_MEMBER_TEMPLATES
#endif /* ! __CINT__ */
#ifdef __CINT__
template<> AliFemtoLorentzVector<double> operator+ (const AliFemtoLorentzVector<double>& v1, const AliFemtoLorentzVector<double>& v2);
template<> AliFemtoLorentzVector<double> operator+ (const AliFemtoLorentzVector<double>& v1, const AliFemtoLorentzVector<float>& v2);
template<> AliFemtoLorentzVector<double> operator+ (const AliFemtoLorentzVector<float>&  v1, const AliFemtoLorentzVector<double>& v2);
template<> AliFemtoLorentzVector<float>  operator+ (const AliFemtoLorentzVector<float>&  v1, const AliFemtoLorentzVector<float>& v2);
template<> AliFemtoLorentzVector<double> operator- (const AliFemtoLorentzVector<double>& v1, const AliFemtoLorentzVector<double>& v2);
template<> AliFemtoLorentzVector<double> operator- (const AliFemtoLorentzVector<double>& v1, const AliFemtoLorentzVector<float>& v2);
template<> AliFemtoLorentzVector<double> operator- (const AliFemtoLorentzVector<float>&  v1, const AliFemtoLorentzVector<double>& v2);
template<> AliFemtoLorentzVector<float>  operator- (const AliFemtoLorentzVector<float>&  v1, const AliFemtoLorentzVector<float>& v2);
template<> AliFemtoLorentzVector<double> operator* (const AliFemtoLorentzVector<double>& v1, const AliFemtoLorentzVector<double>& v2);
template<> AliFemtoLorentzVector<double> operator* (const AliFemtoLorentzVector<double>& v1, const AliFemtoLorentzVector<float>& v2);
template<> AliFemtoLorentzVector<double> operator* (const AliFemtoLorentzVector<float>&  v1, const AliFemtoLorentzVector<double>& v2);
template<> AliFemtoLorentzVector<float>  operator* (const AliFemtoLorentzVector<float>&  v1, const AliFemtoLorentzVector<float>& v2);
template<> AliFemtoLorentzVector<double> operator* (const              double v1, const AliFemtoLorentzVector<double>& v2);
template<> AliFemtoLorentzVector<double> operator* (const              double v1, const AliFemtoLorentzVector<float>&  v2);
template<> AliFemtoLorentzVector<double> operator* (const AliFemtoLorentzVector<double>& v1, const double              v2);
template<> AliFemtoLorentzVector<double> operator* (const AliFemtoLorentzVector<float>&  v1, const double              v2);
template<> AliFemtoLorentzVector<double> operator/ (const AliFemtoLorentzVector<double>& v1, const AliFemtoLorentzVector<double>& v2);
template<> AliFemtoLorentzVector<double> operator/ (const AliFemtoLorentzVector<double>& v1, const AliFemtoLorentzVector<float>& v2);
template<> AliFemtoLorentzVector<double> operator/ (const AliFemtoLorentzVector<float>&  v1, const AliFemtoLorentzVector<double>& v2);
template<> AliFemtoLorentzVector<float>  operator/ (const AliFemtoLorentzVector<float>&  v1, const AliFemtoLorentzVector<float>& v2);
template<> AliFemtoLorentzVector<double> operator/ (const              double v1, const AliFemtoLorentzVector<double>& v2);
template<> AliFemtoLorentzVector<double> operator/ (const              double v1, const AliFemtoLorentzVector<float>&  v2);
template<> AliFemtoLorentzVector<double> operator/ (const AliFemtoLorentzVector<double>& v1, const double              v2);
template<> AliFemtoLorentzVector<double> operator/ (const AliFemtoLorentzVector<float>&  v1, const double              v2);
template<> istream& operator>> (istream& is, const AliFemtoLorentzVector<double>& v);
template<> ostream& operator<< (ostream& os, const AliFemtoLorentzVector<double>& v);
template<> istream& operator>> (istream& is, const AliFemtoLorentzVector<float>& v);
template<> ostream& operator<< (ostream& os, const AliFemtoLorentzVector<float>& v);
template<> double abs(const AliFemtoLorentzVector<double>& v);
template<> float  abs(const AliFemtoLorentzVector<float>& v);
#else
//
//   Non-member operators
//
template<class T, class X>
AliFemtoLorentzVector<T>
operator+ (const AliFemtoLorentzVector<T>& v1, const AliFemtoLorentzVector<X>& v2)
{
    return AliFemtoLorentzVector<T>(v1) += v2;
}

template<class T, class X>
AliFemtoLorentzVector<T>
operator- (const AliFemtoLorentzVector<T>& v1, const AliFemtoLorentzVector<X>& v2)
{
    return AliFemtoLorentzVector<T>(v1) -= v2;
}

template<class T, class X>
T
operator* (const AliFemtoLorentzVector<T>& v1, const AliFemtoLorentzVector<X>& v2)
{
    return v1.t()*v2.t() - v1.vect()*v2.vect();
}

template<class T>
AliFemtoLorentzVector<T>
operator* (const AliFemtoLorentzVector<T>& v, double c)
{
    return AliFemtoLorentzVector<T>(v) *= c;
}

template<class T>
AliFemtoLorentzVector<T> operator* (double c, const AliFemtoLorentzVector<T>& v)
{
    return AliFemtoLorentzVector<T>(v) *= c;
}

template<class T, class X>
AliFemtoLorentzVector<T> operator/ (const AliFemtoLorentzVector<T>& v, X c)
{
    return AliFemtoLorentzVector<T>(v) /= c;
}

template<class T>
ostream& operator<< (ostream& os, const AliFemtoLorentzVector<T>& v)
{
    return os << v.vect() << "\t\t" << v.t();
}

template<class T>
istream&  operator>>(istream& is, AliFemtoLorentzVector<T>& v)
{
    T  x, y, z, t;
    is >> x >> y >> z >> t;
    v.SetX(x);
    v.SetY(y);
    v.SetZ(z);
    v.SetT(t);
    return is;
}

//
//        Non-member functions
//
template<class T>
T abs(const AliFemtoLorentzVector<T>& v) {return v.m();}

#endif /*  __CINT__ */
#endif
