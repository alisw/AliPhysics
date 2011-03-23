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
 * Make AliFmLorentzVector persistent
 *
 * Revision 1.10  2005/07/06 18:49:56  fisyak
 * Replace AliFmHelixD, AliFmLorentzVectorD,AliFmLorentzVectorF,AliFmMatrixD,AliFmMatrixF,AliFmPhysicalHelixD,AliFmThreeVectorD,AliFmThreeVectorF by templated version
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
 * Root Version of AliFmarClassLibrary
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

#include "AliFmThreeVector.h"
template<class T> class AliFmLorentzVector {
public:
    AliFmLorentzVector(T = 0, T = 0, T = 0, T = 0);
    virtual ~AliFmLorentzVector();
    
#ifndef ST_NO_MEMBER_TEMPLATES
    template<class X> AliFmLorentzVector(const AliFmThreeVector<X>&, T);
    template<class X> AliFmLorentzVector(T, const AliFmThreeVector<X>&);   

    template<class X> AliFmLorentzVector(const AliFmLorentzVector<X>&);
    template<class X> AliFmLorentzVector<T>& operator=(const AliFmLorentzVector<X>&);
    // AliFmLorentzVector(const AliFmLorentzVector<T>&);                use default
    // AliFmLorentzVector<T>& operator=(const AliFmLorentzVector<T>&);  use default
#else
    AliFmLorentzVector(const AliFmThreeVector<float>&, T);
    AliFmLorentzVector(T, const AliFmThreeVector<float>&);   
    AliFmLorentzVector(const AliFmLorentzVector<float>&);
    
    AliFmLorentzVector(const AliFmThreeVector<double>&, T);
    AliFmLorentzVector(T, const AliFmThreeVector<double>&);   
    AliFmLorentzVector(const AliFmLorentzVector<double>&);
        
    AliFmLorentzVector<T>& operator=(const AliFmLorentzVector<float>&);
    AliFmLorentzVector<T>& operator=(const AliFmLorentzVector<double>&);
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

    const AliFmThreeVector<T>& vect() const;    
    
    void SetX(T);
    void SetY(T);
    void SetZ(T);
    void SetPx(T);
    void SetPy(T);
    void SetPz(T);
    void SetE(T);
    void SetT(T);
    
#ifndef ST_NO_MEMBER_TEMPLATES
    template <class X> void SetVect(const AliFmThreeVector<X>&);
#else
    void SetVect(const AliFmThreeVector<float>&);
    void SetVect(const AliFmThreeVector<double>&);
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
    template<class X> AliFmLorentzVector<T> boost(const AliFmLorentzVector<X>&) const;
#else
    AliFmLorentzVector<T> boost(const AliFmLorentzVector<float>&) const;
    AliFmLorentzVector<T> boost(const AliFmLorentzVector<double>&) const;
#endif   
    
    AliFmLorentzVector<T>  operator- ();
    AliFmLorentzVector<T>  operator+ ();
    AliFmLorentzVector<T>& operator*= (double);
    AliFmLorentzVector<T>& operator/= (double);

#ifndef ST_NO_MEMBER_TEMPLATES
    template<class X> bool operator == (const AliFmLorentzVector<X>&) const;
    template<class X> bool operator != (const AliFmLorentzVector<X>&) const;
    template<class X> AliFmLorentzVector<T>& operator+= (const AliFmLorentzVector<X>&);
    template<class X> AliFmLorentzVector<T>& operator-= (const AliFmLorentzVector<X>&);
#else    
    bool operator == (const AliFmLorentzVector<float>&) const;
    bool operator != (const AliFmLorentzVector<float>&) const;
    bool operator == (const AliFmLorentzVector<double>&) const;
    bool operator != (const AliFmLorentzVector<double>&) const;

    AliFmLorentzVector<T>& operator+= (const AliFmLorentzVector<float>&);
    AliFmLorentzVector<T>& operator-= (const AliFmLorentzVector<float>&);
    AliFmLorentzVector<T>& operator+= (const AliFmLorentzVector<double>&);
    AliFmLorentzVector<T>& operator-= (const AliFmLorentzVector<double>&);
#endif

protected:
    AliFmThreeVector<T> fThreeVector; // The spatial three-vector
    T	             fX4;             // The fourth components
#ifdef __ROOT__
  ClassDef(AliFmLorentzVector,3)
#endif
};
#ifndef __CINT__
//
//        Implementation of member functions
//
template<class T>
AliFmLorentzVector<T>::AliFmLorentzVector(T ax, T ay, T az, T at)
    : fThreeVector(ax, ay, az), fX4(at) { /* nop */ }

template<class T>
AliFmLorentzVector<T>::~AliFmLorentzVector() { /* nopt */ }    

template<class T>
const AliFmThreeVector<T>& AliFmLorentzVector<T>::vect() const 
{
    return fThreeVector;
}

template<class T>
T AliFmLorentzVector<T>::m2() const
{
    return (fX4*fX4 - fThreeVector*fThreeVector);    
}

template<class T>
T AliFmLorentzVector<T>::Plus() const { return (e() + pz()); }

template<class T>
T AliFmLorentzVector<T>::Minus() const { return (e() - pz()); }

template<class T>
T AliFmLorentzVector<T>::m() const
{
    T mass2 = m2();
    if (mass2 < 0)
	return -::sqrt(-mass2);
    else
	return ::sqrt(mass2);
}

template<class T>
T AliFmLorentzVector<T>::mt2() const
{
    return this->Perp2() + m2();
}

template<class T>
T AliFmLorentzVector<T>::mt() const
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
void AliFmLorentzVector<T>::SetPx(T ax) {fThreeVector.SetX(ax);}

template<class T>
void AliFmLorentzVector<T>::SetPy(T ay) {fThreeVector.SetY(ay);}

template<class T>
void AliFmLorentzVector<T>::SetPz(T az) {fThreeVector.SetZ(az);}

template<class T>
void AliFmLorentzVector<T>::SetX(T ax) {fThreeVector.SetX(ax);}

template<class T>
void AliFmLorentzVector<T>::SetY(T ay) {fThreeVector.SetY(ay);}

template<class T>
void AliFmLorentzVector<T>::SetZ(T az) {fThreeVector.SetZ(az);}

template<class T>
void AliFmLorentzVector<T>::SetT(T at) {fX4 = at;}

template<class T>
void AliFmLorentzVector<T>::SetE(T ae) {fX4 = ae;}

template<class T>
T AliFmLorentzVector<T>::x() const {return fThreeVector.x();}

template<class T>
T AliFmLorentzVector<T>::y() const {return fThreeVector.y();}

template<class T>
T AliFmLorentzVector<T>::z() const {return fThreeVector.z();}

template<class T>
T AliFmLorentzVector<T>::px() const {return fThreeVector.x();}

template<class T>
T AliFmLorentzVector<T>::py() const {return fThreeVector.y();}

template<class T>
T AliFmLorentzVector<T>::pz() const {return fThreeVector.z();}

template<class T>
T AliFmLorentzVector<T>::e() const {return fX4;}

template<class T>
T AliFmLorentzVector<T>::t() const {return fX4;}

template<class T>
T AliFmLorentzVector<T>::Perp() const {return fThreeVector.Perp();}

template<class T>
T AliFmLorentzVector<T>::Perp2() const {return fThreeVector.Perp2();}

template<class T>
T AliFmLorentzVector<T>::PseudoRapidity() const {return fThreeVector.PseudoRapidity();}

template<class T>
T AliFmLorentzVector<T>::Phi() const {return fThreeVector.Phi();}

template<class T>
T AliFmLorentzVector<T>::Theta() const {return fThreeVector.Theta();}

template<class T>
T AliFmLorentzVector<T>::CosTheta() const {return fThreeVector.CosTheta();}

template<class T>
T AliFmLorentzVector<T>::operator() (size_t i) const
{
    if (i < 3)
        return fThreeVector(i);
    else if (i == 3)
        return fX4;
    else {
#ifndef ST_NO_EXCEPTIONS
      throw out_of_range("AliFmLorentzVector<T>::operator(): bad index");  
#else
      cerr << "AliFmLorentzVector<T>::operator(): bad index." << endl;
#endif
      return 0;
    }
}

template<class T>
T& AliFmLorentzVector<T>::operator() (size_t i)
{
    if (i < 3)
        return fThreeVector(i);
    else if (i == 3)
        return fX4;
    else {
#ifndef ST_NO_EXCEPTIONS
      throw out_of_range("AliFmLorentzVector<T>::operator(): bad index");  
#else
      cerr << "AliFmLorentzVector<T>::operator(): bad index." << endl;
#endif
      return fX4;
    }
}

template<class T>
T AliFmLorentzVector<T>::operator[] (size_t i) const
{
    if (i < 3)
        return fThreeVector[i];
    else if (i == 3)
        return fX4;
    else {
#ifndef ST_NO_EXCEPTIONS
      throw out_of_range("AliFmLorentzVector<T>::operator[]: bad index"); 
#else
      cerr << "AliFmLorentzVector<T>::operator[]: bad index." << endl;
#endif
      return 0;
    }
}

template<class T>
T& AliFmLorentzVector<T>::operator[] (size_t i)
{
    if (i < 3)
        return fThreeVector[i];
    else if (i == 3)
        return fX4;
    else {
#ifndef ST_NO_EXCEPTIONS
      throw out_of_range("AliFmLorentzVector<T>::operator[]: bad index"); 
#else
      cerr << "AliFmLorentzVector<T>::operator[]: bad index." << endl;
#endif
      return fX4;
    }
}

template<class T>
T AliFmLorentzVector<T>::Rapidity() const
{
    return 0.5*::log((fX4+fThreeVector.z())/(fX4-fThreeVector.z())+1e-20);
}

template<class T>
AliFmLorentzVector<T> AliFmLorentzVector<T>::operator- ()
{
    return AliFmLorentzVector<T>(-fX4,-fThreeVector);
}

template<class T>
AliFmLorentzVector<T> AliFmLorentzVector<T>::operator+ ()
{
    return *this;
}

template<class T>
AliFmLorentzVector<T>& AliFmLorentzVector<T>::operator*= (double c)
{
    fThreeVector *= c;
    fX4 *= c;
    return *this;
}

template<class T>
AliFmLorentzVector<T>& AliFmLorentzVector<T>::operator/= (double c)
{
    fThreeVector /= c;
    fX4 /= c;
    return *this;
}

#ifndef ST_NO_MEMBER_TEMPLATES
#ifndef WIN32

template<class T>
template<class X>
AliFmLorentzVector<T>::AliFmLorentzVector(const AliFmThreeVector<X> &vec, T at)
	: fThreeVector(vec), fX4(at) { /* nop */ }

template<class T>
template<class X>
AliFmLorentzVector<T>::AliFmLorentzVector(T at, const AliFmThreeVector<X> &vec)
	: fThreeVector(vec), fX4(at) { /* nop */ }

template<class T>
template<class X>
AliFmLorentzVector<T>::AliFmLorentzVector(const AliFmLorentzVector<X> &vec)
	: fThreeVector(vec.vect()), fX4(vec.t()) { /* nop */ }

template<class T>
template<class X>
AliFmLorentzVector<T>
AliFmLorentzVector<T>::boost(const AliFmLorentzVector<X>& pframe) const
{
    T mass               = abs(pframe);
    AliFmThreeVector<T> eta = (-1./mass)*pframe.vect();            // gamma*beta
    T gamma              = fabs(pframe.e())/mass;
    AliFmThreeVector<T> pl  = ((this->vect()*eta)/(eta*eta))*eta;  // longitudinal momentum
    return AliFmLorentzVector<T>(gamma*this->e() - this->vect()*eta,
                              this->vect() + (gamma-1.)*pl - this->e()*eta);
}

template<class T>
template<class X>
void AliFmLorentzVector<T>::SetVect(const AliFmThreeVector<X>& v)
{
    fThreeVector = v;
}

template<class T>
template<class X>
AliFmLorentzVector<T>&
AliFmLorentzVector<T>::operator=(const AliFmLorentzVector<X>& vec)
{
    fThreeVector = vec.vect();
    fX4 = vec.t();
    return *this;
}

template<class T>
template<class X>
bool
AliFmLorentzVector<T>::operator== (const AliFmLorentzVector<X>& v) const
{
    return (fThreeVector == v.vect()) && (fX4 == v.t());
}

template<class T>
template<class X>
bool
AliFmLorentzVector<T>::operator!= (const AliFmLorentzVector<X>& v) const
{
    return !(*this == v);
}

template<class T>
template<class X>
AliFmLorentzVector<T>&
AliFmLorentzVector<T>::operator+= (const AliFmLorentzVector<X>& v)
{
    fThreeVector += v.vect();
    fX4 += v.t();
    return *this;
}

template<class T>
template<class X>
AliFmLorentzVector<T>&
AliFmLorentzVector<T>::operator-= (const AliFmLorentzVector<X>& v)
{
    fThreeVector -= v.vect();
    fX4 -= v.t();
    return *this;
}

#endif 
#else

template<class T>
AliFmLorentzVector<T>::AliFmLorentzVector(const AliFmThreeVector<float> &vec, T t)
	: fThreeVector(vec), fX4(t) { /* nop */ }

template<class T>
AliFmLorentzVector<T>::AliFmLorentzVector(const AliFmThreeVector<double> &vec, T t)
	: fThreeVector(vec), fX4(t) { /* nop */ }

template<class T>
AliFmLorentzVector<T>::AliFmLorentzVector(T t, const AliFmThreeVector<float> &vec)
	: fThreeVector(vec), fX4(t) { /* nop */ }

template<class T>
AliFmLorentzVector<T>::AliFmLorentzVector(T t, const AliFmThreeVector<double> &vec)
	: fThreeVector(vec), fX4(t) { /* nop */ }

template<class T>
AliFmLorentzVector<T>::AliFmLorentzVector(const AliFmLorentzVector<float> &vec)
	: fThreeVector(vec.vect()), fX4(vec.t()) { /* nop */ }
    
template<class T>
AliFmLorentzVector<T>::AliFmLorentzVector(const AliFmLorentzVector<double> &vec)
	: fThreeVector(vec.vect()), fX4(vec.t()) { /* nop */ }
    
template<class T>
AliFmLorentzVector<T>
AliFmLorentzVector<T>::boost(const AliFmLorentzVector<float>& pframe) const
{
    T mass               = abs(pframe);
    AliFmThreeVector<T> eta = (-1./mass)*pframe.vect();            // gamma*beta
    T gamma              = fabs(pframe.e())/mass;
    AliFmThreeVector<T> pl  = ((this->vect()*eta)/(eta*eta))*eta;  // longitudinal momentum
    return AliFmLorentzVector<T>(gamma*this->e() - this->vect()*eta,
                              this->vect() + (gamma-1.)*pl - this->e()*eta);
}

template<class T>
AliFmLorentzVector<T>
AliFmLorentzVector<T>::boost(const AliFmLorentzVector<double>& pframe) const
{
    T mass               = abs(pframe);
    AliFmThreeVector<T> eta = (-1./mass)*pframe.vect();            // gamma*beta
    T gamma              = fabs(pframe.e())/mass;
    AliFmThreeVector<T> pl  = ((this->vect()*eta)/(eta*eta))*eta;  // longitudinal momentum
    return AliFmLorentzVector<T>(gamma*this->e() - this->vect()*eta,
                              this->vect() + (gamma-1.)*pl - this->e()*eta);
}

template<class T>
void AliFmLorentzVector<T>::SetVect(const AliFmThreeVector<float>& v)
{
    fThreeVector = v;
}

template<class T>
void AliFmLorentzVector<T>::SetVect(const AliFmThreeVector<double>& v)
{
    fThreeVector = v;
}

template<class T>
AliFmLorentzVector<T>&
AliFmLorentzVector<T>::operator=(const AliFmLorentzVector<float>& vec)
{
    fThreeVector = vec.vect();
    fX4 = vec.t();
    return *this;
}

template<class T>
AliFmLorentzVector<T>&
AliFmLorentzVector<T>::operator=(const AliFmLorentzVector<double>& vec)
{
    fThreeVector = vec.vect();
    fX4 = vec.t();
    return *this;
}

template<class T>
bool
AliFmLorentzVector<T>::operator== (const AliFmLorentzVector<float>& v) const
{
    return (this->vect() == v.vect()) && (fX4 == v.t());
}

template<class T>
bool
AliFmLorentzVector<T>::operator== (const AliFmLorentzVector<double>& v) const
{
    return (fThreeVector == v.vect()) && (fX4 == v.t());
}

template<class T>
bool
AliFmLorentzVector<T>::operator!= (const AliFmLorentzVector<float>& v) const
{
    return !(*this == v);
}

template<class T>
bool
AliFmLorentzVector<T>::operator!= (const AliFmLorentzVector<double>& v) const
{
    return !(*this == v);
}

template<class T>
AliFmLorentzVector<T>&
AliFmLorentzVector<T>::operator+= (const AliFmLorentzVector<float>& v)
{
    fThreeVector += v.vect();
    fX4 += v.t();
    return *this;
}

template<class T>
AliFmLorentzVector<T>&
AliFmLorentzVector<T>::operator+= (const AliFmLorentzVector<double>& v)
{
    fThreeVector += v.vect();
    fX4 += v.t();
    return *this;
}

template<class T>
AliFmLorentzVector<T>&
AliFmLorentzVector<T>::operator-= (const AliFmLorentzVector<float>& v)
{
    fThreeVector -= v.vect();
    fX4 -= v.t();
    return *this;
}

template<class T>
AliFmLorentzVector<T>&
AliFmLorentzVector<T>::operator-= (const AliFmLorentzVector<double>& v)
{
    fThreeVector -= v.vect();
    fX4 -= v.t();
    return *this;
}

#endif // ST_NO_MEMBER_TEMPLATES
#endif /* ! __CINT__ */
#ifdef __CINT__
template<> AliFmLorentzVector<double> operator+ (const AliFmLorentzVector<double>& v1, const AliFmLorentzVector<double>& v2);
template<> AliFmLorentzVector<double> operator+ (const AliFmLorentzVector<double>& v1, const AliFmLorentzVector<float>& v2);
template<> AliFmLorentzVector<double> operator+ (const AliFmLorentzVector<float>&  v1, const AliFmLorentzVector<double>& v2);
template<> AliFmLorentzVector<float>  operator+ (const AliFmLorentzVector<float>&  v1, const AliFmLorentzVector<float>& v2);
template<> AliFmLorentzVector<double> operator- (const AliFmLorentzVector<double>& v1, const AliFmLorentzVector<double>& v2);
template<> AliFmLorentzVector<double> operator- (const AliFmLorentzVector<double>& v1, const AliFmLorentzVector<float>& v2);
template<> AliFmLorentzVector<double> operator- (const AliFmLorentzVector<float>&  v1, const AliFmLorentzVector<double>& v2);
template<> AliFmLorentzVector<float>  operator- (const AliFmLorentzVector<float>&  v1, const AliFmLorentzVector<float>& v2);
template<> AliFmLorentzVector<double> operator* (const AliFmLorentzVector<double>& v1, const AliFmLorentzVector<double>& v2);
template<> AliFmLorentzVector<double> operator* (const AliFmLorentzVector<double>& v1, const AliFmLorentzVector<float>& v2);
template<> AliFmLorentzVector<double> operator* (const AliFmLorentzVector<float>&  v1, const AliFmLorentzVector<double>& v2);
template<> AliFmLorentzVector<float>  operator* (const AliFmLorentzVector<float>&  v1, const AliFmLorentzVector<float>& v2);
template<> AliFmLorentzVector<double> operator* (const              double v1, const AliFmLorentzVector<double>& v2);
template<> AliFmLorentzVector<double> operator* (const              double v1, const AliFmLorentzVector<float>&  v2);
template<> AliFmLorentzVector<double> operator* (const AliFmLorentzVector<double>& v1, const double              v2);
template<> AliFmLorentzVector<double> operator* (const AliFmLorentzVector<float>&  v1, const double              v2);
template<> AliFmLorentzVector<double> operator/ (const AliFmLorentzVector<double>& v1, const AliFmLorentzVector<double>& v2);
template<> AliFmLorentzVector<double> operator/ (const AliFmLorentzVector<double>& v1, const AliFmLorentzVector<float>& v2);
template<> AliFmLorentzVector<double> operator/ (const AliFmLorentzVector<float>&  v1, const AliFmLorentzVector<double>& v2);
template<> AliFmLorentzVector<float>  operator/ (const AliFmLorentzVector<float>&  v1, const AliFmLorentzVector<float>& v2);
template<> AliFmLorentzVector<double> operator/ (const              double v1, const AliFmLorentzVector<double>& v2);
template<> AliFmLorentzVector<double> operator/ (const              double v1, const AliFmLorentzVector<float>&  v2);
template<> AliFmLorentzVector<double> operator/ (const AliFmLorentzVector<double>& v1, const double              v2);
template<> AliFmLorentzVector<double> operator/ (const AliFmLorentzVector<float>&  v1, const double              v2);
template<> istream& operator>> (istream& is, const AliFmLorentzVector<double>& v);
template<> ostream& operator<< (ostream& os, const AliFmLorentzVector<double>& v);
template<> istream& operator>> (istream& is, const AliFmLorentzVector<float>& v);
template<> ostream& operator<< (ostream& os, const AliFmLorentzVector<float>& v);
template<> double abs(const AliFmLorentzVector<double>& v);
template<> float  abs(const AliFmLorentzVector<float>& v);
#else
//
//   Non-member operators
//
template<class T, class X>
AliFmLorentzVector<T>
operator+ (const AliFmLorentzVector<T>& v1, const AliFmLorentzVector<X>& v2)
{
    return AliFmLorentzVector<T>(v1) += v2;
}

template<class T, class X>
AliFmLorentzVector<T>
operator- (const AliFmLorentzVector<T>& v1, const AliFmLorentzVector<X>& v2)
{
    return AliFmLorentzVector<T>(v1) -= v2;
}

template<class T, class X>
T
operator* (const AliFmLorentzVector<T>& v1, const AliFmLorentzVector<X>& v2)
{
    return v1.t()*v2.t() - v1.vect()*v2.vect();
}

template<class T>
AliFmLorentzVector<T>
operator* (const AliFmLorentzVector<T>& v, double c)
{
    return AliFmLorentzVector<T>(v) *= c;
}

template<class T>
AliFmLorentzVector<T> operator* (double c, const AliFmLorentzVector<T>& v)
{
    return AliFmLorentzVector<T>(v) *= c;
}

template<class T, class X>
AliFmLorentzVector<T> operator/ (const AliFmLorentzVector<T>& v, X c)
{
    return AliFmLorentzVector<T>(v) /= c;
}

template<class T>
ostream& operator<< (ostream& os, const AliFmLorentzVector<T>& v)
{
    return os << v.vect() << "\t\t" << v.t();
}

template<class T>
istream&  operator>>(istream& is, AliFmLorentzVector<T>& v)
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
T abs(const AliFmLorentzVector<T>& v) {return v.m();}

#endif /*  __CINT__ */
#endif
