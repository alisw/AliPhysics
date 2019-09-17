// fjcore -- extracted from FastJet v3.2.1 (http://fastjet.fr)
//
// fjcore constitutes a digest of the main FastJet functionality.
// The files fjcore.hh and fjcore.cc are meant to provide easy access to these 
// core functions, in the form of single files and without the need of a full 
// FastJet installation:
//
//     g++ main.cc fjcore.cc
// 
// with main.cc including fjcore.hh.
//
// A fortran interface, fjcorefortran.cc, is also provided. See the example 
// and the Makefile for instructions.
//
// The results are expected to be identical to those obtained by linking to
// the full FastJet distribution.
//
// NOTE THAT, IN ORDER TO MAKE IT POSSIBLE FOR FJCORE AND THE FULL FASTJET
// TO COEXIST, THE FORMER USES THE "fjcore" NAMESPACE INSTEAD OF "fastjet". 
//
// In particular, fjcore provides:
//
//   - access to all native pp and ee algorithms, kt, anti-kt, C/A.
//     For C/A, the NlnN method is available, while anti-kt and kt
//     are limited to the N^2 one (still the fastest for N < 100k particles)
//   - access to selectors, for implementing cuts and selections
//   - access to all functionalities related to pseudojets (e.g. a jet's
//     structure or user-defined information)
//
// Instead, it does NOT provide:
//
//   - jet areas functionality
//   - background estimation
//   - access to other algorithms via plugins
//   - interface to CGAL
//   - fastjet tools, e.g. filters, taggers
//
// If these functionalities are needed, the full FastJet installation must be
// used. The code will be fully compatible, with the sole replacement of the
// header files and of the fjcore namespace with the fastjet one.
//
// fjcore.hh and fjcore.cc are not meant to be human-readable.
// For documentation, see the full FastJet manual and doxygen at http://fastjet.fr
//
// Like FastJet, fjcore is released under the terms of the GNU General Public
// License version 2 (GPLv2). If you use this code as part of work towards a
// scientific publication, whether directly or contained within another program
// (e.g. Delphes, MadGraph, SpartyJet, Rivet, LHC collaboration software frameworks, 
// etc.), you should include a citation to
// 
//   EPJC72(2012)1896 [arXiv:1111.6097] (FastJet User Manual)
//   and, optionally, Phys.Lett.B641 (2006) 57 [arXiv:hep-ph/0512210]
//
// Copyright (c) 2005-2016, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development and are described in hep-ph/0512210. If you use
//  FastJet as part of work towards a scientific publication, please
//  include a citation to the FastJet paper.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//
#ifndef __FJCORE_HH__
#define __FJCORE_HH__
#define __FJCORE__   // remove all the non-core code (a safekeeper)
#define __FJCORE_DROP_CGAL    // disable CGAL support
#ifndef _INCLUDE_FJCORE_CONFIG_AUTO_H
#define _INCLUDE_FJCORE_CONFIG_AUTO_H 1
#ifndef FJCORE_HAVE_DEMANGLING_SUPPORT 
#endif
#ifndef FJCORE_HAVE_DLFCN_H 
# define FJCORE_HAVE_DLFCN_H  1 
#endif
#ifndef FJCORE_HAVE_EXECINFO_H 
#endif
#ifndef FJCORE_HAVE_GNUCXX_DEPRECATED 
#endif
#ifndef FJCORE_HAVE_INTTYPES_H 
# define FJCORE_HAVE_INTTYPES_H  1 
#endif
#ifndef FJCORE_HAVE_LIBM 
# define FJCORE_HAVE_LIBM  1 
#endif
#ifndef FJCORE_HAVE_MEMORY_H 
# define FJCORE_HAVE_MEMORY_H  1 
#endif
#ifndef FJCORE_HAVE_STDINT_H 
# define FJCORE_HAVE_STDINT_H  1 
#endif
#ifndef FJCORE_HAVE_STDLIB_H 
# define FJCORE_HAVE_STDLIB_H  1 
#endif
#ifndef FJCORE_HAVE_STRINGS_H 
# define FJCORE_HAVE_STRINGS_H  1 
#endif
#ifndef FJCORE_HAVE_STRING_H 
# define FJCORE_HAVE_STRING_H  1 
#endif
#ifndef FJCORE_HAVE_SYS_STAT_H 
# define FJCORE_HAVE_SYS_STAT_H  1 
#endif
#ifndef FJCORE_HAVE_SYS_TYPES_H 
# define FJCORE_HAVE_SYS_TYPES_H  1 
#endif
#ifndef FJCORE_HAVE_UNISTD_H 
# define FJCORE_HAVE_UNISTD_H  1 
#endif
#ifndef FJCORE_LT_OBJDIR 
# define FJCORE_LT_OBJDIR  ".libs/" 
#endif
#ifndef FJCORE_PACKAGE 
# define FJCORE_PACKAGE  "fastjet" 
#endif
#ifndef FJCORE_PACKAGE_BUGREPORT 
# define FJCORE_PACKAGE_BUGREPORT  "" 
#endif
#ifndef FJCORE_PACKAGE_NAME 
# define FJCORE_PACKAGE_NAME  "FastJet" 
#endif
#ifndef FJCORE_PACKAGE_STRING 
# define FJCORE_PACKAGE_STRING  "FastJet 3.2.1" 
#endif
#ifndef FJCORE_PACKAGE_TARNAME 
# define FJCORE_PACKAGE_TARNAME  "fastjet" 
#endif
#ifndef FJCORE_PACKAGE_URL 
# define FJCORE_PACKAGE_URL  "" 
#endif
#ifndef FJCORE_PACKAGE_VERSION 
# define FJCORE_PACKAGE_VERSION  "3.2.1" 
#endif
#ifndef FJCORE_STDC_HEADERS 
# define FJCORE_STDC_HEADERS  1 
#endif
#ifndef FJCORE_VERSION 
# define FJCORE_VERSION  "3.2.1" 
#endif
#ifndef FJCORE_VERSION_MAJOR 
# define FJCORE_VERSION_MAJOR  3 
#endif
#ifndef FJCORE_VERSION_MINOR 
# define FJCORE_VERSION_MINOR  2 
#endif
#ifndef FJCORE_VERSION_NUMBER 
# define FJCORE_VERSION_NUMBER  30201 
#endif
#ifndef FJCORE_VERSION_PATCHLEVEL 
# define FJCORE_VERSION_PATCHLEVEL  1 
#endif
#endif
#ifndef __FJCORE_CONFIG_H__
#define __FJCORE_CONFIG_H__
#endif // __FJCORE_CONFIG_H__
#ifndef __FJCORE_FASTJET_BASE_HH__
#define __FJCORE_FASTJET_BASE_HH__
// TS : enclose the fjcore namespace inside the Pythia8 one
//#define FJCORE_BEGIN_NAMESPACE namespace fjcore {
//#define FJCORE_END_NAMESPACE   }
#define FJCORE_BEGIN_NAMESPACE namespace Pythia8 { namespace fjcore {
#define FJCORE_END_NAMESPACE   }}
#ifdef FJCORE_HAVE_OVERRIDE
# define FJCORE_OVERRIDE  override
#else
# define FJCORE_OVERRIDE  
#endif
#endif // __FJCORE_FASTJET_BASE_HH__
#ifndef __FJCORE_NUMCONSTS__
#define __FJCORE_NUMCONSTS__
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
const double pi = 3.141592653589793238462643383279502884197;
const double twopi = 6.283185307179586476925286766559005768394;
const double pisq  = 9.869604401089358618834490999876151135314;
const double zeta2 = 1.644934066848226436472415166646025189219;
const double zeta3 = 1.202056903159594285399738161511449990765;
const double eulergamma = 0.577215664901532860606512090082402431042;
const double ln2   = 0.693147180559945309417232121458176568076;
FJCORE_END_NAMESPACE
#endif // __FJCORE_NUMCONSTS__
#ifndef __FJCORE_INTERNAL_IS_BASE_HH__
#define __FJCORE_INTERNAL_IS_BASE_HH__
FJCORE_BEGIN_NAMESPACE
template<typename T, T _t>
struct integral_type{
  static const T value = _t;         ///< the value (only member carrying info)
  typedef T value_type;		     ///< a typedef for the type T
  typedef integral_type<T,_t> type;  ///< a typedef for the whole structure
};
template<typename T, T _t>
const T integral_type<T, _t>::value;
typedef integral_type<bool, true>  true_type;  ///< the bool 'true'  value promoted to a type
typedef integral_type<bool, false> false_type; ///< the bool 'false' value promoted to a type
typedef char (&__yes_type)[1]; //< the yes type
typedef char (&__no_type) [2]; //< the no type
template<typename B, typename D>
struct __inheritance_helper{
#if !((_MSC_VER !=0 ) && (_MSC_VER == 1310))   // MSVC 7.1
  template <typename T>
  static __yes_type check_sig(D const volatile *, T);
#else
  static __yes_type check_sig(D const volatile *, long);
#endif
  static __no_type  check_sig(B const volatile *, int);
};
template<typename B, typename D>
struct IsBaseAndDerived{
#if ((_MSC_FULL_VER != 0) && (_MSC_FULL_VER >= 140050000))
#pragma warning(push)
#pragma warning(disable:6334)
#endif
  struct Host{
#if !((_MSC_VER !=0 ) && (_MSC_VER == 1310))
    operator B const volatile *() const;
#else
    operator B const volatile * const&() const;
#endif
    operator D const volatile *();
  };
  static const bool value = ((sizeof(B)!=0) && 
			     (sizeof(D)!=0) && 
			     (sizeof(__inheritance_helper<B,D>::check_sig(Host(), 0)) == sizeof(__yes_type)));
#if ((_MSC_FULL_VER != 0) && (_MSC_FULL_VER >= 140050000))
#pragma warning(pop)
#endif
};
template<class B, class D>
B* cast_if_derived(D* d){
  return IsBaseAndDerived<B,D>::value ? (B*)(d) : 0;
}
FJCORE_END_NAMESPACE
#endif  // __IS_BASE_OF_HH__
#ifndef __FJCORE_FJCORE_DEPRECATED_HH__
#define __FJCORE_FJCORE_DEPRECATED_HH__
#if defined(FJCORE_HAVE_CXX14_DEPRECATED)
# define FJCORE_DEPRECATED               [[deprecated]]
# define FJCORE_DEPRECATED_MSG(message)  [[deprecated(message)]]
#elif defined(FJCORE_HAVE_GNUCXX_DEPRECATED)
# define FJCORE_DEPRECATED               __attribute__((__deprecated__))
# define FJCORE_DEPRECATED_MSG(message)  __attribute__((__deprecated__))
#else
# define FJCORE_DEPRECATED               
# define FJCORE_DEPRECATED_MSG(message) 
#endif
#endif // __FJCORE_FJCORE_DEPRECATED_HH__
#ifndef __FJCORE_SHARED_PTR_HH__
#define __FJCORE_SHARED_PTR_HH__
#include <cstdlib>  // for NULL!!!
#ifdef __FJCORE_USETR1SHAREDPTR
#include <tr1/memory>
#endif // __FJCORE_USETR1SHAREDPTR
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
#ifdef __FJCORE_USETR1SHAREDPTR
template<class T>
class SharedPtr : public std::tr1::shared_ptr<T> {
public:
  SharedPtr() : std::tr1::shared_ptr<T>() {}
  SharedPtr(T * t) : std::tr1::shared_ptr<T>(t) {}
  SharedPtr(const SharedPtr<T> & t) : std::tr1::shared_ptr<T>(t) {}
  #ifdef FJCORE_HAVE_EXPLICIT_FOR_OPERATORS
  explicit
  #endif
  inline operator bool() const {return (this->get()!=NULL);}
  T* operator ()() const{
    return this->get(); // automatically returns NULL when out-of-scope
  }
};
#else // __FJCORE_USETR1SHAREDPTR
template<class T>
class SharedPtr{
public:
  class __SharedCountingPtr;
  SharedPtr() : _ptr(NULL){}
  template<class Y> explicit SharedPtr(Y* ptr){
    _ptr = new __SharedCountingPtr(ptr);
  }
  SharedPtr(SharedPtr const & share) : _ptr(share._get_container()){
    if (_ptr!=NULL) ++(*_ptr);
  }
  ~SharedPtr(){
    if (_ptr==NULL) return;
    _decrease_count();
  }
  void reset(){
    SharedPtr().swap(*this);
  }
  template<class Y> void reset(Y * ptr){
    SharedPtr(ptr).swap(*this);
  }
  template<class Y> void reset(SharedPtr<Y> const & share){
    if (_ptr!=NULL){
      if (_ptr == share._get_container()) return;
      _decrease_count();
    }
    _ptr = share._get_container();  // Note: automatically set it to NULL if share is empty
    if (_ptr!=NULL) ++(*_ptr);
  }
  SharedPtr& operator=(SharedPtr const & share){
    reset(share);
    return *this;
  }
  template<class Y> SharedPtr& operator=(SharedPtr<Y> const & share){
    reset(share);
    return *this;
  }
  FJCORE_DEPRECATED_MSG("Use SharedPtr<T>::get() instead")
  T* operator ()() const{
    if (_ptr==NULL) return NULL;
    return _ptr->get(); // automatically returns NULL when out-of-scope
  }
  inline T& operator*() const{
    return *(_ptr->get());
  }
  inline T* operator->() const{
    if (_ptr==NULL) return NULL;
    return _ptr->get();
  }  
  inline T* get() const{
    if (_ptr==NULL) return NULL;
    return _ptr->get();
  }
  inline bool unique() const{
    return (use_count()==1);
  }
  inline long use_count() const{
    if (_ptr==NULL) return 0;
    return _ptr->use_count(); // automatically returns NULL when out-of-scope
  }
  #ifdef FJCORE_HAVE_EXPLICIT_FOR_OPERATORS
  explicit
  #endif
  inline operator bool() const{
    return (get()!=NULL);
  }
  inline void swap(SharedPtr & share){
    __SharedCountingPtr* share_container = share._ptr;
    share._ptr = _ptr;
    _ptr = share_container;
  }
  void set_count(const long & count){
    if (_ptr==NULL) return;
    _ptr->set_count(count);
  }
  class __SharedCountingPtr{
  public:
    __SharedCountingPtr() : _ptr(NULL), _count(0){}
    template<class Y> explicit __SharedCountingPtr(Y* ptr) : _ptr(ptr), _count(1){}
    ~__SharedCountingPtr(){ 
      if (_ptr!=NULL){ delete _ptr;}
    }
    inline T* get() const {return _ptr;}
    inline long use_count() const {return _count;}
    inline long operator++(){return ++_count;}
    inline long operator--(){return --_count;}
    inline long operator++(int){return _count++;}
    inline long operator--(int){return _count--;}
    void set_count(const long & count){
      _count = count;
    }
  private:
    T *_ptr;              ///< the pointer we're counting the references to
    long _count;  ///< the number of references
  };
private:
  inline __SharedCountingPtr* _get_container() const{
    return _ptr;
  }
  void _decrease_count(){
    (*_ptr)--;
    if (_ptr->use_count()==0)
      delete _ptr; // that automatically deletes the object itself
  }
  __SharedCountingPtr *_ptr;
};
template<class T,class U>
inline bool operator==(SharedPtr<T> const & t, SharedPtr<U> const & u){
  return t.get() == u.get();
}
template<class T,class U>
inline bool operator!=(SharedPtr<T> const & t, SharedPtr<U> const & u){
  return t.get() != u.get();
}
template<class T,class U>
inline bool operator<(SharedPtr<T> const & t, SharedPtr<U> const & u){
  return t.get() < u.get();
}
template<class T>
inline void swap(SharedPtr<T> & a, SharedPtr<T> & b){
  return a.swap(b);
}
template<class T>
inline T* get_pointer(SharedPtr<T> const & t){
  return t.get();
}
#endif // __FJCORE_USETR1SHAREDPTR
FJCORE_END_NAMESPACE      // defined in fastjet/internal/base.hh
#endif   // __FJCORE_SHARED_PTR_HH__
#ifndef __FJCORE_LIMITEDWARNING_HH__
#define __FJCORE_LIMITEDWARNING_HH__
#include <iostream>
#include <string>
#include <list>
FJCORE_BEGIN_NAMESPACE
class LimitedWarning {
public:
  LimitedWarning() : _max_warn(_max_warn_default), _n_warn_so_far(0), _this_warning_summary(0) {}
  LimitedWarning(int max_warn_in) : _max_warn(max_warn_in), _n_warn_so_far(0), _this_warning_summary(0) {}
  void warn(const char * warning) {warn(warning, _default_ostr);}
  void warn(const std::string & warning) {warn(warning.c_str(), _default_ostr);}
  void warn(const char * warning, std::ostream * ostr);
  void warn(const std::string & warning, std::ostream * ostr) {warn(warning.c_str(), ostr);}
  static void set_default_stream(std::ostream * ostr) {
    _default_ostr = ostr;
  }
  static void set_default_max_warn(int max_warn) {
    _max_warn_default = max_warn;
  }
  int max_warn() const {return _max_warn;}
  int n_warn_so_far() const {return _n_warn_so_far;}
  static std::string summary();
private:
  int _max_warn, _n_warn_so_far;
  static int _max_warn_default;
  static std::ostream * _default_ostr;
  typedef std::pair<std::string, unsigned int> Summary;
  static std::list< Summary > _global_warnings_summary;
  Summary * _this_warning_summary;
};
FJCORE_END_NAMESPACE
#endif // __FJCORE_LIMITEDWARNING_HH__
#ifndef __FJCORE_ERROR_HH__
#define __FJCORE_ERROR_HH__
#include<iostream>
#include<string>
#if (!defined(FJCORE_HAVE_EXECINFO_H)) || defined(__FJCORE__)
#endif
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
class Error {
public:
  Error() {}
  Error(const std::string & message);
  virtual ~Error() {}
  std::string message() const {return _message;}
  static void set_print_errors(bool print_errors) {_print_errors = print_errors;}
  static void set_print_backtrace(bool enabled);
  static void set_default_stream(std::ostream * ostr) {
    _default_ostr = ostr;
  }
private:
  std::string _message;                ///< error message
  static bool _print_errors;           ///< do we print anything?
  static bool _print_backtrace;        ///< do we print the backtrace?
  static std::ostream * _default_ostr; ///< the output stream (cerr if not set)
#if (!defined(FJCORE_HAVE_EXECINFO_H)) || defined(__FJCORE__)
  static LimitedWarning _execinfo_undefined;
#endif
};
class InternalError : public Error{
public:
  InternalError(const std::string & message_in) : Error(std::string("*** CRITICAL INTERNAL FASTJET ERROR *** CONTACT THE AUTHORS *** ") + message_in){ }
};
FJCORE_END_NAMESPACE
#endif // __FJCORE_ERROR_HH__
#ifndef __FJCORE_PSEUDOJET_STRUCTURE_BASE_HH__
#define __FJCORE_PSEUDOJET_STRUCTURE_BASE_HH__
#include <vector>
#include <string>
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
class PseudoJet;
class ClusterSequence;
class PseudoJetStructureBase{
public:
  PseudoJetStructureBase(){};
  virtual ~PseudoJetStructureBase(){};
  virtual std::string description() const{ return "PseudoJet with an unknown structure"; }
  virtual bool has_associated_cluster_sequence() const { return false;}
  virtual const ClusterSequence* associated_cluster_sequence() const;
  virtual bool has_valid_cluster_sequence() const {return false;}
  virtual const ClusterSequence * validated_cs() const;
  virtual bool has_partner(const PseudoJet &reference, PseudoJet &partner) const;
  virtual bool has_child(const PseudoJet &reference, PseudoJet &child) const;
  virtual bool has_parents(const PseudoJet &reference, PseudoJet &parent1, PseudoJet &parent2) const;
  virtual bool object_in_jet(const PseudoJet &reference, const PseudoJet &jet) const;
  virtual bool has_constituents() const {return false;}
  virtual std::vector<PseudoJet> constituents(const PseudoJet &reference) const;
  virtual bool has_exclusive_subjets() const {return false;}
  virtual std::vector<PseudoJet> exclusive_subjets(const PseudoJet &reference, const double & dcut) const;
  virtual int n_exclusive_subjets(const PseudoJet &reference, const double & dcut) const;
  virtual std::vector<PseudoJet> exclusive_subjets_up_to (const PseudoJet &reference, int nsub) const;
  virtual double exclusive_subdmerge(const PseudoJet &reference, int nsub) const;
  virtual double exclusive_subdmerge_max(const PseudoJet &reference, int nsub) const;
  virtual bool has_pieces(const PseudoJet & /* reference */) const {
    return false;}
  virtual std::vector<PseudoJet> pieces(const PseudoJet & /* reference */
                                        ) const;
};
FJCORE_END_NAMESPACE
#endif  //  __FJCORE_PSEUDOJET_STRUCTURE_BASE_HH__
#ifndef __FJCORE_PSEUDOJET_HH__
#define __FJCORE_PSEUDOJET_HH__
#include<valarray>
#include<vector>
#include<cassert>
#include<cmath>
#include<iostream>
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
const double MaxRap = 1e5;
const double pseudojet_invalid_phi = -100.0;
const double pseudojet_invalid_rap = -1e200;
class PseudoJet {
 public:
  PseudoJet() : _px(0), _py(0), _pz(0), _E(0) {_finish_init(); _reset_indices();}
  PseudoJet(const double px, const double py, const double pz, const double E);
  template <class L> PseudoJet(const L & some_four_vector);
  PseudoJet(bool /* dummy */) {}
  virtual ~PseudoJet(){};
  inline double E()   const {return _E;}
  inline double e()   const {return _E;} // like CLHEP
  inline double px()  const {return _px;}
  inline double py()  const {return _py;}
  inline double pz()  const {return _pz;}
  inline double phi() const {return phi_02pi();}
  inline double phi_std()  const {
    _ensure_valid_rap_phi();
    return _phi > pi ? _phi-twopi : _phi;}
  inline double phi_02pi() const {
    _ensure_valid_rap_phi();
    return _phi;
  }
  inline double rap() const {
    _ensure_valid_rap_phi();
    return _rap;
  }
  inline double rapidity() const {return rap();} // like CLHEP
  double pseudorapidity() const;
  double eta() const {return pseudorapidity();}
  inline double pt2() const {return _kt2;}
  inline double  pt() const {return sqrt(_kt2);} 
  inline double perp2() const {return _kt2;}  // like CLHEP
  inline double  perp() const {return sqrt(_kt2);}    // like CLHEP
  inline double kt2() const {return _kt2;} // for bkwds compatibility
  inline double  m2() const {return (_E+_pz)*(_E-_pz)-_kt2;}    
  inline double  m() const;    
  inline double mperp2() const {return (_E+_pz)*(_E-_pz);}
  inline double mperp() const {return sqrt(std::abs(mperp2()));}
  inline double mt2() const {return (_E+_pz)*(_E-_pz);}
  inline double mt() const {return sqrt(std::abs(mperp2()));}
  inline double modp2() const {return _kt2+_pz*_pz;}
  inline double modp() const {return sqrt(_kt2+_pz*_pz);}
  inline double Et() const {return (_kt2==0) ? 0.0 : _E/sqrt(1.0+_pz*_pz/_kt2);}
  inline double Et2() const {return (_kt2==0) ? 0.0 : _E*_E/(1.0+_pz*_pz/_kt2);}
  double operator () (int i) const ; 
  inline double operator [] (int i) const { return (*this)(i); }; // this too
  double kt_distance(const PseudoJet & other) const;
  double plain_distance(const PseudoJet & other) const;
  inline double squared_distance(const PseudoJet & other) const {
    return plain_distance(other);}
  inline double delta_R(const PseudoJet & other) const {
    return sqrt(squared_distance(other));
  }
  double delta_phi_to(const PseudoJet & other) const;
  inline double beam_distance() const {return _kt2;}
  std::valarray<double> four_mom() const;
  enum { X=0, Y=1, Z=2, T=3, NUM_COORDINATES=4, SIZE=NUM_COORDINATES };
  PseudoJet & boost(const PseudoJet & prest);
  PseudoJet & unboost(const PseudoJet & prest);
  void operator*=(double);
  void operator/=(double);
  void operator+=(const PseudoJet &);
  void operator-=(const PseudoJet &);
  inline void reset(double px, double py, double pz, double E);
  inline void reset(const PseudoJet & psjet) {
    (*this) = psjet;
  }
  template <class L> inline void reset(const L & some_four_vector) {
    const PseudoJet * pj = fjcore::cast_if_derived<const PseudoJet>(&some_four_vector);
    if (pj){
      (*this) = *pj;
    } else {
      reset(some_four_vector[0], some_four_vector[1],
	    some_four_vector[2], some_four_vector[3]);
    }
  }
  inline void reset_PtYPhiM(double pt_in, double y_in, double phi_in, double m_in=0.0) {
    reset_momentum_PtYPhiM(pt_in, y_in, phi_in, m_in);
    _reset_indices();
  }
  inline void reset_momentum(double px, double py, double pz, double E);
  inline void reset_momentum(const PseudoJet & pj);
  void reset_momentum_PtYPhiM(double pt, double y, double phi, double m=0.0);
  template <class L> inline void reset_momentum(const L & some_four_vector) {
    reset_momentum(some_four_vector[0], some_four_vector[1],
		   some_four_vector[2], some_four_vector[3]);
  }
  void set_cached_rap_phi(double rap, double phi);
  inline int user_index() const {return _user_index;}
  inline void set_user_index(const int index) {_user_index = index;}
  class UserInfoBase{
  public:
    UserInfoBase(){};
    virtual ~UserInfoBase(){}; 
  };
  class InexistentUserInfo : public Error {
  public:
    InexistentUserInfo();
  };
  void set_user_info(UserInfoBase * user_info_in) {
    _user_info.reset(user_info_in);
  }
  template<class L>
  const L & user_info() const{
    if (_user_info.get() == 0) throw InexistentUserInfo();
    return dynamic_cast<const L &>(* _user_info.get());
  }
  bool has_user_info() const{
    return _user_info.get();
  }
  template<class L>
  bool has_user_info() const{
    return _user_info.get() && dynamic_cast<const L *>(_user_info.get());
  }
  const UserInfoBase * user_info_ptr() const{
    return _user_info.get();
  }
  const SharedPtr<UserInfoBase> & user_info_shared_ptr() const{
    return _user_info;
  }
  SharedPtr<UserInfoBase> & user_info_shared_ptr(){
    return _user_info;
  }
  std::string description() const;
  bool has_associated_cluster_sequence() const;
  bool has_associated_cs() const {return has_associated_cluster_sequence();}
  bool has_valid_cluster_sequence() const;
  bool has_valid_cs() const {return has_valid_cluster_sequence();}
  const ClusterSequence* associated_cluster_sequence() const;
  const ClusterSequence* associated_cs() const {return associated_cluster_sequence();}
  inline const ClusterSequence * validated_cluster_sequence() const {
    return validated_cs();
  }
  const ClusterSequence * validated_cs() const;
  void set_structure_shared_ptr(const SharedPtr<PseudoJetStructureBase> &structure);
  bool has_structure() const;
  const PseudoJetStructureBase* structure_ptr() const;
  PseudoJetStructureBase* structure_non_const_ptr();
  const PseudoJetStructureBase* validated_structure_ptr() const;
  const SharedPtr<PseudoJetStructureBase> & structure_shared_ptr() const;
  template<typename StructureType>
  const StructureType & structure() const;
  template<typename TransformerType>
  bool has_structure_of() const;
  template<typename TransformerType>
  const typename TransformerType::StructureType & structure_of() const;
  virtual bool has_partner(PseudoJet &partner) const;
  virtual bool has_child(PseudoJet &child) const;
  virtual bool has_parents(PseudoJet &parent1, PseudoJet &parent2) const;
  virtual bool contains(const PseudoJet &constituent) const;
  virtual bool is_inside(const PseudoJet &jet) const;
  virtual bool has_constituents() const;
  virtual std::vector<PseudoJet> constituents() const;
  virtual bool has_exclusive_subjets() const;
  std::vector<PseudoJet> exclusive_subjets (const double dcut) const;
  int n_exclusive_subjets(const double dcut) const;
  std::vector<PseudoJet> exclusive_subjets (int nsub) const;
  std::vector<PseudoJet> exclusive_subjets_up_to (int nsub) const;
  double exclusive_subdmerge(int nsub) const;
  double exclusive_subdmerge_max(int nsub) const;
  virtual bool has_pieces() const;
  virtual std::vector<PseudoJet> pieces() const;
  inline int cluster_hist_index() const {return _cluster_hist_index;}
  inline void set_cluster_hist_index(const int index) {_cluster_hist_index = index;}
  inline int cluster_sequence_history_index() const {
    return cluster_hist_index();}
  inline void set_cluster_sequence_history_index(const int index) {
    set_cluster_hist_index(index);}
 protected:  
  SharedPtr<PseudoJetStructureBase> _structure;
  SharedPtr<UserInfoBase> _user_info;
 private: 
  double _px,_py,_pz,_E;
  mutable double _phi, _rap;
  double _kt2; 
  int    _cluster_hist_index, _user_index;
  void _finish_init();
  void _reset_indices();
  inline void _ensure_valid_rap_phi() const {
    if (_phi == pseudojet_invalid_phi) _set_rap_phi();
  }
  void _set_rap_phi() const;
  friend PseudoJet operator*(double, const PseudoJet &);
};
PseudoJet operator+(const PseudoJet &, const PseudoJet &);
PseudoJet operator-(const PseudoJet &, const PseudoJet &);
PseudoJet operator*(double, const PseudoJet &);
PseudoJet operator*(const PseudoJet &, double);
PseudoJet operator/(const PseudoJet &, double);
bool operator==(const PseudoJet &, const PseudoJet &);
inline bool operator!=(const PseudoJet & a, const PseudoJet & b) {return !(a==b);}
bool operator==(const PseudoJet & jet, const double val);
inline bool operator==(const double val, const PseudoJet & jet) {return jet == val;}
inline bool operator!=(const PseudoJet & a, const double val)  {return !(a==val);}
inline bool operator!=( const double val, const PseudoJet & a) {return !(a==val);}
inline double dot_product(const PseudoJet & a, const PseudoJet & b) {
  return a.E()*b.E() - a.px()*b.px() - a.py()*b.py() - a.pz()*b.pz();
}
bool have_same_momentum(const PseudoJet &, const PseudoJet &);
PseudoJet PtYPhiM(double pt, double y, double phi, double m = 0.0);
std::vector<PseudoJet> sorted_by_pt(const std::vector<PseudoJet> & jets);
std::vector<PseudoJet> sorted_by_rapidity(const std::vector<PseudoJet> & jets);
std::vector<PseudoJet> sorted_by_E(const std::vector<PseudoJet> & jets);
std::vector<PseudoJet> sorted_by_pz(const std::vector<PseudoJet> & jets);
void sort_indices(std::vector<int> & indices, 
		  const std::vector<double> & values);
template<class T> std::vector<T> objects_sorted_by_values(const std::vector<T> & objects, 
					      const std::vector<double> & values) {
  if (objects.size() != values.size()){
    throw Error("fjcore::objects_sorted_by_values(...): the size of the 'objects' vector must match the size of the 'values' vector");
  }
  std::vector<int> indices(values.size());
  for (size_t i = 0; i < indices.size(); i++) {indices[i] = i;}
  sort_indices(indices, values);
  std::vector<T> objects_sorted(objects.size());
  for (size_t i = 0; i < indices.size(); i++) {
    objects_sorted[i] = objects[indices[i]];
  }
  return objects_sorted;
}
class IndexedSortHelper {
public:
  inline IndexedSortHelper (const std::vector<double> * reference_values) {
    _ref_values = reference_values;
  };
  inline int operator() (const int i1, const int i2) const {
    return  (*_ref_values)[i1] < (*_ref_values)[i2];
  };
private:
  const std::vector<double> * _ref_values;
};
template <class L> inline  PseudoJet::PseudoJet(const L & some_four_vector) {
  reset(some_four_vector);
}
inline void PseudoJet::_reset_indices() { 
  set_cluster_hist_index(-1);
  set_user_index(-1);
  _structure.reset();
  _user_info.reset();
}
inline double PseudoJet::m() const {
  double mm = m2();
  return mm < 0.0 ? -std::sqrt(-mm) : std::sqrt(mm);
}
inline void PseudoJet::reset(double px_in, double py_in, double pz_in, double E_in) {
  _px = px_in;
  _py = py_in;
  _pz = pz_in;
  _E  = E_in;
  _finish_init();
  _reset_indices();
}
inline void PseudoJet::reset_momentum(double px_in, double py_in, double pz_in, double E_in) {
  _px = px_in;
  _py = py_in;
  _pz = pz_in;
  _E  = E_in;
  _finish_init();
}
inline void PseudoJet::reset_momentum(const PseudoJet & pj) {
  _px  = pj._px ;
  _py  = pj._py ;
  _pz  = pj._pz ;
  _E   = pj._E  ;
  _phi = pj._phi;
  _rap = pj._rap;
  _kt2 = pj._kt2;
}
template<typename StructureType>
const StructureType & PseudoJet::structure() const{
  return dynamic_cast<const StructureType &>(* validated_structure_ptr());
}
template<typename TransformerType>
bool PseudoJet::has_structure_of() const{
  if (!_structure) return false;
  return dynamic_cast<const typename TransformerType::StructureType *>(_structure.get()) != 0;
}
template<typename TransformerType>
const typename TransformerType::StructureType & PseudoJet::structure_of() const{
  if (!_structure) 
    throw Error("Trying to access the structure of a PseudoJet without an associated structure");
  return dynamic_cast<const typename TransformerType::StructureType &>(*_structure);
}
PseudoJet join(const std::vector<PseudoJet> & pieces);
PseudoJet join(const PseudoJet & j1);
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2);
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, const PseudoJet & j3);
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, const PseudoJet & j3, const PseudoJet & j4);
FJCORE_END_NAMESPACE
#endif // __FJCORE_PSEUDOJET_HH__
#ifndef __FJCORE_FUNCTION_OF_PSEUDOJET_HH__
#define __FJCORE_FUNCTION_OF_PSEUDOJET_HH__
FJCORE_BEGIN_NAMESPACE
template<typename TOut>
class FunctionOfPseudoJet{
public:
  FunctionOfPseudoJet(){}
  virtual ~FunctionOfPseudoJet(){}
  virtual std::string description() const{ return "";}
  virtual TOut result(const PseudoJet &pj) const = 0;
  TOut operator()(const PseudoJet &pj) const { return result(pj);}
  std::vector<TOut> operator()(const std::vector<PseudoJet> &pjs) const {
    std::vector<TOut> res(pjs.size());
    for (unsigned int i=0; i<pjs.size(); i++)
      res[i] = result(pjs[i]);
    return res;
  }
};
FJCORE_END_NAMESPACE
#endif  // __FJCORE_FUNCTION_OF_PSEUDOJET_HH__
#ifndef __FJCORE_SELECTOR_HH__
#define __FJCORE_SELECTOR_HH__
#include <limits>
#include <cmath>
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
class Selector;
class SelectorWorker {
public:
  virtual ~SelectorWorker() {}
  virtual bool pass(const PseudoJet & jet) const = 0;
  virtual void terminator(std::vector<const PseudoJet *> & jets) const {
    for (unsigned i = 0; i < jets.size(); i++) {
      if (jets[i] && !pass(*jets[i])) jets[i] = NULL;
    }
  }
  virtual bool applies_jet_by_jet() const {return true;}
  virtual std::string description() const {return "missing description";}
  virtual bool takes_reference() const { return false;}
  virtual void set_reference(const PseudoJet & /*reference*/){
    throw Error("set_reference(...) cannot be used for a selector worker that does not take a reference");
  }
  virtual SelectorWorker* copy(){ 
    throw Error("this SelectorWorker has nothing to copy");
  }
  virtual void get_rapidity_extent(double & rapmin, double & rapmax) const {
    rapmax = std::numeric_limits<double>::infinity();
    rapmin = -rapmax; 
  }
  virtual bool is_geometric() const { return false;}
  virtual bool has_finite_area() const;
  virtual bool has_known_area() const { return false;}
  virtual double known_area() const{
    throw Error("this selector has no computable area");
  }
};
class Selector{
public:
  Selector() {}
  Selector(SelectorWorker * worker_in) {_worker.reset(worker_in);}
  virtual ~Selector(){}
  bool pass(const PseudoJet & jet) const {
    if (!validated_worker()->applies_jet_by_jet()) {
      throw Error("Cannot apply this selector to an individual jet");
    }
    return _worker->pass(jet);
  }
  bool operator()(const PseudoJet & jet) const {
    return pass(jet);
  }
  unsigned int count(const std::vector<PseudoJet> & jets) const;
  PseudoJet sum(const std::vector<PseudoJet> & jets) const;
  double scalar_pt_sum(const std::vector<PseudoJet> & jets) const;
  void sift(const std::vector<PseudoJet> & jets,
		  std::vector<PseudoJet> & jets_that_pass,
		  std::vector<PseudoJet> & jets_that_fail) const;
  bool applies_jet_by_jet() const {
    return validated_worker()->applies_jet_by_jet();
  }
  std::vector<PseudoJet> operator()(const std::vector<PseudoJet> & jets) const;
  virtual void nullify_non_selected(std::vector<const PseudoJet *> & jets) const {
    validated_worker()->terminator(jets);
  }
  void get_rapidity_extent(double &rapmin, double &rapmax) const {
    return validated_worker()->get_rapidity_extent(rapmin, rapmax);
  }
  std::string description() const {
    return validated_worker()->description();
  }
  bool is_geometric() const{
    return validated_worker()->is_geometric();
  }
  bool has_finite_area() const{
    return validated_worker()->has_finite_area();
  }
  const SharedPtr<SelectorWorker> & worker() const {return _worker;}
  const SelectorWorker* validated_worker() const {
    const SelectorWorker* worker_ptr = _worker.get();
    if (worker_ptr == 0) throw InvalidWorker();
    return worker_ptr;
  }
  bool takes_reference() const {
    return validated_worker()->takes_reference();
  }
  const Selector & set_reference(const PseudoJet &reference){
    if (! validated_worker()->takes_reference()){
      return *this;
    }
    _copy_worker_if_needed();
    _worker->set_reference(reference);
    return *this;
  }
  class InvalidWorker : public Error {
  public:
    InvalidWorker() : Error("Attempt to use Selector with no valid underlying worker") {}
  };
  class InvalidArea : public Error {
  public:
    InvalidArea() : Error("Attempt to obtain area from Selector for which this is not meaningful") {}
  };
  Selector & operator &=(const Selector & b);
  Selector & operator |=(const Selector & b);
protected:
  void _copy_worker_if_needed(){
    if (_worker.unique()) return;
    _worker.reset(_worker->copy());
  }
private:
  SharedPtr<SelectorWorker> _worker; ///< the underlying worker
};
Selector SelectorIdentity();
Selector operator!(const Selector & s);
Selector operator ||(const Selector & s1, const Selector & s2);
Selector operator&&(const Selector & s1, const Selector & s2);
Selector operator*(const Selector & s1, const Selector & s2);
Selector SelectorPtMin(double ptmin);                    ///< select objects with pt >= ptmin
Selector SelectorPtMax(double ptmax);                    ///< select objects with pt <= ptmax
Selector SelectorPtRange(double ptmin, double ptmax);    ///< select objects with ptmin <= pt <= ptmax
Selector SelectorEtMin(double Etmin);                    ///< select objects with Et >= Etmin
Selector SelectorEtMax(double Etmax);                    ///< select objects with Et <= Etmax
Selector SelectorEtRange(double Etmin, double Etmax);    ///< select objects with Etmin <= Et <= Etmax
Selector SelectorEMin(double Emin);                      ///< select objects with E >= Emin
Selector SelectorEMax(double Emax);                      ///< select objects with E <= Emax
Selector SelectorERange(double Emin, double Emax);       ///< select objects with Emin <= E <= Emax
Selector SelectorMassMin(double Mmin);                      ///< select objects with Mass >= Mmin
Selector SelectorMassMax(double Mmax);                      ///< select objects with Mass <= Mmax
Selector SelectorMassRange(double Mmin, double Mmax);       ///< select objects with Mmin <= Mass <= Mmax
Selector SelectorRapMin(double rapmin);                  ///< select objects with rap >= rapmin
Selector SelectorRapMax(double rapmax);                  ///< select objects with rap <= rapmax
Selector SelectorRapRange(double rapmin, double rapmax); ///< select objects with rapmin <= rap <= rapmax
Selector SelectorAbsRapMin(double absrapmin);                     ///< select objects with |rap| >= absrapmin
Selector SelectorAbsRapMax(double absrapmax);                     ///< select objects with |rap| <= absrapmax
Selector SelectorAbsRapRange(double absrapmin, double absrapmax); ///< select objects with absrapmin <= |rap| <= absrapmax
Selector SelectorEtaMin(double etamin);                  ///< select objects with eta >= etamin
Selector SelectorEtaMax(double etamax);                  ///< select objects with eta <= etamax
Selector SelectorEtaRange(double etamin, double etamax); ///< select objects with etamin <= eta <= etamax
Selector SelectorAbsEtaMin(double absetamin);                     ///< select objects with |eta| >= absetamin
Selector SelectorAbsEtaMax(double absetamax);                     ///< select objects with |eta| <= absetamax
Selector SelectorAbsEtaRange(double absetamin, double absetamax); ///< select objects with absetamin <= |eta| <= absetamax
Selector SelectorPhiRange(double phimin, double phimax); ///< select objects with phimin <= phi <= phimax
Selector SelectorRapPhiRange(double rapmin, double rapmax, double phimin, double phimax);
Selector SelectorNHardest(unsigned int n); 
Selector SelectorCircle(const double radius); 
Selector SelectorDoughnut(const double radius_in, const double radius_out); 
Selector SelectorStrip(const double half_width);
Selector SelectorRectangle(const double half_rap_width, const double half_phi_width);
Selector SelectorPtFractionMin(double fraction);
Selector SelectorIsZero();
FJCORE_END_NAMESPACE      // defined in fastjet/internal/base.hh
#endif // __FJCORE_SELECTOR_HH__
#ifndef __FJCORE_JETDEFINITION_HH__
#define __FJCORE_JETDEFINITION_HH__
#include<cassert>
#include<string>
#include<memory>
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
std::string fastjet_version_string();
enum Strategy {
  N2MHTLazy9AntiKtSeparateGhosts   = -10, 
  N2MHTLazy9   = -7, 
  N2MHTLazy25   = -6, 
  N2MHTLazy9Alt   = -5, 
  N2MinHeapTiled   = -4, 
  N2Tiled     = -3, 
  N2PoorTiled = -2, 
  N2Plain     = -1, 
  N3Dumb      =  0, 
  Best        =  1, 
  NlnN        =  2, 
  NlnN3pi     =  3, 
  NlnN4pi     =  4,
  NlnNCam4pi   = 14,
  NlnNCam2pi2R = 13,
  NlnNCam      = 12, // 2piMultD
  BestFJ30     =  21, 
  plugin_strategy = 999
};
enum JetAlgorithm {
  kt_algorithm=0,
  cambridge_algorithm=1,
  antikt_algorithm=2, 
  genkt_algorithm=3, 
  cambridge_for_passive_algorithm=11,
  genkt_for_passive_algorithm=13, 
  ee_kt_algorithm=50,
  ee_genkt_algorithm=53,
  plugin_algorithm = 99,
  undefined_jet_algorithm = 999
};
typedef JetAlgorithm JetFinder;
const JetAlgorithm aachen_algorithm = cambridge_algorithm;
const JetAlgorithm cambridge_aachen_algorithm = cambridge_algorithm;
enum RecombinationScheme {
  E_scheme=0,
  pt_scheme=1,
  pt2_scheme=2,
  Et_scheme=3,
  Et2_scheme=4,
  BIpt_scheme=5,
  BIpt2_scheme=6,
  WTA_pt_scheme=7,
  WTA_modp_scheme=8,
  external_scheme = 99
};
class ClusterSequence;
class JetDefinition {
public:
  class Plugin;
  class Recombiner;
  JetDefinition(JetAlgorithm jet_algorithm_in, 
                double R_in, 
                RecombinationScheme recomb_scheme_in = E_scheme,
                Strategy strategy_in = Best) {
    *this = JetDefinition(jet_algorithm_in, R_in, recomb_scheme_in, strategy_in, 1);
  }
  JetDefinition(JetAlgorithm jet_algorithm_in, 
                RecombinationScheme recomb_scheme_in = E_scheme,
                Strategy strategy_in = Best) {
    double dummyR = 0.0;
    *this = JetDefinition(jet_algorithm_in, dummyR, recomb_scheme_in, strategy_in, 0);
  }
  JetDefinition(JetAlgorithm jet_algorithm_in, 
                double R_in, 
                double xtra_param_in,
                RecombinationScheme recomb_scheme_in = E_scheme,
                Strategy strategy_in = Best) {
    *this = JetDefinition(jet_algorithm_in, R_in, recomb_scheme_in, strategy_in, 2);
    set_extra_param(xtra_param_in);
  }
  JetDefinition(JetAlgorithm jet_algorithm_in, 
                double R_in, 
                const Recombiner * recombiner_in,
                Strategy strategy_in = Best) {
    *this = JetDefinition(jet_algorithm_in, R_in, external_scheme, strategy_in);
    _recombiner = recombiner_in;
  }
  JetDefinition(JetAlgorithm jet_algorithm_in, 
                const Recombiner * recombiner_in,
                Strategy strategy_in = Best) {
    *this = JetDefinition(jet_algorithm_in, external_scheme, strategy_in);
    _recombiner = recombiner_in;
  }
  JetDefinition(JetAlgorithm jet_algorithm_in, 
                double R_in, 
                double xtra_param_in,
                const Recombiner * recombiner_in,
                Strategy strategy_in = Best) {
    *this = JetDefinition(jet_algorithm_in, R_in, xtra_param_in, external_scheme, strategy_in);
    _recombiner = recombiner_in;
  }
  JetDefinition()  {
    *this = JetDefinition(undefined_jet_algorithm, 1.0);
  }
  JetDefinition(const Plugin * plugin_in) {
    _plugin = plugin_in;
    _strategy = plugin_strategy;
    _Rparam = _plugin->R();
    _extra_param = 0.0; // a dummy value to keep static code checkers happy
    _jet_algorithm = plugin_algorithm;
    set_recombination_scheme(E_scheme);
  }
  JetDefinition(JetAlgorithm jet_algorithm_in, 
                double R_in, 
                RecombinationScheme recomb_scheme_in,
                Strategy strategy_in,
                int nparameters_in);
  FJCORE_DEPRECATED_MSG("This argument ordering is deprecated. Use JetDefinition(alg, R, strategy, scheme[, n_parameters]) instead")
  JetDefinition(JetAlgorithm jet_algorithm_in, 
                double R_in, 
                Strategy strategy_in,
                RecombinationScheme recomb_scheme_in = E_scheme,
                int nparameters_in = 1){
    (*this) = JetDefinition(jet_algorithm_in,R_in,recomb_scheme_in,strategy_in,nparameters_in);
  }
  template <class L> 
  std::vector<PseudoJet> operator()(const std::vector<L> & particles) const;
  static const double max_allowable_R; //= 1000.0;
  void set_recombination_scheme(RecombinationScheme);
  void set_recombiner(const Recombiner * recomb) {
    if (_shared_recombiner) _shared_recombiner.reset(recomb);
    _recombiner = recomb;
    _default_recombiner = DefaultRecombiner(external_scheme);
  }
  void set_recombiner(const JetDefinition &other_jet_def);
  void delete_recombiner_when_unused();
  const Plugin * plugin() const {return _plugin;};
  void delete_plugin_when_unused();
  JetAlgorithm jet_algorithm  () const {return _jet_algorithm  ;}
  JetAlgorithm jet_finder     () const {return _jet_algorithm  ;}
  double    R           () const {return _Rparam      ;}
  double    extra_param () const {return _extra_param ;}
  Strategy  strategy    () const {return _strategy    ;}
  RecombinationScheme recombination_scheme() const {
    return _default_recombiner.scheme();}
  void set_jet_algorithm(JetAlgorithm njf) {_jet_algorithm = njf;}
  void set_jet_finder(JetAlgorithm njf)    {_jet_algorithm = njf;}
  void set_extra_param(double xtra_param) {_extra_param = xtra_param;}
  const Recombiner * recombiner() const {
    return _recombiner == 0 ? & _default_recombiner : _recombiner;}
  bool has_same_recombiner(const JetDefinition &other_jd) const;
  bool is_spherical() const;
  std::string description() const;
  std::string description_no_recombiner() const;
  static std::string algorithm_description(const JetAlgorithm jet_alg);
  static unsigned int n_parameters_for_algorithm(const JetAlgorithm jet_alg);
public:
  class Recombiner {
  public:
    virtual std::string description() const = 0;
    virtual void recombine(const PseudoJet & pa, const PseudoJet & pb, 
                           PseudoJet & pab) const = 0;
    virtual void preprocess(PseudoJet & ) const {};
    virtual ~Recombiner() {};
    inline void plus_equal(PseudoJet & pa, const PseudoJet & pb) const {
      PseudoJet pres; 
      recombine(pa,pb,pres);
      pa = pres;
    }
  };
  class DefaultRecombiner : public Recombiner {
  public:
    DefaultRecombiner(RecombinationScheme recomb_scheme = E_scheme) : 
      _recomb_scheme(recomb_scheme) {}
    virtual std::string description() const FJCORE_OVERRIDE;
    virtual void recombine(const PseudoJet & pa, const PseudoJet & pb, 
                           PseudoJet & pab) const FJCORE_OVERRIDE;
    virtual void preprocess(PseudoJet & p) const FJCORE_OVERRIDE;
    RecombinationScheme scheme() const {return _recomb_scheme;}
  private:
    RecombinationScheme _recomb_scheme;
  };
  class Plugin{
  public:
    virtual std::string description() const = 0;
    virtual void run_clustering(ClusterSequence &) const = 0;
    virtual double R() const = 0;
    virtual bool supports_ghosted_passive_areas() const {return false;}
    virtual void set_ghost_separation_scale(double scale) const;
    virtual double ghost_separation_scale() const {return 0.0;}
    virtual bool exclusive_sequence_meaningful() const {return false;}
    virtual bool is_spherical() const {return false;}
    virtual ~Plugin() {};
  };
private:
  JetAlgorithm _jet_algorithm;
  double    _Rparam;
  double    _extra_param ; ///< parameter whose meaning varies according to context
  Strategy  _strategy  ;
  const Plugin * _plugin;
  SharedPtr<const Plugin> _plugin_shared;
  DefaultRecombiner _default_recombiner;
  const Recombiner * _recombiner;
  SharedPtr<const Recombiner> _shared_recombiner;
};
PseudoJet join(const std::vector<PseudoJet> & pieces, const JetDefinition::Recombiner & recombiner);
PseudoJet join(const PseudoJet & j1, 
	       const JetDefinition::Recombiner & recombiner);
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, 
	       const JetDefinition::Recombiner & recombiner);
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, const PseudoJet & j3, 
	       const JetDefinition::Recombiner & recombiner);
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, const PseudoJet & j3, const PseudoJet & j4, 
	       const JetDefinition::Recombiner & recombiner);
FJCORE_END_NAMESPACE
#endif // __FJCORE_JETDEFINITION_HH__
#ifndef __FJCORE_COMPOSITEJET_STRUCTURE_HH__
#define __FJCORE_COMPOSITEJET_STRUCTURE_HH__
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
class CompositeJetStructure : public PseudoJetStructureBase{
public:
  CompositeJetStructure() : _area_4vector_ptr(0){};
  CompositeJetStructure(const std::vector<PseudoJet> & initial_pieces, 
			const JetDefinition::Recombiner * recombiner = 0);
  virtual ~CompositeJetStructure(){
    if (_area_4vector_ptr) delete _area_4vector_ptr;
  };
  virtual std::string description() const FJCORE_OVERRIDE;
  virtual bool has_constituents() const FJCORE_OVERRIDE;
  virtual std::vector<PseudoJet> constituents(const PseudoJet &jet) const FJCORE_OVERRIDE;
  virtual bool has_pieces(const PseudoJet & /*jet*/) const FJCORE_OVERRIDE {return true;}
  virtual std::vector<PseudoJet> pieces(const PseudoJet &jet) const FJCORE_OVERRIDE;
protected:
  std::vector<PseudoJet> _pieces;  ///< the pieces building the jet
  PseudoJet * _area_4vector_ptr;   ///< pointer to the 4-vector jet area
};
template<typename T> PseudoJet join(const std::vector<PseudoJet> & pieces){
  PseudoJet result(0.0,0.0,0.0,0.0);
  for (unsigned int i=0; i<pieces.size(); i++){
    const PseudoJet it = pieces[i];
    result += it;
  }
  T *cj_struct = new T(pieces);
  result.set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>(cj_struct));
  return result;
}
template<typename T> PseudoJet join(const PseudoJet & j1){
  return join<T>(std::vector<PseudoJet>(1,j1));
}
template<typename T> PseudoJet join(const PseudoJet & j1, const PseudoJet & j2){
  std::vector<PseudoJet> pieces;
  pieces.push_back(j1);
  pieces.push_back(j2);
  return join<T>(pieces);
}
template<typename T> PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, 
				    const PseudoJet & j3){
  std::vector<PseudoJet> pieces;
  pieces.push_back(j1);
  pieces.push_back(j2);
  pieces.push_back(j3);
  return join<T>(pieces);
}
template<typename T> PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, 
				    const PseudoJet & j3, const PseudoJet & j4){
  std::vector<PseudoJet> pieces;
  pieces.push_back(j1);
  pieces.push_back(j2);
  pieces.push_back(j3);
  pieces.push_back(j4);
  return join<T>(pieces);
}
template<typename T> PseudoJet join(const std::vector<PseudoJet> & pieces, 
				    const JetDefinition::Recombiner & recombiner){
  PseudoJet result;
  if (pieces.size()>0){
    result = pieces[0];
    for (unsigned int i=1; i<pieces.size(); i++){
      recombiner.plus_equal(result, pieces[i]);
    }
  }
  T *cj_struct = new T(pieces, &recombiner);
  result.set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>(cj_struct));
  return result;
}
template<typename T> PseudoJet join(const PseudoJet & j1, 
				    const JetDefinition::Recombiner & recombiner){
  return join<T>(std::vector<PseudoJet>(1,j1), recombiner);
}
template<typename T> PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, 
				    const JetDefinition::Recombiner & recombiner){
  std::vector<PseudoJet> pieces;
  pieces.reserve(2);
  pieces.push_back(j1);
  pieces.push_back(j2);
  return join<T>(pieces, recombiner);
}
template<typename T> PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, 
				    const PseudoJet & j3, 
				    const JetDefinition::Recombiner & recombiner){
  std::vector<PseudoJet> pieces;
  pieces.reserve(3);
  pieces.push_back(j1);
  pieces.push_back(j2);
  pieces.push_back(j3);
  return join<T>(pieces, recombiner);
}
template<typename T> PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, 
				    const PseudoJet & j3, const PseudoJet & j4, 
				    const JetDefinition::Recombiner & recombiner){
  std::vector<PseudoJet> pieces;
  pieces.reserve(4);
  pieces.push_back(j1);
  pieces.push_back(j2);
  pieces.push_back(j3);
  pieces.push_back(j4);
  return join<T>(pieces, recombiner);
}
FJCORE_END_NAMESPACE      // defined in fastjet/internal/base.hh
#endif // __FJCORE_MERGEDJET_STRUCTURE_HH__
#ifndef __FJCORE_CLUSTER_SEQUENCE_STRUCTURE_HH__
#define __FJCORE_CLUSTER_SEQUENCE_STRUCTURE_HH__
#include <vector>
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
class ClusterSequenceStructure : public PseudoJetStructureBase{
public:
  ClusterSequenceStructure() : _associated_cs(NULL){}
  ClusterSequenceStructure(const ClusterSequence *cs){
    set_associated_cs(cs);
  };
  virtual ~ClusterSequenceStructure();
  virtual std::string description() const FJCORE_OVERRIDE{
    return "PseudoJet with an associated ClusterSequence";
  }
  virtual bool has_associated_cluster_sequence() const FJCORE_OVERRIDE{ return true;}
  virtual const ClusterSequence* associated_cluster_sequence() const FJCORE_OVERRIDE;
  virtual bool has_valid_cluster_sequence() const FJCORE_OVERRIDE;
  virtual const ClusterSequence * validated_cs() const FJCORE_OVERRIDE;
  virtual void set_associated_cs(const ClusterSequence * new_cs){
    _associated_cs = new_cs;
  }
  virtual bool has_partner(const PseudoJet &reference, PseudoJet &partner) const FJCORE_OVERRIDE;
  virtual bool has_child(const PseudoJet &reference, PseudoJet &child) const FJCORE_OVERRIDE;
  virtual bool has_parents(const PseudoJet &reference, PseudoJet &parent1, PseudoJet &parent2) const FJCORE_OVERRIDE;
  virtual bool object_in_jet(const PseudoJet &reference, const PseudoJet &jet) const FJCORE_OVERRIDE;
  virtual bool has_constituents() const FJCORE_OVERRIDE;
  virtual std::vector<PseudoJet> constituents(const PseudoJet &reference) const FJCORE_OVERRIDE;
  virtual bool has_exclusive_subjets() const FJCORE_OVERRIDE;
  virtual std::vector<PseudoJet> exclusive_subjets(const PseudoJet &reference, const double & dcut) const FJCORE_OVERRIDE;
  virtual int n_exclusive_subjets(const PseudoJet &reference, const double & dcut) const FJCORE_OVERRIDE;
  virtual std::vector<PseudoJet> exclusive_subjets_up_to (const PseudoJet &reference, int nsub) const FJCORE_OVERRIDE;
  virtual double exclusive_subdmerge(const PseudoJet &reference, int nsub) const FJCORE_OVERRIDE;
  virtual double exclusive_subdmerge_max(const PseudoJet &reference, int nsub) const FJCORE_OVERRIDE;
  virtual bool has_pieces(const PseudoJet &reference) const FJCORE_OVERRIDE;
  virtual std::vector<PseudoJet> pieces(const PseudoJet &reference) const FJCORE_OVERRIDE;
protected:
  const ClusterSequence *_associated_cs;
};
FJCORE_END_NAMESPACE
#endif  //  __FJCORE_CLUSTER_SEQUENCE_STRUCTURE_HH__
#ifndef __FJCORE_CLUSTERSEQUENCE_HH__
#define __FJCORE_CLUSTERSEQUENCE_HH__
#include<vector>
#include<map>
#include<memory>
#include<cassert>
#include<iostream>
#include<string>
#include<set>
#include<cmath> // needed to get double std::abs(double)
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
class ClusterSequenceStructure;
class DynamicNearestNeighbours;
class ClusterSequence {
 public: 
  ClusterSequence () : _deletes_self_when_unused(false) {}
  template<class L> ClusterSequence (
			          const std::vector<L> & pseudojets,
				  const JetDefinition & jet_def,
				  const bool & writeout_combinations = false);
  ClusterSequence (const ClusterSequence & cs) : _deletes_self_when_unused(false) {
    transfer_from_sequence(cs);
  }
  ClusterSequence & operator=(const ClusterSequence & cs);
  virtual ~ClusterSequence (); //{}
  std::vector<PseudoJet> inclusive_jets (const double ptmin = 0.0) const;
  int n_exclusive_jets (const double dcut) const;
  std::vector<PseudoJet> exclusive_jets (const double dcut) const;
  std::vector<PseudoJet> exclusive_jets (const int njets) const;
  std::vector<PseudoJet> exclusive_jets_up_to (const int njets) const;
  double exclusive_dmerge (const int njets) const;
  double exclusive_dmerge_max (const int njets) const;
  double exclusive_ymerge (int njets) const {return exclusive_dmerge(njets) / Q2();}
  double exclusive_ymerge_max (int njets) const {return exclusive_dmerge_max(njets)/Q2();}
  int n_exclusive_jets_ycut (double ycut) const {return n_exclusive_jets(ycut*Q2());}
  std::vector<PseudoJet> exclusive_jets_ycut (double ycut) const {
    int njets = n_exclusive_jets_ycut(ycut);
    return exclusive_jets(njets);
  }
  std::vector<PseudoJet> exclusive_subjets (const PseudoJet & jet, 
                                            const double dcut) const;
  int n_exclusive_subjets(const PseudoJet & jet, 
                          const double dcut) const;
  std::vector<PseudoJet> exclusive_subjets (const PseudoJet & jet, 
                                            int nsub) const;
  std::vector<PseudoJet> exclusive_subjets_up_to (const PseudoJet & jet, 
						  int nsub) const;
  double exclusive_subdmerge(const PseudoJet & jet, int nsub) const;
  double exclusive_subdmerge_max(const PseudoJet & jet, int nsub) const;
  double Q() const {return _Qtot;}
  double Q2() const {return _Qtot*_Qtot;}
  bool object_in_jet(const PseudoJet & object, const PseudoJet & jet) const;
  bool has_parents(const PseudoJet & jet, PseudoJet & parent1, 
               PseudoJet & parent2) const;
  bool has_child(const PseudoJet & jet, PseudoJet & child) const;
  bool has_child(const PseudoJet & jet, const PseudoJet * & childp) const;
  bool has_partner(const PseudoJet & jet, PseudoJet & partner) const;
  std::vector<PseudoJet> constituents (const PseudoJet & jet) const;
  void print_jets_for_root(const std::vector<PseudoJet> & jets, 
                           std::ostream & ostr = std::cout) const;
  void print_jets_for_root(const std::vector<PseudoJet> & jets, 
                           const std::string & filename,
			   const std::string & comment = "") const;
  void add_constituents (const PseudoJet & jet, 
			 std::vector<PseudoJet> & subjet_vector) const;
  inline Strategy strategy_used () const {return _strategy;}
  std::string strategy_string () const {return strategy_string(_strategy);}
  std::string strategy_string (Strategy strategy_in) const;
  const JetDefinition & jet_def() const {return _jet_def;}
  void delete_self_when_unused();
  bool will_delete_self_when_unused() const {return _deletes_self_when_unused;}
  void signal_imminent_self_deletion() const;
  double jet_scale_for_algorithm(const PseudoJet & jet) const;
  void plugin_record_ij_recombination(int jet_i, int jet_j, double dij, 
				      int & newjet_k) {
    assert(plugin_activated());
    _do_ij_recombination_step(jet_i, jet_j, dij, newjet_k);
  }
  void plugin_record_ij_recombination(int jet_i, int jet_j, double dij, 
				      const PseudoJet & newjet, 
				      int & newjet_k);
  void plugin_record_iB_recombination(int jet_i, double diB) {
    assert(plugin_activated());
    _do_iB_recombination_step(jet_i, diB);
  }
  class Extras {
  public:
    virtual ~Extras() {}
    virtual std::string description() const {return "This is a dummy extras class that contains no extra information! Derive from it if you want to use it to provide extra information from a plugin jet finder";}
  };
  inline void plugin_associate_extras(Extras * extras_in) {
    _extras.reset(extras_in);
  }
// TS : fully disable use of auto_ptr (already deprecated).
//#ifdef FJCORE_HAVE_AUTO_PTR_INTERFACE
//  FJCORE_DEPRECATED_MSG("Please use ClusterSequence::plugin_associate_extras(Extras * extras_in)) instead")
//  inline void plugin_associate_extras(std::auto_ptr<Extras> extras_in){
//    _extras.reset(extras_in.release());
//  }
//#endif
  inline bool plugin_activated() const {return _plugin_activated;}
  const Extras * extras() const {return _extras.get();}
  template<class GBJ> void plugin_simple_N2_cluster () {
    assert(plugin_activated());
    _simple_N2_cluster<GBJ>();
  }
public:
  struct history_element{
    int parent1; /// index in _history where first parent of this jet
    int parent2; /// index in _history where second parent of this jet
    int child;   /// index in _history where the current jet is
		 /// recombined with another jet to form its child. It
		 /// is Invalid if this jet does not further
		 /// recombine.
    int jetp_index; /// index in the _jets vector where we will find the
    double dij;  /// the distance corresponding to the recombination
		 /// at this stage of the clustering.
    double max_dij_so_far; /// the largest recombination distance seen
			   /// so far in the clustering history.
  };
  enum JetType {Invalid=-3, InexistentParent = -2, BeamJet = -1};
  const std::vector<PseudoJet> & jets()    const;
  const std::vector<history_element> & history() const;
  unsigned int n_particles() const;
  std::vector<int> particle_jet_indices(const std::vector<PseudoJet> &) const;
  std::vector<int> unique_history_order() const;
  std::vector<PseudoJet> unclustered_particles() const;
  std::vector<PseudoJet> childless_pseudojets() const;
  bool contains(const PseudoJet & object) const;
  void transfer_from_sequence(const ClusterSequence & from_seq,
			      const FunctionOfPseudoJet<PseudoJet> * action_on_jets = 0);
  const SharedPtr<PseudoJetStructureBase> & structure_shared_ptr() const{
    return _structure_shared_ptr;
  }
  typedef ClusterSequenceStructure StructureType;
  static void print_banner();
  static void set_fastjet_banner_stream(std::ostream * ostr) {_fastjet_banner_ostr = ostr;}
  static std::ostream * fastjet_banner_stream() {return _fastjet_banner_ostr;}
private:
  static std::ostream * _fastjet_banner_ostr;
protected:
  JetDefinition _jet_def;
  template<class L> void _transfer_input_jets(
                                     const std::vector<L> & pseudojets);
  void _initialise_and_run (const JetDefinition & jet_def,
			    const bool & writeout_combinations);
  void _initialise_and_run_no_decant();
  void _decant_options(const JetDefinition & jet_def,
                       const bool & writeout_combinations);
  void _decant_options_partial();
  void _fill_initial_history();
  void _do_ij_recombination_step(const int jet_i, const int jet_j, 
				 const double dij, int & newjet_k);
  void _do_iB_recombination_step(const int jet_i, const double diB);
  void _set_structure_shared_ptr(PseudoJet & j);
  void _update_structure_use_count();
  Strategy _best_strategy() const;
  class _Parabola {
  public:
    _Parabola(double a, double b, double c) : _a(a), _b(b), _c(c) {}
    inline double operator()(const double R) const {return _c*(_a*R*R + _b*R + 1);}
  private:
    double _a, _b, _c;
  };
  class _Line {
  public:
    _Line(double a, double b) : _a(a), _b(b) {}
    inline double operator()(const double R) const {return _a*R + _b;}
  private:
    double _a, _b;
  };
  std::vector<PseudoJet> _jets;
  std::vector<history_element> _history;
  void get_subhist_set(std::set<const history_element*> & subhist,
                       const  PseudoJet & jet, double dcut, int maxjet) const;
  bool _writeout_combinations;
  int  _initial_n;
  double _Rparam, _R2, _invR2;
  double _Qtot;
  Strategy    _strategy;
  JetAlgorithm  _jet_algorithm;
  SharedPtr<PseudoJetStructureBase> _structure_shared_ptr; //< will actually be of type ClusterSequenceStructure
  int _structure_use_count_after_construction; //< info of use when CS handles its own memory
  mutable bool _deletes_self_when_unused;
 private:
  bool _plugin_activated;
  SharedPtr<Extras> _extras; // things the plugin might want to add
  void _really_dumb_cluster ();
  void _delaunay_cluster ();
  template<class BJ> void _simple_N2_cluster ();
  void _tiled_N2_cluster ();
  void _faster_tiled_N2_cluster ();
  void _minheap_faster_tiled_N2_cluster();
  void _CP2DChan_cluster();
  void _CP2DChan_cluster_2pi2R ();
  void _CP2DChan_cluster_2piMultD ();
  void _CP2DChan_limited_cluster(double D);
  void _do_Cambridge_inclusive_jets();
  void _fast_NsqrtN_cluster();
  void _add_step_to_history( //const int step_number,
                            const int parent1, 
			    const int parent2, const int jetp_index,
			    const double dij);
  void _extract_tree_children(int pos, std::valarray<bool> &, 
		const std::valarray<int> &, std::vector<int> &) const;
  void _extract_tree_parents (int pos, std::valarray<bool> &, 
                const std::valarray<int> &,  std::vector<int> &) const;
  typedef std::pair<int,int> TwoVertices;
  typedef std::pair<double,TwoVertices> DijEntry;
  typedef std::multimap<double,TwoVertices> DistMap;
  void _add_ktdistance_to_map(const int ii, 
			      DistMap & DijMap,
  			      const DynamicNearestNeighbours * DNN);
  static bool _first_time;
  static LimitedWarning _exclusive_warnings;
  static LimitedWarning _changed_strategy_warning;
  struct BriefJet {
    double     eta, phi, kt2, NN_dist;
    BriefJet * NN;
    int        _jets_index;
  };
  class TiledJet {
  public:
    double     eta, phi, kt2, NN_dist;
    TiledJet * NN, *previous, * next; 
    int        _jets_index, tile_index, diJ_posn;
    inline void label_minheap_update_needed() {diJ_posn = 1;}
    inline void label_minheap_update_done()   {diJ_posn = 0;}
    inline bool minheap_update_needed() const {return diJ_posn==1;}
  };
  template <class J> void _bj_set_jetinfo( J * const jet, 
						 const int _jets_index) const;
  void _bj_remove_from_tiles( TiledJet * const jet) const;
  template <class J> double _bj_dist(const J * const jeta, 
			const J * const jetb) const;
  template <class J> double _bj_diJ(const J * const jeta) const;
  template <class J> inline J * _bj_of_hindex(
                          const int hist_index, 
			  J * const head, J * const tail) 
    const {
    J * res;
    for(res = head; res<tail; res++) {
      if (_jets[res->_jets_index].cluster_hist_index() == hist_index) {break;}
    }
    return res;
  }
  template <class J> void _bj_set_NN_nocross(J * const jeta, 
            J * const head, const J * const tail) const;
  template <class J> void _bj_set_NN_crosscheck(J * const jeta, 
            J * const head, const J * const tail) const;
  static const int n_tile_neighbours = 9;
  struct Tile {
    Tile *   begin_tiles[n_tile_neighbours]; 
    Tile **  surrounding_tiles; 
    Tile **  RH_tiles;  
    Tile **  end_tiles; 
    TiledJet * head;    
    bool     tagged;    
  };
  std::vector<Tile> _tiles;
  double _tiles_eta_min, _tiles_eta_max;
  double _tile_size_eta, _tile_size_phi;
  int    _n_tiles_phi,_tiles_ieta_min,_tiles_ieta_max;
  inline int _tile_index (int ieta, int iphi) const {
    return (ieta-_tiles_ieta_min)*_n_tiles_phi
                  + (iphi+_n_tiles_phi) % _n_tiles_phi;
  }
  int  _tile_index(const double eta, const double phi) const;
  void _tj_set_jetinfo ( TiledJet * const jet, const int _jets_index);
  void  _bj_remove_from_tiles(TiledJet * const jet);
  void _initialise_tiles();
  void _print_tiles(TiledJet * briefjets ) const;
  void _add_neighbours_to_tile_union(const int tile_index, 
		 std::vector<int> & tile_union, int & n_near_tiles) const;
  void _add_untagged_neighbours_to_tile_union(const int tile_index, 
		 std::vector<int> & tile_union, int & n_near_tiles);
  struct EEBriefJet {
    double NN_dist;  // obligatorily present
    double kt2;      // obligatorily present == E^2 in general
    EEBriefJet * NN; // must be present too
    int    _jets_index; // must also be present!
    double nx, ny, nz;  // our internal storage for fast distance calcs
  };
  void _simple_N2_cluster_BriefJet();
  void _simple_N2_cluster_EEBriefJet();
};
template<class L> void ClusterSequence::_transfer_input_jets(
                                       const std::vector<L> & pseudojets) {
  _jets.reserve(pseudojets.size()*2);
  for (unsigned int i = 0; i < pseudojets.size(); i++) {
    _jets.push_back(pseudojets[i]);}
}
template<class L> ClusterSequence::ClusterSequence (
			          const std::vector<L> & pseudojets,
				  const JetDefinition & jet_def_in,
				  const bool & writeout_combinations) :
  _jet_def(jet_def_in), _writeout_combinations(writeout_combinations),
  _structure_shared_ptr(new ClusterSequenceStructure(this))
{
  _transfer_input_jets(pseudojets);
  _decant_options_partial();
  _initialise_and_run_no_decant();
}
inline const std::vector<PseudoJet> & ClusterSequence::jets () const {
  return _jets;
}
inline const std::vector<ClusterSequence::history_element> & ClusterSequence::history () const {
  return _history;
}
inline unsigned int ClusterSequence::n_particles() const {return _initial_n;}
#ifndef __CINT__
template<class L>
std::vector<PseudoJet> JetDefinition::operator()(const std::vector<L> & particles) const {
  ClusterSequence * cs = new ClusterSequence(particles, *this);
  std::vector<PseudoJet> jets;
  if (is_spherical()) {
    jets = sorted_by_E(cs->inclusive_jets());
  } else {
    jets = sorted_by_pt(cs->inclusive_jets());
  }
  if (jets.size() != 0) {
    cs->delete_self_when_unused();
  } else {
    delete cs;
  }
  return jets;
}
#endif // __CINT__
template <class J> inline void ClusterSequence::_bj_set_jetinfo(
                            J * const jetA, const int _jets_index) const {
    jetA->eta  = _jets[_jets_index].rap();
    jetA->phi  = _jets[_jets_index].phi_02pi();
    jetA->kt2  = jet_scale_for_algorithm(_jets[_jets_index]);
    jetA->_jets_index = _jets_index;
    jetA->NN_dist = _R2;
    jetA->NN      = NULL;
}
template <class J> inline double ClusterSequence::_bj_dist(
                const J * const jetA, const J * const jetB) const {
#ifndef FJCORE_NEW_DELTA_PHI
  double dphi = std::abs(jetA->phi - jetB->phi);
  double deta = (jetA->eta - jetB->eta);
  if (dphi > pi) {dphi = twopi - dphi;}
#else 
  double dphi = pi-std::abs(pi-std::abs(jetA->phi - jetB->phi));
  double deta = (jetA->eta - jetB->eta);
#endif 
  return dphi*dphi + deta*deta;
}
template <class J> inline double ClusterSequence::_bj_diJ(const J * const jet) const {
  double kt2 = jet->kt2;
  if (jet->NN != NULL) {if (jet->NN->kt2 < kt2) {kt2 = jet->NN->kt2;}}
  return jet->NN_dist * kt2;
}
template <class J> inline void ClusterSequence::_bj_set_NN_nocross(
                 J * const jet, J * const head, const J * const tail) const {
  double NN_dist = _R2;
  J * NN  = NULL;
  if (head < jet) {
    for (J * jetB = head; jetB != jet; jetB++) {
      double dist = _bj_dist(jet,jetB);
      if (dist < NN_dist) {
	NN_dist = dist;
	NN = jetB;
      }
    }
  }
  if (tail > jet) {
    for (J * jetB = jet+1; jetB != tail; jetB++) {
      double dist = _bj_dist(jet,jetB);
      if (dist < NN_dist) {
	NN_dist = dist;
	NN = jetB;
      }
    }
  }
  jet->NN = NN;
  jet->NN_dist = NN_dist;
}
template <class J> inline void ClusterSequence::_bj_set_NN_crosscheck(J * const jet, 
		    J * const head, const J * const tail) const {
  double NN_dist = _R2;
  J * NN  = NULL;
  for (J * jetB = head; jetB != tail; jetB++) {
    double dist = _bj_dist(jet,jetB);
    if (dist < NN_dist) {
      NN_dist = dist;
      NN = jetB;
    }
    if (dist < jetB->NN_dist) {
      jetB->NN_dist = dist;
      jetB->NN = jet;
    }
  }
  jet->NN = NN;
  jet->NN_dist = NN_dist;
}
FJCORE_END_NAMESPACE
#endif // __FJCORE_CLUSTERSEQUENCE_HH__
#ifndef __FJCORE_NNBASE_HH__
#define __FJCORE_NNBASE_HH__
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
class _NoInfo {};
template<class I> class NNInfo {
public:
  NNInfo()         : _info(NULL) {}
  NNInfo(I * info) : _info(info) {}
  template<class BJ> void init_jet(BJ * briefjet, const fjcore::PseudoJet & jet, int index) { briefjet->init(jet, index, _info);}
private:
  I * _info;
};
template<> class NNInfo<_NoInfo>  {
public:
  NNInfo()           {}
  NNInfo(_NoInfo * ) {}
  template<class BJ> void init_jet(BJ * briefjet, const fjcore::PseudoJet & jet, int index) { briefjet->init(jet, index);}
};
template<class I = _NoInfo> class NNBase : public NNInfo<I> {
public:
  NNBase() {}
  NNBase(I * info) : NNInfo<I>(info) {}
  virtual void start(const std::vector<PseudoJet> & jets) = 0;
  virtual double dij_min(int & iA, int & iB) = 0;
  virtual void remove_jet(int iA) = 0;
  virtual void merge_jets(int iA, int iB, const PseudoJet & jet, int jet_index) =  0;
  virtual ~NNBase() {};
};
FJCORE_END_NAMESPACE      // defined in fastjet/internal/base.hh
#endif // __FJCORE_NNBASE_HH__
#ifndef __FJCORE_NNH_HH__
#define __FJCORE_NNH_HH__
FJCORE_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
template<class BJ, class I = _NoInfo> class NNH : public NNBase<I> {
public:
  NNH(const std::vector<PseudoJet> & jets)           : NNBase<I>()     {start(jets);}
  NNH(const std::vector<PseudoJet> & jets, I * info) : NNBase<I>(info) {start(jets);}
  void start(const std::vector<PseudoJet> & jets);
  double dij_min(int & iA, int & iB);
  void remove_jet(int iA);
  void merge_jets(int iA, int iB, const PseudoJet & jet, int jet_index);
  ~NNH() {
    delete[] briefjets;
  }
private:
  class NNBJ; // forward declaration
  void set_NN_crosscheck(NNBJ * jet, NNBJ * begin, NNBJ * end);
  void set_NN_nocross   (NNBJ * jet, NNBJ * begin, NNBJ * end);
  NNBJ * briefjets;
  NNBJ * head, * tail;
  int n;
  std::vector<NNBJ *> where_is;
  class NNBJ : public BJ {
  public:
    void init(const PseudoJet & jet, int index_in) {
      BJ::init(jet);
      other_init(index_in);
    }
    void init(const PseudoJet & jet, int index_in, I * info) {
      BJ::init(jet, info);
      other_init(index_in);
    }
    void other_init(int index_in) {
      _index = index_in;
      NN_dist = BJ::beam_distance();
      NN = NULL;
    }
    int index() const {return _index;}
    double NN_dist;
    NNBJ * NN;
  private:
    int _index;
  };
};
template<class BJ, class I> void NNH<BJ,I>::start(const std::vector<PseudoJet> & jets) {
  n = jets.size();
  briefjets = new NNBJ[n];
  where_is.resize(2*n);
  NNBJ * jetA = briefjets;
  for (int i = 0; i< n; i++) {
    this->init_jet(jetA, jets[i], i);
    where_is[i] = jetA;
    jetA++; // move on to next entry of briefjets
  }
  tail = jetA; // a semaphore for the end of briefjets
  head = briefjets; // a nicer way of naming start
  for (jetA = head + 1; jetA != tail; jetA++) {
    set_NN_crosscheck(jetA, head, jetA);
  }
}
template<class BJ, class I> double NNH<BJ,I>::dij_min(int & iA, int & iB) {
  double diJ_min = briefjets[0].NN_dist;
  int diJ_min_jet = 0;
  for (int i = 1; i < n; i++) {
    if (briefjets[i].NN_dist < diJ_min) {
      diJ_min_jet = i; 
      diJ_min  = briefjets[i].NN_dist;
    }
  }
  NNBJ * jetA = & briefjets[diJ_min_jet];
  iA = jetA->index();
  iB = jetA->NN ? jetA->NN->index() : -1;
  return diJ_min;
}
template<class BJ, class I> void NNH<BJ,I>::remove_jet(int iA) {
  NNBJ * jetA = where_is[iA];
  tail--; n--;
  *jetA = *tail;
  where_is[jetA->index()] = jetA;
  for (NNBJ * jetI = head; jetI != tail; jetI++) {
    if (jetI->NN == jetA) set_NN_nocross(jetI, head, tail);
    if (jetI->NN == tail) {jetI->NN = jetA;}
  }
}
template<class BJ, class I> void NNH<BJ,I>::merge_jets(int iA, int iB, 
					const PseudoJet & jet, int index) {
  NNBJ * jetA = where_is[iA];
  NNBJ * jetB = where_is[iB];
  if (jetA < jetB) std::swap(jetA,jetB);
  this->init_jet(jetB, jet, index);
  if (index >= int(where_is.size())) where_is.resize(2*index);
  where_is[jetB->index()] = jetB;
  tail--; n--;
  *jetA = *tail;
  where_is[jetA->index()] = jetA;
  for (NNBJ * jetI = head; jetI != tail; jetI++) {
    if (jetI->NN == jetA || jetI->NN == jetB) {
      set_NN_nocross(jetI, head, tail);
    } 
    double dist = jetI->distance(jetB);
    if (dist < jetI->NN_dist) {
      if (jetI != jetB) {
	jetI->NN_dist = dist;
	jetI->NN = jetB;
      }
    }
    if (dist < jetB->NN_dist) {
      if (jetI != jetB) {
	jetB->NN_dist = dist;
	jetB->NN      = jetI;
      }
    }
    if (jetI->NN == tail) {jetI->NN = jetA;}
  }
}
template <class BJ, class I> void NNH<BJ,I>::set_NN_crosscheck(NNBJ * jet, 
		    NNBJ * begin, NNBJ * end) {
  double NN_dist = jet->beam_distance();
  NNBJ * NN      = NULL;
  for (NNBJ * jetB = begin; jetB != end; jetB++) {
    double dist = jet->distance(jetB);
    if (dist < NN_dist) {
      NN_dist = dist;
      NN = jetB;
    }
    if (dist < jetB->NN_dist) {
      jetB->NN_dist = dist;
      jetB->NN = jet;
    }
  }
  jet->NN = NN;
  jet->NN_dist = NN_dist;
}
template <class BJ, class I>  void NNH<BJ,I>::set_NN_nocross(
                 NNBJ * jet, NNBJ * begin, NNBJ * end) {
  double NN_dist = jet->beam_distance();
  NNBJ * NN      = NULL;
  if (begin < jet) {
    for (NNBJ * jetB = begin; jetB != jet; jetB++) {
      double dist = jet->distance(jetB);
      if (dist < NN_dist) {
	NN_dist = dist;
	NN = jetB;
      }
    }
  }
  if (end > jet) {
    for (NNBJ * jetB = jet+1; jetB != end; jetB++) {
      double dist = jet->distance (jetB);
      if (dist < NN_dist) {
	NN_dist = dist;
	NN = jetB;
      }
    }
  }
  jet->NN = NN;
  jet->NN_dist = NN_dist;
}
FJCORE_END_NAMESPACE      // defined in fastjet/internal/base.hh
#endif // __FJCORE_NNH_HH__
#endif
