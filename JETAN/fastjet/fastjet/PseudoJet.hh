//STARTHEADER
// $Id: PseudoJet.hh 1510 2009-04-13 08:48:41Z salam $
//
// Copyright (c) 2005-2006, Matteo Cacciari and Gavin Salam
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
//  along with FastJet; if not, write to the Free Software
//  Foundation, Inc.:
//      59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//----------------------------------------------------------------------
//ENDHEADER


#ifndef __FASTJET_PSEUDOJET_HH__
#define __FASTJET_PSEUDOJET_HH__

#include<valarray>
#include<vector>
#include<cassert>
#include<cmath>
#include<iostream>
#include "fastjet/internal/numconsts.hh"

namespace fastjet {      // defined in fastjet/internal/base.hh

//using namespace std;

/// Used to protect against parton-level events where pt can be zero
/// for some partons, giving rapidity=infinity. KtJet fails in those cases.
const double MaxRap = 1e5;

/// Class to contain pseudojets, including minimal information of use to
/// to jet-clustering routines.
class PseudoJet {

 public:
  PseudoJet() {};
  /// construct a pseudojet from explicit components
  PseudoJet(const double px, const double py, const double pz, const double E);
  /// constructor from any object that has px,py,pz,E = some_four_vector[0--3],
  template <class L> PseudoJet(const L & some_four_vector) ;

  // first "const double &" says that result is a reference to the
  // stored value and that we will not change that stored value.
  //
  // second "const" says that "this" will not be modified by these
  // functions.
  inline double E()   const {return _E;}
  inline double e()   const {return _E;} // like CLHEP
  inline double px()  const {return _px;}
  inline double py()  const {return _py;}
  inline double pz()  const {return _pz;}

  /// returns phi (in the range 0..2pi)
  inline double phi() const {return phi_02pi();}

  /// returns phi in the range -pi..pi
  inline double phi_std()  const {
    return _phi > pi ? _phi-twopi : _phi;}

  /// returns phi in the range 0..2pi
  inline double phi_02pi() const {return _phi;}

  /// returns the rapidity or some large value when the rapidity
  /// is infinite
  inline double rap() const {return _rap;}

  /// the same as rap()
  inline double rapidity() const {return _rap;} // like CLHEP

  /// returns the pseudo-rapidity or some large value when the
  /// rapidity is infinite
  double pseudorapidity() const;
  double eta() const {return pseudorapidity();}

  /// returns the squared transverse momentum
  inline double kt2() const {return _kt2;}
  /// returns the squared transverse momentum
  inline double perp2() const {return _kt2;}  // like CLHEP
  /// returns the scalar transverse momentum
  inline double  perp() const {return sqrt(_kt2);}    // like CLHEP
  /// returns the squared invariant mass // like CLHEP
  inline double  m2() const {return (_E+_pz)*(_E-_pz)-_kt2;}    
  /// returns the squared transverse mass = kt^2+m^2
  inline double mperp2() const {return (_E+_pz)*(_E-_pz);}
  /// returns the transverse mass = sqrt(kt^2+m^2)
  inline double mperp() const {return sqrt(std::abs(mperp2()));}
  /// returns the invariant mass 
  /// (If m2() is negative then -sqrt(-m2()) is returned, as in CLHEP)
  inline double  m() const;    
  /// return px^2+py^2+pz^2
  inline double modp2() const {return _kt2+_pz*_pz;}
  /// return the transverse energy
  inline double Et() const {return (_kt2==0) ? 0.0 : _E/sqrt(1.0+_pz*_pz/_kt2);}
  /// return the transverse energy squared
  inline double Et2() const {return (_kt2==0) ? 0.0 : _E*_E/(1.0+_pz*_pz/_kt2);}

  /// returns component i, where X==0, Y==1, Z==2, E==3
  double operator () (int i) const ; 
  /// returns component i, where X==0, Y==1, Z==2, E==3
  inline double operator [] (int i) const { return (*this)(i); }; // this too


  // taken from CLHEP
  enum { X=0, Y=1, Z=2, T=3, NUM_COORDINATES=4, SIZE=NUM_COORDINATES };


  /// transform this jet (given in the rest frame of prest) into a jet
  /// in the lab frame [NOT FULLY TESTED]
  PseudoJet & boost(const PseudoJet & prest);
  /// transform this jet (given in lab) into a jet in the rest
  /// frame of prest  [NOT FULLY TESTED]
  PseudoJet & unboost(const PseudoJet & prest);

  /// return the cluster_hist_index, intended to be used by clustering
  /// routines.
  inline int cluster_hist_index() const {return _cluster_hist_index;}
  /// set the cluster_hist_index, intended to be used by clustering routines.
  inline void set_cluster_hist_index(const int index) {_cluster_hist_index = index;}

  /// alternative name for cluster_hist_index() [perhaps more meaningful]
  inline int cluster_sequence_history_index() const {
    return cluster_hist_index();}
  /// alternative name for set_cluster_hist_index(...) [perhaps more
  /// meaningful]
  inline void set_cluster_sequence_history_index(const int index) {
    set_cluster_hist_index(index);}


  /// return the user_index, intended to allow the user to "add" information
  inline int user_index() const {return _user_index;}
  /// set the user_index, intended to allow the user to "add" information
  inline void set_user_index(const int index) {_user_index = index;}

  /// return a valarray containing the four-momentum (components 0-2
  /// are 3-mom, component 3 is energy).
  std::valarray<double> four_mom() const;

  /// returns kt distance (R=1) between this jet and another
  double kt_distance(const PseudoJet & other) const;

  /// returns squared cylinder (rap-phi) distance between this jet and another
  double plain_distance(const PseudoJet & other) const;
  /// returns squared cylinder (rap-phi) distance between this jet and
  /// another
  inline double squared_distance(const PseudoJet & other) const {
    return plain_distance(other);}

  /// returns other.phi() - this.phi(), constrained to be in 
  /// range -pi .. pi
  double delta_phi_to(const PseudoJet & other) const;

  //// this seemed to compile except if it was used
  //friend inline double 
  //  kt_distance(const PseudoJet & jet1, const PseudoJet & jet2) { 
  //                                      return jet1.kt_distance(jet2);}

  /// returns distance between this jet and the beam
  inline double beam_distance() const {return _kt2;}


  void operator*=(double);
  void operator/=(double);
  void operator+=(const PseudoJet &);
  void operator-=(const PseudoJet &);

  /// reset the 4-momentum according to the supplied components and
  /// put the user and history indices back to their default values
  inline void reset(double px, double py, double pz, double E);
  
  /// reset the PseudoJet to be equal to psjet (including its
  /// indices); NB if the argument is derived from a PseudoJet then
  /// the "reset" used will be the templated version (which does not
  /// know about indices...)
  inline void reset(const PseudoJet & psjet) {
    (*this) = psjet;
  }

  /// reset the 4-momentum according to the supplied generic 4-vector
  /// (accessible via indexing, [0]==px,...[3]==E) and put the user
  /// and history indices back to their default values.
  template <class L> inline void reset(const L & some_four_vector) {
    reset(some_four_vector[0], some_four_vector[1],
          some_four_vector[2], some_four_vector[3]);
  }

 private: 
  // NB: following order must be kept for things to behave sensibly...
  double _px,_py,_pz,_E;
  double _phi, _rap, _kt2; 
  int    _cluster_hist_index, _user_index;
  /// calculate phi, rap, kt2 based on the 4-momentum components
  void _finish_init();
  /// set the indices to default values
  void _reset_indices();

  //vertex_type * vertex0, vertex1;
};


//----------------------------------------------------------------------
// routines for basic binary operations

PseudoJet operator+(const PseudoJet &, const PseudoJet &);
PseudoJet operator-(const PseudoJet &, const PseudoJet &);
PseudoJet operator*(double, const PseudoJet &);
PseudoJet operator*(const PseudoJet &, double);
PseudoJet operator/(const PseudoJet &, double);

inline double dot_product(const PseudoJet & a, const PseudoJet & b) {
  return a.E()*b.E() - a.px()*b.px() - a.py()*b.py() - a.pz()*b.pz();
}

/// returns true if the momenta of the two input jets are identical
bool have_same_momentum(const PseudoJet &, const PseudoJet &);

/// return a pseudojet with the given pt, y, phi and mass
PseudoJet PtYPhiM(double pt, double y, double phi, double m = 0.0);

//----------------------------------------------------------------------
// Routines to do with providing sorted arrays of vectors.

/// return a vector of jets sorted into decreasing transverse momentum
std::vector<PseudoJet> sorted_by_pt(const std::vector<PseudoJet> & jets);

/// return a vector of jets sorted into increasing rapidity
std::vector<PseudoJet> sorted_by_rapidity(const std::vector<PseudoJet> & jets);

/// return a vector of jets sorted into decreasing energy
std::vector<PseudoJet> sorted_by_E(const std::vector<PseudoJet> & jets);

/// return a vector of jets sorted into increasing pz
std::vector<PseudoJet> sorted_by_pz(const std::vector<PseudoJet> & jets);

//----------------------------------------------------------------------
// some code to help sorting

/// sort the indices so that values[indices[0->n-1]] is sorted
/// into increasing order 
void sort_indices(std::vector<int> & indices, 
		  const std::vector<double> & values);

/// given a vector of values with a one-to-one correspondence with the
/// vector of objects, sort objects into an order such that the
/// associated values would be in increasing order (but don't actually
/// touch the values vector in the process).
template<class T> std::vector<T> objects_sorted_by_values(const std::vector<T> & objects, 
					      const std::vector<double> & values);

/// a class that helps us carry out indexed sorting.
class IndexedSortHelper {
public:
  inline IndexedSortHelper (const std::vector<double> * reference_values) {
    _ref_values = reference_values;
  };
  inline int operator() (const int & i1, const int & i2) const {
    return  (*_ref_values)[i1] < (*_ref_values)[i2];
  };
private:
  const std::vector<double> * _ref_values;
};


//----------------------------------------------------------------------
/// constructor from any object that has px,py,pz,E = some_four_vector[0--3],
// NB: do not know if it really needs to be inline, but when it wasn't
//     linking failed with g++ (who knows what was wrong...)
template <class L> inline  PseudoJet::PseudoJet(const L & some_four_vector) {

  _px = some_four_vector[0];
  _py = some_four_vector[1];
  _pz = some_four_vector[2];
  _E  = some_four_vector[3];
  _finish_init();
  // some default values for these two indices
  _reset_indices();
}


//----------------------------------------------------------------------
inline void PseudoJet::_reset_indices() { 
  set_cluster_hist_index(-1);
  set_user_index(-1);
}

//----------------------------------------------------------------------
/// specialization of the "reset" template for case where something
/// is reset to a pseudojet -- it then takes the user and history
/// indices from the psjet
// template<> inline void PseudoJet::reset<PseudoJet>(const PseudoJet & psjet) {
//   (*this) = psjet;
// }

////// fun and games...
////template<class L> class FJVector : public L {
//////  /** Default Constructor: create jet with no constituents */
//////  Vector<L>();
////
////};
////

// taken literally from CLHEP
inline double PseudoJet::m() const {
  double mm = m2();
  return mm < 0.0 ? -std::sqrt(-mm) : std::sqrt(mm);
}


inline void PseudoJet::reset(double px, double py, double pz, double E) {
  _px = px;
  _py = py;
  _pz = pz;
  _E  = E;
  _finish_init();
  _reset_indices();
}


} // fastjet namespace 

#endif // __FASTJET_PSEUDOJET_HH__
