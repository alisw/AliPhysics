//STARTHEADER
// $Id: CircularRange.hh 1502 2009-04-06 10:33:14Z salam $
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


#ifndef __FASTJET_CIRCULARRANGE_HH__
#define __FASTJET_CIRCULARRANGE_HH__

#include "fastjet/RangeDefinition.hh"
#include "fastjet/Error.hh"

namespace fastjet {      // defined in fastjet/internal/base.hh

class CircularRange : public fastjet::RangeDefinition {
public:
  /// constructor
  CircularRange() {_set_invalid_rapphi();}
  
  /// initialise CircularRange with a jet
  CircularRange(const fastjet::PseudoJet & jet, double distance) {
                _distance = distance;
		_rapjet = jet.rap();
		_phijet = jet.phi();
		_total_area = fastjet::pi*_distance*_distance;  }

  /// initialise CircularRange with a (rap,phi) point
  CircularRange(double rap, double phi, double distance) {
                _distance = distance;
		_rapjet = rap;
		_phijet = phi;
		_total_area = fastjet::pi*_distance*_distance;  }

  /// initialise CircularRange with just the radius parameter
  CircularRange(double distance) {
                _set_invalid_rapphi();
                _distance = distance;
		_total_area = fastjet::pi*_distance*_distance;  }
  
  /// destructor
  virtual ~CircularRange() {}
  
  /// return description of range
  virtual inline std::string description() const {
    std::ostringstream ostr;
    ostr << "CircularRange: within distance "<< _distance << " of given jet or point." ;
    return ostr.str(); }

  /// returns true since this range is localizable (i.e. set_position
  /// does something meaningful)
  virtual inline bool is_localizable() const { return true; }
  
  /// return bool according to whether (rap,phi) is in range
  virtual inline bool is_in_range(double rap, double phi) const {
     if (! _rapphi_are_valid()) {
       throw Error("Circular range used without a center having being defined (use set_position())");
     }
     double deltaphi = _phijet - phi;
     if ( deltaphi > pi) { deltaphi -= twopi; }
     else if ( deltaphi < -pi) { deltaphi += twopi; }
     bool inrange = ( (rap-_rapjet)*(rap-_rapjet) +
                deltaphi*deltaphi <= _distance*_distance );
     return inrange; }

  /// return the minimal and maximal rapidity of this range
  virtual inline void get_rap_limits(double & rapmin, double & rapmax) const {
     rapmin = _rapjet - _distance;
     rapmax = _rapjet + _distance; }

private:
  double _distance;

  /// value for phi that marks it as invalid
  const static double _invalid_phi = -1000.0;
  /// set internal phi so as to mark things as invalid
  void _set_invalid_rapphi() {_phijet = _invalid_phi;}
  /// true if rap,phi are valid (tests only phi)
  bool _rapphi_are_valid() const {return _phijet != _invalid_phi;}

};

} // fastjet namespace 

#endif // __FASTJET_CIRCULARRANGE_HH__
