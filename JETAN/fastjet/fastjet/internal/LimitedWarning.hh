//STARTHEADER
// $Id: LimitedWarning.hh 699 2007-06-04 20:20:42Z salam $
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


#ifndef __FASTJET_LIMITEDWARNING_HH__
#define __FASTJET_LIMITEDWARNING_HH__

#include<iostream>
#include<string>

/// class to provide facilities for giving warnings up to some maximum
/// number of times
class LimitedWarning {
public:
  
  /// constructor that provides a default maximum number of warnings
  LimitedWarning() : _max_warn(_max_warn_default), _n_warn_so_far(0) {}

  /// constructor that provides a used-set max number of warnings
  LimitedWarning(int max_warn) : _max_warn(max_warn), _n_warn_so_far(0) {}

  /// output a warning to ostr
  void warn(const std::string & warning, std::ostream & ostr = std::cerr) {
    if (_n_warn_so_far < _max_warn) {
      ostr << "WARNING: ";
      ostr << warning;
      _n_warn_so_far++;
      if (_n_warn_so_far == _max_warn) ostr << " (LAST SUCH WARNING)";
      ostr << std::endl;
    }
  }

private:
  int _max_warn, _n_warn_so_far;
  static const int _max_warn_default = 5;
  
};

#endif // __FASTJET_LIMITEDWARNING_HH__
