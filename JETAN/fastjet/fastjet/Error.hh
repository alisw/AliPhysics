//STARTHEADER
// $Id: Error.hh 389 2006-12-14 18:29:16Z salam $
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



#ifndef __FASTJET_ERROR_HH__
#define __FASTJET_ERROR_HH__

#include<iostream>
#include<string>
#include "fastjet/internal/base.hh"

namespace fastjet {      // defined in fastjet/internal/base.hh

/// class corresponding to errors that will be thrown by fastjet
class Error {
public:
  // constructors
  Error() {;};
  Error(const std::string & message) {
    _message = message; 
    if (_print_errors) std::cerr << "fastjet::Error: "<<message << std::endl;
  };

  std::string message() const {return _message;};

  static void set_print_errors(bool print_errors) {
    _print_errors = print_errors;};

private:
  std::string _message;
  static bool _print_errors;
};


} // fastjet namespace 

#endif // __FASTJET_ERROR_HH__
