//STARTHEADER
// $Id: FjPseudoJet.hh 945 2007-11-09 08:46:04Z salam $
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

#ifndef __FJPSEUDOJET_HH__
#define __FJPSEUDOJET_HH__

#ifndef __BACKWARD_WARNING_V1__
#define __BACKWARD_WARNING_V1__
#warning This file includes at least one deprecated FastJet header from v1. \
All fastjet components (including plugins) should now be accessed by including fastjet/...
#endif // __BACKWARD_WARNING_V1__


#include "fastjet/PseudoJet.hh"

/// typedef which provides backwards compatibility for 
/// user programs based on the v1 interface
typedef fastjet::PseudoJet FjPseudoJet;

/// thought not officially "declared" in the docuementation, 
/// this was used in the v1 fastjet_timing.cc example program,
/// so in order to ensure that it still compiles we define it...
const double twopi = fastjet::twopi;

#endif //__FJPSEUDOJET_HH__
