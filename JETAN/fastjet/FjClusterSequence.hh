//STARTHEADER
// $Id: FjClusterSequence.hh 945 2007-11-09 08:46:04Z salam $
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


#ifndef __FJCLUSTERSEQUENCE_HH__
#define __FJCLUSTERSEQUENCE_HH__


#ifndef __BACKWARD_WARNING_V1__
#define __BACKWARD_WARNING_V1__
#warning This file includes at least one deprecated FastJet header from v1. \
All fastjet components (including plugins) should now be accessed by including fastjet/...
#endif // __BACKWARD_WARNING_V1__

#include "fastjet/ClusterSequence.hh"


/// typedef which provides backwards compatibility for 
/// user programs based on the v1 interface
typedef fastjet::ClusterSequence FjClusterSequence;

/// typedef which provides backwards compatibility for 
/// user programs based on the v1 interface
typedef fastjet::Strategy        FjStrategy;


// below follow redefinitions of all the strategy constants
// to allow a v1 legacy user to access the strategy names
// as before

/// experimental ...
const FjStrategy N2MinHeapTiled   = fastjet::N2MinHeapTiled;
/// fastest from about 50..10^4
const FjStrategy N2Tiled     = fastjet::N2Tiled;
/// legacy
const FjStrategy N2PoorTiled = fastjet::N2PoorTiled;
/// fastest below 50
const FjStrategy N2Plain     = fastjet::N2Plain;
/// worse even than the usual N^3 algorithms
const FjStrategy N3Dumb      = fastjet::N3Dumb;
/// automatic selection of the best (based on N)
const FjStrategy Best        = fastjet::Best;
/// best of the NlnN variants -- best overall for N>10^4
const FjStrategy NlnN        = fastjet::NlnN;
/// legacy N ln N using 3pi coverage of cylinder
const FjStrategy NlnN3pi     = fastjet::NlnN3pi;
/// legacy N ln N using 4pi coverage of cylinder
const FjStrategy NlnN4pi     = fastjet::NlnN4pi;
/// Chan's closest pair method (in a variant with 4pi coverage),
/// for use exclusively with the Cambridge algorithm
const FjStrategy NlnNCam4pi   = fastjet::NlnNCam4pi;
const FjStrategy NlnNCam2pi2R = fastjet::NlnNCam2pi2R;
const FjStrategy NlnNCam      = fastjet::NlnNCam; // 2piMultD


#endif //__FJCLUSTERSEQUENCE_HH__
