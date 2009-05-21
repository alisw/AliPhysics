//STARTHEADER
// $Id: ClusterSequencePassiveArea.hh 1134 2008-03-15 17:05:16Z salam $
//
// Copyright (c) 2005-2007, Matteo Cacciari, Gavin Salam and Gregory Soyez
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

#ifndef __FASTJET_CLUSTERSEQUENCEPASSIVEAREA_HH__
#define __FASTJET_CLUSTERSEQUENCEPASSIVEAREA_HH__


#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence1GhostPassiveArea.hh"
#include<iostream>
#include<vector>

namespace fastjet {      // defined in fastjet/internal/base.hh

//using namespace std;

/// Class that behaves essentially like ClusterSequence except
/// that it also provides access to the area of a jet (which
/// will be a random quantity... Figure out what to do about seeds 
/// later...)
class ClusterSequencePassiveArea : public ClusterSequence1GhostPassiveArea {
public:

  /// constructor based on JetDefinition and PassiveAreaSpec
  template<class L> ClusterSequencePassiveArea
         (const std::vector<L> & pseudojets, 
	  const JetDefinition & jet_def,
	  const GhostedAreaSpec & area_spec,
	  const bool & writeout_combinations = false) ;

  /// return an empty area that's appropriate to the passive area
  /// determination carried out
  virtual double empty_area(const RangeDefinition & range) const;

private:

  /// does the initialisation and running specific to the passive
  /// areas class
  void _initialise_and_run_PA (const JetDefinition & jet_def,
                               const GhostedAreaSpec & area_spec,
                               const bool & writeout_combinations = false);

};




template<class L> ClusterSequencePassiveArea::ClusterSequencePassiveArea 
(const std::vector<L> & pseudojets, 
 const JetDefinition & jet_def,
 const GhostedAreaSpec & area_spec,
 const bool & writeout_combinations) {

  // transfer the initial jets (type L) into our own array
  _transfer_input_jets(pseudojets);

  // run the clustering for passive areas
  _initialise_and_run_PA(jet_def, area_spec, writeout_combinations);

}


  
} // fastjet namespace 

#endif // __FASTJET_CLUSTERSEQUENCEPASSIVEAREA_HH__
