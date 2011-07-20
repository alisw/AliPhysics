//STARTHEADER
// $Id: ClusterSequenceActiveAreaExplicitGhosts.hh 1303 2008-10-29 16:42:09Z salam $
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

#ifndef __FASTJET_CLUSTERSEQUENCEACTIVEAREAEXPLICITGHOSTS_HH_
#define __FASTJET_CLUSTERSEQUENCEACTIVEAREAEXPLICITGHOSTS_HH_ 

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceAreaBase.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/internal/LimitedWarning.hh"
#include<iostream>
#include<vector>
#include <stdio.h>

namespace fastjet {      // defined in fastjet/internal/base.hh

//======================================================================
/// Class that behaves essentially like ClusterSequence except
/// that it also provides access to the area of a jet (which
/// will be a random quantity... Figure out what to do about seeds 
/// later...)
class ClusterSequenceActiveAreaExplicitGhosts : 
  public ClusterSequenceAreaBase {
public:
  /// constructor using a GhostedAreaSpec to specify how the area is
  /// to be measured
  template<class L> ClusterSequenceActiveAreaExplicitGhosts
         (const std::vector<L> & pseudojets, 
          const JetDefinition & jet_def,
	  const GhostedAreaSpec & area_spec,
	  const bool & writeout_combinations = false) 
	   : ClusterSequenceAreaBase() {
           std::vector<L> * ghosts = NULL;
	   _initialise(pseudojets,jet_def,&area_spec,ghosts,0.0,
                       writeout_combinations); }

  template<class L> ClusterSequenceActiveAreaExplicitGhosts
         (const std::vector<L> & pseudojets, 
          const JetDefinition & jet_def,
          const std::vector<L> & ghosts,
          double ghost_area,
	  const bool & writeout_combinations = false) 
	   : ClusterSequenceAreaBase() {
           const GhostedAreaSpec * area_spec = NULL;
	   _initialise(pseudojets,jet_def,area_spec,&ghosts,ghost_area,
                       writeout_combinations); }


  /// does the actual work of initialisation
  template<class L> void _initialise
         (const std::vector<L> & pseudojets, 
          const JetDefinition & jet_def,
	  const GhostedAreaSpec * area_spec,
	  const std::vector<L> * ghosts,
	  double                 ghost_area,
	  const bool & writeout_combinations); 

  //vector<PseudoJet> constituents (const PseudoJet & jet) const;

  /// returns the number of hard particles (i.e. those supplied by the user).
  unsigned int n_hard_particles() const;

  /// returns the area of a jet
  virtual double area (const PseudoJet & jet) const;

  /// returns a four vector corresponding to the sum (E-scheme) of the
  /// ghost four-vectors composing the jet area, normalised such that
  /// for a small contiguous area the p_t of the extended_area jet is
  /// equal to area of the jet.
  virtual PseudoJet area_4vector (const PseudoJet & jet) const;

  /// true if a jet is made exclusively of ghosts
  virtual bool is_pure_ghost(const PseudoJet & jet) const;

  /// true if the entry in the history index corresponds to a
  /// ghost; if hist_ix does not correspond to an actual particle
  /// (i.e. hist_ix < 0), then the result is false.
  bool is_pure_ghost(int history_index) const;

  /// this class does have explicit ghosts
  virtual bool has_explicit_ghosts() const {return true;}

  /// return the total area, up to |y|<maxrap, that consists of
  /// unclustered ghosts
  virtual double empty_area(const RangeDefinition & range) const;

  /// returns the total area under study
  double total_area () const;
  
  /// returns the largest squared transverse momentum among
  /// all ghosts
  double max_ghost_perp2() const {return _max_ghost_perp2;}

  /// returns true if there are any particles whose transverse momentum
  /// if so low that there's a risk of the ghosts having modified the
  /// clustering sequence
  bool has_dangerous_particles() const {return _has_dangerous_particles;}

private:

  int    _n_ghosts;
  double _ghost_area;
  std::vector<bool> _is_pure_ghost;
  std::vector<double> _areas;
  std::vector<PseudoJet> _area_4vectors;
  
  // things related to checks for dangerous particles
  double _max_ghost_perp2;
  bool   _has_dangerous_particles; 
  static LimitedWarning _warnings;

  //static int _n_warn_dangerous_particles;
  //static const int _max_warn_dangerous_particles = 5;

  
  unsigned int _initial_hard_n;

  /// adds the "ghost" momenta, which will be used to estimate
  /// the jet area
  void _add_ghosts(const GhostedAreaSpec & area_spec); 

  /// another way of adding ghosts
  template<class L> void _add_ghosts (
	  const std::vector<L> & ghosts,
	  double                 ghost_area);

  /// routine to be called after the processing is done so as to
  /// establish summary information on all the jets (areas, whether
  /// pure ghost, etc.)
  void _post_process();

};


//----------------------------------------------------------------------
// initialise from some generic type... Has to be made available
// here in order for the template aspect of it to work...
template<class L> void ClusterSequenceActiveAreaExplicitGhosts::_initialise
         (const std::vector<L> & pseudojets, 
          const JetDefinition & jet_def,
	  const GhostedAreaSpec * area_spec,
	  const std::vector<L> * ghosts,
	  double                 ghost_area,
	  const bool & writeout_combinations) {
  // don't reserve space yet -- will be done below

  // insert initial jets this way so that any type L that can be
  // converted to a pseudojet will work fine (basically PseudoJet
  // and any type that has [] subscript access to the momentum
  // components, such as CLHEP HepLorentzVector).
  for (unsigned int i = 0; i < pseudojets.size(); i++) {
    PseudoJet mom(pseudojets[i]);
    //mom.set_user_index(0); // for user's particles (user index now lost...)
    _jets.push_back(mom);
    _is_pure_ghost.push_back(false);
  }

  _initial_hard_n = _jets.size();

  if (area_spec != NULL) {
    _add_ghosts(*area_spec);
  } else {
    _add_ghosts(*ghosts, ghost_area);
  }

  if (writeout_combinations) {
    std::cout << "# Printing particles including ghosts\n";
    for (unsigned j = 0; j < _jets.size(); j++) {
      printf("%5u %20.13f %20.13f %20.13e\n",
	       j,_jets[j].rap(),_jets[j].phi_02pi(),_jets[j].kt2());
    }
    std::cout << "# Finished printing particles including ghosts\n";
  }

  // this will ensure that we can still point to jets without
  // difficulties arising!
  _jets.reserve(_jets.size()*2);

  // run the clustering
  _initialise_and_run(jet_def,writeout_combinations);

  // set up all other information
  _post_process();
}


inline unsigned int ClusterSequenceActiveAreaExplicitGhosts::n_hard_particles() const {return _initial_hard_n;}


//----------------------------------------------------------------------
/// add an explicitly specified bunch of ghosts
template<class L> void ClusterSequenceActiveAreaExplicitGhosts::_add_ghosts (
	  const std::vector<L> & ghosts,
	  double                 ghost_area) {

  
  for (unsigned i = 0; i < ghosts.size(); i++) {
    _is_pure_ghost.push_back(true);
    _jets.push_back(ghosts[i]);
  }
  // and record some info about ghosts
  _ghost_area = ghost_area;
  _n_ghosts   = ghosts.size();
}


} // fastjet namespace 

#endif // __FASTJET_CLUSTERSEQUENCEACTIVEAREAEXPLICITGHOSTS_HH_
