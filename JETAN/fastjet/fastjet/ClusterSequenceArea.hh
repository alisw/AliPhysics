//STARTHEADER
// $Id: ClusterSequenceArea.hh 1303 2008-10-29 16:42:09Z salam $
//
// Copyright (c) 2006-2007, Matteo Cacciari, Gavin Salam and Gregory Soyez
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

#ifndef __FASTJET_CLUSTERSEQUENCEAREA_HH__
#define __FASTJET_CLUSTERSEQUENCEAREA_HH__

#include "fastjet/ClusterSequenceAreaBase.hh"
#include "fastjet/ClusterSequenceActiveArea.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include "fastjet/ClusterSequencePassiveArea.hh"
#include "fastjet/ClusterSequenceVoronoiArea.hh"
#include "fastjet/AreaDefinition.hh"

namespace fastjet {

/// General class for user to obtain ClusterSequence with additional
/// area information.
///
/// Based on the area_def, it automatically dispatches the work to the
/// appropriate actual ClusterSequenceAreaBase-derived-class to do the
/// real work.
class ClusterSequenceArea : public ClusterSequenceAreaBase {
public:
  /// main constructor
  template<class L> ClusterSequenceArea
         (const std::vector<L> & pseudojets, 
	  const JetDefinition & jet_def,
	  const AreaDefinition & area_def_in)  : _area_def(area_def_in) {
    initialize_and_run_cswa(pseudojets, jet_def);
  }

  /// constructor with a GhostedAreaSpec
  template<class L> ClusterSequenceArea
         (const std::vector<L> & pseudojets, 
	  const JetDefinition & jet_def,
	  const GhostedAreaSpec & area_spec)   : _area_def(area_spec){
    initialize_and_run_cswa(pseudojets, jet_def);
  }

  /// constructor with a VoronoiAreaSpec
  template<class L> ClusterSequenceArea
         (const std::vector<L> & pseudojets, 
	  const JetDefinition & jet_def,
	  const VoronoiAreaSpec & area_spec)   : _area_def(area_spec){
    initialize_and_run_cswa(pseudojets, jet_def);
  }

  /// return a reference to the area definition
  const AreaDefinition & area_def() const {return _area_def;}


  /// return the area associated with the given jet
  virtual double area       (const PseudoJet & jet) const {
    return _area_base->area(jet);}

  /// return the error (uncertainty) associated with the determination
  /// of the area of this jet
  virtual double area_error (const PseudoJet & jet) const {
    return _area_base->area_error(jet);}

  /// return the 4-vector area
  virtual PseudoJet area_4vector(const PseudoJet & jet) const {
    return _area_base->area_4vector(jet);}

  // /// return the total area, up to |y|<maxrap, that is free of jets
  // virtual double empty_area(double maxrap) const {
  //   return _area_base->empty_area(maxrap);}
  // 
  // /// return something similar to the number of pure ghost jets
  // /// in the given rapidity range in an active area case.
  // /// For the local implementation we return empty_area/(0.55 pi R^2),
  // /// based on measured properties of ghost jets with kt and cam. Note
  // /// that the number returned is a double.
  // virtual double n_empty_jets(double maxrap) const {
  //   return _area_base->n_empty_jets(maxrap);

  /// return the total area, in the given rap-phi range, that is free of jets
  virtual double empty_area(const RangeDefinition & range) const {
    return _area_base->empty_area(range);}

  /// return something similar to the number of pure ghost jets
  /// in the given rap-phi range in an active area case.
  /// For the local implementation we return empty_area/(0.55 pi R^2),
  /// based on measured properties of ghost jets with kt and cam. Note
  /// that the number returned is a double.
  virtual double n_empty_jets(const RangeDefinition & range) const {
    return _area_base->n_empty_jets(range);
  }

  /// true if a jet is made exclusively of ghosts
  virtual bool is_pure_ghost(const PseudoJet & jet) const {
    return _area_base->is_pure_ghost(jet);
  }

  /// true if this ClusterSequence has explicit ghosts
  virtual bool has_explicit_ghosts() const {
    return _area_base->has_explicit_ghosts();
  }
  

  /// overload version of what's in the ClusterSequenceAreaBase class, which 
  /// additionally checks compatibility between "range" and region in which
  /// ghosts are thrown.
  virtual void get_median_rho_and_sigma(const std::vector<PseudoJet> & all_jets,
					const RangeDefinition & range, 
                                        bool use_area_4vector,
                                        double & median, double & sigma,
                                        double & mean_area,
					bool all_are_incl = false) const {
    _warn_if_range_unsuitable(range);
    ClusterSequenceAreaBase::get_median_rho_and_sigma(
                                 all_jets, range, use_area_4vector,
				 median, sigma, mean_area, all_are_incl);
  }

  /// overload version of what's in the ClusterSequenceAreaBase class,
  /// which actually just does the same thing as the base version (but
  /// since we've overridden the 5-argument version above, we have to
  /// override the 4-argument version too.
  virtual void get_median_rho_and_sigma(const RangeDefinition & range, 
                                        bool use_area_4vector,
                                        double & median, double & sigma) const {
    ClusterSequenceAreaBase::get_median_rho_and_sigma(range,use_area_4vector,
                                                      median,sigma);
  }

  /// overload version of what's in the ClusterSequenceAreaBase class,
  /// which actually just does the same thing as the base version (but
  /// since we've overridden the multi-argument version above, we have to
  /// override the 5-argument version too.
  virtual void get_median_rho_and_sigma(const RangeDefinition & range, 
                                        bool use_area_4vector,
                                        double & median, double & sigma,
					double & mean_area) const {
    ClusterSequenceAreaBase::get_median_rho_and_sigma(range,use_area_4vector,
                                                      median,sigma, mean_area);
  }


  /// overload version of what's in the ClusterSequenceAreaBase class, which 
  /// additionally checks compatibility between "range" and region in which
  /// ghosts are thrown.
  virtual void parabolic_pt_per_unit_area(double & a, double & b, 
                                          const RangeDefinition & range, 
                                          double exclude_above=-1.0, 
                                          bool use_area_4vector=false) const {
    _warn_if_range_unsuitable(range);
    ClusterSequenceAreaBase::parabolic_pt_per_unit_area(
                                a,b,range, exclude_above, use_area_4vector);
  }


private:
  
  /// print a warning if the range is unsuitable for the current
  /// calculation of the area (e.g. because ghosts do not extend
  /// far enough).
  void _warn_if_range_unsuitable(const RangeDefinition & range) const;

  template<class L> void initialize_and_run_cswa (
                                 const std::vector<L> & pseudojets, 
                                 const JetDefinition & jet_def);

  std::auto_ptr<ClusterSequenceAreaBase> _area_base;
  AreaDefinition _area_def;
  static LimitedWarning _range_warnings;

};

//----------------------------------------------------------------------
template<class L> void ClusterSequenceArea::initialize_and_run_cswa(
           const std::vector<L> & pseudojets, 
           const JetDefinition  & jet_def)
 {
  
  ClusterSequenceAreaBase * _area_base_ptr;
  switch(_area_def.area_type()) {
  case active_area:
    _area_base_ptr = new ClusterSequenceActiveArea(pseudojets, 
                                                   jet_def, 
                                                   _area_def.ghost_spec());
    break;
  case active_area_explicit_ghosts:
    _area_base_ptr = new ClusterSequenceActiveAreaExplicitGhosts(pseudojets, 
                                                   jet_def, 
                                                   _area_def.ghost_spec());
    break;
  case voronoi_area:
    _area_base_ptr = new ClusterSequenceVoronoiArea(pseudojets, 
                                                   jet_def, 
                                                   _area_def.voronoi_spec());
    break;
  case one_ghost_passive_area:
    _area_base_ptr = new ClusterSequence1GhostPassiveArea(pseudojets, 
						    jet_def, 
						    _area_def.ghost_spec());
    break;
  case passive_area:
    _area_base_ptr = new ClusterSequencePassiveArea(pseudojets, 
						    jet_def, 
						    _area_def.ghost_spec());
    break;
  default:
    std::cerr << "Error: unrecognized area_type in ClusterSequenceArea:" 
              << _area_def.area_type() << std::endl;
    exit(-1);
  }
  // now copy across the information from the area base class
  _area_base = std::auto_ptr<ClusterSequenceAreaBase>(_area_base_ptr);
  transfer_from_sequence(*_area_base);
}

} // fastjet namespace 

#endif // __FASTJET_CLUSTERSEQUENCEAREA_HH__


