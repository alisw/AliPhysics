// -*- C++ -*-
///////////////////////////////////////////////////////////////////////////////
// File: split_merge.h                                                       //
// Description: header file for splitting/merging (contains the CJet class)  //
// This file is part of the SISCone project.                                 //
// WARNING: this is not the main SISCone trunk but                           //
//          an adaptation to spherical coordinates                           //
// For more details, see http://projects.hepforge.org/siscone                //
//                                                                           //
// Copyright (c) 2006-2008 Gavin Salam and Gregory Soyez                          //
//                                                                           //
// This program is free software; you can redistribute it and/or modify      //
// it under the terms of the GNU General Public License as published by      //
// the Free Software Foundation; either version 2 of the License, or         //
// (at your option) any later version.                                       //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program; if not, write to the Free Software               //
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA //
//                                                                           //
// $Revision:: 268                                                          $//
// $Date:: 2009-03-12 21:24:16 +0100 (Thu, 12 Mar 2009)                     $//
///////////////////////////////////////////////////////////////////////////////

#ifndef __SPH_SPLIT_MERGE_H__
#define __SPH_SPLIT_MERGE_H__

#include <siscone/defines.h>
#include "geom_2d.h"
#include "momentum.h"
#include <stdio.h>
#include <vector>
#include <set>
#include <memory>
#include <string>

namespace siscone_spherical{

/**
 * \class CSphjet
 * real Jet information.
 *
 * This class contains information for one single jet. 
 * That is, first, its momentum carrying information
 * about its centre and pT, and second, its particle
 * contents
 */
class CSphjet{
 public:
  /// default ctor
  CSphjet();

  /// default dtor
  ~CSphjet();

  CSphmomentum v;            ///< jet momentum
  double E_tilde;            ///< sum of E_i [1+sin^2(theta_iJ)]
  int n;                     ///< number of particles inside
  std::vector<int> contents; ///< particle contents (list of indices)

  /// ordering variable used for ordering and overlap in the
  /// split--merge. This variable is automatically set either to
  /// E_tilde or E depending on the siscone parameter.
  ///
  /// Note that the default behaviour is E_tilde and that other
  /// choices may lead to infrared unsafe situations.
  ///
  ///  Note: we use the square of the varible rather than the variable
  /// itself
  ///
  ///  Note 2: for the overlap computation, we shall use the jet energy!
  double sm_var2;

  /// covered range in eta-phi
  CSphtheta_phi_range range;

  /// pass at which the jet has been found
  /// It starts at 0 (first pass), -1 means infinite rapidity
  int pass;
};

/// ordering of jets in pt
///bool jets_pt_less(const CSphjet &j1, const CSphjet &j2);
  
/// ordering of jets in energy (e.g. used in final jets ordering)
bool jets_E_less(const CSphjet &j1, const CSphjet &j2);
  

/// the choices of scale variable that can be used in the split-merge
/// step, both for ordering the protojets and for measuing their
/// overlap; E is defined in E-scheme (4-momentum) recombination;
/// Etilde = \sum_{i\in jet} E_i [1+sin^2(theta_{i,jet})]
///
/// NB: if one changes the order here, one _MUST_ also change the order
///     in the SISCone plugin
enum Esplit_merge_scale {
  SM_E,        ///< Energy (IR unsafe with momentum conservation)
  SM_Etilde    ///< sum_{i \in jet} E_i [1+sin^2(theta_iJ)]
};

/// return the name of the split-merge scale choice
std::string split_merge_scale_name(Esplit_merge_scale sms);

/**
 * \class CSphsplit_merge_ptcomparison
 * a class that allows us to carry out comparisons of pt of jets, using
 * information from exact particle contents where necessary.
 */
class CSphsplit_merge_ptcomparison{
public:
  /// default ctor
  CSphsplit_merge_ptcomparison() : 
    particles(0), split_merge_scale(SM_Etilde){};

  /// return the name corresponding to the SM scale variable
  std::string SM_scale_name() const {
    return split_merge_scale_name(split_merge_scale);}

  std::vector<CSphmomentum> * particles;  ///< pointer to the list of particles
  std::vector<double> * particles_norm2;  ///< pointer to the particles's norm^2

  /// comparison of 2 CSphjet
  bool operator()(const CSphjet &jet1, const CSphjet &jet2) const;

  /**
   * get the difference between 2 jets, calculated such that rounding
   * errors will not affect the result even if the two jets have
   * almost the same content (so that the difference is below the
   * rounding errors)
   *
   * \param j1        first jet
   * \param j2        second jet
   * \param v         jet1-jet2
   * \param E_tilde   jet1-jet2 E_tilde
   */
  void get_difference(const CSphjet &j1, const CSphjet &j2, CSphmomentum *v, double *E_tilde) const;

  /// the following parameter controls the variable we're using for 
  /// the split-merge process i.e. the variable we use for 
  ///  1. ordering jet candidates;
  ///  2. computing the overlap fraction of two candidates.
  /// The default value uses Etilde. The other alternative is E
  /// NOTE: Modifying the default choice can have nasty effects:
  /// - using E is IR-safe for QCD jets but the IR unsafety remains
  ///   for back-to-back jets of unstable narrow-width particles
  ///   (e.g. Higgs).
  /// Therefore, keeping the default value is strongly advised.
  Esplit_merge_scale split_merge_scale;
};


// iterator types
/// iterator definition for the jet candidates structure
typedef std::multiset<siscone_spherical::CSphjet,CSphsplit_merge_ptcomparison>::iterator cjet_iterator;

/// iterator definition for the jet structure
typedef std::vector<siscone_spherical::CSphjet>::iterator jet_iterator;



/**
 * \class CSphsplit_merge
 * Class used to split and merge jets.
 */
class CSphsplit_merge{
 public:
  /// default ctor
  CSphsplit_merge();

  /// default dtor
  ~CSphsplit_merge();


  //////////////////////////////
  // initialisation functions //
  //////////////////////////////

  /**
   * initialisation function
   * \param _particles  list of particles
   * \param protocones  list of protocones (initial jet candidates)
   * \param R2          cone radius (squared)
   * \param Emin        minimal energy allowed for jets
   * \return 0 on success, 1 on error
   */
  int init(std::vector<CSphmomentum> &_particles, std::vector<CSphmomentum> *protocones, double R2, double Emin=0.0);

  /**
   * initialisation function for particle list
   * \param _particles  list of particles
   * \return 0 on success, 1 on error
   */
  int init_particles(std::vector<CSphmomentum> &_particles);

  /**
   * build initial list of left particles
   */
  int init_pleft();

  /**
   * use an energy-dependent boundary for splitting
   * When called with true, the criterium for splitting two protojets 
   * will be to compare D1^2/kt1^2 vs. D2^2/kt2^2, the (anti-)kt-weighted 
   * distance instead of the plain distance D1^2 vs. D2^2.
   * This can be set in order to produce more circular hard jets, 
   * with the same underlying philosophy as for the anti-kt algorithm.
   * We thus expect a behaviour closer to the IterativeCone one. 
   * By default, we use the standard D1^2 vs. D2^2 comparison and this 
   * function is not called.
   */
  inline int set_E_weighted_splitting(bool _use_E_weighted_splitting){
    use_E_weighted_splitting = _use_E_weighted_splitting;
    return 0;
  }

  ////////////////////////
  // cleaning functions //
  ////////////////////////

  /// partial clearance
  int partial_clear();

  /// full clearance
  int full_clear();


  /////////////////////////////////
  // main parts of the algorithm //
  /////////////////////////////////
 
  /**
   * build the list 'p_uncol_hard' from p_remain by clustering
   * collinear particles 
   * note that thins in only used for stable-cone detection 
   * so the parent_index field is unnecessary
   *
   * Note that soft particles are not removed here
   * This is just a remnant from the trunk version
   */
  int merge_collinear_and_remove_soft();

  /**
   * add a list of protocones
   * \param protocones  list of protocones (initial jet candidates)
   * \param R2          cone radius (squared)
   * \param Emin        minimal emergy allowed for jets
   * \return 0 on success, 1 on error
   */
  int add_protocones(std::vector<CSphmomentum> *protocones, double R2, double Emin=0.0);

  /**
   * really do the splitting and merging
   * At the end, the vector jets is filled with the jets found.
   * the 'contents' field of each jets contains the indices
   * of the particles included in that jet. 
   * \param overlap_tshold  threshold for splitting/merging transition
   * \param Emin            minimal energy allowed for jets
   * \return the number of jets is returned
   */
  int perform(double overlap_tshold, double Emin=0.0);


  //////////////////////////////
  // save and debug functions //
  //////////////////////////////

  /// save final jets
  /// \param flux   stream to save the jet contentss
  int save_contents(FILE *flux);

  /// show jets/candidates status
  int show();

  // particle information
  int n;                                  ///< number of particles
  std::vector<CSphmomentum> particles;    ///< list of particles
  std::vector<double> particles_norm2;    ///< norm^2 of the particle (3-vect part)
  int n_left;                             ///< numer of particles that does not belong to any jet
  std::vector<CSphmomentum> p_remain;     ///< list of particles remaining to deal with
  std::vector<CSphmomentum> p_uncol_hard; ///< list of particles remaining with collinear clustering
  int n_pass;                             ///< index of the run

  /// minimal difference in squared distance between a particle and
  /// two overlapping protojets when doing a split (useful when
  /// testing approx. collinear safety)
  double most_ambiguous_split; 

  // jets information
  std::vector<CSphjet> jets;            ///< list of jets

  // working entries
  int *indices;                      ///< maximal size array for indices works
  int idx_size;                      ///< number of elements in indices1

  /// The following flag indicates that identical protocones
  /// are to be merged automatically each time around the split-merge
  /// loop and before anything else happens.
  ///
  /// This flag is only effective if ALLOW_MERGE_IDENTICAL_PROTOCONES
  /// is set in 'defines.h'
  /// Note that this lead to infrared-unsafety so it is disabled
  /// by default
  bool merge_identical_protocones;

  /// member used for detailed comparisons of pt's
  CSphsplit_merge_ptcomparison ptcomparison;

  /// stop split--merge when the SM_var of the hardest protojet 
  /// is below this cut-off. 
  /// This is not collinear-safe so you should not use this
  /// variable unless you really know what you are doing
  /// Note that the cut-off is set on the variable squared.
  double SM_var2_hardest_cut_off;

  /// Energy cutoff for the particles to put in p_uncol_hard 
  /// this is meant to allow removing soft particles in the
  /// stable-cone search.
  double stable_cone_soft_E2_cutoff;

 private:
  /**
   * get the overlap between 2 jets
   * \param j1   first jet
   * \param j2   second jet
   * \param v    returned overlap^2 (determined by the choice of SM variable)
   * \return true if overlapping, false if disjoint
   */
  bool get_overlap(const CSphjet &j1, const CSphjet &j2, double *v);


  /**
   * split the two given jets.
   * during this procedure, the jets j1 & j2 are replaced
   * by 2 new jets. Common particles are associted to the 
   * closest initial jet.
   * \param it_j1  iterator of the first jet in 'candidates'
   * \param it_j2  iterator of the second jet in 'candidates'
   * \param j1     first jet (CSphjet instance)
   * \param j2     second jet (CSphjet instance)
   * \return true on success, false on error
   */
  bool split(cjet_iterator &it_j1, cjet_iterator &it_j2);

  /**
   * merge the two given jet.
   * during this procedure, the jets j1 & j2 are replaced
   * by 1 single jets containing both of them.
   * \param it_j1  iterator of the first jet in 'candidates'
   * \param it_j2  iterator of the second jet in 'candidates'
   * \return true on success, false on error
   */
  bool merge(cjet_iterator &it_j1, cjet_iterator &it_j2);

  /**
   * Check whether or not a jet has to be inserted in the 
   * list of protojets. If it has, set its sm_variable and
   * insert it to the list of protojets.
   * \param jet    jet to insert
   */
  bool insert(CSphjet &jet);

  /**
   * given a 4-momentum and its associated pT, return the 
   * variable tht has to be used for SM
   * \param v          4 momentum of the protojet
   * \param E_tilde    E_tilde of the protojet
   */
  double get_sm_var2(CSphmomentum &v, double &E_tilde);

  /// compute Etilde for a given jet
  void compute_Etilde(CSphjet &j);

  // jet information
  /// list of jet candidates
  std::auto_ptr<std::multiset<CSphjet,CSphsplit_merge_ptcomparison> > candidates;

  /// minimal E
  double E_min;

  /**
   * do we have or not to use the energy-weighted splitting
   * (see description for set_E_weighted_splitting)
   * This will be false by default
   */
  bool use_E_weighted_splitting;

#ifdef ALLOW_MERGE_IDENTICAL_PROTOCONES
  /// checkxor for the candidates (to avoid having twice the same contents)
  std::set<siscone::Creference> cand_refs;
#endif
};

}


#endif
