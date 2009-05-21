// -*- C++ -*-
///////////////////////////////////////////////////////////////////////////////
// File: siscone.h                                                           //
// Description: header file for the main SISCone class                       //
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
// $Revision:: 261                                                          $//
// $Date:: 2008-07-23 17:54:30 +0200 (Wed, 23 Jul 2008)                     $//
///////////////////////////////////////////////////////////////////////////////

#ifndef __SPH_SISCONE_H__
#define __SPH_SISCONE_H__

#include "protocones.h"
#include "split_merge.h"

namespace siscone_spherical{

/**
 * \class CSphsiscone
 * final class: gather everything to compute the jet contents.
 * 
 * This is the class user should use.
 * It computes the jet contents of a list of particles
 * given a cone radius and a threshold for splitting/merging.
 *
 * After the call to 'perform', the vector jets is filled with
 * the jets found. the 'contents' field of each jets contains
 * the indices of the particles included in that jet. 
 */
class CSphsiscone : public CSphstable_cones, public CSphsplit_merge{
 public:
  /// default ctor
  CSphsiscone();

  /// default dtor
  ~CSphsiscone();

  /**
   * compute the jets from a given particle set.
   * We are doing multiple passes such pass n_pass looks for jets among 
   * all particles not put into jets during previous passes.
   * By default the number of passes is infinite (0). 
   * \param _particles   list of particles
   * \param _radius      cone radius
   * \param _f           shared energy threshold for splitting&merging
   * \param _n_pass_max  maximum number of passes (0=full search)
   * \param _Emin        minimum energy of the protojets
   * \param _split_merge_scale    the scale choice for the split-merge procedure
   *        NOTE: SM_Etilde  
   *                    is always IR safe
   *              SM_E
   *                    is IR unsafe for events with mom. conservation
   *
   * \return the number of jets found.
   */
  int compute_jets(std::vector<CSphmomentum> &_particles, double _radius, double _f, 
		   int _n_pass_max=0, double _Emin=0.0,
		   Esplit_merge_scale _split_merge_scale=SM_Etilde);

  /**
   * recompute the jets with a different overlap parameter.
   * we use the same particles and R as in the preceeding call.
   * \param _f           shared energy threshold for splitting&merging
   * \param _Emin        minimum energy of the protojets
   * \param _split_merge_scale    the scale choice for the split-merge procedure
   *                                           split--merge variable
   *        NOTE: using pt leads to IR unsafety for some events with momentum
   *              conservation. So we strongly advise not to change the default
   *              value.
   * \return the number of jets found, -1 if recomputation not allowed.
   */
  int recompute_jets(double _f, double _Emin = 0.0,
		     Esplit_merge_scale _split_merge_scale=SM_Etilde);

  /// list of protocones found pass-by-pass
  std::vector<std::vector<CSphmomentum> > protocones_list;

  // random number initialisation
  static bool init_done;      ///< check random generator initialisation

#ifdef DEBUG_STABLE_CONES
  int nb_hash_cones_total, nb_hash_occupied_total;
#endif

 private:
  bool rerun_allowed;         ///< is recompute_jets allowed ?
};

  
// finally, a bunch of functions to access to 
// basic information (package name, version)
//---------------------------------------------

/** 
 * return SISCone package name.
 * This is nothing but "SISCone", it is a replacement to the
 * PACKAGE_NAME string defined in config.h and which is not
 * public by default.
 * \return the SISCone name as a string
 */
std::string siscone_package_name();

/** 
 * return SISCone version number.
 * \return a string of the form "X.Y.Z" with possible additional tag
 *         (alpha, beta, devel) to mention stability status
 */
std::string siscone_version();

}
#endif
