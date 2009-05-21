// -*- C++ -*-
///////////////////////////////////////////////////////////////////////////////
// File: hash.h                                                              //
// Description: header file for classes hash_element and hash_cones          //
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

#ifndef __SPH_HASH_H__
#define __SPH_HASH_H__

#include "momentum.h"

namespace siscone_spherical{

/**
 * \class sph_hash_element
 * information on store cones candidates.
 *
 * We store in this class the information used to count all 
 * protocones in a first pass. This first pass only count
 * cones with different references and test their stbility
 * with the parent-child particles (border particles).
 */
class sph_hash_element{
 public:
  CSph3vector centre;  ///< centre of the cone
  bool is_stable;      ///< true if stable w.r.t. "border particles"

  sph_hash_element *next;  ///< pointer to the next element
};

/**
 * \class sph_hash_cones
 * list of cones candidates.
 *
 * We store in this class all the hash_elements and give
 * functions to manipulate them.
 */
class sph_hash_cones{
 public:
  /// constructor with initialisation
  /// \param _Np  number of particles
  /// \param _R2  cone radius (squared)
  sph_hash_cones(int _Np, double _R2);

  /// destructor
  ~sph_hash_cones();

  /**
   * insert a new candidate into the hash.
   * \param v       4-momentum of te cone to add
   * \param parent  parent particle defining the cone
   * \param child   child particle defining the cone
   * \param p_io    whether the parent has to belong to the cone or not
   * \param c_io    whether the child has to belong to the cone or not
   * \return 0 on success, 1 on error
   */
  int insert(CSphmomentum *v, CSphmomentum *parent, CSphmomentum *child, bool p_io, bool c_io);

  /**
   * insert a new candidate into the hash.
   * \param v       4-momentum of te cone to add
   * Note, in this case, we assume stability. We also assume
   * that eta and phi are computed for v
   * \return 0 on success, 1 on error
   */
  int insert(CSphmomentum *v);

  /// the cone data itself
  sph_hash_element **hash_array;

  /// number of elements
  int n_cones;

  /// number of occupied cells
#ifdef DEBUG_STABLE_CONES
  int n_occupied_cells;
#endif

  /// number of cells-1
  int mask;

  /// circle radius (squared)
  /// NOTE: need to be set before any call to 'insert'
  double R2;

  /// its squreed tangent
  double tan2R;
};

}
#endif
