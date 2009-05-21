// -*- C++ -*-
///////////////////////////////////////////////////////////////////////////////
// File: vicinity.h                                                          //
// Description: header file for particle vicinity (Cvicinity class)          //
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
// $Revision:: 255                                                          $//
// $Date:: 2008-07-12 17:40:35 +0200 (Sat, 12 Jul 2008)                     $//
///////////////////////////////////////////////////////////////////////////////

#ifndef __SPH_VICINITY_H__
#define __SPH_VICINITY_H__

#include <siscone/vicinity.h>
#include <vector>
#include <list>
#include "momentum.h"
#include <siscone/defines.h>
#ifdef USE_QUADTREE_FOR_STABILITY_TEST
#include <siscone/quadtree.h>
#endif

namespace siscone_spherical{

  
/**
 * \class CSphvicinity_elm
 * \brief element in the vicinity of a parent.
 *
 * class used to manage one points in the vicinity 
 * of a parent point.
 */
class CSphvicinity_elm{
 public:
  /// pointer to the second borderline particle
  CSphmomentum *v;

  /// variable to tell if the particle is inside or outside the cone 
  siscone::Cvicinity_inclusion *is_inside;   

  // centre variables
  CSph3vector centre;         ///< direction of the centre
  double angle;            ///< angle with parent
  bool side;               ///< true if angle on the positive side, false otherwise
  double cocircular_range; ///< amount by which the angle can be varied while
                           ///< maintaining this point within co-circularity margin

  /// list of elements co-circular with this one
  /// NB: empty list uses less mem than vector
  std::list<CSphvicinity_elm * > cocircular;                                          
};

/// ordering pointers to CSphvicinity_elm
bool ve_less(CSphvicinity_elm *ve1, CSphvicinity_elm *ve2);


/**
 * \class CSphvicinity
 * \brief list of element in the vicinity of a parent.
 *
 * class used to manage the points which are in the vicinity 
 * of a parent point.
 */
class CSphvicinity{
 public:
  /// default constructor
  CSphvicinity();

  /// constructor with initialisation (see set_particle_list)
  CSphvicinity(std::vector<CSphmomentum> &_particle_list);

  /// default destructor
  ~CSphvicinity();

  /**
   * set the particle_list
   * \param _particle_list   list of particles (type CSphmomentum)
   */ 
  void set_particle_list(std::vector<CSphmomentum> &_particle_list);

  /**
   * build the vicinity list from the list of points.
   * \param _parent    reference particle
   * \param _VR        vicinity radius
   */
  void build(CSphmomentum *_parent, double _VR);

  // cone kinematical information
  CSphmomentum *parent;      ///< parent vector
  double VR;                 ///< radius of the vicinity
  double VR2;                ///< squared radius of the vicinity
  double cosVR;              ///< cosine of the radius of the vicinity
  double R;                  ///< normal radius
  double R2;                 ///< squared normal radius
  double tan2R;              ///< squared tangent of the normal radius
  double D2_R;               ///< euclidian distance (squared) corresp. to the arc R
  double inv_R_EPS_COCIRC;   ///< R / EPSILON_COCIRCULAR
  double inv_R_2EPS_COCIRC;  ///< R / (2*EPSILON_COCIRCULAR)

  // particle list information
  int n_part;                                 ///< number of particles
  std::vector<CSphmomentum> plist;            ///< the list of particles
  /// the inclusion state of particles
  std::vector<siscone::Cvicinity_inclusion> pincluded; 
  CSphvicinity_elm *ve_list;                  ///< list of vicinity elements built from particle list (size=2*n)
#ifdef USE_QUADTREE_FOR_STABILITY_TEST
  siscone::Cquadtree *quadtree;               ///< quadtree used for final stability tests
#endif

  // vicinity information
  std::vector<CSphvicinity_elm*> vicinity;    ///< list of points in parent's vicinity
  unsigned int vicinity_size;                 ///< number of elements in vicinity

 protected:
  /**
   * append a particle to the 'vicinity' list after
   * having tested it and computed the angular-ordering quantities
   * \param v   vector to test
   */
  void append_to_vicinity(CSphmomentum *v);

  // internal variables
  CSph3vector parent_centre;    ///< parent centre
  CSph3vector angular_dir1;     ///< main direction to measure angles
  CSph3vector angular_dir2;     ///< second direction to measure angles (sign)
};

}

#endif
