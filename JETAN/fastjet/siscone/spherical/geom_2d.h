// -*- C++ -*-
///////////////////////////////////////////////////////////////////////////////
// File: geom_2d.h                                                           //
// Description: header file for two-dimensional geometry tools               //
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

#ifndef __SPH_GEOM_2D_H__
#define __SPH_GEOM_2D_H__

#include <iostream>
#include <math.h>
#include <siscone/defines.h>
#include <siscone/geom_2d.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197
#endif

namespace siscone_spherical{

/** 
 * \class CSphtheta_phi_range
 * \brief class for holding a covering range in eta-phi
 *
 * This class deals with ranges in the eta-phi plane. It
 * implements methods to test if two ranges overlap and
 * to take the union of two overlapping intervals.
 */
class CSphtheta_phi_range{
public:
  /// default ctor
  CSphtheta_phi_range();

  /// ctor with initialisation
  /// we initialise with a centre (in theta,phi) and a radius
  /// \param c_theta   theta coordinate of the centre
  /// \param c_phi   phi coordinate of the centre
  /// \param R       radius
  CSphtheta_phi_range(double c_theta, double c_phi, double R);

  /// assignment of range
  /// \param r   range to assign to current one
  CSphtheta_phi_range& operator = (const CSphtheta_phi_range &r);

  /// add a particle to the range
  /// \param theta  theta coordinate of the particle
  /// \param phi  phi coordinate of the particle
  /// \return 0 on success, 1 on error
  int add_particle(const double theta, const double phi);

  /// theta range as a binary coding of covered cells
  unsigned int theta_range;     

  /// phi range as a binary coding of covered cells
  unsigned int phi_range;     

  /// extremal value for theta
  static double theta_min;  ///< minimal value for theta (set to 0)
  static double theta_max;  ///< maximal value for theta (set to pi)

private:
  /// return the cell index corrsponding to an theta value
  inline unsigned int get_theta_cell(double theta){
    return (unsigned int) (1 << ((int) (32*((theta-theta_min)/(theta_max-theta_min)))));
  }

  /// return the cell index corrsponding to a phi value
  inline unsigned int get_phi_cell(double phi){
    return (unsigned int) (1 << ((int) (32*phi/twopi+16)%32));
  }
};

/// test overlap
/// \param  r1  first range
/// \param  r2  second range
/// \return true if overlap, false otherwise.
bool is_range_overlap(const CSphtheta_phi_range &r1, const CSphtheta_phi_range &r2);

/// compute union
/// Note: we assume that the two intervals overlap
/// \param  r1  first range
/// \param  r2  second range
/// \return union of the two ranges
const CSphtheta_phi_range range_union(const CSphtheta_phi_range &r1, const CSphtheta_phi_range &r2);

}

#endif
