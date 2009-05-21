// -*- C++ -*-
///////////////////////////////////////////////////////////////////////////////
// File: momentum.h                                                          //
// Description: header file for 4-momentum class Cmomentum                   //
// This file is part of the SISCone project.                                 //
// WARNING: this is not the main SISCone trunk but                           //
//          an adaptation to spherical coordinates                           //
// For more details, see http://projects.hepforge.org/siscone                //
//                                                                           //
// Copyright (c) 2006-2008 Gavin Salam and Gregory Soyez                     //
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
// $Revision:: 256                                                          $//
// $Date:: 2008-07-14 13:52:16 +0200 (Mon, 14 Jul 2008)                     $//
///////////////////////////////////////////////////////////////////////////////

#ifndef __SPH_VECTOR_H__
#define __SPH_VECTOR_H__

#include <vector>
#include <math.h>
#include <siscone/reference.h>
#include "geom_2d.h"
#include <siscone/defines.h>

namespace siscone_spherical{

/**
 * \class CSph3vector
 * \brief base class for managing the spatial part of Cmomentum (defined after)
 *
 * This class contains the information for particle or group of
 * particles management.
 * It is adapted to use spherical geometry, where, for our purposes,
 * the only time-consuming operation we need is the computation of 
 * the norm. To compute it once-and-for-all and store it in a local 
 * variable, you should call the 'build_norm' method.
 * On top of that, the angle phi is computed from the x-axis
 * and theta from the "north pole". 
 */
class CSph3vector{
 public:
  /// default ctor
  CSph3vector();

  /// ctor with initialisation
  CSph3vector(double _px, double _py, double _pz);

  /// default dtor
  ~CSph3vector();

  /// assignment of vectors
  CSph3vector& operator = (const CSph3vector &v);

  /// addition of vectors
  /// WARNING= norm is not updated
  const CSph3vector operator + (const CSph3vector &v);

  /// subtraction of vectors
  /// WARNING= norm is not updated
  const CSph3vector operator - (const CSph3vector &v);

  /// division by a constant
  /// WARNING= norm is not updated
  const CSph3vector operator / (const double &r);

  /// incrementation of vectors
  /// WARNING= norm is not updated
  CSph3vector& operator += (const CSph3vector &v);

  /// decrementation of vectors
  /// WARNING= norm is not updated
  CSph3vector& operator -= (const CSph3vector &v);

  /// multiplication by a constant
  /// WARNING= norm is not updated
  CSph3vector& operator *= (const double &r);

  /// division by a constant
  /// WARNING= norm is not updated
  CSph3vector& operator /= (const double &r);

  /// computes pT
  inline double perp() const {return sqrt(perp2());}

  /// computes pT^2
  inline double perp2() const {return px*px+py*py;}

  /// 3-vect norm
  inline double norm() const {return sqrt(px*px+py*py+pz*pz);}

  /// 3-vect norm squared
  inline double norm2() const {return px*px+py*py+pz*pz;}

  /// 3-vect azimuthal angle
  inline double phi() const {return atan2(py, px);}

  /// 3-vect polar angle
  inline double theta() const {return atan2(perp(),pz);}

  /// build the spatial normfrom 4-momentum info
  /// !!!                  WARNING                       !!!
  /// !!! computing the norm is the only time-consuming  !!!
  /// !!! information we need in all computations.       !!!
  /// !!! use this whenever you need repeated access     !!!
  /// !!! to the norm to store it in the local variable  !!!
  void build_norm();

  /// just a useful tool to store theta and phi
  /// locally (in _theta and _phi) in case you need 
  /// repeated access
  void build_thetaphi();

  /// for this direction, compute the two reference directions
  /// used to measure angles
  void get_angular_directions(CSph3vector &angular_dir1, CSph3vector &angular_dir2);

  double px;        ///< x-momentum
  double py;        ///< y-momentum
  double pz;        ///< z-momentum

  double _norm;     ///< particle spatial norm (available ONLY after a call to build_norm)
  double _theta;    ///< particle theta angle (available ONLY after a call to build_thetaphi)
  double _phi;      ///< particle phi angle   (available ONLY after a call to build_thetaphi)

  //////////////////////////////////////////////
  // the following part is used for checksums //
  //////////////////////////////////////////////
  siscone::Creference ref;   ///< reference number for the vector
};

/**
 * \class CSphmomentum
 * \brief base class for dynamic coordinates management
 *
 * This class contains the information for particle or group of
 * particles management.
 * It is adapted to use spherical geometry, where, for our purposes,
 * the only time-consuming operation we need is the computation of 
 * the norm. To compute it once-and-for-all and store it in a local 
 * variable, you should call the 'build_norm' method.
 * On top of that, the angle phi is computed from the x-axis
 * and theta from the "north pole". 
 */
class CSphmomentum : public CSph3vector{
 public:
  /// default ctor
  CSphmomentum();

  /// init from a 3-vect
  CSphmomentum(CSph3vector &init, double E=0.0);

  /// ctor with initialisation
  CSphmomentum(double _px, double _py, double _pz, double _E);

  /// ctor with detailed initialisation
  //CSphmomentum(double _eta, double _phi, siscone::Creference _ref);

  /// default dtor
  ~CSphmomentum();

  /// computes m
  inline double mass() const {return sqrt(mass2());}

  /// computes m^2
  inline double mass2() const {return perpmass2()-perp2();}

  /// transverse mass, mt = sqrt(pt^2+m^2) = sqrt(E^2 - pz^2)
  inline double perpmass() const {return sqrt((E-pz)*(E+pz));}

  /// transverse mass squared, mt^2 = pt^2+m^2 = E^2 - pz^2
  inline double perpmass2() const {return (E-pz)*(E+pz);}

  /// computes transverse energy
  inline double Et() const {return E/sqrt(1.0+pz*pz/perp2());}

  /// computes transverse energy (squared)
  inline double Et2() const {return E*E/(1.0+pz*pz/perp2());}

  /// assignment of vectors
  CSphmomentum& operator = (const CSphmomentum &v);

  /// addition of vectors
  /// !!! WARNING !!! no updating of eta and phi !!!
  const CSphmomentum operator + (const CSphmomentum &v);

  /// incrementation of vectors
  /// !!! WARNING !!! no updating of eta and phi !!!
  CSphmomentum& operator += (const CSphmomentum &v);

  /// decrementation of vectors
  /// !!! WARNING !!! no updating of eta and phi !!!
  CSphmomentum& operator -= (const CSphmomentum &v);

  double E;         ///< energy

  int parent_index; ///< particle number in the parent list
  int index;        ///< internal particle number
};

/// ordering of two vectors
/// this is by default done w.r.t. their references
bool operator < (const CSphmomentum &v1, const CSphmomentum &v2);

/// ordering of vectors in eta (e.g. used in collinear tests)
bool momentum_theta_less(const CSphmomentum &v1, const CSphmomentum &v2);

/// ordering of vectors in pt
bool momentum_pt_less(const CSphmomentum &v1, const CSphmomentum &v2);


//////////////////////////
// some handy utilities //
//////////////////////////

/// square
inline double sqr(double x){return x*x;}

/// dot product for te spatial 3-vect
/// \param v1    first 4-vect
/// \param v2    second 4-vect
inline double dot_product3(const CSph3vector &v1, const CSph3vector &v2){
  //double tmp = v1.px*v2.px + v1.py*v2.py + v1.pz*v2.pz;
  //if (!isfinite(tmp)){
  //  std::cout << "dot_product inf: " << std::endl;
  //  std::cout << "  angles: " << v1._theta << " " << v1._phi << " and " << v2._theta << " " << v2._phi << std::endl; 
  //  std::cout << "  moms  : " << v1.px << " " << v1.py << " " << v1.pz 
  //	      << " and "      << v2.px << " " << v2.py << " " << v2.pz << std::endl;
  //}
  return v1.px*v2.px + v1.py*v2.py + v1.pz*v2.pz;
}

/// cross product for the spatial 3-vect
/// \param v1    first 4-vect
/// \param v2    second 4-vect
inline CSph3vector cross_product3(const CSph3vector &v1, const CSph3vector &v2){
  //CSph3vector tmp;
  //tmp.px = v1.py*v2.pz-v1.pz*v2.py;
  //tmp.py = v1.pz*v2.px-v1.px*v2.pz;
  //tmp.pz = v1.px*v2.py-v1.py*v2.px;
  //return tmp;
  return CSph3vector(v1.py*v2.pz-v1.pz*v2.py,
		  v1.pz*v2.px-v1.px*v2.pz,
		  v1.px*v2.py-v1.py*v2.px);
}

/// squared norm of the cross product for the spatial 3-vect (energy is set to 0)
/// \param v1    first 4-vect
/// \param v2    second 4-vect
inline double norm2_cross_product3(const CSph3vector &v1, const CSph3vector &v2){
  return sqr(v1.py*v2.pz-v1.pz*v2.py) + sqr(v1.pz*v2.px-v1.px*v2.pz) + sqr(v1.px*v2.py-v1.py*v2.px);
}

/// get tangent squared of the spherical distance between two vectors
/// \param v1    vector defining the first point
/// \param v2    vector defining the second point
inline double get_tan2_distance(const CSphmomentum &v1, const CSphmomentum &v2){
  return norm2_cross_product3(v1,v2)/sqr(dot_product3(v1,v2));
}

/// get spherical distance between to vectors
/// \param v1    vector defining the first point
/// \param v2    vector defining the second point
inline double get_distance(const CSph3vector *v1, const CSph3vector *v2){
  return atan2(sqrt(norm2_cross_product3(*v1,*v2)), dot_product3(*v1,*v2));
}

/// return true if the two points are distant by less than get spherical distance between two vectors
/// \param v1      vector defining the first point
/// \param v2      vector defining the second point
/// \param tan2R   tangent squared of the max distance
/// WARNING: using the tangent here is dangerous for R>pi/2.
///          this never happens per se for "regular R" but 
///          it may in the vicinity computation as we're using
///          2R there. 
inline bool is_closer(const CSph3vector *v1, const CSph3vector *v2, const double tan2R){
  double dot = dot_product3(*v1,*v2);
  return (dot>=0) && (norm2_cross_product3(*v1,*v2)<=tan2R*dot*dot);
}

/// return true if the two points are distant by less than  get spherical distance between to vectors
/// \param v1      vector defining the first point
/// \param v2      vector defining the second point
/// \param tan2R   tangent squared of the max distance
/// safer version but computes the norm
inline bool is_closer_safer(const CSph3vector *v1, const CSph3vector *v2, const double cosR){
  return dot_product3(*v1,*v2)>=cosR*sqrt(v1->norm2()*v2->norm2());
  //double dot = dot_product3(*v1,*v2);
  //return (dot>=0) && (norm2_cross_product3(*v1,*v2)<tan2R*dot*dot);
}

/// multiply a vector by a constant
/// WARNING: norm not updated
inline CSph3vector operator * (const double &r, const CSph3vector &v){
  CSph3vector tmp = v;
  return tmp*=r;
}
}
#endif
