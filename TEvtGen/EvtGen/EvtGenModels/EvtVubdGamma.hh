//-----------------------------------------------------------------------
// File and Version Information:
//
// Copyright Information: See EvtGen/COPYRIGHT
//
//
// Description:
//    3                                         2                2
//   d Gamma                 /   _      _ _2  mb           _2  mb
//  ---------- = 12 Gamma   | (1+x-z)(z-x-p ) -- W  + (1-z+p ) -- W 
//         _ 2           0   \                 2  1             2  2
//  dx dz dp                                   2
//                                _   _  _2  mb                 2   \.
//                             + [x(z-x)-p ] -- (W + 2mb W  + mb W ) |
//                                            4   3       4       5 /
//
//   with 
//        2 E           2
//           l    _2   p        2 v.p    _
//   x = ------ , p = --- , z = ------ , x = 1-x
//         mb           2         mb
//                    mb 
//
//   the triple differential decay rate according to
//   hep-ph/9905351 v2
//   
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Sven Menke
//
//-----------------------------------------------------------------------

#ifndef EVTVUBDGAMMA_HH
#define EVTVUBDGAMMA_HH

class EvtVubdGamma {

public:
  
  // Constructors

  EvtVubdGamma(const double &alphas);

  // Destructor

  virtual ~EvtVubdGamma( );

  // Operators

  // Selectors 

  // Modifiers

  // Methods

  double getdGdxdzdp(const double &x, const double &z, const double &p2);

protected:
  
  // Helper functions

  double delta(const double &x, const double &xmin, const double &xmax);

  double getW1nodelta(const double &x, const double &z, const double &p2);

  double getW2nodelta(const double &x, const double &z, const double &p2);

  double getW3nodelta(const double &x, const double &z, const double &p2);

  double getW4nodelta(const double &x, const double &z, const double &p2);

  double getW5nodelta(const double &x, const double &z, const double &p2);

  double getW1delta(const double &x, const double &z);
	      					    
  double getW4plus5delta(const double &x, const double &z);
	      					    
private:

  // Friends
  
  // Data members

  double _alphas;
  double _epsilon1;
  double _epsilon2;
  double _epsilon3;
};


#endif // EVTVUBDGAMMA_HH


