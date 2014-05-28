#ifndef HEPMC_UNITS_H
#define HEPMC_UNITS_H

//--------------------------------------------------------------------------
// Units.h
// Author:  A. Buckley, D. Grellscheid
//
// units used by a GenEvent
// The default units are set by a configure switch at compile time in Units.cc.
//--------------------------------------------------------------------------

#include <iostream>
#include <string>

namespace HepMC {

  ///
  /// \namespace Units
  /// Allow units to be specified within HepMC.
  /// The default units are set at compile time. 
  ///
  namespace Units {

    // Convention: if both types are passed, MomentumUnit always goes first.
    enum MomentumUnit { MEV, GEV };	//!< momentum units
    enum LengthUnit   { MM, CM };	//!< position units
    
    LengthUnit   default_length_unit();		//!< default unit is defined by configure
    MomentumUnit default_momentum_unit();	//!< default unit is defined by configure

    // helper functions
    std::string name( MomentumUnit );	//!< convert enum to string
    std::string name( LengthUnit );	//!< convert enum to string

    /// scaling factor relative to MeV
    double conversion_factor( MomentumUnit from, MomentumUnit to ); 
    double conversion_factor( LengthUnit from, LengthUnit to );

  }	// Units
}	// HepMC

#endif // HEPMC_UNITS_H
