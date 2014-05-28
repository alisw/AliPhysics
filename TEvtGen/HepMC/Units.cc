//--------------------------------------------------------------------------
// Units.cc
// Author:  A. Buckley, D. Grellscheid
//
// units used by a GenEvent
// The default units are set here at compile time.
//--------------------------------------------------------------------------

#include "HepMC/Units.h"

namespace HepMC {

  namespace Units {

   // helper functions
    std::string name(MomentumUnit m) {
      switch (m) {
        case MEV : return "MEV";
        case GEV : return "GEV";
 	default  : return "badValue";
      } 
    }

    std::string name(LengthUnit l) {
      switch (l) {
        case MM : return "MM";
        case CM : return "CM";
  	default : return "badValue";
     } 
    }

    double conversion_factor(MomentumUnit from, MomentumUnit to) 
    {
      if ( from == to )
	return 1.0;
      else if ( from == MEV && to == GEV )
	return 0.001;
      else
	return 1000.0;
    }

    double conversion_factor(LengthUnit from, LengthUnit to) 
    {
      if ( from == to )
	return 1.0;
      else if ( from == MM && to == CM )
	return 0.1;
      else
	return 10.0;
    }

    // if this function fails to compile, rerun configure using --with-length_units
    LengthUnit default_length_unit() {
      //return @HEPMC_DEFAULT_LEN_UNIT@ ;
      return MM ;
    }

    // if this function fails to compile, rerun configure using --with-momentum_units
    MomentumUnit default_momentum_unit() {
      //return @HEPMC_DEFAULT_MOM_UNIT@ ;
      return MEV ;
    }

  }	// Units

}	// HepMC
