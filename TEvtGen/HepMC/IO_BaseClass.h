//--------------------------------------------------------------------------
#ifndef HEPMC_IO_BASECLASS_H
#define HEPMC_IO_BASECLASS_H

//////////////////////////////////////////////////////////////////////////
// Matt.Dobbs@Cern.CH, November 1999, refer to:
// M. Dobbs and J.B. Hansen, "The HepMC C++ Monte Carlo Event Record for
// High Energy Physics", Computer Physics Communications (to be published).
//
// event input/output base class
//////////////////////////////////////////////////////////////////////////
//
// class from which all input/output classes shall inherit from.
// i.e.: if you want to write events to hbook ntuples,
//              then inherit from this class and re-define read_event()
//              and write_event()
//
// (Possible extension: Could make this an input iterator)
//

#include <iostream>
#include "HepMC/GenEvent.h"

namespace HepMC {

    //! all input/output classes inherit from IO_BaseClass

    ///
    /// \class  IO_BaseClass
    /// If you want to write a new IO class, 
    /// then inherit from this class and re-define read_event()
    /// and write_event()
    ///
    class IO_BaseClass {
    public:
        virtual ~IO_BaseClass() {}

	/// write this GenEvent
	virtual void write_event( const GenEvent* ) =0;
	/// fill this GenEvent
	virtual bool fill_next_event( GenEvent* ) =0;
	/// write output to ostr
	virtual void print( std::ostream& ostr = std::cout ) const;
	//
	// the read_next_event() differs from
	// the fill_***() methods in that it creates a new event
	// before calling the  corresponding fill_*** method
	// (they are not intended to be over-ridden)
	GenEvent*    read_next_event();  //!< do not over-ride
	//
	// The overloaded stream operators >>,<< are identical to
	//   read_next_event and write_event methods respectively.
	//   (or read_particle_data_table and write_particle_data_table)
	// the event argument for the overloaded stream operators is a pointer,
	// which is passed by reference.
	//  i.e.  GenEvent* evt; 
	//        io >> evt; 
	// will give the expected result.
	// (note: I don't see any reason to have separate const and non-const
	//  versions of operator<<, but the pedantic ansi standard insists 
	//  on it) 
	/// the same as read_next_event
	virtual       GenEvent*& operator>>( GenEvent*& );
	/// the same as write_event
	virtual const GenEvent*& operator<<( const GenEvent*& );
	/// the same as write_event
	virtual       GenEvent*& operator<<( GenEvent*& );
    };

    //////////////
    // Inlines  //
    //////////////

    inline GenEvent* IO_BaseClass::read_next_event() {
	/// creates a new event and fills it by calling 
	/// the sister method read_next_event( GenEvent* )
	// 
        // 1. create an empty event container
        GenEvent* evt = new GenEvent();
	// 2. fill the evt container - if the read is successful, return the
	//    pointer, otherwise return null and delete the evt
	if ( fill_next_event( evt ) ) return evt;
	// note: the below delete is only reached if read fails
	//       ... thus there is not much overhead in new then delete 
	//       since this statement is rarely reached
	delete evt;
	return 0;
    }

    inline void IO_BaseClass::print( std::ostream& ostr ) const { 
	ostr << "IO_BaseClass: abstract parent I/O class. " <<  std::endl;
    }

    inline GenEvent*& IO_BaseClass::operator>>( GenEvent*& evt ){
	evt = read_next_event();
	return evt;
    }

    inline const GenEvent*& IO_BaseClass::operator<<(
					      const GenEvent*& evt ) {
	write_event( evt );
	return evt;
    }

    inline GenEvent*& IO_BaseClass::operator<<( GenEvent*& evt ) {
	write_event( evt );
	return evt;
    }

} // HepMC

#endif  // HEPMC_IO_BASECLASS_H
//--------------------------------------------------------------------------



