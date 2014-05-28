//--------------------------------------------------------------------------
#ifndef HEPMC_IO_HEPEVT_H
#define HEPMC_IO_HEPEVT_H

//////////////////////////////////////////////////////////////////////////
// Matt.Dobbs@Cern.CH, January 2000, refer to:
// M. Dobbs and J.B. Hansen, "The HepMC C++ Monte Carlo Event Record for
// High Energy Physics", Computer Physics Communications (to be published).
//
// HEPEVT IO class
//////////////////////////////////////////////////////////////////////////
//
// Important note: This class uses HepMC::HEPEVT_Wrapper which is an
//                 interface to the fortran77 HEPEVT common block.
//                 The precision and number of entries in the F77 common 
//                 block can be specified. See HepMC/HEPEVT_Wrapper.h.
//                 You will very likely have to specify these values for your
//                 application.
//
//

#include <map>
#include <vector>
#include "HepMC/IO_BaseClass.h"
#include "HepMC/HEPEVT_Wrapper.h"

namespace HepMC {

    class GenEvent;
    class GenVertex;
    class GenParticle;

    //! HEPEVT IO class

    ///
    /// \class  IO_HEPEVT
    /// IO class for reading the standard HEPEVT common block.
    ///
    class IO_HEPEVT : public IO_BaseClass {
    public:
	IO_HEPEVT();
	virtual           ~IO_HEPEVT();
	bool              fill_next_event( GenEvent* );
	void              write_event( const GenEvent* );
	void              print( std::ostream& ostr = std::cout ) const;
	
	// see comments below for these switches.
	/// default is false
	bool              trust_both_mothers_and_daughters() const;
	/// default is true
	bool              trust_mothers_before_daughters() const;
	/// default is true
	bool              print_inconsistency_errors() const;
	/// default is true
	bool              trust_beam_particles() const;
        /// define mother daughter trust rules
	void              set_trust_mothers_before_daughters( bool b = true );
        /// define mother daughter trust rules
	void              set_trust_both_mothers_and_daughters( bool b = false );
	/// Since HEPEVT has bi-directional pointers, it is possible that
	/// the mother/daughter pointers are inconsistent (though physically
	/// speaking this should never happen). In practise it happens often.
	/// When a conflict occurs (i.e. when mother/daughter pointers are in 
	/// disagreement, where an empty (0) pointer is not considered a 
	/// disagreement) an error is printed. These errors can be turned off 
	/// with:            myio_hepevt.set_print_inconsistency_errors(0);
	/// but it is STRONGLY recommended that you print the HEPEVT 
	/// common and understand the inconsistency BEFORE you turn off the
	/// errors. The messages are there for a reason [remember, there is
	/// no message printed when the information is missing, ... only when
	/// is it inconsistent. User beware.]
	/// You can inspect the HEPEVT common block for inconsistencies with
	///   HEPEVT_Wrapper::check_hepevt_consistency()
	///
	/// There is a switch controlling whether the mother pointers or
	/// the daughters are to be trusted.
	/// For example, in Pythia the mother information is always correctly
	/// included, but the daughter information is often left unfilled: in
	/// this case we want to trust the mother pointers and not necessarily
	/// the daughters. [THIS IS THE DEFAULT]. Unfortunately the reverse
	/// happens for the stdhep(2001) translation of Isajet, so we need
	/// an option to toggle the choices.
	void              set_print_inconsistency_errors( bool b = true );
        /// declare whether or not beam particles exist
	void              set_trust_beam_particles( bool b = true );

    protected: // for internal use only
        /// create a GenParticle
	GenParticle* build_particle( int index );
        /// create a production vertex
	void build_production_vertex( 
	    int i,std::vector<HepMC::GenParticle*>& hepevt_particle, GenEvent* evt );
        /// create an end vertex
	void build_end_vertex( 
	    int i, std::vector<HepMC::GenParticle*>& hepevt_particle, GenEvent* evt );
        /// find this particle in the particle map
	int  find_in_map( 
	    const std::map<HepMC::GenParticle*,int>& m, GenParticle* p) const;

    private: // use of copy constructor is not allowed
	IO_HEPEVT( const IO_HEPEVT& ) : IO_BaseClass() {}

    private: // data members

	bool m_trust_mothers_before_daughters;
	bool m_trust_both_mothers_and_daughters;
	bool m_print_inconsistency_errors; 
	bool m_trust_beam_particles;
    };

    ////////////////////////////
    // INLINES access methods //
    ////////////////////////////
    inline bool IO_HEPEVT::trust_both_mothers_and_daughters() const 
    { return m_trust_both_mothers_and_daughters; }
	
    inline bool IO_HEPEVT::trust_mothers_before_daughters() const 
    { return m_trust_mothers_before_daughters; }

    inline bool IO_HEPEVT::print_inconsistency_errors() const
    { return m_print_inconsistency_errors; }

    inline void IO_HEPEVT::set_trust_both_mothers_and_daughters( bool b )
    { m_trust_both_mothers_and_daughters = b; }

    inline void IO_HEPEVT::set_trust_mothers_before_daughters( bool b )
    { m_trust_mothers_before_daughters = b; }

    inline void IO_HEPEVT::set_print_inconsistency_errors( bool b  )
    { m_print_inconsistency_errors = b; }

    inline bool IO_HEPEVT::trust_beam_particles() const
    { return m_trust_beam_particles; }

    inline void IO_HEPEVT::set_trust_beam_particles( bool b )
    { m_trust_beam_particles = b; }

} // HepMC

#endif  // HEPMC_IO_HEPEVT_H
//--------------------------------------------------------------------------
