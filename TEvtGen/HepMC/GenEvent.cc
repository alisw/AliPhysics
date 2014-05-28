//////////////////////////////////////////////////////////////////////////
// Matt.Dobbs@Cern.CH, September 1999
// Updated: 7.1.2000 iterators complete and working!
// Updated: 10.1.2000 GenEvent::vertex, particle iterators are made 
//                    constant WRT this event ... note that 
//                    GenVertex::***_iterator is not const, since it must
//                    be able to return a mutable pointer to itself.
// Updated: 08.2.2000 the event now holds a set of all attached vertices
//                    rather than just the roots of the graph
// Event record for MC generators (for use at any stage of generation)
//////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include "HepMC/GenEvent.h"
#include "HepMC/GenCrossSection.h"
#include "HepMC/Version.h"
#include "HepMC/StreamHelpers.h"

namespace HepMC {

    GenEvent::GenEvent( int signal_process_id, 
                        int event_number,
			GenVertex* signal_vertex,
			const WeightContainer& weights,
			const std::vector<long>& random_states,
			Units::MomentumUnit mom, 
			Units::LengthUnit len ) :
	m_signal_process_id(signal_process_id), 
	m_event_number(event_number),
	m_mpi(-1),
	m_event_scale(-1), 
	m_alphaQCD(-1), 
	m_alphaQED(-1),
	m_signal_process_vertex(signal_vertex), 
	m_beam_particle_1(0),
	m_beam_particle_2(0),
	m_weights(weights),
	m_random_states(random_states),
	m_vertex_barcodes(),
	m_particle_barcodes(),
	m_cross_section(0), 
	m_heavy_ion(0), 
	m_pdf_info(0),
	m_momentum_unit(mom),
	m_position_unit(len)
    {
        /// This constructor only allows null pointers to HeavyIon and PdfInfo
	///
	/// note: default values for m_event_scale, m_alphaQCD, m_alphaQED
	///       are as suggested in hep-ph/0109068, "Generic Interface..."
	///
    }

    GenEvent::GenEvent( int signal_process_id, int event_number,
			GenVertex* signal_vertex,
			const WeightContainer& weights,
			const std::vector<long>& random_states,
			const HeavyIon& ion, 
			const PdfInfo& pdf,
			Units::MomentumUnit mom, 
			Units::LengthUnit len ) :
	m_signal_process_id(signal_process_id), 
	m_event_number(event_number),
	m_mpi(-1),
	m_event_scale(-1), 
	m_alphaQCD(-1), 
	m_alphaQED(-1),
	m_signal_process_vertex(signal_vertex), 
	m_beam_particle_1(0),
	m_beam_particle_2(0),
	m_weights(weights),
	m_random_states(random_states), 
	m_vertex_barcodes(),
	m_particle_barcodes(),
	m_cross_section(0), 
	m_heavy_ion( new HeavyIon(ion) ), 
	m_pdf_info( new PdfInfo(pdf) ),
	m_momentum_unit(mom),
	m_position_unit(len)
    {
        /// GenEvent makes its own copy of HeavyIon and PdfInfo
	///
	/// note: default values for m_event_scale, m_alphaQCD, m_alphaQED
	///       are as suggested in hep-ph/0109068, "Generic Interface..."
    }

    GenEvent::GenEvent( Units::MomentumUnit mom, 
			Units::LengthUnit len, 
			int signal_process_id, 
                        int event_number,
			GenVertex* signal_vertex,
			const WeightContainer& weights,
			const std::vector<long>& random_states ) :
	m_signal_process_id(signal_process_id), 
	m_event_number(event_number),
	m_mpi(-1),
	m_event_scale(-1), 
	m_alphaQCD(-1), 
	m_alphaQED(-1),
	m_signal_process_vertex(signal_vertex), 
	m_beam_particle_1(0),
	m_beam_particle_2(0),
	m_weights(weights),
	m_random_states(random_states),
	m_vertex_barcodes(),
	m_particle_barcodes(),
	m_cross_section(0), 
	m_heavy_ion(0), 
	m_pdf_info(0),
	m_momentum_unit(mom),
	m_position_unit(len)
    {
        /// constructor requiring units - all else is default
        /// This constructor only allows null pointers to HeavyIon and PdfInfo
	///
	/// note: default values for m_event_scale, m_alphaQCD, m_alphaQED
	///       are as suggested in hep-ph/0109068, "Generic Interface..."
	///
    }

    GenEvent::GenEvent( Units::MomentumUnit mom, 
			Units::LengthUnit len,
			int signal_process_id, int event_number,
			GenVertex* signal_vertex,
			const WeightContainer& weights,
			const std::vector<long>& random_states,
			const HeavyIon& ion, 
			const PdfInfo& pdf ) :
	m_signal_process_id(signal_process_id), 
	m_event_number(event_number),
	m_mpi(-1),
	m_event_scale(-1), 
	m_alphaQCD(-1), 
	m_alphaQED(-1),
	m_signal_process_vertex(signal_vertex), 
	m_beam_particle_1(0),
	m_beam_particle_2(0),
	m_weights(weights),
	m_random_states(random_states), 
	m_vertex_barcodes(),
	m_particle_barcodes(),
	m_cross_section(0), 
	m_heavy_ion( new HeavyIon(ion) ), 
	m_pdf_info( new PdfInfo(pdf) ),
	m_momentum_unit(mom),
	m_position_unit(len)
    {
        /// explicit constructor with units first that takes HeavyIon and PdfInfo
        /// GenEvent makes its own copy of HeavyIon and PdfInfo
	///
	/// note: default values for m_event_scale, m_alphaQCD, m_alphaQED
	///       are as suggested in hep-ph/0109068, "Generic Interface..."
    }

    GenEvent::GenEvent( const GenEvent& inevent ) 
      : m_signal_process_id    ( inevent.signal_process_id() ),
	m_event_number         ( inevent.event_number() ),
	m_mpi                  ( inevent.mpi() ),
	m_event_scale          ( inevent.event_scale() ),
	m_alphaQCD             ( inevent.alphaQCD() ),
	m_alphaQED             ( inevent.alphaQED() ),
	m_signal_process_vertex( /* inevent.m_signal_process_vertex */ ),
	m_beam_particle_1      ( /* inevent.m_beam_particle_1 */ ),
	m_beam_particle_2      ( /* inevent.m_beam_particle_2 */ ),
	m_weights              ( /* inevent.m_weights */ ),
	m_random_states        ( /* inevent.m_random_states */ ),
	m_vertex_barcodes      ( /* inevent.m_vertex_barcodes */ ),
	m_particle_barcodes    ( /* inevent.m_particle_barcodes */ ),
	m_cross_section        ( inevent.cross_section() ? new GenCrossSection(*inevent.cross_section()) : 0 ),
	m_heavy_ion            ( inevent.heavy_ion() ? new HeavyIon(*inevent.heavy_ion()) : 0 ),
	m_pdf_info             ( inevent.pdf_info() ? new PdfInfo(*inevent.pdf_info()) : 0 ),
	m_momentum_unit        ( inevent.momentum_unit() ),
	m_position_unit        ( inevent.length_unit() )
    {
	/// deep copy - makes a copy of all vertices!
	//

	// 1. create a NEW copy of all vertices from inevent
	//    taking care to map new vertices onto the vertices being copied
	//    and add these new vertices to this event.
	//    We do not use GenVertex::operator= because that would copy
	//    the attached particles as well.
	std::map<const GenVertex*,GenVertex*> map_in_to_new;
	for ( GenEvent::vertex_const_iterator v = inevent.vertices_begin();
	      v != inevent.vertices_end(); ++v ) {
	    GenVertex* newvertex = new GenVertex(
	        (*v)->position(), (*v)->id(), (*v)->weights() );
            newvertex->suggest_barcode( (*v)->barcode() );
	    map_in_to_new[*v] = newvertex;
	    add_vertex( newvertex );
	}
	// 2. copy the signal process vertex info.
	if ( inevent.signal_process_vertex() ) {
	    set_signal_process_vertex( 
		map_in_to_new[inevent.signal_process_vertex()] );
	} else set_signal_process_vertex( 0 );
        //
        // 3. create a NEW copy of all particles from inevent
        //    taking care to attach them to the appropriate vertex
	GenParticle* beam1(0);
	GenParticle* beam2(0);
        for ( GenEvent::particle_const_iterator p = inevent.particles_begin();
              p != inevent.particles_end(); ++p ) 
        {
	    GenParticle* oldparticle = *p;
	    GenParticle* newparticle = new GenParticle(*oldparticle);
	    if ( oldparticle->end_vertex() ) {
                map_in_to_new[ oldparticle->end_vertex() ]->
                                         add_particle_in(newparticle);
            }
            if ( oldparticle->production_vertex() ) {
                map_in_to_new[ oldparticle->production_vertex() ]->
                                         add_particle_out(newparticle);
            }
	    if ( oldparticle == inevent.beam_particles().first ) beam1 = newparticle;
	    if ( oldparticle == inevent.beam_particles().second ) beam2 = newparticle;
        }
	set_beam_particles( beam1, beam2 );
	//
	// 4. now that vtx/particles are copied, copy weights and random states
	set_random_states( inevent.random_states() );
	weights() = inevent.weights();
    }

    void GenEvent::swap( GenEvent & other )
    {
        // if a container has a swap method, use that for improved performance
	std::swap(m_signal_process_id    , other.m_signal_process_id    );
	std::swap(m_event_number         , other.m_event_number         );
	std::swap(m_mpi                  , other.m_mpi                  );
	std::swap(m_event_scale          , other.m_event_scale          );
	std::swap(m_alphaQCD             , other.m_alphaQCD             );
	std::swap(m_alphaQED             , other.m_alphaQED             );
	std::swap(m_signal_process_vertex, other.m_signal_process_vertex);
	std::swap(m_beam_particle_1      , other.m_beam_particle_1      );
	std::swap(m_beam_particle_2      , other.m_beam_particle_2      );
	m_weights.swap(           other.m_weights  );
	m_random_states.swap(     other.m_random_states  );
	m_vertex_barcodes.swap(   other.m_vertex_barcodes );
	m_particle_barcodes.swap( other.m_particle_barcodes );
	std::swap(m_cross_section        , other.m_cross_section        );
	std::swap(m_heavy_ion            , other.m_heavy_ion            );
	std::swap(m_pdf_info             , other.m_pdf_info             );
	std::swap(m_momentum_unit       , other.m_momentum_unit       );
	std::swap(m_position_unit       , other.m_position_unit       );
	// must now adjust GenVertex back pointers
	for ( GenEvent::vertex_const_iterator vthis = vertices_begin();
	      vthis != vertices_end(); ++vthis ) {
	    (*vthis)->change_parent_event_( this );
	}
	for ( GenEvent::vertex_const_iterator voth = other.vertices_begin();
	      voth != other.vertices_end(); ++voth ) {
	    (*voth)->change_parent_event_( &other );
	}
    }

    GenEvent::~GenEvent() 
    {
	/// Deep destructor.
	/// deletes all vertices/particles in this GenEvent
	/// deletes the associated HeavyIon and PdfInfo
	delete_all_vertices();
	delete m_cross_section;
	delete m_heavy_ion;
	delete m_pdf_info;
    }

    GenEvent& GenEvent::operator=( const GenEvent& inevent ) 
    {
        /// best practices implementation
	GenEvent tmp( inevent );
	swap( tmp );
	return *this;
    }

    void GenEvent::print( std::ostream& ostr ) const {
	/// dumps the content of this event to ostr
	///   to dump to cout use: event.print();
	///   if you want to write this event to file outfile.txt you could use:
	///      std::ofstream outfile("outfile.txt"); event.print( outfile );
	ostr << "________________________________________"
	     << "________________________________________\n";
	ostr << "GenEvent: #" << event_number() 
	     << " ID=" << signal_process_id() 
	     << " SignalProcessGenVertex Barcode: " 
	     << ( signal_process_vertex() ? signal_process_vertex()->barcode()
		  : 0 )
	     << "\n";
	write_units( ostr );
        write_cross_section(ostr);
	ostr << " Entries this event: " << vertices_size() << " vertices, "
	     << particles_size() << " particles.\n"; 
	if( m_beam_particle_1 && m_beam_particle_2 ) {
	  ostr << " Beam Particle barcodes: " << beam_particles().first->barcode() << " "
	       << beam_particles().second->barcode() << " \n";
	} else {
	  ostr << " Beam Particles are not defined.\n";
	}
	// Random State
	ostr << " RndmState(" << m_random_states.size() << ")=";
	for ( std::vector<long>::const_iterator rs 
		  = m_random_states.begin();
	      rs != m_random_states.end(); ++rs ) { ostr << *rs << " "; }
	ostr << "\n";
	// Weights
	ostr << " Wgts(" << weights().size() << ")=";
	weights().print(ostr);
	ostr << " EventScale " << event_scale() 
	     << " [energy] \t alphaQCD=" << alphaQCD() 
	     << "\t alphaQED=" << alphaQED() << std::endl;
	// print a legend to describe the particle info
	ostr << "                                    GenParticle Legend\n";
 	ostr  << "        Barcode   PDG ID      "
	      << "( Px,       Py,       Pz,     E )"
	      << " Stat  DecayVtx\n";
	ostr << "________________________________________"
	     << "________________________________________\n";
	// Print all Vertices
	for ( GenEvent::vertex_const_iterator vtx = this->vertices_begin();
	      vtx != this->vertices_end(); ++vtx ) {
	    (*vtx)->print(ostr); 
	}
	ostr << "________________________________________"
	     << "________________________________________" << std::endl;
    }

    void GenEvent::print_version( std::ostream& ostr ) const {
        ostr << "---------------------------------------------" << std::endl;
        writeVersion( ostr );
        ostr << "---------------------------------------------" << std::endl;
    }

    bool GenEvent::add_vertex( GenVertex* vtx ) {
	/// returns true if successful - generally will only return false
	/// if the inserted vertex is already included in the event.
	if ( !vtx ) return false;
	// if vtx previously pointed to another GenEvent, remove it from that
	// GenEvent's list
	if ( vtx->parent_event() && vtx->parent_event() != this ) {
	    bool remove_status = vtx->parent_event()->remove_vertex( vtx );
	    if ( !remove_status ) {	       
		std::cerr << "GenEvent::add_vertex ERROR "
			  << "GenVertex::parent_event points to \n"
			  << "an event that does not point back to the "
			  << "GenVertex. \n This probably indicates a deeper "
			  << "problem. " << std::endl;
	    }
	}
	//
	// setting the vertex parent also inserts the vertex into this
	// event
	vtx->set_parent_event_( this );
	return ( m_vertex_barcodes.count(vtx->barcode()) ? true : false );
    }

    bool GenEvent::remove_vertex( GenVertex* vtx ) {
	/// this removes vtx from the event but does NOT delete it.
	/// returns True if an entry vtx existed in the table and was erased
	if ( m_signal_process_vertex == vtx ) m_signal_process_vertex = 0;
	if ( vtx->parent_event() == this ) vtx->set_parent_event_( 0 );
	return ( m_vertex_barcodes.count(vtx->barcode()) ? false : true );
    }

    void GenEvent::clear() 
    {
	/// remove all information from the event
	/// deletes all vertices/particles in this evt
	///
	delete_all_vertices();
	// remove existing objects and set pointers to null
	delete m_cross_section;
	m_cross_section = 0;
	delete m_heavy_ion;
	m_heavy_ion = 0;
	delete m_pdf_info;
	m_pdf_info = 0;
	m_signal_process_id = 0;
        m_beam_particle_1 = 0;
	m_beam_particle_2 = 0;
	m_event_number = 0;
	m_mpi = -1;
	m_event_scale = -1;
	m_alphaQCD = -1;
	m_alphaQED = -1;
	m_weights = std::vector<double>();
	m_random_states = std::vector<long>();
	// resetting unit information
	m_momentum_unit = Units::default_momentum_unit();
	m_position_unit = Units::default_length_unit();
        // error check just to be safe
	if ( m_vertex_barcodes.size() != 0 
	     || m_particle_barcodes.size() != 0 ) {
	    std::cerr << "GenEvent::clear() strange result ... \n"
		      << "either the particle and/or the vertex map isn't empty" << std::endl;
            std::cerr << "Number vtx,particle the event after deleting = "
                      << m_vertex_barcodes.size() << "  " 
		      << m_particle_barcodes.size() << std::endl;
        }
	return;
    }
    
    void GenEvent::delete_all_vertices() {
	/// deletes all vertices in the vertex container
	/// (i.e. all vertices owned by this event)
	/// The vertices are the "owners" of the particles, so as we delete
	///   the vertices, the vertex desctructors are automatically
	///   deleting their particles.

  	// delete each vertex individually (this deletes particles as well)
	while ( !vertices_empty() ) {
	    GenVertex* vtx = ( m_vertex_barcodes.begin() )->second;
            m_vertex_barcodes.erase( m_vertex_barcodes.begin() );
            delete vtx;
 	}
	//
	// Error checking:
	if ( !vertices_empty() || ! particles_empty() ) {
	    std::cerr << "GenEvent::delete_all_vertices strange result ... "
		      << "after deleting all vertices, \nthe particle and "
		      << "vertex maps aren't empty.\n  This probably " 
		      << "indicates deeper problems or memory leak in the "
		      << "code." << std::endl;
            std::cerr << "Number vtx,particle the event after deleting = "
                      << m_vertex_barcodes.size() << "  " 
		      << m_particle_barcodes.size() << std::endl;
	}
    }
    
    bool GenEvent::set_barcode( GenParticle* p, int suggested_barcode )
    {
	if ( p->parent_event() != this ) {
	    std::cerr << "GenEvent::set_barcode attempted, but the argument's"
		      << "\n parent_event is not this ... request rejected."
		      << std::endl;
	    return false;
	}
	// M.Dobbs  Nov 4, 2002
	// First we must check to see if the particle already has a
	// barcode which is different from the suggestion. If yes, we
	// remove it from the particle map.
	if ( p->barcode() != 0 && p->barcode() != suggested_barcode ) {
	    if ( m_particle_barcodes.count(p->barcode()) &&
		 m_particle_barcodes[p->barcode()] == p ) {
		m_particle_barcodes.erase( p->barcode() );
	    }
	    // At this point either the particle is NOT in
	    // m_particle_barcodes, or else it is in the map, but
	    // already with the suggested barcode.
	}
	//
	// First case --- a valid barcode has been suggested
	//     (valid barcodes are numbers greater than zero)
	bool insert_success = true;
	if ( suggested_barcode > 0 ) {
	    if ( m_particle_barcodes.count(suggested_barcode) ) {
		// the suggested_barcode is already used.
		if ( m_particle_barcodes[suggested_barcode] == p ) {
		    // but it was used for this particle ... so everythings ok
		    p->set_barcode_( suggested_barcode );
		    return true;
		}
		insert_success = false;
		suggested_barcode = 0;
	    } else { // suggested barcode is OK, proceed to insert
		m_particle_barcodes[suggested_barcode] = p;
		p->set_barcode_( suggested_barcode );
		return true;
	    }
	}
	//
	// Other possibility -- a valid barcode has not been suggested,
	//    so one is automatically generated
	if ( suggested_barcode < 0 ) insert_success = false;
	if ( suggested_barcode <= 0 ) {
	    if ( !m_particle_barcodes.empty() ) {
		// in this case we find the highest barcode that was used,
		// and increment it by 1
		suggested_barcode = m_particle_barcodes.rbegin()->first;
		++suggested_barcode;
	    }
	    // For the automatically assigned barcodes, the first one
	    //   we use is 10001 ... this is just because when we read 
	    //   events from HEPEVT, we will suggest their index as the
	    //   barcode, and that index has maximum value 10000.
	    //  This way we will immediately be able to recognize the hepevt
	    //   particles from the auto-assigned ones.
	    if ( suggested_barcode <= 10000 ) suggested_barcode = 10001;
	}
	// At this point we should have a valid barcode
	if ( m_particle_barcodes.count(suggested_barcode) ) {
	    std::cerr << "GenEvent::set_barcode ERROR, this should never "
		      << "happen \n report bug to matt.dobbs@cern.ch" 
		      << std::endl;
	}
	m_particle_barcodes[suggested_barcode] = p;
	p->set_barcode_( suggested_barcode );
	return insert_success;
    }

    bool GenEvent::set_barcode( GenVertex* v, int suggested_barcode )
    {
	if ( v->parent_event() != this ) {
	    std::cerr << "GenEvent::set_barcode attempted, but the argument's"
		      << "\n parent_event is not this ... request rejected."
		      << std::endl;
	    return false;
	}
	// M.Dobbs Nov 4, 2002
	// First we must check to see if the vertex already has a
	// barcode which is different from the suggestion. If yes, we
	// remove it from the vertex map.
	if ( v->barcode() != 0 && v->barcode() != suggested_barcode ) {
	    if ( m_vertex_barcodes.count(v->barcode()) &&
		 m_vertex_barcodes[v->barcode()] == v ) {
		m_vertex_barcodes.erase( v->barcode() );
	    }
	    // At this point either the vertex is NOT in
	    // m_vertex_barcodes, or else it is in the map, but
	    // already with the suggested barcode.
	}
	
	//
	// First case --- a valid barcode has been suggested
	//     (valid barcodes are numbers greater than zero)
	bool insert_success = true;
	if ( suggested_barcode < 0 ) {
	    if ( m_vertex_barcodes.count(suggested_barcode) ) {
		// the suggested_barcode is already used.
		if ( m_vertex_barcodes[suggested_barcode] == v ) {
		    // but it was used for this vertex ... so everythings ok
		    v->set_barcode_( suggested_barcode );
		    return true;
		}
		insert_success = false;
		suggested_barcode = 0;
	    } else { // suggested barcode is OK, proceed to insert
		m_vertex_barcodes[suggested_barcode] = v;
		v->set_barcode_( suggested_barcode );
		return true;
	    }
	}
	//
	// Other possibility -- a valid barcode has not been suggested,
	//    so one is automatically generated
	if ( suggested_barcode > 0 ) insert_success = false;
	if ( suggested_barcode >= 0 ) {
	    if ( !m_vertex_barcodes.empty() ) {
		// in this case we find the highest barcode that was used,
		// and increment it by 1, (vertex barcodes are negative)
		suggested_barcode = m_vertex_barcodes.rbegin()->first;
		--suggested_barcode;
	    }
	    if ( suggested_barcode >= 0 ) suggested_barcode = -1;
	}
	// At this point we should have a valid barcode
	if ( m_vertex_barcodes.count(suggested_barcode) ) {
	    std::cerr << "GenEvent::set_barcode ERROR, this should never "
		      << "happen \n report bug to matt.dobbs@cern.ch" 
		      << std::endl;
	}
	m_vertex_barcodes[suggested_barcode] = v;
	v->set_barcode_( suggested_barcode );
	return insert_success;
    }

    /// test to see if we have two valid beam particles
    bool  GenEvent::valid_beam_particles() const {
        bool have1 = false;
        bool have2 = false;
	// first check that both are defined
        if(m_beam_particle_1 && m_beam_particle_2) {
	    // now look for a match with the particle "list"
            for ( particle_const_iterator p = particles_begin();
        	  p != particles_end(); ++p ) {
		if( m_beam_particle_1 == *p ) have1 = true;
		if( m_beam_particle_2 == *p ) have2 = true;
	    }
	}
	if( have1 && have2 ) return true;
	return false;
    }
    
    /// construct the beam particle information using pointers to GenParticle
    /// returns false if either GenParticle* is null
    bool  GenEvent::set_beam_particles(GenParticle* bp1, GenParticle* bp2) {
	m_beam_particle_1 = bp1;
	m_beam_particle_2 = bp2;
	if( m_beam_particle_1 && m_beam_particle_2 ) return true;
	return false;
    }

    /// construct the beam particle information using a std::pair of pointers to GenParticle
    /// returns false if either GenParticle* is null
    bool  GenEvent::set_beam_particles(std::pair<HepMC::GenParticle*, HepMC::GenParticle*> const & bp) {
	return set_beam_particles(bp.first,bp.second);
    }

    void GenEvent::write_units( std::ostream & os ) const {
	os << " Momenutm units:" << std::setw(8) << name(momentum_unit());
	os << "     Position units:" << std::setw(8) << name(length_unit());
	os << std::endl;
    }

    void GenEvent::write_cross_section( std::ostream& os ) const
    {
	// write the GenCrossSection information if the cross section was set
	if( !cross_section() ) return;
	if( cross_section()->is_set() ) {
	    os << " Cross Section: " << cross_section()->cross_section() ;
	    os << " +/- " << cross_section()->cross_section_error() ;
	    os << std::endl;
	}
    }

   bool GenEvent::use_momentum_unit( Units::MomentumUnit newunit ) { 
	// currently not exception-safe. 
	// Easy to fix, though, if needed.
	if ( m_momentum_unit != newunit ) { 
	    const double factor = Units::conversion_factor( m_momentum_unit, newunit );
	    // multiply all momenta by 'factor',  
	    // loop is entered only if particle list is not empty
            for ( GenEvent::particle_iterator p = particles_begin();
                                              p != particles_end(); ++p ) 
            {
		(*p)->convert_momentum(factor);
            }
	    // ... 
	    m_momentum_unit = newunit; 
	}
	return true; 
    }
    
    bool GenEvent::use_length_unit( Units::LengthUnit newunit ) { 
	// currently not exception-safe. 
	// Easy to fix, though, if needed.
	if ( m_position_unit != newunit ) { 
	    const double factor = Units::conversion_factor( m_position_unit, newunit );
	    // multiply all lengths by 'factor', 
	    // loop is entered only if vertex list is not empty
	    for ( GenEvent::vertex_iterator vtx = vertices_begin();
	                                    vtx != vertices_end(); ++vtx ) {
		(*vtx)->convert_position(factor);
	    }
	    // ... 
	    m_position_unit = newunit; 
	} 
	return true; 
    }  

    bool GenEvent::use_momentum_unit( std::string& newunit ) { 
        if     ( newunit == "MEV" ) return use_momentum_unit( Units::MEV );
	else if( newunit == "GEV" ) return use_momentum_unit( Units::GEV );
	else std::cerr << "GenEvent::use_momentum_unit ERROR: use either MEV or GEV\n";
	return false;
    }
    
    bool GenEvent::use_length_unit( std::string& newunit ) { 
        if     ( newunit == "MM" ) return use_length_unit( Units::MM );
	else if( newunit == "CM" ) return use_length_unit( Units::CM );
	else std::cerr << "GenEvent::use_length_unit ERROR: use either MM or CM\n";
	return false;
    }  
    
    void GenEvent::define_units( std::string& new_m, std::string& new_l ) { 

        if     ( new_m == "MEV" ) m_momentum_unit = Units::MEV ;
	else if( new_m == "GEV" ) m_momentum_unit = Units::GEV ;
	else std::cerr << "GenEvent::define_units ERROR: use either MEV or GEV\n";

        if     ( new_l == "MM" ) m_position_unit = Units::MM ;
	else if( new_l == "CM" ) m_position_unit = Units::CM ;
	else std::cerr << "GenEvent::define_units ERROR: use either MM or CM\n";

    }

    bool GenEvent::is_valid() const {
        /// A GenEvent is presumed valid if it has both associated
	/// particles and vertices.   No other information is checked.
        if ( vertices_empty() ) return false;
	if ( particles_empty() ) return false;
	return true;
    }

    std::ostream & GenEvent::write_beam_particles(std::ostream & os, 
                	 std::pair<HepMC::GenParticle *,HepMC::GenParticle *> pr )
    {
	GenParticle* p = pr.first;
	if(!p) {
	   detail::output( os, 0 );
	} else {
	   detail::output( os, p->barcode() );
	}
	p = pr.second;
	if(!p) {
	   detail::output( os, 0 );
	} else {
	   detail::output( os, p->barcode() );
	}

	return os;
    }

    std::ostream & GenEvent::write_vertex(std::ostream & os, GenVertex const * v)
    {
	if ( !v || !os ) {
	    std::cerr << "GenEvent::write_vertex !v||!os, "
		      << "v="<< v << " setting badbit" << std::endl;
	    os.clear(std::ios::badbit); 
	    return os;
	}
	// First collect info we need
	// count the number of orphan particles going into v
	int num_orphans_in = 0;
	for ( GenVertex::particles_in_const_iterator p1
		  = v->particles_in_const_begin();
	      p1 != v->particles_in_const_end(); ++p1 ) {
	    if ( !(*p1)->production_vertex() ) ++num_orphans_in;
	}
	//
	os << 'V';
	detail::output( os, v->barcode() ); // v's unique identifier
	detail::output( os, v->id() );
	detail::output( os, v->position().x() );
	detail::output( os, v->position().y() );
	detail::output( os, v->position().z() );
	detail::output( os, v->position().t() );
	detail::output( os, num_orphans_in );
	detail::output( os, (int)v->particles_out_size() );
	detail::output( os, (int)v->weights().size() );
	for ( WeightContainer::const_iterator w = v->weights().begin(); 
	      w != v->weights().end(); ++w ) {
	    detail::output( os, *w );
	}
	detail::output( os,'\n');
	// incoming particles
	for ( GenVertex::particles_in_const_iterator p2 
		  = v->particles_in_const_begin();
	      p2 != v->particles_in_const_end(); ++p2 ) {
	    if ( !(*p2)->production_vertex() ) {
		write_particle( os, *p2 );
	    }
	}
	// outgoing particles
	for ( GenVertex::particles_out_const_iterator p3 
		  = v->particles_out_const_begin();
	      p3 != v->particles_out_const_end(); ++p3 ) {
	    write_particle( os, *p3 );
	}
	return os;
    }

    std::ostream & GenEvent::write_particle( std::ostream & os, GenParticle const * p )
    {
	if ( !p || !os ) {
	    std::cerr << "GenEvent::write_particle !p||!os, "
		      << "p="<< p << " setting badbit" << std::endl;
	    os.clear(std::ios::badbit); 
	    return os;
	}
	os << 'P';
	detail::output( os, p->barcode() );
	detail::output( os, p->pdg_id() );
	detail::output( os, p->momentum().px() );
	detail::output( os, p->momentum().py() );
	detail::output( os, p->momentum().pz() );
	detail::output( os, p->momentum().e() );
	detail::output( os, p->generated_mass() );
	detail::output( os, p->status() );
	detail::output( os, p->polarization().theta() );
	detail::output( os, p->polarization().phi() );
	// since end_vertex is oftentimes null, this CREATES a null vertex
	// in the map
	detail::output( os,   ( p->end_vertex() ? p->end_vertex()->barcode() : 0 )  );
	os << ' ' << p->flow() << "\n";

	return os;
    }

} // HepMC
