//////////////////////////////////////////////////////////////////////////
// Matt.Dobbs@Cern.CH, September 1999
// GenVertex within an event
//////////////////////////////////////////////////////////////////////////

#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenEvent.h"
#include "HepMC/SearchVector.h"
#include <iomanip>       // needed for formatted output

namespace HepMC {

    GenVertex::GenVertex( const FourVector& position,
			  int id, const WeightContainer& weights ) 
	: m_position(position), m_id(id), m_weights(weights), m_event(0),
	  m_barcode(0)
    {}
    //{
	//s_counter++;
    //}

    GenVertex::GenVertex( const GenVertex& invertex ) 
    : m_position( invertex.position() ),
      m_particles_in(),
      m_particles_out(),
      m_id( invertex.id() ),
      m_weights( invertex.weights() ),
      m_event(0),
      m_barcode(0) 
    {
	/// Shallow copy: does not copy the FULL list of particle pointers.
	/// Creates a copy of  - invertex
	///                    - outgoing particles of invertex, but sets the
	///                      decay vertex of these particles to NULL
	///                    - all incoming particles which do not have a
	///                      creation vertex.
	/// (i.e. it creates copies of all particles which it owns)
	/// (note - impossible to copy the FULL list of particle pointers 
	///         while having the vertex
	///         and particles in/out point-back to one another -- unless you
	///         copy the entire tree -- which we don't want to do)
	//
	for ( particles_in_const_iterator 
		  part1 = invertex.particles_in_const_begin();
	      part1 != invertex.particles_in_const_end(); ++part1 ) {
	    if ( !(*part1)->production_vertex() ) { 
		GenParticle* pin = new GenParticle(**part1);
		add_particle_in(pin);
	    }
	}
	for ( particles_out_const_iterator
		  part2 = invertex.particles_out_const_begin();
	      part2 != invertex.particles_out_const_end(); part2++ ) {
	    GenParticle* pin = new GenParticle(**part2);
	    add_particle_out(pin);
	}
	suggest_barcode( invertex.barcode() );
	//
	//s_counter++;
    }
    
    GenVertex::~GenVertex() {
	//
	// need to delete any particles previously owned by this vertex
	if ( parent_event() ) parent_event()->remove_barcode(this);
	delete_adopted_particles();
	//s_counter--;
    }

    void GenVertex::swap( GenVertex & other)
    {
        m_position.swap( other.m_position );
	m_particles_in.swap( other.m_particles_in );
	m_particles_out.swap( other.m_particles_out );
	std::swap( m_id, other.m_id );
	m_weights.swap( other.m_weights );
	std::swap( m_event, other.m_event );
	std::swap( m_barcode, other.m_barcode );
    }

    GenVertex& GenVertex::operator=( const GenVertex& invertex ) {
	/// Shallow: does not copy the FULL list of particle pointers.
	/// Creates a copy of  - invertex
	///                    - outgoing particles of invertex, but sets the
	///                      decay vertex of these particles to NULL
	///                    - all incoming particles which do not have a
	///                      creation vertex.
	///                    - it does not alter *this's m_event (!)
	/// (i.e. it creates copies of all particles which it owns)
	/// (note - impossible to copy the FULL list of particle pointers 
	///         while having the vertex
	///         and particles in/out point-back to one another -- unless you
	///         copy the entire tree -- which we don't want to do)
	///

        // best practices implementation
	GenVertex tmp( invertex );
	swap( tmp );
	return *this;
    }
    
    bool GenVertex::operator==( const GenVertex& a ) const {
	/// Returns true if the positions and the particles in the lists of a 
	///  and this are identical. Does not compare barcodes.
	/// Note that it is impossible for two vertices to point to the same 
	///  particle's address, so we need to do more than just compare the
	///  particle pointers
	//
	if ( a.position() !=  this->position() ) return false;
	// if the size of the inlist differs, return false.
	if ( a.particles_in_size() !=  this->particles_in_size() ) return false;
	// if the size of the outlist differs, return false.
	if ( a.particles_out_size() !=  this->particles_out_size() ) return false;
	// loop over the inlist and ensure particles are identical
	//   (only do this if the lists aren't empty - we already know
	//   if one isn't, then other isn't either!)
	if ( a.particles_in_const_begin() != a.particles_in_const_end() ) {
	    for ( GenVertex::particles_in_const_iterator 
		      ia = a.particles_in_const_begin(),
		      ib = this->particles_in_const_begin();
		  ia != a.particles_in_const_end(); ia++, ib++ ){
		if ( **ia != **ib ) return false;
	    }
	}
	// loop over the outlist and ensure particles are identical
	//   (only do this if the lists aren't empty - we already know
	//   if one isn't, then other isn't either!)
	if ( a.particles_out_const_begin() != a.particles_out_const_end() ) {
	    for ( GenVertex::particles_out_const_iterator 
		      ia = a.particles_out_const_begin(),
		      ib = this->particles_out_const_begin();
		  ia != a.particles_out_const_end(); ia++, ib++ ){
		if ( **ia != **ib ) return false;
	    }
	}
	return true;
    }

    bool GenVertex::operator!=( const GenVertex& a ) const {
	// Returns true if the positions and lists of a and this are not equal.
	return !( a == *this );
    }  

    void GenVertex::print( std::ostream& ostr ) const {
        // find the current stream state
	std::ios_base::fmtflags orig = ostr.flags();
 	std::streamsize prec = ostr.precision();
	if ( barcode()!=0 ) {
	    if ( position() != FourVector(0,0,0,0) ) {
	        ostr << "Vertex:";
		ostr.width(9);
		ostr << barcode();
		ostr << " ID:";
		ostr.width(5);
		ostr << id();
		ostr << " (X,cT)=";
		ostr.width(9);
                ostr.precision(2);
                ostr.setf(std::ios::scientific, std::ios::floatfield);
		ostr.setf(std::ios_base::showpos);
		ostr << position().x() << ",";
		ostr.width(9);
		ostr << position().y() << ",";
		ostr.width(9);
		ostr << position().z() << ",";
		ostr.width(9);
		ostr << position().t();
                ostr.setf(std::ios::fmtflags(0), std::ios::floatfield);
		ostr.unsetf(std::ios_base::showpos);
	        ostr << std::endl;
	    } else {
		ostr << "GenVertex:";
		ostr.width(9);
		ostr << barcode();
		ostr << " ID:";
		ostr.width(5);
		ostr << id();
		ostr << " (X,cT):0";
	        ostr << std::endl;
	    }
	} else {
	    // If the vertex doesn't have a unique barcode assigned, then
	    //  we print its memory address instead... so that the
	    //  print out gives us a unique tag for the particle.
	    if ( position() != FourVector(0,0,0,0) ) {
	        ostr << "Vertex:";
		ostr.width(9);
		ostr << (void*)this;
		ostr << " ID:";
		ostr.width(5);
		ostr << id();
		ostr << " (X,cT)=";
		ostr.width(9);
                ostr.precision(2);
                ostr.setf(std::ios::scientific, std::ios::floatfield);
		ostr.setf(std::ios_base::showpos);
		ostr << position().x();
		ostr.width(9);
		ostr << position().y();
		ostr.width(9);
		ostr << position().z();
		ostr.width(9);
		ostr << position().t();
                ostr.setf(std::ios::fmtflags(0), std::ios::floatfield);
		ostr.unsetf(std::ios_base::showpos);
	        ostr << std::endl;
	    } else {
	        ostr << "GenVertex:";
		ostr.width(9);
		ostr << (void*)this;
		ostr << " ID:";
		ostr.width(5);
		ostr << id();
		ostr << " (X,cT):0";
	        ostr << std::endl;
	    }
	}

	// print the weights if there are any
	if ( ! weights().empty() ) {
	    ostr << " Wgts(" << weights().size() << ")=";
	    for ( WeightContainer::const_iterator wgt = weights().begin();
		  wgt != weights().end(); wgt++ ) { ostr << *wgt << " "; }
	    ostr << std::endl;
	}
	// print out all the incoming, then outgoing particles
	for ( particles_in_const_iterator part1 = particles_in_const_begin();
	      part1 != particles_in_const_end(); part1++ ) {
	    if ( part1 == particles_in_const_begin() ) {
		ostr << " I:";
		ostr.width(2);
		ostr << m_particles_in.size();
	    } else { ostr << "     "; }
	    //(*part1)->print( ostr );  //uncomment for long debugging printout
	    ostr << **part1 << std::endl;
	}
	for ( particles_out_const_iterator part2 = particles_out_const_begin();
	      part2 != particles_out_const_end(); part2++ ) {
	    if ( part2 == particles_out_const_begin() ) { 
		ostr << " O:";
		ostr.width(2);
		ostr << m_particles_out.size();
	    } else { ostr << "     "; }
	    //(*part2)->print( ostr ); // uncomment for long debugging printout
	    ostr << **part2 << std::endl;
	}
        // restore the stream state
        ostr.flags(orig);
        ostr.precision(prec);
    }

    double GenVertex::check_momentum_conservation() const {
	/// finds the difference between the total momentum out and the total
	/// momentum in vectors, and returns the magnitude of this vector
	/// i.e.         returns | vec{p_in} - vec{p_out} |
	double sumpx = 0, sumpy = 0, sumpz = 0;
	for ( particles_in_const_iterator part1 = particles_in_const_begin();
	      part1 != particles_in_const_end(); part1++ ) {
	    sumpx   += (*part1)->momentum().px();
	    sumpy   += (*part1)->momentum().py();
	    sumpz   += (*part1)->momentum().pz();
	}
	for ( particles_out_const_iterator part2 = particles_out_const_begin();
	      part2 != particles_out_const_end(); part2++ ) {
	    sumpx   -= (*part2)->momentum().px();
	    sumpy   -= (*part2)->momentum().py();
	    sumpz   -= (*part2)->momentum().pz();
	}
	return sqrt( sumpx*sumpx + sumpy*sumpy + sumpz*sumpz );
    }

    void GenVertex::add_particle_in( GenParticle* inparticle ) {
	if ( !inparticle ) return;
	// if inparticle previously had a decay vertex, remove it from that
	// vertex's list
	if ( inparticle->end_vertex() ) {
	    inparticle->end_vertex()->remove_particle_in( inparticle );
	}
	m_particles_in.push_back( inparticle );
	inparticle->set_end_vertex_( this );
    }

    void GenVertex::add_particle_out( GenParticle* outparticle ) {
	if ( !outparticle ) return;
	// if outparticle previously had a production vertex,
	// remove it from that vertex's list
	if ( outparticle->production_vertex() ) {
	    outparticle->production_vertex()->remove_particle_out( outparticle );
	}
	m_particles_out.push_back( outparticle );
	outparticle->set_production_vertex_( this );
    }

    GenParticle* GenVertex::remove_particle( GenParticle* particle ) {
	/// this finds *particle in the in and/or out list and removes it from
	///  these lists ... it DOES NOT DELETE THE PARTICLE or its relations.
	/// you could delete the particle too as follows:
	///      delete vtx->remove_particle( particle );
	/// or if the particle has an end vertex, you could:
	///      delete vtx->remove_particle( particle )->end_vertex();
	/// which would delete the particle's end vertex, and thus would
	/// also delete the particle, since the particle would be 
	/// owned by the end vertex.
	if ( !particle ) return 0;
	if ( particle->end_vertex() == this ) {
	    particle->set_end_vertex_( 0 );
	    remove_particle_in(particle);
	}
	if ( particle->production_vertex() == this ) {
	    particle->set_production_vertex_(0);
	    remove_particle_out(particle);
	}
	return particle;
    }

    void GenVertex::remove_particle_in( GenParticle* particle ) {
	/// this finds *particle in m_particles_in and removes it from that list
	if ( !particle ) return;
	m_particles_in.erase( already_in_vector( &m_particles_in, particle ) );
    }

    void GenVertex::remove_particle_out( GenParticle* particle ) {
	/// this finds *particle in m_particles_out and removes it from that list
	if ( !particle ) return;
	m_particles_out.erase( already_in_vector( &m_particles_out, particle ) );
    }

    void GenVertex::delete_adopted_particles() {
	/// deletes all particles which this vertex owns
	/// to be used by the vertex destructor and operator=
	//
	if ( m_particles_out.empty() && m_particles_in.empty() ) return;
	// 1. delete all outgoing particles which don't have decay vertices.
	//    those that do become the responsibility of the decay vertex
	//    and have their productionvertex pointer set to NULL
	for ( std::vector<GenParticle*>::iterator part1 = m_particles_out.begin();
	      part1 != m_particles_out.end(); ) {
	    if ( !(*part1)->end_vertex() ) {
		delete *(part1++);
	    } else { 
		(*part1)->set_production_vertex_(0);
		++part1;
	    }
	}
	m_particles_out.clear();
	//
	// 2. delete all incoming particles which don't have production
	//    vertices. those that do become the responsibility of the 
	//    production vertex and have their decayvertex pointer set to NULL
	for ( std::vector<GenParticle*>::iterator part2 = m_particles_in.begin();
	      part2 != m_particles_in.end(); ) {
	    if ( !(*part2)->production_vertex() ) { 
		delete *(part2++);
	    } else { 
		(*part2)->set_end_vertex_(0); 
		++part2;
	    }
	}
	m_particles_in.clear();
    }

    bool GenVertex::suggest_barcode( int the_bar_code )
    {
	/// allows a barcode to be suggested for this vertex.
	/// In general it is better to let the event pick the barcode for
	/// you, which is automatic.
	/// Returns TRUE if the suggested barcode has been accepted (i.e. the
	///  suggested barcode has not already been used in the event, 
	///  and so it was used).
	/// Returns FALSE if the suggested barcode was rejected, or if the
	///  vertex is not yet part of an event, such that it is not yet
	///  possible to know if the suggested barcode will be accepted).
	if ( the_bar_code >0 ) {
	    std::cerr << "GenVertex::suggest_barcode WARNING, vertex bar codes"
		      << "\n MUST be negative integers. Positive integers "
		      << "\n are reserved for particles only. Your suggestion "
		      << "\n has been rejected." << std::endl;
	    return false;
	}
	bool success = false;
	if ( parent_event() ) {
	    success = parent_event()->set_barcode( this, the_bar_code );
	} else { set_barcode_( the_bar_code ); }
	return success;
    }

    void GenVertex::set_parent_event_( GenEvent* new_evt ) 
    { 
	GenEvent* orig_evt = m_event;
	m_event = new_evt; 
	//
	// every time a vertex's parent event changes, the map of barcodes
	//   in the new and old parent event needs to be modified to 
	//   reflect this
	if ( orig_evt != new_evt ) {
	    if (new_evt) new_evt->set_barcode( this, barcode() );
	    if (orig_evt) orig_evt->remove_barcode( this );
	    // we also need to loop over all the particles which are owned by 
	    //  this vertex, and remove their barcodes from the old event.
	    for ( particles_in_const_iterator part1=particles_in_const_begin();
		  part1 != particles_in_const_end(); part1++ ) {
		if ( !(*part1)->production_vertex() ) { 
		    if ( orig_evt ) orig_evt->remove_barcode( *part1 );
		    if ( new_evt ) new_evt->set_barcode( *part1, 
							 (*part1)->barcode() );
		}
	    }
	    for ( particles_out_const_iterator
		      part2 = particles_out_const_begin();
		  part2 != particles_out_const_end(); part2++ ) {
		    if ( orig_evt ) orig_evt->remove_barcode( *part2 );
		    if ( new_evt ) new_evt->set_barcode( *part2, 
							 (*part2)->barcode() );
	    }
	}
    }

    void GenVertex::change_parent_event_( GenEvent* new_evt ) 
    { 
	//
	// this method is for use with swap
	// particles and vertices have already been exchanged, 
	// but the backpointer needs to be fixed
	//GenEvent* orig_evt = m_event;
	m_event = new_evt; 
    }

    /////////////
    // Static  //
    /////////////
    //unsigned int GenVertex::counter() { return s_counter; }
    //unsigned int GenVertex::s_counter = 0; 

    /////////////
    // Friends //
    /////////////

    /// send vertex information to ostr for printing
    std::ostream& operator<<( std::ostream& ostr, const GenVertex& vtx ) {
	if ( vtx.barcode()!=0 ) ostr << "BarCode " << vtx.barcode();
	else ostr << "Address " << &vtx;
	ostr << " (X,cT)=";
	if ( vtx.position() != FourVector(0,0,0,0)) {
	    ostr << vtx.position().x() << ","
		 << vtx.position().y() << ","
		 << vtx.position().z() << ","
		 << vtx.position().t();
	} else { ostr << 0; }
	ostr << " #in:" << vtx.particles_in_size()
	     << " #out:" << vtx.particles_out_size();
	return ostr;
    }

    /////////////////////////////
    // edge_iterator           // (protected - for internal use only)
    /////////////////////////////
    // If the user wants the functionality of the edge_iterator, he should
    // use particle_iterator with IteratorRange = family, parents, or children
    //

    GenVertex::edge_iterator::edge_iterator() : m_vertex(0), m_range(family), 
	m_is_inparticle_iter(false), m_is_past_end(true)
    {}

    GenVertex::edge_iterator::edge_iterator( const GenVertex& vtx, 
					  IteratorRange range ) :
	m_vertex(&vtx), m_range(family) 
    {
	// Note: (26.1.2000) the original version of edge_iterator inheritted
	//       from set<GenParticle*>::const_iterator() rather than using
	//       composition as it does now.
	//       The inheritted version suffered from a strange bug, which
	//       I have not fully understood --- it only occurred after many
	//       events were processed and only when I called the delete 
	//       function on past events. I believe it had something to do with
	//       the past the end values, which are now robustly coded in this
	//       version as boolean members.
	//
	// default range is family, only other choices are children/parents
	//    descendants/ancestors not allowed & recasted ot children/parents
	if ( range == descendants || range == children ) m_range = children;
	if ( range == ancestors   || range == parents  ) m_range = parents;
	//
	if ( m_vertex->m_particles_in.empty() &&
	     m_vertex->m_particles_out.empty() ) {
	    // Case: particles_in and particles_out is empty.
	    m_is_inparticle_iter = false;
	    m_is_past_end = true;
	} else if ( m_range == parents && m_vertex->m_particles_in.empty() ){
	    // Case: particles in is empty and parents is requested.
	    m_is_inparticle_iter = true;
	    m_is_past_end = true;
	} else if ( m_range == children && m_vertex->m_particles_out.empty() ){
	    // Case: particles out is empty and children is requested.
	    m_is_inparticle_iter = false;
	    m_is_past_end = true;
	} else if ( m_range == children ) {
	    // Case: particles out is NOT empty, and children is requested
	    m_set_iter = m_vertex->m_particles_out.begin();
	    m_is_inparticle_iter = false;
	    m_is_past_end = false;
	} else if ( m_range == family && m_vertex->m_particles_in.empty() ) {
	    // Case: particles in is empty, particles out is NOT empty,
	    //       and family is requested. Then skip ahead to partilces out.
	    m_set_iter = m_vertex->m_particles_out.begin();
	    m_is_inparticle_iter = false;
	    m_is_past_end = false;
	} else {
	    // Normal scenario: start with the first incoming particle
	    m_set_iter = m_vertex->m_particles_in.begin();
	    m_is_inparticle_iter = true;
	    m_is_past_end = false;
	}
    }

    GenVertex::edge_iterator::edge_iterator( const edge_iterator& p ) {
	*this = p;
    }

    GenVertex::edge_iterator::~edge_iterator() {}

    GenVertex::edge_iterator& GenVertex::edge_iterator::operator=(
	const edge_iterator& p ) {
	m_vertex = p.m_vertex;
	m_range = p.m_range;
	m_set_iter = p.m_set_iter;
	m_is_inparticle_iter = p.m_is_inparticle_iter;
	m_is_past_end = p.m_is_past_end;
	return *this;
    }

    GenParticle* GenVertex::edge_iterator::operator*(void) const {
	if ( !m_vertex || m_is_past_end ) return 0;
	return *m_set_iter;
    }

    GenVertex::edge_iterator& GenVertex::edge_iterator::operator++(void){
	// Pre-fix increment
	//
	// increment the set iterator (unless we're past the end value)
	if ( m_is_past_end ) return *this;
	++m_set_iter;
	// handle cases where m_set_iter points past the end
	if ( m_range == family && m_is_inparticle_iter &&
	     m_set_iter == m_vertex->m_particles_in.end() ) {
	    // at the end on in particle set, and range is family, so move to
	    // out particle set
	    m_set_iter = m_vertex->m_particles_out.begin();
	    m_is_inparticle_iter = false;
	} else if ( m_range == parents && 
		    m_set_iter == m_vertex->m_particles_in.end() ) {
	    // at the end on in particle set, and range is parents only, so
	    // move into past the end state
	    m_is_past_end = true;
	    // might as well bail out now
	    return *this;
	} 
	// are we iterating over input or output particles?
	if( m_is_inparticle_iter ) {
	    // the following is not else if because we might have range=family
	    // with an empty particles_in set.	
	    if ( m_set_iter == m_vertex->m_particles_in.end() ) {
		//whenever out particles end is reached, go into past the end state
		m_is_past_end = true;
	    }
	} else {
	    // the following is not else if because we might have range=family
	    // with an empty particles_out set.	
	    if ( m_set_iter == m_vertex->m_particles_out.end() ) {
		//whenever out particles end is reached, go into past the end state
		m_is_past_end = true;
	    }
	}
	return *this;
    }

    GenVertex::edge_iterator GenVertex::edge_iterator::operator++(int){
	// Post-fix increment
	edge_iterator returnvalue = *this;
	++*this;
	return returnvalue;
    }

    bool GenVertex::edge_iterator::is_parent() const {
	if ( **this && (**this)->end_vertex() == m_vertex ) return true;
	return false;
    }

    bool GenVertex::edge_iterator::is_child() const {
	if ( **this && (**this)->production_vertex() == m_vertex ) return true;
	return false;
    }

    int GenVertex::edges_size( IteratorRange range ) const {
	if ( range == children ) return m_particles_out.size();
	if ( range == parents ) return m_particles_in.size();
	if ( range == family ) return m_particles_out.size()
				      + m_particles_in.size();
	return 0;
    }

    /////////////////////
    // vertex_iterator //
    /////////////////////
    
    GenVertex::vertex_iterator::vertex_iterator() 
	: m_vertex(0), m_range(), m_visited_vertices(0), m_it_owns_set(0), 
	  m_recursive_iterator(0)
    {}

    GenVertex::vertex_iterator::vertex_iterator( GenVertex& vtx_root,
					      IteratorRange range )
	: m_vertex(&vtx_root), m_range(range) 
    {
	// standard public constructor
	//
	m_visited_vertices = new std::set<const GenVertex*>;
	m_it_owns_set = 1;
	m_visited_vertices->insert( m_vertex );
	m_recursive_iterator = 0;
	m_edge = m_vertex->edges_begin( m_range );
	// advance to the first good return value
	if ( !follow_edge_() &&
	     m_edge != m_vertex->edges_end( m_range )) ++*this;
    }

    GenVertex::vertex_iterator::vertex_iterator( GenVertex& vtx_root,
	IteratorRange range, std::set<const HepMC::GenVertex*>& visited_vertices ) :
	m_vertex(&vtx_root), m_range(range), 
	m_visited_vertices(&visited_vertices), m_it_owns_set(0),
	m_recursive_iterator(0) 
    {
	// This constuctor is only to be called internally by this class
        //   for use with its recursive member pointer (m_recursive_iterator).
        // Note: we do not need to insert m_vertex_root in the vertex - that is
        //  the responsibility of this iterator's mother, which is normally
        //  done just before calling this protected constructor.
	m_edge = m_vertex->edges_begin( m_range );
	// advance to the first good return value
	if ( !follow_edge_() &&
	     m_edge != m_vertex->edges_end( m_range )) ++*this;
     }

    GenVertex::vertex_iterator::vertex_iterator( const vertex_iterator& v_iter)
	: m_vertex(0), m_visited_vertices(0), m_it_owns_set(0), 
	  m_recursive_iterator(0) 
    {
        *this = v_iter;
    }

    GenVertex::vertex_iterator::~vertex_iterator() {
	if ( m_recursive_iterator ) delete m_recursive_iterator;
        if ( m_it_owns_set ) delete m_visited_vertices;
    }

    GenVertex::vertex_iterator& GenVertex::vertex_iterator::operator=( 
	const vertex_iterator& v_iter ) 
    {
	// Note: when copying a vertex_iterator that is NOT the owner
	// of its set container, the pointer to the set is copied. Beware!
	// (see copy_with_own_set() if you want a different set pointed to)
	// In practise the user never needs to worry
	// since such iterators are only intended to be used internally.
	//
	// destruct data member pointers
	if ( m_recursive_iterator ) delete m_recursive_iterator;
	m_recursive_iterator = 0;
	if ( m_it_owns_set ) delete m_visited_vertices;
	m_visited_vertices = 0;
	m_it_owns_set = 0;
	// copy the target vertex_iterator to this iterator 
	m_vertex = v_iter.m_vertex;
	m_range = v_iter.m_range;
	if ( v_iter.m_it_owns_set ) {
	    // i.e. this vertex will own its set if v_iter points to any
	    // vertex set regardless of whether v_iter owns the set or not!
	    m_visited_vertices = 
		new std::set<const GenVertex*>(*v_iter.m_visited_vertices);
	    m_it_owns_set = 1;
	} else {
	    m_visited_vertices = v_iter.m_visited_vertices;
	    m_it_owns_set = 0;
	}
	//
        // Note: m_vertex_root is already included in the set of 
        //  tv_iter.m_visited_vertices, we do not need to insert it.
        //
        m_edge = v_iter.m_edge;
	copy_recursive_iterator_( v_iter.m_recursive_iterator );
        return *this;
    }

    GenVertex* GenVertex::vertex_iterator::operator*(void) const {
	// de-reference operator
	//
	// if this iterator has an iterator_node, then we return the iterator
        // node.
        if ( m_recursive_iterator ) return  **m_recursive_iterator;
	//
	// an iterator can only return its m_vertex -- any other vertex
        //  is returned by means of a recursive iterator_node 
        //  (so this is the only place in the iterator code that a vertex 
        //   is returned!) 
        if ( m_vertex ) return m_vertex;
        return 0;
    }

    GenVertex::vertex_iterator& GenVertex::vertex_iterator::operator++(void) {
        // Pre-fix incremental operator
        //
        // check for "past the end condition" denoted by m_vertex=0
        if ( !m_vertex ) return *this;
	// if at the last edge, move into the "past the end condition"
	if ( m_edge == m_vertex->edges_end( m_range ) ) {
	    m_vertex = 0;
	    return *this;
	}
	// check to see if we need to create a new recursive iterator by
	// following the current edge only if a recursive iterator doesn't
	// already exist. If a new recursive_iterator is created, return it.
	if ( follow_edge_() ) {
              return *this;
	}
	//
        // if a recursive iterator already exists, increment it, and return its
        // value (unless the recursive iterator has null root_vertex [its
        // root vertex is set to null if it has already returned its root] 
        // - in which case we delete it) 
        // and return the vertex pointed to by the edge.
        if ( m_recursive_iterator ) {
            ++(*m_recursive_iterator);
            if ( **m_recursive_iterator ) {
		return *this;
	    } else {
                delete m_recursive_iterator;
                m_recursive_iterator = 0;
	    }
	}
	//
        // increment to the next particle edge
        ++m_edge;
	// if m_edge is at the end, then we have incremented through all
	// edges, and it is time to return m_vertex, which we accomplish
	// by returning *this
	if ( m_edge == m_vertex->edges_end( m_range ) ) return *this;
	// otherwise we follow the current edge by recursively ++ing.
	return ++(*this);
    }

    GenVertex::vertex_iterator GenVertex::vertex_iterator::operator++(int) {
	// Post-fix increment
        vertex_iterator returnvalue(*this);
        ++(*this);
        return returnvalue;
    }

    void GenVertex::vertex_iterator::copy_with_own_set( 
	const vertex_iterator& v_iter, 
	std::set<const HepMC::GenVertex*>& visited_vertices ) {
	/// intended for internal use only. (use with care!)
	/// this is the same as the operator= method, but it allows the
	/// user to specify which set container m_visited_vertices points to.
	/// in all cases, this vertex will NOT own its set.
	//
	// destruct data member pointers
	if ( m_recursive_iterator ) delete m_recursive_iterator;
	m_recursive_iterator = 0;
	if ( m_it_owns_set ) delete m_visited_vertices;
	m_visited_vertices = 0;
	m_it_owns_set = false;
	// copy the target vertex_iterator to this iterator 
	m_vertex = v_iter.m_vertex;
	m_range = v_iter.m_range;
	m_visited_vertices = &visited_vertices;
	m_it_owns_set = false;
        m_edge = v_iter.m_edge;
	copy_recursive_iterator_( v_iter.m_recursive_iterator );
    }

    GenVertex* GenVertex::vertex_iterator::follow_edge_() {
	// follows the edge pointed to by m_edge by creating a 
	// recursive iterator for it.
	//
	// if a m_recursive_iterator already exists, 
        // this routine has nothing to do,
        // if there's no m_vertex, there's no point following anything, 
	// also there's no point trying to follow a null edge.
	if ( m_recursive_iterator || !m_vertex || !*m_edge ) return 0;
	//
	// if the range is parents, children, or family (i.e. <= family)
	// then only the iterator which owns the set is allowed to create
	// recursive iterators (i.e. recursivity is only allowed to go one
	// layer deep)
	if ( m_range <= family && m_it_owns_set == 0 ) return 0;
	//
	// M.Dobbs 2001-07-16
	// Take care of the very special-rare case where a particle might
	// point to the same vertex for both production and end
	if ( (*m_edge)->production_vertex() == 
	     (*m_edge)->end_vertex() ) return 0;
	//
	// figure out which vertex m_edge is pointing to
	GenVertex* vtx = ( m_edge.is_parent() ? 
			(*m_edge)->production_vertex() :
			(*m_edge)->end_vertex() );
        // if the pointed to vertex doesn't exist or has already been visited, 
        // then return null
	if ( !vtx || !(m_visited_vertices->insert(vtx).second) ) return 0;
	// follow that edge by creating a recursive iterator
	m_recursive_iterator = new vertex_iterator( *vtx, m_range,
						    *m_visited_vertices);
	// and return the vertex pointed to by m_recursive_iterator 
	return **m_recursive_iterator;
    }
	
    void GenVertex::vertex_iterator::copy_recursive_iterator_( 
	const vertex_iterator* recursive_v_iter ) {
	// to properly copy the recursive iterator, we need to ensure
	// the proper set container is transfered ... then do this 
	// operation .... you guessed it .... recursively!
	//
	if ( !recursive_v_iter ) return;
	m_recursive_iterator = new vertex_iterator();
	m_recursive_iterator->m_vertex = recursive_v_iter->m_vertex;
	m_recursive_iterator->m_range = recursive_v_iter->m_range;
	m_recursive_iterator->m_visited_vertices = m_visited_vertices;
	m_recursive_iterator->m_it_owns_set = 0;
	m_recursive_iterator->m_edge = recursive_v_iter->m_edge;
	m_recursive_iterator->copy_recursive_iterator_( 
	    recursive_v_iter->m_recursive_iterator );
    }

    ///////////////////////////////
    // particle_iterator         //
    ///////////////////////////////

    GenVertex::particle_iterator::particle_iterator() {}

    GenVertex::particle_iterator::particle_iterator( GenVertex& vertex_root,
						     IteratorRange range ) {
	// General Purpose Constructor
	//
	if ( range <= family ) {
	    m_edge = GenVertex::edge_iterator( vertex_root, range ); 
	} else {
	    m_vertex_iterator = GenVertex::vertex_iterator(vertex_root, range);
	    m_edge = GenVertex::edge_iterator( **m_vertex_iterator, 
						  m_vertex_iterator.range() ); 
	}
	advance_to_first_();
    }

    GenVertex::particle_iterator::particle_iterator( 
	const particle_iterator& p_iter ){
	*this = p_iter;
    }

    GenVertex::particle_iterator::~particle_iterator() {}

    GenVertex::particle_iterator& 
    GenVertex::particle_iterator::operator=( const particle_iterator& p_iter )
    {
	m_vertex_iterator = p_iter.m_vertex_iterator;
	m_edge = p_iter.m_edge;
	return *this;
    }

    GenParticle* GenVertex::particle_iterator::operator*(void) const {
	return *m_edge;
    }

    GenVertex::particle_iterator& 
    GenVertex::particle_iterator::operator++(void) {
	//Pre-fix increment 
	//
	if ( *m_edge ) {
	    ++m_edge;
	} else if ( *m_vertex_iterator ) {	// !*m_edge is implicit
	    // past end of edge, but still have more vertices to visit
	    // increment the vertex, checking that the result is valid
	    if ( !*(++m_vertex_iterator) ) return *this;
	    m_edge = GenVertex::edge_iterator( **m_vertex_iterator, 
						  m_vertex_iterator.range() ); 
	} else {	// !*m_edge and !*m_vertex_iterator are implicit
	    // past the end condition: do nothing
	    return *this;
	}
	advance_to_first_();
	return *this;
    }

    GenVertex::particle_iterator GenVertex::particle_iterator::operator++(int){
	//Post-fix increment
	particle_iterator returnvalue(*this);
	++(*this);
	return returnvalue;
    }

    GenParticle* GenVertex::particle_iterator::advance_to_first_() {
	/// if the current edge is not a suitable return value ( because
	/// it is a parent of the vertex root that itself belongs to a 
	/// different vertex ) it advances to the first suitable return value 
	if ( !*m_edge ) return *(++*this);
	// if the range is relatives, we need to uniquely assign each particle
	// to a single vertex so as to guarantee particles are returned
	// exactly once.
	if ( m_vertex_iterator.range() == relatives &&
	     m_edge.is_parent() && 
	     (*m_edge)->production_vertex() ) return *(++*this);
	return *m_edge;
    }

    /// scale the position vector
    /// this method is only for use by GenEvent
    /// convert_position assumes that 4th component of the position vector 
    /// is ctau rather than time and has units of length-time
    void GenVertex::convert_position( const double& f ) {
        m_position = FourVector( f*m_position.x(),
                                 f*m_position.y(),
                                 f*m_position.z(),
                                 f*m_position.t() );
   }

} // HepMC
