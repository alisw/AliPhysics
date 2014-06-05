//--------------------------------------------------------------------------
#ifndef HEPMC_GEN_PARTICLE_H
#define HEPMC_GEN_PARTICLE_H

//////////////////////////////////////////////////////////////////////////
// Matt.Dobbs@Cern.CH, September 1999, refer to:
// M. Dobbs and J.B. Hansen, "The HepMC C++ Monte Carlo Event Record for
// High Energy Physics", Computer Physics Communications (to be published).
//
// particle within an event coming in/out of a vertex
// particle is the basic building block or unit of the event record
//////////////////////////////////////////////////////////////////////////
//
// example:
//      GenParticle* p = new GenParticle( FourVector(1,1,1,3), 11, 1 );
// creates a particle with 4-vector (p,E)=1,1,1,3 - with pdg id 11 (electron)
// and give this particle status =1.
//
// the pointers to end/production vertices can only be set by the
//  vertices themselves - thus to set the production vertex for a particle,
//  you add the particle to that vertex with GenVertex::add_particle_out()
//
// We decide not to have a separate 4 vector for the momentum 
//  at decay time (which MC++ includes to allow dE/dX losses etc). 
//  If you want that, just add a decay vertex with the
//  same particle (modified momentum) going out
//

#include "HepMC/Flow.h"
#include "HepMC/Polarization.h"
#include "HepMC/SimpleVector.h"
#include "HepMC/IteratorRange.h"
#include <iostream>
#ifdef _WIN32
#define hepmc_uint64_t  __int64
#else
#include <stdint.h>	// for uint64_t
#define hepmc_uint64_t   uint64_t
#endif

namespace HepMC {

    class GenVertex;
    class GenEvent; 

    class GenParticleProductionRange;
    class ConstGenParticleProductionRange;
    class GenParticleEndRange;
    class ConstGenParticleEndRange;

    //! The GenParticle class contains information about generated particles

    ///
    /// \class GenParticle 
    /// HepMC::GenParticle 
    /// contains momentum, generated mass, particle ID, decay status, 
    /// flow, polarization, pointers to production and decay vertices
    /// and a unique barcode identfier.
    ///
    class GenParticle {

	friend class GenVertex; // so vertex can set decay/production vertexes
	friend class GenEvent;  // so event can set the barCodes
	/// print particle
	friend std::ostream& operator<<( std::ostream&, const GenParticle& );

    public:
        /// default constructor
        GenParticle(void);
        /// constructor requires momentum and particle ID
 	GenParticle( const FourVector& momentum, int pdg_id,
		     int status = 0, const Flow& itsflow = Flow(),
		     const Polarization& polar = Polarization(0,0) );
	GenParticle( const GenParticle& inparticle ); //!< shallow copy.
	virtual ~GenParticle();

        void swap( GenParticle & other); //!< swap
	GenParticle& operator=( const GenParticle& inparticle ); //!< shallow.
        /// check for equality
	bool         operator==( const GenParticle& ) const;
        /// check for inequality
	bool         operator!=( const GenParticle& ) const;

	/// dump this particle's full info to ostr
	void       print( std::ostream& ostr = std::cout ) const; 

	operator HepMC::FourVector() const; //!< conversion operator

	////////////////////
	// access methods //
	////////////////////

	/// standard 4 momentum
	const FourVector &          momentum() const;
	/// particle ID
	int                  pdg_id() const;
	/// HEPEVT decay status
	int                  status() const;
	/// particle flow
	const Flow &         flow() const;
	/// particle flow index
	int                  flow( int code_index ) const;
        /// polarization information
	const Polarization & polarization() const;
	/// pointer to the production vertex
	GenVertex*           production_vertex() const;
	/// pointer to the decay vertex
	GenVertex*           end_vertex() const;
	/// pointer to the event that owns this particle
	GenEvent*            parent_event() const;

        /// Because of precision issues, the generated mass is not always the 
	/// same as the mass calculated from the momentum 4 vector.
        /// If the generated mass has been set, then generated_mass() 
	/// returns that value.
        /// If the generated mass has not been set, then generated_mass() 
	/// returns the mass calculated from the momentum 4 vector.
        double               generated_mass() const; //!< mass as generated

	/// generatedMass() is included for backwards compatibility with CLHEP HepMC
        double               generatedMass() const { return generated_mass(); }


	///
	/// The barcode is the particle's reference number, every vertex in the
	/// event has a unique barcode. Particle barcodes are positive numbers,
	/// vertex barcodes are negative numbers.
	/// 
	/// Please note that the barcodes are intended for internal use within 
	/// HepMC as a unique identifier for the particles and vertices.
	/// Using the barcode to encode extra information is an abuse of 
	/// the barcode data member and causes confusion among users. 
	/// 
	int                  barcode() const; //!< particle barcode
	
	/// Convenience method.  Returns true if status==1
	bool                 is_undecayed() const;
	/// Convenience method.  Returns true if status==2
	bool                 has_decayed() const;
	/// Convenience method.  Returns true if status==4
	/// Note that using status 4 for beam particles is a new convention which
	/// may not have been implemented by the code originating this GenEvent.
	bool                 is_beam() const;

	/// incoming particle range
	GenParticleProductionRange particles_in( IteratorRange range = relatives );
	/// incoming particle range
	ConstGenParticleProductionRange particles_in( IteratorRange range = relatives ) const;
	/// outgoing particle range
	GenParticleEndRange particles_out( IteratorRange range = relatives );
	/// outgoing particle range
	ConstGenParticleEndRange particles_out( IteratorRange range = relatives ) const;

	/////////////////////
	// mutator methods //
	/////////////////////

	/// In general there is no reason to "suggest_barcode"
	bool                 suggest_barcode( int the_bar_code );

	void   set_momentum( const FourVector& vec4 ); //!< set standard 4 momentum
	void   set_pdg_id( int id ); //!< set particle ID
	void   set_status( int status = 0 ); //!< set decay status
	void   set_flow( const Flow& f ); //!< set particle flow
	void   set_flow( int code_index, int code = 0 ); //!< set particle flow index
	/// set polarization
	void   set_polarization( const Polarization& pol = Polarization(0,0) );
        ///  If you do not call set_generated_mass(), then 
        ///  generated_mass() will simply return the mass calculated from momentum()
        void   set_generated_mass( const double & m ); //!< define the actual generated mass

	///  setGeneratedMass() is included for backwards compatibility with CLHEP HepMC
        void   setGeneratedMass( const double & m )  
	                 { return set_generated_mass(m); }

    protected: // for internal use only by friend GenVertex class

	//static unsigned int counter(); //!< temporary for debugging

        /// set production vertex - for internal use only
	void   set_production_vertex_( GenVertex* productionvertex = 0);
        /// set decay vertex - for internal use only
	void   set_end_vertex_( GenVertex* decayvertex = 0 );
	void   set_barcode_( int the_bar_code ); //!< for use by GenEvent only

        /// scale the momentum vector and generated mass 
        /// this method is only for use by GenEvent
	void convert_momentum( const double& );

    private:
	FourVector       m_momentum;          // momentum vector
	int              m_pdg_id;            // id according to PDG convention
	int              m_status;            // As defined for HEPEVT
	Flow             m_flow;
	Polarization     m_polarization;
	GenVertex*       m_production_vertex; // null if vacuum or beam
	GenVertex*       m_end_vertex;        // null if not-decayed
	int              m_barcode;           // unique identifier in the event
        double           m_generated_mass;    // mass of this particle when it was generated

	//static unsigned int       s_counter;
    };  

    //////////////
    // INLINES  //
    //////////////

    inline GenParticle::operator HepMC::FourVector() const 
    { return m_momentum; }

    inline const FourVector & GenParticle::momentum() const 
    { return m_momentum; }

    inline int GenParticle::pdg_id() const { return m_pdg_id; }

    inline int GenParticle::status() const { return m_status; }

    inline GenVertex* GenParticle::production_vertex() const 
    { return m_production_vertex; }

    inline GenVertex* GenParticle::end_vertex() const { return m_end_vertex; }

    inline const Flow & GenParticle::flow() const { return m_flow; }

    inline int GenParticle::flow( int code_index ) const
    { return m_flow.icode( code_index ); }

    inline const Polarization & GenParticle::polarization() const 
    { return m_polarization; }

    inline void GenParticle::set_momentum( const FourVector& vec4 )
    { m_momentum = vec4; }

    inline void GenParticle::set_pdg_id( int id ) { m_pdg_id = id; }

    inline void GenParticle::set_status( int st ) { m_status = st; }

    inline void GenParticle::set_flow( const Flow& f ) { m_flow = f; }

    inline void GenParticle::set_flow( int code_index, int code ) 
    {
	if ( code == 0 ) { 
	    m_flow.set_unique_icode( code_index );
	} else { 
	    m_flow.set_icode( code_index, code );
	}
    }

    inline void GenParticle::set_polarization( const Polarization& polar )
    { m_polarization = polar; }

    inline int  GenParticle::barcode() const { return m_barcode; }

    inline void GenParticle::set_barcode_( int bc ) { m_barcode = bc; }

    inline bool GenParticle::is_undecayed() const {
        return ( m_status==1 ) ?  true : false;
    }
    inline bool GenParticle::has_decayed() const {
        return ( m_status==2 ) ?  true : false;
    }
    inline bool GenParticle::is_beam() const {
        return ( m_status==4 ) ?  true : false;
    }

} // HepMC

#endif  // HEPMC_GEN_PARTICLE_H
//--------------------------------------------------------------------------

