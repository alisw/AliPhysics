//--------------------------------------------------------------------------
#ifndef HEPMC_FLOW_H
#define HEPMC_FLOW_H

//////////////////////////////////////////////////////////////////////////
// Matt.Dobbs@Cern.CH, January 2000, refer to:
// M. Dobbs and J.B. Hansen, "The HepMC C++ Monte Carlo Event Record for
// High Energy Physics", Computer Physics Communications (to be published).
//
// particle's flow object
// keeps track of an arbitrary number of flow patterns within a graph 
// (i.e. color flow, charge flow, lepton number flow, ...) 
// Flow patterns are coded with an integer, in the same manner as in Herwig.
// Note: 0 is NOT allowed as code index nor as flow code since it
//       is used to indicate null.
//////////////////////////////////////////////////////////////////////////

// This class can be used to keep track of flow patterns within 
//  a graph. An example is color flow. If we have two quarks going through
//  an s-channel gluon to form two more quarks:
//
//  \q1       /q3   then we can keep track of the color flow with the
//   \_______/      HepMC::Flow class as follows: 
//   /   g   \. 
//  /q2       \q4
//
//  lets say the color flows from q2-->g-->q3  and q1-->g-->q4
//  the individual colors are unimportant, but the flow pattern is.
//  We can capture this flow by assigning the first pattern (q2-->g-->q3)
//  a unique (arbitrary) flow code 678 and the second pattern (q1-->g-->q4)
//  flow code 269  ( you can ask HepMC::Flow to choose
//  a unique code for you using Flow::set_unique_icode() ).
//  The first two code indices are reserved for color codes, so we store 
//  these codes with the particles as follows:
//    q2->flow().set_icode(1,678);
//    g->flow().set_icode(1,678);
//    q3->flow().set_icode(1,678);
//    q1->flow().set_icode(1,269);
//    g->flow().set_icode(2,269);
//    q4->flow().set_icode(1,269);
//  later on if we wish to know the color partner of q1 we can ask for a list
//  of all particles connected via this code to q1 which do have less than 
//  2 color partners using:
//    vector<GenParticle*> result=q1->dangling_connected_partners(q1->icode(1),1,2);
//  this will return a list containing q1 and q4.
//    vector<GenParticle*> result=q1->connected_partners(q1->icode(1),1,2);
//  would return a list containing q1, g, and q4.
//

#include <iostream>
#include <map>
#include <vector>

namespace HepMC {

    class GenParticle;

    //! The flow object

    ///
    /// \class  Flow
    /// The particle's flow object
    /// keeps track of an arbitrary number of flow patterns within a graph 
    /// (i.e. color flow, charge flow, lepton number flow, ...) 
    /// Flow patterns are coded with an integer, in the same manner as in Herwig.
    class Flow {

        /// for printing
	friend std::ostream& operator<<( std::ostream& ostr, const Flow& f );
	
    public:
        /// default constructor
	Flow( GenParticle* particle_owner = 0 );
	/// copy
	Flow( const Flow& );
	virtual         ~Flow();
        /// swap
        void swap( Flow & other);
	/// make a copy
	Flow&           operator=( const Flow& );
	/// equality
	bool            operator==( const Flow& a ) const; //compares only flow
	/// inequality
	bool            operator!=( const Flow& a ) const; //patterns not owner

        /// print Flow information to ostr
	void            print( std::ostream& ostr = std::cout ) const;

	/// returns all connected particles which have "code" in any  of the 
	///  num_indices beginning with index code_index.
	std::vector<HepMC::GenParticle*> connected_partners( int code, int code_index =1,
						   int num_indices = 2 ) const;
	/// same as connected_partners, but returns only those particles which
	///  are connected to <=1 other particles (i.e. the flow line "dangles"
	///  at these particles)
	std::vector<HepMC::GenParticle*> dangling_connected_partners( int code, 
			       int code_index = 1, int num_indices = 2 ) const;

	////////////////////
	// access methods //
	////////////////////

	/// find particle owning this Flow
	const GenParticle* particle_owner() const;
	/// flow code
	int             icode( int code_index = 1 ) const;
	/// set flow code
	Flow            set_icode( int code_index, int code );
	/// set unique flow code
	Flow            set_unique_icode( int code_index = 1 );

	//////////////////////
	// container access //
	//////////////////////

        /// return true if there is no flow container
	bool            empty() const;
	/// size of flow pattern container
	int             size() const;
	/// clear flow patterns
        void            clear();
	/// empty flow pattern container
	bool            erase( int code_index );

        /// iterator for flow pattern container
        typedef std::map<int,int>::iterator       iterator;
        /// const iterator for flow pattern container
        typedef std::map<int,int>::const_iterator const_iterator;
	/// beginning of flow pattern container
        iterator            begin();
	/// end of flow pattern container
        iterator            end();
	/// beginning of flow pattern container
        const_iterator      begin() const;
	/// end of flow pattern container
        const_iterator      end() const;

    protected: // intended for internal use only
        /// for internal use only
	void            connected_partners( std::vector<HepMC::GenParticle*>* output, 
					    int code,
					    int code_index,
					    int num_indices ) const;
        /// for internal use only
	void            dangling_connected_partners( std::vector<HepMC::GenParticle*>* 
						     output, 
						     std::vector<HepMC::GenParticle*>*
						     visited_particles, 
						     int code, int code_index, 
						     int num_indices ) const; 
    private:
	GenParticle*         m_particle_owner;
	std::map<int,int> m_icode; // stores flow patterns as(code_index,icode)
    };  

    ///////////////////////////
    // INLINE Access Methods //
    ///////////////////////////

    inline const GenParticle* Flow::particle_owner() const {
	return m_particle_owner;
    }
    inline int Flow::icode( int code_index ) const {
	std::map<int,int>::const_iterator a = m_icode.find(code_index);
	return a==m_icode.end() ? 0 : (*a).second;
    }
    inline Flow Flow::set_icode( int code_index, int code ) {
	m_icode[code_index] = code;
	return *this;
    }
    inline Flow Flow::set_unique_icode( int flow_num ) {
	/// use this method if you want to assign a unique flow code, but
	/// do not want the burden of choosing it yourself
	m_icode[flow_num] = size_t(this);
	return *this;
    }
    inline bool Flow::empty() const { return (bool)m_icode.empty(); }
    inline int Flow::size() const { return (int)m_icode.size(); }
    inline void Flow::clear() { m_icode.clear(); }
    inline bool Flow::erase( int code_index ) {
	// this will return true if the number of elements removed is nonzero
	return m_icode.erase( code_index )==0 ? false : true ;
    }
    inline Flow::iterator Flow::begin() { return m_icode.begin(); }
    inline Flow::iterator Flow::end() { return m_icode.end(); }
    inline Flow::const_iterator Flow::begin() const { return m_icode.begin(); }
    inline Flow::const_iterator Flow::end() const { return m_icode.end(); }

    ///////////////////////////
    // INLINE Operators      //
    ///////////////////////////

    inline bool Flow::operator==( const Flow& a ) const {
	/// equivalent flows have the same flow codes for all flow_numbers 
	/// (i.e. their m_icode maps are identical), but they need not have the
	/// same m_particle owner
	return (m_icode == a.m_icode);
    }
    inline bool Flow::operator!=( const Flow& a ) const {
	return !( *this == a );
    }
    inline Flow& Flow::operator=( const Flow& inflow ) {
	/// copies only the m_icode ... not the particle_owner
	/// this is intuitive behaviour so you can do
	/// oneparticle->flow() = otherparticle->flow()
	//
	m_icode = inflow.m_icode;
	return *this;
    }

} // HepMC

#endif  // HEPMC_FLOW_H
//--------------------------------------------------------------------------

