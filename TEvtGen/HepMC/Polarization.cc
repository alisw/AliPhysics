//////////////////////////////////////////////////////////////////////////
// Matt.Dobbs@Cern.CH, September 1999
//
// Polarization object for a particle. All angles are in radians.
//////////////////////////////////////////////////////////////////////////

#include "HepMC/Polarization.h"

namespace HepMC {

    Polarization::Polarization( )
    : m_theta( 0. ),
      m_phi( 0. ),
      m_defined( false )
    { }

    Polarization::Polarization( double theta, double phi )
    : m_theta( valid_theta(theta) ),
      m_phi  ( valid_phi(phi) ),
      m_defined( true )
    { }

    Polarization::Polarization( const Polarization& inpolar )
    : m_theta( valid_theta( inpolar.theta() ) ),
      m_phi  ( valid_phi(   inpolar.phi()   ) ),
      m_defined( inpolar.is_defined() )
    { }

    Polarization::Polarization( const ThreeVector& vec3in ) 
    : m_theta( valid_theta( vec3in.theta() ) ),
      m_phi  ( valid_phi(   vec3in.phi()   ) ),
      m_defined( true )
    { }

    void Polarization::swap( Polarization & other)
    {
 	std::swap( m_theta, other.m_theta );
 	std::swap( m_phi,   other.m_phi   );
	std::swap( m_defined, other.m_defined );
    }

    Polarization& Polarization::operator=( const Polarization& inpolar ) {
        /// best practices implementation
	Polarization tmp( inpolar );
	swap( tmp ); 
	return *this;
    }

    void Polarization::print( std::ostream& ostr ) const {
	ostr << "Polarization: " << *this << std::endl;
    }

    ////////////////////
    // access methods //
    ////////////////////

    ThreeVector  Polarization::normal3d() const {
	// unit Hep3Vector for easy manipulation
	ThreeVector outvec(0,0,1);      // makes unit vector along Z
	outvec.setTheta( theta() ); // sets phi keeping mag and theta constant
	outvec.setPhi( phi() );     // sets theta keeping mag and phi constant
	return outvec;
    }

    double Polarization::set_theta( double theta ) {
	/// Theta is restricted to be between 0 --> pi
	/// if an out of range value is given, it is translated to this range.
	return m_theta = valid_theta( theta );
    }

    double Polarization::set_phi( double phi ) {
	/// Phi is restricted to be between 0 --> 2pi
	/// if an out of range value is given, it is translated to this range.
	return m_phi = valid_phi( phi );
    }
    
    bool Polarization::is_defined( ) const {
        return m_defined;
    }
    
    void Polarization::set_undefined() {
        m_defined = false;
	m_theta = 0.;
	m_phi = 0.;
    }

    void Polarization::set_theta_phi( double theta, double phi ) {
	set_theta( theta );
	set_phi( phi ) ;
	m_defined = true;
    }

    ThreeVector Polarization::set_normal3d( const ThreeVector& vec3in ) {
	set_theta( vec3in.theta() );
	set_phi( vec3in.phi() );
	m_defined = true;
	return vec3in;
    }

    /////////////////////
    // private methods //
    /////////////////////

    double Polarization::valid_theta( double theta ) {
        // this is just absolute value.
	theta = ( theta>0 ? theta : -theta );
	// translate to 0 < theta < 2pi
	theta = ( theta/(2*HepMC_pi) - int(theta/(2*HepMC_pi)) ) 
		* 2*HepMC_pi;
        // now translate to 0 < theta < pi
	if ( theta > HepMC_pi ) theta = 2*HepMC_pi - theta;
	return theta;
    }

    double Polarization::valid_phi( double phi ) {
	//
	// translate to -2pi < phi < 2pi
	phi = ( phi/(2*HepMC_pi) - int(phi/(2*HepMC_pi)) ) * 2*HepMC_pi;
	// translates to 0 < phi < 2pi
	if ( phi < 0 ) phi = 2*HepMC_pi + phi;
	return phi;
    }

    /////////////
    // Friends //
    /////////////

    /// write theta and phi to the output stream
    std::ostream& operator<<( std::ostream& ostr, const Polarization& polar ) {
	return ostr << "(" << polar.theta() 
		    << ","  << polar.phi() << ")";
    }
    
} // HepMC


