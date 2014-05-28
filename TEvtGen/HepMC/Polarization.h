//--------------------------------------------------------------------------
#ifndef HEPMC_POLARIZATION_H
#define HEPMC_POLARIZATION_H

//////////////////////////////////////////////////////////////////////////
// Matt.Dobbs@Cern.CH, September 1999, refer to:
// M. Dobbs and J.B. Hansen, "The HepMC C++ Monte Carlo Event Record for
// High Energy Physics", Computer Physics Communications (to be published).
//
// Polarization object for a particle. All angles are in radians.
//////////////////////////////////////////////////////////////////////////

#include "HepMC/SimpleVector.h"
#include <iostream>
#include <cmath>

namespace HepMC {

    static const double HepMC_pi = 3.14159265358979323846;  // copy of pi from CLHEP
    
    //! The Polarization class stores theta and phi for a GenParticle

    ///
    /// \class  Polarization
    /// HepMC::Polarization stores a particle's theta and phi in radians.
    /// Use of this information is optional. 
    /// By default, the polarization is set to zero.
    ///
    class Polarization {

        /// print polarization information
	friend std::ostream& operator<<( std::ostream&, const Polarization& );

    public:
        /// default constructor
	Polarization( );
	/// constructor requiring at least one value
	Polarization( double theta, double phi = 0 );
	/// construct from another polarization object
	Polarization( const Polarization& inpolar );
	/// construct using the polar and azimuthal angles from a ThreeVector
	Polarization( const ThreeVector& vec3in );
	virtual       ~Polarization() {}

        /// swap
        void swap( Polarization & other);
	/// make a copy
	Polarization& operator=( const Polarization& inpolar );
	/// equality requires that theta and phi are equal
	bool          operator==( const Polarization& ) const;
	/// inequality results if either theta or phi differ
	bool          operator!=( const Polarization& ) const;

        /// print theta and phi
	void          print( std::ostream& ostr = std::cout ) const;
    
	////////////////////
	// access methods //
	////////////////////
	double        theta() const;    //!< returns polar angle in radians
	double        phi() const;      //!< returns azimuthal angle in radians
	ThreeVector   normal3d() const; //!< unit 3 vector for easy manipulation
	bool          is_defined() const;   //!< returns true if the Polarization has been defined

        /// set polar angle in radians 
	double        set_theta( double theta );
        /// set azimuthal angle in radians 
	double        set_phi( double phi );
	/// set both polar and azimuthal angles in radians
	void          set_theta_phi( double theta, double phi );
	/// sets polarization according to direction of 3 vec
	ThreeVector   set_normal3d( const ThreeVector& vec3in ); 
	/// declares the Polarization as undefined and zeros the values
	void          set_undefined();

    private:
    	/// private method to return a polar angle in the correct range
	double valid_theta( double theta );
    	/// private method to return an azimuthal angle in the correct range
	double valid_phi( double phi );

    private:
	double m_theta; //polar angle of polarization in radians 0< theta <pi
	double m_phi;   //azimuthal angle of polarization in rad. 0< phi <2pi
	bool   m_defined; //used to flag if the Polarization has been defined
    };

    ///////////////////////////
    // INLINE Access Methods //
    ///////////////////////////

    inline double Polarization::theta() const { return m_theta; }
    inline double Polarization::phi() const { return m_phi; }

    ///////////////////////////
    // INLINE Operators      //
    ///////////////////////////

    inline bool Polarization::operator==( const Polarization& a ) const 
    {
	return ( a.theta() == this->theta() && a.phi() == this->phi() && a.is_defined() == this->is_defined() );
    }

    inline bool Polarization::operator!=(const Polarization& a ) const 
    {
	return !( a == *this );
    }

} // HepMC

#endif  // HEPMC_POLARIZATION_H
//--------------------------------------------------------------------------
