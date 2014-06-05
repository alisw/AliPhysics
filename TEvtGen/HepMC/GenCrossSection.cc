//--------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////
// garren@fnal.gov, January 2009
//
// The singleton GenCrossSection class holds run level information, 
// such as the cross section.
// 
//////////////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------

#include <iostream>
#include <sstream>

#include "HepMC/GenCrossSection.h"
#include "HepMC/IO_Exception.h"

namespace HepMC {

GenCrossSection::GenCrossSection( GenCrossSection const & orig )
    : m_cross_section( orig.cross_section() ),
      m_cross_section_error( orig.cross_section_error() ),
      m_is_set( orig.is_set() )
{}

void GenCrossSection::swap( GenCrossSection & other)
{
    std::swap( m_cross_section, other.m_cross_section );
    std::swap( m_cross_section_error, other.m_cross_section_error );
    std::swap( m_is_set, other.m_is_set );
}

GenCrossSection &  GenCrossSection::operator = ( GenCrossSection const & rhs )
{
    GenCrossSection tmp( rhs );
    swap( tmp );
    return *this;
}

bool GenCrossSection::operator==( const GenCrossSection& rhs ) const
{
    if( rhs.cross_section() != this->cross_section() ) return false;
    if( rhs.cross_section_error() != this->cross_section_error() ) return false;
    return true;
}

bool GenCrossSection::operator!=( const GenCrossSection& rhs ) const
{
	return !( rhs == *this );
}


void GenCrossSection::clear() 
{
    m_cross_section       = 0.0;
    m_cross_section_error = 0.0;
    m_is_set              = false;
}

std::ostream & GenCrossSection::write( std::ostream & os ) const
{
    // make sure the stream is valid
    if ( !os ) {
	std::cerr << "GenCrossSection::print !os, setting badbit" << std::endl;
	os.clear(std::ios::badbit); 
	return os;
    }
    // write the GenCrossSection information if the cross section was set
    if( is_set() ) {
	os << "C " << m_cross_section 
	   << " " << m_cross_section_error 
	   << "\n";
    }
    return os;
}

std::istream & GenCrossSection::read( std::istream & is )
{
    // make sure the stream is valid
    if ( !is ) { 
      std::cerr << "GenCrossSection stream input setting badbit." << std::endl;
      is.clear(std::ios::badbit); 
      return is; 
    }
    // check to see if we have a GenCrossSection line
    // This line is optional and may not exist
    if ( is.peek()!='C' ) { 
      return is; 
    }
    // get the GenCrossSection line
    std::string line, firstc;
    std::getline(is,line);
    std::istringstream iline(line);
    // Get first character and throw it away
    iline >> firstc;
    // Now get the numbers
    double xs = 0., xserr = 0.;
    iline >> xs ;
    if(!iline) throw IO_Exception("GenCrossSection::read encounterd invalid data");
    iline >> xserr ;
    if(!iline) throw IO_Exception("GenCrossSection::read encounterd invalid data");
    // set the data members
    set_cross_section( xs, xserr );
    return  is;
}

} // HepMC
