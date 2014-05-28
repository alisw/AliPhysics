//////////////////////////////////////////////////////////////////////////
// Matt.Dobbs@Cern.CH, November 2000, refer to:
// M. Dobbs and J.B. Hansen, "The HepMC C++ Monte Carlo Event Record for
// High Energy Physics", Computer Physics Communications (to be published).
//
// Container for the Weights associated with an event or vertex.
// Basically just an interface to STL vector with extra map-like attributes
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <stdexcept>

#include "HepMC/WeightContainer.h"

namespace HepMC {

WeightContainer::WeightContainer( size_type n, double value ) 
    : m_weights(n,value), m_names()
{ set_default_names(n); }

WeightContainer::WeightContainer( const std::vector<double>& wgts )
    : m_weights(wgts), m_names()
{ set_default_names(size()); }

void WeightContainer::set_default_names( size_type n )
{
    // internal program used by the constructors
    std::ostringstream name;
    for ( size_type count = 0; count<n; ++count ) 
    { 
	name.str(std::string());
	name << count;
	m_names[name.str()] = count;
    }
}

void WeightContainer::push_back( const double& value) 
{ 
    size_type count = m_weights.size();
    m_weights.push_back(value); 
    std::ostringstream name;
    name << count;
    m_names[name.str()] = count;
}

void WeightContainer::pop_back() 
{
    // this needs to remove the last entry in the vector 
    // and ALSO the associated map entry
    size_type vit = size() - 1;
    for ( map_iterator m = m_names.begin(); m != m_names.end(); ++m ) 
    { 
        if( m->second == vit ) { 
	    m_names.erase(m->first); 
	    continue;
	}
    }
    m_weights.pop_back(); 
}

double& WeightContainer::operator[]( const std::string& s ) 
{ 
    const_map_iterator m = m_names.find(s);
    if( m != m_names.end() ) {
        return m_weights[m->second]; 
    }
    // doesn't exist - have to create it
    size_type count = m_weights.size();
    m_weights.push_back(0); 
    m_names[s] = count;
    return m_weights.back(); 
}


const double& WeightContainer::operator[]( const std::string& s ) const
{ 
    const_map_iterator m = m_names.find(s);
    if( m != m_names.end() ) {
        return m_weights[m->second]; 
    }
    // doesn't exist and we cannot create it
    // note that std::map does not support this (const) operator
    // throw an appropriate error, we choose the error thrown by std::vector
    throw std::out_of_range("const WeightContainer::operator[] ERROR: string "+s+" not found in  WeightContainer" );
}

bool WeightContainer::operator==( const WeightContainer & other ) const
{
   if( size() != other.size() ) { return false; }
   if( m_names != other.m_names ) { return false; }
   if( m_weights != other.m_weights ) { return false; }
   return true;
}

bool WeightContainer::operator!=( const WeightContainer & other ) const
{
   return !(*this == other );
}

bool WeightContainer::has_key( const std::string& s ) const
{
    // look up the name in the map
    return m_names.find(s) != m_names.end();
}

void WeightContainer::print( std::ostream& ostr ) const 
{ 
    // print a name, weight pair
    for ( const_map_iterator m = map_begin(); m != map_end(); ++m )
    {
	ostr << "(" << m->first << "," << m_weights[m->second] << ") ";
    }
    ostr << std::endl; 
}

void WeightContainer::write( std::ostream& ostr ) const 
{ 
    size_type count = 0;
    for ( const_iterator w = begin(); w != end(); ++w ) 
    { 
	std::string name;
	for ( const_map_iterator m = map_begin(); m != map_end(); ++m )
	{
	    if( m->second == count ) name = m->first;
	}
	ostr << "Weight " << std::setw(4) << count 
	     << " with name " << std::setw(10) <<  name
	     << " is " << *w << std::endl;
	++count;
    }
}

} // HepMC

